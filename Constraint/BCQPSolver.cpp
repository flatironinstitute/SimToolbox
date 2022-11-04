#include "BCQPSolver.hpp"
#include "Trilinos/TpetraUtil.hpp"
#include "spdlog/spdlog.h"

#include <mpi.h>
#include <omp.h>

#include <limits>
#include <random>

BCQPSolver::BCQPSolver(const Teuchos::RCP<const TOP> &A_, const Teuchos::RCP<const TV> &b_)
    : ARcp(A_), bRcp(b_), mapRcp(b_->getMap()), commRcp(b_->getMap()->getComm()) {
    // make sure A and b match the map and comm specified
    TEUCHOS_TEST_FOR_EXCEPTION(!(ARcp->getDomainMap()->isSameAs(*(bRcp->getMap()))), std::invalid_argument,
                               "A (domain) and b do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(mapRcp->isSameAs(*(bRcp->getMap()))), std::invalid_argument,
                               "map and b do not have the same Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(!(mapRcp->isSameAs(*(ARcp->getDomainMap()))), std::invalid_argument,
                               "map and b do not have the same Map.");
    setDefaultBounds();
}

BCQPSolver::BCQPSolver(int localSize, double diagonal) {
    // set up comm
    commRcp = getMPIWORLDTCOMM();
    // set up row and col maps, contiguous and evenly distributed
    Teuchos::RCP<const TMAP> rowMapRcp = getTMAPFromLocalSize(localSize, commRcp);
    mapRcp = rowMapRcp;

    if (commRcp->getRank() == 0) {
        std::cout << "Total number of processes: " << commRcp->getSize() << std::endl;
        std::cout << "rank: " << commRcp->getRank() << std::endl;
        std::cout << "global size: " << mapRcp->getGlobalNumElements() << std::endl;
        std::cout << "local size: " << mapRcp->getNodeNumElements() << std::endl;
        std::cout << "map: " << mapRcp->description() << std::endl;
    }

    // make sure A and b match the map and comm specified
    // set A and b randomly. maintain SPD of A
    Teuchos::RCP<TV> btemp = Teuchos::rcp(new TV(rowMapRcp, false));
    btemp->randomize(-1, 1);
    bRcp = btemp.getConst();

    // generate a local random matrix

    std::random_device rd;                       // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd());                      // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(-1, 1); // U(-1,1)

    // a random matrix B
    Teuchos::SerialDenseMatrix<int, double> BLocal(localSize, localSize, true); // zeroOut
    for (int i = 0; i < localSize; i++) {
        for (int j = 0; j < localSize; j++) {
            BLocal(i, j) = dis(gen);
        }
    }
    // a random diagonal matrix D. D entries all positive
    Teuchos::SerialDenseMatrix<int, double> ALocal(localSize, localSize, true);
    Teuchos::SerialDenseMatrix<int, double> tempLocal(localSize, localSize, true);
    Teuchos::SerialDenseMatrix<int, double> DLocal(localSize, localSize, true);
    for (int i = 0; i < localSize; i++) {
        DLocal(i, i) = pow(10, dis(gen)); // D in [10^(-1),10]
    }

    // compute B^T D B
    tempLocal.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, DLocal, BLocal, 0.0); // temp = DB
    ALocal.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, BLocal, tempLocal, 0.0);    // A = B^T DB

    // extra additions to the diagonal entries of A
    for (int i = 0; i < localSize; i++) {
        ALocal(i, i) += diagonal;
    }

    // use ALocal as local matrix to fill TCMAT A
    // block diagonal distribution of A
    double droptol = 1e-7;
    Kokkos::View<size_t *> rowCount("rowCount", localSize);
    Kokkos::View<size_t *> rowPointers("rowPointers", localSize + 1);
    for (int i = 0; i < localSize; i++) {
        rowCount[i] = 0;
        for (int j = 0; j < localSize; j++) {
            if (fabs(ALocal(i, j)) > droptol) {
                rowCount[i]++;
            }
        }
    }

    rowPointers[0] = 0;
    for (int i = 1; i < localSize + 1; i++) {
        rowPointers[i] = rowPointers[i - 1] + rowCount[i - 1];
    }
    Kokkos::View<int *> columnIndices("columnIndices", rowPointers[localSize]);
    Kokkos::View<double *> values("values", rowPointers[localSize]);
    int p = 0;
    for (int i = 0; i < localSize; i++) {
        for (int j = 0; j < localSize; j++) {
            if (fabs(ALocal(i, j)) > droptol) {
                columnIndices[p] = j;
                values[p] = ALocal(i, j);
                p++;
            }
        }
    }

    const int myRank = commRcp->getRank();
    const int colIndexCount = rowPointers[localSize];
    std::vector<int> colMapIndex(colIndexCount);
#pragma omp parallel for
    for (int i = 0; i < colIndexCount; i++) {
        colMapIndex[i] = columnIndices[i] + myRank * localSize;
    }

    // sort and unique
    std::sort(colMapIndex.begin(), colMapIndex.end());
    colMapIndex.erase(std::unique(colMapIndex.begin(), colMapIndex.end()), colMapIndex.end());

    Teuchos::RCP<TMAP> colMapRcp = Teuchos::rcp(
        new TMAP(Teuchos::OrdinalTraits<int>::invalid(), colMapIndex.data(), colMapIndex.size(), 0, commRcp));

    // fill matrix Aroot
    Teuchos::RCP<TCMAT> Atemp = Teuchos::rcp(new TCMAT(rowMapRcp, colMapRcp, rowPointers, columnIndices, values));
    Atemp->fillComplete(rowMapRcp, rowMapRcp);
    this->ARcp = Teuchos::rcp_dynamic_cast<const TOP>(Atemp, true);

    spdlog::info("Constraint operator ARcp is " + ARcp->description());

    // dump problem
    dumpTCMAT(Atemp, "Amat");
    dumpTV(bRcp, "bvec");

    generateRandomBounds();
}

int BCQPSolver::solveBBPGD(Teuchos::RCP<TV> &xsolRcp, const double tol, const int iteMax, IteHistory &history) const {
    // map must match
    TEUCHOS_TEST_FOR_EXCEPTION(!this->mapRcp->isSameAs(*(xsolRcp->getMap())), std::invalid_argument,
                               "xsolrcp and A operator do not have the same Map.");

    int mvCount = 0; // count matrix-vector multiplications
    int iteCount = 0;
    spdlog::debug("solving APGD");
    spdlog::debug("Constraint operator ARcp is "+ ARcp->description());

    Teuchos::RCP<TV> xkRcp = Teuchos::rcp(new TV(*xsolRcp, Teuchos::Copy));   // deep copy, xk=x0
    Teuchos::RCP<TV> xkm1Rcp = Teuchos::rcp(new TV(*xsolRcp, Teuchos::Copy)); // deep copy, xkm1=x0

    Teuchos::RCP<TV> gradkRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));   // the grad vector
    Teuchos::RCP<TV> gradkm1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true)); // the grad vector

    Teuchos::RCP<TV> gkdiffRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true)); // gkdiff = gk - gkm1
    Teuchos::RCP<TV> xkdiffRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true)); // xkdiff = xk - xkm1

    // compute grad
    ARcp->apply(*xkm1Rcp, *gradkm1Rcp); // gkm1 = A.dot(xkm1)
    mvCount++;
    gradkm1Rcp->update(1.0, *bRcp, 1.0); // gkm1 = A.dot(xkm1)+b

    // check if initial guess works, use xkdiffRcp as temporary space
    double resPhi = checkProjectionResidual(xkm1Rcp, gradkm1Rcp, xkdiffRcp);
    history.push_back(std::array<double, 6>{{1.0 * iteCount, 0, 0, 0, resPhi, 1.0 * mvCount}});
    if (fabs(resPhi) < tol) {
        // initial guess works, return
        xsolRcp = xkm1Rcp;
        return 0;
    }

    // first step, simple Gradient Descent stepsize = g^T g / g^T A g
    // use xkdiffRcp as temporary space
    // ARcp->apply(*gradkm1Rcp, *xkdiffRcp); // Avec = A * gkm1
    // mvCount++;

    // double gTAg = gradkm1Rcp->dot(*xkdiffRcp);
    // double gTg = pow(gradkm1Rcp->norm2(), 2);

    // if (fabs(gTAg) < 10 * std::numeric_limits<double>::epsilon()) {
    //     gTAg += 10 * std::numeric_limits<double>::epsilon(); // prevent div 0 error
    // }

    // double alpha = gTg / gTAg;

    // first step, Dai&Fletcher2005 Section 5.
    // xkdiffRcp is the prjected vector, after checkProjectionResidual
    double alpha = 1.0 / xkdiffRcp->normInf();

    bool stagFlag = false;

    while (iteCount < iteMax) {
        iteCount++;

        // update xk
        xkRcp->update(-alpha, *gradkm1Rcp, 1.0, *xkm1Rcp, 0.0); // xk = xkm1 - alpha*gkm1
        boundProjection(xkRcp);                                 // Projection xk

        // compute new grad with xk
        ARcp->apply(*xkRcp, *gradkRcp); // gk = A.dot(xk)
        mvCount++;
        gradkRcp->update(1.0, *bRcp, 1.0); // gk = A.dot(xk)+b

        // check convergence, use xkdiffRcp as temporary space
        double resPhi = checkProjectionResidual(xkRcp, gradkRcp, xkdiffRcp);

        // use simple phi tolerance check
        history.push_back(std::array<double, 6>{{1.0 * iteCount, 0, 0, alpha, resPhi, 1.0 * mvCount}});
        if (fabs(resPhi) < tol) {
            break;
        }
        xkdiffRcp->update(1.0, *xkRcp, -1.0, *xkm1Rcp, 0.0);       // xk - xkm1
        gkdiffRcp->update(1.0, *gradkRcp, -1.0, *gradkm1Rcp, 0.0); // gk - gkm1

        double a = 0, b = 0;

        // alternating bb1 and bb2 methods
        if (iteCount % 2 == 0) {
            // Barzilai-Borwein step size Choice 1
            a = pow(xkdiffRcp->norm2(), 2);
            b = xkdiffRcp->dot(*gkdiffRcp);
        } else {
            // Barzilai-Borwein step size Choice 2
            a = xkdiffRcp->dot(*gkdiffRcp);
            b = pow(gkdiffRcp->norm2(), 2);
        }

        if (fabs(b) < 10 * std::numeric_limits<double>::epsilon()) {
            b += 10 * std::numeric_limits<double>::epsilon(); // prevent div 0 error
        }

        alpha = a / b; // new step size

        if (alpha < std::numeric_limits<double>::epsilon() * 10) {
            spdlog::critical("BBPGD Stagnate");
            stagFlag = true;
            break;
        }

        // prepare next iteration
        // swap the contents of pointers directly, be careful
        xkm1Rcp.swap(xkRcp);
        gradkm1Rcp.swap(gradkRcp);
    }

    xsolRcp = xkRcp; // return solution
    if (stagFlag) {
        return 1;
    } else {
        return 0;
    }
}

int BCQPSolver::solveAPGD(Teuchos::RCP<TV> &xsolRcp, const double tol, const int iteMax, IteHistory &history) const {
    // map must match
    TEUCHOS_TEST_FOR_EXCEPTION(!this->mapRcp->isSameAs(*(xsolRcp->getMap())), std::invalid_argument,
                               "xsolrcp and A operator do not have the same Map.")

    int mvCount = 0;
    spdlog::debug("solving APGD");
    spdlog::debug("Constraint operator ARcp is "+ ARcp->description());

    // allocate vectors
    Teuchos::RCP<TV> xkRcp = Teuchos::rcp(new TV(*xsolRcp, Teuchos::Copy)); // deep copy
    Teuchos::RCP<TV> ykRcp = Teuchos::rcp(new TV(*xsolRcp, Teuchos::Copy)); // deep copy, yk=xk

    Teuchos::RCP<TV> xkp1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));
    Teuchos::RCP<TV> ykp1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));

    Teuchos::RCP<TV> gVecRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));

    Teuchos::RCP<TV> tempVecRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true)); // temporary result holder

    Teuchos::RCP<TV> xhatkRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));
    xhatkRcp->putScalar(1.0);

    double thetak = 1;
    double thetakp1 = 1;

    Teuchos::RCP<TV> xkdiffRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));
    xkdiffRcp->update(-1.0, *xhatkRcp, 1.0, *xkRcp, 0.0);

    ARcp->apply(*xkdiffRcp, *tempVecRcp);
    mvCount++;

    const double tempNorm2 = tempVecRcp->norm2();
    const double xkdiffNorm2 = xkdiffRcp->norm2();
    double Lk = (tempNorm2 / xkdiffNorm2);
    double tk = 1.0 / Lk;

    Teuchos::RCP<TV> AxbRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), false));
    Teuchos::RCP<TV> Axbkp1Rcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true));

    history.push_back(std::array<double, 6>{{0, 0, 0, tk, 0, 1.0 * mvCount}});

    // enter main loop
    int iteCount = 0;
    double resmin = std::numeric_limits<double>::max();
    bool stagFlag = false;

    while (iteCount < iteMax) {
        iteCount++;

        // line 7 of Mazhar, 2015, g=A.dot(yk)+b
        ARcp->apply(*ykRcp, *AxbRcp); // Axb = A yk, this does not change in the following Lifshitz loop
        mvCount++;
        gVecRcp->update(1.0, *bRcp, 1.0, *AxbRcp, 0.0);
        // line 8 of Mazhar, 2015
        xkp1Rcp->update(1.0, *ykRcp, -tk, *gVecRcp, 0);
        boundProjection(xkp1Rcp);

        double rightTerm1 = ykRcp->dot(*AxbRcp) * 0.5; // yk.dot(A.dot(yk))*0.5
        double rightTerm2 = ykRcp->dot(*bRcp);         // yk.dot(b)

        while (1) {
            //  xkdiff=xkp1-yk
            xkdiffRcp->update(1.0, *xkp1Rcp, -1.0, *ykRcp, 0.0);

            //  calc Lifshitz condition
            ARcp->apply(*xkp1Rcp, *Axbkp1Rcp);
            mvCount++;
            double leftTerm1 = xkp1Rcp->dot(*Axbkp1Rcp) * 0.5; // xkp1.dot(A.dot(xkp1))*0.5
            double leftTerm2 = xkp1Rcp->dot(*bRcp);            // xkp1.dot(b)

            double rightTerm3 = gVecRcp->dot(*xkdiffRcp);              // g.dot(xkdiff)
            double rightTerm4 = 0.5 * Lk * pow(xkdiffRcp->norm2(), 2); // 0.5*Lk*(xkdiff).dot(xkdiff)
            if ((leftTerm1 + leftTerm2) <= (rightTerm1 + rightTerm2 + rightTerm3 + rightTerm4)) {
                break;
            }
            // line 10 & 11 of Mazhar, 2015
            Lk *= 2;
            tk = 1 / Lk;

            // print Lk and tk for debugging
            // std::cout << Lk << " " << tk << std::endl;

            // line 12 of Mazhar, 2015
            xkp1Rcp->update(1.0, *ykRcp, -tk, *gVecRcp, 0.0);
            boundProjection(xkp1Rcp);
        }

        if (tk < std::numeric_limits<double>::epsilon() * 10) {
            spdlog::critical("APGD Stagnate");
            stagFlag = true;
            break;
        }

        // line 14-16, Mazhar, 2015
        thetakp1 = (-thetak * thetak + thetak * sqrt(4 + thetak * thetak)) / 2;
        double betakp1 = thetak * (1 - thetak) / (thetak * thetak + thetakp1);
        // ykp1=xkp1+betakp1*(xkp1-xk)
        ykp1Rcp->update((1 + betakp1), *xkp1Rcp, -betakp1, *xkRcp, 0);

        // check convergence, line 17, Mazhar, 2015. Replace the metric with the minimum-map function
        Axbkp1Rcp->update(1.0, *bRcp, 1.0); // Axkp1 = A*xkp1 in the Lifshitz loop
        double resPhi = checkProjectionResidual(xkp1Rcp, Axbkp1Rcp, tempVecRcp);
        resPhi = fabs(resPhi);

        // line 18-21, Mazhar, 2015
        if (resPhi < resmin) {
            resmin = resPhi;
            xhatkRcp->update(1.0, *xkp1Rcp, 0.0);
        }

        // line 22-24, Mazhar, 2015
        history.push_back(std::array<double, 6>{{1.0 * iteCount, 0, 0, tk, resPhi, 1.0 * mvCount}});
        if (resPhi < tol) {
            break;
        }

        // line 25-28, Mazhar, 2015
        tempVecRcp->update(1.0, *xkp1Rcp, -1.0, *xkRcp, 0.0);
        if (gVecRcp->dot(*tempVecRcp) > 0) {
            ykp1Rcp->scale(1.0, *xkp1Rcp); // ykp1=xkp1
            thetakp1 = 1;
        }

        // line 29-30, Mazhar, 2015
        Lk *= 0.9;
        tk = 1 / Lk;

        // next iteration
        // swap the contents of pointers directly, be careful
        ykRcp.swap(ykp1Rcp); // yk=ykp1, ykp1 to be updated;
        xkRcp.swap(xkp1Rcp); // xk=xkp1, xkp1 to be updated;
        thetak = thetakp1;
    }
    xsolRcp = xhatkRcp;
    if (stagFlag) {
        return 1;
    } else {
        return 0;
    }
}

int BCQPSolver::selfTest(double tol, int maxIte, int solverChoice) {
    IteHistory history;

    Teuchos::RCP<TV> xsolRcp = Teuchos::rcp(new TV(this->mapRcp.getConst(), true)); // zero initial guess
    prepareSolver();

    // dump problem
    dumpTV(lbRcp, "lbvec");
    dumpTV(ubRcp, "ubvec");

    spdlog::info("START TEST");

    switch (solverChoice) {
    case 1:
        solveAPGD(xsolRcp, tol, maxIte, history);
        dumpTV(xsolRcp, "xsolAPGD");
        break;
    default:
        solveBBPGD(xsolRcp, tol, maxIte, history);
        dumpTV(xsolRcp, "xsolBBPGD");
        break;
    }

    // dump iterative history to csv format
    if (commRcp->getRank() == 0)
        for (const auto &record : history) {
            if (solverChoice == 1) {
                printf("APGD_HISTORY,");
            } else {
                printf("BBPGD_HISTORY,");
            }
            for (const auto &v : record) {
                printf("%.6g, ", v);
            }
            printf("\n");
        }

    return 0;
}

void BCQPSolver::boundProjection(Teuchos::RCP<TV> &vecRcp) const {

    auto vecPtr = vecRcp->getLocalView<Kokkos::HostSpace>(); // LeftLayout
    vecRcp->modify<Kokkos::HostSpace>();
    const int ibound = vecPtr.dimension_0();
    const int c = 0; // vecRcp, lbRcp, ubRcp have only 1 column

    // project to lb
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecRcp->getMap()->isSameAs(*(lbRcp->getMap()))), std::invalid_argument,
                               "vec and lb do not have the same Map.");
    auto lbPtr = lbRcp->getLocalView<Kokkos::HostSpace>();
#pragma omp parallel for
    for (int i = 0; i < ibound; i++) {
        const double temp = vecPtr(i, c);
        vecPtr(i, c) = std::max(temp, lbPtr(i, c));
    }

    // project to ub
    TEUCHOS_TEST_FOR_EXCEPTION(!(vecRcp->getMap()->isSameAs(*(ubRcp->getMap()))), std::invalid_argument,
                               "vec and ub do not have the same Map.");
    auto ubPtr = ubRcp->getLocalView<Kokkos::HostSpace>();
#pragma omp parallel for
    for (int i = 0; i < ibound; i++) {
        const double temp = vecPtr(i, c);
        vecPtr(i, c) = std::min(temp, ubPtr(i, c));
    }

    return;
}

double BCQPSolver::checkProjectionResidual(const Teuchos::RCP<const TV> &XRcp, const Teuchos::RCP<const TV> &YRcp,
                                           const Teuchos::RCP<TV> &QRcp) const {
    const double eps = std::numeric_limits<double>::epsilon() * 100;
    auto xPtr = XRcp->getLocalView<Kokkos::HostSpace>();   // LeftLayout
    auto yPtr = YRcp->getLocalView<Kokkos::HostSpace>();   // LeftLayout
    auto lbPtr = lbRcp->getLocalView<Kokkos::HostSpace>(); // LeftLayout
    auto ubPtr = ubRcp->getLocalView<Kokkos::HostSpace>(); // LeftLayout
    auto qPtr = QRcp->getLocalView<Kokkos::HostSpace>();   // LeftLayout
    QRcp->modify<Kokkos::HostSpace>();
    const int ibound = xPtr.dimension_0();
    const int c = 0; // vecRcp, lbRcp, ubRcp have only 1 column

    bool projectionError = false;
// EQ 2.2 of Dai & Fletcher 2005
#pragma omp parallel for
    for (int i = 0; i < ibound; i++) {
        if (xPtr(i, c) < lbPtr(i, c) + eps) {
            qPtr(i, c) = std::min(yPtr(i, c), 0.0);
        } else if (xPtr(i, c) > ubPtr(i, c) - eps) {
            qPtr(i, c) = std::max(yPtr(i, c), 0.0);
        } else if (xPtr(i, c) > lbPtr(i, c) && xPtr(i, c) < ubPtr(i, c)) {
            qPtr(i, c) = yPtr(i, c);
        } else {
            spdlog::error("projection error occured");
            projectionError = true;
        }
    }

    if (projectionError) {
        dumpTV(XRcp, "XRcp");
        dumpTV(YRcp, "YRcp");
        dumpTV(QRcp, "QRcp");
        std::exit(1);
    }

    return QRcp->normInf();
}

void BCQPSolver::setDefaultBounds() {
    if (!lbSet) {
        const auto &vec = Teuchos::rcp(new TV(bRcp->getMap(), false));
        vec->putScalar(-std::numeric_limits<double>::max() / 10);
        setLowerBound(vec);
    }
    if (!ubSet) {
        const auto &vec = Teuchos::rcp(new TV(bRcp->getMap(), false));
        vec->putScalar(std::numeric_limits<double>::max() / 10);
        setUpperBound(vec);
    }
}

void BCQPSolver::generateRandomBounds() {
    const auto &vec1 = Teuchos::rcp(new TV(bRcp->getMap(), true));
    vec1->randomize(-1, 1);

    const auto &vec2 = Teuchos::rcp(new TV(bRcp->getMap(), true));
    vec2->randomize(-1, 1);

    auto vec1Ptr = vec1->getLocalView<Kokkos::HostSpace>(); // LeftLayout
    auto vec2Ptr = vec2->getLocalView<Kokkos::HostSpace>(); // LeftLayout
    vec1->modify<Kokkos::HostSpace>();
    vec2->modify<Kokkos::HostSpace>();
    const int ibound = vec1Ptr.dimension_0();
#pragma omp parallel for
    for (int i = 0; i < ibound; i++) {
        double a = vec1Ptr(i, 0);
        double b = vec2Ptr(i, 0);
        vec1Ptr(i, 0) = std::min(a, b);
        vec2Ptr(i, 0) = std::max(a, b);
    }

    setLowerBound(vec1);
    setUpperBound(vec2);
}
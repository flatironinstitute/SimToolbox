#include "SylinderNear.hpp"
#include "Util/EigenDef.hpp"
#include <cmath>

void printMat3(const double v[9]) {
    const Emat3 mat = Eigen::Map<const Emat3>(v);
    std::cout << mat << std::endl;
}

void printVec3(const double v[3]) { std::cout << ECmap3(v).transpose() << std::endl; }

void testEpsilon() {
    constexpr double epsilonConst[3][3][3] = {{{0, 0, 0}, {0, 0, 1}, {0, -1, 0}}, //
                                              {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}}, //
                                              {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};
    double epsilon[3][3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                epsilon[i][j][k] = 0.0;
            }
        }
    }
    epsilon[0][1][2] = 1.0;
    epsilon[1][2][0] = 1.0;
    epsilon[2][0][1] = 1.0;
    epsilon[1][0][2] = -1.0;
    epsilon[2][1][0] = -1.0;
    epsilon[0][2][1] = -1.0;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                if (epsilon[i][j][k] != epsilonConst[i][j][k]) {
                    printf("epsilonConst not correct\n");
                    std::exit(1);
                }
            }
        }
    }
}

void testFixedPair() {
    omp_set_num_threads(1);
    Evec3 P0(1, 0, 0);
    Evec3 P1(0, sqrt(3), 0);
    Evec3 Q0(0, 0, 1);
    Evec3 Q1(2, 2 * sqrt(3), 1);

    CalcSylinderNearForce calc;
    calc.conPoolPtr = std::make_shared<ConstraintBlockPool>();
    calc.conPoolPtr->resize(1);

    std::vector<SylinderNearEP> sylinderP(1);
    std::vector<SylinderNearEP> sylinderQ(1);

    { // setup P
        Evec3 center = (P0 + P1) / 2;
        Evec3 direction = (P1 - P0).normalized();
        double length = (P1 - P0).norm();
        sylinderP[0].gid = 0;
        sylinderP[0].globalIndex = 0;
        sylinderP[0].rank = 0;
        sylinderP[0].radiusCollision = 0.4;
        sylinderP[0].lengthCollision = length;
        sylinderP[0].pos[0] = center[0];
        sylinderP[0].pos[1] = center[1];
        sylinderP[0].pos[2] = center[2];
        sylinderP[0].direction[0] = direction[0];
        sylinderP[0].direction[1] = direction[1];
        sylinderP[0].direction[2] = direction[2];
    }
    { // setup Q
        Evec3 center = (Q0 + Q1) / 2;
        Evec3 direction = (Q1 - Q0).normalized();
        double length = (Q1 - Q0).norm();
        sylinderQ[0].gid = 1;
        sylinderQ[0].globalIndex = 1;
        sylinderQ[0].rank = 0;
        sylinderQ[0].radiusCollision = 0.5;
        sylinderQ[0].lengthCollision = length;
        sylinderQ[0].pos[0] = center[0];
        sylinderQ[0].pos[1] = center[1];
        sylinderQ[0].pos[2] = center[2];
        sylinderQ[0].direction[0] = direction[0];
        sylinderQ[0].direction[1] = direction[1];
        sylinderQ[0].direction[2] = direction[2];
    }
    ForceNear fnear;
    calc(sylinderP.data(), 1, sylinderQ.data(), 1, &fnear);
    printf("%zu collisions recorded\n", calc.conPoolPtr->front().size());
    auto stress = calc.conPoolPtr->front().front().stress;
    printf("stress:\n");
    printMat3(stress);
    bool pass = true;
    pass = pass && fabs(stress[0] - 0) < 1e-6;
    pass = pass && fabs(stress[1] - 0) < 1e-6;
    pass = pass && fabs(stress[2] - 0.0160681) < 1e-6;
    pass = pass && fabs(stress[3] - 0) < 1e-6;
    pass = pass && fabs(stress[4] - 0) < 1e-6;
    pass = pass && fabs(stress[5] - 0.0278307) < 1e-6;
    pass = pass && fabs(stress[6] - 0.0160681) < 1e-6;
    pass = pass && fabs(stress[7] - 0.0278307) < 1e-6;
    pass = pass && fabs(stress[8] - 1) < 1e-6;

    if (pass) {
    } else {
        printf("stress wrong\n");
        std::exit(1);
    }
}

void testParallel() {
    omp_set_num_threads(1);
    Evec3 P0(0.9942362769247484, 1.486372506783908, 1.499998650123065);
    Evec3 P1(1.994235986439871, 1.487134714472271, 1.500001668306246);
    Evec3 Q0(1.005687630910662, 1.464364624867677, 1.500001338752708);
    Evec3 Q1(2.005685087440173, 1.462109203504082, 1.4999983414464);

    CalcSylinderNearForce calc;
    calc.conPoolPtr = std::make_shared<ConstraintBlockPool>();
    calc.conPoolPtr->resize(1);

    std::vector<SylinderNearEP> sylinderP(1);
    std::vector<SylinderNearEP> sylinderQ(1);

    { // setup P
        Evec3 center = (P0 + P1) / 2;
        Evec3 direction = (P1 - P0).normalized();
        double length = (P1 - P0).norm();
        sylinderP[0].gid = 0;
        sylinderP[0].globalIndex = 0;
        sylinderP[0].rank = 0;
        sylinderP[0].radiusCollision = 0.4;
        sylinderP[0].lengthCollision = length;
        sylinderP[0].pos[0] = center[0];
        sylinderP[0].pos[1] = center[1];
        sylinderP[0].pos[2] = center[2];
        sylinderP[0].direction[0] = direction[0];
        sylinderP[0].direction[1] = direction[1];
        sylinderP[0].direction[2] = direction[2];
    }
    { // setup Q
        Evec3 center = (Q0 + Q1) / 2;
        Evec3 direction = (Q1 - Q0).normalized();
        double length = (Q1 - Q0).norm();
        sylinderQ[0].gid = 1;
        sylinderQ[0].globalIndex = 1;
        sylinderQ[0].rank = 0;
        sylinderQ[0].radiusCollision = 0.5;
        sylinderQ[0].lengthCollision = length;
        sylinderQ[0].pos[0] = center[0];
        sylinderQ[0].pos[1] = center[1];
        sylinderQ[0].pos[2] = center[2];
        sylinderQ[0].direction[0] = direction[0];
        sylinderQ[0].direction[1] = direction[1];
        sylinderQ[0].direction[2] = direction[2];
    }
    ForceNear fnear;
    calc(sylinderP.data(), 1, sylinderQ.data(), 1, &fnear);
    printf("%zu collisions recorded\n", calc.conPoolPtr->front().size());
    auto stress = calc.conPoolPtr->front().front().stress;
    auto block = calc.conPoolPtr->front().front();
    printf("stress:\n");
    printMat3(stress);
    printf("unscaledForceComI:\n");
    printVec3(block.unscaledForceComI);
    printf("unscaledForceComJ:\n");
    printVec3(block.unscaledForceComJ);
}

void testSphere() {
    omp_set_num_threads(1);
    Evec3 P0(1, 1, 1);
    Evec3 P1 = P0 + Evec3::Random();
    Evec3 Q0 = P0 + Evec3(1, 2., 3.);
    Evec3 Q1 = Q0 + Evec3::Random();

    CalcSylinderNearForce calc;
    calc.conPoolPtr = std::make_shared<ConstraintBlockPool>();
    calc.conPoolPtr->resize(1);

    std::vector<SylinderNearEP> sylinderP(1);
    std::vector<SylinderNearEP> sylinderQ(1);

    { // setup P
        Evec3 center = (P0 + P1) / 2;
        Evec3 direction = (P1 - P0).normalized();
        double length = (P1 - P0).norm();
        sylinderP[0].gid = 0;
        sylinderP[0].globalIndex = 0;
        sylinderP[0].rank = 0;
        sylinderP[0].radiusCollision = length;
        sylinderP[0].lengthCollision = length;
        sylinderP[0].pos[0] = center[0];
        sylinderP[0].pos[1] = center[1];
        sylinderP[0].pos[2] = center[2];
        sylinderP[0].direction[0] = direction[0];
        sylinderP[0].direction[1] = direction[1];
        sylinderP[0].direction[2] = direction[2];
    }
    { // setup Q
        Evec3 center = (Q0 + Q1) / 2;
        Evec3 direction = (Q1 - Q0).normalized();
        double length = (Q1 - Q0).norm();
        sylinderQ[0].gid = 1;
        sylinderQ[0].globalIndex = 1;
        sylinderQ[0].rank = 0;
        sylinderQ[0].radiusCollision = length;
        sylinderQ[0].lengthCollision = length;
        sylinderQ[0].pos[0] = center[0];
        sylinderQ[0].pos[1] = center[1];
        sylinderQ[0].pos[2] = center[2];
        sylinderQ[0].direction[0] = direction[0];
        sylinderQ[0].direction[1] = direction[1];
        sylinderQ[0].direction[2] = direction[2];
    }
    ForceNear fnear;
    calc(sylinderP.data(), 1, sylinderQ.data(), 1, &fnear);
    printf("%zu collisions recorded\n", calc.conPoolPtr->front().size());
    auto stress = calc.conPoolPtr->front().front().stress;
    auto block = calc.conPoolPtr->front().front();
    printf("stress:\n");
    printMat3(stress);
    printf("unscaledForceComI:\n");
    printVec3(block.unscaledForceComI);
    printf("unscaledForceComJ:\n");
    printVec3(block.unscaledForceComJ);

    // check correctness, stress = xf for spheres
    Emat3 stress1;
    block.getStress(stress1);
    Evec3 rij = Emap3(sylinderQ[0].pos) - Emap3(sylinderP[0].pos);
    Evec3 rijunit = rij.normalized();
    Emat3 stress2 = rij * rijunit.transpose();
    Emat3 error = (stress1 - stress2).cwiseAbs();
    if (error.maxCoeff() > 1e-7) {
        printf("sphere-sphere stress wrong\n");
        std::exit(1);
    }
}

void testSylinderSphere() {
    // P is sylinder
    // Q is sphere

    omp_set_num_threads(1);
    Evec3 P0(0, 0, 0);
    Evec3 P1 = P0 + Evec3(0.5, 0, 0);
    Evec3 Q0(0, 0, 1);
    Evec3 Q1 = Q0 + Evec3(1, 0, 0);

    CalcSylinderNearForce calc;
    calc.conPoolPtr = std::make_shared<ConstraintBlockPool>();
    calc.conPoolPtr->resize(1);

    std::vector<SylinderNearEP> sylinderP(1);
    std::vector<SylinderNearEP> sylinderQ(1);

    { // setup P
        Evec3 center = (P0 + P1) / 2;
        Evec3 direction = (P1 - P0).normalized();
        double length = (P1 - P0).norm();
        sylinderP[0].gid = 0;
        sylinderP[0].globalIndex = 0;
        sylinderP[0].rank = 0;
        sylinderP[0].radiusCollision = length;
        sylinderP[0].lengthCollision = length;
        sylinderP[0].pos[0] = center[0];
        sylinderP[0].pos[1] = center[1];
        sylinderP[0].pos[2] = center[2];
        sylinderP[0].direction[0] = direction[0];
        sylinderP[0].direction[1] = direction[1];
        sylinderP[0].direction[2] = direction[2];
    }
    { // setup Q
        Evec3 center = (Q0 + Q1) / 2;
        Evec3 direction = (Q1 - Q0).normalized();
        double length = (Q1 - Q0).norm();
        sylinderQ[0].gid = 1;
        sylinderQ[0].globalIndex = 1;
        sylinderQ[0].rank = 0;
        sylinderQ[0].radiusCollision = 0.4 * length;
        sylinderQ[0].lengthCollision = length;
        sylinderQ[0].pos[0] = center[0];
        sylinderQ[0].pos[1] = center[1];
        sylinderQ[0].pos[2] = center[2];
        sylinderQ[0].direction[0] = direction[0];
        sylinderQ[0].direction[1] = direction[1];
        sylinderQ[0].direction[2] = direction[2];
    }
    ForceNear fnear;
    calc(sylinderP.data(), 1, sylinderQ.data(), 1, &fnear);
    printf("%zu collisions recorded\n", calc.conPoolPtr->front().size());
    auto stress = calc.conPoolPtr->front().front().stress;
    auto block = calc.conPoolPtr->front().front();
    printf("stress:\n");
    printMat3(stress);
    printf("unscaledForceComI:\n");
    printVec3(block.unscaledForceComI);
    printf("unscaledForceComJ:\n");
    printVec3(block.unscaledForceComJ);
}

int main() {
    printf("--------------------testing epsilon tensor\n");
    testEpsilon();
    printf("-------------testing error of a given pair\n");
    testFixedPair();
    printf("---testing stability of parallel sylinders\n");
    testParallel();
    printf("---------------------------------------------\ntesting short sylinders as spheres\n");
    testSphere();
    printf("---------------------------------------------\ntesting sylinder-sphere \n");
    testSylinderSphere();
    return 0;
}
#include "SylinderSystem.hpp"

#include "MPI/CommMPI.hpp"
#include "Util/EquatnHelper.hpp"
#include "Util/GeoUtil.hpp"
#include "Util/IOHelper.hpp"
#include "Util/Logger.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <random>
#include <vector>
#include <cmath>

#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTypeInt32Array.h>
#include <vtkTypeUInt8Array.h>
#include <vtkXMLPPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>

#include <mpi.h>
#include <omp.h>

SylinderSystem::SylinderSystem(const std::string &configFile, const std::string &posFile, int argc, char **argv) {
    initialize(SylinderConfig(configFile), posFile, argc, argv);
}

SylinderSystem::SylinderSystem(const SylinderConfig &runConfig_, const std::string &posFile, int argc, char **argv) {
    initialize(runConfig_, posFile, argc, argv);
}

void SylinderSystem::initialize(const SylinderConfig &runConfig_, const std::string &posFile, int argc, char **argv) {
    runConfig = runConfig_;
    stepCount = 0;
    snapID = 0; // the first snapshot starts from 0 in writeResult

    // store the random seed
    restartRngSeed = runConfig.rngSeed;

    // set MPI
    int mpiflag;
    MPI_Initialized(&mpiflag);
    TEUCHOS_ASSERT(mpiflag);

    Logger::set_level(runConfig.logLevel);
    commRcp = getMPIWORLDTCOMM();

    showOnScreenRank0();

    // TRNG pool must be initialized after mpi is initialized
    rngPoolPtr = std::make_shared<TRngPool>(runConfig.rngSeed);
    conSolverPtr = std::make_shared<ConstraintSolver>();
    conCollectorPtr = std::make_shared<ConstraintCollector>();

    dinfo.initialize(); // init DomainInfo
    setDomainInfo();

    sylinderContainer.initialize();
    sylinderContainer.setAverageTargetNumberOfSampleParticlePerProcess(200); // more sample for better balance

    if (IOHelper::fileExist(posFile)) {
        setInitialFromFile(posFile);
    } else {
        setInitialFromConfig();
    }
    setLinkMapsFromFile(posFile);

    // at this point all sylinders located on rank 0
    commRcp->barrier();
    decomposeDomain();
    exchangeSylinder(); // distribute to ranks, initial domain decomposition

    sylinderNearDataDirectoryPtr = std::make_shared<ZDD<SylinderNearEP>>(sylinderContainer.getNumberOfParticleLocal());

    treeSylinderNumber = 0;
    setTreeSylinder();

    calcVolFrac();

    if (commRcp->getRank() == 0) {
        IOHelper::makeSubFolder("./result"); // prepare the output directory
        writeBox();
    }

    if (!runConfig.sylinderFixed) {
        // 100 NON-B steps to resolve initial configuration collisions
        // no output
        spdlog::warn("Initial Collision Resolution Begin");
        for (int i = 0; i < runConfig.initPreSteps; i++) {
            prepareStep();
            calcVelocityNonCon();
            resolveConstraints();
            saveForceVelocityConstraints();
            sumForceVelocity();
            stepEuler();
        }
        spdlog::warn("Initial Collision Resolution End");
    }

    spdlog::warn("SylinderSystem Initialized. {} local sylinders", sylinderContainer.getNumberOfParticleLocal());
}

void SylinderSystem::reinitialize(const SylinderConfig &runConfig_, const std::string &restartFile, int argc,
                                  char **argv, bool eulerStep) {
    runConfig = runConfig_;

    // Read the timestep information and pvtp filenames from restartFile
    std::string pvtpFileName;
    std::ifstream myfile(restartFile);

    myfile >> restartRngSeed;
    myfile >> stepCount;
    myfile >> snapID;
    myfile >> pvtpFileName;

    // increment the rngSeed forward by one to ensure randomness compared to previous run
    restartRngSeed++;

    // set MPI
    int mpiflag;
    MPI_Initialized(&mpiflag);
    TEUCHOS_ASSERT(mpiflag);

    Logger::set_level(runConfig.logLevel);
    commRcp = getMPIWORLDTCOMM();

    showOnScreenRank0();

    // TRNG pool must be initialized after mpi is initialized
    rngPoolPtr = std::make_shared<TRngPool>(restartRngSeed);
    conSolverPtr = std::make_shared<ConstraintSolver>();
    conCollectorPtr = std::make_shared<ConstraintCollector>();

    dinfo.initialize(); // init DomainInfo
    setDomainInfo();

    sylinderContainer.initialize();
    sylinderContainer.setAverageTargetNumberOfSampleParticlePerProcess(200); // more samples for better balance

    std::string asciiFileName = pvtpFileName;
    auto pos = asciiFileName.find_last_of('.');
    asciiFileName.replace(pos, 5, std::string(".dat")); // replace '.pvtp' with '.dat'
    pos = asciiFileName.find_last_of('_');
    asciiFileName.replace(pos, 1, std::string("Ascii_")); // replace '_' with 'Ascii_'

    std::string baseFolder = getCurrentResultFolder();
    setInitialFromVTKFile(baseFolder + pvtpFileName);

    setLinkMapsFromFile(baseFolder + asciiFileName);

    // VTK data is wrote before the Euler step, thus we need to run one Euler step below
    if (eulerStep)
        stepEuler();

    stepCount++;
    snapID++;

    // at this point all sylinders located on rank 0
    commRcp->barrier();
    applyBoxBC();
    decomposeDomain();
    exchangeSylinder(); // distribute to ranks, initial domain decomposition
    updateSylinderMap();

    sylinderNearDataDirectoryPtr = std::make_shared<ZDD<SylinderNearEP>>(sylinderContainer.getNumberOfParticleLocal());

    treeSylinderNumber = 0;
    setTreeSylinder();
    calcVolFrac();

    spdlog::warn("SylinderSystem Initialized. {} local sylinders", sylinderContainer.getNumberOfParticleLocal());
}

void SylinderSystem::applyForceTo(std::vector<double> forces, std::vector<int> gids) {
    const int size = gids.size();
    assert(size == 3 * forces.size());

    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        auto &sy = sylinderContainer[i];
        for (int j = 0; j < size; j++) {
            if (sy.gid == gids[j]) {
                sy.forceNonB[0] += forces[3 * i + 0];
                sy.forceNonB[1] += forces[3 * i + 1];
                sy.forceNonB[2] += forces[3 * i + 2];
            }
        }
    }
}

void SylinderSystem::setTreeSylinder() {
    // initialize tree
    // always keep tree max_glb_num_ptcl to be twice the global actual particle number.
    const int nGlobal = sylinderContainer.getNumberOfParticleGlobal();
    if (nGlobal > 1.5 * treeSylinderNumber || !treeSylinderNearPtr) {
        // a new larger tree
        treeSylinderNearPtr.reset();
        treeSylinderNearPtr = std::make_unique<TreeSylinderNear>();
        treeSylinderNearPtr->initialize(2 * nGlobal);
        treeSylinderNumber = nGlobal;
    }
}

void SylinderSystem::getOrient(Equatn &orient, const double px, const double py, const double pz, const int threadId) {
    Evec3 pvec;
    if (px < -1 || px > 1) {
        pvec[0] = 2 * rngPoolPtr->getU01(threadId) - 1;
    } else {
        pvec[0] = px;
    }
    if (py < -1 || py > 1) {
        pvec[1] = 2 * rngPoolPtr->getU01(threadId) - 1;
    } else {
        pvec[1] = py;
    }
    if (pz < -1 || pz > 1) {
        pvec[2] = 2 * rngPoolPtr->getU01(threadId) - 1;
    } else {
        pvec[2] = pz;
    }

    // px,py,pz all random, pick uniformly in orientation space
    if (px != pvec[0] && py != pvec[1] && pz != pvec[2]) {
        EquatnHelper::setUnitRandomEquatn(orient, rngPoolPtr->getU01(threadId), rngPoolPtr->getU01(threadId),
                                          rngPoolPtr->getU01(threadId));
        return;
    } else {
        orient = Equatn::FromTwoVectors(Evec3(0, 0, 1), pvec);
    }
}

void SylinderSystem::setInitialFromConfig() {
    // this function init all sylinders on rank 0
    if (runConfig.sylinderLengthSigma > 0) {
        rngPoolPtr->setLogNormalParameters(runConfig.sylinderLength, runConfig.sylinderLengthSigma);
    }

    if (commRcp->getRank() != 0) {
        sylinderContainer.setNumberOfParticleLocal(0);
    } else {
        const double boxEdge[3] = {runConfig.initBoxHigh[0] - runConfig.initBoxLow[0],
                                   runConfig.initBoxHigh[1] - runConfig.initBoxLow[1],
                                   runConfig.initBoxHigh[2] - runConfig.initBoxLow[2]};
        const double minBoxEdge = std::min(std::min(boxEdge[0], boxEdge[1]), boxEdge[2]);
        const double maxLength = minBoxEdge * 0.5;
        const double radius = runConfig.sylinderDiameter / 2;
        const int nSylinderLocal = runConfig.sylinderNumber;
        sylinderContainer.setNumberOfParticleLocal(nSylinderLocal);

#pragma omp parallel
        {
            const int threadId = omp_get_thread_num();
#pragma omp for
            for (int i = 0; i < nSylinderLocal; i++) {
                double length;
                if (runConfig.sylinderLengthSigma > 0) {
                    do { // generate random length
                        length = rngPoolPtr->getLN(threadId);
                    } while (length >= maxLength);
                } else {
                    length = runConfig.sylinderLength;
                }
                double pos[3];
                for (int k = 0; k < 3; k++) {
                    pos[k] = rngPoolPtr->getU01(threadId) * boxEdge[k] + runConfig.initBoxLow[k];
                }
                Equatn orientq;
                getOrient(orientq, runConfig.initOrient[0], runConfig.initOrient[1], runConfig.initOrient[2], threadId);
                double orientation[4];
                Emapq(orientation).coeffs() = orientq.coeffs();
                sylinderContainer[i] = Sylinder(i, radius, radius, length, length, pos, orientation);
                sylinderContainer[i].clear();
            }
        }
    }

    if (runConfig.initCircularX) {
        setInitialCircularCrossSection();
    }
}

void SylinderSystem::setInitialCircularCrossSection() {
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    double radiusCrossSec = 0;            // x, y, z, axis radius
    Evec3 centerCrossSec = Evec3::Zero(); // x, y, z, axis center.
    // x axis
    centerCrossSec = Evec3(0, (runConfig.initBoxHigh[1] - runConfig.initBoxLow[1]) * 0.5 + runConfig.initBoxLow[1],
                           (runConfig.initBoxHigh[2] - runConfig.initBoxLow[2]) * 0.5 + runConfig.initBoxLow[2]);
    radiusCrossSec = 0.5 * std::min(runConfig.initBoxHigh[2] - runConfig.initBoxLow[2],
                                    runConfig.initBoxHigh[1] - runConfig.initBoxLow[1]);
#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            double y = sylinderContainer[i].pos[1];
            double z = sylinderContainer[i].pos[2];
            // replace y,z with position in the circle
            getRandPointInCircle(radiusCrossSec, rngPoolPtr->getU01(threadId), rngPoolPtr->getU01(threadId), y, z);
            sylinderContainer[i].pos[1] = y + centerCrossSec[1];
            sylinderContainer[i].pos[2] = z + centerCrossSec[2];
        }
    }
}

void SylinderSystem::calcVolFrac() {
    // calc volume fraction of sphero cylinders
    // step 1, calc local total volume
    double volLocal = 0;
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
#pragma omp parallel for reduction(+ : volLocal)
    for (int i = 0; i < nLocal; i++) {
        auto &sy = sylinderContainer[i];
        volLocal += 3.1415926535 * (0.25 * sy.length * pow(sy.radius * 2, 2) + pow(sy.radius * 2, 3) / 6);
    }
    double volGlobal = 0;

    Teuchos::reduceAll(*commRcp, Teuchos::SumValueReductionOp<int, double>(), 1, &volLocal, &volGlobal);

    // step 2, reduce to root and compute total volume
    double boxVolume = (runConfig.simBoxHigh[0] - runConfig.simBoxLow[0]) *
                       (runConfig.simBoxHigh[1] - runConfig.simBoxLow[1]) *
                       (runConfig.simBoxHigh[2] - runConfig.simBoxLow[2]);
    spdlog::warn("Volume Sylinder = {:g}", volGlobal);
    spdlog::warn("Volume fraction = {:g}", volGlobal / boxVolume);
}

void SylinderSystem::setInitialFromFile(const std::string &filename) {
    spdlog::warn("Reading file " + filename);

    auto parseSylinder = [&](Sylinder &sy, const std::string &line) {
        std::stringstream liness(line);

        // data variables
        int gid;
        char type;
        double radius;
        double mx, my, mz;
        double px, py, pz;
        int group = -1;

        // Read in a single sylinder from the line 
        liness >> type 
               >> gid 
               >> radius 
               >> mx >> my >> mz 
               >> px >> py >> pz  
               >> group;

        // Set the sylinder object
        Emap3(sy.pos) = Evec3((mx + px), (my + py), (mz + pz)) * 0.5;
        sy.gid = gid;
        sy.group = group;
        sy.isImmovable = (type == 'S') ? true : false;
        sy.radius = radius;
        sy.radiusCollision = radius;
        sy.length = sqrt(pow(px - mx, 2) + pow(py - my, 2) + pow(pz - mz, 2));
        sy.lengthCollision = sy.length;
        if (sy.length > 1e-7) { // Setting direction for sylinder with small length is not useful
            Evec3 direction(px - mx, py - my, pz - mz);
            Emapq(sy.orientation) = Equatn::FromTwoVectors(Evec3(0, 0, 1), direction);
        } else {
            Emapq(sy.orientation) = Equatn::FromTwoVectors(Evec3(0, 0, 1), Evec3(0, 0, 1));
        }
    };

    if (commRcp->getRank() != 0) {
        sylinderContainer.setNumberOfParticleLocal(0);
    } else {
        std::ifstream myfile(filename);
        std::string line;
        std::getline(myfile, line); // read two header lines
        std::getline(myfile, line);

        std::deque<Sylinder> sylinderReadFromFile;
        while (std::getline(myfile, line)) {
            if (line[0] == 'C' || line[0] == 'S') {
                Sylinder sy;
                parseSylinder(sy, line);
                sylinderReadFromFile.push_back(sy);
            }
        }
        myfile.close();

        spdlog::debug("Sylinder number in file {} ", sylinderReadFromFile.size());

        // set on rank 0
        const int nRead = sylinderReadFromFile.size();
        sylinderContainer.setNumberOfParticleLocal(nRead);
#pragma omp parallel for
        for (int i = 0; i < nRead; i++) {
            sylinderContainer[i] = sylinderReadFromFile[i];
            sylinderContainer[i].clear();
        }
    }
}

void SylinderSystem::setLinkMapsFromFile(const std::string &filename) {
    spdlog::warn("Reading file " + filename);

    auto parseTwoLink = [&](Link &link, const std::string &line) {
        std::stringstream liness(line);
        char header;
        liness >> header >> link.prev >> link.next;
    };

    auto parseThreeLink = [&](const std::string &line) {
        std::stringstream liness(line);
        int gidI; 
        int gidJ; 
        int gidK;
        char header;
        liness >> header >> gidI >> gidJ >> gidK;
        assert(header == 'T');
        return std::make_tuple(gidI, gidJ, gidK);
    };

    std::ifstream myfile(filename);
    std::string line;
    std::getline(myfile, line); // read two header lines
    std::getline(myfile, line);

    // Make sure maps are clean
    pinLinkMap.clear();
    pinLinkReverseMap.clear();
    extendLinkMap.clear();
    extendLinkReverseMap.clear();
    bendLinkMap.clear();
    bendLinkReverseMap.clear();
    tribendLinkMap.clear();
    tribendLinkReverseMap.clear();

    while (std::getline(myfile, line)) {
        Link link;
        switch (line[0]) {
            case 'P': {
                parseTwoLink(link, line);
                pinLinkMap.emplace(link.prev, link.next);
                pinLinkReverseMap.emplace(link.next, link.prev);
                break;
            }
            case 'E': {
                parseTwoLink(link, line);
                extendLinkMap.emplace(link.prev, link.next);
                extendLinkReverseMap.emplace(link.next, link.prev);
                break;
            }
            case 'B': {
                parseTwoLink(link, line);
                bendLinkMap.emplace(link.prev, link.next);
                bendLinkReverseMap.emplace(link.next, link.prev);
                break;
            }
            case 'T': {
                auto [gidI, gidJ, gidK] = parseThreeLink(line);
                tribendLinkMap.emplace(gidI, std::make_pair(gidJ, gidK));
                tribendLinkReverseMap.emplace(std::make_pair(gidJ, gidK), gidI);
            break;
            }
            default:
                break;
        }
    }
    myfile.close();

    spdlog::debug("Pin link number in file {} ", pinLinkMap.size());
    spdlog::debug("Extend link number in file {} ", extendLinkMap.size());
    spdlog::debug("Bend link number in file {} ", bendLinkMap.size());
    spdlog::debug("Tribend link number in file {} ", tribendLinkMap.size());
}


void SylinderSystem::setInitialFromVTKFile(const std::string &pvtpFileName) {
    spdlog::warn("Reading file " + pvtpFileName);

    if (commRcp->getRank() != 0) {
        sylinderContainer.setNumberOfParticleLocal(0);
    } else {

        // Read the pvtp file and automatically merge the vtks files into a single polydata
        vtkSmartPointer<vtkXMLPPolyDataReader> reader = vtkSmartPointer<vtkXMLPPolyDataReader>::New();
        reader->SetFileName(pvtpFileName.c_str());
        reader->Update();

        // Extract the polydata (At this point, the polydata is unsorted)
        vtkSmartPointer<vtkPolyData> polydata = reader->GetOutput();
        // geometry data
        vtkSmartPointer<vtkPoints> posData = polydata->GetPoints();
        // Extract the point/cell data
        // int32 types
        vtkSmartPointer<vtkTypeInt32Array> gidData =
            vtkArrayDownCast<vtkTypeInt32Array>(polydata->GetCellData()->GetAbstractArray("gid"));
        vtkSmartPointer<vtkTypeInt32Array> groupData =
            vtkArrayDownCast<vtkTypeInt32Array>(polydata->GetCellData()->GetAbstractArray("group"));
        // unsigned char type
        vtkSmartPointer<vtkTypeUInt8Array> isImmovableData =
            vtkArrayDownCast<vtkTypeUInt8Array>(polydata->GetCellData()->GetAbstractArray("isImmovable"));
        // float/double types
        vtkSmartPointer<vtkDataArray> lengthData = polydata->GetCellData()->GetArray("length");
        vtkSmartPointer<vtkDataArray> lengthCollisionData = polydata->GetCellData()->GetArray("lengthCollision");
        vtkSmartPointer<vtkDataArray> radiusData = polydata->GetCellData()->GetArray("radius");
        vtkSmartPointer<vtkDataArray> radiusCollisionData = polydata->GetCellData()->GetArray("radiusCollision");
        vtkSmartPointer<vtkDataArray> znormData = polydata->GetCellData()->GetArray("znorm");
        vtkSmartPointer<vtkDataArray> velData = polydata->GetCellData()->GetArray("vel");
        vtkSmartPointer<vtkDataArray> omegaData = polydata->GetCellData()->GetArray("omega");

        const int sylinderNumberInFile = posData->GetNumberOfPoints() / 2; // two points per sylinder
        sylinderContainer.setNumberOfParticleLocal(sylinderNumberInFile);
        spdlog::debug("Sylinder number in file {} ", sylinderNumberInFile);

#pragma omp parallel for
        for (int i = 0; i < sylinderNumberInFile; i++) {
            auto &sy = sylinderContainer[i];
            double leftEndpointPos[3] = {0, 0, 0};
            double rightEndpointPos[3] = {0, 0, 0};
            posData->GetPoint(i * 2, leftEndpointPos);
            posData->GetPoint(i * 2 + 1, rightEndpointPos);

            Emap3(sy.pos) = (Emap3(leftEndpointPos) + Emap3(rightEndpointPos)) * 0.5;
            sy.gid = gidData->GetComponent(i, 0);
            sy.group = groupData->GetComponent(i, 0);
            sy.isImmovable = isImmovableData->GetTypedComponent(i, 0) > 0 ? true : false;
            sy.length = lengthData->GetComponent(i, 0);
            sy.lengthCollision = lengthCollisionData->GetComponent(i, 0);
            sy.radius = radiusData->GetComponent(i, 0);
            sy.radiusCollision = radiusCollisionData->GetComponent(i, 0);
            const Evec3 direction(znormData->GetComponent(i, 0), znormData->GetComponent(i, 1),
                                  znormData->GetComponent(i, 2));
            Emapq(sy.orientation) = Equatn::FromTwoVectors(Evec3(0, 0, 1), direction);
            sy.vel[0] = velData->GetComponent(i, 0);
            sy.vel[1] = velData->GetComponent(i, 1);
            sy.vel[2] = velData->GetComponent(i, 2);
            sy.omega[0] = omegaData->GetComponent(i, 0);
            sy.omega[1] = omegaData->GetComponent(i, 1);
            sy.omega[2] = omegaData->GetComponent(i, 2);
        }

        // sort the vector of Sylinders by gid ascending;
        // std::sort(sylinderReadFromFile.begin(), sylinderReadFromFile.end(),
        //           [](const Sylinder &t1, const Sylinder &t2) { return t1.gid < t2.gid; });
    }
    commRcp->barrier();
}

std::string SylinderSystem::getCurrentResultFolder() { return getResultFolderWithID(this->snapID); }

std::string SylinderSystem::getResultFolderWithID(int snapID_) {
    const int num = std::max(400 / commRcp->getSize(), 1); // limit max number of files per folder
    int k = snapID_ / num;
    int low = k * num, high = k * num + num - 1;
    std::string baseFolder =
        "./result/result" + std::to_string(low) + std::string("-") + std::to_string(high) + std::string("/");
    return baseFolder;
}

void SylinderSystem::writeAscii(const std::string &baseFolder) {
    // write a single ascii .dat file
    const int nGlobal = sylinderContainer.getNumberOfParticleGlobal();

    std::string name = baseFolder + std::string("SylinderAscii_") + std::to_string(snapID) + ".dat";
    SylinderAsciiHeader header;
    header.nparticle = nGlobal;
    header.time = stepCount * runConfig.dt;
    sylinderContainer.writeParticleAscii(name.c_str(), header);
    if (commRcp->getRank() == 0) {
        FILE *fptr = fopen(name.c_str(), "a");
        for (const auto &key_value : pinLinkMap) {
            fprintf(fptr, "P %d %d\n", key_value.first, key_value.second);
        }
        for (const auto &key_value : extendLinkMap) {
            fprintf(fptr, "E %d %d\n", key_value.first, key_value.second);
        }
        for (const auto &key_value : bendLinkMap) {
            fprintf(fptr, "B %d %d\n", key_value.first, key_value.second);
        }
        for (const auto &key_value : tribendLinkMap) {
            fprintf(fptr, "T %d %d %d\n", key_value.first, key_value.second.first, key_value.second.second);
        }
        fclose(fptr);
    }
    commRcp->barrier();
}

void SylinderSystem::writeTimeStepInfo(const std::string &baseFolder) {
    if (commRcp->getRank() == 0) {
        // write a single txt file containing timestep and most recent pvtp file names
        std::string name = baseFolder + std::string("../../TimeStepInfo.txt");
        std::string pvtpFileName = std::string("Sylinder_") + std::to_string(snapID) + std::string(".pvtp");

        FILE *restartFile = fopen(name.c_str(), "w");
        fprintf(restartFile, "%u\n", restartRngSeed);
        fprintf(restartFile, "%u\n", stepCount);
        fprintf(restartFile, "%u\n", snapID);
        fprintf(restartFile, "%s\n", pvtpFileName.c_str());
        fclose(restartFile);
    }
}

void SylinderSystem::writeVTK(const std::string &baseFolder) {
    const int rank = commRcp->getRank();
    const int size = commRcp->getSize();
    Sylinder::writeVTP<PS::ParticleSystem<Sylinder>>(sylinderContainer, sylinderContainer.getNumberOfParticleLocal(),
                                                     baseFolder, std::to_string(snapID), rank);
    conCollectorPtr->writeVTP(baseFolder, "", std::to_string(snapID), rank);
    if (rank == 0) {
        Sylinder::writePVTP(baseFolder, std::to_string(snapID), size); // write parallel head
        conCollectorPtr->writePVTP(baseFolder, "", std::to_string(snapID), size);
    }
}

void SylinderSystem::writeBox() {
    FILE *boxFile = fopen("./result/simBox.vtk", "w");
    fprintf(boxFile, "# vtk DataFile Version 3.0\n");
    fprintf(boxFile, "vtk file\n");
    fprintf(boxFile, "ASCII\n");
    fprintf(boxFile, "DATASET RECTILINEAR_GRID\n");
    fprintf(boxFile, "DIMENSIONS 2 2 2\n");
    fprintf(boxFile, "X_COORDINATES 2 float\n");
    fprintf(boxFile, "%g %g\n", runConfig.simBoxLow[0], runConfig.simBoxHigh[0]);
    fprintf(boxFile, "Y_COORDINATES 2 float\n");
    fprintf(boxFile, "%g %g\n", runConfig.simBoxLow[1], runConfig.simBoxHigh[1]);
    fprintf(boxFile, "Z_COORDINATES 2 float\n");
    fprintf(boxFile, "%g %g\n", runConfig.simBoxLow[2], runConfig.simBoxHigh[2]);
    fprintf(boxFile, "CELL_DATA 1\n");
    fprintf(boxFile, "POINT_DATA 8\n");
    fclose(boxFile);
}

void SylinderSystem::writeResult() {
    std::string baseFolder = getCurrentResultFolder();
    IOHelper::makeSubFolder(baseFolder);
    writeAscii(baseFolder);
    writeVTK(baseFolder);
    writeTimeStepInfo(baseFolder);
    snapID++;
}

void SylinderSystem::showOnScreenRank0() {
    if (commRcp->getRank() == 0) {
        printf("-----------SylinderSystem Settings-----------\n");
        runConfig.dump();
    }
}

void SylinderSystem::setDomainInfo() {
    const int pbcX = (runConfig.simBoxPBC[0] ? 1 : 0);
    const int pbcY = (runConfig.simBoxPBC[1] ? 1 : 0);
    const int pbcZ = (runConfig.simBoxPBC[2] ? 1 : 0);
    const int pbcFlag = 100 * pbcX + 10 * pbcY + pbcZ;

    switch (pbcFlag) {
    case 0:
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
        break;
    case 1:
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Z);
        break;
    case 10:
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_Y);
        break;
    case 100:
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_X);
        break;
    case 11:
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_YZ);
        break;
    case 101:
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XZ);
        break;
    case 110:
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XY);
        break;
    case 111:
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
        break;
    }

    PS::F64vec3 rootDomainLow;
    PS::F64vec3 rootDomainHigh;
    for (int k = 0; k < 3; k++) {
        rootDomainLow[k] = runConfig.simBoxLow[k];
        rootDomainHigh[k] = runConfig.simBoxHigh[k];
    }

    dinfo.setPosRootDomain(rootDomainLow, rootDomainHigh); // rootdomain must be specified after PBC
}

void SylinderSystem::decomposeDomain() {
    applyBoxBC();
    dinfo.decomposeDomainAll(sylinderContainer);
}

void SylinderSystem::exchangeSylinder() {
    sylinderContainer.exchangeParticle(dinfo);
    updateSylinderRank();
}

void SylinderSystem::calcMobMatrix() {
    // diagonal hydro mobility operator
    // 3*3 block for translational + 3*3 block for rotational.
    // 3 nnz per row, 18 nnz per tubule

    const double Pi = 3.14159265358979323846;
    const double mu = runConfig.viscosity;

    const int nLocal = sylinderMapRcp->getNodeNumElements();
    TEUCHOS_ASSERT(nLocal == sylinderContainer.getNumberOfParticleLocal());
    const int localSize = nLocal * 6; // local row number

    Kokkos::View<size_t *> rowPointers("rowPointers", localSize + 1);
    rowPointers[0] = 0;
    for (int i = 1; i <= localSize; i++) {
        rowPointers[i] = rowPointers[i - 1] + 3;
    }
    Kokkos::View<int *> columnIndices("columnIndices", rowPointers[localSize]);
    Kokkos::View<double *> values("values", rowPointers[localSize]);

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        const auto &sy = sylinderContainer[i];

        // calculate the Mob Trans and MobRot
        Emat3 MobTrans; //            double MobTrans[3][3];
        Emat3 MobRot;   //            double MobRot[3][3];
        Emat3 qq;
        Emat3 Imqq;
        Evec3 q = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        qq = q * q.transpose();
        Imqq = Emat3::Identity() - qq;

        double dragPara = 0;
        double dragPerp = 0;
        double dragRot = 0;
        sy.calcDragCoeff(mu, dragPara, dragPerp, dragRot);
        const double dragParaInv = sy.isImmovable ? 0.0 : 1 / dragPara;
        const double dragPerpInv = sy.isImmovable ? 0.0 : 1 / dragPerp;
        const double dragRotInv = sy.isImmovable ? 0.0 : 1 / dragRot;

        MobTrans = dragParaInv * qq + dragPerpInv * Imqq;
        MobRot = dragRotInv * qq + dragRotInv * Imqq; // = dragRotInv * Identity
        // MobRot regularized to remove null space.
        // here it becomes identity matrix,
        // no effect on geometric constraints
        // no problem for axissymetric slender body.
        // this simplifies the rotational Brownian calculations.

        // column index is local index
        columnIndices[18 * i] = 6 * i; // line 1 of Mob Trans
        columnIndices[18 * i + 1] = 6 * i + 1;
        columnIndices[18 * i + 2] = 6 * i + 2;
        columnIndices[18 * i + 3] = 6 * i; // line 2 of Mob Trans
        columnIndices[18 * i + 4] = 6 * i + 1;
        columnIndices[18 * i + 5] = 6 * i + 2;
        columnIndices[18 * i + 6] = 6 * i; // line 3 of Mob Trans
        columnIndices[18 * i + 7] = 6 * i + 1;
        columnIndices[18 * i + 8] = 6 * i + 2;
        columnIndices[18 * i + 9] = 6 * i + 3; // line 1 of Mob Rot
        columnIndices[18 * i + 10] = 6 * i + 4;
        columnIndices[18 * i + 11] = 6 * i + 5;
        columnIndices[18 * i + 12] = 6 * i + 3; // line 2 of Mob Rot
        columnIndices[18 * i + 13] = 6 * i + 4;
        columnIndices[18 * i + 14] = 6 * i + 5;
        columnIndices[18 * i + 15] = 6 * i + 3; // line 3 of Mob Rot
        columnIndices[18 * i + 16] = 6 * i + 4;
        columnIndices[18 * i + 17] = 6 * i + 5;

        values[18 * i] = MobTrans(0, 0); // line 1 of Mob Trans
        values[18 * i + 1] = MobTrans(0, 1);
        values[18 * i + 2] = MobTrans(0, 2);
        values[18 * i + 3] = MobTrans(1, 0); // line 2 of Mob Trans
        values[18 * i + 4] = MobTrans(1, 1);
        values[18 * i + 5] = MobTrans(1, 2);
        values[18 * i + 6] = MobTrans(2, 0); // line 3 of Mob Trans
        values[18 * i + 7] = MobTrans(2, 1);
        values[18 * i + 8] = MobTrans(2, 2);
        values[18 * i + 9] = MobRot(0, 0); // line 1 of Mob Rot
        values[18 * i + 10] = MobRot(0, 1);
        values[18 * i + 11] = MobRot(0, 2);
        values[18 * i + 12] = MobRot(1, 0); // line 2 of Mob Rot
        values[18 * i + 13] = MobRot(1, 1);
        values[18 * i + 14] = MobRot(1, 2);
        values[18 * i + 15] = MobRot(2, 0); // line 3 of Mob Rot
        values[18 * i + 16] = MobRot(2, 1);
        values[18 * i + 17] = MobRot(2, 2);
    }

    // mobMat is block-diagonal, so domainMap=rangeMap
    mobilityMatrixRcp =
        Teuchos::rcp(new TCMAT(sylinderMobilityMapRcp, sylinderMobilityMapRcp, rowPointers, columnIndices, values));
    mobilityMatrixRcp->fillComplete(sylinderMobilityMapRcp, sylinderMobilityMapRcp); // domainMap, rangeMap

    spdlog::debug("MobMat Constructed " + mobilityMatrixRcp->description());
}

void SylinderSystem::calcMobOperator() {
    calcMobMatrix();
    mobilityOperatorRcp = mobilityMatrixRcp;
}

void SylinderSystem::calcVelocityNonCon() {
    // velocityNonCon = velocityBrown + velocityPartNonBrown + mobility * forcePartNonBrown
    // if monolayer, set velBrownZ =0, velPartNonBrownZ =0, forcePartNonBrownZ =0
    velocityNonConRcp = Teuchos::rcp<TV>(new TV(sylinderMobilityMapRcp, true)); // allocate and zero out
    auto velNCPtr = velocityNonConRcp->getLocalView<Kokkos::HostSpace>();

    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    TEUCHOS_ASSERT(nLocal * 6 == velocityNonConRcp->getLocalLength());

    if (!forcePartNonBrownRcp.is_null()) {
        // apply mobility
        TEUCHOS_ASSERT(!mobilityOperatorRcp.is_null());
        mobilityOperatorRcp->apply(*forcePartNonBrownRcp, *velocityNonConRcp);
        if (runConfig.monolayer) {
#pragma omp parallel for
            for (int i = 0; i < nLocal; i++) {
                velNCPtr(6 * i + 2, 0) = 0; // vz
                velNCPtr(6 * i + 3, 0) = 0; // omegax
                velNCPtr(6 * i + 4, 0) = 0; // omegay
            }
        }
        // write back to Sylinder members
        auto forcePtr = forcePartNonBrownRcp->getLocalView<Kokkos::HostSpace>();
#pragma omp parallel for
        for (int i = 0; i < nLocal; i++) {
            auto &sy = sylinderContainer[i];
            // torque
            sy.forceNonB[0] = forcePtr(6 * i + 0, 0);
            sy.forceNonB[1] = forcePtr(6 * i + 1, 0);
            sy.forceNonB[2] = forcePtr(6 * i + 2, 0);
            sy.torqueNonB[0] = forcePtr(6 * i + 3, 0);
            sy.torqueNonB[1] = forcePtr(6 * i + 4, 0);
            sy.torqueNonB[2] = forcePtr(6 * i + 5, 0);
        }
    }

    if (!velocityPartNonBrownRcp.is_null()) {
        if (runConfig.monolayer) {
            auto velNBPtr = velocityPartNonBrownRcp->getLocalView<Kokkos::HostSpace>();
#pragma omp parallel for
            for (int i = 0; i < nLocal; i++) {
                velNBPtr(6 * i + 2, 0) = 0; // vz
                velNBPtr(6 * i + 3, 0) = 0; // omegax
                velNBPtr(6 * i + 4, 0) = 0; // omegay
            }
        }
        velocityNonConRcp->update(1.0, *velocityPartNonBrownRcp, 1.0);
    }

    // write back total non Brownian velocity
    // combine and sync the velNonB set in by setForceNonBrown() and setVelocityNonBrown()
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        auto &sy = sylinderContainer[i];
        // velocity
        sy.velNonB[0] = velNCPtr(6 * i + 0, 0);
        sy.velNonB[1] = velNCPtr(6 * i + 1, 0);
        sy.velNonB[2] = velNCPtr(6 * i + 2, 0);
        sy.omegaNonB[0] = velNCPtr(6 * i + 3, 0);
        sy.omegaNonB[1] = velNCPtr(6 * i + 4, 0);
        sy.omegaNonB[2] = velNCPtr(6 * i + 5, 0);
    }

    // add Brownian motion
    if (!velocityBrownRcp.is_null()) {
        if (runConfig.monolayer) {
            auto velBPtr = velocityBrownRcp->getLocalView<Kokkos::HostSpace>();
#pragma omp parallel for
            for (int i = 0; i < nLocal; i++) {
                velBPtr(6 * i + 2, 0) = 0; // vz
                velBPtr(6 * i + 3, 0) = 0; // omegax
                velBPtr(6 * i + 4, 0) = 0; // omegay
            }
        }
        velocityNonConRcp->update(1.0, *velocityBrownRcp, 1.0);
    }
}

void SylinderSystem::sumForceVelocity() {
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        auto &sy = sylinderContainer[i];
        for (int k = 0; k < 3; k++) {
            sy.vel[k] = sy.velNonB[k] + sy.velBrown[k] + sy.velCol[k] + sy.velBi[k];
            sy.omega[k] = sy.omegaNonB[k] + sy.omegaBrown[k] + sy.omegaCol[k] + sy.omegaBi[k];
            sy.force[k] = sy.forceNonB[k] + sy.forceCol[k] + sy.forceBi[k];
            sy.torque[k] = sy.torqueNonB[k] + sy.torqueCol[k] + sy.torqueBi[k];
        }
    }
}

void SylinderSystem::stepEuler() {
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    const double dt = runConfig.dt;

    if (!runConfig.sylinderFixed) {
#pragma omp parallel for
        for (int i = 0; i < nLocal; i++) {
            auto &sy = sylinderContainer[i];
            sy.stepEuler(dt);
        }
    }
}

void SylinderSystem::resolveConstraints() {

    Teuchos::RCP<Teuchos::Time> collectColTimer =
        Teuchos::TimeMonitor::getNewCounter("SylinderSystem::CollectCollision");
    Teuchos::RCP<Teuchos::Time> collectLinkTimer = Teuchos::TimeMonitor::getNewCounter("SylinderSystem::CollectLink");

    spdlog::debug("start collect collisions");
    {
        Teuchos::TimeMonitor mon(*collectColTimer);
        collectPairCollision();
        collectBoundaryCollision();
    }

    spdlog::debug("start collect links");
    {
        Teuchos::TimeMonitor mon(*collectLinkTimer);
        collectPinLinkBilateral();
        collectExtendLinkBilateral();
        collectBendLinkBilateral();
        collectTriBendLinkBilateral();
    }

    // solve collision
    // positive buffer value means collision radius is effectively smaller
    // i.e., less likely to collide
    Teuchos::RCP<Teuchos::Time> solveTimer = Teuchos::TimeMonitor::getNewCounter("SylinderSystem::SolveConstraints");
    {
        Teuchos::TimeMonitor mon(*solveTimer);
        const double buffer = 0;
        spdlog::debug("constraint solver setup");
        conSolverPtr->setup(*conCollectorPtr, mobilityOperatorRcp, velocityNonConRcp, runConfig.dt);
        spdlog::debug("setControl");
        conSolverPtr->setControlParams(runConfig.conResTol, runConfig.conMaxIte, runConfig.conSolverChoice);
        spdlog::debug("solveConstraints");
        conSolverPtr->solveConstraints();
        spdlog::debug("writebackGamma");
        conSolverPtr->writebackGamma();
    }

    saveForceVelocityConstraints();
}

void SylinderSystem::updateSylinderMap() {
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    // setup the new sylinderMap
    sylinderMapRcp = getTMAPFromLocalSize(nLocal, commRcp);
    sylinderMobilityMapRcp = getTMAPFromLocalSize(nLocal * 6, commRcp);

    // setup the globalIndex
    int globalIndexBase = sylinderMapRcp->getMinGlobalIndex(); // this is a contiguous map
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        sylinderContainer[i].globalIndex = i + globalIndexBase;
    }
}

bool SylinderSystem::getIfWriteResultCurrentStep() {
    return (stepCount % static_cast<int>(runConfig.timeSnap / runConfig.dt) == 0);
}

void SylinderSystem::prepareStep() {
    spdlog::warn("CurrentStep {}", stepCount);
    applyBoxBC();

    if (stepCount % 50 == 0) {
        decomposeDomain();
    }

    exchangeSylinder();

    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        auto &sy = sylinderContainer[i];
        sy.clear();
        sy.radiusCollision = sylinderContainer[i].radius * runConfig.sylinderDiameterColRatio;
        sy.lengthCollision = sylinderContainer[i].length * runConfig.sylinderLengthColRatio;
        sy.rank = commRcp->getRank();
        sy.colBuf = runConfig.sylinderColBuf;
    }

    if (runConfig.monolayer) {
        const double monoZ = (runConfig.simBoxHigh[2] + runConfig.simBoxLow[2]) / 2;
#pragma omp parallel for
        for (int i = 0; i < nLocal; i++) {
            auto &sy = sylinderContainer[i];
            sy.pos[2] = monoZ;
            Evec3 drt = Emapq(sy.orientation) * Evec3(0, 0, 1);
            drt[2] = 0;
            drt.normalize();
            Emapq(sy.orientation).setFromTwoVectors(Evec3(0, 0, 1), drt);
        }
    }

    updateSylinderMap();

    buildSylinderNearDataDirectory();

    calcMobOperator();

    conCollectorPtr->clear();

    forcePartNonBrownRcp.reset();
    velocityPartNonBrownRcp.reset();
    velocityNonBrownRcp.reset();
    velocityBrownRcp.reset();
}

void SylinderSystem::setForceNonBrown(const std::vector<double> &forceNonBrown) {
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    TEUCHOS_ASSERT(forceNonBrown.size() == 6 * nLocal);
    TEUCHOS_ASSERT(sylinderMobilityMapRcp->getNodeNumElements() == 6 * nLocal);
    forcePartNonBrownRcp = getTVFromVector(forceNonBrown, commRcp);
}

void SylinderSystem::setVelocityNonBrown(const std::vector<double> &velNonBrown) {
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    TEUCHOS_ASSERT(velNonBrown.size() == 6 * nLocal);
    TEUCHOS_ASSERT(sylinderMobilityMapRcp->getNodeNumElements() == 6 * nLocal);
    velocityPartNonBrownRcp = getTVFromVector(velNonBrown, commRcp);
}

void SylinderSystem::runStep(bool count_flag) {

    if (runConfig.KBT > 0) {
        calcVelocityBrown();
    }

    calcVelocityNonCon();
    resolveConstraints();

    sumForceVelocity();

    if (getIfWriteResultCurrentStep() && count_flag) {
        // write result before moving. guarantee data written is consistent to geometry
        writeResult();
    }

    stepEuler();

    if (count_flag)
        stepCount++;
}

void SylinderSystem::saveForceVelocityConstraints() {
    // save results
    forceUniRcp = conSolverPtr->getForceUni();
    velocityUniRcp = conSolverPtr->getVelocityUni();
    forceBiRcp = conSolverPtr->getForceBi();
    velocityBiRcp = conSolverPtr->getVelocityBi();

    auto velUniPtr = velocityUniRcp->getLocalView<Kokkos::HostSpace>();
    auto velBiPtr = velocityBiRcp->getLocalView<Kokkos::HostSpace>();
    auto forceUniPtr = forceUniRcp->getLocalView<Kokkos::HostSpace>();
    auto forceBiPtr = forceBiRcp->getLocalView<Kokkos::HostSpace>();

    const int sylinderLocalNumber = sylinderContainer.getNumberOfParticleLocal();
    TEUCHOS_ASSERT(velUniPtr.dimension_0() == sylinderLocalNumber * 6);
    TEUCHOS_ASSERT(velUniPtr.dimension_1() == 1);
    TEUCHOS_ASSERT(velBiPtr.dimension_0() == sylinderLocalNumber * 6);
    TEUCHOS_ASSERT(velBiPtr.dimension_1() == 1);

#pragma omp parallel for
    for (int i = 0; i < sylinderLocalNumber; i++) {
        auto &sy = sylinderContainer[i];
        sy.velCol[0] = velUniPtr(6 * i + 0, 0);
        sy.velCol[1] = velUniPtr(6 * i + 1, 0);
        sy.velCol[2] = velUniPtr(6 * i + 2, 0);
        sy.omegaCol[0] = velUniPtr(6 * i + 3, 0);
        sy.omegaCol[1] = velUniPtr(6 * i + 4, 0);
        sy.omegaCol[2] = velUniPtr(6 * i + 5, 0);
        sy.velBi[0] = velBiPtr(6 * i + 0, 0);
        sy.velBi[1] = velBiPtr(6 * i + 1, 0);
        sy.velBi[2] = velBiPtr(6 * i + 2, 0);
        sy.omegaBi[0] = velBiPtr(6 * i + 3, 0);
        sy.omegaBi[1] = velBiPtr(6 * i + 4, 0);
        sy.omegaBi[2] = velBiPtr(6 * i + 5, 0);

        sy.forceCol[0] = forceUniPtr(6 * i + 0, 0);
        sy.forceCol[1] = forceUniPtr(6 * i + 1, 0);
        sy.forceCol[2] = forceUniPtr(6 * i + 2, 0);
        sy.torqueCol[0] = forceUniPtr(6 * i + 3, 0);
        sy.torqueCol[1] = forceUniPtr(6 * i + 4, 0);
        sy.torqueCol[2] = forceUniPtr(6 * i + 5, 0);
        sy.forceBi[0] = forceBiPtr(6 * i + 0, 0);
        sy.forceBi[1] = forceBiPtr(6 * i + 1, 0);
        sy.forceBi[2] = forceBiPtr(6 * i + 2, 0);
        sy.torqueBi[0] = forceBiPtr(6 * i + 3, 0);
        sy.torqueBi[1] = forceBiPtr(6 * i + 4, 0);
        sy.torqueBi[2] = forceBiPtr(6 * i + 5, 0);
    }
}

void SylinderSystem::calcVelocityBrown() {
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    const double Pi = 3.1415926535897932384626433;
    const double mu = runConfig.viscosity;
    const double dt = runConfig.dt;
    const double delta = dt * 0.1; // a small parameter used in RFD algorithm
    const double kBT = runConfig.KBT;
    const double kBTfactor = sqrt(2 * kBT / dt);

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            auto &sy = sylinderContainer[i];
            // constants
            double dragPara = 0;
            double dragPerp = 0;
            double dragRot = 0;
            sy.calcDragCoeff(mu, dragPara, dragPerp, dragRot);
            const double dragParaInv = sy.isImmovable ? 0.0 : 1 / dragPara;
            const double dragPerpInv = sy.isImmovable ? 0.0 : 1 / dragPerp;
            const double dragRotInv = sy.isImmovable ? 0.0 : 1 / dragRot;

            // convert FDPS vec3 to Evec3
            Evec3 direction = Emapq(sy.orientation) * Evec3(0, 0, 1);

            // RFD from Delong, JCP, 2015
            // slender fiber has 0 rot drag, regularize with identity rot mobility
            // trans mobility is this
            Evec3 q = direction;
            Emat3 Nmat = (dragParaInv - dragPerpInv) * (q * q.transpose()) + (dragPerpInv)*Emat3::Identity();
            Emat3 Nmatsqrt = Nmat.llt().matrixL();

            // velocity
            Evec3 Wrot(rngPoolPtr->getN01(threadId), rngPoolPtr->getN01(threadId), rngPoolPtr->getN01(threadId));
            Evec3 Wpos(rngPoolPtr->getN01(threadId), rngPoolPtr->getN01(threadId), rngPoolPtr->getN01(threadId));
            Evec3 Wrfdrot(rngPoolPtr->getN01(threadId), rngPoolPtr->getN01(threadId), rngPoolPtr->getN01(threadId));
            Evec3 Wrfdpos(rngPoolPtr->getN01(threadId), rngPoolPtr->getN01(threadId), rngPoolPtr->getN01(threadId));

            Equatn orientRFD = Emapq(sy.orientation);
            EquatnHelper::rotateEquatn(orientRFD, Wrfdrot, delta);
            q = orientRFD * Evec3(0, 0, 1);
            Emat3 Nmatrfd = (dragParaInv - dragPerpInv) * (q * q.transpose()) + (dragPerpInv)*Emat3::Identity();

            Evec3 vel = kBTfactor * (Nmatsqrt * Wpos);           // Gaussian noise
            vel += (kBT / delta) * ((Nmatrfd - Nmat) * Wrfdpos); // rfd drift. seems no effect in this case
            Evec3 omega = sqrt(dragRotInv) * kBTfactor * Wrot;   // regularized identity rotation drag

            Emap3(sy.velBrown) = vel;
            Emap3(sy.omegaBrown) = omega;
        }
    }

    velocityBrownRcp = Teuchos::rcp<TV>(new TV(sylinderMobilityMapRcp, true));
    auto velocityPtr = velocityBrownRcp->getLocalView<Kokkos::HostSpace>();
    velocityBrownRcp->modify<Kokkos::HostSpace>();

    TEUCHOS_ASSERT(velocityPtr.dimension_0() == nLocal * 6);
    TEUCHOS_ASSERT(velocityPtr.dimension_1() == 1);

#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        const auto &sy = sylinderContainer[i];
        velocityPtr(6 * i, 0) = sy.velBrown[0];
        velocityPtr(6 * i + 1, 0) = sy.velBrown[1];
        velocityPtr(6 * i + 2, 0) = sy.velBrown[2];
        velocityPtr(6 * i + 3, 0) = sy.omegaBrown[0];
        velocityPtr(6 * i + 4, 0) = sy.omegaBrown[1];
        velocityPtr(6 * i + 5, 0) = sy.omegaBrown[2];
    }
}

void SylinderSystem::collectBoundaryCollision() {
    auto collisionPoolPtr = conCollectorPtr->constraintPoolPtr; // shared_ptr
    const int nThreads = collisionPoolPtr->size();
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();

    // process collisions with all boundaries
    for (const auto &bPtr : runConfig.boundaryPtr) {
#pragma omp parallel num_threads(nThreads)
        {
            const int threadId = omp_get_thread_num();
            auto &que = (*collisionPoolPtr)[threadId];
#pragma omp for
            for (int i = 0; i < nLocal; i++) {
                const auto &sy = sylinderContainer[i];
                const Evec3 center = ECmap3(sy.pos);

                // check one point
                auto checkEnd = [&](const Evec3 &Query, const double radius) {
                    double Proj[3], delta[3];
                    bPtr->project(Query.data(), Proj, delta);
                    // if (!bPtr->check(Query.data(), Proj, delta)) {
                    //     printf("boundary projection error\n");
                    // }
                    // if inside boundary, delta = Q-Proj
                    // if outside boundary, delta = Proj-Q
                    double deltanorm = Emap3(delta).norm();
                    Evec3 norm = Emap3(delta) * (1 / deltanorm);
                    Evec3 posI = Query - center;

                    if ((Query - ECmap3(Proj)).dot(ECmap3(delta)) < 0) { // outside boundary
                        que.emplace_back(-deltanorm - radius, 0, sy.gid, sy.gid, sy.globalIndex, sy.globalIndex,
                                         norm.data(), norm.data(), posI.data(), posI.data(), Query.data(), Proj, true,
                                         false, 0.0);
                    } else if (deltanorm <
                               (1 + runConfig.sylinderColBuf * 2) * sy.radiusCollision) { // inside boundary but close
                        que.emplace_back(deltanorm - radius, 0, sy.gid, sy.gid, sy.globalIndex, sy.globalIndex,
                                         norm.data(), norm.data(), posI.data(), posI.data(), Query.data(), Proj, true,
                                         false, 0.0);
                    }
                };

                if (sy.isSphere(true)) {
                    double radius = sy.lengthCollision * 0.5 + sy.radiusCollision;
                    checkEnd(center, radius);
                } else {
                    const Equatn orientation = ECmapq(sy.orientation);
                    const Evec3 direction = orientation * Evec3(0, 0, 1);
                    const double length = sy.lengthCollision;
                    const Evec3 Qm = center - direction * (length * 0.5);
                    const Evec3 Qp = center + direction * (length * 0.5);
                    checkEnd(Qm, sy.radiusCollision);
                    checkEnd(Qp, sy.radiusCollision);
                }
            }
        }
    }
    return;
}

void SylinderSystem::collectPairCollision() {

    CalcSylinderNearForce calcColFtr(conCollectorPtr->constraintPoolPtr, pinLinkMap);

    TEUCHOS_ASSERT(treeSylinderNearPtr);
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    setTreeSylinder();
    treeSylinderNearPtr->calcForceAll(calcColFtr, sylinderContainer, dinfo);
}

std::pair<int, int> SylinderSystem::getMaxGid() {
    int maxGidLocal = 0;
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    for (int i = 0; i < nLocal; i++) {
        maxGidLocal = std::max(maxGidLocal, sylinderContainer[i].gid);
    }

    int maxGidGlobal = maxGidLocal;
    Teuchos::reduceAll(*commRcp, Teuchos::MaxValueReductionOp<int, int>(), 1, &maxGidLocal, &maxGidGlobal);
    spdlog::warn("rank: {}, maxGidLocal: {}, maxGidGlobal {}", commRcp->getRank(), maxGidLocal, maxGidGlobal);

    return std::pair<int, int>(maxGidLocal, maxGidGlobal);
}

void SylinderSystem::calcBoundingBox(double localLow[3], double localHigh[3], double globalLow[3],
                                     double globalHigh[3]) {
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    double lx, ly, lz;
    lx = ly = lz = std::numeric_limits<double>::max();
    double hx, hy, hz;
    hx = hy = hz = std::numeric_limits<double>::min();

    for (int i = 0; i < nLocal; i++) {
        const auto &sy = sylinderContainer[i];
        const Evec3 direction = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        Evec3 pm = ECmap3(sy.pos) - (sy.length * 0.5) * direction;
        Evec3 pp = ECmap3(sy.pos) + (sy.length * 0.5) * direction;
        lx = std::min(std::min(lx, pm[0]), pp[0]);
        ly = std::min(std::min(ly, pm[1]), pp[1]);
        lz = std::min(std::min(lz, pm[2]), pp[2]);
        hx = std::max(std::max(hx, pm[0]), pp[0]);
        hy = std::max(std::max(hy, pm[1]), pp[1]);
        hz = std::max(std::max(hz, pm[2]), pp[2]);
    }

    localLow[0] = lx;
    localLow[1] = ly;
    localLow[2] = lz;
    localHigh[0] = hx;
    localHigh[1] = hy;
    localHigh[2] = hz;

    for (int k = 0; k < 3; k++) {
        globalLow[k] = localLow[k];
        globalHigh[k] = localHigh[k];
    }

    Teuchos::reduceAll(*commRcp, Teuchos::MinValueReductionOp<int, double>(), 3, localLow, globalLow);
    Teuchos::reduceAll(*commRcp, Teuchos::MaxValueReductionOp<int, double>(), 3, localHigh, globalHigh);

    return;
}

void SylinderSystem::updateSylinderRank() {
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    const int rank = commRcp->getRank();
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        sylinderContainer[i].rank = rank;
    }
}

void SylinderSystem::applyBoxBC() { sylinderContainer.adjustPositionIntoRootDomain(dinfo); }

void SylinderSystem::calcConStress() {
    if (runConfig.logLevel > spdlog::level::info)
        return;

    Emat3 sumBiStress = Emat3::Zero();
    Emat3 sumUniStress = Emat3::Zero();
    conCollectorPtr->sumLocalConstraintStress(sumUniStress, sumBiStress, false);

    // scale to nkBT
    const double scaleFactor = 1 / (sylinderMapRcp->getGlobalNumElements() * runConfig.KBT);
    sumBiStress *= scaleFactor;
    sumUniStress *= scaleFactor;
    // mpi reduction
    double uniStressLocal[9];
    double biStressLocal[9];
    double uniStressGlobal[9];
    double biStressGlobal[9];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            uniStressLocal[i * 3 + j] = sumUniStress(i, j);
            uniStressGlobal[i * 3 + j] = 0;
            biStressLocal[i * 3 + j] = sumBiStress(i, j);
            biStressGlobal[i * 3 + j] = 0;
        }
    }

    Teuchos::reduceAll(*commRcp, Teuchos::SumValueReductionOp<int, double>(), 9, uniStressLocal, uniStressGlobal);
    Teuchos::reduceAll(*commRcp, Teuchos::SumValueReductionOp<int, double>(), 9, biStressLocal, biStressGlobal);

    spdlog::info("RECORD: ColXF,{:g},{:g},{:g},{:g},{:g},{:g},{:g},{:g},{:g}", //
                 uniStressGlobal[0], uniStressGlobal[1], uniStressGlobal[2],   //
                 uniStressGlobal[3], uniStressGlobal[4], uniStressGlobal[5],   //
                 uniStressGlobal[6], uniStressGlobal[7], uniStressGlobal[8]);
    spdlog::info("RECORD: BiXF,{:g},{:g},{:g},{:g},{:g},{:g},{:g},{:g},{:g}", //
                 biStressGlobal[0], biStressGlobal[1], biStressGlobal[2],     //
                 biStressGlobal[3], biStressGlobal[4], biStressGlobal[5],     //
                 biStressGlobal[6], biStressGlobal[7], biStressGlobal[8]);
}

void SylinderSystem::calcOrderParameter() {
    if (runConfig.logLevel > spdlog::level::info)
        return;

    double px = 0, py = 0, pz = 0;    // pvec
    double Qxx = 0, Qxy = 0, Qxz = 0; // Qtensor
    double Qyx = 0, Qyy = 0, Qyz = 0; // Qtensor
    double Qzx = 0, Qzy = 0, Qzz = 0; // Qtensor

    const int nLocal = sylinderContainer.getNumberOfParticleLocal();

#pragma omp parallel for reduction(+ : px, py, pz, Qxx, Qxy, Qxz, Qyx, Qyy, Qyz, Qzx, Qzy, Qzz)
    for (int i = 0; i < nLocal; i++) {
        const auto &sy = sylinderContainer[i];
        const Evec3 direction = ECmapq(sy.orientation) * Evec3(0, 0, 1);
        px += direction.x();
        py += direction.y();
        pz += direction.z();
        const Emat3 Q = direction * direction.transpose() - Emat3::Identity() * (1 / 3.0);
        Qxx += Q(0, 0);
        Qxy += Q(0, 1);
        Qxz += Q(0, 2);
        Qyx += Q(1, 0);
        Qyy += Q(1, 1);
        Qyz += Q(1, 2);
        Qzx += Q(2, 0);
        Qzy += Q(2, 1);
        Qzz += Q(2, 2);
    }

    // global average
    const int nGlobal = sylinderContainer.getNumberOfParticleGlobal();
    double pQ[12] = {px, py, pz, Qxx, Qxy, Qxz, Qyx, Qyy, Qyz, Qzx, Qzy, Qzz};
    MPI_Allreduce(MPI_IN_PLACE, pQ, 12, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < 12; i++) {
        pQ[i] *= (1.0 / nGlobal);
    }

    spdlog::info("RECORD: Order P,{:g},{:g},{:g},Q,{:g},{:g},{:g},{:g},{:g},{:g},{:g},{:g},{:g}", //
                 pQ[0], pQ[1], pQ[2],                                                             // pvec
                 pQ[3], pQ[4], pQ[5],                                                             // Qtensor
                 pQ[6], pQ[7], pQ[8],                                                             // Qtensor
                 pQ[9], pQ[10], pQ[11]                                                            // Qtensor
    );
}

std::vector<int> SylinderSystem::addNewSylinder(const std::vector<Sylinder> &newSylinder) {
    // assign unique new gid for sylinders on all ranks
    std::pair<int, int> maxGid = getMaxGid();
    const int maxGidLocal = maxGid.first;
    const int maxGidGlobal = maxGid.second;
    const int newCountLocal = newSylinder.size();

    // collect number of ids from all ranks to rank0
    std::vector<int> newCount(commRcp->getSize(), 0);
    std::vector<int> displ(commRcp->getSize() + 1, 0);
    MPI_Gather(&newCountLocal, 1, MPI_INT, newCount.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> newGid;
    if (commRcp->getRank() == 0) {
        // generate random gid on rank 0
        std::partial_sum(newCount.cbegin(), newCount.cend(), displ.begin() + 1);
        const int newCountGlobal = displ.back();
        newGid.resize(newCountGlobal, 0);
        std::iota(newGid.begin(), newGid.end(), maxGidGlobal + 1);
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(newGid.begin(), newGid.end(), g);
    } else {
        newGid.resize(newCountLocal, 0);
    }

    // scatter from rank 0 to every rank
    std::vector<int> newGidRecv(newCountLocal, 0);
    MPI_Scatterv(newGid.data(), newCount.data(), displ.data(), MPI_INT, //
                 newGidRecv.data(), newCountLocal, MPI_INT, 0, MPI_COMM_WORLD);

    // set new gid
    for (int i = 0; i < newCountLocal; i++) {
        Sylinder sy = newSylinder[i];
        sy.gid = newGidRecv[i];
        sylinderContainer.addOneParticle(sy);
    }

    return newGidRecv;
}

void SylinderSystem::addNewEndLink(const std::vector<Link> &newEndLink) {
    // synchronize newLink to all mpi ranks
    const int newCountLocal = newEndLink.size();
    std::vector<int> newCount(commRcp->getSize(), 0);
    MPI_Allgather(&newCountLocal, 1, MPI_INT, newCount.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<int> displ(commRcp->getSize() + 1, 0);
    std::partial_sum(newCount.cbegin(), newCount.cend(), displ.begin() + 1);
    std::vector<Link> newLinkRecv(displ.back());
    MPI_Allgatherv(newEndLink.data(), newCountLocal, createMPIStructType<Link>(), newLinkRecv.data(), newCount.data(),
                   displ.data(), createMPIStructType<Link>(), MPI_COMM_WORLD);

    // put newLinks into the map, same op on all mpi ranks
    for (const auto &ll : newLinkRecv) {
        extendLinkMap.emplace(ll.prev, ll.next);
        extendLinkReverseMap.emplace(ll.next, ll.prev);
    }
}

void SylinderSystem::addNewCenterLink(const std::vector<Link> &newCenterLink) {
    // synchronize newLink to all mpi ranks
    const int newCountLocal = newCenterLink.size();
    std::vector<int> newCount(commRcp->getSize(), 0);
    MPI_Allgather(&newCountLocal, 1, MPI_INT, newCount.data(), 1, MPI_INT, MPI_COMM_WORLD);
    std::vector<int> displ(commRcp->getSize() + 1, 0);
    std::partial_sum(newCount.cbegin(), newCount.cend(), displ.begin() + 1);
    std::vector<Link> newLinkRecv(displ.back());
    MPI_Allgatherv(newCenterLink.data(), newCountLocal, createMPIStructType<Link>(), newLinkRecv.data(), newCount.data(),
                   displ.data(), createMPIStructType<Link>(), MPI_COMM_WORLD);

    // put newLinks into the map, same op on all mpi ranks
    for (const auto &ll : newLinkRecv) {
        bendLinkMap.emplace(ll.prev, ll.next);
        bendLinkReverseMap.emplace(ll.next, ll.prev);
    }
}

void SylinderSystem::buildSylinderNearDataDirectory() {
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    auto &sylinderNearDataDirectory = *sylinderNearDataDirectoryPtr;
    sylinderNearDataDirectory.gidOnLocal.resize(nLocal);
    sylinderNearDataDirectory.dataOnLocal.resize(nLocal);
#pragma omp parallel for
    for (int i = 0; i < nLocal; i++) {
        sylinderNearDataDirectory.gidOnLocal[i] = sylinderContainer[i].gid;
        sylinderNearDataDirectory.dataOnLocal[i].copyFromFP(sylinderContainer[i]);
    }

    // build index
    sylinderNearDataDirectory.buildIndex();
}

void SylinderSystem::collectPinLinkBilateral() {
    // setup bilateral joint constraints between rod ends
    // need special treatment of periodic boundary conditions

    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    auto &conPool = *(this->conCollectorPtr->constraintPoolPtr);
    if (conPool.size() != omp_get_max_threads()) {
        spdlog::critical("conPool multithread mismatch error");
        std::exit(1);
    }

    // fill the data to find
    auto &gidToFind = sylinderNearDataDirectoryPtr->gidToFind;
    const auto &dataToFind = sylinderNearDataDirectoryPtr->dataToFind;

    std::vector<int> gidDisp(nLocal + 1, 0);
    gidToFind.clear();
    gidToFind.reserve(nLocal);

    // loop over all sylinders
    // if endLinkMap[sy.gid] not empty, find info for all next
    for (int i = 0; i < nLocal; i++) {
        const auto &sy = sylinderContainer[i];
        const auto &range = pinLinkMap.equal_range(sy.gid);
        int count = 0;
        for (auto it = range.first; it != range.second; it++) {
            gidToFind.push_back(it->second); // next
            count++;
        }
        gidDisp[i + 1] = gidDisp[i] + count; // number of links for each local Sylinder
    }

    sylinderNearDataDirectoryPtr->find();

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        auto &conQue = conPool[threadId];
#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            const auto &syI = sylinderContainer[i]; // sylinder
            const int lb = gidDisp[i];
            const int ub = gidDisp[i + 1];

            for (int j = lb; j < ub; j++) {
                const auto &syJ = sylinderNearDataDirectoryPtr->dataToFind[j]; // sylinderNear

                const Evec3 &centerI = ECmap3(syI.pos);
                Evec3 centerJ = ECmap3(syJ.pos);
                // apply PBC on centerJ
                for (int k = 0; k < 3; k++) {
                    if (!runConfig.simBoxPBC[k])
                        continue;
                    double trg = centerI[k];
                    double xk = centerJ[k];
                    findPBCImage(runConfig.simBoxLow[k], runConfig.simBoxHigh[k], xk, trg);
                    centerJ[k] = xk;
                    // error check
                    if (fabs(trg - xk) > 0.5 * (runConfig.simBoxHigh[k] - runConfig.simBoxLow[k])) {
                        spdlog::critical("pbc image error in bilateral links");
                        std::exit(1);
                    }
                }
                // sylinders are not treated as spheres for bilateral constraints
                // constraint is always added between Pp and Qm
                // constraint target length is radiusI + radiusJ + runConfig.endLinkGapcid
                const Evec3 directionI = ECmapq(syI.orientation) * Evec3(0, 0, 1);
                const Evec3 Pp = centerI + directionI * (0.5 * syI.length); // plus end
                const Evec3 directionJ = ECmap3(syJ.direction);
                const Evec3 Qm = centerJ - directionJ * (0.5 * syJ.length);
                const Evec3 Ploc = Pp;
                const Evec3 Qloc = Qm;
                const Evec3 rvec = Ploc - Qloc;
                const Evec3 posI = Ploc - centerI;
                const Evec3 posJ = Qloc - centerJ;

                const std::vector<Evec3> normIVec = {Evec3(1.0, 0.0, 0.0),
                                                    Evec3(0.0, 1.0, 0.0),
                                                    Evec3(0.0, 0.0, 1.0)};

                // We have three pin constraints, one for each dimension. This properly constrains our system.
                for (size_t i = 0; i = 3 ; i++)
                {
                    const double delta0 = rvec[i];
                    const double gammaGuess = 0;
                    const Evec3 unscaledForceComI = normIVec[i];
                    const Evec3 unscaledForceComJ = -unscaledForceComI;
                    const Evec3 unscaledTorqueComI = posI.cross(unscaledForceComI);
                    const Evec3 unscaledTorqueComJ = posJ.cross(unscaledForceComJ);
                    ConstraintBlock conBlock(delta0, gammaGuess,              // current separation, initial guess of gamma
                                            syI.gid, syJ.gid,           //
                                            syI.globalIndex,            //
                                            syJ.globalIndex,            //
                                            unscaledForceComI.data(), unscaledForceComJ.data(), // direction of collision force
                                            unscaledTorqueComI.data(), unscaledTorqueComJ.data(), // location of collision relative to particle center
                                            Ploc.data(), Qloc.data(), // location of collision in lab frame
                                            false, true, 0.0);
                    Emat3 stressIJ;
                    CalcSylinderNearForce::collideStress(directionI, directionJ, centerI, centerJ, syI.length, syJ.length,
                                                        syI.radius, syJ.radius, 1.0, Ploc, Qloc, stressIJ);
                    conBlock.setStress(stressIJ);
                    conQue.push_back(conBlock);
                }
            }
        }
    }
}

void SylinderSystem::collectExtendLinkBilateral() {
    // setup bilateral link constraints
    // need special treatment of periodic boundary conditions

    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    auto &conPool = *(this->conCollectorPtr->constraintPoolPtr);
    if (conPool.size() != omp_get_max_threads()) {
        spdlog::critical("conPool multithread mismatch error");
        std::exit(1);
    }

    // fill the data to find
    auto &gidToFind = sylinderNearDataDirectoryPtr->gidToFind;
    const auto &dataToFind = sylinderNearDataDirectoryPtr->dataToFind;

    std::vector<int> gidDisp(nLocal + 1, 0);
    gidToFind.clear();
    gidToFind.reserve(nLocal);

    // loop over all sylinders
    // if extendlinkMap[sy.gid] not empty, find info for all next
    for (int i = 0; i < nLocal; i++) {
        const auto &sy = sylinderContainer[i];
        const auto &range = extendLinkMap.equal_range(sy.gid);
        int count = 0;
        for (auto it = range.first; it != range.second; it++) {
            gidToFind.push_back(it->second); // next
            count++;
        }
        gidDisp[i + 1] = gidDisp[i] + count; // number of links for each local Sylinder
    }

    sylinderNearDataDirectoryPtr->find();

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        auto &conQue = conPool[threadId];
#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            const auto &syI = sylinderContainer[i]; // sylinder
            const int lb = gidDisp[i];
            const int ub = gidDisp[i + 1];

            for (int j = lb; j < ub; j++) {
                const auto &syJ = sylinderNearDataDirectoryPtr->dataToFind[j]; // sylinderNear

                const Evec3 &centerI = ECmap3(syI.pos);
                Evec3 centerJ = ECmap3(syJ.pos);
                // apply PBC on centerJ
                for (int k = 0; k < 3; k++) {
                    if (!runConfig.simBoxPBC[k])
                        continue;
                    double trg = centerI[k];
                    double xk = centerJ[k];
                    findPBCImage(runConfig.simBoxLow[k], runConfig.simBoxHigh[k], xk, trg);
                    centerJ[k] = xk;
                    // error check
                    if (fabs(trg - xk) > 0.5 * (runConfig.simBoxHigh[k] - runConfig.simBoxLow[k])) {
                        spdlog::critical("pbc image error in bilateral links");
                        std::exit(1);
                    }
                }
                // sylinders are not treated as spheres for bilateral constraints
                // constraint is always added between Pp and Qm
                // constraint target length is radiusI + radiusJ + runConfig.linkGap
                const Evec3 directionI = ECmapq(syI.orientation) * Evec3(0, 0, 1);
                const Evec3 Pp = centerI + directionI * (0.5 * syI.length); // plus end
                const Evec3 directionJ = ECmap3(syJ.direction);
                const Evec3 Qm = centerJ - directionJ * (0.5 * syJ.length);
                const Evec3 Ploc = Pp;
                const Evec3 Qloc = Qm;
                const Evec3 rvec = Qloc - Ploc;
                
                const double delta0 = rvec.norm() - syI.radius - syJ.radius - runConfig.extendLinkGap;
                const double gamma = delta0 < 0 ? -delta0 : 0;
                const Evec3 normI = (Ploc - Qloc).normalized();
                const Evec3 normJ = -normI;
                const Evec3 posI = Ploc - centerI;
                const Evec3 posJ = Qloc - centerJ;
                const Evec3 unscaledForceComI = normI;
                const Evec3 unscaledForceComJ = -unscaledForceComI;
                const Evec3 unscaledTorqueComI = posI.cross(unscaledForceComI);
                const Evec3 unscaledTorqueComJ = posJ.cross(unscaledForceComJ);

                ConstraintBlock conBlock(delta0, gamma,              // current separation, initial guess of gamma
                                         syI.gid, syJ.gid,           //
                                         syI.globalIndex,            //
                                         syJ.globalIndex,            //
                                         unscaledForceComI.data(), unscaledForceComJ.data(), // direction of collision force
                                         unscaledTorqueComI.data(), unscaledTorqueComJ.data(), // location of collision relative to particle center
                                         Ploc.data(), Qloc.data(), // location of collision in lab frame
                                         false, true, runConfig.extendLinkKappa);
                Emat3 stressIJ;
                CalcSylinderNearForce::collideStress(directionI, directionJ, centerI, centerJ, syI.length, syJ.length,
                                                     syI.radius, syJ.radius, 1.0, Ploc, Qloc, stressIJ);
                conBlock.setStress(stressIJ);
                conQue.push_back(conBlock);
            }
        }
    }
}

void SylinderSystem::collectBendLinkBilateral() {
    // setup bilateral angular spring constraints between rod centers
    // need special treatment of periodic boundary conditions
    const double t = runConfig.dt * stepCount;

    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    auto &conPool = *(this->conCollectorPtr->constraintPoolPtr);
    if (conPool.size() != omp_get_max_threads()) {
        spdlog::critical("conPool multithread mismatch error");
        std::exit(1);
    }

    // fill the data to find
    auto &gidToFind = sylinderNearDataDirectoryPtr->gidToFind;
    const auto &dataToFind = sylinderNearDataDirectoryPtr->dataToFind;

    std::vector<int> gidDisp(nLocal + 1, 0);
    gidToFind.clear();
    gidToFind.reserve(nLocal);

    // loop over all sylinders
    // if bendLinkMap[sy.gid] not empty, find info for all next
    for (int i = 0; i < nLocal; i++) {
        const auto &sy = sylinderContainer[i];
        const auto &range = bendLinkMap.equal_range(sy.gid);
        int count = 0;
        for (auto it = range.first; it != range.second; it++) {
            gidToFind.push_back(it->second); // next
            count++;
        }
        gidDisp[i + 1] = gidDisp[i] + count; // number of links for each local Sylinder
    }

    sylinderNearDataDirectoryPtr->find();

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        auto &conQue = conPool[threadId];
#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            const auto &syI = sylinderContainer[i]; // sylinder
            const int lb = gidDisp[i];
            const int ub = gidDisp[i + 1];

            for (int j = lb; j < ub; j++) {
                const auto &syJ = sylinderNearDataDirectoryPtr->dataToFind[j]; // sylinderNear

                const Evec3 &centerI = ECmap3(syI.pos);
                Evec3 centerJ = ECmap3(syJ.pos);
                // apply PBC on centerJ
                for (int k = 0; k < 3; k++) {
                    if (!runConfig.simBoxPBC[k])
                        continue;
                    double trg = centerI[k];
                    double xk = centerJ[k];
                    findPBCImage(runConfig.simBoxLow[k], runConfig.simBoxHigh[k], xk, trg);
                    centerJ[k] = xk;
                    // error check
                    if (fabs(trg - xk) > 0.5 * (runConfig.simBoxHigh[k] - runConfig.simBoxLow[k])) {
                        spdlog::critical("pbc image error in bilateral links");
                        std::exit(1);
                    }
                }
                // sylinders are not treated as spheres for bilateral constraints
                // constraint is always added between Pp and Qm
                // constraint target length is radiusI + radiusJ + runConfig.linkGap
                const Evec3 directionI = ECmapq(syI.orientation) * Evec3(0, 0, 1);
                const Evec3 Pp = centerI + directionI * (0.5 * syI.length); // plus end
                const Evec3 directionJ = ECmap3(syJ.direction);
                const Evec3 Qm = centerJ - directionJ * (0.5 * syJ.length);
                const Evec3 Ploc = Pp;
                const Evec3 Qloc = Qm;

                // We have three bending constraints.
                // They correspond to the satisfaction of kappa_i - kappa_0i - lambda_i / b_i along the triad of directors at the node connecting the two rods.
                // In our case, the triad is uniquely described by the unit quaternion associated with each rod segment. To find the average quaternion at the center of the two nodes, we use slerp. 
                const auto equatI = ECmapq(syI.orientation);
                const auto equatJ = ECmapq(syJ.orientation);
                const auto equatIJ = equatI.slerp(0.5, equatJ).normalized();

                // The triad of directors is given by the following three unit vectors:
                const std::vector<Evec3> dirIJVec = {equatIJ * Evec3(1, 0, 0), 
                                                     equatIJ * Evec3(0, 1, 0), 
                                                     equatIJ * Evec3(0, 0, 1)};

                // The curvature of the two rods is given by kappa = vec[equatI.conjugate() * equatJ - equatI * equatJ.conjugate()]
                // where vec is the vector part of a quaternion.
                const Evec3 curvature = (equatI.conjugate() * equatJ).vec() - (equatI * equatJ.conjugate()).vec();

                for (size_t i = 0; i < 3; i++)
                {
                    // The first constraint is given by kappa_1 - kappa_01 - lambda_1 / b_1 = 0
                    // This constraint is applied along the first director.
                    const double delta0 = curvature[i] - runConfig.preferredCurvature[i];
                    const Evec3 unscaledForceComI(0.0, 0.0, 0.0);
                    const Evec3 unscaledForceComJ(0.0, 0.0, 0.0);
                    const Evec3 unscaledTorqueComI = -dirIJVec[i];
                    const Evec3 unscaledTorqueComJ = -unscaledTorqueComI;
                    const double gammaGuess = 0;
                    ConstraintBlock conBlock(delta0, gammaGuess,        // current separation, initial guess of gamma
                                            syI.gid, syJ.gid,           //
                                            syI.globalIndex,            //
                                            syJ.globalIndex,            //
                                            unscaledForceComI.data(), unscaledForceComJ.data(), // direction of collision force
                                            unscaledTorqueComI.data(), unscaledTorqueComJ.data(), // location of collision relative to particle center
                                            Ploc.data(), Qloc.data(), // location of collision in lab frame
                                            false, true, runConfig.bendLinkKappa[i]);
                    Emat3 stressIJ = Emat3::Zero();
                    conBlock.setStress(stressIJ);
                    conQue.push_back(conBlock);
                }
            }
        }
    }
}

void SylinderSystem::collectTriBendLinkBilateral() {
    // setup bilateral angular spring constraints between three spheres
    // need special treatment of periodic boundary conditions
    const int nLocal = sylinderContainer.getNumberOfParticleLocal();
    auto &conPool = *(this->conCollectorPtr->constraintPoolPtr);
    if (conPool.size() != omp_get_max_threads()) {
        spdlog::critical("conPool multithread mismatch error");
        std::exit(1);
    }

    // fill the data to find
    auto &gidToFind = sylinderNearDataDirectoryPtr->gidToFind;
    const auto &dataToFind = sylinderNearDataDirectoryPtr->dataToFind;

    std::vector<int> gidDisp(nLocal + 1, 0);
    gidToFind.clear();
    gidToFind.reserve(nLocal);

    // loop over all sylinders
    // if tribendLinkMap[sy.gid] not empty, find info for all next
    for (int i = 0; i < nLocal; i++) {
        const auto &sy = sylinderContainer[i];
        const auto& [rangeStart, rangeStop] = tribendLinkMap.equal_range(sy.gid);
        int count = 0;
        for (auto it = rangeStart; it != rangeStop; it++) {
            const auto& [key, pair] = *it;
            const auto& [gidJ, gidK] = pair;

            gidToFind.push_back(gidJ); // left sphere
            gidToFind.push_back(gidK); // right sphere
            count ++;
        }
        gidDisp[i + 1] = gidDisp[i] + count; // number of links for each local Sylinder
    }

    sylinderNearDataDirectoryPtr->find();

#pragma omp parallel
    {
        const int threadId = omp_get_thread_num();
        auto &conQue = conPool[threadId];
#pragma omp for
        for (int i = 0; i < nLocal; i++) {
            const auto &syI = sylinderContainer[i]; // sylinder
            const int lb = gidDisp[i];
            const int ub = gidDisp[i + 1];
            const Evec3 &centerI = ECmap3(syI.pos);

            for (int j = lb; j < ub; j++) {
                const auto &syJ = sylinderNearDataDirectoryPtr->dataToFind[2 * j + 0]; // left sphere
                const auto &syK = sylinderNearDataDirectoryPtr->dataToFind[2 * j + 1]; // right sphere

                Evec3 centerJ = ECmap3(syJ.pos);
                Evec3 centerK = ECmap3(syK.pos);

                // apply PBC on centerJ
                for (int k = 0; k < 3; k++) {
                    if (!runConfig.simBoxPBC[k])
                        continue;
                    double trg = centerI[k];
                    double srcJ = centerJ[k];
                    double srcK = centerK[k];
                    findPBCImage(runConfig.simBoxLow[k], runConfig.simBoxHigh[k], srcJ, trg);
                    findPBCImage(runConfig.simBoxLow[k], runConfig.simBoxHigh[k], srcK, trg);
                    centerJ[k] = srcJ;
                    centerK[k] = srcK;
                    // error check
                    if (fabs(trg - srcJ) > 0.5 * (runConfig.simBoxHigh[k] - runConfig.simBoxLow[k])) {
                        spdlog::critical("pbc image error in bilateral links");
                        std::exit(1);
                    }
                    if (fabs(trg - srcK) > 0.5 * (runConfig.simBoxHigh[k] - runConfig.simBoxLow[k])) {
                        spdlog::critical("pbc image error in bilateral links");
                        std::exit(1);
                    }
                }
                // sylinders are not treated as spheres for bilateral constraints
                // constraint is always added between the center of each rod

                // We have three bending constraints.
                // They correspond to the satisfaction of kappa_i - kappa_0i - lambda_i / b_i along the triad of directors at the node connecting the two rods.
                // In our case, the triad is uniquely described by the unit quaternion associated with each rod segment. To find the average quaternion at the center of the two nodes, we use slerp. 
                Evec3 orientVecJI = centerI - centerJ;
                Evec3 orientVecIK = centerK - centerI;
                double distJI = orientVecJI.norm();
                double distIK = orientVecIK.norm();
                orientVecJI /= distJI;
                orientVecIK /= distIK;

                const auto equatI = Equatn::FromTwoVectors(Evec3(0, 0, 1), orientVecJI);
                const auto equatJ = Equatn::FromTwoVectors(Evec3(0, 0, 1), orientVecIK);
                const auto equatIJ = equatI.slerp(0.5, equatJ).normalized();

                // The triad of directors is given by the following three unit vectors:
                const std::vector<Evec3> dirIJVec = {equatIJ * Evec3(1, 0, 0), 
                                                     equatIJ * Evec3(0, 1, 0), 
                                                     equatIJ * Evec3(0, 0, 1)};

                // The momenet of inertia of the lines connecting the rods
                Emat3 momIntJI = distJI * distJI * (orientVecJI * orientVecJI.transpose() - Emat3::Identity());
                Emat3 momIntIK = distIK * distIK * (orientVecIK * orientVecIK.transpose() - Emat3::Identity());
                Emat3 momIntJIinv = momIntJI.inverse();
                Emat3 momIntIKinv = momIntIK.inverse();

                // The curvature of the three spheres is given by kappa = vec[equatI.conjugate() * equatJ - equatI * equatJ.conjugate()]
                // where vec is the vector part of a quaternion.
                const Evec3 curvature = (equatI.conjugate() * equatJ).vec() - (equatI * equatJ.conjugate()).vec();

                for (size_t i = 0; i < 3; i++)
                {
                    const double delta0 = curvature[i] - runConfig.preferredCurvature[i];
                    const Evec3 unscaledTorqueComBetweenJandI = -dirIJVec[i] ;
                    const Evec3 unscaledTorqueComBetweenIandK = -unscaledTorqueComBetweenJandI;
                    const Evec3 unscaledForceComJ = -momIntJIinv * (distJI*orientVecJI).cross(unscaledTorqueComBetweenJandI);
                    const Evec3 unscaledForceComI = momIntJIinv * (distJI*orientVecJI).cross(unscaledTorqueComBetweenJandI) - momIntIKinv * (distIK*orientVecIK).cross(unscaledTorqueComBetweenIandK);
                    const Evec3 unscaledForceComK = momIntIKinv * (distIK*orientVecIK).cross(unscaledTorqueComBetweenIandK);
                    const Evec3 unscaledTorqueComI(0.0, 0.0, 0.0);
                    const Evec3 unscaledTorqueComJ(0.0, 0.0, 0.0);
                    const Evec3 unscaledTorqueComK(0.0, 0.0, 0.0);

                    const double gammaGuess = 0;
                    ConstraintBlock conBlock(delta0, gammaGuess,        // current separation, initial guess of gamma
                                            syI.gid, syJ.gid, syK.gid,           //
                                            syI.globalIndex, syJ.globalIndex, syK.globalIndex,         //
                                            unscaledForceComI.data(), unscaledForceComJ.data(), unscaledForceComK.data(), // direction of collision force
                                            unscaledTorqueComI.data(), unscaledTorqueComJ.data(), unscaledTorqueComK.data(), // location of collision relative to particle center
                                            centerI.data(), centerJ.data(), centerK.data(), // location of collision in lab frame
                                            false, true, runConfig.tribendLinkKappa[i]);
                    Emat3 stressIJ = Emat3::Zero();
                    conBlock.setStress(stressIJ);
                    conQue.push_back(conBlock);
                }
            }
        }
    }
}


void SylinderSystem::printTimingSummary(const bool zeroOut) {
    if (runConfig.timerLevel <= spdlog::level::info)
        Teuchos::TimeMonitor::summarize();
    if (zeroOut)
        Teuchos::TimeMonitor::zeroOutTimers();
}

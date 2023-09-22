/**
 * @file SylinderSystem_test_api.cpp
 * @author Wen Yan (wenyan4work@gmail.com)
 * @brief This file tests a few api behaviors of SylinderSystem class
 * @version 0.1
 * @date 2021-02-26
 *
 * @copyright Copyright (c) 2021
 *
 */
#include "SylinderSystem.hpp"

#include "MPI/CommMPI.hpp"
#include "Util/Logger.hpp"

void testAddLinks(SylinderSystem &sylinderSystem) {
    /*********************************
     * Test:
     * (1) run 10 steps as given in RunConfig.yaml
     * (2) add one link of 20 segments from each rank
     * (3) run 10 more steps
     *
     * Check:
     * (1) gid is correct
     * (2) linkage is correct
     */

    auto &rngPoolPtr = sylinderSystem.getRngPoolPtr();
    auto &runConfig = sylinderSystem.runConfig;

    // run 10 steps for relaxation
    for (int i = 0; i < 10; i++) {
        sylinderSystem.prepareStep();
        sylinderSystem.runStep();
    }

    // add linked sylinders
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // using rank as group id
    const int nLocalNew = 20;
    std::vector<Sylinder> newSylinder(nLocalNew);
    for (int i = 0; i < nLocalNew; i++) {
        newSylinder[i] = Sylinder(i, runConfig.sylinderDiameter / 2, runConfig.sylinderDiameter / 2,
                                  runConfig.sylinderLength / 2, runConfig.sylinderLength / 2);
    }

    // specify linked sylinders
    newSylinder[0].pos[0] = rngPoolPtr->getU01(0) * runConfig.simBoxHigh[0];
    newSylinder[0].pos[1] = rngPoolPtr->getU01(0) * runConfig.simBoxHigh[1];
    newSylinder[0].pos[2] = rngPoolPtr->getU01(0) * runConfig.simBoxHigh[2];
    Equatn qtemp;
    EquatnHelper::setUnitRandomEquatn(qtemp, rngPoolPtr->getU01(0), rngPoolPtr->getU01(0), rngPoolPtr->getU01(0));
    Emapq(newSylinder[0].orientation).coeffs() = qtemp.coeffs();
    for (int i = 1; i < nLocalNew; i++) {
        const auto &prsy = newSylinder[i - 1];
        Evec3 prvec = ECmapq(prsy.orientation) * Evec3(0, 0, 1);
        Evec3 prend = prvec * 0.5 * prsy.length + ECmap3(prsy.pos);
        auto &sy = newSylinder[i];
        Equatn qtemp;
        EquatnHelper::setUnitRandomEquatn(qtemp, rngPoolPtr->getU01(0), rngPoolPtr->getU01(0), rngPoolPtr->getU01(0));
        Emapq(sy.orientation).coeffs() = qtemp.coeffs();
        Evec3 syvec = Emapq(sy.orientation) * Evec3(0, 0, 1);
        Evec3 sypos = prend + (prsy.radius + sy.radius) * prvec + sy.length * 0.5 * syvec;
        Emap3(sy.pos) = sypos;
    }

    const auto &newGid = sylinderSystem.addNewSylinder(newSylinder);
    assert(nLocalNew == newGid.size());
    // sequentially linked
    std::vector<Link> linkage(nLocalNew - 1);
    for (int i = 0; i < nLocalNew - 1; i++) {
        linkage[i].prev = newGid[i];
        linkage[i].next = newGid[i + 1];
    }
    sylinderSystem.addNewEndLink(linkage);

    auto &linkMap = sylinderSystem.getEndLinkMap();
    for (auto &pn : linkMap) {
        printf("L %d, %d, %d\n", rank, pn.first, pn.second);
    }

    // run 10 steps
    const int nSteps = 10;
    for (int i = 0; i < nSteps; i++) {
        sylinderSystem.prepareStep();
        auto &sylinderContainer = sylinderSystem.getContainer();
        int nLocal = sylinderContainer.getNumberOfParticleLocal();
        std::vector<double> forceNonBrown(nLocal * 6, 0.0);
        for (int i = 0; i < nLocal; i++) {
            const auto &gid = sylinderContainer[i].gid;
            if (linkMap.count(gid) > 0) {
                forceNonBrown[6 * i + 1] = 10; // y
                forceNonBrown[6 * i + 2] = 10; // z
            }
        }
        sylinderSystem.setForceNonBrown(forceNonBrown);
        sylinderSystem.runStep();
    }

    // check added links by external python script
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    Logger::setup_mpi_spdlog();

    {
        auto runConfig = SylinderConfig("RunConfig.yaml");

        SylinderSystem sylinderSystem(runConfig, "posInitial.dat", argc, argv);
        sylinderSystem.setTimer(true);
        testAddLinks(sylinderSystem);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
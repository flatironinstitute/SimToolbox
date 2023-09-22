#include "ConstraintCollector.hpp"
#include "spdlog/spdlog.h"

#include <cstdlib>

ConstraintCollector::ConstraintCollector() {
    const int totalThreads = omp_get_max_threads();
    constraintPoolPtr = std::make_shared<ConstraintBlockPool>();
    constraintPoolPtr->resize(totalThreads);
    for (auto &queue : *constraintPoolPtr) {
        queue.clear();
    }

    spdlog::debug("ConstraintCollector constructed for {} threads", constraintPoolPtr->size());
}

bool ConstraintCollector::valid() const { return constraintPoolPtr->empty(); }

void ConstraintCollector::clear() {
    assert(constraintPoolPtr);
    for (int i = 0; i < constraintPoolPtr->size(); i++) {
        (*constraintPoolPtr)[i].clear();
    }

    // keep the total number of queues
    const int totalThreads = omp_get_max_threads();
    constraintPoolPtr->resize(totalThreads);
}

int ConstraintCollector::getLocalNumberOfConstraints() {
    int sum = 0;
    for (int i = 0; i < constraintPoolPtr->size(); i++) {
        sum += (*constraintPoolPtr)[i].size();
    }
    return sum;
}

void ConstraintCollector::sumLocalConstraintStress(Emat3 &uniStress, Emat3 &biStress, bool withOneSide) const {
    const auto &conPool = *constraintPoolPtr;
    const int poolSize = conPool.size();

    Emat3 biStressTotal = Emat3::Zero();
    Emat3 uniStressTotal = Emat3::Zero();

#pragma omp parallel for
    for (int que = 0; que < poolSize; que++) {
        // reduction of stress, one que for each thread

        Emat3 uniStressSumQue = Emat3::Zero();
        Emat3 biStressSumQue = Emat3::Zero();
        const auto &conQue = conPool[que];
        for (const auto &cBlock : conQue) {
            if (cBlock.oneSide && !withOneSide) {
                // skip counting oneside collision blocks
                continue;
            } else {
                Emat3 stressBlock;
                cBlock.getStress(stressBlock);
                if (cBlock.bilateral) {
                    biStressSumQue = biStressSumQue + stressBlock;
                } else {
                    uniStressSumQue = uniStressSumQue + stressBlock;
                }
            }
        }
#pragma omp critical
        {
            biStressTotal = biStressTotal + biStressSumQue;
            uniStressTotal = uniStressTotal + uniStressSumQue;
        }
    }
    uniStress = uniStressTotal;
    biStress = biStressTotal;
}

void ConstraintCollector::writePVTP(const std::string &folder, const std::string &prefix, const std::string &postfix,
                                    const int nProcs) const {
    std::vector<std::string> pieceNames;

    std::vector<IOHelper::FieldVTU> pointDataFields;
    pointDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "gid");
    pointDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "globalIndex");
    pointDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "unscaledForceComIJ");
    pointDataFields.emplace_back(3, IOHelper::IOTYPE::Float32, "unscaledTorqueComIJ");

    std::vector<IOHelper::FieldVTU> cellDataFields;
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "oneSide");
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Int32, "bilateral");
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "delta0");
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "gamma");
    cellDataFields.emplace_back(1, IOHelper::IOTYPE::Float32, "kappa");
    cellDataFields.emplace_back(9, IOHelper::IOTYPE::Float32, "Stress");

    for (int i = 0; i < nProcs; i++) {
        pieceNames.emplace_back(prefix + std::string("ConBlock_") + std::string("r") + std::to_string(i) + "_" +
                                postfix + ".vtp");
    }

    IOHelper::writePVTPFile(folder + "/" + prefix + std::string("ConBlock_") + postfix + ".pvtp", pointDataFields,
                            cellDataFields, pieceNames);
}

void ConstraintCollector::writeVTP(const std::string &folder, const std::string &prefix, const std::string &postfix,
                                   int rank) const {
    const auto &conPool = *constraintPoolPtr;
    const int conQueNum = conPool.size();
    std::vector<int> conQueSize;
    std::vector<int> conQueIndex;
    buildConIndex(conQueSize, conQueIndex);
    const int cBlockNum = conQueIndex.back();

    // for each block:

    // write VTP for basic data, use float to save some space
    // point and point data
    std::vector<double> pos(6 * cBlockNum); // position always in Float64
    std::vector<int32_t> gid(2 * cBlockNum);
    std::vector<int32_t> globalIndex(2 * cBlockNum);
    std::vector<float> unscaledForceComIJ(6 * cBlockNum);
    std::vector<float> unscaledTorqueComIJ(6 * cBlockNum);

    // point connectivity of line
    std::vector<int32_t> connectivity(2 * cBlockNum);
    std::vector<int32_t> offset(cBlockNum);

    // cell data for ColBlock
    std::vector<int32_t> oneSide(cBlockNum);
    std::vector<int32_t> bilateral(cBlockNum);
    std::vector<float> delta0(cBlockNum);
    std::vector<float> gamma(cBlockNum);
    std::vector<float> kappa(cBlockNum);
    std::vector<float> Stress(9 * cBlockNum);

#pragma omp parallel for
    // connectivity
    for (int i = 0; i < cBlockNum; i++) {
        connectivity[2 * i] = 2 * i;         // index of point 0 in line
        connectivity[2 * i + 1] = 2 * i + 1; // index of point 1 in line
        offset[i] = 2 * i + 2;               // offset is the end of each line. in fortran indexing
    }

#pragma omp parallel for
    // data
    for (int q = 0; q < conQueNum; q++) {
        const int queSize = conQueSize[q];
        const int queIndex = conQueIndex[q];
        for (int c = 0; c < queSize; c++) {
            const int cIndex = queIndex + c;
            const auto &block = conPool[q][c];
            // point position
            pos[6 * cIndex + 0] = block.labI[0];
            pos[6 * cIndex + 1] = block.labI[1];
            pos[6 * cIndex + 2] = block.labI[2];
            pos[6 * cIndex + 3] = block.labJ[0];
            pos[6 * cIndex + 4] = block.labJ[1];
            pos[6 * cIndex + 5] = block.labJ[2];
            // point data
            gid[2 * cIndex + 0] = block.gidI;
            gid[2 * cIndex + 1] = block.gidJ;
            globalIndex[2 * cIndex + 0] = block.globalIndexI;
            globalIndex[2 * cIndex + 1] = block.globalIndexJ;
            unscaledForceComIJ[6 * cIndex + 0] = block.unscaledForceComI[0];
            unscaledForceComIJ[6 * cIndex + 1] = block.unscaledForceComI[1];
            unscaledForceComIJ[6 * cIndex + 2] = block.unscaledForceComI[2];
            unscaledForceComIJ[6 * cIndex + 3] = block.unscaledForceComJ[0];
            unscaledForceComIJ[6 * cIndex + 4] = block.unscaledForceComJ[1];
            unscaledForceComIJ[6 * cIndex + 5] = block.unscaledForceComJ[2];
            unscaledTorqueComIJ[6 * cIndex + 0] = block.unscaledTorqueComI[0];
            unscaledTorqueComIJ[6 * cIndex + 1] = block.unscaledTorqueComI[1];
            unscaledTorqueComIJ[6 * cIndex + 2] = block.unscaledTorqueComI[2];
            unscaledTorqueComIJ[6 * cIndex + 3] = block.unscaledTorqueComJ[0];
            unscaledTorqueComIJ[6 * cIndex + 4] = block.unscaledTorqueComJ[1];
            unscaledTorqueComIJ[6 * cIndex + 5] = block.unscaledTorqueComJ[2];
            // cell data
            oneSide[cIndex] = block.oneSide ? 1 : 0;
            bilateral[cIndex] = block.bilateral ? 1 : 0;
            delta0[cIndex] = block.delta0;
            gamma[cIndex] = block.gamma;
            kappa[cIndex] = block.kappa;
            for (int kk = 0; kk < 9; kk++) {
                Stress[9 * cIndex + kk] = block.stress[kk];
            }
        }
    }

    std::ofstream file(folder + '/' + prefix + std::string("ConBlock_") + "r" + std::to_string(rank) +
                           std::string("_") + postfix + std::string(".vtp"),
                       std::ios::out);

    IOHelper::writeHeadVTP(file);

    // Piece starts
    file << "<Piece NumberOfPoints=\"" << cBlockNum * 2 << "\" NumberOfLines=\"" << cBlockNum << "\">\n";
    // Points
    file << "<Points>\n";
    IOHelper::writeDataArrayBase64(pos, "position", 3, file);
    file << "</Points>\n";
    // cell definition
    file << "<Lines>\n";
    IOHelper::writeDataArrayBase64(connectivity, "connectivity", 1, file);
    IOHelper::writeDataArrayBase64(offset, "offsets", 1, file);
    file << "</Lines>\n";
    // point data
    file << "<PointData Scalars=\"scalars\">\n";
    IOHelper::writeDataArrayBase64(gid, "gid", 1, file);
    IOHelper::writeDataArrayBase64(globalIndex, "globalIndex", 1, file);
    IOHelper::writeDataArrayBase64(unscaledForceComIJ, "unscaledForceComIJ", 3, file);
    IOHelper::writeDataArrayBase64(unscaledTorqueComIJ, "unscaledTorqueComIJ", 3, file);
    file << "</PointData>\n";
    // cell data
    file << "<CellData Scalars=\"scalars\">\n";
    IOHelper::writeDataArrayBase64(oneSide, "oneSide", 1, file);
    IOHelper::writeDataArrayBase64(bilateral, "bilateral", 1, file);
    IOHelper::writeDataArrayBase64(delta0, "delta0", 1, file);
    IOHelper::writeDataArrayBase64(gamma, "gamma", 1, file);
    IOHelper::writeDataArrayBase64(kappa, "kappa", 1, file);
    IOHelper::writeDataArrayBase64(Stress, "Stress", 9, file);
    file << "</CellData>\n";
    // Piece ends
    file << "</Piece>\n";

    IOHelper::writeTailVTP(file);
    file.close();
}

void ConstraintCollector::dumpBlocks() const {
    std::cout << "number of collision queues: " << constraintPoolPtr->size() << std::endl;
    // dump constraint blocks
    for (const auto &blockQue : (*constraintPoolPtr)) {
        std::cout << blockQue.size() << " constraints in this queue" << std::endl;
        for (const auto &block : blockQue) {
            std::cout << block.globalIndexI << " " << block.globalIndexJ << "  delta0:" << block.delta0 << std::endl;
        }
    }
}

int ConstraintCollector::buildConstraintMatrixVector(const Teuchos::RCP<const TMAP> &mobMapRcp, //
                                                     Teuchos::RCP<TCMAT> &DMatTransRcp,         //
                                                     Teuchos::RCP<TV> &delta0Rcp,               //
                                                     Teuchos::RCP<TV> &invKappaRcp,             //
                                                     Teuchos::RCP<TV> &biFlagRcp,               //
                                                     Teuchos::RCP<TV> &gammaGuessRcp) const {
    Teuchos::RCP<const TCOMM> commRcp = mobMapRcp->getComm();

    const auto &conPool = *constraintPoolPtr; // the constraint pool
    const int conQueNum = conPool.size();

    // prepare 1, build the index for block queue
    std::vector<int> conQueSize;
    std::vector<int> conQueIndex;
    buildConIndex(conQueSize, conQueIndex);

    // prepare 2, allocate the map and vectors
    const int localGammaSize = conQueIndex.back();
    Teuchos::RCP<const TMAP> gammaMapRcp = getTMAPFromLocalSize(localGammaSize, commRcp);

    // step 1, count the number of entries to each row
    // each constraint block, correspoding to a gamma, occupies a row
    // 18 entries for two side constraint blocks
    // 9 entries for one side constraint blocks
    Kokkos::View<size_t *> rowPointers("rowPointers", localGammaSize + 1); // last entry is the total nnz in this matrix
    rowPointers[0] = 0;

    std::vector<int> colIndexPool(conQueNum + 1, 0); // The beginning index in crs columnIndices for each constraint queue

    int rowPointerIndex = 0;
    int colIndexCount = 0;
    for (int i = 0; i < conQueNum; i++) {
        const auto &conQue = conPool[i];
        const int jsize = conQue.size();
        for (int j = 0; j < jsize; j++) {
            rowPointerIndex++;
            int cBlockNNZ;
            if (conQue[j].oneSide) {
                cBlockNNZ = 6;
            } else if (conQue[j].gidK == GEO_INVALID_INDEX) {
                cBlockNNZ = 12;
            } else {
                cBlockNNZ = 16;
            }
            rowPointers[rowPointerIndex] = rowPointers[rowPointerIndex - 1] + cBlockNNZ;
            colIndexCount += cBlockNNZ;
        }
        colIndexPool[i + 1] = colIndexCount;
    }

    if (rowPointerIndex != localGammaSize) {
        spdlog::critical("rowPointerIndexError in collision solver");
        std::exit(1);
    }

    // step 2, fill the values to each row
    Kokkos::View<int *> columnIndices("columnIndices", rowPointers[localGammaSize]);
    Kokkos::View<double *> values("values", rowPointers[localGammaSize]);
    // multi-thread filling. nThreads = poolSize, each thread process a queue
#pragma omp parallel for num_threads(conQueNum)
    for (int threadId = 0; threadId < conQueNum; threadId++) {
        // each thread process a queue
        const auto &conQue = conPool[threadId];
        const int conNum = conQue.size();
        int kk = colIndexPool[threadId];

        for (int j = 0; j < conNum; j++) {
            const auto &con = conQue[j];

            // 6 nnz for I
            columnIndices[kk + 0] = 6 * con.globalIndexI + 0;
            columnIndices[kk + 1] = 6 * con.globalIndexI + 1;
            columnIndices[kk + 2] = 6 * con.globalIndexI + 2;
            columnIndices[kk + 3] = 6 * con.globalIndexI + 3;
            columnIndices[kk + 4] = 6 * con.globalIndexI + 4;
            columnIndices[kk + 5] = 6 * con.globalIndexI + 5;
            values[kk + 0] = con.unscaledForceComI[0];
            values[kk + 1] = con.unscaledForceComI[1];
            values[kk + 2] = con.unscaledForceComI[2];
            values[kk + 3] = con.unscaledTorqueComI[0];
            values[kk + 4] = con.unscaledTorqueComI[1];
            values[kk + 5] = con.unscaledTorqueComI[2];
            
            kk += 6;
            if (!con.oneSide) {
                // 6 nnz for J
                columnIndices[kk + 0] = 6 * con.globalIndexJ + 0;
                columnIndices[kk + 1] = 6 * con.globalIndexJ + 1;
                columnIndices[kk + 2] = 6 * con.globalIndexJ + 2;
                columnIndices[kk + 3] = 6 * con.globalIndexJ + 3;
                columnIndices[kk + 4] = 6 * con.globalIndexJ + 4;
                columnIndices[kk + 5] = 6 * con.globalIndexJ + 5;
                values[kk + 0] = con.unscaledForceComJ[0];
                values[kk + 1] = con.unscaledForceComJ[1];
                values[kk + 2] = con.unscaledForceComJ[2];
                values[kk + 3] = con.unscaledTorqueComJ[0];
                values[kk + 4] = con.unscaledTorqueComJ[1];
                values[kk + 5] = con.unscaledTorqueComJ[2];
                kk += 6;

                if (con.gidK != GEO_INVALID_INDEX) {
                    // 6 nnz for K
                    columnIndices[kk + 0] = 6 * con.globalIndexK + 0;
                    columnIndices[kk + 1] = 6 * con.globalIndexK + 1;
                    columnIndices[kk + 2] = 6 * con.globalIndexK + 2;
                    columnIndices[kk + 3] = 6 * con.globalIndexK + 3;
                    columnIndices[kk + 4] = 6 * con.globalIndexK + 4;
                    columnIndices[kk + 5] = 6 * con.globalIndexK + 5;
                    values[kk + 0] = con.unscaledForceComK[0];
                    values[kk + 1] = con.unscaledForceComK[1];
                    values[kk + 2] = con.unscaledForceComK[2];
                    values[kk + 3] = con.unscaledTorqueComK[0];
                    values[kk + 4] = con.unscaledTorqueComK[1];
                    values[kk + 5] = con.unscaledTorqueComK[2];
                    kk += 6;
                }
            }
        }
    }

    // step 3 prepare the partitioned column map
    // Each process own some columns. In the map, processes share entries.
    // 3.1 column map has to cover the contiguous range of the mobility map locally owned
    const int mobMin = mobMapRcp->getMinGlobalIndex();
    const int mobMax = mobMapRcp->getMaxGlobalIndex();
    std::vector<int> colMapIndex(mobMax - mobMin + 1);
#pragma omp parallel for
    for (int i = mobMin; i <= mobMax; i++) {
        colMapIndex[i - mobMin] = i;
    }
    // this is the list of the columns that have nnz entries
    // if the column index is out of [mobMinLID, mobMaxLID], add it to the map
    const int colIndexNum = columnIndices.dimension_0();
    if (colIndexNum != colIndexCount) {
        spdlog::critical("colIndexNum error");
        std::exit(1);
    }
    for (int i = 0; i < colIndexNum; i++) {
        if (columnIndices[i] < mobMin || columnIndices[i] > mobMax)
            colMapIndex.push_back(columnIndices[i]);
    }

    // sort and unique
    std::sort(colMapIndex.begin(), colMapIndex.end());
    auto ip = std::unique(colMapIndex.begin(), colMapIndex.end());
    colMapIndex.resize(std::distance(colMapIndex.begin(), ip));

    // create colMap
    Teuchos::RCP<TMAP> colMapRcp = Teuchos::rcp(
        new TMAP(Teuchos::OrdinalTraits<int>::invalid(), colMapIndex.data(), colMapIndex.size(), 0, commRcp));

    // convert columnIndices from global column index to local column index according to colMap
    auto &colmap = *colMapRcp;
#pragma omp parallel for
    for (int i = 0; i < colIndexNum; i++) {
        columnIndices[i] = colmap.getLocalElement(columnIndices[i]);
    }

    // printf("local number of cols: %d on rank %d\n", colMapRcp->getNodeNumElements(), commRcp->getRank());
    // printf("local number of rows: %d on rank %d\n", mobMapRcp->getNodeNumElements(), commRcp->getRank());
    // dumpTMAP(colMapRcp,"colMap");
    // dumpTMAP(mobMapRcp,"mobMap");

    // step 4, allocate the D^Trans matrix
    DMatTransRcp = Teuchos::rcp(new TCMAT(gammaMapRcp, colMapRcp, rowPointers, columnIndices, values));
    DMatTransRcp->fillComplete(mobMapRcp, gammaMapRcp); // domainMap, rangeMap

    // step 5, fill the delta0, gammaGuess, invKappa, conFlag vectors
    delta0Rcp = Teuchos::rcp(new TV(gammaMapRcp, true));
    invKappaRcp = Teuchos::rcp(new TV(gammaMapRcp, true));
    biFlagRcp = Teuchos::rcp(new TV(gammaMapRcp, true));
    gammaGuessRcp = Teuchos::rcp(new TV(gammaMapRcp, true));
    auto delta0 = delta0Rcp->getLocalView<Kokkos::HostSpace>();
    auto gammaGuess = gammaGuessRcp->getLocalView<Kokkos::HostSpace>();
    auto invKappa = invKappaRcp->getLocalView<Kokkos::HostSpace>();
    auto biFlag = biFlagRcp->getLocalView<Kokkos::HostSpace>();
    delta0Rcp->modify<Kokkos::HostSpace>();
    gammaGuessRcp->modify<Kokkos::HostSpace>();
    invKappaRcp->modify<Kokkos::HostSpace>();
    biFlagRcp->modify<Kokkos::HostSpace>();

#pragma omp parallel for num_threads(conQueNum)
    for (int que = 0; que < conQueNum; que++) {
        const auto &cQue = conPool[que];
        const int cIndexBase = conQueIndex[que];
        const int queSize = cQue.size();
        for (int j = 0; j < queSize; j++) {
            const auto &block = cQue[j];
            const auto idx = cIndexBase + j;
            delta0(idx, 0) = block.delta0;
            gammaGuess(idx, 0) = block.gamma;
            if (block.bilateral) {
                invKappa(idx, 0) = block.kappa > 0 ? 1 / block.kappa : 0;
                biFlag(idx, 0) = 1;
            }
        }
    }

    return 0;
}

int ConstraintCollector::buildConIndex(std::vector<int> &conQueSize, std::vector<int> &conQueIndex) const {
    const auto &conPool = *constraintPoolPtr;
    const int conQueNum = conPool.size();
    conQueSize.resize(conQueNum);
    conQueIndex.resize(conQueNum + 1);
    for (int i = 0; i < conQueNum; i++) {
        conQueSize[i] = conPool[i].size();
    }
    for (int i = 1; i < conQueNum + 1; i++) {
        conQueIndex[i] = conQueSize[i - 1] + conQueIndex[i - 1];
    }
    return 0;
}

int ConstraintCollector::writeBackGamma(const Teuchos::RCP<const TV> &gammaRcp) {
    auto &conPool = *constraintPoolPtr; // the constraint pool
    const int conQueNum = conPool.size();

    std::vector<int> conQueSize;
    std::vector<int> conQueIndex;
    buildConIndex(conQueSize, conQueIndex);

    auto gammaPtr = gammaRcp->getLocalView<Kokkos::HostSpace>();

#pragma omp parallel for num_threads(conQueNum)
    for (int i = 0; i < conQueNum; i++) {
        const int conQueSize = conPool[i].size();
        for (int j = 0; j < conQueSize; j++) {
            conPool[i][j].gamma = gammaPtr(conQueIndex[i] + j, 0);
            for (int k = 0; k < 9; k++) {
                conPool[i][j].stress[k] *= conPool[i][j].gamma;
            }
        }
    }

    return 0;
}
/**
 * @file SylinderSystem.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief System for sylinders
 * @version 1.0
 * @date 2018-12-13
 *
 * @copyright Copyright (c) 2018
 *
 */
#ifndef SYLINDERSYSTEM_HPP_
#define SYLINDERSYSTEM_HPP_

#include "Sylinder.hpp"
#include "SylinderConfig.hpp"
#include "SylinderNear.hpp"

#include "Boundary/Boundary.hpp"
#include "Constraint/ConstraintSolver.hpp"
#include "FDPS/particle_simulator.hpp"
#include "Trilinos/TpetraUtil.hpp"
#include "Trilinos/ZDD.hpp"
#include "Util/TRngPool.hpp"

#include <unordered_map>

// Custom hash for std::pair<int, int>
struct PairHash {
    std::size_t operator () (const std::pair<int, int>& p) const {
        std::hash<int> intHash;
        return intHash(p.first) ^ intHash(p.second);
    }
};

/**
 * @brief A collection of sylinders distributed to multiple MPI ranks.
 *
 */
class SylinderSystem {
    bool enableTimer = false;
    int snapID;                  ///< the current id of the snapshot file to be saved. sequentially numbered from 0
    int stepCount;               ///< timestep Count. sequentially numbered from 0
    unsigned int restartRngSeed; ///< parallel seed used by restarted simulations

    // FDPS stuff
    PS::DomainInfo dinfo; ///< domain size, boundary condition, and decomposition info
    void setDomainInfo();

    PS::ParticleSystem<Sylinder> sylinderContainer;        ///< sylinders
    std::unique_ptr<TreeSylinderNear> treeSylinderNearPtr; ///< short range interaction of sylinders
    int treeSylinderNumber;                                ///< the current max_glb number of treeSylinderNear
    void setTreeSylinder();

    std::unordered_multimap<int, int> endLinkMap;        ///< links prev,next
    std::unordered_multimap<int, int> endLinkReverseMap; ///< links next, prev
    std::unordered_multimap<int, int> centerLinkMap;        ///< links prev,next
    std::unordered_multimap<int, int> centerLinkReverseMap; ///< links next, prev
    std::unordered_multimap<int, std::pair<int,int>> triLinkMap;         ///< links center, (left, right)
    std::unordered_multimap<std::pair<int,int>, int, PairHash> triLinkReverseMap;  ///< links (left, right), center, 

    // Constraint stuff
    std::shared_ptr<ConstraintSolver> conSolverPtr;       ///< pointer to ConstraintSolver
    std::shared_ptr<ConstraintCollector> conCollectorPtr; ///<  pointer to ConstraintCollector
    Teuchos::RCP<const TV> forceUniRcp;                   ///< unilateral constraint force
    Teuchos::RCP<const TV> velocityUniRcp;                ///< unilateral constraint velocity
    Teuchos::RCP<const TV> forceBiRcp;                    ///< bilateral constraint force
    Teuchos::RCP<const TV> velocityBiRcp;                 ///< bilateral constraint velocity

    // computed without knowledge of constraints
    Teuchos::RCP<TV> forcePartNonBrownRcp;    ///< force specified by setForceNonBrown()
    Teuchos::RCP<TV> velocityPartNonBrownRcp; ///< velocity specified by setVelocityNonBrown()
    Teuchos::RCP<TV> velocityNonBrownRcp;     ///< \f$V_{NonBrown} = V_{part,NonBrown}+M F_{part,NonBrown}\f$
    Teuchos::RCP<TV> velocityBrownRcp;        ///< Brownian velocity, generated by calcBrown()
    Teuchos::RCP<TV> velocityNonConRcp;       ///< \f$V_{nc} = V_{Brown}+V_{NonBrown}\f$

    // MPI stuff
    std::shared_ptr<TRngPool> rngPoolPtr;      ///< TRngPool object for thread-safe random number generation
    Teuchos::RCP<const TCOMM> commRcp;         ///< TCOMM, set as a Teuchos::MpiComm object in constrctor
    Teuchos::RCP<TMAP> sylinderMapRcp;         ///< TMAP, contiguous and sequentially ordered 1 dof per sylinder
    Teuchos::RCP<TMAP> sylinderMobilityMapRcp; ///< TMAP, contiguous and sequentially ordered 6 dofs per sylinder
    Teuchos::RCP<TCMAT> mobilityMatrixRcp;     ///< block-diagonal mobility matrix
    Teuchos::RCP<TOP> mobilityOperatorRcp;     ///< full mobility operator (matrix-free), to be implemented

    // Data directory
    std::shared_ptr<ZDD<SylinderNearEP>> sylinderNearDataDirectoryPtr; ///< distributed data directory for sylinder data

    // internal utility functions
    /**
     * @brief generate initial configuration on rank 0 according to runConfig
     *
     */
    void setInitialFromConfig();

    /**
     * @brief set initial configuration as given in the (.dat) file
     *
     * The simBox and BC settings in runConfig are still used
     * @param filename
     */
    void setInitialFromFile(const std::string &filename);

    /**
     * @brief set linkMap from the .dat file
     *
     * Every mpi rank run this simultaneously to set linkMap from the same file
     * @param filename
     */
    void setEndLinkMapFromFile(const std::string &filename);
    void setCenterLinkMapFromFile(const std::string &filename);
    void setTriLinkMapFromFile(const std::string &filename);

    /**
     * @brief set initial configuration as given in the (.dat) file
     *
     * The simBox and BC settings in runConfig are still used
     * @param pvtpFileName
     */
    void setInitialFromVTKFile(const std::string &pvtpFileName);

    /**
     * @brief set initial configuration if runConfig.initCircularX is set
     *
     * This function move the position of all sylinders into a cylindrical tube fit in initBox
     */
    void setInitialCircularCrossSection();

    /**
     * @brief display the configuration on rank 0
     *
     */
    void showOnScreenRank0();

    /**
     * @brief update the sylinderMap and sylinderMobilityMap
     *
     * This function is called in prepareStep(), and no adding/removing/exchanging is allowed before runStep()
     */
    void updateSylinderMap(); ///< update sylindermap and sylinderMobilityMap

    /**
     * @brief write VTK parallel XML file into baseFolder
     *
     * @param baseFolder
     */
    void writeVTK(const std::string &baseFolder);

    /**
     * @brief write Ascii file controlled by FDPS into baseFolder
     *
     * @param baseFolder
     */
    void writeAscii(const std::string &baseFolder);

    /**
     * @brief write a txt file containing timestep and most recent pvtp filenames info into baseFolder
     *
     * @param baseFolder
     */
    void writeTimeStepInfo(const std::string &baseFolder);

    /**
     * @brief write a simple legacy VTK file for simBox
     *
     */
    void writeBox();

    /**
     * @brief Get orientation quaternion with givne px,py,pz
     *
     * component in [px,py,pz] out of range [-1,1] will be randomly generated
     * if all out of range [-1,1], a uniformly random orientation on sphere is generated
     * @param orient
     * @param px
     * @param py
     * @param pz
     * @param threadId openmp thread id for random number generation
     */
    void getOrient(Equatn &orient, const double px, const double py, const double pz, const int threadId);

    /**
     * @brief update the rank data field of sylinder
     *
     */
    void updateSylinderRank();

  public:
    SylinderConfig runConfig; ///< system configuration. Be careful if this is modified on the fly

    /**
     * @brief Construct a new SylinderSystem object
     *
     * initialize() should be called after this constructor
     */
    SylinderSystem() = default;

    /**
     * @brief Construct a new SylinderSystem object
     *
     * This constructor calls initialize() internally
     * @param configFile a yaml file for SylinderConfig
     * @param posFile initial configuration. use empty string ("") for no such file
     * @param argc command line argument
     * @param argv command line argument
     */
    SylinderSystem(const std::string &configFile, const std::string &posFile, int argc, char **argv);

    /**
     * @brief Construct a new SylinderSystem object
     *
     * This constructor calls initialize() internally
     * @param config SylinderConfig object
     * @param posFile initial configuration. use empty string ("") for no such file
     * @param argc command line argument
     * @param argv command line argument
     */
    SylinderSystem(const SylinderConfig &config, const std::string &posFile, int argc, char **argv);

    ~SylinderSystem() = default;

    // forbid copy
    SylinderSystem(const SylinderSystem &) = delete;
    SylinderSystem &operator=(const SylinderSystem &) = delete;

    /**
     * @brief initialize after an empty constructor
     *
     * @param config SylinderConfig object
     * @param posFile initial configuration. use empty string ("") for no such file
     * @param argc command line argument
     * @param argv command line argument
     */
    void initialize(const SylinderConfig &config, const std::string &posFile, int argc, char **argv);

    /**
     * @brief reinitialize from vtk files
     *
     * @param config SylinderConfig object
     * @param restartFile txt file containing timestep and most recent pvtp file names
     * @param argc command line argument
     * @param argv command line argument
     */
    void reinitialize(const SylinderConfig &config, const std::string &restartFile, int argc, char **argv,
                      bool eulerStep = true);

    /**
     * @brief enable the timer in step()
     *
     * @param value
     */
    void setTimer(bool value) { enableTimer = value; }

    /**
     * @brief compute axis-aligned bounding box of sylinders
     *
     * @param localLow
     * @param localHigh
     * @param globalLow
     * @param globalHigh
     */
    void calcBoundingBox(double localLow[3], double localHigh[3], double globalLow[3], double globalHigh[3]);

    /**
     * @brief compute domain decomposition by sampling sylinder distribution
     *
     * domain decomposition must be triggered when particle distribution significantly changes
     */
    void decomposeDomain();

    /**
     * @brief exchange between mpi ranks according to domain decomposition
     *
     * particle exchange must be triggered every timestep:
     */
    void exchangeSylinder();

    /**
     * one-step high level API
     */
    // get information
    /**
     * @brief Get sylinderContainer
     *
     * @return PS::ParticleSystem<Sylinder>&
     */
    const PS::ParticleSystem<Sylinder> &getContainer() { return sylinderContainer; }
    PS::ParticleSystem<Sylinder> &getContainerNonConst() { return sylinderContainer; }

    /**
     * @brief Get the DomainInfo object
     *
     * @return PS::DomainInfo&
     */
    const PS::DomainInfo &getDomainInfo() { return dinfo; }
    PS::DomainInfo &getDomainInfoNonConst() { return dinfo; }

    const std::unordered_multimap<int, int> &getEndLinkMap() { return endLinkMap; }
    const std::unordered_multimap<int, int> &getEndLinkReverseMap() { return endLinkReverseMap; }
    const std::unordered_multimap<int, int> &getCenterLinkMap() { return centerLinkMap; }
    const std::unordered_multimap<int, int> &getCenterLinkReverseMap() { return centerLinkReverseMap; }
    const std::unordered_multimap<int, std::pair<int,int>> &getTriLinkMap() { return triLinkMap; }
    const std::unordered_multimap<std::pair<int,int>, int, PairHash> &getTriLinkReverseMap() { return triLinkReverseMap; }

    /**
     * @brief Get the RngPoolPtr object
     *
     * @return std::shared_ptr<TRngPool>&
     */
    std::shared_ptr<TRngPool> &getRngPoolPtr() { return rngPoolPtr; }

    /**
     * @brief Get the CommRcp object
     *
     * @return Teuchos::RCP<const TCOMM>&
     */
    Teuchos::RCP<const TCOMM> &getCommRcp() { return commRcp; }

    /**
     * @brief prepare a step
     *
     * apply simBox boundary condition
     * decomposeDomain() for every 50 steps
     * exchangeSylinder() at every step
     * clear velocity
     * rebuild map
     * compute mobility matrix&operator
     * between prepareStep() and runStep(), sylinders should not be moved, added, or removed
     */
    void prepareStep();

    /**
     * @brief Set the (optional) forceNonBrownRcp
     *
     * This is optional.
     * this force is added to forceNonB
     * The computed mobility matrix will be applied to this force and the result is added to velNonB
     * @param forceNonBrown
     */
    void setForceNonBrown(const std::vector<double> &forceNonBrown);

    /**
     * @brief Set the (optional) velocityNonBrownRcp
     *
     * This is optional. The result is added to velNonB
     * @param velNonBrown
     */
    void setVelocityNonBrown(const std::vector<double> &velNonBrown);

    ConstraintBlockPool &getConstraintPoolNonConst() { return *(conCollectorPtr->constraintPoolPtr); };

    /**
     * @brief resolve collision with given nonBrownian motion and advance the
     * system configuration
     *
     * @param count_flag: add count step to index
     *
     */
    void runStep(bool count_flag = true);

    // These should run after runStep()
    /**
     * @brief add new Sylinders into the system from all ranks
     *
     * add new sylinders
     * 1. new gids will be randomly generated and assigned to each new sylinder
     * 2. added new sylinders will be appended to the local rank
     *
     * @param newSylinder list of new sylinders.
     * @return the generated new gids of the added new sylinders
     */
    std::vector<int> addNewSylinder(const std::vector<Sylinder> &newSylinder);

    /**
     * @brief add new links into the system from all ranks
     *
     * the newLink will be gathered from all ranks, remove duplication, and synchronized to the linkMap on all ranks
     *
     * @param newLink
     */
    void addNewEndLink(const std::vector<Link> &newEndLink);
    void addNewCenterLink(const std::vector<Link> &newCenterLink);

    /**
     * @brief calculate both Col and Bi stress
     *
     */
    void calcConStress();

    /**
     * @brief calculate polar and nematic order parameter
     *
     * The result is shown on screen
     */
    void calcOrderParameter();

    /**
     * @brief calculate volume fraction
     *
     */
    void calcVolFrac();

    /**
     * detailed low level API
     */
    /**
     * @brief apply periodic boundary condition
     *
     */
    void applyBoxBC();

    // compute non-collision velocity and mobility, before collision resolution
    /**
     * @brief calculate translational and rotational Brownian motion as specified in runConfig
     *
     * write back to sylinder.velBrown/omegaBrown
     */
    void calcVelocityBrown();

    /**
     * @brief calculate known velocity before collision resolution
     *
     * velocityNonCon = velocityBrown + velocityNonBrown + mobility * forceNonBrown
     * velocityNonBrown sums both the values set by setVelocityNonBrown() and directly written to
     * sylinder[i].velNonB/omegaNonB
     * write back to sylinder.velNonB/omegaNonB
     */
    void calcVelocityNonCon();

    /**
     * @brief sum vel = velNonB + velB + velCol + velBi
     *
     */
    void sumForceVelocity();

    /**
     * @brief calculate the mobility matrix (block diagonal)
     *
     */
    void calcMobMatrix();

    /**
     * @brief calculate the mobility operator (full-dense, matrix-free)
     *
     * TODO: to be implemented
     */
    void calcMobOperator();

    /**
     * @brief build the ZDD<SylinderNearEP> object
     *
     */
    void buildSylinderNearDataDirectory();

    /**
     * @brief Get the SylinderNearDataDirectory object
     *
     * @return std::shared_ptr<const ZDD<SylinderNearEP>>&
     */
    std::shared_ptr<ZDD<SylinderNearEP>> &getSylinderNearDataDirectory() { return sylinderNearDataDirectoryPtr; }

    // resolve constraints
    void collectPairCollision();        ///< collect pair collision constraints
    void collectBoundaryCollision();    ///< collect boundary collision constraints
    void collectEndLinkBilateral();     ///< setup link constraints
    void collectCenterLinkBilateral();  ///< setup link constraints
    void collectTriLinkBilateral();     ///< setup link constraints

    void resolveConstraints();           ///< resolve constraints
    void saveForceVelocityConstraints(); ///< write back to sylinder.velCol and velBi

    void stepEuler(); ///< Euler step update position and orientation, with both collision and non-collision velocity

    // write results
    std::string getCurrentResultFolder();          ///< get the current output folder path
    std::string getResultFolderWithID(int snapID); ///< get output folder path with snapID
    bool getIfWriteResultCurrentStep();            ///< check if the current step is writing (set by runConfig)
    int getSnapID() { return snapID; };            ///< get the (sequentially ordered) ID of current snapshot
    int getStepCount() { return stepCount; };      ///< get the (sequentially ordered) count of steps executed
    void writeResult();                            ///< write result regardless of runConfig

    // expose raw vectors and operators
    // non-constraint parts
    Teuchos::RCP<TV> getForcePartNonBrown() const { return forcePartNonBrownRcp; }
    Teuchos::RCP<TV> getVelocityPartNonBrown() const { return velocityPartNonBrownRcp; };
    Teuchos::RCP<TV> getVelocityNonBrown() const { return velocityNonBrownRcp; };
    Teuchos::RCP<TV> getVelocityBrown() const { return velocityBrownRcp; };
    Teuchos::RCP<TV> getVelocityNonCon() const { return velocityNonConRcp; };

    // constraint parts
    Teuchos::RCP<const TV> getForceUni() const { return forceUniRcp; };
    Teuchos::RCP<const TV> getVelocityUni() const { return velocityUniRcp; };
    Teuchos::RCP<const TV> getForceBi() const { return forceBiRcp; };
    Teuchos::RCP<const TV> getVelocityBi() const { return velocityBiRcp; };

    // mobility
    Teuchos::RCP<TCMAT> getMobMatrix() { return mobilityMatrixRcp; };
    Teuchos::RCP<TOP> getMobOperator() { return mobilityOperatorRcp; };

    // get information
    /**
     * @brief Get the local and global max gid for sylinders
     *
     * @return std::pair<int, int> [localMaxGid,globalMaxGid]
     */
    std::pair<int, int> getMaxGid();

    /**
     * @brief
     *
     * @param zeroOut zero out all timing info after printing out
     */
    void printTimingSummary(const bool zeroOut = true);
};

#endif
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

#include "Collision/CollisionSolver.hpp"
#include "FDPS/particle_simulator.hpp"
#include "Trilinos/TpetraUtil.hpp"
#include "Util/TRngPool.hpp"

/**
 * @brief A collection of sylinders distributed to multiple MPI ranks.
 *
 */
class SylinderSystem {
    int snapID;    ///< the current id of the snapshot file to be saved. sequentially numbered from 0
    int stepCount; ///< timestep Count. sequentially numbered from 0

    // FDPS stuff
    PS::DomainInfo dinfo; ///< domain size, boundary condition, and decomposition info
    void setDomainInfo();

    PS::ParticleSystem<Sylinder> sylinderContainer;        ///< sylinders
    std::unique_ptr<TreeSylinderNear> treeSylinderNearPtr; ///< short range interaction of sylinders
    int treeSylinderNumber;                                ///< the current max_glb number of treeSylinderNear
    void setTreeSylinder();

    // Collision stuff
    std::shared_ptr<CollisionSolver> collisionSolverPtr;       ///< pointer to CollisionSolver object
    std::shared_ptr<CollisionCollector> collisionCollectorPtr; ///<  pointer to CollisionCollector object
    Teuchos::RCP<TV> forceNonBrownRcp;                         ///< Non Brownian force, set by user
    Teuchos::RCP<TV> velocityNonBrownRcp;                      ///< Non Brownian velocity, set by user
    Teuchos::RCP<TV> velocityBrownRcp;                         ///< Brownian velocity, generated by calcBrown

    // used in collision solver
    Teuchos::RCP<TV> velocityKnownRcp; ///< \f$V_{known} = V_{Brown}+V_{NonBrown}+M F_{NonBrown}\f$
    Teuchos::RCP<TV> forceColRcp;      ///< collision force solution
    Teuchos::RCP<TV> velocityColRcp;   ///< collision velocity solution

    // MPI stuff
    std::shared_ptr<TRngPool> rngPoolPtr;      ///< TRngPool object for thread-safe random number generation
    Teuchos::RCP<const TCOMM> commRcp;         ///< TCOMM, set as a Teuchos::MpiComm object in constrctor
    Teuchos::RCP<TMAP> sylinderMapRcp;         ///< TMAP, contiguous and sequentially ordered 1 dof per sylinder
    Teuchos::RCP<TMAP> sylinderMobilityMapRcp; ///< TMAP, contiguous and sequentially ordered 6 dofs per sylinder
    Teuchos::RCP<TCMAT> mobilityMatrixRcp;     ///< block-diagonal mobility matrix
    Teuchos::RCP<TOP> mobilityOperatorRcp;     ///< full mobility operator (matrix-free), to be implemented

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
     * @brief write a simple legacy VTK file for simBox
     *
     */
    void writeBox();

    /**
     * @brief directly set the position of sylinders to non-overlap with wall
     *
     * used only for randomly generated initial configuration
     */
    void setPosWithWall();

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

    // /**
    //  * @brief Get random coordinate [x,y] in a circle with radius
    //  *
    //  * @param radius
    //  * @param x
    //  * @param y
    //  * @param threadId openmp thread id for random number generation
    //  */
    // void getRandPointInCircle(const double &radius, double &x, double &y, const int &threadId);

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
    PS::ParticleSystem<Sylinder> &getContainer() { return sylinderContainer; }

    PS::ParticleSystem<Sylinder> *getContainerPtr() { return &sylinderContainer; }

    /**
     * @brief Get the DomainInfo object
     *
     * @return PS::DomainInfo&
     */
    PS::DomainInfo &getDomainInfo() { return dinfo; }

    PS::DomainInfo *getDomainInfoPtr() { return &dinfo; }
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
     * The computed mobility matrix will be applied to this force and the result is added to velKnown
     * @param forceNonBrown
     */
    void setForceNonBrown(const std::vector<double> &forceNonBrown);
    /**
     * @brief Set the (optional) velocityNonBrownRcp
     *
     * This is optional. The result is added to velKnown
     * @param velNonBrown
     */
    void setVelocityNonBrown(const std::vector<double> &velNonBrown);
    /**
     * @brief resolve collision with given nonBrownian motion and advance the system configuration
     *
     */
    void runStep();

    // These should run after runStep()
    /**
     * @brief add new Sylinders into the system from every rank
     *
     * add new sylinders and assign new (unique) gid
     *
     * @param newSylinder
     */
    void addNewSylinder(std::vector<Sylinder> &newSylinder);
    /**
     * @brief calculate collision stress with solved collision problem
     *
     * The result is shown on screen
     */
    void calcColStress();
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
     * velocityKnown = velocityBrown + velocityNonBrown + mobility * forceNonBrown
     * velocityNonBrown sums both the values set by setVelocityNonBrown() and directly written to
     * sylinder[i].velNonB/omegaNonB
     * write back to sylinder.velNonB/omegaNonB
     */
    void calcVelocityKnown();

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

    // resolve collision
    void collectWallCollision();  ///< collect wall collision constraints
    void collectPairCollision();  ///< collect pair collision constraints
    void resolveCollision();      ///< resolve collision
    void saveVelocityCollision(); ///< write back to sylinder.velCol

    void stepEuler(); ///< Euler step update position and orientation, with both collision and non-collision velocity

    // write results
    std::string getCurrentResultFolder();     ///< get the current output folder path
    bool getIfWriteResultCurrentStep();       ///< check if the current step is writing (set by runConfig)
    int getSnapID() { return snapID; };       ///< get the (sequentially ordered) ID of current snapshot
    int getStepCount() { return stepCount; }; ///< get the (sequentially ordered) count of steps executed
    void writeResult();                       ///< write result regardless of runConfig

    // expose raw vectors and operators
    Teuchos::RCP<TV> getForceNonBrown() const { return forceNonBrownRcp; }
    Teuchos::RCP<TV> getVelocityNonBrown() const { return velocityNonBrownRcp; };
    Teuchos::RCP<TV> getVelocityBrown() const { return velocityBrownRcp; };
    Teuchos::RCP<TV> getVelocityKnown() const { return velocityKnownRcp; };
    Teuchos::RCP<TV> getForceCol() const { return forceColRcp; };
    Teuchos::RCP<TV> getVelocityCol() const { return velocityColRcp; };
    Teuchos::RCP<TCMAT> getMobMatrix() { return mobilityMatrixRcp; };
    Teuchos::RCP<TOP> getMobOperator() { return mobilityOperatorRcp; };

    // get information
    /**
     * @brief Get the local and global max gid for sylinders
     *
     * @return std::pair<int, int> [localMaxGid,globalMaxGid]
     */
    std::pair<int, int> getMaxGid();
};

#endif
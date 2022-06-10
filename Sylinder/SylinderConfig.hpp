/**
 * @file SylinderConfig.hpp
 * @author wenyan4work (wenyan4work@gmail.com)
 * @brief read configuration parameters from a yaml file
 * @version 1.0
 * @date 2018-12-13
 *
 * @copyright Copyright (c) 2018
 *
 */
#ifndef SYLINDERCONFIG_HPP_
#define SYLINDERCONFIG_HPP_

#include "Boundary/Boundary.hpp"
#include "Util/GeoCommon.h"

#include <iostream>
#include <memory>

/**
 * @brief read configuration parameters from a yaml file
 *
 */
class SylinderConfig {
  public:
    unsigned int rngSeed; ///< random number seed
    int logLevel;         ///< follows SPDLOG level enum, see Util/Logger.hpp for details
    int timerLevel = 0;   ///< how detailed the timer should be

    // domain setting
    double simBoxHigh[3];   ///< simulation box size
    double simBoxLow[3];    ///< simulation box size
    bool simBoxPBC[3];      ///< flag of true/false of periodic in that direction
    bool monolayer = false; ///< flag for simulating monolayer on x-y plane

    double initBoxHigh[3];      ///< initialize sylinders within this box
    double initBoxLow[3];       ///< initialize sylinders within this box
    double initOrient[3];       ///< initial orientation for each sylinder. >1 <-1 means random
    bool initCircularX = false; ///< set the initial cross-section as a circle in the yz-plane
    int initPreSteps = 100;     ///< number of initial pre steps to resolve potential collisions
    double thermEquilTime = 0;  ///< Time to equilibrate system before running crosslinking steps

    // physical constant
    double viscosity; ///< unit pN/(um^2 s), water ~ 0.0009
    double KBT;       ///< pN.um, 0.00411 at 300K
    double linkKappa; ///< pN/um stiffness of sylinder links
    double linkGap;   ///< um length of gap between sylinder links

    // sylinder settings
    bool sylinderFixed = false; ///< sylinders do not move
    int sylinderNumber;         ///< initial number of sylinders
    double sylinderLength;      ///< sylinder length (mean if sigma>0)
    double sylinderLengthSigma; ///< sylinder length lognormal distribution sigma
    double sylinderDiameter;    ///< sylinder diameter

    // collision radius and diameter
    double sylinderDiameterColRatio; ///< collision diameter = ratio * real diameter
    double sylinderLengthColRatio;   ///< collision length = ratio * real length
    double sylinderColBuf;           ///< threshold for recording possible collision

    // time stepping
    double dt;        ///< timestep size
    double timeTotal; ///< total simulation time
    double timeSnap;  ///< snapshot time. save one group of data for each snapshot

    // constraint solver
    double conResTol;    ///< constraint solver residual
    int conMaxIte;       ///< constraint solver maximum iteration
    int conSolverChoice; ///< choose a iterative solver. 0 for BBPGD, 1 for APGD, etc

    std::vector<std::shared_ptr<Boundary>> boundaryPtr;

    SylinderConfig() = default;
    SylinderConfig(std::string filename);
    ~SylinderConfig() = default;

    void dump() const;
};

#endif
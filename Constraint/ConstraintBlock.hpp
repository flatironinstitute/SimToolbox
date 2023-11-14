/**
 * @file ConstraintBlock.hpp
 * @author Wen Yan (wenyan4work@gmail.com)
 * @brief
 * @version 0.1
 * @date 2019-11-04
 *
 * @copyright Copyright (c) 2019
 *
 */
#ifndef CONSTRAINTBLOCK_HPP_
#define CONSTRAINTBLOCK_HPP_

#include "Util/EigenDef.hpp"
#include "Util/GeoCommon.h"
#include "Util/IOHelper.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <type_traits>
#include <vector>

/**
 * @brief collision constraint information block
 *
 * Each block stores the information for one collision constraint.
 * The blocks are collected by ConstraintCollector and then used to construct the sparse fcTrans matrix
 */
struct ConstraintBlock {
  public:
    double delta0 = 0;                    ///< constraint initial value
    double gamma = 0;                     ///< force magnitude, could be an initial guess
    int gidI = GEO_INVALID_INDEX;         ///< unique global ID of particle I
    int gidJ = GEO_INVALID_INDEX;         ///< unique global ID of particle J
    int gidK = GEO_INVALID_INDEX;         ///< unique global ID of particle K
    int globalIndexI = GEO_INVALID_INDEX; ///< global index of particle I
    int globalIndexJ = GEO_INVALID_INDEX; ///< global index of particle J
    int globalIndexK = GEO_INVALID_INDEX; ///< global index of particle K
    int gcid = GEO_INVALID_INDEX;         ///< unique global ID of constraint
    bool oneSide = false;                 ///< flag for one side constraint. body J does not appear in mobility matrix
    bool bilateral = false;               ///< if this is a bilateral constraint or not
    double kappa = 0;                     ///< spring constant. =0 means no spring
    double labI[3] = {0};                 ///< the labframe location of constraint on particle I
    double labJ[3] = {0};                 ///< the labframe location of constraint on particle J
    double labK[3] = {0};                 ///< the labframe location of constraint on particle K
    double unscaledForceComI[3] = {0};    ///< com force induced by this constraint on particle I for
                                          ///< unit constraint Lagrange multiplier gamma
    double unscaledForceComJ[3] = {0};    ///< com force induced by this constraint on particle J for
                                          ///< unit constraint Lagrange multiplier gamma
    double unscaledForceComK[3] = {0};    ///< com force induced by this constraint on particle K for
                                          ///< unit constraint Lagrange multiplier gamma
    double unscaledTorqueComI[3] = {0};   ///< com torque induced by this constraint on particle I for
                                          ///< unit constraint Lagrange multiplier gamma
    double unscaledTorqueComJ[3] = {0};   ///< com torque induced by this constraint on particle J for
                                          ///< unit constraint Lagrange multiplier gamma
    double unscaledTorqueComK[3] = {0};   ///< com torque induced by this constraint on particle K for
                                          ///< unit constraint Lagrange multiplier gamma
    double stress[9] = {0};              ///< virial stress induced by these constraints
    ///< stress 3x3 matrix (row-major) for unit constraint force gamma

    /**
     * @brief Construct a new empty collision block
     *
     */
    ConstraintBlock() = default;

    /**
     * @brief Construct a new ConstraintBlock object
     *
     */
    ConstraintBlock(double delta0_, double gamma_, 
                    int gidI_, int gidJ_, int globalIndexI_, int globalIndexJ_,
                    const double unscaledForceComI_[3], const double unscaledForceComJ_[3], 
                    const double unscaledTorqueComI_[3], const double unscaledTorqueComJ_[3],
                    const double labI_[3], const double labJ_[3], 
                    bool oneSide_, bool bilateral_, double kappa_, 
                    int gcid_ = GEO_INVALID_INDEX)
        : delta0(delta0_), gamma(gamma_), gidI(gidI_), gidJ(gidJ_), globalIndexI(globalIndexI_),
          globalIndexJ(globalIndexJ_), oneSide(oneSide_), bilateral(bilateral_), kappa(kappa_), gcid(gcid_) {
        for (int d = 0; d < 3; d++) {
            unscaledForceComI[d] = unscaledForceComI_[d];
            unscaledForceComJ[d] = unscaledForceComJ_[d];
            unscaledTorqueComI[d] = unscaledTorqueComI_[d];
            unscaledTorqueComJ[d] = unscaledTorqueComJ_[d];
            labI[d] = labI_[d];
            labJ[d] = labJ_[d];
        }
        std::fill(stress, stress + 9, 0);
    }

    /**
     * @brief Construct a new ConstraintBlock object
     *
     */
    ConstraintBlock(double delta0_, double gamma_, 
                    int gidI_, int gidJ_, int gidK_, 
                    int globalIndexI_, int globalIndexJ_, int globalIndexK_,
                    const double unscaledForceComI_[3], const double unscaledForceComJ_[3], const double unscaledForceComK_[3], 
                    const double unscaledTorqueComI_[3], const double unscaledTorqueComJ_[3], const double unscaledTorqueComK_[3],
                    const double labI_[3], const double labJ_[3], const double labK_[3], 
                    bool oneSide_, bool bilateral_, double kappa_, 
                    int gcid_ = GEO_INVALID_INDEX)
        : delta0(delta0_), gamma(gamma_), 
        gidI(gidI_), gidJ(gidJ_), gidK(gidK_), 
        globalIndexI(globalIndexI_), globalIndexJ(globalIndexJ_), globalIndexK(globalIndexK_),
        oneSide(oneSide_), bilateral(bilateral_), kappa(kappa_), gcid(gcid_) {
        for (int d = 0; d < 3; d++) {
            unscaledForceComI[d] = unscaledForceComI_[d];
            unscaledForceComJ[d] = unscaledForceComJ_[d];
            unscaledForceComK[d] = unscaledForceComK_[d];
            unscaledTorqueComI[d] = unscaledTorqueComI_[d];
            unscaledTorqueComJ[d] = unscaledTorqueComJ_[d];
            unscaledTorqueComK[d] = unscaledTorqueComK_[d];
            labI[d] = labI_[d];
            labJ[d] = labJ_[d];
            labK[d] = labK_[d];
        }
        std::fill(stress, stress + 9, 0);
    }

    void setStress(const Emat3 &stress_) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                stress[i * 3 + j] = stress_(i, j);
            }
        }
    }

    void setStress(const double *stress_) {
        for (int i = 0; i < 9; i++) {
            stress[i] = stress_[i];
        }
    }

    const double *getStress() const { return stress; }

    void getStress(Emat3 &stress_) const {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                stress_(i, j) = stress[i * 3 + j];
            }
        }
    }

    void reverseIJ() {
        std::swap(gidI, gidJ);
        std::swap(globalIndexI, globalIndexJ);
        for (int k = 0; k < 3; k++) {
            std::swap(unscaledForceComI[k], unscaledForceComJ[k]);
            std::swap(unscaledTorqueComI[k], unscaledTorqueComJ[k]);
            std::swap(labI[k], labJ[k]);
        }
    }
};

static_assert(std::is_trivially_copyable<ConstraintBlock>::value, "");
static_assert(std::is_default_constructible<ConstraintBlock>::value, "");

using ConstraintBlockQue = std::deque<ConstraintBlock>;      ///< a queue contains blocks collected by one thread
using ConstraintBlockPool = std::vector<ConstraintBlockQue>; ///< a pool contains queues on different threads

#endif
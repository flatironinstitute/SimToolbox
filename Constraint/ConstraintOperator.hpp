/**
 * @file ConstraintOperator.hpp
 * @author Wen Yan (wenyan4work@gmail.com)
 * @brief
 * @version 0.1
 * @date 2019-10-17
 *
 * @copyright Copyright (c) 2019
 *
 */
#ifndef CONSTRAINTOPERATOR_HPP_
#define CONSTRAINTOPERATOR_HPP_

#include "Trilinos/TpetraUtil.hpp"

#include <array>
#include <deque>
#include <vector>

/**
 * @brief Constraint Operator is a block matrix assembled from four blocks:
 *    [Du^T M Du       Du^T M Db           ]
 *    [Db^T M Du       Db^T M Db  +  K^{-1}]
 * The operator is applied on block vectors: [gammau; gammab]^T
 * Du^T, Db^T, M, and K^{-1} are explicitly constructed before constructing this object
 */
class ConstraintOperator : public TOP {
  public:
    /**
     * @brief Construct a new ConstraintOperator object
     *
     * @param mobOp
     * @param uniDuMat
     * @param biDbMat
     * @param invKappaDiagMat
     */
    ConstraintOperator(Teuchos::RCP<TOP> &mobOp, Teuchos::RCP<TCMAT> &uniDuMatTrans, Teuchos::RCP<TCMAT> &biDbMatTrans,
                       std::vector<double> &invKappaDiagMat);

    /**
     * @brief apply this operator, ensuring the block structure
     *
     * @param X
     * @param Y
     * @param mode if the operator should be applied as transposed
     * @param alpha
     * @param beta
     */
    void apply(const TMV &X, TMV &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

    /**
     * @brief Get the Domain Map object. interface required by Tpetra::Operator
     *
     * @return Teuchos::RCP<const TMAP>
     */
    Teuchos::RCP<const TMAP> getDomainMap() const;

    /**
     * @brief Get the Range Map object. interface required by Tpetra::Operator
     *
     * @return Teuchos::RCP<const TMAP>
     */
    Teuchos::RCP<const TMAP> getRangeMap() const;

    /**
     * @brief return if this operator can be applied as transposed. interface required by Tpetra::Operator
     *
     * @return true
     * @return false
     */
    bool hasTransposeApply() const { return false; }

    Teuchos::RCP<const TV> getForceUni() { return mobForceRcp->getVector(0); }
    Teuchos::RCP<const TV> getForceBi() { return mobForceRcp->getVector(1); }
    Teuchos::RCP<const TV> getVelUni() { return mobVelRcp->getVector(0); }
    Teuchos::RCP<const TV> getVelBi() { return mobVelRcp->getVector(1); }
    Teuchos::RCP<const TMAP> getUniBlockMap(){return gammaUniBlockMapRcp;}
    Teuchos::RCP<const TMAP> getBiBlockMap(){return gammaBiBlockMapRcp;}

  private:
    // comm
    Teuchos::RCP<const TCOMM> commRcp; ///< the mpi communicator
    // constant operators
    Teuchos::RCP<TOP> mobOpRcp;           ///< mobility matrix
    Teuchos::RCP<TCMAT> uniDuMatRcp;      ///< unilateral (collision) constraint geometry matrix D_c
    Teuchos::RCP<TCMAT> uniDuMatTransRcp; ///< explicit transpose of D_c
    Teuchos::RCP<TCMAT> biDbMatRcp;       ///< bilateral (spring) constraint geometry matrix D_b
    Teuchos::RCP<TCMAT> biDbMatTransRcp;  ///< explicit transpose of D_b
    std::vector<double> invKappaDiagMat;  ///< 1/h K^{-1} diagonal matrix, in std::vector format

    // maps
    Teuchos::RCP<const TMAP> mobMapRcp;           ///< map for mobility matrix. 6 DOF per obj
    Teuchos::RCP<const TMAP> gammaMapRcp;         ///< map for combined vector [gammau; gammab]^T
    Teuchos::RCP<const TMAP> gammaUniBlockMapRcp; ///< map for the rows in gammaMapRcp for the gammau block
    Teuchos::RCP<const TMAP> gammaBiBlockMapRcp;  ///< map for the rows in gammaMapRcp for the gammab block

    // working multivectors with 2 columns
    Teuchos::RCP<TMV> mobForceRcp; ///< force & torque vectors = [Du gamma_u, Db gamma_b]
    Teuchos::RCP<TMV> mobVelRcp;   ///< U & Omega vectors = [M Du gamma_u, M Db gamma_b]
    Teuchos::RCP<TMV> deltaUniRcp; ///< changes in constraint vectors = Du^T [M Du gamma_u, M Db gamma_b]
    Teuchos::RCP<TMV> deltaBiRcp;  ///< changes in constraint vectors = Db^T [M Du gamma_u, M Db gamma_b]

    /**
     * @brief build gammaMapRcp
     */
    void buildBlockMaps();
};

#endif
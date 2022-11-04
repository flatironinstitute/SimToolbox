#include "ConstraintSolver.hpp"
#include "Util/Logger.hpp"

void ConstraintSolver::setup(ConstraintCollector &conCollector_, Teuchos::RCP<TOP> &mobOpRcp_,
                             Teuchos::RCP<TV> &velncRcp_, double dt_) {
    reset();

    conCollector = conCollector_;

    dt = dt_;
    mobOpRcp = mobOpRcp_;
    velncRcp = velncRcp_;

    mobMapRcp = mobOpRcp->getDomainMap();

    conCollector.buildConstraintMatrixVector(mobMapRcp, DMatTransRcp, delta0Rcp, invKappaRcp, biFlagRcp, gammaRcp);

    delta0Rcp->scale(1.0 / dt);
    invKappaRcp->scale(1.0 / dt);

    deltancRcp = Teuchos::rcp(new TV(delta0Rcp->getMap(), true));
    DMatTransRcp->apply(*velncRcp, *deltancRcp);

    // the BCQP problem
    MOpRcp = Teuchos::rcp(new ConstraintOperator(mobOpRcp, DMatTransRcp, invKappaRcp));
    qRcp = Teuchos::rcp(new TV(delta0Rcp->getMap(), true));
    qRcp->update(1.0, *delta0Rcp, 1.0, *deltancRcp, 0.0);

    // result
    forcebRcp = Teuchos::rcp(new TV(mobMapRcp, true));
    forceuRcp = Teuchos::rcp(new TV(mobMapRcp, true));
    velbRcp = Teuchos::rcp(new TV(mobMapRcp, true));
    veluRcp = Teuchos::rcp(new TV(mobMapRcp, true));
}

void ConstraintSolver::reset() {
    setControlParams(1e-5, 1000000, 0);

    mobMapRcp.reset(); ///< distributed map for obj mobility. 6 dof per obj
    mobOpRcp.reset();  ///< mobility operator, 6 dof per obj to 6 dof per obj
    forceuRcp.reset(); ///< force vec, 6 dof per obj, due to unilateral constraints
    forcebRcp.reset(); ///< force vec, 6 dof per obj, due to bilateral constraints
    veluRcp.reset();   ///< velocity vec, 6 dof per obj. due to unilateral constraints
    velbRcp.reset();   ///< velocity vec, 6 dof per obj. due to bilateral constraints
    velncRcp.reset();  ///< the non-constraint velocity vel_nc

    // composite vectors and operators
    DMatTransRcp.reset(); ///< D^Trans matrix
    invKappaRcp.reset();  ///< K^{-1} diagonal matrix
    biFlagRcp.reset();    ///< bilateral flag vector
    delta0Rcp.reset();    ///< the current (geometric) delta vector delta_0 = [delta_0u ; delta_0b]
    deltancRcp.reset();   ///< delta_nc = [Du^Trans vel_nc,u ; Db^Trans vel_nc,b]

    // the constraint problem M gamma + q
    MOpRcp.reset();   ///< the operator of BCQP problem. M = [B,C;E,F]
    gammaRcp.reset(); ///< the unknown constraint force magnitude gamma = [gamma_u;gamma_b]
    qRcp.reset();     ///< the constant part of BCQP problem. q = delta_0 + delta_nc
}

void ConstraintSolver::solveConstraints() {
    const auto &commRcp = gammaRcp->getMap()->getComm();
    // solver
    BCQPSolver solver(MOpRcp, qRcp);
    spdlog::debug("solver constructed");

    // the bound of BCQP. 0 for gammau, unbound for gammab.
    Teuchos::RCP<TV> lbRcp = solver.getLowerBound();
    lbRcp->scale(-std::numeric_limits<double>::max() * .1, *biFlagRcp); // 0 if biFlag=0, -inf if biFlag=1
    spdlog::debug("bound constructed");

    // solve
    IteHistory history;
    switch (solverChoice) {
    case 0:
        solver.solveBBPGD(gammaRcp, res * (1.0/ dt), maxIte, history);
        break;
    case 1:
        solver.solveAPGD(gammaRcp, res * (1.0/ dt), maxIte, history);
        break;
    default:
        solver.solveBBPGD(gammaRcp, res * (1.0/ dt), maxIte, history);
        break;
    }

    for (auto it = history.begin(); it != history.end() - 1; it++) {
        auto &p = *it;
        spdlog::debug("RECORD: BCQP history {:g}, {:g}, {:g}, {:g}, {:g}, {:g}", p[0], p[1], p[2], p[3], p[4] * dt, 
                      p[5]);
    }

    auto &p = history.back();
    spdlog::info("RECORD: BCQP residue {:g}, {:g}, {:g}, {:g}, {:g}, {:g}", p[0], p[1], p[2], p[3], p[4] * dt, 
                p[5]);

    // calculate unilateral and bilateral vel/force with solution
    Teuchos::RCP<TCMAT> DMatRcp = MOpRcp->getDMat();
    // bilateral first
    Teuchos::RCP<TV> gammaBiRcp = Teuchos::rcp(new TV(gammaRcp->getMap(), true));
    gammaBiRcp->elementWiseMultiply(1.0, *gammaRcp, *biFlagRcp, 0.0);
    DMatRcp->apply(*gammaBiRcp, *forcebRcp);
    mobOpRcp->apply(*forcebRcp, *velbRcp);
    // unilateral second
    Teuchos::RCP<TV> forceRcp = MOpRcp->getForce();
    Teuchos::RCP<TV> velRcp = MOpRcp->getVel();
    forceuRcp->update(1.0, *forceRcp, -1.0, *forcebRcp, 0.0); // force_u = force - force_b
    veluRcp->update(1.0, *velRcp, -1.0, *velbRcp, 0.0);       // vel_u = vel - vel_b
}

void ConstraintSolver::writebackGamma() { conCollector.writeBackGamma(gammaRcp.getConst()); }
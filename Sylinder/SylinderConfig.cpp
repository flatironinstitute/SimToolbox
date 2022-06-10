#include "SylinderConfig.hpp"

#include "Util/Logger.hpp"
#include "Util/YamlHelper.hpp"

SylinderConfig::SylinderConfig(std::string filename) {

    YAML::Node config = YAML::LoadFile(filename);

    // required parameters
    readConfig(config, VARNAME(rngSeed), rngSeed, "");
    readConfig(config, VARNAME(simBoxLow), simBoxLow, 3, "");
    readConfig(config, VARNAME(simBoxHigh), simBoxHigh, 3, "");
    readConfig(config, VARNAME(simBoxPBC), simBoxPBC, 3, "");

    readConfig(config, VARNAME(viscosity), viscosity, "");
    readConfig(config, VARNAME(KBT), KBT, "");

    readConfig(config, VARNAME(sylinderNumber), sylinderNumber, "");
    readConfig(config, VARNAME(sylinderLength), sylinderLength, "");
    readConfig(config, VARNAME(sylinderDiameter), sylinderDiameter, "");

    readConfig(config, VARNAME(dt), dt, "");
    readConfig(config, VARNAME(timeTotal), timeTotal, "");
    readConfig(config, VARNAME(timeSnap), timeSnap, "");

    readConfig(config, VARNAME(conResTol), conResTol, "");
    readConfig(config, VARNAME(conMaxIte), conMaxIte, "");
    readConfig(config, VARNAME(conSolverChoice), conSolverChoice, "");

    // optional parameters
    logLevel = spdlog::level::info; // default to info
    readConfig(config, VARNAME(logLevel), logLevel, "", true);

    timerLevel = logLevel; // default to info
    readConfig(config, VARNAME(timerLevel), timerLevel, "", true);

    monolayer = false;
    readConfig(config, VARNAME(monolayer), monolayer, "", true);

    std::copy(simBoxLow, simBoxLow + 3, initBoxLow);
    std::copy(simBoxHigh, simBoxHigh + 3, initBoxHigh);
    readConfig(config, VARNAME(initBoxLow), initBoxLow, 3, "", true);
    readConfig(config, VARNAME(initBoxHigh), initBoxHigh, 3, "", true);

    initOrient[0] = initOrient[1] = initOrient[2] = 2;
    readConfig(config, VARNAME(initOrient), initOrient, 3, "", true);

    initCircularX = false;
    readConfig(config, VARNAME(initCircularX), initCircularX, "", true);

    initPreSteps = 100;
    readConfig(config, VARNAME(initPreSteps), initPreSteps, "", true);

    thermEquilTime = 0;
    readConfig(config, VARNAME(thermEquilTime), thermEquilTime, "", true);

    linkKappa = 100;
    linkGap = 0.01;
    readConfig(config, VARNAME(linkKappa), linkKappa, "", true);
    readConfig(config, VARNAME(linkGap), linkGap, "", true);

    sylinderFixed = false;
    readConfig(config, VARNAME(sylinderFixed), sylinderFixed, "", true);
    sylinderLengthSigma = -1;
    readConfig(config, VARNAME(sylinderLengthSigma), sylinderLengthSigma, "", true);
    sylinderDiameterColRatio = 1.0;
    readConfig(config, VARNAME(sylinderDiameterColRatio), sylinderDiameterColRatio, "", true);
    sylinderLengthColRatio = 1.0;
    readConfig(config, VARNAME(sylinderLengthColRatio), sylinderLengthColRatio, "", true);
    sylinderColBuf = 0.3;
    readConfig(config, VARNAME(sylinderColBuf), sylinderColBuf, "", true);

    boundaryPtr.clear();
    if (config["boundaries"]) {
        YAML::Node boundaries = config["boundaries"];
        for (const auto &b : boundaries) {
            std::string name = b["type"].as<std::string>();
            spdlog::debug(name);
            if (name == "wall") {
                boundaryPtr.push_back(std::make_shared<Wall>(b));
            } else if (name == "tube") {
                boundaryPtr.push_back(std::make_shared<Tube>(b));
            } else if (name == "sphere") {
                boundaryPtr.push_back(std::make_shared<SphereShell>(b));
            }
        }
    }
}

void SylinderConfig::dump() const {
    {
        printf("-------------------------------------------\n");
        printf("Run Setting: \n");
        printf("Random number seed: %d\n", rngSeed);
        printf("Log Level: %d\n", logLevel);
        printf("Timer Level: %d\n", timerLevel);
        printf("Simulation box Low: %g,%g,%g\n", simBoxLow[0], simBoxLow[1], simBoxLow[2]);
        printf("Simulation box High: %g,%g,%g\n", simBoxHigh[0], simBoxHigh[1], simBoxHigh[2]);
        printf("Periodicity: %d,%d,%d\n", simBoxPBC[0], simBoxPBC[1], simBoxPBC[2]);
        printf("Initialization box Low: %g,%g,%g\n", initBoxLow[0], initBoxLow[1], initBoxLow[2]);
        printf("Initialization box High: %g,%g,%g\n", initBoxHigh[0], initBoxHigh[1], initBoxHigh[2]);
        printf("Initialization orientation: %g,%g,%g\n", initOrient[0], initOrient[1], initOrient[2]);
        printf("Initialization circular cross: %d\n", initCircularX);
        printf("Initialization pre-steps for collision-resolution: %d\n", initPreSteps);
        printf("Time step size: %g\n", dt);
        printf("Total Time: %g\n", timeTotal);
        printf("Snap Time: %g\n", timeSnap);
        printf("-------------------------------------------\n");
    }
    {
        printf("-------------------------------------------\n");
        printf("For drag and collision: Sylinders with length < diameter are treated as spheres\n");
        printf("-------------------------------------------\n");
    }
    {
        printf("Physical setting: \n");
        printf("viscosity: %g\n", viscosity);
        printf("kBT: %g\n", KBT);
        printf("Link Kappa: %g\n", linkKappa);
        printf("Link Gap: %g\n", linkGap);
        printf("Sylinder Number: %d\n", sylinderNumber);
        printf("Sylinder Length: %g\n", sylinderLength);
        printf("Sylinder Length Sigma: %g\n", sylinderLengthSigma);
        printf("Sylinder Diameter: %g\n", sylinderDiameter);
        printf("Sylinder Length Collision Ratio: %g\n", sylinderLengthColRatio);
        printf("Sylinder Diameter Collision Ratio: %g\n", sylinderDiameterColRatio);
        printf("Sylinder Collision Buffer: %g\n", sylinderColBuf);
        printf("-------------------------------------------\n");
        printf("Constraint Solver Setting:\n");
        printf("Residual Tolerance: %g\n", conResTol);
        printf("Max Iteration: %d\n", conMaxIte);
        printf("Solver Choice: %d\n", conSolverChoice);
        printf("-------------------------------------------\n");
    }
    {
        for (const auto &b : boundaryPtr) {
            b->echo();
        }
    }
}
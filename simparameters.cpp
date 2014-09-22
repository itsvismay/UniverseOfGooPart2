#include "simparameters.h"

SimParameters::SimParameters()
{
    simRunning = false;
    constraint = CH_PENALTY_FORCE;
    timeStep = 0.001;
    NewtonMaxIters = 20;
    NewtonTolerance = 1e-8;
    penaltyStiffness = 1e5;

    activeForces = F_GRAVITY | F_SPRINGS | F_FLOOR | F_DAMPING | F_ELASTIC;
    gravityG = -9.8;
    springStiffness = 100;
    maxSpringStrain = 0.2;
    dampingStiffness = 1.0;

    clickMode = CM_ADDPARTICLE, CM_ADDSAW;
    connector = CT_SPRING;
    particleMass = 1.0;
    particleFixed = false;
    maxSpringDist = 0.25;
    sawRadius= 0.01;

    rodDensity = 2.0;
    rodStretch = 100.0;
    rodBend = 0.05;
    rodSegments = 5;
    ropeDensity = 2.0;
    ropeBend = 0.01;
    ropeSegments = 5;
}

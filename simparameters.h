#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

struct SimParameters
{
    SimParameters();

    enum ConstraintHandler {CH_PENALTY_FORCE, CH_STEP_PROJECT, CH_LAGRANGE};
    enum ClickMode {CM_ADDPARTICLE, CM_ADDSAW};
    enum ConnectorType{ CT_SPRING, CT_RIGID_ROD, CT_FLEXIBLE_ROD, CT_ROPE};
    const static int F_GRAVITY = 1;
    const static int F_SPRINGS = 2;
    const static int F_FLOOR   = 4;
    const static int F_DAMPING = 8;
    const static int F_ELASTIC = 16;

    bool simRunning;
    ConstraintHandler constraint;
    double timeStep;
    double NewtonTolerance;
    int NewtonMaxIters;
    double penaltyStiffness;

    int activeForces;
    double gravityG;
    double springStiffness;
    double maxSpringStrain;
    double dampingStiffness;

    ClickMode clickMode;
    ConnectorType connector;
    double particleMass;
    bool particleFixed;
    double maxSpringDist;
    double sawRadius;
    double rodDensity;
    double rodStretchStiffness;
    double rodBendingStiffness;
    int rodSegments;
    double ropeDensity;
    double ropeBend;
    double ropeSegments;

    //Extra
    bool cloudsOn;
    bool gameModeOn;

};

#endif // SIMPARAMETERS_H

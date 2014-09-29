#ifndef SIMULATION_H
#define SIMULATION_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <set>
#include <QMutex>

typedef Eigen::Triplet<double> Tr;

class SimParameters;

struct Particle
{
public:
    Particle(Eigen::Vector2d pos, double mass, bool isFixed, bool inert) : pos(pos), mass(mass), fixed(isFixed), inert(inert)
    {
        vel.setZero();
        prevpos = pos;
    }
    Particle(Eigen::Vector2d pos, double mass, bool isFixed) : pos(pos), mass(mass), fixed(isFixed)
    {
        vel.setZero();
        prevpos = pos;
    }

    Eigen::Vector2d pos;
    Eigen::Vector2d prevpos;
    Eigen::Vector2d vel;
    double mass;
    bool fixed;
    bool inert;
};


struct Connector
{
public:
    Connector(int p1, int p2, double restlen, double mass) : p1(p1), p2(p2), restlen(restlen), mass(mass) {}
    Connector(int p1, int p2, double restlen) : p1(p1), p2(p2), restlen(restlen) {
        mass = 0;
    }
    int p1;
    int p2;
    double restlen;
    double mass;
};

struct Rod: public Connector
{
public:
    Rod(int p1, int p2, double restlen) : Connector(p1, p2, restlen) {
        isHinge = false;
    }
    Rod(int p1, int p2, double restlen, double mass, bool isHinge) : Connector(p1, p2, restlen, mass), isHinge(isHinge) {
    }
    bool isHinge;
};

struct FlexibleRodHinge
{
public:
    FlexibleRodHinge(int s1, int s2, double stiffness) : s1(s1), s2(s2), stiffness(stiffness) {
    }

    int s1;
    int s2;
    double stiffness;
};

struct RopeHinge
{
public:
    RopeHinge(int s1, int s2, double stiffness) : s1(s1), s2(s2), stiffness(stiffness) {
    }

    int s1;
    int s2;
    double stiffness;
};

struct Spring: public Connector
{
public:
    Spring(int p1, int p2, double stiffness, double restlen) : Connector(p1, p2, restlen), stiffness(stiffness) {
        unsnappable = false;
    }
    Spring(int p1, int p2, double stiffness, double restlen, double mass, bool unsnappable) : Connector(p1, p2, restlen, mass), stiffness(stiffness), unsnappable(unsnappable) {}
    bool unsnappable;
    double stiffness;
};

struct Saw
{
public:
    Saw(Eigen::Vector2d pos, double radius) : pos(pos), radius(radius) {}

    Eigen::Vector2d pos;
    double radius;
};

struct Cloud
{
public:
    Cloud(Eigen::Vector2d pos1, Eigen::Vector2d pos2, Eigen::Vector2d pos3): pos1(pos1), pos2(pos2), pos3(pos3)
    {

    }
    Eigen::Vector2d pos1;
    Eigen::Vector2d pos2;
    Eigen::Vector2d pos3;
};

class Simulation
{
public:
    Simulation(const SimParameters &params);

    void addParticle(double x, double y);
    void addSaw(double x, double y);

    void takeSimulationStep();
    void render();
    void clearScene();

private:
    const SimParameters &params_;
    QMutex renderLock_;

    double time_;
    std::vector<Particle> particles_;
    std::vector<Spring> springs_;
    std::vector<FlexibleRodHinge> flexibleRodHinges_;
    std::vector<RopeHinge> ropeHinges_;
    std::vector<Rod> rods_;
    std::vector<Saw> saws_;
    std::vector<Cloud> clouds_;

    void buildConfiguration(Eigen::VectorXd &q, Eigen::VectorXd &qprev, Eigen::VectorXd &v);
    void unbuildConfiguration(const Eigen::VectorXd &q, const Eigen::VectorXd &v);

    void computeForceAndHessian(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, Eigen::SparseMatrix<double> &H);
    void processGravityForce(Eigen::VectorXd &F);
    void processElasticBendingForce(const Eigen::VectorXd &q, Eigen::VectorXd &F);
    void processSpringForce(const Eigen::VectorXd &q, Eigen::VectorXd &F, std::vector<Tr> &H);
    void processDampingForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Tr> &H);
    void processFloorForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Tr> &H);

    void processPenaltyForce(const Eigen::VectorXd &q, Eigen::VectorXd &F);

    void computeMassInverse(Eigen::SparseMatrix<double> &Minv);
    Eigen::SparseMatrix<double> computeGradGTranspose(const Eigen::VectorXd &q);

    void numericalIntegration(Eigen::VectorXd &q, Eigen::VectorXd &qprev, Eigen::VectorXd &v);
    void computeStepProject(Eigen::VectorXd &q, Eigen::VectorXd &oldq, Eigen::VectorXd &v);
    void computeLagrangeMultipliers(const Eigen::VectorXd &q, const Eigen::VectorXd &F, Eigen::VectorXd &v);

    void pruneOverstrainedSprings();
    void deleteSawedObjects();
    double ptSegmentDist(const Eigen::Vector2d &p, const Eigen::Vector2d &q1, const Eigen::Vector2d &q2);
    void detectSawedSprings(std::set<int> &springsToDelete, std::set<int> &hingesToDelete);
    void detectSawedRods(std::set<int> &rodsToDelete, std::set<int> &ropeHingesToDelete);
    void detectSawedParticles(std::set<int> &particlesToDelete);
};

#endif // SIMULATION_H

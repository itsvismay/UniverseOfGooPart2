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


struct Rod
{
    Rod(int p1, int p2, double restlen) : p1(p1), p2(p2), restlen(restlen) {}

    int p1;
    int p2;
    double restlen;
};

struct Hinge
{
public:
    Hinge(int s1, int s2, double stiffness) : s1(s1), s2(s2), stiffness(stiffness) {}

    int s1;
    int s2;
    double stiffness;
};

struct Spring
{
public:
    Spring(int p1, int p2, double stiffness, double restlen) : p1(p1), p2(p2), stiffness(stiffness), restlen(restlen) {
        mass = 0;
        unsnappable = false;
    }
    Spring(int p1, int p2, double stiffness, double restlen, double mass, bool unsnappable) : p1(p1), p2(p2), stiffness(stiffness), restlen(restlen), mass(mass), unsnappable(unsnappable) {}

    int p1;
    int p2;
    double stiffness;
    double restlen;
    double mass;
    bool unsnappable;
};

struct Saw
{
public:
    Saw(Eigen::Vector2d pos, double radius) : pos(pos), radius(radius) {}

    Eigen::Vector2d pos;
    double radius;
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
    std::vector<Hinge> hinges_;
    std::vector<Rod> rods_;
    std::vector<Saw> saws_;

    void buildConfiguration(Eigen::VectorXd &q, Eigen::VectorXd &qprev, Eigen::VectorXd &v);
    void unbuildConfiguration(const Eigen::VectorXd &q, const Eigen::VectorXd &v);

    void computeForceAndHessian(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, Eigen::SparseMatrix<double> &H);
    void processGravityForce(Eigen::VectorXd &F);
    void processSpringForce(const Eigen::VectorXd &q, Eigen::VectorXd &F, std::vector<Tr> &H);
    void processDampingForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Tr> &H);
    void processFloorForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Tr> &H);

    void processPenaltyForce(const Eigen::VectorXd &q, Eigen::VectorXd &F);

    void computeMassInverse(Eigen::SparseMatrix<double> &Minv);

    void numericalIntegration(Eigen::VectorXd &q, Eigen::VectorXd &qprev, Eigen::VectorXd &v);

    void pruneOverstrainedSprings();
    void deleteSawedObjects();
    double ptSegmentDist(const Eigen::Vector2d &p, const Eigen::Vector2d &q1, const Eigen::Vector2d &q2);
    void detectSawedSprings(std::set<int> &springsToDelete, std::set<int> &hingesToDelete);
    void detectSawedRods(std::set<int> &rodsToDelete);
    void detectSawedParticles(std::set<int> &particlesToDelete);
};

#endif // SIMULATION_H

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
};

struct Spring
{
public:
    Spring(int p1, int p2, double stiffness, double restlen) : p1(p1), p2(p2), stiffness(stiffness), restlen(restlen) {}

    int p1;
    int p2;
    double stiffness;
    double restlen;
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
    std::vector<Saw> saws_;

    void buildConfiguration(Eigen::VectorXd &q, Eigen::VectorXd &qprev, Eigen::VectorXd &v);
    void unbuildConfiguration(const Eigen::VectorXd &q, const Eigen::VectorXd &v);

    void computeForceAndHessian(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, Eigen::SparseMatrix<double> &H);
    void processGravityForce(Eigen::VectorXd &F);
    void processSpringForce(const Eigen::VectorXd &q, Eigen::VectorXd &F, std::vector<Tr> &H);
    void processDampingForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Tr> &H);
    void processFloorForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Tr> &H);

    void computeMassInverse(Eigen::SparseMatrix<double> &Minv);

    void numericalIntegration(Eigen::VectorXd &q, Eigen::VectorXd &qprev, Eigen::VectorXd &v);

    void pruneOverstrainedSprings();
    void deleteSawedObjects();
    double ptSegmentDist(const Eigen::Vector2d &p, const Eigen::Vector2d &q1, const Eigen::Vector2d &q2);
    void detectSawedSprings(std::set<int> &springsToDelete);
    void detectSawedParticles(std::set<int> &particlesToDelete);
};

#endif // SIMULATION_H

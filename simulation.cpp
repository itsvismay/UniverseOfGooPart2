#include "simulation.h"
#include <QGLWidget>
#include "simparameters.h"
#include <iostream>

const double PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

Simulation::Simulation(const SimParameters &params) : params_(params), time_(0)
{
}

void Simulation::render()
{
    double baseradius = 0.02;
    double pulsefactor = 0.1;
    double pulsespeed = 50.0;

    int sawteeth = 20;
    double sawdepth = 0.1;
    double sawangspeed = 10.0;

    double baselinewidth = 0.5;

    int numcirclewedges = 20;

    if(params_.activeForces & SimParameters::F_FLOOR)
    {
        glBegin(GL_TRIANGLES);
        {
            glColor3f(0.3, 1.0, 0.3);

            glVertex2f(-1, -0.5);
            glVertex2f(1, -0.5);
            glVertex2f(-1, -1);

            glVertex2f(-1, -1);
            glVertex2f(1, -0.5);
            glVertex2f(1, -1);
        }
        glEnd();
    }

    renderLock_.lock();
    {
        for(vector<Spring>::iterator it = springs_.begin(); it != springs_.end(); ++it)
        {
            glColor3f(0.0, 0.0, 1.0);
            Vector2d sourcepos = particles_[it->p1].pos;
            Vector2d destpos   = particles_[it->p2].pos;

            double dist = (sourcepos-destpos).norm();

            glLineWidth(baselinewidth/dist);

            glBegin(GL_LINES);
            glVertex2f(sourcepos[0], sourcepos[1]);
            glVertex2f(destpos[0], destpos[1]);
            glEnd();
        }

        for(vector<Particle>::iterator it = particles_.begin(); it != particles_.end(); ++it)
        {
            double radius = baseradius*sqrt(it->mass);
            radius *= (1.0 + pulsefactor*sin(pulsespeed*time_));

            glColor3f(0,0,0);

            if(it->fixed)
            {
                radius = baseradius;
                glColor3f(1.0,0,0);
            }

            glBegin(GL_TRIANGLE_FAN);
            {
                glVertex2f(it->pos[0], it->pos[1]);
                for(int i=0; i<=numcirclewedges; i++)
                {
                    glVertex2f(it->pos[0] + radius * cos(2*PI*i/numcirclewedges),
                               it->pos[1] + radius * sin(2*PI*i/numcirclewedges));
                }
            }
            glEnd();
        }

        for(vector<Saw>::iterator it = saws_.begin(); it != saws_.end(); ++it)
        {
            double outerradius = it->radius;
            double innerradius = (1.0-sawdepth)*outerradius;

            glColor3f(0.5,0.5,0.5);

            glBegin(GL_TRIANGLE_FAN);
            {
                glVertex2f(it->pos[0], it->pos[1]);
                int spokes = 2*sawteeth;
                for(int i=0; i<=spokes; i++)
                {
                    double radius = (i%2==0) ? innerradius : outerradius;
                    glVertex2f(it->pos[0] + radius * cos(2*PI*i/spokes + sawangspeed*time_),
                               it->pos[1] + radius * sin(2*PI*i/spokes + sawangspeed*time_));
                }
            }
            glEnd();
        }
    }
    renderLock_.unlock();
}

void Simulation::takeSimulationStep()
{
    VectorXd q, qprev, v;
    buildConfiguration(q, qprev, v);
    numericalIntegration(q, qprev, v);
    unbuildConfiguration(q, v);

    pruneOverstrainedSprings();
    deleteSawedObjects();
    time_ += params_.timeStep;
}

void Simulation::addParticle(double x, double y)
{
    renderLock_.lock();
    {
        Vector2d newpos(x,y);
        for(int i=0; i<(int)particles_.size(); i++)
        {
            Vector2d pos = particles_[i].pos;
            double dist = (pos-newpos).norm();
            if(dist <= params_.maxSpringDist)
            {
                springs_.push_back(Spring(particles_.size(), i, params_.springStiffness/dist, dist));
            }
        }

        double mass = params_.particleMass;
        if(params_.particleFixed)
            mass = std::numeric_limits<double>::infinity();

        particles_.push_back(Particle(newpos, mass, params_.particleFixed));
    }
    renderLock_.unlock();
}

void Simulation::addSaw(double x, double y)
{
    renderLock_.lock();
        saws_.push_back(Saw(Vector2d(x,y), params_.sawRadius));
    renderLock_.unlock();
}

void Simulation::clearScene()
{
    renderLock_.lock();
    {
        particles_.clear();
        springs_.clear();
        saws_.clear();
    }
    renderLock_.unlock();
}

void Simulation::buildConfiguration(VectorXd &q, VectorXd &qprev, VectorXd &v)
{
    int ndofs = 2*particles_.size();
    q.resize(ndofs);
    qprev.resize(ndofs);
    v.resize(ndofs);

    for(int i=0; i<(int)particles_.size(); i++)
    {
        q.segment<2>(2*i) = particles_[i].pos;
        qprev.segment<2>(2*i) = particles_[i].prevpos;
        v.segment<2>(2*i) = particles_[i].vel;
    }
}

void Simulation::unbuildConfiguration(const VectorXd &q, const VectorXd &v)
{
    int ndofs = q.size();
    assert(ndofs == int(2*particles_.size()));
    for(int i=0; i<ndofs/2; i++)
    {
        particles_[i].prevpos = particles_[i].pos;
        particles_[i].pos = q.segment<2>(2*i);
        particles_[i].vel = v.segment<2>(2*i);
    }
}

void Simulation::computeForceAndHessian(const VectorXd &q, const VectorXd &qprev, Eigen::VectorXd &F, SparseMatrix<double> &H)
{
    F.resize(q.size());
    F.setZero();
    H.resize(q.size(), q.size());
    H.setZero();

    vector<Tr> Hcoeffs;
    if(params_.activeForces & SimParameters::F_GRAVITY)
        processGravityForce(F);
    if(params_.activeForces & SimParameters::F_SPRINGS)
        processSpringForce(q, F, Hcoeffs);
    if(params_.activeForces & SimParameters::F_DAMPING)
        processDampingForce(q, qprev, F, Hcoeffs);
    if(params_.activeForces & SimParameters::F_FLOOR)
        processFloorForce(q, qprev, F, Hcoeffs);

    H.setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());
}

void Simulation::processGravityForce(VectorXd &F)
{
    int nparticles = (int)particles_.size();
    for(int i=0; i<nparticles; i++)
    {
        if(!particles_[i].fixed)
        {
            F[2*i+1] += params_.gravityG*particles_[i].mass;
        }
    }
}

void Simulation::processSpringForce(const VectorXd &q, VectorXd &F, std::vector<Tr> &H)
{
    int nsprings = (int)springs_.size();

    for(int i=0; i<nsprings; i++)
    {
        Vector2d p1 = q.segment<2>(2*springs_[i].p1);
        Vector2d p2 = q.segment<2>(2*springs_[i].p2);
        double dist = (p2-p1).norm();
        Vector2d localF = springs_[i].stiffness*(dist-springs_[i].restlen)/dist * (p2-p1);
        F.segment<2>(2*springs_[i].p1) += localF;
        F.segment<2>(2*springs_[i].p2) -= localF;

        Matrix2d I;
        I << 1, 0, 0, 1;
        Matrix2d localH = springs_[i].stiffness * (1.0 - springs_[i].restlen/dist)*I;
        localH += springs_[i].stiffness*springs_[i].restlen*(p2-p1)*(p2-p1).transpose()/dist/dist/dist;

        for(int j=0; j<2; j++)
            for(int k=0; k<2;k++)
            {
                H.push_back(Tr(2*springs_[i].p1+j, 2*springs_[i].p1+k, localH.coeff(j,k)));
                H.push_back(Tr(2*springs_[i].p2+j, 2*springs_[i].p2+k, localH.coeff(j,k)));
                H.push_back(Tr(2*springs_[i].p1+j, 2*springs_[i].p2+k, -localH.coeff(j,k)));
                H.push_back(Tr(2*springs_[i].p2+j, 2*springs_[i].p1+k, -localH.coeff(j,k)));
            }
    }
}

void Simulation::processDampingForce(const VectorXd &q, const VectorXd &qprev, VectorXd &F, std::vector<Tr> &H)
{
    int nsprings = (int)springs_.size();

    for(int i=0; i<nsprings; i++)
    {
        Vector2d p1 = q.segment<2>(2*springs_[i].p1);
        Vector2d p2 = q.segment<2>(2*springs_[i].p2);
        Vector2d p1prev = qprev.segment<2>(2*springs_[i].p1);
        Vector2d p2prev = qprev.segment<2>(2*springs_[i].p2);

        Vector2d relvel = (p2 - p2prev)/params_.timeStep - (p1 - p1prev)/params_.timeStep;
        Vector2d localF = params_.dampingStiffness*relvel;
        F.segment<2>(2*springs_[i].p1) += localF;
        F.segment<2>(2*springs_[i].p2) -= localF;

        Matrix2d I;
        I << 1, 0, 0, 1;
        Matrix2d localH = params_.dampingStiffness*I/params_.timeStep;

        for(int j=0; j<2; j++)
            for(int k=0; k<2;k++)
            {
                H.push_back(Tr(2*springs_[i].p1+j, 2*springs_[i].p1+k, localH.coeff(j,k)));
                H.push_back(Tr(2*springs_[i].p2+j, 2*springs_[i].p2+k, localH.coeff(j,k)));
                H.push_back(Tr(2*springs_[i].p1+j, 2*springs_[i].p2+k, -localH.coeff(j,k)));
                H.push_back(Tr(2*springs_[i].p2+j, 2*springs_[i].p1+k, -localH.coeff(j,k)));
            }
    }
}

void Simulation::processFloorForce(const VectorXd &q, const VectorXd &qprev, VectorXd &F, std::vector<Tr> &H)
{
    int nparticles = particles_.size();

    double basestiffness = 10000;
    double basedrag = 1000.0;

    for(int i=0; i<nparticles; i++)
    {
        if(q[2*i+1] < -0.5 && ! particles_[i].fixed)
        {
            double vel = (q[2*i+1]-qprev[2*i+1])/params_.timeStep;
            double dist = -0.5 - q[2*i+1];

            F[2*i+1] += basestiffness*dist - basedrag*dist*vel;

            H.push_back(Tr(2*i+1, 2*i+1, basestiffness
                           - 0.5*basedrag/params_.timeStep
                           + basedrag*qprev[2*i+1]/params_.timeStep
                        - 2.0*basedrag*q[2*i+1]/params_.timeStep));
        }
    }
}

void Simulation::computeMassInverse(Eigen::SparseMatrix<double> &Minv)
{
    int ndofs = 2*int(particles_.size());

    Minv.resize(ndofs, ndofs);
    Minv.setZero();

    vector<Tr> Minvcoeffs;
    for(int i=0; i<ndofs/2; i++)
    {
        Minvcoeffs.push_back(Tr(2*i,   2*i,   1.0/particles_[i].mass));
        Minvcoeffs.push_back(Tr(2*i+1, 2*i+1, 1.0/particles_[i].mass));
    }

    Minv.setFromTriplets(Minvcoeffs.begin(), Minvcoeffs.end());
}

void Simulation::numericalIntegration(VectorXd &q, VectorXd &qprev, VectorXd &v)
{
//    VectorXd F;
//    SparseMatrix<double> H;
//    SparseMatrix<double> Minv;

//    computeMassInverse(Minv);

//    switch(params_.integrator)
//    {
//    case SimParameters::TI_EXPLICIT_EULER:
//    {
//        computeForceAndHessian(q, qprev, F, H);
//        q += params_.timeStep*v;
//        v += params_.timeStep*Minv*F;
//        break;
//    }
//    case SimParameters::TI_VELOCITY_VERLET:
//    {
//        VectorXd oldq = q;
//        q += params_.timeStep*v;
//        computeForceAndHessian(q, oldq, F, H);
//        v += params_.timeStep*Minv*F;
//        break;
//    }
//    case SimParameters::TI_IMPLICIT_EULER:
//    {
//        VectorXd guess = q;
//        int iter=0;
//        for(iter = 0; iter < params_.NewtonMaxIters; iter++)
//        {
//            computeForceAndHessian(guess, q, F, H);
//            VectorXd fval = guess - q - params_.timeStep*v - params_.timeStep*params_.timeStep*Minv*F;
//            if(fval.norm() < params_.NewtonTolerance)
//            {
//                break;
//            }
//            SparseMatrix<double> I(q.size(), q.size());
//            I.setIdentity();
//            SparseMatrix<double> Hf = I + params_.timeStep*params_.timeStep*Minv*H;
//            BiCGSTAB<SparseMatrix<double> > solver;
//            solver.compute(Hf);
//            VectorXd deltaguess = solver.solve(-fval);
//            guess += deltaguess;
//        }
//        v = (guess-q)/params_.timeStep;
//        q = guess;
//        break;
//    }
//    case SimParameters::TI_IMPLICIT_MIDPOINT:
//    {
//        VectorXd guess = q;
//        for(int iter = 0; iter < params_.NewtonMaxIters; iter++)
//        {
//            computeForceAndHessian(0.5*(guess+q), 0.5*(q+qprev), F, H);
//            VectorXd fval = guess - q - params_.timeStep*v - params_.timeStep*params_.timeStep*Minv*F/2.0;
//            if(fval.norm() < params_.NewtonTolerance)
//                break;
//            SparseMatrix<double> I(q.size(), q.size());
//            I.setIdentity();
//            SparseMatrix<double> Hf = I + params_.timeStep*params_.timeStep*Minv*H/4.0;
//            BiCGSTAB<SparseMatrix<double> > solver;
//            solver.compute(Hf);
//            VectorXd deltaguess = solver.solve(-fval);
//            guess += deltaguess;
//        }
//        v = 2.0*(guess-q)/params_.timeStep - v;
//        q = guess;
//        break;
//    }
//    }
}

void Simulation::pruneOverstrainedSprings()
{
    int nsprings = springs_.size();

    vector<int> toremove;
    for(int i=0; i<nsprings; i++)
    {
        Vector2d srcpos = particles_[springs_[i].p1].pos;
        Vector2d dstpos = particles_[springs_[i].p2].pos;
        double dist = (dstpos-srcpos).norm();

        double strain = (dist - springs_[i].restlen)/springs_[i].restlen;
        if(strain > params_.maxSpringStrain)
            toremove.push_back(i);
    }

    renderLock_.lock();
    {
        for(vector<int>::reverse_iterator it = toremove.rbegin(); it != toremove.rend(); ++it)
            springs_.erase(springs_.begin() + *it);
    }
    renderLock_.unlock();
}

double Simulation::ptSegmentDist(const Vector2d &p, const Vector2d &q1, const Vector2d &q2)
{
    double t = (p-q1).dot(q2-q1) / (q2-q1).dot(q2-q1);
    double linedistsq = (q1 + t*(q2-q1) - p).squaredNorm();
    double q1dist = (p-q1).squaredNorm();
    double q2dist = (p-q2).squaredNorm();
    double mindistsq = min(linedistsq, min(q1dist, q2dist));
    return sqrt(mindistsq);
}

void Simulation::detectSawedSprings(std::set<int> &springsToDelete)
{
    for(int i=0; i<(int)springs_.size(); i++)
    {
        Vector2d pos1 = particles_[springs_[i].p1].pos;
        Vector2d pos2 = particles_[springs_[i].p2].pos;
        double maxx = max(pos1[0], pos2[0]);
        double minx = min(pos1[0], pos2[0]);
        double maxy = max(pos1[1], pos2[1]);
        double miny = min(pos1[1], pos2[1]);
        for(vector<Saw>::iterator saw = saws_.begin(); saw != saws_.end(); ++saw)
        {
            Vector2d sawpos = saw->pos;
            double sawr = saw->radius;

            if(sawpos[0] - sawr > maxx || sawpos[0] + sawr < minx || sawpos[1] - sawr > maxy || sawpos[1] + sawr < miny)
                continue;

            double sawspringdist = ptSegmentDist(sawpos, pos1, pos2);
            if(sawspringdist <= sawr)
            {
                springsToDelete.insert(i);
                break;
            }
        }
    }
}

void Simulation::detectSawedParticles(std::set<int> &particlesToDelete)
{
    for(int i=0; i<(int)particles_.size(); i++)
    {
        Vector2d partpos = particles_[i].pos;

        if(fabs(partpos[0]) > 2 || fabs(partpos[1]) > 2)
        {
            particlesToDelete.insert(i);
            break;
        }

        for(vector<Saw>::iterator it = saws_.begin(); it != saws_.end(); ++it)
        {
            Vector2d sawpos = it->pos;
            double sqdist = (sawpos-partpos).squaredNorm();
            if(sqdist < it->radius*it->radius)
            {
                particlesToDelete.insert(i);
                break;
            }
        }
    }
}

void Simulation::deleteSawedObjects()
{
    set<int> particlestodelete;
    set<int> springstodelete;
    detectSawedParticles(particlestodelete);
    detectSawedSprings(springstodelete);

    vector<Particle> newparticles;
    vector<Spring> newsprings;
    vector<int> remainingparticlemap;

    if(!particlestodelete.empty())
    {
        for(int i=0; i<(int)springs_.size(); i++)
        {
            if(particlestodelete.count(springs_[i].p1) || particlestodelete.count(springs_[i].p2))
                springstodelete.insert(i);
        }

        for(int i=0; i<(int)particles_.size(); i++)
        {
            if(particlestodelete.count(i) == 0)
            {
                remainingparticlemap.push_back(newparticles.size());
                newparticles.push_back(particles_[i]);
            }
            else
                remainingparticlemap.push_back(-1);
        }
    }
    if(!springstodelete.empty())
    {
        for(int i=0; i<(int)springs_.size(); i++)
        {
            if(springstodelete.count(i) == 0)
            {
                newsprings.push_back(springs_[i]);
            }
        }
    }

    if(!springstodelete.empty() || !particlestodelete.empty())
    {
        renderLock_.lock();
        {
            if(!springstodelete.empty())
                springs_ = newsprings;
            if(!particlestodelete.empty())
            {
                particles_ = newparticles;
                for(vector<Spring>::iterator it = springs_.begin(); it != springs_.end(); ++it)
                {
                    it->p1 = remainingparticlemap[it->p1];
                    it->p2 = remainingparticlemap[it->p2];
                }
            }

        }
        renderLock_.unlock();
    }
}

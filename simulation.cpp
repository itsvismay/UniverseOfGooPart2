#include "simulation.h"
#include <QGLWidget>
#include "simparameters.h"
#include <iostream>
#include <Eigen/Geometry>

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
        // Rendering Rods
        for(vector<Rod>::iterator it = rods_.begin(); it != rods_.end(); ++it)
        {
            glColor3f(0.75, 0.0, 0.75);
            Vector2d sourcepos = particles_[it->p1].pos;
            Vector2d destpos   = particles_[it->p2].pos;

//            double dist = (sourcepos-destpos).norm();

            glLineWidth(8);

            glBegin(GL_LINES);
            glVertex2f(sourcepos[0], sourcepos[1]);
            glVertex2f(destpos[0], destpos[1]);
            glEnd();
        }
        //Rendering Springs
        for(vector<Spring>::iterator it = springs_.begin(); it != springs_.end(); ++it)
        {
            glColor3f(0.0, 0.0, 1.0);
            Vector2d sourcepos = particles_[it->p1].pos;
            Vector2d destpos   = particles_[it->p2].pos;

            double dist = (sourcepos-destpos).norm();

            glLineWidth(baselinewidth/dist);

            if (it->unsnappable)
            {
                glColor3f(0.752941, 0.752941, 0.752941);
                glLineWidth(14);
            }

            glBegin(GL_LINES);
            glVertex2f(sourcepos[0], sourcepos[1]);
            glVertex2f(destpos[0], destpos[1]);
            glEnd();
        }
        //Rendering Particles
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
        //Rendering Saws
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
        Vector2d newParticlePos(x,y);
        double particleMass = params_.particleMass;
        int newParticleIndex = particles_.size();
        if (params_.connector == params_.CT_FLEXIBLE_ROD)
        {
            particles_.push_back(Particle(newParticlePos, particleMass, params_.particleFixed, false));
        }
        for(int i=0; i<(int)particles_.size(); i++)
        {
            if (i == newParticleIndex)
            {
                continue;
            }
            Vector2d pos = particles_[i].pos;
            double dist = (pos-newParticlePos).norm();
            if(!particles_[i].inert && dist <= params_.maxSpringDist)
            {
                if (params_.connector == params_.CT_RIGID_ROD)
                {
                    rods_.push_back(Rod(particles_.size(), i, dist));
                }
                else if (params_.connector == params_.CT_SPRING)
                {
                    springs_.push_back(Spring(particles_.size(), i, params_.springStiffness/dist, dist));
                }
                else if (params_.connector == params_.CT_FLEXIBLE_ROD)
                {
                    int rodSegs = params_.rodSegments;
                    if (rodSegs <= 1)
                    {
                        rodSegs = 2;
                    }
                    double springLength = dist/rodSegs;
                    double springMass = params_.rodDensity * springLength;
                    double springStiffness = params_.rodStretchStiffness/springLength;
                    double hingeStiffness = 0;
                    Vector2d unitVector = (newParticlePos - pos)/dist;
                    Vector2d distanceToMove = unitVector * (dist/rodSegs);
                    Vector2d newInertParticlePos = (distanceToMove * 1) + pos;
                    springs_.push_back(Spring(particles_.size(), i, springStiffness, springLength, springMass, true));
                    particles_.push_back(Particle(newInertParticlePos, springMass, false, true));
                    particles_[i].mass += springMass/2;
                    int j;
                    for (j=2; j<=rodSegs - 1; j++)
                    {
                        newInertParticlePos = (distanceToMove * j) + pos;
                        springs_.push_back(Spring(particles_.size(), particles_.size() - 1, springStiffness, springLength, springMass, true));
                        particles_.push_back(Particle(newInertParticlePos, springMass, false, true));
                        hingeStiffness = (params_.rodBendingStiffness * 2)/(springs_[springs_.size() - 1].restlen + springs_[springs_.size() - 2].restlen);
                        hinges_.push_back(Hinge(springs_.size() - 1, springs_.size() - 2, hingeStiffness));
                    }
                    springs_.push_back(Spring(newParticleIndex, particles_.size() - 1, springStiffness, springLength, springMass, true));
                    hingeStiffness = (params_.rodBendingStiffness * 2)/(springs_[springs_.size() - 1].restlen + springs_[springs_.size() - 2].restlen);
                    hinges_.push_back(Hinge(springs_.size() - 1, springs_.size() - 2, hingeStiffness));
                    particleMass += springMass/2;
                }
            }
        }
        if(params_.particleFixed)
        {
            particleMass = std::numeric_limits<double>::infinity();
        }
        if (params_.connector != params_.CT_FLEXIBLE_ROD)
        {
            particles_.push_back(Particle(newParticlePos, particleMass, params_.particleFixed, false));
        }
        else
        {
            particles_[newParticleIndex].mass = particleMass;
        }
    }
    renderLock_.unlock();
}

void Simulation::addSaw(double x, double y)
{
    renderLock_.lock();
        saws_.push_back(Saw(Vector2d(x,y), params_.sawRadius));
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
    if(params_.activeForces & SimParameters::F_ELASTIC)
        processElasticBendingForce(q, F);
    if(params_.constraint == SimParameters::CH_PENALTY_FORCE)
        processPenaltyForce(q, F);

    H.setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());

}

void Simulation::processElasticBendingForce(const VectorXd &q, VectorXd &F)
{
    Vector3d pi,pj,pk,zUnit;
    zUnit.setZero();
    zUnit[2] = 1;
    int s1Id, s2Id, piIndex, pjIndex, pkIndex;

    for (int i; i<hinges_.size(); i++)
    {
        pi.setZero();
        pj.setZero();
        pk.setZero();
        s1Id = hinges_[i].s1;
        s2Id = hinges_[i].s2;

        if (springs_[s1Id].p1 == springs_[s2Id].p1)
        {
            pj.segment<2>(0) = q.segment<2>(springs_[s1Id].p1*2);
            pi.segment<2>(0) = q.segment<2>(springs_[s1Id].p2*2);
            pk.segment<2>(0) = q.segment<2>(springs_[s2Id].p2*2);
            pjIndex = springs_[s1Id].p1;
            piIndex = springs_[s1Id].p2;
            pkIndex = springs_[s2Id].p2;
        }
        else if (springs_[s1Id].p1 == springs_[s2Id].p2)
        {
            pj.segment<2>(0) = q.segment<2>(springs_[s1Id].p1*2);
            pi.segment<2>(0) = q.segment<2>(springs_[s1Id].p2*2);
            pk.segment<2>(0) = q.segment<2>(springs_[s2Id].p1*2);
            pjIndex = springs_[s1Id].p1;
            piIndex = springs_[s1Id].p2;
            pkIndex = springs_[s2Id].p1;
        }
        else if (springs_[s1Id].p2 == springs_[s2Id].p1)
        {
            pj.segment<2>(0) = q.segment<2>(springs_[s1Id].p2*2);
            pi.segment<2>(0) = q.segment<2>(springs_[s1Id].p1*2);
            pk.segment<2>(0) = q.segment<2>(springs_[s2Id].p2*2);
            pjIndex = springs_[s1Id].p2;
            piIndex = springs_[s1Id].p1;
            pkIndex = springs_[s2Id].p2;
        }
        else if (springs_[s1Id].p2 == springs_[s2Id].p2)
        {
            pj.segment<2>(0) = q.segment<2>(springs_[s1Id].p2*2);
            pi.segment<2>(0) = q.segment<2>(springs_[s1Id].p1*2);
            pk.segment<2>(0) = q.segment<2>(springs_[s2Id].p1*2);
            pjIndex = springs_[s1Id].p2;
            piIndex = springs_[s1Id].p1;
            pkIndex = springs_[s2Id].p1;
        }
        double y = ((pj - pi).cross(pk - pj)).dot(zUnit);
        double x = ((pj - pi).norm() * (pk - pj).norm()) + ((pj - pi).dot(pk - pj));
        double theta = 2 * atan2(y, x);
        Vector3d Fi = (hinges_[i].stiffness * theta * (pj - pi).cross(zUnit)) / (pj - pi).squaredNorm();
        Vector3d Fk = (hinges_[i].stiffness * theta * (pk - pj).cross(zUnit)) / (pk - pj).squaredNorm();
        F.segment<2>(piIndex * 2) += Fi.segment<2>(0);
        F.segment<2>(pkIndex * 2) += Fk.segment<2>(0);
        F.segment<2>(pjIndex * 2) += (-Fi-Fk).segment<2>(0);
    }
}

void Simulation::processPenaltyForce(const Eigen::VectorXd &q, Eigen::VectorXd &F)
{
    for(vector<Rod>::iterator it = rods_.begin(); it != rods_.end(); it++)
    {
        Vector2d particle1 = q.segment<2>(2*it->p1);
        Vector2d particle2 = q.segment<2>(2*it->p2);
        double dist = (particle2 - particle1).norm();
        int i1 = it->p1;
        int i2 = it->p2;

        double localFx = 4 * params_.penaltyStiffness  * (dist*dist - it->restlen*it->restlen) * (particle2[0] - particle1[0]);
        double localFy = 4 * params_.penaltyStiffness  * (dist*dist - it->restlen*it->restlen) * (particle2[1] - particle1[1]);
        F[i1*2] += localFx;
        F[i1*2+1] += localFy;
        F[i2*2] -= localFx;
        F[i2*2+1] -= localFy;
    }
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
        if (springs_[i].unsnappable)
        {
            continue;
        }
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
    VectorXd F;
    SparseMatrix<double> H;
    SparseMatrix<double> Minv;

    computeMassInverse(Minv);

    VectorXd oldq = q;
    if(params_.constraint == SimParameters::CH_PENALTY_FORCE)
    {
        q += params_.timeStep*v;
        computeForceAndHessian(q, oldq, F, H);
        v += params_.timeStep*Minv*F;
    }
    else if(params_.constraint == SimParameters::CH_STEP_PROJECT)
    {
        q += params_.timeStep*v;
        computeForceAndHessian(q, oldq, F, H);
        v += params_.timeStep*Minv*F;
        computeStepProject(q, oldq, v);
    }
    else if(params_.constraint == SimParameters::CH_LAGRANGE)
    {
        q += params_.timeStep*v;
        computeForceAndHessian(q, oldq, F, H);
        computeLagrangeMultipliers(q, F, v);
//        v += params_.timeStep*Minv*F;
    }
}

void Simulation::computeLagrangeMultipliers(const VectorXd &qVV, const VectorXd &F, VectorXd &v)
{
    VectorXd lamGuess(rods_.size());
    lamGuess.setZero();
    SparseMatrix<double> massInverseMatrix(qVV.rows(), qVV.rows());
    massInverseMatrix.setZero();
    computeMassInverse(massInverseMatrix);
    VectorXd fOfx(rods_.size());
    SparseMatrix<double> gradF(rods_.size(), rods_.size());

    VectorXd c = qVV + params_.timeStep*v + params_.timeStep * params_.timeStep * massInverseMatrix * F;
//    cout<<"Here 0"<<endl;
    SparseMatrix<double> gradGTransposeofqVV = computeGradGTranspose(qVV);
//    cout<<"Here 1"<<endl;
    cout<<"\n------------------------"<<endl;
    cout<<"\n------------------------"<<endl;
    cout<<"\n------------------------\n\n"<<endl;
    cout<<"\n\n Tranpose : \n"<<gradGTransposeofqVV;
    for(int newtonIt =0 ; newtonIt< params_.NewtonMaxIters; newtonIt++)
    {
        VectorXd qInside(qVV.rows());
        qInside.setZero();
        fOfx.setZero();
        gradF.setZero();
        cout<<"\n------------------------"<<endl;
        // Compute F of Lambda(i+1)
        qInside = c + (params_.timeStep * params_.timeStep) * massInverseMatrix * gradGTransposeofqVV * lamGuess;

        for(int i=0; i<rods_.size(); i++)
        {
            int p1pos = rods_[i].p1*2;
            int p2pos = rods_[i].p2*2;
            Vector2d p1 = qInside.segment<2>(p1pos);
            Vector2d p2 = qInside.segment<2>(p2pos);
            fOfx[i] += (p2-p1).squaredNorm() - rods_[i].restlen*rods_[i].restlen;
        }
        cout<<"\n F of X : \n\n"<<fOfx;
        if(fOfx.norm() < params_.NewtonTolerance)
        {
            cout<<"\n Not hitting newton iteration end";
            break;
        }
        //Compute Gradient of F(lambda i + 1)
        SparseMatrix<double> leftOfGradF(rods_.size(), qVV.rows());
        leftOfGradF.setZero();
        leftOfGradF = computeGradGTranspose(qInside).transpose();
        SparseMatrix<double> rightOfGradF(qVV.rows(), rods_.size());
        rightOfGradF.setZero();
        rightOfGradF = params_.timeStep * params_.timeStep * massInverseMatrix * gradGTransposeofqVV;

        gradF = leftOfGradF * rightOfGradF;
        cout<<"\n Grad F : \n\n"<<gradF<<endl;
        gradF.makeCompressed();
        SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
        solver.compute(gradF);
        VectorXd deltaLamda = solver.solve(-fOfx);
        lamGuess += deltaLamda;
    }
//    cout<<"Here 2"<<endl;
    cout<<"\n Lam Guess : \n\n"<<lamGuess<<endl;
    v = v + params_.timeStep * massInverseMatrix * F + params_.timeStep * massInverseMatrix * gradGTransposeofqVV * lamGuess;
    cout<<"\n Velocity : \n\n"<<v<<endl;

//        //Compute the Fofx
//        for(int i=0; i< rods_.size(); i++)
//        {
//            int p1pos = rods_[i].p1*2;
//            int p2pos = rods_[i].p2*2;
//            Vector2d p1 = q.segment<2>(p1pos);
//            Vector2d p2 = q.segment<2>(p2pos);

//            qInside.segment<2>(p1pos) += params_.timeStep*params_.timeStep*massInverseMatrix.coeff(p1pos, p1pos)*lamGuess[i]*(p2-p1)*2;
//            qInside.segment<2>(p2pos) += params_.timeStep*params_.timeStep*massInverseMatrix.coeff(p2pos, p2pos)*lamGuess[i]*(p1-p2)*2;
//        }

//        for(int i=0; i<rods_.size(); i++)
//        {
//            int p1pos = rods_[i].p1*2;
//            int p2pos = rods_[i].p2*2;
//            Vector2d p1 = qInside.segment<2>(p1pos);
//            Vector2d p2 = qInside.segment<2>(p2pos);
//            fOfx[i] = (p2-p1).squaredNorm() - rods_[i].restlen*rods_[i].restlen;
//        }
//        cout<<"\n F of x : "<<fOfx<<endl;
//        VectorXd rightSide(rods_.size()*2);
//        for (int i = 0; i < rods_.size(); i++)
//        {
//            int p1pos = rods_[i].p1*2;
//            int p2pos = rods_[i].p2*2;
//            Vector2d p1 = q.segment<2>(p1pos);
//            Vector2d p2 = q.segment<2>(p2pos);

//            Vector2d temp1 = params_.timeStep * params_.timeStep * massInverseMatrix(p1pos, p1pos) * (p2-p1) * 2;
//            Vector2d temp2 = params_.timeStep * params_.timeStep * massInverseMatrix(p2pos, p2pos) * (p1-p2) * 2;
//            // TODO : Complete this


//        }
}

Eigen::SparseMatrix<double> Simulation::computeGradGTranspose(const VectorXd &q)
{
    SparseMatrix<double> gradGTranspose(q.rows(), rods_.size());
    gradGTranspose.setZero();

    for(int i=0; i<rods_.size(); i++)
    {
        int p1index = rods_[i].p1*2;
        int p2index = rods_[i].p2*2;
        Vector2d p1 = q.segment<2>(p1index);
        Vector2d p2 = q.segment<2>(p2index);

        Vector2d gradient1 = (p2 - p1)*-2;
        Vector2d gradient2 = (p1 - p2)*2;

        gradGTranspose.coeffRef(p1index, i) += gradient1[0];
        gradGTranspose.coeffRef(p1index+1, i) += gradient1[1];

        gradGTranspose.coeffRef(p2index, i) += gradient2[0];
        gradGTranspose.coeffRef(p2index+1, i) += gradient2[1];
    }
    return gradGTranspose;

}

void Simulation::computeStepProject(VectorXd &q, VectorXd &oldq, VectorXd &v)
{
    VectorXd qGuess = q;
    VectorXd lamGuess(rods_.size());
    lamGuess.setZero();
    SparseMatrix<double> massInverseMatrix(q.rows(), q.rows());
    massInverseMatrix.setZero();
    computeMassInverse(massInverseMatrix);
    VectorXd forceX(qGuess.rows() + rods_.size());

    for(int i = 0; i < params_.NewtonMaxIters; i++)
    {
        forceX.setZero();
        for(int j = 0; j < rods_.size(); j++)
        {
            int p1pos = rods_[j].p1*2;
            int p2pos = rods_[j].p2*2;
            Vector2d p1 = qGuess.segment<2>(2*rods_[j].p1);
            Vector2d p2 = qGuess.segment<2>(2*rods_[j].p2);
            forceX.segment<2>(rods_[j].p1*2) += lamGuess[j] * massInverseMatrix.coeff(p1pos, p1pos) * (p2 - p1) * 2;
            forceX.segment<2>(rods_[j].p2*2) += lamGuess[j] * massInverseMatrix.coeff(p2pos, p2pos) * (p1 - p2) * 2;
        }

        VectorXd qDifference(qGuess.rows());
        qDifference.setZero();
        qDifference = q - qGuess;
        for(int j =0; j < qDifference.rows(); j++)
        {
            forceX[j] += qDifference[j];
        }

        for(int j = 0; j < rods_.size(); j++)
        {
            Vector2d dstpos(qGuess[rods_[j].p1*2], qGuess[rods_[j].p1*2+1]);
            Vector2d srcpos(qGuess[rods_[j].p2*2], qGuess[rods_[j].p2*2+1]);
            double dist = (dstpos-srcpos).norm();
            forceX[j+qGuess.rows()] += dist*dist - rods_[j].restlen*rods_[j].restlen;
        }
        if(forceX.norm() < params_.NewtonTolerance)
        {
            break;
        }

        // GRADIENT CALCULATION

        SparseMatrix<double> forceGradient(qGuess.rows()+rods_.size(), qGuess.rows()+rods_.size());
        forceGradient.setZero();

        // Calculating Top Left of Force Gradient Matrix
        vector< Triplet<double> > topLeftTriplets;
        topLeftTriplets.reserve(rods_.size()*8);
        double dfxidxi, dfxidxj, dfyidyi, dfyidyj, dfxjdxi, dfxjdxj, dfyjdyi, dfyjdyj;
        for(int j=0; j<rods_.size(); j++)
        {
            int p1pos = rods_[j].p1*2;
            int p2pos = rods_[j].p2*2;
            dfxidxi = -2*lamGuess[j]*massInverseMatrix.coeff(p1pos, p1pos);
            dfxidxj = -dfxidxi;
            dfyidyi = dfxidxi;
            dfyidyj = -dfyidyi;

            dfxjdxi = 2*lamGuess[j]*massInverseMatrix.coeff(p2pos, p2pos);
            dfxjdxj = -dfxjdxi;
            dfyjdyi = dfxjdxi;
            dfyjdyj = -dfyjdyi;
            topLeftTriplets.push_back(Triplet<double>(p1pos, p1pos, dfxidxi));
            topLeftTriplets.push_back(Triplet<double>(p1pos, p2pos, dfxidxj));
            topLeftTriplets.push_back(Triplet<double>(p1pos + 1, p1pos + 1, dfyidyi));
            topLeftTriplets.push_back(Triplet<double>(p1pos + 1, p2pos + 1, dfyidyj));

            topLeftTriplets.push_back(Triplet<double>(p2pos, p1pos, dfxjdxi));
            topLeftTriplets.push_back(Triplet<double>(p2pos, p2pos, dfxjdxj));
            topLeftTriplets.push_back(Triplet<double>(p2pos + 1, p1pos + 1, dfyjdyi));
            topLeftTriplets.push_back(Triplet<double>(p2pos + 1, p2pos + 1, dfyjdyj));
        }
        SparseMatrix<double> topLeftMatrix(qGuess.rows(), qGuess.rows());
        topLeftMatrix.setZero();
        topLeftMatrix.setFromTriplets(topLeftTriplets.begin(), topLeftTriplets.end());
        SparseMatrix<double> identity(qGuess.rows(), qGuess.rows());
        identity.setIdentity();
        topLeftMatrix += identity;
        for(int j = 0; j < topLeftMatrix.rows(); j++)
        {
            for(int k = 0; k < topLeftMatrix.cols(); k++)
            {
                forceGradient.coeffRef(j,k) = topLeftMatrix.coeff(j,k);
            }
        }
        // Calculating Top Right and Bottom Left of Force Gradient Matrix
        for(int j = 0; j < rods_.size(); j++)
        {
            int n = j + qGuess.rows();
            int p1pos = rods_[j].p1 *2;
            int p2pos = rods_[j].p2 *2;
            Vector2d p1Grad = (qGuess.segment<2>(p2pos) - qGuess.segment<2>(p1pos)) * 2;
            Vector2d p2Grad = (qGuess.segment<2>(p1pos) - qGuess.segment<2>(p2pos)) * 2;
            forceGradient.coeffRef(n,p1pos) += p1Grad[0];
            forceGradient.coeffRef(n,p1pos + 1) += p1Grad[1];
            forceGradient.coeffRef(n,p2pos) += p2Grad[0];
            forceGradient.coeffRef(n,p2pos + 1) += p2Grad[1];
            forceGradient.coeffRef(p1pos,n) += p1Grad[0] * massInverseMatrix.coeff(p1pos,p1pos);
            forceGradient.coeffRef(p1pos + 1,n) += p1Grad[1] * massInverseMatrix.coeff(p1pos,p1pos);
            forceGradient.coeffRef(p2pos,n) += p2Grad[0] * massInverseMatrix.coeff(p2pos,p2pos);
            forceGradient.coeffRef(p2pos + 1,n) += p2Grad[1] * massInverseMatrix.coeff(p2pos,p2pos);
        }
        forceGradient.makeCompressed();
        SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;
        solver.compute(forceGradient);
        VectorXd deltaguess = solver.solve(-forceX);
        // cout << "deltaGuess " << deltaguess << endl;
        for(int j = 0; j < qGuess.rows(); j++)
        {
            qGuess[j] -= deltaguess[j];
        }
        int n = qGuess.rows();
        for(int j = 0; j < lamGuess.rows(); j++)
        {
            lamGuess[j] += deltaguess[j+n];
        }
    }
    v = v + (qGuess - q)/params_.timeStep;
    q = qGuess;
}

void Simulation::pruneOverstrainedSprings()
{
    set<int> springstodelete;

    vector<Spring> newsprings;
    vector<int> remainingspringmap;
    int nsprings = springs_.size();
    for(int i=0; i<nsprings; i++)
    {
        Vector2d p1 = particles_[springs_[i].p1].pos;
        Vector2d p2 = particles_[springs_[i].p2].pos;
        double d = (p2-p1).norm();

        double strain = (d - springs_[i].restlen)/springs_[i].restlen;
        if(!springs_[i].unsnappable && strain > params_.maxSpringStrain)
        {
            springstodelete.insert(i);
        }
    }

    renderLock_.lock();
    {
        if(!springstodelete.empty())
        {
            for(int i=0; i<(int)springs_.size(); i++)
            {
                if(springstodelete.count(i) == 0)
                {
                    remainingspringmap.push_back(newsprings.size());
                    newsprings.push_back(springs_[i]);
                }
                else
                {
                    remainingspringmap.push_back(-1);
                }
            }
        }
        if(!springstodelete.empty())
        {
            springs_ = newsprings;
            for(vector<Hinge>::iterator hinge = hinges_.begin(); hinge != hinges_.end(); ++hinge)
            {
                hinge->s1 = remainingspringmap[hinge->s1];
                hinge->s2 = remainingspringmap[hinge->s2];
            }
        }
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

void Simulation::detectSawedSprings(std::set<int> &springsToDelete, std::set<int> &hingesToDelete)
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
                particles_[springs_[i].p1].mass -= springs_[i].mass/2;
                particles_[springs_[i].p2].mass -= springs_[i].mass/2;
                springsToDelete.insert(i);
                for(int j=0; j<(int)hinges_.size(); j++)
                {
                    if (hinges_[j].s1 == i || hinges_[j].s2 == i)
                    {
                        hingesToDelete.insert(j);
                    }
                }
                break;
            }
        }
    }
}

void Simulation::detectSawedRods(std::set<int> &rodsToDelete)
{
    for(int i=0; i<(int)rods_.size(); i++)
    {
        Vector2d pos1 = particles_[rods_[i].p1].pos;
        Vector2d pos2 = particles_[rods_[i].p2].pos;
        double maxx = max(pos1[0], pos2[0]);
        double minx = min(pos1[0], pos2[0]);
        double maxy = max(pos1[1], pos2[1]);
        double miny = min(pos1[1], pos2[1]);
        for(vector<Saw>::iterator saw = saws_.begin(); saw != saws_.end(); ++saw)
        {
            Vector2d sawpos = saw->pos;
            double sawRadius = saw->radius;

            if(sawpos[0] - sawRadius > maxx || sawpos[0] + sawRadius < minx || sawpos[1] - sawRadius > maxy || sawpos[1] + sawRadius < miny)
                continue;

            double sawRodDistance = ptSegmentDist(sawpos, pos1, pos2);
            if(sawRodDistance <= sawRadius)
            {
                rodsToDelete.insert(i);
                break;
            }
        }
    }
}

void Simulation::detectSawedParticles(std::set<int> &particlesToDelete)
{
    for(int i=0; i<(int)particles_.size(); i++)
    {
        Vector2d particlePos = particles_[i].pos;

        if(fabs(particlePos[0]) > 2 || fabs(particlePos[1]) > 2)
        {
            particlesToDelete.insert(i);
            break;
        }

        for(vector<Saw>::iterator it = saws_.begin(); it != saws_.end(); ++it)
        {
            Vector2d sawpos = it->pos;
            double sqdist = (sawpos-particlePos).squaredNorm();
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
    set<int> hingestodelete;
    set<int> rodsToDelete;
    detectSawedParticles(particlestodelete);
    detectSawedSprings(springstodelete, hingestodelete);
    detectSawedRods(rodsToDelete);

    vector<Particle> newparticles;
    vector<Spring> newsprings;
    vector<Hinge> newhinges;
    vector<Rod> newrods;
    vector<int> remainingparticlemap;
    vector<int> remainingspringmap;
    if(!particlestodelete.empty())
    {
        for(int i=0; i<(int)springs_.size(); i++)
        {
            if(particlestodelete.count(springs_[i].p1) || particlestodelete.count(springs_[i].p2))
            {
                springstodelete.insert(i);
                for(int j=0; j<(int)hinges_.size(); j++)
                {
                    if (springstodelete.count(hinges_[j].s1) || springstodelete.count(hinges_[j].s2))
                    {
                        hingestodelete.insert(j);
                    }
                }
            }

        }
        for(int i=0; i<(int)rods_.size(); i++)
        {
            if(particlestodelete.count(rods_[i].p1) || particlestodelete.count(rods_[i].p2))
                rodsToDelete.insert(i);
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
                remainingspringmap.push_back(newsprings.size());
                newsprings.push_back(springs_[i]);
            }
            else
            {
                remainingspringmap.push_back(-1);
            }
        }
    }
    if(!rodsToDelete.empty())
    {
        for(int i=0; i<(int)rods_.size(); i++)
        {
            if(rodsToDelete.count(i) == 0)
            {
                newrods.push_back(rods_[i]);
            }
        }
    }
    if(!hingestodelete.empty())
    {
        for(int i=0; i<(int)hinges_.size(); i++)
        {
            if(hingestodelete.count(i) == 0)
            {
                newhinges.push_back(hinges_[i]);
            }
        }
    }
    if(!springstodelete.empty() || !particlestodelete.empty() || !rodsToDelete.empty() || !hingestodelete.empty())
    {
        renderLock_.lock();
        {
            if(!hingestodelete.empty())
                hinges_ = newhinges;
            if(!rodsToDelete.empty())
                rods_ = newrods;
            if(!springstodelete.empty())
            {
                springs_ = newsprings;
                for(vector<Hinge>::iterator hinge = hinges_.begin(); hinge != hinges_.end(); ++hinge)
                {
                    hinge->s1 = remainingspringmap[hinge->s1];
                    hinge->s2 = remainingspringmap[hinge->s2];
                }
            }
            if(!particlestodelete.empty())
            {
                particles_ = newparticles;
                for(vector<Spring>::iterator it = springs_.begin(); it != springs_.end(); ++it)
                {
                    it->p1 = remainingparticlemap[it->p1];
                    it->p2 = remainingparticlemap[it->p2];
                }
                for(vector<Rod>::iterator it = rods_.begin(); it != rods_.end(); ++it)
                {
                    it->p1 = remainingparticlemap[it->p1];
                    it->p2 = remainingparticlemap[it->p2];
                }
            }
        }
        renderLock_.unlock();
    }
}

void Simulation::clearScene()
{
    renderLock_.lock();
    {
        particles_.clear();
        springs_.clear();
        hinges_.clear();
        saws_.clear();
        rods_.clear();
    }
    renderLock_.unlock();
}

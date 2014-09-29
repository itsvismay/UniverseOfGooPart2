#include "controller.h"
#include "mainwindow.h"
#include "simulation.h"
#include <QDebug>

Controller::Controller(int fps) : QThread(), mw_(NULL), fps_(fps)
{
}

Controller::~Controller()
{
    delete sim_;
}

void Controller::initialize(MainWindow *mw)
{
    mw_ = mw;
    sim_ = new Simulation(params_);
}

void Controller::run()
{
    reset();
    connect(&simtimer_, SIGNAL(timeout()), this, SLOT(simTick()));
    simtimer_.start(1000/fps_);
    exec();
}

void Controller::reset()
{
    params_ = SimParameters();
    QMetaObject::invokeMethod(mw_, "setUIFromParameters", Q_ARG(SimParameters, params_));
    clearScene();
}

void Controller::clearScene()
{
    sim_->clearScene();
}

void Controller::setupGameMode()
{
    sim_->setupGameMode();
}

void Controller::updateParameters(SimParameters params)
{
    params_ = params;
}

void Controller::render()
{
    sim_->render();
}

void Controller::mouseClicked(double x, double y)
{
    switch(params_.clickMode)
    {
    case SimParameters::CM_ADDPARTICLE:
        sim_->addParticle(x, y);
        break;
    case SimParameters::CM_ADDSAW:
        sim_->addSaw(x,y);
        break;
    }
}

void Controller::simTick()
{
    if(params_.simRunning)
    {
        sim_->takeSimulationStep();
    }
}

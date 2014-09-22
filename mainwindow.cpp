#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "simparameters.h"
#include "controller.h"

MainWindow::MainWindow(Controller &cont, int fps, QWidget *parent) :
    QMainWindow(parent),
    cont_(cont),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->GLWidget->setController(&cont);
    simRunning_ = false;
    connect(&renderTimer_, SIGNAL(timeout()), this, SLOT(updateGL()));
    renderTimer_.start(1000/fps);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionExit_triggered()
{
    close();
}

void MainWindow::setParametersFromUI()
{
    SimParameters params;

    params.simRunning = simRunning_;

    if(ui->penaltyForceButton->isChecked())
        params.constraint = SimParameters::CH_PENALTY_FORCE;
    else if(ui->stepAndProjectButton->isChecked())
        params.constraint = SimParameters::CH_STEP_PROJECT;
    else if(ui->lagranceMultiplierButton->isChecked())
        params.constraint = SimParameters::CH_LAGRANGE;


    params.timeStep = ui->timeStepEdit->text().toDouble();
    params.NewtonTolerance = ui->newtonTolEdit->text().toDouble();
    params.NewtonMaxIters = ui->newtonMaxItersEdit->text().toInt();
    params.penaltyStiffness = ui->penaltyStiffnessEdit->text().toDouble();

    params.activeForces = 0;
    if(ui->gravityCheckBox->isChecked())
        params.activeForces |= SimParameters::F_GRAVITY;
    if(ui->springsCheckBox->isChecked())
        params.activeForces |= SimParameters::F_SPRINGS;
    if(ui->floorCheckBox->isChecked())
        params.activeForces |= SimParameters::F_FLOOR;
    if(ui->dampingStiffnessCheckBox->isChecked())
        params.activeForces |= SimParameters::F_DAMPING;
    if(ui->elasticBendingCheckbox->isChecked())
        params.activeForces |= SimParameters::F_ELASTIC;

    params.gravityG = ui->gravityGEdit->text().toDouble();
    params.springStiffness = ui->springStiffnessEdit->text().toDouble();
    params.maxSpringStrain = ui->maxStrainEdit->text().toDouble();
    params.dampingStiffness = ui->dampingStiffnessEdit->text().toDouble();

    if(ui->addParticleButton->isChecked())
        params.clickMode = SimParameters::CM_ADDPARTICLE;
    else if(ui->addSawButton->isChecked())
        params.clickMode = SimParameters::CM_ADDSAW;

    params.particleMass = ui->massEdit->text().toDouble();
    params.maxSpringDist = ui->maxSpringDistEdit->text().toDouble();
    params.particleFixed = ui->isFixedCheckBox->isChecked();

    params.sawRadius = ui->radiusEdit->text().toDouble();

    if(ui->springButton->isChecked())
        params.connector = SimParameters::CT_SPRING;
    else if(ui->rigidRodButton->isChecked())
        params.connector = SimParameters:: CT_RIGID_ROD;
    else if(ui->flexibleRodButton->isChecked())
        params.connector = SimParameters:: CT_FLEXIBLE_ROD;
    else if(ui->ropeButton->isChecked())
        params.connector = SimParameters:: CT_ROPE;

    params.rodDensity = ui->rodDensityEdit->text().toDouble();
    params.rodStretch = ui->rodStretchEdit->text().toDouble();
    params.rodBend = ui->rodBendEdit->text().toDouble();
    params.rodSegments = ui->rodSegmentsEdit->text().toInt();
    params.ropeDensity = ui->ropeDensityEdit->text().toDouble();
    params.ropeBend = ui->ropeBendEdit->text().toDouble();
    params.ropeSegments = ui->rodSegmentsEdit->text().toInt();



    setUIFromParameters(params);
    QMetaObject::invokeMethod(&cont_, "updateParameters", Q_ARG(SimParameters, params));
}

void MainWindow::setUIFromParameters(const SimParameters &params)
{
    if(params.simRunning)
    {
        ui->startSimulationButton->setText(QString("Pause Simulation"));
        simRunning_ = true;
    }
    else
    {
        ui->startSimulationButton->setText(QString("Start Simulation"));
        simRunning_ = false;
    }

    switch(params.constraint)
    {
        case SimParameters::CH_PENALTY_FORCE:
            ui->penaltyForceButton->setChecked(true);
            break;
        case SimParameters::CH_STEP_PROJECT:
            ui->stepAndProjectButton->setChecked(true);
            break;
        case SimParameters::CH_LAGRANGE:
            ui->lagranceMultiplierButton->setChecked(true);
            break;
    }

    ui->timeStepEdit->setText(QString::number(params.timeStep));
    ui->newtonTolEdit->setText(QString::number(params.NewtonTolerance));
    ui->newtonMaxItersEdit->setText(QString::number(params.NewtonMaxIters));
    ui->penaltyStiffnessEdit->setText(QString::number(params.penaltyStiffness));

    ui->gravityCheckBox->setChecked(params.activeForces & SimParameters::F_GRAVITY);
    ui->springsCheckBox->setChecked(params.activeForces & SimParameters::F_SPRINGS);
    ui->floorCheckBox->setChecked(params.activeForces & SimParameters::F_FLOOR);
    ui->dampingStiffnessCheckBox->setChecked(params.activeForces & SimParameters::F_DAMPING);
    ui->elasticBendingCheckbox->setChecked(params.activeForces & SimParameters::F_ELASTIC);

    ui->gravityGEdit->setText(QString::number(params.gravityG));
    ui->springStiffnessEdit->setText(QString::number(params.springStiffness));
    ui->maxStrainEdit->setText(QString::number(params.maxSpringStrain));
    ui->dampingStiffnessEdit->setText(QString::number(params.dampingStiffness));

    if(params.clickMode == SimParameters::CM_ADDPARTICLE)
        ui->addParticleButton->setChecked(true);
    else if(params.clickMode == SimParameters::CM_ADDSAW)
        ui->addSawButton->setChecked(true);

    ui->massEdit->setText(QString::number(params.particleMass));
    ui->maxSpringDistEdit->setText(QString::number(params.maxSpringDist));
    ui->isFixedCheckBox->setChecked(params.particleFixed);
    ui->radiusEdit->setText(QString::number(params.sawRadius));

    if(params.connector == SimParameters:: CT_SPRING)
        ui->springButton->setChecked(true);
    else if(params.connector == SimParameters:: CT_RIGID_ROD)
        ui->rigidRodButton->setChecked(true);
    else if(params.connector == SimParameters:: CT_FLEXIBLE_ROD)
        ui->flexibleRodButton->setChecked(true);
    else if(params.connector == SimParameters:: CT_ROPE)
        ui->ropeButton->setChecked(true);

    ui->rodDensityEdit->setText(QString::number(params.rodDensity));
    ui->rodStretchEdit->setText(QString::number(params.rodStretch));
    ui->rodBendEdit->setText(QString::number(params.rodBend));
    ui->rodSegmentsEdit->setText(QString::number(params.rodSegments));
    ui->ropeDensityEdit->setText(QString::number(params.ropeDensity));
    ui->ropeBendEdit->setText(QString::number(params.ropeBend));
    ui->ropeSegmentsEdit->setText(QString::number(params.ropeSegments));

}

void MainWindow::updateGL()
{
    ui->GLWidget->update();
}

void MainWindow::on_actionReset_Everything_triggered()
{
    QMetaObject::invokeMethod(&cont_, "reset");
}

void MainWindow::on_actionReset_triggered()
{
    QMetaObject::invokeMethod(&cont_, "clearScene");
}

void MainWindow::on_startSimulationButton_clicked()
{
    simRunning_ = !simRunning_;
    setParametersFromUI();
}


void MainWindow::on_timeStepEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_newtonTolEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_newtonMaxItersEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_gravityCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_springsCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_floorCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_dampingStiffnessCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_gravityGEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_springStiffnessEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_maxStrainEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_dampingStiffnessEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_addParticleButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_addSawButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_massEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_maxSpringDistEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_isFixedCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_radiusEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_penaltyForceButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_stepAndProjectButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_lagrangeMultiplierButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_penaltyStiffnessEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_elasticBendingCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_springButton_clicked()
{
    setParametersFromUI();
}
void MainWindow::on_rigidRodButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_flexibleRodButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_ropeButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_rodDensityEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_rodStretchEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_rodBendEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_rodSegmentsEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_ropeDensityEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_ropeBendEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_ropeSegmentsEdit_editingFinished()
{
    setParametersFromUI();
}


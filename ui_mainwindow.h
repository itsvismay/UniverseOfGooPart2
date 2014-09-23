/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Tue Sep 23 15:52:00 2014
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QFrame>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "glpanel.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionExit;
    QAction *actionReset;
    QAction *actionReset_Everything;
    QWidget *centralWidget;
    GLPanel *GLWidget;
    QFrame *parameterFrame;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout;
    QGroupBox *simOptionsBox;
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout;
    QGroupBox *ConstraintBox;
    QRadioButton *penaltyForceButton;
    QRadioButton *stepAndProjectButton;
    QRadioButton *lagrangeMultiplierButton;
    QPushButton *startSimulationButton;
    QGroupBox *SimParametersBox;
    QLabel *timeStepLabel;
    QLabel *newtonTolLabel;
    QLabel *newtonMaxItersLabel;
    QLineEdit *timeStepEdit;
    QLineEdit *newtonTolEdit;
    QLineEdit *newtonMaxItersEdit;
    QLineEdit *penaltyStiffnessEdit;
    QLabel *penaltyStiffness;
    QGroupBox *activeForcesBox;
    QCheckBox *gravityCheckBox;
    QCheckBox *springsCheckBox;
    QCheckBox *floorCheckBox;
    QLabel *gravityGLabel;
    QLineEdit *gravityGEdit;
    QLabel *springStiffnessLabel;
    QLineEdit *springStiffnessEdit;
    QLabel *maxStrainLabel;
    QLineEdit *maxStrainEdit;
    QCheckBox *dampingStiffnessCheckBox;
    QLabel *dampingStiffnessLabel;
    QLineEdit *dampingStiffnessEdit;
    QCheckBox *elasticBendingCheckBox;
    QGroupBox *uiOptionsBox;
    QWidget *horizontalLayoutWidget_2;
    QHBoxLayout *horizontalLayout_2;
    QGroupBox *clickFunctionBox;
    QRadioButton *addParticleButton;
    QRadioButton *addSawButton;
    QCheckBox *isFixedCheckBox;
    QGroupBox *ConnectorTypeBox;
    QRadioButton *springButton;
    QRadioButton *flexibleRodButton;
    QRadioButton *ropeButton;
    QRadioButton *rigidRodButton;
    QGroupBox *particleSettingsBox;
    QLabel *massLabel;
    QLabel *maxSpringDistLabel;
    QLineEdit *massEdit;
    QLineEdit *maxSpringDistEdit;
    QGroupBox *rodSettingsBox;
    QLabel *rodDensityLabel;
    QLabel *rodStretchLabel;
    QLineEdit *rodDensityEdit;
    QLineEdit *rodStretchEdit;
    QLabel *rodBendLabel;
    QLineEdit *rodBendEdit;
    QLineEdit *rodSegmentsEdit;
    QLabel *rodSegmentsLabel;
    QGroupBox *sawSettingsBox;
    QLabel *radiusLabel;
    QLineEdit *radiusEdit;
    QGroupBox *ropeSettingsBox;
    QLabel *ropeDensityLabel;
    QLineEdit *ropeDensityEdit;
    QLabel *ropeBendkLabel;
    QLineEdit *ropeBendEdit;
    QLabel *ropeSegmentsLabel;
    QLineEdit *ropeSegmentsEdit;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuScene;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1200, 800);
        actionExit = new QAction(MainWindow);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        actionReset = new QAction(MainWindow);
        actionReset->setObjectName(QString::fromUtf8("actionReset"));
        actionReset_Everything = new QAction(MainWindow);
        actionReset_Everything->setObjectName(QString::fromUtf8("actionReset_Everything"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        GLWidget = new GLPanel(centralWidget);
        GLWidget->setObjectName(QString::fromUtf8("GLWidget"));
        GLWidget->setGeometry(QRect(10, 0, 731, 731));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(GLWidget->sizePolicy().hasHeightForWidth());
        GLWidget->setSizePolicy(sizePolicy);
        parameterFrame = new QFrame(centralWidget);
        parameterFrame->setObjectName(QString::fromUtf8("parameterFrame"));
        parameterFrame->setGeometry(QRect(749, -1, 441, 731));
        parameterFrame->setFrameShape(QFrame::StyledPanel);
        parameterFrame->setFrameShadow(QFrame::Raised);
        verticalLayoutWidget = new QWidget(parameterFrame);
        verticalLayoutWidget->setObjectName(QString::fromUtf8("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(9, -1, 431, 731));
        verticalLayout = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        simOptionsBox = new QGroupBox(verticalLayoutWidget);
        simOptionsBox->setObjectName(QString::fromUtf8("simOptionsBox"));
        horizontalLayoutWidget = new QWidget(simOptionsBox);
        horizontalLayoutWidget->setObjectName(QString::fromUtf8("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(9, 19, 421, 211));
        horizontalLayout = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        ConstraintBox = new QGroupBox(horizontalLayoutWidget);
        ConstraintBox->setObjectName(QString::fromUtf8("ConstraintBox"));
        penaltyForceButton = new QRadioButton(ConstraintBox);
        penaltyForceButton->setObjectName(QString::fromUtf8("penaltyForceButton"));
        penaltyForceButton->setGeometry(QRect(20, 30, 121, 22));
        stepAndProjectButton = new QRadioButton(ConstraintBox);
        stepAndProjectButton->setObjectName(QString::fromUtf8("stepAndProjectButton"));
        stepAndProjectButton->setGeometry(QRect(20, 50, 181, 22));
        lagrangeMultiplierButton = new QRadioButton(ConstraintBox);
        lagrangeMultiplierButton->setObjectName(QString::fromUtf8("lagrangeMultiplierButton"));
        lagrangeMultiplierButton->setGeometry(QRect(20, 70, 181, 22));
        startSimulationButton = new QPushButton(ConstraintBox);
        startSimulationButton->setObjectName(QString::fromUtf8("startSimulationButton"));
        startSimulationButton->setGeometry(QRect(10, 130, 181, 27));

        horizontalLayout->addWidget(ConstraintBox);

        SimParametersBox = new QGroupBox(horizontalLayoutWidget);
        SimParametersBox->setObjectName(QString::fromUtf8("SimParametersBox"));
        timeStepLabel = new QLabel(SimParametersBox);
        timeStepLabel->setObjectName(QString::fromUtf8("timeStepLabel"));
        timeStepLabel->setGeometry(QRect(10, 30, 81, 21));
        newtonTolLabel = new QLabel(SimParametersBox);
        newtonTolLabel->setObjectName(QString::fromUtf8("newtonTolLabel"));
        newtonTolLabel->setGeometry(QRect(10, 50, 131, 21));
        newtonMaxItersLabel = new QLabel(SimParametersBox);
        newtonMaxItersLabel->setObjectName(QString::fromUtf8("newtonMaxItersLabel"));
        newtonMaxItersLabel->setGeometry(QRect(10, 70, 131, 21));
        timeStepEdit = new QLineEdit(SimParametersBox);
        timeStepEdit->setObjectName(QString::fromUtf8("timeStepEdit"));
        timeStepEdit->setGeometry(QRect(150, 30, 51, 21));
        newtonTolEdit = new QLineEdit(SimParametersBox);
        newtonTolEdit->setObjectName(QString::fromUtf8("newtonTolEdit"));
        newtonTolEdit->setGeometry(QRect(150, 50, 51, 21));
        newtonMaxItersEdit = new QLineEdit(SimParametersBox);
        newtonMaxItersEdit->setObjectName(QString::fromUtf8("newtonMaxItersEdit"));
        newtonMaxItersEdit->setGeometry(QRect(150, 70, 51, 21));
        penaltyStiffnessEdit = new QLineEdit(SimParametersBox);
        penaltyStiffnessEdit->setObjectName(QString::fromUtf8("penaltyStiffnessEdit"));
        penaltyStiffnessEdit->setGeometry(QRect(150, 90, 51, 21));
        penaltyStiffness = new QLabel(SimParametersBox);
        penaltyStiffness->setObjectName(QString::fromUtf8("penaltyStiffness"));
        penaltyStiffness->setGeometry(QRect(10, 90, 131, 21));

        horizontalLayout->addWidget(SimParametersBox);


        verticalLayout->addWidget(simOptionsBox);

        activeForcesBox = new QGroupBox(verticalLayoutWidget);
        activeForcesBox->setObjectName(QString::fromUtf8("activeForcesBox"));
        gravityCheckBox = new QCheckBox(activeForcesBox);
        gravityCheckBox->setObjectName(QString::fromUtf8("gravityCheckBox"));
        gravityCheckBox->setGeometry(QRect(30, 30, 97, 21));
        springsCheckBox = new QCheckBox(activeForcesBox);
        springsCheckBox->setObjectName(QString::fromUtf8("springsCheckBox"));
        springsCheckBox->setGeometry(QRect(30, 50, 97, 21));
        floorCheckBox = new QCheckBox(activeForcesBox);
        floorCheckBox->setObjectName(QString::fromUtf8("floorCheckBox"));
        floorCheckBox->setGeometry(QRect(30, 90, 97, 21));
        gravityGLabel = new QLabel(activeForcesBox);
        gravityGLabel->setObjectName(QString::fromUtf8("gravityGLabel"));
        gravityGLabel->setGeometry(QRect(230, 30, 121, 21));
        gravityGEdit = new QLineEdit(activeForcesBox);
        gravityGEdit->setObjectName(QString::fromUtf8("gravityGEdit"));
        gravityGEdit->setGeometry(QRect(370, 30, 51, 21));
        springStiffnessLabel = new QLabel(activeForcesBox);
        springStiffnessLabel->setObjectName(QString::fromUtf8("springStiffnessLabel"));
        springStiffnessLabel->setGeometry(QRect(230, 50, 121, 21));
        springStiffnessEdit = new QLineEdit(activeForcesBox);
        springStiffnessEdit->setObjectName(QString::fromUtf8("springStiffnessEdit"));
        springStiffnessEdit->setGeometry(QRect(370, 50, 51, 21));
        maxStrainLabel = new QLabel(activeForcesBox);
        maxStrainLabel->setObjectName(QString::fromUtf8("maxStrainLabel"));
        maxStrainLabel->setGeometry(QRect(230, 70, 121, 21));
        maxStrainEdit = new QLineEdit(activeForcesBox);
        maxStrainEdit->setObjectName(QString::fromUtf8("maxStrainEdit"));
        maxStrainEdit->setGeometry(QRect(370, 70, 51, 21));
        dampingStiffnessCheckBox = new QCheckBox(activeForcesBox);
        dampingStiffnessCheckBox->setObjectName(QString::fromUtf8("dampingStiffnessCheckBox"));
        dampingStiffnessCheckBox->setGeometry(QRect(30, 110, 151, 21));
        dampingStiffnessLabel = new QLabel(activeForcesBox);
        dampingStiffnessLabel->setObjectName(QString::fromUtf8("dampingStiffnessLabel"));
        dampingStiffnessLabel->setGeometry(QRect(230, 110, 121, 21));
        dampingStiffnessEdit = new QLineEdit(activeForcesBox);
        dampingStiffnessEdit->setObjectName(QString::fromUtf8("dampingStiffnessEdit"));
        dampingStiffnessEdit->setGeometry(QRect(370, 110, 51, 21));
        elasticBendingCheckBox = new QCheckBox(activeForcesBox);
        elasticBendingCheckBox->setObjectName(QString::fromUtf8("elasticBendingCheckBox"));
        elasticBendingCheckBox->setGeometry(QRect(30, 130, 151, 21));

        verticalLayout->addWidget(activeForcesBox);

        uiOptionsBox = new QGroupBox(verticalLayoutWidget);
        uiOptionsBox->setObjectName(QString::fromUtf8("uiOptionsBox"));
        horizontalLayoutWidget_2 = new QWidget(uiOptionsBox);
        horizontalLayoutWidget_2->setObjectName(QString::fromUtf8("horizontalLayoutWidget_2"));
        horizontalLayoutWidget_2->setGeometry(QRect(10, 20, 421, 211));
        horizontalLayout_2 = new QHBoxLayout(horizontalLayoutWidget_2);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        clickFunctionBox = new QGroupBox(horizontalLayoutWidget_2);
        clickFunctionBox->setObjectName(QString::fromUtf8("clickFunctionBox"));
        addParticleButton = new QRadioButton(clickFunctionBox);
        addParticleButton->setObjectName(QString::fromUtf8("addParticleButton"));
        addParticleButton->setGeometry(QRect(10, 30, 117, 21));
        addSawButton = new QRadioButton(clickFunctionBox);
        addSawButton->setObjectName(QString::fromUtf8("addSawButton"));
        addSawButton->setGeometry(QRect(10, 50, 117, 21));
        isFixedCheckBox = new QCheckBox(clickFunctionBox);
        isFixedCheckBox->setObjectName(QString::fromUtf8("isFixedCheckBox"));
        isFixedCheckBox->setGeometry(QRect(10, 70, 97, 21));
        ConnectorTypeBox = new QGroupBox(clickFunctionBox);
        ConnectorTypeBox->setObjectName(QString::fromUtf8("ConnectorTypeBox"));
        ConnectorTypeBox->setGeometry(QRect(0, 90, 136, 209));
        springButton = new QRadioButton(ConnectorTypeBox);
        springButton->setObjectName(QString::fromUtf8("springButton"));
        springButton->setGeometry(QRect(10, 30, 117, 21));
        flexibleRodButton = new QRadioButton(ConnectorTypeBox);
        flexibleRodButton->setObjectName(QString::fromUtf8("flexibleRodButton"));
        flexibleRodButton->setGeometry(QRect(10, 70, 117, 21));
        ropeButton = new QRadioButton(ConnectorTypeBox);
        ropeButton->setObjectName(QString::fromUtf8("ropeButton"));
        ropeButton->setGeometry(QRect(10, 90, 117, 21));
        rigidRodButton = new QRadioButton(ConnectorTypeBox);
        rigidRodButton->setObjectName(QString::fromUtf8("rigidRodButton"));
        rigidRodButton->setGeometry(QRect(10, 50, 117, 21));

        horizontalLayout_2->addWidget(clickFunctionBox);

        particleSettingsBox = new QGroupBox(horizontalLayoutWidget_2);
        particleSettingsBox->setObjectName(QString::fromUtf8("particleSettingsBox"));
        massLabel = new QLabel(particleSettingsBox);
        massLabel->setObjectName(QString::fromUtf8("massLabel"));
        massLabel->setGeometry(QRect(0, 30, 41, 21));
        maxSpringDistLabel = new QLabel(particleSettingsBox);
        maxSpringDistLabel->setObjectName(QString::fromUtf8("maxSpringDistLabel"));
        maxSpringDistLabel->setGeometry(QRect(0, 50, 111, 21));
        massEdit = new QLineEdit(particleSettingsBox);
        massEdit->setObjectName(QString::fromUtf8("massEdit"));
        massEdit->setGeometry(QRect(70, 30, 51, 21));
        maxSpringDistEdit = new QLineEdit(particleSettingsBox);
        maxSpringDistEdit->setObjectName(QString::fromUtf8("maxSpringDistEdit"));
        maxSpringDistEdit->setGeometry(QRect(70, 70, 51, 21));
        rodSettingsBox = new QGroupBox(particleSettingsBox);
        rodSettingsBox->setObjectName(QString::fromUtf8("rodSettingsBox"));
        rodSettingsBox->setGeometry(QRect(0, 90, 135, 209));
        rodDensityLabel = new QLabel(rodSettingsBox);
        rodDensityLabel->setObjectName(QString::fromUtf8("rodDensityLabel"));
        rodDensityLabel->setGeometry(QRect(0, 30, 61, 21));
        rodStretchLabel = new QLabel(rodSettingsBox);
        rodStretchLabel->setObjectName(QString::fromUtf8("rodStretchLabel"));
        rodStretchLabel->setGeometry(QRect(0, 50, 61, 21));
        rodDensityEdit = new QLineEdit(rodSettingsBox);
        rodDensityEdit->setObjectName(QString::fromUtf8("rodDensityEdit"));
        rodDensityEdit->setGeometry(QRect(70, 30, 51, 21));
        rodStretchEdit = new QLineEdit(rodSettingsBox);
        rodStretchEdit->setObjectName(QString::fromUtf8("rodStretchEdit"));
        rodStretchEdit->setGeometry(QRect(70, 50, 51, 21));
        rodBendLabel = new QLabel(rodSettingsBox);
        rodBendLabel->setObjectName(QString::fromUtf8("rodBendLabel"));
        rodBendLabel->setGeometry(QRect(0, 70, 61, 21));
        rodBendEdit = new QLineEdit(rodSettingsBox);
        rodBendEdit->setObjectName(QString::fromUtf8("rodBendEdit"));
        rodBendEdit->setGeometry(QRect(70, 70, 51, 21));
        rodSegmentsEdit = new QLineEdit(rodSettingsBox);
        rodSegmentsEdit->setObjectName(QString::fromUtf8("rodSegmentsEdit"));
        rodSegmentsEdit->setGeometry(QRect(70, 90, 51, 21));
        rodSegmentsLabel = new QLabel(rodSettingsBox);
        rodSegmentsLabel->setObjectName(QString::fromUtf8("rodSegmentsLabel"));
        rodSegmentsLabel->setGeometry(QRect(0, 90, 61, 21));

        horizontalLayout_2->addWidget(particleSettingsBox);

        sawSettingsBox = new QGroupBox(horizontalLayoutWidget_2);
        sawSettingsBox->setObjectName(QString::fromUtf8("sawSettingsBox"));
        radiusLabel = new QLabel(sawSettingsBox);
        radiusLabel->setObjectName(QString::fromUtf8("radiusLabel"));
        radiusLabel->setGeometry(QRect(0, 30, 51, 21));
        radiusEdit = new QLineEdit(sawSettingsBox);
        radiusEdit->setObjectName(QString::fromUtf8("radiusEdit"));
        radiusEdit->setGeometry(QRect(80, 30, 51, 21));
        ropeSettingsBox = new QGroupBox(sawSettingsBox);
        ropeSettingsBox->setObjectName(QString::fromUtf8("ropeSettingsBox"));
        ropeSettingsBox->setGeometry(QRect(0, 90, 136, 209));
        ropeDensityLabel = new QLabel(ropeSettingsBox);
        ropeDensityLabel->setObjectName(QString::fromUtf8("ropeDensityLabel"));
        ropeDensityLabel->setGeometry(QRect(0, 30, 51, 21));
        ropeDensityEdit = new QLineEdit(ropeSettingsBox);
        ropeDensityEdit->setObjectName(QString::fromUtf8("ropeDensityEdit"));
        ropeDensityEdit->setGeometry(QRect(80, 30, 51, 21));
        ropeBendkLabel = new QLabel(ropeSettingsBox);
        ropeBendkLabel->setObjectName(QString::fromUtf8("ropeBendkLabel"));
        ropeBendkLabel->setGeometry(QRect(0, 50, 51, 21));
        ropeBendEdit = new QLineEdit(ropeSettingsBox);
        ropeBendEdit->setObjectName(QString::fromUtf8("ropeBendEdit"));
        ropeBendEdit->setGeometry(QRect(80, 50, 51, 21));
        ropeSegmentsLabel = new QLabel(ropeSettingsBox);
        ropeSegmentsLabel->setObjectName(QString::fromUtf8("ropeSegmentsLabel"));
        ropeSegmentsLabel->setGeometry(QRect(0, 70, 81, 21));
        ropeSegmentsEdit = new QLineEdit(ropeSettingsBox);
        ropeSegmentsEdit->setObjectName(QString::fromUtf8("ropeSegmentsEdit"));
        ropeSegmentsEdit->setGeometry(QRect(80, 70, 51, 21));

        horizontalLayout_2->addWidget(sawSettingsBox);


        verticalLayout->addWidget(uiOptionsBox);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1200, 25));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuScene = new QMenu(menuBar);
        menuScene->setObjectName(QString::fromUtf8("menuScene"));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuScene->menuAction());
        menuFile->addAction(actionExit);
        menuScene->addAction(actionReset);
        menuScene->addAction(actionReset_Everything);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Universe of Goo", 0, QApplication::UnicodeUTF8));
        actionExit->setText(QApplication::translate("MainWindow", "Exit", 0, QApplication::UnicodeUTF8));
        actionReset->setText(QApplication::translate("MainWindow", "Clear Scene", 0, QApplication::UnicodeUTF8));
        actionReset_Everything->setText(QApplication::translate("MainWindow", "Reset Everything", 0, QApplication::UnicodeUTF8));
        simOptionsBox->setTitle(QApplication::translate("MainWindow", "Simulation Options", 0, QApplication::UnicodeUTF8));
        ConstraintBox->setTitle(QApplication::translate("MainWindow", "Constraint Handling", 0, QApplication::UnicodeUTF8));
        penaltyForceButton->setText(QApplication::translate("MainWindow", "Penalty Force", 0, QApplication::UnicodeUTF8));
        stepAndProjectButton->setText(QApplication::translate("MainWindow", "Step and Project", 0, QApplication::UnicodeUTF8));
        lagrangeMultiplierButton->setText(QApplication::translate("MainWindow", "Lagrange Multiplier", 0, QApplication::UnicodeUTF8));
        startSimulationButton->setText(QApplication::translate("MainWindow", "Start Simulation", 0, QApplication::UnicodeUTF8));
        SimParametersBox->setTitle(QApplication::translate("MainWindow", "Parameters", 0, QApplication::UnicodeUTF8));
        timeStepLabel->setText(QApplication::translate("MainWindow", "Time Step:", 0, QApplication::UnicodeUTF8));
        newtonTolLabel->setText(QApplication::translate("MainWindow", "Newton Tolerance:", 0, QApplication::UnicodeUTF8));
        newtonMaxItersLabel->setText(QApplication::translate("MainWindow", "Newton Max Iters:", 0, QApplication::UnicodeUTF8));
        penaltyStiffness->setText(QApplication::translate("MainWindow", "Penalty Stiffness:", 0, QApplication::UnicodeUTF8));
        activeForcesBox->setTitle(QApplication::translate("MainWindow", "Active Forces", 0, QApplication::UnicodeUTF8));
        gravityCheckBox->setText(QApplication::translate("MainWindow", "Gravity", 0, QApplication::UnicodeUTF8));
        springsCheckBox->setText(QApplication::translate("MainWindow", "Springs", 0, QApplication::UnicodeUTF8));
        floorCheckBox->setText(QApplication::translate("MainWindow", "Floor", 0, QApplication::UnicodeUTF8));
        gravityGLabel->setText(QApplication::translate("MainWindow", "Acceleration:", 0, QApplication::UnicodeUTF8));
        springStiffnessLabel->setText(QApplication::translate("MainWindow", "Base Stiffness:", 0, QApplication::UnicodeUTF8));
        maxStrainLabel->setText(QApplication::translate("MainWindow", "Max Strain:", 0, QApplication::UnicodeUTF8));
        dampingStiffnessCheckBox->setText(QApplication::translate("MainWindow", "Viscous Damping", 0, QApplication::UnicodeUTF8));
        dampingStiffnessLabel->setText(QApplication::translate("MainWindow", "Stiffness:", 0, QApplication::UnicodeUTF8));
        elasticBendingCheckBox->setText(QApplication::translate("MainWindow", "Elastic Bending", 0, QApplication::UnicodeUTF8));
        uiOptionsBox->setTitle(QApplication::translate("MainWindow", "UI Options", 0, QApplication::UnicodeUTF8));
        clickFunctionBox->setTitle(QApplication::translate("MainWindow", "Click Function", 0, QApplication::UnicodeUTF8));
        addParticleButton->setText(QApplication::translate("MainWindow", "Add Particle", 0, QApplication::UnicodeUTF8));
        addSawButton->setText(QApplication::translate("MainWindow", "Add Saw", 0, QApplication::UnicodeUTF8));
        isFixedCheckBox->setText(QApplication::translate("MainWindow", "Is Fixed", 0, QApplication::UnicodeUTF8));
        ConnectorTypeBox->setTitle(QApplication::translate("MainWindow", "Connector Type", 0, QApplication::UnicodeUTF8));
        springButton->setText(QApplication::translate("MainWindow", "Spring", 0, QApplication::UnicodeUTF8));
        flexibleRodButton->setText(QApplication::translate("MainWindow", "FlexibleRod", 0, QApplication::UnicodeUTF8));
        ropeButton->setText(QApplication::translate("MainWindow", "Rope", 0, QApplication::UnicodeUTF8));
        rigidRodButton->setText(QApplication::translate("MainWindow", "Rigid Rod", 0, QApplication::UnicodeUTF8));
        particleSettingsBox->setTitle(QApplication::translate("MainWindow", "Particle Settings", 0, QApplication::UnicodeUTF8));
        massLabel->setText(QApplication::translate("MainWindow", "Mass:", 0, QApplication::UnicodeUTF8));
        maxSpringDistLabel->setText(QApplication::translate("MainWindow", "Max Spring Dist:", 0, QApplication::UnicodeUTF8));
        rodSettingsBox->setTitle(QApplication::translate("MainWindow", "Rod Settings", 0, QApplication::UnicodeUTF8));
        rodDensityLabel->setText(QApplication::translate("MainWindow", "Density:", 0, QApplication::UnicodeUTF8));
        rodStretchLabel->setText(QApplication::translate("MainWindow", "Stretch k", 0, QApplication::UnicodeUTF8));
        rodBendLabel->setText(QApplication::translate("MainWindow", "Bend k", 0, QApplication::UnicodeUTF8));
        rodBendEdit->setText(QString());
        rodSegmentsEdit->setText(QString());
        rodSegmentsLabel->setText(QApplication::translate("MainWindow", "Segments", 0, QApplication::UnicodeUTF8));
        sawSettingsBox->setTitle(QApplication::translate("MainWindow", "Saw Settings", 0, QApplication::UnicodeUTF8));
        radiusLabel->setText(QApplication::translate("MainWindow", "Radius:", 0, QApplication::UnicodeUTF8));
        ropeSettingsBox->setTitle(QApplication::translate("MainWindow", "Rope Settings", 0, QApplication::UnicodeUTF8));
        ropeDensityLabel->setText(QApplication::translate("MainWindow", "Density:", 0, QApplication::UnicodeUTF8));
        ropeBendkLabel->setText(QApplication::translate("MainWindow", "Bend k", 0, QApplication::UnicodeUTF8));
        ropeSegmentsLabel->setText(QApplication::translate("MainWindow", "Segments", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        menuScene->setTitle(QApplication::translate("MainWindow", "Scene", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H

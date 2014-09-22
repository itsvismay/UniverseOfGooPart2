/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created: Mon Sep 22 17:29:47 2014
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "mainwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MainWindow[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      39,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      19,   12,   11,   11, 0x0a,
      54,   11,   11,   11, 0x08,
      65,   11,   11,   11, 0x08,
      91,   11,   11,   11, 0x08,
     129,   11,   11,   11, 0x08,
     156,   11,   11,   11, 0x08,
     191,   11,   11,   11, 0x08,
     225,   11,   11,   11, 0x08,
     260,   11,   11,   11, 0x08,
     300,   11,   11,   11, 0x08,
     329,   11,   11,   11, 0x08,
     358,   11,   11,   11, 0x08,
     385,   11,   11,   11, 0x08,
     423,   11,   11,   11, 0x08,
     457,   11,   11,   11, 0x08,
     498,   11,   11,   11, 0x08,
     533,   11,   11,   11, 0x08,
     575,   11,   11,   11, 0x08,
     606,   11,   11,   11, 0x08,
     632,   11,   11,   11, 0x08,
     662,   11,   11,   11, 0x08,
     701,   11,   11,   11, 0x08,
     730,   11,   11,   11, 0x08,
     762,   11,   11,   11, 0x08,
     794,   11,   11,   11, 0x08,
     828,   11,   11,   11, 0x08,
     866,   11,   11,   11, 0x08,
     908,   11,   11,   11, 0x08,
     944,   11,   11,   11, 0x08,
     970,   11,   11,   11, 0x08,
     998,   11,   11,   11, 0x08,
    1029,   11,   11,   11, 0x08,
    1053,   11,   11,   11, 0x08,
    1089,   11,   11,   11, 0x08,
    1125,   11,   11,   11, 0x08,
    1158,   11,   11,   11, 0x08,
    1195,   11,   11,   11, 0x08,
    1232,   11,   11,   11, 0x08,
    1266,   11,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0\0params\0"
    "setUIFromParameters(SimParameters)\0"
    "updateGL()\0on_actionExit_triggered()\0"
    "on_actionReset_Everything_triggered()\0"
    "on_actionReset_triggered()\0"
    "on_startSimulationButton_clicked()\0"
    "on_timeStepEdit_editingFinished()\0"
    "on_newtonTolEdit_editingFinished()\0"
    "on_newtonMaxItersEdit_editingFinished()\0"
    "on_gravityCheckBox_clicked()\0"
    "on_springsCheckBox_clicked()\0"
    "on_floorCheckBox_clicked()\0"
    "on_dampingStiffnessCheckBox_clicked()\0"
    "on_gravityGEdit_editingFinished()\0"
    "on_springStiffnessEdit_editingFinished()\0"
    "on_maxStrainEdit_editingFinished()\0"
    "on_dampingStiffnessEdit_editingFinished()\0"
    "on_addParticleButton_clicked()\0"
    "on_addSawButton_clicked()\0"
    "on_massEdit_editingFinished()\0"
    "on_maxSpringDistEdit_editingFinished()\0"
    "on_isFixedCheckBox_clicked()\0"
    "on_radiusEdit_editingFinished()\0"
    "on_penaltyForceButton_clicked()\0"
    "on_stepAndProjectButton_clicked()\0"
    "on_lagrangeMultiplierButton_clicked()\0"
    "on_penaltyStiffnessEdit_editingFinished()\0"
    "on_elasticBendingCheckBox_clicked()\0"
    "on_springButton_clicked()\0"
    "on_rigidRodButton_clicked()\0"
    "on_flexibleRodButton_clicked()\0"
    "on_ropeButton_clicked()\0"
    "on_rodDensityEdit_editingFinished()\0"
    "on_rodStretchEdit_editingFinished()\0"
    "on_rodBendEdit_editingFinished()\0"
    "on_rodSegmentsEdit_editingFinished()\0"
    "on_ropeDensityEdit_editingFinished()\0"
    "on_ropeBendEdit_editingFinished()\0"
    "on_ropeSegmentsEdit_editingFinished()\0"
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        MainWindow *_t = static_cast<MainWindow *>(_o);
        switch (_id) {
        case 0: _t->setUIFromParameters((*reinterpret_cast< const SimParameters(*)>(_a[1]))); break;
        case 1: _t->updateGL(); break;
        case 2: _t->on_actionExit_triggered(); break;
        case 3: _t->on_actionReset_Everything_triggered(); break;
        case 4: _t->on_actionReset_triggered(); break;
        case 5: _t->on_startSimulationButton_clicked(); break;
        case 6: _t->on_timeStepEdit_editingFinished(); break;
        case 7: _t->on_newtonTolEdit_editingFinished(); break;
        case 8: _t->on_newtonMaxItersEdit_editingFinished(); break;
        case 9: _t->on_gravityCheckBox_clicked(); break;
        case 10: _t->on_springsCheckBox_clicked(); break;
        case 11: _t->on_floorCheckBox_clicked(); break;
        case 12: _t->on_dampingStiffnessCheckBox_clicked(); break;
        case 13: _t->on_gravityGEdit_editingFinished(); break;
        case 14: _t->on_springStiffnessEdit_editingFinished(); break;
        case 15: _t->on_maxStrainEdit_editingFinished(); break;
        case 16: _t->on_dampingStiffnessEdit_editingFinished(); break;
        case 17: _t->on_addParticleButton_clicked(); break;
        case 18: _t->on_addSawButton_clicked(); break;
        case 19: _t->on_massEdit_editingFinished(); break;
        case 20: _t->on_maxSpringDistEdit_editingFinished(); break;
        case 21: _t->on_isFixedCheckBox_clicked(); break;
        case 22: _t->on_radiusEdit_editingFinished(); break;
        case 23: _t->on_penaltyForceButton_clicked(); break;
        case 24: _t->on_stepAndProjectButton_clicked(); break;
        case 25: _t->on_lagrangeMultiplierButton_clicked(); break;
        case 26: _t->on_penaltyStiffnessEdit_editingFinished(); break;
        case 27: _t->on_elasticBendingCheckBox_clicked(); break;
        case 28: _t->on_springButton_clicked(); break;
        case 29: _t->on_rigidRodButton_clicked(); break;
        case 30: _t->on_flexibleRodButton_clicked(); break;
        case 31: _t->on_ropeButton_clicked(); break;
        case 32: _t->on_rodDensityEdit_editingFinished(); break;
        case 33: _t->on_rodStretchEdit_editingFinished(); break;
        case 34: _t->on_rodBendEdit_editingFinished(); break;
        case 35: _t->on_rodSegmentsEdit_editingFinished(); break;
        case 36: _t->on_ropeDensityEdit_editingFinished(); break;
        case 37: _t->on_ropeBendEdit_editingFinished(); break;
        case 38: _t->on_ropeSegmentsEdit_editingFinished(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData MainWindow::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow,
      qt_meta_data_MainWindow, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MainWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 39)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 39;
    }
    return _id;
}
QT_END_MOC_NAMESPACE

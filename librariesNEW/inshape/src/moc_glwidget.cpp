/****************************************************************************
** Meta object code from reading C++ file 'glwidget.h'
**
** Created: Tue 11. Sep 15:01:18 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "glwidget.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'glwidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.2.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_GLWidget[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
      38,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      10,    9,    9,    9, 0x0a,
      34,   27,    9,    9, 0x0a,
      48,    9,    9,    9, 0x0a,
      63,    9,    9,    9, 0x0a,
      81,    9,    9,    9, 0x0a,
     104,   98,    9,    9, 0x0a,
     122,    9,    9,    9, 0x0a,
     141,    9,    9,    9, 0x0a,
     158,    9,    9,    9, 0x0a,
     174,    9,    9,    9, 0x0a,
     187,    9,    9,    9, 0x0a,
     205,    9,    9,    9, 0x0a,
     214,    9,    9,    9, 0x0a,
     234,    9,    9,    9, 0x0a,
     246,    9,    9,    9, 0x0a,
     257,    9,    9,    9, 0x0a,
     268,    9,    9,    9, 0x0a,
     283,    9,    9,    9, 0x0a,
     295,    9,    9,    9, 0x0a,
     309,    9,    9,    9, 0x0a,
     333,    9,    9,    9, 0x0a,
     351,    9,    9,    9, 0x0a,
     365,    9,    9,    9, 0x0a,
     384,    9,    9,    9, 0x0a,
     395,    9,    9,    9, 0x0a,
     409,    9,    9,    9, 0x0a,
     422,    9,    9,    9, 0x0a,
     438,    9,    9,    9, 0x0a,
     461,  454,    9,    9, 0x0a,
     475,    9,    9,    9, 0x0a,
     491,    9,    9,    9, 0x0a,
     512,  505,    9,    9, 0x0a,
     527,  454,    9,    9, 0x0a,
     549,    9,    9,    9, 0x0a,
     568,    9,    9,    9, 0x0a,
     586,    9,    9,    9, 0x0a,
     599,    9,    9,    9, 0x0a,
     617,    9,    9,    9, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_GLWidget[] = {
    "GLWidget\0\0ToggleLighting()\0iScale\0"
    "SetScale(int)\0SetDrawBoxes()\0"
    "SetDrawPolyLine()\0ToggleDistance()\0"
    "iSize\0SetPointSize(int)\0ToggleSearchTree()\0"
    "ToggleRotation()\0ToggleCircles()\0"
    "ToggleGrid()\0ToggleAnimation()\0Rotate()\0"
    "ToggleTranslation()\0Translate()\0"
    "LogPoint()\0PrintLUB()\0PerformTests()\0"
    "PrintGrid()\0RebuildTree()\0"
    "ChangeSearchMethod(int)\0TransformCurves()\0"
    "AnimateStep()\0ToggleCorrection()\0"
    "DeBugLUB()\0CLUBToggled()\0LUBToggled()\0"
    "DecPointIndex()\0IncPointIndex()\0iDepth\0"
    "SetDepth(int)\0MergeVertices()\0"
    "BruteToggle()\0iPlane\0SetZPlane(int)\0"
    "SetAABBTreeDepth(int)\0SetSpatialMedian()\0"
    "SetObjectMedian()\0SetMinArea()\0"
    "SetObjBVHSimple()\0SetObjBVHMonotonic()\0"
};

const QMetaObject GLWidget::staticMetaObject = {
    { &QGLWidget::staticMetaObject, qt_meta_stringdata_GLWidget,
      qt_meta_data_GLWidget, 0 }
};

const QMetaObject *GLWidget::metaObject() const
{
    return &staticMetaObject;
}

void *GLWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GLWidget))
	return static_cast<void*>(const_cast< GLWidget*>(this));
    return QGLWidget::qt_metacast(_clname);
}

int GLWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ToggleLighting(); break;
        case 1: SetScale((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: SetDrawBoxes(); break;
        case 3: SetDrawPolyLine(); break;
        case 4: ToggleDistance(); break;
        case 5: SetPointSize((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 6: ToggleSearchTree(); break;
        case 7: ToggleRotation(); break;
        case 8: ToggleCircles(); break;
        case 9: ToggleGrid(); break;
        case 10: ToggleAnimation(); break;
        case 11: Rotate(); break;
        case 12: ToggleTranslation(); break;
        case 13: Translate(); break;
        case 14: LogPoint(); break;
        case 15: PrintLUB(); break;
        case 16: PerformTests(); break;
        case 17: PrintGrid(); break;
        case 18: RebuildTree(); break;
        case 19: ChangeSearchMethod((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 20: TransformCurves(); break;
        case 21: AnimateStep(); break;
        case 22: ToggleCorrection(); break;
        case 23: DeBugLUB(); break;
        case 24: CLUBToggled(); break;
        case 25: LUBToggled(); break;
        case 26: DecPointIndex(); break;
        case 27: IncPointIndex(); break;
        case 28: SetDepth((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 29: MergeVertices(); break;
        case 30: BruteToggle(); break;
        case 31: SetZPlane((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 32: SetAABBTreeDepth((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 33: SetSpatialMedian(); break;
        case 34: SetObjectMedian(); break;
        case 35: SetMinArea(); break;
        case 36: SetObjBVHSimple(); break;
        case 37: SetObjBVHMonotonic(); break;
        }
        _id -= 38;
    }
    return _id;
}

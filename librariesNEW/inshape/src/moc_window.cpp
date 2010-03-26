/****************************************************************************
** Meta object code from reading C++ file 'window.h'
**
** Created: Tue 11. Sep 15:01:18 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "window.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'window.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.2.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_Window[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
       8,    7,    7,    7, 0x05,

 // slots: signature, parameters, type, tag, flags
      29,    7,    7,    7, 0x0a,
      36,    7,    7,    7, 0x0a,
      45,    7,    7,    7, 0x0a,
      66,    7,    7,    7, 0x0a,
      84,    7,    7,    7, 0x0a,
     102,    7,    7,    7, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_Window[] = {
    "Window\0\0DebugButtonsChange()\0load()\0"
    "import()\0ToggleDebugButtons()\0"
    "DisableDButtons()\0LoadAnimatedObj()\0"
    "ToggleTreeSlider()\0"
};

const QMetaObject Window::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_Window,
      qt_meta_data_Window, 0 }
};

const QMetaObject *Window::metaObject() const
{
    return &staticMetaObject;
}

void *Window::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Window))
	return static_cast<void*>(const_cast< Window*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int Window::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: DebugButtonsChange(); break;
        case 1: load(); break;
        case 2: import(); break;
        case 3: ToggleDebugButtons(); break;
        case 4: DisableDButtons(); break;
        case 5: LoadAnimatedObj(); break;
        case 6: ToggleTreeSlider(); break;
        }
        _id -= 7;
    }
    return _id;
}

// SIGNAL 0
void Window::DebugButtonsChange()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

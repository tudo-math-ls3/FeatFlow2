/****************************************************************************
**
** Copyright (C) 2005-2007 Trolltech ASA. All rights reserved.
**
** This file is part of the example classes of the Qt Toolkit.
**
** This file may be used under the terms of the GNU General Public
** License version 2.0 as published by the Free Software Foundation
** and appearing in the file LICENSE.GPL included in the packaging of
** this file.  Please review the following information to ensure GNU
** General Public Licensing requirements will be met:
** http://www.trolltech.com/products/qt/opensource.html
**
** If you are unsure which license is appropriate for your use, please
** review the following information:
** http://www.trolltech.com/products/qt/licensing.html or contact the
** sales department at sales@trolltech.com.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
**
****************************************************************************/

#ifndef WINDOW_H
#define WINDOW_H

//===================================================
//				     Includes
//===================================================

#include <QMainWindow>
#include <QWidget>


//////////////* FORWARD DEFINITION *////////////

class QAction;
class QLabel;
class QMenu;
class QScrollArea;
class QSlider;
class QFrame;
class GLWidget;
class QGroupBox;
class QTabWidget;
class GLWidget;

class Window : public QMainWindow
{
    Q_OBJECT

public:
    Window();
    
	

public slots:
	void load();
	void import(void);
    void ToggleDebugButtons();
    void DisableDButtons();
	void LoadAnimatedObj();
	void ToggleTreeSlider();
	

signals:
	void DebugButtonsChange();	
	
private:


	GLWidget *window;

	QWidget    *centralWidget;
    QMenu      *fileMenu;
    QMenu      *helpMenu;
	QFrame     *glFrame;
    GLWidget   *glWidget;
	QAction    *exitAct;
	QGroupBox  *vGroupBox1;
	QGroupBox  *vGroupBox2;
	QGroupBox  *vGroupBox3;
	QGroupBox  *vGroupBox4;
	QGroupBox  *vGroupBox5;
	QGroupBox  *vGroupBox6;
	QGroupBox  *vGroupBox7;
    QTabWidget *TabWidget;
	QSlider    *wSlider;
	QSlider    *gSlider;
	QSlider    *ZSlider;

	bool m_bGB5State;
	bool m_bSLDRState;
	
	
	// Creates the group box for the particle options
	QGroupBox* CreateExclusiveGroup();
	QGroupBox* CreateExclusiveGroup3();

	QGroupBox* createVisualizationGroup(const char* szTitle);
	void createTabWidget();

	QGroupBox* CreateDebugButtons();

	// Creates the animation controls
	QGroupBox* CreateAnimationGroupBox(void);
	// creates the statistics control
	QGroupBox* CreatePushButtonGroup(void);
	//creates a horizontal slider
	QSlider* CreateSlider();
	//create AABBTree options
	QGroupBox* CreateExclusiveGroup2();


	

	
};

#endif

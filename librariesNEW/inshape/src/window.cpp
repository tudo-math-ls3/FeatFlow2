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

#include <QtGui>
#include <QPalette>
#include <QColor>
#include <fileParser.h>

#include <window.h>
#include <grid3d.h>
#include <glwidget.h>


Window::Window()
{

	centralWidget = new QWidget;
	setCentralWidget(centralWidget);

	fileMenu = menuBar()->addMenu(tr("&File"));
	
	//create an action
    exitAct = new QAction(tr("E&xit"), this);
    exitAct->setShortcut(tr("Ctrl+Q"));
    connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));

	QAction *loadAct = new QAction(tr("&Load"), this);
	loadAct->setShortcut(tr("Ctrl+L"));
	connect(loadAct, SIGNAL(triggered()), this, SLOT(load()));

	QAction *ImportAct = new QAction(tr("&Import Opennurbs curves"), this);
	ImportAct->setShortcut(tr("Ctrl+I"));
	connect(ImportAct, SIGNAL(triggered()), this, SLOT(import()));
	
	QAction *aniAct = new QAction(tr("&Import animated obj"), this);
	ImportAct->setShortcut(tr("Ctrl+A"));
	connect(aniAct, SIGNAL(triggered()), this, SLOT(LoadAnimatedObj()));	

	//add the action
	fileMenu->addAction(exitAct);
	fileMenu->addAction(loadAct);
	fileMenu->addAction(ImportAct);
	fileMenu->addAction(aniAct);

	
	//create a frame for the glWidget
	glFrame = new QFrame(this);
	glFrame->setFrameStyle(QFrame::Sunken | QFrame::Panel);
	glFrame->setLineWidth(2);
	
    glWidget = new GLWidget(glFrame);

	window = glWidget;

	QHBoxLayout *frameLayout = new QHBoxLayout(glFrame);
	frameLayout->addWidget(glWidget);

	vGroupBox1 = createVisualizationGroup("VisualizationOptions");

	vGroupBox2 = CreateExclusiveGroup();

	vGroupBox3 = CreateAnimationGroupBox();
	vGroupBox4 = CreatePushButtonGroup();
	vGroupBox5 = CreateDebugButtons();
	
	m_bGB5State = false;
	vGroupBox5->setEnabled(m_bGB5State);	

	vGroupBox6 = CreateExclusiveGroup2();
	vGroupBox7 = CreateExclusiveGroup3();

	createTabWidget();

    QHBoxLayout *mainLayout = new QHBoxLayout;
    mainLayout->addWidget(glFrame);
	mainLayout->addWidget(TabWidget);
	centralWidget->setLayout(mainLayout);

	m_bSLDRState = false;

    setWindowTitle(tr("InShape2D"));
}

QGroupBox *Window::createVisualizationGroup(const char* szTitle)
{
    QGroupBox *groupBox = new QGroupBox(tr(szTitle));

    QCheckBox *checkBox1 = new QCheckBox(tr("&Distance"));
	checkBox1->setChecked(true);
	connect(checkBox1,SIGNAL(clicked()),window,SLOT(ToggleDistance()));
    QCheckBox *checkBox2 = new QCheckBox(tr("&Show Distance Field"));
	connect(checkBox2,SIGNAL(clicked()),window,SLOT(ToggleGrid()));
    checkBox2->setChecked(false);

	QCheckBox *checkBox3 = new QCheckBox(tr("Newton Correction"));
	connect(checkBox3,SIGNAL(clicked()),window,SLOT(ToggleCorrection()));
	checkBox3->setChecked(false);

	QCheckBox *checkBox4 = new QCheckBox(tr("Draw Search Tree"));
	connect(checkBox4,SIGNAL(clicked()),window,SLOT(ToggleSearchTree()));
	connect(checkBox4,SIGNAL(clicked()),this,SLOT(ToggleTreeSlider()));
	checkBox4->setChecked(false);

	QCheckBox *checkBox5 = new QCheckBox(tr("Bounding Circles"));
	connect(checkBox5,SIGNAL(clicked()),window,SLOT(ToggleCircles()));
	checkBox5->setChecked(false);

	QCheckBox *checkBox6 = new QCheckBox(tr("Animation"));
	connect(checkBox6,SIGNAL(clicked()),window,SLOT(ToggleAnimation()));
	checkBox6->setChecked(false);

	QCheckBox *checkBox7 = new QCheckBox(tr("Lighting"));
	connect(checkBox7,SIGNAL(clicked()),window,SLOT(ToggleLighting()));
	checkBox7->setChecked(true);


	wSlider = CreateSlider();
	connect(wSlider, SIGNAL(valueChanged(int)), window, SLOT(SetDepth(int)));
	QLabel *label1 = new QLabel();
	QString qstring1("Object BVH Depth");
	label1->setText(qstring1);
	QString qstring2("Group BVH Depth");
	QLabel *label2 = new QLabel();
	label2->setText(qstring2);
	this->gSlider = CreateSlider();
	connect(gSlider, SIGNAL(valueChanged(int)), window, SLOT(SetAABBTreeDepth(int)));

	QString qstring3("Z Cut-Plane");
	QLabel *label3 = new QLabel();
	label3->setText(qstring3);
	this->ZSlider = CreateSlider();
	ZSlider->setRange(0, GRIDZ-1);
	ZSlider->setSingleStep(1);
	ZSlider->setValue(0);
	ZSlider->setEnabled(true);
	connect(ZSlider, SIGNAL(valueChanged(int)), window, SLOT(SetZPlane(int)));

  
    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(checkBox1);
    vbox->addWidget(checkBox2);
    vbox->addWidget(checkBox3);
    vbox->addWidget(checkBox4);
    vbox->addWidget(checkBox5);
    vbox->addWidget(checkBox6);
	vbox->addWidget(checkBox7);
	vbox->addWidget(label1);
	vbox->addWidget(wSlider);
	vbox->addWidget(label2);
	vbox->addWidget(gSlider);
	vbox->addWidget(label3);
	vbox->addWidget(ZSlider);


    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    return groupBox;
}

void Window::createTabWidget()
{
    TabWidget = new QTabWidget(this);
    TabWidget->setSizePolicy(QSizePolicy::Preferred,
							 QSizePolicy::Ignored);
	
	TabWidget->setFixedWidth(200);

	//==========Visualization tab============
    QWidget *tab1 = new QWidget(this);
    
    QVBoxLayout *tab1hbox = new QVBoxLayout;
    tab1hbox->setMargin(5);
    tab1hbox->addWidget(vGroupBox1);
    tab1->setLayout(tab1hbox);


	//==========Analysis tab============
    QWidget *tab2 = new QWidget(this);

	QPushButton *pushButton1 = new QPushButton(tr("&Animate"));
	connect(pushButton1,SIGNAL(clicked()),window,SLOT(AnimateStep()));
    
    QVBoxLayout *tab2hbox = new QVBoxLayout;
    tab2hbox->setMargin(5);
    tab2hbox->addWidget(vGroupBox2);
	tab2hbox->addWidget(vGroupBox5);
	tab2hbox->addWidget(vGroupBox3);
	tab2hbox->addWidget(pushButton1);
    tab2->setLayout(tab2hbox);

	//==========Options tab============
    QWidget *tab3 = new QWidget(this);
    
    QVBoxLayout *tab3hbox = new QVBoxLayout;
    tab3hbox->setMargin(5);
	tab3hbox->addWidget(vGroupBox6);
	tab3hbox->addWidget(vGroupBox7);
    tab3hbox->addWidget(vGroupBox4);
    tab3->setLayout(tab3hbox);

    TabWidget->addTab(tab1, tr("&Visualization"));
	TabWidget->addTab(tab2, tr("&Analysis"));
	TabWidget->addTab(tab3, tr("Op&tions"));


    
}

// Creates the group box for the particle options
QGroupBox* Window::CreateExclusiveGroup()
{
    QGroupBox *groupBox = new QGroupBox(tr("Search Method"));

    QRadioButton *radio1 = new QRadioButton(tr("&Coherency LBUB"));
	connect(radio1, SIGNAL(clicked()),window, SLOT(CLUBToggled()));
	connect(radio1, SIGNAL(clicked()),this, SLOT(ToggleDebugButtons()));
    QRadioButton *radio2 = new QRadioButton(tr("&LBUB"));
	connect(radio2, SIGNAL(clicked()),window, SLOT(LUBToggled()));
	connect(radio2, SIGNAL(clicked()),this, SLOT(ToggleDebugButtons()));
    QRadioButton *radio3 = new QRadioButton(tr("&DebugLUB"));
	connect(radio3, SIGNAL(clicked()),window, SLOT(DeBugLUB()));
	connect(radio3, SIGNAL(clicked()),this, SLOT(ToggleDebugButtons()));
    QRadioButton *radio4 = new QRadioButton(tr("&Brute Search"));
	connect(radio4, SIGNAL(clicked()),window, SLOT(BruteToggle()));
	connect(radio4, SIGNAL(clicked()),this, SLOT(ToggleDebugButtons()));
	

    radio1->setChecked(true);

    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(radio1);
    vbox->addWidget(radio2);
    vbox->addWidget(radio3);
	vbox->addWidget(radio4);
    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    return groupBox;
}
// Creates the animation controls
QGroupBox* Window::CreateAnimationGroupBox(void)
{
    QGroupBox *groupBox = new QGroupBox(tr("Animation Control"));

    QCheckBox *checkBox1 = new QCheckBox(tr("&Translation"));
    QCheckBox *checkBox2 = new QCheckBox(tr("&Rotation"));
  
    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(checkBox1);
    vbox->addWidget(checkBox2);

	groupBox->setLayout(vbox);

    return groupBox;
}

// creates the statistics control
QGroupBox* Window::CreatePushButtonGroup(void)
{
   QGroupBox *groupBox = new QGroupBox(tr("&Misc Options"));
 

    QPushButton *pushButton0 = new QPushButton(tr("&Benchmark"));
	connect(pushButton0,SIGNAL(clicked()),window,SLOT(PerformTests()));
    QPushButton *pushButton1 = new QPushButton(tr("&Merge Vertices"));
	connect(pushButton1,SIGNAL(clicked()),window,SLOT(MergeVertices()));
	
    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(pushButton0);
    vbox->addWidget(pushButton1);
    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    return groupBox;
}

void Window::load(void)
{

	if(!window->curve.empty())
	{
		for(int i = 0; i < window->m_nNumObjects; i++)
			delete window->curve[i];

		window->curve.clear();
	}//end if
		
	if(window->pACurve)
	{
		for(int i = 0; i < window->m_nNumObjects; i++)
			delete window->ApproxCurves[i];

		window->ApproxCurves.clear();
	}//end if

	QString s( QFileDialog::getOpenFileName(this,
                        tr("Choose a file name"), ".",
                        tr("TEXT (*.txt)")));


	if(s.isEmpty())
	{
		return;
	}//end if

	readNURBSDataFromFile(s.toLatin1(), window->curve);	
	
	window->m_nNumObjects = (int)window->curve.size();

	for(int k = 0; k < window->m_nNumObjects; k++)
	{
		window->curve[k]->genBoxes();
		window->curve[k]->scale(0.75);
		window->curve[k]->subdivideBoxes();
		window->curve[k]->subdivideBoxes();
		window->curve[k]->subdivideBoxes();
		CApproxCurve *myCurve = new CApproxCurve( window->curve[k], 150, static_cast<short>(k) );
		window->ApproxCurves.push_back(myCurve);
	}
	
	cout <<"number of particles: "<<window->m_nNumObjects<<endl;
	
	if(window->m_nNumObjects > 0)
		window->pACurve = window->ApproxCurves[0];
	
	/* position particles in the simulation domain */
	int o = 0;
	for(int k = 0; k < 4 && window->m_nNumObjects > 1; k++)
		window->ApproxCurves[k]->TranslateCurveByVector(VECTOR2(0.2, 0.5));
	
	for(int k = 0; k < 4 && window->m_nNumObjects > 1; k++)
		window->curve[k]->translateCurve(VECTOR2(0.2,0.5));
	
	VECTOR2 vTranslation(0.2,0.5);
	for(int i = 4; i < window->m_nNumObjects; i++)
	{
		if((i % 4) == 0)
			vTranslation.x += 0.6;
				
		window->ApproxCurves[i]->TranslateCurveByVector(vTranslation);
		
	}//end for	


	/* build tree from approximation curves */
	if(window->m_nNumObjects > 1)
	{
		PairIters pair(window->ApproxCurves.begin(), window->ApproxCurves.end());
		window->AABBTree.BuildTree(pair, CAABBTree::MEDIAN_CUT);  //myTree = new CAABBTree(window->ApproxCurves,CAABBTree::MEDIAN_CUT);
	}//end if

}//end load

void Window::import(void)
{

	QString s( QFileDialog::getOpenFileName(this,
                        tr("Choose a file name"), ".",
                        tr("3dm (*.3dm)")));


	if(s.isEmpty())
	{
		return;
	}//end if

	ON::Begin();

	ONX_Model model;

	const char *sFileName = s.toLatin1();

    // open file containing opennurbs archive
    FILE* archive_fp = ON::OpenFile( s.toLatin1(), "rb");
    if ( !archive_fp ) 
    {
     
     cerr<<"Could not open archive "<<sFileName<<endl;
	 exit(0);
    }//end if


	const ON_NurbsCurve* curve = 0;

	ON_BinaryFile archive(ON::read3dm, archive_fp);

    bool rc = model.Read(archive);

    // close the file
    ON::CloseFile( archive_fp );

	


	if(rc)
		cout<<"Archive successfully read."<<endl;
	else
		cout<<"Error during reading"<<endl;

	if(model.IsValid())
		cout<<"Model is valid"<<endl;
	else
		cout<<"Model is not valid"<<endl;


	int nObj = model.m_object_table.Count();

	printf("Number of objects %d\n",model.m_object_table.Count());

	window->CleanUp();	

	window->m_nNumObjects = nObj;

	for(int i = 0; i < nObj; i++)
	{

    	const ONX_Model_Object mo = model.m_object_table[i];
	    curve = ON_NurbsCurve::Cast(mo.m_object);

		const ON_NurbsCurve &rcurve = *curve;

		CONXNurbs *tempCurve = new CONXNurbs(rcurve);
		if(window->m_iObjBVH == GLWidget::SIMPLE)
		{
			tempCurve->CreateSimpleBVH();
		}
		else
		{
			tempCurve->genBoxes(2000);
		}
		tempCurve->subdivideBoxes();
		tempCurve->subdivideBoxes();
		tempCurve->subdivideBoxes();
		window->nxnurbs.push_back(tempCurve);
		CApproxCurve *myCurve = new CApproxCurve( tempCurve, 150, static_cast<short>(i) );
		window->ApproxCurves.push_back(myCurve);
	}//end for
	

	model.Destroy();

	window->pACurve = window->ApproxCurves[0];

	/* build tree from approximation curves */
	if(window->m_nNumObjects > 1)
	{
		PairIters pair(window->ApproxCurves.begin(), window->ApproxCurves.end());
		window->AABBTree.BuildTree(pair, window->m_nBuildMeth);
	}//end if

	int nDepth = 0;

	if(window->m_nNumObjects > 1)
	{
		window->AABBTree.breadthFirstSearch();
		nDepth = window->AABBTree.DepthFirstSearch(window->AABBTree.GetRoot(),0);
		cout<<"Depth: "<<nDepth;
	}//end if

	this->gSlider->setRange(0,nDepth);
	gSlider->setValue(0);

}//end import


/*!
    \fn Window::CreateDebugButtons()
 */
QGroupBox* Window::CreateDebugButtons()
{
	QGroupBox *groupBox = new QGroupBox(tr("Debug Cont&rol"));

/// We create a buttons

	QPushButton *NextButton = new QPushButton(tr("&Next"));
	connect(NextButton,SIGNAL(clicked()),window,SLOT(IncPointIndex()));	
    //QPushButton *PrevButton = new QPushButton(tr("&Previous"));
	//QPushButton *flatButton = new QPushButton(tr("&Flat Button"));

/// Finally, we lay out the widgets vertically, and return the group box that we created:

	QVBoxLayout *vbox = new QVBoxLayout;
	vbox->addWidget(NextButton);
	//vbox->addWidget(PrevButton);
	//vbox->addWidget(flatButton);
	vbox->addStretch(1);
	groupBox->setLayout(vbox);
	return groupBox;	

}


/*!
    \fn Window::ToggleDebugButtons()
 */
void Window::ToggleDebugButtons()
{
	
	if(window->m_Meth == GLWidget::DEBUG)
	{
		m_bGB5State = true;	
		this->vGroupBox5->setEnabled(m_bGB5State);
		this->wSlider->setEnabled(m_bGB5State);
	}//end if
	else
	{
		m_bGB5State = false;	
		this->vGroupBox5->setEnabled(m_bGB5State);
		this->wSlider->setEnabled(m_bGB5State);
	}
}//ToggleDebugButtons


/*!
    \fn Window::DisableDButtons()
 */
void Window::DisableDButtons()
{
    /// @todo implement me
}

/*!
    \fn Window::CreateSlider()
 */
QSlider* Window::CreateSlider()
{
	QSlider *slider = new QSlider(Qt::Horizontal);
	slider->setRange(0, 3);
	slider->setSingleStep(1);
	slider->setValue(3);
	slider->setTickPosition(QSlider::TicksAbove);
	slider->setEnabled(false);
	return slider;
}//CreateSlider


/*!
    \fn Window::LoadAnimatedObj()
 */
void Window::LoadAnimatedObj()
{
	
	if(!window->curve.empty())
	{
		for(int i = 0; i < window->m_nNumObjects; i++)
			delete window->curve[i];

		window->curve.clear();
	}//end if
		
	if(window->pACurve)
	{
		for(int i = 0; i < window->m_nNumObjects; i++)
			delete window->ApproxCurves[i];

		window->ApproxCurves.clear();
	}//end if

	if(window->AABBTree.IsInitialized())
	{
		window->AABBTree.DestroyTree();
	}

	QString s( QFileDialog::getOpenFileName(this,
			   tr("Choose a file name"), ".",
			   tr("3dm (*.3dm)")));


	if(s.isEmpty())
	{
		return;
	}//end if

	ON::Begin();

	ONX_Model model;

	const char *sFileName = s.toLatin1();

    // open file containing opennurbs archive
	FILE* archive_fp = ON::OpenFile( s.toLatin1(), "rb");
	if ( !archive_fp ) 
	{

     
		cerr<<"Could not open archive "<<sFileName<<endl;
		exit(0);
	}//end if


	const ON_NurbsCurve* curve = 0;


	ON_BinaryFile archive(ON::read3dm, archive_fp);

	bool rc = model.Read(archive);

    // close the file
	ON::CloseFile( archive_fp );

	if(rc)
		cout<<"Archive successfully read."<<endl;
	else
		cout<<"Error during reading"<<endl;

	if(model.IsValid())
		cout<<"Model is valid"<<endl;
	else
		cout<<"Model is not valid"<<endl;

	int nObj = model.m_object_table.Count();

	printf("Number of objects %d\n",model.m_object_table.Count());

	window->m_nNumObjects = 1;

	for(int i = 0; i < nObj; i++)
	{

		const ONX_Model_Object mo = model.m_object_table[i];
		curve = ON_NurbsCurve::Cast(mo.m_object);

		const ON_NurbsCurve &rcurve = *curve;

		CONXNurbs *tempCurve = new CONXNurbs(rcurve);

		tempCurve->genBoxes(2000);
		tempCurve->subdivideBoxes();
		tempCurve->subdivideBoxes();
		tempCurve->subdivideBoxes();
		window->nxnurbs.push_back(tempCurve);
		CApproxCurve *myCurve = new CApproxCurve( tempCurve, 150, static_cast<short>(i) );
		window->ApproxCurves.push_back(myCurve);
	}//end for
	
	window->pACurve = window->ApproxCurves[0];
			
	
	for(int j = 1; j <(int)window->ApproxCurves.size(); j++)
	{
		window->ApproxCurves[j]->SetCenter(window->pACurve->GetCenter());
	}//end for

	model.Destroy();
	
}//end LoadAnimatedObj

void Window::ToggleTreeSlider(void)
{
	if(window->m_bSearchTree)
	{
		m_bSLDRState = true;	
		this->gSlider->setEnabled(m_bSLDRState);
	}//end if
	else
	{
		m_bSLDRState = false;	
		this->gSlider->setEnabled(m_bSLDRState);
	}
}

QGroupBox* Window::CreateExclusiveGroup2()
{
    QGroupBox *groupBox = new QGroupBox(tr("BVH build method"));

    QRadioButton *radio1 = new QRadioButton(tr("&Spatial median"));
	connect(radio1, SIGNAL(clicked()),window, SLOT(SetSpatialMedian()));
    QRadioButton *radio2 = new QRadioButton(tr("&Object median"));
	connect(radio2, SIGNAL(clicked()),window, SLOT(SetObjectMedian()));
    QRadioButton *radio3 = new QRadioButton(tr("&Min Area"));
	connect(radio3, SIGNAL(clicked()),window, SLOT(SetMinArea()));

	QPushButton *pushButton0 = new QPushButton(tr("Re&built Search Tree"));
	connect(pushButton0,SIGNAL(clicked()),window,SLOT(RebuildTree()));

    radio1->setChecked(true);

    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(radio1);
    vbox->addWidget(radio2);
    vbox->addWidget(radio3);
	vbox->addWidget(pushButton0);
    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    return groupBox;
}

QGroupBox* Window::CreateExclusiveGroup3()
{
    QGroupBox *groupBox = new QGroupBox(tr("Object BVH build method"));

    QRadioButton *radio1 = new QRadioButton(tr("&Simple subdivision"));
	connect(radio1, SIGNAL(clicked()),window, SLOT(SetObjBVHSimple()));
    QRadioButton *radio2 = new QRadioButton(tr("&MP heuristic"));
	connect(radio2, SIGNAL(clicked()),window, SLOT(SetObjBVHMonotonic()));

    radio2->setChecked(true);

    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(radio1);
    vbox->addWidget(radio2);
    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    return groupBox;
}

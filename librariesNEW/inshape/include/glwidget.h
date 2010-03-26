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

#ifndef GLWIDGET_H
#define GLWIDGET_H


#include <QGLWidget>
#include <Qt>
#include <vector>
#include <QTimer>
#include "nurbs.h"
#include "distops.h"
#include "fileParser.h"
#include "aabbtree.h"
#include "approxcurve.h"
#include "rootfinder.h"
#include "grid.h"
#include "colldetect.h"
#include "obb2.h"
//#include "quadtree.h"

#define NUM_THREADS 2

using namespace std;

class window;

//===================================================
//					DEFINES
//===================================================


#define MAX_TEXTURES  100

#define MAX_PARTICLES 1024

#define MAX_NODES 2048

class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent = 0);
    virtual ~GLWidget();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    void MultithreadNewton();
    void MakeTextures();
    void CalcMaxDistance();
	void CollisionResponse(CollPair &Pair);

	
public slots:

	void ToggleLighting();
	void SetScale(int iScale);
	void SetDrawBoxes();
	void SetDrawPolyLine();
	void ToggleDistance();
	void SetPointSize(int iSize);
	void ToggleSearchTree();
	void ToggleRotation();
	void ToggleCircles();
	void ToggleGrid();
	void ToggleAnimation();
	void Rotate();
	void ToggleTranslation();
	void Translate();
	void LogPoint();
	void PrintLUB();
	void PerformTests();
	void PrintGrid();
	void RebuildTree();
	void ChangeSearchMethod(int);
	void TransformCurves(void);
    void AnimateStep();
    void ToggleCorrection();
    void DeBugLUB();
    void CLUBToggled();
    void LUBToggled();
    void DecPointIndex();
    void IncPointIndex();
    void SetDepth(int iDepth);
	void MergeVertices();
	void BruteToggle();
	void SetZPlane(int iPlane);
	void SetAABBTreeDepth(int iDepth);
	// sets build method to CAABBTree::SPATIAL_MEDIAN
	void SetSpatialMedian();
	// Sets the build method to CAABBTree::MEDIAN_CUT
	void SetObjectMedian();
	// Sets the built method to CAABBTree::MIN_AREA
	void SetMinArea();
	// Sets object build method to simple
	void SetObjBVHSimple();
	// Sets object build method to mp heuristic
	void SetObjBVHMonotonic();
	
protected:
    virtual void initializeGL();
    virtual void paintGL();
    virtual void resizeGL(int w, int h);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);


	
	
    GLuint makeObject();
    void quad(GLdouble x1, GLdouble y1, GLdouble x2, GLdouble y2,
              GLdouble x3, GLdouble y3, GLdouble x4, GLdouble y4);
    void extrude(GLdouble x1, GLdouble y1, GLdouble x2, GLdouble y2);
    void normalizeAngle(int *angle);

	//the Grid object
	CGrid g_Grid;

	//handle to a display list
    GLuint object;

	//maximum distance in the grid
	Real maxDist;
 

	void CoherencySearchMultiple(void);
	void CoherencySearch(void);
	void CoherencySearchSingle(void);
	void LowerBoundUpperBound(void);
	DebugResult LowerBoundUpperBoundSingle(void);
	void RegLubSingle(void);
	void ToggleDrawBoxes(void);
	void ToggleBoundingCircles(void);
	void CollisionCheck(CApproxCurve *pCurve, int n);
	void RebuiltTree(void);


	t_ThreadResource resources[NUM_THREADS];

	CONXNurbs *myNurbs;
	ON_NurbsCurve m_Curve;

///=====================================================	
///        		Light parameters
///      variables for light position, ambient, diffuse
///		 and specular
///=====================================================

	GLfloat light_position[4];
	GLfloat light_ambient[4]; 
	GLfloat light_diffuse[4]; 
	GLfloat light_specular[4];

///===================================================	
///        		 	
///			OpenGL parameters	
///
///===================================================

	GLfloat m_nFieldOfView;
		
///===================================================	
///        		 State variables	
///      These vars control various GUI states	
///===================================================
	
	int     m_Meth;
	int	    m_bCLBUB;
	int	    m_bLBUB;
	int     m_bDebugLUB;
	int  	m_bNewton;
	int		m_nBuildMeth;
	int		m_iObjBVH;
	bool   	m_bTranslation;
	bool	m_bDistanceField;
	bool    m_bCircles;
	bool	m_bSearchTree;
	bool    m_bDrawBoxes;
	float   m_TransZ;
	bool	m_bDrawPolyLine;
	bool    m_bDistance;
	bool    m_bTimer;
	bool    m_NewtonCorrection;
	bool    m_bLighting;

	//controls for 3D visisualization
	int m_iDepth;
	int m_iSTDepth;
	int m_iZPlane;

	//the bounding box of the 2D Grid
	COBB boundingBox;

	//storage for textures
	GLuint g_Textures[MAX_TEXTURES];

	//the velocity vector components of the particle velocities
	double dVelX[MAX_PARTICLES];
	double dVelY[MAX_PARTICLES];

	//the current view translation in X and Y direction
	GLfloat m_TransX, m_TransY;

	//the current rotation values around the X and Y
	GLfloat m_RotX, m_RotY;
	float m_oldy, m_oldx;

	//vector to store the Nurbs curves
	vector<CNurbs*> curve;

	//vector to store the opennurbs curves
	vector<CONXNurbs*> nxnurbs;

	//vector of approximation curves
	vector<CApproxCurve*> ApproxCurves;
	CApproxCurve *pACurve;

	//the number of particles
	int		m_nNumObjects;

	//the aabbtree search data structure
	CAABBTree AABBTree;
 
	//bool init;

	//variables for a debug visualization
	int m_CountI;
	int m_CountJ;

	int m_nWidth ;
	int m_nHeight;



	int iteration;

	// frees memory used by curves and data structures
	void CleanUp(void);

	void MultithreadCoherency(void);
	void BruteSearchSingle(void);	
	void InitializeVelocities(void);

	QTimer animationTimer;
	void PerturbMesh(void);
public:
    QPixmap m_Images[MAX_TEXTURES];
	


///===================================================	
///        		 Enumeration	
///      The enumerations identify program options
///===================================================
	enum
	{
		CLBUB,
		LBUB,
		DEBUG,
		NEWTON,
		BRUTE
	};
	
	enum
	{
		SIMPLE,
		MPHEURISTIC
	};

///===================================================	
///        		  FRIEND CLASSES
///      
///===================================================
	friend class Window;

///========================================================================================	
///      Protected member variables: - m_bMedAxis 		 
///									 - m_bDrawVoronoiEdges	
///========================================================================================
protected:
    bool m_bMedAxis;
    bool m_bDrawVoronoiEdges;

public:
	void GetFictitiousBoundaryValues(void);
	void GetFictitiousBoundaryValuesMultiple(void);
};

#endif

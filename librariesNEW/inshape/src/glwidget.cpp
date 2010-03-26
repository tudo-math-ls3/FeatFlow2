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

#pragma warning(disable:4005)


#include <QtGui>
#include <QtOpenGL>
#include <math.h>
#include <limits>
#include <glwidget.h>
#include <colldetect.h>
#include <matrix4x4.h>
#include <matrix3x3.h>
#include <ObjLoader.h>
#include <model3d.h>
#include <grid3d.h>
#include <feat2ops.h>
#include <distops3.h>
#include <cbvhnode.h>
#include <visualdebug.h>
#include <quaternion.h>
#include <fbvaluesthread.h>

#if defined WIN32
#include <windows.h>
#include <mmsystem.h>
#endif


#ifdef WIN32

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#endif


using namespace std;




DebugSEC info;
CVisualDebug vd;
int g_Counter = 0;



#ifdef WIN32
class BrtDst
{
public:

	int m_End;
	VECTOR3 *vVertices;

	BrtDst(VECTOR3 *vVertices_, int nNumP) : m_End(nNumP), vVertices(vVertices_) {};

	void operator() (const blocked_range<int> &r) const
	{
		CDistOps3 op;
		for(int i = 0; i < m_End ; i++)
		{
			Real dist = op.BruteForceDistance(model, grid3d(i));
			grid3d.m_Grid[i].distance = dist;
		}//end for
	}
};

class BrtInnr
{
public:

	int m_End;

	BrtInnr(int nNumP) : m_End(nNumP){};

	void operator() (const blocked_range<int> &r) const
	{
		CDistOps3 op;
		for(int i = 0; i < m_End ; i++)
		{
			bool inner = op.BruteForceInnerPoints(model,grid3d(i));
			grid3d.m_Grid[i].Inside = inner;
		}//end for
	}//end operator
};

class FBCValues
{
public:

	int m_End;

	CApproxCurve *m_pCurve;

	CGrid *m_pGrid;

	FBCValues(int nNumP, CApproxCurve *pCurve, CGrid *grid) : m_End(nNumP)
	{
		m_pCurve = pCurve;
		m_pGrid = grid;
	};

	void operator() (const blocked_range<int> &r) const
	{

		//class for fictitious boundary operations
		CFeat2Ops op;

		for(int i = 0; i < m_End ; i++)
		{
			Real res = op.IsFictitiousBoundarySimple(m_pCurve, m_pGrid->m_pVertices[i]);
			m_pGrid->SetDistance(i, res);
		}//end for
	}//end operator
};
#endif




GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(parent)
{

	//initialize the light parameters
	light_position[0] = 0.0;
	light_position[1] = 0.0;
	light_position[2] = 10.0;
	light_position[3] = 0.0;

	light_ambient[0] = 0.0;
	light_ambient[1] = 0.0;
	light_ambient[2] = 0.0;
	light_ambient[3] = 1.0;

	light_diffuse[0] = 1.0;
	light_diffuse[1] = 1.0;
	light_diffuse[2] = 1.0;
	light_diffuse[3] = 1.0;

	
	light_specular[0] = 1.0;
	light_specular[1] = 1.0;
	light_specular[2] = 1.0;
	light_specular[3] = 1.0;

	//translation parameters
	m_TransZ = -10.0f;
	m_TransX = 0.0f;
	m_TransY = 0.0f;

	m_RotX   = 0.0f;
	m_RotY   = 0.0f;

	InitializeVelocities();

	//CMatrix3x3<float> mat0(1,2,3,4,5,6,7,8,9);
//	CMatrix3x3<float> mat0(0,1,2,3,2,1,1,1,0);
	//CMatrix3x3<float> mat1(10,16,6,14,12,20,9,7,11);

	/*cout<<mat0<<endl;

	cout<<mat0.Determinate()<<endl;*/



	//cout<<mat1<<endl;

	//CMatrix3x3<float> mat2;

	//mat2 = mat0 * mat1;

	//cout<<mat2<<endl;

	//CMatrix4x4<float> mat(1,2,3,0,0,1,2,-2,-1,-1,3,2,1,1,2,0);

	//cout<<mat<<endl;

	//for(int i = 0; i < 4; i++)
	//	for(int j = 0; j <4; j++)
	//	{
	//		CMatrix3x3<float> subMatrix = mat.GetSubMatrix3x3(i,j);
	//		cout<<subMatrix<<endl;
	//	}
	//

	//cout<<mat.Determinate()<<endl;




	VECTOR3 vBL(-5.0,-5.0,5.0);
	VECTOR3 vTR(5.0,5.0,5.0);
	VECTOR3 vBackTR(5.0,5.0,-5.0);

	//grid3d.InitGrid(vBL, vTR, vBackTR);

	g_Grid.InitGrid(VECTOR2(0,0), VECTOR2(1.0,1.0)); 

	//PerturbMesh();
	//CObjLoader Loader;

	//Loader.ReadModelFromFile("ship.obj");
	//Loader.ReadModelFromFile("ship.obj");

	//model.SetVertices(Loader.GetVertices());
	//model.SetFaces(Loader.GetFaces());
	//model.SetNormals(Loader.m_pVertexNormals);
	//model.BuildTriangles();

	
	//CBVHSubdivider subdiv;
	//BVHTree.BuildTree(model.GetTriangles(), subdiv);
	//BVHTree.BuildTree(model, subdiv);

	//AABBTree3f BVHTree2(BVHTree);

	m_bDrawBoxes        = false;
	m_bDrawPolyLine     = false;
	m_bDistance         = true;
	m_bSearchTree	    = false;
	m_bCircles		    = false;
	m_bDistanceField    = false;
	m_bCLBUB	        = true;
	m_bLBUB			    = false;
	m_bNewton		    = false;
	m_bTimer		    = false;
	m_NewtonCorrection  = false;
	m_bMedAxis          = true;
	m_bDrawVoronoiEdges = true;
	m_bLighting		    = true;
	m_iZPlane           = 0;
	
	iteration = 0;

	m_iDepth = 3;
	m_iSTDepth = 0;
	m_iObjBVH = GLWidget::MPHEURISTIC;
	
   QTime midnight(0, 0, 0);
   srand(midnight.secsTo(QTime::currentTime()));
 	    
   connect(&animationTimer, SIGNAL(timeout()), this, SLOT(TransformCurves()));

	m_CountI = 0;
	m_CountJ = 0;

	myNurbs = NULL;

	pACurve  = NULL;

	/* domain bounding box */
	boundingBox.InitBox(VECTOR2(0.0,0.0), VECTOR2(10.0,10.0));
	
	m_Meth = CLBUB;

	m_nBuildMeth = CAABBTree::SPATIAL_MEDIAN;

}//constructor



GLWidget::~GLWidget()
{
    makeCurrent();

	if(myNurbs)
		delete myNurbs;

	CleanUp();

	ON::End();

}//deconstructor

QSize GLWidget::minimumSizeHint() const
{
    return QSize(50, 50);
}

QSize GLWidget::sizeHint() const
{
    return QSize(400, 400);
}

void GLWidget::initializeGL()
{
 	/* set OpenGL clear color to black */
	qglClearColor( Qt::white );

	MakeTextures();
	
	/* enable smooth shading */
	glShadeModel( GL_SMOOTH );

	float shine =  20.0f;

	//enable lighting
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT , light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE , light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glMaterialfv(GL_FRONT, GL_SPECULAR, light_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS,&shine);
	glEnable(GL_COLOR_MATERIAL);

	//enable z-Buffer
    glEnable(GL_DEPTH_TEST);

	//enable texturing
	//glEnable(GL_TEXTURE_2D);

	glDisable(GL_LIGHTING);
	
}//initializeGL

void GLWidget::paintGL()
{

	static float fRotation = 0.0f;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

	glTranslatef(m_TransX, m_TransY, m_TransZ);

	CQuaternionf qRotX;

	qRotX.AxisAngleToQuat(0,1,0, m_RotX);

	CMatrix4x4f matrix;

	qRotX.CreateMatrix(matrix);

	CQuaternionf qRotY;

	qRotY.AxisAngleToQuat(1,0,0, m_RotY);

	glMultMatrixf(matrix.Get());

	qRotY.CreateMatrix(matrix);

	glMultMatrixf(matrix.Get());


	if(m_bDistance && !ApproxCurves.empty() )
	{
		if(m_nNumObjects > 1)
		{
			switch(m_Meth) {
				
				case CLBUB:
					GetFictitiousBoundaryValuesMultiple();
					//CoherencySearchMultiple();
					break;
				
				case LBUB:
					LowerBoundUpperBound();
					break;
				
				case DEBUG:
					break;
	
				case NEWTON:
					break;
					
			}//end switch
			
		}//end if
		else
		{
			
			DebugResult dr;
			
			switch(m_Meth) {
				
				case CLBUB:
					//CoherencySearchSingle();
					GetFictitiousBoundaryValues();
					break;
				
				case LBUB:
					RegLubSingle();
					break;
				
				case DEBUG:
					dr = LowerBoundUpperBoundSingle();
					break;
	
				case BRUTE:
					BruteSearchSingle();
					break;
					
			}//end switch	
			
			if(m_Meth == DEBUG)
			{
				list<CBCNode*>::iterator iter;
				
				for(iter = dr.candidates[m_iDepth].begin(); iter != dr.candidates[m_iDepth].end(); iter++)
				{
					CBCNode *mIter = *iter;
					
					if(mIter == dr.pBest[m_iDepth]) continue;
					
					if(dr.pBestPoint)
					{
						if(dr.pBestPoint ==mIter)continue;
					}
					
					vd.DrawBCLeaf(mIter, 0,1,0,1);
					
				}//end for
				vd.DrawBCLeaf(dr.pBest[m_iDepth], 1,0,0,1);
				if(dr.pBestPoint != NULL)
					vd.DrawBCLeaf(dr.pBestPoint, 1,0,1,1);
				
				glColor3f(0,0,1.0);
				glBegin(GL_LINES);
				glVertex3f(dr.vPoint.x, dr.vPoint.y, 0.0);
				glVertex3f(dr.vPointOCurve.x, dr.vPointOCurve.y, 0.0);
				glColor3f(0,0.7,0);
				glVertex3f(dr.vPoint.x, dr.vPoint.y, 0.0);
				glVertex3f(dr.vUBounds[m_iDepth].x, dr.vUBounds[m_iDepth].y, 0.0);
				glEnd();
			}//end if		
			
		}//end else
	}//end if Distance

	/* draw approx curves */
	for(int i = 0; i < m_nNumObjects; i++)
	{
		vd.DrawApproxCurve(ApproxCurves[i], 0,0,0);
		//vd.DrawBoundingCircle(ApproxCurves[i]->GetCenter(), ApproxCurves[i]->Radius(), 1,0,0);
// 		vd.DrawApproxCurve(pACurve, 0,0,0);
// 		vd.DrawOBBox(this->nxnurbs[i]->BoundingBox(), 1, 0, 0);
	}//end for

	if(m_bCircles)
	{

		for(int i = 0; i < m_nNumObjects; i++)
		{
			for(int index = 0; index < ApproxCurves[i]->GetNumCircles(); index++)
			{
				vd.DrawCircleTree(ApproxCurves[i]->m_BCTree.GetChild(index), m_iDepth);
			}//end for
		}//end for

	}//end drawBoxes
	
	if(m_bDistanceField)
	{
		CalcMaxDistance();
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, g_Textures[0]);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glTexCoordPointer(2, GL_FLOAT, 0, g_Grid.m_pTCoords);
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(2, GL_FLOAT, 0, g_Grid.m_pVertices);
		vd.DrawtGrid(g_Grid, GRID_POINTS_X, GRID_POINTS_Y);
		glDisable(GL_TEXTURE_2D);
	}
	
	if(m_bSearchTree)
	{
		vd.DrawAABBTree(&AABBTree, m_iSTDepth);
	}//end

// 	if(m_bDrawVoronoiEdges)
// 	{
// 		for(int i = 0; i < m_nNumObjects; i++)
// 		{		
// 
// 			vd.DrawAdjList(ApproxCurves[i]->undirGraph);
// /*			if(ApproxCurves[i]->m_TestVertex)
// 				vd.DrawAxisCircle(ApproxCurves[i]->m_TestVertex, 1,0,0, info);*/
// 
// 		}//end for		
// 	}//end if

	//if(m_bDistance)
	//{
	//	this->LUBSingle3D();
	//	//COLUBSingle3D();
	//	//this->BruteDistance3D();
	//	//BruteInnerPoints();
	//	DoMax();
	//	//glDisable(GL_LIGHTING);
	//	//vd.DrawMeshPoints(grid3d);
	//	//glEnable(GL_LIGHTING);
	//}
	////vd.DrawGrid3D(grid3d);
	//glEnable(GL_TEXTURE_2D);
	//glBindTexture(GL_TEXTURE_2D, g_Textures[0]);
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	//vd.DrawCutPlaneZ(grid3d,m_iZPlane,maxDist);
	//glDisable(GL_TEXTURE_2D);



	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glColor4ub(255, 185, 0, 130);
	//glEnable(GL_BLEND);
	//
	//glEnable(GL_LIGHTING);
	//shader.EnableShader();
	//vd.DrawModel3D(model);
	//shader.DisableShader();
	//glDisable(GL_LIGHTING);
	//glDisable(GL_BLEND);

	//vd.DrawBVHTree(BVHTree);

	//if(!ApproxCurves.empty())
	//{
	//	ApproxCurves[0]->RotateCurve(fRotation);
	//	fRotation+= 0.001f;
	//}
	

	glPointSize(1);
}//end paintGL


void GLWidget::resizeGL(int w, int h)
{
	//establish viewport settings
 	glViewport(0,0,w,h);
	m_nWidth  = w;
	m_nHeight = h;
	
	//set a perspective projection
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluPerspective(60,(float)w/h,0.1,100);
	
	//reset the modelview matrix
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
	m_oldx = event->x();
	m_oldy = event->y();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
	
	if (event->buttons() & Qt::MidButton) {
 
		m_TransZ += (event->y()-m_oldy)*0.035;
	}
	else if (event->buttons() & Qt::RightButton)
	{
		m_TransX += 5*(event->x()-m_oldx)/m_nWidth;
		m_TransY -= 5*(event->y()-m_oldy)/m_nHeight;
	}
	else if (event->buttons() && Qt::LeftButton)
	{
		m_RotX += (event->x()-m_oldx)/5;
		m_RotY += (event->y()-m_oldy)/5;
	}//end else
	m_oldy = event->y();
	m_oldx = event->x();
		
	updateGL();
 
}

GLuint GLWidget::makeObject()
{
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);

    glBegin(GL_QUADS);

    GLdouble x1 = +0.06;
    GLdouble y1 = -0.14;
    GLdouble x2 = +0.14;
    GLdouble y2 = -0.06;
    GLdouble x3 = +0.08;
    GLdouble y3 = +0.00;
    GLdouble x4 = +0.30;
    GLdouble y4 = +0.22;

    quad(x1, y1, x2, y2, y2, x2, y1, x1);
    quad(x3, y3, x4, y4, y4, x4, y3, x3);

    extrude(x1, y1, x2, y2);
    extrude(x2, y2, y2, x2);
    extrude(y2, x2, y1, x1);
    extrude(y1, x1, x1, y1);
    extrude(x3, y3, x4, y4);
    extrude(x4, y4, y4, x4);
    extrude(y4, x4, y3, x3);

    const double Pi = 3.14159265358979323846;
    const int NumSectors = 200;

    for (int i = 0; i < NumSectors; ++i) {
        double angle1 = (i * 2 * Pi) / NumSectors;
        GLdouble x5 = 0.30 * sin(angle1);
        GLdouble y5 = 0.30 * cos(angle1);
        GLdouble x6 = 0.20 * sin(angle1);
        GLdouble y6 = 0.20 * cos(angle1);

        double angle2 = ((i + 1) * 2 * Pi) / NumSectors;
        GLdouble x7 = 0.20 * sin(angle2);
        GLdouble y7 = 0.20 * cos(angle2);
        GLdouble x8 = 0.30 * sin(angle2);
        GLdouble y8 = 0.30 * cos(angle2);

        quad(x5, y5, x6, y6, x7, y7, x8, y8);

        extrude(x6, y6, x7, y7);
        extrude(x8, y8, x5, y5);
    }

    glEnd();

    glEndList();
    return list;
}

void GLWidget::quad(GLdouble x1, GLdouble y1, GLdouble x2, GLdouble y2,
                    GLdouble x3, GLdouble y3, GLdouble x4, GLdouble y4)
{

    glVertex3d(x1, y1, -0.05);
    glVertex3d(x2, y2, -0.05);
    glVertex3d(x3, y3, -0.05);
    glVertex3d(x4, y4, -0.05);

    glVertex3d(x4, y4, +0.05);
    glVertex3d(x3, y3, +0.05);
    glVertex3d(x2, y2, +0.05);
    glVertex3d(x1, y1, +0.05);
}

void GLWidget::extrude(GLdouble x1, GLdouble y1, GLdouble x2, GLdouble y2)
{

    glVertex3d(x1, y1, +0.05);
    glVertex3d(x2, y2, +0.05);
    glVertex3d(x2, y2, -0.05);
    glVertex3d(x1, y1, -0.05);
}

void GLWidget::normalizeAngle(int *angle)
{
    while (*angle < 0)
        *angle += 360 * 16;
    while (*angle > 360 * 16)
        *angle -= 360 * 16;
}

void GLWidget::CoherencySearchMultiple(void)
{

	static int iter = 0;
	DistOps op;
	double dLUB = 1.7E+308;
	CRootFinder rFA;
	SearchResult lRes;
	lRes.dLUB = dLUB;

	VECTOR2 lastPoint(0,0);
	QTime timer;
	timer.start();
	
	int conv = 0;

	for(int i = 0; i < GRID_POINTS_X; i++)
	{
				
		for(int j = 0; j < GRID_POINTS_Y; j++)
		{
			list<CAABBNode*> lCandidates;
 			CApproxCurve *CurvePointer;
			AABBTree.breadthFirstSearchDistance(g_Grid.Point(i,j), lCandidates, dLUB);
			lRes = op.COLUBSearch(lCandidates, g_Grid.Point(i,j), CurvePointer, ApproxCurves, lRes);

			lastPoint = CurvePointer->GetPoint(lRes.data, lRes.point);

			lastPoint = VECTOR2::createVector(lastPoint, g_Grid.Point(i,(j+1)%GRID_POINTS_Y));
			
			g_Grid.SetDistance(i,j,lRes.res);
			
			lRes.dLUB = lastPoint.mag();
			
/*			if(m_NewtonCorrection)
			{
				
				CONXNurbs *Nurbs = nxnurbs[CurvePointer->ID()];
			
				double u = CurvePointer->m_dParamBuckets[result.data][result.point].dP;
				
				if(rFA.ONX_Newton(Nurbs, g_Grid.m_Grid[i][j].coord, u, 1e-7, 30))
				{
				
					conv++;
					VECTOR2 NearestPoint = Nurbs->CoxDeBoor(u);
					double d = VECTOR2::createVector(NearestPoint, g_Grid.m_Grid[i][j].coord).mag();
					if(result.res > 0)
					{
						g_Grid.m_Grid[i][j].newton = d;
				
					}
				}
			}*/			

		}//end for j
	}//end for i

 	cout<<"CPU time: "<<timer.elapsed()<<endl;
// 	cout<<"Converged: "<<conv<<endl;
	iter++;

}//CoherencySearchMultiple

void GLWidget::CoherencySearch(void)
{
}

void GLWidget::CoherencySearchSingle(void)
{

	VECTOR2 lastPoint(0,0);
	double dLUB = 1.7E+308;
	DistOps op;
	CRootFinder rFA;

	int iConv = 0;

	QTime timer;
	timer.start();	

	for(int i = 0; i < GRID_POINTS_X; i++)
	{
		for(int j = 0; j < GRID_POINTS_Y; j++)
		{
			list<CBCNode*> lBFS;
			CBCNode *nBest = NULL;
			SearchResult result = op.COLUBSearch(pACurve, g_Grid.Point(i,j), lBFS, nBest, dLUB);
			
			lastPoint = pACurve->GetPoint(result.data,result.point);

			lastPoint = VECTOR2::createVector(lastPoint, g_Grid.Point(i,(j+1)%GRID_POINTS_Y));
		
			dLUB = lastPoint.mag();

			g_Grid.SetDistance(i,j,result.res);
//			g_Grid.m_Grid[i][j].newton   = result.res;

			
			//double u = pACurve->m_dParamBuckets[result.data][result.point].dP;
			
		}//end for j
	}//end for i
	
	cout<<"CPU time: "<<timer.elapsed()<<endl;

}//end CoherencySingle

void GLWidget::LowerBoundUpperBound(void)
{
	

	DistOps op;

	VECTOR2 lastPoint(0,0);

	QTime timer;
	timer.start();
	
	int conv = 0;

	for(int i = 0; i < GRID_POINTS_X; i++)
	{
				
		for(int j = 0; j < GRID_POINTS_Y; j++)
		{
			list<CAABBNode*> lCandidates;
			CApproxCurve *CurvePointer;
			SearchResult result = op.AABBTreeSearch(&AABBTree, g_Grid.Point(i,j), CurvePointer, ApproxCurves); 
			
			g_Grid.SetDistance(i,j,result.res);
			
		}//end for j
	}//end for i

	cout<<" combined search time: "<<timer.elapsed()<<endl;
	
}

DebugResult GLWidget::LowerBoundUpperBoundSingle(void)
{
	VECTOR2 lastPoint(0,0);
	double dLUB = 1.7E+308;
	DistOps op;
	CRootFinder rFA;
	int iConv = 0;

	list<CBCNode*> lBFS;
	CBCNode *nBest = NULL;
	
	DebugResult result = op.DebugSearch(ApproxCurves[0], g_Grid.Point(m_CountI, m_CountJ), lBFS, nBest);

	g_Grid.SetDistance(m_CountI,m_CountJ,result.res);

	result.vPointOCurve = ApproxCurves[0]->GetPoint(result.data,result.point);

	//m_CountJ = (m_CountJ + 1) % GRID_POINTS_Y;

	return result;
}//end LowerBoundUpperBoundSingle

void GLWidget::RegLubSingle(void)
{
	
	DistOps op;

	VECTOR2 lastPoint(0,0);

	QTime timer;
	timer.start();
	
	int conv = 0;

	int i;
	
	for(i = 0; i < GRID_POINTS_X; i++)
	{
					
		for(int j = 0; j < GRID_POINTS_Y; j++)
		{
			list<CBCNode*> lBFS;
			CBCNode *nBest = NULL;
			SearchResult result = op.LUBSearch(pACurve, g_Grid.Point(i,j), lBFS, nBest);
			g_Grid.SetDistance(i,j,result.res);
		}//end for j
	}//end for i

	cout<<" CPU Time regular LUB: "<<timer.elapsed()<<endl;
	
}

void GLWidget::ToggleDrawBoxes(void)
{
}

void GLWidget::ToggleBoundingCircles(void)
{
}

void GLWidget::ToggleRotation(void)
{
}

void GLWidget::ChangeSearchMethod(int nID)
{

	switch(nID)
	{
		case 0:
			m_bCLBUB = true;
			m_bLBUB  = false;
			m_bNewton = false;
			break;
		case 1:
			m_bCLBUB = false;
			m_bLBUB  = true;
			m_bNewton = false;
			break;
		case 2:
			m_bCLBUB = false;
			m_bLBUB  = false;
			m_bNewton = true;
			break;
		default:
			;
	}//end switch	

}


/* set up scale */
void GLWidget::SetScale(int iScale)
{
	
	m_TransZ = -(GLfloat) iScale;
	updateGL();

}


void GLWidget::SetDrawBoxes()
{
	m_bDrawBoxes = !m_bDrawBoxes;
	updateGL();
}
void GLWidget::SetDrawPolyLine()
{
	m_bDrawPolyLine = !m_bDrawPolyLine;
	updateGL();
}

void GLWidget::ToggleDistance()
{
	m_bDistance = !m_bDistance;
	updateGL();
}

void GLWidget::SetPointSize(int iSize)
{

}

void GLWidget::ToggleSearchTree()
{
	m_bSearchTree = !m_bSearchTree;
	updateGL();
}

void GLWidget::ToggleCircles()
{

	m_bCircles = !m_bCircles;
	updateGL();

}

void GLWidget::ToggleGrid()
{
	m_bDistanceField = !m_bDistanceField;
	updateGL();
}


void GLWidget::Rotate()
{


}//end rotate

void GLWidget::ToggleTranslation()
{
	m_bTranslation = !m_bTranslation;
}

void GLWidget::Translate()
{

}//end translate

void GLWidget::ToggleAnimation()
{
	m_bTimer = !m_bTimer;
	if(m_bTimer)
		animationTimer.start(25);
	else
		animationTimer.stop();
}//end toggleAnimation

void GLWidget::CollisionCheck(CApproxCurve *pCurve, int n)
{
	VECTOR2   vCog = pCurve->GetCenter();
	Real dRad    = pCurve->Radius();
	
	if( fabs(boundingBox.m_Vertices[1].x - vCog.x) <  dRad + 1e-6 )
	{
		dVelX[n] = -dVelX[n];
	}
	if( fabs(boundingBox.m_Vertices[0].x - vCog.x) < dRad + 1e-6 )
	{
		dVelX[n] = -dVelX[n];
	}
	if( fabs(boundingBox.m_Vertices[2].y - vCog.y) < dRad + 1e-6 )
	{
		dVelY[n] = -dVelY[n];
	}
	if( fabs(boundingBox.m_Vertices[0].y - vCog.y) < dRad + 1e-6)
	{
		dVelY[n] = -dVelY[n];
	}

}//end collisionCheck

void GLWidget::RebuildTree()
{

	if(m_nNumObjects <= 1)return;

	if(AABBTree.IsInitialized())
	{
		AABBTree.DestroyTree();
	}//end if

	//Get the built time
#ifdef WIN32
	LONGLONG frequency, lastTime, currentTime;
 
 	double timeElapsed, timeScale, tSum;
 	tSum = 0;
 
 	QueryPerformanceFrequency((LARGE_INTEGER*) &frequency);
 
 	timeScale = 1.0/frequency;

	QueryPerformanceCounter((LARGE_INTEGER*) &lastTime);
#else
	QTime timer;
	timer.start();
#endif
	
	PairIters pair(ApproxCurves.begin(), ApproxCurves.end());
	AABBTree.BuildTree(pair, m_nBuildMeth);

#ifdef WIN32
 	QueryPerformanceCounter((LARGE_INTEGER*) &currentTime);
 	timeElapsed = (currentTime - lastTime) * timeScale;
	cout<<"Built time: "<<timeElapsed<<endl;
#else
	cout<<"Built time: "<<timer.elapsed()<<endl;
#endif

	updateGL();
	
}//end RebuildTree

void GLWidget::LogPoint()
{
	
}//end logPoint

void GLWidget::PrintLUB()
{

}//end printLUB


void GLWidget::PerformTests()
{

}//end PerformTest

void GLWidget::PrintGrid()
{


}//end printGrid


// frees memory used by curves and data structures
void GLWidget::CleanUp(void)
{

	//delete nurbs curves
	if(!curve.empty())
	{
		for(int i = 0; i < m_nNumObjects; i++)
		{
			delete curve[i];
		}

		curve.clear();

	}//end if

	//delete nurbs curves
	if(!this->nxnurbs.empty())
	{
		for(int i = 0; i < m_nNumObjects; i++)
		{
			delete nxnurbs[i];
		}

		nxnurbs.clear();

	}//end if

	if(!ApproxCurves.empty())
	{
		for(int i = 0; i < m_nNumObjects; i++)
			delete ApproxCurves[i];

		ApproxCurves.clear();

	}//end if

	this->m_nNumObjects = 0;

	if(AABBTree.IsInitialized())
	{
		AABBTree.DestroyTree();
	}//end if

}//end cleanup



void GLWidget::MultithreadCoherency(void)
{

	static int times = 0;
	
	DistOps op;
	double dLUB = 1.7E+308;

	int i;
	
	VECTOR2 lastPoint(0,0);
	
	QTime timer;
	
	InitResources(&g_Grid, &AABBTree, ApproxCurves, nxnurbs);
#ifndef WIN32	
/*	pthread_t Threads[MAX_THREADS];
	void* ret[MAX_THREADS];

	pthread_create(&Threads[0], NULL, threadProc1, NULL);
	pthread_create(&Threads[1], NULL, threadProc2, NULL);

	timer.start();

	for( i=0; i < MAX_THREADS; i++)
		pthread_join(Threads[i], &ret[i]);
	
	cout<<"Iteration "<<times<<" Search time: "<<timer.elapsed()<<endl;	*/
#endif		
	DeInitResources();
	times++;


}

void GLWidget::TransformCurves(void)
{

	static int transform = 0;

	int i,j;


	CCollDetect co;
	CollList collisions;
	CollList::iterator lIter;
	
	co.CollAllPairs(ApproxCurves, collisions);
	//co.CollAllBVHs(ApproxCurves, collisions);
	
	bool *bColls = new bool[ApproxCurves.size()];
	for(j = 0; j <(int) ApproxCurves.size(); j++)
		bColls[j] = false;
	
	//compute Collision Response
	for(lIter = collisions.begin(); lIter != collisions.end(); lIter++)
	{
		CollPair cPair = *lIter;
		//if the collisions is not yet handled
		if(!bColls[cPair.first])
		{

//			this->CollisionResponse(cPair);

			dVelX[cPair.first] = -dVelX[cPair.first];
			dVelY[cPair.first] = -dVelY[cPair.first];
			bColls[cPair.first] = true;
			
			dVelX[cPair.second] = -dVelX[cPair.second];
			dVelY[cPair.second] = -dVelY[cPair.second];
			/*bColls[cPair.second] = true;		*/	
			
		}//end if
		
	}//end for

	for(int n = 0; n < (int)ApproxCurves.size(); n++)
	{
		CApproxCurve *pCurve = ApproxCurves[n];
		VECTOR2   vCog = pCurve->GetCenter();
		Real dRad    = pCurve->Radius();
		
		if( fabs(boundingBox.m_Vertices[1].x - vCog.x) <  dRad + 1e-6 )
		{
			if(!bColls[n])
				dVelX[n] = -dVelX[n];
		}
		if( fabs(boundingBox.m_Vertices[0].x - vCog.x) < dRad + 1e-6 )
		{
			if(!bColls[n])
				dVelX[n] = -dVelX[n];
		}
		if( fabs(boundingBox.m_Vertices[2].y - vCog.y) < dRad + 1e-6 )
		{
			if(!bColls[n])
				dVelY[n] = -dVelY[n];
		}
		if( fabs(boundingBox.m_Vertices[0].y - vCog.y) < dRad + 1e-6)
		{
			if(!bColls[n])
				dVelY[n] = -dVelY[n];
		}	
	}//end for

	MATRIX2X2 rotMat((double)cos(0.0), (double)-sin(0.0), (double)sin(0.0), (double)cos(0.0) );
	
	for(i = 0; i < m_nNumObjects; i++)
	{
		VECTOR2 trans(dVelX[i], dVelY[i]);
		ApproxCurves[i]->TransformCurve(rotMat,trans );
	}//end for

	if(AABBTree.IsInitialized())
		AABBTree.updateBBoxes(AABBTree.GetRoot());

	delete[] bColls;
	updateGL();
	
	cout<<"transformation step: "<<transform<<endl;
	transform++;

}//TransformCurves

void GLWidget::InitializeVelocities(void)
{
	for(int i = 0; i < MAX_PARTICLES; i++)
	{

		int x1 = rand()%100;
		int x2 = rand()%100;
		int sign = rand()%2;
		double dx1 = 1e-4 * double(x1);
		double dx2 = 1e-4 * double(x2);
		dx1 = (sign == 0) ? -dx1 : dx1;
		sign = rand()%2;
		dx2 = (sign == 0) ? -dx2 : dx2;
		dVelX[i] = dx1;
		dVelY[i] = dx2;
	}//end for	
}

void GLWidget::PerturbMesh(void)
{


}


/*!
    \fn GLWidget::AnimateStep()
 */
void GLWidget::AnimateStep()
{
	

	TransformCurves();

	//MATRIX2X2 rotMat((double)cos(0.0), (double)-sin(0.0), (double)sin(0.0), (double)cos(0.0));

	//for(int i = 0; i < ApproxCurves.size(); i++)
	//{
	//	VECTOR2 trans(dVelX[i], dVelY[i]);
	//	ApproxCurves[i]->TransformCurve(rotMat,trans );
	//	CollisionCheck(ApproxCurves[i],i);	

	//}//end for

	//
	//
 //	if(myTree)
 //		myTree->updateBBoxes(myTree->getRoot());
	//
	////g_Counter = (g_Counter + 1) % ApproxCurves.size();
	////
	////pACurve = ApproxCurves[g_Counter];

	//updateGL();
	
}//end Animate Step


/*!
    \fn GLWidget::ToggleCorrection()
 */
void GLWidget::ToggleCorrection()
{
	m_NewtonCorrection = !m_NewtonCorrection;
}


/*!
    \fn GLWidget::MultithreadNewton()
 */
void GLWidget::MultithreadNewton()
{
	DistOps op;
	double dLUB = 1.7E+308;

	int i;
	
	VECTOR2 lastPoint(0,0);
	
	QTime timer;
	
	InitResources(&g_Grid, &AABBTree, ApproxCurves, nxnurbs);
#ifndef WIN32	
/*	pthread_t Threads[MAX_THREADS];
	void* ret[MAX_THREADS];

	pthread_create(&Threads[0], NULL, threadProcNewton1, NULL);
	pthread_create(&Threads[1], NULL, threadProcNewton2, NULL);

	timer.start();

	for( i=0; i < MAX_THREADS; i++)
		pthread_join(Threads[i], &ret[i]);
	
	cout<<"Search time: "<<timer.elapsed()<<endl;	*/
#endif		
	DeInitResources();
}


/*!
    \fn GLWidget::MakeTextures()
 */
void GLWidget::MakeTextures()
{
#ifndef WIN32
	QString file("/home/raphael/projects/svn/inshape2dmoi/gradientGMV.bmp");
#else
	QString file("../gradientGMV.bmp");
	//QString file("../woody.bmp");
#endif
	m_Images[0] = QPixmap(file);
	g_Textures[0] = bindTexture(m_Images[0], GL_TEXTURE_2D);

	QImage image = m_Images[0].toImage();

}


/*!
    \fn GLWidget::CalcMaxDistance()
 */
void GLWidget::CalcMaxDistance()
{
	
	g_Grid.MakeTCoords();

}//end CalcMaxDist


/*!
    \fn GLWidget::DeBugLUB()
 */
void GLWidget::DeBugLUB()
{
	m_Meth = DEBUG;
}


/*!
    \fn GLWidget::CLUBToggled()
 */
void GLWidget::CLUBToggled()
{
	m_Meth = CLBUB;
}


/*!
    \fn GLWidget::LUBToggled()
 */
void GLWidget::LUBToggled()
{
	m_Meth = LBUB;
}


/*!
    \fn GLWidget::DecPointIndex()
 */
void GLWidget::DecPointIndex()
{
    /// @todo implement me
}


/*!
    \fn GLWidget::IncPointIndex()
 */
void GLWidget::IncPointIndex()
{
	if(m_CountJ < GRID_POINTS_Y-1)
		this->m_CountJ++;
	else
	{
		m_CountJ = 0;
		if(m_CountI < GRID_POINTS_X-1)
		{
			m_CountI++;
		}//end if
		else
		{
			m_CountI = 0;
		}//End else
	}//end else
   updateGL();
}//end IncPointIndex


/*!
    \fn GLWidget::SetDepth()
 */
void GLWidget::SetDepth(int iDepth)
{
	m_iDepth = iDepth;
	updateGL();
}//end SetDepth


/*!
    \fn GLWidget::MergeVertices()
 */
void GLWidget::MergeVertices()
{
	if(this->m_nNumObjects < 1)
		return;
	
	
	//ApproxCurves[0]->DetermineOptimalEdge(info);
	
	updateGL();
}//end MergeVertices

void GLWidget::CollisionResponse(CollPair &Pair)
{

	double E = 2.00;

	CApproxCurve *pCurve1 = ApproxCurves[Pair.first];
	CApproxCurve *pCurve2 = ApproxCurves[Pair.second];

	//compute vector n from center to center
	VECTOR2 vN = VECTOR2::createVector(pCurve1->GetCenter(), pCurve2->GetCenter());

	vN.Normalize();

	VECTOR2 vT = VECTOR2::CreateNormal(vN);

	VECTOR2 Vel1(dVelX[Pair.first], dVelY[Pair.first]);

	VECTOR2 Vel2(dVelX[Pair.second], dVelY[Pair.second]);

	double vait = Vel1 * vT; 
	double vain = Vel1 * vN; 

	double vbit = Vel2 * vT; 
	double vbin = Vel2 * vN;

	double ma = 20.0;
	double mb = 15.0;

	double vafn = (mb * vbin * (E+1) + vain*(ma - E*mb))/(ma + mb);

	double vbfn = (mb * vain * (E+1) + vbin*(ma - E*mb))/(ma + mb);

	double vaft = vait;

	double vbft = vbit;

	VECTOR2 vAft(vafn, vaft);
	VECTOR2 vBft(vbfn, vbft);

	double xfa = vAft * VECTOR2(vN.x,vT.x);

	double yfa = vAft * VECTOR2(vN.y,vT.y);

	double xfb = vBft * VECTOR2(vN.x,vT.x);

	double yfb = vBft * VECTOR2(vN.y,vT.y);

	dVelX[Pair.first]=xfa;
	dVelY[Pair.first]=yfa;

	dVelX[Pair.second]=xfb;
	dVelY[Pair.second]=yfb;

}//end CollisionResponse

void GLWidget::BruteToggle()
{
	m_Meth = BRUTE;
}

void GLWidget::BruteSearchSingle(void)
{

	DistOps op;

	VECTOR2 lastPoint(0,0);

	QTime timer;
	timer.start();

	for(int i = 0; i < GRID_POINTS_X; i++)
	{
				
		for(int j = 0; j < GRID_POINTS_Y; j++)
		{
			//SearchResult result = op.bruteDistance(pACurve, g_Grid.Point(i,j));
			SearchResult result = op.BruteDistance_SIMD(pACurve, g_Grid.Point(i,j));
			//SearchResult result = op.LUBSearch(pACurve, g_Grid.Point(i,j), lBFS, nBest);
			g_Grid.SetDistance(i,j,result.res);
		}//end for j
	}//end for i

	cout<<" CPU Time Brute Force search: "<<timer.elapsed()<<endl;	
}

void GLWidget::SetAABBTreeDepth(int iDepth)
{
	m_iSTDepth = iDepth;
	updateGL();
}

// sets build method to CAABBTree::SPATIAL_MEDIAN
void GLWidget::SetSpatialMedian()
{
	m_nBuildMeth = CAABBTree::SPATIAL_MEDIAN;
}

// Sets the build method to CAABBTree::MEDIAN_CUT
void GLWidget::SetObjectMedian()
{
	m_nBuildMeth = CAABBTree::MEDIAN_CUT;
}

// Sets the built method to CAABBTree::MIN_AREA
void GLWidget::SetMinArea()
{
	m_nBuildMeth = CAABBTree::MIN_AREA;
}

// Sets object build method to simple
void GLWidget::SetObjBVHSimple()
{
	m_iObjBVH = GLWidget::SIMPLE;
}

// Sets object build method to mp heuristic
void GLWidget::SetObjBVHMonotonic()
{
	m_iObjBVH = GLWidget::MPHEURISTIC;
}

void GLWidget::ToggleLighting()
{

	if(m_bLighting)
	{
		//enable lighting
		glDisable(GL_LIGHT0);		
		glDisable(GL_LIGHTING);		
		glDisable(GL_COLOR_MATERIAL);
	}//end if
	else
	{
		//enable lighting
		glEnable(GL_LIGHT0);		
		glEnable(GL_LIGHTING);		
		glEnable(GL_COLOR_MATERIAL);
	}//end else

	m_bLighting = !m_bLighting;
	updateGL();

}

void GLWidget::SetZPlane(int iPlane)
{
	m_iZPlane = iPlane;
	updateGL();
}//end SetZPlane


void GLWidget::GetFictitiousBoundaryValues(void)
{

	QTime timer;

	timer.start();

	vector<CApproxCurve*>::iterator vIter;

// 	for(vIter = ApproxCurves.begin(); vIter != ApproxCurves.end(); vIter++)
// 	{
// 
// 		CApproxCurve *pCurve = *vIter;
// 
// 		int iMid = ALLPOINTS/2;
// 
// 		CFBValuesThread thread1(0, iMid, pCurve, &g_Grid);
// 
// 		CFBValuesThread thread2(iMid, ALLPOINTS, pCurve, &g_Grid);
// 
// 		thread1.begin();
// 		thread2.begin();
// 
// 		thread1.end();
// 		thread2.end();
// 
// 	}//end for


	cout<<"CPU Time Fictitious Boundary : "<<timer.elapsed()<<endl;

}//end GetFictitiousBoundaryValues

void GLWidget::GetFictitiousBoundaryValuesMultiple(void)
{

	QTime timer;

	timer.start();

	vector<CApproxCurve*>::iterator vIter;

	CFeat2Ops op;

	for(int i = 0; i < ALLPOINTS; i++)
	{
		for(vIter = ApproxCurves.begin(); vIter != ApproxCurves.end(); vIter++)
		{

			CApproxCurve *pCurve = *vIter;
			if(op.IsFictitiousBoundary(pCurve, g_Grid(i)))
			{
				Real res = (Real)7.0;
				g_Grid.SetDistance(i, res);
				break;
			}//end if
			else
			{
				Real res = (Real)10.0;
				g_Grid.SetDistance(i, res);
			}//end else

		}//end for
	}//end for i

	cout<<"CPU Time Fictitious Boundary : "<<timer.elapsed()<<endl;

}//end GetFictitiousBoundaryValuesMultiple
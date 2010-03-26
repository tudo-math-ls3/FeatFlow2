/***************************************************************************
 *   Copyright (C) 2006 by Raphael M�nster   *
 *   raphael@cortez   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#if !defined(_CVISUALDEBUG_H)
#define _CVISUALDEBUG_H
#include "nurbs.h"
#include "boxes.h"
//#include <GL/gl.h>
#include "GL/glut.h"
#include "aabbtree.h"
#include "approxcurve.h"
#include "grid.h"
#include "grid3d.h"
#include "aabb3.h"
#include "cbvhnode.h"


void LoadGLTextures();
class CModel3D;


class CVisualDebug
{
public:
	
	int LOD;
	float m_fPointSize;

	CVisualDebug(int lod = 150): LOD(lod){ m_fPointSize = 2.0; };
	~CVisualDebug(){};
	void DrawHull();
	void DrawControlPolyline(CBSpline* pCurve);
	void DrawCornerCutting(CBSpline* pCurve);
	void DrawCircle(VECTOR2 vCenter, float dRad, float r, float g, float b);
	void DrawCurve(CBSpline* pCurve, float r, float g, float b);
	VECTOR2 DrawCurvePoint(CBSpline* pCurve, float u, float r, float g, float b);
	void DrawPoint(const VECTOR2 &vPoint, float size, float r, float g, float b);
	void DrawShortestLine();
	void DrawONXCurve(const ON_NurbsCurve &curve, float r, float g, float b);
	void DrawSBBNode(CSBBNode *tNode, float r, float g, float b);
	void DrawSBBTree(CSBBTree *pTree);
	void DrawAABB(CAABB &pBox,float r, float g, float b);
	bool DrawSearchBox(float d, const VECTOR2& p);
	void DrawBox(CSBB* pBox, float r, float g, float b);
	void DrawOBBox(COBB* pBox, float r, float g, float b);
	void DrawDerivCurve(CBSpline* pCurve);
	void DrawEdge(const VECTOR2 &vV1 ,const VECTOR2 &vV2);
	void DrawCog(CBSpline* pCurve);
	void DrawVector2f(const VECTOR2 &vV1);
	void DrawAABBTreeNode(CAABBNode *aaBBNode, float r, float g, float b, int d, int iTargetDepth);
	void DrawAABBTree(CAABBTree *pTree, int Depth);
	void setPSize(float fPSize);
	void DrawApproxCurve(CApproxCurve *pCurve, float r, float g, float b);
	void DrawGrid(VECTOR2 Grid[][GRID_POINTS_Y], int GridX, int GridY);
	void DrawtGrid(CGrid &Grid, int GridX, int GridY);
	void DrawBoundingCircle(VECTOR2 vMid, double dRad, float r, float g, float b);
	void DrawBoundingCircle(CSBBNode *tNode, float r, float g, float b);
	void DrawCircleTree(CBCNode *root, int iTrgDepth);
	void DrawPoints(vector<VECTOR2> &points, float fSize, float r, float g, float b);
	float* GetColor(float pressure);
	void DrawBCNode(CBCNode *tNode, float r, float g, float b,int depth, int iTrgDepth);
	
public:
	void DrawEdges(CGrid &grid);
	void DrawBCLeaf(CBCNode *tNode, float r, float g, float b,int depth);
    void DrawMedialAxis(CApproxCurve* pCurve);
    void DrawVoronoiEdges(CApproxCurve *pCurve);
    void DrawVoronoiCircles(CApproxCurve *pCurve);
	void DrawAxisCircle(CMedAxisVertex *pVertex, float r, float g, float b, DebugSEC &info);
//    void DrawAdjList(UndirectedGraph &graph);
	void DrawBVHTree(AABBTree3f &pTree);
	void DrawBVHNode(AABBNode3f *pNode);
	
public:
	void DrawModel3D(CModel3D& pModel);
public:
	void DrawGrid3D(GRID3D& grid);
public:
	void DrawCutPlaneZ(GRID3D& grid,int iSlice, double maxDist);
public:
	int DrawMeshPoints(GRID3D& grid);
public:
	void DrawAABB3(const CAABB3f &bxBox);
public:
	void DrawNormals(CModel3D& model);
};

#endif  //_CVISUALDEBUG_H

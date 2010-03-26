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

#include <visualdebug.h>
#include <obb2.h>
#include <sbbtree.h>
#include <model3d.h>


void CVisualDebug::DrawHull()
{

}



void CVisualDebug::DrawBoundingCircle(VECTOR2 vMid, double dRad, float r, float g, float b)
{

	double twoPi = 2 * 3.1415926535;

	double dT = twoPi / 100.0;
	double T  = 0.0;

	glColor3f(r,g,b);

	glBegin(GL_LINE_STRIP);
	for(int i = 0; i<=100.0; i++) 
	{
		VECTOR2 p1;
		p1.x = vMid.x + dRad * cos(T);
		p1.y = vMid.y + dRad * sin(T);
		glVertex3f(p1.x, p1.y, 0);
		T+=dT;
	}
	glEnd();
	
	glPointSize(5);
	
	glBegin(GL_POINTS);
		glVertex3f(vMid.x, vMid.y, 0.0);
	glEnd();
	
	glPointSize(1);

}//end DrawBoundingCircle

void CVisualDebug::DrawGrid(VECTOR2 Grid[][GRID_POINTS_Y], int GridX, int GridY)
{

	glPointSize(1.0);
	glColor3f(1,0.7,0.3);
	glBegin(GL_POINTS);
	int i,j;
	for(i = 0; i < GridY; i++)
		for(j = 0; j < GridX; j++)
		{
			glVertex3f(Grid[j][i].x, Grid[j][i].y, 0.0);
		}//end for j
	glEnd();

	for (i = 0; i < GridX - 1; i++)
	{
		glBegin(GL_TRIANGLE_STRIP);
		for (j = 0; j < GridY; j++)
		{
			glVertex3f(Grid[i][j].x, Grid[i][j].y, 0);
			glVertex3f(Grid[i+1][j].x, Grid[i+1][j].y, 0);
		}
		glEnd();
	}

}//end DrawGrid

void CVisualDebug::DrawPoints(vector<VECTOR2> &points, float fSize, float r, float g, float b)
{


	glPointSize(fSize);
	glColor3f(r,g,b);
	glBegin(GL_POINTS);
	int i;
	for(i = 0; i < (int)points.size(); i++)
	{
		glVertex3f(points[i].x, points[i].y, 0);
	}//end for
	glEnd();

}

void CVisualDebug::DrawtGrid(CGrid &Grid, int GridX, int GridY)
{
	for(int i = 0; i < GridX-1; i++)
	{
		glDrawElements(GL_TRIANGLE_STRIP, 2*GridY, GL_UNSIGNED_INT, &Grid.m_pIndices[2*i*GridX]);
	}//end for
}

void CVisualDebug::DrawControlPolyline(CBSpline* pCurve)
{
	glPointSize(m_fPointSize);

	// Draw curve hull
	glColor3f(0.8,0.8,0.8);
	
	glBegin(GL_LINE_STRIP);
	for(int i=0;i<pCurve->getNumPoints();++i) {
		glVertex3f( pCurve->GetCP(i).x,
					pCurve->GetCP(i).y, 0 );
	}
	glEnd();

	glColor3f(0,1,0);
	glBegin(GL_POINTS);
	for(int i=0;i!=pCurve->getNumPoints();++i) {
		glVertex3f( pCurve->GetCP(i).x, pCurve->GetCP(i).y, 0 );
	}
	
	
	glEnd();

}


void CVisualDebug::DrawCornerCutting(CBSpline* pCurve)
{

}

void CVisualDebug::DrawAABB(CAABB &pBox,float r, float g, float b)
{
	
	//glColor3f(0.7,0.7,1.0);
	glColor3f(r,g,b);
	glBegin(GL_LINE_STRIP);
	glVertex3f(pBox.m_Vertices[0].x, pBox.m_Vertices[0].y, 0);
	glVertex3f(pBox.m_Vertices[1].x, pBox.m_Vertices[0].y, 0);
	glVertex3f(pBox.m_Vertices[1].x, pBox.m_Vertices[1].y, 0);
	glVertex3f(pBox.m_Vertices[0].x, pBox.m_Vertices[1].y, 0);
	glVertex3f(pBox.m_Vertices[0].x, pBox.m_Vertices[0].y, 0);
	glEnd();

	glColor3f(0,0,1);
	glBegin(GL_POINTS);
	glVertex3f(pBox.m_Vertices[1].x, pBox.m_Vertices[0].y, 0);
	glVertex3f(pBox.m_Vertices[0].x, pBox.m_Vertices[1].y, 0);
	glVertex3f(pBox.m_Vertices[1].x, pBox.m_Vertices[1].y, 0);
	glVertex3f(pBox.m_Vertices[0].x, pBox.m_Vertices[0].y, 0);
	glEnd();

}

void CVisualDebug::DrawAABBTreeNode(CAABBNode *aabbNode, float r, float g, float b, int d, int iTargetDepth)
{

	if(aabbNode != NULL)
	{

		glColor3f(0,0,1);
		DrawAABBTreeNode(aabbNode->m_Children[0],0,0,1, d+1, iTargetDepth);
		glColor3f(1,0,0);
		DrawAABBTreeNode(aabbNode->m_Children[1],1,0,0, d+1, iTargetDepth);
		
		if(d==iTargetDepth)
		{
			CAABB &box = aabbNode->GetBoundingBox();
			DrawAABB(box, r,g,b);
		}
	}//end if

}//end DrawAABBTreeNode

void CVisualDebug::DrawAABBTree(CAABBTree *pTree, int Depth)
{

	DrawAABBTreeNode(pTree->GetRoot(),1,1,0,0,Depth);

}//end DrawAABBTree

void CVisualDebug::DrawApproxCurve(CApproxCurve *pCurve, float r, float g, float b)
{

	VECTOR2 p1;
	VECTOR2 *vSamples = pCurve->GetSamples();
	glColor3f(r,g,b);
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(2, GL_FLOAT, 0, vSamples);
	glDrawArrays(GL_LINE_LOOP,0, pCurve->GetResolution()-1);

	glColor3f(1,0,0);
	glPointSize(m_fPointSize);
	glBegin(GL_POINTS);
		glVertex3f(pCurve->GetCenter().x, pCurve->GetCenter().y, 0);
	glEnd();

}//end DrawApproxCurve

void CVisualDebug::DrawONXCurve(const ON_NurbsCurve &curve, float r, float g, float b)
{

	glColor3f(r,g,b);

	ON_3dPoint p1;
	glColor3f(1,0.7,0.3);

	float dP =(float)( (curve.m_knot[curve.CVCount()] - curve.m_knot[curve.m_order-2]) / LOD);
	float P = curve.m_knot[curve.m_order-2];
	
	//Draw b-spline curve
	glBegin(GL_LINE_STRIP);
	for(int i = 0; i<=LOD; i++) 
	{

		

		if(P > curve.m_knot[curve.CVCount()])
			P=curve.m_knot[curve.CVCount()];
		curve.EvPoint(P, p1);
		
		glVertex3f(p1.x, p1.y, p1.z);
		P+=dP;
	}

	glEnd();


}//end DrawONXCurve

void CVisualDebug::DrawCurve(CBSpline* pCurve, float r, float g, float b)
{	
	VECTOR2 p1;
	glColor3f(1,0.7,0.3);
	float dP =(float)( (pCurve->KnotVector()[pCurve->getNumPoints()]-pCurve->KnotVector()[pCurve->getDeg()]) / LOD);
	float P = pCurve->KnotVector()[pCurve->getDeg()];
	//p1 = deBoor(T,Points, true);
	//Draw b-spline curve
	glBegin(GL_LINE_STRIP);
	for(int i = 0; i<=LOD; i++) 
	{
		if(P > pCurve->KnotVector()[pCurve->getNumPoints()])
			P=pCurve->KnotVector()[pCurve->getNumPoints()];
		p1 = pCurve->CoxDeBoor(P);
		//p1 = deBoor1(P);
		glVertex3f(p1.x, p1.y, 0);
		P+=dP;
	}

	glEnd();
	

}



VECTOR2 CVisualDebug::DrawCurvePoint(CBSpline* pCurve, float u, float r, float g, float b)
{
	VECTOR2 p1 = pCurve->CoxDeBoor(u);
    glPointSize(m_fPointSize);
	glColor3f(r,g,b);
	glBegin(GL_POINTS);
		glVertex3f(p1.x,p1.y,0);
	glEnd();
	return p1;
}


bool CVisualDebug::DrawSearchBox(float d, const VECTOR2& p)
{
	glColor3f(0,0,1);
	glBegin(GL_LINE_STRIP);
		glVertex3f(p.x-d, p.y+d, 0);
		glVertex3f(p.x+d, p.y+d, 0);
		glVertex3f(p.x+d, p.y-d, 0);
		glVertex3f(p.x-d, p.y-d, 0);
		glVertex3f(p.x-d, p.y+d, 0);
	glEnd();
	return true;
}

void CVisualDebug::DrawShortestLine()
{

}

void CVisualDebug::DrawOBBox(COBB* pBox, float r, float g, float b)
{
	glColor3f(r,g,b);
	glBegin(GL_LINE_STRIP);
	glVertex3f(pBox->m_Vertices[3].x, pBox->m_Vertices[3].y, 0);
	glVertex3f(pBox->m_Vertices[2].x, pBox->m_Vertices[2].y, 0);
	glVertex3f(pBox->m_Vertices[1].x, pBox->m_Vertices[1].y, 0);
	glVertex3f(pBox->m_Vertices[0].x, pBox->m_Vertices[0].y, 0);
	glVertex3f(pBox->m_Vertices[3].x, pBox->m_Vertices[3].y, 0);
	glEnd();

	glColor3f(r,g,b);
	glBegin(GL_POINTS);
	glVertex3f(pBox->m_Vertices[3].x, pBox->m_Vertices[3].y, 0);
	glVertex3f(pBox->m_Vertices[2].x, pBox->m_Vertices[2].y, 0);
	glVertex3f(pBox->m_Vertices[1].x, pBox->m_Vertices[1].y, 0);
	glVertex3f(pBox->m_Vertices[0].x, pBox->m_Vertices[0].y, 0);
	glEnd();


}

void CVisualDebug::DrawBox(CSBB* pBox, float r, float g, float b)
{
	glColor3f(r,g,b);
	glBegin(GL_LINE_STRIP);
	glVertex3f(pBox->m_Vertices[3].x, pBox->m_Vertices[3].y, 0);
	glVertex3f(pBox->m_Vertices[2].x, pBox->m_Vertices[2].y, 0);
	glVertex3f(pBox->m_Vertices[1].x, pBox->m_Vertices[1].y, 0);
	glVertex3f(pBox->m_Vertices[0].x, pBox->m_Vertices[0].y, 0);
	glVertex3f(pBox->m_Vertices[3].x, pBox->m_Vertices[3].y, 0);
	glEnd();

	glColor3f(0,0,1);
	glBegin(GL_POINTS);
	glVertex3f(pBox->m_Vertices[3].x, pBox->m_Vertices[3].y, 0);
	glVertex3f(pBox->m_Vertices[2].x, pBox->m_Vertices[2].y, 0);
	glVertex3f(pBox->m_Vertices[1].x, pBox->m_Vertices[1].y, 0);
	glVertex3f(pBox->m_Vertices[0].x, pBox->m_Vertices[0].y, 0);
	glEnd();

	glColor3f(1,1,0);
	glBegin(GL_POINTS);
	glVertex3f(pBox->m_Vertices[pBox->m_CurvePoints[0]].x, pBox->m_Vertices[pBox->m_CurvePoints[0]].y, 0);
	glVertex3f(pBox->m_Vertices[pBox->m_CurvePoints[1]].x, pBox->m_Vertices[pBox->m_CurvePoints[1]].y, 0);
	glEnd();
}

void CVisualDebug::DrawPoint(const VECTOR2 &vPoint, float size, float r, float g, float b)
{
	glPointSize(size);

	glColor3f(r,g,b);

	glBegin(GL_POINTS);
		glVertex3f(vPoint.x, vPoint.y, 0);
	glEnd();
}//end DrawPoint

void CVisualDebug::DrawDerivCurve(CBSpline* pCurve)
{
	VECTOR2* p1;
	glColor3f(1.0,0.7,0.3);
	double dP =(double)( (pCurve->KnotVector()[pCurve->getNumPoints()] - pCurve->KnotVector()[pCurve->getDeg()]) / LOD);
	double P = pCurve->KnotVector()[pCurve->getDeg()];
	//p1 = deBoor(T,Points, true);
	//Draw b-spline curve
	
	glBegin(GL_LINE_STRIP);
	for(int i = 0; i<=LOD; i++) 
	{
		if(P > pCurve->KnotVector()[pCurve->getNumPoints()])
			P=pCurve->KnotVector()[pCurve->getDeg()];
		p1 = pCurve->curveDerivs(P,pCurve->getDeg()-1);
		//p1 = deBoor1(P);
		glVertex3f(p1[1].x, p1[1].y, 0);
		P+=dP;
		delete[] p1;
	}

	glEnd();

}

void CVisualDebug::DrawCog(CBSpline* pCurve)
{
	glPointSize(m_fPointSize);
	glColor3f(1,0,0);
	glBegin(GL_POINTS);
	glVertex3f(pCurve->getCog().x, pCurve->getCog().y, 0);
	glEnd();
}

void CVisualDebug::DrawVector2f(const VECTOR2 &vV1)
{
	glColor3ub(255,120,0);
	glBegin(GL_LINES);
		glVertex3f(0,0,0);
		glVertex3f(vV1.x, vV1.y, 0);
	glEnd();
}

void CVisualDebug::DrawCircle(VECTOR2 vCenter, float dRad, float r, float g, float b)
{
	
	
	int i;

	float twoPi = 2 * 3.1415926535;

	float dT = twoPi / 100.0;
	float T  = 0.0;

	glColor3f(r,g,b);

	glPointSize(3);
	glBegin(GL_POINTS);
	glVertex3f(vCenter.x, vCenter.y, 0.0);
	glEnd();
	glPointSize(1);

	glBegin(GL_LINE_STRIP);
	for(i = 0; i<=100.0; i++) 
	{
		VECTOR2 p1;
		p1.x = vCenter.x + dRad * cos(T);
		p1.y = vCenter.y + dRad * sin(T);
		glVertex3f(p1.x, p1.y, 0);
		T+=dT;
	}
	glEnd();

	
}//end DrawBoundingCircle

void CVisualDebug::DrawEdge(const VECTOR2 &vV1 ,const VECTOR2 &vV2)
{
	glColor3ub(255,120,0);
	glBegin(GL_LINES);
		glVertex3f(vV1.x, vV1.y, 0);
		glVertex3f(vV2.x, vV2.y, 0);
	glEnd();
}

void CVisualDebug::setPSize(float fPSize)
{
	m_fPointSize = fPSize;
}



void CVisualDebug::DrawSBBTree(CSBBTree *pTree)
{

	DrawSBBNode(pTree->getRoot(),1,1,0);

}//end DrawSearchTree

void CVisualDebug::DrawCircleTree(CBCNode *root, int TrgDepth)
{
	DrawBCNode(root, 1, 1, 0,0, TrgDepth);
}//end DrawCircleTree

void CVisualDebug::DrawBCNode(CBCNode *tNode, float r, float g, float b,int depth, int TrgDepth)
{

	//if(!tNode->IsLeaf())
	if(depth < TrgDepth && !tNode->IsLeaf())
	{
		for(int j = 0; j < 2; j++)
		{
			DrawBCNode(tNode->m_Children[j],r,g,b,depth+1, TrgDepth);
		}
	}//end if
	else
	{

		glColor3f(r,g,b);
		VECTOR2 vCenter = tNode->GetCenter();

		double dRad = tNode->GetRad();

		double twoPi = 2 * 3.1415926535;

		double dT = twoPi / 100.0;
		double T  = 0.0;

		glColor3f(r,g,b);

		glBegin(GL_LINE_STRIP);
		for(int i = 0; i<=100.0; i++) 
		{
			VECTOR2 p1;
			p1.x = vCenter.x + dRad * cos(T);
			p1.y = vCenter.y + dRad * sin(T);
			glVertex3f(p1.x, p1.y, 0);
			T+=dT;
		}
		glEnd();

		glPointSize(3.0);
		glBegin(GL_POINTS);
			glVertex3f(vCenter.x, vCenter.y, 0.0);
		glEnd();
		glPointSize(1.0);
	}//end else

}//end DrawBCNode

void CVisualDebug::DrawBoundingCircle(CSBBNode *tNode, float r, float g, float b)
{
	VECTOR2 vCenter = tNode->m_Box->getCenter();

	int i;

	double dRad = -1.7E+308;


	for(i = 0; i < 4; i++)
	{
		VECTOR2 vVertex = tNode->m_Box->m_Vertices[i];

		vVertex  = VECTOR2::createVector(vCenter, vVertex);

		double d = vVertex.mag();

		if(dRad < d)
			dRad = d;

	}//end for

	double twoPi = 2 * 3.1415926535;

	double dT = twoPi / 100.0;
	double T  = 0.0;

	glColor3f(r,g,b);

	glBegin(GL_LINE_STRIP);
	for(int i = 0; i<=100.0; i++) 
	{
		VECTOR2 p1;
		p1.x = vCenter.x + dRad * cos(T);
		p1.y = vCenter.y + dRad * sin(T);
		glVertex3f(p1.x, p1.y, 0);
		T+=dT;
	}
	glEnd();

	
}//end DrawBoundingCircle

void CVisualDebug::DrawSBBNode(CSBBNode *tNode, float r, float g, float b)
{

	if(!tNode->IsLeaf())
	{
		for(int i = 0; i < 2; i++)
		{
			DrawSBBNode(tNode->m_Children[i],r,g,b);
		}
	}//end if
	else
	{

		glColor3f(r,g,b);
		
		glBegin(GL_LINE_STRIP);
		for(int i = 0; i <= 4; i++)
		{
			glVertex3f(tNode->m_Box->m_Vertices[i%4].x, tNode->m_Box->m_Vertices[i%4].y, 0);
		}//end for
		glEnd();

		glPointSize(m_fPointSize);
		glColor3f(0.6,0.6,0.2);
		glBegin(GL_POINTS);
		for(int i = 0; i < 4; i++)
		{
			glVertex3f(tNode->m_Box->m_Vertices[i].x, tNode->m_Box->m_Vertices[i].y, 0);
		}//end for
		glEnd();
	}//end else

}//end DrawSearchTreeNode

float* CVisualDebug::GetColor(float distance){
	float *returnVal = new float[3];
	float percent;
		
	if (distance < 0.0) {
		percent = distance;
		returnVal[0] = 68.0/255.0;
		returnVal[1] = 40/255.0;
		returnVal[2] = 84/255.0;
	}
	else if (distance==0.0)
	{
		returnVal[0] = 0;
		returnVal[1] = 0;
		returnVal[2] = 1;
	}
	else if (distance<=0.1)
	{
		percent = distance;
		returnVal[0] = 0.0;
		returnVal[1] = percent * 97.0f/255;
		returnVal[2] = 1 - percent * 2.0f/255.0f;
	}
	else if (distance<=0.2)
	{
		percent = (distance-0.1);
		returnVal[0] = 0.0;
		returnVal[1] = 105.0f/255 + percent* 113.0f/255;
		returnVal[2] = 1 - percent * 4.0f/255;
	}
	else if (distance<=0.5)
	{
		percent = (distance-0.25);
		returnVal[0] = 0.0;
		returnVal[1] = 210.0f/255 + percent * 45.0f/255;
		returnVal[2] = 251.0f/255 - percent * 43.0f/255;
	}
	else if (distance<=1.0)
	{
		percent = (distance-1.0);
		returnVal[0] = 0.0;
		returnVal[1] = 1;
		returnVal[2] = 208.0f/255 - percent * 96.9f/255; 
	}
	else if (distance<=1.5)
	{
		percent = (distance-1.0);
		returnVal[0] = 0;
		returnVal[1] = 1;
		returnVal[2] = 112.0f/255 - percent * 112.0f/255;
	}
	else if (distance<=1.75)
	{
		percent = (distance-1.5);
		returnVal[0] = percent * 107.0f/255;
		returnVal[1] = 1 - percent * 4.0f/255;
		returnVal[2] = 0;
	}
	else if (distance<=2.0)
	{
		percent = (distance-1.75);
		returnVal[0] = 107.0f/255 + percent * 102.0f/255;
		returnVal[1] = 251.0f/255 - percent * 10.0f/255;
		returnVal[2] = 0.0f;
	}
	else if (distance<=2.5)
	{
		percent = (distance-2.0);
		returnVal[0] = 209.0f/255 + percent * 35.0f/255;
		returnVal[1] = 241.0f/255 - percent * 39.0f/255;
		returnVal[2] = 0.0f;
	}
	else if (distance<=3.0)
	{
		percent = (distance-2.5);
		returnVal[0] = 244.0f/255 + percent * 6.0f/255;
		returnVal[1] = 202.0f/255 - percent * 102.0f/255;
		returnVal[2] = 0.0f;
	}
	else if (distance<=4.0)
	{
		percent = (distance-3.0);
		returnVal[0] = 250.0f/255 + percent * 5.0f/255;
		returnVal[1] = 100.0f/255 - percent * 100.0f/255;
		returnVal[2] = 0;
	}
	else 
	{
		returnVal[0] = 1;
		returnVal[1] = 0;
		returnVal[2] = 0;
	}
	return returnVal;
}


void CVisualDebug::DrawEdges(CGrid &grid)
{
	glLineWidth(2.0);
	for(int j = 0; j < grid.GetX(); j++)
	{
		glBegin(GL_LINES);


		glVertex3f(grid.Point(j,0).x, grid.Point(j,0).y,0.0);
		glVertex3f(grid.Point(j,grid.GetY() - 1).x, grid.Point(j,grid.GetY() - 1).y,0.0);
		

		glEnd();
	}

	for(int j = 0; j < grid.GetX(); j++)
	{
		glBegin(GL_LINES);
		
			glVertex3f(grid.Point(0,j).x, grid.Point(0,j).y,0.0);
			glVertex3f(grid.Point(grid.GetX() - 1, j).x, grid.Point(grid.GetX() - 1, j).y,0.0);

		glEnd();
	}

	glLineWidth(1.0);
	
}


/*!
    \fn CVisualDebug::DrawBCLeaf(CBCNode *tNode, float r, float g, float b,int depth)
 */
void CVisualDebug::DrawBCLeaf(CBCNode *tNode, float r, float g, float b,int depth)
{
	
	glColor3f(r,g,b);
	VECTOR2 vCenter = tNode->GetCenter();

	double dRad = tNode->GetRad();

	double twoPi = 2 * 3.1415926535;

	double dT = twoPi / 100.0;
	double T  = 0.0;

	glColor3f(r,g,b);

	glBegin(GL_LINE_STRIP);
	for(int i = 0; i<=100.0; i++) 
	{
		VECTOR2 p1;
		p1.x = vCenter.x + dRad * cos(T);
		p1.y = vCenter.y + dRad * sin(T);
		glVertex3f(p1.x, p1.y, 0);
		T+=dT;
	}
	glEnd();
	
}//end drawBCLeaves


/*!
    \fn CVisualDebug::DrawMedialAxis(CApproxCurve* pCurve)
 */
void CVisualDebug::DrawMedialAxis(CApproxCurve* pCurve)
{
   
	if(pCurve->m_MedAxis.empty())
		return;
	
	glColor3ub(255,0,150);
	
	//std::vector<VECTOR2> &medAxis = ;
	
	glPointSize(3.0);
	
	glBegin(GL_POINTS);
	for(int i = 0; i < (int)pCurve->m_MedAxis.size(); i++)
		glVertex3f(pCurve->m_MedAxis[i]->m_Circle.m_vCenter.x, pCurve->m_MedAxis[i]->m_Circle.m_vCenter.y, 0.0);	
	glEnd();
	
	glPointSize(2.0);
	
}//end DrawMedialAxis


/*!
    \fn CVisualDebug::DrawVoronoiEdges(CApproxCurve *pCurve)
 */
void CVisualDebug::DrawVoronoiEdges(CApproxCurve *pCurve)
{

	glColor3f(0,0,1);
	
	
	int iSize = (int)pCurve->m_vEdges.size();
	
	glBegin(GL_LINES);
	for(int i = 0; i < iSize; i++)
	{
		tEdge edge = pCurve->m_vEdges[i];
		tVertex vert1 = pCurve->m_vVoronoi[edge.e1]; 
		tVertex vert2 = pCurve->m_vVoronoi[edge.e2]; 
		glVertex3f(vert1.vVertex.x, vert1.vVertex.y, 0.0);
		glVertex3f(vert2.vVertex.x, vert2.vVertex.y, 0.0);
	}//end for
	glEnd();
	
}//end DrawVoronoiEdges


/*!
    \fn CVisualDebug::DrawVoronoiCircles(CApproxCurve *pCurve)
 */
void CVisualDebug::DrawVoronoiCircles(CApproxCurve *pCurve)
{
   
	int iSize = (int)pCurve->m_MedAxis.size();
	
	for(int i = 0; i < iSize; i++)
	{
		//DrawBoundingCircle(pCurve->m_MedAxis[i]->m_Circle.m_vCenter, pCurve->m_Radii[i], 0, 1, 0);
		DrawBoundingCircle(pCurve->m_MedAxis[i]->m_Circle.m_vCenter, pCurve->m_MedAxis[i]->m_Circle.m_dRad, 0, 1, 0);
	}//end for
	
}


/*!
    \fn CVisualDebug::DrawAxisCircle(CMedAxisVertex *pVertex, float r, float g, float b)
 */
void CVisualDebug::DrawAxisCircle(CMedAxisVertex *pVertex, float r, float g, float b, DebugSEC &info)
{
    DrawBoundingCircle(pVertex->m_Circle.m_vCenter, pVertex->m_Circle.GetRadius(), r,g,b);
	
	glColor3f(0,0,0);
	
	glPointSize(2.0);
		
	glBegin(GL_POINTS);
	for(int i = 0; i < (int)pVertex->m_Vertices.size(); i++)
	{
		VECTOR2 vec = *(pVertex->m_Vertices[i]);
		glVertex3f(vec.x, vec.y, 0.0);
	}//end for
	glEnd();
	
	if(info.outcome == 1 || info.outcome == 2)
	{
		glColor3f(0,0,1);
		glBegin(GL_LINES);
			glVertex3f(info.p.x, info.p.y, 0.0);
			glVertex3f(info.q.x, info.q.y, 0.0);
			
/*			glColor3f(1,1,0);
			glVertex3f(info.p.x, info.p.y, 0.0);
			glVertex3f(info.pq.x, info.pq.y, 0.0);*/
			
			
			glColor3f(1,0,0);
			glVertex3f(info.p.x, info.p.y, 0.0);
			glVertex3f(info.r.x, info.r.y, 0.0);
			
			glColor3f(0,1,0);
			glVertex3f(info.q.x, info.q.y, 0.0);
			glVertex3f(info.r.x, info.r.y, 0.0);
			
			glColor3f(1,1,0);
			glVertex3f(info.p1.m_vP0.x, info.p1.m_vP0.y, 0.0);
			glVertex3f(info.p1.m_vP0.x - 2*info.p1.m_vDir.x, info.p1.m_vP0.y - 2*info.p1.m_vDir.y, 0.0);
			
			glColor3f(0,0,0);
			glVertex3f(info.p2.m_vP0.x, info.p2.m_vP0.y, 0.0);
			glVertex3f(info.p2.m_vP0.x - 2*info.p2.m_vDir.x, info.p2.m_vP0.y - 2*info.p2.m_vDir.y, 0.0);
			
			
			
			
		glEnd();
	}//end if
	
	glPointSize(1.0);
	
}//end DrawAxisCircle


/*!
    \fn CVisualDebug::DrawAdjList(UndirectedGraph &graph)
 */
//void CVisualDebug::DrawAdjList(UndirectedGraph &graph)
//{
///*	std::pair<edge_iterator, edge_iterator>
//			edges(const adjacency_list& g);*/
//			
//	///typedef for propertymap
//	typedef property_map<UndirectedGraph, vertex_MedAxis_t>::type MAProperty;
//		
//	//get the property accessor
//	MAProperty prMAV = get(vertex_MedAxis, graph);
//	
//	//shorthand to the data type
//	typedef property_traits<MAProperty>::value_type MAType;
//	
//	//define edge iterator	
//	graph_traits<UndirectedGraph>::edge_iterator out, out_end;
//	
//	//two vertex descriptors for the edges
//	graph_traits<UndirectedGraph>::vertex_descriptor u, v;
//	
//	out = edges(graph).first;
//	out_end = edges(graph).second;
//	
//	glColor3f(0,0,1);
//	//iterate through the edges
//	glBegin(GL_LINES);
//	while(out!=out_end)
//	{
//		u = source(*out,graph);
//		v = target(*out,graph);
//		CMedAxisVertex *VertexProp0 = boost::get(prMAV, u);
//		CMedAxisVertex *VertexProp1 = boost::get(prMAV, v);		
//		glVertex3f(VertexProp0->m_Circle.GetCenter().x, VertexProp0->m_Circle.GetCenter().y, 0.0);
//		glVertex3f(VertexProp1->m_Circle.GetCenter().x, VertexProp1->m_Circle.GetCenter().y, 0.0);
//		out++;
//	}//end while
//	glEnd();
//	
//	glColor3ub(255,0,150);
//
//	glPointSize(3.0);
//	
//	//iterate through the points
//	boost::graph_traits<UndirectedGraph>::vertex_iterator vi;
//	glBegin(GL_POINTS);
//	for (vi = vertices(graph).first; vi != vertices(graph).second; ++vi)
//	{
//		v = *vi;
//		CMedAxisVertex *VertexProp0 = prMAV[v];
//		glVertex3f(VertexProp0->m_Circle.GetCenter().x, VertexProp0->m_Circle.GetCenter().y, 0.0);
//
//	}//end for
//	glEnd();
//	
//	for (vi = vertices(graph).first; vi != vertices(graph).second; ++vi)
//	{
//		v = *vi;
//		CMedAxisVertex *VertexProp0 = prMAV[v];
//		DrawBoundingCircle(VertexProp0->m_Circle.GetCenter(), VertexProp0->m_Circle.GetRadius(), 0.0, 1.0 ,0.0);
//
//	}//end for	
//	
//	glPointSize(1.0);	
//	
//}//end DrawAdjList

void CVisualDebug::DrawModel3D(CModel3D& pModel)
{

	const VEC3Array& pVertices = pModel.GetVertices();

	const FaceArray& pFaces    = pModel.GetFaces();

	const Vec3Array& pNormals = pModel.GetVertexNormals();

	

	glBegin(GL_TRIANGLES);
	//glBegin(GL_LINE_STRIP);

	int i,j;
	
	for(i = 0; i <(int)pFaces.size(); i++)
	{

		for(j = 0; j < 3; j++)
		{
			

			int vertIndex = pFaces[i].VertexIndex[j];

			
			glNormal3f(pNormals[vertIndex].x, pNormals[vertIndex].y, pNormals[vertIndex].z);	
			glVertex3f(pVertices[vertIndex].x, pVertices[vertIndex].y, pVertices[vertIndex].z);
			

		}//end for

	}//end for

	glEnd();


}

void CVisualDebug::DrawGrid3D(GRID3D& grid)
{

	glColor3f(0,0,0);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	for(int i = 0; i < grid.GetSizeX(); i++)
	{
		for(int j = 0; j < grid.GetSizeY(); j++)
		{
			glVertex3f(grid(i,j,0).x, grid(i,j,0).y, grid(i,j,0).z);
			glVertex3f(grid(i,j,grid.GetSizeZ()-1).x, grid(i,j,grid.GetSizeZ()-1).y, grid(i,j,grid.GetSizeZ()-1).z);
		}
	}//end for i
	glEnd();

	glBegin(GL_LINES);
	for(int i = 0; i < grid.GetSizeY(); i++)
	{
		for(int j = 0; j < grid.GetSizeZ(); j++)
		{
			glVertex3f(grid(0,i,j).x, grid(0,i,j).y, grid(0,i,j).z);
			glVertex3f(grid(grid.GetSizeX()-1,i,j).x, grid(grid.GetSizeX()-1,i,j).y, grid(grid.GetSizeX()-1,i,j).z);
		}
	}//end for i
	glEnd();

	glBegin(GL_LINES);
	for(int i = 0; i < grid.GetSizeX(); i++)
	{
		for(int j = 0; j < grid.GetSizeZ(); j++)
		{
			glVertex3f(grid(i,0,j).x, grid(i,0,j).y, grid(i,0,j).z);
			glVertex3f(grid(i,grid.GetSizeY()-1,j).x, grid(i,grid.GetSizeY()-1,j).y, grid(i,grid.GetSizeY()-1,j).z);
		}//end for j
	}//end for i
	glEnd();
	glEnable(GL_LIGHTING);

}//end DrawGrid3D

void CVisualDebug::DrawCutPlaneZ(GRID3D& grid,int iSlice, double maxDist)
{

	int nSlice = grid.GetSizeX() * grid.GetSizeY();
	
	for (int i = 0; i < grid.GetSizeX() - 1; i++)
	{
		glBegin(GL_TRIANGLE_STRIP);
		for (int j = 0; j < grid.GetSizeY(); j++)
		{
			float sc = grid.m_Grid[i*grid.GetSizeX() + j + iSlice * nSlice].distance/maxDist;
			glTexCoord2f(sc, 0.0);
			glVertex3f(grid(i,j,iSlice).x, grid(i,j,iSlice).y, grid(i,j,iSlice).z);
			//sc = grid.Distance(i+1,j,iSlice)/maxDist;
			sc = grid.m_Grid[(i+1)*grid.GetSizeX() + j + iSlice * nSlice].distance/maxDist;
			glTexCoord2f(sc, 0.0);
			glVertex3f(grid(i+1,j,iSlice).x, grid(i+1,j,iSlice).y, grid(i,j,iSlice).z);
		}//end j
		glEnd();
	}//end i	
	
}

int CVisualDebug::DrawMeshPoints(GRID3D& grid)
{

	int nAll = grid.GetSizeX() * grid.GetSizeY() * grid.GetSizeZ();

	glColor3f(1,0,0);
	glPointSize(3);
	glBegin(GL_POINTS);
	for(int i = 0; i < nAll; i++)
	{
		if(grid.m_Grid[i].Inside)
		{
			glColor3f(1,0,0);
			glVertex3f(grid(i).x, grid(i).y, grid(i).z);
		}
		else
		{
			glColor3f(0,1,0);
			glVertex3f(grid(i).x, grid(i).y, grid(i).z);
		}

	}//End for
	glEnd();
	glPointSize(1);
	return 0;
}

void CVisualDebug::DrawAABB3(const CAABB3f& bxBox)
{

	glColor3f(1,0,0);

	glBegin(GL_LINES);
	//1
	glVertex3f(bxBox.GetFBL().x, bxBox.GetFBL().y, bxBox.GetFBL().z); 
	glVertex3f(bxBox.GetFTL().x, bxBox.GetFTL().y, bxBox.GetFTL().z); 
	//2
	glVertex3f(bxBox.GetFTL().x, bxBox.GetFTL().y, bxBox.GetFTL().z); 
	glVertex3f(bxBox.GetFTR().x, bxBox.GetFTR().y, bxBox.GetFTR().z);
	//3
	glVertex3f(bxBox.GetFTR().x, bxBox.GetFTR().y, bxBox.GetFTR().z);
	glVertex3f(bxBox.GetFBR().x, bxBox.GetFBR().y, bxBox.GetFBR().z);
	//4
	glVertex3f(bxBox.GetFBR().x, bxBox.GetFBR().y, bxBox.GetFBR().z);
	glVertex3f(bxBox.GetFBL().x, bxBox.GetFBL().y, bxBox.GetFBL().z); 
	//5
	glVertex3f(bxBox.GetFBR().x, bxBox.GetFBR().y, bxBox.GetFBR().z);
	glVertex3f(bxBox.GetBBR().x, bxBox.GetBBR().y, bxBox.GetBBR().z);
	//6
	glVertex3f(bxBox.GetBBR().x, bxBox.GetBBR().y, bxBox.GetBBR().z);
	glVertex3f(bxBox.GetBTR().x, bxBox.GetBTR().y, bxBox.GetBTR().z);
	//7
	glVertex3f(bxBox.GetBTR().x, bxBox.GetBTR().y, bxBox.GetBTR().z);
	glVertex3f(bxBox.GetFTR().x, bxBox.GetFTR().y, bxBox.GetFTR().z);
	//8
	glVertex3f(bxBox.GetBTR().x, bxBox.GetBTR().y, bxBox.GetBTR().z);
	glVertex3f(bxBox.GetBTL().x, bxBox.GetBTL().y, bxBox.GetBTL().z);
	//9
	glVertex3f(bxBox.GetBBR().x, bxBox.GetBBR().y, bxBox.GetBBR().z);
	glVertex3f(bxBox.GetBBL().x, bxBox.GetBBL().y, bxBox.GetBBL().z);
	//10
	glVertex3f(bxBox.GetBTL().x, bxBox.GetBTL().y, bxBox.GetBTL().z);
	glVertex3f(bxBox.GetBBL().x, bxBox.GetBBL().y, bxBox.GetBBL().z);
	//11
	glVertex3f(bxBox.GetBBL().x, bxBox.GetBBL().y, bxBox.GetBBL().z);
	glVertex3f(bxBox.GetFBL().x, bxBox.GetFBL().y, bxBox.GetFBL().z);
	//12
	glVertex3f(bxBox.GetBTL().x, bxBox.GetBTL().y, bxBox.GetBTL().z);
	glVertex3f(bxBox.GetFTL().x, bxBox.GetFTL().y, bxBox.GetFTL().z);

	glEnd();

}

void CVisualDebug::DrawBVHTree(AABBTree3f &pTree)
{
	DrawBVHNode(pTree.GetRoot());
}//end DrawBVHTree

void CVisualDebug::DrawBVHNode(AABBNode3f *pNode)
{
	if(!pNode->IsLeaf())
	{
		for(int j = 0; j < 2; j++)
		{
			DrawBVHNode(pNode->m_Children[j]);
		}
	}//end if
	else
	{
		DrawAABB3(pNode->GetBV());
	}//end else
}//end DrawBVHTree


void CVisualDebug::DrawNormals(CModel3D& model)
{
	const Vec3Array& vertexNormals = model.GetVertexNormals();
	const VEC3Array& vertices = model.GetVertices();

	for(int i = 0; i <vertices.Size();i++)
	{
		glBegin(GL_LINES);
		glVertex3f(vertices[i].x, vertices[i].y,vertices[i].z);
		glVertex3f(vertices[i].x + vertexNormals[i].x, vertices[i].y + vertexNormals[i].y,vertices[i].z + vertexNormals[i].z);
		glEnd();
	}//end for

}

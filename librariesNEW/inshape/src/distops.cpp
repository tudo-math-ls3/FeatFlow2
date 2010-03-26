/***************************************************************************
 *   Copyright (C) 2006 by Raphael M�nster   *
 *   raphael@Cortez   *
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


#if defined WIN32
#include <windows.h>
#include <mmsystem.h>
#endif




#include <iostream>
#include <math.h>
#include <queue>
//#include <qdatetime.h>
#include <limits>

#if !defined WIN32
#include <pthread.h>
#endif

#include "distops.h"
#include "binstree.h"
#include "aabbtree.h"
#include "approxcurve.h"
#include "rootfinder.h"
#include "matrix2x2.h"
#include "obb2.h"
#include "sbb.h"



t_ThreadResource resources[MAX_THREADS];


void InitResources(CGrid *grid, CAABBTree *pTree, CurveVec &ApproxCurves, ONXVec& rNurbs)
{
	int i,j,k;
	
	vector<CNurbs*> vec;

	for(i = 0; i < MAX_THREADS; i++)
	{
		resources[i].ApproxCurves = ApproxCurves;
		resources[i].NURBS        = rNurbs;
		resources[i].pTree        = new CAABBTree(pTree);
		resources[i].Grid         = grid;
	}//end for    /// @todo implement me
}

void DeInitResources()
{
	
	for(int i = 0; i < MAX_THREADS; i++)
	{
		resources[i].pTree->DestroyTree();
	}//end for	
	
}//DeInitResources


void rotatePoint(double dAngle, CBSpline *pCurve, VECTOR2 &vPoint)
{
	VECTOR2 oldCog = pCurve->getCog();
	VECTOR2 vTrans = oldCog * -1.0;
 	vPoint = vPoint + vTrans;
	
	MATRIX2X2 rotMat((double)cos(dAngle), (double)-sin(dAngle), (double)sin(dAngle), (double)cos(dAngle) );

	vPoint = rotMat * vPoint;

	vPoint = vPoint + oldCog;

}//end rotatePoint



SearchResult DistOps::bruteDistance(CApproxCurve *pCurve, const VECTOR2 &vQuery)
{

	SearchResult distance;
	double d = 1.7E+308;

	int i;
	
	/* initialize distance with default values */
	distance.res  = 1.7E+308;
	distance.data = -1;

	//Real *dValues = new Real[pCurve->GetResolution()];

	if(pCurve->m_bBox->inside(vQuery))
	{
		if(pCurve->IsInElement(vQuery))
		{
			distance.res   = -1.0;
			distance.data  = 0;
			distance.point = 0;
			return distance;
		}
	}


	VECTOR2 c;
	
	VECTOR2 *vSamples = pCurve->GetSamples();
	
	int nNumP = pCurve->GetResolution();

	for(int i = 0; i < nNumP; i++)
	{
		VECTOR2 b = vSamples[i];
		c			= VECTOR2::createVector(b,vQuery);
		d			= c.norm2();
		if(distance.res > d)
		{
			distance.res = d;
			distance.point = i;
		}//end if
	}//end for

	distance.res = sqrt(distance.res);

	return distance;

}//end bruteDistance



////////////////////////////////////////////////////////////////////////////////
//
//computes the closest point on the line from
//vA to vB relative to vPoint
//
//Input:
//	-the points vA and vB define a line segment
//	-a reference point vPoint
//Output:
//	-the closest point on the line segment from
//	vA to vB from vPoints point of view
//
////////////////////////////////////////////////////////////////////////////////

VECTOR2 DistOps::ClosestPointOnLine(VECTOR2 vA, VECTOR2 vB, VECTOR2 vPoint)
{
	/* create the vector from end point vA to vPoint */
	VECTOR2 vVector1 = vPoint - vA;

	/* create a normalized direction vector from end point vA to vB */
	VECTOR2 vVector2 = (vB - vA).GetNormalizedVector();
	
	/* calculate distance AB */
	VECTOR2 vAB = VECTOR2::createVector(vA, vB);
	
    double d = vAB.mag();

	/* the dot product is now: dot(vV2, vV1) = cos(angle(vV1,vV2)) * |vV1| * 1 
	which gives the distance from vA to the projection of vPoint onto line vAB */
	double t = vVector2 * vVector1;

	/* If t is less than or equal to 0, it must be closest to the end point vA */
    if (t <= 0) 
		return vA;

	/* If t is greater than or equal to the length of the line segment it must be closest to the end point vB */
    if (t >= d) 
		return vB;
 
	/* create a vector from vA to the closest point on the line */
    VECTOR2 vVector3 = vVector2 * t;

	/* Addition of vA to this vector yields the closest point on the line*/  
    VECTOR2 vClosestPoint = vA + vVector3;

	/* return the closest point on the line segment */
	return vClosestPoint;
}

double DistOps::distPointBox(VECTOR2 vPoint, CAABB *box)
{

	double d = 0.0;
	double minDist = 1.7E+308;
	VECTOR2 c;
	
	/* store comparison results in boolean */
	bool xGreaterRight;
	//bool yLessTop;
	//bool yGreaterBottom;
	bool yGreaterTop;
	bool yLessBottom;
	bool xLessLeft;

	/* partition 2d space */

	/* there are 9 cases */

	/* case 1 */
	if((xGreaterRight = (vPoint.x >= box->m_Vertices[1].x)) && (vPoint.y <= box->m_Vertices[1].y) && (vPoint.y >= box->m_Vertices[0].y))
	{
		VECTOR2 b(box->m_Vertices[1].x, vPoint.y);
		
		c = VECTOR2::createVector(vPoint, b);
		
		return c.mag();

	}
	/* case 2 */
	else if((yGreaterTop = (vPoint.y > box->m_Vertices[1].y)) && xGreaterRight )
	{
		c = VECTOR2::createVector(vPoint, box->m_Vertices[1]);
		return c.mag();
	}
	/* case 3 */
	else if((yLessBottom = (vPoint.y < box->m_Vertices[0].y)) && xGreaterRight )
	{
		c = VECTOR2::createVector(vPoint, VECTOR2(box->m_Vertices[1].x,box->m_Vertices[0].y));
		return c.mag();
	}
	/* case 4 */
	else if((xLessLeft = (vPoint.x <= box->m_Vertices[0].x)) && (vPoint.y <= box->m_Vertices[1].y) && (vPoint.y >= box->m_Vertices[0].y))
	{
		VECTOR2 b(box->m_Vertices[0].x, vPoint.y);
		c = VECTOR2::createVector(vPoint, b);
		return c.mag();
	}
	/* case 5 */
	else if(xLessLeft && yGreaterTop)
	{
		c = VECTOR2::createVector(vPoint, VECTOR2(box->m_Vertices[0].x,box->m_Vertices[1].y));
		return c.mag();
	}
	/* case 6 */
	else if(xLessLeft && yLessBottom)
	{
		c = VECTOR2::createVector(vPoint, box->m_Vertices[0]);
		
		
		return c.mag();
	}
	/* case 7 */
	else if(!xLessLeft && !xGreaterRight && yLessBottom)
	{
		VECTOR2 b(vPoint.x, box->m_Vertices[0].y);
		
		c = VECTOR2::createVector(vPoint, b);
		
		return c.mag();
	}
	/* case 8 */
	else if(!xLessLeft && !xGreaterRight && yGreaterTop)
	{
		VECTOR2 b(vPoint.x, box->m_Vertices[1].y);
		
		c = VECTOR2::createVector(vPoint, b);
		
		return c.mag();

	}
	else
	{

		VECTOR2 lines[4];
		lines[0] = box->m_Vertices[0];
		lines[1] = VECTOR2(box->m_Vertices[1].x, box->m_Vertices[0].y);
		lines[2] = box->m_Vertices[1];
		lines[3] = VECTOR2(box->m_Vertices[0].x, box->m_Vertices[1].y);
		
		for(int j = 1; j <=4; j++)
		{
			VECTOR2 b = ClosestPointOnLine(lines[(j-1+4)%4], lines[(j+4)%4], vPoint);
			c = VECTOR2::createVector(vPoint, b);
			d = c.mag();
			if(d < minDist)
			{
				minDist   = d;
				
			}
		}

		VECTOR2 vCenter = box->GetCenter();
		c = VECTOR2::createVector(vCenter, vPoint);
		d = c.mag();
		if(d < minDist)
			minDist = d;

	}
	return minDist;	
}


//computes the number of Intersection of a ray and a curve
/////////////////////////////////////////////////////////////////////////////////////////////////
//
// Input:
//	-pncurve: pointer to a CBSpline curve object
//	-pvPoint : a point in the xy-plane
//	the ray R is given implicitly:
//				- R starts at pvPoint.x
//				- R is parallel to the x-Axis
//				- R heads to the right(+inf )
//Output:
//	-the number of intersections between R and the curve
//
/////////////////////////////////////////////////////////////////////////////////////////////////
bool DistOps::numIntersections(CBSpline* pncurve, const VECTOR2& pvPoint)
{
	
	int i;
	int n = (int)pncurve->NumBoxes();
	int nI = 0;
	CRootFinder rFA;
	for(i = 0; i < n; i++)
	{
		CSBB* box = pncurve->GetBoxes()[i];
		//box intersection test with ray originating from pvPoint  to +inf (ray is parallel to x axis)
		if( ((pvPoint.y > box-> m_Vertices[3].y) || (pvPoint.y < box->m_Vertices[0].y) || (pvPoint.x > box->m_Vertices[1].x)) )
		{
			//printf("no hit\n");		
			continue;
		}
		//printf("hit\n");
			
		//test intersection with curve
		double dInsc;
		VECTOR2 vPoint;
		if(box->m_fP1 < box->m_fP2)
		{
			dInsc = rFA.bisectIntersec(pncurve, pvPoint.y, 1e-3, box->m_fP1, box->m_fP2,vPoint);
			if(!(vPoint.x <= pvPoint.x))
				nI++;
		}
		else
		{
			dInsc = rFA.IntersectWrap(pncurve, pvPoint.y, 1e-3, box->m_fP1, box->m_fP2);
			if(!(pncurve->CoxDeBoor(dInsc).x <= pvPoint.x))
				nI++;
		}
	}
	if((nI % 2) == 0)
		return false;
	else 
		return true;
}

bool isInElement(CBSpline* pncurve, const VECTOR2& pvPoint)
{
	
	int i;
	int n = (int)pncurve->NumBoxes();
	int nI = 0;
	CRootFinder rFA;
	for(i = 0; i < n; i++)
	{
		CSBB* box = pncurve->GetBoxes()[i];
		//box intersection test with ray originating from pvPoint  to +inf (ray is parallel to x axis)
		if( ((pvPoint.y > box-> m_Vertices[3].y) || (pvPoint.y < box->m_Vertices[0].y) || (pvPoint.x > box->m_Vertices[1].x)) )
		{
			//printf("no hit\n");		
			continue;
		}
		//printf("hit\n");
			
		//test intersection with curve
		double dInsc;
		if(box->m_fP1 < box->m_fP2)
		{
			//dInsc = rFA.bisectIntersec(pncurve, pvPoint.y, 1e-3, box->m_fP1, box->m_fP2);
			rFA.NewtonIntersec(pncurve, pvPoint.y, 1e-3, box->m_fP1, box->m_fP2, 5);
			if(!(pncurve->CoxDeBoor(dInsc).x <= pvPoint.x))
				nI++;
		}
		else
		{
			dInsc = rFA.IntersectWrap(pncurve, pvPoint.y, 1e-3, box->m_fP1, box->m_fP2);
			if(!(pncurve->CoxDeBoor(dInsc).x <= pvPoint.x))
				nI++;
		}
	}
	if((nI % 2) == 0)
		return false;
	else 
		return true;
}//end isInElement




int DistOps::checkSearchResults(double *dRes1, double *dRes2, const VECTOR2 &vPoint)
{
	if(fabs(dRes1[0] - dRes2[0]) > 0.05)
	{
		//log("Point: \n");
		//log("log.txt","%f\n%f",vPoint.x, vPoint.y);
		//log("difference.txt","Sol1 = %f, Sol2 = %f", dRes1[0], dRes2[0]);
		//log("difference.txt","Difference: %f", (fabs(dRes1[0] - dRes2[0])));
		return 0;
	}
	return 1;	
}

///==================================================
/// 	LUBSEARCH FOR SINGLE OBJECT
///==================================================

SearchResult DistOps::LUBSearch(CApproxCurve *Curve, const VECTOR2 &pvPoint, ListBCND &lCandidates, CBCNode *&pBest)
{

	double d = 1.7E+308;
	double mp = -1.0;
	/* initialize distance with default values */
	SearchResult distance;
	distance.res = 1.7E+308;
	distance.data = -1;
	
	
	ListBCND::iterator lIter;
	pBest = NULL;

	//if(Curve->m_bBox->inside(pvPoint))
	//{
	//	if(Curve->IsInElement(pvPoint))
	//	{
	//		distance.res   = -1.0;
	//		distance.data  = 0;
	//		distance.point = 0;
	//		return distance;
	//	}
	//}
	
	/* breadthFirstSearch */
	breadthFirstSearchDistance(Curve, pvPoint, lCandidates, pBest);

	/* compute closest point in this circle */
	int i;

	int nData = pBest->GetData();

	int nNumP = (int)Curve->m_Buckets[nData].size();

	for(i = 0; i < nNumP; i++)
	{
		VECTOR2 b = Curve->GetPoint(nData,i);
		b			= VECTOR2::createVector(b,pvPoint);
		d			= b.mag();
		if(distance.res > d)
		{
			distance.res   = d;
			distance.data  = nData;
			distance.point = i;
		}//end if
	}//end for

	for(lIter = lCandidates.begin(); lIter != lCandidates.end(); lIter++)
	{
		CBCNode *pBox = *lIter;
		if(pBox == pBest)
			continue;

		if(pBox->GetLowerBound(pvPoint) < distance.res)
		{
			
			nData = pBox->GetData();

			nNumP = (int)Curve->m_Buckets[nData].size();

			/* loop to compute minimum distance for this bounding circle */
			for(i = 0; i < nNumP; i++)
			{
				VECTOR2 b = Curve->GetPoint(nData,i);
				b			= VECTOR2::createVector(b,pvPoint);
				d			= b.mag();
				if(distance.res > d)
				{
					distance.res = d;
					distance.data = static_cast<double>(nData);
					distance.point = static_cast<double>(i);
				}//end if
			}//end for

		}//end if

	}//end for
	return distance;
}//end LUBSearch

void DistOps::breadthFirstSearchDistance(CApproxCurve* pCurve, const VECTOR2 &vQuery,
										 ListBCND &lBFS, CBCNode *&nBest)
{
	ListBCND::iterator lIter;
	double lowerBound  = 1.7E+308;
	double upperBound  = -1.7E+308;
	
	CBCTree &pTree = pCurve->m_BCTree;

	int s = pCurve->GetNumCircles();
	for(int i = 0; i < pTree.GetNumChildren(); i++)
	{
		lBFS.push_back(pTree.GetChild(i));
	}//end for
	
    nBest = NULL;

	int vSize = (int)lBFS.size();

	int nLeafCount = 0;
	int iters = 0;
	/* loop until there are only leaves in the list */
	while(vSize != nLeafCount)
	{
		
		nLeafCount = 0;
		int j = 0;
		lowerBound = 1.7E+308;
		double *dLowerBounds = new double[vSize];
		/* find best upper bound */
		for(lIter = lBFS.begin(); lIter != lBFS.end(); lIter++)
		{
			CBCNode *pNode = *lIter;
			dLowerBounds[j] = pNode->GetLowerBound(vQuery);
			if(lowerBound > dLowerBounds[j])
			{
				lowerBound = dLowerBounds[j];
				nBest = pNode;
			}//end if
			j++;
		}//end for
		
		/* get upper bound for best element */
		upperBound = nBest->GetUpperBound(vQuery);
		
		lIter = lBFS.begin();
		for(int i = 0; i < vSize; i++)
		{
			
			CBCNode *pNode = *lIter;
			if(pNode == nBest)
			{
				if(!pNode->IsLeaf())
				{
					lBFS.push_back(pNode->m_Children[0]);
					lBFS.push_back(pNode->m_Children[1]);
					lIter = lBFS.erase(lIter);
					continue;
				}
				else
				{
					nLeafCount++;
					lIter++;
					continue;
				}
			}//end if
			
			/* Can this subtree be pruned ? */
			if(upperBound > dLowerBounds[i])
			{
				if(!pNode->IsLeaf())
				{
					lBFS.push_back(pNode->m_Children[0]);
					lBFS.push_back(pNode->m_Children[1]);
					lIter = lBFS.erase(lIter);
				}//end if
				else
				{
					nLeafCount++;
					lIter++;
				}
			}//end if
			else
			{
				lIter = lBFS.erase(lIter);
			}//end else
			
		}//end for
		vSize = (int)lBFS.size();
		delete[] dLowerBounds;
		iters++;
	}//end while
	
}//end breadthFirstSearchDistance


///==================================================
/// 	LUBSEARCH FOR MULTIPLE OBJECTS
///==================================================

SearchResult DistOps::LUBSearch(ListAABBND &Candidates, const VECTOR2 &pvPoint, CApproxCurve *&nCurve, CurveVec &ApproxCurves)
{

	SearchResult distance;

	double d = std::numeric_limits<Real>::max();
	double mp = -1.0;
	/* initialize distance with default values */
	distance.res  = std::numeric_limits<Real>::max();
	distance.data = -1;


	/* the list will store the candidate leaves of the bspline curves segment box tree*/
	ListBCND lBFS;
	ListBCND::iterator lIter;
	CBCNode *pBest = NULL;
	
	/* pruning search, pBest will hold the nearest segment box */
	/* lower bound/upper bound search to find the candidate leaves */
	breadthFirstSearchDistance(Candidates, pvPoint, lBFS, pBest);

	/* get the pointer to the curve that the best leaf belongs to */
	nCurve = ApproxCurves[pBest->GetID()];


	if(nCurve->m_bBox->inside(pvPoint))
	{
		//cout<<"inside"<<endl;
		if(nCurve->IsInElement(pvPoint))
		{
		//	cout<<"yeah inside element"<<endl;
			distance.res   = -1.0;
			distance.data  = 0.0;
			distance.point = 0.0;
			return distance;
		}
	}

	
	/* compute closest point in this leaf segment box */
	int i;
	
	int nData = pBest->GetData();
	
	int nNumP = (int)nCurve->m_Buckets[nData].size();

	for(i = 0; i < nNumP; i++)
	{
		VECTOR2 b = nCurve->GetPoint(nData,i);
		b			= VECTOR2::createVector(b,pvPoint);
		d			= b.mag();
		if(distance.res > d)
		{
			distance.res   = d;
			distance.data  = nData;
			distance.point = i;
		}//end if
	}//end for



	
	/* loop through the candidates and check for potential improvement of the minimum distance */
	/* the current minimum distance can be improved if there is a candidate with a lower bound
	that is better than the current minimum distance */
	/* Here use binary tree or priority queue */

	for(lIter = lBFS.begin(); lIter != lBFS.end(); lIter++)
	{
		CBCNode *pBox = *lIter;
		if(pBox == pBest)
			continue;

		if(pBox->GetLowerBound(pvPoint) < distance.res)
		{
			
			nData = pBox->GetData();

			CApproxCurve *Curve = ApproxCurves[pBox->GetID()];

			nNumP = (int)Curve->m_Buckets[nData].size();

			/* loop to compute minimum distance for this bounding circle */
			for(i = 0; i < nNumP; i++)
			{
				VECTOR2 b = Curve->GetPoint(nData,i);
				b			= VECTOR2::createVector(b,pvPoint);
				d			= b.mag();
				if(distance.res > d)
				{
					distance.res   = d;
					distance.data  = nData;
					distance.point = i;
					nCurve         = Curve;
				}//end if
			}//end for

		}//end if

	}//end for
	return distance;
}//end lubSearch


///==================================================
/// 	BAB pruning for multiple objects
///==================================================
void  DistOps::breadthFirstSearchDistance(ListAABBND &Candidates, const VECTOR2 &vQuery, ListBCND &lBFS, CBCNode *&nBest)
{
	
	ListBCND::iterator lIter;
	ListAABBND::iterator Iterator;
	
	/* insert the candidate curve's segment boxes into the list */
	for(Iterator = Candidates.begin(); Iterator != Candidates.end(); Iterator++)
	{
		/* get pointer to the bounding box */
		CAABBNode *pNode = *Iterator;

		/* get pointer to the corresponding curve */
		CApproxCurve *pCurve =  pNode->Boxes().front()->m_paCurve;

		/* get number of segment boxes */
		int s = pCurve->GetNumCircles();

		/* insert segment boxes into the list */
		CBCTree &pTree = pCurve->m_BCTree;
		for(int i = 0; i < pTree.GetNumChildren(); i++)
		{
			lBFS.push_back(pTree.GetChild(i));
		}//end for

	}//end for

	/* init best segment box node with NULL */
	nBest = NULL;

	/* get size of queue */
	int vSize = (int)lBFS.size();
	
	/* count the number of leaves */
	int nLeafCount = 0;

	double lowerBound  =  std::numeric_limits<Real>::max();
	double upperBound  = -std::numeric_limits<Real>::max();

	int iter = 0;
	/* loop until there are only leaves in the list */
	while(vSize != nLeafCount)
	{
		nLeafCount = 0;
		lowerBound = std::numeric_limits<Real>::max();
		double *dLowerBounds = new double[vSize];
		/* find best lower bound */
		/* this has to be improved */
		int j = 0;
		for(lIter = lBFS.begin(); lIter != lBFS.end(); lIter++)
		{
			CBCNode *pNode = *lIter;
			dLowerBounds[j] = pNode->GetLowerBound( vQuery );
			if(lowerBound > dLowerBounds[j])
			{
				lowerBound = dLowerBounds[j];
				nBest = pNode;
			}//end if
			j++;
		}//end for
		
		/* get the best node's upper bound */
		upperBound = nBest->GetUpperBound(  vQuery );
		
		lIter = lBFS.begin();
		for(int i = 0; i < vSize; i++)
		{
			/* get the current node pointer */
			CBCNode *pNode = *lIter;

			/* if this node is the best node */
			if(pNode == nBest)
			{
				if(!pNode->IsLeaf())
				{
					lBFS.push_back(pNode->m_Children[0]);
					lBFS.push_back(pNode->m_Children[1]);
					lIter = lBFS.erase(lIter);
					continue;
				}
				else
				{
					nLeafCount++;
					lIter++;
					continue;
				}
			}//end if
			
			/* Can this subtree be pruned ? */
			else if(upperBound > dLowerBounds[i])
			{
				if(!pNode->IsLeaf())
				{
					lBFS.push_back(pNode->m_Children[0]);
					lBFS.push_back(pNode->m_Children[1]);
					
					lIter = lBFS.erase(lIter);
				}//end if
				else
				{
					nLeafCount++;
					lIter++;
				}
			}//end if
			else
			{
				lIter = lBFS.erase(lIter);
			}//end else
			
		}//end for
		vSize = (int)lBFS.size();
		
		delete[] dLowerBounds;

	}//end while
	
}//end breadthFirstSearch


///==================================================
/// 	LUBSEARCH FOR MULTIPLE OBJECTS
///		STARTS WITH SEARCH IN AABBTREE
///		THEN STARTS SEARCH IN BCTREE
///==================================================
SearchResult DistOps::AABBTreeSearch(CAABBTree *pTree, const VECTOR2 &vPoint, CApproxCurve *&nCurve, CurveVec &ApproxCurves)
{
	ListAABBND lCandidates;
	ListAABBND::iterator lIter;

	//BAB search in the AABBTree
	pTree->breadthFirstSearchDistance(vPoint, lCandidates);

	//BAB search int the Circletrees
	SearchResult result = LUBSearch(lCandidates, vPoint,nCurve, ApproxCurves);
	return result;
}//end AABBTreeSearch




// /* Benchmark code for the Windows version */
// #if defined WIN32
// 
// 
// 	
// 
// 	int i,j,k;
// 	i = j = k = 0;
// 
// 
// 
// 	for(i = 0; i < MAX_THREADS; i++)
// 	{
// 		resources[i].ApproxCurves = ApproxCurves;
// 		resources[i].NURBSCurves  = bcurve;
// 		resources[i].pTree        = new CAABBTree(pTree);
// 
// 		for(k = 0; k < GRIDSIZEX; k++)
// 			for(j = 0; j < GRIDSIZEY; j++)
// 			{
// 				resources[i].Grid[k][j] = grid[k][j];
// 			}
// 
// 	}//end for
// 
// 
// 
// 	LONGLONG frequency, lastTime, currentTime;
// 
// 	double timeElapsed, timeScale, tSum;
// 	tSum = 0;
// 
// 	QueryPerformanceFrequency((LARGE_INTEGER*) &frequency);
// 
// 	timeScale = 1.0/frequency;
// 
// 	double dLUB = 1.7E+308;
// 
// 
// 	VECTOR2 lastPoint(0,0);
// 	VECTOR2 lastRow(0,0);
// 
// 
// 	QueryPerformanceCounter((LARGE_INTEGER*) &lastTime);
// 
// #ifdef _SINGLETHREAD_
// 
// 	//for(int i = 0; i < GRIDSIZEX; i++)
//  //	{
//  //		for(int j = 0; j < GRIDSIZEY; j++)
//  //		{
// 	//		
// 	//		ListAABBND lCandidates;
// 
// 	//		pTree->breadthFirstSearchDistance(grid[i][j], lCandidates, dLUB);
// 
//  //			double *res = COLUBSearch(lCandidates, grid[i][j], nCurve, ApproxCurves,bcurve, dLUB);
//  // 			
// 	//		lastPoint = nCurve->GetPoint(int(res[1]),int(res[2]));
//  //			
//  //			lastPoint = VECTOR2::createVector(lastPoint, grid[i][(j+1)%GRIDSIZEY]);
//  //			
//  //			dLUB = lastPoint.mag();
// 
// 	//		delete[] res;
// 
//  //		}//end for j
//  //	}//end for i
// 
// 
// 	CRootFinder rFA;
// 
// 	int iConv = 0;
// 
// 	double lub[GRIDSIZEX][GRIDSIZEY];
// 	double newton[GRIDSIZEX][GRIDSIZEY];
// 
// 	for(int i = 0; i < GRIDSIZEX; i++)
// 	{
// 		for(int j = 0; j < GRIDSIZEY; j++)
// 		{
// 			list<CAABBNode*> lCandidates;
// 			list<CBCNode*> lBFS;
// 			CBCNode *nBest = NULL;
// 			CApproxCurve *acurve = NULL;
// 			pTree->breadthFirstSearchDistance(grid[i][j], lCandidates, dLUB);
// 			SearchResult result = COLUBSearch(lCandidates, grid[i][j], acurve, ApproxCurves, bcurve, dLUB);
// 			
// 			lastPoint = acurve->GetPoint(result.data, result.point);
// 
// 			lastPoint = VECTOR2::createVector(lastPoint, grid[i][(j+1)%GRIDSIZEY]);
// 		
// 			dLUB = lastPoint.mag();
// 	/*		lub[i][j] = result.res;
// 			newton[i][j] = result.res;
// 			//double u = ApproxCurves[0]->m_dParamBuckets[result.data][result.point].dP;
// 	/*		if(result[0] < ApproxCurves[0]->Radius())
// 			{*/
// //				if(rFA.ONX_Newton(myNurbs, grid[i][j], u, 1e-14, 10))
// //				{
// //					iConv++;
// //					VECTOR2 NearestPoint = myNurbs->CoxDeBoor(u);
// //					double d = VECTOR2::createVector(NearestPoint, grid[i][j]).mag();
// //					if(result.res > 0)
// //    					newton[i][j] = d;
// //				}
// ////			}
// 
// 
// 
// 		}//end for j
// 	}//end for i
// 
// 	//cout<<"converged: "<<iConv<<endl;
// 
// 	//double n = GRIDSIZEX * GRIDSIZEY;
// 	//double dSum = 0.0;
// 	//for(k = 0; k < GRIDSIZEX; k++)
// 	//	for(int l = 0; l < GRIDSIZEY; l++)
// 	//	{
// 	//		dSum += fabs(lub[k][l] - newton[k][l]);
// 	//	}
// 
// 	//dSum/=n;
// 	//cout.precision(16);
// 	//cout<<" Average error: "<<dSum<<endl;
// 
// #else
// 
// 	cout<<"start thread"<<endl;
// 
// 	HANDLE threadHandle1;
// 
// 	threadHandle1 = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)threadProc1, NULL, 0, NULL);
// 
// 	HANDLE threadHandle2;
// 
// 	threadHandle2 = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)threadProc2, NULL, 0, NULL);	
// 
// 	WaitForSingleObject(threadHandle1, INFINITE);
// 
// 	WaitForSingleObject(threadHandle2, INFINITE);
// 
// #endif
// 
// 	QueryPerformanceCounter((LARGE_INTEGER*) &currentTime);
// 	timeElapsed = (currentTime - lastTime) * timeScale;
// 
// 	cout<<"combined search time: "<<timeElapsed<<endl;
// 
// /* Benchmark code for the LINUX/UNIX version */
// #else
// 
// #ifdef USETQ
// 	QTime timer;
// #endif
// 	int i,j;
// 	VECTOR2 lastPoint(0,0);
// 	VECTOR2 lastRow(0,0);
// 
// 	ListBCND Candidates;
// 	CBCNode *Best;
// 	
// 	double dLUB = 1.7E+308;
// 
// /* Benchmark code for the singlethreaded version */
// #ifdef _SINGLETHREAD_
// #ifdef USETQ
// 	timer.start();
// #endif	
// 
// 	for(i = 0; i < GRIDSIZEX; i++)
// 	{
// 		for(j = 0; j < GRIDSIZEY; j++)
// 		{
// 
// // 			ListAABBND lCandidates;
// // 						
// //  			pTree->breadthFirstSearchDistance(grid[i][j], lCandidates, dLUB);
// // 
// //  			double *res = COLUBSearch(lCandidates, grid[i][j], nCurve, ApproxCurves,bcurve, dLUB);
// //  			
// //  			lastPoint = nCurve->GetPoint(int(res[1]),int(res[2]));
// //  			
// //  			lastPoint = VECTOR2::createVector(lastPoint, grid[i][(j+1)%GRIDSIZEY]);
// //  			
// //  			dLUB = lastPoint.mag();
// // 
// // 			delete[] res;
// 
// 		}//end for j
// 	}//end for i	
// 
// 	/* Benchmark code for the multithreaded version */
// #else
// 
// 	
// 	pthread_t Threads[MAX_THREADS];
// 	void* ret[MAX_THREADS];
// 
// 	int nThreads = 4;
// 
// 	pthread_create(&Threads[0], NULL, threadProc1, NULL);
// 	pthread_create(&Threads[1], NULL, threadProc2, NULL);
// 	//pthread_create(&Threads[2], NULL, threadProc3, NULL);
// 	//pthread_create(&Threads[3], NULL, threadProc4, NULL);
// 
// 	cout<<"starting search..."<<endl;
// #ifdef USEQT
// 	timer.start();
// #endif
// 
// #ifdef USEQT	
// 	for( i=0; i < 2; i++)
// 		pthread_join(Threads[i], &ret[i]);
// 	cout<<"Search time: "<<timer.elapsed()<<endl;
// #endif
// 	
// #endif


#ifdef WIN32

void threadProc1()
{

	int i,j;
	VECTOR2 lastPoint(0,0);
	VECTOR2 lastRow(0,0);

	SearchResult lRes;
	
	double dLUB = 1.7E+308;
	
	DistOps op;

	CApproxCurve *nCurve = NULL;

	int GRID2 = GRIDSIZEX/2;

	for(i = 0; i < GRID2;i++)
	{
		for(j = 0; j < GRIDSIZEY; j++)
		{
			
 			ListAABBND lCandidates;
			
			resources[0].pTree->breadthFirstSearchDistance(resources[0].Grid->Point(i,j), lCandidates, lRes.res);
			
 			SearchResult res = op.COLUBSearch(lCandidates, resources[0].Grid->Point(i,j),
										 nCurve, resources[0].ApproxCurves, lRes);

			lastPoint = nCurve->GetPoint(res.data,res.point);
 			
 			lastPoint = VECTOR2::createVector(lastPoint, resources[0].Grid->Point(i,(j+1)%GRIDSIZEY));
 			
 			lRes.dLUB = lastPoint.mag();

		}//end for j
	}//end for i	

}//end threadProc1

void threadProc2()
{

	int i,j;
	VECTOR2 lastPoint(0,0);
	VECTOR2 lastRow(0,0);
	
	double dLUB = 1.7E+308;
	
	DistOps op;

	SearchResult lRes;

	CApproxCurve *nCurve = NULL;

	int iGridCeil = (GRIDSIZEX/2) + 1;

	for(i = iGridCeil; i < GRIDSIZEX; i++)
	{
		for(j = 0; j < GRIDSIZEY; j++)
		{
			//cout<<"i "<<i<<"j "<<j<<endl;
 			ListAABBND lCandidates;
			
 			resources[1].pTree->breadthFirstSearchDistance(resources[1].Grid->Point(i,j), lCandidates, dLUB);
			
 			SearchResult res = op.COLUBSearch(lCandidates, resources[1].Grid->Point(i,j),
										 nCurve, resources[1].ApproxCurves, lRes);

			lastPoint = nCurve->GetPoint( res.data, res.point);
			 			
 			lastPoint = VECTOR2::createVector(lastPoint, resources[1].Grid->Point(i,(j+1)%GRIDSIZEY));
 			
 			dLUB = lastPoint.mag();

		}//end for j
	}//end for i	

}//end threadProc2

#endif

#ifdef LINUX

// void *threadProcNewton1(void *p1)
// {
// 
// 	int i,j;
// 	VECTOR2 lastPoint(0,0);
// 	VECTOR2 lastRow(0,0);
// 	
// 	CRootFinder rFA;
// 	
// 	double dLUB = 1.7E+308;
// 	
// 	DistOps op;
// 
// 	CApproxCurve *nCurve = NULL;
// 	
// 	int GRID2 = GRIDSIZEX/int(MAX_THREADS);
// 
// 	for(i = 0; i < GRID2;i++)
// 	{
// 		for(j = 0; j < GRIDSIZEY; j++)
// 		{
// 			//cout<<"i "<<i<<"j "<<j<<endl;
// 			ListAABBND lCandidates;
// 			
// 			resources[0].pTree->breadthFirstSearchDistance(resources[0].Grid->Point(i,j), lCandidates, dLUB);
// 			
// 			SearchResult res = op.COLUBSearch(lCandidates, resources[0].Grid->Point(i,j),
// 											  nCurve, resources[0].ApproxCurves, dLUB);
// 
// 			lastPoint = nCurve->GetPoint(res.data,res.point);
//  			
// 			lastPoint = VECTOR2::createVector(lastPoint, resources[0].Grid->Point(i,(j+1)%GRIDSIZEY));
// 			
// 			resources[0].Grid->SetDistance(i,j,res.res);
// 			
// 			CONXNurbs *Nurbs = resources[0].NURBS[nCurve->ID()];
// 			
// 			double u = nCurve->m_dParamBuckets[res.data][res.point].dP;
// 				
// 			if(rFA.ONX_Newton(Nurbs, resources[0].Grid->m_Grid[i][j].coord, u, 1e-7, 30))
// 			{
// 				
// 				
// 				VECTOR2 NearestPoint = Nurbs->CoxDeBoor(u);
// 				double d = VECTOR2::createVector(NearestPoint, resources[0].Grid->Point(i,j)).mag();
// 				if(res.res > 0)
// 				{
// 					resources[0].Grid->m_Grid[i][j].newton = d;
// 				
// 				}
// 			}			
//  			
// 			dLUB = lastPoint.mag();
// 
// 		}//end for j
// 	}//end for i	
// 
// 
// }//end threadProc1
// 
// void *threadProcNewton2(void *p2)
// {
// 
// 	int i,j;
// 	VECTOR2 lastPoint(0,0);
// 	VECTOR2 lastRow(0,0);
// 	
// 	CRootFinder rFA;
// 	
// 	double dLUB = 1.7E+308;
// 	
// 	DistOps op;
// 
// 	CApproxCurve *nCurve = NULL;
// 
// 	int iGridCeil = (GRIDSIZEX/(int)MAX_THREADS);
// 	int iEnd      = (GRIDSIZEX * 2/(int)MAX_THREADS); 
// 
// 	for(i = iGridCeil; i < iEnd; i++)
// 	{
// 		for(j = 0; j < GRIDSIZEY; j++)
// 		{
// 			//cout<<"i "<<i<<"j "<<j<<endl;
// 			ListAABBND lCandidates;
// 			
// 			resources[1].pTree->breadthFirstSearchDistance(resources[1].Grid->Point(i,j), lCandidates, dLUB);
// 			
// 			SearchResult res = op.COLUBSearch(lCandidates, resources[1].Grid->Point(i,j),
// 											  nCurve, resources[1].ApproxCurves, dLUB);
// 
// 			lastPoint = nCurve->GetPoint( res.data, res.point);
// 			 			
// 			lastPoint = VECTOR2::createVector(lastPoint, resources[1].Grid->Point(i,(j+1)%GRIDSIZEY));
// 			
// 			resources[1].Grid->SetDistance(i,j,res.res);
// 			
// 			CONXNurbs *Nurbs = resources[1].NURBS[nCurve->ID()];
// 			
// 			double u = nCurve->m_dParamBuckets[res.data][res.point].dP;
// 				
// 			if(rFA.ONX_Newton(Nurbs, resources[1].Grid->m_Grid[i][j].coord, u, 1e-7, 30))
// 			{
// 				
// 				
// 				VECTOR2 NearestPoint = Nurbs->CoxDeBoor(u);
// 				double d = VECTOR2::createVector(NearestPoint, resources[1].Grid->Point(i,j)).mag();
// 				if(res.res > 0)
// 				{
// 					resources[1].Grid->m_Grid[i][j].newton = d;
// 				
// 				}
// 			}						
//  			
// 			dLUB = lastPoint.mag();
// 
// 		}//end for j
// 	}//end for i		
// 
// }//end threadProc2


#endif


///==================================================
/// 	COLUBSEARCH FOR MULTIPLE OBJECTS
///		SEARCHES IN BCTREE
///==================================================
SearchResult DistOps::COLUBSearch(ListAABBND &Candidates, const VECTOR2 &pvPoint, CApproxCurve *&nCurve,
				     CurveVec &ApproxCurves, SearchResult &lRes)
{

	SearchResult distance;
	double d = std::numeric_limits<Real>::max();
	double mp = -1.0;
	/* initialize distance with default values */
	distance.res = std::numeric_limits<Real>::max();
	distance.data = -1;

	double dLUB = lRes.dLUB;

	/* the list will store the candidate leaves of the bspline curves segment box tree*/
	CBCNode *pBest = NULL;
	CBCNode *nextBest = NULL;

	list< NodeBound > myBFS;
	list< NodeBound >::iterator myIter;
	
	/* pruning search, pBest will hold the nearest segment box */
	/* lower bound/upper bound search to find the candidate leaves */
	//breadthFirstSearchDistance(Candidates, pvPoint, lCandidates, pBest);
	/////////////////////////////////////////////////////////////////////
	double lowerBound     =  std::numeric_limits<Real>::max();
	double nextLowerBound =  std::numeric_limits<Real>::max();
	double upperBound	  = -std::numeric_limits<Real>::max();
	
	ListAABBND::iterator Iterator;

	/* init best segment box node with NULL */
	pBest = NULL;

	/* insert the candidate curve's segment boxes into the list */
	for(Iterator = Candidates.begin(); Iterator != Candidates.end(); Iterator++)
	{
		/* get pointer to the bounding box */
		CAABBNode *pNode = *Iterator;
		
		for(int k = 0; k < pNode->GetItemCount();  k++)
		{

			/* get pointer to the corresponding curve */
			CApproxCurve *pCurve =  pNode->Boxes()[k]->m_paCurve;
	
			/* get number of segment boxes */
			int s = pCurve->GetNumCircles();
	
			CBCTree &pTree = pCurve->m_BCTree;
			/* insert segment boxes into the list */
			for(int i = 0; i < s; i++)
			{
				double daLowerBound = pTree.GetChild(i)->GetLowerBound(pvPoint);
				if(daLowerBound <= dLUB)
				{
					if(lowerBound > daLowerBound)
					{
						lowerBound = daLowerBound;
						pBest = pTree.GetChild(i);
					}
					NodeBound tPair;
					tPair.tNode = pTree.GetChild(i);
					tPair.dLBound = daLowerBound;
					myBFS.push_back(tPair);
				}//end if
			}//end for
			
		}//end for
	}//end for

	/* get size of queue */
	int vSize = (int)myBFS.size();
	
	/* count the number of leaves */
	int nLeafCount = 0;
	
	
	int iter = 0;
	/* loop until there are only leaves in the list */
	while(vSize != nLeafCount)
	{
		nLeafCount = 0;

		nextLowerBound = 1.7E+308;

		//DO a best-first-search
		/* get the best node's upper bound */
		double bestUpper = pBest->GetUpperBound(  pvPoint );
		upperBound = (dLUB > bestUpper) ? bestUpper : dLUB;
		
		myIter = myBFS.begin();
		for(int i = 0; i < vSize; i++)
		{
			/* get the current node pointer */
			NodeBound NodeAndBound = *myIter;

			/* if this node is the best node */
			if(NodeAndBound.tNode == pBest)
			{
				if(!NodeAndBound.tNode->IsLeaf())
				{
					NodeBound NodesAndBounds[2];
					NodesAndBounds[0].dLBound = NodeAndBound.tNode->m_Children[0]->GetLowerBound(pvPoint);
					NodesAndBounds[1].dLBound = NodeAndBound.tNode->m_Children[1]->GetLowerBound(pvPoint);
					NodesAndBounds[0].tNode   = NodeAndBound.tNode->m_Children[0];
					NodesAndBounds[1].tNode   = NodeAndBound.tNode->m_Children[1];

					if(NodesAndBounds[0].dLBound < nextLowerBound)
					{
						nextLowerBound = NodesAndBounds[0].dLBound;
						nextBest       = NodesAndBounds[0].tNode;
					}

					if(NodesAndBounds[1].dLBound < nextLowerBound)
					{
						nextLowerBound = NodesAndBounds[1].dLBound;
						nextBest       = NodesAndBounds[1].tNode;
					}

					if(NodesAndBounds[0].dLBound <= upperBound)
					{
						myBFS.push_back(NodesAndBounds[0]);
					}//end if

					if(NodesAndBounds[1].dLBound <= upperBound)
					{
						myBFS.push_back(NodesAndBounds[1]);
					}//end if

					myIter = myBFS.erase(myIter);
					continue;
				}
				else
				{
					nLeafCount++;
					myIter++;
					continue;
				}
			}//end if
			
			/* Can this subtree be pruned ? */
			if(upperBound > NodeAndBound.dLBound )
			{
				if(!NodeAndBound.tNode->IsLeaf())
				{
					NodeBound NodesAndBounds[2];
					NodesAndBounds[0].dLBound = NodeAndBound.tNode->m_Children[0]->GetLowerBound(pvPoint);
					NodesAndBounds[1].dLBound = NodeAndBound.tNode->m_Children[1]->GetLowerBound(pvPoint);
					NodesAndBounds[0].tNode   = NodeAndBound.tNode->m_Children[0];
					NodesAndBounds[1].tNode   = NodeAndBound.tNode->m_Children[1];

					if(NodesAndBounds[0].dLBound < nextLowerBound)
					{
						nextLowerBound = NodesAndBounds[0].dLBound;
						nextBest       = NodesAndBounds[0].tNode;
					}

					if(NodesAndBounds[1].dLBound < nextLowerBound)
					{
						nextLowerBound = NodesAndBounds[1].dLBound;
						nextBest       = NodesAndBounds[1].tNode;
					}
					if(NodesAndBounds[0].dLBound <= upperBound)
						myBFS.push_back(NodesAndBounds[0]);
					if(NodesAndBounds[1].dLBound <= upperBound)
						myBFS.push_back(NodesAndBounds[1]);

					
					myIter = myBFS.erase(myIter);
				}//end if
				else
				{
					//tree.insert(NodeAndBound);
					//myIter = myBFS.erase(myIter);
					nLeafCount++;
					myIter++;
				}
			}//end if
			else
			{
				myIter = myBFS.erase(myIter);
			}//end else
			
		}//end for

		if(nextBest)
			pBest    = nextBest;

		nextBest = NULL;
		vSize = (int)myBFS.size();
			
	}//end while
	
	/////////////////////////////////////////////////////////////////////

//#ifdef _DEBUG
	if(myBFS.empty())
	{
//		cout<<"All nodes were rejected, an error will occur..."<<endl;
		pBest = lRes.pBest;
	}//end if
	else
		lRes.pBest = pBest;
	
//#endif

	/* get the pointer to the curve that the best leaf belongs to */
	nCurve = ApproxCurves[pBest->GetID()];

	if(nCurve->m_bBox->Inside(pvPoint))
	{
		//cout<<"inside"<<endl;
		if(nCurve->IsInElement(pvPoint))
		//if(numIntersections(bcurve[nCurve->ID()], pvPoint))
		{
		//	cout<<"yeah inside element"<<endl;
			distance.res = -1.0;
			distance.data = 0;
			distance.point = 0;
			return distance;
		}
		
	}
	
	/* compute closest point in this leaf segment box */
	int i;
	
	int nData = pBest->GetData();

#ifdef _DEBUG
	if(nData < 0)
	{
		cout<<"Invalid data"<<endl;
	}//end if
#endif
	
	int nNumP = (int)nCurve->m_Buckets[nData].size();

	for(i = 0; i < nNumP; i++)
	{
		VECTOR2 b = nCurve->GetPoint(nData,i);
		b			= VECTOR2::createVector(b,pvPoint);
		d			= b.mag();
		if(distance.res > d)
		{
			distance.res   = d;
			distance.data  = static_cast<double>(nData);
			distance.point = static_cast<double>(i);
		}//end if
	}//end for
	
	
	/* check the remaining canidates for a potential improvement */
	for(myIter = myBFS.begin(); myIter != myBFS.end(); myIter++)
	{
		NodeBound nb  = *myIter;
		CBCNode *pBox = nb.tNode;
		if(pBox == pBest)
			continue;

		/* is an improvement possible ? */
		if(nb.dLBound < distance.res)
		{
			
			nData = pBox->GetData();

			CApproxCurve *Curve = ApproxCurves[pBox->GetID()];

			nNumP = (int)Curve->m_Buckets[nData].size();

			//loop to compute minimum distance for this bounding circle 
			for(i = 0; i < nNumP; i++)
			{
				VECTOR2 b = Curve->GetPoint(nData,i);
				b			= VECTOR2::createVector(b,pvPoint);
				d			= b.mag();
				if(distance.res > d)
				{
					distance.res   = d;
					distance.data  = static_cast<double>(nData);
					distance.point = static_cast<double>(i);
					nCurve      = Curve;
				}//end if
			}//end for

		}//end if

	}//end for
	
	//tree.destroyTree();
	return distance;

}//end lubSearch

SearchResult   DistOps::LUBSearch(ListAABBND &Candidates, const VECTOR2 &pvPoint, CApproxCurve *&nCurve,
			  CurveVec &ApproxCurves, CBCNode *&pBest, ListBCND &lBFS)
{
	SearchResult distance;
	double d = 1.7E+308;
	double mp = -1.0;
	/* initialize distance with default values */
	distance.res  = 1.7E+308;
	distance.data = -1;


	/* the list will store the candidate leaves of the bspline curves segment box tree*/
	
	ListBCND::iterator lIter;
	pBest = NULL;
	
	/* pruning search, pBest will hold the nearest segment box */
	/* lower bound/upper bound search to find the candidate leaves */
	breadthFirstSearchDistance(Candidates, pvPoint, lBFS, pBest);

	/* get the pointer to the curve that the best leaf belongs to */
	nCurve = ApproxCurves[pBest->GetID()];


	if(nCurve->m_bBox->inside(pvPoint))
	{
		//cout<<"inside"<<endl;
		if(nCurve->Inside(pvPoint))
		{
		//	cout<<"yeah inside element"<<endl;
			distance.res   = -1.0;
			distance.data  = 0;
			distance.point = 0;
			return distance;
		}
	}

	
	/* compute closest point in this leaf segment box */
	int i;
	
	int nData = pBest->GetData();
	
	int nNumP = (int)nCurve->m_Buckets[nData].size();

	for(i = 0; i < nNumP; i++)
	{
		VECTOR2 b = nCurve->GetPoint(nData,i);
		b			= VECTOR2::createVector(b,pvPoint);
		d			= b.mag();
		if(distance.res > d)
		{
			distance.res   = d;
			distance.data  = nData;
			distance.point = i;
		}//end if
	}//end for



	
	/* loop through the candidates and check for potential improvement of the minimum distance */
	/* the current minimum distance can be improved if there is a candidate with a lower bound
	that is better than the current minimum distance */
	/* Here use binary tree or priority queue */

	for(lIter = lBFS.begin(); lIter != lBFS.end(); lIter++)
	{
		CBCNode *pBox = *lIter;
		if(pBox == pBest)
			continue;

		if(pBox->GetLowerBound(pvPoint) < distance.res)
		{
			
			nData = pBox->GetData();

			CApproxCurve *Curve = ApproxCurves[pBox->GetID()];

			nNumP = (int)Curve->m_Buckets[nData].size();

			/* loop to compute minimum distance for this bounding circle */
			for(i = 0; i < nNumP; i++)
			{
				VECTOR2 b = Curve->GetPoint(nData,i);
				b			= VECTOR2::createVector(b,pvPoint);
				d			= b.mag();
				if(distance.res > d)
				{
					distance.res   = d;
					distance.data  = nData;
					distance.point = i;
					nCurve      = Curve;
				}//end if
			}//end for

		}//end if

	}//end for
	return distance;
}

///==================================================
/// 	COHERENCY LUBSEARCH FOR SINGLE OBJECT
///		Makes use of coherency information, a global     
///		upper bound is passed via the parameter
///		dLUB
///==================================================
SearchResult DistOps::COLUBSearch(CApproxCurve *pCurve, const VECTOR2 &pvPoint, ListBCND &lBFS, CBCNode *&pBest,double dLUB)
{

	SearchResult distance;
	double d = std::numeric_limits<Real>::max();
	double mp = -1.0;
	/* initialize distance with default values */
	distance.res  = std::numeric_limits<Real>::max();
	distance.data = -1;
	/* the list will store the candidate leaves of the bspline curves segment box tree*/
	ListBCND::iterator lIter;
	pBest = NULL;
	
	int i;

	/* check whether the point is inside the curve */
	if(pCurve->m_bBox->inside(pvPoint))
	{
		
		if(pCurve->IsInElement(pvPoint))
		{
			distance.res   = -1.0;
			distance.data  = 0;
			distance.point = 0;
			return distance;
		}//end if
		
	}//end if

	CBCTree &pTree = pCurve->m_BCTree;
	/* insert segment boxes into the list */
	for(i = 0; i < pTree.GetNumChildren(); i++)
	{
		lBFS.push_back(pTree.GetChild(i));
	}//end for

	
	/* init best segment box node with NULL */
	pBest = NULL;

	/* get size of queue */
	int vSize = (int)lBFS.size();
	
	/* count the number of leaves */
	int nLeafCount = 0;

	double lowerBound  =  std::numeric_limits<Real>::max();
	double upperBound  = -std::numeric_limits<Real>::max();
	
	int iter = 0;
	/* loop until there are only leaves in the list */
	while(vSize != nLeafCount)
	{
		nLeafCount = 0;
		lowerBound = std::numeric_limits<Real>::max();
		double *dLowerBounds = new double[vSize];
		/* find best lower bound */
		int j = 0;
		for(lIter = lBFS.begin(); lIter != lBFS.end(); lIter++)
		{
			CBCNode *pNode = *lIter;
			dLowerBounds[j] = pNode->GetLowerBound( pvPoint );
			if(lowerBound > dLowerBounds[j])
			{
				lowerBound = dLowerBounds[j];
				pBest = pNode;
			}//end if
			j++;
		}//end for
		
		/* get the best node's upper bound */
		double bestUpper = pBest->GetUpperBound(  pvPoint );
		upperBound = (dLUB > bestUpper) ? bestUpper : dLUB;
		
		lIter = lBFS.begin();
		for(i = 0; i < vSize; i++)
		{
			/* get the current node pointer */
			CBCNode *pNode = *lIter;

			/* if this node is the best node */
			if(pNode == pBest)
			{
				if(!pNode->IsLeaf())
				{
					lBFS.push_back(pNode->m_Children[0]);
					lBFS.push_back(pNode->m_Children[1]);
					lIter = lBFS.erase(lIter);
					continue;
				}
				else
				{
					nLeafCount++;
					lIter++;
					continue;
				}
			}//end if
			
			/* Can this subtree be pruned ? */
			if(upperBound > dLowerBounds[i])
			{
				if(!pNode->IsLeaf())
				{
					lBFS.push_back(pNode->m_Children[0]);
					lBFS.push_back(pNode->m_Children[1]);
					
					lIter = lBFS.erase(lIter);
				}//end if
				else
				{
					nLeafCount++;
					lIter++;
				}
			}//end if
			else
			{
				lIter = lBFS.erase(lIter);
			}//end else
			
		}//end for
		vSize = (int)lBFS.size();
		delete[] dLowerBounds;
	}//end while
	
	/* compute closest point in this leaf segment box */
	int nData = pBest->GetData();
	
	int nNumP = (int)pCurve->m_Buckets[nData].size();

	for(i = 0; i < nNumP; i++)
	{
		VECTOR2 b = pCurve->GetPoint(nData,i);
		b			= VECTOR2::createVector(b,pvPoint);
		d			= b.mag();
		if(distance.res > d)
		{
			distance.res   = d;
			distance.data  = nData;
			distance.point = i;
		}//end if
	}//end for

	/* check the other candidates if an improvement is possible */
	for(lIter = lBFS.begin(); lIter != lBFS.end(); lIter++)
	{
		CBCNode *pBox = *lIter;
		if(pBox == pBest)
			continue;

		if(pBox->GetLowerBound(pvPoint) < distance.res)
		{
			nData = pBox->GetData();
			nNumP = (int)pCurve->m_Buckets[nData].size();
			/* loop to compute minimum distance for this bounding circle */
			for(i = 0; i < nNumP; i++)
			{
				VECTOR2 b = pCurve->GetPoint(nData,i);
				b			= VECTOR2::createVector(b,pvPoint);
				d			= b.mag();
				if(distance.res > d)
				{
					distance.res   = d;
					distance.data  = nData;
					distance.point = i;
				}//end if
			}//end for

		}//end if

	}//end for
	return distance;

}//end lubSearch


/*!
    \fn DistOps::DebugSearch(CApproxCurve *Curve, const VECTOR2 &pvPoint, ListBCND &lCandidates, CBCNode *&pBest)
 */
DebugResult DistOps::DebugSearch(CApproxCurve *Curve, const VECTOR2 &pvPoint, ListBCND &lCandidates, CBCNode *&pBest)
{
	double d = 1.7E+308;
	double mp = -1.0;
	/* initialize distance with default values */
	DebugResult distance;
	distance.res = 1.7E+308;
	distance.data = -1;
	
	distance.vPoint = pvPoint;
	
	ListBCND::iterator lIter;
	pBest = NULL;

	distance.pBestPoint = NULL;
	
	if(Curve->m_bBox->inside(pvPoint))
	{
		//cout<<"inside"<<endl;
		if(Curve->IsInElement(pvPoint))
		{
		//	cout<<"yeah inside element"<<endl;
			distance.res   = -1.0;
			distance.data  = 0;
			distance.point = 0;
			return distance;

		}
	}
	
	/* breadthFirstSearch */
	BreadthFirstSearchDebug(Curve, pvPoint, lCandidates, pBest,distance);

	/* compute closest point in this circle */
	int i;

	///* Get the address of the node's data *///
	int nData = pBest->GetData();

	///* Get the number of points contained in the BV *///
	int nNumP = (int)Curve->m_Buckets[nData].size();

	for(i = 0; i < nNumP; i++)
	{
		VECTOR2 b = Curve->GetPoint(nData,i);
		b			= VECTOR2::createVector(b,pvPoint);
		d			= b.mag();
		if(distance.res > d)
		{
			distance.res   = d;
			distance.data  = static_cast<double>(nData);
			distance.point = static_cast<double>(i);
		}//end if
	}//end for

	distance.candidates[3] = lCandidates;
	distance.pBest[3] = pBest;
	
	for(lIter = lCandidates.begin();lIter != lCandidates.end();lIter++)
	{
		CBCNode *pBox = *lIter;
		if(pBox == pBest)
			continue;

		if(pBox->GetLowerBound(pvPoint) < distance.res)
		{
			
			nData = pBox->GetData();

			nNumP = (int)Curve->m_Buckets[nData].size();

			/* loop to compute minimum distance for this bounding circle */
			for(i = 0; i < nNumP; i++)
			{
				VECTOR2 b = Curve->GetPoint(nData,i);
				b			= VECTOR2::createVector(b,pvPoint);
				d			= b.mag();
				if(distance.res > d)
				{
					distance.res = d;
					distance.data = static_cast<double>(nData);
					distance.point = static_cast<double>(i);
					distance.pBestPoint = pBox;
				}//end if
			}//end for

		}//end if

	}//end for
	
	return distance;
	
}//DebugSearch

void DistOps::BreadthFirstSearchDebug(CApproxCurve* pCurve, const VECTOR2 &vQuery, ListBCND &lBFS, CBCNode *&nBest, DebugResult &data)
{
	ListBCND::iterator lIter;
	double lowerBound  = 1.7E+308;
	double upperBound  = -1.7E+308;
	
	CBCTree &pTree = pCurve->m_BCTree;
	for(int i = 0; i < pTree.GetNumChildren(); i++)
		lBFS.push_back(pTree.GetChild(i));
	
    nBest = NULL;

	cout<<"Number of Top Level Nodes: "<<pTree.GetNumChildren()<<endl;
	
	int vSize = (int)lBFS.size();

	int nLeafCount = 0;
	int iters = 0;
	/* loop until there are only leaves in the list */
	while(vSize != nLeafCount)
	{
		int pruned = 0;
		nLeafCount = 0;
		int j = 0;
		lowerBound = 1.7E+308;
		double *dLowerBounds = new double[vSize];
		/* find best upper bound */
		for(lIter = lBFS.begin(); lIter != lBFS.end(); lIter++)
		{
			CBCNode *pNode = *lIter;
			dLowerBounds[j] = pNode->GetLowerBound(vQuery);
			if(lowerBound > dLowerBounds[j])
			{
				lowerBound = dLowerBounds[j];
				nBest = pNode;
				data.pBest[iters] = nBest;
			}//end if
			j++;
		}//end for
		
		cout<<"Elements in queue: "<<vSize<<endl;
		
		/* get upper bound for best element */
		upperBound = nBest->GetUpperBound(vQuery);
		data.vUBounds[iters] = nBest->GetUpperVec();
		
		lIter = lBFS.begin();
		for(int i = 0; i < vSize; i++)
		{
			
			CBCNode *pNode = *lIter;
			if(pNode == nBest)
			{
				if(!pNode->IsLeaf())
				{
					lBFS.push_back(pNode->m_Children[0]);
					lBFS.push_back(pNode->m_Children[1]);
					lIter = lBFS.erase(lIter);
					data.candidates[iters].push_back(pNode);
					continue;
				}
				else
				{
					nLeafCount++;
					lIter++;
					continue;
				}
			}//end if
			
			/* Can this subtree be pruned ? */
			if(upperBound > dLowerBounds[i])
			{
				
				if(!pNode->IsLeaf())
				{
					lBFS.push_back(pNode->m_Children[0]);
					lBFS.push_back(pNode->m_Children[1]);
					data.candidates[iters].push_back(pNode);
					lIter = lBFS.erase(lIter);
				}//end if
				else
				{
					nLeafCount++;
					lIter++;
				}
			}//end if
			else
			{
				pruned++;
				lIter = lBFS.erase(lIter);
			}//end else
			
		}//end for
		vSize = (int)lBFS.size();
		delete[] dLowerBounds;
		cout<<"Iteration: "<<iters<<endl;
		cout<<"Elements pruned: "<<pruned<<endl;
		cout<<"-----------------"<<endl;
		cout<<endl;
		iters++;
	}//end while
	
}//BreadthFirstSearchDebug

SearchResult DistOps::BruteDistance_SIMD(CApproxCurve *pCurve, const VECTOR2 &vQuery)
{

	SearchResult distance;
	
#ifdef WIN32
	
	double d = 1.7E+308;

	int i;
	
	__declspec(align(16)) float vXIJ[4];
	__declspec(align(16)) float vXQ[4];
	__declspec(align(16)) float vResult[4];

	vXQ[0]=vQuery.y;
	vXQ[1]=vQuery.x;
	vXQ[2]=vQuery.y;
	vXQ[3]=vQuery.x;

	/* initialize distance with default values */
	distance.res  = 1.7E+308;
	distance.data = -1;

	//Real *dValues = new Real[pCurve->GetResolution()];

	if(pCurve->m_bBox->inside(vQuery))
	{
		if(pCurve->IsInElement(vQuery))
		{
			distance.res   = -1.0;
			distance.data  = 0;
			distance.point = 0;
			return distance;
		}
	}

	VECTOR2 c;
	
	VECTOR2 *vSamples = pCurve->GetSamples();
	
	int nNumP = pCurve->GetResolution();

	for(int i = 0; i < nNumP; i+=2)
	{
		VECTOR2 vHigh = vSamples[i];
		VECTOR2 vLow  = vSamples[i+1];
		vXIJ[0]       = vLow.y;
		vXIJ[1]       = vLow.x;
		vXIJ[2]       = vHigh.y;
		vXIJ[3]       = vHigh.x;
		_asm
		{
			//move vXIJ into xmm0
			movaps xmm0, vXIJ
			//subtract vXQ from vXIJ
			subps  xmm0, vXQ
			//compute the squared vector components
			//(x0 * x0), (y0 * y0), (x1 * x1), (y1 * y1)
			mulps  xmm0, xmm0
			//copy xmm0 into xmm1
			movaps xmm1, xmm0
			//shuffle the registers
			shufps xmm1, xmm0, SIMD_SHUFFLE(0x02,0x03,0x00,0x01)
			//simply add the registers
			addps xmm0, xmm1
			//transfer the result
			movaps vResult, xmm0
		}
		
		d			  = vResult[1];
		if(distance.res > d)
		{
			distance.res = d;
			distance.point = i+1;
		}//end if

		d			  = vResult[3];
		if(distance.res > d)
		{
			distance.res = d;
			distance.point = i;
		}//end if

	}//end for

	distance.res = sqrt(distance.res);

#endif
	
	return distance;

}//end BruteDistance_SIMD

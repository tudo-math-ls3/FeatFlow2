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
#include <limits>
#include <list>

#include "distops3.h"
#include "intersectorray3tri3.h"



CDistOps3::CDistOps3(void)
{
}

CDistOps3::~CDistOps3(void)
{
}

Real CDistOps3::BruteForceDistance(const CModel3D &model, const VECTOR3 &vQuery) const
{
	using namespace std;
	Real distance;
	double d = std::numeric_limits<double>::max();

	int i;
	
	/* initialize distance with default values */
	distance  = std::numeric_limits<Real>::max();

	VECTOR3 c;
	
	int nNumP = model.NumVertices();

	const VEC3Array& vVertices = model.GetVertices();

	for(int i = 0; i < nNumP; i++)
	{
		VECTOR3 b = vVertices[i];
		c			= VECTOR3::createVector(b,vQuery);
		d			= c.norm2();
		if(distance > d)
		{
			distance = d;
		}//end if
	}//end for

	distance = sqrt(distance);

	return distance;

}//end BruteForceDistance
bool CDistOps3::BruteForceInnerPoints(CModel3D & model, VECTOR3& vQuery)
{

		CRay3f ray3(vQuery,  VECTOR3(1,0,0) );

		CTriangle3f tri3;

		const CAABB3f &rBox = model.GetBox();

		if(!rBox.Inside(vQuery))
			return false;

		int nIntersections = 0;

		for(int i = 0; i < model.GetFaces().size(); i++)
		{
			tri3 = *model.GetTriangle(i);

			CIntersectorRay3Tri3f intersector(ray3, tri3);
			
			if(intersector.Intersection())
				nIntersections++;
		}//end for

		if(nIntersections % 2 == 0)
			return false;
		else
			return true;

}//end BruteForceInnerPoints

Real CDistOps3::SimpleLUB(AABBTree3f &tree, const VECTOR3 &vQuery)
{
	
	double d = std::numeric_limits<double>::max();

	double lowerBound  = std::numeric_limits<double>::max();
	double upperBound  = -std::numeric_limits<double>::max();

	
	list<AABBNode3f*> lBFS;
	list<AABBNode3f*>::iterator lIter;

	for(int i = 0; i < 2; i++)
		lBFS.push_back(tree.GetRoot()->m_Children[i]);
    
	AABBNode3f *nBest = NULL;

	int vSize = (int)lBFS.size();

	int nLeafCount = 0;

	//* loop until there are only leaves in the list */

	while(vSize != nLeafCount)
	{
		
		nLeafCount = 0;
		int j = 0;
		lowerBound = std::numeric_limits<double>::max();
		double *dLowerBounds = new double[vSize];
		/* find best upper bound */
		for(lIter = lBFS.begin(); lIter != lBFS.end(); lIter++)
		{
			AABBNode3f *pNode = *lIter;
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
			
			AABBNode3f *pNode = *lIter;
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

	}//end while

	//* initialize distance with default values */


	///* compute closest point in this circle */

	vector<CTriangle3f*>& vTriangles = nBest->GetTriangles();
	vector<VECTOR3*> &vVerts = nBest->GetVertices();

	int i;

	VECTOR3 vNearest;

	d = DistTriPoint(vVerts, vQuery, vNearest);

	//int nNumP = (int)vVerts.size();

	//for(i = 0; i < nNumP; i++)
	//{

	//	VECTOR3 vVec = *vVerts[i];

	//	Real distance = VECTOR3::createVector(vVec, vQuery).norm2();

	//	if(distance < d)
	//	{
	//		 d = distance;
	//	}//end if

	//}//end for

	for(lIter = lBFS.begin(); lIter != lBFS.end(); lIter++)
	{
		AABBNode3f *pBox = *lIter;
		if(pBox == nBest)
			continue;

		vector<VECTOR3*> &vVertices = pBox->GetVertices();

		Real dist = DistTriPoint(vVertices, vQuery, vNearest);

		if(dist < d)
		{
			d = dist;
		}


	}

	return sqrt(d);

}

Res_t CDistOps3::CoherencyLUB(const AABBTree3f &tree, const VECTOR3 &vQuery, Real rLUB)
{

	Res_t result;

	double d = std::numeric_limits<double>::max();

	double lowerBound  = std::numeric_limits<double>::max();
	double upperBound  = -std::numeric_limits<double>::max();

	list<AABBNode3f*> lBFS;
	list<AABBNode3f*>::iterator lIter;

	for(int i = 0; i < 2; i++)
		lBFS.push_back(tree.GetRoot()->m_Children[i]);
    
	AABBNode3f *nBest = NULL;

	int vSize = (int)lBFS.size();

	int nLeafCount = 0;

	//* loop until there are only leaves in the list */

	while(vSize != nLeafCount)
	{
		
		nLeafCount = 0;
		int j = 0;
		lowerBound = std::numeric_limits<double>::max();
		double *dLowerBounds = new double[vSize];
		/* find best upper bound */
		for(lIter = lBFS.begin(); lIter != lBFS.end(); lIter++)
		{
			AABBNode3f *pNode = *lIter;
			dLowerBounds[j] = pNode->GetLowerBoundSqr(vQuery);
			if(lowerBound > dLowerBounds[j])
			{
				lowerBound = dLowerBounds[j];
				nBest = pNode;
			}//end if
			j++;
		}//end for
		
		/* get upper bound for best element */
		Real bestUpper = nBest->GetUpperBoundSqr(vQuery);
		upperBound = (rLUB > bestUpper) ? bestUpper : rLUB;
		
		lIter = lBFS.begin();
		for(int i = 0; i < vSize; i++)
		{
			
			AABBNode3f *pNode = *lIter;
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

	}//end while

	///* compute closest point in this circle */

	vector<VECTOR3*> &vVerts = nBest->GetVertices();

	VECTOR3 vNearest;

	int i;

	d = DistTriPoint(vVerts, vQuery, vNearest);

	for(lIter = lBFS.begin(); lIter != lBFS.end(); lIter++)
	{
		AABBNode3f *pBox = *lIter;
		if(pBox == nBest)
			continue;

		vector<VECTOR3*> &vVertices = pBox->GetVertices();

		Real dist = DistTriPoint(vVertices, vQuery, vNearest);

		if(dist < d)
		{
			d = dist;
		}//end if

	}//end for

	result.rDistance = d;
	result.vVec      = vNearest;

	return result;
}//end

Real CDistOps3::DistTriPoint(vector<VECTOR3*>& vVerts, const VECTOR3 &vQuery, VECTOR3& vNearest)
{
	int i;

	double d = std::numeric_limits<double>::max();
	
	int nNumP = (int)vVerts.size();

	for(i = 0; i < nNumP; i++)
	{

		VECTOR3 vVec = *vVerts[i];

		Real distance;
		distance = VECTOR3::createVector(vVec, vQuery).norm2();

		if(distance < d)
		{
			 d = distance;
			 vNearest = vVec;
		}//end if

	}//end for

	return d;
}//end DistTriPoint
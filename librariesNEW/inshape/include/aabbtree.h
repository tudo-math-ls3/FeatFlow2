/***************************************************************************
 *   Copyright (C) 2006 by Raphael Mnster   *
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


#if !defined _CAABBNODE_H_
#define _CAABBNODE_H_

//===================================================
//				     Includes
//===================================================

#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>
#include "vector2.h"
#include "aabb2.h"
using namespace std;

//==================================================
//					DEFINES
//==================================================

#define ITEMCOUNT 1

#define MAX_ITEMS 256

#define MAX_NODES 2048

//==================================================
//  		      DECLARATIONS
//==================================================

/* foward declaration */
class CAABBTree;
class DistOps;
class CNurbs;
class COBB;
class CApproxCurve;

typedef std::pair<vector<CApproxCurve*>::const_iterator, vector<CApproxCurve*>::const_iterator> PairIters;

void clearBuckets();

/* class for a node of a AABBTree */
class CAABBNode
{
public:

	/* constructor */
	CAABBNode(const vector<COBB*> &vBoxes, vector<COBB*> *vBucket ,int iBucket);
	CAABBNode(){};
	CAABBNode(CAABBNode *pNode, CAABBTree *pTree);
	/* deconstructor */
	~CAABBNode();

	/* member functions */
	inline int	   GetItemCount(){return (int)m_vBucket[m_iBucket].size();};
 	inline CAABB&  GetBoundingBox(){return m_Box;};
	inline bool    IsLeaf(){return (m_Children[0] == NULL);};
	inline int	   GetAxis(){return (int)m_iSplitAxis;};
	
	vector<COBB*>& Boxes(){return m_vBucket[m_iBucket];};
	inline short   Bucket(){return m_iBucket;};
	
	double		   GetUpperBound(const VECTOR2 &vQuery);
	double		   GetLowerBound(const VECTOR2 &vQuery);
    Real Area();

	Real		   Xmin(){return m_Box.m_Vertices[0].x;};
	Real		   Xmax(){return m_Box.m_Vertices[1].x;};
	Real		   Ymin(){return m_Box.m_Vertices[0].y;};
	Real		   Ymax(){return m_Box.m_Vertices[1].y;};
	
	

	/* public member variables */
	CAABBNode	  *m_Children[2];
	
	/* the lower bound is the nearest curve point on the hull box of the curve */
private:
	/* private member variables */
	CAABB	  m_Box;
	short	  m_iBucket;
	bool	  m_iSplitAxis;
	vector<COBB*> *m_vBucket;
	
	
	
	
	


   	 /* make the tree a friend */
	friend class CAABBTree;

	friend ostream& operator<<(ostream& out, CAABBNode *tNode);

};

class CAABBTree
{
public:
	/* constructor */
	CAABBTree(); 
	CAABBTree(const vector<CNurbs*> &curves);
	CAABBTree(const vector<CApproxCurve*> &curves, int iMethod = 0);
	CAABBTree(CAABBNode *tNode);
	CAABBTree(CAABBTree *pTree);



	/* deconstructor */
	~CAABBTree();

	/* inline member functions */
	inline CAABBNode* GetRoot(){return m_Root;};
	inline CAABBNode* getCurrent(){return m_Current;};
	inline void ClearArea(){m_dArea = 0.0;};
	
	/* member functions */

	void BuildTree(PairIters &pIters, int iMethod = 0);
	void BuildTree(CAABBTree &pTree);

	void       subdivideNode(CAABBNode *pNode, int nBucket);
	void       subdivideNodeApprox(CAABBNode *pNode);
	void	   MedianCut(CAABBNode *pNode, int nBucket);
	void	   MinArea(CAABBNode *pNode, int iBucket);
	void       DeleteSubTree(CAABBNode *pNode);
	void	   copySubTree(CAABBNode *pNode);
	void	   DestroyTree();
	void	   inOrder(CAABBNode *pNode);
	void	   breadthFirstSearch();
	void 	   DepthFirstSearchStack();
	int 	   DepthFirstSearch(CAABBNode *pNode,int iDepth);
	void	   breadthFirstSearchDistance(const VECTOR2 &vQuery, list<CAABBNode*> &lBFS);
	void	   breadthFirstSearchDistance(const VECTOR2 &vQuery, list<CAABBNode*> &lBFS, double dLBoundSuccessor);
	void	   updateBBoxes(CAABBNode *pNode);
	bool       IsInitialized(){return m_bInit;};
	CAABBNode* find(const VECTOR2 &vQuery);
	CAABBNode* findNearestBV(const VECTOR2 &vQuery);
	

	enum
	{
		SPATIAL_MEDIAN,
		MEDIAN_CUT,
		MIN_AREA
	};

	vector<COBB*> Buckets[MAX_NODES];
	int m_nNodes;
	double m_dArea;

private:
	/* private member variables */
	/* pointer to the tree root */
	CAABBNode *m_Root;
	/* pointer to the current node when moving around */
	CAABBNode *m_Current;
	int        m_Count;
	bool       m_bInit;

	



};

#endif






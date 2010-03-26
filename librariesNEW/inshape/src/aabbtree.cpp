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

#include <queue>
#include <list>
#include <stack>
#include <aabbtree.h>
#include <distops.h>
#include <approxcurve.h>
#include <obb2.h>

//vector<COBB*> Buckets[512];


void clearBuckets()
{
// 	for(int i = 0; i < 512; i++)
// 	{
// 		Buckets[i].clear();
// 	}
}

///////////////////////////////////////////
//
//	overloaded output operator
//
//
///////////////////////////////////////////

ostream& operator<<(ostream& out, CAABBNode *b1) 
{
	vector<COBB*>::iterator iter;
	for(iter = b1->Boxes().begin(); iter != b1->Boxes().end(); iter++)
	{
		COBB *box = *iter;
		out << box <<endl;
	}

	return out;
}



///////////////////////////////////////////
//
//	class CAABBNode :
//
//	a node in an AABBTree
//
//	stores bounding boxes and the
//
//  primitives they contain
//
//	
//
///////////////////////////////////////////

/* constructor */
CAABBNode::CAABBNode(const vector<COBB*> &vBoxes, vector<COBB*> *vBucket ,int iBucket) 
{

	m_Box.InitBox(vBoxes);
	//m_Boxes		 = vBoxes;

	m_iSplitAxis = 0;
	
	m_Children[0] = NULL;

	m_Children[1] = NULL;

	m_iBucket = static_cast<short>(iBucket);

	m_vBucket = vBucket;

}//end constructor


CAABBNode::CAABBNode(CAABBNode *pNode, CAABBTree *pTree) : m_Box(pNode->GetBoundingBox().m_Vertices[0], pNode->GetBoundingBox().m_Vertices[1])
{

	

	if(pNode->GetAxis())
		m_iSplitAxis = 1;
	else
		m_iSplitAxis = 0;

	m_Children[0] = NULL;

	m_Children[1] = NULL;

	m_iBucket = pNode->Bucket();
	
	m_vBucket = pTree->Buckets;

	pTree->Buckets[m_iBucket] = pNode->Boxes();

}//end constructor

CAABBNode::~CAABBNode()
{
	if(m_Children[0])
	{
		delete m_Children[0];
		m_Children[0] = NULL;

		delete m_Children[1];
		m_Children[1] = NULL;
	}

}//end deconstructor

/////////////////////////////////////////////
//
//	The function returns the upper bound
//	on the distance between the query point
//	and the bounding box
//
//	Input:
//		- the query point
//
//	Output:
//		- the upper bound of the distance
//
/////////////////////////////////////////////
double CAABBNode::GetUpperBound(const VECTOR2 &vQuery)
{

	///* get bounding box center */
	double x = (m_Box.m_Vertices[1].x - m_Box.m_Vertices[0].x) * 0.5 + m_Box.m_Vertices[0].x;
	double y = (m_Box.m_Vertices[1].y - m_Box.m_Vertices[0].y) * 0.5 + m_Box.m_Vertices[0].y;

	VECTOR2 vCenter(x,y); 

	if( (vQuery.x >= vCenter.x) && (vQuery.y >= vCenter.y) )
	{
		VECTOR2 c = VECTOR2::createVector(vQuery, m_Box.m_Vertices[0]);
		return c.mag();
	}
	else if( (vQuery.x >= vCenter.x) && (vQuery.y < vCenter.y) )
	{
		VECTOR2 c = VECTOR2::createVector(vQuery, VECTOR2(m_Box.m_Vertices[0].x, m_Box.m_Vertices[1].y));
		return c.mag();
	}
	else if( (vQuery.x < vCenter.x) && (vQuery.y >= vCenter.y) )
	{
		VECTOR2 c = VECTOR2::createVector(vQuery, VECTOR2(m_Box.m_Vertices[1].x, m_Box.m_Vertices[0].y));
		return c.mag();
	}
	else
	{
		VECTOR2 c = VECTOR2::createVector(vQuery, m_Box.m_Vertices[1]);
		return c.mag();		
	}

	return -1.0;

}//end GetLowerBound

/////////////////////////////////////////
//
//	The function returns the upper bound
//	on the distance between the box and
//	the query point.
//
//
/////////////////////////////////////////
double CAABBNode::GetLowerBound(const VECTOR2 &vQuery)
{

	VECTOR2 vSol;

	if(vQuery.x < Xmin())
		vSol.x = Xmin()-vQuery.x;
	else if(vQuery.x > Xmax())
		vSol.x = vQuery.x - Xmax();
	else
		vSol.x = 0.0;

	if(vQuery.y < Ymin())
		vSol.y = Ymin()-vQuery.y;
	else if(vQuery.y > Ymax())
		vSol.y = vQuery.y - Ymax();
	else
		vSol.y = 0.0;

	return vSol.mag();

}//end GetLowerBound


//////////////////////////////////////
//
//	class AABBTree
//
//	control class for a hierarchy
//
//	of axis-aligned bounding boxes
//
//////////////////////////////////////

CAABBTree::CAABBTree(CAABBNode *tNode)
{
	m_Root    = tNode;
	m_Current = tNode;
}//end constructor

CAABBTree::CAABBTree(CAABBTree *pTree)
{

	m_Root	= new CAABBNode(pTree->GetRoot(), this);
	
	/* a queue stores the nodes to be subdivided */
	queue<CAABBNode*> qNodesOrig;
																			
	/* initially only the root is in the queue */
	qNodesOrig.push(pTree->GetRoot());
	
	/* a queue stores the nodes to be subdivided */
	queue<CAABBNode*> qNodes;
																			
	/* initially only the root is in the queue */
	qNodes.push(m_Root);
	m_nNodes = pTree->m_nNodes;
	/* build the bounding box hierarchy */
	while(!qNodesOrig.empty())
	{

		int num = (int)qNodes.size();
		/* get front element */
		CAABBNode *origNode = qNodesOrig.front();

		/* remove front element */
		qNodesOrig.pop();
		
		CAABBNode *pNode = qNodes.front();
		
		qNodes.pop();
		
		if(origNode->IsLeaf())
		{
			/* NODE IS A LEAF DO NOTHING */
		
		}//end if
		else
		{
			/* copy nodes */
			pNode->m_Children[0] = new CAABBNode(origNode->m_Children[0], this);
			pNode->m_Children[1] = new CAABBNode(origNode->m_Children[1], this);
		
			/* insert original nodes */
			qNodesOrig.push(origNode->m_Children[0]);
			qNodesOrig.push(origNode->m_Children[1]);
			
			/* insert nodes of the copy */
			qNodes.push(pNode->m_Children[0]);
			qNodes.push(pNode->m_Children[1]);
		

		}//end else
			
	}//end while
	
}//end constructor

void CAABBTree::copySubTree( CAABBNode *pNode )
{

	if(pNode != NULL)
	{
		
	}//end if

}//end copSubTree

/////////////////////////////////////////////////////////
//
//	class AABBTree
//
//	constructor:
//		- Input:
//			-curves: a vector of curves, the tree will be built
//					based on the curves bounding boxes
//			-iMethod: built strategy for the tree, possible values
//					  are SPATIAL_MEDIAN and MEDIAN_CUT, 
//
/////////////////////////////////////////////////////////

CAABBTree::CAABBTree(const vector<CApproxCurve*> &curves, int iMethod)
{
	/* stores the boxes */
	vector<COBB*> vBoxes;

	int iBucket =0;

	/* get the boxes */
	vector<CApproxCurve*>::const_iterator iter;
	cout<<"#Curves: "<<curves.size()<<endl;
	for(iter = curves.begin(); iter != curves.end(); iter++)
	{
		CApproxCurve* curve = *iter;
		COBB *aBox = curve->m_bBox;
		vBoxes.push_back(aBox);
		Buckets[iBucket].push_back(aBox);
	}//end for

	/* create the root node */
	m_Root	  = new CAABBNode(vBoxes, Buckets ,iBucket);
	m_Current = m_Root;

	iBucket++;

	/* a queue stores the nodes to be subdivided */
	queue<CAABBNode*> qNodes;
																			
	/* initially only the root is in the queue */
	qNodes.push(m_Root);

	int nNodes = 1;

	/* build the bounding box hierarchy */
	while(!qNodes.empty())
	{

		int num = (int)qNodes.size();
		/* get front element */
		CAABBNode *aNode = qNodes.front();
		m_dArea+= aNode->Area();
		/* remove front element */
		qNodes.pop();
		
		/* target item count already reached ? */
		if(aNode->GetItemCount() > ITEMCOUNT)
		{
			/* more than one item, subdivide node */
			if(iMethod == SPATIAL_MEDIAN)
			{
				subdivideNode(aNode,iBucket);
				iBucket+=2;
				
			}
			else if(iMethod == MEDIAN_CUT)
			{
				MedianCut(aNode,iBucket);
				iBucket+=2;
			}
			else if(iMethod == MIN_AREA)
			{
				MinArea(aNode, iBucket);
				iBucket+=2;
			}
			else
			{
			}

			/* continue subdivision only if both nodes still contain items */
			if( (aNode->m_Children[0]->GetItemCount() > 0) && (aNode->m_Children[1]->GetItemCount() > 0) )
			{
				/* insert children into queue */
				qNodes.push(aNode->m_Children[0]);
				qNodes.push(aNode->m_Children[1]);
				nNodes+=2;
			}//end if

		}//end if
	
	}//end while
	
	m_nNodes = nNodes;
	cout <<"AABBTree built, number of nodes: "<< nNodes<<endl;

}//end constructor

void CAABBTree::BuildTree(PairIters &pIters, int iMethod)
{

	/* stores the boxes */
	vector<COBB*> vBoxes;

	int iBucket =0;

	/* get the boxes */
	for(;pIters.first != pIters.second; pIters.first++)
	{
		CApproxCurve* curve = *pIters.first;
		COBB *aBox = curve->m_bBox;
		vBoxes.push_back(aBox);
		Buckets[iBucket].push_back(aBox);
	}//end for

	/* create the root node */
	m_Root	  = new CAABBNode(vBoxes, Buckets ,iBucket);
	m_Current = m_Root;

	iBucket++;

	/* a queue stores the nodes to be subdivided */
	queue<CAABBNode*> qNodes;
																			
	/* initially only the root is in the queue */
	qNodes.push(m_Root);

	int nNodes = 1;

	/* build the bounding box hierarchy */
	while(!qNodes.empty())
	{

		int num = (int)qNodes.size();
		/* get front element */
		CAABBNode *aNode = qNodes.front();
		m_dArea+= aNode->Area();
		/* remove front element */
		qNodes.pop();
		
		/* target item count already reached ? */
		if(aNode->GetItemCount() > ITEMCOUNT)
		{
			/* more than one item, subdivide node */
			if(iMethod == SPATIAL_MEDIAN)
			{
				subdivideNode(aNode,iBucket);
				iBucket+=2;
				
			}
			else if(iMethod == MEDIAN_CUT)
			{
				MedianCut(aNode,iBucket);
				iBucket+=2;
			}
			else if(iMethod == MIN_AREA)
			{
				MinArea(aNode, iBucket);
				iBucket+=2;
			}
			else
			{
			}

			/* continue subdivision only if both nodes still contain items */
			if( (aNode->m_Children[0]->GetItemCount() > 0) && (aNode->m_Children[1]->GetItemCount() > 0) )
			{
				/* insert children into queue */
				qNodes.push(aNode->m_Children[0]);
				qNodes.push(aNode->m_Children[1]);
				nNodes+=2;
			}//end if

		}//end if
	
	}//end while
	
	m_nNodes = nNodes;
	m_bInit  = true;
	cout <<"AABBTree built, number of nodes: "<< nNodes<<endl;

}//end BuildTree

CAABBTree::CAABBTree(const vector<CNurbs*> &curves)
{
	
	ClearArea();
	
	/* stores the boxes */
	vector<COBB*> vBoxes;

	/* get the boxes */
	vector<CNurbs*>::const_iterator iter;
	for(iter = curves.begin(); iter != curves.end(); iter++)
	{
		CNurbs* curve = *iter;
		COBB *aBox = curve->BoundingBox();
		vBoxes.push_back(aBox);
	}//end for

	/* create the root node */
	m_Root	  = new CAABBNode(vBoxes, Buckets,0);
    m_Current = m_Root;

	/* a queue stores the nodes to be subdivided */
	queue<CAABBNode*> qNodes;

	
	
	/* initially only the root is in the queue */
	qNodes.push(m_Root);
	/* build the bounding box hierarchy */
	while(!qNodes.empty())
	{
		/* get front element */
		CAABBNode *aNode = qNodes.front();

		/* remove front element */
		qNodes.pop();
		
		/* target item count already reached ? */
		if(aNode->GetItemCount() > 1)
		{
			/* more than one item, subdivide node */
			subdivideNode(aNode,1);
			
			/* stop subdivision if all items were assigned to one child */
			if( (aNode->m_Children[0]->GetItemCount() > 0) && (aNode->m_Children[1]->GetItemCount() > 0) )
			{
				/* insert children into queue */
				qNodes.push(aNode->m_Children[0]);
				qNodes.push(aNode->m_Children[1]);

			}//end if

		}//end if

	}//end while
	
}//end constructor

CAABBTree::CAABBTree() : m_bInit(false)
{
	m_Root    = NULL;
	m_Current = NULL;
}//end constructor

CAABBTree::~CAABBTree()
{
	if(m_bInit)
	{
		DestroyTree();
	}//end if
}//end deconstructor


//TODO: check if m_Root is really = NULL at the end
void CAABBTree::DestroyTree()
{

	if(m_Root != NULL)
	{

		DeleteSubTree(m_Root);
		
	}//end if

	

	for(int i = 0; i < MAX_NODES; i++)
	{
		Buckets[i].clear();
	}//end for

	m_bInit = false;

	if(GetRoot() != NULL)
		cout<<"Root not equal to zero"<<endl;

}//end destroyTree

void CAABBTree::DeleteSubTree(CAABBNode *pNode)
{

	if(pNode != NULL)
	{
		DeleteSubTree( pNode->m_Children[0] );
		pNode->m_Children[0] = NULL;
		DeleteSubTree( pNode->m_Children[1] );
		pNode->m_Children[1] = NULL;
		delete pNode;
		pNode = NULL;
	}//end if

}//end DeleteSubTree

//CAABBNode* CAABBTree::find(const VECTOR2 &vQuery)
//{
//	CAABBNode *t = m_Root;
//
//	int splitAxis;
//	/* while the node pointer is valid */
//	while(!t->IsLeaf())
//	{
//		splitAxis = t->m_iSplitAxis;
//		
//		if(splitAxis == CAABB::xAxis)
//		{
//			if(vQuery.x < t->m_dSplit)
//			{
//				cout<<"choose left tree: "<<endl;
//				cout<<"lower bound left: "<<t->m_Children[0]->GetUpperBound(vQuery)<<endl;
//				cout<<"upper bound right: "<<t->m_Children[1]->GetLowerBound(vQuery)<<endl;
//				t = t->m_Children[0];
//			}
//			else
//			{
//				cout<< "choose right branch" <<endl;
//				cout<<"upper bound left: "<<t->m_Children[0]->GetLowerBound(vQuery)<<endl;
//				cout<<"lower bound right: "<<t->m_Children[1]->GetUpperBound(vQuery)<<endl;
//				t = t->m_Children[1];
//			}
//		}//end if
//		else
//		{
//			if(vQuery.y < t->m_dSplit)
//			{
//				cout<<"choose left tree: "<<endl;
//				cout<<"lower bound left: "<<t->m_Children[0]->GetUpperBound(vQuery)<<endl;
//				cout<<"upper bound right: "<<t->m_Children[1]->GetLowerBound(vQuery)<<endl;
//				t = t->m_Children[0];
//			}
//			else
//			{
//				cout<< "choose right branch" <<endl;
//				cout<<"upper bound left: "<<t->m_Children[0]->GetLowerBound(vQuery)<<endl;
//				cout<<"lower bound right: "<<t->m_Children[1]->GetUpperBound(vQuery)<<endl;
//				t = t->m_Children[1];
//			}
//		}//end else
//
//
//	}//end while
//
//	return t;
//	
//}//end find

CAABBNode* CAABBTree::findNearestBV(const VECTOR2 &vQuery)
{
	CAABBNode *t = m_Root;

	int splitAxis = 0;
	/* while the node pointer is valid */
	while(!t->IsLeaf())
	{
	
	
		cout<<"lower bound left: "<<t->m_Children[0]->GetUpperBound(vQuery)<<endl;
		cout<<"upper bound right: "<<t->m_Children[1]->GetLowerBound(vQuery)<<endl;


	}//end while

	return t;
	
}//end find

//////////////////////////////////////////////
//
//	The function sudivides a node in an AABBTree
//	using the median split strategy and
//	assigns the parents items to the children
//
//////////////////////////////////////////////

void CAABBTree::subdivideNode(CAABBNode *pNode, int nBucket)
{

	/* get the nodes bounding box */
	CAABB &bAABB = pNode->GetBoundingBox();

	/* get the longest axis the bounding volume will be split there */
	bool iAxis    = bAABB.GetLongestAxis();

	/* get the center of the bounding volume */
	VECTOR2 vCenter = bAABB.GetCenter();

	vector<COBB*>::iterator bIter;

	double minDist = 1.7E+308;

	/* here the bounding volume will be split at the longest axis */
	VECTOR2 vSplit;
	
	/* Spatial Median split strategy: */
	/* the bounding box will be split at the OBB vertex that is closest to the bounding box center */
	/* loop over all the items associated with the node to find this vertex */
	for(bIter = pNode->Boxes().begin(); bIter != pNode->Boxes().end(); bIter++)
	{
		
		COBB *box = *bIter;

		/* loop over all box vertices to */
		/* find the box vertex closest to the center of the bounding volume */
		for(int i = 0; i < 4; i++)
		{
			VECTOR2 c = VECTOR2::createVector(box->m_Vertices[i], vCenter);
			double d = c.mag();
			if(minDist > d)
			{
				minDist = d;
				vSplit = box->m_Vertices[i];
			}//end if

		}//end for

	}//end for

	/* assign the parents items to the children */
	vector<COBB*> nBoxes[2];
	
	/* split the items into two buckets relative to the split axis */
	for(bIter = pNode->Boxes().begin(); bIter != pNode->Boxes().end(); bIter++)
	{
		COBB *hBox = *bIter;
	
		/* if longest axis = xAxis */
		if(iAxis == CAABB::xAxis)
		{
			/* assign split axis */
			pNode->m_iSplitAxis = iAxis;

			/* value at that the bounding volume is split along the split axis */
			if(hBox->m_Vertices[0].x < vSplit.x)
				Buckets[nBucket].push_back(hBox);
			else
				Buckets[nBucket+1].push_back(hBox);
		}//end if
		/* if longest axis = yAxis */
		else if(iAxis == CAABB::yAxis)
		{
			/* assign split axis */
			pNode->m_iSplitAxis = iAxis;

			/* value at that the bounding volume is split along the split axis */
			if(hBox->m_Vertices[0].y < vSplit.y)
				Buckets[nBucket].push_back(hBox);
			else
				Buckets[nBucket+1].push_back(hBox);			
		}//end if
		else
		{

		}//end else


	}//end for

	/* create the new children and assign the items to the children */
	pNode->m_Children[0] = new CAABBNode(Buckets[nBucket],Buckets, nBucket);
	pNode->m_Children[1] = new CAABBNode(Buckets[nBucket+1],Buckets, nBucket+1);

}//end subdivideNode

void CAABBTree::MedianCut(CAABBNode *pNode, int nBucket)
{

	/* get the nodes bounding box */
	CAABB &bAABB = pNode->GetBoundingBox();

	/* get the longest axis the bounding volume will be split there */
	bool iAxis    = bAABB.GetLongestAxis();

	/* get the center of the bounding volume */
	VECTOR2 vCenter = bAABB.GetCenter();

	vector<COBB*>::iterator bIter;

	double minDist = 1.7E+308;
		
	VECTOR2 vMedian(0,0);
	
	/* Median Cut strategy: */
	for(bIter = pNode->Boxes().begin(); bIter != pNode->Boxes().end(); bIter++)
	{
		
		COBB *box = *bIter;

		VECTOR2 vBoxCenter = box->getCenter();

		vMedian = vMedian + vBoxCenter;

	}//end for

	double nSize = (double)pNode->Boxes().size();

	if(nSize > 0)
		vMedian /= nSize;

	/* assign the parents items to the children */
	int s = (int) pNode->Boxes().size();
	/* split the items into two buckets relative to the split axis */
	for(bIter = pNode->Boxes().begin(); bIter != pNode->Boxes().end(); bIter++)
	{
		COBB *hBox = *bIter;
	
		/* if longest axis = xAxis */
		if(iAxis == CAABB::xAxis)
		{
			/* assign split axis */
			pNode->m_iSplitAxis = iAxis;

			/* value at that the bounding volume is split along the split axis */
			if(hBox->m_Vertices[0].x < vMedian.x)
			{
				Buckets[nBucket].push_back(hBox);
			}
			else
			{
				Buckets[nBucket+1].push_back(hBox);
			}
		}//end if
		/* if longest axis = yAxis */
		else if(iAxis == CAABB::yAxis)
		{
			/* assign split axis */
			pNode->m_iSplitAxis = iAxis;

			/* value at that the bounding volume is split along the split axis */
			

			if(hBox->m_Vertices[0].y < vMedian.y)
			{
				Buckets[nBucket].push_back(hBox);
			}
			else
			{
				Buckets[nBucket+1].push_back(hBox);
			}
		}//end if
		else
		{

		}//end else


	}//end for

	/* create the new children and assign the items to the children */
	pNode->m_Children[0] = new CAABBNode(Buckets[nBucket], Buckets ,nBucket);
	pNode->m_Children[1] = new CAABBNode(Buckets[nBucket+1], Buckets ,nBucket+1);

}//end MedianCut

void CAABBTree::inOrder(CAABBNode *pNode)
{

	if(pNode != NULL)
	{
		inOrder(pNode->m_Children[0]);
		
		inOrder(pNode->m_Children[1]);
	}

}//end inOrder

void CAABBTree::breadthFirstSearch()
{

	queue<CAABBNode*> qBFS;

	qBFS.push(m_Root);

	cout<<"Begin breadth first search"<<endl;

	int iNodes = 0;

	while(!qBFS.empty())
	{
		CAABBNode *pNode = qBFS.front();
		iNodes++;
		qBFS.pop();
		if(!pNode->IsLeaf())
		{
			qBFS.push(pNode->m_Children[0]);
			qBFS.push(pNode->m_Children[1]);
		}//end if

	}//end while

	cout<<"Number of nodes in AABBTree: "<<iNodes<<endl;

}//end breadthFirstSearch

void CAABBTree::DepthFirstSearchStack()
{

	stack<CAABBNode*> sDFS;

	sDFS.push(m_Root);

	while(!sDFS.empty())
	{
		CAABBNode *pNode = sDFS.top();
		sDFS.pop();
		if(!pNode->IsLeaf())
		{
			sDFS.push(pNode->m_Children[0]);
			sDFS.push(pNode->m_Children[1]);
		}//end if

	}//end while

}//end DepthFirstSearchStack

int CAABBTree::DepthFirstSearch(CAABBNode *pNode,int iDepth)
{

	int dR = 0;
	int dL = 0;
	int iReturn = 0;

	if(!pNode->IsLeaf())
	{
		dL = DepthFirstSearch(pNode->m_Children[0],iDepth+1);
		dR = DepthFirstSearch(pNode->m_Children[1],iDepth+1);
		iReturn = ( dL >= dR ) ? dL : dR;
		return iReturn;
	}//end if
	else
	{
		iReturn = iDepth;
		return iReturn;
	}

}//end DepthFirstSearch

///////////////////////////////////////
//
//	A function that computes the
//	nearest leaf nodes to a grid point.
//
//  Input:
//		- the grid point
//		- a list to filled with the 
//		  nearest nodes
//	Makes use of coherency
//
//
///////////////////////////////////////


void CAABBTree::breadthFirstSearchDistance(const VECTOR2 &vQuery, list<CAABBNode*> &lBFS, double dLBoundSuccessor)
{
	list<CAABBNode*>::iterator lIter;
	double upperBound  = 0.0;
	double lowerBound  = 1.7E+308;

	/* enqueue children of root node */
	CAABBNode *c1 = m_Root->m_Children[0];
	CAABBNode *c2 = m_Root->m_Children[1];

	if(c1)
		lBFS.push_back(c1);
	if(c2)
		lBFS.push_back(c2);

	/* stores the the nearest node */
	CAABBNode *nBest = NULL;

	/* current size of queue */
	int vSize = (int)lBFS.size();
	int leafCount = 0;
	int j;
	do
	{
		leafCount = 0;
		lowerBound = 1.7E+308;
		j = 0;
		double *dLBounds = new double[vSize];
		/* find best upper bound */
		/* maybe try best lower bound or a combination ?? */
		for(lIter = lBFS.begin(); lIter != lBFS.end();lIter++)
		{
			CAABBNode *pNode = *lIter;
			dLBounds[j] = pNode->GetLowerBound(vQuery);
			if(lowerBound > dLBounds[j])
			{
				lowerBound = dLBounds[j];
				nBest = pNode;
			}//end if
			j++;
		}//end for
		
		/* get best nodes upper bound value */
		double bestBound = nBest->GetUpperBound(vQuery);
		upperBound = (dLBoundSuccessor > bestBound ) ? bestBound  : dLBoundSuccessor;
		j = 0;
		/* compare best nodes upper bound to the lower bound of all other nodes in the queue */
		/* if the upper bound is better the subtree can be pruned */
		for(lIter = lBFS.begin(); lIter != lBFS.end();)
		{
			
			CAABBNode *pNode = *lIter;
			
			/* the best node's children make in to the next round */
			if(pNode == nBest)
			{
				if(!pNode->IsLeaf())
				{
					lBFS.push_front(pNode->m_Children[0]);
					lBFS.push_front(pNode->m_Children[1]);
					/* erase parent and go on to the next node */
					lIter = lBFS.erase(lIter);
					j++;
					continue;
					
				}
				else
				{
					leafCount++;
				}
			}//end if
			else if(upperBound > dLBounds[j] || pNode->GetBoundingBox().Inside(vQuery) )
			{
				if(!pNode->IsLeaf())
				{
					lBFS.push_front(pNode->m_Children[0]);
					lBFS.push_front(pNode->m_Children[1]);
					/* erase parent and go on to the next node */
					lIter = lBFS.erase(lIter);
					j++;
					continue;
					
				}//end if
				else
				{
					leafCount++;
				}
			}//end if
			else
			{
				lIter = lBFS.erase(lIter);
				j++;
				continue;
			}//end else
			lIter++;
			j++;
		}//end for
		vSize = (int)lBFS.size();
		delete[] dLBounds;
		/* while not all nodes are leaves */
	}while(vSize != leafCount); //end while
	
}//end breadthFirstSearchDistance


///////////////////////////////////////
//
//	A function that computes the
//	nearest leaf nodes to a grid point.
//
//  Input:
//		- the grid point
//		- a list to filled with the 
//		  nearest nodes
//
///////////////////////////////////////

void CAABBTree::breadthFirstSearchDistance(const VECTOR2 &vQuery, list<CAABBNode*> &lBFS)
{

	list<CAABBNode*>::iterator lIter;
	double upperBound  = 0.0;
	double lowerBound  = 1.7E+308;

	/* enqueue children of root node */
	CAABBNode *c1 = m_Root->m_Children[0];
	CAABBNode *c2 = m_Root->m_Children[1];

	if(c1)
		lBFS.push_back(c1);
	if(c2)
		lBFS.push_back(c2);

	/* stores the the nearest node */
	CAABBNode *nBest = NULL;

	/* current size of queue */
	int vSize = (int)lBFS.size();
	int leafCount = 0;
	int j;
	do
	{
		leafCount = 0;
		lowerBound = 1.7E+308;
		j = 0;
		double *dLBounds = new double[vSize];
		/* find best upper bound */
		/* maybe try best lower bound or a combination ?? */
		for(lIter = lBFS.begin(); lIter != lBFS.end();lIter++)
		{
			CAABBNode *pNode = *lIter;
			dLBounds[j] = pNode->GetLowerBound(vQuery);
			if(lowerBound > dLBounds[j])
			{
				lowerBound = dLBounds[j];
				nBest = pNode;
			}//end if
			j++;
		}//end for
		
		/* get best nodes upper bound value */
		upperBound = nBest->GetUpperBound(vQuery);
		j = 0;
		/* compare best nodes upper bound to the lower bound of all other nodes in the queue */
		/* if the upper bound is better the subtree can be pruned */
		for(lIter = lBFS.begin(); lIter != lBFS.end();)
		{
			
			CAABBNode *pNode = *lIter;
			
			/* the best node's children make in to the next round */
			if(pNode == nBest)
			{
				if(!pNode->IsLeaf())
				{
					lBFS.push_front(pNode->m_Children[0]);
					lBFS.push_front(pNode->m_Children[1]);
					/* erase parent and go on to the next node */
					lIter = lBFS.erase(lIter);
					j++;
					continue;
					
				}
				else
				{
					leafCount++;
				}
			}//end if
			else if(upperBound > dLBounds[j])
			{
				if(!pNode->IsLeaf())
				{
					lBFS.push_front(pNode->m_Children[0]);
					lBFS.push_front(pNode->m_Children[1]);
					/* erase parent and go on to the next node */
					lIter = lBFS.erase(lIter);
					j++;
					continue;
					
				}//end if
				else
				{
					leafCount++;
				}
			}//end if
			else
			{
				lIter = lBFS.erase(lIter);
				j++;
				continue;
			}//end else
			lIter++;
			j++;
		}//end for
		vSize = (int)lBFS.size();
		delete[] dLBounds;
	/* while not all nodes are leaves */
	}while(vSize != leafCount); //end while

}//end breadthFirstSearchDistance

////////////////////////////////////////
//
//	updates the bounding box hierarchy
//	in a bottom up manner after the
//	items have moved
//
////////////////////////////////////////
void CAABBTree::updateBBoxes(CAABBNode *pNode)
{
	m_dArea+=pNode->Area();
	/* if the node is a leaf */
	/* the new bounds are the bounds of the triangles */
	if(pNode->IsLeaf())
	{

		//double dItemCount = static_cast<double>(pNode->GetItemCount());
		//double dArea   = pNode->GetBoundingBox()->Area();
		//double SAH     = dArea/dItemCount; 
		//cout<<"SAH: "<<SAH<<endl;

		pNode->m_Box.CalculateBounds(pNode->Boxes());
		int nItems = pNode->GetItemCount();
		if(nItems == 1)
		{
			COBB* box = pNode->Boxes().front();
			CApproxCurve *pCurve = box->m_paCurve;
			pNode->m_Box.m_Vertices[0] = pCurve->m_bBox->m_Vertices[0];
			pNode->m_Box.m_Vertices[1] = pCurve->m_bBox->m_Vertices[2];
		}//end if
		else if(nItems > 1)
		{

			CAABB *nextBox = &pNode->m_Box;
			for(int i = 0; i < nItems; i++)
			{
				
				CApproxCurve *c0 = pNode->Boxes()[i]->m_paCurve;
				COBB* box0       = c0->m_bBox;

				if(box0->m_Vertices[0].x < nextBox->m_Vertices[0].x)
					nextBox->m_Vertices[0].x = box0->m_Vertices[0].x;

				if(box0->m_Vertices[2].x > nextBox->m_Vertices[1].x)
					nextBox->m_Vertices[1].x = box0->m_Vertices[2].x;

				if(box0->m_Vertices[0].y < nextBox->m_Vertices[0].y)
					nextBox->m_Vertices[0].y = box0->m_Vertices[0].y;

				if(box0->m_Vertices[2].y > nextBox->m_Vertices[2].y)
					nextBox->m_Vertices[2].y = box0->m_Vertices[2].y;

			}//end for

		}//end else
	
	}//end if

	/* if the node is not a leaf */
	/* the bounding box is the bounding box of the two children */
	else
	{
		//double dItemCount = static_cast<double>(pNode->GetItemCount());
		//double dArea   = pNode->GetBoundingBox()->Area();
		//double SAH     = dArea/dItemCount; 
		//cout<<"SAH: "<<SAH<<endl;
		updateBBoxes(pNode->m_Children[0]);
		updateBBoxes(pNode->m_Children[1]);
		pNode->m_Box.CalculateBounds(pNode->m_Children[0]->m_Box, pNode->m_Children[1]->m_Box);
	}//end else

}//end updateBBoxes


/*!
    \fn CAABBNode::Area()
 */
Real CAABBNode::Area()
{
	return m_Box.Area();
}

void CAABBTree::MinArea(CAABBNode *pNode, int iBucket)
{

	int nBucket = iBucket;

	/* get the nodes bounding box */
	CAABB &bAABB = pNode->GetBoundingBox();

	/* get the longest axis the bounding volume will be split there */
	bool iAxis    = bAABB.GetLongestAxis();

	/* get the center of the bounding volume */
	VECTOR2 vCenter = bAABB.GetCenter();

	vector<COBB*>::iterator bIter;

	double minDist = 1.7E+308;

	/* here the bounding volume will be split at the longest axis */
	VECTOR2 vSplit;
	
	/*MinArea strategy:
	* The bounds of the AABBs of the curves are the potential dividing lines
	* test all theses lines in both axes and select the one that minimizes 
	* the sum of the area of the two child nodes
	*/

	/*SELECTION PHASE:
	* select a dividing line
	*/

	int bestAxis   = -1;
	Real bestPos   = -1.0;
	double MArea   = 10000;

	for(int axis = 0; axis < 2; axis++)
	{
		for(bIter = pNode->Boxes().begin(); bIter != pNode->Boxes().end(); bIter++)
		{
			
			COBB *box = *bIter;
			
			//dividing line candidate
			box->m_Vertices[0].x;

			vector<COBB*> boxes1;
			vector<COBB*> boxes2;
			
			//test line
			for(int j = 0; j < pNode->Boxes().size(); j++)
			{

				COBB *hBox = pNode->Boxes()[j];
					if(axis == CAABB::xAxis)
					{
						/* value at that the bounding volume is split along the split axis */
						if(hBox->getCenter().x < box->m_Vertices[0].x)
							boxes1.push_back(hBox);
						else
							boxes2.push_back(hBox);
					}//end if
					/* if longest axis = yAxis */
					else if(axis == CAABB::yAxis)
					{
						/* assign split axis */
						pNode->m_iSplitAxis = iAxis;

						/* value at that the bounding volume is split along the split axis */
						if(hBox->getCenter().y < box->m_Vertices[0].y)
							boxes1.push_back(hBox);
						else
							boxes2.push_back(hBox);			
					}//end if
			}//end for

			//the test assignment is finished
			//now compute the areas

			CAABB test1(boxes1);
			CAABB test2(boxes2);
			Real dArea = test1.Area()+test2.Area();
			if(MArea > dArea)
			{
				bestAxis = axis;
				MArea    = dArea;
				if(axis==CAABB::xAxis)
					bestPos  = box->m_Vertices[0].x;
				else
					bestPos  = box->m_Vertices[0].y;
			}//end if

		}//end for
	}//end for

	/*ASSIGNMENT PHASE:
	* assign the parents items to the children 
	*/
	
	vector<COBB*> nBoxes[2];
	
	/* split the items into two buckets relative to the split axis */
	for(bIter = pNode->Boxes().begin(); bIter != pNode->Boxes().end(); bIter++)
	{
		COBB *hBox = *bIter;
	
		/* if longest axis = xAxis */
		if(bestAxis == CAABB::xAxis)
		{
			/* assign split axis */
			pNode->m_iSplitAxis = bestAxis;

			/* value at that the bounding volume is split along the split axis */
			if(hBox->getCenter().x < bestPos)
				Buckets[nBucket].push_back(hBox);
			else
				Buckets[nBucket+1].push_back(hBox);
		}//end if
		/* if longest axis = yAxis */
		else if(bestAxis == CAABB::yAxis)
		{
			/* assign split axis */
			pNode->m_iSplitAxis = bestAxis;

			/* value at that the bounding volume is split along the split axis */
			if(hBox->getCenter().y < bestPos)
				Buckets[nBucket].push_back(hBox);
			else
				Buckets[nBucket+1].push_back(hBox);			
		}//end if
		else
		{

		}//end else


	}//end for

	/* create the new children and assign the items to the children */
	pNode->m_Children[0] = new CAABBNode(Buckets[nBucket],Buckets, nBucket);
	pNode->m_Children[1] = new CAABBNode(Buckets[nBucket+1],Buckets, nBucket+1);

}//end MinArea

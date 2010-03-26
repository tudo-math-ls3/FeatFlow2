////////////////////////////////////////
//
//	CBCNode:
//			 - a node in a hierarchy
//			   of bounding circles
//
////////////////////////////////////////

#include <queue>
#include <bctree.h>
#include <sbbtree.h>



CBCNode::CBCNode(CSBBNode *pNode)
{
	/* set pointer to NULL */
	m_Children[0] = NULL;
	m_Children[1] = NULL;

	m_vCenter = pNode->m_Box->getCenter();

	int i;

	m_dRad = -1.7E+308;

	/* determine bounding circle radius */
	for(i = 0; i < 4; i++)
	{
		VECTOR2 vVertex = pNode->m_Box->m_Vertices[i];

		vVertex  = VECTOR2::createVector(m_vCenter, vVertex);

		Real d = vVertex.mag();

		if(m_dRad < d)
			m_dRad = d;

	}//end for	

	m_Circle.InitCircle(m_vCenter, m_dRad);

}//end constructor


CBCNode::~CBCNode()
{
	if(m_Children[0])
	{
		delete m_Children[0];
		m_Children[0] = NULL;
	}//end if

	if(m_Children[1])
	{
		delete m_Children[1];
		m_Children[1] = NULL;
	}//end if

}//end deconstructor

void CBCNode::GenerateLowerBound(VECTOR2 *vPoints, int nSize)
{

	Real minDist = 1.7E+308;

	for(int i = 0; i < nSize; i++)
	{
		VECTOR2 vV1 = VECTOR2::createVector(m_vCenter, vPoints[i]);
		Real dist   = vV1.mag();
		if(minDist > dist)
		{
			minDist = dist;
			m_vUpper = vPoints[i];
		}//end 

	}//end for
}//end GenerateLowerBound();

void CBCNode::DeleteSubtree(CBCNode *pNode)
{
	if(pNode != NULL)
	{
		DeleteSubtree(pNode->m_Children[0]);
		pNode->m_Children[0] = NULL;
		DeleteSubtree(pNode->m_Children[1]);
		pNode->m_Children[1] = NULL;
		delete pNode;
		pNode = NULL;
	}//end if

}//end deleteSubree

/*!
    \fn CBCNode::NumChildren()
 */
int CBCNode::NumChildren()
{
    return 2;
}


/*!
    \fn CBCNode::GetChildren()
 */
CBCNode** CBCNode::GetChildren()
{
  return m_Children;
}


/*!
    \fn CBCRoot::NumChildren()
 */
int CBCRoot::NumChildren()
{
    return m_nNumChildren;
}


/*!
    \fn CBCRoot::GetChildren()
 */
CBCNode** CBCRoot::GetChildren()
{
	return m_Successors;
}


/*!
    \fn CBCRoot::InitSuccessors(int nSucc)
 */
void CBCRoot::InitSuccessors(int nSucc)
{
    this->m_Successors = new CBCNode*[nSucc];
}


/*!
    \fn CBCRoot::SetChildren(int i, CBCNode* pNode)
 */
void CBCRoot::SetChildren(int i, CBCNode* pNode)
{
    m_Successors[i] = pNode;
}

//Deconstructor
CBCTree::~CBCTree()
{

	DestroyTree();

}//end deconstructor

void CBCTree::DeleteSubTree(CBCNode *pNode)
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

void CBCTree::CopySubTree(CBCNode *pNode)
{

}//end CopySubTree

//Free node memory
void CBCTree::DestroyTree()
{

	for(int i = 0; i < m_iNumChildren; i++)
	{
		DeleteSubTree(m_Children[i]);
	}//end for

}//end DestroyTree

///==============================================
///	      CBCTree initialisation
///==============================================

//===============================================
//
//	Inits a circle tree for a curve with the
//	quadtree method. That means the circle tree
//	is based on the AABBTree built for the curve.
//
//
//	function parameters:
//	-vBucket: a reference to a vector of integers
//	-vSamples: the curve's sampling points
//	-pTrees: the AABBTree data structure
//	-iNumSamples: the number of samples
//	-iNumChildren: the tree's degree
//	-iID: a handle to the curve
//
//===============================================

void CBCTree::InitTree(vector<int> *&vBucket, VECTOR2 vSamples[], CSBBTree **pTrees, int iNumSamples, int iNumChildren, int iID)
{

	//set the number of children
	m_iNumChildren = iNumChildren;

	//allocate memory for the children of the root
	m_Children = new CBCNode*[m_iNumChildren];

	//a counter that holds the storage information
	int iCount = 0;

	int i;

	//queue data structures for the Top-Down tree construction

	//holds the AABBTree's nodes
	queue<CSBBNode*> qBFS;

	//the circle tree's nodes
	queue<CBCNode*>  qNodes;

	/* create the top level nodes in the hierarchy */
	for(i = 0; i < m_iNumChildren; i++)
	{
		//insert the AABBTrees nodes into the queue
		//and construct circle tree nodes from them
		CSBBNode *SBBNode = pTrees[i]->getRoot();
		m_Children[i] = new CBCNode(SBBNode);
		//set the curve ID
		m_Children[i]->SetID(iID);
		qBFS.push(SBBNode);
		qNodes.push(m_Children[i]);
	}//end for
	
	/* Top-Down build of the tree */
	while(!qBFS.empty())
	{
		CSBBNode *SBBNode = qBFS.front();
		CBCNode  *BCNode  = qNodes.front();

		/* create a lower bound for the node by a heuristic */
		BCNode->GenerateLowerBound(vSamples, iNumSamples);

		qBFS.pop();
		qNodes.pop();
		/* an inner node is subdivided */
		if(!SBBNode->IsLeaf())
		{			
			qBFS.push(SBBNode->m_Children[0]);
			qBFS.push(SBBNode->m_Children[1]);
			BCNode->m_Children[0] = new CBCNode(SBBNode->m_Children[0]);
			BCNode->m_Children[0]->SetID(iID);
			BCNode->m_Children[1] = new CBCNode(SBBNode->m_Children[1]);
			BCNode->m_Children[1]->SetID(iID);
			qNodes.push(BCNode->m_Children[0]);
			qNodes.push(BCNode->m_Children[1]);
		}//end if
		/* for a leaf node the data is assigned */
		else
		{
		//	AssignPointsAndParameters(pCurve,BCNode,nIndex);
			AssignPoints(BCNode, vSamples, vBucket[iCount], iNumSamples, iCount);
			iCount++;
		//	m_Leaves.push_back(BCNode);
		}//end else
	}//end while

}//end InitTree



//=====================================================
//
//	The function assigns the data sampling 
//	points to the node.
//
//	function parameters:
//	-pNode: the node that sampling data is assigned
//	-vSamples: the curves sampling points
//  -vBucket: the vector that holds the point indices
//			  of the points that are inside the node's
//			  bounding circle.
//	-iNumSamples: the number of sampling points
//	-iData: a handle to vBucket
//
//=====================================================

void CBCTree::AssignPoints(CBCNode *pNode, VECTOR2 vSamples[] ,vector<int> &vBucket, int iNumSamples, int iData)
{

	//loop over all sampling points
	for(int i = 0; i < iNumSamples; i++)
	{
		//calculate the distance to the bounding circle's center
		VECTOR2 vV1 = VECTOR2::createVector(pNode->GetCenter(), vSamples[i]);
		double dist   = vV1.mag();

		//if the distance is less than the radius the point is inside
		if(dist <= pNode->GetRad())
		{
			//point is inside, insert it into the vector
			vBucket.push_back(i);
		}//end if
	}//end for
	
	//set the handle to the node's data
	pNode->SetData(iData);

}//end AssignPoints


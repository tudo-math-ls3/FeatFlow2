#include <sbbtree.h>
#include <nurbs.h>

//////////////////////////////////////////
//
//		class CSBBNode
//
//	class in a hierarchy of boxes
//  
//
//////////////////////////////////////////

CSBBNode::CSBBNode(CSBB *pBox)
{

	m_Box   = new CSBB(pBox);

	extern int _ID;
	m_Box->m_ID = _ID;
	_ID++;

	m_Children[0] = NULL;

	m_Children[1] = NULL;
	

}//end constructor

CSBBNode::CSBBNode()
{

	m_Box	= NULL;

	m_Children[0] = NULL;

	m_Children[1] = NULL;

}//end constructor


CSBBNode::~CSBBNode()
{

	if(m_Box)
	{
		delete m_Box;	
		m_Box = NULL;
	}

}//end deconstructor

void CSBBNode::subdivideNode(const VECTOR2 &vPoint,const VECTOR2 &vP1, const VECTOR2 &vP2,  Real dValue)
{
	CSBB *box = new CSBB(vP1, vPoint, m_Box->m_fP1, dValue, m_Box->m_pCurve);
	m_Children[0] = new CSBBNode(box);
	delete box;

	box = new CSBB(vPoint, vP2, dValue, m_Box->m_fP2, m_Box->m_pCurve);
	m_Children[1] = new CSBBNode(box);
	delete box;
}//end subdivideNode

Real CSBBNode::GetLowerBoundDebug(const VECTOR2 &vQuery, VECTOR2 &vFar)
{
	cout<<"Debug Mode: "<<endl;
	
	VECTOR2 vCenter = (m_Box->m_Vertices[0] + m_Box->m_Vertices[1] + m_Box->m_Vertices[2] + m_Box->m_Vertices[3]) * 0.25;

	if( (vQuery.x >= vCenter.x) && (vQuery.y >= vCenter.y) )
	{
		VECTOR2 c = VECTOR2::createVector(vQuery, m_Box->m_Vertices[0]);
		vFar = m_Box->m_Vertices[0];
		return c.mag();
	}
	else if( (vQuery.x >= vCenter.x) && (vQuery.y < vCenter.y) )
	{
		VECTOR2 c = VECTOR2::createVector(vQuery, m_Box->m_Vertices[3]);
		vFar = m_Box->m_Vertices[3];
		return c.mag();
	}
	else if( (vQuery.x < vCenter.x) && (vQuery.y >= vCenter.y) )
	{
		VECTOR2 c = VECTOR2::createVector(vQuery, m_Box->m_Vertices[1]);
		vFar = m_Box->m_Vertices[1];
		return c.mag();
	}
	else
	{
		VECTOR2 c = VECTOR2::createVector(vQuery, m_Box->m_Vertices[2]);
		vFar = m_Box->m_Vertices[2];
		return c.mag();		
	}
	
}//end LowerBoundDebug

//Deprecated
Real CSBBNode::GetUpperBound(const VECTOR2 &vQuery)
{
	return 1.0;
}//end GetUpperBound

/////////////////////////////////////////////////
//
//	Class CSBBTree controls a Hierarchy of
//	segment bounding boxes.
//
//
/////////////////////////////////////////////////

/* constructor */
CSBBTree::CSBBTree()
{
	m_Root = NULL;
}//end constructor

CSBBTree::CSBBTree(CSBBNode *tNode) 
{
	m_Root = tNode;
}//end constructor


void CSBBTree::deleteSubTree(CSBBNode *pNode)
{
	if(pNode != NULL)
	{
		deleteSubTree(pNode->m_Children[0]);
		pNode->m_Children[0] = NULL;
		deleteSubTree(pNode->m_Children[1]);
		pNode->m_Children[1] = NULL;
		delete pNode;
		pNode = NULL;
	}//end if

}//end deleteSubree

void CSBBTree::subdivideNode(CSBBNode *pNode, const VECTOR2 &vPoint,const VECTOR2 &vP1, const VECTOR2 &vP2,  Real dValue)
{
	CSBB *box = new CSBB(vP1, vPoint, pNode->m_Box->m_fP1, dValue, pNode->m_Box->m_pCurve);
	pNode->m_Children[0] = new CSBBNode(box);
	delete box;

	box = new CSBB(vPoint, vP2, dValue, pNode->m_Box->m_fP2, pNode->m_Box->m_pCurve);
	pNode->m_Children[1] = new CSBBNode(box);
	delete box;
}//end subdivideNode

void CSBBTree::translateTree(CSBBNode *pNode, const VECTOR2 &vTrans)
{

	if(pNode != NULL)
	{
		translateTree(pNode->m_Children[0], vTrans);
		pNode->m_Box->translateBox(vTrans);
		translateTree(pNode->m_Children[1], vTrans);
	}//end if

}//end inOrder

/*!
    \fn CSBBTree::SubdivideNode(CSBBNode *pNode, const VECTOR2 &vPoint, Real dValue)
 */
void CSBBTree::SubdivideNode(CSBBNode *pNode, const VECTOR2 &vPoint,  Real dValue)
{
	std::vector<VECTOR2> samples1;
	std::vector<VECTOR2> samples2;
	
	pNode->m_Box->m_pCurve->GenSamples(samples1, 150 ,pNode->m_Box->m_fP1, dValue);
	pNode->m_Box->m_pCurve->GenSamples(samples2, 150 ,dValue, pNode->m_Box->m_fP2);
			
	
	CSBB *box = new CSBB(samples1, pNode->m_Box->m_fP1, dValue, pNode->m_Box->m_pCurve);
	pNode->m_Children[0] = new CSBBNode(box);
	delete box;

	box = new CSBB(samples2, dValue, pNode->m_Box->m_fP2, pNode->m_Box->m_pCurve);
	pNode->m_Children[1] = new CSBBNode(box);
	delete box;
}//SubdivideNode


/*!
    \fn CSBBNode::SubdivideNode(CSBBNode *pNode, Real dValue)
 */
void CSBBNode::SubdivideNode(Real dValue)
{
	std::vector<VECTOR2> samples1;
	std::vector<VECTOR2> samples2;
	
	m_Box->m_pCurve->GenSamples(samples1, 150 ,m_Box->m_fP1, dValue);
	m_Box->m_pCurve->GenSamples(samples2, 150 ,dValue, m_Box->m_fP2);
			
	
	CSBB *box = new CSBB(samples1, m_Box->m_fP1, dValue, m_Box->m_pCurve);
	m_Children[0] = new CSBBNode(box);
	delete box;

	box = new CSBB(samples2, dValue, m_Box->m_fP2, m_Box->m_pCurve);
	m_Children[1] = new CSBBNode(box);
	delete box;
}


#ifndef _BCTREE_H_
#define _BCTREE_H_

///==============================================
///	               INCLUDES
///==============================================

#include <vector>
#include <assert.h>
#include "vector2.h"
#include "matrix2x2.h"
#include "circle.h"

///==============================================
///	          Forward Declarations
///==============================================

class CBasicCurve;
class CBSpline;
class CApproxCurve;
class CSBBNode;
class CSBBTree;
class CBCTree;


///==============================================
///	              CLASS CBCNode
///==============================================

class CBCNode
{
public:
	/* standard constructor */
	CBCNode(){ m_Children[0] = NULL; m_Children[1] = NULL;};
	CBCNode(CBSpline *pCurve);

	//copy data from a sbb node
	CBCNode(CSBBNode *pNode);


	/* deconstructor */
	~CBCNode();

	/* inline member functions */

	//is the node a leaf?
	inline bool		 IsLeaf() { return (m_Children[0]== NULL);};

	//returns the nodes center
	inline VECTOR2   GetCenter(){return m_vCenter;};

	//returns the radius of the bounding circle
	inline Real      GetRad(){return m_dRad;};

	//returns the minimum distance from vQuery to the bounding circle
	inline Real    	 GetLowerBound(const VECTOR2 &vQuery){Real dDist = VECTOR2::createVector(m_vCenter,vQuery).mag(); return (dDist - m_dRad);};

	//sets the data index
	inline void		 SetData(short nData){m_Data = nData;};

	//returns the an upper bound for  distance to the
	inline Real	 	 GetUpperBound(const VECTOR2 &vQuery){return VECTOR2::createVector(m_vUpper,vQuery).mag();};
//	inline Real	 	 GetUpperBound(const VECTOR2 &vQuery){return VECTOR2::createVector(m_vUpper,vQuery).norm2();};
	inline int   	 GetData()
	{
#ifdef _DEBUG
		if(!IsLeaf())
		{
			cout<<"warning request data from inner node"<<endl;
			cout<<m_ID<<endl;
			return -1;
		}
#endif
		return m_Data;
	};
	inline void      SetRadius(Real r){m_dRad = r;};
	inline short	 GetID(){return m_ID;};
	inline void 	 SetID(short sID){m_ID = sID;};
	inline void		 Translate(VECTOR2 vTrans){m_vCenter = m_vCenter + vTrans; m_vUpper += vTrans; };
	inline void      Rotate(MATRIX2X2 &rotMat){m_vCenter = rotMat * m_vCenter; m_vUpper = rotMat * m_vUpper;};
	inline void		 SetCenter(VECTOR2 vCenter){m_vCenter = vCenter;};
	inline VECTOR2 	 GetUpperVec(){return m_vUpper;};

	inline bool DetectSimpleColl(CBCNode* otherNode)
	{

		VECTOR2 vCenterJ = otherNode->GetCenter();
		Real    nRadJ    = otherNode->GetRad();
		Real nDist = VECTOR2::createVector(m_vCenter,vCenterJ).mag();

		if((m_dRad+nRadJ) > nDist )
			return true;
		else
			return false;

	}//end DetectSimpleColl

	/* member functions */
	void   GenerateLowerBound(VECTOR2 *vPoints, int nSize);
	void   DeleteSubtree(CBCNode *pNode);
    virtual int NumChildren();
    virtual CBCNode** GetChildren();
	


	
	/* public member variables */
	CBCNode *m_Children[2];

private:

	/* private member variables */
	VECTOR2 m_vCenter;
	VECTOR2 m_vUpper;
	Real    m_dRad;
	short	  m_Data;
	short	  m_ID;
	CCircler m_Circle;
	
    
};

class CBCTree
{
public:

	//constructor
	CBCTree() : m_iNumChildren(0){};

	//deconstructor
	~CBCTree();

	//member functions
	void InitTree(vector<int> *&vBucket, VECTOR2 vSamples[], CSBBTree **pTrees, int iNumSamples, int iNumChildren, int iID);

	//returns the number of children
	inline int GetNumChildren()
	{
		return m_iNumChildren;
	}//end GetNumChildren

	inline CBCNode* GetChild(int i) {return m_Children[i];};
	
	//Free node memory
	void DestroyTree();

	CBCNode **m_Children;

private:

	void DeleteSubTree(CBCNode *pNode);
	void CopySubTree(CBCNode *pNode);
	void AssignPoints(CBCNode *pNode, VECTOR2 vSamples[] ,vector<int> &vBucket, int iNumSamples, int iData);

	int m_iNumChildren;

};

class CBCRoot :public CBCNode
{
	public:
		/* constructor */
	CBCRoot(){ m_Children[0] = NULL; m_Children[1] = NULL;};
    CBCNode** GetChildren();
    int NumChildren();
    void InitSuccessors(int nSucc);
    void SetChildren(int i, CBCNode* pNode);
	
	
	CBCNode **m_Successors;	
		
	int m_nNumChildren;
    
};







typedef struct {
	CBCNode   *node;
	VECTOR2 point;
}NodePoint;

struct CompareNode 
{
	bool operator() (NodePoint n1, NodePoint n2)
	{
		Real l1, l2;
		l1 = n1.node->GetLowerBound(n1.point);
		l2 = n2.node->GetLowerBound(n2.point);
		return (l1 > l2);
	}
};

#endif
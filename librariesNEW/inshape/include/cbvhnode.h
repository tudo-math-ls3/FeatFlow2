/***************************************************************************
 *   Copyright (C) 2007 by Raphael Münster   *
 *   raphael@reynolds   *
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
#ifndef CBVHNODE_H
#define CBVHNODE_H

#include <queue>
#include <list>
#include "vector3.h"
#include "dynamicarray.h"
#include "Triangle3.h"
#include "aabb3.h"
#include "bvhsubdivider.h"
#include "model3d.h"

/**
	@author Raphael Münster <raphael@reynolds>
	
	A class that holds the node information of
	a hierarchy of arbitrary bounding volumes.
	
*/

///==============================================
///	             DEFINITIONS
///==============================================

using namespace std;


template<class BV, class T>
class CBVHTree;


///==============================================
///	             CLASS CBVHNODE
///==============================================

template<class BV, class T>
class CBVHNode
{
public:
	CBVHNode();
	CBVHNode(CBVHNode<BV, T> *pNode);

    ~CBVHNode();
    T GetLowerBound(const CVector3<T> &vQuery);
    T GetUpperBound(const CVector3<T> &vQuery);

    T GetLowerBoundSqr(const CVector3<T> &vQuery);
    T GetUpperBoundSqr(const CVector3<T> &vQuery);

    /*!
        \fn CBVHNode::SetData(short nData)
     */
    void SetData(int nData)
    {
        m_Data = nData;
    }

    /*!
        \fn CBVHNode::GetData()
     */
    int GetData() const
    {
        return m_Data;
    }

    /*!
        \fn CBVHNode::GetRadX()
     */
    T GetRadX()
    {
        /// @todo implement me
    }

    /*!
        \fn CBVHNode::GetRadY()
     */
    T GetRadY()
    {
        /// @todo implement me
    }

    /*!
        \fn CBVHNode::GetCenter()
     */
    CVector3<T> GetCenter() const
    {
       return m_BV.GetCenter();
    }

    /*!
	\fn CBVHNode::GetID()
	 */
	int GetID() const
	{
       return m_ID;
	}

    /*!
	\fn CBVHNode::SetID(int iID)
	 */
	void SetID(int iID)
	{
       m_ID = iID;
	}

	inline bool IsLeaf()
	{
		return (m_Children[0] == NULL && m_Children[1] == NULL);
	}//end IsLeaf

	void SetTriangles(std::vector<CTriangle3f*> *Tris) {Bucket = Tris;};

	BV& GetBV() {return m_BV;};

	vector<CTriangle3f*>& GetTriangles(){return *Bucket;};

	vector<CVector3<T>*>& GetVertices() { return m_vVertices;};

	void InnerVertices(vector<CVector3<T>*>&);
	void InnerVertices(const CDynamicArray< CVector3<T> >&);

	
///===Member variables===	


	BV m_BV;
	CBVHNode<BV,T> *m_Children[2];

	friend class CBVHTree<BV,T>;
	
	private:

		short   m_Data;
		int     m_ID;
		CVector3<T> m_vCenter;
		vector<CTriangle3f*> *Bucket;
		vector<CVector3<T>*> m_vVertices;
	
};
template<class BV, class T>
CBVHNode<BV, T>::CBVHNode(CBVHNode<BV, T> *pNode)
{

	m_vCenter = pNode->GetCenter();
	Bucket = &pNode->GetTriangles();
	m_vVertices = pNode->GetVertices();
	m_Children[0] = 0; m_Children[1] = 0;

}//end constructor

template<class BV, class T>
CBVHNode<BV, T>::CBVHNode()
{

	m_Children[0] = 0; m_Children[1] = 0;

}//end constructor

template<class BV, class T>
CBVHNode<BV, T>::~CBVHNode()
{

	if(m_Children[0])
	{
		delete m_Children[0];
		m_Children[0] = 0;

		delete m_Children[1];
		m_Children[1] = 0;
	}//end if

}//end deconstructor

/*!
    \fn CBVHNode::GetLowerBound()
 */
template<class BV, class T>
T CBVHNode<BV, T>::GetLowerBound(const CVector3<T> &vQuery) 
{
    return m_BV.MinDistance(vQuery);
}


/*!
    \fn CBVHNode::GetUpperBound()
 */
template<class BV, class T>
T CBVHNode<BV, T>::GetUpperBound(const CVector3<T> &vQuery) 
{
	return m_BV.MaxDistance(vQuery);
}

/*!
    \fn CBVHNode::GetLowerBoundSqr()
 */
template<class BV, class T>
T CBVHNode<BV, T>::GetLowerBoundSqr(const CVector3<T> &vQuery) 
{
    return m_BV.MinDistanceSqr(vQuery);
}


/*!
    \fn CBVHNode::GetUpperBoundSqr()
 */
template<class BV, class T>
T CBVHNode<BV, T>::GetUpperBoundSqr(const CVector3<T> &vQuery) 
{
	return m_BV.MaxDistanceSqr(vQuery);
}

template<class BV, class T>
void CBVHNode<BV, T>::InnerVertices(vector<CVector3<T>*>& vVertices)
{
    
	for(int i = 0; i < vVertices.size(); i++)
	{
		if(m_BV.Inside(*vVertices[i]))
			m_vVertices.push_back(vVertices[i]);
	}//end for

}

template<class BV, class T>
void CBVHNode<BV, T>::InnerVertices(const CDynamicArray< CVector3<T> > &vVertices)
{

	for(int i = 0; i < vVertices.Size(); i++)
	{
		if(m_BV.Inside(vVertices[i]))
			m_vVertices.push_back(&vVertices[i]);
	}//end for

    
}

///==============================================
///	             CLASS CBVHTREE
///==============================================

template<class BV, class T>
class CBVHTree
{
public:
	//constructors
	CBVHTree(const vector<BV*> &vPrimitives, int iMethod);

	CBVHTree(const CBVHTree &pTree);

	CBVHTree(){m_Root = 0;};
	
	//deconstructor
	~CBVHTree();

	//member functions
	void BuildTree(vector< CTriangle3<T>* > &vPrimitives, const CBVHSubdivider &Subdivider);
	void BuildTree(CModel3D &model, const CBVHSubdivider &Subdivider);
	void UpdateTree();
	void DestroyTree();
	void DeleteSubTree(CBVHNode<BV,T> *pNode);
	void CopySubTree(CBVHNode<BV,T> *pNode);
	void SubdivideNode(CBVHNode<BV,T> *pNode);
	bool EvaluateTermiationCrit(CBVHNode<BV,T> *pNode, int iDepth);
	int  DepthFirstSearch(CBVHNode<BV,T> *pNode, int iDepth);

	inline CBVHNode<BV,T>* GetRoot() const {return m_Root;};


	enum
	{
		SPATIAL_MEDIAN,
		MEDIAN_CUT,
		MINAREA
	};

private:
	//pointer to the root node
	CBVHNode<BV,T> *m_Root;

};

template<class BV, class T>
CBVHTree<BV, T>::CBVHTree(const CBVHTree &pTree)
{

	m_Root = new CBVHNode<BV,T>(pTree.GetRoot());
	
	/* a queue stores the nodes to be subdivided */
	queue<CBVHNode<BV,T>*> qNodesOrig;
																			
	/* initially only the root is in the queue */
	qNodesOrig.push(pTree.GetRoot());
	
	/* a queue stores the nodes to be subdivided */
	queue<CBVHNode<BV,T>*> qNodes;
																			
	/* initially only the root is in the queue */
	qNodes.push(m_Root);
	/* build the bounding box hierarchy */
	while(!qNodesOrig.empty())
	{

		int num = (int)qNodes.size();
		/* get front element */
		CBVHNode<BV,T> *origNode = qNodesOrig.front();

		/* remove front element */
		qNodesOrig.pop();
		
		CBVHNode<BV,T> *pNode = qNodes.front();
		
		qNodes.pop();
		
		if(origNode->IsLeaf())
		{
			/* NODE IS A LEAF DO NOTHING */
		
		}//end if
		else
		{
			/* copy nodes */
			pNode->m_Children[0] = new CBVHNode<BV,T>(origNode->m_Children[0]);
			pNode->m_Children[1] = new CBVHNode<BV,T>(origNode->m_Children[1]);
		
			/* insert original nodes */
			qNodesOrig.push(origNode->m_Children[0]);
			qNodesOrig.push(origNode->m_Children[1]);
			
			/* insert nodes of the copy */
			qNodes.push(pNode->m_Children[0]);
			qNodes.push(pNode->m_Children[1]);
		

		}//end else
			
	}//end while

}//end constructor

template<class BV, class T>
CBVHTree<BV, T>::CBVHTree(const vector<BV*> &vPrimitives, int iMethod)
{

	m_Root = new CBVHNode<BV,T>();

	//Top Down built: store the nodes that are subdivided
	//use a queue
	queue< CBVHNode<BV,T>* > qNodes;

	//insert the root into the queue
	qNodes.push(m_Root);

	while(!qNodes.empty())
	{
		//get the first element
		CBVHNode<BV,T> *pNode = qNodes.front();

		//remove it from the queue
		qNodes.pop();

		//should the node be further subdivided?
		if(!EvaluateTermiationCrit(pNode,7))
		{
			//subdivide
			qNodes.push(pNode->m_Children[0]);
			qNodes.push(pNode->m_Children[1]);
		}//end if
		else
		{
			//stop subdivision
		}//end else

	}//end while

}//end constructor

template<class BV, class T>
void CBVHTree<BV, T>::SubdivideNode(CBVHNode<BV,T> *pNode)
{

}//end SubdivideNode

template<class BV, class T>
bool CBVHTree<BV, T>::EvaluateTermiationCrit(CBVHNode<BV,T> *pNode, int iDepth)
{
	if(DepthFirstSearch(m_Root,0) >= 6)
	{
		return true;
	}
	else
		return false;
}//end SubdivideNode

template<class BV, class T>
void CBVHTree<BV, T>::DeleteSubTree(CBVHNode<BV,T> *pNode)
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

template<class BV, class T>
void CBVHTree<BV, T>::DestroyTree()
{

	if(m_Root != NULL)
	{
		DeleteSubTree(m_Root);
	}//end if

}//end DeleteSubTree

template<class BV, class T>
CBVHTree<BV, T>::~CBVHTree()
{

	if(m_Root != NULL)
	{
		DeleteSubTree(m_Root);
	}//end if

}//end deconstructor

template<class BV, class T>
void CBVHTree<BV, T>::BuildTree(vector< CTriangle3<T>* > &vPrimitives, const CBVHSubdivider &Subdivider)
{

	m_Root = new CBVHNode<BV, T>();

	//Top Down built: store the nodes that are subdivided
	//use a queue
	queue< CBVHNode<BV, T>* > qNodes;

	//insert the root into the queue
	qNodes.push(m_Root);
	m_Root->Bucket = &vPrimitives;
	m_Root->GetBV().Init(vPrimitives);

	int iNodes = 0;

	while(!qNodes.empty())
	{
		//get the first element
		CBVHNode<BV, T> *pNode = qNodes.front();

		//remove it from the queue
		qNodes.pop();

		//should the node be further subdivided?
		if(!EvaluateTermiationCrit(pNode, 7))
		{
			//subdivide node(pNode)
			cout<<"Number of triangles: "<<pNode->GetTriangles().size()<<endl;
			Subdivider.SubdivideNode(pNode, pNode->Bucket);
			qNodes.push(pNode->m_Children[0]);
			qNodes.push(pNode->m_Children[1]);
		}//end if
		else
		{
			//stop subdivision
			break;

		}//end else

		iNodes++;

	}//end while

}//end BuildTree

template<class BV, class T>
int CBVHTree<BV, T>::DepthFirstSearch(CBVHNode<BV,T> *pNode, int iDepth)
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

	return iReturn;

}//end DepthFirstSearch

template<class BV, class T>
void CBVHTree<BV, T>::BuildTree(CModel3D &model, const CBVHSubdivider &Subdivider)
{

	m_Root = new CBVHNode<BV, T>();

	//Top Down built: store the nodes that are subdivided
	//use a queue
	queue< CBVHNode<BV, T>* > qNodes;

	//insert the root into the queue
	qNodes.push(m_Root);
	
	m_Root->Bucket = &model.GetTriangles();;
	m_Root->GetBV().Init(model.GetTriangles());
	m_Root->InnerVertices(model.GetVertices());

	int iNodes = 0;

	while(!qNodes.empty())
	{
		//get the first element
		CBVHNode<BV, T> *pNode = qNodes.front();

		//remove it from the queue
		qNodes.pop();

		//should the node be further subdivided?
		if(!EvaluateTermiationCrit(pNode, 7))
		{
			//subdivide node(pNode)
			Subdivider.SubdivideNode(pNode, pNode->Bucket, pNode->GetVertices() );
			qNodes.push(pNode->m_Children[0]);
			qNodes.push(pNode->m_Children[1]);
		}//end if
		else
		{
			//stop subdivision
			break;

		}//end else

		iNodes++;

	}//end while

}//end BuildTree

typedef CBVHTree<CAABB3f,float> AABBTree3f;
typedef CBVHNode<CAABB3f,float> AABBNode3f;

#endif

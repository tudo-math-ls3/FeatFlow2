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

/*
 * This class is used to store the information of a bounding box for a
 * monotonic curve segment of a NURBS or BSpline curve. The vector
 * member variables store the corner vertices of the box. The float member
 * variables are the parameter values of the two box vertices that lie
 * on the curve.
 *
 *
 *
 */
#if !defined(_CSBB_H)
#define _CSBB_H

#include <vector>
#include <list>
#include <iostream>


#include "vector2.h"
#include "matrix2x2.h"
#include "dynamicarray.h"
#include "mathglobals.h"


using namespace std;

/* forward declaration */
class CBasicCurve;
class CBSpline;
class CApproxCurve;
class COBB;

class CAABB
{
public:
	
	/* constructors */
	CAABB(const VECTOR2 &vA, const VECTOR2 &vB);
	CAABB(const vector<COBB*> &vBoxes);
	CAABB(){};
	
	/* member functions */
	bool getLongestAxis();
	void calculateBounds(const vector<COBB*> &vBoxes);
	void calculateBounds(CAABB *b1, CAABB *b2);
	
	virtual bool inside(const VECTOR2& p);
	
	/* inline member functions */	
	virtual inline VECTOR2 getCenter()
	{
		VECTOR2 vCenter;
		vCenter = (m_Vertices[0] + m_Vertices[1] +
				   VECTOR2(m_Vertices[0].x, m_Vertices[1].y) + VECTOR2(m_Vertices[1].x, m_Vertices[0].y)) * 0.25;
		return vCenter;
	};

	inline Real Area(){return (m_Vertices[1].x - m_Vertices[0].x) * (m_Vertices[1].y - m_Vertices[0].y);};

	virtual inline void translateBox(VECTOR2 vTrans)
	{
		m_Vertices[0] = m_Vertices[0] + vTrans;
		m_Vertices[1] = m_Vertices[1] + vTrans;
	};
	
	/* member variables */
	VECTOR2 m_Vertices[4];

	enum
	{
		xAxis,
		yAxis
	};
		
};


class COBB : public CAABB
{
public:
	
	/* constructors */
	COBB(VECTOR2 vA, VECTOR2 vB, CBSpline *pCurve = NULL);
	COBB(CDynamicArray<VECTOR2> &cP, vector<int>& hP, CBSpline *pCurve);
	COBB(COBB* pBox, CApproxCurve *pCurve);
	COBB(const vector<VECTOR2> &cP, const vector<int> &hP, CBSpline *pCurve);
	COBB();

	/* member functions */
	bool    inside(const VECTOR2& p);
	bool    intersection(Real r, Real l, Real t, Real b);
	void    getDimension();
	bool    Overlap(COBB *pBox);

	/* inline member functions */
	inline VECTOR2 getCenter()
	{
		VECTOR2 vCenter;
		vCenter = (m_Vertices[0] + m_Vertices[1] + m_Vertices[2] + m_Vertices[3]) * 0.25;
		return vCenter;
	};
	
	inline void translateBox(VECTOR2 vTrans)
	{
		for(int i = 0; i < 4; i++)
		{
			m_Vertices[i] = m_Vertices[i] + vTrans;
		}//end for
	};

	inline Real  getArea(){return ((m_Vertices[1].x - m_Vertices[0].x)*(m_Vertices[2].y - m_Vertices[1].y));};
	

	/* member variables */
	//make this one an integer ID

	union
	{
		CBSpline *m_pCurve;
		CApproxCurve *m_paCurve;
	};

	friend ostream& operator<<(ostream& out, COBB *b1); 

};


///////////////////////////////
//
//	class segment bounding box
//
//
///////////////////////////////

class CSBB : public COBB
{
public:
	
	/* constructors */
	CSBB(VECTOR2 vA, VECTOR2 vB, Real fA, Real fB, CBSpline *pCurve);
	CSBB(CSBB *pBox);
	CSBB();
	
	/* member functions */
	void    getLUB(const VECTOR2 &vQuery);
	friend ostream& operator<<(ostream& out, CSBB *b1); 

	/* member variables */
	int m_ID;
	/* parameter values */
	Real m_fP1;
	Real m_fP2;
	
	/* 0 = bottom left, 1 = bottom right, 2 = top right, 3 = top left */
	/* stores the two points that are on the curve */
	int m_CurvePoints[2];

};

class CSBBTree;

class CSBBNode
{
public:

	/* constructor */
	CSBBNode(CSBB *pBox);

	CSBBNode();

	/* deconstructor */
	~CSBBNode();

	/* member functions */
	void subdivideNode(const VECTOR2 &vPoint,const VECTOR2 &vP1, const VECTOR2 &vP2,  Real dValue);
	Real getUpperBound(const VECTOR2 &vQuery);
	Real getLowerBoundDebug(const VECTOR2 &vQuery, VECTOR2 &vFar);
	
	
	/* inline member functions */
	inline Real getLowerBound(const VECTOR2 &vQuery)
	{
		Real d = 0.0;
		Real minDist = 1.7E+308;
		for(int j = 0; j < 2; j++)
		{
			VECTOR2 b = m_Box->m_Vertices[m_Box->m_CurvePoints[j]];
			b = VECTOR2::createVector(vQuery, b);
			d = b.mag();
			if(d < minDist)
			{
				minDist   = d;
			}//end if
		}//end for
		return minDist;	
	}//end getLowerBound

	inline bool isLeaf() {return (m_Children[0]==NULL); };


	/* public member variables */
	CSBBNode *m_Children[2];
	CSBB *m_Box;

	/* the tree class is a friend */
	friend class CSBBTree;

};

class CSBBTree
{

public:

	/* constructor */
	CSBBTree();
	CSBBTree(CSBBNode *tNode);
	
	/* deconstructor */
	~CSBBTree();

	/* member functions */
	void deleteSubTree(CSBBNode *pNode);
	void subdivideNode(CSBBNode *pNode, const VECTOR2 &vPoint,const VECTOR2 &vP1, const VECTOR2 &vP2,  Real dValue);
	void translateTree(CSBBNode *pNode, const VECTOR2 &vTrans);
	void destroyTree();

	/* inline member functions */
	inline CSBBNode* getRoot(){return m_Root;};

private:

	CSBBNode *m_Root;

};

typedef struct
{
	VECTOR2 m_vCenter;
	VECTOR2 m_vUpper;
	Real    m_dRad;
	short	  m_Data;
	short	  m_ID;
}tCNode;

class CBCNode
{
public:
	/* constructor */
	CBCNode(){ m_Children[0] = NULL; m_Children[1] = NULL;};
	CBCNode(CBSpline *pCurve);
	CBCNode(CSBBNode *pNode);
	CBCNode(tCNode &tNode);

	/* deconstructor */
	~CBCNode();

	/* inline member functions */
	inline bool		 IsLeaf() { return (m_Children[0]== NULL);};
	inline VECTOR2 GetCenter(){return m_vCenter;};
	inline Real    GetRad(){return m_dRad;};
	inline Real    GetLowerBound(const VECTOR2 &vQuery){Real dDist = VECTOR2::createVector(m_vCenter,vQuery).mag(); return (dDist - m_dRad);};
	inline void		 SetData(short nData){m_Data = nData;};
	inline Real	 GetUpperBound(const VECTOR2 &vQuery){return VECTOR2::createVector(m_vUpper,vQuery).mag();};
	inline int   	 GetData()
	{
#ifdef _DEBUG
		if(!IsLeaf())
		{
			cout<<"warning request data from inner node"<<endl;
			return -1;
		}
#endif
		return m_Data;
	};
	inline short	 GetID(){return m_ID;};
	inline void 	 SetID(short sID){m_ID = sID;};
	inline void		 Translate(VECTOR2 vTrans){m_vCenter = m_vCenter + vTrans; m_vUpper += vTrans; };
	inline void      Rotate(MATRIX2X2 &rotMat){m_vCenter = rotMat * m_vCenter; m_vUpper = rotMat * m_vUpper;};
	inline void		 SetCenter(VECTOR2 vCenter){m_vCenter = vCenter;};
	inline VECTOR2 GetUpperVec(){return m_vUpper;};
	

	/* member functions */
	void   GenerateLowerBound(VECTOR2 *vPoints, int nSize);
	void   DeleteSubtree(CBCNode *pNode);
	
	/* public member variables */
	CBCNode *m_Children[2];

private:

	/* private member variables */
	VECTOR2 m_vCenter;
	VECTOR2 m_vUpper;
	Real    m_dRad;
	short	  m_Data;
	short	  m_ID;
	
    
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

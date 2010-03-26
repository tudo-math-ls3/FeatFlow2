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


#if !defined _CAPPROXCURVE_H_
#define _CAPPROXCURVE_H_

//===================================================
//				     Includes
//===================================================

#include "nurbs.h"
#include "matrix2x2.h"
#include "primops.h"
#include "bctree.h"
#include <cstdlib>
//#include <boost/config.hpp>
#include <iostream>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/property_map.hpp>

using namespace std;
//using namespace boost;

class CMedAxisVertex;

//===================================================
//				Typedefs and structs
//===================================================

typedef struct
{
	double dP;
	double dPStart;
	double dPEnd;
}tParameters;

typedef struct
{
  VECTOR2 vertex;
  int     index;
}VertexIndex;



//enum vertex_MedAxis_t { vertex_MedAxis };
//
//namespace boost {BOOST_INSTALL_PROPERTY(vertex, MedAxis); }
//
////typedef for the custom vertex property
//typedef property<vertex_MedAxis_t, CMedAxisVertex*> MedAxisVP;
//
////typedef for the edge property weight = length with data type double
//typedef property<edge_weight_t, double> EdgeWeight;
//
////An adjacency list for an undirected graph with lists 
//typedef adjacency_list<listS, listS, undirectedS, MedAxisVP, EdgeWeight> UndirectedGraph;
//		
////typedef for vertex property map
//typedef property_map<UndirectedGraph, vertex_MedAxis_t>::type MAProperty;		
//		
////typedef for edge property map		
//typedef property_map<UndirectedGraph, edge_weight_t>::type EWProperty;		
//		
////typedef to refer to the edge_property		
//		
////typename for the three vertex objects
//typedef graph_traits<UndirectedGraph>::vertex_descriptor myVertex;





/*!

///======================CLASS CAPPROXCURVE============================
//
//
//	This class contains sampled curve data to realize efficient
//	minimum distance computation. The curve is approximated by
//	a number of samples. A circle tree data structure is built
//	based on the samples to speed up the computation.
//
//	members: m_ID:  identifies the curve in the global curve array
//			 m_Hierarchy: pointers to the top level hierarchy nodes
//			 m_bBox: the curve's bounding box
//			 m_nResolution: the number of samples
//			 m_Buckets: a vector with length = number of leaves,
//				it contains the items of leaf i at m_Buckets[i]
//
///====================================================================

*/

///==============================================
///	      CLASS CApproxCurve
///==============================================

class CApproxCurve : public CBasicCurve
{

public:

	/* constructor */
	CApproxCurve(CBSpline* pCurve, int Resolution, short ID = 0);
	CApproxCurve();
	CApproxCurve(const char* sFileName);

	/* deconstructor */
	~CApproxCurve();

	/* member functions */
	void ComputeCog();
	void RotateCurve(MATRIX2X2 &rotMat);
	void TranslateCurveByVector(const VECTOR2 &rhs);
	void TranslateCurve(const VECTOR2 &rhs);
	void OutputToFile(char *fileName);
	void CreateHierarchy(CBSpline *pCurve);
	void AssignPoints(CBCNode *&pNode, int n);
	void AssignPointsAndParameters(CBSpline *pCurve, CBCNode *&pNode, int n);
	void TranslateTree(CBCNode *pNode, const VECTOR2 &vTrans);
	void TransformCurve(MATRIX2X2 &rotMat, VECTOR2 &vT);
	void RotateTree(CBCNode *pNode, MATRIX2X2 &rotMat);
	void TransformTree(CBCNode *pNode, MATRIX2X2 &rotMat, VECTOR2 &vTrans1, VECTOR2 &vTrans2);
	bool IsInElement(const VECTOR2 &vQuery);
	bool Inside(const VECTOR2 &vQuery);
	void Dump(char *sFileName);
	void RotateCurve(Real nAngle);

	/// functions for medial axis approximation
	void ApproximateMedialAxis(vector<VECTOR2> &vPoints);
	void ReadEdges(const char* sFileName);
	void FindFormingPoints();
	
	/* inline member function */
	inline double   Radius(){return m_nRadius;};
	inline VECTOR2  GetPoint(int iBucket, int iPos){return m_vSamples[m_Buckets[iBucket][iPos]];};
	inline double   GetParameter(int iBucket, int iPos){return m_dParamBuckets[iBucket][iPos].dP;};
	inline int      GetResolution(){return m_nResolution;};
	inline VECTOR2* GetSamples(){return m_vSamples;};
	inline VECTOR2  GetCenter(){return m_vCenter;};
	inline int      GetNumCircles(){return m_nCircles;};
	inline COBB*    BoundingBox(){return m_bBox;};
	inline short    ID(){return m_ID;};
	inline int      GetBucketSize(int Data){return m_Buckets[Data].size();};

//	void MergeVertices(myVertex &u, myVertex &v, DebugSEC &info);
//	template<typename UndirectedGraph> void BuildBoostGraph();
	//void BuildBoostGraph();
	void MakeInnerEdges(vector<tEdge>	&InnerEdges);
	void DetermineEdge(DebugSEC &info);
//	CCircle GetMergedCircle(myVertex &u, myVertex &v, DebugSEC &info);
	void DetermineOptimalEdge(DebugSEC &info);

    /*!
        \fn CApproxCurve::SetCenter(const VECTOR2 &vCenter)
     */
    void SetCenter(const VECTOR2 &vCenter)
    {
       m_vCenter = vCenter;
    }
 

//	vector<myVertex> GradOne;

	/* public member variables */
	CBCRoot*            m_Root;
	COBB* 				m_bBox;
	vector<int>			*m_Buckets;
	vector<CBCNode*>	m_Leaves;
	vector<tParameters> *m_dParamBuckets;
	vector<tVertex>     m_vVoronoi;
	vector<CMedAxisVertex*>     m_MedAxis;
	vector<tEdge>		m_vEdges;
	double				*m_Radii;
	CMedAxisVertex*     m_TestVertex;
//	UndirectedGraph undirGraph;
	CBCTree             m_BCTree;

private:
	/* member variables */
	double		m_nRadius;
	Real        *m_dParameters;
	VECTOR2	    *m_vSamples;
	VECTOR2     m_vCenter;
   	int			m_nResolution;
	int	        m_nCircles;
	int		    m_nBuckets;
	short	    m_ID;

	
};

#endif

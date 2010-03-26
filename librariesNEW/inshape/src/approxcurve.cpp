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

#include <approxcurve.h>
#include <fileParser.h>
#include <obb2.h>
#include <sbbtree.h>

#include <iostream>
#include <fstream>
#include <queue>
#include <set>
#include <math.h>
#include <limits>

//#ifdef _DEBUG
// #define new DEBUG_NEW
// #define malloc DEBUG_MALLOC
// static char THIS_FILE[] = __FILE__;
//#endif


#define SINGLE -3

int	nIndex = 0;

///==============================================
///	      CApproxCurve constructor
///==============================================

/////////////////////////////////////////////////////////
/*
*
*   CApproxCurve class constructor
*
*	Parameters: - pCurve the corresponding BSpline curve
*			    - Resolution the number of sample points
*				- ID the global curve ID
*/
/////////////////////////////////////////////////////////

CApproxCurve::CApproxCurve(CBSpline *pCurve, int Resolution, short ID) : m_ID(ID)
{


	/* number of samples */
	m_nResolution = Resolution;

	m_Radii = NULL;

	/* allocate memory for samples */
	m_vSamples = new VECTOR2[m_nResolution];
	m_dParameters = new Real[m_nResolution];


	int i;
	
	int nDeg = pCurve->getDeg();

	int nNumP = pCurve->getNumPoints();

	Real end = pCurve->KnotVector()[nNumP-1];

	Real start = pCurve->KnotVector()[nDeg];

	double dT = pCurve->KnotVector()[nNumP] - pCurve->KnotVector()[nDeg];

	if(m_nResolution > 1)
		dT = dT / static_cast<double>(m_nResolution - 1);

	double T = pCurve->KnotVector()[nDeg];

	/* generate sample points */
	for(i = 0; i < m_nResolution; i++)
	{
		
		if(T > pCurve->KnotVector()[nNumP])
			T = pCurve->KnotVector()[nNumP];

		
		m_vSamples[i]    = pCurve->CoxDeBoor(T);
		m_dParameters[i] = T;
		m_vCenter 	    += m_vSamples[i];
		T+=dT;
		
	}//end for

	/* calculate the center of gravity of the curve */
	m_vCenter = m_vCenter * (1.0/double(m_nResolution));
	
	m_nRadius = -1.7E+308;
	
	/* generate a bounding circle */
	for(i = 0; i < m_nResolution; i++)
	{
		VECTOR2 vV1 = VECTOR2::createVector( m_vSamples[i], m_vCenter );
		double d = vV1.mag();
		if(m_nRadius < d)
			m_nRadius = d;
	}//end for
	
	/* get number of root circles */
	m_nCircles = pCurve->getBoxCount();

	/* allocate memory for root circles */
//	m_Hierarchy = new CBCNode*[m_nCircles];

	int branchFactor = 2;

	int nLeaves = static_cast<int>(m_nCircles * pow((double)branchFactor ,int(TARGETDEPTH)));
	m_Buckets = new vector<int>[nLeaves];

	m_dParamBuckets = new vector<tParameters>[nLeaves];
	cout<<"Number of leaf nodes: "<<nLeaves<<endl;

	/* build bounding circle hierarchy */
	m_BCTree.InitTree(m_Buckets, m_vSamples, pCurve->m_Hierarchy, m_nResolution, m_nCircles, m_ID); 
	
	/* get an OBB for the curve */
	m_bBox = new COBB(pCurve->BoundingBox(), this);

	nIndex = 0;
	
	m_dAngle = 0.0;
	
	cPrimOps op;
	
	set<VECTOR2*,comparePoints> vSet;
	
	for(int i = 0; i < m_nResolution;i++)
	{
		vSet.insert(&m_vSamples[i]);
	}//End for
	
	DebugSEC inf;
	
	CCircler circle = op.SEC(vSet,inf);
	
	m_vCenter = circle.GetCenter();
	this->m_nRadius = circle.GetRadius();
	
	
	m_Root = new CBCRoot();
	m_Root->InitSuccessors(m_nCircles);
	for(int k = 0; k < this->m_nCircles; k++)
	{
		m_Root->SetChildren(k, m_BCTree.GetChild(k));
	}//end for
	
	m_Root->SetRadius(circle.GetRadius());
	m_Root->SetCenter(circle.GetCenter());
	
	m_Root->m_nNumChildren = m_nCircles;
	
	
// 	std::vector<VECTOR2> vAll;
// 	readVoronoiVertices("Stern.1.v.node", vAll);
// 	ApproximateMedialAxis(vAll);
// 	ReadEdges("Stern.1.v.edge");
// 	FindFormingPoints();

	

// 	CMedAxisVertex *pVertex; 
// 	CMedAxisVertex *pVertex2;
// 		
// 	
// 	cout<<"# of inner edges: "<<m_vEdges.size()<<endl;
// 	
// 	m_TestVertex = NULL;
// 	
// 	BuildBoostGraph<UndirectedGraph>();
// 	
// 	graph_traits<UndirectedGraph>::vertex_descriptor u, v;
// 	graph_traits<UndirectedGraph>::edge_descriptor e;
// 	
// 	graph_traits<UndirectedGraph>::adjacency_iterator outEdge, out_endEdge;
// // 	graph_traits<UndirectedGraph>::out_edge_iterator outEdge, out_endEdge;
// // 	graph_traits<UndirectedGraph>::in_edge_iterator inEdge, in_endEdge;
// 	graph_traits<UndirectedGraph>::vertex_iterator vIter, vIter_end;
// 
// 	cout << "vertices(g) = "<<num_vertices(undirGraph)<<endl;
// 	cout << "edges(g) = "<<num_edges(undirGraph)<<endl;
// 
// 	
// 	u  = vertex(0, undirGraph);
// 	v  = vertex(1, undirGraph);
// 	bool found;
// 	
// 	tie(e, found) = edge(u,v, undirGraph);
// 	if(found)
// 	{
// 		cout<<"Edge found"<<endl;
// 	}//end if
	
	//typedef for propertymap
/*	typedef property_map<UndirectedGraph, vertex_MedAxis_t>::type MAProperty;
		
	//get the property accessor
	MAProperty prMAV = get(vertex_MedAxis, undirGraph);
	
	typedef property_traits<MAProperty>::value_type MAType;
	
	pVertex   = boost::get(prMAV, u);
	pVertex2  = boost::get(prMAV, v);


 	MergeVertices(pVertex, pVertex2);
	
	int nV = num_vertices(undirGraph);
	int nE = num_edges(undirGraph);*/
	
// 	cout<<"Number of vertices "<<nV<<endl;
// 	cout<<"Number of edges "<<nE<<endl;
	
// 	inEdge = inv_adjacent_vertices(u, undirGraph).first;
// 	in_endEdge = inv_adjacent_vertices(u, undirGraph).second;
	
//	outEdge = adjacent_vertices(u, undirGraph).first;
//	out_endEdge = adjacent_vertices(u, undirGraph).second;
	
// 	outEdge = out_edges(u, undirGraph).first;
// 	out_endEdge = out_edges(u, undirGraph).second;
// 	
// 	inEdge = in_edges(u, undirGraph).first;
// 	in_endEdge = in_edges(u, undirGraph).second;	

/*	clear_vertex(u, undirGraph);
	clear_vertex(v, undirGraph);*/
	

// 	while(out_degree(u, undirGraph) > 0)
// 	{
// 		outEdge = adjacent_vertices(u, undirGraph).first;
// 		remove_edge(u, (*outEdge), undirGraph);
// 	}//end while
// 	
// 	while(out_degree(v, undirGraph) > 0)
// 	{
// 		outEdge = adjacent_vertices(v, undirGraph).first;
// 		remove_edge(v, (*outEdge), undirGraph);
// 	}//end while	
	
  	
//  	remove_vertex(v, undirGraph);		
// 	remove_vertex(u, undirGraph);
	
/*	nV = num_vertices(undirGraph);
	nE = num_edges(undirGraph);
	cout<<"Number of vertices "<<nV<<endl;
	cout<<"Number of edges "<<nE<<endl;*/
	
// 	cout<<"write to file"<<endl;
// 	this->OutputToFile("Stern.node");
	
}//end constructor


///==============================================
///	      CApproxCurve standard constructor
///==============================================

CApproxCurve::CApproxCurve()
{

}//end constructor


//old deprecated function
void CApproxCurve::CreateHierarchy(CBSpline *pCurve)
{

	///* initialize number of buckets */
	//m_nBuckets = 0;

	//int i;
	//queue<CSBBNode*> qBFS;
	//queue<CBCNode*>  qNodes;

	///* create the top level nodes in the hierarchy */
	//for(i = 0; i < pCurve->getBoxCount(); i++)
	//{
	//	CSBBNode *SBBNode = pCurve->m_Hierarchy[i]->getRoot();
	//	m_Hierarchy[i] = new CBCNode(SBBNode);
	//	m_Hierarchy[i]->SetID(m_ID);
	//	qBFS.push(SBBNode);
	//	qNodes.push(m_Hierarchy[i]);
	//}//end for
	//
	///* top down build of the tree */
	//while(!qBFS.empty())
	//{
	//	CSBBNode *SBBNode = qBFS.front();
	//	CBCNode  *BCNode  = qNodes.front();
	//	/* create a lower bound for the node by a heuristic */
	//	BCNode->GenerateLowerBound(m_vSamples, m_nResolution);
	//	qBFS.pop();
	//	qNodes.pop();
	//	/* an inner node is subdivided */
	//	if(!SBBNode->IsLeaf())
	//	{			
	//		qBFS.push(SBBNode->m_Children[0]);
	//		qBFS.push(SBBNode->m_Children[1]);
	//		BCNode->m_Children[0] = new CBCNode(SBBNode->m_Children[0]);
	//		BCNode->m_Children[0]->SetID(m_ID);
	//		BCNode->m_Children[1] = new CBCNode(SBBNode->m_Children[1]);
	//		BCNode->m_Children[1]->SetID(m_ID);
	//		qNodes.push(BCNode->m_Children[0]);
	//		qNodes.push(BCNode->m_Children[1]);
	//	}//end if
	//	/* for a leaf node the data is assigned */
	//	else
	//	{
	//		AssignPointsAndParameters(pCurve,BCNode,nIndex);
	//		m_Leaves.push_back(BCNode);
	//	}//end else
	//}//end while
}//end CreateHierarchy

CApproxCurve::~CApproxCurve()
{

	if(m_vSamples)
	{
		delete[] m_vSamples;
		m_vSamples = NULL;
	}//end if

	if(m_Buckets)
		delete[] m_Buckets;

	if(m_dParamBuckets)
		delete[] m_dParamBuckets;
	
	if(m_dParamBuckets)
		delete[] m_dParameters;
	
	if(m_bBox)
	{
		delete m_bBox;
		m_bBox = NULL;
	}//end if
	
	if(m_Radii)
		delete[] m_Radii;

}//end deconstructor


///////////////////////////////////////////////////////////////
/*
*
*	This function assigns the curve points to the circle tree
*	nodes. Additionally it stores the parameter intervall the
*	current point is in.
*	pNode: this is current circle tree node
*	pCurve: this is the curve for that a circle tree is constructed
*
*/
///////////////////////////////////////////////////////////////

void CApproxCurve::AssignPointsAndParameters(CBSpline *pCurve, CBCNode *&pNode, int n)
{
	
	for(int i = 0; i < m_nResolution; i++)
	{
		VECTOR2 vV1 = VECTOR2::createVector(pNode->GetCenter(), m_vSamples[i]);
		double dist   = vV1.mag();
		if(dist <= pNode->GetRad())
		{
			m_Buckets[n].push_back(i);
			tParameters tParam;
			tParam.dP = m_dParameters[i];
			for(int j = 0; j < (int)pCurve->NumBoxes(); j++)
			{
				CSBB* box = pCurve->GetBoxes()[j];
				if(box->m_fP1 <= tParam.dP && box->m_fP2 >= tParam.dP)
				{
					tParam.dPStart = box->m_fP1;
					tParam.dPEnd   = box->m_fP2;
					break;
				}//end if
				else if(box->m_fP2 <= tParam.dP && box->m_fP1 >= tParam.dP)
				{
					tParam.dPStart = box->m_fP2;
					tParam.dPEnd   = box->m_fP1;
					break;
				}//end if
   
				
			}//end for
			m_dParamBuckets[n].push_back(tParam);
			
		}//end if
	}//end for i
	
	pNode->SetData(nIndex);
	nIndex++;
	
}//end AssignPointsAndParameters

void CApproxCurve::AssignPoints(CBCNode *&pNode, int n)
{

	for(int i = 0; i < m_nResolution; i++)
	{
		VECTOR2 vV1 = VECTOR2::createVector(pNode->GetCenter(), m_vSamples[i]);
		double dist   = vV1.mag();
		if(dist <= pNode->GetRad())
		{
			m_Buckets[n].push_back(i);
			
		}//end if
	}//end for
	
	pNode->SetData(nIndex);

	nIndex++;
}//end AssignPoints

void CApproxCurve::ComputeCog()
{

	m_vCenter.x = 0.0;
	m_vCenter.y = 0.0;

	for(int i = 0; i < m_nResolution; i++)
		m_vCenter = m_vCenter + m_vSamples[i];

	m_vCenter = m_vCenter * (1.0/double(m_nResolution));

}//end ComputeCog

/////////////////////////////////////////////////////////
/*
*
*	The curve is rotatet around its center by a rotation 
*	matrix. The rotation matrix is passed by a parameter
*
*
*/
/////////////////////////////////////////////////////////

void CApproxCurve::RotateCurve(MATRIX2X2 &rotMat)
{

	Real yB =  std::numeric_limits<Real>::max();
	Real yT = -std::numeric_limits<Real>::max();
	Real xL =  std::numeric_limits<Real>::max();
	Real xR = -std::numeric_limits<Real>::max();

	VECTOR2 vTrans = m_vCenter * -1.0;

	for(int i = 0; i < m_nResolution; i++)
	{
		m_vSamples[i] = m_vSamples[i] + vTrans;
		m_vSamples[i] = rotMat * m_vSamples[i];
		m_vSamples[i] = m_vSamples[i] + m_vCenter;

		if(yB > m_vSamples[i].y)
		{
			yB = m_vSamples[i].y;
		}
		if(yT < m_vSamples[i].y)
		{
			yT = m_vSamples[i].y;
		}
		if(xL > m_vSamples[i].x)
		{
			xL = m_vSamples[i].x;
		}
		if(xR < m_vSamples[i].x)
		{
			xR = m_vSamples[i].x;
		}

	}//end for

	/* new bounding box vertices */
	m_bBox->m_Vertices[0].x = xL;
	m_bBox->m_Vertices[0].y = yB;
	m_bBox->m_Vertices[1].x = xR;
	m_bBox->m_Vertices[1].y = yB;
	m_bBox->m_Vertices[2].x = xR; 
	m_bBox->m_Vertices[2].y = yT;
	m_bBox->m_Vertices[3].x = xL;
	m_bBox->m_Vertices[3].y = yT;

	
	for(int j = 0; j < m_nCircles; j++)
	{
		RotateTree(m_BCTree.GetChild(j), rotMat);
	}//end for

}//end RotateCurve

//////////////////////////////////////////////////////
/*
*
*
*	Translates a curve by the vector rhs
*
*
*/
//////////////////////////////////////////////////////

void CApproxCurve::TranslateCurveByVector(const VECTOR2 &rhs)
{

	m_vCenter.x = 0.0;
	m_vCenter.y = 0.0;

	Real yB =  std::numeric_limits<Real>::max();
	Real yT = -std::numeric_limits<Real>::max();
	Real xL =  std::numeric_limits<Real>::max();
	Real xR = -std::numeric_limits<Real>::max();

	int i;
	for(i = 0; i < m_nResolution; i++)
	{
		/* translate samples */
		m_vSamples[i] = m_vSamples[i] + rhs;

		/* calculate new center */
		m_vCenter = m_vCenter + m_vSamples[i];

		if(yB > m_vSamples[i].y)
		{
			yB = m_vSamples[i].y;
		}
		if(yT < m_vSamples[i].y)
		{
			yT = m_vSamples[i].y;
		}
		if(xL > m_vSamples[i].x)
		{
			xL = m_vSamples[i].x;
		}
		if(xR < m_vSamples[i].x)
		{
			xR = m_vSamples[i].x;
		}

	}//end for

	/* new center */
	m_vCenter = m_vCenter * (1.0/static_cast<double>(m_nResolution));

	/* new bounding box vertices */
	m_bBox->m_Vertices[0].x = xL;
	m_bBox->m_Vertices[0].y = yB;
	m_bBox->m_Vertices[1].x = xR;
	m_bBox->m_Vertices[1].y = yB;
	m_bBox->m_Vertices[2].x = xR; 
	m_bBox->m_Vertices[2].y = yT;
	m_bBox->m_Vertices[3].x = xL;
	m_bBox->m_Vertices[3].y = yT;

	for(i = 0; i < m_nCircles; i++)
	{
		TranslateTree(m_BCTree.GetChild(i), rhs);
	}//end for

}//end TranslateCurve

void CApproxCurve::TransformCurve(MATRIX2X2 &rotMat, VECTOR2 &vT)
{
	Real yB = std::numeric_limits<Real>::max();
	Real yT = -std::numeric_limits<Real>::max();
	Real xL = std::numeric_limits<Real>::max();
	Real xR = -std::numeric_limits<Real>::max();

	VECTOR2 vTrans = m_vCenter * -1.0;

	for(int i = 0; i < m_nResolution; i++)
	{
		//translate back to origin
		m_vSamples[i] = m_vSamples[i] + vTrans;
		m_vSamples[i] = rotMat * m_vSamples[i];
		m_vSamples[i] = m_vSamples[i] - vTrans + vT;

		/* find new extreme values */
		if(yB > m_vSamples[i].y)
		{
			yB = m_vSamples[i].y;
		}
		if(yT < m_vSamples[i].y)
		{
			yT = m_vSamples[i].y;
		}
		if(xL > m_vSamples[i].x)
		{
			xL = m_vSamples[i].x;
		}
		if(xR < m_vSamples[i].x)
		{
			xR = m_vSamples[i].x;
		}

	}//end for

	/* calculate new center */
	m_vCenter = m_vCenter + vT;

	m_Root->SetCenter(m_vCenter);

	/* new bounding box vertices */
	m_bBox->m_Vertices[0].x = xL;
	m_bBox->m_Vertices[0].y = yB;
	m_bBox->m_Vertices[1].x = xR;
	m_bBox->m_Vertices[1].y = yB;
	m_bBox->m_Vertices[2].x = xR; 
	m_bBox->m_Vertices[2].y = yT;
	m_bBox->m_Vertices[3].x = xL;
	m_bBox->m_Vertices[3].y = yT;

	
	for(int j = 0; j < m_nCircles; j++)
	{
		TransformTree(m_BCTree.GetChild(j), rotMat, vTrans, vT);
	}//end for

}//end TransformCurve

void CApproxCurve::TranslateTree(CBCNode *pNode, const VECTOR2 &vTrans)
{

	if(pNode != NULL)
	{
		TranslateTree(pNode->m_Children[0], vTrans);
		pNode->Translate(vTrans);
		TranslateTree(pNode->m_Children[1], vTrans);
	}//end if

}//end inOrder

void CApproxCurve::TransformTree(CBCNode *pNode, MATRIX2X2 &rotMat, VECTOR2 &vTrans1, VECTOR2 &vTrans2)
{
	if(pNode != NULL)
	{
		TransformTree(pNode->m_Children[0], rotMat, vTrans1, vTrans2);
		pNode->Translate(vTrans1);
		pNode->Rotate(rotMat);
		VECTOR2 vV1 = vTrans2 - vTrans1;
		pNode->Translate(vV1);
		TransformTree(pNode->m_Children[1], rotMat, vTrans1, vTrans2);
	}//end if
}//end inOrder

void CApproxCurve::RotateTree(CBCNode *pNode, MATRIX2X2 &rotMat)
{

	if(pNode != NULL)
	{
		RotateTree(pNode->m_Children[0], rotMat);
		VECTOR2 vTrans = m_vCenter * -1.0;
		pNode->Translate(vTrans);
		pNode->Rotate(rotMat);
		pNode->Translate(m_vCenter);
		RotateTree(pNode->m_Children[1], rotMat);
	}//end if

}//end inOrder

void CApproxCurve::TranslateCurve(const VECTOR2 &rhs)
{

	for(int i = 0; i < m_nResolution; i++)
	{
		m_vSamples[i] = m_vSamples[i] + rhs;
	}//end for

}//end TranslateCurve

//////////////////////////////////////////////////
/*
*
*	Checks whether the point vQuery in inside
*	the curve. The method is based on the
*	jordan curve theorem
*
*
*/
//////////////////////////////////////////////////


bool CApproxCurve::IsInElement(const VECTOR2 &vQuery)
{

	//the number of intersections
	int nIntersections = 0;

	//loop over all line segments
	for(int i = 1; i < m_nResolution; i++)
	{
		//check if vQuery.x is within the x-Interval of the line segment
		if( (m_vSamples[i-1].x < vQuery.x && vQuery.x < m_vSamples[i].x) ||
			 (m_vSamples[i-1].x > vQuery.x && vQuery.x > m_vSamples[i].x ) )
		{
			//calculate the parameter value of the intersection point
			Real dT = (vQuery.x - m_vSamples[i].x) / (m_vSamples[i-1].x - m_vSamples[i].x);

			//calculate the y coordinate of the intersection point
			Real dY = dT * m_vSamples[i-1].y + (1-dT) * m_vSamples[i].y;
			
			//we are shooting down, so check if vQuery.y is above the point of intersection
			if(vQuery.y >= dY)
				nIntersections++;
		}//end if

	}//end for

	/*
	* if the number of intersections is even, then vQuery is outside the polygon
	* if the number of intersections is odd, vQuery is inside the polygon
	*/

	if(nIntersections % 2 == 0)
		return false;
	else
		return true;

}//end IsInElement

//////////////////////////////////////////////////
/*
*
*	Checks whether the point vQuery in inside
*	the curve. The method is based on the
*	jordan curve theorem
*
*
*/
//////////////////////////////////////////////////

bool CApproxCurve::Inside( const VECTOR2 &vQuery )
{
	
	int nIntersections = 0;	
	
	for(int j = 0; j < (int)m_Leaves.size(); j++)
	{
		CBCNode *pNode = m_Leaves[j];
		double dRad    = pNode->GetRad();
		VECTOR2 vMid = pNode->GetCenter();
		
		if( ((vQuery.x > (vMid.x + dRad) ) || (vQuery.x < (vMid.x - dRad) ) || (vQuery.y > (vMid.y + dRad) )) )
			continue;
		
		/* compute closest point in this leaf segment box */
		int nData = pNode->GetData();
	
		int nNumP = (int)m_Buckets[nData].size();	
		
		for(int i = 1; i < nNumP; i++)
		{
			if( (GetPoint(nData,i-1).x < vQuery.x && vQuery.x < GetPoint(nData,i).x) ||
						  (GetPoint(nData,i-1).x > vQuery.x && vQuery.x > GetPoint(nData,i).x ) )
			{
				double dT = (vQuery.x - GetPoint(nData,i).x) / (GetPoint(nData,i-1).x - GetPoint(nData,i).x);
				double dY = dT * GetPoint(nData,i-1).y + (1-dT) * GetPoint(nData,i).y;

				if(vQuery.y >= dY)
					nIntersections++;
			}//end if

		}//end for
		
	}//end for 
	
	if(nIntersections % 2 == 0)
		return false;
	else
		return true;		

}//end Inside

void CApproxCurve::OutputToFile(char *fileName)
{

	/* create an ofstream object */
	ofstream myFile(fileName,ios::out);
    myFile.precision(16); 

	/* check if file is valid */
    if (! myFile) 
    {
		cout << "Error opening output file: "<<fileName << endl;
        return;
    }

	int i;
	int num = m_nResolution;

	/* write number of points */
	myFile <<num<<"  2 0 0"<<endl;

	for(i = 1; i <= num ; i++)
	{
		myFile <<(i)<<" "<<m_vSamples[i-1].x<<" "<< m_vSamples[i-1].y<<endl;
		//myFile << 0.0 <<endl;
	}//end for i

	/* close file */
    myFile.close();

}//end OutputToFile

void CApproxCurve::Dump(char *sFileName)
{

	ofstream out(sFileName, ios::out);
	out.precision(16);

	int i = 0;
	int j = 0;
	int k = 0;

	if(!out.is_open())
	{
		cerr <<" Error opening "<<sFileName<<endl;
	}///end if

	out<< this->m_nResolution<<endl;
	for(i = 0; i < m_nResolution;i++)
	{
		out<< m_vSamples[i].x<< " " <<m_vSamples[i].y<<endl;
	}//end for
	out<< this->m_nCircles<<endl;
	out<< this->m_nRadius<<endl;
	out<< this->m_vCenter.x<<endl;
	out<< this->m_vCenter.y<<endl;

	/* write out bounding box */
	out << this->m_bBox->m_Vertices[0].x << " " << this->m_bBox->m_Vertices[0].y <<endl;
	out << this->m_bBox->m_Vertices[2].x << " " << this->m_bBox->m_Vertices[2].y <<endl;

	int nNodes = 0;

	int nBranch = 2;

	for(i = 0; i <= (int)TARGETDEPTH; i++)
	{
		nNodes += static_cast<int>(pow((double)nBranch,i));
	}//end for
	
	/* write number of nodes */
	out << nNodes <<endl;

	/* write out each circle tree */
	for(i = 0; i < m_nCircles; i++)
	{
		/* number of tree */
		out << i <<endl;
		queue<CBCNode*>  qNodes;

		qNodes.push(m_BCTree.GetChild(0));

		while(!qNodes.empty())
		{
			CBCNode * pNode = qNodes.front();
			qNodes.pop();
			out << pNode->GetRad() <<endl;
			out << pNode->GetCenter().x <<" "<<pNode->GetCenter().y <<endl;
			out << pNode->GetUpperVec().x <<" "<<pNode->GetUpperVec().y<<endl;
			out << pNode->GetID() <<endl;
			out << pNode->IsLeaf()<<endl;
			out << pNode->GetData()<<endl;
			
			if(!pNode->IsLeaf())
			{
				qNodes.push(pNode->m_Children[0]);
				qNodes.push(pNode->m_Children[1]);
			}//end if

		}//end while

	}//end for

}//end Dump

void CApproxCurve::ApproximateMedialAxis(vector<VECTOR2> &vPoints)
{

	int i;

	for(i = 0; i < (int)vPoints.size(); i++)
	{
		tVertex Vertex;
		Vertex.vVertex = vPoints[i];
		m_vVoronoi.push_back(Vertex);
	}//end for

	cout<<"Number of Voronoi vertices: "<<m_vVoronoi.size()<<endl;

	int iIndex = 0;
	for(i = 0; i < (int)m_vVoronoi.size(); i++)
	{
		
		if(IsInElement(m_vVoronoi[i].vVertex))
		{
			m_vVoronoi[i].inside = true;
			CMedAxisVertex *maVertex = new CMedAxisVertex();
			//m_vVoronoi[i].vVertex);
			maVertex->m_Circle.SetCenter(m_vVoronoi[i].vVertex);
			maVertex->m_Index = i;
			m_MedAxis.push_back(maVertex);
			m_vVoronoi[i].index = iIndex;
			iIndex++;
		}//end if
		else
		{
			m_vVoronoi[i].inside = false;
			m_vVoronoi[i].index = -2;
		}
	}//end for

	cout<<"Medial Axis points: "<<m_MedAxis.size()<<endl;

	vPoints.clear();

}//end ApproximateMedialAxis

void CApproxCurve::ReadEdges(const char* sFileName)
{
	/* open file for reading */
	ifstream fin(sFileName);

	/* quick error check */
	if(!fin.is_open())
	{
		cout<<"error opening file "<<sFileName<<endl;
		exit(0);
	}

	/* stores the number of edges and the number of boundary markers */
	int parameters[2];

	/* count variable */
	int i;
		
	/* read the above parameters */
	for(i = 0; i < 2; i++)
	{
		fin >> parameters[i];
	}//end for

	/* read in all edges */
	for(i = 0; i < parameters[0]; i++)
	{ 
		int num;
		int iStart;
		int iEnd;
		double dX,dY;

		fin >> num;
		fin >> iStart;
		fin >> iEnd;

		/* skip the special case */
		if(iEnd == -1)
		{
			fin >> dX;
			fin >> dY;
			continue;
		}

		iStart -=1;
		iEnd   -=1;

		/* test if both endpoints of the edge are inside */
		if(m_vVoronoi[iStart].inside && m_vVoronoi[iEnd].inside)
		{
			tEdge edge;
			edge.e1 = iStart;
			edge.e2 = iEnd;
			m_vEdges.push_back(edge);		
		}
		
	}//for
	
	fin.close();
		
	cout << "Number of all edges: "<<parameters[0]<<endl;
	cout << "Number of edges between inner vertices: "<<m_vEdges.size()<<endl;

}//end ReadEdges

void CApproxCurve::FindFormingPoints()
{
	
	int i,j;
	
	cout<<"Calculation radii...."<<endl;
	
	std::list<VertexIndex> lPoints;
	std::list<VertexIndex>::iterator lIter;
	std::list<VertexIndex>::iterator lIter0;

	vector<VECTOR2> OrigPoints;

	/* if the voronoi diagramm resolution is not equal to the curve resolution */
	readObjectVertices("Stern.node",OrigPoints);

	/* allocate memory for the leaf nodes */
	m_Radii = new double[m_MedAxis.size()];

		
	for(j = 0; j < (int) this->m_MedAxis.size(); j++)
	{
							
		//tFPoints fPoints;
		VertexIndex vClosest;
				
		lPoints.clear();
			
		for(i = 0; i < m_nResolution; i++)
		{
			VertexIndex vi;
			vi.vertex = m_vSamples[i];
			vi.index  = i;
			lPoints.push_back(vi);		
		}//end for

		for(int k = 0; k < 3; k++)
		{		
		
			double d = 1.7E+308;
					
			for(lIter = lPoints.begin(); lIter != lPoints.end(); lIter++)
			{
				VertexIndex vv = *lIter;
				VECTOR2 vVec = vv.vertex;
				vVec = VECTOR2::createVector(vVec, m_MedAxis[j]->m_Circle.m_vCenter);
				double dist = vVec.mag();
				if(d > dist)
				{
					d = dist;
					lIter0 = lIter;
					vClosest = *lIter;
				}//end if
			}//end for
			
			m_MedAxis[j]->m_Vertices.push_back(&m_vSamples[vClosest.index]);
			//cout<<"Closest"<<vClosest;
			m_MedAxis[j]->m_Circle.m_dRad = d;
			lPoints.erase(lIter0);
			m_Radii[j] = d;
		
		}//end for

	}//end for
	
// 	for(i = 0; i < (int)m_MedAxis.size(); i++)
// 	{
/*		for(j=0; j < (int)m_MedAxis[i]->m_Vertices.size();j++)
		{*/
			//cout<<*m_MedAxis[i]->m_Vertices[j];
// 			cout<<m_MedAxis[i]->m_Vertices[j];
// 		cout<<i<<" radius "<<m_MedAxis[i]->m_Circle.m_dRad<<endl;
// 		}//end for
		
// 	}//end for

	
}//end findFormingPoints	




/*!
    \fn CApproxCurve::MergeVertices(myVertex &u, myVertex &v, DebugSEC &info)
 */
//void CApproxCurve::MergeVertices(myVertex &u, myVertex &v, DebugSEC &info)
//{
//	
//	std::set<VECTOR2*,comparePoints> vSet;
//	std::set<VECTOR2*,comparePoints>::iterator sIter;
//	
//	CMedAxisVertex *vV1; 
//	CMedAxisVertex *vV2;
//	
//	
//	//get the property accessor
//	MAProperty prMAV = get(vertex_MedAxis, undirGraph);
//	
//	vV1   = boost::get(prMAV, u);
//	vV2  = boost::get(prMAV, v);	
//	
//	
//	//Merge the vertex lists
//	int i = 0;
//	for(i = 0; i < (int)vV1->m_Vertices.size(); i++)
//	{
//		vSet.insert(vV1->m_Vertices[i]);
//	}//end for
//	
//	for(i = 0; i < (int)vV2->m_Vertices.size(); i++)
//	{
//		vSet.insert(vV2->m_Vertices[i]);
//	}//end for
//	
//	//compute the smallest enclosing circle of the vertices
//	cPrimOps pO;
//	
//	CCircle circle = pO.SEC(vSet, info);
//
//	cout<<"Radius of smallest circle "<<circle.GetRadius()<<endl;
//	cout<<"Center of new circle "<<circle.GetCenter()<<endl;
//
//	//add the new vertex
//	myVertex newVertex = add_vertex(undirGraph);
//	
//	//create a MedAxisVertex
//	CMedAxisVertex *maVertex = new CMedAxisVertex();
//	
//	//let m_testVertex = maVertex for the purpose of debugging
//	m_TestVertex = maVertex;
//			
//	//assign the new circle
//	maVertex->m_Circle = circle;
//	maVertex->m_Index = -1;
//
//	for(sIter = vSet.begin(); sIter != vSet.end(); sIter++)
//	{
//		VECTOR2 *pVertex = *sIter;
//		maVertex->m_Vertices.push_back(	pVertex );
//	}//end for
//
//	//assign the MedAxis property
//	prMAV[newVertex] = maVertex;
//
//	
//// 	cout<<"aus grad "<<out_degree(u, undirGraph)<<endl;
//// 	cout<<"in grad "<<in_degree(u, undirGraph)<<endl;
//
//	cout<<"before Number of edges "<<num_edges(undirGraph)<<endl;	
//	cout<<"Number of vertices "<<num_vertices(undirGraph)<<endl;
//
//	cout<<"Grad vertex1 : "<<out_degree(u, undirGraph)<<endl;
//
//	cout<<"Grad vertex2 : "<<out_degree(v, undirGraph)<<endl;
//	
//	
//	graph_traits<UndirectedGraph>::adjacency_iterator vIter, vIter_end;
//	
//	remove_edge(u,v,undirGraph);
//	
//	vIter     = adjacent_vertices(u, undirGraph).first;
//	vIter_end = adjacent_vertices(u, undirGraph).second;
//	
//	
//	graph_traits<UndirectedGraph>::edge_descriptor e;	
//	bool found;
//	
//	VECTOR2 vMid = maVertex->m_Circle.GetCenter();
//	
//	typedef UndirectedGraph::edge_property_type Weight;	
//	
//	
//	
//	for(;vIter != vIter_end; ++vIter)
//	{
//		myVertex adjVertex = *vIter;
//		CMedAxisVertex *temp   = prMAV[adjVertex];
//		VECTOR2 vMid2 = temp->m_Circle.GetCenter();
//		double length = (VECTOR2::createVector(vMid, vMid2)).mag();
//		cout<<"Added edge with length: "<<length<<endl;
//		cout<<"("<<vMid.x<<","<<vMid.y<<")--->"<<"("<<vMid2.x<<","<<vMid.y<<")"<<endl;
//		cout<<"================"<<endl;
//		add_edge(newVertex, adjVertex, Weight(length),  undirGraph);
//	}//end for	
//	
//	vIter     = adjacent_vertices(v, undirGraph).first;
//	vIter_end = adjacent_vertices(v, undirGraph).second;
//	
//	for(;vIter != vIter_end; ++vIter)
//	{
//		myVertex adjVertex = *vIter;
//		tie(e, found) = edge(newVertex,adjVertex, undirGraph);
//		if(!found)
//		{
//			CMedAxisVertex *temp   = prMAV[adjVertex];
//			VECTOR2 vMid2 = temp->m_Circle.GetCenter();
//			double length = (VECTOR2::createVector(vMid, vMid2)).mag();	
//			cout<<"Added edge with length: "<<length<<endl;
//			cout<<"("<<vMid.x<<","<<vMid.y<<")--->"<<"("<<vMid2.x<<","<<vMid.y<<")"<<endl;
//			cout<<"================"<<endl;
//			add_edge(newVertex, adjVertex, Weight(length),  undirGraph);			
//		}//end if	
//	}//end for	
//	
//	
//	//delete the old edges
// 	clear_vertex(u, undirGraph);
// 	clear_vertex(v, undirGraph);
//	
//	//delete the old vertices
//  	remove_vertex(v, undirGraph);		
// 	remove_vertex(u, undirGraph);	
//	
////  	//get edge iterator for the graph data structure
//// 	graph_traits<UndirectedGraph>::edge_iterator iE;
//// 	
//// 	//get a property accessor	
//// 	EWProperty prEW  = get(edge_weight, undirGraph);
//// 		
//// 	for(iE = edges(undirGraph).first; iE != edges(undirGraph).second; ++iE)
//// 	{
//// 		
//// 		myVertex vStart,vEnd;
//// 		vStart = source(*iE,undirGraph);
//// 		vEnd   = target(*iE,undirGraph);
//// 		CMedAxisVertex *start   = prMAV[vStart];
//// 		CMedAxisVertex *end     = prMAV[vEnd];
//// 		VECTOR2 vMid3 = start->m_Circle.GetCenter();
//// 		VECTOR2 vMid4 = end->m_Circle.GetCenter();
//// 		prEW[*iE] = (VECTOR2::createVector(vMid3, vMid4)).mag();	
//// 
//// 	}//end for	
//	
//	cout<<"after Number of edges "<<num_edges(undirGraph)<<endl;
//	cout<<"Number of vertices "<<num_vertices(undirGraph)<<endl;	
//	
//	//for(int j = 0; j < GradOne.size(); j++)
//	//{
//	//	cout<<"grad: "<<out_degree(GradOne[j], undirGraph)<<endl;
//	//}//end for
//
//	//free the old MedAxisVertices
//	delete vV1;
//	delete vV2;
//	
//}//end MergeVertices
//
//




/*!
    \fn CApproxCurve::BuildBoostGraph()
 */
//template<typename UndirectedGraph> 
//void CApproxCurve::BuildBoostGraph()
//{
//
//	using namespace boost;
//	
//	//stores the inner edges
//	std::vector<tEdge> vInnerEdges;
//	
//	//find the edges that are between inner vertices
//	MakeInnerEdges(vInnerEdges);
//	
//	//number of inner vertices
//	int V = this->m_MedAxis.size();
//		
//	//initialize the graph data structure
//	undirGraph = UndirectedGraph(V);
//		
//	//get the property accessor
//	MAProperty prMAV = get(vertex_MedAxis, undirGraph);
//	
//	EWProperty prEW  = get(edge_weight, undirGraph);
//	
//	typedef typename property_traits<MAProperty>::value_type MAType;
//		
//	myVertex *v = new myVertex[this->m_MedAxis.size()];
//
//	//assign the MedAxis property to the vertices
//	int i = 0;
//	typename graph_traits<UndirectedGraph>::vertex_iterator vi;
//	for (vi = vertices(undirGraph).first; vi != vertices(undirGraph).second; ++vi)
//	{
//		v[i] = *vi;
//		prMAV[v[i]] = m_MedAxis[i];
//		i++;
//	}//end for
//	
//	//get number of inner edges	
// 	int nEdges = (int)vInnerEdges.size();
//	
//	typename graph_traits<UndirectedGraph>::edge_descriptor e;	
// 	
//	typedef typename UndirectedGraph::edge_property_type Weight;			
//	
// 	//add the inner edges to the adjacency list
//	for(i = 0; i < nEdges; i++)
//	{
//		//add_edge(v[vInnerEdges[i].e1], v[vInnerEdges[i].e2], undirGraph);
//		add_edge(v[vInnerEdges[i].e1], v[vInnerEdges[i].e2], Weight(vInnerEdges[i].length), undirGraph);
//		
//	}//end for
//	
//	cout << "vertices(g) = "<<num_vertices(undirGraph)<<endl;
//
//	int grad = 0;
//	for (vi = vertices(undirGraph).first; vi != vertices(undirGraph).second; ++vi)
//	{
//		if(out_degree(*vi, undirGraph) == 1)
//		{
//			grad++;
//			GradOne.push_back(*vi);
//		}
//	}//end for
//
//	//cout<<"# grad 1 "<<grad<<endl;
//	
//	//now calculate the edge weights
//	
//// 	//get edge iterator for the graph data structure
//// 	typename graph_traits<UndirectedGraph>::edge_iterator iE;
//// 	
//// 	for(iE = edges(undirGraph).first; iE != edges(undirGraph).second; ++iE)
//// 	{
//// 		cout<<"weight "<<prEW[*iE]<<endl;
//// 	}//end for
//	
//	//free memory for the vertex descriptors
//	delete[] v;
//	
//}//end build boost graph
//

/*!
    \fn CApproxCurve::MakeInnerEdges(vector<tEdge>	&InnerEdges)
 */
void CApproxCurve::MakeInnerEdges(vector<tEdge>	&InnerEdges)
{
	
    for(int i = 0; i < m_vEdges.size(); i++)
	{
		
		for(int k = 0; k < m_MedAxis.size(); k++)
		{
			CMedAxisVertex *pVertex = m_MedAxis[k];
			if(m_vEdges[i].e1 == pVertex->m_Index)
			{
				for(int j = k; j < m_MedAxis.size(); j++)
				{
					CMedAxisVertex *pVertex2 = m_MedAxis[j];
					if(m_vEdges[i].e2 == pVertex2->m_Index)
					{
						tEdge anEdge;
						anEdge.e1 = k;
						anEdge.e2 = j;
						VECTOR2 start = pVertex->m_Circle.GetCenter();
						VECTOR2 end   = pVertex2->m_Circle.GetCenter();
						anEdge.length = (VECTOR2::createVector(start, end)).mag();
						InnerEdges.push_back(anEdge);	
					}//end if
				}
			}
			
		}//end for
	}//end for
}///end MakeInnerEdges


/*!
    \fn CApproxCurve::DetermineEdge(DebugSEC &info)
 */
//void CApproxCurve::DetermineEdge(DebugSEC &info)
//{
//	
//	//geet a property accessor
//	EWProperty prEW  = get(edge_weight, undirGraph);
//	
//	//get edge iterator for the graph data structure
//	graph_traits<UndirectedGraph>::edge_iterator iE;
//		
//	//store the minimum weight
//	double dMinWeight = 1.7E+308;
//	
//	
//	
//	//vertices that make up the edge
//	myVertex u,v;
//	
//	//determine shortest edge
//	for(iE = edges(undirGraph).first; iE != edges(undirGraph).second; ++iE)
//	{
//
//		myVertex start = source(*iE, undirGraph);
//		myVertex end = target(*iE, undirGraph);
//
//		if( (out_degree(start, undirGraph) == 1) || (out_degree(end, undirGraph) == 1) )
//			continue;
//
//		double dWeight = prEW[*iE];
//// 		cout<<"weight "<<prEW[*iE]<<endl;
//		if(dWeight < dMinWeight)
//		{
//			dMinWeight = dWeight;
//			u = source(*iE, undirGraph);
//			v = target(*iE, undirGraph);
//		}//end if
//	}//end for*/
//	
//// 	iE = edges(undirGraph).first;	
//// 	u = source(*iE, undirGraph);
//// 	v = target(*iE, undirGraph);	
//	
//	//get the property accessor
//	MAProperty prMAV = get(vertex_MedAxis, undirGraph);	
//	
///*	cout<<"Minimum edge length: "<<dMinWeight<<endl;
//	cout<<"Removing edge: "<<prMAV[u]->m_Circle.GetCenter() <<endl;*/
//	
//	MergeVertices(u,v, info);
//	
//}//end Determine Edge


/*!
    \fn CApproxCurve::GetMergedCircle(myVertex &u, myVertex &v,DebugSEC &info)
 */
//CCircle CApproxCurve::GetMergedCircle(myVertex &u, myVertex &v, DebugSEC &info)
//{
//	
//	std::set<VECTOR2*,comparePoints> vSet;
//	std::set<VECTOR2*,comparePoints>::iterator sIter;
//	
//	CMedAxisVertex *vV1; 
//	CMedAxisVertex *vV2;
//	
//	//get the property accessor
//	MAProperty prMAV = get(vertex_MedAxis, undirGraph);
//	
//	vV1   = boost::get(prMAV, u);
//	vV2  = boost::get(prMAV, v);	
//	
//	
//	//Merge the vertex lists
//	int i = 0;
//	for(i = 0; i < (int)vV1->m_Vertices.size(); i++)
//	{
//		vSet.insert(vV1->m_Vertices[i]);
//	}//end for
//	
//	for(i = 0; i < (int)vV2->m_Vertices.size(); i++)
//	{
//		vSet.insert(vV2->m_Vertices[i]);
//	}//end for
//	
//	//compute the smallest enclosing circle of the vertices
//	cPrimOps pO;
//	
//	CCircle circle = pO.SEC(vSet, info);
//	
//   return circle;
//}


/*!
    \fn CApproxCurve::DetermineOptimalEdge(DebugSEC &info)
 */
//void CApproxCurve::DetermineOptimalEdge(DebugSEC &info)
//{

//	
//	//get edge iterator for the graph data structure
//	graph_traits<UndirectedGraph>::edge_iterator iE;
//		
//	//store the minimum weight
//	double dMinRadius = 1.7E+308;
//	
//	
//	
//	//vertices that make up the edge
//	myVertex u,v;
//	
//	//determine shortest edge
//	for(iE = edges(undirGraph).first; iE != edges(undirGraph).second; ++iE)
//	{
//
//		myVertex start = source(*iE, undirGraph);
//		myVertex end = target(*iE, undirGraph);
//
//		if( (out_degree(start, undirGraph) == 1) || (out_degree(end, undirGraph) == 1) )
//			continue;
//
//		CCircle daCircle = GetMergedCircle(start, end, info);
//		
//		double dRadius =  daCircle.GetRadius();
//
//					
//		if(dRadius < dMinRadius)
//		{
//			dMinRadius = dRadius;
//			u = source(*iE, undirGraph);
//			v = target(*iE, undirGraph);
//		}//end if
//	}//end for
//	
//	info.outcome = 0;
//	
//	MergeVertices(u,v, info);	
//	
//}//DetermineOptimalEdge
//
//

void CApproxCurve::RotateCurve(Real nAngle)
{

	//set up a Rotation matrix
	MATRIX2X2 matRotation;
	matRotation.InitRotationMatrix(nAngle);

	Real yB =  std::numeric_limits<Real>::max();
	Real yT = -std::numeric_limits<Real>::max();
	Real xL =  std::numeric_limits<Real>::max();
	Real xR = -std::numeric_limits<Real>::max();

	VECTOR2 vTrans = m_vCenter * -1.0;

	for(int i = 0; i < m_nResolution; i++)
	{
		m_vSamples[i] = m_vSamples[i] + vTrans;
		m_vSamples[i] = matRotation * m_vSamples[i];
		m_vSamples[i] = m_vSamples[i] + m_vCenter;

		if(yB > m_vSamples[i].y)
		{
			yB = m_vSamples[i].y;
		}
		if(yT < m_vSamples[i].y)
		{
			yT = m_vSamples[i].y;
		}
		if(xL > m_vSamples[i].x)
		{
			xL = m_vSamples[i].x;
		}
		if(xR < m_vSamples[i].x)
		{
			xR = m_vSamples[i].x;
		}

	}//end for

	/* new bounding box vertices */
	m_bBox->m_Vertices[0].x = xL;
	m_bBox->m_Vertices[0].y = yB;
	m_bBox->m_Vertices[1].x = xR;
	m_bBox->m_Vertices[1].y = yB;
	m_bBox->m_Vertices[2].x = xR; 
	m_bBox->m_Vertices[2].y = yT;
	m_bBox->m_Vertices[3].x = xL;
	m_bBox->m_Vertices[3].y = yT;

	
	for(int j = 0; j < m_nCircles; j++)
	{
		RotateTree(m_BCTree.GetChild(j), matRotation);
	}//end for



}//end RotateCurve

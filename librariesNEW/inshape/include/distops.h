//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : InShape2D
//  @ File Name : DistOps.h
//  @ Date : 22.06.2006
//  @ Author : Raphael Mnster
//
//


#if !defined(_DISTOPS_H)
#define _DISTOPS_H

#include <vector>
#include <list>
#include <stdio.h>
#include <stdarg.h>

#include "geomops.h"
#include "nurbs.h"
#include "grid.h"

///===============================================
///					DEFINES
///===============================================

#define FILENAME "log.txt"

#define _SINGLETHREAD_

#ifndef WIN32
#define LINUX
#endif

#define MAX_THREADS 2

#define SIMD_SHUFFLE(srch,src1,desth,dest1) ( (srch << 6) | (src1 << 4) | (desth << 2) | (dest1) )

///===============================================
///				FORWARD DECLARATIONS	
///===============================================


extern void log(const char *strFileName ,char *format, ...);

class CApproxCurve;
class CAABBTree;
class CAABBNode;
class CBCNode;
class CBCTree;
class CAABB;


//===============================================
//
//		Declaration of thread functions in case
//		a multithreaded distance computation is
//		desired.
//
//===============================================

#if defined WIN32
void threadProc1();
void threadProc2();
#endif

#if defined LINUX
void *threadProc1(void *p1);
void *threadProc2(void *p2);
void *threadProcNewton1(void *p1);
void *threadProcNewton2(void *p2);
#endif

using namespace std;


///===============================================
///	Typedefs to provide uncomplicated type names
///===============================================

typedef list<CAABBNode*> ListAABBND;
typedef vector<CApproxCurve*> CurveVec;
typedef list<CBCNode*> ListBCND;
typedef VECTOR2 vGrid[GRIDSIZEX][GRIDSIZEY];
typedef vector<CNurbs*> NURBSVec;
typedef vector<CONXNurbs*> ONXVec;

///===============================================
///					STRUCTS
///===============================================

typedef struct
{
	double* da;
	std::vector<CSBB*> boxes; 
}distInf;

typedef struct
{
	Real res;
	int  data;
	int  point;
	Real dLUB;
	CBCNode *pBest;
}SearchResult;

typedef struct
{
	Real res;
	int data;
	int point;
	VECTOR2  vUBounds[TARGETDEPTH+1];
	ListBCND candidates[TARGETDEPTH+1];
	CBCNode* pBest[TARGETDEPTH+1];
	CBCNode* pBestPoint;
	VECTOR2 vPoint;
	VECTOR2 vPointOCurve;
}DebugResult;

//===============================================
//            FUNCTION DECLARATION
//===============================================
void InitResources(CGrid *grid, CAABBTree *pTree, CurveVec &ApproxCurves, ONXVec& rNurbs);
void DeInitResources();


//==========================================================================
//															
//  This class contains the minimum distance computation functions.
//  The function boxSearch is the main algorithm which computes the distance   
//  between an arbitrary point and a NURBS curve. The functions distPointBox  
//  and ClosestPointOnLine are helper functions in the process of finding     
//  the curve/point distance										     
//															      
//==========================================================================

class DistOps : public GeomOps
{
public:
	VECTOR2	  ClosestPointOnLine(VECTOR2 vA, VECTOR2 vB, VECTOR2 vPoint);
	double    distPointBox(VECTOR2 vPoint, CAABB *box);
	bool      numIntersections(CBSpline* pncurve, const VECTOR2& pvPoint);
	int		  checkSearchResults(double *dRes1, double *dRes2, const VECTOR2 &vPoint);
	//Strange old stuff... potentially dangerous, to be watched by Sch�uble...

	//These functions are better, no dangerous return parameters. get funky Wolfgang...
	SearchResult   bruteDistance(CApproxCurve *pCurve, const VECTOR2 &vQuery);

	SearchResult	  LUBSearch(CApproxCurve *Curve, const VECTOR2 &pvPoint, ListBCND &lCandidates, CBCNode *&pBest);

	void	  breadthFirstSearchDistance(CApproxCurve* pCurve, const VECTOR2 &vQuery, ListBCND &lBFS, CBCNode *&nBest);

	SearchResult   LUBSearch(ListAABBND &Candidates, const VECTOR2 &pvPoint, CApproxCurve *&nCurve, CurveVec &ApproxCurves);

	void      breadthFirstSearchDistance(ListAABBND &Candidates, const VECTOR2 &vQuery, ListBCND &lBFS, CBCNode *&nBest);
	SearchResult AABBTreeSearch(CAABBTree *pTree, const VECTOR2 &vPoint, CApproxCurve *&nCurve, CurveVec &ApproxCurves);

	SearchResult COLUBSearch(ListAABBND &Candidates, const VECTOR2 &pvPoint, CApproxCurve *&nCurve,
				    CurveVec &ApproxCurves,  SearchResult &lRes);
	
	SearchResult   LUBSearch(ListAABBND &Candidates, const VECTOR2 &pvPoint, CApproxCurve *&nCurve,
			  CurveVec &ApproxCurves, CBCNode *&pBest, ListBCND &lBFS);

	SearchResult   COLUBSearch(CApproxCurve *pCurve, const VECTOR2 &pvPoint, ListBCND &lBFS, CBCNode *&pBest,double dLUB);

	DebugResult DebugSearch(CApproxCurve *Curve, const VECTOR2 &pvPoint, ListBCND &lCandidates, CBCNode *&pBest);

	CONXNurbs *myNurbs;

public:
	void BreadthFirstSearchDebug(CApproxCurve* pCurve, const VECTOR2 &vQuery, ListBCND &lBFS, CBCNode *&nBest, DebugResult &data);
public:
	SearchResult BruteDistance_SIMD(CApproxCurve *pCurve, const VECTOR2 &vQuery);
};

typedef struct
{
	vGrid     m_grid;
	CAABBTree *m_pTree;
	VECTOR2   &m_vPoint;
	CApproxCurve *m_nCurve;
	CurveVec  &m_ApproxCurves;
	NURBSVec  &m_bcurve;
}BenchmarkInfo;

typedef struct
{
	CGrid          *Grid;
	CurveVec       ApproxCurves;
	CAABBTree      *pTree;
	ONXVec         NURBS;
}t_ThreadResource;

#endif  //_DISTOPS_H

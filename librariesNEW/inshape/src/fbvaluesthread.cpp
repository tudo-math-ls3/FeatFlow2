#include "fbvaluesthread.h"
#include "approxcurve.h"
#include "feat2ops.h"
#include "grid.h"

CFBValuesThread::CFBValuesThread(int iStart, int iEnd, CApproxCurve *pCurve, CGrid *pGrid) : CBasicThread()
{

	m_pGrid  = pGrid;
	m_pCurve = pCurve;

	m_iStart = iStart;

	m_iEnd = iEnd;

}//end constructor

CFBValuesThread::~CFBValuesThread(void)
{

}//end deconstructor

void CFBValuesThread::CleanUp()
{

}//end CleanUp

void CFBValuesThread::ThreadProc()
{

	CFeat2Ops op;

	for(int i = m_iStart; i < m_iEnd; i++)
	{
		Real res = op.IsFictitiousBoundarySimple(m_pCurve, m_pGrid->m_pVertices[i]);
		m_pGrid->SetDistance(i, res);
	}//end for i

}//End ThreadProc

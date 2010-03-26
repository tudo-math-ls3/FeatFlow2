#pragma once
#include "basicthread.h"

class CGrid;
class CApproxCurve;

class CFBValuesThread :
	public CBasicThread
{
public:

	CFBValuesThread(int iStart, int iEnd, CApproxCurve *pCurve, CGrid *pGrid);
	~CFBValuesThread(void);

	void CleanUp();

	void ThreadProc();

private:

	int m_iStart;
	int m_iEnd;
	CApproxCurve *m_pCurve;
	CGrid        *m_pGrid;


};

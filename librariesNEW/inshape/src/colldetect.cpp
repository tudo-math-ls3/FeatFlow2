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

#include <colldetect.h>
#include <approxcurve.h>
#include <distops.h>

CCollDetect::CCollDetect()
{
}


CCollDetect::~CCollDetect()
{
}


/*!
    \fn CCollDetect::CollAllPairs(CurveVec_t &vCurves, CollList &lCollisions)
 */
void CCollDetect::CollAllPairs(CurveVec_t &vCurves, CollList &lCollisions)
{
	
	int i,j;
	
	//all pairs processing
	for(i = 0; i < (int)vCurves.size(); i++)
	{
		VECTOR2 vCenterI = vCurves[i]->GetCenter();
		Real    nRadI    = vCurves[i]->Radius();
		for(j = i+1; j < (int)vCurves.size(); j++)
		{
			VECTOR2 vCenterJ = vCurves[j]->GetCenter();
			Real    nRadJ    = vCurves[j]->Radius();
			
			Real nDist = VECTOR2::createVector(vCenterI,vCenterJ).mag();
			
			//process a single pair
			if(ProcessPairSimple(nRadI, nRadJ, nDist))
			{
				CollPair collPair(i,j);
				lCollisions.push_back(collPair);	
			}//end if
			
		}//end for
	}//end for i

}//end CollAllPairs


/*!
    \fn CCollDetect::TestBVHs(CApproxCurve *pCurve1,CApproxCurve *pCurve2 )
 */
bool CCollDetect::TestBVHs(CApproxCurve *pCurve1,CApproxCurve *pCurve2 )
{
	
	int i;
	
	VECTOR2 vCenterI = pCurve1->GetCenter();
	Real    dRad1    = pCurve1->Radius();


	VECTOR2 vCenterJ = pCurve2->GetCenter();
	Real    dRad2    = pCurve2->Radius();
	
	Real dDist = VECTOR2::createVector(vCenterI,vCenterJ).mag();	
		
	if(dDist > (dRad1+dRad2) )
		return false;
		
	Real nSmlRad = 0.0;
	VECTOR2 vCenter;
	CApproxCurve *smlCurve = NULL;
	CApproxCurve *lrgCurve = NULL;
	
	if(dRad1 > dRad2)
	{
		
		nSmlRad = dRad2;
		vCenter = vCenterJ;
		smlCurve = pCurve2;
		lrgCurve = pCurve1;
	}//end if
	else
	{
		nSmlRad = dRad1;
		vCenter = vCenterI;
		smlCurve = pCurve1;
		lrgCurve = pCurve2;
	}//end else
	
	int nChildren = lrgCurve->m_Root->NumChildren();
	
	if(DFSColl(smlCurve->m_Root, lrgCurve->m_Root->GetChildren(), nChildren))
		return true;
	else
		return false;

    
}//TestBVHs


/*!
    \fn CCollDetect::DFSColl(CBCNode *pBig, CBCNode *pTargets, int nTargets)
 */
bool CCollDetect::DFSColl(CBCNode *pBig, CBCNode **pTargets, int nTargets)
{
	
	Real nRadBig = pBig->GetRad();
	VECTOR2 vCenter = pBig->GetCenter();
	
	//test all Targets against pBig
	for(int i = 0; i < nTargets; i++)
	{
		Real nRad = pTargets[i]->GetRad();
		VECTOR2 vCenterTar = pTargets[i]->GetCenter();
		Real dDist = VECTOR2::createVector(vCenter,vCenterTar).mag();
		//If there is a collision 
		if(ProcessPairSimple(nRadBig, nRad, dDist))
		{

			if(pBig->IsLeaf() && pTargets[i]->IsLeaf())
			{
				return true;
			}
			else
			{
				if(DFSColl(pTargets[i], pBig->GetChildren(), pBig->NumChildren()))
					return true;
				else
					return false;
			}
		}//end if
	}//end for
		
	//no collision has occured so far, return false
    return false;
}//end DFSColl


/*!
    \fn CCollDetect::CollAllBVHs(CurveVec_t &vCurves, CollList &lCollisions)
 */
void CCollDetect::CollAllBVHs(CurveVec_t &vCurves, CollList &lCollisions)
{
	int i,j;
	
	//all pairs processing
	for(i = 0; i < (int)vCurves.size(); i++)
	{
		VECTOR2 vCenterI = vCurves[i]->GetCenter();
		Real    nRadI    = vCurves[i]->Radius();
		for(j = i+1; j < (int)vCurves.size(); j++)
		{
			VECTOR2 vCenterJ = vCurves[j]->GetCenter();
			Real    nRadJ    = vCurves[j]->Radius();
			
			Real nDist = VECTOR2::createVector(vCenterI,vCenterJ).mag();
			
			//process a single pair
			if(TestBVHs(vCurves[i], vCurves[j]))
			{
				CollPair collPair(i,j);
				lCollisions.push_back(collPair);	
			}//end if
			
		}//end for
	}//end for i
}

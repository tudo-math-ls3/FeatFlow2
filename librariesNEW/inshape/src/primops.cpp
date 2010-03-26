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

#include <math.h>
#include <primops.h>







CMedAxisVertex::CMedAxisVertex()
{
}

CMedAxisVertex::CMedAxisVertex(VECTOR2 vMid)
{
	this->m_Circle.SetCenter(vMid);
}//constructor


CMedAxisVertex::~CMedAxisVertex()
{
}


/*!
    \fn cPrimOps::SEC(set<VECTOR2*,comparePoints> &vSet)
    Rokne's algorithm to find the bounding circle for
	
 */
CCircler cPrimOps::SEC(set<VECTOR2*,comparePoints> &vSet, DebugSEC &info)
{
	double nRad = -1;

	int i;
	VECTOR2 vP;
	double minY = 1.7E+308;
	int iP;
	
	double nMinAngle = 1.7E+308;
	VECTOR2 vQ;
	int iQ;
	
	double nMinPRQ = 1.7E+308;
	VECTOR2 vR;
	int iR;
	
	double pi = 3.1415926535;
	
	vector<VECTOR2> vPts;
	set<VECTOR2*,comparePoints>::iterator sIter;
	
	for(sIter = vSet.begin(); sIter != vSet.end(); sIter++)
	{
		VECTOR2* pV = *sIter;
		VECTOR2 V = *pV;
		vPts.push_back(V);
		
	}//end for
	
	
	
	VECTOR2 vCenter;
	
	int numberOfFPs = (int)vPts.size();
	//================================================
	// first step determine a point P with minimum Py
	//================================================
	for(i = 0; i < (int)vPts.size(); i++)
	{
		if(vPts[i].y < minY)
		{
			vP = vPts[i];
			iP = i;
			minY = vPts[i].y;
		}//end if
	}//end for
	
	info.p = vP;
	
	// get norm of P 
	double nNorm = vP.mag();
	
	//================================================
	// second step find a point Q such that the 
	// angle between vector PQ and X-axis is minimal
	//================================================	
	for(i = 0; i < (int)vPts.size(); i++)
	{
		if(i == iP)
			continue;
		
		//the vector PPi
		VECTOR2 v0 = VECTOR2::createVector(vP, vPts[i]);
		
		//X-Axis
		VECTOR2 vXAxis(1,0);
		
		double nNorm2 = v0.mag();

		double nDot   = v0 * vXAxis;
		
		nDot/= nNorm2;
		
		double nAngle = acos(nDot);
		
		nAngle = nAngle * 180/pi;
		
	
		//set minimum angle and Point Q	
		if(nAngle < nMinAngle)
		{
			vQ = vPts[i];
			iQ = i;
			nMinAngle = nAngle;
		}//end if
		
	}//end for
	
	info.q = vQ;
	
	double nAngP = 100;
	double nAngQ = 100;
	
	//================================================
	// third step: find a point R such that die angle
	// PRQ is minimal
	//================================================	
	while(true)
	{	
		nMinPRQ = 1.7E+308;
		for(i = 0; i < (int)vPts.size(); i++)
		{
			if(i == iP || i == iQ)
				continue;
			
			//calculate angle PRQ
			
			//create vector RQ
			VECTOR2 vRQ = VECTOR2::createVector(vPts[i], vQ);
			//create vector RP
			VECTOR2 vRP = VECTOR2::createVector(vPts[i], vP);
			
			//get the angle
			double nAngle = vRQ.CalcAngle(vRP);
			
			if(nAngle < nMinPRQ)
			{
				nMinPRQ = nAngle;
				vR = vPts[i];
				iR = i;
			}//end if
							
		}//end for
		
		info.r = vR;
		
		VECTOR2 vPR = VECTOR2::createVector(vP,vR);
		VECTOR2 vPQ = VECTOR2::createVector(vP,vQ);
		VECTOR2 vQR = VECTOR2::createVector(vQ, vR);
		VECTOR2 vQP = VECTOR2::createVector(vQ, vP);
		
		
		//calculate angle at P		
		nAngP = vPR.CalcAngle(vPQ);
	
		//calculate angle at Q
		nAngQ = vQR.CalcAngle(vQP);
		
		
		//CASE 1: 
		// the angle at r is obtuse, the center of the circle is (p+q)/2 
		if(nMinPRQ > 90)
		{
			vCenter = vQ + vP;
			vCenter = vCenter * 0.5;
			VECTOR2 vVec1 = vP - vQ;
			vVec1 = vVec1 * 0.5;
			//set Radius of circle
			nRad = vVec1.mag();
			info.outcome = 0;
			break;
		}//end if
		//CASE 2:acute angle
		// angle at PQR < 90, the algorithm is finished... 
		else
		{
			
			if(nAngQ <= 90 && nAngP <= 90)
			{
				
				vPQ = VECTOR2::createVector(vP,vQ);
				// create middle senkrechte for PR 
				VECTOR2 vMPR = vP + (vPR * 0.5);
				VECTOR2 vNPR = VECTOR2::createVector(vP,vR);
				vNPR = VECTOR2::CreateNormal(vNPR);
				VECTOR2 vP2 = vMPR + vNPR;
				CParamLine P1(vMPR, vP2);
				


				// create middle senkrechte for PQ 
				VECTOR2 vMPQ = vP + (vPQ * 0.5);
				VECTOR2 vNPQ = VECTOR2::createVector(vP,vQ);
				vNPQ         = VECTOR2::CreateNormal(vNPQ);
				VECTOR2 vP1 = vMPQ - vNPQ;
				CParamLine P2(vMPQ, vP1);
				info.pq = vPQ;
				
				// create middle senkrechte for QR 
// 				VECTOR2 vMPQ = vQ + (vQR * 0.5);
// 				VECTOR2 vNPQ = VECTOR2::createVector(vQ,vR);
// 				vNPQ.Normalize();
// 				VECTOR2 vP1 = vMPQ - vNPQ;
// 				CParamLine P2(vMPQ, vP1);
								
				info.p1 = P1;
				info.p2 = P2;

				double t1,t2;
				
				//???error P1 and vCenter are the same point (possible???)
				if(CParamLine::IntersectParamLines2D(P1, P2, &t1, &t2))
				{
				
					vCenter = P1.Evaluate(t1);
					VECTOR2 vV0 = VECTOR2::createVector(vCenter, vP);
					//set Radius of circle
					nRad = vV0.mag();
					info.outcome = 1;
					break;
				}
				else
				{
					cout<<"There is some kind of error"<<endl;
					info.outcome = 2;
					break;
				}//
			}//end if
			// check whether the angle at p or q is obtuse and replace this point for r 
			else
			{
				if(nAngQ >= 90)
				{
					vQ = vR;
					iQ = iR;
					info.q = vR;
					
				}
				else if(nAngP >= 90)
				{
					vP = vR;
					iP = iR;
					info.p = vR;
				}
				
			}//end else
			
		}//end else
		
	}//end while

	CCircler circle(vCenter, nRad);
	
	return circle;
	
}

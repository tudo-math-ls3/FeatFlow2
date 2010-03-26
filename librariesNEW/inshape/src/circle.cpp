/***************************************************************************
 *   Copyright (C) 2007 by Raphael Münster   *
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

template<class T>
CCircle<T>::CCircle(const CVector2<T>& vCenter, T dRad) : m_vCenter(vCenter), m_dRad(dRad)
{
}//end constructor


/*!
    \fn CCircle::SmallestEnclosingCircle(vector<VECTOR2> &vPts)
 */
template<class T>
CCircle<T> CCircle<T>::SmallestEnclosingCircle(vector< CVector2<T> > &vPts)
{

	T nRad;

	int i;
	CVector2<T> vP;
	T minY = 1.7E+308;
	int iP;
	
	T nMinAngle = 1.7E+308;
	VECTOR2 vQ;
	int iQ;
	
	T nMinPRQ = 1.7E+308;
	CVector2<T> vR;
	int iR;
	
	T pi = 3.1415926535;
	
	CVector2<T> vCenter;
	
	/// find minimum Py (minimum y value)
	for(i = 0; i < (int)vPts.size(); i++)
	{
		if(vPts[i].y < minY)
		{
			vP = vPts[i];
			iP = i;
			minY = vPts[i].y;
		}//end if
	}//end for
	
//	cout << "Minimum Py "<<vP.x<< " " <<vP.y<<endl;
	
	/// get norm of Py 
	T nNorm = vP.mag();
	
	/// find point Q with minimum angle between vector PQ and X-axis 
	for(i = 0; i < (int)vPts.size(); i++)
	{
		if(i == iP)
			continue;
		
		CVector2<T> v0 = CVector2<T>::createVector(vP, vPts[i]);
		
		CVector2<T> vXAxis(1,0);
		
		T nNorm2 = v0.mag();

		T nDot   = v0 * vXAxis;
		
		nDot/= nNorm2;
		
		T nAngle = acos(nDot);
		
		nAngle = nAngle * 180/pi;
		
	
			
		if(nAngle < nMinAngle)
		{
			vQ = vPts[i];
			iQ = i;
			nMinAngle = nAngle;
		}//end if
		
	}//end for
	
	//cout <<"vQ "<<vQ.x<<" " <<vQ.y<<endl;
	//cout <<"MinAngle: "<< nMinAngle<<endl;
	
	T nAngP = 100;
	T nAngQ = 100;
	
	/// find point R that minimizes angle PRQ 
	while(true)
	{	
		nMinPRQ = 1.7E+308;
		for(i = 0; i < (int)vPts.size(); i++)
		{
			if(i == iP || i == iQ)
				continue;
			
			CVector2<T> vRQ = CVector2<T>::createVector(vPts[i], vQ);
			CVector2<T> vRP = CVector2<T>::createVector(vPts[i], vP);
			
			T nNorm0  = vRQ.mag();
			T nNorm1  = vRP.mag();
			
			T nDot    = vRQ * vRP;
			
			nDot		  /= nNorm0 * nNorm1;
			
			T nAngle = acos(nDot);
			
			nAngle = nAngle * 180/pi;
			
			if(nAngle < nMinPRQ)
			{
				nMinPRQ = nAngle;
				vR = vPts[i];
				iR = i;
			}//end if
							
		}//end for
		
		CVector2<T> vPR = CVector2<T>::createVector(vP,vR);
		CVector2<T> vPQ = CVector2<T>::createVector(vP,vQ);
		CVector2<T> vQR = CVector2<T>::createVector(vQ, vR);
		CVector2<T> vQP = CVector2<T>::createVector(vQ, vP);
		
		nAngP = vPR * vPQ;
		
		T n1 = vPR.mag();
		T n2 = vPQ.mag();
		
		nAngP /= n1 * n2;
		
		nAngP = acos(nAngP);
		
		nAngP = nAngP * 180/pi;
		
		nAngQ = vQR * vQP;
		
		n1 = vQR.mag();
		n2 = vQP.mag();
		
		nAngQ /= n1 * n2;
		
		nAngQ = acos(nAngQ);
		
		nAngQ = nAngQ * 180/pi;
		
		/// the angle at r is obtuse, the center of the circle is (p+q)/2 
		if(nMinPRQ > 90)
		{
			vCenter = vQ + vP;
			vCenter = vCenter * 0.5;
			CVector2<T> vVec1 = vP - vQ;
			vVec1 = vVec1 * 0.5;
			nRad = vVec1.mag();
			//cout<<"Loesung aus 2 Punkten"<<endl;
			break;
		}//end if
		else
		{
			/// angle at r,p,q < 90, the algorithm is finished... 
			if(nAngQ <= 90 && nAngP <= 90)
			{
				/// create middle senkrechte for PR 
				CVector2<T> vMPR = vP + (vPR * 0.5);
				CVector2<T> vNPR = CVector2<T>::createVector(vP,vR);
				vNPR.Normalize();
				CVector2<T> vP2 = vMPR + vNPR;
				//g_NPR = vP2;
				CParamLine P1(vMPR, vP2);
				//g_PR = vMPR;

				/// create middle senkrechte for PQ 
				CVector2<T> vMPQ = vP + (vPQ * 0.5);
				CVector2<T> vNPQ = CVector2<T>::createVector(vP,vQ);
				vNPQ.Normalize();
				CVector2<T> vP1 = vMPQ - vNPQ;
				//g_NPQ = vP1;
				CParamLine P2(vMPQ, vP1);
				//g_PQ = vMPQ;

				//g_P = vP;
				//g_Q = vQ;
				//g_R = vR;
				
				T t1,t2;
				
				if(CParamLine::IntersectParamLines2D(P1, P2, &t1, &t2))
				{
				
				vCenter = P1.Evaluate(t1);
				CVector2<T> vV0 = CVector2<T>::createVector(vCenter, vP);
				nRad = vV0.mag();
				//cout<<"Loesung aus 3 Punkten, Mittelpunkt: "<<vCenter.x<<" "<<vCenter.y<<endl;
				break;
				}
				else
				{
					cout<<"There is some kind of error"<<endl;
					break;
				}//
			}//end if
			/// check whether the angle at p or q is obtuse and replace this point for r 
			else
			{
				if(nAngQ >= 90)
				{
					vQ = vR;
					iQ = iR;
				}
				else if(nAngP >= 90)
				{
					vP = vR;
					iP = iR;
				}
				//cout <<"Winkel an Q: "<<nAngQ<<endl;
				//cout <<"Winkel an P: "<<nAngP<<endl;
				
			}//end else
			
		}//end else
		
	}//end while
		
	//cout <<"Center "<<vCenter.x<<" "<<vCenter.y<<endl;
	//cout <<"vR "<<vR.x<<" " <<vR.y<<endl;
	//cout <<"Min PRQ: "<<nMinPRQ<<endl;
	//return vCenter;

	CCircle<T> circle(vCenter, nRad);
	return circle;
}//end smallestEnclosingCircle


/*!
    \fn CCircle::InitCircle(VECTOR2 vMid, double dRad)
 */
template<class T>
void CCircle<T>::InitCircle(CVector2<T> vMid, T dRad)
{
   m_vCenter = vMid;
   m_dRad = dRad;
}


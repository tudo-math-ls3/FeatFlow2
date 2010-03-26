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

#include <obb2.h>
#include <bspline.h>
#include <approxcurve.h>


//standard constructor
COBB::COBB()
{
	
	for(int i = 0; i < 4; i++)
	{
		m_Vertices[i] = VECTOR2(0,0);
	}//end for
	
	m_pCurve = NULL;
	
}

/////////////////////////////////////////////////////////
//                  Constructor                        //
// Input: Two points that define a box				   //
//        values.                                      //
/////////////////////////////////////////////////////////

COBB::COBB(VECTOR2 vA, VECTOR2 vB, CBSpline *pCurve)
{

	/* bottom left vertex */
	m_Vertices[0].x = vA.x;
	m_Vertices[0].y = vA.y;
	/* bottom right vertex */
	m_Vertices[1].x = vB.x;
	m_Vertices[1].y = vA.y;
	
	/* top right vertex */
	m_Vertices[2].x = vB.x; 
	m_Vertices[2].y = vB.y;
	
	/* top left vertex */
	m_Vertices[3].x = vA.x;
	m_Vertices[3].y = vB.y;
	
	if(pCurve)
		m_pCurve = pCurve;
	else
		pCurve = NULL;

}//end constructor



COBB::COBB(COBB* pBox, CApproxCurve *pCurve)
{
	memcpy(m_Vertices, pBox->m_Vertices, 4*sizeof(VECTOR2));
	m_paCurve = pCurve;
}//end constructor

COBB::COBB(CPArray &cP, vector<int>& hP, CBSpline *pCurve)
{

	m_pCurve = pCurve;

	Real yB =  1.7E+308;
	Real yT = -1.7E+308;
	Real xL =  1.7E+308;
	Real xR = -1.7E+308;

	int num = (int)hP.size();

	for(int i = 0; i < num;i++)
	{
		if(yB > cP[hP[i]].y)
		{
			yB = cP[hP[i]].y;
		}
		if(yT < cP[hP[i]].y)
		{
			yT = cP[hP[i]].y;
		}
		if(xL > cP[hP[i]].x)
		{
			xL = cP[hP[i]].x;
		}
		if(xR < cP[hP[i]].x)
		{
			xR = cP[hP[i]].x;
		}
	}
	
	m_Vertices[0].x = xL;
	m_Vertices[0].y = yB;
	m_Vertices[1].x = xR;
	m_Vertices[1].y = yB;
	m_Vertices[2].x = xR; 
	m_Vertices[2].y = yT;
	m_Vertices[3].x = xL;
	m_Vertices[3].y = yT;

}//end constructor

COBB::COBB(const vector<VECTOR2> &cP, const vector<int> &hP, CBSpline *pCurve)
{

	m_pCurve = pCurve;

	Real yB =  1.7E+308;
	Real yT = -1.7E+308;
	Real xL =  1.7E+308;
	Real xR = -1.7E+308;

	int num = (int)hP.size();

	for(int i = 0; i < num;i++)
	{
		if(yB > cP[hP[i]].y)
			yB = cP[hP[i]].y;
		if(yT < cP[hP[i]].y)
			yT = cP[hP[i]].y;
		if(xL > cP[hP[i]].x)
			xL = cP[hP[i]].x;
		if(xR < cP[hP[i]].x)
			xR = cP[hP[i]].x;
		
	}
	

	m_Vertices[0].x = xL;
	m_Vertices[0].y = yB;
	m_Vertices[1].x = xR;
	m_Vertices[1].y = yB;
	m_Vertices[2].x = xR; 
	m_Vertices[2].y = yT;
	m_Vertices[3].x = xL;
	m_Vertices[3].y = yT;

}//end constructor



ostream& operator<<(ostream& out, COBB *b1) 
{
	for(int i = 0; i < 4; i++)
	{
		out << b1->m_Vertices[i];
	}

	return out;
}

////////////////////////////
//  Simple test:           
//  if p is inside the box
//	Input : Point to be tested
//  Output: true if inside else false
//
////////////////////////////
bool COBB::inside(const VECTOR2& p)
{
	
	
	for(int i = 0; i < 4; i++)
	{
		VECTOR2 vA = VECTOR2::createVector(m_Vertices[i], m_Vertices[(i+1)%4]);
		
		VECTOR2 vB = VECTOR2::CreateNormal(vA);

		vB.Normalize();

		VECTOR2 vC = VECTOR2::createVector(m_Vertices[i], p);

		vC.Normalize();

		Real dDot  = vC * vB;

		if(!((dDot > 0.0) && (dDot < 1.0)))
		{
			
			return false;
		}
	}

	return true;

	//return (((p.x >= m_Vertices[0].x) && (p.x <= m_Vertices[1].x)) && ((p.y >= m_Vertices[0].y) && (p.y <= m_Vertices[3].y)));

}//end inside

//////////////////////////////////////////
// Test if two boxes overlap            
// Input : extreme coordinates of the box
// Output: true if overlap else false
//////////////////////////////////////////
bool COBB::intersection(Real r, Real l, Real t, Real b)
{

	bool xOverlap = false;
	bool yOverlap = false;
	
	// check for xOverlap
	if(r > m_Vertices[1].x)
	{
		if(l < m_Vertices[1].x)
			xOverlap = true;
		else
			return false;
	}
	else
	{
		if(m_Vertices[0].x < r)
			xOverlap = true;
		else
			return false;
	}
	// check for yOverlap
	if(t > m_Vertices[2].y)
	{
		if(b < m_Vertices[2].y)
			yOverlap = true;
		else
			return false;
	}
	else
	{
		if(m_Vertices[1].y < t)
			yOverlap = true;
		else
			return false;
	}
	
	return xOverlap && yOverlap;

}

bool COBB::Overlap(COBB *pBox)
{

	bool bXOverlap = false;
	bool bYOverlap = false;

	if( (m_Vertices[0].x <= pBox->m_Vertices[0].x) && (pBox->m_Vertices[0].x <= m_Vertices[1].x) )
		bXOverlap = true;

	else if( (pBox->m_Vertices[0].x <= m_Vertices[0].x) && (pBox->m_Vertices[1].x >= m_Vertices[0].x) )
		bXOverlap = true;

	if( (m_Vertices[0].y <= pBox->m_Vertices[0].y) && (pBox->m_Vertices[0].y <= m_Vertices[3].y) )
		bYOverlap = true;

	else if( (pBox->m_Vertices[0].y <= m_Vertices[0].y) && (pBox->m_Vertices[3].y >= m_Vertices[0].y) )
		bYOverlap = true;

	return ( bXOverlap && bYOverlap );


}//end overlap

void COBB::getDimension()
{
	Real length = m_Vertices[2].x - m_Vertices[3].x;
	Real height = m_Vertices[3].y - m_Vertices[0].y;
	VECTOR2 vCenter = m_Vertices[2] + m_Vertices[0];
	vCenter = vCenter * 0.5;
	//printf("Length = %f, Height = %f center = (%f,%f)\n", length, height, vCenter.x, vCenter.y);
}//end getDimension

/*!
    \fn COBB::InitBox(VECTOR2 vBL, VECTOR2 vRT)
 */
void COBB::InitBox(VECTOR2 vBL, VECTOR2 vRT)
{
	/* bottom left vertex */
	m_Vertices[0].x = vBL.x;
	m_Vertices[0].y = vBL.y;
	/* bottom right vertex */
	m_Vertices[1].x = vRT.x;
	m_Vertices[1].y = vBL.y;
	
	/* top right vertex */
	m_Vertices[2].x = vRT.x; 
	m_Vertices[2].y = vRT.y;
	
	/* top left vertex */
	m_Vertices[3].x = vBL.x;
	m_Vertices[3].y = vRT.y;
	
	m_pCurve = NULL;


}


/*!
    \fn COBB::COBB(std::vector<VECTOR2> &vPoints)
 */
COBB::COBB(std::vector<VECTOR2> &vPoints)
{

	Real yB =  1.7E+308;
	Real yT = -1.7E+308;
	Real xL =  1.7E+308;
	Real xR = -1.7E+308;

	int num = (int)vPoints.size();

	for(int i = 0; i < num;i++)
	{
		if(yB > vPoints[i].y)
			yB = vPoints[i].y;
		if(yT < vPoints[i].y)
			yT = vPoints[i].y;
		if(xL > vPoints[i].x)
			xL = vPoints[i].x;
		if(xR < vPoints[i].x)
			xR = vPoints[i].x;
		
	}

	m_Vertices[0].x = xL;
	m_Vertices[0].y = yB;
	m_Vertices[1].x = xR;
	m_Vertices[1].y = yB;
	m_Vertices[2].x = xR; 
	m_Vertices[2].y = yT;
	m_Vertices[3].x = xL;
	m_Vertices[3].y = yT;
	
	m_pCurve = NULL;	
	
}//end constructor


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

#include "aabb2.h"
#include "bspline.h"
#include "approxcurve.h"
#include "obb2.h"


CAABB::CAABB(const vector<COBB*> &vBoxes)
{

	Real yT =  -std::numeric_limits<Real>::max();
	Real yB =   std::numeric_limits<Real>::max();
	Real xL =   std::numeric_limits<Real>::max();
	Real xR =  -std::numeric_limits<Real>::max();

	vector<COBB*>::const_iterator bIter;

	for(bIter = vBoxes.begin(); bIter != vBoxes.end(); bIter++)
	{
		COBB *box = *bIter;
		if(box->m_Vertices[0].x < xL)
			xL = box->m_Vertices[0].x;

		if(box->m_Vertices[2].x > xR)
			xR = box->m_Vertices[2].x;

		if(box->m_Vertices[0].y < yB)
			yB = box->m_Vertices[0].y;

		if(box->m_Vertices[2].y > yT)
			yT = box->m_Vertices[2].y;

	}

	m_Vertices[0].x = xL;
	m_Vertices[0].y = yB;
	m_Vertices[1].x = xR;
	m_Vertices[1].y = yT;

	

}//end constructor

void CAABB::InitBox(const vector<COBB*> &vBoxes)
{

	Real yT =  -std::numeric_limits<Real>::max();
	Real yB =   std::numeric_limits<Real>::max();
	Real xL =   std::numeric_limits<Real>::max();
	Real xR =  -std::numeric_limits<Real>::max();

	vector<COBB*>::const_iterator bIter;

	for(bIter = vBoxes.begin(); bIter != vBoxes.end(); bIter++)
	{
		COBB *box = *bIter;
		if(box->m_Vertices[0].x < xL)
			xL = box->m_Vertices[0].x;

		if(box->m_Vertices[2].x > xR)
			xR = box->m_Vertices[2].x;

		if(box->m_Vertices[0].y < yB)
			yB = box->m_Vertices[0].y;

		if(box->m_Vertices[2].y > yT)
			yT = box->m_Vertices[2].y;

	}

	m_Vertices[0].x = xL;
	m_Vertices[0].y = yB;
	m_Vertices[1].x = xR;
	m_Vertices[1].y = yT;


}//end InitBox

CAABB::CAABB(const VECTOR2 &vA, const VECTOR2 &vB)
{
	m_Vertices[0] = vA;
	m_Vertices[1] = vB;
}//end constructor

bool CAABB::GetLongestAxis()
{
	if((m_Vertices[1].x - m_Vertices[0].x) > (m_Vertices[1].y - m_Vertices[0].y))
		return xAxis;
	else
		return yAxis;
}//end getLongestAxis



void CAABB::CalculateBounds(const vector<COBB*> &vBoxes)
{

	Real yT =  -std::numeric_limits<Real>::max();
	Real yB =   std::numeric_limits<Real>::max();
	Real xL =   std::numeric_limits<Real>::max();
	Real xR =  -std::numeric_limits<Real>::max();

	vector<COBB*>::const_iterator bIter;

	for(bIter = vBoxes.begin(); bIter != vBoxes.end(); bIter++)
	{
		COBB *box = *bIter;

		CBasicCurve *curve = box->m_pCurve;

		if(box->m_Vertices[0].x < xL)
			xL = box->m_Vertices[0].x;

		if(box->m_Vertices[2].x > xR)
			xR = box->m_Vertices[2].x;

		if(box->m_Vertices[0].y < yB)
			yB = box->m_Vertices[0].y;

		if(box->m_Vertices[2].y > yT)
			yT = box->m_Vertices[2].y;

	}

	m_Vertices[0].x = xL;
	m_Vertices[0].y = yB;
	m_Vertices[1].x = xR;
	m_Vertices[1].y = yT;

}//end calculateBounds

void CAABB::CalculateBounds(const CAABB &b1, const CAABB &b2)
{
	m_Vertices[0].x = (b1.m_Vertices[0].x < b2.m_Vertices[0].x) ? b1.m_Vertices[0].x : b2.m_Vertices[0].x;

	m_Vertices[1].x = (b1.m_Vertices[1].x > b2.m_Vertices[1].x) ? b1.m_Vertices[1].x : b2.m_Vertices[1].x;

	m_Vertices[0].y = (b1.m_Vertices[0].y < b2.m_Vertices[0].y) ? b1.m_Vertices[0].y : b2.m_Vertices[0].y;

	m_Vertices[1].y = (b1.m_Vertices[1].y > b2.m_Vertices[1].y) ? b1.m_Vertices[1].y : b2.m_Vertices[1].y;

}//end calculateBounds

bool CAABB::Inside(const VECTOR2& p)
{
	return (((p.x >= m_Vertices[0].x) && (p.x <= m_Vertices[1].x)) && ((p.y >= m_Vertices[0].y) && (p.y <= m_Vertices[1].y)));
}//end inside
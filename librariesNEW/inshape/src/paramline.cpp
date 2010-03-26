/***************************************************************************
 *   Copyright (C) 2007 by Raphael M�nster   *
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


#include <paramline.h>

////////////////////////////
/*
*
*	a simple 2d parametric
*	line
*/
////////////////////////////
CParamLine::CParamLine(const VECTOR2 &vP0, const VECTOR2 &vP1)
{
	m_vP0 = vP0;
	m_vP1 = vP1;
	m_vDir = VECTOR2::createVector(vP0, vP1);
}//end constructor

VECTOR2 CParamLine::Evaluate(double t)
{
	VECTOR2 res;
	res = m_vP0 + (m_vDir * t);
	return res;
}//end evalParamLine

int CParamLine::IntersectParamLines2D(CParamLine &p1, CParamLine &p2, double *t1, double *t2)
{

	double nDet = p1.m_vDir.x * p2.m_vDir.y - p1.m_vDir.y * p2.m_vDir.x;
	
	if(nDet <= EPSILON4)
	{
		*t1 = -1;
		*t2 = -1;
		return 0;
	}
	
	*t1 = (p2.m_vDir.x * (p1.m_vP0.y - p2.m_vP0.y) - p2.m_vDir.y * (p1.m_vP0.x - p2.m_vP0.x))/nDet;
	
	*t2 = (p1.m_vDir.x * (p1.m_vP0.y - p2.m_vP0.y) - p1.m_vDir.y * (p1.m_vP0.x - p2.m_vP0.x))/nDet;

	return 1;
		
}//end intersectParamLines2D
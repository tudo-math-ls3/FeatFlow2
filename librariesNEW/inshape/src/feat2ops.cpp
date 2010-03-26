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


#include <feat2ops.h>
#include <approxcurve.h>
#include <obb2.h>


CFeat2Ops::CFeat2Ops(void)
{

}//end constructor

CFeat2Ops::~CFeat2Ops(void)
{

}//end deconstructor

Real CFeat2Ops::IsFictitiousBoundarySimple(CApproxCurve *pCurve, const VECTOR2 &vQuery)
{

	Real result;

	/* check whether the point is inside the interface */
	if(pCurve->m_bBox->Inside(vQuery))
	{
		
		if(pCurve->IsInElement(vQuery))
		{
			result = (Real)7.0;
			return result;
		}//end if
		else
		{
			result = (Real)10.0;
			return result;
		}//end else
		
	}//end if

	result = (Real)10.0;
	return result;

}//end IsFictitiousBoundary

bool CFeat2Ops::IsFictitiousBoundary(CApproxCurve *pCurve, const VECTOR2 &vQuery)
{

	/* check whether the point is inside the interface */
	if(pCurve->m_bBox->Inside(vQuery))
	{
		
		if(pCurve->IsInElement(vQuery))
		{
			return true;
		}//end if
		else
		{
			return false;
		}//end else
		
	}//end if

	return false;

}//end IsFictitiousBoundarySimple
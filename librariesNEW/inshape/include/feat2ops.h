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

#ifndef _FEAT2OPS_H_
#define _FEAT2OPS_H_

#ifdef WIN32
#pragma once
#endif

//===================================================	
//			         INCLUDES
//===================================================

#include "distops.h"

//===================================================	
//			      STRUCTS & TYPES
//===================================================

typedef struct
{
	Real rSolValue;
}t_FBCInfo;

class CFeat2Ops :
	public DistOps
{
public:

	//standard constructor
	CFeat2Ops(void);

	//deconstructor
	~CFeat2Ops(void);

	//test whether vQuery is inside the interface pointed to by pCurve
	Real IsFictitiousBoundarySimple(CApproxCurve *pCurve, const VECTOR2 &vQuery);

	//test whether vQuery is inside the interface pointed to by pCurve, returns true if inside
	bool IsFictitiousBoundary(CApproxCurve *pCurve, const VECTOR2 &vQuery);

};

#endif
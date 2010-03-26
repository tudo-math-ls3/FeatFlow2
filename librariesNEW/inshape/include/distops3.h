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

#ifdef WIN32
#pragma once
#endif

#ifndef _CDISTOPS3_H_
#define _CDISTOPS3_H_


#include "distops.h"
#include "grid3d.h"
#include "model3d.h"
#include "cbvhnode.h"

typedef struct
{
	VECTOR3 vVec;
	Real    rDistance;
}Res_t;

class CDistOps3 :
	public DistOps
{
public:
	CDistOps3(void);
public:
	~CDistOps3(void);

	Real BruteForceDistance(const CModel3D &model, const VECTOR3 &vQuery) const;
	



public:
	bool BruteForceInnerPoints(CModel3D & model, VECTOR3& vQuery);
public:
	Real SimpleLUB(AABBTree3f &tree, const VECTOR3 &vQuery);
	Res_t CoherencyLUB(const AABBTree3f &tree, const VECTOR3 &vQuery, Real rLUB);
	Real DistTriPoint(vector<VECTOR3*>& vVerts, const VECTOR3 &vQuery, VECTOR3& vNearest);
};

#endif

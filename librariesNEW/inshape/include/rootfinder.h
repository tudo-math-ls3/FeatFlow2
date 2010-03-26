/***************************************************************************
 *   Copyright (C) 2006 by Raphael M�nster   *
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

#if !defined(_CROOTFINDINGALG_H)
#define _CROOTFINDINGALG_H

#include <math.h>
#include "vector2.h"

class CBSpline;
class CONXNurbs;

class CRootFinder
{
public:
	double bisect(CBSpline* pCurve, const VECTOR2& pvPoint, double xl , double xr, double e);
	double bisectWrap(CBSpline* pCurve, const VECTOR2& pvPoint, double xl , double xr, double e);
	double bisectX(CBSpline* pCurve, double xl , double xr, double e);
	double bisectY(CBSpline* pCurve, double xl , double xr, double e);
	bool Newton( CBSpline* pCurve, const VECTOR2& pvPoint, double& dX0, double xl, double xr, double de, int imaxIterations);
	double evalF(CBSpline* pCurve, const VECTOR2& pvPoint, double dX, double* y1);
	double NewtonBisect(CBSpline* pCurve, const VECTOR2& pvPoint, double e, double xl, double xr);
	double IntersectWrap(CBSpline* pCurve, double y, double e, double xl, double xr);
	double NewtonIntersec(CBSpline* pCurve, double y, double de, double xl, double xr, int numIter);
	bool ONX_bisect(CONXNurbs* pCurve, const VECTOR2& pvPoint, double xl , double xr, double e, double &u);
	bool ONX_Newton(CONXNurbs* pCurve, const VECTOR2& pvPoint, double& dX0, double de, int imaxIterations);
	double bisectIntersec(CBSpline* pCurve, double y, double e, double xl, double xr, VECTOR2 &vPoint);

	double x;
	int n0;
	
};


#endif


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

#ifdef WIN32
#pragma once
#endif

#if !defined(_CNURBS_H)
#define _CNURBS_H
//#define _USE_QT

#include "bspline.h"


class CNurbs : public CBSpline
{
public:

	CP3Array m_vHomCoords;
	CNurbs(int deg, int numP, Real *kV, Real *w, VECTOR2 *cP, int nWrapped);
	CNurbs(){};
	~CNurbs();
	VECTOR2 deBoor(Real u_v,bool draw = false);
	VECTOR2 CoxDeBoor(Real u);
	void setWeight(Real val, int pos);
	VECTOR3* curveDerivs3(Real u, int n);
	VECTOR2* curveDerivs(Real u, int n);
	VECTOR2 CoxDeBoorClamped(Real u);
	void insertKnot(Real u, int t);
	void scale(Real dScale);
	void translateCurve(VECTOR2 vTrans);
	void rotateCurve(Real dAngle);
	void rotateBack(Real dAngle);
	void rotate(Real dAngle);
	void updateControlPoints();
	void toFile(const char *fileName);
    void GenSamples(std::vector<VECTOR2> &samples, int nS, Real start, Real end);
    

	

	

};



class CONXNurbs : public CBSpline
{

public:

	CONXNurbs(const char *sFileName);
	CONXNurbs(const ON_NurbsCurve& curve);
	VECTOR2 CoxDeBoor(Real u);
	void genBoxes(int nIter);
	void findMonotoneSegments(std::vector<VECTOR2>& vPoints, std::vector<Real>& vRoots, int nIters);
	VECTOR2 curveDerivs(Real u);
	VECTOR2 Ev2Deriv(Real u);
	void createHierarchy(int i = 0);
	double* KnotVector(){return m_Curve.m_knot;};
	void scale(Real dScale);
	void CreateSimpleBVH();
	void GenSamples(std::vector<VECTOR2> &samples, int nS, Real start, Real end);
    void BuiltBvh();
	ON_NurbsCurve m_Curve;

};
//end class CONXNurbs

#endif  //_CNURBS_H

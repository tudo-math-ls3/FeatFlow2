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

#ifndef _CIRCLE_H_
#define _CIRCLE_H_

#include <vector>
#include "vector2.h"
#include "matrix2x2.h"
#include "paramline.h"

using namespace std;

template<class T>
class CCircle
{

public:
	/* constructor */
	CCircle(const CVector2<T>& vCenter, T dRad);
	CCircle(){};

	/* inline member functions */
	inline T GetRadius(){return m_dRad;};
	inline CVector2<T> GetCenter(){return m_vCenter;};

	inline T MinDistance(const CVector2<T> &vQuery){T nDist = (vQuery - m_vCenter).mag(); return nDist - m_dRad;};

	inline T MinDistanceSqr(const CVector2<T> &vQuery){T nDist = (vQuery - m_vCenter).norm2(); return nDist - m_dRad;};

	inline T MaxDistance(const CVector2<T> &vQuery){return (T)0.0;};

	inline T MaxDistanceSqr(const CVector2<T> &vQuery){return (T)0.0;};

	//public member functions
	CCircle SmallestEnclosingCircle(vector< CVector2<T> > &vPts);

    /*!
        \fn CCircle::SetCenter(VECTOR2 vMid)
     */
    inline void SetCenter(const CVector2<T> &vMid)
    {
        m_vCenter = vMid;
    }

	inline void Translate(const CVector2<T> &vVec) {m_vCenter+=vVec;};

	inline void Rotate(T nAngle)
	{
		CMatrix2x2<T> matRotation;
		matRotation.InitRotationMatrix(nAngle);
		m_vCenter = matRotation * m_vCenter;
	};

    void InitCircle(CVector2<T> vMid, T dRad);


	/* private member variables */
	CVector2<T> m_vCenter;
	T    m_dRad;

};

//include the function definitions
#include "circle.cpp"

typedef CCircle<double> CCircled;
typedef CCircle<float>  CCirclef;

typedef CCircle<Real>   CCircler;

#endif
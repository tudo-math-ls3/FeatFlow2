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

#ifdef WIN32
#pragma once
#endif

//===================================================
//					DEFINES
//===================================================

#if !defined _QUATERNION_
#define _QUATERNION_

#define PI 3.14159

//===================================================
//					INCLUDES
//===================================================
#include <iostream>
#include <math.h>
#include "matrix4x4.h"

//===================================================
//	    		FORWARD DEFINITIONS
//===================================================


template<class T>
class CQuaternion
{

public:

	/* default constructor */
	CQuaternion() : x(0), y(0), z(0), w(1) {};

	CQuaternion(T a, T b, T c, T d) : x(a), y(b), z(c), w(d) {};

	/* copy constructor */
	CQuaternion( const CQuaternion &q ) : x(q.x), y(q.y), z(q.z), w(q.w)
	{

	};



	inline void SetValues(T X, T Y, T Z, T W)
	{
		x = X;
		y = Y;
		z = Z;
		w = W;
	}//end SetValues

	inline CQuaternion operator - () const
	{
		return CQuaternion(-x,-y,-z,w);
	};

	inline T norm2()
	{
		return  (x * x) + (y * y) + (z * z) + (w * w);
	}//end  operator

	inline double Mag()
	{
		return sqrt(norm2());
	}//end  Mag

	inline void Normalize()
	{
		double InvMag = 1.0/Mag();
		x*= InvMag;
		y*= InvMag;
		z*= InvMag;
		w*= InvMag;
	}//end Normalize

	inline const CQuaternion& GetNormalized()
	{
		double InvMag = 1.0/Mag();
		x*= InvMag;
		y*= InvMag;
		z*= InvMag;
		w*= InvMag;
		return *this;
	}//end GetNormalized

	inline CQuaternion operator *(const CQuaternion &q) const
	{
		CVector3<T> v1(x,y,z);
		CVector3<T> v2(q.x, q.y, q.z);

		
		T nS = w*q.w - (v1 * v2);

		CVector3<T> q_vec = v2 * w + v1 * q.w + (CVector3<T>::Cross(v1,v2));

		return CQuaternion(q_vec.x, q_vec.y, q_vec.z, nS).Normalize();

	}//end operator *

	friend CQuaternion GetQuaternion(const CVector3<T> &v1, const CVector3<T> &v2);

	void CreateMatrix(CMatrix4x4<T> &pMatrix);

	void AxisAngleToQuat(T X, T Y, T Z, T W);

	union
	{
		T m_Data[4];

		struct
		{
			T x;
			T y;
			T z;
			T w;
		};
	};





};//end class Quaternion

template <class T>
void CQuaternion<T>::CreateMatrix(CMatrix4x4<T> &pMatrix)
{

	T xx = x * x;
	T xy = x * y;
	T xz = x * z;
	T xw = x * w;
	T yy = y * y;
	T yz = y * z;
	T yw = y * w;
	T zz = z * z;
	T zw = z * w;

	pMatrix.m_Entries[0] = 1 - 2 * (yy + zz);
	pMatrix.m_Entries[1] = 2 * ( xy + zw );
	pMatrix.m_Entries[2] = 2 * ( xz - yw );
	pMatrix.m_Entries[3] = 0;

	pMatrix.m_Entries[4] = 2 * (xy - zw);
	pMatrix.m_Entries[5] = 1 - 2 * (xx + zz);
	pMatrix.m_Entries[6] = 2 * (yz + xw);
	pMatrix.m_Entries[7] = 0;

	pMatrix.m_Entries[8]  = 2 * (xz + yw);
	pMatrix.m_Entries[9]  = 2 * (yz - xw);
	pMatrix.m_Entries[10] = 1 - 2 * (xx + yy);
	pMatrix.m_Entries[11] = 0;

	pMatrix.m_Entries[12] = 0;
	pMatrix.m_Entries[13] = 0;
	pMatrix.m_Entries[14] = 0;
	pMatrix.m_Entries[15] = 1;

}//end CreateMatrix

template <class T>
void CQuaternion<T>::AxisAngleToQuat(T X, T Y, T Z, T W)
{

	T p = static_cast<T>(180);

	T angle = static_cast<T>(( W / p) * PI) ;

	T result = static_cast<T>(sin(angle / static_cast<T>(2) ));

	w = static_cast<T>(cos(angle / static_cast<T>(2) ));

	x = static_cast<T> (X * result);
	y = static_cast<T> (Y * result);
	z = static_cast<T> (Z * result);

}//end AxisAngleToQuat

typedef CQuaternion<float> CQuaternionf;
typedef CQuaternion<double> CQuaterniond;

template<class T>
CQuaternion<T> GetQuaternion(const CVector3<T> &v1, const CVector3<T> &v2)
{

	CQuaternion<T> q;

	q.x = v1.y * v2.z - v1.z * v2.y;

	q.y = v1.z * v2.x - v1.x * v2.z;

	q.z = v1.x * v2.y - v1.y * v2.x;

	q.z = (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);

	return q;

}//end constructor


#endif
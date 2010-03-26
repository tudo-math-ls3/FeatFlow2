/***************************************************************************
 *   Copyright (C) 2006 by Raphael Münster   *
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




#if !defined(_CVector3_H)
#define _CVector3_H

//===================================================
//					INCLUDES
//===================================================

#include <iostream>
#include "mathglobals.h"

using namespace std;


template<class T>
class CVector3 {


public:
	/* constructor */
	CVector3(T a, T b, T c): x(a), y(b), z(c) {}

	/* copy constructor */
	CVector3(const CVector3 &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
	}

    /* default constructor */
	CVector3():x(0), y(0), z(0){}
	~CVector3(){};
	
	inline const CVector3& operator=(const CVector3& v)
	{
		
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}//end  operator
	
	inline CVector3 operator+(CVector3 v) const
	{
		return CVector3(x + v.x, y + v.y, z + v.z);
	}//end  operator

	inline CVector3 operator-(CVector3 v) const
	{
		return CVector3(x - v.x, y - v.y, z - v.z);
	}//end  operator

	inline CVector3 operator*(T num) const
	{
		// Return scaled vector
		return CVector3(x * num, y * num, z * num);
	}//end  operator

	inline T operator * (const CVector3 &rhs) const
	{
		return  x * rhs.x +
				y * rhs.y +
				z * rhs.z;
				
	}//end  operator


	inline T norm2()
	{
		return  (x * x) + (y * y) + (z * z);
	}//end  operator

	inline double mag()
	{
		return sqrt(norm2());
	}//end  operator

	inline void Normalize()
	{
		double dInvMag = 1.0/mag();
		x *= dInvMag;
		y *= dInvMag;
		z *= dInvMag;
	}//end Normalize

	inline static CVector3 createVector(const CVector3 &a, const CVector3 &b)
	{
		CVector3 res = b - a;
		return res;
	}//end  operator


	inline const CVector3& operator /= (const T &rhs)
	{
		x /= rhs;
		y /= rhs;
		z /= rhs;
		return *this;
	}//end  operator


	inline const CVector3& operator += (const CVector3 &rhs)
	{

		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return *this;
	}//end  operator

	
	inline const CVector3& operator -= (const CVector3 &rhs)
	{

		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		
		return *this;
	}//end  operator

	
	inline const CVector3& operator *= (const T &d)
	{
		x *= d;
		y *= d;
		z *= d;
		
		return *this;
	}//end  operator

	inline T dot(const CVector3 &a, const CVector3 &b)
	{
		return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
	}//end  operator

		
	inline static CVector3 Cross(CVector3 vVector1, CVector3 vVector2)
	{
		CVector3 vCross;

		vCross.x = ((vVector1.y * vVector2.z) - (vVector1.z * vVector2.y));

		vCross.y = ((vVector1.z * vVector2.x) - (vVector1.x * vVector2.z));

		vCross.z = ((vVector1.x * vVector2.y) - (vVector1.y * vVector2.x));

		return vCross;
	}

	friend ostream& operator<<(ostream& out, const CVector3& v1); 
	
	/* union to allow different access methods */
	union
	{
		T m_dCoords[3];

		struct
		{
			T x;
			T y;
		    T z;
		};
	};
	
		
};


/* typedefs to create float and double vectors */
typedef CVector3<double> CVector3d;
typedef CVector3<float>  CVector3f;

typedef CVector3<Real> VECTOR3;


#endif  //_CVector3_H

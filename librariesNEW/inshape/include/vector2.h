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


#if !defined(_CVector2F_H)
#define _CVector2F_H

/*/////////////////////////////////////////////
//				INCLUDES
*//////////////////////////////////////////////


#include <math.h>
#include <iostream>
#include <xmmintrin.h>
#include "mathglobals.h"


#define R_PI 3.1415926535

using namespace std;

/*//////////////////////////////////////////
*
*	CLASS DEFINITION : a templated vector
*	mainly for use with float and double
*
*///////////////////////////////////////////

template<class T>
class CVector2 {

public:

	/* constructors */
	CVector2(T a, T b) : x(a), y(b) {}

	/* default constructor */
	CVector2(): x(0), y(0){}

	/* copy constructor */
	CVector2(const CVector2 &v)
	{
		x = v.x;
		y = v.y;
	}

	inline const CVector2& operator=(const CVector2& v)
	{
		CVector2 tmp(v);
		x = tmp.x;
		y = tmp.y;
		return *this;
	}//end  operator
	
	inline CVector2 operator+(const CVector2& v) const
	{
		return CVector2(x + v.x, y + v.y);
	}//end  operator

	inline CVector2 operator-(const CVector2& v) const
	{
		return CVector2(x - v.x, y - v.y);
	}//end  operator

	inline CVector2 operator*(T num) const
	{
		// Return scaled vector
		return CVector2(x * num, y * num);
	}//end  operator

	inline T operator * (const CVector2& rhs) const
	{
		return  x * rhs.x +
				y * rhs.y;
				
	}//end  operator


	inline T norm2()
	{
		return  (x * x) + (y * y);
	}//end  operator

	inline double mag()
	{
		return sqrt(norm2());
	}//end  operator

	inline static CVector2 createVector(const CVector2 &a, const CVector2 &b)
	{
		CVector2 res = b - a;
		return res;
	}//end  operator


	inline const CVector2& operator /= (const T &rhs)
	{
		x /= rhs;
		y /= rhs;
		return *this;
	}//end  operator


	inline const CVector2& operator += (const CVector2 &rhs)
	{

		x += rhs.x;
		y += rhs.y;
		return *this;
	}//end  operator

	
	inline const CVector2& operator -= (const CVector2 &rhs)
	{

		x -= rhs.x;
		y -= rhs.y;
		
		return *this;
	}//end  operator

	
	inline const CVector2& operator *= (const T &d)
	{
		x *= d;
		y *= d;
		
		return *this;
	}//end  operator

	inline static CVector2 CreateNormal(const CVector2 &vA)
	{
		return CVector2(-vA.y, vA.x);
	}//end  operator

	inline void Normalize()
	{
		double dInvMag = 1.0/mag();
		x *= dInvMag;
		y *= dInvMag;
	}//end Normalize

	inline CVector2 GetNormalizedVector()
	{
		CVector2 vec;
		double dInvMag = 1.0/mag();
		vec.x = x * dInvMag;
		vec.y = y * dInvMag;
		return vec;
	}//end Normalize


	friend ostream& operator<<(ostream& out, const CVector2& v1); 

	inline T dot(const CVector2 &a, const CVector2 &b)
	{
		return (a.x * b.x) + (a.y * b.y);
	}//end  operator

	inline double CalcAngle(CVector2 v1)
	{
		
		CVector2 v2 = *this;
		
		double nAng = v1 * v2;
		
		double n1 = v1.mag();
		double n2 = v2.mag();
		
		nAng /= n1 * n2;
		
		nAng = acos(nAng);
		
		nAng = nAng * 180.0/double(R_PI);
		
		return nAng;
		
	}//end CalcAngle
	
	~CVector2(){};

	/* union to allow different access methods */
	union
	{
		T m_dCoords[2];

		struct
		{
			T x;
			T y;
		};
	};

};

typedef CVector2<double> CVector2d;
typedef CVector2<float>  CVector2f;

typedef CVector2<Real> VECTOR2;

#endif  //_CVector2F_H



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

#if !defined(_CVector4_H_)
#define _CVector4_H_

//===================================================
//					DEFINITIONS
//===================================================
#ifdef WIN32
#pragma once
#endif



//===================================================
//					INCLUDES
//===================================================
#include <iostream>


using namespace std;

//===================================================
//			CLASS DESCRIPTION:
//			    a templated 3d homogeneous vector  
//				class, mainly for use with float
//				and double data type
//===================================================


template<class T>
class CVector4 {


public:
	/* constructor */
	CVector4(T a, T b, T c, T d): x(a), y(b), z(c), w(d) {}

	CVector4(T a, T b, T c): x(a), y(b), z(c), w(1) {}

	/* copy constructor */
	CVector4(const CVector4 &v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		w = 1;
	}

    /* default constructor */
	CVector4():x(0), y(0), z(0), w(1){}
	~CVector4(){};

//===================================================
//  			Assignment operator		
//===================================================
	
	inline const CVector4& operator=(const CVector4& v)
	{
		
		x = v.x;
		y = v.y;
		z = v.z;
		
		return *this;
	}//end  operator
	
	inline CVector4 operator+(const CVector4 &v) const
	{
		return CVector4(x + v.x, y + v.y, z + v.z, 1);
	}//end  operator

	inline CVector4 operator-(const CVector4 &v) const
	{
		return CVector4(x - v.x, y - v.y, z - v.z, 1);
	}//end  operator

	inline CVector4 operator - () const
	{
		return CVector4(-x,-y,-z);
	}//end operator

	inline CVector4 operator*(T num) const
	{
		// Return scaled vector
		return CVector4(x * num, y * num, z * num);
	}//end  operator

	inline T operator * (const CVector4 &rhs) const
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

	inline static CVector4 createVector(const CVector4 &a, const CVector4 &b) 
	{
		CVector4 res = b - a;
		return res;
	}//end  operator


	inline const CVector4& operator /= (const T &rhs)
	{
		x /= rhs;
		y /= rhs;
		z /= rhs;
		return *this;
	}//end  operator


	inline const CVector4& operator += (const CVector4 &rhs)
	{

		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return *this;
	}//end  operator

	
	inline const CVector4& operator -= (const CVector4 &rhs)
	{

		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		
		return *this;
	}//end  operator

	
	inline const CVector4& operator *= (const T &d)
	{
		x *= d;
		y *= d;
		z *= d;
		
		return *this;
	}//end  operator

	inline T dot(const CVector4 &a, const CVector4 &b)
	{
		return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
	}//end  operator

		
	inline static CVector4 Cross(CVector4 vVector1, CVector4 vVector2)
	{
		CVector4 vCross;

		vCross.x = ((vVector1.y * vVector2.z) - (vVector1.z * vVector2.y));

		vCross.y = ((vVector1.z * vVector2.x) - (vVector1.x * vVector2.z));

		vCross.z = ((vVector1.x * vVector2.y) - (vVector1.y * vVector2.x));

		return vCross;
	}//end Cross

	inline void Normalize()
	{
		double dInvMag = 1.0/mag();
		x *= dInvMag;
		y *= dInvMag;
		z *= dInvMag;
		w  = 1;
	}//end Normalize

	friend ostream& operator<<(ostream& out, const CVector4& v1); 
	
	/* union to allow different access methods */
	union
	{
		T m_dCoords[4];

		struct
		{
			T x;
			T y;
		    T z;
			T w;
		};
	};
	
		
};


/* typedefs to create float and double vectors */
typedef CVector4<double> CVector4d;
typedef CVector4<float>  CVector4f;


#endif  //_CVector4_H_

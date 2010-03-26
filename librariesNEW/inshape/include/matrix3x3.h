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



#if !defined(_CMATRIX3X3_H)
#define _CMATRIX3X3_H

//===================================================
//					INCLUDES
//===================================================
#include <iostream>
#include "vector3.h"

template<class T>
class CMatrix3x3 
{
public:
	CMatrix3x3();
	CMatrix3x3( T m00,  T m01,  T m02,  T m10,  T m11,  T m12,
					    T m20,  T m21,  T m22);
	CMatrix3x3( T entries[]);
	~CMatrix3x3(){};

//===================================================
//  		 Matrix-Vector-Product	   
//===================================================
	CVector3<T> operator*(const CVector3<T> &rhs);

//===================================================
//  		 Matrix-Matrix-Product	   
//===================================================

	CMatrix3x3<T> operator*(const CMatrix3x3<T> &rhs) const;
	


	union
	{
		T m_dEntries[9];

		struct
		{
			T m_d00;
			T m_d01;
			T m_d02;
			T m_d10;
			T m_d11;
			T m_d12;
			T m_d20;
			T m_d21;
			T m_d22;
		};
	};

//===================================================
//  		  Set to the identity matrix
//===================================================
	inline void SetIdentity()
	{
		memset(m_dEntries, 0, 9*sizeof(T));
		m_d00 = 1;
		m_d11 = 1;
		m_d22 = 1;
	}//end SetIdentity

//===================================================
//  		  Set to the zero matrix
//===================================================
	inline void SetZero()
	{
		memset(m_dEntries, 0, 9*sizeof(T));
	}//end SetIdentity

//===================================================
//  		 Transpose the matrix
//===================================================
	inline CMatrix3x3 GetTransposedMatrix()
	{
		return CMatrix3x3(m_d00, m_d10, m_d20,
    				      m_d01, m_d11, m_d21,
	    				  m_d02, m_d12, m_d22);

	}//end GetTransposedMatrix

//===================================================
//  		 Transpose the matrix
//===================================================
	inline void TransposeMatrix()
	{
		CMatrix3x3 mat(m_d00, m_d10, m_d20,
    	               m_d01, m_d11, m_d21,
	    	           m_d02, m_d12, m_d22);

		memcpy(m_dEntries, mat.m_dEntries, 9*sizeof(T));

	}//end TransposeMatrix

//===================================================
//  			   Determinate		
//===================================================

	inline T Determinate()
	{
		T det = m_d00 * (m_d11*m_d22 - m_d21*m_d12) - m_d01 * (m_d10*m_d22 - m_d20*m_d12) + m_d02 * (m_d10*m_d21 - m_d20*m_d11);
		return det;
	}//end Determinate

//===================================================
//  			Assignment operator		
//===================================================
	inline void operator=(const CMatrix3x3 &rhs)
	{
		memcpy(m_dEntries, rhs.m_dEntries, 9 * sizeof(T));
	}//end operator=

//===================================================
//  			output operator		
//===================================================
	friend std::ostream &operator << (std::ostream &out, const CMatrix3x3 &rhs);

};


template<class T>
CMatrix3x3<T>::CMatrix3x3()
{
	
	SetZero();
	
}//end constructor

template<class T>
CMatrix3x3<T>::CMatrix3x3(T m00, T m01, T m02, T m10, T m11, T m12,
					   T m20, T m21, T m22)
: m_d00(m00), m_d01(m01), m_d02(m02), m_d10(m10), m_d11(m11), m_d12(m12), m_d20(m20), m_d21(m21), m_d22(m22) 
{
	

}//end constructor

template<class T>
CMatrix3x3<T>::CMatrix3x3(T entries[])
{
	
	memcpy(m_dEntries, entries, 9 * sizeof(T));
	
}//end constructor

template<class T>
CVector3<T> CMatrix3x3<T>::operator *(const CVector3<T> &rhs)
{
	CVector3<T> res;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			res.m_dCoords[i] += m_dEntries[i*3+j] * rhs.m_dCoords[j];

	return res;
}

//tune in time...
template<class T>
CMatrix3x3<T> CMatrix3x3<T>::operator *(const CMatrix3x3<T> &rhs) const
{
	CMatrix3x3<T> tmp;
	
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j <3; j++)
		{
			for(int k = 0; k < 3; k++)
			{
				tmp.m_dEntries[i*3+j] += m_dEntries[i*3+k] * rhs.m_dEntries[k*3+j];
				
			}//end for k
		}//end for j
	}//end for i

	return tmp;
}//end operator*





typedef CMatrix3x3<double> CMatrix3x3d;
typedef CMatrix3x3<float>  CMatrix3x3f;

typedef CMatrix3x3<Real> MATRIX3X3;

#define IJ3(i,j) (i*3 + j)




#endif
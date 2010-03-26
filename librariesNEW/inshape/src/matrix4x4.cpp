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


#include "matrix4x4.h"

std::ostream &operator << (std::ostream &out, const CMatrix4x4f &rhs)
{
	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			out<< rhs.m_Entries[i*4+j]<<" ";
		}//end for
		out<<endl;
	}//end for
	return out;
}//end operator <<

//bool XMatrix::Inverse(XMatrix& result) const
//{
//	XReal determinant = Determinant();
//	XReal matrixSub3x3[9];
//
//	int i, j, sign;
//
//	// if the determinant is 0 (or as near as makes no difference), we cannot 
//	// calculate the inverse
//	if ( fabs(determinant) < 0.0005 )
//		return false;
//
//	for ( i = 0; i < 4; i++ )
//	{
//		for ( j = 0; j < 4; j++ )
//		{
//			sign = 1 - ( (i +j) % 2 ) * 2;
//
//			SubMatrix3x3( i, j, matrixSub3x3 );
//
//
//			XReal det =   matrixSub3x3[0] * ( matrixSub3x3[4]*matrixSub3x3[8] - matrixSub3x3[7]*matrixSub3x3[5] )
//						- matrixSub3x3[1] * ( matrixSub3x3[3]*matrixSub3x3[8] - matrixSub3x3[6]*matrixSub3x3[5] )
//						+ matrixSub3x3[2] * ( matrixSub3x3[3]*matrixSub3x3[7] - matrixSub3x3[6]*matrixSub3x3[4] );
//
//			result.data[i+j*4] = ( det * sign ) / determinant;
//		}
//	}
//	return true;
//}
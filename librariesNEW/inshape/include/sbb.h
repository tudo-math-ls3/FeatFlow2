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

#include <vector>
#include <limits>
#include "vector2.h"

using namespace std;

class CBSpline;

class CSBB 
{
public:
	
	/* constructors */
	CSBB(VECTOR2 vA, VECTOR2 vB, Real fA, Real fB, CBSpline *pCurve);
	CSBB(CSBB *pBox);
	CSBB(std::vector<VECTOR2> vPoints, Real start, Real end, CBSpline *pCurve);
	CSBB();
	
	/* member functions */
	friend ostream& operator<<(ostream& out, CSBB *b1);

	inline void translateBox(VECTOR2 vTrans)
	{
		for(int i = 0; i < 4; i++)
		{
			m_Vertices[i] = m_Vertices[i] + vTrans;
		}//end for
	};

	/* inline member functions */
	inline VECTOR2 getCenter()
	{
		VECTOR2 vCenter;
		vCenter = (m_Vertices[0] + m_Vertices[1] + m_Vertices[2] + m_Vertices[3]) * 0.25;
		return vCenter;
	};


	/* member variables */
	int m_ID;
	/* parameter values */
	Real m_fP1;
	Real m_fP2;

	/* member variables */
	VECTOR2 m_Vertices[4];
	
	/* 0 = bottom left, 1 = bottom right, 2 = top right, 3 = top left */
	/* stores the two points that are on the curve */
	int m_CurvePoints[2];

	CBSpline *m_pCurve;

};
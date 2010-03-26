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
#include "dynamicarray.h"

using namespace std;

class CBSpline;
class CApproxCurve;


/*
*
*	class for an oriented bounding box
*	vertex numbering:
*	
*	3---2
*	|	|
*	0---1
*
*
*/
class COBB
{
public:
	
	/* constructors */
	COBB(VECTOR2 vA, VECTOR2 vB, CBSpline *pCurve = NULL);
	COBB(CDynamicArray<VECTOR2> &cP, vector<int>& hP, CBSpline *pCurve);
	COBB(COBB* pBox, CApproxCurve *pCurve);
	COBB(const vector<VECTOR2> &cP, const vector<int> &hP, CBSpline *pCurve);
	COBB();

	/* member functions */
	bool    inside(const VECTOR2& p);
	bool    intersection(Real r, Real l, Real t, Real b);
	void    getDimension();
	bool    Overlap(COBB *pBox);

	/* inline member functions */
	inline VECTOR2 getCenter()
	{
		VECTOR2 vCenter;
		vCenter = (m_Vertices[0] + m_Vertices[1] + m_Vertices[2] + m_Vertices[3]) * 0.25;
		return vCenter;
	};
		
	inline void translateBox(VECTOR2 vTrans)
	{
		for(int i = 0; i < 4; i++)
		{
			m_Vertices[i] = m_Vertices[i] + vTrans;
		}//end for
	};

	//test if the point is inside the box
	inline bool Inside(const VECTOR2& p)
	{
		return (((p.x >= m_Vertices[0].x) && (p.x <= m_Vertices[1].x)) && ((p.y >= m_Vertices[0].y) && (p.y <= m_Vertices[3].y)));
	}//end inside

	inline Real  getArea(){return ((m_Vertices[1].x - m_Vertices[0].x)*(m_Vertices[2].y - m_Vertices[1].y));};

	/* member variables */
	//make this one an integer ID

	union
	{
		CBSpline *m_pCurve;
		CApproxCurve *m_paCurve;
	};

	friend ostream& operator<<(ostream& out, COBB *b1); 
    void InitBox(VECTOR2 vBL, VECTOR2 vRT);
    COBB(std::vector<VECTOR2> &vPoints);

	/* member variables */
	VECTOR2 m_Vertices[4];

};
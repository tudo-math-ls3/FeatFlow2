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

class COBB;

class CAABB
{
public:
	
	/* constructors */
	CAABB(const VECTOR2 &vA, const VECTOR2 &vB);
	CAABB(const vector<COBB*> &vBoxes);
	CAABB(){};
	
	/* member functions */
	void InitBox(const vector<COBB*> &vBoxes);
	bool GetLongestAxis();
	void CalculateBounds(const vector<COBB*> &vBoxes);
	void CalculateBounds(const CAABB &b1, const CAABB &b2);
	
	bool Inside(const VECTOR2& p);
	
	/* inline member functions */	
	inline VECTOR2 GetCenter()
	{
		VECTOR2 vCenter;
		vCenter = (m_Vertices[0] + m_Vertices[1] +
				   VECTOR2(m_Vertices[0].x, m_Vertices[1].y) + VECTOR2(m_Vertices[1].x, m_Vertices[0].y)) * 0.25;
		return vCenter;
	};

	inline Real Area(){return (m_Vertices[1].x - m_Vertices[0].x) * (m_Vertices[1].y - m_Vertices[0].y);};

	inline void TranslateBox(VECTOR2 vTrans)
	{
		m_Vertices[0] = m_Vertices[0] + vTrans;
		m_Vertices[1] = m_Vertices[1] + vTrans;
	};
	
	/* member variables */
	VECTOR2 m_Vertices[4];

	enum
	{
		xAxis,
		yAxis
	};
		
};
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
#include <list>
#include <set>
#include "geomops.h"
#include "circle.h"
#include "paramline.h"




using namespace std;

typedef struct
{
	int e1;
	int e2;
	double length;
}tEdge;

typedef struct
{
	VECTOR2 vVertex;
	bool      inside;
	int       index;
	std::vector<VECTOR2> m_FPoints;
}tVertex;


// typedef struct
// {
// 	cVector2f vFPoints[3];
// }tFPoints;



typedef struct
{
	VECTOR2 p;
	VECTOR2 q;
	VECTOR2 r;
	VECTOR2 pq;
	int     outcome;
	CParamLine p1;
	CParamLine p2;
}DebugSEC;




class CMedAxisVertex{
	public:
		CMedAxisVertex();
		CMedAxisVertex(VECTOR2 vMid);

		~CMedAxisVertex();

		CCircler m_Circle;
	
		std::vector<VECTOR2*> m_Vertices;
		
		int m_Index;

};


struct comparePoints
{
	bool operator()(const VECTOR2 *p1, const VECTOR2 *p2) const
	{
		return p1 < p2;
	}
};		

//class AdjListElement
//{
//	public:
//	AdjListElement(){};
//	CMedAxisVertex *m_Vertex;
//	std::list<CMedAxisVertex*> m_Adj;
//	
//};
//
//typedef AdjListElement AdjEle;

class cPrimOps : public GeomOps
{

	void findFormingPoints();

	VECTOR2 smallestEnclosingCircle(vector<VECTOR2> &vPts, float *nRad);

	void mergeCircles();

public:
    //void SEC(CMedAxisVertex vV1, CMedAxisVertex vV2);
	CCircler SEC(set<VECTOR2*,comparePoints> &vSet, DebugSEC &info);
};
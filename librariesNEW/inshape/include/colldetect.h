/***************************************************************************
 *   Copyright (C) 2007 by Raphael Münster   *
 *   raphael@reynolds   *
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
#ifndef CCOLLDETECT_H
#define CCOLLDETECT_H

//===================================================
//				     Includes
//===================================================

#include "vector2.h"
#include <vector>
#include <list>

//===================================================
//				Forward declarations
//===================================================

class CApproxCurve;
class CBCNode;

//===================================================
//				Typedefs and structs
//===================================================

typedef vector<CApproxCurve*> CurveVec_t;

typedef pair<int,int> CollPair;

typedef list<CollPair> CollList;


/**
	@author Raphael Münster <raphael@reynolds>
 */
class CCollDetect{
	public:
		CCollDetect();

		~CCollDetect();
		void CollAllPairs(CurveVec_t &vCurves, CollList &lCollisions);

    /*!
        \fn CCollDetect::ProcessPairSimple(Real dRad1, Real dRad2, Real dDist)
		Returns true, when there is a collision
     */
    bool ProcessPairSimple(Real dRad1, Real dRad2, Real dDist)
    {
		if((dRad1+dRad2) > dDist )
			return true;
		else
			return false;
    }//end ProcessPairSimple
    bool TestBVHs(CApproxCurve *pCurve1,CApproxCurve *pCurve2 );
	bool DFSColl(CBCNode *pBig, CBCNode **pTargets, int nTargets);
    void CollAllBVHs(CurveVec_t &vCurves, CollList &lCollisions);

};

#endif

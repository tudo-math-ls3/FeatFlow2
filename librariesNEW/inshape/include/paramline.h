/***************************************************************************
 *   Copyright (C) 2007 by Raphael Münster   *
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


#ifndef _PARAMLINE_H_
#define _PARAMLINE_H_

#include "vector2.h"

class CParamLine
{
public:
	/* constructor */
	CParamLine(){};
	CParamLine(const VECTOR2 &vP0, const VECTOR2 &vP1);

	/* inline member functions */
	inline VECTOR2 GetDir(){return m_vDir;};
	inline VECTOR2 P0(){return m_vP0;};
	inline VECTOR2 P1(){return m_vP1;};

	VECTOR2 Evaluate(double t);

	static int IntersectParamLines2D(CParamLine &p1, CParamLine &p2, double *t1, double *t2);


	/* private member variables */
	VECTOR2 m_vP0;
	VECTOR2 m_vP1;
	VECTOR2 m_vDir;
};

#endif
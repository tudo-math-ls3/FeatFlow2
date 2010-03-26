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

#ifdef WIN32
#pragma once
#endif

#ifndef  _TORUS_H_
#define  _TORUS_H_

#include "vector3.h"
#include "dynamicarray.h"

typedef CDynamicArray<CVector3f> Vert3Array;

class CTorus
{
public:
	CTorus(void);
	~CTorus(void);

	CTorus(int iNumT, int iNumS);

	void InitTorus(int iNumT, int iNumS);

	inline Vert3Array& GetVertices() {return m_Vertices;}
	inline Vert3Array& GetNormals()  {return m_Normals;}
	inline int         Size()        {return m_Vertices.Size();}

private:

	Vert3Array m_Vertices;
	Vert3Array m_Normals;
	
};

#endif

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

#ifndef _MODEL3D_H_
#define _MODEL3D_H_


#ifdef WIN32
#pragma once
#endif

#include "ObjLoader.h"
#include "Triangle3.h"
#include "aabb3.h"


//typedef vector<CVector3f> VertArray;
//typedef vector<tObjFace>  FaceArray;
//typedef CDynamicArray<CVector3f> Vec3Array;
//typedef CDynamicArray<CVector2f> Vec2Array;

typedef CDynamicArray<VECTOR3> VEC3Array;
typedef std::vector<CTriangle3f*> Tri3Array;

class CModel3D {
public:
	CModel3D(void);
	~CModel3D(void);

	void CleanUp(void);

	const VEC3Array& GetVertices() const {return m_pVertices;};
	const Vec2Array& GetTCoords() const {return m_pTexCoords;};
	const Vec3Array& GetVertexNormals() const {return m_pVertexNormals;};
	const FaceArray& GetFaces() const {return m_pFaces;};

	void SetVertices(const VertArray &pVertices);
	void SetFaces(const FaceArray &pFaces);
	void SetNormals(const Vec3Array &pNormals);
	void SetTCoords(const Vec2Array &pTCoords);
	void BuildTriangles();

	inline bool HasUV() const {return m_bHasUV;}

	inline const tObjFace& GetFace(int i) const {return m_pFaces[i];};

	inline CTriangle3f* GetTriangle(int i)
	{
		return m_trTriangles[i];
	}

	inline Tri3Array& GetTriangles() {return m_trTriangles;};

	inline int  NumVertices() const {return m_pVertices.Size();};

	inline const CAABB3f& GetBox() const {return m_bxBox;};


private:

	inline CTriangle3f* MakeTriangle(int i)
	{
		int vertIndex0 = m_pFaces[i].VertexIndex[0];
		int vertIndex1 = m_pFaces[i].VertexIndex[1];
		int vertIndex2 = m_pFaces[i].VertexIndex[2];
		return new CTriangle3f(m_pVertices[vertIndex0], m_pVertices[vertIndex1], m_pVertices[vertIndex2]);
	}



	
	VEC3Array m_pVertices;

	Vec3Array m_pVertexNormals;

	Vec2Array m_pTexCoords;
	
	FaceArray m_pFaces;

	bool      m_bHasUV;

	CAABB3f   m_bxBox;

	Tri3Array m_trTriangles;

};

#endif
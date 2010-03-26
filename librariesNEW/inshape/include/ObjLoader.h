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

#ifndef _OBJLOADER_H_
#define _OBJLOADER_H_

#include <basicmodelloader.h>
#include <vector>


using namespace std;

typedef struct
{

	int VertexIndex[3];
	int TexIndex[3];

}tObjFace;

typedef vector<CVector3f> VertArray;
typedef vector<tObjFace>  FaceArray;
typedef vector<CVector2f> TexCoordArray;
typedef CDynamicArray<CVector3f> Vec3Array;
typedef CDynamicArray<CVector2f> Vec2Array;

class CObjLoader :
	public CBasicModelLoader
{
public:
	CObjLoader(void);
	~CObjLoader(void);

	/* reads the .obj file specified in strFileName */
	void ReadModelFromFile(char *strFileName);
	
	const VertArray& GetVertices() const;
	const FaceArray& GetFaces() const;
	const Vec3Array& GetNormals() const;

	const TexCoordArray& GetTexCoords(void) const;

	
	bool HasUV(void) const;
	

	Vec3Array m_pVertexNormals;


private:

	void ReadVertex(ifstream &in, char strLine[]);

	void ReadFace(ifstream &in, char strLine[]);

	void ReadTexCoord(ifstream &in, char strLine[]);

	void ReadFaceTex(ifstream &in, char strLine[]);

	/* this function computes vertex normals for the current mesh */
	void CalcVertexNormals();

	void CleanUp(void);

	/* private member variables */

	VertArray m_pVertices;

	Vec3Array m_pNormals;

	TexCoordArray m_pTexCoords;

	FaceArray m_pFaces;

	bool m_bUV;

};

#endif
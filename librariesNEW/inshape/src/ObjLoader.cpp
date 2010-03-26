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

#include "ObjLoader.h"
#include <fstream>
#include <sstream>
#include <string>
#include "aabb3.h"

CObjLoader::CObjLoader(void) : m_bUV(false)
{

	CAABB3<Real> aabb;

}//end constructor

CObjLoader::~CObjLoader(void)
{

}//end deconstructor

void CObjLoader::ReadModelFromFile(char *strFileName) 
{

	ifstream in(strFileName);

	char strLine[256];
	string first;

	while(!in.eof())
	{
		
		in>>first;
		
		if(first == string("#"))
		{
			in.getline(strLine,256);
			continue;
		}
		//case: Vertex
		else if(first == string("v"))
			ReadVertex(in,strLine);
		//case: TexCoord
		else if(first == string("vt"))
		{
			ReadTexCoord(in, strLine);
			m_bUV = true;
		}
		//case: Face
		else if(first == string("f"))
		{
			if(!m_bUV)
				ReadFace(in, strLine);
			else
				ReadFaceTex(in, strLine);
		}
		//default
		else
			in.getline(strLine,256);
		

		
	}//end while

	cout <<"Number of vertices: "<<m_pVertices.size()<<endl;

	m_pNormals.Resize(m_pFaces.size());
	m_pVertexNormals.Resize(m_pVertices.size());

	cout <<"Number of faces: "<<m_pFaces.size()<<endl;

	CalcVertexNormals();

}//end ReadModelFromFile

void CObjLoader::ReadVertex(ifstream &in, char strLine[])
{

	CVector3f vec;
	in >> vec.x;
	in >> vec.z;
	in >> vec.y;
	in.getline(strLine,256);
	m_pVertices.push_back(vec);

}//end ReadVertex

void CObjLoader::ReadFace(ifstream &in, char strLine[])
{

	tObjFace Face;

	for(int i = 0; i < 3; i++)
	{
		in >> Face.VertexIndex[i];
		Face.VertexIndex[i]--;
	}

	in.getline(strLine, 256);
	m_pFaces.push_back(Face);

}//end ReadFace

void CObjLoader::ReadTexCoord(ifstream &in, char strLine[])
{
	
	CVector2f vec;
	in >> vec.x;
	in >> vec.y;
	m_pTexCoords.push_back(vec);
	//cout<<m_pTexCoords.size()<<" "<<vec.x<<" "<<vec.y<<endl;
	in.getline(strLine,256);

}//end ReadTexCoord

void CObjLoader::ReadFaceTex(ifstream &in, char strLine[])
{
	tObjFace Face;

	string s;

	basic_string<char> vertIndex;
	basic_string<char> texIndex;
	int vi;
	int ti;
	int Size = m_pFaces.size();

	for(int i = 0; i < 3; i++)
	{
		
		// Format for a face is vertexIndex/texture index vertexIndex/textureIndex vertexIndex/texture index 
		in >> s;
		
		// find separator 
		basic_string<char>::size_type index = s.find("/");

		// extract indices
		vertIndex = s.substr(0,index);
		texIndex = s.substr(index+1,s.size()-1);

		// convert indices from string to int
		istringstream VertStream(vertIndex);
		istringstream TexStream(texIndex);

		VertStream >> vi;
		TexStream  >> ti;

		//assign the values to the face structure
		Face.VertexIndex[i] = vi-1;
		Face.TexIndex[i]    = ti-1;		
		
	}



	//go to next line
	in.getline(strLine, 256);
	m_pFaces.push_back(Face);

}//end ReadFaceTex

const VertArray& CObjLoader::GetVertices() const
{

	return m_pVertices;

}//end GetVertices

const FaceArray& CObjLoader::GetFaces() const
{

	return m_pFaces;

}//end GetVertices

const Vec3Array& CObjLoader::GetNormals() const
{

	return m_pNormals;

}//end GetNormals

void CObjLoader::CalcVertexNormals()
{
	int i,j;

	for(i = 0; i < (int)m_pFaces.size(); i++)
	{

		int vi0 = m_pFaces[i].VertexIndex[0];
		int vi1 = m_pFaces[i].VertexIndex[1];
		int vi2 = m_pFaces[i].VertexIndex[2];

		CVector3f v0 = m_pVertices[vi0];
		CVector3f v1 = m_pVertices[vi1];
		CVector3f v2 = m_pVertices[vi2];

		CVector3f vV1 = CVector3f::createVector(v0, v1);
		CVector3f vV2 = CVector3f::createVector(v0, v2);

		CVector3f vNormal = CVector3f::Cross(vV1, vV2);
		
		vNormal /= vNormal.mag();

		m_pNormals[i] = vNormal;
		
	}//end for

	for(i = 0; i < (int)m_pVertices.size(); i++)
	{

		CVector3f vSum(0,0,0);
		int count = 0;

		for(j = 0; j < (int) m_pFaces.size(); j++)
		{
			if( m_pFaces[j].VertexIndex[0] == i ||
				m_pFaces[j].VertexIndex[1] == i ||
				m_pFaces[j].VertexIndex[2] == i    )
			{

				vSum+=m_pNormals[j];
				count++;

			}//end if
		}//end for

		vSum/=count;
		m_pVertexNormals[i] = vSum * -1.0;


	}//end for

}//end CalcVertexNormals

const TexCoordArray& CObjLoader::GetTexCoords(void) const
{
	return m_pTexCoords;
}//end GetTexCoords

void CObjLoader::CleanUp(void)
{

	m_pNormals.Destroy();
	m_pVertexNormals.Destroy();

}//end CleanUp
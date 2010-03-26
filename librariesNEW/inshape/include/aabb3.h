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

#ifndef _AABB3_H
#define _AABB3_H

//===================================================
//					DEFINITIONS
//===================================================


//===================================================
//					INCLUDES
//===================================================
#include <iostream>
#include <vector>
#include <limits>
#include <vector3.h>
#include <dynamicarray.h>
#include <Triangle3.h>


using namespace std;


template<class T>
class CAABB3
{
public:

	CAABB3(void){};

	CAABB3(const CDynamicArray< CVector3<T> > &Vec3Array);

	~CAABB3(void){};

	void InitBox(const CDynamicArray< CVector3<T> > &Vec3Array);

	void Init(const vector<CTriangle3f*> &vTriangles);

	bool Inside(const CVector3<T> &vQuery) const;

	int LongestAxis() const;

	T MinDistance(const CVector3<T> &vQuery);

	T MinDistanceSqr(const CVector3<T> &vQuery);

	inline T MaxDistance(const CVector3<T> &vQuery) {return (CVector3<T>::createVector(vQuery,m_vUpper)).mag();};

	inline T MaxDistanceSqr(const CVector3<T> &vQuery) {return (CVector3<T>::createVector(vQuery,m_vUpper)).norm2();};

	//inline methods to access vertices
	inline CVector3<T> GetFBL() const {return m_Verts[0];};

	inline CVector3<T> GetBTR() const {return m_Verts[1];};

	inline CVector3<T> GetFBR() const {return CVector3<T>(m_Verts[1].x, m_Verts[0].y, m_Verts[0].z);};

	inline CVector3<T> GetFTR() const {return CVector3<T>(m_Verts[1].x, m_Verts[1].y, m_Verts[0].z);};

	inline CVector3<T> GetFTL() const {return CVector3<T>(m_Verts[0].x, m_Verts[1].y, m_Verts[0].z);};

	inline CVector3<T> GetBBL() const {return CVector3<T>(m_Verts[0].x, m_Verts[0].y, m_Verts[1].z);};

	inline CVector3<T> GetBBR() const {return CVector3<T>(m_Verts[1].x, m_Verts[0].y, m_Verts[1].z);};

	inline CVector3<T> GetBTL() const {return CVector3<T>(m_Verts[0].x, m_Verts[1].y, m_Verts[1].z);};

	inline T Xmin() const {return m_Verts[0].x;};

	inline T Xmax() const {return m_Verts[1].x;};

	inline T Ymin() const {return m_Verts[0].y;};

	inline T Ymax() const {return m_Verts[1].y;};

	inline T Zmin() const {return m_Verts[0].z;};

	inline T Zmax() const {return m_Verts[1].z;};

	inline CVector3<T> GetCenter() const
	{
		CVector3<T> vCenter = m_Verts[0];
		vCenter.x += 0.5 * (m_Verts[1].x - m_Verts[0].x);
	
		vCenter.y += 0.5 * (m_Verts[1].y - m_Verts[0].y);

		vCenter.z += 0.5 * (m_Verts[1].z - m_Verts[0].z);

		return vCenter;
	};

	

	enum
	{
		XAXIS,
		YAXIS,
		ZAXIS
	};

	CVector3<T> m_Verts[2];

	CVector3<T> m_vUpper;

};
template<class T>
bool CAABB3<T>::Inside(const CVector3<T> &vQuery) const
{
	if(  (Xmin() <= vQuery.x && vQuery.x <= Xmax())
	   &&(Ymin() <= vQuery.y && vQuery.y <= Ymax())
	   && (Zmin() <= vQuery.z && vQuery.z <= Zmax()) )
		return true;
	else
		return false;

}

template<class T>
int CAABB3<T>::LongestAxis() const
{
	T rLength = -std::numeric_limits<T>::max();

	int iAxis = -1;

	T lengths[3];

	lengths[0] = fabs(m_Verts[0].x - m_Verts[1].x);
	lengths[1] = fabs(m_Verts[0].y - m_Verts[1].y);
	lengths[2] = fabs(m_Verts[0].z - m_Verts[1].z);

	for(int i = 0; i < 3; i++)
	{
		if(rLength < lengths[i])
		{
			iAxis = i;
			rLength = lengths[i];
		}//end if
	}//end for

	return iAxis;

}//end LongestAxis

template<class T>
CAABB3<T>::CAABB3(const CDynamicArray< CVector3<T> > &Vec3Array)
{
	T MaxX = std::numeric_limits<T>::min();
	T MinX = std::numeric_limits<T>::max();
	T MaxY = std::numeric_limits<T>::min();
	T MinY = std::numeric_limits<T>::max();
	T MaxZ = std::numeric_limits<T>::min();
	T MinZ = std::numeric_limits<T>::max();

	for(int i = 0; i < Vec3Array.size(); i++)
	{

		if(Vec3Array[i].x < MinX)
		{	//assign min index
			m_Verts[0].x = Vec3Array[i].x;
		}

		if(Vec3Array[i].x < MaxX)
		{	//assign max index
			m_Verts[1].x = Vec3Array[i].x;
		}

		if(Vec3Array[i].y < MinY)
		{	//assign min index
			m_Verts[0].y = Vec3Array[i].y;
		}

		if(Vec3Array[i].y < MaxY)
		{	//assign max index
			m_Verts[1].y = Vec3Array[i].y;
		}

		if(Vec3Array[i].z < MinZ)
		{	//assign min index
			m_Verts[0].z = Vec3Array[i].z;
		}

		if(Vec3Array[i].z < MaxZ)
		{	//assign max index
			m_Verts[1].z = Vec3Array[i].z;
		}

	}//end for

}//end constructor

template<class T>
void CAABB3<T>::Init(const vector<CTriangle3f*> &vTriangles)
{
	T MaxX = -std::numeric_limits<T>::max();
	T MinX = std::numeric_limits<T>::max();
	T MaxY = -std::numeric_limits<T>::max();
	T MinY = std::numeric_limits<T>::max();
	T MaxZ = -std::numeric_limits<T>::max();
	T MinZ = std::numeric_limits<T>::max();

	T MinCenter = std::numeric_limits<T>::max();

	
	CVector3<T> vCenter = this->GetCenter();

	for(int i = 0; i < vTriangles.size(); i++)
	{
		CTriangle3f tri = *vTriangles[i];

		for(int j = 0; j < 3; j++)
		{
			CVector3<T> Vec3 = tri.Get(j);
			if(Vec3.x < MinX)
			{	//assign min index
				MinX = Vec3.x;
			}

			if(Vec3.x > MaxX)
			{	//assign max index
				MaxX = Vec3.x;
			}

			if(Vec3.y < MinY)
			{	//assign min index
				MinY = Vec3.y;
			}

			if(Vec3.y > MaxY)
			{	//assign max index
				MaxY = Vec3.y;
			}

			if(Vec3.z < MinZ)
			{	//assign min index
				MinZ = Vec3.z;
			}

			if(Vec3.z > MaxZ)
			{	//assign max index
				MaxZ = Vec3.z;
			}

			T d = CVector3<T>::createVector(Vec3,vCenter).mag();
			if( d < MinCenter)
			{
				m_vUpper = Vec3;
				MinCenter = d;
			}

		}//end for j
	}//end for

	m_Verts[0].x = MinX;
	m_Verts[0].y = MinY;
	m_Verts[0].z = MinZ;

	m_Verts[1].x = MaxX;
	m_Verts[1].y = MaxY;
	m_Verts[1].z = MaxZ;
}//end InitBox

template<class T>
void CAABB3<T>::InitBox(const CDynamicArray< CVector3<T> > &Vec3Array)
{

	T MaxX = -std::numeric_limits<T>::max();
	T MinX = std::numeric_limits<T>::max();
	T MaxY = -std::numeric_limits<T>::max();
	T MinY = std::numeric_limits<T>::max();
	T MaxZ = -std::numeric_limits<T>::max();
	T MinZ = std::numeric_limits<T>::max();

	for(int i = 0; i < Vec3Array.Size(); i++)
	{

		if(Vec3Array[i].x < MinX)
		{	//assign min index
			MinX = Vec3Array[i].x;
		}

		if(Vec3Array[i].x > MaxX)
		{	//assign max index
			MaxX = Vec3Array[i].x;
		}

		if(Vec3Array[i].y < MinY)
		{	//assign min index
			MinY = Vec3Array[i].y;
		}

		if(Vec3Array[i].y > MaxY)
		{	//assign max index
			MaxY = Vec3Array[i].y;
		}

		if(Vec3Array[i].z < MinZ)
		{	//assign min index
			MinZ = Vec3Array[i].z;
		}

		if(Vec3Array[i].z > MaxZ)
		{	//assign max index
			MaxZ = Vec3Array[i].z;
		}

	}//end for

	m_Verts[0].x = MinX;
	m_Verts[0].y = MinY;
	m_Verts[0].z = MinZ;

	m_Verts[1].x = MaxX;
	m_Verts[1].y = MaxY;
	m_Verts[1].z = MaxZ;



}//end InitBox

template<class T>
T CAABB3<T>::MinDistance(const CVector3<T> &vQuery)
{

	CVector3<T> vSol;

	if(vQuery.x < Xmin())
		vSol.x = Xmin()-vQuery.x;
	else if(vQuery.x > Xmax())
		vSol.x = vQuery.x - Xmax();
	else
		vSol.x = 0.0;

	if(vQuery.y < Ymin())
		vSol.y = Ymin()-vQuery.y;
	else if(vQuery.y > Ymax())
		vSol.y = vQuery.y - Ymax();
	else
		vSol.y = 0.0;

	if(vQuery.z < Zmin())
		vSol.z = Zmin()-vQuery.z;
	else if(vQuery.y > Ymax())
		vSol.z = vQuery.z - Zmax();
	else
		vSol.z = 0.0;

	return vSol.mag();

}//end MinDistance

template<class T>
T CAABB3<T>::MinDistanceSqr(const CVector3<T> &vQuery)
{

	CVector3<T> vSol;

	if(vQuery.x < Xmin())
		vSol.x = Xmin()-vQuery.x;
	else if(vQuery.x > Xmax())
		vSol.x = vQuery.x - Xmax();
	else
		vSol.x = 0.0;

	if(vQuery.y < Ymin())
		vSol.y = Ymin()-vQuery.y;
	else if(vQuery.y > Ymax())
		vSol.y = vQuery.y - Ymax();
	else
		vSol.y = 0.0;

	if(vQuery.z < Zmin())
		vSol.z = Zmin()-vQuery.z;
	else if(vQuery.y > Ymax())
		vSol.z = vQuery.z - Zmax();
	else
		vSol.z = 0.0;

	return vSol.norm2();

}//end MinDistanceSqr

typedef CAABB3<float> CAABB3f;
typedef CAABB3<double> CAABB3d;

#endif
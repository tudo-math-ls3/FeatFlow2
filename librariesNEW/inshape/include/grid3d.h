/****************************************************************************
**
** Copyright (C) 2005-2007 Trolltech ASA. All rights reserved.
**
** This file is part of the example classes of the Qt Toolkit.
**
** This file may be used under the terms of the GNU General Public
** License version 2.0 as published by the Free Software Foundation
** and appearing in the file LICENSE.GPL included in the packaging of
** this file.  Please review the following information to ensure GNU
** General Public Licensing requirements will be met:
** http://www.trolltech.com/products/qt/opensource.html
**
** If you are unsure which license is appropriate for your use, please
** review the following information:
** http://www.trolltech.com/products/qt/licensing.html or contact the
** sales department at sales@trolltech.com.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
**
****************************************************************************/

#ifdef WIN32
#pragma once
#endif

#ifndef _GRID3D_H
#define _GRID3D_H

//===================================================
//					DEFINITIONS
//===================================================

#define GRIDX 15
#define GRIDY 15
#define GRIDZ 15

#define GRIDXYZ (GRIDX * GRIDY * GRIDZ)

//===================================================
//					INCLUDES
//===================================================
#include <iostream>
#include "vector3.h"


//===================================================
//					STRUCTS
//===================================================
typedef struct
{
	Real distance;
	bool Inside;
	Real m_LUB;
}Grid3_t;

template<class T>
class CGrid3D
{
public:
	CGrid3D(void);
	~CGrid3D(void);

	void InitGrid(const CVector3<T> &vBL, const CVector3<T> &vTR, const CVector3<T> &vBackTR);

	inline int GetSizeX() {return m_iSizeX;};
	inline int GetSizeY() {return m_iSizeY;};
	inline int GetSizeZ() {return m_iSizeZ;};

   inline CVector3<T>& operator() (int row, int col, int slice)
   {
       return m_pVertices[m_iSizeX * row + col + slice * m_iSizeX * m_iSizeY];
   }
  
   inline CVector3<T>& operator() (int row, int col, int slice) const
   {
       return m_pVertices[m_iSizeX * row + col + slice * m_iSizeX * m_iSizeY];
   }

  inline CVector3<T>& operator() (int index) const
   {
       return m_pVertices[index];
   }

  inline CVector3<T>* GetVertices() {return m_pVertices;};


    /*!
        \fn CGrid::SetDistance(int i, int j, int k, Real d)
     */
    inline void SetDistance(int i, int j, int k, Real d)
    {
        m_pGrid[i][j][k].distance = d;
    }

    /*!
        \fn CGrid::Distance(int i, int j, int k)
     */
    inline Real Distance(int i, int j, int k) const
    {
        return m_pGrid[i][j][k].distance;
    }

	Grid3_t     m_Grid[GRIDXYZ];

private:
	int m_iSizeX;
	int m_iSizeY;
	int m_iSizeZ;
	CVector3<T> *m_pVertices;
	Grid3_t     m_pGrid[GRIDX][GRIDY][GRIDZ];
	

};


template<class T>
CGrid3D<T>::CGrid3D(void) : m_iSizeX(GRIDX), m_iSizeY(GRIDY), m_iSizeZ(GRIDZ), m_pVertices(0)
{
}

template<class T>
CGrid3D<T>::~CGrid3D(void)
{
	if(m_pVertices)
	{
		delete[] m_pVertices;
		m_pVertices = NULL;
	}

	//if(m_pGrid)
	//{
	//	delete[] m_pGrid;
	//	m_pGrid = NULL;
	//}

}

template<class T>
void CGrid3D<T>::InitGrid(const CVector3<T> &vBL, const CVector3<T> &vTR, const CVector3<T> &vBackTR)
{

	int allVertices = m_iSizeX * m_iSizeY * m_iSizeZ;

	int iSlice = m_iSizeX * m_iSizeY; 

	m_pVertices = new CVector3<T>[allVertices];

	//m_pGrid = new Grid3_t[allVertices];

	Real LengthX = fabs(vTR.x - vBL.x);
	Real LengthY = fabs(vTR.y - vBL.y);
	Real LengthZ = -fabs(vTR.z - vBackTR.z);

	Real xInc = LengthX / Real(m_iSizeX-1);
	Real yInc = LengthY / Real(m_iSizeY-1);
	Real zInc = LengthZ / Real(m_iSizeZ-1);

	int count = 0;
	for(int i = 0; i < m_iSizeX; i++)
	{
		for(int j = 0; j < m_iSizeY; j++)
		{
			for(int z = 0; z <m_iSizeZ; z++)
			{
				m_pVertices[i*m_iSizeX + j + z * iSlice] = CVector3<T>(vBL.x + i*xInc, vBL.y + j*yInc, vBL.z + z * zInc);
				count = i*m_iSizeX + j + z * iSlice;
			}//end for
		}//end for
	}//end for
}//end InitGrid

typedef CGrid3D<Real> GRID3D;

#endif
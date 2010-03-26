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
#ifndef GRID_H
#define GRID_H

/**
The class CGrid holds the global computational grid. 

	@author Raphael Münster <raphael@reynolds>
*/

//===============INCLUDES=================
#include "vector2.h"


#define GRID_POINTS_X 512

#define GRID_POINTS_Y 512

#define GRIDSIZEX GRID_POINTS_X
#define GRIDSIZEY GRID_POINTS_Y

#define ALLPOINTS GRIDSIZEX * GRIDSIZEY

#define DOUBLEALLPOINTS 2 * ALLPOINTS

typedef struct
{
	Real distance;
	Real newton;
	Real m_LUB;
}tGrid;

class CGrid{
public:
    CGrid();

    ~CGrid();
	void InitGrid(const VECTOR2 &vBL, const VECTOR2 &vTR);

    /*!
        \fn CGrid::GetX()
     */
    inline int GetX() const
    {
        return m_SizeX;
    }

    /*!
        \fn CGrid::GetY()
     */
    inline int GetY() const
    {
		return m_SizeY;
    }

    /*!
        \fn CGrid::SetDistance(int i, int j, Real d)
     */
    inline void SetDistance(int i, int j, Real d)
    {
        m_Grid[i*GRID_POINTS_X + j].distance = d;
    }

    inline void SetDistance(int i, Real d)
    {
        m_Grid[i].distance = d;
    }

    /*!
        \fn CGrid::Distance(int i, int j)
     */
    inline Real Distance(int i, int j) const
    {
        return m_Grid[i*GRID_POINTS_X + j].distance;
    }

    /*!
        \fn CGrid::Point(int i, int j)
     */
    inline VECTOR2 Point(int i, int j) 
    {
        return m_pVertices[i*GRID_POINTS_X + j];
    }

   inline CVector2f TCoord(int i, int j) const
   {
        return m_pTCoords[i*GRID_POINTS_X + j];
   }

   inline void SetTCoord(int i, int j, float texCoord)
   {
        m_pTCoords[i*GRID_POINTS_X + j] = CVector2f(texCoord,0.0);
   }

   inline VECTOR2& operator() (int row, int col)
   {
       return m_pVertices[GRID_POINTS_X * row + col];
   }
  
   inline VECTOR2 operator() (int row, int col) const
   {
       return m_pVertices[GRID_POINTS_X * row + col];
   }

   inline VECTOR2& operator() (int iIndex)
   {
	   return m_pVertices[iIndex];
   }//end operator

   inline int TriangleCount()
   {
	   return (m_SizeY-1) * 2 * (m_SizeX-1);
   }

    
	void MakeTCoords(void);
	void BuildIndices(void);
    void PerturbMesh(double dPercent);

	CVector2f	 m_pTCoords[ALLPOINTS];
	unsigned int m_pIndices[DOUBLEALLPOINTS];
	VECTOR2      m_pVertices[ALLPOINTS];
	double m_nMaxDist;
	
private:

	void AddTriangle(int a, int b, int c);

	tGrid		 m_Grid[ALLPOINTS];

	
protected:
	
    int m_SizeX;
	int m_SizeY;
	int m_Count;
	
};
#endif

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
#include <math.h>
#include <grid.h>
#include <mathglobals.h>

CGrid::CGrid()
{
	m_SizeX = GRID_POINTS_X;
	m_SizeY = GRID_POINTS_Y;
	m_Count = 0;

	BuildIndices();

}//end constructor


CGrid::~CGrid()
{
}


/*!
    \fn CGrid::InitGrid()
 */
void CGrid::InitGrid(const VECTOR2 &vBL, const VECTOR2 &vTR)
{

	Real LengthX = fabs(vTR.x - vBL.x);
	Real LengthY = fabs(vTR.y - vBL.y);

	Real xInc = LengthX / Real(GRID_POINTS_X-1);
	Real yInc = LengthY / Real(GRID_POINTS_Y-1);

	for(int i = 0; i < GRID_POINTS_X; i++)
	{
		for(int j = 0; j < GRID_POINTS_Y; j++)
		{
			m_pVertices[i*GRID_POINTS_X + j] = VECTOR2(vBL.x + i*xInc, vBL.y + j*yInc);
			m_Grid[i*GRID_POINTS_X + j].distance = 0.0;
			m_Grid[i*GRID_POINTS_X + j].newton = 0.0;
		}//end for
	}//end for

}//end InitGrid


/*!
    \fn CGrid::PerturbMesh(double dPercent)
 */
void CGrid::PerturbMesh(double dPercent)
{
	Real h = (Point(0,0) - Point(0,1)).mag();

	Real hp = h * dPercent;

	for(int j = 0; j < GRID_POINTS_X; j++)
	{
		for(int i = 0; i < GRID_POINTS_Y; i++)
		{
			VECTOR2 dir(hp*cos(double(rand()%360)), hp*sin(double(rand()%360)));
			m_pVertices[j*GRID_POINTS_X + i]+=dir;

		}
	}


}

void CGrid::MakeTCoords(void)
{
	m_nMaxDist = 0.0;
	
	for(int i = 0; i < GRID_POINTS_X; i++)
	{
				
		for(int j = 0; j < GRID_POINTS_Y; j++)
		{
	
			if(Distance(i,j) > m_nMaxDist) m_nMaxDist = Distance(i,j);

		}//end for j
	}//end for i

	//cout<<"max dist "<<m_nMaxDist<<endl;

	for(int i = 0; i < GRID_POINTS_X; i++)
	{
				
		for(int j = 0; j < GRID_POINTS_Y; j++)
		{
			float sc = Distance(i,j) / m_nMaxDist;
			float t1 = (sc < 0.0001) ? 0.0 : sc;
			SetTCoord(i,j,t1);
		}//end for
	}//end for

}

void CGrid::BuildIndices(void)
{
	int i,j;

	int count = 0;

	for (i = 0; i < m_SizeX - 1; i++)
	{
		for (j = 0; j < m_SizeY; j++)
		{
			m_pIndices[count] = i*m_SizeX+j; 
			count++;
			m_pIndices[count] = (i+1)*m_SizeX+j; 
			count++;
		}//end for
	}//end for

}//end BuildIndices

void CGrid::AddTriangle(int a, int b, int c)
{

	this->m_pIndices[m_Count] = a;
	this->m_pIndices[m_Count+1] = b;
	this->m_pIndices[m_Count+2] = c;

	m_Count+=3;

}//end AddTriangle


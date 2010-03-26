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

#include <math.h>
#include <torus.h>


CTorus::CTorus(void)
{

}

CTorus::CTorus(int iNumT, int iNumS)
{

	InitTorus(iNumT, iNumS);

}

void CTorus::InitTorus(int iNumT, int iNumS)
{

	int i,j,k;
	double s, t, x, y, z, PI2;

	PI2 =  2 * (double)M_PI;

	double R = 0.6;
	double r = 0.2;

	m_Vertices.Resize(iNumS * (iNumT+1) * 2);
	m_Normals.Resize(iNumS * (iNumT+1) * 2);

	int iCount = 0;

	for(i = 0; i < iNumS;i++)
	{

		for(j = 0; j <= iNumT; j++)
		{

			t = static_cast<double>(j)/static_cast<double>(iNumT);

			for(k = 1; k >= 0; k--)
			{
				s = (i + k); 
				s/= iNumS;

				x = (R + r * cos(s*PI2))*cos(t*PI2);
				z = (R + r * cos(s*PI2))*sin(t*PI2);
				y = r * sin(s * PI2);

				double nx, ny, nz;

				nx = R * r * cos(s*PI2) * cos(t*PI2) + r * r * cos(s*PI2) * cos(s*PI2) * cos(t*PI2);
				nz = R * r * cos(s*PI2) * sin(t*PI2) + r * r * cos(s*PI2) * cos(s*PI2) * sin(t*PI2);
				ny = R * r * sin(s*PI2) + r * r * cos(s*PI2) * sin(s*PI2);

				CVector3f vec(nx, ny, nz);
				vec /= vec.mag();

				m_Vertices[iCount] = CVector3f(x,y,z);
				m_Normals[iCount]  = vec;
				
				iCount++;

			}//end for k
		}//end for j
	}//end for i
	
	
}

CTorus::~CTorus(void)
{
}

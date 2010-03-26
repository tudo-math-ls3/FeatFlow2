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

#ifdef WIN32
#pragma once
#endif

#ifndef _RAY3D_H
#define _RAY3D_H

//===================================================
//					DEFINITIONS
//===================================================


//===================================================
//					INCLUDES
//===================================================
#include <iostream>
#include "vector3.h"

template <class T>
class CRay3
{
public:

	/* constructor */
	CRay3(void);

	//pass the origin and direction as parameters
	CRay3(const CVector3<T> &vOrig, const CVector3<T> &vDir);

	/* deconstructor */
	~CRay3(void){};

	//member variables for origin and direction
	CVector3<T> m_vOrig;
	CVector3<T> m_vDir;

};

template<class T>
CRay3<T>::CRay3()
{
}//end constructor

template<class T>
CRay3<T>::CRay3(const CVector3<T> &vOrig, const CVector3<T> &vDir) : m_vOrig(vOrig), m_vDir(vDir) 
{
}//end constructor

typedef CRay3<float> CRay3f;
typedef CRay3<double> CRay3d;

#endif
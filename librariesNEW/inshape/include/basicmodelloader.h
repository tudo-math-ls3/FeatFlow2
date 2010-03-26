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

#ifndef _BASICLOADER_H_
#define _BASICLOADER_H_


#include <cstdlib>
#include <iostream>
#include "vector3.h"
#include "vector2.h"
#include "dynamicarray.h"


class CBasicModelLoader
{
public:

	CBasicModelLoader(void);
	~CBasicModelLoader(void);
	
	virtual void ReadModelFromFile(char *strFileName) = 0;

protected:
	FILE *m_FilePointer;


};

#endif
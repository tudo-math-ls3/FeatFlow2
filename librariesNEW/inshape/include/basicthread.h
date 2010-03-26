/***************************************************************************
 *   Copyright (C) 2007 by Raphael Münster   *
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

#ifndef _BASICTHREAD_H_
#define _BASICTHREAD_H_

#ifdef WIN32
#include <iostream>
#include <windows.h>
#endif

class CBasicThread
{
public:

	//constructor
	CBasicThread(void);

	//deconstructor
	virtual ~CBasicThread(void);

	//starts the thread and executes the threadproc
	void begin();

	//terminates the thread
	void end();

	//this function has to be defined in the derived classes
	virtual void ThreadProc(void);

	//releases the resources used by the thread
	virtual void CleanUp();

	//returns whether the thread is still active
	bool IsActive();

protected:

	bool m_bTerminationFlag;

#ifdef WIN32
	HANDLE m_hThreadHandle;
#endif

};

#endif

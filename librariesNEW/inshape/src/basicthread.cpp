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

#include "basicthread.h"

void WindowsThreadFunction(CBasicThread *pThread)
{

	pThread->ThreadProc();

}//end WindowsThreadFunction

CBasicThread::CBasicThread(void) : m_bTerminationFlag(false), m_hThreadHandle(0)
{
}//end constructor

CBasicThread::~CBasicThread(void)
{
}//end deconstructor

//this function has to be defined in the derived classes
void CBasicThread::ThreadProc(void)
{
}//end ThreadProc

//starts the thread and executes the threadproc
void CBasicThread::begin()
{

	m_hThreadHandle = CreateThread(NULL, NULL, (LPTHREAD_START_ROUTINE)WindowsThreadFunction, this, NULL, NULL);

}//end begin

//terminates the thread
void CBasicThread::end()
{

	m_bTerminationFlag = true;

	WaitForSingleObject(m_hThreadHandle, INFINITE);

}//end end

//releases the resources used by the thread
void CBasicThread::CleanUp()
{
}//end CleanUp

//returns whether the thread is still active
bool CBasicThread::IsActive()
{
	return !m_bTerminationFlag;
}//end IsActive

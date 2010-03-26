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

#ifndef _DYNAMICARRAY_H_
#define _DYNAMICARRAY_H_

#ifdef WIN32
#pragma once
#endif

#include <cstdlib>

template<class T>
class CDynamicArray
{
public:

	/* constructors */
	CDynamicArray(void);
	CDynamicArray(int Size);
	CDynamicArray(const CDynamicArray &DynArray);

	/* public member functions */
	void Resize(int Size);
	void Resize(T *data, int Size);
	int Size() const ;
	T&  Get(int number) const;
	T&  operator[](int number) const;
	CDynamicArray& operator=(const CDynamicArray &DynArray);
	void Destroy();
	bool Empty();

	/* destructor */
	~CDynamicArray(void);



	/* private member variables */
private:
	int m_iSize;
	T   *m_Content;
	

};

template<class T>
void CDynamicArray<T>::Destroy()
{

	if(m_Content)
	{
		delete[] m_Content;
		m_Content = NULL;
	}

}

template<class T>
CDynamicArray<T>::CDynamicArray(void) : m_Content(NULL), m_iSize(0)
{

}//end constructor

template<class T>
CDynamicArray<T>::~CDynamicArray(void)
{
	if(m_Content)
	{
		delete[] m_Content;
		m_Content = NULL;
	}
}//end constructor

template<class T>
CDynamicArray<T>::CDynamicArray(int Size) : m_iSize(Size)
{

	m_Content = new T[Size];
	
}//end constructor

template<class T>
CDynamicArray<T>::CDynamicArray(const CDynamicArray<T> &DynArray)
{

	m_iSize = DynArray.Size();

	m_Content = new T[m_iSize];

	for(int i = 0; i < m_iSize; i++)
	{

		m_Content[i] = DynArray.Get(i);

	}//end for

}//end constructor

template<class T>
int CDynamicArray<T>::Size() const 
{

	return m_iSize;

}//end Size

template<class T>
void CDynamicArray<T>::Resize(int Size)
{

	T *newContents = new T[Size];

	for(int i = 0; i < m_iSize && i < Size;i++)
		newContents[i] = m_Content[i];

	if(m_Content)
	{
		delete[] m_Content;
		m_Content = NULL;
	}

	m_Content = newContents;
	m_iSize = Size;

}//end Resize

template<class T>
void CDynamicArray<T>::Resize(T *data, int Size)
{

	if(m_Content)
	{
		delete[] m_Content;
		m_Content = NULL;
	}

	m_Content = data;
	m_iSize = Size;

}//end Resize


template<class T>
bool CDynamicArray<T>::Empty()
{

	return (m_iSize == 0);

}//end operator

template<class T>
T&  CDynamicArray<T>::Get(int number) const
{

	if(number >= 0 && number < m_iSize)
		return m_Content[number];

	return m_Content[0];

}//end Get

template<class T>
T&  CDynamicArray<T>::operator[](int number) const
{
	return m_Content[number];
}//end operator


template<class T>
CDynamicArray<T>& CDynamicArray<T>::operator=(const CDynamicArray<T> &DynArray)
{

	Destroy();

	m_iSize = DynArray.Size();

	m_Content = new T[m_iSize];

	for(int i = 0; i < m_iSize; i++)
	{

		m_Content[i] = DynArray.Get(i);

	}//end for

	return *this;

}//end operator



#endif

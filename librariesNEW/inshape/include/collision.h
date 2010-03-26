/***************************************************************************
 *   Copyright (C) 2006 by Raphael Münster   *
 *   raphael@cortez   *
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

#include <list>
#include <vector>
#include "cBinSTree.h"


/* forward declaration */
class cApproxCurve;
class cOBB;

class cIntervalTreeNode
{
	public:
	/* constructor */
	cIntervalTreeNode();
	
	/* deconstructor */
	~cIntervalTreeNode();
	
};

class cCollisionDetector
{
	
	public:
	/* constructor */
	cCollisionDetector(vector<cApproxCurve*>& pCurves);

	/* deconstructor */
	~cCollisionDetector();
	
	/* public member functions */
	void BuildIntervalTree();
	
	private:
		vector<double> m_BoxValues;
		
	
};



		



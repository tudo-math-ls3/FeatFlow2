/***************************************************************************
 *   Copyright (C) 2006 by Raphael Mnster   *
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

#if !defined _CQUADTREE_H_
#define _CQUADTREE_H_


#include <iostream>
#include <vector2.h>
#include "bspline.h"
#include "boxes.h"

/* namespace for std lib */
using namespace std;

/* names for the nodes */
enum eNodes
{
	TOP_LEFT,		//0
	TOP_RIGHT,		//1 ...
	BOTTOM_LEFT,
	BOTTOM_RIGHT,
	ON
};

/* forward declaration */

class cQuadTree;

//public:
//	cQuadNode();
//	cSBB	  *box;
//	bool	  isSubdivided;
//	cQuadNode *parent;
//	cQuadNode *children[4];
//	int		  getItemCount();


/* tree node object for a binary tree */

class cQuadTreeNode
{

public:
	

	/* contructor */
	cQuadTreeNode(const cVector2f &bottomLeft, const cVector2f &topRight);
	cQuadTreeNode(const cVector2f &bottomLeft, const cVector2f &topRight, const cVector2f &vKey);

	/* destructor */
	~cQuadTreeNode();

	/* member functions */
	inline bool isLeaf(){return !m_bIsSubdivided;};
	
	void getBoxes(){for(int i = 0; i <(int) m_iBoxes.size();i++)cout<<"Box: "<<m_iBoxes[i]<<endl;};
		
	/* get the child nodes */
	cQuadTreeNode* getTopLeft();
	cQuadTreeNode* getTopRight();
	cQuadTreeNode* getBottomLeft();
	cQuadTreeNode* getBottomRight();
	cQuadTreeNode* getChild(int iChild);
	int			   getItemCount();
	
	cVector2f	   m_vBottomLeft;
	cVector2f	   m_vTopRight;
	cVector2f	   m_vCenter;
	cVector2f	   m_vKey;
	cQuadTreeNode    *m_Children[4];




	/* the main BinSTree class is a friend */
	friend class cQuadTree;

private:
	
	bool			 m_bIsSubdivided;
	vector<int>		 m_iBoxes;
	
};



class cQuadTree
{

public:

	/* constructor */
	cQuadTree();
	cQuadTree(cQuadTreeNode *root);
	cQuadTree(cQuadTreeNode *root, cBSpline *pCurve);

	/* deconstructor */
	~cQuadTree();

	/* member functions */
	void inOrder(cQuadTreeNode *tNode);
	void assignBoxesToNode(cQuadTreeNode *tNode);
	void subdivideNode(cQuadTreeNode *tNode);
	void subdivideNode(cQuadTreeNode *tNode, const cVector2f &vKey);
	void deleteSubTree(cQuadTreeNode *tNode);
	void destroyTree();
	int iCompare(const cVector2f &vQuery, const cVector2f &vKey);
	void insertNode(const cVector2f &vCenter);
	int find(const cVector2f &vQuery);
	cQuadTreeNode* findNode(const cVector2f &vQuery);
	cQuadTreeNode* getRoot(){return m_Root;};

private:
	/* private members */
	cQuadTreeNode *m_Root;
	cQuadTreeNode *m_Current;
	cBSpline *m_pCurve;
	
};




#endif


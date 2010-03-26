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

#include "quadtree.h"


cQuadTreeNode::cQuadTreeNode(const cVector2f &bottomLeft, const cVector2f &topRight)
{
	
	m_vBottomLeft   = bottomLeft;
	m_vTopRight     = topRight;
	m_bIsSubdivided = false;

	double dWidth  = m_vTopRight.x - m_vBottomLeft.x;
	double dHeight = m_vTopRight.y - m_vBottomLeft.y;
	double dWidth2 = dWidth/2.0;
	double dHeight2= dHeight/2.0;

	
	/* calculate center */
	m_vCenter = cVector2f(m_vBottomLeft.x + dWidth2, m_vBottomLeft.y + dHeight2);

	m_vKey    = m_vCenter;

	for(int i = 0; i < 4; i++)
		m_Children[i] = NULL;

}

cQuadTreeNode::cQuadTreeNode(const cVector2f &bottomLeft, const cVector2f &topRight, const cVector2f &vKey)
{
	
	m_vBottomLeft   = bottomLeft;
	m_vTopRight     = topRight;
	m_bIsSubdivided = false;

	double dWidth  = m_vTopRight.x - m_vBottomLeft.x;
	double dHeight = m_vTopRight.y - m_vBottomLeft.y;
	double dWidth2 = dWidth/2.0;
	double dHeight2= dHeight/2.0;

	
	/* calculate center */
	m_vCenter = cVector2f(m_vBottomLeft.x + dWidth2, m_vBottomLeft.y + dHeight2);

	m_vKey    = vKey;

	for(int i = 0; i < 4; i++)
		m_Children[i] = NULL;

}

cQuadTreeNode::~cQuadTreeNode()
{


	for(int i = 0; i < 4; i++)
	{
		if(m_Children[i] != NULL)
		{
			delete m_Children[i];
			m_Children[i] = NULL;
		}
	}

}

cQuadTreeNode* cQuadTreeNode::getChild(int iChild)
{
	if((iChild < 0) || (iChild > 3))
		return NULL;

	if(m_bIsSubdivided)
		return m_Children[iChild];
	else
		return NULL;

}//end getChild

cQuadTreeNode* cQuadTreeNode::getTopLeft()
{
	if(m_bIsSubdivided)
		return m_Children[TOP_LEFT];
	else
		return NULL;
}


cQuadTreeNode* cQuadTreeNode::getTopRight()
{
	if(m_bIsSubdivided)
		return m_Children[TOP_RIGHT];
	else
		return NULL;
}


cQuadTreeNode* cQuadTreeNode::getBottomLeft()
{
	if(m_bIsSubdivided)
		return m_Children[BOTTOM_LEFT];
	else
		return NULL;
}


cQuadTreeNode* cQuadTreeNode::getBottomRight()
{
	if(m_bIsSubdivided)
		return m_Children[BOTTOM_RIGHT];
	else
		return NULL;
}



/* quadtree class definition */
cQuadTree::cQuadTree()
{

	m_Root    = NULL;
	m_Current = NULL;
	m_pCurve  = NULL;

}

cQuadTree::cQuadTree(cQuadTreeNode *root)
{
	m_Root	  = root;
	m_Current = root;
	m_pCurve  = NULL;
}

cQuadTree::cQuadTree(cQuadTreeNode *root, cBSpline *pCurve) : m_pCurve(pCurve) 
{
	m_Root	  = root;
	m_Current = root;
}

cQuadTree::~cQuadTree()
{
	if(m_pCurve)
	{
		delete m_pCurve;
		m_pCurve = NULL;
	}
}

void cQuadTree::destroyTree()
{
	deleteSubTree(m_Root);
}

int cQuadTree::find(const cVector2f &vQuery)
{
	/* search tree for item */
	m_Current = findNode(vQuery);

	/* item found */
	if(m_Current != NULL)
	{
		return 1;
	}
	/* not found */
	else
	{
		return 0;
	}
}

cQuadTreeNode* cQuadTree::findNode(const cVector2f &vQuery)
{

	/* start at root */
	cQuadTreeNode *t = m_Root;
	
	/* parent node */
	cQuadTreeNode *parent = NULL;

	/* while there is a valid node pointer */
	while(t != NULL)
	{
		int direction = iCompare(vQuery, t->m_vKey);
		switch(direction)
		{
		case TOP_RIGHT:
			parent = t;
			t = t->getTopRight();
			break;
		case TOP_LEFT:
			parent = t;
			t = t->getTopLeft();
			break;
		case BOTTOM_RIGHT:
			parent = t;
			t = t->getBottomRight();
			break;
		case BOTTOM_LEFT:
			parent = t;
			t = t->getBottomLeft();
			break;
		case ON:
			parent = t;
			t = NULL;
			break;
		}
		///* check where to go */
		//if((vQuery.x >= t->m_vKey.x) && (vQuery.y >= t->m_vKey.y))
		//{
		//	parent = t;
		//	t = t->getTopRight();
		//}//end if
		//else if((vQuery.x < t->m_vKey.x) && (vQuery.y >= t->m_vKey.y))
		//{
		//	parent = t;
		//	t = t->getTopLeft();
		//}//end if
		//else if((vQuery.x >= t->m_vKey.x) && (vQuery.y < t->m_vKey.y))
		//{
		//	parent = t;
		//	t = t->getBottomRight();
		//}//end if
		//else if((vQuery.x < t->m_vKey.x) && (vQuery.y < t->m_vKey.y))
		//{
		//	parent = t;
		//	t = t->getBottomLeft();
		//}//end if
	}//end while

	/* return node */
    return parent;
	
}

int cQuadTree::iCompare(const cVector2f &vQuery, const cVector2f &vKey)
{
	int xQuery = int(1000 * vQuery.x);
	int yQuery = int(1000 * vQuery.y);
	int xKey   = int(1000 * vKey.x);
	int yKey   = int(1000 * vKey.y);
	/* check where to go */
	if((xQuery > xKey) && (yQuery > yKey))
	{
		return TOP_RIGHT;
	}//end if

	else if((xQuery < xKey) && (yQuery > yKey))
	{
		return	TOP_LEFT;	
	}//end if

	else if((xQuery > xKey) && (yQuery < yKey))
	{
		return BOTTOM_RIGHT;		
	}//end if

	else if((xQuery < xKey) && (yQuery < yKey))
	{
		return BOTTOM_LEFT;
	}//end if

	return ON;
}

void cQuadTree::subdivideNode(cQuadTreeNode *tNode)
{
	if(tNode == NULL)
		return;

	/* node is subdivided */
	tNode->m_bIsSubdivided = true;

	/* calculate center */
	cVector2f center = tNode->m_vCenter;

	/* create new nodes */
	cQuadTreeNode *nTopLeft = new cQuadTreeNode(cVector2f(tNode->m_vBottomLeft.x, center.y),
										        cVector2f(center.x, tNode->m_vTopRight.y));

	cQuadTreeNode *nTopRight = new cQuadTreeNode(center, tNode->m_vTopRight);

	cQuadTreeNode *nBottomRight = new cQuadTreeNode(cVector2f(center.x, tNode->m_vBottomLeft.y),
													cVector2f(tNode->m_vTopRight.x, center.y));

	cQuadTreeNode *nBottomLeft  = new cQuadTreeNode(tNode->m_vBottomLeft, center);

	
	
	/* assign new children */
	tNode->m_Children[TOP_LEFT]     = nTopLeft;
	tNode->m_Children[TOP_RIGHT]    = nTopRight;
	tNode->m_Children[BOTTOM_LEFT]  = nBottomLeft;
	tNode->m_Children[BOTTOM_RIGHT] = nBottomRight;

}//end subdivideNode

void cQuadTree::subdivideNode(cQuadTreeNode *tNode, const cVector2f &vKey) 
{
	if(tNode == NULL)
		return;

	if( (vKey.x < tNode->m_vBottomLeft.x) || (vKey.x > tNode->m_vTopRight.x) ||
		(vKey.y < tNode->m_vBottomLeft.y) || (vKey.y > tNode->m_vTopRight.y) )
	{
		return;
	}
	
	tNode->m_vKey = vKey;	

	/* node is subdivided */
	tNode->m_bIsSubdivided = true;
	
	/* create new nodes */
	/* new node top left */
	cQuadTreeNode *nTopLeft = new cQuadTreeNode(cVector2f(tNode->m_vBottomLeft.x, vKey.y),
												cVector2f(vKey.x, tNode->m_vTopRight.y));
	/* new node top right */
	cQuadTreeNode *nTopRight = new cQuadTreeNode(vKey, tNode->m_vTopRight);

	/* new node bottom right */
	cQuadTreeNode *nBottomRight = new cQuadTreeNode(cVector2f(vKey.x, tNode->m_vBottomLeft.y),
												    cVector2f(tNode->m_vTopRight.x, vKey.y));
	/* new node bottom left */
	cQuadTreeNode *nBottomLeft  = new cQuadTreeNode(tNode->m_vBottomLeft, vKey);
	

	/* assign new children */
	tNode->m_Children[TOP_LEFT]     = nTopLeft;
	tNode->m_Children[TOP_RIGHT]    = nTopRight;
	tNode->m_Children[BOTTOM_LEFT]  = nBottomLeft;
	tNode->m_Children[BOTTOM_RIGHT] = nBottomRight;
}

void cQuadTree::deleteSubTree(cQuadTreeNode *tNode)
{
	if(tNode != NULL)
	{
		/* deallocate the four children */
		for(int i = 0; i < 4; i++)
		{
			deleteSubTree(tNode->m_Children[i]);
			tNode->m_Children[i] = NULL;
		}
		/* delete the current node */
		delete tNode;
		/* set pointer to NULL */
		tNode = NULL;
	}
}//end deleteSubTree

void cQuadTree::insertNode(const cVector2f &vCenter)
{
	cQuadTreeNode *tParent = findNode(vCenter);
	subdivideNode(tParent, vCenter);
}//end insertNode

void cQuadTree::assignBoxesToNode(cQuadTreeNode *tNode)
{

	if(tNode == NULL)
		return;

	for(int c = 0; c < 4; c++)
		assignBoxesToNode(tNode->m_Children[c]);

	if(tNode->isLeaf())
	{
	

		int s =(int) m_pCurve->m_Boxes.size();
		int i,j;
		for(i = 0; i < s; i++)
		{
			
			cSBB* box = m_pCurve->m_Boxes[i];
			//cVector2f bounds[4];
			
			//bounds[0] = tNode->m_vBottomLeft;
			//bounds[1] = cVector2f(tNode->m_vTopRight.x, tNode->m_vBottomLeft.y); 
			//bounds[2] = tNode->m_vTopRight;
			//bounds[3] = cVector2f(tNode->m_vBottomLeft.x, tNode->m_vTopRight.y);

		
			for(j  = 0; j < 4; j++)
			{
				cQuadTreeNode *result = findNode(box->m_Vertices[j]);
				if( (tNode == result) || (box->intersection(tNode->m_vTopRight.x, tNode->m_vBottomLeft.x, tNode->m_vTopRight.y, tNode->m_vBottomLeft.y)) )	
				{
					tNode->m_iBoxes.push_back(i);
					break;
				}
						
			}//end for j
		}//end for i
	}//end if

}//end assignBoxesToNode

void cQuadTree::inOrder(cQuadTreeNode *tNode)
{
	if(tNode != NULL)
	{
		inOrder(tNode->getBottomLeft());
		inOrder(tNode->getBottomRight());
		cout<<"Node Center:"<<tNode->m_vCenter;
		if(tNode == m_Root)
			cout <<"root"<<endl;
		for(int i = 0; i < (int) tNode->m_iBoxes.size(); i++)
		{
			cout<<tNode->m_iBoxes[i]<<", ";
		}
		cout<<endl;
		inOrder(tNode->getTopRight());
		inOrder(tNode->getTopLeft());
	}
}//end inOrder





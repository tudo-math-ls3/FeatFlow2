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

#if !defined _C2DTREE_H_
#define _C2DTREE_H_


#include <cstdlib>
#include <iostream>
#include<vector>
#include "cVector2f.h"
#define AXIS_X 0
#define AXIS_Y 1


using namespace std;

/* forward declaration */
template<class T>
class c2DTree;

/* tree node object for a binary tree */
template<class T>
class c2DTreeNode
{

public:
	T   data;
	int m_iAxis;

	/* contructor */
	c2DTreeNode(const T &item, int axis, cVector2f vKey, c2DTreeNode<T> *lNode = NULL, c2DTreeNode<T> *rNode = NULL);

	/* destructor */
	~c2DTreeNode();
	
	/* get the child nodes */
	c2DTreeNode<T>* getLeft(){return m_nLeft;};
	c2DTreeNode<T>* getRight(){return m_nRight;};
	cVector2f		m_vKey;

	/* the main BinSTree class is a friend */
	friend class c2DTree<T>;

private:
	c2DTreeNode<T> *m_nLeft;
	c2DTreeNode<T> *m_nRight;

};

template<class T>
c2DTreeNode<T>::c2DTreeNode(const T &item, int axis, cVector2f vKey, c2DTreeNode<T> *lNode, c2DTreeNode<T> *rNode): data(item),
m_iAxis(axis), m_nLeft(lNode), m_nRight(rNode), m_vKey(vKey)
{}

template<class T>
c2DTreeNode<T>::~c2DTreeNode()
{
	/* deallocate left */
	if(m_nLeft != NULL)
	{
		delete m_nLeft;
		m_nLeft = NULL;
	}
	/* deallocate right */
	if(m_nRight != NULL)
	{
		delete m_nRight;
		m_nRight = NULL;
	}
}

/* Binary search tree definition */
template<class T>
class c2DTree
{
public:
	
	c2DTree(){m_Root = NULL; m_Current = NULL;};
	c2DTree(c2DTreeNode<T> *root);
	c2DTree(const c2DTree<T> &tree);
	~c2DTree(){};
	
	/* member functions */
	c2DTreeNode<T>* GetTreeNode(T item, int axis, cVector2f vKey, c2DTreeNode<T> *lChild = NULL, c2DTreeNode<T> *rChild = NULL);
	void inOrder(c2DTreeNode<T> *tNode);
	void inOrder(c2DTreeNode<T> *tNode, T tItem, vector<T> &vVals);
	void destroyTree();
	void deleteSubTree(c2DTreeNode<T> *tNode);
	c2DTreeNode<T>* getRoot(){return m_Root;};
	c2DTreeNode<T>* getCurrent(){cout<<"Wert: " <<m_Current->data<<endl;return m_Current;};
	void insert(const T &tItem, const cVector2f &vKey);
	c2DTreeNode<T>* findNode(T &tItem, c2DTreeNode<T> *&parent);
	int find(T &tItem);
	
	

protected:
	/* number of elements in the tree */
	int m_iSize;
	c2DTreeNode<T> *m_Root;
	c2DTreeNode<T> *m_Current;

};

/* class constructor: the parameter is a pointer to the root of the new tree */
template<class T>
c2DTree<T>::c2DTree(c2DTreeNode<T> *root)
{

	m_Root		= root;
	m_Current	= root;
	m_iSize		= 0;

}//end constructor

template<class T>
c2DTreeNode<T>* c2DTree<T>::GetTreeNode(T item, int axis, cVector2f vKey, c2DTreeNode<T> *lChild, c2DTreeNode<T> *rChild)
{
	c2DTreeNode<T> *p;

	p = new c2DTreeNode<T>(item ,axis, vKey, lChild, rChild);

	/* memory used up ? */
	if(p == NULL)
	{
		cerr<<"Memory allocation failure"<<endl;
		exit(1);
	}

	return p;
}//end GetTreeNode

template<class T>
void c2DTree<T>::insert(const T &tItem, const cVector2f &vKey)
{

	/* store a pointer to the current node and one to the parent */
	c2DTreeNode<T> *t = m_Root; 

	c2DTreeNode<T> *parent = NULL;

	c2DTreeNode<T> *nNewNode;
	
	int pAxis;
	double dKey;
	double dData;

	while(t != NULL)
	{
		parent = t;
		
		pAxis = t->m_iAxis;
		if(pAxis == AXIS_X)
		{
			dKey  = t->m_vKey.x;
			dData = vKey.x;
		}
		else
		{
			dKey  = t->m_vKey.y;
			dData = vKey.y;
		}
		if(dData < dKey)
			t = t->getLeft();
		else
			t = t->getRight();
	}//end while

	
	
	/* if parent is NULL, insert as root */
	if(parent == NULL)
	{
		/* create new node */
		nNewNode = GetTreeNode(tItem, AXIS_X, vKey, NULL, NULL);
		m_Root = nNewNode;
	}
	else if(dData < dKey)
	{
		/* calculate child axis value */
		int cAxis = (pAxis + 1) % 2;
		
		/* create new node */
		nNewNode = GetTreeNode(tItem, cAxis, vKey, NULL, NULL);
		
		/* assign new child */
		parent->m_nLeft  = nNewNode;
	}
	else
	{

		/* calculate child axis value */
		int cAxis = (pAxis + 1) % 2;

		/* create new node */
		nNewNode = GetTreeNode(tItem, cAxis, vKey, NULL, NULL);
		
		/* assign new child */
		parent->m_nRight = nNewNode;
	}

	/* assign new node as current */
	m_Current = nNewNode;

	m_iSize++;
	
}//end insert

template<class T>
void c2DTree<T>::inOrder(c2DTreeNode<T> *tNode)
{
	if(tNode != NULL)
	{
		inOrder(tNode->getLeft());
		cout<<"Node: "<<endl;
		if(tNode == m_Root)
			cout<<"root"<<endl;
		cout<<tNode->m_vKey;
		if(tNode->m_iAxis == AXIS_X)
			cout<<"X"<<endl;
		else
			cout<<"Y"<<endl;
		inOrder(tNode->getRight());
	}
}//end inOrder

template<class T>
void c2DTree<T>::destroyTree()
{
	if(m_Root)
	{
		deleteSubTree(m_Root);
		m_Root = NULL;
		m_Current = NULL;
	}

}//end destroyTree

template<class T>
void c2DTree<T>::deleteSubTree(c2DTreeNode<T> *tNode)
{
	if(tNode != NULL)
	{
		deleteSubTree(tNode->getLeft());
		tNode->m_nLeft = NULL;
        deleteSubTree(tNode->getRight());
		tNode->m_nRight = NULL;
		delete tNode;
		tNode = NULL;
	}
}//end deleteSubTree

template<class T>
int c2DTree<T>::find(T &tItem)
{
	c2DTreeNode<T> *parent;

	/* search tree for item */
	m_Current = findNode(tItem, parent);

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

}//end find

template<class T>
c2DTreeNode<T>* c2DTree<T>::findNode(T &tItem, c2DTreeNode<T> *&parent)
{
	/* start at root */
	c2DTreeNode<T> *t = m_Root;

	/* parent of root is NULL */
	parent = NULL;

	/* while there is a valid node pointer */
	while(t != NULL)
	{
		/* stop if item is found */
		if(tItem == t->data)
			break;
		else
		{
			/* update parent */
			parent = t;
			/* check where to go */
			if(tItem < t->data)
				t = t->m_nLeft;
			else
				t = t->m_nRight;
		}//end else
	}//end while

	/* return node */
    return t;

}//end findNode


template<class T>
void c2DTree<T>::inOrder(c2DTreeNode<T> *tNode, T tItem, vector<T> &vVals)
{
	if((tNode != NULL))
	{
		inOrder(tNode->getLeft(), tItem, vVals);
		if(tNode->data <= tItem)
		{
			vVals.push_back(tNode->data);
			inOrder(tNode->getRight(), tItem, vVals);
		}
	}
}//end inOrder


#endif

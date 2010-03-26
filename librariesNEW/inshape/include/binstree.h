/*
*
*
*
*
*
*/

#if !defined _CBINSTREE_H_
#define _CBINSTREE_H_

#include <cstdlib>
#include <iostream>
#include<vector>





template<class T>
class NodeDist
{
public:
	/* constructor */
	NodeDist(double adouble, T* abox);
	
	/* member functions */
	inline bool operator<(NodeDist<T> rhs){return d < rhs.d;};
		
	/* member variables */
	double d;
	T*  box;
	
};

template<class T>
NodeDist<T>::NodeDist(double adouble, T* abox) : d(adouble), box(abox)
{
	
}

struct NodeBound
{
	NodeBound() { tNode = NULL; dLBound = 0.0; };
	NodeBound(CBCNode *Node, double lBound)
	{
		tNode   = Node;
		dLBound = lBound;
	};
	CBCNode *tNode;
	double   dLBound;
	
	inline bool operator <(NodeBound &rhs){return dLBound < rhs.dLBound;};
	
};


/* namespace for std lib */
using namespace std;

/* forward declaration */
template<class T>
class cBinSTree;

/* tree node object for a binary tree */
template<class T>
class cTreeNode
{

public:
	T data;

	/* contructor */
	cTreeNode(const T &item, cTreeNode<T> *lNode = NULL, cTreeNode<T> *rNode = NULL);

	/* destructor */
	~cTreeNode();
	
	/* get the child nodes */
	cTreeNode<T>* getLeft(){return m_nLeft;};
	cTreeNode<T>* getRight(){return m_nRight;};

	/* the main BinSTree class is a friend */
	friend class cBinSTree<T>;

private:
	cTreeNode<T> *m_nLeft;
	cTreeNode<T> *m_nRight;

};

template<class T>
cTreeNode<T>::cTreeNode(const T &item, cTreeNode<T> *lNode, cTreeNode<T> *rNode): data(item), m_nLeft(lNode), m_nRight(rNode)
{}

template<class T>
cTreeNode<T>::~cTreeNode()
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
class cBinSTree
{
public:
	
	cBinSTree(){m_Root = NULL;};
	cBinSTree(cTreeNode<T> *root);
	cBinSTree(const cBinSTree<T> &tree);
	~cBinSTree(){};
	
	/* member functions */
	cTreeNode<T>* GetTreeNode(T item, cTreeNode<T> *lChild = NULL, cTreeNode<T> *rChild = NULL);
	
	void inOrder(cTreeNode<T> *tNode);
	void inOrder(cTreeNode<T> *tNode, T tItem, vector<T> &vVals);
	void destroyTree();
	void deleteSubTree(cTreeNode<T> *tNode);
	cTreeNode<T>* getRoot(){return m_Root;};
	cTreeNode<T>* getCurrent(){cout<<"Wert: " <<m_Current->data<<endl;return m_Current;};
	void insert(T &tItem);
	cTreeNode<T>* findNode(T &tItem, cTreeNode<T> *&parent);
	int find(T &tItem);
	void getLowest();
	

protected:
	/* number of elements in the tree */
	int m_iSize;
	cTreeNode<T> *m_Root;
	cTreeNode<T> *m_Current;

};

/* class constructor: the parameter is a pointer to the root of the new tree */
template<class T>
cBinSTree<T>::cBinSTree(cTreeNode<T> *root)
{

	m_Root		= root;
	m_Current	= root;
	m_iSize		= 0;

}//end constructor

template<class T>
cTreeNode<T>* cBinSTree<T>::GetTreeNode(T item, cTreeNode<T> *lChild, cTreeNode<T> *rChild)
{
	cTreeNode<T> *p;

	p = new cTreeNode<T>(item, lChild, rChild);

	/* memory used up ? */
	if(p == NULL)
	{
		cerr<<"Memory allocation failure"<<endl;
		exit(1);
	}

	return p;
}//end GetTreeNode

template<class T>
void cBinSTree<T>::insert(T &tItem)
{

	/* store a pointer to the current node and one to the parent */
	cTreeNode<T> *t = m_Root; 

	cTreeNode<T> *parent = NULL;

	cTreeNode<T> *nNewNode;

	while(t != NULL)
	{
		parent = t;
		if(tItem < t->data)
			t = t->getLeft();
		else
			t = t->getRight();
	}//end while

	/* create new node */
	nNewNode = GetTreeNode(tItem, NULL, NULL);

	/* if parent is NULL, insert as root */
	if(parent == NULL)
		m_Root = nNewNode;
	else if(tItem < parent->data)
		parent->m_nLeft  = nNewNode;
	else
		parent->m_nRight = nNewNode;

	/* assign new node as current */
	m_Current = nNewNode;

	m_iSize++;
	
}//end insert

template<class T>
void cBinSTree<T>::inOrder(cTreeNode<T> *tNode)
{
	if(tNode != NULL)
	{
		inOrder(tNode->getLeft());
		cout<<tNode->data<<endl;
		inOrder(tNode->getRight());
	}
}//end inOrder

template<class T>
void cBinSTree<T>::destroyTree()
{
	if(m_Root)
	{
		deleteSubTree(m_Root);
		m_Root = NULL;
		m_Current = NULL;
	}

}//end destroyTree

template<class T>
void cBinSTree<T>::deleteSubTree(cTreeNode<T> *tNode)
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
int cBinSTree<T>::find(T &tItem)
{
	cTreeNode<T> *parent;

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
cTreeNode<T>* cBinSTree<T>::findNode(T &tItem, cTreeNode<T> *&parent)
{
	/* start at root */
	cTreeNode<T> *t = m_Root;

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
void cBinSTree<T>::getLowest()
{
	m_Current = m_Root;
	
	while(m_Current->m_nLeft != NULL)
	{
		m_Current = m_Current->m_nLeft;
	}
}

template<class T>
void cBinSTree<T>::inOrder(cTreeNode<T> *tNode, T tItem, vector<T> &vVals)
{
	if((tNode != NULL))
	{
		inOrder(tNode->getLeft(), tItem, vVals);
		if(tNode->data < tItem)
		{
			vVals.push_back(tNode->data);
			inOrder(tNode->getRight(), tItem, vVals);
		}
	}
}//end inOrder

#endif

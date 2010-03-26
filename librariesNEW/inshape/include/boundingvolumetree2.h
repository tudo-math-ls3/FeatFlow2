#pragma once


///==============================================
///	               INCLUDES
///==============================================

#include "vector2.h"
#include <vector>
#include <queue>

using namespace std;


///==============================================
///	          Forward Declarations
///==============================================

template<class BV, class T, class SD>
class CBoundingVolumeTree2;

class CSBBNode;
class CSBBTree;

///==============================================
///	      CLASS CBoundingVolumeNode2
///==============================================

template<class BV, class T>
class CBoundingVolumeNode2
{
public:

	//standard constructor
	CBoundingVolumeNode2();

	//copy constructor
	CBoundingVolumeNode2(CSBBNode *pNode);

	//deconstructor
	~CBoundingVolumeNode2();

	//member functions

	//Is this node a leaf?
	inline bool IsLeaf() {return m_Children[0] == NULL;};

	//returns the center of the node

	//returns the lower bound for the branch and bound calculation
	inline T GetLowerBound(const CVector2<T> &vQuery) {return m_BV.MinDistance(vQuery);};

	//returns the upper bound for the branch and bound calculation
	inline T GetUpperBound(const CVector2<T> &vQuery) {return m_BV.MaxDistance(vQuery);};

	//returns the squared lower bound for the branch and bound calculation
	inline T GetLowerBoundSqr(const CVector2<T> &vQuery) {return m_BV.MinDistanceSqr(vQuery);};

	//returns the squared upper bound for the branch and bound calculation
	inline T GetUpperBoundSqr(const CVector2<T> &vQuery) {return m_BV.MaxDistanceSqr(vQuery)};

	//translate the node by the vector vTranslation
	inline void Translate(const CVector2<T> &vTranslation) {m_BV.Translate(vTranslation;};

	//rotate the node by the angle nAngle
	inline void Rotate(T nAngle) {m_BV.Rotate(nAngle);};

	//set the object id
	inline void SetID(int iID) {m_iID = iID;};

	//return the object id
	inline int GetID() {return m_iID;};

	//set the data index
	inline void SetData(int iData) {m_iData = iData;};

	//return the data index
	inline int GetData() {return m_iData;};

	//returns the center of the bounding volume
	inline CVector2<T> GetCenter() {return m_BV.GetCenter();};

	//public member variables

	//pointer to the children of this node
	CBoundingVolumeNode2<BV,T> *m_Children[2];

	//friends
	friend class CBoundingVolumeTree2<BV,T,class>;

private:

	//the node's bounding volume
	BV m_BV;

	//the data index
	int m_iData;

	//the object reference
	int m_iID;

};


///==============================================
///	      CLASS CBoundingVolumeTree2
///==============================================


template<class BV, class T, class SD>
class CBoundingVolumeTree2
{
public:

	//constructor
	CBoundingVolumeTree2();

	//member functions

	//initializes the tree
	void InitTree(vector<int> *&vBucket, CVector2<T> vSamples[], CSBBTree **pTrees, int iNumSamples, int iNumChildren, int iID);

	//deallocates all nodes etc...
	void DestroyTree();

	//updates the nodes after a translation
	void UpdateTree();

	//delete the subtree with root pNode
	void DeleteSubTree(CBoundingVolumeNode2<BV,T> *pNode);

	//returns the depth of the node pNode
	int GetDepth(CBoundingVolumeNode2<BV,T> *pNode, int iDepth);

	//returns the root
	inline CBoundingVolumeNode2<BV,T>* GetChild(){return m_Children[i];};

	//returns the number of children
	inline int GetNumChildren()
	{
		return m_iNumChildren;
	}//end GetNumChildren

	

	//deconstructor
	virtual ~CBoundingVolumeTree2();

private:

	//private member functions

	//deletes the subtree with root pNode
	void DeleteSubTree(CBoundingVolumeNode2<BV,T> *pNode);

	//copies the subtree with node pNode
	void CopySubTree(CBoundingVolumeNode2<BV,T> *pNode);

	//private member variables
	CBoundingVolumeNode2<BV,T> **m_Children;

	//number of child nodes
	int m_iNumChildren;
	


};

#include "boundingvolumetree2.cpp"
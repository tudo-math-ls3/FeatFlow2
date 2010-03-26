#include "vector2.h"
#include "matrix2x2.h"
#include "dynamicarray.h"
#include "mathglobals.h"
#include "sbb.h"


using namespace std;

/* forward declaration */
class CBasicCurve;
class CBSpline;
class CApproxCurve;










class CSBBTree;

class CSBBNode
{
public:

	/* constructor */
	CSBBNode(CSBB *pBox);

	CSBBNode();

	/* deconstructor */
	~CSBBNode();

	/* member functions */
	void subdivideNode(const VECTOR2 &vPoint,const VECTOR2 &vP1, const VECTOR2 &vP2,  Real dValue);
	Real GetUpperBound(const VECTOR2 &vQuery);
	Real GetLowerBoundDebug(const VECTOR2 &vQuery, VECTOR2 &vFar);
	
	
	/* inline member functions */
	inline Real GetLowerBound(const VECTOR2 &vQuery)
	{
		Real d = 0.0;
		Real minDist = 1.7E+308;
		for(int j = 0; j < 2; j++)
		{
			VECTOR2 b = m_Box->m_Vertices[m_Box->m_CurvePoints[j]];
			b = VECTOR2::createVector(vQuery, b);
			d = b.mag();
			if(d < minDist)
			{
				minDist   = d;
			}//end if
		}//end for
		return minDist;	
	}//end GetLowerBound

	inline bool IsLeaf() {return (m_Children[0]==NULL); };
	void SubdivideNode(Real dValue);


	/* public member variables */
	CSBBNode *m_Children[2];
	CSBB *m_Box;

	/* the tree class is a friend */
	friend class CSBBTree;

};

class CSBBTree
{

public:

	/* constructor */
	CSBBTree();
	CSBBTree(CSBBNode *tNode);
	
	/* deconstructor */
	~CSBBTree();

	/* member functions */
	void deleteSubTree(CSBBNode *pNode);
	void subdivideNode(CSBBNode *pNode, const VECTOR2 &vPoint,const VECTOR2 &vP1, const VECTOR2 &vP2,  Real dValue);
	void translateTree(CSBBNode *pNode, const VECTOR2 &vTrans);
	void destroyTree();

	/* inline member functions */
	inline CSBBNode* getRoot(){return m_Root;};
	void SubdivideNode(CSBBNode *pNode, const VECTOR2 &vPoint,  Real dValue);

private:

	CSBBNode *m_Root;

};
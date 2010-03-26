

///==============================================
///	      CBoundingVolumeNode2 Functions
///==============================================

template<class BV, class T>
CBoundingVolumeNode2<BV, T>::CBoundingVolumeNode2()
{

	m_Children[0] = 0;
	m_Children[1] = 0;

}//end constructor

template<class BV, class T>
CBoundingVolumeNode2<BV, T>::CBoundingVolumeNode2(CSBBNode *pNode) : m_Children[0](0), m_Children[1](0)
{

}//end constructor

template<class BV, class T>
CBoundingVolumeNode2<BV, T>::~CBoundingVolumeNode2()
{

	if(m_Children[0])
	{
		delete m_Children[0];
		m_Children[0] = NULL;
	}//end if

	if(m_Children[1])
	{
		delete m_Children[1];
		m_Children[1] = NULL;
	}//end if

}//end deconstructor

template<class BV, class T>
CBoundingVolumeNode2<BV, T>::CBoundingVolumeNode2()
{
	for(int i = 0; i < 2; i++)
	{
		m_Children[i] = NULL;
	}//end for
}//end constructor

///==============================================
///	      CBoundingVolumeTree2 Functions
///==============================================

template<class BV, class T, class SD>
CBoundingVolumeTree2<BV, T, SD>::CBoundingVolumeTree2() : m_iNumChildren(0)
{

}//end constructor

template<class BV, class T, class SD>
CBoundingVolumeTree2<BV, T, SD>::~CBoundingVolumeTree2()
{

}//end deconstructor

template<class BV, class T, class SD>
void CBoundingVolumeTree2<BV, T, SD>::DeleteSubTree(CBoundingVolumeNode2<BV,T> *pNode)
{

	if(pNode != NULL)
	{
		DeleteSubTree( pNode->m_Children[0] );
		pNode->m_Children[0] = NULL;
		DeleteSubTree( pNode->m_Children[1] );
		pNode->m_Children[1] = NULL;
		delete pNode;
		pNode = NULL;
	}//end if

}//end

template<class BV, class T, class SD>
void CBoundingVolumeTree2<BV, T, SD>::DestroyTree()
{

	for(int i = 0; i < m_iNumChildren; i++)
	{
		DeleteSubTree(m_Children[i]);
	}//end for


}//end destroy tree


template<class BV, class T, class SD>
void CBoundingVolumeTree2<BV, T, SD>::InitTree(vector<int> *&vBucket, CVector2<T> vSamples[], CSBBTree **pTrees, int iNumSamples, int iNumChildren, int iID)
{

	//set the number of children
	m_iNumChildren = iNumChildren;

	////allocate memory for the children of the root
	m_Children = new CBoundingVolumeNode2<BV,T>*[m_iNumChildren];

	////a counter that holds the storage information
	int iCount = 0;

	int i;

	//queue data structures for the Top-Down tree construction

	//holds the AABBTree's nodes
	queue<CSBBNode*> qBFS;

	////the circle tree's nodes
	queue<CBCNode*>  qNodes;

	/* create the top level nodes in the hierarchy */
	for(i = 0; i < m_iNumChildren; i++)
	{
		//insert the AABBTrees nodes into the queue
		//and construct circle tree nodes from them
		CSBBNode *SBBNode = pTrees[i]->getRoot();
		m_Children[i] = new CBoundingVolumeNode2<BV,T>(SBBNode);
		//set the curve ID
		m_Children[i]->SetID(iID);
		qBFS.push(SBBNode);
		qNodes.push(m_Children[i]);
	}//end for
	
	/* Top-Down build of the tree */
	while(!qBFS.empty())
	{
		CSBBNode *SBBNode = qBFS.front();
		CBoundingVolumeNode2<BV,T> *BCNode  = qNodes.front();

		/* create a lower bound for the node by a heuristic */
	//	BCNode->GenerateLowerBound(vSamples, iNumSamples);

		qBFS.pop();
		qNodes.pop();
	//	/* an inner node is subdivided */
		if(!SBBNode->IsLeaf())
		{			
	//		qBFS.push(SBBNode->m_Children[0]);
	//		qBFS.push(SBBNode->m_Children[1]);
	//		BCNode->m_Children[0] = new CBCNode(SBBNode->m_Children[0]);
	//		BCNode->m_Children[0]->SetID(iID);
	//		BCNode->m_Children[1] = new CBCNode(SBBNode->m_Children[1]);
	//		BCNode->m_Children[1]->SetID(iID);
	//		qNodes.push(BCNode->m_Children[0]);
	//		qNodes.push(BCNode->m_Children[1]);
		}//end if
		/* for a leaf node the data is assigned */
		else
		{
	//	//	AssignPointsAndParameters(pCurve,BCNode,nIndex);
	//		AssignPoints(BCNode, vSamples, vBucket[iCount], iNumSamples, iCount);
	//		iCount++;
	//	//	m_Leaves.push_back(BCNode);
		}//end else
	}//end while

}//end InitTree



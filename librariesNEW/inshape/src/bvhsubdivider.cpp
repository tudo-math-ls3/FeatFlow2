#include "bvhsubdivider.h"
#include "cbvhnode.h"


void CBVHSubdivider::SubdivideNode(CBVHNode<CAABB3f,float> *pNode, const vector<CTriangle3f*> *vPrimitives) const
{

	/* get the nodes bounding box */
	const CAABB3f bAABB3 = pNode->GetBV();

	///* get the longest axis the bounding volume will be split there */
	int iAxis    = bAABB3.LongestAxis();

	///* get the center of the bounding volume */
	VECTOR3 vCenter = bAABB3.GetCenter();

	const vector<CTriangle3f*> &vTriangles = *vPrimitives;

	
	//* assign the parents items to the children */
	vector<CTriangle3f*> *v0Triangles = new vector<CTriangle3f*>();
	vector<CTriangle3f*> *v1Triangles = new vector<CTriangle3f*>();

	//* split the items into two buckets relative to the split axis */

	for(int i = 0; i < vTriangles.size(); i++)
	{
		CTriangle3f *Tri = vTriangles[i];
		/* value at that the bounding volume is split along the split axis */
		if(Tri->GetCenter().m_dCoords[iAxis] < vCenter.m_dCoords[iAxis])
		{
			v0Triangles->push_back(Tri);
		}
		else
		{
			v1Triangles->push_back(Tri);
		}
	}//end for

	//* create the new children and assign the items to the children */
	pNode->m_Children[0] = new CBVHNode<CAABB3f, Real>();
	pNode->m_Children[0]->SetTriangles(v0Triangles);
	pNode->m_Children[0]->GetBV().Init(*v0Triangles);
	pNode->m_Children[1] = new CBVHNode<CAABB3f, Real>();
	pNode->m_Children[1]->SetTriangles(v1Triangles);
	pNode->m_Children[1]->GetBV().Init(*v1Triangles);

}//end SubdivideNode

void CBVHSubdivider::SubdivideNode(CBVHNode<CAABB3f,float> *pNode, const vector<CTriangle3f*> *vPrimitives
								   ,vector<VECTOR3*> &vVerts) const
{

	/* get the nodes bounding box */
	const CAABB3f bAABB3 = pNode->GetBV();

	///* get the longest axis the bounding volume will be split there */
	int iAxis    = bAABB3.LongestAxis();

	///* get the center of the bounding volume */
	VECTOR3 vCenter = bAABB3.GetCenter();

	const vector<CTriangle3f*> &vTriangles = *vPrimitives;

	
	//* assign the parents items to the children */
	vector<CTriangle3f*> *v0Triangles = new vector<CTriangle3f*>();
	vector<CTriangle3f*> *v1Triangles = new vector<CTriangle3f*>();

	//* split the items into two buckets relative to the split axis */

	for(int i = 0; i < vTriangles.size(); i++)
	{
		CTriangle3f *Tri = vTriangles[i];
		/* value at that the bounding volume is split along the split axis */
		if(Tri->GetCenter().m_dCoords[iAxis] < vCenter.m_dCoords[iAxis])
		{
			v0Triangles->push_back(Tri);
		}
		else
		{
			v1Triangles->push_back(Tri);
		}
	}//end for

	//* create the new children and assign the items to the children */
	pNode->m_Children[0] = new CBVHNode<CAABB3f, Real>();
	pNode->m_Children[0]->SetTriangles(v0Triangles);
	pNode->m_Children[0]->GetBV().Init(*v0Triangles);
	pNode->m_Children[0]->InnerVertices(vVerts);
	pNode->m_Children[1] = new CBVHNode<CAABB3f, Real>();
	pNode->m_Children[1]->SetTriangles(v1Triangles);
	pNode->m_Children[1]->GetBV().Init(*v1Triangles);
	pNode->m_Children[1]->InnerVertices(vVerts);

}//end SubdivideNode


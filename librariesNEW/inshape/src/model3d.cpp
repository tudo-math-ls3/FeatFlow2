#include "model3d.h"

CModel3D::CModel3D(void) : m_bHasUV(false)
{
	
}

CModel3D::~CModel3D(void)
{

	for(int i = 0; i < m_trTriangles.size(); i++)
	{

		delete m_trTriangles[i];
		m_trTriangles[i] = NULL;

	}//end for

}

void CModel3D::SetVertices(const VertArray &pVertices)
{

	int iSize = (int)pVertices.size();
	VertArray::iterator it;
	m_pVertices.Resize(iSize);
	int i;
	for(i=0; i < iSize; i++)
	{
		m_pVertices[i] = pVertices[i];
	}//end for

	m_bxBox.InitBox(m_pVertices);

}

void CModel3D::SetFaces(const FaceArray &pFaces)
{
	m_pFaces = pFaces;
}

void CModel3D::SetNormals(const Vec3Array &pNormals)
{
	m_pVertexNormals = pNormals;
}

void CModel3D::SetTCoords(const Vec2Array &pTCoords)
{
	m_pTexCoords = pTCoords;
	m_bHasUV = true;
}

void CModel3D::BuildTriangles()
{

	

	for(int i = 0; i < (int)m_pFaces.size(); i++)
	{
		m_trTriangles.push_back(MakeTriangle(i));
	}//end for

}//end BuildTriangles



#include <sbb.h>
#include <bspline.h>

/* class CSBB */

//standard constructor
//bla
CSBB::CSBB()
{
}

/////////////////////////////////////////////////////////
//                  Constructor                        //
// Input: Two points on the curve and their parameter  //
//        values.                                      //
/////////////////////////////////////////////////////////



CSBB::CSBB(VECTOR2 vA, VECTOR2 vB, Real fA, Real fB, CBSpline *pCurve) : m_fP1(fA), m_fP2(fB)
{

	m_pCurve = pCurve;
	
	if((vA.y >= vB.y) && (vA.x >= vB.x))
	{
		m_CurvePoints[0] = 2;
		m_CurvePoints[1] = 0;
		m_Vertices[0]    = vB;
		m_Vertices[1].x  = vA.x;
		m_Vertices[1].y  = vB.y;
		m_Vertices[2]    = vA;
		m_Vertices[3].x  = vB.x;
		m_Vertices[3].y  = vA.y;
	}
	else if((vB.y >= vA.y) && (vB.x >= vA.x))
	{
		m_CurvePoints[0] = 2;
		m_CurvePoints[1] = 0;
		m_Vertices[0]    = vA;
		m_Vertices[1].x  = vB.x;
		m_Vertices[1].y  = vA.y;
		m_Vertices[2]    = vB;
		m_Vertices[3].x  = vA.x;
		m_Vertices[3].y  = vB.y;
	}

	else if((vA.y >= vB.y) && (vB.x >= vA.x))
	{
		m_CurvePoints[0] = 3;
		m_CurvePoints[1] = 1;
		m_Vertices[0].x  = vA.x;
		m_Vertices[0].y  = vB.y;
		m_Vertices[1]    = vB;
		m_Vertices[2].x  = vB.x; 
		m_Vertices[2].y  = vA.y;
		m_Vertices[3]    = vA;
	}
	else
	{
		m_CurvePoints[0] = 3;
		m_CurvePoints[1] = 1;
		m_Vertices[0].x  = vB.x;
		m_Vertices[0].y  = vA.y;
		m_Vertices[1]    = vA;
		m_Vertices[2].x  = vA.x; 
		m_Vertices[2].y  = vB.y;
		m_Vertices[3]    = vB;
		
	}

}//end constructor

CSBB::CSBB(CSBB *pBox)
{

	m_pCurve = pBox->m_pCurve;

	for(int i = 0; i < 4; i++)
	{
		m_Vertices[i] = pBox->m_Vertices[i];
	}//end for

	for(int j = 0; j < 2; j++)
	{
		m_CurvePoints[j] = pBox->m_CurvePoints[j];
	}//end for

	m_fP1 = pBox->m_fP1;
	m_fP2 = pBox->m_fP2;
}//end constructor



//////////////////////////////////////////
// Test if two boxes overlap            
// Input : extreme coordinates of the box
// Output: true if overlap else false
//////////////////////////////////////////



ostream& operator<<(ostream& out, CSBB *b1)
{

	out <<"Begin box: "<<endl;
	for(int i = 0; i < 4; i++)
	{
		out << b1->m_Vertices[i];
	}//end for

	out <<b1->m_fP1<<endl;

	out <<b1->m_fP2<<endl;

	out<<"End box..."<<endl;

	return out;

}//end operator

/*!
    \fn CSBB::CSBB(std::vector<VECTOR2> vPoints, Real start, Real end)
 */
CSBB::CSBB(std::vector<VECTOR2> vPoints, Real start, Real end, CBSpline *pCurve)
{
	
	Real yB =  1.7E+308;
	Real yT = -1.7E+308;
	Real xL =  1.7E+308;
	Real xR = -1.7E+308;

	int num = (int)vPoints.size();

	for(int i = 0; i < num;i++)
	{
		if(yB > vPoints[i].y)
			yB = vPoints[i].y;
		if(yT < vPoints[i].y)
			yT = vPoints[i].y;
		if(xL > vPoints[i].x)
			xL = vPoints[i].x;
		if(xR < vPoints[i].x)
			xR = vPoints[i].x;
		
	}

	this->m_fP1 = start;	
	this->m_fP2 = end;
	
	m_Vertices[0].x = xL;
	m_Vertices[0].y = yB;
	m_Vertices[1].x = xR;
	m_Vertices[1].y = yB;
	m_Vertices[2].x = xR; 
	m_Vertices[2].y = yT;
	m_Vertices[3].x = xL;
	m_Vertices[3].y = yT;
	
	m_pCurve = pCurve;	
}

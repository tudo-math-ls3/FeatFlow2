//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : CBSpline.cpp
//  @ Date : 16.06.2006
//  @ Author : Raphael Mnster
//
//


#include <bspline.h>
#include <graham.h>
#include <exception.h>
#include <matrix2x2.h>
#include <fstream>
#include <queue>
#include <rootfinder.h>
#include <sbbtree.h>
#include <obb2.h>


/////////////////* GLOBALS *////////////////////////

int _ID;


using namespace std;


//////////////////////////////////////////////////////////////////
//                           Contructor                         //
// Input: degree, number of Points, pointer to the beginning of //
//        the array of knots, pointer to the start of the       //
//        control point array.								       //
//					                                              //
//////////////////////////////////////////////////////////////////

CBSpline::CBSpline(int deg, int numP, Real *kV, VECTOR2 *cP, int numWrapped) : m_ideg(deg),
						m_inumPoints(numP), m_sWrap(numWrapped)
{
	int i;
	m_inumKnots		= m_inumPoints + m_ideg;
	m_ControlPoints.Resize(m_inumPoints);
	m_KnotVector.Resize(m_inumKnots);
	m_vCog.x = 0;
	m_vCog.y = 0;
	for(i = 0; i < m_inumPoints; i++)
	{
		m_ControlPoints[i] = cP[i];
	}

	for(i = 0; i < m_inumPoints - numWrapped; i++)
	{
		m_vCog = m_vCog + m_ControlPoints[i];
	}

	m_vCog.x /= m_inumPoints - numWrapped;
	m_vCog.y /= m_inumPoints - numWrapped;

	for(i = 0; i < m_inumKnots; i++)
		m_KnotVector[i] = kV[i];

	m_dAngle = 0.0;
	m_bBox = NULL;



}

CBSpline::~CBSpline()
{
	// We print out this message to help visualize when destructors get called
	// in a program.
 	//printf("Destructor called\n");
	int i;

	for(i = 0; i <(int)m_Boxes.size(); i++)
	{
		delete m_Boxes[i];
	}

	if(m_bBox)
	{
		delete m_bBox;
		m_bBox = NULL;
	}

	if(m_Hierarchy[0])
	{
		destroyHierarchy();
	}

}//end deconstructor



//void CBSpline::drawCornerCutting(vector<VECTOR2>& points, int step)
//{
//	vector<VECTOR2>::iterator Iter;
//
//	if(step==1)
//		glColor3f(1.0,0.0,0.0);
//	else if(step==2)
//		glColor3f(1.0,1.0,0.0);
//	else if(step==3)
//		glColor3f(1.0,0.2,0.84);
//	
//	glBegin(GL_POINTS);
//		for ( Iter = points.begin( ); Iter != points.end( ); Iter++ )
//		{
//			VECTOR2 point = *Iter;
//			glVertex3f(point.x, point.y, 0);
//		}
//	glEnd();
//
//	if(points.size() > 1)
//	{
//		glBegin(GL_LINE_STRIP);
//			for ( Iter = points.begin( ); Iter != points.end( ); Iter++ )
//			{
//				VECTOR2 point = *Iter;
//				glVertex3f(point.x, point.y, 0);
//			}
//		glEnd();
//	}
//}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

VECTOR2 CBSpline::deBoor1(Real u, bool draw)
{
	
	int k = findKnotSpan(u);
	if(k < 0)
		exit(0);
	
	int p = this->m_ideg;
	Real saved;
	Real* left  = new Real[p+1];
	Real* right = new Real[p+1];
	Real* N		= new Real[p+1];

	N[0] = 1.0;

	for(int j = 1; j <= p; j++)
	{
		left[j]  = u - m_KnotVector[k+1-j];
		right[j] = m_KnotVector[k+j] - u;
		saved = 0.0;
		for(int r = 0; r < j; r++)
		{
			Real temp = N[r]/(right[r+1] + left[j-r]);
			N[r]	   = saved + right[r+1] * temp;
			saved	   = left[j-r] * temp;
			//printf("Basis function %f, Index %d, saved = %f\n",N[r],j,saved);
		}
		
		N[j] = saved;
	}

	
	for(int i = 0; i <= p; i++)
	{
		printf("real value %f\n", N[i]);
		
	}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		

	VECTOR2 res(0,0);


	for(int i = 0; i <= p; i++)
	{
		res.x += m_ControlPoints[k-p+i].x * N[i];
		res.y += m_ControlPoints[k-p+i].y * N[i];
	}


	delete left;  
	delete right; 
	delete N;

	return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////
//                                                  //
// This function  evaluates the curve at the        //
// parameter u. The following equation is computed  //
// C(u) = SUM_(0,N) (N_ip(u) * P_i)                      //
//                                                  //
/////////////////////////////////////////////////////



VECTOR2 CBSpline::CoxDeBoor(Real u)
{
	int k = findKnotSpan(u);
	if(k < 0)
	{
		printf( "WARNING INVALID KNOT SPAN FOUND!\n");
		return VECTOR2(0,0);
	}
	
	int p = this->m_ideg;
	Real saved;
	Real* left  = new Real[p+1];
	Real* right = new Real[p+1];
	Real* N		= new Real[p+1];

	N[0] = 1.0;

	for(int j = 1; j <= p; j++)
	{
		left[j]  = u - m_KnotVector[k+1-j];
		right[j] = m_KnotVector[k+j] - u;
		saved = 0.0;
		for(int r = 0; r < j; r++)
		{
			Real temp = N[r]/(right[r+1] + left[j-r]);
			N[r]	   = saved + right[r+1] * temp;
			saved	   = left[j-r] * temp;
			//printf("Basis function %f, Index %d, saved = %f\n",N[r],j,saved);
		}
		
		N[j] = saved;
	}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		


	VECTOR2 res(0,0);


	for(int i = 0; i <= p; i++)
	{
		res.x += m_ControlPoints[k-p+i].x * N[i];
		res.y += m_ControlPoints[k-p+i].y * N[i];
	}


	delete[] left;  
	delete[] right; 
	delete[] N;

	return res;
}


///////////////////////////////////////////////////////////////////////////////////////////////
//                                              
// calculates all derivatives up to the n-th  of the b-spline   
// basis functions N_ip at parameter u.     
//                                              
//////////////////////////////////////////////////////////////////////////////////////////////
Real** CBSpline::derivativeBasis(Real u, int n)
{

	int p = m_ideg;

	int k = findKnotSpan(u);
	if(k < 0)
	{
		printf("k ungueltig\n");
		return new Real*[1];
	}

	if( n > p)
	{
		printf("p zu gross\n");
		exit(0);
	}

	int j,i,r;

	Real saved;
	Real* left  = new Real[p+1];
	Real* right = new Real[p+1];
	
	//stores basis functions and knot differences
	Real** ndu = new Real*[p+1];
	for(i = 0; i < p+1; i++)
		ndu[i] = new Real[p+1];

	//stores kth derivative of N_i-p+j_p
	Real** ders = new Real*[p+1];
	for(i = 0; i < p+1; i++)
		ders[i] = new Real[p+1];

	//stores a_kj and a_k-1j
	Real ** a = new Real*[2];
	for(i = 0; i < 2; i++)
		a[i] = new Real[p+1];

	
//---------------------------------------------------------------------------
	
	
	ndu[0][0] = 1.0;
	for(j = 1; j <= p; j++)
	{
		left[j]  = u - m_KnotVector[k+1-j];
		right[j] = m_KnotVector[k+j] - u;
		saved = 0.0;
		for(r = 0; r < j; r++)
		{
			ndu[j][r]  = right[r+1] + left[j-r];
			Real temp = ndu[r][j-1]/ndu[j][r];

			ndu[r][j]  = saved + right[r+1] * temp;
			saved	   = left[j-r] * temp;
			//printf("Basis function %f, Index %d, saved = %f\n",N[r],j,saved);
		}
		ndu[j][j] = saved;
	}

	
	//load basis functions
	for(j = 0;j<=p;j++)
		ders[0][j] = ndu[j][p];

	//main algorithm
	for(r=0; r<=p; r++)
	{
		int s1  = 0;
		int s2  = 1;
		int j1  = 0;
		int j2  = 0;
		a[0][0] = 1.0f;
		//loop to compute kth derivative
		for(i = 1; i <= n; i++)
		{
			Real d = 0.0;
			int ri = r - i;
			int pi = p - i;
			
			if( r >= i )
			{
				a[s2][0] = a[s1][0]/ndu[pi+1][ri];
				d = a[s2][0]*ndu[ri][pi];
			}

			if(ri >= -1)
				j1 = 1;
			else
				j1 = -ri;

			if(r-1 <= pi)
				j2 = i-1;
			else
				j2 = p - r;

			for(j = j1; j<= j2; j++)
			{
				a[s2][j] = (a[s1][j] - a[s1][j-1])/ndu[pi+1][ri+j];
				d += a[s2][j] * ndu[ri+j][pi];
			}//end for j
			if( r <= pi)
			{
				a[s2][i] = -a[s1][i-1]/ndu[pi+1][r];
				d		+= a[s2][i] * ndu[r][pi];
			}
			ders[i][r] = d;
			//switch rows
			j = s1; s1 = s2; s2 = j;
		}// end for i
	}// end for r

	//multiply by prefactors
	r = p;
	for(i = 1; i <= n; i++)
	{
		for(j = 0; j <=p ; j++)
		{
			ders[i][j] *= r;
		}
		r *= (p-i);
	}
//----------------------------------------------------------------------------------

	for(i = 0; i < p + 1; i++)
		delete[] ndu[i];

	for(i = 0; i < 2; i++)
		delete[] a[i];

	delete[] left;
	delete[] right;
	delete[] ndu;
	delete[] a; 
	return ders;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                               //
// calculates all the derivatives up to the n-th at parameter u   //
// of the b-spline curve and returns an array A that contains    //
// the i-th derivative at A[i].                                  //
//                                                               //
////////////////////////////////////////////////////////////////////////////////////////////////////

VECTOR2* CBSpline::curveDerivs(Real u, int n)
{
	int p = m_ideg;

	int k = findKnotSpan(u);
	if(k < 0)
	{
		EXCEPTION("Invalid knot span");
	}

	if( n > p)
	{
		printf("n zu gross\n");
		exit(0);
	}

		

	VECTOR2* CK = new VECTOR2[n+1];
	
	Real** nders = this->derivativeBasis(u,n);

	for(int i = 0; i<=n ; i++)
	{
		CK[i] = VECTOR2();
		for(int j = 0; j <= p; j++)
		{
			CK[i].x = CK[i].x + nders[i][j] * m_ControlPoints[k-p+j].x;
			CK[i].y = CK[i].y + nders[i][j] * m_ControlPoints[k-p+j].y;
		}
	}

	for(int i = 0; i < p + 1; i++)
		delete[] nders[i];
	delete[] nders;
	return CK;
}

//////////////////////////////////////////////////////////////
//                                                          //
// This function generates bounding boxes for the monotonic //
// segments of the curve.                                   //
//                                                          //
//////////////////////////////////////////////////////////////

void CBSpline::genBoxes()
{

	if(!m_Boxes.empty())
		m_Boxes.clear();
   
	/* these std::vectors store the curve points where the monotonic segments end and begin and the parameter values */
	std::vector<Real>vRoots;
	std::vector<VECTOR2>vPoints;

	/* brute force algorithm to find the monotonic segments */
	findMonotoneSegments(vPoints, vRoots);

	/* here a bounding box for the curve is generated */
	graham aGraham;
	aGraham.getPoints(vPoints, 1000);
	m_Hull = aGraham.compConvexHull();

	if(!m_bBox)
		m_bBox = new COBB(vPoints, m_Hull, this);
	
	/* create bounding boxes for the monotonic segments */
	int s = (int) vPoints.size();
	m_Hierarchy = new CSBBTree*[s];
	m_iBoxes    = s;
	
	_ID = 0;
	for(int i = 0; i < s; i++)
	{
		CSBB *pCSBB = new CSBB(vPoints[i], vPoints[(i+1)%s], vRoots[i], vRoots[(i+1)%s], this);
		pCSBB->m_ID = i;
        m_Boxes.push_back(pCSBB);
		m_Hierarchy[i] = new CSBBTree(new CSBBNode(pCSBB));
	}//end for

	cout <<"Number of parent nodes: "<<s <<endl;

	createHierarchy();

}//end genBoxes



///////////////////////////////////////////////////////////
//
//		a hierarchy of axis-aligned bounding boxes for
//		
//		monotonic curve segments is build here.
//
///////////////////////////////////////////////////////////

void CBSpline::createHierarchy()
{
	/* a queue to store the nodes in the hierarchy */	
	queue<CSBBNode*> qNodes;

	/* insert the parent nodes into the queue */
	for(int i = 0; i < m_iBoxes; i++)
	{
		CSBBNode *boxNode = m_Hierarchy[i]->getRoot();
		qNodes.push(boxNode);
	}//end for
	
	/* the curve segments will be subdivided at this parameter value */
	Real fPmid = 0.0;

	/* leaf count at target depth */
	int nTargetNodeCount = NodeCount(TARGETDEPTH, m_iBoxes);

	/* current node count */
	int nNodeCount = 0;

	/* while depth limit is not yet reached */
	while(!qNodes.empty() && nNodeCount < nTargetNodeCount)
	{
		/* get the front node */
		CSBBNode *bNode = qNodes.front();

		/* release it from the queue */
		qNodes.pop();

		/* get the nodes bounding box */
		CSBB* box = bNode->m_Box;

		/* calculation of the parameter value the segment is split at */
		/* the segments are subdivided in the middle, but special casing  */
		/* is needed when the curve is wrapped around */
		if(box->m_fP2 > box->m_fP1)
		{

			fPmid  = (box->m_fP2 + box->m_fP1) / 2;
				
		}//end if
		/* special case */
		else
		{
			
			Real length = getIntervalLength(box);
			if((length/2 + box->m_fP1) >= m_KnotVector[getNumPoints()])
			{
				fPmid = length/2 + box->m_fP1 - m_KnotVector[getNumPoints()] + m_KnotVector[getDeg()];
			}//end if
			else
			{
				fPmid = box->m_fP1 + length/2;
			}//end else
		}//end else
			VECTOR2 vMid = CoxDeBoor(fPmid);
			VECTOR2 vA   = CoxDeBoor(box->m_fP1);
			VECTOR2 vB   = CoxDeBoor(box->m_fP2);

			///* first new box */
			bNode->subdivideNode(vMid,vA, vB, fPmid);
			qNodes.push(bNode->m_Children[0]);
			qNodes.push(bNode->m_Children[1]);
			
			nNodeCount++;
						
	}//end while
}//end createHierarchy

//////////////////////////////////////////////
//                                          //
// Subdivides the bounding box by splitting //
// the parameter interval in half           //
//                                          //
//////////////////////////////////////////////                         

void CBSpline::subdivideBoxes()
{
	std::vector<CSBB*> boxes;
	int s = (int) m_Boxes.size();
	Real fPmid;
	for(int i = 0; i < s ; i++)
	{
		CSBB* box = m_Boxes[i];

		///* stop subdivision if interval too small */
		//if(getIntervalLength(box) <= 0.02)
		//{
		//	boxes.push_back(new CSBB(box->m_Vertices[box->m_CurvePoints[0]], box->m_Vertices[box->m_CurvePoints[1]], box->m_fP1, box->m_fP2));
		//	continue;
		//}//end if

		///* stop subdivision if area too small */
		//if( box->getArea() <= 1e-3 )
		//{
		//	boxes.push_back(new CSBB(box->m_Vertices[box->m_CurvePoints[0]], box->m_Vertices[box->m_CurvePoints[1]], box->m_fP1, box->m_fP2));
		//	continue;
		//}//end if

		if(box->m_fP2 > box->m_fP1)
		{

			fPmid  = (box->m_fP2 + box->m_fP1) / 2;
			
		}//end if
		else
		{
			Real length = getIntervalLength(box);
			if((length/2 + box->m_fP1) >= m_KnotVector[getNumPoints()])
			{
				 fPmid = length/2 + box->m_fP1 - m_KnotVector[getNumPoints()] + m_KnotVector[getDeg()];
			}//end if
			else
			{
				fPmid = box->m_fP1 + length/2;
			}//end else
		}//end else
			VECTOR2 vB = CoxDeBoor(fPmid);
			VECTOR2 vA = CoxDeBoor(box->m_fP1);

			/* first new box */
			boxes.push_back(new CSBB(vA, vB, box->m_fP1, fPmid, this));
			vA = vB;
			vB = CoxDeBoor(box->m_fP2);
			/* second new box */
			boxes.push_back(new CSBB(vA, vB, fPmid, box->m_fP2, this));
	}//end for
	
	for(int i = 0; i < (int)m_Boxes.size(); i++)
		delete m_Boxes[i];
	m_Boxes.clear();
	m_Boxes.assign(boxes.begin(), boxes.end());

}//end subdivideBoxes

////////////////////////////////////////////////////////
//                                                    //
// A functions that finds the monotonic segments of a // 
// b-spline curve.                                    //
// Input: A vector that stores the parameter values   //
//        at the start and end points of the curve.   // 
//        A vector that stores the curve points at    //
//        the beginning and end of a monotonic        //
//        segment.                                    //
////////////////////////////////////////////////////////

void CBSpline::findMonotoneSegments(std::vector<VECTOR2>& vPoints, std::vector<Real>& vRoots)
{
	
	VECTOR2* deriv;
	Real T  = m_KnotVector[m_ideg];
	bool signX;
	bool signY;
	CRootFinder rFA;

	Real dT = (m_KnotVector[m_inumPoints] - m_KnotVector[m_ideg]) / 2000;
	deriv = curveDerivs(T,2);
	if(deriv[1].y > 0)
		signX = true;
	else
		signX = false;

	if(deriv[1].x > 0)
		signY = true;
	else
		signY = false;

	if(deriv[1].x == 0.0)
	{
		vRoots.push_back(m_KnotVector[m_ideg]);
		vPoints.push_back(CoxDeBoor(m_KnotVector[m_ideg]));
	}

	
	if(deriv[1].y == 0.0)
	{
		vRoots.push_back(m_KnotVector[m_ideg]);
		vPoints.push_back(CoxDeBoor(m_KnotVector[m_ideg]));
	}

	delete[] deriv;
	int count = 0;
	for(int i = 0; i < 2000; i++)
	{
			
		deriv = curveDerivs(T,2);
		if((signX == false) && (deriv[1].y > 0))
		{
			signX		= true;
			Real p = rFA.bisectX(this, T-dT, T, 1e-3);
			vRoots.push_back(p);
			vPoints.push_back(CoxDeBoor(p));
			count++;
		}
		if((signX == true) &&(deriv[1].y < 0))
		{
			signX		= false;
			Real p = rFA.bisectX(this, T-dT, T, 1e-3);
			vRoots.push_back(p);
			vPoints.push_back(CoxDeBoor(p));
			count++;
		}

		if((signY == false) && (deriv[1].x > 0))
		{
			signY		= true;
			Real p = rFA.bisectY(this, T-dT, T, 1e-3);
			vRoots.push_back(p);
			vPoints.push_back(CoxDeBoor(p));
			count++;
		}
		if((signY == true) &&(deriv[1].x < 0))
		{
			signY		= false;
			Real p = rFA.bisectY(this, T-dT, T, 1e-3);
			vRoots.push_back(p);
			vPoints.push_back(CoxDeBoor(p));
			count++;
		}

		T+=dT;
		delete[] deriv;
	}

	
	
}

//////////////////////////////////////////////////////////////
//                                                          //
// A function that inserts a new knot a parameter u t times //
//                                                          //
//////////////////////////////////////////////////////////////

bool CBSpline::insertKnot(Real u, int t)
{

	int p = m_ideg;

	int mp = m_inumPoints + m_ideg + 1;

	int nq = m_inumPoints + t;

	int k = findKnotSpan(u);

	int s = 0;
	if( s = getMultiplicity(u) >= p)
		return false;
	
	//new number of knots
	m_inumKnots+=t;

	//allocate memory for new knot vector and control points
	Real* UQ = new Real[m_inumKnots];
	VECTOR2* QW = new VECTOR2[nq];
	VECTOR2* RW = new VECTOR2[p+1];
	

	//Load new knot vector
	int i,j,l;
	for(i = 0; i <= k ; i++)
		UQ[i] = m_KnotVector[i];
	for(i = 1; i <= t; i++)
		UQ[k+i] = u;
	for(i = k + 1; i < mp; i++)
		UQ[i+t] = m_KnotVector[i];

	//Save unaltered control points
	for(i = 0; i <=k-p; i++)
		QW[i] = m_ControlPoints[i];
	for(i = k-s; i< m_inumPoints; i++)
		QW[i+t] = m_ControlPoints[i];
	for(i = 0; i <= p-s; i++)
		RW[i] = m_ControlPoints[k-p+i];

	//main algorithm : inserts knot value u t times
	for(j = 1; j <= t; j++)
	{
		l = k-p+j;
		for(i = 0; i <= p-j-s; i++)
		{
			Real alpha = (u - m_KnotVector[l+i])/(m_KnotVector[i+k+1] - m_KnotVector[l+i]);
			RW[i].x		 = alpha * RW[i+1].x + (1.0 - alpha) * RW[i].x;
			RW[i].y		 = alpha * RW[i+1].y + (1.0 - alpha) * RW[i].y;
		}
		QW[l] = RW[0];
		QW[k+t-j-s] = RW[p-j-s];
	}

	//load remaining control points
	for(i = l+1; i <k-s; i++)
	{
		QW[i] = RW[i-l];
	}

	m_ControlPoints.Resize(QW,nq);
	m_KnotVector.Resize(UQ, m_inumKnots);
	delete[] RW;
	m_inumPoints +=t;
	//wrapPoints();
	return true;
}

void CBSpline::wrapPoints()
{

	int i;
	for(i = 1; i < m_ideg-1; i++)
	{
		m_ControlPoints[m_inumPoints-m_ideg+i] = m_ControlPoints[i];
		
	}
	for(i = 0; i < m_inumPoints; i++)
		printf("[%f,%f]\n",m_ControlPoints[i].x, m_ControlPoints[i].y);
}

/* recomputes the center of gravity in case the curve was moved or rotated */
void CBSpline::updateCog()
{
	m_vCog.x = 0;
	m_vCog.y = 0;
	for(int i = 0; i < m_inumPoints - m_sWrap; i++)
	{
		m_vCog = m_vCog + m_ControlPoints[i];
	}

	Real invNum =(Real) 1.0/(m_inumPoints - m_sWrap);
	m_vCog = m_vCog * invNum;
}//end updateCog


void CBSpline::genHullBox()
{
	if(m_Hull.empty())
		genHull();

    if(!m_bBox)
		m_bBox = new COBB(m_ControlPoints, m_Hull, this);
}


/* scales the curve by dScale */
void CBSpline::scale(Real dScale)
{

	VECTOR2 vOldCog = getCog();

	for(int i = 0; i < m_inumPoints; i++)
		m_ControlPoints[i] = m_ControlPoints[i] * dScale;
	
	updateCog();
    if(m_bBox)
	{
	  delete m_bBox;
	  m_bBox = NULL;
	}

	genBoxes();

	VECTOR2 vNewCog = getCog();

	translateCurve( (vOldCog - vNewCog) );

}//end scale

/////////////////////////////////////////////////
//
//	translate the curve along the vector vTrans
//
//
/////////////////////////////////////////////////

void CBSpline::translateCurve(VECTOR2 vTrans)
{
	int i;
	//translate control points
	for(i = 0; i < m_inumPoints; i++)
	{
		m_ControlPoints[i] = m_ControlPoints[i] + vTrans;
	}
	
	int m = (int)m_Boxes.size();
	//translate segment boxes
	for(i = 0; i < m; i++)
	{
		m_Boxes[i]->m_Vertices[3] = m_Boxes[i]->m_Vertices[3] + vTrans;
		m_Boxes[i]->m_Vertices[2] = m_Boxes[i]->m_Vertices[2] + vTrans;
		m_Boxes[i]->m_Vertices[0] = m_Boxes[i]->m_Vertices[0] + vTrans;
		m_Boxes[i]->m_Vertices[1] = m_Boxes[i]->m_Vertices[1] + vTrans;
	}
	//translate curve bounding box
	
	
	m_bBox->m_Vertices[0] = m_bBox->m_Vertices[0] + vTrans;
	m_bBox->m_Vertices[1] = m_bBox->m_Vertices[1] + vTrans;

	updateCog();
}//end translateCurve

void CBSpline::rotateBack(Real dAngle)
{
}

////////////////////////////////////////
//
//	rotates the curve dAngle degrees
//  around the center of gravity
//
////////////////////////////////////////
void CBSpline::rotateCurve(Real dAngle)
{
	
	VECTOR2 oldCog = m_vCog;
	VECTOR2 vTrans = m_vCog * -1.0;
 	translateCurve(vTrans);
	int i;
	MATRIX2X2 rotMat((Real)cos(-m_dAngle), (Real)-sin(-m_dAngle), (Real)sin(-m_dAngle), (Real)cos(-m_dAngle) );
 	for(i = 0; i < m_inumPoints; i++)
 	{
 		m_ControlPoints[i] = rotMat * m_ControlPoints[i];
 	}
	

	int m = (int)m_Boxes.size();
	/* rotate segment boxes */
	for(i = 0; i < m; i++)
	{
		m_Boxes[i]->m_Vertices[3] = rotMat * m_Boxes[i]->m_Vertices[3];
		m_Boxes[i]->m_Vertices[2] = rotMat * m_Boxes[i]->m_Vertices[2];
		m_Boxes[i]->m_Vertices[0] = rotMat * m_Boxes[i]->m_Vertices[0];
		m_Boxes[i]->m_Vertices[1] = rotMat * m_Boxes[i]->m_Vertices[1];
	}//end for

	/* rotate curve bounding box */
	m_bBox->m_Vertices[0] = rotMat * m_bBox->m_Vertices[0];
	m_bBox->m_Vertices[1] = rotMat * m_bBox->m_Vertices[1];

	m_dAngle += dAngle;	
	
	MATRIX2X2 rotMat2((Real)cos(m_dAngle), (Real)-sin(m_dAngle), (Real)sin(m_dAngle), (Real)cos(m_dAngle) );
		////myVec = rotMat * myVec;
	for(i = 0; i < m_inumPoints; i++)
	{
		m_ControlPoints[i] = rotMat2 * m_ControlPoints[i];
	}

	/* rotate segment boxes */
	for(i = 0; i < m; i++)
	{
		m_Boxes[i]->m_Vertices[3] = rotMat2 * m_Boxes[i]->m_Vertices[3];
		m_Boxes[i]->m_Vertices[2] = rotMat2 * m_Boxes[i]->m_Vertices[2];
		m_Boxes[i]->m_Vertices[0] = rotMat2 * m_Boxes[i]->m_Vertices[0];
		m_Boxes[i]->m_Vertices[1] = rotMat2 * m_Boxes[i]->m_Vertices[1];
	}

	/* rotate curve bounding box */

	m_bBox->m_Vertices[0] = rotMat2 * m_bBox->m_Vertices[0];
	m_bBox->m_Vertices[1] = rotMat2 * m_bBox->m_Vertices[1];
		
	if(m_dAngle > (2 * 3.1415926535))
		m_dAngle = 0.0;
	
	translateCurve(oldCog);
	
}//end rotateCurve

void CBSpline::destroyHierarchy()
{

	for(int i = 0; i < m_iBoxes; i++)
	{
		
		m_Hierarchy[i]->deleteSubTree(m_Hierarchy[i]->getRoot());

	}//end for

}//end destroy Hierarchy

/////////////////////////////////////////////////////
//
// prints knot vector and control points to a file
//
/////////////////////////////////////////////////////
void CBSpline::toFile(const char *fileName)
{
	// Creates an ofstream object named myFile
	ofstream myFile(fileName,ios::out);
    myFile.precision(16); 

    if (! myFile) // Always test file open
    {
        cout << "Error opening output file" << endl;
        return;
    }

	int i;
	int num = getNumPoints();
	myFile <<num<<endl;
	for(i = 0; i < num ; i++)
	{
		myFile << m_ControlPoints[i].x<<endl;
		myFile << m_ControlPoints[i].y<<endl;
	}//end for i

	
	num = getNumKnots();
	myFile <<num<<endl;
	for(i = 0; i < num; i++)
	{
		myFile <<m_KnotVector[i] <<endl;
	}
    myFile.close();
}


///////////////////////////////////////////////////
//
//	returns the knot span represented
//	by a curve bounding box
//
//
//
///////////////////////////////////////////////////

Real CBSpline::getIntervalLength(CSBB *box)
{

	if(box->m_fP1 <= box->m_fP2)
	{
		return (box->m_fP2 - box->m_fP1);
	}
	else
	{
		return (m_KnotVector[m_inumPoints] - box->m_fP1 + box->m_fP2 - m_KnotVector[m_ideg]);
	}

	return -1.0;

}//end getIntervalLength

int CBSpline::findKnotSpan(Real u_v)
{
	if(u_v == m_KnotVector[m_inumPoints])
		return m_inumPoints - 1;

	for(int i=m_ideg;i<m_inumKnots - m_ideg;i++)
	{
		if(u_v >= m_KnotVector[i] && u_v < m_KnotVector[i+1])
		{
			
			return i;
		}
	}
	return -1;
}//end findKnotSpan

int CBSpline::getMultiplicity(Real u)
{
	int m=0;
	for(int i=0;i<m_inumPoints;i++)
	{
		if(m_KnotVector[i] == u)
			m++;
	}
	return m;
}//end GetMultiplicity

void CBSpline::genHull()
{
	graham grH;
	VECTOR2 *cp = new VECTOR2[m_inumPoints];
	for(int i = 0; i < m_inumPoints; i++)
	{
		cp[i] = m_ControlPoints[i];
	}
	grH.getPoints(cp, m_inumPoints,1000);
	m_Hull = grH.compConvexHull();
	delete[] cp;
}//end getHull

/////////////////////////////////////////////////////
//
//	overloaded operator << to output BSpline data
//
//
/////////////////////////////////////////////////////
ostream& operator<<(ostream& out, CBSpline *pCurve) 
{

	int i;
	int num = pCurve->getNumPoints();
	out <<"CONTROL POINTS:  "<<endl;
	for(i = 0; i < num ; i++)
	{
		out << pCurve->GetCP(i);
	}//end for i

	out <<"KNOT VECTOR:  "<<endl;
	num = pCurve->getNumKnots();
	for(i = 0; i < num; i++)
	{
		out <<pCurve->GetKnot(i) <<endl;
	}

	return out;

}//end operator<<





/*!
    \fn CBSpline::GenSamples(std::vector<VECTOR2> &samples, int nS, Real start, Real end)
 */
void CBSpline::GenSamples(std::vector<VECTOR2> &samples, int nS, Real start, Real end)
{
 
}

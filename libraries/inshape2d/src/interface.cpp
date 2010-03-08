/***************************************************************************
 *   Copyright (C) 2007 by Raphael Mnster   *
 *   raphael@Reynolds   *
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

#include "interface.h"
#include "feat2ops.h"
#include <string>

//================================================
//		    		GLOBALS
//================================================

/* VECTOR OF APPROX CURVES */
vector<CApproxCurve*> ApproxCurves;

/* INITIALIZATION STATE */
bool g_Init = false;

/* NUMBER OF OJECTS */
int g_nObjects = 0;

//Bounding Volume Hierarchy
CAABBTree *BVHTree = NULL;


//================================================
//				FUNCTION DEFINITION
//================================================


//==========================================================================================
//
//		Init functions: 
//
//						1. initgroup: initializes a nurbs group
//		
//
//==========================================================================================

extern "C" void initgroup(int *iNumObj)
{

	/* initialize OpenNurbs */	
	ON::Begin();
	
	/* initialize curves */
	int numPoints = 0;
	
	//name of the nurbs data file

	//two matchstickmen
	//const char *sFileName ="pre/matchstickman2.3dm";

	const char *sFileName ="pre/120particles.3dm";
	//const char *sFileName ="pre/particles_featflow2.3dm";


	
	//load models
	ONX_Model model;

	FILE* archive_fp = ON::OpenFile( sFileName, "rb");
	if ( !archive_fp ) 
	{
		cerr<<"Could not open archive "<<sFileName<<endl;
		exit(0);
	}//end if

	const ON_NurbsCurve* curve = 0;

	ON_BinaryFile archive(ON::read3dm, archive_fp);

	bool rc = model.Read(archive);

    // close the file
	ON::CloseFile( archive_fp );
	
	if(!rc)
		cout<<"Error during reading"<<endl;

	if(!model.IsValid())
		cout<<"Model is not valid"<<endl;

	//set the number of objects
	g_nObjects = model.m_object_table.Count();

	//set the reference parameter
	*iNumObj = g_nObjects;

	//printf("Number of objects %d\n",model.m_object_table.Count());	
 
	for(int i = 0; i < g_nObjects; i++)
	{

		const ONX_Model_Object mo = model.m_object_table[i];

		curve = ON_NurbsCurve::Cast(mo.m_object);
		
		const ON_NurbsCurve &rcurve = *curve;

		CONXNurbs *tempCurve = new CONXNurbs(rcurve);
		
		tempCurve->genBoxes(1500);
		tempCurve->subdivideBoxes();
		tempCurve->subdivideBoxes();
		tempCurve->subdivideBoxes();
		CApproxCurve *myCurve = new CApproxCurve( tempCurve, 200, static_cast<short>(i) );
		ApproxCurves.push_back(myCurve);
	}//end for
	
	/* release model */
	model.Destroy();

	/* Build Object BVH */
	if(g_nObjects > 1)
		BVHTree = new CAABBTree(ApproxCurves, CAABBTree::SPATIAL_MEDIAN);


	/* set initialization variable */
	g_Init = true;	

	cout<<"Initialized a group of nurbs curves..."<<endl;
	cout<<"Number of curves: "<<g_nObjects<<endl;
	cout<<"Name of data file: "<<sFileName<<endl;


}//end initgroup

extern "C" void initcurvefromfile(int *ilength,char *sfile)
{
	/* initialize OpenNurbs */	
	ON::Begin();
	
	/* initialize curves */
	int numPoints = 0;
	int length=*ilength;
	printf("Initializing Nurbs module .\n");
	std::string sfileName(sfile);
	sfileName=sfileName.substr(0,length);

	//load models
	ONX_Model model;

	FILE* archive_fp = ON::OpenFile( sfileName.c_str(), "rb");
	if ( !archive_fp ) 
	{
     
		cerr<<"Could not open archive "<<sfileName<<endl;
		exit(0);
	}//end if

	const ON_NurbsCurve* curve = 0;

	ON_BinaryFile archive(ON::read3dm, archive_fp);

	bool rc = model.Read(archive);

    // close the file
	ON::CloseFile( archive_fp );
	
	if(rc)
		cout<<"Archive successfully read."<<endl;
	else
		cout<<"Error during reading"<<endl;

	if(model.IsValid())
		cout<<"Model is valid"<<endl;
	else
		cout<<"Model is not valid"<<endl;

	//set the number of objects
	g_nObjects = model.m_object_table.Count();

	printf("Number of objects %d\n",model.m_object_table.Count());	
 
	for(int i = 0; i < g_nObjects; i++)
	{

		const ONX_Model_Object mo = model.m_object_table[i];

		curve = ON_NurbsCurve::Cast(mo.m_object);
		
		const ON_NurbsCurve &rcurve = *curve;

		CONXNurbs *tempCurve = new CONXNurbs(rcurve);
		
		tempCurve->genBoxes(1500);
		tempCurve->subdivideBoxes();
		tempCurve->subdivideBoxes();
		tempCurve->subdivideBoxes();
		CApproxCurve *myCurve = new CApproxCurve( tempCurve, 400, static_cast<short>(i) );
		ApproxCurves.push_back(myCurve);
	}//end for
	
	/* release model */
	model.Destroy();

	/* set initialization variable */
	g_Init = true;	

	/* translate the curve so that the COG is at the origin */
	VECTOR2 vCOG=ApproxCurves[0]->GetCenter();
	vCOG=vCOG*-1.0;
	ApproxCurves[0]->TranslateCurveByVector(vCOG);

	printf("\n");

	printf("Nurbs module initialized\n");
	
}//end initcurve


extern "C" void initcurve()
{
	/* initialize OpenNurbs */	
	ON::Begin();
	
	/* initialize curves */
	int numPoints = 0;
	printf("Initializing Nurbs module .\n");

	//const char *sFileName ="pre/bommerholz.3dm";

	//const char *sFileName ="pre/featflowklein4.3dm";
	
	//const char *sFileName ="pre/featcrashSmall.3dm";

	//const char *sFileName ="pre/kreis0502_01_hotzenplotz.3dm";

	//const char *sFileName ="pre/bubble.3dm";

	//const char *sFileName ="pre/matchstickman.3dm";

	const char *sFileName ="pre/bubblesml.3dm";
	
	//load models
	ONX_Model model;

	FILE* archive_fp = ON::OpenFile( sFileName, "rb");
	if ( !archive_fp ) 
	{
     
		cerr<<"Could not open archive "<<sFileName<<endl;
		exit(0);
	}//end if

	const ON_NurbsCurve* curve = 0;

	ON_BinaryFile archive(ON::read3dm, archive_fp);

	bool rc = model.Read(archive);

    // close the file
	ON::CloseFile( archive_fp );
	
	if(rc)
		cout<<"Archive successfully read."<<endl;
	else
		cout<<"Error during reading"<<endl;

	if(model.IsValid())
		cout<<"Model is valid"<<endl;
	else
		cout<<"Model is not valid"<<endl;

	//set the number of objects
	g_nObjects = model.m_object_table.Count();

	printf("Number of objects %d\n",model.m_object_table.Count());	
 
	for(int i = 0; i < g_nObjects; i++)
	{

		const ONX_Model_Object mo = model.m_object_table[i];

		curve = ON_NurbsCurve::Cast(mo.m_object);
		
		const ON_NurbsCurve &rcurve = *curve;

		CONXNurbs *tempCurve = new CONXNurbs(rcurve);
		
		tempCurve->genBoxes(1500);
		tempCurve->subdivideBoxes();
		tempCurve->subdivideBoxes();
		tempCurve->subdivideBoxes();
		CApproxCurve *myCurve = new CApproxCurve( tempCurve, 1000, static_cast<short>(i) );
		ApproxCurves.push_back(myCurve);
	}//end for
	
	/* release model */
	model.Destroy();

	/* set initialization variable */
	g_Init = true;	

	printf("\n");

	printf("Nurbs module initialized\n");
	
}//end initcurve

extern "C" void initobject(int *nID)
{

	/* initialize OpenNurbs */	
	ON::Begin();
	
	/* initialize curves */
	int numPoints = 0;
	printf("Initializing Nurbs module .\n");

	const char *sFileName ="pre/bubblesml.3dm";
	//const char *sFileName ="pre/matchstickman.3dm";
	
	//load models
	ONX_Model model;

	FILE* archive_fp = ON::OpenFile( sFileName, "rb");
	if ( !archive_fp ) 
	{
     
		cerr<<"Could not open archive "<<sFileName<<endl;
		exit(0);
	}//end if

	const ON_NurbsCurve* curve = 0;

	ON_BinaryFile archive(ON::read3dm, archive_fp);

	bool rc = model.Read(archive);

    // close the file
	ON::CloseFile( archive_fp );
	
	if(rc)
		cout<<"Archive successfully read."<<endl;
	else
		cout<<"Error during reading"<<endl;

	if(model.IsValid())
		cout<<"Model is valid"<<endl;
	else
		cout<<"Model is not valid"<<endl;

	//set the number of objects
	g_nObjects = model.m_object_table.Count();

	printf("Number of objects %d\n",model.m_object_table.Count());	
 
	for(int i = 0; i < g_nObjects; i++)
	{

		const ONX_Model_Object mo = model.m_object_table[i];

		curve = ON_NurbsCurve::Cast(mo.m_object);
		
		const ON_NurbsCurve &rcurve = *curve;

		CONXNurbs *tempCurve = new CONXNurbs(rcurve);
		
		tempCurve->genBoxes(1500);
		tempCurve->subdivideBoxes();
		tempCurve->subdivideBoxes();
		tempCurve->subdivideBoxes();
		CApproxCurve *myCurve = new CApproxCurve( tempCurve, 200, static_cast<short>(i) );
		ApproxCurves.push_back(myCurve);
	}//end for
	
	/* release model */
	model.Destroy();

	/* set initialization variable */
	g_Init = true;	

	printf("\n");

	printf("Nurbs module initialized\n");

	int ID =  (int)ApproxCurves.size();

	*nID = ID;

}//end initobject


extern "C" int isingeometry2(double *x, double *y)
{

	double dX = *x;
	double dY = *y;
	VECTOR2 vQuery(dX,dY);

	int res;

	bool result = false;

	if(ApproxCurves[0]->m_bBox->inside(vQuery))
	{
	
		if(ApproxCurves[0]->IsInElement(vQuery))
		{
	
			result = true;
	
		}
	}

	if(result == true)
	{
		return 1;
	}
	else
	{
		return 0;
	}
	
}//end isingeometry

extern "C" int getnumverts()
{

	if(ApproxCurves[0])
		return ApproxCurves[0]->GetResolution();
	else
		return 0;

}//end getnumverts

extern "C" Real getvalx(int *iID)
{
	//cout<<"X-coordinate: "<<ApproxCurves[0]->GetSamples()[*iID].x<<endl;
	return ApproxCurves[0]->GetSamples()[*iID].x;
}//end getvalx

extern "C" Real getvaly(int *iID)
{
	return ApproxCurves[0]->GetSamples()[*iID].y;
}//end getvaly

extern "C" void freememory()
{

	if(!ApproxCurves.empty())
	{

		vector<CApproxCurve*>::iterator vIter;

		for(vIter = ApproxCurves.begin(); vIter != ApproxCurves.end(); vIter++)
		{
			CApproxCurve *pCurve = *vIter;
			delete pCurve;
		}//end for

		ApproxCurves.clear();

	}//end if

	if(BVHTree)
	{
		delete BVHTree;
		BVHTree = NULL;
	}//end if

}//end freememory


extern "C" int isingroup(Real *x, Real *y)
{

	//make a vector object for the grid point
	VECTOR2 vQuery(*x, *y);

	//list of AABB candidates
	list<CAABBNode*> lCandidates;

	//iterator object for the curves
	vector<CApproxCurve*>::iterator vIter;

	//a class to fulfill Feat2 operations
	CFeat2Ops op;

	//loop through all curves
	for(vIter = ApproxCurves.begin(); vIter != ApproxCurves.end(); vIter++)
	{
		//get the current curve
		CApproxCurve *pCurve = *vIter;

		//if the point is inside the curve
		if(op.IsFictitiousBoundary(pCurve, vQuery))
		{
			//point is inside, return...
			return 1;
			cout<<"inside"<<endl;
		}//end if
	}//end for

	//return the result of the Inside/Outside test
	return 0;

}//end isingroup

extern "C" int getnumvertsid(int *iID)
{
	return ApproxCurves[*iID]->GetResolution();
}//end getnumvertsid

extern "C" Real getvalxid(int *iID, int *nX)
{
	//cout<<"X-coordinate: "<<ApproxCurves[*iID]->GetSamples()[*nX].x<<endl;
	return ApproxCurves[*iID]->GetSamples()[*nX].x;
	
}//end getvalxid

extern "C" Real getvalyid(int *iID, int *nY)
{
	return ApproxCurves[*iID]->GetSamples()[*nY].y;
}//end getvalyid

extern "C" void rotateobject(double *dpAngle)
{

	double dAngle = *dpAngle;

	MATRIX2X2 rotMat((double)cos(dAngle), (double)-sin(dAngle), (double)sin(dAngle), (double)cos(dAngle) );

	ApproxCurves[0]->RotateCurve(rotMat);

}//end rotateobject

extern "C" void translateobject(double *dpTransX, double *dpTransY)
{

	VECTOR2 vNewCenter(*dpTransX, *dpTransY);
	
	ApproxCurves[0]->SetCenter(vNewCenter);

	//VECTOR2 vTrans =  vNewCenter - ApproxCurves[0]->GetCenter();

	//ApproxCurves[0]->TranslateCurveByVector(vTrans);

}//End translateobject

extern "C" void getdistance(double *dX, double *dY,double *ddist)
{
	double x = *dX;
	double y = *dY;
	VECTOR2 lastPoint(0,0);
	VECTOR2 vQuery(x,y);
	double dLUB = 1.7E+308;
	DistOps op;
	int iConv = 0;

	list<CBCNode*> lBFS;
	CBCNode *nBest = NULL;
	SearchResult sres=op.LUBSearch(ApproxCurves[0], vQuery, lBFS, nBest);
	*ddist=sres.res;

}

extern "C" void getcenter(double *dX,double *dY)
{
	
	VECTOR2 vCenter=ApproxCurves[0]->GetCenter();
	*dX=vCenter.x;
	*dY=vCenter.y;

}

extern "C" void getvertex(int *nX,double *dX,double *dY)
{
	int index =*nX;
	index-=1;
	VECTOR2 vPoint= ApproxCurves[0]->GetSamples()[*nX];

	*dX=vPoint.x;
	*dY=vPoint.y;

}
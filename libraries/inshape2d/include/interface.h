

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

//===================================================
//					DEFINES
//===================================================


#if !defined(_CPPINTERFACE_H)
#define _CPPINTERFACE_H

//================================================
//		    		INCLUDES
//================================================
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <list>
#include <math.h>
#include "mathglobals.h"
#include "nurbs.h"
#include "distops.h"
#include "aabbtree.h"
#include "fileParser.h"
#include "approxcurve.h"
#include "obb2.h"
//================================================
//		      FUNCTION DECLARATION
//================================================


using namespace std;

//==================================================================
//
//		function names in small letters without underscores
//
//==================================================================

extern "C" void initcurvefromfile(int *ilength,char *sfile);

extern "C" void initcurve();

extern "C" void initobject(int *nID);

extern "C" int isingeometry2(double *x, double *y);

extern "C" void isingeometry(double *x, double *y,int *iin)
{
	*iin=isingeometry2(x,y);
}

extern "C" int getnumverts();

extern "C" Real getvalx(int *iID);

extern "C" Real getvaly(int *iID);

extern "C" int getnumvertsid(int *iID);

extern "C" Real getvalxid(int *iID, int *nX);

extern "C" Real getvalyid(int *iID, int *nX);

extern "C" void getvertex(int *nX,double *dX,double *dY);

extern "C" void rotateobject(double *dpAngle);

extern "C" void initgroup(int *iNumObj);

extern "C" int isingroup(Real *x, Real *y);

extern "C" void translateobject(double *dpTransX, double *dpTransY);

extern "C" void freememory();

extern "C" void getcenter(double *dX,double *dY);

extern "C" void getdistance(double *dX, double *dY,double *ddist);

//==================================================================
//
//		function names in small letters with underscores
//
//==================================================================

extern "C" void getvertex_(int *nX,double *dX,double *dY)
{
	getvertex(nX,dX,dY);
}

extern "C" void getcenter_(double *dX,double *dY)
{
	getcenter(dX,dY);
}

extern "C" void initcurvefromfile_(int *ilength,char *sfile)
{
	initcurvefromfile(ilength,sfile);
}

extern "C" void getdistance_(double *dX, double *dY,double *ddist)
{
	getdistance(dX,dY,ddist);
}

extern "C" void initcurve_()
{
	initcurve();
}

extern "C" void initobject_(int *nID)
{
	initobject(nID);
}

extern "C" void isingeometry_(double *x, double *y,int *iin)
{
	*iin = isingeometry2(x,y);
}

extern "C" int getnumverts_()
{
	return getnumverts();
}

extern "C" Real getvalx_(int *iID)
{
	return getvalx(iID);
}

extern "C" Real getvaly_(int *iID)
{
	return getvaly(iID);
}

extern "C" int getnumvertsid_(int *iID)
{
	return getnumvertsid(iID);
}

extern "C" Real getvalxid_(int *iID, int *nX)
{
	return getvalxid(iID,nX);
}

extern "C" Real getvalyid_(int *iID, int *nX)
{
	return getvalyid(iID, nX);
}

extern "C" void rotateobject_(double *dpAngle)
{
	rotateobject(dpAngle);
}

extern "C" void initgroup_(int *iNumObj)
{
	initgroup(iNumObj);
}

extern "C" int isingroup_(Real *x, Real *y)
{
	return isingroup(x,y);
}

extern "C" void translateobject_(double *dpTransX, double *dpTransY)
{
	translateobject(dpTransX, dpTransY);
}

extern "C" void freememory_()
{
	freememory();
}

//==================================================================
//
//		function names in capital letters with underscores
//
//==================================================================

extern "C" void GETVERTEX_(int *nX,double *dX,double *dY)
{
	getvertex(nX,dX,dY);
}

extern "C" void GETCENTER_(double *dX,double *dY)
{
	getcenter(dX,dY);
}


extern "C" void INITCURVEFROMFILE_(int *ilength,char *sfile)
{
	initcurvefromfile(ilength,sfile);
}

extern "C" void GETDISTANCE_(double *dX, double *dY,double *ddist)
{
	getdistance(dX,dY,ddist);
}

 extern "C" void INITCURVE_()
 {
 	initcurve();
 }
 
 extern "C" void INITOBJECT_(int *nID)
 {
    initobject(nID);
 }
 
 extern "C" void ISINGEOMETRY_(double *x, double *y,int *iin)
 {
 	*iin=isingeometry2(x,y);
 }
 
 extern "C" int GETNUMVERTS_()
 {
 	return getnumverts();
 }
 
 extern "C" Real GETVALX_(int *iID)
 {
 	return getvalx(iID);
 }
 
 extern "C" Real GETVALY_(int *iID)
 {
 	return getvaly(iID);
 }
 
 extern "C" int GETNUMVERTSID_(int *iID)
 {
 	return getnumvertsid(iID);
 }
 
 extern "C" Real GETVALXID_(int *iID, int *nX)
 {
 	return getvalxid(iID,nX);
 }
 
 extern "C" Real GETVALYID_(int *iID, int *nX)
 {
 	return getvalyid(iID, nX);
 }
 
 extern "C" void ROTATEOBJECT_(double *dpAngle)
 {
 	rotateobject(dpAngle);
 }
 
 extern "C" void INITGROUP_(int *iNumObj)
 {
 	initgroup(iNumObj);
 }
 
 extern "C" int ISINGROUP_(Real *x, Real *y)
 {
 	return isingroup(x,y);
 }
 
 extern "C" void TRANSLATEOBJECT_(double *dpTransX, double *dpTransY)
 {
 	translateobject(dpTransX, dpTransY);
 }
 
 extern "C" void FREEMEMORY_()
 {
 	freememory();
 }

//==================================================================
//
//		function names in capital letters without underscores
//
//==================================================================

extern "C" void GETVERTEX(int *nX,double *dX,double *dY)
{
	getvertex(nX,dX,dY);
}

extern "C" void GETCENTER(double *dX,double *dY)
{
	getcenter(dX,dY);
}


extern "C" void INITCURVEFROMFILE(int *ilength,char *sfile)
{
	initcurvefromfile(ilength,sfile);
}

extern "C" void GETDISTANCE(double *dX, double *dY,double *ddist)
{
	getdistance(dX,dY,ddist);
}

extern "C" void INITCURVE()
 {
 	initcurve();
 }
 
 extern "C" void INITOBJECT(int *nID)
 {
 	initobject(nID);
 }
 
 extern "C" void ISINGEOMETRY(double *x, double *y,int *iin)
 {
 	*iin=isingeometry2(x,y);
 }
 
 extern "C" int GETNUMVERTS()
 {
 	return getnumverts();
 }
 
 extern "C" Real GETVALX(int *iID)
 {
 	return getvalx(iID);
 }
 
 extern "C" Real GETVALY(int *iID)
 {
 	return getvaly(iID);
 }
 
 extern "C" int GETNUMVERTSID(int *iID)
 {
 	return getnumvertsid(iID);
 }
 
 extern "C" Real GETVALXID(int *iID, int *nX)
 {
 	return getvalxid(iID,nX);
 }
 
 extern "C" Real GETVALYID(int *iID, int *nX)
 {
 	return getvalyid(iID, nX);
 }
 
 extern "C" void ROTATEOBJECT(double *dpAngle)
 {
 	rotateobject(dpAngle);
 }
 
 extern "C" void INITGROUP(int *iNumObj)
 {
 	initgroup(iNumObj);
 }
 
 extern "C" int ISINGROUP(Real *x, Real *y)
 {
 	return isingroup(x,y);
 }
 
 extern "C" void TRANSLATEOBJECT(double *dpTransX, double *dpTransY)
 {
 	translateobject(dpTransX, dpTransY);
 }
 
 extern "C" void FREEMEMORY()
 {
 	freememory();
 }


//====================================================================
//
//		function names in capital letters without underscores infront
//
//====================================================================

extern "C" void _getvertex(int *nX,double *dX,double *dY)
{
	getvertex(nX,dX,dY);
}

extern "C" void _getcenter(double *dX,double *dY)
{
	getcenter(dX,dY);
}


extern "C" void _initcurvefromfile(int *ilength,char *sfile)
{
	initcurvefromfile(ilength,sfile);
}

extern "C" void _getdistance(double *dX, double *dY,double *ddist)
{
	getdistance(dX,dY,ddist);
}

 extern "C" void _initcurve()
{
	initcurve();
}

extern "C" void _initobject(int *nID)
{
	initobject(nID);
}

extern "C" void _isingeometry(double *x, double *y,int *iin)
{
	*iin=isingeometry2(x,y);
}

extern "C" int _getnumverts()
{
	return getnumverts();
}

extern "C" Real _getvalx(int *iID)
{
	return getvalx(iID);
}

extern "C" Real _getvaly(int *iID)
{
	return getvaly(iID);
}

extern "C" int _getnumvertsid(int *iID)
{
	return getnumvertsid(iID);
}

extern "C" Real _getvalxid(int *iID, int *nX)
{
	return getvalxid(iID,nX);
}

extern "C" Real _getvalyid(int *iID, int *nX)
{
	return getvalyid(iID, nX);
}

extern "C" void _rotateobject(double *dpAngle)
{
	rotateobject(dpAngle);
}

extern "C" void _initgroup(int *iNumObj)
{
	initgroup(iNumObj);
}

extern "C" int _isingroup(Real *x, Real *y)
{
	return isingroup(x,y);
}

extern "C" void _translateobject(double *dpTransX, double *dpTransY)
{
	translateobject(dpTransX, dpTransY);
}

extern "C" void _freememory()
{
	freememory();
}

#endif



***********************************************************************
* This file provides routines to calculate the approximate arc length
* for a fictitious boundary component. The computation is based
* on a volume integration approach with a special trick:
*
* Let S denote the boundary curve of a fictitious boundary component.
* Then the drag- and lift-values are normally calculated by
*
*    int ( sigma * n ) ds
*     S
*
* with the stress tensor sigma=[grad(u)+p] and the normal vector n of
* S. Now replace sigma by n and obtain:
*
*    int ( n * n ) ds = int (1) ds = arclength(s)
*     S                  S
*
* because n has length 1. (To be more exact: we replaced the matrix
* sigma by the 2x2-matrix [n^T] and took the first component...
*                         [n^T]
*
* Using the Gauss formula we can obtain a (more or less) equivalent
* formulation using the volume integration and the characteristic
* function alpha of the fictitious boundary component V:
*
*    arclength(s) = int ( n * n ) ds = int (n * grad(alpha)) dx
*                    S                  V
*
* This is what we evaluate here.
*
* On the other hand the length of the boundary is calculated by
* adding the reconstructed interface of the fictitious boundary
* directly.
***********************************************************************

************************************************************************
* Calculates approximate arc length
*
* Volume integration, constant pressure
*
* In:
*  KVERT,KNPR,KMID,DCORVG,DCORMG - usual geometry data
*  ELE    - Element function used for volume integration
*  IFBC   - Number of the boundary component whose arc length should
*           be computed.
*           =0: compute sum of lengths of all fictitious boundary 
*               components
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
*
* Out:
*  DACLNV - Approximation of the arc length by volume integration
*  DACLNR - Approximation of the arc length by reconstructed interface
************************************************************************

      SUBROUTINE BDARCL (DCORVG,KVERT,DCORMG,KMID,KNPR,NEL,
     *                   NVE,NVT,NBCT,ELE,IFBC,DACLNV,DACLNR,LAMBDA,
     *                   IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'casmbly.inc'
      
      INCLUDE 'cmgadr.inc'

      INCLUDE 'cinidat.inc'
      
      INCLUDE 'ccub.inc'

C parameters
      DOUBLE PRECISION DCORVG(2,*),DCORMG(2,*),DGEOM(*)
      INTEGER KNPR(*),KVERT(NNVE,*),KMID(NNVE,*),IFBC,IGEOM(*)
      INTEGER NEL,NVE,NVT,NBCT
      DOUBLE PRECISION DACLNV,DACLNR,LAMBDA

C externals
      EXTERNAL ELE
      EXTERNAL CAINCP,CAVDCP,CAISCP,CADOCP,CAVDCE
      
C local variables

      INTEGER IINFO(1)
      
C Dummy variables

      DOUBLE PRECISION DU1(1),DU2(1),DP(1)
      
C Handle to DINFO-Array for integration routines.
C
C Structure:
C   double DACLNV - approximation of the arc length
C   double DACLNR - approximation of the arc length by 
C                   reconstructed interface
C   double temp [5], used for calculation
C   double ALPHA [NU]

      INTEGER LDINFO
      
C Allocate memory for DINFO-arrays. We don't need IINFO!
C NU is taken from the COMMON block.

      CALL ZNEW (NU+8,1,LDINFO,'LDINFO')
      
c Use the vector dalpha as a test function in the weak formulation
c (only the difusion part is taken, should take even the
c convection but in featflow it seems to be ok) ... 
c this needs to be checked in timedependent case

C Prepare the ALPHA-vector for all fict. boundary components;
C Use the INALPH-function from DRAGLIFT2.F:

      CALL INALPH (DWORK(L(LDINFO)+7),KNPR,KMID,NEL,NVT,NBCT,NU,0)
      
C Relaxate the alpha

      CALL RXALPH (DWORK(L(LDINFO)+7),DCORVG,KVERT,DCORMG,KMID,NEL,NVT,
     *             NU,IFBC,LAMBDA,IGEOM,DGEOM)
      
C IINFO holds the number of the fictitious boundary component

      IINFO(1) = IFBC
      
C Call the calculation routine; IINFO is a dummy parameter.
C Because this is a special volume integration not using velocity or
C pressure, we are using the pure user defined calculation
C routine here. DU1,DU2 and DP are dummies...
      
      CALL VOLINT (DU1,DU2,DP,KVERT,KMID,
     *    DCORVG,NEL,NVE,ELE,(IELT.EQ.2).OR.(IELT.EQ.3),8,
     *    IINFO, 1, DWORK(L(LDINFO)), NU+8,
     *    CAINCP,CAVDCE,CAVDCP,CAISCP,CADOCP,.TRUE.,0,0)
C     Call wrong!!!
      
C Obtain the values of the arc length

      DACLNV = DWORK(L(LDINFO))
      DACLNR = DWORK(L(LDINFO)+1)
      
C Release memory, finish

      CALL ZDISP (0,LDINFO,'LDINFO')
      
99999 END


************************************************************************
* Volume-integration callback-routines for arc length approximation:
************************************************************************

************************************************************************
* Initialisation
************************************************************************
      SUBROUTINE CAINCP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                   ELE,BNONPR)
      
      IMPLICIT NONE
      
C parameters
      
      INTEGER NIINFO,IINFO(NIINFO),NDINFO
      LOGICAL BNONPR
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      EXTERNAL ELE
      
C initialise everything with 0
      
      DINFO(1) = 0D0
      DINFO(2) = 0D0
      
      END

************************************************************************
* Postprocessing
************************************************************************
      SUBROUTINE CADOCP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                   ELE,BNONPR)
      
      IMPLICIT NONE
      
C parameters
      
      INTEGER NIINFO,IINFO(NIINFO),NDINFO
      LOGICAL BNONPR
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      EXTERNAL ELE

C nothing is done here      

      END
      
************************************************************************
* Calculation on an element
************************************************************************
      SUBROUTINE CAVDCE (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                   IELE,ELE,BNONPR)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
C parameters

      INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELE
      LOGICAL BNONPR
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      EXTERNAL ELE

C local variables

      INTEGER IFBC, IEDGE1, IEDGE2
      DOUBLE PRECISION X1,Y1, X2,Y2, LEN
      LOGICAL BFOUND
      
C Number of fictitious boundary component

      IFBC = IINFO(1)

C This performes the approximation by reconstruction of the interface.
C Call the reconstruction routine, reconstruct the points by 8 
C approximation steps.

      CALL RCFBLI (DWORK(L(LCORVG)), KWORK(L(LVERT)), IELE, 8, 
     *             IFBC, .FALSE., X1,Y1, IEDGE1, X2,Y2, IEDGE2,BFOUND,
     *             IGEOM,DGEOM)
     
C store the length

      IF (BFOUND) THEN
        LEN = DSQRT ((X2-X1)**2+(Y2-Y1)**2)
        DINFO(2) = DINFO(2)+LEN
      END IF
      
      END
      
************************************************************************
* Calculation of (partial) values in a cubature point
************************************************************************
      SUBROUTINE CAVDCP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                   XX,YY,XX1,YY1,ICUBP,IELE,ELE,BNONPR)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'casmbly.inc'
      
C parameters
      
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELE,ICUBP
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      DOUBLE PRECISION XX,YY,XX1,YY1
      LOGICAL BNONPR
      EXTERNAL ELE

C local variables

      DOUBLE PRECISION DAV,DAX,DAY,ALPH
      INTEGER I,IG,IFBC

C The ALPHA-vector is stored in DINFO(6..NU+6)

      DAV=0D0 ! VALUE
      DAX=0D0 ! X-der.
      DAY=0D0 ! Y-der.
      DO I=1,IDFL
        IG=KDFG(I)
        ALPH=DINFO(5+IG)
        DAV=DAV+ALPH*DBAS(KDFL(I),1)
        DAX=DAX+ALPH*DBAS(KDFL(I),2)
        DAY=DAY+ALPH*DBAS(KDFL(I),3)
      END DO

C save these values to DINFO for later use

      DINFO(3) = DAV
      DINFO(4) = DAX
      DINFO(5) = DAY

C Number of fictitious boundary component

      IFBC = IINFO(1)
      
C Try to obtain an approximation of the normal vector in the quadrature
C point and save that in DINFO

      CALL FBDNML (XX,YY,IFBC,DINFO(6),DINFO(7))

      END
      
************************************************************************
* Summing up the values from the cubature point to the values of
* interest.
************************************************************************
      SUBROUTINE CAISCP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                   DETJ,ICUBP,IELE,OMEGA,
     *                   DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'casmbly.inc'
      
C parameters
      
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,ICUBP,IELE
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      DOUBLE PRECISION DETJ,OMEGA
      DOUBLE PRECISION DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV
      
C local variables

      DOUBLE PRECISION DAV,DAX,DAY,DN1,DN2,AH1

C get value/derivative

      DAV = DINFO(3)
      DAX = DINFO(4)
      DAY = DINFO(5)
      
C get the normal vector from the calculation routine above

      DN1 = DINFO(6)
      DN2 = DINFO(7)
      
C Oops, if user defined routine was not able to calculate that,
C we can't calculate the arc length...

      IF ((DN1.EQ.0D0).AND.(DN2.EQ.0D0)) GOTO 99999
      
C form the integrand: n*grad(alpha)

      AH1 = DN1*(-DAX)+DN2*(-DAY)
      
C Sum up to the integral

      DINFO(1)=DINFO(1)+AH1*OMEGA

99999 END
      
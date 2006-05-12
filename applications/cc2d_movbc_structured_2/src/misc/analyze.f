***********************************************************************
* This file contains various functions for analysing different
* mathematical algorithms. The routines in this file are more
* or less independent of the actual CC2D software package,
* i.e. they are not required for a proper operation of the
* simulation. Instead they allow to analyse special aspects
* of different subroutines, i.e. prescribing an analytical
* solution in the calculation of drag- and lift-forces 
* to analyse the accuracy of the calculation routines.
*
* The routines here are not planned of being modular - they are
* just indended to serve for some "quick-tests".
***********************************************************************

***********************************************************************
* User-defined analysis I
*
* Calls the drag/lift-calculation with a prescribed, analytical
* solution on the unit square to test its accuracy.
***********************************************************************

      SUBROUTINE UDAN01
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cnspts.inc'
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'ciniopt.inc'
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfiles.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgadr.inc'
      
      EXTERNAL EXXXNC
      
C local variables

      INTEGER LUP,I
      DOUBLE PRECISION DFWX, DFWY,DARCLN
      INTEGER NESTAT(5)
      INTEGER IFPOST,IFILEN,ITFILM
      
      EXTERNAL UE,UDAAF1,UDAAF2,UDAAF3
      
C Prescribe solutions

      CALL PTSDAT(0D0,DBLE(NY))

      ILEV=NLMAX
      I=2
      CALL SETLEV(I)

      CALL UDAA01 (DWORK(L(LCORMG)),KWORK(L(LMID)),DWORK(KU1),
     *             DWORK(KU2),DWORK(KP))
      
C Evaluate drag/lift
      
      DFWX = 0D0
      DFWY = 0D0
      
C      CALL  BDFVLG(DWORK(L(LU)),DWORK(L(LV)),DWORK(L(LP)),
C     *    DWORK(L(LCORVG)),KWORK(L(LVERT)),
C     *    DWORK(L(LCORMG)),KWORK(L(LMID)),KWORK(L(LNPR)),NEL,NVE,NVT,NBCT,
C     *    EXXXNC,.FALSE.,DFWX,DFWY,DPF(1),DPF(2),NESTAT,0,DVIARX)

C      CALL  BDFVLG(DWORK(L(LU)),DWORK(L(LV)),DWORK(L(LP)),
C     *    DWORK(L(LCORVG)),KWORK(L(LVERT)),
C     *    DWORK(L(LCORMG)),KWORK(L(LMID)),KWORK(L(LNPR)),NEL,NVE,NVT,NBCT,
C     *    EXXXNC,.FALSE.,DFWX,DFWY,1D0,1D0,NESTAT,0,DVIARX)

C Use DPF1=1D0, DPF2=2D0 - then it holds:
C  cD = 2*FD/DPF2 = FD
C  cL = 2*FL/DPF2 = FL
C i.e. the drag/lift coefficients are the drag/lift forces and thus can be
C compared with analytical calculated values.

      CALL  BDFVLG(DWORK(KU1),DWORK(KU2),DWORK(KP),DWORK(L(LCORVG)),
     *    KWORK(L(LVERT)),DWORK(L(LCORMG)),KWORK(L(LMID)),
     *    KWORK(L(LNPR)),NEL,NVE,NVT,NBCT,
     *    EXXXNC,(IELT.EQ.2).OR.(IELT.EQ.3),
     *    DFWX,DFWY,1D0,2D0,NESTAT,0,DVIARX)
      
      PRINT *,'Drag vol-int. :  ',DFWX
      PRINT *,'Lift vol-int. :  ',DFWY
      
      DFWX=0D0

      CALL AFBVOL (KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EXXXNC,8,0,DFWX)
      PRINT *,'AREA vol-int. :  ',DFWX

      DFWX=0D0
      DFWY=0D0

      CALL  BDARCL(DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *          DWORK(L(LCORMG)),KWORK(L(LMID)),KWORK(L(LNPR)),
     *          NEL,NVE,NVT,NBCT,EXXXNC,0,DFWX, DARCLN,DVIARX,0,0)
      ! call not correct!!!

      PRINT *,'ARCL vol-int. :  ',DFWX

      DFWX=0D0
      DFWY=0D0
      
      CALL  BDFRIG(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *    DWORK(L(LCORVG)),KWORK(L(LVERT)),DWORK(L(LCORMG)),
     *    KWORK(L(LMID)),KWORK(L(LNPR)),NEL,NVE,EXXXNC,
     *    (IELT.EQ.2).OR.(IELT.EQ.3),
     *    DFWX,DFWY,1D0,2D0,0)

      PRINT *
      PRINT *,'Drag line-int.:  ',DFWX
      PRINT *,'Lift line-int.:  ',DFWY
      PRINT *,'ARCL line-int.:  ',DARCLN
      
      PRINT *
      
C Analytically compute the 1-function on the interface.
C Should give the arc length:

      CALL UDAN03 (UDAAF1,0,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *          DWORK(L(LCORMG)),KWORK(L(LMID)),KWORK(L(LNPR)),NEL,
     *          NVE,EXXXNC,(IELT.EQ.2).OR.(IELT.EQ.3),DFWX)

      PRINT *,'ARCL anal. l.i.: ',DFWX

C Calculate the drag with the help of a reference
C function. Should differ from the analytical solution only
C by the effects of the cubature formula:

      CALL UDAN03 (UDAAF2,0,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *          DWORK(L(LCORMG)),KWORK(L(LMID)),KWORK(L(LNPR)),NEL,
     *          NVE,EXXXNC,(IELT.EQ.2).OR.(IELT.EQ.3),DFWX)

      PRINT *,'drag anal. l.i.: ',DFWX

C Calculate the lift with the help of a reference
C function. Should differ from the analytical solution only
C by the effects of the cubature formula:

      CALL UDAN03 (UDAAF3,0,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *          DWORK(L(LCORMG)),KWORK(L(LMID)),KWORK(L(LNPR)),NEL,
     *          NVE,EXXXNC,(IELT.EQ.2).OR.(IELT.EQ.3),DFWX)

      PRINT *,'lift anal. l.i.: ',DFWX

      IFPOST=1
      IFILEN=0
      ITFILM=0
      CALL FPOST(IFPOST,IFILEN,ITFILM,UE,MT,DFWX,DFWY,0,'')
      
C      PRINT *,'Drag2: ',DFWX
C      PRINT *,'Lift2: ',DFWY
      
C      CALL UDAN02(DWORK(KU1),DWORK(KU2),DWORK(KP),
C     *            KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
C     *            EXXXNC,DFWX,DFWY)      

C      PRINT *,'Drag appr. l.i. old: ',DFWX
C      PRINT *,'Lift appr. l.i. old: ',DFWY
      
      STOP
      
      END
      
***********************************************************************
* User-defined auxiliary routine I
*
* Set the values of velocity/pressure to analytical values.
*
* Warning: INTPOL overwrites dirichlet boundary values by reference
* solution for GMV-output, so when using this routine, disable
* the overwriting of the boundary values in the call to this
* subroutine - otherwise the boundary looks strange in GMV-files ;)
***********************************************************************

      SUBROUTINE UDAA01 (DCORMG,KMID,U,V,P)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
C parameters
      
      DOUBLE PRECISION U(*),V(*),P(*),DCORMG(2,*)
      INTEGER KMID(NNVE,*)
      
C Reference function on the unit square:
      
      DOUBLE PRECISION X,Y,UF,VF,PF
      
      UF(X,Y) = y*(1D0-x)*(0D0-x)*(1D0-y)*(0D0-y)*exp(x)
C x*(1-x)*exp(x)
C 
      VF(X,Y) = -x*(1D0-x)*(0D0-x)*(1D0-y)*(0D0-y)*exp(y)
C y*(1-y)*exp(y)
C 
      PF(X,Y) = x
      
C local variables
      
      INTEGER IVT,I,IEL
      
      DO IVT=1,NMT
        U(IVT) = UF(DCORMG(1,IVT),DCORMG(2,IVT))
        V(IVT) = VF(DCORMG(1,IVT),DCORMG(2,IVT))
      END DO

      DO IEL=1,NEL
        X=0D0
        Y=0D0
        DO I=1,4
          X = X+DCORMG(1,KMID(I,IEL)-NVT)
          Y = Y+DCORMG(2,KMID(I,IEL)-NVT)
        END DO
        X=0.25D0*X
        Y=0.25D0*Y
        P(IEL) = PF(X,Y)
      END DO
      
      END
      
      
************************************************************************
*    Purpose:  Calculates lift (DFW) and drag (DAW)
*    by integration along reconstructed interface
*-----------------------------------------------------------------------
      SUBROUTINE UDAN02(DU1,DU2,DP,KVERT,KMID,DCORVG,ELE,
     *                  DFW,DAW)
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'casmbly.inc'
      
      INCLUDE 'cnspts.inc'

C parameters

      DOUBLE PRECISION DU1(*),DU2(*),DP(*),DCORVG(2,*),DFW,DAW
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      EXTERNAL ELE

C externals

      INTEGER NDFL
      EXTERNAL NDFL
      
C local variables

      INTEGER IELTYP,IW1,IW2,IVBD,IVT,IVT1,IVT2,ISTOP,II,IVBD1
      INTEGER IVE,JP,JDFL,IEDGE1,IEDGE2
      LOGICAL BFOUND
      DOUBLE PRECISION DLEN,PX1,PX2,PY1,PY2,PXM,PYM,DLH
      DOUBLE PRECISION DTX,DTY,DNX,DNY,DPCONT,XX,YY,DUT

      DFW=0D0
      DAW=0D0
      DLEN=0D0

      IF ((DPF(1).EQ.0D0).OR.(DPF(2).EQ.0D0)) RETURN

      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)

      NCUBP=1
      ICUBP=1

      IW1=1
      IW2=1

      DO 100 IEL=1,NEL
      IF (ISTOP.EQ.1) GOTO 1000

      CALL RCFBLI (DCORVG, KVERT, IEL, 8, 0, .FALSE.,
     *             PX1,PY1, IEDGE1, PX2,PY2, IEDGE2, BFOUND,0,0)
C     Call not correct!!!

      IF (.NOT.BFOUND) GOTO 100

      PXM=0.5D0*(PX1+PX2)
      PYM=0.5D0*(PY1+PY2)
      DLH=DSQRT((PX2-PX1)**2+(PY2-PY1)**2)

      DTX= (PX2-PX1)/DLH
      DTY= (PY2-PY1)/DLH
      DNX=-(PY2-PY1)/DLH
      DNY= (PX2-PX1)/DLH
      DPCONT=DP(IEL)
      DLEN=DLEN+DLH

      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)

      DO 120 IVE = 1,4
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE

      XX=PXM
      YY=PYM

      CALL ELE(0D0,0D0,-2)
      CALL ELE(XX,YY,-3)

      DUT=0
      DO 130 JDFL=1,IDFL
      DUT=DUT+DU1(KDFG(JDFL))*DBAS(KDFL(JDFL),2)*DTX*DNX
     *       +DU2(KDFG(JDFL))*DBAS(KDFL(JDFL),2)*DTY*DNX
     *       +DU1(KDFG(JDFL))*DBAS(KDFL(JDFL),3)*DTX*DNY
     *       +DU2(KDFG(JDFL))*DBAS(KDFL(JDFL),3)*DTY*DNY
130   CONTINUE

      DFW=DFW+DLH*(DPF(1)*DUT*DNY-DPCONT*DNX)
      DAW=DAW-DLH*(DPF(1)*DUT*DNX+DPCONT*DNY)

100   CONTINUE

1000  DFW=2D0*DFW/DPF(2)      
      DAW=2D0*DAW/DPF(2) 
99999 END
      
      
      
************************************************************************
* The following routine is a modification of BDFRIG.
* It performs a line integration with a given function
* along the interface of the fictitious boundary components.
*
* This is done for testing. E.g. if the fictitious boundary is a circle,
* by inserting the function 1 the arc length of the circle should be
* calculated - if the fictitious boundary routines support exact
* reconstruction of the cubature points and weights along the
* interface.
*
* In:
*  CLDERP - Calculate function value in specific point
*  IFBC   - fictitious boundary component along which interface
*           should be integrated
*
* Out:
*  Returns the value of the integral in DVAL.
*
* ELE must be the current element, although it's not used for the
* real calculation - we perform a user-defined-only calculation.
*
* The caller has to provide the following function:
*
* SUBROUTINE CLDERP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
*                    XX,YY,XX1,YY1,ICUBP,IELEM,ELE,BNONPR)
* --------------------------------------------------------------------
* INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELEM,ICUBP
* LOGICAL BNONPR
* DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
* DOUBLE PRECISION XX,YY,XX1,YY1
* EXTERNAL ELE
*  -> Is called once for each cubature point
*  -> Must calculate the function value in the 
*     cubature point (XX,YY) on the real or (XX1,YY1) on the
*     reference element, respectively. 
*  -> ICUBP = number of current cubature point on 
*  -> IELEM = number of current element
*  -> The function has to store the function value in DINFO(1).
*  -> IINFO = number of fictitious boundary component that is
*     currently being analyzed.
*  -> DU1,DU2,DP is not used
*  -> For a description of the other parameters,
*     look at the documentation of the CVDERP function that has to
*     be provided to VOLINT or LININT.
************************************************************************
      
      SUBROUTINE UDAN03 (CLDERP,IFBC,DCORVG,KVERT,DCORMG,KMID,KNPR,NEL,
     *                   NVE,ELE,BNONPR,DVAL)
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'casmbly.inc'
      
      INCLUDE 'cmgadr.inc'

      INCLUDE 'ccub.inc'

C parameters
      DOUBLE PRECISION DCORVG(2,*),DCORMG(2,*)
      INTEGER KNPR(*),KVERT(NNVE,*),KMID(NNVE,*),IFBC,NEL,NVE
      DOUBLE PRECISION DVAL
      LOGICAL BNONPR

C externals - we partially use the routines of BDFRIG

      EXTERNAL ELE
      EXTERNAL UDAAIN,I000,CLDERP,UDAASM,CLGTLN,CLPRBD
      
C local variables

      DOUBLE PRECISION DU1,DU2,DP
      
C Structure of DINFO:
C   double function value
C   double integral value
C   double LineStart [2]
C   double LineEnd [2]
C
C Structure of LINFO
C   int     Current Fictitious boundary component

      DOUBLE PRECISION DINFO(6)
      INTEGER IINFO(1)
      
C In IINFO(1) holds the desired fictitious boundary number:

      IINFO(1) = IFBC
      
      CALL LININT (DU1,DU2,DP,DCORVG,KVERT,DCORMG,KMID,
     *    NEL,NVE,ELE,BNONPR,3,
     *    IINFO, 1, DINFO, 6,
     *    UDAAIN,I000,CLDERP,UDAASM,I000,CLGTLN,CLPRBD,.TRUE.,0,0)
C     Call not correct!!!
      
C Obtain the value:

      DVAL = DINFO(2)
      
99999 END

************************************************************************
* Line-integration callback-routines
************************************************************************

************************************************************************
* Initialisation
************************************************************************
      SUBROUTINE UDAAIN(IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,ELE,BNONPR)
      
      IMPLICIT NONE
      
C parameters
      
      INTEGER NIINFO,IINFO(NIINFO),NDINFO
      LOGICAL BNONPR
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      EXTERNAL ELE
      
C initialise values with 0
      
      DINFO(1) = 0D0
      DINFO(2) = 0D0
      
      END

************************************************************************
* Summing up the values from the cubature point to the values of
* interest.
************************************************************************
      SUBROUTINE UDAASM(IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                  DETJ,ICUBP,IELEM,OMEGA,
     *                  DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'casmbly.inc'
      
C parameters
      
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,ICUBP,IELEM
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      DOUBLE PRECISION DETJ,OMEGA
      DOUBLE PRECISION DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV
      
C Multiply the value with the weight and add it to the integral value

      DINFO(2) = DINFO(2) + DINFO(1)*OMEGA

      END

************************************************************************
* Function to integrate
*
* 1-function - should give the arc length of the fictitious boundary.
************************************************************************

      SUBROUTINE UDAAF1 (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                  XX,YY,XX1,YY1,ICUBP,IELEM,ELE,BNONPR)
      IMPLICIT NONE
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELEM,ICUBP
      LOGICAL BNONPR
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      DOUBLE PRECISION XX,YY,XX1,YY1
      EXTERNAL ELE

      DINFO(1) = 1D0
      
      END
      
************************************************************************
* Function to integrate
*
* Drag-reference-function.
* Should calculate the reference drag-force f_D.
************************************************************************

      SUBROUTINE UDAAF2 (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                   XX,YY,XX1,YY1,ICUBP,IELEM,ELE,BNONPR)
      IMPLICIT NONE
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELEM,ICUBP
      LOGICAL BNONPR
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      DOUBLE PRECISION XX,YY,XX1,YY1
      EXTERNAL ELE

C Reference function on the unit square:
      
      DOUBLE PRECISION X,Y,UF,VF,PF,UFX,VFX,PFX,UFY,VFY,PFY
     
C Reference function prescribed on the unit square, 
C together with their X- any Y-derivative
     
      UF(X,Y) = y*(1D0-x)*(0D0-x)*(1D0-y)*(0D0-y)*exp(x)
      VF(X,Y) = -x*(1D0-x)*(0D0-x)*(1D0-y)*(0D0-y)*exp(y)
      PF(X,Y) = x

      UFX(X,Y) = -y**2*x*(1-y)*exp(x)+y**2*(1-x)*(1-y)*exp(x)+
     *           y**2*(1-x)*x*(1-y)*exp(x)
      VFX(X,Y) = -2*x*(1-x)*(1-y)*y*exp(y)+x**2*(1-y)*y*exp(y)
      PFX(X,Y) = 1

      UFY(X,Y) = 2*y*(1-x)*x*(1-y)*exp(x)-y**2*(1-x)*x*exp(x)
      VFY(X,Y) = x**2*(1-x)*y*exp(y)-x**2*(1-x)*(1-y)*exp(y)-
     *           x**2*(1-x)*(1-y)*y*exp(y)
      PFY(X,Y) = 0
      
C local variables

      DOUBLE PRECISION DN1,DN2,AH1,AH2,DPF1,DPF2,DN
      DOUBLE PRECISION DT1,DT2,DUT
      DOUBLE PRECISION X1,Y1,X2,Y2
      DOUBLE PRECISION DU1V,DU2V,DPV,DU1X,DU1Y,DU2X,DU2Y,DPX,DPY
      INTEGER JDFL
      
C Calculate the function- and derivative values in the current
C cubature point:
      
      DU1V = UF(XX,YY)
      DU2V = VF(XX,YY)
      DPV = PF(XX,YY)

      DU1X = UFX(XX,YY)
      DU2X = VFX(XX,YY)
      DPX = PFX(XX,YY)

      DU1Y = UFY(XX,YY)
      DU2Y = VFY(XX,YY)
      DPY = PFY(XX,YY)
      
C Prescribe DPF1/DPF2 such that the drag coefficient c_D calculates
C to the drag force f_D.

      DPF1 = 1D0
      DPF2 = 2D0
      
C At first build the normal vector of the fictitious boundary component
C with the analytical prescription:

      CALL FBDNML (XX,YY,IINFO(1),DN1,DN2)

C Turn it around to get the (normalised) tangential vector:

      DT1 = DN2
      DT2 = -DN1
      
C Calculate the force vector. This is defined as:
C
C (FD) = int ( sigma * alpha ) dx
C (FL)    V
C
C There are now different possibilities about that sigma.
C The original and most correct formulation is the following, which
C is also known as "deformation formulation":
C
C     sigma = -p*I + dpf1*[ grad(u) + grad(u)^T ]
C
C Unfortunately this gives not the best results in the benchmark
C tests. Another formulation is the following, which is also called
C "gradient formulation":
C
C     sigma = -p*I + dpf1*[ grad(u) ]
C
C This can be used in the case that div(u)=0 - in this case
C the strong formulation of the Navier Stokes Equation, which
C can be derived using this sigma, are the same because of
C                 0 = div_h(u) = grad_h(u)^T 
C for the discretised matrices!
C
C In case of line integrals there's a third possibility how to
C define that integral. This possibility was that one which was 
C proposed in the original paper where the DFG-benchmark was proposed.
C It should be equivalent with the avove formulation, but to be
C honest, we haven't found out why...:
C
C     sigma = -p*I + dpf1 * [ Du_T/Dn_S ] * n
C           = -p*I + dpf1 * <Du*n,t> * t
C
C with n the normal vector and t the tangential vector to the surface.
C 
C Unfortunately using nonconforming finite elements the deformation
C formulation does not give very accurate results, the results with
C the gradient formulation are much better. This is not true anymore
C if other finite elements are used: Using Q2 the deformation
C formulation gives far better results!
C
C The third method was the original method used by FeatFlow for
C line integrals around a non-fictitious boundary components.
C For fictitious boundary components it should give roughly the
C same result as the deformation formulation...

C Implementation of the deformation formulation of the stress tensor:
      
      AH1 = -DPV*DN1 + DPF1*(2*DU1X*DN1+(DU2X*DN2+DU1Y*DN2))
      AH2 = -DPV*DN2 + DPF1*((DU2X*DN1+DU1Y*DN1)+2*DU2Y*DN2)
  
C Implementation of the gradient formulation of the stress tensor:

C      AH1=-DPV*DN1 + DPF1*(DU1X*DN1+DU1Y*DN2)
C      AH2=-DPV*DN2 + DPF1*(DU2X*DN1+DU2Y*DN2)
  
C Implementation of the surface integral with the help of tangential
C vectors:
  
C      DUT = (DU1X*DN1+DU1Y*DN2)*DT1 + (DU2X*DN1+DU2Y*DN2)*DT2
C      AH1 = -DPV*DN1 + DPF1*(DUT*DT1)
C      AH2 = -DPV*DN2 + DPF1*(DUT*DT2)
      
C Store the drag value for the summation

      DINFO(1)=AH1
      
      END

************************************************************************
* Function to integrate
*
* Lift-reference-function.
* Should calculate the reference drag-force f_D.
************************************************************************

      SUBROUTINE UDAAF3 (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                   XX,YY,XX1,YY1,ICUBP,IELEM,ELE,BNONPR)
      IMPLICIT NONE
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELEM,ICUBP
      LOGICAL BNONPR
      DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
      DOUBLE PRECISION XX,YY,XX1,YY1
      EXTERNAL ELE

C Reference function on the unit square:
      
      DOUBLE PRECISION X,Y,UF,VF,PF,UFX,VFX,PFX,UFY,VFY,PFY
     
C Reference function prescribed on the unit square, 
C together with their X- any Y-derivative
     
      UF(X,Y) = y*(1D0-x)*(0D0-x)*(1D0-y)*(0D0-y)*exp(x)
      VF(X,Y) = -x*(1D0-x)*(0D0-x)*(1D0-y)*(0D0-y)*exp(y)
      PF(X,Y) = x

      UFX(X,Y) = -y**2*x*(1-y)*exp(x)+y**2*(1-x)*(1-y)*exp(x)+
     *           y**2*(1-x)*x*(1-y)*exp(x)
      VFX(X,Y) = -2*x*(1-x)*(1-y)*y*exp(y)+x**2*(1-y)*y*exp(y)
      PFX(X,Y) = 1

      UFY(X,Y) = 2*y*(1-x)*x*(1-y)*exp(x)-y**2*(1-x)*x*exp(x)
      VFY(X,Y) = x**2*(1-x)*y*exp(y)-x**2*(1-x)*(1-y)*exp(y)-
     *           x**2*(1-x)*(1-y)*y*exp(y)
      PFY(X,Y) = 0
      
C local variables

      DOUBLE PRECISION DN1,DN2,AH1,AH2,DPF1,DPF2,DN
      DOUBLE PRECISION DT1,DT2,DUT
      DOUBLE PRECISION X1,Y1,X2,Y2
      DOUBLE PRECISION DU1V,DU2V,DPV,DU1X,DU1Y,DU2X,DU2Y,DPX,DPY
      INTEGER JDFL
      
C Calculate the function- and derivative values in the current
C cubature point:
      
      DU1V = UF(XX,YY)
      DU2V = VF(XX,YY)
      DPV = PF(XX,YY)

      DU1X = UFX(XX,YY)
      DU2X = VFX(XX,YY)
      DPX = PFX(XX,YY)

      DU1Y = UFY(XX,YY)
      DU2Y = VFY(XX,YY)
      DPY = PFY(XX,YY)
      
C Prescribe DPF1/DPF2 such that the drag coefficient c_D calculates
C to the drag force f_D.

      DPF1 = 1D0
      DPF2 = 2D0
      
C At first build the normal vector of the fictitious boundary component
C with the analytical prescription:

      CALL FBDNML (XX,YY,IINFO(1),DN1,DN2)

C Turn it around to get the (normalised) tangential vector:

      DT1 = DN2
      DT2 = -DN1
      
C Calculate the force vector. This is defined as:
C
C (FD) = int ( sigma * alpha ) dx
C (FL)    V
C
C There are now different possibilities about that sigma.
C The original and most correct formulation is the following, which
C is also known as "deformation formulation":
C
C     sigma = -p*I + dpf1*[ grad(u) + grad(u)^T ]
C
C Unfortunately this gives not the best results in the benchmark
C tests. Another formulation is the following, which is also called
C "gradient formulation":
C
C     sigma = -p*I + dpf1*[ grad(u) ]
C
C This can be used in the case that div(u)=0 - in this case
C the strong formulation of the Navier Stokes Equation, which
C can be derived using this sigma, are the same because of
C                 0 = div_h(u) = grad_h(u)^T 
C for the discretised matrices!
C
C In case of line integrals there's a third possibility how to
C define that integral. This possibility was that one which was 
C proposed in the original paper where the DFG-benchmark was proposed.
C It should be equivalent with the avove formulation, but to be
C honest, we haven't found out why...:
C
C     sigma = -p*I + dpf1 * [ Du_T/Dn_S ] * n
C           = -p*I + dpf1 * <Du*n,t> * t
C
C with n the normal vector and t the tangential vector to the surface.
C 
C Unfortunately using nonconforming finite elements the deformation
C formulation does not give very accurate results, the results with
C the gradient formulation are much better. This is not true anymore
C if other finite elements are used: Using Q2 the deformation
C formulation gives far better results!
C
C The third method was the original method used by FeatFlow for
C line integrals around a non-fictitious boundary components.
C For fictitious boundary components it should give roughly the
C same result as the deformation formulation...

C Implementation of the deformation formulation of the stress tensor:
      
      AH1 = -DPV*DN1 + DPF1*(2*DU1X*DN1+(DU2X*DN2+DU1Y*DN2))
      AH2 = -DPV*DN2 + DPF1*((DU2X*DN1+DU1Y*DN1)+2*DU2Y*DN2)
  
C Implementation of the gradient formulation of the stress tensor:

C      AH1=-DPV*DN1 + DPF1*(DU1X*DN1+DU1Y*DN2)
C      AH2=-DPV*DN2 + DPF1*(DU2X*DN1+DU2Y*DN2)
  
C Implementation of the surface integral with the help of tangential
C vectors:
  
C      DUT = (DU1X*DN1+DU1Y*DN2)*DT1 + (DU2X*DN1+DU2Y*DN2)*DT2
C      AH1 = -DPV*DN1 + DPF1*(DUT*DT1)
C      AH2 = -DPV*DN2 + DPF1*(DUT*DT2)
      
C Store the drag value for the summation

      DINFO(1)=AH2
      
      END

***********************************************************************
* User-defined analysis 4
*
* Calculates the nodal velocity from two grids. The routine accepts
* two deformed grids TRIA1 and TRIA2, which are snapshots at different
* points in time of the simulation. DT is the time difference between
* the snapshots. UDAN04 will subtract the node positions and divide it
* by DT to get a first order approximation of the velocity of the
* nodes.
*
* In:
*   TRIA1   - STRIA-structure for the first grid
*   TRIA2   - STRIA-structure for the second grid; must be structured
*             the same way as TRIA1.
*   DT      - double
*             Time 
*
* Out:
*   DU       - array [1..NVT] of double
*              Absolute node velocity
*   DUX      - array [1..NVT] of double
*              Node-velocity, X-component
*   DUY      - array [1..NVT] of double
*              Node-velocity, Y-component
*             
* NVT is taken from the grid-structure!
***********************************************************************

      SUBROUTINE UDAN04 (TRIA1, TRIA2, DT,
     *                   DU, DUX, DUY)
      
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      INCLUDE 'cmem.inc'
      
      DOUBLE PRECISION DT,DU(*),DUX(*),DUY(*)
      INTEGER TRIA1(SZTRIA),TRIA2(SZTRIA)
      
      INTEGER I,DCRVG1,DCRVG2
      
      DCRVG1 = L(TRIA1(OLCORVG))
      DCRVG2 = L(TRIA2(OLCORVG))
      DO I=0,TRIA1(ONVT)-1
        DUX (1+I) = (DWORK(DCRVG2+2*I)-DWORK(DCRVG1+2*I)) / DT
        DUY (1+I) = (DWORK(DCRVG2+2*I+1)-DWORK(DCRVG1+2*I+1)) / DT
        DU (1+I) = SQRT(DUX(1+I)**2+DUY(1+I)**2)
      END DO
      
      END
      
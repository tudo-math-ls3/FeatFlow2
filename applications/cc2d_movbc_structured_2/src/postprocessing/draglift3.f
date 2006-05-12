***********************************************************************
* This file provides additional drag-/lift calculation routines like 
* the file draglift2.f. These extended routines take care of the
* tourque-force on objects which allowes to calculate rotation.
***********************************************************************

************************************************************************
* Calculate lift, drag and torque 
*
* Calculates the X- and Y- body-forces as well as the torque for the 
* fictitious boundary component IFBC. If IFBC=0, the forces of all
* fictitious boundary components are calculated simultaneously.
*
* Volume integration, constant pressure
*
* In:
*  DU1/DU2/DP - Solution vector(s) for velocity and pressure
*  DPF1/DPF2 - are correction factors for the drag/lift coefficient.
*              For standard drag/lift-calculation these are provided by
*              the DPF(1)/DPF(2) parameters that are initialised in 
*              INDAT2D.F by the user initialisation routine.
*              If DPF1=1.0, DPF2=0.0, the forces are calculated.
*              Otherwise, the forces are used to calculate the
*              drag/lift coefficients
*  TRIA      - array [1..SZTRIA] of integer
*              Triangulation structure of the underlying mesh
*  DCORVG,
*  KVERT,
*  DCORMG,
*  KMID,
*  KXNPR     - Trangulation arrays; must correspond to TRIA!
*  IFBC      - Number of fictitious boundary component to calculate the
*              body forces for, or =0 to calculate all body forces.
*  ELE       - Element function
*  BNONPR    - true, if ELE is a nonparametric element
*  LAMBDA    - Relaxation parameter in [0..1] for relaxation 
*              of ALPHA vector by reconstructed interface, see RXALPH.
*              0D0=standard=no relaxation
*              1D0=treat diagonal intersection lines like horizontal
*
* Out:
*  DFRCE     - If IFBC>0: array [1..3,1..1] of double
*                DFRCE(1,1) = drag coefficient on f.b. component IFBC
*                DFRCE(2,1) = lift coefficient on f.b. component IFBC
*                DFRCE(3,1) = torque force on f.b. component IFBC
*              If IFBC=0: array [1..3,1..NFBDYC] of double
*                DFRCE(1,i) = drag coefficient on f.b. component i
*                DFRCE(2,1) = lift coefficient on f.b. component i
*                DFRCE(3,1) = torque force on f.b. component i
*              If DPF1=1.0, DPF2<>0, DFRCE(1 / 2,.) returns the drag/lift 
*              force rather than the coefficient!
************************************************************************

      SUBROUTINE BDFVLT (DU1,DU2,DP,NEQU,DCORVG,KVERT,DCORMG,
     *                   KMID,KXNPR,TRIA,
     *                   ELE,BNONPR,DFRCE,DPF1,DPF2,
     *                   IFBC,LAMBDA,IGEOM,DGEOM)
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'ccub.inc'
      
C     parameters

      DOUBLE PRECISION DU1(*),DU2(*),DP(*),DCORVG(2,*),DCORMG(2,*)
      INTEGER KXNPR(2,*),KVERT(NNVE,*),KMID(NNVE,*),IFBC
      INTEGER TRIA(SZTRIA),IGEOM(*)
      INTEGER NEQU
      DOUBLE PRECISION DFRCE(3,*),DPF1,DPF2,LAMBDA,DGEOM(*)
      LOGICAL BNONPR
      
C     externals
      
      EXTERNAL ELE
      EXTERNAL CINCP,CVDCP,CISCP,CDOCP,CVDCE
      INTEGER NDFL,NFBDYC,ISFBDY,TNBC
      EXTERNAL NDFL,NFBDYC,ISFBDY,TNBC
      
C     local variables

      INTEGER LALPHA
      INTEGER IDFL,KDFG(NNCUBP),KDFL(NNCUBP)
      DOUBLE PRECISION DJF(2,2),DCOORD(2,4),OM
      INTEGER IEDG,IPT1,IPT2,IFBED(9)
      DOUBLE PRECISION X,Y,XM,YM
      INTEGER IFIX,IAUX,I,IELTYP,ICUB,IVE,JP,IG
      DOUBLE PRECISION XMASS,YMASS,DMASS,XTORQUE,YTORQUE,ATQ
      DOUBLE PRECISION XI1,XI2,XX,YY,DN1,DN2,AH1,AH2
      DOUBLE PRECISION DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV
      DOUBLE PRECISION DAV,DAX,DAY,ALPH,DINRT
      
C     Allocate memory for the ALPHA-array

      CALL ZNEW (NEQU,1,LALPHA,'DALPHA')
      
C     Use the vector dalpha as a test function in the weak formulation
C     (only the diffusion part is taken, should take even the
C     convection but in featflow it seems to be ok) ... 
C     this needs to be checked in timedependent case

C     Prepare the ALPHA-vector for all fict. boundary components:

      CALL INALPH (DWORK(L(LALPHA)),KXNPR,KMID,
     *             TRIA(ONEL),TRIA(ONVT),TRIA(ONBCT),NEQU,IFBC)
      
C     Relaxate the alpha for nearly-horizontal/vertical lines

      CALL RXALPH (DWORK(L(LALPHA)),DCORVG,KVERT,DCORMG,KMID,
     *             TRIA(ONEL),TRIA(ONVT),
     *             NEQU,IFBC,LAMBDA,IGEOM,DGEOM)
     
C     Initialise the result-vector.
C     Depending on if we should calculate only for one particle
C     or for all, we have to fill it with 0:

      IF (IFBC.NE.0) THEN
        CALL LCL1 (DFRCE,3)
      ELSE
        IAUX = NFBDYC(IGEOM,DGEOM)
        CALL LCL1 (DFRCE,IAUX)
      END IF

C     Ok, let's start with the volume integration.      
      
      IER=0

c     which derivatives of basis functions are needed?

      DO  I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      
c     dummy call of ele sets number of element

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)

      IDFL=NDFL(IELTYP)

C     initialise cubature formula

      ICUB = 8
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
      
C     Prepare parametric element for cubature. 
C     This is done for saving arithmetic operations in later calls.
C
C     In the nonconforming case the element has to be initialised
C     separately on every element, so there's no use in initialising it
C     globally.

      IF (.NOT.BNONPR) THEN
        CALL ELE(0D0,0D0,-2)
      END IF

C     Get direct information about the object if we only have one.

      IF (IFBC.NE.0) THEN
                                                                  
        CALL FBDINF (IFBC,
     *               IFIX,XMASS,YMASS,DMASS,DINRT,
     *               IGEOM,DGEOM)
        
      END IF ! IFBC<>0

C     loop over all elements
C
C     The variable IEL is stored in the COMMON block /ELEM/ in case 
C     that any procedure has to know in which element we are...

      DO IEL=1,TRIA(ONEL)
      
C       Get the degrees of freedom on the current element

        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       Store the coordinates of the corners of the current element in 
C       the COMMON block variables DX/DY (necessary for the element) 
C       as well as in the DCOORD array for the coordinate 
C       transformation.

        DO IVE = 1, TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
          DCOORD (1,IVE) = DCORVG(1,JP)
          DCOORD (2,IVE) = DCORVG(2,JP)
        END DO

C       Prepare nonparametric elements for cubature on current element.

        IF (BNONPR) THEN
          CALL ELE(0D0,0D0,-2)
        END IF

C       Initialise auxiliary Jacobian factors DJF for transformation

        CALL QINIJF (DCOORD,DJF)
        
C       Good. But before we start to calculate the information on the
C       current element, we first have to ficure out, to which
C       fictitious boundary component the current element belongs to.
C       We test the corners and midpoints of the current element
C       to find as many fictitious boundary objects as possible, this 
C       element could give an integral value to!

        IF (IFBC.EQ.0) THEN

C         Loop over the four corners:

          DO IEDG = 1,TRIA(ONVE)
          
C           Calculate the edge midpoint

            IPT1 = KVERT(IEDG,IEL)
            
            X = DCORVG(1,IPT1)
            Y = DCORVG(2,IPT1)
            
C           Calculate in which f.b. component this edge is -
C           if at all. Set to 0 if the f.b.c. is only 
C           of virtual nature and does not describe a real object.

            IAUX = MAX(0,ISFBDY (X,Y,0,IGEOM,DGEOM)-TNBC())
            
C           Make sure we don't have that element already

            DO I=1,IEDG-1
              IF ((IAUX.NE.0).AND.(IAUX.EQ.IFBED(I)))
     *          IAUX = 0
            END DO

            IFBED(IEDG) = IAUX
          
          END DO

C         Loop over the four edges:

          XM = 0.0
          YM = 0.0

          DO IEDG = 1,TRIA(ONVE)
          
C           Calculate the edge midpoint

            IPT1 = KVERT(IEDG,IEL)
            IPT2 = KVERT(MOD(IEDG+1,TRIA(ONVE)),IEL)
            
            X = 0.5D0*(DCORVG(1,IPT1)+DCORVG(1,IPT2))
            Y = 0.5D0*(DCORVG(2,IPT1)+DCORVG(2,IPT2))
            XM = XM+X
            YM = YM+Y
            
C           Calculate in which f.b. component this edge is -
C           if at all. Set to 0 if the f.b.c. is only 
C           of virtual nature and does not describe a real object.

            IAUX = MAX(0,ISFBDY (X,Y,0,IGEOM,DGEOM)-TNBC())
            
C           Make sure we don't have that element already

            DO I=1,TRIA(ONVE)+IEDG-1
              IF ((IAUX.NE.0).AND.(IAUX.EQ.IFBED(I)))
     *          IAUX = 0
            END DO

            IFBED(TRIA(ONVE)+IEDG) = IAUX
          
          END DO
          
C         Calculate the center

          X = 0.25*XM
          Y = 0.25*YM
          
C         Calculate in which f.b. component this edge is -
C         if at all. Set to 0 if the f.b.c. is only 
C         of virtual nature and does not describe a real object.

          IAUX = MAX(0,ISFBDY (X,Y,0,IGEOM,DGEOM)-TNBC())
          
C         Make sure we don't have that element already

          DO I=1,2*TRIA(ONVE)+IEDG-1
            IF ((IAUX.NE.0).AND.(IAUX.EQ.IFBED(I)))
     *        IAUX = 0
          END DO

          IFBED(2*TRIA(ONVE)+IEDG) = IAUX

        END IF ! IFBC=0
       
C       Loop over all cubature points

        DO ICUBP = 1, NCUBP
          
c         Cubature point on the reference element

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
          
C         Calculate Jacobian matrix, determinant and mapping of the 
C         cubature point (XI1,YI1) on the reference element into 
C         (XX,YY) on the "real" element. The Jacobian matrix is stored 
C         in the COMMON block variable DJAC - necessary for the element.

          CALL QTRAF (DCOORD,DJF,DJAC,DETJ,XI1,XI2,XX,YY)

C         Calculate the weighting factor for the current cubature point
C         with the help of the Jacobian determinant

          OM = DOMEGA(ICUBP)*DETJ
          
C         Evaluate the basis functions in the cubature point
C         for the velocities

          IF(BNONPR) THEN
            CALL ELE(XX,YY,-3)
          ELSE
            CALL ELE(XI1,XI2,-3)
          ENDIF
          IF (IER.LT.0) GOTO 99999

C         Evaluate the solution values and derivatives in the 
C         cubature point:

C         X-Component:

          DU1V=0D0 ! value
          DU1X=0D0 ! x dreiv.
          DU1Y=0D0 ! y deriv
          DO I=1,IDFL
            IG=KDFG(I)
            DU1V=DU1V+DU1(IG)*DBAS(KDFL(I),1)
            DU1X=DU1X+DU1(IG)*DBAS(KDFL(I),2)
            DU1Y=DU1Y+DU1(IG)*DBAS(KDFL(I),3)
          END DO
          
C         Y-Component:
          
          DU2V=0D0 ! value
          DU2X=0D0 ! x dreiv.
          DU2Y=0D0 ! y deriv
          DO I=1,IDFL
            IG=KDFG(I)
            DU2V=DU2V+DU2(IG)*DBAS(KDFL(I),1)
            DU2X=DU2X+DU2(IG)*DBAS(KDFL(I),2)
            DU2Y=DU2Y+DU2(IG)*DBAS(KDFL(I),3)
          END DO
          
C         Pressure:
          
          DPV=DP(IEL)

C         Calculate grad(ALPHA) in the cubature point

          DAV=0D0 ! VALUE
          DAX=0D0 ! X-der.
          DAY=0D0 ! Y-der.
          DO I=1,IDFL
            IG=KDFG(I)
            ALPH=DWORK(L(LALPHA)+IG-1)
            DAV=DAV+ALPH*DBAS(KDFL(I),1)
            DAX=DAX+ALPH*DBAS(KDFL(I),2)
            DAY=DAY+ALPH*DBAS(KDFL(I),3)
          END DO

C         Form the integrand. Don't normalise the DN1/DN2, since
C         it is not really a normal vector but a scaled one, 
C         calculated with the help of ALPHA!

          DN1=-DAX
          DN2=-DAY
      
C         Calculate the force vector. This is defined as:
C
C         (FD) = int ( sigma * alpha ) dx
C         (FL)    V
C
C         There are now different possibilities about that sigma.
C         The original and most correct formulation is the following, 
C         which is also known as "deformation formulation":
C
C             sigma = -p*I + dpf1*[ grad(u) + grad(u)^T ]
C
C         Unfortunately this gives not the best results in the 
C         benchmark tests. Another formulation is the following, 
C         which is also called "gradient formulation":
C
C             sigma = -p*I + dpf1*[ grad(u) ]
C
C         This can be used in the case that div(u)=0 - in this case
C         the strong formulation of the Navier Stokes Equation, which
C         can be derived using this sigma, are the same because of
C                         0 = div(u_h) = grad(u_h)^T
C         for the discretised matrices!
C         
C         Unfortunately using nonconforming finite elements the 
C         deformation formulation does not give very accurate results, 
C         the results with the gradient formulation are much better. 
C         This is not true anymore if other finite elements are used: 
C         Using Q2 the deformation formulation gives far better results!
C
C         We give all implementations here - comment the formulation in
C         that you want to use/test!
C
C         For volume integration we use:
C
C         alpha = c*normal vector 
C         -> up to multiplicative constant, so not normalised
C         -> must not be normalised because multiplicative constant
C            is necessary for the integral because of Green's formula!
C            (see construction of this)
              
C         Implementation of the deformation formulation of the stress 
C         tensor:
          
          AH1 = -DPV*DN1 + DPF1*(2D0*DU1X*DN1+(DU2X*DN2+DU1Y*DN2))
          AH2 = -DPV*DN2 + DPF1*((DU2X*DN1+DU1Y*DN1)+2D0*DU2Y*DN2)
          
C         Implementation of the gradient formulation of the stress tensor:

C          AH1=-DPV*DN1+DPF1*(DU1X*DN1+DU1Y*DN2)
C          AH2=-DPV*DN2+DPF1*(DU2X*DN1+DU2Y*DN2)

C         Only one f.b. component?

          IF (IFBC.NE.0) THEN

C           Sum up to the integral.

C           Drag
            DFRCE(1,1) = DFRCE(1,1) + AH1*OM
            
C           Lift
            DFRCE(2,1) = DFRCE(2,1) + AH2*OM

C           Torque
            XTORQUE = XX-XMASS
            YTORQUE = YY-YMASS
            ATQ = XTORQUE*AH2-YTORQUE*AH1         

            DFRCE(3,1) = DFRCE(3,1) + ATQ*OM

          ELSE
          
C           We have to deal with all f.b. components at once.
C           IFBED is a list of up to 9 f.b. components, the current
C           element might belong to.  The calculated integral values
C           have to be added to all of these!

            DO IEDG=1,2*TRIA(ONVE)+1
            
              IAUX = IFBED(IEDG)
              IF (IAUX.NE.0) THEN
              
C               Get information about that special fictitious boundary
C               object:

                CALL FBDINF (IAUX+TNBC(),
     *                       IFIX,XMASS,YMASS,DMASS,DINRT,
     *                       IGEOM,DGEOM)
        
C               Drag
                DFRCE(1,IAUX) = DFRCE(1,IAUX) + AH1*OM
                
C               Lift
                DFRCE(2,IAUX) = DFRCE(2,IAUX) + AH2*OM

C               Torque
                XTORQUE = XX-XMASS
                YTORQUE = YY-YMASS
                ATQ = XTORQUE*AH2-YTORQUE*AH1         

                DFRCE(3,IAUX) = DFRCE(3,IAUX) + ATQ*OM
              
              END IF ! IAUX=0
            
            END DO ! IEDG
          
          END IF ! IFBC<>0

        END DO ! ICUBP
        
      END DO ! IEL

C     Correct the values according to the formulas of the
C     drag/lift coefficients with the help of the information
C     provided by PTSDAT:
      
      IF (DPF2.NE.0D0) THEN
        IF (IFBC.NE.0) THEN
          DFRCE(1,1)=2D0*DFRCE(1,1)/DPF2
          DFRCE(2,1)=2D0*DFRCE(2,1)/DPF2
        ELSE
          DO IAUX = 1,NFBDYC(IGEOM,DGEOM)
            DFRCE(1,IAUX)=2D0*DFRCE(1,IAUX)/DPF2
            DFRCE(2,IAUX)=2D0*DFRCE(2,IAUX)/DPF2
          END DO
        END IF
      END IF      

C     Release memory, finish

      CALL ZDISP (0,LALPHA,'DALPHA')
      
99999 END

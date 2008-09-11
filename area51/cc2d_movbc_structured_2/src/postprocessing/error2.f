************************************************************************
* The routines in this file realize error analyzing routines
* for velocity and pressure components.
************************************************************************

************************************************************************
* Error analysis for velocity components
*
* Parametric version; for element E030,E031
*
* Computes L2-error, H1-error, divergence,... of the velocity vectors
*
* In:
*   DU1    - array [1..*] of double
*            Velocity X-component
*   DU2    - array [1..*] of double
*            Velocity Y-component
*   TRIA   - Current triangulation
*   KVERT,
*   KMID,
*   DCORVG - Usual geometry information; must correspond to TRIA
*   ELE    - SUBROUTINE; Element routine
*   ICUB   - Cubature formula for error analysis
*   UE     - SUBROUTINE; Reference solution
*   UX     - SUBROUTINE; Reference X-derivative
*   UY     - SUBROUTINE; Reference Y-derivative
*   TIMENS - Current simulation time
*   RE     - Reynolds-number of the simulation
*   IPARAM - array [1..*] of integer
*   DPARAM - array [1..*] of double
*            IPARAM/DPARAM are passed to UE,UX,UY as parameter blocks
*            of the solver.
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to boundary
*            routines. Not used in this routine.
*
* Out:
*   ERRL2  - L2-error 
*   ERRH1  - H1-error
*   DIVL2  - L2-norm of divergence
*   DIVLI  - SUP-norm of divergence
*   DNL2   - L2-norm of reference solution
*   DNH1   - H1-norm of reference solution
************************************************************************

      SUBROUTINE ELPQU2(DU1,DU2,TRIA,KVERT,KMID,DCORVG,ELE,ICUB,
     *                  UE,UX,UY,TIMENS,RE,
     *                  ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1,
     *                  IPARAM,DPARAM,IGEOM,DGEOM)

      IMPLICIT NONE

C main COMMON blocks

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'

      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'

C parameters
      
      DOUBLE PRECISION DU1(*), DU2(*),DCORVG(2,*),UE,UX,UY
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),TRIA(SZTRIA),IGEOM(*)
      INTEGER ICUB
      DOUBLE PRECISION TIMENS,RE,DGEOM(*)
      
      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      EXTERNAL ELE
      
C externals

      INTEGER NDFL
      EXTERNAL NDFL
      
C local variables      

      DOUBLE PRECISION ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1
      INTEGER IVE,JP
      DOUBLE PRECISION XX,YY,OM
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4,XI1,XI2
      DOUBLE PRECISION UH1,UH2,UH1X,UH2X,UH1Y,UH2Y
      INTEGER IELTYP,IDER,JDOFE,IEQ,ILO
      DOUBLE PRECISION UVAL1,UVAL2,UVALX1,UVALX2,UVALY1,UVALY2
      
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL
      
      
      SUB='ELPQU'
      IER=0

C     Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IF (IER.NE.0) GOTO 99999
      
C     Get local degrees of freedom
      
      IDFL=NDFL(IELTYP)
      IF (IER.LT.0) GOTO 99999
      
C     Initialize cubature formula
      
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
      
C     Initialize derivative flags; we want to have
C     function value and 1st derivatives:
      
      DO IDER=1,NNDER
        BDER(IDER)=.FALSE.
      END DO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.

      ERRL2= 0D0
      ERRH1= 0D0
      DIVL2= 0D0
      DIVLI=-1D0
      DNL2=  0D0
      DNH1=  0D0
      
C     Call the element top precalculate information in the cubature 
C     points:
      
      CALL ELE(0D0,0D0,-2)

C     Loop over the elements to calculate the information:

      DO IEL=1,TRIA(ONEL)

C       Get the global DOF's on our element

        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       Initialize the element parameters to match to the current element

        DO IVE=1,TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

C       Initialize auxiliary Jacobian factors for the transformaion

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C       Loop over the cubature points to calculate the values there:

        DO ICUBP=1,NCUBP

C         Coordinates of the cubature point on the reference element

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
      
C         Calculate Jacobian of the mapping and Jacobian determinant
            	
          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)

C         And the weight in the cuibature formula
	
          OM=DOMEGA(ICUBP)*DETJ

C         Get the value of the element in the cubature point:

          CALL ELE(XI1,XI2,-3)
          IF (IER.LT.0) GOTO 99999

C         Map the cubature point to the real element
  
          XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *      +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
          YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *      (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2   

C         Calculate value, X-derivative and Y-derivative of
C         DU1 and DU2 in our cubature point by a loop over the
C         DOF's in our element

          UH1 =0D0
          UH2 =0D0
          UH1X=0D0
          UH2X=0D0
          UH1Y=0D0
          UH2Y=0D0
          
          DO JDOFE=1,IDFL
            IEQ=KDFG(JDOFE)
            ILO=KDFL(JDOFE)
            UH1 =UH1 +DU1(IEQ)*DBAS(ILO,1)
            UH2 =UH2 +DU2(IEQ)*DBAS(ILO,1)
            UH1X=UH1X+DU1(IEQ)*DBAS(ILO,2)
            UH2X=UH2X+DU2(IEQ)*DBAS(ILO,2)
            UH1Y=UH1Y+DU1(IEQ)*DBAS(ILO,3)
            UH2Y=UH2Y+DU2(IEQ)*DBAS(ILO,3)
          END DO

C         Subtract calculated value from real solution...
C         Compute L2-error,

          UVAL1 = UE(XX,YY,1,TIMENS,RE,0,TRIA,IPARAM,DPARAM,
     *               IGEOM,DGEOM)
          UVAL2 = UE(XX,YY,2,TIMENS,RE,0,TRIA,IPARAM,DPARAM,
     *               IGEOM,DGEOM)

          ERRL2=ERRL2+OM*((UVAL1-UH1 )**2+(UVAL2-UH2 )**2)
          
C         H1-error,
          
          UVALX1 = UX(XX,YY,1,TIMENS,RE,0,TRIA,IPARAM,DPARAM,
     *               IGEOM,DGEOM)
          UVALX2 = UX(XX,YY,2,TIMENS,RE,0,TRIA,IPARAM,DPARAM,
     *               IGEOM,DGEOM)
          UVALY1 = UY(XX,YY,1,TIMENS,RE,0,TRIA,IPARAM,DPARAM,
     *               IGEOM,DGEOM)
          UVALY2 = UY(XX,YY,2,TIMENS,RE,0,TRIA,IPARAM,DPARAM,
     *               IGEOM,DGEOM)

          ERRH1=ERRH1+OM*((UVALX1-UH1X)**2+(UVALY1-UH1Y)**2)
          ERRH1=ERRH1+OM*((UVALX2-UH2X)**2+(UVALY2-UH2Y)**2)
          
C         divergence in L2-norm
          
          DIVL2=DIVL2+OM*(UH1X+UH2Y)**2
          
C         divergence in SUP-norm,          
          
          DIVLI=MAX(DIVLI,ABS(UH1X+UH2Y))
          
C         as well as L2-norm and H1-norm of the reference solution.

          DNL2=DNL2+OM*(UVAL1**2+UVAL2**2)
          DNH1=DNH1+OM*(UVALX1**2+UVALY1**2)
          DNH1=DNH1+OM*(UVALX2**2+UVALY2**2)

        END DO ! ICUBP
        
      END DO ! IEL

99999 END

************************************************************************
* Error analysis for pressure components
*
* Parametric version; for E010
*
* Computes L2-error, and prints the results onto screen.
*
* In:
*   P      - array [1..*] of double
*            Pressure solution
*   TRIA   - Current triangulation
*   KVERT,
*   KMID,
*   DCORVG - Usual geometry information; must correspond to TRIA
*   ELE    - SUBROUTINE; Element routine
*   ICUB   - Cubature formula for error analysis
*   PE     - SUBROUTINE; Reference solution
*   MFILE  - Handle of file where to write the output to, additionally
*            to the screen
*   TIMENS - Current simulation time
*   RE     - Reynolds-number of the simulation
*   IPARAM - array [1..*] of integer
*   DPARAM - array [1..*] of double
*            IPARAM/DPARAM are passed to PE as parameter blocks
*            of the solver.
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to boundary
*            routines. Not used in this routine.
*
* Out:
*   ERRL2  - L2-error 
*   DNL2   - L2-norm of reference solution
*   PERR   - array [1..*] of double
*     This vector corresponds to the pressure vector P. For each DOF
*     in P, the corresponding entry in this vector receives the
*     absolute error to the reference pressure in this DOF.
************************************************************************

      SUBROUTINE ELPQP2(P,PERR,TRIA,KVERT,KMID,DCORVG,ELE,ICUB,PE,
     *                  TIMENS,RE,ERRL2,DNL2,IPARAM,DPARAM,IGEOM,DGEOM)

      IMPLICIT NONE

C main COMMON blocks

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'

      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'

C parameters
      
      DOUBLE PRECISION P(*),PERR(*),DGEOM(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),TRIA(SZTRIA),IGEOM(*)
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER ICUB
      DOUBLE PRECISION PE
      DOUBLE PRECISION ERRL2,DNL2,TIMENS,RE

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      EXTERNAL ELE
      
C externals

      INTEGER NDFL
      EXTERNAL NDFL

C local variables      

      INTEGER IVE,JP
      DOUBLE PRECISION XX,YY,OM,PH
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4,XI1,XI2
      INTEGER IELTYP,IDER,JDOFE,IEQ,ILO
      DOUBLE PRECISION PVAL
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL

      SUB='ELPQP'
      IER=0

C     Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IF (IER.NE.0) GOTO 99999
      
C     Get local degrees of freedom
      
      IDFL=NDFL(IELTYP)
      IF (IER.LT.0) GOTO 99999
      
C     Initialize cubature formula
      
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
      
C     Initialize derivative flags; we only need function values.
      
      DO IDER=1,NNDER
        BDER(IDER)=.FALSE.
      END DO
      
      BDER(1)=.TRUE.

      ERRL2=0D0
      DNL2=0D0
      
C     Call the element top precalculate information in the cubature 
C     points:

      CALL ELE(0D0,0D0,-2)

C     Loop over the elements to calculate the information:

      DO IEL=1,TRIA(ONEL)

C       Get the global DOF's on our element

        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       Initialize the element parameters to match to the current element

        DO IVE=1,TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

C       Initialize auxiliary Jacobian factors for the transformaion

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C       Loop over the cubature points to calculate the values there:

        DO ICUBP=1,NCUBP

C         Coordinates of the cubature point on the reference element

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
      
C         Calculate Jacobian of the mapping and Jacobian determinant
            	
          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)

C         And the weight in the cuibature formula
	
          OM=DOMEGA(ICUBP)*DETJ

C         Get the value of the element in the cubature point:

          CALL ELE(XI1,XI2,-3)
          IF (IER.LT.0) GOTO 99999

C         Map the cubature point to the real element
  
          XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *      +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
          YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *      (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2   

C         Calculate the pressure in the cubature point:

          PH=0D0
          DO JDOFE=1,IDFL
            IEQ=KDFG(JDOFE)
            ILO=KDFL(JDOFE)
            PH =PH +P(IEQ)*DBAS(ILO,1)
          END DO

C         Subtract calculated value from real solution.
C         Compute L2-error,

          PVAL = PE(XX,YY,0,TIMENS,RE,
     *              0,TRIA,IPARAM,DPARAM,
     *              IGEOM,DGEOM)
      
          ERRL2=ERRL2+OM*(PVAL-PH)**2
          
C         as well as L2-norm of the reference solution.

          DNL2 =DNL2 +OM* PVAL**2
          
C         Calculate the pointwise error of each pressure DOF.
C         Add it equally distributed to all DOF's on the element,
C         as for an arbitrary element, we don't know the coordinates
C         of the DOF's!
          
          DO JDOFE=1,IDFL
            IEQ=KDFG(JDOFE)
            PERR(IEQ) = ABS(PVAL-PH) / DBLE(IDFL)
          END DO
        
        END DO  
        
      END DO
 
99999 END

************************************************************************
* Error analysis for velocity components
*
* Nonparametric version; for element EM30,EM31
*
* Computes L2-error, H1-error, divergence,... and prints the results
* onto screen.
*
* In:
*   DU1    - array [1..NMT] of double
*            Velocity X-component
*   DU2    - array [1..NMT] of double
*            Velocity Y-component
*   TRIA   - Current triangulation
*   KVERT,
*   KMID,
*   DCORVG - Usual geometry information; must correspond to TRIA
*   ELE    - SUBROUTINE; Element routine
*   ICUB   - Cubature formula for error analysis
*   UE     - SUBROUTINE; Reference solution
*   UX     - SUBROUTINE; Reference X-derivative
*   UY     - SUBROUTINE; Reference Y-derivative
*   MFILE  - Handle of file where to write the output to, additionally
*            to the screen
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to boundary
*            routines. Not used in this routine.
*
* Out:
*   ERRL2  - L2-error 
*   ERRH1  - H1-error
*   DIVL2  - L2-norm of divergence
*   DIVLI  - SUP-norm of divergence
*   DNL2   - L2-norm of reference solution
*   DNH1   - H1-norm of reference solution
************************************************************************

      SUBROUTINE ELPQN2(DU1,DU2,TRIA,KVERT,KMID,DCORVG,ELE,ICUB,
     *                  UE,UX,UY,TIMENS,RE,
     *                  ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1,
     *                  IPARAM,DPARAM,IGEOM,DGEOM)

      IMPLICIT NONE

C main COMMON blocks

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'

      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'

C parameters
      
      DOUBLE PRECISION DU1(*), DU2(*),DGEOM(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),TRIA(SZTRIA),IGEOM(*)
      DOUBLE PRECISION DCORVG(2,*),TIMENS,RE
      INTEGER ICUB
      DOUBLE PRECISION UE,UX,UY
      DOUBLE PRECISION ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      EXTERNAL ELE
      
C externals

      INTEGER NDFL
      EXTERNAL NDFL

C local variables      

      DOUBLE PRECISION UH1,UH2,UH1X,UH2X,UH1Y,UH2Y
      INTEGER IVE,JP
      DOUBLE PRECISION XX,YY,OM
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4,XI1,XI2
      INTEGER IELTYP,IDER,JDOFE,IEQ,ILO
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL
      
      DOUBLE PRECISION UVAL1,UVAL2,UVALX1,UVALX2,UVALY1,UVALY2

      SUB='ELPQN'
      IER=0
C     Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IF (IER.NE.0) GOTO 99999
      
C     Get local degrees of freedom
      
      IDFL=NDFL(IELTYP)
      IF (IER.LT.0) GOTO 99999
      
C     Initialize cubature formula
      
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
      
C     Initialize derivative flags; we want to have
C     function value and 1st derivatives:
      
      DO IDER=1,NNDER
        BDER(IDER)=.FALSE.
      END DO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.

      ERRL2= 0D0
      ERRH1= 0D0
      DIVL2= 0D0
      DIVLI=-1D0
      DNL2=  0D0
      DNH1=  0D0
      
C     Loop over the elements to calculate the information:

      DO IEL=1,TRIA(ONEL)

C       Get the global DOF's on our element

        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       Initialize the element parameters to match to the current element

        DO IVE=1,TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

C       Initialize auxiliary Jacobian factors for the transformaion

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C       Call the element to precalculate values in the cubature
C       points of our element:
      
        CALL ELE(0D0,0D0,-2)

C       Loop over the cubature points to calculate the values there:

        DO ICUBP=1,NCUBP

C         Coordinates of the cubature point on the reference element

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
      
C         Calculate Jacobian of the mapping and Jacobian determinant
            	
          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)

C         And the weight in the cuibature formula
	
          OM=DOMEGA(ICUBP)*DETJ

C         Map the cubature point to the real element
  
          XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *      +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
          YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *      (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2   

C         Get the value of the element in the cubature point:

          CALL ELE(XX,YY,-3)
          IF (IER.LT.0) GOTO 99999

C         Calculate value, X-derivative and Y-derivative of
C         DU1 and DU2 in our cubature point by a loop over the
C         DOF's in our element

          UH1 =0D0
          UH2 =0D0
          UH1X=0D0
          UH2X=0D0
          UH1Y=0D0
          UH2Y=0D0
          
          DO JDOFE=1,IDFL
            IEQ=KDFG(JDOFE)
            ILO=KDFL(JDOFE)
            UH1 =UH1 +DU1(IEQ)*DBAS(ILO,1)
            UH2 =UH2 +DU2(IEQ)*DBAS(ILO,1)
            UH1X=UH1X+DU1(IEQ)*DBAS(ILO,2)
            UH2X=UH2X+DU2(IEQ)*DBAS(ILO,2)
            UH1Y=UH1Y+DU1(IEQ)*DBAS(ILO,3)
            UH2Y=UH2Y+DU2(IEQ)*DBAS(ILO,3)
          END DO

C         Subtract calculated value from real solution.
C         Compute L2-error,

          UVAL1 = UE(XX,YY,1,TIMENS,RE,0,TRIA,IPARAM,DPARAM,
     *               IGEOM,DGEOM)
          UVAL2 = UE(XX,YY,2,TIMENS,RE,0,TRIA,IPARAM,DPARAM,
     *               IGEOM,DGEOM)

          ERRL2=ERRL2+OM*((UVAL1-UH1 )**2+(UVAL2-UH2 )**2)
          
C         H1-error,
          
          UVALX1 = UX(XX,YY,1,TIMENS,RE,0,TRIA,IPARAM,DPARAM,
     *               IGEOM,DGEOM)
          UVALX2 = UX(XX,YY,2,TIMENS,RE,0,TRIA,IPARAM,DPARAM,
     *               IGEOM,DGEOM)
          UVALY1 = UY(XX,YY,1,TIMENS,RE,0,TRIA,IPARAM,DPARAM,
     *               IGEOM,DGEOM)
          UVALY2 = UY(XX,YY,2,TIMENS,RE,0,TRIA,IPARAM,DPARAM,
     *               IGEOM,DGEOM)

          ERRH1=ERRH1+OM*((UVALX1-UH1X)**2+(UVALY1-UH1Y)**2)
          ERRH1=ERRH1+OM*((UVALX2-UH2X)**2+(UVALY2-UH2Y)**2)
          
C         divergence in L2-norm
          
          DIVL2=DIVL2+OM*(UH1X+UH2Y)**2
          
C         divergence in SUP-norm,          
          
          DIVLI=MAX(DIVLI,ABS(UH1X+UH2Y))
          
C         as well as L2-norm and H1-norm of the reference solution.

          DNL2=DNL2+OM*(UVAL1**2+UVAL2**2)
          DNH1=DNH1+OM*(UVALX1**2+UVALY1**2)
          DNH1=DNH1+OM*(UVALX2**2+UVALY2**2)

        END DO ! ICUBP
        
      END DO ! IEL

99999 END

************************************************************************
* Error calculation and output routine
*
* This routine performs basic error analysis of the solution vector
* and prints the results onto screen / into file.
*
* In:
*   MFILE  : handle to an file where output is written to
*   IERANA : Cubature formula to use for error analysis
*            =0: don't perform error analysis
*   VECDAT : array [1..SZN2VI] of integer 
*            TNS2DVectorParams-structure, c
*            Defines the form of the RHS vector.
*   NLMIN  : minimum level 
*   NLMAX  : maximum level 
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for all levels.
*   DUP    : array [1..*] of double
*            Current solution vector - on level NLMAX
*   DRHS   : array [1..*] of double
*            RHS vector
*   DAUX   : array [1..*] of double
*            Auxiliary vector
*
*   TIMENS : Current simulation time; can be 0 for stationary
*            simulation
*
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. Are passed to the solution routines UE, UEX,
*            UEY,...
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to boundary
*            routines. Not used in this routine.
*
*   UE,
*   UEX,
*   UEY    : Callback functions; exact (analytical) velocity and
*            derivatives for error analysis
*   PE     : Callback functions; exact (analytical) pressore for error 
*            analysis
************************************************************************
      
      SUBROUTINE ERATRM (MFILE,IERANA,VECDAT,DUP,DRHS,DAUX,
     *                   NLMIN,NLMAX,TRIAS,TIMENS,
     *                   IASMBL,DASMBL,IGEOM,DGEOM,
     *                   UE,PE,UEX,UEY)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cout.inc'
      
      INCLUDE 'stria.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'sassembly.inc'

      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
C     parameters

      INTEGER MFILE,IERANA,IGEOM(*),VECDAT(SZN2VI)
      DOUBLE PRECISION DUP(*),DRHS(*),DAUX(*),DGEOM(*)
      INTEGER IASMBL(*),NLMIN,NLMAX,TRIAS(SZTRIA,NNLEV)
      DOUBLE PRECISION DASMBL(*),TIMENS
      
C     local variables
      
      INTEGER ICUBER,IELT,LVERT,LMID,LCORVG,IELTYP,NEQU,NEQP
      DOUBLE PRECISION ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1
      DOUBLE PRECISION ERRL2P,DNL2P,RE
      
C     Coefficient of exact solution

      DOUBLE PRECISION UE,PE,UEX,UEY
      EXTERNAL UE,PE,UEX,UEY

C     definition of finite elements
      EXTERNAL E030,E031,EM30,EM31,E010
      
C     Activate error analysis at all?
      
      IF (IERANA.LE.0) RETURN
      
C     Get the vector size of velocity and pressure part:

      NEQU = VECDAT(ONU)
      NEQP = VECDAT(ONP)

C     Use cubature formula IERANA for error analysis
      
      ICUBER=IERANA
      
      IELT = IASMBL(OIELEMT)
      LVERT = TRIAS(OLVERT,NLMAX)
      LMID = TRIAS(OLMID,NLMAX)
      LCORVG = TRIAS(OLCORVG,NLMAX)
      
      RE = DASMBL(ORE)
      
      IF (IELT.EQ.0) THEN
        CALL ELPQU2(DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),
     *              KWORK(L(LVERT)),KWORK(L(LMID)),
     *              DWORK(L(LCORVG)),E031,ICUBER,UE,UEX,UEY,
     *              TIMENS,RE,
     *              ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1,
     *              IASMBL,DASMBL,IGEOM,DGEOM)
        IELTYP = 31
      END IF
     
      IF (IELT.EQ.1) THEN
        CALL ELPQU2(DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),
     *              KWORK(L(LVERT)),KWORK(L(LMID)),
     *              DWORK(L(LCORVG)),E030,ICUBER,UE,UEX,UEY,
     *              TIMENS,RE,
     *              ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1,
     *              IASMBL,DASMBL,IGEOM,DGEOM)
        IELTYP = 30
      END IF
     
      IF (IELT.EQ.2) THEN
        CALL ELPQN2(DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),
     *              KWORK(L(LVERT)),KWORK(L(LMID)),
     *              DWORK(L(LCORVG)),EM31,ICUBER,UE,UEX,UEY,
     *              TIMENS,RE,
     *              ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1,
     *              IASMBL,DASMBL,IGEOM,DGEOM)
        IELTYP = -31
      END IF
     
      IF (IELT.EQ.3) THEN
        CALL ELPQN2(DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),
     *              KWORK(L(LVERT)),KWORK(L(LMID)),
     *              DWORK(L(LCORVG)),EM30,ICUBER,UE,UEX,UEY,
     *              TIMENS,RE,
     *              ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1,
     *              IASMBL,DASMBL,IGEOM,DGEOM)
        IELTYP = -30
      END IF

C     Write the absolute error to the reference solution in the
C     pressure part into the auxiliary vector

      CALL ELPQP2(DUP(1+2*NEQU),DAUX(1+2*NEQU),TRIAS(1,NLMAX),
     *           KWORK(L(LVERT)),
     *           KWORK(L(LMID)),DWORK(L(LCORVG)),E010,ICUBER,PE,
     *           TIMENS,RE,ERRL2P,DNL2P,IASMBL,DASMBL,IGEOM,DGEOM)

C     Print out the results.
C
C     Velocity:

      IF (DNL2.LT.1D-15) THEN
        WRITE(MTERM,*)
        WRITE(MTERM,*) '* ELPU * EXACT SOLUTION ZERO !!!'
        WRITE(MFILE,*)
        WRITE(MFILE,*) '* ELPU * EXACT SOLUTION ZERO !!!'
        DNL2=1D0
      ENDIF

      IF (DNH1.LT.1D-15) THEN
        WRITE(MTERM,*)
        WRITE(MTERM,*) '* ELPU * EXACT DERIVATIVE ZERO !!!'
        WRITE(MFILE,*)
        WRITE(MFILE,*) '* ELPU * EXACT DERIVATIVE ZERO !!!'
        DNH1=1D0
      ENDIF

      WRITE(MTERM,*)
      WRITE(MTERM,*) '*ELPU*REL. L2-H1-ERROR',
     *                 SQRT(ERRL2/DNL2),SQRT(ERRH1/DNH1),IELTYP,ICUBER
      WRITE(MTERM,*) '*ELPQU*DIVERGENCE      ',SQRT(DIVL2),DIVLI

      IF (MTERM.NE.MFILE) THEN
        WRITE(MFILE,*)
        WRITE(MFILE,*) '*ELPU*REL. L2-H1-ERROR',
     *                 SQRT(ERRL2/DNL2),SQRT(ERRH1/DNH1),IELTYP,ICUBER
        WRITE(MFILE,*) '*ELPQU*DIVERGENCE      ',SQRT(DIVL2),DIVLI
      ENDIF

C     Additionally to the error analysis with cubature formula
C     IERANA, perform the error analysis also with cubature formula
C     1 = constant interpolation for reference purposes:

      IF (IERANA.NE.1) THEN
      
        ICUBER=1
        
        IF (IELT.EQ.0) 
     *    CALL ELPQU2(DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),KWORK(L(LVERT)),
     *                KWORK(L(LMID)),DWORK(L(LCORVG)),E031,ICUBER,
     *                UE,UEX,UEY,TIMENS,RE,
     *                ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1,
     *                IASMBL,DASMBL,IGEOM,DGEOM)
     
        IF (IELT.EQ.1) 
     *    CALL ELPQU2(DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),KWORK(L(LVERT)),
     *                KWORK(L(LMID)),DWORK(L(LCORVG)),E030,ICUBER,
     *                UE,UEX,UEY,TIMENS,RE,
     *                ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1,
     *                IASMBL,DASMBL,IGEOM,DGEOM)
     
        IF (IELT.EQ.2) 
     *    CALL ELPQN2(DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),KWORK(L(LVERT)),
     *                KWORK(L(LMID)),DWORK(L(LCORVG)),EM31,ICUBER,
     *                UE,UEX,UEY,TIMENS,RE,
     *                ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1,
     *                IASMBL,DASMBL,IGEOM,DGEOM)
     
        IF (IELT.EQ.3) 
     *    CALL ELPQN2(DUP(1),DUP(1+NEQU),TRIAS(1,NLMAX),KWORK(L(LVERT)),
     *                KWORK(L(LMID)),DWORK(L(LCORVG)),EM30,ICUBER,
     *                UE,UEX,UEY,TIMENS,RE,
     *                ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1,
     *                IASMBL,DASMBL,IGEOM,DGEOM)

C       Write the absolute error to the reference solution in the
C       pressure part into the auxiliary vector

        CALL ELPQP2(DUP(1+2*NEQU),DAUX(1+2*NEQU),TRIAS(1,NLMAX),
     *              KWORK(L(LVERT)),
     *              KWORK(L(LMID)),DWORK(L(LCORVG)),E010,ICUBER,PE,
     *              TIMENS,RE,ERRL2P,DNL2P,IASMBL,DASMBL,IGEOM,DGEOM)

C       Print out the results:

        IF (DNL2.LT.1D-15) THEN
          WRITE(MTERM,*)
          WRITE(MTERM,*) '* ELPU * EXACT SOLUTION ZERO !!!'
          WRITE(MFILE,*)
          WRITE(MFILE,*) '* ELPU * EXACT SOLUTION ZERO !!!'
          DNL2=1D0
        ENDIF

        IF (DNH1.LT.1D-15) THEN
          WRITE(MTERM,*)
          WRITE(MTERM,*) '* ELPU * EXACT DERIVATIVE ZERO !!!'
          WRITE(MFILE,*)
          WRITE(MFILE,*) '* ELPU * EXACT DERIVATIVE ZERO !!!'
          DNH1=1D0
        ENDIF

        WRITE(MTERM,*)
        WRITE(MTERM,*) '*ELPU*REL. L2-H1-ERROR',
     *                   SQRT(ERRL2/DNL2),SQRT(ERRH1/DNH1),IELTYP,ICUBER
        WRITE(MTERM,*) '*ELPQU*DIVERGENCE      ',SQRT(DIVL2),DIVLI

        IF (MTERM.NE.MFILE) THEN
          WRITE(MFILE,*)
          WRITE(MFILE,*) '*ELPU*REL. L2-H1-ERROR',
     *                   SQRT(ERRL2/DNL2),SQRT(ERRH1/DNH1),IELTYP,ICUBER
          WRITE(MFILE,*) '*ELPQU*DIVERGENCE      ',SQRT(DIVL2),DIVLI
        ENDIF
        
      ENDIF ! (IERANA.NE.1)
        
      END
      
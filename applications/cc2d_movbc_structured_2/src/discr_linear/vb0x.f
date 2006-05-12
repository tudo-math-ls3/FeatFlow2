************************************************************************
* The routines in this file are extended routines to build a right-
* hand side vector on quadrilateral meshes. They accept the
* extended calling convention, thus passing a user-defined parameter
* block to callback routines.
*
* All routines base on XBV0/VB0-routines of the FEAT library.
************************************************************************

************************************************************************
* Generate RHS-vector, wrapper-routine
*
* Combined version for parametric version and nonparametric elements
*
* This function calls VB0X to generate a RHS-vector, based on a
* user-defined bilinear form. LB is a handle to the vector on the
* heap, that should be overwritten with data. If LB=0, a new handle
* will be allocated and returned in LB.
*
* In:
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation information to use for building the RHS
*   IPARAM - array [1..*] of integer
*            User-defined integer array to pass to the
*            callback-routines
*   DPARAM - array [1..*] of double
*            User-defined double precision array to pass to the
*            callback-routines
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to callback
*            routines. Not used in this routine.
*   NBLOC  - Number of matrix blocks
*   LB     - array [1..NBLOC] pf integer
*            For every matrix block:
*            Handle to an existing RHS-vector on the heap, or
*            =0, if a new vector should be allocated.
*   NEQ    - integer
*            Number of entries in RHS-vector
*   ICLEAR - Clear vectors before calculation; only valid if LB<>0.
*            =0: don't clear vector; add new RHS-vector directly
*                to existing RHS-vector
*            =1: fill RHS-vector with 0 before calculation, thus
*                calculating a complete new RHS-vector
*   ELE    - Element to use for calculation
*   BNONPR - Whether ELE is a nonparametric element
*   COEFF  - DOUBLE PRECISION FUNCTION (X,Y,IA,IBLOC,BFIRST,
*                                       TRIA,IEL,IPARAM,DPARAM,
*                                       IGEOM,DGEOM)
*            Function to determine the value of the RHS at a given
*            point.
*            In: 
*              X,Y    - The point where to evaluate the RHS
*              IA     - descriptor of the current summand in the linear
*                       form
*              IBLOC  - Number of current matrix block; =-1 for dummy-
*                       calls
*              BFIRST - =true on first fall, =false on later calls.
*                        allowes to precalculate information
*              IEL    - Number of current element where the cubature
*                       takes place; =0 for dummy-calls
*              TRIA   - Structure of the current triangulation
*              IPARAM,
*              DPARAM - the user defined arrays from above
*              IGEOM,
*              DGEOM  - the user defined arrays from above
*   KBN     - array [1..NBLOC]
*             For every matrix block: Number of summands in the linear
*             form
*   KB      - array [1..KBN,1..NBLOC] of integer
*             Descriptors of the summands in the linear form
*   BCON    - array [1..NBLOC] of boolean
*             For each matrix block:
*             =true, if the coefficients in the integral are constant
*   ICUB    - Cubature formula to use for integration on quadrilaterals
*   ARR     - array [1..NBLOC] of CHARACTER*(6)
*             Name of RHS-vector; for debug output
************************************************************************

      SUBROUTINE XVB0X(TRIA, IPARAM, DPARAM,IGEOM,DGEOM,
     *                 NBLOC,LB,NEQ,ICLEAR,ELE,BNONPR,
     *                 COEFF,KBN,KB,BCON,ICUB,ARR)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'stria.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'ccub.inc'
      
      INTEGER TRIA(SZTRIA),IPARAM(*),NEQ,ICLEAR,ICUB,IGEOM(*)
      DOUBLE PRECISION DPARAM(*),DGEOM(*)
      INTEGER NBLOC,KBN(NBLOC),KB(NNAB,NBLOC),LB(NBLOC)
      LOGICAL BCON(NBLOC),BNONPR
      CHARACTER ARR(NBLOC)*(6)

      DOUBLE PRECISION COEFF
      EXTERNAL ELE,COEFF

      INTEGER IBLOC,ITYPE,ILEN,LOFF,LOECON

      IER=0

      DO IBLOC=1,NBLOC
      
C       Allocate new handles for the matrix blocks if necessary
      
        IF (LB(IBLOC).EQ.0) THEN
        
          CALL ZNEW(NEQ,1,LB(IBLOC),ARR(IBLOC))
          IF (IER.NE.0) RETURN
          
        ELSE

C        Check input parameters: Type of input array, length,...

          CALL ZTYPE(LB(IBLOC),ITYPE)
          IF (ITYPE.NE.1) THEN
            WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
            CALL WERR(-114,'XVB0  ')
            RETURN
          ENDIF
          CALL ZLEN(LB(IBLOC),ILEN)
          IF (ILEN.LT.NEQ) THEN
            WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
            CALL WERR(-115,'XVB0  ')
            RETURN
          ENDIF
          
C         Clear RHS-vector if desired
          
          IF (ICLEAR.EQ.1) CALL ZCLEAR(LB(IBLOC),ARR(IBLOC))
          
        END IF
      END DO

C     Allocate temporary storage for (perhaps) constant coefficients,
C     offset-array of RHS-vectors in DWORK

      CALL ZNEW(NBLOC,-3,LOFF,'KOFF  ')
      IF (IER.NE.0) RETURN
      
      CALL ZNEW(NBLOC*NNDER,1,LOECON,'COECON')
      IF (IER.NE.0) RETURN

C     Initialize KOFF-array with offsets of RHS-vectors relative
C     to DWORK(1)

      DO IBLOC=1,NBLOC
        KWORK(L(LOFF)+IBLOC-1)=L(LB(IBLOC))-1
      END DO

C     Call the RHS vector generation routine

      CALL VB0X(TRIA,IPARAM,DPARAM,IGEOM,DGEOM,
     *          NBLOC,DWORK(1),KWORK(L(LOFF)),
     *          KWORK(L( TRIA(OLVERT) )), KWORK(L( TRIA(OLMID) )),
     *          DWORK(L( TRIA(OLCORVG) )),
     *          ELE,BNONPR,COEFF,
     *          KBN,KB,BCON,DWORK(L(LOECON)),ICUB)
      IF (IER.NE.0) RETURN
      
C     Release temporary storage
      
      CALL ZDISP(0,LOECON,'COECON')
      CALL ZDISP(0,LOFF,'KOFF  ')

      END

************************************************************************
* Generate RHS-vector
*
* Parametric and nonparametric version
*
* This function generates a RHS-vector, based on a user-defined bilinear
* by integration on quadrilaterals. DB points to the vector to 
* calculate.
*
* In contrast to original FEAT routines, this routine only supports
* one matrix-block, thus the NBLOC-parameter of the original XVB0-
* routine is neglected.
*
* In:
*   TRIA   - current triangulation
*   IPARAM - array [1..*] of integer
*            User-defined integer array to pass to the
*            callback-routines
*   DPARAM - array [1..*] of double
*            User-defined double precision array to pass to the
*            callback-routines
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to callback
*            routines. Not used in this routine.
*   NBLOC  - Number of matrix blocks
*   DB     - array [1..NEQ] of double
*            Right hand side vector where to add the RHS-data to
*   NEQ    - integer
*            Number of entries in RHS-vector
*   KVERT  - array [1..NVE,1..NEL] of integer
*            Array with vertex numbers on each element
*   KMID   - array [1..NVE,1..NEL] of integer
*            Array with on of edges each element
*   DCORVG - array [1..2,1..NVT] of double
*            Coordinates of the vertices
*   ELE    - Element to use for calculation
*   BNONPR - Whether ELE is a nonparametric element
*   COEFF  - DOUBLE PRECISION FUNCTION (X,Y,IA,IBLOC,BFIRST,
*                                       TRIA,IEL,IPARAM,DPARAM,
*                                       IGEOM,DGEOM)
*            Function to determine the value of the RHS at a given
*            point.
*            In: 
*              X,Y    - The point where to evaluate the RHS
*              IA     - descriptor of the current summand in the linear
*                       form
*              IBLOC  - Number of current matrix block; =-1 for dummy-
*                       calls
*              BFIRST - =true on first fall, =false on later calls.
*                        allowes to precalculate information
*              IEL    - Number of current element where the cubature
*                       takes place; =0 for dummy-calls
*              TRIA   - Structure of the current triangulation
*              IPARAM,
*              DPARAM - the user defined arrays from above
*              IGEOM,
*              DGEOM  - the user defined arrays from above
*   KOFF    - array [1..NBLOC] of integer
*             Offset between 1 and the starting index of the RHS-vector
*             of block IBLOC in DB-array. Use this if DB=DWORK(1) to
*             describe the real position of the RHS-arrays in DWORK.
*   KBN     - array [1..NBLOC]
*             For every matrix block: Number of summands in the linear
*             form
*   KB      - array [1..KBN,1..NBLOC] of integer
*             Descriptors of the summands in the linear form
*   BCON    - array [1..NBLOC] of boolean
*             For each matrix block:
*             =true, if the coefficients in the integral are constant
*   COECON  - array [1..NNDER,1..NBLOC] of double
*             Auxiliary array for each matrix block, used to save
*             constant coefficients
*   ICUB    - Cubature formula to use for integration on quadrilaterals
*
* Return:
*   DB      - RHS-vector, filled with data
*
* In case of an error, IER will be set <> 0.
************************************************************************

      SUBROUTINE VB0X(TRIA,IPARAM,DPARAM,IGEOM,DGEOM,
     *                NBLOC,DB,KOFF,KVERT,KMID,DCORVG,
     *                ELE,BNONPR,COEFF,
     *                KBN,KB,BCON,COECON,ICUB)
      
      IMPLICIT NONE
 
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
 
      INCLUDE 'ccub.inc'
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'stria.inc'
      
C parameters
      
      INTEGER IPARAM(*),TRIA(SZTRIA),IGEOM(*)
      DOUBLE PRECISION DPARAM(*),DGEOM(*)
      INTEGER NEQ,ICUB,NBLOC
      DOUBLE PRECISION DB(*),DCORVG(2,*),COECON(NNDER,NBLOC)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      INTEGER KBN(NBLOC),KB(NNAB,NBLOC),KOFF(NBLOC)
      LOGICAL BCON(NBLOC),BNONPR
      
      DOUBLE PRECISION COEFF
      EXTERNAL ELE,COEFF
      
C local variables      
      
      DOUBLE PRECISION DCOORD(2,NNVE)
      DOUBLE PRECISION DJF(2,2),XI1,XI2,XX,YY,AUX,OM
      
      INTEGER IBN,IVE,IELTYP,IBLOC,IDER,IB,JP
      LOGICAL BFIRST
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL,JDOFE,IGLOB
      
      INTEGER NDFL
      EXTERNAL NDFL
      
C Preparation - evaluation of parameters

      IER=0

C Which derivatives of basis functions are needed?
C Set the entries in BDER in the /CUB/ common block properly

      DO IDER=1,NNDER
        BDER(IDER)=.FALSE.
      END DO
      
      DO IBLOC=1,NBLOC
        DO IBN=1,KBN(IBLOC)
          IB=KB(IBN,IBLOC)
          IF (IB.LE.0.OR.IB.GT.NNDER) THEN
            WRITE (CPARAM,'(I15)') IBLOC
            CALL WERR(-117,'VB0X  ')
            RETURN
          END IF
          BDER(IB)=.TRUE.
        END DO
      END DO

C Dummy call of ELE to obtain number of element in IELTYP

      IELTYP=-1
      
      CALL ELE(0D0,0D0,IELTYP)
      IF (IER.NE.0) RETURN
  
C Determine degrees of freedom on the elements
      
      IDFL=NDFL(IELTYP)
      IF (IER.LT.0) RETURN
      
C Initialize cubature formula
      
      CALL CB2Q(ICUB)
      IF (IER.NE.0) RETURN
      
      BFIRST=.TRUE.
      
C Dummy call of COEFF for nonlinear problems
C COEFF must set BDER(IDER)=.TRUE. if derivative IDER is needed

      AUX = COEFF(0D0,0D0,-1,0,BFIRST,0,TRIA,IPARAM,DPARAM,IGEOM,DGEOM)

C If we have constant coefficients, calculate them in advance an save
C them in the auxiliary array. Set BCON0=true, if we have only constant
C coefficients

      DO IBLOC=1,NBLOC
        IF (BCON(IBLOC)) THEN
          DO IBN=1,KBN(IBLOC)
            IB=KB(IBN,IBLOC)
            COECON(IB,IBLOC)=COEFF(0D0,0D0,IB,IBLOC,BFIRST,
     *                             IEL,TRIA,IPARAM,DPARAM,IGEOM,DGEOM)
          END DO
        END IF
      END DO

************************************************************************
C *** Calculation of the linear form
************************************************************************

C     Prepare parametric element for cubature. Save the current cubature
C     formula into the COMMON block variable ICUBP - either for now or
C     for later use. Perform a dummy call to the element in the 
C     conforming case.
C     This is done for saving arithmetic operations in later calls.
C
C     In the nonconforming case the element has to be initialised
C     separately on every element, so there's no use in initialising it
C     globally.

      ICUBP=ICUB
      
      IF (.NOT.BNONPR) THEN
        CALL ELE(0D0,0D0,-2)
      END IF

C Loop over all elements, calculate on each element the information and
C sum them up:

      DO IEL=1,TRIA(ONEL)
      
C       Determine global degrees of freedom. The results are saved in
C       KDFG/KDFL in the common blocks.
      
        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) RETURN

C       Save the corners of the current element in the /CUB/ common block
C       for cubature, as well as in DCOORD for the transformation

        DO IVE=1,TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
          DCOORD(1,IVE) = DCORVG(1,JP)
          DCOORD(2,IVE) = DCORVG(2,JP)
        END DO
        
C       Prepare nonparametric elements for cubature on current element.
C       The COMMON block variable ICUBP was initialised earlier...
C
C       Because the later loop will change ICUBP as it is used for  
C       looping through the cubature points on the element, so we have 
C       to reset it here.

        ICUBP = ICUB
        
        IF (BNONPR) THEN
          CALL ELE(0D0,0D0,-2)
        END IF

C       Calculate auxiliary Jacobian factors

        CALL QINIJF (DCOORD,DJF)

C       Loop over all cubature points

        DO ICUBP=1,NCUBP
        
C         Get the cubature point on the reference element
        
          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)

C         Calculate the transformation of that cubature point to the
C         real element, the Jacobian determinant DETJ and the Jacobian
C         matrix DJAC of the transformation:

          CALL QTRAF (DCOORD,DJF,DJAC,DETJ,XI1,XI2,XX,YY)

C         Calculate the weighting factor for the current cubature 
C         point with the help of the Jacobian determinant

          OM = DOMEGA(ICUBP)*DETJ
          
C         Evaluate the basis functions in the cubature point

          IF(BNONPR) THEN
            CALL ELE(XX,YY,-3)
          ELSE
            CALL ELE(XI1,XI2,-3)
          ENDIF
          IF (IER.LT.0) RETURN
          
C         Summing up over all multiindices
 
          BFIRST=.TRUE.
          
          DO IBLOC=1,NBLOC
 
            DO IBN=1,KBN(IBLOC)
            
              IB=KB(IBN,IBLOC)
              
              IF (.NOT.BCON(IBLOC)) THEN
                AUX=COEFF(XX,YY,IB,IBLOC,BFIRST,TRIA,
     *                    IPARAM,DPARAM,IGEOM,DGEOM)*OM
              ELSE
                AUX=COECON(IB,IBLOC)*OM
              ENDIF
        
              DO JDOFE=1,IDFL
                IGLOB=KDFG(JDOFE)+KOFF(IBLOC)
                DB(IGLOB)=DB(IGLOB)+DBAS(KDFL(JDOFE),IB)*AUX
              END DO
              
            END DO
   
            BFIRST=.FALSE.
            
          END DO
          
        END DO
        
      END DO

      END

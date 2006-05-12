**********************************************************************
* This file contains general evaluation routines for FEM-functions
* as well as routines for producing cutlines through the mesh.
**********************************************************************

**********************************************************************
* Scalar evaluation of FEM function on quadrilaterals
*
* Evaluates a scalar FEM-function in a given point. The return value
* is the value in this point as well as the X- and Y-derivative.
*
* In:
*  NEQ    - length of FEM solution vector 
*  DSOL   - array [1..NEQ] of double
*           FEM solution vector describing the continuous FEM
*           function. The function is assumed to be scalar!
*  ELE    - Element that is used for the discretisation
*  BNONPR - TRUE if ELEM is a nonparametric element,
*           FALSE otherwise
*  DXPOS,
*  DYPOS  - the point where to evaluate
*  IELEM  - The element containing the point (DX,DY), or:
*           -1 = determine the element containing (DX,DY) automatically
*            0 = determine the element containing (DX,DY) automatically
*                and return the element number in IELEM
*  TRIA  - array [1..SZTRIA] of integer
*          Triangulation structure
*  DCORVG,
*  KVERT
*  KMID  - Usual geometry information; must correspond to TRIA!
* 
*
* Out:
*  DVAL   - Value of the function in (DX,DY)
*  DGRX   - X-derivative in this point
*  DGRY   - Y-derivative in this point
*
* if IELEM=0 on call:
*  IELEM  - element containing (DX,DY)
*           
* If the point is not on the element (or if no element was found),
* IER is set to -1 and the calculation is aborted.
**********************************************************************

      SUBROUTINE SCEVLQ (NEQ,DSOL,ELE,BNONPR,DXPOS,DYPOS,IELEM,
     *                   TRIA,DCORVG,KVERT,KMID,
     *                   DVAL, DGRX, DGRY)
     
      IMPLICIT NONE
      
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'stria.inc'
      
C parameters
      
      INTEGER NEQ,IELEM,KVERT(NNVE,*),KMID(NNVE,*),TRIA(SZTRIA)
      DOUBLE PRECISION DSOL(NEQ),DXPOS,DYPOS,DVAL, DGRX, DGRY
      DOUBLE PRECISION DCORVG(2,*),DJF(2,2)
      LOGICAL BNONPR
      EXTERNAL ELE
      
C local variables

      INTEGER IEL1,I,IVE,JP,IELTYP,NEL,NVE
      DOUBLE PRECISION X,Y,DCOORD(2,4)
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL
      
      LOGICAL ISINEL
      INTEGER NDFL
      EXTERNAL ISINEL,NDFL
      
C     Fetch some general information to local variables

      NEL = TRIA(ONEL)
      NVE = TRIA(ONVE)
      
C Do we have an element given?

      IF (IELEM.GT.0) THEN
        IEL1=IELEM
      ELSE
      
C No, we have to search for it. Basic linear search method.

        DO IEL1=1,NEL
          IF (ISINEL (DXPOS,DYPOS,IEL1,DCORVG,KVERT)) THEN
            IF (IELEM.EQ.0) IELEM = IEL1
            GOTO 10
          END IF
        END DO
10      CONTINUE        

      END IF
      
C Quick save-test to be sure that element contains our point
      
      IF ((IEL1.GT.NEL).OR.
     *    (.NOT.ISINEL (DXPOS,DYPOS,IEL1,DCORVG,KVERT))) THEN
        IF (MT.GE.3) THEN
          WRITE (*,'(A,2D24.12,A)') 
     *     'WARNING: Found no element containing node (',DXPOS,DYPOS,')'
        END IF
        IER = -1
        GOTO 99999
      END IF      
      
C On parametric elements we have to transform the point back to the
C reference element

      IF (.NOT.BNONPR) THEN
        
C Build the structure defining the element:

        DO IVE = 1,NVE
          JP=KVERT(IVE,IEL1)
          DCOORD (1,IVE) = DCORVG(1,JP)
          DCOORD (2,IVE) = DCORVG(2,JP)
        END DO
        
C Back-transformation

        CALL QBTRAF (DCOORD, X, Y, DXPOS, DYPOS)
        
      ELSE
      
C Otherwise we can evaluate directly

        X = DXPOS
        Y = DYPOS
        
      END IF

C Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)

C how many basis functions=DOF's do we have on each element?

      IDFL=NDFL(IELTYP)

C Calculate the local and global DOF's on our current element.
C We later have to loop about them...

      CALL NDFGLX(TRIA,IEL1,1,IELTYP,KVERT,KMID, KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999

C     Calculate auxiliary Jacobian factors of the transformation

      CALL QINIJF (DCOORD,DJF)

C     Calculate Jacobian matrix of transformation and its
C     determinant. This is necessary for the element routine ELE
C     to properly calculate derivatives. The result is directly
C     written into the element COMMON block variables DJAC and DETJ.

      CALL QTRDET (DCOORD,DJF,DJAC,DETJ,X,Y)

C From the element we want to have weights for function values
C as well as derivatives:

      DO  I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      
C Initialise the COMMON block of the element to inform the it
C about our current quadrilateral:

      DO IVE = 1, NVE
        JP=KVERT(IVE,IEL1)
        KVE(IVE)=JP
        DX(IVE)=DCORVG(1,JP)
        DY(IVE)=DCORVG(2,JP)
      END DO
      
C Call the element routine to calculate the weights of the basis 
C functions in the current point:

      CALL ELE (X,Y,0)
      
C Is there at all something to do, or is the solution vector empty?!?
      
      IF (NEQ.EQ.0) GOTO 99999
      
C Loop about the local DOF's and sum up the values there together with
C the just calculated weights to obtain the function value in (X,Y):

      DVAL = 0D0
      DGRX = 0D0
      DGRY = 0D0

      DO I=1,IDFL
        JP=KDFG(I)
        DVAL=DVAL+DSOL(JP)*DBAS(KDFL(I),1)
        DGRX=DGRX+DSOL(JP)*DBAS(KDFL(I),2)
        DGRY=DGRY+DSOL(JP)*DBAS(KDFL(I),3)
      END DO
      
99999 END

**********************************************************************
* Scalar evaluation of FEM function on quadrilaterals
* with reconstructed gradients
*
* Evaluates a scalar FEM-function in a given point. The return value
* is the value in this point as well as the X- and Y-derivative.
* In contrast to SCEVLQ, this routine accepts two additional
* arrays DGRADX, DGRADY containing for each node the X- and Y-
* gradient of the function. These arrays can e.g. be calculated with
* reconstructed gradient routines. The routine then returns the
* gradient of the function based on these arrays rather than on the
* original gradient of the FEM-function.
*
* In:
*  NEQ    - length of FEM solution vector 
*  DSOL   - array [1..NEQ] of double
*           FEM solution vector describing the continuous FEM
*           function. The function is assumed to be scalar!
*  DGRADX - array [1..NEQ] of double
*  DGRADX - array [1..NEQ] of double
*           FEM solution vector describing the gradient of the
*           continuous FEM function DSOL.
*  ELE    - Element that is used for the discretisation
*  BNONPR - TRUE if ELEM is a nonparametric element,
*           FALSE otherwise
*  DXPOS,
*  DYPOS  - the point where to evaluate
*  IELEM  - The element containing the point (DX,DY), or:
*           -1 = determine the element containing (DX,DY) automatically
*            0 = determine the element containing (DX,DY) automatically
*                and return the element number in IELEM
*  TRIA  - array [1..SZTRIA] of integer
*          Triangulation structure
*  DCORVG,
*  KVERT
*  KMID  - Usual geometry information; must correspond to TRIA!
* 
*
* Out:
*  DVAL   - Value of the function in (DX,DY)
*  DGRX   - X-derivative in this point
*  DGRY   - Y-derivative in this point
*
* if IELEM=0 on call:
*  IELEM  - element containing (DX,DY)
*           
* If the point is not on the element (or if no element was found),
* IER is set to -1 and the calculation is aborted.
**********************************************************************

      SUBROUTINE SCEVQR (NEQ,DSOL,DGRADX,DGRADY,ELE,BNONPR,
     *                   DXPOS,DYPOS,IELEM,
     *                   TRIA,DCORVG,KVERT,KMID,
     *                   DVAL, DGRX, DGRY)
     
      IMPLICIT NONE
      
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'stria.inc'
      
C parameters
      
      INTEGER NEQ,IELEM,KVERT(NNVE,*),KMID(NNVE,*),TRIA(SZTRIA)
      DOUBLE PRECISION DSOL(NEQ),DGRADX(NEQ), DGRADY(NEQ)
      DOUBLE PRECISION DXPOS,DYPOS,DVAL, DGRX, DGRY
      DOUBLE PRECISION DCORVG(2,*)
      LOGICAL BNONPR
      EXTERNAL ELE
      
C local variables

      INTEGER IEL1,I,IVE,JP,IELTYP,NEL,NVE
      DOUBLE PRECISION X,Y,DCOORD(2,4)
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL
      
      LOGICAL ISINEL
      INTEGER NDFL
      EXTERNAL ISINEL,NDFL
      
C     Fetch some general information to local variables

      NEL = TRIA(ONEL)
      NVE = TRIA(ONVE)
      
C Do we have an element given?

      IF (IELEM.GT.0) THEN
        IEL1=IELEM
      ELSE
      
C No, we have to search for it. Basic linear search method.

        DO IEL1=1,NEL
          IF (ISINEL (DXPOS,DYPOS,IEL1,DCORVG,KVERT)) THEN
            IF (IELEM.EQ.0) IELEM = IEL1
            GOTO 10
          END IF
        END DO
10      CONTINUE        

      END IF
      
C Quick save-test to be sure that element contains our point
      
      IF ((IEL1.GT.NEL).OR.
     *    (.NOT.ISINEL (DXPOS,DYPOS,IEL1,DCORVG,KVERT))) THEN
        IF (MT.GE.3) THEN
          WRITE (*,'(A,2D24.12,A)') 
     *     'WARNING: Found no element containing node (',DXPOS,DYPOS,')'
        END IF
        IER = -1
        GOTO 99999
      END IF      
      
C On parametric elements we have to transform the point back to the
C reference element

      IF (.NOT.BNONPR) THEN
        
C Build the structure defining the element:

        DO IVE = 1,NVE
          JP=KVERT(IVE,IEL1)
          DCOORD (1,IVE) = DCORVG(1,JP)
          DCOORD (2,IVE) = DCORVG(2,JP)
        END DO
        
C Back-transformation

        CALL QBTRAF (DCOORD, X, Y, DXPOS, DYPOS)
        
      ELSE
      
C Otherwise we can evaluate directly

        X = DXPOS
        Y = DYPOS
        
      END IF

C Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)

C how many basis functions=DOF's do we have on each element?

      IDFL=NDFL(IELTYP)

C Calculate the local and global DOF's on our current element.
C We later have to loop about them...

      CALL NDFGLX(TRIA,IEL1,1,IELTYP,KVERT,KMID, KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999

C From the element we want to have weights for function values.
C The values of the derivatice we don't need, as we can calculate
C the derivatives from DGRADX/DGRADY.

      DO  I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      
      BDER(1)=.TRUE.
      
C ---
C In standard evaluation of elements at this point we have to calculate
C the Jacobian and its determinant. But this is only necessary for the
C evaluation of derivatives. So as the derivatives are calculated with
C reconstructed gradients only by evaluation of functionals, we can 
C skip the calculation of the Jacobian with QINIJF / QTRDET here.
C ---      
      
C Initialise the COMMON block of the element to inform the it
C about our current quadrilateral:

      DO IVE = 1, NVE
        JP=KVERT(IVE,IEL1)
        KVE(IVE)=JP
        DX(IVE)=DCORVG(1,JP)
        DY(IVE)=DCORVG(2,JP)
      END DO
      
C Call the element routine to calculate the weights of the basis 
C functions in the current point:

      CALL ELE (X,Y,0)
      
C Is there at all something to do, or is the solution vector empty?!?
      
      IF (NEQ.EQ.0) GOTO 99999
      
C Loop about the local DOF's and sum up the values there together with
C the just calculated weights to obtain the function value in (X,Y):

      DVAL = 0D0
      DGRX = 0D0
      DGRY = 0D0

      DO I=1,IDFL
        JP=KDFG(I)
        DVAL=DVAL+DSOL(JP)*DBAS(KDFL(I),1)
        DGRX=DGRX+DGRADX(JP)*DBAS(KDFL(I),1)
        DGRY=DGRY+DGRADY(JP)*DBAS(KDFL(I),1)
      END DO
      
99999 END

**********************************************************************
* Multiple scalar evaluations of FEM function on quadrilaterals
*
* Evaluates a set of scalar FEM-functions in a given point. This helps
* to save computational time if more than one function is to evaluate.
* The FEM-functions are given as a set of handles into the DWORK-
* array. The caller can decide whether to calculate function value or
* gradient of every function. The subroutine returns the value
* of the appropriate function and/or derivatives in three arrays
* with each element of the array corresponding to a function.
*
* This routine can also be used to calculate the gradient of a
* function using recovered gradients. In this case the two arrays
* containing X- and Y-derivative of a function have to be given
* as a standard function.
*
* In:
*  NSOL   - Number of solutions
*  NEQ    - Length of each FEM solution vector.
*           It's assumed that all solution vectors correspond to the
*           same element type and have the same length.
*  LSOL   - array [1..NSOL] of integer
*           Array of handles to all the NSOL solution vectors
*  FCLC   - array [1..NSOL] of integer
*           For every solution vector: Bitfield. FCLC(I) determines
*           what to calculate for that specific solution LSOL(I):
*             Bit 0: Calculate the function value into DVAL(I)
*             Bit 1: Calculate the X-derivative into DGRADY(I)
*             Bit 2: Calculate the Y-derivative into DGRADY(I)
*  ELE    - Element that is used for the discretisation
*  BNONPR - TRUE if ELEM is a nonparametric element,
*           FALSE otherwise
*  DXPOS,
*  DYPOS  - the point where to evaluate
*  IELEM  - The element containing the point (DX,DY), or:
*           -1 = determine the element containing (DX,DY) automatically
*            0 = determine the element containing (DX,DY) automatically
*                and return the element number in IELEM
*  TRIA  - array [1..SZTRIA] of integer
*          Triangulation structure
*  DCORVG,
*  KVERT
*  KMID  - Usual geometry information; must correspond to TRIA!
*
* Out:
*  DVAL   - array [1..NSOL] of double
*           DVAL(I) contains the value of each function LSOL(I)
*           in the point (DX,DY)
*  DGRX   - array [1..NSOL] of double
*           X-derivative in this point for every function, if to be 
*           calculated
*  DGRY   - array [1..NSOL] of double
*           Y-derivative in this point for every function, if to be 
*           calculated
*
* if IELEM=0 on call:
*  IELEM  - element containing (DX,DY)
*           
* If the point is not on the element (or if no element was found),
* IER is set to -1 and the calculation is aborted.
**********************************************************************

      SUBROUTINE MSCVQR (NSOL,NEQ,LSOL,FCLC,ELE,BNONPR,
     *                   DXPOS,DYPOS,IELEM,
     *                   TRIA,DCORVG,KVERT,KMID,
     *                   DVAL, DGRX, DGRY)
     
      IMPLICIT NONE
      
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'stria.inc'
      
C parameters
      
      INTEGER NSOL,NEQ,IELEM,KVERT(NNVE,*),KMID(NNVE,*),TRIA(SZTRIA)
      INTEGER LSOL(NEQ),FCLC(NSOL)
      DOUBLE PRECISION DXPOS,DYPOS,DVAL(NSOL),DGRX(NSOL),DGRY(NSOL)
      DOUBLE PRECISION DCORVG(2,*)
      LOGICAL BNONPR
      EXTERNAL ELE
      
C local variables

      INTEGER IEL1,I,J,IVE,JP,IELTYP,NEL,NVE
      DOUBLE PRECISION X,Y,DCOORD(2,4),BAS,BASX,BASY,SOL,DJF(2,2)
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL
      
      LOGICAL ISINEL
      INTEGER NDFL
      EXTERNAL ISINEL,NDFL
      
C     Fetch some general information to local variables

      NEL = TRIA(ONEL)
      NVE = TRIA(ONVE)
      
C Do we have an element given?

      IF (IELEM.GT.0) THEN
        IEL1=IELEM
      ELSE
      
C No, we have to search for it. Basic linear search method.

        DO IEL1=1,NEL
          IF (ISINEL (DXPOS,DYPOS,IEL1,DCORVG,KVERT)) THEN
            IF (IELEM.EQ.0) IELEM = IEL1
            GOTO 10
          END IF
        END DO
10      CONTINUE        

      END IF
      
C Quick save-test to be sure that element contains our point
      
      IF ((IEL1.GT.NEL).OR.
     *    (.NOT.ISINEL (DXPOS,DYPOS,IEL1,DCORVG,KVERT))) THEN
        IF (MT.GE.3) THEN
          WRITE (*,'(A,2D24.12,A)') 
     *     'WARNING: Found no element containing node (',DXPOS,DYPOS,')'
        END IF
        IER = -1
        GOTO 99999
      END IF      
      
C On parametric elements we have to transform the point back to the
C reference element

      IF (.NOT.BNONPR) THEN
        
C Build the structure defining the element:

        DO IVE = 1,NVE
          JP=KVERT(IVE,IEL1)
          DCOORD (1,IVE) = DCORVG(1,JP)
          DCOORD (2,IVE) = DCORVG(2,JP)
        END DO
        
C Back-transformation

        CALL QBTRAF (DCOORD, X, Y, DXPOS, DYPOS)
        
      ELSE
      
C Otherwise we can evaluate directly

        X = DXPOS
        Y = DYPOS
        
      END IF

C Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)

C how many basis functions=DOF's do we have on each element?

      IDFL=NDFL(IELTYP)

C Calculate the local and global DOF's on our current element.
C We later have to loop about them...

      CALL NDFGLX(TRIA,IEL1,1,IELTYP,KVERT,KMID, KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999

C From the element we want to have weights for function values.
C The values of the derivatice we don't need, as we can calculate
C the derivatives from DGRADX/DGRADY.

      DO I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      
      DO I=1,NSOL
        IF (IAND(FCLC(I),1).NE.0) BDER(1)=.TRUE.
        IF (IAND(FCLC(I),2).NE.0) BDER(2)=.TRUE.
        IF (IAND(FCLC(I),4).NE.0) BDER(3)=.TRUE.
      END DO
      
C     If derivatives are to be calculated, we need the Jacibian
C     matrix and determinant of the mapping to/from the reference
C     element. If derivatives are not used, we can skip this,
C     as elements don't need information about the Jacobian
C     only for calculating the values of the function.

      IF (BDER(2).OR.BDER(3)) THEN
      
C       Calculate auxiliary Jacobian factors of the transformation

        CALL QINIJF (DCOORD,DJF)

C       Calculate Jacobian matrix of transformation and its
C       determinant. This is necessary for the element routine ELE
C       to properly calculate derivatives. The result is directly
C       written into the element COMMON block variables DJAC and DETJ.

        CALL QTRDET (DCOORD,DJF,DJAC,DETJ,X,Y)
        
      END IF

C Initialise the COMMON block of the element to inform the it
C about our current quadrilateral:

      DO IVE = 1, NVE
        JP=KVERT(IVE,IEL1)
        KVE(IVE)=JP
        DX(IVE)=DCORVG(1,JP)
        DY(IVE)=DCORVG(2,JP)
      END DO
      
C Call the element routine to calculate the weights of the basis 
C functions in the current point:

      CALL ELE (X,Y,0)
      
C Is there at all something to do, or is the solution vector empty?!?
      
      IF (NEQ.EQ.0) GOTO 99999
      
C Clear the return vectors

      DO I=1,NSOL
        DVAL(I) = 0D0
        DGRX(I) = 0D0
        DGRY(I) = 0D0
      END DO

C Loop about the local DOF's and sum up the values there together with
C the just calculated weights to obtain the function value in (X,Y):

      DO I=1,IDFL
        JP=KDFG(I)
        BAS = DBAS(KDFL(I),1)
        BASX = DBAS(KDFL(I),2)
        BASY = DBAS(KDFL(I),3)

C       Loop over all functions, build together function value
C       and derivative

        DO J=1,NSOL
          SOL = DWORK(L(LSOL(J))+JP-1)
        
          IF (IAND(FCLC(J),1).NE.0) THEN
            DVAL(J) = DVAL(J) + SOL*BAS
          END IF
          
          IF (IAND(FCLC(J),2).NE.0) THEN
            DGRX(J) = DGRX(J) + SOL*BASX
          END IF
          
          IF (IAND(FCLC(J),4).NE.0) THEN 
            DGRY(J) = DGRY(J) + SOL*BASY
          END IF
        END DO

      END DO
      
99999 END

**********************************************************************
* Scalar evaluation of FEM function on triangles
*
* Evaluates a scalar FEM-function in a given point. The return value
* is the value in this point as well as the X- and Y-derivative.
*
* In:
*  NEQ    - length of FEM solution vector 
*  DSOL   - array [1..NEQ] of double
*           FEM solution vector describing the continuous FEM
*           function. The function is assumed to be scalar!
*  ELE    - Element that is used for the discretisation
*  BNONPR - TRUE if ELEM is a nonparametric element,
*           FALSE otherwise
*  DXPOS,
*  DYPOS  - the point where to evaluate
*  IELEM  - The element containing the point (DX,DY), or:
*           -1 = determine the element containing (DX,DY) automatically
*            0 = determine the element containing (DX,DY) automatically
*                and return the element number in IELEM
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure
*  DCORVG,
*  KVERT
*  KMID   - Usual geometry information; must correspond to TRIA!
* 
*
* Out:
*  DVAL   - Value of the function in (DX,DY)
*  DGRX   - X-derivative in this point
*  DGRY   - Y-derivative in this point
*
* if IELEM=0 on call:
*  IELEM  - element containing (DX,DY)
*           
* If the point is not on the element (or if no element was found),
* IER is set to -1 and the calculation is aborted.
**********************************************************************

      SUBROUTINE SCEVLT (NEQ,DSOL,ELE,BNONPR,DXPOS,DYPOS,IELEM,
     *                   TRIA,DCORVG,KVERT,KMID,
     *                   DVAL, DGRX, DGRY)
     
      IMPLICIT NONE
      
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'stria.inc'
      
C parameters
      
      INTEGER NEQ,IELEM,KVERT(NNVE,*),KMID(NNVE,*),TRIA(SZTRIA)
      DOUBLE PRECISION DSOL(NEQ),DXPOS,DYPOS,DVAL, DGRX, DGRY
      DOUBLE PRECISION DCORVG(2,*)
      LOGICAL BNONPR
      EXTERNAL ELE
      
C local variables

      INTEGER IEL1,I,J,IVE,JP,IELTYP, IV, NEL,NVE
      DOUBLE PRECISION X1,X2,X3
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL
      DOUBLE PRECISION DCOORD (2,3)
      
      INTEGER NDFL
      EXTERNAL NDFL
      
C     Fetch some general information to local variables

      NEL = TRIA(ONEL)
      NVE = TRIA(ONVE)
      
C Do we have an element given?

      IF (IELEM.GT.0) THEN

        IEL1 = IELEM

C Ok, let's check if the point is really inside of this triangle.
C Calculate the barycentric coordinates:

        DO J=1,3
          IV = KVERT(J,IEL1)
          DCOORD(1,J) = DCORVG(1, IV)
          DCOORD(2,J) = DCORVG(2, IV)
        END DO
        
        CALL TBTRAF (DCOORD, X1, X2, X3, DETJ, DXPOS, DYPOS)
	
        IF ((X1.LT.0D0).OR.(X1.GT.1D0).OR.
     *      (X2.LT.0D0).OR.(X2.GT.1D0).OR.
     *      (X3.LT.0D0).OR.(X3.GT.1D0)) THEN

C Oh, not good, the point is not inside here.
C Abort the further calculation. The caller can call this routine again
C to perform a linear search if desired...

          IER = -1
          GOTO 99999
          
	  ENDIF

      ELSE
	      
C No, we have to search for it. Basic linear search method.

        DO IEL1=1, NEL

C Calculate the barycentric coordinates of the point corresponding to
C the current triangle:
	  
          DO J=1,3
            IV = KVERT(J,IEL1)
            DCOORD(1,J) = DCORVG(1, IV)
            DCOORD(2,J) = DCORVG(2, IV)
          END DO
          
          CALL TBTRAF (DCOORD, X1, X2, X3, DETJ, DXPOS, DYPOS)
          
C Then check if the coordinates are all between 0 and 1;
C if yes, we found our element, if no... go on with the search
          
	    IF ((X1.GE.0D0).AND.(X1.LE.1D0).AND.
     *        (X2.GE.0D0).AND.(X2.LE.1D0).AND.
     *        (X3.GE.0D0).AND.(X3.LE.1D0)) THEN
	
C Is the caller interested in the element number?
        
	      IF (IELEM.EQ.0) IELEM=IEL1
	      GOTO 10
    	  
	    ENDIF
        END DO
10    END IF

      IF(IEL1.GT.NEL) THEN
        IF (MT.GE.3) THEN
	    WRITE (MTERM,'(A,2D24.12)') 
     *          'Warning, found no element containing',DXPOS, DYPOS
          IER=-1
          GOTO 99999
        END IF
      END IF
      
C Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
	
C how many basis functions=DOF's do we have on each element?

      IDFL=NDFL(IELTYP)

C Calculate the local and global DOF's on our current element.
C We later have to loop about them...

      CALL NDFGLX(TRIA,IEL1,1,IELTYP,KVERT,KMID, KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999

C !!!!!!!!!!!!!!!!
C Perhaps at this point the calculation of the Jacobian and its determinant
C is necessary - although this is not necessary if only function values
C are to be calculated. As the routines has never been tested for the
C evaluation of derivatives of FE-functions on triangles, this might
C be an open bug.
C !!!!!!!!!!!!!!!!

C From the element we want to have weights for function values
C as well as derivatives:

      DO  I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      
C Initialise the COMMON block of the element to inform the it
C about our current quadrilateral:
	
      DO IVE = 1, NVE
        JP=KVERT(IVE,IEL1)
        KVE(IVE)=JP
        DX(IVE)=DCORVG(1,JP)
        DY(IVE)=DCORVG(2,JP)
      END DO
      
C Call the element routine to calculate the weights of the basis 
C functions in the current point:

      CALL ELE (X1,X2,X3,0)
      
C Is there at all something to do, or is the solution vector empty?!?
      
      IF (NEQ.EQ.0) GOTO 99999
      
C Loop about the local DOF's and sum up the values there together with
C the just calculated weights to obtain the function value in (X,Y):

      DVAL = 0D0
      DGRX = 0D0
      DGRY = 0D0

      DO I=1,IDFL
        JP=KDFG(I)
        DVAL=DVAL+DSOL(JP)*DBAS(KDFL(I),1)
        DGRX=DGRX+DSOL(JP)*DBAS(KDFL(I),2)
        DGRY=DGRY+DSOL(JP)*DBAS(KDFL(I),3)
      END DO
	
99999 END


**********************************************************************
* Calculate values of a scalar FEM function on a cutline
*
* Evaluates a scalar FEM-function along a given cutline.
* The cutline is given by start- and endpoint and is evaluated
* in equidistant intervals. 
*
* In:
*  NEQ    - length of FEM solution vector 
*  DSOL   - array [1..NEQ] of double
*           FEM solution vector describing the continuous FEM
*           function. The function is assumed to be scalar!
*  ELE    - Element that is used for the discretisation
*  BNONPR - TRUE if ELEM is a nonparametric element,
*           FALSE otherwise
*  BTRI   - TRUE if ELEM is a triangular element,
*           FALSE otherwise
*  DX1, 
*  DY1    - start point of the cutline
*  DX2,
*  DY2    - end point of the cutline
*  NNDS   - Number of nodes on the cutline where the function
*           should be evaluated.
*
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure
*  DCORVG,
*  KVERT
*  KMID   - Usual geometry information; must correspond to TRIA!
* 
*
* Out:
*  DVAL   - array [1..NNDS] of double
*           Receives the function values along the cutline.
*  DGRX   - array [1..NNDS] of double
*           Receives the X-derivatives along the cutline.
*  DGRY   - array [1..NNDS] of double
*           Receives the Y-derivatives values along the cutline.
*  DXPOS  - array [1..NNDS] of double
*           Recevies the X-positions of the points where the
*           function is evaluated
*  DYPOS  - array [1..NNDS] of double
*           Recevies the Y-positions of the points where the
*           function is evaluated
*  IELS   - array [1..NNDS] of integer
*           Recevies the numbers of the elements corresponding to
*           the points in DXPOS/DYPOS
*
* If the point is not on the element (or if no element was found),
* IER is set to -1 and the calculation is aborted.
**********************************************************************

      SUBROUTINE CUTLIN (NEQ,DSOL,ELE,BNONPR, BTRI, 
     *                   DX1,DY1, DX2,DY2, NNDS,
     *                   TRIA,DCORVG,KVERT,KMID,
     *                   DVAL, DGRX, DGRY, DXPOS, DYPOS, IELS)

      IMPLICIT NONE
      
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
C parameters
      
      INTEGER NEQ,NNDS,IELS(NNDS),KVERT(NNVE,*),KMID(NNVE,*)
      INTEGER TRIA(SZTRIA)
      DOUBLE PRECISION DCORVG(2,*), DX1,DY1, DX2,DY2
      DOUBLE PRECISION DSOL(NEQ),DXPOS(NNDS),DYPOS(NNDS)
      DOUBLE PRECISION DVAL(NNDS), DGRX(NNDS), DGRY(NNDS)
      LOGICAL BNONPR, BTRI
      EXTERNAL ELE
      
C local variables

      INTEGER I
      
C Cutlines are postprocessing :)
C Therefore we don't care about runtime...

C The cutline is divided into NNDS-1 intervals. The coordinates
C of the interval endpoints we can calculate directly...

      DO I=1,NNDS
        IELS(I) = 0
        DXPOS(I) = ( DX1 * DBLE(NNDS-I) + DX2 * DBLE(I) ) / DBLE(NNDS)
        DYPOS(I) = ( DY1 * DBLE(NNDS-I) + DY2 * DBLE(I) ) / DBLE(NNDS)
      END DO

C Now loop again through these points and evaluate the function there

      DO I=1,NNDS

C First try to evaluate on the same element as before.

        IF (I.GT.1) IELS(I) = IELS(I-1)
        IF (BTRI) THEN
          CALL SCEVLT (NEQ,DSOL,ELE,BNONPR,DXPOS(I),DYPOS(I),IELS(I),
     *               TRIA,DCORVG,KVERT,KMID,
     *               DVAL(I), DGRX(I), DGRY(I))
        ELSE
          CALL SCEVLQ (NEQ,DSOL,ELE,BNONPR,DXPOS(I),DYPOS(I),IELS(I),
     *               TRIA,DCORVG,KVERT,KMID,
     *               DVAL(I), DGRX(I), DGRY(I))
        END IF
     
C If this doesn't work, search for another element that contains
C the current point.

        IF (IER.LT.0) THEN
          IER = 0
          IELS(I) = 0
          IF (BTRI) THEN
            CALL SCEVLT (NEQ,DSOL,ELE,BNONPR,DXPOS(I),DYPOS(I),IELS(I),
     *               TRIA,DCORVG,KVERT,KMID,
     *               DVAL(I), DGRX(I), DGRY(I))
          ELSE
            CALL SCEVLQ (NEQ,DSOL,ELE,BNONPR,DXPOS(I),DYPOS(I),IELS(I),
     *               TRIA,DCORVG,KVERT,KMID,
     *               DVAL(I), DGRX(I), DGRY(I))
          END IF
        END IF
        
C In case of an error we take the last calculated value.
C Errors e.g. occur if there are two boundary components, so the
C points in the inner boundary component don't have a corresponding
C element in the grid!

        IF (IER.LT.0) THEN
          IER = 0
          IF (I.GT.0) THEN
            DVAL(I) = DVAL(I-1)
            DGRX(I) = DGRX(I-1)
            DGRY(I) = DGRY(I-1)
          END IF
        END IF
      
      END DO
      
      END
      
***********************************************************************
* Compute a cutline and write it to a file
*
* In:
*  DSOL    - solution vector
*  NEQ     - length of solution vector
*  ELE     - Element function
*  BNONPR  - TRUE if the element is nonparametric
*  BTRI   - TRUE if ELEM is a triangular element,
*           FALSE otherwise
*  X1,
*  Y1      - start of the cutline
*  X2,
*  Y2      - end of the cutline
*  NNDCL   - Number of points in the cutline
*  MUNIT   - File handle where to write the data to.
*            The file must be open already!
*
* This file uses the data in the TRIAx-COMMON blocks for the definition
* of the grid.
***********************************************************************

      SUBROUTINE WCUTLN (DSOL,NEQ,ELE,BNONPR,BTRI,X1,Y1,X2,Y2,NNDCL,
     *                   MUNIT)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'stria.inc'
      
      INTEGER NEQ
      DOUBLE PRECISION X1,Y1,X2,Y2,DSOL(NEQ)

      EXTERNAL ELE

      INTEGER MUNIT,NNDCL,ELE
      LOGICAL BNONPR,BTRI

C Some variables for calculating cutlines in the postprocessing

      INTEGER LVAL, LGRX, LGRY, LXPOS, LYPOS, LIELS, TRIA(SZTRIA)
      
      INTEGER I

C Reserve memory for the cutline data

      IER = 0
      
      CALL ZNEW (NNDCL,1,LVAL,'DVAL  ')
      CALL ZNEW (NNDCL,1,LGRX,'DGRY  ')
      CALL ZNEW (NNDCL,1,LGRY,'DGRY  ')
      CALL ZNEW (NNDCL,1,LXPOS,'DXPOS ')
      CALL ZNEW (NNDCL,1,LYPOS,'DYPOS ')
      CALL ZNEW (NNDCL,3,LIELS,'IELS  ')
      
      IF (IER.NE.0) THEN
        WRITE (MTERM,'(A)') 'WCUTLN: not enough memory for cutline'
        RETURN
      END IF

C Fetch the triangulation structure
      
      CALL C2TRIA(TRIA)

C Calculate the cutline

      CALL CUTLIN (NEQ,DSOL,ELE,BNONPR,BTRI, 
     *               X1,Y1, X2,Y2, NNDCL,
     *               TRIA,DWORK(L(TRIA(OLCORVG))),
     *               KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLMID))),
     *               DWORK(L(LVAL)), DWORK(L(LGRX)), DWORK(L(LGRY)),
     *               DWORK(L(LXPOS)), DWORK(L(LYPOS)), KWORK(L(LIELS)))
      
C Write it to the file in GNUPlot-format
      
      WRITE (MUNIT,'(A)') '# X Y Value X-Der. Y-Der. Element-Nr.'

      DO I=1,NNDCL
        WRITE (MUNIT,'(I10,5E24.14,I10)') I,
     *         DWORK(L(LXPOS)+I-1), DWORK(L(LYPOS)+I-1), 
     *         DWORK(L(LVAL)+I-1), DWORK(L(LGRX)+I-1), 
     *         DWORK(L(LGRY)+I-1),KWORK(L(LIELS)+I-1)
      
      END DO

      WRITE (MUNIT,'(A)') ''
      
C Free the memory used
      
      CALL ZDISP (0,LIELS,'IELS  ')
      CALL ZDISP (0,LYPOS,'DYPOS ')
      CALL ZDISP (0,LXPOS,'DXPOS ')
      CALL ZDISP (0,LGRY,'DGRY  ')
      CALL ZDISP (0,LGRX,'DGRY  ')
      CALL ZDISP (0,LVAL,'DVAL  ')
      
      END


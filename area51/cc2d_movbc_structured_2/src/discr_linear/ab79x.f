************************************************************************
* Generate matrix in matrix structure 7/8, wrapper routine
*
* Extended calling convention; quadrilateral version
*
* This routine calls AB7X or ABN7X to generate a the structure of 
* a matrix.
*
* In:
*   LA     : array [1..NBLOC] of integer
*            = 0: Create a new matrix
*            <>0: Handle to an existing matrix which should be
*                 rebuild or updated
*   LCOL   : Handle to the KCOL array; must be the same
*            for all matrix blocks
*   LLD    : Handle to the KLD array; must be the same
*            for all matrix blocks
*   NA     : Number of entries in each matrix block
*   NEQ    : Dimension/Number of rows in the matrix
*   NBLOC  : Number of independent diagonal blocks in the matrix
*   ICLEAR : =0: Add the matrix to an existing one
*            =1: Clear the matrix before building
*   TRIA   : array [1..SZTRIA] of integer
*            Triangulation structure of the underlying mesh
*   ELE    : Finite Element callback routine.
*            Must be quadrilateral, corresponding to
*            the triangulation in TRIA!
*   BNONPR : =false: parametric element ELE
*            =true : nonparametric element ELE
*   COEFF  : Coefficient function for the coefficients in the 
*            bilinear form.
*            DOUBLE PRECISION FUNCTION COEFF(XX,YY,IA,IB,IBLOC,
*                   BFIRST,TRIA,IPARAM,DPARAM)
*   KAB    : array [2,NNAB,KABN] of integer
*            Descriptors for the bilinear form
*   KABN   : Number of additive terms in the bilinear form
*   ICUB   : Cubature formula identifier for CB2Q
*   ISYMM  : =0: don't respect symmetry, create matrix in structure 7
*            =1: matrix is symmetric; create matrix in structure 8
*   ARR    : string[6]
*            Name of the matrix array LA - for output purposes.
*   IPARAM : array [1..*] of integer
*            User defined parameter block, passed to COEFF
*   DPARAM : array [1..*] of double
*            User defined parameter block, passed to COEFF
* Out:
*   If LA=0: LA=handles to the matrices
*   The entries the matrix-array corresponding to LA are filled
*   with data according to the parameters.
************************************************************************

      SUBROUTINE XAB7X(LA,LCOL,LLD,NA,NEQ,NBLOC,ICLEAR,TRIA,ELE,BNONPR,
     *                 COEFF,BCON,KAB,KABN,ICUB,ISYMM,ARR,IPARAM,DPARAM)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'stria.inc'

C     parameters

      INTEGER LA(*),KAB(2,NNAB,*),KABN(*)
      LOGICAL BCON(*)

      CHARACTER ARR*12
      DIMENSION ARR(*)
      
      INTEGER LCOL, LLD, NA, NEQ, NBLOC, ICLEAR, ICUB, ISYMM
      INTEGER TRIA(SZTRIA)

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      DOUBLE PRECISION COEFF
      EXTERNAL COEFF,ELE
      
      LOGICAL BNONPR

C     local variables
      
      INTEGER IBLOC, ITYPE, ILEN, LOFF, LOECON,IELTYP
      INTEGER KVERT,KMID,KCORVG

      SUB='XABM7'
      IF (ICHECK.GE.997) CALL OTRC('XABM7 ','12/08/94')
      
      IER=0

C     Loop through the NBLOC independent blocks of the global
C     matrix:
C         MATRIX =  [ BLOCK1                       ]
C                   [         BLOCK2               ]
C                   [                 ...          ]
C                   [                      BLOCKn  ]
C     All blocks have the same column/row structure.

      DO IBLOC=1,NBLOC
      
C       Create a new matrix block? Check if LA is given:
      
        IF (LA(IBLOC).EQ.0) THEN

C         Yes! Allocate a new matrix.
        
          CALL ZNEW(NA,1,LA(IBLOC),ARR(IBLOC))
          IF (IER.NE.0) GOTO 99999
          
        ELSE
          
C         A matrix is given; check that it has the right dimension
C         and the right data type!
          
          CALL ZTYPE(LA(IBLOC),ITYPE)
          IF (ITYPE.NE.1) THEN
            WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
            CALL WERR(-114,'XABM7 ')
            GOTO 99999
          ENDIF
          
          CALL ZLEN(LA(IBLOC),ILEN)
          IF (ILEN.LT.NA) THEN
            WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
            CALL WERR(-115,'XABM7 ')
            GOTO 99999
          ENDIF
          
C         Clear the array if ICLEAR=1
          
          IF (ICLEAR.EQ.1) CALL ZCLEAR(LA(IBLOC),ARR(IBLOC))
          
        ENDIF
        
      END DO

C     Create an auxiliary array KOFF(NBLOC), which receives
C     for every matrix the starting address of the matrix
C     block in DWORK. 

      CALL ZNEW(NBLOC,-3,LOFF,'KOFF  ')
      IF (IER.NE.0) GOTO 99999
      
      DO IBLOC=1,NBLOC
        KWORK(L(LOFF)+IBLOC-1)=L(LA(IBLOC))-1
      END DO
      
C     Allocate an auxiliary array COECON.
C     In case of a bilinear form with constant coefficients, the
C     coefficients are determined in advance by ABx7X and stored here.
      
      CALL ZNEW(NBLOC*NNDER**2,1,LOECON,'COECON')
      IF (IER.NE.0) GOTO 99999

C     Call AB7X / ABN7X depending on whether we have a parametric
C     or nonparametric element.

      KVERT  = L(TRIA(OLVERT))
      KMID   = L(TRIA(OLMID))
      KCORVG = L(TRIA(OLCORVG))
      
C     Ask the element about its type; switch to nonparametric if
C     the number is < 0 -- this overrides BNONPR in case the caller
C     is wrong. Works only for those nonparametric elements
C     supporting element numbers < 0.
      
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      
      IF ((.NOT.BNONPR).OR.(IELTYP.LE.0)) THEN
      
        CALL AB7X(DWORK(1),KWORK(L(LCOL)),KWORK(L(LLD)),NA,NEQ,
     *          NBLOC,KWORK(L(LOFF)),KWORK(KVERT),KWORK(KMID),
     *          DWORK(KCORVG),TRIA,ELE,COEFF,
     *          BCON,DWORK(L(LOECON)),KAB,KABN,ICUB,ISYMM,
     *          IPARAM,DPARAM)

      ELSE 

        CALL AB7NX(DWORK(1),KWORK(L(LCOL)),KWORK(L(LLD)),NA,NEQ,
     *          NBLOC,KWORK(L(LOFF)),KWORK(KVERT),KWORK(KMID),
     *          DWORK(KCORVG),TRIA,ELE,COEFF,
     *          BCON,DWORK(L(LOECON)),KAB,KABN,ICUB,ISYMM,
     *          IPARAM,DPARAM)
     
      END IF

      IF (IER.NE.0) GOTO 99999
      
C     Release auxiliary arrays.

      CALL ZDISP(0,LOECON,'COECON')
      CALL ZDISP(0,LOFF,'KOFF  ')

99999 END
      
************************************************************************
* Generate matrix in matrix structure 7/8
*
* Extended calling convention; quadrilateral version
*
* This routine creates the entries of a matrix in matrix structure
* 7 or 8 of a bilinear form a(u,v). ELE is used for both, trial and
* test functions. ELE is assumed to be parametric.
*
* In:
*   DA     : array [1..*] of double
*            Matrix which entries are to be generated.
*            Can be <> 0.0, in this case the corresponding new entries
*            are added to the existing ones.
*            Multiple matrix blocks can be calculated at once (providing
*            the matrices have the same structure). The array KOFF
*            describes the starting addresses of the matrices in DA.
*   KCOL   : array [1..NA] of integer
*            Column structure
*   KLD    : array [1..NEQ+1] of integer
*            Row structure
*   NA     : Number of entries in each matrix block
*   NEQ    : Dimension/Number of rows in the matrix
*   NBLOC  : Number of independent diagonal blocks in the matrix
*   KOFF   : array [1..NBLOC] of integer
*            Start indices of the matrices in DA.
*            Matrix I is expected in DA(KOFF(I)..KOFF(I)+NA)
*   DCORVG,
*   KVERT,
*   KMID   : Geometry information arrays of a quadrilateral mesh;
*            must be the arrays given by TRIA, but passed separately
*            for performance reasons.
*   TRIA   : array [1..SZTRIA] of integer
*            Triangulation structure of the underlying mesh
*   ELE    : Finite Element callback routine.
*            Must be quadrilateral, corresponding to
*            the triangulation in TRIA!
*   COEFF  : Coefficient function for the coefficients in the 
*            bilinear form.
*            DOUBLE PRECISION FUNCTION COEFF(XX,YY,IA,IB,IBLOC,
*                   BFIRST,TRIA,IPARAM,DPARAM)
*   BCON   : array [1..NBLOC] of boolean
*            BCON(I)=true,  if the matrix block I has constant 
*                           coefficients
*            BCON(I)=false, otherwise
*   COECON : array [COECON(NNDER,NNDER,NBLOC] of double
*            Auxiliary array, saves the constant coefficients during
*            computation
*   KAB    : array [2,NNAB,KABN] of integer
*            Descriptors for the bilinear form
*   KABN   : Number of additive terms in the bilinear form
*   ICUB   : Cubature formula identifier for CB2Q
*   ISYMM  : =0: don't respect symmetry, create matrix in structure 7
*            =1: matrix is symmetric; create matrix in structure 8
*            =2: matrix is symmetric and in structure 7 or 8;
*                automatically determine matrix structure
*   IPARAM : array [1..*] of integer
*            User defined parameter block, passed to COEFF
*   DPARAM : array [1..*] of double
*            User defined parameter block, passed to COEFF
* Out:
*   The entries of the matrix are added to DA.
*
* Out (COMMON block):
*   IER    = -116 for wrong value in array KAB
************************************************************************

      SUBROUTINE AB7X(DA,KCOL,KLD,NA,NEQ,NBLOC,KOFF,KVERT,KMID,
     *                DCORVG,TRIA,ELE,COEFF,BCON,COECON,
     *                KAB,KABN,ICUB,ISYMM,IPARAM,DPARAM)

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
      
      INTEGER KCOL(*),KVERT(NNVE,*),KMID(NNVE,*),KOFF(*)
      INTEGER KLD(*),KAB(2,NNAB,*),KABN(*),TRIA(SZTRIA)
      INTEGER KENTRY(NNBAS,NNBAS)
      DOUBLE PRECISION DA(*),DCORVG(2,*),COECON(NNDER,NNDER,*)
      LOGICAL BCON(*)
      INTEGER NEQ,NA,ICUB,ISYMM,NBLOC

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      DOUBLE PRECISION COEFF
      EXTERNAL ELE,COEFF

C     externals

      INTEGER NDFL
      EXTERNAL NDFL

C     local variables

      LOGICAL BSYMM,BFIRST,BCON0
      INTEGER I,J,IBLOC,I1,IELTYP,IA,IB,ILD,ICOL,JCOL
      INTEGER IDOFE,JDOFE,JDOFE1,JCOL0,IDFG
      INTEGER IVE,JP,IALBET,JCOLB,IROW,J1,ICOL1,ICOLB,ICOL1B
      DOUBLE PRECISION AUX,DJ1,DJ2,DJ3,DJ4,XI1,XI2,XX,YY,OM,DB
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL

      SUB='ABM7'
      IF (ICHECK.GE.997) CALL OTRC('ABM7  ','12/08/94')

C     Preparation - evaluation of parameters

      IER=0
      
C     Create symmetric matrix in matrix structure 8?      
      
      BSYMM=ISYMM.GE.1
      
C     Which derivatives of basis functions are needed?
C     Check the descriptors of the bilinear form and set BDER
C     according to these.

      DO I = 1,NNDER
        BDER(I)=.FALSE.
      END DO
      
      DO IBLOC=1,NBLOC
        DO I=1,KABN(IBLOC)
          DO J=1,2
          
C           The desriptor KAB gives directly the derivative
C           which is to be computed!

            I1=KAB(J,I,IBLOC)
            
            IF (I1.LE.0.OR.I1.GT.NNDER) THEN
            
C             Oops, invalid descriptor
            
              WRITE (CPARAM,'(I15)') IBLOC
              CALL WERR(-116,'AB07  ')
              GOTO 99999
            ENDIF
            
            BDER(I1)=.TRUE.
          END DO
        END DO
      END DO
      
C     Dummy call of ELE to get the element identifier

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      
C     Get the number of local DOF's:
      
      IDFL=NDFL(IELTYP)
      
C     Initialise the cubature formula
      
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999

C     Dummy call of COEFF for nonlinear problems.
C     COEFF must set BDER(IDER)=.TRUE. in the /ELEM/ COMMON block
C     derivative IDER is needed from the element.
C     This allows to tell the element to calculate additional
C     derivarives to the ones defined by KAB.

      BFIRST=.TRUE.
      AUX=COEFF(0D0,0D0,-1,-1,0,BFIRST,TRIA,IPARAM,DPARAM)
      
C     Check all matrix blocks if their bilinear form has constant
C     coefficients. In that case, get the coefficient and save it
C     in COECON.
C     BCON0 is set to TRUE if all blocks have constant coefficients.      
      
      BCON0=.TRUE.
      DO IBLOC=1,NBLOC
        IF (BCON(IBLOC)) THEN
          DO I=1,KABN(IBLOC)
            IA=KAB(1,I,IBLOC)
            IB=KAB(2,I,IBLOC)
            COECON(IA,IB,IBLOC)=COEFF(0D0,0D0,IA,IB,IBLOC,BFIRST,TRIA,
     *                                IPARAM,DPARAM)    
          END DO
        ELSE
          BCON0=.FALSE.
        ENDIF
      END DO
     
C     Now let's calculate the matrix - storage technique 7 or 8
C
C     Dummy call - ELE may save arithmetic operations.
C     ELE can precalculate the values in the cubature points if
C     it supports that...

      CALL ELE(0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999

C     For the loop through the DOF's later, set JDOFE1=1.
C     For symmetric matrices (structure-8), this will be
C     changed to build only the upper triangular part...

      JDOFE1=1

C     To calculate the matrix, we loop over the elements.

      DO IEL=1,TRIA(ONEL)
      
C       The outstanding feature with finite elements is: A basis
C       function for a DOF on one element has common support only
C       with the DOF's on the same element! E.g. for Q1:
C
C              #. . .#. . .#. . .#
C              .     .     .     .
C              .  *  .  *  .  *  .
C              #-----O-----O. . .#
C              |     |     |     .
C              |     | IEL |  *  .
C              #-----X-----O. . .#
C              |     |     |     .
C              |     |     |  *  .
C              #-----#-----#. . .#
C
C       --> On element IEL, the basis function at "X" only interacts
C           with the basis functions in "O". Elements in the 
C           neighbourhood ("*") have no support, therefore we only have
C           to collect all "O" DOF's.
C
C       Calculate the local DOF's into KDFL.
C       Calculate the global DOF's into KDFG.
C       The local DOF's are of course only in the range 1..IDFL,
C       but are probably unsorted. Each local DOF KDFL(I)=1..IDFL
C       corresponds to a global DOF KDFG(I).
      
        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       For the assembly of the global matrix, we use a "local"
C       approach. At first we build a "local" system matrix according
C       to the current element. This contains all additive
C       contributions of element IEL, which are later added at the
C       right positions to the elements in the global system matrix.
C
C       We have IDFL DOF's per element, therefore there are
C       IDFL*IDFL tupel of basis-/testfunctions (phi_i,psi_j) "active"
C       (i.e. have common support) on our current element, each giving
C       an additive contribution to the system matrix.
C
C       We build a quadratic IDFL*IDFL local matrix:
C       KENTRY(1..IDFL,1..IDFL) receives the position in the global
C         system matrix, where the corresponding value has to be 
C         added to.
C       (The corresponding contrbutions can be saved separately, 
C        but we directly add them to the global matrix in this 
C        approach.)
C
C       Loop through the local matrix to initialise it:

        DO JDOFE=1,IDFL
        
C         Row JDOFE of the local matrix corresponds 
C         to row=global DOF KDFG(JDOFE) in the global matrix.
C         This is the "X" in the above picture.
C         Save the absolute position in the DA-array of the diagonal 
C         into KENTRY:
        
          ILD=KLD(KDFG(JDOFE))
          KENTRY(JDOFE,JDOFE)=ILD
          
C         Now we loop through the other DOF's on the current element
C         (the "O"'s).
C         All these have common support with our current basis function
C         and will therefore give an additive value to the global
C         matrix.
C         In case we have a symmetric matrix, we set JDOFE1 from 1 to
C         JDOFE, which gives only the upper triangular part of the
C         matrix!
          
          IF (BSYMM) JDOFE1=JDOFE
          
C         JCOL0 is set to the position of the first entry in the row:
          
          JCOL0=ILD
          
          DO IDOFE=JDOFE1,IDFL
          
C           Ignore the case IDOFE=JDOFE - this gives the diagonal
C           entry which we already have ("X" = "O").
          
            IF (IDOFE.NE.JDOFE) THEN
            
C             Get the global DOF of the "O" which interacts with 
C             our "X".
            
              IDFG=KDFG(IDOFE)
              
C             Starting in JCOL0 (which points to the beginning of
C             the line initially), loop through the elements in
C             the row to find the position of column IDFG.
C             Jump out of the DO loop if we find the column.
              
              DO JCOL=JCOL0,NA
                IF (KCOL(JCOL).EQ.IDFG) GOTO 113
              END DO
113           CONTINUE

C             Because columns in the global matrix are sorted 
C             ascendingly (except for the diagonal element),
C             the next search can start after the column we just found.
              
              JCOL0=JCOL+1
              
C             Save the position of the matrix entry into the local
C             matrix.
              
              KENTRY(JDOFE,IDOFE)=JCOL
            
            END IF ! IDOFE=JDOFE
            
          END DO ! IDOFE
          
        END DO ! JDOFE

C       To calculate the matrix contributions, we have to evaluate
C       ELE on our current element. For this, we have to tell ELE
C       where we are:
C       - Store the vertex numbers of the current element to KVE
C       - Store the coordinates of the vertices of IEL to DX/DY in 
C         the ELE COMMON block
C       - Calculate the Jacobian of the mapping from the reference
C         element [-1,1]^2 to the real element into DJAC(1..2,1..2)
C         in the /ELEM/ COMMON block
C       - Calculate the Jacobian determinant of DJAC into DETJ in the
C         /ELEM/ COMMON block.
C       The last two points are necessary to correctly calculate
C       derivatives of the element!
C
C       So at first initialise KVE and DX/DY:

        DO IVE = 1, TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO
  
C       Now we want to calculate the Jacobian determinant.
C
C       We want to perform a transformation from the reference quadrilateral
C       [-1,1]x[-1,1] onto a "real" quadrilaterl:
C
C         (-1,1) ___ (1,1)         (x4,y4)   ___  (x3,y3)
C               |  |                        /  \
C               |__|           =>          /____\
C        (-1,-1)    (1,-1)         (x1,y1)        (x2,y2)   
C
C       By theory this can by done with a bilinear mapping, i.e. a mapping
C       Psi:R^2->R^2  of the form:
C
C            Psi (xi1) = ( a1 + a2*xi1 + a3*xi2 + a4*xi1*xi2 )
C                (xi2)   ( b1 + b2*xi1 + b3*xi2 + b4*xi1*xi2 )
C
C       Our special transformation has to map:
C
C        Psi(-1,-1) = (x1,y1)
C        Psi( 1,-1) = (x2,y2)
C        Psi( 1, 1) = (x3,y3)
C        Psi(-1, 1) = (x4,y4)
C
C       This gives the linear system:
C
C        a1 - a2 - a3 - a4 = x1       b1 - b2 - b3 - b4 = y1
C        a1 + a2 - a3 - a4 = x2       b1 + b2 - b3 - b4 = y2
C        a1 + a2 + a3 + a4 = x3       b1 + b2 + b3 + b4 = y3
C        a1 - a2 + a3 - a4 = x4       b1 - b2 + b3 - b4 = y4
C
C       Reorder this to calculate the ai:
C
C        a1 = 1/4 * ( x1 + x2 + x3 + x4)      b1 = 1/4 * ( y1 + y2 + y3 + y4)
C        a2 = 1/4 * (-x1 + x2 + x3 - x4)      b2 = 1/4 * (-y1 + y2 + y3 - y4)
C        a3 = 1/4 * (-x1 - x2 + x3 + x4)      b3 = 1/4 * (-y1 - y2 + y3 + y4)
C        a4 = 1/4 * ( x1 - x2 + x3 - x4)      b4 = 1/4 * ( y1 - y2 + y3 - y4)
C
C       The factors in the brackets in these equations are only dependent on 
C       the corners of the quadrilateral, not of the current point. So they
C       are constant for all points we want to map from the reference
C       element to the real one. We call them here "auxiliary Jacobian factors".
C       They are calculated in advance:
C
C         J1 = 1/2 * (-x1 - x2 + x3 - x4)
C         J2 = 1/2 * ( x1 - x2 + x3 - x4)
C
C         J3 = 1/2 * (-y1 + y2 - y3 + y4)
C         J4 = 1/2 * (-y1 + y2 + y3 - y4)

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
        
C       Using these, one can write:

C         a1  =  1/4 * ( x1 + x2 + x3 + x4)  =  1/2 * (x1 + x2 + J1)
C         a2  =  1/4 * (-x1 + x2 + x3 - x4)  =  1/2 * (x2 - x1 + J2)
C         a3  =  1/4 * (-x1 - x2 + x3 + x4)  =  1/2 * J1 
C         a4  =  1/4 * ( x1 - x2 + x3 - x4)  =  1/2 * J2
C
C         b1  =  1/4 * ( y1 + y2 + y3 + y4)  =  1/2 * (y1 + y3 + J3)
C         b2  =  1/4 * (-y1 + y2 + y3 - y4)  =  1/2 * J4
C         b3  =  1/4 * (-y1 - y2 + y3 + y4)  =  1/2 * (y3 - y1 - J4)
C         b4  =  1/4 * ( y1 - y2 + y3 - y4)  =  1/2 * J3
C
C       Now loop through all cubature points to build the local
C       matrix:

        DO ICUBP = 1,NCUBP

C         Get the coordinates of the cubature point on the reference
C         element from DXI in the /CUB/ COMMON block:

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)

C         Calculate the Jacobian of the mapping from the reference
C         element to the real element, based on:
C
C DPhi ( xi1 )  =  ( a2 + a4*xi2    a3 + a4*xi1 )
C      ( xi2 )     ( b2 + b4*xi2    b3 + b4*xi1 )
C
C               =  ( 1/2*(x2-x1+J2) + 1/2*J2*xi2              1/2*J1 + 1/2*J2*xi1 )
C                  (         1/2*J4 + 1/2*J3*xi2      1/2*(y3-y1-J4) + 1/2*J3*xi1 )

          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1

C         and from that its determinant:
C
C           det DPhi (xi1,xi2)  =  (DPhi[1,1]*DPhi[2,2] - DPhi[1,2]*DPhi[2,1]) (xi1,xi2)

          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)

C         This allows to calculate the current weighting factor in 
C         the cubature formula:

          OM=DOMEGA(ICUBP)*DETJ

C         Also calculate the real coordinates of the cubature point: 
C           (XX,YY) := Psi(XI1,XI2) 
C         We probably need it for the evaluation of non-constant
C         coefficients of the bilinear form.

          XX = 0.5D0*(DX(1)+DX(2)+DJ1)
     *        +0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *        +0.5D0*DJ1*XI2
     *        +0.5D0*DJ2*XI1*XI2
     
          YY = 0.5D0*(DY(1)+DY(3)+DJ3)
     *        +0.5D0*DJ4*XI1
     *        +0.5D0*(DY(3)-DY(1)-DJ4)*XI2
     *        -0.5D0*DJ3*XI1*XI2

C         Call the element subroutine to calculate the values
C         of the basis functions in the cubature point (XI1,XI2).
C
C         The variable ICUBP which we use in the loop tells the
C         element the number of the cubature point where we are
C         at the moment. Therefore, when we call ELE with IPAR=-1,
C         ELE can decide on its own whether to calculate the
C         values on the basis functions again or to use the
C         precalculated ones.

          CALL ELE(XI1,XI2,-3)
          IF (IER.LT.0) GOTO 99999

C         Sum up over all pairs of multiindices.
C         Loop though all matrix blocks:

          BFIRST=.TRUE.
          
          DO IBLOC=1,NBLOC
          
C           Loop through all KABN prescribed combinations of the
C           derivatives of the basis functions.

            DO IALBET=1,KABN(IBLOC)
            
C             Get from KAB the type of the derivatives for the 
C             test and trial functions. The summand we calculate
C             here will be:
C
C             int_... ( phi_i )_IA  *  ( phi_j )_IB
C
C             -> Ix=0: function value, 
C                  =1: first derivative, 
C                  =2: 2nd derivative,...
              
              IA=KAB(1,IALBET,IBLOC)
              IB=KAB(2,IALBET,IBLOC)
              
C             Get the coefficient of the basis function - either
C             from COEFF (for non-constant coefficients) or take
C             the precalculated one from COECON.
C             Multiply it with the cubature weight OM.
              
              IF (.NOT.BCON(IBLOC)) THEN
                AUX=COEFF(XX,YY,IA,IB,IBLOC,BFIRST,TRIA,
     *              IPARAM,DPARAM)*OM
              ELSE
                AUX=COECON(IA,IB,IBLOC)*OM
              ENDIF
              
C             Now loop through all possible combinations of DOF's
C             in the current cubature point. The outer loop
C             loops through the "X" in the above picture:

              DO  JDOFE=1,IDFL
              
C               Get the value of the basis function 
C               phi_i (our "X") in the cubature point:

                DB=DBAS(KDFL(JDOFE),IA)
                
C               Perform an inner loop through the other cubature
C               points (the "O"'s). 
C               For symmetric matrices, change JDOFE1 to JDOFE to 
C               handle only the upper triangular part.
                
                IF (BSYMM) JDOFE1=JDOFE
                
                DO IDOFE=JDOFE1,IDFL

C                 Get the value of the basis function 
C                 phi_j (our "O") in the cubature point. This is
C                 DBAS(KDFL(IDOFE),IB).
C                 Them multiply:
C                    DB * DBAS(..) * AUX
C                 ~= psi_j * phi_j * coefficient * cub.weight
C                 Summing this up gives the integral, so the contribution
C                 to the global matrix. 
C
C                 Simply summing up DB * DBAS(..) * AUX would give
C                 the coefficient of the local matrix. As everything
C                 is only summed up here, we directly incorporate
C                 this value into the global matrix at the position 
C                 given by KENTRY.

                  JCOLB=KENTRY(JDOFE,IDOFE)+KOFF(IBLOC)
                  DA(JCOLB)=DA(JCOLB)+DB*DBAS(KDFL(IDOFE),IB)*AUX
                  
                END DO ! IDOFE
                
              END DO ! JDOFE
              
            END DO ! IALBET

            BFIRST=.FALSE.
            
          END DO ! IBLOC

        END DO ! ICUBP
        
      END DO ! IEL

C     Assembly finished.
C
C     For symmetric matrices we have to do some postprocessing:
C     For ISYMM=1 we have matrix structure 8, saving only the upper
C     triangular part. Nothing has to be done.
C
C     If ISYMM>=2, we have a symmetric structure-7 matrix!
C     In this case, copy the matrix entries of the upper triangulat 
C     part to the lower triangular part!
C     Remark: For matrices in structure 8, nothing will happen here
C     as only elements of the upper triangular part in the matrix
C     exist - so there's no lower triangular part to copy anything
C     to...

      IF (ISYMM.GT.1) THEN
      
C       Loop over the matrix blocks:
      
        DO IBLOC=1,NBLOC

C         Loop through all rows in the matrix

          DO IROW=1,NEQ
          
C           Loop through all offdiagonal entries in row IROW:
          
            DO ICOL=KLD(IROW)+1,KLD(IROW+1)-1
            
C             Only take care of elements below the diagonal.
C             As soon as we reach an element of the upper
C             triangular matrix, cancel the loop.
C             (For matrices in structure 8, this will always
C              be the case!)
            
              J1=KCOL(ICOL)
              IF (J1.GE.IROW) GOTO 401
              
C             Element (IROW,J1) below the diagonal found. This should
C             be the same as (J1,IROW)! To perform 
C               A(IROW,J1) := A(J1,IROW)
C             we have to find element (J1,IROW) in DA!.
C             Look into line J1 to find the element:
              
              DO ICOL1=KLD(J1+1)-1,KLD(J1),-1
              
C               If we find element (J1,IROW),
              
                IF (KCOL(ICOL1).EQ.IROW) THEN
                
C                 Copy the element:   A(IROW,J1) := A(J1,IROW)
C
C                 Of course take care of the starting indices of the 
C                 different matrix blocks in DA!
                
                  ICOLB=ICOL+KOFF(IBLOC)
                  ICOL1B=ICOL1+KOFF(IBLOC)

                  DA(ICOLB)=DA(ICOL1B)
                  
C                 BREAK the loop
                  
                  GOTO 402
                  
                END IF ! KCOL(ICOL1) = IROW
              END DO ! ICOL1
              
402         END DO ! ICOL
401       END DO ! IROW
        END DO ! IBLOC

      ENDIF ! ISYMM > 1

99999 END

************************************************************************
* Generate matrix in matrix structure 7/8
*
* Extended calling convention; quadrilateral version
*
* This routine creates the entries of a matrix in matrix structure
* 7 or 8 of a bilinear form a(u,v). ELE is used for both, trial and
* test functions. ELE is assumed to be nonparametric.
*
* In:
*   DA     : array [1..*] of double
*            Matrix which entries are to be generated.
*            Can be <> 0.0, in this case the corresponding new entries
*            are added to the existing ones.
*            Multiple matrix blocks can be calculated at once (providing
*            the matrices have the same structure). The array KOFF
*            describes the starting addresses of the matrices in DA.
*   KCOL   : array [1..NA] of integer
*            Column structure
*   KLD    : array [1..NEQ+1] of integer
*            Row structure
*   NA     : Number of entries in each matrix block
*   NEQ    : Dimension/Number of rows in the matrix
*   NBLOC  : Number of independent diagonal blocks in the matrix
*   KOFF   : array [1..NBLOC] of integer
*            Start indices of the matrices in DA.
*            Matrix I is expected in DA(KOFF(I)..KOFF(I)+NA)
*   DCORVG,
*   KVERT,
*   KMID   : Geometry information arrays of a quadrilateral mesh;
*            must be the arrays given by TRIA, but passed separately
*            for performance reasons.
*   TRIA   : array [1..SZTRIA] of integer
*            Triangulation structure of the underlying mesh
*   ELE    : Finite Element callback routine.
*            Must be quadrilateral and fit to
*            the triangulation in TRIA!
*   COEFF  : Coefficient function for the coefficients in the 
*            bilinear form.
*            DOUBLE PRECISION FUNCTION COEFF(XX,YY,IA,IB,IBLOC,
*                   BFIRST,TRIA,IPARAM,DPARAM)
*   BCON   : array [1..NBLOC] of boolean
*            BCON(I)=true,  if the matrix block I has constant 
*                           coefficients
*            BCON(I)=false, otherwise
*   COECON : array [COECON(NNDER,NNDER,NBLOC] of double
*            Auxiliary array, saves the constant coefficients during
*            computation
*   KAB    : array [2,NNAB,KABN] of integer
*            Descriptors for the bilinear form
*   KABN   : Number of additive terms in the bilinear form
*   ICUB   : Cubature formula identifier for CB2Q
*   ISYMM  : =0: don't respect symmetry, create matrix in structure 7
*            =1: matrix is symmetric; create matrix in structure 8
*   IPARAM : array [1..*] of integer
*            User defined parameter block, passed to COEFF
*   DPARAM : array [1..*] of double
*            User defined parameter block, passed to COEFF
* Out:
*   The entries of the matrix are added to DA.
*
* Out (COMMON block):
*   IER    = -116 for wrong value in array KAB
************************************************************************

      SUBROUTINE AB7NX(DA,KCOL,KLD,NA,NEQ,NBLOC,KOFF,KVERT,KMID,
     *                DCORVG,TRIA,ELE,COEFF,BCON,COECON,
     *                KAB,KABN,ICUB,ISYMM,IPARAM,DPARAM)

      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'ccub.inc'
      
C parameters
      
      INTEGER KCOL(*),KVERT(NNVE,*),KMID(NNVE,*),KOFF(*)
      INTEGER KLD(*),KAB(2,NNAB,*),KABN(*),TRIA(SZTRIA)
      INTEGER KENTRY(NNBAS,NNBAS)
      DOUBLE PRECISION DA(*),DCORVG(2,*),COECON(NNDER,NNDER,*)
      LOGICAL BCON(*)
      INTEGER NEQ,NA,ICUB,ISYMM,NBLOC

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      DOUBLE PRECISION COEFF
      EXTERNAL ELE,COEFF

C     externals

      INTEGER NDFL
      EXTERNAL NDFL

C     local variables

      LOGICAL BSYMM,BFIRST,BCON0
      INTEGER I,J,IBLOC,I1,IELTYP,IA,IB,ILD,ICOL,JCOL
      INTEGER IDOFE,JDOFE,JDOFE1,JCOL0,IDFG
      INTEGER IVE,JP,IALBET,JCOLB,IROW,J1,ICOL1,ICOLB,ICOL1B
      DOUBLE PRECISION AUX,DJ1,DJ2,DJ3,DJ4,XI1,XI2,XX,YY,OM,DB
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL

      SUB='ABM7'
      IF (ICHECK.GE.997) CALL OTRC('ABM7  ','12/08/94')

C     Preparation - evaluation of parameters

      IER=0
      
C     Create symmetric matrix in matrix structure 8?      
      
      BSYMM=ISYMM.GE.1
      
C     Which derivatives of basis functions are needed?
C     Check the descriptors of the bilinear form and set BDER
C     according to these.


      DO I = 1,NNDER
        BDER(I)=.FALSE.
      END DO
      
      DO IBLOC=1,NBLOC
        DO I=1,KABN(IBLOC)
          DO J=1,2
          
C           The desriptor KAB gives directly the derivative
C           which is to be computed!

            I1=KAB(J,I,IBLOC)
            
            IF (I1.LE.0.OR.I1.GT.NNDER) THEN
            
C             Oops, invalid descriptor
            
              WRITE (CPARAM,'(I15)') IBLOC
              CALL WERR(-116,'AB07  ')
              GOTO 99999
            ENDIF
            
            BDER(I1)=.TRUE.
          END DO
        END DO
      END DO
      
C     Dummy call of ELE to get the element identifier
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      
C     Get the number of local DOF's:
      
      IDFL=NDFL(IELTYP)
      
      
C     Initialise the cubature formula
      
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999

C     Dummy call of COEFF for nonlinear problems.
C     COEFF must set BDER(IDER)=.TRUE. in the /ELEM/ COMMON block
C     derivative IDER is needed from the element.
C     This allows to tell the element to calculate additional
C     derivarives to the ones defined by KAB.

      BFIRST=.TRUE.
      AUX=COEFF(0D0,0D0,-1,-1,0,BFIRST,TRIA,IPARAM,DPARAM)
      
C     Check all matrix blocks if their bilinear form has constant
C     coefficients. In that case, get the coefficient and save it
C     in COECON.
C     BCON0 is set to TRUE if all blocks have constant coefficients.      
      
      BCON0=.TRUE.
      DO IBLOC=1,NBLOC
        IF (BCON(IBLOC)) THEN
          DO I=1,KABN(IBLOC)
            IA=KAB(1,I,IBLOC)
            IB=KAB(2,I,IBLOC)
            COECON(IA,IB,IBLOC)=COEFF(0D0,0D0,IA,IB,IBLOC,BFIRST,TRIA,
     *                                IPARAM,DPARAM)    
          END DO
        ELSE
          BCON0=.FALSE.
        ENDIF
      END DO
     
C     Now let's calculate the matrix - storage technique 7 or 8
C
C     For the loop through the DOF's later, set JDOFE1=1.
C     For symmetric matrices (structure-8), this will be
C     changed to build only the upper triangular part...

      JDOFE1=1

C     To calculate the matrix, we loop over the elements.

      DO IEL=1,TRIA(ONEL)
      
C       The outstanding feature with finite elements is: A basis
C       function for a DOF on one element has common support only
C       with the DOF's on the same element! E.g. for Q1:
C
C              #. . .#. . .#. . .#
C              .     .     .     .
C              .  *  .  *  .  *  .
C              #-----O-----O. . .#
C              |     |     |     .
C              |     | IEL |  *  .
C              #-----X-----O. . .#
C              |     |     |     .
C              |     |     |  *  .
C              #-----#-----#. . .#
C
C       --> On element IEL, the basis function at "X" only interacts
C           with the basis functions in "O". Elements in the 
C           neighbourhood ("*") have no support, therefore we only have
C           to collect all "O" DOF's.
C
C       Calculate the local DOF's into KDFL.
C       Calculate the global DOF's into KDFG.
C       The local DOF's are of course only in the range 1..IDFL,
C       but are probably unsorted. Each local DOF KDFL(I)=1..IDFL
C       corresponds to a global DOF KDFG(I).
      
        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       For the assembly of the global matrix, we use a "local"
C       approach. At first we build a "local" system matrix according
C       to the current element. This contains all additive
C       contributions of element IEL, which are later added at the
C       right positions to the elements in the global system matrix.
C
C       We have IDFL DOF's per element, therefore there are
C       IDFL*IDFL tupel of basis-/testfunctions (phi_i,psi_j) "active"
C       (i.e. have common support) on our current element, each giving
C       an additive contribution to the system matrix.
C
C       We build a quadratic IDFL*IDFL local matrix:
C       KENTRY(1..IDFL,1..IDFL) receives the position in the global
C         system matrix, where the corresponding value has to be 
C         added to.
C       (The corresponding contrbutions can be saved separately, 
C        but we directly add them to the global matrix in this 
C        approach.)
C
C       Loop through the local matrix to initialise it:

        DO JDOFE=1,IDFL
        
C         Row JDOFE of the local matrix corresponds 
C         to row=global DOF KDFG(JDOFE) in the global matrix.
C         This is the "X" in the above picture.
C         Save the absolute position in the DA-array of the diagonal 
C         into KENTRY:
        
          ILD=KLD(KDFG(JDOFE))
          KENTRY(JDOFE,JDOFE)=ILD
          
C         Now we loop through the other DOF's on the current element
C         (the "O"'s).
C         All these have common support with our current basis function
C         and will therefore give an additive value to the global
C         matrix.
C         In case we have a symmetric matrix, we set JDOFE1 from 1 to
C         JDOFE, which gives only the upper triangular part of the
C         matrix!
          
          IF (BSYMM) JDOFE1=JDOFE
          
C         JCOL0 is set to the position of the first entry in the row:
          
          JCOL0=ILD
          
          DO IDOFE=JDOFE1,IDFL
          
C           Ignore the case IDOFE=JDOFE - this gives the diagonal
C           entry which we already have ("X" = "O").
          
            IF (IDOFE.NE.JDOFE) THEN
            
C             Get the global DOF of the "O" which interacts with 
C             our "X".
            
              IDFG=KDFG(IDOFE)
              
C             Starting in JCOL0 (which points to the beginning of
C             the line initially), loop through the elements in
C             the row to find the position of column IDFG.
C             Jump out of the DO loop if we find the column.
              
              DO JCOL=JCOL0,NA
                IF (KCOL(JCOL).EQ.IDFG) GOTO 113
              END DO
113           CONTINUE

C             Because columns in the global matrix are sorted 
C             ascendingly (except for the diagonal element),
C             the next search can start after the column we just found.
              
              JCOL0=JCOL+1
              
C             Save the position of the matrix entry into the local
C             matrix.
              
              KENTRY(JDOFE,IDOFE)=JCOL
            
            END IF ! IDOFE=JDOFE
            
          END DO ! IDOFE
          
        END DO ! JDOFE

C       To calculate the matrix contributions, we have to evaluate
C       ELE on our current element. For this, we have to tell ELE
C       where we are:
C       - Store the vertex numbers of the current element to KVE
C       - Store the coordinates of the vertices of IEL to DX/DY in 
C         the ELE COMMON block
C       - Calculate the Jacobian of the mapping from the reference
C         element [-1,1]^2 to the real element into DJAC(1..2,1..2)
C         in the /ELEM/ COMMON block
C       - Calculate the Jacobian determinant of DJAC into DETJ in the
C         /ELEM/ COMMON block.
C       The last two points are necessary to correctly calculate
C       derivatives of the element!
C
C       So at first initialise KVE and DX/DY:

        DO IVE = 1, TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO
  
C       Now we want to calculate the Jacobian determinant.
C
C       We want to perform a transformation from the reference quadrilateral
C       [-1,1]x[-1,1] onto a "real" quadrilaterl:
C
C         (-1,1) ___ (1,1)         (x4,y4)   ___  (x3,y3)
C               |  |                        /  \
C               |__|           =>          /____\
C        (-1,-1)    (1,-1)         (x1,y1)        (x2,y2)   
C
C       By theory this can by done with a bilinear mapping, i.e. a mapping
C       Psi:R^2->R^2  of the form:
C
C            Psi (xi1) = ( a1 + a2*xi1 + a3*xi2 + a4*xi1*xi2 )
C                (xi2)   ( b1 + b2*xi1 + b3*xi2 + b4*xi1*xi2 )
C
C       Our special transformation has to map:
C
C        Psi(-1,-1) = (x1,y1)
C        Psi( 1,-1) = (x2,y2)
C        Psi( 1, 1) = (x3,y3)
C        Psi(-1, 1) = (x4,y4)
C
C       This gives the linear system:
C
C        a1 - a2 - a3 - a4 = x1       b1 - b2 - b3 - b4 = y1
C        a1 + a2 - a3 - a4 = x2       b1 + b2 - b3 - b4 = y2
C        a1 + a2 + a3 + a4 = x3       b1 + b2 + b3 + b4 = y3
C        a1 - a2 + a3 - a4 = x4       b1 - b2 + b3 - b4 = y4
C
C       Reorder this to calculate the ai:
C
C        a1 = 1/4 * ( x1 + x2 + x3 + x4)      b1 = 1/4 * ( y1 + y2 + y3 + y4)
C        a2 = 1/4 * (-x1 + x2 + x3 - x4)      b2 = 1/4 * (-y1 + y2 + y3 - y4)
C        a3 = 1/4 * (-x1 - x2 + x3 + x4)      b3 = 1/4 * (-y1 - y2 + y3 + y4)
C        a4 = 1/4 * ( x1 - x2 + x3 - x4)      b4 = 1/4 * ( y1 - y2 + y3 - y4)
C
C       The factors in the brackets in these equations are only dependent on 
C       the corners of the quadrilateral, not of the current point. So they
C       are constant for all points we want to map from the reference
C       element to the real one. We call them here "auxiliary Jacobian factors".
C       They are calculated in advance:
C
C         J1 = 1/2 * (-x1 - x2 + x3 - x4)
C         J2 = 1/2 * ( x1 - x2 + x3 - x4)
C
C         J3 = 1/2 * (-y1 + y2 - y3 + y4)
C         J4 = 1/2 * (-y1 + y2 + y3 - y4)

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
        
C       Using these, one can write:

C         a1  =  1/4 * ( x1 + x2 + x3 + x4)  =  1/2 * (x1 + x2 + J1)
C         a2  =  1/4 * (-x1 + x2 + x3 - x4)  =  1/2 * (x2 - x1 + J2)
C         a3  =  1/4 * (-x1 - x2 + x3 + x4)  =  1/2 * J1 
C         a4  =  1/4 * ( x1 - x2 + x3 - x4)  =  1/2 * J2
C
C         b1  =  1/4 * ( y1 + y2 + y3 + y4)  =  1/2 * (y1 + y3 + J3)
C         b2  =  1/4 * (-y1 + y2 + y3 - y4)  =  1/2 * J4
C         b3  =  1/4 * (-y1 - y2 + y3 + y4)  =  1/2 * (y3 - y1 - J4)
C         b4  =  1/4 * ( y1 - y2 + y3 - y4)  =  1/2 * J3
C
C       Dummy call - ELE may save arithmetic operations.
C       ELE can precalculate the values in the cubature points if
C       it supports that...
C
C       Remark: In contrast to parametric elements, this is inside the 
C       loop here!

        CALL ELE(0D0,0D0,-2)
        IF (IER.LT.0) GOTO 99999

C
C       Now loop through all cubature points to build the local
C       matrix:

        DO ICUBP = 1,NCUBP

C         Get the coordinates of the cubature point on the reference
C         element from DXI in the /CUB/ COMMON block:

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)

C         Calculate the Jacobian of the mapping from the reference
C         element to the real element, based on:
C
C DPhi ( xi1 )  =  ( a2 + a4*xi2    a3 + a4*xi1 )
C      ( xi2 )     ( b2 + b4*xi2    b3 + b4*xi1 )
C
C               =  ( 1/2*(x2-x1+J2) + 1/2*J2*xi2              1/2*J1 + 1/2*J2*xi1 )
C                  (         1/2*J4 + 1/2*J3*xi2      1/2*(y3-y1-J4) + 1/2*J3*xi1 )

          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1

C         and from that its determinant:
C
C           det DPhi (xi1,xi2)  =  (DPhi[1,1]*DPhi[2,2] - DPhi[1,2]*DPhi[2,1]) (xi1,xi2)

          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)

C         This allows to calculate the current weighting factor in 
C         the cubature formula:

          OM=DOMEGA(ICUBP)*DETJ

C         Also calculate the real coordinates of the cubature point: 
C           (XX,YY) := Psi(XI1,XI2) 
C         We probably need it for the evaluation of non-constant
C         coefficients of the bilinear form.

          XX = 0.5D0*(DX(1)+DX(2)+DJ1)
     *        +0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *        +0.5D0*DJ1*XI2
     *        +0.5D0*DJ2*XI1*XI2
     
          YY = 0.5D0*(DY(1)+DY(3)+DJ3)
     *        +0.5D0*DJ4*XI1
     *        +0.5D0*(DY(3)-DY(1)-DJ4)*XI2
     *        -0.5D0*DJ3*XI1*XI2

C         Call the element subroutine to calculate the values
C         of the basis functions in the cubature point 
C         in real coordinates (XX,YY).
C
C         The variable ICUBP which we use in the loop tells the
C         element the number of the cubature point where we are
C         at the moment. Therefore, when we call ELE with IPAR=-1,
C         ELE can decide on its own whether to calculate the
C         values on the basis functions again or to use the
C         precalculated ones.

          CALL ELE(XX,YY,-3)
          IF (IER.LT.0) GOTO 99999

C         Sum up over all pairs of multiindices.
C         Loop though all matrix blocks:

          BFIRST=.TRUE.
          
          DO IBLOC=1,NBLOC
          
C           Loop through all KABN prescribed combinations of the
C           derivatives of the basis functions.

            DO IALBET=1,KABN(IBLOC)
            
C             Get from KAB the type of the derivatives for the 
C             test and trial functions. The summand we calculate
C             here will be:
C
C             int_... ( phi_i )_IA  *  ( phi_j )_IB
C
C             -> Ix=0: function value, 
C                  =1: first derivative, 
C                  =2: 2nd derivative,...
              
              IA=KAB(1,IALBET,IBLOC)
              IB=KAB(2,IALBET,IBLOC)
              
C             Get the coefficient of the basis function - either
C             from COEFF (for non-constant coefficients) or take
C             the precalculated one from COECON.
C             Multiply it with the cubature weight OM.
              
              IF (.NOT.BCON(IBLOC)) THEN
                AUX=COEFF(XX,YY,IA,IB,IBLOC,BFIRST,TRIA,
     *              IPARAM,DPARAM)*OM
              ELSE
                AUX=COECON(IA,IB,IBLOC)*OM
              ENDIF
              
C             Now loop through all possible combinations of DOF's
C             in the current cubature point. The outer loop
C             loops through the "X" in the above picture:

              DO  JDOFE=1,IDFL
              
C               Get the value of the basis function 
C               phi_i (our "X") in the cubature point:

                DB=DBAS(KDFL(JDOFE),IA)
                
C               Perform an inner loop through the other cubature
C               points (the "O"'s). 
C               For symmetric matrices, change JDOFE1 to JDOFE to 
C               handle only the upper triangular part.
                
                IF (BSYMM) JDOFE1=JDOFE
                
                DO IDOFE=JDOFE1,IDFL

C                 Get the value of the basis function 
C                 phi_j (our "O") in the cubature point. This is
C                 DBAS(KDFL(IDOFE),IB).
C                 Them multiply:
C                    DB * DBAS(..) * AUX
C                 ~= phi_i * phi_j * coefficient * cub.weight
C                 Summing this up gives the integral, so the contribution
C                 to the global matrix. 
C
C                 Simply summing up DB * DBAS(..) * AUX would give
C                 the coefficient of the local matrix. As everything
C                 is only summed up here, we directly incorporate
C                 this value into the global matrix at the position 
C                 given by KENTRY.

                  JCOLB=KENTRY(JDOFE,IDOFE)+KOFF(IBLOC)
                  DA(JCOLB)=DA(JCOLB)+DB*DBAS(KDFL(IDOFE),IB)*AUX
                  
                END DO ! IDOFE
                
              END DO ! JDOFE
              
            END DO ! IALBET

            BFIRST=.FALSE.
            
          END DO ! IBLOC

        END DO ! ICUBP
        
      END DO ! IEL

C     Assembly finished.
C
C     For symmetric matrices we have to do some postprocessing:
C     For ISYMM=1 we have matrix structure 8, saving only the upper
C     triangular part. Nothing has to be done.
C
C     If ISYMM>=2, we have a symmetric structure-7 matrix!
C     In this case, copy the matrix entries of the upper triangulat 
C     part to the lower triangular part!
C     Remark: For matrices in structure 8, nothing will happen here
C     as only elements of the upper triangular part in the matrix
C     exist - so there's no lower triangular part to copy anything
C     to...

      IF (ISYMM.GT.1) THEN
      
C       Loop over the matrix blocks:
      
        DO IBLOC=1,NBLOC

C         Loop through all rows in the matrix

          DO IROW=1,NEQ
          
C           Loop through all offdiagonal entries in row IROW:
          
            DO ICOL=KLD(IROW)+1,KLD(IROW+1)-1
            
C             Only take care of elements below the diagonal.
C             As soon as we reach an element of the upper
C             triangular matrix, cancel the loop.
C             (For matrices in structure 8, this will always
C              be the case!)
            
              J1=KCOL(ICOL)
              IF (J1.GE.IROW) GOTO 401
              
C             Element (IROW,J1) below the diagonal found. This should
C             be the same as (J1,IROW)! To perform 
C               A(IROW,J1) := A(J1,IROW)
C             we have to find element (J1,IROW) in DA!.
C             Look into line J1 to find the element:
              
              DO ICOL1=KLD(J1+1)-1,KLD(J1),-1
              
C               If we find element (J1,IROW),
              
                IF (KCOL(ICOL1).EQ.IROW) THEN
                
C                 Copy the element:   A(IROW,J1) := A(J1,IROW)
C
C                 Of course take care of the starting indices of the 
C                 different matrix blocks in DA!
                
                  ICOLB=ICOL+KOFF(IBLOC)
                  ICOL1B=ICOL1+KOFF(IBLOC)

                  DA(ICOLB)=DA(ICOL1B)
                  
C                 BREAK the loop
                  
                  GOTO 402
                  
                END IF ! KCOL(ICOL1) = IROW
              END DO ! ICOL1
              
402         END DO ! ICOL
401       END DO ! IROW
        END DO ! IBLOC

      ENDIF ! ISYMM > 1

99999 END

************************************************************************
* Generate matrix in matrix structure 9, wrapper routine
*
* Extended calling convention; quadrilateral version
*
* This routine calls AB7X or ABN7X to generate a the structure of 
* a matrix.
*
* In:
*   LA     : array [1..NBLOC] of integer
*            = 0: Create a new matrix
*            <>0: Handle to an existing matrix which should be
*                 rebuild or updated
*   LCOL   : Handle to the KCOL array; must be the same
*            for all matrix blocks
*   LLD    : Handle to the KLD array; must be the same
*            for all matrix blocks
*   NA     : Number of entries in each matrix block
*   NBLOC  : Number of independent diagonal blocks in the matrix
*   ICLEAR : =0: Add the matrix to an existing one
*            =1: Clear the matrix before building
*   TRIA   : array [1..SZTRIA] of integer
*            Triangulation structure of the underlying mesh
*   ELE1   : Finite Element callback routine of the matrix block
*            which is to generate (e.g. velocity).
*            Must be quadrilateral, corresponding to
*            the triangulation in TRIA!
*   ELE2   : Finite Element callback routine for elements that
*            are coupled to the matrix block which is to generate
*            (e.g. pressure).
*            Must be quadrilateral, corresponding to
*            the triangulation in TRIA!
*   ELE3   : A 3rd Finite Element callback routine which allows
*            to add a nonlinear term consisting of 3 multipliers
*            in the integral form:
*              int_... ... * phi_i * psi_j * eta_k
*                            ^^^^^   ^^^^^   ^^^^^
*                            ELE1    ELE2    ELE3
*            If no nonlinear term is to be added, ELE3 can be set
*            to ELE1 or ELE2. Whether or not nonlinearity is
*            added is decided by COEFF!
*            Must be quadrilateral, corresponding to
*            the triangulation in TRIA!
*   COEFF  : Coefficient function for the coefficients in the 
*            bilinear form.
*            DOUBLE PRECISION FUNCTION COEFF(XX,YY,IA,IB,IBLOC,
*                   BFIRST,TRIA,IPARAM,DPARAM)
*   KAB    : array [2,NNAB,KABN] of integer
*            Descriptors for the bilinear form
*   KABN   : Number of additive terms in the bilinear form
*   ICUB   : Cubature formula identifier for CB2Q
*   ARR    : string[6]
*            Name of the matrix array LA - for output purposes.
*   IPARAM : array [1..*] of integer
*            User defined parameter block, passed to COEFF
*   DPARAM : array [1..*] of double
*            User defined parameter block, passed to COEFF
* Out:
*   If LA=0: LA=handles to the matrices
*   The entries the matrix-array corresponding to LA are filled
*   with data according to the parameters.
*
* Out (COMMON blocks):
*   IER = -114 Type of at least one array is not double prec.  
*         -115 Length of at least one array is < NA
************************************************************************

      SUBROUTINE XAB09X(LA,LCOL,LLD,NA,NBLOC,ICLEAR,TRIA,
     *                 ELE1,ELE2,ELE3,COEFF,BCON,KAB,KABN,
     *                 ICUB,ARR,IPARAM,DPARAM)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'stria.inc'

C parameters

      INTEGER LA(*),KAB(2,NNAB,*),KABN(*)
      LOGICAL BCON(*)

      CHARACTER ARR*12
      DIMENSION ARR(*)
      
      INTEGER LCOL, LLD, NA, NBLOC, ICLEAR, ICUB
      INTEGER TRIA(SZTRIA)
      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      DOUBLE PRECISION COEFF
      EXTERNAL COEFF,ELE1,ELE2,ELE3

C local variables
      
      INTEGER IBLOC, ITYPE, ILEN, LOFF, LOECON
      INTEGER KVERT,KMID,KCORVG

      SUB='XAB09'
      IF (ICHECK.GE.997) CALL OTRC('XAB09 ','12/04/90')
      IER=0

C     Loop through the NBLOC independent blocks of the global
C     matrix:
C         MATRIX =  [ BLOCK1                       ]
C                   [         BLOCK2               ]
C                   [                 ...          ]
C                   [                      BLOCKn  ]
C     All blocks have the same column/row structure.

      DO IBLOC=1,NBLOC
      
C       Create a new matrix block? Check if LA is given:
      
        IF (LA(IBLOC).EQ.0) THEN

C         Yes! Allocate a new matrix.
        
          CALL ZNEW(NA,1,LA(IBLOC),ARR(IBLOC))
          IF (IER.NE.0) GOTO 99999
          
        ELSE
          
C         A matrix is given; check that it has the right dimension
C         and the right data type!
          
          CALL ZTYPE(LA(IBLOC),ITYPE)
          IF (ITYPE.NE.1) THEN
            WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
            CALL WERR(-114,'XAB09 ')
            GOTO 99999
          ENDIF
          
          CALL ZLEN(LA(IBLOC),ILEN)
          IF (ILEN.LT.NA) THEN
            WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
            CALL WERR(-115,'XAB09 ')
            GOTO 99999
          ENDIF
          
C         Clear the array if ICLEAR=1
          
          IF (ICLEAR.EQ.1) CALL ZCLEAR(LA(IBLOC),ARR(IBLOC))
          
        END IF
        
      END DO

C     Create an auxiliary array KOFF(NBLOC), which receives
C     for every matrix the starting address of the matrix
C     block in DWORK. 

      CALL ZNEW(NBLOC,-3,LOFF,'KOFF  ')
      IF (IER.NE.0) GOTO 99999

      DO IBLOC=1,NBLOC
        KWORK(L(LOFF)+IBLOC-1)=L(LA(IBLOC))-1
      END DO
      
C     Allocate an auxiliary array COECON.
C     In case of a bilinear form with constant coefficients, the
C     coefficients are determined in advance by ABx7X and stored here.
      
      CALL ZNEW(NBLOC*NNDER**2,1,LOECON,'COECON')
      IF (IER.NE.0) GOTO 99999

C     Call AB09X to calculate the matrix

      KVERT  = L(TRIA(OLVERT))
      KMID   = L(TRIA(OLMID))
      KCORVG = L(TRIA(OLCORVG))

      CALL AB09X(DWORK(1),KWORK(L(LCOL)),KWORK(L(LLD)),NA,NBLOC,
     *          KWORK(L(LOFF)),KWORK(KVERT),KWORK(KMID),
     *          DWORK(KCORVG),TRIA,ELE1,ELE2,
     *          ELE3,COEFF,BCON,DWORK(L(LOECON)),KAB,KABN,ICUB,
     *          IPARAM,DPARAM)
      IF (IER.NE.0) GOTO 99999

C     Release auxiliary arrays.

      CALL ZDISP(0,LOECON,'COECON')
      CALL ZDISP(0,LOFF,'KOFF  ')

99999 END

************************************************************************
* Generate matrix in matrix structure 9
*
* Extended calling convention; quadrilateral version
*
* This routine creates the entries of a matrix in matrix structure 9.
* ELE1-ELE3 is assumed to be parametric.
*
* In:
*   DA     : array [1..*] of double
*            Matrix which entries are to be generated.
*            Can be <> 0.0, in this case the corresponding new entries
*            are added to the existing ones.
*            Multiple matrix blocks can be calculated at once (providing
*            the matrices have the same structure). The array KOFF
*            describes the starting addresses of the matrices in DA.
*   KCOL   : array [1..NA] of integer
*            Column structure
*   KLD    : array [1..NEQ+1] of integer
*            Row structure
*   NA     : Number of entries in each matrix block
*   NBLOC  : Number of independent diagonal blocks in the matrix
*   KOFF   : array [1..NBLOC] of integer
*            Start indices of the matrices in DA.
*            Matrix I is expected in DA(KOFF(I)..KOFF(I)+NA)
*   DCORVG,
*   KVERT,
*   KMID   : Geometry information arrays of a quadrilateral mesh;
*            must be the arrays given by TRIA, but passed separately
*            for performance reasons.
*   TRIA   : array [1..SZTRIA] of integer
*            Triangulation structure of the underlying mesh
*   ELE1   : Finite Element callback routine of the matrix block
*            which is to generate (e.g. velocity).
*            Must be quadrilateral, corresponding to
*            the triangulation in TRIA!
*   ELE2   : Finite Element callback routine for elements that
*            are coupled to the matrix block which is to generate
*            (e.g. pressure).
*            Must be quadrilateral, corresponding to
*            the triangulation in TRIA!
*   ELE3   : A 3rd Finite Element callback routine which allows
*            to add a nonlinear term consisting of 3 multipliers
*            in the integral form:
*              int_... ... * phi_i * psi_j * eta_k
*                            ^^^^^   ^^^^^   ^^^^^
*                            ELE1    ELE2    ELE3
*            If no nonlinear term is to be added, ELE3 can be set
*            to ELE1 or ELE2. Whether or not nonlinearity is
*            added is decided by COEFF!
*            Must be quadrilateral, corresponding to
*            the triangulation in TRIA!
*            WARNING: Nonlinearity not yet implemented!
*   COEFF  : Coefficient function for the coefficients in the 
*            bilinear form.
*            DOUBLE PRECISION FUNCTION COEFF(XX,YY,IA,IB,IBLOC,
*                   BFIRST,TRIA,IPARAM,DPARAM)
*   BCON   : array [1..NBLOC] of boolean
*            BCON(I)=true,  if the matrix block I has constant 
*                           coefficients
*            BCON(I)=false, otherwise
*   COECON : array [COECON(NNDER,NNDER,NBLOC] of double
*            Auxiliary array, saves the constant coefficients during
*            computation
*   KAB    : array [2,NNAB,KABN] of integer
*            Descriptors for the bilinear form
*   KABN   : Number of additive terms in the bilinear form
*   ICUB   : Cubature formula identifier for CB2Q
*   ISYMM  : =0: don't respect symmetry, create matrix in structure 7
*            =1: matrix is symmetric; create matrix in structure 8
*            =2: matrix is symmetric and in structure 7 or 8;
*                automatically determine matrix structure
*   IPARAM : array [1..*] of integer
*            User defined parameter block, passed to COEFF
*   DPARAM : array [1..*] of double
*            User defined parameter block, passed to COEFF
* Out:
*   The entries of the matrix are added to DA.
*
* Out (COMMON block):
*   IER    = -116 for wrong value in array KAB
************************************************************************

      SUBROUTINE AB09X(DA,KCOL,KLD,NA,NBLOC,KOFF,KVERT,KMID,
     *                DCORVG,TRIA,ELE1,ELE2,ELE3,COEFF,
     *                BCON,COECON,KAB,KABN,ICUB,IPARAM,DPARAM)

      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'ccub.inc'
      
C     parameters
      
      INTEGER KCOL(*),KVERT(NNVE,*),KMID(NNVE,*),KOFF(*)
      INTEGER KLD(*),KAB(2,NNAB,*),KABN(*)
      INTEGER KENTRY(NNBAS,NNBAS),TRIA(SZTRIA)
      DOUBLE PRECISION DA(*),DCORVG(2,*),COECON(NNDER,NNDER,*)
      LOGICAL BCON(*)
      INTEGER NA,ICUB,NBLOC

      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

      DOUBLE PRECISION COEFF
      EXTERNAL ELE1,ELE2,ELE3,COEFF

C     externals

      INTEGER NDFL
      EXTERNAL NDFL

C     local variables

      LOGICAL BFIRST,BCON0
      INTEGER I,J,IBLOC,IA,IB,JCOL,NELE,IBAS
      INTEGER IDOFE,JDOFE,JCOL0,IDFG,IDER,ITYP1(3)
      INTEGER IVE,JP,IALBET,JCOLB,IABN,IAB
      DOUBLE PRECISION AUX,DJ1,DJ2,DJ3,DJ4,XI1,XI2,XX,YY,OM,DB
      DOUBLE PRECISION DBAS1(NNBAS,NNDER,3)
      LOGICAL BDER1(NNDER,3),BELE(3)
      INTEGER KDFG1(NNBAS,3),KDFL1(NNBAS,3),IDFL1(3)

      SUB='AB09'
      IF (ICHECK.GE.997) CALL OTRC('AB09  ','12/04/90')

C     Preparation - evaluation of parameters

      IER=0
      
C     Which derivatives of basis functions are needed?
C     Check the descriptors of the bilinear form and set BDER
C     according to these.
C     We initialise BDER1 here, which is later copied to BDER
C     depending on which element (ELE1, ELE2, ELE3) we call.

      DO IDER=1,NNDER
        DO J=1,2
          BDER1(IDER,J)=.FALSE.
        END DO
      END DO
      
      DO IBLOC=1,NBLOC
        DO IABN=1,KABN(IBLOC)
          DO J=1,2
          
C           The desriptor KAB gives directly the derivative
C           which is to be computed!

            IAB=KAB(J,IABN,IBLOC)

            IF (IAB.LE.0.OR.IAB.GT.NNDER) THEN
            
C             Oops, invalid descriptor
            
              WRITE (CPARAM,'(I15)') IBLOC
              CALL WERR(-116,'AB09  ')
              GOTO 99999
            END IF
            BDER1(IAB,J)=.TRUE.
          END DO
        END DO
      END DO
      
C     Dummy call of all ELEn to get the element identifiers

      DO I=1,3
      
        ITYP1(I)=-1
        
C       Crippled switch-case, depending on I:        
        
        GOTO (41,42,43) I
        
41        CALL ELE1(0D0,0D0,ITYP1(1))
          GOTO 44

42        CALL ELE2(0D0,0D0,ITYP1(2))
          GOTO 44

43        CALL ELE3(0D0,0D0,ITYP1(3))

44      CONTINUE
        IF (IER.NE.0) GOTO 99999

C       Get the number of local DOF's into IDFL1:
      
        IDFL1(I)=NDFL(ITYP1(I))
        
      END DO
      
C     Check which elements between ELE1, ELE2 and ELE3 are the same.
C     NELE=1 if ELE1=ELE2 (e.g. Q1/Q1), 
C     NELE=2 if ELE1<>ELE2 (e.g. Q1~/Q0)
      
      BELE(1)=.TRUE.
      BELE(2)=ITYP1(1).NE.ITYP1(2)
      IF (BELE(2)) THEN
        NELE=2
      ELSE
        NELE=1
      ENDIF
      
C     Set BELE(3)=true if ELE3 <> (ELE1,ELE2)
      
      BELE(3)=ITYP1(3).NE.ITYP1(1).AND.ITYP1(3).NE.ITYP1(2)
      
C     Initialise the cubature formula
      
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999

C     Dummy call of COEFF for nonlinear problems.
C     COEFF must set BDER(IDER)=.TRUE. in the /ELEM/ COMMON block
C     derivative IDER is needed from the element.
C     This allows to tell the element to calculate additional
C     derivarives to the ones defined by KAB.

      BFIRST=.TRUE.
      AUX=COEFF(0D0,0D0,-1,-1,0,BFIRST)

C     Check all matrix blocks if their bilinear form has constant
C     coefficients. In that case, get the coefficient and save it
C     in COECON.
C     BCON0 is set to TRUE if all blocks have constant coefficients.      

      BCON0=.TRUE.
      DO IBLOC=1,NBLOC
        IF (BCON(IBLOC)) THEN
          DO I=1,KABN(IBLOC)
            IA=KAB(1,I,IBLOC)
            IB=KAB(2,I,IBLOC)
            COECON(IA,IB,IBLOC)=COEFF(0D0,0D0,IA,IB,IBLOC,BFIRST,TRIA,
     *                                IPARAM,DPARAM)
          END DO
        ELSE
          BCON0=.FALSE.
        END IF
      END DO
      
C     Now let's calculate the matrix - storage technique 9.
C
C     Dummy call - ELE may save arithmetic operations.
C     ELE can precalculate the values in the cubature points if
C     it supports that...

      ICUBP=ICUB

      DO I=1,NELE
      
C       Initialise every element type only once:
      
        IF (BELE(I)) THEN
        
C         Copy the right BDER1 array to BDER, then call ELE:
        
          DO IDER=1,NNDER
            BDER(IDER)=BDER1(IDER,I)
          END DO
          
C         Crippled switch-case, using GOTO:
          
          GOTO (12,13,14) I
12          CALL ELE1(0D0,0D0,-2)
            GOTO 10
13          CALL ELE2(0D0,0D0,-2)
            GOTO 10
14          CALL ELE3(0D0,0D0,-2)
10        CONTINUE
      
        END IF
      
      END DO

C     To calculate the matrix, we loop over the elements.

      DO IEL=1,TRIA(ONEL)
      
C       The outstanding feature with finite elements is: A basis
C       function for a DOF on one element has common support only
C       with the DOF's on the same element! E.g. for Q1:
C
C              #. . .#. . .#. . .#
C              .     .     .     .
C              .  *  .  *  .  *  .
C              .     .     .     .
C              #-----#-----#. . .#
C              |     | O O |     .
C              |     | IEL |  *  .
C              |     | O O |     .
C              #-----X-----#. . .#
C              |     |     |     .
C              |     |     |  *  .
C              |     |     |     .
C              #-----#-----#. . .#
C              
C       --> On element IEL, the basis function at "X" of element ELE1
C           only interacts with the basis functions in "O" of element 
C           ELE2. Elements in the neighbourhood ("*") have no support, 
C           therefore we only have to collect all "O" DOF's.
C
C       Call NDFGLX to get the global DOF's of both ELE's
C       on our current element (the "X" and all "O"'s). 
C       We don't need the local DOF's, so by setting IPAR=0,
C       the call will only fill KDFG.

        DO I=1,3
          IF (BELE(I)) THEN
            CALL NDFGLX(TRIA,IEL,1,ITYP1(I),KVERT,KMID,
     *                  KDFG1(1,I),KDFL1(1,I))
          END IF
        END DO

C       For the assembly of the global matrix, we use a "local"
C       approach. At first we build a "local" system matrix according
C       to the current element. This contains all additive
C       contributions of element IEL, which are later added at the
C       right positions to the elements in the global system matrix.
C
C       We have IDFL DOF's per element, therefore there are
C       IDFL(1)*IDFL(2) tupels of basis-/testfunctions (phi_i,psi_j) 
C       "active" (i.e. have common support) on our current element, 
C       each giving an additive contribution to the system matrix.
C
C       We build a quadratic IDFL(1)*IDFL(2) local matrix:
C       KENTRY(1..IDFL(1),1..IDFL(2)) receives the position in the 
C         global system matrix, where the corresponding value has 
C         to be  added to.
C       (The corresponding contrbutions can be saved separately, 
C        but we directly add them to the global matrix in this 
C        approach.)
C
C       Loop through the local matrix to initialise it:

        DO JDOFE=1,IDFL1(1)
        
C         Row JDOFE of the local matrix corresponds 
C         to row=global DOF KDFG(JDOFE) in the global matrix.
C         This is the "X" in the above picture.
C         Get the starting position of the corresponding row
C         to JCOL0:

          JCOL0=KLD(KDFG1(JDOFE,1))
          
C         Now we loop through the other DOF's on the current element
C         (the "O"'s).
C         All these have common support with our current basis function
C         and will therefore give an additive value to the global
C         matrix.
          
          DO IDOFE=1,IDFL1(NELE)
            
C           Get the global DOF of the "O" which interacts with 
C           our "X".
            
            IDFG=KDFG1(IDOFE,NELE)
            
C           Starting in JCOL0 (which points to the beginning of
C           the line initially), loop through the elements in
C           the row to find the position of column IDFG.
C           Jump out of the DO loop if we find the column.
            
            DO JCOL=JCOL0,NA
              IF (KCOL(JCOL).EQ.IDFG) GOTO 111
            END DO
111         CONTINUE

C           Because columns in the global matrix are sorted 
C           ascendingly (except for the diagonal element),
C           the next search can start after the column we just found.
            
            JCOL0=JCOL+1
            
C           Save the position of the matrix entry into the local
C           matrix.
            
            KENTRY(JDOFE,IDOFE)=JCOL
            
          END DO ! IDOFE
          
        END DO ! JDOFE
        
C       To calculate the matrix contributions, we have to evaluate
C       ELEn on our current element. For this, we have to tell ELE
C       where we are:
C       - Store the vertex numbers of the current element to KVE
C       - Store the coordinates of the vertices of IEL to DX/DY in 
C         the ELE COMMON block
C       - Calculate the Jacobian of the mapping from the reference
C         element [-1,1]^2 to the real element into DJAC(1..2,1..2)
C         in the /ELEM/ COMMON block
C       - Calculate the Jacobian determinant of DJAC into DETJ in the
C         /ELEM/ COMMON block.
C       The last two points are necessary to correctly calculate
C       derivatives of the element!
C
C       So at first initialise KVE and DX/DY:

        DO IVE = 1, TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO
  
C       Now we want to calculate the Jacobian determinant.
C
C       We want to perform a transformation from the reference quadrilateral
C       [-1,1]x[-1,1] onto a "real" quadrilaterl:
C
C         (-1,1) ___ (1,1)         (x4,y4)   ___  (x3,y3)
C               |  |                        /  \
C               |__|           =>          /____\
C        (-1,-1)    (1,-1)         (x1,y1)        (x2,y2)   
C
C       By theory this can by done with a bilinear mapping, i.e. a mapping
C       Psi:R^2->R^2  of the form:
C
C            Psi (xi1) = ( a1 + a2*xi1 + a3*xi2 + a4*xi1*xi2 )
C                (xi2)   ( b1 + b2*xi1 + b3*xi2 + b4*xi1*xi2 )
C
C       Our special transformation has to map:
C
C        Psi(-1,-1) = (x1,y1)
C        Psi( 1,-1) = (x2,y2)
C        Psi( 1, 1) = (x3,y3)
C        Psi(-1, 1) = (x4,y4)
C
C       This gives the linear system:
C
C        a1 - a2 - a3 - a4 = x1       b1 - b2 - b3 - b4 = y1
C        a1 + a2 - a3 - a4 = x2       b1 + b2 - b3 - b4 = y2
C        a1 + a2 + a3 + a4 = x3       b1 + b2 + b3 + b4 = y3
C        a1 - a2 + a3 - a4 = x4       b1 - b2 + b3 - b4 = y4
C
C       Reorder this to calculate the ai:
C
C        a1 = 1/4 * ( x1 + x2 + x3 + x4)      b1 = 1/4 * ( y1 + y2 + y3 + y4)
C        a2 = 1/4 * (-x1 + x2 + x3 - x4)      b2 = 1/4 * (-y1 + y2 + y3 - y4)
C        a3 = 1/4 * (-x1 - x2 + x3 + x4)      b3 = 1/4 * (-y1 - y2 + y3 + y4)
C        a4 = 1/4 * ( x1 - x2 + x3 - x4)      b4 = 1/4 * ( y1 - y2 + y3 - y4)
C
C       The factors in the brackets in these equations are only dependent on 
C       the corners of the quadrilateral, not of the current point. So they
C       are constant for all points we want to map from the reference
C       element to the real one. We call them here "auxiliary Jacobian factors".
C       They are calculated in advance:
C
C         J1 = 1/2 * (-x1 - x2 + x3 - x4)
C         J2 = 1/2 * ( x1 - x2 + x3 - x4)
C
C         J3 = 1/2 * (-y1 + y2 - y3 + y4)
C         J4 = 1/2 * (-y1 + y2 + y3 - y4)

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
        
C       Using these, one can write:

C         a1  =  1/4 * ( x1 + x2 + x3 + x4)  =  1/2 * (x1 + x2 + J1)
C         a2  =  1/4 * (-x1 + x2 + x3 - x4)  =  1/2 * (x2 - x1 + J2)
C         a3  =  1/4 * (-x1 - x2 + x3 + x4)  =  1/2 * J1 
C         a4  =  1/4 * ( x1 - x2 + x3 - x4)  =  1/2 * J2
C
C         b1  =  1/4 * ( y1 + y2 + y3 + y4)  =  1/2 * (y1 + y3 + J3)
C         b2  =  1/4 * (-y1 + y2 + y3 - y4)  =  1/2 * J4
C         b3  =  1/4 * (-y1 - y2 + y3 + y4)  =  1/2 * (y3 - y1 - J4)
C         b4  =  1/4 * ( y1 - y2 + y3 - y4)  =  1/2 * J3
C
C       Now loop through all cubature points to build the local
C       matrix:

        DO ICUBP = 1, NCUBP

C         Get the coordinates of the cubature point on the reference
C         element from DXI in the /CUB/ COMMON block:

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)

C         Calculate the Jacobian of the mapping from the reference
C         element to the real element, based on:
C
C DPhi ( xi1 )  =  ( a2 + a4*xi2    a3 + a4*xi1 )
C      ( xi2 )     ( b2 + b4*xi2    b3 + b4*xi1 )
C
C               =  ( 1/2*(x2-x1+J2) + 1/2*J2*xi2              1/2*J1 + 1/2*J2*xi1 )
C                  (         1/2*J4 + 1/2*J3*xi2      1/2*(y3-y1-J4) + 1/2*J3*xi1 )

          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1

C         and from that its determinant:
C
C           det DPhi (xi1,xi2)  =  (DPhi[1,1]*DPhi[2,2] - DPhi[1,2]*DPhi[2,1]) (xi1,xi2)

          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)

C         This allows to calculate the current weighting factor in 
C         the cubature formula:

          OM=DOMEGA(ICUBP)*DETJ

C         Call the element subroutines to calculate the values
C         of the basis functions in the cubature point (XI1,XI2).
C
C         The variable ICUBP which we use in the loop tells the
C         element the number of the cubature point where we are
C         at the moment. Therefore, when we call ELE with IPAR=-1,
C         ELE can decide on its own whether to calculate the
C         values on the basis functions again or to use the
C         precalculated ones.

          DO I=1,3
          
C           Evaluate every element only once:
          
            IF (BELE(I)) THEN
            
C             Tell the element what to calculate:
            
              DO IDER=1,NNDER
                BDER(IDER)=BDER1(IDER,I)
              END DO
              
C             Crippled SWITCH-CASE using GOTO.
C             Evaluate the element:
              
              GOTO (22,23,24) I
22              CALL ELE1(XI1,XI2,-3)
                GOTO 25
23              CALL ELE2(XI1,XI2,-3)
                GOTO 25
24              CALL ELE3(XI1,XI2,-3)
25            CONTINUE

C             Copy the calculated values of the basis functions
C             to a local copy DBAS1 - we later use them.

              DO IDER=1,NNDER
                IF (BDER1(IDER,I)) THEN
                  DO IBAS=1,IDFL1(I)
                    DBAS1(IBAS,IDER,I)=DBAS(IBAS,IDER)
                  END DO
                END IF
              END DO
              
            END IF
            
          END DO ! I

C         In case we don't have constant coefficients,
C         also calculate the real coordinates of the cubature point: 
C           (XX,YY) := Psi(XI1,XI2) 

          IF (.NOT.BCON0) THEN
            XX = 0.5D0*(DX(1)+DX(2)+DJ1)
     *          +0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *          +0.5D0*DJ1*XI2
     *          +0.5D0*DJ2*XI1*XI2
     
            YY = 0.5D0*(DY(1)+DY(3)+DJ3)
     *          +0.5D0*DJ4*XI1
     *          +0.5D0*(DY(3)-DY(1)-DJ4)*XI2
     *          -0.5D0*DJ3*XI1*XI2
          END IF

C         Sum up over all pairs of multiindices.
C         Loop though all matrix blocks:

          BFIRST=.TRUE.
          DO IBLOC=1,NBLOC
          
C           Loop through all KABN prescribed combinations of the
C           derivatives of the basis functions.

            DO IALBET=1,KABN(IBLOC)
            
C             Get from KAB the type of the derivatives for the 
C             test and trial functions. The summand we calculate
C             here will be:
C
C             int_... ( phi_i )_IA  *  ( psi_j )_IB
C
C             -> Ix=0: function value, 
C                  =1: first derivative, 
C                  =2: 2nd derivative,...
              
              IA=KAB(1,IALBET,IBLOC)
              IB=KAB(2,IALBET,IBLOC)
              
C             Get the coefficient of the basis function - either
C             from COEFF (for non-constant coefficients) or take
C             the precalculated one from COECON.
C             Multiply it with the cubature weight OM.
              
              IF (.NOT.BCON(IBLOC)) THEN
               AUX=COEFF(XX,YY,IA,IB,IBLOC,BFIRST,TRIA,IPARAM,DPARAM)*OM
              ELSE
               AUX=COECON(IA,IB,IBLOC)*OM
              ENDIF
              
C             Now loop through all possible combinations of DOF's
C             in the current cubature point. The outer loop
C             loops through the "X" in the above picture:

              DO JDOFE=1,IDFL1(1)
              
C               Get the value of the basis function 
C               phi_i (our "X") in the cubature point:
              
                DB=DBAS1(KDFL1(JDOFE,1),IA,1)
                
C               Perform an inner loop through the other cubature
C               points (the "O"'s). 

                DO IDOFE=1,IDFL1(NELE)

C                 Get the value of the basis function 
C                 psi_j (our "O") in the cubature point. This is
C                 DBAS(KDFL(IDOFE),IB).
C                 Them multiply:
C                    DB * DBAS(..) * AUX
C                 ~= phi_i * psi_j * coefficient * cub.weight
C                 Summing this up gives the integral, so the contribution
C                 to the global matrix. 
C
C                 Simply summing up DB * DBAS(..) * AUX would give
C                 the coefficient of the local matrix. As everything
C                 is only summed up here, we directly incorporate
C                 this value into the global matrix at the position 
C                 given by KENTRY.

                  JCOLB=KENTRY(JDOFE,IDOFE)+KOFF(IBLOC)
                  DA(JCOLB)=DA(JCOLB)+
     *                      DB*DBAS1(KDFL1(IDOFE,NELE),IB,NELE)*AUX
     
C                 The nonlinearity is not yet implemented!
     
                END DO ! IDOFE
                
              END DO ! JDOFE
              
            END DO ! IALBET

            BFIRST=.FALSE.
            
          END DO ! IBLOC
        END DO ! ICUBP
      END DO ! IEL
 
C     Assembly finished.
 
99999 END

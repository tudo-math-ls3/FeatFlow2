***********************************************************************
* Recovered gradient technique
*
* This file implements recovered gradient routines in extended calling
* convention, supporting arbitrary grids without COMMON blocks.
*
* For the time being these routines only support ZZ-recovery
* technique on arbitrary quadrilateral grids.
***********************************************************************

**************************************************************************
* Local gradient with transformation to reference element
*
* Computes the local recovered gradient in vertex IVT by a linear 
* approach. For this purpose the values in the center of all elements 
* around the inner vertex JVT are taken.
*
* JVT must be an inner node with at least 3 elements surrounding.
* IVT can be a boundary node. For inner nodes JVT=IVT is usually
* the best choice.
*
* Approach:
*  -Compute the values of the derivatives of the basis functions
*   in the center of all elements surrounding the node JVT
*  -Create a bilinear polynomial with corners in the midpoints of
*   these elements. Prescribe the values of this polynomial in
*   these points.
*  -Take the value of the polynomial in the center (or better to say, 
*   in IVT) and use this as a recovered value.
*
* In:
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure describing the mesh
*  DU     - array [1..*] of double
*           Solution vector
*  DCORVG,
*  KVERT,
*  KMID,
*  KVEL   - Usual geometry information. Must correspond to the handles
*           in the TRIA-structure. Transferred via parameters for
*           easier and quicker access.
*  
*  X,
*  Y      - Coordinates of the point where the reconstructed gradient
*           should be computed
*  IVT    - Vertex number where to evaluate. Should correspond to (X,Y)
*           if possible.
*  JVT    - "Inner reference" node, 1..NVT. Must be a vertex in the inner 
*           of the domain, completely surrounded by elements. JVT should
*           be as close to (X,Y) as possible.
*
*  NVEL   - Maximum number of elements meeting in a vertex; must 
*           correspond to NVEL in the TRIA-structure.
*  NVELJ  - Number of elements meeting in current vertex JVT
*  ELE    - Element that corresponds to solution vector DU; this is
*           also used for reconstruction. Must be conformal, parametric
*           element (Q1,Q2,...).
*  WORK   - array [1..*] of double
*           Workspace vector; will be ignored if LSIZE=-1
*  LWORK  - size of workspace vector, or
*           -1=determine LSIZE, return size of workspace vector in LSIZE
*
* Out:
*  DGRADX,
*  DGRADY - reconstructed gradient in vertex IVT
*
* The routine can be called with a Dummy-call using LSIZE=-1 to determine
* the size of the workspace vector. In this case no computation will be
* done. Instead the routine will return in LSIZE the necessary size of
* the workspace vector WORK. To determine the maximum amount of memory
* necessary for all nodes, the caller should use NVELJ=NNVEL and NVE=NNVE.
*
* The usage of IVT, JVT and (X,Y) is the following:
*
* a) If (X,Y) are the coordinates of a vertex surrounded by at least,
*    3 elementts, IVT=JVT=(X,Y) must be used.
* b) If (X,Y) are the coordinates of a vertex with less than adjacent
*    elements (typically boundary vertex), IVT=(X,Y) must 
*    be used, JVT must be an inner vertex surrounded by elements.
* c) If (X,Y) are the coordinates of an inner node, not corresponding
*    to a vertex (element midpoint e.g.), 
*    - use IVT=JVT=vertex number of an inner vertex surrounded by elements
*    - use JVT=vertex number of an inner vertex surrounded by elements,
*      IVT=vertex number of a neighbour of JVT
* If IVT is not a neighbour of JVT, the routine will find a neighbour
* automatically.
*           
*             *------*------*  *------*------* *------*------*
*             |      |      |  |      |      | |      |      |
*             |      |      |  |      |      | |      |      |
*             |     IVT     |  |      |      | |     IVT     |
*             *-----JVT-----*  *-----JVT-----* *-----JVT-----*
*             |    (X,Y)    |  |      |      | |      |      |
*             |      |      |  |      |      | |      | (X,Y)|
*             |      |      |  |      |      | |      |      |
*             *------*------*  *-----IVT-----* *-----IVT-----*
*                                   (X,Y)

************************************************************************** 

      SUBROUTINE GRDLO2(TRIA,DU,DGRADX,DGRADY,
     *                  DCORVG,KVERT,KMID,KVEL,
     *                  X,Y,IVT,JVT,NVELJ,NVEL,ELE,
     *                  WORK, LWORK)

      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'stria.inc'

C parameters

      INTEGER LWORK,TRIA(SZTRIA)
      DOUBLE PRECISION WORK(*)
      DOUBLE PRECISION DU(*),DGRADX,DGRADY

      DOUBLE PRECISION DCORVG(2,*),X,Y
      INTEGER IVT,JVT,NVELJ,NVEL
      INTEGER KVERT(NNVE,*),KVEL(NVEL,*),KMID(NNVE,*)
      EXTERNAL ELE

C externals      
      
      INTEGER NDFL
      EXTERNAL NDFL,E011

C local variables

      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL,MTSIZE
      DOUBLE PRECISION DCOORD (2,4),DJF(2,2)
      DOUBLE PRECISION XDER,YDER
      INTEGER IVEL,NROWS,NCOLS,LNEED,ROW,NMXDEG,NMNMS
      DOUBLE PRECISION HSCALE,XCENTR,YCENTR
      DOUBLE PRECISION XX,YY,XHAT,YHAT,XI1,XI2
      INTEGER IELTYP,ICUB,INEIG,IVE,NVE,JP,IBAS,IEQ,ILO,I,J,K
      DOUBLE PRECISION DALPHA, DBETA, DDSTI, DXNEIG, DYNEIG
      
C Description of the method:
C
C We want to create a local polynom that interpolates the values of 
C each partial derivative in the midpoints if the four elements
C surrounding JVT:
C
C   |-------|-------|
C   |       |       |
C   |   X.......X   |
C   |   .   |   .   |
C   |---.--JVT--.--IVT
C   |   .   |   .   |
C   |   X.......X   |
C   |       |       |
C   |-------|-------|
C
C This polynom is then evaluated in the point IVT (where for inner 
C vertices IVT=JVT is used), which gives an approximation to the 
C (partial) derivative in that point.
C
C In a more general case (not Q1, but Q2, Q3,...), not the midpoints of
C the elements are used, but rather the Gauss-points for that element.
C For every element we have to use "exactly the correct" Gauss-points,
C otherwise no superconvergence is to expect - i.e. Q1=midpoints,
C Q2=2x2 Gauss,...
C
C So at first we check our element about its type. We only support
C conforming elements here, so in case of another element we can
C immediately stop with an error. Otherwise we have to choose the
C correct cubature formula:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)

      IF (IELTYP.EQ.11) THEN
        ICUB = 1
        NMXDEG = 1
        NCOLS = 4
      ELSE IF (IELTYP.EQ.13) THEN
        ICUB = 4
        NMXDEG = 2
        NCOLS = 9
      ELSE IF (IELTYP.EQ.14) THEN
        ICUB = 8
        NMXDEG = 3
        NCOLS = 16
      ELSE
        WRITE (*,*) 'GRDLO2: Illegal element type: ',IELTYP
      END IF
 
C Number of monoms per variable is 1 larger that maximum degree:

      NMNMS = NMXDEG+1
 
C Initialize the cubature formula for quadrilateral elements 

      CALL CB2Q(ICUB)
      
C Get some parameters from the triangulation:

      NVE = TRIA(ONVE)
      
C Small description of the method:
C
C What we use here is a simple polynomial interpolation. In the
C Gaussian points of the elements surrounding JVT we compute the
C values of the derivatives. Then we create an interpolating polynomial
C that interpolates these values. We evaluate the polynomial in
C IVT to create an approximation to the derivatives of our function
C there.
C
C The polynomial we want to create has the typical form:
C
C    p(x,y) = a + bx + cy + dxy + ...
C
C such that:
C
C    p(Gauss-point) = x-derivative (Gauss-point)
C
C (the same for the Y-derivative) with the number of monoms depending
C on the element type: For Q1 we take a bilinear polynom (4 monoms),
C for Q2 a biquadratic polynom (9 monoms) and so on.
C
C This leads to two linear systems of the form
C
C   (1  X1  Y1  X1Y1  ...)  ( a )   ( DUX (X1,Y1) )
C   (1  X2  Y2  X2Y2  ...)  ( b ) = ( DUX (X2,Y2) )
C   (1  X3  Y3  X3Y3  ...)  ( c )   ( DUX (X3,Y3) )
C   ...
C
C to determine the coefficients a,b,c,... of the polynomial for
C the X-derivative and for the Y-derivative.
C
C The matrix in this system is generally a rectangle matrix.
C For 4 Q1-elements, it's 4x4, for 4 Q2-elements it's 16x9,...
C We overcome the problem of solving this by a Least square
C approach, see later.
C
C The crucial point at the beginning is the memory demand in WORK.
C We organize WORK to take two matrices and a couple of
C vectors. The NCOLS variable set above gives the number of
C columns in the matrix, which is exactly the number of monoms
C in the polynomial. The number of rows are calculated
C by "number of cubature points * number of elements around JVT".
C NCUBP is taken from the /CUB/ Common-block.

      NROWS = NCUBP*MAX(4,NVELJ)
      
C Where we make sure that we have at least 4 rows
C (also in the case that only 3 elements meet in a vertex! There's
C a trick to handle that case, see later).
C
C We have to store at least 2 matrices + 2 vectors for the X- and
C Y-derivatives in all the points. The size of the matrix is

      MTSIZE = NROWS*NCOLS

C And we need totally a memory amound of:

      LNEED = 2*MTSIZE + 2*NROWS
      
C Does the caller want to know the memory?

      IF (LWORK.LT.0) THEN
        LWORK = LNEED
        RETURN
      END IF
      
C Is there enough memory?

      IF (LWORK.LT.LNEED) THEN
C       Should not happen on correct implementation. Therefore we make
C       a direct message to the programmer that something went 
C       wrong here
        WRITE (*,*) 'GRDLO2: Not enough memory!'
        RETURN
      END IF
      
C Calculate local number of degrees of freedom of our element type
      
      IDFL = NDFL(IELTYP)
      
C Ok, our central node is JVT. We will now make a loop over all
C elements surrounding that vertex. On each element we will calculate
C in every Gauss-point the values of the X- and Y-derivative.
C Furthermore we will save all the coordinates for later computation.
C These values will be saved in the WORK-array for later use.
C
C So make a loop over the elements surrounding JVT: 

      DO IVEL = 1, NVELJ

C       On which element are we?
 
        IEL = KVEL(IVEL, JVT)

C       Calculate for that element the global degrees of freedom 
C       into KDFG:

        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL) 

C       Transfer the coordinates of the corners into the COMMON-Block
C       (for call to ELE) as well as to DCOORD for the calculation of
C       the transformation:

        DO IVE=1,NVE           
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP) 
          DCOORD(1,IVE) = DCORVG(1,JP)
          DCOORD(2,IVE) = DCORVG(2,JP)
        END DO
         
C       Calculate auxiliary Jacobian factors of the transformation

        CALL QINIJF (DCOORD,DJF)

C       Loop over the cubature points on our element to calculate the
C       coefficients of the FE-basis functions in that point...

        DO ICUBP = 1, NCUBP

C         Get cubature points on reference element from Common block

          XI1 = DXI(ICUBP,1)    
          XI2 = DXI(ICUBP,2)
            
C         Calculate Jacobian matrix of transformation and its
C         determinant, as well as the real coordinates of our point

          CALL QTRAF (DCOORD,DJF,DJAC,DETJ,XI1,XI2,XX,YY)

C         Evaluate x-and y-derivative at gauss point

          CALL ELE(XI1,XI2,0) 

C         Loop over the basis functions. Sum up the node values
C         weighted by the coefficients of the basis functions to
C         calculate the derivative in the gauss point.

          XDER = 0D0
          YDER = 0D0
            
          DO IBAS = 1,IDFL 
            IEQ=KDFG(IBAS)
            ILO=KDFL(IBAS)
               
C           Sum up X-Derivative

            XDER = XDER + DBAS(ILO,2)*DU(IEQ)

C           Sum up Y-Derivative

            YDER = YDER + DBAS(ILO,3)*DU(IEQ)
          END DO

C         OK, X-derivative is in XDER, Y-derivative in YDER.
C         Save both values as well as the coordinates in our
C         matrix: The X-coordinate in the 2nd column, the Y-coordinate
C         in the 2nd column:

          ROW = (IVEL-1)*NCUBP+ICUBP

          WORK (NROWS + ROW)           = XX
          WORK (2*NROWS + ROW)         = YY
          WORK (2*MTSIZE + ROW)         = XDER
          WORK (2*MTSIZE + NROWS + ROW) = YDER

        END DO
                 
      END DO
      
C The first column of our matrix is typically filled with 1,
C as it represents x^0*y^0:

      DO I=1,NROWS
        WORK(I) = 1D0
      END DO   

C Now we know (in GPTS) all coordinates of all cubature points as well
C as the X- and Y-derivatives there.
C We now want to calculate a matrix A with each line of the form
C    ( 1  X  Y  X*Y ... )
C to calculate the coefficients of our local polynomial, which
C interpolates in the (former cubature-) points the values of the
C derivatives we just collected.
C
C But take care here! Don't use the X/Y-coordinates and try to solve! 
C This can easily lead to numerical instabilities, as it would only 
C consider the case that we are using a standard coordinate system. 
C One example is the Q1-element:
C
C We can of course not directly construct a bilinear polynom of the
C form "a+b*x+c*y+d*x*y" that takes the values in the four midpoints,
C because if the four elements are turned by 45 degrees, the four
C midpoints are forming a "cross":
C 
C     / \
C   /  X  \
C  < X-|-X >
C   \  X  /
C     \ /
C
C and the polynom does not exist in that case (not unisolvent!).
C Therefore we introduce a local coordinate system with midpoint in
C JVT. For this we take one neigbour of JVT and take the
C connecting edge as X-axes. The orthogonal vector of this then
C forms the Y-axes.
C
C For this local coordinate system we will use JVT as center.
C To form it, we need to determine a neighbour of JVT:
C
C Try to use IVT as neighbour if possible to avoid numerical
C instabilities a little bit. If this is not possible (e.g. because
C IVT is a diagonal point like on if the current "corner"-elements
C of a quadrilateral), or not necessary (because we are in the "inner"
C of the domain, so IVT=JVT), search for another arbitrary neighbour.
      
      IF (IVT.NE.JVT) THEN
      
        DO I=1,IVEL
        
          IEL = KVEL(I,JVT)
          
          DO J=1,NNVE
          
            IF (KVERT(J,IEL).EQ.JVT) THEN
            
              INEIG = KVERT(MOD (J,NNVE)+1,IEL)
              IF (INEIG.EQ.IVT) GOTO 100
              INEIG = KVERT(MOD (J+NNVE-1,NNVE)+1,IEL)
              IF (INEIG.EQ.IVT) GOTO 100
              
            END IF
        
          END DO
        
        END DO
        
      END IF
      
C Take an arbitrary element of JVT and look there for a neighbour

      IEL = KVEL(1,JVT)

      DO I=1,NNVE
        IF (KVERT(I,IEL).EQ.JVT) THEN
          INEIG = KVERT(MOD (I,NNVE)+1,IEL)
          GOTO 100
        END IF
      END DO
      WRITE (*,*) 'GRDLOK: KVERT destroyed!'
      STOP
      
100   CONTINUE        

C Now we have a neighbour INEIG of JVT as well as an element IEL
C containing JVT and INEIG.
C
C Determine the coordinates of our central/reference point 
C in the patch. If IVT=JVT we will later reconstruct the gradient there
C directly (or in (X,Y), respectively). 
C In case of a boundary node, then IVT <> JVT and JVT is 
C the center of our local coordinate system.

      XCENTR = DCORVG(1,JVT)
      YCENTR = DCORVG(2,JVT)

C Determine the coordinates of the neighbour and its coordinates:

      DXNEIG = DCORVG(1,INEIG)
      DYNEIG = DCORVG(2,INEIG)

C Determine local X-axes: vector (DALPHA,DBETA) that points from our
C modpoint to our neighbour:
C            CENTER + (DALPHA,DBETA) = NEIGHBOUR
C DDSTI is 1/length of the vector.
C The transposed vector (-DBETA,DALPHA) will be our Y-axes.

      DDSTI = 1D0/DSQRT((DXNEIG - XCENTR)**2 + (DYNEIG - YCENTR)**2)
      DALPHA = (DXNEIG - XCENTR)*DDSTI
      DBETA =  (DYNEIG - YCENTR)*DDSTI

C We need to scale out coordinate system to maintain numerical
C stability. Determine a local h from the element size of the current
C element IEL:

      HSCALE=1D0/
     *         DSQRT((DCORVG(1,KVERT(1,IEL))-DCORVG(1,KVERT(3,IEL)))**2+
     *               (DCORVG(2,KVERT(1,IEL))-DCORVG(2,KVERT(3,IEL)))**2)

C Next we have to build a matrix for the coefficients of the
C polynomial. A standard way to form such a matrix is the structure:
C   (1  X  Y  XY  X^2Y  XY^2  ...)
C   (1  X  Y  XY  X^2Y  XY^2  ...)
C   (1  X  Y  XY  X^2Y  XY^2  ...)
C   ...
C Here X,Y are not the real coordinates of the cubature points,
C but the local coordinates of the cubature points in our local
C coordinate system!
C
C For the time being we saved the real X- and Y-coordinates there.
C So now we have to make a loop over the lines, translate
C the coordinates into the local coordinate system and compute the
C other monoms.
C
C From here on we define the matrix for the coefficients to have a 
C slightly different format than mentioned earlier. We will build
C the following matrix structure:
C
C   (1  X  X^2  X^3 .. Y  XY  X^2Y ...)  ( a )   ( DUX (X1,Y1) )
C   (1  X  X^2  X^3 .. Y  XY  X^2Y ...)  ( b ) = ( DUX (X2,Y2) )
C   (1  X  X^2  X^3 .. Y  XY  X^2Y ...)  ( c )   ( DUX (X3,Y3) )
C
C So the Y-coordinate will be saved not in the 2nd column, but in
C the (degree+1)'th column. This is to have all coefficients
C without the Y at the front, so the polynomials with Y=0 can
C be easier calculated later on.
C
C Make a loop over all "real" rows, ignoring the last line in case 
C that we have only 3 quadrs meeting in JVT (i.e. use NCUBP*NVELJ 
C instead of NROWS as loop count):

      DO ROW = 1,NCUBP*NVELJ

C       Current coordinates:

        XX = WORK (NROWS + ROW)
        YY = WORK (2*NROWS + ROW)
        
C       Determine the local coordinates of the cubature point in our
C       local coordinate system.
 
        XHAT = HSCALE*(DALPHA *(XX - XCENTR) + DBETA* (YY - YCENTR))
        YHAT = HSCALE*(-DBETA *(XX - XCENTR) + DALPHA*(YY - YCENTR))

C       And save these instead of XX/YY into the matrix. The X-
C       coordinate comes again into the 2nd column, the Y-coordinate
C       in the first column of the second monom group:

        WORK (NROWS + ROW)       = XHAT
        WORK (NMNMS*NROWS + ROW) = YHAT

C       Now the next big question: How many additional monoms does our
C       polynom have? Ok, NCOLS-3, but that's not enough to compute
C       them. We have set the maximum degree in the variable NMXDEG:
C       Q1->degree 1 (xy), Q2->degree 2 (x^2y^2),...
C       So now make a loop to calculate the missing monoms. 
C       Calculate in the following order:
C       x^2,x^3,x^4,...,X^NMXDEG
        
        DO J=2,NMXDEG
          WORK(J*NROWS + ROW) = XHAT**J
        END DO
        
C       xy,x^2y,x^3y,...

        DO J=1,NMXDEG
          WORK((NMNMS+J)*NROWS + ROW) = XHAT**J * YHAT
        END DO
        
C       and all the rest: y^2,xy^2,x^2y^2,...
        
        DO I=2,NMXDEG
          DO J=0,NMXDEG
            WORK((I*NMNMS+J)*NROWS + ROW) = XHAT**J * YHAT**I
          END DO
        END DO
        
      END DO
       
C Now a little trick for the case, that we only have 3 Q1-elements
C surrounding our node  JVT. Only in this case we have less
C rows than columns - but we don't have less degrees of freedom
C than necessary!
C
C We have a 4x4-system in WORK for the coefficients of the
C local polynomial p(x,y) = a+b*x+c*y+d*x*y. We want to invert
C this to get out coefficients a,b,c,d. Unfortunately if only
C 3 elements meet in the vertex JVT, our system is indefinite,
C and so we are not able to invert the matrix!
C
C The trick is now: In this case we change our local polynomial.
C We seek for a poplynom in the form p(x,y) = a+b*x+c*y, so
C a linear polynomial on a triangle element!
C This differs from the original polynomial only in that way, that
C the last factor is missing.
C
C To obtain a,b,c in this case, we replace the last column/row of
C the matrix by a unary vector. Then the matrix decomposes into a 
C 3x3-system for a,b,c and a 1x1-system for the "dummy"-coefficient d.
C    xxx0  a   midp1
C    xxx0  b = midp2
C    xxx0  c   midp3
C    0001  d   0
C We can invert that matrix and multiply it with the RHS to obtain
C a,b,c (and d). Furthermore if we set the last component of our RHS-
C vectors to 0 we ensure, that d=0. This way we can use the formula
C p(x,y) = a + b*x + c*y + d*x*y for all cases regardless of if less 
C than 4 elements meet in JVT! 

      IF ((NVELJ.EQ.3).AND.(IELTYP.EQ.11)) THEN
C Modify the matrix
        DO I=1,4
          WORK(3*NROWS + I) = 0D0
          WORK(I*NROWS + 4) = 0D0
        END DO
        WORK(3*NROWS + 4) = 1D0
        
C and the RHS-vectors = X/Y-derivative, which are saved
C behind the two matrix blocks:

        WORK (2*MTSIZE + 4)         = 0D0
        WORK (2*MTSIZE + NROWS + 4) = 0D0
      END IF

C Ok, the matrix is prepared, also the both RHS-vectors are prepared
C behind the two matrix blocks. What is the second matrix block
C good for?
C The problem is the case if more than 4 elements meet in JVT. In this
C case the matrix is not rectangular, i.e. we cannot simply solve
C it to get the coefficients of the polynomial.
C
C A way out of this is the Least-Squares approach:
C Instead of solving a system Ax=b, we solve the system A^TAx=A^Tb.
C This gives a "best approximating" solution X.
C The second matrix block serves as temporary array to compute A^TA.

      IF (NROWS.NE.NCOLS) THEN
      
C       Multiply both RHS by A^T. The last column of the second
C       matrix structure is used as temporary array.

        CALL LCP1 (WORK (2*MTSIZE+1),WORK (2*MTSIZE-NROWS+1),NROWS)
        CALL DGEMV ( 'T', NROWS, NCOLS, 1D0, WORK(1), NROWS, 
     *               WORK (2*MTSIZE-NROWS+1), 1,
     *               0D0, WORK (2*MTSIZE+1), 1 )

        CALL LCP1 (WORK (2*MTSIZE+NROWS+1),WORK (2*MTSIZE-NROWS+1),
     *             NROWS)
        CALL DGEMV ( 'T', NROWS, NCOLS, 1D0, WORK(1), NROWS, 
     *               WORK (2*MTSIZE-NROWS+1), 1,
     *               0D0, WORK (2*MTSIZE+NROWS+1), 1 )
     
C       Build A^TA: Copy the matrix to the 2nd matrix block,
C       multiply by A^T and save back to the first:

        CALL LCP1 (WORK(1),WORK(MTSIZE+1),MTSIZE)
        CALL DGEMM ( 'T', 'N', NCOLS, NCOLS, NROWS, 1D0, 
     *               WORK(MTSIZE+1), NROWS, WORK(MTSIZE+1), NROWS,
     *               0D0, WORK(1), NCOLS )

      END IF

C So we now have a matrix at WORK(1) and two RHS-vectors after both
C matrix blocks. 
C The next part is to solve two linear systems for the
C both RHS's to determine the coefficients of the polynomial.
C
C This we will do by LAPACK-routines. Use DGETRF at first to factorize
C our matrix in-place. DGETRF will write integers into our WORK-array,
C which is a double-precision array. But since integers are generally
C smaller in size, this will destroy nothing in WORK - there's enough
C space.

      CALL DGETRF ( NCOLS, NCOLS, WORK(1), NCOLS, WORK(MTSIZE+1), I )
     
      IF (I.NE.0) THEN
        PRINT *,'GRDLO2: Illegal matrix'
        STOP
      END IF
     
C And solve it two times:

      CALL DGETRS( 'N', NCOLS, 1, WORK(1), NCOLS, WORK(MTSIZE+1), 
     *             WORK(2*MTSIZE+1), NCOLS, I )

      CALL DGETRS( 'N', NCOLS, 1, WORK(1), NCOLS, WORK(MTSIZE+1), 
     *             WORK(2*MTSIZE+NROWS+1), NCOLS, I )

C That's it, we now have the two coefficient vectors in WORK(MTSIZE+1)
C and WORK(MTSIZE+NROWS+1).
C
C Use the calculated coefficients to compute the value of the
C polynomial in the point (X,Y). For this purpose calculate the
C local coordinates of the midpoint - which should be ~ 0, if
C IVT=JVT=(X,Y) is the center of a quadrilateral patch.

      XHAT = HSCALE*(DALPHA *(X - XCENTR) + 
     *                 DBETA* (Y - YCENTR))
      YHAT = HSCALE*(-DBETA *(X - XCENTR) + 
     *                 DALPHA*(Y - YCENTR))
             
C Remark that if IVT=JVT=(X,Y)=midpoint of the patch, XHAT=YHAT=0 and 
C so the value is determined by DRECVX(1) / DRECVY(1)!
C      
C We perform a loop over all monoms and add them together
C to the reconstructed gradient.
 
      DGRADX = 0D0
      DGRADY = 0D0
      
      K=1
      DO I=0,NMXDEG
        DO J=0,NMXDEG
          XX = WORK(2*MTSIZE+K)
          YY = WORK(2*MTSIZE+NROWS+K)
          
          DGRADX = DGRADX + XX * XHAT**DBLE(J) * YHAT**DBLE(I)
          DGRADY = DGRADY + YY * XHAT**DBLE(J) * YHAT**DBLE(I)
          
          K=K+1
        END DO
      END DO

      END

**************************************************************************
* Computation of reconstructed gradient with ZZ-technique,
* Extended calling convenction
*
* This routine calculates for a given solution vector DU based on an
* element ELE the reconstructed gradient DGRADX,DGRADY in each node
* of the grid.
*
* In:
*   TRIA    - array [1..SZTRIA] of integer
*             Triangulation structure of the mesh
*   DU      - array [1..NEQ] of double
*             Solution vector
*   NEQ     - Number of equations in DU 
*   ELE     - Element that corresponds to solution vector DU.
*             Must be conformal (Q1,Q2,...) parametric element.
* 
* Out:
*   DGRADX  - array [1..NEQ] of integer
*             X-Gradient of solution vector
*   DGRADY  - array [1..NEQ] of integer
*             Y-Gradient of solution vector
**************************************************************************

      SUBROUTINE XGRDR2(TRIA,DU,DGRADX,DGRADY,NEQ,ELE)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'stria.inc'
      
C parameters

      INTEGER NEQ,TRIA(SZTRIA)
      DOUBLE PRECISION DU(NEQ),DGRADX(NEQ),DGRADY(NEQ)
      EXTERNAL ELE

      INTEGER NDFL
      EXTERNAL NDFL

C local variables

      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL,IELTYP
      INTEGER IVEL,I,J
      INTEGER IEL1,IEL2
      
      INTEGER KVERT,KMID,KXNPR,KVEL,KCOUNT,KMEL,KCORVG,LWORK,LNEED
      INTEGER IEL,IVT,JVT,NEL,NVE,NVT,IDG
      
      INTEGER LCOUNT,NVEL,KV1,KV2,KNP,IJ(1)
      DOUBLE PRECISION D,X,Y,DGX,DGY,DJ(1)
     
C We need an extended triangulation structure for this routine:

      IF (TRIA(OTRIFL).EQ.0) THEN
        WRITE (*,*) 'XGRDR2: No extended triangulation information!'
        RETURN
      END IF
     
C Clear the destination arrays at first:

      CALL LCL1 (DGRADX,NEQ)
      CALL LCL1 (DGRADY,NEQ)
      
C How do we have to perform the recovering of the gradient now?
C Basically we would have to perform a loop about all vertices
C 1..NVT and call GRDLO2 for them. The result is then written into
C DGRADX/DGRADY.
C
C This works fine if Q1 is used. Unfortunately higher order 
C elements like Q2,Q3,... do not only use vertices, but also edge
C and element midpoints!
C
C What to do there? Well, taking the mean! In case of Q2 and an
C edge midpoint, there are two elements adjacent to it. Calculate
C the recovered gradient on both elements and take the mean:
C 
C *------*------*
C |      |      |
C |   > IMT <   |
C |      |      |
C *------*------*
C      
C An element midpoint is surrounded by 4 corner vertices. In this
C case there are 4 patches we have take into considerarion when taking
C the mean:
C
C *------*------*
C |      |      |
C |  X=======X  |
C |  I       I  |
C *--I   M   I--*
C |  I       I  |
C |  X=======X  |
C |      |      |
C *------*------*
C
C At first we have to reserve some memory for easier computation.
C For calculating the mean we need an array KCOUNT. This saves for
C every DOF how many times the recovered gradient has been calculated
C for that DOF. Afterwards it's used to take the mean.

      CALL ZNEW (NEQ,3,LCOUNT,'KCOUNT')

C Ask the recovery routine about how much memory is necessary:

      LNEED = -1
      CALL GRDLO2(TRIA,DJ,DJ(1),DJ(1),DJ,IJ,IJ,IJ,
     *            0D0,0D0,0,0,NNVE,TRIA(ONVEL),
     *            ELE, DJ, LNEED)

      CALL ZNEW (LNEED,-1,LWORK,'WORK  ')
      
C Transfer some adresses to local variables for easier access

      KCORVG = L(TRIA(OLCORVG))
      KVERT = L(TRIA(OLVERT))
      KMID  = L(TRIA(OLMID))
      KXNPR = L(TRIA(OLXNPR))
      KVEL  = L(TRIA(OLVEL))
      KMEL  = L(TRIA(OLMEL))
      NVEL = TRIA(ONVEL)
      NVE  = TRIA(ONVE)
      NVT = TRIA(ONVT)
      NEL = TRIA(ONEL)
      
      KCOUNT = L(LCOUNT)
      
C Ask the element about its type:

      IELTYP = -1
      CALL ELE (0D0,0D0,IELTYP)

C Calculate local number of degrees of freedom of our element type
      
      IDFL = NDFL(IELTYP)
      
C Now perform a loop over all elements:

      DO IEL = 1,NEL
      
C       On the current element calculate the global DOF's into KDFG:
      
        CALL NDFGLX(TRIA,IEL,1,IELTYP,KWORK(KVERT),KWORK(KMID),
     *              KDFG,KDFL)
        
C       Loop over the local DOF's to calculate the recovered gradient
C       there:

        DO IDG = 1,IDFL
        
C         What's the global DOF to I ?

          IVT = KDFG(IDG)
          
C         If IVT is a corner vertex that has already been tackled,
C         we don't have to do anything. It's enough to calculate
C         that value once :)

          IF ((IVT.LE.NVT).AND.(KWORK(KCOUNT+IVT-1).EQ.0)) THEN

C           Remember the node as being tackled and determine the 
C           coordinates of the vertex.
          
            KWORK(KCOUNT+IVT-1) = KWORK(KCOUNT+IVT-1) + 1
            CALL NDE2XY (IVT,TRIA,X,Y)

C           IVT must be any corner point of the element. So if IVT is
C           not an inner node...

            IF (IVT.GT.NVT) THEN
          
C             ...take an arbitrary corner instead of IVT, but don't change
C             the (X,Y)-coordinates! : IVT=KVERT(1,IEL)

              IVT = KWORK(KVERT+(IEL-1)*NNVE+1)
          
            END IF
          
C           IVT now points to a corner vertex of the current element.
C           Next we have to search for a (is possible: neighbour-) node
C           of IVT which is indide of the domain.
C          
C           Count the number of elements surrounding the current vertex:
 
            DO IVEL = 1, NVEL
              IF(KWORK(KVEL+(IVT-1)*NVEL + IVEL-1).EQ.0) GOTO 10
            ENDDO
10          CONTINUE
            IVEL = IVEL-1
          
C           If there are at least 3 elements surrounding our vertex IVT,
C           that's enough to recover the gradient:
          
            IF (IVEL.GE.3) THEN
              JVT = IVT
              CALL GRDLO2(TRIA,DU,DGX,DGY,
     *                  DWORK(KCORVG),KWORK(KVERT),KWORK(KMID),
     *                  KWORK(KVEL),X,Y,IVT,JVT,
     *                  IVEL,TRIA(ONVEL),ELE,
     *                  DWORK(L(LWORK)), LWORK)
            ELSE

C             Ok, too less elements, not so easy. Let's check if there
C             are at least two elements adjacent to our vertex IVT:

              IEL1 = KWORK(KVEL+(IVT-1)*NVEL)
              IEL2 = KWORK(KVEL+(IVT-1)*NVEL+1)

              IF (IEL2.NE.0) THEN
                 
C               We have two elements. Try to find a neighbour of IVT which
C               is an inner node. Use the KXNPR-array to decide about
C               inner nodes:
                 
                DO I = 1,NVE
                  DO J = 1,NVE
                    KV1 = KWORK(KVERT+(IEL1-1)*NNVE+I-1)
                    KV2 = KWORK(KVERT+(IEL2-1)*NNVE+J-1)
                    KNP = KWORK(KXNPR+2*(KV1-1))
                    IF ((KV1.EQ.KV2).AND.(IAND(KNP,2**8).EQ.0)) THEN
C                     We now know: Vertex KV1 is an inner vertex! (Has been
C                     calculated in the loop above). Take this as the new
C                     patch center JVT: 
                      JVT = KV1  
                      GOTO 20
                    END IF
                  END DO
                END DO

C               Oops, no inner node available! This can only happen on
C               rather coarse grids like the QUAD mesh - but it can 
C               happen. We have no chance, we have to go back to the
C               original evaluation routines to get a value vor the 
C               gradient :(
C
C               Indicate this by setting JVT=0

                JVT = 0

20              CONTINUE

C               ...and continue with the calculation below.
              
              ELSE
              
C               Hmmm, only one element on current node IVT. Loop 
C               through all corners to find a vertex that has at least
C               three adjacent elements

                DO I=1,TRIA(ONVE)
                  JVT = KWORK(KVERT+NNVE*(IEL-1)+I-1)
                  IF (KWORK(KVEL+NVEL*(JVT-1)+2).NE.0) THEN
C                   Node JVT is usable as new center.
                    GOTO 30
                  END IF
                END DO
                
C               Also here: no neighbour vertex with >= 2 elements
C               found. Switch back to standard reconstruction:

                JVT = 0
                
30              CONTINUE

              END IF

C             Is a node JVT found that we can use as a center of our
C             patch?

              IF (JVT.EQ.0) THEN
              
C               Unfortunately not. Then we have to fall back to standard
C               evaluation...
              
                CALL SCEVLQ (NEQ,DU,ELE,.FALSE.,X,Y,IEL,
     *                       TRIA,DWORK(KCORVG),KWORK(KVERT),
     *                       KWORK(KMID),D, DGX, DGY)
              
              ELSE
              
C               Nice, we can continue :)
C               Calculate the number of elements surrounding our
C               node JVT:
              
                DO IVEL = 1, NVEL
                  IF(KWORK(KVEL+(JVT-1)*NVEL + IVEL-1).EQ.0) GOTO 40
                ENDDO
40              CONTINUE
                IVEL = IVEL-1
          
C               And finally call the recovery routine to calculate
C               the gradient, using JVT as center...
              
                CALL GRDLO2(TRIA,DU,DGX,DGY,
     *                  DWORK(KCORVG),KWORK(KVERT),KWORK(KMID),
     *                  KWORK(KVEL),X,Y,IVT,JVT,
     *                  IVEL,TRIA(ONVEL),ELE,
     *                  DWORK(L(LWORK)), LWORK)

              END IF
              
            END IF
            
C           We now have any kind of gradient in DGX/DGY - either
C           calculated by ZZ or (in special cases) by direct
C           evaluation. Add this to the value in the DGRADX/DGRADY
C           vector. Take into consideration that our current node
C           is actually KDFG(IDFG), not IVT!
C           IVT might have been changed to a corner if KDFG(IDFG)
C           points to an inner/edge vertex/node.

            DGRADX (KDFG(IDG)) = DGRADX (KDFG(IDG)) + DGX
            DGRADY (KDFG(IDG)) = DGRADY (KDFG(IDG)) + DGY

          END IF
          
        END DO
        
      END DO
        
C Ok, everything calculated. Now we have to divide the calculated
C numbers by number of how often we tackled each node. This is of
C course only necessary for nodes corresponding to non-corner
C vertices, as corner vertices are only tackled once.

      DO I=NVT+1,NEQ
        D = 1D0/DBLE(KWORK(KCOUNT+I-1))
        DGRADX (I) = DGRADX (I) * D
        DGRADY (I) = DGRADY (I) * D
      END DO

C Release memory, that's it.
      
      CALL ZDISP (0,LCOUNT,'LCOUNT')
      CALL ZDISP (0,LWORK,'WORK  ')
      
      END
      
**************************************************************************
* Computation of reconstructed gradient with Standard-technique,
* Extended calling convenction
*
* This routine calculates for a given solution vector DU based on an
* element ELE the reconstructed gradient DGRADX,DGRADY in each node
* of the grid.
* The standard 1st-order method is used for reconstructing the
* gradient.
*
* In:
*   TRIA    - array [1..SZTRIA] of integer
*             Triangulation structure of the mesh
*   DU      - array [1..NEQ] of double
*             Solution vector
*   NEQ     - Number of equations in DU 
*   ELE     - Element that corresponds to solution vector DU.
*   BNONPR  - =true:  ELE is a nonparametric element
*             =false: ELE is a parametric element
* 
* Out:
*   DGRADX  - array [1..NEQ] of integer
*             X-Gradient of solution vector
*   DGRADY  - array [1..NEQ] of integer
*             Y-Gradient of solution vector
**************************************************************************

      SUBROUTINE XGRDR3(TRIA,DU,DGRADX,DGRADY,NEQ,ELE,BNONPR)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'stria.inc'
      
C parameters

      INTEGER NEQ,TRIA(SZTRIA)
      LOGICAL BNONPR
      DOUBLE PRECISION DU(NEQ),DGRADX(NEQ),DGRADY(NEQ)
      EXTERNAL ELE

      INTEGER NDFL
      EXTERNAL NDFL

C local variables

      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL,IELTYP
      INTEGER I
      
      INTEGER KVERT,KMID,KCOUNT,KCORVG
      INTEGER IEL,IVT,NEL,NVT,IDG
      
      INTEGER LCOUNT
      DOUBLE PRECISION D,X,Y,DGX,DGY
     
C We need an extended triangulation structure for this routine:

      IF (TRIA(OTRIFL).EQ.0) THEN
        WRITE (*,*) 'XGRDR3: No extended triangulation information!'
        RETURN
      END IF
     
C Clear the destination arrays at first:

      CALL LCL1 (DGRADX,NEQ)
      CALL LCL1 (DGRADY,NEQ)
      
C How do we have to perform the recovering of the gradient now?
C Well, the standard implementation of calculating the gradient
C uses the derivative-calculation methods of the element itself,
C adds  everything together and takes the mean about how many
C times each point is tackled.
C
C For calculating the mean we need an array KCOUNT. This saves for
C every DOF how many times the recovered gradient has been calculated
C for that DOF. Afterwards it's used to take the mean.

      CALL ZNEW (NEQ,3,LCOUNT,'KCOUNT')

C Transfer some adresses to local variables for easier access

      KCORVG = L(TRIA(OLCORVG))
      KVERT = L(TRIA(OLVERT))
      KMID  = L(TRIA(OLMID))
      NVT = TRIA(ONVT)
      NEL = TRIA(ONEL)
      
      KCOUNT = L(LCOUNT)
      
C Ask the element about its type:

      IELTYP = -1
      CALL ELE (0D0,0D0,IELTYP)

C Calculate local number of degrees of freedom of our element type
      
      IDFL = NDFL(IELTYP)
      
C Now perform a loop over all elements:

      DO IEL = 1,NEL
      
C       On the current element calculate the global DOF's into KDFG:
      
        CALL NDFGLX(TRIA,IEL,1,IELTYP,KWORK(KVERT),KWORK(KMID),
     *              KDFG,KDFL)
        
C       Loop over the local DOF's to calculate the recovered gradient
C       there:

        DO IDG = 1,IDFL
        
C         What's the global DOF to IDG ?

          IVT = KDFG(IDG)
          
C         Remember the node as being tackled and determine the 
C         coordinates of the vertex.
          
          KWORK(KCOUNT+IVT-1) = KWORK(KCOUNT+IVT-1) + 1
          CALL NDE2XY (IVT,TRIA,X,Y)

C         Use the standard evaluation routine to calculate an
C         approximate gradient in the point (X,Y). The X-derivative
C         will be saved to DGX, the Y-derivative to DGY. The value
C         of the FE-function in this point, which is saved to D,
C         will be ignored in the following.
              
          CALL SCEVLQ (NEQ,DU,ELE,BNONPR,X,Y,IEL,
     *                 TRIA,DWORK(KCORVG),KWORK(KVERT),KWORK(KMID),
     *                 D, DGX, DGY)
              
C         Add the gradient values to the value in the DGRADX/DGRADY
C         vector. 

          DGRADX (IVT) = DGRADX (IVT) + DGX
          DGRADY (IVT) = DGRADY (IVT) + DGY
          
        END DO
        
      END DO
        
C Ok, everything calculated. Now we have to divide the calculated
C numbers by number of how often we tackled each node. 
C In contrast to the ZZ-technique this is necessary for inner
C nodes as well as boundary nodes.

      DO I=1,NEQ
        D = 1D0/DBLE(KWORK(KCOUNT+I-1))
        DGRADX (I) = DGRADX (I) * D
        DGRADY (I) = DGRADY (I) * D
      END DO

C Release memory, that's it.
      
      CALL ZDISP (0,LCOUNT,'LCOUNT')
      
      END

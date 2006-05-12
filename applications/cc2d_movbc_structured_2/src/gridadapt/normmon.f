***********************************************************************
* This file contains for the calculation of weights for the monitor
* function as well as normalising routines.
***********************************************************************

***********************************************************************
* Explaination of the method:
*
* The dynamic approach of the grid adaption method tries to find
* a PHI:Omega->Omega that way, that:
*
*    det grad(PHI(x)) = f(PHI(x))
*
* for a given monitor function f:Omega->Omega. f has to fulfill the
* compatibility equation:
*
*    int_Omega ( 1/f(z) ) dz = |Omega|
*
* Our method here requires an arbitrary monitor function with values
* f(z) in (0,1]. Therefore we have to scale the monitor function
* according to the compatibility equation.
***********************************************************************

***********************************************************************
* Normalize monitor function to Int(1/f)= Domain measure
*
* Wrapper-function.
*
* In:
*  NVT    - Number of vertices
*  LONE   - Handle to: array [1..NVT] of dobule
*           Mean area around each vertex
*  DMSR   - Total measure of the domain
*  LMON   - Handle to monitor function
*
* Out:
*  Normalised monitor function.
***********************************************************************

      SUBROUTINE  XNORM(LONE,NVT,DMSR,LMON)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      DOUBLE PRECISION DMSR
      INTEGER LONE,LMON,NVT

      CALL NORM(DWORK(L(LONE)),NVT,DMSR,DWORK(L(LMON)))
      
      END 

***********************************************************************
* Normalize monitor function to Int(1/f)= Domain measure
*
* This routine performs a quick-and-dirty inexact integration using 
* the mean area around each vertex.
*
* In:
*  NVT    - Number of vertices
*  DONE   - array [1..NVT] of dobule
*           Mean area around each vertex
*  DMSR   - Total measure of the domain
*  DMON   - array [1..NVT] of double
*           Monitor function
*
* Out:
*  DMON   - Normalised monitor function.
***********************************************************************
* Assume an arbitrary function g:Omega->R. The aim is to
* obtain a monitor function f:Omega->R that way, that
*
*    int_Omega 1/f = |Omega|
*
* holds. The approach is to scale g by a constant:    f = c * g
* Let {x} in Omega be the sei of points where our monitor function g
* is given, and let Area(x) be the mean area around each x. Then
* the compatibility equation can be approximated as:
*
*    |Omega| = int_Omega 1/f = sum_x ( 1/f(x) * Area(x) )
*                            = sum_x ( Area(x) / f(x) )
*                            = sum_x ( Area(x) / (c*g(x)) )
*                            = 1/c * sum_x ( Area(x) / g(x) )
*
* i.e.:   c = sum_x ( Area(x) / g(x) )  /  |Omega|
* 
* That formula is implemented below...
***********************************************************************

      SUBROUTINE NORM(DONE,NVT,DMSR,DMON)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      
      DOUBLE PRECISION DMON,DONE,DMSR
      INTEGER NVT
      DIMENSION DMON(NVT),DONE(NVT)
      DOUBLE PRECISION DSUM
      INTEGER IVT
      
      IF (MT.GE.3) THEN
        WRITE (MTERM,FMT='(A,F10.6)') 
     *         'Normalizing monitor function to ',DMSR
      END IF
      
C The old monitor function is g=DMON.
C build:    DSUM = sum_x ( Area(x) / g(x) )
      
      DSUM=0D0
      DO IVT=1,NVT
        DSUM=DSUM+DONE(IVT)/DMON(IVT)
      END DO
      
C build:    c = DSUM / |Omega|
      
      DSUM=DSUM/DMSR
      
C scale g(x) to obtain the actual monitor function f:
C    f(x) = c*g(x)
C Overwrite the old monitor function.
      
      DO IVT=1,NVT
        DMON(IVT)=DMON(IVT)*DSUM
      END DO
!      DSUM=0D0
!      DO 30 IVT=1,NVT
!30    DSUM=DSUM+DONE(IVT)/DMON(IVT)

      END 

***********************************************************************
* Normalize monitor function to Int(f)= Domain measure
*
* This routine performs a quick-and-dirty inexact integration using 
* the mean area around each vertex.
*
* In:
*  NVT    - Number of vertices
*  DONE   - array [1..NVT] of dobule
*           Mean area around each vertex
*  DMSR   - Total measure of the domain
*  DMON   - array [1..NVT] of double
*           Monitor function
*
* Out:
*  DMON   - Normalised monitor function.
*
* Remark:
*    DONE = DMON is allowed.
***********************************************************************
* Assume an arbitrary function g:Omega->R. The aim is to
* obtain a monitor function f:Omega->R that way, that
*
*    int_Omega f = |Omega|
*
* holds. The approach is to scale g by a constant:    f = c * g
* Let {x} in Omega be the set of points where our monitor function g
* is given, and let Area(x) be the mean area around each x. Then
* the compatibility equation can be approximated as:
*
*    |Omega| = int_Omega f   = sum_x ( f(x) * Area(x) )
*                            = sum_x ( (c*g(x)) * Area(x) )
*                            = c * sum_x ( g(x) * Area(x) )
*
* i.e.:   c = |Omega| / sum_x ( Area(x) * g(x) )   
* 
* That formula is implemented below...
***********************************************************************

      SUBROUTINE NORMIF(DONE,NVT,DMSR,DMON)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      
      DOUBLE PRECISION DMON,DONE,DMSR
      INTEGER NVT
      DIMENSION DMON(NVT),DONE(NVT)
      DOUBLE PRECISION DSUM
      INTEGER IVT
      
C The old monitor function is g=DMON.
C build:    DSUM = sum_x ( Area(x) * g(x) )
      
      DSUM=0D0
      DO IVT=1,NVT
        DSUM = DSUM + DONE(IVT)*DMON(IVT)
      END DO
      
C build:    c = |Omega| / DSUM 
      
      DSUM=DMSR / DSUM
      
C scale g(x) to obtain the actual monitor function f:
C    f(x) = c*g(x)
C Overwrite the old monitor function.
      
      DO IVT=1,NVT
        DMON(IVT)=DMON(IVT)*DSUM
      END DO

      END 

**********************************************************************
* Normalise monitor function
*
* Normally the scaled monitor function is calculated around
* each node by FMON_i/AREA_i with AREA_i the average area of the
* elements around node i. For the calculation of the grid deformation
* the monitor function is only needed up to a constant, which is
* unappropriate for visualisation and/or error calculation.
* This function scales the monitor function such that the mean
* value is 1D0 (and every entry is positive).
* (For this purpose the normalising routines of the Neumann filter
*  are used)
*
* In:
*  DMON   - array [1..NVT] of double
*           Given scaled monitor function
*  NVT    - Number of vertices
*
* Out:
*  DMONNM - array [1..NVT] of double
*           Normalised scaled monitor function
*
* DMONNM=DMON is allowed.
**********************************************************************

      SUBROUTINE NRMMON (DMON,NVT,DMONNM)
      
      IMPLICIT NONE
      
      INTEGER NVT
      DOUBLE PRECISION DMON(NVT),DMONNM(NVT)
      
C local variables

      INTEGER I
      DOUBLE PRECISION MEAN,DCLIMV 
      
C Calculate the current mean value
      
      MEAN = DCLIMV (DMON,NVT)
      
C Divide all vector entries by the calculated meanvalue to get an
C all positive vector with mean value 1d0

      DO I=1,NVT
        DMONNM(I)=DMON(I)/MEAN
      END DO
      
      END
      
************************************************************************
* Calculate the weight of a monitor function
*
* This is the maximum ratio between the cell size and the corresponding
* function that determines the grid deformation. THe value is reduced
* by 1 to make it 0-based
*
* Remember: DMON(I) = F(I) / AREA(I)
* With F(I)    = function determining the grid deformation in node I,
*      AREA(I) = mean cell size around node I
* and: 
*  if DMON(I) > 1   =>   D(I) > AREA(I)   =>   ratio = DMON(I)
*  if DMON(I) < 1   =>   D(I) < AREA(I)   =>   ratio = 1/DMON(I)
*
* Example:
*  DMON(I)=100 
*    => Ratio = 100 because Mean cell size is 100 times to large
*  DMON(I)=0.1
*    => Ratio = 1 / 0.1 = 10 because Mean cell size is 10 times to small
*
* In:
*  NVT     - Number of vertices
*  DMON    - array [1..NVT] of double
*            Monitor function to analyse
*
* Out:
*  WEIGHT  - The weight of the monitor function = ratio, >= 1
*  INODE   - The number of the node with the above weight
************************************************************************

      SUBROUTINE MONWGH (DMON,NVT,WEIGHT,INODE)
      
      IMPLICIT NONE
      
      INTEGER NVT, INODE
      DOUBLE PRECISION DMON(NVT),WEIGHT
      
      INTEGER I
      
      DOUBLE PRECISION CWGH
      
      WEIGHT=0D0
      INODE=0
      
      DO I=1,NVT
      
C What is the current weight?
      
        CWGH=DMON(I)
        IF (CWGH.LT.1D0) CWGH=1D0/CWGH
        
        IF ( CWGH.GT.WEIGHT ) THEN
          WEIGHT = CWGH
          INODE=I
        END IF
      END DO

      END 

************************************************************************
* Calculate the cell ratio from the monitor function
*
* Calculates with the help of the current monitor function 
* 1.) the ratio between the smallest and largest cell or
* 2.) the largest monitor function weight
*
* In:
*  NVT     - Number of vertices
*  DMON    - array [1..NVT] of double
*            Monitor function to analyse
*  IMTH    - type of return value:
*            =0: return the ratio between largest and smallest cell
*                in RAT, calculated with the help of the monitor
*                function
*            =1: return largest f_i/area_i (or area_i/f_i if
*                f_i/area_i < 1)
*
* Out:
*  RAT     - ratio
*  INDEMN  - The number of the node with the minimum weight
*  INDEMX  - The number of the node with the maximmum weight
************************************************************************

      SUBROUTINE MONCRT (DMON,NVT,IMTH,RAT,INDEMN,INDEMX)
      
      IMPLICIT NONE
      
      INTEGER NVT, INDEMN, INDEMX, IMTH
      DOUBLE PRECISION DMON(NVT),RAT
      
      INTEGER I
      
      DOUBLE PRECISION MN,MX
      
      RAT = 0D0
      INDEMN = 1
      INDEMX = 1
      MN = DMON(1)
      MX = DMON(1)
 
C Determine minimum and maximum weight = cell-ratio
      
      DO I=1,NVT
        IF ( DMON(I).LT.MN ) THEN
          MN = DMON(I)
          INDEMN = I
        END IF
        IF ( DMON(I).GT.MX ) THEN
          MX = DMON(I)
          INDEMX = I
        END IF
      END DO

      IF (IMTH.EQ.0) THEN

C calculate the ratio

        RAT = MX / MN
        
      ELSE
      
C Calculate the largest monitor function value; invert it
C if it is < 1 to obtain a value > 1!

        IF (MX.LT.1D0) MX=1/MX
        IF (MN.LT.1D0) MX=1/MN
        
        RAT = MAX (MN,MX)
      
      END IF

      END 

************************************************************************
* Monitor function rescaling
*
* Rescales the monitor function. Every entry DMON(I) is taken up to
* a power of PW.
*
* In:
*  NVT     - Number of vertices
*  DMON    - array [1..NVT] of double
*            Monitor function to analyse
*  PW      - Power that should be applied to each entry
*
* Out:
*  DMON    - Rescaled monitor function
*
* The monitor function is rescaled during this process such that the
* minimum value is 1D0!
************************************************************************

      SUBROUTINE MONRSC (DMON,NVT,PW)
      
      IMPLICIT NONE
      
      INTEGER NVT
      DOUBLE PRECISION DMON(NVT),PW
      
      INTEGER I
      
      DO I=1,NVT
        DMON(I) = DMON(I) ** PW
      END DO

      END 

**********************************************************************
* Normalize monitor function to Int(f)= Domain measure
*
* This routine performs exact integration to normalize the monitor
* function. It scales f by multiplication of a constant to enforce
*     int_Omega f = |Omega|
*
* In:
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure
*  DMON   - array [1..NVT] of double
*           Monitor function
*  DMSR   - Total measure of the domain
*  NEQ    - number of elements in DMON
*  ELE    - Element that should be used for the interpretation of
*           DMON as FE-function. Must be a conforming element.
*  ICUB   - Number of cubature formula that should be used
*           for integration
*
* Out:
*  DMON   - Normalised monitor function.
**********************************************************************

      SUBROUTINE NRMM2(TRIA,DMON,DMSR,NEQ,ELE,ICUB)

      IMPLICIT NONE

      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'

      INCLUDE 'cbasictria.inc'

      INCLUDE 'stria.inc'

      INCLUDE 'ccub.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
C     parameters:
      
      INTEGER TRIA(SZTRIA),NEQ,ICUB
      DOUBLE PRECISION DMON(NEQ),DMSR
      
C     local variables:

      INTEGER IELTYP,IDFL,I,IG,JP,IEL1,IVE,NVE
      DOUBLE PRECISION DJF(2,2),OM,DINTF,XI1,XI2
      DOUBLE PRECISION DCOORD(2,4)
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
      INTEGER KVERT, KMID

C     externals

      EXTERNAL ELE
      INTEGER NDFL
      EXTERNAL NDFL
      
      IER = 0
      
C     The formula is as follows. Let g:=c*f denote the scaled monitor
C     function f. Then:
C       int_Omega g = int_Omega c*f = c*int_Omega f = |Omega|
C     So,
C       c = |Omega| / int_Omega f
C
C     At first we have to build int_Omega f.
C
C     Resolve some handles for quicker access:
      
      KVERT = L(TRIA(OLVERT))
      KMID = L(TRIA(OLMID))
      NVE = TRIA(ONVE)
      
C     We need only function values for the integration:

      DO  I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      BDER(1)=.TRUE.

C     Get the element type as well as number of local degrees of
C     freedom:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      
C     initialise cubature formula
        
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999

C     Prepare conforming element for cubature. Save the current cubature
C     formula into the COMMON block variable ICUBP - either for now or
C     for later use. Perform a dummy call to the element in the conforming
C     case.
C     This is done for saving arithmetic operations in later calls.

      ICUBP=ICUB
      
      CALL ELE(0D0,0D0,-2)
        
C     Loop about the elements, calculate the integral into DINTF

      DINTF = 0D0
      
      DO IEL1 = 1,TRIA(ONEL)
      
C       Calculate the local and global DOF's on our current element.
C       We later have to loop about them...

        CALL NDFGL(IEL1,1,IELTYP,KWORK(KVERT),KWORK(KMID), KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999
      
C Store the coordinates of the corners of the current element in 
C the COMMON block variables DX/DY (necessary for the element) 
C as well as in the DCOORD array for the coordinate transformation.

        DO IVE = 1, NVE
          JP=KWORK(KVERT + (IEL1-1)*NNVE+(IVE-1))
          KVE(IVE)=JP
          CALL NDE2XY (JP,TRIA,DCOORD (1,IVE),DCOORD (2,IVE))
          DX(IVE)=DCOORD (1,IVE)
          DY(IVE)=DCOORD (2,IVE)
        END DO
        
C Initialise auxiliary Jacobian factors DJF for transformation

        CALL QINIJF (DCOORD,DJF)

c Loop over all cubature points

        DO ICUBP = 1, NCUBP
          
c Cubature point on the reference element

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
          
C         Calculate Jacobian matrix of transformation and its
C         determinant. This is necessary for the element routine ELE
C         to properly calculate derivatives. The result is directly
C         written into the element COMMON block variables DJAC and DETJ.
 
          CALL QTRDET (DCOORD,DJF,DJAC,DETJ,XI1,XI2)

C         Evaluate the basis functions in the cubature point

          CALL ELE(XI1,XI2,-3)
          IF (IER.LT.0) GOTO 99999
      
C         Calculate the weighting factor for the current cubature point
C         with the help of the Jacobian determinant

          OM = DOMEGA(ICUBP)*DETJ
          
C         Add the values of the basis function to the integral
      
          DO I=1,IDFL
            IG=KDFG(I)
            DINTF = DINTF + OM * DMON(IG)*DBAS(KDFL(I),1)
          END DO
          
        END DO
      
      END DO

C     DINTF now contains the integral. Using c=|Omega|/DINTF,
C     scale the monitor function appropriately.

      IF (DINTF.NE.0D0) THEN

        DINTF = DMSR/DINTF
        DO I=1,NEQ
          DMON(I) = DMON(I)*DINTF
        END DO
        
      END IF
      
C     ... or do nothing if an error occurred

99999 CONTINUE      
      
      END  

**********************************************************************
* Normalize monitor function to Int(1/f)= Domain measure
*
* This routine performs exact integration to normalize the monitor
* function. It scales f by multiplication of a constant to enforce
*     int_Omega 1/f = |Omega|
*
* In:
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure
*  DMON   - array [1..NVT] of double
*           Monitor function
*  DMSR   - Total measure of the domain
*  NEQ    - number of elements in DMON
*  ELE    - Element that should be used for the interpretation of
*           DMON as FE-function. Must be a conforming element.
*  ICUB   - Number of cubature formula that should be used
*           for integration
*
* Out:
*  DMON   - Normalised monitor function.
**********************************************************************

      SUBROUTINE NRMM2R(TRIA,DMON,DMSR,NEQ,ELE,ICUB)
      IMPLICIT NONE

      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'

      INCLUDE 'cbasictria.inc'

      INCLUDE 'stria.inc'

      INCLUDE 'ccub.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
C     parameters:
      
      INTEGER TRIA(SZTRIA),NEQ,ICUB
      DOUBLE PRECISION DMON(NEQ),DMSR
      
C     local variables:

      INTEGER IELTYP,IDFL,I,IG,JP,IEL1,IVE,NVE
      DOUBLE PRECISION DJF(2,2),OM,DINTF,XI1,XI2,DF
      DOUBLE PRECISION DCOORD(2,4)
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
      INTEGER KVERT, KMID

C     externals

      EXTERNAL ELE
      INTEGER NDFL
      EXTERNAL NDFL
      
C     A first idea for this routine is the following:
C     1.) Transform f -> 1/f 
C     2.) Use NRMM2 to normalize 1/f
C     3.) Transform f -> 1/f for obtaining the result
C     This strategy is possible, but slightly wrong! The reason is, that
C     1/f is not a bilinear function, and so normalizing 1/f would
C     not mean for f that int(1/f)=|Omega|.
C
C     Instead we use a different normalization strategy. We use the
C     same integration routine as stated above, but with a slightly
C     different summing in the coefficients of the integral to
C     properly resolve 1/f. 
      
      IER = 0
      
C     The formula is as follows. Let g:=c*f denote the scaled monitor
C     function f. Then:
C       int_Omega g = int_Omega c*f = c*int_Omega f = |Omega|
C     So,
C       c = |Omega| / int_Omega f
C
C     At first we have to build int_Omega f.
C
C     Resolve some handles for quicker access:
      
      KVERT = L(TRIA(OLVERT))
      KMID = L(TRIA(OLMID))
      NVE = TRIA(ONVE)
      
C     We need only function values for the integration:

      DO  I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      BDER(1)=.TRUE.

C     Get the element type as well as number of local degrees of
C     freedom:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      
C     initialise cubature formula
        
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999

C     Prepare conforming element for cubature. Save the current cubature
C     formula into the COMMON block variable ICUBP - either for now or
C     for later use. Perform a dummy call to the element in the conforming
C     case.
C     This is done for saving arithmetic operations in later calls.

      ICUBP=ICUB
      
      CALL ELE(0D0,0D0,-2)
        
C     Loop about the elements, calculate the integral into DINTF

      DINTF = 0D0
      
      DO IEL1 = 1,TRIA(ONEL)
      
C       Calculate the local and global DOF's on our current element.
C       We later have to loop about them...

        CALL NDFGLX(TRIA,IEL1,1,IELTYP,KWORK(KVERT),KWORK(KMID), 
     *              KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999
      
C Store the coordinates of the corners of the current element in 
C the COMMON block variables DX/DY (necessary for the element) 
C as well as in the DCOORD array for the coordinate transformation.

        DO IVE = 1, NVE
          JP=KWORK(KVERT + (IEL1-1)*NNVE+(IVE-1))
          KVE(IVE)=JP
          CALL NDE2XY (JP,TRIA,DCOORD (1,IVE),DCOORD (2,IVE))
          DX(IVE)=DCOORD (1,IVE)
          DY(IVE)=DCOORD (2,IVE)
        END DO
        
C Initialise auxiliary Jacobian factors DJF for transformation

        CALL QINIJF (DCOORD,DJF)

c Loop over all cubature points

        DO ICUBP = 1, NCUBP
          
c Cubature point on the reference element

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
          
C         Calculate Jacobian matrix of transformation and its
C         determinant. This is necessary for the element routine ELE
C         to properly calculate derivatives. The result is directly
C         written into the element COMMON block variables DJAC and DETJ.
 
          CALL QTRDET (DCOORD,DJF,DJAC,DETJ,XI1,XI2)

C         Evaluate the basis functions in the cubature point

          CALL ELE(XI1,XI2,-3)
          IF (IER.LT.0) GOTO 99999
      
C         Add the values of the basis function to calculate the value
C         of the function in the current cubature point.
      
          DF = 0D0
          DO I=1,IDFL
            IG=KDFG(I)
            DF = DF + DMON(IG)*DBAS(KDFL(I),1)
          END DO
          
C         Calculate the weighting factor for the current cubature point
C         with the help of the Jacobian determinant

          OM = DOMEGA(ICUBP)*DETJ
          
C         Now add 1/DF * cubature weight to the integral. This properly
C         calculates int(1/f):

          DINTF = DINTF + OM/DF

        END DO
      
      END DO

C     DINTF now contains the integral. Using c=DINTF/|Omega|,
C     scale the monitor function appropriately:
C
C        int(1/(c*DMON)) = |Omega|
C     => c = 1/|Omega| * int(1/DMON) = DINTF/|Omega|

      IF (DMSR.NE.0D0) THEN

        DINTF = DINTF/DMSR
        DO I=1,NEQ
          DMON(I) = DMON(I)*DINTF
        END DO
        
      END IF
      
C     ... or do nothing if an error occurred

99999 CONTINUE      

      END
      
************************************************************************
* Calculate quality measure
*
* This routine calculates the quality measure of a deformed grid.
* The quality measure is defined by:
*
*      Q := ||q-1||_l2
*
* with q(x) = f(x)/a(x), f(.)=monitor function, a(.)=element area
* distribution on deformed mesh. 
*
* In:
*   TRIA    - array [1..SZTRIA] of integer
*             Triangulation structure of the deformed grid
*   DMON    - array [1..NVT] of double
*             monitor function; must correspond to TRIA
*   DSIZE   - array [1..NVT] of double
*             current area distribution around each vertex
*   DMSR    - Measure of the domain
*
* Out:
*   QMSR    - Quality measure
************************************************************************

      SUBROUTINE GAQMSR (NVT,DMON,DSIZE,DMSR,QMSR)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'stria.inc'
      
      INTEGER NVT
      DOUBLE PRECISION DMON(NVT),DSIZE(NVT),DMSR,QMSR
      
C local variables

      INTEGER LTMP,KTMP,KSIZE,I

C Make a copy of DMON and DSIZE for normalizing

      CALL ZNEW (2*NVT,-1,LTMP,'DTMP  ')
      CALL LCP1 (DMON,DWORK(L(LTMP)),NVT)
      CALL LCP1 (DSIZE,DWORK(L(LTMP)+NVT),NVT)
      
C Normalize monitor function and area distribution such that
C         Int(f) = Int(a) ( = |Omega| )
C This way if f~a, then q(x)~1, so Q=0.
      
      CALL NORMIF(DSIZE,NVT,DMSR,DWORK(L(LTMP)))
      CALL NORMIF(DSIZE,NVT,DMSR,DWORK(L(LTMP)+NVT))
      
C Build q(x)-1:

      KTMP = L(LTMP)
      KSIZE = L(LTMP)+NVT
      
      DO I=0,NVT-1
        DWORK(KTMP+I) = DWORK(KTMP+I)/DWORK(KSIZE+I) - 1D0
      END DO      

C Build the l2-norm of the vector

      CALL LL21(DWORK(KTMP),NVT,QMSR)
      QMSR = QMSR / SQRT(DBLE(NVT))
      
C Release memory, finish
      
      CALL ZDISP (0,LTMP,'DTMP  ')
      
      END
      
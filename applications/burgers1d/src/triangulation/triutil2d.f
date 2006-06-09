************************************************************************
* This file contains various auxiliary routine for maintaining
* grid information.
************************************************************************

************************************************************************
* Average element midpoint coordinates
*
* The follwing function corrects the grid on domains where non-linear/
* curved boundary segments are used. 
* When a boundary segment of the domain is a line, new boundary nodes
* are automatically positioned on the boundary. But when the boundary
* is a curve, a new boundary vertex is not automatically on the
* boundary, but it first has tobe moved to there. This is already
* performed in the refinement routine XSB0X. Unfortunately this
* procedure can lead to very anisotropic elements near the boundary,
* depending on how sharp the curved boundary is.
*
* AVEMPC now tries to reduce these effects of anisotropy. The element
* midpoint of the coarser element (which is at the same time the vertex
* where all the four finer elements meet) is taken as the average of
* the four edge midpoints that arise from natural refinement.
*
* o----------o----------o
* \          \          |
* |\         \          | 
* | \      -> \         |
* |  o---------o--------o   Averaging of the element midpoint by 
* | /      -> /         |   interpolation
* |/         /          |
* /          /          |
* o----------o----------o
*
* In:
*   DCORVG,
*   KVERT,
*   KADJ,
*   NEL     - Gemeotry information at level >= 2. DCORVG must contain
*             the coordinates of the vertices at level 1..current level
*
* Out:
*   DCORVG  - The coordinates of the element corners that correspond
*             to coarse grid element midpoints are corrected.
************************************************************************

      SUBROUTINE AVEMPC(DCORVG,KVERT,KADJ,NEL)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C parameters
      
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER KVERT(NNVE,*),KADJ(NNVE,*)
      INTEGER NEL

C local variables

      INTEGER IADJ3, IADJ4
      INTEGER IVT1, IVT2, IVT3, IVT4
      INTEGER IVTM, IEL
      DOUBLE PRECISION PX1,PX2,PX3,PX4,PY1,PY2,PY3,PY4,PXM,PYM,PX,PY

      DO IEL=1,NEL/4

        IADJ3=KADJ(2,IEL)
        IADJ4=KADJ(3,IEL)
        IVT1=KVERT(2,IEL)
        IVT2=KVERT(4,IEL)
        IVT3=KVERT(2,IADJ3)
        IVT4=KVERT(4,IADJ4)
        IVTM=KVERT(3,IEL)

        PX1=DCORVG(1,IVT1)
        PX2=DCORVG(1,IVT2)
        PX3=DCORVG(1,IVT3)
        PX4=DCORVG(1,IVT4)

        PY1=DCORVG(2,IVT1)
        PY2=DCORVG(2,IVT2)
        PY3=DCORVG(2,IVT3)
        PY4=DCORVG(2,IVT4)

        PXM=DCORVG(1,IVTM)
        PYM=DCORVG(2,IVTM)

        PX=0.25D0*(PX1+PX2+PX3+PX4)
        PY=0.25D0*(PY1+PY2+PY3+PY4)

        DCORVG(1,IVTM)=PX
        DCORVG(2,IVTM)=PY

      END DO

      END

************************************************************************
* Grid disturbance
*
* This routine stochastically disturbes the current grid.
* Every gridpoint is moved DIEPS percent in a - more or less random -
* direction. The grid is assumed to be uniform, the minimum cell
* size is computed with the help of NVT!
*
* In:
*  DCORVG, KNPR, NVT - as usual
*  DIEPS             - Rate of disturbing. 0.2D0 = 20%
*
* The current grid will be modified directly.
************************************************************************

      SUBROUTINE GRDIST(DCORVG,KNPR,NVT,DIEPS)
      
      IMPLICIT NONE
      
C parameters
      
      DOUBLE PRECISION DCORVG(2,*),DIEPS
      INTEGER KNPR(*),NVT
      
C local variables

      DOUBLE PRECISION H, HDIST
      INTEGER IVT
      
      H=1D0/(SQRT(DBLE(NVT))-1)
      HDIST=DIEPS*H

      DO IVT=1,NVT
        IF (KNPR(IVT).EQ.0) THEN
          DCORVG(1,IVT)=DCORVG(1,IVT)+DBLE((-1)**MOD(IVT,17))*HDIST
          DCORVG(2,IVT)=DCORVG(2,IVT)+DBLE((-1)**IVT)*HDIST
        ENDIF
      END DO

      END


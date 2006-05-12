************************************************************************
* This file contains correction routines for the geometry.
* They must be used to ensure a valid grid if the grid is not
* deformed with "standard" parameters
* (e.g. performing adaption on level NLMAX-1, instead of NLMAX,
* recalculation of the midpoints,...)
************************************************************************

************************************************************************
* Recalculate midpoint coordinates
* 
* This routine uses the coarse grid information id DCORVG to build
* up the information about edge midpoints DCORMG.
*
* In:
*  DCORVG - Corner points of current triangulation
*  KVERT  - Array with number of vertices on each element
*  KMID   - Array with number of midpoints on each element
*  NEL    - Number of elements
*  NVT    - Number of vertices
*  NMT    - Number of midpoints
*
* Out:
*  DCORMG - Array with edge midpoints, recalculated from DCORVG
*
* This only works with quadrilaterals.
************************************************************************

      SUBROUTINE RXMPTS (DCORVG,KVERT,KMID,NEL,NVT,NMT,DCORMG)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C parameters
      
      INTEGER NVT,NMT,NEL,KVERT(NNVE,NEL),KMID(NNVE,NEL)
      DOUBLE PRECISION DCORVG(2,NVT),DCORMG(2,NMT)
      
C local variables

      INTEGER IEL,IV
      DOUBLE PRECISION XM,YM,X1,Y1,X2,Y2
      
C loop through all elements, rebuild DCORMG
      
      IF (NMT.LE.0) RETURN
      
      DO IEL=1,NEL
        DO IV=1,4
          X1 = DCORVG(1,KVERT(IV,IEL))
          Y1 = DCORVG(2,KVERT(IV,IEL))
          X2 = DCORVG(1,KVERT(MOD(IV,4)+1,IEL))
          Y2 = DCORVG(2,KVERT(MOD(IV,4)+1,IEL))
          XM = 0.5D0*(X1+X2)
          YM = 0.5D0*(Y1+Y2)
          DCORMG(1,KMID(IV,IEL)-NVT) = XM
          DCORMG(2,KMID(IV,IEL)-NVT) = YM
        END DO
      END DO
      
      END
      
************************************************************************
* Project coarse grid to fine grid
* 
* Recalculates the information about the corner vertices on a regularly
* refined mesh by the vertices on the coarse mesh.
*
* This routine only projects the coordinates of all points from the
* coarse grid to the fine grid. The parameter values of the boundary
* nodes are not calculated.
*
* In:
*  DCRVGC - Corner points of current triangulation, coarse grid
*  KVERTC - Array with number of vertices on each element, coarse grid
*  NELC   - Number of elements, coarse grid
*
*  DCRVGF - Corner points of current triangulation, fine grid
*  KVERTF - Array with number of vertices on each element, fine grid
*  KADJF  - Information about neighbor cells on the fine grid
*  NELF   - Number of elements, fine grid
*
* Out:
*  DCRVGF - New calculated corner points of current triangulation
*           on the fine grid
*
* This only works with quadrilaterals.
************************************************************************

      SUBROUTINE PRC2FG (DCRVGC,KVERTC,NELC,
     *                   DCRVGF,KVERTF,KADJF,NELF)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C parameters
      
      INTEGER NELC,NELF
      INTEGER KVERTC(NNVE,NELC),KVERTF(NNVE,NELF),KADJF(NNVE,*)
      DOUBLE PRECISION DCRVGC(2,*),DCRVGF(2,*)
      
C Externals

      DOUBLE PRECISION PARX,PARY
      EXTERNAL PARX,PARY
      
C local variables

      INTEGER IEL,IV,IVE(4),IVM,IADJ1,IADJ2,IADJ3
      DOUBLE PRECISION XM,YM,X1,Y1,X2,Y2,XC,YC
      
C loop through all elements in the coarse grid, rebuild DCRVGF
      
      DO IEL=1,NELC
        
C Look into the fine grid, compute the vertex numbers for the
C edge midpoints and the center.
C The element number on the coarse grid is the element number
C of the lower left element on the fine grid
        
C Right neighbor
        IADJ1 = KADJF(2,IEL)
C Top neighbor
        IADJ2 = KADJF(3,IEL)
C Top right neighbor
        IADJ3 = KADJF(2,IADJ1)
C Calculate vertex numbers
        IVE(1) = KVERTF(2,IEL)
        IVE(2) = KVERTF(2,IADJ1)
        IVE(3) = KVERTF(2,IADJ3)
        IVE(4) = KVERTF(2,IADJ2)
        IVM = KVERTF(3,IEL)
      
C Correct the edge midpoints
      
        XC = 0D0
        YC = 0D0
        
        DO IV=1,4
          X1 = DCRVGC(1,KVERTC(IV,IEL))
          Y1 = DCRVGC(2,KVERTC(IV,IEL))
          X2 = DCRVGC(1,KVERTC(MOD(IV,4)+1,IEL))
          Y2 = DCRVGC(2,KVERTC(MOD(IV,4)+1,IEL))

C Calculate edge midpoint

          XM = 0.5D0*(X1+X2)
          YM = 0.5D0*(Y1+Y2)

          DCRVGF(1,IVE(IV)) = XM
          DCRVGF(2,IVE(IV)) = YM
          
          XC = XC+XM
          YC = YC+YM
        END DO
        
C Correct the center

        DCRVGF(1,IVM) = 0.25D0*XC
        DCRVGF(2,IVM) = 0.25D0*YC

      END DO
      
C Perform some additional correction of the element midpoints
C (for nonlinear boundaries)

      CALL AVEMPC(DCRVGF,KVERTF,KADJF,NELF)

      END

************************************************************************
* Project coarse grid boundary parameters to fine grid boundary param.
* 
* Recalculates the information about the parametrisation of the
* fine grid boudary. Copies the parameter values of the coarse grid
* boundary nodes to the fine grid boundary nodes. Recalculates
* the new fine grid points by interpolation of the parameters of
* two coarse grid boundary nodes
*
* For now this works only if the grid fulfills the following
* conditions:
* - the line segments are [0,1]-parametrised with equally
*   distributed parametrisation (i.e. something like (1-t^2)*P1+t^2*P2
*   is not allowed as parametrisation of a line)
*
* In:
*  KBCTC   - Pointer vector for vertices on the boundary, coarse grid
*  DVBDPC  - Parameter values of boundary points, coarse grid
*  NBCT    - Number of boundary comonents
*
* Out:
*  KBCTF   - Pointer vector for vertices on the boundary, fine grid
*  DVBDPF  - Parameter values of boundary points, fine grid
*
* This only works with quadrilaterals.
************************************************************************

      SUBROUTINE PRC2FP (KBCTC,DVBDPC,NBCT,
     *                   KBCTF,DVBDPF)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C parameters

      INTEGER NBCT,KBCTC(*),KBCTF(*)
      DOUBLE PRECISION DVBDPC(*),DVBDPF(*)
      
C externals

      DOUBLE PRECISION TMAX
      EXTERNAL TMAX
      
C local variables

      INTEGER IBCT,IVBD,NVBC,NVBF,IFFST,ICFST
      
C loop through the bondary components

      DO IBCT=1,NBCT
      
C loop through the vertices

        NVBC  = KBCTC(IBCT+1)-KBCTC(IBCT)
        NVBF  = KBCTF(IBCT+1)-KBCTF(IBCT)
        IFFST = KBCTF(IBCT)
        ICFST = KBCTC(IBCT)

        DO IVBD=0,NVBC-1
        
C Copy them, leaving a free placeholder inbetween.
C Since we are using regular refinement, every fine grid has
C exactly 2x as many boundary nodes as the coarse grid,
C and there's one new node between two coarse grid nodes.

          DVBDPF(IFFST+2*IVBD) = DVBDPC(ICFST+IVBD)
          
C Recalculate the placeholder by interpolation, except for
C the last node.

          IF (IVBD.GT.0) THEN
            DVBDPF(IFFST+2*IVBD-1) = 0.5D0*
     *        ( DVBDPF(IFFST+2*IVBD-2) + DVBDPF(IFFST+2*IVBD) )
          END IF
        
        END DO
        
C Recalculate the last placeholder. Different treatment as the other 
C points because the last point of the segment is the first of the
C boundary component

        DVBDPF(IFFST+NVBF-1) = 0.5D0*
     *        ( DVBDPF(IFFST+NVBF-2) + TMAX(IBCT) )

      END DO

      END
      
************************************************************************
* Reconstruct parameter values of boundary nodes
* 
* This routine anaylyses the coordinates of vertices on the (real)
* boundary and reconstructs their parameter values.
*
* For now this works only if the grid fulfills the following
* conditions:
* - only line segments, no circles on the boundary
* - the line segments are [0,1]-parametrised with equally
*   distributed parametrisation (i.e. something like (1-t^2)*P1+t^2*P2
*   is not allowed as parametrisation of a line)
* Parameter values that are integer numbers are not corrected; it's
* assumed that these points are fixed anchors of the geometry.
*
* In:
*  DCORVG - Coordinates of the gridpoints
*  KVBD   - Numbers of vertices on the boundary
*  KBCT   - Pointer vector for vertices on the boundary
*  DVBDP  - Parameter values of boundary points
*  NBCT   - Number of boundary comonents
*
* Out:
*  DVBDP  - Parameter values of boundary points, corrected
************************************************************************

      SUBROUTINE RXBNPR (DCORVG,KVBD,KBCT,DVBDP,NBCT)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'

C parameters
      
      INTEGER NBCT,KVBD(*),KBCT(*)
      DOUBLE PRECISION DCORVG(2,*),DVBDP(*)
      
C externals

      DOUBLE PRECISION TMAX
      EXTERNAL TMAX
      
C local variables

      INTEGER IBCT,ISEGM,NSEGM,NVBS,ISI,IEI,I,ICS,ICE,ISV,IEV,IFIRST
      DOUBLE PRECISION DSSTRT,DSEND,SLEN,SPLEN
      
C loop through the bondary components

      DO IBCT=1,NBCT
      
C In each boundary component loop through the boundary segments. One
C segment can contain arbitrary many intermediate segments of the
C refinement, so we first have to find the starting- and ending point
C of the segment. They are marked by the property that:
C - the first point contains parameter value 0D0
C - each corner point of a segment has an integer parameter value
C - the last point is equal to the first point
C The length of the segment we can obtain by using TMAX. Although this
C function returns a double, the returned number is indeed an integer
C because of the [0,1]-parametrization.

C We add a slight correction to TMAX to avoid probable rounding errors
C with AINT...

        NSEGM = AINT(TMAX(IBCT)+0.000001)

C Get the number of vertices on the current boundary segment

        NVBS = KBCT(IBCT+1)-KBCT(IBCT)

C Index of the first node on this boundary:

        IFIRST = KBCT(IBCT)

C ISI and IEI contain the index of the starting- and ending-point
C of the current segment. Initialise them with values that point
C "before" the current boundary component.

        ISI = KBCT(IBCT)-1
        IEI = KBCT(IBCT)
        
C DSSTRT and DSEND contain the parameter value of the starting-
C end ending point of the current segment

        DSSTRT = 0D0
        DSEND  = 0D0

C Loop through the segments

        DO ISEGM=1,NSEGM
        
C Find the starting- and ending point; start the search after the
C last found segment.
C Our last ending point is the new starting point
        
          ISI = IEI
          DSSTRT = DSEND

C Don't do anything if we are on the last segment and
C there's no node on that last boundary segment

          IF (ISI.LT.(NVBS+IFIRST-1)) THEN

C Search the endpoint of the segment

            DO IEI=ISI+1,NVBS+IFIRST-1
              IF (DVBDP(IEI).EQ.AINT(DVBDP(IEI))) GOTO 10
            END DO
10          CONTINUE

C Obtain the parameter value of the ending point.
C Obtain the starting- and ending-index of the points in the "inner"
C of the current boundary segment. Obtain the length of the
C current segment.

            IF (IEI.GT.(NVBS+IFIRST-1)) THEN

C We have left our allowed parameter interval since we have "behind"
C the "last point". The ending point is the starting point of the interval
C since the curve is closed. It has got exactly the maximum
C parameter value!

              DSEND = TMAX(IBCT)
            
C Vertex number of starting- and ending points

              ISV = KVBD (ISI)
              IEV = KVBD (IFIRST)

            ELSE
          
              DSEND = DVBDP(IEI)
            
C Vertex number of starting- and ending points

              ISV = KVBD (ISI)
              IEV = KVBD (IEI)

            END IF

C Starting- and ending index of "inner" nodes

            ICS = ISI+1
            ICE = IEI-1
          
C Length of the line segment:

            SLEN = DSQRT((DCORVG(1,IEV)-DCORVG(1,ISV))**2+
     *                  (DCORVG(2,IEV)-DCORVG(2,ISV))**2)

            IF (SLEN.EQ.0D0) THEN
              WRITE (*,*) 'Error in RXMPTS: SLEN=0 !'
              STOP
            END IF

C Loop through all "inner" points of the segment (ICS..ICE) and
C correct their parameter values. This is done by analysing the
C distance from the point to the starting- and ending-point of
C the line.

            DO I=ICS,ICE
          
C Length of the part of the segment between current node and
C starting node of the segment

              SPLEN = DSQRT((DCORVG(1,KVBD(I))-DCORVG(1,ISV))**2+
     *                      (DCORVG(2,KVBD(I))-DCORVG(2,ISV))**2)
     
C Use this to rebuild the parameter value
     
              DVBDP (I) = DSSTRT + SPLEN / SLEN
            
            END DO

          END IF

        END DO
      
      END DO
      
      END
      
      
************************************************************************
* Reconstruct parameter values of boundary midpoints
* 
* This routine anaylyses the coordinates of midpoits on the (real)
* boundary and reconstructs their parameter values.
* This is done by simple interpolation between the neighboring
* boundary nodes.
*
* For now this works only if the grid fulfills the following
* conditions:
* - the line segments are [0,1]-parametrised with equally
*   distributed parametrisation (i.e. something like (1-t^2)*P1+t^2*P2
*   is not allowed as parametrisation of a line)
*
* In:
*  KVBD   - Numbers of vertices on the boundary
*  KBCT   - Pointer vector for vertices on the boundary
*  DVBDP  - Parameter values of boundary points
*  DMBDP  - Parameter values of boundary midpoints
*  NBCT   - Number of boundary comonents
*
* Out:
*  DMBDP  - Parameter values of boundary midpoints, corrected
************************************************************************

      SUBROUTINE RXBMPR (KBCT,DVBDP,DMBDP,NBCT)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'

C parameters
      
      INTEGER NBCT,KBCT(*)
      DOUBLE PRECISION DVBDP(*),DMBDP(*)
      
C externals

      DOUBLE PRECISION TMAX
      EXTERNAL TMAX
      
C local variables

      INTEGER IBCT,IVBD
      
C loop through the bondary components

      DO IBCT=1,NBCT
      
C loop through the vertices

        DO IVBD=KBCT(IBCT),KBCT(IBCT+1)-2
        
C Rebuild the parameter of the midpoint

          DMBDP(IVBD)=0.5D0*(DVBDP(IVBD)+DVBDP(IVBD+1))
        
        END DO
        
C The "last" boundary midpoint has to be build manually because
C the endpoint of the line is the starting point of the segment

        DMBDP(KBCT(IBCT+1)-1) = 0.5D0*(DVBDP(KBCT(IBCT+1)-1)+TMAX(IBCT))
        
      END DO
      
      END
      
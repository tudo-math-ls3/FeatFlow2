************************************************************************
* This file contains routines to calculate A Posteriori errors.
************************************************************************

************************************************************************
* Calculate A Posteriori H1-error
*
* This routine calculates an approximation to the H1-error using the
* solution of a calculation. 
*
* In:
*   DU     - array [1..NVT] of double
*            Solution vector in the vertices of the grid
*   NVT    - Number of vertices in the grid
*   TRIA   - Current triangulation
*
* Out:
*   ERRH1  - array [1..NVT] of double
*            An approximation to the H1-error in every vertex of the
*            triagulation
************************************************************************

      SUBROUTINE ERPCH1(DU,NVT,TRIA,ERRH1)

      IMPLICIT NONE

      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'

      INTEGER NVT
      DOUBLE PRECISION DU(NVT)
      INTEGER TRIA(SZTRIA)
      DOUBLE PRECISION ERRH1(*)

C     local variables

      INTEGER LGR1X,LGR1Y,LGR2X,LGR2Y,I
      
      EXTERNAL E011

C     Allocate auxiliary vectors for gradient information

      CALL ZNEW (NVT,1,LGR1X,'DGR1X ')
      CALL ZNEW (NVT,1,LGR1Y,'DGR1Y ')
      CALL ZNEW (NVT,1,LGR2X,'DGR2X ')
      CALL ZNEW (NVT,1,LGR2Y,'DGR2Y ')
      
C     The H1-error is defined cellwise as
C
C         | grad(u) - grad(u_h) |_T
C
C     Ok, we don't have the u... instead, we replace grad(u)
C     with grad_2(u_h) with grad_2 being a higher order
C     gradient reconstruction method -- the ZZ-technique.
C     grad(u) we calculate as usual by the first order method.
C
C     Calculate the gradient with standard gradient evaluation

      CALL XGRDR3(TRIA,DU,DWORK(L(LGR1X)),DWORK(L(LGR1Y)),NVT,
     *            E011,.FALSE.)

C     Calculate the gradient with higher-order ZZ-technique

      CALL XGRDR3(TRIA,DU,DWORK(L(LGR2X)),DWORK(L(LGR2Y)),NVT,
     *            E011,.FALSE.)

C     Subtract the solutions

      CALL LLC1 (DWORK(L(LGR1X)),DWORK(L(LGR2X)),NVT,1D0,-1D0)
      CALL LLC1 (DWORK(L(LGR1Y)),DWORK(L(LGR2Y)),NVT,1D0,-1D0)

C     Build the l2-norm in every vertex; this is the approximation
C     to the error

      DO I=0,NVT-1
        ERRH1 (1+I) = SQRT(DWORK(L(LGR1X)+I)**2 + DWORK(L(LGR1Y)+I)**2)
      END DO

C     Release auxiliary vectors

      CALL ZDISP (0,LGR2Y,'DGR2Y ')
      CALL ZDISP (0,LGR2X,'DGR2X ')
      CALL ZDISP (0,LGR1Y,'DGR1Y ')
      CALL ZDISP (0,LGR1X,'DGR1X ')

      END
      
************************************************************************
* Calculate A Posteriori H1-error from velocity vector
*
* This routine calculates an approximation to the H1-error using the
* solution of a calculation. The solution is assumed to be a 2D
* velocity vector, discretized with Q1~. 
*
* In:
*   DU1    - array [1..NU] of double
*            X-velocity
*   DU2    - array [1..NU] of double
*            Y-velocity
*   NU     - Number of velocity components in each velocity
*   TRIA   - Current triangulation
*
* Out:
*   ERRH1  - array [1..NVT] of double
*            An approximation to the H1-error in every vertex of the
*            triagulation
************************************************************************

      SUBROUTINE ERPVH1(DU1,DU2,NU,TRIA,NVT,ERRH1,
     *                  IASMBL,DASMBL,IGEOM,DGEOM)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'

      INCLUDE 'stria.inc'
      
      INCLUDE 'sassembly.inc'
      
      INTEGER NU,NVT
      DOUBLE PRECISION DU1(NU),DU2(NU)
      INTEGER TRIA(SZTRIA)
      DOUBLE PRECISION ERRH1(NVT)
      INTEGER IASMBL(*),IGEOM(*)
      DOUBLE PRECISION DASMBL(*),DGEOM(*)

C     local variables

      INTEGER LH1,LH2,LAUX,I
      INTEGER KCORVG,KXNPR
      DOUBLE PRECISION RE
      EXTERNAL UE
      
C     Allocate auxiliary vectors to calculate the solution in the
C     vertices

      CALL ZNEW (NVT,1,LAUX,'DAUX  ')
      
C     Interpolate the solution to the vertices.
C     This returns handles LH1, LH2.
      
      LH1=0
      LH2=0
      CALL XINTUV (DU1,DU2,TRIA,LH1,LH2)
      
C     Implement boundary conditions, since the interpolation
C     does not correctly handle boundaries:
      
      KCORVG = L(TRIA(OLCORVG))
      KXNPR  = L(TRIA(OLXNPR))
      RE = DASMBL(ORE)
      CALL BDRCOR (DWORK(L(LH1)),DWORK(L(LH2)),
     *             TRIA,DWORK(KCORVG),
     *             KWORK(KXNPR),UE,1D0,RE,
     *             IASMBL,DASMBL,IGEOM,DGEOM)
     
C     Calculate the error of each component

      CALL ERPCH1(DWORK(L(LH1)),NVT,TRIA,ERRH1)
      CALL LCP1 (ERRH1,DWORK(L(LH1)),NVT)
      CALL ERPCH1(DWORK(L(LH2)),NVT,TRIA,ERRH1)
      
C     And add them together using a norm

      DO I=1,NVT
        ERRH1(I) = SQRT(ERRH1(I)**2+DWORK(L(LH1)+I-1)**2)
      END DO
      
C     Release memory

      CALL ZDISP (0,LAUX,'DAUX  ')
      CALL ZDISP (0,LH2,'DH2   ')
      CALL ZDISP (0,LH1,'DH1   ')

      END
            
************************************************************************
* Calculate A Posteriori H1-error from velocity vector
*
* This routine calculates an approximation to the H1-error using the
* solution of a calculation. The solution is assumed to be a 2D
* velocity vector, discretized with Q1. 
*
* In:
*   DU1    - array [1..NVT] of double
*            X-velocity
*   DU2    - array [1..NVT] of double
*            Y-velocity
*   NVT    - Number of velocity components in each velocity
*            = number of vertices in the triangulation
*   TRIA   - Current triangulation
*
* Out:
*   ERRH1  - array [1..NVT] of double
*            An approximation to the H1-error in every vertex of the
*            triagulation
************************************************************************

      SUBROUTINE ERPQH1(DU1,DU2,NVT,TRIA,ERRH1)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'

      INCLUDE 'stria.inc'
      
      INCLUDE 'sassembly.inc'
      
      INTEGER NVT
      DOUBLE PRECISION DU1(NVT),DU2(NVT)
      INTEGER TRIA(SZTRIA)
      DOUBLE PRECISION ERRH1(NVT)

C     local variables

      INTEGER LAUX,I
      
C     Allocate auxiliary vectors to calculate the solution in the
C     vertices

      CALL ZNEW (NVT,1,LAUX,'DAUX  ')
      
C     Calculate the error of each component

      CALL ERPCH1(DU1,NVT,TRIA,ERRH1)
      CALL ERPCH1(DU2,NVT,TRIA,DWORK(L(LAUX)))
      
C     And add them together using a norm

      DO I=1,NVT
        ERRH1(I) = SQRT(ERRH1(I)**2+DWORK(L(LAUX)+I-1)**2)
      END DO
      
C     Release memory

      CALL ZDISP (0,LAUX,'DAUX  ')

      END
                        
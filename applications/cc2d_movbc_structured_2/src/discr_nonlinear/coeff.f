************************************************************************
      DOUBLE PRECISION FUNCTION   COEFST (X,Y,IA,IB,IBLOC,BFIRST,
     *                                    TRIA,IPARAM,DPARAM)
*
*     Coefficient for the Stokes-block.
*     IPARAM/DPARAM points to the assembly structure
************************************************************************
      
      IMPLICIT NONE

C main COMMON blocks
      
      INCLUDE 'sassembly.inc'

C parameters

      INTEGER IA,IB,IBLOC
      DOUBLE PRECISION X,Y
      LOGICAL BFIRST
      
      INTEGER IPARAM(*),TRIA(*)
      DOUBLE PRECISION DPARAM(*)

C local variables

      IF ((IA.EQ.1).AND.(IB.EQ.1)) THEN
       COEFST=1D0
      ELSE
       COEFST=DPARAM(ONY)
      ENDIF
C
      END

************************************************************************
      DOUBLE PRECISION FUNCTION   COEFFM (X,Y,IA,IB,IBLOC,BFIRST,
     *                                    TRIA,IPARAM,DPARAM)
*
*     Coefficient for the Mass matrix
************************************************************************
      
      IMPLICIT NONE

C parameters

      INTEGER IA,IB,IBLOC
      DOUBLE PRECISION X,Y
      LOGICAL BFIRST
      
      INTEGER IPARAM(*),TRIA(*)
      DOUBLE PRECISION DPARAM(*)

C local variables

      IF ((IA.EQ.1).AND.(IB.EQ.1)) THEN
       COEFFM=1D0
      ENDIF
C
      END

************************************************************************
      DOUBLE PRECISION FUNCTION COEFFN(IA,IB,U1L1,U1L2,U2L1,U2L2,
     *                                 A1L,A2L,
     *                                 IDFL, KDFL, KDFG,
     *                                 DELTA,DCMASS,THSTEP,IPRECA,NY,
     *                                 IPARAM,DPARAM)
*
*     Coefficient for the convective-block.
*     Not used in current implementation!
************************************************************************
      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'

C parameters

      DOUBLE PRECISION U1L1(*),U1L2(*),U2L1(*),U2L2(*)
      DOUBLE PRECISION A1L,A2L,DELTA,DCMASS,THSTEP,NY
      INTEGER IA, IB, IPRECA, IDFL, KDFL(*), KDFG(*)
      
      INTEGER IPARAM(*)
      DOUBLE PRECISION DPARAM(*)

C local variables

      DOUBLE PRECISION DU1, DU2, HBAS
      INTEGER JDFL

      DU1=0D0
      DU2=0D0
      
      THSTEP = DPARAM(1)
      IPRECA = IPARAM(1)

      IF (DCMASS.NE.0D0) THEN
C
       DO 10 JDFL=1,IDFL
       HBAS=DBAS(KDFL(JDFL),1)
       DU1=DU1+(A1L*U1L1(KDFG(JDFL))+A2L*U2L1(KDFG(JDFL)))*HBAS
       DU2=DU2+(A1L*U1L2(KDFG(JDFL))+A2L*U2L2(KDFG(JDFL)))*HBAS
10     CONTINUE
C
       IF ((IA.EQ.1).AND.(IB.EQ.2)) COEFFN=DU1
       IF ((IA.EQ.1).AND.(IB.EQ.3)) COEFFN=DU2
       IF ((IA.EQ.2).AND.(IB.EQ.3)) COEFFN=DELTA*DU1*DU2
       IF ((IA.EQ.3).AND.(IB.EQ.2)) COEFFN=DELTA*DU1*DU2
C
       IF ((IA.EQ.1).AND.(IB.EQ.1)) COEFFN=DCMASS/THSTEP
C
       IF (IPRECA.EQ.4) THEN      
        IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFFN=DELTA*DU1**2+NY
        IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFFN=DELTA*DU2**2+NY
       ELSE
        IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFFN=DELTA*DU1**2
        IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFFN=DELTA*DU2**2
       ENDIF
C
      ELSE
C
       COEFFN=0D0
       IF ((IA.EQ.1).AND.(IB.EQ.1)) COEFFN=-1D0/THSTEP
C
      ENDIF
C
C
      END
C
************************************************************************
      DOUBLE PRECISION FUNCTION   COEFFB (X,Y,IA,IB,IBLOC,BFIRST,
     *                                    TRIA,IPARAM,DPARAM)
*
*     Coefficient for the B1/B2-blocks
************************************************************************
      
      IMPLICIT NONE

C parameters

      DOUBLE PRECISION X,Y
      INTEGER IA, IB, IBLOC
      LOGICAL BFIRST
      
      INTEGER IPARAM(*),TRIA(*)
      DOUBLE PRECISION DPARAM(*)

      COEFFB= -1.D0

      END

************************************************************************
* Exact solution - velocity, also for boundary conditions
*
* This routine determines - depending on a point - the exact
* solution value in that point.
*
* In:
*  X,Y    - coordinates of the point
*  IBLOC  - matrix block that specifies the type of the solution to
*           return.
*           = 0: matrix block for x-velocity
*           = 1: matrix block for y-velocity
*  TIMENS - Current time in instationary Navier-Stokes calculation
*           =0 for stationary simulation.
*  RE     - Reynolds number; from the DAT file
*  IPARAM - array [1..*] of integer 
*  DPARAM - array [1..*] of integer 
*           TIntAssembly/TDoubleAssembly assembly structures; gives
*           additional information about the discretization.
*           This is passed to user defined callback routines so that 
*           they can access more detailed information about the
*           problem. Not used in this routine.
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to fictitious boundary
*           routines. Not used in this routine.
*
* "Optional" parameters:
*  INODE  - Number of the node corresponding to (X,Y). Can be 0.
*           If INODE=0, this function must create the necessary 
*              information only by the given (X,Y)-coordinates of 
*              the point
*           If INODE>0, the parameter TRIA must also be defined.
*              In this case INODE is a node number in the triangulation
*              TRIA with coordinates (X,Y). Then this routine can also
*              access additional information about the triangulation
*              (like precomputed values) to compute the result.
*  TRIA   - array [1..SZTRIA] of integer
*           If INODE>0, TRIA defines the triangulation INODE refers to.
*           If INODE=0, this parameter can be a dummy parameter.
************************************************************************

      DOUBLE PRECISION FUNCTION UE (X,Y,IBLOC,TIMENS,RE,
     *                          INODE,TRIA,IPARAM,DPARAM,IGEOM,DGEOM)
      
      IMPLICIT NONE

C parameters

      DOUBLE PRECISION X,Y
      DOUBLE PRECISION DPARAM(*),TIMENS,RE,DGEOM(*)
      INTEGER IBLOC
      LOGICAL BFIRST
      INTEGER INODE,TRIA(*),IPARAM(*),IGEOM(*)
      
C externals
      DOUBLE PRECISION FDATIN
      EXTERNAL FDATIN

      UE=FDATIN(1,IBLOC,X,Y,TIMENS,RE,
     *          INODE,TRIA,IPARAM,DPARAM,IGEOM,DGEOM)
C
      END
C
************************************************************************
*     Exact solution - pressure, only for error analysis
************************************************************************

      DOUBLE PRECISION FUNCTION PE (X,Y, IBLOC,TIMENS,RE,
     *                              INODE,TRIA,IPARAM,DPARAM,
     *                              IGEOM,DGEOM)
      
      IMPLICIT NONE

C parameters

      DOUBLE PRECISION X,Y
      INTEGER IBLOC
      DOUBLE PRECISION DPARAM(*),TIMENS,RE,DGEOM(*)
      INTEGER INODE,TRIA(*),IPARAM(*),IGEOM(*)

C externals
      DOUBLE PRECISION FDATIN
      EXTERNAL FDATIN
      
C local variables    

      PE=FDATIN(4,IBLOC,X,Y,TIMENS,RE,
     *          INODE,TRIA,IPARAM,DPARAM,IGEOM,DGEOM)
C
      END
C
*************************************************************************
      DOUBLE PRECISION FUNCTION UEX(X,Y, IBLOC,TIMENS,RE,
     *                              INODE,TRIA,IPARAM,DPARAM,
     *                              IGEOM,DGEOM)
C
C     x-derivative of exact solution, only for error analysis
*************************************************************************
      
      IMPLICIT NONE

C parameters

      DOUBLE PRECISION X,Y
      INTEGER IBLOC
      DOUBLE PRECISION DPARAM(*),TIMENS,RE,DGEOM(*)
      INTEGER INODE,TRIA(*),IPARAM(*),IGEOM(*)

C externals
      DOUBLE PRECISION FDATIN
      EXTERNAL FDATIN
      
      UEX=FDATIN(2,IBLOC,X,Y,TIMENS,RE,
     *           0,TRIA,IPARAM,DPARAM,IGEOM,DGEOM)
C
      END
C
*************************************************************************
      DOUBLE PRECISION FUNCTION UEY(X,Y, IBLOC,TIMENS,RE,
     *                              INODE,TRIA,IPARAM,DPARAM,
     *                              IGEOM,DGEOM)
C
C     y-derivative of exact solution,, only for error analysis
*************************************************************************
      
      IMPLICIT NONE

C parameters

      DOUBLE PRECISION X,Y
      INTEGER IBLOC
      DOUBLE PRECISION DPARAM(*),TIMENS,RE,DGEOM(*)
      INTEGER INODE,TRIA(*),IPARAM(*),IGEOM(*)

C externals
      DOUBLE PRECISION FDATIN
      EXTERNAL FDATIN
      
      UEY=FDATIN(3,IBLOC,X,Y,TIMENS,RE,
     *           0,TRIA,IPARAM,DPARAM,IGEOM,DGEOM)   
C
      END

************************************************************************
* This file contains initialization routines for the NSDEF2 solver.
************************************************************************

************************************************************************
* Initialize vector structures for Navier-Stokes solver
*
* This routine allocates memory for necessary vectors and
* initializes the VECDAT structures on all levels, which are used 
* for solving the Navier Stokes equation in all cases.
*
* In:
*   NLMIN   : Minimum level in TRIAS
*   NLMAX   : Maximum level in TRIAS
*   TRIAS   : array [1..SZTRIA,1..NNLEV] of integer
*             Triangulation structures on all levels.
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. This tells all assembly-routines how to 
*            set up the nonlinearity in the system matrices, which 
*            cubature formula to use, etc.
*
* Out:
*   VECDAT  : array [1..SZN2VI,1..NNLEV] of integer
*             TNS2DVectorParams structure that defines the vector
*             structure. VECDAT(.,NLMIN..NLMAX) is filled with data.
************************************************************************

      SUBROUTINE ININSV (NLMIN,NLMAX,TRIAS,IASMBL,DASMBL,
     *                   VECDAT)
      
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'sassembly.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
      INCLUDE 'cbasicmg.inc'
      
C     parameters      
      
      INTEGER NLMIN,NLMAX,TRIAS(SZTRIA,NNLEV),IASMBL(*)
      INTEGER VECDAT(SZN2VI,NNLEV)
      DOUBLE PRECISION DASMBL(*)
      
C     Functions to get the number of DOF's

      INTEGER NDFGX
      EXTERNAL NDFGX
            
C     local variables
      
      INTEGER I,NUP

C     Loop through all essential levels:      

      DO I=NLMIN,NLMAX

C       At first clear the structure for level I

        CALL LCL3(VECDAT(1,I),SZN2VI)
        
C       Get the number of unknowns for velocity and pressure.
C       This is done with the central element routine which calculates
C       the global number of DOF's.
        
        IF ((IASMBL(OIELEMT).EQ.0).OR.(IASMBL(OIELEMT).EQ.2)) THEN
          VECDAT(ONU,I) = NDFGX(31,TRIAS(1,I))
        ELSE IF ((IASMBL(OIELEMT).EQ.1).OR.(IASMBL(OIELEMT).EQ.3)) THEN
          VECDAT(ONU,I) = NDFGX(30,TRIAS(1,I))
        END IF
        VECDAT(ONP,I) = NDFGX(10,TRIAS(1,I))
        
C       Summing these values results in the number of unknowns:
        
        NUP = 2*VECDAT(ONU,I)+VECDAT(ONP,I)
        VECDAT(ONEQV,I) = NUP

C       Allocate solution vector

        CALL ZNEW(NUP,1,VECDAT(OLSOL,I),'DSOL  ')

C       Allocate RHS vector

        CALL ZNEW(NUP,1,VECDAT(OLRHS,I),'DRHS  ')
        
C       Allocate auxiliary vector

        CALL ZNEW(NUP,1,VECDAT(OLTMP,I),'DAUX  ')
        
      END DO
      
      END
      
************************************************************************
* Release vector structures for Navier-Stokes solver
*
* This routine releases memory that was allocated by ININSV.
*
* In:
*   NLMIN   : Minimum level in TRIAS
*   NLMAX   : Maximum level in TRIAS
*   VECDAT  : array [1..SZN2VI,1..NNLEV] of integer
*             TNS2DVectorParams structure that defines the vector
*             structure. 
*
* Out:
*  VECDAT(.,NLMIN..NLMAX) is released.
************************************************************************

      SUBROUTINE DONNSV (NLMIN,NLMAX,VECDAT)
      
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'sassembly.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
      INCLUDE 'cbasicmg.inc'
      
C     parameters      
      
      INTEGER NLMIN,NLMAX,TRIAS(SZTRIA,NNLEV)
      INTEGER VECDAT(SZN2VI,NNLEV)
      
C     local variables
      
      INTEGER I

C     Loop through all essential levels:      

      DO I=NLMAX,NLMIN,-1

C       Release auxiliary vector

        CALL ZDISP(0,VECDAT(OLTMP,I),'DAUX  ')
        
C       Release RHS vector

        CALL ZDISP(0,VECDAT(OLRHS,I),'DRHS  ')
        
C       Release solution vector

        CALL ZDISP(0,VECDAT(OLSOL,I),'DSOL  ')

      END DO
      
      END

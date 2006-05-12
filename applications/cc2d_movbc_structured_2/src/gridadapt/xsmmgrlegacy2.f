************************************************************************
* Static grid adaption with correction of multiple levels
*
* Featflow Legacy call
*
* This subroutine performs grid-adaption according to information
* in the COMMON-blocks, using the multigrid COMMON-block structures.
*
* In:
*   IGALEV - Level where the grid adaption should be performed.
*             >0: Perform grid adaption on level IGALEV
*            <=0: Perform grid adaption on level NLMAX+IGALEV
************************************************************************
      
      SUBROUTINE XSMMGL (IGALEV)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
      INCLUDE 'ctria.inc'

C Variables for the static method:
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'sgridadaptmgstaticparams.inc'
      INCLUDE 'cgridadapt.inc'

      INTEGER TRIAS(SZTRIA,NNLEV)
      INTEGER IGALEV

C local variables

      INTEGER II
      DOUBLE PRECISION TGES
      DOUBLE PRECISION TGTOT,TGLSPR,TGLS,TGODE,TGGPR,TGGSM,TGMSM
      DOUBLE PRECISION TGGAPR,TGGRC
      INTEGER NLSYS ,NLSITE,NCODE ,NMFEVL, ILVOLD
      
C     Save the current level
      
      ILVOLD = ILEV
      
C     Clean the TRIAS array - necessary initialization, as the
C     TRIAS is on the stack...

      CALL LCL3 (TRIAS,SZTRIA*NNLEV)
      
C     Build up the TRIAS-array. Copy all level information
C     between level NLMIN and NLMAX to the structure
C     array TRIAS. Start with level 1, not MILV, even if
C     it does perhaps not exist - to make sure to copy
C     as many levels as possible, as they are needed for the
C     multigrid solver in the grid adaption.
    
      II=1
      DO ILEV=NLMIN,NLMAX
      
        CALL SETLEV (II)
        CALL C2TRIA(TRIAS(1,ILEV))

C       Create extended structures

        CALL GENETR (TRIAS(1,ILEV))

      END DO
      
C     Calculate the level of the grid adaption

      II = MIN(IGALEV,NLMAX)
      IF (II.LE.0) II = NLMAX+II
      II = MAX(NLMIN,II)
     
C     Call the grid adaption

      CALL XSMMGW (TRIAS,NLMIN,II,0,0D0,0,0D0)
      
C     Release extended triangulation structures,
C     Recalculate DCORMG if necessary.

      DO ILEV=NLMAX,NLMIN,-1
      
        CALL DISETR (TRIAS(1,ILEV),0)

      END DO
      
C     Restore the old level

      ILEV = MIN(NLMAX,ILVOLD)
      II=1
      CALL SETLEV (II)
      ILEV = ILVOLD

      END 

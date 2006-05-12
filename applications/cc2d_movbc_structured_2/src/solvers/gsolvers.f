************************************************************************
* This file contains maintainance-routines for the general solver
* structures TSolverIParams/TSolverDParams. These routines
* can be used to initialize the structure with standard values.
************************************************************************
      
************************************************************************
* Solver structure initialization
*
* The following routine can be used to initialise the IPARAM/DPARAM
* array structures with default values for the computation. This
* is typically done in a solver-specific initialization routine,
* which initializes solver-specific variables right after 
* this routine initialized the general parameters.
*
* In:
*   -
* Out:
*   IPARAM  - array [1..SZSLVI] of integer
*             Integer parameter structure
*   DPARAM  - array [1..SZSLVI] of integer
*             Double precision parameter structure
************************************************************************

      SUBROUTINE INGSLV (IPARAM, DPARAM)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'

      INTEGER IPARAM(SZSLVI)
      DOUBLE PRECISION DPARAM(SZSLVD)
      
C Clear the structures

      CALL LCL3 (IPARAM,SZSLVI)
      CALL LCL1 (DPARAM,SZSLVD)
      
C Set the non-zero standard values
      
      IPARAM (ONITMIN) = 0
      IPARAM (ONITMAX) = 50
      IPARAM (OTIM  )  = 2
      IPARAM (OIASRLN) = 3
      IPARAM (OMSGTRM) = 2
      IPARAM (OMTRMRS) = 1

      DPARAM (OEPSREL) = 1D-5
      DPARAM (OEPSABS) = 0D0
      DPARAM (ODIVREL) = 1D3
      DPARAM (ODIVABS) = 1D99
      DPARAM (OVECZER) = 1D-12
      
      END

************************************************************************
* General timing auxiliary routine for solvers
*
* This routine can by used by all types of solver to generate timing
* statistics and collect them in a double-array.
* It allowes direct computation of timing information by only
* one call. It accepts two DPARAM blocks:
*  TPARAM - a temporary array storing start-times
*  DPARAM - a destination array for storing final time statistics
* Both parameter blocks must be of the same size.
* If MDE=0, the current system time is saved in TPARAM(OFS). If
* MDE=1, the difference between TPARAM(OFS) and the current system
* time is added to DPARAM(OFS).
*
* In:
*   DPARAM : array [1..*] of double 
*            Destination array for timing statistics
*   TPARAM : array [1..*] of double 
*            Temporary array for storing the start-time
*   OFS    : the offset in DPARAM/TPARAM to the current timing
*            information that should be modified
*   MDE    : the mode.
*            =0: save current time in TPARAM(OFS)
*            =1: add time difference to TPARAM(OFS) in DPARAM(OFS)
*
* Out:
*  TPARAM(OFS) or DPARAM(OFS), depending on MDE
************************************************************************

      SUBROUTINE GTMAUX (TPARAM,DPARAM,OFS,MDE)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'

      DOUBLE PRECISION TIM,DPARAM(*),TPARAM(*)
      INTEGER MDE,OFS

      IF (MDE.EQ.0) THEN
        CALL ZTIME(TPARAM(OFS))
      ELSE
        CALL ZTIME(TIM)
        DPARAM(OFS) = DPARAM(OFS) + TIM - TPARAM(OFS)
      END IF

      END
      

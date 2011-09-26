************************************************************************
* This file describes a structure for collecting statistical information
* about time time that is needed by the different parts of the
* nonlinear solver.
*
* The structure itselves is realized as integer/double precision
* arrays:
*   STIMIN = array [1..SZTIMG] of integer 
*          = TTimings
* where the different entries can be accessed with the offsets defined
* here, e.g.
*   TTUPW = STIMIN(OTTUPW)
*
* The content of the structure is as follows:
*
* TTimings = record
*  
*   // The first part of the structure coincides with the timing
*   // result part of the MG solver structure:
*  
*   TTMG   : Total time for MG solver
*   TTPROL : Time for prolongation in MG
*   TTREST : Time for restriction in MG
*   TTDEF  : Time for defect calculation in MG
*   TTSMTH : Time for smoothing operations in MG
*   TTCGSL : Time for coarse-grid solver in MG
*   TTBC   : Time for boundary conditions in MG
*   TTCGC  : Time for coarse-grid corrections in MG
*
*   // The second part describes the general timing information
*   // for the nonlinear solver:
*
*   TTNL   : Total time for solving the nonlinear equation
*   TTLSOL : Total time for solving of linear equations
*            (as subproblems of the nonlinear equation),
*            independent of the linear solver
*   TTUPW  : Total time for setting up nonlinearity (Upwind/Streamline 
*            diff.) and inclusion of stabilisation terms
*   TTADF  : Total time for setting up linear part of nonlinear system
*            matrices
*   TTBDR  : Total rime for implementing boundary conditions into matrix/
*            vectors
*   TTCOR  : Total time for calculation optimal correction
*   TTPOST : Total time for postprocessing
*   TTLC   : Total time for copying of vectors and restriction of
*            vectors
* end
************************************************************************

      INTEGER OTTMG,OTTPROL,OTTREST,OTTDEF,OTTSMTH,OTTCGSL,OTTBC,OTTCGC 

      INTEGER OTTNL,OTTUPW ,OTTADF ,OTTBDR ,OTTCOR ,OTTPOST,OTTLC
      INTEGER OTTLSOL

      PARAMETER (OTTMG   =  1)
      PARAMETER (OTTPROL =  2)
      PARAMETER (OTTREST =  3)
      PARAMETER (OTTDEF  =  4)
      PARAMETER (OTTSMTH =  5)
      PARAMETER (OTTCGSL =  6)
      PARAMETER (OTTBC   =  7)
      PARAMETER (OTTCGC  =  8)

      PARAMETER (OTTNL   =  9)
      PARAMETER (OTTLSOL = 10)
      PARAMETER (OTTUPW  = 11)
      PARAMETER (OTTADF  = 12) 
      PARAMETER (OTTBDR  = 13)
      PARAMETER (OTTCOR  = 14)
      PARAMETER (OTTPOST = 15)
      PARAMETER (OTTLC   = 16)
      
      INTEGER SZTIMG
      
      PARAMETER (SZTIMG  = 16)
      
************************************************************************
* This file defines the structure that is used as user-defined
* parameter block for the callback routines of M020 when solving
* linear subproblems for NSDEF.
*
* The structure contains all necessary information that is needed
* in the callback routines. It is set up directly before calling the
* linear solver and passed as user-defined parameter block to all
* the callback routines.
*
* The structure contains triangulation and matrix structures as sub-
* structures, therefore STRIA.INC and SMAT2DNS.INC has to be included
* before.
*
* The structure itselves is realized as integer array:
*
*   IPARAM = array [1..SZNLSI] of integer 
*          = TNSDEFILinearSolution
*
* where the different entries can be accessed with the offsets defined
* here, e.g.
*
*   ICNEUM = IPARAM(OICNEUM)
*
* The structure is defined as follows:
*
* TNSDEFILinearSolution = record
*   ICNEUM : During the solution process, this flag receives
*            whether in the current iteration there are Neumann
*            boundary components in the geometry.
*            Solver callback-routines can access this information
*            to adapt to the problem.
*   IMATS  : array [1..SZN2MI,1..NNLEV] of integer 
*            During the solution process, IMATS receives a 
*            TNS2DMatrixParams structure that defines the system
*            matrices of the problem.
*            Solver callback-routines can access this information
*            to adapt to the problem.
*   ITRIAS : array [1..SZTRIA,1..NNLEV] of integer 
*            During the multigrid iteration, ITRIAS receives
*            information about the triangulation of all levels.
*            Solver callback-routines can access this information
*            to adapt to the problem.
* end
*
************************************************************************

      INTEGER OIMATS ,OITRIAS, OICNEUM

      INTEGER SZNLSI 
      
      PARAMETER (OICNEUM = 1)
      PARAMETER (OIMATS  = 2)
      PARAMETER (OITRIAS = OIMATS + NNLEV*SZN2MI)

      PARAMETER (SZNLSI  = OITRIAS + NNLEV*SZTRIA)
      
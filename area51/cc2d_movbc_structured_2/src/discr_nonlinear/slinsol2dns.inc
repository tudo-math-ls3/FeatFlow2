************************************************************************
* This file describes structure for collecting all relevant
* information for solving a linear system in a 2D environment.
************************************************************************
* Each solution vector consists of two velocity and one pressure
* component, attached after each other.
*
* The structure is an extension to the TLinearSolverParams structure
* and adds further variables that allow to access the different
* components of the solution separately. Therefore the file
* SLINSOL.INC must me included prior to this file!
*
* The structure itself is realized as integer array:
*   IPARAM = array [1..SZN2VI] of integer
*          = TNS2DVectorParams
*
* TNS2DVectorParams = record (TLinearSolverParams)
*
*   // The first part is identical to the TLinearSolverParams structure,
*   // as it contains general information about the solution vectors:
*
*   NEQ    : Number of equations in solution vector
*   ...
*   LTMP   : Handle to array [1..NEQ] of double precision.
*
*   // The second part of this structure adds information about the
*   // structure of the solution vectors. Each solution vector
*   // consists of two velocity components of length NU and one
*   // pressure component of length NP:
*
*   NU     : integer
*   NP     : integer
* end
************************************************************************

      INTEGER ONU, ONP, SZN2VI
      
      PARAMETER (ONU    = SZLSLI + 1)
      PARAMETER (ONP    = SZLSLI + 2)
      
      PARAMETER (SZN2VI = SZLSLI + 2)

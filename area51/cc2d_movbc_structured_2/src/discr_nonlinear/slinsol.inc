************************************************************************
* This file describes a structure for collecting all relevant
* information for solving a linear system. This structure only
* represents the very basic information like handles to solution
* vectors, RHS, NEQ,... . If more detailed information is necessary
* (like how a solution vector is to be decomposed into different
* components for velocities into different directions e.g.), the
* structure should be augmented to the specific needs.
*
* The structure itselves is realized as integer array:
*   IPARAM = array [1..SZLSLI] of integer
*          = TLinearSolverParams
*
* TLinearSolverParams = record
*   NEQV   : Number of equations in solution vector
*   LSOL   : Handle to array [1..NEQ] of double precision.
*            Represents a solution vector
*   LRHS   : Handle to array [1..NEQ] of double precision.
*            Represents a right-hand-side vector.
*   LTMP   : Handle to array [1..NEQ] of double precision.
*            Represents an auxiliary defect vector that can be
*            used arbitrarily by the solver.
* end
************************************************************************

      INTEGER ONEQV,OLSOL,OLRHS,OLTMP,SZLSLI
      
      PARAMETER (ONEQV  = 1)
      PARAMETER (OLSOL  = 2)
      PARAMETER (OLRHS  = 3)
      PARAMETER (OLTMP  = 4)
      
      PARAMETER (SZLSLI = 4)

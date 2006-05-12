***********************************************************************
* This file contains initialization routines for the problem of
* solving stationary Navier Stokes equations.
***********************************************************************

***********************************************************************
* Initialize discretisation structures for nonstationary solver
*
* This routine allocates memory for the structures used in the
* nonstationary solver. It initializes all the structures, but does
* not compute anything.
*
* In:
*   NLMIN  : minimum level where calculation is allowed
*   NLMAX  : maximum level where calculation is allowed
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for the underlying meshes
*            on all levels. TRIAS(.,ILEV) represents the triangulation
*            structure on level ILEV in its ummodified form, i.e.
*            without any boundary conditions implemented, probably
*            undeformed (if the user wants to start the grid
*            deformation in every time step from the original grids).
*
* Out:
*   MATDAT : array [1..SZN2MI,1..NNLEV] of integer
*            TNS2DMatrixParams-structures for level NLMIN..NLMAX, 
*            initialized with data
*   VECDAT : array [1..SZN2VI,1..NNLEV] of integer
*            TNS2DVectorParams-structures for level NLMIN..NLMAX. 
*
*   IPARAM : array [1..SZISDI] of integer
*   DPARAM : array [1..SZISDI] of double
*            Integer and double prec. parameter blocks that define the
*            behaviour of the nonstationary solver. 
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer and double prec. parameter block for the stationary 
*            sub-solver NSDEF2.
*   IMGPAR : array [1..SZ020I+3*SZSLVI] of integer
*   DMGPAR : array [1..SZ020D+3*SZSLVD] of double
*            Integer and double parameter blocks for the multigrid 
*            sub-solver M020. 
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. 
*   IADTS  : array [1..SZADTI] of integer
*   DADTS  : array [1..SZADTD] of double
*            Integer and double prec. parameter block that configures
*            the behaviour of the adaptive time stepping.
*
*   LUP    : Handle to array [1..*] of double
*            Iteration vector; not initialized
*   LRHS   : Handle to array [1..*] of double
*            Right hand side vector
*   NEQ    : Length of solution/RHS vector
*
* INMAIN must be called prior to this routine.
* Parametrization and triangulations must be initialized prior to
* calling this routine.
***********************************************************************
      
      SUBROUTINE INISSL (MSHOW,NLMIN,NLMAX,
     *                   TRIAS,
     *                   MATDAT,VECDAT,
     *                   IPARAM,DPARAM,
     *                   ISTPAR,DSTPAR,
     *                   IMGPAR,DMGPAR,
     *                   IASMBL,DASMBL,
     *                   IADTS,DADTS,
     *                   LUP,LRHS,NEQ)
      
      IMPLICIT NONE

      INCLUDE 'stria.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'm020.inc'

      INCLUDE 'stiming.inc'
      INCLUDE 'snsdef.inc'
      INCLUDE 'ststepping.inc'
      INCLUDE 'snonstead.inc'
      
      INCLUDE 'sassembly.inc'
      
      INCLUDE 'sadtstep.inc'

C     parameters

      INTEGER MSHOW,MFILE,NLMIN,NLMAX,TRIAS(SZTRIA,*)
      INTEGER MATDAT(SZN2MI,*),VECDAT(SZN2VI,*),IPARAM(*),ISTPAR(*)
      INTEGER IMGPAR(*),IASMBL(*),IADTS(*)
      INTEGER LUP,LRHS,LAUX,NEQ
      DOUBLE PRECISION DPARAM(*),DSTPAR(*),DMGPAR(*),DASMBL(*),DADTS(*)
      
C     Most things can be initialized with the initialization
C     routine of the stationary solver:

      CALL INSTSL (NLMIN,NLMAX,
     *             TRIAS,
     *             MATDAT,VECDAT,
     *             ISTPAR,DSTPAR,
     *             IMGPAR,DMGPAR,
     *             IASMBL,DASMBL,
     *             LUP,LRHS,NEQ)
     
C     Adaptive time stepping

      CALL INIATS (IADTS,DADTS,1)
      
C     Initialize the parameter block for the nonstationary
C     solver with the appropriate routines:

      CALL INIISD (IPARAM,DPARAM,1)
      
      END
            
***********************************************************************
* Clean up structures for nonstationary solver
*
* This routine cleans up the discretisation structures for the 
* nonstationary solver, releasing all allocated memory.
*
*   MFILE  : handle to an file where output is written to
*   NLMIN  : minimum level where calculation is allowed
*   NLMAX  : maximum level where calculation is allowed
*
*   MATDAT : array [1..SZN2MI,1..NNLEV] of integer
*            TNS2DMatrixParams-structures for level NLMIN..NLMAX, 
*            initialized with data
*   VECDAT : array [1..SZN2VI,1..NNLEV] of integer
*            TNS2DVectorParams-structures for level NLMIN..NLMAX. 
*
*   IPARAM : array [1..SZISDI] of integer
*   DPARAM : array [1..SZISDI] of double
*            Integer and double prec. parameter blocks that define the
*            behaviour of the nonstationary solver. 
*            sub-solver M020. 
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer and double prec. parameter block for the stationary 
*            sub-solver NSDEF2.
*
*   LUP    : Handle to array [1..NUVP] of double
*            Iteration vector
*   LRHS   : Handle to array [1..NUVP] of double
*            Right hand side vector
*
* Out:
*   All handles to vectors / matrices / ... in the above structures
*   are invalid after calling this routine.
***********************************************************************
      
      SUBROUTINE DNISSL (NLMIN,NLMAX,
     *                   MATDAT,VECDAT,
     *                   IPARAM,DPARAM,
     *                   ISTPAR,DSTPAR,
     *                   LUP,LRHS)
      
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'm020.inc'

      INCLUDE 'stiming.inc'
      INCLUDE 'snonstead.inc'

      INTEGER MSHOW,NLMIN,NLMAX
      INTEGER MATDAT(SZN2MI,*),VECDAT(SZN2VI,*),IPARAM(*),ISTPAR(*)
      INTEGER LUP,LRHS
      DOUBLE PRECISION DPARAM(*),DSTPAR(*)
      
C     Release the rest using the clean up routine of the stationary
C     solver

      CALL DNSTSL (NLMIN,NLMAX,
     *             MATDAT,VECDAT,ISTPAR,DSTPAR,
     *             LUP,LRHS)
      
      END
            
************************************************************************
* Generate vector information for the nonstationary solver
*
* This routine calculates all "missing" discretisation information for
* the call to NONST2. After this routine the solver can be called
* to solve the problem.
*
* As this preparation phase already needs information about the problem
* and the discretization, the routine has major access to the parameters
* of the solver, which have to be initialized prior to calling this
* routine by the appropriate initialization routine!
*
* In:
*   MFILE  : handle to an file where output is written to
*   NLMIN  : minimum level where calculation is allowed
*   NLMAX  : maximum level where calculation is allowed
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for the underlying meshes
*            on all levels. TRIAS(.,ILEV) represents the triangulation
*            structure on level ILEV in its ummodified form, i.e.
*            without any boundary conditions implemented, probably
*            undeformed (if the user wants to start the grid
*            deformation in every time step from the original grids).
*   VECDAT : array [1..SZN2VI,1..NNLEV] of integer
*            TNS2DVectorParams-structures for level NLMIN..NLMAX. 
*
*   IPARAM : array [1..*] of integer
*   DPARAM : array [1..*] of double
*            Integer and double prec. parameter block of the solver.
*            For nonstationary solver, this points to the structure of
*            the nonstationary solver. For stationary solver, this points
*            to the structure of the stationary solver and therefore
*            coincides with ISTPAR/DSTPAR!
*   ISTPAR : array [1..SZNSDI] of integer
*   DSTPAR : array [1..SZNSDD] of double
*            Integer and double prec. parameter block for the stationary 
*            sub-solver NSDEF2.
*   IASMBL : array [1..SZASMI] of integer
*   DASMBL : array [1..SZASMD] of double
*            Integer and double prec. parameter block that controls the
*            discretization. 
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to boundary
*            routines. Not used in this routine.
*
*
*   TIMENS : Simulation time; should be 0D0 for preparation of first
*            time step
*   TSTEP  : Length of first time step; can usually be found in the
*            time step control structure.
*            Can be set to 1D0 if IDBTYP=0 for the preparation of the
*            very first timne step of a nonstationary simulation.
*
*   LUP    : Handle to array [1..NUVP] of double
*            Iteration vector; initialized with 0 or preinitialized
*            by reading a solution from disc
*   NEQ    : Length of solution vector
*
*   IBDTYP : Type of boundary conditions to implement into RHS / solution
*            vector. Bitfield.
*            =0: Don't implement boundary conditions;
*                TSTEP is ignored.
*            Bit 0: implement dirichlet boundary conditions
*                   into solution vectors
*
* Out:
*   - Implements boundary conditions into solution if desired
************************************************************************

      SUBROUTINE GNISVC (NLMIN,NLMAX,
     *                   TRIAS,
     *                   VECDAT,
     *                   IPARAM,DPARAM,
     *                   ISTPAR,DSTPAR,
     *                   IASMBL,DASMBL,IGEOM,DGEOM,
     *                   IBDTYP,TIMENS,TSTEP,
     *                   LUP)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'

      INCLUDE 'stria.inc'
      
      INCLUDE 'sassembly.inc'

      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'

      INCLUDE 'cbasicmg.inc'
      INCLUDE 'm020.inc'

      INCLUDE 'stiming.inc'
      INCLUDE 'snsdef.inc'
      INCLUDE 'ststepping.inc'

C     parameters

      INTEGER NLMIN,NLMAX,TRIAS(SZTRIA,*)
      INTEGER VECDAT(SZN2VI,*),ISTPAR(*),IGEOM(*)
      INTEGER IPARAM(*),IASMBL(*)
      INTEGER LUP,IBDTYP
      DOUBLE PRECISION DPARAM(*),DSTPAR(*),DASMBL(*),DGEOM(*)
      DOUBLE PRECISION TIMENS,TSTEP

C     local variables

      INTEGER KU,NEQU
      
C     externals

      DOUBLE PRECISION PARX,PARY,TMAX,UE
      EXTERNAL PARX,PARY,TMAX,UE

C     -----------------------------------------------------------------
C     Step 1: Implement boundary conditions into solution vector
C             on finest level

      KU   = L(LUP)
      NEQU = VECDAT(ONU,NLMAX)

      IF (IAND(IBDTYP,1).NE.0) THEN
        CALL BDRSTX (DWORK(KU),DWORK(KU+NEQU),
     *               KWORK(L(TRIAS(OLXNPR,NLMAX))),
     *               PARX,PARY,UE,TIMENS,
     *               DASMBL(ORE),
     *               TRIAS(1,NLMAX),IASMBL,DASMBL,IGEOM,DGEOM) 
      END IF

      END 

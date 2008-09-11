************************************************************************
* This file contains the standard routine for processing to the next
* time step in the nonstationary solver.
* The user can write its own postprocessing routine and call this
* old one to perform standard postprocessing.
************************************************************************

************************************************************************
* Proceed to next time step in nonstationary solver
*
* This is the default callback routine for proeeding to the next
* time step in the nonstationary Navier-Stokes solver NONST. 
* It's called 
*  - on the beginning of the simulation (with TIMENS=TIMEST and TSTEP=0)
*    and
*  - directly before setting the simulational time from
*    TIMENS to TIMENS+TSTEP.
*
* The basic routine is empty at the moment - this interface is
* only given as a documentation about the parameters. The main
* program normally redefines this routine and provides the solver
* with that modified version.
*
* In:
*   NLMIN  : minimum level 
*   NLMAX  : maximum level 
*   TRIAS  : array [1..SZTRIA,1..NLEV] of integer
*            Triangulation structures for all levels.
*   MATDAT : array [1..SZN2MI,1..NNLEV] of integer
*            TNS2DMatrixParams-structures for level NLMIN..NLMAX, 
*            initialized with data
*   VECDAT : array [1..SZN2VI,1..NNLEV] of integer
*            TNS2DVectorParams-structures for level NLMIN..NLMAX. 
*            This structure array must specify the structure of
*            the vectors on each level. 
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
*   IGEOM  - array [1..*] of integer 
*   DGEOM  - array [1..*] of double 
*            Integer- and double-precision parameter blocks with
*            geometry information. Passed to fictitious boundary
*            routines. Not used in this routine.
*
*   NEQU   : length of velocity vector (number of equations in u)
*   NEQP   : length of pressure vector (numer of equations in p)
*   NUVP   : Total length of solution vector; usually 2*NEQU+NEQP
*   DUP    : array [1..NUVP] of double
*            Current solution vector
*   DRHS   : array [1..NUVP] of double
*            Right hand side for the next time step
*   DAUX   : array [1..NUVP] of double
*            Auxiliary vector.
*
*   TIMEST - Initial simulation time
*   TIMENS - Current simulation time
*   TSTEP  - Current time step size. Can vary from time-step to
*            time-step because of the Fractional-Step scheme
*            or the adaptive time stepping.
*   ICORST - Number of substep in the time stepping scheme
*   NCORST - Total number of substeps in the time stepping scheme
*
* Out:
*   Nothing.
************************************************************************
      
      SUBROUTINE DFPRNT (NLMIN,NLMAX,
     *                   TRIAS,MATDAT,VECDAT,
     *                   IPARAM,DPARAM,ISTPAR,DSTPAR,
     *                   IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *                   NUVP,DUP,DRHS,DAUX,
     *                   TIMEST, TIMENS, TSTEP)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      
      INCLUDE 'stria.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'smat2dns.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
      INCLUDE 'ssolvers.inc'
      INCLUDE 'stiming.inc'
      INCLUDE 'snonstead.inc'
      
      INCLUDE 'sassembly.inc'
      
C parameters

      INTEGER NUVP
      INTEGER IPARAM (*),IMGPAR(*),ISTPAR(*),IASMBL(*),IGEOM(*)
      INTEGER TRIAS(SZTRIA,NNLEV)
      INTEGER MATDAT(SZN2MI,NNLEV),VECDAT(SZN2VI,NNLEV),NLMIN,NLMAX
      DOUBLE PRECISION DPARAM (*),DMGPAR(*),DSTPAR(*),DASMBL(*)
      DOUBLE PRECISION DGEOM(*)
      DOUBLE PRECISION TIMEST, TIMENS, TSTEP
      
      DOUBLE PRECISION DUP(*),DRHS(*),DAUX(*)
      
      END

!##############################################################################
!# ****************************************************************************
!# <name> LinearSolverAutoInitialise </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to initialise a linear solver for a given
!# situation. On one hand, it provides initialisation routines for standard
!# solver settings. On the other hand, there's a parser included which can
!# read user-defined text files from hard disc that configure the setting of
!# a solver.
!#
!# Hint: Reading from a file is not yet implemented.
!# </purpose>
!##############################################################################

MODULE linearsolverautoinitialise

  USE fsystem
  USE linearsolver

  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE linsolinit_CC2DMultigrid(p_rsolverNode, isolverType, nlevels, &
                        niterationsMin, niterationsMax,daccuracyAbs,daccuracyRel,  &
                        nprecSteps, nsmoothingSteps, &
                        nmaxCoarseGridSteps, dcoarseGridAccuracyAbs, &
                        dcoarseGridAccuracyRel)
  
  !<description>
  
  ! This routine builds a solver node for the linear solver for CC2D-like
  ! applications. The following combinations are supported:
  ! 1.) Multigrid solver, VANCA smoother, VANCA coarse grid solver
  ! 2.) BiCGStab-solver, Multigrid preconditioner, 
  !     VANCA smoother, VANCA coarse grid solver
  ! The caller must attach the matrices of all levels to the solver manually
  ! by calling linsol_setMatrices. Afterwards the linear solver can be started.
  
  !</description>
  
  !<input>
  
  ! Type of solver structure which should be build up.
  ! 1 = Multigrid solver, VANCA smoother, VANCA coarse grid solver
  ! 2 = BiCGStab-solver, Multigrid preconditioner, VANCA smoother,
  !     VANCA coarse grid solver
  INTEGER, INTENT(IN)               :: isolverType

  ! Number of levels
  INTEGER, INTENT(IN)               :: nlevels

  ! Minimum number of solver iterations
  INTEGER, INTENT(IN)               :: niterationsMin

  ! Maximum number of solver iterations
  INTEGER, INTENT(IN)               :: niterationsMax

  ! If solver conbination is BiCGStab with MG preconditioning, number of
  ! steps multigrid should perform. Otherwise ignored
  INTEGER, INTENT(IN)               :: nprecSteps

  ! Absolute accuracy of the solver
  REAL(DP), INTENT(IN)                   :: daccuracyAbs

  ! Relativew accuracy of the solver
  REAL(DP), INTENT(IN)                   :: daccuracyRel

  ! Number of pre- and postsmoothing steps
  INTEGER, INTENT(IN)               :: nsmoothingSteps
  
  ! Absolute accuracy on the coarse grid
  REAL(DP), INTENT(IN)                   :: dcoarseGridAccuracyAbs

  ! Relative accuracy on the coarse grid
  REAL(DP), INTENT(IN)                   :: dcoarseGridAccuracyRel
  
  ! Maximum iterations on the coarse grid
  INTEGER, INTENT(IN)               :: nmaxCoarseGridSteps

  ! A list of spatial discretisation structures for the three equations
  ! that are supported by CC2D: x-velocity, y-velocity, pressure
  TYPE(t_spatialDiscretisation), DIMENSION(3) :: RspatialDiscretisation
  
!</input>
  
!<output>
  
  ! The solver node that identifies the solver. The caller must attach
  ! level-dependent data (matrix information) to it with the standard
  ! routines in the module LinearSolver. Afterwards the structure
  ! can be used to solve the problem.
  
  TYPE(t_linsolNode),POINTER :: p_rsolverNode
  
!</output>
  
!</subroutine>

  ! local variables
  TYPE(t_linsolNode),POINTER :: p_rmgSolver
  INTEGER               :: ilevel
  TYPE(t_linsolNode),POINTER :: p_rpreSmoother
  TYPE(t_linsolNode),POINTER :: p_rpostSmoother
  TYPE(t_linsolNode),POINTER :: p_rcoarseGridSolver
  TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo
  TYPE(t_interlevelProjectionBlock) :: rprojection 
  
  ! Create the solver node - either BiCGStab or Multigrid.
  ! If we create BiCGStab, attach multigrid as preconditioner.
  
  CALL linsol_initMultigrid (p_rmgSolver)
  
  IF (isolverType .EQ. 1) THEN
  
    p_rsolverNode => p_rmgSolver
    
  ELSE
  
    CALL linsol_initBiCGStab (p_rsolverNode,p_rmgSolver)
    
    ! Configure MG preconditioning as a fixed number of MG steps
    ! without checking the residuum.
    
    p_rmgSolver%nminIterations = nprecSteps
    p_rmgSolver%nmaxIterations = nprecSteps
    p_rmgSolver%iresCheck      = NO
    
  END IF
  
  ! Configure the solver
  p_rsolverNode%nminIterations = niterationsMin
  p_rsolverNode%nmaxIterations = niterationsMax
  p_rsolverNode%depsRel = daccuracyRel
  p_rsolverNode%depsAbs = daccuracyAbs
  
  ! Initialise a standard interlevel projection structure for all levels.
  CALL mlprj_initProjection (rprojection,RspatialDiscretisation)
  
  ! Continue to configure MG by accessing p_rmgSolver.
  ! Loop through the levels.
  
  DO ilevel = 1,nlevels
  
    ! On the lowest level create a coarse grid solver structure
    IF (ilevel .EQ. 1) THEN
      CALL linsol_initVANCAQ1TP0NS2D (p_rcoarseGridSolver)
      p_rcoarseGridSolver%depsRel = dcoarseGridAccuracyRel
      p_rcoarseGridSolver%depsAbs = dcoarseGridAccuracyAbs
      p_rcoarseGridSolver%nmaxIterations = nmaxCoarseGridSteps
    ELSE
      NULLIFY(p_rcoarseGridSolver)
    END IF
    
    ! Create pre- and postsmoother structure on the current level
    CALL linsol_initVANCAQ1TP0NS2D (p_rpreSmoother)
    CALL linsol_initVANCAQ1TP0NS2D (p_rpostSmoother)
    
    ! Configure the structures to form a smoother. A smoother is a solver
    ! that iterates a finite number of times without respecting the
    ! residuum.
    
    CALL linsol_convertToSmoother (p_rpreSmoother,nsmoothingSteps)
    CALL linsol_convertToSmoother (p_rpostSmoother,nsmoothingSteps)
    
    ! Create the level, attrach it to the solver and proceed to the 
    ! next level
    CALL linsol_addMultigridLevel (p_rlevelInfo,p_rmgSolver, rprojection,&
                    p_rpresmoother,p_rpostsmoother,p_rcoarseGridSolver)
    
  END DO
    
  END SUBROUTINE
  
END MODULE

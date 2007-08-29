!##############################################################################
!# ****************************************************************************
!# <name> heatcond_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic problem structures for the heat conduction
!# problem.
!# </purpose>
!##############################################################################

MODULE heatcond_basic

  USE fsystem
  USE linearsolver
  USE boundary
  USE matrixfilters
  USE vectorfilters
  USE triangulation
  USE spatialdiscretisation
  USE sortstrategy
  USE timestepping
  
  USE collection
  USE paramlist
    
  IMPLICIT NONE
  
  ! Maximum allowed level in this application; must be =9 for 
  ! FEAT 1.x compatibility (still)!
  INTEGER, PARAMETER :: NNLEV = 9

!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (trial/test functions,...)
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! The static matrix (containing Laplace, convection,...) which does not
    ! change in time.
    TYPE(t_matrixBlock) :: rmatrixStatic
    
    ! The mass matrix
    TYPE(t_matrixBlock) :: rmatrixMass

    ! System matrix. May change during the time iteration    
    TYPE(t_matrixBlock) :: rmatrix
    
    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
  
  END TYPE
  
!</typeblock>

!<typeblock description="Configuration block for the time stepping.">

  TYPE t_problem_nonst
  
    ! Configuration block of the time stepping scheme.
    TYPE(t_explicitTimeStepping)        :: rtimestepping
    
    ! Number of current time step
    INTEGER                             :: iiteration
    
    ! Maximum number of time steps
    INTEGER                             :: niterations
    
    ! Start time
    REAL(DP)                            :: dtimemin
    
    ! Current time
    REAL(DP)                            :: dtime
    
    ! Maximum time
    REAL(DP)                            :: dtimemax
  
  END TYPE

!</typeblock>

!<typeblock description="Application-specific type block for heatcond problem">

  TYPE t_problem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    INTEGER :: ilvmin
    
    ! Maximum refinement level
    INTEGER :: ilvmax

    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary

    ! A variable describing the analytic boundary conditions.    
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions

    ! A RHS vector on the finest level used for solving linear systems
    TYPE(t_vectorBlock) :: rrhs
    
    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock) :: rprojection

    ! A filter chain to filter vectors during the solution process
    TYPE(t_filterChain), DIMENSION(1) :: RfilterChain

    ! A parameter block for everything that controls the time dependence.
    TYPE(t_problem_nonst) :: rtimedependence

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by NLMAX!
    TYPE(t_problem_lvl), DIMENSION(NNLEV) :: RlevelInfo
    
    ! A collection object that saves structural data and some 
    ! problem-dependent information which is e.g. passed to 
    ! callback routines.
    TYPE(t_collection) :: rcollection
    
  END TYPE

!</typeblock>

!</types>

END MODULE
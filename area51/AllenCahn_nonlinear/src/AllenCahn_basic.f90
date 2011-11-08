!##############################################################################
!# ****************************************************************************
!# <name> AllenCahn_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic problem structures for the heat conduction
!# problem.
!# </purpose>
!##############################################################################

MODULE AllenCahn_basic

  use fsystem
  use linearsolver
  use boundary
  use matrixfilters
  use vectorfilters
  use triangulation
  use spatialdiscretisation
  use sortstrategy
  use timestepping
  use discretebc
  use linearsystemblock
  use multilevelprojection
  use filtersupport
  
  use collection
  use paramlist
    
  IMPLICIT NONE
  
!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_ACproblem_lvl
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (trial/test functions,...)
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
!~~~~~~~MCai
    ! An object specifying the discretisation (trial/test functions,...)
    TYPE(t_blockDiscretisation)  :: rdiscretisation
!~~~~~
    ! The static matrix (containing Laplace, convection,...) which does not
    ! change in time.
    TYPE(t_matrixBlock) :: rmatrixStatic

    TYPE(t_matrixBlock) :: rmatrixPoly

    TYPE(t_matrixBlock) :: rmatrixJacobian

    TYPE(t_matrixBlock) :: rmatrixConv
        
    ! The mass matrix
    TYPE(t_matrixBlock) :: rmatrixMass

    ! System matrix. May change during the time iteration
    TYPE(t_matrixBlock) :: rmatrix
    
    ! A variable describing the discrete boundary conditions.
    ! The variable points to NULL until the first boundary conditions
    ! are assembled.
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC => NULL()
  
  END TYPE
  
!</typeblock>

!<typeblock description="Configuration block for the time stepping.">

  TYPE t_ACproblem_nonst
  
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

!<typeblock description="Application-specific type block for AllenCahn problem">

  TYPE t_ACproblem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    INTEGER :: NLMIN
    
    ! Maximum refinement level
    INTEGER :: NLMAX

    ! An object for saving the domain:
    TYPE(t_boundary) :: rboundary

    ! A RHS vector on the finest level used for solving linear systems
    TYPE(t_vectorBlock) :: rrhs
    
    ! A solver node that accepts parameters for the linear solver
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock) :: rprojection

    ! A filter chain to filter vectors during the solution process
    TYPE(t_filterChain), DIMENSION(1) :: RfilterChain

    ! A parameter block for everything that controls the time dependence.
    TYPE(t_ACproblem_nonst) :: rtimedependence

    ! An array of t_ACproblem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by NLMAX!
    TYPE(t_ACproblem_lvl), DIMENSION(:), POINTER :: RlevelInfo
    
    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    TYPE(t_collection) :: rcollection
    
  END TYPE

!</typeblock>

!</types>

END MODULE
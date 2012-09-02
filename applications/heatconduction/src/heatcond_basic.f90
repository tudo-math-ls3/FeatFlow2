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

module heatcond_basic

  use fsystem
  use linearsolver
  use boundary
  use matrixfilters
  use vectorfilters
  use triangulation
  use spatialdiscretisation
  use sortstrategybase
  use sortstrategy
  use timestepping
  use discretebc
  use linearsystemscalar
  use linearsystemblock
  use filtersupport
  use multilevelprojection
  
  use collection
  use paramlist
    
  implicit none
  
!<types>

!<typeblock description="Type block defining all information about one level">

  type t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (trial/test functions,...)
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubatureInfo

    ! The static matrix (containing Laplace, convection,...) which does not
    ! change in time.
    type(t_matrixBlock) :: rmatrixStatic
    
    ! The mass matrix
    type(t_matrixBlock) :: rmatrixMass

    ! System matrix. May change during the time iteration
    type(t_matrixBlock) :: rmatrix
    
    ! A variable describing the discrete boundary conditions.
    ! The variable points to NULL until the first boundary conditions
    ! are assembled.
    type(t_discreteBC), pointer :: p_rdiscreteBC => null()
  
    ! Sorting strategy for resorting vectors/matrices.
    type(t_blockSortStrategy) :: rsortStrategy

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1) :: RfilterChain
    
    ! Number of filters in the filter chain
    integer :: nfilters

  end type
  
!</typeblock>

!<typeblock description="Configuration block for the time stepping.">

  type t_problem_nonst
  
    ! Configuration block of the time stepping scheme.
    type(t_explicitTimeStepping)        :: rtimestepping
    
    ! Number of current time step
    integer                             :: iiteration
    
    ! Maximum number of time steps
    integer                             :: niterations
    
    ! Start time
    real(DP)                            :: dtimemin
    
    ! Current time
    real(DP)                            :: dtime
    
    ! Maximum time
    real(DP)                            :: dtimemax
  
  end type

!</typeblock>

!<typeblock description="Application-specific type block for heatcond problem">

  type t_problem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    integer :: ilvmin
    
    ! Maximum refinement level
    integer :: ilvmax

    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! A RHS vector on the finest level used for solving linear systems
    type(t_vectorBlock) :: rrhs
    
    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection

    ! A filter chain to filter vectors during the solution process
    type(t_filterChain), dimension(1) :: RfilterChain

    ! A parameter block for everything that controls the time dependence.
    type(t_problem_nonst) :: rtimedependence
    
    ! Directory for output files
    character(len=SYS_STRLEN) :: sucddir
    
    ! File to the parametrisation and triangulation
    character(len=SYS_STRLEN) :: sprmfile,strifile

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by NLMAX!
    type(t_problem_lvl), dimension(:), pointer :: RlevelInfo
    
    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    type(t_collection) :: rcollection
    
  end type

!</typeblock>

!</types>

end module

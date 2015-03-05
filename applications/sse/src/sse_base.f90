!##############################################################################
!# ****************************************************************************
!# <name> sse_base </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the basic routines and constants.
!# </purpose>
!##############################################################################

module sse_base

  use fsystem
  use genoutput

  use boundary, only : t_boundary
  use collection, only : t_collection
  use discretebc, only : t_discreteBC
  use filtersupport, only : t_filterChain
  use linearsystemblock, only : t_matrixBlock,t_vectorBlock
  use linearsystemscalar, only : t_matrixScalar
  use multilevelprojection, only : t_interlevelProjectionBlock
  use sortstrategybase, only : t_blockSortStrategy
  use spatialdiscretisation, only : t_blockDiscretisation,t_scalarCubatureInfo
  use triangulation, only : t_triangulation
  
  implicit none

  private

  public :: t_problem
  public :: t_problem_lvl
  public sse_getNDIM
  public sse_getNVAR
  public sse_getSection
  public sse_getSolver
  
#ifdef USE_COMPILER_INTEL
  public :: sinh
  public :: cosh
#endif

!<constants>

!<constantblock description="Constants for problem types">

  ! Compute standard Poisson problem in 2D
  integer, parameter, public :: POISSON_SCALAR = 0

  ! Compute standard Poisson problem in 2D as first-order system
  integer, parameter, public :: POISSON_SYSTEM = 1
  
  ! Compute SSE solution from scalar problem in 2D
  integer, parameter, public :: SSE_SCALAR     = 2

  ! Compute SSE solution from first-order system ($\sigma=A\nabla N$) in 2D
  integer, parameter, public :: SSE_SYSTEM1    = 3

  ! Compute SSE solution from first-order system ($\sigma=\nabla N$) in 2D
  integer, parameter, public :: SSE_SYSTEM2    = 4

  ! Compute Corine`s problem in 1D
  integer, parameter, public :: CORINE_1D      = 5

  ! Compute Corine`s problem in 2D
  integer, parameter, public :: CORINE_2D      = 6
!</constantblock>
  
!<constantblock description="Constants for complex numbers">

  ! Real part of complex number
  complex(DP), parameter, public :: creal = cmplx(1.0_DP,0.0_DP)

  ! Imaginary part of complex number
  complex(DP), parameter, public :: cimg = cmplx(0.0_DP,1.0_DP)

!</constantblock>

!<constantblock description="Constants for solvers">

  ! Scalar solver
  integer, parameter, public :: SOLVER_SCALAR       = 0

  ! System solver
  integer, parameter, public :: SOLVER_SYSTEM       = 1
  
  ! Saddle point solver
  integer, parameter, public :: SOLVER_SADDLEPOINT  = 2
  
!</constantblock>
  
!</constants>

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_problem_lvl

    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation

    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo), dimension(3) :: RcubatureInfo

    ! A system matrix for that specific level. The matrix will receive the
    ! discrete system operator.
    type(t_matrixBlock) :: rmatrix

    ! A scalar matrix that will recieve the prolongation matrix for this level.
    type(t_matrixScalar), dimension(6) :: rmatProl,rmatRest

    ! An interlevel projection structure for changing levels (scalar case)
    type(t_interlevelProjectionBlock) :: rprojection

    ! Interlevel projection structures for changing levels (system case)
    type(t_interlevelProjectionBlock) :: rprojectionA
    type(t_interlevelProjectionBlock) :: rprojectionS

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC

    ! Sorting strategy for resorting vectors/matrices.
    type(t_blockSortStrategy) :: rsortStrategy

    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    type(t_filterChain), dimension(1) :: RfilterChain

    ! Number of filters in the filter chain.
    integer :: nfilters

  end type

!</typeblock>

!<typeblock description="Application-specific type block for SSE problem">

  type t_problem

    ! Problem type
    integer :: cproblemtype

    ! Problem subtype
    integer :: cproblemsubtype
    
    ! Minimum refinement level; = Level i in RlevelInfo
    integer :: ilvmin

    ! Maximum refinement level
    integer :: ilvmax

    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! A solution vector and a RHS vector on the finest level.
    type(t_vectorBlock) :: rvector,rrhs

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation.
    type(t_problem_lvl), dimension(:), pointer :: RlevelInfo

    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    type(t_collection) :: rcollection

  end type

!</typeblock>

!</types>
  
contains

  ! ***************************************************************************
  
#ifdef USE_COMPILER_INTEL
!<function>

  elemental function sinh(cx)

!<description>
    ! Complex valued hyperbolic sine functions (available in Fortran 2008)
!</description>

!<input>
    complex(DP), intent(in) :: cx
!</input>

!<result>
    complex(DP) :: sinh
!</result>
!</function>

    sinh = -cmplx(0.0_DP,1.0_DP) * sin(cmplx(0.0_DP,1.0_DP)*cx) 

  end function
#endif
  
  ! ***************************************************************************
  
#ifdef USE_COMPILER_INTEL
!<function>

  elemental function cosh(cx)

!<description>
    ! Complex valued hyperbolic sine functions (available in Fortran 2008)
!</description>

!<input>
    complex(DP), intent(in) :: cx
!</input>

!<result>
    complex(DP) :: cosh
!</result>
!</function>

    cosh = cos(cmplx(0.0_DP,1.0_DP)*cx) 

  end function
#endif

  ! ***************************************************************************

!<function>

  function sse_getNDIM(cproblemtype) result(ndim)

!<description>
    ! This function returns the number of spatial dimensions
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
!</input>

!<result>
    ! Number of spatial dimensions
    integer :: ndim
!</result>
!</function>

    select case(cproblemtype)
    case(POISSON_SCALAR, POISSON_SYSTEM)
      ndim = 2

    case(SSE_SCALAR, SSE_SYSTEM1,SSE_SYSTEM2)
      ndim = 2

    case (CORINE_1D)
      ndim = 1

    case (CORINE_2D)
      ndim = 2

    case default
      ndim = 0
      call output_line("Invalid problem type", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_getNDIM")
      call sys_halt()
    end select
  end function sse_getNDIM
  
  ! ***************************************************************************

!<function>

  function sse_getNVAR(cproblemtype) result(nvar)

!<description>
    ! This function returns the number of variables
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
!</input>

!<result>
    ! Number of variables
    integer :: nvar
!</result>
!</function>

    select case(cproblemtype)
    case(POISSON_SCALAR)
      nvar = 1

    case(POISSON_SYSTEM)
      nvar = 1 + sse_getNDIM(cproblemtype)

    case(SSE_SCALAR)
      nvar = 1 * 2

    case(SSE_SYSTEM1,SSE_SYSTEM2)
      nvar = (1 + sse_getNDIM(cproblemtype)) * 2

    case (CORINE_1D)
      nvar = 6

    case (CORINE_2D)
      nvar = 0

    case default
      nvar = 0

      call output_line("Invalid problem type", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_getNVAR")
      call sys_halt()
    end select
  end function sse_getNVAR

  ! ***************************************************************************

!<function>

  function sse_getSection(cproblemtype) result(sstring)

!<description>
    ! This function returns the section name
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
!</input>

!<result>
    ! Section name
    character(SYS_STRLEN) :: sstring
!</result>
!</function>

    select case(cproblemtype)
    case(POISSON_SCALAR, POISSON_SYSTEM)
      sstring = 'POISSON'

    case(SSE_SCALAR, SSE_SYSTEM1,SSE_SYSTEM2)
      sstring = 'SSE'

    case (CORINE_1D,CORINE_2D)
      sstring = 'CORINE'

    case default
      sstring = ''

      call output_line("Invalid problem type", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_getSection")
      call sys_halt()
    end select
  end function sse_getSection

  ! ***************************************************************************

!<function>

  function sse_getSolver(cproblemtype) result(isolver)

!<description>
    ! This function returns the type of solver
!</description>

!<input>
    ! Problem type
    integer, intent(in) :: cproblemtype
!</input>

!<result>
    ! Section name
    integer :: isolver
!</result>
!</function>

    select case(cproblemtype)
    case(POISSON_SCALAR, SSE_SCALAR)
      isolver = SOLVER_SCALAR

    case(POISSON_SYSTEM, SSE_SYSTEM1, SSE_SYSTEM2)
      ! isolver = SOLVER_SADDLEPOINT !!! not working at the moment
      isolver = SOLVER_SYSTEM

    case (CORINE_1D,CORINE_2D)
      isolver = SOLVER_SCALAR

    case default
      isolver = -1
      
      call output_line("Invalid problem type", &
          OU_CLASS_ERROR,OU_MODE_STD,"sse_getSolver")
      call sys_halt()
    end select
  end function sse_getSolver
  
end module sse_base

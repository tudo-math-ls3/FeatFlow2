!##############################################################################
!# ****************************************************************************
!# <name> filtersupport </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains filter routines that can be used in the solution
!# process of the problem (especially when solving linear subproblems).
!# For using filters in a solver, a ''filter chain'' must be build up.
!# This is a simple array of filter structures that describe the filter
!# that should be applied to a vector one after the other.
!#
!# To build a filter chain, the application must declare an array
!# of type t_filterChain (e.g. of length FILTER_MAXFILTERS) and set the
!# filer type identifier.
!# The routine filter_applyFilterChainVec can be used to apply such a filter
!# chain to a vector. The routine stops when either the end of the array
!# is reached or if a filter is found that contains the tag FILTER_NOFILTER.
!#
!# A filter is always applied to the whole solution/defect vector. One filter
!# is for example the "apply-discrete-boundary-conditions" filter
!# that loops about all subvectors of a given vector to implement the
!# discrete boundary conditions.
!#
!# 'Linear' and 'nonlinear' filters \\
!# -------------------------------- \\
!# There exist two types of filters: 'Linear' and 'nonlinear' ones:
!#
!# a) 'Linear' filters can be applied by a filter chain when solving
!#   a linear system. The routines in this file provide the functionality
!#   to apply such a filter chain to a vector or a matrix.
!#   Examples for such filters are:
!#   - Implement Dirichlet boundary conditions into a matrix/vector,
!#   - Filter a vector to be in <tex>$L^2_0$</tex>.
!#
!# b) 'Nonlinear' filters are usually used inside of a nonlinear loop.
!#   Filters of this type are somehow 'special', as they usually need
!#   extended information. These filters cannot be incorporated in a
!#   filter chain and can therefore only be called manually and not
!#   during the solution process of a linear system.
!#   Examples for such filters are:
!#   - Implementation of slip boundary conditions.
!#
!# The following routines can be found here:
!#
!# 1.) filter_applyFilterChainVec
!#     -> Applies a given filter chain onto a (block) vector.
!#
!# 2.) filter_applyFilterChainMat
!#     -> Applies a given filter chain onto a (block) matrix.
!#
!# How to use filters \\
!# ------------------ \\
!# Using filters is rather easy:
!#
!# a) Declare an array of size FILTER_MAXFILTERS:
!#
!# <code>
!#      TYPE(t_filterChain), DIMENSION(FILTER_MAXFILTERS) :: RfilterChain
!# </code>
!#
!# b) All filters in this chain are automatically initialised to type
!#    FILTER_NOFILTER by the framework. Initialise the first couple of
!#    entries in that array to the filter you need, e.g.:
!#
!# <code>
!#      RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
!#      RfilterChain(2)%ifilterType = FILTER_DISCBCDEFFICT
!#
!#      RfilterChain(3)%ifilterType = FILTER_TOL20
!#      RfilterChain(3)%itoL20component = 3
!# </code>
!#
!#    The above initialisation sets up a filter chain the implementation
!#    of Dirichlet boundary conditions on the real boundary and on
!#    fictitious boundary components to a vector and filters the 3rd
!#    component of a block vector into the space <tex>$L^2_0$</tex>.
!#
!# c) Either apply the filter chain to a matrix or vector, e.g.
!#
!# <code>
!#      filter_applyFilterChainVec (rx, RfilterChain)
!#      filter_applyFilterChainMat (rmatrix, RfilterChain)
!# </code>
!#
!#    or specify it when initialising a linear solver, e.g.
!#
!# <code>
!#      TYPE(t_filterChain), DIMENSION(FILTER_MAXFILTERS), POINTER :: p_filterChain
!#      p_filterChain => RfilterChain
!#      CALL linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_filterChain)
!# </code>
!#
!# Notes:
!# a) The filter chain is processed until FILTER_NOFILTER is reached.
!#    Setting ifilterType=FILTER_DONOTHING deactivates a filter.
!#    Setting ifilterType=FILTER_NOFILTER stops processing the filter chain
!#    at this position.
!#
!# b) There is no memory attached to filters. Therefore, no cleanup is
!#    necessary.
!#
!# </purpose>
!##############################################################################

module filtersupport

!$use omp_lib
  use fsystem
  use genoutput
  use linearsystemblock
  use spatialdiscretisation
  use matrixfilters
  use vectorfilters

  implicit none

  private

  public :: t_filterChain
  public :: filter_applyFilterChainVec
  public :: filter_applyFilterChainMat

!<constants>

!<constantblock description="General constants concerning filters">

  ! A standard length for arrays holding a filter chain.
  integer, parameter, public :: FILTER_MAXFILTERS        =  32

!</constantblock>

!<constantblock description="Filter type flags for t\_filterChain%ifilterType">

  ! Abort-filter. Causes the filter chain to stop filtering immediately.
  integer, parameter, public :: FILTER_ABORT             =  -1

  ! Do-nothing filter; does nothing
  integer, parameter, public :: FILTER_DONOTHING         =  0

  ! Vector filter for imposing discrete boundary conditions of the
  ! real boundary into a solution vector.
  integer, parameter, public :: FILTER_DISCBCSOLREAL     =  1

  ! Vector filter for imposing discrete boundary conditions of the
  ! real boundary into a right-hand-side vector.
  integer, parameter, public :: FILTER_DISCBCRHSREAL     =  2

  ! Vector filter for imposing discrete boundary conditions of the
  ! real boundary into a defect vector.
  integer, parameter, public :: FILTER_DISCBCDEFREAL     =  3

  ! Matrix filter for imposing discrete boundary conditions of the
  ! real boundary into a matrix.
  integer, parameter, public :: FILTER_DISCBCMATREAL     =  4

  ! Vector filter for imposing discrete boundary conditions of
  ! the fictitious boundary into a solution vector.
  integer, parameter, public :: FILTER_DISCBCSOLFICT     =  5

  ! Vector filter for imposing discrete boundary conditions of
  ! the fictitious boundary into a right-hand-side vector.
  integer, parameter, public :: FILTER_DISCBCRHSFICT     =  6

  ! Vector filter for imposing discrete boundary conditions of
  ! the fictitious boundary into a defect vector.
  integer, parameter, public :: FILTER_DISCBCDEFFICT     =  7

  ! Matrix filter for imposing discrete boundary conditions of the
  ! fictitious boundary into a matrix.
  integer, parameter, public :: FILTER_DISCBCMATFICT     =  8

  ! Vector filter for bringing a subvector of a vector to the space <tex>$L^2_0$</tex>.
  integer, parameter, public :: FILTER_TOL20             =  9

  ! Vector filter for bringing the vector sum (small l1-norm) to mean value 0.
  integer, parameter, public :: FILTER_SMALLL1TO0        =  10

!</constantblock>

!</constants>

!<types>

  !<typeblock>

  ! A structure for using filters in a solver. This structure contains all
  ! information that are necessary to apply a filter (level independent!).

  type t_filterChain

    ! Tag that identifies the type of filter to be applied to a vector.
    ! One of the FILTER_XXXX constants
    integer                            :: ifilterType = FILTER_ABORT

    ! Information tag for the TOL20 filter if ifilterType=FILTER_TOL20:
    ! Number of the subvector that should be filtered to be in the
    ! space <tex>$L^2_0$</tex>.
    integer                            :: itoL20component = 0

    ! Information tag for the SMALLL1TOL0 filter if ifilterType=FILTER_SMALLL1TO0:
    ! Number of the subvector that should be filtered.
    integer                            :: ismallL1to0component = 0

  end type

  !</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine filter_applyFilterChainVec (rx, RfilterChain)

!<description>
  ! This routine applies a filter chain on a (block) vector rx. All filters in
  ! the chain are applied one after the other until the end of the array is
  ! reached or the first filter is found with the filter tag FILTER_NOFILTER.
!</description>

!<input>

  ! The filter chain
  type(t_filterChain), dimension(:), intent(in) :: RfilterChain

!</input>

!<inputoutput>

  ! The vector where the filter will be applied to.
  ! This is also the result vector.
  type(t_vectorBlock), intent(inout)            :: rx

!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i, ifilterType

    ! Apply the filters to rx - one after the other
    do i=lbound(RfilterChain,1),ubound(RfilterChain,1)

      ifilterType = RfilterChain(i)%ifilterType

      ! Choose the filter and apply it
      select case (ifilterType)
      case (FILTER_ABORT)
        ! Cancel if we reached the last filter before reaching the end of the
        ! array
        exit

      case (FILTER_DONOTHING)
        ! Do nothing

      case (FILTER_DISCBCSOLREAL)
        ! Impose Dirichlet boundary conditions into the solution vector rx
        call vecfil_discreteBCsol (rx)

      case (FILTER_DISCBCRHSREAL)
        ! Impose Dirichlet boundary conditions into the RHS vector rx
        call vecfil_discreteBCrhs (rx)

      case (FILTER_DISCBCDEFREAL)
        ! Impose Dirichlet boundary conditions into the defect vector rx
        call vecfil_discreteBCdef (rx)

      case (FILTER_DISCBCSOLFICT)
        ! Impose Dirichlet fictitious boundary conditions into the solution vector rx
        call vecfil_discreteFBCsol (rx)

      case (FILTER_DISCBCRHSFICT)
        ! Impose Dirichlet fictitious boundary conditions into the RHS vector rx
        call vecfil_discreteFBCrhs (rx)

      case (FILTER_DISCBCDEFFICT)
        ! Impose Dirichlet fictitious boundary conditions into the defect vector rx
        call vecfil_discreteFBCdef (rx)

      case (FILTER_TOL20)
        ! Bring the subvector itoL20component of rx to the space <tex>$L^2_0$</tex>:
        call vecfil_subvectorToL20 (rx,RfilterChain(i)%itoL20component)

      case (FILTER_SMALLL1TO0)
        ! Bring the subvector itoL20component of rx to the space <tex>$L^2_0$</tex>:
        call vecfil_subvectorSmallL1To0 (rx,RfilterChain(i)%ismallL1to0component)

      case default
        call output_line ('Unknown filter.', OU_CLASS_WARNING, OU_MODE_STD, &
                          'filter_applyFilterChainVec')
        exit

      end select

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine filter_applyFilterChainMat (rmatrix, RfilterChain)

!<description>
  ! This routine applies a filter chain on a (block) matrix rmatrix. All
  ! filters in the chain are applied one after the other until the end of
  ! the array is reached or the first filter is found with the filter
  ! tag FILTER_NOFILTER.
  !
  ! When implementing boundary conditions using filters, note that all
  ! boundary conditions that are implemented with 'nonlinear filters' are
  ! not implemented into the matrix with this routine - they must be
  ! implemented manually! The following boundary conditions are
  ! not implemented with a matrix filter chain:
  ! - Slip boundary conditions
!</description>

!<input>

  ! The filter chain
  type(t_filterChain), dimension(:), intent(in) :: RfilterChain

!</input>

!<inputoutput>

  ! The matrix where the filter will be applied to.
  ! This is also the result matrix.
  type(t_matrixBlock), intent(inout)            :: rmatrix

!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i, ifilterType

    ! Apply the filters to matrixx - one after the other
    do i=lbound(RfilterChain,1),ubound(RfilterChain,1)

      ifilterType = RfilterChain(i)%ifilterType

      ! Choose the filter and apply it
      select case (ifilterType)
      case (FILTER_ABORT)
        ! Cancel if we reached the last filter before reaching the end of the
        ! array
        exit

      case (FILTER_DISCBCSOLREAL,FILTER_DISCBCRHSREAL, &
            FILTER_DISCBCDEFREAL,FILTER_DISCBCMATREAL)
        ! Impose Dirichlet boundary conditions into the matrix rmatrix.
        ! The filter is the same for both, solution and defect filter,
        ! as the matrix modification is teh same (for now).
        call matfil_discreteBC (rmatrix)

      case (FILTER_DISCBCSOLFICT,FILTER_DISCBCRHSFICT, &
            FILTER_DISCBCDEFFICT,FILTER_DISCBCMATFICT)
        ! Impose Dirichlet fictitious boundary conditions into the matrix rmatrix.
        ! The filter is the same for both, solution and defect filter,
        ! as the matrix modification is teh same (for now).
        call matfil_discreteFBC (rmatrix)

      case (FILTER_TOL20,FILTER_DONOTHING)
        ! These filters have no effect for matrices.

      case default
        call output_line ('Unknown filter.', OU_CLASS_WARNING, OU_MODE_STD, &
                          'filter_applyFilterChainMat')
        exit

      end select

    end do

  end subroutine

end module

!##############################################################################
!# ****************************************************************************
!# <name> filtersupport </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains filter routines that can be used in the solution
!# process of the problem (especially when solving linear subproblems).
!# For using filters in a solver, a ""filter chain"" must be build up.
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
!# "Linear" and "nonlinear" filters \\
!# -------------------------------- \\
!# There exist two types of filters: "Linear" and "nonlinear" ones:
!#
!# a) "Linear" filters can be applied by a filter chain when solving
!#   a linear system. The routines in this file provide the functionality
!#   to apply such a filter chain to a vector or a matrix.
!#   Examples for such filters are:
!#   - Implement Dirichlet boundary conditions into a matrix/vector,
!#   - Filter a vector to be in <tex>$L^2_0$</tex>.
!#
!# b) "Nonlinear" filters are usually used inside of a nonlinear loop.
!#   Filters of this type are somehow "special", as they usually need
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
!# The following routines can be used to build up simple filter chains:
!#
!# 1.) filter_initFilterChain
!#     -> Creates a filter chain
!#
!# 1.) filter_clearFilterChain
!#     -> Clears a filter chain
!#
!# 1.) filter_doneFilterChain
!#     -> Releases a filter chain
!#
!# 4.) filter_newFilterDiscBCDef
!#     -> Adds a filter for discrete boundary conditions on the real boundary
!#
!# 5.) filter_newFilterDiscFBCDef
!#     -> Adds a filter for discrete boundary conditions on fict. boundary
!#
!# 6.) filter_newFilterToL20
!#     -> Adds a filter to enforce the L2_0 space for a component
!#
!# 7.) filter_newFilterSmallL1to0
!#     -> Adds a filter to enforce to enforce the l1-norm of a component
!#        to be zero
!#
!# 8.) filter_newFilterOneEntryZero
!#     -> Adds a filter to put one entry to zero
!#
!# 9.) filter_newFilterOverwriteDofs  
!#     -> Adds a filter to set a couple of DOFs to predefined values
!#
!#
!# How to use filters \\
!# ------------------ \\
!# Using filters is rather easy:
!#
!# a) Declare an array of size FILTER_MAXFILTERS:
!#
!# <code>
!#      type(t_filterChain), dimension(FILTER_MAXFILTERS) :: RfilterChain
!       integer :: nsize
!# </code>
!#
!# b) All filters in this chain are automatically initialised to type
!#    FILTER_NOFILTER by the framework. Initialise the first couple of
!#    entries in that array to the filter you need, e.g.:
!#
!# <code>
!#      call filter_newFilterChain (RfilterChain,nsize)
!#      call filter_newFilterDiscBCDef (RfilterChain,nsize,rdiscreteBC)
!#      call filter_newFilterDiscFBCDef (RfilterChain,nsize,rdiscreteFBC)
!#      call filter_newFilterToL20 (RfilterChain,nsize,3)
!#      ...
!#      call filter_doneFilterChain (RfilterChain,nsize)
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
!#      CALL linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
!# </code>
!#
!# Notes:
!# a) The filter chain is processed until FILTER_NOFILTER is reached.
!#    Setting ifilterType=FILTER_DONOTHING deactivates a filter.
!#    Setting ifilterType=FILTER_NOFILTER stops processing the filter chain
!#    at this position.
!#
!# b) There is no memory allocated in filters. Therefore, no cleanup is
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
  
  use discretebc
  use discretefbc

  implicit none

  private

  public :: t_filterChain
  public :: filter_applyFilterChainVec
  public :: filter_applyFilterChainMat
  
  public :: filter_initFilterChain
  public :: filter_clearFilterChain
  public :: filter_doneFilterChain
  
  public :: filter_newFilterDiscBCDef
  public :: filter_newFilterDiscFBCDef
  public :: filter_newFilterToL20
  public :: filter_newFilterSmallL1to0
  public :: filter_newFilterOneEntryZero
  public :: filter_newFilterOverwriteDofs  

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

  ! Vector filter for replacing one entry of a subvector with zero.
  integer, parameter, public :: FILTER_ONEENTRY0         =  11

  ! Vector filter for replacing entries of a subvector by predefined values.
  integer, parameter, public :: FILTER_OVERWRITEDOFS     = 12

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
    
    ! Information tag for the ONEENTRY0 filter if ifilterType=FILTER_ONEENTRY0:
    ! Number of the subvector that should be filtered.
    integer                            :: iblock = 0
    
    ! Information tag for the ONEENTRY0 filter if ifilterType=FILTER_ONEENTRY0:
    ! Number of the entry (row) of the subvector that should be filtered.    
    integer                            :: irow = 0
    
    ! Reference to a boundary condition structure for the Dirichlet boundary 
    ! condition filter or NULL if the predefined boundary conditions in the vector
    ! should be applied.
    type(t_discreteBC), pointer        :: p_rdiscreteBC => null()

    ! Reference to a boundary condition structure for the Dirichlet boundary 
    ! condition filter or NULL if the predefined boundary conditions in the vector
    ! should be applied. Fictitious boundary.
    type(t_discreteFBC), pointer       :: p_rdiscreteFBC => null()
    
    ! Number of DOFs to overwrite with the DOF overwrite filter.
    integer                            :: ndofs = 0

    ! Pointer to a list of DOFs to be overwritten by the DOF overwrite filter.
    integer, dimension(:), pointer     :: p_IdofsOverwrite => null()
    
    ! Pointer to a list of values to overwrite the DOFs with in the DOF overwrite
    ! filter, or NULL() if to overwrite by zero.
    real(DP), dimension(:), pointer     :: p_DdofsOverwrite => null()

  end type

  !</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine filter_initFilterChain (RfilterChain,nsize)

!<description>
  ! Initialises a filter chain.
!</description>

!<output>
  ! The filter chain to be resetted.
  type(t_filterChain), dimension(:), intent(out) :: RfilterChain
  
  ! Number of filters in the filter chain. Is set to zero.
  integer, intent(out) :: nsize
!</output>

!</subroutine>

    ! The "intent(out)" Resets the filter chain.
    ! Just reset nsize if specified.
    nsize = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine filter_clearFilterChain (RfilterChain,nsize)

!<description>
  ! Resets all filters in a filter chain.
!</description>

!<inputoutput>
  ! The filter chain to be resetted.
  type(t_filterChain), dimension(:), intent(inout) :: RfilterChain
  
  ! Number of filters in the filter chain.
  integer, intent(inout) :: nsize
!</inputoutput>

!</subroutine>

    call filter_doneFilterChain (RfilterChain,nsize)
    call filter_initFilterChain (RfilterChain,nsize)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine filter_doneFilterChain (RfilterChain,nsize)

!<description>
  ! Releases a filter chain
!</description>

!<inputoutput>
  ! The filter chain to be cleaned up.
  type(t_filterChain), dimension(:), intent(out) :: RfilterChain
  
  ! Number of filters in the filter chain.
  integer, intent(out) :: nsize
!</inputoutput>

!</subroutine>

    ! Currently, there is no memory allocated by the filter chain.
    ! So just reset it.
    call filter_initFilterChain (RfilterChain,nsize)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine filter_newFilterDiscBCDef (RfilterChain,nsize,rdiscreteBC)

!<description>
  ! Adds a "Discrete boundary condition to defect" filter for the "real boundary".
!</description>

!<input>
  ! OPTIONAL: Boundary condition structure to filter the defect with.
  ! If not specified, the default boundary condition structure in the
  ! vector is used (DEPRECATED!).
  type(t_discreteBC), intent(in), target, optional :: rdiscreteBC
!</input>

!<inputoutput>
  ! The filter chain.
  ! The new filter is appended at position nsize+1.
  ! The array must be large enough!
  type(t_filterChain), dimension(:), intent(inout) :: RfilterChain

  ! Current size of the filter chain.
  ! This number is increased by one and the new filter is appended.
  !
  ! For adding the first filter, this value must be set =0 !
  integer, intent(inout) :: nsize
!</inputoutput>

!</subroutine>

    if (size(RfilterChain) .lt. nsize+1) then
      call output_line ("Filter chain not large enough!", &
          OU_CLASS_ERROR, OU_MODE_STD, "filter_newFilterDiscBCDefReal")
      call sys_halt()
    end if
    
    ! Append the filter.
    nsize = nsize + 1
    RfilterChain(nsize)%ifilterType = FILTER_DISCBCDEFREAL
    nullify(RfilterChain(nsize)%p_rdiscreteBC)
    if (present(rdiscreteBC)) RfilterChain(nsize)%p_rdiscreteBC => rdiscreteBC

#if WARN_DEPREC
    if (.not. present(rdiscreteBC)) then
      call output_line ("Using deprecated feature. Please update your code.", &
          OU_CLASS_WARNING,OU_MODE_STD,"filter_newFilterDiscBCDef")
    end if
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine filter_newFilterDiscFBCDef (RfilterChain,nsize,rdiscreteFBC)

!<description>
  ! Adds a "Discrete boundary condition to defect" filter for the 
  ! "fixtitious boundary".
!</description>

!<input>
  ! OPTIONAL: Boundary condition structure to filter the defect with.
  ! If not specified, the default boundary condition structure in the
  ! vector is used (DEPRECATED!).
  type(t_discreteFBC), intent(in), target, optional :: rdiscreteFBC
!</input>

!<inputoutput>
  ! The filter chain.
  ! The new filter is appended at position nsize+1.
  ! The array must be large enough!
  type(t_filterChain), dimension(:), intent(inout) :: RfilterChain

  ! Current size of the filter chain.
  ! This number is increased by one and the new filter is appended.
  !
  ! For adding the first filter, this value must be set =0 !
  integer, intent(inout) :: nsize
!</inputoutput>

!</subroutine>

    if (size(RfilterChain) .lt. nsize+1) then
      call output_line ("Filter chain not large enough!", &
          OU_CLASS_ERROR, OU_MODE_STD, "filter_newFilterDiscBCDefReal")
      call sys_halt()
    end if
    
    ! Append the filter.
    nsize = nsize + 1
    RfilterChain(nsize)%ifilterType = FILTER_DISCBCDEFFICT
    nullify(RfilterChain(nsize)%p_rdiscreteFBC)
    if (present(rdiscreteFBC)) RfilterChain(nsize)%p_rdiscreteFBC => rdiscreteFBC

#if WARN_DEPREC
    if (.not. present(rdiscreteFBC)) then
      call output_line ("Using deprecated feature. Please update your code.", &
          OU_CLASS_WARNING,OU_MODE_STD,"filter_newFilterDiscFBCDef")
    end if
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine filter_newFilterToL20 (RfilterChain,nsize,itoL20component)

!<description>
  ! Adds a "Vector to L2_0 space" filter to a filter chain.
!</description>

!<input>
  ! Number of the component in the vector which should be filtered.
  integer, intent(in) :: itoL20component
!</input>

!<inputoutput>
  ! The filter chain.
  ! The new filter is appended at position nsize+1.
  ! The array must be large enough!
  type(t_filterChain), dimension(:), intent(inout) :: RfilterChain

  ! Current size of the filter chain.
  ! This number is increased by one and the new filter is appended.
  !
  ! For adding the first filter, this value must be set =0 !
  integer, intent(inout) :: nsize
!</inputoutput>

!</subroutine>

    if (size(RfilterChain) .lt. nsize+1) then
      call output_line ("Filter chain not large enough!", &
          OU_CLASS_ERROR, OU_MODE_STD, "filter_newFilterToL20")
      call sys_halt()
    end if
    
    ! Append the filter.
    nsize = nsize + 1
    RfilterChain(nsize)%ifilterType = FILTER_TOL20
    RfilterChain(nsize)%itoL20component = itoL20component

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine filter_newFilterSmallL1to0 (RfilterChain,nsize,ismallL1to0component)

!<description>
  ! Adds a "l1-norm to zero" filter to a filter chain.
!</description>

!<input>
  ! Number of the component in the vector which should be filtered.
  integer, intent(in) :: ismallL1to0component
!</input>

!<inputoutput>
  ! The filter chain.
  ! The new filter is appended at position nsize+1.
  ! The array must be large enough!
  type(t_filterChain), dimension(:), intent(inout) :: RfilterChain

  ! Current size of the filter chain.
  ! This number is increased by one and the new filter is appended.
  !
  ! For adding the first filter, this value must be set =0 !
  integer, intent(inout) :: nsize
!</inputoutput>

!</subroutine>

    if (size(RfilterChain) .lt. nsize+1) then
      call output_line ("Filter chain not large enough!", &
          OU_CLASS_ERROR, OU_MODE_STD, "filter_newFilterToL20")
      call sys_halt()
    end if
    
    ! Append the filter.
    nsize = nsize + 1
    RfilterChain(nsize)%ifilterType = FILTER_SMALLL1TO0
    RfilterChain(nsize)%ismallL1to0component = ismallL1to0component

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine filter_newFilterOneEntryZero (RfilterChain,nsize,iblock,irow)

!<description>
  ! Adds a "one entry to zero" filter to a filter chain.
!</description>

!<input>
  ! Number of the component in the vector which should be filtered.
  integer, intent(in) :: iblock
  
  ! Entry which should be set to zero
  integer, intent(in) :: irow
!</input>

!<inputoutput>
  ! The filter chain.
  ! The new filter is appended at position nsize+1.
  ! The array must be large enough!
  type(t_filterChain), dimension(:), intent(inout) :: RfilterChain

  ! Current size of the filter chain.
  ! This number is increased by one and the new filter is appended.
  !
  ! For adding the first filter, this value must be set =0 !
  integer, intent(inout) :: nsize
!</inputoutput>

!</subroutine>

    if (size(RfilterChain) .lt. nsize+1) then
      call output_line ("Filter chain not large enough!", &
          OU_CLASS_ERROR, OU_MODE_STD, "filter_newFilterOneEntryZero")
      call sys_halt()
    end if
    
    ! Append the filter.
    nsize = nsize + 1
    RfilterChain(nsize)%ifilterType = FILTER_ONEENTRY0
    RfilterChain(nsize)%iblock = iblock
    RfilterChain(nsize)%irow = irow

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine filter_newFilterOverwriteDofs (RfilterChain,nsize,&
      iblock,ndofs,IdofsOverwrite,DdofsOverwrite)

!<description>
  ! Adds a "DOF overwrite" filter to a filter chain.
!</description>

!<input>
  ! Number of the block in the vector to be filtered.
  integer, intent(in) :: iblock

  ! Number of DOFs in IdofsOverwrite
  integer, intent(in) :: ndofs

  ! List of DOFs to overwrite.
  integer, dimension(:), intent(in), target :: IdofsOverwrite

  ! OPTIONAL: List of DOFs to overwrite the DOFs with.
  ! Of not specified, teh DOFs are overwritten by zero
  REAL(DP), dimension(:), intent(in), target, optional :: DdofsOverwrite
!</input>

!<inputoutput>
  ! The filter chain.
  ! The new filter is appended at position nsize+1.
  ! The array must be large enough!
  type(t_filterChain), dimension(:), intent(inout) :: RfilterChain

  ! Current size of the filter chain.
  ! This number is increased by one and the new filter is appended.
  !
  ! For adding the first filter, this value must be set =0 !
  integer, intent(inout) :: nsize
!</inputoutput>

!</subroutine>

    if (size(RfilterChain) .lt. nsize+1) then
      call output_line ("Filter chain not large enough!", &
          OU_CLASS_ERROR, OU_MODE_STD, "filter_newFilterOneEntryZero")
      call sys_halt()
    end if
    
    ! Append the filter.
    nsize = nsize + 1
    RfilterChain(nsize)%ifilterType = FILTER_OVERWRITEDOFS
    RfilterChain(nsize)%iblock = iblock
    RfilterChain(nsize)%ndofs = ndofs
    RfilterChain(nsize)%p_IdofsOverwrite => IdofsOverwrite
    nullify(RfilterChain(nsize)%p_DdofsOverwrite)
    if (present(DdofsOverwrite)) RfilterChain(nsize)%p_DdofsOverwrite => DdofsOverwrite

  end subroutine

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
        if (associated(RfilterChain(i)%p_rdiscreteBC)) then
          call vecfil_discreteBCsol (rx,RfilterChain(i)%p_rdiscreteBC)
        else
          call vecfil_discreteBCsol (rx)
        end if

      case (FILTER_DISCBCRHSREAL)
        ! Impose Dirichlet boundary conditions into the RHS vector rx
        if (associated(RfilterChain(i)%p_rdiscreteBC)) then
          call vecfil_discreteBCrhs (rx,RfilterChain(i)%p_rdiscreteBC)
        else
          call vecfil_discreteBCrhs (rx)
        end if

      case (FILTER_DISCBCDEFREAL)
        ! Impose Dirichlet boundary conditions into the defect vector rx
        if (associated(RfilterChain(i)%p_rdiscreteBC)) then
          call vecfil_discreteBCdef (rx,RfilterChain(i)%p_rdiscreteBC)
        else
          call vecfil_discreteBCdef (rx)
        end if

      case (FILTER_DISCBCSOLFICT)
        ! Impose Dirichlet fictitious boundary conditions into the solution vector rx
        if (associated(RfilterChain(i)%p_rdiscreteFBC)) then
          call vecfil_discreteFBCsol (rx,RfilterChain(i)%p_rdiscreteFBC)
        else
          call vecfil_discreteFBCsol (rx)
        end if

      case (FILTER_DISCBCRHSFICT)
        ! Impose Dirichlet fictitious boundary conditions into the RHS vector rx
        if (associated(RfilterChain(i)%p_rdiscreteFBC)) then
          call vecfil_discreteFBCrhs (rx,RfilterChain(i)%p_rdiscreteFBC)
        else
          call vecfil_discreteFBCrhs (rx)
        end if

      case (FILTER_DISCBCDEFFICT)
        ! Impose Dirichlet fictitious boundary conditions into the defect vector rx
        if (associated(RfilterChain(i)%p_rdiscreteFBC)) then
          call vecfil_discreteFBCdef (rx,RfilterChain(i)%p_rdiscreteFBC)
        else
          call vecfil_discreteFBCdef (rx)
        end if

      case (FILTER_TOL20)
        ! Bring the subvector itoL20component of rx to the space <tex>$L^2_0$</tex>:
        call vecfil_subvectorToL20 (rx,RfilterChain(i)%itoL20component)

      case (FILTER_SMALLL1TO0)
        ! Bring the subvector itoL20component of rx to the space <tex>$L^2_0$</tex>:
        call vecfil_subvectorSmallL1To0 (rx,RfilterChain(i)%ismallL1to0component)

      case (FILTER_ONEENTRY0)
        ! Replace the entry "irow" of the subvector "iblock" of rx with zero
        call vecfil_OneEntryZero (rx,RfilterChain(i)%iblock,&
                                           RfilterChain(i)%irow)

      case (FILTER_OVERWRITEDOFS)
        ! Overwrite the specified DOFs
        if (associated(RfilterChain(i)%p_DdofsOverwrite)) then
          
          call vecfil_dofOverwrite (rx%RvectorBlock(RfilterChain(i)%iblock),&
              RfilterChain(i)%ndofs,RfilterChain(i)%p_IdofsOverwrite,&
              RfilterChain(i)%p_DdofsOverwrite)
        
        else
        
          call vecfil_dofOverwrite (rx%RvectorBlock(RfilterChain(i)%iblock),&
              RfilterChain(i)%ndofs,RfilterChain(i)%p_IdofsOverwrite)
        
        end if

      case default
        call output_line ("Unknown filter.", OU_CLASS_WARNING, OU_MODE_STD, &
                          "filter_applyFilterChainVec")
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
  ! boundary conditions that are implemented with "nonlinear filters" are
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
        if (associated(RfilterChain(i)%p_rdiscreteBC)) then
          call matfil_discreteBC (rmatrix,RfilterChain(i)%p_rdiscreteBC)
        else
          call matfil_discreteBC (rmatrix)
        end if

      case (FILTER_DISCBCSOLFICT,FILTER_DISCBCRHSFICT, &
            FILTER_DISCBCDEFFICT,FILTER_DISCBCMATFICT)
        ! Impose Dirichlet fictitious boundary conditions into the matrix rmatrix.
        ! The filter is the same for both, solution and defect filter,
        ! as the matrix modification is teh same (for now).
        if (associated(RfilterChain(i)%p_rdiscreteFBC)) then
          call matfil_discreteFBC (rmatrix,RfilterChain(i)%p_rdiscreteFBC)
        else
          call matfil_discreteFBC (rmatrix)
        end if

      case (FILTER_TOL20,FILTER_DONOTHING)
        ! These filters have no effect for matrices.

      case default
        call output_line ("Unknown filter.", OU_CLASS_WARNING, OU_MODE_STD, &
                          "filter_applyFilterChainMat")
        exit

      end select

    end do

  end subroutine

end module

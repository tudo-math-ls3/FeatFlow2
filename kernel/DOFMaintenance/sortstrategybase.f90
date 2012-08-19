!#########################################################################
!# ***********************************************************************
!# <name> sortstrategybase </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains basic structures to support the sorting of
!# matrices and vectors.
!# A sorting strategy is simply a permutation of the numbers 1..NEQ
!# how to permute a scalar vector and its inverse permutation.
!# Both permutations are usually assigned as one large array
!# to a (scalar) vector to indicate how it is resorted.
!#
!# The following routines can be found here:
!#
!# 1.) sstrat_getUnsortedPosInfo
!#     -> Returns a permutation that describes for a sorted vector
!#        the positions in the unsorted vector
!#
!# 2.) sstrat_getSortedPosInfo
!#     -> Returns a permutation that describes for an unsorted vector
!#        the positions in the sorted vector
!#
!# </purpose>
!#########################################################################

module sortstrategybase

  use fsystem
  use storage
  use genoutput

  implicit none

  private

!<constants>

!<constantblock description="Sort strategy identifiers.">

  ! No sort strategy; this must be =0!
  integer, parameter, public :: SSTRAT_UNSORTED     = 0

  ! Cuthill-McKee sort strategy
  integer, parameter, public :: SSTRAT_CM           = 1

  ! Reverse Cuthill-McKee sort strategy
  integer, parameter, public :: SSTRAT_RCM          = 2

  ! Row-wise sorting for point coordinate.
  ! (As calculated by sstrat_initXYZsorting with idirection=0.)
  ! Only for special type of discretisations (<tex>$Q_1$</tex>, <tex>$\tilde Q_1$</tex>), where
  ! the DOF`s can be identified with X/Y/Z coordinates.
  ! Coincides with the sorting strategy of FEAST for simple-type domains
  ! like the unit square.
  integer, parameter, public :: SSTRAT_XYZCOORD     = 3

  ! Column-wise sorting for point coordinate.
  ! (As calculated by sstrat_initXYZsorting with idirection=1.)
  ! Only for special type of discretisations (<tex>$Q_1$</tex>, <tex>$\tilde Q_1$</tex>), where
  ! the DOF`s can be identified with X/Y/Z coordinates.
  integer, parameter, public :: SSTRAT_ZYXCOORD     = 4

  ! General FEAST renumbering.
  ! The DOF`s are numbered rowwise, independent of the geometrical
  ! structure of the domain.
  ! Only for special type of discretisations (<tex>$Q_1$</tex>) and tensor product meshes.
  integer, parameter, public :: SSTRAT_FEAST        = 5

  ! Stochastic renumbering / Random permutation.
  ! The permutation is completely random.
  integer, parameter, public :: SSTRAT_STOCHASTIC   = 6

  ! Hierarchical renumbering: The permutation is calculated by a sequence of meshes,
  ! regularly refined.
  integer, parameter, public :: SSTRAT_HIERARCHICAL = 7

  ! Random permutation.
  ! The permutation is completely random.
  integer, parameter, public :: SSTRAT_RANDOM       = 8

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Encapsules a sorting strategy which can be used
  ! to resort matrices/vectors etc.
  type t_sortStrategy
    
    ! Type of the sorting strategy. One of the SSTRAT_xxxx constants.
    integer :: ctype = SSTRAT_UNSORTED
    
    ! Size of the permutation.
    integer :: nsize = 0
    
    ! TRUE if the structure is a copy of another sorting strategy structure.
    ! Prevents memory from being deallocated more than once.
    logical :: bisCopy = .false.
    
    ! Handle to the permutation. The calculated
    ! permutation Ipermutation given as a parameter to the routines
    ! has length 2*nsize for nsize numbers. It contains the permutation in the
    ! first N and its inverse at the second nsize entries. Here, the
    ! exact meaning of the word "permutation" is as follows:
    !
    !  Ipermutation (position in sorted vector) = position in unsorted vector.
    !  Ipermutation (nsize+position in unsorted vector) = position in sorted vector.    
    integer :: h_Ipermutation = ST_NOHANDLE
    
  end type
  
!</typeblock>

  public :: t_sortStrategy

!<typeblock>

  ! A set of sorting strategies for multiple blocks (e.g., of block vectors).
  type t_blockSortStrategy
    
    ! Number of blocks
    integer :: nblocks = 0
    
    ! List of sorting strategies for all the blocks.
    type(t_sortStrategy), dimension(:), pointer :: p_Rstrategies => null()
    
  end type
  
!</typeblock>

  public :: t_blockSortStrategy

!</types>

  public :: sstrat_getUnsortedPosInfo
  public :: sstrat_getSortedPosInfo

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sstrat_getUnsortedPosInfo (rsortStrategy,p_Ipermutation)

!<description>
  ! Returns a pointer to a permutation with
  !   Ipermutation (position in sorted vector) = position in unsorted vector.
!</description>

!<input>
  ! The sorting strategy structure to be released.
  type(t_sortStrategy), intent(in) :: rsortStrategy
!</input>

!<output>
  ! Pointer to the array.
  integer, dimension(:), pointer :: p_Ipermutation
!</output>

!</subroutine>

    if (rsortStrategy%h_Ipermutation .eq. ST_NOHANDLE) then
      call output_line ("Sort strategy not initialised.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sstrat_getUnsortedPosInfo")
      call sys_halt()
    end if
    
    ! Get the permutation
    call storage_getbase_int (rsortStrategy%h_Ipermutation,p_Ipermutation)
    
    ! Only the first part.
    p_Ipermutation => p_Ipermutation(1:rsortStrategy%nsize)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sstrat_getSortedPosInfo (rsortStrategy,p_Ipermutation)

!<description>
  ! Returns a pointer to a permutation with
  !   Ipermutation (nsize+position in unsorted vector) = position in sorted vector
!</description>

!<input>
  ! The sorting strategy structure to be released.
  type(t_sortStrategy), intent(in) :: rsortStrategy
!</input>

!<output>
  ! Pointer to the array.
  integer, dimension(:), pointer :: p_Ipermutation
!</output>

!</subroutine>

    if (rsortStrategy%h_Ipermutation .eq. ST_NOHANDLE) then
      call output_line ("Sort strategy not initialised.", &
          OU_CLASS_ERROR,OU_MODE_STD,"sstrat_getSortedPosInfo")
      call sys_halt()
    end if

    ! Get the permutation
    call storage_getbase_int (rsortStrategy%h_Ipermutation,p_Ipermutation)
    
    ! Only the second part.
    p_Ipermutation => p_Ipermutation(rsortStrategy%nsize+1:2*rsortStrategy%nsize)

  end subroutine

end module

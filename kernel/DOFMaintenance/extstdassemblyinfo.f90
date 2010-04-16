!##############################################################################
!# ****************************************************************************
!# <name> extstdassemblyinfo </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains additional assembly information that allows to
!# configure the assembly routines. The basic structure t_extScalarAssemblyInfo
!# controls e.g. the cubature rule used during the assembly.
!#
!# A t_extScalarAssemblyInfo structure is associated to a discretisation
!# structure. For the assembly of a matrix or vector with a specific cubature
!# formula e.g., the code can create a default assembly structure, modify the
!# cubature formula and assemble with this:
!#
!# <code>
!#   ! We assume for this example: 2D QUAD mesh
!#
!#   ! Get a structure and modify the cubature formulas.
!#   call easminfo_createDefInfoStructure (rdiscretisation,rassemblyInfo)
!#   rassemblyInfo%p_RinfoBlocks(:)%ccubature = CUB_G4_2D
!#
!#   ! Assemble a matrix based on this.
!#   call bilf_buildMatrixScalar2 (rform,.true.,rmatrix,&
!#       rscalarAssemblyInfo=rassemblyInfo)
!#
!#   ! Release the info structure.
!#   call easminfo_releaseInfoStructure (rassemblyInfo)
!# </code>
!#
!# Routines in this module:
!#
!# 1.) easminfo_createDefInfoStructure
!#     -> Create a default assembly information structure based on a 
!#        discretisation structure.
!#
!# 2.) easminfo_releaseInfoStructure
!#     -> Release an assembly information structure.
!#
!# </purpose>
!##############################################################################

module extstdassemblyinfo

  use fsystem
  use genoutput
  use storage
  use basicgeometry
  use boundary
  use boundaryaux
  use cubature
  use scalarpde
  use spatialdiscretisation
  
  implicit none
  
  private
  
!<types>

!<typeblock>
  
  ! Contains information that configures the assembly of matrices and vectors
  ! for a set of elements.
  type t_extScalarAssemblyInfoBlock
    
    ! Cubature rule to use.
    integer(I32) :: ccubature = 0
    
    ! Id of the element set, this assembly block refers to.
    integer :: ielementDistr = 0
    
    ! Apply the above cubature rule for all elements in the above element set.
    ! =0: Only apply for a limited set of elements specified in the element list
    ! p_IelementList.
    ! =1: Apply the above cubature rule for all elements in element set
    ! ielementDistr.
    integer :: celementListQuantifier = 1
    
    ! Number of elements in this block
    integer :: NEL = 0
    
    ! If celementListQuantifier=0, this specifies handle to a list of elements where to
    ! apply the above cubature rule. All elements shall be in the element set
    ! ielementDistr!
    integer :: h_IelementList = ST_NOHANDLE
    
  end type
  
!</typeblock>

!<typeblock>
  
  ! Contains information that configures the assembly of matrices and vectors.
  type t_extScalarAssemblyInfo
  
    ! Number of assembly information blocks in p_RinfoBlocks.
    integer :: ninfoBlockCount = 0
    
    ! A list of assembly information blocks.
    type(t_extScalarAssemblyInfoBlock), dimension(:), pointer :: p_RinfoBlocks => null()
  
  end type
  
!</typeblock>

!</types>

  public :: t_extScalarAssemblyInfoBlock
  public :: t_extScalarAssemblyInfo
  public :: easminfo_createDefInfoStructure
  public :: easminfo_releaseInfoStructure
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine easminfo_createDefInfoStructure (rdiscretisation, rassemblyInfo, &
      idepCubatureType)
  
!<description>
  ! Creates a default assembly information structure based on a discretisation.
  ! All elements in an element set are assembled with the same cubature formula.
!</description>

!<input>
  ! A discretisation structure, the assembly information structure should be
  ! associated to.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
  
  ! OPTIONAL: Compatibility flag. If present, this identifier allows to transfer
  ! a specific cubature formula from the (deprecated) discretisation structure 
  ! to the assembly structure.
  ! =0: Transfer ccubTypeBilForm.
  ! =1: Transfer ccubTypeLinForm.
  ! =2: Transfer ccubTypeEval.
  integer, intent(in), optional :: idepCubatureType
!</input>

!<output>
  ! Assembly information structure to be created.
  type(t_extScalarAssemblyInfo), intent(out) :: rassemblyInfo
!</output>

!</subroutine>
    ! local variables
    integer :: i

    ! We take as many blocks as element sets.
    rassemblyInfo%ninfoBlockCount = rdiscretisation%inumFESpaces
    allocate(rassemblyInfo%p_RinfoBlocks(rassemblyInfo%ninfoBlockCount))
    
    ! Loop through the element sets and insert the default cubature formula.
    do i = 1,rassemblyInfo%ninfoBlockCount
      rassemblyInfo%p_RinfoBlocks(i)%ielementDistr = i
      
      ! Handle all elements in the same way...
      rassemblyInfo%p_RinfoBlocks(i)%celementListQuantifier = 1
      
      ! Standard cubature formula for that element set.
      rassemblyInfo%p_RinfoBlocks(i)%ccubature = &
          spdiscr_getStdCubature(rdiscretisation%RelementDistr(i)%celement)
          
      rassemblyInfo%p_RinfoBlocks(i)%h_IelementList = ST_NOHANDLE
      
      rassemblyInfo%p_RinfoBlocks(i)%NEL = rdiscretisation%RelementDistr(i)%NEL
      
      ! Probably fetch the DEPRECATED cubature formula from the
      ! discretisation structure.
      if (present(idepCubatureType)) then
        select case (idepCubatureType)
        case (0)
          if (rdiscretisation%RelementDistr(i)%ccubTypeBilForm .ne. 0) then
            rassemblyInfo%p_RinfoBlocks(i)%ccubature = &
                rdiscretisation%RelementDistr(i)%ccubTypeBilForm
          end if

        case (1)
          if (rdiscretisation%RelementDistr(i)%ccubTypeLinForm .ne. 0) then
            rassemblyInfo%p_RinfoBlocks(i)%ccubature = &
                rdiscretisation%RelementDistr(i)%ccubTypeLinForm
          end if

        case (2)
          if (rdiscretisation%RelementDistr(i)%ccubTypeEval .ne. 0) then
            rassemblyInfo%p_RinfoBlocks(i)%ccubature = &
                rdiscretisation%RelementDistr(i)%ccubTypeEval
          end if
        end select
      end if
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine easminfo_releaseInfoStructure (rassemblyInfo)
  
!<description>
  ! Cleans up an assembly information structure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_extScalarAssemblyInfo), intent(inout) :: rassemblyInfo
!</inputoutput>

!</subroutine>
    ! local variables
    integer :: i

    ! Loop through the element sets and release memory
    do i = 1,rassemblyInfo%ninfoBlockCount
      ! Release memory if necessary
      if (rassemblyInfo%p_RinfoBlocks(i)%h_IelementList .ne. ST_NOHANDLE) then
        call storage_free(rassemblyInfo%p_RinfoBlocks(i)%h_IelementList)
      end if
    end do
    
    deallocate(rassemblyInfo%p_RinfoBlocks)

  end subroutine

end module

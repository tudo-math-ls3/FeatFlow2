!##############################################################################
!# ****************************************************************************
!# <name> assemblytemplatesoptc </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains template structures that are used to collect information
!# necessary for the assembly of matrices and vectors in the optimal control
!# problem.
!#
!# The following routines can be found here:
!#
!# 1.) astmplo_createLevelInfoHier
!#     -> Creates a template assembly hierarchy
!#
!# 2.) astmplo_releaseLevelInfoHier
!#     -> Releases a template assembly hierarchy
!#
!# </purpose>
!##############################################################################

module assemblytemplatesoptc

  use fsystem
  use storage
  use boundary
  use triangulation
  use linearsystemscalar
  
  use spatialdiscretisation
  use timediscretisation
  use timescalehierarchy

  use meshhierarchy
  use fespacehierarchybase
  use fespacehierarchy
  use spacetimehierarchy

  implicit none
  
  private
  
  public :: t_staticSpaceAsmTemplatesOptC
  public :: t_staticSpaceAsmHierarchyOptC
  public :: astmplo_createSpaceAsmHier
  public :: astmplo_releaseSpaceAsmHier
  
!<types>

!<typeblock>

  ! A type block specifying all 'static' information in space which are depending
  ! on a discretisation and a triangulation. Such static information can be
  ! precalculated and is valid until the mesh or the FE spaces change.
  ! information in this structure is devoted to the optimal control problem.
  type t_staticSpaceAsmTemplatesOptC

    ! Matrix with a precomputed EOJ stabilisation operator -- if EOJ is active.
    type(t_matrixScalar) :: rmatrixEOJ1
    type(t_matrixScalar) :: rmatrixEOJ2
    
  end type

!</typeblock>

!<typeblock>

  ! A hierarchy of t_staticSpaceAsmTemplatesOptC structures.
  type t_staticSpaceAsmHierarchyOptC

    ! Number of levels in the hierarchy.
    integer :: nlevels = 0
  
    ! The level info structures on all levels.
    type(t_staticSpaceAsmTemplatesOptC), dimension(:), pointer :: p_RasmTemplList => null()
  
  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine astmplo_createSpaceAsmHier (rhierarchy,nlevels)
  
!<description>
  ! Allocates memory for a level hierarchy consisting of nlevels levels
!</description>

!<input>
  ! Number of levels.
  integer :: nlevels
!</input>

!<output>
  ! A t_staticLevelInfoHierarchy to initialise.
  type(t_staticSpaceAsmHierarchyOptC), intent(out) :: rhierarchy
!</output>

!</subroutine>

    ! Allocate memory.
    rhierarchy%nlevels = nlevels
    allocate(rhierarchy%p_RasmTemplList(nlevels))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine astmplo_releaseSpaceAsmHier (rhierarchy)
  
!<description>
  ! Releases a levle info hierarchy.
  !
  ! WARNING: Attached matrices are NOT automatically released!
!</description>

!<inputoutput>
  ! A t_staticLevelInfoHierarchy to clean up.
  type(t_staticSpaceAsmHierarchyOptC), intent(inout) :: rhierarchy
!</inputoutput>

!</subroutine>

    ! Release memory.
    deallocate(rhierarchy%p_RasmTemplList)
    rhierarchy%nlevels = 0

  end subroutine

end module

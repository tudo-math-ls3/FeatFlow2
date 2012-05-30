!##############################################################################
!# ****************************************************************************
!# <name> assemblytemplates </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains general template structures that are used to collect information
!# necessary for the assembly of matrices and vectors.
!#
!# The following routines can be found here:
!#
!# 1.) astmpl_createLevelInfoHier
!#     -> Creates a template assembly hierarchy
!#
!# 2.) astmpl_releaseLevelInfoHier
!#     -> Releases a template assembly hierarchy
!#
!# </purpose>
!##############################################################################

module assemblytemplates

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
  
  public :: t_staticSpaceAsmTemplates
  public :: t_staticSpaceAsmHierarchy
  public :: astmpl_createSpaceAsmHier
  public :: astmpl_releaseSpaceAsmHier
  
!<types>

!<typeblock>

  ! A type block specifying all 'static' information in space which are depending
  ! on a discretisation and a triangulation. Such static information can be
  ! precalculated and is valid until the mesh or the FE spaces change.
  type t_staticSpaceAsmTemplates
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation => null()

    ! A scalar discretisation structure for the velocity space.
    type(t_spatialDiscretisation), pointer :: p_rdiscrVelocity => null()

    ! A scalar discretisation structure for the pressure space.
    type(t_spatialDiscretisation), pointer :: p_rdiscrPressure => null()

    ! A scalar discretisation structure that specifies how to generate
    ! the mass matrix in the velocity FEM space.
    ! May use a different cubature rule that the velocity discretisation
    type(t_spatialDiscretisation), pointer :: p_rdiscrMassVelocity => null()

    ! A scalar discretisation structure that specifies how to generate
    ! the mass matrix in the pressure FEM space.
    type(t_spatialDiscretisation), pointer :: p_rdiscrMassPressure => null()
    
    ! Cubature information structure for Mass-type matrices in the velocity space.
    type(t_scalarCubatureInfo) :: rcubatureInfoMassVelocity

    ! Cubature information structure for Laplace-type matrices.
    type(t_scalarCubatureInfo) :: rcubatureInfoStokes

    ! Cubature information structure for pressure/divergence matrices.
    type(t_scalarCubatureInfo) :: rcubatureInfoDiv

    ! Cubature information structure for Mass-type matrices in the velocity space.
    type(t_scalarCubatureInfo) :: rcubatureInfoMassPressure

    ! Cubature information structure for RHS vectors, momentum equation
    type(t_scalarCubatureInfo) :: rcubatureInfoRHSmomentum

    ! Cubature information structure for RHS vectors, continuity equation
    type(t_scalarCubatureInfo) :: rcubatureInfoRHScontinuity

    ! A template FEM matrix that defines the structure of Laplace/Stokes/...
    ! matrices. The matrix contains only a stucture, no content.
    type(t_matrixScalar) :: rmatrixTemplateFEM

    ! A template FEM matrix that defines the structure of the offdiagonal
    ! matrices. Needed e.g. by Newton and may have different structure
    ! than the diagonal velocity matrices in case EOJ is activated.
    ! The matrix contains only a stucture, no content.
    type(t_matrixScalar) :: rmatrixTemplateFEMOffdiag

    ! A template FEM matrix that defines the structure of the pressure
    ! matrices. The matrix contains only a stucture, no content.
    type(t_matrixScalar) :: rmatrixTemplateFEMPressure

    ! A template FEM matrix that defines the structure of gradient
    ! matrices (B1/B2). The matrix contains only a stucture, no content.
    type(t_matrixScalar) :: rmatrixTemplateGradient

    ! A template FEM matrix that defines the structure of divergence
    ! matrices (D1/D2). The matrix contains only a stucture, no content.
    type(t_matrixScalar) :: rmatrixTemplateDivergence

    ! Precalculated Laplace matrix for that specific level
    type(t_matrixScalar) :: rmatrixLaplace
    
    ! Precalculated B1-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB1

    ! Precalculated B2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB2
    
    ! Precalculated D1-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixD1

    ! Precalculated D2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixD2

    ! Precalculated D1^T-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixD1T

    ! Precalculated D2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixD2T

    ! Precalculated mass matrix for the velocity space.
    type(t_matrixScalar) :: rmatrixMassVelocity

    ! Precalculated mass matrix for the pressure space.
    type(t_matrixScalar) :: rmatrixMassPressure
    
    ! Precalculated EOJ matrix in the primal velocity space.
    ! Calculated with nu=1.
    type(t_matrixScalar) :: rmatrixEOJPrimalVel

    ! Precalculated EOJ matrix in the dual velocity space.
    ! Calculated with nu=1. May share the content with rmatrixEOJPrimal.
    type(t_matrixScalar) :: rmatrixEOJDualVel
  end type

!</typeblock>

!<typeblock>

  ! A hierarchy of t_staticSpaceAsmTemplates structures.
  type t_staticSpaceAsmHierarchy

    ! Number of levels in the hierarchy.
    integer :: nlevels = 0
  
    ! The level info structures on all levels.
    type(t_staticSpaceAsmTemplates), dimension(:), pointer :: p_RasmTemplList => null()
  
  end type

!</typeblock>


!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine astmpl_createSpaceAsmHier (rhierarchy,nlevels)
  
!<description>
  ! Allocates memory for a level hierarchy consisting of nlevels levels
!</description>

!<input>
  ! Number of levels.
  integer :: nlevels
!</input>

!<output>
  ! A t_staticLevelInfoHierarchy to initialise.
  type(t_staticSpaceAsmHierarchy), intent(out) :: rhierarchy
!</output>

!</subroutine>

    ! Allocate memory.
    rhierarchy%nlevels = nlevels
    allocate(rhierarchy%p_RasmTemplList(nlevels))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine astmpl_releaseSpaceAsmHier (rhierarchy)
  
!<description>
  ! Releases a levle info hierarchy.
  !
  ! WARNING: Attached matrices are NOT automatically released!
!</description>

!<inputoutput>
  ! A t_staticLevelInfoHierarchy to clean up.
  type(t_staticSpaceAsmHierarchy), intent(inout) :: rhierarchy
!</inputoutput>

!</subroutine>

    ! Release memory.
    deallocate(rhierarchy%p_RasmTemplList)
    rhierarchy%nlevels = 0

  end subroutine

end module

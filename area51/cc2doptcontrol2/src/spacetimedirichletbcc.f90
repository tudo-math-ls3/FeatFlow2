!##############################################################################
!# ****************************************************************************
!# <name> spacetimedirichletbcc </name>
!# ****************************************************************************
!#
!# <purpose>
!# Encapsules the Dirichlet boundary control boundary conditions in space/time.
!# </purpose>
!##############################################################################

module spacetimedirichletbcc

  use fsystem
  use storage
  use genoutput
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use paramlist
  use timestepping
  use dofmapping
  use discretebc
  use discretefbc
  use derivatives
  
  use collection
  use convection
    
  use constantsoptc
  use structuresoptc
  use user_callback

  !use spacepreconditioner
  !use spacepreconditionerinit
  use spatialbcdef
  use spacediscretisation
  use timediscretisation
  use spacetimevectors
  
  use scalarpde
  use assemblytemplates

  !use timeanalysis
  
  !use spacetimediscretisation

  implicit none
  
  private
  
!<types>
  !<typeblock>
  
  ! Contains the definition of the Neumann boundary conditions for all timesteps
  ! in a space-time discretisation.
  type t_sptiDirichletBCCBoundary
  
    ! Underlying space discretisation
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr => null()

    ! Underlying time-discretisation
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr => null()
    
    ! Number of unknowns in time
    integer :: NEQtime = 0
    
    ! Analytical definition of the Dirichlet control boundary segments
    type(t_boundaryRegionList), dimension(:), pointer :: p_RbdRegion => null()

!    ! Boundary control template matrix for the primal/dual velocity. 
!    ! This matrix contains the structure of an operator
!    ! that works on the Dirichlet control boundary. To save time, it is
!    ! saved as 'differential' matrix to the standard velocity 'template' matrix.
!    ! The matrix can be used for all timesteps as it contains a common structure.
!    type(t_matrixScalar) :: rboudaryOperatorVel
!
!    ! Boundary control template matrix for the dual pressure
!    type(t_matrixScalar) :: rboudaryOperatorGradient
    
  end type
  
  !</typeblock>
!</types>

  public :: t_sptiDirichletBCCBoundary
  public :: stdbcc_createDirichletBCCBd
  public :: stdbcc_releaseDirichletBCCBd
  public :: stdbcc_assembleDirichletBCCBd
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine stdbcc_createDirichletBCCBd (rspaceDiscr,rtimeDiscr,&
      rsptiDirichletBCC)

!<description>
  ! Creates a space-time Dirichlet boundary control boundary definition structure
  ! from a time discretisation.
!</description>

!<input>
  ! Underlying space discretisation.
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr

  ! Underlying time-discretisation.
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr
!</input>

!<output>
  ! Structure to create.
  type(t_sptiDirichletBCCBoundary), intent(out) :: rsptiDirichletBCC
!</output>

!</subroutine>

    ! Initialise the structure.
    rsptiDirichletBCC%p_rtimeDiscr => rtimeDiscr
    rsptiDirichletBCC%p_rspaceDiscr => rspaceDiscr
    rsptiDirichletBCC%NEQtime = rtimeDiscr%nintervals+1
    allocate(rsptiDirichletBCC%p_RbdRegion(rsptiDirichletBCC%NEQtime))
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine stdbcc_releaseDirichletBCCBd (rsptiDirichletBCC)

!<description>
  ! Releases a space-time Dirichlet boundary control boundary definition structure.
!</description>

!<inputoutput>
  ! Structure to release.
  type(t_sptiDirichletBCCBoundary), intent(inout) :: rsptiDirichletBCC
!</inputoutput>

!</subroutine>

    ! Initialise the structure.
    deallocate(rsptiDirichletBCC%p_RbdRegion)
    rsptiDirichletBCC%NEQtime = 0
    nullify(rsptiDirichletBCC%p_rtimeDiscr)
    nullify(rsptiDirichletBCC%p_rspaceDiscr)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine stdbcc_assembleDirichletBCCBd (roptcBDC,rsptiDirichletBCC,rglobalData)

!<description>
  ! Initialises a space-time Neumann boundary definition structure with the
  ! definition of the Neumann boundary based on the analytical
  ! definitions in the callback routines.
!</description>

!<input>
  ! Boundary conditions in the problem.
  type(t_optcBDC), intent(in) :: roptcBDC

  ! Global data, passed to callback routines
  type(t_globalData), intent(inout) :: rglobalData
!</input>

!<inputoutput>
  ! Structure to assemble.
  type(t_sptiDirichletBCCBoundary), intent(inout) :: rsptiDirichletBCC
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,j
    real(DP) :: dtimePrimal, dtimeDual, dtstep
    type(t_bilinearForm) :: rform
    type(t_bdRegionEntry), pointer :: p_rbdRegion
    
    ! Prepare a bilinear form.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff(1:3) = .true.
    rform%Dcoefficients(1:rform%itermCount) = 1.0_DP
    
    ! Loop through all timesteps
    do i=1,rsptiDirichletBCC%NEQtime
      ! Current point in time
      call tdiscr_getTimestep(rsptiDirichletBCC%p_rtimeDiscr,i-1,dtimePrimal,dtstep)
      dtimeDual = dtimePrimal - (1.0_DP-rsptiDirichletBCC%p_rtimeDiscr%dtheta)*dtstep
    
      ! Calculate the BC's.
      call sbc_releaseBoundaryList(rsptiDirichletBCC%p_RbdRegion(i))
      call sbc_assembleBDconditions (roptcBDC,dtimePrimal,dtimeDual,&
            CCSPACE_PRIMALDUAL,rglobalData,SBC_DIRICHLETBCC,&
            rsptiDirichletBCC%p_rtimediscr,rsptiDirichletBCC%p_rspacediscr,&
            rdirichletControlBoundary=rsptiDirichletBCC%p_RbdRegion(i))

    end do
    
  end subroutine

end module

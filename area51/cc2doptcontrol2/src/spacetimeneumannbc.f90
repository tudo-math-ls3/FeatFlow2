!##############################################################################
!# ****************************************************************************
!# <name> spacetimeneumannbc </name>
!# ****************************************************************************
!#
!# <purpose>
!# Encapsules the Neumann boundary conditions in space/time.
!# This type of boundary conditions is saved analytically for each time
!# discretisation as set of boundary segments. During the assembly,
!# of matrices in space, the analytical definition is evaluated and the
!# actual matrix entries are created based on this definition. The
!# definition is therefore independent of the space discretisation!
!# </purpose>
!##############################################################################

module spacetimeneumannbc

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

  !use timeanalysis
  
  !use spacetimediscretisation

  implicit none
  
  private
  
!<types>
  !<typeblock>
  
  ! Contains the definition of the Neumann boundary conditions for all timesteps
  ! in a space-time discretisation.
  type t_sptiNeumannBoundary
    
    ! Underlying space discretisation
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr => null()

    ! Underlying time-discretisation
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr => null()
    
    ! Number of unknowns in time
    integer :: NEQtime = 0
    
    ! Analytical definition of the Neumann boundary segments
    type(t_neumannBoundary), dimension(:), pointer :: p_RneumannBoundary => null()
    
  end type
  
  !</typeblock>
!</types>

  public :: t_sptiNeumannBoundary
  public :: stnm_createNeumannBoundary
  public :: stnm_releaseNeumannBoundary
  public :: stnm_assembleNeumannBoundary
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine stnm_createNeumannBoundary (rspaceDiscr,rtimeDiscr,rsptiNeumannBC)

!<description>
  ! Creates a space-time Neumann boundary definition structure
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
  type(t_sptiNeumannBoundary), intent(out) :: rsptiNeumannBC
!</output>

!</subroutine>

    ! Initialise the structure.
    rsptiNeumannBC%p_rtimeDiscr => rtimeDiscr
    rsptiNeumannBC%p_rspaceDiscr => rspaceDiscr
    rsptiNeumannBC%NEQtime = rtimeDiscr%nintervals+1
    allocate(rsptiNeumannBC%p_RneumannBoundary(rsptiNeumannBC%NEQtime))
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine stnm_releaseNeumannBoundary (rsptiNeumannBC)

!<description>
  ! Releases a space-time Neumann boundary definition structure.
!</description>

!<inputoutput>
  ! Structure to release.
  type(t_sptiNeumannBoundary), intent(inout) :: rsptiNeumannBC
!</inputoutput>

!</subroutine>

    ! Initialise the structure.
    deallocate(rsptiNeumannBC%p_RneumannBoundary)
    rsptiNeumannBC%NEQtime = 0
    nullify(rsptiNeumannBC%p_rtimeDiscr)
    nullify(rsptiNeumannBC%p_rspaceDiscr)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine stnm_assembleNeumannBoundary (roptcBDC,rsptiNeumannBC,rglobalData)

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
  ! Structure to release.
  type(t_sptiNeumannBoundary), intent(inout) :: rsptiNeumannBC
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    real(DP) :: dtimePrimal, dtimeDual, dtstep
    
    ! Loop through all timesteps
    do i=1,rsptiNeumannBC%NEQtime
      ! Current point in time
      call tdiscr_getTimestep(rsptiNeumannBC%p_rtimeDiscr,i-1,dtimePrimal,dtstep)
      dtimeDual = dtimePrimal - (1.0_DP-rsptiNeumannBC%p_rtimeDiscr%dtheta)*dtstep
    
      ! Calculate the Neumann BC's.
      call sbc_releaseNeumannBoundary(rsptiNeumannBC%p_RneumannBoundary(i))
      call sbc_assembleBDconditions (roptcBDC,dtimePrimal,dtimeDual,&
            CCSPACE_PRIMALDUAL,rglobalData,SBC_NEUMANN,&
            rsptiNeumannBC%p_rtimediscr,rsptiNeumannBC%p_rspacediscr,&
            rneumannBoundary=rsptiNeumannBC%p_RneumannBoundary(i))
    end do
  
  end subroutine

end module

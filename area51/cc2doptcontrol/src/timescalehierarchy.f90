!##############################################################################
!# ****************************************************************************
!# <name> timescalehierarchy </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module maintains time scales and hierarchies of time scales.
!#
!# A time scale hierarchy is a hierarchy of time discretisation structures
!# with increasing number of timesteps. Every refinement usually doubles the
!# number of timesteps.
!#
!# The following routines can be found here:
!#
!# </purpose>
!##############################################################################

module timescalehierarchy

  use fsystem
  use genoutput
  use boundary
  use basicgeometry
  use triangulation
  use spatialdiscretisation
  use collection
  
  use timediscretisation
  
  implicit none
  
  private
  
  public :: t_timescaleHierarchy
  public :: tmsh_createHierarchy
  public :: tmsh_releaseHierarchy
  
!<types>

!<typeblock>

  ! A structure that describes a hierarchy of time discretisations.
  type t_timescaleHierarchy
  
    ! Number of levels available in this structure.
    integer :: nlevels = 0
    
    ! Absolute start time of the simulation
    real(DP) :: dtimeInit = 0.0_DP
    
    ! Maximum time of the simulation
    real(DP) :: dtimeMax = 0.0_DP

    ! Level information.
    type(t_timeDiscretisation), dimension(:), pointer :: p_rtimeLevels => null()
    
  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine tmsh_createHierarchy (rcoarseDiscr,rtimeHierarchy,nlevels)

!<description>
  ! Creates a time scale level hierarchy for nlevels refinement levels
  ! based on a coarse time mesh rcoarseDiscr.
  ! The coarse mesh is refined ilevel-1 times.
!</description>
 
!<input>
  ! Definition of the underlying time coarse mesh.
  type(t_timeDiscretisation), intent(in) :: rcoarseDiscr
  
  ! Number of refinement levels in the hierarchy including the coarse mesh.
  integer, intent(in) :: nlevels
!</input>

!<output>
  ! A new time scale level structure.
  type(t_timescaleHierarchy), intent(out) :: rtimeHierarchy
!</output>
  
!</subroutine>

    integer :: i

    ! Basic initialisation
    rtimeHierarchy%nlevels = nlevels
    rtimeHierarchy%dtimeInit       = dtimeInit
    rtimeHierarchy%dtimeMax        = dtimeMax
    allocate(rtimeHierarchy%p_rtimeLevels(nlevels))
    
    ! Initialise the coarse mesh. Copy rcoarseDiscr.
    rtimeHierarchy%p_rtimeLevels(1) = rcoarseDiscr
        
    ! Refine.
    do i=2,nlevels
      rtimeHierarchy%p_rtimeLevels(i) = rtimeHierarchy%p_rtimeLevels(i-1)
      call tdiscr_refineRegular (rtimeHierarchy%p_rtimeLevels(i),1)
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine tmsh_releaseHierarchy (rtimeHierarchy)

!<description>
  ! Releases a time scale hierarchy.
!</description>

!<inputoutput>
  ! The hierarchy structure to release.
  type(t_timescaleHierarchy), intent(inout) :: rtimeHierarchy
!</inputoutput>
  
!</subroutine>
  
    integer :: i
    
    ! Release all levels
    do i=rtimeHierarchy%nlevels,1,-1
      call tdiscr_done (rtimeHierarchy%p_rtimeLevels(i))
    end do
    
    ! Clean up the rest
    deallocate (rtimeHierarchy%p_rtimeLevels)
    
    rtimeHierarchy%nlevels = 0
    
  end subroutine

end module

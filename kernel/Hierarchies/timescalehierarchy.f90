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
!# 1.) tmsh_createHierarchy
!#     -> Creates a time scale hierarchy
!#
!# 2.) tmsh_releaseHierarchy
!#     -> Releases a timescale hierarchy
!#
!# 3.) tmsh_printHierStatistics
!#     -> Prints statistical data of a time scale hierarchy to the terminal
!#
!# </purpose>
!##############################################################################

module timescalehierarchy

!$use omp_lib
  use fsystem
  use genoutput
  use boundary
  use basicgeometry
  use triangulation

  use timediscretisation

  implicit none

  private

  public :: t_timescaleHierarchy
  public :: tmsh_createHierarchy
  public :: tmsh_releaseHierarchy
  public :: tmsh_printHierStatistics

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
    type(t_timeDiscretisation), dimension(:), pointer :: p_RtimeLevels => null()

  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine tmsh_createHierarchy (rcoarseDiscr,rtimeHierarchy,npreref,nlevels,bprint)

!<description>
  ! Creates a time scale level hierarchy for nlevels refinement levels
  ! based on a coarse time mesh rcoarseDiscr.
  ! The coarse mesh is refined ilevel-1 times.
!</description>

!<input>
  ! Definition of the underlying time coarse mesh.
  type(t_timeDiscretisation), intent(in) :: rcoarseDiscr

  ! Number of pre-refinements to get the coarse time mesh from rcoarseDiscr.
  integer, intent(in) :: npreref

  ! Number of refinement levels in the hierarchy including the coarse mesh.
  integer, intent(in) :: nlevels

  ! OPTIONAL: Whether to print the current state of the assembly to
  ! the terminal. If set to TRUE, numbers "2 3 4..:" will be printed
  ! to the terminal for every level currently being processed.
  logical, intent(in), optional :: bprint
!</input>

!<output>
  ! A new time scale level structure.
  type(t_timescaleHierarchy), intent(out) :: rtimeHierarchy
!</output>

!</subroutine>

    integer :: i
    logical :: boutput

    boutput = .false.
    if (present(bprint)) boutput = bprint

    ! Basic initialisation
    rtimeHierarchy%nlevels = nlevels
    rtimeHierarchy%dtimeInit       = rcoarseDiscr%dtimeInit
    rtimeHierarchy%dtimeMax        = rcoarseDiscr%dtimeMax
    allocate(rtimeHierarchy%p_RtimeLevels(nlevels))

    ! Initialise the coarse mesh. Copy rcoarseDiscr.
    rtimeHierarchy%p_RtimeLevels(1) = rcoarseDiscr

    ! Pre-refine the coarse mesh.
    call tdiscr_refineRegular (rtimeHierarchy%p_RtimeLevels(1),npreref)

    ! Refine.
    do i=2,nlevels

      if (present(bprint)) then
        if (bprint) then
          ! Print current state.
          call output_line (" "//trim(sys_siL(i,10)),bnolinebreak=.true.,cdateTimeLogPolicy=OU_DTP_NONE)
        end if
      end if

      rtimeHierarchy%p_RtimeLevels(i) = rtimeHierarchy%p_RtimeLevels(i-1)
      call tdiscr_refineRegular (rtimeHierarchy%p_RtimeLevels(i),1)
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
      call tdiscr_done (rtimeHierarchy%p_RtimeLevels(i))
    end do

    ! Clean up the rest
    deallocate (rtimeHierarchy%p_RtimeLevels)

    rtimeHierarchy%nlevels = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine tmsh_printHierStatistics (rtimeHierarchy)

!<description>
  ! Writes statistics about the mesh hierarchy to the terminal.
!</description>

!<inputoutput>
  ! The hierarchy structure.
  type(t_timescaleHierarchy), intent(inout) :: rtimeHierarchy
!</inputoutput>

!</subroutine>

    integer :: i

    do i=1,rtimeHierarchy%nlevels
      call tdiscr_infoStatistics (rtimeHierarchy%p_RtimeLevels(i),i .eq. 1,i)
    end do

  end subroutine

end module

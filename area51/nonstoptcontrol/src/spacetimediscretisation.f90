!##############################################################################
!# ****************************************************************************
!# <name> spacetimediscretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the space time discretisation.
!# The structure t_ccoptSpaceTimeDiscretisation defines the general
!# discretisation of the underlying problem in space and time and can be
!# initialised and cleaned up with the routines in this module. 
!#
!# The following routines can be found here:
!#
!# 1.) sptidis_initDiscretisation
!#     -> Initialise a space-time discretisation structure.
!#
!# 2.) sptidis_doneDiscretisation
!#     -> Clean up a space-time discretisation structure.
!#
!# 3.) sptidis_infoDiscretisation
!#     -> Print information about the discretisation.
!##############################################################################

module spacetimediscretisation

  use fsystem
  use genoutput
  use spatialdiscretisation
  use timediscretisation
  use dofmapping

  implicit none
  
  private

!<types>

!<typeblock>

  ! Defines the basic shape of the supersystem which realises the coupling
  ! between all timesteps.
  ! The space-time discretisation is always a union of a discretisation in space
  ! and a discretisation in time. No memory is allocated in the structure.
  type t_spaceTimeDiscretisation
  
    ! Pointer to a space discretisation.
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr => null()
  
    ! Pointer to a time discretisation.
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr => null()
    
    ! Number of time-DOF's. For Theta-schemes, this coincides with
    ! the number of intervals, but for more complicated schemes,
    ! this may differ.
    integer :: NEQtime = 0

  end type

!</typeblock>

!</types>

  public :: t_spaceTimeDiscretisation
  public :: sptidis_initDiscretisation
  public :: sptidis_doneDiscretisation
  public :: sptidis_infoDiscretisation

contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine sptidis_initDiscretisation (rspaceDiscr,rtimeDiscr,&
      rspaceTimeDiscr)
  
!<description>
  ! Initialises a space time discretisation structure using
  ! a predefined space and time discretisation
!</description>

!<input>
  ! Space discretisation
  type(t_blockDiscretisation), target :: rspaceDiscr

  ! Time discretisation
  type(t_timeDiscretisation), target :: rtimeDiscr
!</input>

!<output>
  ! Space-time discretisation
  type(t_spaceTimeDiscretisation), intent(out) :: rspaceTimeDiscr
!</output>

!</subroutine>

    ! Save the pointers.
    rspaceTimeDiscr%p_rspaceDiscr => rspaceDiscr
    rspaceTimeDiscr%p_rtimeDiscr => rtimeDiscr

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine sptidis_doneDiscretisation (rspaceTimeDiscr)
  
!<description>
  ! Cleans up a given space time discrtisation structure.
!</description>

!<inputoutput>
  ! Supersystem-structure to be cleaned up.
  type(t_spaceTimeDiscretisation), intent(inout) :: rspaceTimeDiscr
!</inputoutput>

!</subroutine>

    nullify(rspaceTimeDiscr%p_rspaceDiscr)
    nullify(rspaceTimeDiscr%p_rtimeDiscr)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine sptidis_infoDiscretisation (rspaceTimeDiscr)
  
!<description>
  ! Prints statistical information about a time discretisation.
!</description>

!<inputoutput>
  ! Supersystem-structure to be cleaned up.
  type(t_spaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr
!</inputoutput>

!</subroutine>
    
    call output_line ('NEQ in time   : '//sys_siL(rspaceTimeDiscr%NEQTime,10))
    call output_line ('NEQ in space  : '//&
      sys_siL(dof_igetNDofGlobBlock(rspaceTimeDiscr%p_rspaceDiscr),10))

  end subroutine

end module

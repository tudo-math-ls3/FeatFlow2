!##############################################################################
!# ****************************************************************************
!# <name> flagship_signals </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for signal handling
!#
!# The following routines are available:
!#
!# 1.) flagship_SIGINT
!#     -> Performs signal handling for SIGINT
!#
!# 2.) flagship_SIGQUT
!#     -> Performs signal handling for SIGQUIT
!# </purpose>
!##############################################################################

module flagship_signals
  
#include "flagship.h"

!$ use omp_lib
  use fsystem
  use genoutput
  use signals
  
  implicit none

  private
  public :: flagship_SIGINT
  public :: flagship_SIGQUIT
  
contains
  
  !*****************************************************************************
  
!<function>
  
  function flagship_SIGINT(signum) result(sigcount)

!<description>
    ! This subroutine performs signal handling for SIGINT. In essence,
    ! it counts the number if SIGINT`s received and terminates if user
    ! sent SIGINT more than three times.
!</description>

!<input>
    integer, intent(in) :: signum
!</input>

!<result>
    ! signal
    integer :: sigcount
!</result>
!</function>

    ! local variables
    integer, save :: icount = 0
    
    sigcount = icount
    
    if (signum .eq. -1) then
      
      ! Reset counter
      icount = 0
      
    elseif (signum .eq. SIGINT) then
      
      ! Increase counter
      icount = icount+1
      if (icount .ge. 3) then
        call output_line('Simulation terminated due to user interruption (SIGINT)')
        call sys_halt()
      end if
      
    end if
  end function flagship_SIGINT
  
  !*****************************************************************************
  
!<function>

  function flagship_SIGQUIT(signum) result(sigcount)
    
!<description>
    ! This subroutine performs signal handling for SIGQUIT.
!</description>

!<input>
    integer, intent(in) :: signum
!</input>

!<result>
    ! signal
    integer :: sigcount
!</result>
!</function>

    sigcount = 0
    call output_line('Simulation terminated due to user interruption (SIGQUIT)')
    call sys_halt()
  end function flagship_SIGQUIT
  
end module flagship_signals

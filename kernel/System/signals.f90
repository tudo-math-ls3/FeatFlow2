!##############################################################################
!# ****************************************************************************
!# <name> signals </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides basic routines for signal handling in Fortran.
!# Signal handling is rather tricky in Fortran (because the function
!# that is registered as a signal handler is later called by value
!# rather than by reference), so this module provides functions which
!# are used as wrapper for the underlying C routines.
!# </purpose>
!##############################################################################
module signals

!$use omp_lib

  implicit none
  
  private
  public :: fsignal

!<constants>
!<constantblock description="Types of signals">

  integer, parameter, public :: SIGHUP    =  1   ! Hangup (POSIX).
  integer, parameter, public :: SIGINT    =  2   ! Interrupt (ANSI).
  integer, parameter, public :: SIGQUIT   =  3   ! Quit (POSIX).
  integer, parameter, public :: SIGILL    =  4   ! Illegal instruction (ANSI).
  integer, parameter, public :: SIGTRAP   =  5   ! Trace trap (POSIX).
  integer, parameter, public :: SIGABRT   =  6   ! Abort (ANSI).
  integer, parameter, public :: SIGIOT    =  6   ! IOT trap (4.2 BSD).
  integer, parameter, public :: SIGBUS    =  7   ! BUS error (4.2 BSD).
  integer, parameter, public :: SIGFPE    =  8   ! Floating-point exception (ANSI).
  integer, parameter, public :: SIGKILL   =  9   ! Kill, unblockable (POSIX).
  integer, parameter, public :: SIGUSR1   = 10   ! User-defined signal 1 (POSIX).
  integer, parameter, public :: SIGSEGV   = 11   ! Segmentation violation (ANSI).
  integer, parameter, public :: SIGUSR2   = 12   ! User-defined signal 2 (POSIX).
  integer, parameter, public :: SIGPIPE   = 13   ! Broken pipe (POSIX).
  integer, parameter, public :: SIGALRM   = 14   ! Alarm clock (POSIX).
  integer, parameter, public :: SIGTERM   = 15   ! Termination (ANSI).
  integer, parameter, public :: SIGSTKFLT = 16   ! Stack fault.
  integer, parameter, public :: SIGCLD    = 17   ! Same as SIGCHLD (System V).
  integer, parameter, public :: SIGCHLD   = 17   ! Child status has changed (POSIX).
  integer, parameter, public :: SIGCONT   = 18   ! Continue (POSIX).
  integer, parameter, public :: SIGSTOP   = 19   ! Stop, unblockable (POSIX).
  integer, parameter, public :: SIGTSTP   = 20   ! Keyboard stop (POSIX).
  integer, parameter, public :: SIGTTIN   = 21   ! Background read from tty (POSIX).
  integer, parameter, public :: SIGTTOU   = 22   ! Background write to tty (POSIX).
  integer, parameter, public :: SIGURG    = 23   ! Urgent condition on socket (4.2 BSD).
  integer, parameter, public :: SIGXCPU   = 24   ! CPU limit exceeded (4.2 BSD).
  integer, parameter, public :: SIGXFSZ   = 25   ! File size limit exceeded (4.2 BSD).
  integer, parameter, public :: SIGVTALRM = 26   ! Virtual alarm clock (4.2 BSD).
  integer, parameter, public :: SIGPROF   = 27   ! Profiling alarm clock (4.2 BSD).
  integer, parameter, public :: SIGWINCH  = 28   ! Window size change (4.3 BSD, Sun).
  integer, parameter, public :: SIGPOLL   = 29   ! Pollable event occurred (System V).
  integer, parameter, public :: SIGIO     = 29   ! I/O now possible (4.2 BSD).
  integer, parameter, public :: SIGPWR    = 30   ! Power failure restart (System V).
  integer, parameter, public :: SIGSYS    = 31   ! Bad system call.
  integer, parameter, public :: SIGUNUSED = 31

!</constantblock>
!</constants>

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

contains
  
!<subroutine>
  
  subroutine fsignal(sig,func)

!<description>
    ! This subroutine registers a signal handler with a signal. The
    ! second argument is optional. If it is omitted, then the signal
    ! number is deregistered from signal handling.
!</description>

!<input>
    ! signal number
    integer, intent(in) :: sig
    
    ! signal handler
    interface
      function func(sig)
        integer, intent(in) :: sig
        integer :: func
        
      end function func
    end interface
    optional :: func
!</input>
!</subroutine>
    
#if (!WINDOWS)
    
    ! external subroutines written in C
    external :: signal_register
    external :: signal_deregister
    
    ! (De-)register handler
    if (present(func)) then
      call signal_register(sig,func)
    else
      call signal_deregister(sig)
    end if
    
#endif
    
  end subroutine fsignal
end module signals

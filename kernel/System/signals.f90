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

#include "kernel/feat2macros.h"

!$ use omp_lib
  use genoutput

  implicit none

  private
  public :: fsignal

!<constants>

#if FEAT2_PP_OS_IS_LINUX() || FEAT2_PP_OS_IS_SOLARIS()
!<constantblock description="Types of signals available on Linux and SUN Solaris systems">
  integer, parameter, public :: SIGHUP    =  1   ! terminal line hangup
  integer, parameter, public :: SIGINT    =  2   ! interrupt program
  integer, parameter, public :: SIGQUIT   =  3   ! quit program
  integer, parameter, public :: SIGILL    =  4   ! illegal instruction
  integer, parameter, public :: SIGTRAP   =  5   ! trace trap
  integer, parameter, public :: SIGABRT   =  6   ! abort program (formerly SIGIOT)
  integer, parameter, public :: SIGBUS    =  7   ! bus error
  integer, parameter, public :: SIGFPE    =  8   ! floating-point exception
  integer, parameter, public :: SIGKILL   =  9   ! kill program
  integer, parameter, public :: SIGUSR1   = 10   ! user defined signal 1
  integer, parameter, public :: SIGSEGV   = 11   ! segmentation violation
  integer, parameter, public :: SIGUSR2   = 12   ! user defined signal 2
  integer, parameter, public :: SIGPIPE   = 13   ! write on a pipe with no reader
  integer, parameter, public :: SIGALRM   = 14   ! real-time timer expired
  integer, parameter, public :: SIGTERM   = 15   ! software termination signal
  integer, parameter, public :: SIGSTKFLT = 16   ! stack fault
  integer, parameter, public :: SIGCHLD   = 17   ! child status has changed
  integer, parameter, public :: SIGCONT   = 18   ! continue after stop
  integer, parameter, public :: SIGSTOP   = 19   ! stop (cannot be caught or ignored)
  integer, parameter, public :: SIGTSTP   = 20   ! stop signal generated from keyboard
  integer, parameter, public :: SIGTTIN   = 21   ! background read attempted from control terminal
  integer, parameter, public :: SIGTTOU   = 22   ! background write attempted to control terminal
  integer, parameter, public :: SIGURG    = 23   ! urgent condition present on socket
  integer, parameter, public :: SIGXCPU   = 24   ! cpu time limit exceeded (see setrlimit(2))
  integer, parameter, public :: SIGXFSZ   = 25   ! file size limit exceeded (see setrlimit(2))
  integer, parameter, public :: SIGVTALRM = 26   ! virtual time alarm (see setitimer(2))
  integer, parameter, public :: SIGPROF   = 27   ! profiling timer alarm (see setitimer(2))
  integer, parameter, public :: SIGWINCH  = 28   ! window size change
  integer, parameter, public :: SIGIO     = 29   ! I/O is possible on a descriptor (see fcntl(2))
  integer, parameter, public :: SIGPWR    = 30   ! power failure restart
  integer, parameter, public :: SIGSYS    = 31   ! non-existent system call invoked
!</constantblock>

#elif FEAT2_PP_OS_IS_OSX()
!<constantblock description="Types of signals available on Mac OSX systems">
  integer, parameter, public :: SIGHUP    =  1   ! terminal line hangup
  integer, parameter, public :: SIGINT    =  2   ! interrupt program
  integer, parameter, public :: SIGQUIT   =  3   ! quit program
  integer, parameter, public :: SIGILL    =  4   ! illegal instruction
  integer, parameter, public :: SIGTRAP   =  5   ! trace trap
  integer, parameter, public :: SIGABRT   =  6   ! abort program (formerly SIGIOT)
  integer, parameter, public :: SIGEMT    =  7   ! emulate instruction executed
  integer, parameter, public :: SIGFPE    =  8   ! floating-point exception
  integer, parameter, public :: SIGKILL   =  9   ! kill program
  integer, parameter, public :: SIGBUS    = 10   ! bus error
  integer, parameter, public :: SIGSEGV   = 11   ! segmentation violation
  integer, parameter, public :: SIGSYS    = 12   ! non-existent system call invoked
  integer, parameter, public :: SIGPIPE   = 13   ! write on a pipe with no reader
  integer, parameter, public :: SIGALRM   = 14   ! real-time timer expired
  integer, parameter, public :: SIGTERM   = 15   ! software termination signal
  integer, parameter, public :: SIGURG    = 16   ! urgent condition present on socket
  integer, parameter, public :: SIGSTOP   = 17   ! stop (cannot be caught or ignored)
  integer, parameter, public :: SIGTSTP   = 18   ! stop signal generated from keyboard
  integer, parameter, public :: SIGCONT   = 19   ! continue after stop
  integer, parameter, public :: SIGCHLD   = 20   ! child status has changed
  integer, parameter, public :: SIGTTIN   = 21   ! background read attempted from control terminal
  integer, parameter, public :: SIGTTOU   = 22   ! background write attempted to control terminal
  integer, parameter, public :: SIGIO     = 23   ! I/O is possible on a descriptor (see fcntl(2))
  integer, parameter, public :: SIGXCPU   = 24   ! cpu time limit exceeded (see setrlimit(2))
  integer, parameter, public :: SIGXFSZ   = 25   ! file size limit exceeded (see setrlimit(2))
  integer, parameter, public :: SIGVTALRM = 26   ! virtual time alarm (see setitimer(2))
  integer, parameter, public :: SIGPROF   = 27   ! profiling timer alarm (see setitimer(2))
  integer, parameter, public :: SIGWINCH  = 28   ! window size change
  integer, parameter, public :: SIGINFO   = 29   ! status request from keyboard
  integer, parameter, public :: SIGUSR1   = 30   ! user defined signal 1
  integer, parameter, public :: SIGUSR2   = 31   ! user defined signal 2
!</constantblock>

#elif FEAT2_PP_OS_IS_CYGWIN()
!<constantblock description="Types of signals available on Cygwin systems">
  integer, parameter, public :: SIGHUP    =  1   ! hangup
  integer, parameter, public :: SIGINT    =  2   ! interrupt
  integer, parameter, public :: SIGQUIT   =  3   ! quit
  integer, parameter, public :: SIGILL    =  4   ! illegal instruction (not reset when caught)
  integer, parameter, public :: SIGTRAP   =  5   ! trace trap (not reset when caught)
  integer, parameter, public :: SIGABRT   =  6   ! used by abort
  integer, parameter, public :: SIGEMT    =  7   ! EMT instruction
  integer, parameter, public :: SIGFPE    =  8   ! floating point exception
  integer, parameter, public :: SIGKILL   =  9   ! kill (cannot be caught or ignored)
  integer, parameter, public :: SIGBUS    = 10   ! bus error
  integer, parameter, public :: SIGSEGV   = 11   ! segmentation violation
  integer, parameter, public :: SIGSYS    = 12   ! bad argument to system call
  integer, parameter, public :: SIGPIPE   = 13   ! write on a pipe with no one to read it
  integer, parameter, public :: SIGALRM   = 14   ! alarm clock
  integer, parameter, public :: SIGTERM   = 15   ! software termination signal from kill
  integer, parameter, public :: SIGURG    = 16   ! urgent condition on IO channel
  integer, parameter, public :: SIGSTOP   = 17   ! sendable stop signal not from tty
  integer, parameter, public :: SIGTSTP   = 18   ! stop signal from tty
  integer, parameter, public :: SIGCONT   = 19   ! continue a stopped process
  integer, parameter, public :: SIGCHLD   = 20   ! to parent on child stop or exit
  integer, parameter, public :: SIGCLD    = 20   ! System V name for SIGCHLD
  integer, parameter, public :: SIGTTIN   = 21   ! to readers pgrp upon background tty read
  integer, parameter, public :: SIGTTOU   = 22   ! like TTIN for output if (tp->t_local&LTOSTOP)
  integer, parameter, public :: SIGIO     = 23   ! input/output possible
  integer, parameter, public :: SIGPOLL   = 23   ! System V name for SIGIO
  integer, parameter, public :: SIGXCPU   = 24   ! exceeded CPU time limit
  integer, parameter, public :: SIGXFSZ   = 25   ! exceeded file size limit
  integer, parameter, public :: SIGVTALRM = 26   ! virtual time alarm
  integer, parameter, public :: SIGPROF   = 27   ! profiling time alarm
  integer, parameter, public :: SIGWINCH  = 28   ! window changed
  integer, parameter, public :: SIGLOST   = 29   ! resource lost (eg, record-lock lost)
  integer, parameter, public :: SIGPWR    = 29   ! power failure
  integer, parameter, public :: SIGUSR1   = 30   ! user defined signal 1
  integer, parameter, public :: SIGUSR2   = 31   ! user defined signal 2
!</constantblock>
#endif

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

#if FEAT2_PP_OS_IS_UNIX()

    ! external subroutines written in C
    external :: signal_register
    external :: signal_deregister

    ! (De-)register handler
    if (present(func)) then
      call signal_register(sig,func)
    else
      call signal_deregister(sig)
    end if

#else

    call output_line('Signal handling not implemented for OS!',&
        OU_CLASS_WARNING,OU_MODE_STD,'fsignal')
    
#endif

  end subroutine fsignal
end module signals

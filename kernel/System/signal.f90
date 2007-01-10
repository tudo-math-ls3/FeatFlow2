!##############################################################################
!# ****************************************************************************
!# <name> signal </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides basic routines for signal handling in Fortran.
!# Signal handling is rather tricky in Fortran (because the function 
!# that is registered as a signal handler is later called by value 
!# rather than by reference), so this module provides functions which
!# are used as wrapper for the underlying C routines.
!##############################################################################
MODULE signal
  IMPLICIT NONE
  
  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<constants>
!<constantblock description="Types of signals">

  INTEGER, PARAMETER :: SIGHUP    =  1   ! Hangup (POSIX).
  INTEGER, PARAMETER :: SIGINT    =  2   ! Interrupt (ANSI).
  INTEGER, PARAMETER :: SIGQUIT   =  3   ! Quit (POSIX).
  INTEGER, PARAMETER :: SIGILL    =  4   ! Illegal instruction (ANSI).
  INTEGER, PARAMETER :: SIGTRAP   =  5   ! Trace trap (POSIX).
  INTEGER, PARAMETER :: SIGABRT   =  6   ! Abort (ANSI).
  INTEGER, PARAMETER :: SIGIOT    =  6   ! IOT trap (4.2 BSD).
  INTEGER, PARAMETER :: SIGBUS    =  7   ! BUS error (4.2 BSD).
  INTEGER, PARAMETER :: SIGFPE    =  8   ! Floating-point exception (ANSI).
  INTEGER, PARAMETER :: SIGKILL   =  9   ! Kill, unblockable (POSIX).
  INTEGER, PARAMETER :: SIGUSR1   = 10   ! User-defined signal 1 (POSIX).
  INTEGER, PARAMETER :: SIGSEGV   = 11   ! Segmentation violation (ANSI).
  INTEGER, PARAMETER :: SIGUSR2   = 12   ! User-defined signal 2 (POSIX).
  INTEGER, PARAMETER :: SIGPIPE   = 13   ! Broken pipe (POSIX).
  INTEGER, PARAMETER :: SIGALRM   = 14   ! Alarm clock (POSIX).
  INTEGER, PARAMETER :: SIGTERM   = 15   ! Termination (ANSI).
  INTEGER, PARAMETER :: SIGSTKFLT = 16   ! Stack fault.
  INTEGER, PARAMETER :: SIGCLD    = 17   ! Same as SIGCHLD (System V).
  INTEGER, PARAMETER :: SIGCHLD   = 17   ! Child status has changed (POSIX).
  INTEGER, PARAMETER :: SIGCONT   = 18   ! Continue (POSIX).
  INTEGER, PARAMETER :: SIGSTOP   = 19   ! Stop, unblockable (POSIX).
  INTEGER, PARAMETER :: SIGTSTP   = 20   ! Keyboard stop (POSIX).
  INTEGER, PARAMETER :: SIGTTIN   = 21   ! Background read from tty (POSIX).
  INTEGER, PARAMETER :: SIGTTOU   = 22   ! Background write to tty (POSIX).
  INTEGER, PARAMETER :: SIGURG    = 23   ! Urgent condition on socket (4.2 BSD).
  INTEGER, PARAMETER :: SIGXCPU   = 24   ! CPU limit exceeded (4.2 BSD).
  INTEGER, PARAMETER :: SIGXFSZ   = 25   ! File size limit exceeded (4.2 BSD).
  INTEGER, PARAMETER :: SIGVTALRM = 26   ! Virtual alarm clock (4.2 BSD).
  INTEGER, PARAMETER :: SIGPROF   = 27   ! Profiling alarm clock (4.2 BSD).
  INTEGER, PARAMETER :: SIGWINCH  = 28   ! Window size change (4.3 BSD, Sun).
  INTEGER, PARAMETER :: SIGPOLL   = 29   ! Pollable event occurred (System V).
  INTEGER, PARAMETER :: SIGIO     = 29   ! I/O now possible (4.2 BSD).
  INTEGER, PARAMETER :: SIGPWR    = 30   ! Power failure restart (System V).
  INTEGER, PARAMETER :: SIGSYS    = 31   ! Bad system call.
  INTEGER, PARAMETER :: SIGUNUSED = 31

!</constantblock>
!</constants>  

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

CONTAINS
  
!<subroutine>
  
  SUBROUTINE fsignal(sig,func)

!<description>
    ! This subroutine registers a signal handler with a signal. The
    ! second argument is optional. If it is omitted, then the signal
    ! number is deregistered from signal handling.
!</description>

!<input>
    ! signal number
    INTEGER, INTENT(IN) :: sig
    
    ! signal handler
    INTERFACE 
      FUNCTION func(sig)
        INTEGER, INTENT(IN) :: sig
        INTEGER :: func
        
      END FUNCTION func
    END INTERFACE
    OPTIONAL :: func
!</input>
!</subroutine>
    
    ! external subroutines written in C
    EXTERNAL :: signal_register
    EXTERNAL :: signal_deregister
    
    ! (De-)register handler
    IF (PRESENT(func)) THEN
      CALL signal_register(sig,func)
    ELSE
      CALL signal_deregister(sig)
    END IF
    
  END SUBROUTINE fsignal
END MODULE signal

!<description>
! This program is able to calculate
! ||u1 - u2||_\Omega - but u1 and u2 are allowed to live
! on different meshes. In fact, all we need is that \Omega1
! and \Omega2 need to have some overlap. To be more precise,
! u1 can be calculated with elementtype 1
! and u2 can be calculated with elementtype 2
! and u1 can live on mesh 1 and u2 on mesh 2 - all that
! needs to be fulfilled is that mesh 1 and mesh 2 need to
! have an overlap - but if this would not be the case,
! this calculation would be pretty much useless.
! This means that this program is able to cover the most
! general case.
! However, we only compare the velocities, so if you
! need to compare the pressure this does not work at the moment.
!
! If you are familiar with cc2d you will notice that
! many routines look somehow similar - this is because
! I started from a cc2d-based problem, so I copied and
! slightly modified some parts from this application
! (e.g. grid-generation, init of discretisation, ...)
!
! I would like to thank M. Koester and P. Zajac for their
! support when I wrote this program.
!
! M. Schuh (mschuh@mathematik.tu-dortmund.de)
!
!</description>

program ExtFEcomparer

  use ExtFEcomparer_main
  use fparser

  implicit none

  ! For the logfiles
  character(LEN=SYS_STRLEN) :: slogfile,serrorfile,sbenchlogfile

  ! The very first thing in every application:
  ! Initialise system-wide settings:
  call sys_init()

  ! General output init - output to the terminal
  call output_init ()

  ! Get the names of the logfiles
  call ExtFEcomparer_getLogFiles(slogfile,serrorfile,sbenchlogfile)

  ! release the old output to the terminal
  call output_done()

  ! Initialise log file for output.
  call output_init (slogfile,serrorfile,sbenchlogfile)

  ! Initialise the storage management:
  call storage_init(999, 100)

  ! Initialise the parser
  call fparser_init ()

  call output_lbrk ()
  call output_line ("Start to compare the solutions")
  call output_separator (OU_SEP_MINUS)

  ! Now we can start
  call start_ExtFEcomparer_main()

  ! Release the parser
  call fparser_done ()

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display "Handles in use=0" and "Memory in use=0"!
  call output_lbrk ()
  call storage_info(.true.)

  ! release the output
  call output_done()

  ! Clean up the storage management, finish
  call storage_done()

end program

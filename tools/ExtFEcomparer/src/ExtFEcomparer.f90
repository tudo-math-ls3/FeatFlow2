!<description>
! Version 2.1 - 3.11.2014
!
! To sum up the features and present the changelog, we need a notation:
! f_k = Component k of function f, for example if we solve 2D-Navier-Stokes,
! f_1 = Velocity in x direction, f_2 = Velocity in y direction, f_3 = pressure
! f is created by a representing vector of coefficents (i.e. the solution vector
! of the Navier-Stokes-Equations from cc2d), the domain (represented by
! *.tri and *.prm files) and the finite elements that were used, represented
! by a name or (in case of Navier-Stokes) an identifier for the element pair.
! All this has to be present to actually use the tool (for 2 functions).
! Hint: If you want to analyze only 1 function, enter it 2 times as input:
! One time as f, one time as g
!
! Summary of features:
!
! 1) Calculations for 1D and 2D
! 1.1) Calculate ||f_k||_L2(\Omega)
! 1.2) Calculate ||g_k||_L2(\Omega)
! 1.3) Calculate ||f_k - g_l||_L2(\Omega)
! 1.4) Calculate (d^i/dx^i) f_k(x) |i|=0,1 i: multiindex
! 1.5) Calculate (d^j/dx^j) g_k(x) |j|=0,1 j: multiindex
! 1.6) Calculate (d^i/dx^i)f_k(x) - (d^j/dx^j)g_l(x), |i|=0,1; |j| = 0,1 i,j: multiindex
!
! 2) Output
! 2.1) Export the vector of coefficents in an other format
! 2.2) UCD-Output of the mesh
! 2.3) UCD-Output of the input functions
!
! I would like to thank M. Koester and P. Zajac for their
! support when I started writing this program.
!
! M. Schuh (mschuh@mathematik.tu-dortmund.de)
!
!</description>

program ExtFEcomparer

  use ExtFEcomparer_main
  use ExtFEcomparer_parameters

  use fsystem
  use genoutput
  use storage
  use fparser

  implicit none

  ! For the logfiles
  character(LEN=ExtFE_STRLEN) :: slogfile,serrorfile,sbenchlogfile

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
  call storage_init(200, 100)

  ! Initialise the parser
  call fparser_init ()

  call output_lbrk ()
  call output_line ("Starting the Extended Finite Element comparer")

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

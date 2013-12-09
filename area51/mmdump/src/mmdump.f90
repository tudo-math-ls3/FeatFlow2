program mmdump

  use mmdump2d

  implicit none

  type(t_parlist) :: rparam

  ! Initialise system-wide settings:
  call sys_init()

  ! Initialise the output system.
  call output_init ('./out/output.log')

  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Read in parameter list
  call parlst_init(rparam)
  call parlst_readfromfile(rparam, './data/mmdump.dat')

  ! Call the problem to solve. Poisson 1D method 1 - simple:
  call mmdump_2d(rparam)

  ! Release parameter list
  call parlst_done (rparam)

  ! Print out heap statistics
  call output_lbrk ()
  call storage_info(.true.)

  ! Clean up the storage management, finish
  call storage_done()

end program

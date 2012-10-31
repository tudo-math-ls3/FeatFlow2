program fme

  use fsystem
  use genoutput
  use paramlist
  
  use fme2d

  implicit none

  type(t_parlist) :: rparam
  integer :: nDim

  ! Initialise system-wide settings:
  call system_init()

  ! Initialise the output system.
  call output_init()

  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Read in parameter list
  call parlst_init(rparam)
  call parlst_readfromfile(rparam, './data/master.dat')

  ! fetch parameters
  call parlst_getvalue_int(rparam, '', 'NDIM', nDim, 0)
  
  ! call exporter
  select case(nDim)
  case (2)
    call fme_2d(rparam)
  end select

  ! Release parameter list
  call parlst_done (rparam)

  ! Print out heap statistics
  call output_lbrk ()
  call storage_info(.true.)

  ! Clean up the storage management, finish
  call storage_done()

end program

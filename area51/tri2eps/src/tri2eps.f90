program mmdump

  use fsystem
  use genoutput
  use storage
  use paramlist
  use boundary
  use triangulation

  implicit none

  type(t_parlist) :: rparam
  type(t_boundary) :: rbnd
  type(t_triangulation) :: rtria
  character(len=256) :: smesh
  integer :: nLevel
  real(DP), dimension(2) :: Dbox
  real(DP) :: dlineWidth

  ! Initialise system-wide settings:
  call system_init()

  ! Initialise the output system.
  call output_init()

  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Read in parameter list
  call parlst_init(rparam)
  call parlst_readfromfile(rparam, './data/tri2eps.dat')

  ! fetch parameters
  call parlst_getvalue_string(rparam, '', 'SMESH', smesh, 'quad')
  call parlst_getvalue_int(rparam, '', 'NLEVEL', nLevel, 0)
  call parlst_getvalue_double(rparam, '', 'DSIZEX', Dbox(1), 100.0_DP)
  call parlst_getvalue_double(rparam, '', 'DSIZEY', Dbox(2), 100.0_DP)
  call parlst_getvalue_double(rparam, '', 'DLINEWIDTH', dlineWidth, 0.1_DP)

  ! read boundary and mesh
  call output_line('Reading mesh...')
  call boundary_read_prm(rbnd, './mesh/' // trim(smesh) // '.prm')
  call tria_readTriFile2D (rtria, './mesh/' // trim(smesh) // '.tri', rbnd)

  ! refine mesh
  if(nLevel .gt. 0) then
    call output_line('Refining mesh...')
    call tria_quickRefine2LevelOrdering(nLevel, rtria, rbnd)
  end if
  call tria_initStandardMeshFromRaw (rtria,rbnd)

  ! export to eps
  call output_line('Exporting mesh...')
  call tria_exportPostScript(rtria, &
      './eps/' // trim(smesh) // '_lvl' // trim(sys_sil(nLevel,4)) // '.eps', &
      Dbox, dlineWidth = dlineWidth)

  ! cleanup
  call tria_done(rtria)
  call boundary_release(rbnd)

  ! Release parameter list
  call parlst_done (rparam)

  ! Print out heap statistics
  call output_lbrk ()
  call storage_info(.true.)

  ! Clean up the storage management, finish
  call storage_done()

end program

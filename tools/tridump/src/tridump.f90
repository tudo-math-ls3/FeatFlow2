program tridump

  use fparser
  use fsystem
  use genoutput
  use storage
  use boundary
  use triangulation
  use meshgeneration
  use meshmodification
  use ucd

  implicit none

  type(t_boundary) :: rbnd
  type(t_triangulation) :: rtria
  character(len=256) :: sarg, smesh, sname, spredir
  integer :: i, iarg, nref, ndim
  real(DP) :: ddist, dline, dalpha
  real(DP), dimension(2) :: Dbox
  real(DP), dimension(4) :: Drect
  integer, dimension(2) :: Isize
  type(t_ucdExport) :: rucd
  logical :: bgen, bvtk, bgmv, beps, bbnd

  ! Initialise system-wide settings:
  call sys_init()

  ! Initialise the output system.
  call output_init()

  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)
  
  ! Initialise function parser
  call fparser_init()

  if (.not. sys_getenv_string("PREDIR", spredir)) spredir = "./mesh"

  ! print help
  if(sys_ncommandLineArgs() .lt. 1) then
    call output_lbrk()
    call output_line("USAGE: tridump <options>")
    call output_lbrk()
    call output_line("Valid options:")
    call output_line("-read2d <mesh>              Read 2D mesh from <mesh>.tri/prm")
    call output_line("-read3d <mesh>              Read 3D mesh from <mesh>.tri")
    call output_line("-ref <n>                    Refine the mesh <n> times")
    call output_line("-dist <x>                   Disturb mesh by <x>")
    call output_line("-vtk                        Write VTK file")
    call output_line("-gmv                        Write GMV file")
    call output_line("-eps                        Write EPS file")
    call output_line("-box <x> <y>                Set bounding box for EPS export")
    call output_line("-line <w>                   Set line width for EPS export")
    call output_line("-gen-rect                   Generate a rectangular mesh")
    !call output_line("-gen-wave <alpha>           Geneate an wave mesh")
    call output_line("-rect <x0> <x1> <y0> <y1>   Set domain borders for mesh generation")
    call output_line("-size <nx> <ny>             Set sizes for mesh generation")
    call output_lbrk()
    call output_lbrk()
    !                0         1         2        3          4         5         6         7         8
    !                 123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    call output_line("Examples:")
    call output_line("---------")
    call output_line("tridump -vtk -read2d bench1 -ref 4 -dist 0.05")
    call output_line("-> Reads in the 'bench1' mesh, refines it 4 times, disturbs")
    call output_line("   it by a factor of 0.05 and exports it in VTK format.")
    call output_lbrk()
    stop
  end if
  
  ! initialise standard parameters
  nref = 0
  ndim = 0
  bbnd = .false.
  bgen = .false.
  bvtk = .false.
  bgmv = .false.
  beps = .false.
  dalpha = 0.0_DP
  ddist = 0.0_DP
  Dbox = (/ 100.0_DP, 100.0_DP /)
  Drect = (/ 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP /)
  Isize = (/ 1, 1 /)
  dline = 0.1_DP
  smesh = ''
  
  ! loop over all arguments
  iarg = 1
  do while(iarg .lt. sys_ncommandLineArgs())
    ! fetch argument string
    call sys_getcommandLineArg(iarg, sarg)
    iarg = iarg + 1

    ! check argument
    if(sarg .eq. '-vtk') then
      bvtk = .true.
    else if(sarg .eq. '-gmv') then
      bgmv = .true.
    else if(sarg .eq. '-eps') then
      beps = .true.
    else if(sarg .eq. '-ref') then
      call sys_getcommandLineArg(iarg, sarg)
      iarg = iarg + 1
      read (sarg, *) nref
    else if(sarg .eq. '-dist') then
      call sys_getcommandLineArg(iarg, sarg)
      iarg = iarg + 1
      read (sarg, *) ddist
    else if(sarg .eq. '-box') then
      call sys_getcommandLineArg(iarg, sarg)
      iarg = iarg + 1
      read (sarg, *) Dbox(1)
      call sys_getcommandLineArg(iarg, sarg)
      iarg = iarg + 1
      read (sarg, *) Dbox(2)
    else if(sarg .eq. '-line') then
      call sys_getcommandLineArg(iarg, sarg)
      iarg = iarg + 1
      read (sarg, *) dline
    else if(sarg .eq. '-read2d') then
      call sys_getcommandLineArg(iarg, smesh)
      iarg = iarg + 1
      ndim = 2
    else if(sarg .eq. '-read3d') then
      call sys_getcommandLineArg(iarg, smesh)
      iarg = iarg + 1
      ndim = 3
    else if(sarg .eq. '-gen-rect') then
      bgen = .true.
    !else if(sarg .eq. '-gen-wave') then
    !  bgen = .true.
    !  call sys_getcommandLineArg(iarg, sarg)
    !  iarg = iarg + 1
    !  read (sarg, *) dalpha
    else if(sarg .eq. '-size') then
      do i=1,2
        call sys_getcommandLineArg(iarg, sarg)
        iarg = iarg + 1
        read (sarg, *) Isize(i)
      end do
    else if(sarg .eq. '-rect') then
      do i = 1, 4
        call sys_getcommandLineArg(iarg, sarg)
        iarg = iarg + 1
        read (sarg, *) Drect(i)
      end do
    else
      ! unknown parameter
      call output_line("ERROR: unknown parameter '" // trim(sarg) // "'")
      stop
    end if
  end do

  if(.not. bgen) then
    ! read boundary and mesh
    if(ndim .eq. 2) then
      call output_line("Reading mesh from '"//trim(spredir)//"/"//trim(smesh)//".tri'...")
      call boundary_read_prm(rbnd, trim(spredir)//"/"// trim(smesh) // '.prm')
      call tria_readTriFile2D (rtria, trim(spredir)//"/" // trim(smesh) // '.tri', rbnd)
      bbnd = .true.
    else if(ndim .eq. 3) then
      call output_line("Reading mesh from '"//trim(spredir)//"/"//trim(smesh)//".tri'...")
      call tria_readTriFile3D (rtria, trim(spredir)//"/"// trim(smesh) // '.tri')
      bbnd = .false.
    else
      call output_line("ERROR: no input mesh specified")
      stop
    end if
  else
    !if(dalpha .ne. 0.0_DP) then
    !  ! mesh generation
    !  smesh = 'gen-wave'
    !  call meshgen_wave2DQuadMesh (rtria,Drect(1),Drect(2),Drect(3),Drect(4),Isize(1),Isize(2),dalpha)
    !else
      ! mesh generation
      smesh = 'gen-rect'
      call meshgen_rectangular2DQuadMesh (rtria,Drect(1),Drect(2),Drect(3),Drect(4),Isize(1),Isize(2))
    !end if
  end if

  ! refine mesh
  if(nref .gt. 0) then
    call output_line('Refining mesh '//trim(sys_sil(nref,7))//' times...')
    if(bbnd) then
      call tria_quickRefine2LevelOrdering(nref, rtria, rbnd)
    else
      call tria_quickRefine2LevelOrdering(nref, rtria)
    end if
  end if
  
  ! initialise standard mesh
  if(bbnd) then
    call tria_initStandardMeshFromRaw (rtria, rbnd)
  else
    call tria_initStandardMeshFromRaw (rtria)
  end if
  
  ! distort mesh if desired
  if(ddist .ne. 0.0_DP) then
    call output_line("Disturbing mesh by " // trim(sys_sdl(ddist,7)) // "...")
    call meshmod_disturbMesh(rtria, ddist)
  end if

  ! build filename
  if(nref .gt. 0) then
    sname = './out/' // trim(smesh) // '_lvl' // trim(sys_sil(nref,4))
  else
    sname = './out/' // trim(smesh)
  end if
  
  
  ! export to EPS
  if(beps) then
    call output_line("Exporting mesh to '" // trim(sname) // ".eps'")
    call tria_exportPostScript(rtria, trim(sname) // '.eps', Dbox, dlineWidth = dline)
  end if

  ! export to VTK
  if(bvtk) then
    call output_line("Exporting mesh to '" // trim(sname) // ".vtk'")
    call ucd_startVTK (rucd, UCD_FLAG_STANDARD, rtria, trim(sname) //".vtk")
    call ucd_write (rucd)
    call ucd_release (rucd)
  end if

  ! export to GMV
  if(bgmv) then
    call output_line("Exporting mesh to '" // trim(sname) // ".gmv'")
    call ucd_startGMV (rucd, UCD_FLAG_STANDARD, rtria, trim(sname) //".gmv")
    call ucd_write (rucd)
    call ucd_release (rucd)
  end if
  
  ! cleanup
  call tria_done(rtria)
  if(bbnd) &
    call boundary_release(rbnd)

  ! Clean up the storage management, finish
  call storage_done()

end program

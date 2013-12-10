program meshadapt

  use boundary
  use fsystem
  use genoutput
  use hadaptaux
  use hadaptivity
  use io
  use linearsystemscalar
  use storage
  use triangulation

  implicit none
  
  type(t_vectorScalar) :: rindicator
  type(t_boundary) :: rbnd
  type(t_hadapt) :: rhadapt
  type(t_triangulation) :: rtria
  real(DP), dimension(:), pointer :: p_Dindicator
  integer, dimension(13) :: Istat1, Istat2
  character(len=256) :: sarg, smesh, serror, sname, sdata
  integer :: i, iarg, idelay, iel, ilinelen, ios, iunit, nrefmax, ndim
  real(DP) :: dreftol, dcrstol
  logical :: bbnd, bdaemon

  ! Initialise system-wide settings:
  call sys_init()

  ! Initialise the output system.
  call output_init()

  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)
  
  ! Print help
  if(sys_ncommandLineArgs() .lt. 1) then
    call output_lbrk()
    call output_line("USAGE: tridump <options>")
    call output_lbrk()
    call output_line("Valid options:")
    call output_line("-read1d <mesh>              Read 1D mesh from <mesh>.tri")
    call output_line("-read2d <mesh>              Read 2D mesh from <mesh>.tri/prm")
    call output_line("-read3d <mesh>              Read 3D mesh from <mesh>.tri")
    call output_line("-error <file>               Read elementwise error distribution from <file>")
    call output_line("-refmax <n>                 Maximum number of refinement levels")
    call output_line("-reftol <d>                 Tolerance for refinement")
    call output_line("-crstol <d>                 Tolerance for recoarsening")
    call output_line("-daemon <i>                 Run as daemon with <i> seconds delay (stop by Ctrl-D)")
    call output_lbrk()
    call sys_halt()
  end if

  ! Initialise standard parameters
  nrefmax = 0
  ndim = 0
  bbnd = .false.
  bdaemon = .false.
  dreftol = 0.0_DP
  dcrstol = 0.0_DP
  smesh = ''
  serror = ''

  ! Loop over all arguments
  iarg = 1
  do while(iarg .lt. sys_ncommandLineArgs())
    ! fetch argument string
    call sys_getcommandLineArg(iarg, sarg)
    iarg = iarg + 1
    
    ! check argument
    if(sarg .eq. '-read1d') then
      call sys_getcommandLineArg(iarg, smesh)
      iarg = iarg + 1
      ndim = 1
    else if(sarg .eq. '-read2d') then
      call sys_getcommandLineArg(iarg, smesh)
      iarg = iarg + 1
      ndim = 2
    else if(sarg .eq. '-read3d') then
      call sys_getcommandLineArg(iarg, smesh)
      iarg = iarg + 1
      ndim = 3
    else if(sarg .eq. '-error') then
      call sys_getcommandLineArg(iarg, serror)
      iarg = iarg + 1
    else if(sarg .eq. '-refmax') then
      call sys_getcommandLineArg(iarg, sarg)
      iarg = iarg + 1
      read (sarg, *) nrefmax
    else if(sarg .eq. '-reftol') then
      call sys_getcommandLineArg(iarg, sarg)
      iarg = iarg + 1
      read (sarg, *) dreftol
    else if(sarg .eq. '-crstol') then
      call sys_getcommandLineArg(iarg, sarg)
      iarg = iarg + 1
      read (sarg, *) dcrstol
    else if(sarg .eq. '-daemon') then
      call sys_getcommandLineArg(iarg, sarg)
      iarg = iarg + 1
      read (sarg, *) idelay
      bdaemon = .true.
    else
      ! unknown parameter
      call output_line("ERROR: unknown parameter '"//trim(sarg)//"'")
      call sys_halt()
    end if
  end do

  ! Read boundary and mesh
  if(ndim .eq. 1) then
    call output_line("Reading mesh from '"//trim(smesh)//".tri'...")
    call tria_readTriFile1D (rtria, trim(smesh)//'.tri')
    bbnd = .false.
  else if(ndim .eq. 2) then
    inquire(file=trim(smesh)//'.prm', exist=bbnd)
    if (bbnd) then
      call output_line("Reading mesh from '"//trim(smesh)//".tri/prm'...")
      call boundary_read_prm(rbnd, trim(smesh)//'.prm')
      call tria_readTriFile2D (rtria, trim(smesh)//'.tri', rbnd)
    else
      call output_line("Reading mesh from '"//trim(smesh)//".tri'...")
      call tria_readTriFile2D (rtria, trim(smesh)//'.tri')
    end if
  else if(ndim .eq. 3) then
    call output_line("Reading mesh from '"//trim(smesh)//".tri'...")
    call tria_readTriFile3D (rtria, trim(smesh)//'.tri')
    bbnd = .false.
  else
    call output_line("ERROR: no input mesh specified")
    call sys_halt()
  end if
  
  ! Initialise adaptation structure from triangulation
  call hadapt_initFromTriangulation(rhadapt, rtria)

  ! Set some parameters manually
  rhadapt%nsubdividemax        = nrefmax
  rhadapt%iadaptationStrategy  = HADAPT_REDGREEN
  rhadapt%drefinementTolerance = dreftol
  rhadapt%dcoarseningTolerance = dcrstol
  rhadapt%iSpec = ior(rhadapt%iSpec, HADAPT_HAS_PARAMETERS)

  ! Are we in daemon mode?
  if (bdaemon) then
    ! Get status of indicator field file
    call stat(trim(serror), Istat1)
  end if

  ! Infinite loop for potential daemon mode
  daemon: do

    ! Read indicator vector
    call output_line("Reading indicator field from '"//trim(serror)//"'...")
    call io_openFileForReading(trim(serror), iunit, .true.)
    
    ! Read first line from file
    read(iunit, fmt=*) iel
    if (iel .eq. 0) exit daemon
    if (iel .ne. rtria%NEL) then
      call output_line ("Mismatch in number of elements!", &
                        OU_CLASS_ERROR,OU_MODE_STD,"meshadapt")
      call sys_halt()
    end if
    
    ! Create indicator
    call lsyssc_createVector(rindicator, rtria%NEL, .true., ST_DOUBLE)
    call lsyssc_getbase_double(rindicator, p_Dindicator)
    
    do iel=1,rtria%NEL
      call io_readlinefromfile(iunit, sdata, ilinelen, ios)
      p_Dindicator(iel) = sys_str2Double(sdata, "(F20.10)")
    end do
    close(iunit)
    
    ! Perform mesh adaptation
    call hadapt_performAdaptation(rhadapt, rindicator)
    
    ! Release indicator
    call lsyssc_releaseVector(rindicator)
    
    ! Update triangulation structure
    call hadapt_generateRawMesh(rhadapt, rtria)
  
    ! Initialise standard mesh
    if(bbnd) then
      call tria_initStandardMeshFromRaw (rtria, rbnd)
    else
      call tria_initStandardMeshFromRaw (rtria)
    end if
    
    ! Export triangulation structure
    call output_line("Exporting triangulation to '"//trim(smesh)//"_ref.tri'...")
    if (bbnd) then
      call tria_exportTriFile(rtria, trim(smesh)//'_ref.tri', TRI_FMT_STANDARD)
    else
      call tria_exportTriFile(rtria, trim(smesh)//'_ref.tri', TRI_FMT_NOPARAMETRISATION)
    end if
    
    ! Are we in daemon mode?
    if (.not.bdaemon) exit daemon
    
    delay: do
      ! Get status of indicator field file
      call stat(trim(serror), Istat2)
      if (Istat1(10) .ne. Istat2(10)) then
        Istat1 = Istat2
        exit delay
      else
        call sleep(idelay)
      end if
    end do delay
    
  end do daemon
  
  ! Clean up
  call hadapt_releaseAdaptation(rhadapt)
  call tria_done(rtria)
  if(bbnd) call boundary_release(rbnd)
  
  ! Clean up the storage management, finish
  call storage_done()
  
end program meshadapt

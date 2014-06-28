program meshadapt

  use boundary
  use fsystem
  use genoutput
  use hadaptaux
  use hadaptivity
  use io
  use linearsystemscalar
  use signals
  use storage
  use triangulation

  use meshadaptbase

  implicit none
  
  character(len=256) :: sarg,sname
  integer :: i,iarg,nrefmax,ndim
  real(DP) :: dreftol,dcrstol
  logical :: bdaemon

  ! Initialise system-wide settings:
  call sys_init()

  ! Initialise the output system.
  call output_init()

  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)
  
  ! Print help
  if(sys_ncommandLineArgs() .lt. 1) then
    call output_lbrk()
    call output_line("USAGE: meshadapt <options>")
    call output_lbrk()
    call output_line("Valid options:")
    call output_line("-read1d <mesh>              Read 1D mesh from <mesh>.tri")
    call output_line("-read2d <mesh>              Read 2D mesh from <mesh>.tri/prm")
    call output_line("-read3d <mesh>              Read 3D mesh from <mesh>.tri")
    call output_line("-error <file>               Read elementwise error distribution from <file>")
    call output_line("-refmax <n>                 Maximum number of refinement levels")
    call output_line("-reftol <d>                 Tolerance for refinement")
    call output_line("-crstol <d>                 Tolerance for recoarsening")
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
  do while(iarg .le. sys_ncommandLineArgs())
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
      bdaemon = .true.
    else
      ! unknown parameter
      call output_line("ERROR: unknown parameter '"//trim(sarg)//"'")
      call sys_halt()
    end if
  end do

  ! Read boundary and mesh
  select case(ndim)
  case (1)
    call output_line("Reading mesh from './"//trim(smesh)//".tri'...")
    call tria_readTriFile1D (rtria, './'//trim(smesh)//'.tri')
    bbnd = .false.
  case (2)
    inquire(file='./'//trim(smesh)//'.prm', exist=bbnd)
    if (bbnd) then
      call output_line("Reading mesh from './"//trim(smesh)//".tri/prm'...")
      call boundary_read_prm(rbnd, './'//trim(smesh)//'.prm')
      call tria_readTriFile2D (rtria, './'//trim(smesh)//'.tri', rbnd)
    else
      call output_line("Reading mesh from './"//trim(smesh)//".tri'...")
      call tria_readTriFile2D (rtria, './'//trim(smesh)//'.tri')
    end if
  case (3)
    call output_line("Reading mesh from './"//trim(smesh)//".tri'...")
    call tria_readTriFile3D (rtria, './'//trim(smesh)//'.tri')
    bbnd = .false.
  case default
    call output_line("ERROR: no input mesh specified")
    call sys_halt()
  end select
  
  ! Initialise standard mesh
  if(bbnd) then
    call tria_initStandardMeshFromRaw (rtria, rbnd)
  else
    call tria_initStandardMeshFromRaw (rtria)
  end if

  ! Set some parameters manually
  rhadapt%nsubdividemax        = nrefmax
  rhadapt%iadaptationStrategy  = HADAPT_REDGREEN
  rhadapt%drefinementTolerance = dreftol
  rhadapt%dcoarseningTolerance = dcrstol
  rhadapt%iSpec                = ior(rhadapt%iSpec, HADAPT_HAS_PARAMETERS)

  ! Initialise adaptation structure from triangulation
  call hadapt_initFromTriangulation(rhadapt, rtria)

  ! Are we in daemon mode?
  if (bdaemon) then
    call fsignal(SIGINT, perform_adaptation)
    call fsignal(SIGQUIT, perform_adaptation)
    
!    daemon: do
!      call sleep(10)
!    end do daemon
  else
    i = perform_adaptation(SIGINT)
    i = perform_adaptation(SIGQUIT)
  end if

end program meshadapt

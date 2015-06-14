program ferecovery

  use boundary
  use cubature
  use element
  use fsystem
  use genoutput
  use hadaptaux
  use hadaptivity
  use io
  use linearsystemblock
  use linearsystemscalar
  use pprocgradients
  use spatialdiscretisation
  use storage
  use triangulation
  use vectorio

  implicit none

  type(t_blockDiscretisation) :: rblockDiscrSolution, rblockDiscrGradient
  type(t_boundary) :: rbnd
  type(t_triangulation) :: rtria
  type(t_vectorBlock) :: rsolution, rgradient
  real(DP), dimension(:), pointer :: p_Ddata
  character(len=256) :: sarg,scubature,selement,sgrad,smesh,sname,ssol,stoken
  integer ::  i,iarg,iel,ilinelen,ios,itoken,iunit,ntoken,ndim
  integer(I32), dimension(:), allocatable :: Celement
  logical :: bbnd

  ! Initialise system-wide settings:
  call sys_init()

  ! Initialise the output system.
  call output_init()

  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)

  ! Print help
  if(sys_ncommandLineArgs() .lt. 1) then
    call output_lbrk()
    call output_line("USAGE: ferecovery <options>")
    call output_lbrk()
    call output_line("Valid options:")
    call output_line("-read1d <mesh>              Read 1D mesh from <mesh>.tri")
    call output_line("-read2d <mesh>              Read 2D mesh from <mesh>.tri/prm")
    call output_line("-read3d <mesh>              Read 3D mesh from <mesh>.tri")
    call output_line("-sol <file>                 Read nodal solution values from <file>")
    call output_line("-grad <file>                Write nodal gradient values to <file>")
    call output_line("-elem <elementID>           Specify type(s) of element by <elementID>")
    call output_lbrk()
    call sys_halt()
  end if

  ! Initialise standard parameters
  ndim = 0
  bbnd = .false.
  smesh = ''
  sgrad = ''
  ssol = ''

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
    else if(sarg .eq. '-sol') then
      call sys_getcommandLineArg(iarg, ssol)
      iarg = iarg + 1
    else if(sarg .eq. '-grad') then
      call sys_getcommandLineArg(iarg, sgrad)
      iarg = iarg + 1
    else if(sarg .eq. '-elem') then
      call sys_getcommandLineArg(iarg, selement)
      call sys_countTokens(selement, ntoken, ',')
      allocate(Celement(ntoken))
      i = 1
      do itoken = 1,ntoken
        call sys_getNextToken(selement, stoken, i, ',')
        Celement(itoken) = elem_igetID(trim(stoken))
      end do
      iarg = iarg + 1
    else
      ! unknown parameter
      call output_line("ERROR: unknown parameter '" //trim(sarg)//"'")
      call sys_halt()
    end if
  end do

  ! Read boundary and mesh
  select case(ndim)
  case (1)
    call output_line("Reading mesh from '"//"./"//trim(smesh)//".tri'...")
    call tria_readTriFile1D (rtria, "/."//trim(smesh)//'.tri')
    bbnd = .false.
  case (2)
    inquire(file="./"//trim(smesh)//'.prm', exist=bbnd)
    if (bbnd) then
      call output_line("Reading mesh from '"//"./"//trim(smesh)//".tri/prm'...")
      call boundary_read_prm(rbnd, "./"//trim(smesh)//'.prm')
      call tria_readTriFile2D (rtria, "./"//trim(smesh)//'.tri', rbnd)
    else
      call output_line("Reading mesh from '"//"./"//trim(smesh)//".tri'...")
      call tria_readTriFile2D (rtria, "./"//trim(smesh)//'.tri')
    end if
  case (3)
    call output_line("Reading mesh from '"//"./"//trim(smesh)//".tri'...")
    call tria_readTriFile3D (rtria, "./"//trim(smesh)//'.tri')
    bbnd = .false.
  case default
    call output_line("ERROR: no input mesh specified")
    call sys_halt()
  end select

  ! Initialise block and spatial discretisations
  if (bbnd) then
    call spdiscr_initBlockDiscr(rblockDiscrSolution, 1, rtria, rbnd)
    call spdiscr_initBlockDiscr(rblockDiscrGradient, ndim, rtria, rbnd)
    select case(ndim)
    case (1)
      call spdiscr_initDiscr_simple(rblockDiscrSolution%RspatialDiscr(1),&
          Celement(1), rtria, rbnd)
      call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(1),&
          Celement(1), rtria, rbnd)
    case (2)
      if (size(Celement,1) .eq. 1) then
        call spdiscr_initDiscr_simple(rblockDiscrSolution%RspatialDiscr(1),&
            Celement(1), rtria, rbnd)
        call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(1),&
            Celement(1), rtria, rbnd)
        call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(2),&
            Celement(1), rtria, rbnd)
      else
        call spdiscr_initDiscr_simple(rblockDiscrSolution%RspatialDiscr(1),&
            Celement(1), Celement(2), rtria, rbnd)
        call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(1),&
            Celement(1), Celement(2), rtria, rbnd)
        call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(2),&
            Celement(1), Celement(2), rtria, rbnd)
      end if
    case (3)
      call spdiscr_initDiscr_simple(rblockDiscrSolution%RspatialDiscr(1),&
          Celement(1), rtria, rbnd)
      call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(1),&
          Celement(1), rtria, rbnd)
      call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(2),&
          Celement(1), rtria, rbnd)
      call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(3),&
          Celement(1), rtria, rbnd)
    end select
  else
    call spdiscr_initBlockDiscr(rblockDiscrSolution, 1, rtria)
    call spdiscr_initBlockDiscr(rblockDiscrGradient, ndim, rtria)
    select case(ndim)
    case (1)
      call spdiscr_initDiscr_simple(rblockDiscrSolution%RspatialDiscr(1),&
          Celement(1), rtria)
      call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(1),&
          Celement(1), rtria)
    case (2)
      if (size(Celement,1) .eq. 1) then
        call spdiscr_initDiscr_simple(rblockDiscrSolution%RspatialDiscr(1),&
            Celement(1), rtria)
        call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(1),&
            Celement(1), rtria)
        call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(2),&
            Celement(1), rtria)
      else
        call spdiscr_initDiscr_simple(rblockDiscrSolution%RspatialDiscr(1),&
            Celement(1), Celement(2), rtria)
        call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(1),&
            Celement(1), Celement(2), rtria)
        call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(2),&
            Celement(1), Celement(2), rtria)
      end if
    case (3)
      call spdiscr_initDiscr_simple(rblockDiscrSolution%RspatialDiscr(1),&
          Celement(1), rtria)
      call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(1),&
          Celement(1), rtria)
      call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(2),&
          Celement(1), rtria)
      call spdiscr_initDiscr_simple(rblockDiscrGradient%RspatialDiscr(3),&
          Celement(1), rtria)
    end select
  end if

  ! Read solution from file
  call lsysbl_createVector(rblockDiscrSolution, rsolution, .false.)
  call output_line("Reading solution from '"//"./"//trim(ssol)//"'...")
  call vecio_readBlockVectorHR(rsolution, sname, .true., 0, "./"//trim(ssol), .true.)

  ! Recover gradient by superconvergent patch recovery
  call lsysbl_createVector(rblockDiscrGradient, rgradient, .false.)
  call ppgrd_calcGradient(rsolution%RvectorBlock(1), rgradient,&
      PPGRD_ZZTECHNIQUE, PPGRD_NODEPATCH)

  ! Write gradient to file
  call output_line("Writing gradient to '"//"./"//trim(sgrad)//"'...")
  call vecio_writeBlockVectorHR(rgradient, 'grad', .true., 0, "./"//trim(sgrad), '(E20.10)')

  ! Clean up
  deallocate(Celement)
  call lsysbl_releaseVector(rsolution)
  call lsysbl_releaseVector(rgradient)
  call spdiscr_releaseBlockDiscr(rblockDiscrSolution)
  call spdiscr_releaseBlockDiscr(rblockDiscrGradient)
  call tria_done(rtria)
  if(bbnd) call boundary_release(rbnd)

  ! Clean up the storage management, finish
  call storage_done()

end program ferecovery

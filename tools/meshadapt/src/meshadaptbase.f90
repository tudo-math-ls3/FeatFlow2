module meshadaptbase

  use boundary
  use fsystem
  use genoutput
  use hadaptaux
  use hadaptivity
  use io
  use linearalgebra
  use linearsystemscalar
  use signals
  use storage
  use triangulation

  implicit none
  
!<types>
!<typeblock>

  ! A mesh adaptation structure that can be passed between routines.

  type t_meshAdapt

    ! Boundary
    type(t_boundary), pointer :: rboundary => null()

    ! Triangulation
    type(t_triangulation) :: rtriangulation

    ! Adaptation structure
    type(t_hadapt) :: rhadapt    

    ! Maximum refinement level
    integer :: nrefmax = 0

    ! Refinement tolerance
    real(DP) :: dreftol = 0.0_DP

    ! Coarsening tolerance
    real(DP) :: dcrstol = 0.0_DP

  end type t_meshAdapt

!</typeblock>
!</types>

  private
  
  public :: t_meshAdapt

  public :: madapt_alloc
  public :: madapt_dealloc
  public :: madapt_init
  public :: madapt_done
  public :: madapt_step
  public :: madapt_signalhandler

  public :: madapt_getNEL
  public :: madapt_getNVT
  public :: madapt_getNDIM
  public :: madapt_getVertexCoords
  public :: madapt_getNeighboursAtElement
  public :: madapt_getVerticesAtElement

  interface madapt_step
    module procedure madapt_step_direct
    module procedure madapt_step_dble1
    module procedure madapt_step_dble2
    module procedure madapt_step_sngl1
    module procedure madapt_step_sngl2
    module procedure madapt_step_fromfile
  end interface madapt_step

contains
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine madapt_alloc(rmeshAdapt)

!<description>
  ! This subroutine allocates memory for a mesh adaptation structure
!</description>

!<input>
    ! Pointer to mesh adaptation structure
    type(t_meshAdapt), intent(inout), pointer :: rmeshAdapt
!</input>

!</subroutine>

    nullify(rmeshAdapt)
    allocate(rmeshAdapt)

  end subroutine madapt_alloc

  ! ***************************************************************************

!<subroutine>
  
  subroutine madapt_dealloc(rmeshAdapt)

!<description>
  ! This subroutine deallocates memory of a mesh adaptation structure
!</description>

!<inputoutput>
    ! Pointer to mesh adaptation structure
    type(t_meshAdapt), intent(inout), pointer :: rmeshAdapt
!</inputoutput>

!</subroutine>

    if (associated(rmeshAdapt)) deallocate(rmeshAdapt)

  end subroutine madapt_dealloc

  ! ***************************************************************************

!<subroutine>
  
  subroutine madapt_init(rmeshAdapt,ndim,smesh)

!<description>
  ! This subroutine initialises the mesh adaptation structure
!</description>

!<input>
    ! Number of spatial dimension
    integer, intent(in) :: ndim

    ! Name if the mesh file
    character(len=*), intent(in) :: smesh
!</input>

!<output>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(out) :: rmeshAdapt
!</output>

!</subroutine>

    ! local variables
    logical bhasPRMfile

    ! Check if mesh file is not empty
    if (trim(adjustl(smesh)) .eq. '') then
      call output_line("No input mesh specified!",&
          OU_CLASS_ERROR,OU_MODE_STD,"madapt_init")
      call sys_halt()
    end if

    ! Read TRI and PRM file (if it exists)
    select case(ndim)
    case (1)
      call output_line("Reading mesh from './"//trim(smesh)//".tri'...")
      call tria_readTriFile1D(rmeshAdapt%rtriangulation,&
          './'//trim(smesh)//'.tri')
    case (2)
      inquire(file='./'//trim(smesh)//'.prm', exist=bhasPRMfile)
      if (bhasPRMfile) then
        allocate(rmeshAdapt%rboundary)
        call output_line("Reading boundary from './"//trim(smesh)//".prm'...")
        call boundary_read_prm(rmeshAdapt%rboundary, './'//trim(smesh)//'.prm')
        call output_line("Reading mesh from './"//trim(smesh)//".tri'...")
        call tria_readTriFile2D(rmeshAdapt%rtriangulation,&
            './'//trim(smesh)//'.tri', rmeshAdapt%rboundary)
      else
        call output_line("Reading mesh from './"//trim(smesh)//".tri'...")
        call tria_readTriFile2D(rmeshAdapt%rtriangulation,&
            './'//trim(smesh)//'.tri')
      end if
    case (3)
      call output_line("Reading mesh from './"//trim(smesh)//".tri'...")
      call tria_readTriFile3D (rmeshAdapt%rtriangulation,&
          './'//trim(smesh)//'.tri')
    case default
      call output_line("Invalid spatial dimension!",&
          OU_CLASS_ERROR,OU_MODE_STD,"madapt_init")
      call sys_halt()
    end select

    ! Initialise standard mesh
    if(associated(rmeshAdapt%rboundary)) then
      call output_line("Generating standard mesh from raw triangulation...")
      call tria_initStandardMeshFromRaw(rmeshAdapt%rtriangulation,&
          rmeshAdapt%rboundary)
    else
      call tria_initStandardMeshFromRaw(rmeshAdapt%rtriangulation)
    end if

    ! Set some parameters manually
    rmeshAdapt%rhadapt%iadaptationStrategy  = HADAPT_REDGREEN
    rmeshAdapt%rhadapt%iSpec                = ior(rmeshAdapt%rhadapt%iSpec,&
                                                    HADAPT_HAS_PARAMETERS)

    ! Initialise adaptation structure from triangulation
    call output_line("Initialising adaptation data structure...")
    call hadapt_initFromTriangulation(rmeshAdapt%rhadapt, rmeshAdapt%rtriangulation)
    
  end subroutine madapt_init
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine madapt_done(rmeshAdapt)

!<description>
  ! This subroutine finalises the mesh adaptation structure
!</description>

!<inputoutput>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(inout) :: rmeshAdapt
!</inputoutput>

!</subroutine>

    call hadapt_releaseAdaptation(rmeshAdapt%rhadapt)
    call tria_done(rmeshAdapt%rtriangulation)
    if(associated(rmeshAdapt%rboundary)) call boundary_release(rmeshAdapt%rboundary)

  end subroutine madapt_done

  ! ***************************************************************************

!<subroutine>

  subroutine madapt_step_direct(rmeshAdapt,rindicator,nrefmax,dreftol,dcrstol)

!<description>
  ! This subroutine performs a single mesh adaptation step based on
  ! the refinement and re-coarsening tolerances, the maximum
  ! refinement level and the element-wise error indicator given.
  !
  ! The error indicator is required as scalar vectors.
!</description>

!<input>
    ! Element-wise error indicator vector
    type(t_vectorScalar), intent(in) :: rindicator

    ! Maximum refinement level
    integer, intent(in) :: nrefmax

    ! Refinement tolerance
    real(DP), intent(in) :: dreftol

    ! Re-coarsening tolerance
    real(DP), intent(in) :: dcrstol
!</input>

!<inputoutput>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(inout) :: rmeshAdapt
!</inputoutput>

!</subroutine>

    if (nrefmax .lt. rmeshAdapt%rhadapt%nsubdividemax) then
      if (hadapt_getSubdivisonLevel(rmeshAdapt%rhadapt) .le. nrefmax) then
        rmeshAdapt%rhadapt%nsubdividemax = nrefmax
      else
        call output_line("Maximum refinement level cannot be smaller than "//&
            "refinement level in current mesh!",&
            OU_CLASS_ERROR,OU_MODE_STD,"madapt_step_direct")
        call sys_halt()
      end if
    else
      rmeshAdapt%rhadapt%nsubdividemax = nrefmax
    end if

    ! Set parameters
    rmeshAdapt%rhadapt%drefinementTolerance = dreftol
    rmeshAdapt%rhadapt%dcoarseningTolerance = dcrstol

    ! Perform mesh adaptation
    call hadapt_refreshAdaptation(rmeshAdapt%rhadapt, rmeshAdapt%rtriangulation)
    call hadapt_performAdaptation(rmeshAdapt%rhadapt, rindicator)
    
    ! Update triangulation structure
    call hadapt_generateRawMesh(rmeshAdapt%rhadapt, rmeshAdapt%rtriangulation)
    
    ! Initialise standard mesh
    if(associated(rmeshAdapt%rboundary)) then
      call tria_initStandardMeshFromRaw(rmeshAdapt%rtriangulation,&
          rmeshAdapt%rboundary)
    else
      call tria_initStandardMeshFromRaw(rmeshAdapt%rtriangulation)
    end if

  end subroutine madapt_step_direct

  ! ***************************************************************************

!<subroutine>

  subroutine madapt_step_dble1(rmeshAdapt,Dindicator,nrefmax,dreftol,dcrstol)

!<description>
  ! This subroutine performs a single mesh adaptation step based on
  ! the refinement and re-coarsening tolerances, the maximum
  ! refinement level and the element-wise error indicator given.
  !
  ! The error indicator is required as double precision array
!</description>

!<input>
    ! Element-wise error indicator vector
    real(DP), dimension(:), intent(in) :: Dindicator

    ! Maximum refinement level
    integer, intent(in) :: nrefmax

    ! Refinement tolerance
    real(DP), intent(in) :: dreftol

    ! Re-coarsening tolerance
    real(DP), intent(in) :: dcrstol
!</input>

!<inputoutput>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(inout) :: rmeshAdapt
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorScalar) :: rindicator
    real(DP), dimension(:), pointer :: p_Dindicator
        
    ! Create scalar indicator vector
    call lsyssc_createVector(rindicator, size(Dindicator), .true., ST_DOUBLE)
    call lsyssc_getbase_double(rindicator, p_Dindicator)
    call lalg_copyVector(Dindicator, p_Dindicator)

    ! Perform mesh adaptation step
    call madapt_step(rmeshAdapt,rindicator,nrefmax,dreftol,dcrstol)

    ! Release scalar indicator vector
    call lsyssc_releaseVector(rindicator)

  end subroutine madapt_step_dble1

  ! ***************************************************************************

!<subroutine>

  subroutine madapt_step_dble2(rmeshAdapt,nel,Dindicator,nrefmax,dreftol,dcrstol)

!<description>
  ! This subroutine performs a single mesh adaptation step based on
  ! the refinement and re-coarsening tolerances, the maximum
  ! refinement level and the element-wise error indicator given.
  !
  ! The error indicator is required as double precision array
!</description>

!<input>
    ! Number of elements
    integer, intent(in) :: nel

    ! Element-wise error indicator vector
    real(DP), dimension(nel), intent(in) :: Dindicator

    ! Maximum refinement level
    integer, intent(in) :: nrefmax

    ! Refinement tolerance
    real(DP), intent(in) :: dreftol

    ! Re-coarsening tolerance
    real(DP), intent(in) :: dcrstol
!</input>

!<inputoutput>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(inout) :: rmeshAdapt
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorScalar) :: rindicator
    real(DP), dimension(:), pointer :: p_Dindicator

    ! Create scalar indicator vector
    call lsyssc_createVector(rindicator, size(Dindicator), .true., ST_DOUBLE)
    call lsyssc_getbase_double(rindicator, p_Dindicator)
    call lalg_copyVector(Dindicator, p_Dindicator)

    ! Perform mesh adaptation step
    call madapt_step(rmeshAdapt,rindicator,nrefmax,dreftol,dcrstol)

    ! Release scalar indicator vector
    call lsyssc_releaseVector(rindicator)

  end subroutine madapt_step_dble2

  ! ***************************************************************************

!<subroutine>

  subroutine madapt_step_sngl1(rmeshAdapt,Findicator,nrefmax,dreftol,dcrstol)

!<description>
  ! This subroutine performs a single mesh adaptation step based on
  ! the refinement and re-coarsening tolerances, the maximum
  ! refinement level and the element-wise error indicator given.
  !
  ! The error indicator is required as single precision array
!</description>

!<input>
    ! Element-wise error indicator vector
    real(SP), dimension(:), intent(in) :: Findicator

    ! Maximum refinement level
    integer, intent(in) :: nrefmax

    ! Refinement tolerance
    real(DP), intent(in) :: dreftol

    ! Re-coarsening tolerance
    real(DP), intent(in) :: dcrstol
!</input>

!<inputoutput>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(inout) :: rmeshAdapt
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorScalar) :: rindicator
    real(DP), dimension(:), pointer :: p_Dindicator
    
    ! Create scalar indicator vector
    call lsyssc_createVector(rindicator, size(Findicator), .true., ST_DOUBLE)
    call lsyssc_getbase_double(rindicator, p_Dindicator)
    call lalg_copyVector(Findicator, p_Dindicator)

    ! Perform mesh adaptation step
    call madapt_step(rmeshAdapt,rindicator,nrefmax,dreftol,dcrstol)

    ! Release scalar indicator vector
    call lsyssc_releaseVector(rindicator)

  end subroutine madapt_step_sngl1

  ! ***************************************************************************

!<subroutine>

  subroutine madapt_step_sngl2(rmeshAdapt,nel,Findicator,nrefmax,dreftol,dcrstol)

!<description>
  ! This subroutine performs a single mesh adaptation step based on
  ! the refinement and re-coarsening tolerances, the maximum
  ! refinement level and the element-wise error indicator given.
  !
  ! The error indicator is required as single precision array
!</description>

!<input>
    ! Number of elements
    integer, intent(in) :: nel

    ! Element-wise error indicator vector
    real(SP), dimension(nel), intent(in) :: Findicator

    ! Maximum refinement level
    integer, intent(in) :: nrefmax

    ! Refinement tolerance
    real(DP), intent(in) :: dreftol

    ! Re-coarsening tolerance
    real(DP), intent(in) :: dcrstol
!</input>

!<inputoutput>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(inout) :: rmeshAdapt
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorScalar) :: rindicator
    real(DP), dimension(:), pointer :: p_Dindicator
    
    ! Create scalar indicator vector
    call lsyssc_createVector(rindicator, nel, .true., ST_DOUBLE)
    call lsyssc_getbase_double(rindicator, p_Dindicator)
    call lalg_copyVector(Findicator, p_Dindicator)

    ! Perform mesh adaptation step
    call madapt_step(rmeshAdapt,rindicator,nrefmax,dreftol,dcrstol)

    ! Release scalar indicator vector
    call lsyssc_releaseVector(rindicator)

  end subroutine madapt_step_sngl2

  ! ***************************************************************************

!<subroutine>

  subroutine madapt_step_fromfile(rmeshAdapt,sindicator,nrefmax,dreftol,dcrstol)

!<description>
  ! This subroutine performs a single mesh adaptation step based on
  ! the refinement and re-coarsening tolerances, the maximum
  ! refinement level and the element-wise error indicator given.
!</description>

!<input>
    ! Name of the element-wise error indicator file
    character(len=*), intent(in) :: sindicator

    ! Maximum refinement level
    integer, intent(in) :: nrefmax

    ! Refinement tolerance
    real(DP), intent(in) :: dreftol

    ! Re-coarsening tolerance
    real(DP), intent(in) :: dcrstol
!</input>

!<inputoutput>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(inout) :: rmeshAdapt
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorScalar) :: rindicator
    real(DP), dimension(:), pointer :: p_Dindicator
    character(len=SYS_STRLEN) :: sdata
    integer :: iel,ilinelen,ios,iunit,nel

    ! Check if mesh file is not empty
    if (trim(adjustl(sindicator)) .eq. '') then
      call output_line("No indicator specified!",&
          OU_CLASS_ERROR,OU_MODE_STD,"madapt_step_fromfile")
      call sys_halt()
    end if

    ! Read indicator vector
    call output_line("Reading indicator field from './"//trim(sindicator)//"'...")
    call io_openFileForReading('./'//trim(sindicator), iunit, .true.)

    ! Read first line from file
    read(iunit, fmt=*) nel

    ! Create scalar indicator vector
    call lsyssc_createVector(rindicator, nel, .true., ST_DOUBLE)
    call lsyssc_getbase_double(rindicator, p_Dindicator)
    
    do iel=1,nel
      call io_readlinefromfile(iunit, sdata, ilinelen, ios)
      p_Dindicator(iel) = sys_stringToDouble(sdata, "(F20.10)")
    end do
    close(iunit)
    
    ! Perform mesh adaptation step
    call madapt_step(rmeshAdapt,rindicator,nrefmax,dreftol,dcrstol)
    
    ! Release scalar indicator vector
    call lsyssc_releaseVector(rindicator)
    
  end subroutine madapt_step_fromfile

  ! ***************************************************************************

!<function>

  function madapt_signalhandler(isignum) result(iresult)
    
!<description>
  ! This function performs mesh adaptation based on signals received
  ! from an outer daemon process. The following signals are supported:
  !
  ! SIGUSR1: Import mesh from file
  ! SIGUSR2: Export mesh to file
  ! SIGINT : Perform single mesh adaptation step
  ! SIGQUIT: Finalise mesh adaptation and quit
  ! SIGHUP : Finalise mesh adaptation and hangup
  ! SIGTERM: Finalise mesh adaptation and terminate
!</description>

!<input>
    ! Signal number
    integer, intent(in) :: isignum
!</input>

!<result>
    ! Result
    integer :: iresult
!</result>

!</function>

    ! persistent variables
    type(t_meshAdapt), save :: rmeshAdapt
    character(len=SYS_STRLEN), save :: smesh='', sindicator=''
    real(DP), save :: dreftol=0.0_DP, dcrstol=0.0_DP
    integer, save :: nrefmax=1

    ! local variables
    character(len=SYS_STRLEN) :: sarg,sname
    integer :: iarg,istep,ndim
    logical :: bdaemon=.false.
    
    select case(isignum)
    case (SIGUSR1) !----- Import mesh from file --------------------------------

      ! Get arguments from command line
      if(sys_ncommandLineArgs() .lt. 1) then
        call output_lbrk()
        call output_line("USAGE: meshadapt <options>")
        call output_lbrk()
        call output_line("Valid options:")
        call output_line("-read1d <mesh>     Read 1D mesh from <mesh>.tri")
        call output_line("-read2d <mesh>     Read 2D mesh from <mesh>.tri/prm")
        call output_line("-read3d <mesh>     Read 3D mesh from <mesh>.tri")
        call output_line("-indicator <file>  Read elementwise refinement indicator from <file>")
        call output_line("-refmax <n>        Maximum number of refinement levels")
        call output_line("-reftol <d>        Tolerance for refinement")
        call output_line("-crstol <d>        Tolerance for recoarsening")
        call output_lbrk()
        call sys_halt()
      end if
      
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
        else if(sarg .eq. '-indicator') then
          call sys_getcommandLineArg(iarg, sindicator)
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
          call output_line("Unknown parameter '"//trim(sarg)//"'!",&
              OU_CLASS_ERROR,OU_MODE_STD,"madapt_signalhandler")
          call sys_halt()
        end if
      end do
      
      ! Initialise mesh adaptation structure
      call madapt_init(rmeshAdapt,ndim,smesh)
      
      iresult = merge(0,1,bdaemon)

    case (SIGUSR2) !----- Export mesh to file ----------------------------------

      call output_line("Exporting triangulation to './"//&
          trim(smesh)//"_ref.tri'...")
      if (associated(rmeshAdapt%rboundary)) then
        call tria_exportTriFile(rmeshAdapt%rtriangulation, './'//&
            trim(smesh)//'_ref.tri', TRI_FMT_STANDARD)
      else
        call tria_exportTriFile(rmeshAdapt%rtriangulation, './'//&
            trim(smesh)//'_ref.tri', TRI_FMT_NOPARAMETRISATION)
      end if

      iresult = 0

    case (SIGINT) !----- Perform mesh adaptation step --------------------------

      call madapt_step(rmeshAdapt,sindicator,nrefmax,dreftol,dcrstol)
            
      iresult = 0

    case (SIGQUIT,&
          SIGHUP,&
          SIGTERM) !----- Finalise mesh adaptation -----------------------------

      call madapt_done(rmeshAdapt)
      
      ! Clean up the storage management
      call storage_done()
      
      ! Finish
      stop
      
    case default
      
      iresult = -1

    end select

  end function madapt_signalhandler

  ! ***************************************************************************

!<function>

  function madapt_getNEL(rmeshAdapt) result(NEL)

!<description>
    ! Returns the number of elements
!</description>

!<input>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(in) :: rmeshAdapt
!</input>

!<result>
    ! Number of elements
    integer :: NEL
!</result>

!</function>

    NEL = rmeshAdapt%rtriangulation%NEL

  end function madapt_getNEL

  ! ***************************************************************************

!<function>

  function madapt_getNVT(rmeshAdapt) result(NVT)

!<description>
    ! Returns the number of vertices
!</description>

!<input>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(in) :: rmeshAdapt
!</input>

!<result>
    ! Number of elements
    integer :: NVT
!</result>

!</function>

    NVT = rmeshAdapt%rtriangulation%NVT

  end function madapt_getNVT

  ! ***************************************************************************

!<function>

  function madapt_getNDIM(rmeshAdapt) result(NDIM)

!<description>
    ! Returns the number of spatial dimensions
!</description>

!<input>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(in) :: rmeshAdapt
!</input>

!<result>
    ! Number of spatial dimensions
    integer :: NDIM
!</result>

!</function>

    NDIM = rmeshAdapt%rtriangulation%NDIM

  end function madapt_getNDIM

  ! ***************************************************************************

!<subroutine>

  subroutine madapt_getVertexCoords(rmeshAdapt,DvertexCoords)

!<description>
    ! Returns the vertex coordinates
!</description>

!<input>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(in) :: rmeshAdapt
!</input>

!<output>
    ! Array of vertex coordinates
    real(DP), dimension(*), intent(out) :: DvertexCoords
!</output>

!</subroutine>

    ! local variable
    integer :: i,j,n,m
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    
    call storage_getbase_double2d(&
        rmeshAdapt%rtriangulation%h_DvertexCoords, p_DvertexCoords)

    n = size(p_DvertexCoords,1)
    m = size(p_DvertexCoords,2)

    do j=1,m
      do i=1,n
        DvertexCoords(n*(j-1)+i) = p_DvertexCoords(i,j)
      end do
    end do

  end subroutine madapt_getVertexCoords

  ! ***************************************************************************

!<subroutine>

  subroutine madapt_getNeighboursAtElement(rmeshAdapt,IneighboursAtElement)

!<description>
    ! Returns the element-adjacency list
!</description>

!<input>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(in) :: rmeshAdapt
!</input>

!<output>
    ! Element-adjacency array
    integer, dimension(*), intent(out) :: IneighboursAtElement
!</output>

!</subroutine>

!</subroutine>

    ! local variable
    integer :: i,j,n,m
    integer, dimension(:,:), pointer :: p_IneighboursAtElement
    
    call storage_getbase_int2d(&
        rmeshAdapt%rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)

    n = size(p_IneighboursAtElement,1)
    m = size(p_IneighboursAtElement,2)

    do j=1,m
      do i=1,n
        IneighboursAtElement(n*(j-1)+i) = p_IneighboursAtElement(i,j)
      end do
    end do
    
  end subroutine madapt_getNeighboursAtElement

  ! ***************************************************************************

!<subroutine>

  subroutine madapt_getVerticesAtElement(rmeshAdapt,IverticesAtElement)

!<description>
    ! Returns the vertices-at-element list
!</description>

!<input>
    ! Mesh adaptation structure
    type(t_meshAdapt), intent(in) :: rmeshAdapt
!</input>

!<output>
    ! Vertices-at-element array
    integer, dimension(*), intent(out) :: IverticesAtElement
!</output>

!</subroutine>

!</subroutine>

    ! local variable
    integer :: i,j,n,m
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    
    call storage_getbase_int2d(&
        rmeshAdapt%rtriangulation%h_IverticesAtElement, p_IverticesAtElement)

    n = size(p_IverticesAtElement,1)
    m = size(p_IverticesAtElement,2)

    do j=1,m
      do i=1,n
        IverticesAtElement(n*(j-1)+i) = p_IverticesAtElement(i,j)
      end do
    end do
    
  end subroutine madapt_getVerticesAtElement

end module meshadaptbase

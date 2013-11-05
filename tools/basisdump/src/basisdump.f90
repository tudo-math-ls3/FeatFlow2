!##############################################################################
!# ****************************************************************************
!# <name> basisdump </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# </purpose>
!##############################################################################

program basisdump

  use fparser
  use fsystem
  use genoutput
  use storage
  use paramlist
  use basicgeometry
  use triangulation
  use transformation
  use derivatives
  use element
  use elementpreprocessing
  use spdiscprojection
  use ucd
  
  implicit none
  
  type(t_parlist) :: rparam
  type(t_triangulation) :: rtria, rtriaDump
  type(t_evalElementSet) :: reval
  integer :: nref, nbasis, nverts, ibas, ivt, ndone, ntodo, nmaxDer, ider
  integer(I32) :: celement, cshape, cevalTag, ctrafoType
  real(DP), dimension(:,:), pointer :: p_Dvtx
  real(DP), dimension(:,:,:), pointer :: p_Dphi
  integer, dimension(1), parameter :: IelList = (/1/)
  real(DP), dimension(:,:,:,:), pointer :: p_Dbas
  logical, dimension(EL_MAXNDER) :: Bder
  character(len=SYS_STRLEN) :: selement, spredir, sucddir
  type(t_ucdExport) :: rexport

  ! The very first thing in every application:
  ! Initialise system-wide settings:
  call system_init()
  
  ! The very second thing in every program:
  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)
  
  ! Initialise function parser
  call fparser_init()

  if (sys_ncommandLineArgs() .lt. 2) then
    call output_line('USAGE: basisdump <element>')
    call exit
  end if
  
  ! fetch element name
  call getarg(1, selement)

  ! parse the element id
  celement = elem_igetID(selement)
  
  ! fetch the element shape
  cshape = elem_igetShape(celement)

  ! Read in reference mesh
  call output_line("Creating meshes...")
  call createRefMesh(rtria, cshape)
  call tria_initStandardMeshFromRaw(rtria)
  
  ! select number of refinements depending on dimension
  select case(rtria%ndim)
  case (1)
    nref = 9
  case (2)
    nref = 7
  case (3)
    nref = 5
  end select

  ! Refine unto dump level
  call tria_refine2LevelOrdering(rtria, rtriaDump)
  call tria_quickRefine2LevelOrdering (nref-1, rtriaDump)
  call tria_initStandardMeshFromRaw(rtriaDump)
    
  ! Fetch the vertex-coordinates array from the dump triangulation
  call output_line("Preparing for element evaluation...")
  call storage_getbase_double2d(rtriaDump%h_DvertexCoords, p_Dvtx)
  nverts = rtriaDump%NVT
    
  ! fetch the number of basis function
  nbasis = elem_igetNDofLoc(celement)
    
  ! define the evaluation tag
  cevalTag = ior(EL_EVLTAG_REFPOINTS, elem_getEvaluationTag(celement))
    
  ! fetch the trafo type
  ctrafoType = elem_igetTrafoType(celement)
  
  ! fetch the maximum derivative
  nmaxDer = min(elem_getMaxDerivative(celement),1+rtria%ndim)

  ! set up bder
  Bder(:) = .false.
  Bder(1:nmaxDer) = .true.
    
  ! allocate a basis function array
  allocate(p_Dphi(nverts, nmaxDer, nbasis))
    
  ! allocate basis function evaluation array
  allocate(p_Dbas(nbasis,nmaxDer, 100, 1))

    
  call output_line("Evaluating basis functions...")
  ! start evaluating
  ! Notes:
  ! 1. We need to run over the dump-mesh vertices in small blocks, because some old
  !    element implementations cause stack overflows for a sufficiently large number
  !    of evaluation points.
  ! 2. We need to release and reinitialise the element evaluation structure in the
  !    last iteration of the loop, as otherwise the elprep_prepareSetForEvaluation
  !    routine will not recognise that the number of points is less than in the
  !    iterations before.
  ndone = 0
  do while(ndone .lt. nverts)
    
    ntodo = min(100, nverts - ndone)
      
    ! Define an element evaluation set
    call elprep_init(reval)
      
    ! prepare element set
    call elprep_prepareSetForEvaluation (reval, cevalTag, rtria, IelList, ctrafoType, &
        p_Dvtx(:,ndone+1:ndone+ntodo))

    ! evaluate
    call elem_generic_sim2 (celement, reval, Bder, p_Dbas)
      
    ! loop over all basis functions
    do ibas = 1, nbasis
      do ider = 1, nmaxDer
        do ivt = 1, ntodo
          p_Dphi(ndone+ivt, ider, ibas) = p_Dbas(ibas, ider, ivt, 1)
        end do
      end do
    end do
      
    ndone = ndone + ntodo
      
    ! release the element set in the last two iterations
    if(ndone + 100 .gt. nverts) &
      call elprep_releaseElementSet(reval)
      
  end do

  ! Start UCD export to VTK file:
  call output_line("Writing VTK file...")
  if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = "./ucd"
  call ucd_startVTK (rexport, UCD_FLAG_STANDARD, rtriaDump, &
      trim(sucddir)//"/"//trim(selement)  // "_lvl" // trim(sys_sil(nref,4)) //".vtk")
    
  ! loop over all basis functions
  do ider = 1, nmaxDer
    do ibas = 1, nbasis
      ! write UCD variable
      call ucd_addVariableVertexBased (rexport, trim(funcName(ibas,ider)), p_Dphi(:,ider,ibas))
    end do
  end do
    
  ! write and close ucd
  call ucd_write (rexport)
  call ucd_release (rexport)
    
  call output_line("Cleaning up...")
  ! clean up
  deallocate(p_Dphi)
  deallocate(p_Dbas)
  call tria_done(rtriaDump)
  call tria_done(rtria)

  ! Clean up the storage management, finish
  call storage_done()

contains

  character(len=SYS_STRLEN) function funcName(ibas, ider)
  integer, intent(in) :: ibas, ider
  
    select case(ider)
    case (1)
      funcName = "phi_" // trim(sys_si0l(ibas,2))
    case (2)
      funcName = "dx-phi_" // trim(sys_si0l(ibas,2))
    case (3)
      funcName = "dy-phi_" // trim(sys_si0l(ibas,2))
    case (4)
      funcName = "dz-phi_" // trim(sys_si0l(ibas,2))
    end select
  end function

  subroutine createRefMesh(rtria, cshape)
  type(t_triangulation), intent(out) :: rtria
  integer(i32), intent(in) :: cshape

  integer :: i, nvc
  integer, dimension(1) :: Isize1
  integer, dimension(2) :: Isize2
  real(DP), dimension(:,:), pointer :: p_Dvtx, p_Dv2
  integer, dimension(:,:), pointer :: p_Ive
  integer, dimension(:), pointer :: p_Inp
  

    ! decode shape
    rtria%ndim = mod(int(cshape) / 1000000, 100)
    rtria%NVT  = mod(int(cshape)          , 100)
    rtria%NMT  = mod(int(cshape) /     100, 100)
    rtria%NAT  = mod(int(cshape) /   10000, 100)
    rtria%NEL  = 1
    rtria%NNVE = rtria%nvt
    rtria%NNEE = rtria%nmt
    rtria%NNAE = rtria%nat
    rtria%NBCT = 1

    rtria%InelOfType(rtria%NNVE) = 1

    select case(cshape)
    case (BGEOM_SHAPE_TRIA, BGEOM_SHAPE_TETRA)
      nvc = rtria%ndim + 1
    case default
      nvc = rtria%ndim
    end select
    
    Isize2 = (/rtria%ndim,rtria%NNVE/)
    call storage_new('createRefMesh', 'DvertexCoords', Isize2, &
        ST_DOUBLE, rtria%h_DvertexCoords, ST_NEWBLOCK_NOINIT)
    call storage_getbase_double2d(rtria%h_DvertexCoords, p_Dvtx)
    
    Isize2 = (/rtria%NNVE,1/)
    call storage_new('createRefMesh', 'IverticesAtElement', Isize2, &
        ST_INT, rtria%h_IverticesAtElement, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int2d(rtria%h_IverticesAtElement, p_Ive)

    Isize1 = (/rtria%NNVE/)
    call storage_new('createRefMesh', 'InodalProperty', Isize1, &
        ST_INT, rtria%h_InodalProperty, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rtria%h_InodalProperty, p_Inp)
    do i=1, rtria%NNVE
      p_Inp(i) = 1
    end do
    
    if(nvc .eq. rtria%ndim) then
      ! create corner vertices
      call spdp_aux_getCornerRefCoords(p_Dvtx, rtria%ndim, rtria%NNVE)
    else
      allocate(p_Dv2(nvc, rtria%NNVE))
      ! create corner vertices
      call spdp_aux_getCornerRefCoords(p_Dv2, rtria%ndim, rtria%NNVE)
      do i = 1, rtria%NNVE
        p_Dvtx(1:rtria%ndim, i) = p_Dv2(2:rtria%ndim+1, i)
      end do
      deallocate(p_Dv2)
    end if

    do i = 1, rtria%NNVE
      p_Ive(i,1) = i
    end do
    
    call tria_genRawBoundary2D(rtria)
    call tria_initExtendedRawMesh(rtria)
    
  end subroutine

end program

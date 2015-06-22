module stokesdbg_core

use fsystem
use storage
use genoutput
use linearsolver
use boundary
use cubature
use derivatives
use matrixfilters
use vectorfilters
use discretebc
use bcassembly
use triangulation
use meshmodification
use domainintegration
use element
use spatialdiscretisation
use linearsystemscalar
use linearsystemblock
use coarsegridcorrection
use multileveloperators
use multilevelprojection
use spdiscprojection
use filtersupport
use scalarpde
use linearformevaluation
use linearalgebra
use discretebc
use ucd
use collection, only: t_collection
use pprocerror
use bilinearformevaluation
use stdoperators
use paramlist
use quicksolver
use iterationcontrol
use jumpstabilisation

use stokesdbg_aux
use stokesdbg_div_eoj
  
implicit none

!<types>

!<typeblock description="Type block defining all information about one level">

  type t_level
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtria

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscr
    
    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubInfo
    
    ! A system matrix for that specific level.
    type(t_matrixBlock) :: rmatSys

    ! Prolongation matrix for velocity
    type(t_matrixScalar) :: rmatProlVelo

    ! Prolongation matrix for pressure
    type(t_matrixScalar) :: rmatProlPres

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC

    ! Multilevel projection structure
    type(t_interlevelProjectionBlock) :: rprojection
  
  end type
  
!</typeblock>

!<typeblock>

  type t_problem
  
    ! A pointer to the parameter structure
    type(t_parlist), pointer :: p_rparam => null()
  
    ! the problem dirver
    integer :: idriver = 0
    
    ! minimum test level
    integer :: ilevelMin = 0
    
    ! maximum test level
    integer :: ilevelMax = 0
    
    ! coarse mesh level for multigrid
    integer :: ilevelCoarse = 0
    
    ! level array; dimension(ilevelCoarse:ilevelMax)
    type(t_level), dimension(:), pointer :: Rlevels => null()
    
    ! the boundary structure
    type(t_boundary) :: rbnd
    
    ! the dimension of the space
    integer :: ndim = 0
    
    ! number of velocity components
    integer :: ncompVelo = 0
    
    ! stabilisation
    integer :: stabilise = 0
    
    ! eoj parameter
    real(DP) :: deojGamma = 0.0_DP
    
    ! eoj cubature id
    integer(I32) :: ceojCubature = CUB_G3_1D
    
    ! problem parameters
    real(DP) :: dnu = 1.0_DP
    real(DP) :: dalpha = 1.0_DP
    real(DP) :: dbeta = 1.0_DP
    real(DP) :: dgamma = 1.0_DP
    real(DP) :: dsigma = 1.0_DP
    
    ! statistics array; dimension(DSTAT_LENGTH,ilevelMin:ilevelMax)
    real(DP), dimension(:,:), pointer :: p_Dstat
    
    ! statistics array, dimension(ISTAT_LENGTH,ilevelMin:ilevelMax)
    integer, dimension(:,:), pointer :: p_Istat
  
  end type
  
!</typeblock>

!<typeblock>
  
  type t_system
  
    ! A pointer to the underlying problem structure
    type(t_problem), pointer :: p_rproblem => null()
  
    ! the level this system structure is assigned to
    integer :: ilevel = 0
    
    ! the corresponding vectors
    type(t_vectorBlock) :: rvecSol
    type(t_vectorBlock) :: rvecRhs
    type(t_vectorBlock) :: rvecDef
    
    ! intermediate vectors
    type(t_vectorBlock) :: rvecRhsVelo, rvecRhsPres
    type(t_vectorBlock) :: rvecDefVelo, rvecDefPres
    
    ! the filter chain to be used
    type(t_filterChain), dimension(:), pointer :: p_rfilterChain => null()
    
    ! the linear solver to be used
    type(t_linsolNode), pointer :: p_rsolver => null()
    
    ! the iteration control structrue
    type(t_iterationControl) :: riter
  
  end type

!</typeblock>
  
!</types>

!<constants>
  ! index of velocity L2-error
  integer, parameter :: DSTAT_U_L2   = 1
  ! index of velocity H1-error
  integer, parameter :: DSTAT_U_H1   = 2
  ! index of velocity divergence
  integer, parameter :: DSTAT_U_DIV  = 3
  ! index of pressure L2-error
  integer, parameter :: DSTAT_P_L2   = 4
  ! index of initial defect
  integer, parameter :: DSTAT_DEF_I  = 5
  ! index of final defect
  integer, parameter :: DSTAT_DEF_F  = 6
  ! index of convergence rate
  integer, parameter :: DSTAT_CRATE  = 7
  ! length of the p_Dstat array 
  integer, parameter :: DSTAT_LENGTH = 7
  
  ! index of the solver iteration count
  integer, parameter :: ISTAT_NITER  = 1
  ! index of the solver status
  integer, parameter :: ISTAT_STATUS = 2
  ! length of the p_Istat array
  integer, parameter :: ISTAT_LENGTH = 2
!</constants>
contains
  
  ! ***********************************************************************************************
  
  subroutine stdbg_doneProblem(rproblem)
  type(t_problem), intent(inout) :: rproblem
  
  integer :: i
  
    if(associated(rproblem%p_Istat)) &
      deallocate(rproblem%p_Istat)
    if(associated(rproblem%p_Dstat)) &
      deallocate(rproblem%p_Dstat)
  
    if(associated(rproblem%Rlevels)) then
    
      ! loop over all levels
      do i = ubound(rproblem%Rlevels,1), lbound(rproblem%Rlevels, 1), -1
      
        call mlprj_doneProjection(rproblem%Rlevels(i)%rprojection)

        call lsysbl_releaseMatrix(rproblem%Rlevels(i)%rmatSys)
        
        call lsyssc_releaseMatrix(rproblem%Rlevels(i)%rmatProlPres)
        call lsyssc_releaseMatrix(rproblem%Rlevels(i)%rmatProlVelo)
        
        call bcasm_releaseDiscreteBC (rproblem%Rlevels(i)%rdiscreteBC)
        call spdiscr_releaseCubStructure(rproblem%Rlevels(i)%rcubInfo)
        call spdiscr_releaseBlockDiscr(rproblem%Rlevels(i)%rdiscr)
        call tria_done (rproblem%Rlevels(i)%rtria)
        
      end do
      
      deallocate(rproblem%Rlevels)
      
    end if
    
    ! release boundary
    call boundary_release(rproblem%rbnd)

  end subroutine
  
  ! ***********************************************************************************************

  subroutine stdbg_initProblem(rproblem, rparam)
  type(t_problem), intent(out) :: rproblem
  type(t_parlist), target, intent(in) :: rparam

  integer :: ilmin, ilmax, ilcrs
  character(len=64) :: stxt
  
    ! store parameter pointer
    rproblem%p_rparam => rparam
  
    ! fetch parameters
    call parlst_getvalue_int(rparam, '', 'LEVEL_MIN', ilmin, 1)
    call parlst_getvalue_int(rparam, '', 'LEVEL_MAX', ilmax, -1)
    call parlst_getvalue_int(rparam, '', 'LEVEL_COARSE', ilcrs, 0)
    
    ! adjust parameters
    ilmax = max(1, ilmax)
    if(ilmin .le. 0) then
      ilmin = max(1, ilmax+ilmin)
    else
      ilmin = min(ilmin, ilmax)
    end if
  
    if(ilcrs .le. 0) then
      ilcrs = max(1, ilmin+ilcrs)
    else
      ilcrs = min(ilmin, ilcrs)
    end if
    
    ! store parameters
    rproblem%ilevelMin = ilmin
    rproblem%ilevelMax = ilmax
    rproblem%ilevelCoarse = ilcrs
    
    ! fetch stabilisation parameters
    call parlst_getvalue_int(rparam, '', 'STABILISE', rproblem%stabilise, 0)
    call parlst_getvalue_double(rparam, '', 'EOJ_GAMMA', rproblem%deojGamma, 0.01_DP)
    call parlst_getvalue_string(rparam, '', 'EOJ_CUBATURE', stxt, 'G3_1D')
    rproblem%ceojCubature = cub_igetID(stxt)
    
    ! allocate levels
    allocate(rproblem%Rlevels(ilcrs:ilmax))
    
    ! fetch other parameters
    call parlst_getvalue_double(rparam, '', 'DNU', rproblem%dnu, 1.0_DP)
    call parlst_getvalue_double(rparam, '', 'DALPHA', rproblem%dalpha, 1.0_DP)
    call parlst_getvalue_double(rparam, '', 'DBETA', rproblem%dbeta, 1.0_DP)
    call parlst_getvalue_double(rparam, '', 'DGAMMA', rproblem%dgamma, 1.0_DP)
    call parlst_getvalue_double(rparam, '', 'DSIGMA', rproblem%dsigma, 1.0_DP)
    
    ! allocate statistics arrays
    allocate(rproblem%p_Dstat(DSTAT_LENGTH,ilmin:ilmax))
    allocate(rproblem%p_Istat(ISTAT_LENGTH,ilmin:ilmax))
    
  end subroutine
  
  ! ***********************************************************************************************
  
  subroutine stdbg_initTriangulation2D(rproblem, sfile)
  type(t_problem), target, intent(inout) :: rproblem
  character(len=*), intent(in) :: sfile
  
  integer :: i, ilcrs, ilmin, ilmax, idistType
  real(DP) :: ddistFactor
  character(len=SYS_STRLEN) :: spredir
  
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './pre'

    ilcrs = rproblem%ilevelCoarse
    ilmin = rproblem%ilevelMin
    ilmax = rproblem%ilevelMax
    
    ! fetch distortion type
    call parlst_getvalue_int(rproblem%p_rparam, '', 'DIST_TYPE', idistType, 0)
    call parlst_getvalue_double(rproblem%p_rparam, '', 'DIST_FACTOR', ddistFactor, 0.1_DP)

    call output_line("Generating mesh hierarchy...")
    
    ! read boundary
    call boundary_read_prm(rproblem%rbnd, trim(spredir) // '/' // trim(sfile) // '.prm')
    
    ! read coarse mesh
    call tria_readTriFile2D (rproblem%Rlevels(ilcrs)%rtria, &
        trim(spredir) // '/' // trim(sfile) // '.tri', rproblem%rbnd)
    
    ! pre-refine
    call tria_quickRefine2LevelOrdering (ilcrs-1, rproblem%Rlevels(ilcrs)%rtria, rproblem%rbnd)
    
    ! initialise a standard mesh
    call tria_initStandardMeshFromRaw (rproblem%Rlevels(ilcrs)%rtria, rproblem%rbnd)
    
    ! distort coarse mesh if desired
    if(idistType .eq. 1) &
      call meshmod_disturbMesh (rproblem%Rlevels(ilcrs)%rtria, ddistFactor)
    
    ! Now refine the grid for the fine levels.
    do i = ilcrs+1, ilmax
      call tria_refine2LevelOrdering(rproblem%Rlevels(i-1)%rtria, rproblem%Rlevels(i)%rtria, rproblem%rbnd)
      call tria_initStandardMeshFromRaw(rproblem%Rlevels(i)%rtria, rproblem%rbnd)
    end do
    
    ! distort all meshes if desired
    if(idistType .eq. 2) then
      do i = ilmin, ilmax
        call meshmod_disturbMesh (rproblem%Rlevels(i)%rtria, ddistFactor)
      end do
    end if
    
  end subroutine
  
  ! ***********************************************************************************************
  
  subroutine stdbg_initDiscretisation(rproblem)
  type(t_problem), target, intent(inout) :: rproblem
  
  integer(I32):: celemVelo, celemPres, ccubature
  integer :: i, j, ilcrs, ilmax, ndim, ncomp
  type(t_level), pointer :: p_rlvl
  character(len=SYS_STRLEN) :: sname
  
    ilcrs = rproblem%ilevelCoarse
    ilmax = rproblem%ilevelMax
    ndim = rproblem%Rlevels(ilcrs)%rtria%ndim
    rproblem%ndim = ndim
    
    ! read parameters
    call parlst_getvalue_string(rproblem%p_rparam, '', 'ELEMENT_VELOCITY', sname, '')
    celemVelo = elem_igetID(sname)
    call parlst_getvalue_string(rproblem%p_rparam, '', 'ELEMENT_PRESSURE', sname, '')
    celemPres = elem_igetID(sname)
    call parlst_getvalue_string(rproblem%p_rparam, '', 'CUBATURE', sname, '')
    ccubature = cub_igetID(sname)
    
    ! fetch number of velocity FE components
    rproblem%ncompVelo = elem_igetFeDimension(celemVelo)
    
    ! compute total number of components
    if(rproblem%ncompVelo .eq. 1) then
      ! ndim*Velocity + 1*Pressure
      ncomp = ndim + 1
    else
      ! 1*Velocity + 1*Pressure
      ncomp = 2
    end if
    
    ! loop over all problem levels
    do i = ilcrs, ilmax
    
      p_rlvl => rproblem%Rlevels(i)

      ! Create the block discretisation
      if(ndim .eq. 2) then
        call spdiscr_initBlockDiscr (p_rlvl%rdiscr, ncomp, p_rlvl%rtria, rproblem%rbnd)
      else
        call spdiscr_initBlockDiscr (p_rlvl%rdiscr, ncomp, p_rlvl%rtria)
      end if
      
      ! Set up velocity spaces
      if(ndim .eq. 2) then
        call spdiscr_initDiscr_simple (p_rlvl%rdiscr%RspatialDiscr(1), &
            celemVelo, ccubature, p_rlvl%rtria, rproblem%rbnd)
      else
        call spdiscr_initDiscr_simple (p_rlvl%rdiscr%RspatialDiscr(1), &
            celemVelo, ccubature, p_rlvl%rtria)
      end if

      if(rproblem%ncompVelo .eq. 1) then
        ! scalar velocity space: duplicate for all components
        do j = 2, ndim
          call spdiscr_duplicateDiscrSc(p_rlvl%rdiscr%RspatialDiscr(1), p_rlvl%rdiscr%RspatialDiscr(j))
        end do
      end if

      ! Set up pressure space
      call spdiscr_deriveSimpleDiscrSc (p_rlvl%rdiscr%RspatialDiscr(1), &
          celemPres, ccubature, p_rlvl%rdiscr%RspatialDiscr(ncomp))
    
      ! Set up cubature info structure
      call spdiscr_createDefCubStructure(p_rlvl%rdiscr%RspatialDiscr(1), &
          p_rlvl%rcubInfo, ccubature)
      
    end do
  
  end subroutine
  
  ! ***********************************************************************************************

  subroutine stdbg_initProjections(rproblem)
  type(t_problem), target, intent(inout) :: rproblem
  
  integer :: i, j, ilcrs, ilmax, ndim, ncomp
  type(t_level), pointer :: p_rlvlf, p_rlvlc
  
    ilcrs = rproblem%ilevelCoarse
    ilmax = rproblem%ilevelMax
    ndim = rproblem%Rlevels(ilcrs)%rdiscr%ndimension
    ncomp = rproblem%Rlevels(ilcrs)%rdiscr%ncomponents

    ! Initialise multi-level projection for coarse level
    call output_line('Assembling multilevel projections...')
    call mlprj_initProjectionDiscr (rproblem%Rlevels(ilcrs)%rprojection, rproblem%Rlevels(ilcrs)%rdiscr)
    
    do i = ilcrs+1, ilmax
    
      p_rlvlf => rproblem%Rlevels(i)
      p_rlvlc => rproblem%Rlevels(i-1)

      ! Create prolongation matrix structures
      call mlop_create2LvlMatrixStruct(p_rlvlc%rdiscr%RspatialDiscr(1), &
          p_rlvlf%rdiscr%RspatialDiscr(1), LSYSSC_MATRIX9, p_rlvlf%rmatProlVelo)
      call mlop_create2LvlMatrixStruct(p_rlvlc%rdiscr%RspatialDiscr(ncomp), &
          p_rlvlf%rdiscr%RspatialDiscr(ncomp), LSYSSC_MATRIX9, p_rlvlf%rmatProlPres)

      ! Assemble prolongation matrices
      call mlop_build2LvlProlMatrix (p_rlvlc%rdiscr%RspatialDiscr(1),&
          p_rlvlf%rdiscr%RspatialDiscr(1), .true., p_rlvlf%rmatProlVelo)
      call mlop_build2LvlProlMatrix (p_rlvlc%rdiscr%RspatialDiscr(ncomp),&
          p_rlvlf%rdiscr%RspatialDiscr(ncomp), .true., p_rlvlf%rmatProlPres)

      ! Initialise multi-level projection
      call mlprj_initProjectionDiscr (p_rlvlf%rprojection, p_rlvlf%rdiscr)
      do j = 1, ncomp-1
        call mlprj_initMatrixProjection(&
            p_rlvlf%rprojection%RscalarProjection(1,j), p_rlvlf%rmatProlVelo)
      end do
      call mlprj_initMatrixProjection(&
          p_rlvlf%rprojection%RscalarProjection(1,ncomp), p_rlvlf%rmatProlPres)

    end do
    
  end subroutine
  
  ! ***********************************************************************************************

  subroutine stdbg_initMatrices(rproblem, cconstrType)
  type(t_problem), target, intent(inout) :: rproblem
  integer(I32), optional, intent(in) :: cconstrType
  
  integer :: i, j, k, ilcrs, ilmax, ndim, ctype
  type(t_level), pointer :: p_rlvl
  integer, dimension(3) :: Ideriv
  
    ilcrs = rproblem%ilevelCoarse
    ilmax = rproblem%ilevelMax
    ndim = rproblem%Rlevels(ilcrs)%rdiscr%ndimension
    
    if(ndim .eq. 2) then
      Ideriv = (/ DER_DERIV2D_X, DER_DERIV2D_Y, DER_FUNC /)
    else if(ndim .eq. 3) then
      Ideriv = (/ DER_DERIV3D_X, DER_DERIV3D_Y, DER_DERIV3D_Z /)
    end if
    
    ! choose construction type
    select case(rproblem%stabilise)
    case (1,2,3,4)
      ctype = BILF_MATC_EDGEBASED
    case default
      ctype = BILF_MATC_ELEMENTBASED
    end select
    
    do i = ilcrs, ilmax

      call output_line('Assembling Level ' // trim(sys_sil(i,4)) // '...')

      p_rlvl => rproblem%Rlevels(i)

      ! Create block matrix
      call lsysbl_createMatBlockByDiscr (p_rlvl%rdiscr, p_rlvl%rmatSys)

      ! Houston, we have a saddle point problem...
      p_rlvl%rmatSys%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
      
      ! Assemble A-matrix structure
      call bilf_createMatrixStructure(p_rlvl%rdiscr%RspatialDiscr(1),&
          LSYSSC_MATRIX9, p_rlvl%rmatSys%RmatrixBlock(1,1), cconstrType=ctype)

      ! Assemble B-matrix structure
      call bilf_createMatrixStructure (p_rlvl%rdiscr%RspatialDiscr(ndim+1),&
          LSYSSC_MATRIX9, p_rlvl%rmatSys%RmatrixBlock(1,ndim+1),&
          p_rlvl%rdiscr%RspatialDiscr(1))

      ! Transpose B-matrix to obtain D-matrix structure
      call lsyssc_transposeMatrix (p_rlvl%rmatSys%RmatrixBlock(1,ndim+1), &
          p_rlvl%rmatSys%RmatrixBlock(ndim+1,1), LSYSSC_TR_STRUCTURE)

      ! Duplicate matrix structures
      do j = 2, ndim
        call lsyssc_duplicateMatrix (p_rlvl%rmatSys%RmatrixBlock(1,1),&
            p_rlvl%rmatSys%RmatrixBlock(j,j),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

        call lsyssc_assignDiscrDirectMat (p_rlvl%rmatSys%RmatrixBlock(j,j),&
            p_rlvl%rdiscr%RspatialDiscr(j))
        
        call lsyssc_duplicateMatrix (p_rlvl%rmatSys%RmatrixBlock(1,ndim+1),&
            p_rlvl%rmatSys%RmatrixBlock(j,ndim+1),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

        call lsyssc_duplicateMatrix (p_rlvl%rmatSys%RmatrixBlock(ndim+1,1),&
            p_rlvl%rmatSys%RmatrixBlock(ndim+1,j),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      end do
      
      ! duplicate matrix structures for div-div eoj
      select case(rproblem%stabilise)
      case (3,4)
        do j = 1, ndim
          do k = 1, ndim
            if (k .ne. j) then
              call lsyssc_duplicateMatrix (p_rlvl%rmatSys%RmatrixBlock(j,j),&
                p_rlvl%rmatSys%RmatrixBlock(j,k),LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
              call lsyssc_assignDiscrDirectMat (p_rlvl%rmatSys%RmatrixBlock(j,k),&
                p_rlvl%rdiscr%RspatialDiscr(k))
              call lsyssc_allocEmptyMatrix(p_rlvl%rmatSys%RmatrixBlock(j,k), LSYSSC_SETM_ZERO)
            end if
          end do
        end do
      end select

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Matrix Content Assembly
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      ! Assemble laplace matrix for X-velocity
      call stdop_assembleLaplaceMatrix (p_rlvl%rmatSys%RmatrixBlock(1,1),&
          .true., rproblem%dnu, p_rlvl%rcubInfo)

      ! Copy the matrix into the Y-velocity block
      do j = 2, ndim
        call lsyssc_duplicateMatrix (p_rlvl%rmatSys%RmatrixBlock(1,1),&
            p_rlvl%rmatSys%RmatrixBlock(j,j),LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
      end do

      ! Assemble the B-matrices
      do j = 1, ndim
        call stdop_assembleSimpleMatrix (p_rlvl%rmatSys%RmatrixBlock(j,ndim+1),&
          DER_FUNC, Ideriv(j), -1.0_DP, .true., p_rlvl%rcubInfo)
      end do

      ! And transpose the B-matrices to obtain the D-matrices
      do j = 1, ndim
        call lsyssc_transposeMatrix (p_rlvl%rmatSys%RmatrixBlock(j,ndim+1), &
            p_rlvl%rmatSys%RmatrixBlock(ndim+1,j), LSYSSC_TR_CONTENT)
      end do

      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      ! Stabilisation Assembly
      ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      select case(rproblem%stabilise)
      case (1)
        ! reactive jump-stabilisation
        do j = 1, ndim
          call jstab_calcReacJumpStabilisation(p_rlvl%rmatSys%RmatrixBlock(j,j), &
            rproblem%deojGamma, 1.0_DP, rproblem%ceojCubature, rproblem%dnu)
        end do

      case (2)
        ! gradient jump-stabilisation
        do j = 1, ndim
          call jstab_calcUEOJumpStabilisation(p_rlvl%rmatSys%RmatrixBlock(j,j), &
            rproblem%deojGamma, 0.0_DP, 2.0_DP, 1.0_DP, rproblem%ceojCubature, rproblem%dnu)
        end do

      case (3)
        ! divergence jump-stabilisation
        do j = 1, ndim
          call stokesdbg_divEoj(p_rlvl%rmatSys%RmatrixBlock, &
            rproblem%deojGamma, 0.0_DP, 2.0_DP, 1.0_DP, rproblem%ceojCubature, rproblem%dnu)
        end do

      case (4)
        ! normal-flow jump-stabilisation
        do j = 1, ndim
          call stokesdbg_flowEoj(p_rlvl%rmatSys%RmatrixBlock, rproblem%deojGamma, rproblem%ceojCubature)
        end do
      end select
      
    end do
    
  end subroutine
  
  ! ***********************************************************************************************

  subroutine stdbg_initNoSlipBCs(rproblem)
  type(t_problem), target, intent(inout) :: rproblem
  
  type(t_boundaryRegion) :: rrgn
  integer :: i, j, ncomp
  
    if(rproblem%ncompVelo .gt. 1) then
      ncomp = 1
    else
      ncomp = rproblem%ndim
    end if
  
    do i = rproblem%ilevelCoarse, rproblem%ilevelMax

      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(rproblem%Rlevels(i)%rdiscreteBC)

      call boundary_createRegion(rproblem%rbnd,1,0,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      do j = 1, ncomp
        call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,j,rrgn,&
            rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      end do
      
      ! Assign BCs to system matrix
      rproblem%Rlevels(i)%rmatSys%p_rdiscreteBC => rproblem%Rlevels(i)%rdiscreteBC

      ! Filter system matrix
      call matfil_discreteBC(rproblem%Rlevels(i)%rmatSys)

    end do
    
  end subroutine
  
  ! ***********************************************************************************************

  subroutine stdbg_initSlipBCs(rproblem)
  type(t_problem), target, intent(inout) :: rproblem
  
  type(t_boundaryRegion) :: rrgn
  integer :: i
  
    ! ensure that we have a scalar velocity space
    if(rproblem%ncompVelo .ne. 1) then
      call output_line("ERROR: Slip boundary conditions not available "//&
        "for vector-valued Velocity spaces!", OU_CLASS_ERROR, OU_MODE_STD, 'stdbg_initSlipBCs')
      call sys_halt()
    end if
  
    ! loop over all levels
    do i = rproblem%ilevelCoarse, rproblem%ilevelMax

      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(rproblem%Rlevels(i)%rdiscreteBC)

      ! Bottom edge: u2 = 0
      call boundary_createRegion(rproblem%rbnd,1,1,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,2,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)

      ! Right edge: u1 = 0
      call boundary_createRegion(rproblem%rbnd,1,2,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,1,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)

      ! Top edge: u2 = 0
      call boundary_createRegion(rproblem%rbnd,1,3,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,2,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)

      ! Left edge: u1 = 0
      call boundary_createRegion(rproblem%rbnd,1,4,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,1,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      
      ! Assign BCs to system matrix
      rproblem%Rlevels(i)%rmatSys%p_rdiscreteBC => rproblem%Rlevels(i)%rdiscreteBC

      ! Filter system matrix
      call matfil_discreteBC(rproblem%Rlevels(i)%rmatSys)
      
    end do
    
  end subroutine
  
  ! ***********************************************************************************************

  subroutine stdbg_initQuadPoiseulleBCs(rproblem)
  type(t_problem), target, intent(inout) :: rproblem
  
  type(t_boundaryRegion) :: rrgn
  integer :: i, j, ncomp
  
    if(rproblem%ncompVelo .gt. 1) then
      ncomp = 1
    else
      ncomp = rproblem%ndim
    end if
  
    do i = rproblem%ilevelCoarse, rproblem%ilevelMax

      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(rproblem%Rlevels(i)%rdiscreteBC)

      ! Bottom edge: u = 0
      call boundary_createRegion(rproblem%rbnd,1,1,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      do j = 1, ncomp
        call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,j,rrgn,&
            rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      end do

      ! Top edge: u = 0
      call boundary_createRegion(rproblem%rbnd,1,3,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      do j = 1, ncomp
        call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,j,rrgn,&
            rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      end do

      ! Left edge: u1 = par-profile, u2 = 0
      call boundary_createRegion(rproblem%rbnd,1,4,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,1,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcParProfileBC2D)
      if(ncomp .eq. 2) then
        call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,2,rrgn,&
            rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      end if
      
      ! Assign BCs to system matrix
      rproblem%Rlevels(i)%rmatSys%p_rdiscreteBC => rproblem%Rlevels(i)%rdiscreteBC

      ! Filter system matrix
      call matfil_discreteBC(rproblem%Rlevels(i)%rmatSys)
      
    end do
    
  end subroutine

  ! ***********************************************************************************************

  subroutine stdbg_initSingDrivenCavityBCs(rproblem)
  type(t_problem), target, intent(inout) :: rproblem
  
  type(t_boundaryRegion) :: rrgn
  integer :: i, j, ncomp
  
    if(rproblem%ncompVelo .gt. 1) then
      ncomp = 1
    else
      ncomp = rproblem%ndim
    end if
  
    do i = rproblem%ilevelCoarse, rproblem%ilevelMax

      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(rproblem%Rlevels(i)%rdiscreteBC)

      ! Top edge: u1 = 1, u2=0
      call boundary_createRegion(rproblem%rbnd,1,3,rrgn)
      rrgn%iproperties = 0
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,1,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcOneBC2D)
      do j = 2, ncomp
        call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,j,rrgn,&
            rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      end do

      ! Bottom edge: u = 0
      call boundary_createRegion(rproblem%rbnd,1,1,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      do j = 1, ncomp
        call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,j,rrgn,&
            rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      end do

      ! Right edge: u = 0
      call boundary_createRegion(rproblem%rbnd,1,2,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      do j = 1, ncomp
        call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,j,rrgn,&
            rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      end do

      ! Left edge: u = 0
      call boundary_createRegion(rproblem%rbnd,1,4,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      do j = 1, ncomp
        call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,j,rrgn,&
            rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      end do
      
      ! Assign BCs to system matrix
      rproblem%Rlevels(i)%rmatSys%p_rdiscreteBC => rproblem%Rlevels(i)%rdiscreteBC

      ! Filter system matrix
      call matfil_discreteBC(rproblem%Rlevels(i)%rmatSys)
      
    end do
    
  end subroutine

  ! ***********************************************************************************************

  subroutine stdbg_initRegDrivenCavityBCs(rproblem)
  type(t_problem), target, intent(inout) :: rproblem
  
  type(t_boundaryRegion) :: rrgn
  type(t_collection) :: rcollect
  integer :: i, j, ncomp
  
    if(rproblem%ncompVelo .gt. 1) then
      ncomp = 1
    else
      ncomp = rproblem%ndim
    end if
    
    rcollect%DquickAccess(1) = rproblem%dsigma
  
    do i = rproblem%ilevelCoarse, rproblem%ilevelMax

      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(rproblem%Rlevels(i)%rdiscreteBC)

      ! Top edge: u1 = 1, u2=0
      call boundary_createRegion(rproblem%rbnd,1,3,rrgn)
      rrgn%iproperties = 0
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,1,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcRegDrivenCavityBC2D, rcollect)
      do j = 2, ncomp
        call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,j,rrgn,&
            rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      end do

      ! Bottom edge: u = 0
      call boundary_createRegion(rproblem%rbnd,1,1,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      do j = 1, ncomp
        call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,j,rrgn,&
            rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      end do

      ! Right edge: u = 0
      call boundary_createRegion(rproblem%rbnd,1,2,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      do j = 1, ncomp
        call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,j,rrgn,&
            rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      end do

      ! Left edge: u = 0
      call boundary_createRegion(rproblem%rbnd,1,4,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      do j = 1, ncomp
        call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,j,rrgn,&
            rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      end do
      
      ! Assign BCs to system matrix
      rproblem%Rlevels(i)%rmatSys%p_rdiscreteBC => rproblem%Rlevels(i)%rdiscreteBC

      ! Filter system matrix
      call matfil_discreteBC(rproblem%Rlevels(i)%rmatSys)
      
    end do
    
  end subroutine
  ! ***********************************************************************************************
  
  subroutine stdbg_doneSystem(rsystem)
  type(t_system), intent(inout) :: rsystem
  
    if(associated(rsystem%p_rsolver)) then
      call linsol_releaseSolver (rsystem%p_rsolver)
    end if
    
    if(associated(rsystem%p_RfilterChain)) &
      deallocate(rsystem%p_RfilterChain)
    
    if(rsystem%rvecDefPres%NEQ .gt. 0) &
      call lsysbl_releaseVector(rsystem%rvecDefPres)
    if(rsystem%rvecDefVelo%NEQ .gt. 0) &
      call lsysbl_releaseVector(rsystem%rvecDefVelo)
    if(rsystem%rvecRhsPres%NEQ .gt. 0) &
      call lsysbl_releaseVector(rsystem%rvecRhsPres)
    if(rsystem%rvecRhsVelo%NEQ .gt. 0) &
      call lsysbl_releaseVector(rsystem%rvecRhsVelo)
    if(rsystem%rvecDef%NEQ .gt. 0) &
      call lsysbl_releaseVector(rsystem%rvecDef)
    if(rsystem%rvecRhs%NEQ .gt. 0) &
      call lsysbl_releaseVector(rsystem%rvecRhs)
    if(rsystem%rvecSol%NEQ .gt. 0) &
      call lsysbl_releaseVector(rsystem%rvecSol)
  
  end subroutine
  
  ! ***********************************************************************************************

  subroutine stdbg_initSystem(rsystem, rproblem, ilevel)
  type(t_system), intent(out) :: rsystem
  type(t_problem), target, intent(in) :: rproblem
  integer, intent(in) :: ilevel
  
  type(t_level), pointer :: p_rlvl
  
    ! store problem pointer
    rsystem%p_rproblem => rproblem

    ! store system level
    rsystem%ilevel = ilevel
    
    ! fetch problem level
    p_rlvl => rproblem%Rlevels(ilevel)

    ! Create three vectors
    call lsysbl_createVecBlockIndMat (p_rlvl%rmatSys, rsystem%rvecSol, .true.)
    call lsysbl_createVecBlockIndMat (p_rlvl%rmatSys, rsystem%rvecRhs, .true.)
    call lsysbl_createVecBlockIndMat (p_rlvl%rmatSys, rsystem%rvecDef, .true.)
    call lsysbl_createVecBlockIndMat (p_rlvl%rmatSys, rsystem%rvecRhsVelo, .true.)
    call lsysbl_createVecBlockIndMat (p_rlvl%rmatSys, rsystem%rvecRhsPres, .true.)
    call lsysbl_createVecBlockIndMat (p_rlvl%rmatSys, rsystem%rvecDefVelo, .true.)
    call lsysbl_createVecBlockIndMat (p_rlvl%rmatSys, rsystem%rvecDefPres, .true.)

    ! Assemble the RHS vector
    !rform%itermCount = 1
    !rform%Idescriptors(1) = DER_FUNC
    !rcollection%DquickAccess(1) = dnu
    !rcollection%DquickAccess(2) = dc
    !call linf_buildVectorScalar (rform,.true., rsystem%rvecRhs%RvectorBlock(1), &
    !    p_rlvl%rcubInfo, funcRhsX2D, rcollection)
    !call linf_buildVectorScalar (rform,.true., rsystem%rvecRhs%RvectorBlock(2), &
    !    p_rlvl%rcubInfo, funcRhsY2D, rcollection)

    ! Assign BCs to vectors
    !call lsysbl_assignDiscreteBC(rsystem%rvecSol, p_rlvl%rdiscreteBC)
    !call lsysbl_assignDiscreteBC(rsystem%rvecRhs, p_rlvl%rdiscreteBC)
    !call lsysbl_assignDiscreteBC(rsystem%rvecDef, p_rlvl%rdiscreteBC)
    !call lsysbl_assignDiscreteBC(rsystem%rvecRhsVelo, p_rlvl%rdiscreteBC)
    !call lsysbl_assignDiscreteBC(rsystem%rvecRhsPres, p_rlvl%rdiscreteBC)
    !call lsysbl_assignDiscreteBC(rsystem%rvecDefVelo, p_rlvl%rdiscreteBC)
    !call lsysbl_assignDiscreteBC(rsystem%rvecDefPres, p_rlvl%rdiscreteBC)

    ! Filter solution and rhs vectors
    !call vecfil_discreteBCsol(rvecSol)
    !call vecfil_discreteBCrhs(rvecRhs)
    
  end subroutine

  ! ***********************************************************************************************
  
  subroutine stdbg_initMultigrid(rsystem)
  type(t_system), target, intent(inout) :: rsystem

  type(t_linsolNode), pointer :: p_rsmoother, p_rsolver
  type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
  type(t_matrixBlock), dimension(:), pointer :: Rmatrices
  integer :: i, ilvl, ilcrs
  integer :: nsmoothSteps
  real(DP) :: ddamping
  
    ilvl = rsystem%ilevel
    ilcrs = rsystem%p_rproblem%ilevelCoarse
    
    ! fetch MG parameters
    call parlst_getvalue_int(rsystem%p_rproblem%p_rparam, '', 'SMOOTH_STEPS', nsmoothSteps, 4)
    call parlst_getvalue_double(rsystem%p_rproblem%p_rparam, '', 'SMOOTH_DAMP', ddamping, 1.0_DP)
    
    ! Initialise multigrid
    call linsol_initMultigrid2 (rsystem%p_rsolver, ilvl-ilcrs+1, rsystem%p_RfilterChain)
    
    ! Use BiCGStab-Vanka as a coarse grid solver
    call linsol_initVANKA (p_rsmoother)
    call linsol_initBiCGStab (p_rsolver,p_rsmoother,rsystem%p_RfilterChain)
    
    ! The coarse grid in multigrid is always grid 1!
    call linsol_getMultigrid2Level (rsystem%p_rsolver,1,p_rlevelInfo)
    p_rlevelInfo%p_rcoarseGridSolver => p_rsolver

    ! Now set up the other levels...
    do i = ilcrs+1, ilvl
    
      ! Set up the VANKA smoother.
      call linsol_initVANKA (p_rsmoother)
      
      ! We will use 4 smoothing steps with damping parameter 0.7
      call linsol_convertToSmoother(p_rsmoother, nsmoothSteps, ddamping)
      
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (rsystem%p_rsolver, i-ilcrs+1, p_rlevelInfo)
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother

      ! Attach our user-defined projection to the level.
      call linsol_initProjMultigrid2Level(p_rlevelInfo, rsystem%p_rproblem%Rlevels(i)%rprojection)

    end do
    
    ! Set the output level of the solver to 2 for some output
    call parlst_getvalue_int(rsystem%p_rproblem%p_rparam, '', 'IOUTPUT', rsystem%p_rsolver%ioutputLevel, 0)
    rsystem%p_rsolver%nminIterations = 1
    rsystem%p_rsolver%nmaxIterations = 1
    
    ! parse iteration control parameters
    rsystem%riter%nminIterations = 0
    rsystem%riter%nmaxIterations = 1000
    rsystem%riter%dtolRel = 1E-8_DP
    rsystem%riter%dtolAbs = 1E-11_DP
    rsystem%riter%dstagRate = 0.95_DP
    rsystem%riter%nstagIter = 3
    call itc_getParamsFromParlist(rsystem%riter, "ITERCTRL", rsystem%p_rproblem%p_rparam)

    ! Attach the system matrix to the solver.
    allocate(Rmatrices(ilcrs:ilvl))
    do i = ilcrs, ilvl
      call lsysbl_duplicateMatrix (rsystem%p_rproblem%Rlevels(i)%rmatSys,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    call linsol_setMatrices(rsystem%p_rsolver,Rmatrices(ilcrs:ilvl))
    do i=ilcrs,ilvl
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do
    deallocate(Rmatrices)
    
  end subroutine

  ! ***********************************************************************************************
  
  subroutine stdbg_solve(rsystem)
  type(t_system), intent(inout) :: rsystem
    
  type(t_level), pointer :: p_rlvl
  type(t_matrixBlock), pointer :: p_rmatSys
  integer :: ierror, ilogRes, i, ncomp
  real(DP) :: ddef, ddefVelo, ddefPres, ddefDiv
  integer, dimension(4) :: Cnorms = LINALG_NORMEUCLID
  real(DP), dimension(4) :: Dnorms
  
    ! Initialise solver
    !call output_line('Initialising solver...')
    call linsol_initStructure (rsystem%p_rsolver, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Failed to initialise solver structure!", OU_CLASS_ERROR)
      call sys_halt()
    end if
    call linsol_initData (rsystem%p_rsolver, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Failed to initialise solver data!", OU_CLASS_ERROR)
      call sys_halt()
    end if
    
    ! get the problem level
    p_rmatSys => rsystem%p_rproblem%Rlevels(rsystem%ilevel)%rmatSys

    call parlst_getvalue_int(rsystem%p_rproblem%p_rparam, '', 'ILOGRES', ilogRes, 0)
    
    ncomp = p_rmatSys%nblocksPerRow-1
    
    ! Solve...
    if(ilogRes .gt. 0) then
      !                                           1------------------1------------------1------------------1
      call output_line(" Iter  Defect             |f_u - A*u|_2      |D*u|_2            |f_p - B*p|_2")
      call output_line("----------------------------------------------------------------------------------")
    end if
    
    ! start iterating
    call itc_initIteration(rsystem%riter)
    do while(.true.)
    
      ! compute velocity/pressure defects
      call calcDefVelo(p_rmatSys, rsystem%rvecRhsVelo, rsystem%rvecSol, rsystem%rvecDefVelo)
      call calcDefPres(p_rmatSys, rsystem%rvecRhsPres, rsystem%rvecSol, rsystem%rvecDefPres)
      !call calcDef    (p_rmatSys, rsystem%rvecRhs    , rsystem%rvecSol, rsystem%rvecDef)
      call lsysbl_vectorLinearComb(rsystem%rvecDefVelo, rsystem%rvecDefPres, 1.0_DP, 1.0_DP, rsystem%rvecDef)
    
      ! apply filter chains
      call filter_applyFilterChainVec(rsystem%rvecDefVelo, rsystem%p_RfilterChain)
      call filter_applyFilterChainVec(rsystem%rvecDefPres, rsystem%p_RfilterChain)
      call filter_applyFilterChainVec(rsystem%rvecDef    , rsystem%p_RfilterChain)
    
      ! compute defect norms
      ddef     = lsysbl_vectorNorm(rsystem%rvecDef    , LINALG_NORMEUCLID)
      if(iand(ilogRes,1) .ne. 0) then
        call lsysbl_vectorNormBlock(rsystem%rvecDefVelo, Cnorms, Dnorms)
        ddefVelo = 0.0_DP
        do i = 1, ncomp
          ddefVelo = ddefVelo + Dnorms(i)**2
        end do
        ddefVelo = sqrt(ddefVelo)
        ddefDiv  = Dnorms(ncomp+1)
        ddefPres = lsysbl_vectorNorm(rsystem%rvecDefPres, LINALG_NORMEUCLID)
      end if

      ! push new defect
      if(rsystem%riter%cstatus .eq. ITC_STATUS_UNDEFINED) then
        call itc_initResidual(rsystem%riter, ddef)
        rsystem%p_rproblem%p_Dstat(DSTAT_DEF_I, rsystem%ilevel) = ddef
      else
        call itc_pushResidual(rsystem%riter, ddef)
      end if

      ! print current residuals
      if(iand(ilogRes,1) .ne. 0) then
        call output_line(&
          trim(sys_si(rsystem%riter%niterations, 5)) // "  " // &
          trim(sys_sdel(ddef, 12)) // " " // &
          trim(sys_sdel(ddefVelo, 12)) // " " // &
          trim(sys_sdel(ddefDiv, 12)) // " " // &
          trim(sys_sdel(ddefPres, 12)))
      end if
      
      if(rsystem%riter%cstatus .ne. ITC_STATUS_CONTINUE) exit
      
      ! apply preconditioner
      call linsol_precondDefect(rsystem%p_rsolver, rsystem%rvecDef)
      
      ! update solution vector
      call lsysbl_vectorLinearComb(rsystem%rvecDef, rsystem%rvecSol, 1.0_DP, 1.0_DP)
    
    end do
    
    ! store solution statistics
    rsystem%p_rproblem%p_Istat(ISTAT_NITER,rsystem%ilevel) = rsystem%riter%niterations
    rsystem%p_rproblem%p_Istat(ISTAT_STATUS,rsystem%ilevel) = rsystem%riter%cstatus
    rsystem%p_rproblem%p_Dstat(DSTAT_DEF_F, rsystem%ilevel) = ddef
    rsystem%p_rproblem%p_Dstat(DSTAT_CRATE, rsystem%ilevel) = itc_calcConvRateReal(rsystem%riter)

    ! print summary if desired
    if(iand(ilogRes,2) .ne. 0) then
      call output_lbrk()
      call itc_printStatistics(rsystem%riter)
    end if

    ! Release solver data and structure
    call linsol_doneData (rsystem%p_rsolver)
    call linsol_doneStructure (rsystem%p_rsolver)
    
    
  contains
    subroutine calcDefVelo(rmatSys, rvecRhs, rvecSol, rvecDef)
    type(t_matrixBlock), intent(in) :: rmatSys
    type(t_vectorBlock), intent(in) :: rvecRhs, rvecSol
    type(t_vectorBlock), intent(inout) :: rvecDef
    
    integer :: i,n
    
      n = rmatSys%nblocksPerRow-1
    
      call lsysbl_copyVector(rvecRhs, rvecDef)
      do i = 1, n
        ! d_i := b_i - A_ii * u_i
        call lsyssc_matVec(rmatSys%RmatrixBlock(i,i), rvecSol%RvectorBlock(i), &
            rvecDef%RvectorBlock(i), -1.0_DP, 1.0_DP)
        ! d_n := b_n - sum{D_i * u_i}
        call lsyssc_matVec(rmatSys%RmatrixBlock(n+1,i), rvecSol%RvectorBlock(i), &
            rvecDef%RvectorBlock(n+1), -1.0_DP, 1.0_DP)
      end do

    end subroutine

    subroutine calcDefPres(rmatSys, rvecRhs, rvecSol, rvecDef)
    type(t_matrixBlock), intent(in) :: rmatSys
    type(t_vectorBlock), intent(in) :: rvecRhs, rvecSol
    type(t_vectorBlock), intent(inout) :: rvecDef
    
    integer :: i,n
    
      n = rmatSys%nblocksPerRow-1
    
      call lsysbl_copyVector(rvecRhs, rvecDef)
      do i = 1, n
        ! d_i := b_i - B_i * p
        call lsyssc_matVec(rmatSys%RmatrixBlock(i,n+1), rvecSol%RvectorBlock(n+1), &
            rvecDef%RvectorBlock(i), -1.0_DP, 1.0_DP)
      end do

    end subroutine
    
    subroutine calcDef(rmatSys, rvecRhs, rvecSol, rvecDef)
    type(t_matrixBlock), intent(in) :: rmatSys
    type(t_vectorBlock), intent(in) :: rvecRhs, rvecSol
    type(t_vectorBlock), intent(inout) :: rvecDef
    
      call lsysbl_copyVector(rvecRhs, rvecDef)
      call lsysbl_matVec(rmatSys, rvecSol, rvecDef, -1.0_DP, 1.0_DP)

    end subroutine
    
  end subroutine

  ! ***********************************************************************************************

  subroutine stdbg_filterPressureMean(rsystem)
  type(t_system), intent(inout) :: rsystem
  
  type(t_level), pointer :: p_rlvl
  type(t_matrixScalar) :: rmatMass
  type(t_vectorScalar) :: rvecPrim, rvecDual
  real(DP), dimension(:), pointer :: p_Dwork
  real(DP) :: dtol, dint, dvol
  integer :: cinfo, ncomp, niter
  
    ! fetch the level
    p_rlvl => rsystem%p_rproblem%Rlevels(rsystem%ilevel)
    
    ! assemble a matrix structure
    ncomp = p_rlvl%rdiscr%ncomponents
    call bilf_createMatrixStructure (p_rlvl%rdiscr%RspatialDiscr(ncomp), &
        LSYSSC_MATRIX9, rmatMass)
    
    ! assemble a mass matrix
    call stdop_assembleSimpleMatrix (rmatMass,DER_FUNC,DER_FUNC,1.0_DP, .true.,p_rlvl%rcubInfo)
    
    ! assemble dual vector
    call lsyssc_createVecIndMat(rmatMass, rvecDual, .true.)
    call linf_buildSimpleVector(rvecDual, p_rlvl%rcubInfo, stdbg_aux_funcRhsOne, .true., DER_FUNC)
    call lsyssc_duplicateVector(rvecDual, rvecPrim, LSYSSC_DUP_COPY, LSYSSC_DUP_COPY)
    
    ! compute initial euclid norm
    dtol = lsyssc_vectorNorm(rvecPrim, LINALG_NORMEUCLID) * 1E-17_DP
    niter = 1000
    
    ! apply CG-SSOR quick solver
    allocate(p_Dwork(3*rvecPrim%NEQ))
    call qsol_solveCG_SSOR(rmatMass, rvecPrim, p_Dwork, cinfo, niter, dtol, 1.2_DP)
    deallocate(p_Dwork)
    
    ! okay, primal and dual vector computed, now project the pressure
    
    ! compute pressure integral and domain volume
    dint = lsyssc_scalarProduct(rvecDual, rsystem%rvecSol%RvectorBlock(ncomp))
    dvol = lsyssc_scalarProduct(rvecDual, rvecPrim)
    
    ! project pressure
    call lsyssc_vectorLinearComb(rvecPrim, rsystem%rvecSol%RvectorBlock(ncomp), -dint/dvol, 1.0_DP)
    
    !dint = lsyssc_scalarProduct(rvecDual, rsystem%rvecSol%RvectorBlock(ncomp))

    ! release all temporary data
    call lsyssc_releaseVector(rvecPrim)
    call lsyssc_releaseVector(rvecDual)
    call lsyssc_releaseMatrix(rmatMass)
  
  end subroutine
  
  ! ***********************************************************************************************
  
  subroutine stdbg_writeVTK(rsystem)
  type(t_system), intent(inout) :: rsystem
  
  type(t_ucdExport) :: rexport
  character(len=SYS_STRLEN) :: sfilename, sucddir
  real(DP), dimension(:), pointer :: p_Du1, p_Du2, p_Du3, p_Dp
  
    ! fetch VTK filename from parameter list
    call parlst_getvalue_string(rsystem%p_rproblem%p_rparam, '', 'VTKFILE', sfilename, '')
    if(sfilename .eq. '') return
    
    !call output_lbrk()
    !call output_line('Writing VTK output...')

    ! Start UCD export to VTK file:
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './ucd'
    call ucd_startVTK (rexport, UCD_FLAG_STANDARD, rsystem%p_rproblem%Rlevels(rsystem%ilevel)%rtria, &
        trim(sucddir)//'/' // trim(sfilename) // '_' // trim(sys_sil(rsystem%p_rproblem%idriver,4)) // &
        '_lvl' // trim(sys_si0l(rsystem%ilevel,3)) //'.vtk')
    
    ! Allocate temporary memory for projection
    allocate(p_Du1(rsystem%p_rproblem%Rlevels(rsystem%ilevel)%rtria%NVT))
    allocate(p_Du2(rsystem%p_rproblem%Rlevels(rsystem%ilevel)%rtria%NVT))
    allocate(p_Du3(rsystem%p_rproblem%Rlevels(rsystem%ilevel)%rtria%NVT))
    allocate(p_Dp(rsystem%p_rproblem%Rlevels(rsystem%ilevel)%rtria%NEL))
    
    ! project and write pressure
    call spdp_projectToCells(rsystem%rvecSol%RvectorBlock(rsystem%rvecSol%nblocks), p_Dp)
    call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Dp)
    
    ! Project and write velocity
    call spdp_projectToVertices (rsystem%rvecSol%RvectorBlock(1), p_Du1)
    call spdp_projectToVertices (rsystem%rvecSol%RvectorBlock(2), p_Du2)
    if(rsystem%rvecSol%nblocks .gt. 3) then
      call spdp_projectToVertices (rsystem%rvecSol%RvectorBlock(3), p_Du3)
      call ucd_addVarVertBasedVec(rexport, 'velocity', p_Du1, p_Du2, p_Du3)
    else
      call ucd_addVarVertBasedVec(rexport, 'velocity', p_Du1, p_Du2)
    end if
    
    ! Project and write velocity
    call spdp_projectToVertices (rsystem%rvecRhs%RvectorBlock(1), p_Du1)
    call spdp_projectToVertices (rsystem%rvecRhs%RvectorBlock(2), p_Du2)
    if(rsystem%rvecRhs%nblocks .gt. 3) then
      call spdp_projectToVertices (rsystem%rvecRhs%RvectorBlock(3), p_Du3)
      call ucd_addVarVertBasedVec(rexport, 'force (rhs)', p_Du1, p_Du2, p_Du3)
    else
      call ucd_addVarVertBasedVec(rexport, 'force (rhs)', p_Du1, p_Du2)
    end if
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! deallocate temporary arrays
    deallocate(p_Dp)
    deallocate(p_Du3)
    deallocate(p_Du2)
    deallocate(p_Du1)
    
  end subroutine
  
  ! ***********************************************************************************************
  
  subroutine stdbg_printMeshStatistics2D(rproblem)
  type(t_problem), intent(in) :: rproblem
  integer :: i
  
    call output_lbrk()
    call output_line("Mesh Statistics")
    call output_line("===============")
    !                 12345-1234567890-1234567890-1234567890
    call output_line("Level        NVT        NMT        NEL")
    call output_line("--------------------------------------")
    do i = rproblem%ilevelCoarse, rproblem%ilevelMax
      call output_line(&
        trim(sys_si(i,3)) // "   " // &
        trim(sys_si(rproblem%Rlevels(i)%rtria%NVT,10)) // " " // &
        trim(sys_si(rproblem%Rlevels(i)%rtria%NMT,10)) // " " // &
        trim(sys_si(rproblem%Rlevels(i)%rtria%NEL,10)))
    end do
  
  end subroutine
  
  ! ***********************************************************************************************
  
  subroutine stdbg_printMatrixStatistics(rproblem)
  type(t_problem), intent(in) :: rproblem
  integer :: i, n, nV, nP, nA, nB
  
    n = rproblem%Rlevels(rproblem%ilevelMax)%rmatSys%nblocksPerRow - 1
  
    call output_lbrk()
    call output_line("Matrix Statistics")
    call output_line("=================")
    !                 12345-1234567890-1234567890-1234567890-1234567890-1234567890-1234567890
    call output_line("Level   #DOFs(V)   #DOFs(P)      #DOFs    NNZE(A)    NNZE(B)       NNZE")
    call output_line("-----------------------------------------------------------------------")
    do i = rproblem%ilevelCoarse, rproblem%ilevelMax
      nV = rproblem%Rlevels(i)%rmatSys%RmatrixBlock(1,1)%NEQ
      nP = rproblem%Rlevels(i)%rmatSys%RmatrixBlock(1,n+1)%NCOLS
      nA = rproblem%Rlevels(i)%rmatSys%RmatrixBlock(1,1)%NA
      nB = rproblem%Rlevels(i)%rmatSys%RmatrixBlock(1,n+1)%NA
      call output_line(&
        trim(sys_si(i,3)) // "   " // &
        trim(sys_si(nV,10)) // " " // &
        trim(sys_si(nP,10)) // " " // &
        trim(sys_si(n*nV+nP,10)) // " " // &
        trim(sys_si(nA,10)) // " " // &
        trim(sys_si(nB,10)) // " " // &
        trim(sys_si(n*(nA+2*nB),10)))
    end do
  
  end subroutine
  
  ! ***********************************************************************************************
  
  subroutine stdbg_printSolverStatistics(rproblem)
  type(t_problem), intent(in) :: rproblem
  integer :: ilvl
  
    call output_lbrk()
    call output_line("Solver Statistics")
    call output_line("=================")
    !                 12345-123456--123456789012345678-123456789012345678-1234567890123-123456789012
    call output_line("Level  #ITER  Defect (initial)   Defect (final)     Conv-Rate     Status")
    call output_line("---------------------------------------------------------------------------------")
    do ilvl = rproblem%ilevelMin, rproblem%ilevelMax
      call output_line(&
        trim(sys_si(ilvl,3)) // "  " // &
        trim(sys_si(rproblem%p_Istat(ISTAT_NITER,ilvl),7)) // "  " // &
        trim(sys_sdel(rproblem%p_Dstat(DSTAT_DEF_I,ilvl),12)) // " " // &
        trim(sys_sdel(rproblem%p_Dstat(DSTAT_DEF_F,ilvl),12)) // " " // &
        trim(sys_sdel(rproblem%p_Dstat(DSTAT_CRATE,ilvl),7)) // " " // &
        trim(fstatus(rproblem%p_Istat(ISTAT_STATUS,ilvl))))
    end do
  
  contains
    character(len=SYS_STRLEN) function fstatus(cstatus)
    integer, intent(in) :: cstatus
    
      select case(cstatus)
      case (ITC_STATUS_UNDEFINED)
        fstatus = 'undefined'
      case (ITC_STATUS_CONTINUE)
        fstatus = 'continue' ! this should not happen...
      case (ITC_STATUS_CONVERGED)
        fstatus = 'converged'
      case (ITC_STATUS_DIVERGED)
        fstatus = 'diverged'
      case (ITC_STATUS_MAX_ITER)
        fstatus = 'max iter'
      case (ITC_STATUS_STAGNATED)
        fstatus = 'stagnated'
      case (ITC_STATUS_ORBITING)
        fstatus = 'orbiting'
      case (ITC_STATUS_NAN_RES)
        fstatus = 'NaN found'
      case default
        fstatus = 'unknown'
      end select
    end function
  end subroutine

  ! ***********************************************************************************************
  
  subroutine stdbg_printErrorStatistics(rproblem)
  type(t_problem), intent(in) :: rproblem
  
  real(DP), dimension(4) :: Dfactor
  integer :: ilvl,j
    
    ! print statistics
    call output_lbrk()
    call output_line("Error Summary")
    call output_line("=============")
    !                 12345-123456789012345678-123456789012345678-123456789012345678-123456789012345678
    call output_line("Level |u - u_h|_L2       |u - u_h|_H1       |div(u_h)|_L2      |p - p_h|_L2")
    call output_line("---------------------------------------------------------------------------------")
    do ilvl = rproblem%ilevelMin, rproblem%ilevelMax
      call output_line(&
        trim(sys_si(ilvl,3)) // "   " // &
        trim(sys_sdel(rproblem%p_Dstat(DSTAT_U_L2,ilvl),12)) // " " // &
        trim(sys_sdel(rproblem%p_Dstat(DSTAT_U_H1,ilvl),12)) // " " // &
        trim(sys_sdel(rproblem%p_Dstat(DSTAT_U_DIV,ilvl),12)) // " " // &
        trim(sys_sdel(rproblem%p_Dstat(DSTAT_P_L2,ilvl),12)))
    end do
    call output_lbrk()
    
    call output_line("Error Reduction Factors")
    call output_line("-----------------------")
    !                 12345-123456789012345678-123456789012345678-123456789012345678-123456789012345678
    call output_line("Level |u - u_h|_L2       |u - u_h|_H1       |div(u_h)|_L2      |p - p_h|_L2")
    call output_line("---------------------------------------------------------------------------------")
    do ilvl = rproblem%ilevelMin+1, rproblem%ilevelMax
      Dfactor = 0.0_DP
      do j = 1, 4
        if(abs(rproblem%p_Dstat(j,ilvl-1)*SYS_EPSREAL_DP) .lt. abs(rproblem%p_Dstat(j,ilvl))) &
          Dfactor(j) = rproblem%p_Dstat(j,ilvl-1) / rproblem%p_Dstat(j,ilvl)
      end do
      call output_line(&
        trim(sys_si(ilvl,3)) // "   " // &
        trim(sys_sdel(Dfactor(1),12)) // " " // &
        trim(sys_sdel(Dfactor(2),12)) // " " // &
        trim(sys_sdel(Dfactor(3),12)) // " " // &
        trim(sys_sdel(Dfactor(4),12)))
    end do
    call output_lbrk()
    
  end subroutine
  
end module

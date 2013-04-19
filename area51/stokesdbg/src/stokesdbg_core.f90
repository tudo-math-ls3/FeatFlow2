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
use discretebc
use ucd
use collection, only: t_collection
use pprocerror
use bilinearformevaluation
use stdoperators
use paramlist

use stokesdbg_aux
  
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
    
    ! problem parameters
    real(DP) :: dnu = 1.0_DP
    real(DP) :: dalpha = 1.0_DP
    real(DP) :: dbeta = 1.0_DP
    real(DP) :: dgamma = 1.0_DP
    
    ! statistics array; dimension(4,ilevelMin:ileveMax)
    ! Entries:
    ! (1,lvl) -> |u-u_h|_L2
    ! (2,lvl) -> |u-u_h|_H1
    ! (3,lvl) -> |div(u_h)|_L2
    ! (4,lvl) -> |p-p_h|_L2
    real(DP), dimension(:,:), pointer :: p_Dstat
  
  end type
  
!</typeblock>

!<typeblock>
  
  type t_system
  
    ! the level this system structure is assigned to
    integer :: ilevel = 0
    
    ! the corresponding vectors
    type(t_vectorBlock) :: rvecSol
    type(t_vectorBlock) :: rvecRhs
    type(t_vectorBlock) :: rvecTmp
    
    ! the filter chain to be used
    type(t_filterChain), dimension(:), pointer :: p_rfilterChain => null()
    
    ! the linear solver to be used
    type(t_linsolNode), pointer :: p_rsolver => null()
  
  end type

!</typeblock>
  
!</types>

contains
  
  ! ***********************************************************************************************
  
  subroutine stdbg_doneProblem(rproblem)
  type(t_problem), intent(inout) :: rproblem
  
  integer :: i
  
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
  type(t_parlist), intent(in) :: rparam

  integer :: ilmin, ilmax, ilcrs
  
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
    
    ! allocate levels
    allocate(rproblem%Rlevels(ilcrs:ilmax))
    
    ! fetch other parameters
    call parlst_getvalue_double(rparam, '', 'DNU', rproblem%dnu, 1.0_DP)
    call parlst_getvalue_double(rparam, '', 'DALPHA', rproblem%dalpha, 1.0_DP)
    call parlst_getvalue_double(rparam, '', 'DBETA', rproblem%dbeta, 1.0_DP)
    call parlst_getvalue_double(rparam, '', 'DGAMMA', rproblem%dgamma, 1.0_DP)
    
    ! allocate statistics array
    allocate(rproblem%p_Dstat(4,ilmin:ilmax))
    
  end subroutine
  
  ! ***********************************************************************************************
  
  subroutine stdbg_initTriangulation2D(rproblem, rparam, sfile)
  type(t_problem), target, intent(inout) :: rproblem
  type(t_parlist), intent(inout) :: rparam
  character(len=*), intent(in) :: sfile
  
  integer :: i, ilcrs, ilmin, ilmax, idistType
  real(DP) :: ddistFactor
  character(len=SYS_STRLEN) :: spredir
  
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './pre'

    ilcrs = rproblem%ilevelCoarse
    ilmin = rproblem%ilevelMin
    ilmax = rproblem%ilevelMax
    
    ! fetch distortion type
    call parlst_getvalue_int(rparam, '', 'DIST_TYPE', idistType, 0)
    call parlst_getvalue_double(rparam, '', 'DIST_FACTOR', ddistFactor, 0.1_DP)

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
  
  subroutine stdbg_initDiscretisation(rproblem, rparam)
  type(t_problem), target, intent(inout) :: rproblem
  type(t_parlist), intent(in) :: rparam
  
  integer(I32):: celemVelo, celemPres, ccubature
  integer :: i, j, ilcrs, ilmax, ndim
  type(t_level), pointer :: p_rlvl
  character(len=SYS_STRLEN) :: sname
  
    ilcrs = rproblem%ilevelCoarse
    ilmax = rproblem%ilevelMax
    ndim = rproblem%Rlevels(ilcrs)%rtria%ndim
    
    ! read parameters
    call parlst_getvalue_string(rparam, '', 'ELEMENT_VELOCITY', sname, '')
    celemVelo = elem_igetID(sname)
    call parlst_getvalue_string(rparam, '', 'ELEMENT_PRESSURE', sname, '')
    celemPres = elem_igetID(sname)
    call parlst_getvalue_string(rparam, '', 'CUBATURE', sname, '')
    ccubature = cub_igetID(sname)
    
    do i = ilcrs, ilmax
    
      p_rlvl => rproblem%Rlevels(i)

      ! Create the block discretisation
      if(ndim .eq. 2) then
        call spdiscr_initBlockDiscr (p_rlvl%rdiscr, ndim+1, p_rlvl%rtria, rproblem%rbnd)
      else
        call spdiscr_initBlockDiscr (p_rlvl%rdiscr, ndim+1, p_rlvl%rtria)
      end if
      
      ! Set up velocity spaces
      if(ndim .eq. 2) then
        call spdiscr_initDiscr_simple (p_rlvl%rdiscr%RspatialDiscr(1), &
            celemVelo, ccubature, p_rlvl%rtria, rproblem%rbnd)
      else
        call spdiscr_initDiscr_simple (p_rlvl%rdiscr%RspatialDiscr(1), &
            celemVelo, ccubature, p_rlvl%rtria)
      end if

      do j = 2, ndim
        call spdiscr_duplicateDiscrSc(p_rlvl%rdiscr%RspatialDiscr(1), p_rlvl%rdiscr%RspatialDiscr(j))
      end do

      ! Set up pressure space
      call spdiscr_deriveSimpleDiscrSc (p_rlvl%rdiscr%RspatialDiscr(1), &
          celemPres, ccubature, p_rlvl%rdiscr%RspatialDiscr(ndim+1))
    
      ! Set up cubature info structure
      call spdiscr_createDefCubStructure(p_rlvl%rdiscr%RspatialDiscr(1), &
          p_rlvl%rcubInfo, ccubature)
      
    end do
  
  end subroutine
  
  ! ***********************************************************************************************

  subroutine stdbg_initProjections(rproblem)
  type(t_problem), target, intent(inout) :: rproblem
  
  integer :: i, j, ilcrs, ilmax, ndim
  type(t_level), pointer :: p_rlvlf, p_rlvlc
  
    ilcrs = rproblem%ilevelCoarse
    ilmax = rproblem%ilevelMax
    ndim = rproblem%Rlevels(ilcrs)%rdiscr%ndimension

    ! Initialise multi-level projection for coarse level
    call output_line('Assembling multilevel projections...')
    call mlprj_initProjectionDiscr (rproblem%Rlevels(ilcrs)%rprojection, rproblem%Rlevels(ilcrs)%rdiscr)
    
    do i = ilcrs+1, ilmax
    
      p_rlvlf => rproblem%Rlevels(i)
      p_rlvlc => rproblem%Rlevels(i-1)

      ! Create prolongation matrix structures
      call mlop_create2LvlMatrixStruct(p_rlvlc%rdiscr%RspatialDiscr(1), &
          p_rlvlf%rdiscr%RspatialDiscr(1), LSYSSC_MATRIX9, p_rlvlf%rmatProlVelo)
      call mlop_create2LvlMatrixStruct(p_rlvlc%rdiscr%RspatialDiscr(ndim+1), &
          p_rlvlf%rdiscr%RspatialDiscr(ndim+1), LSYSSC_MATRIX9, p_rlvlf%rmatProlPres)

      ! Assemble prolongation matrices
      call mlop_build2LvlProlMatrix (p_rlvlc%rdiscr%RspatialDiscr(1),&
          p_rlvlf%rdiscr%RspatialDiscr(1), .true., p_rlvlf%rmatProlVelo)
      call mlop_build2LvlProlMatrix (p_rlvlc%rdiscr%RspatialDiscr(ndim+1),&
          p_rlvlf%rdiscr%RspatialDiscr(ndim+1), .true., p_rlvlf%rmatProlPres)

      ! Initialise multi-level projection
      call mlprj_initProjectionDiscr (p_rlvlf%rprojection, p_rlvlf%rdiscr)
      do j = 1, ndim
        call mlprj_initMatrixProjection(&
            p_rlvlf%rprojection%RscalarProjection(1,j), p_rlvlf%rmatProlVelo)
      end do
      call mlprj_initMatrixProjection(&
          p_rlvlf%rprojection%RscalarProjection(1,ndim+1), p_rlvlf%rmatProlPres)

    end do
    
  end subroutine
  
  ! ***********************************************************************************************

  subroutine stdbg_initMatrices(rproblem, cconstrType)
  type(t_problem), target, intent(inout) :: rproblem
  integer(I32), optional, intent(in) :: cconstrType
  
  integer :: i, j, ilcrs, ilmax, ndim
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
    
    do i = ilcrs, ilmax

      call output_line('Assembling Level ' // trim(sys_sil(i,4)) // '...')

      p_rlvl => rproblem%Rlevels(i)

      ! Create block matrix
      call lsysbl_createMatBlockByDiscr (p_rlvl%rdiscr, p_rlvl%rmatSys)

      ! Houston, we have a saddle point problem...
      p_rlvl%rmatSys%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
      
      ! Assemble A-matrix structure
      call bilf_createMatrixStructure(p_rlvl%rdiscr%RspatialDiscr(1),&
          LSYSSC_MATRIX9, p_rlvl%rmatSys%RmatrixBlock(1,1))

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
      
    end do
    
  end subroutine
  
  ! ***********************************************************************************************

  subroutine stdbg_initNoSlipBCs(rproblem)
  type(t_problem), target, intent(inout) :: rproblem
  
  type(t_boundaryRegion) :: rrgn
  integer :: i
  
    do i = rproblem%ilevelMin, rproblem%ilevelMax

      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(rproblem%Rlevels(i)%rdiscreteBC)

      call boundary_createRegion(rproblem%rbnd,1,0,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,1,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,2,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      
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
  
    do i = rproblem%ilevelMin, rproblem%ilevelMax

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
  integer :: i
  
    do i = rproblem%ilevelMin, rproblem%ilevelMax

      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(rproblem%Rlevels(i)%rdiscreteBC)

      ! Bottom edge: u = 0
      call boundary_createRegion(rproblem%rbnd,1,1,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,1,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,2,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)

      ! Top edge: u = 0
      call boundary_createRegion(rproblem%rbnd,1,3,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,1,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,2,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)

      ! Left edge: u1 = 0, u2 = par-profile
      call boundary_createRegion(rproblem%rbnd,1,4,rrgn)
      rrgn%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,1,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcZeroBC2D)
      call bcasm_newDirichletBConRealBD (rproblem%Rlevels(i)%rdiscr,1,rrgn,&
          rproblem%Rlevels(i)%rdiscreteBC, stdbg_aux_funcParProfileBC2D)
      
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
      deallocate(rsystem%p_rsolver)
    end if
    
    if(associated(rsystem%p_RfilterChain)) &
      deallocate(rsystem%p_RfilterChain)
    
    if(rsystem%rvecTmp%NEQ .gt. 0) &
      call lsysbl_releaseVector(rsystem%rvecTmp)
    if(rsystem%rvecRhs%NEQ .gt. 0) &
      call lsysbl_releaseVector(rsystem%rvecRhs)
    if(rsystem%rvecSol%NEQ .gt. 0) &
      call lsysbl_releaseVector(rsystem%rvecSol)
  
  end subroutine
  
  ! ***********************************************************************************************

  subroutine stdbg_initSystem(rproblem, rsystem, ilevel)
  type(t_problem), target, intent(in) :: rproblem
  type(t_system), intent(out) :: rsystem
  integer, intent(in) :: ilevel
!  type(t_collection), intent(inout) :: rcollection
  
  type(t_level), pointer :: p_rlvl
  
    rsystem%ilevel = ilevel
    
    p_rlvl => rproblem%Rlevels(ilevel)

    ! Create three vectors
    call lsysbl_createVecBlockIndMat (p_rlvl%rmatSys, rsystem%rvecSol, .true.)
    call lsysbl_createVecBlockIndMat (p_rlvl%rmatSys, rsystem%rvecRhs, .true.)
    call lsysbl_createVecBlockIndMat (p_rlvl%rmatSys, rsystem%rvecTmp, .true.)

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
    call lsysbl_assignDiscreteBC(rsystem%rvecSol, p_rlvl%rdiscreteBC)
    call lsysbl_assignDiscreteBC(rsystem%rvecRhs, p_rlvl%rdiscreteBC)
    call lsysbl_assignDiscreteBC(rsystem%rvecTmp, p_rlvl%rdiscreteBC)

    ! Filter solution and rhs vectors
    !call vecfil_discreteBCsol(rvecSol)
    !call vecfil_discreteBCrhs(rvecRhs)
    
  end subroutine

  ! ***********************************************************************************************
  
  subroutine stdbg_initMultigrid(rproblem, rsystem, rparam)
  type(t_problem), target, intent(in) :: rproblem
  type(t_system), target, intent(inout) :: rsystem
  type(t_parlist), intent(inout) :: rparam

  type(t_linsolNode), pointer :: p_rsmoother, p_rsolver
  type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
  type(t_matrixBlock), dimension(:), pointer :: Rmatrices
  integer :: i, ilvl, ilcrs
  integer :: nsmoothSteps
  real(DP) :: ddamping
  
    ilvl = rsystem%ilevel
    ilcrs = rproblem%ilevelCoarse
    
    ! fetch MG parameters
    call parlst_getvalue_int(rparam, '', 'SMOOTH_STEPS', nsmoothSteps, 4)
    call parlst_getvalue_double(rparam, '', 'SMOOTH_DAMP', ddamping, 1.0_DP)
    
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
      call linsol_initProjMultigrid2Level(p_rlevelInfo, rproblem%Rlevels(i)%rprojection)

    end do
    
    ! Set the output level of the solver to 2 for some output
    call parlst_getvalue_int(rparam, '', 'IOUTPUT', rsystem%p_rsolver%ioutputLevel, 0)
    call parlst_getvalue_int(rparam, '', 'MIN_ITER', rsystem%p_rsolver%nminIterations, 0)
    call parlst_getvalue_int(rparam, '', 'MAX_ITER', rsystem%p_rsolver%nmaxIterations, 1000)
    call parlst_getvalue_double(rparam, '', 'EPSREL', rsystem%p_rsolver%depsRel, 1E-8_DP)
    call parlst_getvalue_double(rparam, '', 'EPSABS', rsystem%p_rsolver%depsAbs, 1E-11_DP)

    ! Attach the system matrix to the solver.
    allocate(Rmatrices(ilcrs:ilvl))
    do i = ilcrs, ilvl
      call lsysbl_duplicateMatrix (rproblem%Rlevels(i)%rmatSys,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    call linsol_setMatrices(rsystem%p_rsolver,Rmatrices(ilcrs:ilvl))
    do i=ilcrs,ilvl
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do
    deallocate(Rmatrices)
    
  end subroutine

  ! ***********************************************************************************************
  
  subroutine stdbg_solve(rproblem, rsystem)
  type(t_problem), intent(inout) :: rproblem
  type(t_system), intent(inout) :: rsystem
    
  integer :: ierror
  
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

    ! Solve...
    !call output_lbrk()
    !call output_line('Solving...')
    call linsol_solveAdaptively (rsystem%p_rsolver, rsystem%rvecSol,rsystem%rvecRhs,rsystem%rvecTmp)

    ! Release solver data and structure
    call linsol_doneData (rsystem%p_rsolver)
    call linsol_doneStructure (rsystem%p_rsolver)
  
  end subroutine

  ! ***********************************************************************************************
  
  subroutine stdbg_writeVTK(rproblem, rsystem, rparam)
  type(t_problem), intent(inout) :: rproblem
  type(t_system), intent(inout) :: rsystem
  type(t_parlist), intent(inout) :: rparam
  
  type(t_ucdExport) :: rexport
  character(len=SYS_STRLEN) :: sfilename, sucddir
  real(DP), dimension(:), pointer :: p_Du1, p_Du2, p_Du3, p_Dp
  
    ! fetch VTK filename from parameter list
    call parlst_getvalue_string(rparam, '', 'VTKFILE', sfilename, '')
    if(sfilename .eq. '') return
    
    !call output_lbrk()
    !call output_line('Writing VTK output...')

    ! Start UCD export to VTK file:
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './ucd'
    call ucd_startVTK (rexport, UCD_FLAG_STANDARD, rproblem%Rlevels(rsystem%ilevel)%rtria, &
        trim(sucddir)//'/' // trim(sfilename) // '_' // trim(sys_sil(rproblem%idriver,4)) // &
        '_lvl' // trim(sys_si0l(rsystem%ilevel,3)) //'.vtk')
    
    ! Allocate temporary memory for projection
    allocate(p_Du1(rproblem%Rlevels(rsystem%ilevel)%rtria%NVT))
    allocate(p_Du2(rproblem%Rlevels(rsystem%ilevel)%rtria%NVT))
    allocate(p_Du3(rproblem%Rlevels(rsystem%ilevel)%rtria%NVT))
    allocate(p_Dp(rproblem%Rlevels(rsystem%ilevel)%rtria%NEL))
    
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
    call output_line("---------------")
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
    call output_line("-----------------")
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
  
  subroutine stdbg_printErrorStatistics(rproblem)
  type(t_problem), intent(in) :: rproblem
  
  real(DP), dimension(4) :: Dfactor
  integer :: ilvl,j
    
    ! print statistics
    call output_lbrk()
    call output_line("Error Summary")
    call output_line("-------------")
    !                 12345-123456789012345678-123456789012345678-123456789012345678-123456789012345678
    call output_line("Level |u - u_h|_L2       |u - u_h|_H1       |div(u_h)|_L2      |p - p_h|_L2")
    call output_line("---------------------------------------------------------------------------------")
    do ilvl = rproblem%ilevelMin, rproblem%ilevelMax
      call output_line(&
        trim(sys_si(ilvl,3)) // "   " // &
        trim(sys_sdel(rproblem%p_Dstat(1,ilvl),12)) // " " // &
        trim(sys_sdel(rproblem%p_Dstat(2,ilvl),12)) // " " // &
        trim(sys_sdel(rproblem%p_Dstat(3,ilvl),12)) // " " // &
        trim(sys_sdel(rproblem%p_Dstat(4,ilvl),12)))
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

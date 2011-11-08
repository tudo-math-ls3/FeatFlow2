!##############################################################################
!# ****************************************************************************
!# <name> elemdbg3d_test1 </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module elemdbg3d_test1

  use fsystem
  use genoutput
  use paramlist
  use storage
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use basicgeometry
  use meshregion
  use triangulation
  use element
  use spatialdiscretisation
  use transformation
  use linearsystemscalar
  use linearsystemblock
  use discretebc
  use scalarpde
  use pprocerror
  use stdoperators
  use meshmodification
  use ucd
  use spdiscprojection
  use convection
  use collection
  use sortstrategy
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  use convection
  use pprocerror
  
  use elemdbg3d_callback
  use disto3d_aux
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine elemdbg3d_1(rparam, sConfigSection, itest)
  type(t_parlist), intent(INOUT) :: rparam
  character(LEN=*), intent(IN) :: sConfigSection
  integer, intent(IN) :: itest
  
!<description>
!</description>

!</subroutine>

  type(t_triangulation) :: rtriangulation
  type(t_blockDiscretisation) :: rdiscretisation
  type(t_linearForm) :: rlinform
  type(t_matrixBlock) :: rmatrix
  type(t_vectorBlock), target :: rvecSol, rvecRhs, rvectmp
  type(t_meshregion) :: rmeshRegion
  type(t_discreteBC), target :: rdiscreteBC
  type(t_linsolNode), pointer :: p_rsolver, p_rprecond
  type(t_matrixBlock), dimension(1) :: Rmatrices
  type(t_errorScVec) :: rerror
  integer :: NLMIN,NLMAX,ierror,ilvl
  real(DP), dimension(:,:), allocatable, target :: Derror
  integer, dimension(:,:), allocatable :: Istat
  integer :: isolver, ioutput, nmaxiter,ccubature,idistType,idistLevel,idistLvl
  integer(I32) :: celement, cshape
  real(DP) :: ddist, depsRel, depsAbs, drelax, daux1, daux2, ddist2
  character(LEN=64) :: selement,scubature
  type(t_bilinearForm) :: rform
  integer :: iwritemesh, h_Ipermute
  type(t_ucdexport) :: rexport
  real(DP), dimension(:), pointer :: p_Ddata
  real(DP) :: dnu,dbeta1,dbeta2,dbeta3,dupsam,dgamma
  integer :: istabil, isolution, ifillin
  !type(t_convStreamlineDiffusion) :: rconfigSD
  type(t_jumpStabilisation) :: rconfigEOJ
  type(t_collection) :: rcollect
  
    h_Ipermute = ST_NOHANDLE

    ! Fetch minimum and maximum levels
    call parlst_getvalue_int(rparam, sConfigSection, 'NLMIN', NLMIN, -1)
    call parlst_getvalue_int(rparam, sConfigSection, 'NLMAX', NLMAX, -1)
    
    if((NLMIN .lt. 1) .or. (NLMAX .lt. NLMIN)) then
      call output_line('Invalid NLMIN/NLMAX parameters',&
           OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg3d_1')
      call sys_halt()
    end if
    
    ! Fetch mesh distortion type
    call parlst_getvalue_int(rparam, sConfigSection, 'IDISTTYPE', idistType, 0)

    ! Fetch mesh distortion level
    call parlst_getvalue_int(rparam, sConfigSection, 'IDISTLEVEL', idistLevel, 0)

    ! Fetch mesh distortion parameters
    call parlst_getvalue_double(rparam, sConfigSection, 'DMESHDISTORT', ddist, 0.1_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DMESHDISTORT2', ddist2, 0.0_DP)
    
    ! Fetch element and cubature rule
    call parlst_getvalue_string(rparam, sConfigSection, 'SELEMENT', selement, '')
    call parlst_getvalue_string(rparam, sConfigSection, 'SCUBATURE', scubature, '')
    
    ! Fetch desired solution
    call parlst_getvalue_int(rparam, sConfigSection, 'ISOLUTION', isolution, 0)
    
    ! Fetch solver type
    call parlst_getvalue_int(rparam, sConfigSection, 'ISOLVER', isolver, 0)
    
    ! Fetch solver output level
    call parlst_getvalue_int(rparam, sConfigSection, 'IOUTPUT', ioutput, 0)
    
    ! Fetch maximum number of iterations
    call parlst_getvalue_int(rparam, sConfigSection, 'NMAXITER', nmaxiter, 100)
    
    ! Fetch problem parameters
    call parlst_getvalue_double(rparam, sConfigSection, 'DNU', dnu, 1.0_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DBETA1', dbeta1, 0.0_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DBETA2', dbeta2, 0.0_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DBETA3', dbeta3, 0.0_DP)
    
    ! Get stabilisation parameters
    call parlst_getvalue_int(rparam, sConfigSection, 'ISTABIL', istabil, 0)
    call parlst_getvalue_double(rparam, sConfigSection, 'DUPSAM', dupsam, 1.0_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DGAMMA', dgamma, 0.01_DP)
    
    ! Fetch solver tolerances
    call parlst_getvalue_double(rparam, sConfigSection, 'DEPSABS', depsAbs, 1E-11_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DEPSREL', depsRel, 1E-8_DP)
    
    ! Fetch relaxation parameter for CG-SSOR
    call parlst_getvalue_double(rparam, sConfigSection, 'DRELAX', drelax, 1.2_DP)

    ! Fetch fill-in level for ILU(k) preconditioner
    call parlst_getvalue_int(rparam, sConfigSection, 'IFILLIN', ifillin, 0)

    ! Writing of the mesh
    call parlst_getvalue_int(rparam, sConfigSection, 'IWRITEMESH', iwritemesh, 0)
    
    ! Parse element and cubature
    celement = elem_igetID(selement)
    ccubature = cub_igetID(scubature)
    
    ! Get the shape of the element
    cshape = elem_igetShape(celement)
    if(cshape .ne. cub_igetShape(ccubature)) then
      call output_line('Element and cubature formula incompatible!', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg3d_1')
      call sys_halt()
    end if
    
    ! Allocate arrays
    allocate(Derror(6,NLMIN:NLMAX))
    allocate(Istat(6,NLMIN:NLMAX))
    
    call output_separator(OU_SEP_STAR)
    call output_line('ELEMENT-DEBUGGER: 3D TEST #1')
    call output_line('============================')
    ! Print out that we are going to do:
    select case(itest)
    case(301)
      call output_line('System.............: L2-projection')
    case(302)
      call output_line('System.............: Poisson')
    case(303)
      call output_line('System.............: Convection-Diffusion')
      call output_line('DNU................: ' // trim(sys_sdEP(dnu,20,12)))
      call output_line('DBETA1.............: ' // trim(sys_sdEP(dbeta1,20,12)))
      call output_line('DBETA2.............: ' // trim(sys_sdEP(dbeta2,20,12)))
      call output_line('DBETA3.............: ' // trim(sys_sdEP(dbeta3,20,12)))
      select case(istabil)
      case(0)
        call output_line('Stabilisation......: none')
      !case(1)
      !  call output_line('Stabilisation......: Streamline-Diffusion')
      !  call output_line('DUPSAM.............: ' // trim(sys_sdEP(dupsam,20,12)))
      case(2)
        call output_line('Stabilisation......: Jump-Stabilisation')
        call output_line('DGAMMA.............: ' // trim(sys_sdEP(dgamma,20,12)))
      case default
        call output_line('Invalid ISTABIL parameter', &
          OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg3d_1')
        call sys_halt()
      end select
    case default
      call output_line('Invalid ITEST parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg3d_1')
      call sys_halt()
    end select
    select case(isolution)
    case(0)
      call output_line('Solution...........: u(x,y,z) = sin(pi*x) * sin(pi*y) * sin(pi*z)')
    case(1)
      call output_line('Solution...........: u(x,y,z) = x')
    case(2)
      call output_line('Solution...........: u(x,y,z) = y')
    case(3)
      call output_line('Solution...........: u(x,y,z) = z')
    case default
      call output_line('Invalid ISOLUTION parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg3d_1')
      call sys_halt()
    end select
    select case(cshape)
    case(BGEOM_SHAPE_HEXA)
      call output_line('Coarse Mesh........: CUBE.tri')
    case(BGEOM_SHAPE_PRISM)
      call output_line('Coarse Mesh........: PRISM.tri')
    case(BGEOM_SHAPE_TETRA)
      call output_line('Tetrahedron elements are not yet supported!', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg3d_1')
      call sys_halt()
    case(BGEOM_SHAPE_PYRA)
      call output_line('Pyramid elements are not yet supported!', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg3d_1')
      call sys_halt()
    case default
      call output_line('Element is not a valid 3D element!', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg3d_1')
      call sys_halt()
    end select
    call output_line('NLMIN..............: ' // trim(sys_siL(NLMIN,4)))
    call output_line('NLMAX..............: ' // trim(sys_siL(NLMAX,4)))
    select case(idistType)
    case(0)
      ! no distortion => nothing else to print
    case(1)
      call output_line('Distortion Type....: index-based')
    case(2)
      call output_line('Distortion Type....: XY-plane-wise index-based')
    case(3)
      call output_line('Distortion Type....: XZ-plane-wise index-based')
    case(4)
      call output_line('Distortion Type....: YZ-plane-wise index-based')
    case default
      call output_line('Invalid IDISTTYPE parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg3d_1')
      call sys_halt()
    end select
    if(idistType .ne. 0) then
      call output_line('Distortion Level...: ' // trim(sys_siL(idistLevel,4)))
      call output_line('Mesh Distortion....: ' // trim(sys_sdL(ddist,8)))
    end if
    if((idistType .ge. 2) .and. (idistType .le. 4)) then
      call output_line('Mesh Distortion 2..: ' // trim(sys_sdL(ddist2,8)))
    end if
    call output_line('Element............: ' // trim(selement))
    call output_line('Cubature rule......: ' // trim(scubature))
    select case(isolver)
    case(0)
      call output_line('Solver.............: UMFPACK4')
    case(1)
      call output_line('Solver.............: CG-SSOR')
      call output_line('Relaxation.........: ' // trim(sys_sdL(drelax,8)))
    case(2)
      call output_line('Solver.............: BiCGStab-ILU(k)')
      if(ifillin .lt. 0) then
        call output_line('IFILLIN must not be less than 0', &
          OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg3d_1')
        call sys_halt()
      end if
      call output_line('Allowed Fill-In....: ' // trim(sys_siL(ifillin,4)))
    case default
      call output_line('Invalid ISOLVER parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg3d_1')
      call sys_halt()
    end select
    call output_line('Output Level.......: ' // trim(sys_siL(ioutput,4)))
    call output_line('Maximum Iter.......: ' // trim(sys_siL(nmaxiter,12)))
    call output_line('Absolute EPS.......: ' // trim(sys_sdEP(depsAbs,20,12)))
    call output_line('Relative EPS.......: ' // trim(sys_sdEP(depsRel,20,12)))
    call output_lbrk()
    
    ! Copy important parameters into quick-access arrays of the collection,
    ! these ones are needed by the callback functions.
    rcollect%IquickAccess(1) = itest
    rcollect%IquickAccess(2) = isolution
    rcollect%DquickAccess(1) = dnu
    rcollect%DquickAccess(2) = dbeta1
    rcollect%DquickAccess(3) = dbeta2
    rcollect%DquickAccess(4) = dbeta3
    
    ! Loop over all levels
    do ilvl = NLMIN, NLMAX
    
      call output_line('Processing Level ' // trim(sys_siL(ilvl,4)) // '...')

      ! Now read in the basic triangulation.
      select case(cshape)
      case(BGEOM_SHAPE_HEXA)
        call tria_readTriFile3D (rtriangulation, './pre/CUBE.tri')
      case(BGEOM_SHAPE_PRISM)
        call tria_readTriFile3D (rtriangulation, './pre/PRISM.tri')
      end select

      if(idistLevel .le. 0) then
        idistLvl = max(1,ilvl-idistLevel)
      else
        idistLvl = idistLevel
      end if
       
      if((idistLvl .le. ilvl) .and. (idistType .ne. 0)) then
      
        ! Refine up to distortion level
        call tria_quickRefine2LevelOrdering (idistLvl-1, rtriangulation)
        
        ! Distort the mesh
        select case(idistType)
        case (1)
          ! index-based distortion
          call meshmod_disturbMesh (rtriangulation, ddist)
        
        case (2)
          ! XY-plane-wise index-based distortion
          call disto3d_distortCubePlaneXY(rtriangulation, ddist, ddist2)
        
        case (3)
          ! XZ-plane-wise index-based distortion
          call disto3d_distortCubePlaneXZ(rtriangulation, ddist, ddist2)
        
        case (4)
          ! YZ-plane-wise index-based distortion
          call disto3d_distortCubePlaneYZ(rtriangulation, ddist, ddist2)
        
        end select
        
        ! Refine up to current level
        if(idistLvl .lt. ilvl) then
          call tria_quickRefine2LevelOrdering (ilvl-idistLvl, rtriangulation)
        end if
      
      else
      
        ! Refine up to current level
        call tria_quickRefine2LevelOrdering (ilvl-1, rtriangulation)
        
      end if
            
      ! And create information about adjacencies and everything one needs from
      ! a triangulation.
      call tria_initStandardMeshFromRaw (rtriangulation)

      ! Set up discretisation
      call spdiscr_initBlockDiscr (rdiscretisation,1,rtriangulation)
      call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
          celement, ccubature, rtriangulation)
                   
      ! Create matrix structure
      call lsysbl_createMatBlockByDiscr(rdiscretisation, rmatrix)
      if((itest .eq. 303) .and. (istabil .eq. 2)) then
        ! Extended stencil for jump-stabilisation:
        ! Don't get confused here - in 3D, the construction type is also
        ! 'BILF_MATC_EDGEBASED'...
        call bilf_createMatrixStructure(rdiscretisation%RspatialDiscr(1),&
                                 LSYSSC_MATRIX9,rmatrix%RmatrixBlock(1,1),&
                                 cconstrType=BILF_MATC_EDGEBASED)
      else
        call bilf_createMatrixStructure(rdiscretisation%RspatialDiscr(1),&
                                 LSYSSC_MATRIX9,rmatrix%RmatrixBlock(1,1))
      end if
      
      ! Create vectors
      call lsysbl_createVecBlockIndMat(rmatrix, rvecSol, .true.)
      call lsysbl_createVecBlockIndMat(rmatrix, rvecRhs, .true.)
      call lsysbl_createVecBlockIndMat(rmatrix, rvecTmp, .true.)
      
      ! Set up discrete BCs
      call bcasm_initDiscreteBC(rdiscreteBC)
      call mshreg_createFromNodalProp(rmeshRegion, rtriangulation, &
                                      MSHREG_IDX_ALL)
      
      ! Describe Dirichlet BCs on that mesh region
      call bcasm_newDirichletBConMR(rdiscretisation, 1, rdiscreteBC, &
                           rmeshRegion, getBoundaryValues3D, rcollect)
      call mshreg_done(rmeshRegion)

      ! Assign BCs
      rmatrix%p_rdiscreteBC => rdiscreteBC
      rvecSol%p_rdiscreteBC => rdiscreteBC
      rvecRhs%p_rdiscreteBC => rdiscreteBC
      
      ! Store statistics
      Istat(1,ilvl) = rmatrix%NEQ
      Istat(2,ilvl) = rmatrix%rmatrixBlock(1,1)%NA
      Istat(3,ilvl) = rtriangulation%NVT
      Istat(4,ilvl) = rtriangulation%NMT
      Istat(5,ilvl) = rtriangulation%NAT
      Istat(6,ilvl) = rtriangulation%NEL

      ! Assemble system matrix
      select case(itest)
      case(301)
        ! Assemble Mass matrix
        call stdop_assembleSimpleMatrix(rmatrix%RmatrixBlock(1,1), &
                                        DER_FUNC3D, DER_FUNC3D)
      
      case(302)
        ! Assemble Laplace matrix
        call stdop_assembleLaplaceMatrix(rmatrix%RmatrixBlock(1,1))
      
      case(303)
        ! Assemble Convection-Diffusion matrix
        rform%itermcount = 6
        rform%Idescriptors(1,1) = DER_DERIV3D_X
        rform%Idescriptors(2,1) = DER_DERIV3D_X
        rform%Idescriptors(1,2) = DER_DERIV3D_Y
        rform%Idescriptors(2,2) = DER_DERIV3D_Y
        rform%Idescriptors(1,3) = DER_DERIV3D_Z
        rform%Idescriptors(2,3) = DER_DERIV3D_Z
        rform%Idescriptors(1,4) = DER_DERIV3D_X
        rform%Idescriptors(2,4) = DER_FUNC3D
        rform%Idescriptors(1,5) = DER_DERIV3D_Y
        rform%Idescriptors(2,5) = DER_FUNC3D
        rform%Idescriptors(1,6) = DER_DERIV3D_Z
        rform%Idescriptors(2,6) = DER_FUNC3D
        rform%Dcoefficients(1) = dnu
        rform%Dcoefficients(2) = dnu
        rform%Dcoefficients(3) = dnu
        rform%Dcoefficients(4) = dbeta1
        rform%Dcoefficients(5) = dbeta2
        rform%Dcoefficients(6) = dbeta3
        call bilf_buildMatrixScalar (rform,.true.,rmatrix%RmatrixBlock(1,1))
      
        if(istabil .eq. 2) then
          ! Assemble the jump stabilisation
          rconfigEOJ%dgamma = dgamma
          rconfigEOJ%dgammastar = dgamma
          rconfigEOJ%dnu = dnu
          rconfigEOJ%ccubType = CUB_G2_2D
          call conv_JumpStabilisation3d (rconfigEOJ, CONV_MODMATRIX, &
                                         rmatrix%RmatrixBlock(1,1))
        end if
      
      end select

      ! Assemble RHS vector
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC3D
      call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
            rlinform,.true.,rvecRhs%RvectorBlock(1),coeff_RHS3D,rcollect)
      
      ! In any case except for L2-projection, filter the system
      if(itest .ne. 301) then
        ! Implement BCs
        call vecfil_discreteBCrhs (rvecRhs)
        call vecfil_discreteBCsol (rvecsol)
        call matfil_discreteBC (rmatrix)
      end if
      
      ! Set up the solver
      select case(isolver)
      case (0)
        ! UMFPACK solver
        call linsol_initUMFPACK4(p_rsolver)
      
      case (1)
        ! CG-SSOR[1.2] solver
        nullify(p_rprecond)
        call linsol_initSSOR(p_rprecond,drelax)
        call linsol_initCG(p_rsolver, p_rprecond)
      
      case (2)
        ! BiCGStab-ILU(k) solver with RCMK
        nullify(p_rprecond)
        call linsol_initMILUs1x1(p_rprecond,ifillin,0.0_DP)
        call linsol_initBiCGStab(p_rsolver, p_rprecond)
        
        ! Calculate a RCMK permutation
        call sstrat_calcRevCuthillMcKee(rmatrix%RmatrixBlock(1,1), h_Ipermute)
        
        ! Permute matrix and vectors
        call lsyssc_sortMatrix(rmatrix%RmatrixBlock(1,1), .true., &
                               SSTRAT_RCM, h_Ipermute)
        call lsyssc_sortVectorInSitu(rvecSol%RvectorBlock(1), &
            rvecTmp%RvectorBlock(1), SSTRAT_RCM, h_Ipermute)
        call lsyssc_sortVectorInSitu(rvecRhs%RvectorBlock(1), &
            rvecTmp%RvectorBlock(1), SSTRAT_RCM, h_Ipermute)

      end select
      
      p_rsolver%ioutputLevel = ioutput
      p_rsolver%depsRel = depsRel
      p_rsolver%depsAbs = depsAbs
      p_rsolver%nmaxiterations = nmaxiter
      
      Rmatrices = (/rmatrix/)
      call linsol_setMatrices(p_rsolver,Rmatrices)
      call linsol_initStructure (p_rsolver, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      call linsol_initData (p_rsolver, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      
      ! Solve the system
      call linsol_solveAdaptively (p_rsolver,rvecSol,rvecRhs,rvecTmp)
      
      ! If we have a permutation, unsort the solution vector
      if(h_Ipermute .ne. ST_NOHANDLE) then
        call lsyssc_sortVectorInSitu(rvecSol%RvectorBlock(1), &
            rvecTmp%RvectorBlock(1), -SSTRAT_RCM)
      end if
      
      ! Calculate the errors to the reference function
      rerror%p_RvecCoeff => rvecSol%RvectorBlock(1:1)
      rerror%p_DerrorL2 => Derror(1:1,ilvl)
      rerror%p_DerrorH1 => Derror(2:2,ilvl)
      call pperr_scalarVec(rerror, getReferenceFunction3D, rcollect)
      
      ! Print the errors
      call output_line('Errors (L2/H1): ' // &
          trim(sys_sdEP(Derror(1,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(2,ilvl),20,12)))

      ! Probably write the mesh to disc
      if (iwritemesh .eq. 1) then
        call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                          'ucd/sol3d_'//TRIM(sys_siL(ilvl,5))//'.gmv')

        ! Project the solution to the vertices
        allocate (p_Ddata(rtriangulation%NVT))
        call spdp_projectToVertices (rvecSol%RvectorBlock(1), p_Ddata)
        call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)

        call ucd_write (rexport)
        call ucd_release (rexport)
        deallocate(p_Ddata)
      else if (iwritemesh .eq. 2) then
        call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                          'ucd/sol3d_'//TRIM(sys_siL(ilvl,5))//'.vtk')

        ! Project the solution to the vertices
        allocate (p_Ddata(rtriangulation%NVT))
        call spdp_projectToVertices (rvecSol%RvectorBlock(1), p_Ddata)
        call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)

        call ucd_write (rexport)
        call ucd_release (rexport)
        deallocate(p_Ddata)
      end if

      ! Clean up this level
      call linsol_doneData (p_rsolver)
      call linsol_doneStructure (p_rsolver)
      call linsol_releaseSolver (p_rsolver)
      call lsysbl_releaseVector (rvecTmp)
      call lsysbl_releaseVector (rvecSol)
      call lsysbl_releaseVector (rvecRhs)
      call lsysbl_releaseMatrix (rmatrix)
      if(h_Ipermute .ne. ST_NOHANDLE) call storage_free(h_Ipermute)
      call bcasm_releaseDiscreteBC (rdiscreteBC)
      call spdiscr_releaseBlockDiscr(rdiscretisation)
      call tria_done (rtriangulation)
    
    end do ! ilvl
    
    ! Print some statistics
    call output_separator(OU_SEP_MINUS)
    call output_line('Level         NEQ        NNZE         NVT' // &
                     '         NMT         NAT         NEL')
    do ilvl = NLMIN, NLMAX
      call output_line(trim(sys_si(ilvl,5)) // &
          trim(sys_si(Istat(1,ilvl),12)) // &
          trim(sys_si(Istat(2,ilvl),12)) // &
          trim(sys_si(Istat(3,ilvl),12)) // &
          trim(sys_si(Istat(4,ilvl),12)) // &
          trim(sys_si(Istat(5,ilvl),12)) // &
          trim(sys_si(Istat(6,ilvl),12)))
    end do ! ilvl

    ! Print out the L2- and H1-errors on each level
    call output_separator(OU_SEP_MINUS)
    call output_line('Level     L2-error            H1-error')
    do ilvl = NLMIN, NLMAX
      call output_line(trim(sys_si(ilvl,5)) // '   ' // &
          trim(sys_sdEP(Derror(1,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(2,ilvl),20,12)))
    end do ! ilvl
    
    ! Print out the L2- and H1- factors for each level pair
    if(NLMAX .gt. NLMIN) then
      call output_separator(OU_SEP_MINUS)
      call output_line('Level     L2-factor           H1-factor')
      do ilvl = NLMIN+1, NLMAX
        
        ! avoid division by zero here
        daux1 = 0.0_DP
        daux2 = 0.0_DP
        if(abs(Derror(1,ilvl)) .gt. SYS_EPSREAL_DP) &
          daux1 = Derror(1,ilvl-1)/Derror(1,ilvl)
        if(abs(Derror(2,ilvl)) .gt. SYS_EPSREAL_DP) &
          daux2 = Derror(2,ilvl-1)/Derror(2,ilvl)
        
        ! print out the factors
        call output_line(trim(sys_si(ilvl,5)) // '   ' // &
            trim(sys_sdEP(daux1,20,12)) // &
            trim(sys_sdEP(daux2,20,12)))
      end do ! ilvl
    end if

    ! Deallocate arrays
    deallocate(Istat)
    deallocate(Derror)
    
  end subroutine

end module

!##############################################################################
!# ****************************************************************************
!# <name> mgrenum2d_test1 </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module mgrenum2d_test1

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
  use linearformevaluation
  use bilinearformevaluation
  use cubature
  use element
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use meshregion
  use spatialdiscretisation
  use scalarpde
  use multileveloperators
  use multilevelprojection
  use stdoperators
  use sortstrategy
  use paramlist
  use random
  use matrixio
  use vectorio
  use collection, only: t_collection
  
  implicit none

  type t_level
    type(t_triangulation) :: rtria
    type(t_blockDiscretisation) :: rdisc
    type(t_matrixBlock) :: rmatSys
    type(t_matrixScalar) :: rmatProl
    type(t_interlevelProjectionBlock) :: rproj
    type(t_discreteBC) :: rdbc
    integer :: h_Isort
  end type

contains

  ! ***************************************************************************

!<subroutine>

  subroutine mgrenum2d_1(rparam)
  type(t_parlist), intent(inout) :: rparam

!</subroutine>

  character(len=256) :: sprmfile, strifile
  character(len=32) :: selement, scubature
  integer(I32) :: celement, ccubature
  integer :: isort
  type(t_level), dimension(:), pointer :: Rlvl
  type(t_boundary) :: rbnd
  type(t_vectorBlock) :: rvecSol,rvecRhs,rvecTmp
  !type(t_boundaryRegion) :: rbndRegion
  type(t_meshRegion) :: rmeshRegion
  type(t_linsolNode), pointer :: p_rsolver,p_rsmoother, p_rprec
  type(t_matrixBlock), dimension(:), pointer :: Rmatrices
  type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
  integer, dimension(:), pointer :: p_Isort
  integer :: NLMIN, NLMAX,ierror
  integer :: i,j,k,n,iseed,cavrgType,cfilterPrjMat,csmoother,csmoothType,nsmoothSteps,&
      cdumpSysMat, cdumpPrjMat, cdumpRhsVec, cdumpSolVec, ccycle, imaxIter, ndim, &
      cCGSolver, cCGPrec, iCGMaxIter
  real(DP) :: ddamping, depsRel

    ! Fetch all necessary parameters
    call parlst_getvalue_int(rparam, '', 'NDIM', ndim, 2)

    call parlst_getvalue_string(rparam, '', 'SPRMFILE', sprmfile, '')
    call parlst_getvalue_string(rparam, '', 'STRIFILE', strifile, '')

    call parlst_getvalue_int(rparam, '', 'NLMIN', NLMIN, 2)
    call parlst_getvalue_int(rparam, '', 'NLMAX', NLMAX, 5)

    call parlst_getvalue_string(rparam, '', 'SELEMENT', selement, '')
    call parlst_getvalue_string(rparam, '', 'SCUBATURE', scubature, '')

    call parlst_getvalue_int(rparam, '', 'ISORT', isort, 0)
    call parlst_getvalue_int(rparam, '', 'ISEED', iseed, RNG_DEFAULT_S)

    call parlst_getvalue_int(rparam, '', 'CAVRGTYPE', cavrgType, 1)
    call parlst_getvalue_int(rparam, '', 'CFILTERPRJMAT', cfilterPrjMat, 0)

    call parlst_getvalue_int(rparam, '', 'CCYCLE', ccycle, 0)
    call parlst_getvalue_double(rparam, '', 'DEPSREL', depsRel, 1E-8_DP)
    call parlst_getvalue_int(rparam, '', 'IMAXITER', imaxIter, 50)

    call parlst_getvalue_int(rparam, '', 'CCGSOLVER', cCGSolver, 0)
    call parlst_getvalue_int(rparam, '', 'CCGPREC', cCGPrec, 0)
    call parlst_getvalue_int(rparam, '', 'ICGMAXITER', iCGMaxIter, 50)

    call parlst_getvalue_int(rparam, '', 'CSMOOTHER', csmoother, 1)
    call parlst_getvalue_int(rparam, '', 'CSMOOTHTYPE', csmoothType, 3)
    call parlst_getvalue_int(rparam, '', 'NSMOOTHSTEPS', nsmoothSteps, 4)
    call parlst_getvalue_double(rparam, '', 'DDAMPING', ddamping, 0.7_DP)

    call parlst_getvalue_int(rparam, '', 'CDUMPSYSMAT', cdumpSysMat, 0)
    call parlst_getvalue_int(rparam, '', 'CDUMPPRJMAT', cdumpPrjMat, 0)
    call parlst_getvalue_int(rparam, '', 'CDUMPRHSVEC', cdumpRhsVec, 0)
    call parlst_getvalue_int(rparam, '', 'CDUMPSOLVEC', cdumpSolVec, 0)

    ! parse element and cubature rule ids
    celement = elem_igetID(selement)
    ccubature = cub_igetID(scubature)

    ! print parameters
    call output_line('Parsed Parameters:')
    call output_line('NDIM           = ' // trim(sys_sil(ndim,4)))
    call output_line('SPRMFILE       = ' // trim(sprmfile))
    call output_line('STRIFILE       = ' // trim(strifile))
    call output_line('NLMIN          = ' // trim(sys_sil(NLMIN,4)))
    call output_line('NLMAX          = ' // trim(sys_sil(NLMAX,4)))
    call output_line('SELEMENT       = ' // trim(elem_getName(celement))) !trim(selement))
    call output_line('SCUBATURE      = ' // trim(cub_getName(ccubature))) !trim(scubature))
    call output_line('ISORT          = ' // trim(sys_sil(isort,4)))
    call output_line('ISEED          = ' // trim(sys_sil(iseed,12)))
    call output_line('CARVGTYPE      = ' // trim(sys_sil(cavrgType,4)))
    call output_line('CFILTERPRJMAT  = ' // trim(sys_sil(cfilterPrjMat,4)))
    call output_line('CCYCLE         = ' // trim(sys_sil(ccycle,4)))
    call output_line('DEPSREL        = ' // trim(sys_sdEP(depsRel, 20, 12)))
    call output_line('IMAXITER       = ' // trim(sys_sil(imaxIter,8)))
    call output_line('CCGSOLVER      = ' // trim(sys_sil(cCGSolver,4)))
    call output_line('CCGPREC        = ' // trim(sys_sil(cCGPrec,4)))
    call output_line('ICGMAXITER     = ' // trim(sys_sil(iCGMaxIter,8)))
    call output_line('CSMOOTHER      = ' // trim(sys_sil(csmoother,4)))
    call output_line('CSMOOTHTYPE    = ' // trim(sys_sil(csmoothType,4)))
    call output_line('NSMOOTHSTEPS   = ' // trim(sys_sil(nsmoothSteps,4)))
    call output_line('DDAMPING       = ' // trim(sys_sdP(ddamping, 8, 6)))
    call output_separator(OU_SEP_STAR)

    ! Allocate memory for all levels
    allocate(Rlvl(NLMIN:NLMAX))

    call output_line('Creating meshes...')
    if(ndim .eq. 2) then

      call boundary_read_prm(rbnd, sprmfile)
      call tria_readTriFile2D (Rlvl(NLMIN)%rtria, strifile, rbnd)
      call tria_quickRefine2LevelOrdering (NLMIN-1, Rlvl(NLMIN)%rtria, rbnd)
      call tria_initStandardMeshFromRaw (Rlvl(NLMIN)%rtria, rbnd)
    
      do i = NLMIN+1, NLMAX
        call tria_refine2LevelOrdering(Rlvl(i-1)%rtria, Rlvl(i)%rtria,rbnd)
        call tria_initStandardMeshFromRaw(Rlvl(i)%rtria, rbnd)
      end do

    else if(ndim .eq. 3) then

      call tria_readTriFile3D (Rlvl(NLMIN)%rtria, strifile)
      call tria_quickRefine2LevelOrdering (NLMIN-1, Rlvl(NLMIN)%rtria)
      call tria_initStandardMeshFromRaw (Rlvl(NLMIN)%rtria)
    
      do i = NLMIN+1, NLMAX
        call tria_refine2LevelOrdering(Rlvl(i-1)%rtria, Rlvl(i)%rtria)
        call tria_initStandardMeshFromRaw(Rlvl(i)%rtria)
      end do

    else

      ! invalid dimension
      call sys_halt()

    end if

    ! create discretisations and system matrices
    call output_line('Creating discretisations...')
    if(ndim .eq. 2) then
      do i = NLMIN, NLMAX
        call spdiscr_initBlockDiscr (Rlvl(i)%rdisc, 1, Rlvl(i)%rtria, rbnd)
        call spdiscr_initDiscr_simple (Rlvl(i)%rdisc%RspatialDiscr(1), &
            celement, ccubature, Rlvl(i)%rtria, rbnd)
      end do
    else if(ndim .eq. 3) then
      do i = NLMIN, NLMAX
        call spdiscr_initBlockDiscr (Rlvl(i)%rdisc, 1, Rlvl(i)%rtria)
        call spdiscr_initDiscr_simple (Rlvl(i)%rdisc%RspatialDiscr(1), &
            celement, ccubature, Rlvl(i)%rtria)
      end do
    end if

    ! create system matrices
    call output_line('Creating system matrices...')
    do i = NLMIN, NLMAX
      call lsysbl_createMatBlockByDiscr(Rlvl(i)%rdisc,Rlvl(i)%rmatSys)
      call bilf_createMatrixStructure (Rlvl(i)%rdisc%RspatialDiscr(1),&
          LSYSSC_MATRIX9, Rlvl(i)%rmatSys%RmatrixBlock(1,1))
      call stdop_assembleLaplaceMatrix(Rlvl(i)%rmatSys%RmatrixBlock(1,1))
    end do

    call output_line('Assembling boundary conditions...')
    do i = NLMIN, NLMAX
      ! Initialise the discrete BC structure
      call bcasm_initDiscreteBC(Rlvl(i)%rdbc)

      ! create a mesh region for the whole boundary
      call mshreg_createFromNodalProp(rmeshRegion, Rlvl(i)%rtria, MSHREG_IDX_ALL)

      ! assemble BCs on the mesh region
      call bcasm_newDirichletBConMR (Rlvl(i)%rdisc, 1, Rlvl(i)%rdbc,&
          rmeshRegion, funcBCZeroMR)

      ! release the mesh region
      call mshreg_done(rmeshRegion)

      ! attach BC to matrix, filter it and detach BC structure again
      Rlvl(i)%rmatSys%p_rdiscreteBC => Rlvl(i)%rdbc
      call matfil_discreteBC (Rlvl(i)%rmatSys)
      Rlvl(i)%rmatSys%p_rdiscreteBC => null()

    end do

    ! print some statistics
    call output_separator(OU_SEP_STAR)
    call output_line('System Statistics:')
    call output_line('------------------')
                    !    012345678901-12345678901-12345678901-12345678901-12345678901-12345678901234
    call output_line('LVL         NVT         NMT         NAT         NEL         NEQ           NNZE')
    do i = NLMIN, NLMAX
      call output_line(trim(sys_si(i,3)) // &
        trim(sys_si(Rlvl(i)%rtria%NVT,12)) // &
        trim(sys_si(Rlvl(i)%rtria%NMT,12)) // &
        trim(sys_si(Rlvl(i)%rtria%NAT,12)) // &
        trim(sys_si(Rlvl(i)%rtria%NEL,12)) // &
        trim(sys_si(Rlvl(i)%rmatSys%RmatrixBlock(1,1)%NEQ,12)) // &
        trim(sys_si(Rlvl(i)%rmatSys%RmatrixBlock(1,1)%NA,15)))
    end do
    call output_separator(OU_SEP_STAR)

    ! create vectors
    call lsysbl_createVecBlockIndMat (Rlvl(NLMAX)%rmatSys,rvecSol, .true.)
    call lsysbl_createVecBlockIndMat (Rlvl(NLMAX)%rmatSys,rvecRhs, .true.)
    call lsysbl_createVecBlockIndMat (Rlvl(NLMAX)%rmatSys,rvecTmp, .true.)

    ! assemble rhs vector
    call linf_buildSimpleVector(rvecRhs%RvectorBlock(1), fcoeff_buildVectorSc_sim=coeffOne)

    ! filter vectors
    rvecRhs%p_rdiscreteBC => Rlvl(NLMAX)%rdbc
    rvecSol%p_rdiscreteBC => Rlvl(NLMAX)%rdbc
    call vecfil_discreteBCrhs (rvecRhs)
    call vecfil_discreteBCsol (rvecSol)
    rvecRhs%p_rdiscreteBC => null()
    rvecSol%p_rdiscreteBC => null()

    ! Loop over all levels except for the coarse-most one
    do i = NLMIN+1, NLMAX
    
      ! Create the matrix structure of the prolongation matrix.
      call mlop_create2LvlMatrixStruct(Rlvl(i-1)%rdisc%RspatialDiscr(1),&
          Rlvl(i)%rdisc%RspatialDiscr(1), LSYSSC_MATRIX9, Rlvl(i)%rmatProl)
      
      ! And assemble the entries of the prolongation matrix.
      call mlop_build2LvlProlMatrix (Rlvl(i-1)%rdisc%RspatialDiscr(1),&
          Rlvl(i)%rdisc%RspatialDiscr(1),.true., Rlvl(i)%rmatProl, cavrgType)

      ! filter prolongation matrix
      if(cfilterPrjMat .ne. 0) &
        call aux_filterProjMat(Rlvl(i)%rmatProl, Rlvl(i-1)%rdbc, Rlvl(i)%rdbc)

    end do

    if(isort .gt. 0) then

      ! calculate sort strategy
      do i = NLMIN, NLMAX

        ! allocate storage
        n = Rlvl(i)%RmatSys%NEQ
        call storage_new('mgrenum2d_1', 'h_Isort', 2*n, ST_INT, Rlvl(i)%h_Isort,ST_NEWBLOCK_NOINIT)
        call storage_getbase_int(Rlvl(i)%h_Isort, p_Isort)

        ! calculate permutation
        select case(isort)
        case (1)
          call sstrat_calcCuthillMcKee(Rlvl(i)%rmatSys%RmatrixBlock(1,1), p_Isort)

        case (2)
          call sstrat_calcXYZsorting(Rlvl(i)%rdisc%RspatialDiscr(1), p_Isort)

        case (3)
          call sstrat_calcStochastic(p_Isort)

        case (4)
          call aux_calcHierarchical(Rlvl(NLMIN:i), p_Isort)

        case (5)
          call sstrat_calcRandom(p_Isort, iseed)

        end select

        ! sort system matrix
        call lsyssc_sortMatrix(Rlvl(i)%rmatSys%RmatrixBlock(1,1), .true.,1, Rlvl(i)%h_Isort)

        ! now manually reset all information concerning the sorting in the matrix
        ! to fool the MG implementation thus avoiding that any unsorting is performed
        ! when prolongating or restricting vectors.
        Rlvl(i)%rmatSys%RmatrixBlock(1,1)%isortStrategy = 0
        Rlvl(i)%rmatSys%RmatrixBlock(1,1)%h_IsortPermutation = ST_NOHANDLE

      end do

      ! sort vectors
      call lsyssc_sortVectorInSitu(rvecSol%RvectorBlock(1), rvecTmp%RvectorBlock(1), 1, Rlvl(NLMAX)%h_Isort)
      call lsyssc_sortVectorInSitu(rvecRhs%RvectorBlock(1), rvecTmp%RvectorBlock(1), 1, Rlvl(NLMAX)%h_Isort)

      ! and remove any sorting information
      rvecSol%RvectorBlock(1)%isortStrategy = 0
      rvecSol%RvectorBlock(1)%h_IsortPermutation = ST_NOHANDLE
      rvecRhs%RvectorBlock(1)%isortStrategy = 0
      rvecRhs%RvectorBlock(1)%h_IsortPermutation = ST_NOHANDLE

    end if

    ! set up projection matrices
    do i = NLMIN+1,NLMAX
      
      ! sort projection matrix
      if(isort .gt. 0) &
        call aux_sortProjMatrix(Rlvl(i)%rmatProl, Rlvl(i-1)%h_Isort, Rlvl(i)%h_Isort)

      ! initialize projection structure
      call mlprj_initProjectionMat (Rlvl(i)%rproj, Rlvl(i)%rmatSys)
      call mlprj_initMatrixProjection(Rlvl(i)%rproj%RscalarProjection(1,1),Rlvl(i)%rmatProl)

    end do

    ! dumping
    do i = NLMIN, NLMAX

      ! dump system matrix to text file?
      if(iand(cdumpSysMat,1) .ne. 0) then
        call matio_writeMatrixHR (Rlvl(i)%rmatSys%RmatrixBlock(1,1), 'System Matrix', &
            .true., 0, 'A_s' // trim(sys_siL(isort,2)) // '_l' // trim(sys_siL(i,2)) // &
            '.txt', '(E20.12)', 1E-12_DP)
      end if

      ! dump system matrix to matlab file?
      if(iand(cdumpSysMat,2) .ne. 0) then
        call matio_spyMatrix('A_s' // trim(sys_siL(isort,2)) // '_l' // trim(sys_siL(i,2)), &
            'System Matrix',Rlvl(i)%rmatSys%RmatrixBlock(1,1),.true.)
      end if

      if(i .gt. NLMIN) then

        ! dump prolongation matrix to text file?
        if(iand(cdumpPrjMat,1) .ne. 0) then
          call matio_writeMatrixHR (Rlvl(i)%rmatProl, 'Prolongation Matrix', &
              .true., 0, 'P_s' // trim(sys_siL(isort,2)) // '_l' // trim(sys_siL(i,2)) // &
              '.txt', '(E20.12)', 1E-12_DP)
        end if

        ! dump prolongation matrix to matlab file?
        if(iand(cdumpSysMat,2) .ne. 0) then
          call matio_spyMatrix('P_s' // trim(sys_siL(isort,2)) // '_l' // trim(sys_siL(i,2)), &
              'Prolongation Matrix',Rlvl(i)%rmatProl,.true.)
        end if

      end if

    end do

    ! dump RHS vector to text file?
    if(iand(cdumpRhsVec,1) .ne. 0) then
      call vecio_writeVectorHR (rvecRhs%RvectorBlock(1), 'RHS', .false., 0, &
        'b_s' // trim(sys_siL(isort,2)) // '.txt', '(E20.12)')
    end if

    ! dump RHS vector to matlab file?
    if(iand(cdumpRhsVec,2) .ne. 0) then
      call vecio_spyVector('b_s' // trim(sys_siL(isort,2)), 'RHS', rvecRhs%RvectorBlock(1), .false.)
    end if

    ! And set up an interlevel projecton structure for the coarse-most level.
    call mlprj_initProjectionMat (Rlvl(NLMIN)%rproj, Rlvl(NLMIN)%rmatSys)

    ! Set up multigrid
    call linsol_initMultigrid2 (p_rsolver,NLMAX-NLMIN+1)
    
    ! set up a coarse grid preconditioner
    if(cCGSolver .ne. 0) then
      select case(cCGPrec)
      case (0)
        nullify(p_rprec)

      case (1)
        call linsol_initJacobi(p_rprec)

      case (2)
        call linsol_initSSOR(p_rprec)

      case (3)
        call linsol_initMILUs1x1 (p_rprec,0,0.0_DP)
      end select
    end if

    ! Set up a coarse grid solver.
    call linsol_getMultigrid2Level (p_rsolver,1,p_rlevelInfo)
    select case(cCGSolver)
    case (0)
      call linsol_initUMFPACK4 (p_rlevelInfo%p_rcoarseGridSolver)

    case (1)
      call linsol_initCG(p_rlevelInfo%p_rcoarseGridSolver, p_rprec)

    case (2)
      call linsol_initBiCGStab(p_rlevelInfo%p_rcoarseGridSolver, p_rprec)

    end select

    p_rlevelInfo%p_rcoarseGridSolver%nmaxIterations = iCGMaxIter

    ! Now set up the other levels...
    do i = NLMIN+1, NLMAX
    
      call linsol_getMultigrid2Level (p_rsolver,i-NLMIN+1,p_rlevelInfo)

      select case(csmoother)
      case (0) ! Jacobi
        call linsol_initJacobi(p_rsmoother)

      case (1) ! SOR
        call linsol_initSOR(p_rsmoother)

      case (2) ! ILU(0)
        call linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)

      end select

      call linsol_convertToSmoother(p_rsmoother, nsmoothSteps, ddamping)
      
      ! Attach smoothers
      if(iand(csmoothType,1) .ne. 0) &
        p_rlevelInfo%p_rpresmoother => p_rsmoother
      if(iand(csmoothType,2) .ne. 0) &
        p_rlevelInfo%p_rpostsmoother => p_rsmoother
      
      ! Attach our user-defined projection to the level.
      call linsol_initProjMultigrid2Level(p_rlevelInfo,Rlvl(i)%rproj)
      
    end do
    
    ! Set the output level of the solver to 2 for some output
    p_rsolver%ioutputLevel = 2

    ! set stopping criterion
    p_rsolver%depsRel = depsRel

    ! set maximum iterations
    p_rsolver%nmaxIterations = imaxIter

    ! set cycle
    p_rsolver%p_rsubnodeMultigrid2%icycle = ccycle
    
    ! Attach the system matrices to the solver.
    allocate(Rmatrices(NLMIN:NLMAX))
    do i = NLMIN, NLMAX
      call lsysbl_duplicateMatrix (Rlvl(i)%rmatSys,Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    call linsol_setMatrices(p_Rsolver,Rmatrices(NLMIN:NLMAX))
    do i=NLMIN,NLMAX
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do
    deallocate(Rmatrices)
    
    ! Solve the system
    call linsol_initStructure (p_rsolver, ierror)
    call linsol_initData (p_rsolver, ierror)
    call linsol_solveAdaptively (p_rsolver,rvecSol,rvecRhs,rvecTmp)
    call linsol_doneData (p_rsolver)
    call linsol_doneStructure (p_rsolver)


    ! dump solution vector to text file?
    if(iand(cdumpSolVec,1) .ne. 0) then
      call vecio_writeVectorHR (rvecSol%RvectorBlock(1), 'SOL', .false., 0, &
        'x_s' // trim(sys_siL(isort,2)) // '.txt', '(E20.12)')
    end if

    ! dump solution vector to matlab file?
    if(iand(cdumpSolVec,2) .ne. 0) then
      call vecio_spyVector('x_s' // trim(sys_siL(isort,2)), 'SOL', rvecSol%RvectorBlock(1), .false.)
    end if

    ! clean up the mess
    call linsol_releaseSolver (p_rsolver)
    do i = NLMAX, NLMIN+1, -1
      call mlprj_doneProjection(Rlvl(i)%rproj)
      call lsyssc_releaseMatrix (Rlvl(i)%rmatProl)
    end do
    call mlprj_doneProjection(Rlvl(NLMIN)%rproj)
    call lsysbl_releaseVector (rvecTmp)
    call lsysbl_releaseVector (rvecRhs)
    call lsysbl_releaseVector (rvecSol)
    do i = NLMAX, NLMIN, -1
      if(Rlvl(i)%h_Isort .ne. ST_NOHANDLE) &
        call storage_free(Rlvl(i)%h_Isort)
      call lsysbl_releaseMatrix (Rlvl(i)%rmatSys)
      call bcasm_releaseDiscreteBC (Rlvl(i)%rdbc)
      call spdiscr_releaseBlockDiscr(Rlvl(i)%rdisc)
      call tria_done (Rlvl(i)%rtria)
    end do
    
    deallocate(Rlvl)
    
    if(ndim .eq. 2) then
      call boundary_release (rbnd)
    end if

  end subroutine

  !************************************************************************************************

  ! callback function for RHS vector assembly
  subroutine coeffOne (rdiscretisation,rform,nelements,npointsPerElement,Dpoints, &
                       IdofsTest,rdomainIntSubset,Dcoefficients,rcollection)
  use basicgeometry
  use domainintegration
  type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
  type(t_linearForm), intent(IN) :: rform
  integer, intent(IN) :: nelements
  integer, intent(IN) :: npointsPerElement
  real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
  integer, dimension(:,:), intent(IN) :: IdofsTest
  type(t_domainIntSubset), intent(IN) :: rdomainIntSubset
  type(t_collection), intent(INOUT), optional :: rcollection
  real(DP), dimension(:,:,:), intent(OUT) :: Dcoefficients
  
    Dcoefficients = 1.0_DP

  end subroutine

  !************************************************************************************************

  ! callback function for homogene Dirichlet BC assembly (boundary-region version)
  subroutine funcBCZeroBR(Icomponents,rdiscretisation,rboundaryRegion,ielement, &
                          cinfoNeeded,iwhere,dwhere, Dvalues, rcollection)
  integer, dimension(:), intent(IN)                           :: Icomponents
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  type(t_boundaryRegion), intent(IN)                          :: rboundaryRegion
  integer(I32), intent(IN)                                    :: ielement
  integer, intent(IN)                                         :: cinfoNeeded
  integer(I32), intent(IN)                                     :: iwhere
  real(DP), intent(IN)                                        :: dwhere
  type(t_collection), intent(INOUT), optional                 :: rcollection
  real(DP), dimension(:), intent(OUT)                         :: Dvalues

    Dvalues = 0.0_DP

  end subroutine

  !************************************************************************************************

  ! callback function for homogene Dirichlet BC assembly (mesh-region version)
  subroutine funcBCZeroMR (Icomponents,rdiscretisation,rmeshRegion,cinfoNeeded,&
                            Iwhere,Dwhere,Dcoords,Dvalues,rcollection)
  integer, dimension(:), intent(in) :: Icomponents
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
  type(t_meshRegion), intent(in) :: rmeshRegion
  integer, intent(in) :: cinfoNeeded
  integer, dimension(4), intent(in) :: Iwhere
  real(DP), dimension(:), intent(in) :: Dwhere
  real(DP), dimension(:), intent(in) :: Dcoords
  type(t_collection), intent(inout), optional :: rcollection
  real(DP), dimension(:), intent(out) :: Dvalues

    Dvalues = 0.0_DP

  end subroutine

  !************************************************************************************************

  ! sorts a prolongation matrix
  subroutine aux_sortProjMatrix(rmat, h_Ipc, h_Ipf)
  type(t_matrixScalar), intent(inout) :: rmat
  integer, intent(inout) :: h_Ipc, h_Ipf

  ! local variables
  integer :: i,j,k,j1,j2
  integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Ipc, p_Ipf, p_Kld2, p_Kcol2
  real(DP), dimension(:), pointer :: p_DA, p_DA2

    ! fetch matrix arrays
    call lsyssc_getbase_Kld(rmat, p_Kld)
    call lsyssc_getbase_Kcol(rmat, p_Kcol)
    call lsyssc_getbase_double(rmat, p_DA)

    ! fetch permutation arrays
    call storage_getbase_int(h_Ipc, p_Ipc)
    call storage_getbase_int(h_Ipf, p_Ipf)

    ! make a backup of the matrix arrays
    allocate(p_Kld2(rmat%NEQ+1))
    allocate(p_Kcol2(rmat%NA))
    allocate(p_DA2(rmat%NA))
    do i = 1, rmat%NEQ+1
      p_Kld2(i) = p_Kld(i)
    end do
    do i = 1, rmat%NA
      p_Kcol2(i) = p_Kcol(i)
      p_DA2(i) = p_DA(i)
    end do

    ! now let's start permuting
    k = 1
    do i = 1, rmat%NEQ
      p_Kld(i) = k
      j1 = p_Kld2(p_Ipf(i))
      j2 = p_Kld2(p_Ipf(i)+1)-1
      do j = j1, j2
        p_Kcol(k) = p_Ipc(rmat%NCOLS+p_Kcol2(j))
        p_DA(k) = p_DA2(j)
        k = k+1
      end do ! j
    end do ! i
    p_Kld(rmat%NEQ+1) = k

    deallocate(p_DA2)
    deallocate(p_Kcol2)
    deallocate(p_Kld2)

    ! re-sort the matrix columns
    call lsyssc_sortMatrixCols(rmat)

  end subroutine

  !************************************************************************************************

  ! filters boundary condition DOF into projection matrix
  subroutine aux_filterProjMat(rmat, rdbc_c, rdbc_f)
  type(t_matrixScalar), intent(inout) :: rmat
  type(t_discreteBC), intent(in) :: rdbc_c, rdbc_f

  integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Idof, p_Iaux
  real(DP), dimension(:), pointer :: p_DA
  integer :: i,j,k,k1,k2

    ! fetch the matrix arrays
    call lsyssc_getbase_Kld(rmat, p_Kld)
    call lsyssc_getbase_Kcol(rmat, p_Kcol)
    call lsyssc_getbase_double(rmat, p_DA)

    ! loop over all fine mesh bcs
    do i = 1, rdbc_f%inumEntriesUsed
      call storage_getbase_int(rdbc_f%p_RdiscBCList(i)%rdirichletBCs%h_IdirichletDOFs, p_Idof)
      do j = 1, rdbc_f%p_RdiscBCList(i)%rdirichletBCs%nDOF
        k1 = p_Kld(p_Idof(j))
        k2 = p_Kld(p_Idof(j)+1)-1
        do k = k1, k2
          p_DA(k) = 0.0_DP
        end do
      end do
    end do

    allocate(p_Iaux(rmat%NCOLS))
    do i = 1, rmat%NCOLS
      p_Iaux(i) = 0
    end do

    do i = 1, rdbc_c%inumEntriesUsed
      call storage_getbase_int(rdbc_c%p_RdiscBCList(i)%rdirichletBCs%h_IdirichletDOFs, p_Idof)
      do j = 1, rdbc_c%p_RdiscBCList(i)%rdirichletBCs%nDOF
        p_Iaux(p_Idof(j)) = 1
      end do
    end do

    do i = 1, rmat%NA
      if(p_Iaux(p_Kcol(i)) .gt. 0) &
        p_DA(i) = 0.0_DP
    end do

    deallocate(p_Iaux)

  end subroutine

  !************************************************************************************************

  ! calculates hierarchical sorting
  subroutine aux_calcHierarchical(Rlvl, p_Isort)
  type(t_level), dimension(:), target, intent(in) :: Rlvl
  integer, dimension(:), intent(out) :: p_Isort

  type(t_spatialDiscretisation), dimension(:), pointer :: Rdisc
  integer :: i, m1, m2

    m1 = lbound(Rlvl,1)
    m2 = ubound(Rlvl,1)

    allocate(Rdisc(m1:m2))
    do i = m1, m2
      call spdiscr_duplicateDiscrSc(Rlvl(i)%rdisc%RspatialDiscr(1), Rdisc(i), .true.)
    end do

    call sstrat_calcHierarchical(Rdisc, p_Isort)

    do i = m2, m1, -1
      call spdiscr_releaseDiscr(Rdisc(i))
    end do

    deallocate(Rdisc)

  end subroutine

end module

!##############################################################################
!# ****************************************************************************
!# <name> elemdbg1d_test1 </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module elemdbg1d_test1

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
  use sortstrategybase
  use sortstrategy
  use bilinearformevaluation
  use linearformevaluation
  use linearsolver
  use convection
  use pprocerror
  use statistics
    
  use elemdbg1d_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine elemdbg1d_1(rparam, sConfigSection, itest)
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
  type(t_discreteBC), target :: rdiscreteBC
  type(t_linsolNode), pointer :: p_rsolver, p_rprecond
  type(t_matrixBlock), dimension(1) :: Rmatrices
  type(t_errorScVec) :: rerror
  integer :: NLMIN,NLMAX,ierror,ilvl,ccubature
  integer :: isolver, ioutput, nmaxiter,idistType,idistLevel,idistLvl
  integer(I32) :: celement, cprimaryelement, cshape
  real(DP), dimension(:,:), allocatable, target :: Derror
  real(DP), dimension(:,:), allocatable, target :: Dtimer
  integer, dimension(:,:), allocatable :: Istat
  real(DP) :: ddist, depsRel, depsAbs, drelax, daux1, daux2, daux3, daux4, daux5
  character(LEN=32) :: selement,sprimaryelement,scubature
  type(t_bilinearForm) :: rform
  integer :: iucd
  type(t_ucdexport) :: rexport
  real(DP) :: dnu,dbeta1,dupsam,dgamma
  integer :: istabil, isolution, ifillin
  !type(t_convStreamlineDiffusion) :: rconfigSD
  type(t_jumpStabilisation) :: rconfigEOJ
  type(t_collection) :: rcollect
  character(len=SYS_STRLEN) :: sucddir
  type(t_blockSortStrategy) :: rsortStrategy
  type(t_scalarCubatureInfo) :: rcubatureInfo
  type(t_timer) :: rtimer

    ! Fetch sytem variables
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './ucd'
    
    ! Fetch minimum and maximum levels
    call parlst_getvalue_int(rparam, sConfigSection, 'NLMIN', NLMIN, -1)
    call parlst_getvalue_int(rparam, sConfigSection, 'NLMAX', NLMAX, -1)
    
    if((NLMIN .lt. 1) .or. (NLMAX .lt. NLMIN)) then
      call output_line('Invalid NLMIN/NLMAX parameters',&
           OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg1d_1')
      call sys_halt()
    end if
    
    ! Fetch mesh distortion type
    call parlst_getvalue_int(rparam, sConfigSection, 'IDISTTYPE', idistType, 0)

    ! Fetch mesh distortion level
    call parlst_getvalue_int(rparam, sConfigSection, 'IDISTLEVEL', idistLevel, 0)

    ! Fetch mesh distortion parameter
    call parlst_getvalue_double(rparam, sConfigSection, 'DMESHDISTORT', ddist, 0.1_DP)
    
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
    
    ! Fetch solver tolerances
    call parlst_getvalue_double(rparam, sConfigSection, 'DEPSABS', depsAbs, 1E-11_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DEPSREL', depsRel, 1E-8_DP)

    ! Fetch problem parameters
    call parlst_getvalue_double(rparam, sConfigSection, 'DNU', dnu, 1.0_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DBETA1', dbeta1, 0.0_DP)
    
    ! Fetch stabilisation parameters
    call parlst_getvalue_int(rparam, sConfigSection, 'ISTABIL', istabil, 0)
    call parlst_getvalue_double(rparam, sConfigSection, 'DUPSAM', dupsam, 1.0_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DGAMMA', dgamma, 0.01_DP)
    
    ! Fetch relaxation parameter for CG-SSOR
    call parlst_getvalue_double(rparam, sConfigSection, 'DRELAX', drelax, 1.2_DP)

    ! Fetch fill-in level for ILU(k) preconditioner
    call parlst_getvalue_int(rparam, sConfigSection, 'IFILLIN', ifillin, 0)

    ! UCD export
    call parlst_getvalue_int(rparam, sConfigSection, 'IUCD', iucd, 0)
    
    ! Parse element and cubature
    celement = elem_igetID(selement)
    ccubature = cub_igetID(scubature)

    ! Get primary element
    cprimaryelement = elem_getPrimaryElement(celement)
    sprimaryelement = elem_getName(cprimaryelement)
    
    ! Get the shape of the element
    cshape = elem_igetShape(celement)
    if(cshape .ne. BGEOM_SHAPE_LINE) then
      call output_line('Element is not a valid 1D element!', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg1d_1')
      call sys_halt()
    end if
    
    ! Allocate arrays
    allocate(Derror(6,NLMIN:NLMAX))
    allocate(Dtimer(5,NLMIN:NLMAX))
    allocate(Istat(4,NLMIN:NLMAX))
    
    call output_separator(OU_SEP_STAR)
    call output_line('ELEMENT-DEBUGGER: 1D TEST #1')
    call output_line('============================')
    
    ! Print out that we are going to do:
    select case(itest)
    case(101)
      call output_line('System.............: L2-projection')
    case(102)
      call output_line('System.............: Poisson')
    case(103)
      call output_line('System.............: Convection-Diffusion')
      call output_line('DNU................: ' // trim(sys_sdEP(dnu,20,12)))
      if(isolution .eq. 2) then
        call output_line('DBETA1.............: ' // trim(sys_sdEP(1.0_DP,20,12)))
      else
        call output_line('DBETA1.............: ' // trim(sys_sdEP(dbeta1,20,12)))
      end if
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
          OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg1d_1')
        call sys_halt()
      end select
    case default
      call output_line('Invalid ITEST parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg1d_1')
      call sys_halt()
    end select
    select case(isolution)
    case(0)
      call output_line('Solution...........: u(x) = sin(pi*x)')
    case(1)
      call output_line('Solution...........: u(x) = x')
    case(2)
      ! Check if the solution is available
      if(itest .eq. 103) then
        ! solution is valid
        call output_line('Solution...........: u(x) = (exp(1/nu) -' // &
                         ' exp(x/nu)) / (exp(1/nu) + 1)')
        dbeta1 = 1.0_DP
      else
        call output_line('ISOLUTION = 2 is only valid for ITEST = 103', &
          OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg1d_1')
        call sys_halt()
      end if
    case default
      call output_line('Invalid ISOLUTION parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg1d_1')
      call sys_halt()
    end select
    call output_line('NLMIN..............: ' // trim(sys_siL(NLMIN,4)))
    call output_line('NLMAX..............: ' // trim(sys_siL(NLMAX,4)))
    select case(idistType)
    case(0)
      ! no distortion => nothing else to print
    case(1)
      call output_line('Distortion Type....: index-based')
    case default
      call output_line('Invalid IDISTTYPE parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg1d_1')
      call sys_halt()
    end select
    if(idistType .ne. 0) then
      call output_line('Distortion Level...: ' // trim(sys_siL(idistLevel,4)))
      call output_line('Mesh Distortion....: ' // trim(sys_sdL(ddist,8)))
    end if
#ifdef USE_LARGEINT
    call output_line('Element............: ' // trim(selement) // &
                     ' (ID=' // trim(sys_siL(int(celement, I64),12)) // ')')
    call output_line('Primary element....: ' // trim(sprimaryelement) // &
                     ' (ID=' // trim(sys_siL(int(cprimaryelement, I64),12)) // ')')
#else
    call output_line('Element............: ' // trim(selement) // &
                     ' (ID=' // trim(sys_siL(celement,12)) // ')')
    call output_line('Primary element....: ' // trim(sprimaryelement) // &
                     ' (ID=' // trim(sys_siL(cprimaryelement,12)) // ')')
#endif
    call output_line('Cubature rule......: ' // trim(scubature) // &
                     ' (ID=' // trim(sys_siL(ccubature,12)) // ')')
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
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg1d_1')
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

    ! Loop over all levels
    do ilvl = NLMIN, NLMAX
    
      call output_line('Processing Level ' // trim(sys_siL(ilvl,4)) // '...')

      ! Start time measurement
      call stat_startTimer(rtimer,STAT_TIMERSHORT)

      ! Now read in the basic triangulation.
      call tria_createRawTria1D(rtriangulation, 0.0_DP, 1.0_DP, 1)
      
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

      ! Stop time measurement
      call stat_stopTimer(rtimer)
      Dtimer(1,ilvl) = rtimer%delapsedReal

      ! Start time measurement
      call stat_clearTimer(rtimer)
      call stat_startTimer(rtimer,STAT_TIMERSHORT)

      ! Set up discretisation
      call spdiscr_initBlockDiscr (rdiscretisation,1,rtriangulation)
      call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
          celement, rtriangulation)

      ! Set up a cubature info structure
      call spdiscr_createDefCubStructure(&  
          rdiscretisation%RspatialDiscr(1), rcubatureInfo, int(ccubature, I32))
      
      ! Create matrix structure
      call lsysbl_createMatBlockByDiscr(rdiscretisation, rmatrix)
      if ((itest .eq. 103) .and. (istabil .eq. 2)) then
        ! Extended matrix stencil
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
      if(itest .ne. 101) then
        call bcasm_initDiscreteBC(rdiscreteBC)
        select case(isolution)
        case(0)
          ! u(x) = sin(pi*x) => homogene Dirichlet BCs
          call bcasm_newDirichletBC_1D(rdiscretisation, rdiscreteBC, 0.0_DP, 0.0_DP)
        case(1)
          ! u(x) = x => 0/1 Dirichlet BCs
          call bcasm_newDirichletBC_1D(rdiscretisation, rdiscreteBC, 0.0_DP, 1.0_DP)
        case(2)
          ! u(x) = (exp(nu)-exp(x*nu)) / (exp(nu)-1) => 1/0 Dirichlet BCs
          call bcasm_newDirichletBC_1D(rdiscretisation, rdiscreteBC, 1.0_DP, 0.0_DP)
        end select

        ! Assign BCs
        rmatrix%p_rdiscreteBC => rdiscreteBC
        rvecSol%p_rdiscreteBC => rdiscreteBC
        rvecRhs%p_rdiscreteBC => rdiscreteBC
      end if
      
      ! Store statistics
      Istat(1,ilvl) = rmatrix%NEQ
      Istat(2,ilvl) = rmatrix%rmatrixBlock(1,1)%NA
      Istat(3,ilvl) = rtriangulation%NVT
      Istat(4,ilvl) = rtriangulation%NEL

      ! Stop time measurement
      call stat_stopTimer(rtimer)
      Dtimer(2,ilvl) = rtimer%delapsedReal

      ! Start time measurement
      call stat_clearTimer(rtimer)
      call stat_startTimer(rtimer,STAT_TIMERSHORT)

      ! Assemble system
      select case(itest)
      case(101)
        ! Assemble Mass matrix
        call stdop_assembleSimpleMatrix(rmatrix%RmatrixBlock(1,1), &
                                        DER_FUNC1D, DER_FUNC1D,&
                                        rcubatureInfo=rcubatureInfo)
      
      case(102)
        ! Assemble Laplace matrix
        call stdop_assembleLaplaceMatrix(rmatrix%RmatrixBlock(1,1),&
                                         rcubatureInfo=rcubatureInfo)
      
      case(103)
        ! Assemble Laplace + Convection
        rform%itermcount = 2
        rform%Idescriptors(1,1) = DER_DERIV1D_X
        rform%Idescriptors(2,1) = DER_DERIV1D_X
        rform%Idescriptors(1,2) = DER_DERIV1D_X
        rform%Idescriptors(2,2) = DER_FUNC1D
        rform%Dcoefficients(1) = dnu
        rform%Dcoefficients(2) = dbeta1
        call bilf_buildMatrixScalar (rform,.false.,rmatrix%RmatrixBlock(1,1),&
                                     rcubatureInfo)

        if(istabil .eq. 2) then
          ! Assemble the jump stabilisation
          rconfigEOJ%dgamma = dgamma
          rconfigEOJ%dgammastar = dgamma
          rconfigEOJ%dnu = dnu
          call conv_JumpStabilisation1d (rconfigEOJ, CONV_MODMATRIX, &
                                         rmatrix%RmatrixBlock(1,1))
        end if
      
      end select
      
      ! Assemble RHS vector
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC1D
      call linf_buildVectorScalar (rlinform,.true.,rvecRhs%RvectorBlock(1),&
                                   rcubatureInfo,coeff_RHS1D,rcollect)

      ! In any case except for L2-projection, filter the system
      if(itest .ne. 101) then
        ! Implement BCs
        call vecfil_discreteBCrhs (rvecRhs)
        call vecfil_discreteBCsol (rvecsol)
        call matfil_discreteBC (rmatrix)
      end if
      
      ! Stop time measurement
      call stat_stopTimer(rtimer)
      Dtimer(3,ilvl) = rtimer%delapsedReal
      
      ! Start time measurement
      call stat_clearTimer(rtimer)
      call stat_startTimer(rtimer,STAT_TIMERSHORT)

      ! Set up the solver
      select case(isolver)
      case (0)
        ! UMFPACK solver
        call linsol_initUMFPACK4(p_rsolver)
      
      case (1)
        ! CG-SSOR solver
        nullify(p_rprecond)
        call linsol_initSSOR(p_rprecond,drelax)
        call linsol_initCG(p_rsolver, p_rprecond)
      
      case (2)
        ! BiCGStab-ILU(k) solver with RCMK
        nullify(p_rprecond)
        call linsol_initMILUs1x1(p_rprecond,ifillin,0.0_DP)
        call linsol_initBiCGStab(p_rsolver, p_rprecond)
        
        ! Calculate a RCMK permutation
        call sstrat_initBlockSorting (rsortStrategy,rdiscretisation)
        call sstrat_initRevCuthillMcKee(rsortStrategy%p_Rstrategies(1),rmatrix%RmatrixBlock(1,1))
        
        ! Attach the sorting strategy
        call lsysbl_setSortStrategy (rmatrix,rsortStrategy,rsortStrategy)
        call lsysbl_setSortStrategy (rvecSol,rsortStrategy)
        call lsysbl_setSortStrategy (rvecRhs,rsortStrategy)
        
        ! Permute matrix and vectors
        call lsysbl_sortMatrix(rmatrix,.true.)
        call lsysbl_sortVector(rvecSol,.true.,rvecTmp%RvectorBlock(1))
        call lsysbl_sortVector(rvecRhs,.true.,rvecTmp%RvectorBlock(1))
      
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
      
      ! If necessary, unsort the solution vector
      call lsyssc_sortVector(rvecSol%RvectorBlock(1), .false., rvecTmp%RvectorBlock(1))
      
      ! Stop time measurement
      call stat_stopTimer(rtimer)
      Dtimer(4,ilvl) = rtimer%delapsedReal

      ! Start time measurement
      call stat_clearTimer(rtimer)
      call stat_startTimer(rtimer,STAT_TIMERSHORT)

      ! Calculate the errors to the reference function
      rerror%p_RvecCoeff => rvecSol%RvectorBlock(1:1)
      rerror%p_DerrorL2 => Derror(1:1,ilvl)
      rerror%p_DerrorH1 => Derror(2:2,ilvl)
      call pperr_scalarVec(rerror, getReferenceFunction1D, rcollect, rcubatureInfo)
      
      ! Print the errors
      call output_line('Errors (L2/H1): ' // &
          trim(sys_sdEP(Derror(1,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(2,ilvl),20,12)))
          
      ! Do we perform UCD output?
      if (iucd .gt. 0) then

        ! What type of output?
        select case(iucd)
        case (1)
          call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
               trim(sucddir) // '/sol1d_' // trim(sys_siL(ilvl,5)) // '.gmv')

        case (2)
          call ucd_startVTK (rexport,UCD_FLAG_STANDARD,rtriangulation,&
             trim(sucddir) // '/sol1d_' // trim(sys_siL(ilvl,5)) // '.vtk')

        end select

        ! Project the solution to the vertices
        call ucd_addVectorByVertex (rexport,'sol',UCD_VAR_STANDARD, rvecSol%RvectorBlock(1))

        ! Write and release ucd
        call ucd_write (rexport)
        call ucd_release (rexport)

      end if
      
      ! Stop time measurement
      call stat_stopTimer(rtimer)
      Dtimer(5,ilvl) = rtimer%delapsedReal

      ! Clean up this level
      call linsol_doneData (p_rsolver)
      call linsol_doneStructure (p_rsolver)
      call linsol_releaseSolver (p_rsolver)
      call lsysbl_releaseVector (rvecTmp)
      call lsysbl_releaseVector (rvecSol)
      call lsysbl_releaseVector (rvecRhs)
      call lsysbl_releaseMatrix (rmatrix)
      call sstrat_doneBlockSorting (rsortStrategy)
      call bcasm_releaseDiscreteBC (rdiscreteBC)
      call spdiscr_releaseBlockDiscr(rdiscretisation)
      call tria_done (rtriangulation)
    
    end do ! ilvl
    
    ! Print some statistics
    call output_separator(OU_SEP_MINUS)
    call output_line('Level         NEQ        NNZE         NVT         NEL')
    do ilvl = NLMIN, NLMAX
      call output_line(trim(sys_si(ilvl,5)) // &
          trim(sys_si(Istat(1,ilvl),12)) // &
          trim(sys_si(Istat(2,ilvl),12)) // &
          trim(sys_si(Istat(3,ilvl),12)) // &
          trim(sys_si(Istat(4,ilvl),12)))
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

#ifdef USE_TIMER
    ! Print out the time measurements for each level
    call output_separator(OU_SEP_MINUS)
    call output_line('Level    CPU_Tria   CPU_Discr   CPU_Assem' // &
                     '   CPU_Solve    CPU_Post')
    
    do ilvl = NLMIN, NLMAX
      call output_line(trim(sys_si(ilvl,5)) // &
          trim(sys_sdEP(Dtimer(1,ilvl),12,4)) // &
          trim(sys_sdEP(Dtimer(2,ilvl),12,4)) // &
          trim(sys_sdEP(Dtimer(3,ilvl),12,4)) // &
          trim(sys_sdEP(Dtimer(4,ilvl),12,4)) // &
          trim(sys_sdEP(Dtimer(5,ilvl),12,4)))
    end do ! ilvl

    ! Print out the timing factors for each level pair
    if(NLMAX .gt. NLMIN) then
      call output_separator(OU_SEP_MINUS)
      call output_line('Level   Tria-fact  Discr-fact  Assem-fact' // &
                       '  Solve-fact   Post-fact')
      do ilvl = NLMIN+1, NLMAX
        
        ! avoid division by zero here
        daux1 = 0.0_DP
        daux2 = 0.0_DP
        daux3 = 0.0_DP
        daux4 = 0.0_DP
        daux5 = 0.0_DP
        
        if(abs(Dtimer(1,ilvl-1)) .gt. SYS_EPSREAL_DP) &
          daux1 = Dtimer(1,ilvl)/Dtimer(1,ilvl-1)
        if(abs(Dtimer(2,ilvl-1)) .gt. SYS_EPSREAL_DP) &
          daux2 = Dtimer(2,ilvl)/Dtimer(2,ilvl-1)
        if(abs(Dtimer(3,ilvl-1)) .gt. SYS_EPSREAL_DP) &
            daux3 = Dtimer(3,ilvl)/Dtimer(3,ilvl-1)
        if(abs(Dtimer(4,ilvl-1)) .gt. SYS_EPSREAL_DP) &
          daux4 = Dtimer(4,ilvl)/Dtimer(4,ilvl-1)
        if(abs(Dtimer(5,ilvl-1)) .gt. SYS_EPSREAL_DP) &
          daux5 = Dtimer(5,ilvl)/Dtimer(5,ilvl-1)

        ! print out the factors
        call output_line(trim(sys_si(ilvl,5)) // &
            trim(sys_sdEP(daux1,12,4)) // &
            trim(sys_sdEP(daux2,12,4)) // &
            trim(sys_sdEP(daux3,12,4)) // &
            trim(sys_sdEP(daux4,12,4)) // &
            trim(sys_sdEP(daux5,12,4)))
      end do ! ilvl
    end if
#endif

    ! Deallocate arrays
    deallocate(Istat)
    deallocate(Dtimer)
    deallocate(Derror)
    
  end subroutine

end module

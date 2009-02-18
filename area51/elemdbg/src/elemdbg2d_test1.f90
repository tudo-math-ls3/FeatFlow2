!##############################################################################
!# ****************************************************************************
!# <name> elemdbg2d_test1 </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module elemdbg2d_test1

  use fsystem
  use genoutput
  use paramlist
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use pprocerror
  use stdoperators
  use meshmodification
    
  use elemdbg2d_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine elemdbg2d_1(rparam, sConfigSection, itest)
  type(t_parlist), intent(INOUT) :: rparam
  character(LEN=*), intent(IN) :: sConfigSection
  integer, intent(IN) :: itest
  
!<description>
!</description>

!</subroutine>

  type(t_boundary) :: rboundary
  type(t_triangulation) :: rtriangulation
  type(t_blockDiscretisation) :: rdiscretisation
  type(t_linearForm) :: rlinform
  type(t_matrixBlock) :: rmatrix
  type(t_vectorBlock), target :: rvecSol, rvecRhs, rvectmp
  type(t_boundaryRegion) :: rboundaryRegion
  type(t_discreteBC), target :: rdiscreteBC
  type(t_linsolNode), pointer :: p_rsolver, p_rprecond
  type(t_matrixBlock), dimension(1) :: Rmatrices
  type(t_errorScVec) :: rerror
  integer :: NLMIN,NLMAX,ierror,ilvl
  real(DP), dimension(:,:), allocatable, target :: Derror
  integer, dimension(:,:), allocatable :: Istat
  integer :: isolver, ioutput, nmaxiter,ccubature
  integer(I32) :: celement, cshape
  real(DP) :: ddist, depsRel, depsAbs, drelax, daux1, daux2
  character(LEN=64) :: selement,scubature
    
    ! Fetch minimum and maximum levels
    call parlst_getvalue_int(rparam, sConfigSection, 'NLMIN', NLMIN, -1)
    call parlst_getvalue_int(rparam, sConfigSection, 'NLMAX', NLMAX, -1)
    
    if((NLMIN .lt. 1) .or. (NLMAX .lt. NLMIN)) then
      call output_line('Invalid NLMIN/NLMAX parameters',&
           OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
      call sys_halt()
    end if

    ! Fetch mesh distortion parameter
    call parlst_getvalue_double(rparam, sConfigSection, 'DMESHDISTORT', ddist, 0.0_DP)
    
    ! Fetch element and cubature rule
    call parlst_getvalue_string(rparam, sConfigSection, 'SELEMENT', selement, '')
    call parlst_getvalue_string(rparam, sConfigSection, 'SCUBATURE', scubature, '')
    
    ! Fetch solver type
    call parlst_getvalue_int(rparam, sConfigSection, 'ISOLVER', isolver, 0)
    
    ! Fetch solver output level
    call parlst_getvalue_int(rparam, sConfigSection, 'IOUTPUT', ioutput, 0)
    
    ! Fetch maximum number of iterations
    call parlst_getvalue_int(rparam, sConfigSection, 'NMAXITER', nmaxiter, 100)
    
    ! Fetch solver tolerances
    call parlst_getvalue_double(rparam, sConfigSection, 'DEPSABS', depsAbs, 1E-11_DP)
    call parlst_getvalue_double(rparam, sConfigSection, 'DEPSREL', depsRel, 1E-8_DP)
    
    ! Fetch relaxation parameter for CG-SSOR
    call parlst_getvalue_double(rparam, sConfigSection, 'DRELAX', drelax, 1.2_DP)
    
    ! Parse element and cubature
    celement = elem_igetID(selement)
    ccubature = cub_igetID(scubature)
    
    ! Get the shape of the element
    cshape = elem_igetShape(celement)
    if(cshape .ne. cub_igetShape(ccubature)) then
      call output_line('Element and cubature formula incompatible', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
      call sys_halt()
    end if
    
    ! Allocate arrays
    allocate(Derror(2,NLMIN:NLMAX))
    allocate(Istat(5,NLMIN:NLMAX))
    
    call output_separator(OU_SEP_STAR)
    call output_line('ELEMENT-DEBUGGER: 2D TEST #1')
    call output_line('============================')
    
    ! Print out that we are going to do:
    select case(itest)
    case(201)
      call output_line('System.........: L2-projection')
    case(202)
      call output_line('System.........: Poisson')
    case default
      call output_line('Invalid ITEST parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
      call sys_halt()
    end select
    select case(cshape)
    case(BGEOM_SHAPE_TRIA)
      call output_line('Coarse Mesh....: TRIA.tri')
    case(BGEOM_SHAPE_QUAD)
      call output_line('Coarse Mesh....: QUAD.tri')
    case default
      call output_line('Element is not a valid 2D element!', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
      call sys_halt()
    end select
    call output_line('NLMIN..........: ' // trim(sys_siL(NLMIN,4)))
    call output_line('NLMAX..........: ' // trim(sys_siL(NLMAX,4)))
    call output_line('Mesh-Distortion: ' // trim(sys_sdL(ddist,8)))
    call output_line('Element........: ' // trim(selement))
    call output_line('Cubature rule..: ' // trim(scubature))
    select case(isolver)
    case(0)
      call output_line('Solver.........: UMFPACK4')
    case(1)
      call output_line('Solver.........: CG-SSOR')
      call output_line('Relaxation.....: ' // trim(sys_sdL(drelax,8)))
    case(2)
      call output_line('Solver.........: BiCGStab-ILU(0)')
    case default
      call output_line('Invalid ISOLVER parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg2d_1')
      call sys_halt()
    end select
    call output_line('Output Level...: ' // trim(sys_siL(ioutput,4)))
    call output_line('Maximum Iter...: ' // trim(sys_siL(nmaxiter,12)))
    call output_line('Absolute EPS...: ' // trim(sys_sdEP(depsAbs,20,12)))
    call output_line('Relative EPS...: ' // trim(sys_sdEP(depsRel,20,12)))
    call output_lbrk()
    
    ! Read in parametrisation
    select case(cshape)
    case(BGEOM_SHAPE_TRIA)
      call boundary_read_prm(rboundary, './pre/TRIA.prm')
    case(BGEOM_SHAPE_QUAD)
      call boundary_read_prm(rboundary, './pre/QUAD.prm')
    end select
    
    ! Loop over all levels
    do ilvl = NLMIN, NLMAX
    
      call output_line('Processing Level ' // trim(sys_siL(ilvl,4)) // '...')

      ! Now read in the basic triangulation.
      select case(cshape)
      case(BGEOM_SHAPE_TRIA)
        call tria_readTriFile2D (rtriangulation, './pre/TRIA.tri', rboundary)
      case(BGEOM_SHAPE_QUAD)
        call tria_readTriFile2D (rtriangulation, './pre/QUAD.tri', rboundary)
      end select
       
      ! Refine it.
      call tria_quickRefine2LevelOrdering (ilvl-1,rtriangulation,rboundary)
      
      ! And create information about adjacencies and everything one needs from
      ! a triangulation.
      call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

      ! Optionally distort the mesh
      if(ddist .ne. 0.0_DP) &
        call meshmod_disturbMesh (rtriangulation, ddist)
      
      ! Set up discretisation
      call spdiscr_initBlockDiscr (rdiscretisation,1,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
          celement, ccubature, rtriangulation, rboundary)
                   
      ! Create matrix structure
      call lsysbl_createMatBlockByDiscr(rdiscretisation, rmatrix)
      call bilf_createMatrixStructure(rdiscretisation%RspatialDiscr(1),&
                               LSYSSC_MATRIX9,rmatrix%RmatrixBlock(1,1))
          
      ! Create vectors
      call lsysbl_createVecBlockIndMat(rmatrix, rvecSol, .true.)
      call lsysbl_createVecBlockIndMat(rmatrix, rvecRhs, .true.)
      call lsysbl_createVecBlockIndMat(rmatrix, rvecTmp, .true.)
      
      ! Set up discrete BCs
      call bcasm_initDiscreteBC(rdiscreteBC)
      call boundary_createRegion(rboundary,1,1,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
           rboundaryRegion,rdiscreteBC,getBoundaryValues2D)
      call boundary_createRegion(rboundary,1,2,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
           rboundaryRegion,rdiscreteBC,getBoundaryValues2D)
      call boundary_createRegion(rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
           rboundaryRegion,rdiscreteBC,getBoundaryValues2D)
      call boundary_createRegion(rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
           rboundaryRegion,rdiscreteBC,getBoundaryValues2D)
                               
      ! Assign BCs
      rmatrix%p_rdiscreteBC => rdiscreteBC
      rvecSol%p_rdiscreteBC => rdiscreteBC
      rvecRhs%p_rdiscreteBC => rdiscreteBC
      
      ! Store statistics
      Istat(1,ilvl) = rmatrix%NEQ
      Istat(2,ilvl) = rmatrix%rmatrixBlock(1,1)%NA
      Istat(3,ilvl) = rtriangulation%NVT
      Istat(4,ilvl) = rtriangulation%NMT
      Istat(5,ilvl) = rtriangulation%NEL

      ! Assemble system
      select case(itest)
      case(201)
        ! Assemble Mass matrix
        call stdop_assembleSimpleMatrix(rmatrix%RmatrixBlock(1,1), &
                                        DER_FUNC2D, DER_FUNC2D)
      
        ! Assemble RHS vector
        rlinform%itermCount = 1
        rlinform%Idescriptors(1) = DER_FUNC2D
        call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
              rlinform,.true.,rvecRhs%RvectorBlock(1),coeff_RHS2D_Mass)
      
      case(202)
        ! Assemble Laplace matrix
        call stdop_assembleLaplaceMatrix(rmatrix%RmatrixBlock(1,1))
        
        ! Assemble RHS vector
        rlinform%itermCount = 1
        rlinform%Idescriptors(1) = DER_FUNC2D
        call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
           rlinform,.true.,rvecRhs%RvectorBlock(1),coeff_RHS2D_Laplace)
        
        ! Implement BCs
        call vecfil_discreteBCrhs (rvecRhs)
        call vecfil_discreteBCsol (rvecsol)
        call matfil_discreteBC (rmatrix)
      
      end select
      
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
        ! BiCGStab-ILU(0) solver
        nullify(p_rprecond)
        call linsol_initMILUs1x1(p_rprecond,0,0.0_DP)
        call linsol_initBiCGStab(p_rsolver, p_rprecond)
      
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
      
      ! Calculate the errors to the reference function
      rerror%p_RvecCoeff => rvecSol%RvectorBlock(1:1)
      rerror%p_DerrorL2 => Derror(1:1,ilvl)
      rerror%p_DerrorH1 => Derror(2:2,ilvl)
      call pperr_scalarVec(rerror, getReferenceFunction2D)

      ! Clean up this level
      call linsol_doneData (p_rsolver)
      call linsol_doneStructure (p_rsolver)
      call linsol_releaseSolver (p_rsolver)
      call lsysbl_releaseVector (rvecTmp)
      call lsysbl_releaseVector (rvecSol)
      call lsysbl_releaseVector (rvecRhs)
      call lsysbl_releaseMatrix (rmatrix)
      call bcasm_releaseDiscreteBC (rdiscreteBC)
      call spdiscr_releaseBlockDiscr(rdiscretisation)
      call tria_done (rtriangulation)
    
    end do ! ilvl
    
    ! Release the domain
    call boundary_release (rboundary)
    
    ! Print some statistics
    call output_separator(OU_SEP_MINUS)
    call output_line('Level         NEQ        NNZE         NVT' // &
                     '         NMT         NEL')
    do ilvl = NLMIN, NLMAX
      call output_line(trim(sys_si(ilvl,5)) // &
          trim(sys_si(Istat(1,ilvl),12)) // &
          trim(sys_si(Istat(2,ilvl),12)) // &
          trim(sys_si(Istat(3,ilvl),12)) // &
          trim(sys_si(Istat(4,ilvl),12)) // &
          trim(sys_si(Istat(5,ilvl),12)))
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
        if(abs(Derror(1,ilvl)) .gt. SYS_EPSREAL) &
          daux1 = Derror(1,ilvl-1)/Derror(1,ilvl)
        if(abs(Derror(2,ilvl)) .gt. SYS_EPSREAL) &
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

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
  use linearsolver
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use meshregion
  use triangulation
  use spatialdiscretisation
  use pprocerror
  use stdoperators
  use meshmodification
    
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
  type(t_vectorBlock) :: rvecSol, rvecRhs, rvectmp
  type(t_discreteBC), target :: rdiscreteBC
  type(t_linsolNode), pointer :: p_rsolver, p_rprecond
  type(t_matrixBlock), dimension(1) :: Rmatrices
  integer :: NLMIN,NLMAX,ierror,ilvl,ccubature
  integer :: isolver, ioutput, nmaxiter
  integer(I32) :: celement
  real(DP), dimension(:,:), allocatable :: Derror
  integer, dimension(:,:), allocatable :: Istat
  real(DP) :: ddist, depsRel, depsAbs, drelax
  character(LEN=64) :: selement,scubature
    
    ! Fetch minimum and maximum levels
    call parlst_getvalue_int(rparam, sConfigSection, 'NLMIN', NLMIN, -1)
    call parlst_getvalue_int(rparam, sConfigSection, 'NLMAX', NLMAX, -1)
    
    if((NLMIN .lt. 1) .or. (NLMAX .lt. NLMIN)) then
      call output_line('Invalid NLMIN/NLMAX parameters',&
           OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg1d_1')
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
    
    ! Allocate arrays
    allocate(Derror(6,NLMIN:NLMAX))
    allocate(Istat(4,NLMIN:NLMAX))
    
    call output_separator(OU_SEP_STAR)
    call output_line('ELEMENT-DEBUGGER: 1D TEST #1')
    call output_line('============================')
    
    ! Print out that we are going to do:
    select case(itest)
    case(101)
      call output_line('System.........: L2-projection')
    case(102)
      call output_line('System.........: Poisson')
    case default
      call output_line('Invalid ITEST parameter', &
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg1d_1')
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
        OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg1d_1')
    end select
    call output_line('Output Level...: ' // trim(sys_siL(ioutput,4)))
    call output_line('Maximum Iter...: ' // trim(sys_siL(nmaxiter,12)))
    call output_line('Absolute EPS...: ' // trim(sys_sdEP(depsAbs,20,12)))
    call output_line('Relative EPS...: ' // trim(sys_sdEP(depsRel,20,12)))
    call output_lbrk()

    ! Loop over all levels
    do ilvl = NLMIN, NLMAX
    
      call output_line('Processing Level ' // trim(sys_siL(ilvl,4)) // '...')

      ! Now read in the basic triangulation.
      call tria_createRawTria1D(rtriangulation, 0.0_DP, 1.0_DP, 1)
      ! Refine it.
      call tria_quickRefine2LevelOrdering (ilvl-1,rtriangulation)
      
      ! And create information about adjacencies and everything one needs from
      ! a triangulation.
      call tria_initStandardMeshFromRaw (rtriangulation)

      ! Optionally distort the mesh
      if(ddist .ne. 0.0_DP) &
        call meshmod_disturbMesh (rtriangulation, ddist)
      
      ! Set up discretisation
      call spdiscr_initBlockDiscr (rdiscretisation,1,rtriangulation)
      call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
          celement, ccubature, rtriangulation)
                   
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
      call bcasm_newDirichletBC_1D(rdiscretisation, rdiscreteBC, 0.0_DP, 0.0_DP)

      ! Assign BCs
      rmatrix%p_rdiscreteBC => rdiscreteBC
      rvecSol%p_rdiscreteBC => rdiscreteBC
      rvecRhs%p_rdiscreteBC => rdiscreteBC
      
      ! Store statistics
      Istat(1,ilvl) = rmatrix%NEQ
      Istat(2,ilvl) = rmatrix%rmatrixBlock(1,1)%NA
      Istat(3,ilvl) = rtriangulation%NVT
      Istat(4,ilvl) = rtriangulation%NEL

      ! Assemble system
      select case(itest)
      case(101)
        ! Assemble Mass matrix
        call stdop_assembleSimpleMatrix(rmatrix%RmatrixBlock(1,1), &
                                        DER_FUNC1D, DER_FUNC1D)
      
        ! Assemble RHS vector
        rlinform%itermCount = 1
        rlinform%Idescriptors(1) = DER_FUNC1D
        call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
              rlinform,.true.,rvecRhs%RvectorBlock(1),coeff_RHS1D_Mass)
      
      case(102)
        ! Assemble Laplace matrix
        call stdop_assembleLaplaceMatrix(rmatrix%RmatrixBlock(1,1))
        
        ! Assemble RHS vector
        rlinform%itermCount = 1
        rlinform%Idescriptors(1) = DER_FUNC1D
        call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
           rlinform,.true.,rvecRhs%RvectorBlock(1),coeff_RHS1D_Laplace)
        
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
      
      case default
        call output_line('Unknown solver',&
             OU_CLASS_ERROR, OU_MODE_STD, 'elemdbg1d_1')
        call sys_halt()
      
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
      
      ! Calculate the error to the reference function.
      call pperr_scalar (rvecSol%RvectorBlock(1),PPERR_L2ERROR,&
                         Derror(1,ilvl), getReferenceFunction1D)

      call pperr_scalar (rvecSol%RvectorBlock(1),PPERR_H1ERROR,&
                         Derror(2,ilvl), getReferenceFunction1D)

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

    ! Print out the results
    call output_separator(OU_SEP_MINUS)
    call output_line('Level     L2-error            H1-error')
    do ilvl = NLMIN, NLMAX
      call output_line(trim(sys_si(ilvl,5)) // '   ' // &
          trim(sys_sdEP(Derror(1,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(2,ilvl),20,12)))
    end do ! ilvl
    
    call output_separator(OU_SEP_MINUS)
    call output_line('Level     L2-factor           H1-factor')
    do ilvl = NLMIN+1, NLMAX
      call output_line(trim(sys_si(ilvl,5)) // '   ' // &
          trim(sys_sdEP(Derror(1,ilvl-1)/Derror(1,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(2,ilvl-1)/Derror(2,ilvl),20,12)))
    end do ! ilvl

    ! Deallocate arrays
    deallocate(Istat)
    deallocate(Derror)
    
  end subroutine

end module

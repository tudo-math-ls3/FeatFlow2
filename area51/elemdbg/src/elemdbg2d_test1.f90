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

  subroutine elemdbg2d_1
  
!<description>
!</description>

!</subroutine>

    type(t_boundary) :: rboundary
    type(t_triangulation) :: rtriangulation
    type(t_blockDiscretisation) :: rdiscretisation
    type(t_linearForm) :: rlinform
    type(t_matrixBlock) :: rmatrix
    type(t_vectorBlock) :: rvecSolM, rvecSolL, rvecRhs, rvectmp
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_discreteBC), target :: rdiscreteBC
    type(t_linsolNode), pointer :: p_rsolver, p_rprecond
    type(t_matrixBlock), dimension(1) :: Rmatrices
    integer :: NLMIN,NLMAX,ierror,ilvl
    real(DP), dimension(:,:), allocatable :: Derror
    integer, dimension(:,:), allocatable :: Istat
    
    
    ! Set up minimum and maximum levels
    NLMIN = 2
    NLMAX = 6
    
    ! Allocate arrays
    allocate(Derror(6,NLMIN:NLMAX))
    allocate(Istat(5,NLMIN:NLMAX))
    
    call output_separator(OU_SEP_STAR)
    call output_line('ELEMENT-DEBUGGER: 2D TEST #1')
    call output_line('============================')
    
    ! Read in parametrisation
    call boundary_read_prm(rboundary, './pre/QUAD.prm')
    
    ! Loop over all levels
    do ilvl = NLMIN, NLMAX
    
      call output_line('Processing Level ' // trim(sys_siL(ilvl,4)) // '...')

      ! Now read in the basic triangulation.
      call tria_readTriFile2D (rtriangulation, './pre/QUAD.tri', rboundary)
       
      ! Refine it.
      call tria_quickRefine2LevelOrdering (ilvl-1,rtriangulation,rboundary)
      
      ! And create information about adjacencies and everything one needs from
      ! a triangulation.
      call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

      ! Optionally distort the mesh
      !call meshmod_disturbMesh (rtriangulation, 0.2_DP)
      
      ! Set up discretisation
      call spdiscr_initBlockDiscr (rdiscretisation,1,rtriangulation, rboundary)
      call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
          EL_Q1_2D, CUB_G5X5, rtriangulation, rboundary)
                   
      ! Create matrix structure
      call lsysbl_createMatBlockByDiscr(rdiscretisation, rmatrix)
      call bilf_createMatrixStructure(rdiscretisation%RspatialDiscr(1),&
                               LSYSSC_MATRIX9,rmatrix%RmatrixBlock(1,1))
      call lsyssc_allocEmptyMatrix(rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_ONE)
          
      ! Create vectors
      call lsysbl_createVecBlockIndMat(rmatrix, rvecSolM, .true.)
      call lsysbl_createVecBlockIndMat(rmatrix, rvecSolL, .true.)
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
      rvecSolM%p_rdiscreteBC => rdiscreteBC
      rvecSolL%p_rdiscreteBC => rdiscreteBC
      rvecRhs%p_rdiscreteBC => rdiscreteBC
      
      ! Store statistics
      Istat(1,ilvl) = rmatrix%NEQ
      Istat(2,ilvl) = rmatrix%rmatrixBlock(1,1)%NA
      Istat(3,ilvl) = rtriangulation%NVT
      Istat(4,ilvl) = rtriangulation%NMT
      Istat(5,ilvl) = rtriangulation%NEL
      
      ! Set up the solver
      nullify(p_rprecond)
      call linsol_initSSOR(p_rprecond,1.2_DP)
      call linsol_initCG(p_rsolver, p_rprecond)
      p_rsolver%depsRel = 1e-8_DP
      p_rsolver%nmaxiterations = 5000
      
      !call linsol_initUMFPACK4(p_rsolver)

      Rmatrices = (/rmatrix/)
      call linsol_setMatrices(p_rsolver,Rmatrices)
      call linsol_initStructure (p_rsolver, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop

      ! -------------------------------------------------------------
      ! MASS SYSTEM
      ! -------------------------------------------------------------

      ! Assemble Mass matrix
      call stdop_assembleSimpleMatrix(rmatrix%RmatrixBlock(1,1), &
                                      DER_FUNC2D, DER_FUNC2D)
      
      ! Assemble RHS vector
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC2D
      call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
            rlinform,.true.,rvecRhs%RvectorBlock(1),coeff_RHS2D_Mass)
      
      ! Implement BCs
      !call vecfil_discreteBCrhs (rvecRhs)
      !call vecfil_discreteBCsol (rvecsolM)
      !call matfil_discreteBC (rmatrix)
      
      call linsol_initData (p_rsolver, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      
      ! Solve the mass system
      call linsol_solveAdaptively (p_rsolver,rvecSolM,rvecRhs,rvecTmp)
      
      ! Release solver data
      call linsol_doneData (p_rsolver)
      
      ! Calculate the error to the reference function.
      call pperr_scalar (rvecSolM%RvectorBlock(1),PPERR_L2ERROR,&
                         Derror(1,ilvl), getReferenceFunction2D)

      call pperr_scalar (rvecSolM%RvectorBlock(1),PPERR_H1ERROR,&
                         Derror(2,ilvl), getReferenceFunction2D)

      ! -------------------------------------------------------------
      ! LAPLACE SYSTEM
      ! -------------------------------------------------------------

      ! Assemble Laplace matrix
      call stdop_assembleLaplaceMatrix(rmatrix%RmatrixBlock(1,1))
      
      ! Assemble RHS vector
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC2D
      call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
         rlinform,.true.,rvecRhs%RvectorBlock(1),coeff_RHS2D_Laplace)
      
      ! Implement BCs
      call vecfil_discreteBCrhs (rvecRhs)
      call vecfil_discreteBCsol (rvecsolL)
      call matfil_discreteBC (rmatrix)
      
      ! Set up the solver
      call linsol_initData (p_rsolver, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      
      ! Solve the mass system
      call linsol_solveAdaptively (p_rsolver,rvecSolL,rvecRhs,rvecTmp)
      
      ! Release solver data
      call linsol_doneData (p_rsolver)
      
      ! Calculate the error to the reference function.
      call pperr_scalar (rvecSolL%RvectorBlock(1),PPERR_L2ERROR,&
                         Derror(3,ilvl), getReferenceFunction2D)

      call pperr_scalar (rvecSolL%RvectorBlock(1),PPERR_H1ERROR,&
                         Derror(4,ilvl), getReferenceFunction2D)

      
      ! Calculate difference between solution of mass and laplace systems
      call lsysbl_vectorLinearComb(rvecSolM,rvecSolL,1.0_DP,-1.0_DP,rvecTmp)
      
      ! Calculate errors of difference
      call pperr_scalar (rvecTmp%RvectorBlock(1),PPERR_L2ERROR,Derror(5,ilvl))
      call pperr_scalar (rvecTmp%RvectorBlock(1),PPERR_H1ERROR,Derror(6,ilvl))

      ! Clean up this level
      call linsol_doneStructure (p_rsolver)
      call linsol_releaseSolver (p_rsolver)
      call lsysbl_releaseVector (rvecTmp)
      call lsysbl_releaseVector (rvecSolL)
      call lsysbl_releaseVector (rvecSolM)
      call lsysbl_releaseVector (rvecRhs)
      call lsysbl_releaseMatrix (rmatrix)
      call bcasm_releaseDiscreteBC (rdiscreteBC)
      call spdiscr_releaseBlockDiscr(rdiscretisation)
      call tria_done (rtriangulation)
    
    end do ! ilvl
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)
    
    ! Print some statistics
    call output_separator(OU_SEP_MINUS)
    call output_line('Level        NEQ        NNZE         NVT         NMT         NEL')
    do ilvl = NLMIN, NLMAX
      call output_line((sys_siL(ilvl,4)) // &
          trim(sys_si(Istat(1,ilvl),12)) // &
          trim(sys_si(Istat(2,ilvl),12)) // &
          trim(sys_si(Istat(3,ilvl),12)) // &
          trim(sys_si(Istat(4,ilvl),12)) // &
          trim(sys_si(Istat(5,ilvl),12)))
    end do ! ilvl

    ! Print out the results
    call output_separator(OU_SEP_MINUS)
    call output_line('Level    L2-error (Mass)     L2-error (Laplace)  L2-error (difference)')
    do ilvl = NLMIN, NLMAX
      call output_line((sys_siL(ilvl,7)) // &
          trim(sys_sdEP(Derror(1,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(3,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(5,ilvl),20,12)))
    end do ! ilvl
    
    call output_separator(OU_SEP_MINUS)
    call output_line('Level    L2-factor (Mass)    L2-factor (Laplace)')
    do ilvl = NLMIN+1, NLMAX
      call output_line((sys_siL(ilvl,7)) // &
          trim(sys_sdEP(Derror(1,ilvl-1)/Derror(1,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(3,ilvl-1)/Derror(3,ilvl),20,12)))
    end do ! ilvl

    ! Print out the results
    call output_separator(OU_SEP_MINUS)
    call output_line('Level    H1-error (Mass)     H1-error (Laplace)  H1-error (difference)')
    do ilvl = NLMIN, NLMAX
      call output_line((sys_siL(ilvl,7)) // &
          trim(sys_sdEP(Derror(2,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(4,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(6,ilvl),20,12)))
    end do ! ilvl
    
    call output_separator(OU_SEP_MINUS)
    call output_line('Level    H1-factor (Mass)    H1-factor (Laplace)')
    do ilvl = NLMIN+1, NLMAX
      call output_line((sys_siL(ilvl,7)) // &
          trim(sys_sdEP(Derror(2,ilvl-1)/Derror(2,ilvl),20,12)) // &
          trim(sys_sdEP(Derror(4,ilvl-1)/Derror(4,ilvl),20,12)))
    end do ! ilvl

    ! Deallocate arrays
    deallocate(Istat)
    deallocate(Derror)
    
  end subroutine

end module

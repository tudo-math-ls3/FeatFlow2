!##############################################################################
!# ****************************************************************************
!# <name> renum2d_method1_mg </name>
!# ****************************************************************************
!#
!# <purpose>
!# MFLOP test suite for testing the speed of a couple of solvers on
!# different levels.
!# </purpose>
!##############################################################################

MODULE renum2d_test2

  USE fsystem
  USE genoutput
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE ucd
  USE pprocerror
  USE genoutput
  USE sortstrategy
  USE statistics
    
  USE renum2d_callback
  
  IMPLICIT NONE

!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_level
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    TYPE(t_blockDiscretisation) :: rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the 
    ! discrete Laplace operator.
    TYPE(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC) :: rdiscreteBC
    
    ! Resorting strategy
    INTEGER :: h_Ipermutation
  
  END TYPE

!</typeblock>

  
!<typeblock description="Solver parameters">

  type t_solverConfig
  
    integer :: NLMIN
    integer :: NLMAX
    integer :: isortStrategy
    integer :: NLCOARSE
    integer :: isolver
    integer :: ipreconditioner
    integer :: ismoother
    integer :: nsm
    integer :: ifillILU
    integer :: ioutputLevel
    integer :: nmaxIterations
    integer :: iresNorm
    integer :: icycle
    real(dp) :: domega
    real(dp) :: depsRel
    integer :: iwriteGMV
    integer :: ielementType
  
  end type
  
!</typeblock>

!</types>

CONTAINS

  subroutine readSolverConfig (rparList,ssection,rsolverConfig)
  
!<description>
  ! Reads a solver configuration from a parameter list.
!</description>  

  type(t_parlist), intent(in) :: rparList
  character(len=*), intent(in) :: ssection
  type(t_solverConfig), intent(out) :: rsolverConfig
  
    ! Read the parameters
    call parlst_getvalue_int (rparlist,ssection,'isortStrategy',rsolverConfig%isortStrategy,0)
    call parlst_getvalue_int (rparlist,ssection,'isolver',rsolverConfig%isolver,0)
    call parlst_getvalue_int (rparlist,ssection,'ipreconditioner',rsolverConfig%ipreconditioner,0)
    call parlst_getvalue_int (rparlist,ssection,'ismoother',rsolverConfig%ismoother,0)
    call parlst_getvalue_int (rparlist,ssection,'nsm',rsolverConfig%nsm,4)
    call parlst_getvalue_int (rparlist,ssection,'ifillILU',rsolverConfig%ifillILU,0)
    call parlst_getvalue_double (rparlist,ssection,'domega',rsolverConfig%domega,1.0_DP)
    
    ! Parameters with a predefined value in the main section

    call parlst_getvalue_int (rparlist, 'PERFORMANCETESTS', 'NLMIN', &
        rsolverConfig%NLMIN, 7)
    call parlst_getvalue_int (rparlist,ssection,'NLMIN',&
        rsolverConfig%NLMIN,rsolverConfig%NLMIN)

    call parlst_getvalue_int (rparlist, 'PERFORMANCETESTS', 'NLMAX', &
        rsolverConfig%NLMAX, 11)
    call parlst_getvalue_int (rparlist,ssection,'NLMAX',&
        rsolverConfig%NLMAX,rsolverConfig%NLMAX)
    
    call parlst_getvalue_int (rparlist, 'PERFORMANCETESTS', 'ioutputLevel', &
        rsolverConfig%ioutputLevel, 0)
    call parlst_getvalue_int (rparlist,ssection,'ioutputLevel',&
        rsolverConfig%ioutputLevel,rsolverConfig%ioutputLevel)
    
    call parlst_getvalue_int (rparlist, 'PERFORMANCETESTS', 'nmaxIterations', &
        rsolverConfig%nmaxIterations, 10000)
    call parlst_getvalue_int (rparlist,ssection,'nmaxIterations',&
        rsolverConfig%nmaxIterations,rsolverConfig%nmaxIterations)
  
    call parlst_getvalue_int (rparlist, 'PERFORMANCETESTS', 'iresNorm', &
        rsolverConfig%iresNorm, 2)
    call parlst_getvalue_int (rparlist,ssection,'iresNorm',&
        rsolverConfig%iresNorm,rsolverConfig%iresNorm)
  
    call parlst_getvalue_double (rparlist, 'PERFORMANCETESTS', 'depsRel', & 
        rsolverConfig%depsRel, 1E-8_DP)
    call parlst_getvalue_double (rparlist,ssection,'iresNorm',&
        rsolverConfig%depsRel,rsolverConfig%depsRel)

    call parlst_getvalue_int (rparlist, 'PERFORMANCETESTS', 'NLCOARSE', &
        rsolverConfig%NLCOARSE, 1)
    call parlst_getvalue_int (rparlist,ssection,'NLCOARSE',&
        rsolverConfig%NLCOARSE,rsolverConfig%NLCOARSE)

    call parlst_getvalue_int (rparlist, 'PERFORMANCETESTS', 'icycle', &
        rsolverConfig%icycle, 0)
    call parlst_getvalue_int (rparlist,ssection,'icycle',&
        rsolverConfig%icycle,rsolverConfig%icycle)

    call parlst_getvalue_int (rparlist, 'PERFORMANCETESTS', 'iwriteGMV', &
        rsolverConfig%iwriteGMV, 0)
    call parlst_getvalue_int (rparlist,ssection,'iwriteGMV',&
        rsolverConfig%iwriteGMV,rsolverConfig%iwriteGMV)

    call parlst_getvalue_int (rparlist, 'PERFORMANCETESTS', 'ielementType', &
        rsolverConfig%ielementType, 0)
    call parlst_getvalue_int (rparlist,ssection,'ielementType',&
        rsolverConfig%ielementType,rsolverConfig%ielementType)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE renum2d_solvertest
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Set up RHS
  ! 4.) Set up matrix
  ! 5.) Create solver structure
  ! 6.) Solve the problem
  ! 7.) Write solution to GMV file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An array of problem levels for the multigrid solver
    TYPE(t_level), DIMENSION(:), TARGET, ALLOCATABLE :: Rlevels

    ! An object for saving the domain:
    TYPE(t_boundary) :: rboundary
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    ! A couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock,rtempBlock2

    ! A variable that is used to specify a region on the boundary.
    TYPE(t_boundaryRegion) :: rboundaryRegion

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rcoarseGridSolver,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(:), ALLOCATABLE :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock) :: rprojection

    ! One level of multigrid
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    ! Some temporary variables
    INTEGER :: i
    
    
    ! Sorting strategy counter
    INTEGER(I32), DIMENSION(:), POINTER :: p_Ipermutation
    TYPE(t_spatialDiscretisation), DIMENSION(:), ALLOCATABLE :: Rdiscretisation
    
    TYPE(t_timer) :: rtimerSolver,rtimerRef,rtimerSort,rtimerMatGen
    
    ! DAT-file parameters
    type(t_parlist) :: rparlist
    integer :: iconfig,nconfigs
    character(len=SYS_STRLEN) :: sconfig
    type(t_solverConfig) :: rsolverConfig
    integer :: ilevel, iperform
    logical :: bexists
    character(len=SYS_STRLEN) :: sprmfile,strifile,sdatafile

    ! Ok, let's start. Get the data from the DAT file.
    call parlst_init(rparlist)

    ! Get the configuration
    sdatafile = './data/renum.dat'
    if (sys_ncommandLineArgs .gt. 0) then
      inquire(file=sys_scommandLineArgs(1,1),exist=bexists)
      if (bexists) sdatafile = sys_scommandLineArgs(1,1)
    end if
    call parlst_readfromfile (rparlist, sdatafile)
    
    ! Start the tests?
    call parlst_getvalue_int (rparlist, 'PERFORMANCETESTS', 'iperform', iperform)
    
    if (iperform .ne. 0) then
    
      call parlst_getvalue_string (rparlist, 'MFLOPTESTS', 'PRMFILE', sprmfile,'./pre/heat_v77.prm')
      call parlst_getvalue_string (rparlist, 'MFLOPTESTS', 'TRIFILE', strifile,'./pre/heat_v77.tri')
    
      ! Read in the parametrisation of the boundary and save
      ! it to rboundary.
      CALL boundary_read_prm(rboundary, sprmfile)
      !CALL boundary_read_prm(rboundary, './pre/QUAD.prm')
          
      ! Loop through all solver configurations
      nconfigs = parlst_querysubstrings(rparList,&
          'PERFORMANCETESTS','configurations')
      
      do iconfig = 1,nconfigs
      
        ! Get the parameters for the configuration
        call parlst_getvalue_string (rparlist, 'PERFORMANCETESTS', &
                                    'configurations', sconfig, '',iconfig)
                                     
        ! Get the configuration
        call readSolverConfig (rparList,trim(sconfig),rsolverConfig)

        ! Loop through the levels
        DO ilevel = rsolverConfig%NLMIN,rsolverConfig%NLMAX
            
          ! Allocate memory for all levels
          ALLOCATE(Rlevels(1:ilevel))
          
          CALL stat_clearTimer (rtimerRef)
          CALL stat_startTimer (rtimerRef)
          
          ! Now read in the basic triangulation into our coarse level.
          CALL tria_readTriFile2D (Rlevels(rsolverConfig%NLCOARSE)%rtriangulation, &
                                  strifile, rboundary)
          !CALL tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
          !                        './pre/QUAD.tri', rboundary)
          
          ! In case we have a triangular element, form a tri mesh
          if (rsolverConfig%ielementType .eq. 1) then
            call tria_rawGridToTri (Rlevels(rsolverConfig%NLCOARSE)%rtriangulation)
          end if
          
          ! Refine it.
          CALL tria_quickRefine2LevelOrdering (rsolverConfig%NLCOARSE-1,&
              Rlevels(rsolverConfig%NLCOARSE)%rtriangulation,rboundary)
          
          ! And create information about adjacencies and everything one needs from
          ! a triangulation.
          CALL tria_initStandardMeshFromRaw (Rlevels(rsolverConfig%NLCOARSE)%rtriangulation,&
              rboundary)
          
          ! Now refine the grid for the fine levels.
          DO i = rsolverConfig%NLCOARSE+1, ilevel

            ! Refine the grid using the 2-Level-Ordering algorithm
            CALL tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
                Rlevels(i)%rtriangulation,rboundary)
            
            ! Create a standard mesh
            CALL tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation,&
              rboundary)
          
          END DO
          
          CALL stat_stopTimer (rtimerRef)

          ! Now we can start to initialise the discretisation. At first, set up
          ! a block discretisation structure that specifies the blocks in the
          ! solution vector. In this simple problem, we only have one block.
          ! Do this for all levels
          DO i = rsolverConfig%NLCOARSE, ilevel
            CALL spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 1, &
                                          Rlevels(i)%rtriangulation, rboundary)
          END DO
          
          ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
          ! structures for every component of the solution vector.
          ! Initialise the first element of the list to specify the element
          ! and cubature rule for this solution component:
          if (rsolverConfig%ielementType .eq. 1) then
            ! P1
            DO i = rsolverConfig%NLCOARSE, ilevel
              CALL spdiscr_initDiscr_simple (&
                  Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
                  EL_E001,CUB_G3_T,Rlevels(i)%rtriangulation, rboundary)
            END DO
          else
            ! Q1
            DO i = rsolverConfig%NLCOARSE, ilevel
              CALL spdiscr_initDiscr_simple (&
                  Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
                  EL_E011,CUB_G2X2,Rlevels(i)%rtriangulation, rboundary)
            END DO
          end if
                        
          CALL stat_clearTimer (rtimerMatGen)
          CALL stat_startTimer (rtimerMatGen)
          
          ! Now as the discretisation is set up, we can start to generate
          ! the structure of the system matrix which is to solve.
          ! We create a scalar matrix, based on the discretisation structure
          ! for our one and only solution component.
          DO i = rsolverConfig%NLCOARSE, ilevel

            ! Initialise the block matrix with default values based on
            ! the discretisation.
            CALL lsysbl_createMatBlockByDiscr (&
                Rlevels(i)%rdiscretisation,Rlevels(i)%rmatrix)    

            ! Now as the discretisation is set up, we can start to generate
            ! the structure of the system matrix which is to solve.
            ! We create that directly in the block (1,1) of the block matrix
            ! using the discretisation structure of the first block.
            CALL bilf_createMatrixStructure ( &
                Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
                LSYSSC_MATRIX9,Rlevels(i)%rmatrix%RmatrixBlock(1,1))
            
            ! Update the structural information of the block matrix, as we manually
            ! changed one of the submatrices:
            CALL lsysbl_updateMatStrucInfo (Rlevels(i)%rmatrix)

            ! And now to the entries of the matrix. For assembling of the entries,
            ! we need a bilinear form, which first has to be set up manually.
            ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
            ! scalar system matrix in 2D.
            rform%itermCount = 2
            rform%Idescriptors(1,1) = DER_DERIV_X
            rform%Idescriptors(2,1) = DER_DERIV_X
            rform%Idescriptors(1,2) = DER_DERIV_Y
            rform%Idescriptors(2,2) = DER_DERIV_Y

            ! In the standard case, we have constant coefficients:
            rform%ballCoeffConstant = .TRUE.
            rform%BconstantCoeff = .TRUE.
            rform%Dcoefficients(1)  = 1.0 
            rform%Dcoefficients(2)  = 1.0 

            ! Now we can build the matrix entries.
            ! We specify the callback function coeff_Laplace for the coefficients.
            ! As long as we use constant coefficients, this routine is not used.
            ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
            ! the framework will call the callback routine to get analytical
            ! data.
            CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                Rlevels(i)%rmatrix%RmatrixBlock(1,1),coeff_Laplace_2D)
            !CALL lsyssc_allocEmptyMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_ZERO)
          
          END DO
          
          CALL stat_stopTimer (rtimerMatGen)
            
          ! Although we could manually create the solution/RHS vector,
          ! the easiest way to set up the vector structure is
          ! to create it by using our matrix as template:
          CALL lsysbl_createVecBlockIndMat (Rlevels(ilevel)%rmatrix,rrhsBlock, .FALSE.)

          ! The vector structure is ready but the entries are missing. 
          ! So the next thing is to calculate the content of that vector.
          !
          ! At first set up the corresponding linear form (f,Phi_j):
          rlinform%itermCount = 1
          rlinform%Idescriptors(1) = DER_FUNC
          
          ! ... and then discretise the RHS to get a discrete version of it.
          ! Again we simply create a scalar vector based on the one and only
          ! discretisation structure.
          ! This scalar vector will later be used as the one and only first
          ! component in a block vector.
          CALL linf_buildVectorScalar (&
              Rlevels(ilevel)%rdiscretisation%RspatialDiscr(1),&
              rlinform,.TRUE.,rrhsBlock%RvectorBlock(1),coeff_RHS_2D)
          
          DO i = rsolverConfig%NLCOARSE, ilevel
          
            ! Initialise the discrete BC structure
            CALL bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)

            ! On edge 1 of boundary component 1 add Dirichlet boundary conditions.      
            CALL boundary_createRegion(rboundary,1,1,rboundaryRegion)
            CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                              rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                              getBoundaryValues_2D)
                                      
            ! Now to the edge 2 of boundary component 1 the domain. 
            CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
            CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                              rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                              getBoundaryValues_2D)
                                      
            ! Edge 3 of boundary component 1.
            CALL boundary_createRegion(rboundary,1,3,rboundaryRegion)
            CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                              rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                              getBoundaryValues_2D)
            
            ! Edge 4 of boundary component 1. That's it.
            CALL boundary_createRegion(rboundary,1,4,rboundaryRegion)
            CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                              rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                              getBoundaryValues_2D)

            ! Edge 1 of boundary component 2. That's it.
            CALL boundary_createRegion(rboundary,2,1,rboundaryRegion)
            CALL bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                              rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                              getBoundaryValues_2D)
            
            ! Hang the pointer into the matrix. That way, these
            ! boundary conditions are always connected to that matrix.
            Rlevels(i)%rmatrix%p_rdiscreteBC => Rlevels(i)%rdiscreteBC
        
            ! Also implement the boundary conditions into the matrix.
            CALL matfil_discreteBC (Rlevels(i)%rmatrix)
            
          END DO

          ! Our right-hand-side also needs to know the boundary conditions.
          rrhsBlock%p_rdiscreteBC => Rlevels(ilevel)%rdiscreteBC
          
          ! Now we have block vectors for the RHS and the matrix. What we
          ! need additionally is a block vector for the solution and
          ! temporary data. Create them using the RHS as template.
          ! Fill the solution vector with 0:
          CALL lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .TRUE.)
          CALL lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .FALSE.)
          CALL lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock2, .FALSE.)
          
          ! Next step is to implement boundary conditions into the RHS,
          ! solution and matrix. This is done using a vector/matrix filter
          ! for discrete boundary conditions.
          ! The discrete boundary conditions are already attached to the
          ! vectors/matrix. Call the appropriate vector/matrix filter that
          ! modifies the vectors/matrix according to the boundary conditions.
          CALL vecfil_discreteBCrhs (rrhsBlock)
          CALL vecfil_discreteBCsol (rvectorBlock)
          
          CALL stat_clearTimer (rtimerSort)
          CALL stat_startTimer (rtimerSort)
          
          ! During the linear solver, the boundary conditions are also
          ! frequently imposed to the vectors. But as the linear solver
          ! does not work with the actual solution vectors but with
          ! defect vectors instead.
          ! So, set up a filter chain that filters the defect vector
          ! during the solution process to implement discrete boundary conditions.
          RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

          ! Resort matrices and vectors
          ALLOCATE(Rdiscretisation(ilevel))
          DO i = rsolverConfig%NLCOARSE, ilevel
          
            ! For stochastic resorting, prepare an array with the discretisation structures.
            CALL spdiscr_duplicateDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
              Rdiscretisation(i), .TRUE.)
        
            ! Prepare a permutation
            CALL storage_new('..', 'Iperm', 2*Rlevels(i)%rmatrix%NEQ, ST_INT, &
                Rlevels(i)%h_Ipermutation,ST_NEWBLOCK_ORDERED)
            CALL storage_getbase_int(Rlevels(i)%h_Ipermutation,p_Ipermutation)
            
            ! Calculate the resorting
            SELECT CASE (rsolverConfig%isortStrategy)
            CASE (0)
              ! Nothing to be done. 2-level ordering
            CASE (1)
              ! Cuthill McKee
              CALL sstrat_calcCuthillMcKee (Rlevels(i)%rmatrix%RmatrixBlock(1,1),p_Ipermutation)

              CALL lsyssc_sortMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),.TRUE.,&
                                      SSTRAT_CM,Rlevels(i)%h_Ipermutation)
              IF (i .EQ. ilevel) THEN
                CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                    rrhsBlock%RvectorBlock(1),rtempBlock2%RvectorBLock(1))
                CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                    rvectorBlock%RvectorBlock(1),rtempBlock2%RvectorBLock(1))
                CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                    rtempBlock%RvectorBlock(1),rtempBlock2%RvectorBlock(1))
              END IF
            CASE (2)
              ! XYZ-Sorting
              CALL sstrat_calcXYZsorting (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
                  p_Ipermutation)

              CALL lsyssc_sortMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),.TRUE.,&
                                      SSTRAT_XYZCOORD,Rlevels(i)%h_Ipermutation)
              IF (i .EQ. ilevel) THEN
                CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                    rrhsBlock%RvectorBlock(1),rtempBlock2%RvectorBLock(1))
                CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                    rvectorBlock%RvectorBlock(1),rtempBlock2%RvectorBLock(1))
                CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                    rtempBlock%RvectorBlock(1),rtempBlock2%RvectorBlock(1))
              END IF

            CASE (3)
              ! Stochastic resorting
              CALL sstrat_calcStochastic (p_Ipermutation)

              CALL lsyssc_sortMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),.TRUE.,&
                                      SSTRAT_STOCHASTIC,Rlevels(i)%h_Ipermutation)
              IF (i .EQ. ilevel) THEN
                CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                    rrhsBlock%RvectorBlock(1),rtempBlock2%RvectorBLock(1))
                CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                    rvectorBlock%RvectorBlock(1),rtempBlock2%RvectorBLock(1))
                CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                    rtempBlock%RvectorBlock(1),rtempBlock2%RvectorBlock(1))
              END IF

            CASE (4)
              ! Hierarchical resorting
              CALL sstrat_calcHierarchical (Rdiscretisation(rsolverConfig%NLCOARSE:i),p_Ipermutation)

              CALL lsyssc_sortMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),.TRUE.,&
                                      SSTRAT_HIERARCHICAL,Rlevels(i)%h_Ipermutation)
              IF (i .EQ. ilevel) THEN
                CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                    rrhsBlock%RvectorBlock(1),rtempBlock2%RvectorBLock(1))
                CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                    rvectorBlock%RvectorBlock(1),rtempBlock2%RvectorBLock(1))
                CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                    rtempBlock%RvectorBlock(1),rtempBlock2%RvectorBlock(1))
              END IF
            END SELECT
            
          END DO
          
          CALL stat_stopTimer (rtimerSort)

          DO i = rsolverConfig%NLCOARSE, ilevel
            CALL spdiscr_releaseDiscr(Rdiscretisation(i))
          END DO
          DEALLOCATE(Rdiscretisation)


          ! Now we have to build up the level information for multigrid.
          !
          ! At first, initialise a standard interlevel projection structure. We
          ! can use the same structure for all levels.
          CALL mlprj_initProjectionMat (rprojection,Rlevels(ilevel)%rmatrix)

          ! Create a Multigrid-solver. Attach the above filter chain
          ! to the solver, so that the solver automatically filters
          ! the vector during the solution process.
          p_RfilterChain => RfilterChain
          
          SELECT CASE (rsolverConfig%isolver)
          
          CASE (0)
          
            CALL linsol_initMultigrid (p_rsolverNode,p_RfilterChain)
            
            ! Set up a coarse grid solver.
            CALL linsol_initUMFPACK4 (p_rcoarsegridSolver)
            
            ! Add the coarse grid level.
            CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverNode,rprojection,&
                                          NULL(), NULL(), p_rcoarseGridSolver)

            ! Now set up the other levels...
            DO i = rsolverConfig%NLCOARSE+1, ilevel
            
              SELECT CASE (rsolverConfig%ismoother)
              CASE (0)
                ! Create a Jacobi smoother
                CALL linsol_initJacobi(p_rsmoother)
              
                ! We will use 4 smoothing steps with damping parameter 0.7
                CALL linsol_convertToSmoother(p_rsmoother, rsolverConfig%nsm, rsolverConfig%domega)
                
              CASE (1)
                ! Create an ILU(s) smoother
                CALL linsol_initMILUs1x1 (p_rsmoother,rsolverConfig%ifillILU,0.0_DP)
                CALL linsol_convertToSmoother(p_rsmoother, rsolverConfig%nsm, rsolverConfig%domega)
              
              END SELECT

              ! And add this multi-grid level. We will use the same smoother
              ! for pre- and post-smoothing.
              CALL linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverNode,rprojection,&
                                            p_rsmoother, p_rsmoother, NULL())
              
            END DO
            
            ! Initialise multigrid-specific parameters
            p_rsolverNode%p_rsubnodeMultigrid%icycle = rsolverConfig%icycle
            
          CASE (1)

            SELECT CASE (rsolverConfig%ipreconditioner)
            CASE (0)
              CALL linsol_initJacobi (p_rsmoother,1.0_DP)
              CALL linsol_initCG (p_rsolverNode,p_rsmoother,p_RfilterChain)
            CASE (1)
              CALL linsol_initMILUs1x1 (p_rsmoother,rsolverConfig%nsm,0.0_DP)
              CALL linsol_initCG (p_rsolverNode,p_rsmoother,p_RfilterChain)
            END SELECT
            
          CASE (2)

            SELECT CASE (rsolverConfig%ipreconditioner)
            CASE (0)
              CALL linsol_initJacobi (p_rsmoother,1.0_DP)
              CALL linsol_initBiCGStab (p_rsolverNode,p_rsmoother,p_RfilterChain)
            CASE (1)
              CALL linsol_initMILUs1x1 (p_rsmoother,rsolverConfig%nsm,0.0_DP)
              CALL linsol_initBiCGStab (p_rsolverNode,p_rsmoother,p_RfilterChain)
            END SELECT
            
          END SELECT
          
          ! Initialise parameters for the solving process
          p_rsolverNode%ioutputlevel = rsolverConfig%ioutputLevel
          p_rsolverNode%nmaxIterations = rsolverConfig%nmaxIterations
          p_rsolverNode%depsRel = rsolverConfig%depsRel
          p_rsolverNode%iresNorm = rsolverConfig%iresNorm
          
          ! Attach the system matrices to the solver.
          !
          ! We copy our matrices to a big matrix array and transfer that
          ! to the setMatrices routines. This intitialises then the matrices
          ! on all levels according to that array. Note that this does not
          ! allocate new memory, we create only 'links' to existing matrices
          ! into Rmatrices(:)!
          ALLOCATE(Rmatrices(rsolverConfig%NLCOARSE:ilevel))
          DO i = rsolverConfig%NLCOARSE, ilevel
            CALL lsysbl_duplicateMatrix (Rlevels(i)%rmatrix,&
                Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          END DO
          
          CALL linsol_setMatrices(p_RsolverNode,Rmatrices(rsolverConfig%NLCOARSE:ilevel))

          ! We can release Rmatrices immediately -- as long as we don't
          ! release Rlevels(i)%rmatrix!
          DO i=rsolverConfig%NLCOARSE,ilevel
            CALL lsysbl_releaseMatrix (Rmatrices(i))
          END DO
          DEALLOCATE(Rmatrices)
          
          ! Initialise structure/data of the solver. This allows the
          ! solver to allocate memory / perform some precalculation
          ! to the problem.
          CALL linsol_initStructure (p_rsolverNode, ierror)
          IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
          CALL linsol_initData (p_rsolverNode, ierror)
          IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
          
          ! Finally solve the system. As we want to solve Ax=b with
          ! b being the real RHS and x being the real solution vector,
          ! we use linsol_solveAdaptively. If b is a defect
          ! RHS and x a defect update to be added to a solution vector,
          ! we would have to use linsol_precondDefect instead.
          
          CALL stat_clearTimer(rtimerSolver)
          CALL stat_startTimer(rtimerSolver)
          CALL linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rrhsBlock,rtempBlock)
          CALL stat_stopTimer(rtimerSolver)
          
          ! That's it, rvectorBlock now contains our solution. We can now
          ! start the postprocessing. 
          ! Start UCD export to GMV file:
          if (rsolverConfig%iwriteGMV .ne. 0) then
            CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,&
                Rlevels(ilevel)%rtriangulation,'gmv/u_'//&
                trim(sys_lowcase(sconfig))//'_lv_'//trim(sys_siL(ilevel,10))//'.gmv')
            
            CALL lsyssc_vectorActivateSorting (rvectorBlock%RvectorBlock(1),.FALSE.)
            CALL lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
            CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
            
            ! Write the file to disc, that's it.
            CALL ucd_write (rexport)
            CALL ucd_release (rexport)
          end if
          
          ! Calculate the error to the reference function.
    !        CALL pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,derror,&
    !                          getReferenceFunction_2D)
    !        CALL output_line ('L2-error: ' // sys_sdEL(derror,10) )
    !
    !        CALL pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
    !                          getReferenceFunction_2D)
    !        CALL output_line ('H1-error: ' // sys_sdEL(derror,10) )
    !        

          WRITE (*,'(A,A,A,I3,A,E18.10,E18.10,A,I5,A,E18.10,A,E18.10,A,E18.10,A,E18.10)') &
            'Solver: ',trim(sys_lowcase(sconfig)),&
            ' Level ',ilevel,&
            ': ',rtimerSolver%delapsedReal,rtimerSolver%delapsedCPU,&
            ' / ',p_rsolverNode%iiterations,&
            ' KRate: ',p_rsolverNode%dconvergenceRate,&
            ' Ref: ',rtimerRef%delapsedReal,&
            ' Sort: ',rtimerSort%delapsedReal,&
            ' MatGen: ',rtimerMatGen%delapsedReal

          ! We are finished - but not completely!
          ! Now, clean up so that all the memory is available again.
          !
          ! Release solver data and structure
          CALL linsol_doneData (p_rsolverNode)
          CALL linsol_doneStructure (p_rsolverNode)
          
          ! Release the interlevel projection structure
          CALL mlprj_doneProjection (rprojection)
          
          ! Release the solver node and all subnodes attached to it (if at all):
          CALL linsol_releaseSolver (p_rsolverNode)
          
          ! Release the block matrix/vectors
          CALL lsysbl_releaseVector (rtempBlock)
          CALL lsysbl_releaseVector (rtempBlock2)
          CALL lsysbl_releaseVector (rvectorBlock)
          CALL lsysbl_releaseVector (rrhsBlock)
          DO i = ilevel, rsolverConfig%NLCOARSE, -1
            CALL lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
            CALL storage_free (Rlevels(i)%h_Ipermutation)
          END DO

          ! Release our discrete version of the boundary conditions
          DO i = ilevel, rsolverConfig%NLCOARSE, -1
            CALL bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
          END DO

          ! Release the discretisation structure and all spatial discretisation
          ! structures in it.
          DO i = ilevel, rsolverConfig%NLCOARSE, -1
            CALL spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
          END DO
          
          ! Release the triangulation. 
          DO i = ilevel, rsolverConfig%NLCOARSE, -1
            IF (ilevel .EQ. rsolverconfig%NLMAX) THEN
              CALL tria_infoStatistics (&
                  Rlevels(i)%rtriangulation,i .EQ. rsolverConfig%NLMAX,i)
            END IF
            CALL tria_done (Rlevels(i)%rtriangulation)
          END DO
          
          DEALLOCATE(Rlevels)
                
        END DO
      
      end do ! iconfig
      
      ! Finally release the domain, that's it.
      CALL boundary_release (rboundary)

    end if

    call parlst_done(rparlist)

  END SUBROUTINE

END MODULE

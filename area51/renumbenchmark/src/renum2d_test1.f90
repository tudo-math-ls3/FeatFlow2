!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# MFLOP test suite for testing the speed of the matrix vector multiplication
!# on a number of levels.
!# </purpose>
!##############################################################################

MODULE renum2d_test1

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use filtersupport
  use linearsystemscalar
  use linearsystemblock
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use scalarpde
  use ucd
  use pprocerror
  use genoutput
  use paramlist
  use statistics
  use sortstrategy
    
  USE renum2d_callback
  
  IMPLICIT NONE

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

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE renum2d_matvec_mfloptest
  
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
    TYPE(t_level), DIMENSION(:), pointer :: Rlevels

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

    ! NLMIN receives the level of the coarse grid.
    INTEGER :: NLMIN,ilevmin

    ! NLMAX receives the level where we want to solve.
    INTEGER :: NLMAX,ilevmax
    
    ! Some temporary variables
    INTEGER :: i
    REAL(DP) :: dmflops
    
    ! Sorting strategy counter
    INTEGER :: isortStrategy, imesh, isortStrategyidx, imeshidx
    INTEGER, DIMENSION(:), POINTER :: p_Ipermutation
    TYPE(t_spatialDiscretisation), DIMENSION(:), POINTER :: Rdiscretisation
    TYPE(t_matrixScalar), POINTER :: p_rmatrix
    
    TYPE(t_timer) :: rtimerSolver
    
    character(len=SYS_STRLEN) :: sprmfile,strifile,sdatafile
    logical :: bexists
    TYPE(t_triangulation) :: rtriangulation

    type(t_parlist) :: rparlist
    integer :: iperform

    ! Ok, let's start. Get the data from the DAT file.
    call parlst_init(rparlist)

    ! Get the configuration
    call sys_getcommandLineArg(1,sdatafile,sdefault='./data/renum.dat')
    call parlst_readfromfile (rparlist, sdatafile)
    
    ! Start the tests?
    call parlst_getvalue_int (rparlist, 'MFLOPTESTS', 'iperform', iperform)
    
    if (iperform .ne. 0) then

      ! Get the test parameters
      call parlst_getvalue_int (rparlist, 'MFLOPTESTS', 'NLMIN', ilevmin)
      call parlst_getvalue_int (rparlist, 'MFLOPTESTS', 'NLMAX', ilevmax)

      ! At first, read in the parametrisation of the boundary and save
      ! it to rboundary.
      call parlst_getvalue_string (rparlist, 'MFLOPTESTS', 'PRMFILE', sprmfile,'./pre/heat_v77.prm')
      call parlst_getvalue_string (rparlist, 'MFLOPTESTS', 'TRIFILE', strifile,'./pre/heat_v77.tri')
      
      CALL boundary_read_prm(rboundary, sprmfile)
          
      ! Ok, let's start. 
      
      NLMIN = 1
      
      DO imeshidx=1,parlst_querysubstrings(rparlist,'MFLOPTESTS','imesh')
      
        ! Get the mesh to test
        call parlst_getvalue_int (rparlist, 'MFLOPTESTS', 'imesh', &
            imesh,iarrayindex=imeshidx)

        DO NLMAX=ilevmin,ilevmax
        
          DO isortStrategyidx=1,parlst_querysubstrings(rparlist,'MFLOPTESTS','isortStrategy')
          
            ! Get the sort strategy to test
            call parlst_getvalue_int (rparlist, 'MFLOPTESTS', 'isortstrategy', &
                isortstrategy,iarrayindex=isortstrategyidx)
        
            ! Allocate memory for all levels
            ALLOCATE(Rlevels(NLMIN:NLMAX))
            
            IF (imesh .EQ. 0) THEN
              ! Now read in the basic triangulation into our coarse level.
              CALL tria_readTriFile2D (rtriangulation, strifile, rboundary)
                                      
              ! Only macro 1
              CALL tria_initStandardMeshFromRaw (rtriangulation,rboundary)
              
              CALL tria_generateSubdomain(rtriangulation,(/1/),&
                  Rlevels(NLMIN)%rtriangulation,rboundary)
             
              CALL tria_done (rtriangulation)

            ELSE
     
              CALL tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
                                      strifile, rboundary)
            END IF
                                    
            ! Refine it.
            CALL tria_quickRefine2LevelOrdering (NLMIN-1,&
                Rlevels(NLMIN)%rtriangulation,rboundary)
            
            ! And create information about adjacencies and everything one needs from
            ! a triangulation.
            CALL tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtriangulation,&
                rboundary)
            
            ! Now refine the grid for the fine levels.
            DO i = NLMIN+1, NLMAX

              ! Refine the grid using the 2-Level-Ordering algorithm
              CALL tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
                  Rlevels(i)%rtriangulation,rboundary)
              
              ! Create a standard mesh
              CALL tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation,&
                rboundary)
            
            END DO

            ! Now we can start to initialise the discretisation. At first, set up
            ! a block discretisation structure that specifies the blocks in the
            ! solution vector. In this simple problem, we only have one block.
            ! Do this for all levels
            DO i = NLMIN, NLMAX
              CALL spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 1, &
                                            Rlevels(i)%rtriangulation, rboundary)
            END DO
            
            ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
            ! structures for every component of the solution vector.
            ! Initialise the first element of the list to specify the element
            ! and cubature rule for this solution component:
            DO i = NLMIN, NLMAX
              CALL spdiscr_initDiscr_simple (&
                  Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
                  EL_E011,CUB_G2X2,Rlevels(i)%rtriangulation, rboundary)
            END DO
                          
            ! Now as the discretisation is set up, we can start to generate
            ! the structure of the system matrix which is to solve.
            ! We create a scalar matrix, based on the discretisation structure
            ! for our one and only solution component.
            DO i = NLMIN, NLMAX

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
            
            END DO
              
            ! Although we could manually create the solution/RHS vector,
            ! the easiest way to set up the vector structure is
            ! to create it by using our matrix as template:
            CALL lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrix,rrhsBlock, .FALSE.)

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
                Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(1),&
                rlinform,.TRUE.,rrhsBlock%RvectorBlock(1),coeff_RHS_2D)
            
            DO i = NLMIN, NLMAX
            
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
            rrhsBlock%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBC
            
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
            
            ! Resort matrices and vectors
            ALLOCATE(Rdiscretisation(NLMAX))
            DO i = NLMIN, NLMAX
            
              ! For stochastic resorting, prepare an array with the discretisation structures.
              CALL spdiscr_duplicateDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
                Rdiscretisation(i), .TRUE.)
          
              ! Prepare a permutation
              CALL storage_new('..', 'Iperm', 2*Rlevels(i)%rmatrix%NEQ, ST_INT, &
                  Rlevels(i)%h_Ipermutation,ST_NEWBLOCK_ORDERED)
              CALL storage_getbase_int(Rlevels(i)%h_Ipermutation,p_Ipermutation)
              
              ! Calculate the resorting
              SELECT CASE (isortStrategy)
              CASE (0)
                ! Nothing to be done. 2-level ordering
              CASE (1)
                ! Cuthill McKee
                CALL sstrat_calcCuthillMcKee (Rlevels(i)%rmatrix%RmatrixBlock(1,1),p_Ipermutation)

                CALL lsyssc_sortMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),.TRUE.,&
                                        SSTRAT_CM,Rlevels(i)%h_Ipermutation)
                IF (i .EQ. NLMAX) THEN
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
                IF (i .EQ. NLMAX) THEN
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
                IF (i .EQ. NLMAX) THEN
                  CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                      rrhsBlock%RvectorBlock(1),rtempBlock2%RvectorBLock(1))
                  CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                      rvectorBlock%RvectorBlock(1),rtempBlock2%RvectorBLock(1))
                  CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                      rtempBlock%RvectorBlock(1),rtempBlock2%RvectorBlock(1))
                END IF

              CASE (4)
                ! Hierarchical resorting
                CALL sstrat_calcHierarchical (Rdiscretisation(NLMIN:i),p_Ipermutation)

                CALL lsyssc_sortMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),.TRUE.,&
                                        SSTRAT_HIERARCHICAL,Rlevels(i)%h_Ipermutation)
                IF (i .EQ. NLMAX) THEN
                  CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                      rrhsBlock%RvectorBlock(1),rtempBlock2%RvectorBLock(1))
                  CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                      rvectorBlock%RvectorBlock(1),rtempBlock2%RvectorBLock(1))
                  CALL lsyssc_synchroniseSortMatVec (Rlevels(i)%rmatrix%RmatrixBlock(1,1),&
                      rtempBlock%RvectorBlock(1),rtempBlock2%RvectorBlock(1))
                END IF
              END SELECT
              
            END DO

            DO i = NLMIN, NLMAX
              CALL spdiscr_releaseDiscr(Rdiscretisation(i))
            END DO
            DEALLOCATE(Rdiscretisation)


            
            ! Finally solve the system. As we want to solve Ax=b with
            ! b being the real RHS and x being the real solution vector,
            ! we use linsol_solveAdaptively. If b is a defect
            ! RHS and x a defect update to be added to a solution vector,
            ! we would have to use linsol_precondDefect instead.
            
            p_rmatrix => Rlevels(NLMAX)%rmatrix%RmatrixBlock(1,1)
            CALL stat_clearTimer(rtimerSolver)
            CALL stat_startTimer(rtimerSolver)
            
            ! Dop a couple of MV-Multiplications
            DO i=1,100
              CALL lsyssc_scalarMatVec (p_rmatrix, &
                  rvectorBlock%RvectorBlock(1), rrhsBlock%RvectorBlock(1), -1.0_DP, 1.0_DP)
            END DO
            
            CALL stat_stopTimer(rtimerSolver)
            
            ! Flops:
            dmflops = 2.0_DP * REAL(i-1,DP) * REAL(p_rmatrix%NA,DP)
            
            ! MFLOPs:
            dmflops = dmflops / (rtimerSolver%delapsedReal * 1000000.0_DP)
            
            ! Release the block matrix/vectors
            CALL lsysbl_releaseVector (rtempBlock)
            CALL lsysbl_releaseVector (rtempBlock2)
            CALL lsysbl_releaseVector (rvectorBlock)
            CALL lsysbl_releaseVector (rrhsBlock)
            DO i = NLMAX, NLMIN, -1
              CALL lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
              CALL storage_free (Rlevels(i)%h_Ipermutation)
            END DO

            ! Release our discrete version of the boundary conditions
            DO i = NLMAX, NLMIN, -1
              CALL bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
            END DO

            ! Release the discretisation structure and all spatial discretisation
            ! structures in it.
            DO i = NLMAX, NLMIN, -1
              CALL spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
            END DO
            
            WRITE (*,'(A,I3,A,I3,A,I3,A,E18.10,A,E18.10,A,F18.10)') &
              'Level ',NLMAX,&
              ', Strategy ',isortStrategy,&
              ', imesh ',imesh,&
              ', Time=',rtimerSolver%delapsedReal,&
              ', MFLOPS= ',dmflops
            
            ! Release the triangulation. 
            DO i = NLMAX, NLMIN, -1
              IF ((NLMAX .EQ. 11) .AND. (isortStrategy .EQ. 4)) THEN
                CALL tria_infoStatistics (Rlevels(i)%rtriangulation,i .EQ. NLMAX,i)
              END IF
              CALL tria_done (Rlevels(i)%rtriangulation)
            END DO
            
            DEALLOCATE(Rlevels)
            
          END DO

        END DO
        
      END DO

      ! Finally release the domain, that's it.
      CALL boundary_release (rboundary)

    end if

    call parlst_done(rparlist)

  END SUBROUTINE

END MODULE

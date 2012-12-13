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

module renum2d_test1

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
  use sortstrategybase
  use sortstrategy
    
  use renum2d_callback
  
  implicit none

  type t_level
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (structure of the
    ! solution, trial/test functions,...)
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the
    ! discrete Laplace operator.
    type(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC) :: rdiscreteBC
    
    ! Resorting strategy
    type(t_blockSortStrategy) :: rsortStrategy
  
  end type

contains

  ! ***************************************************************************

!<subroutine>

  subroutine renum2d_matvec_mfloptest
  
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
    type(t_level), dimension(:), pointer :: Rlevels

    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock,rtempBlock2

    ! A variable that is used to specify a region on the boundary.
    type(t_boundaryRegion) :: rboundaryRegion

    ! NLMIN receives the level of the coarse grid.
    integer :: NLMIN,ilevmin

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX,ilevmax
    
    ! Some temporary variables
    integer :: i
    integer(I32) :: celementtype
    real(DP) :: dmflops
    
    ! Sorting strategy counter
    integer :: isortStrategy, imesh, isortStrategyidx, imeshidx
    integer, dimension(:), pointer :: p_Ipermutation
    type(t_spatialDiscretisation), dimension(:), pointer :: Rdiscretisation
    type(t_matrixScalar), pointer :: p_rmatrix
    
    type(t_timer) :: rtimerSolver
    
    character(len=SYS_STRLEN) :: sprmfile,strifile,sdatafile,selement
    logical :: bexists
    type(t_triangulation) :: rtriangulation

    type(t_parlist) :: rparlist
    integer :: iperform

    ! Ok, let's start. Get the data from the DAT file.
    call parlst_init(rparlist)

    ! Get the configuration
    call sys_getcommandLineArg(1,sdatafile,sdefault='./data/renum.dat')
    call parlst_readfromfile (rparlist, sdatafile)
    
    ! Start the tests?
    call parlst_getvalue_int (rparlist, 'MFLOPTESTS', 'iperform', iperform)
    
    call parlst_getvalue_string (rparlist, 'MFLOPTESTS', 'selement', &
        selement,'EL_Q1')

    ! Get the correct element id
    celementtype = elem_igetID(selement)
    
    if (iperform .ne. 0) then

      ! Get the test parameters
      call parlst_getvalue_int (rparlist, 'MFLOPTESTS', 'NLMIN', ilevmin)
      call parlst_getvalue_int (rparlist, 'MFLOPTESTS', 'NLMAX', ilevmax)

      ! At first, read in the parametrisation of the boundary and save
      ! it to rboundary.
      call parlst_getvalue_string (rparlist, 'MFLOPTESTS', 'PRMFILE', sprmfile,'./pre/heat_v77.prm')
      call parlst_getvalue_string (rparlist, 'MFLOPTESTS', 'TRIFILE', strifile,'./pre/heat_v77.tri')
      
      call boundary_read_prm(rboundary, sprmfile)
          
      ! Ok, let's start.
      
      NLMIN = 1
      
      do imeshidx=1,parlst_querysubstrings(rparlist,'MFLOPTESTS','imesh')
      
        ! Get the mesh to test
        call parlst_getvalue_int (rparlist, 'MFLOPTESTS', 'imesh', &
            imesh,iarrayindex=imeshidx)

        do NLMAX=ilevmin,ilevmax
        
          do isortStrategyidx=1,parlst_querysubstrings(rparlist,'MFLOPTESTS','isortStrategy')
          
            ! Get the sort strategy to test
            call parlst_getvalue_int (rparlist, 'MFLOPTESTS', 'isortstrategy', &
                isortstrategy,iarrayindex=isortstrategyidx)
        
            ! Allocate memory for all levels
            allocate(Rlevels(NLMIN:NLMAX))
            
            if (imesh .eq. 0) then
              ! Now read in the basic triangulation into our coarse level.
              call tria_readTriFile2D (rtriangulation, strifile, rboundary)
                                      
              ! In case we have a triangular element, form a tri mesh
              if (elem_igetNVE(celementtype) .eq. 3) then
                call tria_rawGridToTri (rtriangulation)
              end if

              ! Only macro 1
              call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
              
              call tria_generateSubdomain(rtriangulation,(/1/),&
                  Rlevels(NLMIN)%rtriangulation,rboundary)
             
              call tria_done (rtriangulation)

            else
     
              call tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
                                      strifile, rboundary)
            end if
                                    
            ! Refine it.
            call tria_quickRefine2LevelOrdering (NLMIN-1,&
                Rlevels(NLMIN)%rtriangulation,rboundary)
            
            ! And create information about adjacencies and everything one needs from
            ! a triangulation.
            call tria_initStandardMeshFromRaw (Rlevels(NLMIN)%rtriangulation,&
                rboundary)
            
            ! Now refine the grid for the fine levels.
            do i = NLMIN+1, NLMAX

              ! Refine the grid using the 2-Level-Ordering algorithm
              call tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
                  Rlevels(i)%rtriangulation,rboundary)
              
              ! Create a standard mesh
              call tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation,&
                rboundary)
            
            end do

            ! Now we can start to initialise the discretisation. At first, set up
            ! a block discretisation structure that specifies the blocks in the
            ! solution vector. In this simple problem, we only have one block.
            ! Do this for all levels
            do i = NLMIN, NLMAX
              call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 1, &
                                            Rlevels(i)%rtriangulation, rboundary)
            end do
            
            ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
            ! structures for every component of the solution vector.
            ! Initialise the first element of the list to specify the element
            ! and cubature rule for this solution component:
            do i = NLMIN, NLMAX
              call spdiscr_initDiscr_simple (&
                  Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
                  celementtype,spdiscr_getStdCubature(celementType),&
                  Rlevels(i)%rtriangulation, rboundary)
            end do
                          
            ! Now as the discretisation is set up, we can start to generate
            ! the structure of the system matrix which is to solve.
            ! We create a scalar matrix, based on the discretisation structure
            ! for our one and only solution component.
            do i = NLMIN, NLMAX

              ! Initialise the block matrix with default values based on
              ! the discretisation.
              call lsysbl_createMatBlockByDiscr (&
                  Rlevels(i)%rdiscretisation,Rlevels(i)%rmatrix)

              ! Now as the discretisation is set up, we can start to generate
              ! the structure of the system matrix which is to solve.
              ! We create that directly in the block (1,1) of the block matrix
              ! using the discretisation structure of the first block.
              call bilf_createMatrixStructure ( &
                  Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
                  LSYSSC_MATRIX9,Rlevels(i)%rmatrix%RmatrixBlock(1,1))
              
              ! Update the structural information of the block matrix, as we manually
              ! changed one of the submatrices:
              call lsysbl_updateMatStrucInfo (Rlevels(i)%rmatrix)

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
              rform%ballCoeffConstant = .true.
              rform%BconstantCoeff = .true.
              rform%Dcoefficients(1)  = 1.0
              rform%Dcoefficients(2)  = 1.0

              ! Now we can build the matrix entries.
              ! We specify the callback function coeff_Laplace for the coefficients.
              ! As long as we use constant coefficients, this routine is not used.
              ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
              ! the framework will call the callback routine to get analytical
              ! data.
              call bilf_buildMatrixScalar (rform,.true.,&
                  Rlevels(i)%rmatrix%RmatrixBlock(1,1),coeff_Laplace_2D)
            
            end do
              
            ! Although we could manually create the solution/RHS vector,
            ! the easiest way to set up the vector structure is
            ! to create it by using our matrix as template:
            call lsysbl_createVecBlockIndMat (Rlevels(NLMAX)%rmatrix,rrhsBlock, .false.)

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
            call linf_buildVectorScalar (&
                Rlevels(NLMAX)%rdiscretisation%RspatialDiscr(1),&
                rlinform,.true.,rrhsBlock%RvectorBlock(1),coeff_RHS_2D)
            
            do i = NLMIN, NLMAX
            
              ! Initialise the discrete BC structure
              call bcasm_initDiscreteBC(Rlevels(i)%rdiscreteBC)

              ! On edge 1 of boundary component 1 add Dirichlet boundary conditions.
              call boundary_createRegion(rboundary,1,1,rboundaryRegion)
              call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                                rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                                getBoundaryValues_2D)
                                        
              ! Now to the edge 2 of boundary component 1 the domain.
              call boundary_createRegion(rboundary,1,2,rboundaryRegion)
              call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                                rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                                getBoundaryValues_2D)
                                        
              ! Edge 3 of boundary component 1.
              call boundary_createRegion(rboundary,1,3,rboundaryRegion)
              call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                                rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                                getBoundaryValues_2D)
              
              ! Edge 4 of boundary component 1. That's it.
              call boundary_createRegion(rboundary,1,4,rboundaryRegion)
              call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                                rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                                getBoundaryValues_2D)

              ! Edge 1 of boundary component 2. That's it.
              call boundary_createRegion(rboundary,2,1,rboundaryRegion)
              call bcasm_newDirichletBConRealBD (Rlevels(i)%rdiscretisation,1,&
                                                rboundaryRegion,Rlevels(i)%rdiscreteBC,&
                                                getBoundaryValues_2D)
              
              ! Hang the pointer into the matrix. That way, these
              ! boundary conditions are always connected to that matrix.
              Rlevels(i)%rmatrix%p_rdiscreteBC => Rlevels(i)%rdiscreteBC
          
              ! Also implement the boundary conditions into the matrix.
              call matfil_discreteBC (Rlevels(i)%rmatrix)
              
            end do

            ! Our right-hand-side also needs to know the boundary conditions.
            rrhsBlock%p_rdiscreteBC => Rlevels(NLMAX)%rdiscreteBC
            
            ! Now we have block vectors for the RHS and the matrix. What we
            ! need additionally is a block vector for the solution and
            ! temporary data. Create them using the RHS as template.
            ! Fill the solution vector with 0:
            call lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .true.)
            call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .false.)
            call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock2, .false.)
            
            ! Next step is to implement boundary conditions into the RHS,
            ! solution and matrix. This is done using a vector/matrix filter
            ! for discrete boundary conditions.
            ! The discrete boundary conditions are already attached to the
            ! vectors/matrix. Call the appropriate vector/matrix filter that
            ! modifies the vectors/matrix according to the boundary conditions.
            call vecfil_discreteBCrhs (rrhsBlock)
            call vecfil_discreteBCsol (rvectorBlock)
            
            ! Resort matrices and vectors
            allocate(Rdiscretisation(NLMAX))
            do i = NLMIN, NLMAX
            
              ! For stochastic resorting, prepare an array with the discretisation structures.
              call spdiscr_duplicateDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
                Rdiscretisation(i), .true.)
          
              ! Create a sort strategy structure for our discretisation
              call sstrat_initBlockSorting (Rlevels(i)%rsortStrategy,Rlevels(i)%rdiscretisation)
              
              ! Calculate the resorting
              select case (isortStrategy)
              case (0)
                ! Nothing to be done. 2-level ordering

              case (1)
                ! Cuthill McKee
                call sstrat_initCuthillMcKee (Rlevels(i)%rsortStrategy%p_Rstrategies(1),&
                    Rlevels(i)%rmatrix%RmatrixBlock(1,1))

                ! Attach the sorting strategy to the matrix. The matrix is not yet sorted.
                call lsysbl_setSortStrategy (Rlevels(i)%rmatrix,&
                    Rlevels(i)%rsortStrategy,Rlevels(i)%rsortStrategy)
                
                ! Sort the matrix
                call lsysbl_sortMatrix (Rlevels(i)%rmatrix,.true.)
                
                ! Sort the vectors on the maximum level.
                if (i .eq. NLMAX) then
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rrhsBlock,rtempBlock2%RvectorBlock(1))
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rvectorBlock,rtempBlock2%RvectorBlock(1))
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rtempBlock,rtempBlock2%RvectorBlock(1))
                end if

              case (2)
                ! XYZ-Sorting
                call sstrat_initXYZsorting (Rlevels(i)%rsortStrategy%p_Rstrategies(1),&
                    Rlevels(i)%rdiscretisation%RspatialDiscr(1))

                ! Attach the sorting strategy to the matrix. The matrix is not yet sorted.
                call lsysbl_setSortStrategy (Rlevels(i)%rmatrix,&
                    Rlevels(i)%rsortStrategy,Rlevels(i)%rsortStrategy)
                
                ! Sort the matrix
                call lsysbl_sortMatrix (Rlevels(i)%rmatrix,.true.)
                
                ! Sort the vectors on the maximum level.
                if (i .eq. NLMAX) then
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rrhsBlock,rtempBlock2%RvectorBlock(1))
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rvectorBlock,rtempBlock2%RvectorBlock(1))
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rtempBlock,rtempBlock2%RvectorBlock(1))
                end if

              case (3)
                ! Stochastic resorting
                call sstrat_initStochastic (Rlevels(i)%rsortStrategy%p_Rstrategies(1))

                ! Attach the sorting strategy to the matrix. The matrix is not yet sorted.
                call lsysbl_setSortStrategy (Rlevels(i)%rmatrix,&
                    Rlevels(i)%rsortStrategy,Rlevels(i)%rsortStrategy)
                
                ! Sort the matrix
                call lsysbl_sortMatrix (Rlevels(i)%rmatrix,.true.)
                
                ! Sort the vectors on the maximum level.
                if (i .eq. NLMAX) then
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rrhsBlock,rtempBlock2%RvectorBlock(1))
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rvectorBlock,rtempBlock2%RvectorBlock(1))
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rtempBlock,rtempBlock2%RvectorBlock(1))
                end if

              case (4)
                ! Hierarchical resorting
                call sstrat_initHierarchical (Rlevels(i)%rsortStrategy%p_Rstrategies(1),&
                    Rdiscretisation(NLMIN:i))

                ! Attach the sorting strategy to the matrix. The matrix is not yet sorted.
                call lsysbl_setSortStrategy (Rlevels(i)%rmatrix,&
                    Rlevels(i)%rsortStrategy,Rlevels(i)%rsortStrategy)
                
                ! Sort the matrix
                call lsysbl_sortMatrix (Rlevels(i)%rmatrix,.true.)
                
                ! Sort the vectors on the maximum level.
                if (i .eq. NLMAX) then
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rrhsBlock,rtempBlock2%RvectorBlock(1))
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rvectorBlock,rtempBlock2%RvectorBlock(1))
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rtempBlock,rtempBlock2%RvectorBlock(1))
                end if

              case (5)
                ! Reverse Cuthill McKee
                call sstrat_initRevCuthillMcKee (Rlevels(i)%rsortStrategy%p_Rstrategies(1),&
                    Rlevels(i)%rmatrix%RmatrixBlock(1,1))

                ! Attach the sorting strategy to the matrix. The matrix is not yet sorted.
                call lsysbl_setSortStrategy (Rlevels(i)%rmatrix,&
                    Rlevels(i)%rsortStrategy,Rlevels(i)%rsortStrategy)
                
                ! Sort the matrix
                call lsysbl_sortMatrix (Rlevels(i)%rmatrix,.true.)
                
                ! Sort the vectors on the maximum level.
                if (i .eq. NLMAX) then
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rrhsBlock,rtempBlock2%RvectorBlock(1))
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rvectorBlock,rtempBlock2%RvectorBlock(1))
                  call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rtempBlock,rtempBlock2%RvectorBlock(1))
                end if
              end select
              
            end do

            do i = NLMIN, NLMAX
              call spdiscr_releaseDiscr(Rdiscretisation(i))
            end do
            deallocate(Rdiscretisation)
            
            ! Finally solve the system. As we want to solve Ax=b with
            ! b being the real RHS and x being the real solution vector,
            ! we use linsol_solveAdaptively. If b is a defect
            ! RHS and x a defect update to be added to a solution vector,
            ! we would have to use linsol_precondDefect instead.
            
            p_rmatrix => Rlevels(NLMAX)%rmatrix%RmatrixBlock(1,1)
            call stat_clearTimer(rtimerSolver)
            call stat_startTimer(rtimerSolver)
            
            ! Dop a couple of MV-Multiplications
            do i=1,100
              call lsyssc_scalarMatVec (p_rmatrix, &
                  rvectorBlock%RvectorBlock(1), rrhsBlock%RvectorBlock(1), -1.0_DP, 1.0_DP)
            end do
            
            call stat_stopTimer(rtimerSolver)
            
            ! Flops:
            dmflops = 2.0_DP * real(i-1,DP) * real(p_rmatrix%NA,DP)
            
            ! MFLOPs:
            if (rtimerSolver%delapsedReal .eq. 0.0_DP) then
              dmflops = 0.0
            else
              dmflops = dmflops / (rtimerSolver%delapsedReal * 1000000.0_DP)
            end if
            
            ! Release the block matrix/vectors
            call lsysbl_releaseVector (rtempBlock)
            call lsysbl_releaseVector (rtempBlock2)
            call lsysbl_releaseVector (rvectorBlock)
            call lsysbl_releaseVector (rrhsBlock)
            do i = NLMAX, NLMIN, -1
              call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
              call sstrat_doneBlockSorting (Rlevels(i)%rsortStrategy)
            end do

            ! Release our discrete version of the boundary conditions
            do i = NLMAX, NLMIN, -1
              call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
            end do

            ! Release the discretisation structure and all spatial discretisation
            ! structures in it.
            do i = NLMAX, NLMIN, -1
              call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
            end do
            
            write (*,'(A,I3,A,I3,A,I3,A,E18.10,A,E18.10,A,F18.10)') &
              'Level ',NLMAX,&
              ', Strategy ',isortStrategy,&
              ', imesh ',imesh,&
              ', Time=',rtimerSolver%delapsedReal,&
              ', MFLOPS= ',dmflops
            
            ! Release the triangulation.
            do i = NLMAX, NLMIN, -1
              if ((NLMAX .eq. 11) .and. (isortStrategy .eq. 4)) then
                call tria_infoStatistics (Rlevels(i)%rtriangulation,i .eq. NLMAX,i)
              end if
              call tria_done (Rlevels(i)%rtriangulation)
            end do
            
            deallocate(Rlevels)
            
          end do

        end do
        
      end do

      ! Finally release the domain, that's it.
      call boundary_release (rboundary)

    end if

    call parlst_done(rparlist)

  end subroutine

end module

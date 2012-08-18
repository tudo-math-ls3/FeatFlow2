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

module renum2d_test2

  use fsystem
  use genoutput
  use storage
  use boundary
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
  use multilevelprojection
  use linearsolver
  use bilinearformevaluation
  use linearformevaluation
  use paramlist
  use statistics
  use sortstrategybase
  use sortstrategy
    
  use renum2d_callback
  
  implicit none

!<types>

!<typeblock description="Type block defining all information about one level">

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
    integer(I32) :: ielementType
  
  end type
  
!</typeblock>

!</types>

contains

  subroutine readSolverConfig (rparList,ssection,rsolverConfig)
  
!<description>
  ! Reads a solver configuration from a parameter list.
!</description>

  type(t_parlist), intent(in) :: rparList
  character(len=*), intent(in) :: ssection
  type(t_solverConfig), intent(out) :: rsolverConfig
  
  character(len=SYS_STRLEN) :: selement
  
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

    call parlst_getvalue_string (rparlist, 'PERFORMANCETESTS', 'selement', &
        selement,'EL_Q1')
    call parlst_getvalue_string (rparlist,ssection,'selement',&
        selement,selement)
    
    ! Get the correct element id
    rsolverConfig%ielementType = elem_igetID(selement)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine renum2d_solvertest
  
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
    type(t_level), dimension(:), target, allocatable :: Rlevels

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

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rcoarseGridSolver,p_rsmoother

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), allocatable :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection

    ! One level of multigrid
    type(t_linsolMGLevelInfo), pointer :: p_rlevelInfo
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Some temporary variables
    integer :: i
    
    ! Sorting strategy counter
    type(t_spatialDiscretisation), dimension(:), allocatable :: Rdiscretisation
    
    type(t_timer) :: rtimerSolver,rtimerRef,rtimerSort,rtimerMatGen
    
    ! DAT-file parameters
    type(t_parlist) :: rparlist
    integer :: iconfig,nconfigs
    character(len=SYS_STRLEN) :: sconfig
    type(t_solverConfig) :: rsolverConfig
    integer :: ilevel, iperform
    character(len=SYS_STRLEN) :: sprmfile,strifile,sdatafile

    ! Ok, let's start. Get the data from the DAT file.
    call parlst_init(rparlist)

    ! Get the configuration
    call sys_getcommandLineArg(1,sdatafile,sdefault='./data/renum.dat')
    call parlst_readfromfile (rparlist, sdatafile)
    
    ! Start the tests?
    call parlst_getvalue_int (rparlist, 'PERFORMANCETESTS', 'iperform', iperform)
    
    if (iperform .ne. 0) then
    
      call parlst_getvalue_string (rparlist, 'PERFORMANCETESTS', 'PRMFILE', sprmfile,'./pre/heat_v77.prm')
      call parlst_getvalue_string (rparlist, 'PERFORMANCETESTS', 'TRIFILE', strifile,'./pre/heat_v77.tri')
    
      ! Read in the parametrisation of the boundary and save
      ! it to rboundary.
      call boundary_read_prm(rboundary, sprmfile)
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
        do ilevel = rsolverConfig%NLMIN,rsolverConfig%NLMAX
            
          ! Allocate memory for all levels
          allocate(Rlevels(1:ilevel))
          
          call stat_clearTimer (rtimerRef)
          call stat_startTimer (rtimerRef)
          
          ! Now read in the basic triangulation into our coarse level.
          call tria_readTriFile2D (Rlevels(rsolverConfig%NLCOARSE)%rtriangulation, &
                                   strifile, rboundary)
          !CALL tria_readTriFile2D (Rlevels(NLMIN)%rtriangulation, &
          !                        './pre/QUAD.tri', rboundary)
          
          ! In case we have a triangular element, form a tri mesh
          if (elem_igetNVE(rsolverConfig%ielementType) .eq. 3) then
            call tria_rawGridToTri (Rlevels(rsolverConfig%NLCOARSE)%rtriangulation)
          end if
          
          ! Refine it.
          call tria_quickRefine2LevelOrdering (rsolverConfig%NLCOARSE-1,&
              Rlevels(rsolverConfig%NLCOARSE)%rtriangulation,rboundary)
          
          ! And create information about adjacencies and everything one needs from
          ! a triangulation.
          call tria_initStandardMeshFromRaw (Rlevels(rsolverConfig%NLCOARSE)%rtriangulation,&
              rboundary)
          
          ! Now refine the grid for the fine levels.
          do i = rsolverConfig%NLCOARSE+1, ilevel

            ! Refine the grid using the 2-Level-Ordering algorithm
            call tria_refine2LevelOrdering(Rlevels(i-1)%rtriangulation,&
                Rlevels(i)%rtriangulation,rboundary)
            
            ! Create a standard mesh
            call tria_initStandardMeshFromRaw(Rlevels(i)%rtriangulation,&
              rboundary)
          
          end do
          
          call stat_stopTimer (rtimerRef)

          ! Now we can start to initialise the discretisation. At first, set up
          ! a block discretisation structure that specifies the blocks in the
          ! solution vector. In this simple problem, we only have one block.
          ! Do this for all levels
          do i = rsolverConfig%NLCOARSE, ilevel
            call spdiscr_initBlockDiscr (Rlevels(i)%rdiscretisation, 1, &
                                         Rlevels(i)%rtriangulation, rboundary)
          end do
          
          ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
          ! structures for every component of the solution vector.
          ! Initialise the first element of the list to specify the element
          ! and cubature rule for this solution component:
          do i = rsolverConfig%NLCOARSE, ilevel
            call spdiscr_initDiscr_simple (&
                Rlevels(i)%rdiscretisation%RspatialDiscr(1), &
                rsolverConfig%ielementType,&
                spdiscr_getStdCubature(rsolverConfig%ielementType),&
                Rlevels(i)%rtriangulation, rboundary)
          end do
                        
          call stat_clearTimer (rtimerMatGen)
          call stat_startTimer (rtimerMatGen)
          
          ! Now as the discretisation is set up, we can start to generate
          ! the structure of the system matrix which is to solve.
          ! We create a scalar matrix, based on the discretisation structure
          ! for our one and only solution component.
          do i = rsolverConfig%NLCOARSE, ilevel

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
            !CALL lsyssc_allocEmptyMatrix (Rlevels(i)%rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_ZERO)
          
          end do
          
          call stat_stopTimer (rtimerMatGen)
            
          ! Although we could manually create the solution/RHS vector,
          ! the easiest way to set up the vector structure is
          ! to create it by using our matrix as template:
          call lsysbl_createVecBlockIndMat (Rlevels(ilevel)%rmatrix,rrhsBlock, .false.)

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
              Rlevels(ilevel)%rdiscretisation%RspatialDiscr(1),&
              rlinform,.true.,rrhsBlock%RvectorBlock(1),coeff_RHS_2D)
          
          do i = rsolverConfig%NLCOARSE, ilevel
          
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
          rrhsBlock%p_rdiscreteBC => Rlevels(ilevel)%rdiscreteBC
          
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
          
          call stat_clearTimer (rtimerSort)
          call stat_startTimer (rtimerSort)
          
          ! During the linear solver, the boundary conditions are also
          ! frequently imposed to the vectors. But as the linear solver
          ! does not work with the actual solution vectors but with
          ! defect vectors instead.
          ! So, set up a filter chain that filters the defect vector
          ! during the solution process to implement discrete boundary conditions.
          RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

          ! Resort matrices and vectors
          allocate(Rdiscretisation(ilevel))
          do i = rsolverConfig%NLCOARSE, ilevel
          
            ! For stochastic resorting, prepare an array with the discretisation structures.
            call spdiscr_duplicateDiscrSc (Rlevels(i)%rdiscretisation%RspatialDiscr(1),&
              Rdiscretisation(i), .true.)
        
            ! Create a sort strategy structure for our discretisation
            call sstrat_initBlockSorting (Rlevels(i)%rsortStrategy,Rlevels(i)%rdiscretisation)
            
            ! Calculate the resorting
            select case (rsolverConfig%isortStrategy)
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
              if (i .eq. ilevel) then
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
              if (i .eq. ilevel) then
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
              if (i .eq. ilevel) then
                call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rrhsBlock,rtempBlock2%RvectorBlock(1))
                call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rvectorBlock,rtempBlock2%RvectorBlock(1))
                call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rtempBlock,rtempBlock2%RvectorBlock(1))
              end if

            case (4)
              ! Hierarchical resorting
              call sstrat_initHierarchical (Rlevels(i)%rsortStrategy%p_Rstrategies(1),&
                  Rdiscretisation(rsolverConfig%NLCOARSE:i))

              ! Attach the sorting strategy to the matrix. The matrix is not yet sorted.
              call lsysbl_setSortStrategy (Rlevels(i)%rmatrix,&
                  Rlevels(i)%rsortStrategy,Rlevels(i)%rsortStrategy)
              
              ! Sort the matrix
              call lsysbl_sortMatrix (Rlevels(i)%rmatrix,.true.)
              
              ! Sort the vectors on the maximum level.
              if (i .eq. ilevel) then
                call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rrhsBlock,rtempBlock2%RvectorBlock(1))
                call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rvectorBlock,rtempBlock2%RvectorBlock(1))
                call lsysbl_synchroniseSort (Rlevels(i)%rmatrix,rtempBlock,rtempBlock2%RvectorBlock(1))
              end if
            end select
            
          end do
          
          call stat_stopTimer (rtimerSort)

          do i = rsolverConfig%NLCOARSE, ilevel
            call spdiscr_releaseDiscr(Rdiscretisation(i))
          end do
          deallocate(Rdiscretisation)

          ! Now we have to build up the level information for multigrid.
          !
          ! At first, initialise a standard interlevel projection structure. We
          ! can use the same structure for all levels.
          call mlprj_initProjectionMat (rprojection,Rlevels(ilevel)%rmatrix)

          ! Create a Multigrid-solver. Attach the above filter chain
          ! to the solver, so that the solver automatically filters
          ! the vector during the solution process.
          p_RfilterChain => RfilterChain
          
          select case (rsolverConfig%isolver)
          
          case (0)
          
            call linsol_initMultigrid (p_rsolverNode,p_RfilterChain)
            
            ! Set up a coarse grid solver.
            call linsol_initUMFPACK4 (p_rcoarsegridSolver)
            
            ! Add the coarse grid level.
            call linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverNode,rprojection,&
                                          null(), null(), p_rcoarseGridSolver)

            ! Now set up the other levels...
            do i = rsolverConfig%NLCOARSE+1, ilevel
            
              select case (rsolverConfig%ismoother)
              case (0)
                ! Create a Jacobi smoother
                call linsol_initJacobi(p_rsmoother)
              
                ! We will use 4 smoothing steps with damping parameter 0.7
                call linsol_convertToSmoother(p_rsmoother, rsolverConfig%nsm, rsolverConfig%domega)
                
              case (1)
                ! Create an ILU(s) smoother
                call linsol_initMILUs1x1 (p_rsmoother,rsolverConfig%ifillILU,0.0_DP)
                call linsol_convertToSmoother(p_rsmoother, rsolverConfig%nsm, rsolverConfig%domega)
              
              end select

              ! And add this multi-grid level. We will use the same smoother
              ! for pre- and post-smoothing.
              call linsol_addMultiGridLevel(p_rlevelInfo,p_rsolverNode,rprojection,&
                                            p_rsmoother, p_rsmoother, null())
              
            end do
            
            ! Initialise multigrid-specific parameters
            p_rsolverNode%p_rsubnodeMultigrid%icycle = rsolverConfig%icycle
            
          case (1)

            select case (rsolverConfig%ipreconditioner)
            case (0)
              call linsol_initJacobi (p_rsmoother,1.0_DP)
              call linsol_initCG (p_rsolverNode,p_rsmoother,p_RfilterChain)
            case (1)
              call linsol_initMILUs1x1 (p_rsmoother,rsolverConfig%nsm,0.0_DP)
              call linsol_initCG (p_rsolverNode,p_rsmoother,p_RfilterChain)
            end select
            
          case (2)

            select case (rsolverConfig%ipreconditioner)
            case (0)
              call linsol_initJacobi (p_rsmoother,1.0_DP)
              call linsol_initBiCGStab (p_rsolverNode,p_rsmoother,p_RfilterChain)
            case (1)
              call linsol_initMILUs1x1 (p_rsmoother,rsolverConfig%nsm,0.0_DP)
              call linsol_initBiCGStab (p_rsolverNode,p_rsmoother,p_RfilterChain)
            end select
            
          end select
          
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
          allocate(Rmatrices(rsolverConfig%NLCOARSE:ilevel))
          do i = rsolverConfig%NLCOARSE, ilevel
            call lsysbl_duplicateMatrix (Rlevels(i)%rmatrix,&
                Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          end do
          
          call linsol_setMatrices(p_RsolverNode,Rmatrices(rsolverConfig%NLCOARSE:ilevel))

          ! We can release Rmatrices immediately -- as long as we don't
          ! release Rlevels(i)%rmatrix!
          do i=rsolverConfig%NLCOARSE,ilevel
            call lsysbl_releaseMatrix (Rmatrices(i))
          end do
          deallocate(Rmatrices)
          
          ! Initialise structure/data of the solver. This allows the
          ! solver to allocate memory / perform some precalculation
          ! to the problem.
          call linsol_initStructure (p_rsolverNode, ierror)
          if (ierror .ne. LINSOL_ERR_NOERROR) stop
          call linsol_initData (p_rsolverNode, ierror)
          if (ierror .ne. LINSOL_ERR_NOERROR) stop
          
          ! Finally solve the system. As we want to solve Ax=b with
          ! b being the real RHS and x being the real solution vector,
          ! we use linsol_solveAdaptively. If b is a defect
          ! RHS and x a defect update to be added to a solution vector,
          ! we would have to use linsol_precondDefect instead.
          
          call stat_clearTimer(rtimerSolver)
          call stat_startTimer(rtimerSolver)
          call linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rrhsBlock,rtempBlock)
          call stat_stopTimer(rtimerSolver)
          
          ! That's it, rvectorBlock now contains our solution. We can now
          ! start the postprocessing.
          ! Start UCD export to GMV file:
          if (rsolverConfig%iwriteGMV .ne. 0) then
            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,&
                Rlevels(ilevel)%rtriangulation,'gmv/u_'//&
                trim(sys_lowcase(sconfig))//'_lv_'//trim(sys_siL(ilevel,10))//'.gmv')
            
            call lsysbl_sortVector (rvectorBlock,.false.)
            call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
            call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
            
            ! Write the file to disc, that's it.
            call ucd_write (rexport)
            call ucd_release (rexport)
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

          write (*,'(A,A,A,I3,A,E18.10,E18.10,A,I5,A,E18.10,A,E18.10,A,E18.10,A,E18.10)') &
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
          call linsol_doneData (p_rsolverNode)
          call linsol_doneStructure (p_rsolverNode)
          
          ! Release the interlevel projection structure
          call mlprj_doneProjection (rprojection)
          
          ! Release the solver node and all subnodes attached to it (if at all):
          call linsol_releaseSolver (p_rsolverNode)
          
          ! Release the block matrix/vectors
          call lsysbl_releaseVector (rtempBlock)
          call lsysbl_releaseVector (rtempBlock2)
          call lsysbl_releaseVector (rvectorBlock)
          call lsysbl_releaseVector (rrhsBlock)
          do i = ilevel, rsolverConfig%NLCOARSE, -1
            call lsysbl_releaseMatrix (Rlevels(i)%rmatrix)
            call sstrat_doneBlockSorting (Rlevels(i)%rsortStrategy)
          end do

          ! Release our discrete version of the boundary conditions
          do i = ilevel, rsolverConfig%NLCOARSE, -1
            call bcasm_releaseDiscreteBC (Rlevels(i)%rdiscreteBC)
          end do

          ! Release the discretisation structure and all spatial discretisation
          ! structures in it.
          do i = ilevel, rsolverConfig%NLCOARSE, -1
            call spdiscr_releaseBlockDiscr(Rlevels(i)%rdiscretisation)
          end do
          
          ! Release the triangulation.
          do i = ilevel, rsolverConfig%NLCOARSE, -1
            if (ilevel .eq. rsolverconfig%NLMAX) then
              call tria_infoStatistics (&
                  Rlevels(i)%rtriangulation,i .eq. rsolverConfig%NLMAX,i)
            end if
            call tria_done (Rlevels(i)%rtriangulation)
          end do
          
          deallocate(Rlevels)
                
        end do
      
      end do ! iconfig
      
      ! Finally release the domain, that's it.
      call boundary_release (rboundary)

    end if

    call parlst_done(rparlist)

  end subroutine

end module

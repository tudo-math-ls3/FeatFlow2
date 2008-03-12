!##############################################################################
!# ****************************************************************************
!# <name> heatcond_method2 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a simple heat conduction
!# problem with constant coefficients on a simple domain.
!#
!# The head conduction equation solved in this module is the simplest at all:
!#
!#    d/dt u  - Laplace u  =  f
!#
!# With the RHS f and the boundary conditions being constant in time.
!# The equation is discretised in time by explicit Euler, prividing the
!# following discrete equation:
!#
!#    (1/dt M  +  L) u_{n+1}  =  f  +  1/dt M u_n
!#
!# with M=Mass matrix, L=-Laplace, dt=time step size.
!#
!# The whole solution process is implemented into one routine, using a standard
!# linear solver (not multigrid) for solving the linear subproblems in every
!# timestep. The time domain is [0,T] and discretised by a simple explicit 
!# Euler.
!#
!# The module is an extension of the standard poisson solver.
!# Right hand side and boundary conditions are precomputed and kept constant
!# over time to reduce computational time.
!# </purpose>
!##############################################################################

MODULE heatcond_method2

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
    
  USE heatcond_callback
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE heatcond2
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  !  1.) Read in parametrisation
  !  2.) Read in triangulation
  !  3.) Set up boundary conditions
  !  4.) Set up matrices
  !  5.) Set up RHS
  !  6.) Create solver structure
  !  7.) Loop over all timesteps
  !  8.)  Solve the problem
  !  9.)  Write solution to GMV file
  ! 10.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An object for saving the domain:
    TYPE(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    TYPE(t_blockDiscretisation) :: rdiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    ! A scalar mass matrix, Laplace matrix, temporary matrix and
    ! a vector. The vector accepts the RHS of the problem in scalar form.
    TYPE(t_matrixScalar) :: rmatrixLaplace, rmatrixMass, rmatrix
    TYPE(t_vectorScalar) :: rrhs

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_matrixBlock) :: rmatrixBlock,rmatrixMassBlock
    TYPE(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock,rtimeRhsBlock

    ! A set of variables describing the analytic and discrete boundary
    ! conditions.    
    TYPE(t_boundaryRegion) :: rboundaryRegion
    TYPE(t_discreteBC), TARGET :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
    
    ! NLMAX receives the level where we want to solve.
    INTEGER :: NLMAX
    
    ! Time step size, number of timesteps.
    REAL(DP) :: dtstep
    INTEGER :: ntimesteps
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
    
    ! Time and time step counter
    REAL(DP) :: dtime
    INTEGER :: itimestep
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata

    ! Ok, let's start. 
    !
    ! We want to solve our heat equation problem on level...
    NLMAX = 6
    
    ! Initialise time step size and number of timesteps
    dtstep = 0.01_DP
    ntimesteps = 10

    ! We start at time 0.0.
    dtime = 0.0_DP
    
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    CALL boundary_read_prm(rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rtriangulation, './pre/QUAD.tri', rboundary)
    
    ! Refine it.
    CALL tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    CALL tria_initStandardMeshFromRaw (rtriangulation,rboundary)
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    CALL spdiscr_initBlockDiscr2D (rdiscretisation,1,&
                                   rtriangulation, rboundary)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    CALL spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscretisation(1), &
                                   EL_E011,CUB_TRZ,rtriangulation, rboundary)
                                   !CUB_G2X2

    ! Up to now, everything is 'analytical'.
    ! Let's change that, let's start to discretise!
    !
    ! 1.) Boundary conditions
    !
    ! For implementing boundary conditions, we use a 'filter technique with
    ! discretised boundary conditions'. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    CALL bcasm_initDiscreteBC(rdiscreteBC)
    !
    ! We 'know' already (from the problem definition) that we have four boundary
    ! segments in the domain. Each of these, we want to use for enforcing
    ! some kind of boundary condition.
    !
    ! We ask the bondary routines to create a 'boundary region' - which is
    ! simply a part of the boundary corresponding to a boundary segment.
    ! A boundary region roughly contains the type, the min/max parameter value
    ! and whether the endpoints are inside the region or not.
    CALL boundary_createRegion(rboundary,1,1,rboundaryRegion)
    
    ! We use this boundary region and specify that we want to have Dirichlet
    ! boundary there. The following call does the following:
    ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
    !   We specify icomponent='1' to indicate that we set up the
    !   Dirichlet BC's for the first (here: one and only) component in the 
    !   solution vector.
    ! - Discretise the boundary condition so that the BC's can be applied
    !   to matrices and vectors
    ! - Add the calculated discrete BC's to rdiscreteBC for later use.
    CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)
                              
    ! Now to the edge 2 of boundary component 1 the domain.
    CALL boundary_createRegion(rboundary,1,2,rboundaryRegion)
    CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)
                              
    ! Edge 3 of boundary component 1.
    CALL boundary_createRegion(rboundary,1,3,rboundaryRegion)
    CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)
    
    ! Edge 4 of boundary component 1. That's it.
    CALL boundary_createRegion(rboundary,1,4,rboundaryRegion)
    CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)

    ! 2.) Matrices
    !
    ! At first: create the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    CALL bilf_createMatrixStructure (rdiscretisation%RspatialDiscretisation(1),&
                                     LSYSSC_MATRIX9,rmatrix)
                                     
    ! For the time dependent problem, we need a Laplace matrix and a
    ! Mass matrix. Both have exactly the same structure as rmatrix, but they
    ! have different entries!
    ! We 'copy' the structure of rmatrix to the Mass and Laplace matrix in that
    ! sense, that they 'share' their structure with rmatrix. This means, all
    ! three matrices rmatrix, rmatrixMass and rmatrixLaplace have exactly
    ! the same structure. The arrays in memory defining that structure is
    ! shared among the matrices. This helps to save memory!
    CALL lsyssc_duplicateMatrix (rmatrix,rmatrixLaplace,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)

    CALL lsyssc_duplicateMatrix (rmatrix,rmatrixMass,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
        
    ! Now to the entries of the Laplace matrix. For assembling of the entries,
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
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will use constant coefficients.
    CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrixLaplace)
    
    ! The same way, we create a mass matrix.
    ! We specify the bilinear form (Psi_j, Phi_i) for the
    ! scalar system matrix in 2D.
    
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    rform%Dcoefficients(1)  = 1.0 
    rform%Dcoefficients(2)  = 1.0 

    ! Build the entries.
    CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrixMass)

    ! Allocate memory for a temporary matrix. 
    CALL lsyssc_allocEmptyMatrix(rmatrix,LSYSSC_SETM_ZERO)
    
    ! From the temp matrix and from the mass matrix, we derive 'block' versions.
    ! These are 1x1 block matrices that share their data with the mass matrix
    ! and the temporary matrix, respectively. This helps to deal more easily with
    ! the linear solver later.
    CALL lsysbl_createMatFromScalar (rmatrix,rmatrixBlock,rdiscretisation)
    CALL lsysbl_createMatFromScalar (rmatrixMass,rmatrixMassBlock,rdiscretisation)
    
    ! Hang the pointer of the boudnary conditions into the matrix. That way, these
    ! boundary conditions are always connected to that matrix.
    rmatrixBlock%p_rdiscreteBC => rdiscreteBC
    
    ! To give an overview, we now have (concerning the matrices):
    !
    ! - Laplace matrix -> rmatrixLaplace
    ! - Mass matrix    -> rmatrixMass and rmatrixMassBlock
    ! - Temp matrix    -> rmatrix     and rmatrixBlock (knowing the BC's)
    !
    ! 3.) Vectors
    !
    ! Now we come to the stuff with the vectors. Thís starts with the RHS vector.
    !
    ! For setting up up a RHS vector during the time stepping, we have to
    ! set up the description of the corresponding linear form.
    ! The same has to be done for the right hand side of the problem.
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to get a discrete version of it.
    ! Again we simply create a scalar vector based on the one and only
    ! discretisation structure.
    ! This scalar vector will later be used as the one and only first
    ! component in a block vector.
    CALL linf_buildVectorScalar (rdiscretisation%RspatialDiscretisation(1),&
                                 rlinform,.TRUE.,rrhs,coeff_RHS)

    ! Immediately create a 1-block version of the scalar RHS.
    ! This is then compatible to the 1x1 block matrix rmatrixBlock
    CALL lsysbl_createVecFromScalar (rrhs,rrhsBlock,rdiscretisation)
    
    ! Hang the pointer into the vector. That way, these
    ! boundary conditions are always connected to that vector.
    rrhsBlock%p_rdiscreteBC => rdiscreteBC
                             
    ! We have a block vector for the RHS, but that's not enough. What we
    ! need in total is
    ! - A RHS vector (we have)
    ! - A RHS vector for the solver (changing in every time step, missing)
    ! - A solution vector (changing in every time step, missing)
    ! Allocate memory for the missing vectors. They have exactly the same shape as
    ! our RHS vector, so we can use the RHS to 'derive' the missing vectors.
    ! Note that this will also transfer the connected boudary conditions
    ! to the new vectors.
    ! Btw., the (initial) solution vector is filled with 0.
    CALL lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .TRUE.)
    CALL lsysbl_createVecBlockIndirect (rrhsBlock, rtimeRhsBlock, .FALSE.)
    CALL lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .FALSE.)
    
    ! To give an overview, we now have (concerning the matrices):
    !
    ! - RHS vector "f"            -> rrhs and rrhsBlock (knowing the BC's)
    ! - RHS vector for the solver ->          rtimeRhsBlock (knowing the BC's)
    ! - Solution vector           ->          rvectorBlock (knowing the BC's)
    !
    ! Discretisation Finished!
    
    ! Next step: Set up the solvers!
    !
    ! During the linear solver, the boundary conditions are also
    ! frequently imposed to the vectors. But as the linear solver
    ! does not work with the actual solution vectors but with
    ! defect vectors instead.
    ! So, set up a filter chain that filters the defect vector
    ! during the solution process to implement discrete boundary conditions.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Create a BiCGStab-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    NULLIFY(p_rpreconditioner)
    CALL linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    ! Ok, now it gets interesting. The linear solver needs a matrix.
    ! Using the mass and Laplace matrix, we set up
    !      A = 1/dt M + L
    ! in rmatrix. This is the system matrix, which is valid in all time steps
    ! as we use explicit Euler with a fixed time step!
    !
    ! The boolean parameters in this call say:
    ! - Sum up the entries from rmatrixMass and rmatrixLaplace.
    ! - All matrices have the same structure, even the destination matrix.
    !   So it's not necessary to allocate any memory.
    !
    ! Note that rmatrix shares its data with rmatrixBlock(1,1) so by modifying
    ! rmatrix, we simultaneously modify rmatrixBlock(1,1).
    CALL lsyssc_matrixLinearComb (&
        rmatrixMass,1.0_DP/dtstep, rmatrixLaplace,1.0_DP, rmatrix,&
        .FALSE.,.FALSE.,.TRUE.,.TRUE.)
        
    ! Next step is to implement boundary conditions into the matrix.
    ! This is done using a vector/matrix filter for discrete boundary 
    ! conditions.
    ! The discrete boundary conditions are already attached to the
    ! matrix. Call the appropriate matrix filter that modifies the 
    ! matrix according to the boundary conditions.
    CALL matfil_discreteBC (rmatrixBlock)

    ! Attach the system matrix to the solver.
    ! First create an array with the matrix data (on all levels, but we
    ! only have one level here), then call the initialisation 
    ! routine to attach all these matrices.
    ! Remark: Don't make a call like
    !    CALL linsol_setMatrices(p_RsolverNode,(/rmatrixBlock/))
    ! This doesn't work on all compilers, since the compiler would have
    ! to create a temp array on the stack - which does not always work!
    Rmatrices = (/rmatrixBlock/)
    CALL linsol_setMatrices(p_RsolverNode,Rmatrices)
    
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    CALL linsol_initStructure (p_rsolverNode, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    CALL linsol_initData (p_rsolverNode, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    
    ! Start the time-loop.
    DO itimestep = 1,ntimesteps
      
      ! Next time step.
      dtime = dtime + dtstep
      
      CALL output_separator(OU_SEP_MINUS)
      CALL output_line ('Time step '//TRIM(sys_siL(itimestep,6))// &
                        '     Time '//TRIM(sys_sdL(dtime,5)))
      CALL output_lbrk ()

      ! Create the new RHS for the solver. To do this, calculate:
      !
      !    RHS = f + 1/dt M u

      CALL lsysbl_copyVector (rrhsBlock, rtimeRhsBlock)
      CALL lsysbl_blockMatVec (rmatrixMassBlock, rvectorBlock, rtimeRhsBlock, &
          1.0_DP/dtstep, 1.0_DP)
          
      ! Implement the boundary conditions into the RHS and into the solution vector.
      CALL vecfil_discreteBCrhs (rtimeRhsBlock)
      CALL vecfil_discreteBCsol (rvectorBlock)
      
      ! Solve the system for the new time step:
      !
      !    (M/dt + L) u  =  f + M/dt u_{old}
      !    ^^^^^^^^^^       ^^^^^^^^^^^^^^^^
      !   rmatrixBlock       rtimeRhsBlock
      !
      ! <=>         A x  = b    
      !
      ! As we want to solve Ax=b with b being the real RHS and 
      ! x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
      CALL linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rtimeRhsBlock,rtempBlock)
      
      ! That's it, rvectorBlock now contains our solution. We can now
      ! start the postprocessing. 
      ! Start UCD export to GMV file:
      CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                        'gmv/u2.gmv.'//TRIM(sys_si0L(itimestep,5)))
      
      CALL lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
      CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
      
      ! Write the file to disc, that's it.
      CALL ucd_write (rexport)
      CALL ucd_release (rexport)
      
    END DO
    
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    CALL linsol_doneData (p_rsolverNode)
    CALL linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    CALL lsysbl_releaseVector (rtempBlock)
    CALL lsysbl_releaseVector (rvectorBlock)
    CALL lsysbl_releaseVector (rtimeRhsBlock)
    CALL lsysbl_releaseVector (rrhsBlock)
    CALL lsysbl_releaseMatrix (rmatrixBlock)
    CALL lsysbl_releaseMatrix (rmatrixMassBlock)

    ! Release the matrices -- after releasing the block matrices!
    ! These must exist as long as the block matrices exist, as the block 
    ! matrices/vectors are only 'copies' of the scalar ones, sharing the 
    ! same handles!
    CALL lsyssc_releaseVector (rrhs)
    CALL lsyssc_releaseMatrix (rmatrix)
    CALL lsyssc_releaseMatrix (rmatrixLaplace)
    CALL lsyssc_releaseMatrix (rmatrixMass)
    
    ! Release our discrete version of the boundary conditions
    CALL bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    CALL spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    CALL tria_done (rtriangulation)
    
    ! Finally release the domain, that's it.
    CALL boundary_release (rboundary)
    
  END SUBROUTINE

END MODULE

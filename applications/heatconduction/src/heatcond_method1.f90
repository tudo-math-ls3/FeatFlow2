!##############################################################################
!# ****************************************************************************
!# <name> heatcond_method1 </name>
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
!# The module bases on the standard poisson example and shows which
!# modifications are necessary to create a heat equation solver from a poisson
!# solver. Boundary conditions, matrices and vectors are all reassembled in 
!# every timestep
!# </purpose>
!##############################################################################

MODULE heatcond_method1

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

  SUBROUTINE heatcond1
  
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
    
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    TYPE(t_matrixScalar) :: rmatrix
    TYPE(t_vectorScalar) :: rrhs,rvector

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_matrixBlock) :: rmatrixBlock
    TYPE(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

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
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    ! Time step size, number of timesteps.
    REAL(DP) :: dtstep
    INTEGER :: ntimesteps

    ! Time and time step counter
    REAL(DP) :: dtime
    INTEGER :: itimestep
    
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
                                   EL_E011,CUB_G2X2,rtriangulation, rboundary)
                                   
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    CALL bilf_createMatrixStructure (rdiscretisation%RspatialDiscretisation(1),&
                                     LSYSSC_MATRIX9,rmatrix)
                                     
    ! Use the discretisation to create a solution vector.
    ! Fill it with zero -- that's out initial condition!
    CALL lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscretisation(1),&
        rvector,.TRUE.)             
    
    ! Start the timeloop
    DO itimestep=1,ntimesteps
       
      ! Next time step.
      dtime = dtime + dtstep
      
      CALL output_separator(OU_SEP_MINUS)
      CALL output_line ('Time step '//TRIM(sys_siL(itimestep,6))// &
                        '     Time '//TRIM(sys_sdL(dtime,5)))
      CALL output_lbrk ()

      ! STEP 1: Form the right hand side:  dtimestep*f + M u_{old}
      !
      ! To assemble the basic RHS f, set up the corresponding linear 
      ! form (f,Phi_j):
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC
      
      ! Discretise the RHS. We simply create a scalar vector 
      ! based on the one and only discretisation structure.
      ! The result is rrhs!
      CALL linf_buildVectorScalar (rdiscretisation%RspatialDiscretisation(1),&
                                  rlinform,.TRUE.,rrhs,coeff_RHS)   
                                  

      ! And now to the entries of the mass matrix. 
      ! For assembling of the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (Psi_j, Phi_i) for the
      ! scalar system matrix in 2D.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC
      
      rform%ballCoeffConstant = .TRUE.
      rform%BconstantCoeff = .TRUE.
      rform%Dcoefficients(1)  = 1.0 
      CALL bilf_buildMatrixScalar (rform,.TRUE.,rmatrix)
      
      ! Now form the actual RHS by matrix vector multiplication!
      ! dtimestep*f + M u_{old}
      CALL lsyssc_scalarMatVec(rmatrix,rvector,rrhs,1.0_DP,dtstep)
      
      ! STEP 2: Assemble the system matrix (M + dtimestep*Laplace)
      !
      ! For this purpose, set up the corresponding bilinear form
      ! (grad Psi_j, grad Phi_i):
      rform%itermCount = 2
      rform%Idescriptors(1,1) = DER_DERIV_X
      rform%Idescriptors(2,1) = DER_DERIV_X
      rform%Idescriptors(1,2) = DER_DERIV_Y
      rform%Idescriptors(2,2) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients.
      rform%ballCoeffConstant = .TRUE.
      rform%BconstantCoeff = .TRUE.
      rform%Dcoefficients(1)  = dtstep
      rform%Dcoefficients(2)  = dtstep 

      ! Now we can build the matrix entries.
      ! Note that we set bclear=FALSE in this call, so the Laplace part
      ! is added to the existing (!) mass matrix!
      ! So this results in (M + dtstep*Laplace), as we set 
      ! the coefficient rform%Dcoefficients to dtstep above!
      CALL bilf_buildMatrixScalar (rform,.FALSE.,rmatrix)
      
      ! STEP 3: Create block vectors and boundary conditions.
      !      
      ! The linear solver only works for block matrices/vectors - but above,
      ! we created scalar ones. So the next step is to make a 1x1 block
      ! system from the matrices/vectors above which the linear solver
      ! understands.
      CALL lsysbl_createMatFromScalar (rmatrix,rmatrixBlock,rdiscretisation)
      CALL lsysbl_createVecFromScalar (rrhs,rrhsBlock,rdiscretisation)
      CALL lsysbl_createVecFromScalar (rvector,rvectorBlock,rdiscretisation)
      
      ! Now we have the raw problem. What is missing is the definition of the boudary
      ! conditions.
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

      ! Hang the pointer into the vector and matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      rmatrixBlock%p_rdiscreteBC => rdiscreteBC
      rrhsBlock%p_rdiscreteBC => rdiscreteBC
      rvectorBlock%p_rdiscreteBC => rdiscreteBC
      
      ! Now we have block vectors for the RHS and the matrix. What we
      ! need additionally is a block vector for the solution and
      ! temporary data. Create them using the RHS as template.
      ! Fill the solution vector with 0:
      !CALL lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .TRUE.)
      CALL lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .FALSE.)
      
      ! Next step is to implement boundary conditions into the RHS,
      ! solution and matrix. This is done using a vector/matrix filter
      ! for discrete boundary conditions.
      ! The discrete boundary conditions are already attached to the
      ! vectors/matrix. Call the appropriate vector/matrix filter that
      ! modifies the vectors/matrix according to the boundary conditions.
      CALL vecfil_discreteBCrhs (rrhsBlock)
      CALL vecfil_discreteBCsol (rvectorBlock)
      CALL matfil_discreteBC (rmatrixBlock)
      
      ! STEP 6: Solve the system
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
      
      ! Attach the system matrix to the solver.
      ! First create an array with the matrix data (on all levels, but we
      ! only have one level here), then call the initialisation 
      ! routine to attach all these matrices.
      ! Remark: Don't make a call like
      !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
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
      
      ! Finally solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
      CALL linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rrhsBlock,rtempBlock)
      
      ! STEP 7: Postprocessing
      !
      ! That's it, rvectorBlock now contains our solution. We can now
      ! start the postprocessing. 
      ! Start UCD export to GMV file:
      CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                        'gmv/u1.gmv.'//TRIM(sys_si0L(itimestep,5)))
      
      CALL lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
      CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
      
      ! Write the file to disc, that's it.
      CALL ucd_write (rexport)
      CALL ucd_release (rexport)
      
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
      CALL lsysbl_releaseVector (rrhsBlock)
      CALL lsysbl_releaseMatrix (rmatrixBlock)

      ! Release the scalar matrix/rhs vector which were used to create
      ! the block matrices/vectors. These must exist as long as the
      ! block matrices/vectors exist, as the block matrices/vectors are
      ! only 'copies' of the scalar ones, sharing the same handles!
      CALL lsyssc_releaseVector (rrhs)
      
      ! Release our discrete version of the boundary conditions
      CALL bcasm_releaseDiscreteBC (rdiscreteBC)

    END DO
    
    ! Release the preallocated matrix and the solution vector.
    CALL lsyssc_releaseMatrix (rmatrix)
    CALL lsyssc_releaseVector (rvector)
    
    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    CALL spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    CALL tria_done (rtriangulation)
    
    ! Finally release the domain, that's it.
    CALL boundary_release (rboundary)
    
  END SUBROUTINE

END MODULE

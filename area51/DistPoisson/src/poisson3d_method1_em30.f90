!##############################################################################
!# ****************************************************************************
!# <name> poisson3d_method1_em30 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!# This module is a (provisional) equivalent to the 2D example
!# poisson2d_method7.
!# The element used here is EM30. For GMV output, it is interpolated
!# to Q1 and then written to the output file.
!# </purpose>
!##############################################################################

module poisson3d_method1_em30

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use linearsystemscalar
  use linearsystemblock
  use filtersupport
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use scalarpde
  use element
  use ucd
  use pprocerror
  use genoutput
  use matrixio
  use meshregion
    
  use poisson3d_callback
  use spdiscprojection
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine poisson3d_1_em30
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.) Read in triangulation
  ! 2.) Set up RHS
  ! 3.) Set up matrix
  ! 4.) Create solver structure
  ! 5.) Solve the problem
  ! 6.) Convert the solution vector into GMV`s style
  ! 7.) Write solution to GMV file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let us see...
    !
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! An object for saving the boundary mesh region
    type(t_meshregion) :: rmeshRegion

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar) :: rmatrix
    type(t_vectorScalar) :: rrhs

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlock
    type(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock

    ! A variable describing the discrete boundary conditions.
    type(t_discreteBC), target :: rdiscreteBC

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(1) :: Rmatrices

    ! A filter chain that describes how to filter the matrix/vector
    ! before/during the solution process. The filters usually implement
    ! boundary conditions.
    type(t_filterChain), dimension(1), target :: RfilterChain
    
    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Declarations for projecting a vector to the Q1 space for GMV export
    type(t_vectorBlock) :: rprjVector
    type(t_blockDiscretisation) :: rprjDiscretisation
    type(t_discreteBC), target :: rdiscreteBC_Q1

    ! Ok, let us start.
    !
    ! We want to solve our Poisson problem on level...
    NLMAX = 4
    
    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './pre'

    ! At first, read in the basic triangulation.
    call tria_readTriFile3D (rtriangulation, trim(spredir)//'/CUBE.tri')
    
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation)
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr (rdiscretisation,1,rtriangulation)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   EL_EM30_3D,CUB_G3_3D,rtriangulation)
                 
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrix)
    
    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
    ! scalar system matrix in 3D.
    rform%itermCount = 3
    rform%Idescriptors(1,1) = DER_DERIV3D_X
    rform%Idescriptors(2,1) = DER_DERIV3D_X
    rform%Idescriptors(1,2) = DER_DERIV3D_Y
    rform%Idescriptors(2,2) = DER_DERIV3D_Y
    rform%Idescriptors(1,3) = DER_DERIV3D_Z
    rform%Idescriptors(2,3) = DER_DERIV3D_Z

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    rform%Dcoefficients(1)  = 1.0
    rform%Dcoefficients(2)  = 1.0
    rform%Dcoefficients(3)  = 1.0

    ! Now we can build the matrix entries.
    ! We specify the callback function coeff_Laplace for the coefficients.
    ! As long as we use constant coefficients, this routine is not used.
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will call the callback routine to get analytical
    ! data.
    call bilf_buildMatrixScalar (rform,.true.,rmatrix,coeff_Laplace_3D)
                                  
    ! The same has to be done for the right hand side of the problem.
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC3D
    
    ! ... and then discretise the RHS to get a discrete version of it.
    ! Again we simply create a scalar vector based on the one and only
    ! discretisation structure.
    ! This scalar vector will later be used as the one and only first
    ! component in a block vector.
    call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
                                 rlinform,.true.,rrhs,coeff_RHS_3D)
    
    ! The linear solver only works for block matrices/vectors - but above,
    ! we created scalar ones. So the next step is to make a 1x1 block
    ! system from the matrices/vectors above which the linear solver
    ! understands.
    call lsysbl_createMatFromScalar (rmatrix,rmatrixBlock,rdiscretisation)
    call lsysbl_createVecFromScalar (rrhs,rrhsBlock,rdiscretisation)
    
    ! Now we have the raw problem. What is missing is the definition of the boundary
    ! conditions.
    ! For implementing boundary conditions, we use a `filter technique with
    ! discretised boundary conditions`. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    
    ! In contrast to the 2D examples, we currently do not have an analytic
    ! description of the domain`s boundary, therefore we need a discrete
    ! (mesh-dependent) description of the mesh`s boundary. This can be done
    ! using mesh-regions.
    !
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    call bcasm_initDiscreteBC(rdiscreteBC)
    
    ! Create a mesh region describing the mesh`s boundary based on the
    ! nodal-property-array of the current triangulation.
    call mshreg_createFromNodalProp(rmeshRegion, rtriangulation, &
                                      MSHREG_IDX_ALL)
    
    ! Describe Dirichlet BCs on that mesh region
    call bcasm_newDirichletBConMR(rdiscretisation, 1, rdiscreteBC, rmeshRegion,&
                                  getBoundaryValuesMR_3D)
    
    ! We will not free the mesh region structure here, as we will need it
    ! later for calculating the Q1 solution.

    ! Hang the pointer into the vector and matrix. That way, these
    ! boundary conditions are always connected to that matrix and that
    ! vector.
    rmatrixBlock%p_rdiscreteBC => rdiscreteBC
    rrhsBlock%p_rdiscreteBC => rdiscreteBC
                             
    ! Now we have block vectors for the RHS and the matrix. What we
    ! need additionally is a block vector for the solution and
    ! temporary data. Create them using the RHS as template.
    ! Fill the solution vector with 0:
    call lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .true.)
    call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .false.)
    
    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
    call vecfil_discreteBCrhs (rrhsBlock)
    call vecfil_discreteBCsol (rvectorBlock)
    call matfil_discreteBC (rmatrixBlock)

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
    nullify(p_rpreconditioner)
    !CALL linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
    call linsol_initUMFPACK4 (p_rsolverNode)
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    
    ! Attach the system matrix to the solver.
    ! First create an array with the matrix data (on all levels, but we
    ! only have one level here), then call the initialisation
    ! routine to attach all these matrices.
    ! Remark: Do not make a call like
    !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
    ! This does not work on all compilers, since the compiler would have
    ! to create a temp array on the stack - which does not always work!
    Rmatrices = (/rmatrixBlock/)
    call linsol_setMatrices(p_RsolverNode,Rmatrices)
    
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
    call linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rrhsBlock,rtempBlock)
    
    ! -------------------------------------------------------------------------
    ! Solution vector projection for projecting a solution vector of
    ! an arbitrary element to Q1 such that it can be written into a GMV file
    ! -------------------------------------------------------------------------

    ! Now, Our vector block is off the way GMV`s style. We need to convert the vector
    ! to Q1 as that is the format, GMV`s understands. So the task is to
    ! create a Q1 solution from rvectorBlock.
    !
    ! Step 1: Create a discretisation structure for Q1, based on our
    ! previous discretisation structure:
    call spdiscr_duplicateBlockDiscr(rvectorBlock%p_rblockDiscr,rprjDiscretisation)
    call spdiscr_deriveSimpleDiscrSc (&
                 rvectorBlock%p_rblockDiscr%RspatialDiscr(1), &
                 EL_Q1_3D, CUB_G3_3D, rprjDiscretisation%RspatialDiscr(1))

    ! Step 2: Setup a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)
 
    ! Step 3: Take the original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution(rvectorBlock,rprjVector)

    ! Step 4: Discretise the boundary condition according to the Q1
    ! discretisation:
    call bcasm_initDiscreteBC(rdiscreteBC_Q1)
    
    ! Describe Dirichlet BCs on the mesh region which was created before
    call bcasm_newDirichletBConMR(rprjDiscretisation, 1, rdiscreteBC_Q1, rmeshRegion,&
                                  getBoundaryValuesMR_3D)
    
    ! Free the mesh region structure as we will not need it anymore
    call mshreg_done(rmeshRegion)

    ! Connect the vector to the BC`s
    rprjVector%p_rdiscreteBC => rdiscreteBC_Q1

    ! Step 5: Set up a boundary condition filter for Dirichtley boundary conditions
    ! and pass the vetor through it. This finally implement the Dirichtley
    ! boundary conditions into the output vector
    RfilterChain(1)%ifilterType = FILTER_DISCBCSOLREAL
    call filter_applyFilterChainVec (rprjVector, RfilterChain)

    ! Now we have a Q1 solution in rprjVector.
    !
    ! Step 6: Write the GMV file.
    ! That is it, rvectorBlock now contains our solution in GMV`s style.
    ! We can now start the postprocessing.
      
    ! Call the GMV library to write out a GMV file for our solution.
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       trim(sucddir)//'/u3d_1_em30.gmv')
    
    call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! Step 7: Release our discrete version of the boundary conditions of Q1.
    call bcasm_releaseDiscreteBC(rdiscreteBC_Q1)
    
    ! Release the block matrix/vectors rprjVector
    call lsysbl_releaseVector(rprjVector)
    
    ! Release the discretisation structure.
    call spdiscr_releaseBlockDiscr (rprjDiscretisation)
    
    ! -------------------------------------------------------------------------
    ! Projection and GMV export finished.
    ! -------------------------------------------------------------------------

    ! Calculate the error to the reference function.
    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,derror,&
                       getReferenceFunction_3D)
    call output_line ('L2-error: ' // sys_sdEL(derror,10) )

    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
                       getReferenceFunction_3D)
    call output_line ('H1-error: ' // sys_sdEL(derror,10) )
    
    ! We are finished - but not completely!
    ! Now, clean up so that all the memory is available again.
    !
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the block matrix/vectors
    call lsysbl_releaseVector (rtempBlock)
    call lsysbl_releaseVector (rvectorBlock)
    call lsysbl_releaseVector (rrhsBlock)
    call lsysbl_releaseMatrix (rmatrixBlock)

    ! Release the scalar matrix/rhs vector which were used to create
    ! the block matrices/vectors. These must exist as long as the
    ! block matrices/vectors exist, as the block matrices/vectors are
    ! only 'copies' of the scalar ones, sharing the same handles!
    call lsyssc_releaseVector (rrhs)
    call lsyssc_releaseMatrix (rmatrix)
    
    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation.
    call tria_done (rtriangulation)
    
  end subroutine

end module

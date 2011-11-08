!##############################################################################
!# ****************************************************************************
!# <name> poisson1d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple 1D Poisson
!# problem with constant coefficients on a simple domain.
!# </purpose>
!##############################################################################

module poisson1d_method0_simple

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
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use scalarpde
  use ucd
  use element
  use pprocerror
  use genoutput
  use matrixio
  use vectorio
  use meshregion
  use discretebc
    
  use dg1d_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine poisson1d_0_simple
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.) Create triangulation
  ! 2.) Set up RHS
  ! 3.) Set up matrix
  ! 4.) Create solver structure
  ! 5.) Solve the problem
  ! 6.) Write solution to GMV file
  ! 7.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let us see...
    !
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation
    
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

    ! A set of variables describing the analytic and discrete boundary
    ! conditions.
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
    
    ! Error indicator during initialisation of the solver
    integer :: ierror
    
    ! Error of FE function to reference function
    real(DP) :: derror

    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata

    ! The number of sub-intervals for the discretisation
    integer :: nintervals = 16
    
    ! Ok, let us start.
    ! At first, create the basic triangulation.
    ! Our domain is [0, 1], divided into nintervals sub-intervals.
    call tria_createRawTria1D(rtriangulation, 0.0_DP, 1.0_DP, nintervals)
    
    ! As the tria_createRawTria1D routine always generates a grid
    ! with sub-intervals of equal length, we can optionally disturb
    ! the mesh.
    !CALL meshmod_disturbMesh(rtriangulation, 0.2_DP)

    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation)
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr (rdiscretisation, 1, rtriangulation)
    
    ! In the next step, we will define the element type and the cubature
    ! formula that is to be used. For this 1D poisson-example we currently
    ! have 2 possible element types: linear and quadratic ones.
    ! For linear elements the trapezoidal formula satisfies all our
    ! needs, for quadratic elements we should choose a 3-point Gauss-formula.
    !
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
    ! Setting up a linear element and trapezoidal rule would be...
                                   EL_P1_1D,CUB_TRZ_1D,rtriangulation)
    ! Setting up a quadratic element and 3-point Gauss rule would be...
                                   !EL_P2_1D,CUB_G3_1D,rtriangulation)
    ! Setting up a cubic spline element and 4-point Gauss rule would be...
                                   !EL_S31_1D,CUB_G4_1D,rtriangulation)
                                   
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrix)
    
    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
    ! scalar system matrix in 1D.
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X

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
    call bilf_buildMatrixScalar (rform,.true.,rmatrix,coeff_Laplace_1D)
    
    ! The same has to be done for the right hand side of the problem.
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to get a discrete version of it.
    ! Again we simply create a scalar vector based on the one and only
    ! discretisation structure.
    ! This scalar vector will later be used as the one and only first
    ! component in a block vector.
    call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(1),&
                                 rlinform,.true.,rrhs,coeff_RHS_1D)
    
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
    !
    ! In contrast to the 2D poisson examples, we will directly set the
    ! dirichlet boundary conditions by hand instead of discretising an analytic
    ! boundary condition function using a boundary structure.
    !
    ! Initialise the structure that collects the discrete BC`s:
    call bcasm_initDiscreteBC(rdiscreteBC)

    ! In 1D we have 2 possibilities to describe Dirichlet BCs on the interval
    ! ends. One possibility is to use the bcasm_initDirichletBC_1D routine.
    ! The following call would prescribe 0 on both interval ends:
    !
    ! CALL bcasm_newDirichletBC_1D(rdiscretisation, rdiscreteBC, 0.0_DP, 0.0_DP)
    !
    ! The second possibility is using mesh regions:
    !
    ! Create a mesh region describing the mesh`s boundary based on the
    ! nodal-property-array of the current triangulation.
    call mshreg_createFromNodalProp(rmeshRegion, rtriangulation, &
                                      MSHREG_IDX_ALL)
    
    ! Describe Dirichlet BCs on that mesh region
    call bcasm_newDirichletBConMR(rdiscretisation, 1, rdiscreteBC, rmeshRegion,&
                                  getBoundaryValuesMR_1D)
    
    ! Free the mesh region structure as we will not need it anymore
    call mshreg_done(rmeshRegion)
    
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

    ! Attach the above filter chain to the solver, so that the solver
    ! automatically filters the vector during the solution process!
    !
    ! We now have the option to create a preconditioner for the solver.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Preconditioner Remark:
    ! ----------------------
    ! Please keep in mind that if the 1D grid for the discretisation was
    ! created using the tria_createRawTria1D routine (which is the default
    ! case in this example) and the grid was NOT refined afterwards (which is
    ! also the default case in this example), the resulting matrix will be
    ! tri-diagonal for linear elements.
    ! In this case, a LU-decomposition would produce no fill-in - therefore the
    ! ILU(0) preconditioner computes a complete LU-decomposition of the system
    ! matrix. So do not wonder if you set ILU(0) as a preconditioner and the
    ! solver always converges after 1 iteration... ^_^
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    nullify(p_rpreconditioner)
    ! Setting up a Jacobi-Preconditioner would be..
    !CALL linsol_initJacobi (p_rpreconditioner)
    ! Setting up a Jin-Wei-Tam-Preconditioner would be...
    !CALL linsol_initJinWeiTam (p_rpreconditioner)
    ! Setting up a SOR[1.2]-Preconditioner would be...
    !CALL linsol_initSOR (p_rpreconditioner, 1.2_DP)
    ! Setting up a SSOR[1.2]-Preconditioner would be...
    !CALL linsol_initSSOR (p_rpreconditioner, 1.2_DP)
    ! Setting up a ILU(0)-Preconditioner would be...
    !CALL linsol_initMILUs1x1(p_rpreconditioner, 0, 0.0_DP)
    
    ! We now need to create a solver for the linear system.
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Solver Remark:
    ! --------------
    ! Please keep in mind that the CG solver needs a symmetric preconditioner
    ! and will (most probably) not work with unsymmetric preconditioners as
    ! SOR or (M)ILU(s).
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Setting up a Defect-Correction-Solver would be...
    !CALL linsol_initDefCorr (p_rsolverNode,p_rpreconditioner,RfilterChain)
    ! Setting up a PCG-Solver would be...
    !CALL linsol_initCG (p_rsolverNode,p_rpreconditioner,RfilterChain)
    ! Setting up a BiCGStab-Solver would be...
    call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
    ! Setting up a GMRES(17)-Solver would be...
    !CALL linsol_initGMRES (p_rsolverNode,17,p_rpreconditioner,RfilterChain)
    
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
    
    ! That is it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing.
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       trim(sucddir)//'/u1d_0_simple.gmv')
    
    call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Calculate the error to the reference function.
    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,derror,&
                       getReferenceFunction_1D)
    call output_line ('L2-error: ' // sys_sdEL(derror,10) )
    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
                       getReferenceFunction_1D)
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

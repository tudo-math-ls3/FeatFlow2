!##############################################################################
!# ****************************************************************************
!# <name> poisson2d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Poisson
!# problem with constant coefficients on a simple domain.
!# </purpose>
!##############################################################################

module dg2d_method0_simple

  use fsystem
  use stdoperators
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
  use element
  use ucd
  use pprocerror
  use genoutput
  use cubature
  use dg2d_routines
  
    
  use poisson2d_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine dg2d_0_simple
  
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
    ! We need a couple of variables for this problem. Let us see...
    !
    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! Path to the mesh
    character(len=SYS_STRLEN) :: spredir

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinformconv, rlinformedge
    
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar) :: rmatrixMC, rmatrixCX, rmatrixCY
    type(t_vectorScalar) :: rrhs,rsol,redge,rconv,rsoltemp,rrhstemp,rsolUp

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlock
    type(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock,rsolBlock,redgeBlock,rconvBlock,rsolTempBlock,rsolUpBlock

    ! A set of variables describing the discrete boundary conditions.    
    type(t_boundaryRegion) :: rboundaryRegion
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
    
    real(DP) :: ttime, dt, ttfinal
    
    real(DP), dimension(2) :: vel
    
    integer :: ielementType
    
    dt = 0.01_DP
    ttfinal=0.01_DP
    
    ielementType = EL_DG_T2_2D
    
    
    
        
    vel(1)=1.0_DP
    vel(2)=1.0_DP

    ! Ok, let us start. 
    !
    ! We want to solve our Poisson problem on level...
    NLMAX = 2
    
    ! Get the path $PREDIR from the environment, where to read .prm/.tri files 
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './pre'

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, trim(spredir)//'/QUAD.prm')
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtriangulation, trim(spredir)//'/QUAD.tri', rboundary)
     
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                 rtriangulation, rboundary)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   ielementType,CUB_G3X3,rtriangulation, rboundary)
                 
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrixMC,&
                                     rdiscretisation%RspatialDiscr(1),&
                                     BILF_MATC_EDGEBASED)
    
    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
    ! scalar system matrix in 2D.

    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    
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
    call bilf_buildMatrixScalar (rform,.true.,rmatrixMC,coeff_Laplace_2D)
    
    ! Next we create the Matrices CX and CY 
    ! They represent the diskretisation of the nabla operator

    ! We could do this as we did with the mass matrix MC, but as CX and MC
    ! are of the same structure we can as well create a shared matrix,
    ! that shares the handles for the structure (Kld, Kcol, Kdiagonal)
    ! with the Matrix MC and only creates a new haldle for the entries (Data)
    ! So instead of using
    ! CALL bilf_createMatrixStructure (rdiscretisation%RspatialDiscretisation(1),&
    !                                    LSYSSC_MATRIX9,rmatrixCX)
    ! we use
    call lsyssc_duplicateMatrix(rmatrixMC, rmatrixCX,&
         LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

    if(ielementType .ne. EL_DG_T0_2D) call stdop_assembleSimpleMatrix(rmatrixCX, DER_DERIV_X, DER_FUNC)

    ! Now we do the same for CY
    call lsyssc_duplicateMatrix(rmatrixMC, rmatrixCY,&
         LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

    if(ielementType .ne. EL_DG_T0_2D) call stdop_assembleSimpleMatrix(rmatrixCY, DER_DERIV_Y, DER_FUNC)

    ! The same has to be done for the right hand side of the problem.
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinformedge%itermCount = 1
    rlinformedge%Idescriptors(1) = DER_FUNC2D
    
    
    
    ! Create scalar vectors
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1),rrhs ,.true.,ST_DOUBLE)
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1),redge,.true.,ST_DOUBLE)
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1),rconv,.true.,ST_DOUBLE)
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1),rsol ,.true.,ST_DOUBLE)
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1),rsoltemp ,.true.,ST_DOUBLE)
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1),rrhstemp ,.true.,ST_DOUBLE)
    call lsyssc_createVecByDiscr (rdiscretisation%RspatialDiscr(1),rsolUp ,.true.,ST_DOUBLE)

                                 
    ! Test the new DG edgebased routine                                 
    call linf_dg_buildVectorScalarEdge2d (rlinformedge, CUB_G3_1D, .true.,&
                                              redge,rsol,&
                                              flux_dg_buildVectorScEdge2D_sim)
                                 
    
    ! The linear solver only works for block matrices/vectors - but above,
    ! we created scalar ones. So the next step is to make a 1x1 block
    ! system from the matrices/vectors above which the linear solver
    ! understands.
    call lsysbl_createMatFromScalar (rmatrixMC,rmatrixBlock,rdiscretisation)
    call lsysbl_createVecFromScalar (rrhs,rrhsBlock,rdiscretisation)
    call lsysbl_createVecFromScalar (rsol,rsolBlock,rdiscretisation)
    call lsysbl_createVecFromScalar (rsolUp,rsolUpBlock,rdiscretisation)
    
!    ! Now we have the raw problem. What is missing is the definition of the boundary
!    ! conditions.
!    ! For implementing boundary conditions, we use a `filter technique with
!    ! discretised boundary conditions`. This means, we first have to calculate
!    ! a discrete version of the analytic BC, which we can implement into the
!    ! solution/RHS vectors using the corresponding filter.
!    !
!    ! Create a t_discreteBC structure where we store all discretised boundary
!    ! conditions.
!    call bcasm_initDiscreteBC(rdiscreteBC)
!    !
!    ! We 'know' already (from the problem definition) that we have four boundary
!    ! segments in the domain. Each of these, we want to use for enforcing
!    ! some kind of boundary condition.
!    !
!    ! We ask the bondary routines to create a 'boundary region' - which is
!    ! simply a part of the boundary corresponding to a boundary segment.
!    ! A boundary region roughly contains the type, the min/max parameter value
!    ! and whether the endpoints are inside the region or not.
!    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
!    
!    ! We use this boundary region and specify that we want to have Dirichlet
!    ! boundary there. The following call does the following:
!    ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
!    !   We specify icomponent='1' to indicate that we set up the
!    !   Dirichlet BC`s for the first (here: one and only) component in the 
!    !   solution vector.
!    ! - Discretise the boundary condition so that the BC`s can be applied
!    !   to matrices and vectors
!    ! - Add the calculated discrete BC`s to rdiscreteBC for later use.
!    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)
!                             
!    ! Now to the edge 2 of boundary component 1 the domain.
!    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
!    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)
!                             
!    ! Edge 3 of boundary component 1.
!    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
!    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)
!    
!    ! Edge 4 of boundary component 1. That is it.
!    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
!    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
!                                       rboundaryRegion,rdiscreteBC,&
!                                       getBoundaryValues_2D)
!                             
!    ! Hang the pointer into the vector and matrix. That way, these
!    ! boundary conditions are always connected to that matrix and that
!    ! vector.
!    rmatrixBlock%p_rdiscreteBC => rdiscreteBC
!    rrhsBlock%p_rdiscreteBC => rdiscreteBC
!                             
    ! Now we have block vectors for the RHS and the matrix. What we
    ! need additionally is a block vector for the solution and
    ! temporary data. Create them using the RHS as template.
    ! Fill the solution vector with 0:
    !call lsysbl_createVecBlockIndirect (rrhsBlock, rsolBlock, .true.)
    !call lsysbl_createVecBlockIndirect (rrhsBlock, rsolTempBlock, .true.)
    call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .false.)
!    
!    ! Next step is to implement boundary conditions into the RHS,
!    ! solution and matrix. This is done using a vector/matrix filter
!    ! for discrete boundary conditions.
!    ! The discrete boundary conditions are already attached to the
!    ! vectors/matrix. Call the appropriate vector/matrix filter that
!    ! modifies the vectors/matrix according to the boundary conditions.
!    call vecfil_discreteBCrhs (rrhsBlock)
!    call vecfil_discreteBCsol (rvectorBlock)
!    call matfil_discreteBC (rmatrixBlock)
!
!    ! During the linear solver, the boundary conditions are also
!    ! frequently imposed to the vectors. But as the linear solver
!    ! does not work with the actual solution vectors but with
!    ! defect vectors instead.
!    ! So, set up a filter chain that filters the defect vector
!    ! during the solution process to implement discrete boundary conditions.
!    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
!
    ! Create a BiCGStab-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    nullify(p_rpreconditioner)
    call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,RfilterChain)
    
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
    
    
    
    
    
    ttime = 0.0_DP

    timestepping: do

       ! Compute solution from time step t^n to time step t^{n+1}
       write(*,*)
       write(*,*)
       write(*,*) 'TIME STEP:', ttime
       write(*,*)
       
       
       
       
       
       call lsyssc_copyVector(rsol,rsoltemp)
       
       
       ! Step 1/3
       
       ! Create RHS-Vector
             
       ! First use the dg-function for the edge terms
       call linf_dg_buildVectorScalarEdge2d (rlinformedge, CUB_G3_1D, .true.,&
                                              rrhs,rsolTemp,&
                                              flux_dg_buildVectorScEdge2D_sim)
                                              
       
       call lsyssc_scaleVector (rrhs,-1.0_DP)
       ! Then add the convection terms
       if(ielementType .ne. EL_DG_T0_2D) call lsyssc_scalarMatVec (rmatrixCX, rsolTemp, rrhs, vel(1), 1.0_DP)
       if(ielementType .ne. EL_DG_T0_2D) call lsyssc_scalarMatVec (rmatrixCY, rsolTemp, rrhs, vel(2), 1.0_DP)
       
       ! Solve for solution update
       call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
       
       ! Get new temp solution
       call lsyssc_vectorLinearComb (rsol,rsolUp,1.0_DP,dt,rsoltemp)
       
       
       ! Step 2/3
       
       ! Create RHS-Vector
             
       ! First use the dg-function for the edge terms
       call linf_dg_buildVectorScalarEdge2d (rlinformedge, CUB_G3_1D, .true.,&
                                              rrhs,rsolTemp,&
                                              flux_dg_buildVectorScEdge2D_sim)
       call lsyssc_scaleVector (rrhs,-1.0_DP)
       ! Then add the convection terms
       if(ielementType .ne. EL_DG_T0_2D) call lsyssc_scalarMatVec (rmatrixCX, rsolTemp, rrhs, vel(1), 1.0_DP)
       if(ielementType .ne. EL_DG_T0_2D) call lsyssc_scalarMatVec (rmatrixCY, rsolTemp, rrhs, vel(2), 1.0_DP)
       
       ! Solve for solution update
       call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
       
       ! Get new temp solution
       call lsyssc_vectorLinearComb (rsoltemp,rsolUp,1.0_DP,dt)
       call lsyssc_vectorLinearComb (rsol,rsolUp,0.75_DP,0.25_DP,rsoltemp)


       ! Step 3/3
       
       ! Create RHS-Vector
             
       ! First use the dg-function for the edge terms
       call linf_dg_buildVectorScalarEdge2d (rlinformedge, CUB_G3_1D, .true.,&
                                              rrhs,rsolTemp,&
                                              flux_dg_buildVectorScEdge2D_sim)
       call lsyssc_scaleVector (rrhs,-1.0_DP)
       ! Then add the convection terms
       if(ielementType .ne. EL_DG_T0_2D) call lsyssc_scalarMatVec (rmatrixCX, rsolTemp, rrhs, vel(1), 1.0_DP)
       if(ielementType .ne. EL_DG_T0_2D) call lsyssc_scalarMatVec (rmatrixCY, rsolTemp, rrhs, vel(2), 1.0_DP)
       
       ! Solve for solution update
       call linsol_solveAdaptively (p_rsolverNode,rsolUpBlock,rrhsBlock,rtempBlock)
       
       ! Get new temp solution
       call lsyssc_vectorLinearComb (rsoltemp,rsolUp,1.0_DP,dt)
       call lsyssc_vectorLinearComb (rsolUp,rsol,2.0_DP/3.0_DP,1.0_DP/3.0_DP)       
       
       
    
       ! Go on to the next time step
       ttime = ttime + dt

       ! Leave the time stepping loop if final time is reached
       if (ttime .ge. ttfinal-0.001_DP*dt) exit timestepping

    end do timestepping
    
    ! Output solution to gmv file
    call dg2gmv(rsol,3)


call lsyssc_getbase_double (rsol,p_Ddata)
       write(*,*) p_Ddata
       pause
    
    
!    
!    ! Finally solve the system. As we want to solve Ax=b with
!    ! b being the real RHS and x being the real solution vector,
!    ! we use linsol_solveAdaptively. If b is a defect
!    ! RHS and x a defect update to be added to a solution vector,
!    ! we would have to use linsol_precondDefect instead.
!    call linsol_solveAdaptively (p_rsolverNode,rvectorBlock,rrhsBlock,rtempBlock)
!    
!    ! That is it, rvectorBlock now contains our solution. We can now
!    ! start the postprocessing. 
!    !
!    ! Get the path for writing postprocessing files from the environment variable
!    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
!    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'
!
!    ! Start UCD export to GMV file:
!    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
!                       trim(sucddir)//'/u2d_0_simple.gmv')
!    
!    call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
!    call ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
!    
!    ! Write the file to disc, that is it.
!    call ucd_write (rexport)
!    call ucd_release (rexport)
!    
!    ! Calculate the error to the reference function.
!    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,derror,&
!                       getReferenceFunction_2D)
!    call output_line ('L2-error: ' // sys_sdEL(derror,10) )
!
!    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
!                       getReferenceFunction_2D)
!    call output_line ('H1-error: ' // sys_sdEL(derror,10) )
!    
!    ! We are finished - but not completely!
!    ! Now, clean up so that all the memory is available again.
!    !
!    ! Release solver data and structure
!    call linsol_doneData (p_rsolverNode)
!    call linsol_doneStructure (p_rsolverNode)
!    
!    ! Release the solver node and all subnodes attached to it (if at all):
!    call linsol_releaseSolver (p_rsolverNode)
!    
!    ! Release the block matrix/vectors
!    call lsysbl_releaseVector (rtempBlock)
!    call lsysbl_releaseVector (rvectorBlock)
!    call lsysbl_releaseVector (rrhsBlock)
!    call lsysbl_releaseMatrix (rmatrixBlock)
!
!    ! Release the scalar matrix/rhs vector which were used to create
!    ! the block matrices/vectors. These must exist as long as the
!    ! block matrices/vectors exist, as the block matrices/vectors are
!    ! only 'copies' of the scalar ones, sharing the same handles!
!    call lsyssc_releaseVector (rrhs)
!    call lsyssc_releaseMatrix (rmatrix)
!    
!    ! Release our discrete version of the boundary conditions
!    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that is it.
    call boundary_release (rboundary)
    
  end subroutine

end module

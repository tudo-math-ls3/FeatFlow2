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

module poisson2d_method0_simple
  
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
    
  use poisson2d_callback
  USE statistics

  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine poisson2d_0_simple
  
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
    type(t_linearForm) :: rlinform
    
    ! A scalar matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form.
    type(t_matrixScalar) :: rmatrix, rmatrixLaplace, rmatrixMass
    type(t_vectorScalar) :: rrhs

    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock) :: rmatrixBlock
    type(t_vectorBlock) :: rrhsBlock,rvectorBlock,rtempBlock

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
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror    
    
    ! Error of FE function to reference function
    real(DP) :: derror, Alpha
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir
    real(DP), dimension(:), pointer :: p_Ddata

    	TYPE (t_timer) :: rtimer
       
!*****************************************************************************
! Clear and start the computation time of the subroutine.
!----------------------------------------------------------------------------   
     	CALL stat_clearTimer(rtimer)
     	CALL stat_startTimer(rtimer)

    ! Ok, let's start. 
    !
    ! We want to solve our Poisson problem on level...
    NLMAX = 5
    Alpha=0.0000000001_DP
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
    call spdiscr_initBlockDiscr (rdiscretisation,2,&
                                   rtriangulation, rboundary)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
!     call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
!                                    EL_E011,CUB_G2X2,rtriangulation, rboundary)
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   EL_Q2,CUB_G3X3,rtriangulation, rboundary)

    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(2), &
                                   EL_Q2,CUB_G3X3,rtriangulation, rboundary)
                 
    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rmatrix)
    call  lsyssc_duplicateMatrix(rmatrix,rmatrixLaplace,LSYSSC_DUP_SHARE,&
					LSYSSC_DUP_IGNORE)
    call  lsyssc_duplicateMatrix(rmatrix,rmatrixMass,LSYSSC_DUP_SHARE,&
					LSYSSC_DUP_IGNORE)

  ! Initialise the block matrix with default values based on
    ! the discretisation.
    call lsysbl_createMatBlockByDiscr (rdiscretisation,rmatrixBlock)    
   
      ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create that directly in the block (1,1) of the block matrix
    ! using the discretisation structure of the first block.
    !
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9, rmatrixBlock%RmatrixBlock(1,1))
    call lsyssc_duplicateMatrix (rmatrixBlock%RmatrixBlock(1,1),&
             rmatrixBlock%RmatrixBlock(1,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rmatrixBlock%RmatrixBlock(1,1),&
             rmatrixBlock%RmatrixBlock(2,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    call lsyssc_duplicateMatrix (rmatrixBlock%RmatrixBlock(1,1),&
             rmatrixBlock%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!----------------------------------------------------------------------
! Set up matrixLaplace

    ! And now to the entries of the  Laplace matrix. For assembling of the entries,
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
    call bilf_buildMatrixScalar (rform,.true.,rmatrixLaplace,coeff_Laplace_2D)

   CALL lsyssc_copyMatrix(rmatrixLaplace,rmatrixBlock%RmatrixBlock(1,1))
   CALL lsyssc_copyMatrix(rmatrixLaplace,rmatrixBlock%RmatrixBlock(2,2))
   !---------------------------------------------------------------------- 
! Set up matrixMass

      ! And now to the entries of the mass matrix. 
      ! For assembling of the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (Psi_j, Phi_i) for the
      ! scalar system matrix in 2D.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_FUNC
      
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = 1.0 
      call bilf_buildMatrixScalar (rform,.true.,rmatrixMass) 

   CALL lsyssc_copyMatrix(rmatrixMass,rmatrixBlock%RmatrixBlock(2,1))
   rmatrixBlock%RmatrixBlock(2,1)%dScaleFactor=-1.0_DP
   CALL lsyssc_copyMatrix(rmatrixMass,rmatrixBlock%RmatrixBlock(1,2))
   rmatrixBlock%RmatrixBlock(1,2)%dScaleFactor=1.0_DP/Alpha
!----------------------------------------------------------------------
    ! Update the structural information of the block matrix, as we manually
    ! changed the submatrices:
    call lsysbl_updateMatStrucInfo (rmatrixBlock)

    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    call lsysbl_createVecBlockIndMat (rmatrixBlock,rrhsBlock, .false.)
    call lsysbl_createVecBlockIndMat (rmatrixBlock,rvectorBlock, .false.)
    call lsysbl_createVecBlockIndMat (rmatrixBlock,rtempBlock, .false.)
! ---------------------------------------------------------------------- 
!  Set up RHS_1

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
                                 rlinform,.true.,rrhsBlock%RvectorBlock(1),coeff_RHS_2D1)
! ---------------------------------------------------------------------- 
!  Set up RHS_2

    ! The same has to be done for the right hand side of the problem.
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to get a discrete version of it.
    ! Again we simply create a scalar vector based on the one and only
    ! discretisation structure.
    ! This scalar vector will later be used as the one and only first
    ! component in a block vector.

     call linf_buildVectorScalar (rdiscretisation%RspatialDiscr(2),&
                                 rlinform,.true.,rrhsBlock%RvectorBlock(2),coeff_RHS_2D2)
     CALL lsyssc_scaleVector(rrhsBlock%RvectorBlock(2),-1.0_DP)
!  ---------------------------------------------------------------------- 
    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rvectorBlock)
    call lsysbl_clearVector(rtempBlock)
!  ---------------------------------------------------------------------- 
!     Create BC-Y

    ! Now we have the raw problem. What is missing is the definition of the boudary
    ! conditions.
    ! For implementing boundary conditions, we use a 'filter technique with
    ! discretised boundary conditions'. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! Create a t_discreteBC structure where we store all discretised boundary
    ! conditions.
    call bcasm_initDiscreteBC(rdiscreteBC)
    !
    ! We 'know' already (from the problem definition) that we have four boundary
    ! segments in the domain. Each of these, we want to use for enforcing
    ! some kind of boundary condition.
    !
    ! We ask the bondary routines to create a 'boundary region' - which is
    ! simply a part of the boundary corresponding to a boundary segment.
    ! A boundary region roughly contains the type, the min/max parameter value
    ! and whether the endpoints are inside the region or not.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    
    ! We use this boundary region and specify that we want to have Dirichlet
    ! boundary there. The following call does the following:
    ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
    !   We specify icomponent='1' to indicate that we set up the
    !   Dirichlet BC's for the first (here: one and only) component in the 
    !   solution vector.
    ! - Discretise the boundary condition so that the BC's can be applied
    !   to matrices and vectors
    ! - Add the calculated discrete BC's to rdiscreteBC for later use.
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)

    ! Now to the edge 2 of boundary component 1 the domain.
    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)

    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That's it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D)
!  ---------------------------------------------------------------------   
!     Create BC-P

    ! We 'know' already (from the problem definition) that we have four boundary
    ! segments in the domain. Each of these, we want to use for enforcing
    ! some kind of boundary condition.
    ! We ask the bondary routines to create a 'boundary region' - which is
    ! simply a part of the boundary corresponding to a boundary segment.
    ! A boundary region roughly contains the type, the min/max parameter value
    ! and whether the endpoints are inside the region or not.
    call boundary_createRegion(rboundary,1,1,rboundaryRegion)
    
    ! We use this boundary region and specify that we want to have Dirichlet
    ! boundary there. The following call does the following:
    ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
    !   We specify icomponent='1' to indicate that we set up the
    !   Dirichlet BC's for the first (here: one and only) component in the 
    !   solution vector.
    ! - Discretise the boundary condition so that the BC's can be applied
    !   to matrices and vectors
    ! - Add the calculated discrete BC's to rdiscreteBC for later use.
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D_P)
                             
    ! Now to the edge 2 of boundary component 1 the domain.
    call boundary_createRegion(rboundary,1,2,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D_P)
                             
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D_P)
    
    ! Edge 4 of boundary component 1. That's it.
    call boundary_createRegion(rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                       rboundaryRegion,rdiscreteBC,&
                                       getBoundaryValues_2D_P)
! -----------------------------------------------------------------                   
                             
    ! Hang the pointer into the vector and matrix. That way, these
    ! boundary conditions are always connected to that matrix and that
    ! vector.
    rmatrixBlock%p_rdiscreteBC => rdiscreteBC
    rrhsBlock%p_rdiscreteBC => rdiscreteBC
                             
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
    p_RfilterChain => RfilterChain
    nullify(p_rpreconditioner)
!     call linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
    call linsol_initUMFPACK4 (p_rsolverNode)

    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2
    p_rsolverNode%nmaxIterations=1000
    p_rsolverNode%depsRel = 1E-5_DP
    ! Attach the system matrix to the solver.
    ! First create an array with the matrix data (on all levels, but we
    ! only have one level here), then call the initialisation 
    ! routine to attach all these matrices.
    ! Remark: Don't make a call like
    !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
    ! This doesn't work on all compilers, since the compiler would have
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
    
    ! That's it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing. 
    !
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       TRIM(sucddir)//'/u2d_0_simple.gmv')
    
    call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol_y',UCD_VAR_STANDARD, p_Ddata)
  
    call lsyssc_getbase_double (rvectorBlock%RvectorBlock(2),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol_p',UCD_VAR_STANDARD, p_Ddata)

    ! Calculate the error to the reference function.
    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,derror,&
                       getReferenceFunction_2D)
    call output_line ('L2-error-Y: ' // sys_sdEL(derror,10) )
    call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_H1ERROR,derror,&
                       getReferenceFunction_2D)
    call output_line ('H1-error-Y: ' // sys_sdEL(derror,10) )

!     call pperr_scalar (rvectorBlock%RvectorBlock(2),PPERR_L2ERROR,derror,&
!                        getReferenceFunction_2D)
!     call output_line ('L2-error-P: ' // sys_sdEL(derror,10) )
!         call pperr_scalar (rvectorBlock%RvectorBlock(2),PPERR_H1ERROR,derror,&
!                        getReferenceFunction_2D)
!     call output_line ('H1-error-P: ' // sys_sdEL(derror,10) )   
    
    CALL lsyssc_scaleVector (rvectorBlock%RvectorBlock(2),-1.0_DP/Alpha)
    call lsyssc_getbase_double (rvectorBlock%RvectorBlock(2),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'sol_u',UCD_VAR_STANDARD, p_Ddata)

    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)   
    

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
!     call lsyssc_releaseVector (rrhs)
    call lsyssc_releaseMatrix (rmatrix)
    
    call lsyssc_releaseMatrix (rmatrixLaplace)
    call lsyssc_releaseMatrix (rmatrixMass)
    
    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rdiscreteBC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    call tria_done (rtriangulation)
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)
! *********************************************************************************
! Stop and print the computation time.
!---------------------------------------------------------------------------------      
     	CALL stat_stopTimer(rtimer)
!--------------------------------------------------------------------------------- 
       WRITE(*,*)"Time for Computation :=",rtimer%delapsedReal

    
  end subroutine

end module

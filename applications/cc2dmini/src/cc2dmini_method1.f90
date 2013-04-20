!##############################################################################
!# ****************************************************************************
!# <name> cc2dmini_method1 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a stationary
!# Navier-Stokes problem
!#
!#              $$- \nu Laplace(u) + u*grad(u) + \Nabla p = f $$
!#              $$ \Nable \cdot p = 0$$
!#
!# on a 2D domain for a 2D function $u=(u_1,u_2)$ and a pressure $p$.
!#
!# The routine splits up the tasks of reading the domain, creating
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines.
!#
!# For the nonlinearity, the nonlinear solver is invoked. The
!# defect that is setted up there is preconditioned by a linear Multigrid
!# solver with a simple-VANKA smoother/preconditioner for
!# 2D saddle point problems, Jacobi-Type. As coarse grid solver,
!# UMFPACK is used.
!# </purpose>
!##############################################################################

module cc2dmini_method1

  use fsystem
  use storage
  use linearsolver
  use boundary
  use cubature
  use scalarpde
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use ucd
  use boundarycondition
  use discretebc
  use linearsystemscalar
  use linearsystemblock
  use scalarpde
  use derivatives
  use element
  use bilinearformevaluation
  use linearformevaluation
  use filtersupport
  use multilevelprojection
  
  use collection
  use convection
    
  use cc2dmini_callback
  
  implicit none
  
  ! Maximum allowed level in this application; must be =9 for
  ! FEAT 1.x compatibility (still)!
  integer, parameter :: NNLEV = 9
  
!<types>

!<typeblock description="Type block defining all information about one level">

  type t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...)
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! A system matrix for that specific level.
    type(t_matrixBlock) :: rmatrix

    ! Laplace matrix for that specific level.
    type(t_matrixScalar) :: rmatrixLaplace

    ! B1-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB2

    ! A temporary vector for building the solution when assembling the
    ! matrix on lower levels.
    type(t_vectorBlock) :: rtempVector

    ! A variable describing the discrete boundary conditions fo the velocity
    type(t_discreteBC) :: rdiscreteBC
  
  end type
  
!</typeblock>


!<typeblock description="Application-specific type block for Nav.St. problem">

  type t_problem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    integer :: NLMIN
    
    ! Maximum refinement level
    integer :: NLMAX
    
    ! Viscosity parameter nu = 1/Re
    real(DP) :: dnu

    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! A solution vector and a RHS vector on the finest level.
    type(t_vectorBlock) :: rvector,rrhs

    ! A variable describing the analytic boundary conditions.
    type(t_boundaryConditions), pointer :: p_rboundaryConditions

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by NLMAX!
    type(t_problem_lvl), dimension(NNLEV) :: RlevelInfo
    
    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    type(t_collection) :: rcollection

  end type

!</typeblock>

!</types>
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine c2d1_initParamTriang (ilvmin,ilvmax,rproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
!</description>

!<input>
  ! Minimum refinement level of the mesh; = coarse grid = level 1
  integer, intent(in) :: ilvmin
  
  ! Maximum refinement level
  integer, intent(in) :: ilvmax
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  
    ! Initialise the level in the problem structure
    rproblem%NLMIN = ilvmin
    rproblem%NLMAX = ilvmax

    call collct_setvalue_int(rproblem%rcollection,'NLMIN',ilvmin,.true.)
    call collct_setvalue_int(rproblem%rcollection,'NLMAX',ilvmax,.true.)

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rproblem%rboundary, './pre/bench1.prm')
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation, &
                                                 './pre/bench1.tri', &
        rproblem%rboundary)

    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering(rproblem%NLMIN-1,&
        rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,rproblem%rboundary)

    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (&
        rproblem%RlevelInfo(rproblem%NLMIN)%rtriangulation,rproblem%rboundary)
    
    ! Now, refine to level up to nlmax.
    do i=rproblem%NLMIN+1,rproblem%NLMAX
      call tria_refine2LevelOrdering (rproblem%RlevelInfo(i-1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)
      call tria_initStandardMeshFromRaw (rproblem%RlevelInfo(i)%rtriangulation,&
          rproblem%rboundary)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d1_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: I
  
    ! An object for saving the domain:
    type(t_boundary), pointer :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! An object for the block discretisation on one level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    do i=rproblem%NLMIN,rproblem%NLMAX
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rproblem%rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies 3 blocks in the
      ! solution vector. In this simple problem, we only have one block.
      allocate(p_rdiscretisation)
      call spdiscr_initBlockDiscr (p_rdiscretisation,3,&
                                   p_rtriangulation, p_rboundary)

      ! Save the discretisation structure to our local LevelInfo structure
      ! for later use.
      rproblem%RlevelInfo(i)%p_rdiscretisation => p_rdiscretisation

      ! p_rdiscretisation%Rdiscretisations is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  p_rdiscretisation%RspatialDiscr(1), &
                  EL_EM30,CUB_G2X2, &
                  p_rtriangulation, p_rboundary)
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards.
      call spdiscr_duplicateDiscrSc(p_rdiscretisation%RspatialDiscr(1),&
          p_rdiscretisation%RspatialDiscr(2))
  
      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (p_rdiscretisation%RspatialDiscr(1),&
          EL_Q0,CUB_G2X2,p_rdiscretisation%RspatialDiscr(3))
          
    end do
                                   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d1_initMatVec (rproblem)
  
!<description>
  ! Calculates the system matrix and RHS vector of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_matrixScalar), pointer :: p_rmatrixLaplace
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector,p_rtempVector

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
    do i=rproblem%NLMIN,rproblem%NLMAX
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      ! The global system looks as follows:
      !
      !    ( A         B1 )
      !    (      A    B2 )
      !    ( B1^T B2^T    )
      !
      ! with A = L + nonlinear Convection. We compute in advance
      ! a standard Laplace matrix L which can be added later to the
      ! convection matrix, resulting in the nonlinear system matrix.
      !
      ! Get a pointer to the (scalar) Laplace matrix:
      p_rmatrixLaplace => rproblem%RlevelInfo(i)%rmatrixLaplace
      
      ! and save it to the collection for later use.
      call collct_setvalue_matsca(rproblem%rcollection,'LAPLACE',&
                                  p_rmatrixLaplace,.true.,i)
      
      ! Create the matrix structure of the Laplace matrix:
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrixLaplace)
      
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
      rform%Dcoefficients(1)  = rproblem%dnu
      rform%Dcoefficients(2)  = rproblem%dnu

      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_Stokes for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine,
      ! so the callback routine has access to everything what is
      ! in the collection.
      call bilf_buildMatrixScalar (rform,.true.,&
                                   p_rmatrixLaplace,coeff_Stokes,&
                                   rproblem%rcollection)
      
      ! In the global system, there are two coupling matrices B1 and B2.
      ! Both have the same structure.
      ! Create the matrices structure of the pressure using the 3rd
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(3),LSYSSC_MATRIX9,&
                rproblem%RlevelInfo(i)%rmatrixB1,&
                p_rdiscretisation%RspatialDiscr(1))
                
      ! Duplicate the B1 matrix structure to the B2 matrix, so use
      ! lsyssc_duplicateMatrix to create B2. Share the matrix
      ! structure between B1 and B2 (B1 is the parent and B2 the child).
      ! Don't create a content array yet, it will be created by
      ! the assembly routines later.
      call lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB1,&
                  rproblem%RlevelInfo(i)%rmatrixB2,LSYSSC_DUP_COPY,LSYSSC_DUP_REMOVE)

      ! Build the first pressure matrix B1.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_X

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = -1.0_DP
      
      call bilf_buildMatrixScalar (rform,.true.,&
                                  rproblem%RlevelInfo(i)%rmatrixB1,coeff_Pressure,&
                                  rproblem%rcollection)

      ! Build the second pressure matrix B2.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = -1.0_DP
      
      call bilf_buildMatrixScalar (rform,.true.,&
                                  rproblem%RlevelInfo(i)%rmatrixB2,coeff_Pressure,&
                                  rproblem%rcollection)
                                  
      ! Now let's come to the main system matrix, which is a block matrix.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)
      
      ! Save the system matrix to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rproblem%rcollection,'SYSTEMMAT',p_rmatrix,.true.,i)

      ! Inform the matrix that we build a saddle-point problem.
      ! Normally, imatrixSpec has the value LSYSBS_MSPEC_GENERAL,
      ! but probably some solvers can use the special structure later.
      p_rmatrix%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
      
      ! Let's consider the global system in detail:
      !
      !    ( A         B1 ) = ( A11  A12  A13 )
      !    (      A    B2 )   ( A21  A22  A23 )
      !    ( B1^T B2^T    )   ( A31  A32  A33 )
      !
      ! The matrices A11 and A22 of the global system matrix have exactly
      ! the same structure as the original Laplace matrix from above!
      ! Initialise them with the same structure, i.e. A11, A22 and the
      ! Laplace matrix L share(!) the same structure.
      !
      ! For this purpose, use the "duplicate matric" routine.
      ! The structure of the matrix is shared with the Laplace matrix.
      ! For the content, a new empty array is allocated which will later receive
      ! the entries.
      call lsyssc_duplicateMatrix (p_rmatrixLaplace,&
                  p_rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
                  
      ! The matrix A22 is identical to A11! So mirror A11 to A22 sharing the
      ! structure and the content.
      call lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
                  p_rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      ! Manually change the discretisation structure of the Y-velocity
      ! matrix to the Y-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      call lsyssc_assignDiscretisation (p_rmatrix%RmatrixBlock(2,2),&
          p_rdiscretisation%RspatialDiscr(2))
                                  
      ! The B1/B2 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create copies of the entries. The reason is
      ! that these matrices are modified for boundary conditions later.
      call lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
                                   p_rmatrix%RmatrixBlock(1,3),&
                                   LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      call lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
                                   p_rmatrix%RmatrixBlock(2,3),&
                                   LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
      ! Furthermore, put B1^T and B2^T to the block matrix.
      call lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
                                   p_rmatrix%RmatrixBlock(3,1),&
                                   LSYSSC_TR_VIRTUAL)

      call lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
                                   p_rmatrix%RmatrixBlock(3,2),&
                                   LSYSSC_TR_VIRTUAL)

      ! Now on all levels except for the maximum one, create a temporary
      ! vector on that level, based on the matrix template.
      ! It's used for building the matrices on lower levels.
      if (i .lt. rproblem%NLMAX) then
        p_rtempVector => rproblem%RlevelInfo(i)%rtempVector
        call lsysbl_createVecBlockIndMat (p_rmatrix,p_rtempVector,.false.)
        
        ! Add the temp vector to the collection on level i
        ! for use in the callback routine
        call collct_setvalue_vec(rproblem%rcollection,'RTEMPVEC',p_rtempVector,&
                                .true.,i)
      end if

    end do

    ! (Only) on the finest level, we need to calculate a RHS vector
    ! and to allocate a solution vector.
    
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector

    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    call lsysbl_createVecBlockIndMat (p_rmatrix,p_rrhs, .false.)
    call lsysbl_createVecBlockIndMat (p_rmatrix,p_rvector, .false.)

    ! Save the solution/RHS vector to the collection. Might be used
    ! later (e.g. in nonlinear problems)
    call collct_setvalue_vec(rproblem%rcollection,'RHS',p_rrhs,.true.)
    call collct_setvalue_vec(rproblem%rcollection,'SOLUTION',p_rvector,.true.)
    
    ! The vector structure is ready but the entries are missing.
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the
    ! first block.
    !
    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is
    ! in the collection.
    !
    ! Note that the vector is unsorted after calling this routine!
    !
    ! Discretise the X-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
              p_rrhs%RvectorBlock(1),coeff_RHS_x,&
              rproblem%rcollection)

    ! And the Y-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(2),rlinform,.true.,&
              p_rrhs%RvectorBlock(2),coeff_RHS_y,&
              rproblem%rcollection)

    ! The third subvector must be zero - as it represents the RHS of
    ! the equation "div(u) = 0".
    call lsyssc_clearVector(p_rrhs%RvectorBlock(3))
                                
    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rproblem%rvector)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d1_discretiseBC (rdiscretisation,rdiscreteBC)
  
!<description>
  ! This routine discretises the current boundary conditions for
  ! the discretisation specified by rdiscretisation.
!</description>

!<input>
  ! A discretisation structure specifying the current discretisation.
  type(t_blockDiscretisation), intent(in) :: rdiscretisation
!</input>

!<inputoutput>
  ! A structuree that receives the discretised boundary conditions.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_boundary), pointer :: p_rboundary
    
    p_rboundary => rdiscretisation%p_rboundary

    ! We first set up the boundary conditions for the X-velocity, then those
    ! of the Y-velocity.
    !
    ! We 'know' already (from the problem definition) that we have four boundary
    ! segments in the domain. Each of these, we want to use for enforcing
    ! some kind of boundary condition.
    !
    ! We ask the boundary routines to create a 'boundary region' - which is
    ! simply a part of the boundary corresponding to a boundary segment.
    ! A boundary region roughly contains the type, the min/max parameter value
    ! and whether the endpoints are inside the region or not.
    call boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    
    ! The endpoint of this segment should also be Dirichlet. We set this by
    ! changing the region properties in rboundaryRegion.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
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
                                      getBoundaryValues)
                              
    ! Edge 2 is Neumann boundary, so it's commented out.
    ! CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
    !                                    rboundaryRegion,rdiscreteBC,&
    !                                    getBoundaryValues)
                              
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)
    
    ! Edge 4 of boundary component 1. That's it.
    call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)

    ! The whole 2nd boundary component - if it exists.
    if (boundary_igetNBoundComp(p_rboundary) .ge. 2) then
      call boundary_createRegion(p_rboundary,2,0,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues)
    end if

    ! Now continue with defining the boundary conditions of the Y-velocity:
    !
    ! Define edge 1.
    call boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    
    ! Edge with start- and endpoint.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)
                              
    ! Edge 2 is Neumann boundary, so it's commented out.
    ! CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rdiscretisation,2,&
    !                                    rboundaryRegion,rdiscreteBC,&
    !                                    getBoundaryValues)
                              
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)
    
    ! Edge 4 of boundary component 1. That's it.
    call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)

    ! The whole 2nd boundary component - if it exists.
    if (boundary_igetNBoundComp(p_rboundary) .ge. 2) then
      call boundary_createRegion(p_rboundary,2,0,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d1_initDiscreteBC (rproblem)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

  ! A pointer to the system matrix and the RHS vector as well as
  ! the discretisation
  type(t_matrixBlock), pointer :: p_rmatrix
  type(t_vectorBlock), pointer :: p_rrhs,p_rvector
  type(t_blockDiscretisation), pointer :: p_rdiscretisation

  ! Pointer to structure for saving discrete BC's:
  type(t_discreteBC), pointer :: p_rdiscreteBC
    
    do i=rproblem%NLMIN,rproblem%NLMAX
    
      ! Get our velocity matrix from the problem structure.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => p_rmatrix%p_rblockDiscrTrial
      
      ! For implementing boundary conditions, we use a 'filter technique with
      ! discretised boundary conditions'. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
      call bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%rdiscreteBC)
      
      ! Discretise the boundary conditions.
      call c2d1_discretiseBC (p_rdiscretisation,rproblem%RlevelInfo(i)%rdiscreteBC)
                               
      ! Hang the pointer into the the matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%rdiscreteBC
      
      p_rmatrix%p_rdiscreteBC => p_rdiscreteBC
      
      ! Also hang in the boundary conditions into the temporary vector that is
      ! used for the creation of solutions on lower levels.
      ! This allows us to filter this vector when we create it.
      rproblem%RlevelInfo(i)%rtempVector%p_rdiscreteBC => p_rdiscreteBC
      
    end do

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%NLMAX)%rdiscreteBC
    
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    
    p_rrhs%p_rdiscreteBC => p_rdiscreteBC
    p_rvector%p_rdiscreteBC => p_rdiscreteBC
                
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d1_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,ilvmax
  
    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector
    
    ! Get our the right hand side and solution from the problem structure
    ! on the finest level
    ilvmax = rproblem%NLMAX
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    
    ! Filter the solution and RHS vectors through the boundary-condition-
    ! implementation filter. This implements all discrete boundary
    ! conditions.
    ! Use the RHS-filter for the RHS and the solution-filter for the
    ! solution vector.
    call vecfil_discreteBCrhs (p_rrhs)
    call vecfil_discreteBCsol (p_rvector)

    ! Implement discrete boundary conditions into the matrices on all
    ! levels, too.
    ! In fact, this modifies the B-matrices. The A-matrices are overwritten
    ! later and must then be modified again!
    do i=rproblem%NLMIN ,rproblem%NLMAX
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      call matfil_discreteBC (p_rmatrix)
    end do

  end subroutine

  ! ***************************************************************************
  !<subroutine>
  
    subroutine c2d1_getDefect (ite,rx,rb,rd,p_rcollection)
  
    use fsystem
    use linearsystemblock
    use collection
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration
    ! vector rx and the right hand side vector rb, this routine has to compute the
    ! defect vector rd. The routine accepts a pointer to a collection structure
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<input>
    ! Number of current iteration. 0=build initial defect
    integer, intent(in)                           :: ite

    ! Current iteration vector
    type(t_vectorBlock), intent(in),target        :: rx

    ! Right hand side vector of the equation.
    type(t_vectorBlock), intent(in), target       :: rb
  !</input>
               
  !<inputoutput>
    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    type(t_collection), pointer                   :: p_rcollection

    ! Defect vector b-A(x)x. This must be filled by the callback routine
    ! with data.
    type(t_vectorBlock), intent(inout), target    :: rd
  !</inputoutput>
  
  !</subroutine>

    ! local variables
    integer :: ilvmax
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_matrixScalar), pointer :: p_rmatrixLaplace
    type(t_matrixBlock) :: rmatrixLaplaceBlock
    type(t_convUpwind) :: rupwind

    ! A filter chain to pre-filter the vectors and the matrix.
    type(t_filterChain), dimension(1), target :: RfilterChain

      ! Get minimum/maximum level from the collection
      ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')
      
      ! Get the system and the Laplace matrix on the maximum level
      p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',ilvmax)
      p_rmatrixLaplace => collct_getvalue_matsca (p_rcollection,'LAPLACE',ilvmax)
      
      ! Build a temporary 3x3 block matrix rmatrixLaplace with Laplace
      ! on the main diagonal:
      !
      ! (  L    0   B1 )
      ! (  0    L   B2 )
      ! ( B1^T B2^T 0  )
      !
      call lsysbl_duplicateMatrix (p_rmatrix,rmatrixLaplaceBlock,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          
      call lsyssc_duplicateMatrix (p_rmatrixLaplace,&
          rmatrixLaplaceBlock%RmatrixBlock(1,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      call lsyssc_duplicateMatrix (p_rmatrixLaplace,&
          rmatrixLaplaceBlock%RmatrixBlock(2,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Now, in the first step, we build the linear part of the nonlinear defect:
      !     d_lin = rhs - (-nu * Laplace(.))*solution
      call lsysbl_copyVector (rb,rd)
      call lsysbl_matVec (rmatrixLaplaceBlock, rx, rd, -1.0_DP, 1.0_DP)
      
      ! Release the temporary matrix again.
      call lsysbl_releaseMatrix (rmatrixLaplaceBlock)
      
      ! For the final defect
      !
      !     d = rhs - (-nu * Laplace(.))*solution - u*grad(.)*solution
      !       = d_lin -  u*grad(.)*solution
      !
      ! we need the nonlinearity.
      
      ! Set up the upwind structure for the creation of the defect.
      ! There's not much to do, only initialise the viscosity...
      rupwind%dnu = collct_getvalue_real (p_rcollection,'NU')
      
      ! Set up a filter that modifies the block vectors/matrix
      ! according to boundary conditions.
      ! Initialise the first filter of the filter chain as boundary
      ! implementation filter for defect vectors:
      RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

      ! Apply the filter chain to the defect vector.
      ! As the filter consists only of an implementation filter for
      ! boundary conditions, this implements the boundary conditions
      ! into the defect vector.
      call filter_applyFilterChainVec (rd, RfilterChain)
      
      ! Call the upwind method to calculate the nonlinear defect.
      ! As we calculate only the defect, the matrix is ignored!
      call conv_upwind2d (rx, rx, 1.0_DP, 0.0_DP,&
                          rupwind, CONV_MODDEFECT, &
                          p_rmatrix%RmatrixBlock(1,1), rx, rd)
      
      ! Apply the filter chain to the defect vector again - since the
      ! implementation of the nonlinearity usually changes the Dirichlet
      ! nodes in the vector!
      call filter_applyFilterChainVec (rd, RfilterChain)

      ! That's it
      
    end subroutine
    
  ! ***************************************************************************

  !<subroutine>

    subroutine c2d1_getOptimalDamping (rd,rx,rb,rtemp1,rtemp2,domega,p_rcollection)
  
  !<description>
    ! This subroutine is called inside of the nonlinear loop, to be precise,
    ! inside of c2d1_precondDefect. It calculates an optiman damping parameter
    ! for the nonlinear defect correction.
    !
    ! The nonlinear loop reads:
    !
    !     $$ u_(n+1) = u_n + OMEGA * C^{-1}d_n $$
    !
    ! with $d_n$ the nonlinear defect and $C^{-1}$ a preconditioner (usually
    ! the linearised system).
    ! Based on the current solution $u_n$, the defect vector $d_n$, the RHS
    ! vector $f_n$ and the previous parameter OMEGA, a new
    ! OMEGA=domega value is calculated.
    !
    ! The nonlinear system matrix on the finest level in the collection is
    ! overwritten by $A(u_n+domega_{old}*C^{-1}d_n)$.
  !</description>

  !<input>
    ! Current iteration vector
    type(t_vectorBlock), intent(in)               :: rx

    ! Current RHS vector of the nonlinear equation
    type(t_vectorBlock), intent(in)               :: rb

    ! Defect vector b-A(x)x.
    type(t_vectorBlock), intent(in)               :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    type(t_collection), pointer                   :: p_rcollection
  !</input>

  !<inputoutput>
    ! A temporary vector in the structure of rx
    type(t_vectorBlock), intent(inout)            :: rtemp1

    ! A 2nd temporary vector in the structure of rx
    type(t_vectorBlock), intent(inout)            :: rtemp2

    ! Damping parameter. On entry: an initial value given e.g. by the
    ! previous step.
    ! On return: The new damping parameter.
    real(DP), intent(inout)                       :: domega
  !</inputoutput>
  
  !</subroutine>

    ! local variables
    integer :: ilvmax
    real(DP) :: domegaMin, domegaMax,dskv1,dskv2
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_matrixScalar), pointer :: p_rmatrixLaplace
    type(t_convUpwind) :: rupwind

    ! A filter chain to pre-filter the vectors and the matrix.
    type(t_filterChain), dimension(1), target :: RfilterChain

!    DEBUG!!!:
!    real(dp), dimension(:), pointer :: p_vec,p_def,p_temp1,p_temp2,p_da
!    call lsysbl_getbase_double (rd,p_def)
!    call lsysbl_getbase_double (rx,p_vec)
!    call lsysbl_getbase_double (rtemp1,p_temp1)
!    call lsysbl_getbase_double (rtemp2,p_temp2)
!    ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')
!    p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',ilvmax)
!    call storage_getbase_double (p_rmatrix%RmatrixBlock(1,1)%h_da,p_da)

      ! Get minimum/maximum level from the collection
      ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')
      
      ! Initialise minimum/maximum omega. In a later implementation, this
      ! might be set using a parameter from a DAT-file.
      domegaMin = 0.0_DP
      domegaMax = 2.0_DP
      
      ! Is there anything to do?
      if (domegaMin .ge. domegaMax) then
        ! No - cancel.
        domega = domegaMin
        return
      end if

      ! Get the system and the Laplace matrix on the maximum level
      p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',ilvmax)
      p_rmatrixLaplace => collct_getvalue_matsca (p_rcollection,'LAPLACE',ilvmax)

      ! Set up the upwind structure for the creation of the defect.
      ! There's not much to do, only initialise the viscosity...
      rupwind%dnu = collct_getvalue_real (p_rcollection,'NU')
      
      ! Set up a filter that modifies the block vectors/matrix
      ! according to boundary conditions.
      ! Initialise the first filter of the filter chain as boundary
      ! implementation filter for defect vectors:
      RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
      
      ! We now want to calculate a new OMEGA parameter
      ! with OMGMIN < OMEGA < OMGMAX.
      !
      ! The defect correction for a problem like T(u)u=f has the form
      !
      !       u_(n+1)  =  u_n  +  OMEGA * C * ( f - T(u_n)u_n )
      !                =  u_n  +  OMEGA * d_n
      !
      ! with an appropriate preconditioner C, which we don't care here.
      ! In our case, this iteration system can be written as:
      !
      ! (u1)     (u1)                     ( (f1)   [ A         B1] (u1) )
      ! (u2)  := (u2)  + OMEGA * C^{-1} * ( (f2) - [      A    B2] (u2) )
      ! (p )     (p )                     ( (fp)   [ B1^T B2^T 0 ] (p ) )
      !
      !                                   |------------------------------|
      !                                              = d_n
      !                            |-------------------------------------|
      !                                        = Y = (y1,y2,yp) = rd
      !
      ! with KST1=KST1(u1,u2,p) and Y=rd being the solution from
      ! the Oseen equation with
      !
      !                  [ A         B1 ]
      !    C = T(u_n) =  [      A    B2 ]
      !                  [ B1^T B2^T 0  ]
      !
      ! The parameter OMEGA is calculated as the result of the 1D
      ! minimization problem:
      !
      !   OMEGA = min_omega || T(u^l+omega*Y)*(u^l+omega*Y) - f ||_E
      !
      !           < T(u^l+omegaold*Y)Y , f - T(u^l+omegaold*Y)u^l >
      !        ~= -------------------------------------------------
      !              < T(u^l+omegaold*Y)Y , T(u^l+omegaold*Y)Y >
      !
      ! when choosing omegaold=previous omega, which is a good choice
      ! as one can see by linearization (see p. 170, Turek's book).
      !
      ! Here, ||.||_E denotes the the Euclidian norm to the Euclidian
      ! scalar product <.,.>.
      
      ! ==================================================================
      ! First term of scalar product in the nominator
      !
      ! Calculate the new nonlinear block A at the
      ! point rtemp1 = u_n + omegaold*Y
      ! ==================================================================
      !
      ! At first, calculate the point rtemp1 = u_n+omegaold*Y where
      ! to evaluate the matrix.

      call lsysbl_copyVector(rd,rtemp1)
      call lsysbl_vectorLinearComb (rx,rtemp1,1.0_DP,domega)

      ! Construct the linear part of the nonlinear matrix on the maximum
      ! level.
      !
      ! The system matrix looks like:
      !   (  A    0   B1 )
      !   (  0    A   B2 )
      !   ( B1^T B2^T 0  )
      !
      ! The A-matrix consists of Laplace+Convection.
      ! We build them separately and add together.
      !
      ! So at first, initialise the A-matrix with the Laplace contribution.
      ! We ignore the structure and simply overwrite the content of the
      ! system submatrices with the Laplace matrix.
      call lsyssc_duplicateMatrix (p_rmatrixLaplace,p_rmatrix%RmatrixBlock(1,1),&
                                   LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPY)

      ! Call the upwind method to evaluate nonlinearity part of the matrix
      ! in the point rtemp1.
      call conv_upwind2d (rtemp1, rtemp1, 1.0_DP, 0.0_DP,&
                          rupwind, CONV_MODMATRIX, &
                          p_rmatrix%RmatrixBlock(1,1))
                                    
      ! Apply the filter chain to the matrix.
      ! As the filter consists only of an implementation filter for
      ! boundary conditions, this implements the boundary conditions
      ! into the system matrix.
      call filter_applyFilterChainMat (p_rmatrix, RfilterChain)
        
      ! ==================================================================
      ! Second term of the scalar product in the nominator
      ! Calculate the defect rtemp2 = F-T*u_n.
      ! ==================================================================

      call lsysbl_copyVector (rb,rtemp2)
      call lsysbl_matVec (p_rmatrix, rx, rtemp2, -1.0_DP, 1.0_DP)
      
      ! This is a defect vector - filter it! This e.g. implements boundary
      ! conditions.
      call filter_applyFilterChainVec (rtemp2, RfilterChain)
      
      ! ==================================================================
      ! For all terms in the fraction:
      ! Calculate the value  rtemp1 = T*Y
      ! ==================================================================

      call lsysbl_matVec (p_rmatrix, rd, rtemp1, 1.0_DP, 0.0_DP)
      
      ! This is a defect vector against 0 - filter it! This e.g.
      ! implements boundary conditions.
      call filter_applyFilterChainVec (rtemp1, RfilterChain)
      
      ! ==================================================================
      ! Calculation of the fraction terms.
      ! Calculate nominator:    dskv1:= (T*Y,D)   = (rtemp1,rtemp2)
      ! Calculate denominator:  dskv2:= (T*Y,T*Y) = (rtemp1,rtemp1)
      ! ==================================================================
      
      dskv1 = lsysbl_scalarProduct (rtemp1, rtemp2)
      dskv2 = lsysbl_scalarProduct (rtemp1, rtemp1)
      
      if (dskv2 .lt. 1.0E-40_DP) then
        print *,'Error in c2d1_getOptimalDamping. dskv2 nearly zero.'
        print *,'Optimal damping parameter singular.'
        print *,'Is the triangulation ok??? .tri-file destroyed?'
        stop
      end if
      
      ! Ok, we have the nominator and the denominator. Divide them
      ! by each other to calculate the new OMEGA.
      
      domega = dskv1 / dskv2
      
      ! And make sure it's in the allowed range:
      
      domega = max(domegamin,min(domegamax,domega))
      
      ! That's it, we have our new Omega.
  
    end subroutine

  ! ***************************************************************************

  !<subroutine>

    subroutine c2d1_precondDefect (ite, rd,rx,rb,domega,bsuccess,p_rcollection)
  
    use fsystem
    use linearsystemblock
    use collection
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration
    ! vector rx and the right hand side vector rb, this routine has to compute the
    ! defect vector rd. The routine accepts a pointer to a collection structure
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<inputoutput>
    ! Number of current iteration.
    integer, intent(in)                           :: ite

    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    type(t_vectorBlock), intent(inout), target    :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    type(t_collection), pointer                   :: p_rcollection
    
    ! Damping parameter. Is set to rsolverNode%domega (usually = 1.0_DP)
    ! on the first call to the callback routine.
    ! The callback routine can modify this parameter according to any suitable
    ! algorithm to calculate an 'optimal damping' parameter. The nonlinear loop
    ! will then use this for adding rd to the solution vector:
    ! $$ x_{n+1} = x_n + domega*rd $$
    ! domega will stay at this value until it's changed again.
    real(DP), intent(inout)                       :: domega

    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    logical, intent(inout)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Current iteration vector
    type(t_vectorBlock), intent(in), target       :: rx

    ! Current right hand side of the nonlinear system
    type(t_vectorBlock), intent(in), target       :: rb
  !</input>
  
  !</subroutine>
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_matrixScalar), pointer :: p_rmatrixLaplace
    type(t_vectorScalar), pointer :: p_rvectorTemp,p_rvectorTemp2
    type(t_vectorBlock) :: rtemp1,rtemp2
    integer :: ierror,ilvmax,ilvmin, ilev
    type(t_linsolNode), pointer :: p_rsolverNode
    type(t_convUpwind) :: rupwind
    type(t_vectorBlock), pointer :: p_rvectorFine,p_rvectorCoarse

    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock), pointer :: p_rprojection

    ! A filter chain to pre-filter the vectors and the matrix.
    type(t_filterChain), dimension(1), target :: RfilterChain

!    DEBUG!!!:
!    real(dp), dimension(:), pointer :: p_vec,p_def,p_da
!    call lsysbl_getbase_double (rd,p_def)
!    call lsysbl_getbase_double (rx,p_vec)
!    ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')
!    p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',ilvmax)
!    call storage_getbase_double (p_rmatrix%RmatrixBlock(1,1)%h_da,p_da)

      ! Get minimum and maximum level from the collection
      ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')
      ilvmin = collct_getvalue_int (p_rcollection,'NLMIN')
      
      ! Set up the upwind structure for the creation of the defect.
      ! There's not much to do, only initialise the viscosity...
      rupwind%dnu = collct_getvalue_real (p_rcollection,'NU')
      
      ! Set up a filter that modifies the block vectors/matrix
      ! according to boundary conditions.
      ! Initialise the first filter of the filter chain as boundary
      ! implementation filter for defect vectors:
      RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
      
      ! Get the interlevel projection structure and the temporary vector
      ! from the collection.
      ! Our 'parent' prepared there how to interpolate the solution on the
      ! fine grid to coarser grids.
      p_rprojection => collct_getvalue_ilvp(p_rcollection,'ILVPROJECTION')
      p_rvectorTemp => collct_getvalue_vecsca(p_rcollection,'RTEMPSCALAR')

      ! On all levels, we have to set up the nonlinear system matrix,
      ! so that the linear solver can be applied to it.

      do ilev=ilvmax,ilvmin,-1
      
        ! Get the system matrix and the Laplace matrix
        p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',ilev)
        p_rmatrixLaplace => collct_getvalue_matsca (p_rcollection,'LAPLACE',ilev)
        
        ! On the highest level, we use rx as solution to build the nonlinear
        ! matrix. On lower levels, we have to create a solution
        ! on that level from a fine-grid solution before we can use
        ! it to build the matrix!
        if (ilev .eq. ilvmax) then
          p_rvectorCoarse => rx
        else
          ! Get the temporary vector on level i. Will receive the solution
          ! vector on that level.
          p_rvectorCoarse => collct_getvalue_vec (p_rcollection,'RTEMPVEC',ilev)
          
          ! Get the solution vector on level i+1. This is either the temporary
          ! vector on that level, or the solution vector on the maximum level.
          if (ilev .lt. ilvmax-1) then
            p_rvectorFine => collct_getvalue_vec (p_rcollection,'RTEMPVEC',ilev+1)
          else
            p_rvectorFine => rx
          end if

          ! Interpolate the solution from the finer grid to the coarser grid.
          ! The interpolation is configured in the interlevel projection
          ! structure we got from the collection.
          call mlprj_performInterpolation (p_rprojection,p_rvectorCoarse, &
                                           p_rvectorFine,p_rvectorTemp)

          ! Apply the filter chain to the temp vector.
          ! THis implements the boundary conditions that are attached to it.
          call filter_applyFilterChainVec (p_rvectorCoarse, RfilterChain)

        end if
        
        ! The system matrix looks like:
        !   (  A    0   B1 )
        !   (  0    A   B2 )
        !   ( B1^T B2^T 0  )
        !
        ! The A-matrix consists of Laplace+Convection.
        ! We build them separately and add together.
        !
        ! So at first, initialise the A-matrix with the Laplace contribution.
        ! We ignore the structure and simply overwrite the content of the
        ! system submatrices with the Laplace matrix.
        call lsyssc_duplicateMatrix (p_rmatrixLaplace,p_rmatrix%RmatrixBlock(1,1),&
                                     LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPY)

        ! Call the upwind method to calculate the nonlinear matrix.
        call conv_upwind2d (p_rvectorCoarse, p_rvectorCoarse, 1.0_DP, 0.0_DP,&
                            rupwind, CONV_MODMATRIX, &
                            p_rmatrix%RmatrixBlock(1,1))
                                     
        ! Apply the filter chain to the matrix.
        ! As the filter consists only of an implementation filter for
        ! boundary conditions, this implements the boundary conditions
        ! into the system matrix.
        call filter_applyFilterChainMat (p_rmatrix, RfilterChain)
        
      end do
      
      ! Our 'parent' (the caller of the nonlinear solver) has prepared
      ! a preconditioner node for us (a linear solver with symbolically
      ! factorised matrices). Get this from the collection.
      
      p_rsolverNode => collct_getvalue_linsol(p_rcollection,'LINSOLVER')

      ! Initialise data of the solver. This in fact performs a numeric
      ! factorisation of the matrices in UMFPACK-like solvers.
      call linsol_initData (p_rsolverNode, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      
      ! Finally solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
      call linsol_precondDefect (p_rsolverNode,rd)

      ! Release the numeric factorisation of the matrix.
      ! We don't release the symbolic factorisation, as we can use them
      ! for the next iteration.
      call linsol_doneData (p_rsolverNode)

      ! Finally calculate a new damping parameter domega.
      !
      ! For this purpose, we need two temporary vectors.
      ! On one hand, we have p_rvectorTemp.
      ! Get the second temporary vector from the collection as it was
      ! prepared by our 'parent' that invoked the nonlinear solver.
      p_rvectorTemp2 => collct_getvalue_vecsca(p_rcollection,'RTEMP2SCALAR')
      
      ! Both temp vectors are scalar, but we need block-vectors in the
      ! structure of rx/rb. Derive block vectors in that structure that
      ! share their memory with the scalar temp vectors. Note that the
      ! temp vectors are created large enough by our parent!
      call lsysbl_createVecFromScalar (p_rvectorTemp,rtemp1)
      call lsysbl_enforceStructure (rb,rtemp1)

      call lsysbl_createVecFromScalar (p_rvectorTemp2,rtemp2)
      call lsysbl_enforceStructure (rb,rtemp2)

      ! Calculate the omega
      call c2d1_getOptimalDamping (rd,rx,rb,rtemp1,rtemp2,&
                                   domega,p_rcollection)

      ! Release the temp block vectors. This only cleans up the structure.
      ! The data is not released from heap as it belongs to the
      ! scalar temp vectors.
      call lsysbl_releaseVector (rtemp2)
      call lsysbl_releaseVector (rtemp1)

    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d1_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a nonlinear solver with linear solver
  ! as preconditioner.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilvmin,ilvmax
    integer :: i
    integer :: imaxmem

    ! Error indicator during initialisation of the solver
    integer :: ierror
  
    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    type(t_filterChain), dimension(1), target :: RfilterChain

    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector
    type(t_vectorBlock), target :: rtempBlock
    type(t_vectorScalar), target :: rtempVectorSc,rtempVectorSc2

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rsmoother
    type(t_linsolNode), pointer :: p_rcoarseGridSolver,p_rpreconditioner

    ! The nonlinear solver configuration
    type(t_nlsolNode) :: rnlSol

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(NNLEV) :: Rmatrices
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection

    ! One level of multigrid
    type(t_linsolMGLevelInfo), pointer :: p_rlevelInfo

    ilvmin = rproblem%NLMIN
    ilvmax = rproblem%NLMAX
    
    ! Get our right hand side / solution / matrix on the finest
    ! level from the problem structure.
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    p_rmatrix => rproblem%RlevelInfo(ilvmax)%rmatrix
    
    ! Now we have to build up the level information for multigrid.
    !
    ! At first, initialise a standard interlevel projection structure. We
    ! can use the same structure for all levels. Therefore it's enough
    ! to initialise one structure using the RHS vector on the finest
    ! level to specify the shape of the PDE-discretisation.
    call mlprj_initProjectionVec (rprojection,rproblem%rrhs)
    
    ! During the linear solver, the boundary conditions must
    ! frequently be imposed to the vectors. This is done using
    ! a filter chain. As the linear solver does not work with
    ! the actual solution vectors but with defect vectors instead,
    ! a filter for implementing the real boundary conditions
    ! would be wrong.
    ! Therefore, create a filter chain with one filter only,
    ! which implements Dirichlet-conditions into a defect vector.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid (p_rsolverNode,RfilterChain)
    
    ! Set the output level of the solver for some output
    p_rsolverNode%ioutputLevel = 1

    ! Configure MG to gain only one digit. We can do this, as MG is only
    ! used as preconditioner inside of another solver. Furthermore, perform
    ! at least two, at most 10 steps.
    p_rsolverNode%depsRel = 1E-1_DP
    p_rsolverNode%depsAbs = 0.0_DP
    p_rsolverNode%nminIterations = 2
    p_rsolverNode%nmaxIterations = 10
    
    ! Add the interlevel projection structure to the collection; we can
    ! use it later for setting up nonlinear matrices.
    call collct_setvalue_ilvp(rproblem%rcollection,'ILVPROJECTION',&
                              rprojection,.true.)
    
    ! Then set up smoothers / coarse grid solver:
    imaxmem = 0
    do i=ilvmin,ilvmax
      
      ! On the coarsest grid, set up a coarse grid solver and no smoother
      ! On finer grids, set up a smoother but no coarse grid solver.
      nullify(p_rpreconditioner)
      nullify(p_rsmoother)
      nullify(p_rcoarseGridSolver)
      if (i .eq. ilvmin) then
        ! Set up a BiCGStab solver with VANKA preconditioning as coarse grid solver:
        !CALL linsol_initVANKA (p_rpreconditioner,1.0_DP,LINSOL_VANKA_2DSPQ1TQ0)
        !CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,RfilterChain)
        !p_rcoarseGridSolver%ioutputLevel = 2
        
        ! Setting up UMFPACK coarse grid solver would be:
        call linsol_initUMFPACK4 (p_rcoarseGridSolver)

      else
        ! Set up the VANKA smoother for multigrid with damping parameter 0.7,
        ! 4 smoothing steps:
        call linsol_initVANKA (p_rsmoother,1.0_DP,LINSOL_VANKA_2DNAVST)
        call linsol_convertToSmoother (p_rsmoother,4,0.7_DP)
      end if
    
      ! Add the level.
      call linsol_addMultigridLevel (p_rlevelInfo,p_rsolverNode, rprojection,&
                                     p_rsmoother,p_rsmoother,p_rcoarseGridSolver)

      ! How much memory is necessary for performing the level change?
      ! We ourself must build nonlinear matrices on multiple levels and have
      ! to interpolate the solution vector from finer level to coarser ones.
      ! We need temporary memory for this purpose...
      if (i .gt. ilvmin) then
        ! Pass the system metrices on the coarse/fine grid to
        ! mlprj_getTempMemoryMat to specify the discretisation structures
        ! of all equations in the PDE there.
        imaxmem = max(imaxmem,mlprj_getTempMemoryMat (rprojection,&
                              rproblem%RlevelInfo(i-1)%rmatrix,&
                              rproblem%RlevelInfo(i)%rmatrix))
      end if
    end do
    
    ! Set up a scalar temporary vector that we need for building up nonlinear
    ! matrices. It must be at least as large as MAXMEM and NEQ(finest level),
    ! as we use it for resorting vectors, too.
    call lsyssc_createVector (rtempVectorSc,max(imaxmem,rproblem%rrhs%NEQ),&
                              .false.,ST_DOUBLE)
    call collct_setvalue_vecsca(rproblem%rcollection,'RTEMPSCALAR',&
                                rtempVectorSc,.true.)
    
    ! Set up a second temporary vector that we need for calculating
    ! the optimal defect correction.
    call lsyssc_createVector (rtempVectorSc2,rproblem%rrhs%NEQ,.false.,ST_DOUBLE)
    call collct_setvalue_vecsca(rproblem%rcollection,'RTEMP2SCALAR',&
                                rtempVectorSc2,.true.)
    
    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array.
    Rmatrices(ilvmin:ilvmax) = rproblem%RlevelInfo(ilvmin:ilvmax)%rmatrix
    call linsol_setMatrices(p_rsolverNode,Rmatrices(ilvmin:ilvmax))
    
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initStructure (p_rsolverNode,ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    ! Put the prepared solver node to the collection for later use.
    call collct_setvalue_linsol(rproblem%rcollection,'LINSOLVER',p_rsolverNode,.true.)
    
    ! Create a temporary vector we need for the nonliner iteration.
    call lsysbl_createVecBlockIndirect (rproblem%rrhs, rtempBlock, .false.)

    ! The nonlinear solver structure rnlSol is initialised by the default
    ! initialisation with all necessary information to solve the problem.
    ! We call the nonlinear solver directly. For preconditioning
    ! and defect calculation, we use our own callback routine.
    rnlSol%ioutputLevel = 2
    call nlsol_performSolve(rnlSol,rproblem%rvector,rproblem%rrhs,rtempBlock,&
                            c2d1_getDefect,c2d1_precondDefect,&
                            rcollection=rproblem%rcollection)

    ! Release the temporary vector(s)
    call lsysbl_releaseVector (rtempBlock)
    call lsyssc_releaseVector (rtempVectorSc)
    call lsyssc_releaseVector (rtempVectorSc2)
    
    ! Remove the solver node from the collection - not needed anymore there
    call collct_deletevalue(rproblem%rcollection,'LINSOLVER')
    
    ! Remove the temporary vector from the collection
    call collct_deletevalue(rproblem%rcollection,'RTEMPSCALAR')
    call collct_deletevalue(rproblem%rcollection,'RTEMP2SCALAR')
    
    ! Remove the interlevel projection structure
    call collct_deletevalue(rproblem%rcollection,'ILVPROJECTION')
    
    ! Clean up the linear solver, release all memory, remove the solver node
    ! from memory.
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the multilevel projection structure.
    call mlprj_doneProjection (rprojection)
    
    print *
    print *,'Nonlinear solver statistics'
    print *,'---------------------------'
    print *,'Initial defect: ',rnlSol%DinitialDefect(1)
    print *,'Final defect:  ',rnlSol%DfinalDefect(1)
    print *,'#Iterations:   ',rnlSol%iiterations

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d1_postprocessing (rproblem)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! We need some more variables for postprocessing - i.e. writing
    ! a GMV file.
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2

    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport

    ! A pointer to the solution vector and to the triangulation.
    type(t_vectorBlock), pointer :: p_rvector
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A vector accepting Q1 data
    type(t_vectorBlock) :: rprjVector
    
    ! A discretisation structure for Q1
    type(t_blockDiscretisation) :: rprjDiscretisation
    
    ! Discrete boundary conditions for the output vector
    type(t_discreteBC), target :: rdiscreteBC

    ! A filter chain to pre-filter the vectors and the matrix.
    type(t_filterChain), dimension(1), target :: RfilterChain

    ! Get the solution vector from the problem structure.
    p_rvector => rproblem%rvector
    
    ! The solution vector is probably not in the way, GMV likes it!
    ! GMV for example does not understand Q1~ vectors!
    ! Therefore, we first have to convert the vector to a form that
    ! GMV understands.
    ! GMV understands only Q1 solutions! So the task is now to create
    ! a Q1 solution from p_rvector and write that out.
    !
    ! For this purpose, first create a 'derived' simple discretisation
    ! structure based on Q1 by copying the main guiding block discretisation
    ! structure and modifying the discretisation structures of the
    ! two velocity subvectors:
    
    call spdiscr_duplicateBlockDiscr(p_rvector%p_rblockDiscr,rprjDiscretisation)
    
    call spdiscr_deriveSimpleDiscrSc (&
                 p_rvector%p_rblockDiscr%RspatialDiscr(1), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(1))

    call spdiscr_deriveSimpleDiscrSc (&
                 p_rvector%p_rblockDiscr%RspatialDiscr(2), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(2))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (p_rvector,rprjVector)
    
    ! Discretise the boundary conditions according to the Q1/Q1/Q0
    ! discretisation:
    call c2d1_discretiseBC (rprjDiscretisation,rdiscreteBC)
                            
    ! Connect the vector to the BC's
    rprjVector%p_rdiscreteBC => rdiscreteBC
    
    ! Set up a boundary condition filter for Dirichlet boundary conditions
    ! and pass the vector through it. This finally implements the Dirichlet
    ! boundary conditions into the output vector.
    RfilterChain(1)%ifilterType = FILTER_DISCBCSOLREAL
    call filter_applyFilterChainVec (rprjVector, RfilterChain)
    
    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    !
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      p_rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing.
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing.
    ! Start UCD export to GMV file:
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,'gmv/u1.gmv')
    
    ! Write velocity field
    call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    
    call ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
    call ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)
    
    ! Write pressure
    call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
    call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Release the auxiliary vector
    call lsysbl_releaseVector (rprjVector)
    
    ! Release the discretisation structure.
    call spdiscr_releaseBlockDiscr (rprjDiscretisation)
    
    ! Throw away the discrete BC's - not used anymore.
    call bcasm_releaseDiscreteBC (rdiscreteBC)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d1_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    integer :: i

    ! Release matrices and vectors on all levels
    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Delete the matrix
      call lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrix)

      ! Delete the variables from the collection.
      call collct_deletevalue (rproblem%rcollection,'SYSTEMMAT',i)
      call collct_deletevalue (rproblem%rcollection,'LAPLACE',i)
      
      ! Release Laplace, B1 and B2 matrix
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB2)
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB1)
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixLaplace)
      
      ! Remove the temp vector that was used for interpolating the solution
      ! from higher to lower levels in the nonlinear iteration.
      if (i .lt. rproblem%NLMAX) then
        call lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVector)
        call collct_deletevalue(rproblem%rcollection,'RTEMPVEC',i)
      end if
      
    end do

    ! Delete solution/RHS vector
    call lsysbl_releaseVector (rproblem%rvector)
    call lsysbl_releaseVector (rproblem%rrhs)

    ! Delete the variables from the collection.
    call collct_deletevalue (rproblem%rcollection,'RHS')
    call collct_deletevalue (rproblem%rcollection,'SOLUTION')

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d1_doneBC (rproblem)
  
!<description>
  ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Release our discrete version of the boundary conditions
      call bcasm_releaseDiscreteBC (rproblem%RlevelInfo(i)%rdiscreteBC)
    end do
    
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine c2d1_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Remove the block discretisation structure and all the substructures.
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%p_rdiscretisation)
      
      ! Remove the discretisation from the heap.
      deallocate(rproblem%RlevelInfo(i)%p_rdiscretisation)
    end do
    
  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine c2d1_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    ! Release the triangulation on all levels
    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      call tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    end do
    
    ! Finally release the domain.
    call boundary_release (rproblem%rboundary)

    call collct_deleteValue(rproblem%rcollection,'NLMAX')
    call collct_deleteValue(rproblem%rcollection,'NLMIN')

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc2dmini1
  
!<description>
  ! This is a 'separated' Navier-Stokes solver for solving a Navier-Stokes
  ! problem. The different tasks of the problem are separated into
  ! subroutines. The problem uses a problem-specific structure for the
  ! communication: All subroutines add their generated information to the
  ! structure, so that the other subroutines can work with them.
  ! (This is somehow a cleaner implementation than using a collection!).
  ! For the communication to callback routines of black-box subroutines
  ! (matrix-assembly), a collection is used.
  !
  ! The following tasks are performed by the subroutines:
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

    ! LV receives the level where we want to solve
    integer :: NLMAX,NLMIN
    real(DP) :: dnu
    
    ! A problem structure for our problem
    type(t_problem), pointer :: p_rproblem
    
    integer :: i
    
    ! Ok, let's start.
    !
    ! Allocate the problem structure on the heap -- it's rather large.
    allocate(p_rproblem)
    
    ! NLMIN receives the minimal level where to discretise for supporting
    ! the solution process.
    ! NLMAX receives the level where we want to solve.
    NLMIN = 2
    NLMAX = 4
    
    ! Viscosity parameter:
    dnu = 1.0E0_DP/1000.0E0
    
    p_rproblem%dnu = dnu
    
    ! Initialise the collection
    call collct_init (p_rproblem%rcollection)
    do i=1,NNLEV
      call collct_addlevel (p_rproblem%rcollection)
    end do

    ! Add the (global) viscosity parameter
    call collct_setvalue_real(p_rproblem%rcollection,'NU',dnu,.true.)

    ! So now the different steps - one after the other.
    !
    ! Initialisation
    call c2d1_initParamTriang (NLMIN,NLMAX,p_rproblem)
    call c2d1_initDiscretisation (p_rproblem)
    call c2d1_initMatVec (p_rproblem)
    call c2d1_initDiscreteBC (p_rproblem)
    
    ! Implementation of boundary conditions
    call c2d1_implementBC (p_rproblem)
    
    ! Solve the problem
    call c2d1_solve (p_rproblem)
    
    ! Postprocessing
    call c2d1_postprocessing (p_rproblem)
    
    ! Cleanup
    call c2d1_doneMatVec (p_rproblem)
    call c2d1_doneBC (p_rproblem)
    call c2d1_doneDiscretisation (p_rproblem)
    call c2d1_doneParamTriang (p_rproblem)

    call collct_deletevalue(p_rproblem%rcollection,'NU')

    ! Print some statistical data about the collection - anything forgotten?
    print *
    print *,'Remaining collection statistics:'
    print *,'--------------------------------'
    print *
    call collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection and the problem structure.
    call collct_done (p_rproblem%rcollection)
    
    deallocate(p_rproblem)
    
  end subroutine

end module

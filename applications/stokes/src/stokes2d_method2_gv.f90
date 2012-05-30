!##############################################################################
!# ****************************************************************************
!# <name> stokes2d_method2_gv </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a Stokes
!# problem on a simple domain.
!#
!# The routine splits up the tasks of reading the domain, creating
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines.
!#
!# In contrast to stokes2d_method2_sv, this routine uses the general VANKA
!# for smoothing/preconditioning.
!# </purpose>
!##############################################################################

module stokes2d_method2_gv

  use fsystem
  use storage
  use boundary
  use cubature
  use derivatives
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use spdiscprojection
  use scalarpde
  use bilinearformevaluation
  use linearformevaluation
  use discretebc
  use filtersupport
  use coarsegridcorrection
  use linearsolver
  use ucd
  use pprocerror
  use genoutput
  
  use collection
    
  use stokes2d_callback
  
  implicit none
  
!<types>

!<typeblock description="Type block defining all information about one level">

  type t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...)
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Cubature info structure which encapsules the cubature formula
    type(t_scalarCubatureInfo) :: rcubatureInfo
    
    ! A system matrix for that specific level.
    type(t_matrixBlock) :: rmatrix

    ! B1-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level.
    type(t_matrixScalar) :: rmatrixB2

    ! A variable describing the discrete boundary conditions fo the velocity
    type(t_discreteBC) :: rdiscreteBC
  
  end type
  
!</typeblock>


!<typeblock description="Application-specific type block for Stokes problem">

  type t_problem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    integer :: ilvmin
    
    ! Maximum refinement level
    integer :: ilvmax
    
    ! Viscosity parameter nu = 1/Re
    real(DP) :: dnu

    ! An object for saving the domain:
    type(t_boundary) :: rboundary

    ! A solution vector and a RHS vector on the finest level.
    type(t_vectorBlock) :: rvector,rrhs

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation.
    type(t_problem_lvl), dimension(:), pointer :: RlevelInfo
    
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

  subroutine st2_initParamTriang (ilvmin,ilvmax,rproblem)
  
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
  
  ! Path to the mesh
  character(len=SYS_STRLEN) :: spredir
  
    ! Initialise the level in the problem structure
    rproblem%ilvmin = ilvmin
    rproblem%ilvmax = ilvmax
    allocate(rproblem%RlevelInfo(ilvmin:ilvmax))

    ! Get the path $PREDIR from the environment, where to read .prm/.tri files
    ! from. If that does not exist, write to the directory "./pre".
    if (.not. sys_getenv_string("PREDIR", spredir)) spredir = './pre'

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rproblem%rboundary, trim(spredir)//'/QUAD.prm')
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation, &
        trim(spredir)//'/QUAD.tri', rproblem%rboundary)
    
    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering(rproblem%ilvmin-1,&
        rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,rproblem%rboundary)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (&
        rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,rproblem%rboundary)
    
    ! Now, refine to level up to nlmax.
    do i=rproblem%ilvmin+1,rproblem%ilvmax
      call tria_refine2LevelOrdering (rproblem%RlevelInfo(i-1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary)
      call tria_initStandardMeshFromRaw (rproblem%RlevelInfo(i)%rtriangulation,&
          rproblem%rboundary)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine st2_initDiscretisation (rproblem)
  
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
    
    ! A cubature information structure on one level
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
    
    do i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rproblem%rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Get the cubature information structure
      p_rcubatureInfo => rproblem%RlevelInfo(i)%rcubatureInfo
      
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
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
                  EL_EM30, p_rtriangulation, p_rboundary)
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards.
      call spdiscr_duplicateDiscrSc(p_rdiscretisation%RspatialDiscr(1),&
          p_rdiscretisation%RspatialDiscr(2))
  
      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (p_rdiscretisation%RspatialDiscr(1), &
          EL_Q0, p_rdiscretisation%RspatialDiscr(3))
      
      ! Set up an cubature info structure to tell the code which cubature
      ! formula to use
      ! Create an assembly information structure which tells the code
      ! the cubature formula to use. Standard: Gauss 3x3.
      call spdiscr_createDefCubStructure(&  
          p_rdiscretisation%RspatialDiscr(1),rproblem%RlevelInfo(i)%rcubatureInfo,&
          CUB_GEN_AUTO_G3)
    end do
                                   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine st2_initMatVec (rproblem)
  
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
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
    ! A cubature information structure
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
    
    do i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! Get the cubature information structure
      p_rcubatureInfo => rproblem%RlevelInfo(i)%rcubatureInfo
            
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      call lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)
      
      ! Inform the matrix that we build a saddle-point problem.
      ! Normally, imatrixSpec has the value LSYSBS_MSPEC_GENERAL,
      ! but probably some solvers can use the special structure later.
      p_rmatrix%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
      
      ! Save matrix to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      call collct_setvalue_mat(rproblem%rcollection,'LAPLACE',p_rmatrix,.true.,i)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      !
      ! Create the matrix structure of the X-velocity.
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrix%RmatrixBlock(1,1))

      ! In the Stokes problem, the matrix for the Y-velocity is identical to
      ! the matrix for the X-velocity; both are Laplace-matrices!
      ! Therefore, we can simply make a copy of the matrix for the X-velocity.
      ! This we do later after the entries are created.
      !
      ! In the global system, there are two coupling matrices B1 and B2.
      ! Both have the same structure.
      !
      !    ( A         B1 )
      !    (      A    B2 )
      !    ( B1^T B2^T    )
      !
      ! Create the matrices structure of the pressure using the 3rd
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(3),LSYSSC_MATRIX9,&
                rproblem%RlevelInfo(i)%rmatrixB1,&
                p_rdiscretisation%RspatialDiscr(1))
                
      ! Duplicate the B1 matrix structure to the B2 matrix, so use
      ! lsyssc_duplicateMatrix to create B2. Share the matrix
      ! structure between B1 and B2 (B1 is the parent and B2 the child).
      ! Do not create a content array yet, it will be created by
      ! the assembly routines later.
      call lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB1,&
                  rproblem%RlevelInfo(i)%rmatrixB2,LSYSSC_DUP_COPY,LSYSSC_DUP_REMOVE)
                                       
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
      ! We specify the callback function coeff_Laplace for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine,
      ! so the callback routine has access to everything what is
      ! in the collection.
      !
      ! Build the X-velocity matrix:
      call bilf_buildMatrixScalar (rform,.true.,p_rmatrix%RmatrixBlock(1,1),&
          p_rcubatureInfo,coeff_Stokes_2D,rproblem%rcollection)
      
      ! Duplicate the matrix to the Y-velocity matrix, share structure and
      ! content between them (as the matrices are the same).
      call lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
          p_rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Manually change the discretisation structure of the Y-velocity
      ! matrix to the Y-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      call lsyssc_assignDiscrDirectMat (p_rmatrix%RmatrixBlock(2,2),&
          p_rdiscretisation%RspatialDiscr(2))
                                  
      ! Build the first pressure matrix B1.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_X

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = -1.0_DP
      
      call bilf_buildMatrixScalar (rform,.true.,rproblem%RlevelInfo(i)%rmatrixB1,&
          p_rcubatureInfo,coeff_Pressure_2D,rproblem%rcollection)

      ! Build the second pressure matrix B2.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .true.
      rform%BconstantCoeff = .true.
      rform%Dcoefficients(1)  = -1.0_DP
      
      call bilf_buildMatrixScalar (rform,.true.,rproblem%RlevelInfo(i)%rmatrixB2,&
          p_rcubatureInfo,coeff_Pressure_2D,rproblem%rcollection)
                                  
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
      
      ! Furthermore, put B1^T and B2^T to the block matrix.!
      ! Note that we really create copies of our matrices, as the
      ! general VANCA cannot handle virtually transposed matrices!
      call lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
                                   p_rmatrix%RmatrixBlock(3,1),&
                                   LSYSSC_TR_ALL)

      call lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
                                   p_rmatrix%RmatrixBlock(3,2),&
                                   LSYSSC_TR_ALL)

    end do

    ! (Only) on the finest level, we need to calculate a RHS vector
    ! and to allocate a solution vector.
    
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector

    ! Next step: Create a RHS vector and a solution vector and a temporary
    ! vector. All are filled with zero.
    call lsysbl_createVectorBlock (p_rdiscretisation,p_rrhs,.true.)
    call lsysbl_createVectorBlock (p_rdiscretisation,p_rvector,.true.)

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
    
    ! ... and then discretise the RHS to the first two subvectors of
    ! the block vector using the discretisation structure of the
    ! corresponding blocks.
    !
    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is
    ! in the collection.
    !
    ! Note that the vector is unsorted after calling this routine!
    
    p_rcubatureInfo => rproblem%RlevelInfo(rproblem%ilvmax)%rcubatureInfo
    
    call linf_buildVectorScalar (&
        rlinform,.true.,p_rrhs%RvectorBlock(1),&
        p_rcubatureInfo,coeff_RHS_X_2D,rproblem%rcollection)

    call linf_buildVectorScalar (&
        rlinform,.true.,p_rrhs%RvectorBlock(2),&
        p_rcubatureInfo,coeff_RHS_Y_2D,rproblem%rcollection)
                                
    ! The third subvector must be zero - as it represents the RHS of
    ! the equation "div(u) = 0".
    call lsyssc_clearVector(p_rrhs%RvectorBlock(3))
                                
    ! Clear the solution vector on the finest level.
    call lsysbl_clearVector(rproblem%rvector)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine st2_discretiseBC (rdiscretisation,rdiscreteBC)
  
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
    !   Dirichlet BC`s for the first (here: one and only) component in the
    !   solution vector.
    ! - Discretise the boundary condition so that the BC`s can be applied
    !   to matrices and vectors
    ! - Add the calculated discrete BC`s to rdiscreteBC for later use.
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_2D)
                              
    ! Edge 2 is Neumann boundary, so it is commented out.
    ! CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
    !                                    rboundaryRegion,rdiscreteBC,&
    !                                    getBoundaryValues_2D)
                              
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_2D)

    ! Now continue with defining the boundary conditions of the Y-velocity:
    !
    ! Define edge 1.
    call boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    
    ! Edge with start- and endpoint.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_2D)
                              
    ! Edge 2 is Neumann boundary, so it is commented out.
    ! CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rdiscretisation,2,&
    !                                    rboundaryRegion,rdiscreteBC,&
    !                                    getBoundaryValues_2D)
                              
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That is it.
    call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_2D)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine st2_initDiscreteBC (rproblem)
  
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

  ! Pointer to structure for saving discrete BC`s:
  type(t_discreteBC), pointer :: p_rdiscreteBC
    
    do i=rproblem%ilvmin,rproblem%ilvmax
    
      ! Get our velocity matrix from the problem structure.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => p_rmatrix%p_rblockDiscrTrial
      
      ! For implementing boundary conditions, we use a `filter technique with
      ! discretised boundary conditions`. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
      call bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%rdiscreteBC)
      
      ! Discretise the boundary conditions.
      call st2_discretiseBC (p_rdiscretisation,rproblem%RlevelInfo(i)%rdiscreteBC)
                               
      ! Assign the BC`s to the vectors and the matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%rdiscreteBC

      call lsysbl_assignDiscreteBC(p_rmatrix,p_rdiscreteBC)
      
    end do

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%ilvmax)%rdiscreteBC
    
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    
    call lsysbl_assignDiscreteBC(p_rrhs,p_rdiscreteBC)
    call lsysbl_assignDiscreteBC(p_rvector,p_rdiscreteBC)
                
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine st2_implementBC (rproblem)
  
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
    ilvmax = rproblem%ilvmax
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    
    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
    call vecfil_discreteBCrhs (p_rrhs)
    call vecfil_discreteBCsol (p_rvector)

    ! Implement discrete boundary conditions into the matrices on all
    ! levels, too. Call the appropriate matrix filter to modify
    ! all matrices according to the attached discrete boundary conditions.
    do i=rproblem%ilvmin,rproblem%ilvmax
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      call matfil_discreteBC (p_rmatrix)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine st2_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
    integer :: ilvmin,ilvmax
    integer :: i

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

    ! A solver node that accepts parameters for the linear solver
    type(t_linsolNode), pointer :: p_rsolverNode,p_rsmoother
    type(t_linsolNode), pointer :: p_rcoarseGridSolver,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices
    
    ! One level of multigrid
    type(t_linsolMG2LevelInfo), pointer :: p_rlevelInfo
    
    ilvmin = rproblem%ilvmin
    ilvmax = rproblem%ilvmax
    
    ! Get our right hand side / solution / matrix on the finest
    ! level from the problem structure.
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    p_rmatrix => rproblem%RlevelInfo(ilvmax)%rmatrix
    
    ! Create a temporary vector we need that for some preparation.
    call lsysbl_createVecBlockIndirect (p_rrhs, rtempBlock, .false.)
    
    ! During the linear solver, the boundary conditions must
    ! frequently be imposed to the vectors. This is done using
    ! a filter chain. As the linear solver does not work with
    ! the actual solution vectors but with defect vectors instead,
    ! a filter for implementing the real boundary conditions
    ! would be wrong.
    ! Therefore, create a filter chain with one filter only,
    ! which implements Dirichlet-conditions into a defect vector.
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Now we have to build up the level information for multigrid.
    !
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    call linsol_initMultigrid2 (p_rsolverNode,ilvmax-ilvmin+1,RfilterChain)
    
    ! Then set up smoothers / coarse grid solver:
    do i=ilvmin,ilvmax
      
      ! On the coarsest grid, set up a coarse grid solver and no smoother
      ! On finer grids, set up a smoother but no coarse grid solver.
      nullify(p_rpreconditioner)
      nullify(p_rsmoother)
      nullify(p_rcoarseGridSolver)
      if (i .eq. ilvmin) then
        ! Set up a BiCGStab solver with VANKA preconditioning as coarse grid solver:
        call linsol_initVANKA (p_rpreconditioner)
        call linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,RfilterChain)
        !p_rcoarseGridSolver%ioutputLevel = 2
        
        ! Setting up UMFPACK coarse grid solver would be:
        ! CALL linsol_initUMFPACK4 (p_rcoarseGridSolver)

      else
        ! Set up the VANKA smoother for multigrid with damping parameter 0.7,
        ! 4 smoothing steps:
        call linsol_initVANKA (p_rsmoother)
        call linsol_convertToSmoother (p_rsmoother,4,0.7_DP)
      end if
    
      ! And add this multi-grid level. We will use the same smoother
      ! for pre- and post-smoothing.
      call linsol_getMultigrid2Level (p_rsolverNode,i-ilvmin+1,p_rlevelInfo)
      p_rlevelInfo%p_rcoarseGridSolver => p_rcoarseGridSolver
      p_rlevelInfo%p_rpresmoother => p_rsmoother
      p_rlevelInfo%p_rpostsmoother => p_rsmoother
    end do
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array. Note that this does not
    ! allocate new memory, we create only 'links' to existing matrices
    ! into Rmatrices(:)!
    allocate(Rmatrices(ilvmin:ilvmax))
    do i=ilvmin,ilvmax
      call lsysbl_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end do
    
    call linsol_setMatrices(p_RsolverNode,Rmatrices(ilvmin:ilvmax))
    
    ! We can release Rmatrices immediately -- as long as we do not
    ! release rproblem%RlevelInfo(i)%rmatrix!
    do i=ilvmin,ilvmax
      call lsysbl_releaseMatrix (Rmatrices(i))
    end do
    deallocate(Rmatrices)
    
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    call linsol_initStructure (p_rsolverNode, ierror)
    
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix structure invalid!",OU_CLASS_ERROR)
      call sys_halt()
    end if

    call linsol_initData (p_rsolverNode, ierror)
    
    if (ierror .ne. LINSOL_ERR_NOERROR) then
      call output_line("Matrix singular!",OU_CLASS_ERROR)
      call sys_halt()
    end if
    
    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    call linsol_solveAdaptively (p_rsolverNode,p_rvector,p_rrhs,rtempBlock)
    
    ! Release solver data and structure
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    call linsol_releaseSolver (p_rsolverNode)
    
    ! Release the temporary vector
    call lsysbl_releaseVector (rtempBlock)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine st2_postprocessing (rproblem)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! We need some more variables for postprocessing
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    character(len=SYS_STRLEN) :: sucddir

    ! A pointer to the solution vector and to the triangulation.
    type(t_vectorBlock), pointer :: p_rvector
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A vector accepting Q1 data
    type(t_vectorBlock) :: rprjVector
    
    ! A discretisation structure for Q1
    type(t_blockDiscretisation) :: rprjDiscretisation
    
    ! Discrete boundary conditions for the output vector
    type(t_discreteBC), target :: rdiscreteBC

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
    
    call spdiscr_duplicateBlockDiscr (p_rvector%p_rblockDiscr,rprjDiscretisation)
    
    call spdiscr_deriveSimpleDiscrSc (&
                 p_rvector%p_rblockDiscr%RspatialDiscr(1), &
                 EL_Q1, rprjDiscretisation%RspatialDiscr(1))

    call spdiscr_deriveSimpleDiscrSc (&
                 p_rvector%p_rblockDiscr%RspatialDiscr(2), &
                 EL_Q1, rprjDiscretisation%RspatialDiscr(2))
                 
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
    call bcasm_initDiscreteBC(rdiscreteBC)
    call st2_discretiseBC (rprjDiscretisation,rdiscreteBC)
                            
    ! Connect the vector to the BC`s
    rprjVector%p_rdiscreteBC => rdiscreteBC
    
    ! Send the vector to the boundary-condition implementation filter.
    ! This modifies the vector according to the attached discrete boundary
    ! conditions.
    call vecfil_discreteBCsol (rprjVector)

    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    !
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      p_rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation
    
    ! Get the path for writing postprocessing files from the environment variable
    ! $UCDDIR. If that does not exist, write to the directory "./gmv".
    if (.not. sys_getenv_string("UCDDIR", sucddir)) sucddir = './gmv'

    ! p_rvector now contains our solution. We can now
    ! start the postprocessing.
    ! Start UCD export to VTK file:
    call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
                       trim(sucddir)//'/u2d_2_gv.vtk')

    ! Write velocity field
    call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    
    ! In case we use the VTK exporter, which supports vector output, we will
    ! pass the X- and Y-velocity at once to the ucd module.
    call ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)

    ! If we use the GMV exporter, we might replace the line above by the
    ! following two lines:
    !CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
    !CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)
    
    ! Write pressure
    call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
    call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)
    
    ! Release the auxiliary vector
    call lsysbl_releaseVector (rprjVector)
    
    ! Release the discretisation structure.
    call spdiscr_releaseBlockDiscr (rprjDiscretisation)
    
    ! Throw away the discrete BC`s - not used anymore.
    call bcasm_releaseDiscreteBC (rdiscreteBC)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine st2_calcErrors (rproblem)
  
!<description>
  ! Calculates L2- and H1-errors against analytic solution
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! Error data structures for post-processing
  type(t_errorScVec) :: rerrorU, rerrorP

  ! A cubature information structure
  type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

  ! Error arrays for post-processing
  real(DP), dimension(2), target :: DerrorUL2, DerrorUH1
  real(DP), dimension(1), target :: DerrorPL2

    ! Get the cubature information structure
    p_rcubatureInfo => rproblem%RlevelInfo(rproblem%ilvmax)%rcubatureInfo
    
    ! Store the viscosity parameter nu in the collection's quick access array
    rproblem%rcollection%DquickAccess(1) = rproblem%dnu

    ! Set up the error structure for velocity
    rerrorU%p_RvecCoeff => rproblem%rvector%RvectorBlock(1:2)
    rerrorU%p_DerrorL2 => DerrorUL2
    rerrorU%p_DerrorH1 => DerrorUH1

    ! Set up the error structure for pressure
    rerrorP%p_RvecCoeff => rproblem%rvector%RvectorBlock(3:3)
    rerrorP%p_DerrorL2 => DerrorPL2

    ! Calculate errors of velocity and pressure against analytic solutions.
    call pperr_scalarVec(rerrorU, funcVelocity2D, rproblem%rcollection, rcubatureInfo=p_rcubatureInfo)
    call pperr_scalarVec(rerrorP, funcPressure2D, rproblem%rcollection, rcubatureInfo=p_rcubatureInfo)

    ! Print the errors.
    call output_lbrk()
    call output_line('|u - u_h|_L2 = ' // trim(sys_sdEL(DerrorUL2(1), 10)) &
                                // ' ' // trim(sys_sdEL(DerrorUL2(2), 10)))
    call output_line('|u - u_h|_H1 = ' // trim(sys_sdEL(DerrorUH1(1), 10)) &
                                // ' ' // trim(sys_sdEL(DerrorUH1(2), 10)))
    call output_line('|p - p_h|_L2 = ' // trim(sys_sdEL(derrorPL2(1), 10)))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine st2_doneMatVec (rproblem)
  
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
    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Delete the matrix
      call lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrix)

      ! Delete the variables from the collection.
      call collct_deletevalue (rproblem%rcollection,'LAPLACE',i)
      
      ! Release B1 and B2 matrix
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB2)
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB1)
      
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

  subroutine st2_doneBC (rproblem)
  
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

    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release our discrete version of the boundary conditions
      call bcasm_releaseDiscreteBC (rproblem%RlevelInfo(i)%rdiscreteBC)
    end do
    
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine st2_doneDiscretisation (rproblem)
  
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

    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release the cubature info structure.
      call spdiscr_releaseCubStructure(rproblem%RlevelInfo(i)%rcubatureInfo)

      ! Remove the block discretisation structure and all substructures.
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%p_rdiscretisation)
      
      ! Remove the discretisation from the heap.
      deallocate(rproblem%RlevelInfo(i)%p_rdiscretisation)
    end do
    
  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine st2_doneParamTriang (rproblem)
  
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

    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release the triangulation
      call tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    end do
    
    deallocate(rproblem%RlevelInfo)
    
    ! Finally release the domain.
    call boundary_release (rproblem%rboundary)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine stokes2d_2_gv
  
!<description>
  ! This is a 'separated' stokes solver for solving a stokes
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

    ! NLMIN receives the minimal level where to discretise for supporting
    ! the solution process.
    ! NLMAX receives the level where we want to solve.
    integer :: NLMIN,NLMAX
    real(DP) :: dnu
    
    ! A problem structure for our problem
    type(t_problem), target :: rproblem
    
    integer :: i
    
    ! Ok, let us start.
    ! We want to solve our Laplace problem on level...

    NLMIN = 2
    NLMAX = 7
    
    ! Viscosity parameter:
    dnu = 1E0_DP
    
    rproblem%dnu = dnu
    
    ! Initialise the collection
    call collct_init (rproblem%rcollection)
    do i=1,NLMAX
      call collct_addlevel (rproblem%rcollection)
    end do

    ! So now the different steps - one after the other.
    !
    ! Initialisation
    call st2_initParamTriang (NLMIN,NLMAX,rproblem)
    call st2_initDiscretisation (rproblem)
    call st2_initMatVec (rproblem)
    call st2_initDiscreteBC (rproblem)
    
    ! Implementation of boundary conditions
    call st2_implementBC (rproblem)
    
    ! Solve the problem
    call st2_solve (rproblem)
    
    ! Postprocessing
    call st2_postprocessing (rproblem)
    call st2_calcErrors(rproblem)
    
    ! Cleanup
    call st2_doneMatVec (rproblem)
    call st2_doneBC (rproblem)
    call st2_doneDiscretisation (rproblem)
    call st2_doneParamTriang (rproblem)

    ! Print some statistical data about the collection - anything forgotten?
    print *
    print *,'Remaining collection statistics:'
    print *,'--------------------------------'
    print *
    call collct_printStatistics (rproblem%rcollection)
    
    ! Finally release the collection.
    call collct_done (rproblem%rcollection)
    
  end subroutine

end module

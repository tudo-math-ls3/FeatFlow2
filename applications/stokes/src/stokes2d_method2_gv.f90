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
!# In contrast to stokes2d_method2_sv, this routine uses the general VANCA
!# for smoothing/preconditioning.
!# </purpose>
!##############################################################################

MODULE stokes2d_method2_gv

  USE fsystem
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
  USE coarsegridcorrection
  USE spdiscprojection
  USE ucd
  
  USE collection
    
  USE stokes2d_callback
  
  IMPLICIT NONE
  
!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the block discretisation
    ! (size of subvectors in the solution vector, trial/test functions,...)
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! A system matrix for that specific level. 
    TYPE(t_matrixBlock) :: rmatrix

    ! B1-matrix for that specific level. 
    TYPE(t_matrixScalar) :: rmatrixB1

    ! B2-matrix for that specific level. 
    TYPE(t_matrixScalar) :: rmatrixB2

    ! A variable describing the discrete boundary conditions fo the velocity
    TYPE(t_discreteBC) :: rdiscreteBC
  
  END TYPE
  
!</typeblock>


!<typeblock description="Application-specific type block for stokes problem">

  TYPE t_problem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    INTEGER :: ilvmin
    
    ! Maximum refinement level
    INTEGER :: ilvmax
    
    ! Viscosity parameter nu = 1/Re
    REAL(DP) :: dnu

    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary

    ! A solution vector and a RHS vector on the finest level. 
    TYPE(t_vectorBlock) :: rvector,rrhs

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation.
    TYPE(t_problem_lvl), DIMENSION(:), POINTER :: RlevelInfo
    
    ! A collection object that saves structural data and some 
    ! problem-dependent information which is e.g. passed to 
    ! callback routines.
    TYPE(t_collection) :: rcollection

  END TYPE

!</typeblock>

!</types>
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st2_initParamTriang (ilvmin,ilvmax,rproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
!</description>

!<input>
  ! Minimum refinement level of the mesh; = coarse grid = level 1
  INTEGER, INTENT(IN) :: ilvmin
  
  ! Maximum refinement level
  INTEGER, INTENT(IN) :: ilvmax
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i
  
    ! Initialise the level in the problem structure
    rproblem%ilvmin = ilvmin
    rproblem%ilvmax = ilvmax
    ALLOCATE(rproblem%RlevelInfo(ilvmin:ilvmax))

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    ! Set p_rboundary to NULL() to create a new structure.
    NULLIFY(rproblem%p_rboundary)
    CALL boundary_read_prm(rproblem%p_rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation, &
        './pre/QUAD.tri', rproblem%p_rboundary)
    
    ! Refine the mesh up to the minimum level
    CALL tria_quickRefine2LevelOrdering(rproblem%ilvmin-1,&
        rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,rproblem%p_rboundary)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    CALL tria_initStandardMeshFromRaw (&
        rproblem%RlevelInfo(rproblem%ilvmin)%rtriangulation,rproblem%p_rboundary)
    
    ! Now, refine to level up to nlmax.
    DO i=rproblem%ilvmin+1,rproblem%ilvmax
      CALL tria_refine2LevelOrdering (rproblem%RlevelInfo(i-1)%rtriangulation,&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%p_rboundary)
      CALL tria_initStandardMeshFromRaw (rproblem%RlevelInfo(i)%rtriangulation,&
          rproblem%p_rboundary)
    END DO
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st2_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: I
  
    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! An object for the block discretisation on one level
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

    DO i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rproblem%p_rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies 3 blocks in the
      ! solution vector. In this simple problem, we only have one block.
      ALLOCATE(p_rdiscretisation)
      CALL spdiscr_initBlockDiscr2D (p_rdiscretisation,3,&
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
      CALL spdiscr_initDiscr_simple ( &
                  p_rdiscretisation%RspatialDiscretisation(1), &
                  EL_EM30,CUB_G2X2, &
                  p_rtriangulation, p_rboundary)
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory, 
      ! as both structures will share the same dynamic information afterwards.
      CALL spdiscr_duplicateDiscrSc(p_rdiscretisation%RspatialDiscretisation(1),&
          p_rdiscretisation%RspatialDiscretisation(2))
  
      ! For the pressure (3rd component), we set up a separate discretisation 
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      CALL spdiscr_initDiscr_combined ( &
                  p_rdiscretisation%RspatialDiscretisation(3), &
                  EL_Q0,EL_EM30,CUB_G2X2, &
                  p_rtriangulation, p_rboundary)
    END DO
                                   
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st2_initMatVec (rproblem)
  
!<description>
  ! Calculates the system matrix and RHS vector of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i
  
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector

    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
  
    DO i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)    
      
      ! Inform the matrix that we build a saddle-point problem.
      ! Normally, imatrixSpec has the value LSYSBS_MSPEC_GENERAL,
      ! but probably some solvers can use the special structure later.
      p_rmatrix%imatrixSpec = LSYSBS_MSPEC_SADDLEPOINT
      
      ! Save matrix to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      CALL collct_setvalue_mat(rproblem%rcollection,'LAPLACE',p_rmatrix,.TRUE.,i)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      !
      ! Create the matrix structure of the X-velocity.
      CALL bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscretisation(1),LSYSSC_MATRIX9,&
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
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscretisation.
      CALL bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscretisation(3),LSYSSC_MATRIX9,&
                rproblem%RlevelInfo(i)%rmatrixB1)
                
      ! Duplicate the B1 matrix structure to the B2 matrix, so use
      ! lsyssc_duplicateMatrix to create B2. Share the matrix 
      ! structure between B1 and B2 (B1 is the parent and B2 the child). 
      ! Don't create a content array yet, it will be created by 
      ! the assembly routines later.
      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB1,&
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
      rform%ballCoeffConstant = .TRUE.
      rform%BconstantCoeff = .TRUE.
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
      CALL bilf_buildMatrixScalar (rform,.TRUE.,p_rmatrix%RmatrixBlock(1,1),&
                                   coeff_Stokes_2D,rproblem%rcollection)
      
      ! Duplicate the matrix to the Y-velocity matrix, share structure and
      ! content between them (as the matrices are the same).
      CALL lsyssc_duplicateMatrix (p_rmatrix%RmatrixBlock(1,1),&
                  p_rmatrix%RmatrixBlock(2,2),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Manually change the discretisation structure of the Y-velocity 
      ! matrix to the Y-discretisation structure.
      ! Ok, we use the same discretisation structure for both, X- and Y-velocity,
      ! so this is not really necessary - we do this for sure...
      p_rmatrix%RmatrixBlock(2,2)%p_rspatialDiscretisation => &
        p_rdiscretisation%RspatialDiscretisation(2)
                                  
      ! Build the first pressure matrix B1.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_X

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .TRUE.
      rform%BconstantCoeff = .TRUE.
      rform%Dcoefficients(1)  = -1.0_DP
      
      CALL bilf_buildMatrixScalar (rform,.TRUE.,rproblem%RlevelInfo(i)%rmatrixB1,&
                                   coeff_Pressure_2D,rproblem%rcollection)

      ! Build the second pressure matrix B2.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .TRUE.
      rform%BconstantCoeff = .TRUE.
      rform%Dcoefficients(1)  = -1.0_DP
      
      CALL bilf_buildMatrixScalar (rform,.TRUE.,rproblem%RlevelInfo(i)%rmatrixB2,&
                                   coeff_Pressure_2D,rproblem%rcollection)
                                  
      ! The B1/B2 matrices exist up to now only in our local problem structure.
      ! Put a copy of them into the block matrix.
      !
      ! Note that we share the structure of B1/B2 with those B1/B2 of the
      ! block matrix, while we create copies of the entries. The reason is
      ! that these matrices are modified for bondary conditions later.
      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
                                   p_rmatrix%RmatrixBlock(1,3),&
                                   LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
                                   p_rmatrix%RmatrixBlock(2,3),&
                                   LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
      ! Furthermore, put B1^T and B2^T to the block matrix.!
      ! Note that we really create copies of our matrices, as the
      ! general VANCA cannot handle virtually transposed matrices!
      CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
                                   p_rmatrix%RmatrixBlock(3,1),&
                                   LSYSSC_TR_ALL)

      CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
                                   p_rmatrix%RmatrixBlock(3,2),&
                                   LSYSSC_TR_ALL)

      ! Update the structural information of the block matrix, as we manually
      ! changed the submatrices:
      CALL lsysbl_updateMatStrucInfo (p_rmatrix)
      
    END DO

    ! (Only) on the finest level, we need to calculate a RHS vector
    ! and to allocate a solution vector.
    
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector

    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    CALL lsysbl_createVecBlockIndMat (p_rmatrix,p_rrhs, .FALSE.)
    CALL lsysbl_createVecBlockIndMat (p_rmatrix,p_rvector, .FALSE.)

    ! Save the solution/RHS vector to the collection. Might be used
    ! later (e.g. in nonlinear problems)
    CALL collct_setvalue_vec(rproblem%rcollection,'RHS',p_rrhs,.TRUE.)
    CALL collct_setvalue_vec(rproblem%rcollection,'SOLUTION',p_rvector,.TRUE.)
    
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
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(1),rlinform,.TRUE.,&
              p_rrhs%RvectorBlock(1),coeff_RHS_X_2D,&
              rproblem%rcollection)

    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(1),rlinform,.TRUE.,&
              p_rrhs%RvectorBlock(2),coeff_RHS_Y_2D,&
              rproblem%rcollection)
                                
    ! The third subvector must be zero - as it represents the RHS of
    ! the equation "div(u) = 0".
    CALL lsyssc_clearVector(p_rrhs%RvectorBlock(3))
                                
    ! Clear the solution vector on the finest level.
    CALL lsysbl_clearVector(rproblem%rvector)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st2_discretiseBC (rdiscretisation,rdiscreteBC)
  
!<description>
  ! This routine discretises the current boundary conditions for
  ! the discretisation specified by rdiscretisation.
!</description>

!<input>
  ! A discretisation structure specifying the current discretisation.
  TYPE(t_blockDiscretisation), INTENT(IN) :: rdiscretisation
!</input>

!<inputoutput>
  ! A structuree that receives the discretised boundary conditions.
  TYPE(t_discreteBC), INTENT(INOUT) :: rdiscreteBC
!</inputoutput>

!</subroutine>

    ! local variables
    TYPE(t_boundaryRegion) :: rboundaryRegion    
    TYPE(t_boundary), POINTER :: p_rboundary
    
    p_rboundary => rdiscretisation%p_rboundary

    ! We first set up the boundary conditions for the X-velocity, then those
    ! of the Y-velocity.
    !
    ! We 'know' already (from the problem definition) that we have four boundary
    ! segments in the domain. Each of these, we want to use for inforcing
    ! some kind of boundary condition.
    !
    ! We ask the bondary routines to create a 'boundary region' - which is
    ! simply a part of the boundary corresponding to a boundary segment.
    ! A boundary region roughly contains the type, the min/max parameter value
    ! and whether the endpoints are inside the region or not.
    CALL boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    
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
    CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_2D)
                              
    ! Edge 2 is Neumann boundary, so it's commented out.
    ! CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
    !                                    rboundaryRegion,rdiscreteBC,&
    !                                    getBoundaryValues_2D)
                              
    ! Edge 3 of boundary component 1.
    CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_2D)

    ! Now continue with defining the boundary conditions of the Y-velocity:
    !
    ! Define edge 1.
    CALL boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    
    ! Edge with start- and endpoint.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    CALL bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_2D)
                              
    ! Edge 2 is Neumann boundary, so it's commented out.
    ! CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rdiscretisation,2,&
    !                                    rboundaryRegion,rdiscreteBC,&
    !                                    getBoundaryValues_2D)
                              
    ! Edge 3 of boundary component 1.
    CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    CALL bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_2D)
    
    ! Edge 4 of boundary component 1. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues_2D)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st2_initDiscreteBC (rproblem)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

  ! A pointer to the system matrix and the RHS vector as well as 
  ! the discretisation
  TYPE(t_matrixBlock), POINTER :: p_rmatrix
  TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
  TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

  ! Pointer to structure for saving discrete BC's:
  TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
    
    DO i=rproblem%ilvmin,rproblem%ilvmax
    
      ! Get our velocity matrix from the problem structure.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => p_rmatrix%p_rblockDiscretisation
      
      ! For implementing boundary conditions, we use a 'filter technique with
      ! discretised boundary conditions'. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
      CALL bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%rdiscreteBC)
      
      ! Discretise the boundary conditions.
      CALL st2_discretiseBC (p_rdiscretisation,rproblem%RlevelInfo(i)%rdiscreteBC)
                               
      ! Hang the pointer into the the matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%rdiscreteBC
      
      p_rmatrix%p_rdiscreteBC => p_rdiscreteBC
      
    END DO

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%ilvmax)%rdiscreteBC
    
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector
    
    p_rrhs%p_rdiscreteBC => p_rdiscreteBC
    p_rvector%p_rdiscreteBC => p_rdiscreteBC
                
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st2_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,ilvmax
  
    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
    
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
    CALL vecfil_discreteBCrhs (p_rrhs)
    CALL vecfil_discreteBCsol (p_rvector)

    ! Implement discrete boundary conditions into the matrices on all 
    ! levels, too. Call the appropriate matrix filter to modify
    ! all matrices according to the attached discrete boundary conditions.
    DO i=rproblem%ilvmin,rproblem%ilvmax
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      CALL matfil_discreteBC (p_rmatrix)
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st2_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
    INTEGER :: ilvmin,ilvmax
    INTEGER :: i

    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    
  
    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain

    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
    TYPE(t_vectorBlock), TARGET :: rtempBlock

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rsmoother
    TYPE(t_linsolNode), POINTER :: p_rcoarseGridSolver,p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(:), POINTER :: Rmatrices
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock) :: rprojection

    ! One level of multigrid
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo
    
    ilvmin = rproblem%ilvmin
    ilvmax = rproblem%ilvmax
    
    ! Get our right hand side / solution / matrix on the finest
    ! level from the problem structure.
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector
    p_rmatrix => rproblem%RlevelInfo(ilvmax)%rmatrix
    
    ! Create a temporary vector we need that for some preparation.
    CALL lsysbl_createVecBlockIndirect (p_rrhs, rtempBlock, .FALSE.)
    
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
    ! At first, initialise a standard interlevel projection structure. We
    ! can use the same structure for all levels.
    CALL mlprj_initProjectionMat (rprojection,p_rmatrix)
    
    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    CALL linsol_initMultigrid (p_rsolverNode,p_RfilterChain)
    
    ! Then set up smoothers / coarse grid solver:
    DO i=ilvmin,ilvmax
      
      ! On the coarsest grid, set up a coarse grid solver and no smoother
      ! On finer grids, set up a smoother but no coarse grid solver.
      NULLIFY(p_rpreconditioner)
      NULLIFY(p_rsmoother)
      NULLIFY(p_rcoarseGridSolver)
      IF (i .EQ. ilvmin) THEN
        ! Set up a BiCGStab solver with VANCA preconditioning as coarse grid solver:
        CALL linsol_initVANCA (p_rpreconditioner)
        CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,p_RfilterChain)
        !p_rcoarseGridSolver%ioutputLevel = 2
        
        ! Setting up UMFPACK coarse grid solver would be:
        ! CALL linsol_initUMFPACK4 (p_rcoarseGridSolver)

      ELSE
        ! Set up the VANCA smoother for multigrid with damping parameter 0.7,
        ! 4 smoothing steps:
        CALL linsol_initVANCA (p_rsmoother)
        CALL linsol_convertToSmoother (p_rsmoother,4,0.7_DP)
      END IF
    
      ! Add the level.
      CALL linsol_addMultigridLevel (p_rlevelInfo,p_rsolverNode, rprojection,&
                                     p_rsmoother,p_rsmoother,p_rcoarseGridSolver)
    END DO
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array. Note that this does not
    ! allocate new memory, we create only 'links' to existing matrices
    ! into Rmatrices(:)!
    ALLOCATE(Rmatrices(ilvmin:ilvmax))
    DO i=ilvmin,ilvmax
      CALL lsysbl_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrix,&
          Rmatrices(i),LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    END DO
    
    CALL linsol_setMatrices(p_RsolverNode,Rmatrices(ilvmin:ilvmax))
    
    ! We can release Rmatrices immediately -- as long as we don't
    ! release rproblem%RlevelInfo(i)%rmatrix!
    DO i=ilvmin,ilvmax
      CALL lsysbl_releaseMatrix (Rmatrices(i))
    END DO
    DEALLOCATE(Rmatrices)
    
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    CALL linsol_initStructure (p_rsolverNode,ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    CALL linsol_initData (p_rsolverNode,ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    
    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    CALL linsol_solveAdaptively (p_rsolverNode,p_rvector,p_rrhs,rtempBlock)
    
    ! Release solver data and structure
    CALL linsol_doneData (p_rsolverNode)
    CALL linsol_doneStructure (p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (p_rsolverNode)
    
    ! Release the multilevel projection structure.
    CALL mlprj_doneProjection (rprojection)
    
    ! Release the temporary vector
    CALL lsysbl_releaseVector (rtempBlock)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st2_postprocessing (rproblem)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! We need some more variables for postprocessing
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport

    ! A pointer to the solution vector and to the triangulation.
    TYPE(t_vectorBlock), POINTER :: p_rvector
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! A vector accepting Q1 data
    TYPE(t_vectorBlock) :: rprjVector
    
    ! A discretisation structure for Q1
    TYPE(t_blockDiscretisation) :: rprjDiscretisation
    
    ! Discrete boundary conditions for the output vector
    TYPE(t_discreteBC), TARGET :: rdiscreteBC

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
    
    CALL spdiscr_duplicateBlockDiscr (p_rvector%p_rblockDiscretisation,rprjDiscretisation)
    
    CALL spdiscr_deriveSimpleDiscrSc (&
                 p_rvector%p_rblockDiscretisation%RspatialDiscretisation(1), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscretisation(1))

    CALL spdiscr_deriveSimpleDiscrSc (&
                 p_rvector%p_rblockDiscretisation%RspatialDiscretisation(2), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscretisation(2))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    CALL lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.FALSE.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    CALL spdp_projectSolution (p_rvector,rprjVector)
    
    ! Discretise the boundary conditions according to the Q1/Q1/Q0 
    ! discretisation:
    CALL bcasm_initDiscreteBC(rdiscreteBC)
    CALL st2_discretiseBC (rprjDiscretisation,rdiscreteBC)
                            
    ! Connect the vector to the BC's
    rprjVector%p_rdiscreteBC => rdiscreteBC
    
    ! Send the vector to the boundary-condition implementation filter.
    ! This modifies the vector according to the attached discrete boundary
    ! conditions.
    CALL vecfil_discreteBCsol (rprjVector)

    ! Now we have a Q1/Q1/Q0 solution in rprjVector.
    !
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      p_rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing. 
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
                       'gmv/u2d_2_gv.gmv')

    ! Write velocity field
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    
    ! In case we use the VTK exporter, which supports vector output, we will
    ! pass the X- and Y-velocity at once to the ucd module.
    CALL ucd_addVarVertBasedVec(rexport,'velocity',p_Ddata,p_Ddata2)

    ! If we use the GMV exporter, we might replace the line above by the
    ! following two lines:
    !CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
    !CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)
    
    ! Write pressure
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
    CALL ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
    ! Release the auxiliary vector
    CALL lsysbl_releaseVector (rprjVector)
    
    ! Release the discretisation structure.
    CALL spdiscr_releaseBlockDiscr (rprjDiscretisation)
    
    ! Throw away the discrete BC's - not used anymore.
    CALL bcasm_releaseDiscreteBC (rdiscreteBC)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st2_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    INTEGER :: i

    ! Release matrices and vectors on all levels
    DO i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Delete the matrix
      CALL lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrix)

      ! Delete the variables from the collection.
      CALL collct_deletevalue (rproblem%rcollection,'LAPLACE',i)
      
      ! Release B1 and B2 matrix
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB2)
      CALL lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB1)
      
    END DO

    ! Delete solution/RHS vector
    CALL lsysbl_releaseVector (rproblem%rvector)
    CALL lsysbl_releaseVector (rproblem%rrhs)

    ! Delete the variables from the collection.
    CALL collct_deletevalue (rproblem%rcollection,'RHS')
    CALL collct_deletevalue (rproblem%rcollection,'SOLUTION')

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st2_doneBC (rproblem)
  
!<description>
  ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    DO i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release our discrete version of the boundary conditions
      CALL bcasm_releaseDiscreteBC (rproblem%RlevelInfo(i)%rdiscreteBC)
    END DO
    
  END SUBROUTINE


  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st2_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    DO i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Remove the block discretisation structure and all substructures.
      CALL spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%p_rdiscretisation)
      
      ! Remove the discretisation from the heap.
      DEALLOCATE(rproblem%RlevelInfo(i)%p_rdiscretisation)
    END DO
    
  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st2_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    DO i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release the triangulation
      CALL tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    END DO
    
    DEALLOCATE(rproblem%RlevelInfo)
    
    ! Finally release the domain.
    CALL boundary_release (rproblem%p_rboundary)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE stokes2d_2_gv
  
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
    INTEGER :: NLMIN,NLMAX
    REAL(DP) :: dnu
    
    ! A problem structure for our problem
    TYPE(t_problem), TARGET :: rproblem
    
    INTEGER :: i
    
    ! Ok, let's start. 
    ! We want to solve our Laplace problem on level...

    NLMIN = 2
    NLMAX = 7
    
    ! Viscosity parameter:
    dnu = 1E0_DP
    
    rproblem%dnu = dnu
    
    ! Initialise the collection
    CALL collct_init (rproblem%rcollection)
    DO i=1,NLMAX
      CALL collct_addlevel_all (rproblem%rcollection)
    END DO

    ! So now the different steps - one after the other.
    !
    ! Initialisation
    CALL st2_initParamTriang (NLMIN,NLMAX,rproblem)
    CALL st2_initDiscretisation (rproblem)    
    CALL st2_initMatVec (rproblem)    
    CALL st2_initDiscreteBC (rproblem)
    
    ! Implementation of boundary conditions
    CALL st2_implementBC (rproblem)
    
    ! Solve the problem
    CALL st2_solve (rproblem)
    
    ! Postprocessing
    CALL st2_postprocessing (rproblem)
    
    ! Cleanup
    CALL st2_doneMatVec (rproblem)
    CALL st2_doneBC (rproblem)
    CALL st2_doneDiscretisation (rproblem)
    CALL st2_doneParamTriang (rproblem)

    ! Print some statistical data about the collection - anything forgotten?
    PRINT *
    PRINT *,'Remaining collection statistics:'
    PRINT *,'--------------------------------'
    PRINT *
    CALL collct_printStatistics (rproblem%rcollection)
    
    ! Finally release the collection.
    CALL collct_done (rproblem%rcollection)
    
  END SUBROUTINE

END MODULE

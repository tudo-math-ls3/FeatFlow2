!##############################################################################
!# ****************************************************************************
!# <name> stokes_method1 </name>
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
!# The routine uses the simple-VANCA smoother/preconditioner for
!# 2D saddle point problems, Jacobi-Type.
!# </purpose>
!##############################################################################

MODULE stokes_method1

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
    
  USE stokes_callback
  
  IMPLICIT NONE
  
  ! Maximum allowed level in this application; must be =9 for 
  ! FEAT 1.x compatibility (still)!
  INTEGER, PARAMETER :: NNLEV = 9
  
!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation), POINTER :: p_rtriangulation

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
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
  
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

    ! A variable describing the analytic boundary conditions.    
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by LV!
    TYPE(t_problem_lvl), DIMENSION(NNLEV) :: RlevelInfo
    
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

  SUBROUTINE st1_initParamTriang (ilvmin,ilvmax,rproblem)
  
    INCLUDE 'cout.inc'
    INCLUDE 'cerr.inc'
    INCLUDE 'cmem.inc'
    INCLUDE 'cparametrization.inc'

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
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

  ! local variables
  INTEGER :: i
  
    ! For compatibility to old F77: an array accepting a set of triangulations
    INTEGER, DIMENSION(SZTRIA,NNLEV) :: TRIAS

    ! Variable for a filename:  
    CHARACTER(LEN=60) :: CFILE

    ! Initialise the level in the problem structure
    rproblem%ilvmin = ilvmin
    rproblem%ilvmax = ilvmax

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    ! Set p_rboundary to NULL() to create a new structure.
    NULLIFY(rproblem%p_rboundary)
    CALL boundary_read_prm(rproblem%p_rboundary, './pre/QUAD.prm')
        
    ! Remark that this does not read in the parametrisation for FEAT 1.x.
    ! Unfortunately we still need it for creating the initial triangulation!
    ! Therefore, read the file again wihh FEAT 1.x routines.
    IMESH = 1
    CFILE = './pre/QUAD.prm'
    CALL GENPAR (.TRUE.,IMESH,CFILE)

    ! Now read in the triangulation - in FEAT 1.x syntax.
    ! Refine it to level ilvmin/ilvmax.
    ! This will probably modify ilvmin/ilvmax in case of a level
    ! shift, i.e. if ilvmax > ilvmin+9 !
    ! After this routine, we have to rely on ilvmin/ilvmax in the
    ! problem structure rather than those in the parameters.
    CFILE = './pre/QUAD.tri'
    CALL INMTRI (2,TRIAS,rproblem%ilvmin,rproblem%ilvmax,0,0,CFILE)
    
    ! ... and create a FEAT 2.0 triangulation for that. Until the point where
    ! we recreate the triangulation routines, this method has to be used
    ! to get a triangulation.
    ! Set p_rtriangulation to NULL() to create a new structure on the heap.
    DO i=rproblem%ilvmin,rproblem%ilvmax
      NULLIFY(rproblem%RlevelInfo(i)%p_rtriangulation)
      CALL tria_wrp_tria2Structure(TRIAS(:,i),rproblem%RlevelInfo(i)%p_rtriangulation)
    END DO
    
    ! The TRIAS(,)-array is now part pf the triangulation structure,
    ! we don't need it anymore.
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st1_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

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
      p_rtriangulation => rproblem%RlevelInfo(i)%p_rtriangulation
      
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
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      p_rdiscretisation%RspatialDiscretisation(2) = &
        p_rdiscretisation%RspatialDiscretisation(1)
  
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

  SUBROUTINE st1_initMatVec (rproblem)
  
!<description>
  ! Calculates the system matrix and RHS vector of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

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
      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                   p_rmatrix%RmatrixBlock(1,1),coeff_Stokes,&
                                   rproblem%rcollection)
      
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
      
      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                  rproblem%RlevelInfo(i)%rmatrixB1,coeff_Pressure,&
                                  rproblem%rcollection)

      ! Build the second pressure matrix B2.
      ! Again first set up the bilinear form, then call the matrix assembly.
      rform%itermCount = 1
      rform%Idescriptors(1,1) = DER_FUNC
      rform%Idescriptors(2,1) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients:
      rform%ballCoeffConstant = .TRUE.
      rform%BconstantCoeff = .TRUE.
      rform%Dcoefficients(1)  = -1.0_DP
      
      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                  rproblem%RlevelInfo(i)%rmatrixB2,coeff_Pressure,&
                                  rproblem%rcollection)
                                  
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
      
      ! Furthermore, put B1^T and B2^T to the block matrix.
      CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB1, &
                                   p_rmatrix%RmatrixBlock(3,1),&
                                   LSYSSC_TR_VIRTUAL)

      CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB2, &
                                   p_rmatrix%RmatrixBlock(3,2),&
                                   LSYSSC_TR_VIRTUAL)

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
    
    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the 
    ! first block.
    !
    ! We pass our collection structure as well to this routine, 
    ! so the callback routine has access to everything what is
    ! in the collection.
    !
    ! Note that the vector is unsorted after calling this routine!
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(1),rlinform,.TRUE.,&
              p_rrhs%RvectorBlock(1),coeff_RHS,&
              rproblem%rcollection)
                                
    ! The third subvector must be zero - as it represents the RHS of
    ! the equation "div(u) = 0".
    CALL lsyssc_clearVector(p_rrhs%RvectorBlock(3))
                                
    ! Clear the solution vector on the finest level.
    CALL lsysbl_clearVector(rproblem%rvector)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st1_initAnalyticBC (rproblem)
  
!<description>
  ! This initialises the analytic bonudary conditions of the problem
  ! and saves them to the problem structure.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables

    ! A set of variables describing the analytic boundary conditions.    
    TYPE(t_boundaryRegion) :: rboundaryRegion
    TYPE(t_bcRegion), POINTER :: p_rbcRegion
    
    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! A pointer to the domain
    TYPE(t_boundary), POINTER :: p_rboundary
    
    INTEGER :: i
  
    ! Get the domain from the problem structure
    p_rboundary => rproblem%p_rboundary

    ! For implementing boundary conditions, we use a 'filter technique with
    ! discretised boundary conditions'. This means, we first have to calculate
    ! a discrete version of the analytic BC, which we can implement into the
    ! solution/RHS vectors using the corresponding filter.
    !
    ! At first, we need the analytic description of the boundary conditions.
    ! Initialise a structure for boundary conditions, which accepts this,
    ! on the heap.
    !
    ! We first set up the boundary conditions for the X-velocity, then those
    ! of the Y-velocity.
    !
    ! Set p_rboundaryConditions to NULL() to create a new structure on the heap.
    NULLIFY (rproblem%p_rboundaryConditions)
    CALL bcond_initBC (rproblem%p_rboundaryConditions,p_rboundary)
    
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
    ! boundary there. The following routine adds a new 'boundary condition region'
    ! for the first segment to the boundary condition structure.
    ! The region will be set up as 'Dirichlet boundary'.
    ! We specify icomponent='1' to indicate that we set up the
    ! Dirichlet BC's for the firstcomponent in the solution vector,
    ! the X-velocity.
    ! The routine also returns the created object in p_rbcRegion so that we can
    ! modify it - but accept it as it is, so we can ignore that.
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
                              
    ! Now to the edge 2 of boundary component 1 the domain. We use the
    ! same two routines to add the boundary condition to p_rboundaryConditions.
    !
    ! Edge 2 is Neumann boudary
    !CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    !CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
    !                                   rboundaryRegion,p_rbcRegion)
                              
    ! Edge 3 of boundary component 1.
    CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
    
    ! Edge 4 of boundary component 1. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
      
    ! Now continue with defining the boundary conditions of the Y-velocity:
    !
    ! Define edge 1.
    CALL boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,2,&
                                       rboundaryRegion,p_rbcRegion)
     
    ! Define edge 2 - Neumann boundary                         
    !CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    !CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,2,&
    !                                   rboundaryRegion,p_rbcRegion)
                              
    ! Define Edge 3
    CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,2,&
                                       rboundaryRegion,p_rbcRegion)
    
    ! Define Edge 4. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,2,&
                                       rboundaryRegion,p_rbcRegion)

    ! Install these analytic boundary conditions into all discretisation
    ! structures on all levels.
                               
    DO i=rproblem%ilvmin,rproblem%ilvmax
      
      ! Ask the problem structure to give us the discretisation structure...
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      ! and inform the discretisation which analytic boundary conditions to use:
      p_rdiscretisation%p_rboundaryConditions => rproblem%p_rboundaryConditions

    END DO
    
    ! The pressure does not need boundary conditions.
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st1_initDiscreteBC (rproblem)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

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
      
      ! For the discrete problem, we need a discrete version of the above
      ! boundary conditions. So we have to discretise them.
      ! The following routine gives back p_rdiscreteBC, a pointer to a
      ! discrete version of the boundary conditions. Remark that
      ! the pointer has to be nullified before calling the routine,
      ! otherwise, the routine tries to update the boundary conditions
      ! in p_rdiscreteBC!
      ! getBoundaryValues is a callback routine that specifies the
      ! values on the boundary. We pass our collection structure as well
      ! to this routine, so the callback routine has access to everything what is
      ! in the collection.
      NULLIFY(rproblem%RlevelInfo(i)%p_rdiscreteBC)
      CALL bcasm_discretiseBC (p_rdiscretisation,rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                              .FALSE.,getBoundaryValues,rproblem%rcollection)
                               
      ! Hang the pointer into the the matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%p_rdiscreteBC
      
      p_rmatrix%p_rdiscreteBC => p_rdiscreteBC
      
    END DO

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%ilvmax)%p_rdiscreteBC
    
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector
    
    p_rrhs%p_rdiscreteBC => p_rdiscreteBC
    p_rvector%p_rdiscreteBC => p_rdiscreteBC
                
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st1_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

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

  SUBROUTINE st1_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

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
    TYPE(t_matrixBlock), DIMENSION(NNLEV) :: Rmatrices
    
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

    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    CALL linsol_initMultigrid (p_rsolverNode,p_RfilterChain)
    
    ! Now we have to build up the level information for multigrid.
    !
    ! At first, initialise a standard interlevel projection structure. We
    ! can use the same structure for all levels.
    CALL mlprj_initProjectionMat (rprojection,p_rmatrix)
    
    ! Then set up smoothers / coarse grid solver:
    DO i=ilvmin,ilvmax
      
      ! On the coarsest grid, set up a coarse grid solver and no smoother
      ! On finer grids, set up a smoother but no coarse grid solver.
      NULLIFY(p_rpreconditioner)
      NULLIFY(p_rsmoother)
      NULLIFY(p_rcoarseGridSolver)
      IF (i .EQ. ilvmin) THEN
        ! Set up a BiCGStab solver with VANCA preconditioning as coarse grid solver:
        CALL linsol_initVANCA (p_rpreconditioner,1.0_DP,LINSOL_VANCA_2DSPQ1TQ0)
        CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,p_RfilterChain)
        !p_rcoarseGridSolver%ioutputLevel = 2
        
        ! Setting up UMFPACK coarse grid solver would be:
        ! CALL linsol_initUMFPACK4 (p_rcoarseGridSolver)

      ELSE
        ! Set up the VANCA smoother for multigrid with damping parameter 0.7,
        ! 4 smoothing steps:
        CALL linsol_initVANCA (p_rsmoother,1.0_DP,LINSOL_VANCA_2DSPQ1TQ0)
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
    ! on all levels according to that array.
    Rmatrices(ilvmin:ilvmax) = rproblem%RlevelInfo(ilvmin:ilvmax)%rmatrix
    CALL linsol_setMatrices(p_RsolverNode,Rmatrices(ilvmin:ilvmax))
    
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
    
    ! Release the temporary vector
    CALL lsysbl_releaseVector (rtempBlock)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st1_postprocessing (rproblem)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

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
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC

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
    
    rprjDiscretisation = p_rvector%p_rblockDiscretisation
    
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
    NULLIFY(p_rdiscreteBC)
    CALL bcasm_discretiseBC (rprjDiscretisation,p_rdiscreteBC, &
                            .FALSE.,getBoundaryValues,rproblem%rcollection)
                            
    ! Connect the vector to the BC's
    rprjVector%p_rdiscreteBC => p_rdiscreteBC
    
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
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,'gmv/u1.gmv')

    ! Write velocity field
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    
    CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, p_Ddata2)
    
    ! Write pressure
    CALL lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
    CALL ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
    ! Release the auxiliary vector
    CALL lsysbl_releaseVector (rprjVector)
    
    ! Throw away the discrete BC's - not used anymore.
    CALL bcasm_releaseDiscreteBC (p_rdiscreteBC)
    
    ! Release the auxiliary discretisation structure.
    ! We only release the two substructures we manually created before.
    ! The large structure must not be released - it's a copy of 
    ! another one.
    CALL spdiscr_releaseDiscr (rprjDiscretisation%RspatialDiscretisation(1))
    CALL spdiscr_releaseDiscr (rprjDiscretisation%RspatialDiscretisation(2))
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st1_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
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

  SUBROUTINE st1_doneBC (rproblem)
  
!<description>
  ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    DO i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release our discrete version of the boundary conditions
      CALL bcasm_releaseDiscreteBC (rproblem%RlevelInfo(i)%p_rdiscreteBC)

      ! ...and also the corresponding analytic description.
      CALL bcond_doneBC (rproblem%p_rboundaryConditions)
    END DO
    
  END SUBROUTINE


  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st1_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i
  TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

    DO i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Before we remove the block discretisation structure, remember that
      ! we copied the scalar discretisation structure for the X-velocity
      ! to the Y-velocity.
      ! To prevent errors or wrong deallocation, we manually release the
      ! spatial discretisation structures of each of the components.
      p_rDiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation

      ! Remove spatial discretisation structure of the velocity:
      CALL spdiscr_releaseDiscr(p_rdiscretisation%RspatialDiscretisation(1))
      
      ! Don't remove that of the Y-velocity; there is none :)
      !
      ! Remove the discretisation structure of the pressure.
      CALL spdiscr_releaseDiscr(p_rdiscretisation%RspatialDiscretisation(3))
      
      ! Finally remove the block discretisation structure. Don't release
      ! the substructures again.
      CALL spdiscr_releaseBlockDiscr(p_rdiscretisation,.FALSE.)
      
      ! Remove the discretisation from the heap.
      DEALLOCATE(p_rdiscretisation)
    END DO
    
  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE st1_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    ! For compatibility to old F77: an array accepting a set of triangulations
    INTEGER, DIMENSION(SZTRIA,NNLEV) :: TRIAS


    DO i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release the old FEAT 1.x handles.
      ! Get the old triangulation structure of level ilv from the
      ! FEAT2.0 triangulation:
      TRIAS(:,i) = rproblem%RlevelInfo(i)%p_rtriangulation%Itria
      CALL DNMTRI (i,i,TRIAS)
      
      ! then the FEAT 2.0 stuff...
      CALL tria_done (rproblem%RlevelInfo(i)%p_rtriangulation)
    END DO
    
    ! Finally release the domain.
    CALL boundary_release (rproblem%p_rboundary)
    
    ! Don't forget to throw away the old FEAT 1.0 boundary definition!
    CALL DISPAR

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE stokes1
  
  include 'cmem.inc'
  
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
    CALL st1_initParamTriang (NLMIN,NLMAX,rproblem)
    CALL st1_initDiscretisation (rproblem)    
    CALL st1_initMatVec (rproblem)    
    CALL st1_initAnalyticBC (rproblem)   
    CALL st1_initDiscreteBC (rproblem)
    
    ! Implementation of boundary conditions
    CALL st1_implementBC (rproblem)
    
    ! Solve the problem
    CALL st1_solve (rproblem)
    
    ! Postprocessing
    CALL st1_postprocessing (rproblem)
    
    ! Cleanup
    CALL st1_doneMatVec (rproblem)
    CALL st1_doneBC (rproblem)
    CALL st1_doneDiscretisation (rproblem)
    CALL st1_doneParamTriang (rproblem)

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

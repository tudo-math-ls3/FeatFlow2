!##############################################################################
!# ****************************************************************************
!# <name> codire_method5 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a simple 
!# Convection-Diffusion-Reaction problem with constant coefficients 
!# on a simple domain.
!#
!# The routine splits up the tasks of reading the domain, creating 
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines.
!#
!# The routines here behave similar to codire_method3. In difference,
!# a multigrid-solver with ILU(0) smoother and BiCGStab coarse grid solver
!# is used and matrices/vectors are sorted for Cuthill-McKee.
!# </purpose>
!##############################################################################

MODULE codire_method5

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
  USE sortstrategy
  USE coarsegridcorrection
  USE ucd
  
  USE collection
  USE paramlist
    
  USE codire_callback
  
  IMPLICIT NONE
  
  ! Maximum allowed level in this application; must be =9 for 
  ! FEAT 1.x compatibility (still)!
  INTEGER, PARAMETER :: NNLEV = 9

!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! An object specifying the discretisation (trial/test functions,...)
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! A matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form. The matrix will receive the discrete Laplace operator.
    TYPE(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
  
  END TYPE
  
!</typeblock>


!<typeblock description="Application-specific type block for codire problem">

  TYPE t_problem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    INTEGER :: ilvmin
    
    ! Maximum refinement level
    INTEGER :: ilvmax

    ! A solution vector and a RHS vector on the finest level. 
    TYPE(t_vectorBlock) :: rvector,rrhs

    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary

    ! A variable describing the analytic boundary conditions.    
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by NLMAX!
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

  SUBROUTINE cdrm5_initParamTriang (ilvmin,ilvmax,rproblem)
  
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

  SUBROUTINE cdrm5_initDiscretisation (rproblem)
  
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
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector. In this simple problem, we only have one block.
      ALLOCATE(p_rdiscretisation)
      CALL spdiscr_initBlockDiscr2D (p_rdiscretisation,1,&
                                    p_rtriangulation, p_rboundary)

      ! Save the discretisation structure to our local LevelInfo structure
      ! for later use.
      rproblem%RlevelInfo(i)%p_rdiscretisation => p_rdiscretisation

      ! p_rdiscretisation%Rdiscretisations is a list of scalar 
      ! discretisation structures for every component of the solution vector.
      ! Initialise the first element of the list to specify the element
      ! and cubature rule for this solution component:
      CALL spdiscr_initDiscr_simple ( &
                  p_rdiscretisation%RspatialDiscretisation(1), &
                  EL_E011,CUB_G2X2, &
                  p_rtriangulation, p_rboundary)

    END DO
                                   
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cdrm5_initMatVec (rproblem,rparams)
  
!<description>
  ! Calculates the system matrix and RHS vector of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! A parameter list with informations from the DAT file.
  TYPE(t_parlist), INTENT(IN) :: rparams
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

    ! Arrays for the Cuthill McKee renumbering strategy
    INTEGER, DIMENSION(1) :: H_Iresort 
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iresort

    ! Parameters from the DAT file
    REAL(DP) :: alpha11,alpha12,alpha21,alpha22,beta1,beta2,gamma
    CHARACTER(LEN=10) :: Sstr
  
  
    ! And now to the entries of the matrix. For assembling of the entries,
    ! we need a bilinear form, which first has to be set up manually.
    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
    ! scalar system matrix in 2D.
    
    rform%itermCount = 7
    
    ! alpha * Laplace(u)
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X
    
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_DERIV_X
    
    rform%Idescriptors(1,3) = DER_DERIV_X
    rform%Idescriptors(2,3) = DER_DERIV_Y
    
    rform%Idescriptors(1,4) = DER_DERIV_Y
    rform%Idescriptors(2,4) = DER_DERIV_Y
    
    ! (beta1, beta2)^T * grad(u)
    rform%Idescriptors(1,5) = DER_DERIV_X
    rform%Idescriptors(2,5) = DER_FUNC
    
    rform%Idescriptors(1,6) = DER_DERIV_Y
    rform%Idescriptors(2,6) = DER_FUNC
    
    ! gamma * u
    rform%Idescriptors(1,7) = DER_FUNC       
    rform%Idescriptors(2,7) = DER_FUNC

    ! In the standard case, we have constant coefficients:
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    
    ! get the coefficients from the parameter list
    CALL parlst_getvalue_string (rparams, 'EQUATION', 'ALPHA11', Sstr, '1.0')
    READ(Sstr,*) alpha11
    CALL parlst_getvalue_string (rparams, 'EQUATION', 'ALPHA12', Sstr, '0.0')
    READ(Sstr,*) alpha12
    CALL parlst_getvalue_string (rparams, 'EQUATION', 'ALPHA21', Sstr, '0.0')
    READ(Sstr,*) alpha21
    CALL parlst_getvalue_string (rparams, 'EQUATION', 'ALPHA22', Sstr, '1.0')
    READ(Sstr,*) alpha22
    CALL parlst_getvalue_string (rparams, 'EQUATION', 'BETA1', Sstr, '0.0')
    READ(Sstr,*) beta1
    CALL parlst_getvalue_string (rparams, 'EQUATION', 'BETA2', Sstr, '0.0')
    READ(Sstr,*) beta2
    CALL parlst_getvalue_string (rparams, 'EQUATION', 'GAMMA', Sstr, '0.0')
    READ(Sstr,*) gamma
    
    rform%Dcoefficients(1)  = alpha11
    rform%Dcoefficients(2)  = alpha12
    rform%Dcoefficients(3)  = alpha21
    rform%Dcoefficients(4)  = alpha22
    rform%Dcoefficients(5)  = beta1
    rform%Dcoefficients(6)  = beta2
    rform%Dcoefficients(7)  = gamma

    DO i = rproblem%ilvmin, rproblem%ilvmax
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
    
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
    
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)    

      ! Save matrix and vectors to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      CALL collct_setvalue_mat(rproblem%rcollection,'LAPLACE',p_rmatrix,.TRUE.,i)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      CALL bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscretisation(1),LSYSSC_MATRIX9,&
                p_rmatrix%RmatrixBlock(1,1))
    
      ! Update the structural information of the block matrix, as we manually
      ! changed one of the submatrices:
      CALL lsysbl_updateMatStrucInfo (p_rmatrix)
    
      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_CoDiRe for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine, 
      ! so the callback routine has access to everything what is
      ! in the collection.
      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                   p_rmatrix%RmatrixBlock(1,1),coeff_CoDiRe,&
                                   rproblem%rcollection)

      ! Allocate an array for holding the resorting strategy.
      CALL storage_new ('cdrm5_initMatVec', 'Iresort', &
            p_rmatrix%RmatrixBlock(1,1)%NEQ*2, ST_INT, h_Iresort(1), ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(h_Iresort(1),p_Iresort)
      
      ! Calculate the resorting strategy.
      CALL sstrat_calcCuthillMcKee (p_rmatrix%RmatrixBlock(1,1),p_Iresort)
      
      ! Save the handle of the resorting strategy to the collection.
      CALL collct_setvalue_int(rproblem%rcollection,'LAPLACE-CM',h_Iresort(1),.TRUE.,i)
      
      ! Resort the matrix according to the resorting strategy.
      CALL lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.TRUE.,&
                              SSTRAT_CM,h_Iresort(1))
    
    END DO
    
    ! (Only) on the finest level, we need to calculate a RHS vector
    ! and to allocate a solution vector.
    
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector

    ! Save the solution/RHS vector to the collection. Might be used
    ! later (e.g. in nonlinear problems)
    CALL collct_setvalue_vec(rproblem%rcollection,'RHS',p_rrhs,.TRUE.)
    CALL collct_setvalue_vec(rproblem%rcollection,'SOLUTION',p_rvector,.TRUE.)

    ! Now we want to build up the right hand side. At first we need a block
    ! vector of the right structure. Although we could manually create
    ! that vector, the easiest way to set up the vector structure is
    ! to create it by using our matrix as template:
    CALL lsysbl_createVecBlockIndMat (p_rmatrix,p_rrhs, .FALSE.)
    
    ! The vector structure is done but the entries are missing. 
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
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(1),rlinform,.TRUE.,&
              p_rrhs%RvectorBlock(1),coeff_RHS,&
              rproblem%rcollection)
    
    ! Now we have block vectors for the RHS and the matrix. What we
    ! need additionally is a block vector for the solution. 
    ! Create them using the RHS as template.
    ! Fill the solution vector with 0:
    CALL lsysbl_createVecBlockIndirect (p_rrhs, p_rvector, .TRUE.)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cdrm5_initAnalyticBC (rproblem)
  
!<description>
  ! This initialises the analytic bonudary conditions of the problem
  ! and saves them to the problem structure.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables
  INTEGER :: i

    ! A set of variables describing the analytic boundary conditions.    
    TYPE(t_boundaryRegion) :: rboundaryRegion
    TYPE(t_bcRegion), POINTER :: p_rbcRegion
    
    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! A pointer to the domain
    TYPE(t_boundary), POINTER :: p_rboundary
  
    ! Ask the problem structure to give us the discretisation structure and
    p_rdiscretisation => rproblem%RlevelInfo(1)%p_rdiscretisation
    
    ! Get the domain from the discretisation
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
    ! Set p_rboundaryConditions to NULL() to create a new structure on the heap.
    NULLIFY (rproblem%p_rboundaryConditions)
    CALL bcond_initBC (rproblem%p_rboundaryConditions, p_rboundary)
    
    ! We 'know' already (from the problem definition) that we have four boundary
    ! segments in the domain. Each of these, we want to use for inforcing
    ! some kind of boundary condition.
    !
    ! We ask the bondary routines to create a 'boundary region' - which is
    ! simply a part of the boundary corresponding to a boundary segment.
    ! A boundary region roughly contains the type, the min/max parameter value
    ! and whether the endpoints are inside the region or not.
    CALL boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    
    ! We use this boundary region and specify that we want to have Dirichlet
    ! boundary there. The following routine adds a new 'boundary condition region'
    ! for the first segment to the boundary condition structure.
    ! The region will be set up as 'Dirichlet boundary'.
    ! We specify icomponent='1' to indicate that we set up the
    ! Dirichlet BC's for the first (here: one and only) component in the solution
    ! vector.
    ! The routine also returns the created object in p_rbcRegion so that we can
    ! modify it - but accept it as it is, so we can ignore that.
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
                             
    ! Now to the edge 2 of boundary component 1 the domain. We use the
    ! same two routines to add the boundary condition to p_rboundaryConditions.
    CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
                             
    ! Edge 3 of boundary component 1.
    CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
    
    ! Edge 4 of boundary component 1. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
                             
    ! Install these boundary conditions into all discretisation structures
    ! on all levels.
                               
    DO i=rproblem%ilvmin,rproblem%ilvmax
      
      ! Ask the problem structure to give us the discretisation structure and
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      ! inform the discretisation which analytic boundary conditions to use:
      p_rdiscretisation%p_rboundaryConditions => rproblem%p_rboundaryConditions
      
    END DO
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cdrm5_initDiscreteBC (rproblem)
  
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
    
      ! Get our matrix from the problem structure.
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

  SUBROUTINE cdrm5_implementBC (rproblem)
  
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

  SUBROUTINE cdrm5_solve (rproblem)
  
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
    
    ! Resort the RHS and solution vector according to the resorting
    ! strategy given in the matrix.
    IF (lsysbl_isMatrixSorted(p_rmatrix)) THEN
      ! Use the temporary vector from above to store intermediate data.
      ! The vectors are assumed to know how they are resorted (the strategy
      ! is already attached to them). So call the resorting routines
      ! to resort them as necessary!
      ! We use the first subvector of rtempBlock as temporary data; it's
      ! large enough, as we only have one block.
      CALL lsysbl_sortVectorInSitu (p_rrhs,rtempBlock%RvectorBlock(1),.TRUE.)
      CALL lsysbl_sortVectorInSitu (p_rvector,rtempBlock%RvectorBlock(1),.TRUE.)
    END IF
    
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
        ! Set up a BiCGStab solver with ILU preconditioning as coarse grid solver
        ! would be:
        ! CALL linsol_initMILUs1x1 (p_rpreconditioner,0,0.0_DP)
        ! CALL linsol_initBiCGStab (p_rcoarseGridSolver,p_rpreconditioner,p_RfilterChain)
        
        ! Set up UMFPACK coarse grid solver.
        CALL linsol_initUMFPACK4 (p_rcoarseGridSolver)

      ELSE
        ! Setting up Jacobi smoother for multigrid would be:
        ! CALL linsol_initJacobi (p_rsmoother)

        ! Setting up Jin-Wei-Tam smoother for multigrid would be:
        ! CALL linsol_initJinWeiTam (p_rsmoother)

        ! Setting up VANCA smoother for multigrid would be:
        ! CALL linsol_initVANCA (p_rsmoother)

        ! Set up an ILU smoother for multigrid with damping parameter 0.7,
        ! 4 smoothing steps:
         CALL linsol_initMILUs1x1 (p_rsmoother,0,0.0_DP)
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
    
    ! Unsort the vectors again in case they were resorted before calling 
    ! the solver.
    ! We use the first subvector of rtempBlock as temporary data; it's
    ! large enough, as we only have one block.
    CALL lsysbl_sortVectorInSitu (p_rrhs,rtempBlock%RvectorBlock(1),.FALSE.)
    CALL lsysbl_sortVectorInSitu (p_rvector,rtempBlock%RvectorBlock(1),.FALSE.)
    
    ! Release the temporary vector
    CALL lsysbl_releaseVector (rtempBlock)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cdrm5_postprocessing (rproblem)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables
  
    ! We need some more variables for postprocessing
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport

    ! A pointer to the solution vector and to the triangulation.
    TYPE(t_vectorBlock), POINTER :: p_rvector
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! Get the solution vector from the problem structure.
    p_rvector => rproblem%rvector
    
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      p_rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing. 
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,'gmv/u5.gmv')
    
    CALL lsyssc_getbase_double (p_rvector%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cdrm5_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    INTEGER :: ihandle,i

    ! Release matrices and vectors on all levels
    DO i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Delete the matrix
      CALL lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrix)

      ! Delete the variables from the collection.
      CALL collct_deletevalue (rproblem%rcollection,'LAPLACE',i)
      
      ! Release the permutation for sorting matrix/vectors
      ihandle = collct_getvalue_int (rproblem%rcollection,'LAPLACE-CM',i)
      CALL storage_free (ihandle)
      CALL collct_deletevalue (rproblem%rcollection,'LAPLACE-CM',i)
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

  SUBROUTINE cdrm5_doneBC (rproblem)
  
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

  SUBROUTINE cdrm5_doneDiscretisation (rproblem)
  
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

    DO i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Delete the block discretisation together with the associated
      ! scalar spatial discretisations....
      CALL spdiscr_releaseBlockDiscr(&
                   rproblem%RlevelInfo(i)%p_rdiscretisation, .TRUE.)
      
      ! and remove the allocated block discretisation structure from the heap.
      DEALLOCATE(rproblem%RlevelInfo(i)%p_rdiscretisation)
    END DO
    
  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cdrm5_doneParamTriang (rproblem)
  
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

  SUBROUTINE codire5
  
  include 'cmem.inc'
  
!<description>
  ! This is a 'separated' CoDiRe solver for solving a convection-diffusion-
  ! reaction problem. The different tasks of the problem are separated into
  ! subroutines. The problem uses a problem-specific structure for the 
  ! communication: All subroutines add their generated information to the
  ! structure, so that the other subroutines can work with them.
  ! (THis is somehow a cleaner implementation than using a collection!).
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

    ! A paramlist structure with parameters from the dat file
    TYPE(t_parlist) :: rparams

    ! NLMIN receives the minimal level where to discretise for supporting
    ! the solution process.
    ! NLMAX receives the level where we want to solve
    INTEGER :: NLMIN, NLMAX
    
    ! A problem structure for our problem
    TYPE(t_problem), TARGET :: rproblem
    
    INTEGER :: i
    
    ! Ok, let's start. 
    ! We want to solve our CoDiRe problem on level...
    NLMIN = 1
    NLMAX = 7
    
    ! Initialise the collection.
    CALL collct_init (rproblem%rcollection)
    DO i=1,NLMAX
      CALL collct_addlevel_all (rproblem%rcollection)
    END DO
    
    ! Initialise the parameter list
    CALL parlst_init(rparams)
    
    ! Read the parameters from disc and put a reference to it
    ! to the collection
    CALL parlst_readfromfile(rparams, 'data/codire.dat')
    CALL collct_setvalue_parlst (rproblem%rcollection, 'PARAMS', rparams, .TRUE.)

    ! We want to solve our Laplace problem on level...
    CALL parlst_getvalue_int (rparams, 'GENERAL', 'NLMAX', NLMAX, 7)

    ! So now the different steps - one after the other.
    !
    ! Initialisation
    CALL cdrm5_initParamTriang (NLMIN,NLMAX,rproblem)
    CALL cdrm5_initDiscretisation (rproblem)    
    CALL cdrm5_initMatVec (rproblem,rparams)    
    CALL cdrm5_initAnalyticBC (rproblem)   
    CALL cdrm5_initDiscreteBC (rproblem)
    
    ! Implementation of boundary conditions
    CALL cdrm5_implementBC (rproblem)
    
    ! Solve the problem
    CALL cdrm5_solve (rproblem)
    
    ! Postprocessing
    CALL cdrm5_postprocessing (rproblem)
    
    ! Cleanup
    CALL cdrm5_doneMatVec (rproblem)
    CALL cdrm5_doneBC (rproblem)
    CALL cdrm5_doneDiscretisation (rproblem)
    CALL cdrm5_doneParamTriang (rproblem)
    
    ! Release parameter list
    CALL collct_deletevalue (rproblem%rcollection,'PARAMS')
    CALL parlst_done (rparams)

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
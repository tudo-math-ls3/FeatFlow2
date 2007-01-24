!##############################################################################
!# ****************************************************************************
!# <name> heatcond_method5 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a simple 
!# heat conduction problem with constant coefficients 
!# on a simple domain.
!#
!# The routine splits up the tasks of reading the domain, creating 
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines.
!#
!# The routines here behave similar to heatcond_method3. In difference,
!# a multigrid-solver with ILU(0) smoother and BiCGStab coarse grid solver
!# is used and matrices/vectors are sorted for Cuthill-McKee.
!# </purpose>
!##############################################################################

MODULE heatcond_method5

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
  USE timestepping
  USE genoutput
  
  USE collection
  USE paramlist
    
  USE heatcond_callback
  
  IMPLICIT NONE
  
  ! Maximum allowed level in this application; must be =9 for 
  ! FEAT 1.x compatibility (still)!
  INTEGER, PARAMETER :: NNLEV = 9

!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (trial/test functions,...)
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! The static matrix (containing Laplace, convection,...) which does not
    ! change in time.
    TYPE(t_matrixBlock) :: rmatrixStatic
    
    ! The mass matrix
    TYPE(t_matrixBlock) :: rmatrixMass

    ! System matrix. May change during the time iteration    
    TYPE(t_matrixBlock) :: rmatrix
    
    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
  
  END TYPE
  
!</typeblock>

!<typeblock description="Configuration block for the time stepping.">

  TYPE t_problem_nonst
  
    ! Configuration block of the time stepping scheme.
    TYPE(t_explicitTimeStepping)        :: rtimestepping
    
    ! Number of current time step
    INTEGER                             :: iiteration
    
    ! Maximum number of time steps
    INTEGER                             :: niterations
    
    ! Start time
    REAL(DP)                            :: dtimemin
    
    ! Current time
    REAL(DP)                            :: dtime
    
    ! Maximum time
    REAL(DP)                            :: dtimemax
  
  END TYPE

!</typeblock>

!<typeblock description="Application-specific type block for heatcond problem">

  TYPE t_problem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    INTEGER :: ilvmin
    
    ! Maximum refinement level
    INTEGER :: ilvmax

    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary

    ! A variable describing the analytic boundary conditions.    
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions

    ! A RHS vector on the finest level used for solving linear systems
    TYPE(t_vectorBlock) :: rrhs
    
    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    
    ! A filter chain to filter vectors during the solution process
    TYPE(t_filterChain), DIMENSION(1) :: RfilterChain

    ! A parameter block for everything that controls the time dependence.
    TYPE(t_problem_nonst) :: rtimedependence

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

  SUBROUTINE hc5_initParamTriang (ilvmin,ilvmax,rproblem)
  
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
    DO i=rproblem%ilvmin,rproblem%ilvmax
      CALL tria_wrp_tria2Structure(TRIAS(:,i),rproblem%RlevelInfo(i)%rtriangulation)
    END DO
    
    ! The TRIAS(,)-array is now part pf the triangulation structure,
    ! we don't need it anymore.
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
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
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
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

  SUBROUTINE hc5_initMatVec (rproblem,rparams)
  
!<description>
  ! Calculates the matrices of the linear system
  ! by discretising the problem with the default discretisation structure
  ! in the problem structure.
  ! Sets up a solution vector for the linear system, allocates memory
  ! for a RHS vector.
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
    TYPE(t_bilinearForm) :: rform,rformmass
    TYPE(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    TYPE(t_matrixBlock), POINTER :: p_rmatrixStatic,p_rmatrixMass,p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs

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

    ! For the time dependent problem, we also need a mass matrix. We set up another
    ! bilinear form for that.
    
    rformmass%itermCount = 1
    
    ! One and only term for the mass matrix
    rformmass%Idescriptors(1,1) = DER_FUNC
    rformmass%Idescriptors(2,1) = DER_FUNC

    ! The coefficient in front of the term of the mass matrix
    rformmass%Dcoefficients(1)  = 1.0_DP

    DO i = rproblem%ilvmin, rproblem%ilvmax
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
    
      ! -------------------------------------------------------------
      ! Laplace matrix
    
      p_rmatrixStatic => rproblem%RlevelInfo(i)%rmatrixStatic
    
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrixStatic)    

      ! Save matrix and vectors to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      CALL collct_setvalue_mat(rproblem%rcollection,'LAPLACE',p_rmatrixStatic,.TRUE.,i)

      ! Now as the discretisation is set up, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      CALL bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscretisation(1),LSYSSC_MATRIX9,&
                p_rmatrixStatic%RmatrixBlock(1,1))
    
      ! Update the structural information of the block matrix, as we manually
      ! changed one of the submatrices:
      CALL lsysbl_updateMatStrucInfo (p_rmatrixStatic)
    
      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_heatcond for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine, 
      ! so the callback routine has access to everything what is
      ! in the collection.
      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                   p_rmatrixStatic%RmatrixBlock(1,1),coeff_heatcond,&
                                   rproblem%rcollection)

      ! -------------------------------------------------------------
      ! Mass matrix
      
      p_rmatrixMass => rproblem%RlevelInfo(i)%rmatrixMass
      
      ! Initialise the block matrix with default values based on
      ! the discretisation.
      CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrixMass)    

      ! Save matrix and vectors to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      CALL collct_setvalue_mat(rproblem%rcollection,'MASS',p_rmatrixMass,.TRUE.,i)
      
      ! The structure of the mass matrix is the same as the system matrix.
      ! Initialise the structure as "shared" with the system matrix.
      ! Reserve memory for the entries.
      CALL lsyssc_duplicateMatrix(p_rmatrixStatic%RmatrixBlock(1,1),&
           p_rmatrixMass%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
      ! Update the structural information of the block matrix, as we manually
      ! changed one of the submatrices:
      CALL lsysbl_updateMatStrucInfo (p_rmatrixMass)

      ! Now we can build the matrix entries of the mass matrix.
      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                   p_rmatrixMass%RmatrixBlock(1,1))

      ! -------------------------------------------------------------
      ! System matrix.

      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)    

      ! Save matrix and vectors to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      CALL collct_setvalue_mat(rproblem%rcollection,'SYSTEM',p_rmatrix,.TRUE.,i)
      
      ! The structure of the mass matrix is the same as the system matrix.
      ! Initialise the structure as "shared" with the static matrix.
      ! Reserve memory for the entries.
      CALL lsyssc_duplicateMatrix(p_rmatrixStatic%RmatrixBlock(1,1),&
           p_rmatrix%RmatrixBlock(1,1),LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
      ! Update the structural information of the block matrix, as we manually
      ! changed one of the submatrices:
      CALL lsysbl_updateMatStrucInfo (p_rmatrix)

    END DO

    ! (Only) on the finest level, we need to calculate a RHS vector
    ! and to allocate a solution vector.
    
    p_rrhs    => rproblem%rrhs   
    p_rmatrixStatic => rproblem%RlevelInfo(rproblem%ilvmax)%rmatrixStatic

    ! Save the solution/RHS vector to the collection. Might be used
    ! later (e.g. in nonlinear problems)
    CALL collct_setvalue_vec(rproblem%rcollection,'RHS',p_rrhs,.TRUE.)

    ! Reserve memory for all the vectors on the finest level.
    CALL lsysbl_createVecBlockIndMat (p_rmatrixStatic,p_rrhs, .FALSE.)

    ! -------------------------------------------------------------
    ! Matrix resorting
    
    ! Finally, sort the matrices on all levels. We dfo this after the
    ! creation of the vectors, so the vectors stay unsorted!
    DO i = rproblem%ilvmin, rproblem%ilvmax
    
      p_rmatrixStatic => rproblem%RlevelInfo(i)%rmatrixStatic
      p_rmatrixMass => rproblem%RlevelInfo(i)%rmatrixMass
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! Allocate an array for holding the resorting strategy.
      CALL storage_new ('hc5_initMatVec', 'Iresort', &
            p_rmatrixStatic%RmatrixBlock(1,1)%NEQ*2, ST_INT, h_Iresort(1), &
            ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(h_Iresort(1),p_Iresort)
      
      ! Calculate the resorting strategy.
      CALL sstrat_calcCuthillMcKee (p_rmatrixStatic%RmatrixBlock(1,1),p_Iresort)
      
      ! Save the handle of the resorting strategy to the collection.
      CALL collct_setvalue_int(rproblem%rcollection,'LAPLACE-CM',h_Iresort(1),.TRUE.,i)
      
      ! Resort the matrices according to the sorting strategy.
      ! Note that as we share the structure between all matrices, we first
      ! have to sort the 'child' matrices...
      CALL lsyssc_sortMatrix (p_rmatrixMass%RmatrixBlock(1,1),.TRUE.,&
                              SSTRAT_CM,h_Iresort(1))
      CALL lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.TRUE.,&
                              SSTRAT_CM,h_Iresort(1))

      ! ...before sorting the 'parent' matrix!
      CALL lsyssc_sortMatrix (p_rmatrixStatic%RmatrixBlock(1,1),.TRUE.,&
                              SSTRAT_CM,h_Iresort(1))
                              
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_calcRHS (rproblem,rrhs)
  
!<description>
  ! Calculates the RHS vector at the current point in time.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! A pointer to the RHS vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>

  ! local variables
  INTEGER :: i
  
    ! A linear form describing the analytic problem to solve
    TYPE(t_linearForm) :: rlinform
    
    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

    ! Put the current simulation time as parameter "TIME" into the collection.
    ! Also set Dquickaccess (1) to the simulation time for faster access by the
    ! callback routine.
    rproblem%rcollection%Dquickaccess (1) = rproblem%rtimedependence%dtime 
    CALL collct_setvalue_real(rproblem%rcollection,'TIME',&
         rproblem%rtimedependence%dtime,.TRUE.)

    ! The vector structure is done but the entries are missing. 
    ! So the next thing is to calculate the content of that vector.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    p_rdiscretisation => rproblem%RlevelInfo(rproblem%ilvmax)%p_rdiscretisation
    
    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the 
    ! first block.
    !
    ! We pass our collection structure as well to this routine, 
    ! so the callback routine has access to everything what is
    ! in the collection.
    CALL linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscretisation(1),rlinform,.TRUE.,&
              rrhs%RvectorBlock(1),coeff_RHS,&
              rproblem%rcollection)
    
    ! Remove the "TIME"-parameter from the collection again.
    CALL collct_deletevalue (rproblem%rcollection,'TIME')
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_initAnalyticBC (rproblem)
  
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
    !CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    !CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
    !                                   rboundaryRegion,p_rbcRegion)
    
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

  SUBROUTINE hc5_initDiscreteBC (rproblem, bupdate)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<input>
  ! Whether to update existing boundary conditions.
  ! Should be set to .FALSE. on first call and to .TRUE. for every following
  ! call.
  LOGICAL, INTENT(IN) :: bupdate
!</input>

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
  
    ! Put the current simulation time as parameter "TIME" into the collection.
    ! Also set Dquickaccess (1) to the simulation time for faster access by the
    ! callback routine.
    rproblem%rcollection%Dquickaccess (1) = rproblem%rtimedependence%dtime 
    CALL collct_setvalue_real(rproblem%rcollection,'TIME',&
         rproblem%rtimedependence%dtime,.TRUE.)
    
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
      IF (.NOT. bupdate) NULLIFY(rproblem%RlevelInfo(i)%p_rdiscreteBC)
      CALL bcasm_discretiseBC (p_rdiscretisation,rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                              .FALSE.,getBoundaryValues,rproblem%rcollection)
                               
      ! Hang the pointer into the the matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%p_rdiscreteBC
      
      p_rmatrix%p_rdiscreteBC => p_rdiscreteBC
      
    END DO
    
    ! Remove the "TIME"-parameter from the collection again.
    CALL collct_deletevalue (rproblem%rcollection,'TIME')

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%ilvmax)%p_rdiscreteBC
    
    p_rrhs    => rproblem%rrhs   
    p_rrhs%p_rdiscreteBC => p_rdiscreteBC
                
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_implementBC (rproblem,rvector,rrhs,dtimeWeight)
  
!<description>
  ! Implements boundary conditions into the RHS, a given solution vector
  ! and into the system matrices on all levels specified in rproblem.
!</description>

!<input>
  ! Time stepping weight. Standard is 1.0.
  REAL(DP) :: dtimeWeight
!</input>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! A pointer to the solution vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
  
  ! A pointer to the RHS vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>

  ! local variables
  INTEGER :: i,ilvmax
  
    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
    
    ! Pointer to structure for saving discrete BC's:
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
    
    ! Get our the right hand side and solution from the problem structure
    ! on the finest level
    ilvmax = rproblem%ilvmax
    
    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. Call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%ilvmax)%p_rdiscreteBC
    CALL vecfil_discreteBCrhs (rrhs,dtimeWeight,p_rdiscreteBC)
    CALL vecfil_discreteBCsol (rvector,dtimeWeight,p_rdiscreteBC)

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

  SUBROUTINE hc5_initSolver (rproblem)
  
!<description>
  ! Initialises the linear solver according to the problem rproblem.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
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
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain

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
    
    ! During the linear solver, the boundary conditions must
    ! frequently be imposed to the vectors. This is done using
    ! a filter chain. As the linear solver does not work with 
    ! the actual solution vectors but with defect vectors instead,
    ! a filter for implementing the real boundary conditions 
    ! would be wrong.
    ! Therefore, create a filter chain with one filter only,
    ! which implements Dirichlet-conditions into a defect vector.
    rproblem%RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Create a Multigrid-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => rproblem%RfilterChain
    CALL linsol_initMultigrid (p_rsolverNode,p_RfilterChain)
    
    ! Now we have to build up the level information for multigrid.
    !
    ! At first, initialise a standard interlevel projection structure. We
    ! can use the same structure for all levels.
    CALL mlprj_initProjectionDiscr (rprojection,&
         rproblem%RlevelInfo(ilvmax)%p_rdiscretisation)
    
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
    
    ! Save the solver node in the problem structure, finish
    rproblem%p_rsolverNode => p_rsolverNode

    ! Allocate memory, initialise solver structures according to the
    ! linear system we just attached.
    CALL linsol_initStructure (p_rsolverNode,ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_doneSolver (rproblem)
  
!<description>
  ! Releases the solver from the problem structure.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>
 
    ! Release solver data and structure
    CALL linsol_doneData (rproblem%p_rsolverNode)
    CALL linsol_doneStructure (rproblem%p_rsolverNode)
    
    ! Release the solver node and all subnodes attached to it (if at all):
    CALL linsol_releaseSolver (rproblem%p_rsolverNode)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_timestep (rproblem,rvector,rrhs)
  
!<description>
  ! Performs one time step: $t^n -> t^n+1$. 
  ! Assembles system matrix and RHS vector. 
  ! Solves the corresponding time-step equation and returns the solution vector
  ! at the end of the time step.
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! The current solution vector at time $t^n$. Is replaced by the
  ! solution vector at time $t^{n+1}.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector

  ! The RHS vector at time $t^n$. Is replaced by the RHS at time $t^{n+1}$.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
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
    TYPE(t_vectorBlock), POINTER :: p_rrhs
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
    
    ! We have an equation of the type
    !
    !   d/dt u(x,t)  +  N(u(x,t))  =  f(x,t)
    !
    ! Which is discretised in time with a Theta scheme, leading to
    !
    !   $$ u_{n+1} + w_1*N(u_n+1) 
    !      =   u_n + w_2*N(u_n)  +  w_3*f_{n+1}  +  w_4*f_n $$
    !
    ! with k=time step size, u_{n+1} = u(.,t_{n+1}),etc., c.f. timestepping.f90.
    !
    ! The RHS of that equation therefore contains parts of the solution
    ! u_n, of the old RHS f_n and the new RHS f_{n+1}. At first, we make
    ! a weighted copy of the current RHS f_n to the 'global' RHS vector
    ! according to the time stepping scheme.
    
    ilvmin = rproblem%ilvmin
    ilvmax = rproblem%ilvmax
    
    ! Get our right hand side / solution / matrix on the finest
    ! level from the problem structure.
    p_rmatrix => rproblem%RlevelInfo(ilvmax)%rmatrix
    p_rrhs    => rproblem%rrhs 
    
    ! Create a temporary vector we need that for some preparation.
    CALL lsysbl_createVecBlockIndirect (p_rrhs, rtempBlock, .FALSE.)
    
    ! Set up w_2*N(u_n) + w_4*f_n.
    
    CALL lsysbl_vectorLinearComb(rrhs,p_rrhs,&
         rproblem%rtimedependence%rtimestepping%dweightOldRHS,0.0_DP)
    
    ! Synchronise the sorting of the vectors accoring to the system matrix.
    ! We use the first subvector of rtempBlock as temporary data; it's
    ! large enough, as we have only one block.
    CALL lsysbl_synchroniseSortMatVec (p_rmatrix,p_rrhs,rtempBlock%RvectorBlock(1))
    CALL lsysbl_synchroniseSortMatVec (p_rmatrix,rvector,rtempBlock%RvectorBlock(1))
    
    CALL lsysbl_blockMatVec(rproblem%RlevelInfo(ilvmax)%rmatrixStatic,&
         rvector,p_rrhs,&
         rproblem%rtimedependence%rtimestepping%dweightMatrixRHS,&
         rproblem%rtimedependence%rtimestepping%dweightOldRHS)
         
    ! Add u_n -- or, more precisely, M u_n (with M being the mass matrix),
    ! since the discretisation with finite elements requires that.

    CALL lsysbl_blockMatVec(rproblem%RlevelInfo(ilvmax)%rmatrixMass,&
         rvector,p_rrhs,1.0_DP,1.0_DP)
         
    ! Switch to the next point in time.
    rproblem%rtimedependence%dtime = rproblem%rtimedependence%dtime + &
          rproblem%rtimedependence%rtimestepping%dtstep
          
    ! Generate f_n+1 into the rrhs overwriting the previous RHS.
    CALL hc5_calcRHS (rproblem,rrhs)
    
    ! Add w_3*f_{n+1} to the current RHS. If necessary, unsort p_rrhs back before.
    CALL lsysbl_sortVectorInSitu (p_rrhs,rtempBlock%RvectorBlock(1),.FALSE.)
    
    CALL lsysbl_vectorLinearComb(rrhs,p_rrhs,&
         rproblem%rtimedependence%rtimestepping%dweightNewRHS,1.0_DP)

    ! That's it for the RHS vector.
    !
    ! The LHS "u_{n+1} + w_1*N(u_n+1)" results in the system matrix
    ! "M + w_1 N(.)" for the next linear system to solve. Set up that system
    ! on every level of the discretisation.
    
    DO i = ilvmin,ilvmax
      CALL lsyssc_duplicateMatrix (rproblem%RlevelInfo(i)%rmatrixMass%RmatrixBlock(1,1),&
                                   rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                   LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPY)
      CALL lsyssc_matrixLinearComb (rproblem%RlevelInfo(i)%rmatrixStatic%RmatrixBlock(1,1),&
                                    rproblem%rtimedependence%rtimestepping%dweightMatrixLHS,&
                                    rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                    1.0_DP,&
                                    rproblem%RlevelInfo(i)%rmatrix%RmatrixBlock(1,1),&
                                    .FALSE.,.FALSE.,.TRUE.,.TRUE.)
    END DO
    
    ! Discretise the boundary conditions at the new time instant
    CALL hc5_initDiscreteBC (rproblem,.TRUE.)
          
    ! Implement boundary conditions into the RHS vector, the solution vector
    ! and the current system matrices.
    CALL hc5_implementBC (rproblem,rvector,p_rrhs,1.0_DP)
          
    ! Preparation of the linear system completed!
    !
    ! Attach the system matrices to the solver.
    
    p_rsolverNode => rproblem%p_rsolverNode
    
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array.
    Rmatrices(ilvmin:ilvmax) = rproblem%RlevelInfo(ilvmin:ilvmax)%rmatrix
    CALL linsol_setMatrices(p_rsolverNode,Rmatrices(ilvmin:ilvmax))
    
    ! Initialise data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    CALL linsol_initData (p_rsolverNode,ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    
    ! Synchronise p_rrhs with the matrix so it's compatible to the linear system.
    CALL lsysbl_synchroniseSortMatVec (p_rmatrix,p_rrhs,rtempBlock%RvectorBlock(1))

    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b is a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    CALL linsol_solveAdaptively (p_rsolverNode,rvector,p_rrhs,rtempBlock)
    
    ! rvector is now u_n+1.
    !
    ! Release solver data.
    CALL linsol_doneData (p_rsolverNode)
    
    ! Unsort the vectors again in case they were resorted before calling 
    ! the solver.
    ! We use the first subvector of rtempBlock as temporary data; it's
    ! large enough, as we only have one block.
    CALL lsysbl_sortVectorInSitu (p_rrhs,rtempBlock%RvectorBlock(1),.FALSE.)
    CALL lsysbl_sortVectorInSitu (rvector,rtempBlock%RvectorBlock(1),.FALSE.)
    
    ! Release the temporary vector
    CALL lsysbl_releaseVector (rtempBlock)
    
    ! Finally tell the time stepping scheme that we completed the time step.
    CALL timstp_nextSubstep (rproblem%rtimedependence%rtimestepping)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_postprocessing (rproblem,rvector,iiteration,dtime)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<input>
  ! Number of current iteration
  INTEGER, INTENT(IN) :: iiteration
  
  ! Current simulation time
  REAL(DP), INTENT(IN) :: dtime
!</input>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! The current solution vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>

  ! local variables
  
    ! We need some more variables for postprocessing
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport

    ! A pointer to the solution vector and to the triangulation.
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing. 
    ! Start UCD export to GMV file:
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,&
                       'gmv/u.gmv.'//TRIM(sys_si0L(iiteration,5)))
    CALL ucd_setSimulationTime (rexport,dtime)
    
    CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_doneMatVec (rproblem)
  
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

      ! Delete the variables from the collection.
      CALL collct_deletevalue (rproblem%rcollection,'SYSTEM',i)
      CALL collct_deletevalue (rproblem%rcollection,'LAPLACE',i)
      CALL collct_deletevalue (rproblem%rcollection,'MASS',i)

      ! Delete the system matrix
      CALL lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrix)
      
      ! Delete the mass matrix
      CALL lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixMass)

      ! Delete the matrix
      CALL lsysbl_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixStatic)

      ! Release the permutation for sorting matrix/vectors
      ihandle = collct_getvalue_int (rproblem%rcollection,'LAPLACE-CM',i)
      CALL storage_free (ihandle)
      CALL collct_deletevalue (rproblem%rcollection,'LAPLACE-CM',i)
    END DO

    ! Delete solution/RHS vector
    CALL lsysbl_releaseVector (rproblem%rrhs)

    ! Delete the variables from the collection.
    CALL collct_deletevalue (rproblem%rcollection,'RHS')

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_doneBC (rproblem)
  
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

  SUBROUTINE hc5_doneDiscretisation (rproblem)
  
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

  SUBROUTINE hc5_doneParamTriang (rproblem)
  
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
      TRIAS(:,i) = rproblem%RlevelInfo(i)%rtriangulation%Itria
      CALL DNMTRI (i,i,TRIAS)
      
      ! then the FEAT 2.0 stuff...
      CALL tria_done (rproblem%RlevelInfo(i)%rtriangulation)
    END DO
    
    ! Finally release the domain.
    CALL boundary_release (rproblem%p_rboundary)
    
    ! Don't forget to throw away the old FEAT 1.0 boundary definition!
    CALL DISPAR

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_initparameters (rparams,rproblem)
  
!<description>
  ! Reads the DAT file from disc into the parameter list rparams and
  ! initialises basic variables (number of levels, time stepping technique)
  ! in rproblem according to these settings.
!</description>

!<inputoutput>
  ! A parameter list structure accepting the parameters from the DAT file.
  TYPE(t_parlist), INTENT(INOUT) :: rparams

  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: cscheme,niterations
  REAL(DP) :: dtheta,dtstep,dtimemin,dtimemax

    ! Read the parameters from disc and put a reference to it
    ! to the collection
    CALL parlst_readfromfile(rparams, 'data/heatcond.dat')

    ! We want to solve our Laplace problem on level...
    CALL parlst_getvalue_int (rparams, 'GENERAL', 'NLMIN', rproblem%ilvmin, 7)
    CALL parlst_getvalue_int (rparams, 'GENERAL', 'NLMAX', rproblem%ilvmax, 7)
    
    ! Get the parameters for the time stepping scheme from the parameter list
    CALL parlst_getvalue_int (rparams, 'TIMESTEPPING', 'CSCHEME', cscheme, 0)
    CALL parlst_getvalue_int (rparams, 'TIMESTEPPING', 'NITERATIONS', niterations, 1000)
    CALL parlst_getvalue_double (rparams, 'TIMESTEPPING', 'DTHETA', dtheta, 1.0_DP)
    CALL parlst_getvalue_double (rparams, 'TIMESTEPPING', 'DTSTEP', dtstep, 0.1_DP)
    CALL parlst_getvalue_double (rparams, 'TIMESTEPPING', 'DTIMEMIN', dtimemin, 0.0_DP)
    CALL parlst_getvalue_double (rparams, 'TIMESTEPPING', 'DTIMEMAX', dtimemax, 1.0_DP)
    
    ! Initialise the time stepping in the problem structure
    CALL timstp_init (rproblem%rtimedependence%rtimestepping, &
                      cscheme, dtimemin, dtstep, dtheta)
                     
    rproblem%rtimedependence%niterations = niterations
    
    rproblem%rtimedependence%dtimemin = dtimemin
    rproblem%rtimedependence%dtime = dtimemin
    rproblem%rtimedependence%dtimemax = dtimemax

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_timeloop (rproblem,rvector,rrhs)
  
!<description>
  ! Starts the time discretisation. Proceeds in time until the final time
  ! or the maximum number of time steps is reached.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem

  ! The initial solution vector. Is replaced by the final solution vector.
  ! Must be unsorted.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
  
  ! The initial RHS vector. Is replaced by the final RHS vector.
  ! Must be unsorted and without any boundary conditions implemented.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

  INTEGER :: iiteration
  REAL(DP) :: dtime
  
    ! Let's start the timeloop
  
    iiteration = 1
    rproblem%rtimedependence%dtime = rproblem%rtimedependence%dtimemin
    DO WHILE ((iiteration .LE. rproblem%rtimedependence%niterations) .AND. &
              (rproblem%rtimedependence%dtime .LT. rproblem%rtimedependence%dtimemax))
              
      rproblem%rtimedependence%iiteration = iiteration
      
      CALL output_separator(OU_SEP_MINUS)
      CALL output_line ('Time step '//TRIM(sys_siL(iiteration,6))// &
                        '     Time '//TRIM(sys_sdL(rproblem%rtimedependence%dtime,5)))
      CALL output_lbrk ()
              
      ! Proceed to the next time step
      CALL hc5_timestep (rproblem,rvector,rrhs)
      
      ! Postprocessing. Write out the solution.
      CALL hc5_postprocessing (rproblem,rvector,iiteration,&
           rproblem%rtimedependence%dtime)
           
      iiteration = iiteration+1
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE heatcond5
  
  include 'cmem.inc'
  
!<description>
  ! This is a 'separated' heatcond solver for solving a convection-diffusion-
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

    ! A problem structure for our problem
    TYPE(t_problem), TARGET :: rproblem
    
    ! An initial RHS vector and a solution vector
    TYPE(t_vectorBlock) :: rrhs,rvector
    
    INTEGER :: i
    
    ! Initialise the parameter list
    CALL parlst_init(rparams)
    
    ! Initialise the collection.
    CALL collct_init (rproblem%rcollection)
    CALL collct_setvalue_parlst (rproblem%rcollection, 'PARAMS', rparams, .TRUE.)
    
    ! Read in the parameters from the DAT file and initialise the basic
    ! structures with these.
    CALL hc5_initparameters (rparams,rproblem)

    ! Add space for level information in the collection
    DO i=1,rproblem%ilvmax
      CALL collct_addlevel_all (rproblem%rcollection)
    END DO
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    CALL hc5_initParamTriang (rproblem%ilvmin,rproblem%ilvmax,rproblem)
    CALL hc5_initDiscretisation (rproblem)    
    CALL hc5_initMatVec (rproblem,rparams)    
    CALL hc5_initAnalyticBC (rproblem)   

    ! Use the auxiliary RHS vector on the finest level to create an
    ! initial RHS and solution vector, which we pass later to the timeloop.
    CALL lsysbl_createVecBlockIndirect (rproblem%rrhs,rrhs,.FALSE.)
    CALL lsysbl_createVecBlockIndirect (rproblem%rrhs,rvector,.TRUE.)

    ! Calculate the initial RHS, don't incorporate any BC's.
    CALL hc5_calcRHS (rproblem,rrhs)  
    
    ! Discretise the initial boundary conditions
    CALL hc5_initDiscreteBC (rproblem,.FALSE.)
    
    ! Implement them into the initial solution vector, as we have a zero
    ! vector as initial solution.
    rvector%p_rdiscreteBC => rproblem%rrhs%p_rdiscreteBC
    CALL vecfil_discreteBCsol (rvector)
    
    ! Initialise the solver
    CALL hc5_initSolver (rproblem)
    
    ! Call the timeloop to solve the problem
    CALL hc5_timeloop (rproblem,rvector,rrhs)
    
    ! Release the solver, we donÄt need it anymore
    CALL hc5_doneSolver (rproblem)
    
    ! Cleanup
    CALL hc5_doneMatVec (rproblem)
    CALL hc5_doneBC (rproblem)
    CALL hc5_doneDiscretisation (rproblem)
    CALL hc5_doneParamTriang (rproblem)
    
    ! Release parameter list
    CALL collct_deletevalue (rproblem%rcollection,'PARAMS')
    CALL parlst_done (rparams)
    
    ! Release RHS and solution vector
    CALL lsysbl_releaseVector (rvector)
    CALL lsysbl_releaseVector (rrhs)

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
!##############################################################################
!# ****************************************************************************
!# <name> codire_method3 </name>
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
!# On start of the routine, a data file 'data/codire.dat' is read from
!# disc. The parameters in this file configure the problem to solve.
!# </purpose>
!##############################################################################

MODULE codire_method3

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
  
  USE collection
  USE paramlist
    
  USE codire_callback
  
  IMPLICIT NONE
  
!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! An object specifying the discretisation (trial/test functions,...)
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    
    ! A matrix and vector. The vector accepts the RHS of the problem
    ! in scalar form. The matrix will receive the discrete Laplace operator.
    TYPE(t_matrixBlock) :: rmatrix
    TYPE(t_vectorBlock) :: rvector,rrhs

    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
  
  END TYPE
  
!</typeblock>


!<typeblock description="Application-specific type block for poisson problem">

  TYPE t_problem
  
    ! LV receives the level where we want to solve
    INTEGER :: LV

    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary

    ! A variable describing the analytic boundary conditions.    
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. There is currently
    ! only one level supported, identified by LV!
    TYPE(t_problem_lvl), DIMENSION(1) :: RlevelInfo
    
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

  SUBROUTINE pm2_initParamTriang (ilv,rproblem)
  
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
  ! The level up to where we refine the coarse mesh.
  INTEGER, INTENT(IN) :: ilv
!</input>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

  ! local variables
  
    ! For compatibility to old F77: an array accepting a set of triangulations
    INTEGER, DIMENSION(SZTRIA,20) :: TRIAS

    ! Variable for a filename:  
    CHARACTER(LEN=60) :: CFILE

    ! Initialise the level in the problem structure
    rproblem%LV = ilv

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
    ! Refine it to level LV...
    CFILE = './pre/QUAD.tri'
    CALL INMTRI (2,TRIAS,ilv,ilv,0,CFILE)
    
    ! ... and create a FEAT 2.0 triangulation for that. Until the point where
    ! we recreate the triangulation routines, this method has to be used
    ! to get a triangulation.
    ! Set p_rtriangulation to NULL() to create a new structure on the heap.
    NULLIFY(rproblem%RlevelInfo(1)%p_rtriangulation)
    CALL tria_wrp_tria2Structure(TRIAS(:,ilv),rproblem%RlevelInfo(1)%p_rtriangulation)
    
    ! The TRIAS(,)-array is now part pf the triangulation structure,
    ! we don't need it anymore.
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm2_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

  ! local variables
  
    ! An object for saving the domain:
    TYPE(t_boundary), POINTER :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! Ask the problem structure to give us the boundary and triangulation.
    ! We need it for the discretisation.
    p_rboundary => rproblem%p_rboundary
    p_rtriangulation => rproblem%RlevelInfo(1)%p_rtriangulation
    
    ! Now we can start to initialise the discretisation. Set up
    ! a simple discretisation structure suing the boundary and
    ! triangulation information. Specify the element and cubature rule
    ! to use during the assembly of matrices.
    !
    ! Note that we initialise only one discretisation structure here,
    ! as our solution is scalar. Normally, we have to initialise one
    ! discretisation structure for every component of the solution!
    ! Set p_rdiscretisation to NULL() to create a new structure on the heap.
    NULLIFY(rproblem%RlevelInfo(1)%p_rdiscretisation)
    CALL spdiscr_initDiscr_simple (rproblem%RlevelInfo(1)%p_rdiscretisation,&
                                   EL_E011,CUB_TRZ,&
                                   p_rtriangulation, p_rboundary)
                                   
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm2_initMatVec (rproblem,rparams)
  
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
  
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector

    ! A pointer to the discretisation structure with the data.
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    
    ! Parameters from the DAT file
    REAL(DP) :: alpha11,alpha12,alpha21,alpha22,beta1,beta2,gamma
    CHARACTER(LEN=10) :: Sstr
  
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rproblem%RlevelInfo(1)%p_rdiscretisation
    
    p_rmatrix => rproblem%RlevelInfo(1)%rmatrix
    p_rrhs    => rproblem%RlevelInfo(1)%rrhs   
    p_rvector => rproblem%RlevelInfo(1)%rvector
    
    ! Save matrix and vectors to the collection.
    ! They maybe used later, expecially in nonlinear problems.
    CALL collct_setvalue_vec(rproblem%rcollection,'RHS',p_rrhs,.TRUE.)
    CALL collct_setvalue_vec(rproblem%rcollection,'SOLUTION',p_rvector,.TRUE.)
    CALL collct_setvalue_mat(rproblem%rcollection,'LAPLACE',p_rmatrix,.TRUE.)

    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create that directly in the block (1,1) of the block matrix.
    CALL bilf_createMatrixStructure (p_rdiscretisation,LSYSSC_MATRIX9,&
                                     p_rmatrix%RmatrixBlock(1,1))
    
    ! Update the structural information of the block matrix, as we manually
    ! changed one of the submatrices:
    CALL lsysbl_updateMatStrucInfo (p_rmatrix)
    
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
    rform%Idescriptors(1,5) = DER_FUNC
    rform%Idescriptors(2,5) = DER_DERIV_X
    
    rform%Idescriptors(1,6) = DER_FUNC
    rform%Idescriptors(2,6) = DER_DERIV_Y
    
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
    
    ! Now we can build the matrix entries.
    ! We specify the callback function coeff_Laplace for the coefficients.
    ! As long as we use constant coefficients, this routine is not used.
    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
    ! the framework will call the callback routine to get analytical data.
    !
    ! We pass our collection structure as well to this routine, 
    ! so the callback routine has access to everything what is
    ! in the collection.
    CALL bilf_buildMatrixScalar (p_rdiscretisation,rform,.TRUE.,&
                                 p_rmatrix%RmatrixBlock(1,1),coeff_Laplace,&
                                 rproblem%rcollection)
    
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
    ! the block vector.
    !
    ! We pass our collection structure as well to this routine, 
    ! so the callback routine has access to everything what is
    ! in the collection.
    CALL bilf_buildVectorScalar (p_rdiscretisation,rlinform,.TRUE.,&
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

  SUBROUTINE pm2_initAnalyticBC (rproblem)
  
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
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    
    ! A pointer to the domain
    TYPE(t_boundary), POINTER :: p_rboundary
  
    ! Ask the problem structure to give us the discretisation structure and
    p_rdiscretisation => rproblem%RlevelInfo(1)%p_rdiscretisation
    
    ! Get the domain from the discretisation
    p_rboundary => p_rdiscretisation%p_rdomain
    
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
    CALL scbc_initScalarBC (rproblem%p_rboundaryConditions,&
                            p_rdiscretisation%p_rdomain)
    
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
    ! The region will be set up as 'Dirichlet boundary' with parameter values
    ! from the previous boundary region.
    ! The routine also returns the created object in p_rbcRegion so that we can
    ! modify it - but accept it as it is, so we can ignore that.
    CALL scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,rproblem%p_rboundaryConditions,&
                             rboundaryRegion,p_rbcRegion)
                             
    ! Now to the edge 2 of boundary component 1 the domain. We use the
    ! same two routines to add the boundary condition to p_rboundaryConditions.
    CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    CALL scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,rproblem%p_rboundaryConditions,&
                             rboundaryRegion,p_rbcRegion)
                             
    ! Edge 3 of boundary component 1.
    CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    CALL scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,rproblem%p_rboundaryConditions,&
                             rboundaryRegion,p_rbcRegion)
    
    ! Edge 4 of boundary component 1. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,rproblem%p_rboundaryConditions,&
                             rboundaryRegion,p_rbcRegion)
                             
    ! The boundary conditions are set up, but still the discretisation
    ! does not know about it. So inform the discretisation which
    ! analytic boundary conditions to use:
    p_rdiscretisation%p_rboundaryConditions => rproblem%p_rboundaryConditions
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm2_initDiscreteBC (rproblem)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables

    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    
    ! Pointer to structure for saving discrete BC's:
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC

    ! Get our matrix and right hand side from the problem structure.
    p_rrhs    => rproblem%RlevelInfo(1)%rrhs   
    p_rvector => rproblem%RlevelInfo(1)%rvector
    p_rmatrix => rproblem%RlevelInfo(1)%rmatrix
    
    ! From the matrix or the RHS we have access to the discretisation and the
    ! analytic boundary conditions.
    p_rdiscretisation => p_rmatrix%RmatrixBlock(1,1)%p_rspatialDiscretisation
    
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
    NULLIFY(rproblem%RlevelInfo(1)%p_rdiscreteBC)
    CALL bcasm_discretiseBC (p_rdiscretisation,rproblem%RlevelInfo(1)%p_rdiscreteBC, &
                             .FALSE.,getBoundaryValues,rproblem%rcollection)
                             
    ! Hang the pointer into the vectors and the matrix - more precisely,
    ! to the first block matrix and the first subvector. That way, these
    ! boundary conditions are always connected to that matrix and that
    ! vector.
    p_rdiscreteBC => rproblem%RlevelInfo(1)%p_rdiscreteBC
    
    p_rmatrix%RmatrixBlock(1,1)%p_rdiscreteBC => p_rdiscreteBC
    p_rrhs%RvectorBlock(1)%p_rdiscreteBC => p_rdiscreteBC
    p_rvector%RvectorBlock(1)%p_rdiscreteBC => p_rdiscreteBC
                
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm2_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables
  
    ! A filter chain to pre-filter the vectors and the matrix.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain

    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector

    ! Get our matrix and right hand side from the problem structure.
    p_rrhs    => rproblem%RlevelInfo(1)%rrhs   
    p_rvector => rproblem%RlevelInfo(1)%rvector
    p_rmatrix => rproblem%RlevelInfo(1)%rmatrix
    
    ! The next step is to set up a filter that modifies the block
    ! vectors according to boundary conditions.
    ! Initialise the first filter of the filter chain as boundary
    ! implementation filter:
    RfilterChain(1)%ifilterType = FILTER_DISCBCSOLREAL
    
    ! Apply the filter chain to the matrix and the vectors.
    ! As the filter consists only of an implementation filter for
    ! boundary conditions, this implements the boundary conditions
    ! into the vectors and matrices
    CALL filter_applyFilterChainVec (p_rrhs, RfilterChain)
    CALL filter_applyFilterChainVec (p_rvector, RfilterChain)
    CALL filter_applyFilterChainMat (p_rmatrix, RfilterChain)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm2_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a linear solver.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables
  
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
    TYPE(t_linsolNode), POINTER :: p_rsolverNode, p_rpreconditioner

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices
    
    ! Error indicator during initialisation of the solver
    INTEGER :: ierror    

    ! Get our matrix and right hand side from the problem structure.
    p_rrhs    => rproblem%RlevelInfo(1)%rrhs   
    p_rvector => rproblem%RlevelInfo(1)%rvector
    p_rmatrix => rproblem%RlevelInfo(1)%rmatrix
    
    ! Create a temporary vector for the solver - it needs that.
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

    ! Create a BiCGStab-solver. Attach the above filter chain
    ! to the solver, so that the solver automatically filters
    ! the vector during the solution process.
    p_RfilterChain => RfilterChain
    NULLIFY(p_rpreconditioner)
    CALL linsol_initBiCGStab (p_rsolverNode,p_rpreconditioner,p_RfilterChain)
    
    ! Set the output level of the solver to 2 for some output
    p_rsolverNode%ioutputLevel = 2

    ! Attach the system matrix to the solver.
    ! First create an array with the matrix data (on all levels, but we
    ! only have one level here), then call the initialisation 
    ! routine to attach all these matrices.
    ! Remark: Don't make a call like
    !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
    ! This doesn't work on all compilers, since the compiler would have
    ! to create a temp array on the stack - which does not always work!
    Rmatrices = (/p_rmatrix/)
    CALL linsol_setMatrices(p_RsolverNode,Rmatrices)
    
    ! Initialise structure/data of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    CALL linsol_initStructure (p_rsolverNode,ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    CALL linsol_initData (p_rsolverNode,ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
    
    ! Finally solve the system. As we want to solve Ax=b with
    ! b being the real RHS and x being the real solution vector,
    ! we use linsol_solveAdaptively. If b would be a defect
    ! RHS and x a defect update to be added to a solution vector,
    ! we would have to use linsol_precondDefect instead.
    CALL linsol_solveAdaptively (p_rsolverNode,p_rmatrix,&
                                 p_rvector,p_rrhs,rtempBlock)
    
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

  SUBROUTINE pm2_postprocessing (rproblem)
  
!<description>
  ! Writes the solution into a GMV file.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables
  
    ! We need some more variables for postprocessing - i.e. writing
    ! a GMV file.
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    INTEGER :: NCELLS,NVERTS
    INTEGER :: ihandle

    ! A pointer to the solution vector and to the triangulation.
    TYPE(t_vectorBlock), POINTER :: p_rvector
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! Get the solution vector from the problem structure.
    p_rvector => rproblem%RlevelInfo(1)%rvector
    
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      p_rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing. Call the GMV library to write out
    ! a GMV file for our solution.
    ihandle = sys_getFreeUnit()
    CALL GMVOF0 (ihandle,-2,'gmv/u3.gmv')
    CALL GMVHEA (ihandle)
    CALL GMVTRI (ihandle,p_rtriangulation%Itria,0,NCELLS,NVERTS)
    
    CALL storage_getbase_double (p_rvector%RvectorBlock(1)%h_Ddata,p_Ddata)
    CALL GMVSCA (ihandle,p_rtriangulation%Itria,1,NVERTS,&
                 p_rvector%RvectorBlock(1)%NEQ,p_Ddata,'sol')
    
    CALL GMVFOT (ihandle)
    CLOSE(ihandle)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm2_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! Release matrix and vectors
    CALL lsysbl_releaseVector (rproblem%RlevelInfo(1)%rvector)
    CALL lsysbl_releaseVector (rproblem%RlevelInfo(1)%rrhs)
    CALL lsysbl_releaseMatrix (rproblem%RlevelInfo(1)%rmatrix)

    ! Delete the variables from the collection.
    CALL collct_deletevalue (rproblem%rcollection,'RHS')
    CALL collct_deletevalue (rproblem%rcollection,'SOLUTION')
    CALL collct_deletevalue (rproblem%rcollection,'LAPLACE')

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm2_doneBC (rproblem)
  
!<description>
  ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! Release our discrete version of the boundary conditions
    CALL bcasm_releaseDiscreteBC (rproblem%RlevelInfo(1)%p_rdiscreteBC)

    ! ...and also the corresponding analytic description.
    CALL scbc_doneScalarBC (rproblem%p_rboundaryConditions)
    
  END SUBROUTINE


  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm2_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! Delete the discretisation.
    CALL spdiscr_releaseDiscr(rproblem%RlevelInfo(1)%p_rdiscretisation)
    
  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE pm2_doneParamTriang (rproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! For compatibility to old F77: an array accepting a set of triangulations
    INTEGER, DIMENSION(SZTRIA,20) :: TRIAS

    ! Release the old FEAT 1.x handles.
    ! Get the old triangulation structure of level ilv from the
    ! FEAT2.0 triangulation:
    TRIAS(:,rproblem%LV) = rproblem%RlevelInfo(1)%p_rtriangulation%Itria
    CALL DNMTRI (rproblem%LV,rproblem%LV,TRIAS)
    
    ! then the FEAT 2.0 stuff...
    CALL tria_done (rproblem%RlevelInfo(1)%p_rtriangulation)
    
    ! Finally release the domain.
    CALL boundary_release (rproblem%p_rboundary)
    
    ! Don't forget to throw away the old FEAT 1.0 boundary definition!
    CALL DISPAR

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE codire3
  
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

    ! A paramlist structure with parameters from the dat file
    TYPE(t_parlist) :: rparams

    ! LV receives the level where we want to solve
    INTEGER :: LV
    
    ! A problem structure for our problem
    TYPE(t_problem), TARGET :: rproblem
    
    ! A temporary string
    CHARACTER(LEN=10) :: Sstr
    
    ! Ok, let's start. 
    ! Initialise the collection.
    CALL collct_init (rproblem%rcollection)
    
    ! Initialise the parameter list
    CALL parlst_init(rparams)
    
    ! Read the parameters from disc and put a reference to it
    ! to the collection
    CALL parlst_readfromfile(rparams, 'data/codire.dat')
    CALL collct_setvalue_parlst (rproblem%rcollection, 'PARAMS', rparams, .TRUE.)

    ! We want to solve our Laplace problem on level...

    CALL parlst_getvalue_string (rparams, 'GENERAL', 'NLMAX', Sstr, '7')
    READ(Sstr,*) LV
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    CALL pm2_initParamTriang (LV,rproblem)
    CALL pm2_initDiscretisation (rproblem)    
    CALL pm2_initMatVec (rproblem,rparams)    
    CALL pm2_initAnalyticBC (rproblem)   
    CALL pm2_initDiscreteBC (rproblem)
    
    ! Implementation of boundary conditions
    CALL pm2_implementBC (rproblem)
    
    ! Solve the problem
    CALL pm2_solve (rproblem)
    
    ! Postprocessing
    CALL pm2_postprocessing (rproblem)
    
    ! Cleanup
    CALL pm2_doneMatVec (rproblem)
    CALL pm2_doneBC (rproblem)
    CALL pm2_doneDiscretisation (rproblem)
    CALL pm2_doneParamTriang (rproblem)
    
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

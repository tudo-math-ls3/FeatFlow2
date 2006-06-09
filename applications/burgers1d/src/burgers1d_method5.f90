!##############################################################################
!# ****************************************************************************
!# <name> burgers1d_method5 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a 1D Burgers equation
!# simultaneously in space/time. 
!#
!# The 1D burgers equation in space/time is defined by:
!#
!#    $$  u_t  +  u*u_x  -  \nu u_xx  =  0,   u(x,0) = sin(Pi*x)  $$
!#
!# We solve this equation in the space/time domain
!# $(x,t) \in \Omega=[0,1]x[0,1]$. In this example, we use a direct space-time
!# discretisation, i.e. we don't use a separate time-discretisation.
!# Instead, we replace the $t$ variable by $y$-variable of a usual
!# 2D space discretisation, thus resulting in the formula
!#
!#    $$  u_y  +  u*u_x  -  \nu u_xx  =  0,   u(x,0) = sin(Pi*x)  $$
!#
!# For Solving this in the domain $\Omega$, we follow the usual w3ay:
!# The tasks of reading the domain, creating triangulations, discretisation,
!# solving, postprocessing and creanup into different subroutines. 
!# The communication between these subroutines is done using an 
!# application-specific structure saving problem data.
!#
!# As this problem is nonlinear,we need to invoke a nonlinear solver,
!# which uses a couple of application-spcific callback routines.
!# To provide these routines with necessary data, we build up a collection 
!# structure. This is passed through the solver to the callback routines.
!#
!# The preconditioning in this example is done by a UMFPACK4 Gauss 
!# elimination.
!#
!# </purpose>
!##############################################################################

MODULE burgers1d_method5

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
  USE nonlinearsolver
  
  USE matrixio
  
  USE collection
    
  USE burgers1d_callback
  
  IMPLICIT NONE
  
  ! Maximum allowed level in this application; must be =9 for 
  ! FEAT 1.x compatibility (still)!
  INTEGER, PARAMETER :: NLMAX = 9
  
!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation), POINTER :: p_rtriangulation

    ! An object specifying the discretisation (trial/test functions,...)
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the 
    ! discrete Laplace operator.
    TYPE(t_matrixBlock) :: rmatrix

    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
  
  END TYPE
  
!</typeblock>


!<typeblock description="Application-specific type block for burgers1d problem">

  TYPE t_problem
  
    ! Minimum refinement level; = Level i in RlevelInfo
    INTEGER :: ilvmin
    
    ! Maximum refinement level
    INTEGER :: ilvmax

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
    TYPE(t_problem_lvl), DIMENSION(NLMAX) :: RlevelInfo
    
    ! A collection structure with problem-dependent data
    TYPE(t_collection) :: rcollection
    
  END TYPE

!</typeblock>

!</types>
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_initParamTriang (ilvmin,ilvmax,rproblem)
  
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
    INTEGER, DIMENSION(SZTRIA,NLMAX) :: TRIAS

    ! Variable for a filename:  
    CHARACTER(LEN=60) :: CFILE

    ! Initialise the level in the problem structure
    rproblem%ilvmin = ilvmin
    rproblem%ilvmax = ilvmax
    
    ! Store the min/max level in the collection to be used in 
    ! callback routines
    CALL collct_setvalue_int(rproblem%rcollection,'NLMIN',ilvmin,.TRUE.)
    CALL collct_setvalue_int(rproblem%rcollection,'NLMAX',ilvmax,.TRUE.)

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
    ! problem structure ratzher than those in the parameters.
    CFILE = './pre/QUAD.tri'
    CALL INMTRI (2,TRIAS,rproblem%ilvmin,rproblem%ilvmax,0,CFILE)
    
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

  SUBROUTINE b1d5_initDiscretisation (rproblem)
  
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
    
    DO i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rproblem%p_rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%p_rtriangulation
      
      ! Now we can start to initialise the discretisation. Set up
      ! a simple discretisation structure suing the boundary and
      ! triangulation information. Specify the element and cubature rule
      ! to use during the assembly of matrices.
      !
      ! Note that we initialise only one discretisation structure here,
      ! as our solution is scalar. Normally, we have to initialise one
      ! discretisation structure for every component of the solution!
      ! Set p_rdiscretisation to NULL() to create a new structure on the heap.
      NULLIFY(rproblem%RlevelInfo(i)%p_rdiscretisation)
      CALL spdiscr_initDiscr_simple (rproblem%RlevelInfo(i)%p_rdiscretisation, &
                                    EL_E011,CUB_TRZ,&
                                    p_rtriangulation, p_rboundary)
    END DO
                                   
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_initMatVec (rproblem)
  
!<description>
  ! Calculates the RHS-ector, set up the solution vetor.
  ! Set up the structure of the system matrix/matrices of the linear
  ! system.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

    ! local variables
    INTEGER :: i
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector

    ! A pointer to the discretisation structure with the data.
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
  
    ! Arrays for the Cuthill McKee renumbering strategy
    INTEGER, DIMENSION(1) :: H_Iresort 
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iresort
    
    DO i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! Save matrix to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      CALL collct_setvalue_mat(rproblem%rcollection,'SYSTEMMAT',p_rmatrix,.TRUE.,i)

      ! Now using the discretisation, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix.
      CALL bilf_createMatrixStructure (p_rdiscretisation,LSYSSC_MATRIX9,&
                                       p_rmatrix%RmatrixBlock(1,1))
                                      
      ! Update the structural information of the block matrix, as we manually
      ! changed one of the submatrices:
      CALL lsysbl_updateMatStrucInfo (p_rmatrix)
      
      ! Allocate memory for the matrix, don't calculate the entries.
      ! Remember hat we have a nonlinear matrix, which entries must be build
      ! in evey step of the nonlinear iteration!
      CALL bilf_createEmptyMatrixScalar(p_rmatrix%RmatrixBlock(1,1),.FALSE.)
      
      ! Allocate an array for holding the resorting strategy.
      CALL storage_new ('b1d5_initMatVec', 'Iresort', &
            p_rmatrix%RmatrixBlock(1,1)%NEQ*2, ST_INT, h_Iresort(1), ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(h_Iresort(1),p_Iresort)
      
      ! Calculate the resorting strategy.
      CALL sstrat_calcCuthillMcKee (p_rmatrix%RmatrixBlock(1,1),p_Iresort)
      
      ! Save the handle of the resorting strategy to the collection.
      CALL collct_setvalue_int(rproblem%rcollection,'LAPLACE-CM',h_Iresort(1),.TRUE.,i)
      
      ! We don't resort the matrix yet - this is done later when the entries
      ! are assembled.
      
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
    ! the block vector.
    !
    ! We pass our collection structure as well to this routine, 
    ! so the callback routine has access to everything what is
    ! in the collection.
    !
    ! Note that the vector is unsorted after calling this routine!
    CALL linf_buildVectorScalar(p_rdiscretisation,rlinform,.TRUE.,&
                                p_rrhs%RvectorBlock(1),coeff_RHS,&
                                rproblem%rcollection)
                                
    ! Clear the solution vector on the finest level.
    CALL lsysbl_vectorClear(rproblem%rvector)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_initAnalyticBC (rproblem)
  
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
    ! Set p_rboundaryConditions to NULL() to create a new structure on the heap.
    NULLIFY (rproblem%p_rboundaryConditions)
    CALL scbc_initScalarBC (rproblem%p_rboundaryConditions,p_rboundary)
    
    ! We 'know' already (from the problem definition) that we have four boundary
    ! segments in the domain. Each of these, we want to use for inforcing
    ! some kind of boundary condition.
    ! Our boundary conditions are defined on segment 1,2 and 4 f the QUAD mesh.
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
                              
    ! Ege 3 must be set up as Neumann boundary, which is realised as
    ! simple 'do-$nothing'-boundary conditions. So we don't do anything with edge 3!
    !CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    !CALL scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,rproblem%p_rboundaryConditions,&
    !                        rboundaryRegion,p_rbcRegion)
    
    ! Edge 4 of boundary component 1. That's it.
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,rproblem%p_rboundaryConditions,&
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

  SUBROUTINE b1d5_initDiscreteBC (rproblem)
  
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
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation

  ! Pointer to structure for saving discrete BC's:
  TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
    
    DO i=rproblem%ilvmin,rproblem%ilvmax
    
      ! Get our matrix from the problem structure.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
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
      NULLIFY(rproblem%RlevelInfo(i)%p_rdiscreteBC)
      CALL bcasm_discretiseBC (p_rdiscretisation,rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                              .FALSE.,getBoundaryValues,rproblem%rcollection)
                               
      ! Hang the pointer into the vectors and the matrix - more precisely,
      ! to the first block matrix and the first subvector. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%p_rdiscreteBC
      
      p_rmatrix%RmatrixBlock(1,1)%p_rdiscreteBC => p_rdiscreteBC
      
    END DO

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%ilvmax)%p_rdiscreteBC
    
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector
    
    p_rrhs%RvectorBlock(1)%p_rdiscreteBC => p_rdiscreteBC
    p_rvector%RvectorBlock(1)%p_rdiscreteBC => p_rdiscreteBC
                
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

  ! local variables
  INTEGER :: i,ilvmax
  
    ! A filter chain to pre-filter the vectors and the matrix.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain

    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
    REAL(DP), DIMENSION(:), POINTER :: p_data
    
    ! Set up a filter that modifies the block vectors/matrix
    ! according to boundary conditions.
    ! Initialise the first filter of the filter chain as boundary
    ! implementation filter:
    RfilterChain(1)%ifilterType = FILTER_DISCBCSOLREAL
    
    ! Get our the right hand side and solution from the problem structure
    ! on the finest level
    ilvmax = rproblem%ilvmax
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector
    
    ! Apply the filter chain to the vectors.
    ! As the filter consists only of an implementation filter for
    ! boundary conditions, this implements the boundary conditions
    ! into the vectors and matrices.
    CALL filter_applyFilterChainVec (p_rrhs, RfilterChain)
    CALL filter_applyFilterChainVec (p_rvector, RfilterChain)

  END SUBROUTINE

  ! ***************************************************************************
    SUBROUTINE b1d5_getDefect (ite,rx,rb,rd,p_rcollection)
  
    USE linearsystemblock
    USE collection
    
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
    INTEGER, INTENT(IN)                           :: ite

    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN),TARGET        :: rx

    ! Right hand side vector of the equation.
    TYPE(t_vectorBlock), INTENT(IN)               :: rb
  !</input>
               
  !<inputoutput>
    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    TYPE(t_collection), POINTER                   :: p_rcollection

    ! Defect vector b-A(x)x. This must be filled by the callback routine
    ! with data.
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rd
  !</inputoutput>

    ! local variables
    TYPE(t_bilinearForm) :: rform
    INTEGER :: nlmin,nlmax,i,j
    REAL(DP) :: dnrm
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    REAL(DP), DIMENSION(:), POINTER :: p_data
    INTEGER(I32), DIMENSION(:), POINTER :: p_Kld

    ! A filter chain to pre-filter the vectors and the matrix.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain

      ! Get maximum level from the collection
      nlmax = collct_getvalue_int (p_rcollection,'NLMAX')

      ! Get the system matrix on the maximum level
      p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',nlmax)
      
      ! Put a reference to rx into the collection. This way, we inform the callback
      ! routine of the matrix assembly about the solution vector to use
      ! fot the nonlinear term.
      CALL collct_setvalue_vec(p_rcollection,'RX',rx,.TRUE.)
      
      ! Build the entries with the discretisation routine.
      
      ! For assembling of the entries,
      ! we need a bilinear form, which first has to be set up manually.
      ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
      ! scalar system matrix in 2D.
      
      rform%itermCount = 4
      rform%Idescriptors(1,1) = DER_DERIV_Y   ! u_t
      rform%Idescriptors(2,1) = DER_FUNC
      
      rform%Idescriptors(1,2) = DER_DERIV_X   ! u_x
      rform%Idescriptors(2,2) = DER_FUNC
      
      rform%Idescriptors(1,3) = DER_DERIV_X   ! -u_xx  -> u_x phi_x
      rform%Idescriptors(2,3) = DER_DERIV_X

      ! The 4th last term u_yy is actually not needed. By setting the coefficient
      ! in front of thisterm to 0 (see matrix assembly callback routine), the term
      ! can be switched off. Nevertheless, we add it here for having the
      ! possibility to use it (by setting the cofficient to a value not equal
      ! to 0), which serves as stabilisation for the problem!

      rform%Idescriptors(1,4) = DER_DERIV_Y   ! -u_xx  -> u_x phi_x
      rform%Idescriptors(2,4) = DER_DERIV_Y

      ! In the standard case, we have constant coefficients.
      ! Theoretically, there are some coefficients constant - but for
      ! simplicity, we define them all in the callback routine of the matrix
      ! assembly.
      rform%ballCoeffConstant = .FALSE.
      rform%BconstantCoeff = .FALSE.

      ! Now we can build the matrix entries.
      ! We specify the callback function coeff_Laplace for the coefficients.
      ! As long as we use constant coefficients, this routine is not used.
      ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
      ! the framework will call the callback routine to get analytical data.
      !
      ! We pass our collection structure as well to this routine, 
      ! so the callback routine has access to everything what is
      ! in the collection.
      CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                   p_rmatrix%RmatrixBlock(1,1),coeff_burgers,&
                                   p_rcollection)

      ! Remove the vector from the collection - not necessary anymore.
      CALL collct_deletevalue (p_rcollection,'RX')
      
      ! Set up a filter that modifies the block vectors/matrix
      ! according to boundary conditions.
      ! Initialise the first filter of the filter chain as boundary
      ! implementation filter for defect vectors:
      RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
      
      ! Apply the filter chain to the matrix, too.
      ! This modifies all matrices according to the discrete boundary 
      ! conditions.
      CALL filter_applyFilterChainMat (p_rmatrix, RfilterChain)
      
      ! Build the defect: d=b-Ax
      CALL lsysbl_vectorCopy (rb,rd)
      CALL lsysbl_blockMatVec (p_rmatrix, rx, rd, -1.0_DP, 1.0_DP)
    
      ! Apply the filter chain to the defect vector.
      ! As the filter consists only of an implementation filter for
      ! boundary conditions, this implements the boundary conditions
      ! into the defect vector.
      CALL filter_applyFilterChainVec (rd, RfilterChain)

      ! That's it
  
    END SUBROUTINE
    
  ! ***************************************************************************

    SUBROUTINE b1d5_precondDefect (rd,domega,p_rcollection)
  
    USE linearsystemblock
    USE collection
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration 
    ! vector rx and the right hand side vector rb, this routine has to compute the 
    ! defect vector rd. The routine accepts a pointer to a collection structure 
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<inputoutput>
    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    TYPE(t_collection), POINTER                   :: p_rcollection
    
    ! Damping parameter. Is set to rsolverNode%domega (usually = 1.0_DP)
    ! on call to the callback routine.
    ! The callback routine can modify this parameter according to any suitable
    ! algorithm to calculate an 'optimal damping' parameter. The nonlinear loop
    ! will then use this for adding rd to the solution vector:
    ! $$ x_{n+1} = x_n + domega*rd $$
    REAL(DP), INTENT(INOUT)                       :: domega
  !</inputoutput>
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    INTEGER :: ierror,nlmax
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
  
      ! Get maximum level from the collection
      nlmax = collct_getvalue_int (p_rcollection,'NLMAX')

      ! Get the system matrix on the maximum level
      p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',nlmax)

      ! Our 'parent' (the caller of the nonlinear solver) has prepared
      ! a preconditioner node for us (a linear solver with symbolically
      ! factorised matrices). Get this from the collection.
      
      p_rsolverNode => collct_getvalue_linsol(p_rcollection,'LINSOLVER')

      ! Initialise data of the solver. This in fact performs a numeric
      ! factorisation of the matrices in UMFPACK-like solvers.
      CALL linsol_initData (p_rsolverNode, ierror)
      IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
      
      ! Finally solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
      CALL linsol_precondDefect (p_rsolverNode,rd)

      ! Release the numeric factorisation of the matrix.
      ! We don't release the symbolic factorisation, as we can use them
      ! for the next iteration.
      CALL linsol_doneData (p_rsolverNode)

    END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a nonlinear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

    ! local variables
    TYPE(t_vectorBlock), TARGET :: rtempBlock
    
    TYPE(t_nlsolNode) :: rnlSol
    INTEGER :: ierror,nlmax
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    TYPE(t_matrixBlock), DIMENSION(1) :: Rmatrices
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    
    ! For solving the problem, we need to invoke a nonlinear solver.
    ! This nonlinear solver 
    !  - builds in every step a linearisation of the system matrix
    !  - calls a linear solver for preconditioning
    ! We can save some time in this situation, if we prepare the
    ! linear solver in-advance. This means:
    !  - allocate memory in advance and
    !  - perform a symbolic factorisation of the matrix in advance.
    ! The point is that the entries of the matrix change in every
    ! iteration - but not the matrix structure! So this is something
    ! we can prepare once and use it through the whole solution process!
    !
    ! At first, set up the linear solver as usual:
    CALL linsol_initUMFPACK4 (p_rsolverNode)

    ! Get the system matrix on the finest level...
    p_rmatrix => rproblem%RlevelInfo(rproblem%ilvmax)%rmatrix

    ! And associate it to the solver
    Rmatrices = (/p_rmatrix/)
    CALL linsol_setMatrices(p_rsolverNode,Rmatrices)

    ! Initialise structure of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    ! In fact, solvers like UMFPACK use this for a symbolic factorisation
    ! of the matrix.
    CALL linsol_initStructure (p_rsolverNode, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP
  
    ! Put the prepared solver node to the collection for later use.
    CALL collct_setvalue_linsol(rproblem%rcollection,'LINSOLVER',p_rsolverNode,.TRUE.)
    
    ! Create a temporary vector we need for the nonlinera iteration.
    CALL lsysbl_createVecBlockIndirect (rproblem%rrhs, rtempBlock, .FALSE.)

    ! The nonlinear solver structure rnlSol is initialised by the default
    ! initialisation with all necessary information to solve the problem.
    ! We call the nonlinear solver directly. For preconditioning
    ! and defect calculation, we use our own callback routine.
    !rnlSol%domega = 0.25_DP
    rnlSol%ioutputLevel = 2
    CALL nlsol_performSolve(rnlSol,rproblem%rvector,rproblem%rrhs,rtempBlock,&
                            b1d5_getDefect,b1d5_precondDefect,&
                            rcollection=rproblem%rcollection)

    ! Release the temporary vector
    CALL lsysbl_releaseVector (rtempBlock)
    
    ! Remove the solver node from the collection - not needed anymore there
    CALL collct_deletevalue(rproblem%rcollection,'LINSOLVER')
    
    ! Clean up the linear solver, release all memory, remove the solver node
    ! from memory.
    CALL linsol_releaseSolver (p_rsolverNode)
    
    PRINT *,'Nonlinear solver statistics'
    PRINT *,'---------------------------'
    PRINT *,'Intial defect: ',rnlSol%DinitialDefect(1)
    PRINT *,'Final defect:  ',rnlSol%DfinalDefect(1)
    PRINT *,'#Iterations:   ',rnlSol%iiterations

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_postprocessing (rproblem)
  
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
    p_rvector => rproblem%rvector
    
    ! From the attached discretisation, get the underlying triangulation
    p_rtriangulation => &
      p_rvector%RvectorBlock(1)%p_rspatialDiscretisation%p_rtriangulation
    
    ! p_rvector now contains our solution. We can now
    ! start the postprocessing. Call the GMV library to write out
    ! a GMV file for our solution.
    ihandle = sys_getFreeUnit()
    CALL GMVOF0 (ihandle,-2,'gmv/u5.gmv')
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

  SUBROUTINE b1d5_doneMatVec (rproblem)
  
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
      CALL collct_deletevalue (rproblem%rcollection,'SYSTEMMAT',i)
      
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

  SUBROUTINE b1d5_doneBC (rproblem)
  
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
      CALL scbc_doneScalarBC (rproblem%p_rboundaryConditions)
    END DO
    
  END SUBROUTINE


  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_doneDiscretisation (rproblem)
  
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
      ! Delete the discretisation.
      CALL spdiscr_releaseDiscr(rproblem%RlevelInfo(i)%p_rdiscretisation)
    END DO
    
  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d5_doneParamTriang (rproblem)
  
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
    INTEGER, DIMENSION(SZTRIA,NLMAX) :: TRIAS


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
    
    CALL collct_deleteValue(rproblem%rcollection,'NLMAX')
    CALL collct_deleteValue(rproblem%rcollection,'NLMIN')

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE burgers1d5
  
  include 'cmem.inc'
  
!<description>
  ! This is a 'separated' burgers1d solver for solving a Burgers-1D
  ! problem. The different tasks of the problem are separated into
  ! subroutines. The problem uses a problem-specific structure for the 
  ! communication: All subroutines add their generated information to the
  ! structure, so that the other subroutines can work with them.
  ! For the communication to callback routines of black-box subroutines
  ! (matrix-assembly, preconditioning), a collection is used.
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

    ! LV receives the level where we want to solve
    INTEGER :: LV
    
    ! A problem structure for our problem
    TYPE(t_problem), TARGET :: rproblem
    
    INTEGER :: i
    
    ! Ok, let's start. 
    ! We want to solve our Laplace problem on level...

    LV = 7
    
    ! Initialise the collection
    CALL collct_init (rproblem%rcollection)
    DO i=1,NLMAX
      CALL collct_addlevel_all (rproblem%rcollection)
    END DO

    ! So now the different steps - one after the other.
    !
    ! Initialisation. Pass LV as minimum and maximum level, so
    ! we calculate only on the finest mesh.
    CALL b1d5_initParamTriang (LV,LV,rproblem)
    CALL b1d5_initDiscretisation (rproblem)    
    CALL b1d5_initMatVec (rproblem)    
    CALL b1d5_initAnalyticBC (rproblem)   
    CALL b1d5_initDiscreteBC (rproblem)
    
    ! Implementation of boundary conditions
    CALL b1d5_implementBC (rproblem)
    
    ! Solve the problem
    CALL b1d5_solve (rproblem)
    
    ! Postprocessing
    CALL b1d5_postprocessing (rproblem)
    
    ! Cleanup
    CALL b1d5_doneMatVec (rproblem)
    CALL b1d5_doneBC (rproblem)
    CALL b1d5_doneDiscretisation (rproblem)
    CALL b1d5_doneParamTriang (rproblem)

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

!##############################################################################
!# ****************************************************************************
!# <name> burgers1d_method6 </name>
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
!# As this problem is nonlinear, we need to invoke a nonlinear solver,
!# which uses a couple of application-spcific callback routines.
!# To provide these routines with necessary data, we build up a collection 
!# structure. This is passed through the solver to the callback routines.
!#
!# The preconditioning in this example is done by a Multigrid solver,
!# which performs only a few number of steps, so does not iterate till
!# convergence.
!#
!# </purpose>
!##############################################################################

MODULE burgers1d_method6

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
  USE ucd
  
  USE matrixio
  
  USE collection
    
  USE burgers1d_callback
  
  IMPLICIT NONE
  
!<types>

!<typeblock description="Type block defining all information about one level">

  TYPE t_problem_lvl
  
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation (trial/test functions,...)
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! A system matrix for that specific level. The matrix will receive the 
    ! discrete Laplace operator.
    TYPE(t_matrixBlock) :: rmatrix

    ! A temporary vector for building the solution when assembling the
    ! matrix on lower levels.
    TYPE(t_vectorBlock) :: rtempVector

    ! A variable describing the discrete boundary conditions.    
    TYPE(t_discreteBC) :: rdiscreteBC
  
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

    ! A solver node that accepts parameters for the linear solver    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode

    ! An array of t_problem_lvl structures, each corresponding
    ! to one level of the discretisation. 
    TYPE(t_problem_lvl), DIMENSION(:), POINTER :: RlevelInfo
    
    ! A collection structure with problem-dependent data
    TYPE(t_collection) :: rcollection
    
  END TYPE

!</typeblock>

!</types>
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d6_initParamTriang (ilvmin,ilvmax,rproblem)
  
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
    
    ! Store the min/max level in the collection to be used in 
    ! callback routines
    CALL collct_setvalue_int(rproblem%rcollection,'NLMIN',ilvmin,.TRUE.)
    CALL collct_setvalue_int(rproblem%rcollection,'NLMAX',ilvmax,.TRUE.)

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

  SUBROUTINE b1d6_initDiscretisation (rproblem)
  
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
    
    ! An object for the spatial discretisation
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

    DO i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rproblem%p_rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      
      ! Now we can startto initialise the discretisation. At first, set up
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

  SUBROUTINE b1d6_initMatVec (rproblem)
  
!<description>
  ! Calculates the RHS-ector, set up the solution vetor.
  ! Set up the structure of the system matrix/matrices of the linear
  ! system.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: i
    
    ! A bilinear and linear form describing the analytic problem to solve
    TYPE(t_linearForm) :: rlinform
    
    ! A pointer to the system matrix and the RHS/solution vectors.
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector,p_rtempVector

    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
  
    ! Arrays for the Cuthill McKee renumbering strategy
    INTEGER, DIMENSION(1) :: H_Iresort 
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iresort
    
    DO i=rproblem%ilvmin,rproblem%ilvmax
      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      ! Get the matrix structure; we want to build a template matrix
      ! on the level, which receives the entries later.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix

      ! Initialise the block matrix with default values based on
      ! the discretisation.
      CALL lsysbl_createMatBlockByDiscr (p_rdiscretisation,p_rmatrix)    
      
      ! Save matrix to the collection.
      ! They maybe used later, expecially in nonlinear problems.
      CALL collct_setvalue_mat(rproblem%rcollection,'SYSTEMMAT',p_rmatrix,.TRUE.,i)

      ! Now using the discretisation, we can start to generate
      ! the structure of the system matrix which is to solve.
      ! We create that directly in the block (1,1) of the block matrix
      ! using the discretisation structure of the first block.
      CALL bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscretisation(1),LSYSSC_MATRIX9,&
                p_rmatrix%RmatrixBlock(1,1))

      ! Update the structural information of the block matrix, as we manually
      ! changed one of the submatrices:
      CALL lsysbl_updateMatStrucInfo (p_rmatrix)
      
      ! Allocate memory for the matrix, don't calculate the entries.
      ! Remember hat we have a nonlinear matrix, which entries must be build
      ! in evey step of the nonlinear iteration!
      ! We fill the matrix with 1. This is necessary, as the UMFPACK solver
      ! needs nonzero matrix entries for the symbolic factorisation!
      CALL lsyssc_allocEmptyMatrix(p_rmatrix%RmatrixBlock(1,1),LSYSSC_SETM_ONE)
      
      ! Allocate an array for holding the resorting strategy.
      CALL storage_new ('b1d6_initMatVec', 'Iresort', &
            p_rmatrix%RmatrixBlock(1,1)%NEQ*2, ST_INT, h_Iresort(1), ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(h_Iresort(1),p_Iresort)
      
      ! Calculate the resorting strategy.
      CALL sstrat_calcCuthillMcKee (p_rmatrix%RmatrixBlock(1,1),p_Iresort)
      
      ! Save the handle of the resorting strategy to the collection.
      CALL collct_setvalue_int(rproblem%rcollection,'LAPLACE-CM',h_Iresort(1),.TRUE.,i)
      
      ! Attach the sorting strategy without actually resorting the matrix - 
      ! this is done later when the entries are assembled!
      ! The sorting strategy then holds for the matrix as well as for all
      ! vectors derived from this matrix.
      ! Commenting this line out would lead to totally disabling the
      ! sorting in the whole application...
      CALL lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.TRUE.,&
                              -SSTRAT_CM,h_Iresort(1))

      ! Now on all levels except for the maximum one, create a temporary 
      ! vector on that level, based on the matrix template.
      ! It's used for building the matrices on lower levels.
      IF (i .LT. rproblem%ilvmax) THEN
        p_rtempVector => rproblem%RlevelInfo(i)%rtempVector
        CALL lsysbl_createVecBlockIndMat (p_rmatrix,p_rtempVector,.FALSE.)
        
        ! Add the temp vector to the collection on level i
        ! for use in the callback routine
        CALL collct_setvalue_vec(rproblem%rcollection,'RTEMPVEC',p_rtempVector,&
                                .TRUE.,i)
      END IF
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
                                
    ! Clear the solution vector on the finest level.
    CALL lsysbl_clearVector(rproblem%rvector)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d6_initDiscreteBC (rproblem)
  
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
    TYPE(t_boundaryRegion) :: rboundaryRegion

    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

    ! Pointer to structure for saving discrete BC's:
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
      
    ! A pointer to the domain
    TYPE(t_boundary), POINTER :: p_rboundary

    ! Get the domain from the problem structure
    p_rboundary => rproblem%p_rboundary

    DO i=rproblem%ilvmin,rproblem%ilvmax
    
      ! Get our matrix from the problem structure.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => p_rmatrix%p_rblockDiscretisation
      
      ! Now we have the raw problem. What is missing is the definition of the boudary
      ! conditions.
      ! For implementing boundary conditions, we use a 'filter technique with
      ! discretised boundary conditions'. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
      CALL bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%rdiscreteBC)
      !
      ! We 'know' already (from the problem definition) that we have four boundary
      ! segments in the domain. Each of these, we want to use for enforcing
      ! some kind of boundary condition.
      !
      ! We ask the bondary routines to create a 'boundary region' - which is
      ! simply a part of the boundary corresponding to a boundary segment.
      ! A boundary region roughly contains the type, the min/max parameter value
      ! and whether the endpoints are inside the region or not.
      CALL boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
      
      ! We use this boundary region and specify that we want to have Dirichlet
      ! boundary there. The following call does the following:
      ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
      !   We specify icomponent='1' to indicate that we set up the
      !   Dirichlet BC's for the first (here: one and only) component in the 
      !   solution vector.
      ! - Discretise the boundary condition so that the BC's can be applied
      !   to matrices and vectors
      ! - Add the calculated discrete BC's to rdiscreteBC for later use.
      CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues,rproblem%rcollection)
                               
      ! Now to the edge 2 of boundary component 1 the domain. We use the
      ! same two routines to add the boundary condition to p_rboundaryConditions.
      CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
      CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
         rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
         getBoundaryValues,rproblem%rcollection)
                               
      ! Edge 3 of boundary component 1.
      ! Ege 3 must be set up as Neumann boundary, which is realised as
      ! simple 'do-$nothing'-boundary conditions. So we don't do anything with edge 3
      ! CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
      ! CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
      !     rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
      !     getBoundaryValues,rproblem%rcollection)
      
      ! Edge 4 of boundary component 1. That's it.
      CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
      CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
          rboundaryRegion,rproblem%RlevelInfo(i)%rdiscreteBC,&
          getBoundaryValues,rproblem%rcollection)
                               
      ! Hang the pointer into the vectors and the matrix - more precisely,
      ! to the first block matrix and the first subvector. That way, these
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

  SUBROUTINE b1d6_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: ilvmax
  
    ! A pointer to the RHS/solution vector
    TYPE(t_vectorBlock), POINTER :: p_rrhs,p_rvector
    
    ! Get our the right hand side and solution from the problem structure
    ! on the finest level
    ilvmax = rproblem%ilvmax
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector
    
    ! Next step is to implement boundary conditions into the RHS and
    ! solution. This is done using a vector filter for discrete boundary 
    ! conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors. Call the appropriate vector filter that modifies the vectors
    ! according to the boundary conditions.
    CALL vecfil_discreteBCrhs (p_rrhs)
    CALL vecfil_discreteBCsol (p_rvector)

  END SUBROUTINE

  ! ***************************************************************************
    SUBROUTINE b1d6_getDefect (ite,rx,rb,rd,p_rcollection)
  
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
      INTEGER :: ilvmax
      TYPE(t_matrixBlock), POINTER :: p_rmatrix

      ! Get maximum level from the collection
      ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')

      ! Get the system matrix on the maximum level
      p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',ilvmax)
      
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
      ! We specify the callback function coeff_burgers for the coefficients.
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
      
      ! Implement discrete boundary conditions into the matrix. 
      ! Call the appropriate matrix filter to modify the system matrix 
      ! according to the attached discrete boundary conditions.
      CALL matfil_discreteBC (p_rmatrix)
      
      ! Build the defect: d=b-Ax
      CALL lsysbl_copyVector (rb,rd)
      CALL lsysbl_blockMatVec (p_rmatrix, rx, rd, -1.0_DP, 1.0_DP)
    
      ! Apply the defect-vector filter for discrete boundary conditions
      ! to modify the defect vector according to the (discrete) boundary
      ! conditions.
      CALL vecfil_discreteBCdef (rd)

      ! That's it
  
    END SUBROUTINE
    
  ! ***************************************************************************

    SUBROUTINE b1d6_precondDefect (ite,rd,rx,rb,domega,bsuccess,p_rcollection)
  
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
    ! Number of current iteration. 
    INTEGER, INTENT(IN)                           :: ite

    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    TYPE(t_vectorBlock), INTENT(INOUT)            :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    TYPE(t_collection), POINTER                   :: p_rcollection
    
    ! Damping parameter. Is set to rsolverNode%domega (usually = 1.0_DP)
    ! on the first call to the callback routine.
    ! The callback routine can modify this parameter according to any suitable
    ! algorithm to calculate an 'optimal damping' parameter. The nonlinear loop
    ! will then use this for adding rd to the solution vector:
    ! $$ x_{n+1} = x_n + domega*rd $$
    ! domega will stay at this value until it's changed again.
    REAL(DP), INTENT(INOUT)                       :: domega
  
    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    LOGICAL, INTENT(INOUT)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Current iteration vector
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rx

    ! Current right hand side of the nonlinear system
    TYPE(t_vectorBlock), INTENT(IN), TARGET       :: rb
  !</input>
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    INTEGER :: ierror,ilvmin,ilvmax,i
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    TYPE(t_vectorBlock), POINTER :: p_rvectorFine,p_rvectorCoarse
    TYPE(t_interlevelProjectionBlock), POINTER :: p_rprojection
    TYPE(t_vectorScalar), POINTER :: p_rvectorTemp
    TYPE(t_bilinearForm) :: rform

      ! Our 'parent' (the caller of the nonlinear solver) has prepared
      ! a preconditioner node for us (a linear solver with symbolically
      ! factorised matrices). Get this from the collection.
      p_rsolverNode => collct_getvalue_linsol(p_rcollection,'LINSOLVER')

      ! The matrices are all attached to the solver node, and
      ! symbolic factorisation was already performed. The only thing that
      ! is missing is the numerical factorisation... no wonder,
      ! the matrix entries do not yet exist!
      !
      ! So the task is now: Calculate the matrix entries and perform
      ! a numerical factorisation / initData in the solver!
      !
      ! The matrix on the maximum level is already prepared.
      ! This was done for calculating the nonlinear defect. The matrices
      ! on the other levels are still missing.
      !
      ! Get maximum/minimum level from the collection
      ilvmin = collct_getvalue_int (p_rcollection,'NLMIN')
      ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')
      
      ! Get the interlevel projection structure and the temporary vector
      ! from the collection.
      ! Our 'parent' prepared there how to interpolate the solution on the
      ! fine grid to coarser grids.
      p_rprojection => collct_getvalue_ilvp(p_rcollection,'ILVPROJECTION')
      p_rvectorTemp => collct_getvalue_vecsca(p_rcollection,'RTEMPSCALAR')
      
      ! Prepare the matrix assembly on level < NLMAX.
      !
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
      
      ! Loop through all the levels to set up the matrices.
      ! Note that we loop from the last but one level to the minimum
      ! level, because we must interpolate the solution vector
      ! from the finest to the coarsest for getting the nonlinarity
      ! on all the levels.
      DO i=ilvmax-1,ilvmin,-1
      
        ! Get the destination matrix on that level
        p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',i)
        
        ! Get the temporary vector on level i. Will receive the solution
        ! vector on that level. 
        p_rvectorCoarse => collct_getvalue_vec (p_rcollection,'RTEMPVEC',i)
        
        ! Get the solution vector on level i+1. This is either the temporary
        ! vector on that level, or the solution vector on the maximum level.
        IF (i .LT. ilvmax-1) THEN
          p_rvectorFine => collct_getvalue_vec (p_rcollection,'RTEMPVEC',i+1)
        ELSE
          p_rvectorFine => rx
        END IF
        
        ! Interpolate the solution from the finer grid to the coarser grid.
        ! The interpolation is configured in the interlevel projection
        ! structure we got from the collection.
        CALL mlprj_performInterpolation (p_rprojection,p_rvectorCoarse, &
                                         p_rvectorFine,p_rvectorTemp)
        
        ! Now we have the solution vector on the current level.
        ! Next step is to build the matrix entries with the discretisation routine.
        !
        ! Get the system matrix on the current level.
        p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',i)
        
        ! Put a reference to the solution vector on the current level into 
        ! the collection. This way, we inform the callback
        ! routine of the matrix assembly about the solution vector to use
        ! fot the nonlinear term.
        CALL collct_setvalue_vec(p_rcollection,'RX',p_rvectorCoarse,.TRUE.)
        
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
        
        ! Remove RX from the collection, not needed there anymore.
        CALL collct_deletevalue(p_rcollection,'RX')
        
        ! Implement discrete boundary conditions into the matrix. 
        ! Call the appropriate matrix filter to modify the system matrix
        ! according to the attached discrete boundary conditions.
        CALL matfil_discreteBC (p_rmatrix)
        
        ! Sort the matrix according to the attached sorting strategy -
        ! if there's a sorting strategy attached at all
        ! (this is prepared by the application).
        ! The sorting is important,
        ! - to make the matrices compatible to the vector rd (which is also
        !   resorted later)
        ! - to make the matrices consistent to those we put into the
        !   linear solver; we put sorted matrices into the linear solver!
        CALL lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.TRUE.,&
                                ABS(p_rmatrix%RmatrixBlock(1,1)%isortStrategy))
        
      END DO

      ! Sort the matrix on the maximum level - don't forget this!
      p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',ilvmax)
      CALL lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.TRUE.,&
                              ABS(p_rmatrix%RmatrixBlock(1,1)%isortStrategy))

      ! Ok, system matrices on all levels are assembled now.
      ! Now we turn to invokle the linear solver for preconditioning...
      !
      ! Resort the vector rd before solving the corresponding linear
      ! system/perform the preconditioning. Use p_rvectorTemp as temporary
      ! vector for that purpose - it's prepared large enough.
      CALL lsysbl_sortVectorInSitu (rd,p_rvectorTemp,.TRUE.)
      
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
      ! in the next iteration.
      CALL linsol_doneData (p_rsolverNode)

      ! Unsort our rd again so that it's in the state it was before.
      CALL lsysbl_sortVectorInSitu (rd,p_rvectorTemp,.FALSE.)
      
      ! Unsort the structure of all matrices without unsorting the entries.
      ! This of course means throwing away all matrices, but we don't need
      ! them anymore - they are reassembled in the next sweep.
      DO i=ilvmin,ilvmax
        p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',i)
        CALL lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.FALSE.,&
                                -ABS(p_rmatrix%RmatrixBlock(1,1)%isortStrategy))
      END DO        

    END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d6_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a nonlinear solver.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    TYPE(t_vectorBlock), TARGET :: rtempBlock
    TYPE(t_vectorScalar), TARGET :: rtempVectorSc
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    
    TYPE(t_nlsolNode) :: rnlSol
    INTEGER :: ierror,ilvmin,ilvmax,i
    TYPE(t_linsolNode), POINTER :: p_rsolverNode,p_rpreconditioner,p_rsmoother
    TYPE(t_linsolNode), POINTER :: p_rcoarseGridSolver

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    TYPE(t_matrixBlock), DIMENSION(rproblem%ilvmax) :: Rmatrices

    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain

    ! One level of multigrid
    TYPE(t_linsolMGLevelInfo), POINTER :: p_rlevelInfo
    
    ! An interlevel projection structure for changing levels
    TYPE(t_interlevelProjectionBlock) :: rprojection
    
    INTEGER(PREC_VECIDX) :: imaxmem
    
    ! Min/Max level?
    ilvmin = rproblem%ilvmin
    ilvmax = rproblem%ilvmax
    
    ! Now we have to build up the level information for multigrid.
    !
    ! At first, initialise a standard interlevel projection structure. We
    ! can use the same structure for all levels. Therefore it's enough
    ! to initialise one structure using the RHS vector on the finest
    ! level to specify the shape of the PDE-discretisation.
    CALL mlprj_initProjectionVec (rprojection,rproblem%rrhs)
    
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
    ! At first, set up the linear solver as usual.
    !
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
    CALL collct_setvalue_ilvp(rproblem%rcollection,'ILVPROJECTION',&
                              rprojection,.TRUE.)
    
    ! Then set up smoothers / coarse grid solver:
    imaxmem = 0
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
        ! Set up Jacobi smoother for multigrid would be:
        ! CALL linsol_initJacobi (p_rsmoother)

        ! Set up an ILU smoother for multigrid with damping parameter 0.7,
        ! 2 smoothing steps (pre- and postsmoothing). As the problem is
        ! badly conditioned, we need even ILU(4) to get this problem solved,
        ! otherwise the smoother is diverging!
        ! note that if the trapezoidal rule is used for setting up the matrix,
        ! one would even need ILU(6)!!!
        CALL linsol_initMILUs1x1 (p_rsmoother,4,0.0_DP)
        CALL linsol_convertToSmoother (p_rsmoother,2,0.7_DP)
      END IF
    
      ! Add the level.
      CALL linsol_addMultigridLevel (p_rlevelInfo,p_rsolverNode, rprojection,&
                                     p_rsmoother,p_rsmoother,p_rcoarseGridSolver)
                                     
      ! How much memory is necessary for performing the level change?
      ! We ourself must build nonlinear matrices on multiple levels and have
      ! to interpolate the solution vector from finer level to coarser ones.
      ! We need temporary memory for this purpose...
      IF (i .GT. ilvmin) THEN
        ! Pass the system metrices on the coarse/fine grid to 
        ! mlprj_getTempMemoryMat to specify the discretisation structures
        ! of all equations in the PDE there.
        imaxmem = MAX(imaxmem,mlprj_getTempMemoryMat (rprojection,&
                              rproblem%RlevelInfo(i-1)%rmatrix,&
                              rproblem%RlevelInfo(i)%rmatrix))
      END IF
    END DO
    
    ! Set up a scalar temporary vector that we need for building up nonlinear
    ! matrices. It must be at least as large as MAXMEM and NEQ(finest level),
    ! as we use it for resorting vectors, too.
    CALL lsyssc_createVector (rtempVectorSc,MAX(imaxmem,rproblem%rrhs%NEQ),&
                              .FALSE.,ST_DOUBLE)
    CALL collct_setvalue_vecsca(rproblem%rcollection,'RTEMPSCALAR',rtempVectorSc,.TRUE.)
    
    ! Before attaching the matrices to the solver and the initialisation of
    ! the problem structure, sort the matrix structure on all levels
    ! according top the associated permutation. Don't sort the
    ! entries - there are none!
    DO i=ilvmin,ilvmax
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      CALL lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.FALSE.,&
                              ABS(p_rmatrix%RmatrixBlock(1,1)%isortStrategy))
    END DO        
    
    ! Attach the system matrices to the solver.
    !
    ! We copy our matrices to a big matrix array and transfer that
    ! to the setMatrices routines. This intitialises then the matrices
    ! on all levels according to that array.
    Rmatrices(ilvmin:ilvmax) = rproblem%RlevelInfo(ilvmin:ilvmax)%rmatrix
    CALL linsol_setMatrices(p_RsolverNode,Rmatrices(ilvmin:ilvmax))
    
    ! Initialise structure of the solver. This allows the
    ! solver to allocate memory / perform some precalculation
    ! to the problem.
    ! In fact, solvers like UMFPACK use this for a symbolic factorisation
    ! of the matrix.
    CALL linsol_initStructure (p_rsolverNode, ierror)
    IF (ierror .NE. LINSOL_ERR_NOERROR) STOP

    ! Unsort the matrix structure again. The matrices stay in unsorted form
    ! until the entries are assembled.
    ! Remark: This makes the matrices inconsistent to those attached to the
    !  linear solver! So before invoking the linear solver, the matrices
    !  must be sorted to make their entries consistent!!!
    DO i=ilvmin,ilvmax
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      CALL lsyssc_sortMatrix (p_rmatrix%RmatrixBlock(1,1),.FALSE.,&
                              -ABS(p_rmatrix%RmatrixBlock(1,1)%isortStrategy))
    END DO        
    
    ! Put the prepared solver node to the collection for later use.
    CALL collct_setvalue_linsol(rproblem%rcollection,'LINSOLVER',p_rsolverNode,.TRUE.)
    
    ! Create a temporary vector we need for the nonliner iteration.
    CALL lsysbl_createVecBlockIndirect (rproblem%rrhs, rtempBlock, .FALSE.)

    ! The nonlinear solver structure rnlSol is initialised by the default
    ! initialisation with all necessary information to solve the problem.
    ! We call the nonlinear solver directly. For preconditioning
    ! and defect calculation, we use our own callback routine.
    rnlSol%ioutputLevel = 2
    CALL nlsol_performSolve(rnlSol,rproblem%rvector,rproblem%rrhs,rtempBlock,&
                            b1d6_getDefect,b1d6_precondDefect,&
                            rcollection=rproblem%rcollection)

    ! Release the temporary vector(s)
    CALL lsysbl_releaseVector (rtempBlock)
    CALL lsyssc_releaseVector (rtempVectorSc)
    
    ! Remove the solver node from the collection - not needed anymore there
    CALL collct_deletevalue(rproblem%rcollection,'LINSOLVER')
    
    ! Remove the temporary vector from the collection
    CALL collct_deletevalue(rproblem%rcollection,'RTEMPSCALAR')
    
    ! Remove the interlevel projection structure
    CALL collct_deletevalue(rproblem%rcollection,'ILVPROJECTION')
    
    ! Clean up the linear solver, release all memory, remove the solver node
    ! from memory.
    CALL linsol_releaseSolver (p_rsolverNode)
    
    ! Release the multilevel projection structure.
    CALL mlprj_doneProjection (rprojection)
    
    CALL output_lbrk()
    CALL output_line ('Nonlinear solver statistics')
    CALL output_line ('---------------------------')
    CALL output_line ('Intial defect: '//TRIM(sys_sdEL(rnlSol%DinitialDefect(1),15)))
    CALL output_line ('Final defect:  '//TRIM(sys_sdEL(rnlSol%DfinalDefect(1),15)))
    CALL output_line ('#Iterations:   '//TRIM(sys_siL(rnlSol%iiterations,10)))

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d6_postprocessing (rproblem)
  
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
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,'gmv/u6.gmv')
    
    CALL lsyssc_getbase_double (p_rvector%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'sol',UCD_VAR_STANDARD, p_Ddata)
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d6_doneMatVec (rproblem)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
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
      
      ! Remove the temp vector that was used for interpolating the solution
      ! from higher to lower levels in the nonlinear iteration.
      IF (i .LT. rproblem%ilvmax) THEN
        CALL lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVector)
        CALL collct_deletevalue(rproblem%rcollection,'RTEMPVEC',i)
      END IF
      
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

  SUBROUTINE b1d6_doneBC (rproblem)
  
!<description>
  ! Releases discrete boundary conditions from the heap.
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

  SUBROUTINE b1d6_doneDiscretisation (rproblem)
  
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
      ! Delete the block discretisation together with the associated
      ! scalar spatial discretisations....
      CALL spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%p_rdiscretisation)
     
      ! and remove the allocated block discretisation structure from the heap.
      DEALLOCATE(rproblem%RlevelInfo(i)%p_rdiscretisation)
    END DO
    
  END SUBROUTINE
    
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE b1d6_doneParamTriang (rproblem)
  
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
    
    ! Finally release the domain.
    CALL boundary_release (rproblem%p_rboundary)
    
    CALL collct_deleteValue(rproblem%rcollection,'NLMAX')
    CALL collct_deleteValue(rproblem%rcollection,'NLMIN')

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE burgers1d6
  
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
!</description>

!</subroutine>

    ! NLMIN receives the minimal level where to discretise for supporting
    ! the solution process.
    ! NLMAX receives the level where we want to solve.
    INTEGER :: NLMIN,NLMAX
    
    ! A problem structure for our problem
    TYPE(t_problem), POINTER :: p_rproblem
    
    INTEGER :: i
    
    ! Ok, let's start. 
    !
    ! We want to solve our Laplace problem on level...

    NLMIN = 3
    NLMAX = 7
    
    ! Allocate the problem structure -- it's rather large
    ALLOCATE(p_rproblem)
    ALLOCATE(p_rproblem%RlevelInfo(1:NLMAX))

    ! Initialise the collection
    CALL collct_init (p_rproblem%rcollection)
    DO i=1,NLMAX
      CALL collct_addlevel_all (p_rproblem%rcollection)
    END DO

    ! So now the different steps - one after the other.
    
    ! Initialisation. 
    CALL b1d6_initParamTriang (NLMIN,NLMAX,p_rproblem)
    CALL b1d6_initDiscretisation (p_rproblem)    
    CALL b1d6_initMatVec (p_rproblem)    
    CALL b1d6_initDiscreteBC (p_rproblem)
    
    ! Implementation of boundary conditions
    CALL b1d6_implementBC (p_rproblem)
    
    ! Solve the problem
    CALL b1d6_solve (p_rproblem)
    
    ! Postprocessing
    CALL b1d6_postprocessing (p_rproblem)
    
    ! Cleanup
    CALL b1d6_doneMatVec (p_rproblem)
    CALL b1d6_doneBC (p_rproblem)
    CALL b1d6_doneDiscretisation (p_rproblem)
    CALL b1d6_doneParamTriang (p_rproblem)
    
    ! Print some statistical data about the collection - anything forgotten?
    PRINT *
    PRINT *,'Remaining collection statistics:'
    PRINT *,'--------------------------------'
    PRINT *
    CALL collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection and the problem structure.
    CALL collct_done (p_rproblem%rcollection)
    
    DEALLOCATE(p_rproblem%RlevelInfo)
    DEALLOCATE(p_rproblem)
    
  END SUBROUTINE

END MODULE

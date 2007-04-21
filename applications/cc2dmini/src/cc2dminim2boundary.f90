!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2boundary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the definition of the analytic boundary conditions
!# as well as discretisation routines for the boundary.
!#
!# The following files can be found here:
!#
!# 1.) c2d2_initAnalyticBC
!#     -> Initialise analytic boundary conditions
!#
!# 2.) c2d2_initDiscreteBC
!#     -> Discretise analytic boundary conditions, create discrete
!#        boundary conditions
!#
!# 3.) c2d2_doneBC
!#     -> Release discrete and analytic boundary conditions
!#
!# 4.) c2d2_implementBC
!#     -> Implement discrete boundary conditions into solution/RHS vectors
!#        and matrix on finest level
!#
!# </purpose>
!##############################################################################

MODULE cc2dminim2boundary

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
  USE nonlinearsolver
  USE paramlist
  
  USE collection
  USE convection
    
  USE cc2dminim2basic
  USE cc2dmini_callback
  
  IMPLICIT NONE

CONTAINS
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initAnalyticBC (rproblem)
  
!<description>
  ! This initialises the analytic bonudary conditions of the problem
  ! and saves them to the problem structure.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables

    ! A set of variables describing the analytic boundary conditions.    
    TYPE(t_boundaryRegion) :: rboundaryRegion
    TYPE(t_bcRegion), POINTER :: p_rbcRegion
    
    ! A pointer to the discretisation structure with the data.
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
    
    ! A pointer to the domain
    TYPE(t_boundary), POINTER :: p_rboundary
    
    INTEGER :: i
    INTEGER, DIMENSION(2) :: IvelComp
  
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
    
    ! Edge 4 of boundary component 1. 
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                       rboundaryRegion,p_rbcRegion)
      
    ! The whole 2nd boundary component - if it exists.
    IF (boundary_igetNBoundComp(p_rboundary) .GE. 2) THEN
      CALL boundary_createRegion(p_rboundary,2,0,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,1,&
                                        rboundaryRegion,p_rbcRegion)
    END IF
      
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

    ! The whole 2nd boundary component - if it exists.
    IF (boundary_igetNBoundComp(p_rboundary) .GE. 2) THEN
      CALL boundary_createRegion(p_rboundary,2,0,rboundaryRegion)
      CALL bcond_newDirichletBConRealBD (rproblem%p_rboundaryConditions,2,&
                                        rboundaryRegion,p_rbcRegion)
    END IF
      
    !CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    !IvelComp = (/1,2/)
    !CALL bcond_newPressureDropBConRealBD (rproblem%p_rboundaryConditions,IvelComp,&
    !                                      rboundaryRegion,p_rbcRegion)
      
    ! Install these analytic boundary conditions into all discretisation
    ! structures on all levels.
                               
    DO i=rproblem%NLMIN,rproblem%NLMAX
      
      ! Ask the problem structure to give us the discretisation structure...
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
      ! and inform the discretisation which analytic boundary conditions to use:
      p_rdiscretisation%p_rboundaryConditions => rproblem%p_rboundaryConditions

    END DO
    
    ! The pressure does not need boundary conditions.
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initDiscreteBC (rproblem)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
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
    
    DO i=rproblem%NLMIN,rproblem%NLMAX
    
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
      !
      ! On maximum level, discretrise everything. On lower level, discretise
      ! only for the implementation into the matrices and defect vector. 
      ! That's enough, as the lower levels are only used for preconditioning 
      ! of defect vectors.
      
      NULLIFY(rproblem%RlevelInfo(i)%p_rdiscreteBC)
      IF (i .EQ. rproblem%NLMAX) THEN
        CALL bcasm_discretiseBC (p_rdiscretisation, &
                                 rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                                .FALSE.,getBoundaryValues, &
                                rproblem%rcollection)
      ELSE
        CALL bcasm_discretiseBC (p_rdiscretisation, &
                                 rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                                .FALSE.,getBoundaryValues, &
                                rproblem%rcollection,BCASM_DISCFORDEFMAT)
      END IF
                                       
      ! Hang the pointer into the the matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%p_rdiscreteBC
      
      p_rmatrix%p_rdiscreteBC => p_rdiscreteBC
      
      ! Also hang in the boundary conditions into the temporary vector that is
      ! used for the creation of solutions on lower levels.
      ! This allows us to filter this vector when we create it.
      rproblem%RlevelInfo(i)%rtempVector%p_rdiscreteBC => p_rdiscreteBC
      
    END DO

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%NLMAX)%p_rdiscreteBC
    
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector
    
    p_rrhs%p_rdiscreteBC => p_rdiscreteBC
    p_rvector%p_rdiscreteBC => p_rdiscreteBC
                
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem astructure saving problem-dependent information.
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
    ilvmax = rproblem%NLMAX
    p_rrhs    => rproblem%rrhs   
    p_rvector => rproblem%rvector
    
    ! Implement discrete boundary conditions into RHS vector by 
    ! filtering the vector.
    CALL vecfil_discreteBCrhs (p_rrhs)

    ! Implement discrete boundary conditions into solution vector by
    ! filtering the vector.
    CALL vecfil_discreteBCsol (p_rvector)
    
    ! Implement discrete boundary conditions into the matrices on all 
    ! levels, too.
    ! In fact, this modifies the B-matrices. The A-matrices are overwritten
    ! later and must then be modified again!
    DO i=rproblem%NLMIN ,rproblem%NLMAX
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      CALL matfil_discreteBC (p_rmatrix)
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_doneBC (rproblem)
  
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

    DO i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Release our discrete version of the boundary conditions
      CALL bcasm_releaseDiscreteBC (rproblem%RlevelInfo(i)%p_rdiscreteBC)

      ! ...and also the corresponding analytic description.
      CALL bcond_doneBC (rproblem%p_rboundaryConditions)
    END DO
    
  END SUBROUTINE

END MODULE

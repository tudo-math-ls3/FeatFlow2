!##############################################################################
!# ****************************************************************************
!# <name> heatcond_boundarycondition </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is routines to maintain the boundary conditions -- in analytic
!# and discrete form.
!#
!# The following routines can be found here:
!# 1.) hc5_initDiscreteBC
!#     -> Discretise boundary conditions or update discretisation of boundary
!#        conditions.
!#
!# 2.) hc5_implementBC
!#     -> Implement boundary conditions into matrices/vectors.
!#
!# 3.) hc5_doneBC
!#     -> Release boundary conditions and associated memory.
!# </purpose>
!##############################################################################

MODULE heatcond_boundarycondition

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
  
  USE heatcond_basic
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE hc5_initDiscreteBC (rproblem)
  
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
    TYPE(t_vectorBlock), POINTER :: p_rrhs
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation

    ! Pointer to structure for saving discrete BC's:
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
    
    ! A set of variables describing the analytic boundary conditions.    
    TYPE(t_boundaryRegion) :: rboundaryRegion

    ! A pointer to the domain
    TYPE(t_boundary), POINTER :: p_rboundary
  
    ! Get the domain from the discretisation
    p_rboundary => rproblem%p_rboundary

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
      
      ! For implementing boundary conditions, we use a 'filter technique with
      ! discretised boundary conditions'. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions - if it does not exist. If it exists, clear it to prepare
      ! it for the next step.
      IF (.NOT. ASSOCIATED(rproblem%RlevelInfo(i)%p_rdiscreteBC)) THEN
        ALLOCATE(rproblem%RlevelInfo(i)%p_rdiscreteBC)
        CALL bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%p_rdiscreteBC)
      ELSE
        CALL bcasm_clearDiscreteBC(rproblem%RlevelInfo(i)%p_rdiscreteBC)
      END IF
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
         rboundaryRegion,rproblem%RlevelInfo(i)%p_rdiscreteBC,&
         getBoundaryValues)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
      CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
         rboundaryRegion,rproblem%RlevelInfo(i)%p_rdiscreteBC,&
         getBoundaryValues)
                               
      ! Edge 3 of boundary component 1.
      CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
      CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
         rboundaryRegion,rproblem%RlevelInfo(i)%p_rdiscreteBC,&
         getBoundaryValues)
      
      ! Edge 4 of boundary component 1. That's it.
      CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
      CALL bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
         rboundaryRegion,rproblem%RlevelInfo(i)%p_rdiscreteBC,&
         getBoundaryValues)

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
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
  
  ! A pointer to the solution vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
  
  ! A pointer to the RHS vector.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,ilvmax
  
    ! A pointer to the system matrix and the RHS vector as well as 
    ! the discretisation
    TYPE(t_matrixBlock), POINTER :: p_rmatrix
    
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

  SUBROUTINE hc5_doneBC (rproblem)
  
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
      IF (ASSOCIATED(rproblem%RlevelInfo(i)%p_rdiscreteBC)) THEN
        CALL bcasm_releaseDiscreteBC (rproblem%RlevelInfo(i)%p_rdiscreteBC)
        DEALLOCATE(rproblem%RlevelInfo(i)%p_rdiscreteBC)
      END IF
    END DO
    
  END SUBROUTINE

END MODULE

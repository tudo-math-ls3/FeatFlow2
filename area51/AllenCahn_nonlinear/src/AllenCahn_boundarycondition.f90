!##############################################################################
!# ****************************************************************************
!# <name> AllenCahn_boundarycondition </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is routines to maintain the boundary conditions -- in analytic
!# and discrete form.
!#
!# The following routines can be found here:
!# 1.) AC5_initDiscreteBC
!#     -> Discretise boundary conditions or update discretisation of boundary
!#        conditions.
!#
!# 2.) AC5_implementBC
!#     -> Implement boundary conditions into matrices/vectors.
!#
!# 3.) AC5_doneBC
!#     -> Release boundary conditions and associated memory.
!# </purpose>
!##############################################################################

MODULE AllenCahn_boundarycondition

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use sortstrategy
  use coarsegridcorrection
  use ucd
  use timestepping
  use genoutput
  
  use collection
  use paramlist
    
  use AllenCahn_callback
  use AllenCahn_basic
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE AC5_initDiscreteBC (rACproblem)

!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: i

    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), POINTER :: p_rmatrix
    type(t_vectorBlock), POINTER :: p_rrhs
    type(t_blockDiscretisation), POINTER :: p_rdiscretisation

    ! Pointer to structure for saving discrete BC's:
    type(t_discreteBC), POINTER :: p_rdiscreteBC
    
    ! A set of variables describing the analytic boundary conditions.
    type(t_boundaryRegion) :: rboundaryRegion

    ! A pointer to the domain
    type(t_boundary), POINTER :: p_rboundary
  
    ! Get the domain from the discretisation
    p_rboundary => rACproblem%rboundary

    ! Put the current simulation time as parameter "TIME" into the collection.
    ! Also set Dquickaccess (1) to the simulation time for faster access by the
    ! callback routine.
    rACproblem%rcollection%Dquickaccess (1) = rACproblem%rtimedependence%dtime
    call collct_setvalue_real(rACproblem%rcollection,'TIME',&
         rACproblem%rtimedependence%dtime,.TRUE.)
    
    DO i=rACproblem%NLMIN,rACproblem%NLMAX
    
      ! Get our matrix from the problem structure.
      p_rmatrix => rACproblem%RlevelInfo(i)%rmatrix
      
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => p_rmatrix%p_rblockDiscrTrial
      
      ! For implementing boundary conditions, we use a 'filter technique with
      ! discretised boundary conditions'. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions - if it does not exist. if it exists, clear it to prepare
      ! it for the next step.
      if (.NOT. ASSOCIATED(rACproblem%RlevelInfo(i)%p_rdiscreteBC)) THEN
        ALLOCATE(rACproblem%RlevelInfo(i)%p_rdiscreteBC)
        call bcasm_initDiscreteBC(rACproblem%RlevelInfo(i)%p_rdiscreteBC)
      else
        call bcasm_clearDiscreteBC(rACproblem%RlevelInfo(i)%p_rdiscreteBC)
      end if
      !
      ! We 'know' already (from the problem definition) that we have four boundary
      ! segments in the domain. Each of these, we want to use for enforcing
      ! some kind of boundary condition.
      !

!~~~~~~~~Note that we only have Neumann BC, so we delete all Dirichlet BC part
!~~~~~~~We have changed the Boundary condition to be Neumann type~~~~~~~~~~~~~

      ! We ask the bondary routines to create a 'boundary region' - which is
      ! simply a part of the boundary corresponding to a boundary segment.
      ! A boundary region roughly contains the type, the min/max parameter value
      ! and whether the endpoints are inside the region or not.
      call boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
      
      ! We use this boundary region and specify that we want to have Dirichlet
      ! boundary there. The following call does the following:
      ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
      
      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! Hang the pointer into the the matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rACproblem%RlevelInfo(i)%p_rdiscreteBC
      
      p_rmatrix%p_rdiscreteBC => p_rdiscreteBC
      
    end DO
    
    ! Remove the "TIME"-parameter from the collection again.
    call collct_deletevalue (rACproblem%rcollection,'TIME')

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rACproblem%RlevelInfo(rACproblem%NLMAX)%p_rdiscreteBC
    
    p_rrhs    => rACproblem%rrhs
    p_rrhs%p_rdiscreteBC => p_rdiscreteBC
                
  end SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE AC5_implementBC (rACproblem,rvector,rrhs,dtimeWeight)
  
!<description>
  ! Implements boundary conditions into the RHS, a given solution vector
  ! and into the system matrices on all levels specified in rACproblem.
!</description>

!<input>
  ! Time stepping weight. Standard is 1.0.
  REAL(DP) :: dtimeWeight
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
  
  ! A pointer to the solution vector.
  type(t_vectorBlock), INTENT(INOUT) :: rvector
  
  ! A pointer to the RHS vector.
  type(t_vectorBlock), INTENT(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,NLMAX
  
    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), POINTER :: p_rmatrix
    
    ! Pointer to structure for saving discrete BC's:
    type(t_discreteBC), POINTER :: p_rdiscreteBC
    
    ! Get our the right hand side and solution from the problem structure
    ! on the finest level
    NLMAX = rACproblem%NLMAX
    
    ! Next step is to implement boundary conditions into the RHS,
    ! solution and matrix. This is done using a vector/matrix filter
    ! for discrete boundary conditions.
    ! The discrete boundary conditions are already attached to the
    ! vectors/matrix. call the appropriate vector/matrix filter that
    ! modifies the vectors/matrix according to the boundary conditions.

    p_rdiscreteBC => rACproblem%RlevelInfo(rACproblem%NLMAX)%p_rdiscreteBC
    call vecfil_discreteBCrhs (rrhs, p_rdiscreteBC)
    call vecfil_discreteBCsol (rvector, p_rdiscreteBC)


    ! Implement discrete boundary conditions into the matrices on all
    ! levels, too. call the appropriate matrix filter to modify
    ! all matrices according to the attached discrete boundary conditions.
    DO i=rACproblem%NLMIN,rACproblem%NLMAX
      p_rmatrix => rACproblem%RlevelInfo(i)%rmatrix
      call matfil_discreteBC (p_rmatrix)
    end DO

  end SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE AC5_doneBC (rACproblem)
  
!<description>
  ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_ACproblem), INTENT(INOUT), TARGET :: rACproblem
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i

    DO i=rACproblem%NLMAX,rACproblem%NLMIN,-1
      ! Release our discrete version of the boundary conditions
      if (ASSOCIATED(rACproblem%RlevelInfo(i)%p_rdiscreteBC)) THEN
        call bcasm_releaseDiscreteBC (rACproblem%RlevelInfo(i)%p_rdiscreteBC)
        DEALLOCATE(rACproblem%RlevelInfo(i)%p_rdiscreteBC)
      end if
    end DO
    
  end SUBROUTINE

end MODULE

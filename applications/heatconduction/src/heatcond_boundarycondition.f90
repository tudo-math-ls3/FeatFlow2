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

module heatcond_boundarycondition

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use discretebc
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
    
  use heatcond_callback
  
  use heatcond_basic
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine hc5_initDiscreteBC (rproblem)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Pointer to structure for saving discrete BC`s:
    type(t_discreteBC), pointer :: p_rdiscreteBC
    
    ! A set of variables describing the analytic boundary conditions.
    type(t_boundaryRegion) :: rboundaryRegion

    ! A pointer to the domain
    type(t_boundary), pointer :: p_rboundary
  
    ! Get the domain from the discretisation
    p_rboundary => rproblem%rboundary

    ! Put the current simulation time as parameter "TIME" into the collection.
    ! Also set Dquickaccess (1) to the simulation time for faster access by the
    ! callback routine.
    rproblem%rcollection%Dquickaccess (1) = rproblem%rtimedependence%dtime
    call collct_setvalue_real(rproblem%rcollection,'TIME',&
         rproblem%rtimedependence%dtime,.true.)
    
    do i=rproblem%ilvmin,rproblem%ilvmax
    
      ! Get our matrix from the problem structure.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => p_rmatrix%p_rblockDiscrTrial
      
      ! For implementing boundary conditions, we use a `filter technique with
      ! discretised boundary conditions`. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions - if it does not exist. If it exists, clear it to prepare
      ! it for the next step.
      if (.not. associated(rproblem%RlevelInfo(i)%p_rdiscreteBC)) then
        allocate(rproblem%RlevelInfo(i)%p_rdiscreteBC)
        call bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%p_rdiscreteBC)
      else
        call bcasm_clearDiscreteBC(rproblem%RlevelInfo(i)%p_rdiscreteBC)
      end if
      !
      ! We 'know' already (from the problem definition) that we have four boundary
      ! segments in the domain. Each of these, we want to use for enforcing
      ! some kind of boundary condition.
      !
      ! We ask the boundary routines to create a 'boundary region' - which is
      ! simply a part of the boundary corresponding to a boundary segment.
      ! A boundary region roughly contains the type, the min/max parameter value
      ! and whether the endpoints are inside the region or not.
      call boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
      
      ! We use this boundary region and specify that we want to have Dirichlet
      ! boundary there. The following call does the following:
      ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
      !   We specify icomponent='1' to indicate that we set up the
      !   Dirichlet BC`s for the first (here: one and only) component in the
      !   solution vector.
      ! - Discretise the boundary condition so that the BC`s can be applied
      !   to matrices and vectors
      ! - Add the calculated discrete BC`s to rdiscreteBC for later use.
      call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
         rboundaryRegion,rproblem%RlevelInfo(i)%p_rdiscreteBC,&
         getBoundaryValues)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
         rboundaryRegion,rproblem%RlevelInfo(i)%p_rdiscreteBC,&
         getBoundaryValues)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
         rboundaryRegion,rproblem%RlevelInfo(i)%p_rdiscreteBC,&
         getBoundaryValues)
      
      ! Edge 4 of boundary component 1. That is it.
      call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (p_rdiscretisation,1,&
         rboundaryRegion,rproblem%RlevelInfo(i)%p_rdiscreteBC,&
         getBoundaryValues)

      ! Hang the pointer into the the matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%p_rdiscreteBC
      
      p_rmatrix%p_rdiscreteBC => p_rdiscreteBC
      
    end do
    
    ! Remove the "TIME"-parameter from the collection again.
    call collct_deletevalue (rproblem%rcollection,'TIME')

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%ilvmax)%p_rdiscreteBC
    
    p_rrhs    => rproblem%rrhs
    p_rrhs%p_rdiscreteBC => p_rdiscreteBC
                
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine hc5_implementBC (rproblem,rvector,rrhs,dtimeWeight)
  
!<description>
  ! Implements boundary conditions into the RHS, a given solution vector
  ! and into the system matrices on all levels specified in rproblem.
!</description>

!<input>
  ! Time stepping weight. Standard is 1.0.
  real(DP) :: dtimeWeight
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
  
  ! A pointer to the solution vector.
  type(t_vectorBlock), intent(inout) :: rvector
  
  ! A pointer to the RHS vector.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,ilvmax
  
    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    
    ! Pointer to structure for saving discrete BC`s:
    type(t_discreteBC), pointer :: p_rdiscreteBC
    
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
    call vecfil_discreteBCrhs (rrhs,p_rdiscreteBC)
    call vecfil_discreteBCsol (rvector,p_rdiscreteBC)

    ! Implement discrete boundary conditions into the matrices on all
    ! levels, too. Call the appropriate matrix filter to modify
    ! all matrices according to the attached discrete boundary conditions.
    do i=rproblem%ilvmin,rproblem%ilvmax
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      call matfil_discreteBC (p_rmatrix)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine hc5_doneBC (rproblem)
  
!<description>
  ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    do i=rproblem%ilvmax,rproblem%ilvmin,-1
      ! Release our discrete version of the boundary conditions
      if (associated(rproblem%RlevelInfo(i)%p_rdiscreteBC)) then
        call bcasm_releaseDiscreteBC (rproblem%RlevelInfo(i)%p_rdiscreteBC)
        deallocate(rproblem%RlevelInfo(i)%p_rdiscreteBC)
      end if
    end do
    
  end subroutine

end module

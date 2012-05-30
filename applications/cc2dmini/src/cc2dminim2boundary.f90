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

module cc2dminim2boundary

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
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  
  use collection
  use convection
    
  use cc2dminim2basic
  use cc2dmini_callback
  
  implicit none

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_discretiseBC (rdiscretisation,rdiscreteBC)
  
!<description>
  ! This routine discretises the current boundary conditions for
  ! the discretisation specified by rdiscretisation.
!</description>

!<input>
  ! A discretisation structure specifying the current discretisation.
  type(t_blockDiscretisation), intent(in) :: rdiscretisation
!</input>

!<inputoutput>
  ! A structuree that receives the discretised boundary conditions.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_boundary), pointer :: p_rboundary
    
    p_rboundary => rdiscretisation%p_rboundary

    ! We first set up the boundary conditions for the X-velocity, then those
    ! of the Y-velocity.
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
    
    ! The endpoint of this segment should also be Dirichlet. We set this by
    ! changing the region properties in rboundaryRegion.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    ! We use this boundary region and specify that we want to have Dirichlet
    ! boundary there. The following call does the following:
    ! - Create Dirichlet boundary conditions on the region rboundaryRegion.
    !   We specify icomponent='1' to indicate that we set up the
    !   Dirichlet BC's for the first (here: one and only) component in the
    !   solution vector.
    ! - Discretise the boundary condition so that the BC's can be applied
    !   to matrices and vectors
    ! - Add the calculated discrete BC's to rdiscreteBC for later use.
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)
                              
    ! Edge 2 is Neumann boundary, so it's commented out.
    ! CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rdiscretisation,1,&
    !                                    rboundaryRegion,rdiscreteBC,&
    !                                    getBoundaryValues)
                              
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)
    
    ! Edge 4 of boundary component 1. That's it.
    call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)

    ! The whole 2nd boundary component - if it exists.
    if (boundary_igetNBoundComp(p_rboundary) .ge. 2) then
      call boundary_createRegion(p_rboundary,2,0,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues)
    end if

    ! Now continue with defining the boundary conditions of the Y-velocity:
    !
    ! Define edge 1.
    call boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    
    ! Edge with start- and endpoint.
    rboundaryRegion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
    
    ! As we define the Y-velocity, we now set icomponent=2 in the following call.
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)
                              
    ! Edge 2 is Neumann boundary, so it's commented out.
    ! CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    ! CALL bcasm_newDirichletBConRealBD (rdiscretisation,2,&
    !                                    rboundaryRegion,rdiscreteBC,&
    !                                    getBoundaryValues)
                              
    ! Edge 3 of boundary component 1.
    call boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)
    
    ! Edge 4 of boundary component 1. That's it.
    call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                      rboundaryRegion,rdiscreteBC,&
                                      getBoundaryValues)

    ! The whole 2nd boundary component - if it exists.
    if (boundary_igetNBoundComp(p_rboundary) .ge. 2) then
      call boundary_createRegion(p_rboundary,2,0,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_initDiscreteBC (rproblem)
  
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
  type(t_vectorBlock), pointer :: p_rrhs,p_rvector
  type(t_blockDiscretisation), pointer :: p_rdiscretisation

  ! Pointer to structure for saving discrete BC's:
  type(t_discreteBC), pointer :: p_rdiscreteBC
    
    do i=rproblem%NLMIN,rproblem%NLMAX
    
      ! Get our velocity matrix from the problem structure.
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => p_rmatrix%p_rblockDiscrTrial
      
      ! For implementing boundary conditions, we use a 'filter technique with
      ! discretised boundary conditions'. This means, we first have to calculate
      ! a discrete version of the analytic BC, which we can implement into the
      ! solution/RHS vectors using the corresponding filter.
      !
      ! Create a t_discreteBC structure where we store all discretised boundary
      ! conditions.
      call bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%rdiscreteBC)
      
      ! Discretise the boundary conditions.
      call c2d2_discretiseBC (p_rdiscretisation,rproblem%RlevelInfo(i)%rdiscreteBC)

      ! Hang the pointer into the the matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%rdiscreteBC
      
      p_rmatrix%p_rdiscreteBC => p_rdiscreteBC
      
      ! Also hang in the boundary conditions into the temporary vector that is
      ! used for the creation of solutions on lower levels.
      ! This allows us to filter this vector when we create it.
      rproblem%RlevelInfo(i)%rtempVector%p_rdiscreteBC => p_rdiscreteBC
      
    end do

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%NLMAX)%rdiscreteBC
    
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    
    p_rrhs%p_rdiscreteBC => p_rdiscreteBC
    p_rvector%p_rdiscreteBC => p_rdiscreteBC
                
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_implementBC (rproblem)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,ilvmax
  
    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector
    
    ! Get our the right hand side and solution from the problem structure
    ! on the finest level
    ilvmax = rproblem%NLMAX
    p_rrhs    => rproblem%rrhs
    p_rvector => rproblem%rvector
    
    ! Implement discrete boundary conditions into RHS vector by
    ! filtering the vector.
    call vecfil_discreteBCrhs (p_rrhs)

    ! Implement discrete boundary conditions into solution vector by
    ! filtering the vector.
    call vecfil_discreteBCsol (p_rvector)
    
    ! Implement discrete boundary conditions into the matrices on all
    ! levels, too.
    ! In fact, this modifies the B-matrices. The A-matrices are overwritten
    ! later and must then be modified again!
    do i=rproblem%NLMIN ,rproblem%NLMAX
      p_rmatrix => rproblem%RlevelInfo(i)%rmatrix
      call matfil_discreteBC (p_rmatrix)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_doneBC (rproblem)
  
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

    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Release our discrete version of the boundary conditions
      call bcasm_releaseDiscreteBC (rproblem%RlevelInfo(i)%rdiscreteBC)
    end do
    
  end subroutine

end module

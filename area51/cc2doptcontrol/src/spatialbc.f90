!##############################################################################
!# ****************************************************************************
!# <name> spatialbc </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the definition of the analytic boundary conditions
!# as well as discretisation routines for the boundary.
!#
!# The following files can be found here:
!#
!# 1.) cc_initDiscreteBC
!#     -> Discretise analytic boundary conditions, create discrete
!#        boundary conditions
!#
!# 2.) cc_updateDiscreteBC
!#     -> Reassemble discrete boundary conditions
!#
!# 3.) cc_doneBC
!#     -> Release discrete and analytic boundary conditions
!#
!# 4.) cc_implementBC
!#     -> Implement discrete boundary conditions into solution/RHS vectors
!#        and matrix on finest level
!#
!# </purpose>
!##############################################################################

module spatialbc

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
    
  use basicstructures
  use spatialbcdef
  use user_callback
  
  implicit none

!contains
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_initDiscreteBC (rproblem,dtime,rvector,rrhs)
!
!!<description>
!  ! This calculates the discrete version of the boundary conditions and
!  ! assigns it to the system matrix and RHS vector.
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT), target :: rproblem
!
!  ! A vector structure for the solution vector. The discrete BC structures are
!  ! attached to that.
!  type(t_vectorBlock), intent(INOUT) :: rvector
!
!  ! A vector structure for the RHS vector. The discrete BC structures are
!  ! attached to that.
!  type(t_vectorBlock), intent(INOUT) :: rrhs
!
!  ! Current simulation time.
!  real(dp), intent(in) :: dtime
!!</inputoutput>
!
!!</subroutine>
!
!  ! local variables
!  integer :: i
!
!  ! A pointer to the system matrix and the RHS vector as well as
!  ! the discretisation
!  type(t_blockDiscretisation), pointer :: p_rdiscretisation
!
!  ! Pointer to structure for saving discrete BC's:
!  type(t_discreteBC), pointer :: p_rdiscreteBC
!  type(t_discreteFBC), pointer :: p_rdiscreteFBC
!
!    do i=rproblem%NLMIN,rproblem%NLMAX
!
!      ! From the matrix or the RHS we have access to the discretisation and the
!      ! analytic boundary conditions.
!      p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisation
!
!      ! Initialise the BC/FBC structures.
!      allocate(rproblem%RlevelInfo(i)%p_rdiscreteBC)
!      call bcasm_initDiscreteBC(rproblem%RlevelInfo(i)%p_rdiscreteBC)
!
!      allocate(rproblem%RlevelInfo(i)%p_rdiscreteFBC)
!      call bcasm_initDiscreteFBC(rproblem%RlevelInfo(i)%p_rdiscreteFBC)
!      ! Hang the pointer into the the matrix. That way, these
!      ! boundary conditions are always connected to that matrix and that
!      ! vector.
!      p_rdiscreteBC => rproblem%RlevelInfo(i)%p_rdiscreteBC
!
!      ! Also hang in the boundary conditions into the temporary vector that is
!      ! used for the creation of solutions on lower levels.
!      ! This allows us to filter this vector when we create it.
!      rproblem%RlevelInfo(i)%rtempVector1%p_rdiscreteBC => p_rdiscreteBC
!      rproblem%RlevelInfo(i)%rtempVector2%p_rdiscreteBC => p_rdiscreteBC
!      rproblem%RlevelInfo(i)%rtempVector3%p_rdiscreteBC => p_rdiscreteBC
!
!      ! The same for the fictitious boundary boundary conditions.
!      p_rdiscreteFBC => rproblem%RlevelInfo(i)%p_rdiscreteFBC
!      rproblem%RlevelInfo(i)%rtempVector1%p_rdiscreteBCfict => p_rdiscreteFBC
!      rproblem%RlevelInfo(i)%rtempVector2%p_rdiscreteBCfict => p_rdiscreteFBC
!      rproblem%RlevelInfo(i)%rtempVector3%p_rdiscreteBCfict => p_rdiscreteFBC
!
!    end do
!
!    ! On the finest level, attach the discrete BC also
!    ! to the solution and RHS vector. They need it to be compatible
!    ! to the matrix on the finest level.
!    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%NLMAX)%p_rdiscreteBC
!
!    rrhs%p_rdiscreteBC => p_rdiscreteBC
!    rvector%p_rdiscreteBC => p_rdiscreteBC
!
!    ! The same with the fictitious boundary BC's
!    p_rdiscreteFBC => rproblem%RlevelInfo(rproblem%NLMAX)%p_rdiscreteFBC
!    rrhs%p_rdiscreteBCfict => p_rdiscreteFBC
!    rvector%p_rdiscreteBCfict => p_rdiscreteFBC
!
!    ! Call the update routine to assemble the BC's.
!    call cc_updateDiscreteBC (rproblem,dtime)
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_updateDiscreteBC (rproblem,dtime,rdiscretisation,cbctype,&
!      rdiscreteBC,rdiscreteFBC)
!
!!<description>
!  ! This updates the discrete version of the boundary conditions. The BC's
!  ! are reassembled according to the current situation of the simulation.
!!</description>
!
!!<input>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT), target :: rproblem
!
!  ! Current simulation time.
!  real(dp), intent(in) :: dtime
!
!  ! Discretisation structure identifying the FEM space.
!  type(t_blockDiscretisation), intent(in) :: rdiscretisation
!
!  ! Type of boundary condition to assemble.
!  ! CCDISCBC_PRIMAL: BC for the primal space.
!  ! CCDISCBC_DUAL: BC for the dual space.
!  ! CCDISCBC_PRIMALDUAL: BC for primal and dual space.
!  integer, intent(in) :: cbctype
!
!!</input>
!
!!<inputoutout>
!  ! Structure identifying and receiving the discrete boudary conditions.
!  type(t_discreteBC), intent(inout) :: rdiscreteBC
!
!  ! Structure identifying and receiving the discrete fictitious
!  ! boundary conditions
!  type(t_discreteFBC), intent(inout) :: rdiscreteBC
!!</inputoutput>
!
!!</subroutine>
!
!    ! Initialise the collection for the assembly process with callback routines.
!    ! Basically, this stores the simulation time in the collection if the
!    ! simulation is nonstationary.
!    call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)
!
!    ! Clear the last discretised BC's. We are reassembling them.
!    call bcasm_clearDiscreteBC(rproblem%RlevelInfo(i)%p_rdiscreteBC)
!    call bcasm_clearDiscreteFBC(rproblem%RlevelInfo(i)%p_rdiscreteFBC)
!
!    ! Assemble boundary conditions
!    call cc_assembleBDconditions (rproblem,dtime,p_rdiscretisation,&
!        rproblem%RlevelInfo(i)%p_rdiscreteBC,rproblem%rcollection)
!
!    ! Assemble the boundary conditions of fictitious boundary components
!    call cc_assembleFBDconditions (rproblem,dtime,p_rdiscretisation,&
!        rproblem%RlevelInfo(i)%p_rdiscreteFBC,rproblem%rcollection)
!
!    ! Clean up the collection (as we are done with the assembly, that's it.
!    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_implementBC (rproblem,rvector,rrhs,rdefect)
!
!!<description>
!  ! Implements boundary conditions into the RHS and into a given solution vector.
!  ! Implements boundary conditions into the system matrix on every level
!  ! defined by rproblem.
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT), target :: rproblem
!
!  ! A vector structure of a solution vector. The discrete BC's are implemented
!  ! into that.
!  type(t_vectorBlock), intent(INOUT), optional :: rvector
!
!  ! A vector structure of a RHS vector. The discrete BC's are implamented into that.
!  type(t_vectorBlock), intent(INOUT), optional :: rrhs
!
!  ! A vector structure of a defect vector. The discrete BC's are implamented into that.
!  type(t_vectorBlock), intent(INOUT), optional :: rdefect
!
!!</inputoutput>
!
!!</subroutine>
!
!  ! local variables
!  integer :: ilvmax
!
!    ! Get our the right hand side and solution from the problem structure
!    ! on the finest level
!    ilvmax = rproblem%NLMAX
!
!    if (present(rvector)) then
!
!      ! Implement discrete boundary conditions into solution vector by
!      ! filtering the vector.
!      call vecfil_discreteBCsol (rvector)
!
!      ! Implement discrete boundary conditions of fictitioous boundary comnponents
!      ! into solution vector by filtering the vector.
!      call vecfil_discreteFBCsol (rvector)
!
!    end if
!
!    if (present(rdefect)) then
!
!      ! Implement discrete boundary conditions into solution vector by
!      ! filtering the vector.
!      call vecfil_discreteBCdef (rdefect)
!
!      ! Implement discrete boundary conditions of fictitioous boundary comnponents
!      ! into solution vector by filtering the vector.
!      call vecfil_discreteFBCdef (rdefect)
!
!    end if
!
!    if (present(rrhs)) then
!
!      ! Implement pressure drop boundary conditions into RHS vector
!      ! if there are any.
!      call vecfil_discreteNLPDropBCrhs (rrhs)
!
!      ! Implement discrete boundary conditions into RHS vector by
!      ! filtering the vector.
!      call vecfil_discreteBCrhs (rrhs)
!
!      ! Implement discrete boundary conditions of fictitious boundary components
!      ! into RHS vector by filtering the vector.
!      call vecfil_discreteFBCrhs (rrhs)
!
!    end if
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_doneBC (rproblem)
!
!!<description>
!  ! Releases discrete and analytic boundary conditions from the heap.
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT), target :: rproblem
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables
!    integer :: i
!
!    do i=rproblem%NLMAX,rproblem%NLMIN,-1
!      ! Release our discrete version of the boundary conditions
!      call bcasm_releaseDiscreteBC (rproblem%RlevelInfo(i)%p_rdiscreteBC)
!      deallocate(rproblem%RlevelInfo(i)%p_rdiscreteBC)
!
!      ! as well as the discrete version of the BC's for fictitious boundaries
!      call bcasm_releaseDiscreteFBC (rproblem%RlevelInfo(i)%p_rdiscreteFBC)
!      deallocate(rproblem%RlevelInfo(i)%p_rdiscreteFBC)
!    end do
!
!  end subroutine

end module

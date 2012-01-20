!##############################################################################
!# ****************************************************************************
!# <name> ccboundarycondition </name>
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
!#     -> Release discrete boundary conditions
!#
!# 4.) cc_implementBC
!#     -> Implement discrete boundary conditions into solution/RHS vectors
!#        and matrix on finest level
!#
!# </purpose>
!##############################################################################

module ccboundarycondition

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
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  
  use collection
  use convection
    
  use ccbasic
  use ccboundaryconditionparser
  use cccallback
  
  implicit none

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_initDiscreteBC (rproblem,rvector,rrhs)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
  
  ! A vector structure for the solution vector. The discrete BC structures are 
  ! attached to that.
  type(t_vectorBlock), intent(inout) :: rvector

  ! A vector structure for the RHS vector. The discrete BC structures are 
  ! attached to that.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Pointer to structure for saving discrete BC`s:
    type(t_discreteBC), pointer :: p_rdiscreteBC
    type(t_discreteFBC), pointer :: p_rdiscreteFBC
    
    do i=rproblem%NLMIN,rproblem%NLMAX
    
      ! Hang the pointer into the matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%rdynamicInfo%rdiscreteBC

      ! Also hang in the boundary conditions into the temporary vector that is
      ! used for the creation of solutions on lower levels.
      ! This allows us to filter this vector when we create it.
      rproblem%RlevelInfo(i)%rtempVector%p_rdiscreteBC => p_rdiscreteBC
      
      ! The same for the fictitious boundary boundary conditions.
      p_rdiscreteFBC => rproblem%RlevelInfo(i)%rdynamicInfo%rdiscreteFBC
      rproblem%RlevelInfo(i)%rtempVector%p_rdiscreteBCfict => p_rdiscreteFBC
      
    end do

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%NLMAX)%rdynamicInfo%rdiscreteBC
    
    rrhs%p_rdiscreteBC => p_rdiscreteBC
    rvector%p_rdiscreteBC => p_rdiscreteBC

    ! The same with the fictitious boundary BC`s
    p_rdiscreteFBC => rproblem%RlevelInfo(rproblem%NLMAX)%rdynamicInfo%rdiscreteFBC
    rrhs%p_rdiscreteBCfict => p_rdiscreteFBC
    rvector%p_rdiscreteBCfict => p_rdiscreteFBC
    
    ! Call the update routine to assemble the BC`s.
    call cc_updateDiscreteBC (rproblem)
                
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_updateDiscreteBC (rproblem)
  
!<description>
  ! This updates the discrete version of the boundary conditions. The BC`s
  ! are reassembled according to the current situation of the simulation.
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
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

    do i=rproblem%NLMIN,rproblem%NLMAX
    
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisation
      
      ! Clear the last discretised BC`s. We are reassembling them.
      call bcasm_clearDiscreteBC(rproblem%RlevelInfo(i)%rdynamicInfo%rdiscreteBC)
      call bcasm_clearDiscreteFBC(rproblem%RlevelInfo(i)%rdynamicInfo%rdiscreteFBC)
      
      ! Assemble boundary conditions
      call cc_assembleBDconditions (rproblem,p_rdiscretisation,&
          rproblem%RlevelInfo(i)%rdynamicInfo,rproblem%rcollection)
      
      ! Assemble the boundary conditions of fictitious boundary components
      call cc_assembleFBDconditions (rproblem,p_rdiscretisation,&
          rproblem%RlevelInfo(i)%rdynamicInfo,rproblem%rcollection)

    end do

    ! Clean up the collection (as we are done with the assembly, that is it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_implementBC (rproblem,rvector,rrhs,&
      bsolvector,brhsvector)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! A vector structure for the solution vector. The discrete BC`s are implemented
  ! into that.
  type(t_vectorBlock), intent(inout) :: rvector

  ! A vector structure for the RHS vector. The discrete BC`s are implamented into that.
  type(t_vectorBlock), intent(inout) :: rrhs
  
  ! Whether to implement the BC`s into the solution vector or not
  logical, intent(in) :: bsolvector

  ! Whether to implement the BC`s into the solution vector or not
  logical, intent(in) :: brhsvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ilvmax
  
    ! Get the right hand side and solution from the problem structure
    ! on the finest level
    ilvmax = rproblem%NLMAX
    
    if (bsolvector) then
    
      ! Implement discrete boundary conditions into solution vector by
      ! filtering the vector.
      call vecfil_discreteBCsol (rvector)
      
      ! Implement discrete boundary conditions of fictitioous boundary comnponents
      ! into solution vector by filtering the vector.
      call vecfil_discreteFBCsol (rvector)
      
    end if
    
    if (brhsvector) then
    
      ! Implement pressure drop boundary conditions into RHS vector
      ! if there are any.
      call vecfil_discreteNLPDropBCrhs (rrhs)
      
      ! Implement discrete boundary conditions into RHS vector by 
      ! filtering the vector.
      call vecfil_discreteBCrhs (rrhs)

      ! Implement discrete boundary conditions of fictitious boundary components
      ! into RHS vector by filtering the vector.
      call vecfil_discreteFBCrhs (rrhs)

    end if
    
    ! Implementation of boundary conditions into matrices deactivated.
    ! It is simply not used as the global matrices are actually not used outside
    ! of the nonlinear iteration!
    !
    !  IF (bmatrices) THEN
    !  
    !    ! Implement discrete boundary conditions into the matrices on all 
    !    ! levels, too.
    !    ! In fact, this modifies the B-matrices. The A-matrices are overwritten
    !    ! later and must then be modified again!
    !    DO i=rproblem%NLMIN ,rproblem%NLMAX
    !      p_rmatrix => rproblem%RlevelInfo(i)%rpreallocatedSystemMatrix
    !      CALL matfil_discreteBC (p_rmatrix)  ! standard boundary conditions
    !      CALL matfil_discreteFBC (p_rmatrix)  ! fictitious boundary boundary conditions
    !    END DO
    !    
    !  END IF

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneBC (rproblem)
  
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
      call bcasm_clearDiscreteBC (rproblem%RlevelInfo(i)%rdynamicInfo%rdiscreteBC)
      
      ! as well as the discrete version of the BC`s for fictitious boundaries
      call bcasm_clearDiscreteFBC (rproblem%RlevelInfo(i)%rdynamicInfo%rdiscreteFBC)

      ! Release the Dirichlet edges.      
      if (rproblem%RlevelInfo(i)%rdynamicInfo%hedgesDirichletBC .ne. ST_NOHANDLE) then
        call storage_free (rproblem%RlevelInfo(i)%rdynamicInfo%hedgesDirichletBC)
      end if
      
    end do
    
  end subroutine

end module

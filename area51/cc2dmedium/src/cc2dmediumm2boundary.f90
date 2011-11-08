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
!# 3.) c2d2_updateDiscreteBC
!#     -> Reassemble discrete boundary conditions
!#
!# 4.) c2d2_doneBC
!#     -> Release discrete and analytic boundary conditions
!#
!# 5.) c2d2_implementBC
!#     -> Implement discrete boundary conditions into solution/RHS vectors
!#        and matrix on finest level
!#
!# </purpose>
!##############################################################################

module cc2dmediumm2boundary

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
    
  use cc2dmediumm2basic
  use cc2dmediumm2boundarydef
  use cc2dmedium_callback
  
  implicit none

contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_initAnalyticBC (rproblem)
  
!<description>
  ! This initialises the analytic boundary conditions of the problem
  ! and saves them to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    integer :: i
  
    ! Initialise the boundary condition by parsing the parameter files.
    ! This initialises rproblem%p_rboundaryConditions by parsing DAT-file
    ! parameters in the parameter list rproblem%rparamList
    call c2d2_parseBDconditions (rproblem)

    ! Initialise the boundary conditions of fictitious boundary components
    call c2d2_parseFBDconditions (rproblem)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_initDiscreteBC (rproblem,rvector,rrhs)
  
!<description>
  ! This calculates the discrete version of the boundary conditions and
  ! assigns it to the system matrix and RHS vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! A vector structure for the solution vector. The discrete BC structures are
  ! attached to that.
  type(t_vectorBlock), intent(INOUT) :: rvector

  ! A vector structure for the RHS vector. The discrete BC structures are
  ! attached to that.
  type(t_vectorBlock), intent(INOUT) :: rrhs
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

  ! A pointer to the system matrix and the RHS vector as well as
  ! the discretisation
  type(t_blockDiscretisation), pointer :: p_rdiscretisation

  ! Pointer to structure for saving discrete BC's:
  type(t_discreteBC), pointer :: p_rdiscreteBC
  type(t_discreteFBC), pointer :: p_rdiscreteFBC
    
    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    call c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)

    do i=rproblem%NLMIN,rproblem%NLMAX
    
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
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
      
      nullify(rproblem%RlevelInfo(i)%p_rdiscreteBC)
      if (i .eq. rproblem%NLMAX) then
        call bcasm_discretiseBC (p_rdiscretisation, &
                                 rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                                .false.,getBoundaryValues, &
                                rproblem%rcollection)
      else
        call bcasm_discretiseBC (p_rdiscretisation, &
                                 rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                                .false.,getBoundaryValues, &
                                rproblem%rcollection,BCASM_DISCFORDEFMAT)
      end if
                
      ! The same way, discretise the fictitious boundary conditions and hang
      ! them in into the matrix/vectors.
      nullify(rproblem%RlevelInfo(i)%p_rdiscreteFBC)
      if (i .eq. rproblem%NLMAX) then
        call bcasm_discretiseFBC (p_rdiscretisation,&
                                  rproblem%RlevelInfo(i)%p_rdiscreteFBC,.false., &
                                  getBoundaryValuesFBC,rproblem%rcollection)
      else
        call bcasm_discretiseFBC (p_rdiscretisation,&
                                  rproblem%RlevelInfo(i)%p_rdiscreteFBC,.false., &
                                  getBoundaryValuesFBC,rproblem%rcollection,&
                                  BCASM_DISCFORDEFMAT)
      end if

      ! Hang the pointer into the the matrix. That way, these
      ! boundary conditions are always connected to that matrix and that
      ! vector.
      p_rdiscreteBC => rproblem%RlevelInfo(i)%p_rdiscreteBC

      rproblem%RlevelInfo(i)%rpreallocatedSystemMatrix%p_rdiscreteBC => p_rdiscreteBC
      
      ! Also hang in the boundary conditions into the temporary vector that is
      ! used for the creation of solutions on lower levels.
      ! This allows us to filter this vector when we create it.
      rproblem%RlevelInfo(i)%rtempVector%p_rdiscreteBC => p_rdiscreteBC
      
      ! The same for the fictitious boundary boundary conditions.
      p_rdiscreteFBC => rproblem%RlevelInfo(i)%p_rdiscreteFBC
      rproblem%RlevelInfo(i)%rpreallocatedSystemMatrix%p_rdiscreteBCfict => &
          p_rdiscreteFBC
      rproblem%RlevelInfo(i)%rtempVector%p_rdiscreteBCfict => p_rdiscreteFBC
      
    end do

    ! On the finest level, attach the discrete BC also
    ! to the solution and RHS vector. They need it to be compatible
    ! to the matrix on the finest level.
    p_rdiscreteBC => rproblem%RlevelInfo(rproblem%NLMAX)%p_rdiscreteBC
    
    rrhs%p_rdiscreteBC => p_rdiscreteBC
    rvector%p_rdiscreteBC => p_rdiscreteBC

    ! The same with the fictitious boundary BC's
    p_rdiscreteFBC => rproblem%RlevelInfo(rproblem%NLMAX)%p_rdiscreteFBC
    rrhs%p_rdiscreteBCfict => p_rdiscreteFBC
    rvector%p_rdiscreteBCfict => p_rdiscreteFBC
    
    ! Clean up the collection (as we are done with the assembly, that's it.
    call c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)
                
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_updateDiscreteBC (rproblem, bforce)
  
!<description>
  ! This updates the discrete version of the boundary conditions. The BC's
  ! are reassembled according to the current situation of the simulation.
!</description>

!<input>
  ! Normally, only those BC's are update that are marked as 'nonstatic'
  ! (e.g. for time-dependent simulation). When bforce=TRUE, all discrete
  ! BC's are rebuild, which might be necessary in optimisation environments
  ! when testing multiple stationary configurations.
  logical, intent(IN) :: bforce
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
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
    call c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)

    do i=rproblem%NLMIN,rproblem%NLMAX
    
      ! From the matrix or the RHS we have access to the discretisation and the
      ! analytic boundary conditions.
      p_rdiscretisation => rproblem%RlevelInfo(i)%p_rdiscretisation
      
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
      !
      ! Note that we do not NULLIFY the pointer to the discrete BC's here.
      ! That's the only difference to the initDiscreteBC routine above,
      ! as this instructs bcasm_discretiseBC to update nonstatic BC's.
      
      if (i .eq. rproblem%NLMAX) then
        call bcasm_discretiseBC (p_rdiscretisation, &
                                 rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                                 bforce,getBoundaryValues, &
                                 rproblem%rcollection)
      else
        call bcasm_discretiseBC (p_rdiscretisation, &
                                 rproblem%RlevelInfo(i)%p_rdiscreteBC, &
                                 bforce,getBoundaryValues, &
                                 rproblem%rcollection,BCASM_DISCFORDEFMAT)
      end if
                
      ! The same way, discretise the fictitious boundary conditions and hang
      ! them in into the matrix/vectors.
      if (i .eq. rproblem%NLMAX) then
        call bcasm_discretiseFBC (p_rdiscretisation,&
                                  rproblem%RlevelInfo(i)%p_rdiscreteFBC,bforce, &
                                  getBoundaryValuesFBC,rproblem%rcollection)
      else
        call bcasm_discretiseFBC (p_rdiscretisation,&
                                  rproblem%RlevelInfo(i)%p_rdiscreteFBC,bforce, &
                                  getBoundaryValuesFBC,rproblem%rcollection,&
                                  BCASM_DISCFORDEFMAT)
      end if

    end do

    ! Clean up the collection (as we are done with the assembly, that's it.
    call c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_implementBC (rproblem,rvector,rrhs,&
      bsolvector,brhsvector)
  
!<description>
  ! Implements boundary conditions into the RHS and into a given solution vector.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem

  ! A vector structure for the solution vector. The discrete BC's are implemented
  ! into that.
  type(t_vectorBlock), intent(INOUT) :: rvector

  ! A vector structure for the RHS vector. The discrete BC's are implamented into that.
  type(t_vectorBlock), intent(INOUT) :: rrhs
  
  ! Whether to implement the BC's into the solution vector or not
  logical, intent(IN) :: bsolvector

  ! Whether to implement the BC's into the solution vector or not
  logical, intent(IN) :: brhsvector
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,ilvmax
  
  ! A pointer to the system matrix and the RHS vector as well as
  ! the discretisation
  type(t_matrixBlock), pointer :: p_rmatrix
    
    ! Get our the right hand side and solution from the problem structure
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
    ! It's simply not used as the global matrices are actually not used outside
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

  subroutine c2d2_doneBC (rproblem)
  
!<description>
  ! Releases discrete and analytic boundary conditions from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  type(t_fparser), pointer :: p_rparser

    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Release our discrete version of the boundary conditions
      call bcasm_releaseDiscreteBC (rproblem%RlevelInfo(i)%p_rdiscreteBC)
      
      ! as well as the discrete version of the BC's for fictitious boundaries
      call bcasm_releaseDiscreteFBC (rproblem%RlevelInfo(i)%p_rdiscreteFBC)

      ! ...and also the corresponding analytic description.
      call bcond_doneBC (rproblem%p_rboundaryConditions)
    end do
    
    ! Release the parser object with all the expressions to be evaluated
    ! on the boundary.
    p_rparser => collct_getvalue_pars (rproblem%rcollection, BDC_BDPARSER, &
                                       0, SEC_SBDEXPRESSIONS)
    call fparser_release (p_rparser)
    deallocate (p_rparser)

    ! Remove the boundary value parser from the collection
    call collct_deletevalue (rproblem%rcollection, BDC_BDPARSER)
    
    ! Remove the boundary-expression section we added earlier,
    ! with all their content.
    call collct_deletesection(rproblem%rcollection,SEC_SBDEXPRESSIONS)
    
  end subroutine

end module

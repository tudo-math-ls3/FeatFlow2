!##############################################################################
!# ****************************************************************************
!# <name> spacetimebc </name>
!# ****************************************************************************
!#
!# <purpose>
!# Realises boundary conditions in space and time.
!# </purpose>
!##############################################################################

module spacetimebc

  use fsystem
  use storage
  use genoutput
  use linearalgebra
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use collection
  use vectorio
  use derivatives
  use timediscretisation
  
  use stdoperators
  use spacetimevectors
  use bilinearformevaluation
  use collection
  use vectorfilters
  use discretebc
  use bcassembly
  use boundary

  use physics
  
  use callback
  
  implicit none

  ! Defect BC's
  integer, parameter :: SPOP_DEFECT = 0

  ! RHS BC's
  integer, parameter :: SPOP_RHS = 1

  ! Solution BC's
  integer, parameter :: SPOP_SOLUTION = 2

  ! Encapsules the boundary conditions in space and time.
  type t_spacetimeBC
  
    ! Underlying space discretisation
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr
    
    ! Underlying time discretisation
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    
    ! Physics of the problem.
    type(t_physics), pointer :: p_rphysics

  end type

contains

  ! ***************************************************************************

  subroutine spop_getPrimalDualTime (rtimeDiscr,istep,dprimalTime,ddualTime)

  ! Calculates the primal/dual time corresponding to timestep istep.
  
  ! Underlying time discretisation 
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr

  ! Timestep
  integer, intent(in) :: istep

  ! Corresponding time in the primal solution
  real(DP), intent(out) :: dprimalTime

  ! Corresponding time in the dual solution
  real(DP), intent(out) :: ddualTime
      
    ! local variables
    real(DP) :: dtimeend,dtstep,dtimestart
      
    select case (rtimeDiscr%ctype)
    case (TDISCR_ONESTEPTHETA)
      call tdiscr_getTimestep(rtimediscr,istep-1,dtimeend,dtstep,dtimestart)
      
      ! The primal is always the endpoint
      dprimalTime = dtimeend
      
      ! The dual is shifted by dtheta -- except for in the beginning.
      if (istep .eq. 1) then
        ddualTime = dtimeEnd
      else
        ddualTime = dtimeend*rtimeDiscr%dtheta + dtimestart*(1.0_DP-rtimeDiscr%dtheta)
      end if
      
    end select
      
  end subroutine

  ! ***************************************************************************

  subroutine spop_createBC (rspaceDiscr,rtimeDiscr,rphysics,rbc)

  ! Initialises space-time boundary conditions
  
  ! Underlying spatial discretisation structure
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr

  ! Underlying time discretisation structure
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr

  ! Physics of the problem
  type(t_physics), intent(in), target :: rphysics

  ! Boundary condition structure.
  type(t_spacetimeBC), intent(out) :: rbc
    
    rbc%p_rphysics => rphysics
    rbc%p_rspaceDiscr => rspaceDiscr
    rbc%p_rtimeDiscr => rtimeDiscr
    
  end subroutine

  ! ***************************************************************************

  subroutine spop_releaseBC (rbc)

  ! Cleans up a BC structure.
  
  ! Boundary condition structure to be cleaned up
  type(t_spacetimeBC), intent(inout) :: rbc
    
    nullify(rbc%p_rphysics)
    nullify(rbc%p_rspaceDiscr)
    nullify(rbc%p_rtimeDiscr)
  
  end subroutine
  
  ! ***************************************************************************

  subroutine spop_assembleSpaceBCtime (rbc, dtime, icomponent, rdiscreteBC)

  ! Assemble the boundary conditions at a specific time.
  
  ! The underlying domain.
  
  ! Boundary condition structure 
  type(t_spacetimeBC), intent(in) :: rbc
  
  ! Time where to evaluate the boundary conditions.
  real(DP), intent(in) :: dtime

  ! Component of the equation whose BC's should be assembled.
  integer, intent(in) :: icomponent  
  
  ! Boundary condition structure to place the BC's to.
  ! Must already be initialised. Old existing BC's are not deleted,
  ! the new BC's are appended.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
  
    ! local variables
    integer :: ibc, isegment
    type(t_boundaryRegion) :: rregion
    type(t_collection) :: rcollection
  
    select case (rbc%p_rphysics%cequation)
    case (0)
      do ibc = 1,boundary_igetNBoundComp(rbc%p_rspaceDiscr%p_rboundary)
        do isegment = 1,boundary_igetNsegments(rbc%p_rspaceDiscr%p_rboundary,ibc)
          ! Get the segment and create Dirichlet BC's there.
          call boundary_createRegion (rbc%p_rspaceDiscr%p_rboundary, ibc, isegment, &
              rregion)
                             
          ! IQuickAccess(1) is the number of the primal equation. Here we only have one!
          ! IquickAccess(2) is the equation.
          ! DquickAccess(1) is the time.
          rcollection%IquickAccess(1) = icomponent
          rcollection%IquickAccess(2) = rbc%p_rphysics%cequation
          rcollection%DquickAccess(1) = dtime
          call bcasm_newDirichletBConRealBd (rbc%p_rspaceDiscr, &
              icomponent, rregion, rdiscreteBC, cb_getBoundaryValuesOptC, rcollection)
        end do
      end do
    end select
  
  end subroutine
  
  ! ***************************************************************************

  subroutine spop_assembleSpaceBC (rbc, istep, ctype, rdiscreteBC)

  ! Assemble the boundary conditions at a specific timestep for all components.
  
  ! The underlying domain.
  
  ! Boundary condition structure.
  ! Should be clean before calling this routine.
  type(t_spacetimeBC), intent(in) :: rbc
  
  ! Timestep number
  integer, intent(in) :: istep

  ! Type of BC's. A SPOP_xxxx constant.
  integer, intent(in) :: ctype
  
  ! Boundary condition structure to place the BC's to.
  ! Must already be initialised. Old existing BC's are not deleted,
  ! the new BC's are appended.
  type(t_discreteBC), intent(inout) :: rdiscreteBC

    ! local variables  
    real(DP) :: dprimalTime,ddualTime

    call spop_getPrimalDualTime (rbc%p_rtimeDiscr,istep,dprimalTime,ddualTime)

    select case (ctype)
    case (SPOP_DEFECT,SPOP_SOLUTION)
      call spop_assembleSpaceBCtime (rbc, dprimalTime, 1, rdiscreteBC)
      call spop_assembleSpaceBCtime (rbc, ddualTime, 2, rdiscreteBC)
    case (SPOP_RHS)
      call spop_assembleSpaceBCtime (rbc, ddualTime, 1, rdiscreteBC)
      call spop_assembleSpaceBCtime (rbc, dprimalTime, 2, rdiscreteBC)
    end select
      
  end subroutine
  
  ! ***************************************************************************

  subroutine spop_applyBC (rbc, ctype, rvector)

  ! Applies boundary conditions to a defect vector rvector.
  
  ! Boundary condition structure.
  type(t_spacetimeBC), intent(in) :: rbc
  
  ! Type of BC's. A SPOP_xxxx constant.
  integer, intent(in) :: ctype
  
  ! Vector where to apply boundary conditions to.
  type(t_spaceTimeVector), intent(inout) :: rvector
  
  ! Type of BC's to apply.
  
    ! local variables
    integer :: istep
    type(t_discreteBC) :: rdiscreteBC
    type(t_vectorblock) :: rtempVec
    real(DP), dimension(:), pointer :: p_Ddata
    
    ! Init the BC structure
    call bcasm_initDiscreteBC(rdiscreteBC)
    
    ! Temp vector
    call lsysbl_createVectorBlock (rbc%p_rspaceDiscr,rtempVec,.false.)
    call lsysbl_getbase_double(rtempVec,p_Ddata)
    
    ! Loop through the timesteps.
    do istep = 1,rvector%NEQtime
     
      ! Assemble the BC's. Be careful with the dime discretisation!
      call bcasm_clearDiscreteBC (rdiscreteBC)
      call spop_assembleSpaceBC (rbc, istep, ctype, rdiscreteBC)
      
      ! Implement.
      call sptivec_getTimestepData (rvector, istep, rtempVec)

      select case (ctype)
      case (SPOP_DEFECT)
        call vecfil_discreteBCdef (rtempVec,rdiscreteBC)
        
        ! DEBUG!!!
        if (istep .eq. 1) then
          call lsyssc_clearVector (rtempVec%RvectorBlock(1))
        end if
        if (istep .eq. rvector%NEQtime) then
          call lsyssc_clearVector (rtempVec%RvectorBlock(2))
        end if
      case (SPOP_RHS)
        call vecfil_discreteBCrhs (rtempVec,rdiscreteBC)
      case (SPOP_SOLUTION)
        call vecfil_discreteBCsol (rtempVec,rdiscreteBC)
      end select
        
      call sptivec_setTimestepData (rvector, istep, rtempVec)
    
    end do
    
    ! Release
    call lsysbl_releaseVector(rtempVec)
    call bcasm_releaseDiscreteBC(rdiscreteBC)
  
  end subroutine

end module

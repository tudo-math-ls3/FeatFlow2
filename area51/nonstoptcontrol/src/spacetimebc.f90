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
  
  use analyticprojection

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

  subroutine spop_getPrimalDualTime (rtimeDiscr,istep,dtimePrimal,dtimeDual)

  ! Calculates the primal/dual time corresponding to timestep istep.
  
  ! Underlying time discretisation 
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr

  ! Timestep
  integer, intent(in) :: istep

  ! Corresponding time in the primal solution
  real(DP), intent(out) :: dtimePrimal

  ! Corresponding time in the dual solution
  real(DP), intent(out) :: dtimeDual
      
    ! local variables
    real(DP) :: dtimeend,dtstep,dtimestart
    integer :: ithetaschemetype
    
    ithetaschemetype = rtimediscr%itag
      
    select case (rtimeDiscr%ctype)
    case (TDISCR_ONESTEPTHETA)
      call tdiscr_getTimestep(rtimediscr,istep-1,dtimeend,dtstep,dtimestart)
      
      ! The primal is always the endpoint
      dtimePrimal = dtimeend
      dtimeDual = dtimeend
      
      if (ithetaschemetype .eq. 1) then
        ! The dual is shifted by dtheta -- except for in the beginning.
        if (istep .ne. 1) then
          dtimeDual = dtimeend*rtimeDiscr%dtheta + dtimestart*(1.0_DP-rtimeDiscr%dtheta)
        end if
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
      ! 2D Heat equation
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
          rcollection%IquickAccess(3) = rbc%p_rphysics%creferenceProblem
          rcollection%DquickAccess(1) = dtime
          rcollection%DquickAccess(2) = rbc%p_rphysics%doptControlAlpha
          rcollection%DquickAccess(3) = rbc%p_rphysics%doptControlGamma
          call bcasm_newDirichletBConRealBd (rbc%p_rspaceDiscr, &
              icomponent, rregion, rdiscreteBC, cb_getBoundaryValuesOptC, rcollection)
        end do
      end do

    case (1)
      ! Stokes equation
      do ibc = 1,boundary_igetNBoundComp(rbc%p_rspaceDiscr%p_rboundary)
      
        ! Last segment is Neumann!
        do isegment = 1,boundary_igetNsegments(rbc%p_rspaceDiscr%p_rboundary,ibc)
          if (isegment .ne. 2) then
            ! Get the segment and create Dirichlet BC's there.
            call boundary_createRegion (rbc%p_rspaceDiscr%p_rboundary, ibc, isegment, &
                rregion)
            rregion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
                               
            ! IQuickAccess(1) is the number of the primal equation. Here we only have one!
            ! IquickAccess(2) is the equation.
            ! DquickAccess(1) is the time.
            rcollection%IquickAccess(1) = icomponent
            rcollection%IquickAccess(2) = rbc%p_rphysics%cequation
            rcollection%IquickAccess(3) = rbc%p_rphysics%creferenceProblem
            rcollection%DquickAccess(1) = dtime
            rcollection%DquickAccess(2) = rbc%p_rphysics%doptControlAlpha
            rcollection%DquickAccess(3) = rbc%p_rphysics%doptControlGamma
            call bcasm_newDirichletBConRealBd (rbc%p_rspaceDiscr, &
                icomponent, rregion, rdiscreteBC, cb_getBoundaryValuesOptC, rcollection)
          end if
        end do
      end do

    case (2)
    
      ! 1D Heat equation
    
      select case (rbc%p_rphysics%creferenceProblem)
      case (0)
      case (1)
        ! 1.)
        if (icomponent .eq. 1) then
          call bcasm_newDirichletBC_1D(rbc%p_rspaceDiscr, rdiscreteBC, &
              fct_heatY1 (0.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              fct_heatY1 (1.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              iequation = 1)
        else
          call bcasm_newDirichletBC_1D(rbc%p_rspaceDiscr, rdiscreteBC, &
              fct_heatLambda1 (0.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              fct_heatLambda1 (1.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              iequation = 2)
        end if

      case (2)
        ! 2.)
        if (icomponent .eq. 1) then
          call bcasm_newDirichletBC_1D(rbc%p_rspaceDiscr, rdiscreteBC, &
              fct_heatY2 (0.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              fct_heatY2 (1.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              iequation = 1)
        else
          call bcasm_newDirichletBC_1D(rbc%p_rspaceDiscr, rdiscreteBC, &
              fct_heatLambda2 (0.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha),  &
              fct_heatLambda2 (1.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              iequation = 2)
        end if

      case (3)
        ! 3.)
        if (icomponent .eq. 1) then
          call bcasm_newDirichletBC_1D(rbc%p_rspaceDiscr, rdiscreteBC, &
              fct_heatY3 (0.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              fct_heatY3 (1.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              iequation = 1)
        else
          call bcasm_newDirichletBC_1D(rbc%p_rspaceDiscr, rdiscreteBC, &
              fct_heatLambda3 (0.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              fct_heatLambda3 (1.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              iequation = 2)
        end if
        
      case (4)
        ! 4.)
        if (icomponent .eq. 1) then
          call bcasm_newDirichletBC_1D(rbc%p_rspaceDiscr, rdiscreteBC, &
              fct_heatY4 (0.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              fct_heatY4 (1.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              iequation = 1)
        else
          call bcasm_newDirichletBC_1D(rbc%p_rspaceDiscr, rdiscreteBC, &
              fct_heatLambda4 (0.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              fct_heatLambda4 (1.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              iequation = 2)
        end if
        
      case (5)
        ! 5.)
        if (icomponent .eq. 1) then
          call bcasm_newDirichletBC_1D(rbc%p_rspaceDiscr, rdiscreteBC, &
              fct_heatY5 (0.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              fct_heatY5 (1.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              iequation = 1)
        else
          call bcasm_newDirichletBC_1D(rbc%p_rspaceDiscr, rdiscreteBC, &
              fct_heatLambda5 (0.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              fct_heatLambda5 (1.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              iequation = 2)
        end if
        
      case (6)
        ! 6.)
        if (icomponent .eq. 1) then
          call bcasm_newDirichletBC_1D(rbc%p_rspaceDiscr, rdiscreteBC, &
              fct_heatY6 (0.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              fct_heatY6 (1.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              iequation = 1)
        else
          call bcasm_newDirichletBC_1D(rbc%p_rspaceDiscr, rdiscreteBC, &
              fct_heatLambda6 (0.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              fct_heatLambda6 (1.0_DP,0.0_DP,dtime,rbc%p_rphysics%doptControlAlpha), &
              iequation = 2)
        end if
        
      case default
        call output_line ("Problem not supported in 1D.")
        call sys_halt()
      end select
      
    case default
    
      call output_line ("Problem not supported.")
      call sys_halt()

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
    real(DP) :: dtimePrimal,dtimeDual

    call spop_getPrimalDualTime (rbc%p_rtimeDiscr,istep,dtimePrimal,dtimeDual)

    select case (rbc%p_rphysics%cequation)
    case (0,2)
      ! Heat equation
      select case (ctype)
      case (SPOP_DEFECT,SPOP_SOLUTION)
        call spop_assembleSpaceBCtime (rbc, dtimePrimal, 1, rdiscreteBC)
        call spop_assembleSpaceBCtime (rbc, dtimeDual, 2, rdiscreteBC)
      case (SPOP_RHS)
        call spop_assembleSpaceBCtime (rbc, dtimeDual, 1, rdiscreteBC)
        call spop_assembleSpaceBCtime (rbc, dtimePrimal, 2, rdiscreteBC)
      end select
    
    case (1)
      ! Stokes equation
      select case (ctype)
      case (SPOP_DEFECT,SPOP_SOLUTION)
        call spop_assembleSpaceBCtime (rbc, dtimePrimal, 1, rdiscreteBC)
        call spop_assembleSpaceBCtime (rbc, dtimePrimal, 2, rdiscreteBC)
        call spop_assembleSpaceBCtime (rbc, dtimeDual, 4, rdiscreteBC)
        call spop_assembleSpaceBCtime (rbc, dtimeDual, 5, rdiscreteBC)
      case (SPOP_RHS)
        call spop_assembleSpaceBCtime (rbc, dtimeDual, 1, rdiscreteBC)
        call spop_assembleSpaceBCtime (rbc, dtimeDual, 2, rdiscreteBC)
        call spop_assembleSpaceBCtime (rbc, dtimePrimal, 4, rdiscreteBC)
        call spop_assembleSpaceBCtime (rbc, dtimePrimal, 5, rdiscreteBC)
      end select

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
!        if (istep .eq. 1) then
!          call lsyssc_clearVector (rtempVec%RvectorBlock(1))
!        end if
!        if (istep .eq. rvector%NEQtime) then
!          call lsyssc_clearVector (rtempVec%RvectorBlock(2))
!        end if
      case (SPOP_RHS)
        call vecfil_discreteBCrhs (rtempVec,rdiscreteBC)
      case (SPOP_SOLUTION)
        call vecfil_discreteBCsol (rtempVec,rdiscreteBC)

        ! DEBUG!!!
!        if (istep .eq. rvector%NEQtime) then
!          call lsyssc_clearVector (rtempVec%RvectorBlock(2))
!          call anprj_discrDirect (rtempVec%RvectorBlock(2),&
!              ffunctionTermCond)
!        end if
      end select
        
      call sptivec_setTimestepData (rvector, istep, rtempVec)
    
    end do
    
    ! Release
    call lsysbl_releaseVector(rtempVec)
    call bcasm_releaseDiscreteBC(rdiscreteBC)
  
  end subroutine

end module

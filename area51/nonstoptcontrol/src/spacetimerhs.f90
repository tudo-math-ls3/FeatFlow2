!##############################################################################
!# ****************************************************************************
!# <name> spacetimerhs </name>
!# ****************************************************************************
!#
!# <purpose>
!# Routines ro assemble the RHS in space and space/time.
!# </purpose>
!##############################################################################

module spacetimerhs

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

  use stdoperators
  use scalarpde
  use linearformevaluation
  use domainintegration
  use timediscretisation
  
  use physics
  use spacetimevectors
  use callback
  
  implicit none

contains

  ! ***************************************************************************

  subroutine strhs_assembleSpaceRHS (rphysics,dtime,icomponent,rrhs)
  
  ! Assembles one component of the RHS at a special point in time.

  ! Physics of the problem
  type(t_physics), intent(in), target :: rphysics

  ! Time
  real(DP), intent(in) :: dtime
  
  ! Component to assemble
  integer, intent(in) :: icomponent
  
  ! Destination vector in space
  type(t_vectorScalar), intent(inout) :: rrhs
  
    ! local variables
    type(t_collection) :: rcollection
    type(t_linearForm) :: rlinform
    
    ! IquickAccess(1) = the component
    ! DquickAccess(1) = current time
    rcollection%IquickAccess(1) = icomponent
    rcollection%DquickAccess(1) = dtime
    rcollection%DquickAccess(2) = rphysics%doptControlAlpha
    
    select case (rphysics%cequation)
    case (0)
      ! Heat equation, one primal and one dual component.
      ! icomponent decides upon that.
     
      rlinform%itermcount = 1
      rlinform%Dcoefficients(1) = 1.0_DP
      rlinform%Idescriptors(1) = 1
     
      call linf_buildVectorScalar2 (rlinform, .true., rrhs,&
          rhs_heatEquation, rcollection)
    
    end select

  end subroutine

  ! ***************************************************************************

  subroutine strhs_assembleRHS (rphysics,rrhs)
  
  ! Assembles the space-time RHS.

  ! Physics of the problem
  type(t_physics), intent(in), target :: rphysics

  ! Destination vector in space
  type(t_spacetimeVector), intent(inout) :: rrhs
  
    ! local variables
    integer :: istep
    real(DP) :: dtimePrimal,dtimeDual,dtheta
    real(DP) :: dtimeend,dtstep,dtimestart
    type(t_vectorBlock) :: rrhsSpace
  
    ! Clear the RHS.
    call sptivec_clearVector (rrhs)
    
    ! Temp vector
    call lsysbl_createVectorBlock (rrhs%p_rspaceDiscr,rrhsSpace)
    
    ! Time discretisation scheme
    select case (rrhs%p_rtimeDiscr%ctype)
    case (TDISCR_ONESTEPTHETA)
    
      dtheta = rrhs%p_rtimeDiscr%dtheta
    
      ! Loop through the timesteps
      do istep = 1,rrhs%NEQtime
      
        ! Get the time of the primal and dual equation
        call tdiscr_getTimestep(rrhs%p_rtimeDiscr,istep-1,&
            dtimeend,dtstep,dtimestart)
            
        ! Primal equation lives in the endpoints of the time interval.
        ! Dual equation lives 'inbetween'.
        dtimeprimal = dtheta*dtimeend + (1-dtheta)*dtimestart
        dtimedual = dtimeend
    
        call lsysbl_clearVector (rrhsSpace)
    
        select case (rphysics%cequation)
        case (0)
          ! Heat equation, one primal and one dual component.
          call strhs_assembleSpaceRHS (rphysics,dtimeprimal,1,&
              rrhsSpace%RvectorBlock(1))

          call strhs_assembleSpaceRHS (rphysics,dtimedual,2,&
              rrhsSpace%RvectorBlock(2))
        
        end select
        
        call sptivec_setTimestepData (rrhs, istep, rrhsSpace)
        
      end do
      
    case default
      ! Should not happen
      call sys_halt()
      
    end select

    call lsysbl_releaseVector (rrhsSpace)

  end subroutine

end module

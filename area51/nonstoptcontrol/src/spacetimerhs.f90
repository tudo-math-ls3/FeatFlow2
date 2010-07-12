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

  subroutine strhs_assembleSpaceRHS (rphysics,bclear,dtime,dtstep,&
      ithetaschemetype,dtheta,icomponent,rrhs)
  
  ! Assembles one component of the RHS at a special point in time.

  ! Physics of the problem
  type(t_physics), intent(in), target :: rphysics

  ! Time
  real(DP), intent(in) :: dtime
  
  ! Clear old RHS vector.
  logical, intent(in) :: bclear
  
  ! Timestep
  real(DP), intent(in) :: dtstep

  ! type of the theta scheme
  integer, intent(in) :: ithetaschemetype

  ! Theta scheme identifier
  real(DP), intent(in) :: dtheta
  
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
    rcollection%IquickAccess(2) = rphysics%creferenceProblem
    rcollection%IquickAccess(3) = ithetaschemetype
    rcollection%DquickAccess(1) = dtime
    rcollection%DquickAccess(2) = dtstep
    rcollection%DquickAccess(3) = rphysics%doptControlAlpha
    rcollection%DquickAccess(4) = rphysics%doptControlGamma
    rcollection%DquickAccess(5) = dtheta
    
    select case (rphysics%cequation)
    case (0,2)
      ! Heat equation, one primal and one dual component.
      ! icomponent decides upon that.
     
      rlinform%itermcount = 1
      rlinform%Dcoefficients(1) = 1.0_DP
      rlinform%Idescriptors(1) = 1
     
      call linf_buildVectorScalar2 (rlinform, bclear, rrhs,&
          rhs_heatEquation, rcollection)

    case (1)
      ! Stokes equation.
     
      rlinform%itermcount = 1
      rlinform%Dcoefficients(1) = 1.0_DP
      rlinform%Idescriptors(1) = 1
     
      call linf_buildVectorScalar2 (rlinform, bclear, rrhs,&
          rhs_stokesEquation, rcollection)
    
    case default
    
      call output_line ("Equation not supported.")
      call sys_halt()

    end select
    
    ! DEBUG: Scale the RHS by dtstep.
    !call lsyssc_scaleVector (rrhs,dtstep)

  end subroutine

  ! ***************************************************************************

  subroutine strhs_assembleRHS (rphysics,rrhs)
  
  ! Assembles the space-time RHS.

  ! Physics of the problem
  type(t_physics), intent(in), target :: rphysics

  ! Destination vector in space
  type(t_spacetimeVector), intent(inout) :: rrhs
  
    ! local variables
    integer :: istep,ithetaschemetype
    real(DP) :: dtimePrimal,dtimeDual,dtheta
    real(DP) :: dtimeend,dtstep,dtimestart
    type(t_vectorBlock) :: rrhsSpace
    real(DP), dimension(:), pointer :: p_Ddata
    
    ithetaschemetype = rrhs%p_rtimeDiscr%itag
  
    ! Clear the RHS.
    call sptivec_clearVector (rrhs)
    
    ! Temp vector
    call lsysbl_createVectorBlock (rrhs%p_rspaceDiscr,rrhsSpace)
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rrhsSpace, p_Ddata)
    
    ! Time discretisation scheme
    select case (rrhs%p_rtimeDiscr%ctype)
    case (TDISCR_ONESTEPTHETA)
    
      dtheta = rrhs%p_rtimeDiscr%dtheta
    
      ! Loop through the timesteps
      do istep = 1,rrhs%NEQtime
      
        ! Get the time of the primal and dual equation
        call tdiscr_getTimestep(rrhs%p_rtimeDiscr,istep-1,&
            dtimeend,dtstep,dtimestart)
        
        if (ithetaschemetype .eq. 1) then
          ! Dual rhs lives in the endpoints of the time interval.
          ! Primal equation lives 'inbetween'.
          dtimeprimal = dtheta*dtimeend + (1-dtheta)*dtimestart
          dtimedual = dtimeend
        else
          dtimeprimal = dtimeend
          dtimedual = dtimeend
        end if
    
        call lsysbl_clearVector (rrhsSpace)
    
        select case (rphysics%cequation)
        case (0,2)
          ! Heat equation, one primal and one dual component.
          call strhs_assembleSpaceRHS (rphysics,.true.,dtimeprimal,dtstep,&
              ithetaschemetype,dtheta,1,rrhsSpace%RvectorBlock(1))

          call strhs_assembleSpaceRHS (rphysics,.true.,dtimedual,dtstep,&
              ithetaschemetype,dtheta,2,rrhsSpace%RvectorBlock(2))
              
          if (ithetaschemetype .eq. 1) then
            if (istep .eq. 1) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(2),&
                  1.0_DP-dtheta)
            end if

            if (istep .eq. rrhs%NEQtime) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(2),&
                  (dtheta + dtheta*rphysics%doptControlGamma/dtstep))
            end if
            
          else if (ithetaschemetype .eq. 3) then
            if (istep .eq. rrhs%NEQtime) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(2),&
                  (dtheta + rphysics%doptControlGamma/dtstep))
            end if

            if (istep .lt. rrhs%NEQtime) then
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(2),&
                  dtheta/(1.0_DP-dtheta))

              call strhs_assembleSpaceRHS (rphysics,.false.,dtimedual+dtstep,dtstep,&
                  ithetaschemetype,dtheta,2,rrhsSpace%RvectorBlock(2))
            
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(2),&
                  (1.0_DP-dtheta))
            end if
            
          else
            if (istep .eq. rrhs%NEQtime) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(2),&
                  (1.0_DP + rphysics%doptControlGamma/dtstep))
            end if
          end if

          ! DEBUG!!!
          !call lsysbl_scaleVector (rrhsSpace,rrhs%p_rtimeDiscr%dtstep)
        
        case (1)
          ! Stokes equation.
          call strhs_assembleSpaceRHS (rphysics,.true.,dtimeprimal,dtstep,&
              ithetaschemetype,dtheta,1,rrhsSpace%RvectorBlock(1))

          call strhs_assembleSpaceRHS (rphysics,.true.,dtimedual,dtstep,&
              ithetaschemetype,dtheta,2,rrhsSpace%RvectorBlock(2))

          call strhs_assembleSpaceRHS (rphysics,.true.,dtimedual,dtstep,&
              ithetaschemetype,dtheta,4,rrhsSpace%RvectorBlock(4))

          call strhs_assembleSpaceRHS (rphysics,.true.,dtimedual,dtstep,&
              ithetaschemetype,dtheta,5,rrhsSpace%RvectorBlock(5))

          if (ithetaschemetype .eq. 1) then
            if (istep .eq. 1) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(4),&
                  1.0_DP-dtheta)
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(5),&
                  1.0_DP-dtheta)
            end if
            
            if (istep .eq. rrhs%NEQtime) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(4),&
                  (dtheta + dtheta*rphysics%doptControlGamma/dtstep))

              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(5),&
                  (dtheta + dtheta*rphysics%doptControlGamma/dtstep))
            end if
            
          else

            if (istep .eq. rrhs%NEQtime) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(4),&
                  (1.0_DP + rphysics%doptControlGamma/dtstep))

              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(5),&
                  (1.0_DP + rphysics%doptControlGamma/dtstep))
            end if
            
          end if

        case default
        
          call output_line ("Equation not supported.")
          call sys_halt()

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

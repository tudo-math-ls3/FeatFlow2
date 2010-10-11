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
  ! The assembled RHS is unweighted with respect to the THETA scheme!

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

  ! Theta scheme weight
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
    rcollection%DquickAccess(6) = rphysics%dpar
    rcollection%DquickAccess(7) = rphysics%dcouplePrimalToDual
    rcollection%DquickAccess(8) = rphysics%dcoupleDualToPrimal
    rcollection%DquickAccess(9) = rphysics%dtimeMin
    rcollection%DquickAccess(10) = rphysics%dtimeMax
    
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

  subroutine strhs_assembleSpaceSolRHS (rphysics,bclear,dtime,dtstep,&
      ithetaschemetype,dtheta,icomponent,rrhs)
  
  ! Assembles one component of the solution into the RHS at a special point in time.
  ! The assembled solution-RHS is unweighted with respect to the THETA scheme!

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

  ! Theta scheme weight
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
    rcollection%DquickAccess(6) = rphysics%dpar
    rcollection%DquickAccess(7) = rphysics%dcouplePrimalToDual
    rcollection%DquickAccess(8) = rphysics%dcoupleDualToPrimal
    rcollection%DquickAccess(9) = rphysics%dtimeMin
    rcollection%DquickAccess(10) = rphysics%dtimeMax
    
    select case (rphysics%cequation)
    case (0,2)
      ! Heat equation, one primal and one dual component.
      ! icomponent decides upon that.
     
      rlinform%itermcount = 1
      rlinform%Dcoefficients(1) = 1.0_DP
      rlinform%Idescriptors(1) = 1
     
      call linf_buildVectorScalar2 (rlinform, bclear, rrhs,&
          sol_heatEquation, rcollection)

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
    real(DP) :: dtheta,dtstep
    type(t_vectorBlock) :: rrhsSpace
    type(t_vectorBlock), pointer :: p_rrhsSol,p_rrhsTemp
    type(t_spaceTimeVectorAccess) :: raccessPoolSol,raccessPoolRhs
    real(DP), dimension(:), pointer :: p_Ddata, p_DdataTemp
    
    ithetaschemetype = rrhs%p_rtimeDiscr%itag
  
    ! Clear the RHS.
    call sptivec_clearVector (rrhs)
    
    ! Temp vectors
    call lsysbl_createVectorBlock (rrhs%p_rspaceDiscr,rrhsSpace)
    
    call sptivec_createAccessPool (rrhs%p_rspaceDiscr,raccessPoolRhs,3)
    call sptivec_createAccessPool (rrhs%p_rspaceDiscr,raccessPoolSol,3)
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rrhsSpace, p_Ddata)
    call lsysbl_getbase_double (rrhsSpace, p_DdataTemp)
    
    ! Time discretisation scheme
    select case (rrhs%p_rtimeDiscr%ctype)
    case (TDISCR_ONESTEPTHETA)
    
      dtheta = rrhs%p_rtimeDiscr%dtheta
      dtstep = rrhs%p_rtimeDiscr%dtstep
    
      ! Loop through the timesteps
      do istep = 1,rrhs%NEQtime
      
        ! Create the 'current' rhs.
        call getrhs (rphysics,rrhs%p_rtimediscr,raccesspoolRhs,istep,0,p_rrhsTemp)

        select case (rphysics%cequation)
        case (0,2)

          ! Do some final modifications to the vector, depending on
          ! - the timestep scheme
          ! - the coupling between the primal and dual solutions
          
          select case (ithetaschemetype)
          case (1)
            ! Take the 'current' RHS.
            call lsysbl_copyVector (p_rrhsTemp,rrhsSpace)
          
            if (istep .eq. 1) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(2),&
                  dtheta)
            end if

            if (istep .eq. rrhs%NEQtime-1) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(2),&
                  (1.0_DP + (1.0_DP-dtheta)*rphysics%doptControlGamma/dtstep))
            end if

            if (istep .eq. rrhs%NEQtime) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(2),&
                  (dtheta + dtheta*rphysics%doptControlGamma/dtstep))
            end if
            
          case default
          
            ! Take the mean of the current and last RHS except for if
            ! we are in the 1st timestep.
            call lsysbl_copyVector (p_rrhsTemp,rrhsSpace)
            
            if (istep .eq. rrhs%NEQtime) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(2),&
                  (1.0_DP + rphysics%doptControlGamma/dtstep))
            end if

            ! Forward equation
            if (istep .gt. 1) then
              call getrhs (rphysics,rrhs%p_rtimediscr,raccesspoolRhs,istep-1,0,p_rrhsTemp)
              call lsyssc_vectorLinearComb (p_rrhsTemp%RvectorBlock(1),rrhsSpace%RvectorBlock(1),&
                  1.0_DP-dtheta,dtheta)
            end if

            ! Backward equation
            if (istep .lt. rrhs%NEQtime) then
              call getrhs (rphysics,rrhs%p_rtimediscr,raccesspoolRhs,istep+1,0,p_rrhsTemp)
              call lsyssc_vectorLinearComb (p_rrhsTemp%RvectorBlock(2),rrhsSpace%RvectorBlock(2),&
                  1.0_DP-dtheta,dtheta)
            end if
            
          end select
          
          ! Process the coupling -- or decoupling -- of the equations.
          !

          if ((rphysics%dcoupleDualToPrimal .ne. 1.0_DP) .or. &
              (rphysics%dcouplePrimalToDual .ne. 1.0_DP)) then

            ! Dual decoupling
            select case (ithetaschemetype)
            case (1)

              ! Create the 'current' solution rhs.
              call getrhs (rphysics,rrhs%p_rtimediscr,raccessPoolSol,istep,1,p_rrhsSol)
            
              ! Primal decoupling
              call lsyssc_vectorLinearComb (p_rrhsSol%RvectorBlock(2),&
                  rrhsSpace%RvectorBlock(1),&
                  (1.0_DP-rphysics%dcoupleDualToPrimal)*(-1.0_DP/rphysics%doptControlAlpha),1.0_DP)

              if (istep .eq. 1) then
              
                call lsyssc_vectorLinearComb (p_rrhsSol%RvectorBlock(1),&
                    rrhsSpace%RvectorBlock(2),&
                    (1.0_DP-rphysics%dcouplePrimalToDual)*(1.0_DP-dtheta),1.0_DP)
              
              else if (istep .eq. rrhs%NEQtime-1) then

                call lsyssc_vectorLinearComb (p_rrhsSol%RvectorBlock(1),&
                    rrhsSpace%RvectorBlock(2),&
                    (1.0_DP-rphysics%dcouplePrimalToDual)*&
                      (1.0_DP+(1.0_DP-dtheta)*rphysics%doptControlGamma/dtstep),1.0_DP)

              else if (istep .eq. rrhs%NEQtime) then

                call lsyssc_vectorLinearComb (p_rrhsSol%RvectorBlock(1),&
                    rrhsSpace%RvectorBlock(2),&
                    (1.0_DP-rphysics%dcouplePrimalToDual)*&
                      (rphysics%dcoupleTermCond*dtheta+dtheta*rphysics%doptControlGamma/dtstep),1.0_DP)

              else
              
                ! Standard weights.
                call lsyssc_vectorLinearComb (p_rrhsSol%RvectorBlock(1),&
                    rrhsSpace%RvectorBlock(2),&
                    (1.0_DP-rphysics%dcouplePrimalToDual),1.0_DP)
                    
              end if
            
            case default

              ! Primal decoupling
              if (istep .gt. 1) then
                ! Create the 'previous' solution rhs.
                call getrhs (rphysics,rrhs%p_rtimediscr,raccessPoolSol,istep-1,1,p_rrhsSol)

                call lsyssc_vectorLinearComb (p_rrhsSol%RvectorBlock(2),&
                    rrhsSpace%RvectorBlock(1),&
                    (1.0_DP-dtheta)*(1.0_DP-rphysics%dcoupleDualToPrimal)*(-1.0_DP/rphysics%doptControlAlpha),1.0_DP)

                ! Create the 'current' solution rhs.
                call getrhs (rphysics,rrhs%p_rtimediscr,raccessPoolSol,istep,1,p_rrhsSol)

                call lsyssc_vectorLinearComb (p_rrhsSol%RvectorBlock(2),&
                    rrhsSpace%RvectorBlock(1),&
                    dtheta*(1.0_DP-rphysics%dcoupleDualToPrimal)*(-1.0_DP/rphysics%doptControlAlpha),1.0_DP)
              end if

              ! Dual decoupling
              if (istep .lt. rrhs%NEQtime) then
              
                ! Create the 'previous' solution rhs.
                call getrhs (rphysics,rrhs%p_rtimediscr,raccessPoolSol,istep+1,1,p_rrhsSol)
              
                call lsyssc_vectorLinearComb (p_rrhsSol%RvectorBlock(1),&
                    rrhsSpace%RvectorBlock(2),(1.0_DP-dtheta)*(1.0_DP-rphysics%dcouplePrimalToDual),1.0_DP)

                ! Create the 'current' solution rhs.
                call getrhs (rphysics,rrhs%p_rtimediscr,raccessPoolSol,istep,1,p_rrhsSol)

                call lsyssc_vectorLinearComb (p_rrhsSol%RvectorBlock(1),&
                    rrhsSpace%RvectorBlock(2),dtheta*(1.0_DP-rphysics%dcouplePrimalToDual),1.0_DP)
                    
              else

                ! Create the 'current' solution rhs.
                call getrhs (rphysics,rrhs%p_rtimediscr,raccessPoolSol,istep,1,p_rrhsSol)
              
                call lsyssc_vectorLinearComb (p_rrhsSol%RvectorBlock(1),&
                    rrhsSpace%RvectorBlock(2),&
                    (1.0_DP-rphysics%dcouplePrimalToDual)*&
                        (rphysics%dcoupleTermCond+rphysics%doptControlGamma/dtstep),1.0_DP)
              end if
              
            end select
            
          end if

          ! DEBUG!!!
          !call lsysbl_scaleVector (rrhsSpace,rrhs%p_rtimeDiscr%dtstep)
        
        case (1)
          ! Stokes equation.

          ! Take the 'current' RHS.
          call lsysbl_copyVector (p_rrhsTemp,rrhsSpace)

          if (ithetaschemetype .eq. 1) then
            if (istep .eq. 1) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(4),&
                  1.0_DP-dtheta)
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(5),&
                  1.0_DP-dtheta)
            end if
            
            if (istep .eq. rrhs%NEQtime-1) then
              ! Scale the RHS according to the timestep scheme.
              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(4),&
                  (1.0_DP + (1.0_DP-dtheta)*rphysics%doptControlGamma/dtstep))

              call lsyssc_scaleVector (rrhsSpace%RvectorBlock(5),&
                  (1.0_DP + (1.0_DP-dtheta)*rphysics%doptControlGamma/dtstep))
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

    call sptivec_releaseAccessPool (raccessPoolSol)
    call sptivec_releaseAccessPool (raccessPoolRhs)
    call lsysbl_releaseVector (rrhsSpace)
    
  contains
  
    subroutine getrhs (rphysics,rtimediscr,raccesspool,itimestep,itype,p_rvector)
    
    ! Physics of the problem
    type(t_physics), intent(in), target :: rphysics

    ! Time discretisation
    type(t_timeDiscretisation), intent(in) :: rtimediscr

    ! Vector pool    
    type(t_spaceTimeVectorAccess), intent(inout) :: raccesspool
    
    ! Timestep
    integer, intent(in) :: itimestep
    
    ! Type of RHS vector.
    ! =0: f/z, =1: y/lambda
    integer, intent(in) :: itype
    
    ! Pointer that points to the vector
    type(t_vectorBlock), pointer :: p_rvector
    
      ! local variables
      real(DP) :: dtimeend,dtstep,dtimestart,dtimeprimal,dtimedual
      real(DP) :: dtheta
      integer :: ithetaschemetype
      
      dtheta = rtimediscr%dtheta
      ithetaschemetype = rtimediscr%itag
    
      ! Try to fetch the RHS.
      call sptivec_getVectorFromPool (raccessPool,itimestep,p_rvector)
      
      if (.not. associated(p_rvector)) then
        
        ! Create it.
        call sptivec_getFreeBufferFromPool (raccessPool,itimestep,p_rvector)
        
        ! Get the time of the primal and dual equation
        call tdiscr_getTimestep(rtimediscr,itimestep-1,&
            dtimeend,dtstep,dtimestart)
        
        if (ithetaschemetype .eq. 1) then
          ! Dual rhs lives in the endpoints of the time interval.
          ! Primal rhs lives 'inbetween'.
          dtimeprimal = dtheta*dtimeend + (1-dtheta)*dtimestart
          dtimedual = dtimeend
        else
          dtimeprimal = dtimeend
          dtimedual = dtimeend
        end if
        
        select case (rphysics%cequation)
        case (0,2)
          ! Heat equation, one primal and one dual component.
          select case (itype)
          case (0)
            call strhs_assembleSpaceRHS (rphysics,.true.,dtimeprimal,dtstep,&
                ithetaschemetype,dtheta,1,p_rvector%RvectorBlock(1))

            call strhs_assembleSpaceRHS (rphysics,.true.,dtimedual,dtstep,&
                ithetaschemetype,dtheta,2,p_rvector%RvectorBlock(2))
                
          case (1)

            ! Assemble the solution as RHS
            call strhs_assembleSpaceSolRHS (rphysics,.true.,dtimedual,dtstep,&
                ithetaschemetype,dtheta,1,p_rvector%RvectorBlock(1))

            call strhs_assembleSpaceSolRHS (rphysics,.true.,dtimeprimal,dtstep,&
                ithetaschemetype,dtheta,2,p_rvector%RvectorBlock(2))
                
          end select
              
        case (1)

          ! Stokes equation.
          select case (itype)
          case (0)
            call strhs_assembleSpaceRHS (rphysics,.true.,dtimeprimal,dtstep,&
                ithetaschemetype,dtheta,1,p_rvector%RvectorBlock(1))

            call strhs_assembleSpaceRHS (rphysics,.true.,dtimeprimal,dtstep,&
                ithetaschemetype,dtheta,2,p_rvector%RvectorBlock(2))

            call strhs_assembleSpaceRHS (rphysics,.true.,dtimedual,dtstep,&
                ithetaschemetype,dtheta,3,p_rvector%RvectorBlock(3))

            call strhs_assembleSpaceRHS (rphysics,.true.,dtimedual,dtstep,&
                ithetaschemetype,dtheta,4,p_rvector%RvectorBlock(4))

            call strhs_assembleSpaceRHS (rphysics,.true.,dtimedual,dtstep,&
                ithetaschemetype,dtheta,5,p_rvector%RvectorBlock(5))

            call strhs_assembleSpaceRHS (rphysics,.true.,dtimeprimal,dtstep,&
                ithetaschemetype,dtheta,6,p_rvector%RvectorBlock(6))
            
          case (1)
            call strhs_assembleSpaceSolRHS (rphysics,.true.,dtimedual,dtstep,&
                ithetaschemetype,dtheta,1,p_rvector%RvectorBlock(1))

            call strhs_assembleSpaceSolRHS (rphysics,.true.,dtimedual,dtstep,&
                ithetaschemetype,dtheta,2,p_rvector%RvectorBlock(2))

            call strhs_assembleSpaceSolRHS (rphysics,.true.,dtimeprimal,dtstep,&
                ithetaschemetype,dtheta,3,p_rvector%RvectorBlock(3))

            call strhs_assembleSpaceSolRHS (rphysics,.true.,dtimeprimal,dtstep,&
                ithetaschemetype,dtheta,4,p_rvector%RvectorBlock(4))

            call strhs_assembleSpaceSolRHS (rphysics,.true.,dtimeprimal,dtstep,&
                ithetaschemetype,dtheta,5,p_rvector%RvectorBlock(5))

            call strhs_assembleSpaceSolRHS (rphysics,.true.,dtimedual,dtstep,&
                ithetaschemetype,dtheta,6,p_rvector%RvectorBlock(6))
  
          end select
          
        end select
        
      end if
      
    end subroutine

  end subroutine

end module

!##############################################################################
!# ****************************************************************************
!# <name> optcanalysis </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to perform error analysis for the solution
!# of the optimal control problem. Here, routines can be found that calculate
!# the functional, which is to be minimised by the theory:
!#
!# $$ min J(y,u) = 1/2||y-z||_{L^2} + \gamma/2||y(T)-z(T)||_{L^2} + \alpha/2||u||^2 $$
!#
!# The following routines can be found here:
!#
!# 1.) optcana_stationaryFunctional
!#     -> Calculates the value of the stationary functional which is to be
!#        minimised:
!#     $$ J(y,u) = 1/2||y-z||_{L^2} + \alpha/2||u||^2 $$
!#
!# 2.) optcana_nonstatFunctional
!#     -> Calculates the value of the stationary functional which is to be
!#        minimised:
!#     $$ J(y,u) = 1/2||y-z||^2_{L^2} + \alpha/2||u||^2_{L^2} + gamma/2||y(T)-z(T)||^2_{L^2}$$
!#
!# 3.) optcana_analyticalError
!#     -> Calculates the error ||y-y0|| of a given function y0 to an
!#        analytically given reference function y0.
!# </purpose>
!##############################################################################

module optcanalysis

  use fsystem
  use genoutput
  use basicgeometry
  use linearsystemscalar
  use linearsystemblock
  use pprocerror
  use collection
  use domainintegration
  use spatialdiscretisation
  use fparser
  use pprocerror
  
  use analyticsolution
  
  use timediscretisation
  use timeevaluation
  use spacetimevectors
  use user_callback
  use spacematvecassembly
  use spacetimelinearsystem
  
  use structuresoptc
    
  implicit none
  
  private
  
  !public :: optcana_stationaryFunctional
  public :: optcana_nonstatFunctional
  public :: optcana_analyticalError
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine optcana_evalFunction (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  ! Standard evaluation routine. Evaluates component rcollection%Icollection(1)
  ! of the analytical solution identified by the name "SOL" in the collection.

  integer, intent(IN)                                         :: cderivative
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  integer, intent(IN)                                         :: nelements
  integer, intent(IN)                                         :: npointsPerElement
  real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
  integer, dimension(:,:), intent(IN) :: IdofsTest
  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset
  type(t_collection), intent(INOUT), optional      :: rcollection
  real(DP), dimension(:,:), intent(out) :: Dvalues
  
  integer :: ierror,ctype,icomponent
  
    ! Get data
    icomponent = rcollection%IquickAccess(1)
    ctype = rcollection%IquickAccess(2)
  
    if (ctype .ne. ANSOL_TP_ANALYTICAL) then
      ! Evaluate the reference using precalculated information in the collection
      ! from an analytical function.
      call ansol_evaluate (rcollection,"SOL",icomponent,&
          Dvalues,npointsPerElement,nelements,Dpoints,rdomainIntSubset%p_Ielements,ierror)
    else
      ! This is a function realised by our callback routines.
      ! Call them to get the information. Pass the "connected" collection
      ! which contains postprocessing data.
      call user_coeff_Reference (cderivative,rdiscretisation, &
          nelements,npointsPerElement,Dpoints, &
          IdofsTest,rdomainIntSubset,&
          Dvalues,rcollection%p_rnextCollection)
      ierror = 0
    end if
        
    ! Check that this was ok. If yes, copy the data to the destination.
    if (ierror .ne. 0) then
      call output_line ('Error evaluating RHS function.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'optcana_evalFunction')
      call sys_halt()
    end if

  end subroutine

!!******************************************************************************
!
!!<subroutine>
!
!  subroutine optcana_stationaryFunctional (rglobalData,rsolution,rreference,dalpha,Derror)
!
!!<description>
!  ! This function calculates the value of the functional which is to be
!  ! minimised in the stationary optimal control problem.
!  ! The functional is defined as
!  !   $$ J(y,u) = 1/2||y-z||_{L^2}^2  + \alpha/2||u||^2 $$
!  ! over a spatial domain $\Omega$.
!  !
!  ! For this purpose, the routine evaluates the user-defined callback functions
!  ! user_ffunction_TargetX and user_ffunction_TargetY. The collection must be initialised
!  ! for postprocessing before calling his routine!
!!</description>
!
!!<input>
!  ! Global settings for callback routines.
!  type(t_globalData), intent(inout), target :: rglobalData
!
!  ! Solution vector to compute the norm/error from.
!  type(t_vectorBlock), intent(IN) :: rsolution
!
!  ! Analytic solution defining the reference function z.
!  type(t_anSolution), intent(inout) :: rreference
!
!  ! Regularisation parameter $\alpha$.
!  real(DP), intent(IN) :: dalpha
!!</input>
!
!!<output>
!  ! Returns information about the error.
!  ! Derror(1) = ||y-z||_{L^2}.
!  ! Derror(2) = ||u||.
!  ! Derror(3) = J(y,u) = 1/2||y-z||_{L^2}^2  + \alpha/2||u||^2.
!  ! Norm of the error functional.
!  real(DP), dimension(:), intent(OUT) :: Derror
!!</output>
!
!!</subroutine>
!
!    ! local variables
!    real(DP),dimension(2) :: Derr
!    type(t_collection) :: rcollection
!
!    ! Initialise the collection for the assembly process with callback routines.
!    ! This stores the simulation time in the collection and sets the
!    ! current subvector z for the callback routines.
!    call collct_init(rcollection)
!
!    ! Perform error analysis to calculate and add 1/2||y-z||_{L^2}.
!    if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
!
!      ! Perform error analysis to calculate and add 1/2||y-z||^2_{L^2}.
!      call user_initCollectForAssembly (rglobalData,0.0_DP,rcollection)
!
!      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
!                        user_ffunction_TargetX,rcollection)
!
!      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
!                        user_ffunction_TargetY,rcollection)
!
!      call user_doneCollectForAssembly (rglobalData,rcollection)
!
!    else
!
!      ! Use our standard implementation to evaluate the functional.
!      call ansol_prepareEval (rreference,rcollection,"SOL",0.0_DP)
!
!      ! X-velocity
!      rcollection%IquickAccess(1) = 1
!      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
!          optcana_evalFunction,rcollection)
!
!      ! Y-velocity
!      rcollection%IquickAccess(1) = 2
!      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
!          optcana_evalFunction,rcollection)
!
!      call ansol_doneEval (rcollection,"SOL")
!
!    end if
!
!    Derror(1) = (0.5_DP*(Derr(1)**2+Derr(2)**2))
!    Derror(2) = 0.0_DP
!
!    ! Calculate \alpha/2||u||^2.
!    if (dalpha .ne. 0.0_DP) then
!      call pperr_scalar (rsolution%RvectorBlock(4),PPERR_L2ERROR,Derr(1))
!
!      call pperr_scalar (rsolution%RvectorBlock(5),PPERR_L2ERROR,Derr(2))
!
!      Derror(2) = (0.5_DP*(Derr(1)**2+Derr(2)**2))
!
!      ! Because of u=-lambda/alpha we have:
!      !    alpha/2 ||u||^2 = alpha/2 ||lambda/alpha||^2 = 1/(2*alpha) ||lambda||^2
!    end if
!
!    Derror(3) = 0.5_DP * Derror(1)+(0.5_DP/dalpha) * Derror(2)
!    Derror(1) = sqrt(Derror(1))
!    Derror(2) = sqrt(Derror(2))
!
!    ! Clean up
!    call collct_done(rcollection)
!
!  end subroutine

!******************************************************************************

!<subroutine>

  subroutine optcana_nonstatFunctional (rglobalData,rphysics,rconstraints,rsolution,rreference,&
      dalpha,dgamma,Derror)

!<description>
  ! This function calculates the value of the functional which is to be
  ! minimised in the stationary optimal control problem.
  ! The functional is defined as
  !   $$ J(y,u) = 1/2||y-z||^2_{L^2} + \alpha/2||u||^2_{L^2} + gamma/2||y(T)-z(T)||^2_{L^2}$$
  ! over a spatial domain $\Omega$.
  !
  ! For this purpose, the routine evaluates the user-defined callback functions
  ! user_ffunction_TargetX and user_ffunction_TargetY. The collection must be initialised
  ! for postprocessing before calling his routine!
  !
  ! The function z is given implicitely in the problem structure rproblem
  ! and evaluated in user_ffunction_TargetX and user_ffunction_TargetY!
!</description>
  
!<input>
  ! Physics of the problem.
  type(t_settings_physics) :: rphysics

  ! Solution vector to compute the norm/error from.
  type(t_spacetimeVector), intent(IN) :: rsolution
  
  ! Analytic solution defining the reference function z.
  type(t_anSolution), intent(inout) :: rreference

  ! Constraints in the optimal control problem.
  type(t_optcconstraintsSpaceTime), intent(in) :: rconstraints
  
  ! Regularisation parameter $\alpha$.
  real(DP), intent(IN) :: dalpha

  ! Regularisation parameter $\gamma$.
  real(DP), intent(IN) :: dgamma

  ! Global settings for callback routines.
  type(t_globalData), intent(inout), target :: rglobalData
!</input>

!<output>
  ! Returns information about the error.
  ! Derror(1) = ||y-z||_{L^2}.
  ! Derror(2) = ||u||_{L^2}.
  ! Derror(3) = ||y(T)-z(T)||_{L^2}.
  ! Derror(4) = J(y,u).
  real(DP), dimension(:), intent(OUT) :: Derror
!</output>
  
!</subroutine>
    
    ! local variables
    integer :: isubstep
    real(DP) :: dtstep,dtime
    real(DP),dimension(2) :: Derr
    type(t_collection) :: rcollection
    type(t_vectorBlock) :: rtempVector
    real(dp), dimension(:), pointer :: p_Dx
    type(t_optcconstraintsSpace) :: rconstrSpace
    
    ! Create a temp vector
    call lsysbl_createVectorBlock(rsolution%p_rspaceDiscr,rtempVector,.true.)
    
    call lsysbl_getbase_double (rtempVector,p_Dx)
    
    ! Initialise the collection for the assembly process with callback routines.
    ! This stores the simulation time in the collection and sets the
    ! current subvector z for the callback routines.
    call collct_init(rcollection)
    
    Derror(1:4) = 0.0_DP

    do isubstep = 1,rsolution%NEQtime
    
      ! Current point in time.
      ! For the first iterate, take a look at time interval 1 to
      ! get the length and start point of the interval.
      ! For all further iterates, look at the time interval leading to
      ! that iterate.
      if (isubstep .gt. 1) then
        call tdiscr_getTimestep(rsolution%p_rtimeDiscr,isubstep-1,dtime,dtstep)
      else
        call tdiscr_getTimestep(rsolution%p_rtimeDiscr,isubstep,dtstep=dtstep,dtimestart=dtime)
      end if

      ! Get the solution.
      ! Evaluate the space time function in rvector in the point
      ! in time dtime. Independent of the discretisation in time,
      ! this will give us a vector in space.
      !CALL sptivec_getTimestepData (rsolution, 1+isubstep, rtempVector)
      call tmevl_evaluate(rsolution,dtime,rtempVector)

      select case (rphysics%cequation)
      case (0,1)
        ! Stokes, Navier-Stokes, 2D
        
        ! Compute:
        ! Derror(1) = ||y-z||^2_{L^2}.
        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          
          ! Perform error analysis to calculate and add 1/2||y-z||^2_{L^2}.
          call user_initCollectForVecAssembly (rglobalData,&
              rreference%iid,1,dtime,rcollection)

          call pperr_scalar (rtempVector%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
              user_fct_Target,rcollection)

          call user_doneCollectForAssembly (rglobalData,rcollection)

          call user_initCollectForVecAssembly (rglobalData,&
              rreference%iid,2,dtime,rcollection)

          call pperr_scalar (rtempVector%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
              user_fct_Target,rcollection)
              
          call user_doneCollectForAssembly (rglobalData,rcollection)
        else
          ! Use our standard implementation to evaluate the functional.
          call ansol_prepareEval (rreference,rcollection,"SOL",dtime)

          ! X-velocity
          rcollection%IquickAccess(1) = 1
          call pperr_scalar (rtempVector%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
              optcana_evalFunction,rcollection)

          ! Y-velocity
          rcollection%IquickAccess(1) = 2
          call pperr_scalar (rtempVector%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
              optcana_evalFunction,rcollection)
              
          call ansol_doneEval (rcollection,"SOL")
        end if

        ! We use the summed trapezoidal rule.
        if ((isubstep .eq. 1) .or. (isubstep .eq. rsolution%NEQtime)) then
          Derror(1) = Derror(1) + 0.5_DP*0.5_DP*(Derr(1)**2 + Derr(2)**2) * dtstep
        else
          Derror(1) = Derror(1) + 0.5_DP*(Derr(1)**2 + Derr(2)**2) * dtstep
        end if

        ! Compute:
        ! Derror(3) = ||y(T)-z(T)||^2
        if (isubstep .eq. rsolution%NEQtime) then
          Derror(3) = 0.5_DP*(Derr(1)**2+Derr(2)**2)
        end if
        
        ! Compute:
        ! Derror(2) = ||u|| = ||P[min/max](-1/alpha lambda)||^2_{L^2}.
        ! For that purpose, scale the lambda part and project it if necessary.
        call lsyssc_scaleVector (rtempVector%RvectorBlock(4),-1.0_DP/dalpha)
        call lsyssc_scaleVector (rtempVector%RvectorBlock(5),-1.0_DP/dalpha)
        
        if (rconstraints%ccontrolConstraints .ne. 0) then
          select case (rconstraints%cconstraintsType)
          case (0)
            call smva_projectControlTstepConst (rtempVector%RvectorBlock(4),&
                rconstraints%dumin1,rconstraints%dumax1)
            call smva_projectControlTstepConst (rtempVector%RvectorBlock(5),&
                rconstraints%dumin2,rconstraints%dumax2)
          case (1)
            ! Initialise the space constraints.
            call stlin_initSpaceConstraints (rconstraints,dtime,&
                rsolution%p_rspaceDiscr,rconstrSpace)
            
            ! Implement the constraints
            call smva_projectControlTstepVec (rtempVector%RvectorBlock(4),&
                rconstrSpace%p_rvectorumin%RvectorBlock(1),&
                rconstrSpace%p_rvectorumax%RvectorBlock(1))
            call smva_projectControlTstepVec (rtempVector%RvectorBlock(5),&
                rconstrSpace%p_rvectorumin%RvectorBlock(2),&
                rconstrSpace%p_rvectorumax%RvectorBlock(2))
            
            ! Done.
            call stlin_doneSpaceConstraints (rconstrSpace)
          end select
        end if
        
        !das hier gibt ein falsches Ergebnis1!
        call pperr_scalar (rtempVector%RvectorBlock(4),PPERR_L2ERROR,Derr(1))
        call pperr_scalar (rtempVector%RvectorBlock(5),PPERR_L2ERROR,Derr(2))
              
        ! We use the summed trapezoidal rule.
        if ((isubstep .eq. 1) .or. (isubstep .eq. rsolution%NEQtime)) then
          Derror(2) = Derror(2) + 0.05_DP*0.5_DP*(Derr(1)**2+Derr(2)**2) * dtstep
        else
          Derror(2) = Derror(2) + 0.5_DP*(Derr(1)**2+Derr(2)**2) * dtstep
        end if
        
      end select
      
    end do
    
    call collct_done(rcollection)
      
    ! Normalise...
    ! Derror(1) = Derror(1) / REAL(rsolution%NEQtime,DP)
    ! Derror(2) = Derror(2) / REAL(rsolution%NEQtime,DP)
    
    ! Calculate J(.)
    Derror(4) = 0.5_DP * Derror(1)  +  0.5_DP * dgamma * Derror(3)
    
    if (dalpha .ne. 0.0_DP) then
      ! Calculate:
      !    alpha/2 ||u||^2 = 1/2 ||P(lambda)||^2
      Derror(4) = Derror(4) + 0.5_DP * Derror(2)
      
      ! Calculate ||u|| = 1/alpha ||P(lambda)||
      Derror(2) = 1.0_DP/dalpha * sqrt(Derror(2))
    else
      Derror(2) = 0.0_DP
    end if
    
    ! And the rest
    Derror(3) = sqrt(Derror(3))
    Derror(1) = sqrt(Derror(1))
    
    ! Release temnp vector
    call lsysbl_releaseVector (rtempVector)
    
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine optcana_analyticalError (rglobalData,rconstraints,rsolution,rreference,&
      DerrorU,DerrorP,DerrorLambda,DerrorXi,boutput)

!<description>
  ! Computes the L2-error $||y-y0||_2$ of the given solution y to a reference
  ! solution y0 given as analytical function.
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  type(t_spacetimeVector), intent(IN) :: rsolution
  
  ! Analytic solution defining the reference function z.
  type(t_anSolution), intent(inout) :: rreference

  ! Constraints in the optimal control problem.
  type(t_optcconstraintsSpaceTime), intent(in) :: rconstraints
  
  ! Global settings for callback routines.
  type(t_globalData), intent(inout), target :: rglobalData
  
  ! Flag that determines if the error in each timestep is written to the terminal.
  logical, intent(in) :: boutput
!</input>

!<output>
  ! Returns information about the error. ||y-y0||_{L^2}.
  ! DerrorU(1) = on the interval [0,T]
  ! DerrorU(2) = on the interval [0,T)
  ! DerrorU(3) = on the interval (0,T]
  ! DerrorU(4) = on the interval (0,T)
  real(DP), dimension(:), intent(out) :: DerrorU

  ! Returns information about the error in the pressure. ||p-p0||_{L^2}.
  ! DerrorP(1) = on the interval [0,T]
  ! DerrorP(2) = on the interval [0,T)
  ! DerrorP(3) = on the interval (0,T]
  ! DerrorP(4) = on the interval (0,T)
  real(DP), dimension(:), intent(out) :: DerrorP

  ! Returns information about the error. ||lambda-lambda0||_{L^2}.
  ! DerrorLambda(1) = on the interval [0,T]
  ! DerrorLambda(2) = on the interval [0,T)
  ! DerrorLambda(3) = on the interval (0,T]
  ! DerrorLambda(4) = on the interval (0,T)
  real(DP), dimension(:), intent(out) :: DerrorLambda

  ! Returns information about the error in the pressure. ||xi-xi0||_{L^2}.
  ! DerrorXi(1) = on the interval [0,T]
  ! DerrorXi(2) = on the interval [0,T)
  ! DerrorXi(3) = on the interval (0,T]
  ! DerrorXi(4) = on the interval (0,T)
  real(DP), dimension(:), intent(out) :: DerrorXi
!</output>
  
!</subroutine>
    
    ! local variables
    integer :: isubstep,i
    real(DP) :: dtstep,dtimePrimal,dtimeDual
    real(DP) :: derrU, derrP, derrLambda, derrXi
    real(DP),dimension(6) :: Derr
    type(t_collection) :: rcollection
    type(t_collection), target :: ruserCollection
    type(t_vectorBlock) :: rtempVector
    real(dp), dimension(:), pointer :: p_Ddata
    
    ! Create a temp vector
    call lsysbl_createVectorBlock(rsolution%p_rspaceDiscr,rtempVector,.true.)
    call lsysbl_getbase_double (rtempVector,p_Ddata)

    ! Some basic initialisation.
    ! Initialise the collection for the assembly process with callback routines.
    ! This stores the simulation time in the collection and sets the
    ! current subvector z for the callback routines.
    call collct_init(rcollection)
    
    DerrorU(:) = 0.0_DP
    DerrorP(:) = 0.0_DP
    DerrorLambda(:) = 0.0_DP
    DerrorXi(:) = 0.0_DP

    do isubstep = 1,rsolution%NEQtime
    
      ! Current point in time.
      ! For the first iterate, take a look at time interval 1 to
      ! get the length and start point of the interval.
      ! For all further iterates, look at the time interval leading to
      ! that iterate.
      if (isubstep .gt. 1) then
        call tdiscr_getTimestep(rsolution%p_rtimeDiscr,isubstep-1,dtimePrimal,dtstep)
        dtimeDual = dtimePrimal - (1.0_DP-rsolution%p_rtimeDiscr%dtheta)*dtstep
      else
        call tdiscr_getTimestep(rsolution%p_rtimeDiscr,isubstep,dtstep=dtstep,dtimestart=dtimePrimal)
        dtimeDual = dtimePrimal
      end if

      ! Get the solution.
      ! Evaluate the space time function in rvector in the point
      ! in time dtime. Independent of the discretisation in time,
      ! this will give us a vector in space.
      ! Note 1: The dual solution is shifted by (1-dtheta)*dtstep!
      ! We therefore only have to evaluate once!
      ! Note 2: For Crank-Nicolson, the time discretisation discretises
      ! the primal pressure at the point of the dual velocity and
      ! the dual pressure at the time of the primal velocity.
      ! We compensate for this time shift during the error calculation.
      
      !CALL sptivec_getTimestepData (rsolution, 1+isubstep, rtempVector)
      call tmevl_evaluate(rsolution,dtimePrimal,rtempVector)

      ! Use our standard implementation to evaluate the error.
      if (rreference%ctype .ne. ANSOL_TP_ANALYTICAL) then
        call ansol_prepareEval (rreference,rcollection,"SOL",dtimePrimal)
      else
        ! Prepare the user-defined collection for assembly.
        call collct_init(ruserCollection)
      end if
  
      ! Save the function type to the collection, so the callback knows how
      ! to evaluate.
      rcollection%IquickAccess(2) = rreference%ctype
      
      ! The user-defined collection is the follower of rcollection.
      rcollection%p_rnextCollection => ruserCollection

      ! Perform error analysis to calculate and add 1/2||y-y0||^2_{L^2},...
      ! Primal velocity, dual pressure
      do i=1,2
        rcollection%IquickAccess(1) = i

        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          ! Analytically given data
          call user_initCollectForVecAssembly (rglobalData,&
              rreference%iid,i,dtimePrimal,ruserCollection)
        end if

        call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
            optcana_evalFunction,rcollection)

        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          call user_doneCollectForAssembly (rglobalData,ruserCollection)
        end if
      end do
      
      ! Primal pressure only in the 1st timestep.
      if (isubstep .eq. 0) then
        i=3
        rcollection%IquickAccess(1) = i
        
        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          ! Analytically given data
          call user_initCollectForVecAssembly (rglobalData,&
              rreference%iid,i,dtimePrimal,ruserCollection)
        end if
        
        call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
            optcana_evalFunction,rcollection)

        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          call user_doneCollectForAssembly (rglobalData,ruserCollection)
        end if
      end if
      
      if (isubstep .ne. rsolution%NEQtime) then
        i=6
        rcollection%IquickAccess(1) = i

        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          ! Analytically given data
          call user_initCollectForVecAssembly (rglobalData,&
              rreference%iid,i,dtimePrimal,ruserCollection)
        end if

        call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
            optcana_evalFunction,rcollection)

        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          call user_doneCollectForAssembly (rglobalData,ruserCollection)
        end if
      end if
          
      ! The same for the dual equation.
      ! Dual velocity, primal pressure.
      ! In rtempVector(4..6) is the dual solution at time dtimeDual,
      ! so we don't have to evaluate the function again!

      ! If we have a function in rreference, switch the time for ir.
      if (rreference%ctype .ne. ANSOL_TP_ANALYTICAL) then

        call ansol_doneEval (rcollection,"SOL")
        call ansol_prepareEval (rreference,rcollection,"SOL",dtimeDual)

      end if

      do i=4,5
        rcollection%IquickAccess(1) = i

        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          ! Analytically given data
          call user_initCollectForVecAssembly (rglobalData,&
              rreference%iid,i,dtimeDual,ruserCollection)
        end if

        call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
            optcana_evalFunction,rcollection)

        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          call user_doneCollectForAssembly (rglobalData,ruserCollection)
        end if
            
      end do
      
      if (isubstep .ne. 0) then
        i=3
        rcollection%IquickAccess(1) = i

        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          ! Analytically given data
          call user_initCollectForVecAssembly (rglobalData,&
              rreference%iid,i,dtimeDual,ruserCollection)
        end if

        call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
            optcana_evalFunction,rcollection)

        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          call user_doneCollectForAssembly (rglobalData,ruserCollection)
        end if
      end if
          
      ! Dual pressure only in the last timestep.
      if (isubstep .eq. rsolution%NEQtime) then
        i=6
        rcollection%IquickAccess(1) = i

        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          ! Analytically given data
          call user_initCollectForVecAssembly (rglobalData,&
              rreference%iid,i,dtimeDual,ruserCollection)
        end if

        call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
            optcana_evalFunction,rcollection)

        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          call user_doneCollectForAssembly (rglobalData,ruserCollection)
        end if
      end if

      ! Clean up the collection -- either ours or the user-defined one.
      if (rreference%ctype .ne. ANSOL_TP_ANALYTICAL) then
        call ansol_doneEval (rcollection,"SOL")
      else
        call collct_done (ruserCollection)
      end if

      ! Get the errors in that timestep.
      derrU = sqrt(Derr(1)**2 + Derr(2)**2)
      derrP = Derr(3)
      derrLambda = sqrt(Derr(4)**2 + Derr(5)**2)
      derrXi = Derr(6)

      ! We use the summed trapezoidal rule.
      ! Watch out with the start/end of the time interval when
      ! plugging the values into the error arrays.
      if (isubstep .eq. 1) then
        
        DerrorU(1)      = DerrorU(1)      + 0.5_DP*derrU**2
        DerrorP(1)      = DerrorP(1)      + 0.5_DP*derrP**2
        DerrorLambda(1) = DerrorLambda(1) + 0.5_DP*derrLambda**2
        DerrorXi(1)     = DerrorXi(1)     + 0.5_DP*derrXi**2

        DerrorU(2)      = DerrorU(2)      + 0.5_DP*derrU**2
        DerrorP(2)      = DerrorP(2)      + 0.5_DP*derrP**2
        DerrorLambda(2) = DerrorLambda(2) + 0.5_DP*derrLambda**2
        DerrorXi(2)     = DerrorXi(2)     + 0.5_DP*derrXi**2
      
      end if

      if (isubstep .eq. 2) then
        
        DerrorU(1)      = DerrorU(1)      + derrU**2
        DerrorP(1)      = DerrorP(1)      + derrP**2
        DerrorLambda(1) = DerrorLambda(1) + derrLambda**2
        DerrorXi(1)     = DerrorXi(1)     + derrXi**2

        DerrorU(2)      = DerrorU(2)      + derrU**2
        DerrorP(2)      = DerrorP(2)      + derrP**2
        DerrorLambda(2) = DerrorLambda(2) + derrLambda**2
        DerrorXi(2)     = DerrorXi(2)     + derrXi**2
      
        DerrorU(3)      = DerrorU(3)      + 0.5_DP*derrU**2
        DerrorP(3)      = DerrorP(3)      + 0.5_DP*derrP**2
        DerrorLambda(3) = DerrorLambda(3) + 0.5_DP*derrLambda**2
        DerrorXi(3)     = DerrorXi(3)     + 0.5_DP*derrXi**2

        DerrorU(4)      = DerrorU(4)      + 0.5_DP*derrU**2
        DerrorP(4)      = DerrorP(4)      + 0.5_DP*derrP**2
        DerrorLambda(4) = DerrorLambda(4) + 0.5_DP*derrLambda**2
        DerrorXi(4)     = DerrorXi(4)     + 0.5_DP*derrXi**2
      
      end if

      if ((isubstep .ge. 3) .and. (isubstep .le. rsolution%NEQtime-2)) then
        
        DerrorU(1)      = DerrorU(1)      + derrU**2
        DerrorP(1)      = DerrorP(1)      + derrP**2
        DerrorLambda(1) = DerrorLambda(1) + derrLambda**2
        DerrorXi(1)     = DerrorXi(1)     + derrXi**2

        DerrorU(2)      = DerrorU(2)      + derrU**2
        DerrorP(2)      = DerrorP(2)      + derrP**2
        DerrorLambda(2) = DerrorLambda(2) + derrLambda**2
        DerrorXi(2)     = DerrorXi(2)     + derrXi**2
      
        DerrorU(3)      = DerrorU(3)      + derrU**2
        DerrorP(3)      = DerrorP(3)      + derrP**2
        DerrorLambda(3) = DerrorLambda(3) + derrLambda**2
        DerrorXi(3)     = DerrorXi(3)     + derrXi**2

        DerrorU(4)      = Derroru(4)      + derrU**2
        DerrorP(4)      = DerrorP(4)      + derrP**2
        DerrorLambda(4) = DerrorLambda(4) + derrLambda**2
        DerrorXi(4)     = DerrorXi(4)     + derrXi**2
      
      end if
      
      if (isubstep .eq. rsolution%NEQtime-1) then
        
        DerrorU(1)      = DerrorU(1)      + derrU**2
        DerrorP(1)      = DerrorP(1)      + derrP**2
        DerrorLambda(1) = DerrorLambda(1) + derrLambda**2
        DerrorXi(1)     = DerrorXi(1)     + derrXi**2

        DerrorU(2)      = DerrorU(2)      + 0.5_DP*derrU**2
        DerrorP(2)      = DerrorP(2)      + 0.5_DP*derrP**2
        DerrorLambda(2) = DerrorLambda(2) + 0.5_DP*derrLambda**2
        DerrorXi(2)     = DerrorXi(2)     + 0.5_DP*derrXi**2
      
        DerrorU(3)      = DerrorU(3)      + derrU**2
        DerrorP(3)      = DerrorP(3)      + derrP**2
        DerrorLambda(3) = DerrorLambda(3) + derrLambda**2
        DerrorXi(3)     = DerrorXi(3)     + derrXi**2

        DerrorU(4)      = DerrorU(4)      + 0.5_DP*derrU**2
        DerrorP(4)      = DerrorP(4)      + 0.5_DP*derrP**2
        DerrorLambda(4) = DerrorLambda(4) + 0.5_DP*derrLambda**2
        DerrorXi(4)     = DerrorXi(4)     + 0.5_DP*derrXi**2
      
      end if

      if (isubstep .eq. rsolution%NEQtime) then
        
        DerrorU(1)      = DerrorU(1)      + 0.5_DP*derrU**2
        DerrorP(1)      = DerrorP(1)      + 0.5_DP*derrP**2
        DerrorLambda(1) = DerrorLambda(1) + 0.5_DP*derrLambda**2
        DerrorXi(1)     = DerrorXi(1)     + 0.5_DP*derrXi**2

        DerrorU(3)      = DerrorU(3)      + 0.5_DP*derrU**2
        DerrorP(3)      = DerrorP(3)      + 0.5_DP*derrP**2
        DerrorLambda(3) = DerrorLambda(3) + 0.5_DP*derrLambda**2
        DerrorXi(3)     = DerrorXi(3)     + 0.5_DP*derrXi**2

      end if
      
      if (boutput) then
        call output_line("error("//trim(sys_siL(isubstep,10))//") = "// &
            trim(sys_sdEL(Derr(1),10))//" / "//&
            trim(sys_sdEL(Derr(2),10))//" / "// &
            trim(sys_sdEL(Derr(3),10))//" / "// &
            trim(sys_sdEL(Derr(4),10))//" / "//&
            trim(sys_sdEL(Derr(5),10))//" / "// &
            trim(sys_sdEL(Derr(6),10)))
      end if
          
    end do

    ! Get the error return values.
    DerrorU = sqrt(DerrorU*dtstep)
    DerrorP = sqrt(DerrorP*dtstep)
    DerrorLambda = sqrt(DerrorLambda*dtstep)
    DerrorXi = sqrt(DerrorXi*dtstep)

    ! Release temnp vector and temp data
    call lsysbl_releaseVector (rtempVector)
    call collct_done(rcollection)

  end subroutine

end module

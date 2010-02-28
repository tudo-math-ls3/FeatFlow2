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
  
  use structuresoptc
    
  implicit none
  
  private
  
  public :: optcana_stationaryFunctional
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
  
  integer :: ierror
  
    ! Evaluate
    call ansol_evaluate (rcollection,"SOL",rcollection%IquickAccess(1),&
        Dvalues,npointsPerElement,nelements,Dpoints,rdomainIntSubset%p_Ielements,ierror)
        
    ! Check that this was ok. If yes, copy the data to the destination.
    if (ierror .ne. 0) then
      call output_line ('Error evaluating RHS function.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'optcana_evalFunction')
      call sys_halt()
    end if

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine optcana_stationaryFunctional (rglobalData,rsolution,rreference,dalpha,Derror)

!<description>
  ! This function calculates the value of the functional which is to be
  ! minimised in the stationary optimal control problem.
  ! The functional is defined as
  !   $$ J(y,u) = 1/2||y-z||_{L^2}^2  + \alpha/2||u||^2 $$
  ! over a spatial domain $\Omega$.
  !
  ! For this purpose, the routine evaluates the user-defined callback functions 
  ! user_ffunction_TargetX and user_ffunction_TargetY. The collection must be initialised
  ! for postprocessing before calling his routine!
!</description>
  
!<input>
  ! Global settings for callback routines.
  type(t_globalData), intent(inout), target :: rglobalData

  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(IN) :: rsolution
  
  ! Analytic solution defining the reference function z.
  type(t_anSolution), intent(inout) :: rreference

  ! Regularisation parameter $\alpha$.
  real(DP), intent(IN) :: dalpha
!</input>

!<output>
  ! Returns information about the error.
  ! Derror(1) = ||y-z||_{L^2}.
  ! Derror(2) = ||u||.
  ! Derror(3) = J(y,u) = 1/2||y-z||_{L^2}^2  + \alpha/2||u||^2.
  ! Norm of the error functional.
  real(DP), dimension(:), intent(OUT) :: Derror
!</output>
  
!</subroutine>
    
    ! local variables
    real(DP),dimension(2) :: Derr
    type(t_collection) :: rcollection
    
    ! Initialise the collection for the assembly process with callback routines.
    ! This stores the simulation time in the collection and sets the
    ! current subvector z for the callback routines.
    call collct_init(rcollection)

    ! Perform error analysis to calculate and add 1/2||y-z||_{L^2}.
    if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
      
      ! Perform error analysis to calculate and add 1/2||y-z||^2_{L^2}.
      call user_initCollectForAssembly (rglobalData,0.0_DP,rcollection)

      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
                        user_ffunction_TargetX,rcollection)

      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
                        user_ffunction_TargetY,rcollection)
          
      call user_doneCollectForAssembly (rglobalData,rcollection)
      
    else
    
      ! Use our standard implementation to evaluate the functional.
      call ansol_prepareEval (rreference,rcollection,"SOL",0.0_DP)

      ! X-velocity
      rcollection%IquickAccess(1) = 1
      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
          optcana_evalFunction,rcollection)

      ! Y-velocity
      rcollection%IquickAccess(1) = 2
      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
          optcana_evalFunction,rcollection)
          
      call ansol_doneEval (rcollection,"SOL")
      
    end if
                       
    Derror(1) = (0.5_DP*(Derr(1)**2+Derr(2)**2))
    Derror(2) = 0.0_DP
    
    ! Calculate \alpha/2||u||^2.
    if (dalpha .ne. 0.0_DP) then
      call pperr_scalar (rsolution%RvectorBlock(4),PPERR_L2ERROR,Derr(1))

      call pperr_scalar (rsolution%RvectorBlock(5),PPERR_L2ERROR,Derr(2))
                         
      Derror(2) = (0.5_DP*(Derr(1)**2+Derr(2)**2))
      
      ! Because of u=-lambda/alpha we have:
      !    alpha/2 ||u||^2 = alpha/2 ||lambda/alpha||^2 = 1/(2*alpha) ||lambda||^2
    end if
    
    Derror(3) = 0.5_DP * Derror(1)+(0.5_DP/dalpha) * Derror(2)
    Derror(1) = sqrt(Derror(1))
    Derror(2) = sqrt(Derror(2))
    
    ! Clean up    
    call collct_done(rcollection)
    
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine optcana_nonstatFunctional (rglobalData,rconstraints,rsolution,rreference,&
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

      ! Compute:
      ! Derror(1) = ||y-z||^2_{L^2}.
      if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
        
        ! Perform error analysis to calculate and add 1/2||y-z||^2_{L^2}.
        call user_initCollectForAssembly (rglobalData,dtime,rcollection)

        call pperr_scalar (rtempVector%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
            user_ffunction_TargetX,rcollection)

        call pperr_scalar (rtempVector%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
            user_ffunction_TargetY,rcollection)
            
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
        call smva_projectControlTimestep (rtempVector%RvectorBlock(4),&
            rconstraints%dumin1,rconstraints%dumax1)
        call smva_projectControlTimestep (rtempVector%RvectorBlock(5),&
            rconstraints%dumin2,rconstraints%dumax2)
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
      
    end do
    
    call collct_done(rcollection)
      
    ! Normalise...
    ! Derror(1) = Derror(1) / REAL(rsolution%NEQtime,DP)
    ! Derror(2) = Derror(2) / REAL(rsolution%NEQtime,DP)
    
    ! Calculate J(.)
    Derror(4) = 0.5_DP * Derror(1)  +  0.5_DP * dgamma * Derror(3)
    
    if (dalpha .ne. 0.0_DP) then
      ! Calculate:
      !    alpha/2 ||u||^2
      Derror(4) = Derror(4) + 0.5_DP*dalpha * Derror(2)
      
      ! Calculate ||u|| = 1/alpha ||lambda||
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
      derrorU,derrorP,derrorLambda,derrorXi,boutput)

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
  real(DP), intent(out) :: derrorU

  ! Returns information about the error in the pressure. ||p-p0||_{L^2}.
  real(DP), intent(out) :: derrorP

  ! Returns information about the error. ||lambda-lambda0||_{L^2}.
  real(DP), intent(out) :: derrorLambda

  ! Returns information about the error in the pressure. ||xi-xi0||_{L^2}.
  real(DP), intent(out) :: derrorXi
!</output>
  
!</subroutine>
    
    ! local variables
    integer :: isubstep,i
    real(DP) :: dtstep,dtimePrimal,dtimeDual
    real(DP),dimension(6) :: Derr
    type(t_collection) :: rcollection
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
    
    derrorU = 0.0_DP
    derrorP = 0.0_DP
    derrorLambda = 0.0_DP
    derrorXi = 0.0_DP

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
      call ansol_prepareEval (rreference,rcollection,"SOL",dtimePrimal)

      ! Perform error analysis to calculate and add 1/2||y-y0||^2_{L^2},...
      ! Primal velocity, dual pressure
      do i=1,2
        rcollection%IquickAccess(1) = i
        call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
            optcana_evalFunction,rcollection)
      end do
      
      ! Primal pressure only in the 1st timestep.
      if (isubstep .eq. 0) then
        i=3
        rcollection%IquickAccess(1) = i
        call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
            optcana_evalFunction,rcollection)
      end if
      
      i=6
      rcollection%IquickAccess(1) = i
      call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
          optcana_evalFunction,rcollection)
          
      call ansol_doneEval (rcollection,"SOL")
      
      ! Ignore the solution of the 0th step (initial solution).
      ! They don't contribute to the dual solution!
      
      if (isubstep .gt. 0) then
      
        ! The same for the dual equation.
        ! Dual velocity, primal pressure.
        ! In rtempVector(4..6) is the dual solution at time dtimeDual,
        ! so we don't have to evaluate the flow again!

        call ansol_prepareEval (rreference,rcollection,"SOL",dtimeDual)
        do i=4,5
          rcollection%IquickAccess(1) = i
          call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
              optcana_evalFunction,rcollection)
        end do
        
        i=3
        rcollection%IquickAccess(1) = i
        call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
            optcana_evalFunction,rcollection)
            
        ! Dual pressure only in the last timestep.
        if (isubstep .eq. rsolution%NEQtime) then
          i=6
          rcollection%IquickAccess(1) = i
          call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
              optcana_evalFunction,rcollection)
        end if

        call ansol_doneEval (rcollection,"SOL")
      else
        Derr(4:6) = 0.0_DP
      end if

      ! We use the summed trapezoidal rule.
      if (isubstep .eq. 1) then
        
        derrorU = derrorU + 0.5_DP*(Derr(1)**2 + Derr(2)**2) * dtstep
        derrorP = derrorP + Derr(3)**2 * dtstep

      else if ((isubstep .eq. 2) .or. (isubstep .eq. rsolution%NEQtime)) then

        ! 2nd substep is first solution in the dual.

        derrorU = derrorU + (Derr(1)**2 + Derr(2)**2) * dtstep
        derrorP = derrorP + Derr(3)**2 * dtstep
        
        derrorLambda = derrorLambda + 0.5_DP*(Derr(4)**2 + Derr(5)**2) * dtstep
        derrorXi = derrorXi + Derr(6)**2 * dtstep
      
      else
      
        derrorU = derrorU + (Derr(1)**2 + Derr(2)**2) * dtstep
        derrorP = derrorP + Derr(3)**2 * dtstep

        derrorLambda = derrorLambda + (Derr(4)**2 + Derr(5)**2) * dtstep
        derrorXi = derrorXi + Derr(6)**2 * dtstep
      
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
    derrorU = sqrt(derrorU)
    derrorP = sqrt(derrorP)
    derrorLambda = sqrt(derrorLambda)
    derrorXi = sqrt(derrorXi)

    ! Release temnp vector and temp data
    call lsysbl_releaseVector (rtempVector)
    call collct_done(rcollection)

  end subroutine

end module

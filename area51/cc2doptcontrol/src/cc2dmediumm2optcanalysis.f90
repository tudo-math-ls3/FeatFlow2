!##############################################################################
!# ****************************************************************************
!# <name> cc2dmmediumm2optcanalysis </name>
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
!# 1.) cc_optc_stationaryFunctional
!#     -> Calculates the value of the stationary functional which is to be
!#        minimised:
!#     $$ J(y,u) = 1/2||y-z||_{L^2} + \alpha/2||u||^2 $$
!#
!# 1.) cc_optc_nonstatFunctional
!#     -> Calculates the value of the stationary functional which is to be
!#        minimised:
!#     $$ J(y,u) = 1/2||y-z||^2_{L^2} + \alpha/2||u||^2_{L^2} + gamma/2||y(T)-z(T)||^2_{L^2}$$
!#
!# </purpose>
!##############################################################################

module cc2dmediumm2optcanalysis

  use fsystem
  use linearsystemblock
  use pprocerror
  use spacetimevectors
  use cc2dmedium_callback
  use cc2dmediumm2matvecassembly
    
  implicit none
  
contains

!******************************************************************************

!<subroutine>

  subroutine cc_optc_stationaryFunctional (rsolution,dalpha,Derror,rcollection)

!<description>
  ! This function calculates the value of the functional which is to be
  ! minimised in the stationary optimal control problem.
  ! The functional is defined as
  !   $$ J(y,u) = 1/2||y-z||_{L^2}^2  + \alpha/2||u||^2 $$
  ! over a spatial domain $\Omega$.
  !
  ! For this purpose, the routine evaluates the user-defined callback functions 
  ! ffunction_TargetX and ffunction_TargetY. The collection must be initialised
  ! for postprocessing before calling his routine!
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(IN) :: rsolution
  
  ! Regularisation parameter $\alpha$.
  real(DP), intent(IN) :: dalpha
  
  ! Collection structure of the main application. Is passed to the callback routines.
  type(t_collection), intent(INOUT) :: rcollection
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
    
    ! Perform error analysis to calculate and add 1/2||y-z||_{L^2}.
    call pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
                       ffunction_TargetX,rcollection)

    call pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
                       ffunction_TargetY,rcollection)
                       
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
    
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_optc_nonstatFunctional (rproblem,rsolution,rtempVector,&
      dalpha,dgamma,Derror)

!<description>
  ! This function calculates the value of the functional which is to be
  ! minimised in the stationary optimal control problem.
  ! The functional is defined as
  !   $$ J(y,u) = 1/2||y-z||^2_{L^2} + \alpha/2||u||^2_{L^2} + gamma/2||y(T)-z(T)||^2_{L^2}$$
  ! over a spatial domain $\Omega$.
  !
  ! For this purpose, the routine evaluates the user-defined callback functions 
  ! ffunction_TargetX and ffunction_TargetY. The collection must be initialised
  ! for postprocessing before calling his routine!
  !
  ! The function z is given implicitely in the problem structure rproblem
  ! and evaluated in ffunction_TargetX and ffunction_TargetY!
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  type(t_spacetimeVector), intent(IN) :: rsolution
  
  ! Regularisation parameter $\alpha$.
  real(DP), intent(IN) :: dalpha

  ! Regularisation parameter $\gamma$.
  real(DP), intent(IN) :: dgamma
!</input>

!<inputoutput>
  ! Problem structure defining z and the time discretisation.
  type(t_problem), intent(INOUT) :: rproblem

  ! A block temp vector with size and structure of the subvectors in rsolution.
  type(t_vectorBlock), intent(INOUT) :: rtempVector
!</inputoutput>

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
    real(DP) :: dtstep
    real(DP),dimension(2) :: Derr
    
    Derror(1:4) = 0.0_DP
    dtstep = rsolution%p_rtimeDiscretisation%dtstep

    do isubstep = 0,rsolution%NEQtime-1
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + &
          isubstep*dtstep
      rproblem%rtimedependence%itimestep = isubstep

      ! Get the solution.
      ! Evaluate the space time function in rvector in the point
      ! in time dtime. Independent of the discretisation in time,
      ! this will give us a vector in space.
      !CALL sptivec_getTimestepData (rsolution, 1+isubstep, rtempVector)
      call tmevl_evaluate(rsolution,rproblem%rtimedependence%dtime,rtempVector)

      ! Initialise the collection for the assembly process with callback routines.
      ! This stores the simulation time in the collection and sets the
      ! current subvector z for the callback routines.
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
      
      ! Compute:
      ! Derror(1) = ||y-z||^2_{L^2}.
      
      ! Perform error analysis to calculate and add 1/2||y-z||^2_{L^2}.
      call pperr_scalar (rtempVector%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
                        ffunction_TargetX,rproblem%rcollection)

      call pperr_scalar (rtempVector%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
                        ffunction_TargetY,rproblem%rcollection)

      ! We use the summed trapezoidal rule.
      if ((isubstep .eq. 0) .or. (isubstep .eq. rsolution%NEQtime-1)) then
        Derror(1) = Derror(1) + 0.5_DP*0.5_DP*(Derr(1)**2 + Derr(2)**2) * dtstep
      else
        Derror(1) = Derror(1) + 0.5_DP*(Derr(1)**2 + Derr(2)**2) * dtstep
      end if

      ! Compute:
      ! Derror(3) = ||y(T)-z(T)||^2
      if (isubstep .eq. rsolution%NEQtime-1) then
        Derror(3) = 0.5_DP*(Derr(1)**2+Derr(2)**2)
      end if
      
      ! Compute:
      ! Derror(2) = ||u|| = ||P[min/max](-1/alpha lambda)||^2_{L^2}.
      ! For that purpose, scale the lambda part and project it if necessary.
      call lsyssc_scaleVector (rtempVector%RvectorBlock(4),-1.0_DP/dalpha)
      call lsyssc_scaleVector (rtempVector%RvectorBlock(5),-1.0_DP/dalpha)
      
      call cc_projectControlTimestep (rtempVector%RvectorBlock(4),&
          rproblem%roptControl%dumin1,rproblem%roptControl%dumax1)
      call cc_projectControlTimestep (rtempVector%RvectorBlock(5),&
          rproblem%roptControl%dumin2,rproblem%roptControl%dumax2)
      
      call pperr_scalar (rtempVector%RvectorBlock(4),PPERR_L2ERROR,Derr(1))
      call pperr_scalar (rtempVector%RvectorBlock(5),PPERR_L2ERROR,Derr(2))
            
      ! We use the summed trapezoidal rule.             
      if ((isubstep .eq. 0) .or. (isubstep .eq. rsolution%NEQtime-1)) then
        Derror(2) = Derror(2) + 0.05_DP*0.5_DP*(Derr(1)**2+Derr(2)**2) * dtstep
      else
        Derror(2) = Derror(2) + 0.5_DP*(Derr(1)**2+Derr(2)**2) * dtstep
      end if
      
      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
    end do
    
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
    
  end subroutine

end module

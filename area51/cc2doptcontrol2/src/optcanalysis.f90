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
  use boundary
  use cubature
  use linearsystemscalar
  use linearsystemblock
  use pprocerror
  use collection
  use domainintegration
  use spatialdiscretisation
  use fparser
  use pprocerror
  use derivatives
  use feevaluation
  
  use analyticsolution
  
  use timediscretisation
  use timeevaluation
  use spacetimevectors
  use user_callback
  use spacematvecassembly
  use spacetimelinearsystem
  use newtonderivative
  use spatialbcdef
  use spacetimedirichletbcc
  
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

  ! ***************************************************************************

!<subroutine>

  subroutine fct_diffToTarget (cderivative, rdiscretisation, &
                                   nelements, npointsPerElement, Dpoints, &
                                   IdofsTest, rdomainIntSubset, &
                                   Dvalues, rcollection)
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
  !<description>
    ! Calculates the error "y-z" in the cubature points inside the domain
    ! of interest.
  !</description>
    
  !<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(in) :: cderivative
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
    
  !</input>
  
  !<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
  !</output>
    
  !</subroutine>
  
    integer :: icomponent,i,j
    type(t_vectorBlock), pointer :: p_rvectorBlock
    real(DP), dimension(:,:), allocatable :: DvaluesY
    real(DP) :: dx1,dy1,dx2,dy2,dx,dy
  
    ! Call the user-function to calculate z.
    ! Pass the user-defined collection.
    select case (rcollection%IquickAccess(2))
    case (1)
      ! User definec
      call user_fct_Target (cderivative, rdiscretisation, &
          nelements, npointsPerElement, Dpoints, &
          IdofsTest, rdomainIntSubset, &
          Dvalues, rcollection%p_rnextCollection)
    case (2)
      ! Standard implementation
      call optcana_evalFunction (cderivative, rdiscretisation, &
          nelements, npointsPerElement, Dpoints, &
          IdofsTest, rdomainIntSubset, &
          Dvalues, rcollection%p_rnextCollection)
    end select
        
    ! Fetch some parameters
    icomponent = rcollection%IquickAccess(1)
    p_rvectorBlock => rcollection%p_rvectorQuickAccess1
    dx1 = rcollection%DquickAccess(1)
    dy1 = rcollection%DquickAccess(2)
    dx2 = rcollection%DquickAccess(3)
    dy2 = rcollection%DquickAccess(4)
    
    ! Calculate y-z
    allocate (DvaluesY(ubound(Dvalues,1),ubound(Dvalues,2)))
    
    call fevl_evaluate_sim (p_rvectorBlock%RvectorBlock(icomponent), &
        rdomainIntSubset, DER_FUNC, DvaluesY)
    do i=1,ubound(Dvalues,2)
      do j=1,ubound(Dvalues,1)
      
        ! Check if the point is in the domain of integest.
        dx = Dpoints(1,j,i)
        dy = Dpoints(2,j,i)
        if ((dx .ge. dx1) .and. (dy .ge. dy1) .and. (dx .le. dx2) .and. (dy .le. dy2)) then
          Dvalues(j,i) = (DvaluesY(j,i) - Dvalues(j,i))
        else
          Dvalues(j,i) = 0.0_DP
        end if
        
      end do
    end do
    
    deallocate (DvaluesY)

  end subroutine

!******************************************************************************

    subroutine ffunction_dirichletbcU (cderivative, rdiscretisation, &
                                   DpointsRef, Dpoints, ibct, DpointPar,&
                                   Ielements, Dvalues, rcollection)

    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
  !<description>
    ! Calculates the norm of the control ||u||^2 = ||1/alpha (nu partial_n lambda - xi*n)||^2
    ! in cubature points. Used to evaluate ||u||_boundary
  !</description>
    
  !<input>
    integer, intent(in) :: cderivative
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    real(DP), dimension(:,:,:), intent(in) :: DpointsRef
    real(DP), dimension(:,:,:), intent(in) :: Dpoints
    integer, intent(in) :: ibct
    real(DP), dimension(:,:), intent(in) :: DpointPar
    integer, dimension(:), intent(in) :: Ielements
    type(t_collection), intent(inout), optional :: rcollection
  !</input>
  
  !<output>
    real(DP), dimension(:,:), intent(out) :: Dvalues
  !</output>
    
  !</subroutine>
  
      ! local variables
      real(DP) :: dbetaC,dnx,dny
      type(t_vectorBlock), pointer :: p_rvector
      integer :: iel,ipt
      real(DP), dimension(:,:,:), allocatable :: DvecValues
      
      ! Dummys
      integer, dimension(1,1) :: IdofsTest
      
      ! Get the data for the evaluation
      p_rvector => rcollection%p_rvectorQuickAccess1
      dbetaC = rcollection%DquickAccess(1)
      
      allocate (DvecValues(ubound(Dvalues,1),ubound(Dvalues,2),8))
      
      ! Calculate derivative of lambda and xi.
      call fevl_evaluate_sim (DER_FUNC, DvecValues(:,:,1), p_rvector%RvectorBlock(6), &
          Dpoints, Ielements)
      call fevl_evaluate_sim (DER_DERIV2D_X, DvecValues(:,:,2), p_rvector%RvectorBlock(4), &
          Dpoints, Ielements)
      call fevl_evaluate_sim (DER_DERIV2D_Y, DvecValues(:,:,3), p_rvector%RvectorBlock(4), &
          Dpoints, Ielements)
      call fevl_evaluate_sim (DER_DERIV2D_X, DvecValues(:,:,4), p_rvector%RvectorBlock(5), &
          Dpoints, Ielements)
      call fevl_evaluate_sim (DER_DERIV2D_Y, DvecValues(:,:,5), p_rvector%RvectorBlock(5), &
          Dpoints, Ielements)
          
      ! Calculate normal vectors
      call boundary_getNormalVec2D_sim(rdiscretisation%p_rboundary, ibct, DpointPar, &
          DvecValues(:,:,6), DvecValues(:,:,7), cparType=BDR_PAR_LENGTH)
          
      ! Calculate NU directly into the output array.
      ! NOTE: DOES NOT YET WORK WITH NONCONSTANT COEFFICIENTS!!!
      call smva_calcViscosity (DvecValues(:,:,8),1,&
          rdiscretisation,ubound(Dpoints,3),ubound(Dpoints,2),&
          Dpoints,Ielements,rcollection%p_rnextCollection)
      
      ! Calculate the values
      do iel = 1,size(Ielements)
        do ipt = 1,ubound(DpointPar,1)
            
          Dvalues(ipt,iel) = ( (DvecValues(ipt,iel,8)*DvecValues(ipt,iel,2)*DvecValues(ipt,iel,6) + &
                                DvecValues(ipt,iel,8)*DvecValues(ipt,iel,3)*DvecValues(ipt,iel,7) + &
                                DvecValues(ipt,iel,1)*DvecValues(ipt,iel,6)) ** 2 + &
                               (DvecValues(ipt,iel,8)*DvecValues(ipt,iel,4)*DvecValues(ipt,iel,6) + &
                                DvecValues(ipt,iel,8)*DvecValues(ipt,iel,5)*DvecValues(ipt,iel,7) + &
                                DvecValues(ipt,iel,1)*DvecValues(ipt,iel,7)) ** 2 ) / (dbetaC*dbetaC)
        
        end do
      end do
      
      deallocate(DvecValues)
  
    end subroutine 
    
!******************************************************************************

!<subroutine>

  subroutine optcana_nonstatFunctional (rglobalData,rphysics,rconstraints,roptcBDC,&
      rsolution,rreference,dalphaC,dbetaC,dgammaC,Derror,DobservationArea)

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
  
  ! The boundary conditions
  type(t_optcBDC), intent(in) :: roptcBDC

  ! Regularisation parameter $\alpha$.
  real(DP), intent(IN) :: dalphaC

  ! Regularisation parameter $\beta$.
  real(DP), intent(IN) :: dbetaC

  ! Regularisation parameter $\gamma$.
  real(DP), intent(IN) :: dgammaC

  ! Global settings for callback routines.
  type(t_globalData), intent(inout), target :: rglobalData
  
  ! OPTIONAL: Observation area where to observe the functional.
  real(DP), dimension(:), intent(in), optional :: DobservationArea
!</input>

!<output>
  ! Returns information about the error.
  ! Derror(1) = ||y-z||_{L^2}.
  ! Derror(2) = ||u||_{L^2}.
  ! Derror(3) = ||y(T)-z(T)||_{L^2}.
  ! Derror(4) = J(y,u).
  ! Derror(5) = ||u||_{L^2(Gamma_C)}.
  real(DP), dimension(:), intent(OUT) :: Derror
!</output>
  
!</subroutine>
    
    ! local variables
    integer :: isubstep,i
    real(DP) :: dtstep,dtime
    real(DP),dimension(2) :: Derr
    type(t_collection), target :: rcollection,rlocalcoll
    type(t_vectorBlock), target :: rtempVector, rzeroVector
    real(dp), dimension(:), pointer :: p_Dx
    type(t_optcconstraintsSpace) :: rconstrSpace
    type(t_sptiDirichletBCCBoundary) :: rdirichletBCC
    type(t_bdRegionEntry), pointer :: p_rbdRegionIterator
    
    ! Create a temp vector
    call lsysbl_createVectorBlock(rsolution%p_rspaceDiscr,rtempVector,.true.)
    call lsysbl_createVectorBlock(rsolution%p_rspaceDiscr,rzeroVector,.true.)
    
    call lsysbl_getbase_double (rtempVector,p_Dx)
    
    ! Initialise the collection for the assembly process with callback routines.
    ! This stores the simulation time in the collection and sets the
    ! current subvector z for the callback routines.
    call collct_init(rcollection)

    ! Assemble the dirichlet control boundary conditions
    call stdbcc_createDirichletBCCBd (rsolution%p_rspaceDiscr,rsolution%p_rtimeDiscr,&
        rdirichletBCC)
    call stdbcc_assembleDirichletBCCBd (roptcBDC,rdirichletBCC,rglobalData)
    
    Derror(1:5) = 0.0_DP

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

        ! The local callback function calculates y-z in the domain of
        ! integest. Pass the user-defined callback function
        ! and the necessary parameters.
        rlocalcoll%p_rnextCollection => rcollection
        rlocalcoll%p_rvectorQuickAccess1 => rtempVector
        
        ! Domain of interest
        rlocalcoll%DquickAccess(1) = -SYS_MAXREAL_DP
        rlocalcoll%DquickAccess(2) = -SYS_MAXREAL_DP
        rlocalcoll%DquickAccess(3) = SYS_MAXREAL_DP
        rlocalcoll%DquickAccess(4) = SYS_MAXREAL_DP
        
        if (present(DobservationArea)) then
          ! Put the observation area into the collection
          rlocalcoll%DquickAccess(1) = DobservationArea(1)
          rlocalcoll%DquickAccess(2) = DobservationArea(2)
          rlocalcoll%DquickAccess(3) = DobservationArea(3)
          rlocalcoll%DquickAccess(4) = DobservationArea(4)
        end if

        if (rreference%ctype .eq. ANSOL_TP_ANALYTICAL) then
          
          ! Mark: Reference function given by user_fct_Target.
          rlocalcoll%IquickAccess(2) = 1
          
          ! Perform error analysis to calculate and add 1/2||y-z||^2_{L^2}.
          call user_initCollectForVecAssembly (rglobalData,&
              rreference%iid,1,dtime,rcollection)
          
          ! X-velocity
          rlocalcoll%IquickAccess(1) = 1
          call pperr_scalar (PPERR_L2ERROR,Derr(1),rzeroVector%RvectorBlock(1),&
              fct_diffToTarget,rlocalcoll)

          call user_doneCollectForAssembly (rglobalData,rcollection)

          call user_initCollectForVecAssembly (rglobalData,&
              rreference%iid,2,dtime,rcollection)

          ! Y-velocity
          rlocalcoll%IquickAccess(1) = 2
          call pperr_scalar (PPERR_L2ERROR,Derr(2),rzeroVector%RvectorBlock(2),&
              fct_diffToTarget,rlocalcoll)
              
          call user_doneCollectForAssembly (rglobalData,rcollection)
          
        else
          
          ! Use our standard implementation to evaluate the functional.
          call ansol_prepareEval (rreference,rcollection,"SOL",dtime)

          ! Mark: Reference function given by optcana_evalFunction.
          rlocalcoll%IquickAccess(2) = 2

          ! X-velocity
          rlocalcoll%IquickAccess(1) = 1
          rcollection%IquickAccess(1) = 1
          call pperr_scalar (PPERR_L2ERROR,Derr(1),rzeroVector%RvectorBlock(1),&
              fct_diffToTarget,rlocalcoll)

          ! Y-velocity
          rlocalcoll%IquickAccess(1) = 2
          rcollection%IquickAccess(1) = 2
          call pperr_scalar (PPERR_L2ERROR,Derr(2),rzeroVector%RvectorBlock(2),&
              fct_diffToTarget,rlocalcoll)
              
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
        
        if (dalphaC .gt. 0.0_DP) then
          ! Compute:
          ! Derror(2) = ||u|| = ||P[min/max](-1/alpha lambda)||^2_{L^2}.
          
          ! At first, calculate P(-1/alpha lambda) -- or nothing,
          ! if distriobuted control is deactivated.
          if (rconstraints%ccontrolConstraints .ne. 0) then
            select case (rconstraints%ccontrolConstraintsType)
            case (0)
              call nwder_applyMinMaxProjByDof (1.0_DP,rtempVector%RvectorBlock(4),&
                  -1.0_DP/dalphaC,rtempVector%RvectorBlock(4),&
                  -1.0_DP/dalphaC,rtempVector%RvectorBlock(4),&
                  rconstraints%dumin1,rconstraints%dumax1)

              call nwder_applyMinMaxProjByDof (1.0_DP,rtempVector%RvectorBlock(5),&
                  -1.0_DP/dalphaC,rtempVector%RvectorBlock(5),&
                  -1.0_DP/dalphaC,rtempVector%RvectorBlock(5),&
                  rconstraints%dumin2,rconstraints%dumax2)

            case (1)
              ! Initialise the space constraints.
              call stlin_initSpaceConstraints (rconstraints,dtime,dtime,&
                  rsolution%p_rspaceDiscr,rconstrSpace)
              
              ! Implement the constraints
              call nwder_applyMinMaxProjByDof (1.0_DP,rtempVector%RvectorBlock(4),&
                  -1.0_DP/dalphaC,rtempVector%RvectorBlock(4),&
                  -1.0_DP/dalphaC,rtempVector%RvectorBlock(4),&
                  1.0_DP,1.0_DP,&
                  rconstrSpace%p_rvectorumin%RvectorBlock(1),&
                  rconstrSpace%p_rvectorumax%RvectorBlock(1))

              call nwder_applyMinMaxProjByDof (1.0_DP,rtempVector%RvectorBlock(5),&
                  -1.0_DP/dalphaC,rtempVector%RvectorBlock(5),&
                  -1.0_DP/dalphaC,rtempVector%RvectorBlock(5),&
                  1.0_DP,1.0_DP,&
                  rconstrSpace%p_rvectorumin%RvectorBlock(2),&
                  rconstrSpace%p_rvectorumax%RvectorBlock(2))

              ! Done.
              call stlin_doneSpaceConstraints (rconstrSpace)
            end select
          end if

          !das hier gibt ein falsches Ergebnis1!
          call pperr_scalar (PPERR_L2ERROR,Derr(1),rtempVector%RvectorBlock(4))
          call pperr_scalar (PPERR_L2ERROR,Derr(2),rtempVector%RvectorBlock(5))
                
          ! We use the summed trapezoidal rule.
          if ((isubstep .eq. 1) .or. (isubstep .eq. rsolution%NEQtime)) then
            Derror(2) = Derror(2) + 0.05_DP*0.5_DP*(Derr(1)**2+Derr(2)**2) * dtstep
          else
            Derror(2) = Derror(2) + 0.5_DP*(Derr(1)**2+Derr(2)**2) * dtstep
          end if
          
        end if

        if (dbetaC .gt. 0.0_DP) then
        
          ! Compute on the Dirichlet control boundary:
          ! Derror(2) = ||u||_GammaC = ||-1/alpha (nu*partial_n(u) + xi*n))||^2_{L^2(GammaC)}.
          Derr(1) = 0.0_DP
          
          ! Prepare the collection
          rlocalColl%p_rnextCollection => rcollection
          rlocalColl%DquickAccess(1) = dbetaC
          call smva_prepareViscoAssembly (rphysics,rcollection,rtempVector)
          
          ! Pass rtempVector%RvectorBlock(4) as dummy, filled with 0.
          call lsyssc_clearVector (rtempVector%RvectorBlock(4))
          
          ! Loop over the boundary regions where Dirichlet boundary control
          ! is applied.
          p_rbdRegionIterator => rdirichletBCC%p_RbdRegion(isubstep)%p_rprimalBdHead
          do i = 1,rdirichletBCC%p_RbdRegion(isubstep)%nregionsPrimal
          
            ! Calculate the integral ||u||^2
            call pperr_scalarBoundary2D (PPERR_L2ERROR, CUB_G4_1D, Derr(1),&
                p_rbdRegionIterator%rboundaryRegion, rtempVector%RvectorBlock(4),&
                ffunction_dirichletbcU, rlocalcoll)
                
            ! Next boundary region
            p_rbdRegionIterator => p_rbdRegionIterator%p_nextBdRegion
                
          end do
          
          ! We use the summed trapezoidal rule.
          ! Do not square Derr(1) here, since Derr(1) is already ||u||^2 !
          if ((isubstep .eq. 1) .or. (isubstep .eq. rsolution%NEQtime)) then
            Derror(5) = Derror(5) + 0.05_DP*0.5_DP * Derr(1) * dtstep
          else
            Derror(5) = Derror(5) + 0.5_DP * Derr(1)  * dtstep
          end if
          
        end if
        
      end select
      
    end do

    ! Release the boundary conditions again
    call stdbcc_releaseDirichletBCCBd(rdirichletBCC)

    ! Clean up the collection
    call collct_done(rcollection)
      
    ! Normalise...
    ! Derror(1) = Derror(1) / REAL(rsolution%NEQtime,DP)
    ! Derror(2) = Derror(2) / REAL(rsolution%NEQtime,DP)
    
    ! Calculate J(.)
    Derror(4) = 0.5_DP * Derror(1)  +  0.5_DP * dgammaC * Derror(3)
    
    if (dalphaC .gt. 0.0_DP) then
      ! Calculate:
      !    alpha/2 ||u||^2 = alpha/2 ||P(-1/alpha lambda)||^2
      Derror(4) = Derror(4) + 0.5_DP * dalphaC * Derror(2)
      
      ! Calculate ||u|| = sqrt(||P(-1/alpha lambda)||^2)
      Derror(2) = sqrt(Derror(2))
    else
      Derror(2) = 0.0_DP
    end if

    if (dbetaC .gt. 0.0_DP) then
      ! Calculate:
      !    alpha/2 ||u||^2 = beta/2 ||1/beta lambda)||^2
      Derror(4) = Derror(4) + 0.5_DP * dbetaC * Derror(5)
      
      ! Calculate ||u|| = sqrt(||1/beta lambda||^2)
      Derror(5) = sqrt(Derror(5))
    else
      Derror(5) = 0.0_DP
    end if
    
    ! And the rest
    Derror(3) = sqrt(Derror(3))
    Derror(1) = sqrt(Derror(1))
    
    ! Release temnp vector
    call lsysbl_releaseVector (rtempVector)
    call lsysbl_releaseVector (rzeroVector)
    
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

        call pperr_scalar (PPERR_L2ERROR,Derr(i),rtempVector%RvectorBlock(i),&
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
        
        call pperr_scalar (PPERR_L2ERROR,Derr(i),rtempVector%RvectorBlock(i),&
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

        call pperr_scalar (PPERR_L2ERROR,Derr(i),rtempVector%RvectorBlock(i),&
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

        call pperr_scalar (PPERR_L2ERROR,Derr(i),rtempVector%RvectorBlock(i),&
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

        call pperr_scalar (PPERR_L2ERROR,Derr(i),rtempVector%RvectorBlock(i),&
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

        call pperr_scalar (PPERR_L2ERROR,Derr(i),rtempVector%RvectorBlock(i),&
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

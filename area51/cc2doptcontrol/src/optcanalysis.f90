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

module optcanalysis

  use fsystem
  use linearsystemblock
  use pprocerror
  use spacetimevectors
  use user_callback
  use spacematvecassembly
    
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
    real(DP) :: dtstep,dtime
    real(DP),dimension(2) :: Derr
    real(dp), dimension(:), pointer :: p_Dx
    
    Derror(1:4) = 0.0_DP
    dtstep = rsolution%p_rtimeDiscretisation%dtstep
    
    call lsysbl_getbase_double (rtempVector,p_Dx)

    do isubstep = 0,rsolution%NEQtime-1
      ! Current point in time
      dtime = rproblem%rtimedependence%dtimeInit + isubstep*dtstep

      ! Get the solution.
      ! Evaluate the space time function in rvector in the point
      ! in time dtime. Independent of the discretisation in time,
      ! this will give us a vector in space.
      !CALL sptivec_getTimestepData (rsolution, 1+isubstep, rtempVector)
      call tmevl_evaluate(rsolution,dtime,rtempVector)

      ! Initialise the collection for the assembly process with callback routines.
      ! This stores the simulation time in the collection and sets the
      ! current subvector z for the callback routines.
      call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)
      
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
      
      if (rproblem%roptControl%ccontrolConstraints .ne. 0.0_DP) then
        call cc_projectControlTimestep (rtempVector%RvectorBlock(4),&
            rproblem%roptControl%dumin1,rproblem%roptControl%dumax1)
        call cc_projectControlTimestep (rtempVector%RvectorBlock(5),&
            rproblem%roptControl%dumin2,rproblem%roptControl%dumax2)
      end if
      
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

  ! ***************************************************************************
  
!<subroutine>

  subroutine ffunction_analyticalRef (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! Analytical reference solution given by an expression.
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(IN)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\#local DOF's in trial space,Number of elements)
  integer, dimension(:,:), intent(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(INOUT), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
!</output>
  
!</subroutine>

    real(DP) :: dtime
    integer :: i,j,icomponent
    
    real(DP), dimension(:), allocatable :: DvaluesAct

    type(t_fparser), pointer :: p_rparser
    real(dp), dimension(:,:), allocatable :: p_Dval
    
    ! Get the parser object with the RHS expressions from the collection
    p_rparser => collct_getvalue_pars (rcollection, 'SOLPARSER')
    
    ! Current time
    dtime = rcollection%DquickAccess(1)
    
    ! X/Y component
    icomponent = rcollection%IquickAccess(1)
    
    ! Prepare the array with the values for the function.
    ! X-coordinate, Y-coordinate, time.
    allocate(p_Dval(3,npointsPerElement*nelements))
    do i=1,nelements
      do j=1,npointsPerElement
        p_Dval (1,(i-1)*npointsPerElement+j) = Dpoints(1,j,i)
        p_Dval (2,(i-1)*npointsPerElement+j) = Dpoints(2,j,i)
        p_Dval (3,(i-1)*npointsPerElement+j) = dtime
      end do
    end do
    
    ! Evaluate the 1st expression for the X-rhs
    allocate(DvaluesAct(npointsPerElement*nelements))
    call fparser_evalFunction (p_rparser, icomponent, 2, p_Dval, DvaluesAct)

    ! Reshape the data, that's it.
    do i=0,nelements-1
      do j=1,npointsPerElement
        Dvalues(j,i+1) = DvaluesAct(i*npointsPerElement+j)
      end do
    end do
    
    deallocate(DvaluesAct)
    deallocate(p_Dval)

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_optc_analyticalError (rproblem,rsolution,rtempVector,&
      dalpha,dgamma,&
      ssolutionExpressionY1,ssolutionExpressionY2,ssolutionExpressionP,&
      ssolutionExpressionLAMBDA1,ssolutionExpressionLAMBDA2,ssolutionExpressionXI,&
      derrorU,derrorP,derrorLambda,derrorXi)

!<description>
  ! Computes the L2-error $||y-y0||_2$ of the given solution y to a reference
  ! solution y0 given as expression.
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  type(t_spacetimeVector), intent(IN) :: rsolution
  
  ! Regularisation parameter $\alpha$.
  real(DP), intent(IN) :: dalpha

  ! Regularisation parameter $\gamma$.
  real(DP), intent(IN) :: dgamma
  
  ! Expression specifying the analytical solution y0 / p0, primal/dual.
  character(LEN=SYS_STRLEN) :: ssolutionExpressionY1,ssolutionExpressionY2,ssolutionExpressionP
  character(LEN=SYS_STRLEN) :: ssolutionExpressionLAMBDA1,ssolutionExpressionLAMBDA2
  character(LEN=SYS_STRLEN) :: ssolutionExpressionXI
!</input>

!<inputoutput>
  ! Problem structure defining z and the time discretisation.
  type(t_problem), intent(INOUT) :: rproblem

  ! A block temp vector with size and structure of the subvectors in rsolution.
  type(t_vectorBlock), intent(INOUT) :: rtempVector
!</inputoutput>

!<output>
  ! Returns information about the error. ||y-y0||_{L^2}.
  real(DP), intent(OUT) :: derrorU

  ! Returns information about the error in the pressure. ||p-p0||_{L^2}.
  real(DP), intent(OUT) :: derrorP

  ! Returns information about the error. ||lambda-lambda0||_{L^2}.
  real(DP), intent(OUT) :: derrorLambda

  ! Returns information about the error in the pressure. ||xi-xi0||_{L^2}.
  real(DP), intent(OUT) :: derrorXi
!</output>
  
!</subroutine>
    
    ! local variables
    integer :: isubstep,i
    real(DP) :: dtstep,dtime
    real(DP),dimension(6) :: Derr
    character(LEN=10), dimension(3), parameter :: EXPR_VARIABLES = &
      (/'X    ','Y    ','TIME '/)
    type(t_collection) :: rcollection
    type(t_fparser) :: rsolParser
    real(dp), dimension(:), pointer :: p_Ddata
    
    call lsysbl_getbase_double (rtempVector,p_Ddata)
    
    derrorU = 0.0_DP
    derrorP = 0.0_DP
    derrorLambda = 0.0_DP
    derrorXi = 0.0_DP
    dtstep = rsolution%p_rtimeDiscretisation%dtstep
    
    ! Create a parser and to evaluate the expression.
    call fparser_create (rsolParser,6)
    
    ! Compile the two expressions+
    if (ssolutionExpressionY1 .ne. "") then
      call fparser_parseFunction (rsolParser,1, ssolutionExpressionY1, EXPR_VARIABLES)
    else
      call fparser_parseFunction (rsolParser,1, "0", EXPR_VARIABLES)
    end if
    if (ssolutionExpressionY2 .ne. "") then
      call fparser_parseFunction (rsolParser,2, ssolutionExpressionY2, EXPR_VARIABLES)
    else
      call fparser_parseFunction (rsolParser,2, "0", EXPR_VARIABLES)
    end if
    if (ssolutionExpressionP .ne. "") then
      call fparser_parseFunction (rsolParser,3, ssolutionExpressionP, EXPR_VARIABLES)
    else
      call fparser_parseFunction (rsolParser,3, "0", EXPR_VARIABLES)
    end if
    
    if (ssolutionExpressionLAMBDA1 .ne. "") then
      call fparser_parseFunction (rsolParser,4, ssolutionExpressionLAMBDA1, EXPR_VARIABLES)
    else
      call fparser_parseFunction (rsolParser,4, "0", EXPR_VARIABLES)
    end if
    if (ssolutionExpressionLAMBDA2 .ne. "") then
      call fparser_parseFunction (rsolParser,5, ssolutionExpressionLAMBDA2, EXPR_VARIABLES)
    else
      call fparser_parseFunction (rsolParser,5, "0", EXPR_VARIABLES)
    end if
    if (ssolutionExpressionXi .ne. "") then
      call fparser_parseFunction (rsolParser,6, ssolutionExpressionXi, EXPR_VARIABLES)
    else
      call fparser_parseFunction (rsolParser,6, "0", EXPR_VARIABLES)
    end if

    ! Put the parser to the problem collection to be usable in a callback
    ! routine.
    call collct_init(rcollection)
    call collct_setvalue_pars (rcollection, 'SOLPARSER', &
        rsolParser, .true.)
          
    do isubstep = 0,rsolution%NEQtime-1
      ! Current point in time
      dtime = rproblem%rtimedependence%dtimeInit + isubstep*dtstep

      rcollection%DquickAccess(1) = dtime

      ! Get the solution.
      ! Evaluate the space time function in rvector in the point
      ! in time dtime. Independent of the discretisation in time,
      ! this will give us a vector in space.
      !CALL sptivec_getTimestepData (rsolution, 1+isubstep, rtempVector)
      call tmevl_evaluate(rsolution,dtime,rtempVector)

      ! Compute:
      ! Derror(1) = ||y-z||^2_{L^2}.
      
      ! Perform error analysis to calculate and add 1/2||y-z||^2_{L^2},...
      do i=1,6
        rcollection%IquickAccess(1) = i
        call pperr_scalar (rtempVector%RvectorBlock(i),PPERR_L2ERROR,Derr(i),&
                          ffunction_analyticalRef,rcollection)
      end do

      ! We use the summed trapezoidal rule.
      if ((isubstep .eq. 0) .or. (isubstep .eq. rsolution%NEQtime-1)) then
        derrorU = derrorU + 0.5_DP*0.5_DP*(Derr(1)**2 + Derr(2)**2) * dtstep
        derrorP = derrorP + 0.5_DP*Derr(3)**2 * dtstep

        derrorLambda = derrorLambda + 0.5_DP*0.5_DP*(Derr(4)**2 + Derr(5)**2) * dtstep
        derrorXi = derrorXi + 0.5_DP*Derr(6)**2 * dtstep
      else
        derrorU = derrorU + 0.5_DP*(Derr(1)**2 + Derr(2)**2) * dtstep
        derrorP = derrorP + Derr(3)**2 * dtstep

        derrorLambda = derrorLambda + 0.5_DP*(Derr(4)**2 + Derr(5)**2) * dtstep
        derrorXi = derrorXi + Derr(6)**2 * dtstep
      end if
      call output_line("error("//trim(sys_siL(isubstep+1,10))//") = "// &
        trim(sys_sdEL(Derr(1),10))//" / "//&
        trim(sys_sdEL(Derr(2),10))//" / "// &
        trim(sys_sdEL(Derr(3),10))//" / "// &
        trim(sys_sdEL(Derr(4),10))//" / "//&
        trim(sys_sdEL(Derr(5),10))//" / "// &
        trim(sys_sdEL(Derr(6),10)))

    end do
    
    derrorU = sqrt(derrorU)
    derrorP = sqrt(derrorP)
    derrorLambda = sqrt(derrorLambda)
    derrorXi = sqrt(derrorXi)
    
    ! Clean up
    call collct_deletevalue (rcollection, 'SOLPARSER')
    call collct_done(rcollection)
    
    ! Release the parser
    call fparser_release(rsolParser)
    
  end subroutine

end module

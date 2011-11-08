!##############################################################################
!# ****************************************************************************
!# <name> timeanalysis </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to calculate the time error and time
!# derivative of a solution. The time error is the value
!# of a functional $||u_{n}-u*_{n}||$ in time for a solution
!# $u_{n} = u(t_{n})$, and the norm of the time derivatice of a solution.
!# The value of this (unknown) error is approximated
!# by a time error functional $J(.)$ which can be evaluated by different
!# means, depending on the actual problem to solve.
!#
!# The following routines can be found here:
!#
!# 1.) cc_timeErrorByPredictor
!#     -> Calculate an approximation to the time error using a predictor-like
!#        approach.
!#
!# 2.) cc_timeDerivative
!#     -> Calculate the (1st order) time derivative of a solution.
!# </purpose>
!##############################################################################

module timeanalysis

  use fsystem
  use linearalgebra
  use linearsystemblock
    
  implicit none
  
!<types>

!<typeblock>
  
  ! Time error analysis block. Contains values of different time
  ! error functionals.
  type t_timeError
  
    ! $||rel. error in U||_{L2}$
    real(DP) :: drelUL2

    ! $||rel. error in U||_{\max}$
    real(DP) :: drelUmax

    ! $||rel. error in P||_{L2}$
    real(DP) :: drelPL2

    ! $||rel. error in P||_{\max}$
    real(DP) :: drelPmax
  
  end type
  
!</typeblock>

!<typeblock>
  
  ! Block for norms of the time derivative
  type t_timeDerivatives
  
    ! $||rel. change in U||_{L2}$
    real(DP) :: drelUL2

    ! $||rel. change in U||_{\max}$
    real(DP) :: drelUmax

    ! $||rel. change in P||_{L2}$
    real(DP) :: drelPL2

    ! $||rel. change in P||_{\max}$
    real(DP) :: drelPmax
  
  end type
  
!</typeblock>

!</types>


!<constants>

!<constantblock description="Identifiers for norms of time-dependent solutions.">

  ! $ ||u||_{l2} $
  integer, parameter :: TNRM_L2U    = 1

  ! $ ||u||_{\max} $
  integer, parameter :: TNRM_LMAX   = 2

  ! $ ||p||_{\l2} $
  integer, parameter :: TNRM_P2U    = 3

  ! $ ||p||_{\max} $
  integer, parameter :: TNRM_PMAX   = 4

  ! $ \max ( ||u||_{l2} , ||p||_{l2} ) $
  integer, parameter :: TNRM_L2UP2U     = 5

  ! $ \max ( ||u||_{\max} , ||p||_{\max} ) $
  integer, parameter :: TNRM_L2MAXPMAX  = 6

  ! $ \max ( ||u||_{l2}   , ||p||_{l2} ,
  !          ||u||_{\max} , ||p||_{\max} ) $
  integer, parameter :: TNRM_MAX        = 7

  ! $ \min ( ||u||_{l2}   , ||p||_{l2} ,
  !          ||u||_{\max} , ||p||_{\max} ) $
  integer, parameter :: TNRM_MIN        = 8

!</constantblock>

!</constants>

contains

!******************************************************************************

!<function>

  real(DP) function cc_timeErrorByPredictor (ctimeErrorControl,&
                      rsolution,rpredSolution,rauxVector,&
                      rtimeError) &
           result(dtimeerror)

!<description>
  ! Calculates the value of the time error functional $J(u_{n})$ by using the
  ! explicit time error prediction method. For this purpose, the application
  ! must perform the following tasks:
  ! The time
  !  a) Perform a predictor step to compute a 'predictor solution' u1 at time
  !     $t_{n}$ (e.g. by 1st order method)
  !  b) Perform the real computation to compute a solution u2 at time
  !     $t_{n}$ (e.g. 2nd order method)
  !  c) Use an error functional to compute J(.) to compute a time error
  !     using u1 and u2
  ! Both solution vectors u1 and u2 must be specified to this routine.
  ! The routine uses these to calculate $J(u2)$. This value can then be
  ! used e.g. in the adaptive time stepping routines to calculate a new
  ! time step size.
!</description>

!<input>
  ! Type of norm to use in space/time; former IEPSAD.
  ! =TNRM_L2U      : Calculate dtimeerror=drelUL2; standard
  ! =TNRM_LMAX     : Calculate dtimeerror=drelUmax
  ! =TNRM_P2U      : Calculate dtimeerror=drelPL2
  ! =TNRM_PMAX     : Calculate dtimeerror=drelPmax
  ! =TNRM_L2UP2U   : Calculate dtimeerror=max(drelUL2,drelPL2)
  ! =TNRM_L2MAXPMAX: Calculate dtimeerror=max(drelUmax,drelPmax)
  ! =TNRM_MAX      : Calculate dtimeerror=max(drelUL2,drelPL2,drelUmax,drelPmax)
  ! =TNRM_MIN      : Calculate dtimeerror=min(drelUL2,drelPL2,drelUmax,drelPmax)
  integer                        :: ctimeErrorControl
  
  ! Solution vector u2.
  type(t_vectorBlock), intent(IN) :: rsolution
  
  ! Solution vector u1 of the predictor calculation.
  type(t_vectorBlock), intent(INOUT) :: rpredSolution
!</input>
  
!<inputoutput>
  ! Auxiliary vector; same structure and size as rsolution.
  type(t_vectorBlock), intent(INOUT) :: rauxVector
!</inputoutput>

!<output>
  ! OPTIONAL: Time error analysis block. Returns values of different time error functionals.
  type(t_timeError), intent(inout), target, optional :: rtimeError
!</output>

!<result>
  ! Value of the error functional.
!</result>

!</function>

    ! local variables
    real(DP) :: dtmp
    real(DP), dimension(3) :: Dnorms1,Dnorms2
    integer, dimension(3) :: Cnorms
    type(t_timeError),target :: rtimeErrorLocal
    type(t_timeError), pointer :: p_rtimeError

    ! Write the results of the time error analysis either to the local analysis
    ! block or to the one given as parameter.
    p_rtimeError => rtimeErrorLocal
    if (present(rtimeError)) p_rtimeError => rtimeError

    ! Calculate d:=u2-u1
    call lsysbl_vectorLinearComb (rSolution,rpredSolution,1.0_DP,-1.0_DP,rauxVector)
    
    ! Calculate the different norms of the error functionals for the standard
    ! (Navier-)Stokes equation -- all at once.
    !
    ! ||d||_l2
    Cnorms = LINALG_NORML2
    call lsysbl_vectorNormBlock (rauxVector,Cnorms,Dnorms1)
    call lsysbl_vectorNormBlock (rsolution,Cnorms,Dnorms2)

    ! Compatibility note: For full compatibility to the old CC2D version, one must
    ! test (dtmp .LE. 1.0_DP) everywhere here instead of (dtmp .EQ. 0.0_DP) !

    dtmp = sqrt( 0.5_DP * (Dnorms2(1)**2+Dnorms2(2)**2) )
    if (dtmp .eq. 0.0_DP) dtmp=1.0_DP
    p_rtimeError%drelUL2 = sqrt(0.5_DP * (Dnorms1(1)**2+Dnorms1(2)**2) ) / dtmp

    dtmp = Dnorms2(3)
    if (dtmp .eq. 0.0_DP) dtmp=1.0_DP
    p_rtimeError%drelPL2 = Dnorms1(3) / dtmp
    
    ! ||d||_max
    Cnorms = LINALG_NORMMAX
    call lsysbl_vectorNormBlock (rauxVector,Cnorms,Dnorms1)
    call lsysbl_vectorNormBlock (rsolution,Cnorms,Dnorms2)

    dtmp = max(Dnorms2(1),Dnorms2(2))
    if (dtmp .eq. 0.0_DP) dtmp=1.0_DP
    p_rtimeError%drelUmax = max(Dnorms1(1),Dnorms1(2)) / dtmp

    dtmp = Dnorms2(3)
    if (dtmp .eq. 0.0_DP) dtmp=1.0_DP
    p_rtimeError%drelPmax = Dnorms1(3) / dtmp
    
    ! Get the value of the error functional J(.)
    select case (ctimeErrorControl)
    case DEFAULT
      dtimeerror = p_rtimeError%drelUL2
    case (TNRM_LMAX)
      dtimeerror = p_rtimeError%drelUmax
    case (TNRM_P2U)
      dtimeerror = p_rtimeError%drelPL2
    case (TNRM_PMAX)
      dtimeerror = p_rtimeError%drelPmax
    case (TNRM_L2UP2U)
      dtimeerror = max(p_rtimeError%drelUL2,p_rtimeError%drelPL2)
    case (TNRM_L2MAXPMAX)
      dtimeerror = max(p_rtimeError%drelUmax,p_rtimeError%drelPmax)
    case (TNRM_MAX)
      dtimeerror = max(p_rtimeError%drelUL2,p_rtimeError%drelPL2,&
                       p_rtimeError%drelUmax,p_rtimeError%drelPmax)
    case (TNRM_MIN)
      dtimeerror = min(p_rtimeError%drelUL2,p_rtimeError%drelPL2,&
                       p_rtimeError%drelUmax,p_rtimeError%drelPmax)
    end select
    
  end function

!******************************************************************************

!<function>

  real(DP) function cc_timeDerivative (ctimeErrorControl,&
                      rsolutionNew,rsolutionOld,dtstep,rauxVector,rtimeDerivNorms) &
           result(dtimenorm)

!<description>
  ! Calculates the norm of the time derivative of a solution:
  ! $|| u{n+1} - u_n ||$ / dtstep.
  ! For this purpose, the caller must provide the solution before a time step
  ! and after a time step as well as the time step size.
!</description>

!<input>
  ! Type of norm to use in space/time; former IEPSAD.
  ! =TNRM_L2U      : Calculate dtimeerror=drelUL2; standard
  ! =TNRM_LMAX     : Calculate dtimeerror=drelUmax
  ! =TNRM_P2U      : Calculate dtimeerror=drelPL2
  ! =TNRM_PMAX     : Calculate dtimeerror=drelPmax
  ! =TNRM_L2UP2U   : Calculate dtimeerror=max(drelUL2,drelPL2)
  ! =TNRM_L2MAXPMAX: Calculate dtimeerror=max(drelUmax,drelPmax)
  ! =TNRM_MAX      : Calculate dtimeerror=max(drelUL2,drelPL2,drelUmax,drelPmax)
  ! =TNRM_MIN      : Calculate dtimeerror=min(drelUL2,drelPL2,drelUmax,drelPmax)
  integer                        :: ctimeErrorControl
  
  ! Solution vector $u_{n+1}$ at the end of the time step.
  type(t_vectorBlock), intent(IN) :: rsolutionNew
  
  ! Solution vector $u_n$ at the beginning of the time step.
  type(t_vectorBlock), intent(INOUT) :: rsolutionOld
  
  ! Length of time step
  real(DP), intent(IN) :: dtstep
!</input>
  
!<inputoutput>
  ! Auxiliary vector; same structure and size as rsolution.
  type(t_vectorBlock), intent(INOUT) :: rauxVector
!</inputoutput>

!<output>
  ! OPTIONAL: Time norm analysis block. Returns different norms of the
  ! time derivative.
  type(t_timeDerivatives), intent(INOUT), target, optional :: rtimeDerivNorms
!</output>

!<result>
  ! Norm of the time derivative: $|| u{n+1} - u_n ||$ / dtstep.
!</result>

!</function>

    ! local variables
    real(DP), dimension(3) :: Dnorms1
    integer, dimension(3) :: Cnorms
    type(t_timeDerivatives),target :: rtimeNormLocal
    type(t_timeDerivatives), pointer :: p_rtimeNorm
    integer :: nequ,neqp

    ! Write the results of the time error analysis either to the local analysis
    ! block or to the one given as parameter.
    p_rtimeNorm => rtimeNormLocal
    if (present(rtimeDerivNorms)) p_rtimeNorm => rtimeDerivNorms

    ! Calculate d:=u2-u1
    call lsysbl_vectorLinearComb (rSolutionNew,rsolutionOld,1.0_DP,-1.0_DP,rauxVector)
    
    ! Get the length of the subvectors
    nequ = rSolutionNew%RvectorBlock(1)%NEQ+rSolutionNew%RvectorBlock(2)%NEQ
    neqp = rSolutionNew%RvectorBlock(3)%NEQ
    
    ! Calculate the different norms of the error functionals for the standard
    ! (Navier-)Stokes equation -- all at once.
    !
    ! ||d||_l2 / dtstep
    Cnorms = LINALG_NORML2
    call lsysbl_vectorNormBlock (rauxVector,Cnorms,Dnorms1)
    
    p_rtimeNorm%drelUL2 = sqrt( 0.5_DP * (Dnorms1(1)**2+Dnorms1(2)**2) ) &
        / (sqrt(real(nequ,DP)) * dtstep)
    p_rtimeNorm%drelPL2 = Dnorms1(3) / (sqrt(real(neqp,DP)) * dtstep)

    ! ||d||_max / dtstep
    p_rtimeNorm%drelUmax = max(Dnorms1(1),Dnorms1(2)) / dtstep
    p_rtimeNorm%drelPmax = Dnorms1(3) / dtstep
    
    ! Return the value in the desired norm
    select case (ctimeErrorControl)
    case DEFAULT
      dtimenorm = p_rtimeNorm%drelUL2
    case (TNRM_LMAX)
      dtimenorm = p_rtimeNorm%drelUmax
    case (TNRM_P2U)
      dtimenorm = p_rtimeNorm%drelPL2
    case (TNRM_PMAX)
      dtimenorm = p_rtimeNorm%drelPmax
    case (TNRM_L2UP2U)
      dtimenorm = max(p_rtimeNorm%drelUL2,p_rtimeNorm%drelPL2)
    case (TNRM_L2MAXPMAX)
      dtimenorm = max(p_rtimeNorm%drelUmax,p_rtimeNorm%drelPmax)
    case (TNRM_MAX)
      dtimenorm = max(p_rtimeNorm%drelUL2,p_rtimeNorm%drelPL2,&
                      p_rtimeNorm%drelUmax,p_rtimeNorm%drelPmax)
    case (TNRM_MIN)
      dtimenorm = min(p_rtimeNorm%drelUL2,p_rtimeNorm%drelPL2,&
                      p_rtimeNorm%drelUmax,p_rtimeNorm%drelPmax)
    end select
    
  end function
  
end module

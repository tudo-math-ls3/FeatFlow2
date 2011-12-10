!##############################################################################
!# ****************************************************************************
!# <name> cc2dmmediumm2timeanalysis </name>
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

module cctimeanalysis

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

    ! $||rel. error in VS||_{L2}$
    real(DP) :: drelVSL2

    ! $||rel. error in US||_{\max}$
    real(DP) :: drelVSmax

    ! $||rel. error in V||_{L2}$
    real(DP) :: drelVFL2

    ! $||rel. error in U||_{\max}$
    real(DP) :: drelVFmax


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

    ! $||rel. change in VS||_{L2}$
    real(DP) :: drelVSL2

    ! $||rel. change in VS||_{\max}$
    real(DP) :: drelVSmax

    ! $||rel. change in VF||_{L2}$
    real(DP) :: drelVFL2

    ! $||rel. change in VF||_{\max}$
    real(DP) :: drelVFmax


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
  integer, parameter :: TNRM_LMAXU   = 2

  ! $ ||vS||_{l2} $
  integer, parameter :: TNRM_L2VS    = 3  

  ! $ ||vS||_{\max} $
  integer, parameter :: TNRM_LMAXVS   = 4

  ! $ ||vF||_{l2} $
  integer, parameter :: TNRM_L2VF    = 5  

  ! $ ||vF||_{\max} $
  integer, parameter :: TNRM_LMAXVF   = 6

  ! $ ||p||_{\l2} $
  integer, parameter :: TNRM_P2U    = 7

  ! $ ||p||_{\max} $
  integer, parameter :: TNRM_PMAX   = 8

  ! $ \max ( ||u||_{l2} , ||vS||_{l2} , ||vF||_{l2} , ||p||_{l2} ) $
  integer, parameter :: TNRM_L2UP2U     = 9

  ! $ \max ( ||u||_{\max} ,||vS||_{\max} , ||vF||_{\max} , ||p||_{\max} ) $
  integer, parameter :: TNRM_L2MAXPMAX  = 10

  ! $ \max ( ||u||_{l2},   ||vS||_{l2},    ||vF||_{l2},   ||p||_{l2} , 
  !          ||u||_{\max}, ||vS||_{\max},  ||vF||_{\max}, ||p||_{\max} ) $
  integer, parameter :: TNRM_MAX        = 11

  ! $ \min ( ||u||_{l2},   ||vS||_{l2},    ||vF||_{l2},   ||p||_{l2} , 
  !          ||u||_{\max}, ||vS||_{\max},  ||vF||_{\max}, ||p||_{\max} ) $
  integer, parameter :: TNRM_MIN        = 12

  ! $ \max ( ||u||_{l2} , ||vS||_{l2} , ||vF||_{l2} ) $
  integer, parameter :: TNRM_L2UV     = 13

  ! $ \max ( ||u||_{\max} ,||vS||_{\max},||vF||_{\max} ) $
  integer, parameter :: TNRM_L2MAXUMAXV  = 14

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
  ! =TNRM_L2VS      : Calculate dtimeerror=drelVSL2; standard
  ! =TNRM_L2VF      : Calculate dtimeerror=drelVFL2; standard
  ! =TNRM_LMAXU     : Calculate dtimeerror=drelUmax
  ! =TNRM_LMAXVS     : Calculate dtimeerror=drelVSmax
  ! =TNRM_LMAXVF     : Calculate dtimeerror=drelVFmax
  ! =TNRM_P2U      : Calculate dtimeerror=drelPL2
  ! =TNRM_PMAX     : Calculate dtimeerror=drelPmax
  ! =TNRM_L2UP2U   : Calculate dtimeerror=max(drelUL2,drelVSL2,drelVFL2,drelPL2)
  ! =TNRM_L2MAXPMAX: Calculate dtimeerror=max(drelUmax,drelVSmax,drelVFmax,drelPmax)

  ! =TNRM_MAX      : Calculate dtimeerror=max(drelUL2,drelVSL2,drelVFL2,drelPL2,
!					  drelUmax,drelVSmax,drelVFmax,drelPmax)

  ! =TNRM_MIN      : Calculate dtimeerror=min(drelUL2,drelVSL2,drelVFL2,drelPL2,
!  					  drelUmax,drelVSmax,drelVFmax,drelPmax)

  ! =TNRM_L2UV      : Calculate dtimeerror=max(drelUL2,drelVSL2,drelVFL2,)  
  ! =TNRM_L2MAXUMAXV : Calculate dtimeerror=max(drelUmax,drelVSmax,drelVFmax,) 
  integer                        :: ctimeErrorControl
  
  ! Solution vector y2.
  type(t_vectorBlock), intent(in) :: rsolution
  
  ! Solution vector y1 of the predictor calculation.
  type(t_vectorBlock), intent(inout) :: rpredSolution
!</input>
  
!<inputoutput>
  ! Auxiliary vector; same structure and size as rsolution.
  type(t_vectorBlock), intent(inout) :: rauxVector
!</inputoutput>

!<output>
  ! OPTIONAL: Time error analysis block. Returns values of different time error functionals.
  type(t_timeError), intent(out), target, optional :: rtimeError
!</output>

!<result>
  ! Value of the error functional.
!</result>

!</function>

    ! local variables
    real(DP) :: dtmp
    real(DP), dimension(7) :: Dnorms1,Dnorms2 ! 7: nblocks
    integer, dimension(7) :: Cnorms  ! 7: nblocks
    type(t_timeError),target :: rtimeErrorLocal

    ! Calculate d:=u2-u1
    call lsysbl_vectorLinearComb (rSolution,rpredSolution,1.0_DP,-1.0_DP,rauxVector)
    
    ! Calculate the different norms of the error functionals for the standard
    ! (Navier-)Stokes equation -- all at once.

!  ----------------------------- ||d||_l2  ----------------------------------------

    Cnorms = LINALG_NORML2
    call lsysbl_vectorNormBlock (rauxVector,Cnorms,Dnorms1)
    call lsysbl_vectorNormBlock (rsolution,Cnorms,Dnorms2)

    ! Compatibility note: For full compatibility to the old CC2D version, one must
    ! test (dtmp .LE. 1.0_DP) everywhere here instead of (dtmp .EQ. 0.0_DP) !
! #######################################
!       uS / norm, L2
! #######################################
    dtmp = sqrt( 0.5_DP * (Dnorms2(1)**2+Dnorms2(2)**2) )
    if (dtmp .eq. 0.0_DP) dtmp=1.0_DP
    rtimeErrorLocal%drelUL2 = sqrt(0.5_DP * (Dnorms1(1)**2+Dnorms1(2)**2) ) / dtmp
! #######################################
!       vS / norm, L2
! #######################################
    dtmp = sqrt( 0.5_DP * (Dnorms2(3)**2+Dnorms2(4)**2) )
    if (dtmp .eq. 0.0_DP) dtmp=1.0_DP
    rtimeErrorLocal%drelVSL2 = sqrt(0.5_DP * (Dnorms1(3)**2+Dnorms1(4)**2) ) / dtmp
! #######################################
!       vF / norm, L2
! #######################################
    dtmp = sqrt( 0.5_DP * (Dnorms2(5)**2+Dnorms2(6)**2) )
    if (dtmp .eq. 0.0_DP) dtmp=1.0_DP
    rtimeErrorLocal%drelVFL2 = sqrt(0.5_DP * (Dnorms1(5)**2+Dnorms1(6)**2) ) / dtmp
! #######################################
!       p / norm, L2
! #######################################
    dtmp = Dnorms2(7)
    if (dtmp .eq. 0.0_DP) dtmp=1.0_DP
    rtimeErrorLocal%drelPL2 = Dnorms1(7) / dtmp

!  ----------------------------- ||d||_max  ----------------------------------------
    
    Cnorms = LINALG_NORMMAX
    call lsysbl_vectorNormBlock (rauxVector,Cnorms,Dnorms1)
    call lsysbl_vectorNormBlock (rsolution,Cnorms,Dnorms2)

! *************************
!  uS / norm, max
! **************************
    dtmp = max(Dnorms2(1),Dnorms2(2))
    if (dtmp .eq. 0.0_DP) dtmp=1.0_DP
    rtimeErrorLocal%drelUmax = max(Dnorms1(1),Dnorms1(2)) / dtmp

! *************************
!  vS / norm, max
! **************************
    dtmp = max(Dnorms2(3),Dnorms2(4))
    if (dtmp .eq. 0.0_DP) dtmp=1.0_DP
    rtimeErrorLocal%drelVSmax = max(Dnorms1(3),Dnorms1(4)) / dtmp
! *************************
!  vF / norm, max
! **************************
    dtmp = max(Dnorms2(5),Dnorms2(6))
    if (dtmp .eq. 0.0_DP) dtmp=1.0_DP
    rtimeErrorLocal%drelVFmax = max(Dnorms1(5),Dnorms1(6)) / dtmp
! *************************
!  p / norm, max
! **************************
    dtmp = Dnorms2(7)
    if (dtmp .eq. 0.0_DP) dtmp=1.0_DP
    rtimeErrorLocal%drelPmax = Dnorms1(7) / dtmp
!----------------------------------------------------------------------------------    
    ! Get the value of the error functional J(.)
    select case (ctimeErrorControl)
    case DEFAULT
    case (TNRM_L2U)
      dtimeerror = rtimeErrorLocal%drelUL2
    case (TNRM_L2VS)
      dtimeerror = rtimeErrorLocal%drelVSL2
    case (TNRM_L2VF)
      dtimeerror = rtimeErrorLocal%drelVFL2
    case (TNRM_LMAXU)
      dtimeerror = rtimeErrorLocal%drelUmax
    case (TNRM_LMAXVS)
      dtimeerror = rtimeErrorLocal%drelVSmax
    case (TNRM_LMAXVF)
      dtimeerror = rtimeErrorLocal%drelVFmax
    case (TNRM_P2U)
      dtimeerror = rtimeErrorLocal%drelPL2
    case (TNRM_PMAX)
      dtimeerror = rtimeErrorLocal%drelPmax
    case (TNRM_L2UP2U)
      dtimeerror = max(rtimeErrorLocal%drelUL2,rtimeErrorLocal%drelVSL2,rtimeErrorLocal%drelVFL2,rtimeErrorLocal%drelPL2)
    case (TNRM_L2MAXPMAX)
      dtimeerror = max(rtimeErrorLocal%drelUmax,rtimeErrorLocal%drelVSmax,rtimeErrorLocal%drelVFmax,rtimeErrorLocal%drelPmax)
    case (TNRM_MAX)
      dtimeerror = max(rtimeErrorLocal%drelUL2,rtimeErrorLocal%drelVSL2,rtimeErrorLocal%drelVFL2,rtimeErrorLocal%drelPL2,&
                       rtimeErrorLocal%drelUmax,rtimeErrorLocal%drelVSmax,rtimeErrorLocal%drelVFmax,rtimeErrorLocal%drelPmax)
    case (TNRM_MIN)
      dtimeerror = min(rtimeErrorLocal%drelUL2,rtimeErrorLocal%drelVSL2,rtimeErrorLocal%drelVFL2,rtimeErrorLocal%drelPL2,&
                       rtimeErrorLocal%drelUmax,rtimeErrorLocal%drelVSmax,rtimeErrorLocal%drelVFmax,rtimeErrorLocal%drelPmax)
    case (TNRM_L2UV)
      dtimeerror = max(rtimeErrorLocal%drelUL2,rtimeErrorLocal%drelVSL2,rtimeErrorLocal%drelVFL2)
    case (TNRM_L2MAXUMAXV)
      dtimeerror = max(rtimeErrorLocal%drelUmax,rtimeErrorLocal%drelVSmax,rtimeErrorLocal%drelVFmax)
    end select
    
    if (present(rtimeError)) rtimeError = rtimeErrorLocal
    
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
  ! =TNRM_L2VS      : Calculate dtimeerror=drelVSL2;
  ! =TNRM_L2VF      : Calculate dtimeerror=drelVFL2;
  ! =TNRM_LMAXU     : Calculate dtimeerror=drelUmax
  ! =TNRM_LMAXVS     : Calculate dtimeerror=drelVSmax
  ! =TNRM_LMAXVF     : Calculate dtimeerror=drelVFmax
  ! =TNRM_P2U      : Calculate dtimeerror=drelPL2
  ! =TNRM_PMAX     : Calculate dtimeerror=drelPmax
  ! =TNRM_L2UP2U   : Calculate dtimeerror=max(drelUL2,drelVSL2,drelVFL2,drelPL2)
  ! =TNRM_L2MAXPMAX: Calculate dtimeerror=max(drelUmax,drelVSmax,drelVFmax,drelPmax)
  ! =TNRM_MAX      : Calculate dtimeerror=max(drelUL2,drelVSL2,drelVFL2,drelPL2,drelUmax,drelVSmax,drelVFmax,drelPmax)
  ! =TNRM_MIN      : Calculate dtimeerror=min(drelUL2,drelVSL2,drelVFL2,drelPL2,drelUmax,drelVSmax,drelVSmax,drelPmax) 
 
  integer                        :: ctimeErrorControl
  
  ! Solution vector $u_{n+1}$ at the end of the time step.
  type(t_vectorBlock), intent(in) :: rsolutionNew
  
  ! Solution vector $u_n$ at the beginning of the time step.
  type(t_vectorBlock), intent(inout) :: rsolutionOld
  
  ! Length of time step
  real(DP), intent(in) :: dtstep
!</input>
  
!<inputoutput>
  ! Auxiliary vector; same structure and size as rsolution.
  type(t_vectorBlock), intent(inout) :: rauxVector
!</inputoutput>

!<output>
  ! OPTIONAL: Time norm analysis block. Returns different norms of the
  ! time derivative.
  type(t_timeDerivatives), intent(inout), target, optional :: rtimeDerivNorms
!</output>

!<result>
  ! Norm of the time derivative: $|| y_{n+1} - y_n ||$ / dtstep.
!</result>

!</function>

    ! local variables
    real(DP), dimension(7) :: Dnorms1
    integer, dimension(7) :: Cnorms
    type(t_timeDerivatives),target :: rtimeNormLocal
    type(t_timeDerivatives), pointer :: p_rtimeNorm
    integer :: nequ,neqvS,neqvF,neqp

    ! Write the results of the time error analysis either to the local analysis
    ! block or to the one given as parameter.
    p_rtimeNorm => rtimeNormLocal
    if (present(rtimeDerivNorms)) p_rtimeNorm => rtimeDerivNorms

    ! Calculate d:=y2-y1
    call lsysbl_vectorLinearComb (rSolutionNew,rsolutionOld,1.0_DP,-1.0_DP,rauxVector)
    
    ! Get the length of the subvectors
    ! Get the length of the subvectors
    nequ  = rSolutionNew%RvectorBlock(1)%NEQ+rSolutionNew%RvectorBlock(2)%NEQ
    neqvS = rSolutionNew%RvectorBlock(3)%NEQ+rSolutionNew%RvectorBlock(4)%NEQ
    neqvF = rSolutionNew%RvectorBlock(5)%NEQ+rSolutionNew%RvectorBlock(6)%NEQ
    neqp  = rSolutionNew%RvectorBlock(7)%NEQ
    
    ! Calculate the different norms of the error functionals for the standard
    ! (Navier-)Stokes equation -- all at once.
    !
    ! ||d||_l2 / dtstep
    Cnorms = LINALG_NORML2
    call lsysbl_vectorNormBlock (rauxVector,Cnorms,Dnorms1)
    
    p_rtimeNorm%drelUL2 = sqrt( 0.5_DP * (Dnorms1(1)**2+Dnorms1(2)**2) ) &
        / (sqrt(real(nequ,DP)) * dtstep)
    p_rtimeNorm%drelVSL2 = sqrt( 0.5_DP * (Dnorms1(3)**2+Dnorms1(4)**2) ) &
        / (sqrt(real(neqvS,DP)) * dtstep)
    p_rtimeNorm%drelVFL2 = sqrt( 0.5_DP * (Dnorms1(5)**2+Dnorms1(6)**2) ) &
        / (sqrt(real(neqvF,DP)) * dtstep)
    p_rtimeNorm%drelPL2 = Dnorms1(7) / (sqrt(real(neqp,DP)) * dtstep)

    ! ||d||_max / dtstep
    Cnorms = LINALG_NORMMAX
    call lsysbl_vectorNormBlock (rauxVector,Cnorms,Dnorms1)

    p_rtimeNorm%drelUmax  = max(Dnorms1(1),Dnorms1(2)) / dtstep
    p_rtimeNorm%drelVSmax = max(Dnorms1(3),Dnorms1(4)) / dtstep
    p_rtimeNorm%drelVFmax = max(Dnorms1(5),Dnorms1(6)) / dtstep
    p_rtimeNorm%drelPmax  = Dnorms1(7) / dtstep
    
    ! Return the value in the desired norm
    select case (ctimeErrorControl)
    case DEFAULT
    case (TNRM_L2U)
      dtimenorm = p_rtimeNorm%drelUL2
    case (TNRM_L2VS)
      dtimenorm = p_rtimeNorm%drelVSL2
    case (TNRM_L2VF)
      dtimenorm = p_rtimeNorm%drelVFL2
    case (TNRM_LMAXU)
      dtimenorm = p_rtimeNorm%drelUmax
    case (TNRM_LMAXVS)
      dtimenorm = p_rtimeNorm%drelVSmax
    case (TNRM_LMAXVF)
      dtimenorm = p_rtimeNorm%drelVFmax
    case (TNRM_P2U)
      dtimenorm = p_rtimeNorm%drelPL2
    case (TNRM_PMAX)
      dtimenorm = p_rtimeNorm%drelPmax
    case (TNRM_L2UP2U)
      dtimenorm = max(p_rtimeNorm%drelUL2,p_rtimeNorm%drelVSL2,p_rtimeNorm%drelVFL2,p_rtimeNorm%drelPL2)
    case (TNRM_L2MAXPMAX)
      dtimenorm = max(p_rtimeNorm%drelUmax,p_rtimeNorm%drelVSmax,p_rtimeNorm%drelVFmax,p_rtimeNorm%drelPmax)
    case (TNRM_MAX)
      dtimenorm = max(p_rtimeNorm%drelUL2,p_rtimeNorm%drelVSL2,p_rtimeNorm%drelVFL2,p_rtimeNorm%drelPL2,&
                      p_rtimeNorm%drelUmax,p_rtimeNorm%drelVSmax,p_rtimeNorm%drelVFmax,p_rtimeNorm%drelPmax)
    case (TNRM_MIN)
      dtimenorm = min(p_rtimeNorm%drelUL2,p_rtimeNorm%drelVSL2,p_rtimeNorm%drelVFL2,p_rtimeNorm%drelPL2,&
                      p_rtimeNorm%drelUmax,p_rtimeNorm%drelVSmax,p_rtimeNorm%drelVFmax,p_rtimeNorm%drelPmax)
    case (TNRM_L2UV)
      dtimenorm = max(p_rtimeNorm%drelUL2,p_rtimeNorm%drelVSL2,p_rtimeNorm%drelVFL2)
    case (TNRM_L2MAXUMAXV)
      dtimenorm = max(p_rtimeNorm%drelUmax,p_rtimeNorm%drelVSmax,p_rtimeNorm%drelVFmax)
    end select    
  end function  
end module

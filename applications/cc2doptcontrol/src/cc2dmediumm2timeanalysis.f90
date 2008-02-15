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

MODULE cc2dmediumm2timeanalysis

  USE fsystem
  USE linearsystemblock
    
  IMPLICIT NONE
  
!<types>

!<typeblock>
  
  ! Time error analysis block. Contains values of different time 
  ! error functionals.
  TYPE t_timeError
  
    ! $||rel. error in U||_{L2}$
    REAL(DP) :: drelUL2

    ! $||rel. error in U||_{\max}$
    REAL(DP) :: drelUmax

    ! $||rel. error in P||_{L2}$
    REAL(DP) :: drelPL2

    ! $||rel. error in P||_{\max}$
    REAL(DP) :: drelPmax
  
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! Block for norms of the time derivative 
  TYPE t_timeDerivatives
  
    ! $||rel. change in U||_{L2}$
    REAL(DP) :: drelUL2

    ! $||rel. change in U||_{\max}$
    REAL(DP) :: drelUmax

    ! $||rel. change in P||_{L2}$
    REAL(DP) :: drelPL2

    ! $||rel. change in P||_{\max}$
    REAL(DP) :: drelPmax
  
  END TYPE
  
!</typeblock>

!</types>


!<constants>

!<constantblock description="Identifiers for norms of time-dependent solutions.">

  ! $ ||u||_{l2} $
  INTEGER, PARAMETER :: TNRM_L2U    = 1  

  ! $ ||u||_{\max} $
  INTEGER, PARAMETER :: TNRM_LMAX   = 2

  ! $ ||p||_{\l2} $
  INTEGER, PARAMETER :: TNRM_P2U    = 3

  ! $ ||p||_{\max} $
  INTEGER, PARAMETER :: TNRM_PMAX   = 4

  ! $ \max ( ||u||_{l2} , ||p||_{l2} ) $
  INTEGER, PARAMETER :: TNRM_L2UP2U     = 5

  ! $ \max ( ||u||_{\max} , ||p||_{\max} ) $
  INTEGER, PARAMETER :: TNRM_L2MAXPMAX  = 6

  ! $ \max ( ||u||_{l2}   , ||p||_{l2} , 
  !          ||u||_{\max} , ||p||_{\max} ) $
  INTEGER, PARAMETER :: TNRM_MAX        = 7

  ! $ \min ( ||u||_{l2}   , ||p||_{l2} , 
  !          ||u||_{\max} , ||p||_{\max} ) $
  INTEGER, PARAMETER :: TNRM_MIN        = 8

!</constantblock>

!</constants>

CONTAINS

!******************************************************************************

!<function>

  REAL(DP) FUNCTION cc_timeErrorByPredictor (ctimeErrorControl,&
                      rsolution,rpredSolution,rauxVector,&
                      rtimeError) &
           RESULT(dtimeerror)

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
  ! Type of normto use in space/time; former IEPSAD.
  ! =TNRM_L2U      : Calculate dtimeerror=drelUL2; standard
  ! =TNRM_LMAX     : Calculate dtimeerror=drelUmax
  ! =TNRM_P2U      : Calculate dtimeerror=drelPL2
  ! =TNRM_PMAX     : Calculate dtimeerror=drelPmax
  ! =TNRM_L2UP2U   : Calculate dtimeerror=max(drelUL2,drelPL2)
  ! =TNRM_L2MAXPMAX: Calculate dtimeerror=max(drelUmax,drelPmax)
  ! =TNRM_MAX      : Calculate dtimeerror=max(drelUL2,drelPL2,drelUmax,drelPmax)
  ! =TNRM_MIN      : Calculate dtimeerror=min(drelUL2,drelPL2,drelUmax,drelPmax)  
  INTEGER                        :: ctimeErrorControl
  
  ! Solution vector u2.
  TYPE(t_vectorBlock), INTENT(IN) :: rsolution
  
  ! Solution vector u1 of the predictor calculation.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rpredSolution
!</input>
  
!<inputoutput>
  ! Auxiliary vector; same structure and size as rsolution.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rauxVector
!</inputoutput>

!<output>
  ! OPTIONAL: Time error analysis block. Returns values of different time error functionals.
  TYPE(t_timeError), INTENT(OUT), TARGET, OPTIONAL :: rtimeError
!</output>

!<result>
  ! Value of the error functional.
!</result>

!</function>

    ! local variables
    REAL(DP) :: dtmp
    REAL(DP), DIMENSION(3) :: Dnorms1,Dnorms2
    INTEGER(PREC_VECIDX), DIMENSION(3) :: Cnorms
    TYPE(t_timeError),TARGET :: rtimeErrorLocal
    TYPE(t_timeError), POINTER :: p_rtimeError

    ! Write the results of the time error analysis either to the local analysis
    ! block or to the one given as parameter.
    p_rtimeError => rtimeErrorLocal
    IF (PRESENT(rtimeError)) p_rtimeError => rtimeError

    ! Calculate d:=u2-u1
    CALL lsysbl_vectorLinearComb (rSolution,rpredSolution,1.0_DP,-1.0_DP,rauxVector)
    
    ! Calculate the different norms of the error functionals for the standard
    ! (Navier-)Stokes equation -- all at once.
    !
    ! ||d||_l2
    Cnorms = LINALG_NORML2
    Dnorms1 = lsysbl_vectorNormBlock (rauxVector,Cnorms)
    Dnorms2 = lsysbl_vectorNormBlock (rsolution,Cnorms)

    ! Compatibility note: For full compatibility to the old CC2D version, one must
    ! test (dtmp .LE. 1.0_DP) everywhere here instead of (dtmp .EQ. 0.0_DP) !

    dtmp = SQRT( 0.5_DP * (Dnorms2(1)**2+Dnorms2(2)**2) )
    IF (dtmp .EQ. 0.0_DP) dtmp=1.0_DP
    p_rtimeError%drelUL2 = SQRT(0.5_DP * (Dnorms1(1)**2+Dnorms1(2)**2) ) / dtmp

    dtmp = Dnorms2(3)
    IF (dtmp .EQ. 0.0_DP) dtmp=1.0_DP
    p_rtimeError%drelPL2 = Dnorms1(3) / dtmp
    
    ! ||d||_max
    Cnorms = LINALG_NORMMAX
    Dnorms1 = lsysbl_vectorNormBlock (rauxVector,Cnorms)
    Dnorms2 = lsysbl_vectorNormBlock (rsolution,Cnorms)

    dtmp = MAX(Dnorms2(1),Dnorms2(2))
    IF (dtmp .EQ. 0.0_DP) dtmp=1.0_DP
    p_rtimeError%drelUmax = MAX(Dnorms1(1),Dnorms1(2)) / dtmp

    dtmp = Dnorms2(3)
    IF (dtmp .EQ. 0.0_DP) dtmp=1.0_DP
    p_rtimeError%drelPmax = Dnorms1(3) / dtmp
    
    ! Get the value of the error functional J(.)
    SELECT CASE (ctimeErrorControl)
    CASE DEFAULT
      dtimeerror = p_rtimeError%drelUL2
    CASE (TNRM_LMAX)
      dtimeerror = p_rtimeError%drelUmax
    CASE (TNRM_P2U)
      dtimeerror = p_rtimeError%drelPL2
    CASE (TNRM_PMAX)
      dtimeerror = p_rtimeError%drelPmax
    CASE (TNRM_L2UP2U)
      dtimeerror = MAX(p_rtimeError%drelUL2,p_rtimeError%drelPL2)
    CASE (TNRM_L2MAXPMAX)
      dtimeerror = MAX(p_rtimeError%drelUmax,p_rtimeError%drelPmax)
    CASE (TNRM_MAX)
      dtimeerror = MAX(p_rtimeError%drelUL2,p_rtimeError%drelPL2,&
                       p_rtimeError%drelUmax,p_rtimeError%drelPmax)
    CASE (TNRM_MIN)
      dtimeerror = MIN(p_rtimeError%drelUL2,p_rtimeError%drelPL2,&
                       p_rtimeError%drelUmax,p_rtimeError%drelPmax)
    END SELECT
    
  END FUNCTION

!******************************************************************************

!<function>

  REAL(DP) FUNCTION cc_timeDerivative (ctimeErrorControl,&
                      rsolutionNew,rsolutionOld,dtstep,rauxVector,rtimeDerivNorms) &
           RESULT(dtimenorm)

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
  INTEGER                        :: ctimeErrorControl
  
  ! Solution vector $u_{n+1}$ at the end of the time step.
  TYPE(t_vectorBlock), INTENT(IN) :: rsolutionNew
  
  ! Solution vector $u_n$ at the beginning of the time step.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rsolutionOld
  
  ! Length of time step
  REAL(DP), INTENT(IN) :: dtstep
!</input>
  
!<inputoutput>
  ! Auxiliary vector; same structure and size as rsolution.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rauxVector
!</inputoutput>

!<output>
  ! OPTIONAL: Time norm analysis block. Returns different norms of the 
  ! time derivative.
  TYPE(t_timeDerivatives), INTENT(INOUT), TARGET, OPTIONAL :: rtimeDerivNorms
!</output>

!<result>
  ! Norm of the time derivative: $|| u{n+1} - u_n ||$ / dtstep.
!</result>

!</function>

    ! local variables
    REAL(DP), DIMENSION(3) :: Dnorms1
    INTEGER(PREC_VECIDX), DIMENSION(3) :: Cnorms
    TYPE(t_timeDerivatives),TARGET :: rtimeNormLocal
    TYPE(t_timeDerivatives), POINTER :: p_rtimeNorm
    INTEGER(PREC_VECIDX) :: nequ,neqp

    ! Write the results of the time error analysis either to the local analysis
    ! block or to the one given as parameter.
    p_rtimeNorm => rtimeNormLocal
    IF (PRESENT(rtimeDerivNorms)) p_rtimeNorm => rtimeDerivNorms

    ! Calculate d:=u2-u1
    CALL lsysbl_vectorLinearComb (rSolutionNew,rsolutionOld,1.0_DP,-1.0_DP,rauxVector)
    
    ! Get the length of the subvectors
    nequ = rSolutionNew%RvectorBlock(1)%NEQ+rSolutionNew%RvectorBlock(2)%NEQ
    neqp = rSolutionNew%RvectorBlock(3)%NEQ
    
    ! Calculate the different norms of the error functionals for the standard
    ! (Navier-)Stokes equation -- all at once.
    !
    ! ||d||_l2 / dtstep
    Cnorms = LINALG_NORML2
    Dnorms1 = lsysbl_vectorNormBlock (rauxVector,Cnorms)
    
    p_rtimeNorm%drelUL2 = SQRT( 0.5_DP * (Dnorms1(1)**2+Dnorms1(2)**2) ) &
        / (SQRT(REAL(nequ,DP)) * dtstep)
    p_rtimeNorm%drelPL2 = Dnorms1(3) / (SQRT(REAL(neqp,DP)) * dtstep)

    ! ||d||_max / dtstep
    p_rtimeNorm%drelUmax = MAX(Dnorms1(1),Dnorms1(2)) / dtstep
    p_rtimeNorm%drelPmax = Dnorms1(3) / dtstep
    
    ! Return the value in the desired norm
    SELECT CASE (ctimeErrorControl)
    CASE DEFAULT
      dtimenorm = p_rtimeNorm%drelUL2
    CASE (TNRM_LMAX)
      dtimenorm = p_rtimeNorm%drelUmax
    CASE (TNRM_P2U)
      dtimenorm = p_rtimeNorm%drelPL2
    CASE (TNRM_PMAX)
      dtimenorm = p_rtimeNorm%drelPmax
    CASE (TNRM_L2UP2U)
      dtimenorm = MAX(p_rtimeNorm%drelUL2,p_rtimeNorm%drelPL2)
    CASE (TNRM_L2MAXPMAX)
      dtimenorm = MAX(p_rtimeNorm%drelUmax,p_rtimeNorm%drelPmax)
    CASE (TNRM_MAX)
      dtimenorm = MAX(p_rtimeNorm%drelUL2,p_rtimeNorm%drelPL2,&
                      p_rtimeNorm%drelUmax,p_rtimeNorm%drelPmax)
    CASE (TNRM_MIN)
      dtimenorm = MIN(p_rtimeNorm%drelUL2,p_rtimeNorm%drelPL2,&
                      p_rtimeNorm%drelUmax,p_rtimeNorm%drelPmax)
    END SELECT
    
  END FUNCTION
  
END MODULE

!##############################################################################
!# ****************************************************************************
!# <name> cc2dmmediumm2timeerror </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to calculate the time error, i.e. the value
!# of a functional $||u_{n}-u*_{n}||$ in time for a solution 
!# $u_{n} = u(t_{n})$. The value of this (unknown) error is approximated
!# by a time error functional $J(.)$ which can be evaluated by different
!# means, depending on the actual problem to solve.
!#
!# The following routines can be found here:
!#
!# 1.) c2d2_timeErrorByPredictor
!#     -> Calculate an approximation to the time error using a predictor-like
!#        approach.
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2timeerror

  USE fsystem
  USE linearsystemblock
    
  IMPLICIT NONE
  
!<types>

!<typeblock>
  
  ! Time error analysis block. Contains values of different time 
  ! error functionals.
  TYPE t_timeError
  
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

CONTAINS

!******************************************************************************

!<function>

  REAL(DP) FUNCTION c2d2_timeErrorByPredictor (ctimeErrorControl,&
                      rsolution,rpredSolution,rauxVector,&
                      rtimeError) &
           RESULT(dtimeerror)

!<description>
  ! Calculates the calue of the time error functional $J(u_{n})$ by using the
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
  ! Type of error control in space/time; former IEPSAD.
  ! =1: Calculate dtimeerror=drelUL2; standard
  ! =2: Calculate dtimeerror=drelUmax
  ! =3: Calculate dtimeerror=drelPL2
  ! =4: Calculate dtimeerror=drelPmax
  ! =5: Calculate dtimeerror=max(drelUL2,drelPL2)
  ! =6: Calculate dtimeerror=max(drelUmax,drelPmax)
  ! =7: Calculate dtimeerror=max(drelUL2,drelPL2,drelUmax,drelPmax)
  ! =8: Calculate dtimeerror=min(drelUL2,drelPL2,drelUmax,drelPmax)  
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
!/output>

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

    dtmp = SQRT(Dnorms2(1)**2+Dnorms2(2)**2)
    IF (dtmp .EQ. 0.0_DP) dtmp=1.0_DP
    p_rtimeError%drelUL2 = SQRT(Dnorms1(1)**2+Dnorms1(2)**2) / dtmp

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
    
    ! Get the vlue of the error functional J(.)
    SELECT CASE (ctimeErrorControl)
    CASE DEFAULT
      dtimeerror = p_rtimeError%drelUL2
    CASE (2)
      dtimeerror = p_rtimeError%drelUmax
    CASE (3)
      dtimeerror = p_rtimeError%drelPL2
    CASE (4)
      dtimeerror = p_rtimeError%drelPmax
    CASE (5)
      dtimeerror = MAX(p_rtimeError%drelUL2,p_rtimeError%drelPL2)
    CASE (6)
      dtimeerror = MAX(p_rtimeError%drelUmax,p_rtimeError%drelPmax)
    CASE (7)
      dtimeerror = MAX(p_rtimeError%drelUL2,p_rtimeError%drelPL2,&
                       p_rtimeError%drelUmax,p_rtimeError%drelPmax)
    CASE (8)
      dtimeerror = MIN(p_rtimeError%drelUL2,p_rtimeError%drelPL2,&
                       p_rtimeError%drelUmax,p_rtimeError%drelPmax)
    END SELECT
    
  END FUNCTION
  
END MODULE

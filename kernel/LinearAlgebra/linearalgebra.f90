!##############################################################################
!# ****************************************************************************
!# <name> LinearAlgebra </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic linear algebra routines and serves as a wrapper
!# to BLAS routines. The routines here replace the original LCPx, LSPx,...
!# routines of the old FEAT library.
!# </purpose>
!##############################################################################

MODULE linearalgebra

  USE fsystem

  IMPLICIT NONE
  
CONTAINS

!<subroutine>

  SUBROUTINE lalg_vectorCopy (Dx,Dy)
  
  !<description>
  
  ! Copies vector dx: Dy = Dx
  
  !</description>

  !<input>
  
  ! Source vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dx
  
  !</input>

  !<output>
  
  ! Destination vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dy
  
  !</output>
  
  CALL DCOPY(SIZE(Dx),Dx,1,Dy,1)
  
  !</description>
  
!</subroutine>

  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lalg_vectorScale (Dx,c)
  
  !<description>
  
  ! Scales a vector vector Dx: Dx = c * Dx
  
  !</description>

  !<inputoutput>
  
  ! Source and destination vector
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: Dx
  
  !</inputoutput>

  !<input>

  ! Multiplication factor
  REAL(DP), INTENT(IN) :: c

  !</input>
  
  CALL DSCAL(SIZE(Dx),c,Dx,1)
  
  !</description>
  
!</subroutine>

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorClear (Dx)
  
  !<description>
  
  ! Clears the vector dx: Dx = 0

  !</description>

  !<output>
  
  ! Destination vector to be cleared
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dx
  
  !</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!  DO i = 1,SIZE(Dx)
!    Dx(i) = 0.0_DP
!  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorLinearComb (Dx,Dy,cx,cy)
  
  !<description>
  
  ! Performs a linear combination: Dy = cx * Dx  +  cy * Dy
  
  !</description>

  !<input>
  
  ! First source vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dx
  
  ! Scaling factor for Dx
  REAL(DP), INTENT(IN)               :: cx

  ! Scaling factor for Dy
  REAL(DP), INTENT(IN)               :: cy
  
  !</input>

  !<inputoutput>
  
  ! Second source vector; also receives the result
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: Dy
  
  !</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP) :: c
  
  IF (cy.EQ.0.0_DP) THEN
    CALL DCOPY(SIZE(Dx),Dx,1,Dy,1)
    IF (cx.NE.1.0_DP) CALL DSCAL(SIZE(Dx),cx,Dy,1)
  ELSE IF (cy.EQ.1D0) THEN
    CALL DAXPY(SIZE(Dx),cx,Dx,1,Dy,1)
  ELSE
    c=cx/cy
    CALL DAXPY(SIZE(Dx),c,Dx,1,Dy,1)
    CALL DSCAL(SIZE(Dx),cy,Dy,1)
  ENDIF
  
  END SUBROUTINE
  
END MODULE

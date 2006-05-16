!##############################################################################
!# ****************************************************************************
!# <name> LinearAlgebra </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic linear algebra routines and serves as a wrapper
!# to BLAS routines. The routines here replace the original LCPx, LSPx,...
!# routines of the old FEAT library.
!#
!# The following routines can be found here:
!#
!# 1.) lalg_vectorCopy
!#     -> Copy a vector to another (former LCPx)
!#
!# 2.) lalg_vectorScale
!#     -> Scale a vector (former LSCx)
!#
!# 3.) lalg_vectorClear
!#     -> Clear a vector (former LCLx)
!#
!# 4.) lalg_vectorLinearComb
!#     -> Linear combination of two vectors (former LLCx)
!#
!# 5.) lalg_scalarProduct
!#     -> Calculate the scalar product of two vectors (former LSPx)
!#
!# 6.) lalg_norm
!#     -> Calculate a specific norm of a vector
!#
!# </purpose>
!##############################################################################

MODULE linearalgebra

  USE fsystem

  IMPLICIT NONE
  
!<constants>

!<constantblock description="Constants identifying vector norms">

  ! Sum of the entries
  INTEGER, PARAMETER :: LINALG_NORMSUM    = -1

  ! Euclidian vector norm
  INTEGER, PARAMETER :: LINALG_NORMEUCLID = 0

  ! $l_1$-norm
  INTEGER, PARAMETER :: LINALG_NORML1     = 1
  
  ! $l_2$-norm
  INTEGER, PARAMETER :: LINALG_NORML2     = 2
  
  ! max-norm
  INTEGER, PARAMETER :: LINALG_NORMMAX    = 3
  
!</constantblock>

!</constants>
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorCopyDble (Dx,Dy)
  
!<description>
  ! Copies a double precision vector dx: Dy = Dx
!</description>

!<input>
  
  ! Source vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dx
  
!</input>

!<output>
  
  ! Destination vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dy
  
!</output>
  
!</subroutine>

  CALL DCOPY(SIZE(Dx),Dx,1,Dy,1)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorCopySngl (Sx,Sy)
  
!<description>
  ! Copies a single precision vector: Sy = Sx
!</description>

!<input>
  
  ! Source vector
  REAL(SP), DIMENSION(:), INTENT(IN) :: Sx
  
!</input>

!<output>
  
  ! Destination vector
  REAL(SP), DIMENSION(:), INTENT(IN) :: Sy
  
!</output>
  
!</subroutine>

  CALL SCOPY(SIZE(Sx),Sx,1,Sy,1)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorCopyInt (Ix,Iy)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  INTEGER(I32), DIMENSION(:), INTENT(IN) :: Ix
  
!</input>

!<output>
  
  ! Destination vector
  INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Iy
  
!</output>
  
!</subroutine>

  INTEGER(I32) :: i
  
  ! Does not exist in BLAS!
  DO i=1,SIZE(Ix)
    Iy(i) = Ix(i)
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lalg_vectorScaleDble (Dx,dc)
  
!<description>
  ! Scales a double precision vector: Dx = dc * Dx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: Dx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  REAL(DP), INTENT(IN) :: dc

!</input>
  
!</subroutine>

  CALL DSCAL(SIZE(Dx),dc,Dx,1)
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lalg_vectorScaleSngl (Sx,sc)
  
!<description>
  ! Scales a single precision vector: Dx = sc * Dx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: Sx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  REAL(DP), INTENT(IN) :: sc

!</input>
  
!</subroutine>

  CALL SSCAL(SIZE(Sx),sc,Sx,1)
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lalg_vectorScaleInt (Ix,ic)
  
!<description>
  ! Scales a integer vector: Ix = c * Ix
!</description>

!<inputoutput>
  
  ! Source and destination vector
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: Ix
  
!</inputoutput>

!<input>

  ! Multiplication factor
  INTEGER(I32), INTENT(IN) :: ic

!</input>
  
!</subroutine>

  INTEGER(I32) :: i
  
  ! Does not exist in BLAS
  DO i=1,SIZE(Ix)
    Ix(i) = ic*Ix(i)
  END DO

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorClearDble (Dx)
  
!<description>
  ! Clears a double precision vector: Dx = 0
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
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

  DO i = 1,SIZE(Dx)
    Dx(i) = 0.0_DP
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorClearSngl (Sx)
  
!<description>
  ! Clears a single precision vector: Sx = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  REAL(SP), DIMENSION(:), INTENT(OUT) :: Sx
  
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

  DO i = 1,SIZE(Sx)
    Sx(i) = 0.0_SP
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorClearInt (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

  DO i = 1,SIZE(Ix)
    Ix(i) = 0
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorLinearCombDble (Dx,Dy,dcx,dcy)
  
!<description>
  ! Performs a linear combination: Dy = dcx * Dx  +  dcy * Dy
!</description>

!<input>
  
  ! First source vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dx
  
  ! Scaling factor for Dx
  REAL(DP), INTENT(IN)               :: dcx

  ! Scaling factor for Dy
  REAL(DP), INTENT(IN)               :: dcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: Dy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP) :: c
  
  IF (dcy .EQ. 0.0_DP) THEN
    CALL DCOPY(SIZE(Dx),Dx,1,Dy,1)
    IF (dcx .NE. 1.0_DP) CALL DSCAL(SIZE(Dx),dcx,Dy,1)
  ELSE IF (dcy .EQ. 1D0) THEN
    CALL DAXPY(SIZE(Dx),dcx,Dx,1,Dy,1)
  ELSE
    c=dcx/dcy
    CALL DAXPY(SIZE(Dx),c,Dx,1,Dy,1)
    CALL DSCAL(SIZE(Dx),dcy,Dy,1)
  ENDIF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorLinearCombSngl (Sx,Sy,scx,scy)
  
!<description>
  ! Performs a linear combination: Sy = scx * Sx  +  scy * Sy
!</description>

!<input>
  
  ! First source vector
  REAL(SP), DIMENSION(:), INTENT(IN) :: Sx
  
  ! Scaling factor for Dx
  REAL(SP), INTENT(IN)               :: scx

  ! Scaling factor for Dy
  REAL(SP), INTENT(IN)               :: scy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: Sy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(SP) :: c
  
  IF (scy .EQ. 0.0_SP) THEN
    CALL SCOPY(SIZE(Sx),Sx,1,Sy,1)
    IF (scx .NE. 1.0_SP) CALL SSCAL(SIZE(Sx),scx,Sy,1)
  ELSE IF (scy .EQ. 1D0) THEN
    CALL SAXPY(SIZE(Sx),scx,Sx,1,Sy,1)
  ELSE
    c=scx/scy
    CALL SAXPY(SIZE(Sx),c,Sx,1,Sy,1)
    CALL SSCAL(SIZE(Sx),scy,Sy,1)
  ENDIF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<function>

  REAL(DP) FUNCTION lalg_scalarProductDble (Dx,Dy) RESULT (res)
  
!<description>
  ! Calculates the scalar product of two double precision vectors: 
  ! res = <vector,vector>
!</description>

!<input>
  
  ! First source vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dx
  
  ! Second source vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dy
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  ! local variables
  INTEGER(I32) :: i
  
  res = Dx(1)*Dy(1)
  DO i=2,SIZE(Dx)
    res = res + Dx(i)*Dy(i)
  END DO
  
  END FUNCTION

  ! ***************************************************************************

!<function>

  REAL(SP) FUNCTION lalg_scalarProductSngl (Sx,Sy) RESULT (res)
  
!<description>
  ! Calculates the scalar product of two single precision vectors: 
  ! res = <vector,vector>
!</description>

!<input>
  
  ! First source vector
  REAL(SP), DIMENSION(:), INTENT(IN) :: Sx
  
  ! Second source vector
  REAL(SP), DIMENSION(:), INTENT(IN) :: Sy
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  ! local variables
  INTEGER(I32) :: i
  
  res = Sx(1)*Sy(1)
  DO i=2,SIZE(Sx)
    res = res + Sx(i)*Sy(i)
  END DO
  
  END FUNCTION

  ! ***************************************************************************

!<function>

  REAL(SP) FUNCTION lalg_scalarProductInt (Ix,Iy) RESULT (res)
  
!<description>
  ! Calculates the scalar product of two single precision vectors: 
  ! res = <vector,vector>
!</description>

!<input>
  
  ! First source vector
  INTEGER(I32), DIMENSION(:), INTENT(IN) :: Ix
  
  ! Second source vector
  INTEGER(I32), DIMENSION(:), INTENT(IN) :: Iy
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  ! local variables
  INTEGER(I32) :: i
  
  res = Ix(1)*Iy(1)
  DO i=2,SIZE(Ix)
    res = res + Ix(i)*Iy(i)
  END DO
  
  END FUNCTION

  ! ***************************************************************************

!<subroutine>

  REAL(DP) FUNCTION lalg_normDble (Dx,cnorm,iposMax) RESULT(resnorm)
  
!<description>
  ! Calculates the norm of a double precision vector. cnorm identifies the 
  ! type of norm to calculate.
!</description>

!<input>
  ! Vector to calculate the norm of.
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dx
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  INTEGER, INTENT(IN) :: cnorm
!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  INTEGER(I32), INTENT(OUT), OPTIONAL :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! < 0, if an error occurred (unknown norm).
!</result>
  
!</subroutine>

  INTEGER(I32) :: i,j

  ! Choose the norm to calculate
  SELECT CASE (cnorm)
  CASE (LINALG_NORMSUM)
    ! L1-norm: sum all entries
    resnorm = Dx(1)
    DO i=2,SIZE(Dx)
      resnorm = resnorm + Dx(i)
    END DO

  CASE (LINALG_NORMEUCLID)
    ! Euclidian norm = scalar product <vector,vector>
    resnorm = Dx(1)*Dx(1)
    DO i=2,SIZE(Dx)
      resnorm = resnorm + Dx(i)*Dx(i)
    END DO

  CASE (LINALG_NORML1)
    ! L1-norm: sum all entries, divide by sqrt(vector length).
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = Dx(1)
    DO i=2,SIZE(Dx)
      resnorm = resnorm + Dx(i)
    END DO
    resnorm = resnorm / REAL(SIZE(Dx),DP)

  CASE (LINALG_NORML2)
    ! l_2-norm - like euclidian norm, but divide by vector length.
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = Dx(1)*Dx(1)
    DO i=2,SIZE(Dx)
      resnorm = resnorm + Dx(i)*Dx(i)
    END DO
    resnorm = resnorm / DSQRT(REAL(SIZE(Dx),DP))
    
  CASE (LINALG_NORMMAX)
    ! MAX-norm. Find the absolute largest entry.
    resnorm = ABS(Dx(1))
    j=1
    DO i=2,SIZE(Dx)
      IF (ABS(Dx(i)) .GT. resnorm) THEN
        j = i
        resnorm = ABS(Dx(i))
      END IF
    END DO
    IF (PRESENT(iposMax)) iposMax = j
  CASE DEFAULT
    resnorm = -1.0_DP ! Unknown norm.
  END SELECT
    
  END FUNCTION

  ! ***************************************************************************

!<subroutine>

  REAL(SP) FUNCTION lalg_normSngl (Sx,cnorm,iposMax) RESULT(resnorm)
  
!<description>
  ! Calculates the norm of a single precision vector. cnorm identifies the 
  ! type of norm to calculate.
!</description>

!<input>
  ! Vector to calculate the norm of.
  REAL(SP), DIMENSION(:), INTENT(IN) :: Sx
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  INTEGER, INTENT(IN) :: cnorm
!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  INTEGER(I32), INTENT(OUT), OPTIONAL :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! < 0, if an error occurred (unknown norm).
!</result>
  
!</subroutine>

  INTEGER(I32) :: i,j

  ! Choose the norm to calculate
  SELECT CASE (cnorm)
  CASE (LINALG_NORMSUM)
    ! L1-norm: sum absolute value of all entries
    resnorm = Sx(1)
    DO i=2,SIZE(Sx)
      resnorm = resnorm + ABS(Sx(i))
    END DO

  CASE (LINALG_NORMEUCLID)
    ! Euclidian norm = scalar product <vector,vector>
    resnorm = Sx(1)*Sx(1)
    DO i=2,SIZE(Sx)
      resnorm = resnorm + Sx(i)*Sx(i)
    END DO

  CASE (LINALG_NORML1)
    ! L1-norm: sum sum absolute value of all entries, divide by sqrt(vector length).
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = Sx(1)
    DO i=2,SIZE(Sx)
      resnorm = resnorm + ABS(Sx(i))
    END DO
    resnorm = resnorm / REAL(SIZE(Sx),SP)

  CASE (LINALG_NORML2)
    ! l_2-norm - like euclidian norm, but divide by sqrt(vector length).
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = Sx(1)*Sx(1)
    DO i=2,SIZE(Sx)
      resnorm = resnorm + Sx(i)*Sx(i)
    END DO
    resnorm = resnorm / SQRT(REAL(SIZE(Sx),SP))
    
  CASE (LINALG_NORMMAX)
    ! MAX-norm. Find the absolute largest entry.
    resnorm = ABS(Sx(1))
    j=1
    DO i=2,SIZE(Sx)
      IF (ABS(Sx(i)) .GT. resnorm) THEN
        j = i
        resnorm = ABS(Sx(i))
      END IF
    END DO
    IF (PRESENT(iposMax)) iposMax = j
  CASE DEFAULT
    resnorm = -1.0_SP ! Unknown norm.
  END SELECT
    
  END FUNCTION

END MODULE

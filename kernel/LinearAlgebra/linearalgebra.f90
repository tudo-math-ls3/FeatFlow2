!##############################################################################
!# ****************************************************************************
!# <name> linearalgebra </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic linear algebra routines and serves as a wrapper
!# to BLAS routines. The routines here replace the original LCPx, LSPx,...
!# routines of the old FEAT library.
!#
!# The following routines can be found here:
!#
!# 1.) lalg_copyVectorXXX
!#     -> Copy a vector to another (former LCPx)
!#
!# 2.) lalg_scaleVectorXXX
!#     -> Scale a vector (former LSCx)
!#
!# 3.) lalg_clearVectorXXX
!#     -> Clear a vector (former LCLx)
!#
!# 4.) lalg_setVectorXXX
!#     -> Set a vector to a defined value
!#
!# 5.) lalg_vectorLinearCombXXX
!#     -> Linear combination of two vectors (former LLCx)
!#
!# 6.) lalg_scalarProductXXX
!#     -> Calculate the scalar product of two vectors (former LSPx)
!#
!# 7.) lalg_normXXX
!#     -> Calculate a specific norm of a vector
!#
!# 8.) lalg_vectorSortXXX
!#     -> Resort the entries of a vector according to a given permutation
!#
!# 9.) lalg_vectorAddScalar
!#     -> Adds a scalar to each entry of a vector
!# </purpose>
!##############################################################################

MODULE linearalgebra

  USE fsystem

  IMPLICIT NONE
  
!<constants>

!<constantblock description="Constants identifying vector norms">

  ! Sum of the absolute values of entries
  INTEGER, PARAMETER :: LINALG_NORMSUM    = -1

  ! Euclidian vector norm: (vector,vector)
  INTEGER, PARAMETER :: LINALG_NORMEUCLID = 0

  ! $l_1$-norm: 1/NEQ * sum(abs(entries))
  INTEGER, PARAMETER :: LINALG_NORML1     = 1
  
  ! $l_2$-norm: 1/sqrt(NEQ) * (vector,vector)
  INTEGER, PARAMETER :: LINALG_NORML2     = 2
  
  ! max-norm
  INTEGER, PARAMETER :: LINALG_NORMMAX    = 3
  
!</constantblock>

!</constants>
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_copyVectorDble (Dx,Dy)
  
!<description>
  ! Copies a double precision vector dx: Dy = Dx
!</description>

!<input>
  
  ! Source vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dx
  
!</input>

!<output>
  
  ! Destination vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dy
  
!</output>
  
!</subroutine>

  CALL DCOPY(SIZE(Dx),Dx,1,Dy,1)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_copyVectorSngl (Fx,Fy)
  
!<description>
  ! Copies a single precision vector: Fy = Fx
!</description>

!<input>
  
  ! Source vector
  REAL(SP), DIMENSION(:), INTENT(IN) :: Fx
  
!</input>

!<output>
  
  ! Destination vector
  REAL(SP), DIMENSION(:), INTENT(OUT) :: Fy
  
!</output>
  
!</subroutine>

  CALL SCOPY(SIZE(Fx),Fx,1,Fy,1)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_copyVectorSnglDbl (Fx,Dy)
  
!<description>
  ! Copies single precision vector to double precision vector: Dy = Fx
!</description>

!<input>
  
  ! Source vector
  REAL(SP), DIMENSION(:), INTENT(IN) :: Fx
  
!</input>

!<output>
  
  ! Destination vector
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dy
  
!</output>
  
!</subroutine>
  INTEGER(I32) :: i
  
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
  DO i=1,SIZE(Fx)
    Dy(i) = REAL(Fx(i),DP)
  END DO
!%OMP  end parallel do

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_copyVectorDblSngl (Dx,Fy)
  
!<description>
  ! Copies double precision vector to single precision vector: Fy = Dx
!</description>

!<input>
  
  ! Source vector
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dx
  
!</input>

!<output>
  
  ! Destination vector
  REAL(SP), DIMENSION(:), INTENT(OUT) :: Fy
  
!</output>
  
!</subroutine>
  INTEGER(I32) :: i
  
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
  DO i=1,SIZE(Dx)
    Fy(i) = REAL(Dx(i),SP)
  END DO
!%OMP  end parallel do

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_copyVectorInt (Ix,Iy)
  
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
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
  DO i=1,SIZE(Ix)
    Iy(i) = Ix(i)
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_copyVectorDble2D (Dx,Dy)
  
!<description>
  ! Copies a double precision vector dx: Dy = Dx
!</description>

!<input>
  
  ! Source vector
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dx
  
!</input>

!<output>
  
  ! Destination vector
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dy
  
!</output>
  
!</subroutine>

  CALL DCOPY(SIZE(Dx,1)*SIZE(Dx,2),Dx,1,Dy,1)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_copyVectorSngl2D (Fx,Fy)
  
!<description>
  ! Copies a single precision vector: Fy = Fx
!</description>

!<input>
  
  ! Source vector
  REAL(SP), DIMENSION(:,:), INTENT(IN) :: Fx
  
!</input>

!<output>
  
  ! Destination vector
  REAL(SP), DIMENSION(:,:), INTENT(OUT) :: Fy
  
!</output>
  
!</subroutine>

  CALL SCOPY(SIZE(Fx,1)*SIZE(Fx,2),Fx,1,Fy,1)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_copyVectorSnglDbl2D (Fx,Dy)
  
!<description>
  ! Copies single precision vector to double precision vector: Dy = Fx
!</description>

!<input>
  
  ! Source vector
  REAL(SP), DIMENSION(:,:), INTENT(IN) :: Fx
  
!</input>

!<output>
  
  ! Destination vector
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dy
  
!</output>
  
!</subroutine>
  INTEGER(I32) :: i,j
  
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  DO j=1,SIZE(Fx,2)
    DO i=1,SIZE(Fx,1)
      Dy(i,j) = REAL(Fx(i,j),DP)
    END DO
  END DO
!%OMP  end parallel do

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_copyVectorDblSngl2D (Dx,Fy)
  
!<description>
  ! Copies double precision vector to single precision vector: Fy = Dx
!</description>

!<input>
  
  ! Source vector
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dx
  
!</input>

!<output>
  
  ! Destination vector
  REAL(SP), DIMENSION(:,:), INTENT(OUT) :: Fy
  
!</output>
  
!</subroutine>
  INTEGER(I32) :: i,j
  
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  DO j=1,SIZE(Dx,2)
    DO i=1,SIZE(Dx,1)
      Fy(i,j) = REAL(Dx(i,j),SP)
    END DO
  END DO
!%OMP  end parallel do

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_copyVectorInt2D (Ix,Iy)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: Ix
  
!</input>

!<output>
  
  ! Destination vector
  INTEGER(I32), DIMENSION(:,:), INTENT(OUT) :: Iy
  
!</output>
  
!</subroutine>

  INTEGER(I32) :: i,j
  
  ! Does not exist in BLAS!
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  DO j=1,SIZE(Ix,2)
    DO i=1,SIZE(Ix,1)
      Iy(i,j) = Ix(i,j)
    END DO
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lalg_scaleVectorDble (Dx,dc)
  
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

  SUBROUTINE lalg_scaleVectorSngl (Fx,sc)
  
!<description>
  ! Scales a single precision vector: Dx = sc * Dx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: Fx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  REAL(SP), INTENT(IN) :: sc

!</input>
  
!</subroutine>

  CALL SSCAL(SIZE(Fx),sc,Fx,1)
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lalg_scaleVectorInt (Ix,ic)
  
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
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
  DO i=1,SIZE(Ix)
    Ix(i) = ic*Ix(i)
  END DO
!%OMP  end parallel do

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lalg_scaleVectorDble2D (Dx,dc)
  
!<description>
  ! Scales a double precision vector: Dx = dc * Dx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: Dx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  REAL(DP), INTENT(IN) :: dc

!</input>
  
!</subroutine>

  CALL DSCAL(SIZE(Dx,1)*SIZE(Dx,2),dc,Dx,1)
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lalg_scaleVectorSngl2D (Fx,sc)
  
!<description>
  ! Scales a single precision vector: Dx = sc * Dx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: Fx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  REAL(SP), INTENT(IN) :: sc

!</input>
  
!</subroutine>

  CALL SSCAL(SIZE(Fx,1)*SIZE(Fx,2),sc,Fx,1)
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lalg_scaleVectorInt2D (Ix,ic)
  
!<description>
  ! Scales a integer vector: Ix = c * Ix
!</description>

!<inputoutput>
  
  ! Source and destination vector
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: Ix
  
!</inputoutput>

!<input>

  ! Multiplication factor
  INTEGER(I32), INTENT(IN) :: ic

!</input>
  
!</subroutine>

  INTEGER(I32) :: i,j
  
  ! Does not exist in BLAS
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  DO j=1,SIZE(Ix,2)
    DO i=1,SIZE(Ix,1)
      Ix(i,j) = ic*Ix(i,j)
    END DO
  END DO
!%OMP  end parallel do

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_clearVectorDble (Dx)
  
!<description>
  ! Clears a double precision vector: Dx = 0
!</description>

!<output>
  ! Destination vector to be cleared
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  INTEGER(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
  DO i = 1,SIZE(Dx)
    Dx(i) = 0.0_DP
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_clearVectorSngl (Fx)
  
!<description>
  ! Clears a single precision vector: Fx = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  REAL(SP), DIMENSION(:), INTENT(OUT) :: Fx
  
!</output>
  
!</subroutine>

  ! local variables
  INTEGER(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
  DO i = 1,SIZE(Fx)
    Fx(i) = 0.0_SP
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_clearVectorInt (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  INTEGER(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
  DO i = 1,SIZE(Ix)
    Ix(i) = 0
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_clearVectorDble2D (Dx)
  
!<description>
  ! Clears a double precision vector: Dx = 0
!</description>

!<output>
  ! Destination vector to be cleared
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  INTEGER(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  DO j = 1,SIZE(Dx,2)
    DO i = 1,SIZE(Dx,1)
      Dx(i,j) = 0.0_DP
    END DO
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_clearVectorSngl2D (Fx)
  
!<description>
  ! Clears a single precision vector: Fx = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  REAL(SP), DIMENSION(:,:), INTENT(OUT) :: Fx
  
!</output>
  
!</subroutine>

  ! local variables
  INTEGER(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  DO j = 1,SIZE(Fx,2)
    DO i = 1,SIZE(Fx,1)
      Fx(i,j) = 0.0_SP
    END DO
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_clearVectorInt2D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  INTEGER(I32), DIMENSION(:,:), INTENT(OUT) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  INTEGER(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  DO j = 1,SIZE(Ix,2)
    DO i = 1,SIZE(Ix,1)
      Ix(i,j) = 0
    END DO
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_setVectorDble (Dx,dvalue)
  
!<description>
  ! Sets the vector data to a defined value: Dx = dvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  REAL(DP), INTENT(IN) :: dvalue
!</input>

!<output>
  ! Destination vector to be set
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  INTEGER(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = dvalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
  DO i = 1,SIZE(Dx)
    Dx(i) = dvalue
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_setVectorSngl (Fx,fvalue)
  
!<description>
  ! Sets the vector data to a defined value: Fx = fvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  REAL(SP), INTENT(IN) :: fvalue
!</input>

!<output>
  ! Destination vector to be set
  REAL(SP), DIMENSION(:), INTENT(OUT) :: Fx
!</output>

!</subroutine>

  ! local variables
  INTEGER(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Fx = fvalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
  DO i = 1,SIZE(Fx)
    Fx(i) = fvalue
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_setVectorInt (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  INTEGER(I32), INTENT(IN) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  INTEGER(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = ivalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
  DO i = 1,SIZE(Ix)
    Ix(i) = ivalue
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_setVectorDble2D (Dx,dvalue)
  
!<description>
  ! Sets the vector data to a defined value: Dx = dvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  REAL(DP), INTENT(IN) :: dvalue
!</input>

!<output>
  ! Destination vector to be set
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  INTEGER(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = dvalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  DO j = 1,SIZE(Dx,2)
    DO i = 1,SIZE(Dx,1)
      Dx(i,j) = dvalue
    END DO
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_setVectorSngl2D (Fx,fvalue)
  
!<description>
  ! Sets the vector data to a defined value: Fx = fvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  REAL(SP), INTENT(IN) :: fvalue
!</input>

!<output>
  ! Destination vector to be set
  REAL(SP), DIMENSION(:,:), INTENT(OUT) :: Fx
!</output>

!</subroutine>

  ! local variables
  INTEGER(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Fx = fvalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  DO j = 1,SIZE(Fx,2)
    DO i = 1,SIZE(Fx,1)
      Fx(i,j) = fvalue
    END DO
  END DO
!%OMP  end parallel do
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_setVectorInt2D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  INTEGER(I32), INTENT(IN) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  INTEGER(I32), DIMENSION(:,:), INTENT(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  INTEGER(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = ivalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  DO j = 1,SIZE(Ix,2)
    DO i = 1,SIZE(Ix,1)
      Ix(i,j) = ivalue
    END DO
  END DO
!%OMP  end parallel do
  
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
  ELSE IF (dcy .EQ. 1.0_DP) THEN
    CALL DAXPY(SIZE(Dx),dcx,Dx,1,Dy,1)
  ELSE
    c=dcx/dcy
    CALL DAXPY(SIZE(Dx),c,Dx,1,Dy,1)
    CALL DSCAL(SIZE(Dx),dcy,Dy,1)
  ENDIF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorLinearCombSngl (Fx,Fy,scx,scy)
  
!<description>
  ! Performs a linear combination: Fy = scx * Fx  +  scy * Fy
!</description>

!<input>
  
  ! First source vector
  REAL(SP), DIMENSION(:), INTENT(IN) :: Fx
  
  ! Scaling factor for Dx
  REAL(SP), INTENT(IN)               :: scx

  ! Scaling factor for Dy
  REAL(SP), INTENT(IN)               :: scy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: Fy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(SP) :: c
  
  IF (scy .EQ. 0.0_SP) THEN
    CALL SCOPY(SIZE(Fx),Fx,1,Fy,1)
    IF (scx .NE. 1.0_SP) CALL SSCAL(SIZE(Fx),scx,Fy,1)
  ELSE IF (scy .EQ. 1.0_DP) THEN
    CALL SAXPY(SIZE(Fx),scx,Fx,1,Fy,1)
  ELSE
    c=scx/scy
    CALL SAXPY(SIZE(Fx),c,Fx,1,Fy,1)
    CALL SSCAL(SIZE(Fx),scy,Fy,1)
  ENDIF
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorLinearCombSnglDble (Fx,Dy,scx,dcy)
  
!<description>
  ! Performs a linear combination: Dy = scx * Fx  +  dcy * Dy
!</description>

!<input>
  
  ! First source vector
  REAL(SP), DIMENSION(:), INTENT(IN) :: Fx
  
  ! Scaling factor for Dx
  REAL(SP), INTENT(IN)               :: scx

  ! Scaling factor for Dy
  REAL(DP), INTENT(IN)               :: dcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: Dy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER(I32) :: i
  REAL(DP) :: c
  
  IF (dcy .EQ. 0.0_DP) THEN
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
    DO i=1,SIZE(Fx)
      Dy(i) = Fx(i)
    END DO
!%OMP  end parallel do
    IF (scx .NE. 1.0_SP) THEN
!%OMP  parallel do &
!%OMP& default(shared) & 
!%OMP& private(i)
      DO i=1,SIZE(Fx)
        Dy(i) = scx*Fx(i)
      END DO
!%OMP  end parallel do
    END IF
  ELSE IF (dcy .EQ. 1.0_DP) THEN
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
    DO i=1,SIZE(Fx)
      Dy(i) = Dy(i) + scx*Fx(i)
    END DO
!%OMP  end parallel do
  ELSE
    c=scx/dcy
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
    DO i=1,SIZE(Fx)
      Dy(i) = Dy(i) + c*Fx(i)
    END DO
!%OMP  end parallel do
    CALL DSCAL(SIZE(Dy),dcy,Dy,1)
  ENDIF
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorLinearCombDble2D (Dx,Dy,dcx,dcy)
  
!<description>
  ! Performs a linear combination: Dy = dcx * Dx  +  dcy * Dy
!</description>

!<input>
  
  ! First source vector
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dx
  
  ! Scaling factor for Dx
  REAL(DP), INTENT(IN)                 :: dcx

  ! Scaling factor for Dy
  REAL(DP), INTENT(IN)                 :: dcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: Dy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP) :: c
  
  IF (dcy .EQ. 0.0_DP) THEN
    CALL DCOPY(SIZE(Dx,1)*SIZE(Dx,2),Dx,1,Dy,1)
    IF (dcx .NE. 1.0_DP) CALL DSCAL(SIZE(Dx,1)*SIZE(Dx,2),dcx,Dy,1)
  ELSE IF (dcy .EQ. 1.0_DP) THEN
    CALL DAXPY(SIZE(Dx,1)*SIZE(Dx,2),dcx,Dx,1,Dy,1)
  ELSE
    c=dcx/dcy
    CALL DAXPY(SIZE(Dx,1)*SIZE(Dx,2),c,Dx,1,Dy,1)
    CALL DSCAL(SIZE(Dx,1)*SIZE(Dx,2),dcy,Dy,1)
  ENDIF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorLinearCombSngl2D (Fx,Fy,scx,scy)
  
!<description>
  ! Performs a linear combination: Fy = scx * Fx  +  scy * Fy
!</description>

!<input>
  
  ! First source vector
  REAL(SP), DIMENSION(:,:), INTENT(IN) :: Fx
  
  ! Scaling factor for Dx
  REAL(SP), INTENT(IN)                 :: scx

  ! Scaling factor for Dy
  REAL(SP), INTENT(IN)                 :: scy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: Fy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(SP) :: c
  
  IF (scy .EQ. 0.0_SP) THEN
    CALL SCOPY(SIZE(Fx,1)*SIZE(Fx,2),Fx,1,Fy,1)
    IF (scx .NE. 1.0_SP) CALL SSCAL(SIZE(Fx,1)*SIZE(Fx,2),scx,Fy,1)
  ELSE IF (scy .EQ. 1.0_DP) THEN
    CALL SAXPY(SIZE(Fx,1)*SIZE(Fx,2),scx,Fx,1,Fy,1)
  ELSE
    c=scx/scy
    CALL SAXPY(SIZE(Fx,1)*SIZE(Fx,2),c,Fx,1,Fy,1)
    CALL SSCAL(SIZE(Fx,1)*SIZE(Fx,2),scy,Fy,1)
  ENDIF
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorLinearCombSnglDble2D (Fx,Dy,scx,dcy)
  
!<description>
  ! Performs a linear combination: Dy = scx * Fx  +  dcy * Dy
!</description>

!<input>
  
  ! First source vector
  REAL(SP), DIMENSION(:,:), INTENT(IN) :: Fx
  
  ! Scaling factor for Dx
  REAL(SP), INTENT(IN)                 :: scx

  ! Scaling factor for Dy
  REAL(DP), INTENT(IN)                 :: dcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: Dy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER(I32) :: i,j
  REAL(DP) :: c
  
  IF (dcy .EQ. 0.0_DP) THEN
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
    DO j=1,SIZE(Fx,2)
      DO i=1,SIZE(Fx,1)
        Dy(i,j) = Fx(i,j)
      END DO
    END DO
!%OMP  end parallel do
    IF (scx .NE. 1.0_SP) THEN
!%OMP  parallel do &
!%OMP& default(shared) & 
!%OMP& private(i,j)
      DO j=1,SIZE(Fx,2)
        DO i=1,SIZE(Fx,1)
          Dy(i,j) = scx*Fx(i,j)
        END DO
      END DO
!%OMP  end parallel do
    END IF
  ELSE IF (dcy .EQ. 1.0_DP) THEN
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
    DO j=1,SIZE(Fx,2)
      DO i=1,SIZE(Fx,1)
        Dy(i,j) = Dy(i,j) + scx*Fx(i,j)
      END DO
    END DO
!%OMP  end parallel do
  ELSE
    c=scx/dcy
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
    DO j=1,SIZE(Fx,2)
      DO i=1,SIZE(Fx,1)
        Dy(i,j) = Dy(i,j) + c*Fx(i,j)
      END DO
    END DO
!%OMP  end parallel do
    CALL DSCAL(SIZE(Dy,1)*SIZE(Dy,2),dcy,Dy,1)
  ENDIF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<function>

  REAL(DP) FUNCTION lalg_scalarProductDble (Dx,Dy) RESULT (res)
  
!<description>
  ! Calculates the scalar product of two double precision vectors: 
  ! res = (vector,vector)
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
  !INTEGER(I32) :: i
  REAL(DP) :: DDOT

  res=DDOT(SIZE(Dx),Dx,1,Dy,1)
!!$  res = Dx(1)*Dy(1)
!!$  DO i=2,SIZE(Dx)
!!$    res = res + Dx(i)*Dy(i)
!!$  END DO
  
  END FUNCTION

  ! ***************************************************************************

!<function>

  REAL(SP) FUNCTION lalg_scalarProductSngl (Fx,Fy) RESULT (res)
  
!<description>
  ! Calculates the scalar product of two single precision vectors: 
  ! res = (vector,vector)
!</description>

!<input>
  
  ! First source vector
  REAL(SP), DIMENSION(:), INTENT(IN) :: Fx
  
  ! Second source vector
  REAL(SP), DIMENSION(:), INTENT(IN) :: Fy
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  ! local variables
  !INTEGER(I32) :: i
  REAL(SP) :: SDOT

  res=SDOT(SIZE(Fx),Fx,1,Fy,1)
!!$  res = Fx(1)*Fy(1)
!!$  DO i=2,SIZE(Fx)
!!$    res = res + Fx(i)*Fy(i)
!!$  END DO
  
  END FUNCTION

  ! ***************************************************************************

!<function>

  REAL(SP) FUNCTION lalg_scalarProductInt (Ix,Iy) RESULT (res)
  
!<description>
  ! Calculates the scalar product of two single precision vectors: 
  ! res = (vector,vector)
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

!<function>

  REAL(DP) FUNCTION lalg_scalarProductDble2D (Dx,Dy) RESULT (res)
  
!<description>
  ! Calculates the scalar product of two double precision vectors: 
  ! res = (vector,vector)
!</description>

!<input>
  
  ! First source vector
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dx
  
  ! Second source vector
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dy
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  ! local variables
  !INTEGER(I32) :: i
  REAL(DP) :: DDOT

  res=DDOT(SIZE(Dx,1)*SIZE(Dx,2),Dx,1,Dy,1)
!!$  res = Dx(1)*Dy(1)
!!$  DO i=2,SIZE(Dx)
!!$    res = res + Dx(i)*Dy(i)
!!$  END DO
  
  END FUNCTION

  ! ***************************************************************************

!<function>

  REAL(SP) FUNCTION lalg_scalarProductSngl2D (Fx,Fy) RESULT (res)
  
!<description>
  ! Calculates the scalar product of two single precision vectors: 
  ! res = (vector,vector)
!</description>

!<input>
  
  ! First source vector
  REAL(SP), DIMENSION(:,:), INTENT(IN) :: Fx
  
  ! Second source vector
  REAL(SP), DIMENSION(:,:), INTENT(IN) :: Fy
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  ! local variables
  !INTEGER(I32) :: i
  REAL(SP) :: SDOT

  res=SDOT(SIZE(Fx,1)*SIZE(Fx,2),Fx,1,Fy,1)
!!$  res = Fx(1)*Fy(1)
!!$  DO i=2,SIZE(Fx)
!!$    res = res + Fx(i)*Fy(i)
!!$  END DO
  
  END FUNCTION

  ! ***************************************************************************

!<function>

  REAL(SP) FUNCTION lalg_scalarProductInt2D (Ix,Iy) RESULT (res)
  
!<description>
  ! Calculates the scalar product of two single precision vectors: 
  ! res = (vector,vector)
!</description>

!<input>
  
  ! First source vector
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: Ix
  
  ! Second source vector
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: Iy
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  ! local variables
  INTEGER(I32) :: i,j
  
  res = Ix(1,1)*Iy(1,1)
  DO j=1,SIZE(Ix,2)
    DO i=2,SIZE(Ix,1)
      res = res + Ix(i,j)*Iy(i,j)
    END DO
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
  ! < 0 >1, if an error occurred (unknown norm).
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
    ! Euclidian norm = scalar product (vector,vector)
    resnorm = lalg_scalarProductDble(Dx,Dx)
!!$    resnorm = Dx(1)*Dx(1)
!!$    DO i=2,SIZE(Dx)
!!$      resnorm = resnorm + Dx(i)*Dx(i)
!!$    END DO
    resnorm = SQRT(resnorm)

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
    resnorm = lalg_scalarProductDble(Dx,Dx)
!!$    resnorm = Dx(1)*Dx(1)
!!$    DO i=2,SIZE(Dx)
!!$      resnorm = resnorm + Dx(i)*Dx(i)
!!$    END DO
    resnorm = SQRT(resnorm / REAL(SIZE(Dx),DP))
    
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

  REAL(SP) FUNCTION lalg_normSngl (Fx,cnorm,iposMax) RESULT(resnorm)
  
!<description>
  ! Calculates the norm of a single precision vector. cnorm identifies the 
  ! type of norm to calculate.
!</description>

!<input>
  ! Vector to calculate the norm of.
  REAL(SP), DIMENSION(:), INTENT(IN) :: Fx
  
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
    resnorm = Fx(1)
    DO i=2,SIZE(Fx)
      resnorm = resnorm + ABS(Fx(i))
    END DO

  CASE (LINALG_NORMEUCLID)
    ! Euclidian norm = scalar product (vector,vector)
    resnorm = lalg_scalarProductSngl(Fx,Fx)
!!$    resnorm = Fx(1)*Fx(1)
!!$    DO i=2,SIZE(Fx)
!!$      resnorm = resnorm + Fx(i)*Fx(i)
!!$    END DO
    resnorm = SQRT(resnorm)

  CASE (LINALG_NORML1)
    ! L1-norm: sum sum absolute value of all entries, divide by sqrt(vector length).
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = Fx(1)
    DO i=2,SIZE(Fx)
      resnorm = resnorm + ABS(Fx(i))
    END DO
    resnorm = resnorm / REAL(SIZE(Fx),SP)

  CASE (LINALG_NORML2)
    ! l_2-norm - like euclidian norm, but divide by sqrt(vector length).
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = lalg_scalarProductSNGL(Fx,Fx)
!!$    resnorm = Fx(1)*Fx(1)
!!$    DO i=2,SIZE(Fx)
!!$      resnorm = resnorm + Fx(i)*Fx(i)
!!$    END DO
    resnorm = SQRT(resnorm / REAL(SIZE(Fx),SP))
    
  CASE (LINALG_NORMMAX)
    ! MAX-norm. Find the absolute largest entry.
    resnorm = ABS(Fx(1))
    j=1
    DO i=2,SIZE(Fx)
      IF (ABS(Fx(i)) .GT. resnorm) THEN
        j = i
        resnorm = ABS(Fx(i))
      END IF
    END DO
    IF (PRESENT(iposMax)) iposMax = j
  CASE DEFAULT
    resnorm = -1.0_SP ! Unknown norm.
  END SELECT
    
  END FUNCTION

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorSortDble (Dx, Dd, Itr)
  
  !<description>
    ! Resorts the entries in the vector Dx corresponding to Itr.
    ! In particular, the first SIZE(Itr) entries of Dx are written resortedly
    ! to Dd.
  !</description>
    
  !<input>
    ! Source vector to be sorted
    REAL(DP), DIMENSION(:), INTENT(IN) :: Dx
    
    ! Array with permutation of 1..neq.
    ! Itr(i) defines the number of the entry in Dx that should
    ! move to position i.
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: Itr
  !</input>
    
  !<output>
    ! The resorted vector
    REAL(DP), DIMENSION(:), INTENT(OUT) :: Dd
  !</output>
    
!</subroutine>
    
    ! local variable
    INTEGER(I32) :: ieq
    
    ieq = SIZE(Itr)
    ieq = SIZE(Dx)
    ieq = SIZE(Dd)
    
    DO ieq=1, SIZE(Itr)
      Dd(ieq) = Dx(Itr(ieq))
    END DO
  
  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorSortSngl (Fx, Fd, Itr)
  
  !<description>
    ! Resorts the entries in the vector Fx corresponding to Itr.
    ! In particular, the first SIZE(Itr) entries of Dx are written resortedly
    ! to Dd.
  !</description>
    
  !<input>
    ! Array with permutation of 1..neq
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: Itr

    ! Source vector to be sorted
    REAL(SP), DIMENSION(:), INTENT(IN) :: Fx
  !</input>
    
  !<output>
    ! The resorted vector
    REAL(SP), DIMENSION(:), INTENT(OUT) :: Fd
  !</output>
    
!</subroutine>
    
    ! local variable
    INTEGER(I32) :: ieq
    
    DO ieq=1, SIZE(Itr)
      Fd(ieq) = Fx(Itr(ieq))
    END DO
  
  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorSortInt (Ix, Id, Itr)
  
  !<description>
    ! Resorts the entries in the vector Ix corresponding to Itr.
    ! In particular, the first SIZE(Itr) entries of Dx are written resortedly
    ! to Dd.
  !</description>
    
  !<input>
    ! Array with permutation of 1..neq
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: Itr

    ! Source vector to be sorted
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: Ix
  !</input>
    
  !<output>
    ! The resorted vector
    INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Id
  !</output>
    
!</subroutine>
    
    ! local variable
    INTEGER(I32) :: ieq
    
    DO ieq=1, SIZE(Itr)
      Id(ieq) = Ix(Itr(ieq))
    END DO
  
  END SUBROUTINE 

!<subroutine>

  SUBROUTINE lalg_tensorProductDble(Dx,Dy,Dtensor)

!<description>
    ! Calculates the tensor product of two double precision vectors:
    ! Dtensor = Dx (*) Dy
!</description>

!<input>
    ! First source vector
    REAL(DP), DIMENSION(:), INTENT(IN) :: Dx
  
    ! Second source vector
    REAL(DP), DIMENSION(:), INTENT(IN) :: Dy
!</input>

!<output>
    ! Tensor product
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dtensor
!</output>
!</subroutine>
    
    ! local variables
    INTEGER :: i,j

    DO i=1,SIZE(Dy)
      DO j=1,SIZE(Dx)
        Dtensor(j,i)=Dx(j)*Dy(i)
      END DO
    END DO
  END SUBROUTINE lalg_tensorProductDble

  !****************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorAddScalarDble (Dx,dvalue)
  
!<description>
  ! This routine adds the value dvalue to each entry of the vector Id.
!</description>
  
!<input>
  ! The value to add to every entry.
  REAL(DP) :: dvalue
!</input>

!<inputoutput>
  ! Data array to be modified.
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: Dx
!</inputoutput>

!</subroutine>
    
    REAL(DP) :: dval
    INTEGER(I32) :: i
    
    dval = dvalue
    DO i=1,SIZE(Dx)
      Dx(i) = Dx(i) + dval
    END DO
    
  END SUBROUTINE
  
  !****************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorAddScalarSngl (Fx,fvalue)
  
!<description>
  ! This routine adds the value dvalue to each entry of the vector Fx.
!</description>
  
!<input>
  ! The value to add to every entry.
  REAL(SP) :: fvalue
!</input>

!<inputoutput>
  ! Data array to be modified.
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: Fx
!</inputoutput>

!</subroutine>
    
    REAL(SP) :: fval
    INTEGER(I32) :: i
    
    fval = fvalue
    DO i=1,SIZE(Fx)
      Fx(i) = Fx(i) + fval
    END DO
    
  END SUBROUTINE
  
  !****************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorAddScalarInt (Ix,ivalue)
  
!<description>
  ! This routine adds the value ivalue to each entry of the vector Ix.
!</description>
  
!<input>
  ! The value to add to every entry.
  INTEGER(I32) :: ivalue
!</input>

!<inputoutput>
  ! Data array to be modified.
  INTEGER(I32), DIMENSION(:), INTENT(INOUT) :: Ix
!</inputoutput>

!</subroutine>
    
    INTEGER(I32) :: ival
    INTEGER(I32) :: i
    
    ival = ivalue
    DO i=1,SIZE(Ix)
      Ix(i) = Ix(i) + ival
    END DO
    
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorAddScalarDble2D (Dx,dvalue)
  
!<description>
  ! This routine adds the value dvalue to each entry of the vector Id.
!</description>
  
!<input>
  ! The value to add to every entry.
  REAL(DP) :: dvalue
!</input>

!<inputoutput>
  ! Data array to be modified.
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: Dx
!</inputoutput>

!</subroutine>
    
    REAL(DP) :: dval
    INTEGER(I32) :: i,j
    
    dval = dvalue
    DO j=1,SIZE(Dx,2)
      DO i=1,SIZE(Dx,1)
        Dx(i,j) = Dx(i,j) + dval
      END DO
    END DO
    
  END SUBROUTINE
  
  !****************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorAddScalarSnglD (Fx,fvalue)
  
!<description>
  ! This routine adds the value dvalue to each entry of the vector Fx.
!</description>
  
!<input>
  ! The value to add to every entry.
  REAL(SP) :: fvalue
!</input>

!<inputoutput>
  ! Data array to be modified.
  REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: Fx
!</inputoutput>

!</subroutine>
    
    REAL(SP) :: fval
    INTEGER(I32) :: i,j
    
    fval = fvalue
    DO j=1,SIZE(Fx,2)
      DO i=1,SIZE(Fx,1)
      Fx(i,j) = Fx(i,j) + fval
    END DO
    END DO
    
  END SUBROUTINE
  
  !****************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorAddScalarInt2D (Ix,ivalue)
  
!<description>
  ! This routine adds the value ivalue to each entry of the vector Ix.
!</description>
  
!<input>
  ! The value to add to every entry.
  INTEGER(I32) :: ivalue
!</input>

!<inputoutput>
  ! Data array to be modified.
  INTEGER(I32), DIMENSION(:,:), INTENT(INOUT) :: Ix
!</inputoutput>

!</subroutine>
    
    INTEGER(I32) :: ival
    INTEGER(I32) :: i,j
    
    ival = ivalue
    DO j=1,SIZE(Ix,2)
      DO i=1,SIZE(Ix,1)
        Ix(i,j) = Ix(i,j) + ival
      END DO
    END DO
    
  END SUBROUTINE
  
END MODULE

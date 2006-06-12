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
!# 7.) lalg_vectorSort
!#     -> Resort the entries of a vector according to a given permutation
!#
!# </purpose>
!##############################################################################

MODULE linearalgebra

  USE fsystem

  IMPLICIT NONE
  
!<constants>

!<constantblock description="Constants identifying vector norms">

  ! Sum of the absolute values of entries
  INTEGER, PARAMETER :: LINALG_NORMSUM    = -1

  ! Euclidian vector norm: <vector,vector>
  INTEGER, PARAMETER :: LINALG_NORMEUCLID = 0

  ! $l_1$-norm: 1/NEQ * sum(abs(entries))
  INTEGER, PARAMETER :: LINALG_NORML1     = 1
  
  ! $l_2$-norm: 1/sqrt(SEQ) * <vector,vector>
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
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dy
  
!</output>
  
!</subroutine>

  CALL DCOPY(SIZE(Dx),Dx,1,Dy,1)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorCopySngl (Fx,Fy)
  
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

  SUBROUTINE lalg_vectorCopySnglDbl (Fx,Dy)
  
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
  
  DO i=1,SIZE(Fx)
    Dy(i) = REAL(Fx(i),DP)
  END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lalg_vectorCopyDblSngl (Dx,Fy)
  
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
  
  DO i=1,SIZE(Dx)
    Fy(i) = REAL(Dx(i),SP)
  END DO

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

  SUBROUTINE lalg_vectorScaleSngl (Fx,sc)
  
!<description>
  ! Scales a single precision vector: Dx = sc * Dx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: Fx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  REAL(DP), INTENT(IN) :: sc

!</input>
  
!</subroutine>

  CALL SSCAL(SIZE(Fx),sc,Fx,1)
  
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

  SUBROUTINE lalg_vectorClearSngl (Fx)
  
!<description>
  ! Clears a single precision vector: Fx = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  REAL(SP), DIMENSION(:), INTENT(OUT) :: Fx
  
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

  DO i = 1,SIZE(Fx)
    Fx(i) = 0.0_SP
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
  ELSE IF (scy .EQ. 1D0) THEN
    CALL SAXPY(SIZE(Fx),scx,Fx,1,Fy,1)
  ELSE
    c=scx/scy
    CALL SAXPY(SIZE(Fx),c,Fx,1,Fy,1)
    CALL SSCAL(SIZE(Fx),scy,Fy,1)
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

  REAL(SP) FUNCTION lalg_scalarProductSngl (Fx,Fy) RESULT (res)
  
!<description>
  ! Calculates the scalar product of two single precision vectors: 
  ! res = <vector,vector>
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
  INTEGER(I32) :: i
  
  res = Fx(1)*Fy(1)
  DO i=2,SIZE(Fx)
    res = res + Fx(i)*Fy(i)
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
    resnorm = DSQRT(resnorm)

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
    resnorm = DSQRT(resnorm / REAL(SIZE(Dx),DP))
    
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
    ! Euclidian norm = scalar product <vector,vector>
    resnorm = Fx(1)*Fx(1)
    DO i=2,SIZE(Fx)
      resnorm = resnorm + Fx(i)*Fx(i)
    END DO
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
    resnorm = Fx(1)*Fx(1)
    DO i=2,SIZE(Fx)
      resnorm = resnorm + Fx(i)*Fx(i)
    END DO
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

END MODULE

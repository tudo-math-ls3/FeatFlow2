!##############################################################################
!# ****************************************************************************
!# <name> mprimitives </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of low-level mathematic functions without
!# a big background of finite elements or similar, like linear transformation,
!# parabolic profile, etc.
!#
!# The routines here are usually declared as PURE as they should not have
!# any side effects!
!#
!# The following routines can be found here:
!#
!# 1.) mprim_getParabolicProfile
!#     -> Calculates the value of a parabolic profile along a line.
!#
!# 2.) mprim_invertMatrix
!#     -> Invert a full matrix
!#
!# 3.) mprim_kronecker
!#     -> Compute Kronecker delta symbol
!#
!# 4.) mprim_invert4x4MatrixDirectDble
!#     -> Inverts a 4x4 matrix directly without pivoting.
!#
!# 5.) mprim_invertMatrixPivotDble
!#     -> Inverts a 4x4 matrix directly with pivoting.
!#
!# 6.) mprim_signum
!#     -> Signum function
!# 
!# 7.) mprim_linearRescale
!#     -> Scales a coordinate x linearly from the interval [a,b] to the
!#        interval [c,d]
!#
!# 8.) mprim_quadraticInterpolation
!#     -> Evaluate the quadratic interpolation polynomial of three values.
!#
!# 9.) mprim_SVD_factorise
!#     -> Compute the factorisation for a singular value decomposition
!#
!# 10.) mprim_SVD_backsubst
!#      -> Perform back substitution for a singular value decomposition
!#
!# </purpose>
!##############################################################################

MODULE mprimitives

  USE fsystem
  USE genoutput
  
  IMPLICIT NONE
  
  INTERFACE mprim_signum
    MODULE PROCEDURE mprim_signum_dble
    MODULE PROCEDURE mprim_signum_real
    MODULE PROCEDURE mprim_signum_int32
  END INTERFACE

  ! Alternative name for backward compatibility
  INTERFACE mprim_kronecker
    MODULE PROCEDURE kronecker
  END INTERFACE

CONTAINS

  ! ***************************************************************************

!<function>
  
  ELEMENTAL REAL(DP) FUNCTION mprim_getParabolicProfile (dpos,dlength,dmaxvalue)
  
!<description>
  ! Calculates the value of a parabolic profile along a line.
  ! dpos is a parameter value along a line of length dlength. The return value
  ! is the value of the profile and has its maximum dmax at position 0.5.
!</description>
  
!<input>
  
  ! Position on the line segment where to calculate the velocity
  REAL(DP), INTENT(IN) :: dpos
  
  ! Length of the line segment
  REAL(DP), INTENT(IN) :: dlength
  
  ! Maximum value of the profile 
  REAL(DP), INTENT(IN) :: dmaxvalue
  
!</input>

!<result>
  ! The value of the parabolic profile at position $dpos \in [0,dmax]$.
!</result>
!</function>

    ! Result: Value of the parbolic profile on position dpos on the line segment
    mprim_getParabolicProfile = 4.0_DP*dmaxvalue*dpos*(dlength-dpos)/(dlength*dlength)
    
  END FUNCTION

  ! ***************************************************************************

!<function>
  
  ELEMENTAL REAL(DP) FUNCTION mprim_signum_dble (dval)
  
!<description>
  ! Signum function
!</description>
  
!<input>
  ! Value to be checked
  REAL(DP), INTENT(IN) :: dval
!</input>

!<result>
  ! The value of the parabolic profile at position $dpos \in [0,dmax]$.
!</result>
!</function>

    ! Result: Value of the parbolic profile on position dpos on the line segment
    IF (dval .LT. 0.0_DP) THEN
      mprim_signum_dble = -1.0_DP
    ELSE IF (dval .GT. 0.0_DP) THEN
      mprim_signum_dble = 1.0_DP
    ELSE 
      mprim_signum_dble = 0.0_DP
    END IF

  END FUNCTION

  ! ***************************************************************************

!<function>
  
  ELEMENTAL REAL(SP) FUNCTION mprim_signum_real (fval)
  
!<description>
  ! Signum function.
!</description>
  
!<input>
  ! Value to be checked
  REAL(SP), INTENT(IN) :: fval
!</input>

!<result>
  ! The value of the parabolic profile at position $dpos \in [0,dmax]$.
!</result>
!</function>

    ! Result: Value of the parbolic profile on position dpos on the line segment
    IF (fval .LT. 0.0_SP) THEN
      mprim_signum_real = -1.0_SP
    ELSE IF (fval .GT. 0.0_SP) THEN
      mprim_signum_real = 1.0_SP
    ELSE 
      mprim_signum_real = 0.0_SP
    END IF

  END FUNCTION

  ! ***************************************************************************

!<function>
  
  ELEMENTAL INTEGER(I32) FUNCTION mprim_signum_int32 (ival)
  
!<description>
  ! Signum function.
!</description>
  
!<input>
  ! Value to be checked
  INTEGER(I32), INTENT(IN) :: ival
!</input>

!<result>
  ! The value of the parabolic profile at position $dpos \in [0,dmax]$.
!</result>
!</function>

    ! Result: Value of the parbolic profile on position dpos on the line segment
    SELECT CASE (ival)
      CASE (:-1)
        mprim_signum_int32 = -1_I32
      CASE (0)
        mprim_signum_int32 = 0_I32
      CASE DEFAULT
        mprim_signum_int32 = 1_I32
    END SELECT

  END FUNCTION

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mprim_invertMatrixDble(Da,Df,Dx,ndim,ipar)

!<description>
    ! This subroutine performs the direct inversion of a NxN system.
    ! 
    ! If the parameter ipar=0, then only factorization of matrix A is performed. 
    ! For ipar=1, the vector x is calculated using the factorized matrix A. 
    ! For ipar=2, LAPACK routine DGESV is used to solve the dense linear system Ax=f. 
    ! In addition, for NDIM=2,3,4 explicit formulas are employed to replace the
    ! more expensive LAPACK routine.
!</description>

!<input>
    ! dimension of the matrix
    INTEGER, INTENT(IN) :: ndim

    ! source right-hand side vector
    REAL(DP), DIMENSION(ndim), INTENT(IN) :: Df

    ! What to do?
    ! IPAR = 0 : invert matrix by means of Gaussian elimination with
    !            full pivoting and return the inverted matrix inv(A)
    ! IPAR = 1 : apply inverted matrix to the right-hand side vector
    !            and return x = inv(A)*f
    ! IPAR = 2 : invert matrix and apply it to the right-hand side
    !            vector. Return x = inv(A)*f
    INTEGER, INTENT(IN) :: ipar
!</input>

!<inputoutput> 
    ! source square matrix to be inverted
    REAL(DP), DIMENSION(ndim,ndim), INTENT(INOUT) :: Da
!</inputoutput>

!<output>
    ! destination vector containing inverted matrix times right-hand
    ! side vector
    REAL(DP), DIMENSION(ndim), INTENT(OUT) :: Dx
!</output>

!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(ndim,ndim) :: Db
    REAL(DP), DIMENSION(ndim) :: Dpiv
    INTEGER, DIMENSION(ndim) :: Kindx,Kindy

    REAL(DP) :: dpivot,daux
    INTEGER :: idim1,idim2,ix,iy,indx,indy,info

    SELECT CASE (ipar)
    CASE (0)
      ! Perform factorization of matrix Da

      ! Initialization
      Kindx=0;  Kindy=0

      DO idim1=1,ndim

        ! Determine pivotal element
        dpivot=0

        DO iy=1,ndim
          IF (Kindy(iy) /= 0) CYCLE

          DO ix=1,ndim
            IF (Kindx(ix) /= 0) CYCLE

            IF (ABS(Da(ix,iy)) .LE. ABS(dpivot)) CYCLE
            dpivot=Da(ix,iy);  indx=ix;  indy=iy
          END DO
        END DO

        ! Return if pivotal element is zero
        IF (ABS(dpivot) .LE. 0._DP) RETURN

        Kindx(indx)=indy;  Kindy(indy)=indx;  Da(indx,indy)=1._DP&
            &/dpivot

        DO idim2=1,ndim
          IF (idim2 == indy) CYCLE 
          Da(1:indx-1,idim2)=Da(1:indx-1,idim2)-Da(1:indx-1,  &
              & indy)*Da(indx,idim2)/dpivot
          Da(indx+1:ndim,idim2)=Da(indx+1:ndim,idim2)-Da(indx+1:ndim&
              &,indy)*Da(indx,idim2)/dpivot
        END DO

        DO ix=1,ndim
          IF (ix /= indx) Da(ix,indy)=Da(ix,indy)/dpivot
        END DO

        DO iy=1,ndim
          IF (iy /= indy) Da(indx,iy)=-Da(indx,iy)/dpivot
        END DO
      END DO

      DO ix=1,ndim
        IF (Kindx(ix) == ix) CYCLE

        DO iy=1,ndim
          IF (Kindx(iy) == ix) EXIT
        END DO

        DO idim1=1,ndim
          daux=Da(ix,idim1)
          Da(ix,idim1)=Da(iy,idim1)
          Da(iy,idim1)=daux
        END DO

        Kindx(iy)=Kindx(ix);  Kindx(ix)=ix
      END DO

      DO ix=1,ndim
        IF (Kindy(ix) == ix) CYCLE

        DO iy=1,ndim
          IF (Kindy(iy) == ix) EXIT
        END DO

        DO idim1=1,ndim
          daux=Da(idim1,ix)
          Da(idim1,ix)=Da(idim1,iy)
          Da(idim1,iy)=daux
        END DO
        
        Kindy(iy)=Kindy(ix);  Kindy(ix)=ix
      END DO


    CASE (1)    
      ! Perform inversion of Da to solve the system Da * Dx = Df
      DO idim1=1,ndim
        Dx(idim1)=0
        DO idim2=1,ndim
          Dx(idim1)=Dx(idim1)+Da(idim1,idim2)*Df(idim2)
        END DO
      END DO

    CASE (2)
      ! Solve the dense linear system Ax=f calling LAPACK routine

      SELECT CASE(ndim)
      CASE (2)
        ! Explicit formula for 2x2 system
        Db(1,1)= Da(2,2)
        Db(2,1)=-Da(2,1)
        Db(1,2)=-Da(1,2)
        Db(2,2)= Da(1,1)
        daux=Da(1,1)*Da(2,2)-Da(1,2)*Da(2,1)
        Dx=MATMUL(Db,Df)/daux

      CASE (3)
        ! Explicit formula for 3x3 system
        Db(1,1)=Da(2,2)*Da(3,3)-Da(2,3)*Da(3,2)
        Db(2,1)=Da(2,3)*Da(3,1)-Da(2,1)*Da(3,3)
        Db(3,1)=Da(2,1)*Da(3,2)-Da(2,2)*Da(3,1)
        Db(1,2)=Da(1,3)*Da(3,2)-Da(1,2)*Da(3,3)
        Db(2,2)=Da(1,1)*Da(3,3)-Da(1,3)*Da(3,1)
        Db(3,2)=Da(1,2)*Da(3,1)-Da(1,1)*Da(3,2)
        Db(1,3)=Da(1,2)*Da(2,3)-Da(1,3)*Da(2,2)
        Db(2,3)=Da(1,3)*Da(2,1)-Da(1,1)*Da(2,3)
        Db(3,3)=Da(1,1)*Da(2,2)-Da(1,2)*Da(2,1)
        daux=Da(1,1)*Da(2,2)*Da(3,3)+Da(2,1)*Da(3,2)*Da(1,3)+ Da(3,1)&
            &*Da(1,2)*Da(2,3)-Da(1,1)*Da(3,2)*Da(2,3)- Da(3,1)*Da(2&
            &,2)*Da(1,3)-Da(2,1)*Da(1,2)*Da(3,3)
        Dx=MATMUL(Db,Df)/daux

      CASE (4)
        ! Explicit formula for 4x4 system
        Db(1,1)=Da(2,2)*Da(3,3)*Da(4,4)+Da(2,3)*Da(3,4)*Da(4,2)+Da(2&
            &,4)*Da(3,2)*Da(4,3)- Da(2,2)*Da(3,4)*Da(4,3)-Da(2,3)&
            &*Da(3,2)*Da(4,4)-Da(2,4)*Da(3,3)*Da(4,2)
        Db(2,1)=Da(2,1)*Da(3,4)*Da(4,3)+Da(2,3)*Da(3,1)*Da(4,4)+Da(2&
            &,4)*Da(3,3)*Da(4,1)- Da(2,1)*Da(3,3)*Da(4,4)-Da(2,3)&
            &*Da(3,4)*Da(4,1)-Da(2,4)*Da(3,1)*Da(4,3)
        Db(3,1)=Da(2,1)*Da(3,2)*Da(4,4)+Da(2,2)*Da(3,4)*Da(4,1)+Da(2&
            &,4)*Da(3,1)*Da(4,2)- Da(2,1)*Da(3,4)*Da(4,2)-Da(2,2)&
            &*Da(3,1)*Da(4,4)-Da(2,4)*Da(3,2)*Da(4,1)
        Db(4,1)=Da(2,1)*Da(3,3)*Da(4,2)+Da(2,2)*Da(3,1)*Da(4,3)+Da(2&
            &,3)*Da(3,2)*Da(4,1)- Da(2,1)*Da(3,2)*Da(4,3)-Da(2,2)&
            &*Da(3,3)*Da(4,1)-Da(2,3)*Da(3,1)*Da(4,2)
        Db(1,2)=Da(1,2)*Da(3,4)*Da(4,3)+Da(1,3)*Da(3,2)*Da(4,4)+Da(1&
            &,4)*Da(3,3)*Da(4,2)- Da(1,2)*Da(3,3)*Da(4,4)-Da(1,3)&
            &*Da(3,4)*Da(4,2)-Da(1,4)*Da(3,2)*Da(4,3)
        Db(2,2)=Da(1,1)*Da(3,3)*Da(4,4)+Da(1,3)*Da(3,4)*Da(4,1)+Da(1&
            &,4)*Da(3,1)*Da(4,3)- Da(1,1)*Da(3,4)*Da(4,3)-Da(1,3)&
            &*Da(3,1)*Da(4,4)-Da(1,4)*Da(3,3)*Da(4,1)
        Db(3,2)=Da(1,1)*Da(3,4)*Da(4,2)+Da(1,2)*Da(3,1)*Da(4,4)+Da(1&
            &,4)*Da(3,2)*Da(4,1)- Da(1,1)*Da(3,2)*Da(4,4)-Da(1,2)&
            &*Da(3,4)*Da(4,1)-Da(1,4)*Da(3,1)*Da(4,2)
        Db(4,2)=Da(1,1)*Da(3,2)*Da(4,3)+Da(1,2)*Da(3,3)*Da(4,1)+Da(1&
            &,3)*Da(3,1)*Da(4,2)- Da(1,1)*Da(3,3)*Da(4,2)-Da(1,2)&
            &*Da(3,1)*Da(4,3)-Da(1,3)*Da(3,2)*Da(4,1)
        Db(1,3)=Da(1,2)*Da(2,3)*Da(4,4)+Da(1,3)*Da(2,4)*Da(4,2)+Da(1&
            &,4)*Da(2,2)*Da(4,3)- Da(1,2)*Da(2,4)*Da(4,3)-Da(1,3)&
            &*Da(2,2)*Da(4,4)-Da(1,4)*Da(2,3)*Da(4,2)
        Db(2,3)=Da(1,1)*Da(2,4)*Da(4,3)+Da(1,3)*Da(2,1)*Da(4,4)+Da(1&
            &,4)*Da(2,3)*Da(4,1)- Da(1,1)*Da(2,3)*Da(4,4)-Da(1,3)&
            &*Da(2,4)*Da(4,1)-Da(1,4)*Da(2,1)*Da(4,3)
        Db(3,3)=Da(1,1)*Da(2,2)*Da(4,4)+Da(1,2)*Da(2,4)*Da(4,1)+Da(1&
            &,4)*Da(2,1)*Da(4,2)- Da(1,1)*Da(2,4)*Da(4,2)-Da(1,2)&
            &*Da(2,1)*Da(4,4)-Da(1,4)*Da(2,2)*Da(4,1)
        Db(4,3)=Da(1,1)*Da(2,3)*Da(4,2)+Da(1,2)*Da(2,1)*Da(4,3)+Da(1&
            &,3)*Da(2,2)*Da(4,1)- Da(1,1)*Da(2,2)*Da(4,3)-Da(1,2)&
            &*Da(2,3)*Da(4,1)-Da(1,3)*Da(2,1)*Da(4,2)
        Db(1,4)=Da(1,2)*Da(2,4)*Da(3,3)+Da(1,3)*Da(2,2)*Da(3,4)+Da(1&
            &,4)*Da(2,3)*Da(3,2)- Da(1,2)*Da(2,3)*Da(3,4)-Da(1,3)&
            &*Da(2,4)*Da(3,2)-Da(1,4)*Da(2,2)*Da(3,3)
        Db(2,4)=Da(1,1)*Da(2,3)*Da(3,4)+Da(1,3)*Da(2,4)*Da(3,1)+Da(1&
            &,4)*Da(2,1)*Da(3,3)- Da(1,1)*Da(2,4)*Da(3,3)-Da(1,3)&
            &*Da(2,1)*Da(3,4)-Da(1,4)*Da(2,3)*Da(3,1)
        Db(3,4)=Da(1,1)*Da(2,4)*Da(3,2)+Da(1,2)*Da(2,1)*Da(3,4)+Da(1&
            &,4)*Da(2,2)*Da(3,1)- Da(1,1)*Da(2,2)*Da(3,4)-Da(1,2)&
            &*Da(2,4)*Da(3,1)-Da(1,4)*Da(2,1)*Da(3,2)
        Db(4,4)=Da(1,1)*Da(2,2)*Da(3,3)+Da(1,2)*Da(2,3)*Da(3,1)+Da(1&
            &,3)*Da(2,1)*Da(3,2)- Da(1,1)*Da(2,3)*Da(3,2)-Da(1,2)&
            &*Da(2,1)*Da(3,3)-Da(1,3)*Da(2,2)*Da(3,1)
        daux=Da(1,1)*Da(2,2)*Da(3,3)*Da(4,4)+Da(1,1)*Da(2,3)*Da(3,4)&
            &*Da(4,2)+Da(1,1)*Da(2,4)*Da(3,2)*Da(4,3)+ Da(1,2)*Da(2&
            &,1)*Da(3,4)*Da(4,3)+Da(1,2)*Da(2,3)*Da(3,1)*Da(4,4)+Da(1&
            &,2)*Da(2,4)*Da(3,3)*Da(4,1)+ Da(1,3)*Da(2,1)*Da(3,2)&
            &*Da(4,4)+Da(1,3)*Da(2,2)*Da(3,4)*Da(4,1)+Da(1,3)*Da(2,4)&
            &*Da(3,1)*Da(4,2)+ Da(1,4)*Da(2,1)*Da(3,3)*Da(4,2)+Da(1&
            &,4)*Da(2,2)*Da(3,1)*Da(4,3)+Da(1,4)*Da(2,3)*Da(3,2)*Da(4&
            &,1)- Da(1,1)*Da(2,2)*Da(3,4)*Da(4,3)-Da(1,1)*Da(2,3)&
            &*Da(3,2)*Da(4,4)-Da(1,1)*Da(2,4)*Da(3,3)*Da(4,2)- Da(1&
            &,2)*Da(2,1)*Da(3,3)*Da(4,4)-Da(1,2)*Da(2,3)*Da(3,4)*Da(4&
            &,1)-Da(1,2)*Da(2,4)*Da(3,1)*Da(4,3)- Da(1,3)*Da(2,1)&
            &*Da(3,4)*Da(4,2)-Da(1,3)*Da(2,2)*Da(3,1)*Da(4,4)-Da(1,3)&
            &*Da(2,4)*Da(3,2)*Da(4,1)- Da(1,4)*Da(2,1)*Da(3,2)*Da(4&
            &,3)-Da(1,4)*Da(2,2)*Da(3,3)*Da(4,1)-Da(1,4)*Da(2,3)*Da(3&
            &,1)*Da(4,2)
        Dx=MATMUL(Db,Df)/daux

      CASE DEFAULT
        ! Use LAPACK routine for general NxN system, where N>4
        Dpiv=0; Dx=Df
        CALL DGESV(ndim,1,Da,ndim,Dpiv,Dx,ndim,info)
        
      END SELECT
    END SELECT
  END SUBROUTINE mprim_invertMatrixDble

  ! ***************************************************************************

!<subroutine>

  PURE SUBROUTINE mprim_invert4x4MatrixDirectDble(Da,Db)

!<description>
  ! This subroutine directly inverts a 4x4 system without any pivoting.
  ! 'Da' is a 2-dimensional 4x4 m matrix. The inverse of Da is written
  ! to the 2-dimensional 4x4 matrix Db.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Da and Db are assumed to be 4x4 arrays!
!</description>

!<input>
  ! source square matrix to be inverted
  REAL(DP), DIMENSION(4,4), INTENT(IN) :: Da
!</input>

!<output>
  ! destination square matrix; receives $A^{-1}$.
  REAL(DP), DIMENSION(4,4), INTENT(OUT) :: Db
!</output>

!</subroutine>

    REAL(DP) :: daux

      ! Explicit formula for 4x4 system
      Db(1,1)=Da(2,2)*Da(3,3)*Da(4,4)+Da(2,3)*Da(3,4)*Da(4,2)+Da(2&
          &,4)*Da(3,2)*Da(4,3)- Da(2,2)*Da(3,4)*Da(4,3)-Da(2,3)&
          &*Da(3,2)*Da(4,4)-Da(2,4)*Da(3,3)*Da(4,2)
      Db(2,1)=Da(2,1)*Da(3,4)*Da(4,3)+Da(2,3)*Da(3,1)*Da(4,4)+Da(2&
          &,4)*Da(3,3)*Da(4,1)- Da(2,1)*Da(3,3)*Da(4,4)-Da(2,3)&
          &*Da(3,4)*Da(4,1)-Da(2,4)*Da(3,1)*Da(4,3)
      Db(3,1)=Da(2,1)*Da(3,2)*Da(4,4)+Da(2,2)*Da(3,4)*Da(4,1)+Da(2&
          &,4)*Da(3,1)*Da(4,2)- Da(2,1)*Da(3,4)*Da(4,2)-Da(2,2)&
          &*Da(3,1)*Da(4,4)-Da(2,4)*Da(3,2)*Da(4,1)
      Db(4,1)=Da(2,1)*Da(3,3)*Da(4,2)+Da(2,2)*Da(3,1)*Da(4,3)+Da(2&
          &,3)*Da(3,2)*Da(4,1)- Da(2,1)*Da(3,2)*Da(4,3)-Da(2,2)&
          &*Da(3,3)*Da(4,1)-Da(2,3)*Da(3,1)*Da(4,2)
      Db(1,2)=Da(1,2)*Da(3,4)*Da(4,3)+Da(1,3)*Da(3,2)*Da(4,4)+Da(1&
          &,4)*Da(3,3)*Da(4,2)- Da(1,2)*Da(3,3)*Da(4,4)-Da(1,3)&
          &*Da(3,4)*Da(4,2)-Da(1,4)*Da(3,2)*Da(4,3)
      Db(2,2)=Da(1,1)*Da(3,3)*Da(4,4)+Da(1,3)*Da(3,4)*Da(4,1)+Da(1&
          &,4)*Da(3,1)*Da(4,3)- Da(1,1)*Da(3,4)*Da(4,3)-Da(1,3)&
          &*Da(3,1)*Da(4,4)-Da(1,4)*Da(3,3)*Da(4,1)
      Db(3,2)=Da(1,1)*Da(3,4)*Da(4,2)+Da(1,2)*Da(3,1)*Da(4,4)+Da(1&
          &,4)*Da(3,2)*Da(4,1)- Da(1,1)*Da(3,2)*Da(4,4)-Da(1,2)&
          &*Da(3,4)*Da(4,1)-Da(1,4)*Da(3,1)*Da(4,2)
      Db(4,2)=Da(1,1)*Da(3,2)*Da(4,3)+Da(1,2)*Da(3,3)*Da(4,1)+Da(1&
          &,3)*Da(3,1)*Da(4,2)- Da(1,1)*Da(3,3)*Da(4,2)-Da(1,2)&
          &*Da(3,1)*Da(4,3)-Da(1,3)*Da(3,2)*Da(4,1)
      Db(1,3)=Da(1,2)*Da(2,3)*Da(4,4)+Da(1,3)*Da(2,4)*Da(4,2)+Da(1&
          &,4)*Da(2,2)*Da(4,3)- Da(1,2)*Da(2,4)*Da(4,3)-Da(1,3)&
          &*Da(2,2)*Da(4,4)-Da(1,4)*Da(2,3)*Da(4,2)
      Db(2,3)=Da(1,1)*Da(2,4)*Da(4,3)+Da(1,3)*Da(2,1)*Da(4,4)+Da(1&
          &,4)*Da(2,3)*Da(4,1)- Da(1,1)*Da(2,3)*Da(4,4)-Da(1,3)&
          &*Da(2,4)*Da(4,1)-Da(1,4)*Da(2,1)*Da(4,3)
      Db(3,3)=Da(1,1)*Da(2,2)*Da(4,4)+Da(1,2)*Da(2,4)*Da(4,1)+Da(1&
          &,4)*Da(2,1)*Da(4,2)- Da(1,1)*Da(2,4)*Da(4,2)-Da(1,2)&
          &*Da(2,1)*Da(4,4)-Da(1,4)*Da(2,2)*Da(4,1)
      Db(4,3)=Da(1,1)*Da(2,3)*Da(4,2)+Da(1,2)*Da(2,1)*Da(4,3)+Da(1&
          &,3)*Da(2,2)*Da(4,1)- Da(1,1)*Da(2,2)*Da(4,3)-Da(1,2)&
          &*Da(2,3)*Da(4,1)-Da(1,3)*Da(2,1)*Da(4,2)
      Db(1,4)=Da(1,2)*Da(2,4)*Da(3,3)+Da(1,3)*Da(2,2)*Da(3,4)+Da(1&
          &,4)*Da(2,3)*Da(3,2)- Da(1,2)*Da(2,3)*Da(3,4)-Da(1,3)&
          &*Da(2,4)*Da(3,2)-Da(1,4)*Da(2,2)*Da(3,3)
      Db(2,4)=Da(1,1)*Da(2,3)*Da(3,4)+Da(1,3)*Da(2,4)*Da(3,1)+Da(1&
          &,4)*Da(2,1)*Da(3,3)- Da(1,1)*Da(2,4)*Da(3,3)-Da(1,3)&
          &*Da(2,1)*Da(3,4)-Da(1,4)*Da(2,3)*Da(3,1)
      Db(3,4)=Da(1,1)*Da(2,4)*Da(3,2)+Da(1,2)*Da(2,1)*Da(3,4)+Da(1&
          &,4)*Da(2,2)*Da(3,1)- Da(1,1)*Da(2,2)*Da(3,4)-Da(1,2)&
          &*Da(2,4)*Da(3,1)-Da(1,4)*Da(2,1)*Da(3,2)
      Db(4,4)=Da(1,1)*Da(2,2)*Da(3,3)+Da(1,2)*Da(2,3)*Da(3,1)+Da(1&
          &,3)*Da(2,1)*Da(3,2)- Da(1,1)*Da(2,3)*Da(3,2)-Da(1,2)&
          &*Da(2,1)*Da(3,3)-Da(1,3)*Da(2,2)*Da(3,1)
      daux=Da(1,1)*Da(2,2)*Da(3,3)*Da(4,4)+Da(1,1)*Da(2,3)*Da(3,4)&
          &*Da(4,2)+Da(1,1)*Da(2,4)*Da(3,2)*Da(4,3)+ Da(1,2)*Da(2&
          &,1)*Da(3,4)*Da(4,3)+Da(1,2)*Da(2,3)*Da(3,1)*Da(4,4)+Da(1&
          &,2)*Da(2,4)*Da(3,3)*Da(4,1)+ Da(1,3)*Da(2,1)*Da(3,2)&
          &*Da(4,4)+Da(1,3)*Da(2,2)*Da(3,4)*Da(4,1)+Da(1,3)*Da(2,4)&
          &*Da(3,1)*Da(4,2)+ Da(1,4)*Da(2,1)*Da(3,3)*Da(4,2)+Da(1&
          &,4)*Da(2,2)*Da(3,1)*Da(4,3)+Da(1,4)*Da(2,3)*Da(3,2)*Da(4&
          &,1)- Da(1,1)*Da(2,2)*Da(3,4)*Da(4,3)-Da(1,1)*Da(2,3)&
          &*Da(3,2)*Da(4,4)-Da(1,1)*Da(2,4)*Da(3,3)*Da(4,2)- Da(1&
          &,2)*Da(2,1)*Da(3,3)*Da(4,4)-Da(1,2)*Da(2,3)*Da(3,4)*Da(4&
          &,1)-Da(1,2)*Da(2,4)*Da(3,1)*Da(4,3)- Da(1,3)*Da(2,1)&
          &*Da(3,4)*Da(4,2)-Da(1,3)*Da(2,2)*Da(3,1)*Da(4,4)-Da(1,3)&
          &*Da(2,4)*Da(3,2)*Da(4,1)- Da(1,4)*Da(2,1)*Da(3,2)*Da(4&
          &,3)-Da(1,4)*Da(2,2)*Da(3,3)*Da(4,1)-Da(1,4)*Da(2,3)*Da(3&
          &,1)*Da(4,2)
      Db=Db*(1.0_DP/daux)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  PURE SUBROUTINE mprim_invertMatrixPivotDble(Da,ndim)

!<description>
  ! This subroutine directly inverts a (ndim x ndim) system with pivoting.
  ! 'Da' is a 2-dimensional (ndim x ndim) matrix and will be replaced
  ! by its inverse.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Da and Db are assumed to be 4x4 arrays!
!</description>

!<input>
  ! Dimension of the matrix Da and Db.
  INTEGER, INTENT(IN) :: ndim
!</input>

!<inputoutput>
  ! source square matrix to be inverted
  REAL(DP), DIMENSION(ndim,ndim), INTENT(INOUT) :: Da
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER, DIMENSION(ndim) :: Kindx,Kindy

    REAL(DP) :: dpivot,daux
    INTEGER :: idim1,idim2,ix,iy,indx,indy

    ! Perform factorization of matrix Da

    ! Initialization
    Kindx=0
    Kindy=0

    DO idim1=1,ndim

      ! Determine pivotal element
      dpivot=0

      DO iy=1,ndim
        IF (Kindy(iy) /= 0) CYCLE

        DO ix=1,ndim
          IF (Kindx(ix) /= 0) CYCLE

          IF (ABS(Da(ix,iy)) .LE. ABS(dpivot)) CYCLE
          dpivot=Da(ix,iy);  indx=ix;  indy=iy
        END DO
      END DO

      ! Return if pivotal element is zero
      IF (ABS(dpivot) .LE. 0._DP) RETURN

      Kindx(indx)=indy;  Kindy(indy)=indx;  Da(indx,indy)=1._DP&
          &/dpivot

      DO idim2=1,ndim
        IF (idim2 == indy) CYCLE 
        Da(1:indx-1,idim2)=Da(1:indx-1,idim2)-Da(1:indx-1,  &
            & indy)*Da(indx,idim2)/dpivot
        Da(indx+1:ndim,idim2)=Da(indx+1:ndim,idim2)-Da(indx+1:ndim&
            &,indy)*Da(indx,idim2)/dpivot
      END DO

      DO ix=1,ndim
        IF (ix /= indx) Da(ix,indy)=Da(ix,indy)/dpivot
      END DO

      DO iy=1,ndim
        IF (iy /= indy) Da(indx,iy)=-Da(indx,iy)/dpivot
      END DO
    END DO

    DO ix=1,ndim
      IF (Kindx(ix) == ix) CYCLE

      DO iy=1,ndim
        IF (Kindx(iy) == ix) EXIT
      END DO

      DO idim1=1,ndim
        daux=Da(ix,idim1)
        Da(ix,idim1)=Da(iy,idim1)
        Da(iy,idim1)=daux
      END DO

      Kindx(iy)=Kindx(ix);  Kindx(ix)=ix
    END DO

    DO ix=1,ndim
      IF (Kindy(ix) == ix) CYCLE

      DO iy=1,ndim
        IF (Kindy(iy) == ix) EXIT
      END DO

      DO idim1=1,ndim
        daux=Da(idim1,ix)
        Da(idim1,ix)=Da(idim1,iy)
        Da(idim1,iy)=daux
      END DO
      
      Kindy(iy)=Kindy(ix);  Kindy(ix)=ix
    END DO

  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mprim_invertMatrixSngl(Fa,Ff,Fx,ndim,ipar)

!<description>
    ! This subroutine performs the direct inversion of a NxN system.
    ! 
    ! If the parameter ipar=0, then only factorization of matrix A is performed. 
    ! For ipar=1, the vector x is calculated using the factorized matrix A. 
    ! For ipar=2, LAPACK routine DGESV is used to solve the dense linear system Ax=f. 
    ! In addition, for NDIM=2,3,4 explicit formulas are employed to replace the
    ! more expensive LAPACK routine.
!</description>

!<input>
    ! dimension of the matrix
    INTEGER, INTENT(IN) :: ndim

    ! source right-hand side vector
    REAL(SP), DIMENSION(ndim), INTENT(IN) :: Ff

    ! What to do?
    ! IPAR = 0 : invert matrix by means of Gaussian elimination with
    !            full pivoting and return the inverted matrix inv(A)
    ! IPAR = 1 : apply inverted matrix to the right-hand side vector
    !            and return x = inv(A)*f
    ! IPAR = 2 : invert matrix and apply it to the right-hand side
    !            vector. Return x = inv(A)*f
    INTEGER, INTENT(IN) :: ipar
!</input>

!<inputoutput> 
    ! source square matrix to be inverted
    REAL(SP), DIMENSION(ndim,ndim), INTENT(INOUT) :: Fa
!</inputoutput>

!<output>
    ! destination vector containing inverted matrix times right-hand
    ! side vector
    REAL(SP), DIMENSION(ndim), INTENT(OUT) :: Fx
!</output>
!</subroutine>

    ! local variables
    REAL(SP), DIMENSION(ndim,ndim) :: Fb
    REAL(SP), DIMENSION(ndim) :: Fpiv
    INTEGER, DIMENSION(ndim) :: Kindx,Kindy

    REAL(SP) :: fpivot,faux
    INTEGER :: idim1,idim2,ix,iy,indx,indy,info

    SELECT CASE (ipar)
    CASE (0)
      ! Perform factorization of matrix Da

      ! Initialization
      Kindx=0;  Kindy=0

      DO idim1=1,ndim

        ! Determine pivotal element
        fpivot=0

        DO iy=1,ndim
          IF (Kindy(iy) /= 0) CYCLE

          DO ix=1,ndim
            IF (Kindx(ix) /= 0) CYCLE

            IF (ABS(Fa(ix,iy)) .LE. ABS(fpivot)) CYCLE
            fpivot=Fa(ix,iy);  indx=ix;  indy=iy
          END DO
        END DO

        ! Return if pivotal element is zero
        IF (ABS(fpivot) .LE. 0._DP) RETURN

        Kindx(indx)=indy;  Kindy(indy)=indx;  Fa(indx,indy)=1._SP&
            &/fpivot

        DO idim2=1,ndim
          IF (idim2 == indy) CYCLE 
          Fa(1:indx-1,idim2)=Fa(1:indx-1,idim2)-Fa(1:indx-1,  &
              & indy)*Fa(indx,idim2)/fpivot
          Fa(indx+1:ndim,idim2)=Fa(indx+1:ndim,idim2)-Fa(indx+1:ndim&
              &,indy)*Fa(indx,idim2)/fpivot
        END DO

        DO ix=1,ndim
          IF (ix /= indx) Fa(ix,indy)=Fa(ix,indy)/fpivot
        END DO

        DO iy=1,ndim
          IF (iy /= indy) Fa(indx,iy)=-Fa(indx,iy)/fpivot
        END DO
      END DO

      DO ix=1,ndim
        IF (Kindx(ix) == ix) CYCLE

        DO iy=1,ndim
          IF (Kindx(iy) == ix) EXIT
        END DO

        DO idim1=1,ndim
          faux=Fa(ix,idim1)
          Fa(ix,idim1)=Fa(iy,idim1)
          Fa(iy,idim1)=faux
        END DO

        Kindx(iy)=Kindx(ix);  Kindx(ix)=ix
      END DO

      DO ix=1,ndim
        IF (Kindy(ix) == ix) CYCLE

        DO iy=1,ndim
          IF (Kindy(iy) == ix) EXIT
        END DO

        DO idim1=1,ndim
          faux=Fa(idim1,ix)
          Fa(idim1,ix)=Fa(idim1,iy)
          Fa(idim1,iy)=faux
        END DO
        
        Kindy(iy)=Kindy(ix);  Kindy(ix)=ix
      END DO


    CASE (1)    
      ! Perform inversion of Da to solve the system Da * Dx = Df
      DO idim1=1,ndim
        Fx(idim1)=0
        DO idim2=1,ndim
          Fx(idim1)=Fx(idim1)+Fa(idim1,idim2)*Ff(idim2)
        END DO
      END DO

    CASE (2)
      ! Solve the dense linear system Ax=f calling LAPACK routine

      SELECT CASE(ndim)
      CASE (2)
        ! Explicit formula for 2x2 system
        Fb(1,1)= Fa(2,2)
        Fb(2,1)=-Fa(2,1)
        Fb(1,2)=-Fa(1,2)
        Fb(2,2)= Fa(1,1)
        faux=Fa(1,1)*Fa(2,2)-Fa(1,2)*Fa(2,1)
        Fx=MATMUL(Fb,Ff)/faux

      CASE (3)
        ! Explicit formula for 3x3 system
        Fb(1,1)=Fa(2,2)*Fa(3,3)-Fa(2,3)*Fa(3,2)
        Fb(2,1)=Fa(2,3)*Fa(3,1)-Fa(2,1)*Fa(3,3)
        Fb(3,1)=Fa(2,1)*Fa(3,2)-Fa(2,2)*Fa(3,1)
        Fb(1,2)=Fa(1,3)*Fa(3,2)-Fa(1,2)*Fa(3,3)
        Fb(2,2)=Fa(1,1)*Fa(3,3)-Fa(1,3)*Fa(3,1)
        Fb(3,2)=Fa(1,2)*Fa(3,1)-Fa(1,1)*Fa(3,2)
        Fb(1,3)=Fa(1,2)*Fa(2,3)-Fa(1,3)*Fa(2,2)
        Fb(2,3)=Fa(1,3)*Fa(2,1)-Fa(1,1)*Fa(2,3)
        Fb(3,3)=Fa(1,1)*Fa(2,2)-Fa(1,2)*Fa(2,1)
        faux=Fa(1,1)*Fa(2,2)*Fa(3,3)+Fa(2,1)*Fa(3,2)*Fa(1,3)+ Fa(3,1)&
            &*Fa(1,2)*Fa(2,3)-Fa(1,1)*Fa(3,2)*Fa(2,3)- Fa(3,1)*Fa(2&
            &,2)*Fa(1,3)-Fa(2,1)*Fa(1,2)*Fa(3,3)
        Fx=MATMUL(Fb,Ff)/faux

      CASE (4)
        ! Explicit formula for 4x4 system
        Fb(1,1)=Fa(2,2)*Fa(3,3)*Fa(4,4)+Fa(2,3)*Fa(3,4)*Fa(4,2)+Fa(2&
            &,4)*Fa(3,2)*Fa(4,3)- Fa(2,2)*Fa(3,4)*Fa(4,3)-Fa(2,3)&
            &*Fa(3,2)*Fa(4,4)-Fa(2,4)*Fa(3,3)*Fa(4,2)
        Fb(2,1)=Fa(2,1)*Fa(3,4)*Fa(4,3)+Fa(2,3)*Fa(3,1)*Fa(4,4)+Fa(2&
            &,4)*Fa(3,3)*Fa(4,1)- Fa(2,1)*Fa(3,3)*Fa(4,4)-Fa(2,3)&
            &*Fa(3,4)*Fa(4,1)-Fa(2,4)*Fa(3,1)*Fa(4,3)
        Fb(3,1)=Fa(2,1)*Fa(3,2)*Fa(4,4)+Fa(2,2)*Fa(3,4)*Fa(4,1)+Fa(2&
            &,4)*Fa(3,1)*Fa(4,2)- Fa(2,1)*Fa(3,4)*Fa(4,2)-Fa(2,2)&
            &*Fa(3,1)*Fa(4,4)-Fa(2,4)*Fa(3,2)*Fa(4,1)
        Fb(4,1)=Fa(2,1)*Fa(3,3)*Fa(4,2)+Fa(2,2)*Fa(3,1)*Fa(4,3)+Fa(2&
            &,3)*Fa(3,2)*Fa(4,1)- Fa(2,1)*Fa(3,2)*Fa(4,3)-Fa(2,2)&
            &*Fa(3,3)*Fa(4,1)-Fa(2,3)*Fa(3,1)*Fa(4,2)
        Fb(1,2)=Fa(1,2)*Fa(3,4)*Fa(4,3)+Fa(1,3)*Fa(3,2)*Fa(4,4)+Fa(1&
            &,4)*Fa(3,3)*Fa(4,2)- Fa(1,2)*Fa(3,3)*Fa(4,4)-Fa(1,3)&
            &*Fa(3,4)*Fa(4,2)-Fa(1,4)*Fa(3,2)*Fa(4,3)
        Fb(2,2)=Fa(1,1)*Fa(3,3)*Fa(4,4)+Fa(1,3)*Fa(3,4)*Fa(4,1)+Fa(1&
            &,4)*Fa(3,1)*Fa(4,3)- Fa(1,1)*Fa(3,4)*Fa(4,3)-Fa(1,3)&
            &*Fa(3,1)*Fa(4,4)-Fa(1,4)*Fa(3,3)*Fa(4,1)
        Fb(3,2)=Fa(1,1)*Fa(3,4)*Fa(4,2)+Fa(1,2)*Fa(3,1)*Fa(4,4)+Fa(1&
            &,4)*Fa(3,2)*Fa(4,1)- Fa(1,1)*Fa(3,2)*Fa(4,4)-Fa(1,2)&
            &*Fa(3,4)*Fa(4,1)-Fa(1,4)*Fa(3,1)*Fa(4,2)
        Fb(4,2)=Fa(1,1)*Fa(3,2)*Fa(4,3)+Fa(1,2)*Fa(3,3)*Fa(4,1)+Fa(1&
            &,3)*Fa(3,1)*Fa(4,2)- Fa(1,1)*Fa(3,3)*Fa(4,2)-Fa(1,2)&
            &*Fa(3,1)*Fa(4,3)-Fa(1,3)*Fa(3,2)*Fa(4,1)
        Fb(1,3)=Fa(1,2)*Fa(2,3)*Fa(4,4)+Fa(1,3)*Fa(2,4)*Fa(4,2)+Fa(1&
            &,4)*Fa(2,2)*Fa(4,3)- Fa(1,2)*Fa(2,4)*Fa(4,3)-Fa(1,3)&
            &*Fa(2,2)*Fa(4,4)-Fa(1,4)*Fa(2,3)*Fa(4,2)
        Fb(2,3)=Fa(1,1)*Fa(2,4)*Fa(4,3)+Fa(1,3)*Fa(2,1)*Fa(4,4)+Fa(1&
            &,4)*Fa(2,3)*Fa(4,1)- Fa(1,1)*Fa(2,3)*Fa(4,4)-Fa(1,3)&
            &*Fa(2,4)*Fa(4,1)-Fa(1,4)*Fa(2,1)*Fa(4,3)
        Fb(3,3)=Fa(1,1)*Fa(2,2)*Fa(4,4)+Fa(1,2)*Fa(2,4)*Fa(4,1)+Fa(1&
            &,4)*Fa(2,1)*Fa(4,2)- Fa(1,1)*Fa(2,4)*Fa(4,2)-Fa(1,2)&
            &*Fa(2,1)*Fa(4,4)-Fa(1,4)*Fa(2,2)*Fa(4,1)
        Fb(4,3)=Fa(1,1)*Fa(2,3)*Fa(4,2)+Fa(1,2)*Fa(2,1)*Fa(4,3)+Fa(1&
            &,3)*Fa(2,2)*Fa(4,1)- Fa(1,1)*Fa(2,2)*Fa(4,3)-Fa(1,2)&
            &*Fa(2,3)*Fa(4,1)-Fa(1,3)*Fa(2,1)*Fa(4,2)
        Fb(1,4)=Fa(1,2)*Fa(2,4)*Fa(3,3)+Fa(1,3)*Fa(2,2)*Fa(3,4)+Fa(1&
            &,4)*Fa(2,3)*Fa(3,2)- Fa(1,2)*Fa(2,3)*Fa(3,4)-Fa(1,3)&
            &*Fa(2,4)*Fa(3,2)-Fa(1,4)*Fa(2,2)*Fa(3,3)
        Fb(2,4)=Fa(1,1)*Fa(2,3)*Fa(3,4)+Fa(1,3)*Fa(2,4)*Fa(3,1)+Fa(1&
            &,4)*Fa(2,1)*Fa(3,3)- Fa(1,1)*Fa(2,4)*Fa(3,3)-Fa(1,3)&
            &*Fa(2,1)*Fa(3,4)-Fa(1,4)*Fa(2,3)*Fa(3,1)
        Fb(3,4)=Fa(1,1)*Fa(2,4)*Fa(3,2)+Fa(1,2)*Fa(2,1)*Fa(3,4)+Fa(1&
            &,4)*Fa(2,2)*Fa(3,1)- Fa(1,1)*Fa(2,2)*Fa(3,4)-Fa(1,2)&
            &*Fa(2,4)*Fa(3,1)-Fa(1,4)*Fa(2,1)*Fa(3,2)
        Fb(4,4)=Fa(1,1)*Fa(2,2)*Fa(3,3)+Fa(1,2)*Fa(2,3)*Fa(3,1)+Fa(1&
            &,3)*Fa(2,1)*Fa(3,2)- Fa(1,1)*Fa(2,3)*Fa(3,2)-Fa(1,2)&
            &*Fa(2,1)*Fa(3,3)-Fa(1,3)*Fa(2,2)*Fa(3,1)
        faux=Fa(1,1)*Fa(2,2)*Fa(3,3)*Fa(4,4)+Fa(1,1)*Fa(2,3)*Fa(3,4)&
            &*Fa(4,2)+Fa(1,1)*Fa(2,4)*Fa(3,2)*Fa(4,3)+ Fa(1,2)*Fa(2&
            &,1)*Fa(3,4)*Fa(4,3)+Fa(1,2)*Fa(2,3)*Fa(3,1)*Fa(4,4)+Fa(1&
            &,2)*Fa(2,4)*Fa(3,3)*Fa(4,1)+ Fa(1,3)*Fa(2,1)*Fa(3,2)&
            &*Fa(4,4)+Fa(1,3)*Fa(2,2)*Fa(3,4)*Fa(4,1)+Fa(1,3)*Fa(2,4)&
            &*Fa(3,1)*Fa(4,2)+ Fa(1,4)*Fa(2,1)*Fa(3,3)*Fa(4,2)+Fa(1&
            &,4)*Fa(2,2)*Fa(3,1)*Fa(4,3)+Fa(1,4)*Fa(2,3)*Fa(3,2)*Fa(4&
            &,1)- Fa(1,1)*Fa(2,2)*Fa(3,4)*Fa(4,3)-Fa(1,1)*Fa(2,3)&
            &*Fa(3,2)*Fa(4,4)-Fa(1,1)*Fa(2,4)*Fa(3,3)*Fa(4,2)- Fa(1&
            &,2)*Fa(2,1)*Fa(3,3)*Fa(4,4)-Fa(1,2)*Fa(2,3)*Fa(3,4)*Fa(4&
            &,1)-Fa(1,2)*Fa(2,4)*Fa(3,1)*Fa(4,3)- Fa(1,3)*Fa(2,1)&
            &*Fa(3,4)*Fa(4,2)-Fa(1,3)*Fa(2,2)*Fa(3,1)*Fa(4,4)-Fa(1,3)&
            &*Fa(2,4)*Fa(3,2)*Fa(4,1)- Fa(1,4)*Fa(2,1)*Fa(3,2)*Fa(4&
            &,3)-Fa(1,4)*Fa(2,2)*Fa(3,3)*Fa(4,1)-Fa(1,4)*Fa(2,3)*Fa(3&
            &,1)*Fa(4,2)
        Fx=MATMUL(Fb,Ff)/faux

      CASE DEFAULT
        ! Use LAPACK routine for general NxN system, where N>4
        Fpiv=0; Fx=Ff
        CALL SGESV(ndim,1,Fa,ndim,Fpiv,Fx,ndim,info)
        
      END SELECT
    END SELECT
  END SUBROUTINE mprim_invertMatrixSngl
  
  ! ***************************************************************************

!<function>

  ELEMENTAL FUNCTION kronecker(i,j) RESULT(kron)
    
!<description>
    ! Compute the Kronecker delta symbol 
    ! $\delta_{ij}\left\{\begin{array}{ll}
    !  1 & i=j\\
    !  0 & i\ne j
    !  \end{array}\right.$
!</description>

!<input>
    ! Evaluation points I and J
    INTEGER, INTENT(IN) :: i,j
!</input>

!<result>
    ! Kronecker delta symbol
    INTEGER :: kron
!</result>
!</function>

    kron=MERGE(1,0,i==j)
  END FUNCTION kronecker

  !************************************************************************

!<subroutine>
  
  ELEMENTAL SUBROUTINE mprim_linearRescale(dx,da,db,dc,dd,dy)
  
!<description>
  ! Scales a coordinate x linearly from the interval [a,b] to the
  ! interval [c,d].
!</description>

!<input>
  ! coordinate to be rescaled
  REAL(DP), INTENT(IN) :: dx
  
  ! [a,b] - source interval
  REAL(DP), INTENT(IN) :: da,db
  
  ! [c,d] - destination interval
  REAL(DP), INTENT(IN) :: dc,dd
!</input>

!<output>
  ! Rescaled coordinate
  REAL(DP), INTENT(OUT) :: dy
!</output>

!</subroutine>

    REAL(DP) :: d1,d2,d3

    ! Calculate the coefficients of the transformation 
    !    D1*A+D2 = C, D1*B+D2 = D.
    ! Use them to calculate Y=D1*X+D2.

    IF (dA .EQ. db) THEN
      dy = dc
      RETURN
    END IF

    d3 = 1.0_DP/(da-db)
    d1 = (dc-dd)*d3
    d2 = (-db*dc+da*dd)*d3
    
    dy = d1*dx+d2

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  ELEMENTAL SUBROUTINE mprim_quadraticInterpolation (dx,d1,d2,d3,dy)
  
!<description>
  ! Calculates a quadratic interpolation. dx is a value in the range $[-1,1]$.
  ! The routine calculates the value $dy:=p(dx)$ with $p(.)$ being the quadratic
  ! interpolation polynomial with $p(-1)=d1$, $p(0)=d2$ and $p(1)=d3$.
!<description>
  
!<input>
  ! The parameter value in the range $[-1,1]$ where the polynomial should be evaluated.
  REAL(DP), INTENT(IN) :: dx
  
  ! The value $p(-1)$.
  REAL(DP), INTENT(IN) :: d1

  ! The value $p(0)$.
  REAL(DP), INTENT(IN) :: d2

  ! The value $p(1)$.
  REAL(DP), INTENT(IN) :: d3
!</input>

!<output>
  ! The value $p(dx)$.
  REAL(DP), INTENT(OUT) :: dy
!</output>
  
!</subroutine>

    ! The polynomial p(t) = a + bt + ct^2 has to fulfill:
    !   p(-1)=d1, p(0)=d2, p(1)=d3.
    !
    ! So the polynomial has the shape:
    !   p(t) = d2  +  1/2(d3-d1)t  +  1/2(d3-2d2+d1)t^2
    !
    ! The Horner scheme gives us:
    !   p(t) = 1/2 ( (d3-2d2+d1)t + (d3-d1) ) t + d2
    
    dy = 0.5_DP * ( (d3 - 2.0_DP*d2 + d1)*dx + (d3-d1) ) * dx + d2

  END SUBROUTINE

  !************************************************************************

  SUBROUTINE mprim_SVD_factorise(Da,mdim,ndim,Dd,Db,btransposedOpt)

!<description>
    ! This subroutine computes the factorisation for a singular value
    ! decomposition. Given an ndim-by-mdim matrix Da, the routine
    ! decomposes it into the product $$ A = U * D * B^T$$
    ! where $U$ overwrites the matrix A and is returned in its memory
    ! position, $D$ is a diagonal matrix returned as vector Dd and
    ! $B$ is an n-by-n square matrix which is returned (instead
    ! of its transpose) as matrix Db.
    !
    ! The optional parameter btransposedOpt can be used to indicate 
    ! that matrix A is stored in transposed format. Note that 
    ! matrices D and B are not affected by this fact.
    !
    ! The details of this algorithm are given in numerical recipes in F90
!</description>

!<input>
    ! Dimensions of the rectangular matrix
    INTEGER, INTENT(IN)                           :: ndim,mdim

    ! OPTIONAL: Flag to indicate if the matrix A is transposed
    LOGICAL, INTENT(IN), OPTIONAL                 :: btransposedOpt
!</input>
    
!<inputoutput>
    ! Matrix that should be factorized on input.
    ! Matrix U on output.
    REAL(DP), DIMENSION(mdim,ndim), INTENT(INOUT) :: Da
!</inputoutput>

!<output>
    ! Diagonal matrix
    REAL(DP), DIMENSION(:), INTENT(OUT)           :: Dd

    ! Square matrix
    REAL(DP), DIMENSION(:,:), INTENT(OUT)         :: Db
!</output>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(SIZE(Dd)) :: rv1
    REAL(DP) :: g, scale, anorm, s, f, h, c, x, y, z, dot
    INTEGER  :: its, i, j, jj, k, l, nm, n, m, idot
    INTEGER, PARAMETER :: MAX_ITS = 30
    LOGICAL :: btransposed

    ! Do we have an optional transposed flag?
    IF (PRESENT(btransposedOpt)) THEN
      btransposed = btransposedOpt
    ELSE
      btransposed = .FALSE.
    END IF

    ! Is the matrix transposed?
    IF (btransposed) THEN
      n=mdim; m=ndim
    ELSE
      n=ndim; m=mdim
    END IF

    ! Check if number of equations is larger than number of unknowns
    IF (m < n) THEN
      CALL output_line('Fewer equations than unknowns!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mprim_SVD_factorise')
      CALL sys_halt()
    END IF
    
    
    IF (btransposed) THEN

      ! Householder reduction to bidiagonal form
      g     = 0.0_DP
      scale = 0.0_DP
      anorm = 0.0_DP
      
      DO i = 1, n
        
        l      = i+1
        rv1(i) = scale*g
        g      = 0.0_DP
        scale  = 0.0_DP
        
        IF (i .LE. m) THEN
          scale = SUM(ABS(Da(i,i:m)))
          IF (scale .GT. SYS_EPSREAL) THEN
            Da(i,i:m) = Da(i,i:m)/scale
            s = 0.0_DP
            DO idot = i, m
              s = s+Da(i,idot)*Da(i,idot)
            END DO
            f = Da(i,i)
            g = -SIGN(SQRT(s),f)
            h = f*g-s
            Da(i,i) = f-g
            IF (i .NE. n) THEN
              DO j = l, n
                s = 0.0_DP
                DO k = i, m
                  s = s+Da(i,k)*Da(j,k)
                END DO
                f = s/h
                DO k = i, m
                  Da(j,k) = Da(j,k)+f*Da(i,k)
                END DO
              END DO
            END IF
            DO k = i, m
              Da(i,k) = scale*Da(i,k)
            END DO
          END IF
        END IF
        
        Dd(i) = scale*g
        g     = 0.0_DP
        s     = 0.0_DP
        scale = 0.0_DP

        IF ((i .LE. m) .AND. (i .NE. n)) THEN
          scale = SUM(ABS(Da(l:n,i)))
          IF (scale .NE. 0.0_DP) THEN
            DO k = l, n
              Da(k,i) = Da(k,i)/scale
              s = s+Da(k,i)*Da(k,i)
            END DO
            f = Da(l,i)
            g = -SIGN(SQRT(s),f)
            h = f*g-s
            Da(l,i) = f-g
            DO k = l, n
              rv1(k) = Da(k,i)/h
            END DO
            
            IF (i .NE. m) THEN
              DO j = l, m
                s = 0.0_DP
                DO k = l, n
                  s = s+Da(k,j)*Da(k,i)
                END DO
                DO k = l, n
                  Da(k,j) = Da(k,j)+s*rv1(k)
                END DO
              END DO
            END IF
            
            DO k = l, n
              Da(k,i) = scale*Da(k,i)
            END DO
          END IF
        END IF
        anorm = MAX(anorm,(ABS(Dd(i))+ABS(rv1(i))))
      END DO
      
      ! Accumulation of right-hand transformations
      DO i = n, 1, -1
        IF (i .LT. n) THEN
          IF (g .NE. 0.0_DP) THEN
            DO j = l, n
              Db(j,i) = (Da(j,i)/Da(l,i))/g
            END DO
            DO j = l, n
              s = 0.0_DP
              DO k = l, n
                s = s+Da(k,i)*Db(k,j)
              END DO
              DO k = l, n
                Db(k,j) = Db(k,j)+s*Db(k,i)
              END DO
            END DO
          END IF
          DO j = l, n
            Db(i,j) = 0.0_DP
            Db(j,i) = 0.0_DP
          END DO
        END IF
        Db(i,i) = 1.0_DP
        g = rv1(i)
        l = i
      END DO

      ! Accumulation of left-hand transformations
      DO i = n, 1, -1
        l = i+1
        g = Dd(i)
        IF (i .LT. n) THEN
          DO j = l, n
            Da(j,i) = 0.0_DP
          END DO
        END IF
        IF (g .NE. 0.0_DP) THEN
          g = 1.0_DP/g
          IF (i .NE. n) THEN
            DO j = l, n
              s = 0.0_DP
              DO k = l, m
                s = s+Da(i,k)*Da(j,k)
              END DO
              f = (s/Da(i,i))*g
              DO k = i, m
                Da(j,k) = Da(j,k)+f*Da(i,k)
              END DO
            END DO
          END IF
          DO j = i, m
            Da(i,j) = Da(i,j)*g
          END DO
        ELSE
          DO j = i, m
            Da(i,j) = 0.0_DP
          END DO
        END IF
        Da(i,i) = Da(i,i)+1.0_DP
      END DO

      ! Diagonalization of the bidiagonal form
      DO k = n, 1, -1
        DO its = 1, MAX_ITS
          DO l = k, 1, -1
            nm = l-1
            IF ((ABS(rv1(l))+anorm) .EQ. anorm) GOTO 5
            IF ((ABS(Dd(nm))+anorm) .EQ. anorm) GOTO 4
          END DO
4         CONTINUE
          c = 0.0_DP
          s = 1.0_DP
          DO i = l, k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            IF ((ABS(f)+anorm) .EQ. anorm) GOTO 5
            g = Dd(i)
            h = SQRT(f*f+g*g)
            Dd(i) = h
            h = 1.0_DP/h
            c = g*h
            s = -(f*h)
            DO j = 1, m
              y = Da(nm,j)
              z = Da(i,j)
              Da(nm,j) = (y*c)+(z*s)
              Da(i,j) = -(y*s)+(z*c)
            END DO
          END DO
5         CONTINUE
          z = Dd(k)
          IF (l .EQ. k) THEN
            IF (z .LT. 0.0_DP) THEN
              Dd(k) = -z
              DO j = 1, n
                Db(j,k) = -Db(j,k)
              END DO
            END IF
            GOTO 6
          END IF
          IF (its .EQ. MAX_ITS) THEN
            CALL output_line('Convergence failed!',&
                OU_CLASS_ERROR,OU_MODE_STD,'mprim_SVD_factorise')
            CALL sys_halt()
          END IF
          x = Dd(l)
          nm = k-1
          y = Dd(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_DP*h*y)
          g = SQRT(f*f+1.0_DP)
          f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x
          
          ! Next QR transformation
          c = 1.0_DP
          s = 1.0_DP
          DO j = l, nm
            i = j+1
            g = rv1(i)
            y = Dd(i)
            h = s*g
            g = c*g
            z = SQRT(f*f+h*h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c)+(g*s)
            g = -(x*s)+(g*c)
            h = y*s
            y = y*c
            DO jj = 1, n
              x = Db(jj,j)
              z = Db(jj,i)
              Db(jj,j) = (x*c)+(z*s)
              Db(jj,i) = -(x*s)+(z*c)
            END DO
            z = SQRT(f*f+h*h)
            Dd(j) = z
            IF (z .NE. 0.0_DP) THEN
              z = 1.0_DP/z
              c = f*z
              s = h*z
            END IF
            f = (c*g)+(s*y)
            x = -(s*g)+(c*y)
            DO jj = 1, m
              y = Da(j,jj)
              z = Da(i,jj)
              Da(j,jj) = (y*c)+(z*s)
              Da(i,jj) = -(y*s)+(z*c)
            END DO
          END DO
          rv1(l) = 0.0_DP
          rv1(k) = f
          Dd(k) = x
        END DO
6       CONTINUE
      END DO

    ELSE

      ! Householder reduction to bidiagonal form
      g     = 0.0_DP
      scale = 0.0_DP
      anorm = 0.0_DP
      
      DO i = 1, n
        
        l      = i+1
        rv1(i) = scale*g
        g      = 0.0_DP
        scale  = 0.0_DP
        
        IF (i .LE. m) THEN
          scale = SUM(ABS(Da(i:m,i)))
          IF (scale .GT. SYS_EPSREAL) THEN
            Da(i:m,i) = Da(i:m,i)/scale
            s = 0.0_DP
            DO idot = i, m
              s = s+Da(idot,i)*Da(idot,i)
            END DO
            f = Da(i,i)
            g = -SIGN(SQRT(s),f)
            h = f*g-s
            Da(i,i) = f-g
            IF (i .NE. n) THEN
              DO j = l, n
                s = 0.0_DP
                DO k = i, m
                  s = s+Da(k,i)*Da(k,j)
                END DO
                f = s/h
                DO k = i, m
                  Da(k,j) = Da(k,j)+f*Da(k,i)
                END DO
              END DO
            END IF
            DO k = i, m
              Da(k,i) = scale*Da(k,i)
            END DO
          END IF
        END IF
        
        Dd(i) = scale*g
        g     = 0.0_DP
        s     = 0.0_DP
        scale = 0.0_DP
        
        IF ((i .LE. m) .AND. (i .NE. n)) THEN
          scale = SUM(ABS(Da(i,l:n)))
          IF (scale .NE. 0.0_DP) THEN
            DO k = l, n
              Da(i,k) = Da(i,k)/scale
              s = s+Da(i,k)*Da(i,k)
            END DO
            f = Da(i,l)
            g = -SIGN(SQRT(s),f)
            h = f*g-s
            Da(i,l) = f-g
            DO k = l, n
              rv1(k) = Da(i,k)/h
            END DO
            
            IF (i .NE. m) THEN
              DO j = l, m
                s = 0.0_DP
                DO k = l, n
                  s = s+Da(j,k)*Da(i,k)
                END DO
                DO k = l, n
                  Da(j,k) = Da(j,k)+s*rv1(k)
                END DO
              END DO
            END IF
            
            DO k = l, n
              Da(i,k) = scale*Da(i,k)
            END DO
          END IF
        END IF
        anorm = MAX(anorm,(ABS(Dd(i))+ABS(rv1(i))))
      END DO

      ! Accumulation of right-hand transformations
      DO i = n, 1, -1
        IF (i .LT. n) THEN
          IF (g .NE. 0.0_DP) THEN
            DO j = l, n
              Db(j,i) = (Da(i,j)/Da(i,l))/g
            END DO
            DO j = l, n
              s = 0.0_DP
              DO k = l, n
                s = s+Da(i,k)*Db(k,j)
              END DO
              DO k = l, n
                Db(k,j) = Db(k,j)+s*Db(k,i)
              END DO
            END DO
          END IF
          DO j = l, n
            Db(i,j) = 0.0_DP
            Db(j,i) = 0.0_DP
          END DO
        END IF
        Db(i,i) = 1.0_DP
        g = rv1(i)
        l = i
      END DO

      ! Accumulation of left-hand transformations
      DO i = n, 1, -1
        l = i+1
        g = Dd(i)
        IF (i .LT. n) THEN
          DO j = l, n
            Da(i,j) = 0.0_DP
          END DO
        END IF
        IF (g .NE. 0.0_DP) THEN
          g = 1.0_DP/g
          IF (i .NE. n) THEN
            DO j = l, n
              s = 0.0_DP
              DO k = l, m
                s = s+Da(k,i)*Da(k,j)
              END DO
              f = (s/Da(i,i))*g
              DO k = i, m
                Da(k,j) = Da(k,j)+f*Da(k,i)
              END DO
            END DO
          END IF
          DO j = i, m
            Da(j,i) = Da(j,i)*g
          END DO
        ELSE
          DO j = i, m
            Da(j,i) = 0.0_DP
          END DO
        END IF
        Da(i,i) = Da(i,i)+1.0_DP
      END DO

      ! Diagonalization of the bidiagonal form
      DO k = n, 1, -1
        DO its = 1, MAX_ITS
          DO l = k, 1, -1
            nm = l-1
            IF ((ABS(rv1(l))+anorm) .EQ. anorm) GOTO 2
            IF ((ABS(Dd(nm))+anorm) .EQ. anorm) GOTO 1
          END DO
1         CONTINUE
          c = 0.0_DP
          s = 1.0_DP
          DO i = l, k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            IF ((ABS(f)+anorm) .EQ. anorm) GOTO 2
            g = Dd(i)
            h = SQRT(f*f+g*g)
            Dd(i) = h
            h = 1.0_DP/h
            c = g*h
            s = -(f*h)
            DO j = 1, m
              y = Da(j,nm)
              z = Da(j,i)
              Da(j,nm) = (y*c)+(z*s)
              Da(j,i) = -(y*s)+(z*c)
            END DO
          END DO
2         CONTINUE
          z = Dd(k)
          IF (l .EQ. k) THEN
            IF (z .LT. 0.0_DP) THEN
              Dd(k) = -z
              DO j = 1, n
                Db(j,k) = -Db(j,k)
              END DO
            END IF
            GOTO 3
          END IF
          IF (its .EQ. MAX_ITS) THEN
            CALL output_line('Convergence failed!',&
                OU_CLASS_ERROR,OU_MODE_STD,'mprim_SVD_factorise')
            CALL sys_halt()
          END IF
          x = Dd(l)
          nm = k-1
          y = Dd(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_DP*h*y)
          g = SQRT(f*f+1.0_DP)
          f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x
          
          ! Next QR transformation
          c = 1.0_DP
          s = 1.0_DP
          DO j = l, nm
            i = j+1
            g = rv1(i)
            y = Dd(i)
            h = s*g
            g = c*g
            z = SQRT(f*f+h*h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c)+(g*s)
            g = -(x*s)+(g*c)
            h = y*s
            y = y*c
            DO jj = 1, n
              x = Db(jj,j)
              z = Db(jj,i)
              Db(jj,j) = (x*c)+(z*s)
              Db(jj,i) = -(x*s)+(z*c)
            END DO
            z = SQRT(f*f+h*h)
            Dd(j) = z
            IF (z .NE. 0.0_DP) THEN
              z = 1.0_DP/z
              c = f*z
              s = h*z
            END IF
            f = (c*g)+(s*y)
            x = -(s*g)+(c*y)
            DO jj = 1, m
              y = Da(jj,j)
              z = Da(jj,i)
              Da(jj,j) = (y*c)+(z*s)
              Da(jj,i) = -(y*s)+(z*c)
            END DO
          END DO
          rv1(l) = 0.0_DP
          rv1(k) = f
          Dd(k) = x
        END DO
3       CONTINUE
      END DO
      
    END IF
  END SUBROUTINE mprim_SVD_factorise

  !************************************************************************

  SUBROUTINE mprim_SVD_backsubst(Da,mdim,ndim,Dd,Db,Dx,Df,btransposedOpt)

!<description>
    ! This subroutine solves $A * x = f$ for vector $x$, where the rectangular
    ! matrix $A$ has been decomposed into $Da$, $Dd$ and $Db$ by the routine
    ! mprim_SVD_factorise. The optional parameter btransposedOpt can be used
    ! to indicate that matrix A is stored in transposed format.
!</description>

!<input>
    ! Dimensions of the rectangular matrix
    INTEGER, INTENT(IN)                           :: ndim,mdim

    ! OPTIONAL: Flag to indicate if the matrix A is transposed
    LOGICAL, INTENT(IN), OPTIONAL                 :: btransposedOpt

    ! Factorised matrux U
    REAL(DP), DIMENSION(mdim,ndim), INTENT(IN)    :: Da

    ! Diagonal matrix D
    REAL(DP), DIMENSION(:), INTENT(IN)            :: Dd

    ! Square matrix B
    REAL(DP), DIMENSION(:,:), INTENT(IN)          :: Db

    ! Right-hand side vector
    REAL(DP), DIMENSION(*), INTENT(IN)            :: Df
!</input>

!<output>
    ! Solution vector
    REAL(DP), DIMENSION(:), INTENT(OUT)           :: Dx
!</output>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(SIZE(Dd)) :: Daux
    INTEGER :: i,n,m
    LOGICAL :: btransposed

    ! Do we have an optional transposed flag?
    IF (PRESENT(btransposedOpt)) THEN
      btransposed = btransposedOpt
    ELSE
      btransposed = .FALSE.
    END IF

    ! Is the matrix transposed?
    IF (btransposed) THEN
      n=mdim; m=ndim
    ELSE
      n=ndim; m=mdim
    END IF
    
    ! Check if number of equations is larger than number of unknowns
    IF (m < n) THEN
      CALL output_line('Fewer equations than unknowns!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mprim_SVD_backsubst')
      CALL sys_halt()
    END IF
    
    ! Compute aux = (U^T * f)/D where D_i /= 0
    IF (btransposed) THEN
      CALL DGEMV('n', mdim, ndim, 1.0_DP, Da, mdim, Df, 1, 0.0_DP, Daux, 1)
    ELSE
      CALL DGEMV('t', mdim, ndim, 1.0_DP, Da, mdim, Df, 1, 0.0_DP, Daux, 1)
    END IF

    DO i =1, SIZE(Dd)
      IF (ABS(Dd(i)) .GT. SYS_EPSREAL) THEN
        Daux(i) = Daux(i)/Dd(i)
      ELSE
        Daux(i) = 0.0_DP
      END IF
    END DO
    
    ! Compute x = B * aux
    CALL DGEMV('n', n, n, 1.0_DP, Db, n, Daux, 1, 0.0_DP, Dx, 1)    
  END SUBROUTINE mprim_SVD_backsubst
END MODULE mprimitives

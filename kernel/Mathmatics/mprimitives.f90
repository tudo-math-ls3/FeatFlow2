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
!# </purpose>
!##############################################################################

MODULE mprimitives

  USE fsystem
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<function>
  
  PURE REAL(DP) FUNCTION mprim_getParabolicProfile (dpos,dlength,dmaxvalue)
  
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
    mprim_getParabolicProfile = 4*dmaxvalue*dpos*(dlength-dpos)/(dlength*dlength)
    
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
  
  PURE SUBROUTINE mprim_linearRescale(dx,da,db,dc,dd,dy)
  
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

END MODULE mprimitives

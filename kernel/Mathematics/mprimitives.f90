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
!#  1.) mprim_getParabolicProfile
!#      -> Calculates the value of a parabolic profile along a line.
!#
!#  2.) mprim_invertMatrix
!#      -> Invert a full matrix
!#
!#  3.) mprim_kronecker
!#      -> Compute Kronecker delta symbol
!#
!#  4.) mprim_invert2x2MatrixDirect
!#      -> Inverts a 2x2 matrix directly without pivoting.
!#
!#  5.) mprim_invert3x3MatrixDirect
!#      -> Inverts a 3x3 matrix directly without pivoting.
!#
!#  6.) mprim_invert4x4MatrixDirect
!#      -> Inverts a 4x4 matrix directly without pivoting.
!#
!#  7.) mprim_invert5x5MatrixDirect
!#      -> Inverts a 5x5 matrix directly without pivoting.
!#
!#  8.) mprim_invert6x6MatrixDirect
!#      -> Inverts a 6x6 matrix directly without pivoting.
!#
!#  9.) mprim_invertMatrixPivot
!#      -> Inverts a n x n matrix directly with pivoting.
!#
!# 10.) mprim_signum
!#      -> Signum function
!#
!# 11.) mprim_linearRescale
!#      -> Scales a coordinate x linearly from the interval [a,b] to the
!#         interval [c,d]
!#
!# 12.) mprim_quadraticInterpolation
!#     -> Evaluate the quadratic interpolation polynomial of three values.
!#
!# 13.) mprim_SVD_factorise
!#     -> Compute the factorisation for a singular value decomposition
!#
!# 14.) mprim_SVD_backsubst
!#      -> Perform back substitution for a singular value decomposition
!#
!# 15.) mprim_stdDeviation
!#      -> Calculates the standard deviation of a vector
!#
!# 16.) mprim_meanDeviation
!#      -> Calculates the mean deviation of a vector
!#
!# 17.) mprim_meanValue
!#      -> Calculates the mean value of a vector
!#
!# 18.) mprim_degToRad
!#      -> Converts DEG to RAD
!#
!# 19.) mprim_radToDeg
!#      -> Converts RAD to DEG
!#
!# 20.) mprim_solve2x2Direct
!#      -> Solves a 2x2 matrix directly without pivoting.
!#
!# 21.) mprim_solve3x3Direct
!#      -> Solves a 3x3 matrix directly without pivoting.
!#
!# 22.) mprim_solve2x2BandDiag
!#      -> Solves a 2x2 block system containing only diagonal bands
!#
!# 23.) mprim_minmod
!#      -> Computes the minmod function for two and three parameters
!#
!# 24.) mprim_softmax
!#      -> Computes the soft-maximum function for two parameter
!#
!# 25.) mprim_softmin
!#      -> Computes the soft-minimum function for two parameter
!#
!# 26.) mprim_leastSquaresMin
!#      -> Solves a least-squares minimisation problem
!#
!# 27.) mprim_transposeMatrix
!#      -> Computes the transpose of a matrix
!#
!# 28.) mprim_horner
!#      -> Evaluates a polynomial using the horner scheme.
!#
!# 29.) mprim_hornerd1
!#      -> Evaluates the derivative of a polynomial using the horner scheme.
!#
!# 30.) mprim_hornerd2
!#      -> Evaluates the 2nd derivative of a polynomial using the horner scheme.
!#
!# 31.) mprim_polarToCartesian
!#      -> Converts polar coordinates to cartesian coordinates
!#
!# 32.) mprim_cartesianToPolar
!#      -> Converts cartesian coordinates to polar coordinates
!# </purpose>
!##############################################################################

module mprimitives

!$use omp_lib
  use fsystem
  use genoutput
  use linearalgebra

  implicit none

  private

  interface mprim_signum
    module procedure mprim_signumDP
    module procedure mprim_signumSP
    module procedure mprim_signumInt
  end interface

  public :: mprim_signum

  interface mprim_SVD_factorise
    module procedure mprim_SVD_factoriseDP
    module procedure mprim_SVD_factoriseSP
  end interface mprim_SVD_factorise

  public :: mprim_SVD_factorise
  public :: mprim_SVD_factoriseDP
  public :: mprim_SVD_factoriseSP

  interface mprim_SVD_backsubst
    module procedure mprim_SVD_backsubstDP
    module procedure mprim_SVD_backsubstSP
  end interface

  public :: mprim_SVD_backsubst
  public :: mprim_SVD_backsubstDP
  public :: mprim_SVD_backsubstSP

  ! Alternative name for backward compatibility
  interface mprim_kronecker
    module procedure kronecker
  end interface

  public :: mprim_kronecker

  interface mprim_invertMatrix
    module procedure mprim_invertMatrixDP
    module procedure mprim_invertMatrixSP
  end interface

  public :: mprim_invertMatrix
  public :: mprim_invertMatrixDP
  public :: mprim_invertMatrixSP
  
  interface mprim_invertMatrixPivot
    module procedure mprim_invertMatrixPivotDP
    module procedure mprim_invertMatrixPivotSP
  end interface

  public :: mprim_invertMatrixPivot
  public :: mprim_invertMatrixPivotDP
  public :: mprim_invertMatrixPivotSP
  
  interface mprim_invert2x2MatrixDirect
    module procedure mprim_invert2x2MatrixDirectDP
    module procedure mprim_invert2x2MatrixDirectSP
  end interface

  public :: mprim_invert2x2MatrixDirect
  public :: mprim_invert2x2MatrixDirectDP
  public :: mprim_invert2x2MatrixDirectSP

  interface mprim_invert3x3MatrixDirect
    module procedure mprim_invert3x3MatrixDirectDP
    module procedure mprim_invert3x3MatrixDirectSP
  end interface

  public :: mprim_invert3x3MatrixDirect
  public :: mprim_invert3x3MatrixDirectDP
  public :: mprim_invert3x3MatrixDirectSP

  interface mprim_invert4x4MatrixDirect
    module procedure mprim_invert4x4MatrixDirectDP
    module procedure mprim_invert4x4MatrixDirectSP
  end interface
  
  public :: mprim_invert4x4MatrixDirect
  public :: mprim_invert4x4MatrixDirectDP
  public :: mprim_invert4x4MatrixDirectSP

  interface mprim_invert5x5MatrixDirect
    module procedure mprim_invert5x5MatrixDirectDP
    module procedure mprim_invert5x5MatrixDirectSP
  end interface

  public :: mprim_invert5x5MatrixDirect
  public :: mprim_invert5x5MatrixDirectDP
  public :: mprim_invert5x5MatrixDirectSP

  interface mprim_invert6x6MatrixDirect
    module procedure mprim_invert6x6MatrixDirectDP
    module procedure mprim_invert6x6MatrixDirectSP
  end interface

  public :: mprim_invert6x6MatrixDirect
  public :: mprim_invert6x6MatrixDirectDP
  public :: mprim_invert6x6MatrixDirectSP

  interface mprim_linearRescale
    module procedure mprim_linearRescaleDP
    module procedure mprim_linearRescaleSP
  end interface mprim_linearRescale

  public :: mprim_linearRescale

  interface mprim_quadraticInterpolation
    module procedure mprim_quadraticInterpolationDP
    module procedure mprim_quadraticInterpolationSP
  end interface mprim_quadraticInterpolation

  public :: mprim_quadraticInterpolation  

  interface mprim_stdDeviation
    module procedure mprim_stdDeviationDP
    module procedure mprim_stdDeviationSP
  end interface

  public :: mprim_stdDeviation
  
  interface mprim_meanDeviation
    module procedure mprim_meanDeviationDP
    module procedure mprim_meanDeviationSP
  end interface

  public :: mprim_meanDeviation

  interface mprim_meanValue
    module procedure mprim_meanValueDP
    module procedure mprim_meanValueSP
  end interface

  public :: mprim_meanValue

  interface mprim_solve2x2Direct
    module procedure mprim_solve2x2DirectDP
    module procedure mprim_solve2x2DirectSP
  end interface

  public :: mprim_solve2x2Direct

  interface mprim_solve3x3Direct
    module procedure mprim_solve3x3DirectDP
    module procedure mprim_solve3x3DirectSP
  end interface

  public :: mprim_solve3x3Direct

  interface mprim_solve2x2BandDiag
    module procedure mprim_solve2x2BandDiagDP
    module procedure mprim_solve2x2BandDiagSP
  end interface

  public :: mprim_solve2x2BandDiag

  interface mprim_minmod
    module procedure mprim_minmod2DP
    module procedure mprim_minmod2SP
    module procedure mprim_minmod3DP
    module procedure mprim_minmod3SP
  end interface

  public :: mprim_minmod

  interface mprim_softmax
    module procedure mprim_softmax2DP
    module procedure mprim_softmax3DP
    module procedure mprim_softmax2SP
    module procedure mprim_softmax3SP
  end interface

  public :: mprim_softmax

  interface mprim_softmin
    module procedure mprim_softmin2DP
    module procedure mprim_softmin3DP
    module procedure mprim_softmin2SP
    module procedure mprim_softmin3SP
  end interface

  public :: mprim_softmin

  interface mprim_leastSquaresMin
    module procedure mprim_leastSquaresMinDP
    module procedure mprim_leastSquaresMinSP
  end interface mprim_leastSquaresMin

  public :: mprim_leastSquaresMin
  public :: mprim_leastSquaresMinDP
  public :: mprim_leastSquaresMinSP

  interface mprim_transposeMatrix
    module procedure mprim_transposeMatrix1DP
    module procedure mprim_transposeMatrix1SP
    module procedure mprim_transposeMatrix2DP
    module procedure mprim_transposeMatrix2SP
  end interface mprim_transposeMatrix

  public :: mprim_transposeMatrix

  public :: mprim_getParabolicProfile
  public :: mprim_degToRad
  public :: mprim_radToDeg
  
  public :: mprim_horner
  public :: mprim_hornerd1
  public :: mprim_hornerd2

  interface mprim_horner
    module procedure mprim_hornerSP
    module procedure mprim_hornerDP
  end interface

  interface mprim_hornerd1
    module procedure mprim_hornerd1SP
    module procedure mprim_hornerd1DP
  end interface

  interface mprim_hornerd2
    module procedure mprim_hornerd2SP
    module procedure mprim_hornerd2DP
  end interface
  
  public :: mprim_polarToCartesian
  public :: mprim_cartesianToPolar
  
  interface mprim_polarToCartesian
    module procedure mprim_polarToCartesian2D
    module procedure mprim_polarToCartesian3D
  end interface
  
  interface mprim_cartesianToPolar
    module procedure mprim_cartesianToPolar2D
    module procedure mprim_cartesianToPolar3D
  end interface

contains

  ! ***************************************************************************

!<function>

  elemental real(DP) function mprim_getParabolicProfile (dpos,dlength,dmaxvalue)

!<description>
  ! Calculates the value of a parabolic profile along a line.
  ! dpos is a parameter value along a line of length dlength. The return value
  ! is the value of the profile and has its maximum dmax at position 0.5.
!</description>

!<input>

  ! Position on the line segment where to calculate the velocity
  real(DP), intent(in) :: dpos

  ! Length of the line segment
  real(DP), intent(in) :: dlength

  ! Maximum value of the profile
  real(DP), intent(in) :: dmaxvalue

!</input>

!<result>
  ! The value of the parabolic profile at position <tex>$dpos \in [0,dmax]$</tex>.
!</result>
!</function>

    ! Result: Value of the parbolic profile on position dpos on the line segment
    mprim_getParabolicProfile = 4.0_DP*dmaxvalue*dpos*(dlength-dpos)/(dlength*dlength)

  end function

  ! ***************************************************************************

!<function>

  elemental real(DP) function mprim_signumDP (dval)

!<description>
  ! Signum function
!</description>

!<input>
  ! Value to be checked
  real(DP), intent(in) :: dval
!</input>

!<result>
  ! The value of the parabolic profile at position <tex>$dpos \in [0,dmax]$</tex>.
!</result>
!</function>

    ! Result: Value of the parbolic profile on position dpos on the line segment
    if (dval .lt. 0.0_DP) then
      mprim_signumDP = -1.0_DP
    else if (dval .gt. 0.0_DP) then
      mprim_signumDP = 1.0_DP
    else
      mprim_signumDP = 0.0_DP
    end if

  end function

  ! ***************************************************************************

!<function>

  elemental real(SP) function mprim_signumSP (fval)

!<description>
  ! Signum function.
!</description>

!<input>
  ! Value to be checked
  real(SP), intent(in) :: fval
!</input>

!<result>
  ! The value of the parabolic profile at position <tex>$dpos \in [0,dmax]$</tex>.
!</result>
!</function>

    ! Result: Value of the parbolic profile on position dpos on the line segment
    if (fval .lt. 0.0_SP) then
      mprim_signumSP = -1.0_SP
    else if (fval .gt. 0.0_SP) then
      mprim_signumSP = 1.0_SP
    else
      mprim_signumSP = 0.0_SP
    end if

  end function

  ! ***************************************************************************

!<function>

  elemental integer function mprim_signumInt (ival)

!<description>
  ! Signum function.
!</description>

!<input>
  ! Value to be checked
  integer(I32), intent(in) :: ival
!</input>

!<result>
  ! The value of the parabolic profile at position <tex>$dpos \in [0,dmax]$</tex>.
!</result>
!</function>

    ! Result: Value of the parbolic profile on position dpos on the line segment
    select case (ival)
      case (:-1)
        mprim_signumInt = -1
      case (0)
        mprim_signumInt = 0
      case default
        mprim_signumInt = 1
    end select

  end function

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invertMatrixDP(Da,Df,Dx,ndim,ipar,bsuccess)

!<description>
    ! This subroutine performs the direct inversion of a NxN system.
    !
    ! If the parameter ipar=0, then only factorisation of matrix A is performed.
    ! For ipar=1, the vector x is calculated using the factorised matrix A.
    ! For ipar=2, LAPACK routine DGESV is used to solve the dense linear system Ax=f.
    ! In addition, for NDIM=2,3,4 explicit formulas are employed to replace the
    ! more expensive LAPACK routine.
!</description>

!<input>
    ! dimension of the matrix
    integer, intent(in) :: ndim

    ! source right-hand side vector
    real(DP), dimension(ndim), intent(in) :: Df

    ! What to do?
    ! IPAR = 0 : invert matrix by means of Gaussian elimination with
    !            full pivoting and return the inverted matrix inv(A)
    ! IPAR = 1 : apply inverted matrix to the right-hand side vector
    !            and return x = inv(A)*f
    ! IPAR = 2 : invert matrix and apply it to the right-hand side
    !            vector. Return x = inv(A)*f
    integer, intent(in) :: ipar
!</input>

!<inputoutput>
    ! source square matrix to be inverted
    real(DP), dimension(ndim,ndim), intent(inout) :: Da
!</inputoutput>

!<output>
    ! destination vector containing inverted matrix times right-hand
    ! side vector
    real(DP), dimension(ndim), intent(out) :: Dx

    ! TRUE, if successful. FALSE if the system is indefinite.
    ! If FALSE, Db is undefined.
    logical, intent(out) :: bsuccess
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(ndim,ndim) :: Db
    integer, dimension(ndim) :: Ipiv
    integer, dimension(ndim) :: Kindx,Kindy

    real(DP) :: dpivot,daux
    integer :: idim1,idim2,ix,iy,indx,indy,info

    interface
      pure subroutine DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      use fsystem
      integer, intent(in) :: N,LDA,LDB,NRHS
      integer, intent(inout) :: INFO
      integer, dimension(*), intent(inout) :: IPIV
      real(dp), dimension( LDA, * ), intent(inout) :: A
      real(dp), dimension( LDB, * ), intent(inout) :: B
      end subroutine
    end interface

    bsuccess = .false.

    select case (ipar)
    case (0)
      ! Perform factorisation of matrix Da

      ! Initialization
      Kindx=0;  Kindy=0

      do idim1=1,ndim

        ! Determine pivotal element
        dpivot=0

        do iy=1,ndim
          if (Kindy(iy) /= 0) cycle

          do ix=1,ndim
            if (Kindx(ix) /= 0) cycle

            if (abs(Da(ix,iy)) .le. abs(dpivot)) cycle
            dpivot=Da(ix,iy);  indx=ix;  indy=iy
          end do
        end do

        ! Return if pivotal element is zero
        if (abs(dpivot) .le. 0._DP) return

        Kindx(indx)=indy;  Kindy(indy)=indx;  Da(indx,indy)=1._DP&
            &/dpivot

        do idim2=1,ndim
          if (idim2 == indy) cycle
          Da(1:indx-1,idim2)=Da(1:indx-1,idim2)-Da(1:indx-1,  &
              & indy)*Da(indx,idim2)/dpivot
          Da(indx+1:ndim,idim2)=Da(indx+1:ndim,idim2)-Da(indx+1:ndim&
              &,indy)*Da(indx,idim2)/dpivot
        end do

        do ix=1,ndim
          if (ix /= indx) Da(ix,indy)=Da(ix,indy)/dpivot
        end do

        do iy=1,ndim
          if (iy /= indy) Da(indx,iy)=-Da(indx,iy)/dpivot
        end do
      end do

      do ix=1,ndim
        if (Kindx(ix) == ix) cycle

        do iy=1,ndim
          if (Kindx(iy) == ix) exit
        end do

        do idim1=1,ndim
          daux=Da(ix,idim1)
          Da(ix,idim1)=Da(iy,idim1)
          Da(iy,idim1)=daux
        end do

        Kindx(iy)=Kindx(ix);  Kindx(ix)=ix
      end do

      do ix=1,ndim
        if (Kindy(ix) == ix) cycle

        do iy=1,ndim
          if (Kindy(iy) == ix) exit
        end do

        do idim1=1,ndim
          daux=Da(idim1,ix)
          Da(idim1,ix)=Da(idim1,iy)
          Da(idim1,iy)=daux
        end do

        Kindy(iy)=Kindy(ix);  Kindy(ix)=ix
      end do


    case (1)
      ! Perform inversion of Da to solve the system Da * Dx = Df
      do idim1=1,ndim
        Dx(idim1)=0
        do idim2=1,ndim
          Dx(idim1)=Dx(idim1)+Da(idim1,idim2)*Df(idim2)
        end do
      end do

    case (2)
      ! Solve the dense linear system Ax=f calling LAPACK routine

      select case(ndim)
      case (1)
        if (Da(1,1) .ne. 0.0_DP) then
          Dx(1) = Df(1) / Da(1,1)
        else
          return
        end if

      case (2)
        call mprim_invert2x2MatrixDirectDP(Da,Db,bsuccess)
        if (.not. bsuccess) return
        ! Dx=matmul(Db,Df)
        Dx(1) = Db(1,1)*Df(1) + Db(1,2)*Df(2)
        Dx(2) = Db(2,1)*Df(1) + Db(2,2)*Df(2)

      case (3)
        call mprim_invert3x3MatrixDirectDP(Da,Db, bsuccess)
        if (.not. bsuccess) return
        ! Dx=matmul(Db,Df)
        Dx(1) = Db(1,1)*Df(1) + Db(1,2)*Df(2) + Db(1,3)*Df(3)
        Dx(2) = Db(2,1)*Df(1) + Db(2,2)*Df(2) + Db(2,3)*Df(3)
        Dx(3) = Db(3,1)*Df(1) + Db(3,2)*Df(2) + Db(3,3)*Df(3)

      case (4)
        call mprim_invert4x4MatrixDirectDP(Da,Db,bsuccess)
        if (.not. bsuccess) return
        ! Dx=matmul(Db,Df)
        Dx(1) = Db(1,1)*Df(1) + Db(1,2)*Df(2) &
              + Db(1,3)*Df(3) + Db(1,4)*Df(4)
        Dx(2) = Db(2,1)*Df(1) + Db(2,2)*Df(2) &
              + Db(2,3)*Df(3) + Db(2,4)*Df(4)
        Dx(3) = Db(3,1)*Df(1) + Db(3,2)*Df(2) &
              + Db(3,3)*Df(3) + Db(3,4)*Df(4)
        Dx(4) = Db(4,1)*Df(1) + Db(4,2)*Df(2) &
              + Db(4,3)*Df(3) + Db(4,4)*Df(4)

      case (5)
        call mprim_invert5x5MatrixDirectDP(Da,Db,bsuccess)
        if (.not. bsuccess) return
        ! Dx=matmul(Db,Df)
        Dx(1) = Db(1,1)*Df(1) + Db(1,2)*Df(2) &
              + Db(1,3)*Df(3) + Db(1,4)*Df(4) &
              + Db(1,5)*Df(5)
        Dx(2) = Db(2,1)*Df(1) + Db(2,2)*Df(2) &
              + Db(2,3)*Df(3) + Db(2,4)*Df(4) &
              + Db(2,5)*Df(5)
        Dx(3) = Db(3,1)*Df(1) + Db(3,2)*Df(2) &
              + Db(3,3)*Df(3) + Db(3,4)*Df(4) &
              + Db(3,5)*Df(5)
        Dx(4) = Db(4,1)*Df(1) + Db(4,2)*Df(2) &
              + Db(4,3)*Df(3) + Db(4,4)*Df(4) &
              + Db(4,5)*Df(5)
        Dx(5) = Db(5,1)*Df(1) + Db(5,2)*Df(2) &
              + Db(5,3)*Df(3) + Db(5,4)*Df(4) &
              + Db(5,5)*Df(5)

      case (6)
        call mprim_invert6x6MatrixDirectDP(Da,Db,bsuccess)
        if (.not. bsuccess) return
        ! Dx=matmul(Db,Df)
        Dx(1) = Db(1,1)*Df(1) + Db(1,2)*Df(2) &
              + Db(1,3)*Df(3) + Db(1,4)*Df(4) &
              + Db(1,5)*Df(5) + Db(1,6)*Df(6)
        Dx(2) = Db(2,1)*Df(1) + Db(2,2)*Df(2) &
              + Db(2,3)*Df(3) + Db(2,4)*Df(4) &
              + Db(2,5)*Df(5) + Db(2,6)*Df(6)
        Dx(3) = Db(3,1)*Df(1) + Db(3,2)*Df(2) &
              + Db(3,3)*Df(3) + Db(3,4)*Df(4) &
              + Db(3,5)*Df(5) + Db(3,6)*Df(6)
        Dx(4) = Db(4,1)*Df(1) + Db(4,2)*Df(2) &
              + Db(4,3)*Df(3) + Db(4,4)*Df(4) &
              + Db(4,5)*Df(5) + Db(4,6)*Df(6)
        Dx(5) = Db(5,1)*Df(1) + Db(5,2)*Df(2) &
              + Db(5,3)*Df(3) + Db(5,4)*Df(4) &
              + Db(5,5)*Df(5) + Db(5,6)*Df(6)
        Dx(6) = Db(6,1)*Df(1) + Db(6,2)*Df(2) &
              + Db(6,3)*Df(3) + Db(6,4)*Df(4) &
              + Db(6,5)*Df(5) + Db(6,6)*Df(6)

      case default
        ! Use LAPACK routine for general NxN system, where N > 6
        Ipiv=0; Dx=Df
        call DGESV(ndim,1,Da,ndim,Ipiv,Dx,ndim,info)

        if (info .ne. 0) return

      end select

    end select

    bsuccess = .true.

  end subroutine mprim_invertMatrixDP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invertMatrixSP(Fa,Ff,Fx,ndim,ipar,bsuccess)

!<description>
    ! This subroutine performs the direct inversion of a NxN system.
    !
    ! If the parameter ipar=0, then only factorisation of matrix A is performed.
    ! For ipar=1, the vector x is calculated using the factorised matrix A.
    ! For ipar=2, LAPACK routine DGESV is used to solve the dense linear system Ax=f.
    ! In addition, for NDIM=2,3,4 explicit formulas are employed to replace the
    ! more expensive LAPACK routine.
!</description>

!<input>
    ! dimension of the matrix
    integer, intent(in) :: ndim

    ! source right-hand side vector
    real(SP), dimension(ndim), intent(in) :: Ff

    ! What to do?
    ! IPAR = 0 : invert matrix by means of Gaussian elimination with
    !            full pivoting and return the inverted matrix inv(A)
    ! IPAR = 1 : apply inverted matrix to the right-hand side vector
    !            and return x = inv(A)*f
    ! IPAR = 2 : invert matrix and apply it to the right-hand side
    !            vector. Return x = inv(A)*f
    integer, intent(in) :: ipar
!</input>

!<inputoutput>
    ! source square matrix to be inverted
    real(SP), dimension(ndim,ndim), intent(inout) :: Fa
!</inputoutput>

!<output>
    ! destination vector containing inverted matrix times right-hand
    ! side vector
    real(SP), dimension(ndim), intent(out) :: Fx

    ! TRUE, if successful. FALSE if the system is indefinite.
    ! If FALSE, Fb is undefined.
    logical, intent(out) :: bsuccess
!</output>

!</subroutine>

    ! local variables
    real(SP), dimension(ndim,ndim) :: Fb
    integer, dimension(ndim) :: Ipiv
    integer, dimension(ndim) :: Kindx,Kindy

    real(SP) :: fpivot,faux
    integer :: idim1,idim2,ix,iy,indx,indy,info

    interface
      pure subroutine SGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      use fsystem
      integer, intent(in) :: N,LDA,LDB,NRHS
      integer, intent(inout) :: INFO
      integer, dimension(*), intent(inout) :: IPIV
      real(sp), dimension( LDA, * ), intent(inout) :: A
      real(sp), dimension( LDB, * ), intent(inout) :: B
      end subroutine
    end interface

    bsuccess = .false.

    select case (ipar)
    case (0)
      ! Perform factorisation of matrix Fa

      ! Initialization
      Kindx=0;  Kindy=0

      do idim1=1,ndim

        ! Determine pivotal element
        fpivot=0

        do iy=1,ndim
          if (Kindy(iy) /= 0) cycle

          do ix=1,ndim
            if (Kindx(ix) /= 0) cycle

            if (abs(Fa(ix,iy)) .le. abs(fpivot)) cycle
            fpivot=Fa(ix,iy);  indx=ix;  indy=iy
          end do
        end do

        ! Return if pivotal element is zero
        if (abs(fpivot) .le. 0._SP) return

        Kindx(indx)=indy;  Kindy(indy)=indx;  Fa(indx,indy)=1._SP&
            &/fpivot

        do idim2=1,ndim
          if (idim2 == indy) cycle
          Fa(1:indx-1,idim2)=Fa(1:indx-1,idim2)-Fa(1:indx-1,  &
              & indy)*Fa(indx,idim2)/fpivot
          Fa(indx+1:ndim,idim2)=Fa(indx+1:ndim,idim2)-Fa(indx+1:ndim&
              &,indy)*Fa(indx,idim2)/fpivot
        end do

        do ix=1,ndim
          if (ix /= indx) Fa(ix,indy)=Fa(ix,indy)/fpivot
        end do

        do iy=1,ndim
          if (iy /= indy) Fa(indx,iy)=-Fa(indx,iy)/fpivot
        end do
      end do

      do ix=1,ndim
        if (Kindx(ix) == ix) cycle

        do iy=1,ndim
          if (Kindx(iy) == ix) exit
        end do

        do idim1=1,ndim
          faux=Fa(ix,idim1)
          Fa(ix,idim1)=Fa(iy,idim1)
          Fa(iy,idim1)=faux
        end do

        Kindx(iy)=Kindx(ix);  Kindx(ix)=ix
      end do

      do ix=1,ndim
        if (Kindy(ix) == ix) cycle

        do iy=1,ndim
          if (Kindy(iy) == ix) exit
        end do

        do idim1=1,ndim
          faux=Fa(idim1,ix)
          Fa(idim1,ix)=Fa(idim1,iy)
          Fa(idim1,iy)=faux
        end do

        Kindy(iy)=Kindy(ix);  Kindy(ix)=ix
      end do


    case (1)
      ! Perform inversion of Fa to solve the system Fa * Fx = Ff
      do idim1=1,ndim
        Fx(idim1)=0
        do idim2=1,ndim
          Fx(idim1)=Fx(idim1)+Fa(idim1,idim2)*Ff(idim2)
        end do
      end do

    case (2)
      ! Solve the dense linear system Ax=f calling LAPACK routine

      select case(ndim)
      case (1)
        if (Fa(1,1) .ne. 0.0_SP) then
          Fx(1) = Ff(1) / Fa(1,1)
        else
          return
        end if

      case (2)
        call mprim_invert2x2MatrixDirectSP(Fa,Fb,bsuccess)
        if (.not. bsuccess) return
        ! Fx=matmul(Fb,Ff)
        Fx(1) = Fb(1,1)*Ff(1) + Fb(1,2)*Ff(2)
        Fx(2) = Fb(2,1)*Ff(1) + Fb(2,2)*Ff(2)

      case (3)
        call mprim_invert3x3MatrixDirectSP(Fa,Fb, bsuccess)
        if (.not. bsuccess) return
        ! Fx=matmul(Fb,Ff)
        Fx(1) = Fb(1,1)*Ff(1) + Fb(1,2)*Ff(2) + Fb(1,3)*Ff(3)
        Fx(2) = Fb(2,1)*Ff(1) + Fb(2,2)*Ff(2) + Fb(2,3)*Ff(3)
        Fx(3) = Fb(3,1)*Ff(1) + Fb(3,2)*Ff(2) + Fb(3,3)*Ff(3)

      case (4)
        call mprim_invert4x4MatrixDirectSP(Fa,Fb,bsuccess)
        if (.not. bsuccess) return
        ! Fx=matmul(Fb,Ff)
        Fx(1) = Fb(1,1)*Ff(1) + Fb(1,2)*Ff(2) &
              + Fb(1,3)*Ff(3) + Fb(1,4)*Ff(4)
        Fx(2) = Fb(2,1)*Ff(1) + Fb(2,2)*Ff(2) &
              + Fb(2,3)*Ff(3) + Fb(2,4)*Ff(4)
        Fx(3) = Fb(3,1)*Ff(1) + Fb(3,2)*Ff(2) &
              + Fb(3,3)*Ff(3) + Fb(3,4)*Ff(4)
        Fx(4) = Fb(4,1)*Ff(1) + Fb(4,2)*Ff(2) &
              + Fb(4,3)*Ff(3) + Fb(4,4)*Ff(4)

      case (5)
        call mprim_invert5x5MatrixDirectSP(Fa,Fb,bsuccess)
        if (.not. bsuccess) return
        ! Fx=matmul(Fb,Ff)
        Fx(1) = Fb(1,1)*Ff(1) + Fb(1,2)*Ff(2) &
              + Fb(1,3)*Ff(3) + Fb(1,4)*Ff(4) &
              + Fb(1,5)*Ff(5)
        Fx(2) = Fb(2,1)*Ff(1) + Fb(2,2)*Ff(2) &
              + Fb(2,3)*Ff(3) + Fb(2,4)*Ff(4) &
              + Fb(2,5)*Ff(5)
        Fx(3) = Fb(3,1)*Ff(1) + Fb(3,2)*Ff(2) &
              + Fb(3,3)*Ff(3) + Fb(3,4)*Ff(4) &
              + Fb(3,5)*Ff(5)
        Fx(4) = Fb(4,1)*Ff(1) + Fb(4,2)*Ff(2) &
              + Fb(4,3)*Ff(3) + Fb(4,4)*Ff(4) &
              + Fb(4,5)*Ff(5)
        Fx(5) = Fb(5,1)*Ff(1) + Fb(5,2)*Ff(2) &
              + Fb(5,3)*Ff(3) + Fb(5,4)*Ff(4) &
              + Fb(5,5)*Ff(5)

      case (6)
        call mprim_invert6x6MatrixDirectSP(Fa,Fb,bsuccess)
        if (.not. bsuccess) return
        ! Fx=matmul(Fb,Ff)
        Fx(1) = Fb(1,1)*Ff(1) + Fb(1,2)*Ff(2) &
              + Fb(1,3)*Ff(3) + Fb(1,4)*Ff(4) &
              + Fb(1,5)*Ff(5) + Fb(1,6)*Ff(6)
        Fx(2) = Fb(2,1)*Ff(1) + Fb(2,2)*Ff(2) &
              + Fb(2,3)*Ff(3) + Fb(2,4)*Ff(4) &
              + Fb(2,5)*Ff(5) + Fb(2,6)*Ff(6)
        Fx(3) = Fb(3,1)*Ff(1) + Fb(3,2)*Ff(2) &
              + Fb(3,3)*Ff(3) + Fb(3,4)*Ff(4) &
              + Fb(3,5)*Ff(5) + Fb(3,6)*Ff(6)
        Fx(4) = Fb(4,1)*Ff(1) + Fb(4,2)*Ff(2) &
              + Fb(4,3)*Ff(3) + Fb(4,4)*Ff(4) &
              + Fb(4,5)*Ff(5) + Fb(4,6)*Ff(6)
        Fx(5) = Fb(5,1)*Ff(1) + Fb(5,2)*Ff(2) &
              + Fb(5,3)*Ff(3) + Fb(5,4)*Ff(4) &
              + Fb(5,5)*Ff(5) + Fb(5,6)*Ff(6)
        Fx(6) = Fb(6,1)*Ff(1) + Fb(6,2)*Ff(2) &
              + Fb(6,3)*Ff(3) + Fb(6,4)*Ff(4) &
              + Fb(6,5)*Ff(5) + Fb(6,6)*Ff(6)

      case default
        ! Use LAPACK routine for general NxN system, where N > 6
        Ipiv=0; Fx=Ff
        call SGESV(ndim,1,Fa,ndim,Ipiv,Fx,ndim,info)

        if (info .ne. 0) return

      end select

    end select

    bsuccess = .true.

  end subroutine mprim_invertMatrixSP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invert2x2MatrixDirectDP(Da,Db,bsuccess)

!<description>
  ! This subroutine directly inverts a 2x2 system without any pivoting.
  ! 'Da' is a 2-dimensional 2x2 m matrix. The inverse of Da is written
  ! to the 2-dimensional 2x2 matrix Db.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Da and Db are assumed to be 2x2 arrays!
!</description>

!<input>
  ! source square matrix to be inverted
  real(DP), dimension(2,2), intent(in) :: Da
!</input>

!<output>
  ! destination square matrix; receives <tex>$ A^{-1} $</tex>.
  real(DP), dimension(2,2), intent(out) :: Db

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Db is undefined.
  logical, intent(out) :: bsuccess
!</output>

!</subroutine>

  real(DP) :: daux

    ! Explicit formula for 2x2 system
    Db(1,1)= Da(2,2)
    Db(2,1)=-Da(2,1)
    Db(1,2)=-Da(1,2)
    Db(2,2)= Da(1,1)
    daux=Da(1,1)*Da(2,2)-Da(1,2)*Da(2,1)
    if (daux .ne. 0.0_DP) then
      Db=Db*(1.0_DP/daux)
      bsuccess = .true.
    else
      bsuccess = .false.
    end if

  end subroutine mprim_invert2x2MatrixDirectDP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invert2x2MatrixDirectSP(Fa,Fb,bsuccess)

!<description>
  ! This subroutine directly inverts a 2x2 system without any pivoting.
  ! 'Fa' is a 2-dimensional 2x2 m matrix. The inverse of Fa is written
  ! to the 2-dimensional 2x2 matrix Fb.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Fa and Fb are assumed to be 2x2 arrays!
!</description>

!<input>
  ! source square matrix to be inverted
  real(SP), dimension(2,2), intent(in) :: Fa
!</input>

!<output>
  ! destination square matrix; receives <tex>$ A^{-1} $</tex>.
  real(SP), dimension(2,2), intent(out) :: Fb

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Fb is undefined.
  logical, intent(out) :: bsuccess
!</output>

!</subroutine>

  real(SP) :: faux

    ! Explicit formula for 2x2 system
    Fb(1,1)= Fa(2,2)
    Fb(2,1)=-Fa(2,1)
    Fb(1,2)=-Fa(1,2)
    Fb(2,2)= Fa(1,1)
    faux=Fa(1,1)*Fa(2,2)-Fa(1,2)*Fa(2,1)
    if (faux .ne. 0.0_SP) then
      Fb=Fb*(1.0_SP/faux)
      bsuccess = .true.
    else
      bsuccess = .false.
    end if

  end subroutine mprim_invert2x2MatrixDirectSP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invert3x3MatrixDirectDP(Da,Db,bsuccess)

!<description>
  ! This subroutine directly inverts a 3x3 system without any pivoting.
  ! 'Da' is a 2-dimensional 3x3 m matrix. The inverse of Da is written
  ! to the 2-dimensional 3x3 matrix Db.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Da and Db are assumed to be 3x3 arrays!
!</description>

!<input>
  ! source square matrix to be inverted
  real(DP), dimension(3,3), intent(in) :: Da
!</input>

!<output>
  ! destination square matrix; receives <tex>$ A^{-1} $</tex>.
  real(DP), dimension(3,3), intent(out) :: Db

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Db is undefined.
  logical, intent(out) :: bsuccess
!</output>

!</subroutine>

  real(DP) :: daux

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
    daux = Da(1,1)*Db(1,1)+Da(1,2)*Db(2,1)+Da(1,3)*Db(3,1)
    if (daux .ne. 0.0_DP) then
      !daux=Da(1,1)*Da(2,2)*Da(3,3)+Da(2,1)*Da(3,2)*Da(1,3)+ Da(3,1)&
      !    &*Da(1,2)*Da(2,3)-Da(1,1)*Da(3,2)*Da(2,3)- Da(3,1)*Da(2&
      !    &,2)*Da(1,3)-Da(2,1)*Da(1,2)*Da(3,3)
      Db=Db*(1.0_DP/daux)
      bsuccess = .true.
    else
      bsuccess = .false.
    end if

  end subroutine mprim_invert3x3MatrixDirectDP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invert3x3MatrixDirectSP(Fa,Fb,bsuccess)

!<description>
  ! This subroutine directly inverts a 3x3 system without any pivoting.
  ! 'Fa' is a 2-dimensional 3x3 m matrix. The inverse of Fa is written
  ! to the 2-dimensional 3x3 matrix Fb.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Fa and Fb are assumed to be 3x3 arrays!
!</description>

!<input>
  ! source square matrix to be inverted
  real(SP), dimension(3,3), intent(in) :: Fa
!</input>

!<output>
  ! destination square matrix; receives <tex>$ A^{-1} $</tex>.
  real(SP), dimension(3,3), intent(out) :: Fb

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Fb is undefined.
  logical, intent(out) :: bsuccess
!</output>

!</subroutine>

  real(SP) :: faux

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
    faux = Fa(1,1)*Fb(1,1)+Fa(1,2)*Fb(2,1)+Fa(1,3)*Fb(3,1)
    if (faux .ne. 0.0_SP) then
      !faux=Fa(1,1)*Fa(2,2)*Fa(3,3)+Fa(2,1)*Fa(3,2)*Fa(1,3)+ Fa(3,1)&
      !    &*Fa(1,2)*Fa(2,3)-Fa(1,1)*Fa(3,2)*Fa(2,3)- Fa(3,1)*Fa(2&
      !    &,2)*Fa(1,3)-Fa(2,1)*Fa(1,2)*Fa(3,3)
      Fb=Fb*(1.0_SP/faux)
      bsuccess = .true.
    else
      bsuccess = .false.
    end if

  end subroutine mprim_invert3x3MatrixDirectSP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invert4x4MatrixDirectDP(Da,Db,bsuccess)

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
  real(DP), dimension(4,4), intent(in) :: Da
!</input>

!<output>
  ! destination square matrix; receives <tex>$ A^{-1} $</tex>.
  real(DP), dimension(4,4), intent(out) :: Db

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Db is undefined.
  logical, intent(out) :: bsuccess
!</output>

!</subroutine>

  ! auxiliary variables
  real(DP) :: det,daux
  real(DP), dimension(6) :: W

    ! 2x2 determinants of rows 3-4
    W(1)=Da(3,1)*Da(4,2)-Da(3,2)*Da(4,1)
    W(2)=Da(3,1)*Da(4,3)-Da(3,3)*Da(4,1)
    W(3)=Da(3,1)*Da(4,4)-Da(3,4)*Da(4,1)
    W(4)=Da(3,2)*Da(4,3)-Da(3,3)*Da(4,2)
    W(5)=Da(3,2)*Da(4,4)-Da(3,4)*Da(4,2)
    W(6)=Da(3,3)*Da(4,4)-Da(3,4)*Da(4,3)
    ! pre-calculate first column of inverse
    Db(1,1)= Da(2,2)*W(6)-Da(2,3)*W(5)+Da(2,4)*W(4)
    Db(2,1)=-Da(2,1)*W(6)+Da(2,3)*W(3)-Da(2,4)*W(2)
    Db(3,1)= Da(2,1)*W(5)-Da(2,2)*W(3)+Da(2,4)*W(1)
    Db(4,1)=-Da(2,1)*W(4)+Da(2,2)*W(2)-Da(2,3)*W(1)
    daux = (Da(1,1)*Db(1,1)+Da(1,2)*Db(2,1)&
            +Da(1,3)*Db(3,1)+Da(1,4)*Db(4,1))
    if (daux .ne. 0.0_DP) then
      ! calculate determinant of A
      det = 1.0_DP / daux
      ! update first column of inverse
      Db(1,1)=Db(1,1)*det
      Db(2,1)=Db(2,1)*det
      Db(3,1)=Db(3,1)*det
      Db(4,1)=Db(4,1)*det
      ! calculate second column of inverse
      Db(1,2)=det*(-Da(1,2)*W(6)+Da(1,3)*W(5)-Da(1,4)*W(4))
      Db(2,2)=det*( Da(1,1)*W(6)-Da(1,3)*W(3)+Da(1,4)*W(2))
      Db(3,2)=det*(-Da(1,1)*W(5)+Da(1,2)*W(3)-Da(1,4)*W(1))
      Db(4,2)=det*( Da(1,1)*W(4)-Da(1,2)*W(2)+Da(1,3)*W(1))
      ! 2x2 determinants of rows 1-2
      W(1)=Da(1,1)*Da(2,2)-Da(1,2)*Da(2,1)
      W(2)=Da(1,1)*Da(2,3)-Da(1,3)*Da(2,1)
      W(3)=Da(1,1)*Da(2,4)-Da(1,4)*Da(2,1)
      W(4)=Da(1,2)*Da(2,3)-Da(1,3)*Da(2,2)
      W(5)=Da(1,2)*Da(2,4)-Da(1,4)*Da(2,2)
      W(6)=Da(1,3)*Da(2,4)-Da(1,4)*Da(2,3)
      ! calculate third column of inverse
      Db(1,3)=det*( Da(4,2)*W(6)-Da(4,3)*W(5)+Da(4,4)*W(4))
      Db(2,3)=det*(-Da(4,1)*W(6)+Da(4,3)*W(3)-Da(4,4)*W(2))
      Db(3,3)=det*( Da(4,1)*W(5)-Da(4,2)*W(3)+Da(4,4)*W(1))
      Db(4,3)=det*(-Da(4,1)*W(4)+Da(4,2)*W(2)-Da(4,3)*W(1))
      ! calculate fourth column of inverse
      Db(1,4)=det*(-Da(3,2)*W(6)+Da(3,3)*W(5)-Da(3,4)*W(4))
      Db(2,4)=det*( Da(3,1)*W(6)-Da(3,3)*W(3)+Da(3,4)*W(2))
      Db(3,4)=det*(-Da(3,1)*W(5)+Da(3,2)*W(3)-Da(3,4)*W(1))
      Db(4,4)=det*( Da(3,1)*W(4)-Da(3,2)*W(2)+Da(3,3)*W(1))

      ! 'old' implementation follows

!      real(DP) :: daux
!
!        ! Explicit formula for 4x4 system
!        Db(1,1)=Da(2,2)*Da(3,3)*Da(4,4)+Da(2,3)*Da(3,4)*Da(4,2)+Da(2&
!            &,4)*Da(3,2)*Da(4,3)- Da(2,2)*Da(3,4)*Da(4,3)-Da(2,3)&
!            &*Da(3,2)*Da(4,4)-Da(2,4)*Da(3,3)*Da(4,2)
!        Db(2,1)=Da(2,1)*Da(3,4)*Da(4,3)+Da(2,3)*Da(3,1)*Da(4,4)+Da(2&
!            &,4)*Da(3,3)*Da(4,1)- Da(2,1)*Da(3,3)*Da(4,4)-Da(2,3)&
!            &*Da(3,4)*Da(4,1)-Da(2,4)*Da(3,1)*Da(4,3)
!        Db(3,1)=Da(2,1)*Da(3,2)*Da(4,4)+Da(2,2)*Da(3,4)*Da(4,1)+Da(2&
!            &,4)*Da(3,1)*Da(4,2)- Da(2,1)*Da(3,4)*Da(4,2)-Da(2,2)&
!            &*Da(3,1)*Da(4,4)-Da(2,4)*Da(3,2)*Da(4,1)
!        Db(4,1)=Da(2,1)*Da(3,3)*Da(4,2)+Da(2,2)*Da(3,1)*Da(4,3)+Da(2&
!            &,3)*Da(3,2)*Da(4,1)- Da(2,1)*Da(3,2)*Da(4,3)-Da(2,2)&
!            &*Da(3,3)*Da(4,1)-Da(2,3)*Da(3,1)*Da(4,2)
!        Db(1,2)=Da(1,2)*Da(3,4)*Da(4,3)+Da(1,3)*Da(3,2)*Da(4,4)+Da(1&
!            &,4)*Da(3,3)*Da(4,2)- Da(1,2)*Da(3,3)*Da(4,4)-Da(1,3)&
!            &*Da(3,4)*Da(4,2)-Da(1,4)*Da(3,2)*Da(4,3)
!        Db(2,2)=Da(1,1)*Da(3,3)*Da(4,4)+Da(1,3)*Da(3,4)*Da(4,1)+Da(1&
!            &,4)*Da(3,1)*Da(4,3)- Da(1,1)*Da(3,4)*Da(4,3)-Da(1,3)&
!            &*Da(3,1)*Da(4,4)-Da(1,4)*Da(3,3)*Da(4,1)
!        Db(3,2)=Da(1,1)*Da(3,4)*Da(4,2)+Da(1,2)*Da(3,1)*Da(4,4)+Da(1&
!            &,4)*Da(3,2)*Da(4,1)- Da(1,1)*Da(3,2)*Da(4,4)-Da(1,2)&
!            &*Da(3,4)*Da(4,1)-Da(1,4)*Da(3,1)*Da(4,2)
!        Db(4,2)=Da(1,1)*Da(3,2)*Da(4,3)+Da(1,2)*Da(3,3)*Da(4,1)+Da(1&
!            &,3)*Da(3,1)*Da(4,2)- Da(1,1)*Da(3,3)*Da(4,2)-Da(1,2)&
!            &*Da(3,1)*Da(4,3)-Da(1,3)*Da(3,2)*Da(4,1)
!        Db(1,3)=Da(1,2)*Da(2,3)*Da(4,4)+Da(1,3)*Da(2,4)*Da(4,2)+Da(1&
!            &,4)*Da(2,2)*Da(4,3)- Da(1,2)*Da(2,4)*Da(4,3)-Da(1,3)&
!            &*Da(2,2)*Da(4,4)-Da(1,4)*Da(2,3)*Da(4,2)
!        Db(2,3)=Da(1,1)*Da(2,4)*Da(4,3)+Da(1,3)*Da(2,1)*Da(4,4)+Da(1&
!            &,4)*Da(2,3)*Da(4,1)- Da(1,1)*Da(2,3)*Da(4,4)-Da(1,3)&
!            &*Da(2,4)*Da(4,1)-Da(1,4)*Da(2,1)*Da(4,3)
!        Db(3,3)=Da(1,1)*Da(2,2)*Da(4,4)+Da(1,2)*Da(2,4)*Da(4,1)+Da(1&
!            &,4)*Da(2,1)*Da(4,2)- Da(1,1)*Da(2,4)*Da(4,2)-Da(1,2)&
!            &*Da(2,1)*Da(4,4)-Da(1,4)*Da(2,2)*Da(4,1)
!        Db(4,3)=Da(1,1)*Da(2,3)*Da(4,2)+Da(1,2)*Da(2,1)*Da(4,3)+Da(1&
!            &,3)*Da(2,2)*Da(4,1)- Da(1,1)*Da(2,2)*Da(4,3)-Da(1,2)&
!            &*Da(2,3)*Da(4,1)-Da(1,3)*Da(2,1)*Da(4,2)
!        Db(1,4)=Da(1,2)*Da(2,4)*Da(3,3)+Da(1,3)*Da(2,2)*Da(3,4)+Da(1&
!            &,4)*Da(2,3)*Da(3,2)- Da(1,2)*Da(2,3)*Da(3,4)-Da(1,3)&
!            &*Da(2,4)*Da(3,2)-Da(1,4)*Da(2,2)*Da(3,3)
!        Db(2,4)=Da(1,1)*Da(2,3)*Da(3,4)+Da(1,3)*Da(2,4)*Da(3,1)+Da(1&
!            &,4)*Da(2,1)*Da(3,3)- Da(1,1)*Da(2,4)*Da(3,3)-Da(1,3)&
!            &*Da(2,1)*Da(3,4)-Da(1,4)*Da(2,3)*Da(3,1)
!        Db(3,4)=Da(1,1)*Da(2,4)*Da(3,2)+Da(1,2)*Da(2,1)*Da(3,4)+Da(1&
!            &,4)*Da(2,2)*Da(3,1)- Da(1,1)*Da(2,2)*Da(3,4)-Da(1,2)&
!            &*Da(2,4)*Da(3,1)-Da(1,4)*Da(2,1)*Da(3,2)
!        Db(4,4)=Da(1,1)*Da(2,2)*Da(3,3)+Da(1,2)*Da(2,3)*Da(3,1)+Da(1&
!            &,3)*Da(2,1)*Da(3,2)- Da(1,1)*Da(2,3)*Da(3,2)-Da(1,2)&
!            &*Da(2,1)*Da(3,3)-Da(1,3)*Da(2,2)*Da(3,1)
!        daux=Da(1,1)*Da(2,2)*Da(3,3)*Da(4,4)+Da(1,1)*Da(2,3)*Da(3,4)&
!            &*Da(4,2)+Da(1,1)*Da(2,4)*Da(3,2)*Da(4,3)+ Da(1,2)*Da(2&
!            &,1)*Da(3,4)*Da(4,3)+Da(1,2)*Da(2,3)*Da(3,1)*Da(4,4)+Da(1&
!            &,2)*Da(2,4)*Da(3,3)*Da(4,1)+ Da(1,3)*Da(2,1)*Da(3,2)&
!            &*Da(4,4)+Da(1,3)*Da(2,2)*Da(3,4)*Da(4,1)+Da(1,3)*Da(2,4)&
!            &*Da(3,1)*Da(4,2)+ Da(1,4)*Da(2,1)*Da(3,3)*Da(4,2)+Da(1&
!            &,4)*Da(2,2)*Da(3,1)*Da(4,3)+Da(1,4)*Da(2,3)*Da(3,2)*Da(4&
!            &,1)- Da(1,1)*Da(2,2)*Da(3,4)*Da(4,3)-Da(1,1)*Da(2,3)&
!            &*Da(3,2)*Da(4,4)-Da(1,1)*Da(2,4)*Da(3,3)*Da(4,2)- Da(1&
!            &,2)*Da(2,1)*Da(3,3)*Da(4,4)-Da(1,2)*Da(2,3)*Da(3,4)*Da(4&
!            &,1)-Da(1,2)*Da(2,4)*Da(3,1)*Da(4,3)- Da(1,3)*Da(2,1)&
!            &*Da(3,4)*Da(4,2)-Da(1,3)*Da(2,2)*Da(3,1)*Da(4,4)-Da(1,3)&
!            &*Da(2,4)*Da(3,2)*Da(4,1)- Da(1,4)*Da(2,1)*Da(3,2)*Da(4&
!            &,3)-Da(1,4)*Da(2,2)*Da(3,3)*Da(4,1)-Da(1,4)*Da(2,3)*Da(3&
!            &,1)*Da(4,2)
!        Db=Db*(1.0_DP/daux)

        bsuccess = .true.

      else
        bsuccess = .false.
      end if

  end subroutine mprim_invert4x4MatrixDirectDP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invert4x4MatrixDirectSP(Fa,Fb,bsuccess)

!<description>
  ! This subroutine directly inverts a 4x4 system without any pivoting.
  ! 'Fa' is a 2-dimensional 4x4 m matrix. The inverse of Fa is written
  ! to the 2-dimensional 4x4 matrix Fb.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Fa and Fb are assumed to be 4x4 arrays!
!</description>

!<input>
  ! source square matrix to be inverted
  real(SP), dimension(4,4), intent(in) :: Fa
!</input>

!<output>
  ! destination square matrix; receives <tex>$ A^{-1} $</tex>.
  real(SP), dimension(4,4), intent(out) :: Fb

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Fb is undefined.
  logical, intent(out) :: bsuccess
!</output>

!</subroutine>

  ! auxiliary variables
  real(SP) :: det,faux
  real(SP), dimension(6) :: W

    ! 2x2 determinants of rows 3-4
    W(1)=Fa(3,1)*Fa(4,2)-Fa(3,2)*Fa(4,1)
    W(2)=Fa(3,1)*Fa(4,3)-Fa(3,3)*Fa(4,1)
    W(3)=Fa(3,1)*Fa(4,4)-Fa(3,4)*Fa(4,1)
    W(4)=Fa(3,2)*Fa(4,3)-Fa(3,3)*Fa(4,2)
    W(5)=Fa(3,2)*Fa(4,4)-Fa(3,4)*Fa(4,2)
    W(6)=Fa(3,3)*Fa(4,4)-Fa(3,4)*Fa(4,3)
    ! pre-calculate first column of inverse
    Fb(1,1)= Fa(2,2)*W(6)-Fa(2,3)*W(5)+Fa(2,4)*W(4)
    Fb(2,1)=-Fa(2,1)*W(6)+Fa(2,3)*W(3)-Fa(2,4)*W(2)
    Fb(3,1)= Fa(2,1)*W(5)-Fa(2,2)*W(3)+Fa(2,4)*W(1)
    Fb(4,1)=-Fa(2,1)*W(4)+Fa(2,2)*W(2)-Fa(2,3)*W(1)
    faux = (Fa(1,1)*Fb(1,1)+Fa(1,2)*Fb(2,1)&
            +Fa(1,3)*Fb(3,1)+Fa(1,4)*Fb(4,1))
    if (faux .ne. 0.0_SP) then
      ! calculate determinant of A
      det = 1.0_SP / faux
      ! update first column of inverse
      Fb(1,1)=Fb(1,1)*det
      Fb(2,1)=Fb(2,1)*det
      Fb(3,1)=Fb(3,1)*det
      Fb(4,1)=Fb(4,1)*det
      ! calculate second column of inverse
      Fb(1,2)=det*(-Fa(1,2)*W(6)+Fa(1,3)*W(5)-Fa(1,4)*W(4))
      Fb(2,2)=det*( Fa(1,1)*W(6)-Fa(1,3)*W(3)+Fa(1,4)*W(2))
      Fb(3,2)=det*(-Fa(1,1)*W(5)+Fa(1,2)*W(3)-Fa(1,4)*W(1))
      Fb(4,2)=det*( Fa(1,1)*W(4)-Fa(1,2)*W(2)+Fa(1,3)*W(1))
      ! 2x2 determinants of rows 1-2
      W(1)=Fa(1,1)*Fa(2,2)-Fa(1,2)*Fa(2,1)
      W(2)=Fa(1,1)*Fa(2,3)-Fa(1,3)*Fa(2,1)
      W(3)=Fa(1,1)*Fa(2,4)-Fa(1,4)*Fa(2,1)
      W(4)=Fa(1,2)*Fa(2,3)-Fa(1,3)*Fa(2,2)
      W(5)=Fa(1,2)*Fa(2,4)-Fa(1,4)*Fa(2,2)
      W(6)=Fa(1,3)*Fa(2,4)-Fa(1,4)*Fa(2,3)
      ! calculate third column of inverse
      Fb(1,3)=det*( Fa(4,2)*W(6)-Fa(4,3)*W(5)+Fa(4,4)*W(4))
      Fb(2,3)=det*(-Fa(4,1)*W(6)+Fa(4,3)*W(3)-Fa(4,4)*W(2))
      Fb(3,3)=det*( Fa(4,1)*W(5)-Fa(4,2)*W(3)+Fa(4,4)*W(1))
      Fb(4,3)=det*(-Fa(4,1)*W(4)+Fa(4,2)*W(2)-Fa(4,3)*W(1))
      ! calculate fourth column of inverse
      Fb(1,4)=det*(-Fa(3,2)*W(6)+Fa(3,3)*W(5)-Fa(3,4)*W(4))
      Fb(2,4)=det*( Fa(3,1)*W(6)-Fa(3,3)*W(3)+Fa(3,4)*W(2))
      Fb(3,4)=det*(-Fa(3,1)*W(5)+Fa(3,2)*W(3)-Fa(3,4)*W(1))
      Fb(4,4)=det*( Fa(3,1)*W(4)-Fa(3,2)*W(2)+Fa(3,3)*W(1))

      ! 'old' implementation follows

!      real(SP) :: faux
!
!        ! Explicit formula for 4x4 system
!        Fb(1,1)=Fa(2,2)*Fa(3,3)*Fa(4,4)+Fa(2,3)*Fa(3,4)*Fa(4,2)+Fa(2&
!            &,4)*Fa(3,2)*Fa(4,3)- Fa(2,2)*Fa(3,4)*Fa(4,3)-Fa(2,3)&
!            &*Fa(3,2)*Fa(4,4)-Fa(2,4)*Fa(3,3)*Fa(4,2)
!        Fb(2,1)=Fa(2,1)*Fa(3,4)*Fa(4,3)+Fa(2,3)*Fa(3,1)*Fa(4,4)+Fa(2&
!            &,4)*Fa(3,3)*Fa(4,1)- Fa(2,1)*Fa(3,3)*Fa(4,4)-Fa(2,3)&
!            &*Fa(3,4)*Fa(4,1)-Fa(2,4)*Fa(3,1)*Fa(4,3)
!        Fb(3,1)=Fa(2,1)*Fa(3,2)*Fa(4,4)+Fa(2,2)*Fa(3,4)*Fa(4,1)+Fa(2&
!            &,4)*Fa(3,1)*Fa(4,2)- Fa(2,1)*Fa(3,4)*Fa(4,2)-Fa(2,2)&
!            &*Fa(3,1)*Fa(4,4)-Fa(2,4)*Fa(3,2)*Fa(4,1)
!        Fb(4,1)=Fa(2,1)*Fa(3,3)*Fa(4,2)+Fa(2,2)*Fa(3,1)*Fa(4,3)+Fa(2&
!            &,3)*Fa(3,2)*Fa(4,1)- Fa(2,1)*Fa(3,2)*Fa(4,3)-Fa(2,2)&
!            &*Fa(3,3)*Fa(4,1)-Fa(2,3)*Fa(3,1)*Fa(4,2)
!        Fb(1,2)=Fa(1,2)*Fa(3,4)*Fa(4,3)+Fa(1,3)*Fa(3,2)*Fa(4,4)+Fa(1&
!            &,4)*Fa(3,3)*Fa(4,2)- Fa(1,2)*Fa(3,3)*Fa(4,4)-Fa(1,3)&
!            &*Fa(3,4)*Fa(4,2)-Fa(1,4)*Fa(3,2)*Fa(4,3)
!        Fb(2,2)=Fa(1,1)*Fa(3,3)*Fa(4,4)+Fa(1,3)*Fa(3,4)*Fa(4,1)+Fa(1&
!            &,4)*Fa(3,1)*Fa(4,3)- Fa(1,1)*Fa(3,4)*Fa(4,3)-Fa(1,3)&
!            &*Fa(3,1)*Fa(4,4)-Fa(1,4)*Fa(3,3)*Fa(4,1)
!        Fb(3,2)=Fa(1,1)*Fa(3,4)*Fa(4,2)+Fa(1,2)*Fa(3,1)*Fa(4,4)+Fa(1&
!            &,4)*Fa(3,2)*Fa(4,1)- Fa(1,1)*Fa(3,2)*Fa(4,4)-Fa(1,2)&
!            &*Fa(3,4)*Fa(4,1)-Fa(1,4)*Fa(3,1)*Fa(4,2)
!        Fb(4,2)=Fa(1,1)*Fa(3,2)*Fa(4,3)+Fa(1,2)*Fa(3,3)*Fa(4,1)+Fa(1&
!            &,3)*Fa(3,1)*Fa(4,2)- Fa(1,1)*Fa(3,3)*Fa(4,2)-Fa(1,2)&
!            &*Fa(3,1)*Fa(4,3)-Fa(1,3)*Fa(3,2)*Fa(4,1)
!        Fb(1,3)=Fa(1,2)*Fa(2,3)*Fa(4,4)+Fa(1,3)*Fa(2,4)*Fa(4,2)+Fa(1&
!            &,4)*Fa(2,2)*Fa(4,3)- Fa(1,2)*Fa(2,4)*Fa(4,3)-Fa(1,3)&
!            &*Fa(2,2)*Fa(4,4)-Fa(1,4)*Fa(2,3)*Fa(4,2)
!        Fb(2,3)=Fa(1,1)*Fa(2,4)*Fa(4,3)+Fa(1,3)*Fa(2,1)*Fa(4,4)+Fa(1&
!            &,4)*Fa(2,3)*Fa(4,1)- Fa(1,1)*Fa(2,3)*Fa(4,4)-Fa(1,3)&
!            &*Fa(2,4)*Fa(4,1)-Fa(1,4)*Fa(2,1)*Fa(4,3)
!        Fb(3,3)=Fa(1,1)*Fa(2,2)*Fa(4,4)+Fa(1,2)*Fa(2,4)*Fa(4,1)+Fa(1&
!            &,4)*Fa(2,1)*Fa(4,2)- Fa(1,1)*Fa(2,4)*Fa(4,2)-Fa(1,2)&
!            &*Fa(2,1)*Fa(4,4)-Fa(1,4)*Fa(2,2)*Fa(4,1)
!        Fb(4,3)=Fa(1,1)*Fa(2,3)*Fa(4,2)+Fa(1,2)*Fa(2,1)*Fa(4,3)+Fa(1&
!            &,3)*Fa(2,2)*Fa(4,1)- Fa(1,1)*Fa(2,2)*Fa(4,3)-Fa(1,2)&
!            &*Fa(2,3)*Fa(4,1)-Fa(1,3)*Fa(2,1)*Fa(4,2)
!        Fb(1,4)=Fa(1,2)*Fa(2,4)*Fa(3,3)+Fa(1,3)*Fa(2,2)*Fa(3,4)+Fa(1&
!            &,4)*Fa(2,3)*Fa(3,2)- Fa(1,2)*Fa(2,3)*Fa(3,4)-Fa(1,3)&
!            &*Fa(2,4)*Fa(3,2)-Fa(1,4)*Fa(2,2)*Fa(3,3)
!        Fb(2,4)=Fa(1,1)*Fa(2,3)*Fa(3,4)+Fa(1,3)*Fa(2,4)*Fa(3,1)+Fa(1&
!            &,4)*Fa(2,1)*Fa(3,3)- Fa(1,1)*Fa(2,4)*Fa(3,3)-Fa(1,3)&
!            &*Fa(2,1)*Fa(3,4)-Fa(1,4)*Fa(2,3)*Fa(3,1)
!        Fb(3,4)=Fa(1,1)*Fa(2,4)*Fa(3,2)+Fa(1,2)*Fa(2,1)*Fa(3,4)+Fa(1&
!            &,4)*Fa(2,2)*Fa(3,1)- Fa(1,1)*Fa(2,2)*Fa(3,4)-Fa(1,2)&
!            &*Fa(2,4)*Fa(3,1)-Fa(1,4)*Fa(2,1)*Fa(3,2)
!        Fb(4,4)=Fa(1,1)*Fa(2,2)*Fa(3,3)+Fa(1,2)*Fa(2,3)*Fa(3,1)+Fa(1&
!            &,3)*Fa(2,1)*Fa(3,2)- Fa(1,1)*Fa(2,3)*Fa(3,2)-Fa(1,2)&
!            &*Fa(2,1)*Fa(3,3)-Fa(1,3)*Fa(2,2)*Fa(3,1)
!        faux=Fa(1,1)*Fa(2,2)*Fa(3,3)*Fa(4,4)+Fa(1,1)*Fa(2,3)*Fa(3,4)&
!            &*Fa(4,2)+Fa(1,1)*Fa(2,4)*Fa(3,2)*Fa(4,3)+ Fa(1,2)*Fa(2&
!            &,1)*Fa(3,4)*Fa(4,3)+Fa(1,2)*Fa(2,3)*Fa(3,1)*Fa(4,4)+Fa(1&
!            &,2)*Fa(2,4)*Fa(3,3)*Fa(4,1)+ Fa(1,3)*Fa(2,1)*Fa(3,2)&
!            &*Fa(4,4)+Fa(1,3)*Fa(2,2)*Fa(3,4)*Fa(4,1)+Fa(1,3)*Fa(2,4)&
!            &*Fa(3,1)*Fa(4,2)+ Fa(1,4)*Fa(2,1)*Fa(3,3)*Fa(4,2)+Fa(1&
!            &,4)*Fa(2,2)*Fa(3,1)*Fa(4,3)+Fa(1,4)*Fa(2,3)*Fa(3,2)*Fa(4&
!            &,1)- Fa(1,1)*Fa(2,2)*Fa(3,4)*Fa(4,3)-Fa(1,1)*Fa(2,3)&
!            &*Fa(3,2)*Fa(4,4)-Fa(1,1)*Fa(2,4)*Fa(3,3)*Fa(4,2)- Fa(1&
!            &,2)*Fa(2,1)*Fa(3,3)*Fa(4,4)-Fa(1,2)*Fa(2,3)*Fa(3,4)*Fa(4&
!            &,1)-Fa(1,2)*Fa(2,4)*Fa(3,1)*Fa(4,3)- Fa(1,3)*Fa(2,1)&
!            &*Fa(3,4)*Fa(4,2)-Fa(1,3)*Fa(2,2)*Fa(3,1)*Fa(4,4)-Fa(1,3)&
!            &*Fa(2,4)*Fa(3,2)*Fa(4,1)- Fa(1,4)*Fa(2,1)*Fa(3,2)*Fa(4&
!            &,3)-Fa(1,4)*Fa(2,2)*Fa(3,3)*Fa(4,1)-Fa(1,4)*Fa(2,3)*Fa(3&
!            &,1)*Fa(4,2)
!        Fb=Fb*(1.0_SP/faux)

        bsuccess = .true.

      else
        bsuccess = .false.
      end if

  end subroutine mprim_invert4x4MatrixDirectSP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invert5x5MatrixDirectDP(Da,Db,bsuccess)

!<description>
  ! This subroutine directly inverts a 5x5 system without any pivoting.
  ! 'Da' is a 2-dimensional 5x5 matrix. The inverse of Da is written
  ! to the 2-dimensional 5x5 matrix Db.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Da and Db are assumed to be 5x5 arrays!
!</description>

!<input>
  ! source square matrix to be inverted
  real(DP), dimension(5,5), intent(in) :: Da
!</input>

!<output>
  ! destination square matrix; receives <tex>$ A^{-1} $</tex>.
  real(DP), dimension(5,5), intent(out) :: Db

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Db is undefined.
  logical, intent(out) :: bsuccess
!</output>

!</subroutine>

  ! auxiliary variables
  real(DP) :: det,daux
  real(DP), dimension(10) :: V,W

    ! 2x2 determinants of rows 4-5
    V( 1) = Da(4,1)*Da(5,2)-Da(4,2)*Da(5,1)
    V( 2) = Da(4,1)*Da(5,3)-Da(4,3)*Da(5,1)
    V( 3) = Da(4,1)*Da(5,4)-Da(4,4)*Da(5,1)
    V( 4) = Da(4,1)*Da(5,5)-Da(4,5)*Da(5,1)
    V( 5) = Da(4,2)*Da(5,3)-Da(4,3)*Da(5,2)
    V( 6) = Da(4,2)*Da(5,4)-Da(4,4)*Da(5,2)
    V( 7) = Da(4,2)*Da(5,5)-Da(4,5)*Da(5,2)
    V( 8) = Da(4,3)*Da(5,4)-Da(4,4)*Da(5,3)
    V( 9) = Da(4,3)*Da(5,5)-Da(4,5)*Da(5,3)
    V(10) = Da(4,4)*Da(5,5)-Da(4,5)*Da(5,4)
    ! 3x3 determinants of rows 3-4-5
    W( 1) = Da(3,1)*V( 5)-Da(3,2)*V( 2)+Da(3,3)*V( 1)
    W( 2) = Da(3,1)*V( 6)-Da(3,2)*V( 3)+Da(3,4)*V( 1)
    W( 3) = Da(3,1)*V( 7)-Da(3,2)*V( 4)+Da(3,5)*V( 1)
    W( 4) = Da(3,1)*V( 8)-Da(3,3)*V( 3)+Da(3,4)*V( 2)
    W( 5) = Da(3,1)*V( 9)-Da(3,3)*V( 4)+Da(3,5)*V( 2)
    W( 6) = Da(3,1)*V(10)-Da(3,4)*V( 4)+Da(3,5)*V( 3)
    W( 7) = Da(3,2)*V( 8)-Da(3,3)*V( 6)+Da(3,4)*V( 5)
    W( 8) = Da(3,2)*V( 9)-Da(3,3)*V( 7)+Da(3,5)*V( 5)
    W( 9) = Da(3,2)*V(10)-Da(3,4)*V( 7)+Da(3,5)*V( 6)
    W(10) = Da(3,3)*V(10)-Da(3,4)*V( 9)+Da(3,5)*V( 8)
    ! pre-calculate first column of inverse
    Db(1,1) = Da(2,2)*W(10)-Da(2,3)*W( 9)+Da(2,4)*W( 8)-Da(2,5)*W( 7)
    Db(2,1) =-Da(2,1)*W(10)+Da(2,3)*W( 6)-Da(2,4)*W( 5)+Da(2,5)*W( 4)
    Db(3,1) = Da(2,1)*W( 9)-Da(2,2)*W( 6)+Da(2,4)*W( 3)-Da(2,5)*W( 2)
    Db(4,1) =-Da(2,1)*W( 8)+Da(2,2)*W( 5)-Da(2,3)*W( 3)+Da(2,5)*W( 1)
    Db(5,1) = Da(2,1)*W( 7)-Da(2,2)*W( 4)+Da(2,3)*W( 2)-Da(2,4)*W( 1)
    daux = (Da(1,1)*Db(1,1)+Da(1,2)*Db(2,1)&
            +Da(1,3)*Db(3,1)+Da(1,4)*Db(4,1)+Da(1,5)*Db(5,1))
    if (daux .ne. 0.0_DP) then
      ! calculate determinant of A
      det = 1.0_DP / daux
      ! update first column of inverse
      Db(1,1)=Db(1,1)*det
      Db(2,1)=Db(2,1)*det
      Db(3,1)=Db(3,1)*det
      Db(4,1)=Db(4,1)*det
      Db(5,1)=Db(5,1)*det
      ! calculate second column of inverse
      Db(1,2) = det*(-Da(1,2)*W(10)+Da(1,3)*W( 9)-Da(1,4)*W( 8)+Da(1,5)*W( 7))
      Db(2,2) = det*( Da(1,1)*W(10)-Da(1,3)*W( 6)+Da(1,4)*W( 5)-Da(1,5)*W( 4))
      Db(3,2) = det*(-Da(1,1)*W( 9)+Da(1,2)*W( 6)-Da(1,4)*W( 3)+Da(1,5)*W( 2))
      Db(4,2) = det*( Da(1,1)*W( 8)-Da(1,2)*W( 5)+Da(1,3)*W( 3)-Da(1,5)*W( 1))
      Db(5,2) = det*(-Da(1,1)*W( 7)+Da(1,2)*W( 4)-Da(1,3)*W( 2)+Da(1,4)*W( 1))
      ! 3x3 determinants of rows 2-4-5
      W( 1) = Da(2,1)*V( 5)-Da(2,2)*V( 2)+Da(2,3)*V( 1)
      W( 2) = Da(2,1)*V( 6)-Da(2,2)*V( 3)+Da(2,4)*V( 1)
      W( 3) = Da(2,1)*V( 7)-Da(2,2)*V( 4)+Da(2,5)*V( 1)
      W( 4) = Da(2,1)*V( 8)-Da(2,3)*V( 3)+Da(2,4)*V( 2)
      W( 5) = Da(2,1)*V( 9)-Da(2,3)*V( 4)+Da(2,5)*V( 2)
      W( 6) = Da(2,1)*V(10)-Da(2,4)*V( 4)+Da(2,5)*V( 3)
      W( 7) = Da(2,2)*V( 8)-Da(2,3)*V( 6)+Da(2,4)*V( 5)
      W( 8) = Da(2,2)*V( 9)-Da(2,3)*V( 7)+Da(2,5)*V( 5)
      W( 9) = Da(2,2)*V(10)-Da(2,4)*V( 7)+Da(2,5)*V( 6)
      W(10) = Da(2,3)*V(10)-Da(2,4)*V( 9)+Da(2,5)*V( 8)
      ! calculate third column of inverse
      Db(1,3) = det*( Da(1,2)*W(10)-Da(1,3)*W( 9)+Da(1,4)*W( 8)-Da(1,5)*W( 7))
      Db(2,3) = det*(-Da(1,1)*W(10)+Da(1,3)*W( 6)-Da(1,4)*W( 5)+Da(1,5)*W( 4))
      Db(3,3) = det*( Da(1,1)*W( 9)-Da(1,2)*W( 6)+Da(1,4)*W( 3)-Da(1,5)*W( 2))
      Db(4,3) = det*(-Da(1,1)*W( 8)+Da(1,2)*W( 5)-Da(1,3)*W( 3)+Da(1,5)*W( 1))
      Db(5,3) = det*( Da(1,1)*W( 7)-Da(1,2)*W( 4)+Da(1,3)*W( 2)-Da(1,4)*W( 1))
      ! 2x2 determinants of rows 1-2
      V( 1) = Da(1,1)*Da(2,2)-Da(1,2)*Da(2,1)
      V( 2) = Da(1,1)*Da(2,3)-Da(1,3)*Da(2,1)
      V( 3) = Da(1,1)*Da(2,4)-Da(1,4)*Da(2,1)
      V( 4) = Da(1,1)*Da(2,5)-Da(1,5)*Da(2,1)
      V( 5) = Da(1,2)*Da(2,3)-Da(1,3)*Da(2,2)
      V( 6) = Da(1,2)*Da(2,4)-Da(1,4)*Da(2,2)
      V( 7) = Da(1,2)*Da(2,5)-Da(1,5)*Da(2,2)
      V( 8) = Da(1,3)*Da(2,4)-Da(1,4)*Da(2,3)
      V( 9) = Da(1,3)*Da(2,5)-Da(1,5)*Da(2,3)
      V(10) = Da(1,4)*Da(2,5)-Da(1,5)*Da(2,4)
      ! 3x3 determinants of rows 1-2-3
      W( 1) = Da(3,1)*V( 5)-Da(3,2)*V( 2)+Da(3,3)*V( 1)
      W( 2) = Da(3,1)*V( 6)-Da(3,2)*V( 3)+Da(3,4)*V( 1)
      W( 3) = Da(3,1)*V( 7)-Da(3,2)*V( 4)+Da(3,5)*V( 1)
      W( 4) = Da(3,1)*V( 8)-Da(3,3)*V( 3)+Da(3,4)*V( 2)
      W( 5) = Da(3,1)*V( 9)-Da(3,3)*V( 4)+Da(3,5)*V( 2)
      W( 6) = Da(3,1)*V(10)-Da(3,4)*V( 4)+Da(3,5)*V( 3)
      W( 7) = Da(3,2)*V( 8)-Da(3,3)*V( 6)+Da(3,4)*V( 5)
      W( 8) = Da(3,2)*V( 9)-Da(3,3)*V( 7)+Da(3,5)*V( 5)
      W( 9) = Da(3,2)*V(10)-Da(3,4)*V( 7)+Da(3,5)*V( 6)
      W(10) = Da(3,3)*V(10)-Da(3,4)*V( 9)+Da(3,5)*V( 8)
      ! calculate fourth column of inverse
      Db(1,4) = det*( Da(5,2)*W(10)-Da(5,3)*W( 9)+Da(5,4)*W( 8)-Da(5,5)*W( 7))
      Db(2,4) = det*(-Da(5,1)*W(10)+Da(5,3)*W( 6)-Da(5,4)*W( 5)+Da(5,5)*W( 4))
      Db(3,4) = det*( Da(5,1)*W( 9)-Da(5,2)*W( 6)+Da(5,4)*W( 3)-Da(5,5)*W( 2))
      Db(4,4) = det*(-Da(5,1)*W( 8)+Da(5,2)*W( 5)-Da(5,3)*W( 3)+Da(5,5)*W( 1))
      Db(5,4) = det*( Da(5,1)*W( 7)-Da(5,2)*W( 4)+Da(5,3)*W( 2)-Da(5,4)*W( 1))
      ! calculate fifth column of inverse
      Db(1,5) = det*(-Da(4,2)*W(10)+Da(4,3)*W( 9)-Da(4,4)*W( 8)+Da(4,5)*W( 7))
      Db(2,5) = det*( Da(4,1)*W(10)-Da(4,3)*W( 6)+Da(4,4)*W( 5)-Da(4,5)*W( 4))
      Db(3,5) = det*(-Da(4,1)*W( 9)+Da(4,2)*W( 6)-Da(4,4)*W( 3)+Da(4,5)*W( 2))
      Db(4,5) = det*( Da(4,1)*W( 8)-Da(4,2)*W( 5)+Da(4,3)*W( 3)-Da(4,5)*W( 1))
      Db(5,5) = det*(-Da(4,1)*W( 7)+Da(4,2)*W( 4)-Da(4,3)*W( 2)+Da(4,4)*W( 1))

      bsuccess = .true.
    else
      bsuccess = .false.
    end if

  end subroutine mprim_invert5x5MatrixDirectDP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invert5x5MatrixDirectSP(Fa,Fb,bsuccess)

!<description>
  ! This subroutine directly inverts a 5x5 system without any pivoting.
  ! 'Fa' is a 2-dimensional 5x5 matrix. The inverse of Fa is written
  ! to the 2-dimensional 5x5 matrix Fb.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Fa and Fb are assumed to be 5x5 arrays!
!</description>

!<input>
  ! source square matrix to be inverted
  real(SP), dimension(5,5), intent(in) :: Fa
!</input>

!<output>
  ! destination square matrix; receives <tex>$ A^{-1} $</tex>.
  real(SP), dimension(5,5), intent(out) :: Fb

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Fb is undefined.
  logical, intent(out) :: bsuccess
!</output>

!</subroutine>

  ! auxiliary variables
  real(SP) :: det,faux
  real(SP), dimension(10) :: V,W

    ! 2x2 determinants of rows 4-5
    V( 1) = Fa(4,1)*Fa(5,2)-Fa(4,2)*Fa(5,1)
    V( 2) = Fa(4,1)*Fa(5,3)-Fa(4,3)*Fa(5,1)
    V( 3) = Fa(4,1)*Fa(5,4)-Fa(4,4)*Fa(5,1)
    V( 4) = Fa(4,1)*Fa(5,5)-Fa(4,5)*Fa(5,1)
    V( 5) = Fa(4,2)*Fa(5,3)-Fa(4,3)*Fa(5,2)
    V( 6) = Fa(4,2)*Fa(5,4)-Fa(4,4)*Fa(5,2)
    V( 7) = Fa(4,2)*Fa(5,5)-Fa(4,5)*Fa(5,2)
    V( 8) = Fa(4,3)*Fa(5,4)-Fa(4,4)*Fa(5,3)
    V( 9) = Fa(4,3)*Fa(5,5)-Fa(4,5)*Fa(5,3)
    V(10) = Fa(4,4)*Fa(5,5)-Fa(4,5)*Fa(5,4)
    ! 3x3 determinants of rows 3-4-5
    W( 1) = Fa(3,1)*V( 5)-Fa(3,2)*V( 2)+Fa(3,3)*V( 1)
    W( 2) = Fa(3,1)*V( 6)-Fa(3,2)*V( 3)+Fa(3,4)*V( 1)
    W( 3) = Fa(3,1)*V( 7)-Fa(3,2)*V( 4)+Fa(3,5)*V( 1)
    W( 4) = Fa(3,1)*V( 8)-Fa(3,3)*V( 3)+Fa(3,4)*V( 2)
    W( 5) = Fa(3,1)*V( 9)-Fa(3,3)*V( 4)+Fa(3,5)*V( 2)
    W( 6) = Fa(3,1)*V(10)-Fa(3,4)*V( 4)+Fa(3,5)*V( 3)
    W( 7) = Fa(3,2)*V( 8)-Fa(3,3)*V( 6)+Fa(3,4)*V( 5)
    W( 8) = Fa(3,2)*V( 9)-Fa(3,3)*V( 7)+Fa(3,5)*V( 5)
    W( 9) = Fa(3,2)*V(10)-Fa(3,4)*V( 7)+Fa(3,5)*V( 6)
    W(10) = Fa(3,3)*V(10)-Fa(3,4)*V( 9)+Fa(3,5)*V( 8)
    ! pre-calculate first column of inverse
    Fb(1,1) = Fa(2,2)*W(10)-Fa(2,3)*W( 9)+Fa(2,4)*W( 8)-Fa(2,5)*W( 7)
    Fb(2,1) =-Fa(2,1)*W(10)+Fa(2,3)*W( 6)-Fa(2,4)*W( 5)+Fa(2,5)*W( 4)
    Fb(3,1) = Fa(2,1)*W( 9)-Fa(2,2)*W( 6)+Fa(2,4)*W( 3)-Fa(2,5)*W( 2)
    Fb(4,1) =-Fa(2,1)*W( 8)+Fa(2,2)*W( 5)-Fa(2,3)*W( 3)+Fa(2,5)*W( 1)
    Fb(5,1) = Fa(2,1)*W( 7)-Fa(2,2)*W( 4)+Fa(2,3)*W( 2)-Fa(2,4)*W( 1)
    faux = (Fa(1,1)*Fb(1,1)+Fa(1,2)*Fb(2,1)&
            +Fa(1,3)*Fb(3,1)+Fa(1,4)*Fb(4,1)+Fa(1,5)*Fb(5,1))
    if (faux .ne. 0.0_SP) then
      ! calculate determinant of A
      det = 1.0_SP / faux
      ! update first column of inverse
      Fb(1,1)=Fb(1,1)*det
      Fb(2,1)=Fb(2,1)*det
      Fb(3,1)=Fb(3,1)*det
      Fb(4,1)=Fb(4,1)*det
      Fb(5,1)=Fb(5,1)*det
      ! calculate second column of inverse
      Fb(1,2) = det*(-Fa(1,2)*W(10)+Fa(1,3)*W( 9)-Fa(1,4)*W( 8)+Fa(1,5)*W( 7))
      Fb(2,2) = det*( Fa(1,1)*W(10)-Fa(1,3)*W( 6)+Fa(1,4)*W( 5)-Fa(1,5)*W( 4))
      Fb(3,2) = det*(-Fa(1,1)*W( 9)+Fa(1,2)*W( 6)-Fa(1,4)*W( 3)+Fa(1,5)*W( 2))
      Fb(4,2) = det*( Fa(1,1)*W( 8)-Fa(1,2)*W( 5)+Fa(1,3)*W( 3)-Fa(1,5)*W( 1))
      Fb(5,2) = det*(-Fa(1,1)*W( 7)+Fa(1,2)*W( 4)-Fa(1,3)*W( 2)+Fa(1,4)*W( 1))
      ! 3x3 determinants of rows 2-4-5
      W( 1) = Fa(2,1)*V( 5)-Fa(2,2)*V( 2)+Fa(2,3)*V( 1)
      W( 2) = Fa(2,1)*V( 6)-Fa(2,2)*V( 3)+Fa(2,4)*V( 1)
      W( 3) = Fa(2,1)*V( 7)-Fa(2,2)*V( 4)+Fa(2,5)*V( 1)
      W( 4) = Fa(2,1)*V( 8)-Fa(2,3)*V( 3)+Fa(2,4)*V( 2)
      W( 5) = Fa(2,1)*V( 9)-Fa(2,3)*V( 4)+Fa(2,5)*V( 2)
      W( 6) = Fa(2,1)*V(10)-Fa(2,4)*V( 4)+Fa(2,5)*V( 3)
      W( 7) = Fa(2,2)*V( 8)-Fa(2,3)*V( 6)+Fa(2,4)*V( 5)
      W( 8) = Fa(2,2)*V( 9)-Fa(2,3)*V( 7)+Fa(2,5)*V( 5)
      W( 9) = Fa(2,2)*V(10)-Fa(2,4)*V( 7)+Fa(2,5)*V( 6)
      W(10) = Fa(2,3)*V(10)-Fa(2,4)*V( 9)+Fa(2,5)*V( 8)
      ! calculate third column of inverse
      Fb(1,3) = det*( Fa(1,2)*W(10)-Fa(1,3)*W( 9)+Fa(1,4)*W( 8)-Fa(1,5)*W( 7))
      Fb(2,3) = det*(-Fa(1,1)*W(10)+Fa(1,3)*W( 6)-Fa(1,4)*W( 5)+Fa(1,5)*W( 4))
      Fb(3,3) = det*( Fa(1,1)*W( 9)-Fa(1,2)*W( 6)+Fa(1,4)*W( 3)-Fa(1,5)*W( 2))
      Fb(4,3) = det*(-Fa(1,1)*W( 8)+Fa(1,2)*W( 5)-Fa(1,3)*W( 3)+Fa(1,5)*W( 1))
      Fb(5,3) = det*( Fa(1,1)*W( 7)-Fa(1,2)*W( 4)+Fa(1,3)*W( 2)-Fa(1,4)*W( 1))
      ! 2x2 determinants of rows 1-2
      V( 1) = Fa(1,1)*Fa(2,2)-Fa(1,2)*Fa(2,1)
      V( 2) = Fa(1,1)*Fa(2,3)-Fa(1,3)*Fa(2,1)
      V( 3) = Fa(1,1)*Fa(2,4)-Fa(1,4)*Fa(2,1)
      V( 4) = Fa(1,1)*Fa(2,5)-Fa(1,5)*Fa(2,1)
      V( 5) = Fa(1,2)*Fa(2,3)-Fa(1,3)*Fa(2,2)
      V( 6) = Fa(1,2)*Fa(2,4)-Fa(1,4)*Fa(2,2)
      V( 7) = Fa(1,2)*Fa(2,5)-Fa(1,5)*Fa(2,2)
      V( 8) = Fa(1,3)*Fa(2,4)-Fa(1,4)*Fa(2,3)
      V( 9) = Fa(1,3)*Fa(2,5)-Fa(1,5)*Fa(2,3)
      V(10) = Fa(1,4)*Fa(2,5)-Fa(1,5)*Fa(2,4)
      ! 3x3 determinants of rows 1-2-3
      W( 1) = Fa(3,1)*V( 5)-Fa(3,2)*V( 2)+Fa(3,3)*V( 1)
      W( 2) = Fa(3,1)*V( 6)-Fa(3,2)*V( 3)+Fa(3,4)*V( 1)
      W( 3) = Fa(3,1)*V( 7)-Fa(3,2)*V( 4)+Fa(3,5)*V( 1)
      W( 4) = Fa(3,1)*V( 8)-Fa(3,3)*V( 3)+Fa(3,4)*V( 2)
      W( 5) = Fa(3,1)*V( 9)-Fa(3,3)*V( 4)+Fa(3,5)*V( 2)
      W( 6) = Fa(3,1)*V(10)-Fa(3,4)*V( 4)+Fa(3,5)*V( 3)
      W( 7) = Fa(3,2)*V( 8)-Fa(3,3)*V( 6)+Fa(3,4)*V( 5)
      W( 8) = Fa(3,2)*V( 9)-Fa(3,3)*V( 7)+Fa(3,5)*V( 5)
      W( 9) = Fa(3,2)*V(10)-Fa(3,4)*V( 7)+Fa(3,5)*V( 6)
      W(10) = Fa(3,3)*V(10)-Fa(3,4)*V( 9)+Fa(3,5)*V( 8)
      ! calculate fourth column of inverse
      Fb(1,4) = det*( Fa(5,2)*W(10)-Fa(5,3)*W( 9)+Fa(5,4)*W( 8)-Fa(5,5)*W( 7))
      Fb(2,4) = det*(-Fa(5,1)*W(10)+Fa(5,3)*W( 6)-Fa(5,4)*W( 5)+Fa(5,5)*W( 4))
      Fb(3,4) = det*( Fa(5,1)*W( 9)-Fa(5,2)*W( 6)+Fa(5,4)*W( 3)-Fa(5,5)*W( 2))
      Fb(4,4) = det*(-Fa(5,1)*W( 8)+Fa(5,2)*W( 5)-Fa(5,3)*W( 3)+Fa(5,5)*W( 1))
      Fb(5,4) = det*( Fa(5,1)*W( 7)-Fa(5,2)*W( 4)+Fa(5,3)*W( 2)-Fa(5,4)*W( 1))
      ! calculate fifth column of inverse
      Fb(1,5) = det*(-Fa(4,2)*W(10)+Fa(4,3)*W( 9)-Fa(4,4)*W( 8)+Fa(4,5)*W( 7))
      Fb(2,5) = det*( Fa(4,1)*W(10)-Fa(4,3)*W( 6)+Fa(4,4)*W( 5)-Fa(4,5)*W( 4))
      Fb(3,5) = det*(-Fa(4,1)*W( 9)+Fa(4,2)*W( 6)-Fa(4,4)*W( 3)+Fa(4,5)*W( 2))
      Fb(4,5) = det*( Fa(4,1)*W( 8)-Fa(4,2)*W( 5)+Fa(4,3)*W( 3)-Fa(4,5)*W( 1))
      Fb(5,5) = det*(-Fa(4,1)*W( 7)+Fa(4,2)*W( 4)-Fa(4,3)*W( 2)+Fa(4,4)*W( 1))

      bsuccess = .true.
    else
      bsuccess = .false.
    end if

  end subroutine mprim_invert5x5MatrixDirectSP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invert6x6MatrixDirectDP(Da,Db,bsuccess)

!<description>
  ! This subroutine directly inverts a 6x6 system without any pivoting.
  ! 'Da' is a 2-dimensional 6x6 m matrix. The inverse of Da is written
  ! to the 2-dimensional 6x6 matrix Db.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Da and Db are assumed to be 6x6 arrays!
!</description>

!<input>
  ! source square matrix to be inverted
  real(DP), dimension(6,6), intent(in) :: Da
!</input>

!<output>
  ! destination square matrix; receives <tex>$ A^{-1} $</tex>.
  real(DP), dimension(6,6), intent(out) :: Db

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Db is undefined.
  logical, intent(out) :: bsuccess
!</output>

!</subroutine>

  ! auxiliary variables
  real(DP) :: det,daux
  real(DP), dimension(15) :: U, W
  real(DP), dimension(20) :: V

    ! 2x2 determinants of rows 5-6
    U(1) = Da(5,1)*Da(6,2)-Da(5,2)*Da(6,1)
    U(2) = Da(5,1)*Da(6,3)-Da(5,3)*Da(6,1)
    U(3) = Da(5,1)*Da(6,4)-Da(5,4)*Da(6,1)
    U(4) = Da(5,1)*Da(6,5)-Da(5,5)*Da(6,1)
    U(5) = Da(5,1)*Da(6,6)-Da(5,6)*Da(6,1)
    U(6) = Da(5,2)*Da(6,3)-Da(5,3)*Da(6,2)
    U(7) = Da(5,2)*Da(6,4)-Da(5,4)*Da(6,2)
    U(8) = Da(5,2)*Da(6,5)-Da(5,5)*Da(6,2)
    U(9) = Da(5,2)*Da(6,6)-Da(5,6)*Da(6,2)
    U(10) = Da(5,3)*Da(6,4)-Da(5,4)*Da(6,3)
    U(11) = Da(5,3)*Da(6,5)-Da(5,5)*Da(6,3)
    U(12) = Da(5,3)*Da(6,6)-Da(5,6)*Da(6,3)
    U(13) = Da(5,4)*Da(6,5)-Da(5,5)*Da(6,4)
    U(14) = Da(5,4)*Da(6,6)-Da(5,6)*Da(6,4)
    U(15) = Da(5,5)*Da(6,6)-Da(5,6)*Da(6,5)
    ! 3x3 determinants of rows 4-5-6
    V(1) = Da(4,1)*U(6)-Da(4,2)*U(2)+Da(4,3)*U(1)
    V(2) = Da(4,1)*U(7)-Da(4,2)*U(3)+Da(4,4)*U(1)
    V(3) = Da(4,1)*U(8)-Da(4,2)*U(4)+Da(4,5)*U(1)
    V(4) = Da(4,1)*U(9)-Da(4,2)*U(5)+Da(4,6)*U(1)
    V(5) = Da(4,1)*U(10)-Da(4,3)*U(3)+Da(4,4)*U(2)
    V(6) = Da(4,1)*U(11)-Da(4,3)*U(4)+Da(4,5)*U(2)
    V(7) = Da(4,1)*U(12)-Da(4,3)*U(5)+Da(4,6)*U(2)
    V(8) = Da(4,1)*U(13)-Da(4,4)*U(4)+Da(4,5)*U(3)
    V(9) = Da(4,1)*U(14)-Da(4,4)*U(5)+Da(4,6)*U(3)
    V(10) = Da(4,1)*U(15)-Da(4,5)*U(5)+Da(4,6)*U(4)
    V(11) = Da(4,2)*U(10)-Da(4,3)*U(7)+Da(4,4)*U(6)
    V(12) = Da(4,2)*U(11)-Da(4,3)*U(8)+Da(4,5)*U(6)
    V(13) = Da(4,2)*U(12)-Da(4,3)*U(9)+Da(4,6)*U(6)
    V(14) = Da(4,2)*U(13)-Da(4,4)*U(8)+Da(4,5)*U(7)
    V(15) = Da(4,2)*U(14)-Da(4,4)*U(9)+Da(4,6)*U(7)
    V(16) = Da(4,2)*U(15)-Da(4,5)*U(9)+Da(4,6)*U(8)
    V(17) = Da(4,3)*U(13)-Da(4,4)*U(11)+Da(4,5)*U(10)
    V(18) = Da(4,3)*U(14)-Da(4,4)*U(12)+Da(4,6)*U(10)
    V(19) = Da(4,3)*U(15)-Da(4,5)*U(12)+Da(4,6)*U(11)
    V(20) = Da(4,4)*U(15)-Da(4,5)*U(14)+Da(4,6)*U(13)
    ! 4x4 determinants of rows 3-4-5-6
    W(1) = Da(3,1)*V(11)-Da(3,2)*V(5)+Da(3,3)*V(2)-Da(3,4)*V(1)
    W(2) = Da(3,1)*V(12)-Da(3,2)*V(6)+Da(3,3)*V(3)-Da(3,5)*V(1)
    W(3) = Da(3,1)*V(13)-Da(3,2)*V(7)+Da(3,3)*V(4)-Da(3,6)*V(1)
    W(4) = Da(3,1)*V(14)-Da(3,2)*V(8)+Da(3,4)*V(3)-Da(3,5)*V(2)
    W(5) = Da(3,1)*V(15)-Da(3,2)*V(9)+Da(3,4)*V(4)-Da(3,6)*V(2)
    W(6) = Da(3,1)*V(16)-Da(3,2)*V(10)+Da(3,5)*V(4)-Da(3,6)*V(3)
    W(7) = Da(3,1)*V(17)-Da(3,3)*V(8)+Da(3,4)*V(6)-Da(3,5)*V(5)
    W(8) = Da(3,1)*V(18)-Da(3,3)*V(9)+Da(3,4)*V(7)-Da(3,6)*V(5)
    W(9) = Da(3,1)*V(19)-Da(3,3)*V(10)+Da(3,5)*V(7)-Da(3,6)*V(6)
    W(10) = Da(3,1)*V(20)-Da(3,4)*V(10)+Da(3,5)*V(9)-Da(3,6)*V(8)
    W(11) = Da(3,2)*V(17)-Da(3,3)*V(14)+Da(3,4)*V(12)-Da(3,5)*V(11)
    W(12) = Da(3,2)*V(18)-Da(3,3)*V(15)+Da(3,4)*V(13)-Da(3,6)*V(11)
    W(13) = Da(3,2)*V(19)-Da(3,3)*V(16)+Da(3,5)*V(13)-Da(3,6)*V(12)
    W(14) = Da(3,2)*V(20)-Da(3,4)*V(16)+Da(3,5)*V(15)-Da(3,6)*V(14)
    W(15) = Da(3,3)*V(20)-Da(3,4)*V(19)+Da(3,5)*V(18)-Da(3,6)*V(17)
    ! pre-calculate first column of inverse
    Db(1,1) = Da(2,2)*W(15)-Da(2,3)*W(14)+Da(2,4)*W(13)-Da(2,5)*W(12)+Da(2,6)*W(11)
    Db(2,1) =-Da(2,1)*W(15)+Da(2,3)*W(10)-Da(2,4)*W(9)+Da(2,5)*W(8)-Da(2,6)*W(7)
    Db(3,1) = Da(2,1)*W(14)-Da(2,2)*W(10)+Da(2,4)*W(6)-Da(2,5)*W(5)+Da(2,6)*W(4)
    Db(4,1) =-Da(2,1)*W(13)+Da(2,2)*W(9)-Da(2,3)*W(6)+Da(2,5)*W(3)-Da(2,6)*W(2)
    Db(5,1) = Da(2,1)*W(12)-Da(2,2)*W(8)+Da(2,3)*W(5)-Da(2,4)*W(3)+Da(2,6)*W(1)
    Db(6,1) =-Da(2,1)*W(11)+Da(2,2)*W(7)-Da(2,3)*W(4)+Da(2,4)*W(2)-Da(2,5)*W(1)

    daux = (Da(1,1)*Db(1,1)+Da(1,2)*Db(2,1)+Da(1,3)*Db(3,1)+&
            Da(1,4)*Db(4,1)+Da(1,5)*Db(5,1)+Da(1,6)*Db(6,1))

    if (daux .ne. 0.0_DP) then

      ! calculate determinant of A
      det = 1.0_DP / daux
      ! update first column of inverse
      Db(1,1) = det*Db(1,1)
      Db(2,1) = det*Db(2,1)
      Db(3,1) = det*Db(3,1)
      Db(4,1) = det*Db(4,1)
      Db(5,1) = det*Db(5,1)
      Db(6,1) = det*Db(6,1)
      ! calculate second column of inverse
      Db(1,2) = det*(-Da(1,2)*W(15)+Da(1,3)*W(14)-Da(1,4)*W(13)+Da(1,5)*W(12)-Da(1,6)*W(11))
      Db(2,2) = det*( Da(1,1)*W(15)-Da(1,3)*W(10)+Da(1,4)*W(9)-Da(1,5)*W(8)+Da(1,6)*W(7))
      Db(3,2) = det*(-Da(1,1)*W(14)+Da(1,2)*W(10)-Da(1,4)*W(6)+Da(1,5)*W(5)-Da(1,6)*W(4))
      Db(4,2) = det*( Da(1,1)*W(13)-Da(1,2)*W(9)+Da(1,3)*W(6)-Da(1,5)*W(3)+Da(1,6)*W(2))
      Db(5,2) = det*(-Da(1,1)*W(12)+Da(1,2)*W(8)-Da(1,3)*W(5)+Da(1,4)*W(3)-Da(1,6)*W(1))
      Db(6,2) = det*( Da(1,1)*W(11)-Da(1,2)*W(7)+Da(1,3)*W(4)-Da(1,4)*W(2)+Da(1,5)*W(1))
      ! 3x3 determinants of rows 2-5-6
      V(1) = Da(2,1)*U(6)-Da(2,2)*U(2)+Da(2,3)*U(1)
      V(2) = Da(2,1)*U(7)-Da(2,2)*U(3)+Da(2,4)*U(1)
      V(3) = Da(2,1)*U(8)-Da(2,2)*U(4)+Da(2,5)*U(1)
      V(4) = Da(2,1)*U(9)-Da(2,2)*U(5)+Da(2,6)*U(1)
      V(5) = Da(2,1)*U(10)-Da(2,3)*U(3)+Da(2,4)*U(2)
      V(6) = Da(2,1)*U(11)-Da(2,3)*U(4)+Da(2,5)*U(2)
      V(7) = Da(2,1)*U(12)-Da(2,3)*U(5)+Da(2,6)*U(2)
      V(8) = Da(2,1)*U(13)-Da(2,4)*U(4)+Da(2,5)*U(3)
      V(9) = Da(2,1)*U(14)-Da(2,4)*U(5)+Da(2,6)*U(3)
      V(10) = Da(2,1)*U(15)-Da(2,5)*U(5)+Da(2,6)*U(4)
      V(11) = Da(2,2)*U(10)-Da(2,3)*U(7)+Da(2,4)*U(6)
      V(12) = Da(2,2)*U(11)-Da(2,3)*U(8)+Da(2,5)*U(6)
      V(13) = Da(2,2)*U(12)-Da(2,3)*U(9)+Da(2,6)*U(6)
      V(14) = Da(2,2)*U(13)-Da(2,4)*U(8)+Da(2,5)*U(7)
      V(15) = Da(2,2)*U(14)-Da(2,4)*U(9)+Da(2,6)*U(7)
      V(16) = Da(2,2)*U(15)-Da(2,5)*U(9)+Da(2,6)*U(8)
      V(17) = Da(2,3)*U(13)-Da(2,4)*U(11)+Da(2,5)*U(10)
      V(18) = Da(2,3)*U(14)-Da(2,4)*U(12)+Da(2,6)*U(10)
      V(19) = Da(2,3)*U(15)-Da(2,5)*U(12)+Da(2,6)*U(11)
      V(20) = Da(2,4)*U(15)-Da(2,5)*U(14)+Da(2,6)*U(13)
      ! 4x4 determinants of rows 1-2-5-6
      W(1) = Da(1,1)*V(11)-Da(1,2)*V(5)+Da(1,3)*V(2)-Da(1,4)*V(1)
      W(2) = Da(1,1)*V(12)-Da(1,2)*V(6)+Da(1,3)*V(3)-Da(1,5)*V(1)
      W(3) = Da(1,1)*V(13)-Da(1,2)*V(7)+Da(1,3)*V(4)-Da(1,6)*V(1)
      W(4) = Da(1,1)*V(14)-Da(1,2)*V(8)+Da(1,4)*V(3)-Da(1,5)*V(2)
      W(5) = Da(1,1)*V(15)-Da(1,2)*V(9)+Da(1,4)*V(4)-Da(1,6)*V(2)
      W(6) = Da(1,1)*V(16)-Da(1,2)*V(10)+Da(1,5)*V(4)-Da(1,6)*V(3)
      W(7) = Da(1,1)*V(17)-Da(1,3)*V(8)+Da(1,4)*V(6)-Da(1,5)*V(5)
      W(8) = Da(1,1)*V(18)-Da(1,3)*V(9)+Da(1,4)*V(7)-Da(1,6)*V(5)
      W(9) = Da(1,1)*V(19)-Da(1,3)*V(10)+Da(1,5)*V(7)-Da(1,6)*V(6)
      W(10) = Da(1,1)*V(20)-Da(1,4)*V(10)+Da(1,5)*V(9)-Da(1,6)*V(8)
      W(11) = Da(1,2)*V(17)-Da(1,3)*V(14)+Da(1,4)*V(12)-Da(1,5)*V(11)
      W(12) = Da(1,2)*V(18)-Da(1,3)*V(15)+Da(1,4)*V(13)-Da(1,6)*V(11)
      W(13) = Da(1,2)*V(19)-Da(1,3)*V(16)+Da(1,5)*V(13)-Da(1,6)*V(12)
      W(14) = Da(1,2)*V(20)-Da(1,4)*V(16)+Da(1,5)*V(15)-Da(1,6)*V(14)
      W(15) = Da(1,3)*V(20)-Da(1,4)*V(19)+Da(1,5)*V(18)-Da(1,6)*V(17)
      ! calculate third column of inverse
      Db(1,3) = det*( Da(4,2)*W(15)-Da(4,3)*W(14)+Da(4,4)*W(13)-Da(4,5)*W(12)+Da(4,6)*W(11))
      Db(2,3) = det*(-Da(4,1)*W(15)+Da(4,3)*W(10)-Da(4,4)*W(9)+Da(4,5)*W(8)-Da(4,6)*W(7))
      Db(3,3) = det*( Da(4,1)*W(14)-Da(4,2)*W(10)+Da(4,4)*W(6)-Da(4,5)*W(5)+Da(4,6)*W(4))
      Db(4,3) = det*(-Da(4,1)*W(13)+Da(4,2)*W(9)-Da(4,3)*W(6)+Da(4,5)*W(3)-Da(4,6)*W(2))
      Db(5,3) = det*( Da(4,1)*W(12)-Da(4,2)*W(8)+Da(4,3)*W(5)-Da(4,4)*W(3)+Da(4,6)*W(1))
      Db(6,3) = det*(-Da(4,1)*W(11)+Da(4,2)*W(7)-Da(4,3)*W(4)+Da(4,4)*W(2)-Da(4,5)*W(1))
      ! calculate fourth column of inverse
      Db(1,4) = det*(-Da(3,2)*W(15)+Da(3,3)*W(14)-Da(3,4)*W(13)+Da(3,5)*W(12)-Da(3,6)*W(11))
      Db(2,4) = det*( Da(3,1)*W(15)-Da(3,3)*W(10)+Da(3,4)*W(9)-Da(3,5)*W(8)+Da(3,6)*W(7))
      Db(3,4) = det*(-Da(3,1)*W(14)+Da(3,2)*W(10)-Da(3,4)*W(6)+Da(3,5)*W(5)-Da(3,6)*W(4))
      Db(4,4) = det*( Da(3,1)*W(13)-Da(3,2)*W(9)+Da(3,3)*W(6)-Da(3,5)*W(3)+Da(3,6)*W(2))
      Db(5,4) = det*(-Da(3,1)*W(12)+Da(3,2)*W(8)-Da(3,3)*W(5)+Da(3,4)*W(3)-Da(3,6)*W(1))
      Db(6,4) = det*( Da(3,1)*W(11)-Da(3,2)*W(7)+Da(3,3)*W(4)-Da(3,4)*W(2)+Da(3,5)*W(1))
      ! 2x2 determinants of rows 3-4
      U(1) = Da(3,1)*Da(4,2)-Da(3,2)*Da(4,1)
      U(2) = Da(3,1)*Da(4,3)-Da(3,3)*Da(4,1)
      U(3) = Da(3,1)*Da(4,4)-Da(3,4)*Da(4,1)
      U(4) = Da(3,1)*Da(4,5)-Da(3,5)*Da(4,1)
      U(5) = Da(3,1)*Da(4,6)-Da(3,6)*Da(4,1)
      U(6) = Da(3,2)*Da(4,3)-Da(3,3)*Da(4,2)
      U(7) = Da(3,2)*Da(4,4)-Da(3,4)*Da(4,2)
      U(8) = Da(3,2)*Da(4,5)-Da(3,5)*Da(4,2)
      U(9) = Da(3,2)*Da(4,6)-Da(3,6)*Da(4,2)
      U(10) = Da(3,3)*Da(4,4)-Da(3,4)*Da(4,3)
      U(11) = Da(3,3)*Da(4,5)-Da(3,5)*Da(4,3)
      U(12) = Da(3,3)*Da(4,6)-Da(3,6)*Da(4,3)
      U(13) = Da(3,4)*Da(4,5)-Da(3,5)*Da(4,4)
      U(14) = Da(3,4)*Da(4,6)-Da(3,6)*Da(4,4)
      U(15) = Da(3,5)*Da(4,6)-Da(3,6)*Da(4,5)
      ! 3x3 determinants of rows 2-3-4
      V(1) = Da(2,1)*U(6)-Da(2,2)*U(2)+Da(2,3)*U(1)
      V(2) = Da(2,1)*U(7)-Da(2,2)*U(3)+Da(2,4)*U(1)
      V(3) = Da(2,1)*U(8)-Da(2,2)*U(4)+Da(2,5)*U(1)
      V(4) = Da(2,1)*U(9)-Da(2,2)*U(5)+Da(2,6)*U(1)
      V(5) = Da(2,1)*U(10)-Da(2,3)*U(3)+Da(2,4)*U(2)
      V(6) = Da(2,1)*U(11)-Da(2,3)*U(4)+Da(2,5)*U(2)
      V(7) = Da(2,1)*U(12)-Da(2,3)*U(5)+Da(2,6)*U(2)
      V(8) = Da(2,1)*U(13)-Da(2,4)*U(4)+Da(2,5)*U(3)
      V(9) = Da(2,1)*U(14)-Da(2,4)*U(5)+Da(2,6)*U(3)
      V(10) = Da(2,1)*U(15)-Da(2,5)*U(5)+Da(2,6)*U(4)
      V(11) = Da(2,2)*U(10)-Da(2,3)*U(7)+Da(2,4)*U(6)
      V(12) = Da(2,2)*U(11)-Da(2,3)*U(8)+Da(2,5)*U(6)
      V(13) = Da(2,2)*U(12)-Da(2,3)*U(9)+Da(2,6)*U(6)
      V(14) = Da(2,2)*U(13)-Da(2,4)*U(8)+Da(2,5)*U(7)
      V(15) = Da(2,2)*U(14)-Da(2,4)*U(9)+Da(2,6)*U(7)
      V(16) = Da(2,2)*U(15)-Da(2,5)*U(9)+Da(2,6)*U(8)
      V(17) = Da(2,3)*U(13)-Da(2,4)*U(11)+Da(2,5)*U(10)
      V(18) = Da(2,3)*U(14)-Da(2,4)*U(12)+Da(2,6)*U(10)
      V(19) = Da(2,3)*U(15)-Da(2,5)*U(12)+Da(2,6)*U(11)
      V(20) = Da(2,4)*U(15)-Da(2,5)*U(14)+Da(2,6)*U(13)
      ! 4x4 determinants of rows 1-2-3-4
      W(1) = Da(1,1)*V(11)-Da(1,2)*V(5)+Da(1,3)*V(2)-Da(1,4)*V(1)
      W(2) = Da(1,1)*V(12)-Da(1,2)*V(6)+Da(1,3)*V(3)-Da(1,5)*V(1)
      W(3) = Da(1,1)*V(13)-Da(1,2)*V(7)+Da(1,3)*V(4)-Da(1,6)*V(1)
      W(4) = Da(1,1)*V(14)-Da(1,2)*V(8)+Da(1,4)*V(3)-Da(1,5)*V(2)
      W(5) = Da(1,1)*V(15)-Da(1,2)*V(9)+Da(1,4)*V(4)-Da(1,6)*V(2)
      W(6) = Da(1,1)*V(16)-Da(1,2)*V(10)+Da(1,5)*V(4)-Da(1,6)*V(3)
      W(7) = Da(1,1)*V(17)-Da(1,3)*V(8)+Da(1,4)*V(6)-Da(1,5)*V(5)
      W(8) = Da(1,1)*V(18)-Da(1,3)*V(9)+Da(1,4)*V(7)-Da(1,6)*V(5)
      W(9) = Da(1,1)*V(19)-Da(1,3)*V(10)+Da(1,5)*V(7)-Da(1,6)*V(6)
      W(10) = Da(1,1)*V(20)-Da(1,4)*V(10)+Da(1,5)*V(9)-Da(1,6)*V(8)
      W(11) = Da(1,2)*V(17)-Da(1,3)*V(14)+Da(1,4)*V(12)-Da(1,5)*V(11)
      W(12) = Da(1,2)*V(18)-Da(1,3)*V(15)+Da(1,4)*V(13)-Da(1,6)*V(11)
      W(13) = Da(1,2)*V(19)-Da(1,3)*V(16)+Da(1,5)*V(13)-Da(1,6)*V(12)
      W(14) = Da(1,2)*V(20)-Da(1,4)*V(16)+Da(1,5)*V(15)-Da(1,6)*V(14)
      W(15) = Da(1,3)*V(20)-Da(1,4)*V(19)+Da(1,5)*V(18)-Da(1,6)*V(17)
      ! calculate fifth column of inverse
      Db(1,5) = det*( Da(6,2)*W(15)-Da(6,3)*W(14)+Da(6,4)*W(13)-Da(6,5)*W(12)+Da(6,6)*W(11))
      Db(2,5) = det*(-Da(6,1)*W(15)+Da(6,3)*W(10)-Da(6,4)*W(9)+Da(6,5)*W(8)-Da(6,6)*W(7))
      Db(3,5) = det*( Da(6,1)*W(14)-Da(6,2)*W(10)+Da(6,4)*W(6)-Da(6,5)*W(5)+Da(6,6)*W(4))
      Db(4,5) = det*(-Da(6,1)*W(13)+Da(6,2)*W(9)-Da(6,3)*W(6)+Da(6,5)*W(3)-Da(6,6)*W(2))
      Db(5,5) = det*( Da(6,1)*W(12)-Da(6,2)*W(8)+Da(6,3)*W(5)-Da(6,4)*W(3)+Da(6,6)*W(1))
      Db(6,5) = det*(-Da(6,1)*W(11)+Da(6,2)*W(7)-Da(6,3)*W(4)+Da(6,4)*W(2)-Da(6,5)*W(1))
      ! calculate sixth column of inverse
      Db(1,6) = det*(-Da(5,2)*W(15)+Da(5,3)*W(14)-Da(5,4)*W(13)+Da(5,5)*W(12)-Da(5,6)*W(11))
      Db(2,6) = det*( Da(5,1)*W(15)-Da(5,3)*W(10)+Da(5,4)*W(9)-Da(5,5)*W(8)+Da(5,6)*W(7))
      Db(3,6) = det*(-Da(5,1)*W(14)+Da(5,2)*W(10)-Da(5,4)*W(6)+Da(5,5)*W(5)-Da(5,6)*W(4))
      Db(4,6) = det*( Da(5,1)*W(13)-Da(5,2)*W(9)+Da(5,3)*W(6)-Da(5,5)*W(3)+Da(5,6)*W(2))
      Db(5,6) = det*(-Da(5,1)*W(12)+Da(5,2)*W(8)-Da(5,3)*W(5)+Da(5,4)*W(3)-Da(5,6)*W(1))
      Db(6,6) = det*( Da(5,1)*W(11)-Da(5,2)*W(7)+Da(5,3)*W(4)-Da(5,4)*W(2)+Da(5,5)*W(1))

      bsuccess = .true.

    else

      bsuccess = .false.

    end if

  end subroutine mprim_invert6x6MatrixDirectDP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invert6x6MatrixDirectSP(Fa,Fb,bsuccess)

!<description>
  ! This subroutine directly inverts a 6x6 system without any pivoting.
  ! 'Fa' is a 2-dimensional 6x6 m matrix. The inverse of Fa is written
  ! to the 2-dimensional 6x6 matrix Fb.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Fa and Fb are assumed to be 6x6 arrays!
!</description>

!<input>
  ! source square matrix to be inverted
  real(SP), dimension(6,6), intent(in) :: Fa
!</input>

!<output>
  ! destination square matrix; receives <tex>$ A^{-1} $</tex>.
  real(SP), dimension(6,6), intent(out) :: Fb

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Fb is undefined.
  logical, intent(out) :: bsuccess
!</output>

!</subroutine>

  ! auxiliary variables
  real(SP) :: det,faux
  real(SP), dimension(15) :: U, W
  real(SP), dimension(20) :: V

    ! 2x2 determinants of rows 5-6
    U(1) = Fa(5,1)*Fa(6,2)-Fa(5,2)*Fa(6,1)
    U(2) = Fa(5,1)*Fa(6,3)-Fa(5,3)*Fa(6,1)
    U(3) = Fa(5,1)*Fa(6,4)-Fa(5,4)*Fa(6,1)
    U(4) = Fa(5,1)*Fa(6,5)-Fa(5,5)*Fa(6,1)
    U(5) = Fa(5,1)*Fa(6,6)-Fa(5,6)*Fa(6,1)
    U(6) = Fa(5,2)*Fa(6,3)-Fa(5,3)*Fa(6,2)
    U(7) = Fa(5,2)*Fa(6,4)-Fa(5,4)*Fa(6,2)
    U(8) = Fa(5,2)*Fa(6,5)-Fa(5,5)*Fa(6,2)
    U(9) = Fa(5,2)*Fa(6,6)-Fa(5,6)*Fa(6,2)
    U(10) = Fa(5,3)*Fa(6,4)-Fa(5,4)*Fa(6,3)
    U(11) = Fa(5,3)*Fa(6,5)-Fa(5,5)*Fa(6,3)
    U(12) = Fa(5,3)*Fa(6,6)-Fa(5,6)*Fa(6,3)
    U(13) = Fa(5,4)*Fa(6,5)-Fa(5,5)*Fa(6,4)
    U(14) = Fa(5,4)*Fa(6,6)-Fa(5,6)*Fa(6,4)
    U(15) = Fa(5,5)*Fa(6,6)-Fa(5,6)*Fa(6,5)
    ! 3x3 determinants of rows 4-5-6
    V(1) = Fa(4,1)*U(6)-Fa(4,2)*U(2)+Fa(4,3)*U(1)
    V(2) = Fa(4,1)*U(7)-Fa(4,2)*U(3)+Fa(4,4)*U(1)
    V(3) = Fa(4,1)*U(8)-Fa(4,2)*U(4)+Fa(4,5)*U(1)
    V(4) = Fa(4,1)*U(9)-Fa(4,2)*U(5)+Fa(4,6)*U(1)
    V(5) = Fa(4,1)*U(10)-Fa(4,3)*U(3)+Fa(4,4)*U(2)
    V(6) = Fa(4,1)*U(11)-Fa(4,3)*U(4)+Fa(4,5)*U(2)
    V(7) = Fa(4,1)*U(12)-Fa(4,3)*U(5)+Fa(4,6)*U(2)
    V(8) = Fa(4,1)*U(13)-Fa(4,4)*U(4)+Fa(4,5)*U(3)
    V(9) = Fa(4,1)*U(14)-Fa(4,4)*U(5)+Fa(4,6)*U(3)
    V(10) = Fa(4,1)*U(15)-Fa(4,5)*U(5)+Fa(4,6)*U(4)
    V(11) = Fa(4,2)*U(10)-Fa(4,3)*U(7)+Fa(4,4)*U(6)
    V(12) = Fa(4,2)*U(11)-Fa(4,3)*U(8)+Fa(4,5)*U(6)
    V(13) = Fa(4,2)*U(12)-Fa(4,3)*U(9)+Fa(4,6)*U(6)
    V(14) = Fa(4,2)*U(13)-Fa(4,4)*U(8)+Fa(4,5)*U(7)
    V(15) = Fa(4,2)*U(14)-Fa(4,4)*U(9)+Fa(4,6)*U(7)
    V(16) = Fa(4,2)*U(15)-Fa(4,5)*U(9)+Fa(4,6)*U(8)
    V(17) = Fa(4,3)*U(13)-Fa(4,4)*U(11)+Fa(4,5)*U(10)
    V(18) = Fa(4,3)*U(14)-Fa(4,4)*U(12)+Fa(4,6)*U(10)
    V(19) = Fa(4,3)*U(15)-Fa(4,5)*U(12)+Fa(4,6)*U(11)
    V(20) = Fa(4,4)*U(15)-Fa(4,5)*U(14)+Fa(4,6)*U(13)
    ! 4x4 determinants of rows 3-4-5-6
    W(1) = Fa(3,1)*V(11)-Fa(3,2)*V(5)+Fa(3,3)*V(2)-Fa(3,4)*V(1)
    W(2) = Fa(3,1)*V(12)-Fa(3,2)*V(6)+Fa(3,3)*V(3)-Fa(3,5)*V(1)
    W(3) = Fa(3,1)*V(13)-Fa(3,2)*V(7)+Fa(3,3)*V(4)-Fa(3,6)*V(1)
    W(4) = Fa(3,1)*V(14)-Fa(3,2)*V(8)+Fa(3,4)*V(3)-Fa(3,5)*V(2)
    W(5) = Fa(3,1)*V(15)-Fa(3,2)*V(9)+Fa(3,4)*V(4)-Fa(3,6)*V(2)
    W(6) = Fa(3,1)*V(16)-Fa(3,2)*V(10)+Fa(3,5)*V(4)-Fa(3,6)*V(3)
    W(7) = Fa(3,1)*V(17)-Fa(3,3)*V(8)+Fa(3,4)*V(6)-Fa(3,5)*V(5)
    W(8) = Fa(3,1)*V(18)-Fa(3,3)*V(9)+Fa(3,4)*V(7)-Fa(3,6)*V(5)
    W(9) = Fa(3,1)*V(19)-Fa(3,3)*V(10)+Fa(3,5)*V(7)-Fa(3,6)*V(6)
    W(10) = Fa(3,1)*V(20)-Fa(3,4)*V(10)+Fa(3,5)*V(9)-Fa(3,6)*V(8)
    W(11) = Fa(3,2)*V(17)-Fa(3,3)*V(14)+Fa(3,4)*V(12)-Fa(3,5)*V(11)
    W(12) = Fa(3,2)*V(18)-Fa(3,3)*V(15)+Fa(3,4)*V(13)-Fa(3,6)*V(11)
    W(13) = Fa(3,2)*V(19)-Fa(3,3)*V(16)+Fa(3,5)*V(13)-Fa(3,6)*V(12)
    W(14) = Fa(3,2)*V(20)-Fa(3,4)*V(16)+Fa(3,5)*V(15)-Fa(3,6)*V(14)
    W(15) = Fa(3,3)*V(20)-Fa(3,4)*V(19)+Fa(3,5)*V(18)-Fa(3,6)*V(17)
    ! pre-calculate first column of inverse
    Fb(1,1) = Fa(2,2)*W(15)-Fa(2,3)*W(14)+Fa(2,4)*W(13)-Fa(2,5)*W(12)+Fa(2,6)*W(11)
    Fb(2,1) =-Fa(2,1)*W(15)+Fa(2,3)*W(10)-Fa(2,4)*W(9)+Fa(2,5)*W(8)-Fa(2,6)*W(7)
    Fb(3,1) = Fa(2,1)*W(14)-Fa(2,2)*W(10)+Fa(2,4)*W(6)-Fa(2,5)*W(5)+Fa(2,6)*W(4)
    Fb(4,1) =-Fa(2,1)*W(13)+Fa(2,2)*W(9)-Fa(2,3)*W(6)+Fa(2,5)*W(3)-Fa(2,6)*W(2)
    Fb(5,1) = Fa(2,1)*W(12)-Fa(2,2)*W(8)+Fa(2,3)*W(5)-Fa(2,4)*W(3)+Fa(2,6)*W(1)
    Fb(6,1) =-Fa(2,1)*W(11)+Fa(2,2)*W(7)-Fa(2,3)*W(4)+Fa(2,4)*W(2)-Fa(2,5)*W(1)

    faux = (Fa(1,1)*Fb(1,1)+Fa(1,2)*Fb(2,1)+Fa(1,3)*Fb(3,1)+&
            Fa(1,4)*Fb(4,1)+Fa(1,5)*Fb(5,1)+Fa(1,6)*Fb(6,1))

    if (faux .ne. 0.0_SP) then

      ! calculate determinant of A
      det = 1.0_SP / faux
      ! update first column of inverse
      Fb(1,1) = det*Fb(1,1)
      Fb(2,1) = det*Fb(2,1)
      Fb(3,1) = det*Fb(3,1)
      Fb(4,1) = det*Fb(4,1)
      Fb(5,1) = det*Fb(5,1)
      Fb(6,1) = det*Fb(6,1)
      ! calculate second column of inverse
      Fb(1,2) = det*(-Fa(1,2)*W(15)+Fa(1,3)*W(14)-Fa(1,4)*W(13)+Fa(1,5)*W(12)-Fa(1,6)*W(11))
      Fb(2,2) = det*( Fa(1,1)*W(15)-Fa(1,3)*W(10)+Fa(1,4)*W(9)-Fa(1,5)*W(8)+Fa(1,6)*W(7))
      Fb(3,2) = det*(-Fa(1,1)*W(14)+Fa(1,2)*W(10)-Fa(1,4)*W(6)+Fa(1,5)*W(5)-Fa(1,6)*W(4))
      Fb(4,2) = det*( Fa(1,1)*W(13)-Fa(1,2)*W(9)+Fa(1,3)*W(6)-Fa(1,5)*W(3)+Fa(1,6)*W(2))
      Fb(5,2) = det*(-Fa(1,1)*W(12)+Fa(1,2)*W(8)-Fa(1,3)*W(5)+Fa(1,4)*W(3)-Fa(1,6)*W(1))
      Fb(6,2) = det*( Fa(1,1)*W(11)-Fa(1,2)*W(7)+Fa(1,3)*W(4)-Fa(1,4)*W(2)+Fa(1,5)*W(1))
      ! 3x3 determinants of rows 2-5-6
      V(1) = Fa(2,1)*U(6)-Fa(2,2)*U(2)+Fa(2,3)*U(1)
      V(2) = Fa(2,1)*U(7)-Fa(2,2)*U(3)+Fa(2,4)*U(1)
      V(3) = Fa(2,1)*U(8)-Fa(2,2)*U(4)+Fa(2,5)*U(1)
      V(4) = Fa(2,1)*U(9)-Fa(2,2)*U(5)+Fa(2,6)*U(1)
      V(5) = Fa(2,1)*U(10)-Fa(2,3)*U(3)+Fa(2,4)*U(2)
      V(6) = Fa(2,1)*U(11)-Fa(2,3)*U(4)+Fa(2,5)*U(2)
      V(7) = Fa(2,1)*U(12)-Fa(2,3)*U(5)+Fa(2,6)*U(2)
      V(8) = Fa(2,1)*U(13)-Fa(2,4)*U(4)+Fa(2,5)*U(3)
      V(9) = Fa(2,1)*U(14)-Fa(2,4)*U(5)+Fa(2,6)*U(3)
      V(10) = Fa(2,1)*U(15)-Fa(2,5)*U(5)+Fa(2,6)*U(4)
      V(11) = Fa(2,2)*U(10)-Fa(2,3)*U(7)+Fa(2,4)*U(6)
      V(12) = Fa(2,2)*U(11)-Fa(2,3)*U(8)+Fa(2,5)*U(6)
      V(13) = Fa(2,2)*U(12)-Fa(2,3)*U(9)+Fa(2,6)*U(6)
      V(14) = Fa(2,2)*U(13)-Fa(2,4)*U(8)+Fa(2,5)*U(7)
      V(15) = Fa(2,2)*U(14)-Fa(2,4)*U(9)+Fa(2,6)*U(7)
      V(16) = Fa(2,2)*U(15)-Fa(2,5)*U(9)+Fa(2,6)*U(8)
      V(17) = Fa(2,3)*U(13)-Fa(2,4)*U(11)+Fa(2,5)*U(10)
      V(18) = Fa(2,3)*U(14)-Fa(2,4)*U(12)+Fa(2,6)*U(10)
      V(19) = Fa(2,3)*U(15)-Fa(2,5)*U(12)+Fa(2,6)*U(11)
      V(20) = Fa(2,4)*U(15)-Fa(2,5)*U(14)+Fa(2,6)*U(13)
      ! 4x4 determinants of rows 1-2-5-6
      W(1) = Fa(1,1)*V(11)-Fa(1,2)*V(5)+Fa(1,3)*V(2)-Fa(1,4)*V(1)
      W(2) = Fa(1,1)*V(12)-Fa(1,2)*V(6)+Fa(1,3)*V(3)-Fa(1,5)*V(1)
      W(3) = Fa(1,1)*V(13)-Fa(1,2)*V(7)+Fa(1,3)*V(4)-Fa(1,6)*V(1)
      W(4) = Fa(1,1)*V(14)-Fa(1,2)*V(8)+Fa(1,4)*V(3)-Fa(1,5)*V(2)
      W(5) = Fa(1,1)*V(15)-Fa(1,2)*V(9)+Fa(1,4)*V(4)-Fa(1,6)*V(2)
      W(6) = Fa(1,1)*V(16)-Fa(1,2)*V(10)+Fa(1,5)*V(4)-Fa(1,6)*V(3)
      W(7) = Fa(1,1)*V(17)-Fa(1,3)*V(8)+Fa(1,4)*V(6)-Fa(1,5)*V(5)
      W(8) = Fa(1,1)*V(18)-Fa(1,3)*V(9)+Fa(1,4)*V(7)-Fa(1,6)*V(5)
      W(9) = Fa(1,1)*V(19)-Fa(1,3)*V(10)+Fa(1,5)*V(7)-Fa(1,6)*V(6)
      W(10) = Fa(1,1)*V(20)-Fa(1,4)*V(10)+Fa(1,5)*V(9)-Fa(1,6)*V(8)
      W(11) = Fa(1,2)*V(17)-Fa(1,3)*V(14)+Fa(1,4)*V(12)-Fa(1,5)*V(11)
      W(12) = Fa(1,2)*V(18)-Fa(1,3)*V(15)+Fa(1,4)*V(13)-Fa(1,6)*V(11)
      W(13) = Fa(1,2)*V(19)-Fa(1,3)*V(16)+Fa(1,5)*V(13)-Fa(1,6)*V(12)
      W(14) = Fa(1,2)*V(20)-Fa(1,4)*V(16)+Fa(1,5)*V(15)-Fa(1,6)*V(14)
      W(15) = Fa(1,3)*V(20)-Fa(1,4)*V(19)+Fa(1,5)*V(18)-Fa(1,6)*V(17)
      ! calculate third column of inverse
      Fb(1,3) = det*( Fa(4,2)*W(15)-Fa(4,3)*W(14)+Fa(4,4)*W(13)-Fa(4,5)*W(12)+Fa(4,6)*W(11))
      Fb(2,3) = det*(-Fa(4,1)*W(15)+Fa(4,3)*W(10)-Fa(4,4)*W(9)+Fa(4,5)*W(8)-Fa(4,6)*W(7))
      Fb(3,3) = det*( Fa(4,1)*W(14)-Fa(4,2)*W(10)+Fa(4,4)*W(6)-Fa(4,5)*W(5)+Fa(4,6)*W(4))
      Fb(4,3) = det*(-Fa(4,1)*W(13)+Fa(4,2)*W(9)-Fa(4,3)*W(6)+Fa(4,5)*W(3)-Fa(4,6)*W(2))
      Fb(5,3) = det*( Fa(4,1)*W(12)-Fa(4,2)*W(8)+Fa(4,3)*W(5)-Fa(4,4)*W(3)+Fa(4,6)*W(1))
      Fb(6,3) = det*(-Fa(4,1)*W(11)+Fa(4,2)*W(7)-Fa(4,3)*W(4)+Fa(4,4)*W(2)-Fa(4,5)*W(1))
      ! calculate fourth column of inverse
      Fb(1,4) = det*(-Fa(3,2)*W(15)+Fa(3,3)*W(14)-Fa(3,4)*W(13)+Fa(3,5)*W(12)-Fa(3,6)*W(11))
      Fb(2,4) = det*( Fa(3,1)*W(15)-Fa(3,3)*W(10)+Fa(3,4)*W(9)-Fa(3,5)*W(8)+Fa(3,6)*W(7))
      Fb(3,4) = det*(-Fa(3,1)*W(14)+Fa(3,2)*W(10)-Fa(3,4)*W(6)+Fa(3,5)*W(5)-Fa(3,6)*W(4))
      Fb(4,4) = det*( Fa(3,1)*W(13)-Fa(3,2)*W(9)+Fa(3,3)*W(6)-Fa(3,5)*W(3)+Fa(3,6)*W(2))
      Fb(5,4) = det*(-Fa(3,1)*W(12)+Fa(3,2)*W(8)-Fa(3,3)*W(5)+Fa(3,4)*W(3)-Fa(3,6)*W(1))
      Fb(6,4) = det*( Fa(3,1)*W(11)-Fa(3,2)*W(7)+Fa(3,3)*W(4)-Fa(3,4)*W(2)+Fa(3,5)*W(1))
      ! 2x2 determinants of rows 3-4
      U(1) = Fa(3,1)*Fa(4,2)-Fa(3,2)*Fa(4,1)
      U(2) = Fa(3,1)*Fa(4,3)-Fa(3,3)*Fa(4,1)
      U(3) = Fa(3,1)*Fa(4,4)-Fa(3,4)*Fa(4,1)
      U(4) = Fa(3,1)*Fa(4,5)-Fa(3,5)*Fa(4,1)
      U(5) = Fa(3,1)*Fa(4,6)-Fa(3,6)*Fa(4,1)
      U(6) = Fa(3,2)*Fa(4,3)-Fa(3,3)*Fa(4,2)
      U(7) = Fa(3,2)*Fa(4,4)-Fa(3,4)*Fa(4,2)
      U(8) = Fa(3,2)*Fa(4,5)-Fa(3,5)*Fa(4,2)
      U(9) = Fa(3,2)*Fa(4,6)-Fa(3,6)*Fa(4,2)
      U(10) = Fa(3,3)*Fa(4,4)-Fa(3,4)*Fa(4,3)
      U(11) = Fa(3,3)*Fa(4,5)-Fa(3,5)*Fa(4,3)
      U(12) = Fa(3,3)*Fa(4,6)-Fa(3,6)*Fa(4,3)
      U(13) = Fa(3,4)*Fa(4,5)-Fa(3,5)*Fa(4,4)
      U(14) = Fa(3,4)*Fa(4,6)-Fa(3,6)*Fa(4,4)
      U(15) = Fa(3,5)*Fa(4,6)-Fa(3,6)*Fa(4,5)
      ! 3x3 determinants of rows 2-3-4
      V(1) = Fa(2,1)*U(6)-Fa(2,2)*U(2)+Fa(2,3)*U(1)
      V(2) = Fa(2,1)*U(7)-Fa(2,2)*U(3)+Fa(2,4)*U(1)
      V(3) = Fa(2,1)*U(8)-Fa(2,2)*U(4)+Fa(2,5)*U(1)
      V(4) = Fa(2,1)*U(9)-Fa(2,2)*U(5)+Fa(2,6)*U(1)
      V(5) = Fa(2,1)*U(10)-Fa(2,3)*U(3)+Fa(2,4)*U(2)
      V(6) = Fa(2,1)*U(11)-Fa(2,3)*U(4)+Fa(2,5)*U(2)
      V(7) = Fa(2,1)*U(12)-Fa(2,3)*U(5)+Fa(2,6)*U(2)
      V(8) = Fa(2,1)*U(13)-Fa(2,4)*U(4)+Fa(2,5)*U(3)
      V(9) = Fa(2,1)*U(14)-Fa(2,4)*U(5)+Fa(2,6)*U(3)
      V(10) = Fa(2,1)*U(15)-Fa(2,5)*U(5)+Fa(2,6)*U(4)
      V(11) = Fa(2,2)*U(10)-Fa(2,3)*U(7)+Fa(2,4)*U(6)
      V(12) = Fa(2,2)*U(11)-Fa(2,3)*U(8)+Fa(2,5)*U(6)
      V(13) = Fa(2,2)*U(12)-Fa(2,3)*U(9)+Fa(2,6)*U(6)
      V(14) = Fa(2,2)*U(13)-Fa(2,4)*U(8)+Fa(2,5)*U(7)
      V(15) = Fa(2,2)*U(14)-Fa(2,4)*U(9)+Fa(2,6)*U(7)
      V(16) = Fa(2,2)*U(15)-Fa(2,5)*U(9)+Fa(2,6)*U(8)
      V(17) = Fa(2,3)*U(13)-Fa(2,4)*U(11)+Fa(2,5)*U(10)
      V(18) = Fa(2,3)*U(14)-Fa(2,4)*U(12)+Fa(2,6)*U(10)
      V(19) = Fa(2,3)*U(15)-Fa(2,5)*U(12)+Fa(2,6)*U(11)
      V(20) = Fa(2,4)*U(15)-Fa(2,5)*U(14)+Fa(2,6)*U(13)
      ! 4x4 determinants of rows 1-2-3-4
      W(1) = Fa(1,1)*V(11)-Fa(1,2)*V(5)+Fa(1,3)*V(2)-Fa(1,4)*V(1)
      W(2) = Fa(1,1)*V(12)-Fa(1,2)*V(6)+Fa(1,3)*V(3)-Fa(1,5)*V(1)
      W(3) = Fa(1,1)*V(13)-Fa(1,2)*V(7)+Fa(1,3)*V(4)-Fa(1,6)*V(1)
      W(4) = Fa(1,1)*V(14)-Fa(1,2)*V(8)+Fa(1,4)*V(3)-Fa(1,5)*V(2)
      W(5) = Fa(1,1)*V(15)-Fa(1,2)*V(9)+Fa(1,4)*V(4)-Fa(1,6)*V(2)
      W(6) = Fa(1,1)*V(16)-Fa(1,2)*V(10)+Fa(1,5)*V(4)-Fa(1,6)*V(3)
      W(7) = Fa(1,1)*V(17)-Fa(1,3)*V(8)+Fa(1,4)*V(6)-Fa(1,5)*V(5)
      W(8) = Fa(1,1)*V(18)-Fa(1,3)*V(9)+Fa(1,4)*V(7)-Fa(1,6)*V(5)
      W(9) = Fa(1,1)*V(19)-Fa(1,3)*V(10)+Fa(1,5)*V(7)-Fa(1,6)*V(6)
      W(10) = Fa(1,1)*V(20)-Fa(1,4)*V(10)+Fa(1,5)*V(9)-Fa(1,6)*V(8)
      W(11) = Fa(1,2)*V(17)-Fa(1,3)*V(14)+Fa(1,4)*V(12)-Fa(1,5)*V(11)
      W(12) = Fa(1,2)*V(18)-Fa(1,3)*V(15)+Fa(1,4)*V(13)-Fa(1,6)*V(11)
      W(13) = Fa(1,2)*V(19)-Fa(1,3)*V(16)+Fa(1,5)*V(13)-Fa(1,6)*V(12)
      W(14) = Fa(1,2)*V(20)-Fa(1,4)*V(16)+Fa(1,5)*V(15)-Fa(1,6)*V(14)
      W(15) = Fa(1,3)*V(20)-Fa(1,4)*V(19)+Fa(1,5)*V(18)-Fa(1,6)*V(17)
      ! calculate fifth column of inverse
      Fb(1,5) = det*( Fa(6,2)*W(15)-Fa(6,3)*W(14)+Fa(6,4)*W(13)-Fa(6,5)*W(12)+Fa(6,6)*W(11))
      Fb(2,5) = det*(-Fa(6,1)*W(15)+Fa(6,3)*W(10)-Fa(6,4)*W(9)+Fa(6,5)*W(8)-Fa(6,6)*W(7))
      Fb(3,5) = det*( Fa(6,1)*W(14)-Fa(6,2)*W(10)+Fa(6,4)*W(6)-Fa(6,5)*W(5)+Fa(6,6)*W(4))
      Fb(4,5) = det*(-Fa(6,1)*W(13)+Fa(6,2)*W(9)-Fa(6,3)*W(6)+Fa(6,5)*W(3)-Fa(6,6)*W(2))
      Fb(5,5) = det*( Fa(6,1)*W(12)-Fa(6,2)*W(8)+Fa(6,3)*W(5)-Fa(6,4)*W(3)+Fa(6,6)*W(1))
      Fb(6,5) = det*(-Fa(6,1)*W(11)+Fa(6,2)*W(7)-Fa(6,3)*W(4)+Fa(6,4)*W(2)-Fa(6,5)*W(1))
      ! calculate sixth column of inverse
      Fb(1,6) = det*(-Fa(5,2)*W(15)+Fa(5,3)*W(14)-Fa(5,4)*W(13)+Fa(5,5)*W(12)-Fa(5,6)*W(11))
      Fb(2,6) = det*( Fa(5,1)*W(15)-Fa(5,3)*W(10)+Fa(5,4)*W(9)-Fa(5,5)*W(8)+Fa(5,6)*W(7))
      Fb(3,6) = det*(-Fa(5,1)*W(14)+Fa(5,2)*W(10)-Fa(5,4)*W(6)+Fa(5,5)*W(5)-Fa(5,6)*W(4))
      Fb(4,6) = det*( Fa(5,1)*W(13)-Fa(5,2)*W(9)+Fa(5,3)*W(6)-Fa(5,5)*W(3)+Fa(5,6)*W(2))
      Fb(5,6) = det*(-Fa(5,1)*W(12)+Fa(5,2)*W(8)-Fa(5,3)*W(5)+Fa(5,4)*W(3)-Fa(5,6)*W(1))
      Fb(6,6) = det*( Fa(5,1)*W(11)-Fa(5,2)*W(7)+Fa(5,3)*W(4)-Fa(5,4)*W(2)+Fa(5,5)*W(1))

      bsuccess = .true.

    else

      bsuccess = .false.

    end if

  end subroutine mprim_invert6x6MatrixDirectSP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invertMatrixPivotDP(Da,ndim,bsuccess)

!<description>
  ! This subroutine directly inverts a (ndim x ndim) system with pivoting.
  ! 'Da' is a 2-dimensional (ndim x ndim) matrix and will be replaced
  ! by its inverse.
!</description>

!<input>
  ! Dimension of the matrix Da.
  integer, intent(in) :: ndim
!</input>

!<inputoutput>
  ! source square matrix to be inverted
  real(DP), dimension(ndim,ndim), intent(inout) :: Da

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Db is undefined.
  logical, intent(out) :: bsuccess
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(ndim) :: Kindx,Kindy

    real(DP) :: dpivot,daux
    integer :: idim1,idim2,ix,iy,indx,indy

    ! Perform factorisation of matrix Da

    bsuccess = .false.

    ! Initialisation
    Kindx=0
    Kindy=0

    do idim1=1,ndim

      ! Determine pivotal element
      dpivot=0

      do iy=1,ndim
        if (Kindy(iy) /= 0) cycle

        do ix=1,ndim
          if (Kindx(ix) /= 0) cycle

          if (abs(Da(ix,iy)) .le. abs(dpivot)) cycle
          dpivot=Da(ix,iy);  indx=ix;  indy=iy
        end do
      end do

      ! Return if pivotal element is zero
      if (dpivot .eq. 0.0_DP) return

      Kindx(indx)=indy;  Kindy(indy)=indx;  Da(indx,indy)=1._DP&
          &/dpivot

      do idim2=1,ndim
        if (idim2 == indy) cycle
        Da(1:indx-1,idim2)=Da(1:indx-1,idim2)-Da(1:indx-1,  &
            & indy)*Da(indx,idim2)/dpivot
        Da(indx+1:ndim,idim2)=Da(indx+1:ndim,idim2)-Da(indx+1:ndim&
            &,indy)*Da(indx,idim2)/dpivot
      end do

      do ix=1,ndim
        if (ix /= indx) Da(ix,indy)=Da(ix,indy)/dpivot
      end do

      do iy=1,ndim
        if (iy /= indy) Da(indx,iy)=-Da(indx,iy)/dpivot
      end do
    end do

    do ix=1,ndim
      if (Kindx(ix) == ix) cycle

      do iy=1,ndim
        if (Kindx(iy) == ix) exit
      end do

      do idim1=1,ndim
        daux=Da(ix,idim1)
        Da(ix,idim1)=Da(iy,idim1)
        Da(iy,idim1)=daux
      end do

      Kindx(iy)=Kindx(ix);  Kindx(ix)=ix
    end do

    do ix=1,ndim
      if (Kindy(ix) == ix) cycle

      do iy=1,ndim
        if (Kindy(iy) == ix) exit
      end do

      do idim1=1,ndim
        daux=Da(idim1,ix)
        Da(idim1,ix)=Da(idim1,iy)
        Da(idim1,iy)=daux
      end do

      Kindy(iy)=Kindy(ix);  Kindy(ix)=ix
    end do

    bsuccess = .true.

  end subroutine mprim_invertMatrixPivotDP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_invertMatrixPivotSP(Fa,ndim,bsuccess)

!<description>
  ! This subroutine directly inverts a (ndim x ndim) system with pivoting.
  ! 'Fa' is a 2-dimensional (ndim x ndim) matrix and will be replaced
  ! by its inverse.
!</description>

!<input>
  ! Dimension of the matrix Fa.
  integer, intent(in) :: ndim
!</input>

!<inputoutput>
  ! source square matrix to be inverted
  real(SP), dimension(ndim,ndim), intent(inout) :: Fa

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Db is undefined.
  logical, intent(out) :: bsuccess
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(ndim) :: Kindx,Kindy

    real(SP) :: fpivot,faux
    integer :: idim1,idim2,ix,iy,indx,indy

    ! Perform factorisation of matrix Fa

    bsuccess = .false.

    ! Initialisation
    Kindx=0
    Kindy=0

    do idim1=1,ndim

      ! Determine pivotal element
      fpivot=0

      do iy=1,ndim
        if (Kindy(iy) /= 0) cycle

        do ix=1,ndim
          if (Kindx(ix) /= 0) cycle

          if (abs(Fa(ix,iy)) .le. abs(fpivot)) cycle
          fpivot=Fa(ix,iy);  indx=ix;  indy=iy
        end do
      end do

      ! Return if pivotal element is zero
      if (fpivot .eq. 0.0_SP) return

      Kindx(indx)=indy;  Kindy(indy)=indx;  Fa(indx,indy)=1._SP&
          &/fpivot

      do idim2=1,ndim
        if (idim2 == indy) cycle
        Fa(1:indx-1,idim2)=Fa(1:indx-1,idim2)-Fa(1:indx-1,  &
            & indy)*Fa(indx,idim2)/fpivot
        Fa(indx+1:ndim,idim2)=Fa(indx+1:ndim,idim2)-Fa(indx+1:ndim&
            &,indy)*Fa(indx,idim2)/fpivot
      end do

      do ix=1,ndim
        if (ix /= indx) Fa(ix,indy)=Fa(ix,indy)/fpivot
      end do

      do iy=1,ndim
        if (iy /= indy) Fa(indx,iy)=-Fa(indx,iy)/fpivot
      end do
    end do

    do ix=1,ndim
      if (Kindx(ix) == ix) cycle

      do iy=1,ndim
        if (Kindx(iy) == ix) exit
      end do

      do idim1=1,ndim
        faux=Fa(ix,idim1)
        Fa(ix,idim1)=Fa(iy,idim1)
        Fa(iy,idim1)=faux
      end do

      Kindx(iy)=Kindx(ix);  Kindx(ix)=ix
    end do

    do ix=1,ndim
      if (Kindy(ix) == ix) cycle

      do iy=1,ndim
        if (Kindy(iy) == ix) exit
      end do

      do idim1=1,ndim
        faux=Fa(idim1,ix)
        Fa(idim1,ix)=Fa(idim1,iy)
        Fa(idim1,iy)=faux
      end do

      Kindy(iy)=Kindy(ix);  Kindy(ix)=ix
    end do

    bsuccess = .true.

  end subroutine mprim_invertMatrixPivotSP

  ! ***************************************************************************

!<function>

  elemental function kronecker(i,j) result(kron)

!<description>
    ! Compute the Kronecker delta symbol
    ! <tex> $$
    !  \delta_{ij}\left\{\begin{array}{ll}
    !  1 & i=j\\
    !  0 & i\ne j
    !  \end{array}\right. $$</tex>
!</description>

!<input>
    ! Evaluation points I and J
    integer, intent(in) :: i,j
!</input>

!<result>
    ! Kronecker delta symbol
    integer :: kron
!</result>
!</function>

    kron=merge(1,0,i==j)
  end function kronecker

  !************************************************************************

!<subroutine>

  elemental subroutine mprim_linearRescaleDP(dx,da,db,dc,dd,dy)

!<description>
  ! Scales a coordinate x linearly from the interval [a,b] to the
  ! interval [c,d].
!</description>

!<input>
  ! coordinate to be rescaled
  real(DP), intent(in) :: dx

  ! [a,b] - source interval
  real(DP), intent(in) :: da,db

  ! [c,d] - destination interval
  real(DP), intent(in) :: dc,dd
!</input>

!<output>
  ! Rescaled coordinate
  real(DP), intent(out) :: dy
!</output>

!</subroutine>

    real(DP) :: d1,d2,d3

    ! Calculate the coefficients of the transformation
    !    D1*A+D2 = C, D1*B+D2 = D.
    ! Use them to calculate Y=D1*X+D2.

    if (dA .eq. db) then
      dy = dc
      return
    end if

    d3 = 1.0_DP/(da-db)
    d1 = (dc-dd)*d3
    d2 = (-db*dc+da*dd)*d3

    dy = d1*dx+d2

  end subroutine mprim_linearRescaleDP

  !************************************************************************

!<subroutine>

  elemental subroutine mprim_linearRescaleSP(fx,fa,fb,fc,fd,fy)

!<description>
  ! Scales a coordinate x linearly from the interval [a,b] to the
  ! interval [c,d].
!</description>

!<input>
  ! coordinate to be rescaled
  real(SP), intent(in) :: fx

  ! [a,b] - source interval
  real(SP), intent(in) :: fa,fb

  ! [c,d] - destination interval
  real(SP), intent(in) :: fc,fd
!</input>

!<output>
  ! Rescaled coordinate
  real(SP), intent(out) :: fy
!</output>

!</subroutine>

    real(SP) :: f1,f2,f3

    ! Calculate the coefficients of the transformation
    !    F1*A+F2 = C, F1*B+F2 = D.
    ! Use them to calculate Y=F1*X+F2.

    if (fA .eq. fb) then
      fy = fc
      return
    end if

    f3 = 1.0_SP/(fa-fb)
    f1 = (fc-fd)*f3
    f2 = (-fb*fc+fa*fd)*f3

    fy = f1*fx+f2

  end subroutine mprim_linearRescaleSP

  ! ***************************************************************************

!<subroutine>

  elemental subroutine mprim_quadraticInterpolationDP (dx,d1,d2,d3,dy)

!<description>
  ! <tex>
  ! Calculates a quadratic interpolation. dx is a value in the range $[-1,1]$.
  ! The routine calculates the value $dy:=p(dx)$ with $p(.)$ being the quadratic
  ! interpolation polynomial with $p(-1)=d1$, $p(0)=d2$ and $p(1)=d3$.
  ! </tex>
!</description>

!<input>
  ! The parameter value in the range <tex>$ [-1,1] $</tex> where the polynomial should be evaluated.
  real(DP), intent(in) :: dx

  ! The value <tex>$ p(-1) $</tex>.
  real(DP), intent(in) :: d1

  ! The value <tex>$ p(0) $</tex>.
  real(DP), intent(in) :: d2

  ! The value $ p(1) $</tex>.
  real(DP), intent(in) :: d3
!</input>

!<output>
  ! The value <tex>$ p(dx) $</tex>.
  real(DP), intent(out) :: dy
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

  end subroutine mprim_quadraticInterpolationDP

  ! ***************************************************************************

!<subroutine>

  elemental subroutine mprim_quadraticInterpolationSP (fx,f1,f2,f3,fy)

!<description>
  ! <tex>
  ! Calculates a quadratic interpolation. fx is a value in the range $[-1,1]$.
  ! The routine calculates the value $fy:=p(fx)$ with $p(.)$ being the quadratic
  ! interpolation polynomial with $p(-1)=f1$, $p(0)=f2$ and $p(1)=f3$.
  ! </tex>
!</description>

!<input>
  ! The parameter value in the range <tex>$ [-1,1] $</tex> where the polynomial should be evaluated.
  real(SP), intent(in) :: fx

  ! The value <tex>$ p(-1) $</tex>.
  real(SP), intent(in) :: f1

  ! The value <tex>$ p(0) $</tex>.
  real(SP), intent(in) :: f2

  ! The value $ p(1) $</tex>.
  real(SP), intent(in) :: f3
!</input>

!<output>
  ! The value <tex>$ p(fx) $</tex>.
  real(SP), intent(out) :: fy
!</output>

!</subroutine>

    ! The polynomial p(t) = a + bt + ct^2 has to fulfill:
    !   p(-1)=f1, p(0)=f2, p(1)=f3.
    !
    ! So the polynomial has the shape:
    !   p(t) = f2  +  1/2(f3-f1)t  +  1/2(f3-2f2+f1)t^2
    !
    ! The Horner scheme gives us:
    !   p(t) = 1/2 ( (f3-2f2+f1)t + (f3-f1) ) t + f2

    fy = 0.5_SP * ( (f3 - 2.0_SP*f2 + f1)*fx + (f3-f1) ) * fx + f2

  end subroutine mprim_quadraticInterpolationSP

  !************************************************************************

!<subroutine>

  subroutine mprim_SVD_factoriseDP(Da,Dd,Db,ndim1,ndim2,ndim3,btransposedOpt)

!<description>
    ! <tex>
    ! This subroutine computes the factorisation for a singular value
    ! decomposition. Given an ndim2-by-ndim1 matrix Da, the routine
    ! decomposes it into the product $ A = U * D * B^T $
    ! where $U$ overwrites the matrix A and is returned in its memory
    ! position, $D$ is a diagonal matrix returned as vector Dd and
    ! $B$ is an n-by-n square matrix which is returned (instead
    ! of its transpose) as matrix Db.
    !
    ! The optional parameter btransposedOpt can be used to indicate
    ! that matrix A is stored in transposed format. Note that
    ! matrices D and B are not affected by this fact.
    !
    ! The details of this algorithm are given in numerical recipes in F90.
    ! </tex>
!</description>

!<input>
    ! Dimensions of the rectangular matrix Da
    integer, intent(in) :: ndim1,ndim2

    ! Dimension of the diagonal matrix Dd and the suare matrix Db
    integer, intent(in) :: ndim3

    ! OPTIONAL: Flag to indicate if the matrix A is transposed
    logical, intent(in), optional :: btransposedOpt
!</input>

!<inputoutput>
    ! Matrix that should be factorised on input.
    ! Matrix U on output.
    real(DP), dimension(ndim1,ndim2), intent(inout) :: Da
!</inputoutput>

!<output>
    ! Diagonal matrix
    real(DP), dimension(ndim3), intent(out) :: Dd

    ! Square matrix
    real(DP), dimension(ndim3,ndim3), intent(out) :: Db
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(:), allocatable :: rv1
    real(DP) :: g, scale, anorm, s, f, h, c, x, y, z
    integer  :: its, i, j, jj, k, l, nm, n, m, idot
    integer, parameter :: MAX_ITS = 30
    logical :: btransposed

    ! Do we have an optional transposed flag?
    if (present(btransposedOpt)) then
      btransposed = btransposedOpt
    else
      btransposed = .false.
    end if

    ! Is the matrix transposed?
    if (btransposed) then
      n=ndim1; m=ndim2
    else
      n=ndim2; m=ndim1
    end if

    ! Check if number of equations is larger than number of unknowns
    if (m < n) then
      call output_line('Fewer equations than unknowns!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mprim_SVD_factoriseDP')
      call sys_halt()
    end if

    ! Allocate internal memory
    allocate(rv1(ndim3))

    if (btransposed) then

      ! Householder reduction to bidiagonal form
      g     = 0.0_DP
      scale = 0.0_DP
      anorm = 0.0_DP

      do i = 1, n

        l      = i+1
        rv1(i) = scale*g
        g      = 0.0_DP
        scale  = 0.0_DP

        if (i .le. m) then
          scale = sum(abs(Da(i,i:m)))
          if (scale .ne. 0.0_DP) then
            Da(i,i:m) = Da(i,i:m)/scale
            s = 0.0_DP
            do idot = i, m
              s = s+Da(i,idot)*Da(i,idot)
            end do
            f = Da(i,i)
            g = -sign(sqrt(s),f)
            h = f*g-s
            Da(i,i) = f-g
            if (i .ne. n) then
              do j = l, n
                s = 0.0_DP
                do k = i, m
                  s = s+Da(i,k)*Da(j,k)
                end do
                f = s/h
                do k = i, m
                  Da(j,k) = Da(j,k)+f*Da(i,k)
                end do
              end do
            end if
            do k = i, m
              Da(i,k) = scale*Da(i,k)
            end do
          end if
        end if

        Dd(i) = scale*g
        g     = 0.0_DP
        s     = 0.0_DP
        scale = 0.0_DP

        if ((i .le. m) .and. (i .ne. n)) then
          scale = sum(abs(Da(l:n,i)))
          if (scale .ne. 0.0_DP) then
            do k = l, n
              Da(k,i) = Da(k,i)/scale
              s = s+Da(k,i)*Da(k,i)
            end do
            f = Da(l,i)
            g = -sign(sqrt(s),f)
            h = f*g-s
            Da(l,i) = f-g
            do k = l, n
              rv1(k) = Da(k,i)/h
            end do

            if (i .ne. m) then
              do j = l, m
                s = 0.0_DP
                do k = l, n
                  s = s+Da(k,j)*Da(k,i)
                end do
                do k = l, n
                  Da(k,j) = Da(k,j)+s*rv1(k)
                end do
              end do
            end if

            do k = l, n
              Da(k,i) = scale*Da(k,i)
            end do
          end if
        end if
        anorm = max(anorm,(abs(Dd(i))+abs(rv1(i))))
      end do

      ! Accumulation of right-hand transformations
      do i = n, 1, -1
        if (i .lt. n) then
          if (g .ne. 0.0_DP) then
            do j = l, n
              Db(j,i) = (Da(j,i)/Da(l,i))/g
            end do
            do j = l, n
              s = 0.0_DP
              do k = l, n
                s = s+Da(k,i)*Db(k,j)
              end do
              do k = l, n
                Db(k,j) = Db(k,j)+s*Db(k,i)
              end do
            end do
          end if
          do j = l, n
            Db(i,j) = 0.0_DP
            Db(j,i) = 0.0_DP
          end do
        end if
        Db(i,i) = 1.0_DP
        g = rv1(i)
        l = i
      end do

      ! Accumulation of left-hand transformations
      do i = n, 1, -1
        l = i+1
        g = Dd(i)
        if (i .lt. n) then
          do j = l, n
            Da(j,i) = 0.0_DP
          end do
        end if
        if (g .ne. 0.0_DP) then
          g = 1.0_DP/g
          if (i .ne. n) then
            do j = l, n
              s = 0.0_DP
              do k = l, m
                s = s+Da(i,k)*Da(j,k)
              end do
              f = (s/Da(i,i))*g
              do k = i, m
                Da(j,k) = Da(j,k)+f*Da(i,k)
              end do
            end do
          end if
          do j = i, m
            Da(i,j) = Da(i,j)*g
          end do
        else
          do j = i, m
            Da(i,j) = 0.0_DP
          end do
        end if
        Da(i,i) = Da(i,i)+1.0_DP
      end do

      ! Diagonalisation of the bidiagonal form
      do k = n, 1, -1
        do its = 1, MAX_ITS
          do l = k, 1, -1
            nm = l-1
            if ((abs(rv1(l))+anorm) .eq. anorm) goto 5
            if ((abs(Dd(nm))+anorm) .eq. anorm) goto 4
          end do
4         continue
          c = 0.0_DP
          s = 1.0_DP
          do i = l, k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            if ((abs(f)+anorm) .eq. anorm) goto 5
            g = Dd(i)
            h = sqrt(f*f+g*g)
            Dd(i) = h
            h = 1.0_DP/h
            c = g*h
            s = -(f*h)
            do j = 1, m
              y = Da(nm,j)
              z = Da(i,j)
              Da(nm,j) = (y*c)+(z*s)
              Da(i,j) = -(y*s)+(z*c)
            end do
          end do
5         continue
          z = Dd(k)
          if (l .eq. k) then
            if (z .lt. 0.0_DP) then
              Dd(k) = -z
              do j = 1, n
                Db(j,k) = -Db(j,k)
              end do
            end if
            goto 6
          end if
          if (its .eq. MAX_ITS) then
            call output_line('Convergence failed!',&
                OU_CLASS_ERROR,OU_MODE_STD,'mprim_SVD_factoriseDP')
            call sys_halt()
          end if
          x = Dd(l)
          nm = k-1
          y = Dd(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_DP*h*y)
          g = sqrt(f*f+1.0_DP)
          f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x

          ! Next QR transformation
          c = 1.0_DP
          s = 1.0_DP
          do j = l, nm
            i = j+1
            g = rv1(i)
            y = Dd(i)
            h = s*g
            g = c*g
            z = sqrt(f*f+h*h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c)+(g*s)
            g = -(x*s)+(g*c)
            h = y*s
            y = y*c
            do jj = 1, n
              x = Db(jj,j)
              z = Db(jj,i)
              Db(jj,j) = (x*c)+(z*s)
              Db(jj,i) = -(x*s)+(z*c)
            end do
            z = sqrt(f*f+h*h)
            Dd(j) = z
            if (z .ne. 0.0_DP) then
              z = 1.0_DP/z
              c = f*z
              s = h*z
            end if
            f = (c*g)+(s*y)
            x = -(s*g)+(c*y)
            do jj = 1, m
              y = Da(j,jj)
              z = Da(i,jj)
              Da(j,jj) = (y*c)+(z*s)
              Da(i,jj) = -(y*s)+(z*c)
            end do
          end do
          rv1(l) = 0.0_DP
          rv1(k) = f
          Dd(k) = x
        end do
6       continue
      end do

    else

      ! Householder reduction to bidiagonal form
      g     = 0.0_DP
      scale = 0.0_DP
      anorm = 0.0_DP

      do i = 1, n

        l      = i+1
        rv1(i) = scale*g
        g      = 0.0_DP
        scale  = 0.0_DP

        if (i .le. m) then
          scale = sum(abs(Da(i:m,i)))
          if (scale .ne. 0.0_DP) then
            Da(i:m,i) = Da(i:m,i)/scale
            s = 0.0_DP
            do idot = i, m
              s = s+Da(idot,i)*Da(idot,i)
            end do
            f = Da(i,i)
            g = -sign(sqrt(s),f)
            h = f*g-s
            Da(i,i) = f-g
            if (i .ne. n) then
              do j = l, n
                s = 0.0_DP
                do k = i, m
                  s = s+Da(k,i)*Da(k,j)
                end do
                f = s/h
                do k = i, m
                  Da(k,j) = Da(k,j)+f*Da(k,i)
                end do
              end do
            end if
            do k = i, m
              Da(k,i) = scale*Da(k,i)
            end do
          end if
        end if

        Dd(i) = scale*g
        g     = 0.0_DP
        s     = 0.0_DP
        scale = 0.0_DP

        if ((i .le. m) .and. (i .ne. n)) then
          scale = sum(abs(Da(i,l:n)))
          if (scale .ne. 0.0_DP) then
            do k = l, n
              Da(i,k) = Da(i,k)/scale
              s = s+Da(i,k)*Da(i,k)
            end do
            f = Da(i,l)
            g = -sign(sqrt(s),f)
            h = f*g-s
            Da(i,l) = f-g
            do k = l, n
              rv1(k) = Da(i,k)/h
            end do

            if (i .ne. m) then
              do j = l, m
                s = 0.0_DP
                do k = l, n
                  s = s+Da(j,k)*Da(i,k)
                end do
                do k = l, n
                  Da(j,k) = Da(j,k)+s*rv1(k)
                end do
              end do
            end if

            do k = l, n
              Da(i,k) = scale*Da(i,k)
            end do
          end if
        end if
        anorm = max(anorm,(abs(Dd(i))+abs(rv1(i))))
      end do

      ! Accumulation of right-hand transformations
      do i = n, 1, -1
        if (i .lt. n) then
          if (g .ne. 0.0_DP) then
            do j = l, n
              Db(j,i) = (Da(i,j)/Da(i,l))/g
            end do
            do j = l, n
              s = 0.0_DP
              do k = l, n
                s = s+Da(i,k)*Db(k,j)
              end do
              do k = l, n
                Db(k,j) = Db(k,j)+s*Db(k,i)
              end do
            end do
          end if
          do j = l, n
            Db(i,j) = 0.0_DP
            Db(j,i) = 0.0_DP
          end do
        end if
        Db(i,i) = 1.0_DP
        g = rv1(i)
        l = i
      end do

      ! Accumulation of left-hand transformations
      do i = n, 1, -1
        l = i+1
        g = Dd(i)
        if (i .lt. n) then
          do j = l, n
            Da(i,j) = 0.0_DP
          end do
        end if
        if (g .ne. 0.0_DP) then
          g = 1.0_DP/g
          if (i .ne. n) then
            do j = l, n
              s = 0.0_DP
              do k = l, m
                s = s+Da(k,i)*Da(k,j)
              end do
              f = (s/Da(i,i))*g
              do k = i, m
                Da(k,j) = Da(k,j)+f*Da(k,i)
              end do
            end do
          end if
          do j = i, m
            Da(j,i) = Da(j,i)*g
          end do
        else
          do j = i, m
            Da(j,i) = 0.0_DP
          end do
        end if
        Da(i,i) = Da(i,i)+1.0_DP
      end do

      ! Diagonalisation of the bidiagonal form
      do k = n, 1, -1
        do its = 1, MAX_ITS
          do l = k, 1, -1
            nm = l-1
            if ((abs(rv1(l))+anorm) .eq. anorm) goto 2
            if ((abs(Dd(nm))+anorm) .eq. anorm) goto 1
          end do
1         continue
          c = 0.0_DP
          s = 1.0_DP
          do i = l, k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            if ((abs(f)+anorm) .eq. anorm) goto 2
            g = Dd(i)
            h = sqrt(f*f+g*g)
            Dd(i) = h
            h = 1.0_DP/h
            c = g*h
            s = -(f*h)
            do j = 1, m
              y = Da(j,nm)
              z = Da(j,i)
              Da(j,nm) = (y*c)+(z*s)
              Da(j,i) = -(y*s)+(z*c)
            end do
          end do
2         continue
          z = Dd(k)
          if (l .eq. k) then
            if (z .lt. 0.0_DP) then
              Dd(k) = -z
              do j = 1, n
                Db(j,k) = -Db(j,k)
              end do
            end if
            goto 3
          end if
          if (its .eq. MAX_ITS) then
            call output_line('Convergence failed!',&
                OU_CLASS_ERROR,OU_MODE_STD,'mprim_SVD_factoriseDP')
            call sys_halt()
          end if
          x = Dd(l)
          nm = k-1
          y = Dd(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_DP*h*y)
          g = sqrt(f*f+1.0_DP)
          f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x

          ! Next QR transformation
          c = 1.0_DP
          s = 1.0_DP
          do j = l, nm
            i = j+1
            g = rv1(i)
            y = Dd(i)
            h = s*g
            g = c*g
            z = sqrt(f*f+h*h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c)+(g*s)
            g = -(x*s)+(g*c)
            h = y*s
            y = y*c
            do jj = 1, n
              x = Db(jj,j)
              z = Db(jj,i)
              Db(jj,j) = (x*c)+(z*s)
              Db(jj,i) = -(x*s)+(z*c)
            end do
            z = sqrt(f*f+h*h)
            Dd(j) = z
            if (z .ne. 0.0_DP) then
              z = 1.0_DP/z
              c = f*z
              s = h*z
            end if
            f = (c*g)+(s*y)
            x = -(s*g)+(c*y)
            do jj = 1, m
              y = Da(jj,j)
              z = Da(jj,i)
              Da(jj,j) = (y*c)+(z*s)
              Da(jj,i) = -(y*s)+(z*c)
            end do
          end do
          rv1(l) = 0.0_DP
          rv1(k) = f
          Dd(k) = x
        end do
3       continue
      end do

    end if

    ! Deallocate internal memory
    deallocate(rv1)

  end subroutine mprim_SVD_factoriseDP

  !************************************************************************

!<subroutine>

  subroutine mprim_SVD_factoriseSP(Fa,Fd,Fb,ndim1,ndim2,ndim3,btransposedOpt)

!<description>
    ! <tex>
    ! This subroutine computes the factorisation for a singular value
    ! decomposition. Given an ndim2-by-ndim1 matrix Fa, the routine
    ! decomposes it into the product $ A = U * D * B^T $
    ! where $U$ overwrites the matrix A and is returned in its memory
    ! position, $D$ is a diagonal matrix returned as vector Fd and
    ! $B$ is an n-by-n square matrix which is returned (instead
    ! of its transpose) as matrix Fb.
    !
    ! The optional parameter btransposedOpt can be used to indicate
    ! that matrix A is stored in transposed format. Note that
    ! matrices D and B are not affected by this fact.
    !
    ! The details of this algorithm are given in numerical recipes in F90.
    ! </tex>
!</description>

!<input>
    ! Dimensions of the rectangular matrix
    integer, intent(in) :: ndim1,ndim2,ndim3

    ! OPTIONAL: Flag to indicate if the matrix A is transposed
    logical, intent(in), optional :: btransposedOpt
!</input>

!<inputoutput>
    ! Matrix that should be factorised on input.
    ! Matrix U on output.
    real(SP), dimension(ndim1,ndim2), intent(inout) :: Fa
!</inputoutput>

!<output>
    ! Diagonal matrix
    real(SP), dimension(ndim3), intent(out) :: Fd

    ! Square matrix
    real(SP), dimension(ndim3,ndim3), intent(out) :: Fb
!</output>
!</subroutine>

    ! local variables
    real(SP), dimension(:), allocatable :: rv1
    real(SP) :: g, scale, anorm, s, f, h, c, x, y, z
    integer  :: its, i, j, jj, k, l, nm, n, m, idot
    integer, parameter :: MAX_ITS = 30
    logical :: btransposed

    ! Do we have an optional transposed flag?
    if (present(btransposedOpt)) then
      btransposed = btransposedOpt
    else
      btransposed = .false.
    end if

    ! Is the matrix transposed?
    if (btransposed) then
      n=ndim1; m=ndim2
    else
      n=ndim2; m=ndim1
    end if

    ! Check if number of equations is larger than number of unknowns
    if (m < n) then
      call output_line('Fewer equations than unknowns!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mprim_SVD_factoriseSP')
      call sys_halt()
    end if

    ! Allocate internal memory
    allocate(rv1(ndim3))

    if (btransposed) then

      ! Householder reduction to bidiagonal form
      g     = 0.0_SP
      scale = 0.0_SP
      anorm = 0.0_SP

      do i = 1, n

        l      = i+1
        rv1(i) = scale*g
        g      = 0.0_SP
        scale  = 0.0_SP

        if (i .le. m) then
          scale = sum(abs(Fa(i,i:m)))
          if (scale .ne. 0.0_SP) then
            Fa(i,i:m) = Fa(i,i:m)/scale
            s = 0.0_SP
            do idot = i, m
              s = s+Fa(i,idot)*Fa(i,idot)
            end do
            f = Fa(i,i)
            g = -sign(sqrt(s),f)
            h = f*g-s
            Fa(i,i) = f-g
            if (i .ne. n) then
              do j = l, n
                s = 0.0_SP
                do k = i, m
                  s = s+Fa(i,k)*Fa(j,k)
                end do
                f = s/h
                do k = i, m
                  Fa(j,k) = Fa(j,k)+f*Fa(i,k)
                end do
              end do
            end if
            do k = i, m
              Fa(i,k) = scale*Fa(i,k)
            end do
          end if
        end if

        Fd(i) = scale*g
        g     = 0.0_SP
        s     = 0.0_SP
        scale = 0.0_SP

        if ((i .le. m) .and. (i .ne. n)) then
          scale = sum(abs(Fa(l:n,i)))
          if (scale .ne. 0.0_SP) then
            do k = l, n
              Fa(k,i) = Fa(k,i)/scale
              s = s+Fa(k,i)*Fa(k,i)
            end do
            f = Fa(l,i)
            g = -sign(sqrt(s),f)
            h = f*g-s
            Fa(l,i) = f-g
            do k = l, n
              rv1(k) = Fa(k,i)/h
            end do

            if (i .ne. m) then
              do j = l, m
                s = 0.0_SP
                do k = l, n
                  s = s+Fa(k,j)*Fa(k,i)
                end do
                do k = l, n
                  Fa(k,j) = Fa(k,j)+s*rv1(k)
                end do
              end do
            end if

            do k = l, n
              Fa(k,i) = scale*Fa(k,i)
            end do
          end if
        end if
        anorm = max(anorm,(abs(Fd(i))+abs(rv1(i))))
      end do

      ! Accumulation of right-hand transformations
      do i = n, 1, -1
        if (i .lt. n) then
          if (g .ne. 0.0_SP) then
            do j = l, n
              Fb(j,i) = (Fa(j,i)/Fa(l,i))/g
            end do
            do j = l, n
              s = 0.0_SP
              do k = l, n
                s = s+Fa(k,i)*Fb(k,j)
              end do
              do k = l, n
                Fb(k,j) = Fb(k,j)+s*Fb(k,i)
              end do
            end do
          end if
          do j = l, n
            Fb(i,j) = 0.0_SP
            Fb(j,i) = 0.0_SP
          end do
        end if
        Fb(i,i) = 1.0_SP
        g = rv1(i)
        l = i
      end do

      ! Accumulation of left-hand transformations
      do i = n, 1, -1
        l = i+1
        g = Fd(i)
        if (i .lt. n) then
          do j = l, n
            Fa(j,i) = 0.0_SP
          end do
        end if
        if (g .ne. 0.0_SP) then
          g = 1.0_SP/g
          if (i .ne. n) then
            do j = l, n
              s = 0.0_SP
              do k = l, m
                s = s+Fa(i,k)*Fa(j,k)
              end do
              f = (s/Fa(i,i))*g
              do k = i, m
                Fa(j,k) = Fa(j,k)+f*Fa(i,k)
              end do
            end do
          end if
          do j = i, m
            Fa(i,j) = Fa(i,j)*g
          end do
        else
          do j = i, m
            Fa(i,j) = 0.0_SP
          end do
        end if
        Fa(i,i) = Fa(i,i)+1.0_SP
      end do

      ! Diagonalisation of the bidiagonal form
      do k = n, 1, -1
        do its = 1, MAX_ITS
          do l = k, 1, -1
            nm = l-1
            if ((abs(rv1(l))+anorm) .eq. anorm) goto 5
            if ((abs(Fd(nm))+anorm) .eq. anorm) goto 4
          end do
4         continue
          c = 0.0_SP
          s = 1.0_SP
          do i = l, k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            if ((abs(f)+anorm) .eq. anorm) goto 5
            g = Fd(i)
            h = sqrt(f*f+g*g)
            Fd(i) = h
            h = 1.0_SP/h
            c = g*h
            s = -(f*h)
            do j = 1, m
              y = Fa(nm,j)
              z = Fa(i,j)
              Fa(nm,j) = (y*c)+(z*s)
              Fa(i,j) = -(y*s)+(z*c)
            end do
          end do
5         continue
          z = Fd(k)
          if (l .eq. k) then
            if (z .lt. 0.0_SP) then
              Fd(k) = -z
              do j = 1, n
                Fb(j,k) = -Fb(j,k)
              end do
            end if
            goto 6
          end if
          if (its .eq. MAX_ITS) then
            call output_line('Convergence failed!',&
                OU_CLASS_ERROR,OU_MODE_STD,'mprim_SVD_factoriseSP')
            call sys_halt()
          end if
          x = Fd(l)
          nm = k-1
          y = Fd(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_SP*h*y)
          g = sqrt(f*f+1.0_SP)
          f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x

          ! Next QR transformation
          c = 1.0_SP
          s = 1.0_SP
          do j = l, nm
            i = j+1
            g = rv1(i)
            y = Fd(i)
            h = s*g
            g = c*g
            z = sqrt(f*f+h*h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c)+(g*s)
            g = -(x*s)+(g*c)
            h = y*s
            y = y*c
            do jj = 1, n
              x = Fb(jj,j)
              z = Fb(jj,i)
              Fb(jj,j) = (x*c)+(z*s)
              Fb(jj,i) = -(x*s)+(z*c)
            end do
            z = sqrt(f*f+h*h)
            Fd(j) = z
            if (z .ne. 0.0_SP) then
              z = 1.0_SP/z
              c = f*z
              s = h*z
            end if
            f = (c*g)+(s*y)
            x = -(s*g)+(c*y)
            do jj = 1, m
              y = Fa(j,jj)
              z = Fa(i,jj)
              Fa(j,jj) = (y*c)+(z*s)
              Fa(i,jj) = -(y*s)+(z*c)
            end do
          end do
          rv1(l) = 0.0_SP
          rv1(k) = f
          Fd(k) = x
        end do
6       continue
      end do

    else

      ! Householder reduction to bidiagonal form
      g     = 0.0_SP
      scale = 0.0_SP
      anorm = 0.0_SP

      do i = 1, n

        l      = i+1
        rv1(i) = scale*g
        g      = 0.0_SP
        scale  = 0.0_SP

        if (i .le. m) then
          scale = sum(abs(Fa(i:m,i)))
          if (scale .ne. 0.0_SP) then
            Fa(i:m,i) = Fa(i:m,i)/scale
            s = 0.0_SP
            do idot = i, m
              s = s+Fa(idot,i)*Fa(idot,i)
            end do
            f = Fa(i,i)
            g = -sign(sqrt(s),f)
            h = f*g-s
            Fa(i,i) = f-g
            if (i .ne. n) then
              do j = l, n
                s = 0.0_SP
                do k = i, m
                  s = s+Fa(k,i)*Fa(k,j)
                end do
                f = s/h
                do k = i, m
                  Fa(k,j) = Fa(k,j)+f*Fa(k,i)
                end do
              end do
            end if
            do k = i, m
              Fa(k,i) = scale*Fa(k,i)
            end do
          end if
        end if

        Fd(i) = scale*g
        g     = 0.0_SP
        s     = 0.0_SP
        scale = 0.0_SP

        if ((i .le. m) .and. (i .ne. n)) then
          scale = sum(abs(Fa(i,l:n)))
          if (scale .ne. 0.0_SP) then
            do k = l, n
              Fa(i,k) = Fa(i,k)/scale
              s = s+Fa(i,k)*Fa(i,k)
            end do
            f = Fa(i,l)
            g = -sign(sqrt(s),f)
            h = f*g-s
            Fa(i,l) = f-g
            do k = l, n
              rv1(k) = Fa(i,k)/h
            end do

            if (i .ne. m) then
              do j = l, m
                s = 0.0_SP
                do k = l, n
                  s = s+Fa(j,k)*Fa(i,k)
                end do
                do k = l, n
                  Fa(j,k) = Fa(j,k)+s*rv1(k)
                end do
              end do
            end if

            do k = l, n
              Fa(i,k) = scale*Fa(i,k)
            end do
          end if
        end if
        anorm = max(anorm,(abs(Fd(i))+abs(rv1(i))))
      end do

      ! Accumulation of right-hand transformations
      do i = n, 1, -1
        if (i .lt. n) then
          if (g .ne. 0.0_SP) then
            do j = l, n
              Fb(j,i) = (Fa(i,j)/Fa(i,l))/g
            end do
            do j = l, n
              s = 0.0_SP
              do k = l, n
                s = s+Fa(i,k)*Fb(k,j)
              end do
              do k = l, n
                Fb(k,j) = Fb(k,j)+s*Fb(k,i)
              end do
            end do
          end if
          do j = l, n
            Fb(i,j) = 0.0_SP
            Fb(j,i) = 0.0_SP
          end do
        end if
        Fb(i,i) = 1.0_SP
        g = rv1(i)
        l = i
      end do

      ! Accumulation of left-hand transformations
      do i = n, 1, -1
        l = i+1
        g = Fd(i)
        if (i .lt. n) then
          do j = l, n
            Fa(i,j) = 0.0_SP
          end do
        end if
        if (g .ne. 0.0_SP) then
          g = 1.0_SP/g
          if (i .ne. n) then
            do j = l, n
              s = 0.0_SP
              do k = l, m
                s = s+Fa(k,i)*Fa(k,j)
              end do
              f = (s/Fa(i,i))*g
              do k = i, m
                Fa(k,j) = Fa(k,j)+f*Fa(k,i)
              end do
            end do
          end if
          do j = i, m
            Fa(j,i) = Fa(j,i)*g
          end do
        else
          do j = i, m
            Fa(j,i) = 0.0_SP
          end do
        end if
        Fa(i,i) = Fa(i,i)+1.0_SP
      end do

      ! Diagonalisation of the bidiagonal form
      do k = n, 1, -1
        do its = 1, MAX_ITS
          do l = k, 1, -1
            nm = l-1
            if ((abs(rv1(l))+anorm) .eq. anorm) goto 2
            if ((abs(Fd(nm))+anorm) .eq. anorm) goto 1
          end do
1         continue
          c = 0.0_SP
          s = 1.0_SP
          do i = l, k
            f = s*rv1(i)
            rv1(i) = c*rv1(i)
            if ((abs(f)+anorm) .eq. anorm) goto 2
            g = Fd(i)
            h = sqrt(f*f+g*g)
            Fd(i) = h
            h = 1.0_SP/h
            c = g*h
            s = -(f*h)
            do j = 1, m
              y = Fa(j,nm)
              z = Fa(j,i)
              Fa(j,nm) = (y*c)+(z*s)
              Fa(j,i) = -(y*s)+(z*c)
            end do
          end do
2         continue
          z = Fd(k)
          if (l .eq. k) then
            if (z .lt. 0.0_SP) then
              Fd(k) = -z
              do j = 1, n
                Fb(j,k) = -Fb(j,k)
              end do
            end if
            goto 3
          end if
          if (its .eq. MAX_ITS) then
            call output_line('Convergence failed!',&
                OU_CLASS_ERROR,OU_MODE_STD,'mprim_SVD_factoriseSP')
            call sys_halt()
          end if
          x = Fd(l)
          nm = k-1
          y = Fd(nm)
          g = rv1(nm)
          h = rv1(k)
          f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0_SP*h*y)
          g = sqrt(f*f+1.0_SP)
          f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x

          ! Next QR transformation
          c = 1.0_SP
          s = 1.0_SP
          do j = l, nm
            i = j+1
            g = rv1(i)
            y = Fd(i)
            h = s*g
            g = c*g
            z = sqrt(f*f+h*h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c)+(g*s)
            g = -(x*s)+(g*c)
            h = y*s
            y = y*c
            do jj = 1, n
              x = Fb(jj,j)
              z = Fb(jj,i)
              Fb(jj,j) = (x*c)+(z*s)
              Fb(jj,i) = -(x*s)+(z*c)
            end do
            z = sqrt(f*f+h*h)
            Fd(j) = z
            if (z .ne. 0.0_SP) then
              z = 1.0_SP/z
              c = f*z
              s = h*z
            end if
            f = (c*g)+(s*y)
            x = -(s*g)+(c*y)
            do jj = 1, m
              y = Fa(jj,j)
              z = Fa(jj,i)
              Fa(jj,j) = (y*c)+(z*s)
              Fa(jj,i) = -(y*s)+(z*c)
            end do
          end do
          rv1(l) = 0.0_SP
          rv1(k) = f
          Fd(k) = x
        end do
3       continue
      end do

    end if

    ! Deallocate internal memory
    deallocate(rv1)

  end subroutine mprim_SVD_factoriseSP

  !************************************************************************

!<subroutine>

  subroutine mprim_SVD_backsubstDP(Da,Dd,Db,Dx,Df,ndim1,ndim2,ndim3,btransposedOpt)

!<description>
    ! This subroutine solves <tex>$ A * x = f $</tex> for vector <tex>$ x $</tex>,
    ! where the rectangular
    ! matrix $A$ has been decomposed into <tex>$ Da $</tex>, <tex>$ Dd $</tex>
    ! and <tex>$ Db $</tex> by the routine
    ! mprim_SVD_factorise. The optional parameter btransposedOpt can be used
    ! to indicate that matrix A is stored in transposed format.
!</description>

!<input>
    ! Dimensions of the rectangular matrix
    integer, intent(in) :: ndim1,ndim2

    ! Dimension of the diagonal matrix Dd and the suare matrix Db
    integer, intent(in) :: ndim3

    ! OPTIONAL: Flag to indicate if the matrix A is transposed
    logical, intent(in), optional :: btransposedOpt

    ! Factorised matrix U
    real(DP), dimension(ndim1,ndim2), intent(in) :: Da

    ! Diagonal matrix D
    real(DP), dimension(ndim3), intent(in) :: Dd

    ! Square matrix B
    real(DP), dimension(ndim3,ndim3), intent(in) :: Db

    ! Right-hand side vector
    real(DP), dimension(ndim3), intent(in) :: Df
!</input>

!<output>
    ! Solution vector
    real(DP), dimension(ndim3), intent(out) :: Dx
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(:), allocatable :: Daux
    integer :: i,n,m
    logical :: btransposed

    ! Do we have an optional transposed flag?
    if (present(btransposedOpt)) then
      btransposed = btransposedOpt
    else
      btransposed = .false.
    end if

    ! Is the matrix transposed?
    if (btransposed) then
      n=ndim1; m=ndim2
    else
      n=ndim2; m=ndim1
    end if

    ! Check if number of equations is larger than number of unknowns
    if (m < n) then
      call output_line('Fewer equations than unknowns!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mprim_SVD_backsubstDP')
      call sys_halt()
    end if

    ! Allocate internal memory
    allocate(Daux(ndim3))

    ! Compute aux = (U^T * f)/D where D_i /= 0
    if (btransposed) then
      call DGEMV('n', ndim1, ndim2, 1.0_DP, Da, ndim1, Df, 1, 0.0_DP, Daux, 1)
    else
      call DGEMV('t', ndim1, ndim2, 1.0_DP, Da, ndim1, Df, 1, 0.0_DP, Daux, 1)
    end if

    do i =1, size(Dd)
      if (Dd(i) .ne. 0.0_DP) then
        Daux(i) = Daux(i)/Dd(i)
      else
        Daux(i) = 0.0_DP
      end if
    end do

    ! Compute x = B * aux
    call DGEMV('n', n, n, 1.0_DP, Db, n, Daux, 1, 0.0_DP, Dx, 1)

    ! Deallocate internal memory
    deallocate(Daux)
  end subroutine mprim_SVD_backsubstDP

  !************************************************************************

!<subroutine>

  subroutine mprim_SVD_backsubstSP(Fa,Fd,Fb,Fx,Ff,ndim1,ndim2,ndim3,btransposedOpt)

!<description>
    ! This subroutine solves <tex>$ A * x = f $</tex> for vector <tex>$ x $</tex>,
    ! where the rectangular
    ! matrix $A$ has been decomposed into <tex>$ Fa $</tex>, <tex>$ Fd $</tex>
    ! and <tex>$ Fb $</tex> by the routine
    ! mprim_SVD_factorise. The optional parameter btransposedOpt can be used
    ! to indicate that matrix A is stored in transposed format.
!</description>

!<input>
    ! Dimensions of the rectangular matrix
    integer, intent(in) :: ndim1,ndim2

    ! Dimension of the diagonal matrix Dd and the suare matrix Db
    integer, intent(in) :: ndim3

    ! OPTIONAL: Flag to indicate if the matrix A is transposed
    logical, intent(in), optional :: btransposedOpt

    ! Factorised matrix U
    real(SP), dimension(ndim1,ndim2), intent(in) :: Fa

    ! Diagonal matrix D
    real(SP), dimension(ndim3), intent(in) :: Fd

    ! Square matrix B
    real(SP), dimension(ndim3,ndim3), intent(in) :: Fb

    ! Right-hand side vector
    real(SP), dimension(ndim3), intent(in) :: Ff
!</input>

!<output>
    ! Solution vector
    real(SP), dimension(ndim3), intent(out) :: Fx
!</output>
!</subroutine>

    ! local variables
    real(SP), dimension(:), allocatable :: Faux
    integer :: i,n,m
    logical :: btransposed

    ! Do we have an optional transposed flag?
    if (present(btransposedOpt)) then
      btransposed = btransposedOpt
    else
      btransposed = .false.
    end if

    ! Is the matrix transposed?
    if (btransposed) then
      n=ndim1; m=ndim2
    else
      n=ndim2; m=ndim1
    end if

    ! Check if number of equations is larger than number of unknowns
    if (m < n) then
      call output_line('Fewer equations than unknowns!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mprim_SVD_backsubstSP')
      call sys_halt()
    end if

    ! Allocate internal memory
    allocate(Faux(ndim3))

    ! Compute aux = (U^T * f)/D where D_i /= 0
    if (btransposed) then
      call SGEMV('n', ndim1, ndim2, 1.0_SP, Fa, ndim1, Ff, 1, 0.0_SP, Faux, 1)
    else
      call SGEMV('t', ndim1, ndim2, 1.0_SP, Fa, ndim1, Ff, 1, 0.0_SP, Faux, 1)
    end if

    do i =1, size(Fd)
      if (Fd(i) .ne. 0.0_SP) then
        Faux(i) = Faux(i)/Fd(i)
      else
        Faux(i) = 0.0_SP
      end if
    end do

    ! Compute x = B * aux
    call SGEMV('n', n, n, 1.0_SP, Fb, n, Faux, 1, 0.0_SP, Fx, 1)

    ! Deallocate internal memory
    deallocate(Faux)
  end subroutine mprim_SVD_backsubstSP

  !************************************************************************

!<function>

  pure function mprim_stdDeviationDP(Dval) result(stdDev)

!<description>
    ! This function calculates the standard deviation of the given data
!</description>

!<input>
    ! double data
    real(DP), dimension(:), intent(in) :: Dval
!</input>

!<result>
    ! standard deviation
    real(DP) :: stdDev
!</result>
!</function>

    ! local variable
    real(DP)     :: mean
    integer(I32) :: i

    ! Compute mean
    mean = Dval(1)

    do i = 2, size(Dval)
      mean = mean + Dval(i)
    end do

    mean = mean/real(size(Dval), DP)

    ! Compute standard deviation
    stdDev = (Dval(1)-mean)*(Dval(1)-mean)

    do i = 2, size(Dval)
      stdDev = stdDev + (Dval(i)-mean)*(Dval(i)-mean)
    end do

    stdDev = sqrt(stdDev/real(size(Dval), DP))

  end function mprim_stdDeviationDP

  !************************************************************************

!<function>

  pure function mprim_stdDeviationSP(Fval) result(stdDev)

!<description>
    ! This function calculates the standard deviation of the given data
!</description>

!<input>
    ! single data
    real(SP), dimension(:), intent(in) :: Fval
!</input>

!<result>
    ! standard deviation
    real(SP) :: stdDev
!</result>
!</function>

    ! local variable
    real(SP)     :: mean
    integer(I32) :: i

    ! Compute mean
    mean = Fval(1)

    do i = 2, size(Fval)
      mean = mean + Fval(i)
    end do

    mean = mean/real(size(Fval), SP)

    ! Compute standard deviation
    stdDev = (Fval(1)-mean)*(Fval(1)-mean)

    do i = 2, size(Fval)
      stdDev = stdDev + (Fval(i)-mean)*(Fval(i)-mean)
    end do

    stdDev = sqrt(stdDev/real(size(Fval), SP))

  end function mprim_stdDeviationSP

  !************************************************************************

!<function>

  pure function mprim_stdDeviationInt(Ival) result(stdDev)

!<description>
    ! This function calculates the standard deviation of the given data
!</description>

!<input>
    ! integer data
    integer, dimension(:), intent(in) :: Ival
!</input>

!<result>
    ! standard deviation
    real(DP) :: stdDev
!</result>
!</function>

    ! local variable
    real(DP)     :: mean
    integer(I32) :: i

    ! Compute mean
    mean = real(Ival(1), DP)

    do i = 2, size(Ival)
      mean = mean + real(Ival(i), DP)
    end do

    mean = mean/real(size(Ival), DP)

    ! Compute standard deviation
    stdDev = (Ival(1)-mean)*(Ival(1)-mean)

    do i = 2, size(Ival)
      stdDev = stdDev + (Ival(i)-mean)*(Ival(i)-mean)
    end do

    stdDev = sqrt(stdDev/real(size(Ival), DP))

  end function mprim_stdDeviationInt

  !************************************************************************

!<function>

  pure function mprim_meanDeviationDP(Dval) result(meanDev)

!<description>
    ! This function calculates the mean deviation of the given data
!</description>

!<input>
    ! double data
    real(DP), dimension(:), intent(in) :: Dval
!</input>

!<result>
    ! mean deviation
    real(DP) :: meanDev
!</result>
!</function>

    ! local variable
    real(DP)     :: mean
    integer(I32) :: i

    ! Compute mean
    mean = Dval(1)

    do i = 2, size(Dval)
      mean = mean + Dval(i)
    end do

    mean = mean/real(size(Dval), DP)

    ! Compute mean deviation
    meanDev = abs(Dval(1)-mean)

    do i = 2, size(Dval)
      meanDev = meanDev + abs(Dval(i)-mean)
    end do

    meandev = meanDev/real(size(Dval), DP)

  end function mprim_meanDeviationDP

  !************************************************************************

!<function>

  pure function mprim_meanDeviationSP(Fval) result(meanDev)

!<description>
    ! This function calculates the mean deviation of the given data
!</description>

!<input>
    ! single data
    real(SP), dimension(:), intent(in) :: Fval
!</input>

!<result>
    ! mean deviation
    real(SP) :: meanDev
!</result>
!</function>

    ! local variable
    real(SP)     :: mean
    integer(I32) :: i

    ! Compute mean
    mean = Fval(1)

    do i = 2, size(Fval)
      mean = mean + Fval(i)
    end do

    mean = mean/real(size(Fval), SP)

    ! Compute mean deviation
    meanDev = abs(Fval(1)-mean)

    do i = 2, size(Fval)
      meanDev = meanDev + abs(Fval(i)-mean)
    end do

    meanDev = meanDev/real(size(Fval), SP)

  end function mprim_meanDeviationSP

  !************************************************************************

!<function>

  pure function mprim_meanDeviationInt(Ival) result(meanDev)

!<description>
    ! This function calculates the mean deviation of the given data
!</description>

!<input>
    ! integer data
    integer, dimension(:), intent(in) :: Ival
!</input>

!<result>
    ! mean deviation
    real(DP) :: meanDev
!</result>
!</function>

    ! local variable
    real(DP)     :: mean
    integer(I32) :: i

    ! Compute mean
    mean = real(Ival(1), DP)

    do i = 2, size(Ival)
      mean = mean + real(Ival(i), DP)
    end do

    mean = mean/real(size(Ival), DP)

    ! Compute mean deviation
    meanDev = abs(Ival(1)-mean)

    do i = 2, size(Ival)
      meanDev = meanDev + abs(Ival(i)-mean)
    end do

    meanDev = meanDev/real(size(Ival), DP)

  end function mprim_meanDeviationInt

  !************************************************************************

!<function>

  pure function mprim_meanValueDP(Dval) result(meanVal)

!<description>
    ! This function calculates the mean value of the given data
!</description>

!<input>
    ! double data
    real(DP), dimension(:), intent(in) :: Dval
!</input>

!<result>
    ! mean value
    real(DP) :: meanVal
!</result>
!</function>

    ! local variable
    integer(I32) :: i

    ! Compute mean value
    meanVal = Dval(1)

    do i = 2, size(Dval)
      meanVal = meanVal + Dval(i)
    end do

    meanVal = meanVal/real(size(Dval), DP)
  end function mprim_meanValueDP

  !************************************************************************

!<function>

  pure function mprim_meanValueSP(Fval) result(meanVal)

!<description>
    ! This function calculates the mean value of the given data
!</description>

!<input>
    ! single data
    real(SP), dimension(:), intent(in) :: Fval
!</input>

!<result>
    ! mean value
    real(SP) :: meanVal
!</result>
!</function>

    ! local variable
    integer(I32) :: i

    ! Compute mean value
    meanVal = Fval(1)

    do i = 2, size(Fval)
      meanVal = meanVal + Fval(i)
    end do

    meanVal = meanVal/real(size(Fval), SP)
  end function mprim_meanValueSP

  !************************************************************************

!<function>

  pure function mprim_meanValueInt(Ival) result(meanVal)

!<description>
    ! This function calculates the mean value of the given data
!</description>

!<input>
    ! integer data
    integer, dimension(:), intent(in) :: Ival
!</input>

!<result>
    ! mean value
    integer :: meanVal
!</result>
!</function>

    ! local variable
    integer(I32) :: i

    ! Compute mean value
    meanVal = Ival(1)

    do i = 2, size(Ival)
      meanVal = meanVal + Ival(i)
    end do

    meanVal = meanVal/size(Ival)
  end function mprim_meanValueInt

  ! *****************************************************************************

!<function>

  elemental function mprim_degToRad(d) result(r)

!<description>
    ! This function converts DEG to RAD
!</description>

!<input>
    ! DEG
    real(DP), intent(in) :: d
!</input>

!<result>
    ! RAD
    real(DP) :: r
!</result>
!</function>

    r = d * (SYS_PI / 180._DP)
  end function mprim_degToRad

  ! *****************************************************************************

!<function>

  elemental function mprim_radToDeg(r) result(d)

!<description>
    ! This function converts RAD to DEG
!</description>

!<input>
    ! RAD
    real(DP), intent(in) :: r
!</input>

!<result>
    ! DEG
    real(DP) :: d
!</result>
!</function>

    d = r * (180._DP / SYS_PI)

  end function mprim_radToDeg

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_solve2x2DirectDP(Da,Db)

!<description>
  ! This subroutine directly solves a 2x2 system without any pivoting.
  ! 'Da' is a 2x2 matrix, Db a 2-tupel with the right hand
  ! side which is replaced by the solution.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Da and Db are assumed to be 2x2 arrays!
!</description>

!<input>
  ! Matrix A.
  real(DP), dimension(2,2), intent(in) :: Da
!</input>

!<inputoutput>
  ! On entry: RHS vector b.
  ! On exit: Solution x with Ax=b.
  real(DP), dimension(2), intent(out) :: Db
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP) :: ddet,x1,x2

    ! Explicit formula for 2x2 systems, computed with Maple.
    ddet = 1.0_DP / (Da(1,1)*Da(2,2)-Da(1,2)*Da(2,1))
    x1 = Db(1)
    x2 = Db(2)
    Db(1) = -(Da(1,2)*x2 - Da(2,2)*x1) * ddet
    Db(2) =  (Da(1,1)*x2 - Da(2,1)*x1) * ddet

  end subroutine mprim_solve2x2DirectDP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_solve2x2DirectSP(Fa,Fb)

!<description>
  ! This subroutine directly solves a 2x2 system without any pivoting.
  ! 'Fa' is a 2x2 matrix, Fb a 2-tupel with the right hand
  ! side which is replaced by the solution.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Fa and Fb are assumed to be 2x2 arrays!
!</description>

!<input>
  ! Matrix A.
  real(SP), dimension(2,2), intent(in) :: Fa
!</input>

!<inputoutput>
  ! On entry: RHS vector b.
  ! On exit: Solution x with Ax=b.
  real(SP), dimension(2), intent(out) :: Fb
!</inputoutput>

!</subroutine>

    ! local variables
    real(SP) :: fdet,x1,x2

    ! Explicit formula for 2x2 systems, computed with Maple.
    fdet = 1.0_SP / (Fa(1,1)*Fa(2,2)-Fa(1,2)*Fa(2,1))
    x1 = Fb(1)
    x2 = Fb(2)
    Fb(1) = -(Fa(1,2)*x2 - Fa(2,2)*x1) * fdet
    Fb(2) =  (Fa(1,1)*x2 - Fa(2,1)*x1) * fdet

  end subroutine mprim_solve2x2DirectSP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_solve3x3DirectDP(Da,Db)

!<description>
  ! This subroutine directly solves a 3x3 system without any pivoting.
  ! 'Da' is a 3x3 matrix, Db a 3-tupel with the right hand
  ! side which is replaced by the solution.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Da and Db are assumed to be 2x2 arrays!
!</description>

!<input>
  ! Matrix A.
  real(DP), dimension(3,3), intent(in) :: Da
!</input>

!<inputoutput>
  ! On entry: RHS vector b.
  ! On exit: Solution x with Ax=b.
  real(DP), dimension(3), intent(inout) :: Db
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP) :: ddet,x1,x2,x3

    ! Explicit formula for 3x3 systems, computed with Maple.
    ddet = 1.0_DP / &
        (Da(1,1)*Da(2,2)*Da(3,3) &
        -Da(1,1)*Da(3,2)*Da(2,3) &
        -Da(2,1)*Da(1,2)*Da(3,3) &
        +Da(3,2)*Da(2,1)*Da(1,3) &
        -Da(2,2)*Da(3,1)*Da(1,3) &
        +Da(3,1)*Da(1,2)*Da(2,3))
    x1 = Db(1)
    x2 = Db(2)
    x3 = Db(3)
    Db(1) = ddet * &
        (Da(1,2)*Da(2,3)*x3 &
        -Da(1,2)*x2*Da(3,3) &
        +Da(1,3)*Da(3,2)*x2 &
        -Da(1,3)*Da(2,2)*x3 &
        +x1*Da(2,2)*Da(3,3) &
        -x1*Da(3,2)*Da(2,3))
    Db(2) = - ddet * &
        (Da(1,1)*Da(2,3)*x3 &
        -Da(1,1)*x2*Da(3,3) &
        -Da(2,1)*Da(1,3)*x3 &
        -Da(2,3)*Da(3,1)*x1 &
        +x2*Da(3,1)*Da(1,3) &
        +Da(2,1)*x1*Da(3,3))
    Db(3) = ddet * &
        (Da(3,2)*Da(2,1)*x1 &
        -Da(1,1)*Da(3,2)*x2 &
        +Da(1,1)*Da(2,2)*x3 &
        -Da(2,2)*Da(3,1)*x1 &
        -Da(2,1)*Da(1,2)*x3 &
        +Da(3,1)*Da(1,2)*x2)

  end subroutine mprim_solve3x3DirectDP

  ! ***************************************************************************

!<subroutine>

  pure subroutine mprim_solve3x3DirectSP(Fa,Fb)

!<description>
  ! This subroutine directly solves a 3x3 system without any pivoting.
  ! 'Fa' is a 3x3 matrix, Fb a 3-tupel with the right hand
  ! side which is replaced by the solution.
  !
  ! Warning: For speed reasons, there is no array bounds checking
  ! activated in this routine! Fa and Fb are assumed to be 2x2 arrays!
!</description>

!<input>
  ! Matrix A.
  real(SP), dimension(3,3), intent(in) :: Fa
!</input>

!<inputoutput>
  ! On entry: RHS vector b.
  ! On exit: Solution x with Ax=b.
  real(SP), dimension(3), intent(inout) :: Fb
!</inputoutput>

!</subroutine>

    ! local variables
    real(SP) :: fdet,x1,x2,x3

    ! Explicit formula for 3x3 systems, computed with Maple.
    fdet = 1.0_SP / &
        (Fa(1,1)*Fa(2,2)*Fa(3,3) &
        -Fa(1,1)*Fa(3,2)*Fa(2,3) &
        -Fa(2,1)*Fa(1,2)*Fa(3,3) &
        +Fa(3,2)*Fa(2,1)*Fa(1,3) &
        -Fa(2,2)*Fa(3,1)*Fa(1,3) &
        +Fa(3,1)*Fa(1,2)*Fa(2,3))
    x1 = Fb(1)
    x2 = Fb(2)
    x3 = Fb(3)
    Fb(1) = fdet * &
        (Fa(1,2)*Fa(2,3)*x3 &
        -Fa(1,2)*x2*Fa(3,3) &
        +Fa(1,3)*Fa(3,2)*x2 &
        -Fa(1,3)*Fa(2,2)*x3 &
        +x1*Fa(2,2)*Fa(3,3) &
        -x1*Fa(3,2)*Fa(2,3))
    Fb(2) = - fdet * &
        (Fa(1,1)*Fa(2,3)*x3 &
        -Fa(1,1)*x2*Fa(3,3) &
        -Fa(2,1)*Fa(1,3)*x3 &
        -Fa(2,3)*Fa(3,1)*x1 &
        +x2*Fa(3,1)*Fa(1,3) &
        +Fa(2,1)*x1*Fa(3,3))
    Fb(3) = fdet * &
        (Fa(3,2)*Fa(2,1)*x1 &
        -Fa(1,1)*Fa(3,2)*x2 &
        +Fa(1,1)*Fa(2,2)*x3 &
        -Fa(2,2)*Fa(3,1)*x1 &
        -Fa(2,1)*Fa(1,2)*x3 &
        +Fa(3,1)*Fa(1,2)*x2)

  end subroutine mprim_solve3x3DirectSP

  ! ************************************************************************

!<subroutine>

  pure subroutine mprim_solve2x2BandDiagDP (neqA,Da,Db,Dd,Dc,Dvec1,Dvec2)

!<description>
  ! This routine solves to a 2x2 block matrix with all blocks consisting
  ! of only diagonal bands.
  !
  ! <!--
  !
  ! The system that is to be solved here is assumed to have the following shape:
  !
  !   (  A             B             ) ( U1 ) = ( F1 )
  !   (       ..            ..       ) ( U1 )   ( F1 )
  !   (             A             B  ) ( U1 )   ( F1 )
  !   (  D             C             ) ( U2 )   ( F2 )
  !   (       ..            ..       ) ( U2 )   ( F2 )
  !   (             D             C  ) ( U2 )   ( F2 )
  !
  ! or in short:
  !
  !   ( A B ) = (U) = (F)
  !   ( D C )   (U)   (F)
  !
  ! with all matrices diagonal matrices.
  !
  ! -->
!</description>

!<input>
  ! Dimension of the matrix blocks
  integer, intent(in) :: neqA

  ! Submatrix A, only diagonal entries.
  real(DP), dimension(*), intent(in) :: Da

  ! Diagonal entries in the submatrix B.
  real(DP), dimension(*), intent(in) :: Db

  ! Diagonal entries in the submatrix D.
  real(DP), dimension(*), intent(in) :: Dd

  ! Diagonal elements of the local system matrix C
  real(DP), dimension(*), intent(in) :: Dc
!</input>

!<inputoutput>
  ! On entry: Right hand side F1.
  ! On exit: Solution vector U1.
  real(DP), dimension(*), intent(inout) :: Dvec1

  ! On entry: Right hand side F2.
  ! On exit: Solution vector U2.
  real(DP), dimension(*), intent(inout) :: Dvec2
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    real(DP) :: ddet,a11,a12,a21,a22,x1,x2

    ! Such a system can be solved by reducing it to neqA*neqA 2x2 systems
    ! as there is only minimal coupling present between the entries.
    !
    ! Loop over the entries of the A matrix.
    do i=1,neqA

      ! Fetch a 2x2 system from the big matrix
      a11 = Da(i)
      a21 = Dd(i)
      a12 = Db(i)
      a22 = Dc(i)
      x1 = Dvec1(i)
      x2 = Dvec2(i)

      ! Solve the system, overwrite the input vector.
      ddet = 1.0_DP / (a11*a22 - a12*a21)
      Dvec1(i) = -(a12*x2 - a22*x1) * ddet
      Dvec2(i) =  (a11*x2 - a21*x1) * ddet

    end do

  end subroutine mprim_solve2x2BandDiagDP

  ! ************************************************************************

!<subroutine>

  pure subroutine mprim_solve2x2BandDiagSP (neqA,Fa,Fb,Dd,Fc,Fvec1,Fvec2)

!<description>
  ! This routine solves to a 2x2 block matrix with all blocks consisting
  ! of only diagonal bands.
  !
  ! <!--
  !
  ! The system that is to be solved here is assumed to have the following shape:
  !
  !   (  A             B             ) ( U1 ) = ( F1 )
  !   (       ..            ..       ) ( U1 )   ( F1 )
  !   (             A             B  ) ( U1 )   ( F1 )
  !   (  D             C             ) ( U2 )   ( F2 )
  !   (       ..            ..       ) ( U2 )   ( F2 )
  !   (             D             C  ) ( U2 )   ( F2 )
  !
  ! or in short:
  !
  !   ( A B ) = (U) = (F)
  !   ( D C )   (U)   (F)
  !
  ! with all matrices diagonal matrices.
  !
  ! -->
!</description>

!<input>
  ! Dimension of the matrix blocks
  integer, intent(in) :: neqA

  ! Submatrix A, only diagonal entries.
  real(SP), dimension(*), intent(in) :: Fa

  ! Diagonal entries in the submatrix B.
  real(SP), dimension(*), intent(in) :: Fb

  ! Diagonal entries in the submatrix D.
  real(SP), dimension(*), intent(in) :: Dd

  ! Diagonal elements of the local system matrix C
  real(SP), dimension(*), intent(in) :: Fc
!</input>

!<inputoutput>
  ! On entry: Right hand side F1.
  ! On exit: Solution vector U1.
  real(SP), dimension(*), intent(inout) :: Fvec1

  ! On entry: Right hand side F2.
  ! On exit: Solution vector U2.
  real(SP), dimension(*), intent(inout) :: Fvec2
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    real(SP) :: fdet,a11,a12,a21,a22,x1,x2

    ! Such a system can be solved by reducing it to neqA*neqA 2x2 systems
    ! as there is only minimal coupling present between the entries.
    !
    ! Loop over the entries of the A matrix.
    do i=1,neqA

      ! Fetch a 2x2 system from the big matrix
      a11 = Fa(i)
      a21 = Dd(i)
      a12 = Fb(i)
      a22 = Fc(i)
      x1 = Fvec1(i)
      x2 = Fvec2(i)

      ! Solve the system, overwrite the input vector.
      fdet = 1.0_SP / (a11*a22 - a12*a21)
      Fvec1(i) = -(a12*x2 - a22*x1) * fdet
      Fvec2(i) =  (a11*x2 - a21*x1) * fdet

    end do

  end subroutine mprim_solve2x2BandDiagSP

  !*****************************************************************************

!<function>

  elemental function mprim_minmod2DP(a,b) result (c)

!<description>
    ! The minmod functions returns zero if the two arguments a and b
    ! have different sign and the argument with the smallest absolute
    ! value otherwise.
!</description>

!<input>
    real(DP), intent(in) :: a,b
!</input>

!<result>
    real(DP) :: c
!</result>
!</function>

    if (a*b .le. 0.0_DP) then
      c = 0.0_DP
    else
      c = sign(min(abs(a), abs(b)), a)
    end if
  end function

  !*****************************************************************************

!<function>

  elemental function mprim_minmod2SP(a,b) result (c)

!<description>
    ! The minmod functions returns zero if the two arguments a and b
    ! have different sign and the argument with the smallest absolute
    ! value otherwise.
!</description>

!<input>
    real(SP), intent(in) :: a,b
!</input>

!<result>
    real(SP) :: c
!</result>
!</function>

    if (a*b .le. 0.0_DP) then
      c = 0.0_DP
    else
      c = sign(min(abs(a), abs(b)), a)
    end if
  end function

  !*****************************************************************************

!<function>

  elemental function mprim_minmod3DP(a,b,c) result (d)

!<description>
    ! The minmod functions returns zero if the two arguments a and b
    ! have different sign and the scaling parameter d by which the
    ! third argument c has to be scaled to obtain minmod(a,b) otherwise.
!</description>

!<input>
    real(DP), intent(in) :: a,b,c
!</input>

!<result>
    real(DP) :: d
!</result>
!</function>

    if (a*b*abs(c) .le. 0.0_DP) then
      d = 0.0_DP
    elseif (abs(a) .lt. abs(b)) then
      d = a/c
    else
      d = b/c
    end if
  end function

  !*****************************************************************************

!<function>

  elemental function mprim_minmod3SP(a,b,c) result (d)

!<description>
    ! The minmod functions returns zero if the two arguments a and b
    ! have different sign and the scaling parameter d by which the
    ! third argument c hass to be scaled to obtain minmod(a,b) otherwise.
!</description>

!<input>
    real(SP), intent(in) :: a,b,c
!</input>

!<result>
    real(SP) :: d
!</result>
!</function>

    if (a*b*abs(c) .le. 0.0_DP) then
      d = 0.0_DP
    elseif (abs(a) .lt. abs(b)) then
      d = a/c
    else
      d = b/c
    end if
  end function

  !*****************************************************************************

!<function>

  elemental function mprim_softmax2DP(a,b) result (c)

!<description>
    ! The softmax function computes the maximum value of the given
    ! data a and b but with smooth transition, that is, it does not
    ! switch abruptly between values a and b. This implementation is
    ! based on the idea of John D. Cook, c.f.
    !   http://www.johndcook.com/blog/2010/01/13/soft-maximum/
    ! This implementation is taken from
    !   http://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/
!</description>

!<input>
    real(DP), intent(in) :: a,b
!</input>

!<result>
    real(DP) :: c
!</result>
!</function>

    ! local variables
    real(DP) :: dminval,dmaxval

    dminval = min(a,b)
    dmaxval = max(a,b)

    c = dmaxval + log(1.0_DP+exp(dminval-dmaxval))

  end function

  !*****************************************************************************

!<function>

  elemental function mprim_softmax3DP(a,b,c) result (d)

!<description>
    ! The softmax function computes the maximum value of the given
    ! data a and b but with smooth transition, that is, it does not
    ! switch abruptly between values a and b. This implementation is
    ! based on the idea of John D. Cook, c.f.
    !   http://www.johndcook.com/blog/2010/01/13/soft-maximum/
    ! The third parameter c is used to control the sharpness of the
    ! soft-maximum function. This implementation is taken from
    !   http://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/
!</description>

!<input>
    real(DP), intent(in) :: a,b,c
!</input>

!<result>
    real(DP) :: d
!</result>
!</function>

    ! local variables
    real(DP) :: dminval,dmaxval

    dminval = min(a,b)
    dmaxval = max(a,b)

    d = dmaxval + log(1.0_DP+exp(c*(dminval-dmaxval)))/c

  end function

  !*****************************************************************************

!<function>

  elemental function mprim_softmax2SP(a,b) result (c)

!<description>
    ! The softmax function computes the maximum value of the given
    ! data a and b but with smooth transition, that is, it does not
    ! switch abruptly between values a and b. This implementation is
    ! based on the idea of John D. Cook, c.f.
    !   http://www.johndcook.com/blog/2010/01/13/soft-maximum/
    ! This implementation is taken from
    !   http://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/
!</description>

!<input>
    real(SP), intent(in) :: a,b
!</input>

!<result>
    real(SP) :: c
!</result>
!</function>

    ! local variables
    real(SP) :: fminval,fmaxval

    fminval = min(a,b)
    fmaxval = max(a,b)

    c = fmaxval + log(1.0_SP+exp(fminval-fmaxval))

  end function

  !*****************************************************************************

!<function>

  elemental function mprim_softmax3SP(a,b,c) result (d)

!<description>
    ! The softmax function computes the maximum value of the given
    ! data a and b but with smooth transition, that is, it does not
    ! switch abruptly between values a and b. This implementation is
    ! based on the idea of John D. Cook, c.f.
    !   http://www.johndcook.com/blog/2010/01/13/soft-maximum/
    ! The third parameter c is used to control the sharpness of the
    ! soft-maximum function. This implementation is taken from
    !   http://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/
!</description>

!<input>
    real(SP), intent(in) :: a,b,c
!</input>

!<result>
    real(SP) :: d
!</result>
!</function>

    ! local variables
    real(SP) :: fminval,fmaxval

    fminval = min(a,b)
    fmaxval = max(a,b)

    d = fmaxval + log(1.0_SP+exp(c*(fminval-fmaxval)))/c

  end function

!*****************************************************************************

!<function>

  elemental function mprim_softmin2DP(a,b) result (c)

!<description>
    ! The softmin function computes the minimum value of the given
    ! data a and b but with smooth transition, that is, it does not
    ! switch abruptly between values a and b. This implementation is
    ! based on the idea of John D. Cook, c.f.
    !   http://www.johndcook.com/blog/2010/01/13/soft-maximum/
    ! This implementation is taken from
    !   http://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/
!</description>

!<input>
    real(DP), intent(in) :: a,b
!</input>

!<result>
    real(DP) :: c
!</result>
!</function>

    ! local variables
    real(DP) :: dminval,dmaxval

    dminval = min(a,b)
    dmaxval = max(a,b)

    c = dminval - log(1.0_DP+exp(dminval-dmaxval))

  end function

  !*****************************************************************************

!<function>

  elemental function mprim_softmin3DP(a,b,c) result (d)

!<description>
    ! The softmin function computes the minimum value of the given
    ! data a and b but with smooth transition, that is, it does not
    ! switch abruptly between values a and b. This implementation is
    ! based on the idea of John D. Cook, c.f.
    !   http://www.johndcook.com/blog/2010/01/13/soft-maximum/
    ! The third parameter c is used to control the sharpness of the
    ! soft-maximum function. This implementation is taken from
    !   http://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/
!</description>

!<input>
    real(DP), intent(in) :: a,b,c
!</input>

!<result>
    real(DP) :: d
!</result>
!</function>

    ! local variables
    real(DP) :: dminval,dmaxval

    dminval = min(a,b)
    dmaxval = max(a,b)

    d = dminval - log(1.0_DP+exp(c*(dminval-dmaxval)))/c

  end function

  !*****************************************************************************

!<function>

  elemental function mprim_softmin2SP(a,b) result (c)

!<description>
    ! The softmin function computes the minimum value of the given
    ! data a and b but with smooth transition, that is, it does not
    ! switch abruptly between values a and b. This implementation is
    ! based on the idea of John D. Cook, c.f.
    !   http://www.johndcook.com/blog/2010/01/13/soft-maximum/
    ! This implementation is taken from
    !   http://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/
!</description>

!<input>
    real(SP), intent(in) :: a,b
!</input>

!<result>
    real(SP) :: c
!</result>
!</function>

    ! local variables
    real(SP) :: fminval,fmaxval

    fminval = min(a,b)
    fmaxval = max(a,b)

    c = fminval - log(1.0_SP+exp(fminval-fmaxval))

  end function

  !*****************************************************************************

!<function>

  elemental function mprim_softmin3SP(a,b,c) result (d)

!<description>
    ! The softmin function computes the minimum value of the given
    ! data a and b but with smooth transition, that is, it does not
    ! switch abruptly between values a and b. This implementation is
    ! based on the idea of John D. Cook, c.f.
    !   http://www.johndcook.com/blog/2010/01/13/soft-maximum/
    ! The third parameter c is used to control the sharpness of the
    ! soft-maximum function. This implementation is taken from
    !   http://www.johndcook.com/blog/2010/01/20/how-to-compute-the-soft-maximum/
!</description>

!<input>
    real(SP), intent(in) :: a,b,c
!</input>

!<result>
    real(SP) :: d
!</result>
!</function>

    ! local variables
    real(SP) :: fminval,fmaxval

    fminval = min(a,b)
    fmaxval = max(a,b)

    d = fminval - log(1.0_SP+exp(c*(fminval-fmaxval)))/c

  end function

  !*****************************************************************************

!<subroutine>

  subroutine mprim_leastSquaresMinDP(Da,Db,Dx,ndim1,ndim2,ndim3,ndim4,&
      btransposedOpt,drcondOpt)

!<description>
    ! This subroutine solves the least-squares minimisation problem
    ! $$ min_x \|Ax-b \|_2 $$
    ! for given matrix $A$ and vector $b$. If matrix A has full rank
    ! the least-squares minimisation problem is solved using QR or LQ
    ! factorisation. If this approach fails then the least-squares
    ! minimisation problem is solved using a complete orthorgonal
    ! factorisation of matrix A which works also for rank-deficient
    ! matrices.
!</description>

!<input>
    ! Dimension of matrix
    integer, intent(in) :: ndim1,ndim2

    ! Dimension of right-hand side vector
    integer, intent(in) :: ndim3

    ! Dimension of solution vector
    integer, intent(in) :: ndim4

    ! Matrix for the least-squares minimisation
    real(DP), dimension(ndim1,ndim2), intent(in) :: Da

    ! Vector for the least-squares minimisation
    real(DP), dimension(ndim3), intent(in) :: Db

    ! OPTIONAL: Flag which indicates whether matrix A or its transpose
    ! should be used in the least-squares minimisation problem
    logical, intent(in), optional :: btransposedOpt

    ! OPTIONAL: Tolerance to determine the effective rank of matrix
    ! Da, which is defined as the order of the largest leading
    ! triangular submatrix in the QR factorisation with pivoting of A
    ! whose estimated condition number < 1/drcondOpt.
    real(DP), intent(in), optional :: drcondOpt
!</input>

!<output>
    ! Solution vector from the least-squares minimisation
    real(DP), dimension(ndim4), intent(out) :: Dx
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_Da
    real(DP), dimension(:), pointer :: p_Db,p_Dwork
    integer, dimension(:), pointer :: p_Ipvt
    integer :: nwork,naux,info,irank,n,m
    logical :: btransposed

    ! Do we have an optional transposed flag?
    if (present(btransposedOpt)) then
      btransposed = btransposedOpt
    else
      btransposed = .false.
    end if

    ! Determine size of working memory
    nwork = ndim1+ndim2

    ! Allocate working memory
    allocate(p_Dwork(nwork))

    ! Make a copy of the input matrix A which will be overwritten,
    ! e.g., be the QR or LQ factorisation
    allocate(p_Da(ndim1,ndim2))
    call lalg_copyVector(Da,p_Da,ndim1,ndim2)

    ! Determine size of source/destination vector
    naux = max(ndim3,ndim4)

    ! Make a copy of the input vector b which will be overwritten be
    ! the solution vector
    allocate(p_Db(naux))
    call lalg_copyVector(Db,p_Db,ndim3)
    
    ! Solve the least-squares minimisation problem using QR or LQ
    ! factorisation. It is assumed that matrix A has full rank.
    call dgels(merge('T','N',btransposed), ndim1, ndim2, 1, p_Da, ndim1,&
               p_Db, naux, p_Dwork, nwork, info)

    ! Check if solution procedure terminated without errors
    if (info .eq. 0) then

      ! Copy solution vector to output array
      call lalg_copyVector(p_Db,Dx,ndim4)

      ! Release internal memory
      deallocate(p_Da,p_Db,p_Dwork)

      ! That is it
      return
    end if
    
    ! Otherwise, ...
    if (info .lt. 0) then
      call output_line('Invalid parameter!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mprim_leastSquareMinDP')
      call sys_halt()
    else
      
      ! Matrix A does not have full rank. In this case, we need to
      ! solve the least-squares minimisation problem using SVD.

      ! Release internal memory
      deallocate(p_Dwork)

      ! Again, make a copy of the input vector b which will be
      ! overwritten be the solution vector
      call lalg_copyVector(Db,p_Db,ndim3)

      ! Transposition of matrix A cannot be handled implicitly. In this
      ! case, we must transpose matrix A by hand before performing SVD.
      if (btransposed) then

        n = ndim1
        m = ndim2

        ! Transpose matrix A by hand
        deallocate(p_Da)
        allocate(p_Da(m,n))
        call mprim_transposeMatrix(Da,p_Da)

      else

        m = ndim1
        n = ndim2

        ! Make a copy of the input matrix A which will be overwritten.
        call lalg_copyVector(Da,p_Da,ndim1,ndim2)

      end if

      ! Allocate internal memory
      allocate(p_Ipvt(n))
      call lalg_clearVector(p_Ipvt)
      
      ! Determine size of working memory
      naux  = min(m,n)
      nwork = max( naux+3*n+1, 2*naux+1 )
      allocate(p_Dwork(nwork))
      
      ! Do we have an optional tolerance paramter?
      if (present(drcondOpt)) then
        ! Solve the least-squares minimisation problem using SVD.
        call dgelsy(m, n, 1, p_Da, m, p_Db, m, p_Ipvt,&
                    drcondOpt, irank, p_Dwork, nwork, info)
      else
        ! Solve the least-squares minimisation problem using SVD.
        call dgelsy(m, n, 1, p_Da, m, p_Db, m, p_Ipvt,&
                    1e-12_DP, irank, p_Dwork, nwork, info)
      end if

      ! Check if solution procedure terminated without errors
      if (info .eq. 0) then
        
        ! Copy solution vector to output array
        call lalg_copyVector(p_Db,Dx,ndim4)
        
        ! Release internal memory
        deallocate(p_Da,p_Db,p_Dwork,p_Ipvt)
        
        ! That is it
        return
      else
        call output_line('Unable to solve least-squares minimisation problem!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mprim_leastSquareMinDP')
        call sys_halt()       
      end if
      
    end if

  end subroutine mprim_leastSquaresMinDP

  !*****************************************************************************

!<subroutine>

  subroutine mprim_leastSquaresMinSP(Fa,Fb,Fx,ndim1,ndim2,ndim3,ndim4,&
      btransposedOpt,frcondOpt)

!<description>
    ! This subroutine solves the least-squares minimisation problem
    ! $$ min_x \|Ax-b \|_2 $$
    ! for given matrix $A$ and vector $b$. If matrix A has full rank
    ! the least-squares minimisation problem is solved using QR or LQ
    ! factorisation. If this approach fails then the least-squares
    ! minimisation problem is solved using a complete orthorgonal
    ! factorisation of matrix A which works also for rank-deficient
    ! matrices.
!</description>

!<input>
    ! Dimension of matrix and vectors
    integer, intent(in) :: ndim1,ndim2

    ! Dimension of right-hand side vector
    integer, intent(in) :: ndim3

    ! Dimension of solution vector
    integer, intent(in) :: ndim4

    ! Matrix for the least-squares minimisation
    real(SP), dimension(ndim1,ndim2), intent(in) :: Fa

    ! Vector for the least-squares minimisation
    real(SP), dimension(ndim3), intent(in) :: Fb

    ! OPTIONAL: Flag which indicates whether matrix A or its transpose
    ! should be used in the least-squares minimisation problem
    logical, intent(in), optional :: btransposedOpt

    ! OPTIONAL: Tolerance to determine the effective rank of matrix
    ! Da, which is defined as the order of the largest leading
    ! triangular submatrix in the QR factorisation with pivoting of A
    ! whose estimated condition number < 1/frcondOpt.
    real(SP), intent(in), optional :: frcondOpt
!</input>

!<output>
    ! Solution vector from the least-squares minimisation
    real(SP), dimension(ndim4), intent(out) :: Fx
!</output>

!</subroutine>

    ! local variables
    real(SP), dimension(:,:), pointer :: p_Fa
    real(SP), dimension(:), pointer :: p_Fb,p_Fwork
    integer, dimension(:), pointer :: p_Ipvt
    integer :: nwork,naux,info,irank,n,m
    logical :: btransposed

    ! Do we have an optional transposed flag?
    if (present(btransposedOpt)) then
      btransposed = btransposedOpt
    else
      btransposed = .false.
    end if

    ! Determine size of working memory
    nwork = ndim1+ndim2

    ! Allocate working memory
    allocate(p_Fwork(nwork))

    ! Make a copy of the input matrix A which will be overwritten,
    ! e.g., be the QR or LQ factorisation
    allocate(p_Fa(ndim1,ndim2))
    call lalg_copyVector(Fa,p_Fa,ndim1,ndim2)

    ! Determine size of source/destination vector
    naux = max(ndim3,ndim4)

    ! Make a copy of the input vector b which will be overwritten be
    ! the solution vector
    allocate(p_Fb(naux))
    call lalg_copyVector(Fb,p_Fb,ndim3)
    
    ! Solve the least-squares minimisation problem using QR or LQ
    ! factorisation. It is assumed that matrix A has full rank.
    call dgels(merge('T','N',btransposed), ndim1, ndim2, 1, p_Fa, ndim1,&
               p_Fb, naux, p_Fwork, nwork, info)

    ! Check if solution procedure terminated without errors
    if (info .eq. 0) then

      ! Copy solution vector to output array
      call lalg_copyVector(p_Fb,Fx,ndim4)

      ! Release internal memory
      deallocate(p_Fa,p_Fb,p_Fwork)

      ! That is it
      return
    end if
    
    ! Otherwise, ...
    if (info .lt. 0) then
      call output_line('Invalid parameter!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mprim_leastSquareMinSP')
      call sys_halt()
    else
      
      ! Matrix A does not have full rank. In this case, we need to
      ! solve the least-squares minimisation problem using SVD.

      ! Release internal memory
      deallocate(p_Fwork)

      ! Again, make a copy of the input vector b which will be
      ! overwritten be the solution vector
      call lalg_copyVector(Fb,p_Fb,ndim3)

      ! Transposition of matrix A cannot be handled implicitly. In this
      ! case, we must transpose matrix A by hand before performing SVD.
      if (btransposed) then

        n = ndim1
        m = ndim2

        ! Transpose matrix A by hand
        deallocate(p_Fa)
        allocate(p_Fa(m,n))
        call mprim_transposeMatrix(Fa,p_Fa)

      else

        m = ndim1
        n = ndim2

        ! Make a copy of the input matrix A which will be overwritten.
        call lalg_copyVector(Fa,p_Fa,ndim1,ndim2)

      end if

      ! Allocate internal memory
      allocate(p_Ipvt(n))
      call lalg_clearVector(p_Ipvt)
      
      ! Determine size of working memory
      naux = min(m,n)
      nwork = max( naux+3*n+1, 2*naux+1 )
      allocate(p_Fwork(nwork))

      ! Do we have an optional tolerance paramter?
      if (present(frcondOpt)) then
        ! Solve the least-squares minimisation problem using SVD.
        call dgelsy(m, n, 1, p_Fa, m, p_Fb, m, p_Ipvt,&
                    frcondOpt, irank, p_Fwork, nwork, info)
      else
        ! Solve the least-squares minimisation problem using SVD.
        call dgelsy(m, n, 1, p_Fa, m, p_Fb, m, p_Ipvt,&
                    1e-12_SP, irank, p_Fwork, nwork, info)
      end if

      ! Check if solution procedure terminated without errors
      if (info .eq. 0) then
        
        ! Copy solution vector to output array
        call lalg_copyVector(p_Fb,Fx,ndim4)
        
        ! Release internal memory
        deallocate(p_Fa,p_Fb,p_Fwork,p_Ipvt)
        
        ! That is it
        return
      else
        call output_line('Unable to solve least-squares minimisation problem!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mprim_leastSquareMinSP')
        call sys_halt()       
      end if
      
    end if

  end subroutine mprim_leastSquaresMinSP

  !*****************************************************************************

!<subroutine>

  subroutine mprim_transposeMatrix1DP(Da)

!<description>
    ! This subroutine computes the transpose of the source matrix A in-place.
!</description>

!<inputoutput>
    ! On input: source matrix that should be transposed.
    ! On output: transposed matrix.
    real(DP), dimension(:,:), intent(inout) :: Da
!</inputoutput>
!</subroutine>

#if !defined(USE_INTEL_MKL)
    real(DP), dimension(:), allocatable :: Daux
    integer :: i,j,n
#endif

    ! Check if matrix is square matrix
    if (size(Da,1) .ne. size(Da,2)) then
      call output_line('Matrix must be square matrix for in-place transposition!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mprim_transposeMatrix1DP')
      call sys_halt()
    end if
    
#ifdef USE_INTEL_MKL
    ! Call in-place transposition from Intel MKL library
    call mkl_dimatcopy('C', 'T', size(Da,1), size(Da,2), 1.0_DP, Da, size(Da,1), size(Da,2))
#else
    n = size(Da,1)
    allocate(Daux(n))

    ! Loop over all columns of the square matrix
    do i=1,n
      ! Make a copy of the i-th column
      Daux(i:n) = Da(i:n,i)
      
      ! Swap content of i-th row and column
      do j=i,n
        Da(j,i) = Da(i,j)
        Da(i,j) = Daux(j)
      end do
    end do

    deallocate(Daux)
#endif

  end subroutine mprim_transposeMatrix1DP

  !*****************************************************************************

!<subroutine>

  subroutine mprim_transposeMatrix1SP(Fa)

!<description>
    ! This subroutine computes the transpose of the source matrix A in-place.
!</description>

!<inputoutput>
    ! On input: source matrix that should be transposed.
    ! On output: transposed matrix.
    real(SP), dimension(:,:), intent(inout) :: Fa
!</inputoutput>
!</subroutine>

#if !defined(USE_INTEL_MKL)
    real(SP), dimension(:), allocatable :: Faux
    integer :: i,j,n
#endif

    ! Check if matrix is square matrix
    if (size(Fa,1) .ne. size(Fa,2)) then
      call output_line('Matrix must be square matrix for in-place transposition!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mprim_transposeMatrix1SP')
      call sys_halt()
    end if
    
#ifdef USE_INTEL_MKL
    ! Call in-place transposition from Intel MKL library
    call mkl_simatcopy('C', 'T', size(Fa,1), size(Fa,2), 1.0_SP, Fa, size(Fa,1), size(Fa,2))
#else
    n = size(Fa,1)
    allocate(Faux(n))

    ! Loop over all columns of the square matrix
    do i=1,n
      ! Make a copy of the i-th column
      Faux(i:n) = Fa(i:n,i)
      
      ! Swap content of i-th row and column
      do j=i,n
        Fa(j,i) = Fa(i,j)
        Fa(i,j) = Faux(j)
      end do
    end do

    deallocate(Faux)
#endif

  end subroutine mprim_transposeMatrix1SP

  !*****************************************************************************

!<subroutine>

  subroutine mprim_transposeMatrix2DP(DaSrc,DaDest)

!<description>
    ! This subroutine computes the transpose of the source matrix A out-of-place.
!</description>

!<input>
    ! Source matrix that should be transposed.
    real(DP), dimension(:,:), intent(in) :: DaSrc
!</input>

!<output>
    ! Destination matrix containing the transposed of matrix A on output
    real(DP), dimension(:,:), intent(out) :: DaDest
!</output>
!</subroutine>

#if !defined(USE_INTEL_MKL)
    integer :: i,j
#endif

    ! Check if matrix dimensions are compatible
    if ((size(DaSrc,1) .ne. size(DaDest,2)) .or.&
        (size(DaSrc,2) .ne. size(DaDest,1))) then
      call output_line('Matrix dimensions mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mprim_transposeMatrix2DP')
      call sys_halt()
    end if

#ifdef USE_INTEL_MKL
    ! Call out-of-place transposition from Intel MKL library
    call mkl_domatcopy('C', 'T', size(DaSrc,1), size(DaSrc,2), 1.0_DP,&
                       DaSrc, size(DaSrc,1), DaDest, size(DaDest,1))
#else
    ! Loop over all columns of the square matrix
    !$omp parallel do private(j) default(shared)
    do i=1,size(DaSrc,1)
      ! Swap content of i-th row and column
      do j=1,size(DaSrc,2)
        DaDest(j,i) = DaSrc(i,j)
      end do
    end do
    !$omp end parallel do
#endif

  end subroutine mprim_transposeMatrix2DP

  !*****************************************************************************

!<subroutine>

  subroutine mprim_transposeMatrix2SP(FaSrc,FaDest)

!<description>
    ! This subroutine computes the transpose of the source matrix A out-of-place.
!</description>

!<input>
    ! Source matrix that should be transposed.
    real(SP), dimension(:,:), intent(in) :: FaSrc
!</input>

!<output>
    ! Destination matrix containing the transposed of matrix A on output
    real(SP), dimension(:,:), intent(out) :: FaDest
!</output>
!</subroutine>

#if !defined(USE_INTEL_MKL)
    integer :: i,j
#endif

    ! Check if matrix dimensions are compatible
    if ((size(FaSrc,1) .ne. size(FaDest,2)) .or.&
        (size(FaSrc,2) .ne. size(FaDest,1))) then
      call output_line('Matrix dimensions mismatch!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mprim_transposeMatrix2SP')
      call sys_halt()
    end if

#ifdef USE_INTEL_MKL
    ! Call out-of-place transposition from Intel MKL library
    call mkl_somatcopy('C', 'T', size(FaSrc,1), size(FaSrc,2), 1.0_SP,&
                       FaSrc, size(FaSrc,1), FaDest, size(FaDest,1))
#else
    ! Loop over all columns of the square matrix
    !$omp parallel do private(j) default(shared)
    do i=1,size(FaSrc,1)
      ! Swap content of i-th row and column
      do j=1,size(FaSrc,2)
        FaDest(j,i) = FaSrc(i,j)
      end do
    end do
    !$omp end parallel do
#endif

  end subroutine mprim_transposeMatrix2SP

  !*****************************************************************************

!<function>

  pure real(DP) function mprim_hornerDP(dx,Dai)
  
!<description>
  ! Applies the horner scheme to evaluate a polynomial
  !   p(x) = a_0 + a_1 x + a_2 x^2 + ...
  ! given as a list of coefficients.
!</description>

!<input>
  ! Point where to evaluate
  real(DP), intent(in) :: dx

  ! List of coefficients
  real(DP), dimension(:), intent(in) :: Dai
!</input>

!<result>
  ! p(x)
!</result>

!</function>
    
    integer :: i
    
    mprim_hornerDP = 0.0_DP
    do i = ubound(Dai,1),1,-1
      mprim_hornerDP = mprim_hornerDP * dx + Dai(i)
    end do
    
  end function

  !*****************************************************************************

!<function>

  pure real(DP) function mprim_hornerd1DP(dx,Dai)
  
!<description>
  ! Applies the extended horner scheme to evaluate the derivative p'(x) of a polynomial
  !   p(x) = a_0 + a_1 x + a_2 x^2 + ...
  ! given as a list of coefficients.
!</description>

!<input>
  ! Point where to evaluate
  real(DP), intent(in) :: dx

  ! List of coefficients
  real(DP), dimension(:), intent(in) :: Dai
!</input>

!<result>
  ! p'(x)
!</result>

!</function>

    integer :: i
    real(DP) :: dy
    
    dy = 0.0_DP
    mprim_hornerd1DP = 0.0_DP
    do i = ubound(Dai,1),1,-1
      mprim_hornerd1DP = mprim_hornerd1DP * dx + dy
      dy = dy * dx + Dai(i)
    end do
    
  end function
  
  !*****************************************************************************

!<function>

  pure real(DP) function mprim_hornerd2DP(dx,Dai)
  
!<description>
  ! Applies the extended horner scheme to evaluate the 2nd derivative p''(x) of a polynomial
  !   p(x) = a_0 + a_1 x + a_2 x^2 + ...
  ! given as a list of coefficients.
!</description>

!<input>
  ! Point where to evaluate
  real(DP), intent(in) :: dx

  ! List of coefficients
  real(DP), dimension(:), intent(in) :: Dai
!</input>

!<result>
  ! p''(x)
!</result>

!</function>

    integer :: i
    real(DP) :: dy, dy1
    
    dy = 0.0_DP
    dy1 = 0.0_DP
    mprim_hornerd2DP = 0.0_DP
    do i = ubound(Dai,1),1,-1
      mprim_hornerd2DP = mprim_hornerd2DP * dx + dy1
      dy1 = dy1 * dx + dy
      dy = dy * dx + Dai(i)
    end do
    
  end function
  
  !*****************************************************************************

!<function>

  pure real(DP) function mprim_hornerSP(fx,Fai)
  
!<description>
  ! Applies the horner scheme to evaluate a polynomial
  !   p(x) = a_0 + a_1 x + a_2 x^2 + ...
  ! given as a list of coefficients.
!</description>

!<input>
  ! Point where to evaluate
  real(SP), intent(in) :: fx

  ! List of coefficients
  real(SP), dimension(:), intent(in) :: Fai
!</input>

!<result>
  ! p(x)
!</result>

!</function>
    
    integer :: i
    
    mprim_hornerSP = 0.0_DP
    do i = ubound(Fai,1),1,-1
      mprim_hornerSP = mprim_hornerSP * fx + Fai(i)
    end do
    
  end function

  !*****************************************************************************

!<function>

  pure real(DP) function mprim_hornerd1SP(fx,Fai)
  
!<description>
  ! Applies the extended horner scheme to evaluate the derivative p'(x) of a polynomial
  !   p(x) = a_0 + a_1 x + a_2 x^2 + ...
  ! given as a list of coefficients.
!</description>

!<input>
  ! Point where to evaluate
  real(SP), intent(in) :: fx

  ! List of coefficients
  real(SP), dimension(:), intent(in) :: Fai
!</input>

!<result>
  ! p'(x)
!</result>

!</function>

    integer :: i
    real(SP) :: fy
    
    fy = 0.0_DP
    mprim_hornerd1SP = 0.0_DP
    do i = ubound(Fai,1),1,-1
      mprim_hornerd1SP = mprim_hornerd1SP * fx + fy
      fy = fy * fx + Fai(i)
    end do
    
  end function
  
  !*****************************************************************************

!<function>

  pure real(DP) function mprim_hornerd2SP(fx,Fai)
  
!<description>
  ! Applies the extended horner scheme to evaluate the 2nd derivative p''(x) of a polynomial
  !   p(x) = a_0 + a_1 x + a_2 x^2 + ...
  ! given as a list of coefficients.
!</description>

!<input>
  ! Point where to evaluate
  real(SP), intent(in) :: fx

  ! List of coefficients
  real(SP), dimension(:), intent(in) :: Fai
!</input>

!<result>
  ! p''(x)
!</result>

!</function>

    integer :: i
    real(SP) :: fy, fy1
    
    fy = 0.0_DP
    fy1 = 0.0_DP
    mprim_hornerd2SP = 0.0_DP
    do i = ubound(Fai,1),1,-1
      mprim_hornerd2SP = mprim_hornerd2SP * fx + fy1
      fy1 = fy1 * fx + fy
      fy = fy * fx + Fai(i)
    end do
    
  end function
  
  ! ***************************************************************************

!<subroutine>

  elemental subroutine mprim_polarToCartesian2D(dr,dphi,dx,dy)

!<description>
  ! Calculates 2D cartesian coordinates from polar coordinates
!</description>

!<input>
  ! Radius and angle of a point.
  real(DP), intent(in) :: dr, dphi
!</input>

!<inputoutput>
  ! X/Y cartesian coordinate of the point
  real(DP), intent(out) :: dx,dy
!</inputoutput>

!</subroutine>

    dx = dr * cos(dphi)
    dy = dr * sin(dphi)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  elemental subroutine mprim_polarToCartesian3D(dr,dphi,dpsi,dx,dy,dz)

!<description>
  ! Calculates 3D cartesian coordinates from polar coordinates
!</description>

!<input>
  ! Radius, Y-angle and Z-angle of a point.
  real(DP), intent(in) :: dr, dphi, dpsi
!</input>

!<inputoutput>
  ! X/Y/Z cartesian coordinate of the point.
  real(DP), intent(out) :: dx,dy,dz
!</inputoutput>

!</subroutine>

    dx = dr * cos(dphi)
    dy = dr * sin(dphi)
    dz = dr * sin(dpsi)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  elemental subroutine mprim_cartesianToPolar2D(dx,dy,dr,dphi)

!<description>
  ! Calculates 2D cartesian coordinates from polar coordinates
!</description>

!<input>
  ! X/Y cartesian coordinate of a point
  real(DP), intent(in) :: dx,dy
!</input>

!<inputoutput>
  ! Radius and angle of a point.
  real(DP), intent(out) :: dr, dphi
!</inputoutput>

!</subroutine>

    dr = sqrt(dx**2 + dy**2)
    dphi = atan2 (dy,dx)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  elemental subroutine mprim_cartesianToPolar3D(dx,dy,dz,dr,dphi,dpsi)

!<description>
  ! Calculates 3D cartesian coordinates from polar coordinates
!</description>

!<input>
  ! X/Y/Z cartesian coordinate of a point
  real(DP), intent(in) :: dx,dy,dz
!</input>

!<inputoutput>
  ! Radius, Y-angle and Z-angle of a point.
  real(DP), intent(out) :: dr, dphi, dpsi
!</inputoutput>

!</subroutine>

    dr = sqrt(dx**2 + dy**2 + dz**2)
    dphi = atan2 (dy,dx)
    
    if (dr .gt. SYS_EPSREAL_DP) then
      dpsi = asin (dy/dr)
    else
      dpsi = 0.0_DP
    end if

  end subroutine

end module mprimitives

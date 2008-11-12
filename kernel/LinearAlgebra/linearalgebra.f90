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
!# 8.) lalg_errorNormXXX
!#     -> Calculate a specific norm from the difference of two vectors
!#
!# 9.) lalg_vectorSortXXX
!#     -> Resort the entries of a vector according to a given permutation
!#
!# 10.) lalg_vectorAddScalarXXX
!#      -> Adds a scalar to each entry of a vector
!#
!# 11.) lalg_vectorCompMultXXX
!#      -> Multiply two vectors componentwise
!# </purpose>
!##############################################################################

module linearalgebra

  use fsystem

  implicit none
  
  ! Remark:
  ! As currently quad precision is equal to double precision, the quad
  ! precision routines are not part of the interfaces to avoid that the
  ! compiler complains about ambiguous interfaces...
  
  interface lalg_copyVectorInt
    module procedure lalg_copyVectorI8
    module procedure lalg_copyVectorI8I16
    module procedure lalg_copyVectorI8I32
    module procedure lalg_copyVectorI8I64
    module procedure lalg_copyVectorI16
    module procedure lalg_copyVectorI16I8
    module procedure lalg_copyVectorI16I32
    module procedure lalg_copyVectorI16I64
    module procedure lalg_copyVectorI32
    module procedure lalg_copyVectorI32I8
    module procedure lalg_copyVectorI32I16
    module procedure lalg_copyVectorI32I64
    module procedure lalg_copyVectorI64
    module procedure lalg_copyVectorI64I8
    module procedure lalg_copyVectorI64I16
    module procedure lalg_copyVectorI64I32
  end interface
  
  interface lalg_copyVectorInt2D
    module procedure lalg_copyVectorI8_2D
    module procedure lalg_copyVectorI8I16_2D
    module procedure lalg_copyVectorI8I32_2D
    module procedure lalg_copyVectorI8I64_2D
    module procedure lalg_copyVectorI16_2D
    module procedure lalg_copyVectorI16I8_2D
    module procedure lalg_copyVectorI16I32_2D
    module procedure lalg_copyVectorI16I64_2D
    module procedure lalg_copyVectorI32_2D
    module procedure lalg_copyVectorI32I8_2D
    module procedure lalg_copyVectorI32I16_2D
    module procedure lalg_copyVectorI32I64_2D
    module procedure lalg_copyVectorI64_2D
    module procedure lalg_copyVectorI64I8_2D
    module procedure lalg_copyVectorI64I16_2D
    module procedure lalg_copyVectorI64I32_2D
  end interface  

  interface lalg_copyVectorReal
    module procedure lalg_copyVectorSngl
    module procedure lalg_copyVectorDble
    !module procedure lalg_copyVectorQuad
    module procedure lalg_copyVectorSnglDbl
    !module procedure lalg_copyVectorSnglQuad
    module procedure lalg_copyVectorDblSngl
    !module procedure lalg_copyVectorDblQuad
    !module procedure lalg_copyVectorQuadSngl
    !module procedure lalg_copyVectorQuadDbl
  end interface

  interface lalg_copyVectorReal2D
    module procedure lalg_copyVectorSngl2D
    module procedure lalg_copyVectorDble2D
    !module procedure lalg_copyVectorQuad2D
    module procedure lalg_copyVectorSnglDbl2D
    !module procedure lalg_copyVectorSnglQuad2D
    module procedure lalg_copyVectorDblSngl2D
    !module procedure lalg_copyVectorDblQuad2D
    !module procedure lalg_copyVectorQuadSngl2D
    !module procedure lalg_copyVectorQuadDbl2D
  end interface

  interface lalg_copyVector
    module procedure lalg_copyVectorSngl
    module procedure lalg_copyVectorDble
    !module procedure lalg_copyVectorQuad
    module procedure lalg_copyVectorSnglDbl
    !module procedure lalg_copyVectorSnglQuad
    module procedure lalg_copyVectorDblSngl
    !module procedure lalg_copyVectorDblQuad
    !module procedure lalg_copyVectorQuadSngl
    !module procedure lalg_copyVectorQuadDbl
    module procedure lalg_copyVectorI8
    module procedure lalg_copyVectorI8I16
    module procedure lalg_copyVectorI8I32
    module procedure lalg_copyVectorI8I64
    module procedure lalg_copyVectorI16
    module procedure lalg_copyVectorI16I8
    module procedure lalg_copyVectorI16I32
    module procedure lalg_copyVectorI16I64
    module procedure lalg_copyVectorI32
    module procedure lalg_copyVectorI32I8
    module procedure lalg_copyVectorI32I16
    module procedure lalg_copyVectorI32I64
    module procedure lalg_copyVectorI64
    module procedure lalg_copyVectorI64I8
    module procedure lalg_copyVectorI64I16
    module procedure lalg_copyVectorI64I32
    module procedure lalg_copyVectorLogical
    module procedure lalg_copyVectorChar
    module procedure lalg_copyVectorSngl2D
    module procedure lalg_copyVectorDble2D
    !module procedure lalg_copyVectorQuad2D
    module procedure lalg_copyVectorSnglDbl2D
    !module procedure lalg_copyVectorSnglQuad2D
    module procedure lalg_copyVectorDblSngl2D
    !module procedure lalg_copyVectorDblQuad2D
    !module procedure lalg_copyVectorQuadSngl2D
    !module procedure lalg_copyVectorQuadDbl2D
    module procedure lalg_copyVectorI8_2D
    module procedure lalg_copyVectorI8I16_2D
    module procedure lalg_copyVectorI8I32_2D
    module procedure lalg_copyVectorI8I64_2D
    module procedure lalg_copyVectorI16_2D
    module procedure lalg_copyVectorI16I8_2D
    module procedure lalg_copyVectorI16I32_2D
    module procedure lalg_copyVectorI16I64_2D
    module procedure lalg_copyVectorI32_2D
    module procedure lalg_copyVectorI32I8_2D
    module procedure lalg_copyVectorI32I16_2D
    module procedure lalg_copyVectorI32I64_2D
    module procedure lalg_copyVectorI64_2D
    module procedure lalg_copyVectorI64I8_2D
    module procedure lalg_copyVectorI64I16_2D
    module procedure lalg_copyVectorI64I32_2D
    module procedure lalg_copyVectorLogical2D
    module procedure lalg_copyVectorChar2D
  end interface

  interface lalg_scaleVector
    module procedure lalg_scaleVectorSngl
    module procedure lalg_scaleVectorDble
    !module procedure lalg_scaleVectorQuad
    module procedure lalg_scaleVectorSngl2D
    module procedure lalg_scaleVectorDble2D
    !module procedure lalg_scaleVectorQuad2D
  end interface

  interface lalg_clearVectorInt
    module procedure lalg_clearVectorI8
    module procedure lalg_clearVectorI16
    module procedure lalg_clearVectorI32
    module procedure lalg_clearVectorI64
  end interface

  interface lalg_clearVectorInt2D
    module procedure lalg_clearVectorI8_2D
    module procedure lalg_clearVectorI16_2D
    module procedure lalg_clearVectorI32_2D
    module procedure lalg_clearVectorI64_2D
  end interface

  interface lalg_clearVector
    module procedure lalg_clearVectorSngl
    module procedure lalg_clearVectorDble
    !module procedure lalg_clearVectorQuad
    module procedure lalg_clearVectorI8
    module procedure lalg_clearVectorI16
    module procedure lalg_clearVectorI32
    module procedure lalg_clearVectorI64
    module procedure lalg_clearVectorSngl2D
    module procedure lalg_clearVectorDble2D
    !module procedure lalg_clearVectorQuad2D
    module procedure lalg_clearVectorI8_2D
    module procedure lalg_clearVectorI16_2D
    module procedure lalg_clearVectorI32_2D
    module procedure lalg_clearVectorI64_2D
  end interface

  interface lalg_setVectorInt
    module procedure lalg_setVectorI8
    module procedure lalg_setVectorI16
    module procedure lalg_setVectorI32
    module procedure lalg_setVectorI64
  end interface

  interface lalg_setVectorInt2D
    module procedure lalg_setVectorI8_2D
    module procedure lalg_setVectorI16_2D
    module procedure lalg_setVectorI32_2D
    module procedure lalg_setVectorI64_2D
  end interface

  interface lalg_setVector
    module procedure lalg_setVectorSngl
    module procedure lalg_setVectorDble
    !module procedure lalg_setVectorQuad
    module procedure lalg_setVectorI8
    module procedure lalg_setVectorI16
    module procedure lalg_setVectorI32
    module procedure lalg_setVectorI64
    module procedure lalg_setVectorLogical
    module procedure lalg_setVectorChar
    module procedure lalg_setVectorSngl2D
    module procedure lalg_setVectorDble2D
    !module procedure lalg_setVectorQuad2D
    module procedure lalg_setVectorI8_2D
    module procedure lalg_setVectorI16_2D
    module procedure lalg_setVectorI32_2D
    module procedure lalg_setVectorI64_2D
    module procedure lalg_setVectorLogical2D
    module procedure lalg_setVectorChar2D
  end interface

  interface lalg_vectorLinearComb
    module procedure lalg_vectorLinearCombSngl
    module procedure lalg_vectorLinearCombDble
    !module procedure lalg_vectorLinearCombQuad
    module procedure lalg_vectorLinearCombSnglDble
    !module procedure lalg_vectorLinearCombSnglQuad
    !module procedure lalg_vectorLinearCombDblQuad
    module procedure lalg_vectorLinearCombSngl2D
    module procedure lalg_vectorLinearCombDble2D
    !module procedure lalg_vectorLinearCombQuad2D
    module procedure lalg_vectorLinearCombSnglDble2D
    !module procedure lalg_vectorLinearCombSnglQuad2D
    !module procedure lalg_vectorLinearCombDblQuad2D
  end interface

  interface lalg_scalarProduct
    module procedure lalg_scalarProductSngl
    module procedure lalg_scalarProductDble
    !module procedure lalg_scalarProductQuad
    module procedure lalg_scalarProductSngl2D
    module procedure lalg_scalarProductDble2D
    !module procedure lalg_scalarProductQuad2D
  end interface

  interface lalg_norm
    module procedure lalg_normSngl
    module procedure lalg_normDble
    !module procedure lalg_normQuad
  end interface

  interface lalg_errorNorm
    module procedure lalg_errorNormSngl
    module procedure lalg_errorNormDble
    !module procedure lalg_errorNormQuad
  end interface

  interface lalg_vectorSortInt    
    module procedure lalg_vectorSortI32
    module procedure lalg_vectorSortI64
  end interface

  interface lalg_vectorSort
    module procedure lalg_vectorSortSngl
    module procedure lalg_vectorSortDble
    !module procedure lalg_vectorSortQuad
    module procedure lalg_vectorSortI32
    module procedure lalg_vectorSortI64
  end interface

  interface lalg_vectorAddScalarInt
    module procedure lalg_vectorAddScalarI8
    module procedure lalg_vectorAddScalarI16
    module procedure lalg_vectorAddScalarI32
    module procedure lalg_vectorAddScalarI64
  end interface

  interface lalg_vectorAddScalar
    module procedure lalg_vectorAddScalarSngl
    module procedure lalg_vectorAddScalarDble
    !module procedure lalg_vectorAddScalarQuad
    module procedure lalg_vectorAddScalarSngl2D
    module procedure lalg_vectorAddScalarDble2D
    !module procedure lalg_vectorAddScalarQuad2D
    module procedure lalg_vectorAddScalarI8
    module procedure lalg_vectorAddScalarI16
    module procedure lalg_vectorAddScalarI32
    module procedure lalg_vectorAddScalarI64
  end interface

!<constants>

!<constantblock description="Constants identifying vector norms">

  ! Sum of the absolute values of entries
  integer, parameter :: LINALG_NORMSUM    = -1

  ! Euclidian vector norm: (vector,vector)
  integer, parameter :: LINALG_NORMEUCLID = 0

  ! $l_1$-norm: 1/NEQ * sum(abs(entries))
  integer, parameter :: LINALG_NORML1     = 1
  
  ! $l_2$-norm: 1/sqrt(NEQ) * (vector,vector)
  integer, parameter :: LINALG_NORML2     = 2
  
  ! max-norm
  integer, parameter :: LINALG_NORMMAX    = 3
  
!</constantblock>

!</constants>
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSngl (Fx,Fy,n)
  
!<description>
  ! Copies a single precision vector: Fy = Fx
!</description>

!<input>
  
  ! Source vector
  real(SP), dimension(:), intent(IN) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:), intent(OUT) :: Fy
  
!</output>
  
!</subroutine>

    if (.not. present(n)) then
      call SCOPY(size(Fx),Fx,1,Fy,1)
    else
      call SCOPY(n,Fx,1,Fy,1)
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDble (Dx,Dy,n)
  
!<description>
  ! Copies a double precision vector dx: Dy = Dx
!</description>

!<input>
  
  ! Source vector
  real(DP), dimension(:), intent(IN) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:), intent(OUT) :: Dy
  
!</output>
  
!</subroutine>

    if (.not. present(n)) then
      call DCOPY(size(Dx),Dx,1,Dy,1)
    else
      call DCOPY(n,Dx,1,Dy,1)
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQuad (Qx,Qy,n)
  
!<description>
  ! Copies a quad precision vector dx: Qy = Qx
!</description>

!<input>
  
  ! Source vector
  real(QP), dimension(:), intent(IN) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:), intent(OUT) :: Qy
  
!</output>
  
!</subroutine>

  integer :: i

    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Qx)
        Qy(i) = Qx(i)
      end do
      !$omp end parallel do
      
      !call QCOPY(size(Qx),Qx,1,Qy,1)
      
    else

      !$omp parallel do
      do i = 1, n
        Qy(i) = Qx(i)
      end do
      !$omp end parallel do
      
      !call QCOPY(n,Qx,1,Qy,1)

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSnglDbl (Fx,Dy,n)
  
!<description>
  ! Copies single precision vector to double precision vector: Dy = Fx
!</description>

!<input>
  
  ! Source vector
  real(SP), dimension(:), intent(IN) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:), intent(OUT) :: Dy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i=1,size(Fx)
        Dy(i) = real(Fx(i),DP)
      end do
      !$omp end parallel do
    
    else

      !$omp parallel do
      do i=1,n
        Dy(i) = real(Fx(i),DP)
      end do
      !$omp end parallel do
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSnglQuad (Fx,Qy,n)
  
!<description>
  ! Copies single precision vector to quad precision vector: Qy = Fx
!</description>

!<input>
  
  ! Source vector
  real(SP), dimension(:), intent(IN) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:), intent(OUT) :: Qy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i=1,size(Fx)
        Qy(i) = real(Fx(i),QP)
      end do
      !$omp end parallel do
    
    else

      !$omp parallel do
      do i=1,n
        Qy(i) = real(Fx(i),QP)
      end do
      !$omp end parallel do
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDblSngl (Dx,Fy,n)
  
!<description>
  ! Copies double precision vector to single precision vector: Fy = Dx
!</description>

!<input>
  
  ! Source vector
  real(DP), dimension(:), intent(IN) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:), intent(OUT) :: Fy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then

      !$omp parallel do
      do i=1,size(Dx)
        Fy(i) = real(Dx(i),SP)
      end do
      !$omp end parallel do

    else
    
      !$omp parallel do
      do i=1,n
        Fy(i) = real(Dx(i),SP)
      end do
      !$omp end parallel do

    end if  

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDblQuad (Dx,Qy,n)
  
!<description>
  ! Copies double precision vector to quad precision vector: Qy = Dx
!</description>

!<input>
  
  ! Source vector
  real(DP), dimension(:), intent(IN) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:), intent(OUT) :: Qy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then

      !$omp parallel do
      do i=1,size(Dx)
        Qy(i) = real(Dx(i),QP)
      end do
      !$omp end parallel do

    else
    
      !$omp parallel do
      do i=1,n
        Qy(i) = real(Dx(i),QP)
      end do
      !$omp end parallel do

    end if  

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQuadSngl (Qx,Fy,n)
  
!<description>
  ! Copies quad precision vector to single precision vector: Fy = Qx
!</description>

!<input>
  
  ! Source vector
  real(QP), dimension(:), intent(IN) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:), intent(OUT) :: Fy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then

      !$omp parallel do
      do i=1,size(Qx)
        Fy(i) = real(Qx(i),SP)
      end do
      !$omp end parallel do

    else
    
      !$omp parallel do
      do i=1,n
        Fy(i) = real(Qx(i),SP)
      end do
      !$omp end parallel do

    end if  

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQuadDbl (Qx,Dy,n)
  
!<description>
  ! Copies quad precision vector to double precision vector: Dy = Qx
!</description>

!<input>
  
  ! Source vector
  real(QP), dimension(:), intent(IN) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:), intent(OUT) :: Dy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then

      !$omp parallel do
      do i=1,size(Qx)
        Dy(i) = real(Qx(i),DP)
      end do
      !$omp end parallel do

    else
    
      !$omp parallel do
      do i=1,n
        Dy(i) = real(Qx(i),DP)
      end do
      !$omp end parallel do

    end if  

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI8 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I8), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = Ix(i)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = Ix(i)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI8I16 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I8), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I16)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = int(Ix(i),I16)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI8I32 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I8), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I32)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = int(Ix(i),I32)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI8I64 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I8), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I64)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = int(Ix(i),I64)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI16 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I16), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = Ix(i)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = Ix(i)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI16I8 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I16), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I8)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = int(Ix(i),I8)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI16I32 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I16), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I32)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = int(Ix(i),I32)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI16I64 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I16), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I64)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = int(Ix(i),I64)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI32 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = Ix(i)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = Ix(i)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI32I8 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I8)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = int(Ix(i),I8)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI32I16 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I16)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = int(Ix(i),I16)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI32I64 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I64)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = int(Ix(i),I64)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI64 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I64), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = Ix(i)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = Ix(i)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI64I8 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I64), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I8)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = int(Ix(i),I8)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI64I16 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I64), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I16)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = int(Ix(i),I16)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI64I32 (Ix,Iy,n)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I64), dimension(:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I32)
      end do
      !$omp end parallel do
      
    else

      !$omp parallel do
      do i = 1, n
        Iy(i) = int(Ix(i),I32)
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorLogical (Lx,Ly,n)
  
!<description>
  ! Copies a logical vector Lx: Ly = Lx
!</description>

!<input>
  
  ! Source vector
  logical, dimension(:), intent(IN) :: Lx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  logical, dimension(:), intent(OUT) :: Ly
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i=1,size(Lx)
        Ly(i) = Lx(i)
      end do
      !$omp end parallel do
    
    else

      !$omp parallel do
      do i=1,n
        Ly(i) = Lx(i)
      end do
      !$omp end parallel do
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorChar (Sx,Sy,n)
  
!<description>
  ! Copies a character vector sx: Sy = Sx
!</description>

!<input>
  
  ! Source vector
  character, dimension(:), intent(IN) :: Sx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  
  ! Destination vector
  character, dimension(:), intent(OUT) :: Sy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i=1,size(Sx)
        Sy(i) = Sx(i)
      end do
      !$omp end parallel do
    
    else

      !$omp parallel do
      do i=1,n
        Sy(i) = Sx(i)
      end do
      !$omp end parallel do
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSngl2D (Fx,Fy,n,m)
  
!<description>
  ! Copies a single precision vector: Fy = Fx
!</description>

!<input>
  
  ! Source vector
  real(SP), dimension(:,:), intent(IN) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:,:), intent(OUT) :: Fy
  
!</output>
  
!</subroutine>

    if (present(n) .and. present(m)) then
      call SCOPY(n*m,Fx,1,Fy,1)
    else
      call SCOPY(size(Fx,1)*size(Fx,2),Fx,1,Fy,1)
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDble2D (Dx,Dy,n,m)
  
!<description>
  ! Copies a double precision vector: Dy = Dx
!</description>

!<input>
  
  ! Source vector
  real(DP), dimension(:,:), intent(IN) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m
  
!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:,:), intent(OUT) :: Dy
  
!</output>
  
!</subroutine>

    if (present(n) .and. present(m)) then
      call DCOPY(n*m,Dx,1,Dy,1)
    else
      call DCOPY(size(Dx,1)*size(Dx,2),Dx,1,Dy,1)
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQuad2D (Qx,Qy,n,m)
  
!<description>
  ! Copies a quad precision vector: Qy = Qx
!</description>

!<input>
  
  ! Source vector
  real(QP), dimension(:,:), intent(IN) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m
  
!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:,:), intent(OUT) :: Qy
  
!</output>
  
!</subroutine>

  integer :: i,j

    if (present(n) .and. present(m)) then
    
      !$omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Qy(i,j) = Qx(i,j)
        end do
      end do
      !$omp end parallel do
      
      !call QCOPY(n*m,Qx,1,Qy,1)
      
    else

      !$omp parallel do private(i)
      do j = 1, size(Qx,2)
        do i = 1, size(Qx,1)
          Qy(i,j) = Qx(i,j)
        end do
      end do
      !$omp end parallel do

      !call QCOPY(size(Qx,1)*size(Qx,2),Qx,1,Qy,1)
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSnglDbl2D (Fx,Dy,n,m)
  
!<description>
  ! Copies single precision vector to double precision vector: Dy = Fx
!</description>

!<input>
  
  ! Source vector
  real(SP), dimension(:,:), intent(IN) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:,:), intent(OUT) :: Dy
  
!</output>
  
!</subroutine>
  integer :: i,j
  
    if (present(n) .and. present(m)) then

      !$omp parallel do private(i)
      do j=1,m
        do i=1,n
          Dy(i,j) = real(Fx(i,j),DP)
        end do
      end do
     !$omp end parallel do

    else

      !$omp parallel do private(i)
      do j=1,size(Fx,2)
        do i=1,size(Fx,1)
          Dy(i,j) = real(Fx(i,j),DP)
        end do
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSnglQuad2D (Fx,Qy,n,m)
  
!<description>
  ! Copies single precision vector to quad precision vector: Qy = Fx
!</description>

!<input>
  
  ! Source vector
  real(SP), dimension(:,:), intent(IN) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:,:), intent(OUT) :: Qy
  
!</output>
  
!</subroutine>
  integer :: i,j
  
    if (present(n) .and. present(m)) then

      !$omp parallel do private(i)
      do j=1,m
        do i=1,n
          Qy(i,j) = real(Fx(i,j),QP)
        end do
      end do
     !$omp end parallel do

    else

      !$omp parallel do private(i)
      do j=1,size(Fx,2)
        do i=1,size(Fx,1)
          Qy(i,j) = real(Fx(i,j),QP)
        end do
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDblSngl2D (Dx,Fy,n,m)
  
!<description>
  ! Copies double precision vector to single precision vector: Fy = Dx
!</description>

!<input>
  
  ! Source vector
  real(DP), dimension(:,:), intent(IN) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:,:), intent(OUT) :: Fy
  
!</output>
  
!</subroutine>
  integer :: i,j
  
    if (present(n) .and. present(m)) then

      !$omp parallel do private(i)
      do j=1,m
        do i=1,n
          Fy(i,j) = real(Dx(i,j),SP)
        end do
      end do
      !$omp end parallel do

    else

      !$omp parallel do private(i)
      do j=1,size(Dx,2)
        do i=1,size(Dx,1)
          Fy(i,j) = real(Dx(i,j),SP)
        end do
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDblQuad2D (Dx,Qy,n,m)
  
!<description>
  ! Copies double precision vector to quad precision vector: Qy = Dx
!</description>

!<input>
  
  ! Source vector
  real(DP), dimension(:,:), intent(IN) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:,:), intent(OUT) :: Qy
  
!</output>
  
!</subroutine>
  integer :: i,j
  
    if (present(n) .and. present(m)) then

      !$omp parallel do private(i)
      do j=1,m
        do i=1,n
          Qy(i,j) = real(Dx(i,j),QP)
        end do
      end do
      !$omp end parallel do

    else

      !$omp parallel do private(i)
      do j=1,size(Dx,2)
        do i=1,size(Dx,1)
          Qy(i,j) = real(Dx(i,j),QP)
        end do
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQuadSngl2D (Qx,Fy,n,m)
  
!<description>
  ! Copies quad precision vector to single precision vector: Fy = Qx
!</description>

!<input>
  
  ! Source vector
  real(QP), dimension(:,:), intent(IN) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:,:), intent(OUT) :: Fy
  
!</output>
  
!</subroutine>
  integer :: i,j
  
    if (present(n) .and. present(m)) then

      !$omp parallel do private(i)
      do j=1,m
        do i=1,n
          Fy(i,j) = real(Qx(i,j),SP)
        end do
      end do
      !$omp end parallel do

    else

      !$omp parallel do private(i)
      do j=1,size(Qx,2)
        do i=1,size(Qx,1)
          Fy(i,j) = real(Qx(i,j),SP)
        end do
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQuadDbl2D (Qx,Dy,n,m)
  
!<description>
  ! Copies quad precision vector to double precision vector: Dy = Qx
!</description>

!<input>
  
  ! Source vector
  real(QP), dimension(:,:), intent(IN) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:,:), intent(OUT) :: Dy
  
!</output>
  
!</subroutine>
  integer :: i,j
  
    if (present(n) .and. present(m)) then

      !$omp parallel do private(i)
      do j=1,m
        do i=1,n
          Dy(i,j) = real(Qx(i,j),DP)
        end do
      end do
     !$omp end parallel do

    else

      !$omp parallel do private(i)
      do j=1,size(Qx,2)
        do i=1,size(Qx,1)
          Dy(i,j) = real(Qx(i,j),DP)
        end do
      end do
      !$omp end parallel do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI8_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I8), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = Ix(i,j)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = Ix(i,j)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI8I16_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I8), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I16)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I16)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI8I32_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I8), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I32)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I32)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI8I64_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I8), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I64)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I64)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI16_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I16), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = Ix(i,j)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = Ix(i,j)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI16I8_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I16), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I8)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I8)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI16I32_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I16), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I32)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I32)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI16I64_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I16), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I64)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I64)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI32_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = Ix(i,j)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = Ix(i,j)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI32I8_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I8)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I8)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI32I16_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I16)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I16)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI32I64_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I64)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I64)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI64_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I64), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = Ix(i,j)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = Ix(i,j)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI64I8_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I64), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I8)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I8)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI64I16_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I64), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I16)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I16)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI64I32_2D (Ix,Iy,n,m)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I64), dimension(:,:), intent(IN) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      !%omp parallel do private(i)
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I32)
        end do
      end do
      !%omp end parallel do
    
    else
    
      !%omp parallel do private(i)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I32)
        end do
      end do
      !%omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorLogical2D (Lx,Ly,n,m)
  
!<description>
  ! Copies a logical vector Lx: Ly = Lx
!</description>

!<input>
  
  ! Source vector
  logical, dimension(:,:), intent(IN) :: Lx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  logical, dimension(:,:), intent(OUT) :: Ly
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if (present(n) .and. present(m)) then

      !$omp parallel do private(i)
      do j=1,m
        do i=1,n
          Ly(i,j) = Lx(i,j)
        end do
      end do
      !$omp end parallel do

    else
      
      !$omp parallel do private(i)
      do j=1,size(Lx,2)
        do i=1,size(Lx,1)
          Ly(i,j) = Lx(i,j)
        end do
      end do
      !$omp end parallel do

    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorChar2D (Sx,Sy,n,m)
  
!<description>
  ! Copies a character vector Sx: Sy = Sx
!</description>

!<input>
  
  ! Source vector
  character, dimension(:,:), intent(IN) :: Sx
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  character, dimension(:,:), intent(OUT) :: Sy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if (present(n) .and. present(m)) then

      !$omp parallel do private(i)
      do j=1,m
        do i=1,n
          Sy(i,j) = Sx(i,j)
        end do
      end do
      !$omp end parallel do

    else

      !$omp parallel do private(i)
      do j=1,size(Sx,2)
        do i=1,size(Sx,1)
          Sy(i,j) = Sx(i,j)
        end do
      end do
      !$omp end parallel do

    end if
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lalg_scaleVectorSngl (Fx,sc,n)
  
!<description>
  ! Scales a single precision vector: Dx = sc * Dx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  real(SP), dimension(:), intent(INOUT) :: Fx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(SP), intent(IN) :: sc

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>
  
!</subroutine>

    if (.not. present(n)) then
      if(sc .eq. 0.0_SP) then
        call lalg_clearVectorSngl(Fx)
      else if(sc .ne. 1.0_SP) then
        call SSCAL(size(Fx),sc,Fx,1)
      end if
    else
      if(sc .eq. 0.0_SP) then
        call lalg_clearVectorSngl(Fx,n)
      else if(sc .ne. 1.0_SP) then
        call SSCAL(n,sc,Fx,1)
      end if
    end if
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lalg_scaleVectorDble (Dx,dc,n)
  
!<description>
  ! Scales a double precision vector: Dx = dc * Dx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  real(DP), dimension(:), intent(INOUT) :: Dx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(DP), intent(IN) :: dc

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>
  
!</subroutine>

    if (.not. present(n)) then
      if(dc .eq. 0.0_DP) then
        call lalg_clearVectorDble(Dx)
      else if(dc .ne. 1.0_DP) then
        call DSCAL(size(Dx),dc,Dx,1)
      end if
    else
      if(dc .eq. 0.0_DP) then
        call lalg_clearVectorDble(Dx,n)
      else if(dc .ne. 1.0_DP) then
        call DSCAL(n,dc,Dx,1)
      end if
    end if
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lalg_scaleVectorQuad (Qx,qc,n)
  
!<description>
  ! Scales a quad precision vector: Qx = qc * Qx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  real(QP), dimension(:), intent(INOUT) :: Qx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(QP), intent(IN) :: qc

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
      if(qc .eq. 0.0_QP) then
        call lalg_clearVectorQuad(Qx)
      else if(qc .ne. 1.0_QP) then
      
        !$omp parallel do
        do i = 1, size(Qx)
          Qx(i) = qc * Qx(i)
        end do
        !$omp end parallel do
      
        !call QSCAL(size(Qx),qc,Qx,1)
      end if
    else
      if(qc .eq. 0.0_QP) then
        call lalg_clearVectorQuad(Qx,n)
      else if(qc .ne. 1.0_QP) then

        !$omp parallel do
        do i = 1, n
          Qx(i) = qc * Qx(i)
        end do
        !$omp end parallel do

        !call QSCAL(n,qc,Qx,1)
      end if
    end if
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lalg_scaleVectorSngl2D (Fx,sc)
  
!<description>
  ! Scales a single precision vector: Dx = sc * Dx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  real(SP), dimension(:,:), intent(INOUT) :: Fx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(SP), intent(IN) :: sc

!</input>
  
!</subroutine>

    if(sc .eq. 0.0_DP) then
      call lalg_clearVectorSngl2D(Fx)
    else if(sc .ne. 1.0_SP) then
      call SSCAL(size(Fx,1)*size(Fx,2),sc,Fx,1)
    end if
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lalg_scaleVectorDble2D (Dx,dc)
  
!<description>
  ! Scales a double precision vector: Dx = dc * Dx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  real(DP), dimension(:,:), intent(INOUT) :: Dx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(DP), intent(IN) :: dc

!</input>
  
!</subroutine>

    if(dc .eq. 0.0_DP) then
      call lalg_clearVectorDble2D(Dx)
    else if(dc .ne. 1.0_DP) then
      call DSCAL(size(Dx,1)*size(Dx,2),dc,Dx,1)
    end if
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lalg_scaleVectorQuad2D (Qx,qc)
  
!<description>
  ! Scales a quad precision vector: Qx = dc * Qx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  real(QP), dimension(:,:), intent(INOUT) :: Qx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(QP), intent(IN) :: qc

!</input>
  
!</subroutine>

  integer :: i,j

    if(qc .eq. 0.0_QP) then
      call lalg_clearVectorQuad2D(Qx)
    else if(qc .ne. 1.0_QP) then
    
      !$omp parallel do private(i)
      do j = 1, size(Qx,2)
        do i = 1, size(Qx,1)
          Qx(i,j) = qc*Qx(i,j)
        end do
      end do
      !$omp end parallel do
    
      !call QSCAL(size(Qx,1)*size(Qx,2),qc,Qx,1)
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorSngl (Fx,n)
  
!<description>
  ! Clears a single precision vector: Fx = 0
!</description>

!<input>
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  real(SP), dimension(:), intent(OUT) :: Fx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1,size(Fx)
        Fx(i) = 0.0_SP
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,n
        Fx(i) = 0.0_SP
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorDble (Dx,n)
  
!<description>
  ! Clears a double precision vector: Dx = 0
!</description>

!<input>
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  real(DP), dimension(:), intent(OUT) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1,size(Dx)
        Dx(i) = 0.0_DP
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,n
        Dx(i) = 0.0_DP
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorQuad (Qx,n)
  
!<description>
  ! Clears a quad precision vector: Qx = 0
!</description>

!<input>
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  real(QP), dimension(:), intent(OUT) :: Qx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1,size(Qx)
        Qx(i) = 0.0_QP
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,n
        Qx(i) = 0.0_QP
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI8 (Ix,n)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<input>
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  integer(I8), dimension(:), intent(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1,size(Ix)
        Ix(i) = 0_I8
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,size(Ix)
        Ix(i) = 0_I8
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI16 (Ix,n)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<input>
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  integer(I16), dimension(:), intent(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1,size(Ix)
        Ix(i) = 0_I16
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,size(Ix)
        Ix(i) = 0_I16
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI32 (Ix,n)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<input>
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  integer(I32), dimension(:), intent(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1,size(Ix)
        Ix(i) = 0_I32
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,size(Ix)
        Ix(i) = 0_I32
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI64 (Ix,n)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<input>
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  integer(I64), dimension(:), intent(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1,size(Ix)
        Ix(i) = 0_I64
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,size(Ix)
        Ix(i) = 0_I64
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorSngl2D (Fx)
  
!<description>
  ! Clears a single precision vector: Fx = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  real(SP), dimension(:,:), intent(OUT) :: Fx
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)
    do j = 1,size(Fx,2)
      do i = 1,size(Fx,1)
        Fx(i,j) = 0.0_SP
      end do
    end do
    !$omp end parallel do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorDble2D (Dx)
  
!<description>
  ! Clears a double precision vector: Dx = 0
!</description>

!<output>
  ! Destination vector to be cleared
  real(DP), dimension(:,:), intent(OUT) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)  
    do j = 1,size(Dx,2)
      do i = 1,size(Dx,1)
        Dx(i,j) = 0.0_DP
      end do
    end do
    !$omp end parallel do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorQuad2D (Qx)
  
!<description>
  ! Clears a quad precision vector: Qx = 0
!</description>

!<output>
  ! Destination vector to be cleared
  real(QP), dimension(:,:), intent(OUT) :: Qx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)  
    do j = 1,size(Qx,2)
      do i = 1,size(Qx,1)
        Qx(i,j) = 0.0_QP
      end do
    end do
    !$omp end parallel do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI8_2D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I8), dimension(:,:), intent(OUT) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)  
    do j = 1,size(Ix,2)
      do i = 1,size(Ix,1)
        Ix(i,j) = 0_I8
      end do
    end do
    !$omp end parallel do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI16_2D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I16), dimension(:,:), intent(OUT) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)  
    do j = 1,size(Ix,2)
      do i = 1,size(Ix,1)
        Ix(i,j) = 0_I16
      end do
    end do
    !$omp end parallel do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI32_2D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I32), dimension(:,:), intent(OUT) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)  
    do j = 1,size(Ix,2)
      do i = 1,size(Ix,1)
        Ix(i,j) = 0_I32
      end do
    end do
    !$omp end parallel do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI64_2D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I64), dimension(:,:), intent(OUT) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)  
    do j = 1,size(Ix,2)
      do i = 1,size(Ix,1)
        Ix(i,j) = 0_I32
      end do
    end do
    !$omp end parallel do
  
  end subroutine
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorSngl (Fx,fvalue,n)
  
!<description>
  ! Sets the vector data to a defined value: Fx = fvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  real(SP), intent(IN) :: fvalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  real(SP), dimension(:), intent(OUT) :: Fx
!</output>

!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      !$omp parallel do
      do i = 1,size(Fx)
        Fx(i) = fvalue
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,n
        Fx(i) = fvalue
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorDble (Dx,dvalue,n)
  
!<description>
  ! Sets the vector data to a defined value: Dx = dvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  real(DP), intent(IN) :: dvalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  real(DP), dimension(:), intent(OUT) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then

      !$omp parallel do
      do i = 1,size(Dx)
        Dx(i) = dvalue
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,size(Dx)
        Dx(i) = dvalue
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorQuad (Qx,qvalue,n)
  
!<description>
  ! Sets the vector data to a defined value: Qx = qvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  real(QP), intent(IN) :: qvalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  real(QP), dimension(:), intent(OUT) :: Qx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then

      !$omp parallel do
      do i = 1,size(Qx)
        Qx(i) = qvalue
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,size(Qx)
        Qx(i) = qvalue
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI8 (Ix,ivalue,n)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I8), intent(IN) :: ivalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  integer(I8), dimension(:), intent(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
   
      !$omp parallel do
      do i = 1,size(Ix)
        Ix(i) = ivalue
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,n
        Ix(i) = ivalue
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI16 (Ix,ivalue,n)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I16), intent(IN) :: ivalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  integer(I16), dimension(:), intent(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
   
      !$omp parallel do
      do i = 1,size(Ix)
        Ix(i) = ivalue
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,n
        Ix(i) = ivalue
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI32 (Ix,ivalue,n)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I32), intent(IN) :: ivalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  integer(I32), dimension(:), intent(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
   
      !$omp parallel do
      do i = 1,size(Ix)
        Ix(i) = ivalue
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,n
        Ix(i) = ivalue
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI64 (Ix,ivalue,n)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I64), intent(IN) :: ivalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  integer(I64), dimension(:), intent(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
   
      !$omp parallel do
      do i = 1,size(Ix)
        Ix(i) = ivalue
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,n
        Ix(i) = ivalue
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorLogical (Lx,lvalue,n)
  
!<description>
  ! Sets the vector data to a defined value: Lx = lvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  logical, intent(IN) :: lvalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  logical, dimension(:), intent(OUT) :: Lx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1,size(Lx)
        Lx(i) = lvalue
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,n
        Lx(i) = lvalue
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorChar (Sx,svalue,n)
  
!<description>
  ! Sets the vector data to a defined value: Sx = svalue
!</description>

!<input>
  ! The value, the vector should be set to.
  character, intent(IN) :: svalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  character, dimension(:), intent(OUT) :: Sx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      !$omp parallel do
      do i = 1,size(Sx)
        Sx(i) = svalue
      end do
      !$omp end parallel do
      
    else
    
      !$omp parallel do
      do i = 1,n
        Sx(i) = svalue
      end do
      !$omp end parallel do
      
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorSngl2D (Fx,fvalue)
  
!<description>
  ! Sets the vector data to a defined value: Fx = fvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  real(SP), intent(IN) :: fvalue
!</input>

!<output>
  ! Destination vector to be set
  real(SP), dimension(:,:), intent(OUT) :: Fx
!</output>

!</subroutine>

  ! local variables
  integer :: i,j
  
    !$omp parallel do private(i)
    do j = 1,size(Fx,2)
      do i = 1,size(Fx,1)
        Fx(i,j) = fvalue
      end do
    end do
    !$omp end parallel do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorDble2D (Dx,dvalue)
  
!<description>
  ! Sets the vector data to a defined value: Dx = dvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  real(DP), intent(IN) :: dvalue
!</input>

!<output>
  ! Destination vector to be set
  real(DP), dimension(:,:), intent(OUT) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)  
    do j = 1,size(Dx,2)
      do i = 1,size(Dx,1)
        Dx(i,j) = dvalue
      end do
    end do
    !$omp end parallel do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorQuad2D (Qx,qvalue)
  
!<description>
  ! Sets the vector data to a defined value: Qx = qvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  real(QP), intent(IN) :: qvalue
!</input>

!<output>
  ! Destination vector to be set
  real(QP), dimension(:,:), intent(OUT) :: Qx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)  
    do j = 1,size(Qx,2)
      do i = 1,size(Qx,1)
        Qx(i,j) = qvalue
      end do
    end do
    !$omp end parallel do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI8_2D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I8), intent(IN) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  integer(I8), dimension(:,:), intent(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)  
    do j = 1,size(Ix,2)
      do i = 1,size(Ix,1)
        Ix(i,j) = ivalue
      end do
    end do
    !$omp end parallel do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI16_2D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I16), intent(IN) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  integer(I16), dimension(:,:), intent(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)  
    do j = 1,size(Ix,2)
      do i = 1,size(Ix,1)
        Ix(i,j) = ivalue
      end do
    end do
    !$omp end parallel do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI32_2D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I32), intent(IN) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  integer(I32), dimension(:,:), intent(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)  
    do j = 1,size(Ix,2)
      do i = 1,size(Ix,1)
        Ix(i,j) = ivalue
      end do
    end do
    !$omp end parallel do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI64_2D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I64), intent(IN) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  integer(I64), dimension(:,:), intent(OUT) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)  
    do j = 1,size(Ix,2)
      do i = 1,size(Ix,1)
        Ix(i,j) = ivalue
      end do
    end do
    !$omp end parallel do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorLogical2D (Lx,lvalue)
  
!<description>
  ! Sets the vector data to a defined value: Lx = lvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  logical, intent(IN) :: lvalue
!</input>

!<output>
  ! Destination vector to be set
  logical, dimension(:,:), intent(OUT) :: Lx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j
  
    !$omp parallel do private(i)
    do j = 1,size(Lx,2)
      do i = 1,size(Lx,1)
        Lx(i,j) = lvalue
      end do
    end do
    !$omp end parallel do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorChar2D (Sx,svalue)
  
!<description>
  ! Sets the vector data to a defined value: Sx = svalue
!</description>

!<input>
  ! The value, the vector should be set to.
  character, intent(IN) :: svalue
!</input>

!<output>
  ! Destination vector to be set
  character, dimension(:,:), intent(OUT) :: Sx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    !$omp parallel do private(i)
    do j = 1,size(Sx,2)
      do i = 1,size(Sx,1)
        Sx(i,j) = svalue
      end do
    end do
    !$omp end parallel do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSngl (Fx,Fy,scx,scy,n)
  
!<description>
  ! Performs a linear combination: Fy = scx * Fx  +  scy * Fy
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:), intent(IN) :: Fx
  
  ! Scaling factor for Dx
  real(SP), intent(IN)               :: scx

  ! Scaling factor for Dy
  real(SP), intent(IN)               :: scy

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(SP), dimension(:), intent(INOUT) :: Fy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(SP) :: c
  
    if (.not. present(n)) then

      if (scy .eq. 0.0_SP) then
        call SCOPY(size(Fx),Fx,1,Fy,1)
        if (scx .ne. 1.0_SP) call SSCAL(size(Fx),scx,Fy,1)
      else if (scy .eq. 1.0_DP) then
        call SAXPY(size(Fx),scx,Fx,1,Fy,1)
      else
        c=scx/scy
        call SAXPY(size(Fx),c,Fx,1,Fy,1)
        call SSCAL(size(Fx),scy,Fy,1)
      endif
      
    else
    
      if (scy .eq. 0.0_SP) then
        call SCOPY(n,Fx,1,Fy,1)
        if (scx .ne. 1.0_SP) call SSCAL(size(Fx),scx,Fy,1)
      else if (scy .eq. 1.0_DP) then
        call SAXPY(n,scx,Fx,1,Fy,1)
      else
        c=scx/scy
        call SAXPY(n,c,Fx,1,Fy,1)
        call SSCAL(n,scy,Fy,1)
      endif

    end if  
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombDble (Dx,Dy,dcx,dcy,n)
  
!<description>
  ! Performs a linear combination: Dy = dcx * Dx  +  dcy * Dy
!</description>

!<input>
  
  ! First source vector
  real(DP), dimension(:), intent(IN) :: Dx
  
  ! Scaling factor for Dx
  real(DP), intent(IN)               :: dcx

  ! Scaling factor for Dy
  real(DP), intent(IN)               :: dcy
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(DP), dimension(:), intent(INOUT) :: Dy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(DP) :: c
  
    if (.not. present(n)) then
    
      if (dcy .eq. 0.0_DP) then
        call DCOPY(size(Dx),Dx,1,Dy,1)
        if (dcx .ne. 1.0_DP) call DSCAL(size(Dx),dcx,Dy,1)
      else if (dcy .eq. 1.0_DP) then
        call DAXPY(size(Dx),dcx,Dx,1,Dy,1)
      else
        c=dcx/dcy
        call DAXPY(size(Dx),c,Dx,1,Dy,1)
        call DSCAL(size(Dx),dcy,Dy,1)
      endif
      
    else
    
      if (dcy .eq. 0.0_DP) then
        call DCOPY(n,Dx,1,Dy,1)
        if (dcx .ne. 1.0_DP) call DSCAL(size(Dx),dcx,Dy,1)
      else if (dcy .eq. 1.0_DP) then
        call DAXPY(n,dcx,Dx,1,Dy,1)
      else
        c=dcx/dcy
        call DAXPY(n,c,Dx,1,Dy,1)
        call DSCAL(n,dcy,Dy,1)
      endif
      
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombQuad (Qx,Qy,qcx,qcy,n)
  
!<description>
  ! Performs a linear combination: Qy = qcx * Qx  +  qcy * Qy
!</description>

!<input>
  
  ! First source vector
  real(QP), dimension(:), intent(IN) :: Qx
  
  ! Scaling factor for Dx
  real(QP), intent(IN)               :: qcx

  ! Scaling factor for Dy
  real(QP), intent(IN)               :: qcy
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:), intent(INOUT) :: Qy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  ! local variables
  integer :: i,k
  
    if (.not. present(n)) then
      k = size(Qx)
    else
      k = n
    end if
  
    if(qcx .eq. 0.0_QP) then
    
      if(qcy .eq. 1.0_QP) then
        
        ! Nothing to do
        return
      
      else if(qcy .eq. 0.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = 0.0_QP
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = qcy*Qy(i)
        end do
        !$omp end parallel do
      
      end if
    
    else if(qcx .eq. 1.0_QP) then
    
      if(qcy .eq. 1.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = Qy(i) + Qx(i)
        end do
        !$omp end parallel do
      
      else if(qcy .eq. 0.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = Qx(i)
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = qcy*Qy(i) + Qx(i)
        end do
        !$omp end parallel do

      end if
    
    else
    
      if(qcy .eq. 1.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = Qy(i) + qcx*Qx(i)
        end do
        !$omp end parallel do
      
      else if(qcy .eq. 0.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = qcx*Qx(i)
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = qcy*Qy(i) + qcx*Qx(i)
        end do
        !$omp end parallel do

      end if
    
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSnglDble (Fx,Dy,scx,dcy,n)
  
!<description>
  ! Performs a linear combination: Dy = scx * Fx  +  dcy * Dy
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:), intent(IN) :: Fx
  
  ! Scaling factor for Dx
  real(SP), intent(IN)               :: scx

  ! Scaling factor for Dy
  real(DP), intent(IN)               :: dcy
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(DP), dimension(:), intent(INOUT) :: Dy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i,k
  real(DP) :: dcx
  
    if (.not. present(n)) then
      k = size(Fx)
    else
      k = n
    end if
  
    if(scx .eq. 0.0_SP) then
    
      if(dcy .eq. 1.0_DP) then
        
        ! Nothing to do
        return
      
      else if(dcy .eq. 0.0_DP) then
      
        ! Simply clear Dy
        !$omp parallel do
        do i = 1, k
          Dy(i) = 0.0_DP
        end do
        !$omp end parallel do
      
      else
      
        ! Call DSCAL
        call DSCAL(k,dcy,Dy,1)
      
      end if
    
    else if(scx .eq. 1.0_SP) then
    
      if(dcy .eq. 1.0_DP) then
      
        !$omp parallel do
        do i = 1, k
          Dy(i) = Dy(i) + real(Fx(i),DP)
        end do
        !$omp end parallel do
      
      else if(dcy .eq. 0.0_DP) then
      
        !$omp parallel do
        do i = 1, k
          Dy(i) = real(Fx(i),DP)
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do
        do i = 1, k
          Dy(i) = dcy*Dy(i) + real(Fx(i),DP)
        end do
        !$omp end parallel do

      end if
    
    else
    
      ! Convert scx to double precision
      dcx = real(scx,DP)

      if(dcy .eq. 1.0_DP) then
      
        !$omp parallel do
        do i = 1, k
          Dy(i) = Dy(i) + dcx*real(Fx(i),DP)
        end do
        !$omp end parallel do
      
      else if(dcy .eq. 0.0_DP) then
      
        !$omp parallel do
        do i = 1, k
          Dy(i) = dcx*real(Fx(i),DP)
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do
        do i = 1, k
          Dy(i) = dcy*Dy(i) + dcx*real(Fx(i),DP)
        end do
        !$omp end parallel do

      end if
    
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSnglQuad (Fx,Qy,scx,qcy,n)
  
!<description>
  ! Performs a linear combination: Qy = scx * Fx  +  qcy * Qy
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:), intent(IN) :: Fx
  
  ! Scaling factor for Dx
  real(SP), intent(IN)               :: scx

  ! Scaling factor for Qy
  real(QP), intent(IN)               :: qcy
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:), intent(INOUT) :: Qy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i,k
  real(QP) :: qcx
  
    if (.not. present(n)) then
      k = size(Fx)
    else
      k = n
    end if
  
    if(scx .eq. 0.0_SP) then
    
      if(qcy .eq. 1.0_QP) then
        
        ! Nothing to do
        return
      
      else if(qcy .eq. 0.0_QP) then
      
        ! Simply clear Qy
        !$omp parallel do
        do i = 1, k
          Qy(i) = 0.0_QP
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = qcy*Qy(i)
        end do
        !$omp end parallel do
      
      end if
    
    else if(scx .eq. 1.0_SP) then
    
      if(qcy .eq. 1.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = Qy(i) + real(Fx(i),QP)
        end do
        !$omp end parallel do
      
      else if(qcy .eq. 0.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = real(Fx(i),QP)
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = qcy*Qy(i) + real(Fx(i),QP)
        end do
        !$omp end parallel do

      end if
    
    else
    
      ! Convert scx to quad precision
      qcx = real(scx,QP)

      if(qcy .eq. 1.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = Qy(i) + qcx*real(Fx(i),QP)
        end do
        !$omp end parallel do
      
      else if(qcy .eq. 0.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = qcx*real(Fx(i),QP)
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = qcy*Qy(i) + qcx*real(Fx(i),QP)
        end do
        !$omp end parallel do

      end if
    
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombDblQuad (Dx,Qy,dcx,qcy,n)
  
!<description>
  ! Performs a linear combination: Qy = dcx * Dx  +  qcy * Qy
!</description>

!<input>
  
  ! First source vector
  real(DP), dimension(:), intent(IN) :: Dx
  
  ! Scaling factor for Dx
  real(DP), intent(IN)               :: dcx

  ! Scaling factor for Qy
  real(QP), intent(IN)               :: qcy
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:), intent(INOUT) :: Qy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i,k
  real(QP) :: qcx
  
    if (.not. present(n)) then
      k = size(Dx)
    else
      k = n
    end if
  
    if(dcx .eq. 0.0_DP) then
    
      if(qcy .eq. 1.0_QP) then
        
        ! Nothing to do
        return
      
      else if(qcy .eq. 0.0_QP) then
      
        ! Simply clear Qy
        !$omp parallel do
        do i = 1, k
          Qy(i) = 0.0_QP
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = qcy*Qy(i)
        end do
        !$omp end parallel do
      
      end if
    
    else if(dcx .eq. 1.0_DP) then
    
      if(qcy .eq. 1.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = Qy(i) + real(Dx(i),QP)
        end do
        !$omp end parallel do
      
      else if(qcy .eq. 0.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = real(Dx(i),QP)
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = qcy*Qy(i) + real(Dx(i),QP)
        end do
        !$omp end parallel do

      end if
    
    else
    
      ! Convert dcx to quad precision
      qcx = real(dcx,QP)

      if(qcy .eq. 1.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = Qy(i) + qcx*real(Dx(i),QP)
        end do
        !$omp end parallel do
      
      else if(qcy .eq. 0.0_QP) then
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = qcx*real(Dx(i),QP)
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do
        do i = 1, k
          Qy(i) = qcy*Qy(i) + qcx*real(Dx(i),QP)
        end do
        !$omp end parallel do

      end if
    
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSngl2D (Fx,Fy,scx,scy)
  
!<description>
  ! Performs a linear combination: Fy = scx * Fx  +  scy * Fy
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:,:), intent(IN) :: Fx
  
  ! Scaling factor for Dx
  real(SP), intent(IN)                 :: scx

  ! Scaling factor for Dy
  real(SP), intent(IN)                 :: scy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(SP), dimension(:,:), intent(INOUT) :: Fy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(SP) :: c
  
    if (scy .eq. 0.0_SP) then
      call SCOPY(size(Fx,1)*size(Fx,2),Fx,1,Fy,1)
      if (scx .ne. 1.0_SP) call SSCAL(size(Fx,1)*size(Fx,2),scx,Fy,1)
    else if (scy .eq. 1.0_DP) then
      call SAXPY(size(Fx,1)*size(Fx,2),scx,Fx,1,Fy,1)
    else
      c=scx/scy
      call SAXPY(size(Fx,1)*size(Fx,2),c,Fx,1,Fy,1)
      call SSCAL(size(Fx,1)*size(Fx,2),scy,Fy,1)
    endif
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombDble2D (Dx,Dy,dcx,dcy)
  
!<description>
  ! Performs a linear combination: Dy = dcx * Dx  +  dcy * Dy
!</description>

!<input>
  
  ! First source vector
  real(DP), dimension(:,:), intent(IN) :: Dx
  
  ! Scaling factor for Dx
  real(DP), intent(IN)                 :: dcx

  ! Scaling factor for Dy
  real(DP), intent(IN)                 :: dcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(DP), dimension(:,:), intent(INOUT) :: Dy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(DP) :: c
  
    if (dcy .eq. 0.0_DP) then
      call DCOPY(size(Dx,1)*size(Dx,2),Dx,1,Dy,1)
      if (dcx .ne. 1.0_DP) call DSCAL(size(Dx,1)*size(Dx,2),dcx,Dy,1)
    else if (dcy .eq. 1.0_DP) then
      call DAXPY(size(Dx,1)*size(Dx,2),dcx,Dx,1,Dy,1)
    else
      c=dcx/dcy
      call DAXPY(size(Dx,1)*size(Dx,2),c,Dx,1,Dy,1)
      call DSCAL(size(Dx,1)*size(Dx,2),dcy,Dy,1)
    endif
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombQuad2D (Qx,Qy,qcx,qcy)
  
!<description>
  ! Performs a linear combination: Qy = qcx * Qx  +  qcy * Qy
!</description>

!<input>
  
  ! First source vector
  real(QP), dimension(:,:), intent(IN) :: Qx
  
  ! Scaling factor for Qx
  real(QP), intent(IN)                 :: qcx

  ! Scaling factor for Qy
  real(QP), intent(IN)                 :: qcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:,:), intent(INOUT) :: Qy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i,j
  
    !$omp parallel do private(i)
    do j = 1, size(Qx,2)
      do i = 1, size(Qx,1)
        Qy(i,j) = qcy*Qy(i,j) + qcx*Qx(i,j)
      end do
    end do
    !$omp end parallel do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSnglDble2D (Fx,Dy,scx,dcy)
  
!<description>
  ! Performs a linear combination: Dy = scx * Fx  +  dcy * Dy
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:,:), intent(IN) :: Fx
  
  ! Scaling factor for Dx
  real(SP), intent(IN)                 :: scx

  ! Scaling factor for Dy
  real(DP), intent(IN)                 :: dcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(DP), dimension(:,:), intent(INOUT) :: Dy
  
!</inputoutput>
  
!</subroutine>

  integer :: i,j
  real(DP) :: dcx
  
    if(scx .eq. 0.0_SP) then
    
      if(dcy .eq. 1.0_DP) then
        
        ! Nothing to do
        return
      
      else if(dcy .eq. 0.0_DP) then
      
        ! Simply clear Dy
        !$omp parallel do private(i)
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = 0.0_DP
          end do
        end do
        !$omp end parallel do
      
      else
      
        ! Call DSCAL
        call DSCAL(size(Dy,1)*size(Dy,2),dcy,Dy,1)
      
      end if
    
    else if(scx .eq. 1.0_SP) then
    
      if(dcy .eq. 1.0_DP) then
      
        !$omp parallel do private(i)
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = Dy(i,j) + real(Fx(i,j),DP)
          end do
        end do
        !$omp end parallel do
      
      else if(dcy .eq. 0.0_DP) then
      
        !$omp parallel do private(i)
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = real(Fx(i,j),DP)
          end do
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do private(i)
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = dcy*Dy(i,j) + real(Fx(i,j),DP)
          end do
        end do
        !$omp end parallel do

      end if
    
    else
    
      ! Convert scx to double precision
      dcx = real(scx,DP)

      if(dcy .eq. 1.0_DP) then
      
        !$omp parallel do private(i)
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = Dy(i,j) + dcx*real(Fx(i,j),DP)
          end do
        end do
        !$omp end parallel do
      
      else if(dcy .eq. 0.0_DP) then
      
        !$omp parallel do private(i)
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = dcx*real(Fx(i,j),DP)
          end do
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do private(i)
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = dcy*Dy(i,j) + dcx*real(Fx(i,j),DP)
          end do
        end do
        !$omp end parallel do

      end if
    
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSnglQuad2D (Fx,Qy,scx,qcy)
  
!<description>
  ! Performs a linear combination: Qy = scx * Fx  +  qcy * Qy
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:,:), intent(IN) :: Fx
  
  ! Scaling factor for Dx
  real(SP), intent(IN)                 :: scx

  ! Scaling factor for Qy
  real(QP), intent(IN)                 :: qcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:,:), intent(INOUT) :: Qy
  
!</inputoutput>
  
!</subroutine>

  integer :: i,j
  real(QP) :: qcx
  
    if(scx .eq. 0.0_SP) then
    
      if(qcy .eq. 1.0_QP) then
        
        ! Nothing to do
        return
      
      else if(qcy .eq. 0.0_QP) then
      
        ! Simply clear Qy
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = 0.0_QP
          end do
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcy*Qy(i,j)
          end do
        end do
        !$omp end parallel do
      
      end if
    
    else if(scx .eq. 1.0_SP) then
    
      if(qcy .eq. 1.0_QP) then
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = Qy(i,j) + real(Fx(i,j),QP)
          end do
        end do
        !$omp end parallel do
      
      else if(qcy .eq. 0.0_QP) then
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = real(Fx(i,j),QP)
          end do
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcy*Qy(i,j) + real(Fx(i,j),QP)
          end do
        end do
        !$omp end parallel do

      end if
    
    else
    
      ! Convert scx to quad precision
      qcx = real(scx,QP)

      if(qcy .eq. 1.0_QP) then
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = Qy(i,j) + qcx*real(Fx(i,j),QP)
          end do
        end do
        !$omp end parallel do
      
      else if(qcy .eq. 0.0_QP) then
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcx*real(Fx(i,j),QP)
          end do
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcy*Qy(i,j) + qcx*real(Fx(i,j),QP)
          end do
        end do
        !$omp end parallel do

      end if
    
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombDblQuad2D (Dx,Qy,dcx,qcy)
  
!<description>
  ! Performs a linear combination: Qy = dcx * Dx  +  qcy * Qy
!</description>

!<input>
  
  ! First source vector
  real(DP), dimension(:,:), intent(IN) :: Dx
  
  ! Scaling factor for Dx
  real(DP), intent(IN)                 :: dcx

  ! Scaling factor for Qy
  real(QP), intent(IN)                 :: qcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:,:), intent(INOUT) :: Qy
  
!</inputoutput>
  
!</subroutine>

  integer :: i,j
  real(QP) :: qcx
  
    if(dcx .eq. 0.0_DP) then
    
      if(qcy .eq. 1.0_QP) then
        
        ! Nothing to do
        return
      
      else if(qcy .eq. 0.0_QP) then
      
        ! Simply clear Qy
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = 0.0_QP
          end do
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcy*Qy(i,j)
          end do
        end do
        !$omp end parallel do
      
      end if
    
    else if(dcx .eq. 1.0_DP) then
    
      if(qcy .eq. 1.0_QP) then
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = Qy(i,j) + real(Dx(i,j),QP)
          end do
        end do
        !$omp end parallel do
      
      else if(qcy .eq. 0.0_QP) then
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = real(Dx(i,j),QP)
          end do
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcy*Qy(i,j) + real(Dx(i,j),QP)
          end do
        end do
        !$omp end parallel do

      end if
    
    else
    
      ! Convert dcx to quad precision
      qcx = real(dcx,QP)

      if(qcy .eq. 1.0_QP) then
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = Qy(i,j) + qcx*real(Dx(i,j),QP)
          end do
        end do
        !$omp end parallel do
      
      else if(qcy .eq. 0.0_QP) then
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcx*real(Dx(i,j),QP)
          end do
        end do
        !$omp end parallel do
      
      else
      
        !$omp parallel do private(i)
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcy*Qy(i,j) + qcx*real(Dx(i,j),QP)
          end do
        end do
        !$omp end parallel do

      end if
    
    end if
  
  end subroutine

  ! ***************************************************************************

!<function>

  real(SP) function lalg_scalarProductSngl (Fx,Fy,n) result (res)
  
!<description>
  ! Calculates the scalar product of two single precision vectors: 
  ! res = (vector,vector)
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:), intent(IN) :: Fx
  
  ! Second source vector
  real(SP), dimension(:), intent(IN) :: Fy

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  real(SP) :: SDOT

    if (.not. present(n)) then
      res = SDOT(size(Fx),Fx,1,Fy,1)
    else
      res = SDOT(size(Fx),Fx,1,Fy,1)
    end if
  
  end function

  ! ***************************************************************************

!<function>

  real(DP) function lalg_scalarProductDble (Dx,Dy,n) result (res)
  
!<description>
  ! Calculates the scalar product of two double precision vectors: 
  ! res = (vector,vector)
!</description>

!<input>
  
  ! First source vector
  real(DP), dimension(:), intent(IN) :: Dx
  
  ! Second source vector
  real(DP), dimension(:), intent(IN) :: Dy

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  real(DP) :: DDOT

    if (.not. present(n)) then
      res = DDOT(size(Dx),Dx,1,Dy,1)
    else
      res = DDOT(size(Dx),Dx,1,Dy,1)
    end if
  
  end function
  
  ! ***************************************************************************

!<function>

  real(QP) function lalg_scalarProductQuad (Qx,Qy,n) result (res)
  
!<description>
  ! Calculates the scalar product of two quad precision vectors: 
  ! res = (vector,vector)
!</description>

!<input>
  
  ! First source vector
  real(QP), dimension(:), intent(IN) :: Qx
  
  ! Second source vector
  real(QP), dimension(:), intent(IN) :: Qy
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  !real(DP) :: QDOT
  integer :: i

    res = 0.0_QP
    if (.not. present(n)) then
    
      !$omp parallel do reduction(+:res)
      do i = 1, size(Qx)
        res = res + Qx(i)*Qy(i)
      end do
      !$omp end parallel do
      
      !res = QDOT(size(Qx),Qx,1,Qy,1)
    else
    
      !$omp parallel do reduction(+:res)
      do i = 1, n
        res = res + Qx(i)*Qy(i)
      end do
      !$omp end parallel do
      
      !res = QDOT(n,Qx,1,Qy,1)
    end if
  
  end function

  ! ***************************************************************************

!<function>

  real(SP) function lalg_scalarProductSngl2D (Fx,Fy) result (res)
  
!<description>
  ! Calculates the scalar product of two single precision vectors: 
  ! res = (vector,vector)
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:,:), intent(IN) :: Fx
  
  ! Second source vector
  real(SP), dimension(:,:), intent(IN) :: Fy
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  real(SP) :: SDOT

    res = SDOT(size(Fx,1)*size(Fx,2),Fx,1,Fy,1)
  
  end function

  ! ***************************************************************************

!<function>

  real(DP) function lalg_scalarProductDble2D (Dx,Dy) result (res)
  
!<description>
  ! Calculates the scalar product of two double precision vectors: 
  ! res = (vector,vector)
!</description>

!<input>
  
  ! First source vector
  real(DP), dimension(:,:), intent(IN) :: Dx
  
  ! Second source vector
  real(DP), dimension(:,:), intent(IN) :: Dy
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  real(DP) :: DDOT

    res = DDOT(size(Dx,1)*size(Dx,2),Dx,1,Dy,1)
  
  end function

  ! ***************************************************************************

!<function>

  real(DP) function lalg_scalarProductQuad2D (Qx,Qy) result (res)
  
!<description>
  ! Calculates the scalar product of two quad precision vectors: 
  ! res = (vector,vector)
!</description>

!<input>
  
  ! First source vector
  real(QP), dimension(:,:), intent(IN) :: Qx
  
  ! Second source vector
  real(QP), dimension(:,:), intent(IN) :: Qy
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  !real(QP) :: QDOT
  integer :: i,j
  
    res = 0.0_QP
  
    !$omp parallel do private(i) reduction(+:res)
    do j = 1, size(Qx,2)
      do i = 1, size(Qx,1)
        res = res + Qx(i,j)*Qy(i,j)
      end do
    end do
    !$omp end parallel do

    !res = QDOT(size(Qx,1)*size(Qx,2),Qx,1,Qy,1)
  
  end function

  ! ***************************************************************************

!<subroutine>

  real(SP) function lalg_normSngl (Fx,cnorm,iposMax,n) result(resnorm)
  
!<description>
  ! Calculates the norm of a single precision vector. cnorm identifies the 
  ! type of norm to calculate.
!</description>

!<input>
  ! Vector to calculate the norm of.
  real(SP), dimension(:), intent(IN) :: Fx
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(IN) :: cnorm

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(OUT), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! -1, if an error occurred (unknown norm).
!</result>
  
!</subroutine>

  real(SP) :: SASUM, SNRM2
  integer :: ISAMAX

  integer :: isize, iposmaxlocal
  
    isize = size(Fx)
    if (present(n)) isize = n

    ! Choose the norm to calculate
    select case (cnorm)
    case (LINALG_NORMSUM)
      ! L1-norm: sum absolute value of all entries
      resnorm = SASUM(isize,Fx,1)

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      resnorm = SNRM2(isize,Fx,1)

    case (LINALG_NORML1)
      ! L1-norm: sum sum absolute value of all entries, divide by sqrt(vector length).
      ! So, scale such that the vektor (1111...) to has norm = 1.
      resnorm = SASUM(isize,Fx,1) / real(isize,SP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by sqrt(vector length).
      ! So, scale such that the vektor (1111...) to has norm = 1.
      resnorm = SNRM2(isize,Fx,1) / sqrt(real(isize,SP))
      
    case (LINALG_NORMMAX)
      ! MAX-norm. Find the absolute largest entry.
      ! With the BLAS routine, calculate the position. Then get the entry.
      iposmaxlocal = ISAMAX(isize,Fx,1)
      resnorm = abs(Fx(iposmaxlocal))
      
      if(present(iposMax)) then
        iposMax = iposmaxlocal
      end if

    case default
      ! Unknown norm
      resnorm = -1.0_SP
      
    end select
      
  end function

  ! ***************************************************************************

!<subroutine>

  real(DP) function lalg_normDble (Dx,cnorm,iposMax,n) result(resnorm)
  
!<description>
  ! Calculates the norm of a double precision vector. cnorm identifies the 
  ! type of norm to calculate.
!</description>

!<input>
  ! Vector to calculate the norm of.
  real(DP), dimension(:), intent(IN) :: Dx
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(IN) :: cnorm

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(OUT), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! -1, if an error occurred (unknown norm).
!</result>
  
!</subroutine>

  real(DP) :: DASUM,DNRM2
  integer :: IDAMAX

  integer :: isize,iposmaxlocal
  
    isize = size(Dx)
    if (present(n)) isize = n

    ! Choose the norm to calculate
    select case (cnorm)
    case (LINALG_NORMSUM)
      ! L1-norm: sum all entries
      resnorm = DASUM(isize,Dx,1)

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      resnorm = DNRM2(isize,Dx,1)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vektor (1111...) to has norm = 1.
      resnorm = DASUM(isize,Dx,1) / real(isize,DP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vektor (1111...) to has norm = 1.
      resnorm = DNRM2(isize,Dx,1) / sqrt(real(isize,DP))
      
    case (LINALG_NORMMAX)
      ! MAX-norm. Find the absolute largest entry.
      ! With the BLAS routine, calculate the position. Then get the entry.
      iposmaxlocal = IDAMAX(isize,Dx,1)
      resnorm = abs(Dx(iposmaxlocal))
      
      if(present(iposMax)) then
        iposMax = iposmaxlocal
      end if
        
    case default
      ! Unknown norm
      resnorm = -1.0_DP
      
    end select
    
  end function

  ! ***************************************************************************

!<subroutine>

  real(QP) function lalg_normQuad (Qx,cnorm,iposMax,n) result(resnorm)
  
!<description>
  ! Calculates the norm of a quad precision vector. cnorm identifies the 
  ! type of norm to calculate.
!</description>

!<input>
  ! Vector to calculate the norm of.
  real(QP), dimension(:), intent(IN) :: Qx
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(IN) :: cnorm

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(OUT), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! -1, if an error occurred (unknown norm).
!</result>
  
!</subroutine>

  !real(QP) :: QASUM,QNRM2,IQAMAX

  integer :: i,isize
  
    isize = size(Qx)
    if (present(n)) isize = n
    
    resnorm = 0.0_QP

    ! Choose the norm to calculate
    select case (cnorm)
    case (LINALG_NORMSUM)
      ! L1-norm: sum all entries
      !$omp parallel do reduction(+:resnorm)
      do i = 1, isize
        resnorm = resnorm + abs(Qx(i))
      end do
      !$omp end parallel do
      !resnorm = QASUM(isize,Qx,1)

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      !$omp parallel do reduction(+:resnorm)
      do i = 1, isize
        resnorm = resnorm + Qx(i)*Qx(i)
      end do
      !$omp end parallel do
      resnorm = sqrt(resnorm)
      !resnorm = QNRM2(isize,Qx,1)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vektor (1111...) to has norm = 1.
      !$omp parallel do reduction(+:resnorm)
      do i = 1, isize
        resnorm = resnorm + abs(Qx(i))
      end do
      !$omp end parallel do
      resnorm = resnorm / real(isize,QP)
      !resnorm = QASUM(isize,Dx,1) / real(isize,QP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vektor (1111...) to has norm = 1.
      !$omp parallel do reduction(+:resnorm)
      do i = 1, isize
        resnorm = resnorm + Qx(i)*Qx(i)
      end do
      !$omp end parallel do
      resnorm = sqrt(resnorm / real(isize,QP))
      !resnorm = QNRM2(isize,Qx,1) / sqrt(real(isize,QP))
      
    case (LINALG_NORMMAX)
      ! MAX-norm. Find the absolute largest entry.
      if(present(iposMax)) then
        ! iposMax is needed - calculate norm by hand instead of calling
        ! the BLAS routine.
        resnorm = abs(Qx(1))
        iposMax = 1
        do i=2,isize
          if (abs(Qx(i)) .gt. resnorm) then
            iposMax = i
            resnorm = abs(Qx(i))
          end if
        end do
      else
        !$omp parallel do reduction(max:resnorm)
        do i = 1, isize
          resnorm = max(resnorm,abs(Qx(i)))
        end do
        !$omp end parallel do
      end if
        
    case default
      ! Unknown norm
      resnorm = -1.0_QP
      
    end select
    
  end function

  ! ***************************************************************************

!<subroutine>

  real(SP) function lalg_errorNormSngl (Fx,Fy,cnorm,iposMax,n) result(resnorm)
  
!<description>
  ! Calculates the norm of two double precision vectors, !!Fx-Fy!!
  ! cnorm identifies the type of norm to calculate.
!</description>

!<input>
  ! Vectors to calculate the norm of their difference
  real(SP), dimension(:), intent(IN) :: Fx,Fy
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(IN) :: cnorm

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(OUT), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! < 0 >1, if an error occurred (unknown norm).
!</result>
  
!</subroutine>

  real(SP) :: stemp
  integer :: i,j,isize
  
    isize = size(Fx)
    if (present(n)) isize = n
    
    resnorm = 0.0_SP

    ! Choose the norm to calculate
    select case (cnorm)
    case (LINALG_NORMSUM)
      ! L1-norm: sum all entries
      !$omp parallel do reduction(+:resnorm)
      do i = 1, isize
        resnorm = resnorm + abs(Fx(i)-Fy(i))
      end do
      !$omp end parallel do

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      !$omp parallel do reduction(+:resnorm)
      do i = 1, isize
        resnorm = resnorm + (Fx(i)-Fy(i))*(Fx(i)-Fy(i))
      end do
      !$omp end parallel do
      resnorm = sqrt(resnorm)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vector (1111...) to has norm = 1.
      !$omp parallel do reduction(+:resnorm)
      do i = 1, isize
        resnorm = resnorm + abs(Fx(i)-Fy(i))
      end do
      !$omp end parallel do
      resnorm = resnorm / real(isize,SP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vektor (1111...) to has norm = 1.
      !$omp parallel do reduction(+:resnorm)
      do i = 1, isize
        resnorm = resnorm + (Fx(i)-Fy(i))*(Fx(i)-Fy(i))
      end do
      !$omp end parallel do
      resnorm = sqrt(resnorm / real(isize,SP))
      
    case (LINALG_NORMMAX)
      ! MAX-norm. Find the absolute largest entry.
      resnorm = abs(Fx(1)-Fy(1))
      j=1
      do i=2,isize
        stemp = abs(Fx(i)-Fy(i))
        if (stemp .gt. resnorm) then
          j = i
          resnorm = stemp
        end if
      end do
      if (present(iposMax)) iposMax = j
      
    case default
      resnorm = -1.0_DP ! Unknown norm.
    end select
    
  end function

  ! ***************************************************************************

!<subroutine>

  real(DP) function lalg_errorNormDble (Dx,Dy,cnorm,iposMax,n,Dw) result(resnorm)
  
!<description>
  ! Calculates the norm of two double precision vectors, !!Dx-Dy!!
  ! cnorm identifies the type of norm to calculate.
  ! The optional parameter Dw can be used to specify a weighting
  ! vector, e.g., the lumped mass matrix, by which each component
  ! is scaled before computing the sum of contributions.
!</description>

!<input>
  ! Vectors to calculate the norm of their difference
  real(DP), dimension(:), intent(IN) :: Dx,Dy
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(IN) :: cnorm

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

  ! OPTIONAL: Weighting vector
  real(DP), dimension(:), intent(IN), optional :: Dw
!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(OUT), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! -1, if an error occurred (unknown norm).
!</result>
  
!</subroutine>

  real(DP) :: dtemp
  integer :: i,j,isize
  
    isize = size(Dx)
    if (present(n)) isize = n
    
    resnorm = 0.0_DP

    ! Choose the norm to calculate
    select case (cnorm)
    case (LINALG_NORMSUM)
      ! L1-norm: sum all entries
      if (present(Dw)) then
        !$omp parallel do reduction(+:resnorm)
        do i = 1, isize
          resnorm = resnorm + Dw(i)*abs(Dx(i)-Dy(i))
        end do
        !$omp end parallel do
      else
        !$omp parallel do reduction(+:resnorm)
        do i = 1, isize
          resnorm = resnorm + abs(Dx(i)-Dy(i))
        end do
        !$omp end parallel do
      end if

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      if (present(Dw)) then
        !$omp parallel do reduction(+:resnorm)
        do i = 1, isize
          resnorm = resnorm + Dw(i)*(Dx(i)-Dy(i))*(Dx(i)-Dy(i))
        end do
        !$omp end parallel do
      else
        !$omp parallel do reduction(+:resnorm)
        do i = 1, isize
          resnorm = resnorm + (Dx(i)-Dy(i))*(Dx(i)-Dy(i))
        end do
        !$omp end parallel do
      end if
      resnorm = sqrt(resnorm)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vektor (1111...) to has norm = 1.
      !$omp parallel do reduction(+:resnorm)
      do i = 1, isize
        resnorm = resnorm + abs(Dx(i)-Dy(i))
      end do
      !$omp end parallel do
      resnorm = resnorm / real(isize,DP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vektor (1111...) to has norm = 1.
      !$omp parallel do reduction(+:resnorm)
      do i = 1, isize
        resnorm = resnorm + (Dx(i)-Dy(i))*(Dx(i)-Dy(i))
      end do
      !$omp end parallel do
      resnorm = sqrt(resnorm / real(isize,DP))
      
    case (LINALG_NORMMAX)
      ! MAX-norm. Find the absolute largest entry.
      resnorm = abs(Dx(1)-Dy(1))
      j=1
      do i=2,isize
        dtemp = abs(Dx(i)-Dy(i))
        if (dtemp .gt. resnorm) then
          j = i
          resnorm = dtemp
        end if
      end do
      if (present(iposMax)) iposMax = j
      
    case default
      resnorm = -1.0_DP ! Unknown norm.
    end select
    
  end function

  ! ***************************************************************************

!<subroutine>

  real(QP) function lalg_errorNormQuad (Qx,Qy,cnorm,iposMax,n,Qw) result(resnorm)
  
!<description>
  ! Calculates the norm of two quad precision vectors, !!Qx-Qy!!
  ! cnorm identifies the type of norm to calculate.
  ! The optional parameter Qw can be used to specify a weighting
  ! vector, e.g., the lumped mass matrix, by which each component
  ! is scaled before computing the sum of contributions.
!</description>

!<input>
  ! Vectors to calculate the norm of their difference
  real(QP), dimension(:), intent(IN) :: Qx,Qy
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(IN) :: cnorm

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

  ! OPTIONAL: Weighting vector
  real(QP), dimension(:), intent(IN), optional :: Qw
!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(OUT), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! -1, if an error occurred (unknown norm).
!</result>
  
!</subroutine>

  real(QP) :: qtemp
  integer :: i,j,isize
  
    isize = size(Qx)
    if (present(n)) isize = n
    
    resnorm = 0.0_DP

    ! Choose the norm to calculate
    select case (cnorm)
    case (LINALG_NORMSUM)
      ! L1-norm: sum all entries
      if (present(Qw)) then
        !$omp parallel do reduction(+:resnorm)
        do i = 1, isize
          resnorm = resnorm + Qw(i)*abs(Qx(i)-Qy(i))
        end do
        !$omp end parallel do
      else
        !$omp parallel do reduction(+:resnorm)
        do i = 1, isize
          resnorm = resnorm + abs(Qx(i)-Qy(i))
        end do
        !$omp end parallel do
      end if

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      if (present(Qw)) then
        !$omp parallel do reduction(+:resnorm)
        do i = 1, isize
          resnorm = resnorm + Qw(i)*(Qx(i)-Qy(i))*(Qx(i)-Qy(i))
        end do
        !$omp end parallel do
      else
        !$omp parallel do reduction(+:resnorm)
        do i = 1, isize
          resnorm = resnorm + (Qx(i)-Qy(i))*(Qx(i)-Qy(i))
        end do
        !$omp end parallel do
      end if
      resnorm = sqrt(resnorm)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vektor (1111...) to has norm = 1.
      !$omp parallel do reduction(+:resnorm)
      do i = 1, isize
        resnorm = resnorm + abs(Qx(i)-Qy(i))
      end do
      !$omp end parallel do
      resnorm = resnorm / real(isize,DP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vektor (1111...) to has norm = 1.
      !$omp parallel do reduction(+:resnorm)
      do i = 1, isize
        resnorm = resnorm + (Qx(i)-Qy(i))*(Qx(i)-Qy(i))
      end do
      !$omp end parallel do
      resnorm = sqrt(resnorm / real(isize,DP))
      
    case (LINALG_NORMMAX)
      ! MAX-norm. Find the absolute largest entry.
      resnorm = abs(Qx(1)-Qy(1))
      j=1
      do i=2,isize
        qtemp = abs(Qx(i)-Qy(i))
        if (qtemp .gt. resnorm) then
          j = i
          resnorm = qtemp
        end if
      end do
      if (present(iposMax)) iposMax = j
      
    case default
      resnorm = -1.0_QP ! Unknown norm.
    end select
    
  end function

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorSortSngl (Fx, Fd, Itr)
  
!<description>
  ! Resorts the entries in the vector Fx corresponding to Itr.
  ! In particular, the first SIZE(Itr) entries of Dx are written resortedly
  ! to Dd.
!</description>
  
!<input>
  ! Array with permutation of 1..neq
  integer, dimension(:), intent(IN) :: Itr

  ! Source vector to be sorted
  real(SP), dimension(:), intent(IN) :: Fx
!</input>
  
!<output>
  ! The resorted vector
  real(SP), dimension(:), intent(OUT) :: Fd
!</output>
    
!</subroutine>
    
  ! local variable
  integer :: i
  
    do i = 1, size(Itr)
      Fd(i) = Fx(Itr(i))
    end do
  
  end subroutine 

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorSortDble (Dx, Dd, Itr)
  
!<description>
  ! Resorts the entries in the vector Dx corresponding to Itr.
  ! In particular, the first SIZE(Itr) entries of Dx are written resortedly
  ! to Dd.
!</description>
  
!<input>
  ! Source vector to be sorted
  real(DP), dimension(:), intent(IN) :: Dx
  
  ! Array with permutation of 1..neq.
  ! Itr(i) defines the number of the entry in Dx that should
  ! move to position i.
  integer, dimension(:), intent(IN) :: Itr
!</input>
  
!<output>
  ! The resorted vector
  real(DP), dimension(:), intent(OUT) :: Dd
!</output>
    
!</subroutine>
    
  ! local variable
  integer :: i
  
    do i = 1, size(Itr)
      Dd(i) = Dx(Itr(i))
    end do
  
  end subroutine 

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorSortQuad (Qx, Qd, Itr)
  
!<description>
  ! Resorts the entries in the vector Qx corresponding to Itr.
  ! In particular, the first SIZE(Itr) entries of Qx are written resortedly
  ! to Qd.
!</description>
  
!<input>
  ! Source vector to be sorted
  real(QP), dimension(:), intent(IN) :: Qx
  
  ! Array with permutation of 1..neq.
  ! Itr(i) defines the number of the entry in Dx that should
  ! move to position i.
  integer, dimension(:), intent(IN) :: Itr
!</input>
  
!<output>
  ! The resorted vector
  real(QP), dimension(:), intent(OUT) :: Qd
!</output>
    
!</subroutine>
    
  ! local variable
  integer :: i
  
    do i = 1, size(Itr)
      Qd(i) = Qx(Itr(i))
    end do
  
  end subroutine 

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorSortI32 (Ix, Id, Itr)
  
!<description>
  ! Resorts the entries in the vector Ix corresponding to Itr.
  ! In particular, the first SIZE(Itr) entries of Dx are written resortedly
  ! to Dd.
!</description>
  
!<input>
  ! Array with permutation of 1..neq
  integer, dimension(:), intent(IN) :: Itr

  ! Source vector to be sorted
  integer(I32), dimension(:), intent(IN) :: Ix
!</input>
  
!<output>
  ! The resorted vector
  integer(I32), dimension(:), intent(OUT) :: Id
!</output>
    
!</subroutine>
    
  ! local variable
  integer :: i
  
    do i = 1, size(Itr)
      Id(i) = Ix(Itr(i))
    end do
  
  end subroutine 

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorSortI64 (Ix, Id, Itr)
  
!<description>
  ! Resorts the entries in the vector Ix corresponding to Itr.
  ! In particular, the first SIZE(Itr) entries of Dx are written resortedly
  ! to Dd.
!</description>
  
!<input>
  ! Array with permutation of 1..neq
  integer, dimension(:), intent(IN) :: Itr

  ! Source vector to be sorted
  integer(I64), dimension(:), intent(IN) :: Ix
!</input>
  
!<output>
  ! The resorted vector
  integer(I64), dimension(:), intent(OUT) :: Id
!</output>
    
!</subroutine>
    
  ! local variable
  integer :: i
  
    do i = 1, size(Itr)
      Id(i) = Ix(Itr(i))
    end do
  
  end subroutine 

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_tensorProductSngl(Sx,Sy,Stensor)

!<description>
  ! Calculates the tensor product of two single precision vectors:
  ! Stensor = Sx (*) Sy
!</description>

!<input>
  ! First source vector
  real(SP), dimension(:), intent(IN) :: Sx

  ! Second source vector
  real(SP), dimension(:), intent(IN) :: Sy
!</input>

!<output>
  ! Tensor product
  real(SP), dimension(:,:), intent(OUT) :: Stensor
!</output>
!</subroutine>
    
  ! local variables
  integer :: i,j

    !$omp parallel do private(j)
    do i = 1, size(Sy)
      do j = 1, size(Sx)
        Stensor(j,i) = Sx(j)*Sy(i)
      end do
    end do
    !$omp end parallel do
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_tensorProductDble(Dx,Dy,Dtensor)

!<description>
  ! Calculates the tensor product of two double precision vectors:
  ! Dtensor = Dx (*) Dy
!</description>

!<input>
  ! First source vector
  real(DP), dimension(:), intent(IN) :: Dx

  ! Second source vector
  real(DP), dimension(:), intent(IN) :: Dy
!</input>

!<output>
  ! Tensor product
  real(DP), dimension(:,:), intent(OUT) :: Dtensor
!</output>
!</subroutine>
    
  ! local variables
  integer :: i,j

    !$omp parallel do private(j)
    do i = 1, size(Dy)
      do j = 1, size(Dx)
        Dtensor(j,i) = Dx(j)*Dy(i)
      end do
    end do
    !$omp end parallel do
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_tensorProductQuad(Qx,Qy,Qtensor)

!<description>
  ! Calculates the tensor product of two quad precision vectors:
  ! Qtensor = Qx (*) Qy
!</description>

!<input>
  ! First source vector
  real(QP), dimension(:), intent(IN) :: Qx

  ! Second source vector
  real(QP), dimension(:), intent(IN) :: Qy
!</input>

!<output>
  ! Tensor product
  real(QP), dimension(:,:), intent(OUT) :: Qtensor
!</output>
!</subroutine>
    
  ! local variables
  integer :: i,j

    !$omp parallel do private(j)
    do i = 1, size(Qy)
      do j = 1, size(Qx)
        Qtensor(j,i) = Qx(j)*Qy(i)
      end do
    end do
    !$omp end parallel do
      
  end subroutine
  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarSngl (Fx,fvalue,n)
  
!<description>
  ! This routine adds the value dvalue to each entry of the vector Fx.
!</description>
  
!<input>
  ! The value to add to every entry.
  real(SP) :: fvalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  real(SP), dimension(:), intent(INOUT) :: Fx
!</inputoutput>

!</subroutine>
    
  integer :: i
  
    ! Nothing to do?
    if(fvalue .eq. 0.0_SP) return
    
    if (.not. present(n)) then
      
      !$omp parallel do
      do i=1,size(Fx)
        Fx(i) = Fx(i) + fvalue
      end do
      !$omp end parallel do
      
    else
      
      !$omp parallel do
      do i=1,n
        Fx(i) = Fx(i) + fvalue
      end do
      !$omp end parallel do
      
    end if
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarDble (Dx,dvalue,n)
  
!<description>
  ! This routine adds the value dvalue to each entry of the vector Id.
!</description>
  
!<input>
  ! The value to add to every entry.
  real(DP) :: dvalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  real(DP), dimension(:), intent(INOUT) :: Dx
!</inputoutput>

!</subroutine>
    
  integer :: i
  
    ! Nothing to do?
    if(dvalue .eq. 0.0_DP) return
    
    if (.not. present(n)) then
    
      !$omp parallel do
      do i=1,size(Dx)
        Dx(i) = Dx(i) + dvalue
      end do
      !$omp parallel do
      
    else
    
      !$omp parallel do
      do i=1,n
        Dx(i) = Dx(i) + dvalue
      end do
      !$omp parallel do
      
    end if
    
  end subroutine
  

  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarQuad (Qx,qvalue,n)
  
!<description>
  ! This routine adds the value qvalue to each entry of the vector Id.
!</description>
  
!<input>
  ! The value to add to every entry.
  real(QP) :: qvalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  real(QP), dimension(:), intent(INOUT) :: Qx
!</inputoutput>

!</subroutine>
    
  integer :: i
  
    ! Nothing to do?
    if(qvalue .eq. 0.0_DP) return
    
    if (.not. present(n)) then
    
      !$omp parallel do
      do i=1,size(Qx)
        Qx(i) = Qx(i) + qvalue
      end do
      !$omp parallel do
      
    else
    
      !$omp parallel do
      do i=1,n
        Qx(i) = Qx(i) + qvalue
      end do
      !$omp parallel do
      
    end if
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarSngl2D (Fx,fvalue)
  
!<description>
  ! This routine adds the value dvalue to each entry of the vector Fx.
!</description>
  
!<input>
  ! The value to add to every entry.
  real(SP) :: fvalue
!</input>

!<inputoutput>
  ! Data array to be modified.
  real(SP), dimension(:,:), intent(INOUT) :: Fx
!</inputoutput>

!</subroutine>
    
  integer :: i,j
  
    ! Nothing to do?
    if(fvalue .eq. 0.0_SP) return
    
    !$omp parallel do private(i)
    do j=1,size(Fx,2)
      do i=1,size(Fx,1)
        Fx(i,j) = Fx(i,j) + fvalue
      end do
    end do
    !$omp end parallel do
    
  end subroutine
  
  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarDble2D (Dx,dvalue)
  
!<description>
  ! This routine adds the value dvalue to each entry of the vector Dx.
!</description>
  
!<input>
  ! The value to add to every entry.
  real(DP) :: dvalue
!</input>

!<inputoutput>
  ! Data array to be modified.
  real(DP), dimension(:,:), intent(INOUT) :: Dx
!</inputoutput>

!</subroutine>
    
  integer :: i,j
  
    ! Nothing to do?
    if(dvalue .eq. 0.0_DP) return
    
    !$omp parallel do private(i)
    do j = 1, size(Dx,2)
      do i = 1, size(Dx,1)
        Dx(i,j) = Dx(i,j) + dvalue
      end do
    end do
    !$omp end parallel do
    
  end subroutine
  
  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarQuad2D (Qx,qvalue)
  
!<description>
  ! This routine adds the value qvalue to each entry of the vector Qx.
!</description>
  
!<input>
  ! The value to add to every entry.
  real(QP) :: qvalue
!</input>

!<inputoutput>
  ! Data array to be modified.
  real(QP), dimension(:,:), intent(INOUT) :: Qx
!</inputoutput>

!</subroutine>
    
  integer :: i,j
  
    ! Nothing to do?
    if(qvalue .eq. 0.0_QP) return
    
    !$omp parallel do private(i)
    do j = 1, size(Qx,2)
      do i = 1, size(Qx,1)
        Qx(i,j) = Qx(i,j) + qvalue
      end do
    end do
    !$omp end parallel do
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarI8 (Ix,ivalue,n)
  
!<description>
  ! This routine adds the value ivalue to each entry of the vector Ix.
!</description>
  
!<input>
  ! The value to add to every entry.
  integer(I8) :: ivalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  integer(I8), dimension(:), intent(INOUT) :: Ix
!</inputoutput>

!</subroutine>
    
    integer :: i
    
    if (.not. present(n)) then
      !$omp parallel do
      do i = 1, size(Ix)
        Ix(i) = Ix(i) + ivalue
      end do
      !$omp end parallel do
    else
      !$omp parallel do
      do i = 1, n
        Ix(i) = Ix(i) + ivalue
      end do
      !$omp end parallel do
    end if
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarI16 (Ix,ivalue,n)
  
!<description>
  ! This routine adds the value ivalue to each entry of the vector Ix.
!</description>
  
!<input>
  ! The value to add to every entry.
  integer(I16) :: ivalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  integer(I16), dimension(:), intent(INOUT) :: Ix
!</inputoutput>

!</subroutine>
    
    integer :: i
    
    if (.not. present(n)) then
      !$omp parallel do
      do i = 1, size(Ix)
        Ix(i) = Ix(i) + ivalue
      end do
      !$omp end parallel do
    else
      !$omp parallel do
      do i = 1, n
        Ix(i) = Ix(i) + ivalue
      end do
      !$omp end parallel do
    end if
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarI32 (Ix,ivalue,n)
  
!<description>
  ! This routine adds the value ivalue to each entry of the vector Ix.
!</description>
  
!<input>
  ! The value to add to every entry.
  integer(I32) :: ivalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  integer(I32), dimension(:), intent(INOUT) :: Ix
!</inputoutput>

!</subroutine>
    
    integer :: i
    
    if (.not. present(n)) then
      !$omp parallel do
      do i = 1, size(Ix)
        Ix(i) = Ix(i) + ivalue
      end do
      !$omp end parallel do
    else
      !$omp parallel do
      do i = 1, n
        Ix(i) = Ix(i) + ivalue
      end do
      !$omp end parallel do
    end if
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarI64 (Ix,ivalue,n)
  
!<description>
  ! This routine adds the value ivalue to each entry of the vector Ix.
!</description>
  
!<input>
  ! The value to add to every entry.
  integer(I64) :: ivalue

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  integer(I64), dimension(:), intent(INOUT) :: Ix
!</inputoutput>

!</subroutine>
    
    integer :: i
    
    if (.not. present(n)) then
      !$omp parallel do
      do i = 1, size(Ix)
        Ix(i) = Ix(i) + ivalue
      end do
      !$omp end parallel do
    else
      !$omp parallel do
      do i = 1, n
        Ix(i) = Ix(i) + ivalue
      end do
      !$omp end parallel do
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultSngl (Fx,Fy,sc,n)
  
!<description>
  ! Performs componentwise multiplication: Fy = sc * Fx * Fy
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:), intent(IN) :: Fx
  
  ! Scaling factor
  real(SP), intent(IN)               :: sc

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  real(SP), dimension(:), intent(INOUT) :: Fy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      !$omp parallel do
      do i = 1, size(Fy,1)
        Fy(i) = sc*Fx(i)*Fy(i)
      end do
      !$omp end parallel do

    else

      !$omp parallel do
      do i = 1, n
        Fy(i) = sc*Fx(i)*Fy(i)
      end do
      !$omp end parallel do

    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultDble (Dx,Dy,dc,n)
  
!<description>
  ! Performs componentwise multiplication: Dy = dc * Dx * Dy
!</description>

!<input>
  
  ! First source vector
  real(DP), dimension(:), intent(IN) :: Dx
  
  ! Scaling factor
  real(DP), intent(IN)               :: dc

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  real(DP), dimension(:), intent(INOUT) :: Dy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      !$omp parallel do
      do i = 1, size(Dy,1)
        Dy(i) = dc*Dx(i)*Dy(i)
      end do
      !$omp end parallel do

    else

      !$omp parallel do
      do i = 1, n
        Dy(i) = dc*Dx(i)*Dy(i)
      end do
      !$omp end parallel do

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultQuad (Qx,Qy,qc,n)
  
!<description>
  ! Performs componentwise multiplication: Qy = qc * Qx * Qy
!</description>

!<input>
  
  ! First source vector
  real(QP), dimension(:), intent(IN) :: Qx
  
  ! Scaling factor
  real(QP), intent(IN)               :: qc

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:), intent(INOUT) :: Qy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      !$omp parallel do
      do i = 1, size(Qy,1)
        Qy(i) = qc*Qx(i)*Qy(i)
      end do
      !$omp end parallel do

    else

      !$omp parallel do
      do i = 1, n
        Qy(i) = qc*Qx(i)*Qy(i)
      end do
      !$omp end parallel do

    end if

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultDbleSngl (Fx,Dy,dc,n)
  
!<description>
  ! Performs componentwise multiplication: Dy = dc * Fx * Dy
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:), intent(IN) :: Fx
  
  ! Scaling factor
  real(DP), intent(IN)               :: dc

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  real(DP), dimension(:), intent(INOUT) :: Dy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      !$omp parallel do
      do i = 1, size(Dy,1)
        Dy(i) = dc*Fx(i)*Dy(i)
      end do
      !$omp end parallel do

    else

      !$omp parallel do
      do i = 1, n
        Dy(i) = dc*Fx(i)*Dy(i)
      end do
      !$omp end parallel do

    end if

  end subroutine
  
! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultQuadSngl (Fx,Qy,Qc,n)
  
!<description>
  ! Performs componentwise multiplication: Qy = qc * Fx * Qy
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:), intent(IN) :: Fx
  
  ! Scaling factor
  real(QP), intent(IN)               :: qc

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:), intent(INOUT) :: Qy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      !$omp parallel do
      do i = 1, size(Qy,1)
        Qy(i) = qc*real(Fx(i),QP)*Qy(i)
      end do
      !$omp end parallel do

    else

      !$omp parallel do
      do i = 1, n
        Qy(i) = qc*real(Fx(i),QP)*Qy(i)
      end do
      !$omp end parallel do

    end if

  end subroutine
 
! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultQuadDble (Dx,Qy,Qc,n)
  
!<description>
  ! Performs componentwise multiplication: Qy = qc * Dx * Qy
!</description>

!<input>
  
  ! First source vector
  real(DP), dimension(:), intent(IN) :: Dx
  
  ! Scaling factor
  real(QP), intent(IN)               :: qc

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:), intent(INOUT) :: Qy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      !$omp parallel do
      do i = 1, size(Qy,1)
        Qy(i) = qc*real(Dx(i),QP)*Qy(i)
      end do
      !$omp end parallel do

    else

      !$omp parallel do
      do i = 1, n
        Qy(i) = qc*real(Dx(i),QP)*Qy(i)
      end do
      !$omp end parallel do

    end if

  end subroutine
 
end module

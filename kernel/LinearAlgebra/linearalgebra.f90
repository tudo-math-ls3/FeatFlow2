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
  
  private
  
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
  
  public :: lalg_copyVectorInt
  
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
  
  public :: lalg_copyVectorInt2D

    interface lalg_copyVectorInt3D
    module procedure lalg_copyVectorI8_3D
    module procedure lalg_copyVectorI8I16_3D
    module procedure lalg_copyVectorI8I32_3D
    module procedure lalg_copyVectorI8I64_3D
    module procedure lalg_copyVectorI16_3D
    module procedure lalg_copyVectorI16I8_3D
    module procedure lalg_copyVectorI16I32_3D
    module procedure lalg_copyVectorI16I64_3D
    module procedure lalg_copyVectorI32_3D
    module procedure lalg_copyVectorI32I8_3D
    module procedure lalg_copyVectorI32I16_3D
    module procedure lalg_copyVectorI32I64_3D
    module procedure lalg_copyVectorI64_3D
    module procedure lalg_copyVectorI64I8_3D
    module procedure lalg_copyVectorI64I16_3D
    module procedure lalg_copyVectorI64I32_3D
  end interface  
  
  public :: lalg_copyVectorInt3D

  interface lalg_copyVector
    module procedure lalg_copyVectorSngl
    module procedure lalg_copyVectorDble
    module procedure lalg_copyVectorSnglDbl
    module procedure lalg_copyVectorDblSngl
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
    module procedure lalg_copyVectorSnglDbl2D
    module procedure lalg_copyVectorDblSngl2D
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
    module procedure lalg_copyVectorSngl3D
    module procedure lalg_copyVectorDble3D
    module procedure lalg_copyVectorSnglDbl3D
    module procedure lalg_copyVectorDblSngl3D
    module procedure lalg_copyVectorI8_3D
    module procedure lalg_copyVectorI8I16_3D
    module procedure lalg_copyVectorI8I32_3D
    module procedure lalg_copyVectorI8I64_3D
    module procedure lalg_copyVectorI16_3D
    module procedure lalg_copyVectorI16I8_3D
    module procedure lalg_copyVectorI16I32_3D
    module procedure lalg_copyVectorI16I64_3D
    module procedure lalg_copyVectorI32_3D
    module procedure lalg_copyVectorI32I8_3D
    module procedure lalg_copyVectorI32I16_3D
    module procedure lalg_copyVectorI32I64_3D
    module procedure lalg_copyVectorI64_3D
    module procedure lalg_copyVectorI64I8_3D
    module procedure lalg_copyVectorI64I16_3D
    module procedure lalg_copyVectorI64I32_3D
    module procedure lalg_copyVectorLogical3D
    module procedure lalg_copyVectorChar3D
#ifdef ENABLE_QUADPREC
    module procedure lalg_copyVectorQuad
    module procedure lalg_copyVectorSnglQuad
    module procedure lalg_copyVectorDblQuad
    module procedure lalg_copyVectorQuadSngl
    module procedure lalg_copyVectorQuadDbl
    module procedure lalg_copyVectorQuad2D
    module procedure lalg_copyVectorSnglQuad2D
    module procedure lalg_copyVectorDblQuad2D
    module procedure lalg_copyVectorQuadSngl2D
    module procedure lalg_copyVectorQuadDbl2D
    module procedure lalg_copyVectorQuad3D
    module procedure lalg_copyVectorSnglQuad3D
    module procedure lalg_copyVectorDblQuad3D
    module procedure lalg_copyVectorQuadSngl3D
    module procedure lalg_copyVectorQuadDbl3D
#endif
  end interface
  
  public :: lalg_copyVector

  public :: lalg_copyVectorSngl
  public :: lalg_copyVectorDble
  public :: lalg_copyVectorQuad
  public :: lalg_copyVectorSnglDbl
  public :: lalg_copyVectorSnglQuad
  public :: lalg_copyVectorDblSngl
  public :: lalg_copyVectorDblQuad
  public :: lalg_copyVectorQuadSngl
  public :: lalg_copyVectorQuadDbl
  public :: lalg_copyVectorLogical
  public :: lalg_copyVectorChar

  public :: lalg_copyVectorSngl2D
  public :: lalg_copyVectorDble2D
  public :: lalg_copyVectorQuad2D
  public :: lalg_copyVectorSnglDbl2D
  public :: lalg_copyVectorSnglQuad2D
  public :: lalg_copyVectorDblSngl2D
  public :: lalg_copyVectorDblQuad2D
  public :: lalg_copyVectorQuadSngl2D
  public :: lalg_copyVectorQuadDbl2D
  public :: lalg_copyVectorLogical2D
  public :: lalg_copyVectorChar2D

  public :: lalg_copyVectorSngl3D
  public :: lalg_copyVectorDble3D
  public :: lalg_copyVectorQuad3D
  public :: lalg_copyVectorSnglDbl3D
  public :: lalg_copyVectorSnglQuad3D
  public :: lalg_copyVectorDblSngl3D
  public :: lalg_copyVectorDblQuad3D
  public :: lalg_copyVectorQuadSngl3D
  public :: lalg_copyVectorQuadDbl3D
  public :: lalg_copyVectorLogical3D
  public :: lalg_copyVectorChar3D

  interface lalg_scaleVector
    module procedure lalg_scaleVectorSngl
    module procedure lalg_scaleVectorDble
    module procedure lalg_scaleVectorSngl2D
    module procedure lalg_scaleVectorDble2D
    module procedure lalg_scaleVectorSngl3D
    module procedure lalg_scaleVectorDble3D
#ifdef ENABLE_QUADPREC
    module procedure lalg_scaleVectorQuad
    module procedure lalg_scaleVectorQuad2D
    module procedure lalg_scaleVectorQuad3D
#endif
  end interface
  
  public :: lalg_scaleVector
  public :: lalg_scaleVectorSngl
  public :: lalg_scaleVectorDble
  public :: lalg_scaleVectorQuad
  public :: lalg_scaleVectorSngl2D
  public :: lalg_scaleVectorDble2D
  public :: lalg_scaleVectorQuad2D
  public :: lalg_scaleVectorSngl3D
  public :: lalg_scaleVectorDble3D
  public :: lalg_scaleVectorQuad3D

  interface lalg_clearVectorInt
    module procedure lalg_clearVectorI8
    module procedure lalg_clearVectorI16
    module procedure lalg_clearVectorI32
    module procedure lalg_clearVectorI64
  end interface
  
  public :: lalg_clearVectorInt

  interface lalg_clearVectorInt2D
    module procedure lalg_clearVectorI8_2D
    module procedure lalg_clearVectorI16_2D
    module procedure lalg_clearVectorI32_2D
    module procedure lalg_clearVectorI64_2D
  end interface
  
  public :: lalg_clearVectorInt2D

  interface lalg_clearVectorInt3D
    module procedure lalg_clearVectorI8_3D
    module procedure lalg_clearVectorI16_3D
    module procedure lalg_clearVectorI32_3D
    module procedure lalg_clearVectorI64_3D
  end interface
  
  public :: lalg_clearVectorInt3D

  interface lalg_clearVector
    module procedure lalg_clearVectorSngl
    module procedure lalg_clearVectorDble
    module procedure lalg_clearVectorI8
    module procedure lalg_clearVectorI16
    module procedure lalg_clearVectorI32
    module procedure lalg_clearVectorI64
    module procedure lalg_clearVectorSngl2D
    module procedure lalg_clearVectorDble2D
    module procedure lalg_clearVectorI8_2D
    module procedure lalg_clearVectorI16_2D
    module procedure lalg_clearVectorI32_2D
    module procedure lalg_clearVectorI64_2D
    module procedure lalg_clearVectorSngl3D
    module procedure lalg_clearVectorDble3D
    module procedure lalg_clearVectorI8_3D
    module procedure lalg_clearVectorI16_3D
    module procedure lalg_clearVectorI32_3D
    module procedure lalg_clearVectorI64_3D
#ifdef ENABLE_QUADPREC
    module procedure lalg_clearVectorQuad
    module procedure lalg_clearVectorQuad2D
    module procedure lalg_clearVectorQuad3D
#endif
  end interface
  
  public :: lalg_clearVector
  public :: lalg_clearVectorSngl
  public :: lalg_clearVectorDble
  public :: lalg_clearVectorQuad
  public :: lalg_clearVectorSngl2D
  public :: lalg_clearVectorDble2D
  public :: lalg_clearVectorQuad2D
  public :: lalg_clearVectorSngl3D
  public :: lalg_clearVectorDble3D
  public :: lalg_clearVectorQuad3D

  interface lalg_setVectorInt
    module procedure lalg_setVectorI8
    module procedure lalg_setVectorI16
    module procedure lalg_setVectorI32
    module procedure lalg_setVectorI64
  end interface
  
  public :: lalg_setVectorInt

  interface lalg_setVectorInt2D
    module procedure lalg_setVectorI8_2D
    module procedure lalg_setVectorI16_2D
    module procedure lalg_setVectorI32_2D
    module procedure lalg_setVectorI64_2D
  end interface
  
  public :: lalg_setVectorInt2D

   interface lalg_setVectorInt3D
    module procedure lalg_setVectorI8_3D
    module procedure lalg_setVectorI16_3D
    module procedure lalg_setVectorI32_3D
    module procedure lalg_setVectorI64_3D
  end interface
  
  public :: lalg_setVectorInt3D

  interface lalg_setVector
    module procedure lalg_setVectorSngl
    module procedure lalg_setVectorDble
    module procedure lalg_setVectorI8
    module procedure lalg_setVectorI16
    module procedure lalg_setVectorI32
    module procedure lalg_setVectorI64
    module procedure lalg_setVectorLogical
    module procedure lalg_setVectorChar
    module procedure lalg_setVectorSngl2D
    module procedure lalg_setVectorDble2D
    module procedure lalg_setVectorI8_2D
    module procedure lalg_setVectorI16_2D
    module procedure lalg_setVectorI32_2D
    module procedure lalg_setVectorI64_2D
    module procedure lalg_setVectorLogical2D
    module procedure lalg_setVectorChar2D
    module procedure lalg_setVectorSngl3D
    module procedure lalg_setVectorDble3D
    module procedure lalg_setVectorI8_3D
    module procedure lalg_setVectorI16_3D
    module procedure lalg_setVectorI32_3D
    module procedure lalg_setVectorI64_3D
    module procedure lalg_setVectorLogical3D
    module procedure lalg_setVectorChar3D
#ifdef ENABLE_QUADPREC
    module procedure lalg_setVectorQuad
    module procedure lalg_setVectorQuad2D
    module procedure lalg_setVectorQuad3D
#endif
  end interface
  
  public :: lalg_setVector

  public :: lalg_setVectorQuad
  public :: lalg_setVectorQuad2D
  public :: lalg_setVectorQuad3D
  public :: lalg_setVectorSngl
  public :: lalg_setVectorSngl2D
  public :: lalg_setVectorSngl3D
  public :: lalg_setVectorDble
  public :: lalg_setVectorDble2D
  public :: lalg_setVectorDble3D
  public :: lalg_setVectorLogical
  public :: lalg_setVectorLogical2D
  public :: lalg_setVectorLogical3D
  public :: lalg_setVectorChar
  public :: lalg_setVectorChar2D
  public :: lalg_setVectorChar3D

  interface lalg_vectorLinearComb
    module procedure lalg_vectorLinearCombSngl
    module procedure lalg_vectorLinearCombDble
    module procedure lalg_vectorLinearCombSnglDble
    module procedure lalg_vectorLinearCombSngl2D
    module procedure lalg_vectorLinearCombDble2D
    module procedure lalg_vectorLinearCombSnglDble2D
#ifdef ENABLE_QUADPREC
    module procedure lalg_vectorLinearCombQuad
    module procedure lalg_vectorLinearCombSnglQuad
    module procedure lalg_vectorLinearCombDblQuad
    module procedure lalg_vectorLinearCombQuad2D
    module procedure lalg_vectorLinearCombSnglQuad2D
    module procedure lalg_vectorLinearCombDblQuad2D
#endif
  end interface
  
  public :: lalg_vectorLinearComb
  public :: lalg_vectorLinearCombSngl
  public :: lalg_vectorLinearCombDble
  public :: lalg_vectorLinearCombQuad
  public :: lalg_vectorLinearCombSnglDble
  public :: lalg_vectorLinearCombSnglQuad
  public :: lalg_vectorLinearCombDblQuad
  public :: lalg_vectorLinearCombSngl2D
  public :: lalg_vectorLinearCombDble2D
  public :: lalg_vectorLinearCombQuad2D
  public :: lalg_vectorLinearCombSnglDble2D
  public :: lalg_vectorLinearCombSnglQuad2D
  public :: lalg_vectorLinearCombDblQuad2D

  interface lalg_scalarProduct
    module procedure lalg_scalarProductSngl
    module procedure lalg_scalarProductDble
    module procedure lalg_scalarProductSngl2D
    module procedure lalg_scalarProductDble2D
#ifdef ENABLE_QUADPREC
    module procedure lalg_scalarProductQuad
    module procedure lalg_scalarProductQuad2D
#endif
  end interface
  
  public :: lalg_scalarProduct
  public :: lalg_scalarProductSngl
  public :: lalg_scalarProductDble
  public :: lalg_scalarProductQuad
  public :: lalg_scalarProductSngl2D
  public :: lalg_scalarProductDble2D
  public :: lalg_scalarProductQuad2D

  interface lalg_norm
    module procedure lalg_normSngl
    module procedure lalg_normDble
#ifdef ENABLE_QUADPREC
    module procedure lalg_normQuad
#endif
  end interface
  
  public :: lalg_norm,lalg_normSngl,lalg_normDble,lalg_normQuad

  interface lalg_errorNorm
    module procedure lalg_errorNormSngl
    module procedure lalg_errorNormDble
#ifdef ENABLE_QUADPREC
    module procedure lalg_errorNormQuad
#endif
  end interface
  
  public :: lalg_errorNorm
  public :: lalg_errorNormSngl
  public :: lalg_errorNormDble
  public :: lalg_errorNormQuad

  interface lalg_vectorSortInt    
    module procedure lalg_vectorSortI32
    module procedure lalg_vectorSortI64
  end interface
  
  public :: lalg_vectorSortInt

  interface lalg_vectorSort
    module procedure lalg_vectorSortSngl
    module procedure lalg_vectorSortDble
    module procedure lalg_vectorSortI32
    module procedure lalg_vectorSortI64
#ifdef ENABLE_QUADPREC
    module procedure lalg_vectorSortQuad
#endif
  end interface
  
  public :: lalg_vectorSort
  public :: lalg_vectorSortSngl
  public :: lalg_vectorSortDble
  public :: lalg_vectorSortQuad
  public :: lalg_vectorSortI32
  public :: lalg_vectorSortI64

  interface lalg_vectorAddScalarInt
    module procedure lalg_vectorAddScalarI8
    module procedure lalg_vectorAddScalarI16
    module procedure lalg_vectorAddScalarI32
    module procedure lalg_vectorAddScalarI64
  end interface
  
  public :: lalg_vectorAddScalarInt

  interface lalg_vectorAddScalar
    module procedure lalg_vectorAddScalarSngl
    module procedure lalg_vectorAddScalarDble
    module procedure lalg_vectorAddScalarSngl2D
    module procedure lalg_vectorAddScalarDble2D
    module procedure lalg_vectorAddScalarI8
    module procedure lalg_vectorAddScalarI16
    module procedure lalg_vectorAddScalarI32
    module procedure lalg_vectorAddScalarI64
#ifdef ENABLE_QUADPREC
    module procedure lalg_vectorAddScalarQuad
    module procedure lalg_vectorAddScalarQuad2D
#endif
  end interface
  
  public :: lalg_vectorAddScalar
  public :: lalg_vectorAddScalarSngl
  public :: lalg_vectorAddScalarDble
  public :: lalg_vectorAddScalarQuad
  public :: lalg_vectorAddScalarSngl2D
  public :: lalg_vectorAddScalarDble2D
  public :: lalg_vectorAddScalarQuad2D
  public :: lalg_vectorAddScalarI8
  public :: lalg_vectorAddScalarI16
  public :: lalg_vectorAddScalarI32
  public :: lalg_vectorAddScalarI64

  interface lalg_vectorCompMult
    module procedure lalg_vectorCompMultSngl
    module procedure lalg_vectorCompMultDble
    module procedure lalg_vectorCompMultDbleSngl
#ifdef ENABLE_QUADPREC
    module procedure lalg_vectorCompMultQuad
    module procedure lalg_vectorCompMultQuadSngl
    module procedure lalg_vectorCompMultQuadDble
#endif
  end interface

  public :: lalg_vectorCompMult
  public :: lalg_vectorCompMultSngl
  public :: lalg_vectorCompMultDble
  public :: lalg_vectorCompMultQuad
  public :: lalg_vectorCompMultDbleSngl
  public :: lalg_vectorCompMultQuadSngl
  public :: lalg_vectorCompMultQuadDble

!<constants>

!<constantblock description="Constants identifying vector norms">

  ! Sum of the absolute values of entries
  integer, parameter, public :: LINALG_NORMSUM    = -1

  ! Euclidian vector norm: (vector,vector)
  integer, parameter, public :: LINALG_NORMEUCLID = 0

  ! $l_1$-norm: 1/NEQ * sum(abs(entries))
  integer, parameter, public :: LINALG_NORML1     = 1
  
  ! $l_2$-norm: 1/sqrt(NEQ) * (vector,vector)
  integer, parameter, public :: LINALG_NORML2     = 2
  
  ! max-norm
  integer, parameter, public :: LINALG_NORMMAX    = 3
  
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
  real(SP), dimension(:), intent(in) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:), intent(out) :: Fy
  
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
  real(DP), dimension(:), intent(in) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:), intent(out) :: Dy
  
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
  real(QP), dimension(:), intent(in) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:), intent(out) :: Qy
  
!</output>
  
!</subroutine>

  integer :: i

    if (.not. present(n)) then
    
      do i = 1, size(Qx)
        Qy(i) = Qx(i)
      end do
      
      !call QCOPY(size(Qx),Qx,1,Qy,1)
      
    else

      do i = 1, n
        Qy(i) = Qx(i)
      end do
      
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
  real(SP), dimension(:), intent(in) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:), intent(out) :: Dy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Fx)
        Dy(i) = real(Fx(i),DP)
      end do
    
    else

      do i = 1, n
        Dy(i) = real(Fx(i),DP)
      end do
    
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
  real(SP), dimension(:), intent(in) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:), intent(out) :: Qy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Fx)
        Qy(i) = real(Fx(i),QP)
      end do
    
    else

      do i = 1, n
        Qy(i) = real(Fx(i),QP)
      end do
    
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
  real(DP), dimension(:), intent(in) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:), intent(out) :: Fy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then

      do i = 1, size(Dx)
        Fy(i) = real(Dx(i),SP)
      end do

    else
    
      do i = 1, n
        Fy(i) = real(Dx(i),SP)
      end do

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
  real(DP), dimension(:), intent(in) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:), intent(out) :: Qy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then

      do i = 1, size(Dx)
        Qy(i) = real(Dx(i),QP)
      end do

    else
    
      do i = 1, n
        Qy(i) = real(Dx(i),QP)
      end do

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
  real(QP), dimension(:), intent(in) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:), intent(out) :: Fy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then

      do i = 1, size(Qx)
        Fy(i) = real(Qx(i),SP)
      end do

    else
    
      do i = 1, n
        Fy(i) = real(Qx(i),SP)
      end do

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
  real(QP), dimension(:), intent(in) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:), intent(out) :: Dy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then

      do i = 1, size(Qx)
        Dy(i) = real(Qx(i),DP)
      end do

    else
    
      do i = 1, n
        Dy(i) = real(Qx(i),DP)
      end do

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
  integer(I8), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = Ix(i)
      end do
      
    else

      do i = 1, n
        Iy(i) = Ix(i)
      end do
      
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
  integer(I8), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I16)
      end do
      
    else

      do i = 1, n
        Iy(i) = int(Ix(i),I16)
      end do
      
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
  integer(I8), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I32)
      end do
      
    else

      do i = 1, n
        Iy(i) = int(Ix(i),I32)
      end do
      
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
  integer(I8), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I64)
      end do
      
    else

      do i = 1, n
        Iy(i) = int(Ix(i),I64)
      end do
      
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
  integer(I16), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = Ix(i)
      end do
      
    else

      do i = 1, n
        Iy(i) = Ix(i)
      end do
      
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
  integer(I16), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I8)
      end do
      
    else

      do i = 1, n
        Iy(i) = int(Ix(i),I8)
      end do
      
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
  integer(I16), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I32)
      end do
      
    else

      do i = 1, n
        Iy(i) = int(Ix(i),I32)
      end do
      
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
  integer(I16), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I64)
      end do
      
    else

      do i = 1, n
        Iy(i) = int(Ix(i),I64)
      end do
      
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
  integer(I32), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = Ix(i)
      end do
      
    else

      do i = 1, n
        Iy(i) = Ix(i)
      end do
      
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
  integer(I32), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I8)
      end do
      
    else

      do i = 1, n
        Iy(i) = int(Ix(i),I8)
      end do
      
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
  integer(I32), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I16)
      end do
      
    else

      do i = 1, n
        Iy(i) = int(Ix(i),I16)
      end do
      
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
  integer(I32), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I64)
      end do
      
    else

      do i = 1, n
        Iy(i) = int(Ix(i),I64)
      end do
      
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
  integer(I64), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = Ix(i)
      end do
      
    else

      do i = 1, n
        Iy(i) = Ix(i)
      end do
      
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
  integer(I64), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I8)
      end do
      
    else

      do i = 1, n
        Iy(i) = int(Ix(i),I8)
      end do
      
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
  integer(I64), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I16)
      end do
      
    else

      do i = 1, n
        Iy(i) = int(Ix(i),I16)
      end do
      
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
  integer(I64), dimension(:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I32)
      end do
      
    else

      do i = 1, n
        Iy(i) = int(Ix(i),I32)
      end do
      
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
  logical, dimension(:), intent(in) :: Lx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  logical, dimension(:), intent(out) :: Ly
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Lx)
        Ly(i) = Lx(i)
      end do
    
    else

      do i = 1, n
        Ly(i) = Lx(i)
      end do
    
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
  character, dimension(:), intent(in) :: Sx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  
  ! Destination vector
  character, dimension(:), intent(out) :: Sy
  
!</output>
  
!</subroutine>
  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Sx)
        Sy(i) = Sx(i)
      end do
    
    else

      do i = 1, n
        Sy(i) = Sx(i)
      end do
    
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
  real(SP), dimension(:,:), intent(in) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:,:), intent(out) :: Fy
  
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
  real(DP), dimension(:,:), intent(in) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m
  
!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:,:), intent(out) :: Dy
  
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
  real(QP), dimension(:,:), intent(in) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m
  
!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:,:), intent(out) :: Qy
  
!</output>
  
!</subroutine>

  integer :: i,j

    if (present(n) .and. present(m)) then
    
      do j = 1, m
        do i = 1, n
          Qy(i,j) = Qx(i,j)
        end do
      end do
      
      !call QCOPY(n*m,Qx,1,Qy,1)
      
    else

      do j = 1, size(Qx,2)
        do i = 1, size(Qx,1)
          Qy(i,j) = Qx(i,j)
        end do
      end do

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
  real(SP), dimension(:,:), intent(in) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:,:), intent(out) :: Dy
  
!</output>
  
!</subroutine>
  integer :: i,j
  
    if (present(n) .and. present(m)) then

      do j = 1, m
        do i = 1, n
          Dy(i,j) = real(Fx(i,j),DP)
        end do
      end do

    else

      do j = 1, size(Fx,2)
        do i = 1, size(Fx,1)
          Dy(i,j) = real(Fx(i,j),DP)
        end do
      end do
      
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
  real(SP), dimension(:,:), intent(in) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:,:), intent(out) :: Qy
  
!</output>
  
!</subroutine>
  integer :: i,j
  
    if (present(n) .and. present(m)) then

      do j = 1, m
        do i = 1, n
          Qy(i,j) = real(Fx(i,j),QP)
        end do
      end do

    else

      do j = 1, size(Fx,2)
        do i = 1, size(Fx,1)
          Qy(i,j) = real(Fx(i,j),QP)
        end do
      end do
      
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
  real(DP), dimension(:,:), intent(in) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:,:), intent(out) :: Fy
  
!</output>
  
!</subroutine>
  integer :: i,j
  
    if (present(n) .and. present(m)) then

      do j = 1, m
        do i = 1, n
          Fy(i,j) = real(Dx(i,j),SP)
        end do
      end do

    else

      do j = 1, size(Dx,2)
        do i = 1, size(Dx,1)
          Fy(i,j) = real(Dx(i,j),SP)
        end do
      end do
      
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
  real(DP), dimension(:,:), intent(in) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:,:), intent(out) :: Qy
  
!</output>
  
!</subroutine>
  integer :: i,j
  
    if (present(n) .and. present(m)) then

      do j = 1, m
        do i = 1, n
          Qy(i,j) = real(Dx(i,j),QP)
        end do
      end do

    else

      do j = 1, size(Dx,2)
        do i = 1, size(Dx,1)
          Qy(i,j) = real(Dx(i,j),QP)
        end do
      end do
      
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
  real(QP), dimension(:,:), intent(in) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:,:), intent(out) :: Fy
  
!</output>
  
!</subroutine>
  integer :: i,j
  
    if (present(n) .and. present(m)) then

      do j = 1, m
        do i = 1, n
          Fy(i,j) = real(Qx(i,j),SP)
        end do
      end do

    else

      do j = 1, size(Qx,2)
        do i = 1, size(Qx,1)
          Fy(i,j) = real(Qx(i,j),SP)
        end do
      end do
      
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
  real(QP), dimension(:,:), intent(in) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:,:), intent(out) :: Dy
  
!</output>
  
!</subroutine>
  integer :: i,j
  
    if (present(n) .and. present(m)) then

      do j = 1, m
        do i = 1, n
          Dy(i,j) = real(Qx(i,j),DP)
        end do
      end do

    else

      do j = 1, size(Qx,2)
        do i = 1, size(Qx,1)
          Dy(i,j) = real(Qx(i,j),DP)
        end do
      end do
      
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
  integer(I8), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = Ix(i,j)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = Ix(i,j)
        end do
      end do

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
  integer(I8), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I16)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I16)
        end do
      end do

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
  integer(I8), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I32)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I32)
        end do
      end do

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
  integer(I8), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I64)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I64)
        end do
      end do

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
  integer(I16), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = Ix(i,j)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = Ix(i,j)
        end do
      end do

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
  integer(I16), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I8)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I8)
        end do
      end do

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
  integer(I16), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I32)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I32)
        end do
      end do

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
  integer(I16), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I64)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I64)
        end do
      end do

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
  integer(I32), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = Ix(i,j)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = Ix(i,j)
        end do
      end do

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
  integer(I32), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I8)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I8)
        end do
      end do

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
  integer(I32), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I16)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I16)
        end do
      end do

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
  integer(I32), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I64)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I64)
        end do
      end do

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
  integer(I64), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = Ix(i,j)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = Ix(i,j)
        end do
      end do

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
  integer(I64), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I8)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I8)
        end do
      end do

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
  integer(I64), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I16)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I16)
        end do
      end do

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
  integer(I64), dimension(:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if(present(m) .and. present(n)) then
    
      do j = 1, m
        do i = 1, n
          Iy(i,j) = int(Ix(i,j),I32)
        end do
      end do
    
    else
    
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Iy(i,j) = int(Ix(i,j),I32)
        end do
      end do

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
  logical, dimension(:,:), intent(in) :: Lx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  logical, dimension(:,:), intent(out) :: Ly
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if (present(n) .and. present(m)) then

      do j = 1, m
        do i = 1, n
          Ly(i,j) = Lx(i,j)
        end do
      end do

    else
      
      do j = 1, size(Lx,2)
        do i = 1, size(Lx,1)
          Ly(i,j) = Lx(i,j)
        end do
      end do

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
  character, dimension(:,:), intent(in) :: Sx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m

!</input>

!<output>
  
  ! Destination vector
  character, dimension(:,:), intent(out) :: Sy
  
!</output>
  
!</subroutine>

  integer :: i,j
  
    if (present(n) .and. present(m)) then

      do j = 1, m
        do i = 1, n
          Sy(i,j) = Sx(i,j)
        end do
      end do

    else

      do j = 1, size(Sx,2)
        do i = 1, size(Sx,1)
          Sy(i,j) = Sx(i,j)
        end do
      end do

    end if
  
  end subroutine

   ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSngl3D (Fx,Fy,n,m,o)
  
!<description>
  ! Copies a single precision vector: Fy = Fx
!</description>

!<input>
  
  ! Source vector
  real(SP), dimension(:,:,:), intent(in) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:,:,:), intent(out) :: Fy
  
!</output>
  
!</subroutine>

    if (present(n) .and. present(m) .and. present(o)) then
      call SCOPY(n*m*o,Fx,1,Fy,1)
    else
      call SCOPY(size(Fx,1)*size(Fx,2)*size(Fx,3),Fx,1,Fy,1)
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDble3D (Dx,Dy,n,m,o)
  
!<description>
  ! Copies a double precision vector: Dy = Dx
!</description>

!<input>
  
  ! Source vector
  real(DP), dimension(:,:,:), intent(in) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o
  
!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:,:,:), intent(out) :: Dy
  
!</output>
  
!</subroutine>

    if (present(n) .and. present(m) .and. present(o)) then
      call DCOPY(n*m*o,Dx,1,Dy,1)
    else
      call DCOPY(size(Dx,1)*size(Dx,2)*size(Dx,3),Dx,1,Dy,1)
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQuad3D (Qx,Qy,n,m,o)
  
!<description>
  ! Copies a quad precision vector: Qy = Qx
!</description>

!<input>
  
  ! Source vector
  real(QP), dimension(:,:,:), intent(in) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o
  
!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:,:,:), intent(out) :: Qy
  
!</output>
  
!</subroutine>

  integer :: i,j,k

    if (present(n) .and. present(m) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Qy(i,j,k) = Qx(i,j,k)
          end do
        end do
      end do
      
      !call QCOPY(n*m*o,Qx,1,Qy,1)
      
    else

      do k = 1, size(Qx,3)
        do j = 1, size(Qx,2)
          do i = 1, size(Qx,1)
            Qy(i,j,k) = Qx(i,j,k)
          end do
        end do
      end do

      !call QCOPY(size(Qx,1)*size(Qx,2)*size(Qx,3),Qx,1,Qy,1)
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSnglDbl3D (Fx,Dy,n,m,o)
  
!<description>
  ! Copies single precision vector to double precision vector: Dy = Fx
!</description>

!<input>
  
  ! Source vector
  real(SP), dimension(:,:,:), intent(in) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:,:,:), intent(out) :: Dy
  
!</output>
  
!</subroutine>
  integer :: i,j,k
  
    if (present(n) .and. present(m) .and. present(o)) then

      do k = 1, o
        do j = 1, m
          do i = 1, n
            Dy(i,j,k) = real(Fx(i,j,k),DP)
          end do
        end do
      end do

    else

      do k = 1, size(Fx,3)
        do j = 1, size(Fx,2)
          do i = 1, size(Fx,1)
            Dy(i,j,k) = real(Fx(i,j,k),DP)
          end do
        end do
      end do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSnglQuad3D (Fx,Qy,n,m,o)
  
!<description>
  ! Copies single precision vector to quad precision vector: Qy = Fx
!</description>

!<input>
  
  ! Source vector
  real(SP), dimension(:,:,:), intent(in) :: Fx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:,:,:), intent(out) :: Qy
  
!</output>
  
!</subroutine>
  integer :: i,j,k
  
    if (present(n) .and. present(m) .and. present(o)) then

      do k = 1, o
        do j = 1, m
          do i = 1, n
            Qy(i,j,k) = real(Fx(i,j,k),QP)
          end do
        end do
      end do

    else

      do k = 1, size(Fx,3)
        do j = 1, size(Fx,2)
          do i = 1, size(Fx,1)
            Qy(i,j,k) = real(Fx(i,j,k),QP)
          end do
        end do
      end do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDblSngl3D (Dx,Fy,n,m,o)
  
!<description>
  ! Copies double precision vector to single precision vector: Fy = Dx
!</description>

!<input>
  
  ! Source vector
  real(DP), dimension(:,:,:), intent(in) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:,:,:), intent(out) :: Fy
  
!</output>
  
!</subroutine>
  integer :: i,j,k
  
    if (present(n) .and. present(m) .and. present(o)) then

      do k = 1, o
        do j = 1, m
          do i = 1, n
            Fy(i,j,k) = real(Dx(i,j,k),SP)
          end do
        end do
      end do

    else

      do k = 1, size(Dx,3)
        do j = 1, size(Dx,2)
          do i = 1, size(Dx,1)
            Fy(i,j,k) = real(Dx(i,j,k),SP)
          end do
        end do
      end do
      
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDblQuad3D (Dx,Qy,n,m,o)
  
!<description>
  ! Copies double precision vector to quad precision vector: Qy = Dx
!</description>

!<input>
  
  ! Source vector
  real(DP), dimension(:,:,:), intent(in) :: Dx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  real(QP), dimension(:,:,:), intent(out) :: Qy
  
!</output>
  
!</subroutine>
  integer :: i,j,k
  
    if (present(n) .and. present(m) .and. present(o)) then

      do k = 1, o
        do j = 1, m
          do i = 1, n
            Qy(i,j,k) = real(Dx(i,j,k),QP)
          end do
        end do
      end do

    else

      do k = 1, size(Dx,3)
        do j = 1, size(Dx,2)
          do i = 1, size(Dx,1)
            Qy(i,j,k) = real(Dx(i,j,k),QP)
          end do
        end do
      end do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQuadSngl3D (Qx,Fy,n,m,o)
  
!<description>
  ! Copies quad precision vector to single precision vector: Fy = Qx
!</description>

!<input>
  
  ! Source vector
  real(QP), dimension(:,:,:), intent(in) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:,:,:), intent(out) :: Fy
  
!</output>
  
!</subroutine>
  integer :: i,j,k
  
    if (present(n) .and. present(m) .and. present(o)) then

      do k = 1, o
        do j = 1, m
          do i = 1, n
            Fy(i,j,k) = real(Qx(i,j,k),SP)
          end do
        end do
      end do

    else

      do k = 1, size(Qx,3)
        do j = 1, size(Qx,2)
          do i = 1, size(Qx,1)
            Fy(i,j,k) = real(Qx(i,j,k),SP)
          end do
        end do
      end do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQuadDbl3D (Qx,Dy,n,m,o)
  
!<description>
  ! Copies quad precision vector to double precision vector: Dy = Qx
!</description>

!<input>
  
  ! Source vector
  real(QP), dimension(:,:,:), intent(in) :: Qx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:,:,:), intent(out) :: Dy
  
!</output>
  
!</subroutine>
  integer :: i,j,k
  
    if (present(n) .and. present(m) .and. present(o)) then

      do k = 1, o
        do j = 1, m
          do i = 1, n
            Dy(i,j,k) = real(Qx(i,j,k),DP)
          end do
        end do
      end do

    else

      do k = 1, size(Qx,3)
        do j = 1, size(Qx,2)
          do i = 1, size(Qx,1)
            Dy(i,j,k) = real(Qx(i,j,k),DP)
          end do
        end do
      end do
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI8_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I8), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = Ix(i,j,k)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = Ix(i,j,k)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI8I16_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I8), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = int(Ix(i,j,k),I16)
          end do
        end do
      end do

    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = int(Ix(i,j,k),I16)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI8I32_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I8), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = int(Ix(i,j,k),I32)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = int(Ix(i,j,k),I32)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI8I64_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I8), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = int(Ix(i,j,k),I64)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = int(Ix(i,j,k),I64)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI16_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I16), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = Ix(i,j,k)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = Ix(i,j,k)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI16I8_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I16), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = int(Ix(i,j,k),I8)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = int(Ix(i,j,k),I8)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI16I32_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I16), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = int(Ix(i,j,k),I32)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = int(Ix(i,j,k),I32)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI16I64_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I16), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = int(Ix(i,j,k),I64)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = int(Ix(i,j,k),I64)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI32_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = Ix(i,j,k)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = Ix(i,j,k)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI32I8_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = int(Ix(i,j,k),I8)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = int(Ix(i,j,k),I8)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI32I16_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = int(Ix(i,j,k),I16)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = int(Ix(i,j,k),I16)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI32I64_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = int(Ix(i,j,k),I64)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = int(Ix(i,j,k),I64)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI64_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I64), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I64), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = Ix(i,j,k)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = Ix(i,j,k)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI64I8_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I64), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I8), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = int(Ix(i,j,k),I8)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = int(Ix(i,j,k),I8)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI64I16_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I64), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I16), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = int(Ix(i,j,k),I16)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = int(Ix(i,j,k),I16)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorI64I32_3D (Ix,Iy,n,m,o)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I64), dimension(:,:,:), intent(in) :: Ix
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:,:), intent(out) :: Iy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if(present(m) .and. present(n) .and. present(o)) then
    
      do k = 1, o
        do j = 1, m
          do i = 1, n
            Iy(i,j,k) = int(Ix(i,j,k),I32)
          end do
        end do
      end do
    
    else
    
      do k = 1, size(Ix,3)
        do j = 1, size(Ix,2)
          do i = 1, size(Ix,1)
            Iy(i,j,k) = int(Ix(i,j,k),I32)
          end do
        end do
      end do

    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorLogical3D (Lx,Ly,n,m,o)
  
!<description>
  ! Copies a logical vector Lx: Ly = Lx
!</description>

!<input>
  
  ! Source vector
  logical, dimension(:,:,:), intent(in) :: Lx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  logical, dimension(:,:,:), intent(out) :: Ly
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if (present(n) .and. present(m)) then

      do k = 1, o
        do j = 1, m
          do i = 1, n
            Ly(i,j,k) = Lx(i,j,k)
          end do
        end do
      end do

    else
      
      do k = 1, size(Lx,3)
        do j = 1, size(Lx,2)
          do i = 1, size(Lx,1)
            Ly(i,j,k) = Lx(i,j,k)
          end do
        end do
      end do

    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorChar3D (Sx,Sy,n,m,o)
  
!<description>
  ! Copies a character vector Sx: Sy = Sx
!</description>

!<input>
  
  ! Source vector
  character, dimension(:,:,:), intent(in) :: Sx
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n,m,o

!</input>

!<output>
  
  ! Destination vector
  character, dimension(:,:,:), intent(out) :: Sy
  
!</output>
  
!</subroutine>

  integer :: i,j,k
  
    if (present(n) .and. present(m)) then

      do k = 1, o
        do j = 1, m
          do i = 1, n
            Sy(i,j,k) = Sx(i,j,k)
          end do
        end do
      end do

    else

      do k = 1, size(Sx,3)
        do j = 1, size(Sx,2)
          do i = 1, size(Sx,1)
            Sy(i,j,k) = Sx(i,j,k)
          end do
        end do
      end do

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
  real(SP), dimension(:), intent(inout) :: Fx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(SP), intent(in) :: sc

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

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
  real(DP), dimension(:), intent(inout) :: Dx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(DP), intent(in) :: dc

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

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
  real(QP), dimension(:), intent(inout) :: Qx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(QP), intent(in) :: qc

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>
  
!</subroutine>

  integer :: i
  
    if (.not. present(n)) then
      if(qc .eq. 0.0_QP) then
        call lalg_clearVectorQuad(Qx)
      else if(qc .ne. 1.0_QP) then
      
        do i = 1, size(Qx)
          Qx(i) = qc * Qx(i)
        end do
      
        !call QSCAL(size(Qx),qc,Qx,1)
      end if
    else
      if(qc .eq. 0.0_QP) then
        call lalg_clearVectorQuad(Qx,n)
      else if(qc .ne. 1.0_QP) then

        do i = 1, n
          Qx(i) = qc * Qx(i)
        end do

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
  real(SP), dimension(:,:), intent(inout) :: Fx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(SP), intent(in) :: sc

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
  real(DP), dimension(:,:), intent(inout) :: Dx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(DP), intent(in) :: dc

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
  real(QP), dimension(:,:), intent(inout) :: Qx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(QP), intent(in) :: qc

!</input>
  
!</subroutine>

  integer :: i,j

    if(qc .eq. 0.0_QP) then
      call lalg_clearVectorQuad2D(Qx)
    else if(qc .ne. 1.0_QP) then
    
      do j = 1, size(Qx,2)
        do i = 1, size(Qx,1)
          Qx(i,j) = qc*Qx(i,j)
        end do
      end do
    
      !call QSCAL(size(Qx,1)*size(Qx,2),qc,Qx,1)
    end if
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lalg_scaleVectorSngl3D (Fx,sc)
  
!<description>
  ! Scales a single precision vector: Dx = sc * Dx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  real(SP), dimension(:,:,:), intent(inout) :: Fx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(SP), intent(in) :: sc

!</input>
  
!</subroutine>

    if(sc .eq. 0.0_DP) then
      call lalg_clearVectorSngl3D(Fx)
    else if(sc .ne. 1.0_SP) then
      call SSCAL(size(Fx,1)*size(Fx,2)*size(Fx,3),sc,Fx,1)
    end if
  
  end subroutine

! ***************************************************************************
  
!<subroutine>

  subroutine lalg_scaleVectorDble3D (Dx,dc)
  
!<description>
  ! Scales a double precision vector: Dx = dc * Dx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  real(DP), dimension(:,:,:), intent(inout) :: Dx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(DP), intent(in) :: dc

!</input>
  
!</subroutine>

    if(dc .eq. 0.0_DP) then
      call lalg_clearVectorDble3D(Dx)
    else if(dc .ne. 1.0_DP) then
      call DSCAL(size(Dx,1)*size(Dx,2)*size(Dx,3),dc,Dx,1)
    end if
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lalg_scaleVectorQuad3D (Qx,qc)
  
!<description>
  ! Scales a quad precision vector: Qx = dc * Qx
!</description>

!<inputoutput>
  
  ! Source and destination vector
  real(QP), dimension(:,:,:), intent(inout) :: Qx
  
!</inputoutput>

!<input>

  ! Multiplication factor
  real(QP), intent(in) :: qc

!</input>
  
!</subroutine>

  integer :: i,j,k

    if(qc .eq. 0.0_QP) then
      call lalg_clearVectorQuad3D(Qx)
    else if(qc .ne. 1.0_QP) then
    
      do k = 1, size(Qx,3)
        do j = 1, size(Qx,2)
          do i = 1, size(Qx,1)
            Qx(i,j,k) = qc*Qx(i,j,k)
          end do
        end do
      end do
    
      !call QSCAL(size(Qx,1)*size(Qx,2)*size(Qx,3),qc,Qx,1)
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
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  real(SP), dimension(:), intent(out) :: Fx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Fx)
        Fx(i) = 0.0_SP
      end do
      
    else
    
      do i = 1, n
        Fx(i) = 0.0_SP
      end do
      
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
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  real(DP), dimension(:), intent(out) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then
    
      do i = 1, size(Dx)
        Dx(i) = 0.0_DP
      end do
      
    else
    
      do i = 1, n
        Dx(i) = 0.0_DP
      end do
      
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
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  real(QP), dimension(:), intent(out) :: Qx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then
    
      do i = 1, size(Qx)
        Qx(i) = 0.0_QP
      end do
      
    else
    
      do i = 1, n
        Qx(i) = 0.0_QP
      end do
      
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
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  integer(I8), dimension(:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Ix(i) = 0_I8
      end do
      
    else
    
      do i = 1, size(Ix)
        Ix(i) = 0_I8
      end do
      
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
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  integer(I16), dimension(:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Ix(i) = 0_I16
      end do
      
    else
    
      do i = 1, size(Ix)
        Ix(i) = 0_I16
      end do
      
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
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  integer(I32), dimension(:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Ix(i) = 0_I32
      end do
      
    else
    
      do i = 1, size(Ix)
        Ix(i) = 0_I32
      end do
      
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
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be cleared
  integer(I64), dimension(:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Ix)
        Ix(i) = 0_I64
      end do
      
    else
    
      do i = 1, size(Ix)
        Ix(i) = 0_I64
      end do
      
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
  real(SP), dimension(:,:), intent(out) :: Fx
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Fx,2)
      do i = 1, size(Fx,1)
        Fx(i,j) = 0.0_SP
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorDble2D (Dx)
  
!<description>
  ! Clears a double precision vector: Dx = 0
!</description>

!<output>
  ! Destination vector to be cleared
  real(DP), dimension(:,:), intent(out) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Dx,2)
      do i = 1, size(Dx,1)
        Dx(i,j) = 0.0_DP
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorQuad2D (Qx)
  
!<description>
  ! Clears a quad precision vector: Qx = 0
!</description>

!<output>
  ! Destination vector to be cleared
  real(QP), dimension(:,:), intent(out) :: Qx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Qx,2)
      do i = 1, size(Qx,1)
        Qx(i,j) = 0.0_QP
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI8_2D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I8), dimension(:,:), intent(out) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Ix,2)
      do i = 1, size(Ix,1)
        Ix(i,j) = 0_I8
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI16_2D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I16), dimension(:,:), intent(out) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Ix,2)
      do i = 1, size(Ix,1)
        Ix(i,j) = 0_I16
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI32_2D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I32), dimension(:,:), intent(out) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Ix,2)
      do i = 1, size(Ix,1)
        Ix(i,j) = 0_I32
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI64_2D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I64), dimension(:,:), intent(out) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Ix,2)
      do i = 1, size(Ix,1)
        Ix(i,j) = 0_I32
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorSngl3D (Fx)
  
!<description>
  ! Clears a single precision vector: Fx = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  real(SP), dimension(:,:,:), intent(out) :: Fx
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Fx,3)
      do j = 1, size(Fx,2)
        do i = 1, size(Fx,1)
          Fx(i,j,k) = 0.0_SP
        end do
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorDble3D (Dx)
  
!<description>
  ! Clears a double precision vector: Dx = 0
!</description>

!<output>
  ! Destination vector to be cleared
  real(DP), dimension(:,:,:), intent(out) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Dx,3)
      do j = 1, size(Dx,2)
        do i = 1, size(Dx,1)
          Dx(i,j,k) = 0.0_DP
        end do
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorQuad3D (Qx)
  
!<description>
  ! Clears a quad precision vector: Qx = 0
!</description>

!<output>
  ! Destination vector to be cleared
  real(QP), dimension(:,:,:), intent(out) :: Qx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Qx,3)
      do j = 1, size(Qx,2)
        do i = 1, size(Qx,1)
          Qx(i,j,k) = 0.0_QP
        end do
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI8_3D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I8), dimension(:,:,:), intent(out) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Ix,3)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Ix(i,j,k) = 0_I8
        end do
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI16_3D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I16), dimension(:,:,:), intent(out) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Ix,3)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Ix(i,j,k) = 0_I16
        end do
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI32_3D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I32), dimension(:,:,:), intent(out) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Ix,3)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Ix(i,j,k) = 0_I32
        end do
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorI64_3D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I64), dimension(:,:,:), intent(out) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Ix,3)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Ix(i,j,k) = 0_I32
        end do
      end do
    end do
  
  end subroutine
 
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorSngl (Fx,fvalue,n)
  
!<description>
  ! Sets the vector data to a defined value: Fx = fvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  real(SP), intent(in) :: fvalue

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  real(SP), dimension(:), intent(out) :: Fx
!</output>

!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      do i = 1, size(Fx)
        Fx(i) = fvalue
      end do
      
    else
    
      do i = 1, n
        Fx(i) = fvalue
      end do
      
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
  real(DP), intent(in) :: dvalue

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  real(DP), dimension(:), intent(out) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then

      do i = 1, size(Dx)
        Dx(i) = dvalue
      end do
      
    else
    
      do i = 1, size(Dx)
        Dx(i) = dvalue
      end do
      
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
  real(QP), intent(in) :: qvalue

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  real(QP), dimension(:), intent(out) :: Qx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then

      do i = 1, size(Qx)
        Qx(i) = qvalue
      end do
      
    else
    
      do i = 1, size(Qx)
        Qx(i) = qvalue
      end do
      
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
  integer(I8), intent(in) :: ivalue

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  integer(I8), dimension(:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
   
      do i = 1, size(Ix)
        Ix(i) = ivalue
      end do
      
    else
    
      do i = 1, n
        Ix(i) = ivalue
      end do
      
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
  integer(I16), intent(in) :: ivalue

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  integer(I16), dimension(:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
   
      do i = 1, size(Ix)
        Ix(i) = ivalue
      end do
      
    else
    
      do i = 1, n
        Ix(i) = ivalue
      end do
      
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
  integer(I32), intent(in) :: ivalue

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  integer(I32), dimension(:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
   
      do i = 1, size(Ix)
        Ix(i) = ivalue
      end do
      
    else
    
      do i = 1, n
        Ix(i) = ivalue
      end do
      
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
  integer(I64), intent(in) :: ivalue

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  integer(I64), dimension(:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
   
      do i = 1, size(Ix)
        Ix(i) = ivalue
      end do
      
    else
    
      do i = 1, n
        Ix(i) = ivalue
      end do
      
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
  logical, intent(in) :: lvalue

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  logical, dimension(:), intent(out) :: Lx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Lx)
        Lx(i) = lvalue
      end do
      
    else
    
      do i = 1, n
        Lx(i) = lvalue
      end do
      
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
  character, intent(in) :: svalue

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
!</input>

!<output>
  ! Destination vector to be set
  character, dimension(:), intent(out) :: Sx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i
  
    if (.not. present(n)) then
    
      do i = 1, size(Sx)
        Sx(i) = svalue
      end do
      
    else
    
      do i = 1, n
        Sx(i) = svalue
      end do
      
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
  real(SP), intent(in) :: fvalue
!</input>

!<output>
  ! Destination vector to be set
  real(SP), dimension(:,:), intent(out) :: Fx
!</output>

!</subroutine>

  ! local variables
  integer :: i,j
  
    do j = 1, size(Fx,2)
      do i = 1, size(Fx,1)
        Fx(i,j) = fvalue
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorDble2D (Dx,dvalue)
  
!<description>
  ! Sets the vector data to a defined value: Dx = dvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  real(DP), intent(in) :: dvalue
!</input>

!<output>
  ! Destination vector to be set
  real(DP), dimension(:,:), intent(out) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Dx,2)
      do i = 1, size(Dx,1)
        Dx(i,j) = dvalue
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorQuad2D (Qx,qvalue)
  
!<description>
  ! Sets the vector data to a defined value: Qx = qvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  real(QP), intent(in) :: qvalue
!</input>

!<output>
  ! Destination vector to be set
  real(QP), dimension(:,:), intent(out) :: Qx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Qx,2)
      do i = 1, size(Qx,1)
        Qx(i,j) = qvalue
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI8_2D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I8), intent(in) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  integer(I8), dimension(:,:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Ix,2)
      do i = 1, size(Ix,1)
        Ix(i,j) = ivalue
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI16_2D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I16), intent(in) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  integer(I16), dimension(:,:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Ix,2)
      do i = 1, size(Ix,1)
        Ix(i,j) = ivalue
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI32_2D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I32), intent(in) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  integer(I32), dimension(:,:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Ix,2)
      do i = 1, size(Ix,1)
        Ix(i,j) = ivalue
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI64_2D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I64), intent(in) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  integer(I64), dimension(:,:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Ix,2)
      do i = 1, size(Ix,1)
        Ix(i,j) = ivalue
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorLogical2D (Lx,lvalue)
  
!<description>
  ! Sets the vector data to a defined value: Lx = lvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  logical, intent(in) :: lvalue
!</input>

!<output>
  ! Destination vector to be set
  logical, dimension(:,:), intent(out) :: Lx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j
  
    do j = 1, size(Lx,2)
      do i = 1, size(Lx,1)
        Lx(i,j) = lvalue
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorChar2D (Sx,svalue)
  
!<description>
  ! Sets the vector data to a defined value: Sx = svalue
!</description>

!<input>
  ! The value, the vector should be set to.
  character, intent(in) :: svalue
!</input>

!<output>
  ! Destination vector to be set
  character, dimension(:,:), intent(out) :: Sx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j

    do j = 1, size(Sx,2)
      do i = 1, size(Sx,1)
        Sx(i,j) = svalue
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorSngl3D (Fx,fvalue)
  
!<description>
  ! Sets the vector data to a defined value: Fx = fvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  real(SP), intent(in) :: fvalue
!</input>

!<output>
  ! Destination vector to be set
  real(SP), dimension(:,:,:), intent(out) :: Fx
!</output>

!</subroutine>

  ! local variables
  integer :: i,j,k
  
    do k = 1, size(Fx,3)
      do j = 1, size(Fx,2)
        do i = 1, size(Fx,1)
          Fx(i,j,k) = fvalue
        end do
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorDble3D (Dx,dvalue)
  
!<description>
  ! Sets the vector data to a defined value: Dx = dvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  real(DP), intent(in) :: dvalue
!</input>

!<output>
  ! Destination vector to be set
  real(DP), dimension(:,:,:), intent(out) :: Dx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Dx,3)
      do j = 1, size(Dx,2)
        do i = 1, size(Dx,1)
          Dx(i,j,k) = dvalue
        end do
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorQuad3D (Qx,qvalue)
  
!<description>
  ! Sets the vector data to a defined value: Qx = qvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  real(QP), intent(in) :: qvalue
!</input>

!<output>
  ! Destination vector to be set
  real(QP), dimension(:,:,:), intent(out) :: Qx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Qx,3)
      do j = 1, size(Qx,2)
        do i = 1, size(Qx,1)
          Qx(i,j,k) = qvalue
        end do
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI8_3D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I8), intent(in) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  integer(I8), dimension(:,:,:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Ix,3)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Ix(i,j,k) = ivalue
        end do
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI16_3D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I16), intent(in) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  integer(I16), dimension(:,:,:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Ix,3)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Ix(i,j,k) = ivalue
        end do
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI32_3D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I32), intent(in) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  integer(I32), dimension(:,:,:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Ix,3)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Ix(i,j,k) = ivalue
        end do
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorI64_3D (Ix,ivalue)
  
!<description>
  ! Sets the vector data to a defined value: Ix = ivalue
!</description>

!<input>
  ! The value, the vector should be set to.
  integer(I64), intent(in) :: ivalue
!</input>

!<output>
  ! Destination vector to be set
  integer(I64), dimension(:,:,:), intent(out) :: Ix
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Ix,3)
      do j = 1, size(Ix,2)
        do i = 1, size(Ix,1)
          Ix(i,j,k) = ivalue
        end do
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorLogical3D (Lx,lvalue)
  
!<description>
  ! Sets the vector data to a defined value: Lx = lvalue
!</description>

!<input>
  ! The value, the vector should be set to.
  logical, intent(in) :: lvalue
!</input>

!<output>
  ! Destination vector to be set
  logical, dimension(:,:,:), intent(out) :: Lx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k
  
    do k = 1, size(Lx,3)
      do j = 1, size(Lx,2)
        do i = 1, size(Lx,1)
          Lx(i,j,k) = lvalue
        end do
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorChar3D (Sx,svalue)
  
!<description>
  ! Sets the vector data to a defined value: Sx = svalue
!</description>

!<input>
  ! The value, the vector should be set to.
  character, intent(in) :: svalue
!</input>

!<output>
  ! Destination vector to be set
  character, dimension(:,:,:), intent(out) :: Sx
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,j,k

    do k = 1, size(Sx,3)
      do j = 1, size(Sx,2)
        do i = 1, size(Sx,1)
          Sx(i,j,k) = svalue
        end do
      end do
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSngl (Fx,Fy,scx,scy,n)
  
!<description>
  ! Performs a linear combination: Fy = scx * Fx  +  scy * Fy
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:), intent(in) :: Fx
  
  ! Scaling factor for Dx
  real(SP), intent(in)               :: scx

  ! Scaling factor for Dy
  real(SP), intent(in)               :: scy

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(SP), dimension(:), intent(inout) :: Fy
  
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
  real(DP), dimension(:), intent(in) :: Dx
  
  ! Scaling factor for Dx
  real(DP), intent(in)               :: dcx

  ! Scaling factor for Dy
  real(DP), intent(in)               :: dcy
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(DP), dimension(:), intent(inout) :: Dy
  
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
  real(QP), dimension(:), intent(in) :: Qx
  
  ! Scaling factor for Dx
  real(QP), intent(in)               :: qcx

  ! Scaling factor for Dy
  real(QP), intent(in)               :: qcy
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:), intent(inout) :: Qy
  
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
      
        do i = 1, k
          Qy(i) = 0.0_QP
        end do
      
      else
      
        do i = 1, k
          Qy(i) = qcy*Qy(i)
        end do
      
      end if
    
    else if(qcx .eq. 1.0_QP) then
    
      if(qcy .eq. 1.0_QP) then
      
        do i = 1, k
          Qy(i) = Qy(i) + Qx(i)
        end do
      
      else if(qcy .eq. 0.0_QP) then
      
        do i = 1, k
          Qy(i) = Qx(i)
        end do
      
      else
      
        do i = 1, k
          Qy(i) = qcy*Qy(i) + Qx(i)
        end do

      end if
    
    else
    
      if(qcy .eq. 1.0_QP) then
      
        do i = 1, k
          Qy(i) = Qy(i) + qcx*Qx(i)
        end do
      
      else if(qcy .eq. 0.0_QP) then
      
        do i = 1, k
          Qy(i) = qcx*Qx(i)
        end do
      
      else
      
        do i = 1, k
          Qy(i) = qcy*Qy(i) + qcx*Qx(i)
        end do

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
  real(SP), dimension(:), intent(in) :: Fx
  
  ! Scaling factor for Dx
  real(SP), intent(in)               :: scx

  ! Scaling factor for Dy
  real(DP), intent(in)               :: dcy
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(DP), dimension(:), intent(inout) :: Dy
  
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
        do i = 1, k
          Dy(i) = 0.0_DP
        end do
      
      else
      
        ! Call DSCAL
        call DSCAL(k,dcy,Dy,1)
      
      end if
    
    else if(scx .eq. 1.0_SP) then
    
      if(dcy .eq. 1.0_DP) then
      
        do i = 1, k
          Dy(i) = Dy(i) + real(Fx(i),DP)
        end do
      
      else if(dcy .eq. 0.0_DP) then
      
        do i = 1, k
          Dy(i) = real(Fx(i),DP)
        end do
      
      else
      
        do i = 1, k
          Dy(i) = dcy*Dy(i) + real(Fx(i),DP)
        end do

      end if
    
    else
    
      ! Convert scx to double precision
      dcx = real(scx,DP)

      if(dcy .eq. 1.0_DP) then
      
        do i = 1, k
          Dy(i) = Dy(i) + dcx*real(Fx(i),DP)
        end do
      
      else if(dcy .eq. 0.0_DP) then
      
        do i = 1, k
          Dy(i) = dcx*real(Fx(i),DP)
        end do
      
      else
      
        do i = 1, k
          Dy(i) = dcy*Dy(i) + dcx*real(Fx(i),DP)
        end do

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
  real(SP), dimension(:), intent(in) :: Fx
  
  ! Scaling factor for Dx
  real(SP), intent(in)               :: scx

  ! Scaling factor for Qy
  real(QP), intent(in)               :: qcy
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:), intent(inout) :: Qy
  
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
        do i = 1, k
          Qy(i) = 0.0_QP
        end do
      
      else
      
        do i = 1, k
          Qy(i) = qcy*Qy(i)
        end do
      
      end if
    
    else if(scx .eq. 1.0_SP) then
    
      if(qcy .eq. 1.0_QP) then
      
        do i = 1, k
          Qy(i) = Qy(i) + real(Fx(i),QP)
        end do
      
      else if(qcy .eq. 0.0_QP) then
      
        do i = 1, k
          Qy(i) = real(Fx(i),QP)
        end do
      
      else
      
        do i = 1, k
          Qy(i) = qcy*Qy(i) + real(Fx(i),QP)
        end do

      end if
    
    else
    
      ! Convert scx to quad precision
      qcx = real(scx,QP)

      if(qcy .eq. 1.0_QP) then
      
        do i = 1, k
          Qy(i) = Qy(i) + qcx*real(Fx(i),QP)
        end do
      
      else if(qcy .eq. 0.0_QP) then
      
        do i = 1, k
          Qy(i) = qcx*real(Fx(i),QP)
        end do
      
      else
      
        do i = 1, k
          Qy(i) = qcy*Qy(i) + qcx*real(Fx(i),QP)
        end do

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
  real(DP), dimension(:), intent(in) :: Dx
  
  ! Scaling factor for Dx
  real(DP), intent(in)               :: dcx

  ! Scaling factor for Qy
  real(QP), intent(in)               :: qcy
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:), intent(inout) :: Qy
  
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
        do i = 1, k
          Qy(i) = 0.0_QP
        end do
      
      else
      
        do i = 1, k
          Qy(i) = qcy*Qy(i)
        end do
      
      end if
    
    else if(dcx .eq. 1.0_DP) then
    
      if(qcy .eq. 1.0_QP) then
      
        do i = 1, k
          Qy(i) = Qy(i) + real(Dx(i),QP)
        end do
      
      else if(qcy .eq. 0.0_QP) then
      
        do i = 1, k
          Qy(i) = real(Dx(i),QP)
        end do
      
      else
      
        do i = 1, k
          Qy(i) = qcy*Qy(i) + real(Dx(i),QP)
        end do

      end if
    
    else
    
      ! Convert dcx to quad precision
      qcx = real(dcx,QP)

      if(qcy .eq. 1.0_QP) then
      
        do i = 1, k
          Qy(i) = Qy(i) + qcx*real(Dx(i),QP)
        end do
      
      else if(qcy .eq. 0.0_QP) then
      
        do i = 1, k
          Qy(i) = qcx*real(Dx(i),QP)
        end do
      
      else
      
        do i = 1, k
          Qy(i) = qcy*Qy(i) + qcx*real(Dx(i),QP)
        end do

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
  real(SP), dimension(:,:), intent(in) :: Fx
  
  ! Scaling factor for Dx
  real(SP), intent(in)                 :: scx

  ! Scaling factor for Dy
  real(SP), intent(in)                 :: scy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(SP), dimension(:,:), intent(inout) :: Fy
  
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
  real(DP), dimension(:,:), intent(in) :: Dx
  
  ! Scaling factor for Dx
  real(DP), intent(in)                 :: dcx

  ! Scaling factor for Dy
  real(DP), intent(in)                 :: dcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(DP), dimension(:,:), intent(inout) :: Dy
  
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
  real(QP), dimension(:,:), intent(in) :: Qx
  
  ! Scaling factor for Qx
  real(QP), intent(in)                 :: qcx

  ! Scaling factor for Qy
  real(QP), intent(in)                 :: qcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:,:), intent(inout) :: Qy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i,j
  
    do j = 1, size(Qx,2)
      do i = 1, size(Qx,1)
        Qy(i,j) = qcy*Qy(i,j) + qcx*Qx(i,j)
      end do
    end do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSnglDble2D (Fx,Dy,scx,dcy)
  
!<description>
  ! Performs a linear combination: Dy = scx * Fx  +  dcy * Dy
!</description>

!<input>
  
  ! First source vector
  real(SP), dimension(:,:), intent(in) :: Fx
  
  ! Scaling factor for Dx
  real(SP), intent(in)                 :: scx

  ! Scaling factor for Dy
  real(DP), intent(in)                 :: dcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(DP), dimension(:,:), intent(inout) :: Dy
  
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
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = 0.0_DP
          end do
        end do
      
      else
      
        ! Call DSCAL
        call DSCAL(size(Dy,1)*size(Dy,2),dcy,Dy,1)
      
      end if
    
    else if(scx .eq. 1.0_SP) then
    
      if(dcy .eq. 1.0_DP) then
      
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = Dy(i,j) + real(Fx(i,j),DP)
          end do
        end do
      
      else if(dcy .eq. 0.0_DP) then
      
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = real(Fx(i,j),DP)
          end do
        end do
      
      else
      
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = dcy*Dy(i,j) + real(Fx(i,j),DP)
          end do
        end do

      end if
    
    else
    
      ! Convert scx to double precision
      dcx = real(scx,DP)

      if(dcy .eq. 1.0_DP) then
      
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = Dy(i,j) + dcx*real(Fx(i,j),DP)
          end do
        end do
      
      else if(dcy .eq. 0.0_DP) then
      
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = dcx*real(Fx(i,j),DP)
          end do
        end do
      
      else
      
        do j = 1, size(Dy,2)
          do i = 1, size(Dy,1)
            Dy(i,j) = dcy*Dy(i,j) + dcx*real(Fx(i,j),DP)
          end do
        end do

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
  real(SP), dimension(:,:), intent(in) :: Fx
  
  ! Scaling factor for Dx
  real(SP), intent(in)                 :: scx

  ! Scaling factor for Qy
  real(QP), intent(in)                 :: qcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:,:), intent(inout) :: Qy
  
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
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = 0.0_QP
          end do
        end do
      
      else
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcy*Qy(i,j)
          end do
        end do
      
      end if
    
    else if(scx .eq. 1.0_SP) then
    
      if(qcy .eq. 1.0_QP) then
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = Qy(i,j) + real(Fx(i,j),QP)
          end do
        end do
      
      else if(qcy .eq. 0.0_QP) then
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = real(Fx(i,j),QP)
          end do
        end do
      
      else
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcy*Qy(i,j) + real(Fx(i,j),QP)
          end do
        end do

      end if
    
    else
    
      ! Convert scx to quad precision
      qcx = real(scx,QP)

      if(qcy .eq. 1.0_QP) then
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = Qy(i,j) + qcx*real(Fx(i,j),QP)
          end do
        end do
      
      else if(qcy .eq. 0.0_QP) then
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcx*real(Fx(i,j),QP)
          end do
        end do
      
      else
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcy*Qy(i,j) + qcx*real(Fx(i,j),QP)
          end do
        end do

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
  real(DP), dimension(:,:), intent(in) :: Dx
  
  ! Scaling factor for Dx
  real(DP), intent(in)                 :: dcx

  ! Scaling factor for Qy
  real(QP), intent(in)                 :: qcy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:,:), intent(inout) :: Qy
  
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
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = 0.0_QP
          end do
        end do
      
      else
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcy*Qy(i,j)
          end do
        end do
      
      end if
    
    else if(dcx .eq. 1.0_DP) then
    
      if(qcy .eq. 1.0_QP) then
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = Qy(i,j) + real(Dx(i,j),QP)
          end do
        end do
      
      else if(qcy .eq. 0.0_QP) then
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = real(Dx(i,j),QP)
          end do
        end do
      
      else
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcy*Qy(i,j) + real(Dx(i,j),QP)
          end do
        end do

      end if
    
    else
    
      ! Convert dcx to quad precision
      qcx = real(dcx,QP)

      if(qcy .eq. 1.0_QP) then
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = Qy(i,j) + qcx*real(Dx(i,j),QP)
          end do
        end do
      
      else if(qcy .eq. 0.0_QP) then
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcx*real(Dx(i,j),QP)
          end do
        end do
      
      else
      
        do j = 1, size(Qy,2)
          do i = 1, size(Qy,1)
            Qy(i,j) = qcy*Qy(i,j) + qcx*real(Dx(i,j),QP)
          end do
        end do

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
  real(SP), dimension(:), intent(in) :: Fx
  
  ! Second source vector
  real(SP), dimension(:), intent(in) :: Fy

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
  
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
  real(DP), dimension(:), intent(in) :: Dx
  
  ! Second source vector
  real(DP), dimension(:), intent(in) :: Dy

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
  
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
  real(QP), dimension(:), intent(in) :: Qx
  
  ! Second source vector
  real(QP), dimension(:), intent(in) :: Qy
  
  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  !real(DP) :: QDOT
  integer :: i

    res = 0.0_QP
    if (.not. present(n)) then
    
      do i = 1, size(Qx)
        res = res + Qx(i)*Qy(i)
      end do
      
      !res = QDOT(size(Qx),Qx,1,Qy,1)
    else
    
      do i = 1, n
        res = res + Qx(i)*Qy(i)
      end do
      
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
  real(SP), dimension(:,:), intent(in) :: Fx
  
  ! Second source vector
  real(SP), dimension(:,:), intent(in) :: Fy
  
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
  real(DP), dimension(:,:), intent(in) :: Dx
  
  ! Second source vector
  real(DP), dimension(:,:), intent(in) :: Dy
  
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
  real(QP), dimension(:,:), intent(in) :: Qx
  
  ! Second source vector
  real(QP), dimension(:,:), intent(in) :: Qy
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  !real(QP) :: QDOT
  integer :: i,j
  
    res = 0.0_QP
  
    do j = 1, size(Qx,2)
      do i = 1, size(Qx,1)
        res = res + Qx(i,j)*Qy(i,j)
      end do
    end do

    !res = QDOT(size(Qx,1)*size(Qx,2),Qx,1,Qy,1)
  
  end function

  ! ***************************************************************************

!<function>

  real(SP) function lalg_normSngl (Fx,cnorm,iposMax,n) result(resnorm)
  
!<description>
  ! Calculates the norm of a single precision vector. cnorm identifies the 
  ! type of norm to calculate.
!</description>

!<input>
  ! Vector to calculate the norm of.
  real(SP), dimension(:), intent(in) :: Fx
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(in) :: cnorm

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(out), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! -1, if an error occurred (unknown norm).
!</result>
  
!</function>

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
      ! So, scale such that the vector (1111...) has norm = 1.
      resnorm = SASUM(isize,Fx,1) / real(isize,SP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by sqrt(vector length).
      ! So, scale such that the vector (1111...) has norm = 1.
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

!<function>

  real(DP) function lalg_normDble (Dx,cnorm,iposMax,n) result(resnorm)
  
!<description>
  ! Calculates the norm of a double precision vector. cnorm identifies the 
  ! type of norm to calculate.
!</description>

!<input>
  ! Vector to calculate the norm of.
  real(DP), dimension(:), intent(in) :: Dx
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(in) :: cnorm

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(out), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! -1, if an error occurred (unknown norm).
!</result>
  
!</function>

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
      ! So, scale such that the vector (1111...) has norm = 1.
      resnorm = DASUM(isize,Dx,1) / real(isize,DP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vector (1111...) has norm = 1.
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

!<function>

  real(QP) function lalg_normQuad (Qx,cnorm,iposMax,n) result(resnorm)
  
!<description>
  ! Calculates the norm of a quad precision vector. cnorm identifies the 
  ! type of norm to calculate.
!</description>

!<input>
  ! Vector to calculate the norm of.
  real(QP), dimension(:), intent(in) :: Qx
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(in) :: cnorm

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(out), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! -1, if an error occurred (unknown norm).
!</result>
  
!</function>

  !real(QP) :: QASUM,QNRM2,IQAMAX

  integer :: i,isize
  
    isize = size(Qx)
    if (present(n)) isize = n
    
    resnorm = 0.0_QP

    ! Choose the norm to calculate
    select case (cnorm)
    case (LINALG_NORMSUM)
      ! L1-norm: sum all entries
      do i = 1, isize
        resnorm = resnorm + abs(Qx(i))
      end do
      !resnorm = QASUM(isize,Qx,1)

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      do i = 1, isize
        resnorm = resnorm + Qx(i)*Qx(i)
      end do
      resnorm = sqrt(resnorm)
      !resnorm = QNRM2(isize,Qx,1)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vector (1111...) has norm = 1.
      do i = 1, isize
        resnorm = resnorm + abs(Qx(i))
      end do
      resnorm = resnorm / real(isize,QP)
      !resnorm = QASUM(isize,Dx,1) / real(isize,QP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vector (1111...) has norm = 1.
      do i = 1, isize
        resnorm = resnorm + Qx(i)*Qx(i)
      end do
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
        do i = 1, isize
          resnorm = max(resnorm,abs(Qx(i)))
        end do
      end if
        
    case default
      ! Unknown norm
      resnorm = -1.0_QP
      
    end select
    
  end function

  ! ***************************************************************************

!<function>

  real(SP) function lalg_errorNormSngl (Fx,Fy,cnorm,iposMax,n) result(resnorm)
  
!<description>
  ! Calculates the norm of two double precision vectors, !!Fx-Fy!!
  ! cnorm identifies the type of norm to calculate.
!</description>

!<input>
  ! Vectors to calculate the norm of their difference
  real(SP), dimension(:), intent(in) :: Fx,Fy
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(in) :: cnorm

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(out), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! < 0 >1, if an error occurred (unknown norm).
!</result>
  
!</function>

  real(SP) :: stemp
  integer :: i,j,isize
  
    isize = size(Fx)
    if (present(n)) isize = n
    
    resnorm = 0.0_SP

    ! Choose the norm to calculate
    select case (cnorm)
    case (LINALG_NORMSUM)
      ! L1-norm: sum all entries
      do i = 1, isize
        resnorm = resnorm + abs(Fx(i)-Fy(i))
      end do

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      do i = 1, isize
        resnorm = resnorm + (Fx(i)-Fy(i))*(Fx(i)-Fy(i))
      end do
      resnorm = sqrt(resnorm)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vector (1111...) to has norm = 1.
      do i = 1, isize
        resnorm = resnorm + abs(Fx(i)-Fy(i))
      end do
      resnorm = resnorm / real(isize,SP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vector (1111...) has norm = 1.
      do i = 1, isize
        resnorm = resnorm + (Fx(i)-Fy(i))*(Fx(i)-Fy(i))
      end do
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

!<function>

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
  real(DP), dimension(:), intent(in) :: Dx,Dy
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(in) :: cnorm

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

  ! OPTIONAL: Weighting vector
  real(DP), dimension(:), intent(in), optional :: Dw
!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(out), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! -1, if an error occurred (unknown norm).
!</result>
  
!</function>

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
        do i = 1, isize
          resnorm = resnorm + Dw(i)*abs(Dx(i)-Dy(i))
        end do
      else
        do i = 1, isize
          resnorm = resnorm + abs(Dx(i)-Dy(i))
        end do
      end if

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      if (present(Dw)) then
        do i = 1, isize
          resnorm = resnorm + Dw(i)*(Dx(i)-Dy(i))*(Dx(i)-Dy(i))
        end do
      else
        do i = 1, isize
          resnorm = resnorm + (Dx(i)-Dy(i))*(Dx(i)-Dy(i))
        end do
      end if
      resnorm = sqrt(resnorm)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vector (1111...) has norm = 1.
      do i = 1, isize
        resnorm = resnorm + abs(Dx(i)-Dy(i))
      end do
      resnorm = resnorm / real(isize,DP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vector (1111...) has norm = 1.
      do i = 1, isize
        resnorm = resnorm + (Dx(i)-Dy(i))*(Dx(i)-Dy(i))
      end do
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

!<function>

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
  real(QP), dimension(:), intent(in) :: Qx,Qy
  
  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(in) :: cnorm

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n

  ! OPTIONAL: Weighting vector
  real(QP), dimension(:), intent(in), optional :: Qw
!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(out), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! -1, if an error occurred (unknown norm).
!</result>
  
!</function>

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
        do i = 1, isize
          resnorm = resnorm + Qw(i)*abs(Qx(i)-Qy(i))
        end do
      else
        do i = 1, isize
          resnorm = resnorm + abs(Qx(i)-Qy(i))
        end do
      end if

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      if (present(Qw)) then
        do i = 1, isize
          resnorm = resnorm + Qw(i)*(Qx(i)-Qy(i))*(Qx(i)-Qy(i))
        end do
      else
        do i = 1, isize
          resnorm = resnorm + (Qx(i)-Qy(i))*(Qx(i)-Qy(i))
        end do
      end if
      resnorm = sqrt(resnorm)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vector (1111...) has norm = 1.
      do i = 1, isize
        resnorm = resnorm + abs(Qx(i)-Qy(i))
      end do
      resnorm = resnorm / real(isize,DP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vector (1111...) has norm = 1.
      do i = 1, isize
        resnorm = resnorm + (Qx(i)-Qy(i))*(Qx(i)-Qy(i))
      end do
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
  integer, dimension(:), intent(in) :: Itr

  ! Source vector to be sorted
  real(SP), dimension(:), intent(in) :: Fx
!</input>
  
!<output>
  ! The resorted vector
  real(SP), dimension(:), intent(out) :: Fd
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
  real(DP), dimension(:), intent(in) :: Dx
  
  ! Array with permutation of 1..neq.
  ! Itr(i) defines the number of the entry in Dx that should
  ! move to position i.
  integer, dimension(:), intent(in) :: Itr
!</input>
  
!<output>
  ! The resorted vector
  real(DP), dimension(:), intent(out) :: Dd
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
  real(QP), dimension(:), intent(in) :: Qx
  
  ! Array with permutation of 1..neq.
  ! Itr(i) defines the number of the entry in Dx that should
  ! move to position i.
  integer, dimension(:), intent(in) :: Itr
!</input>
  
!<output>
  ! The resorted vector
  real(QP), dimension(:), intent(out) :: Qd
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
  integer, dimension(:), intent(in) :: Itr

  ! Source vector to be sorted
  integer(I32), dimension(:), intent(in) :: Ix
!</input>
  
!<output>
  ! The resorted vector
  integer(I32), dimension(:), intent(out) :: Id
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
  integer, dimension(:), intent(in) :: Itr

  ! Source vector to be sorted
  integer(I64), dimension(:), intent(in) :: Ix
!</input>
  
!<output>
  ! The resorted vector
  integer(I64), dimension(:), intent(out) :: Id
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
  real(SP), dimension(:), intent(in) :: Sx

  ! Second source vector
  real(SP), dimension(:), intent(in) :: Sy
!</input>

!<output>
  ! Tensor product
  real(SP), dimension(:,:), intent(out) :: Stensor
!</output>
!</subroutine>
    
  ! local variables
  integer :: i,j

    do i = 1, size(Sy)
      do j = 1, size(Sx)
        Stensor(j,i) = Sx(j)*Sy(i)
      end do
    end do
      
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
  real(DP), dimension(:), intent(in) :: Dx

  ! Second source vector
  real(DP), dimension(:), intent(in) :: Dy
!</input>

!<output>
  ! Tensor product
  real(DP), dimension(:,:), intent(out) :: Dtensor
!</output>
!</subroutine>
    
  ! local variables
  integer :: i,j

    do i = 1, size(Dy)
      do j = 1, size(Dx)
        Dtensor(j,i) = Dx(j)*Dy(i)
      end do
    end do
      
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
  real(QP), dimension(:), intent(in) :: Qx

  ! Second source vector
  real(QP), dimension(:), intent(in) :: Qy
!</input>

!<output>
  ! Tensor product
  real(QP), dimension(:,:), intent(out) :: Qtensor
!</output>
!</subroutine>
    
  ! local variables
  integer :: i,j

    do i = 1, size(Qy)
      do j = 1, size(Qx)
        Qtensor(j,i) = Qx(j)*Qy(i)
      end do
    end do
      
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
  integer, intent(in), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  real(SP), dimension(:), intent(inout) :: Fx
!</inputoutput>

!</subroutine>
    
  integer :: i
  
    ! Nothing to do?
    if(fvalue .eq. 0.0_SP) return
    
    if (.not. present(n)) then
      
      do i = 1, size(Fx)
        Fx(i) = Fx(i) + fvalue
      end do
      
    else
      
      do i = 1, n
        Fx(i) = Fx(i) + fvalue
      end do
      
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
  integer, intent(in), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  real(DP), dimension(:), intent(inout) :: Dx
!</inputoutput>

!</subroutine>
    
  integer :: i
  
    ! Nothing to do?
    if(dvalue .eq. 0.0_DP) return
    
    if (.not. present(n)) then
    
      do i = 1, size(Dx)
        Dx(i) = Dx(i) + dvalue
      end do
      
    else
    
      do i = 1, n
        Dx(i) = Dx(i) + dvalue
      end do
      
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
  integer, intent(in), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  real(QP), dimension(:), intent(inout) :: Qx
!</inputoutput>

!</subroutine>
    
  integer :: i
  
    ! Nothing to do?
    if(qvalue .eq. 0.0_DP) return
    
    if (.not. present(n)) then
    
      do i = 1, size(Qx)
        Qx(i) = Qx(i) + qvalue
      end do
      
    else
    
      do i = 1, n
        Qx(i) = Qx(i) + qvalue
      end do
      
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
  real(SP), dimension(:,:), intent(inout) :: Fx
!</inputoutput>

!</subroutine>
    
  integer :: i,j
  
    ! Nothing to do?
    if(fvalue .eq. 0.0_SP) return
    
    do j = 1, size(Fx,2)
      do i = 1, size(Fx,1)
        Fx(i,j) = Fx(i,j) + fvalue
      end do
    end do
    
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
  real(DP), dimension(:,:), intent(inout) :: Dx
!</inputoutput>

!</subroutine>
    
  integer :: i,j
  
    ! Nothing to do?
    if(dvalue .eq. 0.0_DP) return
    
    do j = 1, size(Dx,2)
      do i = 1, size(Dx,1)
        Dx(i,j) = Dx(i,j) + dvalue
      end do
    end do
    
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
  real(QP), dimension(:,:), intent(inout) :: Qx
!</inputoutput>

!</subroutine>
    
  integer :: i,j
  
    ! Nothing to do?
    if(qvalue .eq. 0.0_QP) return
    
    do j = 1, size(Qx,2)
      do i = 1, size(Qx,1)
        Qx(i,j) = Qx(i,j) + qvalue
      end do
    end do
    
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
  integer, intent(in), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  integer(I8), dimension(:), intent(inout) :: Ix
!</inputoutput>

!</subroutine>
    
    integer :: i
    
    if (.not. present(n)) then
      do i = 1, size(Ix)
        Ix(i) = Ix(i) + ivalue
      end do
    else
      do i = 1, n
        Ix(i) = Ix(i) + ivalue
      end do
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
  integer, intent(in), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  integer(I16), dimension(:), intent(inout) :: Ix
!</inputoutput>

!</subroutine>
    
    integer :: i
    
    if (.not. present(n)) then
      do i = 1, size(Ix)
        Ix(i) = Ix(i) + ivalue
      end do
    else
      do i = 1, n
        Ix(i) = Ix(i) + ivalue
      end do
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
  integer, intent(in), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  integer(I32), dimension(:), intent(inout) :: Ix
!</inputoutput>

!</subroutine>
    
    integer :: i
    
    if (.not. present(n)) then
      do i = 1, size(Ix)
        Ix(i) = Ix(i) + ivalue
      end do
    else
      do i = 1, n
        Ix(i) = Ix(i) + ivalue
      end do
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
  integer, intent(in), optional :: n
!</input>

!<inputoutput>
  ! Data array to be modified.
  integer(I64), dimension(:), intent(inout) :: Ix
!</inputoutput>

!</subroutine>
    
    integer :: i
    
    if (.not. present(n)) then
      do i = 1, size(Ix)
        Ix(i) = Ix(i) + ivalue
      end do
    else
      do i = 1, n
        Ix(i) = Ix(i) + ivalue
      end do
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
  real(SP), dimension(:), intent(in) :: Fx
  
  ! Scaling factor
  real(SP), intent(in)               :: sc

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  real(SP), dimension(:), intent(inout) :: Fy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      do i = 1, size(Fy,1)
        Fy(i) = sc*Fx(i)*Fy(i)
      end do

    else

      do i = 1, n
        Fy(i) = sc*Fx(i)*Fy(i)
      end do

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
  real(DP), dimension(:), intent(in) :: Dx
  
  ! Scaling factor
  real(DP), intent(in)               :: dc

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  real(DP), dimension(:), intent(inout) :: Dy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      do i = 1, size(Dy,1)
        Dy(i) = dc*Dx(i)*Dy(i)
      end do

    else

      do i = 1, n
        Dy(i) = dc*Dx(i)*Dy(i)
      end do

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
  real(QP), dimension(:), intent(in) :: Qx
  
  ! Scaling factor
  real(QP), intent(in)               :: qc

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:), intent(inout) :: Qy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      do i = 1, size(Qy,1)
        Qy(i) = qc*Qx(i)*Qy(i)
      end do

    else

      do i = 1, n
        Qy(i) = qc*Qx(i)*Qy(i)
      end do

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
  real(SP), dimension(:), intent(in) :: Fx
  
  ! Scaling factor
  real(DP), intent(in)               :: dc

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  real(DP), dimension(:), intent(inout) :: Dy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      do i = 1, size(Dy,1)
        Dy(i) = dc*Fx(i)*Dy(i)
      end do

    else

      do i = 1, n
        Dy(i) = dc*Fx(i)*Dy(i)
      end do

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
  real(SP), dimension(:), intent(in) :: Fx
  
  ! Scaling factor
  real(QP), intent(in)               :: qc

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:), intent(inout) :: Qy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      do i = 1, size(Qy,1)
        Qy(i) = qc*real(Fx(i),QP)*Qy(i)
      end do

    else

      do i = 1, n
        Qy(i) = qc*real(Fx(i),QP)*Qy(i)
      end do

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
  real(DP), dimension(:), intent(in) :: Dx
  
  ! Scaling factor
  real(QP), intent(in)               :: qc

  ! OPTIONAL: Size of the vector
  integer, intent(in), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  real(QP), dimension(:), intent(inout) :: Qy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i

    if (.not. present(n)) then

      do i = 1, size(Qy,1)
        Qy(i) = qc*real(Dx(i),QP)*Qy(i)
      end do

    else

      do i = 1, n
        Qy(i) = qc*real(Dx(i),QP)*Qy(i)
      end do

    end if

  end subroutine
 
end module

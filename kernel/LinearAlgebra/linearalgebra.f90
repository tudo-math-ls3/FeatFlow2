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
!#
!# 12.) lalg_initPerfConfig
!#       -> Initialises the global performance configuration
!# </purpose>
!##############################################################################

module linearalgebra

!$ use omp_lib
  use fsystem
  use perfconfig

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
    module procedure lalg_copyVectorSP
    module procedure lalg_copyVectorDP
    module procedure lalg_copyVectorSPDP
    module procedure lalg_copyVectorDPSP
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
    module procedure lalg_copyVectorSP2D
    module procedure lalg_copyVectorDP2D
    module procedure lalg_copyVectorSPDP2D
    module procedure lalg_copyVectorDPSP2D
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
    module procedure lalg_copyVectorSP3D
    module procedure lalg_copyVectorDP3D
    module procedure lalg_copyVectorSPDP3D
    module procedure lalg_copyVectorDPSP3D
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
    module procedure lalg_copyVectorQP
    module procedure lalg_copyVectorSPQP
    module procedure lalg_copyVectorDPQP
    module procedure lalg_copyVectorQPSP
    module procedure lalg_copyVectorQPDP
    module procedure lalg_copyVectorQP2D
    module procedure lalg_copyVectorSPQP2D
    module procedure lalg_copyVectorDPQP2D
    module procedure lalg_copyVectorQPSP2D
    module procedure lalg_copyVectorQPDP2D
    module procedure lalg_copyVectorQP3D
    module procedure lalg_copyVectorSPQP3D
    module procedure lalg_copyVectorDPQP3D
    module procedure lalg_copyVectorQPSP3D
    module procedure lalg_copyVectorQPDP3D
#endif
  end interface

  public :: lalg_initPerfConfig

  public :: lalg_copyVector

  public :: lalg_copyVectorSP
  public :: lalg_copyVectorDP
  public :: lalg_copyVectorQP
  public :: lalg_copyVectorSPDP
  public :: lalg_copyVectorSPQP
  public :: lalg_copyVectorDPSP
  public :: lalg_copyVectorDPQP
  public :: lalg_copyVectorQPSP
  public :: lalg_copyVectorQPDP
  public :: lalg_copyVectorLogical
  public :: lalg_copyVectorChar

  public :: lalg_copyVectorSP2D
  public :: lalg_copyVectorDP2D
  public :: lalg_copyVectorQP2D
  public :: lalg_copyVectorSPDP2D
  public :: lalg_copyVectorSPQP2D
  public :: lalg_copyVectorDPSP2D
  public :: lalg_copyVectorDPQP2D
  public :: lalg_copyVectorQPSP2D
  public :: lalg_copyVectorQPDP2D
  public :: lalg_copyVectorLogical2D
  public :: lalg_copyVectorChar2D

  public :: lalg_copyVectorSP3D
  public :: lalg_copyVectorDP3D
  public :: lalg_copyVectorQP3D
  public :: lalg_copyVectorSPDP3D
  public :: lalg_copyVectorSPQP3D
  public :: lalg_copyVectorDPSP3D
  public :: lalg_copyVectorDPQP3D
  public :: lalg_copyVectorQPSP3D
  public :: lalg_copyVectorQPDP3D
  public :: lalg_copyVectorLogical3D
  public :: lalg_copyVectorChar3D

  interface lalg_scaleVector
    module procedure lalg_scaleVectorSP
    module procedure lalg_scaleVectorDP
    module procedure lalg_scaleVectorSP2D
    module procedure lalg_scaleVectorDP2D
    module procedure lalg_scaleVectorSP3D
    module procedure lalg_scaleVectorDP3D
#ifdef ENABLE_QUADPREC
    module procedure lalg_scaleVectorQP
    module procedure lalg_scaleVectorQP2D
    module procedure lalg_scaleVectorQP3D
#endif
  end interface

  public :: lalg_scaleVector
  public :: lalg_scaleVectorSP
  public :: lalg_scaleVectorDP
  public :: lalg_scaleVectorQP
  public :: lalg_scaleVectorSP2D
  public :: lalg_scaleVectorDP2D
  public :: lalg_scaleVectorQP2D
  public :: lalg_scaleVectorSP3D
  public :: lalg_scaleVectorDP3D
  public :: lalg_scaleVectorQP3D

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
    module procedure lalg_clearVectorSP
    module procedure lalg_clearVectorDP
    module procedure lalg_clearVectorI8
    module procedure lalg_clearVectorI16
    module procedure lalg_clearVectorI32
    module procedure lalg_clearVectorI64
    module procedure lalg_clearVectorSP2D
    module procedure lalg_clearVectorDP2D
    module procedure lalg_clearVectorI8_2D
    module procedure lalg_clearVectorI16_2D
    module procedure lalg_clearVectorI32_2D
    module procedure lalg_clearVectorI64_2D
    module procedure lalg_clearVectorSP3D
    module procedure lalg_clearVectorDP3D
    module procedure lalg_clearVectorI8_3D
    module procedure lalg_clearVectorI16_3D
    module procedure lalg_clearVectorI32_3D
    module procedure lalg_clearVectorI64_3D
#ifdef ENABLE_QUADPREC
    module procedure lalg_clearVectorQP
    module procedure lalg_clearVectorQP2D
    module procedure lalg_clearVectorQP3D
#endif
  end interface

  public :: lalg_clearVector
  public :: lalg_clearVectorSP
  public :: lalg_clearVectorDP
  public :: lalg_clearVectorQP
  public :: lalg_clearVectorSP2D
  public :: lalg_clearVectorDP2D
  public :: lalg_clearVectorQP2D
  public :: lalg_clearVectorSP3D
  public :: lalg_clearVectorDP3D
  public :: lalg_clearVectorQP3D

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
    module procedure lalg_setVectorSP
    module procedure lalg_setVectorDP
    module procedure lalg_setVectorI8
    module procedure lalg_setVectorI16
    module procedure lalg_setVectorI32
    module procedure lalg_setVectorI64
    module procedure lalg_setVectorLogical
    module procedure lalg_setVectorChar
    module procedure lalg_setVectorSP2D
    module procedure lalg_setVectorDP2D
    module procedure lalg_setVectorI8_2D
    module procedure lalg_setVectorI16_2D
    module procedure lalg_setVectorI32_2D
    module procedure lalg_setVectorI64_2D
    module procedure lalg_setVectorLogical2D
    module procedure lalg_setVectorChar2D
    module procedure lalg_setVectorSP3D
    module procedure lalg_setVectorDP3D
    module procedure lalg_setVectorI8_3D
    module procedure lalg_setVectorI16_3D
    module procedure lalg_setVectorI32_3D
    module procedure lalg_setVectorI64_3D
    module procedure lalg_setVectorLogical3D
    module procedure lalg_setVectorChar3D
#ifdef ENABLE_QUADPREC
    module procedure lalg_setVectorQP
    module procedure lalg_setVectorQP2D
    module procedure lalg_setVectorQP3D
#endif
  end interface

  public :: lalg_setVector

  public :: lalg_setVectorQP
  public :: lalg_setVectorQP2D
  public :: lalg_setVectorQP3D
  public :: lalg_setVectorSP
  public :: lalg_setVectorSP2D
  public :: lalg_setVectorSP3D
  public :: lalg_setVectorDP
  public :: lalg_setVectorDP2D
  public :: lalg_setVectorDP3D
  public :: lalg_setVectorLogical
  public :: lalg_setVectorLogical2D
  public :: lalg_setVectorLogical3D
  public :: lalg_setVectorChar
  public :: lalg_setVectorChar2D
  public :: lalg_setVectorChar3D

  interface lalg_vectorLinearComb
    module procedure lalg_vectorLinearCombSP
    module procedure lalg_vectorLinearCombDP
    module procedure lalg_vectorLinearCombSPDP
    module procedure lalg_vectorLinearCombSP2D
    module procedure lalg_vectorLinearCombDP2D
    module procedure lalg_vectorLinearCombSPDP2D
    module procedure lalg_vectorLinearCombSP3D
    module procedure lalg_vectorLinearCombDP3D
    module procedure lalg_vectorLinearCombSPDP3D
#ifdef ENABLE_QUADPREC
    module procedure lalg_vectorLinearCombQP
    module procedure lalg_vectorLinearCombSPQP
    module procedure lalg_vectorLinearCombDPQP
    module procedure lalg_vectorLinearCombQP2D
    module procedure lalg_vectorLinearCombSPQP2D
    module procedure lalg_vectorLinearCombDPQP2D
    module procedure lalg_vectorLinearCombQP3D
    module procedure lalg_vectorLinearCombSPQP3D
    module procedure lalg_vectorLinearCombDPQP3D
#endif
  end interface

  public :: lalg_vectorLinearComb
  public :: lalg_vectorLinearCombSP
  public :: lalg_vectorLinearCombDP
  public :: lalg_vectorLinearCombQP
  public :: lalg_vectorLinearCombSPDP
  public :: lalg_vectorLinearCombSPQP
  public :: lalg_vectorLinearCombDPQP
  public :: lalg_vectorLinearCombSP2D
  public :: lalg_vectorLinearCombDP2D
  public :: lalg_vectorLinearCombQP2D
  public :: lalg_vectorLinearCombSPDP2D
  public :: lalg_vectorLinearCombSPQP2D
  public :: lalg_vectorLinearCombDPQP2D
  public :: lalg_vectorLinearCombSP3D
  public :: lalg_vectorLinearCombDP3D
  public :: lalg_vectorLinearCombQP3D
  public :: lalg_vectorLinearCombSPDP3D
  public :: lalg_vectorLinearCombSPQP3D
  public :: lalg_vectorLinearCombDPQP3D

  interface lalg_scalarProduct
    module procedure lalg_scalarProductSP
    module procedure lalg_scalarProductDP
    module procedure lalg_scalarProductSP2D
    module procedure lalg_scalarProductDP2D
#ifdef ENABLE_QUADPREC
    module procedure lalg_scalarProductQP
    module procedure lalg_scalarProductQP2D
#endif
  end interface

  public :: lalg_scalarProduct
  public :: lalg_scalarProductSP
  public :: lalg_scalarProductDP
  public :: lalg_scalarProductQP
  public :: lalg_scalarProductSP2D
  public :: lalg_scalarProductDP2D
  public :: lalg_scalarProductQP2D

  interface lalg_norm
    module procedure lalg_normSP
    module procedure lalg_normDP
#ifdef ENABLE_QUADPREC
    module procedure lalg_normQP
#endif
  end interface

  public :: lalg_norm,lalg_normSP,lalg_normDP,lalg_normQP

  interface lalg_errorNorm
    module procedure lalg_errorNormSP
    module procedure lalg_errorNormDP
#ifdef ENABLE_QUADPREC
    module procedure lalg_errorNormQP
#endif
  end interface

  public :: lalg_errorNorm
  public :: lalg_errorNormSP
  public :: lalg_errorNormDP
  public :: lalg_errorNormQP

  interface lalg_vectorSortInt
    module procedure lalg_vectorSortI32
    module procedure lalg_vectorSortI64
  end interface

  public :: lalg_vectorSortInt

  interface lalg_vectorSort
    module procedure lalg_vectorSortSP
    module procedure lalg_vectorSortDP
    module procedure lalg_vectorSortI32
    module procedure lalg_vectorSortI64
#ifdef ENABLE_QUADPREC
    module procedure lalg_vectorSortQP
#endif
  end interface

  public :: lalg_vectorSort
  public :: lalg_vectorSortSP
  public :: lalg_vectorSortDP
  public :: lalg_vectorSortQP
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
    module procedure lalg_vectorAddScalarSP
    module procedure lalg_vectorAddScalarDP
    module procedure lalg_vectorAddScalarSP2D
    module procedure lalg_vectorAddScalarDP2D
    module procedure lalg_vectorAddScalarI8
    module procedure lalg_vectorAddScalarI16
    module procedure lalg_vectorAddScalarI32
    module procedure lalg_vectorAddScalarI64
#ifdef ENABLE_QUADPREC
    module procedure lalg_vectorAddScalarQP
    module procedure lalg_vectorAddScalarQP2D
#endif
  end interface

  public :: lalg_vectorAddScalar
  public :: lalg_vectorAddScalarSP
  public :: lalg_vectorAddScalarDP
  public :: lalg_vectorAddScalarQP
  public :: lalg_vectorAddScalarSP2D
  public :: lalg_vectorAddScalarDP2D
  public :: lalg_vectorAddScalarQP2D
  public :: lalg_vectorAddScalarI8
  public :: lalg_vectorAddScalarI16
  public :: lalg_vectorAddScalarI32
  public :: lalg_vectorAddScalarI64

  interface lalg_vectorCompMult
    module procedure lalg_vectorCompMultSP
    module procedure lalg_vectorCompMultDP
    module procedure lalg_vectorCompMultDPSP
#ifdef ENABLE_QUADPREC
    module procedure lalg_vectorCompMultQP
    module procedure lalg_vectorCompMultQPSP
    module procedure lalg_vectorCompMultQPDP
#endif
  end interface

  public :: lalg_vectorCompMult
  public :: lalg_vectorCompMultSP
  public :: lalg_vectorCompMultDP
  public :: lalg_vectorCompMultQP
  public :: lalg_vectorCompMultDPSP
  public :: lalg_vectorCompMultQPSP
  public :: lalg_vectorCompMultQPDP

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

  !*****************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: lalg_perfconfig

  !*****************************************************************************

contains

  ! ****************************************************************************

!<subroutine>

  subroutine lalg_initPerfConfig(rperfconfig)

!<description>
  ! This routine initialises the global performance configuration
!</description>

!<input>
  ! OPTIONAL: performance configuration that should be used to initialise
  ! the global performance configuration. If not present, the values of
  ! the legacy constants is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
!</subroutine>

    if (present(rperfconfig)) then
      lalg_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(lalg_perfconfig)
    end if

  end subroutine lalg_initPerfConfig

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSP (Fx,Fy,n)

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

  subroutine lalg_copyVectorDP (Dx,Dy,n)

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

  subroutine lalg_copyVectorQP (Qx,Qy,n)

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

    if (.not. present(n)) then
      call QCOPY(size(Qx),Qx,1,Qy,1)
    else
      call QCOPY(n,Qx,1,Qy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSPDP (Fx,Dy,n)

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

    if (.not. present(n)) then
      call SDCOPY(size(Fx),Fx,1,Dy,1)
    else
      call SDCOPY(n,Fx,1,Dy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSPQP (Fx,Qy,n)

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

    if (.not. present(n)) then
      call SQCOPY(size(Fx),Fx,1,Qy,1)
    else
      call SQCOPY(n,Fx,1,Qy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDPSP (Dx,Fy,n)

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

    if (.not. present(n)) then
      call DSCOPY(size(Dx),Dx,1,Fy,1)
    else
      call DSCOPY(n,Dx,1,Fy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDPQP (Dx,Qy,n)

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

    if (.not. present(n)) then
      call DQCOPY(size(Dx),Dx,1,Qy,1)
    else
      call DQCOPY(n,Dx,1,Qy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQPSP (Qx,Fy,n)

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

    if (.not. present(n)) then
      call QSCOPY(size(Qx),Qx,1,Fy,1)
    else
      call QSCOPY(n,Qx,1,Fy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQPDP (Qx,Dy,n)

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

    if (.not. present(n)) then
      call QDCOPY(size(Qx),Qx,1,Dy,1)
    else
      call QDCOPY(n,Qx,1,Dy,1)
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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = Ix(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = Ix(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I16)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = int(Ix(i),I16)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I32)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = int(Ix(i),I32)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I64)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = int(Ix(i),I64)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = Ix(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = Ix(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I8)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
      do i = 1, n
        Iy(i) = int(Ix(i),I8)
      end do
      !$omp end parallel do simd
#else
      do i = 1, n
        Iy(i) = int(Ix(i),I8)
      end do
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I32)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = int(Ix(i),I32)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I64)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = int(Ix(i),I64)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = Ix(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = Ix(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I8)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = int(Ix(i),I8)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I16)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = int(Ix(i),I16)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I64)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = int(Ix(i),I64)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = Ix(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = Ix(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I8)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = int(Ix(i),I8)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I16)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = int(Ix(i),I16)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Iy(i) = int(Ix(i),I32)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Iy(i) = int(Ix(i),I32)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Lx) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Lx)
        Ly(i) = Lx(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Ly(i) = Lx(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Sx) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Sx)
        Sy(i) = Sx(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Sy(i) = Sx(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSP2D (Fx,Fy,n,m)

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

  subroutine lalg_copyVectorDP2D (Dx,Dy,n,m)

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

  subroutine lalg_copyVectorQP2D (Qx,Qy,n,m)

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

    if (present(n) .and. present(m)) then
      call QCOPY(n*m,Qx,1,Qy,1)
    else
      call QCOPY(size(Qx,1)*size(Qx,2),Qx,1,Qy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSPDP2D (Fx,Dy,n,m)

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

    if (present(n) .and. present(m)) then
      call SDCOPY(n*m,Fx,1,Dy,1)
    else
      call SDCOPY(size(Fx,1)*size(Fx,2),Fx,1,Dy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSPQP2D (Fx,Qy,n,m)

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

    if (present(n) .and. present(m)) then
      call SQCOPY(n*m,Fx,1,Qy,1)
    else
      call SQCOPY(size(Fx,1)*size(Fx,2),Fx,1,Qy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDPSP2D (Dx,Fy,n,m)

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

    if (present(n) .and. present(m)) then
      call DSCOPY(n*m,Dx,1,Fy,1)
    else
      call DSCOPY(size(Dx,1)*size(Dx,2),Dx,1,Fy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDPQP2D (Dx,Qy,n,m)

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

    if (present(n) .and. present(m)) then
      call DQCOPY(n*m,Dx,1,Qy,1)
    else
      call DQCOPY(size(Dx,1)*size(Dx,2),Dx,1,Qy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQPSP2D (Qx,Fy,n,m)

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

    if (present(n) .and. present(m)) then
      call QSCOPY(n*m,Qx,1,Fy,1)
    else
      call QSCOPY(size(Qx,1)*size(Qx,2),Qx,1,Fy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQPDP2D (Qx,Dy,n,m)

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

    if (present(n) .and. present(m)) then
      call QDCOPY(n*m,Qx,1,Dy,1)
    else
      call QDCOPY(size(Qx,1)*size(Qx,2),Qx,1,Dy,1)
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

  subroutine lalg_copyVectorSP3D (Fx,Fy,n,m,o)

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

  subroutine lalg_copyVectorDP3D (Dx,Dy,n,m,o)

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

  subroutine lalg_copyVectorQP3D (Qx,Qy,n,m,o)

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

    if (present(n) .and. present(m) .and. present(o)) then
      call QCOPY(n*m*o,Qx,1,Qy,1)
    else
      call QCOPY(size(Qx,1)*size(Qx,2)*size(Qx,3),Qx,1,Qy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSPDP3D (Fx,Dy,n,m,o)

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

    if (present(n) .and. present(m) .and. present(o)) then
      call SDCOPY(n*m*o,Fx,1,Dy,1)
    else
      call SDCOPY(size(Fx,1)*size(Fx,2)*size(Fx,3),Fx,1,Dy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSPQP3D (Fx,Qy,n,m,o)

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

    if (present(n) .and. present(m) .and. present(o)) then
      call SQCOPY(n*m*o,Fx,1,Qy,1)
    else
      call SQCOPY(size(Fx,1)*size(Fx,2)*size(Fx,3),Fx,1,Qy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDPSP3D (Dx,Fy,n,m,o)

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

    if (present(n) .and. present(m) .and. present(o)) then
      call DSCOPY(n*m*o,Dx,1,Fy,1)
    else
      call DSCOPY(size(Dx,1)*size(Dx,2)*size(Dx,3),Dx,1,Fy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDPQP3D (Dx,Qy,n,m,o)

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

    if (present(n) .and. present(m) .and. present(o)) then
      call DQCOPY(n*m*o,Dx,1,Qy,1)
    else
      call DQCOPY(size(Dx,1)*size(Dx,2)*size(Dx,3),Dx,1,Qy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQPSP3D (Qx,Fy,n,m,o)

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

    if (present(n) .and. present(m) .and. present(o)) then
      call QSCOPY(n*m*o,Qx,1,Fy,1)
    else
      call QSCOPY(size(Qx,1)*size(Qx,2)*size(Qx,3),Qx,1,Fy,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorQPDP3D (Qx,Dy,n,m,o)

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

    if (present(n) .and. present(m) .and. present(o)) then
      call QDCOPY(n*m*o,Qx,1,Dy,1)
    else
      call QDCOPY(size(Qx,1)*size(Qx,2)*size(Qx,3),Qx,1,Dy,1)
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

  subroutine lalg_scaleVectorSP (Fx,sc,n)

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
        call lalg_clearVectorSP(Fx)
      else if(sc .ne. 1.0_SP) then
        call SSCAL(size(Fx),sc,Fx,1)
      end if
    else
      if(sc .eq. 0.0_SP) then
        call lalg_clearVectorSP(Fx,n)
      else if(sc .ne. 1.0_SP) then
        call SSCAL(n,sc,Fx,1)
      end if
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_scaleVectorDP (Dx,dc,n)

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
        call lalg_clearVectorDP(Dx)
      else if(dc .ne. 1.0_DP) then
        call DSCAL(size(Dx),dc,Dx,1)
      end if
    else
      if(dc .eq. 0.0_DP) then
        call lalg_clearVectorDP(Dx,n)
      else if(dc .ne. 1.0_DP) then
        call DSCAL(n,dc,Dx,1)
      end if
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_scaleVectorQP (Qx,qc,n)

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

    if (.not. present(n)) then
      if(qc .eq. 0.0_QP) then
        call lalg_clearVectorQP(Qx)
      else if(qc .ne. 1.0_QP) then
        call QSCAL(size(Qx),qc,Qx,1)
      end if
    else
      if(qc .eq. 0.0_QP) then
        call lalg_clearVectorQP(Qx,n)
      else if(qc .ne. 1.0_QP) then
        call QSCAL(n,qc,Qx,1)
      end if
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_scaleVectorSP2D (Fx,sc)

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
      call lalg_clearVectorSP2D(Fx)
    else if(sc .ne. 1.0_SP) then
      call SSCAL(size(Fx,1)*size(Fx,2),sc,Fx,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_scaleVectorDP2D (Dx,dc)

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
      call lalg_clearVectorDP2D(Dx)
    else if(dc .ne. 1.0_DP) then
      call DSCAL(size(Dx,1)*size(Dx,2),dc,Dx,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_scaleVectorQP2D (Qx,qc)

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

    if(qc .eq. 0.0_QP) then
      call lalg_clearVectorQP2D(Qx)
    else if(qc .ne. 1.0_QP) then
      call QSCAL(size(Qx,1)*size(Qx,2),qc,Qx,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_scaleVectorSP3D (Fx,sc)

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
      call lalg_clearVectorSP3D(Fx)
    else if(sc .ne. 1.0_SP) then
      call SSCAL(size(Fx,1)*size(Fx,2)*size(Fx,3),sc,Fx,1)
    end if

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine lalg_scaleVectorDP3D (Dx,dc)

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
      call lalg_clearVectorDP3D(Dx)
    else if(dc .ne. 1.0_DP) then
      call DSCAL(size(Dx,1)*size(Dx,2)*size(Dx,3),dc,Dx,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_scaleVectorQP3D (Qx,qc)

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

    if(qc .eq. 0.0_QP) then
      call lalg_clearVectorQP3D(Qx)
    else if(qc .ne. 1.0_QP) then
      call QSCAL(size(Qx,1)*size(Qx,2)*size(Qx,3),qc,Qx,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorSP (Fx,n)

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

    if (.not. present(n)) then
      call SSET(size(Fx),0.0_SP,Fx,1)
    else
      call SSET(n,0.0_SP,Fx,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorDP (Dx,n)

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

    if (.not. present(n)) then
      call DSET(size(Dx),0.0_DP,Dx,1)
    else
      call DSET(n,0.0_DP,Dx,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorQP (Qx,n)

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

    if (.not. present(n)) then
      call QSET(size(Qx),0.0_QP,Qx,1)
    else
      call QSET(n,0.0_QP,Qx,1)
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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = 0_I8
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = 0_I8
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = 0_I16
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = 0_I16
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = 0_I32
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = 0_I32
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = 0_I64
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = 0_I64
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorSP2D (Fx)

!<description>
  ! Clears a single precision vector: Fx = 0
!</description>

!<output>

  ! Destination vector to be cleared
  real(SP), dimension(:,:), intent(out) :: Fx

!</output>

!</subroutine>

    call SSET(size(Fx,1)*size(Fx,2),0.0_SP,Fx,1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorDP2D (Dx)

!<description>
  ! Clears a double precision vector: Dx = 0
!</description>

!<output>
  ! Destination vector to be cleared
  real(DP), dimension(:,:), intent(out) :: Dx
!</output>

!</subroutine>

    call DSET(size(Dx,1)*size(Dx,2),0.0_DP,Dx,1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorQP2D (Qx)

!<description>
  ! Clears a quad precision vector: Qx = 0
!</description>

!<output>
  ! Destination vector to be cleared
  real(QP), dimension(:,:), intent(out) :: Qx
!</output>

!</subroutine>

    call QSET(size(Qx,1)*size(Qx,2),0.0_QP,Qx,1)

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

  subroutine lalg_clearVectorSP3D (Fx)

!<description>
  ! Clears a single precision vector: Fx = 0
!</description>

!<output>

  ! Destination vector to be cleared
  real(SP), dimension(:,:,:), intent(out) :: Fx

!</output>

!</subroutine>

    call SSET(size(Fx,1)*size(Fx,2)*size(Fx,3),0.0_SP,Fx,1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorDP3D (Dx)

!<description>
  ! Clears a double precision vector: Dx = 0
!</description>

!<output>
  ! Destination vector to be cleared
  real(DP), dimension(:,:,:), intent(out) :: Dx
!</output>

!</subroutine>

    call DSET(size(Dx,1)*size(Dx,2)*size(Dx,3),0.0_DP,Dx,1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorQP3D (Qx)

!<description>
  ! Clears a quad precision vector: Qx = 0
!</description>

!<output>
  ! Destination vector to be cleared
  real(QP), dimension(:,:,:), intent(out) :: Qx
!</output>

!</subroutine>

    call QSET(size(Qx,1)*size(Qx,2)*size(Qx,3),0.0_QP,Qx,1)

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

  subroutine lalg_setVectorSP (Fx,fvalue,n)

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

    if (.not. present(n)) then
      call SSET(size(Fx),fvalue,Fx,1)
    else
      call SSET(n,fvalue,Fx,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorDP (Dx,dvalue,n)

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

    if (.not. present(n)) then
      call DSET(size(Dx),dvalue,Dx,1)
    else
      call DSET(n,dvalue,Dx,1)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorQP (Qx,qvalue,n)

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

    if (.not. present(n)) then
      call QSET(size(Qx),qvalue,Qx,1)
    else
      call QSET(n,qvalue,Qx,1)
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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Ix(i) = ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Ix(i) = ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Ix(i) = ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Ix(i) = ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Lx) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Lx)
        Lx(i) = lvalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Lx(i) = lvalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Sx) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Sx)
        Sx(i) = svalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Sx(i) = svalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorSP2D (Fx,fvalue)

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

    call SSET(size(Fx,1)*size(Fx,2),fvalue,Fx,1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorDP2D (Dx,dvalue)

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

    call DSET(size(Dx,1)*size(Dx,2),dvalue,Dx,1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorQP2D (Qx,qvalue)

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

    call QSET(size(Qx,1)*size(Qx,2),qvalue,Qx,1)

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

  subroutine lalg_setVectorSP3D (Fx,fvalue)

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

    call SSET(size(Fx,1)*size(Fx,2)*size(Fx,3),fvalue,Fx,1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorDP3D (Dx,dvalue)

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

    call DSET(size(Dx,1)*size(Dx,2)*size(Dx,3),dvalue,Dx,1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorQP3D (Qx,qvalue)

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

    call QSET(size(Qx,1)*size(Qx,2)*size(Qx,3),qvalue,Qx,1)

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

  subroutine lalg_vectorLinearCombSP (Fx,Fy,scx,scy,n)

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
        if (scx .ne. 1.0_SP) call SSCAL(size(Fy),scx,Fy,1)
      else if (scy .eq. 1.0_SP) then
        call SAXPY(size(Fx),scx,Fx,1,Fy,1)
      else
        c=scx/scy
        call SAXPY(size(Fx),c,Fx,1,Fy,1)
        call SSCAL(size(Fy),scy,Fy,1)
      endif

    else

      if (scy .eq. 0.0_SP) then
        call SCOPY(n,Fx,1,Fy,1)
        if (scx .ne. 1.0_SP) call SSCAL(n,scx,Fy,1)
      else if (scy .eq. 1.0_SP) then
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

  subroutine lalg_vectorLinearCombDP (Dx,Dy,dcx,dcy,n)

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
        if (dcx .ne. 1.0_DP) call DSCAL(size(Dy),dcx,Dy,1)
      else if (dcy .eq. 1.0_DP) then
        call DAXPY(size(Dx),dcx,Dx,1,Dy,1)
      else
        c=dcx/dcy
        call DAXPY(size(Dx),c,Dx,1,Dy,1)
        call DSCAL(size(Dy),dcy,Dy,1)
      endif

    else

      if (dcy .eq. 0.0_DP) then
        call DCOPY(n,Dx,1,Dy,1)
        if (dcx .ne. 1.0_DP) call DSCAL(n,dcx,Dy,1)
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

  subroutine lalg_vectorLinearCombQP (Qx,Qy,qcx,qcy,n)

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
  real(QP) :: c

    if (.not. present(n)) then

      if (qcy .eq. 0.0_QP) then
        call QCOPY(size(Qx),Qx,1,Qy,1)
        if (qcx .ne. 1.0_QP) call QSCAL(size(Qy),qcx,Qy,1)
      else if (qcy .eq. 1.0_QP) then
        call QAXPY(size(Qx),qcx,Qx,1,Qy,1)
      else
        c=qcx/qcy
        call QAXPY(size(Qx),c,Qx,1,Qy,1)
        call QSCAL(size(Qy),qcy,Qy,1)
      endif

    else

      if (qcy .eq. 0.0_QP) then
        call QCOPY(n,Qx,1,Qy,1)
        if (qcx .ne. 1.0_QP) call QSCAL(n,qcx,Qy,1)
      else if (qcy .eq. 1.0_QP) then
        call QAXPY(n,qcx,Qx,1,Qy,1)
      else
        c=qcx/qcy
        call QAXPY(n,c,Qx,1,Qy,1)
        call QSCAL(n,qcy,Qy,1)
      endif

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSPDP (Fx,Dy,scx,dcy,n)

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
  real(SP) :: c

    if (.not. present(n)) then

      if (dcy .eq. 0.0_DP) then
        call SDCOPY(size(Fx),Fx,1,Dy,1)
        if (scx .ne. 1.0_SP) call DSCAL(size(Dy),real(scx,DP),Dy,1)
      else if (dcy .eq. 1.0_DP) then
        call SDAXPY(size(Fx),scx,Fx,1,Dy,1)
      else
        c=real(scx/dcy,SP)
        call SDAXPY(size(Fx),c,Fx,1,Dy,1)
        call DSCAL(size(Dy),dcy,Dy,1)
      endif

    else

      if (dcy .eq. 0.0_DP) then
        call SDCOPY(n,Fx,1,Dy,1)
        if (scx .ne. 1.0_SP) call DSCAL(size(Dy),real(scx,DP),Dy,1)
      else if (dcy .eq. 1.0_DP) then
        call SDAXPY(n,scx,Fx,1,Dy,1)
      else
        c=real(scx/dcy,SP)
        call SDAXPY(n,c,Fx,1,Dy,1)
        call DSCAL(n,dcy,Dy,1)
      endif

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSPQP (Fx,Qy,scx,qcy,n)

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
  real(SP) :: c

    if (.not. present(n)) then

      if (qcy .eq. 0.0_QP) then
        call SQCOPY(size(Fx),Fx,1,Qy,1)
        if (scx .ne. 1.0_SP) call QSCAL(size(Qy),real(scx,QP),Qy,1)
      else if (qcy .eq. 1.0_QP) then
        call SQAXPY(size(Fx),scx,Fx,1,Qy,1)
      else
        c=real(scx/qcy,SP)
        call SQAXPY(size(Fx),c,Fx,1,Qy,1)
        call QSCAL(size(Qy),qcy,Qy,1)
      endif

    else

      if (qcy .eq. 0.0_QP) then
        call SQCOPY(n,Fx,1,Qy,1)
        if (scx .ne. 1.0_SP) call QSCAL(size(Qy),real(scx,QP),Qy,1)
      else if (qcy .eq. 1.0_QP) then
        call SQAXPY(n,scx,Fx,1,Qy,1)
      else
        c=real(scx/qcy,SP)
        call SQAXPY(n,c,Fx,1,Qy,1)
        call QSCAL(n,qcy,Qy,1)
      endif

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombDPQP (Dx,Qy,dcx,qcy,n)

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
  real(DP) :: c

    if (.not. present(n)) then

      if (qcy .eq. 0.0_QP) then
        call DQCOPY(size(Dx),Dx,1,Qy,1)
        if (dcx .ne. 1.0_DP) call QSCAL(size(Qy),real(dcx,QP),Qy,1)
      else if (qcy .eq. 1.0_QP) then
        call DQAXPY(size(Dx),dcx,Dx,1,Qy,1)
      else
        c=real(dcx/qcy,DP)
        call DQAXPY(size(Dx),c,Dx,1,Qy,1)
        call QSCAL(size(Qy),qcy,Qy,1)
      endif

    else

      if (qcy .eq. 0.0_QP) then
        call DQCOPY(n,Dx,1,Qy,1)
        if (dcx .ne. 1.0_DP) call QSCAL(size(Qy),real(dcx,QP),Qy,1)
      else if (qcy .eq. 1.0_QP) then
        call DQAXPY(n,dcx,Dx,1,Qy,1)
      else
        c=real(dcx/qcy,DP)
        call DQAXPY(n,c,Dx,1,Qy,1)
        call QSCAL(n,qcy,Qy,1)
      endif

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSP2D (Fx,Fy,scx,scy)

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
  integer :: n

  n = size(Fx,1)*size(Fx,2)

    if (scy .eq. 0.0_SP) then
      call SCOPY(n,Fx,1,Fy,1)
      if (scx .ne. 1.0_SP) call SSCAL(n,scx,Fy,1)
    else if (scy .eq. 1.0_SP) then
      call SAXPY(n,scx,Fx,1,Fy,1)
    else
      c=scx/scy
      call SAXPY(n,c,Fx,1,Fy,1)
      call SSCAL(n,scy,Fy,1)
    endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombDP2D (Dx,Dy,dcx,dcy)

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
  integer :: n

  n = size(Dx,1)*size(Dx,2)

    if (dcy .eq. 0.0_DP) then
      call DCOPY(n,Dx,1,Dy,1)
      if (dcx .ne. 1.0_DP) call DSCAL(n,dcx,Dy,1)
    else if (dcy .eq. 1.0_DP) then
      call DAXPY(n,dcx,Dx,1,Dy,1)
    else
      c=dcx/dcy
      call DAXPY(n,c,Dx,1,Dy,1)
      call DSCAL(n,dcy,Dy,1)
    endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombQP2D (Qx,Qy,qcx,qcy)

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
  real(QP) :: c
  integer :: n

  n = size(Qx,1)*size(Qx,2)

    if (qcy .eq. 0.0_QP) then
      call QCOPY(n,Qx,1,Qy,1)
      if (qcx .ne. 1.0_QP) call QSCAL(n,qcx,Qy,1)
    else if (qcy .eq. 1.0_QP) then
      call QAXPY(n,qcx,Qx,1,Qy,1)
    else
      c=qcx/qcy
      call QAXPY(n,c,Qx,1,Qy,1)
      call QSCAL(n,qcy,Qy,1)
    endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSPDP2D (Fx,Dy,scx,dcy)

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

  ! local variables
  real(SP) :: c
  integer :: n

  n = size(Fx,1)*size(Fx,2)

    if (dcy .eq. 0.0_DP) then
      call SDCOPY(n,Fx,1,Dy,1)
      if (scx .ne. 1.0_SP) call DSCAL(n,real(scx,DP),Dy,1)
    else if (dcy .eq. 1.0_DP) then
      call SDAXPY(n,scx,Fx,1,Dy,1)
    else
      c=real(scx/dcy,SP)
      call SDAXPY(n,c,Fx,1,Dy,1)
      call DSCAL(n,dcy,Dy,1)
    endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSPQP2D (Fx,Qy,scx,qcy)

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

  ! local variables
  real(SP) :: c
  integer :: n

  n = size(Fx,1)*size(Fx,2)

    if (qcy .eq. 0.0_QP) then
      call SQCOPY(n,Fx,1,Qy,1)
      if (scx .ne. 1.0_SP) call QSCAL(n,real(scx,QP),Qy,1)
    else if (qcy .eq. 1.0_QP) then
      call SQAXPY(n,scx,Fx,1,Qy,1)
    else
      c=real(scx/qcy,SP)
      call SQAXPY(n,c,Fx,1,Qy,1)
      call QSCAL(n,qcy,Qy,1)
    endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombDPQP2D (Dx,Qy,dcx,qcy)

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

  ! local variables
  real(DP) :: c
  integer :: n

  n = size(Dx,1)*size(Dx,2)

    if (qcy .eq. 0.0_QP) then
      call DQCOPY(n,Dx,1,Qy,1)
      if (dcx .ne. 1.0_DP) call QSCAL(n,real(dcx,QP),Qy,1)
    else if (qcy .eq. 1.0_QP) then
      call DQAXPY(n,dcx,Dx,1,Qy,1)
    else
      c=real(dcx/qcy,DP)
      call DQAXPY(n,c,Dx,1,Qy,1)
      call QSCAL(n,qcy,Qy,1)
    endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSP3D (Fx,Fy,scx,scy)

!<description>
  ! Performs a linear combination: Fy = scx * Fx  +  scy * Fy
!</description>

!<input>

  ! First source vector
  real(SP), dimension(:,:,:), intent(in) :: Fx

  ! Scaling factor for Dx
  real(SP), intent(in)                   :: scx

  ! Scaling factor for Dy
  real(SP), intent(in)                   :: scy

!</input>

!<inputoutput>

  ! Second source vector; also receives the result
  real(SP), dimension(:,:,:), intent(inout) :: Fy

!</inputoutput>

!</subroutine>

  ! local variables
  real(SP) :: c
  integer :: n

  n = size(Fx,1)*size(Fx,2)*size(Fx,3)

    if (scy .eq. 0.0_SP) then
      call SCOPY(n,Fx,1,Fy,1)
      if (scx .ne. 1.0_SP) call SSCAL(n,scx,Fy,1)
    else if (scy .eq. 1.0_SP) then
      call SAXPY(n,scx,Fx,1,Fy,1)
    else
      c=scx/scy
      call SAXPY(n,c,Fx,1,Fy,1)
      call SSCAL(n,scy,Fy,1)
    endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombDP3D (Dx,Dy,dcx,dcy)

!<description>
  ! Performs a linear combination: Dy = dcx * Dx  +  dcy * Dy
!</description>

!<input>

  ! First source vector
  real(DP), dimension(:,:,:), intent(in) :: Dx

  ! Scaling factor for Dx
  real(DP), intent(in)                   :: dcx

  ! Scaling factor for Dy
  real(DP), intent(in)                   :: dcy

!</input>

!<inputoutput>

  ! Second source vector; also receives the result
  real(DP), dimension(:,:,:), intent(inout) :: Dy

!</inputoutput>

!</subroutine>

  ! local variables
  real(DP) :: c
  integer :: n

  n = size(Dx,1)*size(Dx,2)*size(Dx,3)

    if (dcy .eq. 0.0_DP) then
      call DCOPY(n,Dx,1,Dy,1)
      if (dcx .ne. 1.0_DP) call DSCAL(n,dcx,Dy,1)
    else if (dcy .eq. 1.0_DP) then
      call DAXPY(n,dcx,Dx,1,Dy,1)
    else
      c=dcx/dcy
      call DAXPY(n,c,Dx,1,Dy,1)
      call DSCAL(n,dcy,Dy,1)
    endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombQP3D (Qx,Qy,qcx,qcy)

!<description>
  ! Performs a linear combination: Qy = qcx * Qx  +  qcy * Qy
!</description>

!<input>

  ! First source vector
  real(QP), dimension(:,:,:), intent(in) :: Qx

  ! Scaling factor for Qx
  real(QP), intent(in)                   :: qcx

  ! Scaling factor for Qy
  real(QP), intent(in)                   :: qcy

!</input>

!<inputoutput>

  ! Second source vector; also receives the result
  real(QP), dimension(:,:,:), intent(inout) :: Qy

!</inputoutput>

!</subroutine>

  ! local variables
  real(QP) :: c
  integer :: n

  n = size(Qx,1)*size(Qx,2)*size(Qx,3)

    if (qcy .eq. 0.0_QP) then
      call QCOPY(n,Qx,1,Qy,1)
      if (qcx .ne. 1.0_QP) call QSCAL(n,qcx,Qy,1)
    else if (qcy .eq. 1.0_QP) then
      call QAXPY(n,qcx,Qx,1,Qy,1)
    else
      c=qcx/qcy
      call QAXPY(n,c,Qx,1,Qy,1)
      call QSCAL(n,qcy,Qy,1)
    endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSPDP3D (Fx,Dy,scx,dcy)

!<description>
  ! Performs a linear combination: Dy = scx * Fx  +  dcy * Dy
!</description>

!<input>

  ! First source vector
  real(SP), dimension(:,:,:), intent(in) :: Fx

  ! Scaling factor for Dx
  real(SP), intent(in)                   :: scx

  ! Scaling factor for Dy
  real(DP), intent(in)                   :: dcy

!</input>

!<inputoutput>

  ! Second source vector; also receives the result
  real(DP), dimension(:,:,:), intent(inout) :: Dy

!</inputoutput>

!</subroutine>

  ! local variables
  real(SP) :: c
  integer :: n

  n = size(Fx,1)*size(Fx,2)*size(Fx,3)

    if (dcy .eq. 0.0_DP) then
      call SDCOPY(n,Fx,1,Dy,1)
      if (scx .ne. 1.0_SP) call DSCAL(n,real(scx,DP),Dy,1)
    else if (dcy .eq. 1.0_DP) then
      call SDAXPY(n,scx,Fx,1,Dy,1)
    else
      c=real(scx/dcy,SP)
      call SDAXPY(n,c,Fx,1,Dy,1)
      call DSCAL(n,dcy,Dy,1)
    endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombSPQP3D (Fx,Qy,scx,qcy)

!<description>
  ! Performs a linear combination: Qy = scx * Fx  +  qcy * Qy
!</description>

!<input>

  ! First source vector
  real(SP), dimension(:,:,:), intent(in) :: Fx

  ! Scaling factor for Dx
  real(SP), intent(in)                   :: scx

  ! Scaling factor for Qy
  real(QP), intent(in)                   :: qcy

!</input>

!<inputoutput>

  ! Second source vector; also receives the result
  real(QP), dimension(:,:,:), intent(inout) :: Qy

!</inputoutput>

!</subroutine>

  ! local variables
  real(SP) :: c
  integer :: n

  n = size(Fx,1)*size(Fx,2)*size(Fx,3)

    if (qcy .eq. 0.0_QP) then
      call SQCOPY(n,Fx,1,Qy,1)
      if (scx .ne. 1.0_SP) call QSCAL(n,real(scx,QP),Qy,1)
    else if (qcy .eq. 1.0_QP) then
      call SQAXPY(n,scx,Fx,1,Qy,1)
    else
      c=real(scx/qcy,SP)
      call SQAXPY(n,c,Fx,1,Qy,1)
      call QSCAL(n,qcy,Qy,1)
    endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorLinearCombDPQP3D (Dx,Qy,dcx,qcy)

!<description>
  ! Performs a linear combination: Qy = dcx * Dx  +  qcy * Qy
!</description>

!<input>

  ! First source vector
  real(DP), dimension(:,:,:), intent(in) :: Dx

  ! Scaling factor for Dx
  real(DP), intent(in)                   :: dcx

  ! Scaling factor for Qy
  real(QP), intent(in)                   :: qcy

!</input>

!<inputoutput>

  ! Second source vector; also receives the result
  real(QP), dimension(:,:,:), intent(inout) :: Qy

!</inputoutput>

!</subroutine>

  ! local variables
  real(DP) :: c
  integer :: n

  n = size(Dx,1)*size(Dx,2)*size(Dx,3)

    if (qcy .eq. 0.0_QP) then
      call DQCOPY(n,Dx,1,Qy,1)
      if (dcx .ne. 1.0_DP) call QSCAL(n,real(dcx,QP),Qy,1)
    else if (qcy .eq. 1.0_QP) then
      call DQAXPY(n,dcx,Dx,1,Qy,1)
    else
      c=real(dcx/qcy,DP)
      call DQAXPY(n,c,Dx,1,Qy,1)
      call QSCAL(n,qcy,Qy,1)
    endif

  end subroutine

  ! ***************************************************************************

!<function>

  real(SP) function lalg_scalarProductSP (Fx,Fy,n) result (res)

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

  real(DP) function lalg_scalarProductDP (Dx,Dy,n) result (res)

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

  real(QP) function lalg_scalarProductQP (Qx,Qy,n) result (res)

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Qx) > lalg_perfconfig%NITEMMIN_OMP) reduction(+:res)
#endif
      do i = 1, size(Qx)
        res = res + Qx(i)*Qy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

      !res = QDOT(size(Qx),Qx,1,Qy,1)
    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP) reduction(+:res)
#endif
      do i = 1, n
        res = res + Qx(i)*Qy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

      !res = QDOT(n,Qx,1,Qy,1)
    end if

  end function

  ! ***************************************************************************

!<function>

  real(SP) function lalg_scalarProductSP2D (Fx,Fy) result (res)

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

  real(DP) function lalg_scalarProductDP2D (Dx,Dy) result (res)

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

  real(DP) function lalg_scalarProductQP2D (Qx,Qy) result (res)

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

  real(SP) function lalg_normSP (Fx,cnorm,iposMax,n) result(resnorm)

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

  real(DP) function lalg_normDP (Dx,cnorm,iposMax,n) result(resnorm)

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

  real(QP) function lalg_normQP (Qx,cnorm,iposMax,n) result(resnorm)

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

!  real(QP) :: QASUM,QNRM2
!  integer :: IQAMAX

  integer :: isize,iposmaxlocal

    isize = size(Qx)
    if (present(n)) isize = n

    ! Choose the norm to calculate
    select case (cnorm)
    case (LINALG_NORMSUM)
      ! L1-norm: sum all entries
      resnorm = QASUM(isize,Qx,1)

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      resnorm = QNRM2(isize,Qx,1)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vector (1111...) has norm = 1.
      resnorm = QASUM(isize,Qx,1) / real(isize,QP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vector (1111...) has norm = 1.
      resnorm = QNRM2(isize,Qx,1) / sqrt(real(isize,QP))

    case (LINALG_NORMMAX)
      ! MAX-norm. Find the absolute largest entry.
      ! With the BLAS routine, calculate the position. Then get the entry.
      iposmaxlocal = IQAMAX(isize,Qx,1)
      resnorm = abs(Qx(iposmaxlocal))

      if(present(iposMax)) then
        iposMax = iposmaxlocal
      end if

      case default
      ! Unknown norm
      resnorm = -1.0_QP

    end select

  end function

  ! ***************************************************************************

!<function>

  real(SP) function lalg_errorNormSP (Fx,Fy,cnorm,iposMax,n,Fw) result(resnorm)

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

  ! OPTIONAL: Weighting vector
  real(DP), dimension(:), intent(in), optional :: Fw
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
      if (present(Fw)) then
#ifdef HAS_OPENMP40
        !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
        do i = 1, isize
          resnorm = resnorm + Fw(i)*abs(Fx(i)-Fy(i))
        end do
#ifdef HAS_OPENMP40
        !$omp end parallel do simd
#endif
      else
#ifdef HAS_OPENMP40
        !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
        do i = 1, isize
          resnorm = resnorm + abs(Fx(i)-Fy(i))
        end do
#ifdef HAS_OPENMP40
        !$omp end parallel do simd
#endif
      end if

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      if (present(Fw)) then
#ifdef HAS_OPENMP40
        !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
        do i = 1, isize
          resnorm = resnorm + Fw(i)*(Fx(i)-Fy(i))*(Fx(i)-Fy(i))
        end do
#ifdef HAS_OPENMP40
        !$omp end parallel do simd
#endif
      else
#ifdef HAS_OPENMP40
        !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
        do i = 1, isize
          resnorm = resnorm + (Fx(i)-Fy(i))*(Fx(i)-Fy(i))
        end do
#ifdef HAS_OPENMP40
        !$omp end parallel do simd
#endif
      end if
      resnorm = sqrt(resnorm)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vector (1111...) to has norm = 1.
#ifdef HAS_OPENMP40
      !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
      do i = 1, isize
        resnorm = resnorm + abs(Fx(i)-Fy(i))
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif
      resnorm = resnorm / real(isize,SP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vector (1111...) has norm = 1.
#ifdef HAS_OPENMP40
      !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
      do i = 1, isize
        resnorm = resnorm + (Fx(i)-Fy(i))*(Fx(i)-Fy(i))
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif
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

  real(DP) function lalg_errorNormDP (Dx,Dy,cnorm,iposMax,n,Dw) result(resnorm)

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
#ifdef HAS_OPENMP40
        !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
        do i = 1, isize
          resnorm = resnorm + Dw(i)*abs(Dx(i)-Dy(i))
        end do
#ifdef HAS_OPENMP40
        !$omp end parallel do simd
#endif
      else
#ifdef HAS_OPENMP40
        !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
        do i = 1, isize
          resnorm = resnorm + abs(Dx(i)-Dy(i))
        end do
#ifdef HAS_OPENMP40
        !$omp end parallel do simd
#endif
      end if

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      if (present(Dw)) then
#ifdef HAS_OPENMP40
        !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
        do i = 1, isize
          resnorm = resnorm + Dw(i)*(Dx(i)-Dy(i))*(Dx(i)-Dy(i))
        end do
#ifdef HAS_OPENMP40
        !$omp end parallel do simd
#endif
      else
#ifdef HAS_OPENMP40
        !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
        do i = 1, isize
          resnorm = resnorm + (Dx(i)-Dy(i))*(Dx(i)-Dy(i))
        end do
#ifdef HAS_OPENMP40
        !$omp end parallel do simd
#endif
      end if
      resnorm = sqrt(resnorm)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vector (1111...) has norm = 1.
#ifdef HAS_OPENMP40
      !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
      do i = 1, isize
        resnorm = resnorm + abs(Dx(i)-Dy(i))
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif
      resnorm = resnorm / real(isize,DP)

    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vector (1111...) has norm = 1.
#ifdef HAS_OPENMP40
      !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
      do i = 1, isize
        resnorm = resnorm + (Dx(i)-Dy(i))*(Dx(i)-Dy(i))
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif
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

  real(QP) function lalg_errorNormQP (Qx,Qy,cnorm,iposMax,n,Qw) result(resnorm)

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
#ifdef HAS_OPENMP40
        !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
        do i = 1, isize
          resnorm = resnorm + Qw(i)*abs(Qx(i)-Qy(i))
        end do
#ifdef HAS_OPENMP40
        !$omp end parallel do simd
#endif
      else
#ifdef HAS_OPENMP40
        !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
        do i = 1, isize
          resnorm = resnorm + abs(Qx(i)-Qy(i))
        end do
#ifdef HAS_OPENMP40
        !$omp end parallel do simd
#endif
      end if

    case (LINALG_NORMEUCLID)
      ! Euclidian norm = scalar product (vector,vector)
      if (present(Qw)) then
#ifdef HAS_OPENMP40
        !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
        do i = 1, isize
          resnorm = resnorm + Qw(i)*(Qx(i)-Qy(i))*(Qx(i)-Qy(i))
        end do
#ifdef HAS_OPENMP40
        !$omp end parallel do simd
#endif
      else
#ifdef HAS_OPENMP40
        !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
        do i = 1, isize
          resnorm = resnorm + (Qx(i)-Qy(i))*(Qx(i)-Qy(i))
        end do
#ifdef HAS_OPENMP40
        !$omp end parallel do simd
#endif
      end if
      resnorm = sqrt(resnorm)

    case (LINALG_NORML1)
      ! L1-norm: sum all entries, divide by sqrt(vector length).
      ! So, scale such that the vector (1111...) has norm = 1.
#ifdef HAS_OPENMP40
      !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
      do i = 1, isize
        resnorm = resnorm + abs(Qx(i)-Qy(i))
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif
      resnorm = resnorm / real(isize,DP)
      
    case (LINALG_NORML2)
      ! l_2-norm - like euclidian norm, but divide by vector length.
      ! So, scale such that the vector (1111...) has norm = 1.
#ifdef HAS_OPENMP40
      !$omp parallel do simd if(isize > lalg_perfconfig%NITEMMIN_OMP) reduction(+:resnorm)
#endif
      do i = 1, isize
        resnorm = resnorm + (Qx(i)-Qy(i))*(Qx(i)-Qy(i))
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif
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

  subroutine lalg_vectorSortSP (Fx, Fd, Itr)

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

#ifdef HAS_OPENMP40
    !$omp parallel do simd if(size(Itr) > lalg_perfconfig%NITEMMIN_OMP)
#endif
    do i = 1, size(Itr)
      Fd(i) = Fx(Itr(i))
    end do
#ifdef HAS_OPENMP40
    !$omp end parallel do simd
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorSortDP (Dx, Dd, Itr)

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

#ifdef HAS_OPENMP40
    !$omp parallel do simd if(size(Itr) > lalg_perfconfig%NITEMMIN_OMP)
#endif
    do i = 1, size(Itr)
      Dd(i) = Dx(Itr(i))
    end do
#ifdef HAS_OPENMP40
    !$omp end parallel do simd
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorSortQP (Qx, Qd, Itr)

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

#ifdef HAS_OPENMP40
    !$omp parallel do simd if(size(Itr) > lalg_perfconfig%NITEMMIN_OMP)
#endif
    do i = 1, size(Itr)
      Qd(i) = Qx(Itr(i))
    end do
#ifdef HAS_OPENMP40
    !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
    !$omp parallel do simd if(size(Itr) > lalg_perfconfig%NITEMMIN_OMP)
#endif
    do i = 1, size(Itr)
      Id(i) = Ix(Itr(i))
    end do
#ifdef HAS_OPENMP40
    !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
    !$omp parallel do simd if(size(Itr) > lalg_perfconfig%NITEMMIN_OMP)
#endif
    do i = 1, size(Itr)
      Id(i) = Ix(Itr(i))
    end do
#ifdef HAS_OPENMP40
    !$omp end parallel do simd
#endif

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_tensorProductSP(Sx,Sy,Stensor)

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

  subroutine lalg_tensorProductDP(Dx,Dy,Dtensor)

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

  subroutine lalg_tensorProductQP(Qx,Qy,Qtensor)

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

  subroutine lalg_vectorAddScalarSP (Fx,fvalue,n)

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Fx) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Fx)
        Fx(i) = Fx(i) + fvalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Fx(i) = Fx(i) + fvalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarDP (Dx,dvalue,n)

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Dx) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Dx)
        Dx(i) = Dx(i) + dvalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Dx(i) = Dx(i) + dvalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine


  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarQP (Qx,qvalue,n)

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Qx) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Qx)
        Qx(i) = Qx(i) + qvalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Qx(i) = Qx(i) + qvalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarSP2D (Fx,fvalue)

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

  subroutine lalg_vectorAddScalarDP2D (Dx,dvalue)

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

  subroutine lalg_vectorAddScalarQP2D (Qx,qvalue)

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = Ix(i) + ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Ix(i) = Ix(i) + ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = Ix(i) + ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Ix(i) = Ix(i) + ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = Ix(i) + ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Ix(i) = Ix(i) + ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Ix) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Ix)
        Ix(i) = Ix(i) + ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Ix(i) = Ix(i) + ivalue
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultSP (Fx,Fy,sc,n)

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Fy) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Fy,1)
        Fy(i) = sc*Fx(i)*Fy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Fy(i) = sc*Fx(i)*Fy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultDP (Dx,Dy,dc,n)

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Dy) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Dy,1)
        Dy(i) = dc*Dx(i)*Dy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Dy(i) = dc*Dx(i)*Dy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultQP (Qx,Qy,qc,n)

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Qy) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Qy,1)
        Qy(i) = qc*Qx(i)*Qy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Qy(i) = qc*Qx(i)*Qy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultDPSP (Fx,Dy,dc,n)

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Dy) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Dy,1)
        Dy(i) = dc*Fx(i)*Dy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Dy(i) = dc*Fx(i)*Dy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultQPSP (Fx,Qy,Qc,n)

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Qy) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Qy,1)
        Qy(i) = qc*real(Fx(i),QP)*Qy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Qy(i) = qc*real(Fx(i),QP)*Qy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultQPDP (Dx,Qy,Qc,n)

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

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(size(Qy) > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, size(Qy,1)
        Qy(i) = qc*real(Dx(i),QP)*Qy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    else

#ifdef HAS_OPENMP40
      !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
      do i = 1, n
        Qy(i) = qc*real(Dx(i),QP)*Qy(i)
      end do
#ifdef HAS_OPENMP40
      !$omp end parallel do simd
#endif

    end if

  end subroutine

!*******************************************************************************
! ! !               PRIVATE SUBROUTINE FOR INTERNAL USE ONLY               ! ! !
!*******************************************************************************

!<subroutine>
  subroutine qscal(n,qa,qx,incx)

!<description>
  ! Scales a quad precision vector by a constant qx := qa * qx
  !
  ! REMARK: This subroutine is a port of subroutine DSCAL from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<inputoutput>

  ! Source and destination vector
  real(QP), dimension(*), intent(inout) :: qx

!</inputoutput>

!<input>

  ! Length of the vector
  integer, intent(in) :: n

  ! Increment
  integer, intent(in) :: incx

  ! Multiplication factor
  real(QP), intent(in) :: qa

!</input>
!</subroutine>

  ! local variables
  integer :: i,m,mp1,nincx

  if( n.le.0 .or. incx.le.0 )return
  if(incx.eq.1)go to 20

  ! code for increment not equal to 1
  nincx = n*incx
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = 1,nincx,incx
    qx(i) = qa*qx(i)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  return

  ! code for increment equal to 1
20 m = mod(n,5)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    qx(i) = qa*qx(i)
  end do
  if( n .lt. 5 ) return
40 mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = mp1,n,5
    qx(i) = qa*qx(i)
    qx(i + 1) = qa*qx(i + 1)
    qx(i + 2) = qa*qx(i + 2)
    qx(i + 3) = qa*qx(i + 3)
    qx(i + 4) = qa*qx(i + 4)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end subroutine qscal

! ***************************************************************************

!<subroutine>
  subroutine qcopy(n,qx,incx,qy,incy)

!<description>
  ! Copies a quad precision vector to another vector qy := qx
  !
  ! REMARK: This subroutine is a port of subroutine DCOPY from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Destination vector
  real(QP), dimension(*), intent(out) :: qy

!</output>

!<input>

  ! Source vector
  real(QP), dimension(*), intent(in) :: qx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increments
  integer, intent(in) :: incx,incy

!</input>
!</subroutine>

  ! local variables
  integer :: i,ix,iy,m,mp1

  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20

  ! code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
    qy(iy) = qx(ix)
    ix = ix + incx
    iy = iy + incy
  end do
  return

  ! code for both increments equal to 1
20 m = mod(n,7)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    qy(i) = qx(i)
  end do
  if( n .lt. 7 ) return
40 mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = mp1,n,7
    qy(i) = qx(i)
    qy(i + 1) = qx(i + 1)
    qy(i + 2) = qx(i + 2)
    qy(i + 3) = qx(i + 3)
    qy(i + 4) = qx(i + 4)
    qy(i + 5) = qx(i + 5)
    qy(i + 6) = qx(i + 6)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end subroutine qcopy

! ***************************************************************************

!<subroutine>
  subroutine qaxpy(n,qa,qx,incx,qy,incy)

!<description>
  ! Computes constant times a quad precision vector plus a vector.
  !
  ! REMARK: This subroutine is a port of subroutine DAXPY from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Second source and destination vector
  real(QP), dimension(*), intent(out) :: qy

!</output>

!<input>

  ! First source vector
  real(QP), dimension(*), intent(in) :: qx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increments
  integer, intent(in) :: incx,incy

  ! Multiplication factor
  real(QP), intent(in) :: qa

!</input>
!</subroutine>

  ! local variables
  integer :: i,ix,iy,m,mp1

  if(n.le.0)return
  if (qa .eq. 0.0_QP) return
  if(incx.eq.1.and.incy.eq.1)go to 20

  ! code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
    qy(iy) = qy(iy) + qa*qx(ix)
    ix = ix + incx
    iy = iy + incy
  end do
  return

  ! code for both increments equal to 1
20 m = mod(n,4)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    qy(i) = qy(i) + qa*qx(i)
  end do
  if( n .lt. 4 ) return
40 mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = mp1,n,4
    qy(i) = qy(i) + qa*qx(i)
    qy(i + 1) = qy(i + 1) + qa*qx(i + 1)
    qy(i + 2) = qy(i + 2) + qa*qx(i + 2)
    qy(i + 3) = qy(i + 3) + qa*qx(i + 3)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end subroutine qaxpy

! ***************************************************************************

!<subroutine>
  function qasum(n,qx,incx)

!<description>
  ! Takes the sum of the absolute values.
  !
  ! REMARK: This subroutine is a port of subroutine DASUM from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<input>

  ! First source vector
  real(QP), dimension(*), intent(in) :: qx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increment
  integer, intent(in) :: incx

!</input>

!<result>
  ! Sum of the absolute values
  real(QP) :: qasum
!</result>
!</function>

  ! local variables
  real(QP) :: qtemp
  integer i,m,mp1,nincx

  qasum = 0.0_QP
  qtemp = 0.0_QP
  if( n.le.0 .or. incx.le.0 )return
  if(incx.eq.1)go to 20

  ! code for increment not equal to 1
  nincx = n*incx
  do i = 1,nincx,incx
    qtemp = qtemp + abs(qx(i))
  end do
  qasum = qtemp
  return

  ! code for increment equal to 1
20 m = mod(n,6)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    qtemp = qtemp + abs(qx(i))
  end do
  if( n .lt. 6 ) go to 60
40 mp1 = m + 1
  do i = mp1,n,6
    qtemp = qtemp + abs(qx(i)) + abs(qx(i + 1)) + abs(qx(i + 2))&
          + abs(qx(i + 3)) + abs(qx(i + 4)) + abs(qx(i + 5))
  end do
60 qasum = qtemp
  return

  end function qasum

! ***************************************************************************

!<subroutine>
  function qnrm2(n,qx,incx)

!<description>
  ! Computes the euclidean norm of quad precision vector
  !
  ! REMARK: This subroutine is a port of subroutine DNRM2 from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<input>

  ! First source vector
  real(QP), dimension(*), intent(in) :: qx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increment
  integer, intent(in) :: incx

!</input>

!<result>
  ! Euclidean norm
  real(QP) :: qnrm2
!</result>
!</function>

  ! local variables
  real(QP), parameter :: ONE = 1.0_QP, ZERO = 0.0_QP
  real(QP) :: absxi, norm, scale, ssq
  integer :: ix

  if( n.lt.1 .or. incx.lt.1 )then
    norm  = zero
  else if( n.eq.1 )then
    norm  = abs( qx( 1 ) )
  else
    scale = zero
    ssq   = one
    do ix = 1, 1 + ( n - 1 )*incx, incx
      if( qx( ix ).ne.zero )then
        absxi = abs( qx( ix ) )
        if( scale.lt.absxi )then
          ssq   = one   + ssq*( scale/absxi )**2
          scale = absxi
        else
          ssq   = ssq   +     ( absxi/scale )**2
        end if
      end if
    end do
    norm  = scale * sqrt( ssq )
  end if

  qnrm2 = norm

  end function qnrm2

! ***************************************************************************

!<subroutine>
  function iqamax(n,qx,incx)

!<description>
  ! Finds the index of element having max. absolute value.
  !
  ! REMARK: This subroutine is a port of subroutine IDAMAX from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<input>

  ! First source vector
  real(QP), dimension(*), intent(in) :: qx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increments
  integer, intent(in) :: incx

!</input>

!<result>
  ! Index of element having max. absolute value
  integer :: iqamax
!</result>
!</function>

  ! local variables
  real(QP) :: qmax
  integer :: i,ix

  iqamax = 0
  if( n.lt.1 .or. incx.le.0 ) return
  iqamax = 1
  if(n.eq.1)return
  if(incx.eq.1)go to 20

  ! code for increment not equal to 1
  ix = 1
  qmax = abs(qx(1))
  ix = ix + incx
  do i = 2,n
    if(abs(qx(ix)).le.qmax) then
      ix = ix + incx
      cycle
    end if
    iqamax = i
    qmax = abs(qx(ix))
    ix = ix + incx
  end do
  return

  ! code for increment equal to 1
20 qmax = abs(qx(1))
  do i = 2,n
    if(abs(qx(i)).le.qmax) cycle
    iqamax = i
    qmax = abs(qx(i))
  end do

  end function iqamax

! ***************************************************************************

!<subroutine>
  subroutine sdcopy(n,sx,incx,dy,incy)

!<description>
  ! Copies a single precision vector to a double precision vector dy := sx
  !
  ! REMARK: This subroutine is a port of subroutine DCOPY from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Destination vector
  real(DP), dimension(*), intent(out) :: dy

!</output>

!<input>

  ! Source vector
  real(SP), dimension(*), intent(in) :: sx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increments
  integer, intent(in) :: incx,incy

!</input>
!</subroutine>

  ! local variables
  integer :: i,ix,iy,m,mp1

  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20

  ! code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
    dy(iy) = real(sx(ix),DP)
    ix = ix + incx
    iy = iy + incy
  end do
  return

  ! code for both increments equal to 1
20 m = mod(n,7)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    dy(i) = real(sx(i),DP)
  end do
  if( n .lt. 7 ) return
40 mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = mp1,n,7
    dy(i) = real(sx(i),DP)
    dy(i + 1) = real(sx(i + 1),DP)
    dy(i + 2) = real(sx(i + 2),DP)
    dy(i + 3) = real(sx(i + 3),DP)
    dy(i + 4) = real(sx(i + 4),DP)
    dy(i + 5) = real(sx(i + 5),DP)
    dy(i + 6) = real(sx(i + 6),DP)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end subroutine sdcopy

! ***************************************************************************

!<subroutine>
  subroutine dscopy(n,dx,incx,sy,incy)

!<description>
  ! Copies a single precision vector to a double precision vector sy := dx
  !
  ! REMARK: This subroutine is a port of subroutine DCOPY from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Destination vector
  real(SP), dimension(*), intent(out) :: sy

!</output>

!<input>

  ! Source vector
  real(DP), dimension(*), intent(in) :: dx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increments
  integer, intent(in) :: incx,incy

!</input>
!</subroutine>

  ! local variables
  integer :: i,ix,iy,m,mp1

  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20

  ! code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
    sy(iy) = real(dx(ix),SP)
    ix = ix + incx
    iy = iy + incy
  end do
  return

  ! code for both increments equal to 1
20 m = mod(n,7)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    sy(i) = real(dx(i),SP)
  end do
  if( n .lt. 7 ) return
40 mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = mp1,n,7
    sy(i) = real(dx(i),SP)
    sy(i + 1) = real(dx(i + 1),SP)
    sy(i + 2) = real(dx(i + 2),SP)
    sy(i + 3) = real(dx(i + 3),SP)
    sy(i + 4) = real(dx(i + 4),SP)
    sy(i + 5) = real(dx(i + 5),SP)
    sy(i + 6) = real(dx(i + 6),SP)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
end subroutine dscopy

! ***************************************************************************

!<subroutine>
  subroutine sqcopy(n,sx,incx,qy,incy)

!<description>
  ! Copies a single precision vector to a quad precision vector qy := sx
  !
  ! REMARK: This subroutine is a port of subroutine DCOPY from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Destination vector
  real(QP), dimension(*), intent(out) :: qy

!</output>

!<input>

  ! Source vector
  real(SP), dimension(*), intent(in) :: sx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increments
  integer, intent(in) :: incx,incy

!</input>
!</subroutine>

  ! local variables
  integer :: i,ix,iy,m,mp1

  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20

  ! code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
    qy(iy) = real(sx(ix),QP)
    ix = ix + incx
    iy = iy + incy
  end do
  return

  ! code for both increments equal to 1
20 m = mod(n,7)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    qy(i) = real(sx(i),QP)
  end do
  if( n .lt. 7 ) return
40 mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = mp1,n,7
    qy(i) = real(sx(i),QP)
    qy(i + 1) = real(sx(i + 1),QP)
    qy(i + 2) = real(sx(i + 2),QP)
    qy(i + 3) = real(sx(i + 3),QP)
    qy(i + 4) = real(sx(i + 4),QP)
    qy(i + 5) = real(sx(i + 5),QP)
    qy(i + 6) = real(sx(i + 6),QP)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end subroutine sqcopy

! ***************************************************************************

!<subroutine>
  subroutine qscopy(n,qx,incx,sy,incy)

!<description>
  ! Copies a quad precision vector to a single precision vector sy := qx
  !
  ! REMARK: This subroutine is a port of subroutine DCOPY from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Destination vector
  real(SP), dimension(*), intent(out) :: sy

!</output>

!<input>

  ! Source vector
  real(QP), dimension(*), intent(in) :: qx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increments
  integer, intent(in) :: incx,incy

!</input>
!</subroutine>

  ! local variables
  integer :: i,ix,iy,m,mp1

  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20

  ! code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
    sy(iy) = real(qx(ix),SP)
    ix = ix + incx
    iy = iy + incy
  end do
  return

  ! code for both increments equal to 1
20 m = mod(n,7)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    sy(i) = real(qx(i),SP)
  end do
  if( n .lt. 7 ) return
40 mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = mp1,n,7
    sy(i) = real(qx(i),SP)
    sy(i + 1) = real(qx(i + 1),SP)
    sy(i + 2) = real(qx(i + 2),SP)
    sy(i + 3) = real(qx(i + 3),SP)
    sy(i + 4) = real(qx(i + 4),SP)
    sy(i + 5) = real(qx(i + 5),SP)
    sy(i + 6) = real(qx(i + 6),SP)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
end subroutine qscopy

! ***************************************************************************

!<subroutine>
  subroutine dqcopy(n,dx,incx,qy,incy)

!<description>
  ! Copies a double precision vector to a quad precision vector qy := dx
  !
  ! REMARK: This subroutine is a port of subroutine DCOPY from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Destination vector
  real(QP), dimension(*), intent(out) :: qy

!</output>

!<input>

  ! Source vector
  real(DP), dimension(*), intent(in) :: dx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increments
  integer, intent(in) :: incx,incy

!</input>
!</subroutine>

  ! local variables
  integer :: i,ix,iy,m,mp1

  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20

  ! code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
    qy(iy) = real(dx(ix),QP)
    ix = ix + incx
    iy = iy + incy
  end do
  return

  ! code for both increments equal to 1
20 m = mod(n,7)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    qy(i) = real(dx(i),QP)
  end do
  if( n .lt. 7 ) return
40 mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = mp1,n,7
    qy(i) = real(dx(i),QP)
    qy(i + 1) = real(dx(i + 1),QP)
    qy(i + 2) = real(dx(i + 2),QP)
    qy(i + 3) = real(dx(i + 3),QP)
    qy(i + 4) = real(dx(i + 4),QP)
    qy(i + 5) = real(dx(i + 5),QP)
    qy(i + 6) = real(dx(i + 6),QP)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end subroutine dqcopy

! ***************************************************************************

!<subroutine>
  subroutine qdcopy(n,qx,incx,dy,incy)

!<description>
  ! Copies a quad precision vector to a double precision vector dy := qx
  !
  ! REMARK: This subroutine is a port of subroutine DCOPY from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Destination vector
  real(DP), dimension(*), intent(out) :: dy

!</output>

!<input>

  ! Source vector
  real(QP), dimension(*), intent(in) :: qx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increments
  integer, intent(in) :: incx,incy

!</input>
!</subroutine>

  ! local variables
  integer :: i,ix,iy,m,mp1

  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20

  ! code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
    dy(iy) = real(qx(ix),DP)
    ix = ix + incx
    iy = iy + incy
  end do
  return

  ! code for both increments equal to 1
20 m = mod(n,7)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    dy(i) = real(qx(i),DP)
  end do
  if( n .lt. 7 ) return
40 mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = mp1,n,7
    dy(i) = real(qx(i),DP)
    dy(i + 1) = real(qx(i + 1),DP)
    dy(i + 2) = real(qx(i + 2),DP)
    dy(i + 3) = real(qx(i + 3),DP)
    dy(i + 4) = real(qx(i + 4),DP)
    dy(i + 5) = real(qx(i + 5),DP)
    dy(i + 6) = real(qx(i + 6),DP)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end subroutine qdcopy

! ***************************************************************************

!<subroutine>
  subroutine sdaxpy(n,sa,sx,incx,dy,incy)

!<description>
  ! Computes constant times a single vector plus a double precision vector.
  !
  ! REMARK: This subroutine is a port of subroutine DAXPY from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Second source and destination vector
  real(DP), dimension(*), intent(out) :: dy

!</output>

!<input>

  ! First source vector
  real(SP), dimension(*), intent(in) :: sx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increments
  integer, intent(in) :: incx,incy

  ! Multiplication factor
  real(SP), intent(in) :: sa

!</input>
!</subroutine>

  ! local variables
  integer :: i,ix,iy,m,mp1

  if(n.le.0)return
  if (sa .eq. 0.0_SP) return
  if(incx.eq.1.and.incy.eq.1)go to 20

  ! code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
    dy(iy) = dy(iy) + sa*sx(ix)
    ix = ix + incx
    iy = iy + incy
  end do
  return

  ! code for both increments equal to 1
20 m = mod(n,4)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    dy(i) = dy(i) + sa*sx(i)
  end do
  if( n .lt. 4 ) return
40 mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = mp1,n,4
    dy(i) = dy(i) + sa*sx(i)
    dy(i + 1) = dy(i + 1) + sa*sx(i + 1)
    dy(i + 2) = dy(i + 2) + sa*sx(i + 2)
    dy(i + 3) = dy(i + 3) + sa*sx(i + 3)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end subroutine sdaxpy

! ***************************************************************************

!<subroutine>
  subroutine sqaxpy(n,sa,sx,incx,qy,incy)

!<description>
  ! Computes constant times a single vector plus a quad precision vector.
  !
  ! REMARK: This subroutine is a port of subroutine DAXPY from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Second source and destination vector
  real(QP), dimension(*), intent(out) :: qy

!</output>

!<input>

  ! First source vector
  real(SP), dimension(*), intent(in) :: sx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increments
  integer, intent(in) :: incx,incy

  ! Multiplication factor
  real(SP), intent(in) :: sa

!</input>
!</subroutine>

  ! local variables
  integer :: i,ix,iy,m,mp1

  if(n.le.0)return
  if (sa .eq. 0.0_SP) return
  if(incx.eq.1.and.incy.eq.1)go to 20

  ! code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
    qy(iy) = qy(iy) + sa*sx(ix)
    ix = ix + incx
    iy = iy + incy
  end do
  return

  ! code for both increments equal to 1
20 m = mod(n,4)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    qy(i) = qy(i) + sa*sx(i)
  end do
  if( n .lt. 4 ) return
40 mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = mp1,n,4
    qy(i) = qy(i) + sa*sx(i)
    qy(i + 1) = qy(i + 1) + sa*sx(i + 1)
    qy(i + 2) = qy(i + 2) + sa*sx(i + 2)
    qy(i + 3) = qy(i + 3) + sa*sx(i + 3)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end subroutine sqaxpy

! ***************************************************************************

!<subroutine>
  subroutine dqaxpy(n,da,dx,incx,qy,incy)

!<description>
  ! Computes constant times a double vector plus a quad precision vector.
  !
  ! REMARK: This subroutine is a port of subroutine DAXPY from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Second source and destination vector
  real(DP), dimension(*), intent(out) :: qy

!</output>

!<input>

  ! First source vector
  real(DP), dimension(*), intent(in) :: dx

  ! Length of the vector
  integer, intent(in) :: n

  ! Increments
  integer, intent(in) :: incx,incy

  ! Multiplication factor
  real(DP), intent(in) :: da

!</input>
!</subroutine>

  ! local variables
  integer :: i,ix,iy,m,mp1

  if(n.le.0)return
  if (da .eq. 0.0_DP) return
  if(incx.eq.1.and.incy.eq.1)go to 20

  ! code for unequal increments or equal increments not equal to 1
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do i = 1,n
    qy(iy) = qy(iy) + da*dx(ix)
    ix = ix + incx
    iy = iy + incy
  end do
  return

  ! code for both increments equal to 1
20 m = mod(n,4)
  if( m .eq. 0 ) go to 40
  do i = 1,m
    qy(i) = qy(i) + da*dx(i)
  end do
  if( n .lt. 4 ) return
40 mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
  do i = mp1,n,4
    qy(i) = qy(i) + da*dx(i)
    qy(i + 1) = qy(i + 1) + da*dx(i + 1)
    qy(i + 2) = qy(i + 2) + da*dx(i + 2)
    qy(i + 3) = qy(i + 3) + da*dx(i + 3)
  end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end subroutine dqaxpy

! ***************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure&
#endif
  subroutine sset(n,sa,sx,incx)

!<description>
  ! Sets a single precision vector equal to a constant sx := sa
  !
  ! REMARK: This subroutine is an adaptation of subroutine DSCAL from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Destination vector
  real(SP), dimension(*), intent(out) :: sx

!</output>

!<input>

  ! Length of the vector
  integer, intent(in) :: n

  ! Increment
  integer, intent(in) :: incx

  ! Constant
  real(SP), intent(in) :: sa

!</input>
!</subroutine>

  ! local variables
  integer :: i,m,mp1,nincx

  if (n.le.0 .or. incx.le.0) return
  if (incx .eq. 1) then
    ! code for increment equal to 1
    m = mod(n,5)
    if (m .ne. 0) then
      do i = 1,m
        sx(i) = sa
      end do
      if (n .lt. 5) return
    end if
    mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
    do i = mp1,n,5
      sx(i) = sa
      sx(i + 1) = sa
      sx(i + 2) = sa
      sx(i + 3) = sa
      sx(i + 4) = sa
    end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  else
    ! code for increment not equal to 1
    nincx = n*incx
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
    do i = 1,nincx,incx
      sx(i) = sa
    end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end if

  end subroutine sset

! ***************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure&
#endif
  subroutine dset(n,da,dx,incx)

!<description>
  ! Sets a double precision vector equal to a constant dx := da
  !
  ! REMARK: This subroutine is an adaptation of subroutine DSCAL from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Destination vector
  real(DP), dimension(*), intent(out) :: dx

!</output>

!<input>

  ! Length of the vector
  integer, intent(in) :: n

  ! Increment
  integer, intent(in) :: incx

  ! Constant
  real(DP), intent(in) :: da

!</input>
!</subroutine>

  ! local variables
  integer :: i,m,mp1,nincx

  if (n.le.0 .or. incx.le.0) return
  if (incx .eq. 1) then
    ! code for increment equal to 1
    m = mod(n,5)
    if (m .ne. 0) then
      do i = 1,m
        dx(i) = da
      end do
      if (n .lt. 5) return
    end if
    mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
    do i = mp1,n,5
      dx(i) = da
      dx(i + 1) = da
      dx(i + 2) = da
      dx(i + 3) = da
      dx(i + 4) = da
    end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  else
    ! code for increment not equal to 1
    nincx = n*incx
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
    do i = 1,nincx,incx
      dx(i) = da
    end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end if

  end subroutine dset

! ***************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure&
#endif
  subroutine qset(n,qa,qx,incx)

!<description>
  ! Sets a quad precision vector equal to a constant qx := qa
  !
  ! REMARK: This subroutine is an adaptation of subroutine DSCAL from the BLAS
  !         reference implementation by Jack Dongarra, Linpack, 3/11/78.
!</description>

!<output>

  ! Destination vector
  real(QP), dimension(*), intent(out) :: qx

!</output>

!<input>

  ! Length of the vector
  integer, intent(in) :: n

  ! Increment
  integer, intent(in) :: incx

  ! Constant
  real(QP), intent(in) :: qa

!</input>
!</subroutine>

  ! local variables
  integer :: i,m,mp1,nincx

  if (n.le.0 .or. incx.le.0) return
  if (incx .eq. 1) then
    ! code for increment equal to 1
    m = mod(n,5)
    if (m .ne. 0) then
      do i = 1,m
        qx(i) = qa
      end do
      if (n .lt. 5) return
    end if
    mp1 = m + 1
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
    do i = mp1,n,5
      qx(i) = qa
      qx(i + 1) = qa
      qx(i + 2) = qa
      qx(i + 3) = qa
      qx(i + 4) = qa
    end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  else
    ! code for increment not equal to 1
    nincx = n*incx
#ifdef HAS_OPENMP40
  !$omp parallel do simd if(n > lalg_perfconfig%NITEMMIN_OMP)
#else
  !$omp parallel do if(n > lalg_perfconfig%NITEMMIN_OMP)
#endif
    do i = 1,nincx,incx
      qx(i) = qa
    end do
#ifdef HAS_OPENMP40
  !$omp end parallel do simd
#else
  !$omp end parallel do
#endif
  end if

  end subroutine qset

end module

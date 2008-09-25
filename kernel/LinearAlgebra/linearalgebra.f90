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

  interface lalg_copyVector
    module procedure lalg_copyVectorDble
    module procedure lalg_copyVectorSngl
    module procedure lalg_copyVectorSnglDbl
    module procedure lalg_copyVectorDblSngl
    module procedure lalg_copyVectorInt
    module procedure lalg_copyVectorLogical
    module procedure lalg_copyVectorChar
    module procedure lalg_copyVectorDble2D
    module procedure lalg_copyVectorSngl2D
    module procedure lalg_copyVectorSnglDbl2D
    module procedure lalg_copyVectorDblSngl2D
    module procedure lalg_copyVectorInt2D
    module procedure lalg_copyVectorLogical2D
    module procedure lalg_copyVectorChar2D
  end interface

  interface lalg_scaleVector
    module procedure lalg_scaleVectorDble
    module procedure lalg_scaleVectorSngl
    module procedure lalg_scaleVectorInt
    module procedure lalg_scaleVectorDble2D
    module procedure lalg_scaleVectorSngl2D
    module procedure lalg_scaleVectorInt2D
  end interface

  interface lalg_clearVector
    module procedure lalg_clearVectorDble
    module procedure lalg_clearVectorSngl
    module procedure lalg_clearVectorInt
    module procedure lalg_clearVectorDble2D
    module procedure lalg_clearVectorSngl2D
    module procedure lalg_clearVectorInt2D
  end interface

  interface lalg_setVector
    module procedure lalg_setVectorDble
    module procedure lalg_setVectorSngl
    module procedure lalg_setVectorInt
    module procedure lalg_setVectorLogical
    module procedure lalg_setVectorChar
    module procedure lalg_setVectorDble2D
    module procedure lalg_setVectorSngl2D
    module procedure lalg_setVectorInt2D
    module procedure lalg_setVectorLogical2D
    module procedure lalg_setVectorChar2D
  end interface

  interface lalg_vectorLinearComb
    module procedure lalg_vectorLinearCombDble
    module procedure lalg_vectorLinearCombSngl
    module procedure lalg_vectorLinearCombSnglDble
    module procedure lalg_vectorLinearCombDble2D
    module procedure lalg_vectorLinearCombSngl2D
    module procedure lalg_vectorLinearCombSnglDble2D
  end interface

  interface lalg_scalarProduct
    module procedure lalg_scalarProductDble
    module procedure lalg_scalarProductSngl
    module procedure lalg_scalarProductInt
    module procedure lalg_scalarProductDble2D
    module procedure lalg_scalarProductSngl2D
    module procedure lalg_scalarProductInt2D
  end interface

  interface lalg_norm
    module procedure lalg_normDble
    module procedure lalg_normSngl
  end interface

  interface lalg_errorNorm
    module procedure lalg_errorNormDble
    module procedure lalg_errorNormSngl
  end interface
    
  interface lalg_vectorSort
    module procedure lalg_vectorSortDble
    module procedure lalg_vectorSortSngl
    module procedure lalg_vectorSortInt
  end interface

  interface lalg_vectorAddScalar
    module procedure lalg_vectorAddScalarDble
    module procedure lalg_vectorAddScalarSngl
    module procedure lalg_vectorAddScalarInt
    module procedure lalg_vectorAddScalarDble2D
    module procedure lalg_vectorAddScalarSngl2D
    module procedure lalg_vectorAddScalarInt2D
  end interface

  interface lalg_vectorCompMult
    module procedure lalg_vectorCompMultDble
    module procedure lalg_vectorCompMultSngl
    module procedure lalg_vectorCompMultInt
    module procedure lalg_vectorCompMultDbleSngl
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
    call SCOPY(size(Fx),Fx,1,Fy,1)
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
  integer(I32) :: i
  
  if (.not. present(n)) then
  
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i=1,size(Fx)
      Dy(i) = real(Fx(i),DP)
    end do
  !%OMP  end parallel do
  
  else

  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i=1,n
      Dy(i) = real(Fx(i),DP)
    end do
  !%OMP  end parallel do
  
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
  integer(I32) :: i
  
  if (.not. present(n)) then

  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i=1,size(Dx)
      Fy(i) = real(Dx(i),SP)
    end do
  !%OMP  end parallel do

  else
  
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i=1,n
      Fy(i) = real(Dx(i),SP)
    end do
  !%OMP  end parallel do

  end if  

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorInt (Ix,Iy,n)
  
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

  integer(I32) :: i
  
  if (.not. present(n)) then
  
    ! Does not exist in BLAS!
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i=1,size(Ix)
      Iy(i) = Ix(i)
    end do
  !%OMP  end parallel do

  else
  
    ! Does not exist in BLAS!
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i=1,n
      Iy(i) = Ix(i)
    end do
  !%OMP  end parallel do

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
  integer(I32) :: i
  
  if (.not. present(n)) then
  
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i=1,size(Lx)
      Ly(i) = Lx(i)
    end do
  !%OMP  end parallel do
  
  else

  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i=1,n
      Ly(i) = Lx(i)
    end do
  !%OMP  end parallel do
  
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
  integer(I32) :: i
  
  if (.not. present(n)) then
  
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i=1,size(Sx)
      Sy(i) = Sx(i)
    end do
  !%OMP  end parallel do
  
  else

  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i=1,n
      Sy(i) = Sx(i)
    end do
  !%OMP  end parallel do
  
  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDble2D (Dx,Dy)
  
!<description>
  ! Copies a double precision vector dx: Dy = Dx
!</description>

!<input>
  
  ! Source vector
  real(DP), dimension(:,:), intent(IN) :: Dx
  
!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:,:), intent(OUT) :: Dy
  
!</output>
  
!</subroutine>

  call DCOPY(size(Dx,1)*size(Dx,2),Dx,1,Dy,1)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSngl2D (Fx,Fy)
  
!<description>
  ! Copies a single precision vector: Fy = Fx
!</description>

!<input>
  
  ! Source vector
  real(SP), dimension(:,:), intent(IN) :: Fx
  
!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:,:), intent(OUT) :: Fy
  
!</output>
  
!</subroutine>

  call SCOPY(size(Fx,1)*size(Fx,2),Fx,1,Fy,1)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorSnglDbl2D (Fx,Dy)
  
!<description>
  ! Copies single precision vector to double precision vector: Dy = Fx
!</description>

!<input>
  
  ! Source vector
  real(SP), dimension(:,:), intent(IN) :: Fx
  
!</input>

!<output>
  
  ! Destination vector
  real(DP), dimension(:,:), intent(OUT) :: Dy
  
!</output>
  
!</subroutine>
  integer(I32) :: i,j
  
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j=1,size(Fx,2)
    do i=1,size(Fx,1)
      Dy(i,j) = real(Fx(i,j),DP)
    end do
  end do
!%OMP  end parallel do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorDblSngl2D (Dx,Fy)
  
!<description>
  ! Copies double precision vector to single precision vector: Fy = Dx
!</description>

!<input>
  
  ! Source vector
  real(DP), dimension(:,:), intent(IN) :: Dx
  
!</input>

!<output>
  
  ! Destination vector
  real(SP), dimension(:,:), intent(OUT) :: Fy
  
!</output>
  
!</subroutine>
  integer(I32) :: i,j
  
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j=1,size(Dx,2)
    do i=1,size(Dx,1)
      Fy(i,j) = real(Dx(i,j),SP)
    end do
  end do
!%OMP  end parallel do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorInt2D (Ix,Iy)
  
!<description>
  ! Copies an integer vector Ix: Iy = Ix
!</description>

!<input>
  
  ! Source vector
  integer(I32), dimension(:,:), intent(IN) :: Ix
  
!</input>

!<output>
  
  ! Destination vector
  integer(I32), dimension(:,:), intent(OUT) :: Iy
  
!</output>
  
!</subroutine>

  integer(I32) :: i,j
  
  ! Does not exist in BLAS!
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j=1,size(Ix,2)
    do i=1,size(Ix,1)
      Iy(i,j) = Ix(i,j)
    end do
  end do
!%OMP  end parallel do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorLogical2D (Lx,Ly)
  
!<description>
  ! Copies a logical vector Lx: Ly = Lx
!</description>

!<input>
  
  ! Source vector
  logical, dimension(:,:), intent(IN) :: Lx
  
!</input>

!<output>
  
  ! Destination vector
  logical, dimension(:,:), intent(OUT) :: Ly
  
!</output>
  
!</subroutine>

  integer(I32) :: i,j
  
  ! Does not exist in BLAS!
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j=1,size(Lx,2)
    do i=1,size(Lx,1)
      Ly(i,j) = Lx(i,j)
    end do
  end do
!%OMP  end parallel do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_copyVectorChar2D (Sx,Sy)
  
!<description>
  ! Copies a character vector Sx: Sy = Sx
!</description>

!<input>
  
  ! Source vector
  character, dimension(:,:), intent(IN) :: Sx
  
!</input>

!<output>
  
  ! Destination vector
  character, dimension(:,:), intent(OUT) :: Sy
  
!</output>
  
!</subroutine>

  integer(I32) :: i,j
  
  ! Does not exist in BLAS!
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j=1,size(Sx,2)
    do i=1,size(Sx,1)
      Sy(i,j) = Sx(i,j)
    end do
  end do
!%OMP  end parallel do
  
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
    call DSCAL(size(Dx),dc,Dx,1)
  else
    call DSCAL(n,dc,Dx,1)
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
    call SSCAL(size(Fx),sc,Fx,1)
  else
    call SSCAL(n,sc,Fx,1)
  end if
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine lalg_scaleVectorInt (Ix,ic,n)
  
!<description>
  ! Scales a integer vector: Ix = c * Ix
!</description>

!<inputoutput>
  
  ! Source and destination vector
  real(DP), dimension(:), intent(INOUT) :: Ix
  
!</inputoutput>

!<input>

  ! Multiplication factor
  integer(I32), intent(IN) :: ic

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n

!</input>
  
!</subroutine>

  integer(I32) :: i
  
  if (.not. present(n)) then
    ! Does not exist in BLAS
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i=1,size(Ix)
      Ix(i) = ic*Ix(i)
    end do
  !%OMP  end parallel do
  else
    ! Does not exist in BLAS
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i=1,n
      Ix(i) = ic*Ix(i)
    end do
  !%OMP  end parallel do
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

  call DSCAL(size(Dx,1)*size(Dx,2),dc,Dx,1)
  
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

  call SSCAL(size(Fx,1)*size(Fx,2),sc,Fx,1)
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine lalg_scaleVectorInt2D (Ix,ic)
  
!<description>
  ! Scales a integer vector: Ix = c * Ix
!</description>

!<inputoutput>
  
  ! Source and destination vector
  real(DP), dimension(:,:), intent(INOUT) :: Ix
  
!</inputoutput>

!<input>

  ! Multiplication factor
  integer(I32), intent(IN) :: ic

!</input>
  
!</subroutine>

  integer(I32) :: i,j
  
  ! Does not exist in BLAS
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j=1,size(Ix,2)
    do i=1,size(Ix,1)
      Ix(i,j) = ic*Ix(i,j)
    end do
  end do
!%OMP  end parallel do

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
  integer(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

  if (.not. present(n)) then
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,size(Dx)
      Dx(i) = 0.0_DP
    end do
  !%OMP  end parallel do
  else
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,n
      Dx(i) = 0.0_DP
    end do
  !%OMP  end parallel do
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
  integer(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

  if (.not. present(n)) then
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,size(Fx)
      Fx(i) = 0.0_SP
    end do
  !%OMP  end parallel do
  else
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,n
      Fx(i) = 0.0_SP
    end do
  !%OMP  end parallel do
  end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorInt (Ix,n)
  
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
  integer(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

  if (.not. present(n)) then
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,size(Ix)
      Ix(i) = 0
    end do
  !%OMP  end parallel do
  else
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,size(Ix)
      Ix(i) = 0
    end do
  !%OMP  end parallel do
  end if
  
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
  integer(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j = 1,size(Dx,2)
    do i = 1,size(Dx,1)
      Dx(i,j) = 0.0_DP
    end do
  end do
!%OMP  end parallel do
  
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
  integer(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j = 1,size(Fx,2)
    do i = 1,size(Fx,1)
      Fx(i,j) = 0.0_SP
    end do
  end do
!%OMP  end parallel do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_clearVectorInt2D (Ix)
  
!<description>
  ! Clears an integer vector: Ix = 0
!</description>

!<output>
  
  ! Destination vector to be cleared
  integer(I32), dimension(:,:), intent(OUT) :: Ix
  
!</output>
  
!</subroutine>

  ! local variables
  integer(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = 0
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j = 1,size(Ix,2)
    do i = 1,size(Ix,1)
      Ix(i,j) = 0
    end do
  end do
!%OMP  end parallel do
  
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
  integer(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = dvalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

  if (.not. present(n)) then
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,size(Dx)
      Dx(i) = dvalue
    end do
  !%OMP  end parallel do
  else
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,size(Dx)
      Dx(i) = dvalue
    end do
  !%OMP  end parallel do
  end if
  
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
  integer(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Fx = fvalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

  if (.not. present(n)) then
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,size(Fx)
      Fx(i) = fvalue
    end do
  !%OMP  end parallel do
  else
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,n
      Fx(i) = fvalue
    end do
  !%OMP  end parallel do
  end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorInt (Ix,ivalue,n)
  
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
  integer(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Ix = ivalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

  if (.not. present(n)) then
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,size(Ix)
      Ix(i) = ivalue
    end do
  !%OMP  end parallel do
  else
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,n
      Ix(i) = ivalue
    end do
  !%OMP  end parallel do
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
  integer(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Lx = lvalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

  if (.not. present(n)) then
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,size(Lx)
      Lx(i) = lvalue
    end do
  !%OMP  end parallel do
  else
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,n
      Lx(i) = lvalue
    end do
  !%OMP  end parallel do
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
  integer(I32) :: i
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Sx = svalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

  if (.not. present(n)) then
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,size(Sx)
      Sx(i) = svalue
    end do
  !%OMP  end parallel do
  else
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
    do i = 1,n
      Sx(i) = svalue
    end do
  !%OMP  end parallel do
  end if
  
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
  integer(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Dx = dvalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j = 1,size(Dx,2)
    do i = 1,size(Dx,1)
      Dx(i,j) = dvalue
    end do
  end do
!%OMP  end parallel do
  
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
  integer(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Fx = fvalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j = 1,size(Fx,2)
    do i = 1,size(Fx,1)
      Fx(i,j) = fvalue
    end do
  end do
!%OMP  end parallel do
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lalg_setVectorInt2D (Ix,ivalue)
  
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
  integer(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Ix = ivalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j = 1,size(Ix,2)
    do i = 1,size(Ix,1)
      Ix(i,j) = ivalue
    end do
  end do
!%OMP  end parallel do
  
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
  integer(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Lx = lvalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j = 1,size(Lx,2)
    do i = 1,size(Lx,1)
      Lx(i,j) = lvalue
    end do
  end do
!%OMP  end parallel do
  
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
  integer(I32) :: i,j
  
  ! We trust the internal loop unrolling functions of the compiler...
  ! normally we can use:
  
  !Lx = lvalue
  
  ! But if the compiler does not support that, maybe we have to go back
  ! to the standard DO loop...

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
  do j = 1,size(Sx,2)
    do i = 1,size(Sx,1)
      Sx(i,j) = svalue
    end do
  end do
!%OMP  end parallel do
  
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
  integer(I32) :: i
  real(DP) :: c
  
  if (.not. present(n)) then
  
    if (dcy .eq. 0.0_DP) then
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
      do i=1,size(Fx)
        Dy(i) = Fx(i)
      end do
  !%OMP  end parallel do
      if (scx .ne. 1.0_SP) then
  !%OMP  parallel do &
  !%OMP& default(shared) & 
  !%OMP& private(i)
        do i=1,size(Fx)
          Dy(i) = scx*Fx(i)
        end do
  !%OMP  end parallel do
      end if
    else if (dcy .eq. 1.0_DP) then
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
      do i=1,size(Fx)
        Dy(i) = Dy(i) + scx*Fx(i)
      end do
  !%OMP  end parallel do
    else
      c=scx/dcy
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
      do i=1,size(Fx)
        Dy(i) = Dy(i) + c*Fx(i)
      end do
  !%OMP  end parallel do
      call DSCAL(size(Dy),dcy,Dy,1)
    endif
    
  else
  
    if (dcy .eq. 0.0_DP) then
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
      do i=1,n
        Dy(i) = Fx(i)
      end do
  !%OMP  end parallel do
      if (scx .ne. 1.0_SP) then
  !%OMP  parallel do &
  !%OMP& default(shared) & 
  !%OMP& private(i)
        do i=1,n
          Dy(i) = scx*Fx(i)
        end do
  !%OMP  end parallel do
      end if
    else if (dcy .eq. 1.0_DP) then
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
      do i=1,n
        Dy(i) = Dy(i) + scx*Fx(i)
      end do
  !%OMP  end parallel do
    else
      c=scx/dcy
  !%OMP  parallel do &
  !%OMP& default(shared) &
  !%OMP& private(i)
      do i=1,n
        Dy(i) = Dy(i) + c*Fx(i)
      end do
  !%OMP  end parallel do
      call DSCAL(n,dcy,Dy,1)
    endif
    
  end if
  
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

  ! local variables
  integer(I32) :: i,j
  real(DP) :: c
  
  if (dcy .eq. 0.0_DP) then
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
    do j=1,size(Fx,2)
      do i=1,size(Fx,1)
        Dy(i,j) = Fx(i,j)
      end do
    end do
!%OMP  end parallel do
    if (scx .ne. 1.0_SP) then
!%OMP  parallel do &
!%OMP& default(shared) & 
!%OMP& private(i,j)
      do j=1,size(Fx,2)
        do i=1,size(Fx,1)
          Dy(i,j) = scx*Fx(i,j)
        end do
      end do
!%OMP  end parallel do
    end if
  else if (dcy .eq. 1.0_DP) then
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
    do j=1,size(Fx,2)
      do i=1,size(Fx,1)
        Dy(i,j) = Dy(i,j) + scx*Fx(i,j)
      end do
    end do
!%OMP  end parallel do
  else
    c=scx/dcy
!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i,j)
    do j=1,size(Fx,2)
      do i=1,size(Fx,1)
        Dy(i,j) = Dy(i,j) + c*Fx(i,j)
      end do
    end do
!%OMP  end parallel do
    call DSCAL(size(Dy,1)*size(Dy,2),dcy,Dy,1)
  endif
  
  end subroutine
  
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

  ! local variables
  !INTEGER(I32) :: i
  real(DP) :: DDOT

  if (.not. present(n)) then
    res=DDOT(size(Dx),Dx,1,Dy,1)
!!$  res = Dx(1)*Dy(1)
!!$  DO i=2,SIZE(Dx)
!!$    res = res + Dx(i)*Dy(i)
!!$  END DO
  else
    res=DDOT(n,Dx,1,Dy,1)
!!$  res = Dx(1)*Dy(1)
!!$  DO i=2,n
!!$    res = res + Dx(i)*Dy(i)
!!$  END DO
  end if
  
  end function

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

  ! local variables
  !INTEGER(I32) :: i
  real(SP) :: SDOT

  if (.not. present(n)) then
    res=SDOT(size(Fx),Fx,1,Fy,1)
  !!$  res = Fx(1)*Fy(1)
  !!$  DO i=2,SIZE(Fx)
  !!$    res = res + Fx(i)*Fy(i)
  !!$  END DO
  else
    res=SDOT(size(Fx),Fx,1,Fy,1)
  !!$  res = Fx(1)*Fy(1)
  !!$  DO i=2,SIZE(Fx)
  !!$    res = res + Fx(i)*Fy(i)
  !!$  END DO
  end if
  
  end function

  ! ***************************************************************************

!<function>

  real(SP) function lalg_scalarProductInt (Ix,Iy,n) result (res)
  
!<description>
  ! Calculates the scalar product of two single precision vectors: 
  ! res = (vector,vector)
!</description>

!<input>
  
  ! First source vector
  integer(I32), dimension(:), intent(IN) :: Ix
  
  ! Second source vector
  integer(I32), dimension(:), intent(IN) :: Iy
  
  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  ! local variables
  integer(I32) :: i
  
  if (.not. present(n)) then
    res = Ix(1)*Iy(1)
    do i=2,size(Ix)
      res = res + Ix(i)*Iy(i)
    end do
  else
    res = Ix(1)*Iy(1)
    do i=2,n
      res = res + Ix(i)*Iy(i)
    end do
  end if
  
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

  ! local variables
  !INTEGER(I32) :: i
  real(DP) :: DDOT

  res=DDOT(size(Dx,1)*size(Dx,2),Dx,1,Dy,1)
!!$  res = Dx(1)*Dy(1)
!!$  DO i=2,SIZE(Dx)
!!$    res = res + Dx(i)*Dy(i)
!!$  END DO
  
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

  ! local variables
  !INTEGER(I32) :: i
  real(SP) :: SDOT

  res=SDOT(size(Fx,1)*size(Fx,2),Fx,1,Fy,1)
!!$  res = Fx(1)*Fy(1)
!!$  DO i=2,SIZE(Fx)
!!$    res = res + Fx(i)*Fy(i)
!!$  END DO
  
  end function

  ! ***************************************************************************

!<function>

  real(SP) function lalg_scalarProductInt2D (Ix,Iy) result (res)
  
!<description>
  ! Calculates the scalar product of two single precision vectors: 
  ! res = (vector,vector)
!</description>

!<input>
  
  ! First source vector
  integer(I32), dimension(:,:), intent(IN) :: Ix
  
  ! Second source vector
  integer(I32), dimension(:,:), intent(IN) :: Iy
  
!</input>

!<result>
  ! The scalar product of the two vectors.
!</result>

!</function>

  ! local variables
  integer(I32) :: i,j
  
  res = Ix(1,1)*Iy(1,1)
  do j=1,size(Ix,2)
    do i=2,size(Ix,1)
      res = res + Ix(i,j)*Iy(i,j)
    end do
  end do
  
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
  integer(I32), intent(OUT), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! < 0 >1, if an error occurred (unknown norm).
!</result>
  
!</subroutine>

  integer(I32) :: i,j
  integer :: isize
  
  isize = size(Dx)
  if (present(n)) isize = n

  ! Choose the norm to calculate
  select case (cnorm)
  case (LINALG_NORMSUM)
    ! L1-norm: sum all entries
    resnorm = abs(Dx(1))
    do i=2,isize
      resnorm = resnorm + abs(Dx(i))
    end do

  case (LINALG_NORMEUCLID)
    ! Euclidian norm = scalar product (vector,vector)
    resnorm = lalg_scalarProductDble(Dx,Dx,isize)
!!$    resnorm = Dx(1)*Dx(1)
!!$    DO i=2,isize
!!$      resnorm = resnorm + Dx(i)*Dx(i)
!!$    END DO
    resnorm = sqrt(resnorm)

  case (LINALG_NORML1)
    ! L1-norm: sum all entries, divide by sqrt(vector length).
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = abs(Dx(1))
    do i=2,isize
      resnorm = resnorm + abs(Dx(i))
    end do
    resnorm = resnorm / real(isize,DP)

  case (LINALG_NORML2)
    ! l_2-norm - like euclidian norm, but divide by vector length.
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = lalg_scalarProductDble(Dx,Dx,isize)
!!$    resnorm = Dx(1)*Dx(1)
!!$    DO i=2,isize
!!$      resnorm = resnorm + Dx(i)*Dx(i)
!!$    END DO
    resnorm = sqrt(resnorm / real(isize,DP))
    
  case (LINALG_NORMMAX)
    ! MAX-norm. Find the absolute largest entry.
    resnorm = abs(Dx(1))
    j=1
    do i=2,isize
      if (abs(Dx(i)) .gt. resnorm) then
        j = i
        resnorm = abs(Dx(i))
      end if
    end do
    if (present(iposMax)) iposMax = j
  case DEFAULT
    resnorm = -1.0_DP ! Unknown norm.
  end select
    
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
  integer(I32), intent(OUT), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! < 0, if an error occurred (unknown norm).
!</result>
  
!</subroutine>

  integer(I32) :: i,j
  integer :: isize
  
  isize = size(Fx)
  if (present(n)) isize = n

  ! Choose the norm to calculate
  select case (cnorm)
  case (LINALG_NORMSUM)
    ! L1-norm: sum absolute value of all entries
    resnorm = Fx(1)
    do i=2,isize
      resnorm = resnorm + abs(Fx(i))
    end do

  case (LINALG_NORMEUCLID)
    ! Euclidian norm = scalar product (vector,vector)
    resnorm = lalg_scalarProductSngl(Fx,Fx)
!!$    resnorm = Fx(1)*Fx(1)
!!$    DO i=2,isize
!!$      resnorm = resnorm + Fx(i)*Fx(i)
!!$    END DO
    resnorm = sqrt(resnorm)

  case (LINALG_NORML1)
    ! L1-norm: sum sum absolute value of all entries, divide by sqrt(vector length).
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = Fx(1)
    do i=2,isize
      resnorm = resnorm + abs(Fx(i))
    end do
    resnorm = resnorm / real(isize,SP)

  case (LINALG_NORML2)
    ! l_2-norm - like euclidian norm, but divide by sqrt(vector length).
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = lalg_scalarProductSNGL(Fx,Fx)
!!$    resnorm = Fx(1)*Fx(1)
!!$    DO i=2,isize
!!$      resnorm = resnorm + Fx(i)*Fx(i)
!!$    END DO
    resnorm = sqrt(resnorm / real(isize,SP))
    
  case (LINALG_NORMMAX)
    ! MAX-norm. Find the absolute largest entry.
    resnorm = abs(Fx(1))
    j=1
    do i=2,isize
      if (abs(Fx(i)) .gt. resnorm) then
        j = i
        resnorm = abs(Fx(i))
      end if
    end do
    if (present(iposMax)) iposMax = j
  case DEFAULT
    resnorm = -1.0_SP ! Unknown norm.
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
  integer(I32), intent(OUT), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! < 0 >1, if an error occurred (unknown norm).
!</result>
  
!</subroutine>

  real(DP) :: dtemp
  integer(I32) :: i,j
  integer :: isize
  
  isize = size(Dx)
  if (present(n)) isize = n

  ! Choose the norm to calculate
  select case (cnorm)
  case (LINALG_NORMSUM)
    ! L1-norm: sum all entries
    if (present(Dw)) then
      resnorm = Dw(1)*abs(Dx(1)-Dy(1))
      do i=2,isize
        resnorm = resnorm + Dw(i)*abs(Dx(i)-Dy(i))
      end do
    else
      resnorm = abs(Dx(1)-Dy(1))
      do i=2,isize
        resnorm = resnorm + abs(Dx(i)-Dy(i))
      end do
    end if

  case (LINALG_NORMEUCLID)
    ! Euclidian norm = scalar product (vector,vector)
    if (present(Dw)) then
      resnorm = Dw(1)*(Dx(1)-Dy(1))*(Dx(1)-Dy(1))
      do i=2,isize
        resnorm = resnorm + Dw(i)*(Dx(i)-Dy(i))*(Dx(i)-Dy(i))
      end do
    else
      resnorm = (Dx(1)-Dy(1))*(Dx(1)-Dy(1))
      do i=2,isize
        resnorm = resnorm + (Dx(i)-Dy(i))*(Dx(i)-Dy(i))
      end do
    end if
    resnorm = sqrt(resnorm)

  case (LINALG_NORML1)
    ! L1-norm: sum all entries, divide by sqrt(vector length).
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = abs(Dx(1)-Dy(1))
    do i=2,isize
      resnorm = resnorm + abs(Dx(i)-Dy(i))
    end do
    resnorm = resnorm / real(isize,DP)

  case (LINALG_NORML2)
    ! l_2-norm - like euclidian norm, but divide by vector length.
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = (Dx(1)-Dy(1))*(Dx(1)-Dy(1))
    do i=2,isize
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
  case DEFAULT
    resnorm = -1.0_DP ! Unknown norm.
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
  integer(I32), intent(OUT), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! < 0 >1, if an error occurred (unknown norm).
!</result>
  
!</subroutine>

  real(SP) :: stemp
  integer(I32) :: i,j
  integer :: isize
  
  isize = size(Fx)
  if (present(n)) isize = n

  ! Choose the norm to calculate
  select case (cnorm)
  case (LINALG_NORMSUM)
    ! L1-norm: sum all entries
    resnorm = abs(Fx(1)-Fy(1))
    do i=2,isize
      resnorm = resnorm + abs(Fx(i)-Fy(i))
    end do

  case (LINALG_NORMEUCLID)
    ! Euclidian norm = scalar product (vector,vector)
    resnorm = (Fx(1)-Fy(1))*(Fx(1)-Fy(1))
    do i=2,isize
      resnorm = resnorm + (Fx(i)-Fy(i))*(Fx(i)-Fy(i))
    end do
    resnorm = sqrt(resnorm)

  case (LINALG_NORML1)
    ! L1-norm: sum all entries, divide by sqrt(vector length).
    ! So, scale such that the vector (1111...) to has norm = 1.
    resnorm = abs(Fx(1)-Fy(1))
    do i=2,isize
      resnorm = resnorm + abs(Fx(i)-Fy(i))
    end do
    resnorm = resnorm / real(isize,SP)

  case (LINALG_NORML2)
    ! l_2-norm - like euclidian norm, but divide by vector length.
    ! So, scale such that the vektor (1111...) to has norm = 1.
    resnorm = (Fx(1)-Fy(1))*(Fx(1)-Fy(1))
    do i=2,isize
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
  case DEFAULT
    resnorm = -1.0_DP ! Unknown norm.
  end select
    
  end function

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
    integer(I32), dimension(:), intent(IN) :: Itr
  !</input>
    
  !<output>
    ! The resorted vector
    real(DP), dimension(:), intent(OUT) :: Dd
  !</output>
    
!</subroutine>
    
    ! local variable
    integer(I32) :: ieq
    
    ieq = size(Itr)
    ieq = size(Dx)
    ieq = size(Dd)
    
    do ieq=1, size(Itr)
      Dd(ieq) = Dx(Itr(ieq))
    end do
  
  end subroutine 

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
    integer(I32), dimension(:), intent(IN) :: Itr

    ! Source vector to be sorted
    real(SP), dimension(:), intent(IN) :: Fx
  !</input>
    
  !<output>
    ! The resorted vector
    real(SP), dimension(:), intent(OUT) :: Fd
  !</output>
    
!</subroutine>
    
    ! local variable
    integer(I32) :: ieq
    
    do ieq=1, size(Itr)
      Fd(ieq) = Fx(Itr(ieq))
    end do
  
  end subroutine 

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorSortInt (Ix, Id, Itr)
  
  !<description>
    ! Resorts the entries in the vector Ix corresponding to Itr.
    ! In particular, the first SIZE(Itr) entries of Dx are written resortedly
    ! to Dd.
  !</description>
    
  !<input>
    ! Array with permutation of 1..neq
    integer(I32), dimension(:), intent(IN) :: Itr

    ! Source vector to be sorted
    integer(I32), dimension(:), intent(IN) :: Ix
  !</input>
    
  !<output>
    ! The resorted vector
    integer(I32), dimension(:), intent(OUT) :: Id
  !</output>
    
!</subroutine>
    
    ! local variable
    integer(I32) :: ieq
    
    do ieq=1, size(Itr)
      Id(ieq) = Ix(Itr(ieq))
    end do
  
  end subroutine 

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

    do i=1,size(Dy)
      do j=1,size(Dx)
        Dtensor(j,i)=Dx(j)*Dy(i)
      end do
    end do
  end subroutine lalg_tensorProductDble

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
    
    real(DP) :: dval
    integer(I32) :: i
    
    if (.not. present(n)) then
      dval = dvalue
      do i=1,size(Dx)
        Dx(i) = Dx(i) + dval
      end do
    else
      dval = dvalue
      do i=1,n
        Dx(i) = Dx(i) + dval
      end do
    end if
    
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
    
    real(SP) :: fval
    integer(I32) :: i
    
    if (.not. present(n)) then
      fval = fvalue
      do i=1,size(Fx)
        Fx(i) = Fx(i) + fval
      end do
    else
      fval = fvalue
      do i=1,n
        Fx(i) = Fx(i) + fval
      end do
    end if
    
  end subroutine
  
  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarInt (Ix,ivalue,n)
  
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
    
    integer(I32) :: ival
    integer(I32) :: i
    
    if (.not. present(n)) then
      ival = ivalue
      do i=1,size(Ix)
        Ix(i) = Ix(i) + ival
      end do
    else
      ival = ivalue
      do i=1,n
        Ix(i) = Ix(i) + ival
      end do
    end if
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarDble2D (Dx,dvalue)
  
!<description>
  ! This routine adds the value dvalue to each entry of the vector Id.
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
    
    real(DP) :: dval
    integer(I32) :: i,j
    
    dval = dvalue
    do j=1,size(Dx,2)
      do i=1,size(Dx,1)
        Dx(i,j) = Dx(i,j) + dval
      end do
    end do
    
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
    
    real(SP) :: fval
    integer(I32) :: i,j
    
    fval = fvalue
    do j=1,size(Fx,2)
      do i=1,size(Fx,1)
      Fx(i,j) = Fx(i,j) + fval
    end do
    end do
    
  end subroutine
  
  !****************************************************************************

!<subroutine>

  subroutine lalg_vectorAddScalarInt2D (Ix,ivalue)
  
!<description>
  ! This routine adds the value ivalue to each entry of the vector Ix.
!</description>
  
!<input>
  ! The value to add to every entry.
  integer(I32) :: ivalue
!</input>

!<inputoutput>
  ! Data array to be modified.
  integer(I32), dimension(:,:), intent(INOUT) :: Ix
!</inputoutput>

!</subroutine>
    
    integer(I32) :: ival
    integer(I32) :: i,j
    
    ival = ivalue
    do j=1,size(Ix,2)
      do i=1,size(Ix,1)
        Ix(i,j) = Ix(i,j) + ival
      end do
    end do
    
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
  integer(I32) :: i

  if (.not. present(n)) then

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
    do i = 1, size(Dy,1)
      Dy(i) = dc*Dx(i)*Dy(i)
    end do
!%OMP  end parallel do

  else

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
    do i = 1, n
      Dy(i) = dc*Dx(i)*Dy(i)
    end do
!%OMP  end parallel do

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
  integer(I32) :: i

  if (.not. present(n)) then

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
    do i = 1, size(Fy,1)
      Fy(i) = sc*Fx(i)*Fy(i)
    end do
!%OMP  end parallel do

  else

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
    do i = 1, n
      Fy(i) = sc*Fx(i)*Fy(i)
    end do
!%OMP  end parallel do

  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lalg_vectorCompMultInt (Ix,Iy,ic,n)
  
!<description>
  ! Performs componentwise multiplication: Iy = ic * Ix * Iy
!</description>

!<input>
  
  ! First source vector
  integer(I32), dimension(:), intent(IN) :: Ix
  
  ! Scaling factor
  integer(I32), intent(IN)               :: ic

  ! OPTIONAL: Size of the vector
  integer, intent(IN), optional :: n
  
!</input>
 
!<inputoutput>
  
  ! Second source vector; also receives the result
  integer(I32), dimension(:), intent(INOUT) :: Iy
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer(I32) :: i

  if (.not. present(n)) then

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
    do i = 1, size(Iy,1)
      Iy(i) = ic*Ix(i)*Iy(i)
    end do
!%OMP  end parallel do

  else

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
    do i = 1, n
      Iy(i) = ic*Ix(i)*Iy(i)
    end do
!%OMP  end parallel do

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
  integer(I32) :: i

  if (.not. present(n)) then

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
    do i = 1, size(Dy,1)
      Dy(i) = dc*Fx(i)*Dy(i)
    end do
!%OMP  end parallel do

  else

!%OMP  parallel do &
!%OMP& default(shared) &
!%OMP& private(i)
    do i = 1, n
      Dy(i) = dc*Fx(i)*Dy(i)
    end do
!%OMP  end parallel do

  end if

  end subroutine
 
end module

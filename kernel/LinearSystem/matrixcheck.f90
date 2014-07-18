!##############################################################################
!# ****************************************************************************
!# <name> matrixcheck </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of basic testing routines for matrices.
!#
!# The following routines can be found in this module:
!#
!# 1.) mchk_isDiagDominant
!#     -> Checks if the matrix is diagonally dominant
!#
!# 2.) mchk_isMMatrix
!#     -> Checks if the matrix satisfies sufficient conditions of an M-matrix.
!#
!# 3.) mchk_isZMatrix
!#     -> Checks if the matrix satisfies sufficient conditions of an Z-matrix.
!#
!# The following auxiliary routines can be found in this module:
!#
!# 1.) mchk_calcStrongConnComp
!#     -> Calculates the strongly connected components of a directed graph
!# </purpose>
!##############################################################################

module matrixcheck

!$ use omp_lib
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemscalar
  use stackInt
  use storage

  implicit none

  private

!<constants>

!<constantblock description="Global constants for matrix checks">

  ! Check if matrix is diagonally dominant
  integer, parameter, public :: MCHK_DIAGDOMINANT             = 0

  ! Check if matrix is strictly diagonally dominant
  integer, parameter, public :: MCHK_DIAGDOMINANT_STRICT      = 1

  ! Check if matrix is irreducibly diagonally dominant
  integer, parameter, public :: MCHK_DIAGDOMINANT_IRREDUCIBLE = 2

!</constantblock>

!</constants>

  public :: mchk_isDiagDominant
  public :: mchk_isMMatrix
  public :: mchk_isZMatrix

contains

  ! ***************************************************************************

!<subroutine>

  subroutine mchk_isDiagDominant(rmatrix, ccheckType, dtreshold, bresult)

!<description>
    ! This subroutine checks if the matrix is (strictly/irreducibly)
    ! diagonally dominant depending on the parameter ccheckType. If
    ! the optional parameter bresult is present, then the result is
    ! returned. Otherwise, the outcome of the test is printed.
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! One of the MCHK_DIAGONALDOMINANT_XXX constants
    integer, intent(in) :: ccheckType

    ! OPTIONAL: entries which are smaller than the given treshold
    ! value are not considered in the check.
    real(DP), intent(in), optional :: dtreshold
!</input>

!<output>
    ! OPTIONAL:: If present, the result of the check is returned.
    logical, intent(out), optional :: bresult
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_Da
    real(SP), dimension(:), pointer :: p_Fa
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Kdiagonal
    integer, dimension(:), pointer :: p_IsccIdx, p_Iscc
    integer :: nscc
    logical :: bresult1

    ! Check if matrix is square matrix
    if (rmatrix%NCOLS .ne. rmatrix%NEQ) then
      if (present(bresult)) then
        bresult = .false.
      else
        call output_line('Matrix is not square matrix!',&
            OU_CLASS_MSG,OU_MODE_STD,'mchk_isDiagDominant')
      end if
      return
    end if

    ! What matrix type are we?
    select case(rmatrix%cmatrixFormat)

    case(LSYSSC_MATRIX7)

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      ! What data type are we?
      select case(rmatrix%cdataType)

      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call check_ddMat79DP(rmatrix%NEQ,&
            p_Kld, p_Kcol, p_Kld, p_Da, ccheckType, bresult1)

      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call check_ddMat79SP(rmatrix%NEQ,&
            p_Kld, p_Kcol, p_Kld, p_Fa, ccheckType, bresult1)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mchk_isDiagDominant')
        call sys_halt()
      end select

    case(LSYSSC_MATRIX9)

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

      ! What data type are we?
      select case(rmatrix%cdataType)

      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call check_ddMat79DP(rmatrix%NEQ,&
            p_Kld, p_Kcol, p_Kdiagonal, p_Da, ccheckType, bresult1)

      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call check_ddMat79SP(rmatrix%NEQ,&
            p_Kld, p_Kcol, p_Kdiagonal, p_Fa, ccheckType, bresult1)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mchk_isDiagDominant')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mchk_isDiagDominant')
      call sys_halt()
    end select

    ! Special check if matrix is irreducible
    if (bresult1 .and.&
        ccheckType .eq. MCHK_DIAGDOMINANT_IRREDUCIBLE) then
      call mchk_calcStrongConnComp(rmatrix, .true., nscc, dtreshold, p_IsccIdx, p_Iscc)
      bresult1 = (nscc .eq. 1)
    end if

    ! Report outcome of check?
    if (present(bresult)) then
      bresult = bresult1
    else
      select case(ccheckType)
      case (MCHK_DIAGDOMINANT)
        if (bresult1) then
          call output_line('Matrix is diagonally dominant!',&
              OU_CLASS_MSG,OU_MODE_STD,'mchk_isDiagDominant')
        else
          call output_line('Matrix is not diagonally dominant!',&
              OU_CLASS_MSG,OU_MODE_STD,'mchk_isDiagDominant')
        end if
      case (MCHK_DIAGDOMINANT_STRICT)
        if (bresult1) then
          call output_line('Matrix is strictly diagonally dominant!',&
              OU_CLASS_MSG,OU_MODE_STD,'mchk_isDiagDominant')
        else
          call output_line('Matrix is not strictly diagonally dominant!',&
              OU_CLASS_MSG,OU_MODE_STD,'mchk_isDiagDominant')
        end if
      case (MCHK_DIAGDOMINANT_IRREDUCIBLE)
        if (bresult1) then
          call output_line('Matrix is irreducibly diagonally dominant!',&
              OU_CLASS_MSG,OU_MODE_STD,'mchk_isDiagDominant')
        else
          call output_line('Matrix is not irreducibly diagonally dominant!',&
              OU_CLASS_MSG,OU_MODE_STD,'mchk_isDiagDominant')
        end if
      end select
    end if

  contains

    !***************************************************************************
    ! Diagonal dominance check for double-valued matrix stored in format 7/9

    pure subroutine check_ddMat79DP(NEQ, Kld, Kcol, Kdiagonal, Da,&
                                    ccheckType, bresult)

      real(DP), dimension(:), intent(in) :: Da
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal
      integer, intent(in) :: NEQ, ccheckType
      logical, intent(out) :: bresult

      ! local variables
      real(DP) :: dtmp,ddiag
      integer :: ieq,ia
      logical :: ballDiagDominance,ballStrictDiagDominance
      logical :: banyDiagDominance,banyStrictDiagDominance

      ! Initialisation
      ballDiagDominance = .true.
      banyDiagDominance = .false.
      ballStrictDiagDominance = .true.
      banyStrictDiagDominance = .false.

      !$omp parallel do default(shared) private(dtmp,ddiag,ia)&
      !$omp reduction(.and. : ballDiagDominance,ballStrictDiagDominance)&
      !$omp reduction(.or.  : banyDiagDominance,banyStrictDiagDominance)
      do ieq = 1,NEQ

        ! Compute absolute row sum
        dtmp = 0.0_DP
        do ia = Kld(ieq), Kld(ieq+1)-1
          dtmp = dtmp + abs(Da(ia))
        end do

        ! Compute absolute value of diagonal coefficient
        ddiag = 2.0_DP * abs(Da(Kdiagonal(ieq)))

        ! Check for diagonal dominance
        ballDiagDominance =  ballDiagDominance .and. (ddiag .ge. dtmp)
        banyDiagDominance =  banyDiagDominance .or.  (ddiag .ge. dtmp)

        ! Check for strict diagonal dominance
        ballStrictDiagDominance =  ballStrictDiagDominance .and. (ddiag .gt. dtmp)
        banyStrictDiagDominance =  banyStrictDiagDominance .or.  (ddiag .gt. dtmp)
      end do
      !$omp end parallel do

      ! What check should be applied?
      select case(ccheckType)
      case (MCHK_DIAGDOMINANT)
        bresult = ballDiagDominance
      case (MCHK_DIAGDOMINANT_STRICT)
        bresult = ballStrictDiagDominance
      case (MCHK_DIAGDOMINANT_IRREDUCIBLE)
        bresult = ballDiagDominance .and. banyStrictDiagDominance
      end select
    end subroutine check_ddMat79DP

    !***************************************************************************
    ! Diagonal dominance check for single-valued matrix stored in format 7/9

    pure subroutine check_ddMat79SP(NEQ, Kld, Kcol, Kdiagonal, Fa,&
                                    ccheckType, bresult)

      real(SP), dimension(:), intent(in) :: Fa
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal
      integer, intent(in) :: NEQ, ccheckType
      logical, intent(out) :: bresult

      ! local variables
      real(SP) :: ftmp,fdiag
      integer :: ieq,ia
      logical :: ballDiagDominance,ballStrictDiagDominance
      logical :: banyDiagDominance,banyStrictDiagDominance

      ! Initialisation
      ballDiagDominance = .true.
      banyDiagDominance = .false.
      ballStrictDiagDominance = .true.
      banyStrictDiagDominance = .false.

      !$omp parallel do default(shared) private(ftmp,fdiag,ia)&
      !$omp reduction(.and. : ballDiagDominance,ballStrictDiagDominance)&
      !$omp reduction(.or.  : banyDiagDominance,banyStrictDiagDominance)
      do ieq = 1,NEQ

        ! Compute absolute row sum
        ftmp = 0.0_SP
        do ia = Kld(ieq), Kld(ieq+1)-1
          ftmp = ftmp + abs(Fa(ia))
        end do

        ! Compute absolute value of diagonal coefficient
        fdiag = 2.0_SP * abs(Fa(Kdiagonal(ieq)))

        ! Check for diagonal dominance
        ballDiagDominance =  ballDiagDominance .and. (fdiag .ge. ftmp)
        banyDiagDominance =  banyDiagDominance .or.  (fdiag .ge. ftmp)

        ! Check for strict diagonal dominance
        ballStrictDiagDominance =  ballStrictDiagDominance .and. (fdiag .gt. ftmp)
        banyStrictDiagDominance =  banyStrictDiagDominance .or.  (fdiag .gt. ftmp)
      end do
      !$omp end parallel do

      ! What check should be applied?
      select case(ccheckType)
      case (MCHK_DIAGDOMINANT)
        bresult = ballDiagDominance
      case (MCHK_DIAGDOMINANT_STRICT)
        bresult = ballStrictDiagDominance
      case (MCHK_DIAGDOMINANT_IRREDUCIBLE)
        bresult = ballDiagDominance .and. banyStrictDiagDominance
      end select
    end subroutine check_ddMat79SP  

  end subroutine mchk_isDiagDominant

  ! ***************************************************************************

!<subroutine>

  subroutine mchk_isMMatrix(rmatrix, dtreshold, bresult)

!<description>
    ! This subroutine checks if the matrix satisfies the sufficient
    ! conditions of an M-matrix, that is, all diagonal entries are
    ! strictly positive, all off-diagonal entries are non-positive and
    ! the matrix is irreducibly diagonally dominant.
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! OPTIONAL: entries which are smaller than the given treshold
    ! value are not considered in the check.
    real(DP), intent(in), optional :: dtreshold
!</input>

!<output>
    ! OPTIONAL:: If present, the result of the check is returned.
    logical, intent(out), optional :: bresult
!</output>
!</subroutine>

    ! local variable
    type(t_matrixScalar) :: rmatrixDense
    real(DP), dimension(:), pointer :: p_Da
    real(SP), dimension(:), pointer :: p_Fa
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Kdiagonal
    integer, dimension(:), allocatable :: Iperm
    real(DP) :: dtrhold
    logical :: bresult1

    dtrhold = 0.0_DP
    if (present(dtreshold)) dtrhold = dtreshold

    ! --------------------------------------------------------------------------
    ! Sufficient conditions for square matrices to be M-matrix are as follows:
    !
    ! #1: A is Z-matrix (a_{ij} <= 0 for all j/=i) and
    !     A is strictly diagonally dominant and a_{ii} > 0 for all i
    !
    ! #2: A is Z-matrix (a_{ij} <= 0 for all j/=i) and
    !     A is irreducibly diagonally dominant and a_{ii} > 0 for all i
    !
    ! #3: A is Z-matrix (a_{ij} <= 0 for all j/=i) and
    !     all leading principal minors are positive
    ! --------------------------------------------------------------------------

    ! Check if matrix is square matrix
    bresult1 = rmatrix%NEQ .eq. rmatrix%NCOLS

    if (.not.bresult1) then
      if (present(bresult)) then
        bresult = bresult1
      else
        call output_line('Matrix is not an M-matrix!',&
            OU_CLASS_MSG,OU_MODE_STD,'mchk_isMMatrix')
      end if
      return
    end if

    ! Check if matrix is Z-matrix
    call mchk_isZmatrix(rmatrix, dtrhold, bresult1)

    if (.not.bresult1) then
      if (present(bresult)) then
        bresult = bresult1
      else
        call output_line('Matrix is not an M-matrix!',&
            OU_CLASS_MSG,OU_MODE_STD,'mchk_isMMatrix')
      end if
      return
    end if

    ! Check if all diagonal entries are strictly positive (cf. #1 and #2)

    ! What matrix type are we?
    select case(rmatrix%cmatrixFormat)

    case(LSYSSC_MATRIX7)

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      ! What data type are we?
      select case(rmatrix%cdataType)

      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call check_signMat79DP(rmatrix%NEQ, rmatrix%dscaleFactor,&
            dtrhold, p_Kld, p_Kcol, p_Da, bresult1)

      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call check_signMat79SP(rmatrix%NEQ, real(rmatrix%dscaleFactor,SP),&
            real(dtrhold,SP), p_Kld, p_Kcol, p_Fa, bresult1)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mchk_isMMatrix')
        call sys_halt()
      end select

    case(LSYSSC_MATRIX9)

      ! Set pointers
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

      ! What data type are we?
      select case(rmatrix%cdataType)

      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call check_signMat79DP(rmatrix%NEQ, rmatrix%dscaleFactor,&
            dtrhold, p_Kdiagonal, p_Kcol, p_Da, bresult1)

      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call check_signMat79SP(rmatrix%NEQ, real(rmatrix%dscaleFactor,SP),&
            real(dtrhold,SP), p_Kdiagonal, p_Kcol, p_Fa, bresult1)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mchk_isMMatrix')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mchk_isMMatrix')
      call sys_halt()
    end select

    ! If all diagonal entries are strictky positive then it suffices
    ! to check if matrix A is strictly or irreducibly diagonally
    ! dominant to detect that matrix A is(!) an M-matrix. If this
    ! check fails this does not mean that A is not an M-matrix.

    if (bresult1) then
      ! Check if matrix is strictly diagonally dominant (cf. #1)
      call mchk_isDiagDominant(rmatrix,&
          MCHK_DIAGDOMINANT_STRICT, dtreshold, bresult1)

      if (bresult1) then
        ! Report outcome of check?
        if (present(bresult)) then
          bresult = bresult1
        else
          call output_line('Matrix is an M-matrix!',&
              OU_CLASS_MSG,OU_MODE_STD,'mchk_isMMatrix')
        end if
        return
      end if

      ! Check if matrix is irreducibly diagonally dominant (cf. #2)
      call mchk_isDiagDominant(rmatrix,&
          MCHK_DIAGDOMINANT_IRREDUCIBLE, dtreshold, bresult1)

      if (bresult1) then
        ! Report outcome of check?
        if (present(bresult)) then
          bresult = bresult1
        else
          call output_line('Matrix is an M-matrix!',&
              OU_CLASS_MSG,OU_MODE_STD,'mchk_isMMatrix')
        end if
        return
      end if
    end if

    ! If we end up here, then matrix A does not satisfy the sufficient
    ! condition of an M-matrix but nontheless A can be an M-matrix.
    call output_line('Check did not return a definite answer. Thus, '//&
        'we permform an extensive check (time consuming)!',&
        OU_CLASS_MSG,OU_MODE_STD,'mchk_isMMatrix')

    ! Duplicate matrix and convert it to full matrix
    call lsyssc_duplicateMatrix(rmatrix, rmatrixDense,&
        LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
    call lsyssc_convertMatrix(rmatrixDense, LSYSSC_MATRIX1, .true.)

    ! Allocate memory for permutation vector
    allocate(Iperm(rmatrixDense%NEQ))

    select case(rmatrixDense%cdataType)
    case(ST_DOUBLE)
      call lsyssc_getbase_double(rmatrixDense, p_Da)
      bresult1 = check_InverseMat1DP(rmatrix%NEQ, rmatrix%NCOLS,&
          p_Da, Iperm, dtrhold)

    case(ST_SINGLE)
      call lsyssc_getbase_single(rmatrixDense, p_Fa)
      bresult1 = check_InverseMat1SP(rmatrix%NEQ, rmatrix%NCOLS,&
          p_Fa, Iperm, real(dtrhold,SP))

    case default
      call output_line('Unsupported data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mchk_isMMatrix')
      call sys_halt()
    end select

    ! Release temporal memory
    call lsyssc_releaseMatrix(rmatrixDense)

    ! Deallocate temporal memory
    deallocate(Iperm)

    ! Report outcome of check?
    if (present(bresult)) then
      bresult = bresult1
    else
      if (bresult1) then
        call output_line('Matrix is an M-matrix!',&
            OU_CLASS_MSG,OU_MODE_STD,'mchk_isMMatrix')
      else
        call output_line('Matrix is not an M-matrix!',&
            OU_CLASS_MSG,OU_MODE_STD,'mchk_isMMatrix')
      end if
    end if

  contains

    ! Here, the real working routines start

    !***************************************************************************
    ! Sign check for double-valued matrix stored in format 7/9

    pure subroutine check_signMat79DP(NEQ, dscale, dtreshold, Kdiagonal,&
                                      Kcol, Da, bresult)

      integer, intent(in) :: NEQ
      real(DP), intent(in) :: dscale,dtreshold
      real(DP), dimension(:), intent(in) :: Da
      integer, dimension(:), intent(in) :: Kdiagonal, Kcol
      logical, intent(out) :: bresult

      ! local variables
      integer :: ieq

      ! Initialisation
      bresult = .true.

      if (dscale .eq. 1.0_DP) then
        !$omp parallel do default(shared) reduction(.and. : bresult)
        do ieq = 1,NEQ
          ! Check if diagonal entry is strictly positive
          if (abs(Da(Kdiagonal(ieq))) .le. dtreshold) cycle
          bresult = (bresult .and. Da(Kdiagonal(ieq)) .gt. 0.0_DP)
        end do
        !$omp end parallel do
      else
        !$omp parallel do default(shared) reduction(.and. : bresult)
        do ieq = 1,NEQ
          ! Check if diagonal entry is strictly positive
          if (abs(dscale*Da(Kdiagonal(ieq))) .le. dtreshold) cycle
          bresult = (bresult .and. dscale*Da(Kdiagonal(ieq)) .gt. 0.0_DP)
        end do
        !$omp end parallel do
      end if

    end subroutine check_signMat79DP

    !***************************************************************************
    ! Sign check for single-valued matrix stored in format 7/9

    pure subroutine check_signMat79SP(NEQ, fscale, ftreshold, Kdiagonal,&
                                      Kcol, Fa, bresult)

      integer, intent(in) :: NEQ
      real(SP), intent(in) :: fscale,ftreshold
      real(SP), dimension(:), intent(in) :: Fa
      integer, dimension(:), intent(in) :: Kdiagonal, Kcol
      logical, intent(out) :: bresult

      ! local variables
      integer :: ieq

      ! Initialisation
      bresult = .true.

      if (fscale .eq. 1.0_SP) then
        !$omp parallel do default(shared) reduction(.and. : bresult)
        do ieq = 1,NEQ
          ! Check if diagonal entry is strictly positive
          if (abs(Fa(Kdiagonal(ieq))) .le. ftreshold) cycle
          bresult = (bresult .and. Fa(Kdiagonal(ieq)) .gt. 0.0_SP)
        end do
        !$omp end parallel do
      else
        !$omp parallel do default(shared) reduction(.and. : bresult)
        do ieq = 1,NEQ
          ! Check if diagonal entry is strictly positive
          if (abs(fscale*Fa(Kdiagonal(ieq))) .le. ftreshold) cycle
          bresult = (bresult .and. fscale*Fa(Kdiagonal(ieq)) .gt. 0.0_SP)
        end do
        !$omp end parallel do
      end if

    end subroutine check_signMat79SP

    !***************************************************************************
    ! Compute LU-factorisation of double-valued dense matrix and check
    ! if inverse matrix exists and all of its entries are non-negative

    function check_InverseMat1DP(NEQ, NCOLS, DA, Iperm, dtreshold)&
        result(bresult)

      integer, intent(in) :: NEQ,NCOLS
      real(DP), intent(in) :: dtreshold

      real(DP), dimension(NEQ,NCOLS), intent(inout) :: Da
      integer, dimension(:), intent(inout) :: Iperm

      logical :: bresult

      ! local variables
      real(DP), dimension(:), allocatable :: Db,Dx
      integer :: iaux,icol,ieq,ipiv,jeq,iuv
      real(DP) :: daux


      ! Initialisation
      bresult = .true.

      ! Compute the LU-decomposition of matrix using Gaussian
      ! elimination with partial pivoting

      do ieq = 1,NEQ
        Iperm(ieq) = ieq
      end do

      ! Perform Gaussian elimination
      row: do ieq = 1,NEQ

        ! Row-wise pivoting
        ipiv = ieq

        do jeq = ieq,NEQ
          if (abs(Da(Iperm(jeq),ieq)) .gt. abs(Da(Iperm(ipiv),ieq))) then
            ipiv = jeq
          end if
        end do

        ! Swap rows IPIV <-> IEQ
        if (Iperm(ipiv) .ne. Iperm(ieq)) then
          iaux        = Iperm(ieq)
          Iperm(ieq)  = Iperm(ipiv)
          Iperm(ipiv) = iaux
        end if

        ! Get pivot element
        daux = Da(Iperm(ieq),ieq)

        ! Check for zero determinant
        if (abs(daux) .le. SYS_EPSREAL_DP) then
          ! Matrix is not invertible
          bresult = .false.
          return
        end if

        ! Elimination step
        col: do jeq = ieq+1,NEQ
          Da(Iperm(jeq),ieq) = Da(Iperm(jeq),ieq) / daux
          do icol = ieq+1,NCOLS
            Da(Iperm(jeq),icol) = Da(Iperm(jeq),icol)&
                                - Da(Iperm(ieq),icol)*Da(Iperm(jeq),ieq)
          end do
        end do col
      end do row

      ! Allocate temporal memory
      allocate(Db(NEQ),Dx(NCOLS))

      ! Loop over all unit vectors and compute inverse matrix
      ! column-by-column
      do iuv = 1,NEQ

        ! Initialise i-th unit vector
        Db = 0.0_DP
        Db(iuv) = 1.0_DP

        ! Forward substitution
        do ieq = 1,NEQ-1
          do jeq = ieq+1,NEQ
            Db(Iperm(jeq)) = Db(Iperm(jeq))-Db(Iperm(ieq))*Da(Iperm(jeq),ieq)
          end do
        end do

        ! Initialise solution
        Dx = 0.0_DP

        ! Backward substitution
        do ieq = NEQ,1,-1
          Dx(Iperm(ieq)) = Db(Iperm(ieq))
          do jeq = ieq+1,NEQ
            Dx(Iperm(ieq)) = Dx(Iperm(ieq)) - Da(Iperm(ieq),jeq)*Dx(Iperm(jeq))
          end do
          Dx(Iperm(ieq)) = Dx(Iperm(ieq)) / Da(Iperm(ieq),ieq)

          ! Check if entry is negative
          if (Dx(Iperm(ieq)) .lt. -dtreshold) then
            bresult = .false.
            deallocate(Db,Dx)
            return
          end if
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Db,Dx)

    end function check_InverseMat1DP

    !***************************************************************************
    ! Compute LU-factorisation of single-valued dense matrix and check
    ! if inverse matrix exists and all of its entries are non-negative

    function check_InverseMat1SP(NEQ, NCOLS, FA, Iperm, ftreshold)&
        result(bresult)

      integer, intent(in) :: NEQ,NCOLS
      real(SP), intent(in) :: ftreshold

      real(SP), dimension(NEQ,NCOLS), intent(inout) :: Fa
      integer, dimension(:), intent(inout) :: Iperm

      logical :: bresult

      ! local variables
      real(SP), dimension(:), allocatable :: Fb,Fx
      integer :: iaux,icol,ieq,ipiv,jeq,iuv
      real(SP) :: faux


      ! Initialisation
      bresult = .true.

      ! Compute the LU-decomposition of matrix using Gaussian
      ! elimination with partial pivoting

      do ieq = 1,NEQ
        Iperm(ieq) = ieq
      end do

      ! Perform Gaussian elimination
      row: do ieq = 1,NEQ

        ! Row-wise pivoting
        ipiv = ieq

        do jeq = ieq,NEQ
          if (abs(Fa(Iperm(jeq),ieq)) .gt. abs(Fa(Iperm(ipiv),ieq))) then
            ipiv = jeq
          end if
        end do

        ! Swap rows IPIV <-> IEQ
        if (Iperm(ipiv) .ne. Iperm(ieq)) then
          iaux        = Iperm(ieq)
          Iperm(ieq)  = Iperm(ipiv)
          Iperm(ipiv) = iaux
        end if

        ! Get pivot element
        faux = Fa(Iperm(ieq),ieq)

        ! Check for zero determinant
        if (abs(faux) .le. SYS_EPSREAL_SP) then
          ! Matrix is not invertible
          bresult = .false.
          return
        end if

        ! Elimination step
        col: do jeq = ieq+1,NEQ
          Fa(Iperm(jeq),ieq) = Fa(Iperm(jeq),ieq) / faux
          do icol = ieq+1,NCOLS
            Fa(Iperm(jeq),icol) = Fa(Iperm(jeq),icol)&
                                - Fa(Iperm(ieq),icol)*Fa(Iperm(jeq),ieq)
          end do
        end do col
      end do row

      ! Allocate temporal memory
      allocate(Fb(NEQ),Fx(NCOLS))

      ! Loop over all unit vectors and compute inverse matrix
      ! column-by-column
      do iuv = 1,NEQ

        ! Initialise i-th unit vector
        Fb = 0.0_SP
        Fb(iuv) = 1.0_SP

        ! Forward substitution
        do ieq = 1,NEQ-1
          do jeq = ieq+1,NEQ
            Fb(Iperm(jeq)) = Fb(Iperm(jeq))-Fb(Iperm(ieq))*Fa(Iperm(jeq),ieq)
          end do
        end do

        ! Initialise solution
        Fx = 0.0_SP

        ! Backward substitution
        do ieq = NEQ,1,-1
          Fx(Iperm(ieq)) = Fb(Iperm(ieq))
          do jeq = ieq+1,NEQ
            Fx(Iperm(ieq)) = Fx(Iperm(ieq)) - Fa(Iperm(ieq),jeq)*Fx(Iperm(jeq))
          end do
          Fx(Iperm(ieq)) = Fx(Iperm(ieq)) / Fa(Iperm(ieq),ieq)

          ! Check if entry is negative
          if (Fx(Iperm(ieq)) .lt. -dtreshold) then
            bresult = .false.
            deallocate(Fb,Fx)
            return
          end if
        end do
      end do

      ! Deallocate temporal memory
      deallocate(Fb,Fx)

    end function check_InverseMat1SP

  end subroutine mchk_isMMatrix

  ! ***************************************************************************

!<subroutine>

  subroutine mchk_isZMatrix(rmatrix, dtreshold, bresult)

!<description>
    ! This subroutine checks if the matrix is a Z-matrix.
    ! That is, all off-diagonal entries are non-positive.
    !
    ! $$ A=(a_{ij}),\quad a_{ij}\le 0\, \forall j\ne i $$
    !
    ! If the optional parameter dtreshold is present, then the above
    ! condition is relaxed as follows:
    !
    ! $$ A=(a_{ij}),\quad a_{ij}\le -|dtreshold|\, \forall j\ne i $$
    !
    ! If the optional parameter bresult is present, then the result
    ! is returned. Otherwise, the outcome of the test is printed.
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! OPTIONAL: entries which are smaller than the given treshold
    ! value are not considered in the check.
    real(DP), intent(in), optional :: dtreshold
!</input>

!<output>
    ! OPTIONAL:: If present, the result of the check is returned.
    logical, intent(out), optional :: bresult
!</output>
!</subroutine>

    ! local variable
    real(DP), dimension(:), pointer :: p_Da
    real(SP), dimension(:), pointer :: p_Fa
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Kdiagonal
    real(DP) :: dtrhold
    logical :: bresult1

    dtrhold = 0.0_DP
    if (present(dtreshold)) dtrhold = dtreshold

    ! Check if all off-diagonal entries are non-positive

    ! What matrix type are we?
    select case(rmatrix%cmatrixFormat)

    case(LSYSSC_MATRIX7)

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      ! What data type are we?
      select case(rmatrix%cdataType)

      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call check_signMat7DP(rmatrix%NEQ, rmatrix%dscaleFactor,&
            dtrhold, p_Kld, p_Kcol, p_Da, bresult1)

      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call check_signMat7SP(rmatrix%NEQ, real(rmatrix%dscaleFactor,SP),&
            real(dtrhold,SP), p_Kld, p_Kcol, p_Fa, bresult1)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mchk_isMMatrix')
        call sys_halt()
      end select

    case(LSYSSC_MATRIX9)

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

      ! What data type are we?
      select case(rmatrix%cdataType)

      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call check_signMat9DP(rmatrix%NEQ, rmatrix%dscaleFactor,&
            dtrhold, p_Kld, p_Kcol, p_Kdiagonal, p_Da, bresult1)

      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call check_signMat9SP(rmatrix%NEQ, real(rmatrix%dscaleFactor,SP),&
            real(dtrhold,SP),p_Kld, p_Kcol, p_Kdiagonal, p_Fa, bresult1)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mchk_isZMatrix')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mchk_isZMatrix')
      call sys_halt()
    end select

    ! Report outcome of check?
    if (present(bresult)) then
      bresult = bresult1
    else
      if (bresult1) then
        call output_line('Matrix is a Z-matrix!',&
            OU_CLASS_MSG,OU_MODE_STD,'mchk_isZMatrix')
      else
        call output_line('Matrix is not a Z-matrix!',&
            OU_CLASS_MSG,OU_MODE_STD,'mchk_isZMatrix')
      end if
    end if

  contains

    ! Here, the real working routines start

    !***************************************************************************
    ! Sign check for double-valued matrix stored in format 7

    pure subroutine check_signMat7DP(NEQ, dscale, dtreshold, Kld, Kcol,&
                                     Da, bresult)

      integer, intent(in) :: NEQ
      real(DP), intent(in) :: dscale,dtreshold
      real(DP), dimension(:), intent(in) :: Da
      integer, dimension(:), intent(in) :: Kld, Kcol
      logical, intent(out) :: bresult

      ! local variables
      integer :: ieq,ia

      ! Initialisation
      bresult = .true.

      if (dscale .eq. 1.0_DP) then
        !$omp parallel do default(shared) private(ia) reduction(.and. : bresult)
        do ieq = 1,NEQ
          ! Check if off-diagonal entries are non-positive
          do ia = Kld(ieq)+1, Kld(ieq+1)-1
            if (abs(Da(ia)) .le. dtreshold) cycle
            bresult = (bresult .and. Da(ia) .le. 0.0_DP)
          end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do default(shared) private(ia) reduction(.and. : bresult)
        do ieq = 1,NEQ
          ! Check if off-diagonal entries are non-positive
          do ia = Kld(ieq)+1, Kld(ieq+1)-1
            if (abs(dscale*Da(ia)) .le. dtreshold) cycle
            bresult = (bresult .and. dscale*Da(ia) .le. 0.0_DP)
          end do
        end do
        !$omp end parallel do
      end if

    end subroutine check_signMat7DP

    !***************************************************************************
    ! Sign check for single-valued matrix stored in format 7

    pure subroutine check_signMat7SP(NEQ, fscale, ftreshold, Kld, Kcol,&
                                     Fa, bresult)

      integer, intent(in) :: NEQ
      real(SP), intent(in) :: fscale,ftreshold
      real(SP), dimension(:), intent(in) :: Fa
      integer, dimension(:), intent(in) :: Kld, Kcol
      logical, intent(out) :: bresult

      ! local variables
      integer :: ieq,ia

      ! Initialisation
      bresult = .true.

      if (fscale .eq. 1.0_SP) then
        !$omp parallel do default(shared) private(ia) reduction(.and. : bresult)
        do ieq = 1,NEQ
          ! Check if off-diagonal entries are non-positive
          do ia = Kld(ieq)+1, Kld(ieq+1)-1
            if (abs(Fa(ia)) .le. ftreshold) cycle
            bresult = (bresult .and. Fa(ia) .le. 0.0_SP)
          end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do default(shared) private(ia) reduction(.and. : bresult)
        do ieq = 1,NEQ
          ! Check if off-diagonal entries are non-positive
          do ia = Kld(ieq)+1, Kld(ieq+1)-1
            if (abs(fscale*Fa(ia)) .le. ftreshold) cycle
            bresult = (bresult .and. fscale*Fa(ia) .le. 0.0_SP)
          end do
        end do
        !$omp end parallel do
      end if

    end subroutine check_signMat7SP

    !***************************************************************************
    ! Sign check for double-valued matrix stored in format 9

    pure subroutine check_signMat9DP(NEQ, dscale, dtreshold, Kld, Kcol,&
                                     Kdiagonal, Da, bresult)

      integer, intent(in) :: NEQ
      real(DP), intent(in) :: dscale,dtreshold
      real(DP), dimension(:), intent(in) :: Da
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal
      logical, intent(out) :: bresult

      ! local variables
      integer :: ieq,ia

      ! Initialisation
      bresult = .true.

      if (dscale .eq. 1.0_DP) then
        !$omp parallel do default(shared) private(ia) reduction(.and. : bresult)
        do ieq = 1,NEQ
          ! Check if off-diagonal entries are non-positive
          do ia = Kld(ieq), Kdiagonal(ieq)-1
            if (abs(Da(ia)) .le. dtreshold) cycle
            bresult = (bresult .and. Da(ia) .le. 0.0_DP)
          end do
          do ia = Kdiagonal(ieq)+1, Kld(ieq+1)-1
            if (abs(Da(ia)) .le. dtreshold) cycle
            bresult = (bresult .and. Da(ia) .le. 0.0_DP)
          end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do default(shared) private(ia) reduction(.and. : bresult)
        do ieq = 1,NEQ
          ! Check if off-diagonal entries are non-positive
          do ia = Kld(ieq), Kdiagonal(ieq)-1
            if (abs(dscale*Da(ia)) .le. dtreshold) cycle
            bresult = (bresult .and. dscale*Da(ia) .le. 0.0_DP)
          end do
          do ia = Kdiagonal(ieq)+1, Kld(ieq+1)-1
            if (abs(dscale*Da(ia)) .le. dtreshold) cycle
            bresult = (bresult .and. dscale*Da(ia) .le. 0.0_DP)
          end do
        end do
        !$omp end parallel do
      end if

    end subroutine check_signMat9DP

    !***************************************************************************
    ! Sign check for single-valued matrix stored in format 9

    pure subroutine check_signMat9SP(NEQ, fscale, ftreshold, Kld, Kcol,&
                                     Kdiagonal, Fa, bresult)

      integer, intent(in) :: NEQ
      real(SP), intent(in) :: fscale,ftreshold
      real(SP), dimension(:), intent(in) :: Fa
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal
      logical, intent(out) :: bresult

      ! local variables
      integer :: ieq,ia

      ! Initialisation
      bresult = .true.

      if (fscale .eq. 1.0_SP) then
        !$omp parallel do default(shared) private(ia) reduction(.and. : bresult)
        do ieq = 1,NEQ
          ! Check if off-diagonal entries are non-positive
          do ia = Kld(ieq), Kdiagonal(ieq)-1
            if (abs(Fa(ia)) .le. ftreshold) cycle
            bresult = (bresult .and. Fa(ia) .le. 0.0_SP)
          end do
          do ia = Kdiagonal(ieq)+1, Kld(ieq+1)-1
            if (abs(Fa(ia)) .le. ftreshold) cycle
            bresult = (bresult .and. Fa(ia) .le. 0.0_SP)
          end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do default(shared) private(ia) reduction(.and. : bresult)
        do ieq = 1,NEQ
          ! Check if off-diagonal entries are non-positive
          do ia = Kld(ieq), Kdiagonal(ieq)-1
            if (abs(fscale*Fa(ia)) .le. ftreshold) cycle
            bresult = (bresult .and. fscale*Fa(ia) .le. 0.0_SP)
          end do
          do ia = Kdiagonal(ieq)+1, Kld(ieq+1)-1
            if (abs(fscale*Fa(ia)) .le. ftreshold) cycle
            bresult = (bresult .and. fscale*Fa(ia) .le. 0.0_SP)
          end do
        end do
        !$omp end parallel do
      end if

    end subroutine check_signMat9SP

  end subroutine mchk_isZMatrix

  ! ***************************************************************************

!<subroutine>

  subroutine mchk_calcStrongConnComp(rmatrix, bdata, nscc, dtreshold,&
                                     p_IsccIdx, p_Iscc)

!<description>
    ! This subroutine calculates the strongly connected components of
    ! the graph that corresponds to to matrix rmatrix. This algorithm
    ! is known as Tarjan's arlgorithm.
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! Whether to consider the real data of the matrix or its sparsity
    ! pattern
    logical, intent(in) :: bdata

    ! OPTIONAL: If bdata=true, then matrix coefficients with absolute
    ! value smaller than the given treshold are considered zero
    real(DP), intent(in), optional :: dtreshold
!</input>

!<output>
    ! Number of strongly connected components
    integer, intent(out) :: nscc

    ! OPTIONAL: Pointer to integer array storing the strongly
    ! connected components. If not present, then this routine returns
    ! only the number of strongly connected components (nscc).
    integer, dimension(:), pointer, optional :: p_IsccIdx, p_Iscc
!</output>
!</subroutine>

    ! locale variables
    type(t_stackInt) :: rstack
    integer, dimension(:), pointer :: p_Kld,p_Kcol
    integer, dimension(:), allocatable :: Index,LowLink
    real(DP), dimension(:), pointer :: p_Da
    real(SP), dimension(:), pointer :: p_Fa
    real(DP) :: dtrhold
    integer :: ieq,idx

    dtrhold = 0.0_DP
    if (present(dtreshold)) dtrhold = dtreshold

    ! Initialisation
    call stack_create(rstack, rmatrix%NEQ)

    allocate(index(rmatrix%NEQ), LowLink(rmatrix%NEQ))
    call lalg_clearVector(Index)
    idx  = 1
    nscc = 0

    ! What matrix type are we?
    select case(rmatrix%cmatrixFormat)

    case(LSYSSC_MATRIX7, LSYSSC_MATRIX9)

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      if (bdata) then

        ! What data type are we?
        select case(rmatrix%cdataType)

        case (ST_DOUBLE)
          call lsyssc_getbase_double(rmatrix, p_Da)

          ! Tarjan`s algorithm based on the double-valued matrix
          ! coefficients
          do ieq = 1, rmatrix%NEQ
            if (index(ieq) .eq. 0)&
                call strongconnectMat79DP(rstack, p_Kld, p_Kcol, p_Da, ieq,&
                                          Index, LowLink, idx, nscc)
          end do

          ! Generate list of strongly connected components?
          if (present(p_IsccIdx) .and. present(p_Iscc)) then
            allocate(p_IsccIdx(nscc+1), p_Iscc(rmatrix%NEQ))

            call stack_clear(rstack)
            call lalg_clearVector(Index)
            idx  = 1
            nscc = 0
            p_IsccIdx(1) = 1

            ! Tarjan`s algorithm based on the double-valued matrix
            ! coefficients
            do ieq = 1, rmatrix%NEQ
              if (index(ieq) .eq. 0)&
                  call strongconnectMat79DP(rstack, p_Kld, p_Kcol, p_Da, ieq,&
                                            Index, LowLink, idx, nscc,&
                                            p_IsccIdx, p_Iscc)
            end do
          end if

        case (ST_SINGLE)
          call lsyssc_getbase_single(rmatrix, p_Fa)

          ! Tarjan`s algorithm based on the double-valued matrix
          ! coefficients
          do ieq = 1, rmatrix%NEQ
            if (index(ieq) .eq. 0)&
                call strongconnectMat79SP(rstack, p_Kld, p_Kcol, p_Fa, ieq,&
                                          Index, LowLink, idx, nscc)
          end do

          ! Generate list of strongly connected components?
          if (present(p_IsccIdx) .and. present(p_Iscc)) then
            allocate(p_IsccIdx(nscc+1), p_Iscc(rmatrix%NEQ))

            call stack_clear(rstack)
            call lalg_clearVector(Index)
            idx  = 1
            nscc = 0
            p_IsccIdx(1) = 1

            ! Tarjan`s algorithm based on the double-valued matrix
            ! coefficients
            do ieq = 1, rmatrix%NEQ
              if (index(ieq) .eq. 0)&
                  call strongconnectMat79SP(rstack, p_Kld, p_Kcol, p_Fa, ieq,&
                                            Index, LowLink, idx, nscc,&
                                            p_IsccIdx, p_Iscc)
            end do
          end if

        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mchk_calcStrongConnComp')
          call sys_halt()
        end select

      else

        ! Tarjan`s algorithm based on the sparsity graph
        do ieq = 1, rmatrix%NEQ
          if (index(ieq) .eq. 0)&
              call strongconnectMat79(rstack, p_Kld, p_Kcol, ieq,&
                                      Index, LowLink, idx, nscc)
        end do

        ! Generate list of strongly connected components?
        if (present(p_IsccIdx) .and. present(p_Iscc)) then
          allocate(p_IsccIdx(nscc+1), p_Iscc(rmatrix%NEQ))

          call stack_clear(rstack)
          call lalg_clearVector(Index)
          idx  = 1
          nscc = 0
          p_IsccIdx(1) = 1

          ! Tarjan`s algorithm based on the sparsity graph
          do ieq = 1, rmatrix%NEQ
            if (index(ieq) .eq. 0)&
                call strongconnectMat79(rstack, p_Kld, p_Kcol, ieq,&
                                        Index, LowLink, idx, nscc,&
                                        p_IsccIdx, p_Iscc)
          end do
        end if

      end if

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mchk_calcStrongConnComp')
      call sys_halt()
    end select

    ! Finalisation
    call stack_release(rstack)
    deallocate(Index,LowLink)

  contains

    ! Here, the real working routine follows

    !***************************************************************************
    ! Compute SCC based on the sparsity graph for matrix format 7 and 9

    recursive subroutine strongconnectMat79(rstack, Kld, Kcol, ieq, Index,&
                                            LowLink, idx, nscc, IsccIdx, Iscc)

      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: ieq

      type(t_stackInt), intent(inout) :: rstack
      integer, dimension(:), intent(inout) :: Index,LowLink
      integer, intent(inout) :: idx,nscc

      integer, dimension(:), optional :: IsccIdx, Iscc

      ! local variables
      integer :: ia,jeq,iidx

      ! Set the depth index for IEQ to the smallest unused index
      index(ieq)   = idx
      LowLink(ieq) = idx
      idx = idx+1
      call stack_push(rstack, ieq)

      ! Consider successor of IEQ
      do ia = Kld(ieq), Kld(ieq+1)-1
        jeq = Kcol(ia)

        if (index(jeq) .eq. 0) then
          ! Successor JEQ has not yet been visited; recurse on it
          call strongconnectMat79(rstack, Kld, Kcol, jeq, Index,&
                                  LowLink, idx, nscc, IsccIdx, Iscc)

          LowLink(ieq) = min(LowLink(ieq), LowLink(jeq))
        elseif (stack_contains(rstack, jeq)) then
          ! Successor JEQ is in the stack and hence in the current
          ! strongly connected component
          LowLink(ieq) = min(LowLink(ieq), index(jeq))
        end if
      end do

      ! If IEQ is a root node, pop the stack and generate a new
      ! strongly connected component; two versions!

      if (present(IsccIdx) .and. present(Iscc)) then

        if (LowLink(ieq) .eq. index(ieq)) then
          nscc = nscc+1
          iidx = IsccIdx(nscc)
          scc1: do
            call stack_pop(rstack, jeq)
            Iscc(iidx) = jeq
            iidx = iidx+1
            if (jeq .eq. ieq) exit scc1
          end do scc1
          IsccIdx(nscc+1) = iidx
        end if

      else

        if (LowLink(ieq) .eq. index(ieq)) then
          nscc = nscc+1
          scc2: do
            call stack_pop(rstack, jeq)
            if (jeq .eq. ieq) exit scc2
          end do scc2
        end if

      end if

    end subroutine strongconnectMat79

    !***************************************************************************
    ! Compute SCC based on double-based data for matrix format 7 and 9

    recursive subroutine strongconnectMat79DP(rstack, Kld, Kcol, Da, ieq, Index,&
                                              LowLink, idx, nscc, IsccIdx, Iscc)

      real(DP), dimension(:), intent(in) :: Da
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: ieq

      type(t_stackInt), intent(inout) :: rstack
      integer, dimension(:), intent(inout) :: Index,LowLink
      integer, intent(inout) :: idx,nscc

      integer, dimension(:), optional :: IsccIdx, Iscc

      ! local variables
      integer :: ia,jeq,iidx

      ! Set the depth index for IEQ to the smallest unused index
      index(ieq)   = idx
      LowLink(ieq) = idx
      idx = idx+1
      call stack_push(rstack, ieq)

      ! Consider successor of IEQ
      do ia = Kld(ieq), Kld(ieq+1)-1
        if (abs(Da(ia)) .le. dtrhold) cycle
        jeq = Kcol(ia)

        if (index(jeq) .eq. 0) then
          ! Successor JEQ has not yet been visited; recurse on it
          call strongconnectMat79DP(rstack, Kld, Kcol, Da, jeq, Index,&
                                    LowLink, idx, nscc, IsccIdx, Iscc)

          LowLink(ieq) = min(LowLink(ieq), LowLink(jeq))
        elseif (stack_contains(rstack, jeq)) then
          ! Successor JEQ is in the stack and hence in the current
          ! strongly connected component
          LowLink(ieq) = min(LowLink(ieq), index(jeq))
        end if
      end do

      ! If IEQ is a root node, pop the stack and generate a new
      ! strongly connected component; two versions!

      if (present(IsccIdx) .and. present(Iscc)) then

        if (LowLink(ieq) .eq. index(ieq)) then
          nscc = nscc+1
          iidx = IsccIdx(nscc)
          scc1: do
            call stack_pop(rstack, jeq)
            Iscc(iidx) = jeq
            iidx = iidx+1
            if (jeq .eq. ieq) exit scc1
          end do scc1
          IsccIdx(nscc+1) = iidx
        end if

      else

        if (LowLink(ieq) .eq. index(ieq)) then
          nscc = nscc+1
          scc2: do
            call stack_pop(rstack, jeq)
            if (jeq .eq. ieq) exit scc2
          end do scc2
        end if

      end if

    end subroutine strongconnectMat79DP

    !***************************************************************************
    ! Compute SCC based on single-based data for matrix format 7 and 9

    recursive subroutine strongconnectMat79SP(rstack, Kld, Kcol, Fa, ieq, Index,&
                                              LowLink, idx, nscc, IsccIdx, Iscc)

      real(SP), dimension(:), intent(in) :: Fa
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: ieq

      type(t_stackInt), intent(inout) :: rstack
      integer, dimension(:), intent(inout) :: Index,LowLink
      integer, intent(inout) :: idx,nscc

      integer, dimension(:), optional :: IsccIdx, Iscc

      ! local variables
      integer :: ia,jeq,iidx

      ! Set the depth index for IEQ to the smallest unused index
      index(ieq)   = idx
      LowLink(ieq) = idx
      idx = idx+1
      call stack_push(rstack, ieq)

      ! Consider successor of IEQ
      do ia = Kld(ieq), Kld(ieq+1)-1
        if (abs(Fa(ia)) .le. real(dtrhold,SP)) cycle
        jeq = Kcol(ia)

        if (index(jeq) .eq. 0) then
          ! Successor JEQ has not yet been visited; recurse on it
          call strongconnectMat79SP(rstack, Kld, Kcol, Fa, jeq, Index,&
                                    LowLink, idx, nscc, IsccIdx, Iscc)

          LowLink(ieq) = min(LowLink(ieq), LowLink(jeq))
        elseif (stack_contains(rstack, jeq)) then
          ! Successor JEQ is in the stack and hence in the current
          ! strongly connected component
          LowLink(ieq) = min(LowLink(ieq), index(jeq))
        end if
      end do

      ! If IEQ is a root node, pop the stack and generate a new
      ! strongly connected component; two versions!

      if (present(IsccIdx) .and. present(Iscc)) then

        if (LowLink(ieq) .eq. index(ieq)) then
          nscc = nscc+1
          iidx = IsccIdx(nscc)
          scc1: do
            call stack_pop(rstack, jeq)
            Iscc(iidx) = jeq
            iidx = iidx+1
            if (jeq .eq. ieq) exit scc1
          end do scc1
          IsccIdx(nscc+1) = iidx
        end if

      else

        if (LowLink(ieq) .eq. index(ieq)) then
          nscc = nscc+1
          scc2: do
            call stack_pop(rstack, jeq)
            if (jeq .eq. ieq) exit scc2
          end do scc2
        end if

      end if

    end subroutine strongconnectMat79SP

  end subroutine mchk_calcStrongConnComp

end module matrixcheck


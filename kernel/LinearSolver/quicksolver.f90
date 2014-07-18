!##############################################################################
!# ****************************************************************************
!# <name> quicksolver </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a hand full of 'fire-and-forget' linear solvers.
!# The routines in this module are all one-call solvers which do not need
!# any initialisation and do not allocate memory.
!# The solvers implemented in this module are not meant as an alternative to
!# the ones in the linearsolver module, they are meant as 'low-level' solvers
!# for other kernel modules which cannot use the linearsolver module for
!# whatever reason.
!#
!# As the routines do not allocate memory, the caller is responsible to pass
!# a working array to the solver which is used to store temporary vectors.
!#
!# The following routines can be found in this module:
!#
!#  1.) qsol_solveCG
!#      -> Solves a linear system using an unpreconditioned CG solver.
!#         Supported Matrix types: Format-7 / Format-9 CSR
!#         Working Memory: 3*n
!#
!#  2.) qsol_solveCG_SSOR
!#      -> Solves a linear system using a CG solver with built-in SSOR
!#         preconditioner.
!#         Supported Matrix types: Format-9 CSR
!#         Working memory: 3*n
!#
!#  3.) qsol_solveSOR
!#      -> Solves a linear system using a SOR solver.
!#         Supported Matrix types: Format-9 CSR
!#         Working memory: n
!#
!#  4.) qsol_precSOR
!#      -> Applies the SOR preconditioner onto a defect vector.
!#         Supported Matrix types: Format-9 CSR
!#         Working memory: none
!#
!#  5.) qsol_precSSOR
!#      -> Applies the SSOR preconditioner onto a defect vector.
!#         Supported Matrix types: Format-9 CSR
!#         Working memory: none
!#
!# 20.) qsol_solveDiagSchurComp
!#      -> Solves a 2x2 block system with diagonal diagonal blocks
!#         using a Schur complement approach.
!#
!# 21.) qsol_solveTridiag
!#      -> Solves a general tridiagonal systems with an extended Thomas
!#         algorithm.
!#
!# </purpose>
!##############################################################################

module quicksolver

!$ use omp_lib
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemscalar
  use mprimitives
  use perfconfig

  implicit none

  private

  public :: qsol_solveCG
  public :: qsol_solveCG_lsyssc
  public :: qsol_solveCG_double
  public :: qsol_solveCG_SSOR
  public :: qsol_solveCG_SSOR_lsyssc
  public :: qsol_solveCG_SSOR_double
  public :: qsol_solveSOR
  public :: qsol_solveSOR_lsyssc
  public :: qsol_solveSOR_double
  public :: qsol_precSOR
  public :: qsol_precSOR_lsyssc
  public :: qsol_precSOR_double
  public :: qsol_precSSOR
  public :: qsol_precSSOR_lsyssc
  public :: qsol_precSSOR_double

  public :: qsol_solveDiagSchurComp
  public :: qsol_solveTridiag

!<constants>

!<constantblock>
  ! Operation completed successfully
  integer, parameter, public :: QSOL_INFO_SUCCESS           = 0

  ! Maximum iterations reached
  integer, parameter, public :: QSOL_INFO_MAX_ITER          = 1

  ! Internal solver error
  integer, parameter, public :: QSOL_INFO_INTERNAL_ERROR    = 2
!</constantblock>

!</constants>

!************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: qsol_perfconfig

!************************************************************************

  interface qsol_solveCG
    module procedure qsol_solveCG_lsyssc
    module procedure qsol_solveCG_double
  end interface

  interface qsol_solveCG_SSOR
    module procedure qsol_solveCG_SSOR_lsyssc
    module procedure qsol_solveCG_SSOR_double
  end interface

  interface qsol_solveSOR
    module procedure qsol_solveSOR_lsyssc
    module procedure qsol_solveSOR_double
  end interface

  interface qsol_precSOR
    module procedure qsol_precSOR_lsyssc
    module procedure qsol_precSOR_double
  end interface

  interface qsol_precSSOR
    module procedure qsol_precSSOR_lsyssc
    module procedure qsol_precSSOR_double
  end interface

contains

  ! ***************************************************************************

!<subroutine>

  subroutine qsol_solveCG_lsyssc(rmat, rvec, Dwork, cinfo, niter, dtol)

!<description>
  ! This routine is a one-call CG solver.
!</description>

!<input>
  ! The system matrix.
  type(t_matrixScalar), intent(in) :: rmat
!</input>

!<output>
  ! One of the QSOL_INFO_XXXX constants describing the result.
  integer, intent(out) :: cinfo
!</output>

!<inputoutput>
  ! On entry, the vector containing the right-hand-side of the linear system.
  ! On exit, the vector containing the solution of the linear system.
  type(t_vectorScalar), intent(inout) :: rvec

  ! Work array. Its length must be at least 3*n.
  real(DP), dimension(:), target, intent(inout) :: Dwork

  ! On entry, the maximum number of allowed iterations. Must be > 0.
  ! On exit, the total number of performed iterations.
  integer, intent(inout) :: niter

  ! On entry, the absolute tolerance in euclidian norm that is to be reached.
  !   If set to a value <= 0, no residual checks are performed, and the solver
  !   always performs niter iterations, unless an error occurs.
  ! On exit,
  !   if dtol was  > 0 on entry, the absolute final defect,
  !   if dtol was <= 0 on entry, dtol is left unchanged.
  real(DP), intent(inout) :: dtol
!</inputoutput>

!</suboutine>

  ! local variables
  integer, dimension(:), pointer :: p_Kld, p_Kcol
  real(DP), dimension(:), pointer :: p_Da, p_Dx

    ! Get the arrays
    call lsyssc_getbase_Kld(rmat, p_Kld)
    call lsyssc_getbase_Kcol(rmat, p_Kcol)

    call lsyssc_getbase_double(rmat, p_Da)
    call lsyssc_getbase_double(rvec, p_Dx)

    ! Call double version
    call qsol_solveCG_double(rmat%NEQ, p_Kld, p_Kcol, p_Da, p_Dx, Dwork, &
                             cinfo, niter, dtol)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine qsol_solveCG_double(n, Kld, Kcol, Da, Dx, Dwork, cinfo, niter, dtol)

!<description>
  ! This routine is a one-call CG solver.
!</description>

!<input>
  ! The size of the linear system.
  integer, intent(in) :: n

  ! The Kld array of the matrix.
  integer, dimension(:), intent(in) :: Kld

  ! The Kcol array of the matrix.
  integer, dimension(:), intent(in) :: Kcol

  ! The data array of the matrix.
  real(DP), dimension(:), intent(in) :: Da
!</input>

!<output>
  ! One of the QSOL_INFO_XXXX constants describing the result.
  integer, intent(out) :: cinfo
!</output>

!<inputoutput>
  ! On entry, the vector containing the right-hand-side of the linear system.
  ! On exit, the vector containing the solution of the linear system.
  real(DP), dimension(:), intent(inout) :: Dx

  ! Work array. Its length must be at least 3*n.
  real(DP), dimension(:), target, intent(inout) :: Dwork

  ! On entry, the maximum number of allowed iterations. Must be > 0.
  ! On exit, the total number of performed iterations.
  integer, intent(inout) :: niter

  ! On entry, the absolute tolerance in euclidian norm that is to be reached.
  !   If set to a value <= 0, no residual checks are performed, and the solver
  !   always performs niter iterations, unless an error occurs.
  ! On exit,
  !   if dtol was  > 0 on entry, the absolute final defect,
  !   if dtol was <= 0 on entry, dtol is left unchanged.
  real(DP), intent(inout) :: dtol
!</inputoutput>

!</suboutine>

  ! local CG variables
  real(DP) :: dalpha, dbeta, dgamma, dtol2
  integer :: ite

  ! temporary sub-vectors
  real(DP), dimension(:), pointer :: p_Ddef, p_Ddir, p_Dtmp

    ! Calculate squared tolerance
    if(dtol .gt. 0.0_DP) then
      dtol2 = dtol**2
    else
      dtol2 = 0.0_DP
    end if

    ! Get the sub-vectors
    p_Ddef => Dwork(    1 :   n)  ! defect (/gradient) vector
    p_Ddir => Dwork(  n+1 : 2*n)  ! descend direction vector
    p_Dtmp => Dwork(2*n+1 : 3*n)  ! temporary vector

    ! First of all, copy Dx to Ddef
    call lalg_copyVector(Dx,p_Ddef,n)

    ! And clear the solution vector
    call lalg_clearVector(Dx,n)

    ! copy Ddef to Ddir
    call lalg_copyVector(p_Ddef,p_Ddir,n)

    ! Calculate initial gamma
    dgamma = lalg_scalarProduct(p_Ddef,p_Ddef,n)

    ! Check against tolerance if desired
    if(dgamma .le. dtol2) then

      ! No iterations where performed
      niter = 0

      ! Store final defect
      dtol = sqrt(dgamma)

      ! Success
      cinfo = QSOL_INFO_SUCCESS

      ! And return here
      return

    end if

    ! Okay, start the CG iteration
    do ite = 1, niter

      ! Calculate:
      ! 1. Dtmp := A * Ddir
      ! 2. dalpha := < Ddir, Dtmp >
      call qsol_mvmult_CG(n,Kld,Kcol,Da,p_Ddir,p_Dtmp,dalpha)

      ! Calculate alpha
      if(dalpha .eq. 0.0_DP) then

        ! Internal solver error
        cinfo = QSOL_INFO_INTERNAL_ERROR

        ! Store number of iterations
        niter = ite

        ! Store last defect
        if(dtol .gt. 0.0_DP) then
          dtol = sqrt(dgamma)
        end if

        return

      end if
      dalpha = dgamma / dalpha

      ! Calculate Dx = Dx + alpha*Ddir
      call lalg_vectorLinearComb(p_Ddir,Dx,dalpha,1.0_DP,n)

      ! Calculate Ddef = Ddef - alpha*Dtmp
      call lalg_vectorLinearComb(p_Dtmp,p_Ddef,-dalpha,1.0_DP,n)

      ! Calculate new gamma and beta
      dbeta = dgamma
      dgamma = lalg_scalarProduct(p_Ddef,p_Ddef,n)
      dbeta = dgamma / dbeta

      ! Check against tolerance if desired
      if(dgamma .le. dtol2) then

        ! No iterations where performed
        niter = ite

        ! Store final defect
        dtol = sqrt(dgamma)

        ! Success
        cinfo = QSOL_INFO_SUCCESS

        ! And return here
        return

      end if

      ! Calculate Ddir = Ddef + beta*Ddir
      call lalg_vectorLinearComb(p_Ddef,p_Ddir,1.0_DP,dbeta,n)

    end do

    ! If we come out here, then the maximum amount of iterations was performed.
    cinfo = QSOL_INFO_MAX_ITER

    ! Store final defect
    if(dtol .gt. 0.0_DP) then
      dtol = sqrt(dgamma)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine qsol_solveCG_SSOR_lsyssc(rmat, rvec, Dwork, cinfo, niter, dtol, &
                                      drelax)

!<description>
  ! This routine is a one-call CG solver.
!</description>

!<input>
  ! The system matrix.
  type(t_matrixScalar), intent(in) :: rmat

  ! OPTIONAL: The relaxation parameter of the SSOR preconditioner.
  ! If given, must be in range (0,2). If not given, 1 is used.
  real(DP), optional, intent(in) :: drelax
!</input>

!<output>
  ! One of the QSOL_INFO_XXXX constants describing the result.
  integer, intent(out) :: cinfo
!</output>

!<inputoutput>
  ! On entry, the vector containing the right-hand-side of the linear system.
  ! On exit, the vector containing the solution of the linear system.
  type(t_vectorScalar), intent(inout) :: rvec

  ! Work array. Its length must be at least 3*n.
  real(DP), dimension(:), target, intent(inout) :: Dwork

  ! On entry, the maximum number of allowed iterations. Must be > 0.
  ! On exit, the total number of performed iterations.
  integer, intent(inout) :: niter

  ! On entry, the absolute tolerance in euclidian norm that is to be reached.
  !   If set to a value <= 0, no residual checks are performed, and the solver
  !   always performs niter iterations, unless an error occurs.
  ! On exit,
  !   if dtol was  > 0 on entry, the absolute final defect,
  !   if dtol was <= 0 on entry, dtol is left unchanged.
  real(DP), intent(inout) :: dtol
!</inputoutput>

!</suboutine>

  ! local variables
  integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiag
  real(DP), dimension(:), pointer :: p_Da, p_Dx
  real(DP) :: drlx

    if(present(drelax)) then
      drlx = drelax
    else
      drlx = 1.0_DP
    end if

    ! Get the arrays
    call lsyssc_getbase_Kld(rmat, p_Kld)
    call lsyssc_getbase_Kcol(rmat, p_Kcol)
    call lsyssc_getbase_Kdiagonal(rmat, p_Kdiag)

    call lsyssc_getbase_double(rmat, p_Da)
    call lsyssc_getbase_double(rvec, p_Dx)

    ! Call double version
    call qsol_solveCG_SSOR_double(rmat%NEQ, p_Kld, p_Kcol, p_Kdiag, p_Da, &
                                  p_Dx, Dwork, cinfo, niter, dtol, drlx)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine qsol_solveCG_SSOR_double(n, Kld, Kcol, Kdiag, Da, Dx, Dwork, &
                                      cinfo, niter, dtol, drelax)

!<description>
  ! This routine is a one-call CG solver with built-in SSOR preconditioner.
!</description>

!<input>
  ! The size of the linear system.
  integer, intent(in) :: n

  ! The Kld array of the matrix.
  integer, dimension(:), intent(in) :: Kld

  ! The Kcol array of the matrix.
  integer, dimension(:), intent(in) :: Kcol

  ! The Kdiagonal array of the matrix.
  integer, dimension(:), intent(in) :: Kdiag

  ! The data array of the matrix.
  real(DP), dimension(:), intent(in) :: Da

  ! OPTIONAL: The relaxation parameter of the SSOR preconditioner.
  ! If given, must be in range (0,2). If not given, 1 is used.
  real(DP), optional, intent(in) :: drelax
!</input>

!<output>
  ! One of the QSOL_INFO_XXXX constants describing the result.
  integer, intent(out) :: cinfo
!</output>

!<inputoutput>
  ! On entry, the vector containing the right-hand-side of the linear system.
  ! On exit, the vector containing the solution of the linear system.
  real(DP), dimension(:), intent(inout) :: Dx

  ! Work array. Its length must be at least 3*n.
  real(DP), dimension(:), target, intent(inout) :: Dwork

  ! On entry, the maximum number of allowed iterations. Must be > 0.
  ! On exit, the total number of performed iterations.
  integer, intent(inout) :: niter

  ! On entry, the absolute tolerance in euclidian norm that is to be reached.
  !   If set to a value <= 0, no residual checks are performed, and the solver
  !   always performs niter iterations, unless an error occurs.
  ! On exit,
  !   if dtol was  > 0 on entry, the absolute final defect,
  !   if dtol was <= 0 on entry, dtol is left unchanged.
  real(DP), intent(inout) :: dtol
!</inputoutput>

!</suboutine>

  ! local CG variables
  real(DP) :: dalpha, dbeta, dgamma, ddef, drlx
  integer :: ite

  ! temporary sub-vectors
  real(DP), dimension(:), pointer :: p_Ddef, p_Ddir, p_Dtmp

    ! Choose relaxation parameter.
    if(present(drelax)) then
      drlx = drelax
    else
      drlx = 1.0_DP
    end if

    ! Get the sub-vectors
    p_Ddef => Dwork(    1 :   n)  ! defect (/gradient) vector
    p_Ddir => Dwork(  n+1 : 2*n)  ! descend direction vector
    p_Dtmp => Dwork(2*n+1 : 3*n)  ! temporary vector

    ! First of all, copy Dx to Ddef
    call lalg_copyVector(Dx,p_Ddef,n)

    ! And clear the solution vector
    call lalg_clearVector(Dx,n)

    ! copy Ddef to Ddir
    call lalg_copyVector(p_Ddef,p_Ddir,n)

    ! Check defect?
    if(dtol .gt. 0.0_DP) then

      ! Calculate defect then
      ddef = lalg_norm(p_Ddef, LINALG_NORMEUCLID, n=n)

      ! Check against tolerance if desired
      if(ddef .le. dtol) then

        ! No iterations where performed
        niter = 0

        ! Store final defect
        dtol = ddef

        ! Success
        cinfo = QSOL_INFO_SUCCESS

        ! And return here
        return

      end if

    end if

    ! Call the preconditioner
    call qsol_precSSOR(n,Kld,Kcol,Kdiag,Da,p_Ddir,drlx)

    ! Calculate initial gamma
    dgamma = lalg_scalarProduct(p_Ddef,p_Ddir,n)

    ! Okay, start the CG iteration
    do ite = 1, niter

      ! Make sure gamma is not zero
      if(dgamma .eq. 0.0_DP) then

        ! Internal solver error
        cinfo = QSOL_INFO_INTERNAL_ERROR

        ! Store number of iterations
        niter = ite

        ! Store last defect
        if(dtol .gt. 0.0_DP) then
          dtol = ddef
        end if

        return

      end if

      ! Calculate:
      ! 1. Dtmp := A * Ddir
      ! 2. dalpha := < Ddir, Dtmp >
      call qsol_mvmult_CG(n,Kld,Kcol,Da,p_Ddir,p_Dtmp,dalpha)

      ! Calculate alpha
      if(dalpha .eq. 0.0_DP) then

        ! Internal solver error
        cinfo = QSOL_INFO_INTERNAL_ERROR

        ! Store number of iterations
        niter = ite

        ! Store last defect
        if(dtol .gt. 0.0_DP) then
          dtol = ddef
        end if

        return

      end if
      dalpha = dgamma / dalpha

      ! Calculate Dx = Dx + alpha*Ddir
      call lalg_vectorLinearComb(p_Ddir,Dx,dalpha,1.0_DP,n)

      ! Calculate Ddef = Ddef - alpha*Dtmp
      call lalg_vectorLinearComb(p_Dtmp,p_Ddef,-dalpha,1.0_DP,n)

      ! Check defect?
      if(dtol .gt. 0.0_DP) then

        ! Calculate defect then
        ddef = lalg_norm(p_Ddef, LINALG_NORMEUCLID, n=n)

        ! Check against tolerance if desired
        if(ddef .le. dtol) then

          ! Store number of iterations
          niter = ite

          ! Store final defec
          dtol = ddef

          ! Success
          cinfo = QSOL_INFO_SUCCESS

          ! And return here
          return

        end if

      end if

      ! Copy Ddef to Dtmp
      call lalg_copyVector(p_Ddef,p_Dtmp,n)

      ! Apply preconditioner onto Dtmp
      call qsol_precSSOR(n,Kld,Kcol,Kdiag,Da,p_Dtmp,drlx)

      ! Calculate new gamma and beta
      dbeta = dgamma
      dgamma = lalg_scalarProduct(p_Ddef,p_Dtmp,n)
      dbeta = dgamma / dbeta

      ! Calculate Ddir = Dtmp + beta*Ddir
      call lalg_vectorLinearComb(p_Dtmp,p_Ddir,1.0_DP,dbeta,n)

    end do

    ! If we come out here, then the maximum amount of iterations was performed.
    cinfo = QSOL_INFO_MAX_ITER

    ! Store final defect
    if(dtol .gt. 0.0_DP) then
      dtol = ddef
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine qsol_mvmult_CG(n,Kld,Kcol,Da,Dx,Dy,dalpha)

!<description>
  ! PRIVATE AUXILIARY ROUTINE:
  ! This routine performs two tasks at once:
  ! 1. Calculate <tex>$  y := A * x                      $</tex>
  ! 2. Calculate <tex>$  alpha := < x, A*x > = < x, y >  $</tex>
  !
  ! This routine is used by the CG solvers in this module.
!</description>

!<input>
  ! The size of Dx
  integer, intent(in) :: n

  ! The Kld array of the matrix.
  integer, dimension(*), intent(in) :: Kld

  ! The Kcol array of the matrix.
  integer, dimension(*), intent(in) :: Kcol

  ! The DA array of the matrix.
  real(DP), dimension(*), intent(in) :: Da

  ! The input vector that is to be multiplied by the matrix.
  real(DP), dimension(*), intent(in) :: Dx
!</input>

!<output>
  ! The output vector that recieves A*x
  real(DP), dimension(*), intent(out) :: Dy

  ! The result of the scalar product < x, A*x >
  real(DP), intent(out) :: dalpha
!</output>

!</subroutine>

  ! local variables
  integer :: i,j
  real(DP) :: dt

    dalpha = 0.0_DP

    !$omp parallel do default(shared) private(j,dt) &
    !$omp reduction(+:dalpha) if (n > qsol_perfconfig%NEQMIN_OMP)
    do i = 1, n
      dt = 0.0_DP
      do j = Kld(i), Kld(i+1)-1
        dt = dt + Da(j)*Dx(Kcol(j))
      end do ! j
      dalpha = dalpha + Dx(i)*dt
      Dy(i) = dt
    end do ! i
    !$omp end parallel do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine qsol_solveSOR_lsyssc(rmat, rvec, Dwork, cinfo, niter, dtol, drelax)

!<description>
  ! This routine is a one-call SOR solver.
!</description>

!<input>
  ! The system matrix.
  type(t_matrixScalar), intent(in) :: rmat

  ! OPTIONAL: The relaxation parameter of the SOR solver.
  ! If given, must be in range (0,2). If not given, 1 is used.
  real(DP), optional, intent(in) :: drelax
!</input>

!<output>
  ! One of the QSOL_INFO_XXXX constants describing the result.
  integer, intent(out) :: cinfo
!</output>

!<inputoutput>
  ! On entry, the vector containing the right-hand-side of the linear system.
  ! On exit, the vector containing the solution of the linear system.
  type(t_vectorScalar), intent(inout) :: rvec

  ! Work array. Its length must be at least n.
  real(DP), dimension(:), intent(inout) :: Dwork

  ! On entry, the maximum number of allowed iterations. Must be > 0.
  ! On exit, the total number of performed iterations.
  integer, intent(inout) :: niter

  ! On entry, the absolute tolerance in euclidian norm that is to be reached.
  !   If set to a value <= 0, no residual checks are performed, and the solver
  !   always performs niter iterations, unless an error occurs.
  ! On exit,
  !   if dtol was  > 0 on entry, the absolute final defect,
  !   if dtol was <= 0 on entry, dtol is left unchanged.
  real(DP), intent(inout) :: dtol
!</inputoutput>

!</suboutine>

  ! local variables
  integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiag
  real(DP), dimension(:), pointer :: p_Da, p_Dx
  real(DP) :: drlx

    if(present(drelax)) then
      drlx = drelax
    else
      drlx = 1.0_DP
    end if

    ! Make sure the work array has sufficient length.
    if(size(Dwork,1) .lt. rmat%NEQ) then
      call output_line ('ERROR: Insufficient work array!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'qsol_solveSOR_lsyssc')
      call sys_halt()
    end if

    ! Get the arrays
    call lsyssc_getbase_Kld(rmat, p_Kld)
    call lsyssc_getbase_Kcol(rmat, p_Kcol)
    call lsyssc_getbase_Kdiagonal(rmat, p_Kdiag)

    call lsyssc_getbase_double(rmat, p_Da)
    call lsyssc_getbase_double(rvec, p_Dx)

    ! Call double version
    call qsol_solveSOR_double(rmat%NEQ, p_Kld, p_Kcol, p_Kdiag, p_Da, &
                              p_Dx, Dwork, cinfo, niter, dtol, drlx)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine qsol_solveSOR_double(n, Kld, Kcol, Kdiag, Da, Dx, Dwork, cinfo, &
                                  niter, dtol, drelax)

!<description>
  ! This routine is a one-call SOR solver.
!</description>

!<input>
  ! The size of the linear system.
  integer, intent(in) :: n

  ! The Kld array of the matrix.
  integer, dimension(*), intent(in) :: Kld

  ! The Kcol array of the matrix.
  integer, dimension(*), intent(in) :: Kcol

  ! The Kdiagonal array of the matrix.
  integer, dimension(*), intent(in) :: Kdiag

  ! The data array of the matrix.
  real(DP), dimension(*), intent(in) :: Da

  ! OPTIONAL: The relaxation parameter of the SOR solver.
  ! If given, must be in range (0,2). If not given, 1 is used.
  real(DP), optional, intent(in) :: drelax
!</input>

!<output>
  ! One of the QSOL_INFO_XXXX constants describing the result.
  integer, intent(out) :: cinfo
!</output>

!<inputoutput>
  ! On entry, the vector containing the right-hand-side of the linear system.
  ! On exit, the vector containing the solution of the linear system.
  real(DP), dimension(*), intent(inout) :: Dx

  ! Work array. Its length must be at least n.
  real(DP), dimension(*), intent(inout) :: Dwork

  ! On entry, the maximum number of allowed iterations. Must be > 0.
  ! On exit, the total number of performed iterations.
  integer, intent(inout) :: niter

  ! On entry, the absolute tolerance in euclidian norm that is to be reached.
  !   If set to a value <= 0, no residual checks are performed, and the solver
  !   always performs niter iterations, unless an error occurs.
  ! On exit,
  !   if dtol was  > 0 on entry, the absolute final defect,
  !   if dtol was <= 0 on entry, dtol is left unchanged.
  real(DP), intent(inout) :: dtol
!</inputoutput>

!</suboutine>

  ! local variables
  integer :: ite,i,j
  real(DP) :: drlx,daux,ddef,dtol2

    ! Choose relaxation parameter
    if(present(drelax)) then
      drlx = drelax
    else
      drlx = 1.0_DP
    end if

    ! Calculate squared tolerance
    if(dtol .gt. 0.0_DP) then
      dtol2 = dtol**2
    else
      dtol2 = 0.0_DP
    end if

    ! Remark:
    ! This SOR implementation does not check the real defect against the
    ! specified tolerance.
    ! The real defect vector in the k-th iteration would be:
    !
    !  d_k := b - A * x_k
    !
    ! Instead we check the norm of the following vector:
    !
    !  q_k := b - L * x_k - (D + U) * x_{k-1}
    !
    ! But as we assume that SOR converges monotonely, it should hold that:
    !
    !  || d_k || <= || q_k ||
    !
    ! So (hopefully) there should not be any problems here.

    ! Copy the rhs to work array and reset iteration vector
    do i = 1, n
      Dwork(i) = Dx(i)
      Dx(i) = 0.0_DP
    end do

    ! Launch the SOR solver.
    do ite = 1, niter

      ! Reset the defect
      ddef = 0.0_DP

      ! Loop over all matrix rows
      do i = 1, n

        ! Get b(i)
        daux = Dwork(i)

        ! aux = b(i) - A(i,.) * x(.)
        do j = Kld(i), Kld(i+1)-1
          daux = daux - Da(j)*Dx(Kcol(j))
        end do

        ! Calculate new x(i)
        Dx(i) = Dx(i) + drlx*daux / Da(Kdiag(i))

        ! Update defect
        ddef = ddef + daux*daux

      end do

      ! Check defect
      if(ddef .le. dtol2) then

        ! Okay, we are done
        cinfo = QSOL_INFO_SUCCESS

        ! Store 'final defect'
        dtol = sqrt(ddef)

        ! Store number of iterations
        niter = ite

        ! Get out
        return

      end if

    end do

    ! Maximum number of iterations reached
    cinfo = QSOL_INFO_MAX_ITER

    ! Store final defect if desired
    if(dtol .gt. 0.0_DP) dtol = sqrt(ddef)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine qsol_precSOR_lsyssc(rmat,rvec,drelax)

!<description>
  ! Applies the SOR preconditioner onto a defect vector.
!</description>

!<input>
  ! The system matrix.
  type(t_matrixScalar), intent(in) :: rmat

  ! OPTIONAL: The relaxation parameter of the SOR preconditioner.
  ! If given, must be in range (0,2). If not given, 1 is used.
  real(DP), optional, intent(in) :: drelax
!<input>

!<inputoutput>
  ! The vector that is to be preconditioned.
  type(t_vectorScalar), intent(inout) :: rvec
!</inputoutput>

!</subroutine>

  ! local variables
  integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiag
  real(DP), dimension(:), pointer :: p_Da, p_Dx
  real(DP) :: drlx

    if(present(drelax)) then
      drlx = drelax
    else
      drlx = 1.0_DP
    end if

    ! Get the arrays
    call lsyssc_getbase_Kld(rmat, p_Kld)
    call lsyssc_getbase_Kcol(rmat, p_Kcol)
    call lsyssc_getbase_Kdiagonal(rmat, p_Kdiag)

    call lsyssc_getbase_double(rmat, p_Da)
    call lsyssc_getbase_double(rvec, p_Dx)

    ! Call double version
    call qsol_precSOR_double(rmat%NEQ, p_Kld, p_Kcol, p_Kdiag, p_Da, &
                              p_Dx, drlx)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine qsol_precSOR_double(n,Kld,Kcol,Kdiag,Da,Dx,drelax)

!<description>
  ! Applies the SOR preconditioner onto a defect vector.
!</description>

!<input>
  ! The dimension of the matrix.
  integer, intent(in) :: n

  ! The Kld array of the matrix
  integer, dimension(*), intent(in) :: Kld

  ! The Kcol array of the matrix
  integer, dimension(*), intent(in) :: Kcol

  ! The Kdiagonal array of the matrix
  integer, dimension(*), intent(in) :: Kdiag

  ! The data array of the matrix
  real(DP), dimension(*), intent(in) :: Da

  ! OPTIONAL: The relaxation parameter of the SOR preconditioner.
  ! If given, must be in range (0,2). If not given, 1 is used.
  real(DP), optional, intent(in) :: drelax
!<input>

!<inputoutput>
  ! The vector that is to be preconditioned.
  real(DP), dimension(*), intent(inout) :: Dx
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,j,k
  real(DP) :: dt, drlx

    ! Choose relaxation parameter
    if(present(drelax)) then
      drlx = drelax
    else
      drlx = 1.0_DP
    end if

    if(drlx .ne. 1.0_DP) then

      ! Forward insertion
      do i = 1, n
        dt = 0.0_DP
        k = Kdiag(i)
        do j = Kld(i), k-1
          dt = dt + Da(j)*Dx(Kcol(j))
        end do
        Dx(i) = (Dx(i) - drlx*dt) / Da(k)
      end do

    else ! Unrelaxed GS preconditioner

      ! Forward insertion
      do i = 1, n
        dt = 0.0_DP
        k = Kdiag(i)
        do j = Kld(i), k-1
          dt = dt + Da(j)*Dx(Kcol(j))
        end do
        Dx(i) = (Dx(i) - dt) / Da(k)
      end do

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine qsol_precSSOR_lsyssc(rmat,rvec,drelax)

!<description>
  ! Applies the SSOR preconditioner onto a defect vector.
!</description>

!<input>
  ! The system matrix.
  type(t_matrixScalar), intent(in) :: rmat

  ! OPTIONAL: The relaxation parameter of the SSOR preconditioner.
  ! If given, must be in range (0,2). If not given, 1 is used.
  real(DP), optional, intent(in) :: drelax
!<input>

!<inputoutput>
  ! The vector that is to be preconditioned.
  type(t_vectorScalar), intent(inout) :: rvec
!</inputoutput>

!</subroutine>

  ! local variables
  integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiag
  real(DP), dimension(:), pointer :: p_Da, p_Dx
  real(DP) :: drlx

    if(present(drelax)) then
      drlx = drelax
    else
      drlx = 1.0_DP
    end if

    ! Get the arrays
    call lsyssc_getbase_Kld(rmat, p_Kld)
    call lsyssc_getbase_Kcol(rmat, p_Kcol)
    call lsyssc_getbase_Kdiagonal(rmat, p_Kdiag)

    call lsyssc_getbase_double(rmat, p_Da)
    call lsyssc_getbase_double(rvec, p_Dx)

    ! Call double version
    call qsol_precSSOR_double(rmat%NEQ, p_Kld, p_Kcol, p_Kdiag, p_Da, &
                              p_Dx, drlx)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine qsol_precSSOR_double(n,Kld,Kcol,Kdiag,Da,Dx,drelax)

!<description>
  ! Applies the SSOR preconditioner onto a defect vector.
!</description>

!<input>
  ! The dimension of the matrix.
  integer, intent(in) :: n

  ! The Kld array of the matrix
  integer, dimension(*), intent(in) :: Kld

  ! The Kcol array of the matrix
  integer, dimension(*), intent(in) :: Kcol

  ! The Kdiagonal array of the matrix
  integer, dimension(*), intent(in) :: Kdiag

  ! The data array of the matrix
  real(DP), dimension(*), intent(in) :: Da

  ! OPTIONAL: The relaxation parameter of the SSOR preconditioner.
  ! If given, must be in range (0,2). If not given, 1 is used.
  real(DP), optional, intent(in) :: drelax
!<input>

!<inputoutput>
  ! The vector that is to be preconditioned.
  real(DP), dimension(*), intent(inout) :: Dx
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,j,k
  real(DP) :: dt, drlx

    ! Choose relaxation parameter
    if(present(drelax)) then
      drlx = drelax
    else
      drlx = 1.0_DP
    end if

    if(drlx .ne. 1.0_DP) then

      ! Forward insertion
      do i = 1, n
        dt = 0.0_DP
        k = Kdiag(i)
        do j = Kld(i), k-1
          dt = dt + Da(j)*Dx(Kcol(j))
        end do
        Dx(i) = (Dx(i) - drlx*dt) / Da(k)
      end do

      ! Backward insertion
      do i = n, 1, -1
        dt = 0.0_DP
        k = Kdiag(i)
        do j = Kld(i+1)-1, k+1, -1
          dt = dt + Da(j)*Dx(Kcol(j))
        end do
        Dx(i) = Dx(i) - ((drlx*dt) / Da(k))
      end do

    else ! Unrelaxed SGS preconditioner

      ! Forward insertion
      do i = 1, n
        dt = 0.0_DP
        k = Kdiag(i)
        do j = Kld(i), k-1
          dt = dt + Da(j)*Dx(Kcol(j))
        end do
        Dx(i) = (Dx(i) - dt) / Da(k)
      end do

      ! Backward insertion
      do i = n, 1, -1
        dt = 0.0_DP
        k = Kdiag(i)
        do j = Kld(i+1)-1, k+1, -1
          dt = dt + Da(j)*Dx(Kcol(j))
        end do
        Dx(i) = Dx(i) - (dt / Da(k))
      end do

    end if

  end subroutine

  ! ************************************************************************

!<subroutine>

  pure subroutine qsol_solveDiagSchurComp (ndimA,ndimC,Du,Df,Da,Db,Dd,Dc,bsuccess)

!<description>
  ! This routine applies a Schur Complement decomposition to a 2x2 saddle point
  ! matrix where only entries on the diagonals of the diagonal blocks exist.
  !
  ! The system that is to be solved here is assumed to have the following shape:
  !
  ! <verb>
  !   (  A             +----------+  ) ( U ) = ( F )
  !   (       ..       |    B     |  ) ( U )   ( F )
  !   (             A  +----------+  ) ( U )   ( F )
  !   (  +----------+  C             ) ( U )   ( F )
  !   (  |    D     |       ..       ) ( U )   ( F )
  !   (  +----------+             C  ) ( U )   ( F )
  ! </verb>
  !
  ! or in short:
  !
  ! <verb>
  !   ( A B ) = (U) = (F)
  !   ( D C )   (U)   (F)
  ! </verb>
!</description>

!<input>
  ! Dimension of the A-matrix
  integer, intent(in) :: ndimA

  ! Dimension of the C-matrix
  integer, intent(in) :: ndimC

  ! Submatrix A, only diagonal entries.
  real(DP), dimension(*), intent(in) :: Da

  ! Entries in the submatrix B.
  real(DP), dimension(ndimA,*), intent(in) :: Db

  ! Entries in the submatrix D.
  real(DP), dimension(ndimC,*), intent(in) :: Dd

  ! Diagonal elements of the local system matrix C
  real(DP), dimension(*), intent(in) :: Dc

  ! Local RHS vector.
  real(DP), dimension(*), intent(in) :: Df
!</input>

!<output>
  ! SOlution vector.
  real(DP), dimension(*), intent(out) :: Du

  ! TRUE, if successful. FALSE if the system is indefinite.
  ! If FALSE, Du is undefined.
  logical, intent(out) :: bsuccess
!</output>

!</subroutine>

    ! local variables
    integer :: i,j,k,info
    integer, dimension(ndimC) :: Ipiv

    ! A^-1
    real(DP), dimension(ndimA) :: Dainv

    ! DA^-1, saved transposed for efficiency reasons
    real(DP), dimension(ndimA,ndimC) :: Ddainv

    ! S and S^-1
    real(DP), dimension(ndimC,ndimC) :: Ds, Dsinv

    ! Temporary RHS vector
    real(DP), dimension(ndimA) :: Dftemp1
    real(DP), dimension(ndimC) :: Dftemp2

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

    ! The system can be written as:
    !   ( A  B ) = (u) = (f)
    !   ( D  C )   (p)   (g)

    ! To solve with the Schur Complement, we have to solve in two steps.
    !
    ! 1.) Solve:  (C - D A^-1 B) p = g - D A^-1 f
    !         <=>                p = (C - D A^-1 B)^-1 (g - D A^-1 f)
    !         <=>                p = S^-1 (g - D A^-1 f)
    !
    ! 2.) Solve:               A u = f - B p
    !         <=>                u = A^-1 (f - B p)
    !
    ! * Solving A^-1 is easy as A is a diagonal matrix.
    ! * Solving S^-1 involves an inversion of a full matrix using LAPACK
    !   or a direct inversion.
    ! * D A^-1 can be easily computed by scaling D by A.
    !
    ! At first, compute A^-1 and copy f to Dftemp.
    do i=1,ndimA
      Dainv(i) = 1.0_DP/Da(i)
    end do

    do i=1,ndimA
      Dftemp1(i) = Df(i)
    end do

    do i=1,ndimC
      Dftemp2(i) = Df(ndimA+i)
    end do

    do i=1,ndimC
      do j=1,ndimA
        ! Compute D A^-1.
        ! D A^-1 is saved transposed for efficiency reasons!
        Ddainv(j,i) = Dd(i,j)*Dainv(j)

        ! Compute (g - D A^-1 f)
        Dftemp2(i) = Dftemp2(i) - Ddainv(j,i) * Df(j)
      end do
    end do

    ! Compute S = (C - D A^-1 B)
    ! (Order of loops changed for efficiency reasons...)

    do j=1,ndimC
      ! Compute (- (D A^-1) B)
      do i=1,ndimC
        Ds(i,j) = 0.0_DP
        do k=1,ndimA
          Ds(i,j) = Ds(i,j) - Ddainv(k,i) * Db(k,j)
        end do
      end do

      ! Add C; consists only of the diagonal.
      Ds(j,j) = Ds(j,j) + Dc(j)
    end do

    ! Invert S and compute p.
    !
    ! call mprim_invertMatrix(Ds, &
    !     Dftemp2,Du(ndimA+1:ndimA+ndimC),ndimC,2)

    select case(ndimC)
    case (1)
      if (Ds(1,1) .le. SYS_EPSREAL_DP) return
      Du(ndimA+1) = Dftemp2(1) / Ds(1,1)

    case (2)
      call mprim_invert2x2MatrixDirect(Ds,Dsinv,bsuccess)
      if (.not. bsuccess) return
      Du(ndimA+1) = Dsinv(1,1)*Dftemp2(1) + Dsinv(1,2)*Dftemp2(2)
      Du(ndimA+2) = Dsinv(2,1)*Dftemp2(1) + Dsinv(2,2)*Dftemp2(2)

    case (3)
      call mprim_invert3x3MatrixDirect(Ds,Dsinv,bsuccess)
      if (.not. bsuccess) return
      Du(ndimA+1) = Dsinv(1,1)*Dftemp2(1) + Dsinv(1,2)*Dftemp2(2) + Dsinv(1,3)*Dftemp2(3)
      Du(ndimA+2) = Dsinv(2,1)*Dftemp2(1) + Dsinv(2,2)*Dftemp2(2) + Dsinv(2,3)*Dftemp2(3)
      Du(ndimA+3) = Dsinv(3,1)*Dftemp2(1) + Dsinv(3,2)*Dftemp2(2) + Dsinv(3,3)*Dftemp2(3)

    case (4)
      call mprim_invert4x4MatrixDirect(Ds,Dsinv,bsuccess)
      if (.not. bsuccess) return
      Du(ndimA+1) = Dsinv(1,1)*Dftemp2(1) + Dsinv(1,2)*Dftemp2(2) &
                  + Dsinv(1,3)*Dftemp2(3) + Dsinv(1,4)*Dftemp2(4)
      Du(ndimA+2) = Dsinv(2,1)*Dftemp2(1) + Dsinv(2,2)*Dftemp2(2) &
                  + Dsinv(2,3)*Dftemp2(3) + Dsinv(2,4)*Dftemp2(4)
      Du(ndimA+3) = Dsinv(3,1)*Dftemp2(1) + Dsinv(3,2)*Dftemp2(2) &
                  + Dsinv(3,3)*Dftemp2(3) + Dsinv(3,4)*Dftemp2(4)
      Du(ndimA+4) = Dsinv(4,1)*Dftemp2(1) + Dsinv(4,2)*Dftemp2(2) &
                  + Dsinv(4,3)*Dftemp2(3) + Dsinv(4,4)*Dftemp2(4)

    case (5)
      call mprim_invert5x5MatrixDirect(Ds,Dsinv,bsuccess)
      if (.not. bsuccess) return
      Du(ndimA+1) = Dsinv(1,1)*Dftemp2(1) + Dsinv(1,2)*Dftemp2(2) &
                  + Dsinv(1,3)*Dftemp2(3) + Dsinv(1,4)*Dftemp2(4) &
                  + Dsinv(1,5)*Dftemp2(5)
      Du(ndimA+2) = Dsinv(2,1)*Dftemp2(1) + Dsinv(2,2)*Dftemp2(2) &
                  + Dsinv(2,3)*Dftemp2(3) + Dsinv(2,4)*Dftemp2(4) &
                  + Dsinv(2,5)*Dftemp2(5)
      Du(ndimA+3) = Dsinv(3,1)*Dftemp2(1) + Dsinv(3,2)*Dftemp2(2) &
                  + Dsinv(3,3)*Dftemp2(3) + Dsinv(3,4)*Dftemp2(4) &
                  + Dsinv(3,5)*Dftemp2(5)
      Du(ndimA+4) = Dsinv(4,1)*Dftemp2(1) + Dsinv(4,2)*Dftemp2(2) &
                  + Dsinv(4,3)*Dftemp2(3) + Dsinv(4,4)*Dftemp2(4) &
                  + Dsinv(4,5)*Dftemp2(5)
      Du(ndimA+5) = Dsinv(5,1)*Dftemp2(1) + Dsinv(5,2)*Dftemp2(2) &
                  + Dsinv(5,3)*Dftemp2(3) + Dsinv(5,4)*Dftemp2(4) &
                  + Dsinv(5,5)*Dftemp2(5)

    case (6)
      call mprim_invert6x6MatrixDirect(Ds,Dsinv,bsuccess)
      if (.not. bsuccess) return
      Du(ndimA+1) = Dsinv(1,1)*Dftemp2(1) + Dsinv(1,2)*Dftemp2(2) &
                  + Dsinv(1,3)*Dftemp2(3) + Dsinv(1,4)*Dftemp2(4) &
                  + Dsinv(1,5)*Dftemp2(5) + Dsinv(1,6)*Dftemp2(6)
      Du(ndimA+2) = Dsinv(2,1)*Dftemp2(1) + Dsinv(2,2)*Dftemp2(2) &
                  + Dsinv(2,3)*Dftemp2(3) + Dsinv(2,4)*Dftemp2(4) &
                  + Dsinv(2,5)*Dftemp2(5) + Dsinv(2,6)*Dftemp2(6)
      Du(ndimA+3) = Dsinv(3,1)*Dftemp2(1) + Dsinv(3,2)*Dftemp2(2) &
                  + Dsinv(3,3)*Dftemp2(3) + Dsinv(3,4)*Dftemp2(4) &
                  + Dsinv(3,5)*Dftemp2(5) + Dsinv(3,6)*Dftemp2(6)
      Du(ndimA+4) = Dsinv(4,1)*Dftemp2(1) + Dsinv(4,2)*Dftemp2(2) &
                  + Dsinv(4,3)*Dftemp2(3) + Dsinv(4,4)*Dftemp2(4) &
                  + Dsinv(4,5)*Dftemp2(5) + Dsinv(4,6)*Dftemp2(6)
      Du(ndimA+5) = Dsinv(5,1)*Dftemp2(1) + Dsinv(5,2)*Dftemp2(2) &
                  + Dsinv(5,3)*Dftemp2(3) + Dsinv(5,4)*Dftemp2(4) &
                  + Dsinv(5,5)*Dftemp2(5) + Dsinv(5,6)*Dftemp2(6)
      Du(ndimA+6) = Dsinv(6,1)*Dftemp2(1) + Dsinv(6,2)*Dftemp2(2) &
                  + Dsinv(6,3)*Dftemp2(3) + Dsinv(6,4)*Dftemp2(4) &
                  + Dsinv(6,5)*Dftemp2(5) + Dsinv(6,6)*Dftemp2(6)

    case default
      ! Use LAPACK routine for general NxN system, where N > 6
      Ipiv=0
      call DGESV(ndimC,1,Ds,ndimC,Ipiv,Dftemp2,ndimC,info)
      if (info .ne. 0) return
      Du(ndimA+1:ndimA+ndimC) = Dftemp2(1:ndimC)

    end select

    ! Calculate u = A^-1 (f - B p)
    do i=1,ndimA
      do j=1,ndimC
        Dftemp1(i) = Dftemp1(i) - Db(i,j) * Du(ndimA+j)
      end do
      Du(i) = Dainv(i) * Dftemp1(i)
    end do

    ! That's it.

    bsuccess = .true.

  end subroutine

  ! ************************************************************************

!<subroutine>

  pure subroutine qsol_solveTridiag (neq,iidxoffdiag,Dvec,Da,Db,Dd,Dbtemp)

!<description>
  ! This routine applies a modified Thomas algorithm to solve
  ! a tridiagonal system specified as a set of diagonals.
  ! The algorithm support not only 'simple' tridiagonal
  ! systems (with the offdiagonals directly above/below the diagonal)
  ! but also systems with the offdiagonals with some rows/columns
  ! distance from the diagonal.
  !
  ! The system that is to be solved here is assumed to have the following shape:
  !
  ! <verb>
  !   (  A      B         ) ( U ) = ( F )
  !   (     ..     ..     ) ( U )   ( F )
  !   (  D      A      B  ) ( U )   ( F )
  !   (     ..     ..     ) ( U )   ( F )
  !   (         D      A  ) ( U )   ( F )
  ! </verb>
  !
  ! The basic algorithm used here can be found e.g. at
  !   [ http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm ]
!</description>

!<input>
  ! Dimension of the matrix
  integer, intent(in) :: neq

  ! Index of the offdiagonal B/D matrix, i.e. the difference in the row/column
  ! number between the start of A and B/D. E.g.: An index of 1 indicates a
  ! trisiagonal matrix with B/D being the 1st offdiagonal. An index of 0 is
  ! not allowed as this would result in a diagonal matrix with vanishing B/D.
  integer, intent(in) :: iidxoffdiag

  ! Entries in the diagonal A.
  ! real(DP), dimension(neq), intent(inout) :: Da
  real(DP), dimension(:), intent(inout) :: Da

  ! Entries in the diagonal B.
  ! real(DP), dimension(neq-iidxoffdiag), intent(inout) :: Db
  real(DP), dimension(:), intent(inout) :: Db

  ! Entries in the diagonal D.
  ! real(DP), dimension(neq-iidxoffdiag), intent(inout) :: Dd
  real(DP), dimension(:), intent(inout) :: Dd
!</input>

!<inputoutput>
  ! Temporary vector of size 2*neq.
  ! real(DP), dimension(2*neq), intent(inout) :: Dbtemp
  real(DP), dimension(:), intent(inout) :: Dbtemp

  ! On Entry: RHS vector.
  ! On exit: Solution vector.
  ! real(DP), dimension(neq), intent(out) :: Dvec
  real(DP), dimension(:), intent(out) :: Dvec
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    real(DP) :: dx

    ! Forward sweep. Modify the coefficients to factorise the matrix.
    !
    ! We have e.g.
    !
    !  ( a1          b1             | f1 )    ( a1                  b1     | f1 )
    !  (     a2          b2         | f2 )    (     a2                  b2 | f2 )
    !  (         a3          b3     | f3 )    (         a3                 | f3 )
    !  ( d1          a4          b4 | f4 ) or (             a4             | f4 )
    !  (     d2          a5         | f5 )    (                 a5         | f5 )
    !  (         d3          a6     | f6 )    ( d1                  a6     | f6 )
    !  (             d4          a7 | f7 )    (     d2                  a7 | f7 )
    !
    ! a) Divide down the B-matrix.
    do i=1,min(iidxoffdiag,neq-iidxoffdiag)
      Dbtemp(i) = Db(i) / Da(i)
      Dbtemp(neq+i) = Dvec(i) / Da(i)
    end do

    ! Now we have in the first rows:
    !
    !  ( 1           c1             | g1 )    ( 1                   c1     | g1 )
    !  (     1           c2         | g2 )    (     1                   c2 | g2 )
    !  (         1           c3     | g3 )    (         a3                 | f3 )
    !  ( d1          a4          b4 | f4 ) or (             a4             | f4 )
    !  (     d2          a5         | f5 )    (                 a5         | f5 )
    !  (         d3          a6     | f6 )    ( d1                  a6     | f6 )
    !  (             d4          a7 | f7 )    (     d2                  a7 | f7 )

    ! There may be some entries left that are not touched by the offdiagonal.
    do i=min(iidxoffdiag,neq-iidxoffdiag)+1,iidxoffdiag
      Dbtemp(neq+i) = Dvec(i) / Da(i)
    end do

    ! leading to:
    !
    !  ( 1           c1             | g1 )    ( 1                   c1     | g1 )
    !  (     1           c2         | g2 )    (     1                   c2 | g2 )
    !  (         1           c3     | g3 )    (         1                  | g3 )
    !  ( d1          a4          b4 | f4 ) or (             1              | g4 )
    !  (     d2          a5         | f5 )    (                 1          | g5 )
    !  (         d3          a6     | f6 )    ( d1                  a6     | f6 )
    !  (             d4          a7 | f7 )    (     d2                  a7 | f7 )
    !
    ! b) Factorise the B matrix and eliminate D.
    do i = iidxoffdiag+1,neq-iidxoffdiag
      ! Calc redundant factor
      dx = Da(i) - Dd(i-iidxoffdiag) * Dbtemp(i-iidxoffdiag)

      ! Factorise B, result is written to Dbtemp(1..neq-iidxoffdiag)
      Dbtemp(i) = Db(i) / dx

      ! Modify the RHS appropriately. Result is written to
      ! Dbtemp(neq+1..2*neq)
      Dbtemp(neq+i) = &
          (Dvec(i) - Dd(i-iidxoffdiag) * Dbtemp(neq+i-iidxoffdiag)) / dx
    end do

    ! So we have:
    !
    !  ( 1           c1             | g1 )    ( 1                   c1     | g1 )
    !  (     1           c2         | g2 )    (     1                   c2 | g2 )
    !  (         1           c3     | g3 )    (         1                  | g3 )
    !  ( 0           1           c4 | g4 ) or (             1              | g4 )
    !  (     d2          a5         | f5 )    (                 1          | g5 )
    !  (         d3          a6     | f6 )    ( 0                   1      | g6 )
    !  (             d4          a7 | f7 )    (     0                   1  | g7 )
    !
    ! We have to continue the modification of the RHS for the remaining
    ! entries of D.
    do i = max(iidxoffdiag,neq-iidxoffdiag)+1,neq
      ! Calc redundant factor
      dx = Da(i) - Dd(i-iidxoffdiag) * Dbtemp(i-iidxoffdiag)

      ! Modify the RHS appropriately. Result is written to
      ! Dbtemp(neq+1..2*neq)
      Dbtemp(neq+i) = &
          (Dvec(i) - Dd(i-iidxoffdiag) * Dbtemp(neq+i-iidxoffdiag)) / dx
    end do

    ! So the matrix we have now is:
    !
    !  ( 1           c1             | g1 )    ( 1                   c1     | g1 )
    !  (     1           c2         | g2 )    (     1                   c2 | g2 )
    !  (         1           c3     | g3 )    (         1                  | g3 )
    !  ( 0           1           c4 | g4 ) or (             1              | g4 )
    !  (     0           1          | g5 )    (                 1          | g5 )
    !  (         0           1      | g6 )    ( 0                   1      | g6 )
    !  (             0           1  | g7 )    (     0                   1  | g7 )

    ! c) Back substitution to canculate the solution.
    ! At first take those elements of the factorised RHS which are
    ! already computed.
    do i=neq-iidxoffdiag+1,neq
      Dvec(i) = Dbtemp(neq+i)
    end do

    ! That means on the example:
    !
    !  ( 1           c1             | g1      )    ( 1                   c1     | g1      )
    !  (     1           c2         | g2      )    (     1                   c2 | g2      )
    !  (         1           c3     | g3      )    (         1                  | g3 = x3 )
    !  ( 0           1           c4 | g4      ) or (             1              | g4 = x4 )
    !  (     0           1          | g5 = x5 )    (                 1          | g5 = x5 )
    !  (         0           1      | g6 = x6 )    ( 0                   1      | g6 = x6 )
    !  (             0           1  | g7 = x7 )    (     0                   1  | g7 = x7 )

    ! Now do the back-substitution for the remaining entries.
    do i = neq-iidxoffdiag,1,-1
      Dvec(i) = Dbtemp(neq+i) - Dbtemp(i)*Dvec(i+iidxoffdiag)
    end do

  end subroutine

end module

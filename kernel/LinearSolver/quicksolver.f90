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
!# the linearsolver module, they are meant as 'low-level' solvers for module
!# which cannot use the linearsolver module for whatever reason, e.g. the
!# multilevelprojection module, which itself is included by the linearsolver
!# module.
!#
!# As the routines do not allocate memory, the caller is responsible to pass
!# a working array to the solver which is used to store temporary vectors.
!#
!# The following routines can be found in this module:
!#
!#  1.) qsol_solveCG
!#      -> Solves a linear system using an unpreconditioned CG solver.
!#         This solver needs 3*n working memory.
!#
!#  2.) qsol_solveCG_SSOR
!#      -> Solves a linear system using a CG solver with built-in SSOR
!#         preconditioner. This solver needs 3*n working memory.
!#
!#  3.) qsol_solveSOR
!#      -> Solves a linear system using a SOR solver.
!#         This solver needs 1*n working memory.
!#
!#  4.) qsol_precSOR
!#      -> Applies the SOR preconditioner onto a defect vector.
!#
!#  5.) qsol_precSSOR
!#      -> Applies the SSOR preconditioner onto a defect vector.
!#
!# </purpose>
!##############################################################################

module quicksolver

use fsystem
use linearalgebra

implicit none

private

public :: qsol_solveCG
public :: qsol_solveCG_SSOR
public :: qsol_solveSOR
public :: qsol_precSOR
public :: qsol_precSSOR

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

contains

  ! ***************************************************************************

!<subroutine>

  subroutine qsol_solveCG(n, Kld, Kcol, Da, Dx, Dwork, cinfo, niter, dtol)

!<description>
  ! This routine is a one-call CG solver.
!</description>

!<input>
  ! The size of the linear system.
  integer, intent(IN) :: n
  
  ! The Kld array of the matrix.
  integer, dimension(:), intent(IN) :: Kld
  
  ! The Kcol array of the matrix.
  integer, dimension(:), intent(IN) :: Kcol
  
  ! The data array of the matrix.
  real(DP), dimension(:), intent(IN) :: Da
!</input>

!<output>
  ! One of the QSOL_INFO_XXXX constants describing the result.
  integer, intent(INOUT) :: cinfo
!</output>

!<inputoutput>
  ! On entry, the vector containing the right-hand-side of the linear system.
  ! On exit, the vector containing the solution of the linear system.
  real(DP), dimension(:), intent(INOUT) :: Dx
  
  ! Work array. Its length must be at least 3*n.
  real(DP), dimension(:), target, intent(INOUT) :: Dwork
  
  ! On entry, the maximum number of allowed CG iterations. Must be > 0.
  ! On exit, the total number of performed CG iterations.
  integer, intent(INOUT) :: niter
  
  ! On entry, the absolute tolerance in euclidian norm that is to be reached.
  !   If set to a value <= 0, no residual checks are performed, and the solver
  !   always performs niter iterations, unless an error occurs.
  ! On exit,
  !   if dtol was  > 0 on entry, the absolute final defect,
  !   if dtol was <= 0 on entry, dtol is left unchanged.
  real(DP), intent(INOUT) :: dtol
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
    call lalg_copyVectorDble(Dx,p_Ddef,n)
    
    ! And clear the solution vector
    call lalg_clearVectorDble(Dx,n)
    
    ! copy Ddef to Ddir
    call lalg_copyVectorDble(p_Ddef,p_Ddir,n)
    
    ! Calculate initial gamma
    dgamma = lalg_scalarProductDble(p_Ddef,p_Ddef,n)

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
      call lalg_vectorLinearCombDble(p_Ddir,Dx,dalpha,1.0_DP,n)
      
      ! Calculate Ddef = Ddef - alpha*Dtmp
      call lalg_vectorLinearCombDble(p_Dtmp,p_Ddef,-dalpha,1.0_DP,n)
      
      ! Calculate new gamma and beta
      dbeta = dgamma
      dgamma = lalg_scalarProductDble(p_Ddef,p_Ddef,n)
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
      call lalg_vectorLinearCombDble(p_Ddef,p_Ddir,1.0_DP,dbeta,n)
    
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

  subroutine qsol_solveCG_SSOR(n, Kld, Kcol, Kdiag, Da, Dx, Dwork, cinfo, &
                               niter, dtol, drelax)

!<description>
  ! This routine is a one-call CG solver with built-in SSOR preconditioner.
!</description>

!<input>
  ! The size of the linear system.
  integer, intent(IN) :: n
  
  ! The Kld array of the matrix.
  integer, dimension(:), intent(IN) :: Kld
  
  ! The Kcol array of the matrix.
  integer, dimension(:), intent(IN) :: Kcol
  
  ! The Kdiagonal array of the matrix.
  integer, dimension(:), intent(IN) :: Kdiag
  
  ! The data array of the matrix.
  real(DP), dimension(:), intent(IN) :: Da

  ! The relaxation parameter of the SSOR preconditioner.
  ! If given, must be in range (0,2). If not given, 1 is used.
  real(DP), optional, intent(IN) :: drelax
!</input>

!<output>
  ! One of the QSOL_INFO_XXXX constants describing the result.
  integer, intent(INOUT) :: cinfo
!</output>

!<inputoutput>
  ! On entry, the vector containing the right-hand-side of the linear system.
  ! On exit, the vector containing the solution of the linear system.
  real(DP), dimension(:), intent(INOUT) :: Dx
  
  ! Work array. Its length must be at least 3*n.
  real(DP), dimension(:), target, intent(INOUT) :: Dwork
  
  ! On entry, the maximum number of allowed CG iterations. Must be > 0.
  ! On exit, the total number of performed CG iterations.
  integer, intent(INOUT) :: niter
  
  ! On entry, the absolute tolerance in euclidian norm that is to be reached.
  !   If set to a value <= 0, no residual checks are performed, and the solver
  !   always performs niter iterations, unless an error occurs.
  ! On exit,
  !   if dtol was  > 0 on entry, the absolute final defect,
  !   if dtol was <= 0 on entry, dtol is left unchanged.
  real(DP), intent(INOUT) :: dtol
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
    call lalg_copyVectorDble(Dx,p_Ddef,n)
    
    ! And clear the solution vector
    call lalg_clearVectorDble(Dx,n)
    
    ! copy Ddef to Ddir
    call lalg_copyVectorDble(p_Ddef,p_Ddir,n)
    
    ! Check defect?
    if(dtol .gt. 0.0_DP) then
    
      ! Calculate defect then
      ddef = lalg_normDble(p_Ddef, LINALG_NORMEUCLID, n=n)

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
    dgamma = lalg_scalarProductDble(p_Ddef,p_Ddir,n)
    
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
      call lalg_vectorLinearCombDble(p_Ddir,Dx,dalpha,1.0_DP,n)
      
      ! Calculate Ddef = Ddef - alpha*Dtmp
      call lalg_vectorLinearCombDble(p_Dtmp,p_Ddef,-dalpha,1.0_DP,n)
      
      ! Check defect?
      if(dtol .gt. 0.0_DP) then
      
        ! Calculate defect then
        ddef = lalg_normDble(p_Ddef, LINALG_NORMEUCLID, n=n)

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
      call lalg_copyVectorDble(p_Ddef,p_Dtmp,n)
      
      ! Apply preconditioner onto Dtmp
      call qsol_precSSOR(n,Kld,Kcol,Kdiag,Da,p_Dtmp,drlx)
      
      ! Calculate new gamma and beta
      dbeta = dgamma
      dgamma = lalg_scalarProductDble(p_Ddef,p_Dtmp,n)
      dbeta = dgamma / dbeta
      
      ! Calculate Ddir = Dtmp + beta*Ddir
      call lalg_vectorLinearCombDble(p_Dtmp,p_Ddir,1.0_DP,dbeta,n)
    
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
  
  pure subroutine qsol_mvmult_CG(n,Kld,Kcol,Da,Dx,Dy,dalpha)
  
!<description>
  ! PRIVATE AUXILIARY ROUTINE:
  ! This routine performs two tasks at once:
  ! 1. Calculate y := A * x
  ! 2. Calculate alpha := < x, A*x > = < x, y >
  !
  ! This routine is used by the CG solvers in this module.
!</description>

!<input>
  ! The size of Dx
  integer, intent(IN) :: n
  
  ! The Kld array of the mass matrix.
  integer, dimension(*), intent(IN) :: Kld
  
  ! The Kcol array of the mass matrix.
  integer, dimension(*), intent(IN) :: Kcol

  ! The DA array of the mass matrix.
  real(DP), dimension(*), intent(IN) :: Da
  
  ! The input vector that is to be multiplied by the matrix.
  real(DP), dimension(*), intent(IN) :: Dx
!</input>

!<output>
  ! The output vector that recieves A*x
  real(DP), dimension(*), intent(OUT) :: Dy
  
  ! The result of the scalar product < x, A*x >
  real(DP), intent(OUT) :: dalpha
!</output>

!</subroutine>

  ! local variables
  integer :: i,j
  real(DP) :: dt
  
    dalpha = 0.0_DP
    
    !$omp parallel do default(shared) private(j,dt) reduction(+:dalpha)
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

  subroutine qsol_solveSOR(n, Kld, Kcol, Kdiag, Da, Dx, Dwork, cinfo, &
                           niter, dtol, drelax)

!<description>
  ! This routine is a one-call SOR solver.
!</description>

!<input>
  ! The size of the linear system.
  integer, intent(IN) :: n
  
  ! The Kld array of the matrix.
  integer, dimension(*), intent(IN) :: Kld
  
  ! The Kcol array of the matrix.
  integer, dimension(*), intent(IN) :: Kcol
  
  ! The Kdiagonal array of the matrix.
  integer, dimension(*), intent(IN) :: Kdiag
  
  ! The data array of the matrix.
  real(DP), dimension(*), intent(IN) :: Da

  ! The relaxation parameter of the SOR solver.
  ! If given, must be in range (0,2). If not given, 1 is used.
  real(DP), optional, intent(IN) :: drelax
!</input>

!<output>
  ! One of the QSOL_INFO_XXXX constants describing the result.
  integer, intent(INOUT) :: cinfo
!</output>

!<inputoutput>
  ! On entry, the vector containing the right-hand-side of the linear system.
  ! On exit, the vector containing the solution of the linear system.
  real(DP), dimension(*), intent(INOUT) :: Dx
  
  ! Work array. Its length must be at least 1*n.
  real(DP), dimension(*), target, intent(INOUT) :: Dwork
  
  ! On entry, the maximum number of allowed CG iterations. Must be > 0.
  ! On exit, the total number of performed CG iterations.
  integer, intent(INOUT) :: niter
  
  ! On entry, the absolute tolerance in euclidian norm that is to be reached.
  !   If set to a value <= 0, no residual checks are performed, and the solver
  !   always performs niter iterations, unless an error occurs.
  ! On exit,
  !   if dtol was  > 0 on entry, the absolute final defect,
  !   if dtol was <= 0 on entry, dtol is left unchanged.
  real(DP), intent(INOUT) :: dtol
!</inputoutput>

!</suboutine>

  ! local variables
  integer :: ite,i,j,k
  real(DP) :: drlx,daux,daii,ddef,dtol2
  
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
    ! As we assume that SOR converges monotonely, is should hold that:
    !
    !  || d_k || <= || q_k ||
    !
    ! So (hopefully) there shouldn't be any problems here.

    ! Copy the rhs to work array and set iteration vector to D^-1 * b
    !$omp parallel do
    do i = 1, n
      Dwork(i) = Dx(i)
      Dx(i) = Dx(i) / Da(Kdiag(i))
    end do
    !$omp end parallel do
    
    ! Launch the SOR solver
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
        Dx(i) = Dx(i) + drlx * (daux / Da(Kdiag(i)))

        ! Update defect
        ddef = ddef + daux*daux
        
      end do
      
      ! Check defect
      if(ddef .le. dtol2) then
        
        ! Okay, we're done
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

  pure subroutine qsol_precSOR(n,Kld,Kcol,Kdiag,Da,Dx,drelax)

!<description>
  ! Applies the SOR preconditioner onto a defect vector.
!</description>

!<input>
  ! The dimension of the matrix.
  integer, intent(IN) :: n
  
  ! The Kld array of the matrix
  integer, dimension(*), intent(IN) :: Kld
  
  ! The Kcol array of the matrix
  integer, dimension(*), intent(IN) :: Kcol
  
  ! The Kdiagonal array of the matrix
  integer, dimension(*), intent(IN) :: Kdiag
  
  ! The data array of the matrix
  real(DP), dimension(*), intent(IN) :: Da
  
  ! OPTIONAL: The relaxation parameter of SOR. If not given, 1 is used.
  real(DP), optional, intent(IN) :: drelax
!<input>

!<inputoutput>
  ! The vector that is to be preconditioned.
  real(DP), dimension(*), intent(INOUT) :: Dx
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

  pure subroutine qsol_precSSOR(n,Kld,Kcol,Kdiag,Da,Dx,drelax)

!<description>
  ! Applies the SSOR preconditioner onto a defect vector.
!</description>

!<input>
  ! The dimension of the matrix.
  integer, intent(IN) :: n
  
  ! The Kld array of the matrix
  integer, dimension(*), intent(IN) :: Kld
  
  ! The Kcol array of the matrix
  integer, dimension(*), intent(IN) :: Kcol
  
  ! The Kdiagonal array of the matrix
  integer, dimension(*), intent(IN) :: Kdiag
  
  ! The data array of the matrix
  real(DP), dimension(*), intent(IN) :: Da
  
  ! OPTIONAL: The relaxation parameter of SSOR. If not given, 1 is used.
  real(DP), optional, intent(IN) :: drelax
!<input>

!<inputoutput>
  ! The vector that is to be preconditioned.
  real(DP), dimension(*), intent(INOUT) :: Dx
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
  
end module
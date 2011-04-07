!##############################################################################
!# ****************************************************************************
!# <name> solverlinear </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides basic routines for solving linear system of the form
!#
!#   $$Au=f$$
!#
!# where $A$ is a block matrix and $f$ denotes a block vector. In addition,
!# a wrapper-routine for scalar equations is provided. On the highest level
!# of the linear solver hierarchy, this problem is solved by means of a
!# geometric multigrid approach. If only one multigrid-level is available than
!# the coarse grid solver is directly employed.
!#
!# The following routines are available:
!#
!# 1.) linsol_solveMultigrid = linsol_solveMultigridScalar /
!#                             linsol_solveMultigridBlock
!#     -> Solve the linear system Au=f by the multigrid method
!#
!# 2.) linsol_solveTwogrid
!#     -> Perform one linear two-grid step
!#
!# 3.) linsol_solveSinglegrid
!#     -> Solve the linear system Au=f by means of a single-grid solver
!#
!# 4.) linsol_solveUMFPACK
!#     -> Umfpack solver
!#
!# 5.) linsol_solveJacobi
!#     -> Jacobi solver
!#
!# 6.) linsol_solveSSOR
!#     -> (S)SOR solver
!#
!# 7.) linsol_solveBicgstab
!#     -> BiCGSTAB method
!#
!# 8.) linsol_solveFgmres
!#     -> flexible GMRES(m) method
!#
!# 9.) linsol_precond
!#     -> Perform preconditioning $C^{-1}Au=C^{-1}f$
!#
!# 10.) linsol_precondJacobi
!#      -> Jacobi preconditioner
!#
!# 11.) linsol_precondSSOR
!#      -> (S)SOR preconditioner
!#
!# 12.) linsol_precondILU
!#      -> ILU preconditioner
!#
!# 13.) linsol_smooth
!#      -> Perform smoothing of $Au=f$
!#
!# 14.) linsol_smoothJacobi
!#      -> Jacobi smoother
!#
!# 15.) linsol_smoothSSOR
!#      -> (S)SOR smoother
!#
!# 16.) linsol_smoothILU
!#      -> ILU smoother
!# </purpose>
!##############################################################################

module solverlinear

  use boundaryfilter
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solveraux
  use sortstrategy
  use storage

  implicit none

  private
  public :: linsol_solveMultigrid

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface linsol_solveMultigrid
    module procedure linsol_solveMultigridScalar
    module procedure linsol_solveMultigridBlock
  end interface

contains

!<subroutine>

  subroutine linsol_solveMultigridScalar(rproblemLevel, rsolver, ru, rf)

!<description>
    ! This subroutine solves the linear system
    !   $$Au=f$$
    ! for a scalar matrix $A$ and a scalar right-hand side vector $f$
    ! by means of the geometric multigrid method.
    ! NOTE: This subroutine serves as a wrapper for scalar equations.
    ! The scalar quantities are converted to their 1-block counterparts
    ! and the block version of the geometric multigrid method is called.
!</description>

!<input>
    ! Multigrid structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! right-hand side vector
    type(t_vectorScalar), intent(in) :: rf
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! Solution vector
    type(t_vectorScalar), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: ruBlock,rfBlock

    ! Convert scalar vectors to block vectors
    call lsysbl_createVecFromScalar(ru, ruBlock)
    call lsysbl_createVecFromScalar(rf, rfBlock)

    ! Perform geometric multigrid for the 1-block vector u
    call linsol_solveMultigridBlock(rproblemLevel, rsolver, ruBlock, rfBlock)

    ! Release auxiliary block vectors
    call lsysbl_releaseVector(ruBlock)
    call lsysbl_releaseVector(rfBlock)
  end subroutine linsol_solveMultigridScalar

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_solveMultigridBlock(rproblemLevel, rsolver, ru, rf)

!<description>
    ! This subroutine solves the linear system
    !   $$Au=f$$
    ! by means of the geometric multigrid method. If only one multigrid-level
    ! is available the coarse grid solver is directly employed.
    ! Otherwise, the the linear two-grid solver is called recursively starting
    ! at the highest multigrid level.
!</description>

!<input>
    ! Multigrid structure
    type(t_problemLevel) :: rproblemLevel

    ! right-hand side vector
    type(t_vectorBlock), intent(in) :: rf
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! Solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solverMultigrid), pointer :: p_solverMultigrid
    type(t_solver), pointer :: p_solverSinglegrid
    type(t_solver) :: rsolverTemp
    type(t_vectorBlock), pointer :: p_rres
    real(DP) :: doldDefect
    integer :: imgstep,mgcycle

    ! Check compatibility
    call lsysbl_isVectorCompatible(ru, rf)

    ! Check for zero right-hand side
    if (rsolver%drhsZero > 0.0_DP) then
      rsolver%dinitialRHS = lsysbl_vectorNorm(rf, LINALG_NORMMAX)

      if (rsolver%dinitialRHS .le. rsolver%drhsZero) then
        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
          call output_line('!!! Zero initial right-hand side '//trim(&
                           sys_sdEL(rsolver%dinitialRHS,5))//' !!!')
        end if

        ! Clear solution vector and adjust solver status
        call lsysbl_clearVector(ru)
        rsolver%istatus          = SV_ZERO_RHS
        rsolver%dconvergenceRate = 0.0_DP

        ! That is it, return.
        return
      end if
    end if

    ! What kind of solver are we?
    select case(rsolver%csolverType)
    case (SV_LINEAR)
      !-------------------------------------------------------------------------
      ! Single grid solver: Au=b

      call linsol_solveSinglegrid(rsolver, ru, rf)

    case (SV_LINEARMG)

      ! Check if multigrid solver exists
      if (.not.associated(rsolver%p_solverMultigrid)) then
        call output_line('Multigrid solver does not exists!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linsol_solveMultigridBlock')
        call sys_halt()
      end if

      ! Set pointer
      p_solverMultigrid => rsolver%p_solverMultigrid

      ! Check if multigrid solver is called for only one grid level
      if (p_solverMultigrid%nlmin .eq. p_solverMultigrid%nlmax) then

        !-----------------------------------------------------------------------
        ! Single grid solver: Au=b

        ! Check if single-grid solver exists
        if (.not.associated(p_solverMultigrid%p_solverCoarsegrid)) then
          call output_line('Coarsegrid solver does not exists!',&
              OU_CLASS_ERROR,OU_MODE_STD,'linsol_solveMultigridBlock')
          call sys_halt()
        end if

        ! Set pointer
        p_solverSinglegrid => p_solverMultigrid%p_solverCoarsegrid

        ! The outer nonlinear solver may have modified the linear multigrid
        ! solver, e.g., Newton`s method may have defined an absolute tolerance
        ! for the linear solver to satisfy some sufficient decrease condition.
        ! Hence, copy the original input parameters of the linear single grid
        ! solver to the temporal structure rsolverTemp and impose the new values.
        call solver_copySolver(p_solverSinglegrid, rsolverTemp, .true., .false.)
        call solver_copySolver(rsolver, p_solverSinglegrid, .true., .false.)

        ! Perform single-grid solution
        call linsol_solveSinglegrid(p_solverSinglegrid, ru, rf)

        ! Compute statistical data from the single-grid solver
        call solver_statistics(p_solverSinglegrid, p_solverSinglegrid%iiterations)

        ! Restore the old data from the temporal structure rsolverTemp
        ! and adopt the output parameters from the single grid solver
        ! Restore the old data from the linear coarse grid solver
        call solver_copySolver(rsolverTemp, p_solverSinglegrid, .true., .false.)
        call solver_copySolver(p_solverSinglegrid, rsolver, .false., .true.)


        select case(p_solverSinglegrid%isolver)
        case(LINSOL_SOLVER_JACOBI,&
             LINSOL_SOLVER_SOR,&
             LINSOL_SOLVER_SSOR)
          ! Jacobi- or (S)SOR solver
          if (rsolver%ioutputLevel .ge. SV_IOLEVEL_INFO) then
            call output_lbrk()
            call output_separator(OU_SEP_TILDE)
            call output_line('Single-grid solution         '//solver_getstatus(rsolver))
            call output_line('Number of linear iterations: '//trim(sys_siL(rsolver%iiterations,5)))
            call output_separator(OU_SEP_TILDE)
            call output_lbrk()
          end if

        case (LINSOL_SOLVER_BICGSTAB,&
              LINSOL_SOLVER_GMRES,&
              LINSOL_SOLVER_UMFPACK4)
          ! BiCGSTAB or GMRES solver or direct solver
          if (rsolver%ioutputLevel .ge. SV_IOLEVEL_INFO) then
            call output_lbrk()
            call output_separator(OU_SEP_TILDE)
            call output_line('Single-grid solution         '//solver_getstatus(rsolver))
            call output_line('Number of linear iterations: '//trim(sys_siL(rsolver%iiterations,5)))
            call output_line('Convergence rate:            '//trim(sys_sdEL(rsolver%dconvergenceRate,5)))
            call output_separator(OU_SEP_TILDE)
            call output_lbrk()
          end if

        case default
          call output_line('Unsupported single-grid solver!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'linsol_solveMultigridBlock')
          call sys_halt()
        end select

      else

        !-----------------------------------------------------------------------
        ! Multigrid solver: Au=b

        ! Set pointer for residual vector
        p_rres => p_solverMultigrid%rtempVectors(3*rproblemLevel%ilev-&
                                                 2*p_solverMultigrid%nlmin)

        ! Compute the initial linear residual: res=f-A*u
        call lsysbl_copyVector(rf, p_rres)
        call lsysbl_blockMatVec(p_solverMultigrid%rmatrix(rproblemLevel%ilev),&
                                ru, p_rres, -1.0_DP, 1.0_DP)

        ! Compute norm of initial defect
        rsolver%dinitialDefect = lsysbl_vectorNorm(p_rres, rsolver%iresNorm)
        rsolver%dfinalDefect   = rsolver%dinitialDefect
        doldDefect             = rsolver%dinitialDefect

        ! Check if initial residual is too large ...
        if (solver_testDivergence(rsolver)) then
          if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
            call output_line('!!! Norm of initial residual is too large '//&
                             trim(sys_sdEL(rsolver%dinitialDefect,5))//' !!!')
          end if

          ! Clear solution vector and adjust solver status
          call lsysbl_clearVector(ru)
          rsolver%istatus = SV_INF_DEF

          ! That is it, return.
          return

        elseif (rsolver%dinitialDefect .le. rsolver%ddefZero) then
          ! ... or if it satisfies the desired tolerance already
          if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
            call output_line('!!! Zero initial residual '//&
                             trim(sys_sdEL(rsolver%dinitialDefect,5))//' !!!')
          end if

          ! Clear solution vector and adjust solver status
          call lsysbl_clearVector(ru)
          rsolver%istatus = SV_ZERO_DEF

          ! That is it, return.
          return
        end if


        ! Perform prescribed number of multigrid steps
        mgstep: do imgstep = 1, p_solverMultigrid%ilmax

          ! Perform one linear two-grid step
          mgcycle = merge(1, 2, p_solverMultigrid%icycle .eq. 1)
          call linsol_solveTwogrid(rproblemLevel, rsolver, ru, rf, mgcycle)

          ! Compute the new linear defect: res=f-A*u
          call lsysbl_copyVector(rf, p_rres)
          call lsysbl_blockMatVec(p_solverMultigrid%rmatrix(rproblemLevel%ilev),&
                                  ru, p_rres, -1.0_DP, 1.0_DP)

          ! Compute norm of new linear defect
          rsolver%dfinalDefect = lsysbl_vectorNorm(p_rres, rsolver%iresNorm)

          if (rsolver%ioutputLevel .ge. SV_IOLEVEL_VERBOSE) then
            call output_lbrk()
            call output_separator(OU_SEP_TILDE)
            call output_line('Linear multigrid step:     '//trim(sys_siL(imgstep,5)))
            call output_line('Norm of residual:          '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
            call output_line('Norm of previous residual: '//trim(sys_sdEL(doldDefect,5)))
            call output_line('Variation of residual:     '//trim(sys_sdEL(&
                             rsolver%dfinalDefect/max(SYS_EPSREAL_DP, doldDefect),5)))
            call output_line('Improvement of residual:   '//trim(sys_sdEL(&
                             rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
            call output_separator(OU_SEP_TILDE)
            call output_lbrk()
          end if

          ! Check if residual increased too much
          if (solver_testDivergence(rsolver)) then
            if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
              call output_line('!!! Residual increased by factor '//trim(sys_sdEL(&
                  rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
            end if

            ! Clear solution vector and adjust solver status
            call lsysbl_clearVector(ru)
            rsolver%istatus          = SV_INCR_DEF
            rsolver%dconvergenceRate = 1.0_DP

            ! That is it, return.
            return
          end if

          ! Check convergence criteria
          if (imgstep .ge. p_solverMultigrid%ilmin) then
            if (solver_testConvergence(rsolver)) then
              rsolver%istatus = SV_CONVERGED
              exit mgstep
            end if
            if (solver_testStagnation(rsolver, doldDefect)) then
              rsolver%istatus = SV_STAGNATED
              exit mgstep
            end if
          end if

          ! Save norm of old defect
          doldDefect = rsolver%dfinalDefect
        end do mgstep

        ! Multigrid convergence rates
        call solver_statistics(rsolver, imgstep)

        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_INFO) then
          call output_lbrk()
          call output_separator(OU_SEP_PERC)
          call output_line('Linear multigrid solution '//solver_getstatus(rsolver))
          call output_line('Number of multigrid steps: '//trim(sys_siL(rsolver%iiterations,5)))
          call output_line('Norm of final residual:    '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
          call output_line('Norm of initial residual:  '//trim(sys_sdEL(rsolver%dinitialDefect,5)))
          call output_line('Improvement of residual:   '//trim(sys_sdEL(&
                           rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
          call output_line('Convergence rate:          '//trim(sys_sdEL(rsolver%dconvergenceRate,5)))
          call output_separator(OU_SEP_PERC)
          call output_lbrk()
        end if

      end if


    case default
      call output_line(' Invalid solver!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linsol_solveMultigridBlock')
      call sys_halt()
    end select
  end subroutine linsol_solveMultigridBlock

  ! ***************************************************************************

!<subroutine>

  recursive subroutine linsol_solveTwogrid(rproblemLevel, rsolver, ruf, rff, mgcycle)

!<description>
    ! This subroutine performs one linear two-grid step of geometric multigrid.
    ! Given the initial solution $U^f$ and the right-hand side $F^f$ on the fine grid,
    ! a prescribed number of pre-smoothing steps is performed. Afterward, the
    ! linear residual $R^f=F^f-A^f*u^f$ is evaluated and restricted to the coarse
    ! grid, i.e. $R^c=I_f^c R^f$. The solution on the coarse grid is initialised
    ! by zero and the coarse grid correction is computed by calling the two-grid
    ! algorithm recursively. Finally, the fine grid solution is corrected as
    ! $U^f:=U^f+I_c^f U^c$ and a prescribed number of post-smoothing steps is performed.
!</description>

!<input>
    ! Multigrid level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! Right-hand side vector on fine grid
    type(t_vectorBlock), intent(in) :: rff
!</input>

!<inputoutput>
    ! Multigrid solver structure
    type(t_solver), intent(inout) :: rsolver

    ! Solution vector in fine grid
    type(t_vectorBlock), intent(inout) :: ruf

    ! Identifier for multigrid cycle
    integer, intent(inout) :: mgcycle
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_raux,p_rresf,p_rresc
    integer :: icycle,ioffset,nsmooth


    ! Check if current level is already the coarse grid level
    if (rproblemLevel%ilev .eq. rsolver%p_solverMultigrid%nlmin) then

      !-------------------------------------------------------------------------
      ! Solve coarse grid problem

      ! Set pointers
      p_raux => rsolver%p_solverMultigrid%rtempVectors(rproblemLevel%ilev)

      ! Solve on the coarsest grid
      call lsysbl_clearVector(p_raux)
      call linsol_solveSinglegrid(rsolver%p_solverMultigrid%p_solverCoarsegrid,&
                                  p_raux, rff)

      ! Apply increment: u:=u+omega*aux
      if (rsolver%domega < 0.0_DP) then
        call lsysbl_vectorLinearComb(p_raux, ruf, 1.0_DP, 1.0_DP)
      else
        call lsysbl_vectorLinearComb(p_raux, ruf, rsolver%domega, 1.0_DP)
      end if

      ! Adjust multigrid cycle
      if (rsolver%p_solverMultigrid%icycle .eq. 0) mgcycle = 1

    else

      !-------------------------------------------------------------------------
      ! Recursive application of two-grid algorithm

      ! Set pointers
      ioffset =  -2*rsolver%p_solverMultigrid%nlmin+3*(rproblemLevel%ilev-1)
      p_raux  => rsolver%p_solverMultigrid%rtempVectors(ioffset+1)
      p_rresc => rsolver%p_solverMultigrid%rtempVectors(ioffset+2)
      p_rresf => rsolver%p_solverMultigrid%rtempVectors(ioffset+3)

      call lsysbl_clearVector(p_raux)
      call lsysbl_clearVector(p_rresc)

      ! Apply presmoothing steps
      if (rsolver%p_solverMultigrid%npresmooth > 0) then
        nsmooth = rsolver%p_solverMultigrid%npresmooth *&
                  rsolver%p_solverMultigrid%nsmoothfactor**&
                  (rsolver%p_solverMultigrid%nlmax-rproblemLevel%ilev)
        call linsol_smooth(rproblemLevel,&
                           rsolver%p_solverMultigrid%p_smoother(rproblemLevel%ilev),&
                           ruf, rff, nsmooth)
      end if

      ! Compute the linear residual: r=f-A*u
      call lsysbl_copyVector(rff, p_rresf)
      call lsysbl_blockMatVec(rsolver%p_solverMultigrid%rmatrix(rproblemLevel%ilev),&
                              ruf, p_rresf, -1.0_DP, 1.0_DP)

      ! Restrict the residual: resc=R(f-A*u)
      call solver_restrictionBlock(rproblemLevel%rtriangulation,&
                                   rproblemLevel%p_rproblemLevelCoarse%rtriangulation,&
                                   p_rresf, p_rresc)

      ! Filter boundary conditions
      call bdrf_filterVectorByValue(rsolver%rboundaryCondition, p_rresc, 0.0_DP)

      ! Compute the coarse grid correction recursively, whereby
      ! raux is initialised by zeros (see above)
      do icycle = 1, merge(1, mgcycle, rproblemLevel%p_rproblemLevelCoarse%ilev .eq.&
                                       rsolver%p_solverMultigrid%nlmin)
        call linsol_solveTwogrid(rproblemLevel%p_rproblemLevelCoarse,&
                                 rsolver, p_raux, p_rresc, mgcycle)
      end do

      ! Set pointers
      ioffset = -2*rsolver%p_solverMultigrid%nlmin+3*(rproblemLevel%ilev-1)
      p_raux  => rsolver%p_solverMultigrid%rtempVectors(ioffset+1)
      p_rresc => rsolver%p_solverMultigrid%rtempVectors(ioffset+2)
      p_rresf => rsolver%p_solverMultigrid%rtempVectors(ioffset+3)

      ! Perform prolongation: resf=P(raux)
      call lsysbl_clearVector(p_rresf)
      call solver_prolongationBlock(rproblemLevel%p_rproblemLevelCoarse%rtriangulation,&
                                    rproblemLevel%rtriangulation, p_raux, p_rresf)

      ! Filter boundary conditions
      call bdrf_filterVectorByValue(rsolver%rboundaryCondition, p_rresf, 0.0_DP)

      ! Update the solution: u:=u+omega*raux
      if (rsolver%domega < 0.0_DP) then
        call lsysbl_vectorLinearComb(p_rresf, ruf, 1.0_DP, 1.0_DP)
      else
        call lsysbl_vectorLinearComb(p_rresf, ruf, rsolver%domega, 1.0_DP)
      end if

      ! Apply postsmoothing steps
      if (rsolver%p_solverMultigrid%npostsmooth > 0) then
        nsmooth = rsolver%p_solverMultigrid%npostsmooth *&
                  rsolver%p_solverMultigrid%nsmoothfactor**&
                  (rsolver%p_solverMultigrid%nlmax-rproblemLevel%ilev)
        call linsol_smooth(rproblemLevel,&
                           rsolver%p_solverMultigrid%p_smoother(rproblemLevel%ilev),&
                           ruf, rff, nsmooth)
      end if

      ! Adjust multigrid cycle
      if ((rproblemLevel%ilev .eq. rsolver%p_solverMultigrid%nlmax) .and.&
          (rsolver%p_solverMultigrid%icycle .eq. 0)) mgcycle=2
    end if
  end subroutine linsol_solveTwogrid

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_solveSinglegrid(rsolver, ru, rf)

!<description>
    ! This subroutine solves the linear system
    !   $$Au=f$$
    ! by means of a linear single-grid solver
!</description>

!<input>
    ! Right-hand side vector
    type(t_vectorBlock), intent(in) :: rf
!</input>

!<inputoutput>
    ! Linear solver structure
    type(t_solver), intent(inout) :: rsolver

    ! Solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! Check if solver is correct
    if (rsolver%csolverType .ne. SV_LINEAR) then
      call output_line('Invalid solver type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linsol_solveSinglegrid')
      call sys_halt()
    end if

    ! Waht type of solver are we?
    select case(rsolver%isolver)
    case(LINSOL_SOLVER_UMFPACK4)
      ! Solve directly: UMFPACK4 (only for scalar problems)
      call linsol_solveUMFPACK(rsolver, ru, rf)


    case(LINSOL_SOLVER_JACOBI)
      ! Solve iteratively: Jacobi solver
      call linsol_solveJacobi(rsolver, ru, rf)


    case(LINSOL_SOLVER_SOR,LINSOL_SOLVER_SSOR)
      ! Solve iteratively: Gauss-Seidel/SOR solver
      call linsol_solveSSOR(rsolver, ru, rf)


    case(LINSOL_SOLVER_BICGSTAB)
      ! Solve iteratively: BiCGSTAB
      call linsol_solveBicgstab(rsolver, ru, rf)


    case(LINSOL_SOLVER_GMRES)
      ! Solve iteratively: Flexible GMRES
      call linsol_solveFgmres(rsolver, ru, rf)


    case default
      call output_line('Invalid linear solver!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linsol_solveSinglegrid')
      call sys_halt()
    end select
  end subroutine linsol_solveSinglegrid

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_solveUMFPACK(rsolver, ru, rf)

!<description>
    ! This subroutine solves the linear system
    !   $$Au=f$$
    ! directly by means of UMFPACK
!</description>

!<input>
    ! r.h.s. vector
    type(t_vectorBlock), intent(in) :: rf
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solverUMFPACK), pointer :: p_solver
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rr
    real(DP), dimension(:), pointer :: p_Du,p_Df
    real(DP), dimension(90) :: Dinfo
    integer :: ksys = 1

    ! Check compatibility
    call lsysbl_isVectorCompatible(ru, rf)

    ! Check if vector comes from a scalar one
    if (ru%nblocks .ne. 1 .or. rf%nblocks .ne. 1) then
      call output_line('UMFPACK4 can only be used for scalar equations!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linsol_solveUMFPACK')
      call sys_halt()
    end if

    ! Check if vector is stored in interleave format
    if ((ru%RvectorBlock(1)%NVAR .ne. 1) .or. &
        (rf%RvectorBlock(1)%NVAR .ne. 1)) then
      call output_line('Scalar vectors must not contain multiple components!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linsol_solveUMFPACK')
      call sys_halt()
    end if

    if (rsolver%drhsZero > 0.0_DP) then
      rsolver%dinitialRHS = lsysbl_vectorNorm(rf, LINALG_NORMMAX)

      if (rsolver%dinitialRHS .le. rsolver%drhsZero) then
        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
          call output_line('!!! Zero initial right-hand side '//trim(&
                           sys_sdEL(rsolver%dinitialRHS,5))//' !!!')
        end if

        ! Clear solution vector and adjust solver status
        call lsysbl_clearVector(ru)
        rsolver%istatus          = SV_ZERO_RHS
        rsolver%dconvergenceRate = 0.0_DP

        ! That is it, return.
        return
      end if
    end if

    ! Set pointers
    p_solver  => rsolver%p_solverUMFPACK
    p_rmatrix => p_solver%rmatrix
    p_rr      => p_solver%rtempVector
    call lsyssc_getbase_double(ru%RvectorBlock(1), p_Du)
    call lsyssc_getbase_double(rf%RvectorBlock(1), p_Df)

    ! Solve the system directly
    call UMF4SOL(ksys, p_Du, p_Df, p_solver%inumeric, p_solver%Dcontrol, Dinfo)

    ! Check Dinfo(1) if there is an error
    select case(int(Dinfo(1)))
    case (0)
      ! ok.

    case (1)
      ! Singular matrix
      call output_line('Matrix is singular!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linsol_solveUMFPACK')
      call sys_halt()

    case (-1)
      ! Insufficient memory
      call output_line('Insufficient memory!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linsol_solveUMFPACK')
      call sys_halt()

    case default
      ! Unknown error
      call output_line('Internal error!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linsol_solveUMFPACK')
    end select

    ! Compute linear residual
    call lsysbl_copyVector(rf, p_rr)
    call lsysbl_blockMatVec(p_rmatrix, ru, p_rr, -1.0_DP, 1.0_DP)
    rsolver%dfinalDefect = lsysbl_vectorNorm(p_rr, rsolver%iresNorm)

    ! Compute solver statistics
    call solver_statistics(rsolver, 1)
  end subroutine linsol_solveUMFPACK

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_solveJacobi(rsolver, ru, rf)

!<description>
    ! This subroutine solves the linear system Au=f
    ! by means of Jacobi iterations
!</description>

!<input>
    ! Right-hand side vector
    type(t_vectorBlock), intent(in) :: rf
!</input>
!<inputoutput>
    ! Linear solver structure
    type(t_solver) :: rsolver

    ! Solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixBlock),  pointer :: p_rmatrix
    type(t_vectorBlock),  pointer :: p_rres
    real(DP) :: doldDefect
    integer :: iiterations


    ! Check for zero right-hand side
    if (rsolver%drhsZero > 0.0_DP) then
      rsolver%dinitialRHS = lsysbl_vectorNorm(rf, LINALG_NORMMAX)

      if (rsolver%dinitialRHS .le. rsolver%drhsZero) then
        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
          call output_line('!!! Zero initial right-hand side '//trim(&
                           sys_sdEL(rsolver%dinitialRHS,5))//' !!!')
        end if

        ! Clear solution vector and adjust solver status
        call lsysbl_clearVector(ru)
        rsolver%istatus          = SV_ZERO_RHS
        rsolver%dconvergenceRate = 0.0_DP

        ! That is it, return.
        return
      end if
    end if

    ! Set pointers
    p_rmatrix => rsolver%p_solverJacobi%rmatrix
    p_rres    => rsolver%p_solverJacobi%rtempVector

    ! Check compatibility
    call lsysbl_isVectorCompatible(ru, rf)
    call lsysbl_isVectorCompatible(ru, p_rres)
    call lsysbl_isMatrixCompatible(ru, p_rmatrix, .false.)

    ! Compute initial residual
    call lsysbl_copyVector(rf, p_rres)
    call lsysbl_blockMatVec(p_rmatrix, ru, p_rres, -1.0_DP, 1.0_DP)

    ! Compute norm of initial defect
    rsolver%dinitialDefect   = lsysbl_vectorNorm(p_rres, rsolver%iresNorm)
    rsolver%dfinalDefect     = rsolver%dinitialDefect
    doldDefect               = rsolver%dinitialDefect

    ! Check if initial residual is too large ...
    if (solver_testDivergence(rsolver)) then
      if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
        call output_line('!!! Norm of initial residual is too large '//&
                         trim(sys_sdEL(rsolver%dinitialDefect,5))//' !!!')
      end if

      ! Clear solution vector and adjust solver status
      call lsysbl_clearVector(ru)
      rsolver%istatus = SV_INF_DEF

      ! That is it, return.
      return

    elseif (rsolver%dinitialDefect .le. rsolver%ddefZero) then
      ! ... or if it satisfies the desired tolerance already
      if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
        call output_line('!!! Zero initial residual '//&
                         trim(sys_sdEL(rsolver%dinitialDefect,5))//' !!!')
      end if

      ! Clear solution vector and adjust solver status
      call lsysbl_clearVector(ru)
      rsolver%istatus = SV_ZERO_DEF

      ! That is it, return.
      return
    end if


    ! Iterative correction
    correction: do iiterations = 1, rsolver%nmaxIterations

      ! Precondition the linear residual
      call linsol_precondJacobi(rsolver, p_rres)

      ! Update solution
      call lsysbl_vectorLinearComb(p_rres, ru, 1.0_DP, 1.0_DP)

      ! Compute residual
      call lsysbl_copyVector(rf, p_rres)
      call lsysbl_blockMatVec(p_rmatrix, ru, p_rres, -1.0_DP, 1.0_DP)

      ! Compute norm of residual
      rsolver%dfinalDefect = lsysbl_vectorNorm(p_rres, rsolver%iresNorm)

      if (rsolver%ioutputLevel .ge. SV_IOLEVEL_VERBOSE) then
        call output_lbrk()
        call output_separator(OU_SEP_TILDE)
        call output_line('Jacobi step:               '//trim(sys_siL(iiterations,5)))
        call output_line('Norm of residual:          '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
        call output_line('Norm of previous residual: '//trim(sys_sdEL(doldDefect,5)))
        call output_line('Variation of residual:     '//trim(sys_sdEL(&
                         rsolver%dfinalDefect/max(SYS_EPSREAL_DP, doldDefect),5)))
        call output_line('Improvement of residual:   '//trim(sys_sdEL(&
                         rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
        call output_separator(OU_SEP_TILDE)
        call output_lbrk()
      end if

      ! Check if residual increased too much
      if (solver_testDivergence(rsolver)) then
        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
          call output_line('!!! Residual increased by factor '//trim(sys_sdEL(&
                           rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
        end if

        ! Clear solution vector and adjust solver status
        call lsysbl_clearVector(ru)
        rsolver%istatus          = SV_INCR_DEF
        rsolver%dconvergenceRate = 1.0_DP

        ! That is it, return.
        return
      end if

      ! Check convergence criteria
      if (iiterations .ge. rsolver%nminIterations) then
        if (solver_testConvergence(rsolver)) then
          rsolver%istatus = SV_CONVERGED
          exit correction
        end if
        if (solver_testStagnation(rsolver, doldDefect)) then
          rsolver%istatus = SV_STAGNATED
          exit correction
        end if
      end if

      ! Save norm of old defect
      doldDefect = rsolver%dfinalDefect
    end do correction

    ! Compute convergence rate
    call solver_statistics(rsolver, iiterations)

    if (rsolver%ioutputLevel .ge. SV_IOLEVEL_INFO) then
      call output_lbrk()
      call output_separator(OU_SEP_PERC)
      call output_line('Jacobi solution            '//solver_getstatus(rsolver))
      call output_line('Number of Jacobi steps:    '//trim(sys_siL(rsolver%iiterations,5)))
      call output_line('Norm of final residual:    '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
      call output_line('Norm of initial residual:  '//trim(sys_sdEL(rsolver%dinitialDefect,5)))
      call output_line('Improvement of residual:   '//trim(sys_sdEL(&
                       rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
      call output_line('Convergence rate:          '//trim(sys_sdEL(rsolver%dconvergenceRate,5)))
      call output_separator(OU_SEP_PERC)
      call output_lbrk()
    end if

  end subroutine linsol_solveJacobi

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_solveSSOR(rsolver, ru, rf)

!<description>
    ! This subroutine solves the linear system
    !   $$Au=f$$
    ! by means of (S)SOR iterations
!</description>

!<input>
    ! Right-hand side vector
    type(t_vectorBlock), intent(in) :: rf
!</input>
!<inputoutput>
    ! Linear solver structure
    type(t_solver) :: rsolver

    ! Solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock),  pointer :: p_rres
    real(DP) :: doldDefect
    integer :: iiterations


    ! Check for zero right-hand side
    if (rsolver%drhsZero > 0.0_DP) then
      rsolver%dinitialRHS = lsysbl_vectorNorm(rf, LINALG_NORMMAX)

      if (rsolver%dinitialRHS .le. rsolver%drhsZero) then
        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
          call output_line('!!! Zero initial right-hand side '//trim(&
                           sys_sdEL(rsolver%dinitialRHS,5))//' !!!')
        end if

        ! Clear solution vector and adjust solver status
        call lsysbl_clearVector(ru)
        rsolver%istatus          = SV_ZERO_RHS
        rsolver%dconvergenceRate = 0.0_DP

        ! That is it, return.
        return
      end if
    end if

    ! Set pointers
    p_rmatrix => rsolver%p_solverSSOR%rmatrix
    p_rres    => rsolver%p_solverSSOR%rtempVector

    ! Check compatibility
    call lsysbl_isVectorCompatible(ru, rf)
    call lsysbl_isVectorCompatible(ru, p_rres)
    call lsysbl_isMatrixCompatible(ru, p_rmatrix, .false.)

    ! Compute initial residual
    call lsysbl_copyVector(rf, p_rres)
    call lsysbl_blockMatVec(p_rmatrix, ru, p_rres, -1.0_DP, 1.0_DP)

    ! Compute norm of initial defect
    rsolver%dinitialDefect = lsysbl_vectorNorm(p_rres, rsolver%iresNorm)
    rsolver%dfinalDefect   = rsolver%dinitialDefect
    doldDefect             = rsolver%dinitialDefect

    ! Check if initial residual is too large ...
    if (solver_testDivergence(rsolver)) then
      if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
        call output_line('!!! Norm of initial residual is too large '//&
                         trim(sys_sdEL(rsolver%dinitialDefect,5))//' !!!')
      end if

      ! Clear solution vector and adjust solver status
      call lsysbl_clearVector(ru)
      rsolver%istatus = SV_INF_DEF

      ! That is it, return.
      return

    elseif (rsolver%dinitialDefect .le. rsolver%ddefZero) then
      ! ... or if it satisfies the desired tolerance already
      if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
        call output_line('!!! Zero initial residual '//&
                         trim(sys_sdEL(rsolver%dinitialDefect,5))//' !!!')
      end if

      ! Clear solution vector and adjust solver status
      call lsysbl_clearVector(ru)
      rsolver%istatus = SV_ZERO_DEF

      ! That is it, return.
      return
    end if


    ! Iterative correction
    correction: do iiterations = 1, rsolver%nmaxIterations

      ! Precondition the linear residual
      call linsol_precondSSOR(rsolver, p_rres)

      ! Update solution
      call lsysbl_vectorLinearComb(p_rres, ru, 1.0_DP, 1.0_DP)

      ! Compute residual
      call lsysbl_copyVector(rf, p_rres)
      call lsysbl_blockMatVec(p_rmatrix, ru, p_rres, -1.0_DP, 1.0_DP)

      ! Compute norm of residual
      rsolver%dfinalDefect = lsysbl_vectorNorm(p_rres, rsolver%iresNorm)

      if (rsolver%ioutputLevel .ge. SV_IOLEVEL_VERBOSE) then
        call output_lbrk()
        call output_separator(OU_SEP_TILDE)
        call output_line('(S)SOR step:               '//trim(sys_siL(iiterations,5)))
        call output_line('Norm of residual:          '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
        call output_line('Norm of previous residual: '//trim(sys_sdEL(doldDefect,5)))
        call output_line('Variation of residual:     '//trim(sys_sdEL(&
                         rsolver%dfinalDefect/max(SYS_EPSREAL_DP, doldDefect),5)))
        call output_line('Improvement of residual:   '//trim(sys_sdEL(&
                         rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
        call output_separator(OU_SEP_TILDE)
        call output_lbrk()
      end if

      ! Check if residual increased too much
      if (solver_testDivergence(rsolver)) then
        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
          call output_line('!!! Residual increased by factor '//trim(sys_sdEL(&
                           rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
        end if

        ! Clear solution vector and adjust solver status
        call lsysbl_clearVector(ru)
        rsolver%istatus          = SV_INCR_DEF
        rsolver%dconvergenceRate = 1.0_DP

        ! That is it, return.
        return
      end if

      ! Check convergence criteria
      if (iiterations .ge. rsolver%nminIterations) then
        if (solver_testConvergence(rsolver)) then
          rsolver%istatus = SV_CONVERGED
          exit correction
        end if
        if (solver_testStagnation(rsolver, doldDefect)) then
          rsolver%istatus = SV_STAGNATED
          exit correction
        end if
      end if

      ! Save norm of old defect
      doldDefect = rsolver%dfinalDefect
    end do correction

    ! Compute convergence rate
    call solver_statistics(rsolver, iiterations)

    if (rsolver%ioutputLevel .ge. SV_IOLEVEL_INFO) then
      call output_lbrk()
      call output_separator(OU_SEP_PERC)
      call output_line('(S)SOR solution            '//solver_getstatus(rsolver))
      call output_line('Number of (S)SOR steps:    '//trim(sys_siL(rsolver%iiterations,5)))
      call output_line('Norm of final residual:    '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
      call output_line('Norm of initial residual:  '//trim(sys_sdEL(rsolver%dinitialDefect,5)))
      call output_line('Improvement of residual:   '//trim(sys_sdEL(&
                       rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
      call output_line('Convergence rate:          '//trim(sys_sdEL(rsolver%dconvergenceRate,5)))
      call output_separator(OU_SEP_PERC)
      call output_lbrk()
    end if

  end subroutine linsol_solveSSOR

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_solveBicgstab(rsolver, ru, rf)

!<description>
    ! This subroutine solves the linear system
    !   $$Au=f$ $
    ! by means of the BiCGSTAB method
!</description>

!<input>
    ! Right-hand side vector
    type(t_vectorBlock), intent(in) :: rf
!</input>

!<inputoutput>
    ! Linear solver structure
    type(t_solver) :: rsolver

    ! Solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! restart parameter: if the denominator in the computation of
    ! omega is too small, then the BiCGStab algorithm is restarted
    ! whereby the computed solution vector u_j is used as u_0
    real(DP), parameter :: EPS_RESTART = 1E-5_DP

    ! local variables
    type(t_solverBiCGSTAB), pointer :: p_solverBiCGSTAB
    type(t_solver), pointer :: p_solverPrecond
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rr,p_rr0,p_rr1,p_rp,p_rpa,p_rpa1,p_rsa,p_rsa1
    real(DP) :: alpha,beta,omega,rho0,rho1,nrm0,doldDefect
    integer :: neq,iiterations,iiterations0
    logical :: bprec


    ! Check for zero right-hand side
    if (rsolver%drhsZero > 0.0_DP) then
      rsolver%dinitialRHS = lsysbl_vectorNorm(rf, LINALG_NORMMAX)

      if (rsolver%dinitialRHS .le. rsolver%drhsZero) then
        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
          call output_line('!!! Zero initial right-hand side '//trim(&
                           sys_sdEL(rsolver%dinitialRHS,5))//' !!!')
        end if

        ! Clear solution vector and adjust solver status
        call lsysbl_clearVector(ru)
        rsolver%istatus          = SV_ZERO_RHS
        rsolver%dconvergenceRate = 0.0_DP

        ! That is it, return.
        return
      end if
    end if


    ! Set pointers
    p_solverBiCGSTAB => rsolver%p_solverBiCGSTAB
    p_rmatrix        => p_solverBiCGSTAB%rmatrix
    p_rr             => p_solverBiCGSTAB%rtempVectors(1)
    p_rr0            => p_solverBiCGSTAB%rtempVectors(2)
    p_rr1            => p_solverBiCGSTAB%rtempVectors(3)
    p_rp             => p_solverBiCGSTAB%rtempVectors(4)
    p_rpa            => p_solverBiCGSTAB%rtempVectors(5)
    p_rpa1           => p_solverBiCGSTAB%rtempVectors(6)
    p_rsa            => p_solverBiCGSTAB%rtempVectors(7)
    p_rsa1           => p_solverBiCGSTAB%rtempVectors(8)

    ! Check for preconditioner
    if (associated(p_solverBiCGSTAB%p_precond)) then
      bprec = .true.
      p_solverPrecond => p_solverBiCGSTAB%p_precond
    else
      bprec = .false.
    end if

    ! Check compatibility
    call lsysbl_isVectorCompatible(ru, rf)
    call lsysbl_isMatrixCompatible(ru, p_rmatrix, .false.)

    ! Initialisation
    neq         = ru%NEQ
    iiterations = 1

    ! Compute initial residual r=C(-1)*(f-Au)
    call lsysbl_copyVector(rf, p_rr)
    call lsysbl_blockMatVec(p_rmatrix, ru, p_rr, -1.0_DP, 1.0_DP)

    ! Compute norm of initial defect
    rsolver%dinitialDefect   = lsysbl_vectorNorm(p_rr, rsolver%iresNorm)
    rsolver%dfinalDefect     = rsolver%dinitialDefect
    doldDefect               = rsolver%dinitialDefect

    ! Initialisation
100 rho0  = 1.0_DP
    alpha = 0.0_DP
    omega = 1.0_DP
    iiterations0  = iiterations

    ! ... and apply preconditioner if required
    if (bprec) then
      call lsysbl_copyVector(p_rr, p_rr1)
      call linsol_precond(p_solverPrecond, p_rr)
    end if

    ! Check if initial residual is too large ...
    if (solver_testDivergence(rsolver)) then
      if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
        call output_line('!!! Norm of initial residual is too large '//&
                         trim(sys_sdEL(rsolver%dinitialDefect,5))//' !!!')
      end if

      ! Clear solution vector and adjust solver status
      call lsysbl_clearVector(ru)
      rsolver%istatus = SV_INF_DEF

      ! That is it, return.
      return

    elseif (rsolver%dinitialDefect .le. rsolver%ddefZero) then
      ! ... or if it satisfies the desired tolerance already
      if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
        call output_line('!!! Zero initial residual '//&
                         trim(sys_sdEL(rsolver%dinitialDefect,5))//' !!!')
      end if

      ! Clear solution vector and adjust solver status
      call lsysbl_clearVector(ru)
      rsolver%istatus = SV_ZERO_DEF

      ! That is it, return.
      return
    end if


    ! Store initial residual
    call lsysbl_copyVector(p_rr, p_rr0)
    nrm0 = lsysbl_vectorNorm(p_rr, rsolver%iresNorm)

    ! Iterative correction
    correction: do iiterations = iiterations0, rsolver%nmaxIterations

      rho1 = lsysbl_scalarProduct(p_rr, p_rr0)

      if (iiterations .eq. 1) then
        ! Set p:=r
        call lsysbl_copyVector(p_rr, p_rp)
      else
        beta = (rho1/rho0)*(alpha/omega)
        rho0 = rho1

        ! Compute p:=r+beta*p-beta*omega*pa
        call lsysbl_vectorLinearComb(p_rr, p_rp, 1.0_DP, beta)
        call lsysbl_vectorLinearComb(p_rpa, p_rp, -beta*omega, 1.0_DP)
      end if

      ! Compute pa:=C(-1)*A*p
      call lsysbl_blockMatVec(p_rmatrix, p_rp, p_rpa, 1.0_DP, 0.0_DP)
      if (bprec) then
        call lsysbl_copyVector(p_rpa, p_rpa1)
        call linsol_precond(p_solverPrecond, p_rpa)
      end if

      ! Compute alpha:=(r,r0)/(pa,r0)
      alpha = lsysbl_scalarProduct(p_rpa, p_rr0)

      ! Check if algorithm needs to be restarted
      if ((iiterations > iiterations0) .and. &
          (abs(alpha) .le. lsysbl_vectorNorm(p_rpa,&
                           rsolver%iresNorm)*nrm0*EPS_RESTART)) then
        call lsysbl_copyVector(rf, p_rr)
        call lsysbl_blockMatVec(p_rmatrix, ru, p_rr, -1.0_DP, 1.0_DP)
        rsolver%dfinalDefect = lsysbl_vectorNorm(p_rr, rsolver%iresNorm)
        goto 100
      end if

      alpha = rho1/alpha

      ! Compute r:=r-alpha*pa
      call lsysbl_vectorLinearComb(p_rpa, p_rr, -alpha, 1.0_DP)
      if (bprec) call lsysbl_vectorLinearComb(p_rpa1, p_rr1, -alpha, 1.0_DP)

      ! Compute sa:=C(-1)*A*r
      call lsysbl_blockMatVec(p_rmatrix, p_rr, p_rsa, 1.0_DP, 0.0_DP)
      if (bprec) then
        call lsysbl_copyVector(p_rsa, p_rsa1)
        call linsol_precond(p_solverPrecond, p_rsa)
      end if

      ! Compute omega:=(sa,r)/(sa,sa)
      omega = lsysbl_scalarProduct(p_rsa,p_rr)/lsysbl_scalarProduct(p_rsa,p_rsa)

      ! Compute u:=u+alpha*p+omega*r
      call lsysbl_vectorLinearComb(p_rp, ru, alpha, 1.0_DP)
      call lsysbl_vectorLinearComb(p_rr, ru, omega, 1.0_DP)

      ! Compute r:=r-omega*sa
      call lsysbl_vectorLinearComb(p_rsa, p_rr, -omega, 1.0_DP)
      if (bprec) then
        call lsysbl_vectorLinearComb(p_rsa1, p_rr1, -omega, 1.0_DP)
        rsolver%dfinalDefect = lsysbl_vectorNorm(p_rr1, rsolver%iresNorm)
      else
        rsolver%dfinalDefect = lsysbl_vectorNorm(p_rr, rsolver%iresNorm)
      end if

      if (rsolver%ioutputLevel .ge. SV_IOLEVEL_VERBOSE) then
        call output_lbrk()
        call output_separator(OU_SEP_TILDE)
        call output_line('BiCGSTAB step:             '//trim(sys_siL(iiterations,5)))
        call output_line('Norm of residual:          '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
        call output_line('Norm of previous residual: '//trim(sys_sdEL(doldDefect,5)))
        call output_line('Variation of residual:     '//trim(sys_sdEL(&
                         rsolver%dfinalDefect/max(SYS_EPSREAL_DP, doldDefect),5)))
        call output_line('Improvement of residual:   '//trim(sys_sdEL(&
                         rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
        call output_separator(OU_SEP_TILDE)
        call output_lbrk()
      end if

      ! Check if residual increased too much
      if (solver_testDivergence(rsolver)) then
        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
          call output_line('!!! Residual increased by factor '//trim(sys_sdEL(&
                           rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
        end if

        ! Clear solution vector and adjust solver status
        call lsysbl_clearVector(ru)
        rsolver%istatus          = SV_INCR_DEF
        rsolver%dconvergenceRate = 1.0_DP

        ! That is it, return.
        return
      end if

      ! Check convergence criteria
      if (iiterations .ge. rsolver%nminIterations) then
        if (solver_testConvergence(rsolver)) then
          rsolver%istatus = SV_CONVERGED
          exit correction
        end if
        if (solver_testStagnation(rsolver, doldDefect)) then
          rsolver%istatus = SV_STAGNATED
          exit correction
        end if
      end if

      ! Save norm of old defect
      doldDefect = rsolver%dfinalDefect
    end do correction

    ! Compute convergence rate
    call solver_statistics(rsolver, iiterations)

    if (rsolver%ioutputLevel .ge. SV_IOLEVEL_INFO) then
      call output_lbrk()
      call output_separator(OU_SEP_PERC)
      call output_line('BiCGSTAB solution          '//solver_getstatus(rsolver))
      call output_line('Number of BiCGSTAB steps:  '//trim(sys_siL(rsolver%iiterations,5)))
      call output_line('Norm of final residual:    '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
      call output_line('Norm of initial residual:  '//trim(sys_sdEL(rsolver%dinitialDefect,5)))
      call output_line('Improvement of residual:   '//trim(sys_sdEL(&
                       rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
      call output_line('Convergence rate:          '//trim(sys_sdEL(rsolver%dconvergenceRate,5)))
      call output_separator(OU_SEP_PERC)
      call output_lbrk()
    end if

  end subroutine linsol_solveBicgstab

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_solveFgmres(rsolver, ru, rf)

!<description>
    ! This subroutine solves the linear system $Au=f$
    ! by means of the flexible GMRES(m) method
!</description>

!<input>
    ! Right-hand side vector
    type(t_vectorBlock), intent(in) :: rf
!</input>

!<inputoutput>
    ! Linear solver structure
    type(t_solver) :: rsolver

    ! Solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solverGMRES), pointer :: p_solverGMRES
    type(t_solver), pointer :: p_solverPrecond
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), dimension(:), pointer :: p_rv
    type(t_vectorBlock), dimension(:), pointer :: p_rz
    real(DP), dimension(:,:), pointer :: p_Dh
    real(DP), dimension(:), pointer :: p_Dc,p_Ds,p_Dq
    real(DP) :: dtmp,doldDefect
    integer :: neq,i,k,l,iiterations
    logical :: bprec


    ! Check for zero right-hand side
    if (rsolver%drhsZero > 0.0_DP) then
      rsolver%dinitialRHS = lsysbl_vectorNorm(rf, LINALG_NORMMAX)

      if (rsolver%dinitialRHS .le. rsolver%drhsZero) then
        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
          call output_line('!!! Zero initial right-hand side '//trim(&
                           sys_sdEL(rsolver%dinitialRHS,5))//' !!!')
        end if

        ! Clear solution vector and adjust solver status
        call lsysbl_clearVector(ru)
        rsolver%istatus          = SV_ZERO_RHS
        rsolver%dconvergenceRate = 0.0_DP

        ! That is it, return.
        return
      end if
    end if

    ! Set pointers
    p_solverGMRES  => rsolver%p_solverGMRES
    p_rmatrix      => p_solverGMRES%rmatrix
    p_rv           => p_solverGMRES%rv
    p_rz           => p_solverGMRES%rz; l=size(p_rz)
    call storage_getbase_double2D(p_solverGMRES%h_Dh, p_Dh)
    call storage_getbase_double(p_solverGMRES%h_Dc, p_Dc)
    call storage_getbase_double(p_solverGMRES%h_Ds, p_Ds)
    call storage_getbase_double(p_solverGMRES%h_Dq, p_Dq)

    ! Check for preconditioner
    if (associated(p_solverGMRES%p_precond)) then
      bprec = .true.
      p_solverPrecond => p_solverGMRES%p_precond
    else
      bprec = .false.
    end if

    ! Check compatibility
    call lsysbl_isVectorCompatible(ru, rf)
    call lsysbl_isMatrixCompatible(ru, p_rmatrix, .false.)

    ! Initialisation
    neq = ru%NEQ

    ! Compute initial residual v(1)=f-A*u
    call lsysbl_copyVector(rf, p_rv(1))
    call lsysbl_blockMatVec(p_rmatrix, ru, p_rv(1), -1.0_DP, 1.0_DP)

    ! Compute norm of initial defect
    rsolver%dinitialDefect = lsysbl_vectorNorm(p_rv(1), rsolver%iresNorm)
    rsolver%dfinalDefect   = rsolver%dinitialDefect
    doldDefect             = rsolver%dinitialDefect

    ! Check if initial residual is too large ...
    if (solver_testDivergence(rsolver)) then
      if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
        call output_line('!!! Norm of initial residual is too large '//&
                         trim(sys_sdEL(rsolver%dinitialDefect,5))//' !!!')
      end if

      ! Clear solution vector and adjust solver status
      call lsysbl_clearVector(ru)
      rsolver%istatus = SV_INF_DEF

      ! That is it, return.
      return

    elseif (rsolver%dinitialDefect .le. rsolver%ddefZero) then
      ! ... or if it satisfies the desired tolerance already
      if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
        call output_line('!!! Zero initial residual '//&
                         trim(sys_sdEL(rsolver%dinitialDefect,5))//' !!!')
      end if

      ! Clear solution vector and adjust solver status
      call lsysbl_clearVector(ru)
      rsolver%istatus = SV_ZERO_DEF

      ! That is it, return.
      return
    end if

    iiterations = 0
    outer: do while(iiterations .lt. rsolver%nmaxIterations)
      if (rsolver%dfinalDefect .lt. rsolver%depsAbs) exit outer

      ! Construct the elementary vector e1 scaled be residual norm
      call lalg_clearVectorDble(p_Dq)
      p_Dq(1) = rsolver%dfinalDefect

      ! Compute the first column of v(1)=r/!!r!!
      dtmp = 1.0_DP/p_Dq(1)
      call lsysbl_scaleVector(p_rv(1), dtmp)

      inner: do i = 1, l

        ! Increase number of iterations
        iiterations = iiterations+1

        ! Compute r=A*z(i) where z(i)=M(-1)*v(i)
        call lsysbl_copyVector(p_rv(i), p_rz(i))

        if (bprec) then
          call linsol_precond(p_solverPrecond, p_rz(i))
        end if

        ! Set v(i+1):=A*z(i)
        call lsysbl_copyVector(p_rv(i), p_rv(i+1))
        call lsysbl_blockMatVec(p_rmatrix, p_rz(i), p_rv(i+1), 1.0_DP, 0.0_DP)

        ! Construct i-th column of orthonormal basis h using
        ! a modified Gram-Schmidt algorithm
        do k = 1, i
          p_Dh(k,i) = lsysbl_scalarProduct(p_rv(i+1), p_rv(k))
          call lsysbl_vectorLinearComb(p_rv(k), p_rv(i+1), -p_Dh(k,i), 1.0_DP)
        end do
        P_Dh(i+1,i) = lsysbl_vectorNorm(p_rv(i+1), rsolver%iresNorm)

        ! Re-orthogonalize
        do k = 1, i
          dtmp = lsysbl_scalarProduct(p_rv(i+1), p_rv(k))
          p_Dh(k,i) = p_Dh(k,i)+dtmp
          call lsysbl_vectorLinearComb(p_rv(k), p_rv(i+1), -dtmp, 1.0_DP)
        end do
        p_Dh(i+1,i) = lsysbl_vectorNorm(p_rv(i+1), rsolver%iresNorm)

        if (abs(p_Dh(i+1,i)) > SYS_EPSREAL_DP) then
          dtmp = 1.0_DP/p_Dh(i+1,i)
          call lsysbl_scaleVector(p_rv(i+1), dtmp)
        end if

        ! Apply Givens rotations to the i-th column of h which
        ! renders the Hessenberg matrix an upper triangular matrix
        do k = 1, i-1
          dtmp = p_Dh(k,i)
          p_Dh(k,i)   = p_Dc(k)*dtmp - p_Ds(k)*p_Dh(k+1,i)
          p_Dh(k+1,i) = p_Ds(k)*dtmp + p_Dc(k)*p_Dh(k+1,i)
        end do

        ! Get next plane rotation
        dtmp = max(sqrt(p_Dh(i,i)**2 + p_Dh(i+1,i)**2), SYS_EPSREAL_DP)
        p_Dc(i) =  p_Dh(i,i)/dtmp
        p_Ds(i) = -p_Dh(i+1,i)/dtmp

        dtmp = p_Dq(i)
        p_Dq(i)   = p_Dc(i)*dtmp - p_Ds(i)*p_Dq(i+1)
        p_Dq(i+1) = p_Ds(i)*dtmp + p_Ds(i)*p_Dq(i+1)

        ! Determine residual norm and check convergence
        p_Dh(i,i) = p_Dc(i)*p_Dh(i,i) - p_Ds(i)*p_Dh(i+1,i)
        p_Dh(i+1,i) = 0.0_DP

        ! Check early convergence
        if (abs(p_Dq(i+1)) .le. rsolver%depsAbs) exit inner
      end do inner
      i = min(i, l)

      ! Solve the triangular system q:=h(-1)*q
      call DTRSV('U', 'N', 'N', i, p_Dh, l+1, p_Dq, 1)

      ! Update the solution u=u+p*z
      do k = 1, i
        call lsysbl_vectorLinearComb(p_rz(k), ru, p_Dq(k), 1.0_DP)
      end do

      ! Compute residual v(1)=f-A*u
      call lsysbl_copyVector(rf, p_rv(1))
      call lsysbl_blockMatVec(p_rmatrix, ru, p_rv(1), -1.0_DP, 1.0_DP)

      ! The norm of the residual is implicitly given by p(i+1) but
      ! may be highly inaccurate. Therefore compute norm explicitly
      rsolver%dfinalDefect = lsysbl_vectorNorm(p_rv(1), rsolver%iresNorm)

      if (rsolver%ioutputLevel .ge. SV_IOLEVEL_VERBOSE) then
        call output_lbrk()
        call output_separator(OU_SEP_TILDE)
        call output_line('FGMRES step:               '//trim(sys_siL(iiterations,5)))
        call output_line('Norm of residual:          '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
        call output_line('Norm of previous residual: '//trim(sys_sdEL(doldDefect,5)))
        call output_line('Variation of residual:     '//trim(sys_sdEL(&
                         rsolver%dfinalDefect/max(SYS_EPSREAL_DP, doldDefect),5)))
        call output_line('Improvement of residual:   '//trim(sys_sdEL(&
                         rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
        call output_separator(OU_SEP_TILDE)
        call output_lbrk()
      end if

      ! Check if residual increased too much
      if (solver_testDivergence(rsolver)) then
        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
          call output_line('!!! Residual increased by factor '//trim(sys_sdEL(&
                           rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
        end if

        ! Clear solution vector and adjust solver status
        call lsysbl_clearVector(ru)
        rsolver%istatus          = SV_INCR_DEF
        rsolver%dconvergenceRate = 1.0_DP

        ! That is it, return.
        return
      end if

      ! Check convergence criteria
      if (iiterations .ge. rsolver%nminIterations) then
        if (solver_testConvergence(rsolver)) then
          rsolver%istatus = SV_CONVERGED
          exit outer
        end if
        if (solver_testStagnation(rsolver, doldDefect)) then
          rsolver%istatus = SV_STAGNATED
          exit outer
        end if
      end if

      ! Save norm of old defect
      doldDefect = rsolver%dfinalDefect
    end do outer

    ! Compute convergence rate
    call solver_statistics(rsolver, iiterations)

    if (rsolver%ioutputLevel .ge. SV_IOLEVEL_INFO) then
      call output_lbrk()
      call output_separator(OU_SEP_PERC)
      call output_line('FGMRES solution            '//solver_getstatus(rsolver))
      call output_line('Number of FGMRES steps:    '//trim(sys_siL(rsolver%iiterations,5)))
      call output_line('Norm of final residual:    '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
      call output_line('Norm of initial residual:  '//trim(sys_sdEL(rsolver%dinitialDefect,5)))
      call output_line('Improvement of residual:   '//trim(sys_sdEL(&
                       rsolver%dfinalDefect/max(SYS_EPSREAL_DP, rsolver%dinitialDefect),5)))
      call output_line('Convergence rate:          '//trim(sys_sdEL(rsolver%dconvergenceRate,5)))
      call output_separator(OU_SEP_PERC)
      call output_lbrk()
    end if

  end subroutine linsol_solveFgmres

  ! *****************************************************************************

!<subroutine>

  subroutine linsol_precond(rsolver, ru)

!<description>
    ! This subroutine performs linear preconditioning of the given vector $u$.
    ! The given solver is applied as preconditioner.
!</description>

!<input>
    ! Linear Solver structure
    type(t_solver), intent(in) :: rsolver
!</input>

!<inputoutput>
    ! Vector to be preconditioned
    type(t_vectorBlock) :: ru
!</inputoutput>
!</subroutine>

    ! Check of solver is correct
    if (rsolver%csolverType .ne. SV_LINEAR) then
      call output_line('Unsupported solver type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linsol_precond')
      call sys_halt()
    end if

    ! What kind of solver should be applied?
    select case (rsolver%isolver)
    case(LINSOL_SOLVER_NONE)
      ! do nothing

    case(LINSOL_SOLVER_JACOBI)
      ! Jacobi preconditioner
      call linsol_precondJacobi(rsolver, ru)


    case(LINSOL_SOLVER_SOR, LINSOL_SOLVER_SSOR)
      ! (S)SOR preconditioner
      call linsol_precondSSOR(rsolver, ru)


    case(LINSOL_SOLVER_ILU)
      ! ILU preconditioner
      call linsol_precondILU(rsolver, ru)


    case default
      call output_line('Unsupported preconditioner type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linsol_precond')
      call sys_halt()
    end select
  end subroutine linsol_precond

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_precondJacobi(rsolver, ru)

!<description>
    ! This subroutine performs preconditioning by means of Jacobi iterations
!</description>

!<input>
    ! solver
    type(t_solver), intent(in) :: rsolver
!</input>

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixBlock),  pointer :: p_rmatrixBlock
    type(t_matrixScalar), pointer :: p_rmatrix
    type(t_vectorScalar), pointer :: p_rvector
    real(DP), dimension(:), pointer :: p_DA,p_Du
    integer, dimension(:), pointer :: p_Kdiagonal
    real(DP) :: domega
    integer :: iblock

    ! Set pointers
    p_rmatrixBlock  => rsolver%p_solverJacobi%rmatrix

    ! Set relaxation parameter
    if (rsolver%domega < 0.0_DP) then
      domega = 1.0_DP
    else
      domega = rsolver%domega
    end if

    ! Loop over all diagonal blocks
    do iblock = 1, ru%nblocks

      ! Set pointers
      p_rmatrix => p_rmatrixBlock%RmatrixBlock(iblock,iblock)
      p_rvector => ru%RvectorBlock(iblock)

      ! What kind of matrix are we?
      select case(p_rmatrix%cmatrixFormat)

      case (LSYSSC_MATRIXD)
        call lsyssc_getbase_double(p_rmatrix, p_DA)
        call lsyssc_getbase_double(p_rvector, p_Du)
        call jacobi_MatD_double(p_DA, p_Du, p_rvector%NEQ, domega)

      case (LSYSSC_MATRIX7)
        call lsyssc_getbase_Kld   (p_rmatrix, p_Kdiagonal)
        call lsyssc_getbase_double(p_rmatrix, p_DA)
        call lsyssc_getbase_double(p_rvector, p_Du)

        call jacobi_Mat79_double(p_Kdiagonal, p_DA, p_Du, p_rvector%NEQ, domega)

      case (LSYSSC_MATRIX9)
        call lsyssc_getbase_Kdiagonal(p_rmatrix, p_Kdiagonal)
        call lsyssc_getbase_double(p_rmatrix, p_DA)
        call lsyssc_getbase_double(p_rvector, p_Du)

        call jacobi_Mat79_double(p_Kdiagonal, p_DA, p_Du, p_rvector%NEQ, domega)

      case (LSYSSC_MATRIX7INTL)
        call lsyssc_getbase_Kld   (p_rmatrix, p_Kdiagonal)
        call lsyssc_getbase_double(p_rmatrix, p_DA)
        call lsyssc_getbase_double(p_rvector, p_Du)

        ! What kind of interleave matrix format are we?
        select case(p_rmatrix%cinterleavematrixFormat)

        case (LSYSSC_MATRIXD)
          call jacobi_Mat79IntlD_double(p_Kdiagonal,&
              p_rvector%NVAR, p_DA, p_Du, p_rvector%NEQ, domega)

        case (LSYSSC_MATRIX1)
          call jacobi_Mat79Intl1_double(p_Kdiagonal,&
              p_rvector%NVAR, p_DA, p_Du, p_rvector%NEQ, domega)

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'linsol_precondJacobi')
          call sys_halt()
        end select

      case (LSYSSC_MATRIX9INTL)
        call lsyssc_getbase_Kdiagonal(p_rmatrix, p_Kdiagonal)
        call lsyssc_getbase_double(p_rmatrix, p_DA)
        call lsyssc_getbase_double(p_rvector, p_Du)

        ! What kind of interleave matrix format are we?
        select case(p_rmatrix%cinterleavematrixFormat)

        case (LSYSSC_MATRIXD)
          call jacobi_Mat79IntlD_double(p_Kdiagonal,&
              p_rvector%NVAR, p_DA, p_Du, p_rvector%NEQ, domega)

        case (LSYSSC_MATRIX1)
          call jacobi_Mat79Intl1_double(p_Kdiagonal,&
              p_rvector%NVAR, p_DA, p_Du, p_rvector%NEQ, domega)

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'linsol_precondJacobi')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linsol_precondJacobi')
        call sys_halt()
      end select
    end do

  contains

    ! Here, the real working routines follow

    !*************************************************************
    ! This subroutine performs one Jacobi step for a double
    ! precision scalar matrix which is given as diagonal matrix

    subroutine jacobi_MatD_double(Da, Du, neq, domega)

      real(DP), dimension(:), intent(in) :: Da
      real(DP), dimension(:), intent(inout) :: Du
      integer, intent(in) :: neq
      real(DP), intent(in) :: domega

      integer :: ieq

      !$omp parallel do default(shared)
      do ieq = 1, neq
        Du(ieq) = domega*Du(ieq)/Da(ieq)
      end do
      !$omp end parallel do
    end subroutine jacobi_MatD_double

    !*************************************************************
    ! This subroutine performs one Jacobi step for a double
    ! precision scalar matrix which is stored in format 7 or 9

    subroutine jacobi_Mat79_double(Kdiagonal, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kdiagonal
      real(DP), dimension(:), intent(in) :: Da
      real(DP), dimension(:), intent(inout) :: Du
      integer, intent(in) :: neq
      real(DP), intent(in) :: domega

      integer :: ieq

      !$omp parallel do default(shared)
      do ieq = 1, neq
        Du(ieq) = domega*Du(ieq)/Da(Kdiagonal(ieq))
      end do
      !$omp end parallel do
    end subroutine jacobi_Mat79_double

    !*************************************************************
    ! This subroutine performs one Jacobi step for a double
    ! precision scalar interleave matrix which is stored in format
    ! 7 or 9, whereby each interleave matrix is a diagonal one

    subroutine jacobi_Mat79IntlD_double(Kdiagonal, nvar, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kdiagonal
      real(DP), dimension(nvar,*), intent(in) :: Da
      real(DP), dimension(nvar,*), intent(inout) :: Du
      integer, intent(in) :: neq
      integer, intent(in) :: nvar
      real(DP), intent(in) :: domega

      integer :: ieq,ivar

      !$omp parallel do default(shared) private(ivar)
      do ieq=1,neq
        do ivar=1,nvar
          Du(ivar,ieq)=domega*Du(ivar,ieq)/Da(ivar,Kdiagonal(ieq))
        end do
      end do
      !$omp end parallel do
    end subroutine jacobi_Mat79IntlD_double

    !*************************************************************
    ! This subroutine performs one Jacobi step for a double
    ! precision scalar interleave matrix which is stored in format
    ! 7 or 9, whereby each interleave matrix is a full one

    subroutine jacobi_Mat79Intl1_double(Kdiagonal,nvar,Da,Du,neq,domega)

      integer, dimension(:), intent(in) :: Kdiagonal
      real(DP), dimension(nvar,nvar,*), intent(in) :: Da
      real(DP), dimension(nvar,*), intent(inout) :: Du
      integer, intent(in) :: neq
      integer, intent(in) :: nvar
      real(DP), intent(in) :: domega

      integer :: ieq,ivar

      !$omp parallel do default(shared) private(ivar)
      do ieq=1,neq
        do ivar=1,nvar
          Du(ivar,ieq)=domega*Du(ivar,ieq)/Da(ivar,ivar,Kdiagonal(ieq))
        end do
      end do
      !$omp end parallel do
    end subroutine jacobi_Mat79Intl1_double
  end subroutine linsol_precondJacobi

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_precondSSOR(rsolver, ru)

!<description>
    ! This subroutine performs preconditioning by means of (S)SOR
    ! iterations
!</description>

!<input>
    ! solver
    type(t_solver), intent(in) :: rsolver
!</input>

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixBlock), pointer :: p_rmatrixBlock
    type(t_matrixScalar), pointer :: p_rmatrix
    type(t_vectorScalar), pointer :: p_rvector
    real(DP), dimension(:), pointer :: p_DA,p_Du
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    integer, dimension(:), pointer :: p_Kcol

    real(DP) :: domega
    integer :: iblock
    logical :: bsymmetric

    ! Initialisation
    bsymmetric = (rsolver%isolver .eq. LINSOL_SOLVER_SSOR)

    ! Set relaxation parameter
    if (rsolver%domega < 0.0_DP) then
      domega = 1.0_DP
    else
      domega = rsolver%domega
    end if

    ! Set pointers
    p_rmatrixBlock => rsolver%p_solverSSOR%rmatrix

    ! Loop over all blocks
    do iblock = 1, ru%nblocks

      ! Set pointers
      p_rmatrix => p_rmatrixBlock%RmatrixBlock(iblock,iblock)
      p_rvector => ru%RvectorBlock(iblock)

      ! What kind of matrix are we?
      select case(p_rmatrix%cmatrixFormat)

      case (LSYSSC_MATRIX7)
        call lsyssc_getbase_Kld   (p_rmatrix, p_Kld)
        call lsyssc_getbase_Kcol  (p_rmatrix, p_Kcol)
        call lsyssc_getbase_double(p_rmatrix, p_DA)
        call lsyssc_getbase_double(p_rvector, p_Du)

        if (bsymmetric) then
          call ssor_Mat7_double(p_Kld, p_Kcol,&
              p_DA, p_Du, p_rvector%NEQ, domega)
        else
          call sor_Mat7_double(p_Kld, p_Kcol,&
              p_DA, p_Du, p_rvector%NEQ, domega)
        end if

      case (LSYSSC_MATRIX9)
        call lsyssc_getbase_Kdiagonal(p_rmatrix, p_Kdiagonal)
        call lsyssc_getbase_Kld   (p_rmatrix, p_Kld)
        call lsyssc_getbase_Kcol  (p_rmatrix, p_Kcol)
        call lsyssc_getbase_double(p_rmatrix, p_DA)
        call lsyssc_getbase_double(p_rvector, p_Du)

        if (bsymmetric) then
          call ssor_Mat9_double(p_Kdiagonal, p_Kld,&
              p_Kcol, p_DA, p_Du, p_rvector%NEQ, domega)
        else
          call sor_Mat9_double(p_Kdiagonal, p_Kld,&
              p_Kcol, p_DA, p_Du, p_rvector%NEQ, domega)
        end if

      case (LSYSSC_MATRIX7INTL)
        call lsyssc_getbase_Kld   (p_rmatrix, p_Kld)
        call lsyssc_getbase_Kcol  (p_rmatrix, p_Kcol)
        call lsyssc_getbase_double(p_rmatrix, p_DA)
        call lsyssc_getbase_double(p_rvector, p_Du)

        ! What kind of interleave matrix are we?
        select case(p_rmatrix%cinterleavematrixFormat)

        case (LSYSSC_MATRIXD)
          if (bsymmetric) then
            call ssor_Mat7IntlD_double(p_Kld, p_Kcol,&
                p_rvector%NVAR, p_DA, p_Du, p_rvector%NEQ, domega)
          else
            call sor_Mat7IntlD_double(p_Kld, p_Kcol,&
                p_rvector%NVAR, p_DA, p_Du, p_rvector%NEQ, domega)
          end if

        case (LSYSSC_MATRIX1)
          if (bsymmetric) then
            call ssor_Mat7Intl1_double(p_Kld, p_Kcol,&
                p_rvector%NVAR, p_DA, p_Du, p_rvector%NEQ, domega)
          else
            call sor_Mat7Intl1_double(p_Kld, p_Kcol,&
                p_rvector%NVAR, p_DA, p_Du, p_rvector%NEQ, domega)
          end if

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'linsol_precondSSOR')
          call sys_halt()
        end select

      case (LSYSSC_MATRIX9INTL)
        call lsyssc_getbase_Kdiagonal(p_rmatrix, p_Kdiagonal)
        call lsyssc_getbase_Kld   (p_rmatrix, p_Kld)
        call lsyssc_getbase_Kcol  (p_rmatrix, p_Kcol)
        call lsyssc_getbase_double(p_rmatrix, p_DA)
        call lsyssc_getbase_double(p_rvector, p_Du)

        ! What kind of interleave matrix are we?
        select case(p_rmatrix%cinterleavematrixFormat)

        case (LSYSSC_MATRIXD)
          if (bsymmetric) then
            call ssor_Mat9IntlD_double(p_Kdiagonal, p_Kld, p_Kcol,&
                p_rvector%NVAR, p_DA, p_Du, p_rvector%NEQ, domega)
          else
            call sor_Mat9IntlD_double(p_Kdiagonal, p_Kld, p_Kcol,&
                p_rvector%NVAR, p_DA, p_Du, p_rvector%NEQ, domega)
          end if

        case (LSYSSC_MATRIX1)
          if (bsymmetric) then
            call ssor_Mat9Intl1_double(p_Kdiagonal, p_Kld, p_Kcol,&
                p_rvector%NVAR, p_DA, p_Du, p_rvector%NEQ, domega)
          else
            call sor_Mat9Intl1_double(p_Kdiagonal, p_Kld, p_Kcol,&
                p_rvector%NVAR, p_DA, p_Du, p_rvector%NEQ, domega)
          end if

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'linsol_precondSSOR')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linsol_precondSSOR')
        call sys_halt()
      end select
    end do

  contains

    ! Here, the real working routines follow

    !*************************************************************
    ! This subroutine performs one SOR step for a double
    ! precision scalar matrix which is stored in format 7

    subroutine sor_Mat7_double(Kld, Kcol, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), dimension(:), intent(in) :: Da
      real(DP), dimension(:), intent(inout) :: Du
      integer, intent(in) :: neq
      real(DP), intent(in) :: domega

      real(DP) :: daux
      integer :: ieq,icol,ild

      ! Process complete matrix from top-left to bottom-right
      do ieq = 1, neq
        daux = 0.0_DP

        ! Sum up row
        do ild = Kld(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)
          daux = daux+Da(ild)*Du(icol)
        end do

        ! Update solution vector
        Du(ieq)=(Du(ieq)-daux*domega)/Da(Kld(ieq))
      end do
    end subroutine sor_Mat7_double

    !*************************************************************
    ! This subroutine performs one SOR step for a double
    ! precision scalar matrix which is stored in format 9

    subroutine sor_Mat9_double(Kdiagonal, Kld, Kcol, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), dimension(:), intent(in) :: Da
      real(DP), dimension(:), intent(inout) :: Du
      integer, intent(in) :: neq
      real(DP), intent(in) :: domega

      real(DP) :: daux
      integer :: ieq,icol,ild

      ! Process complete matrix from top-left to bottom-right
      do ieq = 1, neq
        daux = 0.0_DP

        ! Sum up left part of row
        do ild = Kld(ieq), Kdiagonal(ieq)-1
          icol = Kcol(ild)
          daux = daux+Da(ild)*Du(icol)
        end do

        ! Sum up right part of row
        do ild = Kdiagonal(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)
          daux = daux+Da(ild)*Du(icol)
        end do

        ! Update solution vector
        Du(ieq)=(Du(ieq)-daux*domega)/Da(Kdiagonal(ieq))
      end do
    end subroutine sor_Mat9_double

    !*************************************************************
    ! This subroutine performs one SSOR step for a double
    ! precision scalar matrix which is stored in format 7

    subroutine ssor_Mat7_double(Kld, Kcol, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), dimension(:), intent(in) :: Da
      real(DP), dimension(:), intent(inout) :: Du
      integer, intent(in) :: neq
      real(DP), intent(in) :: domega

      real(DP) :: daux
      integer :: ieq,icol,ild

      ! Process lower left triangular matrix
      do ieq = 1, neq
        daux = 0.0_DP

        ! Sum up left part of row
        do ild = Kld(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)
          if (icol .ge. ieq) exit
          daux = daux+Da(ild)*Du(icol)
        end do

        ! Update solution vector
        Du(ieq) = (Du(ieq)-daux*domega)/Da(Kld(ieq))
      end do

      ! Process upper right triangular matrix
      do ieq = neq-1, 1, -1
        daux = 0.0_DP

        ! Sum up right part of row
        do ild = Kld(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)
          if (icol .le. ieq) cycle
          daux = daux+Da(ild)*Du(icol)
        end do

        ! Update solution vector
        Du(ieq) = Du(ieq)-daux*domega/Da(Kld(ieq))
      end do
    end subroutine ssor_Mat7_double

    !*************************************************************
    ! This subroutine performs one SSOR step for a double
    ! precision scalar matrix which is stored in format 9

    subroutine ssor_Mat9_double(Kdiagonal, Kld, Kcol, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), dimension(:), intent(in) :: Da
      real(DP), dimension(:), intent(inout) :: Du
      integer, intent(in) :: neq
      real(DP), intent(in) :: domega

      real(DP) :: daux
      integer :: ieq,icol,ild

      ! Process lower left triangular matrix
      do ieq = 1, neq
        daux = 0.0_DP

        ! Sum up left part of row
        do ild = Kld(ieq), Kdiagonal(ieq)-1
          icol = Kcol(ild)
          daux = daux+Da(ild)*Du(icol)
        end do

        ! Update solution vector
        Du(ieq) = (Du(ieq)-daux*domega)/Da(Kdiagonal(ieq))
      end do

      ! Process upper right triangular matrix
      do ieq = neq-1, 1, -1
        daux = 0.0_DP

        ! Sum up right part of row
        do ild = Kdiagonal(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)
          daux = daux+Da(ild)*Du(icol)
        end do

        ! Update solution vector
        Du(ieq)=Du(ieq)-daux*domega/Da(Kdiagonal(ieq))
      end do
    end subroutine ssor_Mat9_double

    !*************************************************************
    ! This subroutine performs one SOR step for a double
    ! precision scalar interleave matrix which is stored in format 7,
    ! whereby each intereave matrix is a diagonal one

    subroutine sor_Mat7IntlD_double(Kld, Kcol, nvar, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), dimension(nvar,*), intent(in) :: Da
      real(DP), dimension(nvar,*), intent(inout) :: Du
      integer, intent(in) :: neq
      integer, intent(in) :: nvar
      real(DP), intent(in) :: domega

      real(DP), dimension(nvar) :: daux
      integer :: ieq,icol,ild

      ! Process complete matrix from top-left to bottom-right
      do ieq = 1, neq
        daux = 0.0_DP

        ! Sum up row
        do ild = Kld(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)
          daux = daux+Da(:,ild)*Du(:,icol)
        end do

        ! Update solution vector
        Du(:,ieq) = (Du(:,ieq)-daux*domega)/Da(:,Kld(ieq))
      end do
    end subroutine sor_Mat7IntlD_double

    !*************************************************************
    ! This subroutine performs one SOR step for a double
    ! precision scalar interleave matrix which is stored in format 9,
    ! whereby each intereave matrix is a diagonal one

    subroutine sor_Mat9IntlD_double(Kdiagonal, Kld, Kcol, nvar, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), dimension(nvar,*), intent(in) :: Da
      real(DP), dimension(nvar,*), intent(inout) :: Du
      integer, intent(in) :: neq
      integer, intent(in) :: nvar
      real(DP), intent(in) :: domega

      real(DP), dimension(nvar) :: daux
      integer :: ieq,icol,ild

      ! Process complete matrix from top-left to bottom-right
      do ieq = 1,neq
        daux = 0.0_DP

        ! Sum up left part of row
        do ild = Kld(ieq), Kdiagonal(ieq)-1
          icol = Kcol(ild)
          daux = daux+Da(:,ild)*Du(:,icol)
        end do

        ! Sum up right part of row
        do ild = Kdiagonal(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)
          daux = daux+Da(:,ild)*Du(:,icol)
        end do

        ! Update solution vector
        Du(:,ieq)=(Du(:,ieq)-daux*domega)/Da(:,Kdiagonal(ieq))
      end do
    end subroutine sor_Mat9IntlD_double

    !*************************************************************
    ! This subroutine performs one SSOR step for a double
    ! precision scalar interleave matrix which is stored in format 7,
    ! whereby each interleave matrix is a diagonal one

    subroutine ssor_Mat7IntlD_double(Kld, Kcol, nvar, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), dimension(nvar,*), intent(in) :: Da
      real(DP), dimension(nvar,*), intent(inout) :: Du
      integer, intent(in) :: neq
      integer, intent(in) :: nvar
      real(DP), intent(in) :: domega

      real(DP), dimension(nvar) :: daux
      integer :: ieq,icol,ild

      ! Process lower left triangular matrix
      do ieq = 1, neq
        daux = 0.0_DP

        ! Sum up left part of row
        do ild = Kld(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)
          if (icol .ge. ieq) exit
          daux = daux+Da(:,ild)*Du(:,icol)
        end do

        ! Update solution vector
        Du(:,ieq)=(Du(:,ieq)-daux*domega)/Da(:,Kld(ieq))
      end do

      ! Process upper right triangular matrix
      do ieq = neq-1, 1, -1
        daux = 0.0_DP

        ! Sum up right part of row
        do ild = Kld(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)
          if (icol .le. ieq) cycle
          daux = daux+Da(:,ild)*Du(:,icol)
        end do

        ! Update solution vector
        Du(:,ieq) = Du(:,ieq)-daux*domega/Da(:,Kld(ieq))
      end do
    end subroutine ssor_Mat7IntlD_double

    !*************************************************************
    ! This subroutine performs one SSOR step for a double
    ! precision scalar interleave matrix which is stored in format 9,
    ! whereby each interleave matrix is a diagonal one

    subroutine ssor_Mat9IntlD_double(Kdiagonal, Kld, Kcol, nvar, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), dimension(nvar,*), intent(in) :: Da
      real(DP), dimension(nvar,*), intent(inout) :: Du
      integer, intent(in) :: neq
      integer, intent(in) :: nvar
      real(DP), intent(in) :: domega

      real(DP), dimension(nvar) :: daux
      integer :: ieq,icol,ild

      ! Process lower left triangular matrix
      do ieq = 1,neq
        daux = 0.0_DP

        ! Sum up left part of row
        do ild = Kld(ieq), Kdiagonal(ieq)-1
          icol = Kcol(ild)
          daux = daux+Da(:,ild)*Du(:,icol)
        end do

        ! Update solution vector
        Du(:,ieq) = (Du(:,ieq)-daux*domega)/Da(:,Kdiagonal(ieq))
      end do

      ! Process upper right triangular matrix
      do ieq = neq-1, 1, -1
        daux = 0.0_DP

        ! Sum up right part of row
        do ild = Kdiagonal(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)
          daux = daux+Da(:,ild)*Du(:,icol)
        end do

        ! Update solution vector
        Du(:,ieq) = Du(:,ieq)-daux*domega/Da(:,Kdiagonal(ieq))
      end do
    end subroutine ssor_Mat9IntlD_double

    !*************************************************************
    ! This subroutine performs one SOR step for a double
    ! precision scalar interleave matrix which is stored in format 7,
    ! whereby each intereave matrix is a full one

    subroutine sor_Mat7Intl1_double(Kld, Kcol, nvar, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), dimension(nvar,nvar,*), intent(in) :: Da
      real(DP), dimension(nvar,*), intent(inout) :: Du
      integer, intent(in) :: neq
      integer, intent(in) :: nvar
      real(DP), intent(in) :: domega

      real(DP), dimension(nvar) :: daux
      integer :: ieq,icol,ild,ivar,jvar

      ! Process complete matrix from top-left to bottom-right
      do ieq = 1, neq
        daux = 0.0_DP

        ! Sum up row
        do ild = Kld(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)

          do ivar = 1, nvar
            do jvar = 1, nvar
              daux(ivar) = daux(ivar)+Da(ivar,jvar,ild)*Du(jvar,icol)
            end do
          end do
        end do

        ! Update solution vector
        do ivar = 1,nvar
          Du(ivar,ieq) = (Du(ivar,ieq)-daux(ivar)*domega)/Da(ivar,ivar,Kld(ieq))
        end do
      end do
    end subroutine sor_Mat7Intl1_double

    !*************************************************************
    ! This subroutine performs one SOR step for a double
    ! precision scalar interleave matrix which is stored in format 9,
    ! whereby each intereave matrix is a full one

    subroutine sor_Mat9Intl1_double(Kdiagonal, Kld, Kcol, nvar, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), dimension(nvar,nvar,*), intent(in) :: Da
      real(DP), dimension(nvar,*), intent(inout) :: Du
      integer, intent(in) :: neq
      integer, intent(in) :: nvar
      real(DP), intent(in) :: domega

      real(DP), dimension(nvar) :: daux
      integer :: ieq,icol,ild,ivar,jvar

      ! Process complete matrix from top-left to bottom-right
      do ieq = 1, neq
        daux = 0.0_DP

        ! Sum up left part of row
        do ild = Kld(ieq), Kdiagonal(ieq)-1
          icol = Kcol(ild)

          do ivar=1,nvar
            do jvar=1,nvar
              daux(ivar)=daux(ivar)+Da(ivar,jvar,ild)*Du(jvar,icol)
            end do
          end do
        end do

        ! Sum up right part of row
        do ild = Kdiagonal(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)

          do ivar=1,nvar
            do jvar=1,nvar
              daux(ivar)=daux(ivar)+Da(ivar,jvar,ild)*Du(jvar,icol)
            end do
          end do
        end do

        ! Update solution vector
        do ivar = 1, nvar
          Du(ivar,ieq) = (Du(ivar,ieq)-daux(ivar)*domega)/Da(ivar,ivar,Kdiagonal(ieq))
        end do
      end do
    end subroutine sor_Mat9Intl1_double

    !*************************************************************
    ! This subroutine performs one SSOR step for a double
    ! precision scalar interleave matrix which is stored in format 7,
    ! whereby each interleave matrix is a full one

    subroutine ssor_Mat7Intl1_double(Kld, Kcol, nvar, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), dimension(nvar,nvar,*), intent(in) :: Da
      real(DP), dimension(nvar,*), intent(inout) :: Du
      integer, intent(in) :: neq
      integer, intent(in) :: nvar
      real(DP), intent(in) :: domega

      real(DP), dimension(nvar) :: daux
      integer :: ieq,icol,ild,ivar,jvar

      ! Process lower left triangular matrix
      do ieq = 1, neq
        daux = 0.0_DP

        ! Phase 1: Process diagonal blocks as in the scalar case
        do ild = Kld(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)
          if (icol .ge. ieq) exit

          do ivar = 1, nvar
            daux(ivar) = daux(ivar)+Da(ivar,ivar,ild)*Du(ivar,icol)
          end do
        end do

        ! Phase 2: Include sum of lower-left blocks
        do ild = Kld(ieq), Kld(ieq+1)-1
           icol = Kcol(ild)

          ! Loop over all row-blocks
          do ivar = 2, nvar

            ! Loop over all lower-left column-blocks
            do jvar = 1, ivar-1
              daux(ivar) = daux(ivar)+Da(ivar,jvar,ild)*Du(jvar,icol)
            end do
          end do
        end do

        ! Update solution vector
        do ivar = 1, nvar
          Du(ivar,ieq) = (Du(ivar,ieq)-daux(ivar)*domega)/Da(ivar,ivar,Kld(ieq))
        end do
      end do

      ! Process upper right triangular matrix
      do ieq = neq-1, 1, -1
        daux = 0.0_DP

        ! Phase 1: Process diagonal blocks as in the scalar case
        do ild = Kld(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)
          if (icol .le. ieq) cycle

          do ivar = 1,nvar
            daux(ivar) = daux(ivar)+Da(ivar,ivar,ild)*Du(ivar,icol)
          end do
        end do

        ! Phase 2: Include sum of upper-right blocks
        do ild = Kld(ieq),Kld(ieq+1)-1
          icol = Kcol(ieq)

          ! Loop over all row-blocks
          do ivar = 1, nvar

            ! Loop over all upper-right column-blocks
            do jvar = ivar+1, nvar
              daux(ivar) = daux(ivar)+Da(ivar,jvar,ild)*Du(jvar,icol)
            end do
          end do
        end do

        ! Update solution vector
        do ivar = 1, nvar
          Du(ivar,ieq) = Du(ivar,ieq)-daux(ivar)*domega/Da(ivar,ivar,Kld(ieq))
        end do
      end do
    end subroutine ssor_Mat7Intl1_double

    !*************************************************************
    ! This subroutine performs one SSOR step for a double
    ! precision scalar interleave matrix which is stored in format 9,
    ! whereby each interleave matrix is a full one

    subroutine ssor_Mat9Intl1_double(Kdiagonal, Kld, Kcol, nvar, Da, Du, neq, domega)

      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), dimension(nvar,nvar,*), intent(in) :: Da
      real(DP), dimension(nvar,*), intent(inout) :: Du
      integer, intent(in) :: neq
      integer, intent(in) :: nvar
      real(DP), intent(in) :: domega

      real(DP), dimension(nvar) :: daux
      integer :: ieq,icol,ild,ivar,jvar

      ! Process lower left triangular matrix
      do ieq = 1, neq
        daux = 0.0_DP

        ! Phase 1: Process diagonal blocks as in the scalar case
        do ild = Kld(ieq), Kdiagonal(ieq)-1
          icol = Kcol(ild)

          do ivar = 1,nvar
            daux(ivar) = daux(ivar)+Da(ivar,ivar,ild)*Du(ivar,icol)
          end do
        end do

        ! Phase 2: Include sum of lower-left blocks
        do ild = Kld(ieq), Kld(ieq+1)-1
           icol = Kcol(ild)

          ! Loop over all row-blocks
          do ivar = 2, nvar

            ! Loop over all lower-left column-blocks
            do jvar = 1, ivar-1
              daux(ivar) = daux(ivar)+Da(ivar,jvar,ild)*Du(jvar,icol)
            end do
          end do
        end do

        ! Update solution vector
        do ivar = 1, nvar
          Du(ivar,ieq) = (Du(ivar,ieq)-daux(ivar)*domega)/Da(ivar,ivar,Kdiagonal(ieq))
        end do
      end do

      ! Process upper right triangular matrix
      do ieq = neq-1, 1, -1
        daux = 0.0_DP

        ! Phase 1: Process diagonal blocks as in the scalar case
        do ild = Kdiagonal(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ild)

          do ivar = 1, nvar
            daux(ivar) = daux(ivar)+Da(ivar,ivar,ild)*Du(ivar,icol)
          end do
        end do

        ! Phase 2: Include sum of upper-right blocks
        do ild = Kld(ieq), Kld(ieq+1)-1
          icol = Kcol(ieq)

          ! Loop over all row-blocks
          do ivar = 1, nvar

            ! Loop over all upper-right column-blocks
            do jvar = ivar+1, nvar
              daux(ivar) = daux(ivar)+Da(ivar,jvar,ild)*Du(jvar,icol)
            end do
          end do
        end do

        ! Update solution vector
        do ivar = 1, nvar
          Du(ivar,ieq) = Du(ivar,ieq)-daux(ivar)*domega/Da(ivar,ivar,Kdiagonal(ieq))
        end do
      end do
    end subroutine ssor_Mat9Intl1_double
  end subroutine linsol_precondSSOR

  ! ***************************************************************************

!<subroutine>

  recursive subroutine linsol_precondILU(rsolver, ru)

!<description>
    ! This subroutine performs preconditioning by means of ILU
!</description>

!<input>
    ! solver
    type(t_solver), intent(in) :: rsolver
!</input>

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! Declare our LUSOLT-routine from SPLIB as interface to be sure,
    ! parameter interfaces are checked by compiler
    interface
      subroutine lusolt(n,x,lu,jlu,uptr)
        use fsystem

        integer(I32) jlu(*),uptr(*),n
        real(DP)     x(n)

        ! Note that we changed the interface here in contrast to the original
        ! LUSOLT routine - to make it possible to pass an integer array as
        ! double precision array. Bad practise, but SPLIB is set up this way.
        integer(I32) :: lu(*)
      end subroutine lusolt
    end interface

    ! local variables
    type(t_solverILU), pointer :: p_solver
    integer :: iblock


    ! Set pointer
    p_solver => rsolver%p_solverILU

    ! Do we have a block vector?
    if (ru%nblocks .ne. 1) then

      ! Check if vector and preconditioner are compatible
      if (associated(p_solver%p_rsolverBlockILU)) then
        if (size(p_solver%p_rsolverBlockILU,1) .ne. ru%nblocks) then
          call output_line('Block-ILU preconditioner and vector are not compatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'linsol_precondILU')
        end if
      else
        call output_line('Block-ILU preconditioner is not initialised!',&
            OU_CLASS_ERROR,OU_MODE_STD,'linsol_precondILU')
        call sys_halt()
      end if

      ! Loop over all blocks
      do iblock = 1, ru%nblocks
        call do_scalarILU(p_solver%p_rsolverBlockILU(iblock), ru%RvectorBlock(iblock))
      end do

    else

      call do_scalarILU(p_solver, ru%RvectorBlock(1))

    end if

  contains

    ! Here, the working routine follows

    !**************************************************************
    ! Perform (M)ILU(s) preconditioning for a scalar vector

    subroutine do_scalarILU(rsolver, ru)

      type(t_solverILU), intent(in) :: rsolver
      type(t_vectorScalar), intent(inout) :: ru

      type(t_matrixScalar), pointer :: p_rmatrix
      real(DP), dimension(:), pointer :: p_Du,p_Ddata,p_Da
      integer, dimension(:), pointer :: p_Kld
      integer, dimension(:), pointer :: p_Kcol,p_Kdiagonal
      integer(I32), dimension(:), pointer :: p_Idata,p_lu,p_jlu,p_ilup


      ! What kind of ILU algorithm should be used?
      if (rsolver%ifill .le. 0) then

        ! Set pointer to scalar matrix
        p_rmatrix  => rsolver%rmatrix%RmatrixBlock(1,1)

        ! What kind of matrix are we?
        select case(p_rmatrix%cmatrixFormat)
        case (LSYSSC_MATRIX7)
          call lsyssc_getbase_Kld    (p_rmatrix, p_Kld)
          call lsyssc_getbase_Kcol   (p_rmatrix, p_Kcol)
          call lsyssc_getbase_double (ru, p_Du)
          call storage_getbase_double(rsolver%h_Ddata, p_Ddata)

          call do_Mat7MILU0(p_Ddata, p_Kld, p_Kcol,&
                            p_rmatrix%NEQ, p_rmatrix%NVAR, p_Du)


        case (LSYSSC_MATRIX9)
          call lsyssc_getbase_Kld      (p_rmatrix, p_Kld)
          call lsyssc_getbase_Kcol     (p_rmatrix, p_Kcol)
          call lsyssc_getbase_Kdiagonal(p_rmatrix, p_Kdiagonal)
          call lsyssc_getbase_double   (ru, p_Du)
          call storage_getbase_double  (rsolver%h_Ddata, p_Ddata)

          call do_Mat9MILU0(p_Ddata, p_Kld, p_Kcol, p_Kdiagonal,&
                            p_rmatrix%NEQ, p_rmatrix%NVAR, p_Du)


        case (LSYSSC_MATRIX7INTL)

          ! What kind of interleave format are we?
          select case(p_rmatrix%cinterleavematrixFormat)

          case (LSYSSC_MATRIXD)
            call lsyssc_getbase_Kld    (p_rmatrix, p_Kld)
            call lsyssc_getbase_Kcol   (p_rmatrix, p_Kcol)
            call lsyssc_getbase_double (ru, p_Du)
            call storage_getbase_double(rsolver%h_Ddata, p_Ddata)

            call do_Mat7MILU0(p_Ddata, p_Kld, p_Kcol,&
                              p_rmatrix%NEQ, p_rmatrix%NVAR, p_Du)

          case (LSYSSC_MATRIX1)
            call lsyssc_getbase_Kld    (p_rmatrix, p_Kld)
            call lsyssc_getbase_Kcol   (p_rmatrix, p_Kcol)
            call lsyssc_getbase_double (p_rmatrix, p_Da)
            call lsyssc_getbase_double (ru, p_Du)
            call storage_getbase_double(rsolver%h_Ddata, p_Ddata)

            call do_Mat7Intl1BILU0(p_Da, p_Ddata, p_Kld, p_Kcol,&
                                   p_rmatrix%NEQ, p_rmatrix%NVAR, p_Du)


          case default
            call output_line('Unsupported interleave matrix format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'linsol_precondILU')
            call sys_halt()
          end select


        case (LSYSSC_MATRIX9INTL)

          ! What kind of interleave format are we?
          select case(p_rmatrix%cinterleavematrixFormat)

          case (LSYSSC_MATRIXD)
            call lsyssc_getbase_Kld      (p_rmatrix, p_Kld)
            call lsyssc_getbase_Kcol     (p_rmatrix, p_Kcol)
            call lsyssc_getbase_Kdiagonal(p_rmatrix, p_Kdiagonal)
            call lsyssc_getbase_double   (ru, p_Du)
            call storage_getbase_double  (rsolver%h_Ddata, p_Ddata)

            call do_Mat9MILU0(p_Ddata, p_Kld, p_Kcol, p_Kdiagonal,&
                              p_rmatrix%NEQ, p_rmatrix%NVAR, p_Du)


          case (LSYSSC_MATRIX1)
            call lsyssc_getbase_Kld      (p_rmatrix, p_Kld)
            call lsyssc_getbase_Kcol     (p_rmatrix, p_Kcol)
            call lsyssc_getbase_Kdiagonal(p_rmatrix, p_Kdiagonal)
            call lsyssc_getbase_double   (p_rmatrix, p_Da)
            call lsyssc_getbase_double   (ru, p_Du)
            call storage_getbase_double  (rsolver%h_Ddata, p_Ddata)

            call do_Mat9Intl1BILU0(p_Da, p_Ddata, p_Kld, p_Kcol,&
                                   p_Kdiagonal, p_rmatrix%NEQ, p_rmatrix%NVAR, p_Du)


          case default
            call output_line('Unsupported interleave matrix format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'linsol_precondILU')
            call sys_halt()
          end select


        case default
          call output_line('Unsupported matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'linsol_precondILU')
          call sys_halt()
        end select

      else

        ! Set pointers
        call lsyssc_getbase_double(ru, p_Du)
        call storage_getbase_int(rsolver%h_Idata, p_Idata)
        p_lu   => p_Idata(rsolver%lu:)
        p_jlu  => p_Idata(rsolver%jlu:)
        p_ilup => p_Idata(rsolver%ilup:)

        ! Solve the system by calling SPLIB
        call LUSOLT(int(size(p_Du), I32), p_Du, p_lu, p_jlu, p_ilup)

      end if
    end subroutine do_scalarILU


    !**************************************************************
    ! Perform (M)ILU(0) preconditioning for a matrix in format CSR7

    subroutine do_Mat7MILU0(Da, Kld, Kcol, neq, nvar, Du)

      real(DP), dimension(nvar,*), intent(in) :: Da
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      integer, intent(in) :: neq
      integer, intent(in) :: nvar

      real(DP), dimension(nvar,*), intent(inout) :: Du

      real(DP), dimension(nvar) :: Daux
      integer :: ia,ild,ieq,icol


      ! Forward substitution: L*y = b
      forward: do ieq = 1, neq

        Daux = 0.0_DP

        fwdin: do ia = Kld(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ia)
          if (icol .ge. ieq) exit fwdin
          Daux = Daux+Da(:,ia)*Du(:,icol)
        end do fwdin

        Du(:,ieq) = Du(:,ieq)-Daux

      end do forward

      ! Scale last entry by diagonal
      Daux = Da(:,Kld(neq))
      Du(:,neq) = Du(:,neq)/Daux

      ! Backward substitution: U*x = y
      backward: do ieq = neq-1, 1, -1

        ild = Kld(ieq)
        Daux = 0.0_DP

        bwdin: do ia = Kld(ieq+1)-1, ild+1, -1
          icol = Kcol(ia)
          if (icol .le. ieq) exit bwdin
          Daux = Daux+Da(:,ia)*Du(:,icol)
        end do bwdin

        Du(:,ieq)=(Du(:,ieq)-Daux)/Da(:,ild)

      end do backward
    end subroutine do_Mat7MILU0


    !**************************************************************
    ! Perform (M)ILU(0) peconditioning for a matrix in format CSR9

    subroutine do_Mat9MILU0(Da, Kld, Kcol, Kdiagonal, neq, nvar, Du)

      real(DP), dimension(nvar,*), intent(in) :: Da
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: neq
      integer, intent(in) :: nvar

      real(DP), dimension(nvar,*), intent(inout) :: Du

      real(DP), dimension(nvar) :: Daux
      integer :: ia,ild,ieq,icol


      ! Forward substitution
      forward: do ieq = 1, neq

        Daux = 0.0_DP

        fwdin: do ia = Kld(ieq), Kdiagonal(ieq)-1
          icol = Kcol(ia)
          Daux = Daux+Da(:,ia)*Du(:,icol)
        end do fwdin

        Du(:,ieq) = Du(:,ieq)-Daux

      end do forward

      ! Scale last entry by diagonal
      Daux = Da(:,Kdiagonal(neq))
      Du(:,neq) = Du(:,neq)/Daux

      ! Backward substitution
      backward: do ieq = neq-1, 1, -1

        ild = Kdiagonal(ieq)
        Daux = 0.0_DP

        bwdin: do ia = Kld(ieq+1)-1, Kdiagonal(ieq)+1, -1
          icol = Kcol(ia)
          Daux = Daux+Da(:,ia)*Du(:,icol)
        end do bwdin

        Du(:,ieq) = (Du(:,ieq)-Daux)/Da(:,ild)

      end do backward
    end subroutine do_Mat9MILU0


    !**************************************************************
    ! Perform BILU(0) peconditioning for a matrix in format CSR7,
    ! whereby the local blocks are full matrices

    subroutine do_mat7Intl1BILU0(Da, DaDiag, Kld, Kcol, neq, nvar, Du)

      real(DP), dimension(nvar,nvar,*), intent(in) :: Da
      real(DP), dimension(nvar,nvar,neq), intent(in) :: DaDiag
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      integer, intent(in) :: neq
      integer, intent(in) :: nvar

      real(DP), dimension(nvar,*), intent(inout) :: Du

      real(DP), dimension(nvar) :: DauxBlock
      real(DP) :: daux
      integer :: ia,ild,ieq,icol,ivar,jvar,ii

      ! Forward block-substitution: [D+Lower(A)]*y = b,
      ! whereby D is the LU decomposition of the diagonal blocks
      ! and Lower(A) stands for the lower triangular parts of A

      ! Compute the inverse of the first diagonal block
      ii = 0
      do ivar = 1, nvar
        daux = Du(ivar,1)
        if (ii .ne. 0) then
          do jvar = ii, ivar-1
            daux = daux-DaDiag(ivar,jvar,1)*Du(jvar,1)
          end do
        elseif (daux.ne.0.0_DP) then
          ii = ivar
        end if
        Du(ivar,1) = daux
      end do

      do ivar = nvar, 1, -1
        daux = Du(ivar,1)
        do jvar = ivar+1, nvar
          daux = daux-DaDiag(ivar,jvar,1)*Du(jvar,1)
        end do
        Du(ivar,1) = daux/DaDiag(ivar,ivar,1)
      end do

      ! Loop over rows 2,3,...,NEQ
      do ieq = 2, neq

        DauxBlock = 0.0_DP
        fwdin: do ia = Kld(ieq)+1, Kld(ieq+1)-1
          icol = Kcol(ia)
          if (icol .ge. ieq) exit fwdin
          DauxBlock = DauxBlock+matmul(Da(:,:,ia),Du(:,icol))
        end do fwdin

        Du(:,ieq)= Du(:,ieq)-DauxBlock

        ! Compute the inverse of the diagonal block
        ii = 0
        do ivar = 1, nvar
          daux = Du(ivar,ieq)
          if (ii .ne. 0) then
            do jvar = ii, ivar-1
              daux = daux-DaDiag(ivar,jvar,ieq)*Du(jvar,ieq)
            end do
          elseif (daux.ne.0.0_DP) then
            ii = ivar
          end if
          Du(ivar,ieq) = daux
        end do

        do ivar = nvar, 1, -1
          daux = Du(ivar,ieq)
          do jvar = ivar+1, nvar
            daux = daux-DaDiag(ivar,jvar,ieq)*Du(jvar,ieq)
          end do
          Du(ivar,ieq) = daux/DaDiag(ivar,ivar,ieq)
        end do
      end do
    end subroutine do_mat7Intl1BILU0

    !**************************************************************
    ! Perform BILU(0) peconditioning for a matrix in format CSR9,
    ! whereby the local blocks are full matrices

    subroutine do_mat9Intl1BILU0(Da, DaDiag, Kld, Kcol, Kdiagonal, neq, nvar, Du)

      real(DP), dimension(nvar,nvar,*), intent(in) :: Da
      real(DP), dimension(nvar,nvar,neq), intent(in) :: DaDiag
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, intent(in) :: neq
      integer, intent(in) :: nvar

      real(DP), dimension(nvar,*), intent(inout) :: Du

      real(DP), dimension(nvar) :: DauxBlock
      real(DP) :: daux
      integer :: ia,ild,ieq,icol,ivar,jvar,ii

      ! Forward block-substitution: [D+Lower(A)]*y = b,
      ! whereby D is the LU decomposition of the diagonal blocks
      ! and Lower(A) stands for the lower triangular parts of A

      ! Compute the inverse of the first diagonal block
      ii = 0
      do ivar = 1, nvar
        daux = Du(ivar,1)
        if (ii .ne. 0) then
          do jvar = ii, ivar-1
            daux = daux-DaDiag(ivar,jvar,1)*Du(jvar,1)
          end do
        elseif (daux.ne.0.0_DP) then
          ii = ivar
        end if
        Du(ivar,1) = daux
      end do

      do ivar = nvar, 1, -1
        daux = Du(ivar,1)
        do jvar = ivar+1, nvar
          daux = daux-DaDiag(ivar,jvar,1)*Du(jvar,1)
        end do
        Du(ivar,1) = daux/DaDiag(ivar,ivar,1)
      end do

      ! Loop over rows 2,3,...,NEQ
      do ieq = 2, neq

        DauxBlock = 0.0_DP
        fwdin: do ia = Kld(ieq), Kdiagonal(ieq)-1
          icol = Kcol(ia)
          DauxBlock = DauxBlock+matmul(Da(:,:,ia),Du(:,icol))
        end do fwdin

        Du(:,ieq)= Du(:,ieq)-DauxBlock

        ! Compute the inverse of the diagonal block
        ii = 0
        do ivar = 1, nvar
          daux = Du(ivar,ieq)
          if (ii .ne. 0) then
            do jvar = ii, ivar-1
              daux = daux-DaDiag(ivar,jvar,ieq)*Du(jvar,ieq)
            end do
          elseif (daux.ne.0.0_DP) then
            ii = ivar
          end if
          Du(ivar,ieq) = daux
        end do

        do ivar = nvar, 1, -1
          daux = Du(ivar,ieq)
          do jvar = ivar+1, nvar
            daux = daux-DaDiag(ivar,jvar,ieq)*Du(jvar,ieq)
          end do
          Du(ivar,ieq) = daux/DaDiag(ivar,ivar,ieq)
        end do
      end do
    end subroutine do_mat9Intl1BILU0
  end subroutine linsol_precondILU

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_smooth(rproblemLevel, rsolver, ru, rf, nsmooth)

!<description>
    ! Given an initial solution $U$ and a right-hand side $F$, this
    ! subroutine perform a prescribed number of smoothing steps.
!</description>

!<input>
    ! multigrid level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! right-hand side vector
    type(t_vectorBlock), intent(in) :: rf

    ! number of smoothing steps
    integer, intent(in) :: nsmooth
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! Check if solver is correct
    if (rsolver%csolverType .ne. SV_LINEAR) then
      call output_line('Unsupported solver type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linsol_smootherBlock')
      call sys_halt()
    end if

    ! What kind of smoother should be applied?
    select case(rsolver%isolver)
    case(LINSOL_SOLVER_JACOBI)
      ! Jacobi smoother
      rsolver%nminIterations = nsmooth
      rsolver%nmaxIterations = nsmooth

      call linsol_smoothJacobi(rsolver, ru, rf)

      rsolver%iiterations = nsmooth
      rsolver%niterations = rsolver%niterations + nsmooth


    case(LINSOL_SOLVER_SOR,&
         LINSOL_SOLVER_SSOR)
      ! (S)SOR smoother
      rsolver%nminIterations = nsmooth
      rsolver%nmaxIterations = nsmooth

      call linsol_smoothSSOR(rsolver, ru, rf)

      rsolver%iiterations = nsmooth
      rsolver%niterations = rsolver%niterations + nsmooth

    case(LINSOL_SOLVER_ILU)
      ! ILU(0) smoother
      rsolver%nminIterations = nsmooth
      rsolver%nmaxIterations = nsmooth

      call linsol_smoothILU(rsolver, ru, rf)

      rsolver%iiterations = nsmooth
      rsolver%niterations = rsolver%niterations + nsmooth

    case default
      call output_line('Unsupported smoother type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'linsol_smooth')
      call sys_halt()
    end select
  end subroutine linsol_smooth

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_smoothJacobi(rsolver, ru, rf)

!<description>
    ! This subroutine performs smoothing by means of Jacobi iterations
!</description>

!<input>
    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! r.h.s. vector
    type(t_vectorBlock), intent(in) :: rf
!</input>

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solverJacobi), pointer :: rsolverJacobi
    type(t_matrixBlock),  pointer :: rmatrix
    type(t_vectorBlock),  pointer :: rres
    integer :: iiterations

    ! Set pointers
    rsolverJacobi => rsolver%p_solverJacobi
    rmatrix       => rsolverJacobi%rmatrix
    rres          => rsolverJacobi%rtempVector

    ! Compute initial residual
    call lsysbl_copyVector(rf, rres)
    call lsysbl_blockMatVec(rmatrix, ru, rres, -1.0_DP, 1.0_DP)

    ! Iterative correction
    do iiterations = 1, rsolver%nmaxIterations

      ! Precondition the linear residual
      call linsol_precondJacobi(rsolver, rres)

      ! Update solution
      call lsysbl_vectorLinearComb(rres, ru, 1.0_DP, 1.0_DP)

      ! Compute residual
      call lsysbl_copyVector(rf, rres)
      call lsysbl_blockMatVec(rmatrix, ru, rres, -1.0_DP, 1.0_DP)

    end do
  end subroutine linsol_smoothJacobi

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_smoothSSOR(rsolver, ru, rf)

!<description>
    ! This subroutine performs smoothing by means of (S)SOR iterations
!</description>

!<input>
    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! r.h.s. vector
    type(t_vectorBlock), intent(in) :: rf
!</input>

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solverSSOR), pointer :: rsolverSSOR
    type(t_matrixBlock),  pointer :: rmatrix
    type(t_vectorBlock),  pointer :: rres
    integer :: iiterations

    ! Set pointers
    rsolverSSOR => rsolver%p_solverSSOR
    rmatrix     => rsolverSSOR%rmatrix
    rres        => rsolverSSOR%rtempVector

    ! Compute initial residual
    call lsysbl_copyVector(rf, rres)
    call lsysbl_blockMatVec(rmatrix, ru, rres, -1.0_DP, 1.0_DP)

    ! Iterative correction
    do iiterations = 1, rsolver%nmaxIterations

      ! Precondition the linear residual
      call linsol_precondSSOR(rsolver, rres)

      ! Update solution
      call lsysbl_vectorLinearComb(rres, ru, 1.0_DP, 1.0_DP)

      ! Compute residual
      call lsysbl_copyVector(rf,rres)
      call lsysbl_blockMatVec(rmatrix, ru, rres, -1.0_DP, 1.0_DP)

    end do
  end subroutine linsol_smoothSSOR

  ! ***************************************************************************

!<subroutine>

  subroutine linsol_smoothILU(rsolver, ru, rf)

!<description>
    ! This subroutine performs smoothing by means of ILU factorisation
!</description>

!<input>
    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! r.h.s. vector
    type(t_vectorBlock), intent(in) :: rf
!</input>

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(inout) :: ru
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solverILU), pointer :: rsolverILU
    type(t_matrixBlock), pointer :: rmatrix
    type(t_vectorBlock),  pointer :: raux
    integer :: iiterations

    ! Set pointers
    rsolverILU => rsolver%p_solverILU
    rmatrix    => rsolverILU%rmatrix
    raux       => rsolverILU%rtempVector

    ! Check compatibility
    call lsysbl_isMatrixCompatible(ru, rmatrix, .false.)
    call lsysbl_isVectorCompatible(ru, rf)
    call lsysbl_isVectorCompatible(ru, raux)

    do iiterations = 1, rsolver%nmaxIterations

      call lsysbl_copyVector (rf, raux)
      call lsysbl_blockMatVec(rmatrix, ru, raux, -1.0_DP, 1.0_DP)
      call linsol_precondILU (rsolver, raux)
      call lsysbl_vectorLinearComb(raux, ru, rsolver%domega, 1.0_DP)

    end do
  end subroutine linsol_smoothILU
end module solverlinear

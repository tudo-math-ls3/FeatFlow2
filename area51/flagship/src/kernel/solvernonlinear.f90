!##############################################################################
!# ****************************************************************************
!# <name> solvernonlinear </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides basic routines for solving nonlinear
!# algebraic systems of the generic form
!#
!# $$F(U;\bar U)=g(U),$$
!#
!# where the given solution $\bar U$ can be either an initial guess
!# $U_0$ or the solution from the last time step $U^n$.
!# In the latter case, the new solution $U$ stands for the solution 
!# at time $t^{n+1}=t^n+\Delta t$.
!# The quantities $\bar U$ and $U$ can be the nodal solution vectors for
!# a scalar quantity, e.g., $U=[u]^T$, or for a vector-valued
!# quantity, e.g., $U=[u_1,u_2,\dots,u_N]^T$.
!#
!# A) The nonlinear system is solved by means of a fixed-point iteration
!#
!# $$U^{(m+1)}=U^{(m)}+\omega^{(m)}C(U^{(m)})^{-1}G(U^{(m)};\bar U)$$
!#
!# where $G(U^{(m)};\bar U)$ denotes the nonlinear residual vector
!# (see above) and $C(U^{(m)})$ stands for the preconditioner of the
!# nonlinear iteration.
!#  There are three different types of nonlinear preconditioners:
!#
!# 1.) Segregated/block-diagonal approach:
!#     The preconditioner is a block-diagonal matrix so that all
!#     nonlinear equations can be solved independently.
!#
!# 2.) Defect-correction approach
!#     The preconditioner is a full block-matrix so that the complete
!#     system of equations is solved simultaneously.
!#
!# 3.) inexact Newton approach
!#     The preconditioner is a full block-matrix so that the complete
!#     system of equations is solved simultaneously. In addition to
!#     the standard defect-correction approach, the provisional
!#     solution increment may be scaled so as to satisfy the
!#     sufficient decrease condition.
!#
!# B) The nonlinear system is solved by means of a nonlinear multigrid iteration
!#
!# The following routines are available:
!#
!# 1.) nlsol_solveMultigrid = nlsol_solveMultigridScalar /
!#                            nlsol_solveMultigridBlock
!#     -> Solve the nonlinear system by means of the nonlinear multigrid approach
!#
!# 2.) nlsol_solveFixedpoint = nlsol_solveFixedpointScalar / 
!#                             nlsol_solveFixedpointBlock
!#     -> Solve the nonlinear system by means of a fixed-point approach
!#        on a single grid level
!#
!# 3.) nlsol_backtracking
!#     -> Perform backtracking to globalize inexact Newton iteration
!#
!# 4.) nlsol_calcForcingTerm
!#     -> Compute the forcing term for inexact Newton iteration
!#
!# 5.) nlsol_checkStagnation
!#     -> Check if nonlinear solver stagnates
!#
!# 6.) nlsol_relaxTimestep
!#     -> Relax the current time step.
!# </purpose>
!##############################################################################

module solvernonlinear

  use boundaryfilter
  use collection
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use solver
  use solverlinear
  use storage
 
  implicit none

  private
  public :: nlsol_solveMultigrid
  public :: nlsol_solveFixedpoint

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

  interface nlsol_solveMultigrid
    module procedure nlsol_solveMultigridScalar
    module procedure nlsol_solveMultigridBlock
  end interface

  interface nlsol_solveFixedpoint
    module procedure nlsol_solveFixedpointScalar
    module procedure nlsol_solveFixedpointBlock
  end interface

contains

  ! *****************************************************************************
  ! *****************************************************************************
  ! *****************************************************************************

!<subroutine>

  subroutine nlsol_solveMultigridScalar(rproblemLevel, rtimestep, rsolver,&
                                        rsolution, rsolutionInitial,& 
                                        fcb_calcResidual, fcb_calcJacobian,&
                                        fcb_applyJacobian, fcb_setBoundary,&
                                        rcollection, rrhs)

!<description>
    ! This subroutine solves the nonlinear system
    ! $$  F(U;\bar U)=g(U)  $$
    ! by means of the nonlinear multigrid approach.
    ! Note that this subroutine serves as wrapper for scalar equations.
!</description>

!<input>
    ! Initial solution vector
    type(t_vectorScalar), intent(IN) :: rsolutionInitial

    ! Callback routines
    include 'intf_nlsolcallback.inc'

    ! OPTIONAL: right-hand side vector
    type(t_vectorScalar), intent(IN), optional :: rrhs
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver
    
    ! solution vector
    type(t_vectorScalar), intent(INOUT) :: rsolution

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rsolutionBlock
    type(t_vectorBlock) :: rsolutionInitialBlock
    type(t_vectorBlock) :: rrhsBlock

    if (present(rrhs)) then

      ! Convert scalar vectors into block vectors
      call lsysbl_createVecFromScalar(rrhs, rrhsBlock)
      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)
      call lsysbl_createVecFromScalar(rsolutionInitial, rsolutionInitialBlock)
      
      ! Apply block-version of nonlinear multigrid
      call nlsol_solveMultigridBlock(rproblemLevel, rtimestep, rsolver,&
                                     rsolutionBlock, rsolutionInitialBlock,&
                                     fcb_calcResidual, fcb_calcJacobian,&
                                     fcb_applyJacobian, fcb_setBoundary,&
                                     rcollection, rrhsBlock)
      
      ! Release temporal block vectors
      call lsysbl_releaseVector(rrhsBlock)
      call lsysbl_releaseVector(rsolutionBlock)
      call lsysbl_releaseVector(rsolutionInitialBlock)
      
    else
      
      ! Convert scalar vectors into block vectors
      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)
      call lsysbl_createVecFromScalar(rsolutionInitial, rsolutionInitialBlock)
      
      ! Apply block-version of nonlinear multigrid
      call nlsol_solveMultigridBlock(rproblemLevel, rtimestep, rsolver,&
                                     rsolutionBlock, rsolutionInitialBlock,&
                                     fcb_calcResidual, fcb_calcJacobian,&
                                     fcb_applyJacobian, fcb_setBoundary,&
                                     rcollection)
      
      ! Release temporal block vectors
      call lsysbl_releaseVector(rsolutionBlock)
      call lsysbl_releaseVector(rsolutionInitialBlock)

    end if
    
  end subroutine nlsol_solveMultigridScalar

  ! *****************************************************************************

!<subroutine>

  subroutine nlsol_solveMultigridBlock(rproblemLevel, rtimestep, rsolver,&
                                       rsolution, rsolutionInitial,&
                                       fcb_calcResidual, fcb_calcJacobian,&
                                       fcb_applyJacobian, fcb_setBoundary,&
                                       rcollection, rrhs)

!<description>
    ! This subroutine solves the nonlinear system
    ! $$  F(U;\bar U)=g(U)  $$
    ! by means of the nonlinear multigrid approach.
!</description>

!<input>
    ! Initial solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionInitial

    ! Callback routines
    include 'intf_nlsolcallback.inc'

    ! OPTIONAL: right-hand side vector
    type(t_vectorBlock), intent(IN), optional :: rrhs
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver
    
    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution
    
    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_solverMultigrid), pointer :: p_solverMultigrid
    type(t_solver), pointer :: p_solverNonlinear
    integer :: imgstep, mgcycle

    ! What kind of solver are we?
    select case(rsolver%csolverType)
      
    case (SV_NONLINEAR)
      !-------------------------------------------------------------------------
      ! Single grid solver: G(u)=f
      !-------------------------------------------------------------------------
      
      select case(rsolver%isolver)

      case (NLSOL_SOLVER_FIXEDPOINT)
        call nlsol_solveFixedPoint(rproblemLevel, rtimestep, rsolver,&
                                   rsolution, rsolutionInitial,&
                                   fcb_calcResidual, fcb_calcJacobian,&
                                   fcb_applyJacobian, fcb_setBoundary,&
                                   rcollection, rrhs)
        
      case DEFAULT
        call output_line('Invalid solver type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'nlsol_solveMultigridBlock')
        call sys_halt()
      end select
      

    case (SV_NONLINEARMG)
      
      ! Check if multigrid solver exists
      if (.not.associated(rsolver%p_solverMultigrid)) then
        call output_line('Multigrid solver does not exists!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'nlsol_solveMultigridBlock')
        call sys_halt()
      end if
      
      ! Set pointer
      p_solverMultigrid => rsolver%p_solverMultigrid

      ! Check if multigrid solver is called for only one grid level
      if (p_solverMultigrid%nlmin .eq. p_solverMultigrid%nlmax) then

        !-----------------------------------------------------------------------
        ! Single grid solver: G(u)=f
        !-----------------------------------------------------------------------

        ! Check if single-grid solver exists
        if (.not.associated(p_solverMultigrid%p_solverCoarsegrid)) then
          call output_line('Coarsegrid solver does not exists!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'nlsol_solveMultigridBlock')
          call sys_halt()
        end if

        ! Set pointer
        p_solverNonlinear => p_solverMultigrid%p_solverCoarsegrid

        ! What kind of solver are we?
        select case(p_solverNonlinear%isolver)

        case (NLSOL_SOLVER_FIXEDPOINT)
          call nlsol_solveFixedPoint(rproblemLevel, rtimestep, p_solverNonlinear,&
                                     rsolution, rsolutionInitial,&
                                     fcb_calcResidual, fcb_calcJacobian,&
                                     fcb_applyJacobian, fcb_setBoundary,&
                                     rcollection, rrhs)
          
        case DEFAULT
          call output_line('Invalid solver type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'nlsol_solveMultigridBlock')
          call sys_halt()
        end select
        

      else
        
        !-----------------------------------------------------------------------
        ! Multigrid grid solver: G(u)=f
        !-----------------------------------------------------------------------
        
        ! Perform prescribed number of multigrid steps
        mgstep: do imgstep = 1, p_solverMultigrid%ilmax
          
          ! Perform one nonlinear two-grid step
          mgcycle = merge(1, 2, p_solverMultigrid%icycle .eq. 1)
          call nlsol_solveTwogrid(rproblemLevel, rtimestep, rsolver,&
                                  rsolution, rsolutionInitial,&
                                  fcb_calcResidual, fcb_calcJacobian,&
                                  fcb_applyJacobian, fcb_setBoundary,&
                                  rcollection, rrhs)
          
          if (rsolver%ioutputLevel .ge. SV_IOLEVEL_VERBOSE) then
            call output_lbrk()
            call output_separator(OU_SEP_TILDE)
            call output_line('Nonlinear multigrid step: '//trim(sys_siL(imgstep,5)))
            call output_line('Norm of residual:         '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
            call output_line('Improvement of residual:  '//trim(sys_sdEL(rsolver%dfinalDefect/&
                                                          max(SYS_EPSREAL, rsolver%dinitialDefect),5)))
            call output_separator(OU_SEP_TILDE)
            call output_lbrk()
          end if

          ! Write nonlinear defect to file?
          if (rsolver%ioutputLevel .ge. SV_IOLEVEL_FILE) then
            write(UNIT=rsolver%iunitLogfile, FMT='("(01),",I10,2(",",E16.8E3))')&
                imgstep, rsolver%dfinalDefect,&
                rsolver%dfinalDefect/max(SYS_EPSREAL, rsolver%dinitialDefect)
          end if
          

          ! Check if residual increased too much
          if (rsolver%dfinalDefect > rsolver%ddivRel*rsolver%dinitialDefect .or.&
              rsolver%dfinalDefect > rsolver%ddivAbs) then
            if (rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
              call output_line('!!! Residual increased by factor '//trim(sys_sdEL(&
                  rsolver%dfinalDefect/max(SYS_EPSREAL, rsolver%dinitialDefect),5))//' !!!')
            end if

            ! Adjust solver status
            rsolver%istatus = SV_INCR_DEF

            ! That's it, return
            return
          end if

          ! Check convergence criteria
          if (imgstep .ge. p_solverMultigrid%ilmin) then
            if ((rsolver%dfinalDefect .le. rsolver%depsAbs) .or.&
                (rsolver%dfinalDefect .le. rsolver%depsRel*rsolver%dinitialDefect)) then
              exit mgstep
            end if
          end if
          
        end do mgstep

        ! Multigrid convergence rates
        call solver_statistics(rsolver, imgstep)

        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_INFO) then
          call output_lbrk()
          call output_separator(OU_SEP_TILDE)
          call output_line('Nonlinear multigrid solution '//solver_getstatus(rsolver))
          call output_line('Number of multigrid steps:   '//trim(sys_siL(rsolver%iiterations,5)))
          call output_line('Norm of final residual:      '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
          call output_line('Norm of initial residual:    '//trim(sys_sdEL(rsolver%dinitialDefect,5)))
          call output_line('Residual improvement:        '//trim(sys_sdEL(rsolver%dfinalDefect/&
                                                            max(SYS_EPSREAL, rsolver%dinitialDefect),5)))
          call output_line('Convergence rate:            '//trim(sys_sdEL(rsolver%dconvergenceRate,5)))
          call output_separator(OU_SEP_TILDE)
          call output_lbrk()
        end if
        
        ! Write nonlinear convergence rate to logfile?
        if (rsolver%ioutputLevel .ge. SV_IOLEVEL_FILE) then
          write(UNIT=rsolver%iunitLogfile,FMT='("(02),",I10,4(",",E16.8E3))')&
              rsolver%iiterations, rsolver%dfinalDefect, rsolver%dinitialDefect,&
              rsolver%dfinalDefect/max(SYS_EPSREAL, rsolver%dinitialDefect),&
              rsolver%dconvergenceRate
        end if
      end if
      

    case DEFAULT
      call output_line('Invalid solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'nlsol_solveMultigridBlock')
      call sys_halt()
    end select
  end subroutine nlsol_solveMultigridBlock

  ! *****************************************************************************

!<subroutine>

  recursive subroutine nlsol_solveTwogrid(rproblemLevel, rtimestep, rsolver,&
                                          rsolution, rsolutionInitial,&
                                          fcb_calcResidual, fcb_calcJacobian,&
                                          fcb_applyJacobian, fcb_setBoundary,&
                                          rcollection, rrhs)

!<description>
    ! This subroutine performs one nonlinear two-grid step of nonlinear multigrid.
!</description>

!<input>
    ! Initial solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionInitial

    ! Callback routines
    include 'intf_nlsolcallback.inc'

    ! OPTIONAL: right-hand side
    type(t_vectorBlock), intent(IN), optional :: rrhs
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver
    
    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    call output_line('Nonlinear two-grid solver is not implemented yet!',&
                     OU_CLASS_ERROR,OU_MODE_STD,'nlsol_solveTwogrid')
    call sys_halt()
  end subroutine nlsol_solveTwogrid
    
  ! *****************************************************************************

!<subroutine>

  subroutine nlsol_solveFixedpointScalar(rproblemLevel, rtimestep, rsolver,&
                                         rsolution, rsolutionInitial,&
                                         fcb_calcResidual, fcb_calcJacobian,&
                                         fcb_applyJacobian, fcb_setBoundary,&
                                         rcollection, rrhs)

!<description>
    ! This subroutine performs one nonlinear step to solve the system
    ! $$  dU/dT+F(U;\bar U)=g(U)  $$
    ! for a scalar quantity $U$. The solution at state $\bar U$ must be
    ! given. The solution at the new state $U$ is computed either by
    ! means of a fixed-point defect correction procedure or via an
    ! inexact Newton algorithm.
    !
    ! NOTE: This subroutine serves as a wrapper for scalar problems.
    ! The scalar vector $\bar U$ is converted into a 1-block vector and
    ! the block version if this subroutine is called.
!</description>

!<input>
    ! Initial solution vector
    type(t_vectorScalar), intent(IN) :: rsolutionInitial

    ! Callback routines
    include 'intf_nlsolcallback.inc'

    ! OPTIONAL: right-hand side
    type(t_vectorScalar), intent(IN), optional :: rrhs
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver
    
    ! solution vector
    type(t_vectorScalar), intent(INOUT) :: rsolution

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rrhsBlock
    type(t_vectorBlock) :: rsolutionBlock
    type(t_vectorBlock) :: rsolutionInitialBlock
    
    
    if (present(rrhs)) then

      ! Convert scalar vectors into block vectors
      call lsysbl_createVecFromScalar(rrhs, rrhsBlock)    
      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)
      call lsysbl_createVecFromScalar(rsolutionInitial, rsolutionInitialBlock)
      
      ! Call block version
      call nlsol_solveFixedpointBlock(rproblemLevel, rtimestep, rsolver,&
                                      rsolutionBlock, rsolutionInitialBlock,&
                                      fcb_calcResidual, fcb_calcJacobian,&
                                      fcb_applyJacobian, fcb_setBoundary,&
                                      rcollection, rrhsBlock)

      ! Release temporal block vectors
      call lsysbl_releaseVector(rrhsBlock)
      call lsysbl_releaseVector(rsolutionBlock)
      call lsysbl_releaseVector(rsolutionInitialBlock)
      
    else

      ! Convert scalar vectors into block vectors
      call lsysbl_createVecFromScalar(rsolution, rsolutionBlock)
      call lsysbl_createVecFromScalar(rsolutionInitial, rsolutionInitialBlock)
      
      ! Call block version
      call nlsol_solveFixedpointBlock(rproblemLevel, rtimestep, rsolver,&
                                      rsolutionBlock, rsolutionInitialBlock,&
                                      fcb_calcResidual, fcb_calcJacobian,&
                                      fcb_applyJacobian, fcb_setBoundary,&
                                      rcollection)

      ! Release temporal block vectors
      call lsysbl_releaseVector(rsolutionBlock)
      call lsysbl_releaseVector(rsolutionInitialBlock)
    end if

  end subroutine nlsol_solveFixedpointScalar

  ! ***************************************************************************

!<subroutine>
  
  subroutine nlsol_solveFixedpointBlock(rproblemLevel, rtimestep, rsolver,&
                                        rsolution, rsolutionInitial, &
                                        fcb_calcResidual, fcb_calcJacobian,&
                                        fcb_applyJacobian, fcb_setBoundary,&
                                        rcollection, rrhs)

!<description>
    ! This subroutine performs one nonlinear step to solve the system
    ! $$  dU/dT+F(U;\bar U)=g(U)  $$
    ! for a vector-valued quantity $U$. The solution at state $\bar U$
    ! must be given. The solution at the new state $U$ is computed
    ! either by means of a fixed-point defect correction procedure or
    ! via an inexact Newton algorithm.
!</description>

!<input>
    ! Initial solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionInitial

    ! Callback routines
    include 'intf_nlsolcallback.inc'

    ! OPTIONAL: right-hand side vector
    type(t_vectorBlock), intent(IN), optional :: rrhs
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep

    ! solver structure
    type(t_solver), intent(INOUT), target :: rsolver
        
    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_vectorBlock), pointer :: p_rrhs
    type(t_vectorBlock), pointer :: p_rres
    type(t_vectorBlock), pointer :: p_raux
    type(t_vectorBlock), pointer :: p_rufs
    type(t_vectorBlock), pointer :: p_rresfs
    type(t_solver), pointer :: p_rsolver,p_rsolverLinear
    type(t_solver) :: rsolverTemp
    real(DP) :: eta,drtlm,drtjs,redfac,doldDefect
    integer :: iiterations,iblock,ibacktrackingsteps,istatus
    logical :: bstagnate
    logical :: bcompatible


    ! Set pointer to nonlinear solver
    p_rsolver => solver_getNextSolverByType(rsolver, SV_NONLINEAR)
    
    if (.not. associated(p_rsolver)) then
      call output_line('Unsupported/invalid solver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'nlsol_solveFixedpointBlock')
      call sys_halt()
    end if
    
    ! Check if fixed-point iteration should be performed
    if (p_rsolver%isolver .ne. NLSOL_SOLVER_FIXEDPOINT) then
      call output_line('Invalid solver for fixed-point iteration!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'nlsol_solveFixedpointBlock')
      call sys_halt()
    end if

    ! Set pointer to linear solver
    p_rsolverLinear => p_rsolver%p_solverSubnode
    if (.not. associated(p_rsolverLinear)) then
      call output_line('No subsolver is associated!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'nlsol_solveFixedpointBlock')
      call sys_halt()
    end if

    ! Check if subsolver is a linear solver
    if ((p_rsolverLinear%csolverType .ne. SV_LINEARMG) .and.&
        (p_rsolverLinear%csolverType .ne. SV_LINEAR  )) then
      call output_line('Unsupported/invalid subsolver type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'nlsol_solveFixedpointBlock')
      call sys_halt()
    end if

    
    !---------------------------------------------------------------------------
    ! Initialization
    !---------------------------------------------------------------------------  
    
    select case(p_rsolver%iprecond)

    case(NLSOL_PRECOND_BLOCKD,&
         NLSOL_PRECOND_DEFCOR,&
         NLSOL_PRECOND_NEWTON_FAILED)

      ! Set pointers
      p_rrhs => p_rsolver%p_solverDefcor%rTempVectors(1)
      p_rres => p_rsolver%p_solverDefcor%rTempVectors(2)
      p_raux => p_rsolver%p_solverDefcor%rTempVectors(3)

      ! Check if compatible vectors are available.
      call lsysbl_isVectorCompatible(rsolution, p_rres, bcompatible)
      if (.not.bcompatible) then
        call lsysbl_resizeVectorBlock(p_rres, rsolution, .false.)
        call lsysbl_resizeVectorBlock(p_rrhs, rsolution, .false.)
        call lsysbl_resizeVectorBlock(p_raux, rsolution, .false.)
      end if

      ! Compute the initial residual and the constant right-hand side
      call fcb_calcResidual(rproblemLevel, rtimestep, p_rsolver,&
                            rsolution, rsolutionInitial, p_rrhs, p_rres, 0, rcollection)

      ! Apply given right-hand side vector
      if (present(rrhs)) call lsysbl_vectorLinearComb(rrhs, p_rres, 1.0_DP, 1.0_DP)
      
      ! Impose boundary conditions
      call fcb_setBoundary(rproblemLevel, rtimestep, p_rsolver,&
                           rsolution, p_rres, rsolutionInitial, rcollection)
      if (rcollection%Iquickaccess(1) < 0) then
        p_rsolver%istatus = SV_INF_DEF
        return
      end if

    
    case(NLSOL_PRECOND_NEWTON)
      
      ! Set pointers
      p_rrhs   => p_rsolver%p_solverNewton%rTempVectors(1)
      p_rres   => p_rsolver%p_solverNewton%rTempVectors(2)
      p_raux   => p_rsolver%p_solverNewton%rTempVectors(3)
      p_rufs   => p_rsolver%p_solverNewton%rTempVectors(4)
      p_rresfs => p_rsolver%p_solverNewton%rTempVectors(5)

      ! Check if compatible vectors are available.
      call lsysbl_isVectorCompatible(rsolution, p_rres, bcompatible)
      if (.not.bcompatible) then
        call lsysbl_resizeVectorBlock(p_rres,   rsolution, .false.)
        call lsysbl_resizeVectorBlock(p_rrhs,   rsolution, .false.)
        call lsysbl_resizeVectorBlock(p_raux,   rsolution, .false.)
        call lsysbl_resizeVectorBlock(p_rufs,   rsolution, .false.)
        call lsysbl_resizeVectorBlock(p_rresfs, rsolution, .false.)
      end if

      ! Compute the initial residual and the constant right-hand side
      call fcb_calcResidual(rproblemLevel, rtimestep, p_rsolver,&
                            rsolution, rsolutionInitial, p_rrhs, p_rres, 0, rcollection)

      ! Apply given right-hand side vector
      if (present(rrhs)) call lsysbl_vectorLinearComb(rrhs, p_rres, 1.0_DP, 1.0_DP)
   
      ! Assemble Jacobian matrix
      call fcb_calcJacobian(rproblemLevel, rtimestep, p_rsolver,&
                            rsolution, rsolutionInitial, .false., rcollection)
      
      ! Impose boundary conditions
      call fcb_setBoundary(rproblemLevel, rtimestep, p_rsolver,&
                           rsolution, p_rres, rsolutionInitial, rcollection)
      if (rcollection%Iquickaccess(1) < 0) then
        p_rsolver%istatus = SV_INF_DEF
        return
      end if


    case DEFAULT
      call output_line('Invalid nonlinear preconditioner!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'nlsol_solveFixedpointBlock')
      call sys_halt()
    end select


    ! Compute norm of initial residual
    p_rsolver%dinitialDefect = lsysbl_vectorNorm(p_rres, p_rsolver%iresNorm)
    p_rsolver%dfinalDefect   = p_rsolver%dinitialDefect
    doldDefect               = p_rsolver%dinitialDefect
    p_rsolver%istatus        = SV_CONVERGED
    
    ! Check if initial residual is too large ...
    if (p_rsolver%dinitialDefect > p_rsolver%ddivAbs) then
      if (p_rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
        call output_line('!!! Norm of initial residual is too large '//&
                         trim(sys_sdEL(p_rsolver%dinitialDefect,5))//' !!!')
      end if
      
      ! Adjust solver status
      p_rsolver%istatus = SV_INF_DEF
      
      ! Also adjust solver status of top-most solver
      rsolver%istatus = SV_INF_DEF
      
      ! That's it, return.
      return
      
    elseif (p_rsolver%dinitialDefect .le. p_rsolver%ddefZero) then
      ! ... or if it satisfies the desired tolerance already
      if (p_rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
        call output_line('!!! Zero initial residual '//&
                         trim(sys_sdEL(p_rsolver%dinitialDefect,5))//' !!!')
      end if
      
      ! Adjust solver status
      p_rsolver%istatus = SV_ZERO_DEF
      
      ! Also adjust solver status of top-most solver
      rsolver%istatus = SV_ZERO_DEF
      
      ! That's it, return.
      return
    end if

    ! Initialize check for stagnation
    if (p_rsolver%depsStag > 0.0_DP)&
        bstagnate = nlsol_checkStagnation(p_rsolver, 0)

    ! Perform relaxation of the time step
    call nlsol_relaxTimestep(rtimestep, p_rsolver)
    
    !----------------------------------------------------------------
    ! Iterative fixed-point correction
    !----------------------------------------------------------------
    
    fixpoint: do iiterations = 1, p_rsolver%nmaxIterations
            
      ! Store old defect values
      doldDefect = p_rsolver%dfinalDefect
      
      ! Compute solution increment
1     select case(p_rsolver%iprecond)
        
      case(NLSOL_PRECOND_BLOCKD)
        
        !---------------------------------------------------------------------
        ! Segregated/Block-diagonal approach
        !
        ! Apply block-diagonal preconditioning, that is, solve
        !   A[iblock,iblock]*aux[iblock]=res[iblock]
        ! for each iblock=1,nblock separately.
        !---------------------------------------------------------------------
        
        call lsysbl_clearVector(p_raux)
        do iblock = 1, rsolution%nblocks
          call linsol_solveMultigrid(rproblemLevel, p_rsolverLinear, &
                                     p_raux%RvectorBlock(iblock),&
                                     p_rres%RvectorBlock(iblock))
        end do
        
        ! Check status of linear solver: Early good return
        if ((p_rsolverLinear%istatus .eq. SV_ZERO_RHS) .or. &
            (p_rsolverLinear%istatus .eq. SV_ZERO_DEF)) exit fixpoint
        
        ! Add increment to solution vector
        if (p_rsolver%domega < 0.0_DP) then
          call lsysbl_vectorLinearComb(p_raux, rsolution, 1.0_DP, 1.0_DP)
        else
          call lsysbl_vectorLinearComb(p_raux, rsolution, p_rsolver%domega, 1.0_DP)
        end if
        
        ! Compute new defect
        call fcb_calcResidual(rproblemLevel, rtimestep, p_rsolver,&
                              rsolution, rsolutionInitial, p_rrhs, p_rres, iiterations, rcollection)
        
        ! Apply given right-hand side vector
        if (present(rrhs)) call lsysbl_vectorLinearComb(rrhs, p_rres, 1.0_DP, 1.0_DP)

        ! Impose boundary conditions
        call fcb_setBoundary(rproblemLevel, rtimestep, p_rsolver,&
                             rsolution, p_rres, rsolutionInitial, rcollection)
        if (rcollection%Iquickaccess(1) < 0) then
          p_rsolver%istatus = SV_INF_DEF
          return
        end if

        ! Compute norm of new defect
        p_rsolver%dfinalDefect = lsysbl_vectorNorm(p_rres, p_rsolver%iresNorm)

        if (p_rsolver%ioutputLevel .ge. SV_IOLEVEL_VERBOSE) then
          call output_lbrk()
          call output_separator(OU_SEP_TILDE)
          call output_line('Block-diagonal solution step '//trim(sys_siL(iiterations,5)))
          call output_line('Norm of residual:            '//trim(sys_sdEL(p_rsolver%dfinalDefect,5)))
          call output_line('Norm of previous residual:   '//trim(sys_sdEL(doldDefect,5)))
          call output_line('Variation of residual:       '//trim(sys_sdEL(&
                           p_rsolver%dfinalDefect/max(SYS_EPSREAL, doldDefect),5)))
          call output_line('Improvement of residual:     '//trim(sys_sdEL(&
                           p_rsolver%dfinalDefect/max(SYS_EPSREAL, p_rsolver%dinitialDefect),5)))
          call output_separator(OU_SEP_TILDE)
          call output_lbrk()
        end if
        
        ! Write nonlinear defect to logfile?
        if (p_rsolver%ioutputLevel .ge. SV_IOLEVEL_FILE) then
          write(UNIT=p_rsolver%iunitLogfile,FMT='("(01),",I10,3(",",E16.8E3))')&
              iiterations, p_rsolver%dfinalDefect,&
              p_rsolver%dfinalDefect/max(SYS_EPSREAL, doldDefect),&
              p_rsolver%dfinalDefect/max(SYS_EPSREAL, p_rsolver%dinitialDefect)
        end if
        

      case(NLSOL_PRECOND_DEFCOR,&
           NLSOL_PRECOND_NEWTON_FAILED)

        !---------------------------------------------------------------------
        ! Defect-correction
        !---------------------------------------------------------------------
        
        ! Call linear multigrid to solve A*aux=r
        call lsysbl_clearVector(p_raux)
        call linsol_solveMultigrid(rproblemLevel, p_rsolverLinear, p_raux, p_rres)
        
        ! Check status of linear solver: Early good return
        if ((p_rsolverLinear%istatus .eq. SV_ZERO_RHS) .or. &
            (p_rsolverLinear%istatus .eq. SV_ZERO_DEF)) exit fixpoint
                
        ! Add increment to solution vector
        if (p_rsolver%domega < 0.0_DP) then
          call lsysbl_vectorLinearComb(p_raux, rsolution, 1.0_DP, 1.0_DP)
        else
          call lsysbl_vectorLinearComb(p_raux, rsolution, p_rsolver%domega, 1.0_DP)
        end if
        
        ! Compute new defect
        call fcb_calcResidual(rproblemLevel, rtimestep, p_rsolver,&
                              rsolution, rsolutionInitial, p_rrhs, p_rres, iiterations, rcollection)

        ! Apply given right-hand side vector
        if (present(rrhs)) call lsysbl_vectorLinearComb(rrhs, p_rres, 1.0_DP, 1.0_DP)

        ! Impose boundary conditions
        call fcb_setBoundary(rproblemLevel, rtimestep, p_rsolver,&
                             rsolution, p_rres, rsolutionInitial, rcollection)
        if (rcollection%Iquickaccess(1) < 0) then
          p_rsolver%istatus = SV_INF_DEF
          return
        end if

        ! Compute norm of new defect
        p_rsolver%dfinalDefect = lsysbl_vectorNorm(p_rres, p_rsolver%iresNorm)

        if (p_rsolver%ioutputLevel .ge. SV_IOLEVEL_VERBOSE) then
          call output_lbrk()
          call output_separator(OU_SEP_TILDE)
          call output_line('Defect correction step       '//trim(sys_siL(iiterations,5)))
          call output_line('Norm of residual:            '//trim(sys_sdEL(p_rsolver%dfinalDefect,5)))
          call output_line('Norm of previous residual:   '//trim(sys_sdEL(doldDefect,5)))
          call output_line('Variation of residual:       '//trim(sys_sdEL(&
                           p_rsolver%dfinalDefect/max(SYS_EPSREAL, doldDefect),5)))
          call output_line('Improvement of residual:     '//trim(sys_sdEL(&
                           p_rsolver%dfinalDefect/max(SYS_EPSREAL, p_rsolver%dinitialDefect),5)))
          call output_separator(OU_SEP_TILDE)
          call output_lbrk()
        end if
        
        ! Write nonlinear defect to logfile?
        if (p_rsolver%ioutputLevel .ge. SV_IOLEVEL_FILE) then
          write(UNIT=p_rsolver%iunitLogfile,FMT='("(01),",I10,3(",",E16.8E3))')&
              iiterations, p_rsolver%dfinalDefect,&
              p_rsolver%dfinalDefect/max(SYS_EPSREAL, doldDefect),&
              p_rsolver%dfinalDefect/max(SYS_EPSREAL, p_rsolver%dinitialDefect)
        end if
                    

      case(NLSOL_PRECOND_NEWTON)

        !---------------------------------------------------------------------
        ! inexact Newton method
        !---------------------------------------------------------------------

        ! Compute forcing term (reset forcing term in the first iteration)
        call nlsol_calcForcingTerm(p_rsolver%dfinalDefect, doldDefect,&
                                   p_rsolverLinear%dfinalDefect,&
                                   drtlm, redfac, eta, iiterations,&
                                   p_rsolver%p_solverNewton%dforcingStrategy)
        
        ! Modify the solver structure for inexact Newton: The linear
        ! system needs to be solved such that the inexact Newton
        ! condition is satisfied. Hence, prescibe absolute tolerances.
        call solver_copySolver(p_rsolverLinear, rsolverTemp, .true., .false.)
        p_rsolverLinear%depsAbs = max(p_rsolverLinear%depsAbs, eta*p_rsolver%dfinalDefect)

        ! Call linear multigrid to solve A*aux=r
        call lsysbl_clearVector(p_raux)
        call linsol_solveMultigrid(rproblemLevel, p_rsolverLinear, p_raux, p_rres)

        ! Reset tolerance for linear solver
        call solver_copySolver(rsolverTemp, p_rsolverLinear, .true., .false.)

        ! Check status of linear solver: Early good return
        if ((p_rsolverLinear%istatus .eq. SV_ZERO_RHS) .or. &
            (p_rsolverLinear%istatus .eq. SV_ZERO_DEF)) exit fixpoint
        
        ! Compute the values G(U).TRANSP.*F(U)+J(U)*dU
        call lsysbl_copyVector(p_rres, p_rresfs)
        call fcb_applyJacobian(rproblemLevel, p_raux, p_rresfs, 1.0_DP, 1.0_DP, rcollection)
        drtlm = lsysbl_scalarProduct(p_rres, p_rresfs)

        ! Comput the value G(U).T.*J(U)*dU
        if (p_rsolver%iresNorm .eq. LINALG_NORMEUCLID) then
          drtjs = drtlm-p_rsolver%dfinalDefect**2
        else
          drtjs = drtlm-lsysbl_scalarProduct(p_rres, p_rres)
        end if

        ! Backtracking
        ibacktrackingSteps = nlsol_backtracking(rproblemLevel, rtimestep, p_rsolver,&
                                                rsolution, rsolutionInitial, p_rufs, p_raux,&
                                                p_rrhs, p_rres, p_rresfs, fcb_calcResidual,&
                                                fcb_setBoundary, doldDefect, drtjs, eta,&
                                                redfac, iiterations, rcollection, rrhs)
        
        ! Perform failsave defect correction if required
        if (ibacktrackingSteps .ne. 0) then
          
          call fcb_calcJacobian(rproblemLevel, rtimestep, p_rsolver,&
                                rsolution, rsolutionInitial, .true., rcollection)
          p_rsolver%iprecond = NLSOL_PRECOND_NEWTON_FAILED

          ! Impose boundary conditions
          call fcb_setBoundary(rproblemLevel, rtimestep, p_rsolver,&
                               rsolution, p_rres, rsolutionInitial, rcollection)
          if (rcollection%Iquickaccess(1) < 0) then
            p_rsolver%istatus = SV_INF_DEF
            return
          end if
          
          goto 1  ! ugly but quickest way

        else
          
          ! Assemble Jacobian matrix
          if (mod(iiterations-1, max(1, p_rsolver%p_solverNewton%iupdateFrequency)) .eq. 0)&
              call fcb_calcJacobian(rproblemLevel, rtimestep, p_rsolver, rsolution, rsolutionInitial, .false., rcollection)
          
          ! Impose boundary conditions
          call fcb_setBoundary(rproblemLevel, rtimestep, p_rsolver,&
                               rsolution, p_rres, rsolutionInitial, rcollection)
          if (rcollection%Iquickaccess(1) < 0) then
            p_rsolver%istatus = SV_INF_DEF
            return
          end if
          
        end if

        if (p_rsolver%ioutputLevel .ge. SV_IOLEVEL_VERBOSE) then
          call output_lbrk()
          call output_separator(OU_SEP_TILDE)
          call output_line('Inexact Newton step          '//trim(sys_siL(iiterations,5)))
          call output_line('Norm of residual:            '//trim(sys_sdEL(p_rsolver%dfinalDefect,5)))
          call output_line('Norm of previous residual:   '//trim(sys_sdEL(doldDefect,5)))
          call output_line('Variation of residual:       '//trim(sys_sdEL(&
                           p_rsolver%dfinalDefect/max(SYS_EPSREAL, doldDefect),5)))
          call output_line('Improvement of residual:     '//trim(sys_sdEL(&
                           p_rsolver%dfinalDefect/max(SYS_EPSREAL, p_rsolver%dinitialDefect),5)))
          call output_separator(OU_SEP_TILDE)
          call output_lbrk()
        end if
       
        ! Write nonlinear defect to logfile?
        if (p_rsolver%ioutputLevel .ge. SV_IOLEVEL_FILE) then
          write(UNIT=p_rsolver%iunitLogfile,FMT='("(01),",I10,3(",",E16.8E3))')&
              iiterations, p_rsolver%dfinalDefect,&
              p_rsolver%dfinalDefect/max(SYS_EPSREAL, doldDefect),&
              p_rsolver%dfinalDefect/max(SYS_EPSREAL, p_rsolver%dinitialDefect)
        end if
        
        
      case DEFAULT
        call output_line('Invalid nonlinear preconditioner!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'nlsol_solveFixedpointBlock')
        call sys_halt()
      end select
      
      ! Reset preconditioner from failed value (if required)
      p_rsolver%iprecond = abs(p_rsolver%iprecond)
      
      
      !--------------------------------------------------------------
      ! Convergence analysis
      !--------------------------------------------------------------
      
      ! Check if residual increased too much
      if (p_rsolver%dfinalDefect > p_rsolver%ddivRel*p_rsolver%dinitialDefect .or. &
          p_rsolver%dfinalDefect > p_rsolver%ddivAbs) then
        if (p_rsolver%ioutputLevel .ge. SV_IOLEVEL_WARNING) then
          call output_line('!!! Residual increased by factor '//trim(sys_sdEL(&
              p_rsolver%dfinalDefect/max(SYS_EPSREAL, p_rsolver%dinitialDefect),5))//' !!!')
        end if
        
        ! Adjust solver status
        p_rsolver%istatus = SV_INCR_DEF
        
        ! Also adjust solver status of top-most solver
        rsolver%istatus = SV_INCR_DEF
        
        ! That's it, return.
        return
      end if

      ! Check convergence criteria
      if (iiterations .ge. p_rsolver%nminIterations) then
        if ((p_rsolver%dfinalDefect .le. p_rsolver%depsAbs) .or.&
            (p_rsolver%dfinalDefect .le. p_rsolver%depsRel*p_rsolver%dinitialDefect)) then
          exit fixpoint
        elseif(p_rsolver%depsStag > 0.0_DP) then
          if (nlsol_checkStagnation(p_rsolver, iiterations)) exit fixpoint
        end if
      end if
      
    end do fixpoint
    
    ! Compute convergence rates and copy them to the top-most solver
    call solver_statistics(p_rsolver, iiterations)
    call solver_copySolver(p_rsolver, rsolver, .false., .true.)
    
    if (rsolver%ioutputLevel .ge. SV_IOLEVEL_INFO) then
      call output_lbrk()
      call output_separator(OU_SEP_HASH)
      call output_line('Nonlinear solution              '//solver_getstatus(rsolver))
      call output_line('Number of nonlinear iterations: '//trim(sys_siL(rsolver%iiterations,5)))
      call output_line('Norm of final residual:         '//trim(sys_sdEL(rsolver%dfinalDefect,5)))
      call output_line('Norm of initial residual:       '//trim(sys_sdEL(rsolver%dinitialDefect,5)))
      call output_line('Improvement of residual:        '//trim(sys_sdEL(&
                       rsolver%dfinalDefect/max(SYS_EPSREAL, rsolver%dinitialDefect),5)))
      call output_line('Convergence rate:               '//trim(sys_sdEL(rsolver%dconvergenceRate,5)))
      call output_separator(OU_SEP_HASH)
      call output_lbrk()
    end if

    ! Write nonlinear convergence rate to logfile?
    if (p_rsolver%ioutputLevel .ge. SV_IOLEVEL_FILE) then
      write(UNIT=p_rsolver%iunitLogfile,FMT='("(02),",I10,4(",",E16.8E3))')&
          rsolver%iiterations, rsolver%dfinalDefect, rsolver%dinitialDefect,&
          rsolver%dfinalDefect/max(SYS_EPSREAL, rsolver%dinitialDefect),&
          rsolver%dconvergenceRate
    end if
  end subroutine nlsol_solveFixedpointBlock

  !*****************************************************************************

!<function>

  function nlsol_backtracking(rproblemLevel, rtimestep, rsolver, ru, ru0,&
                              rufs, raux, rrhs, rres, rresfs, fcb_calcResidual,&
                              fcb_setBoundary, doldDefect, drtjs, eta, redfac,&
                              iiterations, rcollection, rf) result(istatus)

!<description>
    ! This function performs backtracking for the (inexact) Newton
    ! iteration and returns nonzero if backtracking failed
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(IN) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! norm of the old defect vector
    real(DP), intent(IN) :: doldDefect

    ! value of R(transposed)*(J*s)
    real(DP), intent(IN) :: drtjs

    ! iteration number
    integer, intent(IN) :: iiterations

    ! OPTIONAL: right-hand side vector
    type(t_vectorBlock), intent(IN), optional :: rf

    ! Callback routines
    include 'intf_nlsolcallback.inc'
!</input>

!<inputoutput>
    ! multigrid structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! nonlinear solver structure
    type(t_solver), intent(INOUT) :: rsolver
    
    ! last/new solution vector
    type(t_vectorBlock), intent(INOUT) :: ru
    
    ! failsave copy of last solution vector
    type(t_vectorBlock), intent(INOUT) :: ruFS

    ! constant right-hand side
    type(t_vectorBlock), intent(INOUT) :: rrhs
    
    ! new residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! failsave copy of last residual vector
    type(t_vectorBlock), intent(INOUT) :: rresfs
    
    ! (trial) inexact Newton step
    type(t_vectorBlock), intent(INOUT) :: raux

    ! (modified) forcing term
    real(DP), intent(INOUT) :: eta

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>

!<output>
    ! reduction factor
    real(DP), intent(OUT) :: redfac
!</output>

!<result>
    ! Status of backtracking
    integer :: istatus
!</result>
!</function>
    
    ! T - parameter for sufficient decrease condition. The default value is 1E-4.
    real(DP), parameter :: t = 1E-4_DP
    
    ! THETA_MIN - when backtracking occurs, this is the smallest reduction factor
    !             that will be applied to the current step in a single backtracking
    !             reduction. The default value is 0.1. Valid values are in the
    !             range [0.0, THETA_MAX] (see below).
    real(DP), parameter :: THETA_MIN = 0.1_DP
    
    ! THETA_MAX - when backtracking occurs, this is the largest reduction factor
    !             that will be applied to the current step in a single backtracking
    !             reduction. The default value is 0.5. Valid values are in the
    !             range [THETA_MIN,1.0) (see above).
    real(DP), parameter :: THETA_MAX = 0.5_DP

    
    ! local variables
    real(DP) :: theta
    integer :: ibacktrackingsteps


    ! Initialization
    redfac  = 1.0_DP
    istatus = 0
    call lsysbl_copyVector(ru,   ruFS)
    call lsysbl_copyVector(rres, rresFS)
    
    ! Backtracking loop
    backtrack: do ibacktrackingsteps = 1, rsolver%p_solverNewton%nmaxBacktrackingSteps
      
      ! Add increment to solution vector: u=u_k+aux
      if (rsolver%domega < 0.0_DP) then
        call lsysbl_vectorLinearComb(raux, ruFS, 1.0_DP, 1.0_DP, ru)
      else
        call lsysbl_vectorLinearComb(raux, ruFS, rsolver%domega, 1.0_DP, ru)
      end if
      
      ! Compute new provisional defect
      call fcb_calcResidual(rproblemLevel, rtimestep, rsolver,&
                            ru, ru0, rrhs, rres, iiterations, rcollection)
      
      ! Apply given right-hand side vector
      if (present(rf)) call lsysbl_vectorLinearComb(rf, rres, 1.0_DP, 1.0_DP)

      ! Impose boundary conditions
      call fcb_setBoundary(rproblemLevel, rtimestep, rsolver, ru, rres, ru0, rcollection)
      if (rcollection%Iquickaccess(1) < 0) exit backtrack

      ! Compute norm of nonlinear residual R(u_k+s_k)
      rsolver%dfinalDefect = lsysbl_vectorNorm(rres, rsolver%iresNorm)

      ! Return if sufficient decrease condition is satisfied or no checks need to be performed
      if ((rsolver%p_solverNewton%icheckSufficientDecrease .eq. 0) .or. &
          (rsolver%dfinalDefect .le. (1.0_DP-T*(1.0_DP-eta))*doldDefect)) then
        return
      end if
      
      ! Otherwise, compute a new step-length theta
      theta = rsolver%dfinalDefect*rsolver%dfinalDefect -&
              doldDefect*doldDefect - 2*drtjs*redfac
      if (abs(theta) > SYS_EPSREAL) then
        theta = -drtjs*redfac/theta
        theta = max(theta_min, min(theta, theta_max))
      else
        theta = theta_max
      end if
      
      ! Reduce inexact Newton step
      call lsysbl_scaleVector(raux, theta)
      eta    = 1.0_DP-theta*(1.0_DP-eta)
      redfac = theta*redfac
    end do backtrack

    ! If we end up here, then no acceptable Newton step was found
    call lsysbl_copyVector(ruFS, ru)
    call lsysbl_copyVector(rresFS, rres)
    istatus = -1
  end function nlsol_backtracking

  !*****************************************************************************

!<subroutine>
  
  subroutine nlsol_calcForcingTerm(dDefect, doldDefect, dlinearDefect, drtlm,&
                                   redfac, eta, iiterations, dforcingStrategy)

!<description>
    ! This subroutine computes the forcing term for the inexact Newton
    ! algorithm. The following stratagies are available:
    !
    ! (1) $(\|R(u_k)|-\|R(u_{k-1})+J(u_{k-1})*s_{k-1}\})/\|R(u_{k-1})\|$
    !
    ! (2) $\gamma*(\|R(u_k)\|/\|R(u_{k-1})\|)^\alpha$
    !     for user supplied $\gamma\in(0,1]$ and $\alpha\in(1,2]$
    !
    ! (3) $1/2^{k+1}$
    !
    ! (4) $\min(1/(k+2),\|R(u_k)\|)$
    !
    ! (otherwise) constant forcing term specified in the input file
    !
    ! The first two expressions are from S.C. Eisenstat and H.F.
    ! Walker, "Choosing the forcing term in an inexact Newton method",
    ! SIAM J. Scientific Computing, 17 (1996), pp. 16-32. The first
    ! gives convergence that is q-superlinear and of r-order 
    ! (1+sqrt(5))/2. The second gives convergence that is of q-order
    ! $\alpha$ when $\gamma<1$ and, when $\gamma=1$, of r-order 
    ! $\alpha$ and q-order $p$ for every $p\in[1,\alpha)$. In case
    ! $\alpha=2$ and $\gamma=1$, the second choice yields convergence
    ! that is r-quadratic and of q-order $p$ for every $p\in[1,2)$
    ! The third gives q-linear convergence with asymptotic rate
    ! constant $\alpha$ in a certain norm; see R.S. Dembo, S.C.
    ! Eisenstat, amd T. Steihaug, "Inexact Newton methods", SIAM J.
    ! Numer. Anal., 18 (1982), pp.400-408.
!</description>

!<input>
    ! current nonlinear residual R(u_k)
    real(DP), intent(IN) :: dDefect

    ! last nonlinear residual R(u_{k-1})
    real(DP), intent(IN) :: doldDefect

    ! current residual of the linear model
    !  R(u_{k-1})+J(u_{k-1})*s_{k-1}
    real(DP), intent(IN) :: dlinearDefect

    ! value of R(transposed)*(linear model)
    real(DP), intent(IN) :: drtlm

    ! net reduction factor of globalization
    real(DP), intent(IN) :: redfac

     ! number of current iteration
    integer, intent(IN) :: iiterations

    ! Strategy for choosing the forcing term
    real(DP), intent(IN) :: dforcingStrategy
!</input>

!<inputoutput>
    ! forcing term
    real(DP) :: eta
!</inputoutput>
!</subroutine>
    
    ! CHOICE1_EXP - parameter used in the update of the forcing term.
    !               This is the exponent for determining the ETA_MIN safeguard.
    !               The default value is CHOICE1_EXP = (1+SQRT(5))/2. A larger
    !               value woll allow ETA to decrease more rapidly, while a
    !               smaller value will result in a larger valuve for the safeguard.
    real(DP), parameter :: choice1_exp = 1.61803398874989_DP
    
    ! CHOICE2_EXP - parameter used in the update of the forcing term.
    !               This is the exponent ALPHA in the expression
    !               gamma*(||fcur||/||fprev||)**alpha; it is also used to
    !               determine the ETA_MIN safeguard. The default value is 2.0.
    !               Valid values are in the range (1.0, 2.0].
    real(DP), parameter :: choice2_exp = 2.0_DP
    
    ! CHOICE2_COEFF - parameter used in the update of the forcing term.
    !                 This is the coefficient GAMMA used in the expression
    !                 gamma*(||fcur||/||fprev||)**alpha; it is also used to
    !                 determine the ETA_MIN safeguard. The default value is 1.0.
    !                 Valid values are in the range (0.0, 1.0]
    real(DP), parameter :: choice2_coeff = 1.0_DP
    
    ! ETA_CUTOFF - parameter used to determine when to disable safeguarding the 
    !              update of the forcing term. The default value is 0.1. A value
    !              of 0.0 will enable safeguarding always; a value of 1.0 will 
    !              disable safeguarding always.
    real(DP), parameter :: eta_cutoff = 0.1_DP
    
    ! ETA_MAX - parameter used to provide an upper bound on the forcing term.
    !           This is necessary to ensure convergence of the inexact Newton
    !           iterates and is imposed whenever ETA would be too large otherwise.
    !           The default value of ETA_MAX is 1.0-1E-4. When backtracking 
    !           occurs several times during a nonlinear solution process, the
    !           forcing term can remain near ETA_MAX for several nonlinear steps
    !           and cause the nonlinear iterations to nearly stagnate. In such cases
    !           a smaller value of ETA_MAX may prevent this. Valid values are
    !           in the range (0.0, 1.0).
    real(DP), parameter :: eta_max = 1.0_DP-1E-4_DP
    
    ! ETA_MIN - parameter used to provide a lower bound on the forcing term.
    !           The default value is 0.1.
    real(DP), parameter :: eta_min = 0.1_DP
    
    ! ETA_0 - parameter used to initialize ETA in the first iteration. The default
    !         value is 0.5. If a used-supplied fixed ETA should be used, then
    !         this value is neglected
    real(DP), parameter :: eta_0 = 0.5_DP
    
    ! local variables
    real(DP) :: alpha,gamma,etamin
    real(DP) :: temp1,temp2
    
    
    ! Initialize quantities and return in the first iteration
    if (iiterations .eq. 1) then
      select case(int(dForcingStrategy))
      case (1, 2, 3, 4)
        eta = eta_0

      case DEFAULT
        eta = dforcingStrategy
      end select
      return
    end if

    
    select case(int(dForcingStrategy))
    case (1)
      ! Choice No. 1 by Eisenstat and Walker
      !
      ! eta_k=abs(||R(u_k)|| - ||R(u_{k-1})+J(u_{k-1})*s_{k-1}||) / ||R(u_{k-1})||
      !
      ! Safeguard:
      ! eta_k <- MAX(eta_k,eta_{k-1}^(1+sqrt(5))/2 whenever eta_{k-1}^(1+sqrt(5))/2 > 0.1
      alpha  = choice1_exp
      etamin = eta**alpha
      if (etamin .le. ETA_CUTOFF) etamin = 0.0_DP

      temp1  = 1.0_DP-redfac
      temp2  = sqrt( (temp1*doldDefect)**2+2*redfac*temp1*drtlm+(redfac*dlinearDefect)**2 )
      eta    = abs(dDefect-temp2)/doldDefect

      
    case (2)
      ! Choice No. 2 by Eisenstat and Walker
      !
      ! eta=GAMMA*( ||R(u_k)||/||R(u_{k-1})|| )^ALPHA
      !
      ! Safeguard:
      ! eta_k <- MAX(eta_k,GAMMA*eta_{k-1}^ALPHA) whenever GAMMA*eta_{k-1}^ALPHA > 0.1
      alpha  = choice2_exp
      gamma  = choice2_coeff
      etamin = gamma*eta*alpha
      if (etamin .le. ETA_CUTOFF) etamin = 0.0_DP

      eta = gamma*(dDefect/doldDefect)**alpha

    case (3)
      ! Decreasing forcing term by Brown and Saad which results in
      ! local q-superlinear convergence
      !
      ! eta_k=1/2^(k+1)
      etamin = ETA_MIN
      eta=1.0_DP/real(2**(min(30, iiterations+1)),DP)


    case (4)
      ! Choice by Dembo and Steihaug which results in local 
      ! q-quadratic convergence
      !
      ! eta_k=MIN(1/(k+2), ||R(u_k)||)
      etamin = ETA_MIN
      eta=min(1.0_DP/real(iiterations+2, DP), dDefect)


    case DEFAULT
      ! Constant forcing term specified in the input file which
      ! results in local q-linear convergence
      etamin = SYS_EPSREAL
      eta    = dforcingStrategy
    end select
    
    ! Apply global safeguards and store values for next iteration
    eta = max(etamin, min(eta, ETA_MAX))
  end subroutine nlsol_calcForcingTerm

  !*****************************************************************************

!<function>
  
  function nlsol_checkStagnation(rsolver, iiterations) result (stagnate)

!<description>
    ! This function checks if the nonlinear residual stagnates. Due
    ! to the application of algebraic flux limiting, the residual
    ! oscillates from one iteration to another so that the ratio of
    ! two consecutive residuals will never approach unity. To
    ! overcome this problem, stagnation of the nonlinear residual is
    ! monitored as follows: First, a linear approximation to the last
    ! couple of residuals is obtained by means of least sqaures fit.
    ! Then the coefficient of the linear term yields the slope of the
    ! approximate residual evolution. If the slope tends to zero,
    ! then the nonlinear iteration is considered to stagnate.
!</description>

!<input>
    ! current iteration number
    integer, intent(IN) :: iiterations
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver
!</inputoutput>

!<result>
    ! indicates if the residual stagnates
    logical :: stagnate
!</result>
!</function>
    
    ! local parameters
    real(DP), dimension(2,3) :: X3=reshape((/1D0,1D0,1D0,2D0,1D0,3D0/),(/2,3/))

    ! local variables
    real(DP), dimension(3), save :: Ddef
    real(DP), dimension(2) :: Dy
    real(DP) :: c1

    ! Initialization
    stagnate=.false.

    select case(iiterations)
    case (0)
      Ddef(1) = log10(rsolver%dfinalDefect)


    case (1)
      Ddef(2) = log10(rsolver%dfinalDefect)


    case DEFAULT
      Ddef(3) = log10(rsolver%dfinalDefect)
      Dy      = matmul(X3,Ddef)
      Ddef    = cshift(Ddef,1)
      c1      = -Dy(1)+0.5_DP*Dy(2)
      
      if (abs(c1) < rsolver%depsStag) then
        stagnate        = .true.
        rsolver%istatus = SV_STAGNATED
      end if
    end select
  end function nlsol_checkStagnation

  ! *****************************************************************************

!<subroutine>

  subroutine nlsol_relaxTimestep(rtimestep, rsolver)

!<description>
    ! This subroutine relaxes the time step subject to the given solver.
    ! This routine is not to be confused with the checking routine for
    ! the time step. This routine is used in the nonlinear solver, whereby
    ! the checking and possible adaption of the time step is done before.
!</description>

!<input>
    ! solver
    type(t_solver), intent(IN) :: rsolver
!</input>

!<inputoutput>
    ! time-stepping algorithm
    type(t_timestep), intent(INOUT) :: rtimestep
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_serController), pointer :: p_rserController
    real(DP) :: dStepOpt

    select case(rtimestep%iadaptTimestep)

    case (SV_TIMESTEP_SERADAPT)
      !-------------------------------------------------------------------------
      ! Switched evolution relaxation (SER)
      !
      ! From: W. Mulder and B. Van Leer. "Experiments with typical upwind 
      ! methods for the Euler equations", J. Comp. Phys., 59(1985), 232-246.
      !
      ! The new pseudo time step is computed from the following formula
      !
      !    Dt_opt = Dt * !! Defect(t^{n-1}) !! / !! Defect(t^n) !!
      !
      !-------------------------------------------------------------------------

      ! Set pointer
      p_rserController => rtimestep%p_rserController

      ! Store norm of stationary defect
      p_rserController%dsteadyDefect = rsolver%dinitialDefect

      ! Check if information about the stationary residual
      ! from the previous time step is available
      if (p_rserController%dsteadyDefect1 < SYS_EPSREAL) then

        ! Then adopt time step from previous iteration
        dStepOpt = rtimestep%dStep
      else
        
        ! Compute new time step size
        dStepOpt = rtimestep%dStep * p_rserController%dsteadyDefect1/&
                                     p_rserController%dsteadyDefect

        ! Limit the growth and reduction of the time step
        dStepOpt = max(p_rserController%dDecreaseFactor*rtimestep%dStep,&
                       min(p_rserController%dIncreaseFactor*rtimestep%dStep,&
                           dStepOpt))

        ! Impose upper/lower bounds on absolute values
        dStepOpt = max(rtimestep%dminStep, min(rtimestep%dmaxStep, dStepOpt))

        ! Adjust the simulation time accordingly
        rtimestep%dTime = rtimestep%dTime - rtimestep%dStep + dStepOpt
        rtimestep%dStep = dStepOpt
      end if
      
    end select
  end subroutine nlsol_relaxTimestep
  
end module solvernonlinear

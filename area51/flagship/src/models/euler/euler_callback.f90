!##############################################################################
!# ****************************************************************************
!# <name> euler_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to solve the
!# compressible Euler/Navier Stokes equations in arbitrary spatial dimensions.
!#
!# The following callback functions are available:
!#
!# 1.) euler_nlsolverCallback
!#     -> Callback routine for the nonlinear solver
!#
!# ****************************************************************************
!#
!# The following auxiliary routines are available:
!#
!# 1.) euler_calcPrecondThetaScheme
!#     -> Calculates the nonlinear preconditioner
!#        used in the two-level theta-scheme
!#
!# 2.) euler_calcJacobianThetaScheme
!#     -> Calculates the Jacobian matrix
!#        used in the two-level theta-scheme
!#
!# 3.) euler_calcResidualThetaScheme
!#     -> Calculates the nonlinear residual vector
!#        used in the two-level theta-scheme
!#
!# 4.) euler_calcRhsThetaScheme
!#     -> Calculates the explicit right-hand side vector
!#        used in the two-level theta-scheme
!#
!# 5.) euler_calcRhsRungeKuttaScheme
!#     -> Calculates the right-hand side vector
!#        used in the explicit Runge-Kutta scheme
!#
!# 6.) euler_setBoundaryConditions
!#     -> Imposes boundary conditions for nonlinear solver
!#        by filtering the system matrix and the solution/residual
!#        vector explicitly (i.e. strong boundary conditions)
!#
!# 7.) euler_calcLinearisedFCT
!#     -> Calculates the linearised FCT correction
!#
!# 8.) euler_calcFluxFCT
!#     -> Calculates the raw antidiffusive fluxes for FCT algorithm
!#
!# 9.) euler_calcCorrectionFCT
!#     -> Calculates the contribution of the antidiffusive fluxes
!#        limited by the FCT algorithm and applies them to the residual
!#
!# Frequently asked questions?
!#
!# 1.) What is the magic behind subroutine 'transp_nlsolverCallback'?
!#
!#     -> This is the main callback routine which is called by the
!#        nonlinear solution algorithm. The specifier ioperationSpec
!#        specifies the type of operations that should be performed,
!#        e.g., assemble the preconditioner, the residual and impose
!#        boundary values. If you want to implement a special routine
!#        for assembling the preconditioner, than you have to call it
!#        from transp_nlsolverCallback instead if the standard routine
!#        transp_calcResidualThetaScheme. Note that changing the
!#        interface of this callback routine would require to update
!#        ALL models. Hence, the interface should only be changed if
!#        really necessary.
!#
!# </purpose>
!##############################################################################

module euler_callback

  use afcstabilisation
  use basicgeometry
  use boundaryfilter
  use collection
  use euler_basic
  use euler_callback1d
  use euler_callback2d
  use euler_callback3d
  use flagship_basic
  use fsystem
  use genoutput
  use groupfemsystem
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use solveraux
  use statistics
  use storage
  use timestepaux

  implicit none

  private
  public :: euler_calcJacobianThetaScheme
  public :: euler_calcPrecondThetaScheme
  public :: euler_calcResidualThetaScheme
  public :: euler_calcRhsRungeKuttaScheme
  public :: euler_calcRhsThetaScheme
  public :: euler_nlsolverCallback
  public :: euler_setBoundaryConditions
  public :: euler_calcLinearisedFCT
  public :: euler_calcFluxFCT
  public :: euler_calcCorrectionFCT

contains

  !*****************************************************************************

!<subroutine>

  subroutine euler_nlsolverCallback(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rrhs, rres, istep,&
      ioperationSpec, rcollection, istatus, rsource)

!<description>
    ! This subroutine is called by the nonlinear solver and it is responsible
    ! to assemble preconditioner, right-hand side vector, residual vector, etc.
!</description>

!<input>
    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolution0

    ! number of solver step
    integer, intent(in) :: istep

    ! specifier for operations
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: given source vector
    type(t_vectorBlock), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(inout) :: rtimestep

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! right-hand side vector
    type(t_vectorBlock), intent(inout) :: rrhs

    ! residual vector
    type(t_vectorBlock), intent(inout) :: rres

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! status flag
    integer, intent(out) :: istatus
!</output>
!</subroutine>


    ! Do we have to calculate the preconditioner?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) then

      ! Compute the preconditioner
      call euler_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rcollection)
    end if


    ! Do we have to calculate the constant right-hand side?
    ! --------------------------------------------------------------------------
    if ((iand(ioperationSpec, NLSOL_OPSPEC_CALCRHS)  .ne. 0)) then

      ! Compute the right-hand side
      call euler_calcRhsRungeKuttaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rrhs, istep,&
          rcollection)
    end if


    ! Do we have to calculate the residual?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      if (istep .eq. 0) then
        ! Compute the constant right-hand side
        call euler_calcRhsThetaScheme(rproblemLevel, rtimestep,&
            rsolver, rsolution0, rrhs, rcollection, rsource)
      end if

      ! Compute the residual
      call euler_calcResidualThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rrhs, rres, istep,&
          rcollection)
    end if


    ! Do we have to impose boundary conditions?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      ! Impose boundary conditions
      call euler_setBoundaryConditions(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rres, rcollection)
    end if


    ! Set status flag
    istatus = 0

  end subroutine euler_nlsolverCallback

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rcollection)

!<description>
    ! This subroutine calculates the nonlinear preconditioner and
    ! configures the linear solver structure accordingly. Depending on
    ! the nonlinear solver, the low-order evolution operator or the
    ! Jacobian operator is adopted as nonlinear preconditioner.
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    integer :: systemMatrix, lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ, inviscidAFC
    integer :: isystemCoupling, isystemPrecond, isystemFormat, imasstype, ivar, jvar

    ! Start time measurement for matrix evaluation
    p_rtimer => collct_getvalue_timer(rcollection, 'rtimerAssemblyMatrix')
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'systemmatrix', systemMatrix)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'coeffMatrix_CX', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'coeffMatrix_CY', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'coeffMatrix_CZ', coeffMatrix_CZ)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'inviscidAFC', inviscidAFC)

    !---------------------------------------------------------------------------
    ! Check if fully explicit time-stepping is used
    !---------------------------------------------------------------------------
    if (rtimestep%theta .le. SYS_EPSREAL) then

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'isystemformat', isystemFormat)
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'imasstype', imasstype)
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'lumpedmassmatrix', lumpedMassMatrix)
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'consistentmassmatrix', consistentMassMatrix)

      select case(isystemFormat)
      case (SYSTEM_INTERLEAVEFORMAT)

        select case(imasstype)
        case (MASS_LUMPED)
          call lsyssc_spreadMatrix(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(systemMatrix))
        case (MASS_CONSISTENT)
          call lsyssc_spreadMatrix(&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rproblemLevel%Rmatrix(systemMatrix))
        case DEFAULT
          call output_line('Empty system matrix is invalid!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
          call sys_halt()
        end select

      case (SYSTEM_BLOCKFORMAT)

        select case(imasstype)
        case (MASS_LUMPED)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, euler_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(&
                rproblemLevel%Rmatrix(lumpedMassMatrix), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                .false., .false., .true., .true.)
          end do

        case (MASS_CONSISTENT)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, euler_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(&
                rproblemLevel%Rmatrix(consistentMassMatrix), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                .false., .false., .true., .true.)
          end do

        case DEFAULT
          call output_line('Empty system matrix is invalid!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
          call sys_halt()
        end select

      case DEFAULT
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
        call sys_halt()
      end select

      ! That is it
      return
    end if

    !---------------------------------------------------------------------------
    ! Assemble divergence operator:  $ -\nabla\cdot\mathbf{F}(U) $
    !
    ! Remark: In future versions this routine will assemble both the primal
    !         and the dual preconditioner depending on the variables
    !         primaldual which has to be extracted from the collection.
    !---------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'isystemCoupling', isystemCoupling)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'isystemPrecond', isystemPrecond)

    ! What kind of coupling is applied?
    select case(isystemCoupling)

    case (SYSTEM_SEGREGATED)

      !-------------------------------------------------------------------------
      ! Assemble block-diagonal divergence operator
      !-------------------------------------------------------------------------

      ! What kind of preconditioner is applied?
      select case(isystemPrecond)

      case (DISSIPATION_ZERO)

        ! Assemble divergence operator without dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonalDiag1d,&
              euler_calcMatrixGalerkinDiag1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonalDiag2d,&
              euler_calcMatrixGalerkinDiag2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonalDiag3d,&
              euler_calcMatrixGalerkinDiag3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonalDiag1d,&
              euler_calcMatrixScalarDissDiag1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonalDiag2d,&
              euler_calcMatrixScalarDissDiag2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonalDiag3d,&
              euler_calcMatrixScalarDissDiag3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecond')
          call sys_halt()
        end select


      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonalDiag1d,&
              euler_calcMatrixTensorDissDiag1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonalDiag2d,&
              euler_calcMatrixTensorDissDiag2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonalDiag3d,&
              euler_calcMatrixTensorDissDiag3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with the Rusanov dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonalDiag1d,&
              euler_calcMatrixRusanovDiag1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonalDiag2d,&
              euler_calcMatrixRusanovDiag2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonalDiag3d,&
              euler_calcMatrixRusanovDiag3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
        call sys_halt()
      end select


    case (SYSTEM_ALLCOUPLED)

      !-------------------------------------------------------------------------
      ! Assemble full block transport operator
      !-------------------------------------------------------------------------

      ! What kind of preconditioner is applied?
      select case(isystemPrecond)

      case (DISSIPATION_ZERO)

        ! Assemble divergence operator without dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonal1d, euler_calcMatrixGalerkin1d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonal2d, euler_calcMatrixGalerkin2d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonal3d, euler_calcMatrixGalerkin3d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonal1d,&
              euler_calcMatrixScalarDiss1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonal2d,&
              euler_calcMatrixScalarDiss2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonal3d,&
              euler_calcMatrixScalarDiss3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonal1d,&
              euler_calcMatrixTensorDiss1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonal2d,&
              euler_calcMatrixTensorDiss2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonal3d,&
              euler_calcMatrixTensorDiss3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with the Rusanov dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonal1d, euler_calcMatrixRusanov1d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonal2d, euler_calcMatrixRusanov2d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcMatrixDiagonal3d, euler_calcMatrixRusanov3d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
        call sys_halt()
      end select


    case DEFAULT
      call output_line('Invalid type of flow coupling!',&
          OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Assemble the global system operator
    !---------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'isystemformat', isystemFormat)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'imasstype', imasstype)

    select case(isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)

      !-------------------------------------------------------------------------
      ! Assemble global system operator in interleave matrix format
      !-------------------------------------------------------------------------

      select case(imasstype)
      case (MASS_LUMPED)

        !-----------------------------------------------------------------------
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_L)-theta*dt*L $
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(1), 'lumpedmassmatrix', lumpedMassMatrix)

        call lsyssc_MatrixLinearComb(&
            rproblemLevel%Rmatrix(lumpedMassMatrix), 1.0_DP,&
            rproblemLevel%Rmatrix(systemMatrix),&
            rtimestep%theta*rtimestep%dStep,&
            rproblemLevel%Rmatrix(systemMatrix),&
            .false., .false., .true., .true.)

      case (MASS_CONSISTENT)

        !-----------------------------------------------------------------------
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_C)-theta*dt*L $
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(1), 'consistentmassmatrix', consistentMassMatrix)

        call lsyssc_MatrixLinearComb(&
            rproblemLevel%Rmatrix(consistentMassMatrix), 1.0_DP,&
            rproblemLevel%Rmatrix(systemMatrix),&
            rtimestep%theta*rtimestep%dStep,&
            rproblemLevel%Rmatrix(systemMatrix),&
            .false., .false., .true., .true.)

      case DEFAULT

        !-----------------------------------------------------------------------
        ! Use the global operator for steady-state flow
        !
        !   $ A = -L $
        !-----------------------------------------------------------------------

      end select


    case (SYSTEM_BLOCKFORMAT)

      !-------------------------------------------------------------------------
      ! Assemble global system operator in block matrix format
      !-------------------------------------------------------------------------

      select case(imasstype)
      case (MASS_LUMPED)

        !-----------------------------------------------------------------------
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_L)-theta*dt*L $
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(1), 'lumpedmassmatrix', lumpedMassMatrix)

        do ivar = 1, euler_getNVAR(rproblemLevel)
          do jvar = 1, euler_getNVAR(rproblemLevel)

            if (ivar .eq. jvar) then

              call lsyssc_MatrixLinearComb(&
                  rproblemLevel%Rmatrix(lumpedMassMatrix), 1.0_DP,&
                  rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                  rtimestep%theta*rtimestep%dStep,&
                  rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                  .false., .false., .true., .true.)

            elseif (isystemCoupling .eq. SYSTEM_ALLCOUPLED) then

              call lsyssc_scaleMatrix(&
                  rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,jvar),&
                  rtimestep%theta*rtimestep%dStep)

            end if

          end do ! jvar
        end do ! ivar


      case (MASS_CONSISTENT)

        !-----------------------------------------------------------------------
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_C)-theta*dt*L $
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(1), 'consistentmassmatrix', consistentMassMatrix)

        do ivar = 1, euler_getNVAR(rproblemLevel)
          do jvar = 1, euler_getNVAR(rproblemLevel)

            if (ivar .eq. jvar) then

              call lsyssc_MatrixLinearComb(&
                  rproblemLevel%Rmatrix(consistentMassMatrix), 1.0_DP,&
                  rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                  rtimestep%theta*rtimestep%dStep,&
                  rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                  .false., .false., .true.,.true.)

            elseif (isystemCoupling .eq. SYSTEM_ALLCOUPLED) then

              call lsyssc_scaleMatrix(&
                  rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,jvar),&
                  rtimestep%theta*rtimestep%dStep)

            end if

          end do ! jvar
        end do ! ivar


      case DEFAULT

        !-----------------------------------------------------------------------
        ! Use the global operator for steady-state flow
        !
        !   $ A = -L $
        !-----------------------------------------------------------------------

      end select


    case DEFAULT
      call output_line('Invalid system format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrecondThetaScheme')
      call sys_halt()
    end select


    ! Impose boundary conditions
    call bdrf_filterMatrix(rsolver%rboundaryCondition, &
        rproblemLevel%RmatrixBlock(systemMatrix), 1.0_DP)

    ! Ok, we updated the (nonlinear) system operator successfully. Now we still
    ! have to link it to the solver hierarchy. This is done recursively.
    call flagship_updateSolverMatrix(rproblemLevel, rsolver,&
        systemMatrix, isystemFormat, UPDMAT_ALL,&
        rproblemLevel%ilev, rproblemLevel%ilev)

    ! Finally, we have to update the content of the solver hierarchy
    call solver_updateContent(rsolver)

    ! Stop time measurement for global operator
    call stat_stopTimer(p_rtimer)

  end subroutine euler_calcPrecondThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcJacobianThetaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rcollection)

!<description>
    ! This callback subroutine computes the Jacobian matrix for the
    !  compressible Euler/Navier-Stokes equations
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolution0
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout), target :: rproblemLevel

    ! nonlinear solver structure
    type(t_solver), intent(inout) :: rsolver

    ! current solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    print *, "!!! The calculation of the Jacobian matrix for the compressible !!!"
    print *, "!!! Euler/Navier Stokes equations has yet not been implemented  !!!"
    stop

  end subroutine euler_calcJacobianThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcRhsThetaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rrhs, rcollection, rsource)

!<description>
    ! This subroutine computes the constant right-hand side
    !
    !  $$ rhs = [M + (1-\theta)\Delta t K^n]U^n + S^n$$
    !
    ! where the (scaled) source term is optional.
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! OPTIONAL: source vector
    type(t_vectorBlock), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! right-hand side vector
    type(t_vectorBlock), intent(inout) :: rrhs

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    real(DP) :: dscale
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: consistentMassMatrix, lumpedMassMatrix, massMatrix
    integer :: inviscidAFC, imasstype, idissipationtype
    integer :: iblock


    ! Start time measurement for residual/rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection, 'rtimerAssemblyVector')
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'coeffmatrix_cx', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'coeffmatrix_cy', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'coeffmatrix_cz', coeffMatrix_CZ)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'imasstype', imasstype)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'idissipationtype', idissipationtype)

    ! Do we have some kind of mass matrix?
    select case(imasstype)
    case (MASS_LUMPED, MASS_CONSISTENT)

      ! Do we have some explicit part?
      if (rtimestep%theta .lt. 1.0_DP) then

        ! Compute scaling parameter
        dscale = (1.0_DP-rtimestep%theta) * rtimestep%dStep

        ! What type if stabilization is applied?
        select case(rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation)

        case (AFCSTAB_GALERKIN)

          !---------------------------------------------------------------------
          ! Compute the initial high-order right-hand side
          !
          !   $$ rhs = (1-theta)*dt*K(U^n)*U^n $$
          !---------------------------------------------------------------------

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildDivVector(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                euler_calcFluxGalerkin1d, dscale, .true., rrhs)

          case (NDIM2D)
            call gfsys_buildDivVector(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                euler_calcFluxGalerkin2d, dscale, .true., rrhs)

          case (NDIM3D)
            call gfsys_buildDivVector(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                euler_calcFluxGalerkin3d, dscale, .true., rrhs)
          end select

        case (AFCSTAB_UPWIND,&
              AFCSTAB_FEMFCT_CLASSICAL,&
              AFCSTAB_FEMFCT_ITERATIVE,&
              AFCSTAB_FEMFCT_IMPLICIT,&
              AFCSTAB_FEMFCT_LINEARISED)

          !---------------------------------------------------------------------
          ! Compute the initial low-order right-hand side
          !
          !   $$ rhs = (1-theta)*dt*L(U^n)*U^n $$
          !---------------------------------------------------------------------

          ! What type of dissipation is applied?
          select case(idissipationtype)

          case (DISSIPATION_ZERO)

            ! Assemble divergence operator without dissipation

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxGalerkin1d, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxGalerkin2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxGalerkin3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_SCALAR)

            ! Assemble divergence operator with scalar dissipation

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxScalarDiss1d, dscale, .true. , rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxScalarDiss2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxScalarDiss3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_SCALAR_DSPLIT)

            ! Assemble divergence operator with scalar dissipation
            ! adopting dimensional splitting

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxScalarDiss1d, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxScalarDissDiSp2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxScalarDissDiSp3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_TENSOR)

            ! Assemble divergence operator with tensorial dissipation

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxTensorDiss1d, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxTensorDiss2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxTensorDiss3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_TENSOR_DSPLIT)

            ! Assemble divergence operator with tensorial dissipation
            ! adopting dimensional splitting

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxTensorDiss1d, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxTensorDissDiSp2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxTensorDissDiSp3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_RUSANOV)

            ! Assemble divergence operator with Rusanov flux

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxRusanov1d, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxRusanov2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxRusanov3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_RUSANOV_DSPLIT)

            ! Assemble divergence operator with Rusanov flux
            ! adopting dimensional splitting

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxRusanov1d, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxRusanovDiSp2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxRusanovDiSp3d, dscale, .true., rrhs)
            end select

          case DEFAULT
            call output_line('Invalid type of dissipation!',&
                OU_CLASS_ERROR,OU_MODE_STD,'euler_calcRhsThetaScheme')
            call sys_halt()
          end select

        case (AFCSTAB_FEMTVD)

          !---------------------------------------------------------------------
          ! Compute the initial low-order right-hand side + FEM-TVD stabilisation
          !
          !   $$ rhs = (1-theta)dt*L(U^n)*U^n + F(U^n) $$
          !---------------------------------------------------------------------

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildDivVectorTVD(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                euler_calcFluxGalerkinNoBdr1d,&
                euler_calcCharacteristics1d, dscale, .true., rrhs)

          case (NDIM2D)
            call gfsys_buildDivVectorTVD(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                euler_calcFluxGalerkinNoBdr2d,&
                euler_calcCharacteristics2d, dscale, .true., rrhs)

          case (NDIM3D)
            call gfsys_buildDivVectorTVD(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                euler_calcFluxGalerkinNoBdr3d,&
                euler_calcCharacteristics3d, dscale, .true., rrhs)
          end select

        case DEFAULT
          call output_line('Invalid type of stabilisation!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcRhsThetaScheme')
          call sys_halt()
        end select

        !-----------------------------------------------------------------------
        ! Compute the transient term
        !
        !   $$ rhs := M*U^n + rhs $$
        !-----------------------------------------------------------------------

        massMatrix = merge(lumpedMassMatrix,&
            consistentMassMatrix, imasstype .eq. MASS_LUMPED)

        do iblock = 1, rsolution%nblocks
          call lsyssc_scalarMatVec(&
              rproblemLevel%Rmatrix(massMatrix),&
              rsolution%RvectorBlock(iblock),&
              rrhs%RvectorBlock(iblock), 1.0_DP , 1.0_DP)
        end do

      else ! theta = 1

        !-----------------------------------------------------------------------
        ! Compute the transient term
        !
        !   $$ rhs = M*U^n $$
        !-----------------------------------------------------------------------

        massMatrix = merge(lumpedMassMatrix,&
            consistentMassMatrix, imasstype .eq. MASS_LUMPED)

        do iblock = 1, rsolution%nblocks
          call lsyssc_scalarMatVec(&
              rproblemLevel%Rmatrix(massMatrix),&
              rsolution%RvectorBlock(iblock),&
              rrhs%RvectorBlock(iblock), 1.0_DP , 0.0_DP)
        end do

      end if ! theta

    case DEFAULT

      !-------------------------------------------------------------------------
      ! Initialize the constant right-hand side by zeros
      !
      !   $$ rhs = 0 $$
      !-------------------------------------------------------------------------

      ! Clear vector
      call lsysbl_clearVector(rrhs)

    end select

    ! Apply the source vector to the right-hand side (if any)
    if (present(rsource))&
        call lsysbl_vectorLinearComb(rsource, rrhs, 1.0_DP, 1.0_DP)

    ! Stop time measurement for rhs evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine euler_calcRhsThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcResidualThetaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rrhs, rres, ite, rcollection,&
      rsource)

!<description>
    ! This subroutine computes the nonlinear residual vector
    !
    ! $$ res^{(m)} = rhs - [M-\theta\Delta t K^{(m)}]U^{(m)} - S^{(m)} $$
    !
    ! for the standard two-level theta-scheme, whereby the (scaled)
    ! source term $s^{(m)}$ is optional. The constant right-hand side
    !
    !  $$ rhs = [M + (1-\theta)\Delta t K^n]U^n + S^n$$
    !
    ! must be provided via the precomputed vector rrhs.
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolution0

    ! number of nonlinear iteration
    integer, intent(in) :: ite

    ! OPTIONAL:  source vector
    type(t_vectorBlock), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! right-hand side vector
    type(t_vectorBlock), intent(inout) :: rrhs

    ! residual vector
    type(t_vectorBlock), intent(inout) :: rres

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    type(t_vectorBlock), pointer :: p_rpredictor
    real(DP) :: dscale
    integer(I32) :: ioperationSpec
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: consistentMassMatrix, lumpedMassMatrix, massMatrix
    integer :: inviscidAFC, imasstype, idissipationtype
    integer :: iblock

    ! Start time measurement for residual/rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection, 'rtimerAssemblyVector')
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'coeffmatrix_cx', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'coeffmatrix_cy', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'coeffmatrix_cz', coeffMatrix_CZ)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'imasstype', imasstype)

    ! Do we have some kind of mass matrix?
    select case(imasstype)
    case (MASS_LUMPED, MASS_CONSISTENT)

      !-------------------------------------------------------------------------
      ! Initialize the residual for transient flows
      !
      !   $$ res = rhs-M*U^{(m)} $$
      !-------------------------------------------------------------------------

      ! Apply constant right-hand side
      call lsysbl_copyVector(rrhs, rres)

      massMatrix = merge(lumpedMassMatrix,&
          consistentMassMatrix, imasstype .eq. MASS_LUMPED)

      ! Apply mass matrix
      do iblock = 1, rsolution%nblocks
        call lsyssc_scalarMatVec(&
            rproblemLevel%Rmatrix(massMatrix),&
            rsolution%RvectorBlock(iblock),&
            rres%RvectorBlock(iblock) , -1._DP, 1.0_DP)
      end do

    case DEFAULT

      !-----------------------------------------------------------------------
      ! Initialize the residual for stationary flows zeros
      !
      !   $$ res = rhs $$
      !-----------------------------------------------------------------------

      ! Apply constant right-hand side
      call lsysbl_copyVector(rrhs, rres)

    end select

    !---------------------------------------------------------------------------
    ! Update the residual vector
    !---------------------------------------------------------------------------

    ! Compute scaling parameter
    dscale = rtimestep%theta*rtimestep%dStep

    ! What type if stabilization is applied?
    select case(rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation)

    case (AFCSTAB_GALERKIN)

      !-------------------------------------------------------------------------
      ! Compute the high-order residual
      !
      !   $$ res := res + dt*theta*K(U^{(m)})*U^(m) $$
      !-------------------------------------------------------------------------

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkin1d, dscale, .false., rres)

      case (NDIM2D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkin2d, dscale, .false., rres)

      case (NDIM3D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkin3d, dscale, .false., rres)
      end select

    case (AFCSTAB_UPWIND,&
          AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT,&
          AFCSTAB_FEMFCT_LINEARISED)

      !-------------------------------------------------------------------------
      ! Compute the low-order residual
      !
      !   $$ res := res + dt*theta*L(U^{(m)})*U^(m) $$
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'idissipationtype', idissipationtype)

      ! What type of dissipation is applied?
      select case(idissipationtype)

      case (DISSIPATION_ZERO)

        ! Assemble divergence operator without dissipation
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxGalerkin1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxGalerkin2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxGalerkin3d, dscale, .false., rres)
        end select

      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss3d, dscale, .false., rres)
        end select

      case (DISSIPATION_SCALAR_DSPLIT)

        ! Assemble divergence operator with scalar dissipation
        ! adopting dimensional splitting

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDissDiSp2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDissDiSp3d, dscale, .false., rres)
        end select

      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss3d, dscale, .false., rres)
        end select

      case (DISSIPATION_TENSOR_DSPLIT)

        ! Assemble divergence operator with tensorial dissipation
        ! adopting dimensional splitting

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDissDiSp2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDissDiSp3d, dscale, .false., rres)
        end select

      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with Rusanov flux

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov3d, dscale, .false., rres)
        end select

      case (DISSIPATION_RUSANOV_DSPLIT)

        ! Assemble divergence operator with Rusanov flux adopting
        ! dimensional splitting

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanovDiSp2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanovDiSp3d, dscale, .false., rres)
        end select

      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_calcResidualThetaScheme')
        call sys_halt()
      end select

    case (AFCSTAB_FEMTVD)

      !-------------------------------------------------------------------------
      ! Compute the low-order residual + FEM-TVD stabilization
      !
      !   $$ res = res + dt*theta*L(U^{(m)})*U^{(m)} + F(U^{(m)}) $$
      !-------------------------------------------------------------------------

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkinNoBdr1d,&
            euler_calcCharacteristics1d, dscale, .false., rres)

      case (NDIM2D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkinNoBdr2d,&
            euler_calcCharacteristics2d, dscale, .false., rres)

      case (NDIM3D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkinNoBdr3d,&
            euler_calcCharacteristics3d, dscale , .false., rres)
      end select

    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'euler_calcResidualThetaScheme')
      call sys_halt()
    end select


    !-------------------------------------------------------------------------
    ! Perform algebraic flux correction for the inviscid term
    !
    !   $$ res = res + f^*(u^(m),u^n) $$
    !-------------------------------------------------------------------------

    ! What type if stabilization is applied?
    select case(rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation)
    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT)


      ! Set pointer to predictor
      p_rpredictor => rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorPredictor

      ! Compute low-order predictor ...
      if (ite .eq. 0) then
        ! ... only in the zeroth iteration
        if (rtimestep%theta .ne. 1.0_DP) then
          call lsysbl_invertedDiagMatVec(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rrhs, 1.0_DP, p_rpredictor)
        else
          call lsysbl_copyVector(rsolution, p_rpredictor)
        end if
      elseif (rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation&
              .eq. AFCSTAB_FEMFCT_ITERATIVE) then
        ! ... in each iteration for iterative limiting
        call lsysbl_invertedDiagMatVec(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rrhs, 1.0_DP, p_rpredictor)
      end if

      ! Assemble the raw antidiffusive fluxes
      call euler_calcFluxFCT(rproblemLevel, rsolution, rsolution,&
          rtimestep%theta, rtimestep%dStep, 1.0_DP, (ite .eq. 0), rcollection)

      !-------------------------------------------------------------------------
      ! Set operation specifier
      !-------------------------------------------------------------------------

      if (ite .eq. 0) then
        ! Perform standard flux correction in zeroth iteration
        ioperationSpec = AFCSTAB_FCTALGO_STANDARD
      else
        select case(rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation)
        case (AFCSTAB_FEMFCT_CLASSICAL)
          ! Perform standard flux correction without recomputing bounds
          ioperationSpec = AFCSTAB_FCTALGO_STANDARD-&
                           AFCSTAB_FCTALGO_BOUNDS

        case (AFCSTAB_FEMFCT_IMPLICIT)
          ! Perform semi-implicit flux correction
          ioperationSpec = AFCSTAB_FCTALGO_INITALPHA+&
                           AFCSTAB_FCTALGO_LIMITEDGE+&
                           AFCSTAB_FCTALGO_CORRECT+&
                           AFCSTAB_FCTALGO_CONSTRAIN

        case (AFCSTAB_FEMFCT_ITERATIVE)
          ! Perform standard flux correction
          ioperationSpec = AFCSTAB_FCTALGO_STANDARD
        end select
      end if

      ! Apply FEM-FCT algorithm
      call euler_calcCorrectionFCT(rproblemLevel, p_rpredictor,&
          rtimestep%dStep, .false., ioperationSpec, rres, rcollection)

      ! Subtract corrected antidiffusion from right-hand side
      if (rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation&
          .eq. AFCSTAB_FEMFCT_ITERATIVE) then
        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            p_rpredictor, rtimestep%dStep, .false.,&
            AFCSTAB_FCTALGO_CORRECT, rrhs)
      end if
    end select

    ! Apply the source vector to the residual  (if any)
    if (present(rsource)) then
      call lsysbl_vectorLinearComb(rsource, rres, -1.0_DP, 1.0_DP)
    end if

    ! Stop time measurement for residual/rhs evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine euler_calcResidualThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcRhsRungeKuttaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rrhs, istep, rcollection)

!<description>
    ! This subroutine computes the right-hand side vector
    ! used in the explicit Lax-Wendroff time-stepping scheme
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolution0

    ! number of explicit step
    integer, intent(in) :: istep
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! right-hand side vector
    type(t_vectorBlock), intent(inout) :: rrhs

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    real(DP) :: dscale
    integer :: lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: imasstype, inviscidAFC, idissipationtype
    integer :: iblock


    print *, "WARNING: This subroutine has not been tested!"
    stop

    ! Start time measurement for residual/rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection, 'rtimerAssemblyVector')
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'coeffMatrix_CX', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'coeffMatrix_CY', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'coeffMatrix_CZ', coeffMatrix_CZ)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'inviscidAFC', inviscidAFC)

    !---------------------------------------------------------------------------
    ! Initialize the right-hand side vector
    !---------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'imasstype', imasstype)

    select case(imasstype)
    case (MASS_LUMPED)

      !-------------------------------------------------------------------------
      ! Initialize the right-hand side vector
      !
      !  $ M_L*U $
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'lumpedmassmatrix', lumpedMassMatrix)

      do iblock = 1, rsolution%nblocks
        call lsyssc_scalarMatVec(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rsolution%RvectorBlock(iblock),&
            rrhs%RvectorBlock(iblock), 1.0_DP, 0.0_DP)
      end do

    case(MASS_CONSISTENT)

      !-------------------------------------------------------------------------
      ! Initialize the right-hand side vector
      !
      !  $ M_C*U $
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'consistentmassmatrix',&
          consistentMassMatrix)

      do iblock = 1, rsolution%nblocks
        call lsyssc_scalarMatVec(&
            rproblemLevel%Rmatrix(consistentMassMatrix),&
            rsolution%RvectorBlock(iblock),&
            rrhs%RvectorBlock(iblock), 1.0_DP, 0.0_DP)
      end do

    case DEFAULT

      ! Initialize the right-hand side vector by zeros
      call lsysbl_clearVector(rrhs)
    end select


    !---------------------------------------------------------------------------
    ! Compute the divergence term of the right-hand side
    !---------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'inviscidAFC', inviscidAFC)

    !---------------------------------------------------------------------------
    ! Compute the scaling parameter
    !
    !   $ weight*(1-theta)*dt $
    !---------------------------------------------------------------------------

    dscale = rtimestep%DmultistepWeights(istep)*(1.0_DP-rtimestep%theta)*rtimestep%dStep

    select case(rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation)

    case (AFCSTAB_GALERKIN)

      !-------------------------------------------------------------------------
      ! Compute the high-order right-hand side
      !
      !   $$ rhs = M*U+weight*(1-theta)*dt*K(U)*U $$
      !-------------------------------------------------------------------------

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkin1d, dscale, .false., rrhs)

      case (NDIM2D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkin2d, dscale, .false., rrhs)

      case (NDIM3D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkin3d, dscale, .false., rrhs)
      end select

    case (AFCSTAB_UPWIND,&
          AFCSTAB_FEMFCT_LINEARISED)

      !-----------------------------------------------------------------------
      ! Compute the low-order right-hand side
      !
      !   $$ rhs = M*U+weight*(1-theta)*dt*L(U)*U $$
      !-----------------------------------------------------------------------


      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'idissipationtype', idissipationtype)

      select case(idissipationtype)

      case (DISSIPATION_ZERO)

        ! Assemble divergence operator without dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxGalerkin1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxGalerkin2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxGalerkin3d, dscale, .false., rrhs)
        end select

      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss3d, dscale, .false., rrhs)

        end select

      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss3d, dscale, .false., rrhs)

        end select

      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with Rusanov flux

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov3d, dscale, .false., rrhs)

        end select

      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_calcRhsRungeKuttaScheme')
        call sys_halt()
      end select

    case (AFCSTAB_FEMTVD)

      !-----------------------------------------------------------------------
      ! Compute the low-order right-hand side + FEM-TVD stabilization
      !
      !   $$ rhs = M*U+weight*(1-theta)*dt*L(U)*U + F(U) $$
      !-----------------------------------------------------------------------

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkinNoBdr1d,&
            euler_calcCharacteristics1d, dscale, .false., rrhs)

      case (NDIM2D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkinNoBdr2d,&
            euler_calcCharacteristics2d, dscale, .false., rrhs)

      case (NDIM3D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkinNoBdr3d,&
            euler_calcCharacteristics3d, dscale, .false., rrhs)

      end select

    case DEFAULT
      call output_line('Invalid type of stabilization!',&
          OU_CLASS_ERROR,OU_MODE_STD,'euler_calcRhsRungeKuttaScheme')
      call sys_halt()
    end select

    ! Stop time measurement for global operator
    call stat_stopTimer(p_rtimer)

  end subroutine euler_calcRhsRungeKuttaScheme

  !*****************************************************************************

!<subroutine>

  subroutine euler_setBoundaryConditions(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rres, rcollection)

!<description>
    ! This subroutine imposes the nonlinear boundary conditions.
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! nonlinear solver structure
    type(t_solver), intent(in) :: rsolver

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolution0
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! residual vector
    type(t_vectorBlock), intent(inout) :: rres

    ! collection structure to provide additional
    ! information to the boundary setting routine
    type(t_collection), intent(InOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    integer :: imatrix, istatus


    ! Get parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

    ! What type of nonlinear preconditioner are we?
    select case(rsolver%iprecond)
    case (NLSOL_PRECOND_BLOCKD,&
          NLSOL_PRECOND_DEFCOR,&
          NLSOL_PRECOND_NEWTON_FAILED)

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'systemmatrix', imatrix)

    case (NLSOL_PRECOND_NEWTON)

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'jacobianmatrix', imatrix)

    case DEFAULT
      call output_line('Invalid nonlinear preconditioner!',&
          OU_CLASS_ERROR,OU_MODE_STD,'euler_setBoundaryConditions')
      call sys_halt()
    end select


    ! Impose boundary conditions for the solution vector and impose
    ! zeros in the residual vector and the off-diagonal positions
    ! of the system matrix which is obtained from the collection

    select case(rproblemLevel%rtriangulation%ndim)
    case (NDIM1D)
      call bdrf_filterSolution(&
          rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixBlock(imatrix),&
          rsolution, rres, rsolution0, rtimestep%dTime,&
          euler_calcBoundaryvalues1d, istatus)

    case (NDIM2D)
      call bdrf_filterSolution(&
          rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixBlock(imatrix),&
          rsolution, rres, rsolution0, rtimestep%dTime,&
          euler_calcBoundaryvalues2d, istatus)

    case (NDIM3D)
      call bdrf_filterSolution(&
          rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixBlock(imatrix),&
          rsolution, rres, rsolution0, rtimestep%dTime,&
          euler_calcBoundaryvalues3d, istatus)

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'euler_setBoundaryConditions')
      call sys_halt()
    end select

  end subroutine euler_setBoundaryConditions

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcLinearisedFCT(rbdrCond, rproblemLevel,&
      rtimestep, rsolver, rsolution, rcollection, rsource)

!<description>
    ! This subroutine calculates the linearised FCT correction
!</description>

!<input>
    ! boundary condition structure
    type(t_boundaryCondition), intent(in) :: rbdrCond

    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! OPTIONAL: source vector
    type(t_vectorBlock), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_timestep) :: rtimestepAux
    type(t_vectorBlock), pointer :: p_rpredictor
    type(t_vectorBlock) :: rlowBound, rupBound
    type(t_vectorScalar) :: rbeta, rvectorScalar
    type(t_parlist), pointer :: p_rparlist
    real(DP), dimension(:), pointer :: p_ML, p_Dalpha, p_Dbeta, p_Dflux
    real(DP), dimension(:), pointer :: p_Ddata, p_DlowBound, p_DupBound
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    character(len=SYS_STRLEN) :: sfailsafevariable
    integer :: inviscidAFC, lumpedMassMatrix, isystemFormat
    integer :: ifailsafe, nfailsafe, ivariable, nvariable

    ! Get parameters from parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

    ! Get more parameter from parameter list
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'inviscidAFC', inviscidAFC)

    ! Do we have to apply linearised FEM-FCT?
    if (rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation&
        .ne. AFCSTAB_FEMFCT_LINEARISED) return

    ! Get more parameters from parameter list
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'lumpedmassmatrix', lumpedmassmatrix)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'nfailsafe', nfailsafe)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'isystemformat', isystemFormat)
    
    !---------------------------------------------------------------------------
    ! Prepare failsafe flux correction
    !---------------------------------------------------------------------------
    if (nfailsafe .gt. 0) then
 
      ! Get number of failsafe variables
      nvariable = max(1,&
          parlst_querysubstrings(p_rparlist,&
          rcollection%SquickAccess(1), 'sfailsafevariable'))
      
      ! Create block vectors for the upper and lower bounds
      call lsysbl_createVectorBlock(rlowBound,&
          rproblemLevel%Rafcstab(inviscidAFC)%NEQ, nvariable, .false.)
      call lsysbl_createVectorBlock(rupBound,&
          rproblemLevel%Rafcstab(inviscidAFC)%NEQ, nvariable, .false.)
      
      ! Initialise lower bounds from low-order solution
      do ivariable = 1, nvariable
        
        ! Get variable declaration string
        call parlst_getvalue_string(p_rparlist,&
            rcollection%SquickAccess(1), 'sfailsafevariable',&
            sfailsafevariable, isubstring=ivariable)
        
        ! Get variable data from low-order solution
        call euler_getVariable(rsolution, trim(sfailsafevariable),&
            rlowBound%RvectorBlock(ivariable))
      end do
      
      ! Copy lower bounds to upper bounds
      call lsysbl_copyVector(rlowBound, rupBound)
    end if
    
    !---------------------------------------------------------------------------
    ! Linearised FEM-FCT algorithm
    !---------------------------------------------------------------------------

    ! Initialize dummy timestep
    rtimestepAux%dStep = 1.0_DP
    rtimestepAux%theta = 0.0_DP

    ! Set pointer to predictor
    p_rpredictor => rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorPredictor
    
    ! Compute low-order "right-hand side" without theta parameter
    call euler_calcRhsThetaScheme(rproblemLevel, rtimestepAux,&
        rsolver, rsolution, p_rpredictor, rcollection, rsource)

    ! Compute low-order predictor
    call lsysbl_invertedDiagMatVec(&
        rproblemLevel%Rmatrix(lumpedMassMatrix),&
        p_rpredictor, 1.0_DP, p_rpredictor)

    ! Compute the raw antidiffusive fluxes
    call euler_calcFluxFCT(rproblemLevel, p_rpredictor,&
        rsolution, 0.0_DP, 1.0_DP, 1.0_DP, .true., rcollection)
    
    ! Apply FEM-FCT correction
    call euler_calcCorrectionFCT(rproblemLevel,&
        rsolution, rtimestep%dStep, .false.,&
        AFCSTAB_FCTALGO_STANDARD+&
        AFCSTAB_FCTALGO_SCALEBYMASS,&
        rsolution, rcollection)
    
    !---------------------------------------------------------------------------
    ! Perform failsafe flux correction
    !---------------------------------------------------------------------------
    if (nfailsafe .gt. 0) then
      
      ! Set pointers
      call afcstab_getbase_IverticesAtEdge(&
          rproblemLevel%Rafcstab(inviscidAFC), p_IverticesAtEdge)
      
      ! Compute upper and lower bounds from low-order solution
      do ivariable = 1, nvariable
        
        ! Make a copy of the scalar subvector
        call lsyssc_copyVector(rlowBound%RvectorBlock(ivariable), rvectorScalar)

        ! Set pointers
        call lsyssc_getbase_double(rlowBound%RvectorBlock(ivariable), p_DlowBound)
        call lsyssc_getbase_double(rupBound%RvectorBlock(ivariable), p_DupBound)
        call lsyssc_getbase_double(rvectorScalar, p_Ddata)
        
        ! Compute bounds for variable
        call computeBounds(p_IverticesAtEdge, p_Ddata, p_DlowBound, p_DupBound)
      end do
      
      ! Make a copy of the flux corrected predictor
      call lsysbl_copyVector(rsolution, p_rpredictor)
      
      ! Initialise the edgewise correction factors
      call lsyssc_createVector(rbeta,&
          rproblemLevel%Rafcstab(inviscidAFC)%NEDGE, .true.)
      
      ! Set pointers
      call lsyssc_getbase_double(rbeta, p_Dbeta)
      call lsyssc_getbase_double(&
          rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorAlpha, p_Dalpha)
      call lsyssc_getbase_double(&
          rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorFlux, p_Dflux)
      call lsyssc_getbase_double(&
          rproblemLevel%Rmatrix(lumpedMassMatrix), p_ML)

      ! Failsafe steps
      do ifailsafe = 1, nfailsafe

        ! Loop over failsafe variables
        do ivariable = 1, nvariable
          
          ! Get variable declaration string
          call parlst_getvalue_string(p_rparlist,&
              rcollection%SquickAccess(1), 'sfailsafevariable',&
              sfailsafevariable, isubstring=ivariable)
          
          ! Get variable data flux corrected solution
          call euler_getVariable(rsolution, trim(sfailsafevariable), rvectorScalar)
          
          ! Set pointers
          call lsyssc_getbase_double(rlowBound%RvectorBlock(ivariable), p_DlowBound)
          call lsyssc_getbase_double(rupBound%RvectorBlock(ivariable), p_DupBound)
          call lsyssc_getbase_double(rvectorScalar, p_Ddata)
          
          ! Compute failsafe correction factors
          call computeFailsafe(p_IverticesAtEdge,&
              rtimestep%dStep*real(ifailsafe, DP)/real(nfailsafe, DP),&
              p_Ddata, p_DlowBound, p_DupBound, 1e-8_DP, p_Dbeta)
        end do
        
        ! Restore flux correction solution
        call lsysbl_copyVector(p_rpredictor, rsolution)
        call lsysbl_getbase_double(rsolution, p_Ddata)
        
        ! Apply correction factors to solution vector
        if (isystemFormat .eq. SYSTEM_INTERLEAVEFORMAT) then
          call applyFailsafeInterleaveFormat(p_IverticesAtEdge,&
              rproblemLevel%Rafcstab(inviscidAFC)%NEDGE,&
              rproblemLevel%Rafcstab(inviscidAFC)%NEQ,&
              rproblemLevel%Rafcstab(inviscidAFC)%NVAR,&
              p_ML, p_Dalpha, p_Dbeta, p_Dflux, p_Ddata)
        else
          call applyFailsafeBlockFormat(p_IverticesAtEdge,&
              rproblemLevel%Rafcstab(inviscidAFC)%NEDGE,&
              rproblemLevel%Rafcstab(inviscidAFC)%NEQ,&
              rproblemLevel%Rafcstab(inviscidAFC)%NVAR,&
              p_ML, p_Dalpha, p_Dbeta, p_Dflux, p_Ddata)
        end if
      end do

      ! Release temporal vector
      call lsysbl_releaseVector(rlowBound)
      call lsysbl_releaseVector(rupBound)
      call lsyssc_releaseVector(rbeta)
      call lsyssc_releaseVector(rvectorScalar)
    end if
    

    ! Impose boundary conditions for the solution vector
    select case(rproblemLevel%rtriangulation%ndim)
    case (NDIM1D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dTime, euler_calcBoundaryvalues1d)

    case (NDIM2D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dTime, euler_calcBoundaryvalues2d)

    case (NDIM3D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dTime, euler_calcBoundaryvalues3d)
    end select

  contains
    
    ! Here are some internal working routines 
    
    !**************************************************************
    ! Compute the upper and lower bounds
    
    subroutine computeBounds(IverticesAtEdge,&
        Dx, DlowBound, DupBound)

      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      real(DP), dimension(:), intent(in) :: Dx
      
      real(DP), dimension(:), intent(inout) :: DlowBound, DupBound

      ! local variables
      integer :: iedge,i,j

      ! Loop over all edges
      do iedge = 1, size(IverticesAtEdge,2)

        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        
        DlowBound(i) = min(DlowBound(i), Dx(j))
        DlowBound(j) = min(DlowBound(j), Dx(i))
        DupBound(i)  = max(DupBound(i),  Dx(j))
        DupBound(j)  = max(DupBound(j),  Dx(i))
      end do

    end subroutine computeBounds

    !**************************************************************
    ! Compute the failsafe correction factors
    
    subroutine computeFailsafe(IverticesAtEdge,&
        dscale, Dx, DlowBound, DupBound, dtolerance, Dbeta)
      
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      real(DP), dimension(:), intent(in) :: Dx, Dlowbound, DupBound
      real(DP), intent(in) :: dscale, dtolerance
      
      real(DP), dimension(:), intent(inout) :: Dbeta
      
      ! local variables
      integer :: iedge,i,j
      
      ! Loop over all edges
      !$omp parallel do private(i,j)
      do iedge = 1, size(IverticesAtEdge,2)
        
        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        if ((Dx(i) .lt. DlowBound(i)-dtolerance) .or.&
            (Dx(j) .lt. DlowBound(j)-dtolerance) .or.&
            (Dx(i) .gt. DupBound(i)+dtolerance) .or.&
            (Dx(j) .gt. DupBound(j)+dtolerance)) Dbeta(iedge) = dscale
      end do
      !$omp end parallel do

    end subroutine computeFailsafe

    !**************************************************************
    ! Apply the failsafe correction factors
    ! The data vecto is stored in interleave format

    subroutine applyFailsafeInterleaveFormat(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, ML, Dalpha, Dbeta, Dflux, Dx)
      
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      real(DP), dimension(:), intent(in) :: ML, Dalpha, Dbeta
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      integer, intent(in) :: NEDGE, NEQ, NVAR
      
      real(DP), dimension(NVAR,NEQ), intent(inout) :: Dx
      
      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,i,j

      ! Loop over all edges
      do iedge = 1, size(IverticesAtEdge,2)
        
        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Compute portion of corrected antidiffusive flux
        F_ij = Dbeta(iedge) * Dalpha(iedge) * Dflux(:,iedge)

        ! Remove flux from solution
        Dx(:,i) = Dx(:,i) - F_ij/ML(i)
        Dx(:,j) = Dx(:,j) + F_ij/ML(j)
      end do
      
    end subroutine applyFailsafeInterleaveFormat

    !**************************************************************
    ! Apply the failsafe correction factors
    ! The data vecto is stored in interleave format

    subroutine applyFailsafeBlockFormat(IverticesAtEdge,&
        NEDGE, NEQ, NVAR, ML, Dalpha, Dbeta, Dflux, Dx)
      
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      real(DP), dimension(:), intent(in) :: ML, Dalpha, Dbeta
      real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux
      integer, intent(in) :: NEDGE, NEQ, NVAR
      
      real(DP), dimension(NEQ,NVAR), intent(inout) :: Dx
      
      ! local variables
      real(DP), dimension(NVAR) :: F_ij
      integer :: iedge,i,j

      ! Loop over all edges
      do iedge = 1, size(IverticesAtEdge,2)
        
        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Compute portion of corrected antidiffusive flux
        F_ij = Dbeta(iedge) * Dalpha(iedge) * Dflux(:,iedge)

        ! Remove flux from solution
        Dx(i,:) = Dx(i,:) - F_ij/ML(i)
        Dx(j,:) = Dx(j,:) + F_ij/ML(j)
      end do
      
    end subroutine applyFailsafeBlockFormat

  end subroutine euler_calcLinearisedFCT

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcFluxFCT(rproblemLevel, rsolution1, rsolution2,&
      theta, tstep, dscale, binit, rcollection)

!<description>
    ! This subroutine calculates the raw antidiffusive fluxes for
    ! the different FCT algorithms
!</description>

!<input>
    ! solution vectors
    type(t_vectorBlock), intent(in) :: rsolution1, rsolution2

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! Switch for flux assembly
    ! TRUE  : assemble the initial antidiffusive flux
    ! FALSE : assemble the antidiffusive flux using some initial values
    logical, intent(in) :: binit
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    integer :: inviscidAFC, lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: idissipationtype, imassantidiffusiontype

    ! Get parameters from parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'coeffMatrix_CX', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'coeffMatrix_CY', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'coeffMatrix_CZ', coeffMatrix_CZ)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'lumpedmassmatrix', lumpedmassmatrix)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'consistentmassmatrix', consistentmassmatrix)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'imassantidiffusiontype', imassantidiffusiontype)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'idissipationtype', idissipationtype)


    ! What type of dissipation is applied?
    select case(idissipationtype)

    case (DISSIPATION_SCALAR)

      ! Assemble raw antidiffusive fluxes using scalar dissipation

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTScalarDiss1d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTScalarDiss1d,&
              theta, tstep, dscale, binit)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTScalarDiss2d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTScalarDiss2d,&
              theta, tstep, dscale, binit)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTScalarDiss3d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTScalarDiss3d,&
              theta, tstep, dscale, binit)
        end if
      end select


    case (DISSIPATION_TENSOR)

      ! Assemble raw antidiffusive fluxes using tensorial dissipation

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTTensorDiss1d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTTensorDiss1d,&
              theta, tstep, dscale, binit)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTTensorDiss2d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTTensorDiss2d,&
              theta, tstep, dscale, binit)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTTensorDiss3d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTTensorDiss3d,&
              theta, tstep, dscale, binit)
        end if
      end select


    case (DISSIPATION_RUSANOV)

      ! Assemble raw antidiffusive fluxes using th Rusanov dissipation

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTRusanov1d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTRusanov1d,&
              theta, tstep, dscale, binit)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTRusanov2d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTRusanov2d,&
              theta, tstep, dscale, binit)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTRusanov3d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, euler_calcFluxFCTRusanov3d,&
              theta, tstep, dscale, binit)
        end if
      end select


    case DEFAULT
      call output_line('Invalid type of dissipation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'euler_calcFluxFCT')
      call sys_halt()
    end select

  end subroutine euler_calcFluxFCT

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcCorrectionFCT(rproblemLevel, rsolution, &
      dscale, bclear, ioperationSpec, rresidual, rcollection)

!<description>
    ! This subroutine calculates the raw antidiffusive fluxes for
    ! the different FCT algorithms
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_FCT_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! residual vector
    type(t_vectorBlock), intent(inout) :: rresidual

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    character(len=SYS_STRLEN) :: slimitingvariable
    integer(I32) :: iopSpec
    integer :: inviscidAFC, lumpedMassMatrix, ivariable, nvariable


    ! Get parameters from parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'lumpedmassmatrix', lumpedmassmatrix)

    ! Get number of limiting variables
    nvariable = max(1,&
        parlst_querysubstrings(p_rparlist,&
        rcollection%SquickAccess(1), 'slimitingvariable'))

    ! Copy operation specifier and disable the correction step
    ! if sequential/multiplicative flux correction is performed
    if (nvariable .gt. 1) then
      iopSpec = iand(ioperationSpec, not(AFCSTAB_FCTALGO_CORRECT))
    else
      iopSpec = ioperationSpec
    end if

    ! Loop over items in the list of variables that should
    ! be limited sequentially, i.e., in multiplicative way
    do ivariable = 1, nvariable

      ! Get variable declaration string
      call parlst_getvalue_string(p_rparlist,&
          rcollection%SquickAccess(1), 'slimitingvariable',&
          slimitingvariable, isubstring=ivariable)

      ! What type of flux transformation is applied?
      if (trim(slimitingvariable) .eq. 'density') then

        ! Apply FEM-FCT algorithm for density fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxDensity1d, euler_trafoDiffDensity1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxDensity2d, euler_trafoDiffDensity2d)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxDensity3d, euler_trafoDiffDensity3d)
        end select

      elseif (trim(slimitingvariable) .eq. 'energy') then

        ! Apply FEM-FCT algorithm for energy fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxEnergy1d, euler_trafoDiffEnergy1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxEnergy2d, euler_trafoDiffEnergy2d)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxEnergy3d, euler_trafoDiffEnergy3d)
        end select

      elseif (trim(slimitingvariable) .eq. 'pressure') then

        ! Apply FEM-FCT algorithm for pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxPressure1d, euler_trafoDiffPressure1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxPressure2d, euler_trafoDiffPressure2d)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxPressure3d, euler_trafoDiffPressure3d)
        end select

      elseif (trim(slimitingvariable) .eq. 'velocity_x') then

        ! Apply FEM-FCT algorithm for x-velocity fluxes
        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rsolution, dscale, bclear, iopSpec, rresidual,&
            euler_trafoFluxVelocity1d, euler_trafoDiffVelocity1d)
        
      elseif (trim(slimitingvariable) .eq. 'velocity_y') then

        ! Apply FEM-FCT algorithm for y-velocity fluxes
        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rsolution, dscale, bclear, iopSpec, rresidual,&
            euler_trafoFluxVelocity2d, euler_trafoDiffVelocity2d)

      elseif (trim(slimitingvariable) .eq. 'velocity_z') then

        ! Apply FEM-FCT algorithm for z-velocity fluxes
        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rsolution, dscale, bclear, iopSpec, rresidual,&
            euler_trafoFluxVelocity3d, euler_trafoDiffVelocity3d)
        
      elseif (trim(slimitingvariable) .eq. 'momentum_x') then

        ! Apply FEM-FCT algorithm for x-momentum fluxes
        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rsolution, dscale, bclear, iopSpec, rresidual,&
            euler_trafoFluxMomentum1d, euler_trafoDiffMomentum1d)

      elseif (trim(slimitingvariable) .eq. 'momentum_y') then

        ! Apply FEM-FCT algorithm for y-momentum fluxes
        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rsolution, dscale, bclear, iopSpec, rresidual,&
            euler_trafoFluxMomentum2d, euler_trafoDiffMomentum2d)

      elseif (trim(slimitingvariable) .eq. 'momentum_z') then

        ! Apply FEM-FCT algorithm for z-momentum fluxes
        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rsolution, dscale, bclear, iopSpec, rresidual,&
            euler_trafoFluxMomentum3d, euler_trafoDiffMomentum3d)
        
      elseif (trim(slimitingvariable) .eq. 'density,energy') then

        ! Apply FEM-FCT algorithm for density and energy fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxDenEng1d, euler_trafoDiffDenEng1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxDenEng2d, euler_trafoDiffDenEng2d)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxDenEng3d, euler_trafoDiffDenEng3d)
        end select

      elseif (trim(slimitingvariable) .eq. 'density,pressure') then

        ! Apply FEM-FCT algorithm for density and pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxDenPre1d, euler_trafoDiffDenPre1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxDenPre2d, euler_trafoDiffDenPre2d)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxDenPre3d, euler_trafoDiffDenPre3d)
        end select

      elseif (trim(slimitingvariable) .eq. 'density,energy,momentum') then

        ! Apply FEM-FCT algorithm for full conservative fluxes
        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rsolution, dscale, bclear, iopSpec, rresidual)

      elseif (trim(slimitingvariable) .eq. 'density,pressure,velocity') then

        ! Apply FEM-FCT algorithm for density, velocity and pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxDenPreVel1d, euler_trafoDiffDenPreVel1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxDenPreVel2d, euler_trafoDiffDenPreVel2d)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              euler_trafoFluxDenPreVel3d, euler_trafoDiffDenPreVel3d)
        end select

      elseif (trim(slimitingvariable) .eq. 'none') then
        
        ! Apply raw antidiffusive fluxes without correction
        iopSpec = ioperationSpec
        iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_ADINCREMENTS))
        iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_BOUNDS))
        iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_LIMITNODAL))
        iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_LIMITEDGE))
        
        ! Enforce existence of edgewise correction factors
        rproblemLevel%Rafcstab(inviscidAFC)%iSpec = ior(&
            rproblemLevel%Rafcstab(inviscidAFC)%iSpec, AFCSTAB_HAS_EDGELIMITER)

        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rsolution, dscale, bclear, iopSpec, rresidual)
        
        ! Nothing more needs to be done
        return

      else
        call output_line('Invalid type of flux transformation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_calcCorrectionFCT')
        call sys_halt()
      end if

      ! Disable the initialisation of edgewise correction factors
      ! in all but the first iteration over variables
      iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_INITALPHA))
    end do

    ! Perform the correction step separately (if required)
    if (nvariable .gt. 1) then
      
      ! Copy original specifier
      iopSpec = ioperationSpec

      ! Remove all tasks which might have been performed before
      iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_INITALPHA))
      iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_ADINCREMENTS))
      iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_BOUNDS))
      iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_LIMITNODAL))
      iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_LIMITEDGE))
      
      call gfsys_buildDivVectorFCT(&
          rproblemLevel%Rmatrix(lumpedMassMatrix),&
          rproblemLevel%Rafcstab(inviscidAFC),&
          rsolution, dscale, bclear, iopSpec, rresidual)
    end if

  end subroutine euler_calcCorrectionFCT

end module euler_callback

!##############################################################################
!# ****************************************************************************
!# <name> eulerlagrange_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to solve the
!# compressible Euler/Navier Stokes equations in arbitrary spatial dimensions.
!#
!# The following callback functions are available:
!#
!# 1.) eulerlagrange_nlsolverCallback
!#     -> Callback routine for the nonlinear solver
!#
!# ****************************************************************************
!#
!# The following auxiliary routines are available:
!#
!# 1.) eulerlagrange_calcPrecondThetaScheme
!#     -> Calculates the nonlinear preconditioner
!#        used in the two-level theta-scheme
!#
!# 2.) eulerlagrange_calcJacobianThetaScheme
!#     -> Calculates the Jacobian matrix
!#        used in the two-level theta-scheme
!#
!# 3.) eulerlagrange_calcResidualThetaScheme
!#     -> Calculates the nonlinear residual vector
!#        used in the two-level theta-scheme
!#
!# 4.) eulerlagrange_calcRhsThetaScheme
!#     -> Calculates the explicit right-hand side vector
!#        used in the two-level theta-scheme
!#
!# 5.) eulerlagrange_calcRhsRungeKuttaScheme
!#     -> Calculates the right-hand side vector
!#        used in the explicit Runge-Kutta scheme
!#
!# 6.) eulerlagrange_setBoundaryConditions
!#     -> Imposes boundary conditions for nonlinear solver
!#        by filtering the system matrix and the solution/residual
!#        vector explicitly (i.e. strong boundary conditions)
!#
!# 7.) eulerlagrange_calcLinearisedFCT
!#     -> Calculates the linearised FCT correction
!#
!# 8.) eulerlagrange_calcFluxFCT
!#     -> Calculates the raw antidiffusive fluxes for FCT algorithm
!#
!# 9.) eulerlagrange_calcCorrectionFCT
!#     -> Calculates the contribution of the antidiffusive fluxes
!#        limited by the FCT algorithm and applies them to the residual
!#
!# 10.) eulerlagrange_limitEdgewiseVelocity
!#      -> Performs synchronized flux correction for the velocity
!#
!# 11.) eulerlagrange_limitEdgewiseMomentum
!#      -> Performs synchronized flux correction for the momentum
!#
!# 12.) eulerlagrange_coeffVectorFE
!#      -> Callback routine for the evaluation of linear forms
!#         using a given FE-solution for interpolation
!#
!# 13.) eulerlagrange_coeffVectorAnalytic
!#      -> Callback routine for the evaluation of linear forms
!#         using a given FE-solution for interpolation

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

module eulerlagrange_callback

  use afcstabilisation
  use basicgeometry
  use boundaryfilter
  use collection
  use eulerlagrange_basic
  use eulerlagrange_callback1d
  use eulerlagrange_callback2d
  use eulerlagrange_callback3d
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
  public :: eulerlagrange_calcJacobianThetaScheme
  public :: eulerlagrange_calcPrecondThetaScheme
  public :: eulerlagrange_calcResidualThetaScheme
  public :: eulerlagrange_calcRhsRungeKuttaScheme
  public :: eulerlagrange_calcRhsThetaScheme
  public :: eulerlagrange_nlsolverCallback
  public :: eulerlagrange_setBoundaryConditions
  public :: eulerlagrange_calcLinearisedFCT
  public :: eulerlagrange_calcFluxFCT
  public :: eulerlagrange_calcCorrectionFCT
  public :: eulerlagrange_limitEdgewiseVelocity
  public :: eulerlagrange_limitEdgewiseMomentum
  public :: eulerlagrange_coeffVectorFE
  public :: eulerlagrange_coeffVectorAnalytic

contains

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_nlsolverCallback(rproblemLevel, rtimestep,&
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
      call eulerlagrange_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rcollection)
    end if


    ! Do we have to calculate the constant right-hand side?
    ! --------------------------------------------------------------------------
    if ((iand(ioperationSpec, NLSOL_OPSPEC_CALCRHS)  .ne. 0)) then

      ! Compute the right-hand side
      call eulerlagrange_calcRhsRungeKuttaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rrhs, istep,&
          rcollection)
    end if


    ! Do we have to calculate the residual?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      if (istep .eq. 0) then
        ! Compute the constant right-hand side
        call eulerlagrange_calcRhsThetaScheme(rproblemLevel, rtimestep,&
            rsolver, rsolution0, rrhs, rcollection, rsource)
      end if

      ! Compute the residual
      call eulerlagrange_calcResidualThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rrhs, rres, istep,&
          rcollection)
    end if


    ! Do we have to impose boundary conditions?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      ! Impose boundary conditions
      call eulerlagrange_setBoundaryConditions(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rres, rcollection)
    end if


    ! Set status flag
    istatus = 0

  end subroutine eulerlagrange_nlsolverCallback

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
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
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
          call sys_halt()
        end select

      case (SYSTEM_BLOCKFORMAT)

        select case(imasstype)
        case (MASS_LUMPED)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, eulerlagrange_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(&
                rproblemLevel%Rmatrix(lumpedMassMatrix), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                .false., .false., .true., .true.)
          end do

        case (MASS_CONSISTENT)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, eulerlagrange_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(&
                rproblemLevel%Rmatrix(consistentMassMatrix), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                .false., .false., .true., .true.)
          end do

        case DEFAULT
          call output_line('Empty system matrix is invalid!',&
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
          call sys_halt()
        end select

      case DEFAULT
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
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
              eulerlagrange_calcMatrixDiagonalDiag1d,&
              eulerlagrange_calcMatrixGalerkinDiag1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonalDiag2d,&
              eulerlagrange_calcMatrixGalerkinDiag2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonalDiag3d,&
              eulerlagrange_calcMatrixGalerkinDiag3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonalDiag1d,&
              eulerlagrange_calcMatrixScalarDissDiag1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonalDiag2d,&
              eulerlagrange_calcMatrixScalarDissDiag2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonalDiag3d,&
              eulerlagrange_calcMatrixScalarDissDiag3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecond')
          call sys_halt()
        end select


      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonalDiag1d,&
              eulerlagrange_calcMatrixTensorDissDiag1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonalDiag2d,&
              eulerlagrange_calcMatrixTensorDissDiag2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonalDiag3d,&
              eulerlagrange_calcMatrixTensorDissDiag3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with the Rusanov dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonalDiag1d,&
              eulerlagrange_calcMatrixRusanovDiag1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonalDiag2d,&
              eulerlagrange_calcMatrixRusanovDiag2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonalDiag3d,&
              eulerlagrange_calcMatrixRusanovDiag3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
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
              eulerlagrange_calcMatrixDiagonal1d, eulerlagrange_calcMatrixGalerkin1d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonal2d, eulerlagrange_calcMatrixGalerkin2d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonal3d, eulerlagrange_calcMatrixGalerkin3d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonal1d,&
              eulerlagrange_calcMatrixScalarDiss1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonal2d,&
              eulerlagrange_calcMatrixScalarDiss2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonal3d,&
              eulerlagrange_calcMatrixScalarDiss3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonal1d,&
              eulerlagrange_calcMatrixTensorDiss1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonal2d,&
              eulerlagrange_calcMatrixTensorDiss2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonal3d,&
              eulerlagrange_calcMatrixTensorDiss3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with the Rusanov dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonal1d, eulerlagrange_calcMatrixRusanov1d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonal2d, eulerlagrange_calcMatrixRusanov2d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcMatrixDiagonal3d, eulerlagrange_calcMatrixRusanov3d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
        call sys_halt()
      end select


    case DEFAULT
      call output_line('Invalid type of flow coupling!',&
          OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
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

        do ivar = 1, eulerlagrange_getNVAR(rproblemLevel)
          do jvar = 1, eulerlagrange_getNVAR(rproblemLevel)

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

        do ivar = 1, eulerlagrange_getNVAR(rproblemLevel)
          do jvar = 1, eulerlagrange_getNVAR(rproblemLevel)

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
          OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcPrecondThetaScheme')
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

  end subroutine eulerlagrange_calcPrecondThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcJacobianThetaScheme(rproblemLevel, rtimestep,&
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

  end subroutine eulerlagrange_calcJacobianThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcRhsThetaScheme(rproblemLevel, rtimestep,&
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
                eulerlagrange_calcFluxGalerkin1d, dscale, .true., rrhs)

          case (NDIM2D)
            call gfsys_buildDivVector(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                eulerlagrange_calcFluxGalerkin2d, dscale, .true., rrhs)

          case (NDIM3D)
            call gfsys_buildDivVector(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                eulerlagrange_calcFluxGalerkin3d, dscale, .true., rrhs)
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
                  eulerlagrange_calcFluxGalerkin1d, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxGalerkin2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxGalerkin3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_SCALAR)

            ! Assemble divergence operator with scalar dissipation

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxScalarDiss1d, dscale, .true. , rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxScalarDiss2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxScalarDiss3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_SCALAR_DSPLIT)

            ! Assemble divergence operator with scalar dissipation
            ! adopting dimensional splitting

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxScalarDiss1d, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxScalarDissDiSp2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxScalarDissDiSp3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_TENSOR)

            ! Assemble divergence operator with tensorial dissipation

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxTensorDiss1d, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxTensorDiss2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxTensorDiss3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_TENSOR_DSPLIT)

            ! Assemble divergence operator with tensorial dissipation
            ! adopting dimensional splitting

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxTensorDiss1d, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxTensorDissDiSp2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxTensorDissDiSp3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_RUSANOV)

            ! Assemble divergence operator with Rusanov flux

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxRusanov1d, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxRusanov2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxRusanov3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_RUSANOV_DSPLIT)

            ! Assemble divergence operator with Rusanov flux
            ! adopting dimensional splitting

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxRusanov1d, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxRusanovDiSp2d, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  eulerlagrange_calcFluxRusanovDiSp3d, dscale, .true., rrhs)
            end select

          case DEFAULT
            call output_line('Invalid type of dissipation!',&
                OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcRhsThetaScheme')
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
                eulerlagrange_calcFluxGalerkinNoBdr1d,&
                eulerlagrange_calcCharacteristics1d, dscale, .true., rrhs)

          case (NDIM2D)
            call gfsys_buildDivVectorTVD(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                eulerlagrange_calcFluxGalerkinNoBdr2d,&
                eulerlagrange_calcCharacteristics2d, dscale, .true., rrhs)

          case (NDIM3D)
            call gfsys_buildDivVectorTVD(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                eulerlagrange_calcFluxGalerkinNoBdr3d,&
                eulerlagrange_calcCharacteristics3d, dscale, .true., rrhs)
          end select

        case DEFAULT
          call output_line('Invalid type of stabilisation!',&
              OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcRhsThetaScheme')
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

  end subroutine eulerlagrange_calcRhsThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcResidualThetaScheme(rproblemLevel, rtimestep,&
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
            eulerlagrange_calcFluxGalerkin1d, dscale, .false., rres)

      case (NDIM2D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            eulerlagrange_calcFluxGalerkin2d, dscale, .false., rres)

      case (NDIM3D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            eulerlagrange_calcFluxGalerkin3d, dscale, .false., rres)
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
              eulerlagrange_calcFluxGalerkin1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxGalerkin2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxGalerkin3d, dscale, .false., rres)
        end select

      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxScalarDiss1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxScalarDiss2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxScalarDiss3d, dscale, .false., rres)
        end select

      case (DISSIPATION_SCALAR_DSPLIT)

        ! Assemble divergence operator with scalar dissipation
        ! adopting dimensional splitting

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxScalarDiss1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxScalarDissDiSp2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxScalarDissDiSp3d, dscale, .false., rres)
        end select

      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxTensorDiss1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxTensorDiss2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxTensorDiss3d, dscale, .false., rres)
        end select

      case (DISSIPATION_TENSOR_DSPLIT)

        ! Assemble divergence operator with tensorial dissipation
        ! adopting dimensional splitting

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxTensorDiss1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxTensorDissDiSp2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxTensorDissDiSp3d, dscale, .false., rres)
        end select

      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with Rusanov flux

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxRusanov1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxRusanov2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxRusanov3d, dscale, .false., rres)
        end select

      case (DISSIPATION_RUSANOV_DSPLIT)

        ! Assemble divergence operator with Rusanov flux adopting
        ! dimensional splitting

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxRusanov1d, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxRusanovDiSp2d, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxRusanovDiSp3d, dscale, .false., rres)
        end select

      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcResidualThetaScheme')
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
            eulerlagrange_calcFluxGalerkinNoBdr1d,&
            eulerlagrange_calcCharacteristics1d, dscale, .false., rres)

      case (NDIM2D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            eulerlagrange_calcFluxGalerkinNoBdr2d,&
            eulerlagrange_calcCharacteristics2d, dscale, .false., rres)

      case (NDIM3D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            eulerlagrange_calcFluxGalerkinNoBdr3d,&
            eulerlagrange_calcCharacteristics3d, dscale , .false., rres)
      end select

    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcResidualThetaScheme')
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
      call eulerlagrange_calcFluxFCT(rproblemLevel, rsolution, rsolution,&
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
      call eulerlagrange_calcCorrectionFCT(rproblemLevel, p_rpredictor,&
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

  end subroutine eulerlagrange_calcResidualThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcRhsRungeKuttaScheme(rproblemLevel, rtimestep,&
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
            eulerlagrange_calcFluxGalerkin1d, dscale, .false., rrhs)

      case (NDIM2D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            eulerlagrange_calcFluxGalerkin2d, dscale, .false., rrhs)

      case (NDIM3D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            eulerlagrange_calcFluxGalerkin3d, dscale, .false., rrhs)
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
              eulerlagrange_calcFluxGalerkin1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxGalerkin2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxGalerkin3d, dscale, .false., rrhs)
        end select

      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxScalarDiss1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxScalarDiss2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxScalarDiss3d, dscale, .false., rrhs)

        end select

      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxTensorDiss1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxTensorDiss2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxTensorDiss3d, dscale, .false., rrhs)

        end select

      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with Rusanov flux

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxRusanov1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxRusanov2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              eulerlagrange_calcFluxRusanov3d, dscale, .false., rrhs)

        end select

      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcRhsRungeKuttaScheme')
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
            eulerlagrange_calcFluxGalerkinNoBdr1d,&
            eulerlagrange_calcCharacteristics1d, dscale, .false., rrhs)

      case (NDIM2D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            eulerlagrange_calcFluxGalerkinNoBdr2d,&
            eulerlagrange_calcCharacteristics2d, dscale, .false., rrhs)

      case (NDIM3D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            eulerlagrange_calcFluxGalerkinNoBdr3d,&
            eulerlagrange_calcCharacteristics3d, dscale, .false., rrhs)

      end select

    case DEFAULT
      call output_line('Invalid type of stabilization!',&
          OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcRhsRungeKuttaScheme')
      call sys_halt()
    end select

    ! Stop time measurement for global operator
    call stat_stopTimer(p_rtimer)

  end subroutine eulerlagrange_calcRhsRungeKuttaScheme

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_setBoundaryConditions(rproblemLevel, rtimestep,&
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
          OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_setBoundaryConditions')
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
          eulerlagrange_calcBoundaryvalues1d, istatus)

    case (NDIM2D)
      call bdrf_filterSolution(&
          rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixBlock(imatrix),&
          rsolution, rres, rsolution0, rtimestep%dTime,&
          eulerlagrange_calcBoundaryvalues2d, istatus)

    case (NDIM3D)
      call bdrf_filterSolution(&
          rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixBlock(imatrix),&
          rsolution, rres, rsolution0, rtimestep%dTime,&
          eulerlagrange_calcBoundaryvalues3d, istatus)

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_setBoundaryConditions')
      call sys_halt()
    end select

  end subroutine eulerlagrange_setBoundaryConditions

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcLinearisedFCT(rbdrCond, rproblemLevel,&
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
        call eulerlagrange_getVariable(rsolution, trim(sfailsafevariable),&
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
    call eulerlagrange_calcRhsThetaScheme(rproblemLevel, rtimestepAux,&
        rsolver, rsolution, p_rpredictor, rcollection, rsource)
    
    ! Compute low-order predictor
    call lsysbl_invertedDiagMatVec(&
        rproblemLevel%Rmatrix(lumpedMassMatrix),&
        p_rpredictor, 1.0_DP, p_rpredictor)

    ! Compute the raw antidiffusive fluxes
    call eulerlagrange_calcFluxFCT(rproblemLevel, p_rpredictor,&
        rsolution, 0.0_DP, 1.0_DP, 1.0_DP, .true., rcollection)
    
    ! Apply FEM-FCT correction
    call eulerlagrange_calcCorrectionFCT(rproblemLevel,&
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
          call eulerlagrange_getVariable(rsolution, trim(sfailsafevariable), rvectorScalar)
          
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
          rtimestep%dTime, eulerlagrange_calcBoundaryvalues1d)

    case (NDIM2D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dTime, eulerlagrange_calcBoundaryvalues2d)

    case (NDIM3D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dTime, eulerlagrange_calcBoundaryvalues3d)
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

  end subroutine eulerlagrange_calcLinearisedFCT

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcFluxFCT(rproblemLevel, rsolution1, rsolution2,&
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
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTScalarDiss1d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTScalarDiss1d,&
              theta, tstep, dscale, binit)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTScalarDiss2d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTScalarDiss2d,&
              theta, tstep, dscale, binit)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTScalarDiss3d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTScalarDiss3d,&
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
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTTensorDiss1d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTTensorDiss1d,&
              theta, tstep, dscale, binit)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTTensorDiss2d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTTensorDiss2d,&
              theta, tstep, dscale, binit)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTTensorDiss3d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTTensorDiss3d,&
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
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTRusanov1d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTRusanov1d,&
              theta, tstep, dscale, binit)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTRusanov2d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTRusanov2d,&
              theta, tstep, dscale, binit)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTRusanov3d,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, rsolution2, eulerlagrange_calcFluxFCTRusanov3d,&
              theta, tstep, dscale, binit)
        end if
      end select


    case DEFAULT
      call output_line('Invalid type of dissipation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcFluxFCT')
      call sys_halt()
    end select

  end subroutine eulerlagrange_calcFluxFCT

  !*****************************************************************************

!<subroutine>

  subroutine eulerlagrange_calcCorrectionFCT(rproblemLevel, rsolution, &
      dscale, bclear, ioperationSpec, rresidual, rcollection, rafcstab)

!<description>
    ! This subroutine calculates the raw antidiffusive fluxes for
    ! the different FCT algorithms
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

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
    ! residual vector
    type(t_vectorBlock), intent(inout) :: rresidual

    ! collection structure
    type(t_collection), intent(inout) :: rcollection

    ! OPTIONAL: stabilisation structure
    type(t_afcstab), intent(inout), optional, target :: rafcstab
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_afcstab), pointer :: p_rafcstab
    character(len=SYS_STRLEN) :: slimitingvariable
    integer(I32) :: iopSpec
    integer :: inviscidAFC, lumpedMassMatrix
    integer :: ivariable, nvariable, nvartransformed


    ! Get parameters from parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1),&
        'lumpedmassmatrix', lumpedmassmatrix)
    
    ! Set pointer to stabilisation structure
    if (present(rafcstab)) then
      p_rafcstab => rafcstab
    else
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1),&
          'inviscidAFC', inviscidAFC)
      p_rafcstab => rproblemLevel%Rafcstab(inviscidAFC)
    end if

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

      ! Get number of variables to be limited simultaneously
      nvartransformed = eulerlagrange_getNVARtransformed(rproblemLevel, slimitingvariable)

      ! What type of flux transformation is applied?
      if (trim(slimitingvariable) .eq. 'density') then

        ! Apply FEM-FCT algorithm for density fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxDensity1d, eulerlagrange_trafoDiffDensity1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxDensity2d, eulerlagrange_trafoDiffDensity2d)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxDensity3d, eulerlagrange_trafoDiffDensity3d)
        end select

      elseif (trim(slimitingvariable) .eq. 'energy') then

        ! Apply FEM-FCT algorithm for energy fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxEnergy1d, eulerlagrange_trafoDiffEnergy1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxEnergy2d, eulerlagrange_trafoDiffEnergy2d)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxEnergy3d, eulerlagrange_trafoDiffEnergy3d)
        end select

      elseif (trim(slimitingvariable) .eq. 'pressure') then

        ! Apply FEM-FCT algorithm for pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxPressure1d, eulerlagrange_trafoDiffPressure1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxPressure2d, eulerlagrange_trafoDiffPressure2d)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxPressure3d, eulerlagrange_trafoDiffPressure3d)
        end select

      elseif (trim(slimitingvariable) .eq. 'velocity') then

        ! Apply FEM-FCT algorithm for velocity fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxVelocity1d, eulerlagrange_trafoDiffVelocity1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxVelocity2d, eulerlagrange_trafoDiffVelocity2d,&
              fcb_limitEdgewise=eulerlagrange_limitEdgewiseVelocity)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxVelocity3d, eulerlagrange_trafoDiffVelocity3d,&
              fcb_limitEdgewise=eulerlagrange_limitEdgewiseVelocity)
        end select

      elseif (trim(slimitingvariable) .eq. 'momentum') then

        ! Apply FEM-FCT algorithm for momentum fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxMomentum1d, eulerlagrange_trafoDiffMomentum1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxMomentum2d, eulerlagrange_trafoDiffMomentum2d,&
              fcb_limitEdgewise=eulerlagrange_limitEdgewiseMomentum)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxMomentum3d, eulerlagrange_trafoDiffMomentum3d,&
              fcb_limitEdgewise=eulerlagrange_limitEdgewiseMomentum)
        end select
        
      elseif (trim(slimitingvariable) .eq. 'density,energy') then

        ! Apply FEM-FCT algorithm for density and energy fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxDenEng1d, eulerlagrange_trafoDiffDenEng1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxDenEng2d, eulerlagrange_trafoDiffDenEng2d)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxDenEng3d, eulerlagrange_trafoDiffDenEng3d)
        end select

      elseif (trim(slimitingvariable) .eq. 'density,pressure') then

        ! Apply FEM-FCT algorithm for density and pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxDenPre1d, eulerlagrange_trafoDiffDenPre1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxDenPre2d, eulerlagrange_trafoDiffDenPre2d)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxDenPre3d, eulerlagrange_trafoDiffDenPre3d)
        end select

      elseif (trim(slimitingvariable) .eq. 'density,energy,momentum') then

        nvartransformed = eulerlagrange_getNVARtransformed(rproblemLevel, slimitingvariable)

        ! Apply FEM-FCT algorithm for full conservative fluxes
        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
            rsolution, dscale, bclear, iopSpec, rresidual)

      elseif (trim(slimitingvariable) .eq. 'density,pressure,velocity') then

        nvartransformed = eulerlagrange_getNVARtransformed(rproblemLevel, slimitingvariable)

        ! Apply FEM-FCT algorithm for density, velocity and pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxDenPreVel1d, eulerlagrange_trafoDiffDenPreVel1d)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxDenPreVel2d, eulerlagrange_trafoDiffDenPreVel2d)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              eulerlagrange_trafoFluxDenPreVel3d, eulerlagrange_trafoDiffDenPreVel3d)
        end select

      elseif (trim(slimitingvariable) .eq. 'none') then
        
        ! Apply raw antidiffusive fluxes without correction
        iopSpec = ioperationSpec
        iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_ADINCREMENTS))
        iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_BOUNDS))
        iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_LIMITNODAL))
        iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_LIMITEDGE))
        
        ! Enforce existence of edgewise correction factors
        p_rafcstab%iSpec = ior(p_rafcstab%iSpec, AFCSTAB_HAS_EDGELIMITER)

        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
            rsolution, dscale, bclear, iopSpec, rresidual)
        
        ! Nothing more needs to be done
        return

      else
        call output_line('Invalid type of flux transformation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_calcCorrectionFCT')
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
          rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
          rsolution, dscale, bclear, iopSpec, rresidual)
    end if

  end subroutine eulerlagrange_calcCorrectionFCT

  ! ***************************************************************************

!<subroutine>

  subroutine eulerlagrange_limitEdgewiseVelocity(IverticesAtEdge, NEDGE, NEQ,&
      NVAR, NVARtransformed, ndim1, ndim2, Dx, Dflux, Drp, Drm, Dalpha,&
      fcb_calcFluxTransformation, Dflux0, rcollection)

!<description>
    ! This subroutine computes the edgewise correction factors
    ! for the velocity vector in synchronized fashion.
    ! Note that this subroutine is designed for vectors in
    ! interleave and block format, whereby the concrete format
    ! is determined by means of the variables ndim1 and ndim2.
!</description>

!<input>
    ! Number of edges
    integer, intent(in) :: NEDGE
    
    ! Number of nodes
    integer, intent(in) :: NEQ
    
    ! Number of solution variables
    integer, intent(IN) :: NVAR

    ! Number of transformed variables
    integer, intent(IN) :: NVARtransformed

    ! Dimensions of the solution vector
    integer, intent(in) :: ndim1, ndim2

    ! Solution used for flux transformation
    real(DP), dimension(ndim1,ndim2), intent(in) :: Dx

    ! Raw antidiffusive flux
    real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux

    ! Nodal correction factors
    real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! OPTIONAL: callback function to compute variable transformation
    include '../../../../../kernel/PDEOperators/intf_gfsyscallback.inc'
    optional :: fcb_calcFluxTransformation

    ! OPTIONAL: Antidiffusive flux for constraining
    real(DP), dimension(NVAR,NEDGE), intent(in), optional :: Dflux0
!</intput>

!<inputoutput>
    ! Edgewise correction factors
    real(DP), dimension(:), intent(inout) :: Dalpha

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(NVARtransformed) :: F_ij,F_ji,R_ij,R_ji,Ui,Uj,Uij
    real(DP) :: alpha_ij
    integer :: iedge,i,j
    
    ! Do we have to use the explicit fluxes as constraints?
    if (present(Dflux0)) then
      print *, "Not implemented yet"
      stop
      
    else

      if ((ndim1 .eq. NVAR) .and. (ndim2 .eq. NEQ)) then

        !-----------------------------------------------------------------------
        ! The vector is given in interleave format
        !-----------------------------------------------------------------------

        ! Loop over all edges
        do iedge = 1, NEDGE
          
          ! Get node numbers and matrix positions
          i  = IverticesAtEdge(1, iedge)
          j  = IverticesAtEdge(2, iedge)
          
          ! Compute transformed velocity fluxes
          call fcb_calcFluxTransformation(&
              Dx(:,i), Dx(:,j), Dflux(:,iedge), F_ij, F_ji)
          
          ! Compute nodal correction factors
          R_ij = merge(Drp(:,i), Drm(:,i), F_ij .ge. 0.0_DP)
          R_ji = merge(Drp(:,j), Drm(:,j), F_ji .ge. 0.0_DP)
          
          ! Compute edgewise correction factors
          R_ij = min(R_ij, R_ji)
          
          ! Compute velocity average
          Uij = 0.5_DP*(Dx(2:NVARtransformed+1,i)/Dx(1,i)+&
                        Dx(2:NVARtransformed+1,j)/Dx(1,j))
          
          ! Compute correction factor
          alpha_ij = sum(R_ij*Uij*Uij)/(sum(Uij*Uij)+SYS_EPSREAL)
          
          ! Compute multiplicative correction factor
          Dalpha(iedge) = Dalpha(iedge) *alpha_ij
        end do

      elseif ((ndim1 .eq. NEQ) .and. (ndim2 .eq. NVAR)) then

        !-----------------------------------------------------------------------
        ! The vector is given in block format
        !-----------------------------------------------------------------------

        ! Loop over all edges
        do iedge = 1, NEDGE
          
          ! Get node numbers and matrix positions
          i  = IverticesAtEdge(1, iedge)
          j  = IverticesAtEdge(2, iedge)
          
          ! Compute transformed velocity fluxes
          call fcb_calcFluxTransformation(&
              Dx(i,:), Dx(j,:), Dflux(:,iedge), F_ij, F_ji)
          
          ! Compute nodal correction factors
          R_ij = merge(Drp(:,i), Drm(:,i), F_ij .ge. 0.0_DP)
          R_ji = merge(Drp(:,j), Drm(:,j), F_ji .ge. 0.0_DP)
          
          ! Compute edgewise correction factors
          R_ij = min(R_ij, R_ji)
          
          ! Compute velocity average
          Uij = 0.5_DP*(Dx(i,2:NVARtransformed+1)/Dx(i,1)+&
                        Dx(j,2:NVARtransformed+1)/Dx(j,1))
          
          ! Compute correction factor
          alpha_ij = sum(R_ij*Uij*Uij)/(sum(Uij*Uij)+SYS_EPSREAL)
          
          ! Compute multiplicative correction factor
          Dalpha(iedge) = Dalpha(iedge) *alpha_ij
        end do

      else
        
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_limitEdgewiseVelocity')
        call sys_halt()
        
      end if
    end if
    
  end subroutine eulerlagrange_limitEdgewiseVelocity

  ! ***************************************************************************

!<subroutine>

  subroutine eulerlagrange_limitEdgewiseMomentum(IverticesAtEdge, NEDGE, NEQ,&
      NVAR, NVARtransformed, ndim1, ndim2, Dx, Dflux, Drp, Drm, Dalpha,&
      fcb_calcFluxTransformation, Dflux0, rcollection)

!<description>
    ! This subroutine computes the edgewise correction factors
    ! for the momentum vector in synchronized fashion.
    ! Note that this subroutine is designed for vectors in
    ! interleave and block format, whereby the concrete format
    ! is determined by means of the variables ndim1 and ndim2.
!</description>

!<input>
    ! Number of edges
    integer, intent(in) :: NEDGE
    
    ! Number of nodes
    integer, intent(in) :: NEQ
    
    ! Number of solution variables
    integer, intent(IN) :: NVAR

    ! Number of transformed variables
    integer, intent(IN) :: NVARtransformed

    ! Dimensions of the solution vector
    integer, intent(in) :: ndim1, ndim2

    ! Solution used for flux transformation
    real(DP), dimension(ndim1,ndim2), intent(in) :: Dx

    ! Raw antidiffusive flux
    real(DP), dimension(NVAR,NEDGE), intent(in) :: Dflux

    ! Nodal correction factors
    real(DP), dimension(NVARtransformed,NEQ), intent(in) :: Drp,Drm

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! OPTIONAL: callback function to compute variable transformation
    include '../../../../../kernel/PDEOperators/intf_gfsyscallback.inc'
    optional :: fcb_calcFluxTransformation

    ! OPTIONAL: Antidiffusive flux for constraining
    real(DP), dimension(NVAR,NEDGE), intent(in), optional :: Dflux0
!</intput>

!<inputoutput>
    ! Edgewise correction factors
    real(DP), dimension(:), intent(inout) :: Dalpha

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(NVARtransformed) :: F_ij,F_ji,R_ij,R_ji,Ui,Uj,Uij
    real(DP) :: alpha_ij
    integer :: iedge,i,j
    
    ! Do we have to use the explicit fluxes as constraints?
    if (present(Dflux0)) then
      print *, "Not implemented yet"
      stop
      
    else

      if ((ndim1 .eq. NVAR) .and. (ndim2 .eq. NEQ)) then

        !-----------------------------------------------------------------------
        ! The vector is given in interleave format
        !-----------------------------------------------------------------------

        ! Loop over all edges
        do iedge = 1, NEDGE
          
          ! Get node numbers and matrix positions
          i  = IverticesAtEdge(1, iedge)
          j  = IverticesAtEdge(2, iedge)
          
          ! Compute transformed velocity fluxes
          call fcb_calcFluxTransformation(&
              Dx(:,i), Dx(:,j), Dflux(:,iedge), F_ij, F_ji)
          
          ! Compute nodal correction factors
          R_ij = merge(Drp(:,i), Drm(:,i), F_ij .ge. 0.0_DP)
          R_ji = merge(Drp(:,j), Drm(:,j), F_ji .ge. 0.0_DP)
          
          ! Compute edgewise correction factors
          R_ij = min(R_ij, R_ji)
          
          ! Compute momentum average
          Uij = 0.5_DP*(Dx(2:NVARtransformed+1,i)+&
                        Dx(2:NVARtransformed+1,j))
          
          ! Compute correction factor
          alpha_ij = sum(R_ij*Uij*Uij)/(sum(Uij*Uij)+SYS_EPSREAL)
          
          ! Compute multiplicative correction factor
          Dalpha(iedge) = Dalpha(iedge) *alpha_ij
        end do

      elseif ((ndim1 .eq. NEQ) .and. (ndim2 .eq. NVAR)) then

        !-----------------------------------------------------------------------
        ! The vector is goven in block format
        !-----------------------------------------------------------------------

        ! Loop over all edges
        do iedge = 1, NEDGE
          
          ! Get node numbers and matrix positions
          i  = IverticesAtEdge(1, iedge)
          j  = IverticesAtEdge(2, iedge)
          
          ! Compute transformed velocity fluxes
          call fcb_calcFluxTransformation(&
              Dx(i,:), Dx(j,:), Dflux(:,iedge), F_ij, F_ji)
          
          ! Compute nodal correction factors
          R_ij = merge(Drp(:,i), Drm(:,i), F_ij .ge. 0.0_DP)
          R_ji = merge(Drp(:,j), Drm(:,j), F_ji .ge. 0.0_DP)
          
          ! Compute edgewise correction factors
          R_ij = min(R_ij, R_ji)
          
          ! Compute momentum average
          Uij = 0.5_DP*(Dx(i,2:NVARtransformed+1)+&
                        Dx(j,2:NVARtransformed+1))
          
          ! Compute correction factor
          alpha_ij = sum(R_ij*Uij*Uij)/(sum(Uij*Uij)+SYS_EPSREAL)
          
          ! Compute multiplicative correction factor
          Dalpha(iedge) = Dalpha(iedge) *alpha_ij
        end do

      else

        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_limitEdgewiseMomentum')
        call sys_halt()
        
      end if
    end if

  end subroutine eulerlagrange_limitEdgewiseMomentum

  ! ***************************************************************************

!<subroutine>

  subroutine eulerlagrange_coeffVectorFE(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use collection
    use derivatives
    use domainintegration
    use feevaluation
    use fsystem
    use scalarpde
    use spatialdiscretisation
    use triangulation

!<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points.
    !
    ! The following data must be passed to this routine in the collection in order 
    ! to work correctly:
    !
    ! IquickAccess(1) = systemFormat
    ! IquickAccess(2) = ivar
    ! p_rvectorQuickAccess1 => evaluation solution vector
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements test in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional :: rcollection
!</input>
  
!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(DP), dimension(:), pointer :: Ddata
    integer :: isystemFormat,ivar,iel,ipoint
    
    if (.not. present(rcollection)) then
      Dcoefficients(:,:,:) = 0.0_DP
      return
    end if

    ! Get the parameters from the collection
    isystemFormat = rcollection%IquickAccess(1)
    ivar = rcollection%IquickAccess(2)
    p_rvector => rcollection%p_rvectorQuickAccess1
    
    ! What type of system format are we?
    select case(isystemFormat)
   
    case(SYSTEM_INTERLEAVEFORMAT)
      
      print *, "Not available"
      stop
      
    case(SYSTEM_BLOCKFORMAT)
      
      ! Allocate temporal array
      allocate(Ddata(npointsPerElement))
      
      ! Loop over all elements
      do iel = 1, nelements
        
        ! Evaluate solution in cubature points
        call fevl_evaluate(DER_FUNC, Ddata,&
            p_rvector%RvectorBlock(ivar), Dpoints(:,:,iel))
        
        ! Loop over all cubature points
        do ipoint = 1, npointsPerElement
          Dcoefficients(1,ipoint,iel) = Ddata(ipoint)
        end do
      end do
      
      ! Deallocate temporal array
      deallocate(Ddata)

    case default
      call output_line ('Invalid system format!', &
          OU_CLASS_ERROR,OU_MODE_STD,'eulerlagrange_coeffVectorFE')
      call sys_halt()
    end select
    
  end subroutine eulerlagrange_coeffVectorFE

  ! ***************************************************************************

!<subroutine>

  subroutine eulerlagrange_coeffVectorAnalytic(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fparser
    use scalarpde
    use spatialdiscretisation
    use triangulation

!<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points.
    !
    ! The following data must be passed to this routine in the collection in order 
    ! to work correctly:
    !
    ! DquickAccess(1)          = dtime % simulation time
    ! IquickAccess(itermCount) = icomp % number of the function to be evaluated
    ! SquickAccess(1)          = name of the function parser
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in) :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements test in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional :: rcollection
!</input>
  
!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
!</output>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dtime
    integer :: itermCount, ipoint, iel, ndim, icomp
    
    
    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(1)))
    
    ! This subroutine assumes that the first quick access double value
    ! holds the simulation time
    dtime  = rcollection%DquickAccess(1)

    ! Loop over all components of the linear form
    do itermCount = 1, ubound(Dcoefficients,1)

      ! Moreover, this subroutine assumes the quick access integer
      ! 'itermCount' holds the number of the function to be evaluated
      icomp = rcollection%IquickAccess(itermCount)

      if (dtime < 0.0) then
        
        ! Evaluate all coefficients using the function parser
        do iel = 1, nelements
          call fparser_evalFunction(p_rfparser, icomp, 2,&
              Dpoints(:,:,iel), Dcoefficients(itermCount,:,iel))
        end do

      else

        ! Initialize values
        Dvalue = 0.0_DP
        Dvalue(NDIM3D+1) = dtime
        
        ! Set number of spatial dimensions
        ndim = size(Dpoints, 1)
        
        do iel = 1, nelements
          do ipoint = 1, npointsPerElement
            
            ! Set values for function parser
            Dvalue(1:ndim) = Dpoints(:, ipoint, iel)
            
            ! Evaluate function parser
            call fparser_evalFunction(p_rfparser, icomp, Dvalue,&
                Dcoefficients(itermCount,ipoint,iel))
          end do
        end do

      end if

    end do ! itermCount
    
  end subroutine eulerlagrange_coeffVectorAnalytic

  !*****************************************************************************

!<subroutine>

  pure subroutine eulerlagrange_calcMassFluxFCT(&
      U1_i, U1_j, U2_i, U2_j, C_ij, C_ji,&
      i, j, dscale1, dscale2, F_ij)

!<description>
    ! This subroutine computes the raw antidiffusive mass fluxes
    ! for FCT algorithms in arbitrary dimension
!</description>

!<input>
    ! local solution at nodes I and J
    real(DP), dimension(:), intent(in) :: U1_i,U1_j,U2_i,U2_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(in) :: C_ij,C_ji

    ! scaling coefficients
    real(DP), intent(in) :: dscale1,dscale2

    ! node numbers
    integer, intent(in) :: i, j
!</input>

!<output>
    ! raw antidiffusive flux
    real(DP), dimension(:), intent(out) :: F_ij
!</output>
!</subroutine>
  
    ! Compute raw antidiffusive mass fluxes
    F_ij = dscale1*(U1_i-U1_j)
    
  end subroutine eulerlagrange_calcMassFluxFCT

end module eulerlagrange_callback

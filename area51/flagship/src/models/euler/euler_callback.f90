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
    integer :: iblock, imassantidiffusiontype

    
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
            call gfsys_buildResidual(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                euler_calcFluxGalerkin1d, dscale, .true., rrhs)
            
          case (NDIM2D)
            call gfsys_buildResidual(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                euler_calcFluxGalerkin2d, dscale, .true., rrhs)
            
          case (NDIM3D)
            call gfsys_buildResidual(&
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
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxGalerkin1d, dscale, .true., rrhs)
              
            case (NDIM2D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxGalerkin2d, dscale, .true., rrhs)
              
            case (NDIM3D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxGalerkin3d, dscale, .true., rrhs)
            end select
            
          case (DISSIPATION_SCALAR)
            
            ! Assemble divergence operator with scalar dissipation
            
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxScalarDiss1d, dscale, .true. , rrhs)
              
            case (NDIM2D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxScalarDiss2d, dscale, .true., rrhs)
              
            case (NDIM3D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxScalarDiss3d, dscale, .true., rrhs)
            end select
            
          case (DISSIPATION_SCALAR_DSPLIT)
            
            ! Assemble divergence operator with scalar dissipation
            ! adopting dimensional splitting
            
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxScalarDiss1d, dscale, .true., rrhs)
              
            case (NDIM2D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxScalarDissDiSp2d, dscale, .true., rrhs)
              
            case (NDIM3D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxScalarDissDiSp3d, dscale, .true., rrhs)
            end select
            
          case (DISSIPATION_TENSOR)

            ! Assemble divergence operator with tensorial dissipation
            
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxTensorDiss1d, dscale, .true., rrhs)
              
            case (NDIM2D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxTensorDiss2d, dscale, .true., rrhs)
              
            case (NDIM3D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxTensorDiss3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_TENSOR_DSPLIT)
            
            ! Assemble divergence operator with tensorial dissipation
            ! adopting dimensional splitting
          
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxTensorDiss1d, dscale, .true., rrhs)
              
            case (NDIM2D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxTensorDissDiSp2d, dscale, .true., rrhs)
              
            case (NDIM3D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxTensorDissDiSp3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_RUSANOV)
            
            ! Assemble divergence operator with Rusanov flux
            
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxRusanov1d, dscale, .true., rrhs)
              
            case (NDIM2D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxRusanov2d, dscale, .true., rrhs)
              
            case (NDIM3D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxRusanov3d, dscale, .true., rrhs)
            end select

          case (DISSIPATION_RUSANOV_DSPLIT)

            ! Assemble divergence operator with Rusanov flux 
            ! adopting dimensional splitting
          
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxRusanov1d, dscale, .true., rrhs)
              
            case (NDIM2D)
              call gfsys_buildResidual(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  euler_calcFluxRusanovDiSp2d, dscale, .true., rrhs)
              
            case (NDIM3D)
              call gfsys_buildResidual(&
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
            call gfsys_buildResidualTVD(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                euler_calcFluxGalerkinNoBdr1d,&
                euler_calcCharacteristics1d, dscale, .true., rrhs)
            
          case (NDIM2D)
            call gfsys_buildResidualTVD(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                euler_calcFluxGalerkinNoBdr2d,&
                euler_calcCharacteristics2d, dscale, .true., rrhs)
            
          case (NDIM3D)
            call gfsys_buildResidualTVD(&
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


    ! What type if stabilization is applied?
    select case(rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation)
    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT)
      
      !-------------------------------------------------------------------------
      ! Initialise the raw antidiffusive fluxes
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1),&
          'imassantidiffusiontype',&
          imassantidiffusiontype)
      
      ! Should we apply consistent mass antidiffusion?
      if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
        call gfsys_buildFluxFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rsolution, rsolution, euler_calcFluxFCTScalarDiss1d,&
            rtimestep%theta, rtimestep%dStep, .true.,&
            rproblemLevel%Rmatrix(consistentMassMatrix))
      else
        call gfsys_buildFluxFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rsolution, rsolution, euler_calcFluxFCTScalarDiss1d,&
            rtimestep%theta, rtimestep%dStep, .true.)
      end if
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

    ! right-hand side vector
    type(t_vectorBlock), intent(in) :: rrhs

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

    ! residual vector
    type(t_vectorBlock), intent(inout) :: rres

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
    integer :: iblock, imassantidiffusiontype
    
    type(t_vectorBlock) :: rulow

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
        call gfsys_buildResidual(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkin1d, dscale, .false., rres)
        
      case (NDIM2D)
        call gfsys_buildResidual(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkin2d, dscale, .false., rres)
        
      case (NDIM3D)
        call gfsys_buildResidual(&
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
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxGalerkin1d, dscale, .false., rres)
          
        case (NDIM2D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxGalerkin2d, dscale, .false., rres)
          
        case (NDIM3D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxGalerkin3d, dscale, .false., rres)
        end select
        
      case (DISSIPATION_SCALAR)
        
        ! Assemble divergence operator with scalar dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss1d, dscale, .false., rres)
          
        case (NDIM2D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss2d, dscale, .false., rres)
          
        case (NDIM3D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss3d, dscale, .false., rres)
        end select
        
      case (DISSIPATION_SCALAR_DSPLIT)
        
        ! Assemble divergence operator with scalar dissipation
        ! adopting dimensional splitting
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss1d, dscale, .false., rres)
          
        case (NDIM2D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDissDiSp2d, dscale, .false., rres)
          
        case (NDIM3D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDissDiSp3d, dscale, .false., rres)
        end select
        
      case (DISSIPATION_TENSOR)
        
        ! Assemble divergence operator with tensorial dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss1d, dscale, .false., rres)
          
        case (NDIM2D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss2d, dscale, .false., rres)
          
        case (NDIM3D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss3d, dscale, .false., rres)
        end select
        
      case (DISSIPATION_TENSOR_DSPLIT)
        
        ! Assemble divergence operator with tensorial dissipation
        ! adopting dimensional splitting
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss1d, dscale, .false., rres)
          
        case (NDIM2D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDissDiSp2d, dscale, .false., rres)
          
        case (NDIM3D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDissDiSp3d, dscale, .false., rres)
        end select

      case (DISSIPATION_RUSANOV)
        
        ! Assemble divergence operator with Rusanov flux
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov1d, dscale, .false., rres)
          
        case (NDIM2D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov2d, dscale, .false., rres)
          
        case (NDIM3D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov3d, dscale, .false., rres)
        end select
        
      case (DISSIPATION_RUSANOV_DSPLIT)
        
        ! Assemble divergence operator with Rusanov flux adopting
        ! dimensional splitting
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov1d, dscale, .false., rres)
          
        case (NDIM2D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanovDiSp2d, dscale, .false., rres)
          
        case (NDIM3D)
          call gfsys_buildResidual(&
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
        call gfsys_buildResidualTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkinNoBdr1d,&
            euler_calcCharacteristics1d, dscale, .false., rres)

      case (NDIM2D)
        call gfsys_buildResidualTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkinNoBdr2d,&
            euler_calcCharacteristics2d, dscale, .false., rres)
        
      case (NDIM3D)
        call gfsys_buildResidualTVD(&
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

    ! What type if stabilization is applied?
    select case(rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation)
    case (AFCSTAB_FEMFCT_CLASSICAL,&
          AFCSTAB_FEMFCT_ITERATIVE,&
          AFCSTAB_FEMFCT_IMPLICIT)
      
      !-------------------------------------------------------------------------
      ! Update the raw antidiffusive fluxes and apply them to the residual
      !-------------------------------------------------------------------------
      
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1),&
          'imassantidiffusiontype', imassantidiffusiontype)
      
      ! Should we apply consistent mass antidiffusion?
      if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
        call gfsys_buildFluxFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rsolution, rsolution, euler_calcFluxFCTScalarDiss1d,&
            rtimestep%theta, rtimestep%dStep, .false.,&
            rproblemLevel%Rmatrix(consistentMassMatrix))
      else
        call gfsys_buildFluxFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rsolution, rsolution, euler_calcFluxFCTScalarDiss1d,&
            rtimestep%theta, rtimestep%dStep, .false.)
      end if
      
      if (rtimestep%theta .lt. 1.0_DP) then

        call lsysbl_duplicateVector(rrhs, rulow,&
            LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_EMPTY)
        
        call lsysbl_invertedDiagMatVec(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rrhs, 1.0_DP, rulow)
        
        call gfsys_buildResidualFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rulow, rtimestep%dStep, .false.,&
            AFCSTAB_FCTALGO_STANDARD,&
            rres, euler_calcTrafoPrimitive1d)

        call lsysbl_releaseVector(rulow)
      else
        if (ite .eq. 0) then
          call gfsys_buildResidualFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution0, rtimestep%dStep, .false.,&
              AFCSTAB_FCTALGO_STANDARD,&
              rres, euler_calcTrafoPrimitive1d)
        else
          call gfsys_buildResidualFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution0, rtimestep%dStep, .false.,&
              AFCSTAB_FCTALGO_STANDARD -&
              AFCSTAB_FCTALGO_BOUNDS,&
              rres, euler_calcTrafoPrimitive1d)
        end if
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
        call gfsys_buildResidual(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkin1d, dscale, .false., rrhs)

      case (NDIM2D)
        call gfsys_buildResidual(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkin2d, dscale, .false., rrhs)

      case (NDIM3D)
        call gfsys_buildResidual(&
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
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxGalerkin1d, dscale, .false., rrhs)
          
        case (NDIM2D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxGalerkin2d, dscale, .false., rrhs)
          
        case (NDIM3D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxGalerkin3d, dscale, .false., rrhs)
        end select
        
      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxScalarDiss3d, dscale, .false., rrhs)

        end select

      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxTensorDiss3d, dscale, .false., rrhs)

        end select

      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with Rusanov flux

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildResidual(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              euler_calcFluxRusanov2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildResidual(&
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
        call gfsys_buildResidualTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkinNoBdr1d,&
            euler_calcCharacteristics1d, dscale, .false., rrhs)

      case (NDIM2D)
        call gfsys_buildResidualTVD(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            euler_calcFluxGalerkinNoBdr2d,&
            euler_calcCharacteristics2d, dscale, .false., rrhs)

      case (NDIM3D)
        call gfsys_buildResidualTVD(&
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

end module euler_callback

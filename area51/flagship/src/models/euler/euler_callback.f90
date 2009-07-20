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
!# 2.) euler_calcPreconditioner
!#     -> Calculates the nonlinear preconditioner
!#
!# 3.) euler_calcJacobian
!#     -> Calculates the Jacobian matrix
!#
!# 4.) euler_applyJacobian
!#     -> Applies the Jacobian matrix to a given vector
!#
!# 5.) euler_calcResidual
!#     -> Calculates the nonlinear residual vector
!#
!# 6.) euler_calcRHS
!#     -> Calculates the right-hand side vector
!#
!# 7.) euler_setBoundary
!#     -> Imposes boundary conditions for nonlinear solver
!#
!# 8.) euler_calcLinearizedFCT
!#     -> Calculates the linearized FCT correction
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
  public :: euler_nlsolverCallback
  public :: euler_calcPreconditioner
  public :: euler_calcJacobian
  public :: euler_applyJacobian
  public :: euler_calcResidual
  public :: euler_calcRHS
  public :: euler_setBoundary
  public :: euler_calcLinearizedFCT

contains

  !*****************************************************************************

!<subroutine>

  subroutine euler_nlsolverCallback(rproblemLevel, rtimestep,&
      rsolver,rsolution, rsolutionInitial, rrhs, rres, istep,&
      ioperationSpec, rcollection, istatus, rb)

!<description>
    ! This subroutine is called by the nonlinear solver and it is responsible
    ! to assemble preconditioner, right-hand side vector, residual vector, etc.
!</description>

!<input>
    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolutionInitial
    
    ! number of solver step
    integer, intent(in) :: istep
    
    ! specifier for operations
    integer(I32), intent(in) :: ioperationSpec

    ! OPTIONAL: constant right-hand side vector
    type(t_vectorBlock), intent(in), optional :: rb
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
      
      call euler_calcPreconditioner(rproblemLevel, rtimestep,&
          rsolver, rsolution, rcollection)
    end if
    
    
    ! Do we have to calculate the constant right-hand side?
    ! --------------------------------------------------------------------------
    if ((iand(ioperationSpec, NLSOL_OPSPEC_CALCRHS)  .ne. 0)) then

      call euler_calcRhs(rproblemLevel, rtimestep, rsolver,&
          rsolution, rsolutionInitial, rrhs, istep, rcollection)
    end if


    ! Do we have to calculate the residual
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then
      
      call euler_calcResidual(rproblemLevel, rtimestep, rsolver,&
          rsolution, rsolutionInitial, rrhs, rres, istep, rcollection)
    end if
    
    
    ! Do we have to impose boundary conditions?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then
      
      call euler_setBoundary(rproblemLevel, rtimestep, rsolver,&
          rsolution, rsolutionInitial, rres, rcollection)
    end if
    
    
    ! Set status flag
    istatus = 0
    
  end subroutine euler_nlsolverCallback

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcPreconditioner(rproblemLevel, rtimestep,&
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
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
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
          call lsyssc_spreadMatrix(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                   rproblemLevel%Rmatrix(systemMatrix))
        case (MASS_CONSISTENT)
          call lsyssc_spreadMatrix(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                   rproblemLevel%Rmatrix(systemMatrix))
        case DEFAULT
          call output_line('Empty system matrix is invalid!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
          call sys_halt()
        end select

      case (SYSTEM_BLOCKFORMAT)
        
        select case(imasstype)
        case (MASS_LUMPED)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, euler_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(rproblemLevel%Rmatrix(lumpedMassMatrix), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                .false., .false., .true., .true.)
          end do
          
        case (MASS_CONSISTENT)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, euler_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(rproblemLevel%Rmatrix(consistentMassMatrix), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                .false., .false., .true., .true.)
          end do

        case DEFAULT
          call output_line('Empty system matrix is invalid!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
          call sys_halt()
        end select
          
      case DEFAULT
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
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

      case (AFCSTAB_GALERKIN)

        ! Assemble divergence operator for standard Galerin scheme

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcMatrixDiagonalDiag1d,&
              euler_calcMatrixGalerkinDiag1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcMatrixDiagonalDiag2d,&
              euler_calcMatrixGalerkinDiag2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcMatrixDiagonalDiag3d,&
              euler_calcMatrixGalerkinDiag3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
          call sys_halt()
        end select
        
        
      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcMatrixDiagonalDiag1d,&
              euler_calcMatrixScalarDissDiag1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcMatrixDiagonalDiag2d,&
              euler_calcMatrixScalarDissDiag2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcMatrixDiagonalDiag3d,&
              euler_calcMatrixScalarDissDiag3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_TENSOR)
        
        ! Assemble divergence operator with tensorial dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcMatrixDiagonalDiag1d,&
              euler_calcMatrixTensorDissDiag1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcMatrixDiagonalDiag2d,&
              euler_calcMatrixTensorDissDiag2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))
          
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcMatrixDiagonalDiag3d,&
              euler_calcMatrixTensorDissDiag3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)
        
        ! Assemble divergence operator with the Rusanov dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcMatrixDiagonalDiag1d,&
              euler_calcMatrixRusanovDiag1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcMatrixDiagonalDiag2d,&
              euler_calcMatrixRusanovDiag2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))
          
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcMatrixDiagonalDiag3d,&
              euler_calcMatrixRusanovDiag3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))
          
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
          call sys_halt()
        end select


      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
        call sys_halt()
      end select
      
      
    case (SYSTEM_ALLCOUPLED)

      !-------------------------------------------------------------------------
      ! Assemble full block transport operator
      !-------------------------------------------------------------------------

      ! What kind of preconditioner is applied?
      select case(isystemPrecond)
        
      case (AFCSTAB_GALERKIN)
        
        ! Assemble divergence operator for standard Galerin scheme
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcMatrixDiagonal1d, euler_calcMatrixGalerkin1d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcMatrixDiagonal2d, euler_calcMatrixGalerkin2d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcMatrixDiagonal3d, euler_calcMatrixGalerkin3d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_SCALAR)
        
        ! Assemble divergence operator with scalar dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcMatrixDiagonal1d,&
              euler_calcMatrixScalarDiss1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcMatrixDiagonal2d,&
              euler_calcMatrixScalarDiss2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcMatrixDiagonal3d,&
              euler_calcMatrixScalarDiss3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
          call sys_halt()
        end select

        
      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcMatrixDiagonal1d,&
              euler_calcMatrixTensorDiss1d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcMatrixDiagonal2d,&
              euler_calcMatrixTensorDiss2d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcMatrixDiagonal3d,&
              euler_calcMatrixTensorDiss3d, 1.0_DP, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with the Rusanov dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcMatrixDiagonal1d, euler_calcMatrixRusanov1d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcMatrixDiagonal2d, euler_calcMatrixRusanov2d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcMatrixDiagonal3d, euler_calcMatrixRusanov3d,&
              1.0_DP, .true., rproblemLevel&
              %RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
          call sys_halt()
        end select
          
        
      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
        call sys_halt()
      end select


    case DEFAULT
      call output_line('Invalid type of flow coupling!',&
          OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
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

        call lsyssc_MatrixLinearComb(rproblemLevel&
            %Rmatrix(lumpedMassMatrix), 1.0_DP, rproblemLevel&
            %Rmatrix(systemMatrix), rtimestep%theta*rtimestep%dStep,&
            rproblemLevel%Rmatrix(systemMatrix), .false., .false.,&
            .true., .true.)

      case (MASS_CONSISTENT)

        !-----------------------------------------------------------------------
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_C)-theta*dt*L $
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(1), 'consistentmassmatrix', consistentMassMatrix)

        call lsyssc_MatrixLinearComb(rproblemLevel&
            %Rmatrix(consistentMassMatrix), 1.0_DP, rproblemLevel&
            %Rmatrix(systemMatrix), rtimestep%theta*rtimestep%dStep,&
            rproblemLevel%Rmatrix(systemMatrix), .false., .false.,&
            .true., .true.)
        
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

              call lsyssc_MatrixLinearComb(rproblemLevel&
                  %Rmatrix(lumpedMassMatrix), 1.0_DP, rproblemLevel&
                  %RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar)&
                  , rtimestep%theta*rtimestep%dStep, rproblemLevel&
                  %RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                  .false., .false., .true., .true.)

            elseif (isystemCoupling .eq. SYSTEM_ALLCOUPLED) then

              call lsyssc_scaleMatrix(rproblemLevel&
                  %RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,jvar),&
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

              call lsyssc_MatrixLinearComb(rproblemLevel&
                  %Rmatrix(consistentMassMatrix), 1.0_DP,&
                  rproblemLevel%RmatrixBlock(systemMatrix)&
                  %RmatrixBlock(ivar,ivar), rtimestep%theta*rtimestep&
                  %dStep, rproblemLevel%RmatrixBlock(systemMatrix)&
                  %RmatrixBlock(ivar,ivar), .false., .false., .true.,.true.)

            elseif (isystemCoupling .eq. SYSTEM_ALLCOUPLED) then

              call lsyssc_scaleMatrix(rproblemLevel&
                  %RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,jvar),&
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
          OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
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
    
  end subroutine euler_calcPreconditioner

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcJacobian(rproblemLevel, rtimestep, rsolver,&
      rsolution, rsolutionInitial, rcollection)

!<description>
    ! This callback subroutine computes the Jacobian matrix for the
    !  compressible Euler/Navier-Stokes equations
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolutionInitial
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

  end subroutine euler_calcJacobian

  !*****************************************************************************

!<subroutine>

  subroutine euler_applyJacobian(rproblemLevel, rx, ry, cx, cy, rcollection)

!<description>
    ! This subroutine applies the (scaled) Jacobian matrix to 
    ! the vector rx and adds the result to the vector ry.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! vector to which the Jacobian should be applied to
    type(t_vectorBlock), intent(in) :: rx

    ! factor by which rx should be scaled
    real(DP), intent(in) :: cx

    ! factor by which ry should be scaled
    real(DP), intent(in) :: cy
!</input>

!<inputoutput>
    ! vector to which the result should be added
    type(t_vectorBlock), intent(inout) :: ry

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    integer :: jacobianMatrix

    ! Get parameters from parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'jacobianMatrix', jacobianMatrix)

    ! We know where the Jacobian matrix is stored and can apply it
    ! by means of matrix vector multiplication
    call lsysbl_blockMatVec(rproblemLevel&
        %RmatrixBlock(jacobianMatrix), rx, ry, cx, cy)

  end subroutine euler_applyJacobian

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcResidual(rproblemLevel, rtimestep, rsolver,&
      rsolution, rsolutionInitial, rrhs, rres, ite, rcollection)

!<description>
    ! This subroutine computes the nonlinear residual vector and the
    ! constant right-hand side (only in the first iteration).
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolutionInitial

    ! number of nonlinear iteration
    integer, intent(in) :: ite
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

    ! residual vector
    type(t_vectorBlock), intent(inout) :: rres

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: consistentMassMatrix, lumpedMassMatrix
    integer :: inviscidAFC, imasstype, idissipationtype
    integer :: iblock

    
    ! Start time measurement for residual/rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection, 'rtimerAssemblyVector')
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Are we in the zeroth iteration?
    if (ite .eq. 0) then

      !-------------------------------------------------------------------------
      ! In the first nonlinear iteration update we compute
      !
      ! - the residual $ r=dt*L(U^n)*U^n+f $ and
      !
      ! - the right hand side $ b=[M+(1-theta)*dt*L(U^n)]*U^n $
      !
      !-------------------------------------------------------------------------
      
      ! Get parameters from parameter list which are required unconditionally
      p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'coeffmatrix_cx', coeffMatrix_CX)
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'coeffmatrix_cy', coeffMatrix_CY)
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'coeffmatrix_cz', coeffMatrix_CZ)
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'inviscidAFC', inviscidAFC)

      select case(rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation)
        
      case (AFCSTAB_GALERKIN)

        !-----------------------------------------------------------------------
        ! Compute the initial high-order residual
        !
        !   $ res = dt*K(U^n)*U^n $
        !-----------------------------------------------------------------------

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcFluxGalerkin1d, rtimestep%dStep, .true., rres)

        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcFluxGalerkin2d, rtimestep%dStep, .true., rres)
          
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcFluxGalerkin3d, rtimestep%dStep, .true., rres)

        end select
        
      case (AFCSTAB_UPWIND,&
            AFCSTAB_FEMFCT_LINEARIZED)

        !-----------------------------------------------------------------------
        ! Compute the initial low-order residual
        !
        !   $ res = dt*L(U^n)*U^n $
        !-----------------------------------------------------------------------
        
        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(1), 'idissipationtype', idissipationtype)

        select case(idissipationtype)

        case (DISSIPATION_SCALAR)

          ! Assemble divergence operator with scalar dissipation

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                euler_calcFluxScalarDiss1d, rtimestep%dStep, .true., rres)

          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                euler_calcFluxScalarDiss2d, rtimestep%dStep, .true., rres)

          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                euler_calcFluxScalarDiss3d, rtimestep%dStep, .true.,&
                rres)

          end select

        case (DISSIPATION_SCALAR_DSPLIT)

          ! Assemble divergence operator with scalar dissipation
          !  adopting dimensional splitting

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                euler_calcFluxScalarDiss1d, rtimestep%dStep, .true., rres)
            
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                euler_calcFluxDSplitScalarDiss2d, rtimestep%dStep, .true., rres)

          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                euler_calcFluxDSplitScalarDiss3d, rtimestep%dStep, .true., rres)

          end select

        case (DISSIPATION_TENSOR)

          ! Assemble divergence operator with tensorial dissipation
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                euler_calcFluxTensorDiss1d, rtimestep%dStep, .true., rres)

          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                euler_calcFluxTensorDiss2d, rtimestep%dStep, .true., rres)

          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                euler_calcFluxTensorDiss3d, rtimestep%dStep, .true., rres)

          end select

        case (DISSIPATION_TENSOR_DSPLIT)

          ! Assemble divergence operator with tensorial dissipation
          !  adopting dimensional splitting
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                euler_calcFluxTensorDiss1d, rtimestep%dStep, .true., rres)

          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                euler_calcFluxDSplitTensorDiss2d, rtimestep%dStep, .true., rres)

          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                euler_calcFluxDSplitTensorDiss3d, rtimestep%dStep, .true., rres)

          end select

        case (DISSIPATION_RUSANOV)

          ! Assemble divergence operator with Rusanov flux
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                euler_calcFluxRusanov1d, rtimestep%dStep, .true., rres)

          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                euler_calcFluxRusanov2d, rtimestep%dStep, .true., rres)

          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                euler_calcFluxRusanov3d, rtimestep%dStep, .true., rres)

          end select

        case (DISSIPATION_RUSANOV_DSPLIT)

          ! Assemble divergence operator with Rusanov flux adopting
          !  dimensional splitting
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                euler_calcFluxRusanov1d, rtimestep%dStep, .true., rres)

          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                euler_calcFluxDSplitRusanov2d, rtimestep%dStep, .true., rres)

          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                euler_calcFluxDSplitRusanov3d, rtimestep%dStep, .true., rres)

          end select

        case DEFAULT
          call output_line('Invalid type of dissipation!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcResidual')
          call sys_halt()
        end select

      case (AFCSTAB_FEMTVD)

        !-----------------------------------------------------------------------
        ! Compute the initial low-order residual + FEM-TVD stabilization
        !
        !   $ res = dt*L(U^n)*U^n + F(U^n) $
        !-----------------------------------------------------------------------

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidualTVD(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcFluxGalerkinNoBdr1d,&
              euler_calcCharacteristics1d, rproblemLevel&
              %Rafcstab(inviscidAFC), rtimestep%dStep, .true., rres)

        case (NDIM2D)
          call gfsys_buildResidualTVD(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcFluxGalerkinNoBdr2d,&
              euler_calcCharacteristics2d, rproblemLevel&
              %Rafcstab(inviscidAFC), rtimestep%dStep, .true., rres)

        case (NDIM3D)
          call gfsys_buildResidualTVD(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcFluxGalerkinNoBdr3d,&
              euler_calcCharacteristics3d, rproblemLevel&
              %Rafcstab(inviscidAFC), rtimestep%dStep, .true., rres)

        end select

      case DEFAULT
        call output_line('Invalid type of stabilisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_calcResidual')
        call sys_halt()
      end select

      !-------------------------------------------------------------------------
      ! Compute the constant right-hand side for (pseudo-)transient flows
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'imasstype', imasstype)

      select case(imasstype)
      case (MASS_LUMPED)

        !-----------------------------------------------------------------------
        ! Compute the constant right-hand side
        !
        !   $ rhs = M_L*U^n + (1-theta)*res $
        !-----------------------------------------------------------------------
        
        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(1), 'lumpedmassmatrix', lumpedMassMatrix)
        
        call lsysbl_vectorLinearComb(rres, rrhs, 1.0_DP, 0.0_DP)
        
        do iblock = 1, rsolution%nblocks
          call lsyssc_scalarMatVec(rproblemLevel&
              %Rmatrix(lumpedMassMatrix), rsolution&
              %RvectorBlock(iblock), rrhs%RvectorBlock(iblock),&
              1._DP , 1.0_DP-rtimestep%theta)
        end do

      case (MASS_CONSISTENT)

        !-----------------------------------------------------------------------
        ! Compute the constant right-hand side
        !
        !   $ rhs = M_C*U^n + (1-theta)*res $
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(1), 'consistentmassmatrix', consistentMassMatrix)

        call lsysbl_vectorLinearComb(rres, rrhs, 1.0_DP, 0.0_DP)
        
        do iblock = 1, rsolution%nblocks
          call lsyssc_scalarMatVec(rproblemLevel&
              %Rmatrix(consistentMassMatrix), rsolution&
              %RvectorBlock(iblock), rrhs%RvectorBlock(iblock),&
              1._DP , 1.0_DP-rtimestep%theta)
        end do
      end select

    else   ! ite > 0

      !-------------------------------------------------------------------------
      ! In all subsequent nonlinear iterations only the residual vector is
      ! updated, using the right-hand side vector from the first iteration
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'imasstype', imasstype)
      
      select case(imasstype)
      case (MASS_LUMPED)

        !-----------------------------------------------------------------------
        ! Adopt the constant right-hand side cetor
        !
        !   $ res = rhs-M_L*U^{(m)} $
        !-----------------------------------------------------------------------
        
        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(1), 'lumpedmassmatrix', lumpedMassMatrix)

        call lsysbl_vectorLinearComb(rrhs, rres, 1.0_DP, 0.0_DP)
        
        do iblock = 1, rsolution%nblocks
          call lsyssc_scalarMatVec(rproblemLevel&
              %Rmatrix(lumpedMassMatrix), rsolution&
              %RvectorBlock(iblock), rres%RvectorBlock(iblock), &
              -1._DP, 1.0_DP)
        end do

      case (MASS_CONSISTENT)

        !-----------------------------------------------------------------------
        ! Adopt the constant right-hand side vector
        !
        !   $ res = rhs-M_C*U^{(m)} $
        !-----------------------------------------------------------------------
        
        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(1), 'consistentmassmatrix', consistentMassMatrix)

        call lsysbl_vectorLinearComb(rrhs, rres, 1.0_DP, 0.0_DP)
        
        do iblock = 1, rsolution%nblocks
          call lsyssc_scalarMatVec(rproblemLevel&
              %Rmatrix(consistentMassMatrix), rsolution&
              %RvectorBlock(iblock), rres%RvectorBlock(iblock), &
              -1._DP, 1.0_DP)
        end do

      case DEFAULT

        !-----------------------------------------------------------------------
        ! Initialize the residual vector by zeros
        !
        !   $ res = 0 $
        !-----------------------------------------------------------------------

        call lsysbl_clearVector(rres)
      end select

      !-------------------------------------------------------------------------
      ! Update the residual vector
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'coeffmatrix_cx', coeffMatrix_CX)
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'coeffmatrix_cy', coeffMatrix_CY)
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'coeffmatrix_cz', coeffMatrix_CZ)
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'inviscidAFC', inviscidAFC)

      select case(rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation)
        
      case (AFCSTAB_GALERKIN)

        !-----------------------------------------------------------------------
        ! Compute the high-order residual
        !
        !   $ res = res + dt*theta*K(U^{(m)})*U^(m) $
        !-----------------------------------------------------------------------
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcFluxGalerkin1d, rtimestep%theta*rtimestep&
              %dStep, .false., rres)

        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcFluxGalerkin2d, rtimestep%theta*rtimestep&
              %dStep, .false., rres)

        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcFluxGalerkin3d, rtimestep%theta*rtimestep&
              %dStep, .false., rres)

        end select

      case (AFCSTAB_UPWIND,&
            AFCSTAB_FEMFCT_LINEARIZED)

        !-----------------------------------------------------------------------
        ! Compute the low-order residual
        !
        !   $ res = res + dt*theta*L(U^{(m)})*U^(m) $
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(1), 'idissipationtype', idissipationtype)

        select case(idissipationtype)

        case (DISSIPATION_SCALAR)

          ! Assemble divergence operator with scalar dissipation

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                euler_calcFluxScalarDiss1d, rtimestep%theta*rtimestep&
                %dStep, .false., rres)

          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                euler_calcFluxScalarDiss2d, rtimestep%theta*rtimestep&
                %dStep, .false., rres)

          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                euler_calcFluxScalarDiss3d, rtimestep%theta*rtimestep&
                %dStep, .false., rres)

          end select

        case (DISSIPATION_SCALAR_DSPLIT)

          ! Assemble divergence operator with scalar dissipation
          !  adopting dimensional splitting

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                euler_calcFluxScalarDiss1d, rtimestep%theta*rtimestep&
                %dStep, .false., rres)

          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                euler_calcFluxDSplitScalarDiss2d, rtimestep%theta&
                *rtimestep%dStep, .false., rres)

          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                euler_calcFluxDSplitScalarDiss3d, rtimestep%theta&
                *rtimestep%dStep, .false., rres)

          end select

        case (DISSIPATION_TENSOR)

          ! Assemble divergence operator with tensorial dissipation

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                euler_calcFluxTensorDiss1d, rtimestep%theta*rtimestep&
                %dStep, .false., rres)

          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                euler_calcFluxTensorDiss2d, rtimestep%theta*rtimestep&
                %dStep, .false., rres)

          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                euler_calcFluxTensorDiss3d, rtimestep%theta*rtimestep&
                %dStep, .false., rres)

          end select

        case (DISSIPATION_TENSOR_DSPLIT)

          ! Assemble divergence operator with tensorial dissipation
          !  adopting dimensional splitting

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                euler_calcFluxTensorDiss1d, rtimestep%theta*rtimestep&
                %dStep, .false., rres)

          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                euler_calcFluxDSplitTensorDiss2d, rtimestep%theta&
                *rtimestep%dStep, .false., rres)

          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                euler_calcFluxDSplitTensorDiss3d, rtimestep%theta&
                *rtimestep%dStep, .false., rres)

          end select

        case (DISSIPATION_RUSANOV)

          ! Assemble divergence operator with Rusanov flux

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                euler_calcFluxRusanov1d, rtimestep%theta*rtimestep&
                %dStep, .false., rres)

          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                euler_calcFluxRusanov2d, rtimestep%theta*rtimestep&
                %dStep, .false., rres)

          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                euler_calcFluxRusanov3d, rtimestep%theta*rtimestep&
                %dStep, .false., rres)

          end select

        case (DISSIPATION_RUSANOV_DSPLIT)

          ! Assemble divergence operator with Rusanov flux adopting
          !  dimensional splitting

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                euler_calcFluxRusanov1d, rtimestep%theta*rtimestep&
                %dStep, .false., rres)

          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                euler_calcFluxDSplitRusanov2d, rtimestep%theta&
                *rtimestep%dStep, .false., rres)

          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                euler_calcFluxDSplitRusanov3d, rtimestep%theta&
                *rtimestep%dStep, .false., rres)

          end select

        case DEFAULT
          call output_line('Invalid type of dissipation!',&
              OU_CLASS_ERROR,OU_MODE_STD,'euler_calcResidual')
          call sys_halt()
        end select
          
      case (AFCSTAB_FEMTVD)

        !-----------------------------------------------------------------------
        ! Compute the low-order residual + FEM-TVD stabilization
        !
        !   $ res = res + dt*theta*L(U^{(m)})*U^{(m)} + F(U^{(m)}) $
        !-----------------------------------------------------------------------

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidualTVD(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcFluxGalerkinNoBdr1d,&
              euler_calcCharacteristics1d, rproblemLevel&
              %Rafcstab(inviscidAFC),&
              rtimestep%theta*rtimestep%dStep, .false., rres)

        case (NDIM2D)
          call gfsys_buildResidualTVD(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcFluxGalerkinNoBdr2d,&
              euler_calcCharacteristics2d, rproblemLevel&
              %Rafcstab(inviscidAFC),&
              rtimestep%theta*rtimestep%dStep , .false., rres)

        case (NDIM3D)
          call gfsys_buildResidualTVD(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcFluxGalerkinNoBdr3d,&
              euler_calcCharacteristics3d, rproblemLevel&
              %Rafcstab(inviscidAFC),&
              rtimestep%theta*rtimestep%dStep , .false., rres)

        end select
        
      case DEFAULT
        call output_line('Invalid type of stabilisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_calcResidual')
        call sys_halt()
      end select
      
    end if
    
    ! Stop time measurement for residual/rhs evaluation
    call stat_stopTimer(p_rtimer)
    
  end subroutine euler_calcResidual

  !*****************************************************************************
  
!<subroutine>

  subroutine euler_calcRHS(rproblemLevel, rtimestep, rsolver,&
      rsolution, rsolutionInitial, rrhs, istep, rcollection)

!<description>
    ! This subroutine computes the right-hand side vector
    ! used in the explicit Lax-Wendroff time-stepping scheme
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolutionInitial

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
        call lsyssc_scalarMatVec(rproblemLevel&
            %Rmatrix(lumpedMassMatrix), rsolution&
            %RvectorBlock(iblock), rrhs%RvectorBlock(iblock), 1.0_DP, 0.0_DP)
      end do

    case(MASS_CONSISTENT)
      
      !-------------------------------------------------------------------------
      ! Initialize the right-hand side vector
      !
      !  $ M_C*U $
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'consistentmassmatrix', consistentMassMatrix)

      do iblock = 1, rsolution%nblocks
        call lsyssc_scalarMatVec(rproblemLevel&
            %Rmatrix(consistentMassMatrix), rsolution&
            %RvectorBlock(iblock), rrhs%RvectorBlock(iblock), 1.0_DP, 0.0_DP)
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
      !   $ rhs = M*U+weight*(1-theta)*dt*K(U)*U $
      !-------------------------------------------------------------------------
      
      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildResidual(rproblemLevel&
            %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
            euler_calcFluxGalerkin1d, dscale, .false., rrhs)

      case (NDIM2D)
        call gfsys_buildResidual(rproblemLevel&
            %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
            euler_calcFluxGalerkin2d, dscale, .false., rrhs)

      case (NDIM3D)
        call gfsys_buildResidual(rproblemLevel&
            %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
            euler_calcFluxGalerkin3d, dscale, .false., rrhs)

      end select

    case (AFCSTAB_UPWIND)

      !-----------------------------------------------------------------------
      ! Compute the low-order right-hand side
      !
      !   $ rhs = M*U+weight*(1-theta)*dt*L(U)*U $
      !-----------------------------------------------------------------------


      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'idissipationtype', idissipationtype)

      select case(idissipationtype)
        
      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcFluxScalarDiss1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcFluxScalarDiss2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcFluxScalarDiss3d, dscale, .false., rrhs)

        end select

      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcFluxTensorDiss1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcFluxTensorDiss2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcFluxTensorDiss3d, dscale, .false., rrhs)

        end select

      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with Rusanov flux

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              euler_calcFluxRusanov1d, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              euler_calcFluxRusanov2d, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              euler_calcFluxRusanov3d, dscale, .false., rrhs)

        end select

      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'euler_calcRHS')
        call sys_halt()
      end select

    case (AFCSTAB_FEMTVD)

      !-----------------------------------------------------------------------
      ! Compute the low-order right-hand side + FEM-TVD stabilization
      !
      !   $ rhs = M*U+weight*(1-theta)*dt*L(U)*U + F(U) $
      !-----------------------------------------------------------------------

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildResidualTVD(rproblemLevel&
            %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
            euler_calcFluxGalerkinNoBdr1d,&
            euler_calcCharacteristics1d, rproblemLevel&
            %Rafcstab(inviscidAFC), dscale, .false., rrhs)

      case (NDIM2D)
        call gfsys_buildResidualTVD(rproblemLevel&
            %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
            euler_calcFluxGalerkinNoBdr2d,&
            euler_calcCharacteristics2d, rproblemLevel&
            %Rafcstab(inviscidAFC), dscale, .false., rrhs)

      case (NDIM3D)
        call gfsys_buildResidualTVD(rproblemLevel&
            %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
            euler_calcFluxGalerkinNoBdr3d,&
            euler_calcCharacteristics3d, rproblemLevel&
            %Rafcstab(inviscidAFC), dscale, .false., rrhs)

      end select

    case DEFAULT
      call output_line('Invalid type of stabilization!',&
          OU_CLASS_ERROR,OU_MODE_STD,'euler_calcRHS')
      call sys_halt()
    end select

    ! Stop time measurement for global operator
    call stat_stopTimer(p_rtimer)
    
  end subroutine euler_calcRHS

  !*****************************************************************************

!<subroutine>

  subroutine euler_setBoundary(rproblemLevel, rtimestep, rsolver,&
      rsolution, rsolutionInitial, rres, rcollection)

!<description>
    ! This subroutine imposes the nonlinear boundary conditions.
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! nonlinear solver structure
    type(t_solver), intent(in) :: rsolver

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolutionInitial
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
      
      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'systemmatrix', imatrix)
      
    case (NLSOL_PRECOND_NEWTON)
      
      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'jacobianmatrix', imatrix)

    case DEFAULT
      call output_line('Invalid nonlinear preconditioner!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_setBoundary')
      call sys_halt()
    end select
      

    ! Impose boundary conditions for the solution vector and impose
    ! zeros in the residual vector and the off-diagonal positions 
    ! of the system matrix which is obtained from the collection
    
    select case(rproblemLevel%rtriangulation%ndim)
    case (NDIM1D)
      call bdrf_filterSolution(rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixBlock(imatrix), rsolution, rres,&
          rsolutionInitial, rtimestep%dTime,&
          euler_calcBoundaryvalues1d, istatus)

    case (NDIM2D)
      call bdrf_filterSolution(rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixBlock(imatrix), rsolution, rres,&
          rsolutionInitial, rtimestep%dTime,&
          euler_calcBoundaryvalues2d, istatus)

    case (NDIM3D)
      call bdrf_filterSolution(rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixBlock(imatrix), rsolution, rres,&
          rsolutionInitial, rtimestep%dTime,&
          euler_calcBoundaryvalues3d, istatus)

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'euler_setBoundary')
      call sys_halt()
    end select
    
  end subroutine euler_setBoundary

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcLinearizedFCT(rbdrCond, rproblemLevel,&
      rtimestep, rsolution, rcollection)

!<description>
    ! This subroutine calculates the linearized FCT correction
!</description>

!<input>
    ! boundary condition structure
    type(t_boundaryCondition), intent(in) :: rbdrCond

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep
!</input>

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixScalar), pointer :: p_rmatrix
    type(t_parlist), pointer :: p_rparlist
    type(t_vectorScalar) :: rflux0, rflux, ralpha
    type(t_vectorBlock) :: rdata
    real(DP), dimension(:), pointer :: p_MC, p_ML, p_Cx, p_Cy, p_u, p_flux0, p_flux, p_data, p_alpha
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal, p_Ksep
    integer :: h_Ksep, templatematrix, lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, nedge

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'templatematrix', templateMatrix)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'coeffMatrix_CX', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'coeffMatrix_CY', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'lumpedmassmatrix', lumpedMassMatrix)
    
    ! Set pointers to template matrix
    p_rmatrix => rproblemLevel%Rmatrix(templatematrix)
    call lsyssc_getbase_Kld(p_rmatrix, p_Kld)
    call lsyssc_getbase_Kcol(p_rmatrix, p_Kcol)
    call lsyssc_getbase_Kdiagonal(p_rmatrix, p_Kdiagonal)
    
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(consistentMassMatrix), p_MC)
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(lumpedMassMatrix), p_ML)
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(coeffMatrix_CX), p_Cx)
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(coeffMatrix_CY), p_Cy)

    ! Create diagonal separator
    h_Ksep = ST_NOHANDLE
    call storage_copy(p_rmatrix%h_Kld, h_Ksep)
    call storage_getbase_int(h_Ksep, p_Ksep, p_rmatrix%NEQ+1)

    ! Compute number of edges
    nedge = int(0.5*(p_rmatrix%NA-p_rmatrix%NEQ))

    ! Create auxiliary vectors
    call lsyssc_createVector(rflux0, nedge, NVAR2D, .true., ST_DOUBLE)
    call lsyssc_createVector(rflux,  nedge, NVAR2D, .true., ST_DOUBLE)
    call lsyssc_createVector(ralpha, nedge, .false., ST_DOUBLE)
    call lsysbl_createVectorBlock(rsolution, rdata, .true.)
    
    ! Set pointers
    call lsysbl_getbase_double(rsolution, p_u)
    call lsysbl_getbase_double(rdata, p_data)
    call lsyssc_getbase_double(rflux, p_flux)
    call lsyssc_getbase_double(rflux0, p_flux0)
    call lsyssc_getbase_double(ralpha, p_alpha)

    ! Initialize alpha with ones
    p_alpha = 1.0_DP
    
    ! Build the flux
    call buildFluxCons2d(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
        p_rmatrix%NEQ, NVAR2D, nedge, p_u, rtimestep%dStep, p_MC,&
        p_ML, p_Cx, p_Cy, p_data, p_flux, p_flux0)

    ! Build the correction
    call buildCorrectionCons(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
        p_rmatrix%NEQ, NVAR2D, nedge, p_ML, p_flux, p_flux0, 1,&
        p_alpha, p_u)
    call buildCorrectionCons(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
        p_rmatrix%NEQ, NVAR2D, nedge, p_ML, p_flux, p_flux0, 4,&
        p_alpha, p_u)

    ! Apply correction to low-order solution
    call applyCorrection(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
        p_rmatrix%NEQ, NVAR2D, nedge, p_ML, p_flux, p_alpha, p_data,&
        p_u)

    ! Set boundary conditions explicitly
    call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
        rtimestep%dTime, euler_calcBoundaryvalues2d)

    ! Release flux vectors
    call storage_free(h_Ksep)
    call lsyssc_releaseVector(rflux0)
    call lsyssc_releaseVector(rflux)
    call lsyssc_releaseVector(ralpha)
    call lsysbl_releaseVector(rdata)

  contains

    !***************************************************************************

    subroutine buildFluxCons2d(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        NEDGE, u, dscale, MC, ML, Cx, Cy, troc, flux, flux0)

      real(DP), dimension(NVAR,NEQ), intent(in) :: u
      real(DP), dimension(:), intent(in) :: MC,ML,Cx,Cy
      real(DP), intent(in) :: dscale
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NVAR,NEDGE
      
      integer, dimension(:), intent(inout) :: Ksep
      real(DP), dimension(NVAR,NEDGE), intent(inout) :: flux0,flux

      real(DP), dimension(NVAR,NEQ), intent(out) :: troc     

      real(DP), parameter :: c_d = 1.0

      ! local variables
      real(DP), dimension(NVAR) :: K_ij,K_ji,D_ij,Diff,F_ij,F_ji
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      integer :: ij,ji,i,j,iedge

      ! Initialize time rate of change
      call lalg_clearVector(troc)
      
      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1; iedge = iedge+1
          
          ! Compute coefficients
          C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
          C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)
          
          ! Calculate low-order flux
          call euler_calcFluxRusanov2d(u(:,i), u(:,j), C_ij, C_ji, i, j, dscale, F_ij, F_ji)
          
          ! Update the time rate of change vector
          troc(:,i) = troc(:,i) + F_ij
          troc(:,j) = troc(:,j) + F_ji

          ! Calculate diffusion coefficient
          call euler_calcMatrixRusanovDiag2d(u(:,i), u(:,j), C_ij, C_ji, i, j, dscale, K_ij, K_ji, D_ij)
          
          ! Compute solution difference
          Diff = u(:,i)-u(:,j)         

          ! Compute the raw antidiffusive flux
          flux0(:,iedge) = D_ij*Diff
          
        end do
      end do

      ! Scale the time rate of change by the lumped mass matrix
      do i = 1, NEQ
        troc(:,i) = troc(:,i)/ML(i)
      end do

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1

        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1
          
          ! Apply mass antidiffusion
          flux(:,iedge) = flux0(:,iedge) + MC(ij)*(troc(:,i)-troc(:,j))
          
          ! Update edge counter
          iedge = iedge-1
          
        end do
      end do

    end subroutine buildFluxCons2d
    
    !***************************************************************************

    subroutine  buildCorrectionCons(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        NVAR, NEDGE, ML, flux, flux0, ivar, alpha, u)

      real(DP), dimension(NVAR,NEDGE), intent(in) :: flux0
      real(DP), dimension(:), intent(in) :: ML
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NVAR,NEDGE,ivar
      
      real(DP), dimension(NVAR,NEDGE), intent(inout) :: flux
      real(DP), dimension(NVAR,NEQ), intent(inout) :: u
      real(DP), dimension(:), intent(inout) :: alpha
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(:), allocatable :: pp,pm,qp,qm,rp,rm
      real(DP) :: f_ij,f0_ij,diff,aux,p_ij,p_ji,p0_ij,p0_ji,u_i,u_j,v_i,v_j,r_i,r_j
      integer :: ij,ji,i,j,iedge

      allocate(pp(neq), pm(neq), qp(neq), qm(neq), rp(neq), rm(neq))
      
      pp = 0.0_DP; pm = 0.0_DP
      qp = 0.0_DP; qm = 0.0_DP
      rp = 1.0_DP; rm = 1.0_DP

      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; iedge = iedge+1
          
!!$          ! Flux correction in primitive variables
!!$          select case(ivar)
!!$          case default

            ! Flux correction in conservative variables
            f_ij  = flux(ivar,iedge)
            f0_ij = flux0(ivar,iedge)
            diff  = u(ivar,j)-u(ivar,i)
            
            ! MinMod prelimiting of antidiffusive fluxes
            if (f_ij > SYS_EPSREAL .and. f0_ij > SYS_EPSREAL) then
              aux = min(f_ij, 2*f0_ij)
              alpha(iedge) = min(alpha(iedge), aux/f_ij)
              f_ij = aux
            elseif (f_ij < - SYS_EPSREAL .and. f0_ij < -SYS_EPSREAL) then
              aux = max(f_ij, 2*f0_ij)
              alpha(iedge) = min(alpha(iedge), aux/f_ij)
              f_ij = aux
            else
              f_ij = 0.0; alpha(iedge) = 0.0
            end if
            
            ! Sums of raw antidiffusive fluxes
            pp(i) = pp(i) + max(0.0_DP,  f_ij)
            pp(j) = pp(j) + max(0.0_DP, -f_ij)
            pm(i) = pm(i) + min(0.0_DP,  f_ij)
            pm(j) = pm(j) + min(0.0_DP, -f_ij)
            
            ! Sums of admissible edge contributions
            qp(i) = max(qp(i),  diff)
            qp(j) = max(qp(j), -diff)
            qm(i) = min(qm(i),  diff)
            qm(j) = min(qm(j), -diff)
          
            
!!$          case (4)
!!$            ! Velocities
!!$            u_i = u(2,i)/u(1,i);   v_i = u(3,i)/u(1,i)
!!$            u_j = u(2,j)/u(1,j);   v_j = u(3,j)/u(1,j)
!!$
!!$            ! Pressure variables
!!$            p_ij = (GAMMA-1) * (flux(4, iedge) + &
!!$                         0.5 * (u_i*u_i + v_i*v_i)*flux(1, iedge) -&
!!$                                u_i*flux(2, iedge) - v_i*flux(3, iedge) )
!!$            p_ji = -(GAMMA-1) * (flux(4, iedge) + &
!!$                          0.5 * (u_j*u_j + v_j*v_j)*flux(1, iedge) -&
!!$                                 u_j*flux(2, iedge) - v_j*flux(3, iedge) )
!!$            p0_ij = (GAMMA-1) * (flux0(4, iedge) + &
!!$                          0.5 * (u_i*u_i + v_i*v_i)*flux0(1, iedge) -&
!!$                                 u_i*flux0(2, iedge) - v_i*flux0(3, iedge) )
!!$            p0_ji = -(GAMMA-1) * (flux0(4, iedge) + &
!!$                           0.5 * (u_j*u_j + v_j*v_j)*flux0(1, iedge) -&
!!$                                  u_j*flux0(2, iedge) - v_j*flux0(3, iedge) )
!!$
!!$            ! Solution differences
!!$            diff = (GAMMA-1) * (u(4,j) + &
!!$                         0.5 * (u_j*u_j + v_j*v_j)*u(1,j) -&
!!$                                u_j*u(2,j) - v_j*u(3,j) )&
!!$                 - (GAMMA-1) * (u(4,i) + &
!!$                         0.5 * (u_i*u_i + v_i*v_i)*u(1,i) -&
!!$                                u_i*u(2,i) - v_i*u(3,i) )
!!$
!!$            ! Pressure variables
!!$            p_ij =   (GAMMA-1) * ( u(4,i)*flux(1,iedge) - u(2,i)*flux(2,iedge) &
!!$                                  -u(3,i)*flux(3,iedge) + u(1,i)*flux(4,iedge) )
!!$            p_ji =  -(GAMMA-1) * ( u(4,j)*flux(1,iedge) - u(2,j)*flux(2,iedge) &
!!$                                  -u(3,j)*flux(3,iedge) + u(1,j)*flux(4,iedge) )
!!$            p0_ij =  (GAMMA-1) * ( u(4,i)*flux0(1,iedge) - u(2,i)*flux0(2,iedge) &
!!$                                  -u(3,i)*flux0(3,iedge) + u(1,i)*flux0(4,iedge) )
!!$            p0_ji = -(GAMMA-1) * ( u(4,j)*flux0(1,iedge) - u(2,j)*flux0(2,iedge) &
!!$                                  -u(3,j)*flux0(3,iedge) + u(1,j)*flux0(4,iedge) )
!!$
!!$            ! Solution differences
!!$            diff =  (GAMMA-1) * ( u(4,j)*u(1,j) - u(2,j)*u(2,j) &
!!$                                 -u(3,j)*u(3,j) + u(1,j)*u(4,j) )&
!!$                   -(GAMMA-1) * ( u(4,i)*u(1,i) - u(2,i)*u(2,i) &
!!$                                 -u(3,i)*u(3,i) + u(1,i)*u(4,i) )
!!$            
!!$            ! MinMod prelimiting of antidiffusive fluxes
!!$            if ((p_ij >  SYS_EPSREAL .and. p0_ij >  SYS_EPSREAL) .or.&
!!$                (p_ij < -SYS_EPSREAL .and. p0_ij < -SYS_EPSREAL)) then
!!$              aux = min(p_ij, 2*p0_ij)
!!$              alpha(iedge) = min(alpha(iedge), aux/p_ij)
!!$              p_ij = aux
!!$            else
!!$              p_ij = 0.0; alpha(iedge) = 0.0
!!$            end if
!!$
!!$            if ((p_ji >  SYS_EPSREAL .and. p0_ji >  SYS_EPSREAL) .or.&
!!$                (p_ji < -SYS_EPSREAL .and. p0_ji < -SYS_EPSREAL)) then
!!$              aux = min(p_ji, 2*p0_ji)
!!$              alpha(iedge) = min(alpha(iedge), aux/p_ji)
!!$              p_ji = aux
!!$            else
!!$              p_ji = 0.0; alpha(iedge) = 0.0
!!$            end if
!!$            
!!$            ! Sums of raw antidiffusive fluxes
!!$            pp(i) = pp(i) + max(0.0_DP, p_ij)
!!$            pp(j) = pp(j) + max(0.0_DP, p_ji)
!!$            pm(i) = pm(i) + min(0.0_DP, p_ij)
!!$            pm(j) = pm(j) + min(0.0_DP, p_ji)
!!$            
!!$            ! Sums of admissible edge contributions
!!$            qp(i) = max(qp(i),  diff)
!!$            qp(j) = max(qp(j), -diff)
!!$            qm(i) = min(qm(i),  diff)
!!$            qm(j) = min(qm(j), -diff)
!!$          end select
            
        end do
      end do

      ! Compute nodal correction factors
      do i = 1, NEQ
        qp(i) = qp(i)*ML(i)
        qm(i) = qm(i)*ML(i)

        if (pp(i) > qp(i) + SYS_EPSREAL) rp(i) = qp(i)/pp(i)
        if (pm(i) < qm(i) - SYS_EPSREAL) rm(i) = qm(i)/pm(i)
      end do
      
      ! Loop over all rows (backward)
      do i = NEQ, 1, -1
        
        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

!!$          ! Flux correction in primitive variables
!!$          select case(ivar)
!!$          case default
            
            ! Flux correction in conservative variables
            f_ij = flux(ivar,iedge)
            
            ! Limit conservative fluxes
            if (f_ij > SYS_EPSREAL) then
              alpha(iedge) = min(alpha(iedge), rp(i), rm(j))
            elseif (f_ij < -SYS_EPSREAL) then
              alpha(iedge) = min(alpha(iedge), rm(i), rp(j))
            end if

!!$          case (4)
!!$            
!!$            ! Flux correction in primitive variables
!!$            
!!$            ! Velocities
!!$            u_i = u(2,i)/u(1,i);   v_i = u(3,i)/u(1,i)
!!$            u_j = u(2,j)/u(1,j);   v_j = u(3,j)/u(1,j)
!!$
!!$            ! Pressure variables
!!$            p_ij = (GAMMA-1) * (flux(4, iedge) + &
!!$                         0.5 * (u_i*u_i + v_i*v_i)*flux(1, iedge) -&
!!$                                u_i*flux(2, iedge) - v_i*flux(3, iedge) )
!!$            p_ji = -(GAMMA-1) * (flux(4, iedge) + &
!!$                          0.5 * (u_j*u_j + v_j*v_j)*flux(1, iedge) -&
!!$                                 u_j*flux(2, iedge) - v_j*flux(3, iedge) )
!!$
!!$            ! Pressure variables
!!$            p_ij = (GAMMA-1) * ( u(4,i)*flux(1,iedge) - u(2,i)*flux(2,iedge) &
!!$                                -u(3,i)*flux(3,iedge) + u(1,i)*flux(4,iedge) )
!!$            
!!$            p_ji = -(GAMMA-1) * ( u(4,j)*flux(1,iedge) - u(2,j)*flux(2,iedge) &
!!$                                 -u(3,j)*flux(3,iedge) + u(1,j)*flux(4,iedge) )
!!$
!!$            if (p_ij > SYS_EPSREAL) then
!!$              r_i = rp(i)
!!$            elseif (p_ij < -SYS_EPSREAL) then
!!$              r_i = rm(i)
!!$            else
!!$              r_i = 1.0
!!$            end if
!!$
!!$            if (p_ji > SYS_EPSREAL) then
!!$              r_j = rp(j)
!!$            elseif (p_ji < -SYS_EPSREAL) then
!!$              r_j = rm(j)
!!$            else
!!$              r_j = 1.0
!!$            end if
!!$
!!$            alpha(iedge) = min(alpha(iedge), r_i, r_j)
!!$            
!!$          end select
          
          ! Update edge counter
          iedge = iedge-1
          
        end do
      end do

    end subroutine buildCorrectionCons

    !***************************************************************************

    subroutine applyCorrection(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR,&
        NEDGE, ML, flux, alpha, data, u)
      
      real(DP), dimension(NVAR,NEDGE), intent(in) :: flux
      real(DP), dimension(:), intent(in) :: ML,alpha
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NVAR,NEDGE
      
      real(DP), dimension(NVAR,NEQ), intent(inout) :: u,data
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(NVAR) :: f_ij
      integer :: ij,ji,i,j,iedge

      ! Initialize correction
      call lalg_clearVector(data)

      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J. That is, explore the
        ! upper triangular matrix
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; iedge = iedge+1

          ! Limit raw antidiffusive flux
          f_ij = alpha(iedge)*flux(:,iedge)
          
          ! Apply correction
          data(:,i) = data(:,i) + f_ij
          data(:,j) = data(:,j) - f_ij
        end do
      end do

      ! Loop over all rows
      do i = 1, NEQ
        u(:,i) = u(:,i) + data(:,i)/ML(i)
      end do

      ! Just to be sure that this routine can be called repeatedly
      Ksep = Kld
      
    end subroutine applyCorrection

    !***************************************************************************

    pure elemental function minmod(a,b)
      real(DP), intent(in) :: a,b
      real(DP) :: minmod

      if (a > 0 .and. b > 0) then
        minmod = min(a,b)
      elseif (a < 0 .and. b < 0) then
        minmod = max(a,b)
      else
        minmod = 0
      end if
    end function minmod

  end subroutine euler_calcLinearizedFCT

end module euler_callback

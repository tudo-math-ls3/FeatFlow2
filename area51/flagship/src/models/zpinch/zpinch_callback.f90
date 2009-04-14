!##############################################################################
!# ****************************************************************************
!# <name> zpinch_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to solve the 
!# compressible Euler/Navier Stokes equations in arbitrary spatial dimensions.
!#
!# The following callback functions are available:
!#
!# 1.) zpinch_nlsolverCallback
!#     -> Callback routine for the nonlinear solver
!#
!# 2.) zpinch_calcPreconditioner
!#     -> Calculates the nonlinear preconditioner
!#
!# 3.) zpinch_calcJacobian
!#     -> Calculates the Jacobian matrix
!#
!# 4.) zpinch_applyJacobian
!#     -> Applies the Jacobian matrix to a given vector
!#
!# 5.) zpinch_calcResidual
!#     -> Calculates the nonlinear residual vector
!#
!# 6.) zpinch_calcRHS
!#     -> Calculates the right-hand side vector
!#
!# 7.) zpinch_setBoundary
!#     -> Imposes boundary conditions for nonlinear solver
!#
!# 8.) zpinch_calcLinearizedFCT
!#     -> Calculates the linearized FCT correction
!#
!# 9.) zpinch_calcSourceTerm
!#     -> Calculates the source term
!#
!# </purpose>
!##############################################################################

module zpinch_callback

  use afcstabilisation
  use boundaryfilter
  use collection
  use zpinch_basic
  use zpinch_callback1d
  use zpinch_callback2d
  use zpinch_callback3d
  use flagship_basic
  use fsystem
  use genoutput
  use groupfemsystem
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solveraux
  use statistics
  use storage
  use timestepaux

  implicit none

  private
  public :: zpinch_nlsolverCallback
  public :: zpinch_calcPreconditioner
  public :: zpinch_calcJacobian
  public :: zpinch_applyJacobian
  public :: zpinch_calcResidual
  public :: zpinch_calcRHS
  public :: zpinch_setBoundary
  public :: zpinch_calcLinearizedFCT
  public :: zpinch_calcSourceTerm

contains

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_nlsolverCallback(rproblemLevel, rtimestep, rsolver,&
                                    rsolution, rsolutionInitial,&
                                    rrhs, rres, istep, ioperationSpec,&
                                    rcollection, istatus, rb)

!<description>
    ! This subroutine is called by the nonlinear solver and it is responsible
    ! to assemble preconditioner, right-hand side vector, residual vector, etc.
!</description>

!<input>
    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionInitial
    
    ! number of solver step
    integer, intent(IN) :: istep
    
    ! specifier for operations
    integer(I32), intent(IN) :: ioperationSpec

    ! OPTIONAL: constant right-hand side vector
    type(t_vectorBlock), intent(IN), optional :: rb
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
    
    ! right-hand side vector
    type(t_vectorBlock), intent(INOUT) :: rrhs

    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres
        
    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>

!<output>
    ! status flag
    integer, intent(OUT) :: istatus
!</output>
!</subroutine>


    ! Do we have to calculate the preconditioner?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) then
      
      call zpinch_calcPreconditioner(rproblemLevel, rtimestep, rsolver,&
                                     rsolution, rcollection)
    end if
    
    
!!$    ! Do we have to calculate the constant right-hand side?
!!$    ! --------------------------------------------------------------------------
!!$    if ((iand(ioperationSpec, NLSOL_OPSPEC_CALCRHS)  .ne. 0)) then
!!$
!!$      call zpinch_calcrhs(rproblemLevel, rtimestep, rsolver,&
!!$                          rsolution, rvector, rcollection)
!!$    end if

    ! Do we have to calculate the residual
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then
      
      call zpinch_calcResidual(rproblemLevel, rtimestep, rsolver,&
                               rsolution, rsolutionInitial,&
                               rrhs, rres, istep, rcollection)
    end if
    
    
    ! Do we have to impose boundary conditions?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then
      
      call zpinch_setBoundary(rproblemLevel, rtimestep, rsolver,&
                              rsolution, rsolutionInitial, rres, rcollection)
    end if
    
    
    ! Set status flag
    istatus = 0
    
  end subroutine zpinch_nlsolverCallback

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcPreconditioner(rproblemLevel, rtimestep, rsolver,&
                                       rsolution, rcollection)

!<description>
    ! This subroutine calculates the nonlinear preconditioner and
    ! configures the linear solver structure accordingly. Depending on
    ! the nonlinear solver, the low-order evolution operator or the
    ! Jacobian operator is adopted as nonlinear preconditioner.
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(IN) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(IN) :: rsolution
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_timer), pointer :: rtimer
    integer :: systemMatrix, lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: isystemCoupling, isystemPrecond, isystemFormat, imasstype, ivar, jvar

    ! Start time measurement for matrix evaluation
    rtimer => collct_getvalue_timer(rcollection, 'timerAssemblyMatrix')
    call stat_startTimer(rtimer, STAT_TIMERSHORT)

    ! Get parameters from collection which are required unconditionally
    systemMatrix   = collct_getvalue_int(rcollection, 'systemmatrix')
    coeffMatrix_CX = collct_getvalue_int(rcollection, 'coeffMatrix_CX')
    coeffMatrix_CY = collct_getvalue_int(rcollection, 'coeffMatrix_CY')
    coeffMatrix_CZ = collct_getvalue_int(rcollection, 'coeffMatrix_CZ')

    !---------------------------------------------------------------------------
    ! Check if fully explicit time-stepping is used
    !---------------------------------------------------------------------------
    if (rtimestep%theta .le. SYS_EPSREAL) then
      
      isystemFormat        = collct_getvalue_int(rcollection, 'isystemformat')
      imasstype            = collct_getvalue_int(rcollection, 'imasstype')
      lumpedMassMatrix     = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
      consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')

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
                           OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
          call sys_halt()
        end select

      case (SYSTEM_BLOCKFORMAT)
        
        select case(imasstype)
        case (MASS_LUMPED)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, zpinch_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(rproblemLevel%Rmatrix(lumpedMassMatrix), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                .false., .false., .true., .true.)
          end do
          
        case (MASS_CONSISTENT)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, zpinch_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(rproblemLevel%Rmatrix(consistentMassMatrix), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar), 1.0_DP,&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                .false., .false., .true., .true.)
          end do

        case DEFAULT
          call output_line('Empty system matrix is invalid!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
          call sys_halt()
        end select
          
      case DEFAULT
        call output_line('Invalid system format!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
        call sys_halt()
      end select

      ! That's it
      return
    end if

    !---------------------------------------------------------------------------
    ! Assemble divergence operator:  $ -\nabla\cdot\mathbf{F}(U) $
    !
    ! Remark: In future versions this routine will assemble both the primal
    !         and the dual preconditioner depending on the variables
    !         primaldual which has to be extracted from the collection.
    !---------------------------------------------------------------------------

    isystemCoupling = collct_getvalue_int(rcollection, 'isystemCoupling')
    isystemPrecond = collct_getvalue_int(rcollection, 'isystemPrecond')

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
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsolution, zpinch_calcMatrixDiagonalDiag1d,&
                                      zpinch_calcMatrixGalerkinDiag1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, zpinch_calcMatrixDiagonalDiag2d,&
                                      zpinch_calcMatrixGalerkinDiag2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, zpinch_calcMatrixDiagonalDiag3d,&
                                      zpinch_calcMatrixGalerkinDiag3d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
          call sys_halt()
        end select
        
        
      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsolution, zpinch_calcMatrixDiagonalDiag1d,&
                                      zpinch_calcMatrixScalarDissDiag1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, zpinch_calcMatrixDiagonalDiag2d,&
                                      zpinch_calcMatrixScalarDissDiag2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, zpinch_calcMatrixDiagonalDiag3d,&
                                      zpinch_calcMatrixScalarDissDiag3d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_TENSOR)
        
        ! Assemble divergence operator with tensorial dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsolution, zpinch_calcMatrixDiagonalDiag1d,&
                                      zpinch_calcMatrixTensorDissDiag1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, zpinch_calcMatrixDiagonalDiag2d,&
                                      zpinch_calcMatrixTensorDissDiag2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, zpinch_calcMatrixDiagonalDiag3d,&
                                      zpinch_calcMatrixTensorDissDiag3d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)
        
        ! Assemble divergence operator with the Rusanov dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsolution, zpinch_calcMatrixDiagonalDiag1d,&
                                      zpinch_calcMatrixRusanovDiag1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, zpinch_calcMatrixDiagonalDiag2d,&
                                      zpinch_calcMatrixRusanovDiag2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, zpinch_calcMatrixDiagonalDiag3d,&
                                      zpinch_calcMatrixRusanovDiag3d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
          call sys_halt()
        end select


      case DEFAULT
        call output_line('Invalid type of dissipation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
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
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsolution, zpinch_calcMatrixDiagonal1d,&
                                      zpinch_calcMatrixGalerkin1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, zpinch_calcMatrixDiagonal2d,&
                                      zpinch_calcMatrixGalerkin2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, zpinch_calcMatrixDiagonal3d,&
                                      zpinch_calcMatrixGalerkin3d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_SCALAR)
        
        ! Assemble divergence operator with scalar dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsolution, zpinch_calcMatrixDiagonal1d,&
                                      zpinch_calcMatrixScalarDiss1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, zpinch_calcMatrixDiagonal2d,&
                                      zpinch_calcMatrixScalarDiss2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, zpinch_calcMatrixDiagonal3d,&
                                      zpinch_calcMatrixScalarDiss3d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
          call sys_halt()
        end select

        
      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsolution, zpinch_calcMatrixDiagonal1d,&
                                      zpinch_calcMatrixTensorDiss1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, zpinch_calcMatrixDiagonal2d,&
                                      zpinch_calcMatrixTensorDiss2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, zpinch_calcMatrixDiagonal3d,&
                                      zpinch_calcMatrixTensorDiss3d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with the Rusanov dissipation
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsolution, zpinch_calcMatrixDiagonal1d,&
                                      zpinch_calcMatrixRusanov1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, zpinch_calcMatrixDiagonal2d,&
                                      zpinch_calcMatrixRusanov2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, zpinch_calcMatrixDiagonal3d,&
                                      zpinch_calcMatrixRusanov3d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
          call sys_halt()
        end select
          
        
      case DEFAULT
        call output_line('Invalid type of dissipation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
        call sys_halt()
      end select


    case DEFAULT
      call output_line('Invalid type of flow coupling!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
      call sys_halt()
    end select

    
    !---------------------------------------------------------------------------
    ! Assemble the global system operator
    !---------------------------------------------------------------------------

    isystemFormat = collct_getvalue_int(rcollection, 'isystemformat')
    imasstype     = collct_getvalue_int(rcollection, 'imasstype')

    select case(isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)
      
      !-------------------------------------------------------------------------
      ! Assemble global system operator in interleave matrix format
      !-------------------------------------------------------------------------

      select case(imasstype)
      case (MASS_LUMPED)
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_L)-theta*dt*L $

        lumpedMassMatrix = collct_getvalue_int(rcollection, 'lumpedmassmatrix')

        call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(lumpedMassMatrix), 1.0_DP,&
            rproblemLevel%Rmatrix(systemMatrix),&
                                     rtimestep%theta*rtimestep%dStep,&
                                     rproblemLevel%Rmatrix(systemMatrix),&
                                     .false., .false., .true., .true.)

      case (MASS_CONSISTENT)
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_C)-theta*dt*L $
        
        consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')

        call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(consistentMassMatrix), 1.0_DP,&
                                     rproblemLevel%Rmatrix(systemMatrix),&
                                     rtimestep%theta*rtimestep%dStep,&
                                     rproblemLevel%Rmatrix(systemMatrix),&
                                     .false., .false., .true., .true.)
        
      case DEFAULT
        ! Use the global operator for steady-state flow
        !
        !   $ A = -L $
        
      end select
      
      
    case (SYSTEM_BLOCKFORMAT)

      !-------------------------------------------------------------------------
      ! Assemble global system operator in block matrix format
      !-------------------------------------------------------------------------

      select case(imasstype)
      case (MASS_LUMPED)
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_L)-theta*dt*L $

        lumpedMassMatrix = collct_getvalue_int(rcollection, 'lumpedmassmatrix')

        do ivar = 1, zpinch_getNVAR(rproblemLevel)
          do jvar = 1, zpinch_getNVAR(rproblemLevel)

            if (ivar .eq. jvar) then
              call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(lumpedMassMatrix), 1.0_DP,&
                                           rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                                           rtimestep%theta*rtimestep%dStep,&
                                           rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                                           .false., .false., .true., .true.)
            elseif (isystemCoupling .eq. SYSTEM_ALLCOUPLED) then
              call lsyssc_scaleMatrix(rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,jvar),&
                                      rtimestep%theta*rtimestep%dStep)
            end if

          end do
        end do

      
      case (MASS_CONSISTENT)
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_C)-theta*dt*L $

        consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')

        do ivar = 1, zpinch_getNVAR(rproblemLevel)
          do jvar = 1, zpinch_getNVAR(rproblemLevel)

            if (ivar .eq. jvar) then
              call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(consistentMassMatrix), 1.0_DP,&
                                           rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                                           rtimestep%theta*rtimestep%dStep,&
                                           rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                                           .false., .false., .true., .true.)
            elseif (isystemCoupling .eq. SYSTEM_ALLCOUPLED) then
              call lsyssc_scaleMatrix(rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,jvar),&
                                      rtimestep%theta*rtimestep%dStep)
            end if

          end do
        end do


      case DEFAULT
        ! Use the global operator for steady-state flow
        !
        !   $ A = -L $

      end select


    case DEFAULT
      call output_line('Invalid system format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcPreconditioner')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Impose boundary conditions
    !---------------------------------------------------------------------------

    call bdrf_filterMatrix(rsolver%rboundaryCondition, rproblemLevel%rtriangulation,&
                           rproblemLevel%RmatrixBlock(systemMatrix), 1.0_DP)

    ! Ok, we updated the (nonlinear) system operator successfully. Now we still 
    ! have to link it to the solver hierarchy. This is done recursively.
    call flagship_updateSolverMatrix(rproblemLevel, rsolver, systemMatrix,&
                                     isystemFormat, UPDMAT_ALL,&
                                     rproblemLevel%ilev, rproblemLevel%ilev)

    ! Finally, we have to update the content of the solver hierarchy
    call solver_updateContent(rsolver)

    ! Stop time measurement for global operator
    call stat_stopTimer(rtimer)

  end subroutine zpinch_calcPreconditioner

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcJacobian(rproblemLevel, rtimestep, rsolver,&
                                 rsolution, rsolutionInitial, rcollection)

!<description>
    ! This callback subroutine computes the Jacobian matrix for the
    !  compressible Euler/Navier-Stokes equations
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionInitial
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT), target :: rproblemLevel

    ! nonlinear solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! current solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    print *, "!!! The calculation of the Jacobian matrix for the compressible !!!"
    print *, "!!! Euler/Navier Stokes equations has yet not been implemented  !!!"
    stop

  end subroutine zpinch_calcJacobian

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_applyJacobian(rproblemLevel, rx, ry, cx, cy, rcollection)

!<description>
    ! This subroutine applies the (scaled) Jacobian matrix to 
    ! the vector rx and adds the result to the vector ry.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! vector to which the Jacobian should be applied to
    type(t_vectorBlock), intent(IN) :: rx

    ! factor by which rx should be scaled
    real(DP), intent(IN) :: cx

    ! factor by which ry should be scaled
    real(DP), intent(IN) :: cy
!</input>

!<inputoutput>
    ! vector to which the result should be added
    type(t_vectorBlock), intent(INOUT) :: ry

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: jacobianMatrix

    jacobianMatrix = collct_getvalue_int(rcollection, 'jacobianMatrix')

    ! We know where the Jacobian matrix is stored and can apply it
    ! by means of matrix vector multiplication
    call lsysbl_blockMatVec(rproblemLevel%RmatrixBlock(jacobianMatrix), rx, ry, cx, cy)

  end subroutine zpinch_applyJacobian

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcResidual(rproblemLevel, rtimestep, rsolver,&
                                rsolution, rsolutionInitial,&
                                rrhs, rres, ite, rcollection)

!<description>
    ! This subroutine computes the nonlinear residual vector and the
    ! constant right-hand side (only in the first iteration).
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(IN) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionInitial

    ! number of nonlinear iteration
    integer, intent(IN) :: ite
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver
    
    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! right-hand side vector
    type(t_vectorBlock), intent(INOUT) :: rrhs

    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_timer), pointer :: rtimer
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: consistentMassMatrix, lumpedMassMatrix
    integer :: inviscidAFC, imasstype, idissipationtype
    integer :: iblock


    ! Update the global system operator
    call zpinch_calcPreconditioner(rproblemLevel, rtimestep, rsolver, rsolution, rcollection)

    ! Start time measurement for residual/rhs evaluation
    rtimer => collct_getvalue_timer(rcollection, 'timerAssemblyVector')
    call stat_startTimer(rtimer, STAT_TIMERSHORT)

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
      
      coeffMatrix_CX = collct_getvalue_int(rcollection, 'coeffmatrix_cx')
      coeffMatrix_CY = collct_getvalue_int(rcollection, 'coeffmatrix_cy')
      coeffMatrix_CZ = collct_getvalue_int(rcollection, 'coeffmatrix_cz')
      inviscidAFC    = collct_getvalue_int(rcollection, 'inviscidAFC')

      select case(rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation)
        
      case (AFCSTAB_GALERKIN)

        ! Compute the initial high-order residual
        !
        !   $ res = dt*K(U^n)*U^n $

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                   rsolution, zpinch_calcFluxGalerkin1d,&
                                   rtimestep%dStep, .true., rres)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsolution, zpinch_calcFluxGalerkin2d,&
                                   rtimestep%dStep, .true., rres)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsolution, zpinch_calcFluxGalerkin3d,&
                                   rtimestep%dStep, .true., rres)
        end select

        !-----------------------------------------------------------------------
        
      case (AFCSTAB_UPWIND,&
            AFCSTAB_FEMFCT_LINEARIZED)

        ! Compute the initial low-order residual
        !
        !   $ res = dt*L(U^n)*U^n $

        idissipationtype = collct_getvalue_int(rcollection, 'idissipationtype')

        select case(idissipationtype)

        case (DISSIPATION_SCALAR)

          ! Assemble divergence operator with scalar dissipation

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, zpinch_calcFluxScalarDiss1d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, zpinch_calcFluxScalarDiss2d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, zpinch_calcFluxScalarDiss3d,&
                                     rtimestep%dStep, .true., rres)
          end select

        case (DISSIPATION_SCALAR_DSPLIT)

          ! Assemble divergence operator with scalar dissipation adopting dimensional splitting

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, zpinch_calcFluxScalarDiss1d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, zpinch_calcFluxDSplitScalarDiss2d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, zpinch_calcFluxDSplitScalarDiss3d,&
                                     rtimestep%dStep, .true., rres)
          end select

        case (DISSIPATION_TENSOR)

          ! Assemble divergence operator with tensorial dissipation
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, zpinch_calcFluxTensorDiss1d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, zpinch_calcFluxTensorDiss2d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, zpinch_calcFluxTensorDiss3d,&
                                     rtimestep%dStep, .true., rres)
          end select

        case (DISSIPATION_TENSOR_DSPLIT)

          ! Assemble divergence operator with tensorial dissipation adopting dimensional splitting
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, zpinch_calcFluxTensorDiss1d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, zpinch_calcFluxDSplitTensorDiss2d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, zpinch_calcFluxDSplitTensorDiss3d,&
                                     rtimestep%dStep, .true., rres)
          end select

        case (DISSIPATION_RUSANOV)

          ! Assemble divergence operator with Rusanov flux
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, zpinch_calcFluxRusanov1d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, zpinch_calcFluxRusanov2d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, zpinch_calcFluxRusanov3d,&
                                     rtimestep%dStep, .true., rres)
          end select

        case (DISSIPATION_RUSANOV_DSPLIT)

          ! Assemble divergence operator with Rusanov flux adopting dimensional splitting
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, zpinch_calcFluxRusanov1d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, zpinch_calcFluxDSplitRusanov2d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, zpinch_calcFluxDSplitRusanov3d,&
                                     rtimestep%dStep, .true., rres)
          end select

        case DEFAULT
          call output_line('Invalid type of dissipation!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcResidual')
          call sys_halt()
        end select

        !-----------------------------------------------------------------------

      case (AFCSTAB_FEMTVD)

        ! Compute the initial low-order residual + FEM-TVD stabilization
        !
        !   $ res = dt*L(U^n)*U^n + F(U^n) $

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsolution, zpinch_calcFluxGalerkinNoBdr1d,&
                                      zpinch_calcCharacteristics1d,&
                                      rproblemLevel%Rafcstab(inviscidAFC), rtimestep%dStep, .true., rres)
        case (NDIM2D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, zpinch_calcFluxGalerkinNoBdr2d,&
                                      zpinch_calcCharacteristics2d,&
                                      rproblemLevel%Rafcstab(inviscidAFC), rtimestep%dStep, .true., rres)
        case (NDIM3D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, zpinch_calcFluxGalerkinNoBdr3d,&
                                      zpinch_calcCharacteristics3d,&
                                      rproblemLevel%Rafcstab(inviscidAFC), rtimestep%dStep, .true., rres)
        end select

        !-----------------------------------------------------------------------

      case DEFAULT
        call output_line('Invalid type of stabilisation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcResidual')
        call sys_halt()
      end select

      !-------------------------------------------------------------------------
      ! Compute the constant right-hand side for (pseudo-)transient flows
      !-------------------------------------------------------------------------

      imasstype = collct_getvalue_int(rcollection, 'imasstype')

      select case(imasstype)
      case (MASS_LUMPED)

        ! Compute the constant right-hand side
        !
        !   $ rhs = M_L*U^n + (1-theta)*res $

        lumpedMassMatrix = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
        call lsysbl_vectorLinearComb(rres, rrhs, 1.0_DP, 0.0_DP)
        
        do iblock = 1, rsolution%nblocks
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                   rsolution%RvectorBlock(iblock),&
                                   rrhs%RvectorBlock(iblock),&
                                   1._DP, 1.0_DP-rtimestep%theta)
        end do

      case (MASS_CONSISTENT)

        ! Compute the constant right-hand side
        !
        !   $ rhs = M_C*U^n + (1-theta)*res $

        consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
        call lsysbl_vectorLinearComb(rres, rrhs, 1.0_DP, 0.0_DP)
        
        do iblock = 1, rsolution%nblocks
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                   rsolution%RvectorBlock(iblock),&
                                   rrhs%RvectorBlock(iblock),&
                                   1._DP, 1.0_DP-rtimestep%theta)
        end do
      end select

    else   ! ite > 0

      !-------------------------------------------------------------------------
      ! In all subsequent nonlinear iterations only the residual vector is
      ! updated, using the right-hand side vector from the first iteration
      !-------------------------------------------------------------------------

      imasstype = collct_getvalue_int(rcollection, 'imasstype')
      
      select case(imasstype)
      case (MASS_LUMPED)

        ! Adopt the constant right-hand side cetor
        !
        !   $ res = rhs-M_L*U^{(m)} $
        
        lumpedMassMatrix = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
        call lsysbl_vectorLinearComb(rrhs, rres, 1.0_DP, 0.0_DP)
        
        do iblock = 1, rsolution%nblocks
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                   rsolution%RvectorBlock(iblock),&
                                   rres%RvectorBlock(iblock), -1._DP, 1.0_DP)
        end do

      case (MASS_CONSISTENT)

        ! Adopt the constant right-hand side vector
        !
        !   $ res = rhs-M_C*U^{(m)} $
        
        consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
        call lsysbl_vectorLinearComb(rrhs, rres, 1.0_DP, 0.0_DP)
        
        do iblock = 1, rsolution%nblocks
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                   rsolution%RvectorBlock(iblock),&
                                   rres%RvectorBlock(iblock), -1._DP, 1.0_DP)
        end do

      case DEFAULT

        ! Initialize the residual vector by zeros
        !
        !   $ res = 0 $

        call lsysbl_clearVector(rres)
      end select

      !-------------------------------------------------------------------------
      ! Update the residual vector
      !-------------------------------------------------------------------------

      coeffMatrix_CX = collct_getvalue_int(rcollection, 'coeffmatrix_cx')
      coeffMatrix_CY = collct_getvalue_int(rcollection, 'coeffmatrix_cy')
      coeffMatrix_CZ = collct_getvalue_int(rcollection, 'coeffmatrix_cz')
      inviscidAFC    = collct_getvalue_int(rcollection, 'inviscidAFC')

      select case(rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation)
        
      case (AFCSTAB_GALERKIN)
        
        ! Compute the high-order residual
        !
        !   $ res = res + dt*theta*K(U^{(m)})*U^(m) $
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                   rsolution, zpinch_calcFluxGalerkin1d,&
                                   rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsolution, zpinch_calcFluxGalerkin2d,&
                                   rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsolution, zpinch_calcFluxGalerkin3d,&
                                   rtimestep%theta*rtimestep%dStep, .false., rres)
        end select

        !-----------------------------------------------------------------------

      case (AFCSTAB_UPWIND,&
            AFCSTAB_FEMFCT_LINEARIZED)

        ! Compute the low-order residual
        !
        !   $ res = res + dt*theta*L(U^{(m)})*U^(m) $

        idissipationtype = collct_getvalue_int(rcollection, 'idissipationtype')

        select case(idissipationtype)

        case (DISSIPATION_SCALAR)

          ! Assemble divergence operator with scalar dissipation

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, zpinch_calcFluxScalarDiss1d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, zpinch_calcFluxScalarDiss2d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, zpinch_calcFluxScalarDiss3d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          end select

        case (DISSIPATION_SCALAR_DSPLIT)

          ! Assemble divergence operator with scalar dissipation adopting dimensional splitting

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, zpinch_calcFluxScalarDiss1d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, zpinch_calcFluxDSplitScalarDiss2d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, zpinch_calcFluxDSplitScalarDiss3d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          end select

        case (DISSIPATION_TENSOR)

          ! Assemble divergence operator with tensorial dissipation

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, zpinch_calcFluxTensorDiss1d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, zpinch_calcFluxTensorDiss2d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, zpinch_calcFluxTensorDiss3d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          end select

        case (DISSIPATION_TENSOR_DSPLIT)

          ! Assemble divergence operator with tensorial dissipation adopting dimensional splitting

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, zpinch_calcFluxTensorDiss1d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, zpinch_calcFluxDSplitTensorDiss2d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, zpinch_calcFluxDSplitTensorDiss3d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          end select

        case (DISSIPATION_RUSANOV)

          ! Assemble divergence operator with Rusanov flux

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, zpinch_calcFluxRusanov1d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, zpinch_calcFluxRusanov2d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, zpinch_calcFluxRusanov3d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          end select

        case (DISSIPATION_RUSANOV_DSPLIT)

          ! Assemble divergence operator with Rusanov flux adopting dimensional splitting

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, zpinch_calcFluxRusanov1d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, zpinch_calcFluxDSplitRusanov2d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, zpinch_calcFluxDSplitRusanov3d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          end select

        case DEFAULT
          call output_line('Invalid type of dissipation!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcResidual')
          call sys_halt()
        end select

        !-----------------------------------------------------------------------
          
      case (AFCSTAB_FEMTVD)

        ! Compute the low-order residual + FEM-TVD stabilization
        !
        !   $ res = res + dt*theta*L(U^{(m)})*U^{(m)} + F(U^{(m)}) $

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsolution, zpinch_calcFluxGalerkinNoBdr1d,&
                                      zpinch_calcCharacteristics1d,&
                                      rproblemLevel%Rafcstab(inviscidAFC),&
                                      rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM2D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, zpinch_calcFluxGalerkinNoBdr2d,&
                                      zpinch_calcCharacteristics2d,&
                                      rproblemLevel%Rafcstab(inviscidAFC),&
                                      rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM3D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, zpinch_calcFluxGalerkinNoBdr3d,&
                                      zpinch_calcCharacteristics3d,&
                                      rproblemLevel%Rafcstab(inviscidAFC),&
                                      rtimestep%theta*rtimestep%dStep, .false., rres)
        end select

        !-----------------------------------------------------------------------

      case DEFAULT
        call output_line('Invalid type of stabilisation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcResidual')
        call sys_halt()
      end select

    end if
        
  end subroutine zpinch_calcResidual

  !*****************************************************************************
  
!<subroutine>

  subroutine zpinch_calcRHS(rproblemLevel, rtimestep, rsolver,&
                            rsolution, rsolutionInitial,&
                            rrhs, istep, rcollection)

!<description>
    ! This subroutine computes the right-hand side vector
    ! used in the explicit Lax-Wendroff time-stepping scheme
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(IN) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionInitial

    ! number of explicit step
    integer, intent(IN) :: istep
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver
    
    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! right-hand side vector
    type(t_vectorBlock), intent(INOUT) :: rrhs

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_timer), pointer :: rtimer
    real(DP) :: dscale
    integer :: lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: imasstype, inviscidAFC, idissipationtype
    integer :: iblock

    ! Start time measurement for residual/rhs evaluation
    rtimer => collct_getvalue_timer(rcollection, 'timerAssemblyVector')
    call stat_startTimer(rtimer, STAT_TIMERSHORT)

    ! Get parameters from collection which are required unconditionally
    coeffMatrix_CX = collct_getvalue_int(rcollection, 'coeffMatrix_CX')
    coeffMatrix_CY = collct_getvalue_int(rcollection, 'coeffMatrix_CY')
    coeffMatrix_CZ = collct_getvalue_int(rcollection, 'coeffMatrix_CZ')
    inviscidAFC    = collct_getvalue_int(rcollection, 'inviscidAFC')


    !---------------------------------------------------------------------------
    ! Initialize the right-hand side vector
    !---------------------------------------------------------------------------
    imasstype = collct_getvalue_int(rcollection, 'imasstype')

    select case(imasstype)
    case (MASS_LUMPED)
      
      ! Initialize the right-hand side vector
      !
      !  $ M_L*U

      lumpedMassMatrix = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
      do iblock = 1, rsolution%nblocks
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                 rsolution%RvectorBlock(iblock),&
                                 rrhs%RvectorBlock(iblock), 1.0_DP, 0.0_DP)
      end do

    case(MASS_CONSISTENT)
      
      ! Initialize the right-hand side vector
      !
      !  $ M_C*U

      consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
      do iblock = 1, rsolution%nblocks
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(consistentMassMatrix),&
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
    
    inviscidAFC = collct_getvalue_int(rcollection, 'inviscidAFC')

    ! Compute the scaling parameter
    !
    !   $ weight*(1-theta)*dt

    dscale = rtimestep%DmultistepWeights(istep)*(1.0_DP-rtimestep%theta)*rtimestep%dStep

    select case(rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation)      
      
    case (AFCSTAB_GALERKIN)

      ! Compute the high-order right-hand side
      !
      !   $ rhs = M*U+weight*(1-theta)*dt*K(U)*U
      
      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                 rsolution, zpinch_calcFluxGalerkin1d, dscale, .false., rrhs)
      case (NDIM2D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                 rsolution, zpinch_calcFluxGalerkin2d, dscale, .false., rrhs)
      case (NDIM3D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                 rsolution, zpinch_calcFluxGalerkin3d, dscale, .false., rrhs)
      end select

      !-----------------------------------------------------------------------

    case (AFCSTAB_UPWIND)

      ! Compute the low-order right-hand side
      !
      !   $ rhs = M*U+weight*(1-theta)*dt*L(U)*U
      
      idissipationtype = collct_getvalue_int(rcollection, 'idissipationtype')

      select case(idissipationtype)
        
      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                   rsolution, zpinch_calcFluxScalarDiss1d, dscale, .false., rrhs)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsolution, zpinch_calcFluxScalarDiss2d, dscale, .false., rrhs)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsolution, zpinch_calcFluxScalarDiss3d, dscale, .false., rrhs)
        end select

      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                   rsolution, zpinch_calcFluxTensorDiss1d, dscale, .false., rrhs)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsolution, zpinch_calcFluxTensorDiss2d, dscale, .false., rrhs)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsolution, zpinch_calcFluxTensorDiss3d, dscale, .false., rrhs)
        end select

      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with Rusanov flux

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                   rsolution, zpinch_calcFluxRusanov1d, dscale, .false., rrhs)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsolution, zpinch_calcFluxRusanov2d, dscale, .false., rrhs)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsolution, zpinch_calcFluxRusanov3d, dscale, .false., rrhs)
        end select

      case DEFAULT
        call output_line('Invalid type of dissipation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcRHS')
        call sys_halt()
      end select

      !-----------------------------------------------------------------------

    case (AFCSTAB_FEMTVD)

      ! Compute the low-order right-hand side + FEM-TVD stabilization
      !
      !   $ rhs = M*U+weight*(1-theta)*dt*L(U)*U + F(U)

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                    rsolution, zpinch_calcFluxGalerkinNoBdr1d,&
                                    zpinch_calcCharacteristics1d,&
                                    rproblemLevel%Rafcstab(inviscidAFC), dscale, .false., rrhs)
      case (NDIM2D)
        call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                    rsolution, zpinch_calcFluxGalerkinNoBdr2d,&
                                    zpinch_calcCharacteristics2d,&
                                    rproblemLevel%Rafcstab(inviscidAFC), dscale, .false., rrhs)
      case (NDIM3D)
        call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                    rsolution, zpinch_calcFluxGalerkinNoBdr3d,&
                                    zpinch_calcCharacteristics3d,&
                                    rproblemLevel%Rafcstab(inviscidAFC), dscale, .false., rrhs)
      end select

      !-----------------------------------------------------------------------

    case DEFAULT
      call output_line('Invalid type of stabilization!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcRHS')
      call sys_halt()
    end select

    ! Stop time measurement for global operator
    call stat_stopTimer(rtimer)
    
  end subroutine zpinch_calcRHS

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_setBoundary(rproblemLevel, rtimestep, rsolver,&
                                rsolution, rsolutionInitial, rres, rcollection)

!<description>
    ! This subroutine imposes the nonlinear boundary conditions.
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! nonlinear solver structure
    type(t_solver), intent(IN) :: rsolver

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionInitial
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! collection structure to provide additional
    ! information to the boundary setting routine
    type(t_collection), intent(InOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: imatrix, istatus
    

    ! What type of nonlinear preconditioner are we?
    select case(rsolver%iprecond)
    case (NLSOL_PRECOND_BLOCKD,&
          NLSOL_PRECOND_DEFCOR,&
          NLSOL_PRECOND_NEWTON_FAILED)
      
      imatrix = collct_getvalue_int(rcollection, 'SystemMatrix')
      
    case (NLSOL_PRECOND_NEWTON)
      
      imatrix = collct_getvalue_int(rcollection, 'JacobianMatrix')

    case DEFAULT
      call output_line('Invalid nonlinear preconditioner!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'zpinch_setBoundary')
      call sys_halt()
    end select
      

    ! Impose boundary conditions for the solution vector and impose
    ! zeros in the residual vector and the off-diagonal positions 
    ! of the system matrix which is obtained from the collection
    
    select case(rproblemLevel%rtriangulation%ndim)
    case (NDIM1D)
      call bdrf_filterSolution(rsolver%rboundaryCondition,&
                               rproblemLevel%rtriangulation,&
                               rproblemLevel%RmatrixBlock(imatrix),&
                               rsolution, rres, rsolutionInitial,&
                               rtimestep%dTime,&
                               rproblemLevel%p_rproblem%rboundary,&
                               zpinch_calcBoundaryvalues1d, istatus)
    case (NDIM2D)
      call bdrf_filterSolution(rsolver%rboundaryCondition,&
                               rproblemLevel%rtriangulation,&
                               rproblemLevel%RmatrixBlock(imatrix),&
                               rsolution, rres, rsolutionInitial,&
                               rtimestep%dTime,&
                               rproblemLevel%p_rproblem%rboundary,&
                               zpinch_calcBoundaryvalues2d, istatus)
    case (NDIM3D)
      call bdrf_filterSolution(rsolver%rboundaryCondition,&
                               rproblemLevel%rtriangulation,&
                               rproblemLevel%RmatrixBlock(imatrix),&
                               rsolution, rres, rsolutionInitial,&
                               rtimestep%dTime,&
                               rproblemLevel%p_rproblem%rboundary,&
                               zpinch_calcBoundaryvalues3d, istatus)
    case DEFAULT
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'zpinch_setBoundary')
      call sys_halt()
    end select
    
  end subroutine zpinch_setBoundary

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcLinearizedFCT(rbdrCond, rproblemLevel, rtimestep, rsolution, rcollection)

!<description>
    ! This subroutine calculates the linearized FCT correction
!</description>

!<input>
    ! boundary condition structure
    type(t_boundaryCondition), intent(IN) :: rbdrCond

    ! problem level structure
    type(t_problemLevel), intent(IN) :: rproblemLevel

    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep
!</input>

!<inputoutput>
    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixScalar), pointer :: p_rmatrix
    type(t_vectorScalar) :: rflux0, rflux, ralpha
    type(t_vectorBlock) :: rdata
    real(DP), dimension(:), pointer :: p_MC, p_ML, p_Cx, p_Cy, p_u, p_flux0, p_flux, p_data, p_alpha
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal, p_Ksep
    integer :: h_Ksep, templatematrix, lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, nedge

    templateMatrix       = collct_getvalue_int(rcollection, 'templatematrix')
    consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentMassMatrix')
    lumpedMassMatrix     = collct_getvalue_int(rcollection, 'lumpedMassMatrix')
    coeffMatrix_CX       = collct_getvalue_int(rcollection, 'coeffMatrix_CX')
    coeffMatrix_CY       = collct_getvalue_int(rcollection, 'coeffMatrix_CY')
    
    p_rmatrix => rproblemLevel%Rmatrix(templatematrix)

    ! Set pointers
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

    p_alpha = 1.0_DP
        
    ! Build the flux
    call buildFluxCons2d(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                         NVAR2D, nedge, p_u, rtimestep%dStep,&
                         p_MC, p_ML, p_Cx, p_Cy, p_data, p_flux, p_flux0)

    ! Build the correction
    call buildCorrectionCons(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                             NVAR2D, nedge, p_ML, p_flux, p_flux0, 1, p_alpha, p_u)
    call buildCorrectionCons(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                             NVAR2D, nedge, p_ML, p_flux, p_flux0, 4, p_alpha, p_u)
    call buildCorrectionCons(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                             NVAR2D, nedge, p_ML, p_flux, p_flux0, 5, p_alpha, p_u)

    ! Apply correction to low-order solution
    call applyCorrection(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                         NVAR2D, nedge, p_ML, p_flux, p_alpha, p_data, p_u)

    ! Set boundary conditions explicitly
    call bdrf_filterVectorExplicit(rbdrCond, rproblemLevel%rtriangulation,&
                                   rsolution, rtimestep%dTime,&
                                   rproblemLevel%p_rproblem%rboundary,&
                                   zpinch_calcBoundaryvalues2d)

    ! Release flux vectors
    call storage_free(h_Ksep)
    call lsyssc_releaseVector(rflux0)
    call lsyssc_releaseVector(rflux)
    call lsyssc_releaseVector(ralpha)
    call lsysbl_releaseVector(rdata)

  contains
    
    !***************************************************************************
    ! Build the raw antidiffusive fluxes

    subroutine buildFluxCons2d(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, NEDGE, u,&
                               dscale, MC, ML, Cx, Cy, troc, flux, flux0)

      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: MC,ML,Cx,Cy
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR,NEDGE
      
      integer, dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(NVAR,NEDGE), intent(INOUT) :: flux0,flux

      real(DP), dimension(NVAR,NEQ), intent(OUT) :: troc     

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
          call zpinch_calcFluxRusanov2d(u(:,i), u(:,j), C_ij, C_ji, dscale, F_ij, F_ji)
          
          ! Update the time rate of change vector
          troc(:,i) = troc(:,i) + F_ij
          troc(:,j) = troc(:,j) + F_ji

          ! Calculate diffusion coefficient
          call zpinch_calcMatrixRusanovDiag2d(u(:,i), u(:,j), C_ij, C_ji, dscale, K_ij, K_ji, D_ij)
          
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

    subroutine  buildCorrectionCons(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, NEDGE,&
                                    ML, flux, flux0, ivar, alpha, u)

      real(DP), dimension(NVAR,NEDGE), intent(IN) :: flux0
      real(DP), dimension(:), intent(IN) :: ML
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR,NEDGE,ivar
      
      real(DP), dimension(NVAR,NEDGE), intent(INOUT) :: flux
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: u
      real(DP), dimension(:), intent(INOUT) :: alpha
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(:), allocatable :: pp,pm,qp,qm,rp,rm
      real(DP) :: f_ij,f0_ij,diff,aux,u_ij,v_ij,p
      integer :: ij,ji,i,j,iedge

      allocate(pp(neq), pm(neq), qp(neq), qm(neq), rp(neq), rm(neq))
      
      pp = 0.0_DP; pm = 0.0_DP
      qp = 0.0_DP; qm = 0.0_DP
      rp = 0.0_DP; rm = 0.0_DP

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
          
          if (ivar .eq. 4) then

            u_ij = 0.5_DP * (u(2,i)/u(1,i)+u(2,j)/u(1,j))
            v_ij = 0.5_DP * (u(3,i)/u(1,i)+u(3,j)/u(1,j))

            f_ij = (GAMMA-1) * (flux(4, iedge) + &
                                0.5*(u_ij*u_ij + v_ij*v_ij)*flux(1, iedge) -&
                                u_ij*flux(2, iedge) - v_ij*flux(3, iedge))

            diff = (GAMMA-1) * ((u(4,j)-u(4,i)) + &
                                0.5*(u_ij*u_ij + v_ij*v_ij)*(u(1,j)-u(1,i)) -&
                                u_ij*(u(2,j)-u(2,i)) - v_ij*(u(3,j)-u(3,i)))

          else
            
            f_ij = flux(ivar,iedge)
            diff = u(ivar,j)-u(ivar,i)
            
          end if

          ! Ignore fluxes which are too small
          if (abs(f_ij) .le. sqrt(SYS_EPSREAL)) f_ij = 0.0_DP
          if (abs(diff) .le. sqrt(SYS_EPSREAL)) diff = 0.0_DP

          ! Prelimit antidiffusive fluxes
          if (f_ij * diff < -SYS_EPSREAL) then
            alpha(iedge) = 0.0_DP; f_ij = 0.0_DP
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

        end do
      end do

      ! Compute nodal correction factors
      do i = 1, NEQ
        if (pp(i) >  sqrt(SYS_EPSREAL)) rp(i) = min(1.0_DP, ML(i)*qp(i)/pp(i))
        if (pm(i) < -sqrt(SYS_EPSREAL)) rm(i) = min(1.0_DP, ML(i)*qm(i)/pm(i))
      end do
      
      ! Loop over all rows (backward)
      do i = NEQ, 1, -1
        
        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          if (ivar .eq. 4) then

            u_ij = 0.5_DP * (u(2,i)/u(1,i)+u(2,j)/u(1,j))
            v_ij = 0.5_DP * (u(3,i)/u(1,i)+u(3,j)/u(1,j))

            f_ij = (GAMMA-1) * (flux(4, iedge) + &
                                0.5*(u_ij*u_ij + v_ij*v_ij)*flux(1, iedge) -&
                                u_ij*flux(2, iedge) - v_ij*flux(3, iedge))
            
          else
            
            f_ij = flux(ivar,iedge)
            
          end if

          ! Limit conservative fluxes
          if (f_ij > sqrt(SYS_EPSREAL)) then
            alpha(iedge) = min(alpha(iedge),rp(i), rm(j))
          elseif (f_ij < -sqrt(SYS_EPSREAL)) then
            alpha(iedge) = min(alpha(iedge),rm(i), rp(j))
          end if
          
          ! Update edge counter
          iedge = iedge-1
          
        end do
      end do

    end subroutine buildCorrectionCons
    
    !***************************************************************************

    subroutine applyCorrection(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, NEDGE,&
                               ML, flux, alpha, data, u)

      
      real(DP), dimension(NVAR,NEDGE), intent(IN) :: flux
      real(DP), dimension(:), intent(IN) :: ML,alpha
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR,NEDGE
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: u,data
      integer, dimension(:), intent(INOUT) :: Ksep
      
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
          where (abs(f_ij) .le. SYS_EPSREAL) f_ij = 0.0
          
          ! Apply correction
          data(:,i) = data(:,i) + f_ij
          data(:,j) = data(:,j) - f_ij
        end do
      end do

      do i = 1, NEQ
        u(:,i) = u(:,i) + data(:,i)/ML(i)
      end do
      
    end subroutine applyCorrection

    !***************************************************************************

    pure elemental function minmod(a,b)
      real(DP), intent(IN) :: a,b
      real(DP) :: minmod

      if (a > 0 .and. b > 0) then
        minmod = min(a,b)
      elseif (a < 0 .and. b < 0) then
        minmod = max(a,b)
      else
        minmod = 0
      end if
    end function minmod

  end subroutine zpinch_calcLinearizedFCT

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcSourceTerm(rproblemLevel, rtimestep, rsolution, rcollection)

!<description>
    ! Calculates the source term
!</description>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel
    
    ! time-stepping structure
    type(t_timestep), intent(INOUT) :: rtimestep
        
    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution
    
    ! collection structure
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>


    ! local variable
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata, p_MC, p_ML
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    integer :: neq, nvar, isystemFormat, consistentMassMatrix, lumpedMassMatrix
    
    ! Get lumped and consistent mass matrix
    lumpedMassMatrix     = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
    consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(lumpedMassMatrix), p_ML)
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(consistentMassMatrix), p_MC)

    call lsyssc_getbase_Kld(rproblemLevel%Rmatrix(consistentMassMatrix), p_Kld)
    call lsyssc_getbase_Kcol(rproblemLevel%Rmatrix(consistentMassMatrix), p_Kcol)
    
    ! Set pointer to global solution vector
    call lsysbl_getbase_double(rsolution, p_Ddata)
    
    ! Set pointer to the vertex coordinates
    call storage_getbase_double2D(&
        rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)
    
    ! Set dimensions
    nvar = zpinch_getNVAR(rproblemLevel)
    neq  = rsolution%NEQ/nvar
        
    isystemFormat = collct_getvalue_int(rcollection, 'isystemformat')
    select case(isystemFormat)
    case (SYSTEM_INTERLEAVEFORMAT)
      call calcSourceTermInterleaveFormat(rtimestep%dTime, rtimestep%dStep, neq, nvar,&
                                          p_DvertexCoords, p_Kld, p_Kcol, p_MC, p_ML, p_Ddata)
    case DEFAULT
      call output_line('Invalid system format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcSourceTerm')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the real working routines follow
    
    !**************************************************************
    
    subroutine calcSourceTermInterleaveFormat(dtime, dstep, neq, nvar, DvertexCoords,&
                                              Kld, Kcol, MC, ML, Ddata)

      real(DP), dimension(:,:), intent(IN) :: DvertexCoords
      real(DP), dimension(:), intent(IN) :: MC, ML
      real(DP), intent(IN) :: dtime, dstep
      integer, dimension(:), intent(IN) :: Kld, Kcol
      integer, intent(IN) :: neq, nvar
      
      real(DP), dimension(nvar,neq), intent(INOUT) :: Ddata
      
      ! local variables
      real(DP), dimension(:,:), allocatable :: DsourceTerm
      real(DP) :: drad, dang, daux, dscale, x1, x2, p, rq
      integer :: i,j,ij


      ! Compute the scaling parameter
      dscale = -dstep * 12.0 * (1.0-dtime**4) * dtime**2
!!$      dscale = -dstep * 12.0 * dtime*dtime

      allocate(DsourceTerm(2,neq)); DsourceTerm = 0.0_DP

      ! Loop over all rows
      do i = 1, neq

        ! Loop over all columns
        do ij = Kld(i), Kld(i+1)-1

          ! Get columns number
          j = Kcol(ij)
        
          ! Get coodrinates at node j
          x1 = DvertexCoords(1, j)
          x2 = DvertexCoords(2, j)
          
          ! Compute polar coordinates
          drad = sqrt(x1*x1 + x2*x2)
          dang = atan2(x2, x1)
          
          ! Compute unit vector into origin
          x1 = cos(dang)
          x2 = sin(dang)
          
          ! Compute source term
          if (Ddata(nvar,j) > sqrt(SYS_EPSREAL)) then
            daux = dscale * MC(ij) * Ddata(nvar,j) / max(drad, 1.0e-4_DP)
          else
            daux = 0.0_DP
          end if
                    
          ! Impose source values
          DsourceTerm(1, i) = DsourceTerm(1, i) + daux * x1
          DsourceTerm(2, i) = DsourceTerm(2, i) + daux * x2

        end do
      end do

      do i = 1, neq

        ! Compute kinetic energy from momentum values without source term
        rq = 0.5 * ( Ddata(2, i)*Ddata(2, i) +&
                     Ddata(3, i)*Ddata(3, i) ) / Ddata(1, i)

        ! Compute pressure value
        p = Ddata(4, i) - rq

        ! Update momentum equations
        Ddata(2, i) = Ddata(2, i) + DsourceTerm(1, i)/ML(i)
        Ddata(3, i) = Ddata(3, i) + DsourceTerm(2, i)/ML(i)

        ! Compute kinetic energy from momentum values with source term
        rq = 0.5 * ( Ddata(2, i)*Ddata(2, i) +&
                     Ddata(3, i)*Ddata(3, i) ) / Ddata(1, i)

        ! Update total energy equation
        Ddata(4, i) = p + rq

      end do
      
      deallocate(DsourceTerm)

    end subroutine calcSourceTermInterleaveFormat

  end subroutine zpinch_calcSourceTerm

end module zpinch_callback

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
!# 1.) euler_calcPreconditioner
!#     -> Calculates the nonlinear preconditioner
!#
!# 2.) euler_calcJacobian
!#     -> Calculates the Jacobian matrix
!#
!# 3.) euler_applyJacobian
!#     -> Applies the Jacobian matrix to a given vector
!#
!# 4.) euler_calcResidual
!#     -> Calculates the nonlinear residual vector
!#
!# 5.) euler_calcRHS
!#     -> Calculates the right-hand side vector
!#
!# </purpose>
!##############################################################################

module euler_callback

  use afcstabilisation
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
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solver
  use statistics
  use storage

  implicit none

  private
  public :: euler_calcPreconditioner
  public :: euler_calcJacobian
  public :: euler_applyJacobian
  public :: euler_calcResidual
  public :: euler_calcRHS

contains

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcPreconditioner(rproblemLevel, rtimestep, rsolver,&
                                      rsol, rcollection)

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
    type(t_vectorBlock), intent(IN) :: rsol
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
    integer :: isystemCoupling, iprecond, isystemFormat, imasstype, ivar, jvar

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
    iprecond = collct_getvalue_int(rcollection, 'iprecond')

    ! What kind of coupling is applied?
    select case(isystemCoupling)
      
    case (SYSTEM_SEGREGATED)
      
      !-------------------------------------------------------------------------
      ! Assemble block-diagonal divergence operator
      !-------------------------------------------------------------------------
      
      ! What kind of preconditioner is applied?
      select case(iprecond)

      case (AFCSTAB_GALERKIN)

        ! Assemble divergence operator for standard Galerin scheme

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsol, euler_calcMatrixGalerkinDiag1d, -1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsol, euler_calcMatrixGalerkinDiag2d, -1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsol, euler_calcMatrixGalerkinDiag3d, -1.0_DP, .true.,&
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
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsol, euler_calcMatrixScalarDissDiag1d, -1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsol, euler_calcMatrixScalarDissDiag2d, -1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsol, euler_calcMatrixScalarDissDiag3d, -1.0_DP, .true.,&
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
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsol, euler_calcMatrixTensorDissDiag1d, -1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsol, euler_calcMatrixTensorDissDiag2d, -1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsol, euler_calcMatrixTensorDissDiag3d, -1.0_DP, .true.,&
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
      select case(iprecond)
        
      case (AFCSTAB_GALERKIN)
        
        ! Assemble divergence operator for standard Galerin scheme
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsol, euler_calcMatrixGalerkin1d, -1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsol, euler_calcMatrixGalerkin2d, -1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsol, euler_calcMatrixGalerkin3d, -1.0_DP, .true.,&
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
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsol, euler_calcMatrixScalarDiss1d, -1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsol, euler_calcMatrixScalarDiss2d, -1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsol, euler_calcMatrixScalarDiss3d, -1.0_DP, .true.,&
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
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsol, euler_calcMatrixTensorDiss1d, -1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsol, euler_calcMatrixTensorDiss2d, -1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsol, euler_calcMatrixTensorDiss3d, -1.0_DP, .true.,&
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


    case DEFAULT
      call output_line('Invalid type of flow coupling!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
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

        do ivar = 1, euler_getNVAR(rproblemLevel)
          do jvar = 1, euler_getNVAR(rproblemLevel)

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

        do ivar = 1, euler_getNVAR(rproblemLevel)
          do jvar = 1, euler_getNVAR(rproblemLevel)

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
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPreconditioner')
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

  end subroutine euler_calcPreconditioner

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcJacobian(rproblemLevel, rtimestep, rsolver,&
                                rsol, rsol0, bfailure, rcollection)

!<description>
    ! This callback subroutine computes the Jacobian matrix for the
    !  compressible Euler/Navier-Stokes equations
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: rsol0

    ! Newton subiteration failed, return to defect correction
    logical, intent(IN) :: bfailure
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT), target :: rproblemLevel

    ! nonlinear solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! current solution vector
    type(t_vectorBlock), intent(INOUT) :: rsol

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
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

  end subroutine euler_applyJacobian

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcResidual(rproblemLevel, rtimestep, rsolver,&
                                rsol, rsol0, rrhs, rres, ite, rcollection)

!<description>
    ! This subroutine computes the nonlinear residual vector and the
    ! constant right-hand side (only in the first iteration).
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(IN) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: rsol0

    ! number of nonlinear iteration
    integer, intent(IN) :: ite
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver
    
    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsol

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
    call euler_calcPreconditioner(rproblemLevel, rtimestep, rsolver, rsol, rcollection)

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
                                   rsol, euler_calcFluxGalerkin1d, rtimestep%dStep, .true., rres)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsol, euler_calcFluxGalerkin2d, rtimestep%dStep, .true., rres)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsol, euler_calcFluxGalerkin3d, rtimestep%dStep, .true., rres)
        end select

        !-----------------------------------------------------------------------
        
      case (AFCSTAB_UPWIND)

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
                                     rsol, euler_calcFluxScalarDiss1d, rtimestep%dStep, .true., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsol, euler_calcFluxScalarDiss2d, rtimestep%dStep, .true., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsol, euler_calcFluxScalarDiss3d, rtimestep%dStep, .true., rres)
          end select

        case (DISSIPATION_TENSOR)

          ! Assemble divergence operator with tensorial dissipation
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsol, euler_calcFluxTensorDiss1d, rtimestep%dStep, .true., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsol, euler_calcFluxTensorDiss2d, rtimestep%dStep, .true., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsol, euler_calcFluxTensorDiss3d, rtimestep%dStep, .true., rres)
          end select

        case DEFAULT
          call output_line('Invalid type of dissipation!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_calcResidual')
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
                                      rsol, euler_calcFluxTVD1d, euler_calcCharacteristics1d,&
                                      rproblemLevel%Rafcstab(inviscidAFC), rtimestep%dStep, .true., rres)
        case (NDIM2D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsol, euler_calcFluxTVD2d, euler_calcCharacteristics2d,&
                                      rproblemLevel%Rafcstab(inviscidAFC), rtimestep%dStep, .true., rres)
        case (NDIM3D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsol, euler_calcFluxTVD3d, euler_calcCharacteristics3d,&
                                      rproblemLevel%Rafcstab(inviscidAFC), rtimestep%dStep, .true., rres)
        end select

        !-----------------------------------------------------------------------

      case DEFAULT
        call output_line('Invalid type of stabilisation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_calcResidual')
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
        
        do iblock = 1, rsol%nblocks
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                   rsol%RvectorBlock(iblock),&
                                   rrhs%RvectorBlock(iblock),&
                                   1._DP, 1.0_DP-rtimestep%theta)
        end do

      case (MASS_CONSISTENT)

        ! Compute the constant right-hand side
        !
        !   $ rhs = M_C*U^n + (1-theta)*res $

        consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
        call lsysbl_vectorLinearComb(rres, rrhs, 1.0_DP, 0.0_DP)
        
        do iblock = 1, rsol%nblocks
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                   rsol%RvectorBlock(iblock),&
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
        
        do iblock = 1, rsol%nblocks
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                   rsol%RvectorBlock(iblock),&
                                   rres%RvectorBlock(iblock), -1._DP, 1.0_DP)
        end do

      case (MASS_CONSISTENT)

        ! Adopt the constant right-hand side vector
        !
        !   $ res = rhs-M_C*U^{(m)} $
        
        consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
        call lsysbl_vectorLinearComb(rrhs, rres, 1.0_DP, 0.0_DP)
        
        do iblock = 1, rsol%nblocks
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                   rsol%RvectorBlock(iblock),&
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
                                   rsol, euler_calcFluxGalerkin1d,&
                                   rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsol, euler_calcFluxGalerkin2d,&
                                   rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsol, euler_calcFluxGalerkin3d,&
                                   rtimestep%theta*rtimestep%dStep, .false., rres)
        end select

        !-----------------------------------------------------------------------

      case (AFCSTAB_UPWIND)

        ! Compute the low-order residual
        !
        !   $ res = res + dt*theta*L(U^{(m)})*U^(m) $

        idissipationtype = collct_getvalue_int(rcollection, 'idissipationtype')

        select case(idissipationtype)

        case (DISSIPATION_SCALAR)



          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsol, euler_calcFluxScalarDiss1d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsol, euler_calcFluxScalarDiss2d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsol, euler_calcFluxScalarDiss3d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          end select

        case (DISSIPATION_TENSOR)

          ! Assemble divergence operator with tensorial dissipation

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsol, euler_calcFluxTensorDiss1d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsol, euler_calcFluxTensorDiss2d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsol, euler_calcFluxTensorDiss3d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          end select

        case DEFAULT
          call output_line('Invalid type of dissipation!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_calcResidual')
          call sys_halt()
        end select

        !-----------------------------------------------------------------------
          
      case (AFCSTAB_FEMTVD)

        ! Compute th elow-order residual + FEM-TVD stabilization
        !
        !   $ res = res + dt*theta*L(U^{(m)})*U^{(m)} + F(U^{(m)}) $

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsol, euler_calcFluxTVD1d, euler_calcCharacteristics1d,&
                                      rproblemLevel%Rafcstab(inviscidAFC),&
                                      rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM2D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsol, euler_calcFluxTVD2d, euler_calcCharacteristics2d,&
                                      rproblemLevel%Rafcstab(inviscidAFC),&
                                      rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM3D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsol, euler_calcFluxTVD3d, euler_calcCharacteristics3d,&
                                      rproblemLevel%Rafcstab(inviscidAFC),&
                                      rtimestep%theta*rtimestep%dStep, .false., rres)
        end select

        !-----------------------------------------------------------------------

      case DEFAULT
        call output_line('Invalid type of stabilisation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_calcResidual')
        call sys_halt()
      end select

    end if
        
  end subroutine euler_calcResidual

  !*****************************************************************************
  
!<subroutine>

  subroutine euler_calcRHS(rproblemLevel, rtimestep, rsolver,&
                           rsol, rsol0, rrhs, istep, rcollection)

!<description>
    ! This subroutine computes the right-hand side vector
    ! used in the explicit Lax-Wendroff time-stepping scheme
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(IN) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: rsol0

    ! number of explicit step
    integer, intent(IN) :: istep
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver
    
    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsol

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
      do iblock = 1, rsol%nblocks
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                 rsol%RvectorBlock(iblock),&
                                 rrhs%RvectorBlock(iblock), 1.0_DP, 0.0_DP)
      end do

    case(MASS_CONSISTENT)
      
      ! Initialize the right-hand side vector
      !
      !  $ M_C*U

      consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
      do iblock = 1, rsol%nblocks
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                 rsol%RvectorBlock(iblock),&
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
                                 rsol, euler_calcFluxGalerkin1d, dscale, .false., rrhs)
      case (NDIM2D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                 rsol, euler_calcFluxGalerkin2d, dscale, .false., rrhs)
      case (NDIM3D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                 rsol, euler_calcFluxGalerkin3d, dscale, .false., rrhs)
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
                                   rsol, euler_calcFluxScalarDiss1d, dscale, .false., rrhs)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsol, euler_calcFluxScalarDiss2d, dscale, .false., rrhs)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsol, euler_calcFluxScalarDiss3d, dscale, .false., rrhs)
        end select

      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                   rsol, euler_calcFluxTensorDiss1d, dscale, .false., rrhs)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsol, euler_calcFluxTensorDiss2d, dscale, .false., rrhs)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsol, euler_calcFluxTensorDiss3d, dscale, .false., rrhs)
        end select

      case DEFAULT
        call output_line('Invalid type of dissipation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_calcRHS')
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
                                    rsol, euler_calcFluxTVD1d, euler_calcCharacteristics1d,&
                                    rproblemLevel%Rafcstab(inviscidAFC), dscale, .false., rrhs)
      case (NDIM2D)
        call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                    rsol, euler_calcFluxTVD2d, euler_calcCharacteristics2d,&
                                    rproblemLevel%Rafcstab(inviscidAFC), dscale, .false., rrhs)
      case (NDIM3D)
        call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                    rsol, euler_calcFluxTVD3d, euler_calcCharacteristics3d,&
                                    rproblemLevel%Rafcstab(inviscidAFC), dscale, .false., rrhs)
      end select

      !-----------------------------------------------------------------------

    case DEFAULT
      call output_line('Invalid type of stabilization!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_calcRHS')
      call sys_halt()
    end select

    ! Stop time measurement for global operator
    call stat_stopTimer(rtimer)
    
  end subroutine euler_calcRHS

end module euler_callback

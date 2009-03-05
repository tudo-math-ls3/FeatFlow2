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

  subroutine euler_nlsolverCallback(rproblemLevel, rtimestep, rsolver,&
                                    rsolution, rsolutionInitial, rrhs, rres,&
                                    istep, ioperationSpec, rcollection, istatus)

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
      
      call euler_calcPreconditioner(rproblemLevel, rtimestep, rsolver,&
                                    rsolution, rcollection)
    end if
    
    
!!$    ! Do we have to calculate the constant right-hand side?
!!$    ! --------------------------------------------------------------------------
!!$    if ((iand(ioperationSpec, NLSOL_OPSPEC_CALCRHS)  .ne. 0)) then
!!$
!!$      call euler_calcrhs(rproblemLevel, rtimestep, rsolver,&
!!$                         rsolution, rvector, rcollection)
!!$    end if

    ! Do we have to calculate the residual
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then
      
      call euler_calcResidual(rproblemLevel, rtimestep, rsolver,&
                              rsolution, rsolutionInitial,&
                              rrhs, rres, istep, rcollection)
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

  subroutine euler_calcPreconditioner(rproblemLevel, rtimestep, rsolver,&
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
                                      rsolution, euler_calcMatrixGalerkinDiag1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, euler_calcMatrixGalerkinDiag2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, euler_calcMatrixGalerkinDiag3d, 1.0_DP, .true.,&
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
                                      rsolution, euler_calcMatrixScalarDissDiag1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, euler_calcMatrixScalarDissDiag2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, euler_calcMatrixScalarDissDiag3d, 1.0_DP, .true.,&
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
                                      rsolution, euler_calcMatrixTensorDissDiag1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, euler_calcMatrixTensorDissDiag2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, euler_calcMatrixTensorDissDiag3d, 1.0_DP, .true.,&
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
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                      rsolution, euler_calcMatrixGalerkin1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, euler_calcMatrixGalerkin2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, euler_calcMatrixGalerkin3d, 1.0_DP, .true.,&
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
                                      rsolution, euler_calcMatrixScalarDiss1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, euler_calcMatrixScalarDiss2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, euler_calcMatrixScalarDiss3d, 1.0_DP, .true.,&
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
                                      rsolution, euler_calcMatrixTensorDiss1d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM2D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, euler_calcMatrixTensorDiss2d, 1.0_DP, .true.,&
                                      rproblemLevel%RmatrixBlock(systemMatrix))
        case (NDIM3D)
          call gfsys_buildDivOperator(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, euler_calcMatrixTensorDiss3d, 1.0_DP, .true.,&
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
    call euler_calcPreconditioner(rproblemLevel, rtimestep, rsolver, rsolution, rcollection)

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
                                   rsolution, euler_calcFluxGalerkin1d, rtimestep%dStep, .true., rres)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsolution, euler_calcFluxGalerkin2d, rtimestep%dStep, .true., rres)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsolution, euler_calcFluxGalerkin3d, rtimestep%dStep, .true., rres)
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
                                     rsolution, euler_calcFluxScalarDiss1d, rtimestep%dStep, .true., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, euler_calcFluxScalarDiss2d, rtimestep%dStep, .true., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, euler_calcFluxScalarDiss3d, rtimestep%dStep, .true., rres)
          end select

        case (DISSIPATION_TENSOR)

          ! Assemble divergence operator with tensorial dissipation
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, euler_calcFluxTensorDiss1d, rtimestep%dStep, .true., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, euler_calcFluxTensorDiss2d, rtimestep%dStep, .true., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, euler_calcFluxTensorDiss3d, rtimestep%dStep, .true., rres)
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
                                      rsolution, euler_calcFluxTVD1d, euler_calcCharacteristics1d,&
                                      rproblemLevel%Rafcstab(inviscidAFC), rtimestep%dStep, .true., rres)
        case (NDIM2D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, euler_calcFluxTVD2d, euler_calcCharacteristics2d,&
                                      rproblemLevel%Rafcstab(inviscidAFC), rtimestep%dStep, .true., rres)
        case (NDIM3D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, euler_calcFluxTVD3d, euler_calcCharacteristics3d,&
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
                                   rsolution, euler_calcFluxGalerkin1d,&
                                   rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsolution, euler_calcFluxGalerkin2d,&
                                   rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsolution, euler_calcFluxGalerkin3d,&
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
                                     rsolution, euler_calcFluxScalarDiss1d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, euler_calcFluxScalarDiss2d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, euler_calcFluxScalarDiss3d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          end select

        case (DISSIPATION_TENSOR)

          ! Assemble divergence operator with tensorial dissipation

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, euler_calcFluxTensorDiss1d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, euler_calcFluxTensorDiss2d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                     rsolution, euler_calcFluxTensorDiss3d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          end select

        case DEFAULT
          call output_line('Invalid type of dissipation!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_calcResidual')
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
                                      rsolution, euler_calcFluxTVD1d, euler_calcCharacteristics1d,&
                                      rproblemLevel%Rafcstab(inviscidAFC),&
                                      rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM2D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                      rsolution, euler_calcFluxTVD2d, euler_calcCharacteristics2d,&
                                      rproblemLevel%Rafcstab(inviscidAFC),&
                                      rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM3D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                      rsolution, euler_calcFluxTVD3d, euler_calcCharacteristics3d,&
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
                                 rsolution, euler_calcFluxGalerkin1d, dscale, .false., rrhs)
      case (NDIM2D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                 rsolution, euler_calcFluxGalerkin2d, dscale, .false., rrhs)
      case (NDIM3D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                 rsolution, euler_calcFluxGalerkin3d, dscale, .false., rrhs)
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
                                   rsolution, euler_calcFluxScalarDiss1d, dscale, .false., rrhs)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsolution, euler_calcFluxScalarDiss2d, dscale, .false., rrhs)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsolution, euler_calcFluxScalarDiss3d, dscale, .false., rrhs)
        end select

      case (DISSIPATION_TENSOR)

        ! Assemble divergence operator with tensorial dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                   rsolution, euler_calcFluxTensorDiss1d, dscale, .false., rrhs)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsolution, euler_calcFluxTensorDiss2d, dscale, .false., rrhs)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                   rsolution, euler_calcFluxTensorDiss3d, dscale, .false., rrhs)
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
                                    rsolution, euler_calcFluxTVD1d, euler_calcCharacteristics1d,&
                                    rproblemLevel%Rafcstab(inviscidAFC), dscale, .false., rrhs)
      case (NDIM2D)
        call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                    rsolution, euler_calcFluxTVD2d, euler_calcCharacteristics2d,&
                                    rproblemLevel%Rafcstab(inviscidAFC), dscale, .false., rrhs)
      case (NDIM3D)
        call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                                    rsolution, euler_calcFluxTVD3d, euler_calcCharacteristics3d,&
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

  !*****************************************************************************

!<subroutine>

  subroutine euler_setBoundary(rproblemLevel, rtimestep, rsolver,&
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
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_setBoundary')
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
                               euler_calcBoundaryvalues1d, istatus)
    case (NDIM2D)
      call bdrf_filterSolution(rsolver%rboundaryCondition,&
                               rproblemLevel%rtriangulation,&
                               rproblemLevel%RmatrixBlock(imatrix),&
                               rsolution, rres, rsolutionInitial,&
                               rtimestep%dTime,&
                               rproblemLevel%p_rproblem%rboundary,&
                               euler_calcBoundaryvalues2d, istatus)
    case (NDIM3D)
      call bdrf_filterSolution(rsolver%rboundaryCondition,&
                               rproblemLevel%rtriangulation,&
                               rproblemLevel%RmatrixBlock(imatrix),&
                               rsolution, rres, rsolutionInitial,&
                               rtimestep%dTime,&
                               rproblemLevel%p_rproblem%rboundary,&
                               euler_calcBoundaryvalues3d, istatus)
    case DEFAULT
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_setBoundary')
      call sys_halt()
    end select
    
  end subroutine euler_setBoundary

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcLinearizedFCT(rbdrCond, rproblemLevel, rtimestep, rsolution, rcollection)

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

    real(DP), dimension(NDIM2D) :: XDir2D = (/ 1.0_DP, 0.0_DP /)
    real(DP), dimension(NDIM2D) :: YDir2D = (/ 0.0_DP, 1.0_DP /)
    

    ! local variables
    type(t_matrixScalar), pointer :: p_rmatrix
    type(t_vectorScalar) :: rflux0X, rflux0Y, rfluxX, rfluxY
    real(DP), dimension(:), pointer :: p_MC, p_ML, p_Cx, p_Cy, p_u, p_flux0X, p_flux0Y, p_fluxX, p_fluxY
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal, p_Ksep
    integer :: h_Ksep, templatematrix, lumpedMassMatrix, consistentMassMatrix, coeffMatrix_CX, coeffMatrix_CY, nedge

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

    call lsysbl_getbase_double(rsolution, p_u)

    ! Create diagonal separator
    h_Ksep = ST_NOHANDLE
    call storage_copy(p_rmatrix%h_Kld, h_Ksep)
    call storage_getbase_int(h_Ksep, p_Ksep, p_rmatrix%NEQ+1)

    ! Compute number of edges
    nedge = int(0.5*(p_rmatrix%NA-p_rmatrix%NEQ))

    call lsyssc_createVector(rflux0X, nedge, NVAR2D, .true., ST_DOUBLE)
    call lsyssc_createVector(rflux0Y, nedge, NVAR2D, .true., ST_DOUBLE)
    call lsyssc_createVector(rfluxX, nedge, NVAR2D, .true., ST_DOUBLE)
    call lsyssc_createVector(rfluxY, nedge, NVAR2D, .true., ST_DOUBLE)

    ! Create flux vectors
    call lsyssc_getbase_double(rflux0X, p_flux0X)
    call lsyssc_getbase_double(rflux0Y, p_flux0Y)
    call lsyssc_getbase_double(rfluxX, p_fluxX)
    call lsyssc_getbase_double(rfluxY, p_fluxY)

    ! Build the flux vectors
    call buildFluxVectors2D(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                            NVAR2D, nedge, p_u, rtimestep%dStep, p_MC, p_ML, p_Cx, p_Cy,&
                            p_flux0X, p_flux0Y, p_fluxX, p_fluxY)

    call limitSolutionVector(Xdir2D, p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                             NVAR2D, nedge, p_ML, p_u, p_flux0X, p_fluxX)

    
    ! Set boundary conditions explicitly
    call bdrf_filterVectorExplicit(rbdrCond, rproblemLevel%rtriangulation,&
                                   rsolution, rtimestep%dTime,&
                                   rproblemLevel%p_rproblem%rboundary,&
                                   euler_calcBoundaryvalues2d)

    call limitSolutionVector(YDir2D, p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                             NVAR2D, nedge, p_ML, p_u, p_flux0X, p_fluxX)

    ! Set boundary conditions explicitly
    call bdrf_filterVectorExplicit(rbdrCond, rproblemLevel%rtriangulation,&
                                   rsolution, rtimestep%dTime,&
                                   rproblemLevel%p_rproblem%rboundary,&
                                   euler_calcBoundaryvalues2d)

    ! Release flux vectors
    call storage_free(h_Ksep)
    call lsyssc_releaseVector(rflux0X)
    call lsyssc_releaseVector(rfluxX)
    call lsyssc_releaseVector(rflux0Y)
    call lsyssc_releaseVector(rfluxY)

  contains

    subroutine buildFluxVectors2D(Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, NEDGE, u, dscale,&
                                  MC, ML, Cx, Cy, flux0X, flux0Y, fluxX, fluxY)
      
      real(DP), dimension(NVAR,NEQ), intent(IN) :: u
      real(DP), dimension(:), intent(IN) :: Cx,Cy,MC,ML
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR,NEDGE
      
      real(DP), dimension(NVAR,NEDGE), intent(INOUT) :: flux0X,flux0Y,fluxX,fluxY
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      integer :: ij,ji,i,j, iedge

      real(DP), dimension(:,:), allocatable :: rX,rY
      real(DP), dimension(NDIM2D) :: a_ij, s_ij
      real(DP), dimension(NVAR*NVAR) :: R_ij
      real(DP), dimension(NVAR) :: Lbd_ij, W_ij, Daux, Daux1, Daux2, Wp_ij, Wm_ij

      allocate(rX(nvar,neq), rY(nvar,neq)); rX=0; rY=0

      ! Reset Ksep
      Ksep = Kld; iedge = 0
      
      do i = 1, NEQ
        
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1; iedge = iedge+1

          ! Compute skew-symmetric and symmetric coefficient
          a_ij(1) = 0.5_DP*(Cx(ij)-Cx(ji));   a_ij(2) = 0.5_DP*(Cy(ij)-Cy(ji))
          s_ij(1) = 0.5_DP*(Cx(ij)+Cx(ji));   s_ij(2) = 0.5_DP*(Cy(ij)+Cy(ji))

          
          ! Calculate characteristic variables for x-direction
          call euler_calcCharacteristics2d(u(:,i), u(:,j), XDir2D, W_ij, Lbd_ij, R_ij)
          Daux = a_ij(1)*Lbd_ij; Wm_ij = min(0.0_DP, Daux)*W_ij; Wp_ij = max(0.0_DP, Daux)*W_ij

          call DGEMV('n', NVAR, NVAR, 1.0_DP, R_ij, NVAR, Wm_ij, 1, 0.0_DP, Daux1, 1)
          rX(:,i) = rX(:,i) - 2*dscale*Daux1
          
          call DGEMV('n', NVAR, NVAR, 1.0_DP, R_ij, NVAR, Wp_ij, 1, 0.0_DP, Daux1, 1)
          rX(:,j) = rX(:,j) - 2*dscale*Daux1

          Daux1 = max(0.0_DP, min(abs(Daux)+s_ij(1)*Lbd_ij,2*abs(Daux)))*W_ij
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, Daux1, 1, 0.0_DP, Daux2, 1)
          fluxX(:,iedge) = -Daux2

          
          ! Calculate characteristic variables for y-direction
          call euler_calcCharacteristics2d(u(:,i), u(:,j), YDir2D, W_ij, Lbd_ij, R_ij)
          Daux = a_ij(2)*Lbd_ij; Wm_ij = min(0.0_DP, Daux)*W_ij; Wp_ij = max(0.0_DP, Daux)*W_ij

          call DGEMV('n', NVAR, NVAR, 1.0_DP, R_ij, NVAR, Wm_ij, 1, 0.0_DP, Daux, 1)
          rY(:,i) = rY(:,i) - 2*dscale*Daux
          
          call DGEMV('n', NVAR, NVAR, 1.0_DP, R_ij, NVAR, Wp_ij, 1, 0.0_DP, Daux, 1)
          rY(:,j) = rY(:,j) - 2*dscale*Daux

          Daux1 = max(0.0_DP, min(abs(Daux)+s_ij(2)*Lbd_ij,2*abs(Daux)))*W_ij
          call DGEMV('n', NVAR, NVAR, dscale, R_ij, NVAR, Daux1, 1, 0.0_DP, Daux2, 1)
          fluxY(:,iedge) = -Daux2

        end do
      end do

      do i = 1, NEQ
        rX(:,i) = rX(:,i)/ML(i)
        rY(:,i) = rY(:,i)/ML(i)
      end do

      flux0X = fluxX
      flux0Y = fluxY

      ! Reset Ksep
      Ksep = Kld;  iedge = 0
      
      ! Add contribution of consistent mass matrix
      do i = 1, NEQ
        
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; iedge = iedge+1

          fluxX(:,iedge) = MC(ij)*(rX(:,i)-rX(:,j))+fluxX(:,iedge)
          fluxY(:,iedge) = MC(ij)*(rY(:,i)-rY(:,j))+fluxY(:,iedge)
          
        end do
      end do

      deallocate(rX, rY)

    end subroutine buildFluxVectors2D


    subroutine limitSolutionVector(Dweight, Kld, Kcol, Kdiagonal, Ksep, NEQ, NVAR, NEDGE, ML, u, flux0, flux)

      include 'intf_gfsyscallback.inc'

      real(DP), dimension(NVAR,NEDGE), intent(INOUT) :: flux0,flux
      real(DP), dimension(:), intent(IN) :: Dweight, ML
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NVAR,NEDGE
      
      real(DP), dimension(NVAR,NEQ), intent(INOUT) :: u
      integer, dimension(:), intent(INOUT) :: Ksep

      ! local variables
      real(DP), dimension(:,:), allocatable :: pp,pm,qp,qm,rp,rm,r
      real(DP), dimension(NVAR*NVAR) :: R_ij, L_ij
      real(DP), dimension(NVAR) :: Lbd_ij, W_ij, Daux, F_ij, F0_ij
      integer :: ij,ji,i,j,iedge,ivar

      allocate(pp(nvar,neq), pm(nvar,neq), qp(nvar,neq), qm(nvar,neq), rp(nvar,neq), rm(nvar,neq), r(nvar, neq))
      
      pp=0; pm=0; qp=0; qm=0; rp=1; rm=1; r=0
      
      ! Reset Ksep
      Ksep = Kld; iedge = 0
      
      ! Add contribution of consistent mass matrix
      do i = 1, NEQ
        
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; iedge = iedge+1

          call euler_calcCharacteristics2d(u(:,i), u(:,j), Dweight, W_ij, Lbd_ij, R_ij, L_ij)

          call DGEMV('n', NVAR, NVAR, 1.0_DP, L_ij, NVAR, flux(:,iedge), 1, 0.0_DP, F_ij, 1)
          call DGEMV('n', NVAR, NVAR, 1.0_DP, L_ij, NVAR, flux0(:,iedge), 1, 0.0_DP, F0_ij, 1)

          F_ij = minmod(F_ij, F0_ij);   flux(:,iedge) = F_ij

          ! Sums of raw antidiffusive fluxes
          pp(:,i) = pp(:,i) + max(0.0_DP,  F_ij)
          pp(:,j) = pp(:,j) + max(0.0_DP, -F_ij)
          pm(:,i) = pm(:,i) + min(0.0_DP,  F_ij)
          pm(:,j) = pm(:,j) + min(0.0_DP, -F_ij)

          ! Sums of admissible edge contributions
          qp(:,i) = max(qp(:,i),  F0_ij)
          qp(:,j) = max(qp(:,j), -F0_ij)
          qm(:,i) = min(qm(:,i),  F0_ij)
          qm(:,j) = min(qm(:,j), -F0_ij)

        end do
      end do


      ! Compute nodal correction factors
      do i = 1, NEQ
        do ivar = 1, nvar
          if (pp(ivar,i) >  SYS_EPSREAL) rp(ivar,i) = min(1.0_DP, ML(i)*qp(ivar,i)/pp(ivar,i))
          if (pm(ivar,i) < -SYS_EPSREAL) rm(ivar,i) = min(1.0_DP, ML(i)*qm(ivar,i)/pm(ivar,i))
        end do
      end do


      ! Reset Ksep
      Ksep = Kld; iedge = 0

      ! Add contribution of consistent mass matrix
      do i = 1, NEQ
        
        do ij = Kdiagonal(i)+1, Kld(i+1)-1
          
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; iedge = iedge+1

          F_ij = flux(:,iedge)

          where (F_ij > 0)
            F_ij = min(rp(:,i), rm(:,j))*F_ij
          elsewhere
            F_ij = min(rm(:,i), rp(:,j))*F_ij
          end where

          call euler_calcCharacteristics2d(u(:,i), u(:,j), Dweight, W_ij, Lbd_ij, R_ij)

          call DGEMV('n', NVAR, NVAR, 1.0_DP, R_ij, NVAR, F_ij, 1, 0.0_DP, Daux, 1)

          r(:,i) = r(:,i) + Daux
          r(:,j) = r(:,j) - Daux

        end do
      end do

      do i = 1, NEQ
        u(:,i) = u(:,i) + r(:,i)/ML(i)
      end do

      deallocate(pp,pm,qp,qm,rp,rm,r)
      
    end subroutine limitSolutionVector


    pure function minmod(a,b)
      real(DP), dimension(:), intent(IN) :: a,b
      real(DP), dimension(size(a)) :: minmod

      where (a > 0 .and. b > 0)
        minmod = min(a,b)
      elsewhere (a < 0 .and. b < 0)
        minmod = max(a,b)
      elsewhere
        minmod = 0
      end where
    end function minmod

  end subroutine euler_calcLinearizedFCT

end module euler_callback

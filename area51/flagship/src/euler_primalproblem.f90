!##############################################################################
!# ****************************************************************************
!# <name> euler_primalproblem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to solve the primal 
!# problem for the compressible Euler/Navier-Stokes equations
!#
!# The following routines are available:
!#
!# 1.) euler_calcPrimalPrecond
!#     -> update the global preconditioner for the solution of the
!#        primal problem for the compressible Navier-Stokes equations
!#
!# 2.) fcb_calcPrimalRHS
!#     -> compute the right-hand side vector for the primal problem
!#
!# 3.) fcb_calcPrimalResidual
!#     -> compute the defect vector for the primal problem
!#
!# 4.) fcb_calcPrimalJacobian
!#     -> compute the Jacobian matrix for the primal problem
!#
!# 5.) fcb_applyPrimalJacobian
!#     -> apply the Jacobain matrix for the primal problem
!#
!# </purpose>
!##############################################################################

module euler_primalproblem

  use afcstabilisation
  use fsystem
  use genoutput
  use groupfemsystem
  use linearsystemblock
  use linearsystemscalar
  use statistics
  use storage

  use boundaryfilter
  use euler_basic
  use euler_callback1d
  use euler_callback2d
  use euler_callback3d
  use problem
  use solver

  implicit none

  private
  public :: euler_calcPrimalPrecond
  public :: fcb_calcPrimalRHS
  public :: fcb_calcPrimalResidual
  public :: fcb_calcPrimalJacobian
  public :: fcb_applyPrimalJacobian

contains

  !*****************************************************************************

!<subroutine>

  subroutine euler_calcPrimalPrecond(rproblemLevel, rtimestep, rsolver,&
                                  ru, imode, nlminOpt, nlmaxOpt)

!<description>
    ! This subroutine updates the global system operator for the compressible
    ! Navier-Stokes equations and configures the linear solver structure accordingly.
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! current solution vector
    type(t_vectorBlock), intent(IN), target :: ru

    ! mode: (1) primal or (2) dual problem
    integer, intent(IN) :: imode

    ! OPTIONAL: minimal multigrid level
    integer, intent(IN), optional :: nlminOpt
    
    ! OPTIONAL: maximal multigrid level
    integer, intent(IN), optional :: nlmaxOpt
!</input>

!<inputoutput>
    ! finest multigrid level tostart with
    type(t_problemLevel), intent(INOUT), target :: rproblemLevel

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_problemLevel), pointer :: rproblemLevelTmp
    type(t_matrixBlock), pointer :: p_rmatrixA
    type(t_vectorBlock), pointer :: p_ru
    type(t_vectorBlock), target :: ruFine
    type(t_vectorBlock) :: ruCoarse   
    integer :: nlmin,nlmax,iblock,jblock
    
    ! Start time measurement for global operator
    call stat_startTimer(rtimer_assembly_matrix, STAT_TIMERSHORT)

    ! Check if vector has correct dimensions
    select case(isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)
      if ((ru%nblocks .ne. 1) .or.&
          (ru%RvectorBlock(1)%NVAR .ne. euler_getNVAR(rproblemLevel))) then
        call output_line('Invalid format of solution vector!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
        call sys_halt()
      end if

    case (SYSTEM_BLOCKFORMAT)
      if (ru%nblocks .ne. rproblemLevel%rtriangulation%ndim) then
        call output_line('Invalid format of solution vector!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
        call sys_halt()
      else
        do iblock = 1, ru%nblocks
          if (ru%RvectorBlock(iblock)%NVAR .ne. 1) then
            call output_line('Invalid format of solution vector!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
            call sys_halt()
          end if
        end do
      end if

    case DEFAULT
      call output_line('Unsupported system format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
      call sys_halt()
    end select
    
    
    ! Set minimal/maximal levels: If given by the user, adopt these values.
    ! Otherwise, use the values of the given multigrid solver structure.
    if (present(nlminOpt)) then
      nlmin = nlminOpt
    else
      nlmin = solver_getNLMIN(rsolver,rproblemLevel%ilev)
    end if

    if (present(nlmaxOpt)) then
      nlmax = nlmaxOpt
    else
      nlmax = solver_getNLMAX(rsolver,rproblemLevel%ilev)
    end if
    
    
    ! Update the matrices on all required multigrid levels NLMIN <= ILEV <= NLMAX
    rproblemLevelTmp => rproblemLevel
    p_ru   => ru
    do while(associated(rproblemLevelTmp))
      
      ! Do we have to assemble operators for this level?
      if (rproblemLevelTmp%ilev > nlmax) then

        ! Do we have to restrict the given solution vector?
        if (nlmin < nlmax) then
          call lsysbl_createVectorBlock(rproblemLevelTmp%p_rproblemLevelCoarse%rdiscretisation, &
                                        ruCoarse, .true., ST_DOUBLE)
          call solver_restrictionBlock(rproblemLevelTmp%rtriangulation,&
                                       rproblemLevelTmp%p_rproblemLevelCoarse%rtriangulation,&
                                       p_ru, ruCoarse)
          call lsysbl_swapVectors(ruCoarse, ruFine)
          call lsysbl_releaseVector(ruCoarse)
          p_ru => ruFine
        end if
        
        rproblemLevelTmp => rproblemLevelTmp%p_rproblemLevelCoarse
        cycle
      elseif (rproblemLevelTmp%ilev < nlmin) then
        exit
      end if
      
      
      ! Set pointer
      p_rmatrixA => rproblemLevelTmp%RmatrixBlock(CNSE_MATRIX_A)
      
      ! Check if discrete operator needs to be initialized
      if (rproblemLevelTmp%Rafcstab(CNSE_AFCSTAB_INVISCID)%iSpec .eq. AFCSTAB_UNDEFINED)&
          call gfsys_initStabilisation(p_rmatrixA,&
                                       rproblemLevelTmp%Rafcstab(CNSE_AFCSTAB_INVISCID))
      
      !-------------------------------------------------------------------------
      ! Assemble divergence operator
      !-------------------------------------------------------------------------

      ! What kind of coupling is applied?
      select case(icoupled)
        
      ! Assemble block-diagonal divergence operator
      case (FLOW_SEGREGATED)
        
        ! What kind of preconditioner is applied?
        select case(iprecond)
        case (AFCSTAB_GALERKIN)
          ! Assemble divergence operator for standard Galerin scheme
          select case(rproblemLevelTmp%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                        p_ru, fcb_calcMatrixGalerkinDiag1d,&
                                        1._DP, .true., p_rmatrixA)
          case (NDIM2D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                        p_ru, fcb_calcMatrixGalerkinDiag2d,&
                                        1._DP, .true., p_rmatrixA)
          case (NDIM3D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                        p_ru, fcb_calcMatrixGalerkinDiag3d,&
                                        1._DP, .true., p_rmatrixA)
          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
            call sys_halt()
          end select

          
        case (AFCSTAB_SCALARDISSIPATION)
          ! Assemble divergence operator with scalar dissipation
          select case(rproblemLevelTmp%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                        p_ru, fcb_calcMatrixScalarDissDiag1d,&
                                        1._DP, .true., p_rmatrixA)
          case (NDIM2D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                        p_ru, fcb_calcMatrixScalarDissDiag2d,&
                                        1._DP, .true., p_rmatrixA)
          case (NDIM3D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                        p_ru, fcb_calcMatrixScalarDissDiag3d,&
                                        1._DP, .true., p_rmatrixA)
          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
            call sys_halt()
          end select

          
        case (AFCSTAB_TENSORDISSIPATION)
          ! Assemble divergence operator with tensorial dissipation
          select case(rproblemLevelTmp%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                        p_ru, fcb_calcMatrixTensorDissDiag1d,&
                                        1._DP, .true., p_rmatrixA)
          case (NDIM2D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                        p_ru, fcb_calcMatrixTensorDissDiag2d,&
                                        1._DP, .true., p_rmatrixA)
          case (NDIM3D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                        p_ru, fcb_calcMatrixTensorDissDiag3d,&
                                        1._DP, .true., p_rmatrixA)
          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
            call sys_halt()
          end select

          
        case DEFAULT
          call output_line('Invalid type of dissipation!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
          call sys_halt()
        end select
        
        
      ! Assemble full block transport operator
      case (FLOW_ALLCOUPLED)
          
        ! What kind of preconditioner is applied?
        select case(iprecond)         
        case (AFCSTAB_GALERKIN)
          ! Assemble divergence operator for standard Galerin scheme
          select case(rproblemLevelTmp%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                        p_ru, fcb_calcMatrixGalerkin1d,&
                                        1._DP, .true., p_rmatrixA)
          case (NDIM2D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                        p_ru, fcb_calcMatrixGalerkin2d,&
                                        1._DP, .true., p_rmatrixA)
          case (NDIM3D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                        p_ru, fcb_calcMatrixGalerkin3d,&
                                        1._DP, .true., p_rmatrixA)
          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
            call sys_halt()
          end select

          
        case (AFCSTAB_SCALARDISSIPATION)             
          ! Assemble divergence operator with scalar dissipation
          select case(rproblemLevelTmp%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                        p_ru, fcb_calcMatrixScalarDiss1d,&
                                        1._DP, .true., p_rmatrixA)
          case (NDIM2D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                        p_ru, fcb_calcMatrixScalarDiss2d,&
                                        1._DP, .true., p_rmatrixA)
          case (NDIM3D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                        p_ru, fcb_calcMatrixScalarDiss3d,&
                                        1._DP, .true., p_rmatrixA)
          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
            call sys_halt()
          end select


        case (AFCSTAB_TENSORDISSIPATION)
          ! Assemble divergence operator with tensorial dissipation
          select case(rproblemLevelTmp%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                        p_ru, fcb_calcMatrixTensorDiss1d,&
                                        1._DP, .true., p_rmatrixA)
          case (NDIM2D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                        p_ru, fcb_calcMatrixTensorDiss2d,&
                                        1._DP, .true., p_rmatrixA)
          case (NDIM3D)
            call gfsys_buildDivOperator(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                        p_ru, fcb_calcMatrixTensorDiss3d,&
                                        1._DP, .true., p_rmatrixA)
          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
            call sys_halt()
          end select
          
          
        case DEFAULT
          call output_line('Invalid type of dissipation!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
          call sys_halt()
        end select
          

      case DEFAULT
        call output_line('Invalid flow coupling!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
        call sys_halt()
      end select
          
      
      !-------------------------------------------------------------------------
      ! Assemble the global system operator
      !-------------------------------------------------------------------------
      
      if (isystemFormat .eq. SYSTEM_INTERLEAVEFORMAT) then

        ! What type of flow are we?
        select case(iflowtype)
        case (FLOW_TRANSIENT,&
              FLOW_PSEUDOTRANSIENT)
          ! Compute the global operator for transient flows
          !
          !     A = blockdiag(M_L)-theta*dt*L
          !
          if (rtimestep%theta .gt. 0.0_DP) then
            call lsyssc_MatrixLinearComb(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_ML), 1._DP,&
                                         rproblemLevelTmp%Rmatrix(CNSE_MATRIX_A),&
                                         -rtimestep%theta*rtimestep%dStep,&
                                         rproblemLevelTmp%Rmatrix(CNSE_MATRIX_A),&
                                         .false., .false., .true., .true.)
          else
            call lsyssc_spreadMatrix(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_ML),&
                                     rproblemLevelTmp%Rmatrix(CNSE_MATRIX_A))
          end if
          
        case (FLOW_STEADYSTATE)
          ! Compute the global operator for steady-state flow
          !
          !     A = -L
          !
          call lsyssc_scaleMatrix(rproblemLevelTmp%Rmatrix(CNSE_MATRIX_A), -1._DP)
          
        case DEFAULT
          call output_line('Invalid flow type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
          call sys_halt()
        end select

      else
        
        ! What type of flow are we?
        select case(iflowtype)
        case (FLOW_TRANSIENT, FLOW_PSEUDOTRANSIENT)
          ! Compute the global operator for transient flows
          !
          !     A = diag(M_L)-theta*dt*L
          !
          do iblock = 1, euler_getNVAR(rproblemLevelTmp)
            call lsyssc_MatrixLinearComb(&
                rproblemLevelTmp%Rmatrix(CNSE_MATRIX_ML), 1._DP,&
                rproblemLevelTmp%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(iblock,iblock),&
                -rtimestep%theta*rtimestep%dStep,&
                rproblemLevelTmp%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(iblock,iblock),&
                .false., .false., .true., .true.)
          end do

          ! Scale the off-diagonal blocks by "-theta*dt" if required
          if (icoupled .eq. FLOW_ALLCOUPLED) then
            do iblock = 1, euler_getNVAR(rproblemLevelTmp)
              do jblock = 1, euler_getNVAR(rproblemLevelTmp)
                if (iblock .eq. jblock) cycle
                call lsyssc_scaleMatrix(&
                    rproblemLevelTmp%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(iblock,jblock),&
                    -rtimestep%theta*rtimestep%dStep)
              end do
            end do
          end if
          
        case (FLOW_STEADYSTATE)
          ! Compute the global operator for steady-state flow
          !
          !     A = -L
          !
          select case(icoupled)

          case (FLOW_SEGREGATED)
            do iblock = 1, euler_getNVAR(rproblemLevelTmp)
              call lsyssc_scaleMatrix(&
                  rproblemLevelTmp%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(iblock,iblock),&
                  -1.0_DP)
            end do

          case (FLOW_ALLCOUPLED)
            do iblock = 1, euler_getNVAR(rproblemLevelTmp)
              do jblock = 1, euler_getNVAR(rproblemLevelTmp)
                call lsyssc_scaleMatrix(&
                    rproblemLevelTmp%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(iblock,jblock),&
                    -1.0_DP)
              end do
            enddo

          case DEFAULT
            call output_line('Unsupported block matrix format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
          call sys_halt()
          end select
            
        case DEFAULT
          call output_line('Invalid flow type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalPrecond')
          call sys_halt()
        end select

      end if

      !-------------------------------------------------------------------------
      ! Impose boundary conditions
      !-------------------------------------------------------------------------
      
      call bdrf_filterMatrix(rsolver%rboundaryCondition,&
                             rproblemLevelTmp%rtriangulation, p_rmatrixA, 1.0_DP)


      ! Switch to next coarser multigrid level
      rproblemLevelTmp => rproblemLevelTmp%p_rproblemLevelCoarse

      ! Do we have to restrict the given solution vector?
      if (associated(rproblemLevelTmp) .and. nlmin < nlmax) then
        call lsysbl_createVectorBlock(rproblemLevelTmp%rdiscretisation, &
                                      ruCoarse, .true., ST_DOUBLE)
        call solver_restrictionBlock(rproblemLevelTmp%p_rproblemLevelFine%rtriangulation,&
                                     rproblemLevelTmp%rtriangulation, p_ru, ruCoarse)
        call lsysbl_swapVectors(ruCoarse, ruFine)
        call lsysbl_releaseVector(ruCoarse)
        p_ru => ruFine
      end if
    end do
    
    ! Do we have to release the temporal vector
    if (nlmin < nlmax) call lsysbl_releaseVector(ruFine)

    ! Ok, we updated the nonlinear system operator successfully. Now we still 
    ! have to link it to the solver hierarchy. This is done recursively.
    call euler_updateSolverMatrix(rproblemLevel, rsolver, CNSE_MATRIX_A,&
                               UPDMAT_ALL, nlmin, nlmax)

    ! Finally, we have to update the content of the solver hierarchy
    call solver_updateContent(rsolver)

    ! Stop time measurement for global operator
    call stat_stopTimer(rtimer_assembly_matrix)

  end subroutine euler_calcPrimalPrecond

  !*****************************************************************************

!<subroutine>

  subroutine fcb_calcPrimalRHS(rproblemLevel, rtimestep, rsolver,&
                               ru, ru0, rrhs, istep, imode, istatus)

!<description>
    ! This subroutine computes the right-hand side vector to be
    ! used in the explicit Lax-Wendroff time stepping scheme.
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! number of explicit step
    integer, intent(IN) :: istep

    ! mode: (1) primal or (2) dual problem
    integer, intent(IN) :: imode
!</input>

!<inputoutput>
    ! finest multigrid level to start with
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! nonlinear solver structure
    type(t_solver), intent(INOUT) :: rsolver
    
    ! current solution vector
    type(t_vectorBlock), intent(INOUT) :: ru

    ! constant right-hand side
    type(t_vectorBlock), intent(INOUT) :: rrhs

    ! OPTIONAL: status of the callback function
    integer, intent(INOUT), optional :: istatus
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP) :: dscale
    integer :: iblock

    ! Start time measurement for residual/rhs evaluation
    call stat_startTimer(rtimer_assembly_resrhs, STAT_TIMERSHORT)

    ! Initialize the right-hand side vector
    select case (iflowtype)
    case (FLOW_TRANSIENT,&
          FLOW_PSEUDOTRANSIENT)
      do iblock = 1, ru%nblocks
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                 ru%RvectorBlock(iblock),&
                                 rrhs%RvectorBlock(iblock), 1._DP, 0._DP)
      end do

    case (FLOW_STEADYSTATE)
      call lsysbl_clearVector(rrhs)

    case DEFAULT
      call output_line('Invalid flow type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalRHS')
      call sys_halt()
    end select

    
    ! Comput the right-hand side vector
    select case(rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID)%ctypeAFCstabilisation)      
    case (AFCSTAB_GALERKIN)
      ! Compute weights
      dscale = rtimestep%DmultistepWeights(istep)*&
               rtimestep%dStep*(1._DP-rtimestep%theta)
      
      ! Assemble high-order right-hand side
      !
      !   rhs = weight*(1-theta)*dt*K(U)*U
      !
      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                 ru, fcb_calcFluxGalerkin1d,&
                                 dscale, .false., rrhs)
      case (NDIM2D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                 ru, fcb_calcFluxGalerkin2d,&
                                 dscale, .false., rrhs)
      case (NDIM3D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                 ru, fcb_calcFluxGalerkin3d,&
                                 dscale, .false., rrhs)
      case DEFAULT
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalRHS')
        call sys_halt()
      end select
      
  
    case (AFCSTAB_UPWIND)
      ! Compute weights
      dscale = rtimestep%DmultistepWeights(istep)*&
               rtimestep%dStep*(1._DP-rtimestep%theta)

      ! Assemble low-order right-hand side
      !
      !   rhs = weight*(1-theta)*dt*L(U)*U
      !
      select case(rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID)%idissipation)
      case (AFCSTAB_SCALARDISSIPATION)
        
        ! Scalar dissipation
        !
        !   D = d_{ij}
        !
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                   ru, fcb_calcFluxScalarDiss1d,&
                                   dscale, .false., rrhs)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                   ru, fcb_calcFluxScalarDiss2d,&
                                   dscale, .false., rrhs)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                   ru, fcb_calcFluxScalarDiss3d,&
                                   dscale, .false., rrhs)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalRHS')
          call sys_halt()
        end select
          
        
      case (AFCSTAB_TENSORDISSIPATION)
        
        ! Tensorial dissipation
        !
        !   D = D_{ij}
        !
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                   ru, fcb_calcFluxTensorDiss1d,&
                                   dscale, .false., rrhs)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                   ru, fcb_calcFluxTensorDiss2d,&
                                   dscale, .false., rrhs)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                   ru, fcb_calcFluxTensorDiss3d,&
                                   dscale, .false., rrhs)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalRHS')
          call sys_halt()
        end select
          
        
      case DEFAULT
        call output_line('Invalid type of dissipation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalRHS')
        call sys_halt()
      end select

      
    case (AFCSTAB_FEMFCT)
      ! Compute weights
      dscale = rtimestep%DmultistepWeights(istep)*&
               rtimestep%dStep*(1._DP-rtimestep%theta)

      ! Assemble low-order right-hand side + FEM-FCT stabilization
      !
      !   rhs = weight*(1-theta)*dt*L(U)*U + F(U)
      !
      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                 ru, fcb_calcFluxScalarDiss1d,&
                                 dscale, .false., rrhs)
      case (NDIM2D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                 ru, fcb_calcFluxScalarDiss2d,&
                                 dscale, .false., rrhs)
      case (NDIM3D)
        call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                 ru, fcb_calcFluxScalarDiss3d,&
                                 dscale, .false., rrhs)
      case DEFAULT
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalRHS')
        call sys_halt()
      end select
        

      ! Compute explicit low-order predictor 
      !
      !   u^star = u + M_L^{-1})*[(1-theta)*dt*L(U)*U]
      !
      do iblock = 1, rrhs%nblocks
        call lsyssc_invertedDiagMatVec(rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                       rrhs%RvectorBlock(iblock), 1._DP,&
                                       rrhs%RvectorBlock(iblock))
      end do
      call lsysbl_vectorLinearComb(ru, rrhs, 1._DP, 1._DP)

      ! Compute FEM-FCT stabilization
      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildResidualFCT(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                    rproblemLevel%Rmatrix(CNSE_MATRIX_ML), ru,&
                                    fcb_calcFluxScalarDiss1d,&
                                    fcb_calcRawFluxFCT1d,&
                                    fcb_calcCharacteristics1d,&
                                    rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                    dscale, rtimestep%theta,&
                                    rtimestep%DmultistepWeights(istep)*rtimestep%dStep,&
                                    .false., rrhs, rrhs, rproblemLevel%Rmatrix(CNSE_MATRIX_MC))
      case (NDIM2D)
        call gfsys_buildResidualFCT(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                    rproblemLevel%Rmatrix(CNSE_MATRIX_ML), ru,&
                                    fcb_calcFluxScalarDiss2d,&
                                    fcb_calcRawFluxFCT2d,&
                                    fcb_calcCharacteristics2d,&
                                    rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                    dscale, rtimestep%theta,&
                                    rtimestep%DmultistepWeights(istep)*rtimestep%dStep,&
                                    .false., rrhs, rrhs, rproblemLevel%Rmatrix(CNSE_MATRIX_MC))
      case (NDIM3D)
        call gfsys_buildResidualFCT(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                    rproblemLevel%Rmatrix(CNSE_MATRIX_ML), ru,&
                                    fcb_calcFluxScalarDiss3d,&
                                    fcb_calcRawFluxFCT3d,&
                                    fcb_calcCharacteristics3d,&
                                    rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                    dscale, rtimestep%theta,&
                                    rtimestep%DmultistepWeights(istep)*rtimestep%dStep,&
                                    .false., rrhs, rrhs, rproblemLevel%Rmatrix(CNSE_MATRIX_MC))
      case DEFAULT
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalRHS')
        call sys_halt()
      end select


    case (AFCSTAB_FEMTVD)
      ! Compute weights
      dscale = rtimestep%DmultistepWeights(istep)*&
               rtimestep%dStep*(1._DP-rtimestep%theta)

      ! Assemble low-order right-hand side + FEM-TVD stabilization
      !
      !   rhs = weight*(1-theta)*dt*L(U)*U + F(U)
      !
      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                    ru, fcb_calcFluxTVD1d,&
                                    fcb_calcCharacteristics1d,&
                                    rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                    dscale, .false., rrhs)
      case (NDIM2D)
        call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                    ru, fcb_calcFluxTVD2d,&
                                    fcb_calcCharacteristics2d,&
                                    rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                    dscale, .false., rrhs)
      case (NDIM3D)
        call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                    ru, fcb_calcFluxTVD3d,&
                                    fcb_calcCharacteristics3d,&
                                    rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                    dscale, .false., rrhs)
      case DEFAULT
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'euler_calcPrimalRHS')
        call sys_halt()
      end select

      
      case DEFAULT
        call output_line('Invalid type of dissipation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalRHS')
        call sys_halt()
      end select

      ! Stop time measurement for residual/rhs evaluation
      call stat_stopTimer(rtimer_assembly_resrhs)

    end subroutine fcb_calcPrimalRHS

  !*****************************************************************************

!<subroutine>

    subroutine fcb_calcPrimalResidual(rproblemLevel, rtimestep, rsolver,&
                                      ru, ru0, rrhs, rres, ite, imode, istatus)

!<description>
    ! This subroutine computes the defect for 
    ! algebraic flux correction schemes
    !   $$ r = \Delta t L^n u^n $$
    ! and the constant right-hand 
    !   $$ b^n = M_L u^n+(1-\theta)\Delta t L^n u^n$$
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! number of nonlinear iteration
    integer, intent(IN) :: ite

    ! mode: (1) primal or (2) dual problem
    integer, intent(IN) :: imode
!</input>

!<inputoutput>
    ! multigrid level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! nonlinear solver structure
    type(t_solver), intent(INOUT) :: rsolver
    
    ! current solution vector
    type(t_vectorBlock), intent(INOUT) :: ru

    ! constant right-hand side
    type(t_vectorBlock), intent(INOUT) :: rrhs

    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! OPTIONAL: status of the callback function
    integer, intent(INOUT), optional :: istatus
!</inputoutput>
!</subroutine>

    ! local variable
    integer :: iblock

    ! We don't have to check vectors for compatibility.
    ! This will be done in the update procedure for the NSE-operator
    ! Update the global system operator if no Newton method is used
    call euler_calcPrimalPrecond(rproblemLevel, rtimestep, rsolver, ru, imode)

    ! Start time measurement for residual/rhs evaluation
    call stat_startTimer(rtimer_assembly_resrhs, STAT_TIMERSHORT)

    if (ite .eq. 0) then
      
      !-----------------------------------------------------
      ! In the first nonlinear iteration update the residual 
      ! $r=dt*L*U^n+f$ and the constant right hand side 
      ! $b=[M_L+(1-theta)*dt*L]*U^n$ are computed
      !-----------------------------------------------------
      
      ! Comput the initial residual vector
      select case(rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID)%ctypeAFCstabilisation)
        
      case (AFCSTAB_GALERKIN)

        ! Assemble high-order residual
        !
        !   res = dt*K*U^n
        !
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                   ru, fcb_calcFluxGalerkin1d,&
                                   rtimestep%dStep, .true., rres)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                   ru, fcb_calcFluxGalerkin2d,&
                                   rtimestep%dStep, .true., rres)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                   ru, fcb_calcFluxGalerkin3d,&
                                   rtimestep%dStep, .true., rres)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
          call sys_halt()
        end select

        
      case (AFCSTAB_UPWIND)

        ! Assemble low-order redisual
        !
        !   res = dt*LU^n
        !
        select case(rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID)%idissipation)
        case (AFCSTAB_SCALARDISSIPATION)
          
          ! Scalar dissipation
          !
          !   D = d_{ij}
          !
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                     ru, fcb_calcFluxScalarDiss1d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                     ru, fcb_calcFluxScalarDiss2d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                     ru, fcb_calcFluxScalarDiss3d,&
                                     rtimestep%dStep, .true., rres)
          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
            call sys_halt()
          end select

          
        case (AFCSTAB_TENSORDISSIPATION)
          
          ! Tensorial dissipation
          !
          !   D = D_{ij}
          !
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                     ru, fcb_calcFluxTensorDiss1d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                     ru, fcb_calcFluxTensorDiss2d,&
                                     rtimestep%dStep, .true., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                     ru, fcb_calcFluxTensorDiss3d,&
                                     rtimestep%dStep, .true., rres)
          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
            call sys_halt()
          end select
          
          
        case DEFAULT
          call output_line('Invalid type of dissipation!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
          call sys_halt()
        end select

        
      case (AFCSTAB_FEMFCT)

        ! Assemble low-order residual + FEM-FCT stabilization
        !
        !   res = dt*L(U)*U + F(U)
        
        ! Compute low-order residual with scalar dissipation
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                   ru, fcb_calcFluxScalarDiss1d,&
                                   (1._DP-rtimestep%dStep)*rtimestep%dStep,&
                                   .true., rrhs)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                   ru, fcb_calcFluxScalarDiss2d,&
                                   (1._DP-rtimestep%dStep)*rtimestep%dStep,&
                                   .true., rrhs)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                   ru, fcb_calcFluxScalarDiss3d,&
                                   (1._DP-rtimestep%dStep)*rtimestep%dStep,&
                                   .true., rrhs)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
          call sys_halt()
        end select
         

        ! Compute explicit low-order predictor
        !
        !   u^star = u + M_L^{-1})*[(1-theta)*dt*L(U)*U]
        !
        do iblock = 1, rrhs%nblocks
          call lsyssc_invertedDiagMatVec(rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                         rrhs%RvectorBlock(iblock), 1._DP,&
                                         rrhs%RvectorBlock(iblock))
        end do
        call lsysbl_vectorLinearComb(ru, rrhs, 1._DP, 1._DP)

        ! Initialize FCT limiter and compute flux correction
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidualFCT(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                      rproblemLevel%Rmatrix(CNSE_MATRIX_ML), ru,&
                                      fcb_calcFluxScalarDiss1d,&
                                      fcb_calcRawFluxFCT1d,&
                                      fcb_calcCharacteristics1d,&
                                      rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                      rtimestep%dStep, rtimestep%theta,&
                                      rtimestep%dStep, .true., rres, rrhs,&
                                      rproblemLevel%Rmatrix(CNSE_MATRIX_MC))
        case (NDIM2D)
          call gfsys_buildResidualFCT(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                      rproblemLevel%Rmatrix(CNSE_MATRIX_ML), ru,&
                                      fcb_calcFluxScalarDiss2d,&
                                      fcb_calcRawFluxFCT2d,&
                                      fcb_calcCharacteristics2d,&
                                      rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                      rtimestep%dStep, rtimestep%theta,&
                                      rtimestep%dStep, .true., rres, rrhs,&
                                      rproblemLevel%Rmatrix(CNSE_MATRIX_MC))
        case (NDIM3D)
          call gfsys_buildResidualFCT(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                      rproblemLevel%Rmatrix(CNSE_MATRIX_ML), ru,&
                                      fcb_calcFluxScalarDiss3d,&
                                      fcb_calcRawFluxFCT3d,&
                                      fcb_calcCharacteristics3d,&
                                      rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                      rtimestep%dStep, rtimestep%theta,&
                                      rtimestep%dStep, .true., rres, rrhs,&
                                      rproblemLevel%Rmatrix(CNSE_MATRIX_MC))
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
          call sys_halt()
        end select
        
        
      case (AFCSTAB_FEMTVD)
        
        ! Assemble low-order residual + FEM-TVD stabilization
        !
        !   res = dt*L(U)*U + F(U)
        !
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                      ru, fcb_calcFluxTVD1d,&
                                      fcb_calcCharacteristics1d,&
                                      rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                      rtimestep%dStep, .true., rres)
        case (NDIM2D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                      ru, fcb_calcFluxTVD2d,&
                                      fcb_calcCharacteristics2d,&
                                      rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                      rtimestep%dStep, .true., rres)
        case (NDIM3D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                      ru, fcb_calcFluxTVD3d,&
                                      fcb_calcCharacteristics3d,&
                                      rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                      rtimestep%dStep, .true., rres)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
          call sys_halt()
        end select

        
      case DEFAULT
        call output_line('Invalid type of stabilisation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
        call sys_halt()
      end select
      
      ! Compute the constant right-hand side for (pseudo-)transient flows
      ! 
      !   rhs = M_L*U^n + (1-theta)*res
      !
      select case(iflowtype)
      case (FLOW_TRANSIENT,&
            FLOW_PSEUDOTRANSIENT)
        
        do iblock = 1, ru%nblocks
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                   ru%RvectorBlock(iblock),&
                                   rrhs%RvectorBlock(iblock), 1._DP, 0._DP)
        end do
        call lsysbl_vectorLinearComb(rres, rrhs, (1._DP-rtimestep%theta), 1._DP)

      end select

    else   ! ite > 0
      
      !--------------------------------------------------------------
      ! In all subsequent nonlinear iterations update only the
      ! residual vector and make use of the constant r.h.s
      !--------------------------------------------------------------

      ! What type of flow are we?
      select case(iflowtype)
        
      case (FLOW_TRANSIENT,&
            FLOW_PSEUDOTRANSIENT)
        ! Adopt constant right-hand
        !
        !   res = rhs-ML*u
        !
        call lsysbl_vectorLinearComb(rrhs, rres, 1._DP, 0._DP)
        do iblock = 1, ru%nblocks
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                   ru%RvectorBlock(iblock),&
                                   rres%RvectorBlock(iblock), -1._DP, 1._DP)
        end do
        
      case (FLOW_STEADYSTATE)
        
        ! Clear residual vector
        !
        !   res = 0
        !
        call lsysbl_clearVector(rres)

      case DEFAULT
        call output_line('Invalid flow type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
      end select
      
      ! Update residual vector
      select case(rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID)%ctypeAFCstabilisation)

      case (AFCSTAB_GALERKIN)
        !
        ! Assemble high-order residual
        !
        !   res = res + dt*theta*K*U^(m)
        !
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                   ru, fcb_calcFluxGalerkin1d,&
                                   rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM2D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                   ru, fcb_calcFluxGalerkin2d,&
                                   rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM3D)
          call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                   ru, fcb_calcFluxGalerkin3d,&
                                   rtimestep%theta*rtimestep%dStep, .false., rres)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
          call sys_halt()
        end select
          
        
      case (AFCSTAB_UPWIND)
        !
        ! Assemble low-order redisual
        !
        !   res = res + dt*theta*L*U^(m)
        !
        select case(rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID)%idissipation)
        case (AFCSTAB_SCALARDISSIPATION)
          
          ! Scalar dissipation
          !
          !   D = d_{ij}
          !
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                     ru, fcb_calcFluxScalarDiss1d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                     ru, fcb_calcFluxScalarDiss2d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                     ru, fcb_calcFluxScalarDiss3d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
            call sys_halt()
          end select
          
          
        case (AFCSTAB_TENSORDISSIPATION)
          
          ! Tensorial dissipation
          !
          !   D = D_{ij}
          !
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                     ru, fcb_calcFluxTensorDiss1d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM2D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                     ru, fcb_calcFluxTensorDiss2d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case (NDIM3D)
            call gfsys_buildResidual(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                     ru, fcb_calcFluxTensorDiss3d,&
                                     rtimestep%theta*rtimestep%dStep, .false., rres)
          case DEFAULT
            call output_line('Invalid spatial dimension!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
            call sys_halt()
          end select
          
          
        case DEFAULT
          call output_line('Invalid type of dissipation!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
          call sys_halt()
        end select
        

      case (AFCSTAB_FEMFCT)
        
        ! Assemble low-order residual + FEM-FCT stabilization
        !
        !   res = res + F(U)
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidualFCT(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                      rproblemLevel%Rmatrix(CNSE_MATRIX_ML), ru,&
                                      fcb_calcFluxScalarDiss1d,&
                                      fcb_calcRawFluxFCT1d,&
                                      fcb_calcCharacteristics1d,&
                                      rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                      rtimestep%theta*rtimestep%dStep,&
                                      rtimestep%theta, rtimestep%dStep, .false.,&
                                      rres, rmatrixMC=rproblemLevel%Rmatrix(CNSE_MATRIX_MC))
        case (NDIM2D)
          call gfsys_buildResidualFCT(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                      rproblemLevel%Rmatrix(CNSE_MATRIX_ML), ru,&
                                      fcb_calcFluxScalarDiss2d,&
                                      fcb_calcRawFluxFCT2d,&
                                      fcb_calcCharacteristics2d,&
                                      rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                      rtimestep%theta*rtimestep%dStep,&
                                      rtimestep%theta, rtimestep%dStep, .false.,&
                                      rres, rmatrixMC=rproblemLevel%Rmatrix(CNSE_MATRIX_MC))
        case (NDIM3D)
          call gfsys_buildResidualFCT(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                      rproblemLevel%Rmatrix(CNSE_MATRIX_ML), ru,&
                                      fcb_calcFluxScalarDiss3d,&
                                      fcb_calcRawFluxFCT3d,&
                                      fcb_calcCharacteristics3d,&
                                      rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                      rtimestep%theta*rtimestep%dStep,&
                                      rtimestep%theta, rtimestep%dStep, .false.,&
                                      rres, rmatrixMC=rproblemLevel%Rmatrix(CNSE_MATRIX_MC))
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
          call sys_halt()
        end select
        

      case (AFCSTAB_FEMTVD)
        !
        ! Assemble low-order residual + FEM-TVD stabilization
        !
        !   res = res + dt*theta*L(U)*U + F(U)
        !
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                      ru, fcb_calcFluxTVD1d,&
                                      fcb_calcCharacteristics1d,&
                                      rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                      rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM2D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                      ru, fcb_calcFluxTVD2d,&
                                      fcb_calcCharacteristics2d,&
                                      rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                      rtimestep%theta*rtimestep%dStep, .false., rres)
        case (NDIM3D)
          call gfsys_buildResidualTVD(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                      ru, fcb_calcFluxTVD3d,&
                                      fcb_calcCharacteristics3d,&
                                      rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID),&
                                      rtimestep%theta*rtimestep%dStep, .false., rres)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
          call sys_halt()
        end select

        
      case DEFAULT
        call output_line('Invalid type of stabilisation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalResidual')
        call sys_halt()
      end select
      
    end if

    ! Stop time measurement for residual/rhs evaluation
    call stat_stopTimer(rtimer_assembly_resrhs)

  end subroutine fcb_calcPrimalResidual
  
  !*****************************************************************************

!<subroutine>

  subroutine fcb_calcPrimalJacobian(rproblemLevel, rtimestep, rsolver,&
                                    ru, ru0, imode, bfailure, istatus)

!<description>
    ! This callback subroutine computes the Jacobian matrix for the
    !  compressible Euler/Navier-Stokes equations
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: ru0

    ! mode: (1) primal or (2) dual problem
    integer, intent(IN) :: imode

    ! OPTIONAL: Newton subiteration failed, return to defect correction
    logical, intent(IN), optional :: bfailure
!</input>

!<inputoutput>
    ! multigrid level structure
    type(t_problemLevel), intent(INOUT), target :: rproblemLevel

    ! nonlinear solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! current solution vector
    type(t_vectorBlock), intent(INOUT) :: ru

    ! OPTIONAL: status of the callback function
    integer, intent(INOUT), optional :: istatus
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixBlock), pointer :: p_rmatrixA
    real(DP) :: hstep
    integer :: iblock,jblock
    
    ! Start time measurement for matrix evaluation
    call stat_startTimer(rtimer_assembly_matrix, STAT_TIMERSHORT)

    !----------------------------------------------------------------
    ! If Newton iteration failed completely, return to
    ! low-order preconditioner C=ML-theta*tstep*L(u)
    !----------------------------------------------------------------
    if (present(bfailure)) then
      if (bfailure) then
        ! Update the global system operator
        call euler_calcPrimalPrecond(rproblemLevel, rtimestep, rsolver, ru, imode)

        ! What type of flow are we?
        select case(iflowtype)
        case (FLOW_TRANSIENT,&
              FLOW_PSEUDOTRANSIENT)
          call euler_updateSolverMatrix(rproblemLevel, rsolver, CNSE_MATRIX_A,&
                                     UPDMAT_JAC_TRANSIENT, rproblemLevel%ilev, rproblemLevel%ilev)
          
        case (FLOW_STEADYSTATE)
          call euler_updateSolverMatrix(rproblemLevel, rsolver, CNSE_MATRIX_A,&
                                     UPDMAT_JAC_STEADY, rproblemLevel%ilev, rproblemLevel%ilev)

        case DEFAULT
          call output_line('Invalid flow type!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalJacobian')
        end select
        
        ! Update the content of the solver structure
        call solver_updateContent(rsolver)

        ! That's it, return
        return
      end if
    end if

    !---------------------------------------------------------------------------
    ! If Newton iteration should be performed, generate the Jacobian
    ! matrix and configure the matrices for the linear solver
    !---------------------------------------------------------------------------
    
    ! Compute step lenth
    select case(int(rsolver%p_solverNewton%dperturbationStrategy))
    case (PERTURB_NITSOL)
      ! Choice h=[(1+|u|)*EPS]^{1/(1+p)} by Pernice et al. in
      ! M. Pernice, H.F. Walker, NITSOL: a Newton iterative solver
      ! for nonlinear systems, SIAM J. Sci. Comput. 19 (1998) 302-318.
      hstep = ( (1+&
          lsysbl_vectorNorm(ru, LINALG_NORMEUCLID))*SYS_EPSREAL )**(1._DP/3._DP)
      
    case (PERTURB_SQRTEPS)
      hstep= sqrt(SYS_EPSREAL)

    case DEFAULT
      hstep=max(SYS_EPSREAL, rsolver%p_solverNewton%dperturbationStrategy)
    end select

    ! Check if vector has correct dimensions
    select case(isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)
      if ((ru%nblocks  .ne. 1) .or.&
          (ru0%nblocks .ne. 1) .or.&
          (ru%RvectorBlock(1)%NVAR  .ne. euler_getNVAR(rproblemLevel)) .or. &
          (ru0%RvectorBlock(1)%NVAR .ne. euler_getNVAR(rproblemLevel))) then
        call output_line('Invalid format of solution vector!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcJacobian')
        call sys_halt()
      end if

    case (SYSTEM_BLOCKFORMAT)
      if ((ru%nblocks  .ne. euler_getNVAR(rproblemLevel)) .or.&
          (ru0%nblocks .ne. euler_getNVAR(rproblemLevel))) then
        call output_line('Invalid format of solution vector!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcJacobian')
        call sys_halt()
      else
        do iblock = 1, ru%nblocks
          if ((ru%RvectorBlock(iblock)%NVAR  .ne. 1) .or.&
              (ru0%RvectorBlock(iblock)%NVAR .ne. 1)) then
            call output_line('Invalid format of solution vector!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcJacobian')
            call sys_halt()
          end if
        end do
      end if

    case DEFAULT
      call output_line('Unsupported system format!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcJacobian')
      call sys_halt()
    end select

    ! Set pointer
    p_rmatrixA => rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)

    ! Check if discrete operator needs to be initialized
    if (rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID)%iSpec .eq. AFCSTAB_UNDEFINED)&
        call gfsys_initStabilisation(p_rmatrixA,&
                                     rproblemLevel%Rafcstab(CNSE_AFCSTAB_INVISCID))

    !---------------------------------------------------------------------------
    ! Assemble inviscid part
    !---------------------------------------------------------------------------
    
    ! What kind of coupling is applied?
    select case(icoupled)
      
    ! Assemble block-diagonal divergence operator
    case (FLOW_SEGREGATED)
      
      ! What kind of preconditioner is applied?
      select case(iprecond)
      case (AFCSTAB_GALERKIN)
        ! Assemble Jacobian operator for standard Galerin scheme
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                      ru, fcb_calcMatrixGalerkinDiag1d,&
                                      hstep, .true., p_rmatrixA)
        case (NDIM2D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                      ru, fcb_calcMatrixGalerkinDiag2d,&
                                      hstep, .true., p_rmatrixA)
        case (NDIM3D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                      ru, fcb_calcMatrixGalerkinDiag3d,&
                                      hstep, .true., p_rmatrixA)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalJacobian')
          call sys_halt()
        end select

        
      case (AFCSTAB_SCALARDISSIPATION)
        ! Assemble Jacobian operator with scalar dissipation
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                      ru, fcb_calcMatrixScalarDissDiag1d,&
                                      hstep, .true., p_rmatrixA)
        case (NDIM2D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                      ru, fcb_calcMatrixScalarDissDiag2d,&
                                      hstep, .true., p_rmatrixA)
        case (NDIM3D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                      ru, fcb_calcMatrixScalarDissDiag3d,&
                                      hstep, .true., p_rmatrixA)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalJacobian')
          call sys_halt()
        end select
        
        
        
      case (AFCSTAB_TENSORDISSIPATION)
        ! Assemble Jacobian operator with tensorial dissipation
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                      ru, fcb_calcMatrixTensorDissDiag1d,&
                                      hstep, .true., p_rmatrixA)
        case (NDIM2D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                      ru, fcb_calcMatrixTensorDissDiag2d,&
                                      hstep, .true., p_rmatrixA)
        case (NDIM3D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                      ru, fcb_calcMatrixTensorDissDiag3d,&
                                      hstep, .true., p_rmatrixA)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalJacobian')
          call sys_halt()
        end select


      case DEFAULT
        call output_line('Invalid type of dissipation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalJacobian')
        call sys_halt()
      end select

        
    ! Assemble full block transport operator
    case (FLOW_ALLCOUPLED)
      
      ! What kind of preconditioner is applied?
      select case(iprecond)         
      case (AFCSTAB_GALERKIN)
        ! Assemble Jacobian operator for standard Galerin scheme
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                      ru, fcb_calcMatrixGalerkin1d,&
                                     hstep, .true., p_rmatrixA)
        case (NDIM2D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                      ru, fcb_calcMatrixGalerkin2d,&
                                     hstep, .true., p_rmatrixA)
        case (NDIM3D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                      ru, fcb_calcMatrixGalerkin3d,&
                                      hstep, .true., p_rmatrixA)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalJacobian')
          call sys_halt()
        end select

        
      case (AFCSTAB_SCALARDISSIPATION)
        ! Assemble Jacobian operator with scalar dissipation
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                      ru, fcb_calcMatrixScalarDiss1d,&
                                      hstep, .true., p_rmatrixA)
        case (NDIM2D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                      ru, fcb_calcMatrixScalarDiss2d,&
                                      hstep, .true., p_rmatrixA)
        case (NDIM3D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                      ru, fcb_calcMatrixScalarDiss3d,&
                                      hstep, .true., p_rmatrixA)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalJacobian')
          call sys_halt()
        end select


      case (AFCSTAB_TENSORDISSIPATION)
        ! Assemble Jacobian operator with tensorial dissipation
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CX),&
                                      ru, fcb_calcMatrixTensorDiss1d,&
                                      hstep, .true., p_rmatrixA)
        case (NDIM2D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CY),&
                                      ru, fcb_calcMatrixTensorDiss2d,&
                                      hstep, .true., p_rmatrixA)
        case (NDIM3D)
          call gfsys_buildDivJacobian(rproblemLevel%Rmatrix(CNSE_MATRIX_CX:CNSE_MATRIX_CZ),&
                                      ru, fcb_calcMatrixTensorDiss3d,&
                                      hstep, .true., p_rmatrixA)
        case DEFAULT
          call output_line('Invalid spatial dimension!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalJacobian')
          call sys_halt()
        end select


      case DEFAULT
        call output_line('Invalid type of dissipation!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalJacobian')
        call sys_halt()
      end select


    case DEFAULT
      call output_line('Invalid flow coupling!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalJacobian')
      call sys_halt()
    end select

    
    !---------------------------------------------------------------------------
    ! Assemble the global Jacobian operator
    !---------------------------------------------------------------------------
    
    if (isystemFormat .eq. SYSTEM_INTERLEAVEFORMAT) then
      
      ! What type of flow are we?
      select case(iflowtype)
      case (FLOW_TRANSIENT,&
            FLOW_PSEUDOTRANSIENT)
        ! Compute the global operator for transient flows
        !
        !     A = blockdiag(M_L)-theta*dt*J
        !
        if (rtimestep%theta .gt. 0.0_DP) then
          call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(CNSE_MATRIX_ML), 1._DP,&
                                       rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                       -rtimestep%theta*rtimestep%dStep,&
                                       rproblemLevel%Rmatrix(CNSE_MATRIX_A),&
                                       .false., .false., .true., .true.)
        else
          call lsyssc_spreadMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_ML),&
                                   rproblemLevel%Rmatrix(CNSE_MATRIX_A))
        end if
        
      case (FLOW_STEADYSTATE)
        ! Compute the global operator for steady-state flow
        !
        !     A = -J
        !
        call lsyssc_scaleMatrix(rproblemLevel%Rmatrix(CNSE_MATRIX_A), -1._DP)
        
      case DEFAULT
        call output_line('Invalid flow type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'calcPrimalJacobian')
        call sys_halt()
      end select
      
    else
      
      ! What type of flow are we?
      select case(iflowtype)
      case (FLOW_TRANSIENT,&
            FLOW_PSEUDOTRANSIENT)
        ! Compute the global operator for transient flows
        !
        !     A = diag(M_L)-theta*dt*J
        !
        do iblock = 1, euler_getNVAR(rproblemLevel)
          call lsyssc_MatrixLinearComb(&
              rproblemLevel%Rmatrix(CNSE_MATRIX_ML), 1._DP,&
              rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(iblock,iblock),&
              -rtimestep%theta*rtimestep%dStep,&
              rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(iblock,iblock),&
              .false., .false., .true., .true.)
        end do
        
        ! Scale the off-diagonal blocks by "-theta*dt" if required
        if (icoupled .eq. FLOW_ALLCOUPLED) then
          do iblock = 1, euler_getNVAR(rproblemLevel)
            do jblock = 1, euler_getNVAR(rproblemLevel)
              if (iblock .eq. jblock) cycle
              call lsyssc_scaleMatrix(&
                  rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(iblock,jblock),&
                  -rtimestep%theta*rtimestep%dStep)
            end do
          end do
        end if
        
      case (FLOW_STEADYSTATE)
        ! Compute the global operator for steady-state flow
        !
        !     A = -J
        !
        select case(icoupled)
          
        case (FLOW_SEGREGATED)
          do iblock = 1, euler_getNVAR(rproblemLevel)
            call lsyssc_scaleMatrix(&
                rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(iblock,iblock), -1.0_DP)
          end do
          
        case (FLOW_ALLCOUPLED)
          do iblock = 1, euler_getNVAR(rproblemLevel)
            do jblock = 1, euler_getNVAR(rproblemLevel)
              call lsyssc_scaleMatrix(&
                  rproblemLevel%RmatrixBlock(CNSE_MATRIX_A)%RmatrixBlock(iblock,jblock), -1.0_DP)
            end do
          enddo
          
        case DEFAULT
          call output_line('Unsupported block matrix format!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalJacobian')
          call sys_halt()
        end select
        
      case DEFAULT
        call output_line('Invalid flow type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'calcPrimalJacobian')
        call sys_halt()
      end select
      
    end if

    ! !!!!
    ! NOTE
    ! !!!!

    ! If one day the Jacobian matrix for the stabilization can be assembled
    ! then it will go here, and moreover, the extended sparsity pattern
    ! has to be generated and adopted for matrix J.   
   
    !---------------------------------------------------------------------------
    ! Impose boundary conditions
    !---------------------------------------------------------------------------
    
    call bdrf_filterMatrix(rsolver%rboundaryCondition, rproblemLevel%rtriangulation,&
                           rproblemLevel%Rmatrix(CNSE_MATRIX_J), 1.0_DP)

    ! Ok, we updated the Jacobian matrix successfully. Now we stil have to
    ! link it to the solver hierarch. This is done recursively.

    ! What type of flow are we?
    select case(iflowtype)
    case (FLOW_TRANSIENT,&
          FLOW_PSEUDOTRANSIENT)
      call euler_updateSolverMatrix(rproblemLevel, rsolver, CNSE_MATRIX_J,&
                                 UPDMAT_JAC_TRANSIENT, rproblemLevel%ilev, rproblemLevel%ilev)

    case (FLOW_STEADYSTATE)
      call euler_updateSolverMatrix(rproblemLevel, rsolver, CNSE_MATRIX_J, &
                                 UPDMAT_JAC_STEADY, rproblemLevel%ilev, rproblemLevel%ilev)

    case DEFAULT
      call output_line('Invalid flow type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'fcb_calcPrimalJacobian')
      call sys_halt()
    end select

    ! Finally, we have to update the content of the solver hierarchy
    call solver_updateContent(rsolver)

    ! Stop time measurement for matrix evaluation
    call stat_stopTimer(rtimer_assembly_matrix)

  end subroutine fcb_calcPrimalJacobian

  !*****************************************************************************

!<subroutine>

  subroutine fcb_applyPrimalJacobian(rproblemLevel, rx, ry, cx, cy, istatus)

!<description>
    ! This subroutine applies the (scaled) Jacobian matrix 
    ! to the vector rx and adds the result to the vector ry
!</description>

!<input>
    ! multigrid level structure
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

    ! OPTIONAL: status of the callback function
    integer, intent(INOUT), optional :: istatus
!</inputoutput>
!</subroutine>

    ! We know where the Jacobian matrix is stored and can apply it
    ! by means of matrix vector multiplication
    call lsysbl_blockMatVec(rproblemLevel%RmatrixBlock(CNSE_MATRIX_A), rx, ry, cx, cy)

  end subroutine fcb_applyPrimalJacobian

end module euler_primalproblem

!##############################################################################
!# ****************************************************************************
!# <name> codire_problem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to solve the primal
!# problem for a scalar conservation law in arbitrary spatial dimensions.
!#
!# The following routines are available:
!#
!# 1.) codire_calcPrecond
!#     -> calculate the preconditioner
!#
!# 2.) codire_calcRHS
!#     -> calculate the right-hand side vector
!#
!# 4.) codire_calcJacobian
!#     -> calculate the Jacobian matrix
!#
!# 5.) codire_applyJacobian
!#     -> apply the Jacobian matrix
!#
!# 6.) codire_calcGalerkinResidual
!#     -> calculate the Galerkin residual
!#
!# The following auxiliary routine are available:
!#
!# 1.) codire_setVelocity
!#     -> set the velocity vector internally
!#
!# 2.) codire_updateVelocity
!#     -> update the velocity vector
!#
!# 3.) codire_calcConvectionConst1d
!#     -> compute convection coefficients for constant velocity in 1D
!#
!# 4.) codire_calcConvectionConst2d
!#     -> compute convection coefficients for constant velocity in 2D
!#
!# 5.) codire_calcConvectionConst3d
!#     -> compute convection coefficients for constant velocity in 3D
!#
!# 6.) codire_calcConvectionBurgersSpT2d
!#     -> compute convection coefficients for 1D Burgers' equation in space-time
!#
!# 7.) codire_calcConvection_BuckLevSpT2d
!#     -> compute convection coefficients for 1D Buckley-Leverett equation in space-time
!#
!# 8.) codire_calcConvection_Burgers1d
!#     -> compute convection coefficients for Burgers' equation in 1D
!#
!# 9.) codire_calcConvection_Burgers2d
!#     -> compute convection coefficients for Burgers' equation in quasi-2D
!#
!# 10.) codire_calcConvection_BuckLev1d
!#     -> compute convection coefficients for Buckley-Leverett equation in 1D
!#
!# </purpose>
!##############################################################################

module codire_problem

  use afcstabilisation
  use boundaryfilter
  use codire_basic
  use collection
  use fparser
  use fsystem
  use genoutput
  use groupfemscalar
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solver
  use statistics
  use storage

  implicit none

  private

!  public :: codire_calcPrecond
!  public :: codire_calcRHS
!  public :: codire_calcJacobian
!  public :: codire_applyJacobian
!  public :: codire_calcGalerkinResidual

!<globals>

!!$  !*****************************************************************
!!$  ! Pointers to X-, Y- and Z-component of velocity vector
!!$  real(DP), dimension(:), pointer, save :: p_DvelocityX => null()
!!$  real(DP), dimension(:), pointer, save :: p_DvelocityY => null()
!!$  real(DP), dimension(:), pointer, save :: p_DvelocityZ => null()

!</globals>

  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

contains
  
!!$  !*****************************************************************************
!!$
!!$!<subroutine>
!!$
!!$  subroutine codire_calcPrecond(rproblemLevel, rtimestep, rsolver, ru, rcollection)
!!$
!!$!<description>
!!$    ! This subroutine calculates the nonlinear preconditioner and
!!$    ! configures the linear solver structure accordingly. Depending on
!!$    ! the nonlinear solver, the low-order evolution operator or the
!!$    ! Jacobian operator is adopted as nonlinear preconditioner.
!!$!</description>
!!$
!!$!<input>
!!$    ! time-stepping structure
!!$    type(t_timestep), intent(IN) :: rtimestep
!!$
!!$    ! solution vector
!!$    type(t_vectorBlock), intent(IN) :: ru
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! problem level structure
!!$    type(t_problemLevel), intent(INOUT) :: rproblemLevel
!!$
!!$    ! solver structure
!!$    type(t_solver), intent(INOUT) :: rsolver
!!$
!!$    ! collection
!!$    type(t_collection), intent(INOUT) :: rcollection
!!$!</inputoutput>
!!$!</subroutine>
!!$
!!$    ! local variables
!!$    type(t_problemLevel), pointer :: rproblemLevelTmp
!!$    type(t_vectorBlock), pointer :: p_ru
!!$    type(t_timer), pointer :: rtimer
!!$    type(t_vectorBlock), target :: ruFine
!!$    type(t_vectorBlock) :: ruCoarse
!!$    integer :: nlmin,nlmax
!!$    logical :: bStabilize,bcompatible,bisDivFree
!!$
!!$    integer :: systemMatrix
!!$    integer :: lumpedMassMatrix, consistentMassMatrix
!!$    integer :: transportMatrix
!!$    integer :: coeffMatrix_CX,coeffMatrix_CY,coeffMatrix_CZ
!!$    integer :: coeffMatrix_S
!!$    integer :: imasstype,iconvectiontype,idiffusiontype
!!$    integer :: convectionAFC, diffusionAFC
!!$
!!$    integer :: imode = 1
!!$
!!$    
!!$    ! Start time measurement for global operator
!!$    rtimer => collct_getvalue_timer(rcollection, 'timerAssemblyMatrix')
!!$    call stat_startTimer(rtimer, STAT_TIMERSHORT)
!!$
!!$stop
!!$    
!!$    ! Update velocity vector and return if velocity has not changed. This can be
!!$    ! done, because the system operator can only change due to velocity variations.
!!$    if (.not.codire_updateVelocity(rproblemLevel%p_rproblem, rtimestep%dTime,&
!!$                                   rproblemLevel%ilev, rproblemLevel%ilev)) return
!!$    
!!$    ! Get parameters from collection which are required unconditionally
!!$    transportMatrix = collct_getvalue_int(rcollection, 'transportmatrix')
!!$    coeffMatrix_S   = collct_getvalue_int(rcollection, 'coeffMatrix_S')
!!$
!!$    !---------------------------------------------------------------------------
!!$    ! Assemble diffusion operator
!!$    !---------------------------------------------------------------------------
!!$    
!!$    idiffusiontype = collct_getvalue_int(rcollection, 'idiffusiontype')
!!$    
!!$    select case(idiffusiontype)
!!$    case (DIFFUSION_ZERO)
!!$      ! zero diffusion, clear the transport matrix
!!$      call lsyssc_clearMatrix(rproblemLevel%Rmatrix(transportMatrix))
!!$      
!!$    case (DIFFUSION_ISOTROPIC)
!!$      ! Isotropic diffusion
!!$      call gfsc_buildDiffusionOperator(rproblemLevel%Rmatrix(coeffMatrix_S),&
!!$                                       1._DP , .false., .true.,&
!!$                                       rmatrixDest=rproblemLevel%Rmatrix(transportMatrix))
!!$    
!!$    case (DIFFUSION_ANISOTROPIC)
!!$      ! Anisotropic diffusion
!!$      diffusionAFC = collct_getvalue_int(rcollection, 'diffusionAFC')
!!$
!!$      ! Check if stabilisation structure is initialised
!!$      if (rproblemLevel%Rafcstab(diffusionAFC)%iSpec .eq. AFCSTAB_UNDEFINED) then
!!$        call gfsc_initStabilisation(rproblemLevel%Rmatrix(coeffMatrix_S),&
!!$                                    rproblemLevel%Rafcstab(diffusionAFC))
!!$      end if
!!$      
!!$      ! What kind of stabilisation should be applied?
!!$      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
!!$      case (AFCSTAB_SYMMETRIC)
!!$        call gfsc_buildDiffusionOperator(rproblemLevel%Rmatrix(coeffMatrix_S),&
!!$                                         1._DP, .true., .true.,&
!!$                                         rproblemLevel%Rafcstab(diffusionAFC),&
!!$                                         rproblemLevel%Rmatrix(transportMatrix))
!!$      case DEFAULT
!!$        call gfsc_buildDiffusionOperator(rproblemLevel%Rmatrix(coeffMatrix_S),&
!!$                                         1._DP, .false., .true.,&
!!$                                         rmatrixDest=rproblemLevel%Rmatrix(transportMatrix))
!!$      end select
!!$
!!$    case (DIFFUSION_VARIABLE)
!!$      print *, "Variable diffusion matrices are yet not implemented!"
!!$      stop
!!$      
!!$    case DEFAULT
!!$      call output_line('Invalid type of diffusion!',&
!!$                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcPrecond')
!!$      call sys_halt()
!!$    end select
!!$    
!!$    
!!$    !---------------------------------------------------------------------------
!!$    ! Assemble convective operator
!!$    !---------------------------------------------------------------------------
!!$      
!!$    ! Set velocity vector for current level
!!$    call codire_setVelocity(rproblemLevel)
!!$    
!!$    ! Check if stabilisation structure is initialised
!!$    if (rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%iSpec .eq. AFCSTAB_UNDEFINED)&
!!$        call gfsc_initStabilisation(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
!!$                                    rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$
!!$    ! Do we have to perform stabilization of AFC-type?
!!$    bStabilize = (AFCSTAB_GALERKIN .ne.&
!!$                  rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation)
!!$    
!!$    ! Check if velocity is assumed to be discretely divergence free
!!$    bisDivFree = (ivelocitytype .gt. 0)
!!$    
!!$    if (bStabilize) then
!!$      
!!$      ! Assemble stabilized transport operator on finest level
!!$      select case(abs(ivelocitytype))
!!$      case (VELOCITY_ZERO)
!!$        ! zero velocity, do nothing
!!$        
!!$      case (VELOCITY_CONSTANT,&
!!$            VELOCITY_TIMEDEP)
!!$        ! linear velocity
!!$        select case(rproblemLevel%rtriangulation%ndim)
!!$        case (NDIM1D)
!!$          select case(imode)
!!$          case (1)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                              ru, codire_calcPrimalConvConst1d, bisDivFree,&
!!$                                              bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                              rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$          case (2)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                              ru, codire_calcDualConvConst1d, bisDivFree,&
!!$                                              bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                              rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$          end select
!!$          
!!$        case (NDIM2D)
!!$          select case(imode)
!!$          case (1)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                              ru, codire_calcPrimalConvConst2d, bisDivFree,&
!!$                                              bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                              rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$          case (2)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                              ru, codire_calcDualConvConst2d, bisDivFree,&
!!$                                              bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                              rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$          end select
!!$
!!$        case (NDIM3D)
!!$          select case(imode)
!!$          case (1)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CZ),&
!!$                                              ru, codire_calcPrimalConvConst3d, bisDivFree,&
!!$                                              bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                              rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$          case (2)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CZ),&
!!$                                              ru, codire_calcDualConvConst3d, bisDivFree,&
!!$                                              bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                              rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$          end select
!!$          
!!$        case DEFAULT
!!$          call output_line('Invalid spatial dimension!',&
!!$                           OU_CLASS_ERROR,OU_MODE_STD,'codire_calcPrecond')
!!$          call sys_halt()
!!$        end select
!!$        
!!$
!!$      case (VELOCITY_BURGERS_SPACETIME)
!!$        ! nonlinear Burgers' equation in space-time
!!$        call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                          ru, codire_calcConvectionBurgersSpT2d, bisDivFree,&
!!$                                          bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                          rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$        
!!$      case (VELOCITY_BUCKLEV_SPACETIME)
!!$        ! nonlinear Buckley-Leverett equation in space-time
!!$        call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                          ru, codire_calcConvectionBuckLevSpT2d, bisDivFree,&
!!$                                          bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                          rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$        
!!$      case (VELOCITY_BURGERS1D)
!!$        ! nonlinear Burgers' equation in 1D
!!$        call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                          ru, codire_calcConvectionBurgers1d, bisDivFree,&
!!$                                          bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                          rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$        
!!$      case (VELOCITY_BURGERS2D)
!!$        ! nonlinear Burgers' equation in 2D
!!$        call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                          ru, codire_calcConvectionBurgers2d, bisDivFree,&
!!$                                          bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                          rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$
!!$      case (VELOCITY_BUCKLEV1D)
!!$        ! nonlinear Buckley-Leverett equation in 1D
!!$        call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                          ru, codire_calcConvectionBuckLev1d, bisDivFree,&
!!$                                          bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                          rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$
!!$      case DEFAULT
!!$        call output_line('Invalid velocity profile!',&
!!$                         OU_CLASS_ERROR,OU_MODE_STD,'codire_calcPrecond')
!!$        call sys_halt()
!!$      end select
!!$        
!!$    else
!!$      
!!$      ! Assemble high-/low-order transport operator on coarser levels
!!$      select case(abs(ivelocitytype))
!!$      case (VELOCITY_ZERO)
!!$        ! zero velocity, do nothing
!!$        
!!$      case (VELOCITY_CONSTANT,&
!!$            VELOCITY_TIMEDEP)
!!$        ! linear velocity
!!$        select case(rproblemLevelTmp%rtriangulation%ndim)
!!$        case (NDIM1D)
!!$          select case(imode)
!!$          case (1)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                              ru, codire_calcPrimalConvConst1d, bisDivFree,&
!!$                                              bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$          case (2)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                              ru, codire_calcDualConvConst1d, bisDivFree,&
!!$                                              bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$          end select
!!$          
!!$        case (NDIM2D)
!!$          select case(imode)
!!$          case (1)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                              ru, codire_calcPrimalConvConst2d, bisDivFree,&
!!$                                              bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$          case (2)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                              ru, codire_calcDualConvConst2d, bisDivFree,&
!!$                                              bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$          end select
!!$          
!!$        case (NDIM3D)
!!$          select case(imode)
!!$          case (1)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CZ),&
!!$                                              ru, codire_calcPrimalConvConst3d, bisDivFree,&
!!$                                              bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$          case (2)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CZ),&
!!$                                              ru, codire_calcDualConvConst3d, bisDivFree,&
!!$                                              bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$          end select
!!$
!!$        case DEFAULT
!!$          call output_line('Invalid spatial dimension!',&
!!$                           OU_CLASS_ERROR,OU_MODE_STD,'codire_calcPrecond')
!!$          call sys_halt()
!!$        end select
!!$        
!!$          
!!$      case (VELOCITY_BURGERS_SPACETIME)
!!$        ! nonlinear Burgers' equation in space-time
!!$        call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                          ru, codire_calcConvectionBurgersSpT2d, bisDivFree,&
!!$                                          bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$          
!!$      case (VELOCITY_BUCKLEV_SPACETIME)
!!$        ! nonlinear Buckley-Leverett equation in space-time
!!$        call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                          ru, codire_calcConvectionBuckLevSpT2d, bisDivFree,&
!!$                                          bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$
!!$      case (VELOCITY_BURGERS1D)
!!$        ! nonlinear Burgers' equation in 1D
!!$        call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                          ru, codire_calcConvectionBurgers1d, bisDivFree,&
!!$                                          bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$        
!!$      case (VELOCITY_BURGERS2D)
!!$        ! nonlinear Burgers' equation in 2D
!!$        call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                          ru, codire_calcConvectionBurgers2d, bisDivFree,&
!!$                                          bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$        
!!$      case (VELOCITY_BUCKLEV1D)
!!$        ! nonlinear Buckley-Leverett equation in 1D
!!$        call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                          ru, codire_calcConvectionBuckLev1d, bisDivFree,&
!!$                                          bStabilize, .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$        
!!$      case DEFAULT
!!$        call output_line('Invalid velocity profile!',&
!!$                         OU_CLASS_ERROR,OU_MODE_STD,'codire_calcPrecond1d')
!!$        call sys_halt()
!!$      end select
!!$    end if
!!$    
!!$    !---------------------------------------------------------------------------
!!$    ! Assemble the global system operator
!!$    !---------------------------------------------------------------------------
!!$    
!!$    imasstype = collct_getvalue_int(rcollection, 'imasstype')
!!$
!!$    select case(imasstype)
!!$    case (MASS_LUMPED)
!!$      ! Compute the global operator for transient flow
!!$      !
!!$      !   $ A = ML-theta*dt*L $
!!$
!!$      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(lumpedMassMatrix), 1._DP,&
!!$                                   rproblemLevel%Rmatrix(transportMatrix),&
!!$                                   -rtimestep%theta*rtimestep%dStep,&
!!$                                   rproblemLevel%Rmatrix(systemMatrix),&
!!$                                   .false., .false., .true., .true.)
!!$    case (MASS_CONSISTENT)
!!$      ! Compute the global operator for transient flow
!!$      !
!!$      !   $ A = MC-theta*dt*L $
!!$
!!$      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(consistentMassMatrix), 1._DP,&
!!$                                   rproblemLevel%Rmatrix(transportMatrix),&
!!$                                   -rtimestep%theta*rtimestep%dStep,&
!!$                                   rproblemLevel%Rmatrix(systemMatrix),&
!!$                                   .false., .false., .true., .true.)
!!$
!!$    case DEFAULT
!!$      ! Compute the global operator for steady-state flow
!!$      !
!!$      !     A = -L
!!$      !
!!$      call lsyssc_copyMatrix(rproblemLevel%Rmatrix(transportMatrix),&
!!$                             rproblemLevel%Rmatrix(systemMatrix))
!!$      call lsyssc_scaleMatrix(rproblemLevel%Rmatrix(systemMatrix), -1._DP)
!!$    
!!$    end select
!!$    
!!$    !---------------------------------------------------------------------------
!!$    ! Impose boundary conditions
!!$    !---------------------------------------------------------------------------
!!$      
!!$    call bdrf_filterMatrix(rsolver%rboundaryCondition,&
!!$                           rproblemLevel%rtriangulation,&
!!$                           rproblemLevel%Rmatrix(systemMatrix), 1.0_DP)
!!$        
!!$    ! Ok, we updated the (nonlinear) system operator successfully. Now we still 
!!$    ! have to link it to the solver hierarchy. This is done recursively.
!!$    call codire_updateSolverMatrix(rproblemLevel, rsolver, systemMatrix, UPDMAT_ALL,&
!!$                                   rproblemLevel%ilev, rproblemLevel%ilev)
!!$
!!$    ! Finally, we have to update the content of the solver hierarchy
!!$    call solver_updateContent(rsolver)
!!$
!!$    ! Stop time measurement for global operator
!!$    call stat_stopTimer(rtimer)
!!$    
!!$  end subroutine codire_calcPrecond

!!$  !*****************************************************************************
!!$  
!!$!<subroutine>
!!$
!!$  subroutine codire_calcRHS(rproblemLevel, rtimestep, rsolver,&
!!$                            ru, ru0, rrhs, istep, rcollection)
!!$
!!$!<description>
!!$    ! This subroutine computes the right-hand side vector to be
!!$    ! used in the explicit Lax-Wendroff time stepping scheme.
!!$!</description>
!!$
!!$!<input>
!!$    ! time-stepping structure
!!$    type(t_timestep), intent(IN) :: rtimestep
!!$
!!$    ! initial solution vector
!!$    type(t_vectorBlock), intent(IN) :: ru0
!!$
!!$    ! number of explicit step
!!$    integer, intent(IN) :: istep
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! problem level structure
!!$    type(t_problemLevel), intent(INOUT) :: rproblemLevel
!!$
!!$    ! solver structure
!!$    type(t_solver), intent(INOUT) :: rsolver
!!$    
!!$    ! solution vector
!!$    type(t_vectorBlock), intent(INOUT) :: ru
!!$
!!$    ! right-hand side vector
!!$    type(t_vectorBlock), intent(INOUT) :: rrhs
!!$
!!$    ! collection
!!$    type(t_collection), intent(INOUT) :: rcollection
!!$!</inputoutput>
!!$!</subroutine>
!!$    
!!$    ! local variables
!!$    real(DP) :: dscale
!!$    logical :: bafc
!!$
!!$    ! Start time measurement for residual/rhs evaluation
!!$    call stat_startTimer(rtimer_assembly_resrhs, STAT_TIMERSHORT)
!!$
!!$    ! For nonlinear conservation laws, the global system operator
!!$    ! needs to be updated in each nonlinear iteration. The update 
!!$    ! routine is written such that it first determines if the problem
!!$    ! is nonlinear and returns without matrix update otherwise.
!!$    call codire_calcPrecond(rproblemLevel, rtimestep, rsolver, ru, rcollection)
!!$
!!$    select case(iflowtype)
!!$    case (FLOW_TRANSIENT,&
!!$          FLOW_PSEUDOTRANSIENT)
!!$
!!$       ! Compute the right-hand side
!!$       !
!!$       !     rhs = weight*(1-theta)*dt*L(u)*u
!!$       !
!!$       ! where weight is given by the time-stepping scheme
!!$       dscale = rtimestep%DmultistepWeights(istep)*&
!!$                rtimestep%dStep*(1._DP-rtimestep%theta)
!!$       call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                ru%rvectorBlock(1), rrhs%RvectorBlock(1), dscale, 0._DP)
!!$
!!$    case (FLOW_STEADYSTATE)
!!$
!!$       ! Compute the right-hand side
!!$       !
!!$       !     rhs = weight*(1-theta)*dt*L(u)*u
!!$       !
!!$       ! where weight is given by the time-stepping scheme
!!$       dscale = rtimestep%DmultistepWeights(istep)*&
!!$                rtimestep%dStep*(1._DP-rtimestep%theta)
!!$       call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(CDEQ_MATRIX_L), ru%rvectorBlock(1),&
!!$                                rrhs%RvectorBlock(1), dscale, 0._DP)
!!$       
!!$    case DEFAULT
!!$       call output_line('Invalid flow type!',&
!!$                        OU_CLASS_ERROR,OU_MODE_STD,'codire_calcRHS')
!!$       call sys_halt()
!!$    end select
!!$    
!!$
!!$    ! Apply stabilization of AFC type for convective part?
!!$    bafc = (AFCSTAB_GALERKIN .ne.&
!!$            rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation .and.&
!!$            AFCSTAB_UPWIND   .ne.&
!!$            rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation)
!!$
!!$    ! Perform algebraic flux correction (if required)
!!$    !
!!$    !     rhs = rhs + f^*(u^n+1,u^n)
!!$    !
!!$    if (bafc) then
!!$       
!!$       ! What kind of stabilisation should be applied?
!!$       select case(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation)
!!$          
!!$       case (AFCSTAB_FEMFCT,&
!!$             AFCSTAB_FEMFCT_EXP)
!!$          dscale = rtimestep%DmultistepWeights(istep)*rtimestep%dStep
!!$          call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
!!$                                     rproblemLevel%Rmatrix(CDEQ_MATRIX_ML),&
!!$                                     ru, rtimestep%theta, dscale, .true., rrhs,&
!!$                                     rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$          
!!$       case (AFCSTAB_FEMTVD,&
!!$             AFCSTAB_FEMGP)
!!$          call gfsc_buildResidualTVD(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),&
!!$                                     ru, ru0, rtimestep%theta, dscale, rrhs,&
!!$                                     rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION))
!!$       end select
!!$    end if
!!$
!!$    
!!$    ! Apply stabilization of AFC type for diffusive part?
!!$    bafc = (AFCSTAB_GALERKIN .ne.&
!!$            rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%ctypeAFCstabilisation .and.&
!!$            AFCSTAB_UPWIND   .ne.&
!!$            rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%ctypeAFCstabilisation)
!!$
!!$    ! Perform algebraic flux correction (if required)
!!$    !
!!$    !     rhs = rhs + g^*(u^n+1,u^n)
!!$    !
!!$    if (bafc) then
!!$       
!!$       ! What kind of stabilisation should be applied?
!!$       select case(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%ctypeAFCstabilisation)
!!$          
!!$       case (AFCSTAB_SYMMETRIC)
!!$          call gfsc_buildResidualSymm(ru, 1._DP, rrhs,&
!!$                                      rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION))
!!$       end select
!!$    end if
!!$    
!!$    ! Stop time measurement for residual/rhs evaluation
!!$    call stat_stopTimer(rtimer_assembly_resrhs)
!!$
!!$  end subroutine codire_calcRHS

  

!!$  !*****************************************************************************
!!$
!!$!<subroutine>
!!$
!!$  subroutine codire_calcJacobian(rproblemLevel, rtimestep, rsolver,&
!!$                                 ru, ru0, bfailure, rcollection)
!!$
!!$!<description>
!!$    ! This callback subroutine computes the Jacobian matrix.
!!$!</description>
!!$
!!$!<input>
!!$    ! time-stepping algorithm
!!$    type(t_timestep), intent(IN) :: rtimestep
!!$
!!$    ! initial solution vector
!!$    type(t_vectorBlock), intent(IN) :: ru0
!!$
!!$    ! Newton subiteration failed, return to defect correction
!!$    logical, intent(IN) :: bfailure
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! problem level structure
!!$    type(t_problemLevel), intent(INOUT) :: rproblemLevel
!!$
!!$    ! solver structure
!!$    type(t_solver), intent(INOUT) :: rsolver
!!$
!!$    ! solution vector
!!$    type(t_vectorBlock), intent(INOUT) :: ru
!!$
!!$    ! collection
!!$    type(t_collection), intent(INOUT) :: rcollection
!!$!</inputoutput>
!!$!</subroutine>
!!$
!!$    ! local variables
!!$    real(DP) :: hstep
!!$    logical :: bStabilize, bisExactStructure
!!$
!!$    integer :: imode = 1
!!$
!!$    ! Start time measurement for matrix evaluation
!!$    call stat_startTimer(rtimer_assembly_matrix, STAT_TIMERSHORT)
!!$
!!$    !---------------------------------------------------------------------------
!!$    ! If Newton iteration failed completely, return to
!!$    ! low-order preconditioner C=ML-theta*tstep*L(u)
!!$    !---------------------------------------------------------------------------
!!$    if (bfailure) then
!!$      ! What type of flow are we?
!!$      select case(iflowtype)
!!$      case (FLOW_TRANSIENT,&
!!$            FLOW_PSEUDOTRANSIENT)
!!$        call codire_updateSolverMatrix(rproblemLevel, rsolver, CDEQ_MATRIX_A,&
!!$                                       UPDMAT_JAC_TRANSIENT, rproblemLevel%ilev, rproblemLevel%ilev)
!!$        
!!$      case (FLOW_STEADYSTATE)
!!$        call codire_updateSolverMatrix(rproblemLevel, rsolver, CDEQ_MATRIX_A,&
!!$                                       UPDMAT_JAC_STEADY, rproblemLevel%ilev, rproblemLevel%ilev)
!!$        
!!$      case DEFAULT
!!$        call output_line('Invalid flow type!',&
!!$                         OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
!!$      end select
!!$      
!!$      ! Update the content of the solver structure
!!$      call solver_updateContent(rsolver)
!!$      
!!$      ! That's it, return
!!$      return
!!$    end if
!!$
!!$    !---------------------------------------------------------------------------
!!$    ! If Newton iteration should be performed, generate the Jacobian
!!$    ! matrix and configure the matrices for the linear solver
!!$    !---------------------------------------------------------------------------
!!$    
!!$    ! The Jacobian matrix for the low-order transport operator needs
!!$    ! to be generated only in case of nonlinear governing equations.
!!$    ! In this case, the corresponding transport operator L has to be
!!$    ! updated in each nonlinear iteration. Hence, the object
!!$    ! CDEQ_MATRIX_L can be used to store the Jacobian matrix.
!!$    
!!$
!!$    ! Compute step lenth
!!$    select case(int(rsolver%p_solverNewton%dperturbationStrategy))
!!$    case (PERTURB_NITSOL)
!!$      ! Choice h=[(1+|u|)*EPS]^{1/(1+p)} by Pernice et al. in
!!$      ! M. Pernice, H.F. Walker, NITSOL: a Newton iterative solver
!!$      ! for nonlinear systems, SIAM J. Sci. Comput. 19 (1998) 302-318.
!!$      hstep = ( (1+&
!!$          lsysbl_vectorNorm(ru, LINALG_NORMEUCLID))*SYS_EPSREAL )**(1._DP/3._DP)
!!$      
!!$    case (PERTURB_SQRTEPS)
!!$      hstep= sqrt(SYS_EPSREAL)
!!$
!!$    case DEFAULT
!!$      hstep=max(SYS_EPSREAL, rsolver%p_solverNewton%dperturbationStrategy)
!!$    end select
!!$
!!$
!!$    !---------------------------------------------------------------------------
!!$    ! Assemble diffusive contribution
!!$    !---------------------------------------------------------------------------
!!$
!!$    select case(idiffusiontype)
!!$    case (DIFFUSION_ZERO)
!!$      ! zero diffusion, clear the system matrix
!!$      call lsyssc_clearMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$
!!$      
!!$    case (DIFFUSION_ISOTROPIC)
!!$      ! Isotropic diffusion
!!$      call gfsc_buildDiffusionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
!!$                                       1._DP , .false., .true.,&
!!$                                       rmatrixDest=rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$
!!$      
!!$    case (DIFFUSION_ANISOTROPIC)
!!$      ! Anisotropic diffusion
!!$      bStabilize = (AFCSTAB_GALERKIN .ne.&
!!$                    rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%ctypeAFCstabilisation)
!!$
!!$      call gfsc_buildDiffusionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
!!$                                       1._DP, bStabilize, .true.,&
!!$                                       rmatrixDest=rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$      
!!$      
!!$    case DEFAULT
!!$      call output_line('Invalid type of diffusion!',&
!!$                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
!!$      call sys_halt()
!!$    end select
!!$
!!$
!!$    !---------------------------------------------------------------------------
!!$    ! Assemble convective contribution
!!$    !---------------------------------------------------------------------------
!!$
!!$    ! Set the velocity for current level
!!$    call codire_setVelocity(rproblemLevel)
!!$
!!$    ! Do we have to perform stabilization of AFC-type?
!!$    bStabilize = (AFCSTAB_GALERKIN .ne.&
!!$                  rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation)
!!$    
!!$    select case (abs(ivelocitytype))
!!$    case (VELOCITY_ZERO)
!!$      ! zero velocity, do nothing
!!$      
!!$    case (VELOCITY_CONSTANT,&
!!$          VELOCITY_TIMEDEP) 
!!$      ! linear velocity
!!$      select case(rproblemLevel%rtriangulation%ndim)
!!$      case (NDIM1D)
!!$        select case(imode)
!!$        case (1)
!!$          call gfsc_buildConvectionJacobian(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),ru,&
!!$                                            codire_calcPrimalConvConst1d, hstep, bStabilize,&
!!$                                            .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$        case (2)
!!$          call gfsc_buildConvectionJacobian(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),ru,&
!!$                                            codire_calcDualConvConst1d, hstep, bStabilize,&
!!$                                            .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$        end select
!!$
!!$      case (NDIM2D)
!!$        select case(imode)
!!$        case (1)
!!$          call gfsc_buildConvectionJacobian(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),ru,&
!!$                                            codire_calcPrimalConvConst2d, hstep, bStabilize,&
!!$                                            .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$        case (2)
!!$          call gfsc_buildConvectionJacobian(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),ru,&
!!$                                            codire_calcDualConvConst2d, hstep, bStabilize,&
!!$                                            .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$        end select
!!$
!!$      case (NDIM3D)
!!$        select case(imode)
!!$        case (1)
!!$          call gfsc_buildConvectionJacobian(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CZ),ru,&
!!$                                            codire_calcPrimalConvConst3d, hstep, bStabilize,&
!!$                                            .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$        case (2)
!!$          call gfsc_buildConvectionJacobian(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CZ),ru,&
!!$                                            codire_calcDualConvConst3d, hstep, bStabilize,&
!!$                                            .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$        end select
!!$
!!$      case DEFAULT
!!$        call output_line('Invalid spatial dimension!',&
!!$                         OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
!!$        call sys_halt()
!!$      end select
!!$
!!$      
!!$    case (VELOCITY_BURGERS_SPACETIME)
!!$      ! nonlinear Burgers' equation in space-time
!!$      call gfsc_buildConvectionJacobian(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY), ru,&
!!$                                        codire_calcConvectionBurgersSpT2d, hstep, bStabilize,&
!!$                                        .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$      
!!$    case (VELOCITY_BUCKLEV_SPACETIME)
!!$      ! nonlinear Buckley-Leverett equation in space-time
!!$      call gfsc_buildConvectionJacobian(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY), ru,&
!!$                                        codire_calcConvectionBuckLevSpT2d, hstep, bStabilize,&
!!$                                        .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$      
!!$    case (VELOCITY_BURGERS1D)
!!$      ! nonlinear Burgers' equation in 1D
!!$      call gfsc_buildConvectionJacobian(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX), ru,&
!!$                                        codire_calcConvectionBurgers1d, hstep, bStabilize,&
!!$                                        .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$
!!$    case (VELOCITY_BURGERS2D)
!!$      ! nonlinear Burgers' equation in 2D
!!$      call gfsc_buildConvectionJacobian(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY), ru,&
!!$                                        codire_calcConvectionBurgers2d, hstep, bStabilize,&
!!$                                        .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$
!!$    case (VELOCITY_BUCKLEV1D)
!!$      ! nonlinear Buckley-Leverett equation in 1D
!!$      call gfsc_buildConvectionJacobian(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX), ru,&
!!$                                        codire_calcConvectionBuckLev1d, hstep, bStabilize,&
!!$                                        .false., rproblemLevel%Rmatrix(CDEQ_MATRIX_L))
!!$      
!!$    case DEFAULT
!!$      call output_line('Unsupported velocity type!',&
!!$                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
!!$      call sys_halt()
!!$    end select
!!$    
!!$    
!!$    ! Check if Jacobian matrix has the same structure as the system matrix
!!$    bisExactStructure = .true.
!!$
!!$    ! Check for convective stabilization
!!$    if (rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%iextendedJacobian .ne. 0) then
!!$      select case(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation)
!!$      case (AFCSTAB_FEMGP,&
!!$            AFCSTAB_FEMTVD)
!!$        bisExactStructure = .false.
!!$      end select
!!$    end if
!!$        
!!$    ! Check for diffusive stabilization
!!$    if (rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%iextendedJacobian .ne. 0) then
!!$      select case(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%ctypeAFCstabilisation)
!!$      case (AFCSTAB_SYMMETRIC)
!!$        bisExactStructure = .false.
!!$      end select
!!$    end if
!!$
!!$
!!$    ! Check if the global system operator needs to be generated
!!$    if (rproblemLevel%Rmatrix(CDEQ_MATRIX_J)%cmatrixFormat .eq. LSYSSC_MATRIXUNDEFINED) then
!!$      
!!$      if (bisExactStructure) then
!!$        ! Adopt the standard sparsity pattern
!!$        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
!!$                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
!!$                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
!!$      else
!!$        ! Extend the standard sparsity pattern
!!$        call afcstab_generateExtSparsity(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
!!$                                         rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$      end if
!!$    end if
!!$
!!$    
!!$    !---------------------------------------------------------------------------
!!$    ! Assemble the global system operator for the high-/low-order contribution
!!$    !---------------------------------------------------------------------------
!!$    
!!$    ! What type of flow are we?
!!$    select case(iflowtype)
!!$    case (FLOW_TRANSIENT,&
!!$          FLOW_PSEUDOTRANSIENT)
!!$      ! Compute the global Jacobian for transient flow
!!$      !
!!$      !     J = ML-theta*dt*L
!!$      !
!!$      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(CDEQ_MATRIX_L),&
!!$                                   -rtimestep%theta*rtimestep%dStep,&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_ML), 1._DP,&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
!!$                                   .false., .false., .true., bisExactStructure)
!!$
!!$    case (FLOW_STEADYSTATE)
!!$      ! Compute the global Jacobian for steady-state flow
!!$      !
!!$      !     J = -L
!!$      !
!!$      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(CDEQ_MATRIX_L), -1._DP,&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J), 0._DP,&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
!!$                                   .false., .false., .true., bisExactStructure)
!!$
!!$    case DEFAULT
!!$      call output_line('Invalid flow type!',&
!!$                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
!!$    end select
!!$    
!!$
!!$    ! What kind of diffusion are we?
!!$    select case(idiffusiontype)
!!$    case (DIFFUSION_ZERO, DIFFUSION_ISOTROPIC)
!!$      ! zero diffusion or isotropic diffusion, do nothing
!!$      
!!$    case (DIFFUSION_ANISOTROPIC)
!!$      ! Anisotropic diffusion
!!$      select case(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION)%ctypeAFCstabilisation)
!!$
!!$      case (AFCSTAB_SYMMETRIC)
!!$        call gfsc_buildJacobianSymm(ru, 1._DP, hstep, .false.,&
!!$                                    rproblemLevel%Rafcstab(CDEQ_AFCSTAB_DIFFUSION),&
!!$                                    rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$      end select
!!$
!!$      
!!$    case DEFAULT
!!$      call output_line('Invalid type of diffusion!',&
!!$                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
!!$      call sys_halt()
!!$    end select
!!$
!!$
!!$    ! What kind of velocity are we?
!!$    select case(abs(ivelocitytype))
!!$    case (VELOCITY_ZERO)
!!$      ! zero velocity, do nothing
!!$      
!!$    case(VELOCITY_CONSTANT,&
!!$         VELOCITY_TIMEDEP) 
!!$      ! linear velocity
!!$      select case(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation)
!!$
!!$      case (AFCSTAB_FEMFCT)
!!$        call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC), ru,&
!!$                                   rtimestep%theta, rtimestep%dStep, hstep, .false.,&
!!$                                   rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$        
!!$      case (AFCSTAB_FEMTVD,&
!!$            AFCSTAB_FEMGP)
!!$        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(CDEQ_MATRIX_MC),ru, ru0,&
!!$                                   rtimestep%theta, rtimestep%dStep, hstep, .false.,&
!!$                                   rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$      end select
!!$      
!!$
!!$    case(VELOCITY_BURGERS_SPACETIME)
!!$      ! nonlinear Burgers' equation in space-time
!!$      select case(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation)
!!$
!!$      case (AFCSTAB_FEMFCT)
!!$        call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_MC), ru,&
!!$                                   codire_calcConvectionBurgersSpT2d,&
!!$                                   rtimestep%theta, rtimestep%dStep, hstep, .false.,&
!!$                                   rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$
!!$      case (AFCSTAB_FEMTVD,&
!!$            AFCSTAB_FEMGP)
!!$        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_MC), ru, ru0,&
!!$                                   codire_calcConvectionBurgersSpT2d,&
!!$                                   rtimestep%theta, rtimestep%dStep, hstep, .false.,&
!!$                                   rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$      end select
!!$        
!!$      
!!$    case(VELOCITY_BUCKLEV_SPACETIME)
!!$      ! nonlinear Buckley-Leverett equation in space-time
!!$      select case(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation)
!!$
!!$      case (AFCSTAB_FEMFCT)
!!$        call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_MC), ru,&
!!$                                   codire_calcConvectionBuckLevSpT2d,&
!!$                                   rtimestep%theta, rtimestep%dStep, hstep, .false.,&
!!$                                   rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$
!!$      case (AFCSTAB_FEMTVD,&
!!$            AFCSTAB_FEMGP)
!!$        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_MC), ru, ru0,&
!!$                                   codire_calcConvectionBuckLevSpT2d,&
!!$                                   rtimestep%theta, rtimestep%dStep, hstep, .false.,&
!!$                                   rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$      end select
!!$
!!$      
!!$    case(VELOCITY_BURGERS1D)
!!$      ! nonlinear Burgers' equation in 1D
!!$      select case(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation)
!!$
!!$      case (AFCSTAB_FEMFCT)
!!$        call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_MC), ru,&
!!$                                   codire_calcConvectionBurgers1d,&
!!$                                   rtimestep%theta, rtimestep%dStep, hstep, .false.,&
!!$                                   rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$
!!$      case (AFCSTAB_FEMTVD,&
!!$            AFCSTAB_FEMGP)
!!$        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_MC), ru, ru0,&
!!$                                   codire_calcConvectionBurgers1d,&
!!$                                   rtimestep%theta, rtimestep%dStep, hstep, .false.,&
!!$                                   rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$      end select
!!$
!!$
!!$    case(VELOCITY_BURGERS2D)
!!$      ! nonlinear Burgers' equation in 2D
!!$      select case(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation)
!!$
!!$      case (AFCSTAB_FEMFCT)
!!$        call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_MC), ru,&
!!$                                   codire_calcConvectionBurgers2d,&
!!$                                   rtimestep%theta, rtimestep%dStep, hstep, .false.,&
!!$                                   rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$
!!$      case (AFCSTAB_FEMTVD,&
!!$            AFCSTAB_FEMGP)
!!$        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_MC), ru, ru0,&
!!$                                   codire_calcConvectionBurgers2d,&
!!$                                   rtimestep%theta, rtimestep%dStep, hstep, .false.,&
!!$                                   rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$      end select
!!$
!!$
!!$    case(VELOCITY_BUCKLEV1D)
!!$      ! nonlinear Buckley-Leverett equation in 1D
!!$      select case(rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION)%ctypeAFCstabilisation)
!!$
!!$      case (AFCSTAB_FEMFCT)
!!$        call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_MC), ru,&
!!$                                   codire_calcConvectionBuckLev1d,&
!!$                                   rtimestep%theta, rtimestep%dStep, hstep, .false.,&
!!$                                   rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$
!!$      case (AFCSTAB_FEMTVD,&
!!$            AFCSTAB_FEMGP)
!!$        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_MC), ru, ru0,&
!!$                                   codire_calcConvectionBuckLev1d,&
!!$                                   rtimestep%theta, rtimestep%dStep, hstep, .false.,&
!!$                                   rproblemLevel%Rafcstab(CDEQ_AFCSTAB_CONVECTION),&
!!$                                   rproblemLevel%Rmatrix(CDEQ_MATRIX_J))
!!$      end select
!!$
!!$
!!$    case DEFAULT
!!$      call output_line('Unsupported velocity type!',&
!!$                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
!!$      call sys_halt()
!!$    end select
!!$    
!!$    !---------------------------------------------------------------------------
!!$    ! Impose boundary conditions
!!$    !---------------------------------------------------------------------------
!!$    
!!$    call bdrf_filterMatrix(rsolver%rboundaryCondition, rproblemLevel%rtriangulation,&
!!$                           rproblemLevel%Rmatrix(CDEQ_MATRIX_J), 1.0_DP)
!!$
!!$    ! Ok, we updated the Jacobian matrix successfully. Now we still have to
!!$    ! link it to the solver hierarchy. This is done recursively.
!!$
!!$    ! What type of flow are we?
!!$    select case(iflowtype)
!!$    case (FLOW_TRANSIENT,&
!!$          FLOW_PSEUDOTRANSIENT)
!!$      call codire_updateSolverMatrix(rproblemLevel, rsolver, CDEQ_MATRIX_J, &
!!$                                  UPDMAT_JAC_TRANSIENT, rproblemLevel%ilev, rproblemLevel%ilev)
!!$
!!$    case (FLOW_STEADYSTATE)
!!$      call codire_updateSolverMatrix(rproblemLevel, rsolver, CDEQ_MATRIX_J, &
!!$                                  UPDMAT_JAC_STEADY, rproblemLevel%ilev, rproblemLevel%ilev)
!!$
!!$    case DEFAULT
!!$      call output_line('Invalid flow type!',&
!!$                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
!!$      call sys_halt()
!!$    end select
!!$
!!$    ! Finally, we have to update the content of the solver hierarchy
!!$    call solver_updateContent(rsolver)
!!$    
!!$    ! Stop time measurement for matrix evaluation
!!$    call stat_stopTimer(rtimer_assembly_matrix)
!!$
!!$  end subroutine codire_calcJacobian

!!$  !*****************************************************************************
!!$
!!$!<subroutine>
!!$
!!$  subroutine codire_applyJacobian(rproblemLevel, rx, ry, cx, cy, rcollection)
!!$
!!$!<description>
!!$    ! This subroutine applies the (scaled) Jacobian matrix to 
!!$    ! the vector rx and adds the result to the vector ry.
!!$!</description>
!!$
!!$!<input>
!!$    ! multigrid level structure
!!$    type(t_problemLevel), intent(IN) :: rproblemLevel
!!$
!!$    ! vector to which the Jacobian should be applied to
!!$    type(t_vectorBlock), intent(IN) :: rx
!!$
!!$    ! factor by which rx should be scaled
!!$    real(DP), intent(IN) :: cx
!!$
!!$    ! factor by which ry should be scaled
!!$    real(DP), intent(IN) :: cy
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! vector to which the result should be added
!!$    type(t_vectorBlock), intent(INOUT) :: ry
!!$
!!$    ! collection
!!$    type(t_collection), intent(INOUT) :: rcollection
!!$!</inputoutput>
!!$!</subroutine>
!!$
!!$    ! We know where the Jacobian matrix is stored and can apply it
!!$    ! by means of matrix vector multiplication
!!$    call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(CDEQ_MATRIX_J),&
!!$                             rx%RvectorBlock(1),&
!!$                             ry%RvectorBlock(1), cx, cy)
!!$
!!$  end subroutine codire_applyJacobian

!!$  !*****************************************************************************
!!$
!!$!<subroutine>
!!$
!!$  subroutine codire_calcGalerkinResidual(rproblemLevel, rsolution,&
!!$                                         rresidual, rcollection, rf)
!!$
!!$!<description>
!!$    ! This subroutine calculates the Galerkin residual for a
!!$    ! (converged) steady-state solution to the primal problem
!!$!</description>
!!$
!!$!<input>
!!$    ! solution vector
!!$    type(t_vectorBlock), intent(IN) :: rsolution
!!$
!!$    ! OPTIONAL: right-hand side vector
!!$    type(t_vectorBlock), intent(IN), optional :: rf
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! problem level structure
!!$    type(t_problemLevel), intent(INOUT) :: rproblemLevel
!!$
!!$    ! residual vector
!!$    type(t_vectorBlock), intent(INOUT) :: rresidual
!!$
!!$    ! collection
!!$    type(t_collection), intent(INOUT) :: rcollection
!!$!</inputoutput>
!!$!</subroutine>
!!$
!!$    ! local variables
!!$    type(t_matrixScalar) :: rmatrix
!!$    logical :: bcompatible,bisDivFree
!!$
!!$    integer :: imode = 1
!!$
!!$    ! Check if vectors are compatible
!!$    call lsysbl_isVectorCompatible(rsolution, rresidual, bcompatible)
!!$    
!!$    ! Re-create residual vector if not compatible
!!$    if (.not.bcompatible) then
!!$      call lsysbl_releaseVector(rresidual)
!!$      call lsysbl_duplicateVector(rsolution, rresidual,&
!!$                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
!!$    end if
!!$
!!$    ! Create empty matrix for the global system operator
!!$    call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
!!$                                rmatrix, LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
!!$
!!$    !------------------------------------------------------------
!!$    ! Assemble diffusive operator
!!$    !------------------------------------------------------------
!!$    
!!$    select case(idiffusiontype)
!!$    case (DIFFUSION_ZERO)
!!$      ! zero diffusion, clear the system matrix
!!$      call lsyssc_clearMatrix(rmatrix)
!!$
!!$      
!!$    case (DIFFUSION_ISOTROPIC)
!!$      ! Isotropic diffusion
!!$      call gfsc_buildDiffusionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
!!$                                       1._DP , .false., .true.,&
!!$                                       rmatrixDest=rmatrix)
!!$      
!!$      
!!$    case (DIFFUSION_ANISOTROPIC)
!!$      ! Assemble high-order transport operator on coarser levels
!!$      call gfsc_buildDiffusionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_S),&
!!$                                       1._DP, .false., .true.,&
!!$                                       rmatrixDest=rmatrix)
!!$      
!!$    case DEFAULT
!!$      call output_line('Invalid type of diffusion!',&
!!$                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcGalerkinResidual')
!!$      call sys_halt()
!!$    end select
!!$
!!$    !---------------------------------------------------------------------------
!!$    ! Assemble convective operator
!!$    !---------------------------------------------------------------------------
!!$    
!!$    ! Set velocity vector for current level
!!$    call codire_setVelocity(rproblemLevel)
!!$
!!$    ! Check if velocity is assumed to be discretely divergence free
!!$    bisDivFree = (ivelocitytype .gt. 0)
!!$
!!$    ! Assemble high-order transport operator on coarser levels
!!$    select case(abs(ivelocitytype))
!!$    case (VELOCITY_ZERO)
!!$      ! zero velocity, do nothing
!!$      
!!$    case (VELOCITY_CONSTANT,&
!!$          VELOCITY_TIMEDEP)
!!$      ! linear velocity
!!$      select case(rproblemLevel%rtriangulation%ndim)
!!$      case (NDIM1D)
!!$        select case(imode)
!!$        case (1)
!!$          call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                            rsolution, codire_calcPrimalConvConst1d,&
!!$                                            bisDivFree, .false., .false., rmatrix)
!!$        case (2)
!!$          call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                            rsolution, codire_calcDualConvConst1d,&
!!$                                            bisDivFree, .false., .false., rmatrix)
!!$        end select
!!$          
!!$      case (NDIM2D)
!!$        select case(imode)
!!$          case (1)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                              rsolution, codire_calcPrimalConvConst2d,&
!!$                                              bisDivFree, .false., .false., rmatrix)
!!$          case (2)
!!$            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                              rsolution, codire_calcDualConvConst2d,&
!!$                                              bisDivFree, .false., .false., rmatrix)
!!$          end select
!!$          
!!$      case (NDIM3D)
!!$        select case(imode)
!!$        case (1)
!!$          call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CZ),&
!!$                                            rsolution, codire_calcPrimalConvConst3d,&
!!$                                            bisDivFree, .false., .false., rmatrix)
!!$        case (2)
!!$          call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CZ),&
!!$                                            rsolution, codire_calcDualConvConst3d,&
!!$                                            bisDivFree, .false., .false., rmatrix)
!!$        end select
!!$
!!$      case DEFAULT
!!$        call output_line('Invalid spatial dimension!',&
!!$                         OU_CLASS_ERROR,OU_MODE_STD,'codire_calcGalerkinResidual')
!!$        call sys_halt()
!!$      end select
!!$      
!!$      
!!$    case (VELOCITY_BURGERS_SPACETIME)
!!$      ! nonlinear Burgers' equation in space-time
!!$      call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                        rsolution, codire_calcConvectionBurgersSpT2d,&
!!$                                        bisDivFree, .false., .false., rmatrix)
!!$      
!!$    case (VELOCITY_BUCKLEV_SPACETIME)
!!$      ! nonlinear Buckley-Leverett equation in space-time
!!$      call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                        rsolution, codire_calcConvectionBuckLevSpT2d,&
!!$                                        bisDivFree, .false., .false., rmatrix)
!!$      
!!$    case (VELOCITY_BURGERS1D)
!!$      ! nonlinear Burgers' equation in 1D
!!$      call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                        rsolution, codire_calcConvectionBurgers1d,&
!!$                                        bisDivFree, .false., .false., rmatrix)
!!$
!!$    case (VELOCITY_BURGERS2D)
!!$      ! nonlinear Burgers' equation in 2D
!!$      call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
!!$                                        rsolution, codire_calcConvectionBurgers2d,&
!!$                                        bisDivFree, .false., .false., rmatrix)
!!$      
!!$    case (VELOCITY_BUCKLEV1D)
!!$      ! nonlinear Buckley-Leverett equation in 1D
!!$      call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
!!$                                        rsolution, codire_calcConvectionBuckLev1d,&
!!$                                        bisDivFree, .false., .false., rmatrix)
!!$      
!!$    case DEFAULT
!!$      call output_line('Invalid velocity profile!',&
!!$                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcGalerkinResidual')
!!$      call sys_halt()
!!$    end select
!!$    
!!$    ! Compute the Galerkin residual
!!$    if (present(rf)) then
!!$      call lsysbl_copyVector(rf, rresidual)
!!$      call lsyssc_scalarMatVec(rmatrix, rsolution%RvectorBlock(1),&
!!$                               rresidual%RvectorBlock(1), -1._DP, 1._DP)
!!$    else
!!$      call lsyssc_scalarMatVec(rmatrix, rsolution%RvectorBlock(1),&
!!$                               rresidual%RvectorBlock(1), -1._DP, 0._DP)
!!$    end if
!!$
!!$    ! Release memory
!!$    call lsyssc_releaseMatrix(rmatrix)
!!$
!!$  end subroutine codire_calcGalerkinResidual
  
  !*****************************************************************************
  !*****************************************************************************
  !*****************************************************************************

!!$!<subroutine>
!!$
!!$  subroutine codire_setVelocity(rproblemLevel)
!!$
!!$!<description>
!!$    ! This subroutine sets the internal velocity vector 
!!$    ! to the level of the given multigrid structure.
!!$!</description>
!!$
!!$!<input>
!!$    ! multigrid structure
!!$    type(t_problemLevel), intent(IN) :: rproblemLevel
!!$!</input>
!!$!</subroutine>
!!$
!!$    ! What spatial dimension are we?
!!$    select case(rproblemLevel%rtriangulation%ndim)
!!$    case (NDIM1D)
!!$      call lsyssc_getbase_double(rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%RvectorBlock(1), p_DvelocityX)
!!$
!!$    case (NDIM2D)
!!$      call lsyssc_getbase_double(rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%RvectorBlock(1), p_DvelocityX)
!!$      call lsyssc_getbase_double(rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%RvectorBlock(2), p_DvelocityY)
!!$
!!$    case (NDIM3D)
!!$      call lsyssc_getbase_double(rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%RvectorBlock(1), p_DvelocityX)
!!$      call lsyssc_getbase_double(rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%RvectorBlock(2), p_DvelocityY)
!!$      call lsyssc_getbase_double(rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%RvectorBlock(3), p_DvelocityZ)
!!$
!!$    case DEFAULT
!!$      call output_line('Invalid spatial dimension!',&
!!$                       OU_CLASS_ERROR,OU_MODE_STD,'codire_setVelocity')
!!$      call sys_halt()
!!$    end select
!!$    
!!$  end subroutine codire_setVelocity

!!$  !*****************************************************************************
!!$
!!$!<function>
!!$
!!$  function codire_updateVelocity(rproblem, ttime, nlminOpt, nlmaxOpt) result(bmodified)
!!$
!!$!<description>
!!$    ! This subroutine updates the internal velocity vector for all levels.
!!$    ! If the optional parameters NLMIN and NLMAX are given, then only those
!!$    ! multigrid levels are updated; otherwise, all levels are updated.
!!$!</description>
!!$    
!!$!<input>
!!$    ! problem data structure
!!$    type(t_problem), intent(IN) :: rproblem
!!$
!!$    ! simulation time
!!$    real(DP), intent(IN) :: ttime
!!$
!!$    ! OPTIONAL: minimum multigrid level
!!$    integer, intent(IN), optional :: nlminOpt
!!$
!!$    ! OPTIONAL: maximum multigrid level
!!$    integer, intent(IN), optional :: nlmaxOpt    
!!$!</input>
!!$
!!$!<result>
!!$    ! is TRUE if velocity vector has been modified
!!$    logical :: bmodified
!!$!</result>
!!$!</function>
!!$    
!!$    ! local variables
!!$    type(t_problemLevel), pointer :: rproblemLevel
!!$    real(DP), dimension(:,:), pointer :: p_Dcoords
!!$    real(DP), dimension(1) :: Dvalue
!!$    integer(I32), dimension(2) :: Isize
!!$    integer :: NEQ,nlmin,nlmax
!!$
!!$    ! Set minimal/maximal levels: If given by the user, adopt these values.
!!$    ! Otherwise, initialize velocity vector on all levels of the multigrid
!!$    ! structure.
!!$    if (present(nlminOpt)) then
!!$      nlmin = nlminOpt
!!$    else
!!$      nlmin = rproblem%p_rproblemLevelMin%ilev
!!$    end if
!!$    
!!$    if (present(nlmaxOpt)) then
!!$      nlmax = nlmaxOpt
!!$    else
!!$      nlmax = rproblem%p_rproblemLevelMax%ilev
!!$    end if
!!$
!!$    ! Set time value
!!$    Dvalue = ttime
!!$    
!!$    ! Check if updated of velocity is forced. Note that the variable 
!!$    ! bvelocityUpdate is global and can be modified from external routines.
!!$    if (bvelocityUpdate) then
!!$
!!$      ! Ok, unmark velocity update and set modified to TRUE.
!!$      bvelocityUpdate = .false.
!!$      bmodified       = .true.
!!$
!!$      ! Process all multigrid levels
!!$      rproblemLevel => rproblem%p_rproblemLevelMax
!!$      do while(associated(rproblemLevel))
!!$        
!!$        ! Do we have to initialize the velocity vector for this level?
!!$        if (rproblemLevel%ilev < nlmin .or. rproblemLevel%ilev > nlmax) then
!!$          rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
!!$          cycle
!!$        end if
!!$        
!!$        ! Get total number of vertices [NVT,NVT]
!!$        Isize = rproblemLevel%rtriangulation%NVT
!!$        
!!$        ! Check if velocity vector exists 
!!$        if (associated(rproblemLevel%rvectorBlock)) then
!!$          
!!$          ! Check if velocity vector needs to be resized
!!$          NEQ = rproblemLevel%rvectorBlock(CDEQ_VELOCITY)%NEQ
!!$          if (NEQ .eq. 0) then
!!$            ! Create two-dimensional velocity vector
!!$            call lsysbl_createVectorBlock(rproblemLevel%rvectorBlock(CDEQ_VELOCITY),&
!!$                                          Isize, .true.)
!!$          elseif(NEQ .ne. sum(Isize)) then
!!$            ! Resize two-dimensional velocity vector
!!$            call lsysbl_resizeVectorBlock(rproblemLevel%rvectorBlock(CDEQ_VELOCITY),&
!!$                                          Isize, .true.)
!!$          end if
!!$          
!!$        else
!!$          
!!$          ! Allocate block vector for velocity
!!$          allocate(rproblemLevel%rvectorBlock(CDEQ_VELOCITY))
!!$          
!!$          ! Create two-dimensional velocity vector
!!$          call lsysbl_createVectorBlock(rproblemLevel%rvectorBlock(CDEQ_VELOCITY),&
!!$                                        Isize,.true.)
!!$          
!!$        end if
!!$        
!!$        ! Are we on the finest level? 
!!$        if (rproblemLevel%ilev .eq. nlmax) then
!!$          
!!$          ! Evaluate function parser
!!$          call codire_setVelocity(rproblemLevel)
!!$          call storage_getbase_double2d(rproblemLevel%rtriangulation%h_DvertexCoords,  p_Dcoords)
!!$          p_Dcoords => p_Dcoords(:,1:rproblemLevel%rtriangulation%NVT)
!!$
!!$          ! What spatial dimension are we?
!!$          select case (rproblemLevel%rtriangulation%ndim)
!!$          case (NDIM1D)
!!$            call fparser_evalFunction(rvelocityParser, 1, 2, p_Dcoords, p_DvelocityX, Dvalue)
!!$
!!$          CASE (NDIM2D)
!!$            call fparser_evalFunction(rvelocityParser, 1, 2, p_Dcoords, p_DvelocityX, Dvalue)
!!$            call fparser_evalFunction(rvelocityParser, 2, 2, p_Dcoords, p_DvelocityY, Dvalue)
!!$
!!$          CASE (NDIM3D)
!!$            call fparser_evalFunction(rvelocityParser, 1, 2, p_Dcoords, p_DvelocityX, Dvalue)
!!$            call fparser_evalFunction(rvelocityParser, 2, 2, p_Dcoords, p_DvelocityY, Dvalue)
!!$            call fparser_evalFunction(rvelocityParser, 3, 2, p_Dcoords, p_DvelocityZ, Dvalue)
!!$
!!$          case DEFAULT
!!$            call output_line('Invalid spatial dimension!',&
!!$                             OU_CLASS_ERROR,OU_MODE_STD,'codire_updateVelocity')
!!$            call sys_halt()
!!$          end select
!!$          
!!$        else
!!$          
!!$          ! Restrict velocity vector from finer grid
!!$          call solver_restrictionBlock(&
!!$              rproblemLevel%p_rproblemLevelFine%rtriangulation,&
!!$              rproblemLevel%rtriangulation,&
!!$              rproblemLevel%p_rproblemLevelFine%rvectorBlock(CDEQ_VELOCITY),&
!!$              rproblemLevel%rvectorBlock(CDEQ_VELOCITY))
!!$          
!!$        end if
!!$        
!!$        ! Switch to next coarser level
!!$        rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
!!$      end do
!!$
!!$    else
!!$
!!$      ! What kind of velocity are we?
!!$      select case(abs(ivelocitytype))
!!$      case (VELOCITY_ZERO,&
!!$            VELOCITY_CONSTANT)
!!$        bmodified = .false.
!!$        
!!$      case (VELOCITY_TIMEDEP)
!!$        bmodified = .true.
!!$        
!!$        ! Process all multigrid levels
!!$        rproblemLevel => rproblem%p_rproblemLevelMax
!!$        do while(associated(rproblemLevel))
!!$          
!!$          ! Do we have to initialize the velocity vector for this level?
!!$          if (rproblemLevel%ilev < nlmin .or. rproblemLevel%ilev > nlmax) then
!!$            rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
!!$            cycle
!!$          end if
!!$          
!!$          ! Are we on the finest level? 
!!$          if (rproblemLevel%ilev .eq. nlmax) then
!!$            
!!$            ! Call function parser
!!$            call codire_setVelocity(rproblemLevel)
!!$            call storage_getbase_double2d(rproblemLevel%rtriangulation%h_DvertexCoords, p_Dcoords)
!!$            p_Dcoords => p_Dcoords(:,1:rproblemLevel%rtriangulation%NVT)
!!$
!!$            ! What spatial dimension are we?
!!$            select case (rproblemLevel%rtriangulation%ndim)
!!$            case (NDIM1D)
!!$              call fparser_evalFunction(rvelocityParser, 1, 2, p_Dcoords, p_DvelocityX, Dvalue)
!!$              
!!$            CASE (NDIM2D)
!!$              call fparser_evalFunction(rvelocityParser, 1, 2, p_Dcoords, p_DvelocityX, Dvalue)
!!$              call fparser_evalFunction(rvelocityParser, 2, 2, p_Dcoords, p_DvelocityY, Dvalue)
!!$              
!!$            CASE (NDIM3D)
!!$              call fparser_evalFunction(rvelocityParser, 1, 2, p_Dcoords, p_DvelocityX, Dvalue)
!!$              call fparser_evalFunction(rvelocityParser, 2, 2, p_Dcoords, p_DvelocityY, Dvalue)
!!$              call fparser_evalFunction(rvelocityParser, 3, 2, p_Dcoords, p_DvelocityZ, Dvalue)
!!$              
!!$            case DEFAULT
!!$              call output_line('Invalid spatial dimension!',&
!!$                               OU_CLASS_ERROR,OU_MODE_STD,'codire_updateVelocity')
!!$              call sys_halt()
!!$            end select
!!$            
!!$          else
!!$            
!!$            ! Restrict velocity vector from finer grid
!!$            call solver_restrictionBlock(&
!!$                rproblemLevel%p_rproblemLevelFine%rtriangulation,&
!!$                rproblemLevel%rtriangulation,&
!!$                rproblemLevel%p_rproblemLevelFine%rvectorBlock(CDEQ_VELOCITY),&
!!$                rproblemLevel%rvectorBlock(CDEQ_VELOCITY))
!!$            
!!$          end if
!!$          
!!$          ! Switch to next coarser level
!!$          rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
!!$        end do
!!$        
!!$      case (VELOCITY_BURGERS_SPACETIME:VELOCITY_BURGERS2D)
!!$        bmodified = .true.
!!$        
!!$        ! Process all multigrid levels
!!$        rproblemLevel => rproblem%p_rproblemLevelMax
!!$        do while(associated(rproblemLevel))
!!$          
!!$          ! Do we have to initialize the velocity vector for this level?
!!$          if (rproblemLevel%ilev < nlmin .or. rproblemLevel%ilev > nlmax) then
!!$            rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
!!$            cycle
!!$          end if
!!$          
!!$          ! Are we on the finest level? 
!!$          if (rproblemLevel%ilev .eq. nlmax) then
!!$            
!!$            ! Call function parser
!!$            call codire_setVelocity(rproblemLevel)
!!$            call storage_getbase_double2d(rproblemLevel%rtriangulation%h_DvertexCoords, p_Dcoords)
!!$            p_Dcoords => p_Dcoords(:,1:rproblemLevel%rtriangulation%NVT)
!!$
!!$
!!$            call fparser_evalFunction(rvelocityParser, 2, 2, p_Dcoords, p_DvelocityY, Dvalue)
!!$            
!!$          else
!!$            
!!$            ! Restrict velocity vector from finer grid
!!$            call solver_restrictionBlock(&
!!$                rproblemLevel%p_rproblemLevelFine%rtriangulation,&
!!$                rproblemLevel%rtriangulation,&
!!$                rproblemLevel%p_rproblemLevelFine%rvectorBlock(CDEQ_VELOCITY),&
!!$                rproblemLevel%rvectorBlock(CDEQ_VELOCITY))
!!$            
!!$          end if
!!$          
!!$          ! Switch to next coarser level
!!$          rproblemLevel => rproblemLevel%p_rproblemLevelCoarse
!!$        end do
!!$        
!!$      case DEFAULT
!!$        call output_line('Invalid type of velocity!',&
!!$                         OU_CLASS_ERROR,OU_MODE_STD,'codire_updateVelocity')
!!$        call sys_halt()
!!$      end select
!!$    end if
!!$  end function codire_updateVelocity
!!$
!!$  !*****************************************************************************
!!$
!!$!<subroutine>
!!$
!!$  pure subroutine codire_calcPrimalConvConst1d(u_i, u_j, C_ij, C_ji, i, j,&
!!$                                               k_ij, k_ji, istatus)
!!$
!!$!<description>
!!$    ! This subroutine computes the convective matrix coefficients
!!$    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
!!$    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 1D
!!$!</description>
!!$    
!!$!<input>
!!$    ! solution vector
!!$    real(DP), intent(IN) :: u_i,u_j
!!$
!!$    ! coefficients from spatial discretization
!!$    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!!$
!!$    ! nodal indices
!!$    integer, intent(IN) :: i,j
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! OPTIONAL: status of the callback function
!!$    integer, intent(INOUT), optional :: istatus
!!$!</inputoutput>
!!$
!!$!<output>
!!$    ! convective coefficients
!!$    real(DP), intent(OUT) :: k_ij,k_ji
!!$!</output>
!!$!</subroutine>
!!$
!!$    k_ij = -p_DvelocityX(j)*C_ij(1)
!!$    k_ji = -p_DvelocityX(i)*C_ji(1)
!!$
!!$  end subroutine codire_calcPrimalConvConst1d
!!$
!!$  !*****************************************************************************
!!$
!!$!<subroutine>
!!$
!!$  pure subroutine codire_calcDualConvConst1d(u_i, u_j, C_ij, C_ji, i, j,&
!!$                                          k_ij, k_ji, istatus)
!!$
!!$!<description>
!!$    ! This subroutine computes the convective matrix coefficients
!!$    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
!!$    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 1D
!!$!</description>
!!$    
!!$!<input>
!!$    ! solution vector
!!$    real(DP), intent(IN) :: u_i,u_j
!!$
!!$    ! coefficients from spatial discretization
!!$    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!!$
!!$    ! nodal indices
!!$    integer, intent(IN) :: i,j
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! OPTIONAL: status of the callback function
!!$    integer, intent(INOUT), optional :: istatus
!!$!</inputoutput>
!!$
!!$!<output>
!!$    ! convective coefficients
!!$    real(DP), intent(OUT) :: k_ij,k_ji
!!$!</output>
!!$!</subroutine>
!!$
!!$    k_ij = -p_DvelocityX(j)*C_ji(1)
!!$    k_ji = -p_DvelocityX(i)*C_ij(1)
!!$
!!$  end subroutine codire_calcDualConvConst1d
!!$
!!$  !*****************************************************************************
!!$
!!$!<subroutine>
!!$
!!$  pure subroutine codire_calcPrimalConvConst2d(u_i, u_j, C_ij, C_ji, i, j,&
!!$                                            k_ij, k_ji, istatus)
!!$
!!$!<description>
!!$    ! This subroutine computes the convective matrix coefficients
!!$    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
!!$    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D
!!$!</description>
!!$    
!!$!<input>
!!$    ! solution vector
!!$    real(DP), intent(IN) :: u_i,u_j
!!$
!!$    ! coefficients from spatial discretization
!!$    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!!$
!!$    ! nodal indices
!!$    integer, intent(IN) :: i,j
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! OPTIONAL: status of the callback function
!!$    integer, intent(INOUT), optional :: istatus
!!$!</inputoutput>
!!$
!!$!<output>
!!$    ! convective coefficients
!!$    real(DP), intent(OUT) :: k_ij,k_ji
!!$!</output>
!!$!</subroutine>
!!$
!!$    k_ij = -p_DvelocityX(j)*C_ij(1)-p_DvelocityY(j)*C_ij(2)
!!$    k_ji = -p_DvelocityX(i)*C_ji(1)-p_DvelocityY(i)*C_ji(2)
!!$
!!$  end subroutine codire_calcPrimalConvConst2d
!!$
!!$  !*****************************************************************************
!!$
!!$!<subroutine>
!!$
!!$  pure subroutine codire_calcDualConvConst2d(u_i, u_j, C_ij, C_ji, i, j,&
!!$                                          k_ij, k_ji, istatus)
!!$
!!$!<description>
!!$    ! This subroutine computes the convective matrix coefficients
!!$    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
!!$    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 2D
!!$!</description>
!!$    
!!$!<input>
!!$    ! solution vector
!!$    real(DP), intent(IN) :: u_i,u_j
!!$
!!$    ! coefficients from spatial discretization
!!$    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!!$
!!$    ! nodal indices
!!$    integer, intent(IN) :: i,j
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! OPTIONAL: status of the callback function
!!$    integer, intent(INOUT), optional :: istatus
!!$!</inputoutput>
!!$
!!$!<output>
!!$    ! convective coefficients
!!$    real(DP), intent(OUT) :: k_ij,k_ji
!!$!</output>
!!$!</subroutine>
!!$
!!$    k_ij = -p_DvelocityX(j)*C_ji(1)-p_DvelocityY(j)*C_ji(2)
!!$    k_ji = -p_DvelocityX(i)*C_ij(1)-p_DvelocityY(i)*C_ij(2)
!!$
!!$  end subroutine codire_calcDualConvConst2d
!!$
!!$  !*****************************************************************************
!!$
!!$!<subroutine>
!!$
!!$  pure subroutine codire_calcPrimalConvConst3d(u_i, u_j, C_ij, C_ji, i, j,&
!!$                                            k_ij, k_ji, istatus)
!!$
!!$!<description>
!!$    ! This subroutine computes the convective matrix coefficients
!!$    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
!!$    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 3D
!!$!</description>
!!$    
!!$!<input>
!!$    ! solution vector
!!$    real(DP), intent(IN) :: u_i,u_j
!!$
!!$    ! coefficients from spatial discretization
!!$    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!!$
!!$    ! nodal indices
!!$    integer, intent(IN) :: i,j
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! OPTIONAL: status of the callback function
!!$    integer, intent(INOUT), optional :: istatus
!!$!</inputoutput>
!!$
!!$!<output>
!!$    ! convective coefficients
!!$    real(DP), intent(OUT) :: k_ij,k_ji
!!$!</output>
!!$!</subroutine>
!!$
!!$    k_ij = -p_DvelocityX(j)*C_ij(1)-p_DvelocityY(j)*C_ij(2)-p_DvelocityZ(j)*C_ij(3)
!!$    k_ji = -p_DvelocityX(i)*C_ji(1)-p_DvelocityY(i)*C_ji(2)-p_DvelocityZ(i)*C_ji(3)
!!$
!!$  end subroutine codire_calcPrimalConvConst3d
!!$
!!$  !*****************************************************************************
!!$
!!$!<subroutine>
!!$
!!$  pure subroutine codire_calcDualConvConst3d(u_i, u_j, C_ij, C_ji, i, j,&
!!$                                          k_ij, k_ji, istatus)
!!$
!!$!<description>
!!$    ! This subroutine computes the convective matrix coefficients
!!$    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
!!$    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 3D
!!$!</description>
!!$    
!!$!<input>
!!$    ! solution vector
!!$    real(DP), intent(IN) :: u_i,u_j
!!$
!!$    ! coefficients from spatial discretization
!!$    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!!$
!!$    ! nodal indices
!!$    integer, intent(IN) :: i,j
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! OPTIONAL: status of the callback function
!!$    integer, intent(INOUT), optional :: istatus
!!$!</inputoutput>
!!$
!!$!<output>
!!$    ! convective coefficients
!!$    real(DP), intent(OUT) :: k_ij,k_ji
!!$!</output>
!!$!</subroutine>
!!$
!!$    k_ij = -p_DvelocityX(j)*C_ji(1)-p_DvelocityY(j)*C_ji(2)-p_DvelocityZ(j)*C_ji(3)
!!$    k_ji = -p_DvelocityX(i)*C_ij(1)-p_DvelocityY(i)*C_ij(2)-p_DvelocityZ(i)*C_ij(3)
!!$
!!$  end subroutine codire_calcDualConvConst3d
!!$
!!$  !*****************************************************************************
!!$    
!!$!<subroutine>
!!$
!!$  pure subroutine codire_calcConvectionBurgersSpT2d(u_i, u_j, C_ij, C_ji, i, j,&
!!$                                                 k_ij, k_ji, istatus)
!!$
!!$!<description>
!!$    ! This subroutine computes the convective matrix coefficients
!!$    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the 
!!$    ! one-dimensional Burgers equation $du/dt+df(u)/dx=0$, whereby
!!$    ! the flux function is given by $f(u)=0.5*u^2$.
!!$!</description>
!!$   
!!$!<input>
!!$    ! solution vector
!!$    real(DP), intent(IN) :: u_i,u_j
!!$
!!$    ! coefficients from spatial discretization
!!$    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!!$
!!$    ! nodal indices
!!$    integer, intent(IN) :: i,j
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! OPTIONAL: status of the callback function
!!$    integer, intent(INOUT), optional :: istatus
!!$!</inputoutput>
!!$
!!$!<output>
!!$    ! convective coefficients
!!$    real(DP), intent(OUT) :: k_ij,k_ji
!!$!</output>
!!$!</subroutine>
!!$
!!$    k_ij = -0.5_DP*(u_i+u_j)*C_ij(1)-p_DvelocityY(j)*C_ij(2)
!!$    k_ji = -0.5_DP*(u_i+u_j)*C_ji(1)-p_DvelocityY(i)*C_ji(2)
!!$
!!$  end subroutine codire_calcConvectionBurgersSpT2d
!!$
!!$  !*****************************************************************************
!!$  
!!$!<subroutine>
!!$
!!$  pure subroutine codire_calcConvectionBuckLevSpT2d(u_i, u_j, C_ij, C_ji, i, j,&
!!$                                                 k_ij, k_ji, istatus)
!!$
!!$!<description>
!!$    ! This subroutine computes the convective matrix coefficients
!!$    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the 
!!$    ! Buckley-Leverett equation $du/dt+df(u)/dx=0$, whereby the
!!$    ! flux function is given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
!!$    !
!!$    ! Here, the characteristic velocity $a(u)=f^\prime(u)$ is given
!!$    ! by $a(u)=\frac{4u(1-u)}{(3u^2-2u+1)^2}$.
!!$!</description>
!!$   
!!$!<input>
!!$    ! solution vector
!!$    real(DP), intent(IN) :: u_i,u_j
!!$
!!$    ! coefficients from spatial discretization
!!$    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!!$
!!$    ! nodal indices
!!$    integer, intent(IN) :: i,j
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! OPTIONAL: status of the callback function
!!$    integer, intent(INOUT), optional :: istatus
!!$!</inputoutput>
!!$
!!$!<output>
!!$    ! convective coefficients
!!$    real(DP), intent(OUT) :: k_ij,k_ji
!!$!</output>
!!$!</subroutine>
!!$
!!$    ! local variables
!!$    real(DP) :: v_i,v_j
!!$    
!!$    v_i = 4*u_i*(1-u_i)/(3*u_i*u_i-2*u_i+1)**2
!!$    v_j = 4*u_j*(1-u_j)/(3*u_j*u_j-2*u_j+1)**2
!!$
!!$    k_ij = -v_j*C_ij(1)-p_DvelocityY(j)*C_ij(2)
!!$    k_ji = -v_i*C_ji(1)-p_DvelocityY(i)*C_ji(2)
!!$        
!!$  end subroutine codire_calcConvectionBuckLevSpT2d
!!$
!!$  !*****************************************************************************
!!$    
!!$!<subroutine>
!!$
!!$  pure subroutine codire_calcConvectionBurgers1d(u_i, u_j, C_ij, C_ji, i, j,&
!!$                                              k_ij, k_ji, istatus)
!!$
!!$!<description>
!!$    ! This subroutine computes the convective matrix coefficients
!!$    ! $k_{ij}$ and $k_{ji}$ for Burgers' equation in 1D.
!!$!</description>
!!$   
!!$!<input>
!!$    ! solution vector
!!$    real(DP), intent(IN) :: u_i,u_j
!!$
!!$    ! coefficients from spatial discretization
!!$    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!!$
!!$    ! nodal indices
!!$    integer, intent(IN) :: i,j
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! OPTIONAL: status of the callback function
!!$    integer, intent(INOUT), optional :: istatus
!!$!</inputoutput>
!!$
!!$!<output>
!!$    ! convective coefficients
!!$    real(DP), intent(OUT) :: k_ij,k_ji
!!$!</output>
!!$!</subroutine>
!!$
!!$    k_ij = -0.5_DP*(u_i+u_j)*C_ij(1)
!!$    k_ji = -0.5_DP*(u_i+u_j)*C_ji(1)
!!$
!!$  end subroutine codire_calcConvectionBurgers1d
!!$
!!$  !*****************************************************************************
!!$    
!!$!<subroutine>
!!$
!!$  pure subroutine codire_calcConvectionBurgers2d(u_i, u_j, C_ij, C_ji, i, j,&
!!$                                              k_ij, k_ji, istatus)
!!$
!!$!<description>
!!$    ! This subroutine computes the convective matrix coefficients
!!$    ! $k_{ij}$ and $k_{ji}$ for Burgers' equation in 2D.
!!$!</description>
!!$   
!!$!<input>
!!$    ! solution vector
!!$    real(DP), intent(IN) :: u_i,u_j
!!$
!!$    ! coefficients from spatial discretization
!!$    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!!$
!!$    ! nodal indices
!!$    integer, intent(IN) :: i,j
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! OPTIONAL: status of the callback function
!!$    integer, intent(INOUT), optional :: istatus
!!$!</inputoutput>
!!$
!!$!<output>
!!$    ! convective coefficients
!!$    real(DP), intent(OUT) :: k_ij,k_ji
!!$!</output>
!!$!</subroutine>
!!$
!!$    k_ij = -0.5_DP*(u_i+u_j)*(C_ij(1)+C_ij(2))
!!$    k_ji = -0.5_DP*(u_i+u_j)*(C_ji(1)+C_ji(2))
!!$
!!$  end subroutine codire_calcConvectionBurgers2d
!!$  
!!$  !*****************************************************************************
!!$  
!!$!<subroutine>
!!$
!!$  pure subroutine codire_calcConvectionBuckLev1d(u_i, u_j, C_ij, C_ji, i, j,&
!!$                                              k_ij, k_ji, istatus)
!!$
!!$!<description>
!!$    ! This subroutine computes the convective matrix coefficients
!!$    ! $k_{ij}$ and $k_{ji}$ for the Buckley-Leverett equation 
!!$    ! $du/dt+df(u)/dx=0$ in 1D, whereby the flux function is 
!!$    ! given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
!!$    !
!!$    ! Here, the characteristic velocity $a(u)=f^\prime(u)$ is given
!!$    ! by $a(u)=\frac{4u(1-u)}{(3u^2-2u+1)^2}$.
!!$!</description>
!!$   
!!$!<input>
!!$    ! solution vector
!!$    real(DP), intent(IN) :: u_i,u_j
!!$
!!$    ! coefficients from spatial discretization
!!$    real(DP), dimension(:), intent(IN) :: C_ij,C_ji
!!$
!!$    ! nodal indices
!!$    integer, intent(IN) :: i,j
!!$!</input>
!!$
!!$!<inputoutput>
!!$    ! OPTIONAL: status of the callback function
!!$    integer, intent(INOUT), optional :: istatus
!!$!</inputoutput>
!!$
!!$!<output>
!!$    ! convective coefficients
!!$    real(DP), intent(OUT) :: k_ij,k_ji
!!$!</output>
!!$!</subroutine>
!!$
!!$    k_ij = -(4*u_j*(1-u_j)/(3*u_j*u_j-2*u_j+1)**2)*C_ij(1)
!!$    k_ji = -(4*u_i*(1-u_i)/(3*u_i*u_i-2*u_i+1)**2)*C_ji(1)
!!$    
!!$  end subroutine codire_calcConvectionBuckLev1d
  
end module codire_problem

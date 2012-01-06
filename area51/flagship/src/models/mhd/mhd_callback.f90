!##############################################################################
!# ****************************************************************************
!# <name> mhd_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the compressible MHD equations in arbitrary spatial dimensions.
!#
!# The following callback functions are available:
!#
!# 1.) mhd_nlsolverCallback
!#     -> Callback routine for the nonlinear solver
!#
!# ****************************************************************************
!#
!# The following auxiliary routines are available:
!#
!# 1.) mhd_calcPrecondThetaScheme
!#     -> Calculates the nonlinear preconditioner
!#        used in the two-level theta-scheme
!#
!# 2.) mhd_calcJacobianThetaScheme
!#     -> Calculates the Jacobian matrix
!#        used in the two-level theta-scheme
!#
!# 3.) mhd_calcResidualThetaScheme
!#     -> Calculates the nonlinear residual vector
!#        used in the two-level theta-scheme
!#
!# 4.) mhd_calcRhsThetaScheme
!#     -> Calculates the explicit right-hand side vector
!#        used in the two-level theta-scheme
!#
!# 5.) mhd_calcRhsRungeKuttaScheme
!#     -> Calculates the right-hand side vector
!#        used in the explicit Runge-Kutta scheme
!#
!# 6.) mhd_calcLinearisedFCT
!#     -> Calculates the linearised FCT correction
!#
!# 7.) mhd_calcFluxFCT
!#     -> Calculates the raw antidiffusive fluxes for FCT algorithm
!#
!# 8.) mhd_calcCorrectionFCT
!#     -> Calculates the contribution of the antidiffusive fluxes
!#        limited by the FCT algorithm and applies them to the residual
!#
!# 9.) mhd_limitEdgewiseVelocity
!#     -> Performs synchronised flux correction for the velocity
!#
!# 10.) mhd_limitEdgewiseMomentum
!#      -> Performs synchronised flux correction for the momentum
!#
!# 11.) mhd_limitEdgewiseMagfield
!#      -> Performs synchronised flux correction for the magnetic field
!#
!# 12.) mhd_calcDivergenceVector
!#      -> Calculates the divergence vector.
!#
!# 13.) hydro_calcTimeDerivative
!#      -> Cacluates the approximate time derivative
!#
!# 14.) mhd_coeffVectorFE
!#      -> Callback routine for the evaluation of linear forms
!#         using a given FE-solution for interpolation
!#
!# 15.) mhd_coeffVectorAnalytic
!#      -> Callback routine for the evaluation of linear forms
!#         using a given FE-solution for interpolation
!#
!# 16.) mhd_parseBoundaryCondition
!#      -> Callback routine for the treatment of boundary conditions
!#
!# 17.) mhd_setBoundaryCondition
!#      -> Imposes boundary conditions for nonlinear solver
!#         by filtering the system matrix and the solution/residual
!#         vector explicitly (i.e. strong boundary conditions)
!#
!#
!# Frequently asked questions?
!#
!# 1.) What is the magic behind subroutine 'mhd_nlsolverCallback'?
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

module mhd_callback

#ifdef MHD_NDIM
#undef MHD_NDIM
#endif
#include "mhd.h"

  use afcstabbase
  use afcstabsystem
  use basicgeometry
  use boundary
  use boundarycondaux
  use boundaryfilter
  use cubature
  use collection
  use derivatives
  use domainintegration
  use feevaluation
  use fparser
  use mhd_basic
  use mhd_callback1d
  use mhd_callback2d
  use mhd_callback3d
  use flagship_basic
  use fsystem
  use genoutput
  use groupfemsystem
  use linearalgebra
  use linearformevaluation
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use scalarpde
  use solveraux
  use spatialdiscretisation
  use statistics
  use storage
  use timestepaux
  use triangulation

  implicit none

  private
  public :: mhd_nlsolverCallback
  public :: mhd_calcPrecondThetaScheme
  public :: mhd_calcJacobianThetaScheme
  public :: mhd_calcResidualThetaScheme
  public :: mhd_calcRhsThetaScheme
  public :: mhd_calcRhsRungeKuttaScheme
  public :: mhd_setBoundaryCondition
  public :: mhd_calcLinearisedFCT
  public :: mhd_calcFluxFCT
  public :: mhd_calcCorrectionFCT
  public :: mhd_calcDivergenceVector
  public :: mhd_calcTimeDerivative
  public :: mhd_limitEdgewiseVelocity
  public :: mhd_limitEdgewiseMomentum
  public :: mhd_limitEdgewiseMagfield
  public :: mhd_coeffVectorFE
  public :: mhd_coeffVectorAnalytic
  public :: mhd_parseBoundaryCondition

contains

  !*****************************************************************************

!<subroutine>

  subroutine mhd_nlsolverCallback(rproblemLevel, rtimestep,&
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

    ! local variables
    character(len=SYS_STRLEN) :: ssectionName
    
    ! Get section name
    call collct_getvalue_string(rcollection,&
        'ssectionname', ssectionName)

    ! Do we have to calculate the preconditioner?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) then

      ! Compute the preconditioner
      call mhd_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, ssectionName, rcollection)
    end if


    ! Do we have to calculate the constant right-hand side?
    ! --------------------------------------------------------------------------
    if ((iand(ioperationSpec, NLSOL_OPSPEC_CALCRHS)  .ne. 0)) then

      ! Compute the right-hand side
      call mhd_calcRhsRungeKuttaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rrhs, istep, ssectionName,&
	  rcollection)
    end if


    ! Do we have to calculate the residual?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then
      if (istep .eq. 0) then
        ! Compute the constant right-hand side
        call mhd_calcRhsThetaScheme(rproblemLevel, rtimestep,&
            rsolver, rsolution0, rrhs, ssectionName, rcollection,&
	    rsource)
      end if

      ! Compute the residual
      call mhd_calcResidualThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rrhs, rres, istep,&
	  ssectionName, rcollection)
    end if


    ! Do we have to impose boundary conditions?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      ! Impose boundary conditions
      call mhd_setBoundaryCondition(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rres, ssectionName, rcollection)
    end if


    ! Set status flag
    istatus = 0

  end subroutine mhd_nlsolverCallback

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, ssectionName, rcollection)

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

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    real(DP) :: dscale
    integer :: systemMatrix,lumpedMassMatrix,consistentMassMatrix
    integer :: inviscidAFC,inviscidGFEM
    integer :: isystemCoupling,isystemPrecond,isystemFormat,imasstype,ivar

    ! Start time measurement for matrix evaluation
    p_rtimer => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyMatrix', ssectionName=ssectionName)
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
          'rparlist', ssectionName=ssectionName)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'systemmatrix', systemMatrix)

    !---------------------------------------------------------------------------
    ! Check if fully explicit time-stepping is used
    !---------------------------------------------------------------------------
    if (rtimestep%theta .eq. 0.0_DP) then

      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'isystemformat', isystemFormat)
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'imasstype', imasstype)
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'consistentmassmatrix', consistentMassMatrix)

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
        case default
          call output_line('Empty system matrix is invalid!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecondThetaScheme')
          call sys_halt()
        end select

      case (SYSTEM_BLOCKFORMAT)

        select case(imasstype)
        case (MASS_LUMPED)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, mhd_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(&
                rproblemLevel%Rmatrix(lumpedMassMatrix),&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                1.0_DP, 1.0_DP, .false., .false., .true., .true.)
          end do

        case (MASS_CONSISTENT)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, mhd_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(&
                rproblemLevel%Rmatrix(consistentMassMatrix),&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                1.0_DP, 1.0_DP, .false., .false., .true., .true.)
          end do

        case default
          call output_line('Empty system matrix is invalid!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecondThetaScheme')
          call sys_halt()
        end select

      case default
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecondThetaScheme')
        call sys_halt()
      end select

      ! That is it
      return
    end if

    !---------------------------------------------------------------------------
    ! Assemble divergence operator (for the left-hand side):
    !
    ! (1) without integration by parts:
    !     $$ \int_\Omega w \nabla \cdot {\bf f}(u) {\rm d}{\bf x} $$
    !
    !     then boundary conditions need to be imposed in strong sense
    !     by filtering the system matrix, the solution vector and/or
    !     the residual explicitly.
    !
    ! (2) with integration by parts:
    !     $$ -\int_\Omega \nabla w \cdot {\bf f}(u) {\rm d}{\bf x} $$
    !
    !     with weakly imposed boundary conditions
    !
    !     $$ \int_{\Gamma_-} w {\bf f}(u_0) \cdot {\bf n} {\rm d}{\bf s} $$
    !
    !     imposed at some part of the boundary. At the remaining part
    !     of the boundary nothing is prescribed and the corresponding
    !     boundary integral is built into the bilinear form
    !
    !     $$ \int_{\Gamma_+} w {\bf f}(u) \cdot {\bf n} {\rm d}{\bf s} $$
    !
    !
    ! Remark: By convention, call-back routines assume that the
    !         divergence operator is used on the right-hand
    !         side. Therefore, we habe to multiply it by -1.
    !
    !
    ! Remark: In future versions this routine will assemble both the primal
    !         and the dual preconditioner depending on the variables
    !         primaldual which has to be extracted from the collection.
    !---------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'isystemCoupling', isystemCoupling)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'isystemPrecond', isystemPrecond)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'isystemformat', isystemFormat)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'imasstype', imasstype)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'inviscidGFEM', inviscidGFEM, inviscidAFC)

    ! Compute scaling parameter
    select case (imasstype)
    case (MASS_LUMPED, MASS_CONSISTENT)
      dscale = -rtimestep%theta*rtimestep%dStep
    case default
      dscale = -1.0_DP
    end select

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
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatGalMatD1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiagMatD1d_sim, rcollection=rcollection)

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatGalMatD2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiagMatD2d_sim, rcollection=rcollection)

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatGalMatD3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiagMatD3d_sim, rcollection=rcollection)

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatScDissMatD1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiagMatD1d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatScDissMatD2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiagMatD2d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatScDissMatD3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiagMatD3d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecond')
          call sys_halt()
        end select


      case (DISSIPATION_ROE)

        ! Assemble divergence operator with Roe-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatRoeDissMatD1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiagMatD1d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatRoeDissMatD2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiagMatD2d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatRoeDissMatD3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiagMatD3d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with the Rusanov-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatRusDissMatD1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiagMatD1d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatRusDissMatD2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiagMatD2d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatRusDissMatD3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiagMatD3d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case default
        ! Clear system matrix and apply (lumped) mass matrix only
        call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
      end select


    case (SYSTEM_ALLCOUPLED)

      !-------------------------------------------------------------------------
      ! Assemble full block divergence operator
      !-------------------------------------------------------------------------

      ! What kind of preconditioner is applied?
      select case(isystemPrecond)

      case (DISSIPATION_ZERO)

        ! Assemble divergence operator without dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatGal1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiag1d_sim, rcollection=rcollection)

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatGal2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiag2d_sim, rcollection=rcollection)

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatGal3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiag3d_sim, rcollection=rcollection)

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatScDiss1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiag1d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatScDiss2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiag2d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatScDiss3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiag3d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_ROE)

        ! Assemble divergence operator with Roe-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatRoeDiss1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiag1d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatRoeDiss2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiag2d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatRoeDiss3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiag3d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with the Rusanov-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatRusDiss1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiag1d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatRusDiss2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiag2d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, mhd_calcMatRusDiss3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              mhd_calcMatDiag3d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case default
        ! Clear system matrix and apply (lumped) mass matrix only
        call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
      end select


    case default
      call output_line('Invalid type of flow coupling!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecondThetaScheme')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Assemble the global system operator
    !---------------------------------------------------------------------------

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
        !
        ! Since we have assembled "-L" it suffices to multiply it by "theta*dt"
        ! and add it to the lumped mass matrix.
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)

        call lsyssc_MatrixLinearComb(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rmatrix(systemMatrix),&
            1.0_DP, 1.0_DP, .false., .false., .true., .true.)

      case (MASS_CONSISTENT)

        !-----------------------------------------------------------------------
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_C)-theta*dt*L $
        !
        ! Since we have assembled "-L" it suffices to multiply it by "theta*dt"
        ! and add it to the consistent mass matrix.
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'consistentmassmatrix', consistentMassMatrix)

        call lsyssc_MatrixLinearComb(&
            rproblemLevel%Rmatrix(consistentMassMatrix),&
            rproblemLevel%Rmatrix(systemMatrix),&
            1.0_DP, 1.0_DP, .false., .false., .true., .true.)

      case default

        !-----------------------------------------------------------------------
        ! Use the global operator for steady-state flow
        !
        !   $ A = -L $
        !
        ! Since we have assembled "-L" nothing needs to be done.
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
        !
        ! Since we have assembled "-L" it suffices to multiply it by "theta*dt"
        ! and add it to the lumped mass matrix.
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)

        do ivar = 1, mhd_getNVAR(rproblemLevel)
          call lsyssc_MatrixLinearComb(&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
              1.0_DP, 1.0_DP, .false., .false., .true., .true.)
        end do


      case (MASS_CONSISTENT)

        !-----------------------------------------------------------------------
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_C)-theta*dt*L $
        !
        ! Since we have assembled "-L" it suffices to multiply it by "theta*dt"
        ! and add it to the consistent mass matrix.
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'consistentmassmatrix', consistentMassMatrix)

        do ivar = 1, mhd_getNVAR(rproblemLevel)
          call lsyssc_MatrixLinearComb(&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
              1.0_DP, 1.0_DP, .false., .false., .true., .true.)
        end do


      case default

        !-----------------------------------------------------------------------
        ! Use the global operator for steady-state flow
        !
        !   $ A = -L $
        !
        ! Since we have assembled "-L" nothing needs to be done.
        !-----------------------------------------------------------------------

      end select


    case default
      call output_line('Invalid system format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcPrecondThetaScheme')
      call sys_halt()
    end select


    ! Impose boundary conditions in strong sence (if any)
    call bdrf_filterMatrix(rsolver%rboundaryCondition, &
        rproblemLevel%RmatrixBlock(systemMatrix))

    ! Ok, we updated the (nonlinear) system operator successfully. Now we still
    ! have to link it to the solver hierarchy. This is done recursively.
    call flagship_updateSolverMatrix(rproblemLevel, rsolver,&
        systemMatrix, isystemFormat, UPDMAT_ALL,&
        rproblemLevel%ilev, rproblemLevel%ilev)

    ! Finally, we have to update the content of the solver hierarchy
    call solver_updateContent(rsolver)

    ! Stop time measurement for global operator
    call stat_stopTimer(p_rtimer)

  end subroutine mhd_calcPrecondThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcJacobianThetaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, ssectionName, rcollection)

!<description>
    ! This callback subroutine computes the Jacobian matrix for the
    !  compressible Euler/Navier-Stokes equations
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolution0

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout), target :: rproblemLevel

    ! nonlinear solver structure
    type(t_solver), intent(inout) :: rsolver

    ! current solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    print *, "!!! The calculation of the Jacobian matrix for the !!!"
    print *, "!!! MHD equations has yet not been implemented     !!!"
    stop

  end subroutine mhd_calcJacobianThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcRhsThetaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rrhs, ssectionName, rcollection, rsource)

!<description>
    ! This subroutine computes the constant right-hand side
    !
    !  $$ rhs = M*U^n + (1-\theta) * \Delta t * div F(U^n) + S(U^n) + b.c.`s  $$
    !
    ! where the source term is optional.
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

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

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    type(t_vectorBlock), pointer :: p_rpredictor
    real(DP) :: dscale
    integer :: consistentMassMatrix, lumpedMassMatrix, massMatrix
    integer :: imasstype, iblock, inviscidAFC


    ! Start time measurement for residual/rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyVector', ssectionName=ssectionName)
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection,&
          'rparlist', ssectionName=ssectionName)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'imasstype', imasstype)

    ! Do we have some kind of mass matrix?
    select case(imasstype)
    case (MASS_LUMPED, MASS_CONSISTENT)

      ! Do we have an explicit part?
      if (rtimestep%theta .ne. 1.0_DP) then

        ! Compute scaling parameter
        dscale = (1.0_DP-rtimestep%theta) * rtimestep%dStep

        !-----------------------------------------------------------------------
        ! Compute the divergence operator for the right-hand side
        ! evaluated at the solution from the previous(!) iteration
        !
        !   $$ rhs = (1-\theta) * \Delta t * [div F(U^n) + geomSource(U^n) $$
        !-----------------------------------------------------------------------

        call mhd_calcDivergenceVector(rproblemLevel,&
            rsolver%rboundaryCondition, rsolution,&
            rtimestep%dTime-rtimestep%dStep, dscale, .true.,&
            rrhs, ssectionName, rcollection)

!!$        call mhd_calcGeometricSourceterm(p_rparlist, ssectionName,&
!!$            rproblemLevel, rsolution, dscale, .false., rrhs, rcollection)

        !-----------------------------------------------------------------------
        ! Compute the transient term
        !
        !   $$ rhs := rhs + M*U^n $$
        !-----------------------------------------------------------------------

        ! What type of mass matrix should be used?
        massMatrix = merge(lumpedMassMatrix,&
            consistentMassMatrix, imasstype .eq. MASS_LUMPED)

        ! Apply mass matrix to solution vector
        do iblock = 1, rsolution%nblocks
          call lsyssc_scalarMatVec(&
              rproblemLevel%Rmatrix(massMatrix),&
              rsolution%RvectorBlock(iblock),&
              rrhs%RvectorBlock(iblock), 1.0_DP , 1.0_DP)
        end do

        !-----------------------------------------------------------------------
        ! Perform preparation tasks for algebraic flux correction schemes
        ! of FCT-type which are based on a low-order predictor
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'inviscidAFC', inviscidAFC, 0)

        if (inviscidAFC > 0) then

          ! What type of stabilisation are we?
          select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)

          case (AFCSTAB_NLINFCT_EXPLICIT,&
                AFCSTAB_NLINFCT_ITERATIVE,&
                AFCSTAB_NLINFCT_IMPLICIT)

            ! Compute the low-order predictor based on the right-hand side
            ! and assemble the explicit part of the raw-antidiffusive fluxes

            ! Set pointer to predictor
            p_rpredictor => rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorPredictor

            ! Compute $\tilde u = (M_L)^{-1}*b^n$
            call lsysbl_invertedDiagMatVec(&
                rproblemLevel%Rmatrix(lumpedMassMatrix),&
                rrhs, 1.0_DP, p_rpredictor)

            ! Set specifier
            rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec =&
                ior(rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec,&
                    AFCSTAB_HAS_PREDICTOR)

            ! Assemble explicit part of the raw-antidiffusive fluxes
            call mhd_calcFluxFCT(rproblemLevel, rsolution,&
                rtimestep%theta, rtimestep%dStep, 1.0_DP, .true., .true.,&
                AFCSTAB_FCTFLUX_EXPLICIT, ssectionName, rcollection,&
                rsolutionPredictor=p_rpredictor)
          end select
        end if

      else ! theta = 1

        !-----------------------------------------------------------------------
        ! Compute the transient term
        !
        !   $$ rhs = M*U^n $$
        !-----------------------------------------------------------------------

        ! What type of mass matrix should be used?
        massMatrix = merge(lumpedMassMatrix,&
            consistentMassMatrix, imasstype .eq. MASS_LUMPED)

        ! Apply mass matrix to solution vector
        do iblock = 1, rsolution%nblocks
          call lsyssc_scalarMatVec(&
              rproblemLevel%Rmatrix(massMatrix),&
              rsolution%RvectorBlock(iblock),&
              rrhs%RvectorBlock(iblock), 1.0_DP , 0.0_DP)
        end do

        !-----------------------------------------------------------------------
        ! Perform preparation tasks for algebraic flux correction schemes
        ! of FCT-type which are based on a low-order predictor
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'inviscidAFC', inviscidAFC, 0)

        if (inviscidAFC > 0) then

          ! What type of stabilisation are we?
          select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)

          case (AFCSTAB_NLINFCT_EXPLICIT,&
                AFCSTAB_NLINFCT_ITERATIVE,&
                AFCSTAB_NLINFCT_IMPLICIT)

            ! Compute the low-order predictor based on the right-hand side
            ! and assemble the explicit part of the raw-antidiffusive fluxes

            ! Set pointer to predictor
            p_rpredictor => rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorPredictor

            ! Compute $\tilde u = (M_L)^{-1}*b^n = u^n$
            call lsysbl_copyVector(rsolution, p_rpredictor)

            ! Set specifier
            rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec =&
                ior(rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec,&
                    AFCSTAB_HAS_PREDICTOR)

            ! Assemble explicit part of the raw-antidiffusive fluxes
            call mhd_calcFluxFCT(rproblemLevel, rsolution,&
                rtimestep%theta, rtimestep%dStep, 1.0_DP, .true., .true.,&
                AFCSTAB_FCTFLUX_EXPLICIT, ssectionName, rcollection,&
                rsolutionPredictor=p_rpredictor)
          end select
        end if

      end if ! theta
      
    case default

      !-------------------------------------------------------------------------
      ! Initialise the constant right-hand side by zeros
      !
      !   $$ rhs = 0 $$
      !
      ! Note that there is no explicit part from algebraic flux corretion
      !-------------------------------------------------------------------------

      ! Clear right-hand side vector
      call lsysbl_clearVector(rrhs)
    end select

    ! Apply the source vector to the right-hand side (if any)
    if (present(rsource)) then
      if (rsource%NEQ .gt. 0)&
        call lsysbl_vectorLinearComb(rsource, rrhs, 1.0_DP, 1.0_DP)
    end if

    ! Stop time measurement for rhs evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine mhd_calcRhsThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcResidualThetaScheme(rproblemLevel,&
      rtimestep, rsolver, rsolution, rsolution0, rrhs, rres,&
      ite, ssectionName, rcollection, rsource)

!<description>
    ! This subroutine computes the nonlinear residual vector
    !
    ! $$ res^{(m)} = rhs-[M*U^{(m)}-\theta\Delta t div F(U^{(m)})-S^{(m)}-b.c.`s $$
    !
    ! for the standard two-level theta-scheme, whereby the  source
    ! term $S^{(m)}$ is optional. The constant right-hand side
    !
    !  $$ rhs = [M*U^n + (1-\theta)\Delta t div F(U^n) + S^n + b.c.`s $$
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

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

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

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    type(t_vectorBlock), pointer :: p_rpredictor
    real(DP) :: dscale
    integer(I32) :: ioperationSpec
    integer :: consistentMassMatrix, lumpedMassMatrix, massMatrix
    integer :: massafc, inviscidAFC, viscousAFC, imasstype, iblock


    ! Start time measurement for residual/rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyVector', ssectionName=ssectionName)
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection,&
          'rparlist', ssectionName=ssectionName)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'imasstype', imasstype)

    !-------------------------------------------------------------------------
    ! Initialise the residual by the constant right-hand side
    !
    !   $$ res := rhs $$
    !-------------------------------------------------------------------------
    call lsysbl_copyVector(rrhs, rres)

    ! Do we have some kind of mass matrix?
    select case(imasstype)
    case (MASS_LUMPED, MASS_CONSISTENT)

      ! Compute scaling parameter
      dscale = rtimestep%theta*rtimestep%dStep
      
      !-----------------------------------------------------------------------
      ! Compute the transient term
      !
      !   $$ res := res - M*U^{(m)} $$
      !-----------------------------------------------------------------------
      
      ! What type of mass matrix should be used?
      massMatrix = merge(lumpedMassMatrix,&
          consistentMassMatrix, imasstype .eq. MASS_LUMPED)

      ! Apply mass matrix to solution vector
      do iblock = 1, rsolution%nblocks
        call lsyssc_scalarMatVec(&
            rproblemLevel%Rmatrix(massMatrix),&
            rsolution%RvectorBlock(iblock),&
            rres%RvectorBlock(iblock) , -1._DP, 1.0_DP)
      end do
      
    case default
      
      ! Set scaling parameter
      dscale = 1.0_DP
      
    end select

    !---------------------------------------------------------------------------
    ! Update the residual vector
    !
    !   $$ res := res + dscale * [div F(U^{(m)}) + S(U^{(m)})]$$
    !
    ! where
    !
    !   $dscale = \theta * \Delta t$ for transient flows
    !   $dscale = 1$                 for steady-state flows
    !---------------------------------------------------------------------------

    ! Do we have an implicit part?
    if (dscale .ne. 0.0_DP) then
      ! Compute the implicit part of the divergence term
      call mhd_calcDivergenceVector(rproblemLevel,&
          rsolver%rboundaryCondition, rsolution, rtimestep%dTime,&
          dscale, .false., rres, ssectionName, rcollection)
      
      ! Build the geometric source term (if any)
!!$      call mhd_calcGeometricSourceterm(p_rparlist, ssectionName,&
!!$          rproblemLevel, rsolution, dscale, .false., rres, rcollection)
    end if

    !-------------------------------------------------------------------------
    ! Perform algebraic flux correction for the mass term (if required)
    !
    !   $$ res := res + dscale*fmass(u^(m),u^n) $$
    !-------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'massAFC', massAFC, 0)

    if (massAFC > 0) then
      
      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(massAFC)%cafcstabType)

      case (AFCSTAB_NLINLPT_MASS)
        print *, "AFCSTAB_NLINLPT_MASS not implemented yet"
        stop

      end select
    end if

    !-------------------------------------------------------------------------
    ! Perform algebraic flux correction for the inviscid term (if required)
    !
    !   $$ res := res + dscale*finviscid(u^(m),u^n) $$
    !-------------------------------------------------------------------------

    ! What type if stabilisation is applied?
    select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)
    case (AFCSTAB_NLINFCT_EXPLICIT,&
          AFCSTAB_NLINFCT_ITERATIVE,&
          AFCSTAB_NLINFCT_IMPLICIT)


      ! Set pointer to the predictor vector
      p_rpredictor => rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorPredictor

      ! Set operation specifier
      ioperationSpec = AFCSTAB_FCTFLUX_IMPLICIT
      if (ite .gt. 0)&
          ! This has only influence on iterative FCT algorithm
          ioperationSpec = ioperationSpec + AFCSTAB_FCTFLUX_REJECTED

      ! Assemble implicit part of the raw-antidiffusive fluxes
      call mhd_calcFluxFCT(rproblemLevel, rsolution, rtimestep%theta,&
          rtimestep%dStep, 1.0_DP, .true., .true., ioperationSpec,&
          ssectionName, rcollection, rsolutionPredictor=p_rpredictor)
      
      ! Set operation specifier
      if (ite .eq. 0) then
        ! Perform standard flux correction in zeroth iteration
        ioperationSpec = AFCSTAB_FCTALGO_STANDARD
      else
        select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)
        case (AFCSTAB_NLINFCT_EXPLICIT)
          ! Perform standard flux correction without recomputing bounds
          ioperationSpec = AFCSTAB_FCTALGO_STANDARD-&
                           AFCSTAB_FCTALGO_BOUNDS

        case (AFCSTAB_NLINFCT_IMPLICIT)
          ! Perform semi-implicit flux correction
          ioperationSpec = AFCSTAB_FCTALGO_INITALPHA+&
                           AFCSTAB_FCTALGO_LIMITEDGE+&
                           AFCSTAB_FCTALGO_CORRECT+&
                           AFCSTAB_FCTALGO_CONSTRAIN

        case (AFCSTAB_NLINFCT_ITERATIVE)
          ! Perform standard flux correction
          ioperationSpec = AFCSTAB_FCTALGO_STANDARD
        end select
      end if

      ! Perform flux correction
      call mhd_calcCorrectionFCT(rproblemLevel, p_rpredictor,&
          rtimestep%dStep, .false., ioperationSpec, rres,&
          ssectionName, rcollection)

      ! Special treatment for iterative FCT-algorithm
      if (rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType .eq.&
          AFCSTAB_NLINFCT_ITERATIVE) then
        ! Subtract corrected antidiffusion from right-hand side
        call afcsys_buildVectorFCT(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            p_rpredictor, rtimestep%dStep, .false.,&
            AFCSTAB_FCTALGO_CORRECT, rrhs,&
            rcollection=rcollection)

        ! Recompute the low-order predictor for the next limiting step
        call lsysbl_invertedDiagMatVec(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rrhs, 1.0_DP, p_rpredictor)
      end if

      !-------------------------------------------------------------------------
      ! Remark: Some other algebraic flux correction algorithms which
      ! are not based on the computation of an auxiliary low-order
      ! predictor are implemented in subroutine hydro_calcDivergenceVector
      !-------------------------------------------------------------------------

    end select

    !-------------------------------------------------------------------------
    ! Perform algebraic flux correction for the viscous term (if required)
    !
    !   $$ res = res + dscale*fviscous(u^{(m)},u^n) $$
    !-------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'viscousAFC', viscousAFC, 0)

    if (viscousAFC > 0) then
      
      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(viscousAFC)%cafcstabType)
        
      case (AFCSTAB_NLINLPT_SYMMETRIC)
        print *, "AFCSTAB_NLINLPT_SYMMETRIC not implemented yet"
        stop

      end select
    end if

    ! Apply the source vector to the residual  (if any)
    if (present(rsource)) then
      if (rsource%NEQ .gt. 0)&
      call lsysbl_vectorLinearComb(rsource, rres, -1.0_DP, 1.0_DP)
    end if

    ! Stop time measurement for residual/rhs evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine mhd_calcResidualThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcRhsRungeKuttaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rrhs, istep, ssectionName,&
      rcollection, rsource)

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

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

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

    ! right-hand side vector
    type(t_vectorBlock), intent(inout) :: rrhs

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    type(t_vectorBlock), pointer :: p_rpredictor
    real(DP) :: dscale
    integer :: lumpedMassMatrix, consistentMassMatrix, massMatrix
    integer :: imasstype, iblock, massAFC, inviscidAFC, viscousAFC


    ! Start time measurement for residual/rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyVector', ssectionName=ssectionName)
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection,&
          'rparlist', ssectionName=ssectionName)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'imasstype', imasstype)

    !---------------------------------------------------------------------------
    ! Compute the scaling parameter
    !
    !   $ dscale = weight * \Delta t $
    !---------------------------------------------------------------------------
    
    dscale = rtimestep%DmultistepWeights(istep)*rtimestep%dStep

    !---------------------------------------------------------------------------
    ! Compute the divergence operator for the right-hand side
    ! evaluated at the solution from the previous(!) iteration
    !
    !   $$ rhs = dscale * [div F(U^n) + geomSource(U^n)]$$
    !---------------------------------------------------------------------------

    if (dscale .ne. 0.0_DP) then

      ! Compute the explicit part of the divergence term
      call mhd_calcDivergenceVector(rproblemLevel,&
          rsolver%rboundaryCondition, rsolution,&
          rtimestep%dTime-rtimestep%dStep, dscale, .true.,&
          rrhs, ssectionName, rcollection)

      ! Build the geometric source term (if any)
!!$      call mhd_calcGeometricSourceterm(p_rparlist, ssectionName,&
!!$          rproblemLevel, rsolution, dscale, .false., rrhs, rcollection)
    end if
      
    select case(imasstype)
    case (MASS_LUMPED, MASS_CONSISTENT)

      !-------------------------------------------------------------------------
      ! Compute the transient term
      !
      !   $$ rhs := rhs + M*U^n $$
      !--------------------------------------------------------------------------

      ! What type of mass matrix should be used?
      massMatrix = merge(lumpedMassMatrix,&
          consistentMassMatrix, imasstype .eq. MASS_LUMPED)
      
      ! Apply mass matrix to solution vector
      do iblock = 1, rsolution%nblocks
        call lsyssc_scalarMatVec(&
            rproblemLevel%Rmatrix(massMatrix),&
            rsolution%RvectorBlock(iblock),&
            rrhs%RvectorBlock(iblock), 1.0_DP , 1.0_DP)
      end do
    end select

    !---------------------------------------------------------------------------
    ! Perform algebraic flux correction for the mass term (if required)
    !
    !   $$ rhs := rhs + weight*dt*fmass(u^n+1,u^n) $$
    !--------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'massAFC', massAFC, 0)
    
    if (massAFC > 0) then

      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(massAFC)%cafcstabType)

      case (AFCSTAB_NLINLPT_MASS)
        print *, "AFCSTAB_NLINLPT_MASS not implemented yet"
        stop
      end select
    end if

    !---------------------------------------------------------------------------
    ! Perform algebraic flux correction for the inviscid term (if required)
    !
    !   $$ rhs := rhs + weight*dt*finviscid(u^n+1,u^n) $$
    !---------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'inviscidAFC', inviscidAFC, 0)

    if (inviscidAFC > 0) then

      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)
        
      case (AFCSTAB_NLINFCT_EXPLICIT,&
            AFCSTAB_NLINFCT_IMPLICIT,&
            AFCSTAB_NLINFCT_ITERATIVE)

        ! Set pointer to predictor
        p_rpredictor => rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorPredictor
        
        ! Compute $\tilde u = (M_L)^{-1}*b^n$
        call lsysbl_invertedDiagMatVec(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rrhs, 1.0_DP, p_rpredictor)
        
        ! Set specifier
        rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec =&
            ior(rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec,&
            AFCSTAB_HAS_PREDICTOR)

        ! Assemble explicit part of the raw-antidiffusive fluxes
        call mhd_calcFluxFCT(rproblemLevel, rsolution,&
            rtimestep%theta, rtimestep%dStep, 1.0_DP, .true., .true.,&
            AFCSTAB_FCTFLUX_EXPLICIT, ssectionName, rcollection,&
            rsolutionPredictor=p_rpredictor)

        ! Perform flux correction
        call mhd_calcCorrectionFCT(rproblemLevel, p_rpredictor,&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD, rrhs,&
            ssectionName, rcollection)
      end select
    end if

    !---------------------------------------------------------------------------
    ! Perform algebraic flux correction for the viscous term (if required)
    !
    !   $$ rhs := rhs + weight*dt*fviscous(u^n+1,u^n) $$
    !---------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'viscousAFC', viscousAFC, 0)

    if (viscousAFC > 0) then
      
      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(viscousAFC)%cafcstabType)

      case (AFCSTAB_NLINLPT_SYMMETRIC)
        print *, "AFCSTAB_NLINLPT_SYMMETRIC not implemented yet"
        stop
      end select
    end if

    ! Apply the source vector to the right-hand side (if any)
    if (present(rsource)) then
      if (rsource%NEQ .gt. 0)&
          call lsysbl_vectorLinearComb(rsource, rrhs, 1.0_DP, 1.0_DP)
    end if

    ! Stop time measurement for global operator
    call stat_stopTimer(p_rtimer)

  end subroutine mhd_calcRhsRungeKuttaScheme

  !*****************************************************************************

!<subroutine>

  subroutine mhd_setBoundaryCondition(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rres, ssectionName, rcollection)

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

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName
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
    p_rparlist => collct_getvalue_parlst(rcollection,&
          'rparlist', ssectionName=ssectionName)

    ! What type of nonlinear preconditioner are we?
    select case(rsolver%iprecond)
    case (NLSOL_PRECOND_BLOCKD,&
          NLSOL_PRECOND_DEFCOR,&
          NLSOL_PRECOND_NEWTON_FAILED)

      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'systemmatrix', imatrix)

    case (NLSOL_PRECOND_NEWTON)

      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'jacobianmatrix', imatrix)

    case default
      call output_line('Invalid nonlinear preconditioner!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_setBoundaryCondition')
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
          mhd_calcBoundaryvalues1d, istatus)

    case (NDIM2D)
      call bdrf_filterSolution(&
          rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixBlock(imatrix),&
          rsolution, rres, rsolution0, rtimestep%dTime,&
          mhd_calcBoundaryvalues2d, istatus)

    case (NDIM3D)
      call bdrf_filterSolution(&
          rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixBlock(imatrix),&
          rsolution, rres, rsolution0, rtimestep%dTime,&
          mhd_calcBoundaryvalues3d, istatus)

    case default
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_setBoundaryCondition')
      call sys_halt()
    end select

  end subroutine mhd_setBoundaryCondition

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcLinearisedFCT(rproblemLevel, rtimestep,&
      rsolver, rsolution, ssectionName, rcollection, rsource,&
      rvector1, rvector2, rvector3)

!<description>
    ! This subroutine calculates the linearised FCT correction
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

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

    ! OPTIONAL: auxiliary vectors used to compute the approximation to
    ! the time derivative (if not present, then temporal memory is allocated)
    type(t_vectorBlock), intent(inout), target, optional :: rvector1
    type(t_vectorBlock), intent(inout), target, optional :: rvector2
    type(t_vectorBlock), intent(inout), target, optional :: rvector3
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_vectorBlock), pointer :: p_rvector1
    character(len=SYS_STRLEN), dimension(:), pointer :: SfailsafeVariables
    integer :: inviscidAFC,nfailsafe,ivariable,nvariable
    integer :: imassantidiffusiontype,lumpedMassMatrix
    logical :: bisAccepted


    ! Get parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'inviscidAFC', inviscidAFC, 0)

    !---------------------------------------------------------------------------
    ! Linearised FEM-FCT algorithm for the inviscid term (if any)

    if (inviscidAFC > 0) then

      ! What type of stabilisation are we?
      select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)
        
        case (AFCSTAB_LINFCT)
        ! Get parameters from parameter list
        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)
        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'lumpedmassmatrix', lumpedmassmatrix)
        call parlst_getvalue_int(p_rparlist, ssectionName,&
            'nfailsafe', nfailsafe)
        
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          
          ! Set up vector for computing the approximate time derivative
          if (present(rvector1)) then
            p_rvector1 => rvector1
          else
            allocate(p_rvector1)
          end if
          
          ! Compute approximate time derivative
          call mhd_calcTimeDerivative(rproblemLevel, rtimestep,&
              rsolver, rsolution, ssectionName, rcollection, p_rvector1,&
              rsource, rvector2, rvector3)
          
          ! Build the raw antidiffusive fluxes and include
          ! contribution from the consistent mass matrix
          call mhd_calcFluxFCT(rproblemLevel, rsolution, 0.0_DP,&
              1.0_DP, 1.0_DP, .true., .true., AFCSTAB_FCTFLUX_EXPLICIT,&
              ssectionName, rcollection, rsolutionTimeDeriv=p_rvector1)
          
          ! Release temporal memory
          if (.not.present(rvector1)) then
            call lsysbl_releaseVector(p_rvector1)
            deallocate(p_rvector1)
          end if
      
        else
          
          ! Build the raw antidiffusive fluxes without including
          ! the contribution from consistent mass matrix
          call mhd_calcFluxFCT(rproblemLevel, rsolution, 0.0_DP,&
              1.0_DP, 1.0_DP, .true., .true., AFCSTAB_FCTFLUX_EXPLICIT,&
              ssectionName, rcollection)
        end if
    
        !-----------------------------------------------------------------------
        ! Perform failsafe flux correction (if required)
        !-----------------------------------------------------------------------
        
        if (nfailsafe .gt. 0) then
          
          ! Get number of failsafe variables
          nvariable = max(1,&
              parlst_querysubstrings(p_rparlist,&
              ssectionName, 'sfailsafevariable'))
          
          ! Allocate character array that stores all failsafe variable names
          allocate(SfailsafeVariables(nvariable))
          
          ! Initialise character array with failsafe variable names
          do ivariable = 1, nvariable
            call parlst_getvalue_string(p_rparlist,&
                ssectionName, 'sfailsafevariable',&
                Sfailsafevariables(ivariable), isubstring=ivariable)
          end do
          
          ! Compute FEM-FCT correction
          call mhd_calcCorrectionFCT(rproblemLevel,&
              rsolution, rtimestep%dStep, .false.,&
              AFCSTAB_FCTALGO_STANDARD-&
              AFCSTAB_FCTALGO_CORRECT,&
              rsolution, ssectionName, rcollection)
          
          ! Apply failsafe flux correction
          call afcsys_failsafeFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, rtimestep%dStep, 1e-8_DP,&
              AFCSTAB_FAILSAFEALGO_STANDARD, bisAccepted,&
              nsteps=nfailsafe, CvariableNames=SfailsafeVariables,&
              fcb_extractVariable=mhd_getVariable,&
              rcollection=rcollection)
          
          ! Deallocate temporal memory
          deallocate(SfailsafeVariables)
          
        else
          
          ! Apply linearised FEM-FCT correction
          call mhd_calcCorrectionFCT(rproblemLevel,&
              rsolution, rtimestep%dStep, .false.,&
              AFCSTAB_FCTALGO_STANDARD+&
              AFCSTAB_FCTALGO_SCALEBYMASS,&
              rsolution, ssectionName, rcollection)
        end if

      end select
    end if


    ! Impose boundary conditions for the solution vector
    select case(rproblemLevel%rtriangulation%ndim)
    case (NDIM1D)
      call bdrf_filterVectorExplicit(rsolver%rboundaryCondition,&
          rsolution, rtimestep%dTime, mhd_calcBoundaryvalues1d)

    case (NDIM2D)
      call bdrf_filterVectorExplicit(rsolver%rboundaryCondition,&
          rsolution, rtimestep%dTime, mhd_calcBoundaryvalues2d)

    case (NDIM3D)
      call bdrf_filterVectorExplicit(rsolver%rboundaryCondition,&
          rsolution, rtimestep%dTime, mhd_calcBoundaryvalues3d)
    end select

  end subroutine mhd_calcLinearisedFCT

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcFluxFCT(rproblemLevel, rsolution, theta, tstep, dscale,&
      bclear, bquickAssembly, ioperationSpec, ssectionName, rcollection,&
      rsolutionTimeDeriv, rsolutionPredictor)

!<description>
    ! This subroutine calculates the raw antidiffusive fluxes for
    ! the different FCT algorithms
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! implicitness parameter
    real(DP), intent(in) :: theta

    ! time step size
    real(DP), intent(in) :: tstep

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! Switch for flux assembly
    ! TRUE  : destination flux is cleared before assembly
    ! FALSE : destination flux is no cleared before assembly
    logical, intent(in) :: bclear

    ! Switch for flux assembly
    ! TRUE  : fluxes are not modified externally so that
    !         quicker assembly procedures may be feasible
    ! FALSE : fluxes are truely assembled even if this
    !         leads to an expensive addition of zeros
    logical, intent(in) :: bquickAssembly

    ! Operation specification tag. This is a bitfield coming from an OR
    ! combination of different AFCSTAB_FCTFLUX_xxxx constants and specifies
    ! which operations need to be performed by this subroutine.
    integer(I32), intent(in) :: ioperationSpec

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: approximate time derivative of solution vector
    type(t_vectorBlock), intent(in), optional :: rsolutionTimeDeriv

    ! OPTIONAL: low-order predictor used for prelimiting
    type(t_vectorBlock), intent(in), optional :: rsolutionPredictor
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
    integer :: inviscidAFC, inviscidGFEM, lumpedMassMatrix, consistentMassMatrix
    integer :: idissipationtype, imassantidiffusiontype

    ! Get parameters from parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
          'rparlist', ssectionName=ssectionName)

    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'inviscidGFEM', inviscidGFEM, inviscidAFC)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'lumpedmassmatrix', lumpedmassmatrix)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'consistentmassmatrix', consistentmassmatrix)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'imassantidiffusiontype', imassantidiffusiontype)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'idissipationtype', idissipationtype)


    ! What type of dissipation is applied?
    select case(abs(idissipationtype))

    case (DISSIPATION_SCALAR)

      ! Assemble raw antidiffusive fluxes using scalar dissipation

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTScDiss1d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTScDiss1d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rcollection=rcollection)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTScDiss2d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTScDiss2d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rcollection=rcollection)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTScDiss3d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTScDiss3d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rcollection=rcollection)
        end if
      end select


    case (DISSIPATION_ROE)

      ! Assemble raw antidiffusive fluxes using Roe-type dissipation

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTRoeDiss1d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTRoeDiss1d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rcollection=rcollection)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTRoeDiss2d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTRoeDiss2d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rcollection=rcollection)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTRoeDiss3d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTRoeDiss3d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rcollection=rcollection)
        end if
      end select


    case (DISSIPATION_RUSANOV)

      ! Assemble raw antidiffusive fluxes using th Rusanov-type dissipation

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTRusDiss1d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTRusDiss1d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTRusDiss2d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTRusDiss2d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTRusDiss3d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              theta, tstep, dscale, bclear, bquickAssembly, ioperationSpec,&
              mhd_calcFluxFCTRusDiss3d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        end if
      end select


    case default
      call output_line('Invalid type of dissipation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcFluxFCT')
      call sys_halt()
    end select

  end subroutine mhd_calcFluxFCT

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcCorrectionFCT(rproblemLevel, rsolution, &
      dscale, bclear, ioperationSpec, rresidual, ssectionName,&
      rcollection, rafcstab, slimitingvariableName)

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

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName
    
    ! OPTIONAL: Parameter name of limiting variables in parameter list
    ! If not present, then the default string 'slimitingvariable' is used
    character(len=*), intent(in), optional :: slimitingvariableName
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
    p_rparlist => collct_getvalue_parlst(rcollection,&
          'rparlist', ssectionName=ssectionName)

    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'lumpedmassmatrix', lumpedmassmatrix)
    
    ! Set pointer to stabilisation structure
    if (present(rafcstab)) then
      p_rafcstab => rafcstab
    else
      call parlst_getvalue_int(p_rparlist, ssectionName,&
          'inviscidAFC', inviscidAFC)
      p_rafcstab => rproblemLevel%Rafcstab(inviscidAFC)
    end if
    
    ! Get number of limiting variables
    if (present(slimitingvariableName)) then
      nvariable = max(1, parlst_querysubstrings(p_rparlist,&
          ssectionName, slimitingvariableName))
    else
      nvariable = max(1, parlst_querysubstrings(p_rparlist,&
          ssectionName, 'slimitingvariable'))
    end if
    
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
      if (present(slimitingvariableName)) then
        call parlst_getvalue_string(p_rparlist,&
            ssectionName, slimitingvariableName,&
            slimitingvariable, isubstring=ivariable)
      else
        call parlst_getvalue_string(p_rparlist,&
            ssectionName, 'slimitingvariable',&
            slimitingvariable, isubstring=ivariable)
      end if

      ! Get number of variables to be limited simultaneously
      nvartransformed = mhd_getNVARtransformed(rproblemLevel, slimitingvariable)

      ! What type of flux transformation is applied?
      if (trim(slimitingvariable) .eq. 'density') then

        ! Apply FEM-FCT algorithm for density fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxDensity1d_sim, mhd_trafoDiffDensity1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxDensity2d_sim, mhd_trafoDiffDensity2d_sim,&
              rcollection=rcollection)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxDensity3d_sim, mhd_trafoDiffDensity3d_sim,&
              rcollection=rcollection)
        end select

      elseif (trim(slimitingvariable) .eq. 'energy') then

        ! Apply FEM-FCT algorithm for energy fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxEnergy1d_sim, mhd_trafoDiffEnergy1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxEnergy2d_sim, mhd_trafoDiffEnergy2d_sim,&
              rcollection=rcollection)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxEnergy3d_sim, mhd_trafoDiffEnergy3d_sim,&
              rcollection=rcollection)
        end select

      elseif (trim(slimitingvariable) .eq. 'pressure') then

        ! Apply FEM-FCT algorithm for pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxPressure1d_sim, mhd_trafoDiffPressure1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxPressure2d_sim, mhd_trafoDiffPressure2d_sim,&
              rcollection=rcollection)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxPressure3d_sim, mhd_trafoDiffPressure3d_sim,&
              rcollection=rcollection)
        end select

      elseif (trim(slimitingvariable) .eq. 'velocity') then

        ! Apply FEM-FCT algorithm for velocity fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxVelocity1d_sim, mhd_trafoDiffVelocity1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxVelocity2d_sim, mhd_trafoDiffVelocity2d_sim,&
              rcollection=rcollection,&
              fcb_limitEdgewise=mhd_limitEdgewiseVelocity)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxVelocity3d_sim, mhd_trafoDiffVelocity3d_sim,&
              rcollection=rcollection,&
              fcb_limitEdgewise=mhd_limitEdgewiseVelocity)
        end select

      elseif (trim(slimitingvariable) .eq. 'momentum') then

        ! Apply FEM-FCT algorithm for momentum fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxMomentum1d_sim, mhd_trafoDiffMomentum1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxMomentum2d_sim, mhd_trafoDiffMomentum2d_sim,&
              rcollection=rcollection,&
              fcb_limitEdgewise=mhd_limitEdgewiseMomentum)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxMomentum3d_sim, mhd_trafoDiffMomentum3d_sim,&
	      rcollection=rcollection,&
              fcb_limitEdgewise=mhd_limitEdgewiseMomentum)
        end select

      elseif (trim(slimitingvariable) .eq. 'magneticfield') then

        ! Apply FEM-FCT algorithm for magnetic field fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxMagfield1d_sim, mhd_trafoDiffMagfield1d_sim,&
	      rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxMagfield2d_sim, mhd_trafoDiffMagfield2d_sim,&
	      rcollection=rcollection,&
              fcb_limitEdgewise=mhd_limitEdgewiseMomentum)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxMagfield3d_sim, mhd_trafoDiffMagfield3d_sim,&
	      rcollection=rcollection,&
              fcb_limitEdgewise=mhd_limitEdgewiseMomentum)
        end select
        
      elseif (trim(slimitingvariable) .eq. 'density,energy') then

        ! Apply FEM-FCT algorithm for density and energy fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxDenEng1d_sim, mhd_trafoDiffDenEng1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxDenEng2d_sim, mhd_trafoDiffDenEng2d_sim,&
              rcollection=rcollection)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxDenEng3d_sim, mhd_trafoDiffDenEng3d_sim,&
              rcollection=rcollection)
        end select

      elseif (trim(slimitingvariable) .eq. 'density,pressure') then

        ! Apply FEM-FCT algorithm for density and pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxDenPre1d_sim, mhd_trafoDiffDenPre1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxDenPre2d_sim, mhd_trafoDiffDenPre2d_sim,&
              rcollection=rcollection)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxDenPre3d_sim, mhd_trafoDiffDenPre3d_sim,&
              rcollection=rcollection)
        end select

      elseif (trim(slimitingvariable) .eq. 'density,energy,momentum') then

        nvartransformed = mhd_getNVARtransformed(rproblemLevel, slimitingvariable)

        ! Apply FEM-FCT algorithm for full conservative fluxes
        call afcsys_buildVectorFCT(&
            p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rsolution, dscale, bclear, iopSpec, rresidual, rcollection=rcollection)

      elseif (trim(slimitingvariable) .eq. 'density,pressure,velocity') then

        nvartransformed = mhd_getNVARtransformed(rproblemLevel, slimitingvariable)

        ! Apply FEM-FCT algorithm for density, velocity and pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxDenPreVel1d_sim, mhd_trafoDiffDenPreVel1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxDenPreVel2d_sim, mhd_trafoDiffDenPreVel2d_sim,&
              rcollection=rcollection)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              mhd_trafoFluxDenPreVel3d_sim, mhd_trafoDiffDenPreVel3d_sim,&
              rcollection=rcollection)
        end select

      elseif (trim(slimitingvariable) .eq. 'none') then
        
        if (nvariable .eq. 1) then
          ! Apply raw antidiffusive fluxes without correction
          iopSpec = ioperationSpec
          iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_ADINCREMENTS))
          iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_BOUNDS))
          iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_LIMITNODAL))
          iopSpec = iand(iopSpec, not(AFCSTAB_FCTALGO_LIMITEDGE))
          
          ! Enforce existence of edgewise correction factors
          p_rafcstab%istabilisationSpec =&
              ior(p_rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGELIMITER)
          
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, rcollection=rcollection)
          
          ! Nothing more needs to be done
          return
        else
          cycle
        end if

      else
        call output_line('Invalid type of flux transformation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcCorrectionFCT')
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
      
      call afcsys_buildVectorFCT(&
          p_rafcstab, rproblemLevel%Rmatrix(lumpedMassMatrix),&
          rsolution, dscale, bclear, iopSpec, rresidual, rcollection=rcollection)
    end if

  end subroutine mhd_calcCorrectionFCT

  ! ***************************************************************************

!<subroutine>

  subroutine mhd_limitEdgewiseVelocity(IedgeListIdx, IedgeList,&
      NEDGE, NEQ, NVAR, NVARtransformed, ndim1, ndim2, Dx, Dflux, Dalpha,&
      Drp, Drm, fcb_calcFluxTransformation_sim, DfluxConstr, rcollection)

!<description>
    ! This subroutine computes the edgewise correction factors
    ! for the velocity vector in synchronised fashion.
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

    ! Index pointer for edge data structure
    integer, dimension(:), intent(in) :: IedgeListIdx

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IedgeList

    ! OPTIONAL: callback function to compute variable transformation
    include '../../../../../kernel/PDEOperators/intf_calcFluxTransformation_sim.inc'
    optional :: fcb_calcFluxTransformation_sim

    ! OPTIONAL: Antidiffusive flux for constraining
    real(DP), dimension(NVAR,NEDGE), intent(in), optional :: DfluxConstr
!</intput>

!<inputoutput>
    ! Edgewise correction factors
    real(DP), dimension(:), intent(inout) :: Dalpha

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! auxiliary arras
    real(DP), dimension(:,:,:), pointer :: DdataAtEdge
    real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge

    ! local variables
    real(DP), dimension(NVARtransformed) :: R_ij,R_ji,Uij
    real(DP) :: alpha_ij
    integer :: idx,IEDGEset,IEDGEmax,i,j,iedge
    
    ! Do we have to use the explicit fluxes as constraints?
    if (present(DfluxConstr)) then
      print *, "Not implemented yet"
      stop
      
    else

      if ((ndim1 .eq. NVAR) .and. (ndim2 .eq. NEQ)) then

        !-----------------------------------------------------------------------
        ! The vector is given in interleave format
        !-----------------------------------------------------------------------

        ! Allocate temporal memory
        allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
        allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

        ! Loop over the edges
        do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM
          
          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          ! Use callback function to compute transformed fluxes
          call fcb_calcFluxTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
              Dflux(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
              DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1))

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Compute nodal correction factors
            R_ij = merge(Drp(:,i), Drm(:,i),&
                         DtransformedFluxesAtEdge(:,1,idx) .ge. 0.0_DP)
            R_ji = merge(Drp(:,j), Drm(:,j),&
                         DtransformedFluxesAtEdge(:,2,idx) .ge. 0.0_DP)
            
            ! Compute edgewise correction factors
            R_ij = min(R_ij, R_ji)
            
            ! Compute velocity average
            Uij = 0.5_DP*(Dx(2:NVARtransformed+1,i)/Dx(1,i)+&
                          Dx(2:NVARtransformed+1,j)/Dx(1,j))
            
            ! Compute correction factor
            alpha_ij = sum(R_ij*Uij*Uij)/(sum(Uij*Uij)+SYS_EPSREAL_DP)
            
            ! Compute multiplicative correction factor
            Dalpha(iedge) = Dalpha(iedge) *alpha_ij
          end do
        end do

        ! Deallocate temporal memory
        deallocate(DdataAtEdge)
        deallocate(DtransformedFluxesAtEdge)
        
      elseif ((ndim1 .eq. NEQ) .and. (ndim2 .eq. NVAR)) then

        !-----------------------------------------------------------------------
        ! The vector is given in block format
        !-----------------------------------------------------------------------

        ! Allocate temporal memory
        allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
        allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

        ! Loop over the edges
        do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM
          
          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          ! Use callback function to compute transformed fluxes
          call fcb_calcFluxTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
              Dflux(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
              DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1))
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Compute nodal correction factors
            R_ij = merge(Drp(:,i), Drm(:,i),&
                         DtransformedFluxesAtEdge(:,1,idx) .ge. 0.0_DP)
            R_ji = merge(Drp(:,j), Drm(:,j),&
                         DtransformedFluxesAtEdge(:,2,idx) .ge. 0.0_DP)
            
            ! Compute edgewise correction factors
            R_ij = min(R_ij, R_ji)

            ! Compute velocity average
            Uij = 0.5_DP*(Dx(i,2:NVARtransformed+1)/Dx(i,1)+&
                          Dx(j,2:NVARtransformed+1)/Dx(j,1))
            
            ! Compute correction factor
            alpha_ij = sum(R_ij*Uij*Uij)/(sum(Uij*Uij)+SYS_EPSREAL_DP)
            
            ! Compute multiplicative correction factor
            Dalpha(iedge) = Dalpha(iedge) *alpha_ij
          end do
        end do
        
        ! Deallocate temporal memory
        deallocate(DdataAtEdge)
        deallocate(DtransformedFluxesAtEdge)
        
      else
        
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mhd_limitEdgewiseVelocity')
        call sys_halt()
        
      end if
    end if
    
  end subroutine mhd_limitEdgewiseVelocity

  ! ***************************************************************************

!<subroutine>

  subroutine mhd_limitEdgewiseMomentum(IedgeListIdx, IedgeList,&
      NEDGE, NEQ, NVAR, NVARtransformed, ndim1, ndim2, Dx, Dflux, Dalpha,&
      Drp, Drm, fcb_calcFluxTransformation_sim, DfluxConstr, rcollection)

!<description>
    ! This subroutine computes the edgewise correction factors
    ! for the momentum vector in synchronised fashion.
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

    ! Index pointer for edge data structure
    integer, dimension(:), intent(in) :: IedgeListIdx

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IedgeList

    ! OPTIONAL: callback function to compute variable transformation
    include '../../../../../kernel/PDEOperators/intf_calcFluxTransformation_sim.inc'
    optional :: fcb_calcFluxTransformation_sim

    ! OPTIONAL: Antidiffusive flux for constraining
    real(DP), dimension(NVAR,NEDGE), intent(in), optional :: DfluxConstr
!</intput>

!<inputoutput>
    ! Edgewise correction factors
    real(DP), dimension(:), intent(inout) :: Dalpha

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! auxiliary arras
    real(DP), dimension(:,:,:), pointer :: DdataAtEdge
    real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge
    
    ! local variables
    real(DP), dimension(NVARtransformed) :: R_ij,R_ji,Uij
    real(DP) :: alpha_ij
    integer :: idx,IEDGEset,IEDGEmax,i,j,iedge
    
    ! Do we have to use the explicit fluxes as constraints?
    if (present(DfluxConstr)) then
      print *, "Not implemented yet"
      stop
      
    else

      if ((ndim1 .eq. NVAR) .and. (ndim2 .eq. NEQ)) then

        !-----------------------------------------------------------------------
        ! The vector is given in interleave format
        !-----------------------------------------------------------------------

        ! Allocate temporal memory
        allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
        allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

        ! Loop over the edges
        do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM
          
          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          ! Use callback function to compute transformed fluxes
          call fcb_calcFluxTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
              Dflux(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
              DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1))

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Compute nodal correction factors
            R_ij = merge(Drp(:,i), Drm(:,i),&
                         DtransformedFluxesAtEdge(:,1,idx) .ge. 0.0_DP)
            R_ji = merge(Drp(:,j), Drm(:,j),&
                         DtransformedFluxesAtEdge(:,2,idx) .ge. 0.0_DP)
            
            ! Compute edgewise correction factors
            R_ij = min(R_ij, R_ji)
            
            ! Compute momentum average
            Uij = 0.5_DP*(Dx(2:NVARtransformed+1,i)+&
                          Dx(2:NVARtransformed+1,j))
            
            ! Compute correction factor
            alpha_ij = sum(R_ij*Uij*Uij)/(sum(Uij*Uij)+SYS_EPSREAL_DP)
            
            ! Compute multiplicative correction factor
            Dalpha(iedge) = Dalpha(iedge) *alpha_ij
          end do
        end do

        ! Deallocate temporal memory
        deallocate(DdataAtEdge)
        deallocate(DtransformedFluxesAtEdge)

      elseif ((ndim1 .eq. NEQ) .and. (ndim2 .eq. NVAR)) then

        !-----------------------------------------------------------------------
        ! The vector is goven in block format
        !-----------------------------------------------------------------------

        ! Allocate temporal memory
        allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
        allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

        ! Loop over the edges
        do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM
          
          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          ! Use callback function to compute transformed fluxes
          call fcb_calcFluxTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
              Dflux(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
              DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1))
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Compute nodal correction factors
            R_ij = merge(Drp(:,i), Drm(:,i),&
                         DtransformedFluxesAtEdge(:,1,idx) .ge. 0.0_DP)
            R_ji = merge(Drp(:,j), Drm(:,j),&
                         DtransformedFluxesAtEdge(:,2,idx) .ge. 0.0_DP)
            
            ! Compute edgewise correction factors
            R_ij = min(R_ij, R_ji)

            ! Compute momentum average
            Uij = 0.5_DP*(Dx(i,2:NVARtransformed+1)+&
                          Dx(j,2:NVARtransformed+1))
            
            ! Compute correction factor
            alpha_ij = sum(R_ij*Uij*Uij)/(sum(Uij*Uij)+SYS_EPSREAL_DP)
            
            ! Compute multiplicative correction factor
            Dalpha(iedge) = Dalpha(iedge) *alpha_ij
          end do
        end do
        
        ! Deallocate temporal memory
        deallocate(DdataAtEdge)
        deallocate(DtransformedFluxesAtEdge)

      else

        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mhd_limitEdgewiseMomentum')
        call sys_halt()
        
      end if
    end if

  end subroutine mhd_limitEdgewiseMomentum

  ! ***************************************************************************

!<subroutine>

  subroutine mhd_limitEdgewiseMagfield(IedgeListIdx, IedgeList,&
      NEDGE, NEQ, NVAR, NVARtransformed, ndim1, ndim2, Dx, Dflux, Dalpha,&
      Drp, Drm, fcb_calcFluxTransformation_sim, DfluxConstr, rcollection)

!<description>
    ! This subroutine computes the edgewise correction factors
    ! for the magnetic field in synchronised fashion.
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

    ! Index pointer for edge data structure
    integer, dimension(:), intent(in) :: IedgeListIdx

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IedgeList

    ! OPTIONAL: callback function to compute variable transformation
    include '../../../../../kernel/PDEOperators/intf_calcFluxTransformation_sim.inc'
    optional :: fcb_calcFluxTransformation_sim

    ! OPTIONAL: Antidiffusive flux for constraining
    real(DP), dimension(NVAR,NEDGE), intent(in), optional :: DfluxConstr
!</intput>

!<inputoutput>
    ! Edgewise correction factors
    real(DP), dimension(:), intent(inout) :: Dalpha

    ! OPTIONAL: collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>
!</subroutine>

    ! auxiliary arras
    real(DP), dimension(:,:,:), pointer :: DdataAtEdge
    real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge
    
    ! local variables
    real(DP), dimension(NVARtransformed) :: R_ij,R_ji,Uij
    real(DP) :: alpha_ij
    integer :: idx,IEDGEset,IEDGEmax,i,j,iedge
    
    ! Do we have to use the explicit fluxes as constraints?
    if (present(DfluxConstr)) then
      print *, "Not implemented yet"
      stop
      
    else

      if ((ndim1 .eq. NVAR) .and. (ndim2 .eq. NEQ)) then

        !-----------------------------------------------------------------------
        ! The vector is given in interleave format
        !-----------------------------------------------------------------------

        ! Allocate temporal memory
        allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
        allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

        ! Loop over the edges
        do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM
          
          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(:,IedgeList(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IedgeList(2,iedge))
          end do

          ! Use callback function to compute transformed fluxes
          call fcb_calcFluxTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
              Dflux(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
              DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1))

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Compute nodal correction factors
            R_ij = merge(Drp(:,i), Drm(:,i),&
                         DtransformedFluxesAtEdge(:,1,idx) .ge. 0.0_DP)
            R_ji = merge(Drp(:,j), Drm(:,j),&
                         DtransformedFluxesAtEdge(:,2,idx) .ge. 0.0_DP)
            
            ! Compute edgewise correction factors
            R_ij = min(R_ij, R_ji)
            
            ! Compute average of magnetic field
            Uij = 0.5_DP*(Dx(5:NVARtransformed+4,i)+&
                          Dx(5:NVARtransformed+4,j))
            
            ! Compute correction factor
            alpha_ij = sum(R_ij*Uij*Uij)/(sum(Uij*Uij)+SYS_EPSREAL_DP)
            
            ! Compute multiplicative correction factor
            Dalpha(iedge) = Dalpha(iedge) *alpha_ij
          end do
        end do

        ! Deallocate temporal memory
        deallocate(DdataAtEdge)
        deallocate(DtransformedFluxesAtEdge)

      elseif ((ndim1 .eq. NEQ) .and. (ndim2 .eq. NVAR)) then

        !-----------------------------------------------------------------------
        ! The vector is goven in block format
        !-----------------------------------------------------------------------

        ! Allocate temporal memory
        allocate(DdataAtEdge(NVAR,2,GFSYS_NEDGESIM))
        allocate(DtransformedFluxesAtEdge(NVARtransformed,2,GFSYS_NEDGESIM))

        ! Loop over the edges
        do IEDGEset = 1, NEDGE, GFSYS_NEDGESIM
          
          ! We always handle GFSYS_NEDGESIM edges simultaneously.
          ! How many edges have we actually here?
          ! Get the maximum edge number, such that we handle
          ! at most GFSYS_NEDGESIM edges simultaneously.
          
          IEDGEmax = min(NEDGE, IEDGEset-1+GFSYS_NEDGESIM)
          
          ! Loop through all edges in the current set
          ! and prepare the auxiliary arrays
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Fill auxiliary arrays
            DdataAtEdge(:,1,idx) = Dx(IedgeList(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IedgeList(2,iedge),:)
          end do

          ! Use callback function to compute transformed fluxes
          call fcb_calcFluxTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1), &
              Dflux(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
              DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1))
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IedgeList(1,iedge)
            j = IedgeList(2,iedge)

            ! Compute nodal correction factors
            R_ij = merge(Drp(:,i), Drm(:,i),&
                         DtransformedFluxesAtEdge(:,1,idx) .ge. 0.0_DP)
            R_ji = merge(Drp(:,j), Drm(:,j),&
                         DtransformedFluxesAtEdge(:,2,idx) .ge. 0.0_DP)
            
            ! Compute edgewise correction factors
            R_ij = min(R_ij, R_ji)

            ! Compute average of magnetic field
            Uij = 0.5_DP*(Dx(i,5:NVARtransformed+4)+&
                          Dx(j,5:NVARtransformed+4))
            
            ! Compute correction factor
            alpha_ij = sum(R_ij*Uij*Uij)/(sum(Uij*Uij)+SYS_EPSREAL_DP)
            
            ! Compute multiplicative correction factor
            Dalpha(iedge) = Dalpha(iedge) *alpha_ij
          end do
        end do
        
        ! Deallocate temporal memory
        deallocate(DdataAtEdge)
        deallocate(DtransformedFluxesAtEdge)

      else

        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mhd_limitEdgewiseMomentum')
        call sys_halt()
        
      end if
    end if

  end subroutine mhd_limitEdgewiseMagfield

  ! ***************************************************************************

!<subroutine>

  subroutine mhd_coeffVectorFE(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

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
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_coeffVectorFE')
      call sys_halt()
    end select
    
  end subroutine mhd_coeffVectorFE

  ! ***************************************************************************

!<subroutine>

  subroutine mhd_coeffVectorAnalytic(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

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
    
    
    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))
    
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

        ! Initialise values
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
    
  end subroutine mhd_coeffVectorAnalytic

  ! *****************************************************************************

!<subroutine>

  subroutine mhd_parseBoundaryCondition(cbdrCondType, ndimension,&
      ibdrCondType, nexpressions)

!<description>
    ! This subroutine parses the boundarry condition defined in the
    ! parameter file. That is, the string cbdrCondType is tansformed
    ! into an integer value ibdrCondType and the number of
    ! mathematical expressions corresponding to the given boundary
    ! type and the spatial dimension ndimension are returned.
!</description>

!<input>
    ! character string: type of boundary conditions
    character(len=*), intent(in) :: cbdrCondType

    ! number of spatial dimensions
    integer, intent(in) :: ndimension
!</input>

!<output>
    ! type of boundary condition
    integer, intent(out) :: ibdrCondType

    ! number of mathematical expressions
    integer, intent(out) :: nexpressions
!</outpu>
!</subroutine>

    ! Determine type of boundary condition in numeral form
    select case (sys_upcase(cbdrCondType))

    ! Supersonic outlet boundary conditions
    case ('SUPEROUTLET_STRONG')
      ibdrCondType = BDRC_SUPEROUTLET
      ! No strong boundary conditions are prescribed
    case ('SUPEROUTLET_WEAK')
      ibdrCondType = BDRC_SUPEROUTLET + BDRC_WEAK
    case ('SUPEROUTLET_WEAK_LUMPED')
      ibdrCondType = BDRC_SUPEROUTLET + BDRC_WEAK + BDRC_LUMPED
      
    ! Periodic boundary conditions
    case ('PERIODIC_STRONG')
      ibdrCondType = BDRC_PERIODIC + BDRC_STRONG
    case ('PERIODIC_WEAK')
      ibdrCondType = BDRC_PERIODIC + BDRC_WEAK
    case ('PERIODIC_WEAK_LUMPED')
      ibdrCondType = BDRC_PERIODIC + BDRC_WEAK + BDRC_LUMPED

    ! Antiperiodic boundary conditions
    case ('ANTIPERIODIC_STRONG')
      ibdrCondType = BDRC_ANTIPERIODIC + BDRC_STRONG
    case ('ANTIPERIODIC_WEAK')
      ibdrCondType = BDRC_ANTIPERIODIC + BDRC_WEAK
    case ('ANTIPERIODIC_WEAK_LUMPED')
      ibdrCondType = BDRC_ANTIPERIODIC + BDRC_WEAK + BDRC_LUMPED

    case default
      call output_line('Invalid type of boundary conditions!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_parseBoundaryCondition')
      call sys_halt()
    end select

    
    ! Determine number of mathematical expressions
    select case (iand(ibdrCondType, BDRC_TYPEMASK))

    case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      nexpressions = 2 ! penalty parameter $\epsilon$
                       ! switching parameter $\gamma$

    case default
      nexpressions = 0
    end select
    
  end subroutine mhd_parseBoundaryCondition

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcDivergenceVector(rproblemLevel, rboundaryCondition,&
      rsolution, dtime, dscale, bclear, rvector, ssectionName, rcollection)

!<description>
    ! This subroutine computes the discrete divergence vector
    !
    !  $$ div F(u) + b.c.`s  $$
    !
    ! where the (scaled) source term is optional.
!</description>

!<input>
    ! boundary condition
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! Switch for vector assembly
    ! TRUE  : clear vector before assembly
    ! FLASE : assemble vector in an additive way
    logical, intent(in) :: bclear

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! destination vector for the divergence vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    integer :: inviscidAFC,inviscidGFEM,idissipationtype

    ! Set pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)

    ! Get parameter from parameter list
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'inviscidAFC', inviscidAFC, 0)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'inviscidGFEM', inviscidGFEM, inviscidAFC)

    ! Do we have a zero scling parameter?
    if (dscale .eq. 0.0_DP) then
      if (bclear) call lsysbl_clearVector(rvector)
    else
      
      ! Check if group finite element structure and stabilisation
      ! structure are both available
      if ((inviscidGFEM .le. 0) .or. (inviscidAFC .le. 0)) return

      ! What type if stabilisation is applied?
      select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)
        
      case (AFCSTAB_GALERKIN)
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildVectorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution,&
              mhd_calcFluxGal1d_sim, dscale, bclear, rvector, rcollection)
          
        case (NDIM2D)
          call gfsys_buildVectorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution,&
              mhd_calcFluxGal2d_sim, dscale, bclear, rvector, rcollection)
          
        case (NDIM3D)
          call gfsys_buildVectorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution,&
              mhd_calcFluxGal3d_sim, dscale, bclear, rvector, rcollection)
        end select

        !-----------------------------------------------------------------------
        
      case (AFCSTAB_UPWIND,&
            AFCSTAB_NLINFCT_EXPLICIT,&
            AFCSTAB_NLINFCT_ITERATIVE,&
            AFCSTAB_NLINFCT_IMPLICIT,&
            AFCSTAB_LINFCT,&
            AFCSTAB_NLINLPT_UPWINDBIASED,&
            AFCSTAB_LINLPT_UPWINDBIASED)

        ! Get parameter from parameter list
        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'idissipationtype', idissipationtype)
        
        ! What type of dissipation is applied?
        select case(idissipationtype)
          
        case (DISSIPATION_ZERO)
          
          ! Assemble divergence of flux without dissipation
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxGal1d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxGal2d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxGal3d_sim, dscale, bclear, rvector, rcollection)
          end select
          
          !---------------------------------------------------------------------

        case (DISSIPATION_SCALAR)
          
          ! Assemble divergence of flux with scalar dissipation
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxScDiss1d_sim, dscale, bclear , rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxScDiss2d_sim, dscale, bclear, rvector, rcollection)

          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxScDiss3d_sim, dscale, bclear, rvector, rcollection)
          end select
          
          !---------------------------------------------------------------------

        case (DISSIPATION_SCALAR_DSPLIT)
          
          ! Assemble divergence of flux with scalar dissipation
          ! adopting dimensional splitting
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxScDiss1d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxScDissDiSp2d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxScDissDiSp3d_sim, dscale, bclear, rvector, rcollection)
          end select
          
          !---------------------------------------------------------------------

        case (DISSIPATION_ROE)
          
          ! Assemble divergence of flux with Roe-type dissipation
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxRoeDiss1d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxRoeDiss2d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxRoeDiss3d_sim, dscale, bclear, rvector, rcollection)
          end select
          
          !---------------------------------------------------------------------

        case (DISSIPATION_ROE_DSPLIT)
          
          ! Assemble divergence of flux with Roe-type dissipation
          ! adopting dimensional splitting
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxRoeDiss1d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxRoeDissDiSp2d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxRoeDissDiSp3d_sim, dscale, bclear, rvector, rcollection)
          end select
          
          !---------------------------------------------------------------------

        case (DISSIPATION_RUSANOV)
          
          ! Assemble divergence of flux with Rusanov-type flux
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxRusDiss1d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxRusDiss2d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxRusDiss3d_sim, dscale, bclear, rvector, rcollection)
          end select
          
          !---------------------------------------------------------------------

        case (DISSIPATION_RUSANOV_DSPLIT)
          
          ! Assemble divergence of flux with Rusanov-type flux
          ! adopting dimensional splitting
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxRusDiss1d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxRusDissDiSp2d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                mhd_calcFluxRusDissDiSp3d_sim, dscale, bclear, rvector, rcollection)
          end select
          
        case default
          call output_line('Invalid type of dissipation!',&
              OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcDivergenceVector')
          call sys_halt()
        end select
        
        !-----------------------------------------------------------------------

      case (AFCSTAB_TVD)
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorTVD(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, NDIM1D,&
              mhd_calcFluxGalNoBdr1d_sim,&
              mhd_calcCharacteristics1d_sim, dscale, bclear, rvector, rcollection)
          
        case (NDIM2D)
          call afcsys_buildVectorTVD(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, NDIM2D,&
              mhd_calcFluxGalNoBdr2d_sim,&
              mhd_calcCharacteristics2d_sim, dscale, bclear, rvector, rcollection)
          
        case (NDIM3D)
          call afcsys_buildVectorTVD(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, NDIM3D,&
              mhd_calcFluxGalNoBdr3d_sim,&
              mhd_calcCharacteristics3d_sim, dscale, bclear, rvector, rcollection)
        end select
        
        !-----------------------------------------------------------------------

      case default
        call output_line('Invalid type of stabilisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcDivergenceVector')
        call sys_halt()
      end select
      
      !-------------------------------------------------------------------------
      ! Evaluate linear form for boundary integral (if any)
      !-------------------------------------------------------------------------
      
      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call mhd_calcLinfBdrCond1D(rproblemLevel, rboundaryCondition,&
            rsolution, dtime, -dscale, ssectionName, mhd_coeffVectorBdr1d_sim,&
            rvector, rcollection)
        
      case (NDIM2D)
        call mhd_calcLinfBdrCond2D(rproblemLevel, rboundaryCondition,&
            rsolution, dtime, -dscale, ssectionName, mhd_coeffVectorBdr2d_sim,&
            rvector, rcollection)
        
      case (NDIM3D)
!!$        call mhd_calcLinfBdrCond3D(rproblemLevel, rboundaryCondition,&
!!$            rsolution, dtime, -dscale, ssectionName, mhd_coeffVectorBdr3d_sim,&
!!$            rvector, rcollection)
        print *, "Boundary conditions in 3D have not been implemented yet!"
        stop
      end select
      
    end if

  end subroutine mhd_calcDivergenceVector

  !*****************************************************************************

!<subroutine>

  subroutine mhd_calcTimeDerivative(rproblemLevel, rtimestep,&
      rsolver, rsolution, ssectionName, rcollection, rvector,&
      rsource, rvector1, rvector2)

!<description>
    ! This subroutine calculates the approximate time derivative
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: source vector
    type(t_vectorBlock), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! collection structure
    type(t_collection), intent(inout) :: rcollection

    ! destination vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! OPTIONAL: auxiliary vectors used to compute the approximation to
    ! the time derivative (if not present, then temporal memory is allocated)
    type(t_vectorBlock), intent(inout), target, optional :: rvector1
    type(t_vectorBlock), intent(inout), target, optional :: rvector2
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_vectorBlock), pointer :: p_rvector1, p_rvector2
    real(DP) :: dnorm0, dnorm
    real(DP) :: depsAbsApproxTimeDerivative,depsRelApproxTimeDerivative
    integer :: inviscidAFC,viscousAFC
    integer :: iblock,iapproxtimederivativetype
    integer :: lumpedMassMatrix,consistentMassMatrix
    integer :: cafcstabTypeInviscid
    integer :: cafcstabTypeViscous
    integer :: ite,nmaxIterationsApproxTimeDerivative
    integer(I32) :: istabilisationSpecInviscid
    integer(I32) :: istabilisationSpecViscous
    logical :: bcompatible

    ! Set pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    
    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'inviscidAFC', inviscidAFC, 0)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'viscousAFC', viscousAFC, 0)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'iapproxtimederivativetype', iapproxtimederivativetype)

    ! Check if rvector is compatible to the solution vector;
    ! otherwise create new vector as a duplicate of the solution vector
    call lsysbl_isVectorCompatible(rvector, rsolution, bcompatible)
    if (.not.bcompatible)&
        call lsysbl_duplicateVector(rsolution, rvector,&
        LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

    !---------------------------------------------------------------------------
    
    ! How should we compute the approximate time derivative?
    select case(iapproxtimederivativetype)
      
    case(AFCSTAB_GALERKIN)
      
      ! Get more parameters from parameter list
      call parlst_getvalue_double(p_rparlist,&
          ssectionName, 'depsAbsApproxTimeDerivative',&
          depsAbsApproxTimeDerivative, 1e-4_DP)
      call parlst_getvalue_double(p_rparlist,&
          ssectionName, 'depsRelApproxTimeDerivative',&
          depsRelApproxTimeDerivative, 1e-2_DP)
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'nmaxIterationsApproxTimeDerivative',&
          nmaxIterationsApproxTimeDerivative, 5)

      ! Set up vector1 for computing the approximate time derivative
      if (present(rvector1)) then
        p_rvector1 => rvector1
      else
        allocate(p_rvector1)
      end if
      
      ! Check if rvector1 is compatible to the solution vector;
      ! otherwise create new vector as a duplicate of the solution vector
      call lsysbl_isVectorCompatible(p_rvector1, rsolution, bcompatible)
      if (.not.bcompatible)&
          call lsysbl_duplicateVector(rsolution, p_rvector1,&
          LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      
      ! Set up vector2 for computing the approximate time derivative
      if (present(rvector2)) then
        p_rvector2 => rvector2
      else
        allocate(p_rvector2)
      end if
      
      ! Check if rvector2 is compatible to the solution vector;
      ! otherwise create new vector as a duplicate of the solution vector
      call lsysbl_isVectorCompatible(p_rvector2, rsolution, bcompatible)
      if (.not.bcompatible)&
          call lsysbl_duplicateVector(rsolution, p_rvector2,&
          LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

      ! Make a backup copy of the stabilisation types because we
      ! have to overwrite them to enforce using the standard
      ! Galerkin scheme; this implies that their specification flags
      ! are changed, so make a backup copy of them, too
      if (inviscidAFC > 0) then
        cafcstabTypeInviscid&
            = rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType
        istabilisationSpecInviscid&
            = rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec
        rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType&
            = AFCSTAB_GALERKIN
      end if

      if (viscousAFC > 0) then
        cafcstabTypeViscous&
            = rproblemLevel%Rafcstab(viscousAFC)%cafcstabType
        istabilisationSpecViscous&
            = rproblemLevel%Rafcstab(viscousAFC)%istabilisationSpec
        rproblemLevel%Rafcstab(viscousAFC)%cafcstabType&
            = AFCSTAB_GALERKIN
      end if

      ! Compute $K(u^L)*u^L$ and store the result in rvector1
      call mhd_calcDivergenceVector(rproblemLevel,&
          rsolver%rboundaryCondition, rsolution, rtimestep%dTime,&
          1.0_DP, .true., p_rvector1, ssectionName, rcollection)
      
      ! Build the geometric source term (if any)
!!$      call mhd_calcGeometricSourceterm(p_rparlist, ssectionName,&
!!$          rproblemLevel, rsolution, 1.0_DP, .false., p_rvector1, rcollection)
      
      ! Apply the source vector to the residual (if any)
      if (present(rsource)) then
        if (rsource%NEQ .gt. 0)&
            call lsysbl_vectorLinearComb(rsource, p_rvector1, 1.0_DP, 1.0_DP)
      end if
      
      ! Reset stabilisation structures to their original configuration
      if (inviscidAFC > 0) then
        rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType&
            = cafcstabTypeInviscid
        rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec&
            = istabilisationSpecInviscid
      end if

      if (viscousAFC > 0) then
        rproblemLevel%Rafcstab(viscousAFC)%cafcstabType&
            = cafcstabTypeViscous
        rproblemLevel%Rafcstab(viscousAFC)%istabilisationSpec&
            = istabilisationSpecViscous
      end if

      ! Scale rvector1 by the inverse of the lumped mass matrix and store
      ! the result in rvector; this is the solution of the lumped version
      call lsysbl_invertedDiagMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
          p_rvector1, 1.0_DP, rvector)

      ! Store norm of the initial guess from the lumped version
      dnorm0 = lsysbl_vectorNorm(rvector1, LINALG_NORML2)
      
      richardson: do ite = 1, nmaxIterationsApproxTimeDerivative
        ! Initialise rvector2 by the constant right-hand side
        call lsysbl_copyVector(p_rvector1, p_rvector2)
        
        ! Compute the residual $rhs-M_C*u$ and store the result in rvector3
        do iblock = 1,rsolution%nblocks
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(consistentMassMatrix),&
              rvector%RvectorBlock(iblock), p_rvector2%RvectorBlock(iblock),&
              -1.0_DP, 1.0_DP)
        end do
          
        ! Scale rvector2 by the inverse of the lumped mass matrix
        call lsysbl_invertedDiagMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
            p_rvector2, 1.0_DP, p_rvector2)
        
        ! Apply solution increment (rvector2) to the previous solution iterate
        call lsysbl_vectorLinearComb(p_rvector2, rvector, 1.0_DP, 1.0_DP)
        
        ! Check for convergence
        dnorm = lsysbl_vectorNorm(p_rvector2, LINALG_NORML2)
        if ((dnorm .le. depsAbsApproxTimeDerivative) .or.&
            (dnorm .le. depsRelApproxTimeDerivative*dnorm0)) exit richardson
      end do richardson
      
      ! Release temporal memory
      if (.not.present(rvector1)) then
        call lsysbl_releaseVector(p_rvector1)
        deallocate(p_rvector1)
      end if
      if (.not.present(rvector2)) then
        call lsysbl_releaseVector(p_rvector2)
        deallocate(p_rvector2)
      end if
      
      !-----------------------------------------------------------------------
      
    case(AFCSTAB_UPWIND)

      ! Make a backup copy of the stabilisation types because we
      ! have to overwrite them to enforce using the standard
      ! Galerkin scheme; this implies that their specification flags
      ! are changed, so make a backup copy of them, too
      if (inviscidAFC > 0) then
        cafcstabTypeInviscid&
            = rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType
        istabilisationSpecInviscid&
            = rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec
        rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType&
            = AFCSTAB_UPWIND
      end if

      if (viscousAFC > 0) then
        cafcstabTypeViscous&
            = rproblemLevel%Rafcstab(viscousAFC)%cafcstabType
        istabilisationSpecViscous&
            = rproblemLevel%Rafcstab(viscousAFC)%istabilisationSpec
        rproblemLevel%Rafcstab(viscousAFC)%cafcstabType&
            = AFCSTAB_DMP
      end if

      ! Compute $L(u^L)*u^L$ and store the result in rvector
      call mhd_calcDivergenceVector(rproblemLevel,&
          rsolver%rboundaryCondition, rsolution, rtimestep%dTime,&
          1.0_DP, .true., rvector, ssectionName, rcollection)
      
      ! Build the geometric source term (if any)
!!$      call mhd_calcGeometricSourceterm(p_rparlist, ssectionName,&
!!$          rproblemLevel, rsolution, 1.0_DP, .false., rvector, rcollection)
      
      ! Apply the source vector to the residual (if any)
      if (present(rsource)) then
        if (rsource%NEQ .gt. 0)&
            call lsysbl_vectorLinearComb(rsource, rvector, 1.0_DP, 1.0_DP)
      end if

      ! Reset stabilisation structures to their original configuration
      if (inviscidAFC > 0) then
        rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType&
            = cafcstabTypeInviscid
        rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec&
            = istabilisationSpecInviscid
      end if

      if (viscousAFC > 0) then
        rproblemLevel%Rafcstab(viscousAFC)%cafcstabType&
            = cafcstabTypeViscous
        rproblemLevel%Rafcstab(viscousAFC)%istabilisationSpec&
            = istabilisationSpecViscous
      end if

      ! Scale it by the inverse of the lumped mass matrix
      call lsysbl_invertedDiagMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
          rvector, 1.0_DP, rvector)

    case default
      call output_line('Unsupported type of divergence term!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mhd_calcTimeDerivative')
      call sys_halt()
    end select

  end subroutine mhd_calcTimeDerivative

end module mhd_callback

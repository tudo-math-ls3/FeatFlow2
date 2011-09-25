!##############################################################################
!# ****************************************************************************
!# <name> hydro_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to solve the
!# compressible Euler/Navier Stokes equations in arbitrary spatial dimensions.
!#
!# The following callback functions are available:
!#
!# 1.) hydro_nlsolverCallback
!#     -> Callback routine for the nonlinear solver
!#
!# ****************************************************************************
!#
!# The following auxiliary routines are available:
!#
!# 1.) hydro_calcPrecondThetaScheme
!#     -> Calculates the nonlinear preconditioner
!#        used in the two-level theta-scheme
!#
!# 2.) hydro_calcJacobianThetaScheme
!#     -> Calculates the Jacobian matrix
!#        used in the two-level theta-scheme
!#
!# 3.) hydro_calcResidualThetaScheme
!#     -> Calculates the nonlinear residual vector
!#        used in the two-level theta-scheme
!#
!# 4.) hydro_calcRhsThetaScheme
!#     -> Calculates the explicit right-hand side vector
!#        used in the two-level theta-scheme
!#
!# 5.) hydro_calcRhsRungeKuttaScheme
!#     -> Calculates the right-hand side vector
!#        used in the explicit Runge-Kutta scheme
!#
!# 6.) hydro_setBoundaryCondition
!#     -> Imposes boundary conditions for nonlinear solver
!#        by filtering the system matrix and the solution/residual
!#        vector explicitly (i.e. strong boundary conditions)
!#
!# 7.) hydro_calcLinearisedFCT
!#     -> Calculates the linearised FCT correction
!#
!# 8.) hydro_calcFluxFCT
!#     -> Calculates the raw antidiffusive fluxes for FCT algorithm
!#
!# 9.) hydro_calcCorrectionFCT
!#     -> Calculates the contribution of the antidiffusive fluxes
!#        limited by the FCT algorithm and applies them to the residual
!#
!# 10.) hydro_limitEdgewiseVelocity
!#      -> Performs synchronised flux correction for the velocity
!#
!# 11.) hydro_limitEdgewiseMomentum
!#      -> Performs synchronised flux correction for the momentum
!#
!# 12.) hydro_coeffVectorFE
!#      -> Callback routine for the evaluation of linear forms
!#         using a given FE-solution for interpolation
!#
!# 13.) hydro_coeffVectorAnalytic
!#      -> Callback routine for the evaluation of linear forms
!#         using a given FE-solution for interpolation
!#
!# 14.) hydro_parseBoundaryCondition
!#      -> Callback routine for the treatment of boundary conditions
!#
!# 15.) hydro_calcBilfBdrCond = hydro_calcBilfBdrCond1D /
!#                              hydro_calcBilfBdrCond2D
!#      -> Calculates the bilinear form arising from the weak
!#        imposition of boundary conditions
!#
!# 16.) hydro_calcLinfBdrCond = hydro_calcLinfBdrCond1D /
!#                              hydro_calcLinfBdrCond2D
!#      -> Calculates the linear form arising from the weak
!#         imposition of boundary conditions
!#
!# 17.) hydro_calcGeometricSourceterm
!#      -> Calculates the geometric source term for axi-symmetric,
!#         cylindrical or sperical symmetric coordinate systems.
!#
!# Frequently asked questions?
!#
!# 1.) What is the magic behind subroutine 'hydro_nlsolverCallback'?
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

module hydro_callback

#include "hydro.h"
  use afcstabilisation
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
  use hydro_basic
  use hydro_callback1d
  use hydro_callback2d
  use hydro_callback3d
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
  public :: hydro_nlsolverCallback
  public :: hydro_calcPrecondThetaScheme
  public :: hydro_calcJacobianThetaScheme
  public :: hydro_calcResidualThetaScheme
  public :: hydro_calcRhsThetaScheme
  public :: hydro_calcRhsRungeKuttaScheme
  public :: hydro_setBoundaryCondition
  public :: hydro_calcLinearisedFCT
  public :: hydro_calcFluxFCT
  public :: hydro_calcCorrectionFCT
  public :: hydro_calcGeometricSourceTerm
  public :: hydro_limitEdgewiseVelocity
  public :: hydro_limitEdgewiseMomentum
  public :: hydro_coeffVectorFE
  public :: hydro_coeffVectorAnalytic
  public :: hydro_parseBoundaryCondition
  public :: hydro_calcBilfBdrCond1D
  public :: hydro_calcBilfBdrCond2D
  public :: hydro_calcLinfBdrCond1D
  public :: hydro_calcLinfBdrCond2D
  !*****************************************************************************

!<constants>

!<constantblock>
  ! Minimum number of equations for OpenMP parallelisation: If the number of
  ! equations is below this value, then no parallelisation is performed.
#ifndef HYDRO_GEOMSOURCE_NEQMIN_OMP
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: HYDRO_GEOMSOURCE_NEQMIN_OMP = 10000
#else
  integer, public            :: HYDRO_GEOMSOURCE_NEQMIN_OMP = 10000
#endif
#endif

!</constantblock>

!</constants>

contains

  !*****************************************************************************

!<subroutine>

  subroutine hydro_nlsolverCallback(rproblemLevel, rtimestep,&
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
      call hydro_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, ssectionName, rcollection)
    end if


    ! Do we have to calculate the constant right-hand side?
    ! --------------------------------------------------------------------------
    if ((iand(ioperationSpec, NLSOL_OPSPEC_CALCRHS)  .ne. 0)) then

      ! Compute the right-hand side
      call hydro_calcRhsRungeKuttaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rrhs, istep, ssectionName,&
          rcollection)
    end if


    ! Do we have to calculate the residual?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      if (istep .eq. 0) then
        ! Compute the constant right-hand side
        call hydro_calcRhsThetaScheme(rproblemLevel, rtimestep,&
            rsolver, rsolution0, rrhs, ssectionName, rcollection,&
            rsource)
      end if

      ! Compute the residual
      call hydro_calcResidualThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rrhs, rres, istep,&
          ssectionName, rcollection)
    end if


    ! Do we have to impose boundary conditions?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      ! Impose boundary conditions
      call hydro_setBoundaryCondition(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rres, ssectionName, rcollection)
    end if


    ! Set status flag
    istatus = 0

  end subroutine hydro_nlsolverCallback

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
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
    integer :: systemMatrix, lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ, inviscidAFC
    integer :: isystemCoupling, isystemPrecond, isystemFormat, imasstype, ivar

    ! Start time measurement for matrix evaluation
    p_rtimer => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyMatrix', ssectionName=ssectionName)
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'systemmatrix', systemMatrix)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'coeffMatrix_CX', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'coeffMatrix_CY', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'coeffMatrix_CZ', coeffMatrix_CZ)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'inviscidAFC', inviscidAFC)

    !---------------------------------------------------------------------------
    ! Check if fully explicit time-stepping is used
    !---------------------------------------------------------------------------
    if (rtimestep%theta .le. SYS_EPSREAL_DP) then

      call parlst_getvalue_int(p_rparlist, ssectionName,&
          'isystemformat', isystemFormat)
      call parlst_getvalue_int(p_rparlist, ssectionName,&
          'imasstype', imasstype)
      call parlst_getvalue_int(p_rparlist, ssectionName,&
          'lumpedmassmatrix', lumpedMassMatrix)
      call parlst_getvalue_int(p_rparlist, ssectionName,&
          'consistentmassmatrix', consistentMassMatrix)

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
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecondThetaScheme')
          call sys_halt()
        end select

      case (SYSTEM_BLOCKFORMAT)

        select case(imasstype)
        case (MASS_LUMPED)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, hydro_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(&
                rproblemLevel%Rmatrix(lumpedMassMatrix),&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                1.0_DP, 1.0_DP, .false., .false., .true., .true.)
          end do

        case (MASS_CONSISTENT)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, hydro_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(&
                rproblemLevel%Rmatrix(consistentMassMatrix),&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                1.0_DP, 1.0_DP, .false., .false., .true., .true.)
          end do

        case DEFAULT
          call output_line('Empty system matrix is invalid!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecondThetaScheme')
          call sys_halt()
        end select

      case DEFAULT
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecondThetaScheme')
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
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiagMatD1d_sim, hydro_calcMatGalMatD1d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiagMatD2d_sim, hydro_calcMatGalMatD2d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiagMatD3d_sim, hydro_calcMatGalMatD3d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiagMatD1d_sim, hydro_calcMatScDissMatD1d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiagMatD2d_sim, hydro_calcMatScDissMatD2d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiagMatD3d_sim, hydro_calcMatScDissMatD3d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecond')
          call sys_halt()
        end select


      case (DISSIPATION_ROE)

        ! Assemble divergence operator with Roe-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiagMatD1d_sim, hydro_calcMatRoeDissMatD1d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiagMatD2d_sim, hydro_calcMatRoeDissMatD2d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiagMatD3d_sim, hydro_calcMatRoeDissMatD3d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with the Rusanov-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiagMatD1d_sim, hydro_calcMatRusDissMatD1d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiagMatD2d_sim, hydro_calcMatRusDissMatD2d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiagMatD3d_sim, hydro_calcMatRusDissMatD3d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case DEFAULT
        ! Clear system matrix and apply (lumped) mass matrix only
        call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
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
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiag1d_sim, hydro_calcMatGal1d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiag2d_sim, hydro_calcMatGal2d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiag3d_sim, hydro_calcMatGal3d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiag1d_sim, hydro_calcMatScDiss1d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiag2d_sim, hydro_calcMatScDiss2d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiag3d_sim, hydro_calcMatScDiss3d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_ROE)

        ! Assemble divergence operator with Roe-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiag1d_sim, hydro_calcMatRoeDiss1d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiag2d_sim, hydro_calcMatRoeDiss2d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiag3d_sim, hydro_calcMatRoeDiss3d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with the Rusanov-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiag1d_sim, hydro_calcMatRusDiss1d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM2D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiag2d_sim, hydro_calcMatRusDiss2d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case (NDIM3D)
          call gfsys_buildDivOperator(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcMatDiag3d_sim, hydro_calcMatRusDiss3d_sim,&
              dscale, .true., rproblemLevel%RmatrixBlock(systemMatrix))

        case DEFAULT
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecondThetaScheme')
          call sys_halt()
        end select


      case DEFAULT
        ! Clear system matrix and apply (lumped) mass matrix only
        call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
      end select


    case DEFAULT
      call output_line('Invalid type of flow coupling!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecondThetaScheme')
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

      case DEFAULT

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

        do ivar = 1, hydro_getNVAR(rproblemLevel)
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

        do ivar = 1, hydro_getNVAR(rproblemLevel)
          call lsyssc_MatrixLinearComb(&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
              1.0_DP, 1.0_DP, .false., .false., .true., .true.)
        end do


      case DEFAULT

        !-----------------------------------------------------------------------
        ! Use the global operator for steady-state flow
        !
        !   $ A = -L $
        !
        ! Since we have assembled "-L" nothing needs to be done.
        !-----------------------------------------------------------------------

      end select


    case DEFAULT
      call output_line('Invalid system format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecondThetaScheme')
      call sys_halt()
    end select


    ! Impose boundary conditions in strong sence (if any)
    call bdrf_filterMatrix(rsolver%rboundaryCondition,&
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

  end subroutine hydro_calcPrecondThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcJacobianThetaScheme(rproblemLevel, rtimestep,&
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

    print *, "!!! The calculation of the Jacobian matrix for the compressible !!!"
    print *, "!!! Euler/Navier Stokes equations has yet not been implemented  !!!"
    stop

  end subroutine hydro_calcJacobianThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcRhsThetaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rrhs, ssectionName, rcollection, rsource)

!<description>
    ! This subroutine computes the constant right-hand side
    !
    !  $$ rhs = [M + (1-\theta)\Delta t K^n]U^n + S^n + b.c.`s  $$
    !
    ! where the (scaled) source term is optional.
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
    real(DP) :: dscale
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: consistentMassMatrix, lumpedMassMatrix, massMatrix
    integer :: inviscidAFC, imasstype, idissipationtype
    integer :: iblock


    ! Start time measurement for residual/rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyVector', ssectionName=ssectionName)
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'coeffmatrix_cx', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'coeffmatrix_cy', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'coeffmatrix_cz', coeffMatrix_CZ)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'imasstype', imasstype)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'idissipationtype', idissipationtype)

    ! Do we have some kind of mass matrix?
    select case(imasstype)
    case (MASS_LUMPED, MASS_CONSISTENT)

      ! Do we have an explicit part?
      if (rtimestep%theta .lt. 1.0_DP) then

        ! Compute scaling parameter
        dscale = (1.0_DP-rtimestep%theta) * rtimestep%dStep

        ! What type if stabilisation is applied?
        select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)

          !---------------------------------------------------------------------
          ! Compute the initial high-order right-hand side
          !
          !   $$ rhs = (1-theta)*dt*K(U^n)*U^n $$
          !---------------------------------------------------------------------

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildDivVector(&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                hydro_calcFluxGal1d_sim, dscale, .true., rrhs)

          case (NDIM2D)
            call gfsys_buildDivVector(&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                hydro_calcFluxGal2d_sim, dscale, .true., rrhs)

          case (NDIM3D)
            call gfsys_buildDivVector(&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxGal3d_sim, dscale, .true., rrhs)
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

            ! Assemble divergence of flux without dissipation

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxGal1d_sim, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxGal2d_sim, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxGal3d_sim, dscale, .true., rrhs)
            end select

          case (DISSIPATION_SCALAR)

            ! Assemble divergence of flux with scalar dissipation

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxScDiss1d_sim, dscale, .true. , rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxScDiss2d_sim, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxScDiss3d_sim, dscale, .true., rrhs)
            end select

          case (DISSIPATION_SCALAR_DSPLIT)

            ! Assemble divergence of flux with scalar dissipation
            ! adopting dimensional splitting

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxScDiss1d_sim, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxScDissDiSp2d_sim, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxScDissDiSp3d_sim, dscale, .true., rrhs)
            end select

          case (DISSIPATION_ROE)

            ! Assemble divergence of flux with Roe-type dissipation

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxRoeDiss1d_sim, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxRoeDiss2d_sim, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxRoeDiss3d_sim, dscale, .true., rrhs)
            end select

          case (DISSIPATION_ROE_DSPLIT)

            ! Assemble divergence of flux with Roe-type dissipation
            ! adopting dimensional splitting

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxRoeDiss1d_sim, dscale, .true., rrhs)
              
            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxRoeDissDiSp2d_sim, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxRoeDissDiSp3d_sim, dscale, .true., rrhs)
            end select

          case (DISSIPATION_RUSANOV)

            ! Assemble divergence of flux with Rusanov-type flux

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxRusDiss1d_sim, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxRusDiss2d_sim, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxRusDiss3d_sim, dscale, .true., rrhs)
            end select

          case (DISSIPATION_RUSANOV_DSPLIT)

            ! Assemble divergence of flux with Rusanov-type flux
            ! adopting dimensional splitting

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxRusDiss1d_sim, dscale, .true., rrhs)

            case (NDIM2D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxRusDissDiSp2d_sim, dscale, .true., rrhs)

            case (NDIM3D)
              call gfsys_buildDivVector(&
                  rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
                  hydro_calcFluxRusDissDiSp3d_sim, dscale, .true., rrhs)
            end select

          case DEFAULT
            call output_line('Invalid type of dissipation!',&
                OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcRhsThetaScheme')
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
                rproblemLevel%Rafcstab(inviscidAFC), rsolution, NDIM1D,&
                hydro_calcFluxGalNoBdr1d_sim,&
                hydro_calcCharacteristics1d_sim, dscale, .true., rrhs)

          case (NDIM2D)
            call gfsys_buildDivVectorTVD(&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution, NDIM2D,&
                hydro_calcFluxGalNoBdr2d_sim,&
                hydro_calcCharacteristics2d_sim, dscale, .true., rrhs)

          case (NDIM3D)
            call gfsys_buildDivVectorTVD(&
                rproblemLevel%Rafcstab(inviscidAFC), rsolution, NDIM3D,&
                hydro_calcFluxGalNoBdr3d_sim,&
                hydro_calcCharacteristics3d_sim, dscale, .true., rrhs)
          end select

        case DEFAULT
          call output_line('Invalid type of stabilisation!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcRhsThetaScheme')
          call sys_halt()
        end select

        !-----------------------------------------------------------------------
        ! Evaluate linear form for boundary integral (if any)
        !-----------------------------------------------------------------------

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call hydro_calcLinfBdrCond1D(rproblemLevel, rsolver,&
              rsolution, rtimestep%dTime-rtimestep%dStep, -dscale,&
              hydro_coeffVectorBdr1d_sim, rrhs, ssectionName, rcollection)       

        case (NDIM2D)
          call hydro_calcLinfBdrCond2D(rproblemLevel, rsolver,&
              rsolution, rtimestep%dTime-rtimestep%dStep, -dscale,&
              hydro_coeffVectorBdr2d_sim, rrhs, ssectionName, rcollection)

        case (NDIM3D)
          print *, "Not implemented yet"
          stop

        end select

        !-----------------------------------------------------------------------
        ! Compute the transient term
        !
        !   $$ rhs := M*U^n + rhs $$
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

        ! Build the geometric source term (if any)
        call hydro_calcGeometricSourceterm(p_rparlist, ssectionName,&
            rproblemLevel, rsolution, dscale, .false., rrhs, rcollection)

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
    if (present(rsource)) then
      if (rsource%NEQ .gt. 0)&
        call lsysbl_vectorLinearComb(rsource, rrhs, 1.0_DP, 1.0_DP)
    end if

    ! Stop time measurement for rhs evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine hydro_calcRhsThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcResidualThetaScheme(rproblemLevel,&
      rtimestep, rsolver, rsolution, rsolution0, rrhs, rres,&
      ite, ssectionName, rcollection, rsource)

!<description>
    ! This subroutine computes the nonlinear residual vector
    !
    ! $$ res^{(m)} = rhs - [M-\theta\Delta t K^{(m)}]U^{(m)} - S^{(m)} - b.c.`s $$
    !
    ! for the standard two-level theta-scheme, whereby the (scaled)
    ! source term $s^{(m)}$ is optional. The constant right-hand side
    !
    !  $$ rhs = [M + (1-\theta)\Delta t K^n]U^n + S^n + b.c.`s $$
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
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: consistentMassMatrix, lumpedMassMatrix, massMatrix
    integer :: inviscidAFC, imasstype, idissipationtype
    integer :: iblock

    ! Start time measurement for residual/rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyVector', ssectionName=ssectionName)
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'coeffmatrix_cx', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'coeffmatrix_cy', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'coeffmatrix_cz', coeffMatrix_CZ)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'imasstype', imasstype)

    ! Compute scaling parameter
    dscale = rtimestep%theta*rtimestep%dStep

    ! Do we have some kind of mass matrix?
    select case(imasstype)
    case (MASS_LUMPED, MASS_CONSISTENT)

      !-------------------------------------------------------------------------
      ! Initialize the residual for transient flows
      !
      !   $$ res = rhs-M*U^{(m)} - b.c.`s $$
      !-------------------------------------------------------------------------

      ! Apply constant right-hand side
      call lsysbl_copyVector(rrhs, rres)

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

      !-------------------------------------------------------------------------
      ! Evaluate linear form for boundary integral
      ! and geometric source term (if any)
      !-------------------------------------------------------------------------
      
      ! Do we have an implicit part?
      if (rtimestep%theta .gt. 0.0_DP) then
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call hydro_calcLinfBdrCond1D(rproblemLevel, rsolver,&
              rsolution, rtimestep%dTime, -dscale,&
              hydro_coeffVectorBdr1d_sim, rres, ssectionName, rcollection)
          
        case (NDIM2D)
          call hydro_calcLinfBdrCond2D(rproblemLevel, rsolver,&
              rsolution, rtimestep%dTime, -dscale,&
              hydro_coeffVectorBdr2d_sim, rres, ssectionName, rcollection)
          
        case (NDIM3D)
          print *, "Not implemented yet"
          stop
          
        end select

        ! Build the geometric source term (if any)
        call hydro_calcGeometricSourceterm(p_rparlist, ssectionName,&
            rproblemLevel, rsolution, dscale, .false., rres, rcollection)
      end if
    case DEFAULT

      !-----------------------------------------------------------------------
      ! Initialize the residual for stationary flows zeros
      !
      !   $$ res = rhs - b.c.`s $$
      !-----------------------------------------------------------------------

      ! Apply constant right-hand side
      call lsysbl_copyVector(rrhs, rres)

      ! Evaluate linear form for boundary integral (if any)
      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call hydro_calcLinfBdrCond1D(rproblemLevel,&
            rsolver, rsolution, rtimestep%dTime, -1.0_DP,&
            hydro_coeffVectorBdr1d_sim, rres, ssectionName, rcollection)

      case (NDIM2D)
        call hydro_calcLinfBdrCond2D(rproblemLevel,&
            rsolver, rsolution, rtimestep%dTime, -1.0_DP,&
            hydro_coeffVectorBdr2d_sim, rres, ssectionName, rcollection)
        
      case (NDIM3D)
        print *, "Not implemented yet!"
        stop
      end select

      ! Build the geometric source term (if any)
      call hydro_calcGeometricSourceterm(p_rparlist, ssectionName,&
          rproblemLevel, rsolution, 1.0_DP, .false., rres, rcollection)
      
    end select

    !---------------------------------------------------------------------------
    ! Update the residual vector
    !---------------------------------------------------------------------------

    ! What type if stabilisation is applied?
    select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)

    case (AFCSTAB_GALERKIN)

      !-------------------------------------------------------------------------
      ! Compute the high-order residual
      !
      !   $$ res := res + dt*theta*K(U^{(m)})*U^(m) $$
      !-------------------------------------------------------------------------

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            hydro_calcFluxGal1d_sim, dscale, .false., rres)

      case (NDIM2D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            hydro_calcFluxGal2d_sim, dscale, .false., rres)

      case (NDIM3D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            hydro_calcFluxGal3d_sim, dscale, .false., rres)
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
          ssectionName, 'idissipationtype', idissipationtype)

      ! What type of dissipation is applied?
      select case(idissipationtype)

      case (DISSIPATION_ZERO)

        ! Assemble divergence of flux without dissipation
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxGal1d_sim, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxGal2d_sim, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxGal3d_sim, dscale, .false., rres)
        end select

      case (DISSIPATION_SCALAR)

        ! Assemble divergence of flux with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxScDiss1d_sim, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxScDiss2d_sim, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxScDiss3d_sim, dscale, .false., rres)
        end select

      case (DISSIPATION_SCALAR_DSPLIT)

        ! Assemble divergence of flux with scalar dissipation
        ! adopting dimensional splitting

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxScDiss1d_sim, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxScDissDiSp2d_sim, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxScDissDiSp3d_sim, dscale, .false., rres)
        end select

      case (DISSIPATION_ROE)

        ! Assemble divergence of flux with Roe-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRoeDiss1d_sim, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRoeDiss2d_sim, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRoeDiss3d_sim, dscale, .false., rres)
        end select

      case (DISSIPATION_ROE_DSPLIT)

        ! Assemble divergence of flux with Roe-type dissipation
        ! adopting dimensional splitting

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRoeDiss1d_sim, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRoeDissDiSp2d_sim, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRoeDissDiSp3d_sim, dscale, .false., rres)
        end select

      case (DISSIPATION_RUSANOV)

        ! Assemble divergence of flux with Rusanov-type flux

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRusDiss1d_sim, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRusDiss2d_sim, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRusDiss3d_sim, dscale, .false., rres)
        end select

      case (DISSIPATION_RUSANOV_DSPLIT)

        ! Assemble divergence of flux with Rusanov-type flux adopting
        ! dimensional splitting

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRusDiss1d_sim, dscale, .false., rres)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRusDissDiSp2d_sim, dscale, .false., rres)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRusDissDiSp3d_sim, dscale, .false., rres)
        end select

      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcResidualThetaScheme')
        call sys_halt()
      end select

    case (AFCSTAB_FEMTVD)

      !-------------------------------------------------------------------------
      ! Compute the low-order residual + FEM-TVD stabilisation
      !
      !   $$ res = res + dt*theta*L(U^{(m)})*U^{(m)} + F(U^{(m)}) $$
      !-------------------------------------------------------------------------

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution, NDIM1D,&
            hydro_calcFluxGalNoBdr1d_sim,&
            hydro_calcCharacteristics1d_sim, dscale, .false., rres)

      case (NDIM2D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution, NDIM2D,&
            hydro_calcFluxGalNoBdr2d_sim,&
            hydro_calcCharacteristics2d_sim, dscale, .false., rres)

      case (NDIM3D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution, NDIM3D,&
            hydro_calcFluxGalNoBdr3d_sim,&
            hydro_calcCharacteristics3d_sim, dscale , .false., rres)
      end select

    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcResidualThetaScheme')
      call sys_halt()
    end select


    !-------------------------------------------------------------------------
    ! Perform algebraic flux correction for the inviscid term
    !
    !   $$ res = res + f^*(u^(m),u^n) $$
    !-------------------------------------------------------------------------

    ! What type if stabilisation is applied?
    select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)
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
      elseif (rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType&
              .eq. AFCSTAB_FEMFCT_ITERATIVE) then
        ! ... in each iteration for iterative limiting
        call lsysbl_invertedDiagMatVec(&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rrhs, 1.0_DP, p_rpredictor)
      end if

      ! Assemble the raw antidiffusive fluxes
      call hydro_calcFluxFCT(rproblemLevel, rsolution, rtimestep%theta,&
          rtimestep%dStep, 1.0_DP, (ite .eq. 0), rsolution, ssectionName,&
          rcollection)

      !-------------------------------------------------------------------------
      ! Set operation specifier
      !-------------------------------------------------------------------------

      if (ite .eq. 0) then
        ! Perform standard flux correction in zeroth iteration
        ioperationSpec = AFCSTAB_FCTALGO_STANDARD
      else
        select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)
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
      call hydro_calcCorrectionFCT(rproblemLevel, p_rpredictor,&
          rtimestep%dStep, .false., ioperationSpec, rres,&
          ssectionName, rcollection)

      ! Subtract corrected antidiffusion from right-hand side
      if (rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType&
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
      if (rsource%NEQ .gt. 0)&
      call lsysbl_vectorLinearComb(rsource, rres, -1.0_DP, 1.0_DP)
    end if

    ! Stop time measurement for residual/rhs evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine hydro_calcResidualThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcRhsRungeKuttaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rrhs, istep, ssectionName, rcollection)

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
    real(DP) :: dscale
    integer :: lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: imasstype, inviscidAFC, idissipationtype
    integer :: iblock


    print *, "WARNING: This subroutine has not been tested!"
    stop

    ! Start time measurement for residual/rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyVector', ssectionName=ssectionName)
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'coeffMatrix_CX', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'coeffMatrix_CY', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'coeffMatrix_CZ', coeffMatrix_CZ)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'inviscidAFC', inviscidAFC)

    !---------------------------------------------------------------------------
    ! Initialize the right-hand side vector
    !---------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'imasstype', imasstype)

    select case(imasstype)
    case (MASS_LUMPED)

      !-------------------------------------------------------------------------
      ! Initialize the right-hand side vector
      !
      !  $ M_L*U $
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)

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
          ssectionName, 'consistentmassmatrix',&
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
        ssectionName, 'inviscidAFC', inviscidAFC)

    !---------------------------------------------------------------------------
    ! Compute the scaling parameter
    !
    !   $ weight*(1-theta)*dt $
    !---------------------------------------------------------------------------

    dscale = rtimestep%DmultistepWeights(istep)*(1.0_DP-rtimestep%theta)*rtimestep%dStep

    select case(rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType)

    case (AFCSTAB_GALERKIN)

      !-------------------------------------------------------------------------
      ! Compute the high-order right-hand side
      !
      !   $$ rhs = M*U+weight*(1-theta)*dt*K(U)*U $$
      !-------------------------------------------------------------------------

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            hydro_calcFluxGal1d_sim, dscale, .false., rrhs)

      case (NDIM2D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            hydro_calcFluxGal2d_sim, dscale, .false., rrhs)

      case (NDIM3D)
        call gfsys_buildDivVector(&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
            hydro_calcFluxGal3d_sim, dscale, .false., rrhs)
      end select

    case (AFCSTAB_UPWIND,&
          AFCSTAB_FEMFCT_LINEARISED)

      !-----------------------------------------------------------------------
      ! Compute the low-order right-hand side
      !
      !   $$ rhs = M*U+weight*(1-theta)*dt*L(U)*U $$
      !-----------------------------------------------------------------------


      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'idissipationtype', idissipationtype)

      select case(idissipationtype)

      case (DISSIPATION_ZERO)

        ! Assemble divergence of flux without dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxGal1d_sim, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxGal2d_sim, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxGal3d_sim, dscale, .false., rrhs)
        end select

      case (DISSIPATION_SCALAR)

        ! Assemble divergence of flux with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxScDiss1d_sim, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxScDiss2d_sim, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxScDiss3d_sim, dscale, .false., rrhs)
        end select

      case (DISSIPATION_ROE)

        ! Assemble divergence of flux with Roe-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRoeDiss1d_sim, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRoeDiss2d_sim, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRoeDiss3d_sim, dscale, .false., rrhs)
        end select

      case (DISSIPATION_RUSANOV)

        ! Assemble divergence of flux with Rusanov-type flux

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRusDiss1d_sim, dscale, .false., rrhs)

        case (NDIM2D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRusDiss2d_sim, dscale, .false., rrhs)

        case (NDIM3D)
          call gfsys_buildDivVector(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              hydro_calcFluxRusDiss3d_sim, dscale, .false., rrhs)
        end select

      case DEFAULT
        call output_line('Invalid type of dissipation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcRhsRungeKuttaScheme')
        call sys_halt()
      end select

    case (AFCSTAB_FEMTVD)

      !-----------------------------------------------------------------------
      ! Compute the low-order right-hand side + FEM-TVD stabilisation
      !
      !   $$ rhs = M*U+weight*(1-theta)*dt*L(U)*U + F(U) $$
      !-----------------------------------------------------------------------

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution, NDIM1D,&
            hydro_calcFluxGalNoBdr1d_sim,&
            hydro_calcCharacteristics1d_sim, dscale, .false., rrhs)

      case (NDIM2D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution, NDIM2D,&
            hydro_calcFluxGalNoBdr2d_sim,&
            hydro_calcCharacteristics2d_sim, dscale, .false., rrhs)

      case (NDIM3D)
        call gfsys_buildDivVectorTVD(&
            rproblemLevel%Rafcstab(inviscidAFC), rsolution, NDIM3D,&
            hydro_calcFluxGalNoBdr3d_sim,&
            hydro_calcCharacteristics3d_sim, dscale, .false., rrhs)

      end select

    case DEFAULT
      call output_line('Invalid type of stabilisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcRhsRungeKuttaScheme')
      call sys_halt()
    end select

    ! Stop time measurement for global operator
    call stat_stopTimer(p_rtimer)

  end subroutine hydro_calcRhsRungeKuttaScheme

  !*****************************************************************************

!<subroutine>

  subroutine hydro_setBoundaryCondition(rproblemLevel, rtimestep,&
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

    case DEFAULT
      call output_line('Invalid nonlinear preconditioner!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_setBoundaryCondition')
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
          hydro_calcBoundaryvalues1d, istatus)

    case (NDIM2D)
      call bdrf_filterSolution(&
          rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixBlock(imatrix),&
          rsolution, rres, rsolution0, rtimestep%dTime,&
          hydro_calcBoundaryvalues2d, istatus)

    case (NDIM3D)
      call bdrf_filterSolution(&
          rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixBlock(imatrix),&
          rsolution, rres, rsolution0, rtimestep%dTime,&
          hydro_calcBoundaryvalues3d, istatus)

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_setBoundaryCondition')
      call sys_halt()
    end select

  end subroutine hydro_setBoundaryCondition
  
  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcLinearisedFCT(rbdrCond, rproblemLevel,&
      rtimestep, rsolver, rsolution, ssectionName, rcollection, rsource)

!<description>
    ! This subroutine calculates the linearised FCT correction
!</description>

!<input>
    ! boundary condition structure
    type(t_boundaryCondition), intent(in) :: rbdrCond

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
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_timestep) :: rtimestepAux
    type(t_vectorBlock), pointer :: p_rpredictor
    type(t_parlist), pointer :: p_rparlist
    character(len=SYS_STRLEN), dimension(:), pointer :: SfailsafeVariables
    integer :: inviscidAFC,lumpedMassMatrix,imassantidiffusiontype
    integer :: nfailsafe,ivariable,nvariable

    ! Set pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)

    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'inviscidAFC', inviscidAFC)

    ! Do we have to apply linearised FEM-FCT?
    if (inviscidAFC .le. 0) return
    if (rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType&
        .ne. AFCSTAB_FEMFCT_LINEARISED) return

    ! Get more parameters from parameter list
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'lumpedmassmatrix', lumpedmassmatrix)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'imassantidiffusiontype', imassantidiffusiontype)
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'nfailsafe', nfailsafe)
    
    !---------------------------------------------------------------------------
    ! Linearised FEM-FCT algorithm
    !---------------------------------------------------------------------------

    ! Should we apply consistent mass antidiffusion?
    if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
    ! Initialize dummy timestep
    rtimestepAux%dStep = 1.0_DP
    rtimestepAux%theta = 0.0_DP

    ! Set pointer to predictor
    p_rpredictor => rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorPredictor
    
    ! Compute low-order "right-hand side" without theta parameter
    call hydro_calcRhsThetaScheme(rproblemLevel, rtimestepAux,&
          rsolver, rsolution, p_rpredictor, ssectionName, rcollection,&
          rsource)

    ! Compute low-order predictor
    call lsysbl_invertedDiagMatVec(&
        rproblemLevel%Rmatrix(lumpedMassMatrix),&
        p_rpredictor, 1.0_DP, p_rpredictor)

      ! Build the raw antidiffusive fluxes with contribution from
      ! consistent mass matrix
    call hydro_calcFluxFCT(rproblemLevel, rsolution, 0.0_DP,&
          1.0_DP, 1.0_DP, .true., p_rpredictor, ssectionName, rcollection)
    else
      ! Build the raw antidiffusive fluxes without contribution from
      ! consistent mass matrix
      call hydro_calcFluxFCT(rproblemLevel, rsolution, 0.0_DP,&
          1.0_DP, 1.0_DP, .true., rsolution, ssectionName, rcollection)
    end if
    
    !---------------------------------------------------------------------------
    ! Perform failsafe flux correction (if required)
    !---------------------------------------------------------------------------

    if (nfailsafe .gt. 0) then

      ! Get number of failsafe variables
      nvariable = max(1,&
          parlst_querysubstrings(p_rparlist,&
          ssectionName, 'sfailsafevariable'))

      ! Allocate character array that stores all failsafe variable names
      allocate(SfailsafeVariables(nvariable))
      
      ! Initialize character array with failsafe variable names
      do ivariable = 1, nvariable
        call parlst_getvalue_string(p_rparlist,&
            ssectionName, 'sfailsafevariable',&
            Sfailsafevariables(ivariable), isubstring=ivariable)
      end do

      ! Compute FEM-FCT correction
      call hydro_calcCorrectionFCT(rproblemLevel,&
          rsolution, rtimestep%dStep, .false.,&
          AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_CORRECT,&
          rsolution, ssectionName, rcollection)
      
      ! Apply failsafe flux correction
      call afcstab_failsafeLimiting(&
          rproblemLevel%Rafcstab(inviscidAFC),&
          rproblemLevel%Rmatrix(lumpedMassMatrix),&
          SfailsafeVariables, rtimestep%dStep, nfailsafe,&
          hydro_getVariable, rsolution, p_rpredictor)

      ! Deallocate temporal memory
      deallocate(SfailsafeVariables)

    else
      
      ! Apply linearised FEM-FCT correction
      call hydro_calcCorrectionFCT(rproblemLevel,&
          rsolution, rtimestep%dStep, .false.,&
          AFCSTAB_FCTALGO_STANDARD+&
          AFCSTAB_FCTALGO_SCALEBYMASS,&
          rsolution, ssectionName, rcollection)
    end if

    ! Impose boundary conditions for the solution vector
    select case(rproblemLevel%rtriangulation%ndim)
    case (NDIM1D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dTime, hydro_calcBoundaryvalues1d)

    case (NDIM2D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dTime, hydro_calcBoundaryvalues2d)

    case (NDIM3D)
      call bdrf_filterVectorExplicit(rbdrCond, rsolution,&
          rtimestep%dTime, hydro_calcBoundaryvalues3d)
    end select

  end subroutine hydro_calcLinearisedFCT

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcFluxFCT(rproblemLevel, rsolution1, theta,&
      tstep, dscale, binit, rsolution2, ssectionName, rcollection)

!<description>
    ! This subroutine calculates the raw antidiffusive fluxes for
    ! the different FCT algorithms
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution1

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

    ! OPTIONAL: second solution vector
    type(t_vectorBlock), intent(in) :: rsolution2

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName
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
    integer :: idissipationtype, imassantidiffusiontype

    ! Get parameters from parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)

    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'inviscidAFC', inviscidAFC)
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
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTScDiss1d_sim,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rsolution2, rcollection)
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTScDiss1d_sim,&
              theta, tstep, dscale, binit, rcollection=rcollection)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTScDiss2d_sim,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rsolution2, rcollection)
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTScDiss2d_sim,&
              theta, tstep, dscale, binit, rcollection=rcollection)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTScDiss3d_sim,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rsolution2, rcollection)
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTScDiss3d_sim,&
              theta, tstep, dscale, binit, rcollection=rcollection)
        end if
      end select


    case (DISSIPATION_ROE)

      ! Assemble raw antidiffusive fluxes using Roe-type dissipation

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTRoeDiss1d_sim,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rsolution2, rcollection)
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTRoeDiss1d_sim,&
              theta, tstep, dscale, binit, rcollection=rcollection)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTRoeDiss2d_sim,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rsolution2, rcollection)
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTRoeDiss2d_sim,&
              theta, tstep, dscale, binit, rcollection=rcollection)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTRoeDiss3d_sim,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rsolution2, rcollection)
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTRoeDiss3d_sim,&
              theta, tstep, dscale, binit, rcollection=rcollection)
        end if
      end select


    case (DISSIPATION_RUSANOV)

      ! Assemble raw antidiffusive fluxes using th Rusanov-type dissipation

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTRusDiss1d_sim,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rsolution2, rcollection)
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTRusDiss1d_sim,&
              theta, tstep, dscale, binit, rcollection=rcollection)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTRusDiss2d_sim,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rsolution2, rcollection)
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTRusDiss2d_sim,&
              theta, tstep, dscale, binit, rcollection=rcollection)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTRusDiss3d_sim,&
              theta, tstep, dscale, binit,&
              rproblemLevel%Rmatrix(consistentMassMatrix),&
              rsolution2, rcollection)
        else
          call gfsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rsolution1, hydro_calcFluxFCTRusDiss3d_sim,&
              theta, tstep, dscale, binit, rcollection=rcollection)
        end if
      end select


    case DEFAULT
      call output_line('Invalid type of dissipation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcFluxFCT')
      call sys_halt()
    end select

  end subroutine hydro_calcFluxFCT

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcCorrectionFCT(rproblemLevel, rsolution,&
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
      nvartransformed = hydro_getNVARtransformed(rproblemLevel, slimitingvariable)

      ! What type of flux transformation is applied?
      if (trim(slimitingvariable) .eq. 'density') then

        ! Apply FEM-FCT algorithm for density fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDensity1d_sim, hydro_trafoDiffDensity1d_sim)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDensity2d_sim, hydro_trafoDiffDensity2d_sim)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDensity3d_sim, hydro_trafoDiffDensity3d_sim)
        end select

      elseif (trim(slimitingvariable) .eq. 'energy') then

        ! Apply FEM-FCT algorithm for energy fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxEnergy1d_sim, hydro_trafoDiffEnergy1d_sim)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxEnergy2d_sim, hydro_trafoDiffEnergy2d_sim)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxEnergy3d_sim, hydro_trafoDiffEnergy3d_sim)
        end select

      elseif (trim(slimitingvariable) .eq. 'pressure') then

        ! Apply FEM-FCT algorithm for pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxPressure1d_sim, hydro_trafoDiffPressure1d_sim)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxPressure2d_sim, hydro_trafoDiffPressure2d_sim)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxPressure3d_sim, hydro_trafoDiffPressure3d_sim)
        end select

      elseif (trim(slimitingvariable) .eq. 'velocity') then

        ! Apply FEM-FCT algorithm for velocity fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxVelocity1d_sim, hydro_trafoDiffVelocity1d_sim)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxVelocity2d_sim, hydro_trafoDiffVelocity2d_sim,&
              fcb_limitEdgewise=hydro_limitEdgewiseVelocity)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxVelocity3d_sim, hydro_trafoDiffVelocity3d_sim,&
              fcb_limitEdgewise=hydro_limitEdgewiseVelocity)
        end select

      elseif (trim(slimitingvariable) .eq. 'momentum') then

        ! Apply FEM-FCT algorithm for momentum fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxMomentum1d_sim, hydro_trafoDiffMomentum1d_sim)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxMomentum2d_sim, hydro_trafoDiffMomentum2d_sim,&
              fcb_limitEdgewise=hydro_limitEdgewiseMomentum)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxMomentum3d_sim, hydro_trafoDiffMomentum3d_sim,&
              fcb_limitEdgewise=hydro_limitEdgewiseMomentum)
        end select
        
      elseif (trim(slimitingvariable) .eq. 'density,energy') then

        ! Apply FEM-FCT algorithm for density and energy fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenEng1d_sim, hydro_trafoDiffDenEng1d_sim)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenEng2d_sim, hydro_trafoDiffDenEng2d_sim)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenEng3d_sim, hydro_trafoDiffDenEng3d_sim)
        end select

      elseif (trim(slimitingvariable) .eq. 'density,pressure') then

        ! Apply FEM-FCT algorithm for density and pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenPre1d_sim, hydro_trafoDiffDenPre1d_sim)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenPre2d_sim, hydro_trafoDiffDenPre2d_sim)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenPre3d_sim, hydro_trafoDiffDenPre3d_sim)
        end select

      elseif (trim(slimitingvariable) .eq. 'density,energy,momentum') then

        nvartransformed = hydro_getNVARtransformed(rproblemLevel, slimitingvariable)

        ! Apply FEM-FCT algorithm for full conservative fluxes
        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
            rsolution, dscale, bclear, iopSpec, rresidual)

      elseif (trim(slimitingvariable) .eq. 'density,pressure,velocity') then

        nvartransformed = hydro_getNVARtransformed(rproblemLevel, slimitingvariable)

        ! Apply FEM-FCT algorithm for density, velocity and pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenPreVel1d_sim, hydro_trafoDiffDenPreVel1d_sim)
        case (NDIM2D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenPreVel2d_sim, hydro_trafoDiffDenPreVel2d_sim)
        case (NDIM3D)
          call gfsys_buildDivVectorFCT(&
              rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenPreVel3d_sim, hydro_trafoDiffDenPreVel3d_sim)
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

        call gfsys_buildDivVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrix), p_rafcstab,&
            rsolution, dscale, bclear, iopSpec, rresidual)
        
        ! Nothing more needs to be done
        return
        else
          cycle
        end if

      else
        call output_line('Invalid type of flux transformation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcCorrectionFCT')
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

  end subroutine hydro_calcCorrectionFCT

  ! ***************************************************************************

!<subroutine>

  subroutine hydro_limitEdgewiseVelocity(IverticesAtEdge, NEDGE, NEQ,&
      NVAR, NVARtransformed, ndim1, ndim2, Dx, Dflux, Drp, Drm, Dalpha,&
      fcb_calcFluxTransformation_sim, Dflux0, rcollection)

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

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! OPTIONAL: callback function to compute variable transformation
    include '../../../../../kernel/PDEOperators/intf_gfsyscallback.inc'
    optional :: fcb_calcFluxTransformation_sim

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

    ! auxiliary arras
    real(DP), dimension(:,:,:), pointer :: DdataAtEdge
    real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge

    ! local variables
    real(DP), dimension(NVARtransformed) :: R_ij,R_ji,Uij
    real(DP) :: alpha_ij
    integer :: idx,IEDGEset,IEDGEmax,i,j,iedge
    
    ! Do we have to use the explicit fluxes as constraints?
    if (present(Dflux0)) then
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
            DdataAtEdge(:,1,idx) = Dx(:,IverticesAtEdge(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IverticesAtEdge(2,iedge))
          end do

          ! Use callback function to compute transformed fluxes
          call fcb_calcFluxTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              Dflux(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
              DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1))

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IverticesAtEdge(1,iedge)
            j = IverticesAtEdge(2,iedge)

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
            DdataAtEdge(:,1,idx) = Dx(IverticesAtEdge(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IverticesAtEdge(2,iedge),:)
          end do

          ! Use callback function to compute transformed fluxes
          call fcb_calcFluxTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              Dflux(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
              DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1))
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IverticesAtEdge(1,iedge)
            j = IverticesAtEdge(2,iedge)

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
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_limitEdgewiseVelocity')
        call sys_halt()
        
      end if
    end if
    
  end subroutine hydro_limitEdgewiseVelocity

  ! ***************************************************************************

!<subroutine>

  subroutine hydro_limitEdgewiseMomentum(IverticesAtEdge, NEDGE, NEQ,&
      NVAR, NVARtransformed, ndim1, ndim2, Dx, Dflux, Drp, Drm, Dalpha,&
      fcb_calcFluxTransformation_sim, Dflux0, rcollection)

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

    ! Edge data structure
    integer, dimension(:,:), intent(in) :: IverticesAtEdge

    ! OPTIONAL: callback function to compute variable transformation
    include '../../../../../kernel/PDEOperators/intf_gfsyscallback.inc'
    optional :: fcb_calcFluxTransformation_sim

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

    ! auxiliary arras
    real(DP), dimension(:,:,:), pointer :: DdataAtEdge
    real(DP), dimension(:,:,:), pointer :: DtransformedFluxesAtEdge
    
    ! local variables
    real(DP), dimension(NVARtransformed) :: R_ij,R_ji,Uij
    real(DP) :: alpha_ij
    integer :: idx,IEDGEset,IEDGEmax,i,j,iedge
    
    ! Do we have to use the explicit fluxes as constraints?
    if (present(Dflux0)) then
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
            DdataAtEdge(:,1,idx) = Dx(:,IverticesAtEdge(1,iedge))
            DdataAtEdge(:,2,idx) = Dx(:,IverticesAtEdge(2,iedge))
          end do

          ! Use callback function to compute transformed fluxes
          call fcb_calcFluxTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              Dflux(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
              DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1))

          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IverticesAtEdge(1,iedge)
            j = IverticesAtEdge(2,iedge)

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
            DdataAtEdge(:,1,idx) = Dx(IverticesAtEdge(1,iedge),:)
            DdataAtEdge(:,2,idx) = Dx(IverticesAtEdge(2,iedge),:)
          end do

          ! Use callback function to compute transformed fluxes
          call fcb_calcFluxTransformation_sim(&
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
              Dflux(:,IEDGEset:IEDGEmax), IEDGEmax-IEDGEset+1,&
              DtransformedFluxesAtEdge(:,:,1:IEDGEmax-IEDGEset+1))
          
          ! Loop through all edges in the current set
          ! and scatter the entries to the global vector
          do idx = 1, IEDGEmax-IEDGEset+1
            
            ! Get actual edge number
            iedge = idx+IEDGEset-1
            
            ! Get position of nodes
            i = IverticesAtEdge(1,iedge)
            j = IverticesAtEdge(2,iedge)

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
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_limitEdgewiseMomentum')
        call sys_halt()
        
      end if
    end if

  end subroutine hydro_limitEdgewiseMomentum

  ! ***************************************************************************

!<subroutine>

  subroutine hydro_coeffVectorFE(rdiscretisation, rform,&
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
      call output_line ('Invalid system format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_coeffVectorFE')
      call sys_halt()
    end select
    
  end subroutine hydro_coeffVectorFE

  ! ***************************************************************************

!<subroutine>

  subroutine hydro_coeffVectorAnalytic(rdiscretisation, rform,&
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
    
  end subroutine hydro_coeffVectorAnalytic

  ! *****************************************************************************

!<subroutine>

  subroutine hydro_parseBoundaryCondition(cbdrCondType, ndimension,&
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

    case ('FREESLIP_STRONG')
      ibdrCondType = BDRC_FREESLIP + BDRC_STRONG

    case ('FREESLIP_WEAK')
      ibdrCondType = BDRC_FREESLIP + BDRC_WEAK

    case ('RLXFREESLIP_STRONG')
      ibdrCondType = BDRC_RLXFREESLIP + BDRC_STRONG

    case ('RLXFREESLIP_WEAK')
      ibdrCondType = BDRC_RLXFREESLIP + BDRC_WEAK

    case ('VISCOUSWALL_STRONG')
      ibdrCondType = BDRC_VISCOUSWALL + BDRC_STRONG

    case ('VISCOUSWALL_WEAK')
      ibdrCondType = BDRC_VISCOUSWALL + BDRC_WEAK

    case ('SUPEROUTLET_STRONG')
      ibdrCondType = BDRC_SUPEROUTLET
      ! No strong boundary conditions are prescribed

    case ('SUPEROUTLET_WEAK')
      ibdrCondType = BDRC_SUPEROUTLET + BDRC_WEAK

    case ('SUBOUTLET_STRONG')
      ibdrCondType = BDRC_SUBOUTLET + BDRC_STRONG

    case ('SUBOUTLET_WEAK')
      ibdrCondType = BDRC_SUBOUTLET + BDRC_WEAK

    case ('MASSOUTLET_STRONG')
      ibdrCondType = BDRC_MASSOUTLET + BDRC_STRONG

    case ('MASSOUTLET_WEAK')
      ibdrCondType = BDRC_MASSOUTLET + BDRC_WEAK
      
    case ('FREESTREAM_STRONG')
      ibdrCondType = BDRC_FREESTREAM + BDRC_STRONG
      
    case ('FREESTREAM_WEAK')
      ibdrCondType = BDRC_FREESTREAM + BDRC_WEAK

    case ('SUPERINLET_STRONG')
      ibdrCondType = BDRC_SUPERINLET + BDRC_STRONG

    case ('SUPERINLET_WEAK')
      ibdrCondType = BDRC_SUPERINLET + BDRC_WEAK

    case ('SUBINLET_STRONG')
      ibdrCondType = BDRC_SUBINLET + BDRC_STRONG

    case ('SUBINLET_WEAK')
      ibdrCondType = BDRC_SUBINLET + BDRC_WEAK

    case ('MASSINLET_STRONG')
      ibdrCondType = BDRC_MASSINLET + BDRC_STRONG

    case ('MASSINLET_WEAK')
      ibdrCondType = BDRC_MASSINLET + BDRC_WEAK

    case ('PERIODIC_STRONG')
      ibdrCondType = BDRC_PERIODIC + BDRC_STRONG
      
    case ('PERIODIC_WEAK')
      ibdrCondType = BDRC_PERIODIC + BDRC_WEAK
      
    case ('ANTIPERIODIC_STRONG')
      ibdrCondType = BDRC_ANTIPERIODIC + BDRC_STRONG
      
    case ('ANTIPERIODIC_WEAK')
      ibdrCondType = BDRC_ANTIPERIODIC + BDRC_WEAK

    case default
      read(cbdrCondType, '(I3)') ibdrCondType
    end select

    
    ! Determine number of mathematical expressions
    select case (iand(ibdrCondType, BDRC_TYPEMASK))
      
    case (BDRC_FREESLIP, BDRC_VISCOUSWALL, BDRC_SUPEROUTLET)
      nexpressions = 0

    case (BDRC_SUBOUTLET, BDRC_MASSOUTLET, BDRC_RLXFREESLIP)
      nexpressions = 1
     
    case (BDRC_MASSINLET)
      nexpressions = 2

    case (BDRC_SUBINLET)
      nexpressions = 3

    case (BDRC_FREESTREAM, BDRC_SUPERINLET)
      nexpressions = ndimension+2
      
    case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      nexpressions = -1

    case default
      nexpressions = 0
    end select

  end subroutine hydro_parseBoundaryCondition

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcBilfBdrCond1D(rproblemLevel, rsolver,&
      rsolution, dtime, dscale, fcoeff_buildMatrixScBdr1D_sim,&
      rmatrix, ssectionName, rcollection)

!<description>
    ! This subroutine computes the bilinear form arising from the weak
    ! imposition of boundary conditions in 1D.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(in), target :: rsolution

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! callback routine for nonconstant coefficient matrices.
    include '../../../../../kernel/DOFMaintenance/intf_coefficientMatrixScBdr1D.inc'
!</intput>

!<inputoutput>
    ! matrix
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! At the moment, nothing is done in this subroutine and it should
    ! not be called. It may be necessary to assemble some bilinear
    ! forms at the boundary in future.

  end subroutine hydro_calcBilfBdrCond1D

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcBilfBdrCond2D(rproblemLevel, rsolver,&
      rsolution, dtime, dscale, fcoeff_buildMatrixScBdr2D_sim,&
      rmatrix, ssectionName, rcollection, cconstrType)

!<description>
    ! This subroutine computes the bilinear form arising from the weak
    ! imposition of boundary conditions in 2D.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(in), target :: rsolution

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! callback routine for nonconstant coefficient matrices.
    include '../../../../../kernel/DOFMaintenance/intf_coefficientMatrixScBdr2D.inc'

    ! OPTIONAL: One of the BILF_MATC_xxxx constants that allow to
    ! specify the matrix construction method. If not specified,
    ! BILF_MATC_ELEMENTBASED is used.
    integer, intent(in), optional :: cconstrType
!</intput>

!<inputoutput>
    ! matrix
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! At the moment, nothing is done in this subroutine and it should
    ! not be called. It may be necessary to assemble some bilinear
    ! forms at the boundary in future.

  end subroutine hydro_calcBilfBdrCond2D

!*****************************************************************************

!<subroutine>

  subroutine hydro_calcLinfBdrCond1D(rproblemLevel, rsolver,&
      rsolution, dtime, dscale, fcoeff_buildVectorBlBdr1D_sim,&
      rvector, ssectionName, rcollection)

!<description>
    ! This subroutine computes the linear form arising from the weak
    ! imposition of boundary conditions in 1D.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(in), target :: rsolution

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! callback routine for nonconstant coefficient vectors.
    include '../../../../../kernel/DOFMaintenance/intf_coefficientVectorBlBdr1D.inc'
!</intput>

!<inputoutput>
    ! residual/right-hand side vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_boundaryCondition), pointer :: p_rboundaryCondition
    type(t_collection) :: rcollectionTmp
    type(t_linearForm) :: rform
    integer, dimension(:), pointer :: p_IbdrCondType
    integer :: ibct

    ! Evaluate linear form for boundary integral and return if
    ! there are no weak boundary conditions available
    p_rboundaryCondition => rsolver%rboundaryCondition
    if (.not.p_rboundaryCondition%bWeakBdrCond) return

    ! Check if we are in 1D
    if (rproblemLevel%rtriangulation%ndim .ne. NDIM1D) then
      call output_line('Spatial dimension must be 1D!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcLinfBdrCond1D')
      call sys_halt()
    end if

    ! Initialize temporal collection structure
    call collct_init(rcollectionTmp)

    ! Prepare quick access arrays of temporal collection structure
    rcollectionTmp%SquickAccess(1) = ''
    rcollectionTmp%SquickAccess(2) = 'rfparser'
    rcollectionTmp%DquickAccess(1) = dtime
    rcollectionTmp%DquickAccess(2) = dscale

    ! Attach user-defined collection structure to temporal collection
    ! structure (may be required by the callback function)
    rcollectionTmp%p_rnextCollection => rcollection

    ! Attach solution vector to temporal collection structure
    rcollectionTmp%p_rvectorQuickAccess1 => rsolution

    ! Attach function parser from boundary conditions to collection
    ! structure and specify its name in quick access string array
    call collct_setvalue_pars(rcollectionTmp, 'rfparser',&
        p_rboundaryCondition%rfparser, .true.)

    
    ! Set pointers
    call storage_getbase_int(p_rboundaryCondition%h_IbdrCondType, p_IbdrCondType)
    
    ! Loop over all boundary components
    do ibct = 1, p_rboundaryCondition%iboundarycount
      
      ! Check if this component has weak boundary conditions
      if (iand(p_IbdrCondType(ibct), BDRC_WEAK) .ne. BDRC_WEAK) cycle

      ! Prepare further quick access arrays of temporal collection
      ! structure with boundary component, type and maximum expressions
      rcollectionTmp%IquickAccess(1) = p_IbdrCondType(ibct)
      rcollectionTmp%IquickAccess(2) = ibct
      rcollectionTmp%IquickAccess(3) = p_rboundaryCondition%nmaxExpressions
      
      ! Initialize the linear form
      rform%itermCount = 1
      rform%Idescriptors(1) = DER_FUNC
      
      ! Assemble the linear form
      if (rvector%nblocks .eq. 1) then
        call linf_buildVecIntlScalarBdr1d(rform, .false.,&
            rvector%RvectorBlock(1), fcoeff_buildVectorBlBdr1D_sim,&
            ibct, rcollectionTmp)
      else
        call linf_buildVectorBlockBdr1d(rform, .false.,&
            rvector, fcoeff_buildVectorBlBdr1D_sim,&
            ibct, rcollectionTmp)
      end if
      
    end do ! ibct
    
    ! Release temporal collection structure
    call collct_done(rcollectionTmp)

  end subroutine hydro_calcLinfBdrCond1D

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcLinfBdrCond2D(rproblemLevel, rsolver,&
      rsolution, dtime, dscale, fcoeff_buildVectorBlBdr2D_sim,&
      rvector, ssectionName, rcollection)

!<description>
    ! This subroutine computes the linear form arising from the weak
    ! imposition of boundary conditions in 2D.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(in), target :: rsolution

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! callback routine for nonconstant coefficient vectors.
    include '../../../../../kernel/DOFMaintenance/intf_coefficientVectorBlBdr2D.inc'
!</intput>

!<inputoutput>
    ! residual/right-hand side vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_boundaryCondition), pointer :: p_rboundaryCondition
    type(t_collection) :: rcollectionTmp
    type(t_boundaryRegion) :: rboundaryRegion,rboundaryRegionMirror,rregion
    type(t_linearForm) :: rform
    integer, dimension(:), pointer :: p_IbdrCondCpIdx, p_IbdrCondType
    integer, dimension(:), pointer :: p_IbdrCompPeriodic, p_IbdrCondPeriodic
    integer :: ibct, isegment, ccubTypeBdr

    ! Evaluate linear form for boundary integral and return if
    ! there are no weak boundary conditions available
    p_rboundaryCondition => rsolver%rboundaryCondition
    if (.not.p_rboundaryCondition%bWeakBdrCond) return

    ! Check if we are in 2D
    if (rproblemLevel%rtriangulation%ndim .ne. NDIM2D) then
      call output_line('Spatial dimension must be 2D!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcLinfBdrCond2D')
      call sys_halt()
    end if

    ! Get pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)

    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'ccubTypeBdr', ccubTypeBdr)

    ! Initialize temporal collection structure
    call collct_init(rcollectionTmp)

    ! Prepare quick access arrays of temporal collection structure
    rcollectionTmp%SquickAccess(1) = ''
    rcollectionTmp%SquickAccess(2) = 'rfparser'
    rcollectionTmp%DquickAccess(1) = dtime
    rcollectionTmp%DquickAccess(2) = dscale
    rcollectionTmp%IquickAccess(4) = ccubTypeBdr

    ! Attach user-defined collection structure to temporal collection
    ! structure (may be required by the callback function)
    rcollectionTmp%p_rnextCollection => rcollection

    ! Attach solution vector to first quick access vector of the
    ! temporal collection structure
    rcollectionTmp%p_rvectorQuickAccess1 => rsolution

    ! Attach function parser from boundary conditions to collection
    ! structure and specify its name in quick access string array
    call collct_setvalue_pars(rcollectionTmp, 'rfparser',&
        p_rboundaryCondition%rfparser, .true.)

    
    ! Set pointers
    call storage_getbase_int(p_rboundaryCondition%h_IbdrCondCpIdx,&
        p_IbdrCondCpIdx)
    call storage_getbase_int(p_rboundaryCondition%h_IbdrCondType,&
        p_IbdrCondType)

    ! Set additional pointers for periodic boundary conditions
    if (p_rboundaryCondition%bPeriodic) then
      call storage_getbase_int(p_rboundaryCondition%h_IbdrCompPeriodic,&
          p_IbdrCompPeriodic)
      call storage_getbase_int(p_rboundaryCondition%h_IbdrCondPeriodic,&
          p_IbdrCondPeriodic)
    end if
    
    ! Loop over all boundary components
    do ibct = 1, p_rboundaryCondition%iboundarycount
      
      ! Loop over all boundary segments
      do isegment = p_IbdrCondCpIdx(ibct), p_IbdrCondCpIdx(ibct+1)-1
        
        ! Check if this segment has weak boundary conditions
        if (iand(p_IbdrCondType(isegment), BDRC_WEAK) .ne. BDRC_WEAK) cycle
        
        ! Prepare further quick access arrays of temporal collection
        ! structure with boundary component, type and maximum expressions
        rcollectionTmp%IquickAccess(1) = p_IbdrCondType(isegment)
        rcollectionTmp%IquickAccess(2) = isegment
        rcollectionTmp%IquickAccess(3) = p_rboundaryCondition%nmaxExpressions
        
        ! Initialize the linear form
        rform%itermCount = 1
        rform%Idescriptors(1) = DER_FUNC
        
        ! Create boundary segment in 01-parametrisation
        call bdrc_createRegion(p_rboundaryCondition, ibct,&
            isegment-p_IbdrCondCpIdx(ibct)+1, rboundaryRegion)
        
        ! Check if special treatment of mirror boundary condition is required
        if ((iand(p_IbdrCondType(isegment), BDRC_TYPEMASK) .eq. BDRC_PERIODIC) .or.&
            (iand(p_IbdrCondType(isegment), BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC)) then
          
          ! Create boundary region for mirror boundary in 01-parametrisation
          call bdrc_createRegion(p_rboundaryCondition, p_IbdrCompPeriodic(isegment),&
              p_IbdrCondPeriodic(isegment)-p_IbdrCondCpIdx(p_IbdrCompPeriodic(isegment))+1,&
              rboundaryRegionMirror)
          
          ! Attach boundary regin to temporal collection structure
          call collct_setvalue_bdreg(rcollectionTmp, 'rboundaryRegionMirror',&
              rboundaryRegionMirror, .true.)

          ! In the callback-function, the minimum/maximum parameter
          ! values of the boundary region and its mirrored
          ! counterpartqq are required in length parametrisation to
          ! determine the parameter values of the mirrored cubature
          ! points. Therefore, we make a copy of both boundary
          ! regions, convert them to length parametrisation and attach
          ! the minimum/maximum parameter values to the quick access
          ! arrays of the temporal collection structure.
          rregion = rboundaryRegion
          call boundary_convertRegion(&
              rvector%RvectorBlock(1)%p_rspatialDiscr%p_rboundary,&
              rregion, BDR_PAR_LENGTH)
          
          ! Prepare quick access array of temporal collection structure
          rcollectionTmp%DquickAccess(3) = rregion%dminParam
          rcollectionTmp%DquickAccess(4) = rregion%dmaxParam

          rregion = rboundaryRegionMirror
          call boundary_convertRegion(&
              rvector%RvectorBlock(1)%p_rspatialDiscr%p_rboundary,&
              rregion, BDR_PAR_LENGTH)

          ! Prepare quick access array of temporal collection structure
          rcollectionTmp%DquickAccess(5) = rregion%dminParam
          rcollectionTmp%DquickAccess(6) = rregion%dmaxParam          
        end if
        ! Assemble the linear form
        if (rvector%nblocks .eq. 1) then
          call linf_buildVecIntlScalarBdr2d(rform, ccubTypeBdr, .false.,&
              rvector%RvectorBlock(1), fcoeff_buildVectorBlBdr2D_sim,&
              rboundaryRegion, rcollectionTmp)
        else
          call linf_buildVectorBlockBdr2d(rform, ccubTypeBdr, .false.,&
              rvector, fcoeff_buildVectorBlBdr2D_sim,&
              rboundaryRegion, rcollectionTmp)
        end if
        
      end do ! isegment
    end do ! ibct
    
    ! Release temporal collection structure
    call collct_done(rcollectionTmp)
    
  end subroutine hydro_calcLinfBdrCond2D

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcGeometricSourceterm(rparlist, ssectionName,&
      rproblemLevel, rsolution, dscale, bclear, rsource, rcollection)

!<description>
    ! This subroutine calculates the geometric source term for
    ! axi-symmetric, cylindrically symmetric or sperically symmetric
    ! coordinate system.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! Clear source vector?
    logical, intent(in) :: bclear
!</input>

!<inputoutput>
    ! source vector to be assembled
    type(t_vectorBlock), intent(inout) :: rsource

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:,:), pointer :: p_Dcoords
    real(DP), dimension(:), pointer :: p_DdataSolution
    real(DP), dimension(:), pointer :: p_DdataSource
    real(DP), dimension(:), pointer :: p_DdataMassMatrix
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    character(LEN=SYS_STRLEN) :: smode
    real(DP) :: deffectiveRadius
    integer :: neq, nvar, icoordsystem, igeometricsourcetype
    integer :: isystemFormat, massmatrix

    ! Check of source and solution vector are compatible
    call lsysbl_isVectorCompatible(rsolution, rsource)


    ! Get parameters from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'mode', smode)
    
    ! Are we in primal or dual mode?
    if (trim(smode) .eq. 'primal') then
    
      !-------------------------------------------------------------------------
      ! We are in primal mode
      !-------------------------------------------------------------------------
      
      ! Get further parameters from parameter list
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'isystemformat', isystemformat)
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'icoordsystem', icoordsystem, COORDS_CARTESIAN)
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'igeometricsourcetype', igeometricsourcetype, MASS_LUMPED)
      call parlst_getvalue_double(rparlist,&
          ssectionName, 'deffectiveRadius', deffectiveRadius, 1e-4_DP)
      
      ! Get pointers
      call lsysbl_getbase_double(rsolution, p_DdataSolution)
      call lsysbl_getbase_double(rsource, p_DdataSource)

      ! Get number of equations and variables
      nvar = hydro_getNVAR(rproblemLevel)
      neq  = rsolution%NEQ/nvar
      
      ! Get coordinates of the triangulation
      ! NOTE: This implementation only works for linear and bilinear finite
      !       elements where the nodal degrees of freedom are located at the
      !       vertices of the triangulation. For higher-order finite elements
      !       the linearform needs to be assembled by numerical integration
      call storage_getbase_double2d(&
          rproblemLevel%rtriangulation%h_DvertexCoords, p_Dcoords)
      
      ! What type of coordinate system are we?
      select case(icoordsystem)
      case (COORDS_CARTESIAN)
        ! No geometric source term required, clear vector (if required)
        if (bclear) call lsysbl_clearVector(rsource)
        

      case (COORDS_AXIALSYMMETRY)
        !-----------------------------------------------------------------------
        ! Axi-symmetric (dalpha=1) flow (2D approximation to 3D flow)
        !-----------------------------------------------------------------------
        
        select case(igeometricsourcetype)
        case (MASS_LUMPED)
          ! Use lumped mass matrix
          call parlst_getvalue_int(rparlist,&
              ssectionName, 'lumpedmassmatrix', massmatrix)
          call lsyssc_getbase_double(&
              rproblemLevel%Rmatrix(massmatrix), p_DdataMassMatrix)

          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            call doSource2DIntlLumped(dscale, deffectiveRadius,&
                neq, nvar, bclear, p_Dcoords, p_DdataMassMatrix,&
                p_DdataSolution, p_DdataSource)

          case (SYSTEM_BLOCKFORMAT)
            call doSource2DBlockLumped(dscale, deffectiveRadius,&
                neq, nvar, bclear, p_Dcoords, p_DdataMassMatrix,&
                p_DdataSolution, p_DdataSource)

          case default
            call output_line('Invalid system format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcGeometricSourceterm')
            call sys_halt()
          end select
          
        case (MASS_CONSISTENT)
          ! Use consistent mass matrix
          call parlst_getvalue_int(rparlist,&
              ssectionName, 'consistentmassmatrix', massmatrix)
          call lsyssc_getbase_double(&
              rproblemLevel%Rmatrix(massmatrix), p_DdataMassMatrix)
          call lsyssc_getbase_Kld(rproblemLevel%Rmatrix(massmatrix), p_Kld)
          call lsyssc_getbase_Kcol(rproblemLevel%Rmatrix(massmatrix), p_Kcol)
          
          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            call doSource2DIntlConsistent(dscale, deffectiveRadius,&
                neq, nvar, bclear, p_Dcoords, p_DdataMassMatrix, p_Kld,&
                p_Kcol, p_DdataSolution, p_DdataSource)

          case (SYSTEM_BLOCKFORMAT)
            call doSource2DBlockConsistent(dscale, deffectiveRadius,&
                neq, nvar, bclear, p_Dcoords, p_DdataMassMatrix, p_Kld,&
                p_Kcol, p_DdataSolution, p_DdataSource)

          case default
            call output_line('Invalid system format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcGeometricSourceterm')
            call sys_halt()
          end select
          
        case default
          call output_line('Unsupported geometric source type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcGeometricSourceterm')
          call sys_halt()
        end select
        

      case (COORDS_CYLINDRICALSYMMETRY)
        !-----------------------------------------------------------------------
        ! Cylindrically symmetric (dalpha=1) flow (1D approximation to 2D flow)
        !-----------------------------------------------------------------------
        
        select case(igeometricsourcetype)
        case (MASS_LUMPED)
          ! Use lumped mass matrix
          call parlst_getvalue_int(rparlist,&
              ssectionName, 'lumpedmassmatrix', massmatrix)
          call lsyssc_getbase_double(&
              rproblemLevel%Rmatrix(massmatrix), p_DdataMassMatrix)
          
          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            call doSource1DIntlLumped(dscale, deffectiveRadius,&
                neq, nvar, bclear, p_Dcoords, p_DdataMassMatrix,&
                p_DdataSolution, p_DdataSource)
            
          case (SYSTEM_BLOCKFORMAT)
            call doSource1DBlockLumped(dscale, deffectiveRadius,&
                neq, nvar, bclear, p_Dcoords, p_DdataMassMatrix,&
                p_DdataSolution, p_DdataSource)

          case default
            call output_line('Invalid system format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcGeometricSourceterm')
            call sys_halt()
          end select
          
        case (MASS_CONSISTENT)
          ! Use consistent mass matrix
          call parlst_getvalue_int(rparlist,&
              ssectionName, 'consistentmassmatrix', massmatrix)
          call lsyssc_getbase_double(&
              rproblemLevel%Rmatrix(massmatrix), p_DdataMassMatrix)
          call lsyssc_getbase_Kld(rproblemLevel%Rmatrix(massmatrix), p_Kld)
          call lsyssc_getbase_Kcol(rproblemLevel%Rmatrix(massmatrix), p_Kcol)

          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            call doSource1DIntlConsistent(dscale, deffectiveRadius,&
                neq, nvar, bclear, p_Dcoords, p_DdataMassMatrix, p_Kld,&
                p_Kcol, p_DdataSolution, p_DdataSource)

          case (SYSTEM_BLOCKFORMAT)
            call doSource1DBlockConsistent(dscale, deffectiveRadius,&
                neq, nvar, bclear, p_Dcoords, p_DdataMassMatrix, p_Kld,&
                p_Kcol, p_DdataSolution, p_DdataSource)
            
          case default
            call output_line('Invalid system format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcGeometricSourceterm')
            call sys_halt()
          end select
          
        case default
          call output_line('Unsupported geometric source type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcGeometricSourceterm')
          call sys_halt()
        end select
        
        
      case (COORDS_SPHERICALSYMMETRY)
        !-----------------------------------------------------------------------
        ! Spherically symmetric (dalpha=2) flow (1D approximation to 3D flow)
        !-----------------------------------------------------------------------
        
        select case(igeometricsourcetype)
        case (MASS_LUMPED)
          ! Use lumped mass matrix
          call parlst_getvalue_int(rparlist,&
              ssectionName, 'lumpedmassmatrix', massmatrix)
          call lsyssc_getbase_double(&
              rproblemLevel%Rmatrix(massmatrix), p_DdataMassMatrix)

          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            call doSource1DIntlLumped(2*dscale, deffectiveRadius,&
                neq, nvar, bclear, p_Dcoords, p_DdataMassMatrix,&
                p_DdataSolution, p_DdataSource)

          case (SYSTEM_BLOCKFORMAT)
            call doSource1DBlockLumped(2*dscale, deffectiveRadius,&
                neq, nvar, bclear, p_Dcoords, p_DdataMassMatrix,&
                p_DdataSolution, p_DdataSource)

          case default
            call output_line('Invalid system format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcGeometricSourceterm')
            call sys_halt()
          end select
           
        case (MASS_CONSISTENT)
          ! Use consistent mass matrix
          call parlst_getvalue_int(rparlist,&
              ssectionName, 'consistentmassmatrix', massmatrix)
          call lsyssc_getbase_double(&
              rproblemLevel%Rmatrix(massmatrix), p_DdataMassMatrix)
          call lsyssc_getbase_Kld(rproblemLevel%Rmatrix(massmatrix), p_Kld)
          call lsyssc_getbase_Kcol(rproblemLevel%Rmatrix(massmatrix), p_Kcol)

          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            call doSource1DIntlConsistent(2*dscale, deffectiveRadius,&
                neq, nvar, bclear, p_Dcoords, p_DdataMassMatrix, p_Kld,&
                p_Kcol, p_DdataSolution, p_DdataSource)

          case (SYSTEM_BLOCKFORMAT)
            call doSource1DBlockConsistent(2*dscale, deffectiveRadius,&
                neq, nvar, bclear, p_Dcoords, p_DdataMassMatrix, p_Kld,&
                p_Kcol, p_DdataSolution, p_DdataSource)

          case default
            call output_line('Invalid system format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcGeometricSourceterm')
            call sys_halt()
          end select
          
        case default
          call output_line('Unsupported geometric source type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcGeometricSourceterm')
          call sys_halt()
        end select
        
      case default
        call output_line('Invalid coordinate system!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcGeometricSourceterm')
        call sys_halt()
      end select
      
    elseif (trim(smode) .eq. 'dual') then
      
      !-------------------------------------------------------------------------
      ! We are in dual mode
      !-------------------------------------------------------------------------
      
      print *, "Dual mode not implemented yet"
      stop
      
    else
      call output_line('Invalid mode!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcGeometricSourceterm')
      call sys_halt()
    end if

  contains

    ! Here, the real working routines follow

    !**************************************************************
    ! Calculate the geometric source term for cylindrically (dalpha=1)
    ! or spherically symmetric (dalpha=2) flow in 1D. This routine
    ! assembles the geometric source term for systems stored in
    ! interleaved format using the lumped mass matrix.
    
    subroutine doSource1DIntlLumped(deffectiveScale,&
        deffectiveRadius, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar
      
      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:,:), intent(in) :: Dcoords

      ! Lumped mass matrix
      real(DP), dimension(:), intent(in) :: DdataMassMatrix

      ! Solution vector
      real(DP), dimension(nvar,neq), intent(in) :: DdataSolution

      
      ! Source vector
      real(DP), dimension(nvar,neq), intent(inout) :: DdataSource


      ! local variables
      real(DP) :: daux, dradius, dvel, dpre
      integer :: ieq

      
      ! Do we have to clear the source vector?
      if (bclear) then

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(1,ieq); dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = X_VELOCITY_1T_FROM_CONSVAR(DdataSolution, NVAR1D, ieq)
          dpre = PRESSURE_1T_FROM_CONSVAR_1D(DdataSolution, NVAR1D, ieq)

          ! Overwrite the geometric source term
          DdataSource(1,ieq) = daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR1D, ieq)
          DdataSource(2,ieq) = daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR1D, ieq) * dvel
          DdataSource(3,ieq) = daux * DdataMassMatrix(ieq) *&
                               (TOTAL_ENERGY_1T_FROM_CONSVAR(DdataSolution, NVAR1D, ieq)+dpre)*dvel
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(1,ieq); dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = X_VELOCITY_1T_FROM_CONSVAR(DdataSolution, NVAR1D, ieq)
          dpre = PRESSURE_1T_FROM_CONSVAR_1D(DdataSolution, NVAR1D, ieq)

          ! Update the geometric source term
          DdataSource(1,ieq) = DdataSource(1,ieq) + daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR1D, ieq)
          DdataSource(2,ieq) = DdataSource(2,ieq) + daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR1D, ieq) * dvel
          DdataSource(3,ieq) = DdataSource(3,ieq) + daux * DdataMassMatrix(ieq) *&
                               (TOTAL_ENERGY_1T_FROM_CONSVAR(DdataSolution, NVAR1D, ieq)+dpre)*dvel
        end do
        !$omp end parallel do

      end if

    end subroutine doSource1DIntlLumped

    !**************************************************************
    ! Calculate the geometric source term for cylindrically (dalpha=1)
    ! or spherically symmetric (dalpha=2) flow in 1D. This routine
    ! assembles the geometric source term for systems stored in
    ! block format using the lumped mass matrix.
    
    subroutine doSource1DBlockLumped(deffectiveScale,&
        deffectiveRadius, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar
      
      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:,:), intent(in) :: Dcoords

      ! Lumped mass matrix
      real(DP), dimension(:), intent(in) :: DdataMassMatrix

      ! Solution vector
      real(DP), dimension(neq,nvar), intent(in) :: DdataSolution

      
      ! Source vector
      real(DP), dimension(neq,nvar), intent(inout) :: DdataSource


      ! local variables
      real(DP) :: daux, dradius, dvel, dpre
      integer :: ieq

      
      ! Do we have to clear the source vector?
      if (bclear) then

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(1,ieq); dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = X_VELOCITY_1L_FROM_CONSVAR(DdataSolution, NVAR1D, ieq)
          dpre = PRESSURE_1L_FROM_CONSVAR_1D(DdataSolution, NVAR1D, ieq)

          ! Overwrite the geometric source term
          DdataSource(ieq,1) = daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR1D, ieq)
          DdataSource(ieq,2) = daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR1D, ieq) * dvel
          DdataSource(ieq,3) = daux * DdataMassMatrix(ieq) *&
                               (TOTAL_ENERGY_1L_FROM_CONSVAR(DdataSolution, NVAR1D, ieq)+dpre)*dvel
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(1,ieq); dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = X_VELOCITY_1L_FROM_CONSVAR(DdataSolution, NVAR1D, ieq)
          dpre = PRESSURE_1L_FROM_CONSVAR_1D(DdataSolution, NVAR1D, ieq)

          ! Update the geometric source term
          DdataSource(ieq,1) = DdataSource(ieq,1) + daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR1D, ieq)
          DdataSource(ieq,2) = DdataSource(ieq,2) + daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR1D, ieq) * dvel
          DdataSource(ieq,3) = DdataSource(ieq,3) + daux * DdataMassMatrix(ieq) *&
                               (TOTAL_ENERGY_1L_FROM_CONSVAR(DdataSolution, NVAR1D, ieq)+dpre)*dvel
        end do
        !$omp end parallel do

      end if

    end subroutine doSource1DBlockLumped    
    
    !**************************************************************
    ! Calculate the geometric source term for cylindrically (dalpha=1)
    ! or spherically symmetric (dalpha=2) flow in 1D. This routine
    ! assembles the geometric source term for systems stored in
    ! interleaved format using the consistent mass matrix.
    
    subroutine doSource1DIntlConsistent(deffectiveScale,&
        deffectiveRadius, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, Kld, Kcol, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar

      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:,:), intent(in) :: Dcoords

      ! Consistent mass matrix
      real(DP), dimension(:), intent(in) :: DdataMassMatrix

      ! Sparsity pattern of the mass matrix
      integer, dimension(:), intent(in) :: Kld,Kcol

      ! Solution vector
      real(DP), dimension(nvar,neq), intent(in) :: DdataSolution

      
      ! Source vector
      real(DP), dimension(nvar,neq), intent(inout) :: DdataSource


      ! local variables
      real(DP), dimension(NVAR1D) :: Ddata
      real(DP) :: daux, dradius, dvel, dpre
      integer :: ieq,ia,jeq

      
      ! Do we have to clear the source vector?
      if (bclear) then

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(Ddata,daux,dradius,ia,jeq,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)
            
            ! Get the r-coordinate and compute the radius
            daux = Dcoords(1,jeq); dradius = max(abs(daux), deffectiveRadius)
            
            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
            
            ! Compute the radial velocity and pressure
            dvel = X_VELOCITY_1T_FROM_CONSVAR(DdataSolution, NVAR1D, jeq)
            dpre = PRESSURE_1T_FROM_CONSVAR_1D(DdataSolution, NVAR1D, jeq)

            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR1D, jeq)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR1D, jeq) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       (TOTAL_ENERGY_1T_FROM_CONSVAR(DdataSolution, NVAR1D, jeq)+dpre)*dvel
          end do
          
          ! Overwrite the geometric source term
          DdataSource(:,ieq) = Ddata
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(Ddata,daux,dradius,ia,jeq,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)
            
            ! Get the r-coordinate and compute the radius
            daux = Dcoords(1,jeq); dradius = max(abs(daux), deffectiveRadius)

            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
           
            ! Compute the radial velocity and pressure
            dvel = X_VELOCITY_1T_FROM_CONSVAR(DdataSolution, NVAR1D, jeq)
            dpre = PRESSURE_1T_FROM_CONSVAR_1D(DdataSolution, NVAR1D, jeq)
 
            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR1D, jeq)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR1D, jeq) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       (TOTAL_ENERGY_1T_FROM_CONSVAR(DdataSolution, NVAR1D, jeq)+dpre)*dvel
          end do
          
          ! Update the geometric source term
          DdataSource(:,ieq) = DdataSource(:,ieq) + Ddata
        end do
        !$omp end parallel do

      end if
      
    end subroutine doSource1DIntlConsistent
    
    !**************************************************************
    ! Calculate the geometric source term for cylindrically (dalpha=1)
    ! or spherically symmetric (dalpha=2) flow in 1D. This routine
    ! assembles the geometric source term for systems stored in
    ! block format using the consistent mass matrix.
    
    subroutine doSource1DBlockConsistent(deffectiveScale,&
        deffectiveRadius, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, Kld, Kcol, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar

      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:,:), intent(in) :: Dcoords

      ! Consistent mass matrix
      real(DP), dimension(:), intent(in) :: DdataMassMatrix

      ! Sparsity pattern of the mass matrix
      integer, dimension(:), intent(in) :: Kld,Kcol

      ! Solution vector
      real(DP), dimension(neq,nvar), intent(in) :: DdataSolution

      
      ! Source vector
      real(DP), dimension(neq,nvar), intent(inout) :: DdataSource


      ! local variables
      real(DP), dimension(NVAR1D) :: Ddata
      real(DP) :: daux, dradius, dvel, dpre
      integer :: ieq,ia,jeq

      
      ! Do we have to clear the source vector?
      if (bclear) then

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(Ddata,daux,dradius,ia,jeq,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)
            
            ! Get the r-coordinate and compute the radius
            daux = Dcoords(1,jeq); dradius = max(abs(daux), deffectiveRadius)
            
            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
            
            ! Compute the radial velocity and pressure
            dvel = X_VELOCITY_1L_FROM_CONSVAR(DdataSolution, NVAR1D, jeq)
            dpre = PRESSURE_1L_FROM_CONSVAR_1D(DdataSolution, NVAR1D, jeq)

            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR1D, jeq)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR1D, jeq) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       (TOTAL_ENERGY_1L_FROM_CONSVAR(DdataSolution, NVAR1D, jeq)+dpre)*dvel
          end do
          
          ! Overwrite the geometric source term
          DdataSource(ieq,:) = Ddata
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(Ddata,daux,dradius,ia,jeq,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)
            
            ! Get the r-coordinate and compute the radius
            daux = Dcoords(1,jeq); dradius = max(abs(daux), deffectiveRadius)

            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
           
            ! Compute the radial velocity and pressure
            dvel = X_VELOCITY_1L_FROM_CONSVAR(DdataSolution, NVAR1D, jeq)
            dpre = PRESSURE_1L_FROM_CONSVAR_1D(DdataSolution, NVAR1D, jeq)
 
            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR1D, jeq)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR1D, jeq) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       (TOTAL_ENERGY_1L_FROM_CONSVAR(DdataSolution, NVAR1D, jeq)+dpre)*dvel
          end do
          
          ! Update the geometric source term
          DdataSource(ieq,:) = DdataSource(ieq,:) + Ddata
        end do
        !$omp end parallel do

      end if
      
    end subroutine doSource1DBlockConsistent

    !**************************************************************
    ! Calculate the geometric source term for axi-symmetric (dalpha=1)
    ! flow in 2D. This routine assembles the geometric source term for
    ! systems stored in interleaved format using the lumped mass matrix.
    
    subroutine doSource2DIntlLumped(deffectiveScale,&
        deffectiveRadius, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar
      
      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:,:), intent(in) :: Dcoords

      ! Lumped mass matrix
      real(DP), dimension(:), intent(in) :: DdataMassMatrix

      ! Solution vector
      real(DP), dimension(nvar,neq), intent(in) :: DdataSolution

      
      ! Source vector
      real(DP), dimension(nvar,neq), intent(inout) :: DdataSource


      ! local variables
      real(DP) :: daux, dradius, dvel, dpre
      integer :: ieq

      
      ! Do we have to clear the source vector?
      if (bclear) then

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(1,ieq); dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = X_VELOCITY_1T_FROM_CONSVAR(DdataSolution, NVAR2D, ieq)
          dpre = PRESSURE_1T_FROM_CONSVAR_1D(DdataSolution, NVAR2D, ieq)

          ! Overwrite the geometric source term
          DdataSource(1,ieq) = daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR2D, ieq)
          DdataSource(2,ieq) = daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR2D, ieq) * dvel
          DdataSource(3,ieq) = daux * DdataMassMatrix(ieq) *&
                               Y_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR2D, ieq) * dvel
          DdataSource(4,ieq) = daux * DdataMassMatrix(ieq) *&
                               (TOTAL_ENERGY_1T_FROM_CONSVAR(DdataSolution, NVAR2D, ieq)+dpre)*dvel
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(1,ieq); dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = X_VELOCITY_1T_FROM_CONSVAR(DdataSolution, NVAR2D, ieq)
          dpre = PRESSURE_1T_FROM_CONSVAR_1D(DdataSolution, NVAR2D, ieq)

          ! Update the geometric source term
          DdataSource(1,ieq) = DdataSource(1,ieq) + daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR2D, ieq)
          DdataSource(2,ieq) = DdataSource(2,ieq) + daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR2D, ieq) * dvel
          DdataSource(3,ieq) = DdataSource(3,ieq) + daux * DdataMassMatrix(ieq) *&
                               Y_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR2D, ieq) * dvel
          DdataSource(4,ieq) = DdataSource(4,ieq) + daux * DdataMassMatrix(ieq) *&
                               (TOTAL_ENERGY_1T_FROM_CONSVAR(DdataSolution, NVAR2D, ieq)+dpre)*dvel
        end do
        !$omp end parallel do

      end if

    end subroutine doSource2DIntlLumped

    !**************************************************************
    ! Calculate the geometric source term for axi-symmetric (dalpha=1)
    ! flow in 2D. This routine assembles the geometric source term for
    ! systems stored in block format using the lumped mass matrix.
    
    subroutine doSource2DBlockLumped(deffectiveScale,&
        deffectiveRadius, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar
      
      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:,:), intent(in) :: Dcoords

      ! Lumped mass matrix
      real(DP), dimension(:), intent(in) :: DdataMassMatrix

      ! Solution vector
      real(DP), dimension(neq,nvar), intent(in) :: DdataSolution

      
      ! Source vector
      real(DP), dimension(neq,nvar), intent(inout) :: DdataSource


      ! local variables
      real(DP) :: daux, dradius, dvel, dpre
      integer :: ieq

      
      ! Do we have to clear the source vector?
      if (bclear) then

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(1,ieq); dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = X_VELOCITY_1L_FROM_CONSVAR(DdataSolution, NVAR2D, ieq)
          dpre = PRESSURE_1L_FROM_CONSVAR_1D(DdataSolution, NVAR2D, ieq)

          ! Overwrite the geometric source term
          DdataSource(ieq,1) = daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR2D, ieq)
          DdataSource(ieq,2) = daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR2D, ieq) * dvel
          DdataSource(ieq,3) = daux * DdataMassMatrix(ieq) *&
                               Y_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR2D, ieq) * dvel
          DdataSource(ieq,4) = daux * DdataMassMatrix(ieq) *&
                               (TOTAL_ENERGY_1L_FROM_CONSVAR(DdataSolution, NVAR2D, ieq)+dpre)*dvel
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(1,ieq); dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = X_VELOCITY_1L_FROM_CONSVAR(DdataSolution, NVAR2D, ieq)
          dpre = PRESSURE_1L_FROM_CONSVAR_1D(DdataSolution, NVAR2D, ieq)

          ! Update the geometric source term
          DdataSource(ieq,1) = DdataSource(ieq,1) + daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR2D, ieq)
          DdataSource(ieq,2) = DdataSource(ieq,2) + daux * DdataMassMatrix(ieq) *&
                               X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR2D, ieq) * dvel
          DdataSource(ieq,3) = DdataSource(ieq,3) + daux * DdataMassMatrix(ieq) *&
                               Y_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR2D, ieq) * dvel
          DdataSource(ieq,4) = DdataSource(ieq,4) + daux * DdataMassMatrix(ieq) *&
                               (TOTAL_ENERGY_1L_FROM_CONSVAR(DdataSolution, NVAR2D, ieq)+dpre)*dvel
        end do
        !$omp end parallel do

      end if

    end subroutine doSource2DBlockLumped
    
    !**************************************************************
    ! Calculate the geometric source term for axi-symmetric (dalpha=1)
    ! flow in 2D. This routine assembles the geometric source term for
    ! systems stored in interleaved format using the consistent mass matrix.
    
    subroutine doSource2DIntlConsistent(deffectiveScale,&
        deffectiveRadius, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, Kld, Kcol, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar

      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:,:), intent(in) :: Dcoords

      ! Consistent mass matrix
      real(DP), dimension(:), intent(in) :: DdataMassMatrix

      ! Sparsity pattern of the mass matrix
      integer, dimension(:), intent(in) :: Kld,Kcol

      ! Solution vector
      real(DP), dimension(nvar,neq), intent(in) :: DdataSolution

      
      ! Source vector
      real(DP), dimension(nvar,neq), intent(inout) :: DdataSource


      ! local variables
      real(DP), dimension(NVAR2D) :: Ddata
      real(DP) :: daux, dradius, dvel, dpre
      integer :: ieq,ia,jeq

      
      ! Do we have to clear the source vector?
      if (bclear) then

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(Ddata,daux,dradius,ia,jeq,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)
            
            ! Get the r-coordinate and compute the radius
            daux = Dcoords(1,jeq); dradius = max(abs(daux), deffectiveRadius)
            
            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
            
            ! Compute the radial velocity and pressure
            dvel = X_VELOCITY_1T_FROM_CONSVAR(DdataSolution, NVAR2D, jeq)
            dpre = PRESSURE_1T_FROM_CONSVAR_1D(DdataSolution, NVAR2D, jeq)

            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR2D, jeq)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR2D, jeq) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       Y_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR2D, jeq) * dvel
            Ddata(4) = Ddata(4) + daux * DdataMassMatrix(ia) *&
                       (TOTAL_ENERGY_1T_FROM_CONSVAR(DdataSolution, NVAR2D, jeq)+dpre)*dvel
          end do
          
          ! Overwrite the geometric source term
          DdataSource(:,ieq) = Ddata
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(Ddata,daux,dradius,ia,jeq,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)
            
            ! Get the r-coordinate and compute the radius
            daux = Dcoords(1,jeq); dradius = max(abs(daux), deffectiveRadius)

            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
           
            ! Compute the radial velocity and pressure
            dvel = X_VELOCITY_1T_FROM_CONSVAR(DdataSolution, NVAR2D, jeq)
            dpre = PRESSURE_1T_FROM_CONSVAR_1D(DdataSolution, NVAR2D, jeq)
 
            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR2D, jeq)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR2D, jeq) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       Y_MOMENTUM_1T_FROM_CONSVAR(DdataSolution, NVAR2D, jeq) * dvel
            Ddata(4) = Ddata(4) + daux * DdataMassMatrix(ia) *&
                       (TOTAL_ENERGY_1T_FROM_CONSVAR(DdataSolution, NVAR2D, jeq)+dpre)*dvel
          end do
          
          ! Update the geometric source term
          DdataSource(:,ieq) = DdataSource(:,ieq) + Ddata
        end do
        !$omp end parallel do

      end if
      
    end subroutine doSource2DIntlConsistent
    
    !**************************************************************
    ! Calculate the geometric source term for axi-symmetric (dalpha=1)
    ! flow in 2D. This routine assembles the geometric source term for
    ! systems stored in block format using the consistent mass matrix.
    
    subroutine doSource2DBlockConsistent(deffectiveScale,&
        deffectiveRadius, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, Kld, Kcol, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar

      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:,:), intent(in) :: Dcoords

      ! Consistent mass matrix
      real(DP), dimension(:), intent(in) :: DdataMassMatrix

      ! Sparsity pattern of the mass matrix
      integer, dimension(:), intent(in) :: Kld,Kcol

      ! Solution vector
      real(DP), dimension(neq,nvar), intent(in) :: DdataSolution

      
      ! Source vector
      real(DP), dimension(neq,nvar), intent(inout) :: DdataSource


      ! local variables
      real(DP), dimension(NVAR2D) :: Ddata
      real(DP) :: daux, dradius, dvel, dpre
      integer :: ieq,ia,jeq

      
      ! Do we have to clear the source vector?
      if (bclear) then

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(Ddata,daux,dradius,ia,jeq,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)
            
            ! Get the r-coordinate and compute the radius
            daux = Dcoords(1,jeq); dradius = max(abs(daux), deffectiveRadius)
            
            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
            
            ! Compute the radial velocity and pressure
            dvel = X_VELOCITY_1L_FROM_CONSVAR(DdataSolution, NVAR2D, jeq)
            dpre = PRESSURE_1L_FROM_CONSVAR_1D(DdataSolution, NVAR2D, jeq)

            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR2D, jeq)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR2D, jeq) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       Y_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR2D, jeq) * dvel
            Ddata(4) = Ddata(4) + daux * DdataMassMatrix(ia) *&
                       (TOTAL_ENERGY_1L_FROM_CONSVAR(DdataSolution, NVAR2D, jeq)+dpre)*dvel
          end do
          
          ! Overwrite the geometric source term
          DdataSource(ieq,:) = Ddata
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(Ddata,daux,dradius,ia,jeq,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)
            
            ! Get the r-coordinate and compute the radius
            daux = Dcoords(1,jeq); dradius = max(abs(daux), deffectiveRadius)

            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
           
            ! Compute the radial velocity and pressure
            dvel = X_VELOCITY_1L_FROM_CONSVAR(DdataSolution, NVAR2D, jeq)
            dpre = PRESSURE_1L_FROM_CONSVAR_1D(DdataSolution, NVAR2D, jeq)
 
            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR2D, jeq)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       X_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR2D, jeq) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       Y_MOMENTUM_1L_FROM_CONSVAR(DdataSolution, NVAR2D, jeq) * dvel
            Ddata(4) = Ddata(4) + daux * DdataMassMatrix(ia) *&
                       (TOTAL_ENERGY_1L_FROM_CONSVAR(DdataSolution, NVAR2D, jeq)+dpre)*dvel
          end do
          
          ! Update the geometric source term
          DdataSource(ieq,:) = DdataSource(ieq,:) + Ddata
        end do
        !$omp end parallel do

      end if
      
    end subroutine doSource2DBlockConsistent

  end subroutine hydro_calcGeometricSourceterm
end module hydro_callback

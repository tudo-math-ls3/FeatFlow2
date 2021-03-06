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
!# 1.) hydro_calcPreconditioner
!#     -> Calculates the nonlinear preconditioner
!#
!# 2.) hydro_calcJacobian
!#     -> Calculates the Jacobian matrix
!#
!# 3.) hydro_calcResidual
!#     -> Calculates the nonlinear residual vector
!#
!# 4.) hydro_calcRhs
!#     -> Calculates the explicit right-hand side vector
!#
!# 5.) hydro_calcLinearisedFCT
!#     -> Calculates the linearised FCT correction
!#
!# 6.) hydro_calcFluxFCT
!#     -> Calculates the raw antidiffusive fluxes for FCT algorithm
!#
!# 7.) hydro_calcCorrectionFCT
!#     -> Calculates the contribution of the antidiffusive fluxes
!#        limited by the FCT algorithm and applies them to the residual
!#
!# 8.) hydro_limitEdgewiseVelocity
!#     -> Performs synchronised flux correction for the velocity
!#
!# 9.) hydro_limitEdgewiseMomentum
!#     -> Performs synchronised flux correction for the momentum
!#
!# 10.) hydro_calcGeometricSourceterm
!#      -> Calculates the geometric source term for axi-symmetric,
!#         cylindrical or sperical symmetric coordinate systems.
!#
!# 11.) hydro_calcDivergenceVector
!#      -> Calculates the divergence vector.
!#
!# 12.) hydro_calcTimeDerivative
!#      -> Calcuates the approximate time derivative
!#
!# 13.) hydro_coeffVectorFE
!#      -> Callback routine for the evaluation of linear forms
!#         using a given FE-solution for interpolation
!#
!# 14.) hydro_coeffVectorAnalytic
!#      -> Callback routine for the evaluation of linear forms
!#         using a given FE-solution for interpolation
!#
!# 15.) hydro_parseBoundaryCondition
!#      -> Callback routine for the treatment of boundary conditions
!#
!# 16.) hydro_setBoundaryCondition
!#      -> Imposes boundary conditions for nonlinear solver
!#         by filtering the system matrix and the solution/residual
!#         vector explicitly (i.e. strong boundary conditions)
!#
!# 17.) hydro_calcCFLnumber
!#      -> Calculates the CFL-number
!#
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
!#        transp_calcResidual. Note that changing the
!#        interface of this callback routine would require to update
!#        ALL models. Hence, the interface should only be changed if
!#        really necessary.
!#
!# </purpose>
!##############################################################################

module hydro_callback

#include "flagship.h"
#include "hydro.h"
#include "hydro_callback.h"

!$ use omp_lib
  use afcstabbase
  use afcstabsystem
  use afcstabsystemfct
  use afcstabsystemtvd
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
  use solverbase
  use spatialdiscretisation
  use statistics
  use storage
  use timestepbase
  use triangulation

  ! Modules from hydrodynamic model
  use hydro_basic
  use hydro_callback1d
  use hydro_callback2d
  use hydro_callback3d

  implicit none

  private
  public :: hydro_nlsolverCallback
  public :: hydro_calcPreconditioner
  public :: hydro_calcJacobian
  public :: hydro_calcResidual
  public :: hydro_calcRhs
  public :: hydro_setBoundaryCondition
  public :: hydro_calcLinearisedFCT
  public :: hydro_calcFluxFCT
  public :: hydro_calcCorrectionFCT
  public :: hydro_calcGeometricSourceTerm
  public :: hydro_calcDivergenceVector
  public :: hydro_calcTimeDerivative
  public :: hydro_limitEdgewiseVelocity
  public :: hydro_limitEdgewiseMomentum
  public :: hydro_coeffVectorFE
  public :: hydro_coeffVectorAnalytic
  public :: hydro_parseBoundaryCondition
  public :: hydro_calcCFLnumber

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

#ifdef ENABLE_COPROCESSOR_SUPPORT
    ! Copy solution vector to coprocessor device; in the first
    ! nonlinear step, the solution vector from the previous time step
    ! is also required on the device. Thus, we need to transfer it
    call lsysbl_copyH2D_Vector(rsolution, .false., .false., 0_I64)
    if (istep .eq. 0)&
        call lsysbl_copyH2D_Vector(rsolution0, .false., .false., 0_I64)
#endif

    ! Get section name
    call collct_getvalue_string(rcollection,&
        'ssectionname', ssectionName)


    ! Do we have to calculate the preconditioner?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) then

      ! Compute the preconditioner
      call hydro_calcPreconditioner(rproblemLevel, rtimestep,&
          rsolver, rsolution, ssectionName, rcollection)
    end if

    
    ! Do we have to calculate the residual?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      if (istep .eq. 0) then
        ! Compute the constant right-hand side
        call hydro_calcRhs(rproblemLevel, rtimestep,&
            rsolver, rsolution0, rrhs, ssectionName, rcollection,&
            rsource)
      end if

      ! Compute the residual
      call hydro_calcResidual(rproblemLevel, rtimestep,&
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

  subroutine hydro_calcPreconditioner(rproblemLevel, rtimestep,&
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
    ! Check if fully explicit time-stepping is used. Then the system
    ! matrix equals the (lumped/consistent) mass matrix.
    !---------------------------------------------------------------------------
    if (rtimestep%dscaleImplicit .eq. 0.0_DP) then

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
          call lsyssc_spreadDiagMatrix(&
              rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rproblemLevel%RmatrixScalar(systemMatrix))
        case (MASS_CONSISTENT)
          call lsyssc_spreadMatrix(&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rproblemLevel%RmatrixScalar(systemMatrix))
        case default
          call output_line('Empty system matrix is invalid!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPreconditioner')
          call sys_halt()
        end select

      case (SYSTEM_BLOCKFORMAT)

        select case(imasstype)
        case (MASS_LUMPED)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, hydro_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(&
                rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                1.0_DP, 1.0_DP, .false., .false., .true., .true.)
          end do

        case (MASS_CONSISTENT)
          call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
          do ivar = 1, hydro_getNVAR(rproblemLevel)
            call lsyssc_matrixLinearComb(&
                rproblemLevel%RmatrixScalar(consistentMassMatrix),&
                rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
                1.0_DP, 1.0_DP, .false., .false., .true., .true.)
          end do

        case default
          call output_line('Empty system matrix is invalid!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPreconditioner')
          call sys_halt()
        end select

      case default
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPreconditioner')
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
      dscale = -rtimestep%dscaleImplicit*rtimestep%dStep
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

        call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))

        ! Assemble divergence operator without dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatGalMatD1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiagMatD1d_sim, rcollection=rcollection)

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatGalMatD2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiagMatD2d_sim, rcollection=rcollection&
              COPROC_FCB_CALCOPERATOREDGESYS(hydro_calcDivMatGalMatD2d_cuda))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatGalMatD3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiagMatD3d_sim, rcollection=rcollection)

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatScDissMatD1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiagMatD1d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatScDissMatD2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiagMatD2d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))
          
        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatScDissMatD3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiagMatD3d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC)&
              COPROC_FCB_CALCOPERATOREDGESYS(hydro_calcDivMatScDissMatD2d_cuda))

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPrecond')
          call sys_halt()
        end select


      case (DISSIPATION_ROE)

        ! Assemble divergence operator with Roe-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatRoeDissMatD1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiagMatD1d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatRoeDissMatD2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiagMatD2d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC)&
              COPROC_FCB_CALCOPERATOREDGESYS(hydro_calcDivMatRoeDissMatD2d_cuda))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatRoeDissMatD3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiagMatD3d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with the Rusanov-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatRusDissMatD1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiagMatD1d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatRusDissMatD2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiagMatD2d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC)&
              COPROC_FCB_CALCOPERATOREDGESYS(hydro_calcDivMatRusDissMatD2d_cuda))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatRusDissMatD3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiagMatD3d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPreconditioner')
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
              rsolution, hydro_calcMatGalerkin1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiag1d_sim, rcollection=rcollection)

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatGalerkin2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiag2d_sim, rcollection=rcollection&
              COPROC_FCB_CALCOPERATOREDGESYS(hydro_calcDivMatGalerkin2d_cuda))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatGalerkin3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiag3d_sim, rcollection=rcollection)
      
        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_SCALAR)

        ! Assemble divergence operator with scalar dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatScDiss1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiag1d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatScDiss2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiag2d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC)&
              COPROC_FCB_CALCOPERATOREDGESYS(hydro_calcDivMatScDiss2d_cuda))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatScDiss3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiag3d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_ROE)

        ! Assemble divergence operator with Roe-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatRoeDiss1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiag1d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatRoeDiss2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiag2d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC)&
              COPROC_FCB_CALCOPERATOREDGESYS(hydro_calcDivMatRoeDiss2d_cuda))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatRoeDiss3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiag3d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPreconditioner')
          call sys_halt()
        end select


      case (DISSIPATION_RUSANOV)

        ! Assemble divergence operator with the Rusanov-type dissipation

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatRusDiss1d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiag1d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case (NDIM2D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatRusDiss2d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiag2d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC)&
              COPROC_FCB_CALCOPERATOREDGESYS(hydro_calcDivMatRusDiss2d_cuda))

        case (NDIM3D)
          call gfsys_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, hydro_calcMatRusDiss3d_sim, dscale, .true.,&
              rproblemLevel%RmatrixBlock(systemMatrix),&
              hydro_calcMatDiag3d_sim, rcollection=rcollection,&
              rafcstab=rproblemLevel%Rafcstab(inviscidAFC))

        case default
          call output_line('Invalid spatial dimension!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPreconditioner')
          call sys_halt()
        end select


      case default
        ! Clear system matrix and apply (lumped) mass matrix only
        call lsysbl_clearMatrix(rproblemLevel%RmatrixBlock(systemMatrix))
      end select


    case default
      call output_line('Invalid type of flow coupling!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPreconditioner')
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
        !   $ A = blockdiag(M_L) - scaleImplicit * dt * L $
        !
        ! Since we have assembled "-scaleImplicit * dt * L" it
        ! suffices add it to the lumped mass matrix.
        ! -----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)

        call lsyssc_MatrixLinearComb(&
            rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
            rproblemLevel%RmatrixScalar(systemMatrix),&
            1.0_DP, 1.0_DP, .false., .false., .true., .true.)

      case (MASS_CONSISTENT)

        !-----------------------------------------------------------------------
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_C) - scaleImplicit * dt * L $
        !
        ! Since we have assembled "-scaleImplicit * dt * L" it
        ! suffices to add it to the consistent mass matrix.
        ! -----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'consistentmassmatrix', consistentMassMatrix)

        call lsyssc_MatrixLinearComb(&
            rproblemLevel%RmatrixScalar(consistentMassMatrix),&
            rproblemLevel%RmatrixScalar(systemMatrix),&
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
        !   $ A = blockdiag(M_L) - scaleImplicit * dt * L $
        !
        ! Since we have assembled "-scale * dt * L" it suffices to add
        ! it to the lumped mass matrix.
        ! -----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)

        do ivar = 1, hydro_getNVAR(rproblemLevel)
          call lsyssc_MatrixLinearComb(&
              rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rproblemLevel%RmatrixBlock(systemMatrix)%RmatrixBlock(ivar,ivar),&
              1.0_DP, 1.0_DP, .false., .false., .true., .true.)
        end do


      case (MASS_CONSISTENT)

        !-----------------------------------------------------------------------
        ! Compute the global operator for transient flows
        !
        !   $ A = blockdiag(M_C) - scaleExplicit * dt * L $
        !
        ! Since we have assembled "-scaleImplicit * dt * L" it
        ! suffices to add it to the consistent mass matrix.
        ! -----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'consistentmassmatrix', consistentMassMatrix)

        do ivar = 1, hydro_getNVAR(rproblemLevel)
          call lsyssc_MatrixLinearComb(&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
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
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcPreconditioner')
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

  end subroutine hydro_calcPreconditioner

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcJacobian(rproblemLevel, rtimestep,&
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

  end subroutine hydro_calcJacobian

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcRhs(rproblemLevel, rtimestep,&
      rsolver, rsolution, rrhs, ssectionName, rcollection, rsource)

!<description>
    ! This subroutine computes the constant right-hand side
    !
    !  $$ rhs = M*U^n + scaleExplicit * dt * div F(U^n) + S(U^n) + b.c.`s  $$
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
      if (rtimestep%dscaleExplicit .ne. 0.0_DP) then

        ! Compute scaling parameter
        dscale = rtimestep%dscaleExplicit*rtimestep%dStep

        !-----------------------------------------------------------------------
        ! Compute the divergence operator for the right-hand side
        ! evaluated at the solution from the previous(!) iteration
        !
        !   $$ rhs = scaleExplicit * dt * [div F(U^n) + geomSource(U^n) $$
        !-----------------------------------------------------------------------

        call hydro_calcDivergenceVector(rproblemLevel,&
            rsolver%rboundaryCondition, rsolution,&
            rtimestep%dTime-rtimestep%dStep, dscale, .true.,&
            rrhs, ssectionName, rcollection)

        call hydro_calcGeometricSourceterm(p_rparlist, ssectionName,&
            rproblemLevel, rsolution, dscale, .false., rrhs, rcollection)

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
          call lsyssc_matVec(&
              rproblemLevel%RmatrixScalar(massMatrix),&
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
                rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
                rrhs, 1.0_DP, p_rpredictor)
            
            ! Set specifier
            rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec =&
                ior(rproblemLevel%Rafcstab(inviscidAFC)%istabilisationSpec,&
                    AFCSTAB_HAS_PREDICTOR)
            
            ! Assemble explicit part of the raw-antidiffusive fluxes
            call hydro_calcFluxFCT(rproblemLevel, rsolution,&
                rtimestep%dStep, rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                1.0_DP, .true., .true., AFCSTAB_FCTFLUX_EXPLICIT,&
                ssectionName, rcollection, rsolutionPredictor=p_rpredictor)
          end select
        end if

      else ! dscaleExplicit == 0

        !-----------------------------------------------------------------------
        ! Compute the transient term
        !
        !   $$ rhs = M * U^n $$
        !-----------------------------------------------------------------------

        ! What type of mass matrix should be used?
        massMatrix = merge(lumpedMassMatrix,&
            consistentMassMatrix, imasstype .eq. MASS_LUMPED)

        ! Apply mass matrix to solution vector
        do iblock = 1, rsolution%nblocks
          call lsyssc_matVec(&
              rproblemLevel%RmatrixScalar(massMatrix),&
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
            call hydro_calcFluxFCT(rproblemLevel, rsolution,&
                rtimestep%dStep, rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                1.0_DP, .true., .true., AFCSTAB_FCTFLUX_EXPLICIT,&
                ssectionName, rcollection, rsolutionPredictor=p_rpredictor)
          end select
        end if
        
      end if ! dscaleExplicit == 0

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

  end subroutine hydro_calcRhs

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcResidual(rproblemLevel,&
      rtimestep, rsolver, rsolution, rsolution0, rrhs, rres,&
      ite, ssectionName, rcollection, rsource)

!<description>
    ! This subroutine computes the nonlinear residual vector
    !
    ! $$ res^{(m)} = rhs - [M*U^{(m)} - scaleImplicit * dt * div F(U^{(m)})-S^{(m)}-b.c.`s $$
    !
    ! for any two-level time integration scheme, whereby the source
    ! term $S^{(m)}$ is optional. The constant right-hand side
    !
    !  $$ rhs = [M*U^n + scaleExplicit * dt * div F(U^n) + S^n + b.c.`s $$
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
      dscale = rtimestep%dscaleImplicit*rtimestep%dStep
      
      !-----------------------------------------------------------------------
      ! Compute the transient term
      !
      !   $$ res := res - M * U^{(m)} $$
      !-----------------------------------------------------------------------
      
      ! What type of mass matrix should be used?
      massMatrix = merge(lumpedMassMatrix,&
          consistentMassMatrix, imasstype .eq. MASS_LUMPED)

      ! Apply mass matrix to solution vector
      do iblock = 1, rsolution%nblocks
        call lsyssc_matVec(&
            rproblemLevel%RmatrixScalar(massMatrix),&
            rsolution%RvectorBlock(iblock),&
            rres%RvectorBlock(iblock), -1.0_DP, 1.0_DP)
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
    !   $dscale = scaleImplicit * dt$ for transient flows
    !   $dscale = 1$                  for steady-state flows
    !---------------------------------------------------------------------------

    ! Do we have an implicit part?
    if (dscale .ne. 0.0_DP) then
      ! Compute the implicit part of the divergence term
      call hydro_calcDivergenceVector(rproblemLevel,&
          rsolver%rboundaryCondition, rsolution, rtimestep%dTime,&
          dscale, .false., rres, ssectionName, rcollection)
      
      ! Build the geometric source term (if any)
      call hydro_calcGeometricSourceterm(p_rparlist, ssectionName,&
          rproblemLevel, rsolution, dscale, .false., rres, rcollection)
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
      call hydro_calcFluxFCT(rproblemLevel, rsolution,&
          rtimestep%dStep, rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
          1.0_DP, .true., .true., ioperationSpec, ssectionName,&
          rcollection, rsolutionPredictor=p_rpredictor)
      
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
      call hydro_calcCorrectionFCT(rproblemLevel, p_rpredictor,&
          rtimestep%dStep, .false., ioperationSpec, rres,&
          ssectionName, rcollection)

      ! Special treatment for iterative FCT-algorithm
      if (rproblemLevel%Rafcstab(inviscidAFC)%cafcstabType .eq.&
          AFCSTAB_NLINFCT_ITERATIVE) then
        ! Subtract corrected antidiffusion from right-hand side
        call afcsys_buildVectorFCT(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
            p_rpredictor, rtimestep%dStep, .false.,&
            AFCSTAB_FCTALGO_CORRECT, rrhs,&
            rcollection=rcollection)

        ! Recompute the low-order predictor for the next limiting step
        call lsysbl_invertedDiagMatVec(&
            rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
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

  end subroutine hydro_calcResidual

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
    select case(rsolver%ipreconditioner)
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

    case default
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_setBoundaryCondition')
      call sys_halt()
    end select

  end subroutine hydro_setBoundaryCondition
  
  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcLinearisedFCT(rproblemLevel, rtimestep,&
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
        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'nfailsafe', nfailsafe)
        
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          
          ! Set up vector for computing the approximate time derivative
          if (present(rvector1)) then
            p_rvector1 => rvector1
          else
            allocate(p_rvector1)
          end if
          
          ! Compute approximate time derivative
          call hydro_calcTimeDerivative(rproblemLevel, rtimestep,&
              rsolver, rsolution, ssectionName, rcollection, p_rvector1,&
              rsource, rvector2, rvector3)
          
          ! Build the raw antidiffusive fluxes and include
          ! contribution from the consistent mass matrix
          call hydro_calcFluxFCT(rproblemLevel, rsolution,&
              1.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, .true., .true.,&
              AFCSTAB_FCTFLUX_EXPLICIT, ssectionName,&
              rcollection, rsolutionTimeDeriv=p_rvector1)
          
          ! Release temporal memory
          if (.not.present(rvector1)) then
            call lsysbl_releaseVector(p_rvector1)
            deallocate(p_rvector1)
          end if
      
        else
          
          ! Build the raw antidiffusive fluxes without including
          ! the contribution from consistent mass matrix
          call hydro_calcFluxFCT(rproblemLevel, rsolution,&
              1.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, .true., .true.,&
              AFCSTAB_FCTFLUX_EXPLICIT, ssectionName, rcollection)
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
          
          ! Compute FEM-FCT correction without applying the
          ! antidiffusive correction term to the low-order solution
          call hydro_calcCorrectionFCT(rproblemLevel,&
              rsolution, rtimestep%dStep, .false.,&
              AFCSTAB_FCTALGO_STANDARD-&
              AFCSTAB_FCTALGO_CORRECT,&
              rsolution, ssectionName, rcollection)
          
          ! Apply failsafe flux correction
          call afcsys_failsafeFCT(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, rtimestep%dStep, 1e-8_DP,&
              AFCSTAB_FAILSAFEALGO_STANDARD, bisAccepted,&
              nsteps=nfailsafe, CvariableNames=SfailsafeVariables,&
              fcb_extractVariable=hydro_getVariable,&
              rcollection=rcollection)
                    
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

      end select
    end if

        
    ! Impose boundary conditions for the solution vector
    select case(rproblemLevel%rtriangulation%ndim)
    case (NDIM1D)
      call bdrf_filterVectorExplicit(rsolver%rboundaryCondition,&
          rsolution, rtimestep%dTime, hydro_calcBoundaryvalues1d)

    case (NDIM2D)
      call bdrf_filterVectorExplicit(rsolver%rboundaryCondition,&
          rsolution, rtimestep%dTime, hydro_calcBoundaryvalues2d)

    case (NDIM3D)
      call bdrf_filterVectorExplicit(rsolver%rboundaryCondition,&
          rsolution, rtimestep%dTime, hydro_calcBoundaryvalues3d)
    end select

  end subroutine hydro_calcLinearisedFCT

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcFluxFCT(rproblemLevel, rsolution, tstep,&
      dscaleExplicit, dscaleImplicit, dscale, bclear, bquickAssembly,&
      ioperationSpec, ssectionName, rcollection,&
      rsolutionTimeDeriv, rsolutionPredictor)

!<description>
    ! This subroutine calculates the raw antidiffusive fluxes for
    ! the different FCT algorithms
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! time step size
    real(DP), intent(in) :: tstep

    ! scaling factors
    real(DP), intent(in) :: dscaleExplicit
    real(DP), intent(in) :: dscaleImplicit
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
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTScDiss1d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTScDiss1d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rcollection=rcollection)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTScDiss2d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTScDiss2d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rcollection=rcollection)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTScDiss3d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTScDiss3d_sim,&
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
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTRoeDiss1d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTRoeDiss1d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rcollection=rcollection)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTRoeDiss2d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTRoeDiss2d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rcollection=rcollection)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTRoeDiss3d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTRoeDiss3d_sim,&
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
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTRusDiss1d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTRusDiss1d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        end if

      case (NDIM2D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTRusDiss2d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTRusDiss2d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        end if

      case (NDIM3D)
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTRusDiss3d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rxTimeDeriv=rsolutionTimeDeriv,&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        else
          call afcsys_buildFluxFCT(&
              rproblemLevel%Rafcstab(inviscidAFC), rsolution,&
              tstep, dscaleExplicit, dscaleImplicit, dscale, bclear,&
              bquickAssembly, ioperationSpec, hydro_calcFluxFCTRusDiss3d_sim,&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rxPredictor=rsolutionPredictor,&
              rcollection=rcollection)
        end if
      end select


    case default
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
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix), &
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDensity1d_sim, hydro_trafoDiffDensity1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDensity2d_sim, hydro_trafoDiffDensity2d_sim,&
              rcollection=rcollection)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDensity3d_sim, hydro_trafoDiffDensity3d_sim,&
              rcollection=rcollection)
        end select

      elseif (trim(slimitingvariable) .eq. 'energy') then

        ! Apply FEM-FCT algorithm for energy fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxEnergy1d_sim, hydro_trafoDiffEnergy1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxEnergy2d_sim, hydro_trafoDiffEnergy2d_sim,&
              rcollection=rcollection)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxEnergy3d_sim, hydro_trafoDiffEnergy3d_sim,&
              rcollection=rcollection)
        end select

      elseif (trim(slimitingvariable) .eq. 'pressure') then

        ! Apply FEM-FCT algorithm for pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxPressure1d_sim, hydro_trafoDiffPressure1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxPressure2d_sim, hydro_trafoDiffPressure2d_sim,&
              rcollection=rcollection)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxPressure3d_sim, hydro_trafoDiffPressure3d_sim,&
              rcollection=rcollection)
        end select

      elseif (trim(slimitingvariable) .eq. 'velocity') then

        ! Apply FEM-FCT algorithm for velocity fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxVelocity1d_sim, hydro_trafoDiffVelocity1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxVelocity2d_sim, hydro_trafoDiffVelocity2d_sim,&
              rcollection=rcollection,&
              fcb_limitEdgewise=hydro_limitEdgewiseVelocity)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxVelocity3d_sim, hydro_trafoDiffVelocity3d_sim,&
              rcollection=rcollection,&
              fcb_limitEdgewise=hydro_limitEdgewiseVelocity)
        end select

      elseif (trim(slimitingvariable) .eq. 'momentum') then

        ! Apply FEM-FCT algorithm for momentum fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxMomentum1d_sim, hydro_trafoDiffMomentum1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxMomentum2d_sim, hydro_trafoDiffMomentum2d_sim,&
              rcollection=rcollection,&
              fcb_limitEdgewise=hydro_limitEdgewiseMomentum)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxMomentum3d_sim, hydro_trafoDiffMomentum3d_sim,&
              rcollection=rcollection,&
              fcb_limitEdgewise=hydro_limitEdgewiseMomentum)
        end select
        
      elseif (trim(slimitingvariable) .eq. 'density,energy') then

        ! Apply FEM-FCT algorithm for density and energy fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenEng1d_sim, hydro_trafoDiffDenEng1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenEng2d_sim, hydro_trafoDiffDenEng2d_sim,&
              rcollection=rcollection)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenEng3d_sim, hydro_trafoDiffDenEng3d_sim,&
              rcollection=rcollection)
        end select

      elseif (trim(slimitingvariable) .eq. 'density,pressure') then

        ! Apply FEM-FCT algorithm for density and pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenPre1d_sim, hydro_trafoDiffDenPre1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenPre2d_sim, hydro_trafoDiffDenPre2d_sim,&
              rcollection=rcollection)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenPre3d_sim, hydro_trafoDiffDenPre3d_sim,&
              rcollection=rcollection)
        end select

      elseif (trim(slimitingvariable) .eq. 'density,energy,momentum') then

        nvartransformed = hydro_getNVARtransformed(rproblemLevel, slimitingvariable)

        ! Apply FEM-FCT algorithm for full conservative fluxes
        call afcsys_buildVectorFCT(&
            p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
            rsolution, dscale, bclear, iopSpec, rresidual, rcollection=rcollection)

      elseif (trim(slimitingvariable) .eq. 'density,pressure,velocity') then

        nvartransformed = hydro_getNVARtransformed(rproblemLevel, slimitingvariable)

        ! Apply FEM-FCT algorithm for density, velocity and pressure fluxes
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenPreVel1d_sim, hydro_trafoDiffDenPreVel1d_sim,&
              rcollection=rcollection)
        case (NDIM2D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenPreVel2d_sim, hydro_trafoDiffDenPreVel2d_sim,&
              rcollection=rcollection)
        case (NDIM3D)
          call afcsys_buildVectorFCT(&
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual, nvartransformed,&
              hydro_trafoFluxDenPreVel3d_sim, hydro_trafoDiffDenPreVel3d_sim,&
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
              p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rsolution, dscale, bclear, iopSpec, rresidual,&
              rcollection=rcollection)
          
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
      
      call afcsys_buildVectorFCT(&
          p_rafcstab, rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
          rsolution, dscale, bclear, iopSpec, rresidual,&
          rcollection=rcollection)
    end if

  end subroutine hydro_calcCorrectionFCT

  !***************************************************************************

!<subroutine>

  subroutine hydro_limitEdgewiseVelocity(IedgeListIdx, IedgeList,&
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
    include '../../../../../kernel/PDEOperators/intf_calcFluxTransform_sim.inc'
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
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
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
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
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
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_limitEdgewiseVelocity')
        call sys_halt()
        
      end if
    end if
    
  end subroutine hydro_limitEdgewiseVelocity

  !***************************************************************************

!<subroutine>

  subroutine hydro_limitEdgewiseMomentum(IedgeListIdx, IedgeList,&
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
    include '../../../../../kernel/PDEOperators/intf_calcFluxTransform_sim.inc'
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
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
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
              DdataAtEdge(:,:,1:IEDGEmax-IEDGEset+1),&
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
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_limitEdgewiseMomentum')
        call sys_halt()
        
      end if
    end if

  end subroutine hydro_limitEdgewiseMomentum

  !***************************************************************************

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

  !***************************************************************************

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
    
  end subroutine hydro_coeffVectorAnalytic

  !*****************************************************************************

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

    ! Free-selip boundary conditions
    case ('FREESLIP_STRONG')
      ibdrCondType = BDRC_FREESLIP + BDRC_STRONG
    case ('FREESLIP_WEAK')
      ibdrCondType = BDRC_FREESLIP + BDRC_WEAK
    case ('FREESLIP_WEAK_LUMPED')
      ibdrCondType = BDRC_FREESLIP + BDRC_WEAK + BDRC_LUMPED
      
    ! Relaxed free-slip boundary conditions
    case ('RLXFREESLIP_STRONG')
      ibdrCondType = BDRC_RLXFREESLIP + BDRC_STRONG
    case ('RLXFREESLIP_WEAK')
      ibdrCondType = BDRC_RLXFREESLIP + BDRC_WEAK
    case ('RLXFREESLIP_WEAK_LUMPED')
      ibdrCondType = BDRC_RLXFREESLIP + BDRC_WEAK + BDRC_LUMPED

    ! Viscous wall boundary conditions
    case ('VISCOUSWALL_STRONG')
      ibdrCondType = BDRC_VISCOUSWALL + BDRC_STRONG
    case ('VISCOUSWALL_WEAK')
      ibdrCondType = BDRC_VISCOUSWALL + BDRC_WEAK
    case ('VISCOUSWALL_WEAK_LUMPED')
      ibdrCondType = BDRC_VISCOUSWALL + BDRC_WEAK + BDRC_LUMPED

    ! Freestream boundary conditions
    case ('FREESTREAM_STRONG')
      ibdrCondType = BDRC_FREESTREAM + BDRC_STRONG     
    case ('FREESTREAM_WEAK')
      ibdrCondType = BDRC_FREESTREAM + BDRC_WEAK
    case ('FREESTREAM_WEAK_LUMPED')
      ibdrCondType = BDRC_FREESTREAM + BDRC_WEAK + BDRC_LUMPED

    ! Subsonic inlet boundary conditions
    case ('SUBINLET_STRONG')
      ibdrCondType = BDRC_SUBINLET + BDRC_STRONG
    case ('SUBINLET_WEAK')
      ibdrCondType = BDRC_SUBINLET + BDRC_WEAK
    case ('SUBINLET_WEAK_LUMPED')
      ibdrCondType = BDRC_SUBINLET + BDRC_WEAK + BDRC_LUMPED

    ! Subsonic outlet boundary conditions
    case ('SUBOUTLET_STRONG')
      ibdrCondType = BDRC_SUBOUTLET + BDRC_STRONG
    case ('SUBOUTLET_WEAK')
      ibdrCondType = BDRC_SUBOUTLET + BDRC_WEAK
    case ('SUBOUTLET_WEAK_LUMPED')
      ibdrCondType = BDRC_SUBOUTLET + BDRC_WEAK + BDRC_LUMPED

    ! Massflow inlet boundary conditions
    case ('MASSINLET_STRONG')
      ibdrCondType = BDRC_MASSINLET + BDRC_STRONG
    case ('MASSINLET_WEAK')
      ibdrCondType = BDRC_MASSINLET + BDRC_WEAK
    case ('MASSINLET_WEAK_LUMPED')
      ibdrCondType = BDRC_MASSINLET + BDRC_WEAK + BDRC_LUMPED

    ! Massflow outlet boundary conditions
    case ('MASSOUTLET_STRONG')
      ibdrCondType = BDRC_MASSOUTLET + BDRC_STRONG
    case ('MASSOUTLET_WEAK')
      ibdrCondType = BDRC_MASSOUTLET + BDRC_WEAK
    case ('MASSOUTLET_WEAK_LUMPED')
      ibdrCondType = BDRC_MASSOUTLET + BDRC_WEAK + BDRC_LUMPED
      
    ! Mach outlet boundary conditions
    case ('MACHOUTLET_STRONG')
      ibdrCondType = BDRC_MACHOUTLET + BDRC_STRONG
    case ('MACHOUTLET_WEAK')
      ibdrCondType = BDRC_MACHOUTLET + BDRC_WEAK
    case ('MACHOUTLET_WEAK_LUMPED')
      ibdrCondType = BDRC_MACHOUTLET + BDRC_WEAK + BDRC_LUMPED

    ! Supersonic inlet boundary conditions
    case ('SUPERINLET_STRONG')
      ibdrCondType = BDRC_SUPERINLET + BDRC_STRONG
    case ('SUPERINLET_WEAK')
      ibdrCondType = BDRC_SUPERINLET + BDRC_WEAK
    case ('SUPERINLET_WEAK_LUMPED')
      ibdrCondType = BDRC_SUPERINLET + BDRC_WEAK + BDRC_LUMPED

    ! Supersonic outlet boundary conditions
    case ('SUPEROUTLET_STRONG')
      ibdrCondType = BDRC_SUPEROUTLET
      ! No strong boundary conditions are prescribed
    case ('SUPEROUTLET_WEAK')
      ibdrCondType = BDRC_SUPEROUTLET + BDRC_WEAK
    case ('SUPEROUTLET_WEAK_LUMPED')
      ibdrCondType = BDRC_SUPEROUTLET + BDRC_WEAK + BDRC_LUMPED

    ! Open boundary conditions
    case ('OPEN_STRONG')
      ibdrCondType = BDRC_OPEN
      ! No strong boundary conditions are prescribed
    case ('OPEN_WEAK')
      ibdrCondType = BDRC_OPEN + BDRC_WEAK
    case ('OPEN_WEAK_LUMPED')
      ibdrCondType = BDRC_OPEN + BDRC_WEAK + BDRC_LUMPED
   
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
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_parseBoundaryCondition')
      call sys_halt()
    end select

    
    ! Determine number of mathematical expressions
    select case (iand(ibdrCondType, BDRC_TYPEMASK))
    case (BDRC_OPEN)
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
      nexpressions = 2 ! penalty parameter $\epsilon$
                       ! switching parameter $\gamma$

    case default
      nexpressions = 0
    end select

  end subroutine hydro_parseBoundaryCondition

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
    real(DP), dimension(:), pointer :: p_DdofCoords
    real(DP), dimension(:), pointer :: p_DdataSolution
    real(DP), dimension(:), pointer :: p_DdataSource
    real(DP), dimension(:), pointer :: p_DdataMassMatrix
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    character(LEN=SYS_STRLEN) :: smode
    real(DP) :: deffectiveRadius
    integer :: neq, nvar, icoordsystem, igeometricsourcetype
    integer :: isystemFormat, massmatrix, dofCoords

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
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'dofCoords', dofCoords, 0)
      
      ! Get pointers
      call lsysbl_getbase_double(rsolution, p_DdataSolution)
      call lsysbl_getbase_double(rsource, p_DdataSource)

      ! Get number of equations and variables
      nvar = hydro_getNVAR(rproblemLevel)
      neq  = rsolution%NEQ/nvar
      
      ! Get coordinates of the global DOF`s
      call lsysbl_getbase_double(&
          rproblemLevel%RvectorBlock(dofCoords), p_DdofCoords)
      
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
              rproblemLevel%RmatrixScalar(massmatrix), p_DdataMassMatrix)

          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            call doSource2DIntlLumped(dscale, deffectiveRadius,&
                rproblemLevel%rtriangulation%ndim, neq, nvar, bclear,&
                p_DdofCoords, p_DdataMassMatrix,&
                p_DdataSolution, p_DdataSource)

          case (SYSTEM_BLOCKFORMAT)
            call doSource2DBlockLumped(dscale, deffectiveRadius,&
                rproblemLevel%rtriangulation%ndim, neq, nvar, bclear,&
                p_DdofCoords, p_DdataMassMatrix,&
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
              rproblemLevel%RmatrixScalar(massmatrix), p_DdataMassMatrix)
          call lsyssc_getbase_Kld(rproblemLevel%RmatrixScalar(massmatrix), p_Kld)
          call lsyssc_getbase_Kcol(rproblemLevel%RmatrixScalar(massmatrix), p_Kcol)
          
          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            call doSource2DIntlConsistent(dscale, deffectiveRadius,&
                rproblemLevel%rtriangulation%ndim, neq, nvar, bclear,&
                p_DdofCoords, p_DdataMassMatrix, p_Kld, p_Kcol,&
                p_DdataSolution, p_DdataSource)

          case (SYSTEM_BLOCKFORMAT)
            call doSource2DBlockConsistent(dscale, deffectiveRadius,&
                rproblemLevel%rtriangulation%ndim, neq, nvar, bclear,&
                p_DdofCoords, p_DdataMassMatrix, p_Kld, p_Kcol,&
                p_DdataSolution, p_DdataSource)

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
              rproblemLevel%RmatrixScalar(massmatrix), p_DdataMassMatrix)
          
          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            call doSource1DIntlLumped(dscale, deffectiveRadius,&
                rproblemLevel%rtriangulation%ndim, neq, nvar, bclear,&
                p_DdofCoords, p_DdataMassMatrix,&
                p_DdataSolution, p_DdataSource)
            
          case (SYSTEM_BLOCKFORMAT)
            call doSource1DBlockLumped(dscale, deffectiveRadius,&
                rproblemLevel%rtriangulation%ndim, neq, nvar, bclear,&
                p_DdofCoords, p_DdataMassMatrix,&
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
              rproblemLevel%RmatrixScalar(massmatrix), p_DdataMassMatrix)
          call lsyssc_getbase_Kld(rproblemLevel%RmatrixScalar(massmatrix), p_Kld)
          call lsyssc_getbase_Kcol(rproblemLevel%RmatrixScalar(massmatrix), p_Kcol)

          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            call doSource1DIntlConsistent(dscale, deffectiveRadius,&
                rproblemLevel%rtriangulation%ndim, neq, nvar, bclear,&
                p_DdofCoords, p_DdataMassMatrix, p_Kld, p_Kcol,&
                p_DdataSolution, p_DdataSource)

          case (SYSTEM_BLOCKFORMAT)
            call doSource1DBlockConsistent(dscale, deffectiveRadius,&
                rproblemLevel%rtriangulation%ndim, neq, nvar, bclear,&
                p_DdofCoords, p_DdataMassMatrix, p_Kld, p_Kcol,&
                p_DdataSolution, p_DdataSource)
            
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
              rproblemLevel%RmatrixScalar(massmatrix), p_DdataMassMatrix)

          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            call doSource1DIntlLumped(2*dscale, deffectiveRadius,&
                rproblemLevel%rtriangulation%ndim, neq, nvar, bclear,&
                p_DdofCoords, p_DdataMassMatrix,&
                p_DdataSolution, p_DdataSource)

          case (SYSTEM_BLOCKFORMAT)
            call doSource1DBlockLumped(2*dscale, deffectiveRadius,&
                rproblemLevel%rtriangulation%ndim, neq, nvar, bclear,&
                p_DdofCoords, p_DdataMassMatrix,&
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
              rproblemLevel%RmatrixScalar(massmatrix), p_DdataMassMatrix)
          call lsyssc_getbase_Kld(rproblemLevel%RmatrixScalar(massmatrix), p_Kld)
          call lsyssc_getbase_Kcol(rproblemLevel%RmatrixScalar(massmatrix), p_Kcol)

          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)
            call doSource1DIntlConsistent(2*dscale, deffectiveRadius,&
                rproblemLevel%rtriangulation%ndim, neq, nvar, bclear,&
                p_DdofCoords, p_DdataMassMatrix, p_Kld, p_Kcol,&
                p_DdataSolution, p_DdataSource)

          case (SYSTEM_BLOCKFORMAT)
            call doSource1DBlockConsistent(2*dscale, deffectiveRadius,&
                rproblemLevel%rtriangulation%ndim, neq, nvar, bclear,&
                p_DdofCoords, p_DdataMassMatrix, p_Kld, p_Kcol,&
                p_DdataSolution, p_DdataSource)

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
        deffectiveRadius, ndim, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Spatial dimension
      integer, intent(in) :: ndim

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar
      
      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:), intent(in) :: Dcoords

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
          daux = Dcoords(ndim*(ieq-1)+1)
          dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = XVELOCITY2_1D(DdataSolution,IDX2_FORWARD,ieq,_,_)
          dpre = PRESSURE2_1D(DdataSolution,IDX2_FORWARD,ieq,_,_)

          ! Overwrite the geometric source term
          DdataSource(1,ieq) = daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_1D(DdataSolution,IDX2_FORWARD,ieq,_,_)
          DdataSource(2,ieq) = daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_1D(DdataSolution,IDX2_FORWARD,ieq,_,_) * dvel
          DdataSource(3,ieq) = daux * DdataMassMatrix(ieq) *&
                               (TOTALENERGY2_1D(DdataSolution,IDX2_FORWARD,ieq,_,_)+dpre)*dvel
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(ndim*(ieq-1)+1)
          dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = XVELOCITY2_1D(DdataSolution,IDX2_FORWARD,ieq,_,_)
          dpre = PRESSURE2_1D(DdataSolution,IDX2_FORWARD,ieq,_,_)

          ! Update the geometric source term
          DdataSource(1,ieq) = DdataSource(1,ieq) + daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_1D(DdataSolution,IDX2_FORWARD,ieq,_,_)
          DdataSource(2,ieq) = DdataSource(2,ieq) + daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_1D(DdataSolution,IDX2_FORWARD,ieq,_,_) * dvel
          DdataSource(3,ieq) = DdataSource(3,ieq) + daux * DdataMassMatrix(ieq) *&
                               (TOTALENERGY2_1D(DdataSolution,IDX2_FORWARD,ieq,_,_)+dpre)*dvel
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
        deffectiveRadius, ndim, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Spatial dimension
      integer, intent(in) :: ndim

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar
      
      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:), intent(in) :: Dcoords

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
          daux = Dcoords(ndim*(ieq-1)+1)
          dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = XVELOCITY2_1D(DdataSolution,IDX2_REVERSE,ieq,_,_)
          dpre = PRESSURE2_1D(DdataSolution,IDX2_REVERSE,ieq,_,_)

          ! Overwrite the geometric source term
          DdataSource(ieq,1) = daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_1D(DdataSolution,IDX2_REVERSE,ieq,_,_)
          DdataSource(ieq,2) = daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_1D(DdataSolution,IDX2_REVERSE,ieq,_,_) * dvel
          DdataSource(ieq,3) = daux * DdataMassMatrix(ieq) *&
                               (TOTALENERGY2_1D(DdataSolution,IDX2_REVERSE,ieq,_,_)+dpre)*dvel
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(ndim*(ieq-1)+1)
          dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = XVELOCITY2_1D(DdataSolution,IDX2_REVERSE,ieq,_,_)
          dpre = PRESSURE2_1D(DdataSolution,IDX2_REVERSE,ieq,_,_)

          ! Update the geometric source term
          DdataSource(ieq,1) = DdataSource(ieq,1) + daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_1D(DdataSolution,IDX2_REVERSE,ieq,_,_)
          DdataSource(ieq,2) = DdataSource(ieq,2) + daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_1D(DdataSolution,IDX2_REVERSE,ieq,_,_) * dvel
          DdataSource(ieq,3) = DdataSource(ieq,3) + daux * DdataMassMatrix(ieq) *&
                               (TOTALENERGY2_1D(DdataSolution,IDX2_REVERSE,ieq,_,_)+dpre)*dvel
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
        deffectiveRadius, ndim, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, Kld, Kcol, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Spatial dimension
      integer, intent(in) :: ndim

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar

      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:), intent(in) :: Dcoords

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
            daux = Dcoords(ndim*(jeq-1)+1)
            dradius = max(abs(daux), deffectiveRadius)
            
            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
            
            ! Compute the radial velocity and pressure
            dvel = XVELOCITY2_1D(DdataSolution,IDX2_FORWARD,jeq,_,_)
            dpre = PRESSURE2_1D(DdataSolution,IDX2_FORWARD,jeq,_,_)

            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_1D(DdataSolution,IDX2_FORWARD,jeq,_,_)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_1D(DdataSolution,IDX2_FORWARD,jeq,_,_) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       (TOTALENERGY2_1D(DdataSolution,IDX2_FORWARD,jeq,_,_)+dpre)*dvel
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
            daux = Dcoords(ndim*(jeq-1)+1)
            dradius = max(abs(daux), deffectiveRadius)

            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
           
            ! Compute the radial velocity and pressure
            dvel = XVELOCITY2_1D(DdataSolution,IDX2_FORWARD,jeq,_,_)
            dpre = PRESSURE2_1D(DdataSolution,IDX2_FORWARD,jeq,_,_)
 
            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_1D(DdataSolution,IDX2_FORWARD,jeq,_,_)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_1D(DdataSolution,IDX2_FORWARD,jeq,_,_) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       (TOTALENERGY2_1D(DdataSolution,IDX2_FORWARD,jeq,_,_)+dpre)*dvel
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
        deffectiveRadius, ndim, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, Kld, Kcol, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Spatial dimension
      integer, intent(in) :: ndim
      
      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar

      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:), intent(in) :: Dcoords

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
            daux = Dcoords(ndim*(jeq-1)+1)
            dradius = max(abs(daux), deffectiveRadius)
            
            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
            
            ! Compute the radial velocity and pressure
            dvel = XVELOCITY2_1D(DdataSolution,IDX2_REVERSE,jeq,_,_)
            dpre = PRESSURE2_1D(DdataSolution,IDX2_REVERSE,jeq,_,_)

            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_1D(DdataSolution,IDX2_REVERSE,jeq,_,_)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_1D(DdataSolution,IDX2_REVERSE,jeq,_,_) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       (TOTALENERGY2_1D(DdataSolution,IDX2_REVERSE,jeq,_,_)+dpre)*dvel
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
            daux = Dcoords(ndim*(jeq-1)+1)
            dradius = max(abs(daux), deffectiveRadius)

            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
           
            ! Compute the radial velocity and pressure
            dvel = XVELOCITY2_1D(DdataSolution,IDX2_REVERSE,jeq,_,_)
            dpre = PRESSURE2_1D(DdataSolution,IDX2_REVERSE,jeq,_,_)
 
            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_1D(DdataSolution,IDX2_REVERSE,jeq,_,_)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_1D(DdataSolution,IDX2_REVERSE,jeq,_,_) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       (TOTALENERGY2_1D(DdataSolution,IDX2_REVERSE,jeq,_,_)+dpre)*dvel
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
        deffectiveRadius, ndim, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius
      
      ! Spatial dimension
      integer, intent(in) :: ndim
      
      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar
      
      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:), intent(in) :: Dcoords

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
          daux = Dcoords(ndim*(ieq-1)+1)
          dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = XVELOCITY2_2D(DdataSolution,IDX2_FORWARD,ieq,_,_)
          dpre = PRESSURE2_2D(DdataSolution,IDX2_FORWARD,ieq,_,_)

          ! Overwrite the geometric source term
          DdataSource(1,ieq) = daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_2D(DdataSolution,IDX2_FORWARD,ieq,_,_)
          DdataSource(2,ieq) = daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_2D(DdataSolution,IDX2_FORWARD,ieq,_,_) * dvel
          DdataSource(3,ieq) = daux * DdataMassMatrix(ieq) *&
                               YMOMENTUM2_2D(DdataSolution,IDX2_FORWARD,ieq,_,_) * dvel
          DdataSource(4,ieq) = daux * DdataMassMatrix(ieq) *&
                               (TOTALENERGY2_2D(DdataSolution,IDX2_FORWARD,ieq,_,_)+dpre)*dvel
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(ndim*(ieq-1)+1)
          dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = XVELOCITY2_2D(DdataSolution,IDX2_FORWARD,ieq,_,_)
          dpre = PRESSURE2_2D(DdataSolution,IDX2_FORWARD,ieq,_,_)

          ! Update the geometric source term
          DdataSource(1,ieq) = DdataSource(1,ieq) + daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_2D(DdataSolution,IDX2_FORWARD,ieq,_,_)
          DdataSource(2,ieq) = DdataSource(2,ieq) + daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_2D(DdataSolution,IDX2_FORWARD,ieq,_,_) * dvel
          DdataSource(3,ieq) = DdataSource(3,ieq) + daux * DdataMassMatrix(ieq) *&
                               YMOMENTUM2_2D(DdataSolution,IDX2_FORWARD,ieq,_,_) * dvel
          DdataSource(4,ieq) = DdataSource(4,ieq) + daux * DdataMassMatrix(ieq) *&
                               (TOTALENERGY2_2D(DdataSolution,IDX2_FORWARD,ieq,_,_)+dpre)*dvel
        end do
        !$omp end parallel do

      end if

    end subroutine doSource2DIntlLumped

    !**************************************************************
    ! Calculate the geometric source term for axi-symmetric (dalpha=1)
    ! flow in 2D. This routine assembles the geometric source term for
    ! systems stored in block format using the lumped mass matrix.
    
    subroutine doSource2DBlockLumped(deffectiveScale,&
        deffectiveRadius, ndim, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Spatial dimension
      integer, intent(in) :: ndim

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar
      
      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:), intent(in) :: Dcoords

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
          daux = Dcoords(ndim*(ieq-1)+1)
          dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = XVELOCITY2_2D(DdataSolution,IDX2_REVERSE,ieq,_,_)
          dpre = PRESSURE2_2D(DdataSolution,IDX2_REVERSE,ieq,_,_)

          ! Overwrite the geometric source term
          DdataSource(ieq,1) = daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_2D(DdataSolution,IDX2_REVERSE,ieq,_,_)
          DdataSource(ieq,2) = daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_2D(DdataSolution,IDX2_REVERSE,ieq,_,_) * dvel
          DdataSource(ieq,3) = daux * DdataMassMatrix(ieq) *&
                               YMOMENTUM2_2D(DdataSolution,IDX2_REVERSE,ieq,_,_) * dvel
          DdataSource(ieq,4) = daux * DdataMassMatrix(ieq) *&
                               (TOTALENERGY2_2D(DdataSolution,IDX2_REVERSE,ieq,_,_)+dpre)*dvel
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius,dpre,dvel)&
        !$omp if(neq > HYDRO_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(ndim*(ieq-1)+1)
          dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Compute the radial velocity and pressure
          dvel = XVELOCITY2_2D(DdataSolution,IDX2_REVERSE,ieq,_,_)
          dpre = PRESSURE2_2D(DdataSolution,IDX2_REVERSE,ieq,_,_)

          ! Update the geometric source term
          DdataSource(ieq,1) = DdataSource(ieq,1) + daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_2D(DdataSolution,IDX2_REVERSE,ieq,_,_)
          DdataSource(ieq,2) = DdataSource(ieq,2) + daux * DdataMassMatrix(ieq) *&
                               XMOMENTUM2_2D(DdataSolution,IDX2_REVERSE,ieq,_,_) * dvel
          DdataSource(ieq,3) = DdataSource(ieq,3) + daux * DdataMassMatrix(ieq) *&
                               YMOMENTUM2_2D(DdataSolution,IDX2_REVERSE,ieq,_,_) * dvel
          DdataSource(ieq,4) = DdataSource(ieq,4) + daux * DdataMassMatrix(ieq) *&
                               (TOTALENERGY2_2D(DdataSolution,IDX2_REVERSE,ieq,_,_)+dpre)*dvel
        end do
        !$omp end parallel do

      end if

    end subroutine doSource2DBlockLumped
    
    !**************************************************************
    ! Calculate the geometric source term for axi-symmetric (dalpha=1)
    ! flow in 2D. This routine assembles the geometric source term for
    ! systems stored in interleaved format using the consistent mass matrix.
    
    subroutine doSource2DIntlConsistent(deffectiveScale,&
        deffectiveRadius, ndim, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, Kld, Kcol, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Spatial dimension
      integer, intent(in) :: ndim

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar

      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:), intent(in) :: Dcoords

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
            daux = Dcoords(ndim*(jeq-1)+1)
            dradius = max(abs(daux), deffectiveRadius)
            
            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
            
            ! Compute the radial velocity and pressure
            dvel = XVELOCITY2_2D(DdataSolution,IDX2_FORWARD,jeq,_,_)
            dpre = PRESSURE2_2D(DdataSolution,IDX2_FORWARD,jeq,_,_)

            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_2D(DdataSolution,IDX2_FORWARD,jeq,_,_)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_2D(DdataSolution,IDX2_FORWARD,jeq,_,_) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       YMOMENTUM2_2D(DdataSolution,IDX2_FORWARD,jeq,_,_) * dvel
            Ddata(4) = Ddata(4) + daux * DdataMassMatrix(ia) *&
                       (TOTALENERGY2_2D(DdataSolution,IDX2_FORWARD,jeq,_,_)+dpre)*dvel
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
            daux = Dcoords(ndim*(jeq-1)+1)
            dradius = max(abs(daux), deffectiveRadius)

            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
           
            ! Compute the radial velocity and pressure
            dvel = XVELOCITY2_2D(DdataSolution,IDX2_FORWARD,jeq,_,_)
            dpre = PRESSURE2_2D(DdataSolution,IDX2_FORWARD,jeq,_,_)
 
            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_2D(DdataSolution,IDX2_FORWARD,jeq,_,_)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_2D(DdataSolution,IDX2_FORWARD,jeq,_,_) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       YMOMENTUM2_2D(DdataSolution,IDX2_FORWARD,jeq,_,_) * dvel
            Ddata(4) = Ddata(4) + daux * DdataMassMatrix(ia) *&
                       (TOTALENERGY2_2D(DdataSolution,IDX2_FORWARD,jeq,_,_)+dpre)*dvel
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
        deffectiveRadius, ndim, neq, nvar, bclear, Dcoords,&
        DdataMassMatrix, Kld, Kcol, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Spatial dimension
      integer, intent(in) :: ndim

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq
      
      ! Number of variables
      integer, intent(in) :: nvar

      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:), intent(in) :: Dcoords

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
            daux = Dcoords(ndim*(jeq-1)+1)
            dradius = max(abs(daux), deffectiveRadius)
            
            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
            
            ! Compute the radial velocity and pressure
            dvel = XVELOCITY2_2D(DdataSolution,IDX2_REVERSE,jeq,_,_)
            dpre = PRESSURE2_2D(DdataSolution,IDX2_REVERSE,jeq,_,_)

            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_2D(DdataSolution,IDX2_REVERSE,jeq,_,_)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_2D(DdataSolution,IDX2_REVERSE,jeq,_,_) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       YMOMENTUM2_2D(DdataSolution,IDX2_REVERSE,jeq,_,_) * dvel
            Ddata(4) = Ddata(4) + daux * DdataMassMatrix(ia) *&
                       (TOTALENERGY2_2D(DdataSolution,IDX2_REVERSE,jeq,_,_)+dpre)*dvel
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
            daux = Dcoords(ndim*(jeq-1)+1)
            dradius = max(abs(daux), deffectiveRadius)

            ! Compute unit vector into the origin, scale it be the user-
            ! defined scaling parameter dscale and devide it by the radius
            daux = -sign(1.0_DP, daux) * deffectiveScale / dradius
           
            ! Compute the radial velocity and pressure
            dvel = XVELOCITY2_2D(DdataSolution,IDX2_REVERSE,jeq,_,_)
            dpre = PRESSURE2_2D(DdataSolution,IDX2_REVERSE,jeq,_,_)
 
            ! Update the geometric source term
            Ddata(1) = Ddata(1) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_2D(DdataSolution,IDX2_REVERSE,jeq,_,_)
            Ddata(2) = Ddata(2) + daux * DdataMassMatrix(ia) *&
                       XMOMENTUM2_2D(DdataSolution,IDX2_REVERSE,jeq,_,_) * dvel
            Ddata(3) = Ddata(3) + daux * DdataMassMatrix(ia) *&
                       YMOMENTUM2_2D(DdataSolution,IDX2_REVERSE,jeq,_,_) * dvel
            Ddata(4) = Ddata(4) + daux * DdataMassMatrix(ia) *&
                       (TOTALENERGY2_2D(DdataSolution,IDX2_REVERSE,jeq,_,_)+dpre)*dvel
          end do
          
          ! Update the geometric source term
          DdataSource(ieq,:) = DdataSource(ieq,:) + Ddata
        end do
        !$omp end parallel do

      end if
      
    end subroutine doSource2DBlockConsistent

  end subroutine hydro_calcGeometricSourceterm

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcDivergenceVector(rproblemLevel, rboundaryCondition,&
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
              hydro_calcFluxGalerkin1d_sim, dscale, bclear, rvector, rcollection)
          
        case (NDIM2D)
          call gfsys_buildVectorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution,&
              hydro_calcFluxGalerkin2d_sim, dscale, bclear, rvector, rcollection)
          
        case (NDIM3D)
          call gfsys_buildVectorEdge(&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution,&
              hydro_calcFluxGalerkin3d_sim, dscale, bclear, rvector, rcollection)
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
                hydro_calcFluxGalerkin1d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxGalerkin2d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecGalerkin2d_cuda))
            
          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxGalerkin3d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecGalerkin2d_cuda))
          end select
          
          !---------------------------------------------------------------------
          
        case (DISSIPATION_SCALAR)
          
          ! Assemble divergence of flux with scalar dissipation
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxScDiss1d_sim, dscale, bclear , rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxScDiss2d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecScDiss2d_cuda))

          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxScDiss3d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecScDiss3d_cuda))
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
                hydro_calcFluxScDiss1d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxScDissDiSp2d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecScDissDiSp2d_cuda))

          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxScDissDiSp3d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecScDissDiSp3d_cuda))
          end select
          
          !---------------------------------------------------------------------

        case (DISSIPATION_ROE)
          
          ! Assemble divergence of flux with Roe-type dissipation
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxRoeDiss1d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxRoeDiss2d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecRoeDiss2d_cuda))

          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxRoeDiss3d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecRoeDiss3d_cuda))
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
                hydro_calcFluxRoeDiss1d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxRoeDissDiSp2d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecRoeDissDiSp2d_cuda))

          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxRoeDissDiSp3d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecRoeDissDiSp3d_cuda))
          end select
          
          !---------------------------------------------------------------------

        case (DISSIPATION_RUSANOV)

          ! Assemble divergence of flux with Rusanov-type flux
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxRusDiss1d_sim, dscale, bclear, rvector, rcollection)

          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxRusDiss2d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecRusDiss2d_cuda))

          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxRusDiss3d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecRusDiss3d_cuda))
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
                hydro_calcFluxRusDiss1d_sim, dscale, bclear, rvector, rcollection)
            
          case (NDIM2D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxRusDissDiSp2d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecRusDissDiSp2d_cuda))

          case (NDIM3D)
            call gfsys_buildVectorEdge(&
                rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
                rsolution,&
                hydro_calcFluxRusDissDiSp3d_sim, dscale, bclear, rvector, rcollection&
                COPROC_FCB_CALCVECTOREDGESYS(hydro_calcDivVecRusDissDiSp3d_cuda))
          end select
          
        case default
          call output_line('Invalid type of dissipation!',&
              OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivergenceVector')
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
              hydro_calcFluxGalNoBdr1d_sim,&
              hydro_calcCharacteristics1d_sim, dscale, bclear, rvector, rcollection)
          
        case (NDIM2D)
          call afcsys_buildVectorTVD(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, NDIM2D,&
              hydro_calcFluxGalNoBdr2d_sim,&
              hydro_calcCharacteristics2d_sim, dscale, bclear, rvector, rcollection)
          
        case (NDIM3D)
          call afcsys_buildVectorTVD(&
              rproblemLevel%Rafcstab(inviscidAFC),&
              rproblemLevel%RgroupFEMBlock(inviscidGFEM)%RgroupFEMBlock(1),&
              rsolution, NDIM3D,&
              hydro_calcFluxGalNoBdr3d_sim,&
              hydro_calcCharacteristics3d_sim, dscale, bclear, rvector, rcollection)
        end select

        !-----------------------------------------------------------------------

      case default
        call output_line('Invalid type of stabilisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcDivergenceVector')
        call sys_halt()
      end select
      
      !-------------------------------------------------------------------------
      ! Evaluate linear form for boundary integral (if any)
      !-------------------------------------------------------------------------

      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        call hydro_calcLinfBdrCond1D(rproblemLevel, rboundaryCondition,&
            rsolution, dtime, -dscale, ssectionName, hydro_coeffVectorBdr1d_sim,&
            rvector, rcollection)
        
      case (NDIM2D)
        call hydro_calcLinfBdrCond2D(rproblemLevel, rboundaryCondition,&
            rsolution, dtime, -dscale, ssectionName, hydro_coeffVectorBdr2d_sim,&
            rvector, rcollection)
        
      case (NDIM3D)
!!$        call hydro_calcLinfBdrCond3D(rproblemLevel, rboundaryCondition,&
!!$            rsolution, dtime, -dscale, ssectionName, hydro_coeffVectorBdr3d_sim,&
!!$            rvector, rcollection)
        print *, "Boundary conditions in 3D have not been implemented yet!"
        stop
      end select

    end if

  end subroutine hydro_calcDivergenceVector

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcTimeDerivative(rproblemLevel, rtimestep,&
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
      call hydro_calcDivergenceVector(rproblemLevel,&
          rsolver%rboundaryCondition, rsolution, rtimestep%dTime,&
          1.0_DP, .true., p_rvector1, ssectionName, rcollection)
      
      ! Build the geometric source term (if any)
      call hydro_calcGeometricSourceterm(p_rparlist, ssectionName,&
          rproblemLevel, rsolution, 1.0_DP, .false., p_rvector1, rcollection)
      
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
      call lsysbl_invertedDiagMatVec(rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
          p_rvector1, 1.0_DP, rvector)

      ! Store norm of the initial guess from the lumped version
      dnorm0 = lsysbl_vectorNorm(rvector1, LINALG_NORML2)
      
      richardson: do ite = 1, nmaxIterationsApproxTimeDerivative
        ! Initialise rvector2 by the constant right-hand side
        call lsysbl_copyVector(p_rvector1, p_rvector2)
        
        ! Compute the residual $rhs-M_C*u$ and store the result in rvector3
        do iblock = 1,rsolution%nblocks
          call lsyssc_matVec(rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rvector%RvectorBlock(iblock), p_rvector2%RvectorBlock(iblock),&
              -1.0_DP, 1.0_DP)
        end do
          
        ! Scale rvector2 by the inverse of the lumped mass matrix
        call lsysbl_invertedDiagMatVec(rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
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
      call hydro_calcDivergenceVector(rproblemLevel,&
          rsolver%rboundaryCondition, rsolution, rtimestep%dTime,&
          1.0_DP, .true., rvector, ssectionName, rcollection)
      
      ! Build the geometric source term (if any)
      call hydro_calcGeometricSourceterm(p_rparlist, ssectionName,&
          rproblemLevel, rsolution, 1.0_DP, .false., rvector, rcollection)
      
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
      call lsysbl_invertedDiagMatVec(rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
          rvector, 1.0_DP, rvector)

    case default
      call output_line('Unsupported type of divergence term!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcTimeDerivative')
      call sys_halt()
    end select

  end subroutine hydro_calcTimeDerivative

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcCFLnumber(rproblemLevel, rtimestep, rsolution,&
      ssectionName, rcollection, dCFLnumber)

!<description>
    ! This subroutine calculates the local CFL-number of the given solution
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<inputoutput>
    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>

!<output>
    ! local CFL-number
    real(DP), intent(out) :: dCFLnumber
!</output>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    real(DP), dimension(:), pointer :: p_Dx,p_Ddata
    integer :: isystemFormat,lumpedMassMatrix

    ! Get parameters from parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'isystemformat', isystemFormat)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)

    ! Set pointer to lumped mass matrix
    call lsyssc_getbase_double(&
        rproblemLevel%RmatrixScalar(lumpedMassMatrix), p_Dx)

    ! Set pointer to solution vector
    call lsysbl_getbase_double(rsolution, p_Ddata)

    select case(rproblemLevel%rtriangulation%ndim)
    case (NDIM1D)
      select case(isystemFormat)
      case (SYSTEM_INTERLEAVEFORMAT)
        call calcCFDnumber1DIntl(rsolution%NEQ, rtimestep%dStep,&
            p_Ddata, p_Dx, dCFLnumber)
      case (SYSTEM_BLOCKFORMAT)
        call calcCFDnumber1DBlock(rsolution%NEQ, rtimestep%dStep,&
            p_Ddata, p_Dx, dCFLnumber)
      case default
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcCFLnumber')
        call sys_halt()
      end select
      
    case (NDIM2D)
      select case(isystemFormat)
      case (SYSTEM_INTERLEAVEFORMAT)
        call calcCFDnumber2DIntl(rsolution%NEQ, rtimestep%dStep,&
            p_Ddata, p_Dx, dCFLnumber)
      case (SYSTEM_BLOCKFORMAT)
        call calcCFDnumber2DBlock(rsolution%NEQ, rtimestep%dStep,&
            p_Ddata, p_Dx, dCFLnumber)
      case default
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcCFLnumber')
        call sys_halt()
      end select
      
    case (NDIM3D)
      select case(isystemFormat)
      case (SYSTEM_INTERLEAVEFORMAT)
        call calcCFDnumber3DIntl(rsolution%NEQ, rtimestep%dStep,&
            p_Ddata, p_Dx, dCFLnumber)
      case (SYSTEM_BLOCKFORMAT)
        call calcCFDnumber3DBlock(rsolution%NEQ, rtimestep%dStep,&
            p_Ddata, p_Dx, dCFLnumber)
      case default
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcCFLnumber')
        call sys_halt()
      end select
      
    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'hydro_calcCFLnumber')
      call sys_halt()
    end select

  contains

    ! Here, the real working routines follow

    !**************************************************************
    ! Calculates the CFL-number for the one-dimensional solution
    ! vector stored in interleaved format

    subroutine calcCFDnumber1DIntl(neq, dStep, Ddata, Dx, dCFLnumber)

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq

      ! Global time step
      real(DP), intent(in) :: dStep

      ! Solution data
      real(DP), dimension(NVAR1D,neq), intent(in) :: Ddata

      ! Lumped mass matrix
      real(DP), dimension(neq), intent(in) :: Dx

      ! Local CFD number
      real(DP), intent(out) :: dCFLnumber

      ! local variables
      real(DP) :: dlambda
      integer :: ieq

      dCFLnumber = 0.0_DP
      do ieq=1,neq
        dlambda = VELMAGNITUDE2_1D(Ddata,IDX2_FORWARD,ieq,_,_)+&
                    SOUNDSPEED2_1D(Ddata,IDX2_FORWARD,ieq,_,_)
        dCFLnumber = max(dCFLnumber, dStep*dlambda/sqrt(Dx(ieq)))
      end do
      
    end subroutine calcCFDnumber1DIntl

    !**************************************************************
    ! Calculates the CFL-number for the one-dimensional solution
    ! vector stored in block format

    subroutine calcCFDnumber1DBlock(neq, dStep, Ddata, Dx, dCFLnumber)

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq

      ! Global time step
      real(DP), intent(in) :: dStep

      ! Solution data
      real(DP), dimension(neq,NVAR1D), intent(in) :: Ddata

      ! Lumped mass matrix
      real(DP), dimension(neq), intent(in) :: Dx

      ! Local CFD number
      real(DP), intent(out) :: dCFLnumber

      ! local variables
      real(DP) :: dlambda
      integer :: ieq

      dCFLnumber = 0.0_DP
      do ieq=1,neq
        dlambda = VELMAGNITUDE2_1D(Ddata,IDX2_REVERSE,ieq,_,_)+&
                    SOUNDSPEED2_1D(Ddata,IDX2_REVERSE,ieq,_,_)
        dCFLnumber = max(dCFLnumber, dStep*dlambda/sqrt(Dx(ieq)))
      end do
      
    end subroutine calcCFDnumber1DBlock

    !**************************************************************
    ! Calculates the CFL-number for the two-dimensional solution
    ! vector stored in interleaved format

    subroutine calcCFDnumber2DIntl(neq, dStep, Ddata, Dx, dCFLnumber)

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq

      ! Global time step
      real(DP), intent(in) :: dStep

      ! Solution data
      real(DP), dimension(NVAR2D,neq), intent(in) :: Ddata

      ! Lumped mass matrix
      real(DP), dimension(neq), intent(in) :: Dx

      ! Local CFD number
      real(DP), intent(out) :: dCFLnumber

      ! local variables
      real(DP) :: dlambda
      integer :: ieq

      dCFLnumber = 0.0_DP
      do ieq=1,neq
        dlambda = VELMAGNITUDE2_2D(Ddata,IDX2_FORWARD,ieq,_,_)+&
                    SOUNDSPEED2_2D(Ddata,IDX2_FORWARD,ieq,_,_)
        dCFLnumber = max(dCFLnumber, dStep*dlambda/sqrt(Dx(ieq)))
      end do
      
    end subroutine calcCFDnumber2DIntl

    !**************************************************************
    ! Calculates the CFL-number for the two-dimensional solution
    ! vector stored in block format

    subroutine calcCFDnumber2DBlock(neq, dStep, Ddata, Dx, dCFLnumber)

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq

      ! Global time step
      real(DP), intent(in) :: dStep

      ! Solution data
      real(DP), dimension(neq,NVAR2D), intent(in) :: Ddata

      ! Lumped mass matrix
      real(DP), dimension(neq), intent(in) :: Dx

      ! Local CFD number
      real(DP), intent(out) :: dCFLnumber

      ! local variables
      real(DP) :: dlambda
      integer :: ieq

      dCFLnumber = 0.0_DP
      do ieq=1,neq
        dlambda = VELMAGNITUDE2_2D(Ddata,IDX2_REVERSE,ieq,_,_)+&
                    SOUNDSPEED2_2D(Ddata,IDX2_REVERSE,ieq,_,_)
        dCFLnumber = max(dCFLnumber, dStep*dlambda/sqrt(Dx(ieq)))
      end do
      
    end subroutine calcCFDnumber2DBlock

    !**************************************************************
    ! Calculates the CFL-number for the three-dimensional solution
    ! vector stored in interleaved format

    subroutine calcCFDnumber3DIntl(neq, dStep, Ddata, Dx, dCFLnumber)

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq

      ! Global time step
      real(DP), intent(in) :: dStep

      ! Solution data
      real(DP), dimension(NVAR3D,neq), intent(in) :: Ddata

      ! Lumped mass matrix
      real(DP), dimension(neq), intent(in) :: Dx

      ! Local CFD number
      real(DP), intent(out) :: dCFLnumber

      ! local variables
      real(DP) :: dlambda
      integer :: ieq

      dCFLnumber = 0.0_DP
      do ieq=1,neq
        dlambda = VELMAGNITUDE2_3D(Ddata,IDX2_FORWARD,ieq,_,_)+&
                    SOUNDSPEED2_3D(Ddata,IDX2_FORWARD,ieq,_,_)
        dCFLnumber = max(dCFLnumber, dStep*dlambda/sqrt(Dx(ieq)))
      end do
      
    end subroutine calcCFDnumber3DIntl

    !**************************************************************
    ! Calculates the CFL-number for the three-dimensional solution
    ! vector stored in block format

    subroutine calcCFDnumber3DBlock(neq, dStep, Ddata, Dx, dCFLnumber)

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq

      ! Global time step
      real(DP), intent(in) :: dStep

      ! Solution data
      real(DP), dimension(neq,NVAR3D), intent(in) :: Ddata

      ! Lumped mass matrix
      real(DP), dimension(neq), intent(in) :: Dx

      ! Local CFD number
      real(DP), intent(out) :: dCFLnumber

      ! local variables
      real(DP) :: dlambda
      integer :: ieq

      dCFLnumber = 0.0_DP
      do ieq=1,neq
        dlambda = VELMAGNITUDE2_3D(Ddata,IDX2_REVERSE,ieq,_,_)+&
                    SOUNDSPEED2_3D(Ddata,IDX2_REVERSE,ieq,_,_)
        dCFLnumber = max(dCFLnumber, dStep*dlambda/sqrt(Dx(ieq)))
      end do
      
    end subroutine calcCFDnumber3DBlock

  end subroutine hydro_calcCFLnumber

end module hydro_callback

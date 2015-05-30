!##############################################################################
!# ****************************************************************************
!# <Name> transport_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve scalar conservation laws in arbitrary spatial dimensions.
!#
!# The following callback functions are available:
!#
!# 1.) transp_nlsolverCallback
!#     -> Callback routine for the nonlinear solver
!#
!# ****************************************************************************
!#
!# The following auxiliary routines are available:
!#
!# 1.) transp_calcPreconditioner
!#     -> Calculates the nonlinear preconditioner
!#
!# 2.) transp_calcJacobian
!#     -> Calculates the Jacobian matrix
!#
!# 3.) transp_calcResidual
!#     -> Calculates the nonlinear residual vector
!#
!# 4.) transp_calcRhs
!#     -> Calculates the explicit right-hand side vector
!#
!# 5.) transp_calcVelocityField
!#     -> Calculates the velocity field
!#
!# 6.) transp_calcLinearisedFCT
!#     -> Calculates the linearised FCT correction
!#
!# 7.) transp_calcLinearisedLPT
!#     -> Calculates the linearised linearity-preserving correction
!#
!# 8.) transp_calcBoundsAtEdgesLPT
!#     -> Calculates the bounds at edges for linearity-preserving
!#        flux correction
!#
!# 9.) transp_calcGeometricSourceterm
!#     -> Calculates the geometric source term for axi-symmetric,
!#        cylindrical or sperical symmetric coordinate systems.
!#
!# 10.) transp_calcTransportOperator
!#      -> Calculates the discrete transport operator.
!#
!# 11.) transp_calcTimeDerivative
!#      -> Cacluates the approximate time derivative
!#
!# 12.) transp_coeffVectorFE
!#      -> Callback routine for the evaluation of linear forms
!#         using a given FE-solution for interpolation
!#
!# 13.) transp_coeffVectorAnalytic
!#      -> Callback routine for the evaluation of linear forms
!#         using an analytic expression for the load-vector
!#
!# 14.) transp_refFuncAnalytic
!#      -> Callback routine for the evaluation of the reference
!#         target function for goal-oriented error estimation
!#
!# 15.) transp_weightFuncAnalytic
!#      -> Callback routine for the evaluation of the weights in
!#         the target functional for goal-oriented error estimation
!#
!# 16.) transp_weightFuncGHill
!#      -> Callback routine for the evaluation of the weights in
!#         the rotation of a Gaussian hill benchmark
!#
!# 17.) transp_parseBoundaryCondition
!#      -> Callback routine for the treatment of boundary conditions
!#
!# 18.) transp_setBoundaryCondition
!#      -> Imposes boundary conditions for nonlinear solver
!#         by filtering the system matrix and the solution/residual
!#         vector explicitly (i.e. strong boundary conditions)
!#
!# 19.) transp_calcBilfBdrCondition
!#      -> Wrapper routine for the calculation of the bilinear form
!#         arising from the weak imposition of boundary conditions;
!#         this routine calls transp_calcBilfBdrCondXd with the
!#         correct callback routines depending on the type of velocity
!#         and the mode, i.e. primal or dual
!#
!# 20.) transp_calcLinfBdrCondition
!#      -> Wrapper routine for the calculation of the linear form
!#         arising from the weak imposition of boundary conditions;
!#         this routine calls transp_calcLinfBdrCondXd with the
!#         correct callback routines depending on the type of velocity
!#         and the mode, i.e. primal or dual
!#
!# 21.) transp_calcOperatorBdrCondition
!#      -> Calculates the discrete operator at the boundary by the
!#         group finite element formulation
!#
!# 22.) transp_calcVectorBdrCondition
!#      -> Calculates the discrete vector at the boundary by the
!#         group finite element formulation
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
!#        transp_calcResidual. Note that changing the
!#        interface of this callback routine would require to update
!#        ALL models. Hence, the interface should only be changed if
!#        really necessary.
!#
!# 2.) Where do I have to implement 'what' if I need another type of velocity?
!#
!#     -> In essence, you have to add some code in each subroutine
!#        which does a SELECT CASE on IVELOCITYTYPE. To be more
!#        precise you should search for the tag @FAQ2: step-by-step
!#        and create a new CASE that corresponds to the new type.
!# </purpose>
!##############################################################################

module transport_callback

#include "flagship.h"

!$ use omp_lib
  use afcstabbase
  use afcstabscalar
  use afcstabscalarfct
  use afcstabscalargp
  use afcstabscalarlpt
  use afcstabscalarsymm
  use afcstabscalartvd
  use basicgeometry
  use bilinearformevaluation
  use boundary
  use boundarycondaux
  use boundaryfilter
  use collection
  use cubature
  use derivatives
  use dofmapping
  use domainintegration
  use feevaluation
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use groupfembase
  use groupfemscalar
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

  ! Modules from transport model
  use transport_basic
  use transport_callback1d
  use transport_callback2d
  use transport_callback3d

  implicit none

  private
  public :: transp_nlsolverCallback
  public :: transp_calcPreconditioner
  public :: transp_calcJacobian
  public :: transp_calcResidual
  public :: transp_calcRhs
  public :: transp_setBoundaryCondition
  public :: transp_calcVelocityField
  public :: transp_calcLinearisedFCT
  public :: transp_calcLinearisedLPT
  public :: transp_calcGeometricSourceterm
  public :: transp_calcTransportOperator
  public :: transp_calcTimeDerivative
  public :: transp_coeffVectorFE
  public :: transp_coeffVectorAnalytic
  public :: transp_refFuncAnalytic
  public :: transp_refFuncAnalyticBdr2D
  public :: transp_weightFuncAnalytic
  public :: transp_weightFuncGHill
  public :: transp_parseBoundaryCondition
  public :: transp_calcOperatorBdrCondition
  public :: transp_calcVectorBdrCondition

  !*****************************************************************************

!<constants>

!<constantblock>
  ! Minimum number of equations for OpenMP parallelisation: If the number of
  ! equations is below this value, then no parallelisation is performed.
#ifndef TRANSP_GEOMSOURCE_NEQMIN_OMP
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: TRANSP_GEOMSOURCE_NEQMIN_OMP = 10000
#else
  integer, public            :: TRANSP_GEOMSOURCE_NEQMIN_OMP = 10000
#endif
#endif
  
!</constantblock>

!</constants>

contains

  !*****************************************************************************

!<subroutine>

  subroutine transp_nlsolverCallback(rproblemLevel, rtimestep,&
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
    type(t_parlist), pointer :: p_rparlist
    character(len=SYS_STRLEN) :: ssectionName
    integer(i32) :: iSpec
    integer :: jacobianMatrix


    ! Get section name
    call collct_getvalue_string(rcollection,&
        'ssectionname', ssectionName)

    !###########################################################################
    ! REMARK: The order in which the operations are performed is
    ! essential. This is due to the fact that the calculation of the
    ! residual/rhs requires the discrete transport operator to be
    ! initialised which is assembled in the calculation of the
    ! preconditioner. To prevent the re-assembly of the
    ! preconditioner twice, we remove the specifier
    ! NLSOL_OPSPEC_CALCPRECOND if the residual/rhs vector is built.
    !###########################################################################

    ! Make a local copy
    iSpec = ioperationSpec

    
    ! Do we have to calculate the residual?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      if (istep .eq. 0) then
        ! If the preconditioner and/or the transport operator have not
        ! been initialised (e.g., in the very first solution step or
        ! if the mesh has been regenerated) they have to be generated
        ! explicitly for the old solution from the previous time step.
        if ((iand(rproblemLevel%iproblemSpec, TRANSP_PRECOND_INIT) .ne. 0) .or.&
            (iand(rproblemLevel%iproblemSpec, TRANSP_TROPER_INIT)  .ne. 0))&
            call transp_calcPreconditioner(rproblemLevel, rtimestep,&
            rsolver, rsolution0, ssectionName, rcollection)

        ! Compute the constant right-hand side only in the very first
        ! step of each nonlinear solution loop
        call transp_calcRhs(rproblemLevel, rtimestep,&
            rsolver, rsolution0, rrhs, ssectionName, rcollection, rsource)
      end if

      ! Compute the preconditioner
      call transp_calcPreconditioner(rproblemLevel, rtimestep,&
          rsolver, rsolution, ssectionName, rcollection)

      ! Compute the residual
      call transp_calcResidual(rproblemLevel, rtimestep, rsolver,&
          rsolution, rsolution0, rrhs, rres, istep, ssectionName, rcollection)

      ! Remove specifier for the preconditioner (if any)
      iSpec = iand(iSpec, not(NLSOL_OPSPEC_CALCPRECOND))
    end if


    ! Do we have to calculate the preconditioner?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) then

      ! Compute the preconditioner
      call transp_calcPreconditioner(rproblemLevel, rtimestep,&
          rsolver, rsolution, ssectionName, rcollection)
    end if


    ! Do we have to calculate the Jacobian operator?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_CALCJACOBIAN) .ne. 0) then

      ! Compute the Jacobian matrix
      call transp_calcJacobian(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, ssectionName, rcollection)
    end if


    ! Do we have to impose boundary conditions?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      ! Impose boundary conditions
      call transp_setBoundaryCondition(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rres, ssectionName, rcollection)
    end if


    ! Do we have to apply the Jacobian operator?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_APPLYJACOBIAN) .ne. 0) then

      p_rparlist => collct_getvalue_parlst(rcollection,&
          'rparlist', ssectionName=ssectionName)
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'jacobianMatrix', jacobianMatrix)

      ! Apply Jacobian matrix
      call lsyssc_matVec(rproblemLevel%RmatrixScalar(jacobianMatrix),&
          rsolution%RvectorBlock(1), rres%RvectorBlock(1), 1.0_DP, 1.0_DP)
    end if


    ! Set status flag
    istatus = 0

  end subroutine transp_nlsolverCallback

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcPreconditioner(rproblemLevel,&
      rtimestep, rsolver, rsolution, ssectionName, rcollection,&
      fcb_calcMatrixDiagPrimal_sim, fcb_calcMatrixPrimal_sim,&
      fcb_calcMatrixDiagDual_sim, fcb_calcMatrixDual_sim,&
      fcb_coeffMatBdrPrimal1d_sim, fcb_coeffMatBdrDual1d_sim,&
      fcb_coeffMatBdrPrimal2d_sim, fcb_coeffMatBdrDual2d_sim,&
      fcb_coeffMatBdrPrimal3d_sim, fcb_coeffMatBdrDual3d_sim,&
      bforceUpdate)

!<description>
    ! This subroutine calculates the nonlinear preconditioner for the
    ! primal/dual problem and configures the linear solver structure
    ! accordingly. The details of this subroutine are explained in code.
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCalcMatrix.inc'
    optional :: fcb_calcMatrixDiagPrimal_sim
    optional :: fcb_calcMatrixPrimal_sim
    optional :: fcb_calcMatrixDiagDual_sim
    optional :: fcb_calcMatrixDual_sim

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCoeffMatBdr.inc'
    optional :: fcb_coeffMatBdrPrimal1d_sim
    optional :: fcb_coeffMatBdrPrimal2d_sim
    optional :: fcb_coeffMatBdrPrimal3d_sim
    optional :: fcb_coeffMatBdrDual1d_sim
    optional :: fcb_coeffMatBdrDual2d_sim
    optional :: fcb_coeffMatBdrDual3d_sim

    ! OPTIONAL: flag to overwrite update notifier
    ! If true, then the preconditioner is updated even if the
    ! update notifier does not indicate the need of an update
    logical, intent(in), optional :: bforceUpdate
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    integer :: lumpedMassMatrix, consistentMassMatrix
    integer :: systemMatrix, transportMatrix, imasstype
    logical :: breturn


    breturn = .true.
    if (present(bforceUpdate)) breturn = .not.(bforceUpdate)

    ! Check if the preconditioner and/or the discrete transport
    ! operator have to be updated and return otherwise.
    if ((iand(rproblemLevel%iproblemSpec, TRANSP_PRECOND_INIT)   .eq. 0) .and.&
        (iand(rproblemLevel%iproblemSpec, TRANSP_TROPER_INIT)    .eq. 0) .and.&
        (iand(rproblemLevel%iproblemSpec, TRANSP_PRECOND_UPDATE) .eq. 0) .and.&
        (iand(rproblemLevel%iproblemSpec, TRANSP_TROPER_UPDATE)  .eq. 0) .and.&
        breturn) return

    ! Start time measurement for matrix evaluation
    p_rtimer => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyMatrix', ssectionName=ssectionName)
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Remove init/update notifier from the preconditioner. Depending
    ! on the velocity, diffusion type it will be re-activited below.
    rproblemLevel%iproblemSpec = iand(rproblemLevel%iproblemSpec,&
        not(TRANSP_PRECOND_INIT+TRANSP_PRECOND_UPDATE))

    ! Set pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)

    ! Get parameter from parameter list
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'imasstype', imasstype)
    
    ! What type of mass matrix are we?
    select case(imasstype)
    case (MASS_LUMPED)
      
      !-------------------------------------------------------------------------
      ! Compute the global operator for transient flow
      !
      !   $ A = M_L - scaleImplicit * dt * L $
      !-------------------------------------------------------------------------

      ! Get additional parameters from parameter list
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'systemmatrix', systemMatrix)
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'transportmatrix', transportMatrix)

      ! Assemble the discrete transport operator
      call transp_calcTransportOperator(rproblemlevel,&
        rsolver%rboundaryCondition, rsolution, rtimestep%dTime, 1.0_DP,&
        rproblemLevel%RmatrixScalar(transportMatrix), ssectionName, rcollection,&
        fcb_calcMatrixDiagPrimal_sim, fcb_calcMatrixPrimal_sim,&
        fcb_calcMatrixDiagDual_sim,  fcb_calcMatrixDual_sim,&
        fcb_coeffMatBdrPrimal1d_sim, fcb_coeffMatBdrDual1d_sim,&
        fcb_coeffMatBdrPrimal2d_sim, fcb_coeffMatBdrDual2d_sim,&
        fcb_coeffMatBdrPrimal3d_sim, fcb_coeffMatBdrDual3d_sim,&
        bforceUpdate)
     
      ! Build global system operator
      call lsyssc_MatrixLinearComb(&
          rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
          rproblemLevel%RmatrixScalar(transportMatrix),&
          1.0_DP, -rtimestep%dscaleImplicit*rtimestep%dStep,&
          .false., .false., .true., .true.,&
          rproblemLevel%RmatrixScalar(systemMatrix))

    case (MASS_CONSISTENT)

      !-------------------------------------------------------------------------
      ! Compute the global operator for transient flow
      !
      !   $ A = M_C - scaleImplicit * dt * L $
      !-------------------------------------------------------------------------

      ! Get additional parameters from parameter list
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'consistentmassmatrix', consistentMassMatrix)
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'systemmatrix', systemMatrix)
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'transportmatrix', transportMatrix)

      ! Assemble the discrete transport operator
      call transp_calcTransportOperator(rproblemlevel,&
        rsolver%rboundaryCondition, rsolution, rtimestep%dTime, 1.0_DP,&
        rproblemLevel%RmatrixScalar(transportMatrix), ssectionName, rcollection,&
        fcb_calcMatrixDiagPrimal_sim, fcb_calcMatrixPrimal_sim,&
        fcb_calcMatrixDiagDual_sim,  fcb_calcMatrixDual_sim,&
        fcb_coeffMatBdrPrimal1d_sim, fcb_coeffMatBdrDual1d_sim,&
        fcb_coeffMatBdrPrimal2d_sim, fcb_coeffMatBdrDual2d_sim,&
        fcb_coeffMatBdrPrimal3d_sim, fcb_coeffMatBdrDual3d_sim,&
        bforceUpdate)

      ! Build global system operator
      call lsyssc_MatrixLinearComb(&
          rproblemLevel%RmatrixScalar(consistentMassMatrix),&
          rproblemLevel%RmatrixScalar(transportMatrix),&
          1.0_DP, -rtimestep%dscaleImplicit*rtimestep%dStep,&
          .false., .false., .true., .true.,&
          rproblemLevel%RmatrixScalar(systemMatrix))

    case default

      !-------------------------------------------------------------------------
      ! Compute the global operator for steady-state flow
      !
      !   $ A = -L $
      !-------------------------------------------------------------------------

      ! Get additional parameters from parameter list
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'transportmatrix', transportMatrix)
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'systemmatrix', systemMatrix)
    
      ! Build global system operator
      call transp_calcTransportOperator(rproblemlevel,&
        rsolver%rboundaryCondition, rsolution, rtimestep%dTime, 1.0_DP,&
        rproblemLevel%RmatrixScalar(transportMatrix), ssectionName, rcollection,&
        fcb_calcMatrixDiagPrimal_sim, fcb_calcMatrixPrimal_sim,&
        fcb_calcMatrixDiagDual_sim,  fcb_calcMatrixDual_sim,&
        fcb_coeffMatBdrPrimal1d_sim, fcb_coeffMatBdrDual1d_sim,&
        fcb_coeffMatBdrPrimal2d_sim, fcb_coeffMatBdrDual2d_sim,&
        fcb_coeffMatBdrPrimal3d_sim, fcb_coeffMatBdrDual3d_sim,&
        bforceUpdate)

      call lsyssc_copyMatrix(rproblemLevel%RmatrixScalar(transportMatrix),&
          rproblemLevel%RmatrixScalar(systemMatrix))
      call lsyssc_scaleMatrix(rproblemLevel%RmatrixScalar(systemMatrix), -1.0_DP)
      
    end select

    ! Impose boundary conditions in strong sence (if any)
    if (rsolver%rboundaryCondition%bStrongBdrCond) then
      call bdrf_filterMatrix(rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixScalar(systemMatrix), 1.0_DP)
    end if

    ! Set update notifier for the preconditioner if the update
    ! notifier for the discret transport operator has been set
    if (iand(rproblemLevel%iproblemSpec, TRANSP_TROPER_UPDATE) .ne. 0) then
      rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                       TRANSP_PRECOND_UPDATE)
    end if

    ! Ok, we updated the (nonlinear) system operator successfully. Now we still
    ! have to link it to the solver hierarchy. This is done recursively.
    call flagship_updateSolverMatrix(rproblemLevel, rsolver,&
        systemMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL,&
        rproblemlevel%ilev, rproblemLevel%ilev)

    ! Finally, we have to update the content of the solver hierarchy
    call solver_updateContent(rsolver)

    ! Stop time measurement for global operator
    call stat_stopTimer(p_rtimer)

  end subroutine transp_calcPreconditioner

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcJacobian(rproblemLevel,&
      rtimestep, rsolver, rsolution, rsolution0, ssectionName, rcollection,&
      fcb_calcMatrixDiagPrimal_sim, fcb_calcMatrixPrimal_sim,&
      fcb_calcMatrixDiagDual_sim, fcb_calcMatrixDual_sim,&
      fcb_coeffMatBdrPrimal1d_sim, fcb_coeffMatBdrDual1d_sim,&
      fcb_coeffMatBdrPrimal2d_sim, fcb_coeffMatBdrDual2d_sim,&
      fcb_coeffMatBdrPrimal3d_sim, fcb_coeffMatBdrDual3d_sim)

!<description>
    ! This callback subroutine computes the Jacobian matrix.
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolution0

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCalcMatrix.inc'
    optional :: fcb_calcMatrixDiagPrimal_sim
    optional :: fcb_calcMatrixPrimal_sim
    optional :: fcb_calcMatrixDiagDual_sim
    optional :: fcb_calcMatrixDual_sim

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCoeffMatBdr.inc'
    optional :: fcb_coeffMatBdrPrimal1d_sim
    optional :: fcb_coeffMatBdrPrimal2d_sim
    optional :: fcb_coeffMatBdrPrimal3d_sim
    optional :: fcb_coeffMatBdrDual1d_sim
    optional :: fcb_coeffMatBdrDual2d_sim
    optional :: fcb_coeffMatBdrDual3d_sim
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    type(t_collection) :: rcollectionTmp
    real(DP) :: hstep
    character(LEN=SYS_STRLEN) :: smode
    integer :: transportMatrix, jacobianMatrix
    integer :: consistentMassMatrix, lumpedMassMatrix
    integer :: convectionGFEM, coeffMatrix_S
    integer :: convectionAFC, diffusionAFC, ijacobianFormat
    integer :: imasstype, imassantidiffusiontype, ivelocitytype, idiffusiontype
    integer :: velocityfield
    logical :: bisExactStructure, bisExtendedSparsity

    !###########################################################################
    ! REMARK: If you want to add a new type of velocity/diffusion,
    ! then search for the tag @FAQ2: in this subroutine and create a
    ! new CASE which performs the corresponding task for the new type
    ! of velocity/ diffusion.
    !###########################################################################

    ! Start time measurement for matrix evaluation
    p_rtimer => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyMatrix', ssectionName=ssectionName)
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'transportmatrix', transportMatrix)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'jacobianmatrix', jacobianMatrix)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'coeffMatrix_S', coeffMatrix_S)

    ! The Jacobian matrix for the low-order transport operator needs
    ! to be generated only in case of nonlinear governing equations.
    ! In this case, the corresponding transport operator L has to be
    ! updated in each nonlinear iteration. Hence, we can temporarily
    ! build the Jacobian matrix using the memory of the operator L.

    ! Compute step lenth of the solution perturbation
    select case(int(rsolver%p_rsolverNewton%dperturbationStrategy))
    case (PERTURB_NITSOL)
      ! Choice h=[(1+|u|)*EPS]^{1/(1+p)} by Pernice et al. in
      ! M. Pernice, H.F. Walker, NITSOL: a Newton iterative solver
      ! for nonlinear systems, SIAM J. Sci. Comput. 19 (1998) 302-318.
      hstep = ( (1+lsysbl_vectorNorm(rsolution,&
                   LINALG_NORMEUCLID))*SYS_EPSREAL_DP )**(1.0_DP/3._DP)

    case (PERTURB_SQRTEPS)
      hstep= sqrt(SYS_EPSREAL_DP)

    case default
      hstep = max(SYS_EPSREAL_DP,&
                  rsolver%p_rsolverNewton%dperturbationStrategy)
    end select


    !---------------------------------------------------------------------------
    ! Assemble diffusion operator:
    !
    ! $$ \int_\Omega \nabla w \cdot (D \nabla u) {\rm d}{\bf x} $$
    !
    ! The diffusion operator is symmetric so that it is the same for
    ! the primal and the dual problem. If no diffusion is present,
    ! i.e. $D \equiv 0$, then the transport operator is initialised by
    ! zeros. If there is anisotropic diffusion, i.e. $D=D({\bf x},t)$,
    ! then we may also initialise some stabilisation structure.
    !
    ! The bilinear form for the diffusion operator consists of the
    ! volume integral (see above) only.
    !
    ! Non-homogeneous Neumann boundary conditions
    !
    ! $$ {\bf n} \cdot \nabla u = h \qquad \mbox{on} \quad \Gamma_N $$
    !
    ! are built into the right-hand side vector.
    !---------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'idiffusiontype', idiffusiontype)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'diffusionAFC', diffusionAFC)

    ! @FAQ2: What type of diffusion are we?
    select case(idiffusiontype)
    case (DIFFUSION_ZERO)
      ! zero diffusion, clear the system matrix
      call lsyssc_clearMatrix(rproblemLevel%RmatrixScalar(transportMatrix))

    case (DIFFUSION_ISOTROPIC,&
          DIFFUSION_ANISOTROPIC)
      ! Isotropic diffusion
      call lsyssc_duplicateMatrix(&
          rproblemLevel%RmatrixScalar(coeffMatrix_S),&
          rproblemLevel%RmatrixScalar(transportMatrix),&
          LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPY)

    case (DIFFUSION_VARIABLE)
      call output_line('Variable diffusion matrices are yet not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
      call sys_halt()

    case default
      call output_line('Invalid type of diffusion!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Assemble convection operator:
    !
    ! $$ \int_\Omega w \nabla \cdot ({\bf v} u) {\rm d}{\bf x} $$
    !
    ! with weakly imposed Dirichlet boundary conditions (if any):
    !
    ! $$ \int_{\Gamma_D} w (u_D-u) {\bf v} \cdot {\bf n} {\rm d}{\bf s} $$
    !
    ! The given Dirichlet boundary data $u_D$ is built into the right
    ! -hand side/residual vector whereas the surface integral for the
    ! unknown solution $u$ is added to the transport operator.
    !
    ! The convective operator is skew-symmetric so that we have to
    ! distinguish between the primal and the dual problem.
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(p_rparlist,&
        ssectionName, 'mode', smode)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'convectionAFC', convectionAFC)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'convectionGFEM', convectionGFEM, convectionAFC)

    ! Attach user-defined collection structure to temporal collection
    ! structure (may be required by the callback function)
    rcollectionTmp%p_rnextCollection => rcollection
    
    ! Set vector for velocity field (if any)
    if (transp_hasVelocityVector(ivelocityType)) then
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'velocityfield', velocityfield)
      rcollectionTmp%p_rvectorQuickAccess1 => rproblemLevel%RvectorBlock(velocityfield)
    end if

    ! Are we in primal or dual mode?
    if (trim(smode) .eq. 'primal') then

      !-------------------------------------------------------------------------
      ! We are in primal mode which means that we have to build two
      ! bilinear forms: one for the volume integral and one for the
      ! surface integral which is used to weakly impose Dirichlet
      ! boundary conditions. Note that lumping is performed for the
      ! boundary term to prevent the generation of oscillations.
      !-------------------------------------------------------------------------

      ! @FAQ2: What type of velocity are we?
      select case (abs(ivelocitytype))
      case default
        ! The user-defined callback function for matrix coefficients
        ! is used if present; otherwise an error is thrown
        if (present(fcb_calcMatrixPrimal_sim)) then

          select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
            
          case (AFCSTAB_GALERKIN)
            
            ! Apply standard Galerkin discretisation
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, fcb_calcMatrixPrimal_sim, hstep, 1.0_DP,&
                .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                fcb_calcMatrixDiagPrimal_sim, rcollection=rcollectionTmp)

          case default

            ! Apply low-order discretisation
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, fcb_calcMatrixPrimal_sim, hstep, 1.0_DP,&
                .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                fcb_calcMatrixDiagPrimal_sim, rcollection=rcollectionTmp)

          end select
          
        else ! callback function not present

          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
          call sys_halt()

        end if

        !-----------------------------------------------------------------------

      case (VELOCITY_ZERO)
        ! zero velocity, do nothing

        !-----------------------------------------------------------------------

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP)
        ! linear velocity

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)
          
          ! Apply standard Galerkin discretisation
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatGalConvP1d_sim, hstep, 1.0_DP,&
                .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                transp_calcMatDiagConvP1d_sim, rcollection=rcollectionTmp)
            
          case (NDIM2D)
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatGalConvP2d_sim, hstep, 1.0_DP,&
                .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                transp_calcMatDiagConvP2d_sim, rcollection=rcollectionTmp)
            
          case (NDIM3D)
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatGalConvP3d_sim, hstep, 1.0_DP,&
                .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                transp_calcMatDiagConvP3d_sim, rcollection=rcollectionTmp)
          end select
          
        case default
          
          ! Apply low-order discretisation
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwConvP1d_sim, hstep, 1.0_DP,&
                .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                transp_calcMatDiagConvP1d_sim, rcollection=rcollectionTmp)
            
          case (NDIM2D)
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwConvP2d_sim, hstep, 1.0_DP,&
                .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                transp_calcMatDiagConvP2d_sim, rcollection=rcollectionTmp)
            
          case (NDIM3D)
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwConvP3d_sim, hstep, 1.0_DP,&
                .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                transp_calcMatDiagConvP3d_sim, rcollection=rcollectionTmp)
          end select

        end select

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers` equation in space-time

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
          
        case (AFCSTAB_GALERKIN)
          
          ! Apply standard Galerkin discretisation
          call gfsc_buildJacobian(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatGalSTBurgP2d_sim, hstep, 1.0_DP,&
              .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
              transp_calcMatDiagSTBurgP2d_sim, rcollection=rcollectionTmp)

        case default
          
          ! Apply low-order discretisation
          call gfsc_buildJacobian(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwSTBurgP2d_sim, hstep, 1.0_DP,&
              .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
              transp_calcMatDiagSTBurgP2d_sim, rcollection=rcollectionTmp)

        end select

        !-----------------------------------------------------------------------

      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
          
        case (AFCSTAB_GALERKIN)
          
          ! Apply standard Galerkin discretisation
          call gfsc_buildJacobian(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatGalSTBLevP2d_sim, hstep, 1.0_DP,&
              .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
              transp_calcMatDiagSTBLevP2d_sim, rcollection=rcollectionTmp)

        case default
          
          ! Apply low-order discretisation
          call gfsc_buildJacobian(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwSTBLevP2d_sim, hstep, 1.0_DP,&
              .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
              transp_calcMatDiagSTBLevP2d_sim, rcollection=rcollectionTmp)
          
        end select
        
        !-----------------------------------------------------------------------
        
      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers` equation in 1D

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)

          ! Apply standard Galerkin discretisation
          call gfsc_buildJacobian(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatGalBurgP1d_sim, hstep, 1.0_DP,&
              .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
              transp_calcMatDiagBurgP1d_sim, rcollection=rcollectionTmp)

        case default

          ! Apply low-order discretisation
          call gfsc_buildJacobian(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwBurgP1d_sim, hstep, 1.0_DP,&
              .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
              transp_calcMatDiagBurgP1d_sim, rcollection=rcollectionTmp)

        end select
          
        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers` equation in 2D

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)

          ! Apply standard Galerkin discretisation
          call gfsc_buildJacobian(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatGalBurgP2d_sim, hstep, 1.0_DP,&
              .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
              transp_calcMatDiagBurgP2d_sim, rcollection=rcollectionTmp)

        case default

          ! Apply low-order discretisation
          call gfsc_buildJacobian(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwBurgP2d_sim, hstep, 1.0_DP,&
              .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
              transp_calcMatDiagBurgP2d_sim, rcollection=rcollectionTmp)

        end select

        !-----------------------------------------------------------------------

      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)

          ! Apply standard Galerkin discretisation
          call gfsc_buildJacobian(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatGalBLevP1d_sim, hstep, 1.0_DP,&
              .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
              transp_calcMatDiagBLevP1d_sim, rcollection=rcollectionTmp)

        case default
          
          ! Apply low-order discretisation
          call gfsc_buildJacobian(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwBLevP1d_sim, hstep, 1.0_DP,&
              .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
              transp_calcMatDiagBLevP1d_sim, rcollection=rcollectionTmp)

        end select

      end select

      !-------------------------------------------------------------------------

      ! Evaluate bilinear form for boundary integral (if any)
      call transp_calcBilfBdrCondition(rproblemLevel,&
          rsolver%rboundaryCondition, rsolution, ssectionName,&
          smode, ivelocitytype, rtimestep%dTime, 1.0_DP, .false.,&
          rproblemLevel%RmatrixScalar(transportMatrix), rcollection,&
          fcb_coeffMatBdrPrimal1d_sim, fcb_coeffMatBdrDual1d_sim,&
          fcb_coeffMatBdrPrimal2d_sim, fcb_coeffMatBdrDual2d_sim,&
          fcb_coeffMatBdrPrimal3d_sim, fcb_coeffMatBdrDual3d_sim)


    elseif (trim(smode) .eq. 'dual') then

      !-------------------------------------------------------------------------
      ! We are in dual mode which means that we have to build two
      ! bilinear forms: one for the volume integral and one for the
      ! surface integral which is used to weakly impose Dirichlet
      ! boundary conditions. Note that lumping is performed for the
      ! boundary term to prevent the generation of oscillations.
      !-------------------------------------------------------------------------

      ! @FAQ2: What type of velocity are we?
      select case (abs(ivelocitytype))
      case default
        ! The user-defined callback function for matrix coefficients
        ! is used if present; otherwise an error is thrown
        if (present(fcb_calcMatrixDual_sim)) then

          select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
            
          case (AFCSTAB_GALERKIN)
            
            ! Apply standard Galerkin discretisation
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, fcb_calcMatrixDual_sim, hstep, 1.0_DP,&
                .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                fcb_calcMatrixDiagDual_sim, rcollection=rcollectionTmp)

          case default
            
            ! Apply low-order discretisation
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, fcb_calcMatrixDual_sim, hstep, 1.0_DP,&
                .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                fcb_calcMatrixDiagDual_sim, rcollection=rcollectionTmp)
            
          end select
          
        else ! callback function not present

          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
          call sys_halt()

        end if

      !-------------------------------------------------------------------------

      case (VELOCITY_ZERO)
        ! zero velocity, do nothing

      !-------------------------------------------------------------------------

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP)
        ! linear velocity

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)
          
          ! Apply standard Galerkin discretisation
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatGalConvD1d_sim, hstep, 1.0_DP,&
                .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                transp_calcMatDiagConvD1d_sim, rcollection=rcollectionTmp)
            
          case (NDIM2D)
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatGalConvD2d_sim, hstep, 1.0_DP,&
                .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                transp_calcMatDiagConvD2d_sim, rcollection=rcollectionTmp)
            
          case (NDIM3D)
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatGalConvD3d_sim, hstep, 1.0_DP,&
                .false., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                transp_calcMatDiagConvD3d_sim, rcollection=rcollectionTmp)
          end select

        case default

          ! Apply low-order discretisation
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwConvD1d_sim, hstep, 1.0_DP,&
                .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                transp_calcMatDiagConvD1d_sim, rcollection=rcollectionTmp)
            
          case (NDIM2D)
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwConvD2d_sim, hstep, 1.0_DP,&
                .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                transp_calcMatDiagConvD2d_sim, rcollection=rcollectionTmp)
            
          case (NDIM3D)
            call gfsc_buildJacobian(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwConvD3d_sim, hstep, 1.0_DP,&
                .true., .false., rproblemLevel%RmatrixScalar(transportMatrix),&
                transp_calcMatDiagConvD3d_sim, rcollection=rcollectionTmp)
          end select

        end select

        ! @TODO: The dual mode has only been implemented for linear
        ! convection. If you need to compute the dual problem for
        ! some other velocity type, then you have to add the
        ! implementation of the dual transport operator below!

      end select

      !-------------------------------------------------------------------------

      ! Evaluate bilinear form for boundary integral (if any)
      call transp_calcBilfBdrCondition(rproblemLevel,&
          rsolver%rboundaryCondition, rsolution, ssectionName,&
          smode, ivelocitytype, rtimestep%dTime, 1.0_DP, .false.,&
          rproblemLevel%RmatrixScalar(transportMatrix), rcollection,&
          fcb_coeffMatBdrPrimal1d_sim, fcb_coeffMatBdrDual1d_sim,&
          fcb_coeffMatBdrPrimal2d_sim, fcb_coeffMatBdrDual2d_sim,&
          fcb_coeffMatBdrPrimal3d_sim, fcb_coeffMatBdrDual3d_sim)

    else
      call output_line('Invalid mode!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
      call sys_halt()
    end if


    ! Check if the Jacobian operator has extended sparsity pattern
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'ijacobianFormat', ijacobianFormat)
    if (ijacobianFormat .eq. 0) then
      bisExactStructure   = .true.
      bisExtendedSparsity = .false.
    else
      bisExactStructure   = .false.
      bisExtendedSparsity = .true.
    end if


    !---------------------------------------------------------------------------
    ! Assemble the global system operator for the high-/low-order contribution
    !---------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'imasstype', imasstype)

    select case(imasstype)
    case (MASS_LUMPED)

      !-------------------------------------------------------------------------
      ! Compute the global Jacobian for transient flow
      !
      !   $ J = M_L - scaleImplicit * dt * L $
      !-------------------------------------------------------------------------

      ! Get additional parameters from parameter list
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)

      ! Build global Jacobian matrix
      call lsyssc_MatrixLinearComb(&
          rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
          rproblemLevel%RmatrixScalar(transportMatrix),&
          1.0_DP, -rtimestep%dscaleImplicit*rtimestep%dStep,&
          .false., .false., .true., bisExactStructure,&
          rproblemLevel%RmatrixScalar(jacobianMatrix))

    case (MASS_CONSISTENT)

      !-------------------------------------------------------------------------
      ! Compute the global Jacobian for transient flow
      !
      !   $ J = M_C - scaleImplicit * dt * L $
      !-------------------------------------------------------------------------

      ! Get additional parameters from parameter list
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'consistentmassmatrix', consistentMassMatrix)

      ! Build global Jacobian matrix
      call lsyssc_MatrixLinearComb(&
          rproblemLevel%RmatrixScalar(consistentMassMatrix),&
          rproblemLevel%RmatrixScalar(transportMatrix),&
          1.0_DP, -rtimestep%dscaleImplicit*rtimestep%dStep,&
          .false., .false., .true., bisExactStructure,&
          rproblemLevel%RmatrixScalar(jacobianMatrix))

    case default

      !-------------------------------------------------------------------------
      ! Compute the global Jacobian for steady-state flow
      !
      !   $ J = -L $
      !-------------------------------------------------------------------------

      call lsyssc_MatrixLinearComb(&
          rproblemLevel%RmatrixScalar(transportMatrix),&
          rproblemLevel%RmatrixScalar(jacobianMatrix),&
          -1.0_DP, 0.0_DP, .false., .false., .true., bisExactStructure)

    end select

    !---------------------------------------------------------------------------
    ! Assemble the Jacobian matrix for the diffusion operator
    !---------------------------------------------------------------------------

    ! FAQ2: What kind of diffusion are we?
    select case(idiffusiontype)
    case (DIFFUSION_ZERO, DIFFUSION_ISOTROPIC)
      ! zero diffusion or isotropic diffusion, do nothing

    case (DIFFUSION_ANISOTROPIC)
      ! Anisotropic diffusion
      select case(rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType)
      case (AFCSTAB_SYMMETRIC)
        call afcsc_buildJacobianSymm(&
            rsolution, 1.0_DP, hstep, .false.,&
            rproblemLevel%Rafcstab(diffusionAFC),&
            rproblemLevel%RmatrixScalar(jacobianMatrix),&
            bisExtendedSparsity)
      end select

    case (DIFFUSION_VARIABLE)
      call output_line('Variable diffusion matrices are yet not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
      call sys_halt()

    case default
      call output_line('Invalid type of diffusion!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Assemble the Jacobian matrix for the convection operator
    !---------------------------------------------------------------------------

    ! Are we in primal or dual mode?
    if (trim(smode) .eq. 'primal') then
      
      !-------------------------------------------------------------------------
      ! We are in primal mode
      !-------------------------------------------------------------------------

      ! @FAQ2: What type of velocity are we?
      select case(abs(ivelocitytype))
      case default
        ! The user-defined callback function for matrix coefficients
        ! is used if present; otherwise an error is thrown
        if (present(fcb_calcMatrixPrimal_sim)) then

          select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
          case (AFCSTAB_NLINFCT_EXPLICIT,&
                AFCSTAB_NLINFCT_IMPLICIT,&
                AFCSTAB_NLINFCT_ITERATIVE)

            call parlst_getvalue_int(p_rparlist,&
                ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

            ! Should we apply consistent mass antidiffusion?
            if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
              call afcsc_buildJacobianFCT(&
                  rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                  rsolution, fcb_calcMatrixPrimal_sim,&
                  rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%RmatrixScalar(jacobianMatrix),&
                  rproblemLevel%RmatrixScalar(consistentMassMatrix))
            else
              call afcsc_buildJacobianFCT(&
                  rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                  rsolution, fcb_calcMatrixPrimal_sim,&
                  rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%RmatrixScalar(jacobianMatrix))
            end if

          case (AFCSTAB_TVD)
            call afcsc_buildJacobianTVD(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, fcb_calcMatrixPrimal_sim,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
            
          case (AFCSTAB_GP)
            
            call parlst_getvalue_int(p_rparlist,&
                ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

            ! Should we apply consistent mass antidiffusion?
            if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
              call afcsc_buildJacobianGP(&
                  rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                  rproblemLevel%RmatrixScalar(consistentMassMatrix),&
                  rsolution, rsolution0, fcb_calcMatrixPrimal_sim,&
                  rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                  rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%RmatrixScalar(jacobianMatrix),&
                  bisExtendedSparsity)
            else
              call afcsc_buildJacobianTVD(&
                  rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                  rsolution, fcb_calcMatrixPrimal_sim,&
                  rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%RmatrixScalar(jacobianMatrix),&
                  bisExtendedSparsity)
            end if

          end select

        else ! callback function not present

          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
          call sys_halt()

        end if


      case (VELOCITY_ZERO)
        ! zero velocity, do nothing


      case(VELOCITY_CONSTANT,&
           VELOCITY_TIMEDEP)
        ! linear velocity

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
        case (AFCSTAB_NLINFCT_EXPLICIT,&
              AFCSTAB_NLINFCT_IMPLICIT,&
              AFCSTAB_NLINFCT_ITERATIVE)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianFCT(rsolution,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                rproblemLevel%RmatrixScalar(consistentMassMatrix))
          else
            call afcsc_buildJacobianFCT(rsolution,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix))
          end if


        case (AFCSTAB_TVD)
          call afcsc_buildJacobianTVD(rsolution,&
              rtimestep%dStep, hstep, .false.,&
              rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%RmatrixScalar(jacobianMatrix),&
              bisExtendedSparsity)

        case (AFCSTAB_GP)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianGP(&
                rproblemLevel%RmatrixScalar(consistentMassMatrix),&
                rsolution, rsolution0,&
                rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          else
            call afcsc_buildJacobianTVD(rsolution,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          end if
        end select


      case(VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers` equation in space-time

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
        case (AFCSTAB_NLINFCT_EXPLICIT,&
              AFCSTAB_NLINFCT_IMPLICIT,&
              AFCSTAB_NLINFCT_ITERATIVE)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianFCT(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwSTBurgP2d_sim,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                rproblemLevel%RmatrixScalar(consistentMassMatrix))
          else
            call afcsc_buildJacobianFCT(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwSTBurgP2d_sim,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix))
          end if

        case (AFCSTAB_TVD)
          call afcsc_buildJacobianTVD(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwSTBurgP2d_sim,&
              rtimestep%dStep, hstep, .false.,&
              rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%RmatrixScalar(jacobianMatrix),&
              bisExtendedSparsity)

        case (AFCSTAB_GP)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianGP(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rproblemLevel%RmatrixScalar(consistentMassMatrix),&
                rsolution, rsolution0, transp_calcMatUpwSTBurgP2d_sim,&
                rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          else
            call afcsc_buildJacobianTVD(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwSTBurgP2d_sim,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          end if
        end select


      case(VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
        case (AFCSTAB_NLINFCT_EXPLICIT,&
              AFCSTAB_NLINFCT_IMPLICIT,&
              AFCSTAB_NLINFCT_ITERATIVE)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianFCT(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwSTBLevP2d_sim,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                rproblemLevel%RmatrixScalar(consistentMassMatrix))
          else
            call afcsc_buildJacobianFCT(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwSTBLevP2d_sim,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix))
          end if

        case (AFCSTAB_TVD)
          call afcsc_buildJacobianTVD(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwSTBLevP2d_sim,&
              rtimestep%dStep, hstep, .false.,&
              rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%RmatrixScalar(jacobianMatrix),&
              bisExtendedSparsity)

        case (AFCSTAB_GP)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianGP(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rproblemLevel%RmatrixScalar(consistentMassMatrix),&
                rsolution, rsolution0, transp_calcMatUpwSTBLevP2d_sim,&
                rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          else
            call afcsc_buildJacobianTVD(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwSTBLevP2d_sim,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          end if
        end select


      case(VELOCITY_BURGERS1D)
        ! nonlinear Burgers` equation in 1D

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
        case (AFCSTAB_NLINFCT_EXPLICIT,&
              AFCSTAB_NLINFCT_IMPLICIT,&
              AFCSTAB_NLINFCT_ITERATIVE)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianFCT(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwBurgP1d_sim,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                rproblemLevel%RmatrixScalar(consistentMassMatrix))
          else
            call afcsc_buildJacobianFCT(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwBurgP1d_sim,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix))
          end if

        case (AFCSTAB_TVD)
          call afcsc_buildJacobianTVD(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwBurgP1d_sim,&
              rtimestep%dStep, hstep, .false.,&
              rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%RmatrixScalar(jacobianMatrix),&
              bisExtendedSparsity)

        case (AFCSTAB_GP)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianGP(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rproblemLevel%RmatrixScalar(consistentMassMatrix),&
                rsolution, rsolution0, transp_calcMatUpwBurgP1d_sim,&
                rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          else
            call afcsc_buildJacobianTVD(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwBurgP1d_sim,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          end if
        end select


      case(VELOCITY_BURGERS2D)
        ! nonlinear Burgers` equation in 2D

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
        case (AFCSTAB_NLINFCT_EXPLICIT,&
              AFCSTAB_NLINFCT_IMPLICIT,&
              AFCSTAB_NLINFCT_ITERATIVE)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianFCT(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwBurgP2d_sim,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                rproblemLevel%RmatrixScalar(consistentMassMatrix))
          else
            call afcsc_buildJacobianFCT(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwBurgP2d_sim,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix))
          end if

        case (AFCSTAB_TVD)
          call afcsc_buildJacobianTVD(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwBurgP2d_sim,&
              rtimestep%dStep, hstep, .false.,&
              rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%RmatrixScalar(jacobianMatrix),&
              bisExtendedSparsity)

        case (AFCSTAB_GP)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianGP(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rproblemLevel%RmatrixScalar(consistentMassMatrix),&
                rsolution, rsolution0, transp_calcMatUpwBurgP2d_sim,&
                rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix), bisExtendedSparsity)
          else
            call afcsc_buildJacobianTVD(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwBurgP2d_sim,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          end if
        end select


      case(VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
        case (AFCSTAB_NLINFCT_EXPLICIT,&
              AFCSTAB_NLINFCT_IMPLICIT,&
              AFCSTAB_NLINFCT_ITERATIVE)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianFCT(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwBLevP1d_sim,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                rproblemLevel%RmatrixScalar(consistentMassMatrix))
          else
            call afcsc_buildJacobianFCT(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwBLevP1d_sim,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix))
          end if

        case (AFCSTAB_TVD)
          call afcsc_buildJacobianTVD(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwBLevP1d_sim,&
              rtimestep%dStep, hstep, .false.,&
              rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%RmatrixScalar(jacobianMatrix),&
              bisExtendedSparsity)

        case (AFCSTAB_GP)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianGP(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rproblemLevel%RmatrixScalar(consistentMassMatrix),&
                rsolution, rsolution0, transp_calcMatUpwBLevP1d_sim,&
                rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          else
            call afcsc_buildJacobianTVD(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwBLevP1d_sim,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          end if
        end select

      end select


    elseif (trim(smode) .eq. 'dual') then

      !-------------------------------------------------------------------------
      ! We are in dual mode
      !-------------------------------------------------------------------------

      ! @FAQ2: What type of velocity are we?
      select case(abs(ivelocitytype))
      case default
        ! The user-defined callback function for matrix coefficients
        ! is used if present; otherwise an error is thrown
        if (present(fcb_calcMatrixDual_sim)) then

          select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
          case (AFCSTAB_NLINFCT_EXPLICIT,&
                AFCSTAB_NLINFCT_IMPLICIT,&
                AFCSTAB_NLINFCT_ITERATIVE)

            call parlst_getvalue_int(p_rparlist,&
                ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

            ! Should we apply consistent mass antidiffusion?
            if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
              call afcsc_buildJacobianFCT(&
                  rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                  rsolution, fcb_calcMatrixDual_sim,&
                  rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%RmatrixScalar(jacobianMatrix),&
                  rproblemLevel%RmatrixScalar(consistentMassMatrix))
            else
              call afcsc_buildJacobianFCT(&
                  rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                  rsolution, fcb_calcMatrixDual_sim,&
                  rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%RmatrixScalar(jacobianMatrix))
            end if

          case (AFCSTAB_TVD)
            call afcsc_buildJacobianTVD(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, fcb_calcMatrixDual_sim,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
            
          case (AFCSTAB_GP)
            
            call parlst_getvalue_int(p_rparlist,&
                ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

            ! Should we apply consistent mass antidiffusion?
            if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
              call afcsc_buildJacobianGP(&
                  rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                  rproblemLevel%RmatrixScalar(consistentMassMatrix),&
                  rsolution, rsolution0, fcb_calcMatrixDual_sim,&
                  rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                  rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%RmatrixScalar(jacobianMatrix),&
                  bisExtendedSparsity)
            else
              call afcsc_buildJacobianTVD(&
                  rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                  rsolution, fcb_calcMatrixDual_sim,&
                  rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%RmatrixScalar(jacobianMatrix),&
                  bisExtendedSparsity)
            end if
          end select

        else ! callback function not present

          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
          call sys_halt()

        end if


      case (VELOCITY_ZERO)
        ! zero velocity, do nothing


      case(VELOCITY_CONSTANT,&
           VELOCITY_TIMEDEP)
        ! linear velocity

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
        case (AFCSTAB_NLINFCT_EXPLICIT,&
              AFCSTAB_NLINFCT_IMPLICIT,&
              AFCSTAB_NLINFCT_ITERATIVE)

          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianFCT(rsolution,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                rproblemLevel%RmatrixScalar(consistentMassMatrix))
          else
            call afcsc_buildJacobianFCT(rsolution,&
                rtimestep%p_rthetaScheme%theta, rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix))
          end if


        case (AFCSTAB_TVD)
          call afcsc_buildJacobianTVD(rsolution,&
              rtimestep%dStep, hstep, .false.,&
              rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%RmatrixScalar(jacobianMatrix),&
              bisExtendedSparsity)

        case (AFCSTAB_GP)
          call parlst_getvalue_int(p_rparlist,&
              ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call afcsc_buildJacobianGP(&
                rproblemLevel%RmatrixScalar(consistentMassMatrix),&
                rsolution, rsolution0,&
                rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          else
            call afcsc_buildJacobianTVD(rsolution,&
                rtimestep%dStep, hstep, .false.,&
                rproblemLevel%Rafcstab(convectionAFC),&
                rproblemLevel%RmatrixScalar(jacobianMatrix),&
                bisExtendedSparsity)
          end if
        end select

      end select

    else
      call output_line('Invalid mode!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
      call sys_halt()
    end if


    ! Impose boundary conditions in strong sence (if any)
    if (rsolver%rboundaryCondition%bStrongBdrCond) then
      call bdrf_filterMatrix(rsolver%rboundaryCondition,&
          rproblemLevel%RmatrixScalar(jacobianMatrix), 1.0_DP)
    end if

    ! Ok, we updated the Jacobian matrix successfully. Now we still have to
    ! link it to the solver hierarchy. This is done recursively.

    ! What type of flow are we?
    select case(imasstype)
    case (MASS_LUMPED,&
          MASS_CONSISTENT)
      call flagship_updateSolverMatrix(rproblemLevel, rsolver,&
          jacobianMatrix, SYSTEM_INTERLEAVEFORMAT,&
          UPDMAT_JAC_TRANSIENT, rproblemLevel%ilev, rproblemLevel%ilev)

    case default
      call flagship_updateSolverMatrix(rproblemLevel, rsolver,&
          jacobianMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_JAC_STEADY,&
          rproblemLevel%ilev, rproblemLevel%ilev)
    end select

    ! Finally, we have to update the content of the solver hierarchy
    call solver_updateContent(rsolver)

    ! Stop time measurement for matrix evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine transp_calcJacobian

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcRhs(rproblemLevel, rtimestep,&
      rsolver, rsolution, rrhs, ssectionName, rcollection, rsource,&
      fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
      fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
      fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim)

!<description>
    ! This subroutine computes the constant right-hand side
    !
    !  $$ rhs = [M + scaleExplicit * dt * T^n]u^n + s^n + b.c.`s$$
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

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCoeffVecBdr.inc'
    optional :: fcb_coeffVecBdrPrimal1d_sim
    optional :: fcb_coeffVecBdrPrimal2d_sim
    optional :: fcb_coeffVecBdrPrimal3d_sim
    optional :: fcb_coeffVecBdrDual1d_sim
    optional :: fcb_coeffVecBdrDual2d_sim
    optional :: fcb_coeffVecBdrDual3d_sim
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! right-hand side vector
    type(t_vectorBlock), intent(inout) :: rrhs

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rpredictor
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    character(LEN=SYS_STRLEN) :: smode
    real(DP) :: dscale
    integer :: consistentMassMatrix, lumpedMassMatrix
    integer :: transportMatrix, massMatrix
    integer :: imasstype, imassantidiffusiontype, ivelocitytype
    integer :: convectionAFC, diffusionAFC


    !---------------------------------------------------------------------------
    ! We calculate the constant right-hand side vector based on the
    ! transport operator. If boundary conditions are imposed in weak
    ! sense, then additional surface integrals are applied to the
    ! right-hand side of the primal problem. In the dual problem,
    ! weakly imposed boundary conditions are built into the target
    ! functional which is supplied via the vector 'rb'.
    ! ---------------------------------------------------------------------------

    ! Start time measurement for rhs evaluation
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
        ssectionName, 'transportmatrix', transportMatrix)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'imasstype', imasstype)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_string(p_rparlist,&
        ssectionName, 'mode', smode)

    ! Do we have some kind of mass matrix?
    select case(imasstype)
    case (MASS_LUMPED, MASS_CONSISTENT)

      !-------------------------------------------------------------------------
      ! Compute the constant right-hand side
      !
      !   $$ rhs = M*u^n + scaleExplicit * dt * T(u^n) * u^n + b.c.`s $$
      !-------------------------------------------------------------------------

      ! Do we have an explicit part?
      if (rtimestep%dscaleExplicit .ne. 0.0_DP) then

        ! Compute scaling parameter
        dscale = rtimestep%dscaleExplicit * rtimestep%dStep

        ! Build transport term $scaleExplicit * dt * T(u^n) * u^n$, where
        ! $T(u^n)$ denotes the discrete transport operator
        call lsyssc_matVec(&
            rproblemLevel%RmatrixScalar(transportMatrix),&
            rsolution%rvectorBlock(1),&
            rrhs%RvectorBlock(1), dscale, 0.0_DP)
        
        ! Build transient term $M_L*u^n$ or $M_C*u^n$
        massMatrix = merge(lumpedMassMatrix, consistentMassMatrix,&
                           imasstype .eq. MASS_LUMPED)

        call lsyssc_matVec(rproblemLevel%RmatrixScalar(massMatrix),&
            rsolution%RvectorBlock(1), rrhs%RvectorBlock(1), 1.0_DP, 1.0_DP)

        ! Evaluate linear form for the boundary integral (if any)
        call transp_calcLinfBdrCondition(rproblemLevel,&
            rsolver%rboundaryCondition, rsolution, ssectionName, smode,&
            ivelocitytype, rtimestep%dTime-rtimestep%dStep,&
            dscale, .false.,rrhs, rcollection,&
            fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
            fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
            fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim)

        ! Build the geometric source term (if any)
        call transp_calcGeometricSourceterm(p_rparlist, ssectionName,&
            rproblemLevel, rsolution, dscale, .false., rrhs, rcollection)
        
        !-----------------------------------------------------------------------
        ! Apply algebraic flux correction for the convection term (if required).
        !
        !   $$ rhs = rhs + scaleExplicit * dt * fconv(u^n)
        !
        ! Here, only the explicit part is built into the constant
        ! right-hand side vector for those limiting schemes which allow
        ! for a strict separation of explicit and implicit contributions.
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'convectionAFC', convectionAFC, 0)
        
        if (convectionAFC .gt. 0) then
          
          ! What type of stabilisation are we?
          select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
            
          case (AFCSTAB_TVD)
            call afcsc_buildVectorTVD(&
                rproblemLevel%Rafcstab(convectionAFC),&
                rsolution, dscale, .false., AFCSTAB_TVDALGO_STANDARD, rrhs)

            !-------------------------------------------------------------------

          case (AFCSTAB_NLINLPT_UPWINDBIASED)

            ! Check if slope-based local-extremum diminishing bounds
            ! at edges are available and compute them otherwise
            if (iand(rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec,&
                AFCSTAB_HAS_EDGEBOUNDS) .eq. 0)&
                call transp_calcBoundsAtEdgesLPT(p_rparlist, ssectionName,&
                rproblemLevel, rproblemLevel%Rafcstab(convectionAFC))

            ! Build the raw antidiffusive fluxes
            call afcsc_buildFluxLPT(rproblemLevel%Rafcstab(convectionAFC),&
                rsolution, 1.0_DP, .true., AFCSTAB_LPTFLUX+AFCSTAB_LPTFLUX_BOUNDS)

            ! Apply explicit part of the FEM-LPT correction to right-hand side
            call afcsc_buildVectorLPT(&
                rproblemLevel%Rafcstab(convectionAFC),&
                rsolution, dscale, .false.,&
                AFCSTAB_LPTALGO_STANDARD, rrhs)

            !-------------------------------------------------------------------
            ! Some limiters require preparation tasks to be performed
            !-------------------------------------------------------------------
          case (AFCSTAB_NLINFCT_EXPLICIT,&
                AFCSTAB_NLINFCT_ITERATIVE,&
                AFCSTAB_NLINFCT_IMPLICIT)

            ! Compute the low-order predictor based on the right-hand side
            ! and assemble the explicit part of the raw-antidiffusive fluxes

            ! Set pointer to predictor
            p_rpredictor => rproblemLevel%Rafcstab(convectionAFC)%p_rvectorPredictor

            ! Compute $\tilde u = (M_L)^{-1}*b^n$
            call lsysbl_invertedDiagMatVec(&
                rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
                rrhs, 1.0_DP, p_rpredictor)

            ! Set specifier
            rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec =&
                ior(rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec,&
                    AFCSTAB_HAS_PREDICTOR)

            ! Should we apply consistent mass antidiffusion?
            call parlst_getvalue_int(p_rparlist,&
                ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)
            
            if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
              ! Assemble explicit part of the raw-antidiffusive fluxes
              ! with the contribution of the consistent mass matrix
              call afcsc_buildFluxFCT(&
                  rproblemLevel%Rafcstab(convectionAFC), rsolution,&
                  rtimestep%dStep, rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                  1.0_DP, .true., .true., AFCSTAB_FCTFLUX_EXPLICIT,&
                  rmatrix=rproblemLevel%RmatrixScalar(consistentMassMatrix),&
                  rxPredictor=p_rpredictor)
            else
              ! Assemble explicit part of the raw-antidiffusive fluxes
              ! without the contribution of the consistent mass matrix
              call afcsc_buildFluxFCT(&
                  rproblemLevel%Rafcstab(convectionAFC), rsolution,&
                  1.0_DP, rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                  1.0_DP, .true., .true., AFCSTAB_FCTFLUX_EXPLICIT,&
                  rxPredictor=p_rpredictor)
            end if
            
          end select

        end if
        
        !-----------------------------------------------------------------------
        ! Apply algebraic flux correction for the diffusion term (if required).
        !
        !   $$ rhs = rhs + scaleExplicit * dt * fdiff(u^n)
        !
        ! Here, only the explicit part is built into the constant
        ! right-hand side vector for those limiting schemes which allow
        ! for a strict separation of explicit and implicit contributions.
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'diffusionAFC', diffusionAFC, 0)
        
        if (diffusionAFC .gt. 0) then
          
          ! What kind of stabilisation should be applied?
          select case(rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType)
            
          case (AFCSTAB_SYMMETRIC)
            call afcsc_buildVectorSymm(&
                rproblemLevel%Rafcstab(diffusionAFC), rsolution, dscale, rrhs)

            !-------------------------------------------------------------------
            
          case (AFCSTAB_NLINLPT_SYMMETRIC)

            ! Check if slope-based local-extremum diminishing bounds
            ! at edges are available and compute them otherwise
            if (iand(rproblemLevel%Rafcstab(diffusionAFC)%istabilisationSpec,&
                AFCSTAB_HAS_EDGEBOUNDS) .eq. 0)&
                call transp_calcBoundsAtEdgesLPT(p_rparlist, ssectionName,&
                rproblemLevel, rproblemLevel%Rafcstab(diffusionAFC))
            
            ! Build the raw antidiffusive fluxes
            call afcsc_buildFluxLPT(rproblemLevel%Rafcstab(diffusionAFC),&
                rsolution, 1.0_DP, .true., AFCSTAB_LPTFLUX+AFCSTAB_LPTFLUX_BOUNDS)

            ! Apply explicit part of the FEM-LPT correction to right-hand side
            call afcsc_buildVectorLPT(&
                rproblemLevel%Rafcstab(diffusionAFC),&
                rsolution, dscale, .false.,&
                AFCSTAB_LPTALGO_STANDARD, rrhs)

          end select
          
        end if
        
      else ! diffusionAFC .le. 0

        ! Build transient term $M_L*u^n$ or $M_C*u^n$
        massMatrix = merge(lumpedMassMatrix, consistentMassMatrix,&
                           imasstype .eq. MASS_LUMPED)

        call lsyssc_matVec(rproblemLevel%RmatrixScalar(massMatrix),&
            rsolution%RvectorBlock(1), rrhs%RvectorBlock(1), 1.0_DP, 0.0_DP)

        !-----------------------------------------------------------------------
        ! Some limiters require preparation tasks to be performed even
        ! though the explicit part of raw-antidiffusive fluxes vanishes
        !-----------------------------------------------------------------------

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'convectionAFC', convectionAFC, 0)
        
        if (convectionAFC .gt. 0) then
          
          ! What type of stabilisation are we?
          select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
            
          case (AFCSTAB_NLINFCT_EXPLICIT,&
                AFCSTAB_NLINFCT_ITERATIVE,&
                AFCSTAB_NLINFCT_IMPLICIT)

            ! Compute the low-order predictor based on the right-hand side
            ! and assemble the explicit part of the raw-antidiffusive fluxes
            
            ! Set pointer to predictor
            p_rpredictor => rproblemLevel%Rafcstab(convectionAFC)%p_rvectorPredictor
            
            ! Compute $\tilde u = (M_L)^{-1}*b^n = u^n$
            call lsysbl_copyVector(rsolution, p_rpredictor)
            
            ! Set specifier
            rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec =&
                ior(rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec,&
                    AFCSTAB_HAS_PREDICTOR)

            ! Should we apply consistent mass antidiffusion?
            call parlst_getvalue_int(p_rparlist,&
                ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

            if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
              ! Assemble explicit part of the raw-antidiffusive fluxes
              ! with the contribution of the consistent mass matrix
              call afcsc_buildFluxFCT(&
                  rproblemLevel%Rafcstab(convectionAFC), rsolution,&
                  rtimestep%dStep, rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                  1.0_DP, .true., .true., AFCSTAB_FCTFLUX_EXPLICIT,&
                  rmatrix=rproblemLevel%RmatrixScalar(consistentMassMatrix),&
                  rxPredictor=p_rpredictor)
            else
              ! Assemble explicit part of the raw-antidiffusive fluxes
              ! without the contribution of the consistent mass matrix
              call afcsc_buildFluxFCT(&
                  rproblemLevel%Rafcstab(convectionAFC), rsolution,&
                  1.0_DP, rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
                  1.0_DP, .true., .true., AFCSTAB_FCTFLUX_EXPLICIT)
            end if

          end select
        end if

      end if ! diffusionAFC


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

  end subroutine transp_calcRhs

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcResidual(rproblemLevel,&
      rtimestep, rsolver, rsolution, rsolution0, rrhs, rres,&
      ite, ssectionName, rcollection, rsource,&
      fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
      fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
      fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim)

!<description>
    ! This subroutine computes the nonlinear residual vector
    !
    !   $$ res^{(m)} = rhs - [M - scaleImplicit * dt T^{(m)}] * u^{(m)} - s^{(m)} + b.c.`s $$
    !
    ! whereby the (scaled) source term is optional. The constant
    ! right-hand side
    !
    !  $$ rhs = [M + scaleExplicit * dt * T^n] * u^n + s^n + b.c.`s $$
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

    ! iteration number
    integer, intent(in) :: ite

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: source vector
    type(t_vectorBlock), intent(in), optional :: rsource

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCoeffVecBdr.inc'
    optional :: fcb_coeffVecBdrPrimal1d_sim
    optional :: fcb_coeffVecBdrPrimal2d_sim
    optional :: fcb_coeffVecBdrPrimal3d_sim
    optional :: fcb_coeffVecBdrDual1d_sim
    optional :: fcb_coeffVecBdrDual2d_sim
    optional :: fcb_coeffVecBdrDual3d_sim
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
    type(t_vectorBlock), pointer :: p_rpredictor
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    character(LEN=SYS_STRLEN) :: smode
    real(DP) :: dscale
    integer(I32) :: ioperationSpec
    integer :: transportMatrix, massMatrix
    integer :: consistentMassMatrix, lumpedMassMatrix
    integer :: massAFC, convectionAFC, diffusionAFC
    integer :: imasstype, imassantidiffusiontype, ivelocitytype


    !---------------------------------------------------------------------------
    ! We calculate the nonlinear residual vector based on the
    ! transport operator. If boundary conditions are imposed in weak
    ! sense, then they give rise to additional surface integrals in
    ! the bilinear form, i.e., the transport operator so that no
    ! additional linear forms for boundary conditions are required
    !---------------------------------------------------------------------------

    ! Start time measurement for residual evaluation
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
        ssectionName, 'transportmatrix', transportMatrix)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'imasstype', imasstype)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_string(p_rparlist,&
        ssectionName, 'mode', smode)

    !---------------------------------------------------------------------------
    ! Initialise the residual by the constant right-hand side
    !
    !   $$ res := rhs $$
    !---------------------------------------------------------------------------
    call lsysbl_copyVector(rrhs, rres)

    ! Do we have some kind of mass matrix?
    select case(imasstype)
    case (MASS_LUMPED, MASS_CONSISTENT)

      !-------------------------------------------------------------------------
      ! Compute the residual for transient flows
      !
      !   $$ res := res - [M - scaleImplicit * dt * K(u^{(m)})] * u^{(m)} + b.c.`s $$
      !-------------------------------------------------------------------------

      ! Do we have an implicit part?
      if (rtimestep%dscaleImplicit .ne. 0.0_DP) then
        
        ! Compute scaling parameter
        dscale = rtimestep%dscaleImplicit*rtimestep%dStep
        
        ! Apply transport operator to the solution vector
        call lsyssc_matVec(rproblemLevel%RmatrixScalar(transportMatrix),&
            rsolution%rvectorBlock(1), rres%RvectorBlock(1), dscale, 1.0_DP)
        
        ! Evaluate linear form for the boundary integral (if any)
        call transp_calcLinfBdrCondition(rproblemLevel,&
            rsolver%rboundaryCondition, rsolution, ssectionName, smode,&
            ivelocitytype, rtimestep%dTime, dscale, .false., rres, rcollection,&
            fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
            fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
            fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim)

        ! Build the geometric source term (if any)
        call transp_calcGeometricSourceterm(p_rparlist, ssectionName,&
            rproblemLevel, rsolution, dscale, .false., rres, rcollection)

      else

        ! Set scaling parameter to zero
        dscale = 0.0_DP
        
      end if

      ! What type of mass matrix should be used?
      massMatrix = merge(lumpedMassMatrix, consistentMassMatrix,&
                         imasstype .eq. MASS_LUMPED)

      ! Apply mass matrix to the solution vector
      call lsyssc_matVec(rproblemLevel%RmatrixScalar(massMatrix),&
          rsolution%RvectorBlock(1), rres%RvectorBlock(1), -1.0_DP, 1.0_DP)
      

    case default

      !-------------------------------------------------------------------------
      ! Compute the residual for stationary flows
      !
      !   $$ res := res + K(u^{(m)})*u^{(m)} + b.c.`s $$
      !-------------------------------------------------------------------------

      ! Set scaling parameter
      dscale = 1.0_DP

      ! Apply transport operator to the solution vector
      call lsyssc_matVec(rproblemLevel%RmatrixScalar(transportMatrix),&
          rsolution%rvectorBlock(1), rres%RvectorBlock(1), dscale, 1.0_DP)

      ! Evaluate linear form for boundary integral (if any)
      call transp_calcLinfBdrCondition(rproblemLevel,&
          rsolver%rboundaryCondition, rsolution, ssectionName, smode,&
          ivelocitytype, rtimestep%dTime, dscale, .false., rres, rcollection,&
          fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
          fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
          fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim)

      ! Build the geometric source term (if any)
      call transp_calcGeometricSourceterm(p_rparlist, ssectionName,&
          rproblemLevel, rsolution, dscale, .false., rres, rcollection)
      
    end select

    !-------------------------------------------------------------------------
    ! Perform algebraic flux correction for the mass term (if required)
    !
    !   $$ res := res + dstep*fmass(u^(m),u^n) $$
    !-------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'massAFC', massAFC, 0)
    
    if (massAFC .gt. 0) then
      
      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(massAFC)%cafcstabType)
        
      case (AFCSTAB_NLINLPT_MASS)

        ! Check if slope-based local-extremum diminishing bounds
        ! at edges are available and compute them otherwise
        if (iand(rproblemLevel%Rafcstab(massAFC)%istabilisationSpec,&
            AFCSTAB_HAS_EDGEBOUNDS) .eq. 0)&
            call transp_calcBoundsAtEdgesLPT(p_rparlist, ssectionName,&
            rproblemLevel, rproblemLevel%Rafcstab(massAFC))

        ! Set pointer to predictor
        p_rpredictor => rproblemLevel%Rafcstab(massAFC)%p_rvectorPredictor

        ! Compute approximate time derivative
        call lsysbl_copyVector(rsolution0, p_rpredictor)
        call lsysbl_vectorLinearComb(rsolution, p_rpredictor,&
            1.0_DP/rtimestep%dStep, -1.0_DP/rtimeStep%dStep)

        ! Build the raw antidiffusive fluxes
        call afcsc_buildFluxLPT(rproblemLevel%Rafcstab(massAFC),&
            p_rpredictor, 1.0_DP, .true., AFCSTAB_LPTFLUX+AFCSTAB_LPTFLUX_BOUNDS,&
            rproblemLevel%RmatrixScalar(consistentMassMatrix))

        ! Apply FEM-LPT correction to the residual vector
        call afcsc_buildVectorLPT(&
            rproblemLevel%Rafcstab(massAFC),&
            p_rpredictor, rtimeStep%dStep, .false.,&
            AFCSTAB_LPTALGO_STANDARD, rres)

      end select
      
    end if

    !-------------------------------------------------------------------------
    ! Perform algebraic flux correction for the convective term (if required)
    !
    !   $$ res = res + dscale*fconv(u^(m),u^n) $$
    !
    ! or
    !
    !   $$ res = res + dscale*fconv_impl(u^(m),u^n) $$
    !   $$ rhs = rhs + dscale*fconv_expl(u^(m),u^n) $$
    !
    ! whereby the explicit part is built into the right-hand side.
    !-------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'convectionAFC', convectionAFC, 0)
    
    if (convectionAFC .gt. 0) then
      
      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
        
      case (AFCSTAB_NLINFCT_EXPLICIT,&
            AFCSTAB_NLINFCT_ITERATIVE,&
            AFCSTAB_NLINFCT_IMPLICIT)

        ! Set pointer to predictor
        p_rpredictor => rproblemLevel%Rafcstab(convectionAFC)%p_rvectorPredictor

        ! Should we apply consistent mass antidiffusion?
        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)
        
        ! Set operation specifier
        ioperationSpec = AFCSTAB_FCTFLUX_IMPLICIT
        if (ite .gt. 0)&
            ! This has only influence on iterative FCT algorithm
            ioperationSpec = ioperationSpec + AFCSTAB_FCTFLUX_REJECTED

        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          ! Assemble implicit part of the raw-antidiffusive fluxes
          ! with the contribution of the consistent mass matrix
          call afcsc_buildFluxFCT(&
              rproblemLevel%Rafcstab(convectionAFC), rsolution,&
              rtimestep%dStep, rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
              1.0_DP, .true., .true., ioperationSpec,&
              rmatrix=rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rxPredictor=p_rpredictor)
        else
          ! Assemble implicit part of the raw-antidiffusive fluxes
          ! without the contribution of the consistent mass matrix
          call afcsc_buildFluxFCT(&
              rproblemLevel%Rafcstab(convectionAFC), rsolution,&
              rtimestep%dStep, rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
              1.0_DP, .true., .true., ioperationSpec,&
              rxPredictor=p_rpredictor)
        end if
        
        ! Set operation specifier
        if (ite .eq. 0) then
          ! Perform standard flux correction in zeroth iteration
          ioperationSpec = AFCSTAB_FCTALGO_STANDARD
        else
          select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
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
        call afcsc_buildVectorFCT(&
            rproblemLevel%Rafcstab(convectionAFC),&
            rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
            p_rpredictor, rtimestep%dStep, .false.,&
            ioperationSpec, rres)

        ! Special treatment for iterative FCT-algorithm
        if (rproblemLevel%Rafcstab(convectionAFC)%cafcstabType .eq.&
            AFCSTAB_NLINFCT_ITERATIVE) then
          ! Subtract corrected antidiffusion from right-hand side
          call afcsc_buildVectorFCT(&
              rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              p_rpredictor, rtimestep%dStep, .false.,&
              AFCSTAB_FCTALGO_CORRECT, rrhs)

          ! Recompute the low-order predictor for the next limiting step
          call lsysbl_invertedDiagMatVec(&
              rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
              rrhs, 1.0_DP, p_rpredictor)
        end if

        !-----------------------------------------------------------------------

      case (AFCSTAB_TVD)
        call afcsc_buildVectorTVD(&
            rproblemLevel%Rafcstab(convectionAFC),&
            rsolution, dscale, .false., AFCSTAB_TVDALGO_STANDARD, rres)

        !-----------------------------------------------------------------------

      case (AFCSTAB_GP)

        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call afcsc_buildVectorGP(&
              rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rsolution, rsolution0,&
              rtimestep%dscaleExplicit, rtimestep%dscaleImplicit,&
              rtimestep%dStep, .false., AFCSTAB_TVDALGO_STANDARD, rres)
        else
          call afcsc_buildVectorTVD(&
              rproblemLevel%Rafcstab(convectionAFC),&
              rsolution, dscale, .false., AFCSTAB_TVDALGO_STANDARD, rres)
        end if

        !-----------------------------------------------------------------------

      case (AFCSTAB_NLINLPT_UPWINDBIASED)
        
        ! Check if slope-based local-extremum diminishing bounds
        ! at edges are available and compute them otherwise
        if (iand(rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec,&
            AFCSTAB_HAS_EDGEBOUNDS) .eq. 0)&
            call transp_calcBoundsAtEdgesLPT(p_rparlist, ssectionName,&
            rproblemLevel, rproblemLevel%Rafcstab(convectionAFC))
        
        ! Build the raw antidiffusive fluxes
        call afcsc_buildFluxLPT(rproblemLevel%Rafcstab(convectionAFC),&
            rsolution, 1.0_DP, .true., AFCSTAB_LPTFLUX+AFCSTAB_LPTFLUX_BOUNDS)
        
        ! Apply FEM-LPT correction to the residual vector
        call afcsc_buildVectorLPT(&
            rproblemLevel%Rafcstab(convectionAFC),&
            rsolution, dscale, .false.,&
            AFCSTAB_LPTALGO_STANDARD, rres)
                
      end select

    end if

    !-------------------------------------------------------------------------
    ! Perform algebraic flux correction for the diffusive term (if required)
    !
    !   $$ res = res + dscale*fdiff(u^{(m)},u^n) $$
    !-------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'diffusionAFC', diffusionAFC, 0)

    if (diffusionAFC .gt. 0) then

      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType)
        
      case (AFCSTAB_SYMMETRIC)
        call afcsc_buildVectorSymm(&
            rproblemLevel%Rafcstab(diffusionAFC), rsolution, dscale, rres)
        
        !-----------------------------------------------------------------------
        
      case (AFCSTAB_NLINLPT_SYMMETRIC)
        
        ! Check if slope-based local-extremum diminishing bounds
        ! at edges are available and compute them otherwise
        if (iand(rproblemLevel%Rafcstab(diffusionAFC)%istabilisationSpec,&
            AFCSTAB_HAS_EDGEBOUNDS) .eq. 0)&
            call transp_calcBoundsAtEdgesLPT(p_rparlist, ssectionName,&
            rproblemLevel, rproblemLevel%Rafcstab(diffusionAFC))
        
        ! Build the raw antidiffusive fluxes
        call afcsc_buildFluxLPT(rproblemLevel%Rafcstab(diffusionAFC),&
            rsolution, 1.0_DP, .true., AFCSTAB_LPTFLUX+AFCSTAB_LPTFLUX_BOUNDS)
        
        ! Apply explicit part of the FEM-LPT correction to right-hand side
        call afcsc_buildVectorLPT(&
            rproblemLevel%Rafcstab(diffusionAFC),&
            rsolution, dscale, .false.,&
            AFCSTAB_LPTALGO_STANDARD, rres)
        
      end select
      
    end if

    ! Apply the source vector to the residual (if any)
    if (present(rsource)) then
      if (rsource%NEQ .gt. 0)&
          call lsysbl_vectorLinearComb(rsource, rres, -1.0_DP, 1.0_DP)
    end if

    ! Stop time measurement for residual evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine transp_calcResidual

  !*****************************************************************************

!<subroutine>

  subroutine transp_setBoundaryCondition(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rres, ssectionName, rcollection,&
      rmatrix, rboundaryCondition)

!<description>
    ! This subroutine imposes the Dirichlet boundary conditions in
    ! strong sense by filtering the system matrix, the solution
    ! vector and/or the residual vector explicitly.
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! nonlinear solver structure
    type(t_solver), intent(in) :: rsolver

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolution0

    ! OPTIONAL: boundary condition; if not present, the boundary
    ! condition associated with the solver structure is used
    type(t_boundaryCondition), intent(in), target, optional :: rboundaryCondition
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! residual vector
    type(t_vectorBlock), intent(inout) :: rres

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! collection structure to provide additional
    ! information to the boundary setting routine
    type(t_collection), intent(InOUT) :: rcollection

    ! OPTIONAL: system matrix; if not present, the system matrix
    ! associated with the solver structure is used
    type(t_matrixScalar), intent(inout), target, optional :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_matrixScalar), pointer :: p_rmatrix
    type(t_boundaryCondition), pointer :: p_rboundaryCondition
    integer :: imatrix


    ! What system matrix should be used?
    if (present(rmatrix)) then
      p_rmatrix => rmatrix
    else
      
      ! Get parameter list
      p_rparlist => collct_getvalue_parlst(rcollection,&
          'rparlist', ssectionName=ssectionName)
      
      ! What type of preconditioner are we?
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
            OU_CLASS_ERROR, OU_MODE_STD,'transp_setBoundaryCondition')
        call sys_halt()
      end select
      
      p_rmatrix => rproblemLevel%RmatrixScalar(imatrix)
      
    end if

    ! What boundary conditions should be used?
    if (present(rboundaryCondition)) then
      p_rboundaryCondition => rboundaryCondition
    else
      p_rboundaryCondition => rsolver%rboundaryCondition
    end if


    ! Impose boundary conditions for the solution vector and impose
    ! zeros in the residual vector and the off-diagonal positions of
    ! the system matrix which depends on the nonlinear solver
    call bdrf_filterSolution(p_rboundaryCondition,&
        p_rmatrix , rsolution, rres, rsolution0, rtimestep%dTime)

  end subroutine transp_setBoundaryCondition
 
  !*****************************************************************************

!<subroutine>

  subroutine transp_calcBilfBdrCondition(rproblemLevel,&
      rboundaryCondition, rsolution, ssectionName, smode,&
      ivelocitytype, dtime, dscale, bclear, rmatrix, rcollection,&
      fcb_coeffMatBdrPrimal1d_sim, fcb_coeffMatBdrDual1d_sim,&
      fcb_coeffMatBdrPrimal2d_sim, fcb_coeffMatBdrDual2d_sim,&
      fcb_coeffMatBdrPrimal3d_sim, fcb_coeffMatBdrDual3d_sim)

!<description>
    ! This subroutine is a shortcut for building the bilinear form
    ! arising from the weak imposition of boundary conditions. It
    ! calls the more general routine transp_calcBilfBdrCondXd
    ! using the corresponding callback routines.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! boundary condition
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! problem mode (primal, dual)
    character(LEN=*), intent(in) :: smode

    ! type of velocity
    integer, intent(in) :: ivelocitytype

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! Whether to clear the matrix before calculating the entries.
    ! If .FALSE., the new matrix entries are added to the existing entries.
    logical, intent(in) :: bclear

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCoeffMatBdr.inc'
    optional :: fcb_coeffMatBdrPrimal1d_sim
    optional :: fcb_coeffMatBdrPrimal2d_sim
    optional :: fcb_coeffMatBdrPrimal3d_sim
    optional :: fcb_coeffMatBdrDual1d_sim
    optional :: fcb_coeffMatBdrDual2d_sim
    optional :: fcb_coeffMatBdrDual3d_sim
!</intput>

!<inputoutput>
    ! scalar matrix where to store the bilinear form
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    integer :: primalBdrGFEM,dualBdrGFEM
    
    
    ! Set pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    
    
    !###########################################################################
    ! REMARK: If you want to add a new type of velocity/diffusion,
    ! then search for the tag @FAQ2: in this subroutine and create a
    ! new CASE which performs the corresponding task for the new type
    ! of velocity/ diffusion.
    !###########################################################################

    ! Are we in primal or dual mode?
    if (trim(smode) .eq. 'primal') then
    
      ! Check if group finite element formulation is applicable
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'primalbdrGFEM', primalBdrGFEM, 0)
  
      ! @FAQ2: What type of velocity are we?
      select case(abs(ivelocitytype))
        
      case default
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          ! The user-defined callback function is used if present;
          ! otherwise an error is throws
          if (present(fcb_coeffMatBdrPrimal1d_sim)) then
            call transp_calcBilfBdrCond1D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, fcb_coeffMatBdrPrimal1d_sim,&
                bclear, rmatrix, rcollection)
          else ! callback function not present
            call output_line('Missing user-defined callback function!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcBilfBdrCondition')
            call sys_halt()
          end if
        
        case (NDIM2D)
          ! The user-defined callback function is used if present;
          ! otherwise an error is throws
          if (present(fcb_coeffMatBdrPrimal2d_sim)) then
            call transp_calcBilfBdrCond2D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, fcb_coeffMatBdrPrimal2d_sim,&
                bclear, rmatrix, rcollection)
          else ! callback function not present
            call output_line('Missing user-defined callback function!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcBilfBdrCondition')
            call sys_halt()
          end if

        case (NDIM3D)
          ! The user-defined callback function is used if present;
          ! otherwise an error is throws
          if (present(fcb_coeffMatBdrPrimal3d_sim)) then
!!$            call transp_calcBilfBdrCond3D(rproblemLevel,&
!!$                rboundaryContidion, rsolution, ssectionName,&
!!$                dtime, dscale, fcb_coeffMatBdrPrimal3d_sim,&
!!$                bclear, rmatrix, rcollection)
          else ! callback function not present
            call output_line('Missing user-defined callback function!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcBilfBdrCondition')
            call sys_halt()
          end if
        end select

        !-----------------------------------------------------------------------

      case (VELOCITY_ZERO, VELOCITY_CONSTANT, VELOCITY_TIMEDEP)
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          ! linear velocity in 1D
          if (primalBdrGFEM .gt. 0) then
            call transp_calcOperatorBdrCondition(rproblemLevel,&
                rproblemLevel%RgroupFEMBlock(primalBdrGFEM),&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_calcMatBdrConvP1d_sim,&
                bclear, rmatrix, rcollection)
          else
            call transp_calcBilfBdrCond1D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_coeffMatBdrConvP1d_sim,&
                bclear, rmatrix, rcollection)
          end if

        case (NDIM2D)
          ! linear velocity in 2D
          if (primalBdrGFEM .gt. 0) then
            call transp_calcOperatorBdrCondition(rproblemLevel,&
                rproblemLevel%RgroupFEMBlock(primalBdrGFEM),&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_calcMatBdrConvP2d_sim,&
                bclear, rmatrix, rcollection)
          else
            call transp_calcBilfBdrCond2D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_coeffMatBdrConvP2d_sim,&
                bclear, rmatrix, rcollection)
          end if
          
        case (NDIM3D)
          ! linear velocity in 3D
          if (primalBdrGFEM .gt. 0) then
            call transp_calcOperatorBdrCondition(rproblemLevel,&
                rproblemLevel%RgroupFEMBlock(primalBdrGFEM),&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_calcMatBdrConvP3d_sim,&
                bclear, rmatrix, rcollection)
          else
!!$          call transp_calcBilfBdrCond3D(rproblemLevel,&
!!$              rboundaryCondition, rsolution, ssectionName,&
!!$              dtime, dscale, transp_coeffMatBdrConvP3d_sim,&
!!$              bclear, rmatrix, , rcollection)
          end if
        end select

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers` equation in space-time
        call transp_calcBilfBdrCond2D(rproblemLevel,&
            rboundaryCondition, rsolution, ssectionName,&
            dtime, dscale, transp_coeffMatBdrSTBurgP2d_sim,&
            bclear, rmatrix, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time
        call transp_calcBilfBdrCond2D(rproblemLevel,&
            rboundaryCondition, rsolution, ssectionName,&
            dtime, dscale, transp_coeffMatBdrSTBLevP2d_sim,&
            bclear, rmatrix, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers` equation in 1D
        call transp_calcBilfBdrCond1D(rproblemLevel,&
            rboundaryCondition, rsolution, ssectionName,&
            dtime, dscale, transp_coeffMatBdrBurgP1d_sim,&
            bclear, rmatrix, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers` equation in 2D
        call transp_calcBilfBdrCond2D(rproblemLevel,&
            rboundaryCondition, rsolution, ssectionName,&
            dtime, dscale, transp_coeffMatBdrBurgP2d_sim,&
            bclear, rmatrix, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D
        call transp_calcBilfBdrCond1D(rproblemLevel,&
            rboundaryCondition, rsolution, ssectionName,&
            dtime, dscale, transp_coeffMatBdrBLevP1d_sim,&
            bclear, rmatrix, rcollection)

      end select

      !-----------------------------------------------------------------------

    elseif (trim(smode) .eq. 'dual') then

      ! Check if group finite element formulation is applicable
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'dualbdrGFEM', dualBdrGFEM, 0)

      ! @FAQ2: What type of velocity are we?
      select case(abs(ivelocitytype))
        
      case default
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          ! The user-defined callback function is used if present;
          ! otherwise an error is throws
          if (present(fcb_coeffMatBdrDual1d_sim)) then
            call transp_calcBilfBdrCond1D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, fcb_coeffMatBdrDual1d_sim,&
                bclear, rmatrix, rcollection)
          else ! callback function not present
            call output_line('Missing user-defined callback function!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcBilfBdrCondition')
            call sys_halt()
          end if

        case (NDIM2D)
          ! The user-defined callback function is used if present;
          ! otherwise an error is throws
          if (present(fcb_coeffMatBdrDual2d_sim)) then
            call transp_calcBilfBdrCond2D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, fcb_coeffMatBdrDual2d_sim,&
                bclear, rmatrix, rcollection)
          else ! callback function not present
            call output_line('Missing user-defined callback function!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcBilfBdrCondition')
            call sys_halt()
          end if

        case (NDIM3D)
          ! The user-defined callback function is used if present;
          ! otherwise an error is throws
          if (present(fcb_coeffMatBdrDual3d_sim)) then
!!$            call transp_calcBilfBdrCond3D(rproblemLevel,&
!!$                rboundaryCondition, rsolution, ssectionName,&
!!$                dtime, dscale, fcb_coeffMatBdrDual3d_sim,&
!!$                bclear, rmatrix, rcollection)
          else ! callback function not present
            call output_line('Missing user-defined callback function!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcBilfBdrCondition')
            call sys_halt()
          end if
        end select

        !-----------------------------------------------------------------------
        
      case (VELOCITY_ZERO, VELOCITY_CONSTANT, VELOCITY_TIMEDEP)
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          ! linear velocity in 1D
          if (dualBdrGFEM .gt. 0) then
            call transp_calcOperatorBdrCondition(rproblemLevel,&
                rproblemLevel%RgroupFEMBlock(dualBdrGFEM),&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_calcMatBdrConvD1d_sim,&
                bclear, rmatrix, rcollection)
          else
            call transp_calcBilfBdrCond1D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_coeffMatBdrConvD1d_sim,&
                bclear, rmatrix, rcollection)
          end if

        case (NDIM2D)
          ! linear velocity in 2D
          if (dualBdrGFEM .gt. 0) then
            call transp_calcOperatorBdrCondition(rproblemLevel,&
                rproblemLevel%RgroupFEMBlock(dualBdrGFEM),&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_calcMatBdrConvD2d_sim,&
                bclear, rmatrix, rcollection)
          else
            call transp_calcBilfBdrCond2D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_coeffMatBdrConvD2d_sim,&
                bclear, rmatrix, rcollection)
          end if

        case (NDIM3D)
          ! linear velocity in 3D
          if (dualBdrGFEM .gt. 0) then
            call transp_calcOperatorBdrCondition(rproblemLevel,&
                rproblemLevel%RgroupFEMBlock(dualBdrGFEM),&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_calcMatBdrConvD3d_sim,&
                bclear, rmatrix, rcollection)
          else
!!$          call transp_calcBilfBdrCond3D(rproblemLevel,&
!!$              rboundaryCondition, rsolution, ssectionName,&
!!$              dtime, dscale, transp_coeffMatBdrConvD3d_sim,&
!!$              bclear, rmatrix, rcollection)
          end if
        end select

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers` equation in space-time
!!$        call transp_calcBilfBdrCond2D(rproblemLevel,&
!!$            rboundaryCondition, rsolution, ssectionName,&
!!$            dtime, dscale, transp_coeffMatBdrSTBurgersD2d_sim,&
!!$            bclear, rmatrix, , rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time
!!$        call transp_calcBilfBdrCond2D(rproblemLevel,&
!!$            rboundaryCondition, rsolution, ssectionName,&
!!$            dtime, dscale, transp_coeffMatBdrSTBLevD2d_sim,&
!!$            bclear, rmatrix, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers` equation in 1D
!!$        call transp_calcBilfBdrCond1D(rproblemLevel,&
!!$            rboundaryCondition, rsolution, ssectionName,&
!!$            dtime, dscale, transp_coeffMatBdrBurgersD1d_sim,&
!!$            bclear, rmatrix, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers` equation in 2D
!!$        call transp_calcBilfBdrCond2D(rproblemLevel,&
!!$            rboundaryCondition, rsolution, ssectionName,&
!!$            dtime, dscale, transp_coeffMatBdrBurgersD2d_sim,&
!!$            bclear, rmatrix, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D
!!$        call transp_calcBilfBdrCond1D(rproblemLevel,&
!!$            rboundaryCondition, rsolution, ssectionName,&
!!$            dtime, dscale, transp_coeffMatBdrBLevD1d_sim,&
!!$            bclear, rmatrix, rcollection)

      end select

    else
      call output_line('Invalid mode!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcBilfBdrCondition')
      call sys_halt()
    end if
    
  end subroutine transp_calcBilfBdrCondition

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcLinfBdrCondition(rproblemLevel,&
      rboundaryCondition, rsolution, ssectionName, smode,&
      ivelocitytype, dtime, dscale, bclear, rvector, rcollection,&
      fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
      fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
      fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim)

!<description>
    ! This subroutine is a shortcut for building the linear form
    ! arising from the weak imposition of boundary conditions. It
    ! calls the more general routine transp_calcLinfBdrCondXd
    ! using the corresponding callback routines.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! boundary condition
    type(t_boundaryCondition), intent(in) :: rboundaryCondition
    
    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! problem mode (primal, dual)
    character(LEN=*), intent(in) :: smode

    ! type of velocity
    integer, intent(in) :: ivelocitytype

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! Whether to clear the vector before calculating the entries.
    ! If .FALSE., the new vector entries are added to the existing entries.
    logical, intent(in) :: bclear

    ! user-defined callback functions
    include 'intf_transpCoeffVecBdr.inc'
    optional :: fcb_coeffVecBdrPrimal1d_sim
    optional :: fcb_coeffVecBdrPrimal2d_sim
    optional :: fcb_coeffVecBdrPrimal3d_sim
    optional :: fcb_coeffVecBdrDual1d_sim
    optional :: fcb_coeffVecBdrDual2d_sim
    optional :: fcb_coeffVecBdrDual3d_sim
!</intput>

!<inputoutput>
    ! vector where to store the linear form
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

    ! local variable
    type(t_parlist), pointer :: p_rparlist
    integer :: primalBdrGFEM,dualBdrGFEM
    
    
    ! Set pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)


    !###########################################################################
    ! REMARK: If you want to add a new type of velocity/diffusion,
    ! then search for the tag @FAQ2: in this subroutine and create a
    ! new CASE which performs the corresponding task for the new type
    ! of velocity/ diffusion.
    !###########################################################################

    ! Check if group finite element formulation is applicable
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'primalbdrGFEM', primalBdrGFEM, 0)

    ! Are we in primal or dual mode?
    if (trim(smode) .eq. 'primal') then
      
      ! @FAQ2: What type of velocity are we?
      select case(abs(ivelocitytype))
        
      case default
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          ! The user-defined callback function is used if present;
          ! otherwise an error is throws
          if (present(fcb_coeffVecBdrPrimal1d_sim)) then
            call transp_calcLinfBdrCond1D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, fcb_coeffVecBdrPrimal1d_sim,&
                bclear, rvector, rcollection)
          else ! callback function not present
            call output_line('Missing user-defined callback function!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcLinfBdrCondition')
            call sys_halt()
          end if
          
        case (NDIM2D)
          ! The user-defined callback function is used if present;
          ! otherwise an error is throws
          if (present(fcb_coeffVecBdrPrimal2d_sim)) then
            call transp_calcLinfBdrCond2D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, fcb_coeffVecBdrPrimal2d_sim,&
                bclear, rvector, rcollection)
          else ! callback function not present
            call output_line('Missing user-defined callback function!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcLinfBdrCondition')
            call sys_halt()
          end if
          
        case (NDIM3D)
          ! The user-defined callback function is used if present;
          ! otherwise an error is throws
          if (present(fcb_coeffVecBdrPrimal3d_sim)) then
!!$            call transp_calcLinfBdrCond3D(rproblemLevel,&
!!$                rboundaryCondition, rsolution, ssectionName,&
!!$                dtime, dscale, fcb_coeffVecBdrPrimal3d_sim,&
!!$                bclear, rvector, rcollection)
          else ! callback function not present
            call output_line('Missing user-defined callback function!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcLinfBdrCondition')
            call sys_halt()
          end if
        end select

        !-----------------------------------------------------------------------

      case (VELOCITY_ZERO, VELOCITY_CONSTANT, VELOCITY_TIMEDEP)
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          ! linear velocity in 1D
          if (primalBdrGFEM .gt. 0) then
            call transp_calcVectorBdrCondition(rproblemLevel,&
                rproblemLevel%RgroupFEMBlock(primalBdrGFEM),&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_calcVecBdrConvP1d_sim,&
                bclear, rvector, rcollection)
          else
            call transp_calcLinfBdrCond1D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_coeffVecBdrConvP1d_sim,&
                bclear, rvector, rcollection)
          end if

        case (NDIM2D)
          ! linear velocity in 2D
          if (primalBdrGFEM .gt. 0) then
            call transp_calcVectorBdrCondition(rproblemLevel,&
                rproblemLevel%RgroupFEMBlock(primalBdrGFEM),&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_calcVecBdrConvP2d_sim,&
                bclear, rvector, rcollection)
          else
            call transp_calcLinfBdrCond2D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_coeffVecBdrConvP2d_sim,&
                bclear, rvector, rcollection)
          end if
          
        case (NDIM3D)
          ! linear velocity in 3D
          if (primalBdrGFEM .gt. 0) then
            call transp_calcVectorBdrCondition(rproblemLevel,&
                rproblemLevel%RgroupFEMBlock(primalBdrGFEM),&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_calcVecBdrConvP3d_sim,&
                bclear, rvector, rcollection)
          else
!!$          call transp_calcLinfBdrCond3D(rproblemLevel,&
!!$              rboundaryCondition, rsolution, ssectionName,&
!!$              dtime, dscale, transp_coeffVecBdrConvP3d_sim,&
!!$              bclear, rvector, rcollection)
          end if
        end select

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers` equation in space-time
        call transp_calcLinfBdrCond2D(rproblemLevel,&
            rboundaryCondition, rsolution, ssectionName,&
            dtime, dscale, transp_coeffVecBdrSTBurgP2d_sim,&
            bclear, rvector, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time
        call transp_calcLinfBdrCond2D(rproblemLevel,&
            rboundaryCondition, rsolution, ssectionName,&
            dtime, dscale, transp_coeffVecBdrSTBLevP2d_sim,&
            bclear, rvector, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers` equation in 1D
        call transp_calcLinfBdrCond1D(rproblemLevel,&
            rboundaryCondition, rsolution, ssectionName,&
            dtime, dscale, transp_coeffVecBdrBurgP1d_sim,&
            bclear, rvector, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers` equation in 2D
        call transp_calcLinfBdrCond2D(rproblemLevel,&
            rboundaryCondition, rsolution, ssectionName,&
            dtime, dscale, transp_coeffVecBdrBurgP2d_sim,&
            bclear, rvector, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D
        call transp_calcLinfBdrCond1D(rproblemLevel,&
            rboundaryCondition, rsolution, ssectionName,&
            dtime, dscale, transp_coeffVecBdrBLevP1d_sim,&
            bclear, rvector, rcollection)

      end select

      !-------------------------------------------------------------------------

    elseif (trim(smode) .eq. 'dual') then

      ! Check if group finite element formulation is applicable
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'dualbdrGFEM', dualBdrGFEM, 0)

      ! @FAQ2: What type of velocity are we?
      select case(abs(ivelocitytype))
        
      case default
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          ! The user-defined callback function is used if present;
          ! otherwise an error is throws
          if (present(fcb_coeffVecBdrDual1d_sim)) then
            call transp_calcLinfBdrCond1D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, fcb_coeffVecBdrDual1d_sim,&
                bclear, rvector, rcollection)
          else ! callback function not present
            call output_line('Missing user-defined callback function!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcLinfBdrCondition')
            call sys_halt()
          end if
          
        case (NDIM2D)
          ! The user-defined callback function is used if present;
          ! otherwise an error is throws
          if (present(fcb_coeffVecBdrDual2d_sim)) then
            call transp_calcLinfBdrCond2D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, fcb_coeffVecBdrDual2d_sim,&
                bclear, rvector, rcollection)
          else ! callback function not present
            call output_line('Missing user-defined callback function!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcLinfBdrCondition')
            call sys_halt()
          end if

        case (NDIM3D)
          ! The user-defined callback function is used if present;
          ! otherwise an error is throws
          if (present(fcb_coeffVecBdrDual3d_sim)) then
!!$            call transp_calcLinfBdrCond3D(rproblemLevel,&
!!$                rboundaryCondition, rsolution, ssectionName,&
!!$                dtime, dscale, fcb_coeffVecBdrDual3d_sim,&
!!$                bclear, rvector, rcollection)
          else ! callback function not present
            call output_line('Missing user-defined callback function!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcLinfBdrCondition')
            call sys_halt()
          end if
        end select

        !-----------------------------------------------------------------------

      case (VELOCITY_ZERO, VELOCITY_CONSTANT, VELOCITY_TIMEDEP)
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          ! linear velocity in 1D
          if (dualBdrGFEM .gt. 0) then
            call transp_calcVectorBdrCondition(rproblemLevel,&
                rproblemLevel%RgroupFEMBlock(dualBdrGFEM),&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_calcVecBdrConvD1d_sim,&
                bclear, rvector, rcollection)
          else
            call transp_calcLinfBdrCond1D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_coeffVecBdrConvD1d_sim,&
                bclear, rvector, rcollection)
          end if

        case (NDIM2D)
          ! linear velocity in 2D
          if (dualBdrGFEM .gt. 0) then
            call transp_calcVectorBdrCondition(rproblemLevel,&
                rproblemLevel%RgroupFEMBlock(dualBdrGFEM),&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_calcVecBdrConvD2d_sim,&
                bclear, rvector, rcollection)
          else
            call transp_calcLinfBdrCond2D(rproblemLevel,&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_coeffVecBdrConvD2d_sim,&
                bclear, rvector, rcollection)
          end if

        case (NDIM3D)
          ! linear velocity in 3D
          if (dualBdrGFEM .gt. 0) then
            call transp_calcVectorBdrCondition(rproblemLevel,&
                rproblemLevel%RgroupFEMBlock(dualBdrGFEM),&
                rboundaryCondition, rsolution, ssectionName,&
                dtime, dscale, transp_calcVecBdrConvD3d_sim,&
                bclear, rvector, rcollection)
          else
!!$          call transp_calcLinfBdrCond3D(rproblemLevel,&
!!$              rboundaryCondition, rsolution, ssectionName,&
!!$              dtime, dscale, transp_coeffVecBdrConvD3d_sim,&
!!$              bclear, rvector, rcollection)
          end if
        end select

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers` equation in space-time
!!$        call transp_calcLinfBdrCond2D(rproblemLevel,&
!!$            rboundaryCondition, rsolution, ssectionName,&
!!$            dtime, dscale, transp_coeffVecBdrSTBurgersD2d_sim,&
!!$            bclear, rvector, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time
!!$        call transp_calcLinfBdrCond2D(rproblemLevel,&
!!$            rboundaryCondition, rsolution, ssectionName,&
!!$            dtime, dscale, transp_coeffVecBdrSTBLevD2d_sim,&
!!$            bclear, rvector, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers` equation in 1D
!!$        call transp_calcLinfBdrCond1D(rproblemLevel,&
!!$            rboundaryCondition, rsolution, ssectionName,&
!!$            dtime, dscale, transp_coeffVecBdrBurgersD1d_sim,&
!!$            bclear, rvector, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers` equation in 2D
!!$        call transp_calcLinfBdrCond2D(rproblemLevel,&
!!$            rboundaryCondition, rsolution, ssectionName,&
!!$            dtime, dscale, transp_coeffVecBdrBurgersD2d_sim,&
!!$            bclear, rvector, rcollection)

        !-----------------------------------------------------------------------

      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D
!!$        call transp_calcLinfBdrCond1D(rproblemLevel,&
!!$            rboundaryCondition, rsolution, ssectionName,&
!!$            dtime, dscale, transp_coeffVecBdrBLevD1d_sim,&
!!$            bclear, rvector, rcollection)

      end select

    else
      call output_line('Invalid mode!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcLinfBdrCondition')
      call sys_halt()
    end if

  end subroutine transp_calcLinfBdrCondition

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcOperatorBdrCondition(rproblemLevel,&
      rgroupFEMBlock, rboundaryCondition, rsolution, ssectionName,&
      dtime, dscale, fcb_calcMatBdr_sim, bclear, rmatrix, rcollection)
    
!<description>
    ! This subroutine calculates the discrete operator at the boundary
    ! using the group finite element formulation.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! block of group finite element set
    type(t_groupFEMBlock), intent(IN) :: rgroupFEMBlock

    ! boundary condition
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! Whether to clear the matrix before calculating the entries.
    ! If .FALSE., the new matrix entries are added to the existing entries.
    logical, intent(in) :: bclear
    
    ! user-defined callback functions
    include 'intf_transpCalcMatBdr.inc'
!</intput>

!<inputoutput>
    ! scalar matrix where to store the bilinear form
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_collection) :: rcollectionTmp
    integer, dimension(:), pointer :: p_IbdrCondCpIdx,p_IbdrCondType
    integer :: ibdc,isegment,ivelocitytype,velocityfield,dofcoords
    

    ! Return if there are no weak boundary conditions
    if (.not.rboundaryCondition%bWeakBdrCond) return
    
    ! Set pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    
    ! Attach user-defined collection structure to temporal collection
    ! structure (may be required by the callback function)
    rcollectionTmp%p_rnextCollection => rcollection
    
    ! Prepare quick access arrays of temporal collection structure
    rcollectionTmp%DquickAccess(1) = dtime
    
    ! Set coordinates of the degrees of freedom
    call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'dofCoords', dofCoords, 0)
    if (dofCoords .gt. 0) then
      rcollectionTmp%p_rvectorQuickAccess1 => rproblemLevel%RvectorBlock(dofCoords)
    else
      nullify(rcollectionTmp%p_rvectorQuickAccess1)
    end if

    ! Set vector for velocity field (if any)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    if (transp_hasVelocityVector(ivelocityType)) then
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'velocityfield', velocityfield)
      rcollectionTmp%p_rvectorQuickAccess2 => rproblemLevel%RvectorBlock(velocityfield)
    else
      nullify(rcollectionTmp%p_rvectorQuickAccess2)
    end if

    ! Clear matrix?
    if (bclear) call lsyssc_clearMatrix(rmatrix)

    ! Set pointers
    call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx,&
        p_IbdrCondCpIdx)
    call storage_getbase_int(rboundaryCondition%h_IbdrCondType,&
        p_IbdrCondType)
    
    ! Loop over all boundary components
    do ibdc = 1, rboundaryCondition%iboundarycount
      
      ! Loop over all boundary segments
      do isegment = p_IbdrCondCpIdx(ibdc), p_IbdrCondCpIdx(ibdc+1)-1
        
        ! Check if this segment has weak boundary conditions
        if (iand(p_IbdrCondType(isegment), BDRC_WEAK) .ne. BDRC_WEAK) cycle

        ! Prepare further quick access arrays of temporal collection
        ! structure with boundary component and type
        rcollectionTmp%IquickAccess(1) = p_IbdrCondType(isegment)
        rcollectionTmp%IquickAccess(2) = isegment
        
        ! What type of boundary conditions are we?
        select case(iand(p_IbdrCondType(isegment), BDRC_TYPEMASK))
          
        case (BDRC_ROBIN)
          ! Do nothing since boundary conditions are build into the
          ! linear form and the bilinear form has no boundary term
          
        case (BDRC_HOMNEUMANN, BDRC_INHOMNEUMANN,&
              BDRC_FLUX, BDRC_PERIODIC, BDRC_ANTIPERIODIC)
          
          ! Assemble the discrete operator consistently
          call gfsc_buildOperatorNode(rgroupFEMBlock%RgroupFEMBlock(isegment),&
              rsolution, fcb_calcMatBdr_sim, dscale, .false., rmatrix,&
              cconstrType=GFEM_MATC_CONSISTENT, rcollection=rcollectionTmp)
          
        case (BDRC_DIRICHLET)
          ! Assemble the discrete operator with lumping
          call gfsc_buildOperatorNode(rgroupFEMBlock%RgroupFEMBlock(isegment),&
              rsolution, fcb_calcMatBdr_sim, dscale, .false., rmatrix,&
              cconstrType=GFEM_MATC_LUMPED, rcollection=rcollectionTmp)
          
        case default
          call output_line('Unsupported type of boundary conditions!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcOperatorBdrCondition')
          call sys_halt()
          
        end select
        
      end do ! isegment
      
    end do ! ibdc

  end subroutine transp_calcOperatorBdrCondition

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcVectorBdrCondition(rproblemLevel, rgroupFEMBlock,&
      rboundaryCondition, rsolution, ssectionName, dtime, dscale,&
      fcb_calcVecBdr_sim, bclear, rvector, rcollection)

!<description>
    ! This subroutine calculates the discrete operator at the boundary
    ! using the group finite element formulation.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! block of group finite element sets
    type(t_groupFEMBlock), intent(IN) :: rgroupFEMBlock

    ! boundary condition
    type(t_boundaryCondition), intent(in) :: rboundaryCondition

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! Whether to clear the vector before calculating the entries.
    ! If .FALSE., the new vector entries are added to the existing entries.
    logical, intent(in) :: bclear
    
    ! user-defined callback functions
    include 'intf_transpCalcVecBdr.inc'
!</intput>

!<inputoutput>
    ! vector where to store the linear form
    type(t_vectorBlock), intent(inout) :: rvector

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_collection) :: rcollectionTmp
    integer, dimension(:), pointer :: p_IbdrCondCpIdx, p_IbdrCondType
    integer, dimension(:), pointer :: p_IbdrCompPeriodic, p_IbdrCondPeriodic
    integer :: ibdc,isegment,velocityfield,dofcoords,ivelocitytype
    

    ! Return if there are no weak boundary conditions
    if (.not.rboundaryCondition%bWeakBdrCond) return
    
    ! Set pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)

    ! Initialise temporal collection structure
    call collct_init(rcollectionTmp)

    ! Attach user-defined collection structure to temporal collection
    ! structure (may be required by the callback function)
    rcollectionTmp%p_rnextCollection => rcollection

    ! Attach function parser from boundary conditions to collection
    ! structure and specify its name in quick access string array
    call collct_setvalue_pars(rcollectionTmp, 'rfparser',&
        rboundaryCondition%rfparser, .true.)

    ! Prepare quick access arrays of temporal collection structure
    rcollectionTmp%SquickAccess(1) = ''
    rcollectionTmp%SquickAccess(2) = 'rfparser'
    rcollectionTmp%DquickAccess(1) = dtime

    ! Set coordinates of the degrees of freedom
    call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'dofCoords', dofCoords, 0)
    if (dofCoords .gt. 0) then
      rcollectionTmp%p_rvectorQuickAccess1 => rproblemLevel%RvectorBlock(dofCoords)
    else
      nullify(rcollectionTmp%p_rvectorQuickAccess1)
    end if

    ! Set vector for velocity field (if any)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    if (transp_hasVelocityVector(ivelocityType)) then
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'velocityfield', velocityfield)
      rcollectionTmp%p_rvectorQuickAccess2 => rproblemLevel%RvectorBlock(velocityfield)
    else
      nullify(rcollectionTmp%p_rvectorQuickAccess2)
    end if

    ! Clear vector?
    if (bclear) call lsysbl_clearVector(rvector)

    ! Set pointers
    call storage_getbase_int(rboundaryCondition%h_IbdrCondCpIdx,&
        p_IbdrCondCpIdx)
    call storage_getbase_int(rboundaryCondition%h_IbdrCondType,&
        p_IbdrCondType)

    ! Set additional pointers for periodic boundary conditions
    if (rboundaryCondition%bPeriodic) then
      call storage_getbase_int(rboundaryCondition%h_IbdrCompPeriodic,&
          p_IbdrCompPeriodic)
      call storage_getbase_int(rboundaryCondition%h_IbdrCondPeriodic,&
          p_IbdrCondPeriodic)
    end if
    
    ! Loop over all boundary components
    do ibdc = 1, rboundaryCondition%iboundarycount
      
      ! Loop over all boundary segments
      do isegment = p_IbdrCondCpIdx(ibdc), p_IbdrCondCpIdx(ibdc+1)-1
        
        ! Check if this segment has weak boundary conditions
        if (iand(p_IbdrCondType(isegment), BDRC_WEAK) .ne. BDRC_WEAK) cycle
        
        ! Prepare further quick access arrays of temporal collection
        ! structure with boundary component and type
        rcollectionTmp%IquickAccess(1) = p_IbdrCondType(isegment)
        rcollectionTmp%IquickAccess(2) = isegment
        
        ! What type of boundary conditions are we?
        select case(iand(p_IbdrCondType(isegment), BDRC_TYPEMASK))
          
        case (BDRC_HOMNEUMANN)
          ! Do nothing for homogeneous Neumann boundary conditions
          ! since the boundary integral vanishes by construction
          
        case (BDRC_INHOMNEUMANN, BDRC_ROBIN,&
              BDRC_FLUX, BDRC_DIRICHLET,&
              BDRC_PERIODIC, BDRC_ANTIPERIODIC)

          ! Check if special treatment of mirror boundary condition is required
          if ((iand(p_IbdrCondType(isegment), BDRC_TYPEMASK) .eq. BDRC_PERIODIC) .or.&
              (iand(p_IbdrCondType(isegment), BDRC_TYPEMASK) .eq. BDRC_ANTIPERIODIC)) then
            call output_line('Mirror boundaries are not yet supported!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcVectorBdrCondition')
            call sys_halt()
          end if

          ! Assemble the discrete vector
          call gfsc_buildVectorNode(rgroupFEMBlock%RgroupFEMBlock(isegment),&
              rsolution, fcb_calcVecBdr_sim, dscale, .false., rvector,&
              rcollectionTmp)
          
        case default
          call output_line('Unsupported type of boundary copnditions!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcVectorBdrCondition')
          call sys_halt()
          
        end select
        
      end do ! isegment
      
    end do ! ibdc
    
    ! Release temporal collection structure
    call collct_done(rcollectionTmp)
    
  end subroutine transp_calcVectorBdrCondition

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcVelocityField(rparlist, ssectionName,&
      rproblemLevel, dtime, rcollection, nlminOpt)
    
!<description>
    ! This subroutine calculates the velocity fields from the function
    ! parser. The result is stored separately for each problem level.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! simulation time
    real(DP), intent(in) :: dtime

    ! OPTIONAL: minimum problem level
    integer, intent(in), optional :: nlminOpt
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout), target :: rproblemLevel

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_problemLevel), pointer :: p_rproblemLevel
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
    real(DP), dimension(:), pointer :: p_Ddata,p_DdofCoords
    character(LEN=SYS_STRLEN) :: svelocityname
    integer :: neq,idim,ndim,nlmin
    integer :: ivelocitytype,velocityfield,discretisation,dofCoords

    ! Get parameters from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype, VELOCITY_ZERO)

    select case(ivelocitytype)
    case (VELOCITY_CONSTANT, VELOCITY_TIMEDEP)

      ! Get parameter from parameter list
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'velocityfield', velocityfield, 0)
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'discretisation', discretisation, 0)
      call parlst_getvalue_int(rparlist,&
          ssectionName, 'dofCoords', dofCoords, 0)
      
      ! Get function parser from collection
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)
      
      ! Set minimum problem level
      nlmin = rproblemLevel%ilev
      if (present(nlminOpt)) nlmin = nlminOpt
      
      ! Loop over all problem levels
      p_rproblemLevel => rproblemLevel
      do while(associated(p_rproblemLevel))
        
        if ((velocityfield  .gt. 0) .and.&
            (discretisation .gt. 0) .and.&
            (dofCoords      .gt. 0)) then
          
          ! Get number of degrees of freedom and spatial dimension
          p_rspatialDiscr =>&
              p_rproblemLevel%RblockDiscretisation(discretisation)%RspatialDiscr(1)
          neq  = dof_igetNDofGlob(p_rspatialDiscr)
          ndim = p_rspatialDiscr%ndimension
          
          ! Create/resize velocity vector if required
          if (p_rproblemLevel%RvectorBlock(velocityfield)%NEQ .eq. 0) then
            call lsysbl_createVectorBlock(&
                p_rproblemLevel%rvectorBlock(velocityfield), neq, ndim, .true.)
          elseif (p_rproblemLevel%RvectorBlock(velocityfield)%NEQ .ne. neq*ndim) then
            call lsysbl_resizeVectorBlock(&
                p_rproblemLevel%rvectorBlock(velocityfield), neq, .true.)
          end if
          
          ! Get coordinates of the DOF`s of the current problem level
          call lsysbl_getbase_double(&
              p_rproblemLevel%RvectorBlock(dofCoords), p_DdofCoords)
          
          ! Loop over all spatial dimensions
          do idim = 1, ndim
            
            ! Attach discretisation structure
            p_rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(idim)&
                %p_rspatialDiscr => p_rspatialDiscr

            ! Get scalar subvector
            call lsyssc_getbase_double(&
                p_rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(idim), p_Ddata)
            
            ! Retrieve function name from parameter list
            call parlst_getvalue_string(rparlist, ssectionName,&
                'svelocityname', svelocityname, isubString=idim)
            
            ! Evaluate all coefficients of the scalar subvector
            call fparser_evalFuncBlDbleFixName(p_rfparser, svelocityname,&
                ndim, neq, p_DdofCoords, neq, p_Ddata, (/dtime/))
          end do
          
        else
          call output_line('Coordinates of DOFs not available!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcVelocityField')
          call sys_halt()
        end if
        
        ! Set update notifier for transport operator in problem level structure
        p_rproblemLevel%iproblemSpec = ior(p_rproblemLevel%iproblemSpec,&
                                           TRANSP_TROPER_UPDATE)
        
        ! Proceed to coarser problem level if minimum level has not been reached
        if (p_rproblemLevel%ilev .le. nlmin) exit
        p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      end do

    case default
      return
    end select
    
  end subroutine transp_calcVelocityField

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcLinearisedFCT(rproblemLevel, rtimestep,&
      rsolver, rsolution, ssectionName, rcollection, rsource,&
      fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
      fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
      fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim,&
      rmatrix, rvector1, rvector2, rvector3)

!<description>
    ! This subroutine calculates the linearised FCT correction
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: source vector
    type(t_vectorBlock), intent(in), optional :: rsource

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCoeffVecBdr.inc'
    optional :: fcb_coeffVecBdrPrimal1d_sim
    optional :: fcb_coeffVecBdrPrimal2d_sim
    optional :: fcb_coeffVecBdrPrimal3d_sim
    optional :: fcb_coeffVecBdrDual1d_sim
    optional :: fcb_coeffVecBdrDual2d_sim
    optional :: fcb_coeffVecBdrDual3d_sim
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection

    ! OPTIONAL: matrix to store the discrete transport operator used
    ! to compute the approximation to the time derivative (if not
    ! present, then temporal memory is allocated)
    type(t_matrixScalar), intent(inout), target, optional :: rmatrix

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
    integer :: convectionAFC,imassantidiffusiontype
    integer :: lumpedMassMatrix,consistentMassMatrix


    ! Get parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'convectionAFC', convectionAFC, 0)
    
    !---------------------------------------------------------------------------
    ! Linearised FEM-FCT algorithm for the convective term (if any)
    !---------------------------------------------------------------------------
    
    if (convectionAFC .gt. 0) then

      ! What type of stabilisation are we?
      select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

      case (AFCSTAB_LINFCT)
        ! Get parameters from parameter list
        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'imassantidiffusiontype', imassantidiffusiontype)
        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'consistentmassmatrix', consistentmassmatrix)
        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'lumpedmassmatrix', lumpedmassmatrix)
        
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          
          ! Set up vector for computing the approximate time derivative
          if (present(rvector1)) then
            p_rvector1 => rvector1
          else
            allocate(p_rvector1)
          end if

          ! Compute approximate time derivative
          call transp_calcTimeDerivative(rproblemLevel, rtimestep, rsolver,&
              rsolution, ssectionName, rcollection, p_rvector1, rsource,&
              fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
              fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
              fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim,&
              rmatrix, rvector2, rvector3)
          
          ! Build the raw antidiffusive fluxes and include
          ! contribution from the consistent mass matrix
          call afcsc_buildFluxFCT(&
              rproblemLevel%Rafcstab(convectionAFC),&
              rsolution, 1.0_DP, 0.0_DP, 1.0_DP, 1.0_DP,&
              .true., .true., AFCSTAB_FCTFLUX_EXPLICIT,&
              rmatrix=rproblemLevel%RmatrixScalar(consistentMassMatrix),&
              rxTimeDeriv=p_rvector1)
          
          ! Release temporal memory
          if (.not.present(rvector1)) then
            call lsysbl_releaseVector(p_rvector1)
            deallocate(p_rvector1)
          end if
          
        else
          
          ! Build the raw antidiffusive fluxes without including
          ! the contribution from consistent mass matrix
          call afcsc_buildFluxFCT(&
              rproblemLevel%Rafcstab(convectionAFC),&
              rsolution, 1.0_DP, 0.0_DP, 1.0_DP, 1.0_DP,&
              .true., .true., AFCSTAB_FCTFLUX_EXPLICIT)
        end if
        
        ! Apply linearised FEM-FCT correction
        call afcsc_buildVectorFCT(&
            rproblemLevel%Rafcstab(convectionAFC),&
            rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
            rsolution, rtimestep%dStep, .false.,&
            AFCSTAB_FCTALGO_STANDARD+&
            AFCSTAB_FCTALGO_SCALEBYMASS, rsolution)
        
      end select
    end if
    
    ! Impose boundary conditions for the solution vector
    call bdrf_filterVectorExplicit(rsolver%rboundaryCondition,&
        rsolution, rtimestep%dTime)
    
  end subroutine transp_calcLinearisedFCT

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcLinearisedLPT(rproblemLevel, rtimestep,&
      rsolver, rsolution, ssectionName, rcollection, rsource,&
      fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
      fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
      fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim,&
      rmatrix, rvector1, rvector2, rvector3)

!<description>
    ! This subroutine calculates the linearised correction from
    ! linearity-preserving flux correction
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: source vector
    type(t_vectorBlock), intent(in), optional :: rsource

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCoeffVecBdr.inc'
    optional :: fcb_coeffVecBdrPrimal1d_sim
    optional :: fcb_coeffVecBdrPrimal2d_sim
    optional :: fcb_coeffVecBdrPrimal3d_sim
    optional :: fcb_coeffVecBdrDual1d_sim
    optional :: fcb_coeffVecBdrDual2d_sim
    optional :: fcb_coeffVecBdrDual3d_sim
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection

    ! OPTIONAL: matrix to store the discrete transport operator used
    ! to compute the approximation to the time derivative (if not
    ! present, then temporal memory is allocated)
    type(t_matrixScalar), intent(inout), target, optional :: rmatrix

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
    integer :: massAFC,convectionAFC,diffusionAFC
    integer :: lumpedMassMatrix,consistentMassMatrix

    ! Get parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)
    
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'massAFC', massAFC, 0)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'convectionAFC', convectionAFC, 0)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'diffusionAFC', diffusionAFC, 0)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'lumpedMassMatrix', lumpedMassMatrix)

    !---------------------------------------------------------------------------
    ! Linearised FEM-LPT algorithm for the convective term (if any)
    !---------------------------------------------------------------------------
    
    if (convectionAFC .gt. 0) then

      ! What type of stabilisation are we?
      select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
        
      case (AFCSTAB_LINLPT_UPWINDBIASED)
        
        ! Check if slope-based local-extremum diminishing bounds
        ! at edges are available and compute them otherwise
        if (iand(rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec,&
            AFCSTAB_HAS_EDGEBOUNDS) .eq. 0)&
            call transp_calcBoundsAtEdgesLPT(p_rparlist, ssectionName,&
            rproblemLevel, rproblemLevel%Rafcstab(convectionAFC))
        
        ! Build the raw antidiffusive fluxes
        call afcsc_buildFluxLPT(rproblemLevel%Rafcstab(convectionAFC),&
            rsolution, 1.0_DP, .true., AFCSTAB_LPTFLUX+AFCSTAB_LPTFLUX_BOUNDS)
        
        ! Apply linearised FEM-LPT correction
        call afcsc_buildVectorLPT(&
          rproblemLevel%Rafcstab(convectionAFC),&
          rsolution, rtimestep%dStep, .false.,&
          AFCSTAB_LPTALGO_STANDARD + AFCSTAB_LPTALGO_SCALEBYMASS,&
          rsolution, rproblemLevel%RmatrixScalar(lumpedMassMatrix))
      end select
    end if
    
    !---------------------------------------------------------------------------
    ! Linearised FEM-LPT algorithm for the diffusion term (if any)
    !---------------------------------------------------------------------------
    
    if (diffusionAFC .gt. 0) then
      
      ! What type of stabilisation are we?
      select case(rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType)
        
      case (AFCSTAB_LINLPT_SYMMETRIC)

        ! Check if slope-based local-extremum diminishing bounds
        ! at edges are available and compute them otherwise
        if (iand(rproblemLevel%Rafcstab(diffusionAFC)%istabilisationSpec,&
            AFCSTAB_HAS_EDGEBOUNDS) .eq. 0)&
            call transp_calcBoundsAtEdgesLPT(p_rparlist, ssectionName,&
            rproblemLevel, rproblemLevel%Rafcstab(diffusionAFC))
        
        ! Build the raw antidiffusive fluxes
        call afcsc_buildFluxLPT(rproblemLevel%Rafcstab(diffusionAFC),&
            rsolution, 1.0_DP, .true., AFCSTAB_LPTFLUX+AFCSTAB_LPTFLUX_BOUNDS)
        
        ! Apply linearised FEM-LPT correction
        call afcsc_buildVectorLPT(&
            rproblemLevel%Rafcstab(diffusionAFC),&
            rsolution, rtimestep%dStep, .false.,&
            AFCSTAB_LPTALGO_STANDARD + AFCSTAB_LPTALGO_SCALEBYMASS,&
            rsolution, rproblemLevel%RmatrixScalar(lumpedMassMatrix))
      end select
    end if

    !---------------------------------------------------------------------------
    ! Linearised FEM-LPT algorithm for the mass term (if any)
    !---------------------------------------------------------------------------

    if (massAFC .gt. 0) then
      
      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(massAFC)%cafcstabType)
        
      case (AFCSTAB_LINLPT_MASS)
        
        call parlst_getvalue_int(p_rparlist,&
            ssectionName, 'consistentMassMatrix', consistentMassMatrix)
        
        ! Check if slope-based local-extremum diminishing bounds
        ! at edges are available and compute them otherwise
        if (iand(rproblemLevel%Rafcstab(massAFC)%istabilisationSpec,&
            AFCSTAB_HAS_EDGEBOUNDS) .eq. 0)&
            call transp_calcBoundsAtEdgesLPT(p_rparlist, ssectionName,&
            rproblemLevel, rproblemLevel%Rafcstab(massAFC))
        
        ! Set up vector for computing the approximate time derivative
        if (present(rvector1)) then
          p_rvector1 => rvector1
        else
          allocate(p_rvector1)
        end if
        
        ! Compute approximate time derivative
        call transp_calcTimeDerivative(rproblemLevel, rtimestep, rsolver,&
            rsolution, ssectionName, rcollection, p_rvector1, rsource,&
            fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
            fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
            fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim,&
            rmatrix, rvector2, rvector3)

        ! Build the raw antidiffusive fluxes based on the approximate
        ! time derivative stored in vector p_rvector1
        call afcsc_buildFluxLPT(rproblemLevel%Rafcstab(massAFC),&
            p_rvector1, 1.0_DP, .true., AFCSTAB_LPTFLUX+AFCSTAB_LPTFLUX_BOUNDS,&
            rproblemLevel%RmatrixScalar(consistentMassMatrix))

        ! Apply linearised FEM-LPT correction
        call afcsc_buildVectorLPT(&
          rproblemLevel%Rafcstab(massAFC),&
          p_rvector1, rtimestep%dStep, .false.,&
          AFCSTAB_LPTALGO_STANDARD + AFCSTAB_LPTALGO_SCALEBYMASS,&
          rsolution, rproblemLevel%RmatrixScalar(lumpedMassMatrix))
      
        ! Release temporal memory
        if (.not.present(rvector1)) then
          call lsysbl_releaseVector(p_rvector1)
          deallocate(p_rvector1)
        end if
      end select
    end if
    
    ! Impose boundary conditions for the solution vector
    call bdrf_filterVectorExplicit(rsolver%rboundaryCondition,&
        rsolution, rtimestep%dTime)
          
  end subroutine transp_calcLinearisedLPT

  ! ***************************************************************************

!<subroutine>

  subroutine transp_calcBoundsAtEdgesLPT(rparlist, ssectionName,&
      rproblemLevel, rafcstab)

!<description>
    ! This subroutine calculates the LED bounds at the edges for
    ! linearity-preserving flux correction
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(len=*), intent(in) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel
!</input>

!<inputoutput>
    ! stabilisation structure
    type(t_afcstab), intent(inout) :: rafcstab
!</inputoutput>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), pointer :: p_DdofCoords
    integer :: lumpedmassmatrix,dofCoords
    integer :: coeffMatrix_Cx,coeffMatrix_Cy,coeffMatrix_Cz
    
    ! Get parameters from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'dofCoords', dofCoords, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedmassmatrix, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cx', coeffMatrix_Cx, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cy', coeffMatrix_Cy, 0)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'coeffmatrix_cz', coeffMatrix_Cz, 0)
    
    ! Get coordinates of the degrees of freedom
    if (dofCoords .gt. 0) then
      call lsysbl_getbase_double(&
          rproblemLevel%RvectorBlock(dofCoords), p_DdofCoords)
    else
      call output_line('Problem does not provide coordinates of DOFs!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_calcBoundsAtEdgesLPT')
      call sys_halt()
    end if

    ! Compute bounds for linearity-preserving limiter
    select case(rproblemLevel%rtriangulation%ndim)
      
    case(NDIM1D)
      call afcstab_buildBoundsLPT(rafcstab,&
          rproblemLevel%RmatrixScalar(coeffMatrix_CX),&
          rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
          p_DdofCoords, 2.0_DP)
      
    case(NDIM2D)
      call afcstab_buildBoundsLPT(rafcstab,&
          rproblemLevel%RmatrixScalar(coeffMatrix_CX),&
          rproblemLevel%RmatrixScalar(coeffMatrix_CY),&
          rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
          p_DdofCoords, 2.0_DP)
      
    case(NDIM3D)
      call afcstab_buildBoundsLPT(rafcstab,&
          rproblemLevel%RmatrixScalar(coeffMatrix_CX),&
          rproblemLevel%RmatrixScalar(coeffMatrix_CY),&
          rproblemLevel%RmatrixScalar(coeffMatrix_CZ),&
          rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
          p_DdofCoords, 2.0_DP)
      
    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_calcBoundsAtEdgesLPT')
      call sys_halt()
    end select
    
    ! Set specifier for bounds at edges
    rafcstab%istabilisationSpec =&
        ior(rafcstab%istabilisationSpec, AFCSTAB_HAS_EDGEBOUNDS)
    
  end subroutine transp_calcBoundsAtEdgesLPT

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVectorFE(rdiscretisation, rform,&
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

    ! An array accepting the DOF`s on all elements in the test space.
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
    integer :: iel,ipoint
    
    if (.not. present(rcollection)) then
      Dcoefficients(:,:,:) = 0.0_DP
      return
    end if

    ! Get the parameters from the collection
    p_rvector => rcollection%p_rvectorQuickAccess1

    ! Allocate temporal array
    allocate(Ddata(npointsPerElement))
    
    ! Loop over all elements
    do iel = 1, nelements
      
      ! Evaluate solution in cubature points
      call fevl_evaluate(DER_FUNC, Ddata,&
          p_rvector%RvectorBlock(1), Dpoints(:,:,iel))
      
      ! Loop over all cubature points
      do ipoint = 1, npointsPerElement
        Dcoefficients(1,ipoint,iel) = Ddata(ipoint)
      end do
    end do
    
    ! Deallocate temporal array
    deallocate(Ddata)

  end subroutine transp_coeffVectorFE

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVectorAnalytic(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

!<description>
    ! This subroutine is called during the vector assembly. It has to
    ! compute the coefficients in front of the terms of the linear
    ! form. This routine can be used universaly for arbitrary linear
    ! forms for which the coefficients are evaluated analytically
    ! using a function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
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

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   DquickAccess(1):            simulation time
    !   DquickAccess(1:ntermCount): number of the analytic function(s)
    !   SquickAccess(1):            section name in the collection
    !   SquickAccess(2):            string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</inputoutput>

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
    real(DP), dimension(1) :: dtime
    integer :: itermCount, iel, icomp
    
    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))
    
    ! This subroutine assumes that the first quick access double
    ! value holds the simulation time
    dtime(1) = rcollection%DquickAccess(1)

    ! Loop over all components of the linear form
    do itermCount = 1, ubound(Dcoefficients,1)

      ! Moreover, this subroutine assumes the quick access integer
      ! 'itermCount' holds the number of the function to be evaluated
      icomp = rcollection%IquickAccess(itermCount)

      if (dtime(1) .lt. 0.0_DP) then

        ! Evaluate all coefficients using the function parser
        !$omp parallel do
        do iel = 1, nelements
          call fparser_evalFunction(p_rfparser, icomp, 2,&
              Dpoints(:,:,iel), Dcoefficients(itermCount,:,iel))
        end do
        !$omp end parallel do

      else

        ! Evaluate all coefficients using the function parser
        !$omp parallel do
        do iel = 1, nelements
          call fparser_evalFunction(p_rfparser, icomp, 2,&
              Dpoints(:,:,iel), Dcoefficients(itermCount,:,iel), dtime)
        end do
        !$omp end parallel do

      end if

    end do ! itermCount

  end subroutine transp_coeffVectorAnalytic

  !*****************************************************************************

!<subroutine>

  subroutine transp_refFuncAnalytic(cderivative, rdiscretisation,&
      nelements, npointsPerElement, Dpoints, IdofsTest,&
      rdomainIntSubset, Dvalues, rcollection)

!<description>
    ! This subroutine is called during the calculation of errors. It
    ! has to compute the (analytical) values of a function in a couple
    ! of points on a couple of elements. These values are compared to
    ! those of a computed FE function and used to calculate an error.
    ! This routine can be used universally for arbitrary reference
    ! functions evaluated analytically using a function parser which
    ! is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to
    ! compute simultaneously for all these points.
!</description>

!<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(in) :: cderivative

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   DquickAccess(1): simulation time
    !   SquickAccess(1): section name in the collection
    !   SquickAccess(2): string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(1) :: dtime
    integer :: iel, icomp

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access double
    ! value holds the simulation time
    dtime(1) = rcollection%DquickAccess(1)

    ! Get number of the analytic reference function
    icomp = rcollection%IquickAccess(cderivative)

    ! Evaluate all values using the function parser
    !$omp parallel do
    do iel = 1, nelements
      call fparser_evalFunction(p_rfparser, icomp, 2,&
          Dpoints(:,:,iel), Dvalues(:,iel), dtime)
    end do
    !$omp end parallel do
    
  end subroutine transp_refFuncAnalytic

  !*****************************************************************************

!<subroutine>

  subroutine transp_refFuncAnalyticBdr2D(cderivative, rdiscretisation,&
      DpointsRef, Dpoints, ibct, DpointPar, Ielements, Dvalues, rcollection)

!<description>
    ! This subroutine is called during the calculation of errors with boundary
    ! integrals. It has to compute the (analytical) values of a
    ! function in a couple of points Dpoints on a couple of elements
    ! Ielements. These values are compared to those of a computed FE
    ! function to calculate an error.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real and reference coordinates.
    ! It has to to simultaneously compute the desired values for all these points.
!</description>

!<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(in) :: cderivative

    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates correspond to the reference
    ! element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: DpointsRef

    ! This is an array of all points on all the elements where coefficients
    ! are needed. It specifies the coordinates of the points where
    ! information is needed. These coordinates are world coordinates,
    ! i.e. on the real element.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! This is the number of the boundary component that contains the
    ! points in Dpoint. All points are on the same boundary component.
    integer, intent(in) :: ibct

    ! For every point under consideration, this specifies the parameter
    ! value of the point on the boundary component. The parameter value
    ! is calculated in LENGTH PARAMETRISATION!
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(in) :: DpointPar

    ! This is a list of elements (corresponding to Dpoints) where information
    ! is needed. To an element iel=Ielements(i), the array Dpoints(:,:,i)
    ! specifies the points where information is needed.
    ! DIMENSION(nelements)
    integer, dimension(:), intent(in) :: Ielements

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(1) :: dtime
    integer :: iel, icomp

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access double
    ! value holds the simulation time
    dtime(1) = rcollection%DquickAccess(1)

    ! Get number of the analytic reference function
    icomp = rcollection%IquickAccess(cderivative)

    ! Evaluate all values using the function parser
    !$omp parallel do
    do iel = 1, size(Dpoints,3)
      call fparser_evalFunction(p_rfparser, icomp, 2,&
          Dpoints(:,:,iel), Dvalues(:,iel), dtime)
    end do
    !$omp end parallel do
    
  end subroutine transp_refFuncAnalyticBdr2D

  !*****************************************************************************

!<subroutine>

  subroutine transp_weightFuncAnalytic(rdiscretisation, nelements,&
      npointsPerElement, Dpoints, IdofsTest, rdomainIntSubset,&
      Dvalues, rcollection)

!<description>
    ! This subroutine is called during the calculation of errors. It
    ! has to compute the values of a weighting function in a couple of
    ! points on a couple of elements. These values are multiplied by
    ! the calculated error.  This routine can be used universally for
    ! arbitrary weighting functions evaluated analytically using a
    ! function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to
    ! compute simultaneously for all these points.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   DquickAccess(1): simulation time
    !   SquickAccess(1): section name in the collection
    !   SquickAccess(2): string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(1) :: dtime
    integer :: iel, icomp

    ! This subroutine assumes that the first and second quick access
    ! string values hold the section name and the name of the function
    ! parser in the collection, respectively.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(2)),&
        ssectionName=trim(rcollection%SquickAccess(1)))

    ! This subroutine assumes that the first quick access double
    ! value holds the simulation time
    dtime(1)  = rcollection%DquickAccess(1)

    ! Get number of the analytic reference function
    icomp = rcollection%IquickAccess(2)

    ! Evaluate all values by the function parser
    !$omp parallel do
    do iel = 1, nelements
      call fparser_evalFunction(p_rfparser, icomp, 2,&
          Dpoints(:,:,iel), Dvalues(:,iel), dtime)
    end do
    !$omp end parallel do

  end subroutine transp_weightFuncAnalytic

  !*****************************************************************************

!<subroutine>

  subroutine transp_weightFuncGHill(rdiscretisation, nelements,&
      npointsPerElement, Dpoints, IdofsTest, rdomainIntSubset,&
      Dvalues, rcollection)

!<description>
    ! This subroutine is called during the calculation of errors. It
    ! has to compute the values of a weighting function in a couple of
    ! points on a couple of elements. These values are multiplied by
    ! the calculated error.  This routine can be used universally for
    ! arbitrary weighting functions evaluated analytically using a
    ! function parser which is passed using the collection.
    !
    ! The routine accepts a set of elements and a set of points on
    ! these elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to
    ! compute simultaneously for all these points.
!</description>

!<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation

    ! Number of elements, where the coefficients must be computed.
    integer, intent(in) :: nelements

    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement

    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! An array accepting the DOF`s on all elements in the test space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    ! This subroutine assumes the following data:
    !   DquickAccess(1): simulation time
    !   SquickAccess(1): section name in the collection
    !   SquickAccess(2): string identifying the function parser
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>
!</subroutine>

    ! local variables
    real(DP) :: dxhat,dyhat
    integer :: iel,ipoint

    select case(rcollection%IquickAccess(1))
    case(1)
      ! weight for x_hat
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          Dvalues(ipoint,iel) = -Dpoints(1,ipoint,iel)
        end do
      end do

    case(2) ! weight for y_hat
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          Dvalues(ipoint,iel) = -Dpoints(2,ipoint,iel)
        end do
      end do

    case(3) ! sigma
      dxhat = rcollection%DquickAccess(1)
      dyhat = rcollection%DquickAccess(2)
      
      do iel = 1, nelements
        do ipoint = 1, npointsPerElement
          Dvalues(ipoint,iel) = -(Dpoints(1,ipoint,iel)-dxhat)**2 -&
                                 (Dpoints(2,ipoint,iel)-dyhat)**2
        end do
      end do
      
    case default
      Dvalues = 0
    end select
    
  end subroutine transp_weightFuncGHill

  ! *****************************************************************************

!<subroutine>

  subroutine transp_parseBoundaryCondition(cbdrCondType, ndimension,&
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

    ! Homogeneous Neumann boundary conditions
    case ('HOMNEUMANN_STRONG')
      ibdrCondType = BDRC_HOMNEUMANN
      ! No strong boundary conditions are prescribed
    case ('HOMNEUMANN_WEAK')
      ibdrCondType = BDRC_HOMNEUMANN + BDRC_WEAK
    case ('HOMNEUMANN_WEAK_LUMPED')
      ibdrCondType = BDRC_HOMNEUMANN + BDRC_WEAK + BDRC_LUMPED
      
    ! Inhomogeneous Neumann boundary conditions
    case ('INHOMNEUMANN_STRONG')
      ibdrCondType = BDRC_INHOMNEUMANN
      ! No strong boundary conditions are prescribed
    case ('INHOMNEUMANN_WEAK')
      ibdrCondType = BDRC_INHOMNEUMANN + BDRC_WEAK
    case ('INHOMNEUMANN_WEAK_LUMPED')
      ibdrCondType = BDRC_INHOMNEUMANN + BDRC_WEAK + BDRC_LUMPED

    ! Dirichlet boundary conditions
    case ('DIRICHLET_STRONG')
      ibdrCondType = BDRC_DIRICHLET + BDRC_STRONG
    case ('DIRICHLET_WEAK')
      ibdrCondType = BDRC_DIRICHLET + BDRC_WEAK
    case ('DIRICHLET_WEAK_LUMPED')
      ibdrCondType = BDRC_DIRICHLET + BDRC_WEAK + BDRC_LUMPED
      

    ! Robin boundary conditions
    case ('ROBIN_STRONG')
      ibdrCondType = BDRC_ROBIN + BDRC_STRONG
    case ('ROBIN_WEAK')
      ibdrCondType = BDRC_ROBIN + BDRC_WEAK
    case ('ROBIN_LUMPED_WEAK')
      ibdrCondType = BDRC_ROBIN + BDRC_WEAK + BDRC_LUMPED

    ! Flux boundary conditions
    case ('FLUX_STRONG')
      ibdrCondType = BDRC_FLUX + BDRC_STRONG
    case ('FLUX_WEAK')
      ibdrCondType = BDRC_FLUX + BDRC_WEAK
    case ('FLUX_LUMPED_LUMPED')
      ibdrCondType = BDRC_FLUX + BDRC_WEAK + BDRC_LUMPED

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
          OU_CLASS_ERROR,OU_MODE_STD,'transp_parseBoundaryCondition')
      call sys_halt()
    end select

    
    ! Determine number of mathematical expressions
    select case (iand(ibdrCondType, BDRC_TYPEMASK))

    case (BDRC_HOMNEUMANN)
      nexpressions = 0
      
    case (BDRC_INHOMNEUMANN, BDRC_FLUX)
      nexpressions = 1 ! prescribed boundary value $g$

    case (BDRC_ROBIN)
      nexpressions = 2 ! prescribed boundary value $g$
                       ! multiplication factor $\alpha$

    case (BDRC_DIRICHLET)
      if (iand(ibdrCondType, BDRC_STRONG) .eq. BDRC_STRONG) then
        nexpressions = 1 ! prescribed boundary value $g$
      else
        nexpressions = 3 ! prescribed boundary value $g$
                         ! penalty parameter $\epsilon$
                         ! switching parameter $\gamma$
      end if

    case (BDRC_PERIODIC, BDRC_ANTIPERIODIC)
      nexpressions = 2 ! penalty parameter $\epsilon$
                       ! switching parameter $\gamma$

    case default
      nexpressions = 0
    end select

  end subroutine transp_parseBoundaryCondition

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcGeometricSourceterm(rparlist, ssectionName,&
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
    real(DP), dimension(:), pointer :: p_DdataVelocity
    real(DP), dimension(:), pointer :: p_DdataSource
    real(DP), dimension(:), pointer :: p_DdataMassMatrix
    integer, dimension(:), pointer :: p_Kld,p_Kcol
    character(LEN=SYS_STRLEN) :: smode
    real(DP) :: deffectiveRadius
    integer :: icoordsystem,massmatrix,velocityfield
    integer :: ivelocitytype,igeometricsourcetype,dofCoords

    ! Check of source and solution vector are compatible
    call lsysbl_isVectorCompatible(rsolution, rsource)


    ! Get parameters from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'mode', smode)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype, VELOCITY_ZERO)
    
    ! Are we in primal or dual mode?
    if (trim(smode) .eq. 'primal') then
    
      !-------------------------------------------------------------------------
      ! We are in primal mode
      !-------------------------------------------------------------------------
      
      ! @FAQ2: Which type of velocity are we?
      select case(abs(ivelocitytype))
        
      case (VELOCITY_ZERO, VELOCITY_BURGERS1D, VELOCITY_BURGERS2D,&
            VELOCITY_BUCKLEV1D, VELOCITY_BURGERS_SPACETIME,&
            VELOCITY_BUCKLEV_SPACETIME)
        ! zero velocity, Burgers or Buckley-Leverett equation
        ! do nothing (i.e. clear the source vector if required)
        if (bclear) call lsysbl_clearVector(rsource)
        
      case (VELOCITY_EXTERNAL,VELOCITY_CONSTANT,VELOCITY_TIMEDEP)
        
        !-----------------------------------------------------------------------
        ! problem structure has an explicit velocity field
        !-----------------------------------------------------------------------
        
        ! Get further parameters from parameter list
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'icoordsystem', icoordsystem, COORDS_CARTESIAN)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'igeometricsourcetype', igeometricsourcetype, MASS_LUMPED)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'velocityfield', velocityfield, 0)
        call parlst_getvalue_double(rparlist,&
            ssectionName, 'deffectiveRadius', deffectiveRadius, 1e-4_DP)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'dofCoords', dofCoords, 0)
        
        if ((velocityfield .gt. 0) .and. (dofCoords .gt. 0)) then

          ! Get pointers
          call lsysbl_getbase_double(rsolution, p_DdataSolution)
          call lsysbl_getbase_double(rsource, p_DdataSource)
          call lsysbl_getbase_double(&
              rproblemLevel%RvectorBlock(velocityfield), p_DdataVelocity)
          
          ! Get coordinates of the global DOF`s
          call lsysbl_getbase_double(&
              rproblemLevel%RvectorBlock(dofCoords), p_DdofCoords)
          
          ! What type of coordinate system are we?
          select case(icoordsystem)
          case (COORDS_CARTESIAN)
            ! No geometric source term required, clear vector (if required)
            if (bclear) call lsysbl_clearVector(rsource)
            
            
          case (COORDS_AXIALSYMMETRY,COORDS_CYLINDRICALSYMMETRY)
            ! Axi-symmetric (dalpha=1) flow (2D approximation to 3D flow) or
            ! cylindrically symmetric (dalpha=1) flow (1D approximation to 2D flow)
            
            select case(igeometricsourcetype)
            case (MASS_LUMPED)
              call parlst_getvalue_int(rparlist,&
                  ssectionName, 'lumpedmassmatrix', massmatrix)
              call lsyssc_getbase_double(&
                  rproblemLevel%RmatrixScalar(massmatrix), p_DdataMassMatrix)
              call doSourceVelocityLumped(dscale, deffectiveRadius,&
                  rproblemLevel%rtriangulation%ndim, rsource%NEQ,&
                  bclear, p_DdofCoords, p_DdataMassMatrix,&
                  p_DdataVelocity, p_DdataSolution, p_DdataSource)
              
            case (MASS_CONSISTENT)
              call parlst_getvalue_int(rparlist,&
                  ssectionName, 'consistentmassmatrix', massmatrix)
              call lsyssc_getbase_double(&
                  rproblemLevel%RmatrixScalar(massmatrix), p_DdataMassMatrix)
              call lsyssc_getbase_Kld(rproblemLevel%RmatrixScalar(massmatrix), p_Kld)
              call lsyssc_getbase_Kcol(rproblemLevel%RmatrixScalar(massmatrix), p_Kcol)
              call doSourceVelocityConsistent(dscale, deffectiveRadius,&
                  rproblemLevel%rtriangulation%ndim, rsource%NEQ,&
                  bclear, p_DdofCoords, p_DdataMassMatrix, p_Kld, p_Kcol,&
                  p_DdataVelocity, p_DdataSolution, p_DdataSource)
              
            case default
              call output_line('Unsupported geometric source type!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'transp_calcGeometricSourceterm')
              call sys_halt()
            end select
            
            
          case (COORDS_SPHERICALSYMMETRY)
            ! Spherically symmetric (dalpha=2) flow (1D approximation to 3D flow)
            
            select case(igeometricsourcetype)
            case (MASS_LUMPED)
              call parlst_getvalue_int(rparlist,&
                  ssectionName, 'lumpedmassmatrix', massmatrix)
              call lsyssc_getbase_double(&
                  rproblemLevel%RmatrixScalar(massmatrix), p_DdataMassMatrix)
              call doSourceVelocityLumped(2*dscale, deffectiveRadius,&
                  rproblemLevel%rtriangulation%ndim, rsource%NEQ,&
                  bclear, p_DdofCoords, p_DdataMassMatrix,&
                  p_DdataVelocity, p_DdataSolution, p_DdataSource)
              
            case (MASS_CONSISTENT)
              call parlst_getvalue_int(rparlist,&
                  ssectionName, 'lumpedmassmatrix', massmatrix)
              call lsyssc_getbase_double(&
                  rproblemLevel%RmatrixScalar(massmatrix), p_DdataMassMatrix)
              call lsyssc_getbase_Kld(rproblemLevel%RmatrixScalar(massmatrix), p_Kld)
              call lsyssc_getbase_Kcol(rproblemLevel%RmatrixScalar(massmatrix), p_Kcol)
              call doSourceVelocityConsistent(2*dscale, deffectiveRadius,&
                  rproblemLevel%rtriangulation%ndim, rsource%NEQ,&
                  bclear, p_DdofCoords, p_DdataMassMatrix, p_Kld,&
                  p_Kcol, p_DdataVelocity, p_DdataSolution, p_DdataSource)
              
            case default
              call output_line('Unsupported geometric source type!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'transp_calcGeometricSourceterm')
              call sys_halt()
            end select
            
          case default
            call output_line('Invalid coordinate system!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcGeometricSourceterm')
            call sys_halt()
          end select
          
        else
          call output_line('Coordinates of DOFs not available!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcGeometricSourceterm')
          call sys_halt()
        end if
        
      case default
        call output_line('Unsupported velocity type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_calcGeometricSourceterm')
        call sys_halt()
      end select
      
    elseif (trim(smode) .eq. 'dual') then

      !-------------------------------------------------------------------------
      ! We are in dual mode
      !-------------------------------------------------------------------------

      ! @FAQ2: Which type of velocity are we?
      select case(abs(ivelocitytype))

      case (VELOCITY_ZERO, VELOCITY_BURGERS1D, VELOCITY_BURGERS2D,&
            VELOCITY_BUCKLEV1D, VELOCITY_BURGERS_SPACETIME,&
            VELOCITY_BUCKLEV_SPACETIME)
        ! zero velocity, 1D-Burgers or 1D-Buckley-Leverett equation
        ! do nothing (i.e. clear the source vector if required)
        if (bclear) call lsysbl_clearVector(rsource)

      case (VELOCITY_EXTERNAL,VELOCITY_CONSTANT,VELOCITY_TIMEDEP)

        !-----------------------------------------------------------------------
        ! problem structure has an explicit velocity field
        !-----------------------------------------------------------------------
        
        ! Get further parameters from parameter list
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'icoordsystem', icoordsystem, COORDS_CARTESIAN)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'igeometricsourcetype', igeometricsourcetype, MASS_LUMPED)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'velocityfield', velocityfield, VELOCITY_ZERO)
        call parlst_getvalue_double(rparlist,&
            ssectionName, 'deffectiveRadius', deffectiveRadius, 1e-4_DP)
        call parlst_getvalue_int(rparlist,&
            ssectionName, 'dofCoords', dofCoords, 0)
        
        if ((velocityfield .gt. 0) .and. (dofCoords .gt. 0)) then
          
          ! Get pointers
          call lsysbl_getbase_double(rsolution, p_DdataSolution)
          call lsysbl_getbase_double(rsource, p_DdataSource)
          call lsysbl_getbase_double(&
              rproblemLevel%RvectorBlock(velocityfield), p_DdataVelocity)
          
          ! Get coordinates of the global DOF`s
          call lsysbl_getbase_double(&
              rproblemLevel%RvectorBlock(dofCoords), p_DdofCoords)
          
          ! What type of coordinate system are we?
          select case(icoordsystem)
          case (COORDS_CARTESIAN)
            ! No geometric source term required, clear vector (if required)
            if (bclear) call lsysbl_clearVector(rsource)
            
          case default
            call output_line('Invalid coordinate system!',&
                OU_CLASS_ERROR,OU_MODE_STD,'transp_calcGeometricSourceterm')
            call sys_halt()
          end select
          
        else
          call output_line('Coordinates of DOFs not available!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcGeometricSourceterm')
          call sys_halt()
        end if
        
      case default
        call output_line('Unsupported velocity type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_calcGeometricSourceterm')
        call sys_halt()
      end select
      
    else
      call output_line('Invalid mode!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcGeometricSourceterm')
      call sys_halt()
    end if
    
  contains

    ! Here, the real working routines follow

    !**************************************************************
    ! Calculate the geometric source term for axi-symmetric flow
    ! (dalpha=1) in 2D and cylindrically (dalpha=1) or spherically
    ! symmetric (dalpha=2) flow in 1D. This routine assembles the
    ! geometric source term for a given velocity field and solution
    ! vector using the lumped mass matrix.
    
    subroutine doSourceVelocityLumped(deffectiveScale,&
        deffectiveRadius, ndim, neq, bclear, Dcoords, DdataMassMatrix,&
        DdataVelocity, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Spatial dimension
      integer, intent(in) :: ndim

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq

      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:), intent(in) :: Dcoords

      ! Lumped mass matrix
      real(DP), dimension(:), intent(in) :: DdataMassMatrix

      ! Velocity vector
      real(DP), dimension(:), intent(in) :: DdataVelocity

      ! Solution vector
      real(DP), dimension(:), intent(in) :: DdataSolution

      
      ! Source vector
      real(DP), dimension(:), intent(inout) :: DdataSource


      ! local variables
      real(DP) :: daux, dradius
      integer :: ieq

      
      ! Do we have to clear the source vector?
      if (bclear) then

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius) if(neq .gt. TRANSP_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(ndim*(ieq-1)+1)
          dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Overwrite the geometric source term
          DdataSource(ieq) = daux * DdataMassMatrix(ieq) *&
                             DdataVelocity(ieq) * DdataSolution(ieq)
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,dradius) if(neq .gt. TRANSP_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Get the r-coordinate and compute the radius
          daux = Dcoords(ndim*(ieq-1)+1)
          dradius = max(abs(daux), deffectiveRadius)

          ! Compute unit vector into the origin, scale it be the user-
          ! defined scaling parameter dscale and devide it by the radius
          daux = -sign(1.0_DP, daux) * deffectiveScale / dradius

          ! Update the geometric source term
          DdataSource(ieq) = DdataSource(ieq) + daux * DdataMassMatrix(ieq) *&
                             DdataVelocity(ieq) * DdataSolution(ieq)
        end do
        !$omp end parallel do

      end if

    end subroutine doSourceVelocityLumped

    !**************************************************************
    ! Calculate the geometric source term for axi-symmetric flow
    ! (dalpha=1) in 2D and cylindrically (dalpha=1) or spherically
    ! symmetric (dalpha=2) flow in 1D. This routine assembles the
    ! geometric source term for a given velocity field and solution
    ! vector using the consistent mass matrix.
    
    subroutine doSourceVelocityConsistent(deffectiveScale,&
        deffectiveRadius, ndim, neq, bclear, Dcoords, DdataMassMatrix,&
        Kld, Kcol, DdataVelocity, DdataSolution, DdataSource)

      ! Effective scaling parameter (dalpha * dscale)
      real(DP), intent(in) :: deffectiveScale

      ! Effectiive radius
      real(DP), intent(in) :: deffectiveRadius

      ! Spatial dimension
      integer, intent(in) :: ndim

      ! Number of equation (nodal degrees of freedom)
      integer, intent(in) :: neq

      ! Clear source vector?
      logical, intent(in) :: bclear

      ! Coordinates of the nodal degrees of freedom
      real(DP), dimension(:), intent(in) :: Dcoords

      ! Consistent mass matrix
      real(DP), dimension(:), intent(in) :: DdataMassMatrix

      ! Sparsity pattern of the mass matrix
      integer, dimension(:), intent(in) :: Kld,Kcol

      ! Velocity vector
      real(DP), dimension(:), intent(in) :: DdataVelocity

      ! Solution vector
      real(DP), dimension(:), intent(in) :: DdataSolution

      
      ! Source vector
      real(DP), dimension(:), intent(inout) :: DdataSource


      ! local variables
      real(DP) :: ddata, daux, dradius
      integer :: ieq,ia,jeq

      
      ! Do we have to clear the source vector?
      if (bclear) then

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,ddata,dradius,ia,jeq) if(neq .gt. TRANSP_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Clear temporal data
          ddata = 0.0_DP

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
            
            ! Update the geometric source term
            ddata = ddata + daux * DdataMassMatrix(ia) *&
                             DdataVelocity(jeq) * DdataSolution(jeq)
          end do
          
          ! Overwrite the geometric source term
          DdataSource(ieq) = ddata
        end do
        !$omp end parallel do

      else ! bclear

        ! Loop over all degrees of freedom
        !$omp parallel do default(shared)&
        !$omp private(daux,ddata,dradius,ia,jeq) if(neq .gt. TRANSP_GEOMSOURCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Clear temporal data
          ddata = 0.0_DP

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
            
            ! Update the geometric source term
            ddata = ddata + daux * DdataMassMatrix(ia) *&
                             DdataVelocity(jeq) * DdataSolution(jeq)
          end do
          
          ! Update the geometric source term
          DdataSource(ieq) = DdataSource(ieq) + ddata
        end do
        !$omp end parallel do

      end if
      
    end subroutine doSourceVelocityConsistent

  end subroutine transp_calcGeometricSourceterm

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcTransportOperator(rproblemLevel, rboundaryCondition,&
      rsolution, dtime, dscale, rmatrix, ssectionName, rcollection,&
      fcb_calcMatrixDiagPrimal_sim, fcb_calcMatrixPrimal_sim,&
      fcb_calcMatrixDiagDual_sim,  fcb_calcMatrixDual_sim,&
      fcb_coeffMatBdrPrimal1d_sim, fcb_coeffMatBdrDual1d_sim,&
      fcb_coeffMatBdrPrimal2d_sim, fcb_coeffMatBdrDual2d_sim,&
      fcb_coeffMatBdrPrimal3d_sim, fcb_coeffMatBdrDual3d_sim,&
      bforceUpdate)

!<description>
    ! This subroutine calculates the discrete transport operator for
    ! the primal and dual problem, respectively. The resulting
    ! operator includes convective and diffusive terms as well as
    ! contributions from boundary terms due to weakly imposed boundary
    ! condition.

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

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCalcMatrix.inc'
    optional :: fcb_calcMatrixDiagPrimal_sim
    optional :: fcb_calcMatrixPrimal_sim
    optional :: fcb_calcMatrixDiagDual_sim
    optional :: fcb_calcMatrixDual_sim

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCoeffMatBdr.inc'
    optional :: fcb_coeffMatBdrPrimal1d_sim
    optional :: fcb_coeffMatBdrPrimal2d_sim
    optional :: fcb_coeffMatBdrPrimal3d_sim
    optional :: fcb_coeffMatBdrDual1d_sim
    optional :: fcb_coeffMatBdrDual2d_sim
    optional :: fcb_coeffMatBdrDual3d_sim

    ! OPTIONAL: flag to overwrite update notifier
    ! If true, then the transport operator is updated even if the
    ! update notifier does not indicate the need of an update
    logical, intent(in), optional :: bforceUpdate
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! destination matrix for the transport operator
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_collection) :: rcollectionTmp
    character(LEN=SYS_STRLEN) :: smode
    integer :: ivelocitytype,idiffusiontype,coeffMatrix_S
    integer :: convectionAFC,diffusionAFC,velocityfield
    integer :: convectionGFEM,diffusionGFEM
    logical :: breturn

    !###########################################################################
    ! REMARK: If you want to add a new type of velocity/diffusion,
    ! then search for the tag @FAQ2: in this subroutine and create a
    ! new CASE which performs the corresponding task for the new type
    ! of velocity/ diffusion.
    !###########################################################################

    breturn = .true.
    if (present(bforceUpdate)) breturn = .not.(bforceUpdate)

    ! Check if the discrete transport operator has to be updated and
    ! return otherwise.
    if ((iand(rproblemLevel%iproblemSpec, TRANSP_TROPER_INIT)   .eq. 0) .and.&
        (iand(rproblemLevel%iproblemSpec, TRANSP_TROPER_UPDATE) .eq. 0) .and.&
        breturn) return
    
    ! Set pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)

    ! Get positions of coefficient matrices from parameter list
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'coeffmatrix_S', coeffMatrix_S)

    ! Remove update notifier from the transport operator. Depending on
    ! the velocity, diffusion type it will be re-activited below.
    rproblemLevel%iproblemSpec = iand(rproblemLevel%iproblemSpec,&
        not(TRANSP_TROPER_INIT+TRANSP_TROPER_UPDATE))

    !---------------------------------------------------------------------------
    ! Assemble diffusion operator (for the right-hand side):
    !
    ! After integration by parts it is given by
    !
    ! $$ Su = -\int_\Omega \nabla w \cdot (D \nabla u) {\rm d}{\bf x} $$
    !
    ! The diffusion operator is symmetric so that it is the same for
    ! the primal and the dual problem. If no diffusion is present,
    ! i.e. $D \equiv 0$, then the transport operator is initialised by
    ! zeros. If there is anisotropic diffusion, i.e. $D=D({\bf x},t)$,
    ! then we may also initialise some stabilisation structure.
    !
    ! The bilinear form for the diffusion operator consists of the
    ! volume integral (see above) only.
    !
    ! Non-homogeneous Neumann boundary conditions
    !
    ! $$ q_N = {\bf n} \cdot \nabla u = h \qquad \mbox{on} \quad \Gamma_N $$
    !
    ! are built into the right-hand side vector.
    !---------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'idiffusiontype', idiffusiontype)

    ! Primal and dual mode are equivalent
    ! @FAQ2: Which type of diffusion are we?
    select case(idiffusiontype)
    case (DIFFUSION_ZERO)
      ! zero diffusion, clear the transport matrix
      call lsyssc_clearMatrix(rmatrix)

      !-------------------------------------------------------------------------

    case (DIFFUSION_ISOTROPIC)
      ! Isotropic diffusion
      call lsyssc_duplicateMatrix(&
          rproblemLevel%RmatrixScalar(coeffMatrix_S),&
          rmatrix, LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPY)
      if (dscale .ne. 1.0_DP)&
          call lsyssc_scaleMatrix(rmatrix, dscale)

      !-------------------------------------------------------------------------

    case (DIFFUSION_ANISOTROPIC)
      ! Anisotropic diffusion
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'diffusionAFC', diffusionAFC)
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'diffusionGFEM', diffusionGFEM, diffusionAFC)

      if (diffusionAFC .gt. 0) then

        ! What kind of stabilisation should be applied?
        select case(rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType)

        case (AFCSTAB_DMP, AFCSTAB_SYMMETRIC,&
              AFCSTAB_LINLPT_SYMMETRIC,&
              AFCSTAB_NLINLPT_SYMMETRIC)
          ! Satisfy discrete maximum principle
          ! and assemble stabilisation structure
          call gfsc_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(diffusionGFEM)%RgroupFEMBlock(1),&
              dscale, .true., rmatrix, rafcstab=rproblemLevel%Rafcstab(diffusionAFC))

        case default
          ! Compute the standard Galerkin approximation
          call lsyssc_duplicateMatrix(&
          rproblemLevel%RmatrixScalar(coeffMatrix_S),&
          rmatrix, LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPY)
          if (dscale .ne. 1.0_DP)&
              call lsyssc_scaleMatrix(rmatrix, dscale)
        end select

      else ! diffusionAFC .lt. 0

        ! Compute the standard Galerkin approximation
        call lsyssc_duplicateMatrix(&
          rproblemLevel%RmatrixScalar(coeffMatrix_S),&
          rmatrix, LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPY)
        if (dscale .ne. 1.0_DP)&
            call lsyssc_scaleMatrix(rmatrix, dscale)

      end if   ! diffusionAFC

      !-------------------------------------------------------------------------

    case (DIFFUSION_VARIABLE)
      call output_line('Variable diffusion matrices are yet not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcTransportOperator')
      call sys_halt()

      ! Set update notifier for transport operator in problem level structure
      rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                       TRANSP_TROPER_UPDATE)

    case default
      call output_line('Invalid type of diffusion!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcTransportOperator')
      call sys_halt()
    end select

    
    !---------------------------------------------------------------------------
    ! Assemble convective operator (for the right-hand side):
    !
    ! (1) without integration by parts:
    !     $$ -\int_\Omega w \nabla \cdot {\bf f}(u) {\rm d}{\bf x} $$
    !
    !     then boundary conditions need to be imposed in strong sense
    !     by filtering the system matrix, the solution vector and/or
    !     the residual explicitly.
    !
    ! (2) with integration by parts:
    !     $$ \int_\Omega \nabla w \cdot {\bf f}(u) {\rm d}{\bf x} $$
    !
    !     with weakly imposed boundary conditions
    !
    !     $$ -\int_{\Gamma_-} w {\bf f}(u_0) \cdot {\bf n} {\rm d}{\bf s} $$
    !
    !     imposed at some part of the boundary. At the remaining part
    !     of the boundary nothing is prescribed and the corresponding
    !     boundary integral is built into the bilinear form
    !
    !     $$ -\int_{\Gamma_+} w {\bf f}(u) \cdot {\bf n} {\rm d}{\bf s} $$
    !
    ! The convective operator is skew-symmetric so that we have to
    ! distinguish between the primal and the dual problem.
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(p_rparlist,&
        ssectionName, 'mode', smode)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'convectionAFC', convectionAFC)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'convectionGFEM', convectionGFEM, convectionAFC)

    ! Attach user-defined collection structure to temporal collection
    ! structure (may be required by the callback function)
    rcollectionTmp%p_rnextCollection => rcollection

    ! Set simulation time
    rcollectionTmp%DquickAccess(1) = dtime

    ! Set vector for velocity field (if any)
    if (transp_hasVelocityVector(ivelocityType)) then
      call parlst_getvalue_int(p_rparlist,&
          ssectionName, 'velocityfield', velocityfield)
      rcollectionTmp%p_rvectorQuickAccess1 => rproblemLevel%RvectorBlock(velocityfield)
    end if


    ! Are we in primal or dual mode?
    if (trim(smode) .eq. 'primal') then

      !-------------------------------------------------------------------------
      ! We are in primal mode which means that we have to build two
      ! bilinear forms: one for the volume integral and one for the
      ! surface integral which is used to weakly impose Dirichlet
      ! boundary conditions. Note that lumping is performed for the
      ! boundary term to prevent the generation of oscillations.
      !-------------------------------------------------------------------------

      ! @FAQ2: Which type of velocity are we?
      select case(abs(ivelocitytype))

      case default
        ! The user-defined callback function for matrix coefficients
        ! is used if present; otherwise an error is thrown
        if (present(fcb_calcMatrixPrimal_sim) .and.&
            present(fcb_calcMatrixDiagPrimal_sim)) then

          select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
            
          case (AFCSTAB_GALERKIN)
            
            ! Apply standard Galerkin discretisation
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, fcb_calcMatrixPrimal_sim, dscale, .false., rmatrix,&
                fcb_calcMatrixDiagPrimal_sim, rcollection=rcollectionTmp)
            
          case default
            
            ! Apply low-order discretisation
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, fcb_calcMatrixPrimal_sim, dscale, .false., rmatrix,&
                fcb_calcMatrixDiagPrimal_sim, rcollection=rcollectionTmp,&
                rafcstab=rproblemLevel%Rafcstab(convectionAFC))
            
          end select

        else ! callback function not present

          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcTransportOperator')
          call sys_halt()

        end if

        !-----------------------------------------------------------------------

      case (VELOCITY_ZERO)
        ! zero velocity, do nothing

        !-----------------------------------------------------------------------

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP)
        ! linear velocity

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)

          ! Apply standard Galerkin discretisation
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatGalConvP1d_sim, dscale, .false., rmatrix,&
                transp_calcMatDiagConvP1d_sim, rcollection=rcollectionTmp)

          case (NDIM2D)
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatGalConvP2d_sim, dscale, .false., rmatrix,&
                transp_calcMatDiagConvP2d_sim, rcollection=rcollectionTmp)

          case (NDIM3D)
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatGalConvP3d_sim, dscale, .false., rmatrix,&
                transp_calcMatDiagConvP3d_sim, rcollection=rcollectionTmp)
          end select

        case default

          ! Apply low-order discretisation
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwConvP1d_sim, dscale, .false., rmatrix,&
                transp_calcMatDiagConvP1d_sim, rcollection=rcollectionTmp,&
                rafcstab=rproblemLevel%Rafcstab(convectionAFC))

          case (NDIM2D)
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwConvP2d_sim, dscale, .false., rmatrix,&
                transp_calcMatDiagConvP2d_sim, rcollection=rcollectionTmp,&
                rafcstab=rproblemLevel%Rafcstab(convectionAFC))

          case (NDIM3D)
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwConvP3d_sim, dscale, .false., rmatrix,&
                transp_calcMatDiagConvP3d_sim, rcollection=rcollectionTmp,&
                rafcstab=rproblemLevel%Rafcstab(convectionAFC))
          end select

        end select

        if (abs(ivelocitytype) .eq. VELOCITY_TIMEDEP) then
          ! Set update notifier for transport operator in problem level structure
          rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                           TRANSP_TROPER_UPDATE)
        end if

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers` equation in space-time

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)

          ! Apply standard Galerkin discretisation
          call gfsc_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatGalSTBurgP2d_sim, dscale, .false., rmatrix,&
              transp_calcMatDiagSTBurgP2d_sim, rcollection=rcollectionTmp)

        case default

          ! Apply low-order discretisation
          call gfsc_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwSTBurgP2d_sim, dscale, .false., rmatrix,&
              transp_calcMatDiagSTBurgP2d_sim, rcollection=rcollectionTmp,&
              rafcstab=rproblemLevel%Rafcstab(convectionAFC))

        end select

        ! Set update notifier for transport operator in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         TRANSP_TROPER_UPDATE)

        !-----------------------------------------------------------------------

      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)

          ! Apply standard Galerkin discretisation
          call gfsc_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatGalSTBLevP2d_sim, dscale, .false., rmatrix,&
              transp_calcMatDiagSTBLevP2d_sim, rcollection=rcollectionTmp)

        case default

          ! Apply low-order discretisation
          call gfsc_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwSTBLevP2d_sim, dscale, .false., rmatrix,&
              transp_calcMatDiagSTBLevP2d_sim, rcollection=rcollectionTmp,&
              rafcstab=rproblemLevel%Rafcstab(convectionAFC))

        end select

        ! Set update notifier for transport operator in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         TRANSP_TROPER_UPDATE)

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers` equation in 1D

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)

          ! Apply standard Galerkin discretisation
          call gfsc_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatGalBurgP1d_sim, dscale, .false., rmatrix,&
              transp_calcMatDiagBurgP1d_sim, rcollection=rcollectionTmp)

        case default

          ! Apply low-order discretisation
          call gfsc_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwBurgP1d_sim, dscale, .false., rmatrix,&
              transp_calcMatDiagBurgP1d_sim, rcollection=rcollectionTmp,&
              rafcstab=rproblemLevel%Rafcstab(convectionAFC))

        end select

        ! Set update notifier for transport operator in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         TRANSP_TROPER_UPDATE)

        !-----------------------------------------------------------------------

      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers` equation in 2D

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)

          ! Apply standard Galerkin discretisation
          call gfsc_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatGalBurgP2d_sim, dscale, .false., rmatrix,&
              transp_calcMatDiagBurgP2d_sim, rcollection=rcollectionTmp)

        case default

          ! Apply low-order discretisation
          call gfsc_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwBurgP2d_sim, dscale, .false., rmatrix,&
              transp_calcMatDiagBurgP2d_sim, rcollection=rcollectionTmp,&
              rafcstab=rproblemLevel%Rafcstab(convectionAFC))

        end select

        ! Set update notifier for transport operator in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         TRANSP_TROPER_UPDATE)

        !-----------------------------------------------------------------------

      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)

          ! Apply standard Galerkin discretisation
          call gfsc_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatGalBLevP1d_sim, dscale, .false., rmatrix,&
              transp_calcMatDiagBLevP1d_sim, rcollection=rcollectionTmp)

        case default

          ! Apply low-order discretisation
          call gfsc_buildOperatorEdge(&
              rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
              rsolution, transp_calcMatUpwBLevP1d_sim, dscale, .false., rmatrix,&
              transp_calcMatDiagBLevP1d_sim, rcollection=rcollectionTmp,&
              rafcstab=rproblemLevel%Rafcstab(convectionAFC))

        end select

        ! Set update notifier for transport operator in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         TRANSP_TROPER_UPDATE)
      end select

      !-------------------------------------------------------------------------
           
      ! Evaluate bilinear form for boundary integral (if any)
      call transp_calcBilfBdrCondition(rproblemLevel,&
          rboundaryCondition, rsolution, ssectionName, smode,&
          ivelocitytype, dTime, dscale, .false., rmatrix, rcollection,&
          fcb_coeffMatBdrPrimal1d_sim, fcb_coeffMatBdrDual1d_sim,&
          fcb_coeffMatBdrPrimal2d_sim, fcb_coeffMatBdrDual2d_sim,&
          fcb_coeffMatBdrPrimal3d_sim, fcb_coeffMatBdrDual3d_sim)

    elseif (trim(smode) .eq. 'dual') then

      !-------------------------------------------------------------------------
      ! We are in dual mode which means that we have to build two
      ! bilinear forms: one for the volume integral and one for the
      ! surface integral which is used to weakly impose Dirichlet
      ! boundary conditions. Note that lumping is performed for the
      ! boundary term to prevent the generation of oscillations.
      !-------------------------------------------------------------------------

      ! @FAQ2: Which type of velocity are we?
      select case(abs(ivelocitytype))
      case default
        ! The user-defined callback function for matrix coefficients
        ! is used if present; otherwise an error is thrown
        if (present(fcb_calcMatrixDual_sim) .and.&
            present(fcb_calcMatrixDiagDual_sim)) then

          select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)
            
          case (AFCSTAB_GALERKIN)
            
            ! Apply standard Galerkin discretisation
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, fcb_calcMatrixDual_sim, dscale, .false., rmatrix,&
                fcb_calcMatrixDiagDual_sim, rcollection=rcollectionTmp)
            
          case default
            
            ! Apply low-order discretisation
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, fcb_calcMatrixDual_sim, dscale, .false., rmatrix,&
                fcb_calcMatrixDiagDual_sim, rcollection=rcollectionTmp,&
                rafcstab=rproblemLevel%Rafcstab(convectionAFC))
            
          end select

        else ! callback function not present

          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcTransportOperator')
          call sys_halt()

        end if

        !-----------------------------------------------------------------------

      case (VELOCITY_ZERO)
        ! zero velocity, do nothing

        !-----------------------------------------------------------------------

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP)
        ! linear velocity

        select case(rproblemLevel%Rafcstab(convectionAFC)%cafcstabType)

        case (AFCSTAB_GALERKIN)

          ! Apply standard Galerkin discretisation
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatGalConvD1d_sim, dscale, .false., rmatrix,&
                transp_calcMatDiagConvD1d_sim, rcollection=rcollectionTmp)

          case (NDIM2D)
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatGalConvD2d_sim, dscale, .false., rmatrix,&
                transp_calcMatDiagConvD2d_sim, rcollection=rcollectionTmp)

          case (NDIM3D)
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatGalConvD3d_sim, dscale, .false., rmatrix,&
                transp_calcMatDiagConvD3d_sim, rcollection=rcollectionTmp)
          end select
          
        case default

          ! Apply low-order discretisation
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwConvD1d_sim, dscale, .false., rmatrix,&
                transp_calcMatDiagConvD1d_sim, rcollection=rcollectionTmp,&
                rafcstab=rproblemLevel%Rafcstab(convectionAFC))

          case (NDIM2D)
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwConvD2d_sim, dscale, .false., rmatrix,&
                transp_calcMatDiagConvD2d_sim, rcollection=rcollectionTmp,&
                rafcstab=rproblemLevel%Rafcstab(convectionAFC))

          case (NDIM3D)
            call gfsc_buildOperatorEdge(&
                rproblemLevel%RgroupFEMBlock(convectionGFEM)%RgroupFEMBlock(1),&
                rsolution, transp_calcMatUpwConvD3d_sim, dscale, .false., rmatrix,&
                transp_calcMatDiagConvD3d_sim, rcollection=rcollectionTmp,&
                rafcstab=rproblemLevel%Rafcstab(convectionAFC))
          end select

        end select

        if (abs(ivelocitytype) .eq. VELOCITY_TIMEDEP) then
          ! Set update notifier for transport operator in problem level structure
          rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                       TRANSP_TROPER_UPDATE)
        end if

        !-----------------------------------------------------------------------

        ! @TODO: The dual mode has only been implemented for linear
        ! convection. If you need to compute the dual problem for
        ! some other velocity type, then you have to add the
        ! implementation of the dual transport operator below!

      end select

      ! Evaluate bilinear form for boundary integral (if any)
      call transp_calcBilfBdrCondition(rproblemLevel,&
          rboundaryCondition, rsolution, ssectionName, smode,&
          ivelocitytype, dTime, dscale, .false., rmatrix, rcollection,&
          fcb_coeffMatBdrPrimal1d_sim, fcb_coeffMatBdrDual1d_sim,&
          fcb_coeffMatBdrPrimal2d_sim, fcb_coeffMatBdrDual2d_sim,&
          fcb_coeffMatBdrPrimal3d_sim, fcb_coeffMatBdrDual3d_sim)

    else
      call output_line('Invalid mode!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcTransportOperator')
      call sys_halt()
    end if

  end subroutine transp_calcTransportOperator

!*****************************************************************************

!<subroutine>

  subroutine transp_calcTimeDerivative(rproblemLevel, rtimestep,&
      rsolver, rsolution, ssectionName, rcollection, rvector, rsource,&
      fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
      fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
      fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim,&
      rmatrix, rvector1, rvector2)

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

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCoeffVecBdr.inc'
    optional :: fcb_coeffVecBdrPrimal1d_sim
    optional :: fcb_coeffVecBdrPrimal2d_sim
    optional :: fcb_coeffVecBdrPrimal3d_sim
    optional :: fcb_coeffVecBdrDual1d_sim
    optional :: fcb_coeffVecBdrDual2d_sim
    optional :: fcb_coeffVecBdrDual3d_sim
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! collection structure
    type(t_collection), intent(inout) :: rcollection

    ! destination vector
    type(t_vectorBlock), intent(inout) :: rvector

    ! OPTIONAL: matrix to store the discrete transport operator used
    ! to compute the approximation to the time derivative (if not
    ! present, then temporal memory is allocated)
    type(t_matrixScalar), intent(inout), target, optional :: rmatrix

    ! OPTIONAL: auxiliary vectors used to compute the approximation to
    ! the time derivative (if not present, then temporal memory is allocated)
    type(t_vectorBlock), intent(inout), target, optional :: rvector1
    type(t_vectorBlock), intent(inout), target, optional :: rvector2
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_matrixScalar), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rvector1, p_rvector2
    character(LEN=SYS_STRLEN) :: smode
    real(DP) :: dnorm0, dnorm
    real(DP) :: depsAbsApproxTimeDerivative,depsRelApproxTimeDerivative
    integer :: convectionAFC,diffusionAFC,ivelocitytype
    integer :: iapproxtimederivativetype
    integer :: lumpedMassMatrix,consistentMassMatrix
    integer :: cafcstabTypeConvection
    integer :: cafcstabTypeDiffusion
    integer :: ite,nmaxIterationsApproxTimeDerivative
    integer(I32) :: istabilisationSpecConvection
    integer(I32) :: istabilisationSpecDiffusion
    logical :: bcompatible, bforceUpdate

    ! Get parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)

    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'convectionAFC', convectionAFC, 0)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'diffusionAFC', diffusionAFC, 0)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'iapproxtimederivativetype', iapproxtimederivativetype)
    call parlst_getvalue_int(p_rparlist,&
        ssectionName, 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_string(p_rparlist,&
        ssectionName, 'mode', smode)

    ! Check if rvector is compatible to the solution vector;
    ! otherwise create new vector as a duplicate of the solution vector
    call lsysbl_isVectorCompatible(rvector, rsolution, bcompatible)
    if (.not.bcompatible)&
        call lsysbl_duplicateVector(rsolution, rvector,&
        LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
    
    ! Set up matrix for discrete transport operator
    if (present(rmatrix)) then
      p_rmatrix => rmatrix
    else
      allocate(p_rmatrix)
    end if

    ! Check if matrix is compatible to consistent mass matrix;
    ! otherwise create matrix as a duplicate of the consistent mass matrix
    call lsyssc_isMatrixCompatible(p_rmatrix,&
        rproblemLevel%RmatrixScalar(consistentMassMatrix), bcompatible)
    if (.not.bcompatible) then
      call lsyssc_duplicateMatrix(&
          rproblemLevel%RmatrixScalar(consistentMassMatrix),&
          p_rmatrix, LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      bforceUpdate = .true.
    else
      bforceUpdate = .false.
    end if
    
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
      if (convectionAFC .gt. 0) then
        cafcstabTypeConvection&
            = rproblemLevel%Rafcstab(convectionAFC)%cafcstabType
        istabilisationSpecConvection&
            = rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec
        rproblemLevel%Rafcstab(convectionAFC)%cafcstabType&
            = AFCSTAB_GALERKIN
      end if
      
      if (diffusionAFC .gt. 0) then
        cafcstabTypeDiffusion&
            = rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType
        istabilisationSpecDiffusion&
            = rproblemLevel%Rafcstab(diffusionAFC)%istabilisationSpec
        rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType&
            = AFCSTAB_GALERKIN
      end if
      
      ! Create discrete transport operator of high-order (Galerkin method)
      call transp_calcTransportOperator(rproblemLevel,&
          rsolver%rboundaryCondition, rsolution, rtimestep%dTime,&
          1.0_DP, p_rmatrix, ssectionName, rcollection,&
          bforceUpdate=bforceUpdate)
      
      ! Reset stabilisation structures to their original configuration
      if (convectionAFC .gt. 0) then
        rproblemLevel%Rafcstab(convectionAFC)%cafcstabType&
            = cafcstabTypeConvection
        rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec&
            = istabilisationSpecConvection
      end if
      
      if (diffusionAFC .gt. 0) then
        rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType&
            = cafcstabTypeDiffusion
        rproblemLevel%Rafcstab(diffusionAFC)%istabilisationSpec&
            = istabilisationSpecDiffusion
      end if
      
      ! Compute $K(u^L)*u^L$ and store the result in rvector1
      call lsyssc_matVec(p_rmatrix,&
          rsolution%rvectorBlock(1),&
          p_rvector1%RvectorBlock(1), 1.0_DP, 0.0_DP)
      
      ! Evaluate linear form for the boundary integral (if any)
      call transp_calcLinfBdrCondition(&
          rproblemLevel, rsolver%rboundaryCondition,&
          rsolution, ssectionName, smode, ivelocitytype,&
          rtimestep%dTime, 1.0_DP, .false., p_rvector1, rcollection,&
          fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
          fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
          fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim)
      
      ! Build the geometric source term (if any)
      call transp_calcGeometricSourceterm(p_rparlist, ssectionName,&
          rproblemLevel, rsolution, 1.0_DP, .false., p_rvector1, rcollection)
      
      ! Apply the source vector to the residual (if any)
      if (present(rsource)) then
        if (rsource%NEQ .gt. 0)&
            call lsysbl_vectorLinearComb(rsource, p_rvector1, 1.0_DP, 1.0_DP)
      end if
      
      ! Scale rvector1 by the inverse of the lumped mass matrix and store
      ! the result in rvector; this is the solution of the lumped version
      call lsysbl_invertedDiagMatVec(rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
          p_rvector1, 1.0_DP, rvector)
      
      ! Store norm of the initial guess from the lumped version
      dnorm0 = lsyssc_vectorNorm(rvector%RvectorBlock(1), LINALG_NORML2)
      
      richardson: do ite = 1, nmaxIterationsApproxTimeDerivative
        ! Initialise rvector2 by the constant right-hand side
        call lsysbl_copyVector(p_rvector1, p_rvector2)
        
        ! Compute the residual $rhs-M_C*u$ and store the result in rvector2
        call lsyssc_matVec(rproblemLevel%RmatrixScalar(consistentMassMatrix),&
            rvector%RvectorBlock(1), p_rvector2%RvectorBlock(1), -1.0_DP, 1.0_DP)
        
        ! Scale rvector2 by the inverse of the lumped mass matrix
        call lsysbl_invertedDiagMatVec(rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
            p_rvector2, 1.0_DP, p_rvector2)
        
        ! Apply solution increment (rvector2) to the previous solution iterate
        call lsysbl_vectorLinearComb(p_rvector2, rvector, 1.0_DP, 1.0_DP)
        
        ! Check for convergence
        dnorm = lsyssc_vectorNorm(p_rvector2%RvectorBlock(1), LINALG_NORML2)
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

      !-------------------------------------------------------------------------
      
    case(AFCSTAB_UPWIND)
      
      ! Make a backup copy of the stabilisation types because we
      ! have to overwrite them to enforce using the standard
      ! Galerkin scheme; this implies that their specification flags
      ! are changed, so make a backup copy of them, too
      if (convectionAFC .gt. 0) then
        cafcstabTypeConvection&
            = rproblemLevel%Rafcstab(convectionAFC)%cafcstabType
        istabilisationSpecConvection&
            = rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec
        rproblemLevel%Rafcstab(convectionAFC)%cafcstabType&
            = AFCSTAB_UPWIND
      end if
      
      if (diffusionAFC .gt. 0) then
        cafcstabTypeDiffusion&
            = rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType
        istabilisationSpecDiffusion&
            = rproblemLevel%Rafcstab(diffusionAFC)%istabilisationSpec
        rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType&
            = AFCSTAB_DMP
      end if
      
      ! Create discrete transport operator of low-order
      call transp_calcTransportOperator(rproblemLevel,&
          rsolver%rboundaryCondition, rsolution, rtimestep%dTime,&
          1.0_DP, p_rmatrix, ssectionName, rcollection,&
          bforceUpdate=bforceUpdate)
      
      ! Reset stabilisation structures for further usage
      if (convectionAFC .gt. 0) then
        rproblemLevel%Rafcstab(convectionAFC)%cafcstabType&
            = cafcstabTypeConvection
        rproblemLevel%Rafcstab(convectionAFC)%istabilisationSpec&
            = istabilisationSpecConvection
      end if
      
      if (diffusionAFC .gt. 0) then
        rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType&
            = cafcstabTypeDiffusion
        rproblemLevel%Rafcstab(diffusionAFC)%istabilisationSpec&
            = istabilisationSpecDiffusion
      end if
      
      ! Compute $L(u^L)*u^L$ and store the result in rvector
      call lsyssc_matVec(p_rmatrix,&
          rsolution%rvectorBlock(1),&
          rvector%RvectorBlock(1), 1.0_DP, 0.0_DP)
      
      ! Evaluate linear form for the boundary integral (if any)
      call transp_calcLinfBdrCondition(&
          rproblemLevel, rsolver%rboundaryCondition,&
          rsolution, ssectionName, smode, ivelocitytype,&
          rtimestep%dTime, 1.0_DP, .false., rvector, rcollection,&
          fcb_coeffVecBdrPrimal1d_sim, fcb_coeffVecBdrDual1d_sim,&
          fcb_coeffVecBdrPrimal2d_sim, fcb_coeffVecBdrDual2d_sim,&
          fcb_coeffVecBdrPrimal3d_sim, fcb_coeffVecBdrDual3d_sim)
      
      ! Build the geometric source term (if any)
      call transp_calcGeometricSourceterm(p_rparlist, ssectionName,&
          rproblemLevel, rsolution, 1.0_DP, .false., rvector, rcollection)
      
      ! Apply the source vector to the residual (if any)
      if (present(rsource)) then
        if (rsource%NEQ .gt. 0)&
            call lsysbl_vectorLinearComb(rsource, rvector, 1.0_DP, 1.0_DP)
      end if
      
      ! Scale it by the inverse of the lumped mass matrix
      call lsysbl_invertedDiagMatVec(rproblemLevel%RmatrixScalar(lumpedMassMatrix),&
          rvector, 1.0_DP, rvector)
      
    case default
      call output_line('Unsupported type of transport operator!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcTimeDerivative')
      call sys_halt()
    end select
       
    ! Release temporal memory
    if (.not.present(rmatrix)) then
      call lsyssc_releaseMatrix(rmatrix)
      deallocate(p_rmatrix)
    end if
    
  end subroutine transp_calcTimeDerivative

end module transport_callback

!##############################################################################
!# ****************************************************************************
!# <name> transport_callback </name>
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
!# 1.) transp_calcPrecondThetaScheme
!#     -> Calculates the nonlinear preconditioner
!#        used in the two-level theta-scheme
!#
!# 2.) transp_calcJacobianThetaScheme
!#     -> Calculates the Jacobian matrix
!#        used in the two-level theta-scheme
!#
!# 3.) transp_calcResidualThetaScheme
!#     -> Calculates the nonlinear residual vector
!#        used in the two-level theta-scheme
!#
!# 4.) transp_calcRhsThetaScheme
!#     -> Calculates the explicit right-hand side vector
!#        used in the two-level theta-scheme
!#
!# 5.) transp_calcRhsRungeKuttaScheme
!#     -> Calculates the right-hand side vector
!#        used in the explicit Runge-Kutta scheme
!#
!# 6.) transp_setBoundaryConditions
!#     -> Imposes boundary conditions for nonlinear solver
!#        by filtering the system matrix and the solution/residual
!#        vector explicitly (i.e. strong boundary conditions)
!#
!# 7.) transp_calcBilfBoundaryConditions
!#     -> Calculates the bilinear form arising from the weak
!#        imposition of boundary conditions
!#
!# 8.) transp_calcLinfBoundaryConditions
!#     -> Calculates the linear form arising from the weak
!#        imposition of boundary conditions 
!#
!# 9.) transp_calcVelocityField
!#     -> Calculates the velocity field
!#
!# 10.) transp_setVelocityField
!#      -> Sets the velocity field internally
!#
!# 11.) transp_calcLinearizedFCT
!#      -> Calculates the linearized FCT correction
!#
!# 12.) transp_coeffVectorAnalytic
!#      -> Callback routine for the evaluation of linear forms
!#         using an analytic expression for the load-vector
!#
!# 13.) transp_refFuncAnalytic
!#      -> Callback routine for the evaluation of the reference
!#         target function for goal-oriented error estimation
!#
!# 14.) transp_weightFuncAnalytic
!#      -> Callback routine for the evaluation of the weights in
!#         the target functional for goal-oriented error estimation
!#
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
!#        transp_calcResidualThetaScheme. Note that changing the
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

  use afcstabilisation
  use basicgeometry
  use bilinearformevaluation
  use boundary
  use boundaryfilter
  use collection
  use cubature
  use derivatives
  use dofmapping
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use groupfemscalar
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

  use transport_basic
  use transport_callback1d
  use transport_callback2d
  use transport_callback3d

  implicit none

  private
  public :: transp_calcLinfBoundaryConditions
  public :: transp_calcJacobianThetaScheme
  public :: transp_calcLinearizedFCT
  public :: transp_calcPrecondThetaScheme
  public :: transp_calcRhsRungeKuttaScheme
  public :: transp_calcRhsThetaScheme
  public :: transp_calcResidualThetaScheme
  public :: transp_calcVelocityField
  public :: transp_coeffVectorAnalytic
  public :: transp_nlsolverCallback
  public :: transp_refFuncAnalytic
  public :: transp_setBoundaryConditions
  public :: transp_setVelocityField
  public :: transp_weightFuncAnalytic

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
    integer(i32) :: iSpec
    integer :: jacobianMatrix
    
    
    !###########################################################################
    ! REMARK: The order in which the operations are performed is
    ! essential. This is due to the fact that the calculation of the
    ! residual/rhs requires the discrete transport operator to be
    ! initialized which is assembled in the calculation of the
    ! preconditioner. To prevent the re-assembly of the
    ! preconditioner twice, we remove the specifier
    ! NLSOL_OPSPEC_CALCPRECOND if the residual/rhs vector is built.
    !###########################################################################

    ! Make a local copy
    iSpec = ioperationSpec
    
    ! Do we have to calculate the constant right-hand side?
    ! --------------------------------------------------------------------------
    if ((iand(iSpec, NLSOL_OPSPEC_CALCRHS)  .ne. 0)) then

      ! Compute the preconditioner
      call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rcollection)

      ! Compute the right-hand side
      call transp_calcRhsRungeKuttaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rrhs, istep, rcollection,&
          rsource)

      ! Remove specifier for the preconditioner (if any)
      iSpec = iand(iSpec, not(NLSOL_OPSPEC_CALCPRECOND))
    end if
    
    
    ! Do we have to calculate the residual?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      if (istep .eq. 0) then
        ! Compute the constant right-hand side
        call transp_calcRhsThetaScheme(rproblemLevel, rtimestep,&
            rsolver, rsolution0, rrhs, rcollection, rsource)
      end if
      
      ! Compute the preconditioner
      call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rcollection)
      
      ! Compute the residual
      call transp_calcResidualThetaScheme(rproblemLevel, rtimestep, rsolver,&
          rsolution, rsolution0, rrhs, rres, istep, rcollection)

      ! Remove specifier for the preconditioner (if any)
      iSpec = iand(iSpec, not(NLSOL_OPSPEC_CALCPRECOND))
    end if
    

    ! Do we have to calculate the preconditioner?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) then
     
      ! Compute the preconditioner
      call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rcollection)
    end if

    
    ! Do we have to calculate the Jacobian operator?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_CALCJACOBIAN) .ne. 0) then
      
      ! Compute the Jacobian matrix
      call transp_calcJacobianThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rcollection)
    end if
    
    
    ! Do we have to impose boundary conditions?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then
      
      ! Impose boundary conditions
      call transp_setBoundaryConditions(rproblemLevel, rtimestep,&
          rsolver, rsolution, rsolution0, rres, rcollection)
    end if
    

    ! Do we have to apply the Jacobian operator?
    ! --------------------------------------------------------------------------
    if (iand(iSpec, NLSOL_OPSPEC_APPLYJACOBIAN) .ne. 0) then

      p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'jacobianMatrix', jacobianMatrix)

      ! Apply Jacobian matrix
      call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(jacobianMatrix),&
          rsolution%RvectorBlock(1), rres%RvectorBlock(1), 1.0_DP,&
          1.0_DP)
    end if
    
    
    ! Set status flag
    istatus = 0
    
  end subroutine transp_nlsolverCallback
  
  !*****************************************************************************

!<subroutine>

  subroutine transp_calcPrecondThetaScheme(rproblemLevel,&
      rtimestep, rsolver, rsolution, rcollection,&
      fcb_calcMatrixPrimal, fcb_calcMatrixDual,&
      fcb_coeffMatBdrPrimal_sim, fcb_coeffMatBdrDual_sim)

!<description>
    ! This subroutine calculates the nonlinear preconditioner for the
    !  primal/dual problem and configures the linear solver structure
    !  accordingly. The details of this subroutine are explained in code.
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! user-defined callback functions
    include 'intf_transpCalcMatrix.inc'
    optional :: fcb_calcMatrixPrimal
    optional :: fcb_calcMatrixDual

    ! user-defined callback functions
    include 'intf_transpCoeffMatBdr.inc'
    optional :: fcb_coeffMatBdrPrimal_sim
    optional :: fcb_coeffMatBdrDual_sim
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
    character(LEN=SYS_STRLEN) :: smode
    logical :: bStabilize, bbuildAFC, bconservative
    integer :: systemMatrix, transportMatrix, lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ, coeffMatrix_S
    integer :: imasstype, ivelocitytype, idiffusiontype
    integer :: convectionAFC, diffusionAFC, velocityfield

    !###########################################################################
    ! REMARK: If you want to add a new type of velocity/diffusion,
    ! then search for the tag @FAQ2: in this subroutine and create a
    ! new CASE which performs the corresponding task for the new type
    ! of velocity/ diffusion.
    !###########################################################################
    
    ! Check if the preconditioner has to be updated and return otherwise.
    if (iand(rproblemLevel%iproblemSpec, PROBLEV_MSPEC_UPDATE) .eq. 0) return
    
    ! Start time measurement for matrix evaluation
    p_rtimer => collct_getvalue_timer(rcollection, 'rtimerAssemblyMatrix')
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)
    
    ! Remove update notifier for further calls. Depending on the
    ! velocity, diffusion type it will be re-activited below.
    rproblemLevel%iproblemSpec = iand(rproblemLevel%iproblemSpec,&
                                      not(PROBLEV_MSPEC_UPDATE))

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'transportmatrix', transportMatrix)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'coeffMatrix_CX', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist,&
    rcollection%SquickAccess(1), 'coeffMatrix_CY', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'coeffMatrix_CZ', coeffMatrix_CZ)
    call parlst_getvalue_int(p_rparlist,&
        Rcollection%squickaccess(1), 'Coeffmatrix_s', coeffMatrix_S)
    
    !---------------------------------------------------------------------------
    ! Assemble diffusion operator:
    !
    ! $$ \int_\Omega \nabla w \cdot (D \nabla u) {\rm d}{\bf x} $$
    ! 
    ! The diffusion operator is symmetric so that it is the same for
    ! the primal and the dual problem. If no diffusion is present,
    ! i.e. $D \equiv 0$, then the transport operator is initialized by
    ! zeros. If there is anisotropic diffusion, i.e. $D=D({\bf x},t)$,
    ! then we may also initialize some stabilization structure.
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
        rcollection%SquickAccess(1), 'idiffusiontype', idiffusiontype)
    
    ! Primal and dual mode are equivalent
    ! @FAQ2: Which type of diffusion are we?
    select case(idiffusiontype)
    case (DIFFUSION_ZERO)
      ! zero diffusion, clear the transport matrix
      call lsyssc_clearMatrix(rproblemLevel%Rmatrix(transportMatrix))
      
    case (DIFFUSION_ISOTROPIC)
      ! Isotropic diffusion
      call gfsc_buildDiffusionOperator(rproblemLevel&
          %Rmatrix(coeffMatrix_S), .false., .true., rproblemLevel&
          %Rmatrix(transportMatrix))

    case (DIFFUSION_ANISOTROPIC)
      ! Anisotropic diffusion
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'diffusionAFC', diffusionAFC)

      if (diffusionAFC > 0) then
        
        ! What kind of stabilisation should be applied?
        select case(rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation)

        case (AFCSTAB_DMP)
          ! Satisfy discrete maximum principle
          call gfsc_buildDiffusionOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_S), .true., .true., rproblemLevel&
              %Rmatrix(transportMatrix))
          
        case (AFCSTAB_SYMMETRIC)
          ! Satisfy discrete maximum principle and assemble stabilization structure
          call gfsc_buildDiffusionOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_S), .true., .true., rproblemLevel&
              %Rmatrix(transportMatrix), rproblemLevel%Rafcstab(diffusionAFC))

        case default
          ! Compute the standard Galerkin approximation
          call gfsc_buildDiffusionOperator(rproblemLevel&
              %Rmatrix(coeffMatrix_S), .false., .true., rproblemLevel&
              %Rmatrix(transportMatrix))
        end select

      else   ! diffusionAFC < 0

        ! Compute the standard Galerkin approximation
        call gfsc_buildDiffusionOperator(rproblemLevel&
            %Rmatrix(coeffMatrix_S), .false., .true., rproblemLevel&
            %Rmatrix(transportMatrix))

      end if   ! diffusionAFC


    case (DIFFUSION_VARIABLE)
      call output_line('Variable diffusion matrices are yet not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcPrecondThetaScheme')
      call sys_halt()

      ! Set update notification in problem level structure
      rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                       PROBLEV_MSPEC_UPDATE)
      
    case default
      call output_line('Invalid type of diffusion!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcPrecondThetaScheme')
      call sys_halt()
    end select
    
    
    !---------------------------------------------------------------------------
    ! Assemble convective operator:
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
        rcollection%SquickAccess(1), 'mode', smode)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'convectionAFC', convectionAFC)

    if (convectionAFC > 0) then

      ! Check if stabilization should be applied
      bStabilize = AFCSTAB_GALERKIN .ne.&
                   rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation

      ! Check if stabilization structure should be built
      bbuildAFC = bStabilize .and. AFCSTAB_UPWIND .ne.&
                  rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation
      
    else   ! convectionAFC < 0

      bStabilize = .false.
      bbuildAFC  = .false.

    end if   ! convectionAFC

    ! Check if conservative or non-conservative formulation should be applied
    bconservative = (ivelocitytype .gt. 0)

    ! Set velocity vector (if any)
    if (transp_hasVelocityVector(ivelocityType)) then
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'velocityfield', velocityfield)
      call transp_setVelocityField(rproblemLevel%RvectorBlock(velocityfield))
    end if

    ! Are we in primal or dual mode?
    select case(trim(smode))
    case('primal')
      
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
        if (present(fcb_calcMatrixPrimal)) then
          
          if (bbuildAFC) then

            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsc_buildConvectionOperator(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rsolution, fcb_calcMatrixPrimal, .true., .false.,&
                  rproblemLevel%Rmatrix(transportMatrix),&
                  rproblemLevel%Rafcstab(convectionAFC), bconservative)
              
            case (NDIM2D)
              call gfsc_buildConvectionOperator(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rsolution, fcb_calcMatrixPrimal, .true., .false.,&
                  rproblemLevel%Rmatrix(transportMatrix),&
                  rproblemLevel%Rafcstab(convectionAFC), bconservative)
              
            case (NDIM3D)
              call gfsc_buildConvectionOperator(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rsolution, fcb_calcMatrixPrimal, .true., .false.,&
                  rproblemLevel%Rmatrix(transportMatrix),&
                  rproblemLevel%Rafcstab(convectionAFC), bconservative)
            end select

          else   ! bbuildAFC = false
            
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsc_buildConvectionOperator(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rsolution, fcb_calcMatrixPrimal, bStabilize,&
                  .false., rproblemLevel%Rmatrix(transportMatrix),&
                  bisConservative = bconservative)
              
            case (NDIM2D)
              call gfsc_buildConvectionOperator(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rsolution, fcb_calcMatrixPrimal, bStabilize,&
                  .false., rproblemLevel%Rmatrix(transportMatrix),&
                  bisConservative = bconservative)
              
            case (NDIM3D)
              call gfsc_buildConvectionOperator(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rsolution, fcb_calcMatrixPrimal, bStabilize,&
                  .false., rproblemLevel%Rmatrix(transportMatrix),&
                  bisConservative = bconservative)
            end select
            
          end if
          
        else ! callback function not present
          
          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcPrecondThetaScheme')
          call sys_halt()
          
        end if

        ! The user-defined callback function for matrix coefficients
        ! is used if present; otherwise an error is throws
        if (present(fcb_coeffMatBdrPrimal_sim)) then
          
          ! Evaluate bilinear form for boundary integral (if any)
          call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
              rtimestep%dTime, 1.0_DP, fcb_coeffMatBdrPrimal_sim,&
              rproblemLevel%Rmatrix(transportMatrix), rcollection,&
              BILF_MATC_LUMPED)

        else ! callback function not present
          
          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcPrecondThetaScheme')
          call sys_halt()
          
        end if

        
      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP)
        ! linear velocity

        if (bbuildAFC) then
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                rsolution, transp_calcMatrixPrimalConst1d,&
                .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
                rproblemLevel%Rafcstab(convectionAFC), bconservative)

          case (NDIM2D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rsolution, transp_calcMatrixPrimalConst2d,&
                .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
                rproblemLevel%Rafcstab(convectionAFC), bconservative)

          case (NDIM3D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                rsolution, transp_calcMatrixPrimalConst3d,&
                .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
                rproblemLevel%Rafcstab(convectionAFC), bconservative)
          end select
          
        else   ! bbuildAFC = false
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                rsolution, transp_calcMatrixPrimalConst1d,&
                bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix),&
                bisConservative = bconservative)

          case (NDIM2D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rsolution, transp_calcMatrixPrimalConst2d,&
                bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix),&
                bisConservative = bconservative)

          case (NDIM3D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                rsolution, transp_calcMatrixPrimalConst3d,&
                bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix),&
                bisConservative = bconservative)
          end select
          
        end if
        
        ! Evaluate bilinear form for boundary integral (if any)
        call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
            rtimestep%dTime, 1.0_DP, transp_coeffMatBdrPrimalConst2d,&
            rproblemLevel%Rmatrix(transportMatrix), rcollection,&
            BILF_MATC_LUMPED)

        if (abs(ivelocitytype) .eq. VELOCITY_TIMEDEP) then
          ! Set update notification in problem level structure
          rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                           PROBLEV_MSPEC_UPDATE)
        end if
        

      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers` equation in space-time

        if (bbuildAFC) then
          
          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsolution, transp_calcMatrixPrimalBurgersSpT2d,&
              .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
              rproblemLevel%Rafcstab(convectionAFC), bconservative)

        else   ! bbuildAFC = no

          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsolution, transp_calcMatrixPrimalBurgersSpT2d,&
              bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix),&
              bisConservative = bconservative)

        end if
        
        ! Evaluate bilinear form for boundary integral (if any)
        call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
            rtimestep%dTime, 1.0_DP,&
            transp_coeffMatBdrPrimalBurgersSpT2d, rproblemLevel&
            %Rmatrix(transportMatrix), rcollection, BILF_MATC_LUMPED)

        ! Set update notification in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         PROBLEV_MSPEC_UPDATE)


      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time

        if (bbuildAFC) then

          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsolution, transp_calcMatrixPrimalBuckLevSpT2d,&
              .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
              rproblemLevel%Rafcstab(convectionAFC), bconservative)

        else   ! bbuildAFC = no

          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsolution, transp_calcMatrixPrimalBuckLevSpT2d,&
              bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix),&
              bisConservative = bconservative)

        end if

        ! Evaluate bilinear form for boundary integral (if any)
        call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
            rtimestep%dTime, 1.0_DP,&
            transp_coeffMatBdrPrimalBuckLevSpT2d, rproblemLevel&
            %Rmatrix(transportMatrix), rcollection, BILF_MATC_LUMPED)

        ! Set update notification in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         PROBLEV_MSPEC_UPDATE)


      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers` equation in 1D

        if (bbuildAFC) then

          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rsolution, transp_calcMatrixPrimalBurgers1d,&
              .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
              rproblemLevel%Rafcstab(convectionAFC), bconservative)

        else   ! bbuildAFC = no

          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rsolution, transp_calcMatrixPrimalBurgers1d,&
              bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix),&
              bisConservative = bconservative)

        end if

        ! @TODO: Weak boundary conditions
        
        ! Set update notification in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         PROBLEV_MSPEC_UPDATE)
        

      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers` equation in 2D

        if (bbuildAFC) then

          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsolution, transp_calcMatrixPrimalBurgers2d,&
              .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
              rproblemLevel%Rafcstab(convectionAFC), bconservative)

        else   ! bbuildAFC = no

          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsolution, transp_calcMatrixPrimalBurgers2d,&
              bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix),&
              bisConservative = bconservative)

        end if

        ! Evaluate bilinear form for boundary integral (if any)
        call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
            rtimestep%dTime, 1.0_DP,&
            transp_coeffMatBdrPrimalBurgers2d, rproblemLevel&
            %Rmatrix(transportMatrix), rcollection, BILF_MATC_LUMPED)

        ! Set update notification in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         PROBLEV_MSPEC_UPDATE)


      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D

        if (bbuildAFC) then

          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rsolution, transp_calcMatrixPrimalBuckLev1d,&
              .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
              rproblemLevel%Rafcstab(convectionAFC), bconservative)

        else   ! bbuildAFC = no

          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rsolution, transp_calcMatrixPrimalBuckLev1d,&
              bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix),&
              bisConservative = bconservative)

        end if

        ! @TODO: Weak boundary conditions

        ! Set update notification in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         PROBLEV_MSPEC_UPDATE)
      end select

      
    case('dual')

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
        !  is used if present; otherwise an error is thrown
        if (present(fcb_calcMatrixDual)) then
          
          if (bbuildAFC) then
            
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsc_buildConvectionOperator(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rsolution, fcb_calcMatrixDual, .true., .false.,&
                  rproblemLevel%Rmatrix(transportMatrix),&
                  rproblemLevel%Rafcstab(convectionAFC), bconservative)
              
            case (NDIM2D)
              call gfsc_buildConvectionOperator(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rsolution, fcb_calcMatrixDual, .true., .false.,&
                  rproblemLevel%Rmatrix(transportMatrix),&
                  rproblemLevel%Rafcstab(convectionAFC), bconservative)
              
            case (NDIM3D)
              call gfsc_buildConvectionOperator(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rsolution, fcb_calcMatrixDual, .true., .false.,&
                  rproblemLevel%Rmatrix(transportMatrix),&
                  rproblemLevel%Rafcstab(convectionAFC), bconservative)
            end select

          else   ! bbuildAFC = false
            
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsc_buildConvectionOperator(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rsolution, fcb_calcMatrixDual, bStabilize,&
                  .false., rproblemLevel%Rmatrix(transportMatrix),&
                  bisConservative = bconservative)
              
            case (NDIM2D)
              call gfsc_buildConvectionOperator(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rsolution, fcb_calcMatrixDual, bStabilize,&
                  .false., rproblemLevel%Rmatrix(transportMatrix),&
                  bisConservative = bconservative)
              
            case (NDIM3D)
              call gfsc_buildConvectionOperator(&
                  rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rsolution, fcb_calcMatrixDual, bStabilize,&
                  .false., rproblemLevel%Rmatrix(transportMatrix),&
                  bisConservative = bconservative)
            end select
            
          end if
          
        else ! callback function not present
          
          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcPrecondThetaScheme')
          call sys_halt()
          
        end if

        ! The user-defined callback function for matrix coefficients
        !  is used if present; otherwise an error is throws
        if (present(fcb_coeffMatBdrDual_sim)) then
          
          ! Evaluate bilinear form for boundary integral (if any)
          call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
              rtimestep%dTime, -1.0_DP, fcb_coeffMatBdrDual_sim,&
              rproblemLevel%Rmatrix(transportMatrix), rcollection,&
              BILF_MATC_LUMPED)

        else ! callback function not present
          
          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcPrecondThetaScheme')
          call sys_halt()
          
        end if


      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP)
        ! linear velocity

        if (bbuildAFC) then

          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                rsolution, transp_calcMatrixDualConst1d,&
                .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
                rproblemLevel%Rafcstab(convectionAFC), bconservative)

          case (NDIM2D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rsolution, transp_calcMatrixDualConst2d,&
                .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
                rproblemLevel%Rafcstab(convectionAFC), bconservative)

          case (NDIM3D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                rsolution, transp_calcMatrixDualConst3d,&
                .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
                rproblemLevel%Rafcstab(convectionAFC), bconservative)
          end select

        else   ! bbuildAFC = false
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                rsolution, transp_calcMatrixDualConst1d,&
                bstabilize, .false., rproblemLevel%Rmatrix(transportMatrix),&
                bisConservative = bconservative)

          case (NDIM2D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rsolution, transp_calcMatrixDualConst2d,&
                bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix),&
                bisConservative = bconservative)

          case (NDIM3D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                rsolution, transp_calcMatrixDualConst3d,&
                bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix),&
                bisConservative = bconservative)
          end select

        end if

        ! Evaluate bilinear form for boundary integral (if any)
        call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
            rtimestep%dTime, -1.0_DP, transp_coeffMatBdrDualConst2d,&
            rproblemLevel%Rmatrix(transportMatrix), rcollection,&
            BILF_MATC_LUMPED)
        
        if (abs(ivelocitytype) .eq. VELOCITY_TIMEDEP) then
          ! Set update notification in problem level structure
          rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                       PROBLEV_MSPEC_UPDATE)
        end if
        
        ! @TODO: The dual mode has only been implemented for linear
        ! convection. If you need to compute the dual problem for
        ! some other velocity type, then you have to add the
        ! implementation of the dual transport operator below!

      end select

    case DEFAULT
      call output_line('Invalid mode!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcPrecondThetaScheme')
      call sys_halt()
    end select
        
    !---------------------------------------------------------------------------
    ! Assemble the global system operator
    !---------------------------------------------------------------------------
    
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'systemmatrix', systemMatrix)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'imasstype', imasstype)

    select case(imasstype)
    case (MASS_LUMPED)

      !-------------------------------------------------------------------------
      ! Compute the global operator for transient flow
      !
      !   $ A = ML-theta*dt*L $
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1),&
          'lumpedmassmatrix', lumpedMassMatrix)

      call lsyssc_MatrixLinearComb(rproblemLevel&
          %Rmatrix(lumpedMassMatrix), 1.0_DP, rproblemLevel&
          %Rmatrix(transportMatrix), -rtimestep%theta*rtimestep%dStep,&
          rproblemLevel%Rmatrix(systemMatrix), .false., .false.,&
          .true., .true.)

    case (MASS_CONSISTENT)

      !-------------------------------------------------------------------------
      ! Compute the global operator for transient flow
      !
      !   $ A = MC-theta*dt*L $
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1),&
          'consistentmassmatrix', consistentMassMatrix)

      call lsyssc_MatrixLinearComb(rproblemLevel&
          %Rmatrix(consistentMassMatrix), 1.0_DP, rproblemLevel&
          %Rmatrix(transportMatrix), -rtimestep%theta*rtimestep%dStep,&
          rproblemLevel%Rmatrix(systemMatrix), .false., .false.,&
          .true., .true.)

    case DEFAULT

      !-------------------------------------------------------------------------
      ! Compute the global operator for steady-state flow
      !
      !   $ A = -L $
      !-------------------------------------------------------------------------
      
      call lsyssc_copyMatrix(rproblemLevel%Rmatrix(transportMatrix),&
          rproblemLevel%Rmatrix(systemMatrix))
      call lsyssc_scaleMatrix(rproblemLevel%Rmatrix(systemMatrix), -1.0_DP)
      
    end select
    
    ! Impose boundary conditions in strong sence (if any)
    if (rsolver%rboundaryCondition%bStrongBdrCond) then
      call bdrf_filterMatrix(rsolver%rboundaryCondition,&
          rproblemLevel%Rmatrix(systemMatrix), 1.0_DP)
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
    
  end subroutine transp_calcPrecondThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcJacobianThetaScheme(rproblemLevel,&
      rtimestep, rsolver, rsolution, rsolution0, rcollection,&
      fcb_calcMatrixPrimal, fcb_calcMatrixDual,&
      fcb_coeffMatBdrPrimal_sim, fcb_coeffMatBdrDual_sim)

!<description>
    ! This callback subroutine computes the Jacobian matrix.
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolution0

    ! user-defined callback functions
    include 'intf_transpCalcMatrix.inc'
    optional :: fcb_calcMatrixPrimal
    optional :: fcb_calcMatrixDual

    ! user-defined callback functions
    include 'intf_transpCoeffMatBdr.inc'
    optional :: fcb_coeffMatBdrPrimal_sim
    optional :: fcb_coeffMatBdrDual_sim
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    real(DP) :: hstep
    character(LEN=SYS_STRLEN) :: smode
    integer :: transportMatrix, jacobianMatrix
    integer :: consistentMassMatrix, lumpedMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ, coeffMatrix_S
    integer :: convectionAFC, diffusionAFC, ijacobianFormat
    integer :: imasstype, imassantidiffusiontype, ivelocitytype, idiffusiontype
    integer :: velocityfield
    logical :: bStabilize, bisExactStructure, bisExtendedSparsity

    !###########################################################################
    ! REMARK: If you want to add a new type of velocity/diffusion,
    ! then search for the tag @FAQ2: in this subroutine and create a
    ! new CASE which performs the corresponding task for the new type
    ! of velocity/ diffusion.
    !###########################################################################
    
    ! Start time measurement for matrix evaluation
    p_rtimer => collct_getvalue_timer(rcollection, 'rtimerAssemblyMatrix')
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'transportmatrix', transportMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'jacobianmatrix', jacobianMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'coeffMatrix_CX', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'coeffMatrix_CY', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'coeffMatrix_CZ', coeffMatrix_CZ)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'coeffMatrix_S', coeffMatrix_S)
    
    ! The Jacobian matrix for the low-order transport operator needs
    ! to be generated only in case of nonlinear governing equations.
    ! In this case, the corresponding transport operator L has to be
    ! updated in each nonlinear iteration. Hence, we can temporarily
    ! build the Jacobian matrix using the memory of the operator L.

    ! Compute step lenth of the solution perturbation
    select case(int(rsolver%p_solverNewton%dperturbationStrategy))
    case (PERTURB_NITSOL)
      ! Choice h=[(1+|u|)*EPS]^{1/(1+p)} by Pernice et al. in
      ! M. Pernice, H.F. Walker, NITSOL: a Newton iterative solver
      ! for nonlinear systems, SIAM J. Sci. Comput. 19 (1998) 302-318.
      hstep = ( (1+lsysbl_vectorNorm(rsolution,&
                   LINALG_NORMEUCLID))*SYS_EPSREAL )**(1.0_DP/3._DP)
      
    case (PERTURB_SQRTEPS)
      hstep= sqrt(SYS_EPSREAL)

    case DEFAULT
      hstep = max(SYS_EPSREAL,&
                  rsolver%p_solverNewton%dperturbationStrategy)
    end select
    

    !---------------------------------------------------------------------------
    ! Assemble diffusion operator:
    !
    ! $$ \int_\Omega \nabla w \cdot (D \nabla u) {\rm d}{\bf x} $$
    ! 
    ! The diffusion operator is symmetric so that it is the same for
    ! the primal and the dual problem. If no diffusion is present,
    ! i.e. $D \equiv 0$, then the transport operator is initialized by
    ! zeros. If there is anisotropic diffusion, i.e. $D=D({\bf x},t)$,
    ! then we may also initialize some stabilization structure.
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
        rcollection%SquickAccess(1), 'idiffusiontype', idiffusiontype)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'diffusionAFC', diffusionAFC)

    ! @FAQ2: What type of diffusion are we?
    select case(idiffusiontype)
    case (DIFFUSION_ZERO)
      ! zero diffusion, clear the system matrix
      call lsyssc_clearMatrix(rproblemLevel%Rmatrix(transportMatrix))
      
    case (DIFFUSION_ISOTROPIC,&
          DIFFUSION_ANISOTROPIC)
      ! Isotropic diffusion
      call gfsc_buildDiffusionOperator(rproblemLevel&
          %Rmatrix(coeffMatrix_S), .false., .true., rproblemLevel&
          %Rmatrix(transportMatrix))
      
    case (DIFFUSION_VARIABLE)
      print *, "Variable diffusion matrices are yet not implemented!"
      stop
      
    case DEFAULT
      call output_line('Invalid type of diffusion!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobianThetaScheme')
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
        rcollection%SquickAccess(1), 'mode', smode)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'convectionAFC', convectionAFC)

    if (convectionAFC > 0) then

      ! Check if stabilization should be applied
      bStabilize = (AFCSTAB_GALERKIN .ne.&
                    rproblemLevel%Rafcstab(convectionAFC)&
                    %ctypeAFCstabilisation)

    else   ! convectionAFC < 0

      bStabilize = .false.

    end if   ! convectionAFC

    ! Set velocity vector (if any)
    if (transp_hasVelocityVector(ivelocityType)) then
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'velocityfield', velocityfield)
      call transp_setVelocityField(rproblemLevel%RvectorBlock(velocityfield))
    end if

    ! Are we in primal or dual mode?
    select case(trim(smode))
    case('primal')

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
        if (present(fcb_calcMatrixPrimal)) then
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildConvectionJacobian(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                fcb_calcMatrixPrimal, hstep, bStabilize,&
                .false., rproblemLevel%Rmatrix(transportMatrix))
            
          case (NDIM2D)
            call gfsc_buildConvectionJacobian(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                fcb_calcMatrixPrimal, hstep, bStabilize,&
                .false., rproblemLevel%Rmatrix(transportMatrix))
            
          case (NDIM3D)
            call gfsc_buildConvectionJacobian(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                fcb_calcMatrixPrimal, hstep, bStabilize,&
                .false., rproblemLevel%Rmatrix(transportMatrix))
          end select
          
        else ! callback function not present
          
          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobianThetaScheme')
          call sys_halt()

        end if

        ! The user-defined callback function for matrix coefficients
        ! is used if present; otherwise an error is throws
        if (present(fcb_coeffMatBdrPrimal_sim)) then
          
          ! Evaluate bilinear form for boundary integral (if any)
          call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
              rtimestep%dTime, 1.0_DP, fcb_coeffMatBdrPrimal_sim,&
              rproblemLevel%Rmatrix(transportMatrix), rcollection,&
              BILF_MATC_LUMPED)

        else ! callback function not present
          
          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobianThetaScheme')
          call sys_halt()
          
        end if


      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP) 
        ! linear velocity
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsc_buildConvectionJacobian(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              transp_calcMatrixPrimalConst1d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))

        case (NDIM2D)
          call gfsc_buildConvectionJacobian(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              transp_calcMatrixPrimalConst2d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))

        case (NDIM3D)
          call gfsc_buildConvectionJacobian(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              transp_calcMatrixPrimalConst3d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))
        end select

        ! Evaluate bilinear form for boundary integral (if any)
        call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
            rtimestep%dTime, 1.0_DP, transp_coeffMatBdrPrimalConst2d,&
            rproblemLevel%Rmatrix(transportMatrix), rcollection,&
            BILF_MATC_LUMPED)
      

      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers` equation in space-time
        call gfsc_buildConvectionJacobian(rproblemLevel&
            %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
            transp_calcMatrixPrimalBurgersSpT2d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))
        
        ! Evaluate bilinear form for boundary integral (if any)
        call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
            rtimestep%dTime, 1.0_DP,&
            transp_coeffMatBdrPrimalBurgersSpT2d, rproblemLevel&
            %Rmatrix(transportMatrix), rcollection, BILF_MATC_LUMPED)


      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time
        call gfsc_buildConvectionJacobian(rproblemLevel&
            %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
            transp_calcMatrixPrimalBuckLevSpT2d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))

        ! Evaluate bilinear form for boundary integral (if any)
        call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
            rtimestep%dTime, 1.0_DP,&
            transp_coeffMatBdrPrimalBuckLevSpT2d, rproblemLevel&
            %Rmatrix(transportMatrix), rcollection, BILF_MATC_LUMPED)
        

      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers` equation in 1D
        call gfsc_buildConvectionJacobian(rproblemLevel&
            %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
            transp_calcMatrixPrimalBurgers1d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))

        ! @TODO: Add weak boundary conditions
        

      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers` equation in 2D
        call gfsc_buildConvectionJacobian(rproblemLevel&
            %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
            transp_calcMatrixPrimalBurgers2d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))

        ! Evaluate bilinear form for boundary integral (if any)
        call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
            rtimestep%dTime, 1.0_DP,&
            transp_coeffMatBdrPrimalBurgers2d, rproblemLevel&
            %Rmatrix(transportMatrix), rcollection, BILF_MATC_LUMPED)

        
      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D
        call gfsc_buildConvectionJacobian(rproblemLevel&
            %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
            transp_calcMatrixPrimalBuckLev1d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))

        ! @TODO: Add weak boundary conditions

      end select


    case('dual')
      
      ! @FAQ2: What type of velocity are we?
      select case (abs(ivelocitytype))
      case default
        ! The user-defined callback function for matrix coefficients
        ! is used if present; otherwise an error is thrown
        if (present(fcb_calcMatrixDual)) then
          
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildConvectionJacobian(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                fcb_calcMatrixDual, hstep, bStabilize,&
                .false., rproblemLevel%Rmatrix(transportMatrix))
            
          case (NDIM2D)
            call gfsc_buildConvectionJacobian(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                fcb_calcMatrixDual, hstep, bStabilize,&
                .false., rproblemLevel%Rmatrix(transportMatrix))
            
          case (NDIM3D)
            call gfsc_buildConvectionJacobian(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                fcb_calcMatrixDual, hstep, bStabilize,&
                .false., rproblemLevel%Rmatrix(transportMatrix))
          end select
          
        else ! callback function not present
          
          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobianThetaScheme')
          call sys_halt()
          
        end if

        ! The user-defined callback function for matrix coefficients
        ! is used if present; otherwise an error is throws
        if (present(fcb_coeffMatBdrDual_sim)) then
          
          ! Evaluate bilinear form for boundary integral (if any)
          call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
              rtimestep%dTime, 1.0_DP, fcb_coeffMatBdrDual_sim,&
              rproblemLevel%Rmatrix(transportMatrix), rcollection,&
              BILF_MATC_LUMPED)

        else ! callback function not present
          
          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobianThetaScheme')
          call sys_halt()
          
        end if


      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        
      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP) 
        ! linear velocity
        
        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsc_buildConvectionJacobian(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              transp_calcMatrixDualConst1d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))

        case (NDIM2D)
          call gfsc_buildConvectionJacobian(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              transp_calcMatrixDualConst2d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))

        case (NDIM3D)
          call gfsc_buildConvectionJacobian(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
              transp_calcMatrixDualConst3d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))
        end select

        ! Evaluate bilinear form for boundary integral (if any)
        call transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
            rtimestep%dTime, 1.0_DP, transp_coeffMatBdrDualConst2d,&
            rproblemLevel%Rmatrix(transportMatrix), rcollection,&
            BILF_MATC_LUMPED)
        
        ! @TODO: The dual mode has only been implemented for linear
        ! convection. If you need to compute the dual problem for
        ! some other velocity type, then you have to add the
        ! implementation of the dual transport operator below!

      end select


    case default
      call output_line('Invalid mode!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobianThetaScheme')
      call sys_halt()
    end select
    
    
    ! Check if the Jacobian operator has extended sparsity pattern
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'ijacobianFormat', ijacobianFormat)
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
        rcollection%SquickAccess(1), 'imasstype', imasstype)

    select case(imasstype)
    case (MASS_LUMPED)
      
      !-------------------------------------------------------------------------
      ! Compute the global Jacobian for transient flow
      !
      !   $ J = ML-theta*dt*L $
      !-------------------------------------------------------------------------
      
      call lsyssc_MatrixLinearComb(rproblemLevel&
          %Rmatrix(transportMatrix), -rtimestep%theta*rtimestep%dStep,&
          rproblemLevel%Rmatrix(lumpedMassMatrix), 1.0_DP,&
          rproblemLevel%Rmatrix(jacobianMatrix), .false., .false.,&
          .true., bisExactStructure)

    case (MASS_CONSISTENT)

      !-------------------------------------------------------------------------
      ! Compute the global Jacobian for transient flow
      !
      !   $ J = MC-theta*dt*L $
      !-------------------------------------------------------------------------
      
      call lsyssc_MatrixLinearComb(rproblemLevel&
          %Rmatrix(transportMatrix), -rtimestep%theta*rtimestep%dStep,&
          rproblemLevel%Rmatrix(consistentMassMatrix), 1.0_DP,&
          rproblemLevel%Rmatrix(jacobianMatrix), .false., .false.,&
          .true., bisExactStructure)

    case DEFAULT

      !-------------------------------------------------------------------------
      ! Compute the global Jacobian for steady-state flow
      !
      !   $ J = -L $
      !-------------------------------------------------------------------------
      
      call lsyssc_MatrixLinearComb(rproblemLevel&
          %Rmatrix(transportMatrix), -1.0_DP, rproblemLevel&
          %Rmatrix(jacobianMatrix), 0.0_DP, rproblemLevel&
          %Rmatrix(jacobianMatrix), .false., .false., .true.,&
          bisExactStructure)

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
      select case(rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_SYMMETRIC)
        call gfsc_buildJacobianSymm(rsolution, 1.0_DP, hstep, .false.,&
            rproblemLevel%Rafcstab(diffusionAFC), rproblemLevel&
            %Rmatrix(jacobianMatrix), bisExtendedSparsity)
      end select

    case (DIFFUSION_VARIABLE)
      print *, "Variable diffusion matrices are yet not implemented!"
      stop
      
    case DEFAULT
      call output_line('Invalid type of diffusion!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobianThetaScheme')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Assemble the Jacobian matrix for the convection operator
    !---------------------------------------------------------------------------
    
    ! Are we in primal or dual mode?
    select case(trim(smode))
    case('primal')
      
      !-------------------------------------------------------------------------
      ! We are in primal mode
      !-------------------------------------------------------------------------

      ! @FAQ2: What type of velocity are we?
      select case(abs(ivelocitytype))
      case default
        ! The user-defined callback function for matrix coefficients
        ! is used if present; otherwise an error is thrown
        if (present(fcb_calcMatrixPrimal)) then
          
          select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
          case (AFCSTAB_FEMFCT,&
                AFCSTAB_FEMFCT_CLASSICAL)
            
            call parlst_getvalue_int(p_rparlist,&
                rcollection%SquickAccess(1),&
                'imassantidiffusiontype', imassantidiffusiontype)
            
            ! Should we apply consistent mass antidiffusion?
            if (imassantidiffusiontype .eq. MASS_CONSISTENT) then

              select case(rproblemLevel%rtriangulation%ndim)
              case (NDIM1D)
                call gfsc_buildJacobianFCT(rproblemLevel&
                    %Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                    rsolution, fcb_calcMatrixPrimal, rtimestep&
                    %theta, rtimestep%dStep, hstep, .false.,&
                    rproblemLevel%Rafcstab(convectionAFC),&
                    rproblemLevel%Rmatrix(jacobianMatrix),&
                    rproblemLevel%Rmatrix(consistentMassMatrix))

              case (NDIM2D)
                call gfsc_buildJacobianFCT(rproblemLevel&
                    %Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                    rsolution, fcb_calcMatrixPrimal, rtimestep&
                    %theta, rtimestep%dStep, hstep, .false.,&
                    rproblemLevel%Rafcstab(convectionAFC),&
                    rproblemLevel%Rmatrix(jacobianMatrix),&
                    rproblemLevel%Rmatrix(consistentMassMatrix))

              case (NDIM3D)
                call gfsc_buildJacobianFCT(rproblemLevel&
                    %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                    rsolution, fcb_calcMatrixPrimal, rtimestep&
                    %theta, rtimestep%dStep, hstep, .false.,&
                    rproblemLevel%Rafcstab(convectionAFC),&
                    rproblemLevel%Rmatrix(jacobianMatrix),&
                    rproblemLevel%Rmatrix(consistentMassMatrix))
              end select
              
            else
              
              select case(rproblemLevel%rtriangulation%ndim)
              case (NDIM1D)
                call gfsc_buildJacobianFCT(rproblemLevel&
                    %Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                    rsolution, fcb_calcMatrixPrimal, rtimestep&
                    %theta, rtimestep%dStep, hstep, .false.,&
                    rproblemLevel%Rafcstab(convectionAFC),&
                    rproblemLevel%Rmatrix(jacobianMatrix))

              case (NDIM2D)
                call gfsc_buildJacobianFCT(rproblemLevel&
                    %Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                    rsolution, fcb_calcMatrixPrimal, rtimestep&
                    %theta, rtimestep%dStep, hstep, .false.,&
                    rproblemLevel%Rafcstab(convectionAFC),&
                    rproblemLevel%Rmatrix(jacobianMatrix))
                
              case (NDIM3D)
                call gfsc_buildJacobianFCT(rproblemLevel&
                    %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                    rsolution, fcb_calcMatrixPrimal, rtimestep&
                    %theta, rtimestep%dStep, hstep, .false.,&
                    rproblemLevel%Rafcstab(convectionAFC),&
                    rproblemLevel%Rmatrix(jacobianMatrix))
              end select

            end if
            
          case (AFCSTAB_FEMTVD)
            
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsc_buildJacobianTVD(rproblemLevel&
                  %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                  fcb_calcMatrixPrimal, rtimestep%dStep, hstep,&
                  .false., rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%Rmatrix(jacobianMatrix),&
                  bisExtendedSparsity)

            case (NDIM2D)
              call gfsc_buildJacobianTVD(rproblemLevel&
                  %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                  fcb_calcMatrixPrimal, rtimestep%dStep, hstep,&
                  .false., rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%Rmatrix(jacobianMatrix),&
                  bisExtendedSparsity)
              
            case (NDIM3D)
              call gfsc_buildJacobianTVD(rproblemLevel&
                  %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                  fcb_calcMatrixPrimal, rtimestep%dStep, hstep,&
                  .false., rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%Rmatrix(jacobianMatrix),&
                  bisExtendedSparsity)
            end select
            
          case (AFCSTAB_FEMGP)
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsc_buildJacobianGP(rproblemLevel&
                  %Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                  rsolution, rsolution0, fcb_calcMatrixPrimal,&
                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%Rmatrix(jacobianMatrix),&
                  bisExtendedSparsity)
              
            case (NDIM2D)
              call gfsc_buildJacobianGP(rproblemLevel&
                  %Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                  rsolution, rsolution0, fcb_calcMatrixPrimal,&
                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%Rmatrix(jacobianMatrix),&
                  bisExtendedSparsity)
              
            case (NDIM3D)
              call gfsc_buildJacobianGP(rproblemLevel&
                  %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                  rsolution, rsolution0, fcb_calcMatrixPrimal,&
                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%Rmatrix(jacobianMatrix),&
                  bisExtendedSparsity)
            end select
            
          end select

        else ! callback function not present
          
          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobianThetaScheme')
          call sys_halt()
          
        end if


      case (VELOCITY_ZERO)
        ! zero velocity, do nothing

      
      case(VELOCITY_CONSTANT,&
           VELOCITY_TIMEDEP) 
        ! linear velocity
        
        select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
        case (AFCSTAB_FEMFCT,&
             AFCSTAB_FEMFCT_CLASSICAL)
          
          call parlst_getvalue_int(p_rparlist,&
              rcollection%SquickAccess(1),&
              'imassantidiffusiontype', imassantidiffusiontype)
          
          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call gfsc_buildJacobianFCT(rsolution, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel&
                %Rmatrix(jacobianMatrix), rproblemLevel&
                %Rmatrix(consistentMassMatrix))
          else
            call gfsc_buildJacobianFCT(rsolution, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel%Rmatrix(jacobianMatrix))
          end if
          
          
        case (AFCSTAB_FEMTVD)
          call gfsc_buildJacobianTVD(rsolution, rtimestep%dStep, hstep,&
              .false., rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%Rmatrix(jacobianMatrix),&
              bisExtendedSparsity)
          
        case (AFCSTAB_FEMGP)
          call gfsc_buildJacobianGP(rproblemLevel&
              %Rmatrix(consistentMassMatrix), rsolution,&
              rsolution0, rtimestep%theta, rtimestep%dStep, hstep,&
              .false., rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%Rmatrix(jacobianMatrix),&
              bisExtendedSparsity)
        end select
        
        
      case(VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers` equation in space-time
        
        select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
        case (AFCSTAB_FEMFCT,&
              AFCSTAB_FEMFCT_CLASSICAL)
          
          call parlst_getvalue_int(p_rparlist,&
              rcollection%SquickAccess(1),&
              'imassantidiffusiontype', imassantidiffusiontype)
          
          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call gfsc_buildJacobianFCT(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                transp_calcMatrixPrimalBurgersSpT2d, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel&
                %Rmatrix(jacobianMatrix), rproblemLevel&
                %Rmatrix(consistentMassMatrix))
          else
            call gfsc_buildJacobianFCT(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                transp_calcMatrixPrimalBurgersSpT2d, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel&
                %Rmatrix(jacobianMatrix))
          end if
          
        case (AFCSTAB_FEMTVD)
          call gfsc_buildJacobianTVD(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              transp_calcMatrixPrimalBurgersSpT2d, rtimestep%dStep,&
              hstep, .false., rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%Rmatrix(jacobianMatrix),&
              bisExtendedSparsity)
          
        case (AFCSTAB_FEMGP)
          call gfsc_buildJacobianGP(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rproblemLevel&
              %Rmatrix(consistentMassMatrix), rsolution,&
              rsolution0, transp_calcMatrixPrimalBurgersSpT2d,&
              rtimestep%theta, rtimestep%dStep, hstep, .false.,&
              rproblemLevel%Rafcstab(convectionAFC), rproblemLevel&
              %Rmatrix(jacobianMatrix), bisExtendedSparsity)
        end select
        
        
      case(VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time
        
        select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
        case (AFCSTAB_FEMFCT,&
              AFCSTAB_FEMFCT_CLASSICAL)
          
          call parlst_getvalue_int(p_rparlist,&
              rcollection%SquickAccess(1),&
              'imassantidiffusiontype', imassantidiffusiontype)
          
          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call gfsc_buildJacobianFCT(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                transp_calcMatrixPrimalBuckLevSpT2d, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel&
                %Rmatrix(jacobianMatrix), rproblemLevel&
                %Rmatrix(consistentMassMatrix))
          else
            call gfsc_buildJacobianFCT(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                transp_calcMatrixPrimalBuckLevSpT2d, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel&
                %Rmatrix(jacobianMatrix))
          end if
          
        case (AFCSTAB_FEMTVD)
          call gfsc_buildJacobianTVD(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              transp_calcMatrixPrimalBuckLevSpT2d, rtimestep%dStep,&
              hstep, .false., rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%Rmatrix(jacobianMatrix),&
              bisExtendedSparsity)
          
        case (AFCSTAB_FEMGP)
          call gfsc_buildJacobianGP(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rproblemLevel&
              %Rmatrix(consistentMassMatrix), rsolution,&
              rsolution0, transp_calcMatrixPrimalBuckLevSpT2d,&
              rtimestep%theta, rtimestep%dStep, hstep, .false.,&
              rproblemLevel%Rafcstab(convectionAFC), rproblemLevel&
              %Rmatrix(jacobianMatrix), bisExtendedSparsity)
        end select
        
        
      case(VELOCITY_BURGERS1D)
        ! nonlinear Burgers` equation in 1D
        
        select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
        case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)
          
          call parlst_getvalue_int(p_rparlist,&
              rcollection%SquickAccess(1),&
              'imassantidiffusiontype', imassantidiffusiontype)
          
          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call gfsc_buildJacobianFCT(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                transp_calcMatrixPrimalBurgers1d, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel&
                %Rmatrix(jacobianMatrix), rproblemLevel&
                %Rmatrix(consistentMassMatrix))
          else
            call gfsc_buildJacobianFCT(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                transp_calcMatrixPrimalBurgers1d, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel&
                %Rmatrix(jacobianMatrix))
          end if
          
        case (AFCSTAB_FEMTVD)
          call gfsc_buildJacobianTVD(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              transp_calcMatrixPrimalBurgers1d, rtimestep%dStep, hstep,&
              .false., rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%Rmatrix(jacobianMatrix),&
              bisExtendedSparsity)
          
        case (AFCSTAB_FEMGP)
          call gfsc_buildJacobianGP(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rproblemLevel&
              %Rmatrix(consistentMassMatrix), rsolution,&
              rsolution0, transp_calcMatrixPrimalBurgers1d,&
              rtimestep%theta, rtimestep%dStep, hstep, .false.,&
              rproblemLevel%Rafcstab(convectionAFC), rproblemLevel&
              %Rmatrix(jacobianMatrix), bisExtendedSparsity)
        end select
        
        
      case(VELOCITY_BURGERS2D)
        ! nonlinear Burgers` equation in 2D
        
        select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
        case (AFCSTAB_FEMFCT,&
              AFCSTAB_FEMFCT_CLASSICAL)
          
          call parlst_getvalue_int(p_rparlist,&
              rcollection%SquickAccess(1),&
              'imassantidiffusiontype', imassantidiffusiontype)
          
          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call gfsc_buildJacobianFCT(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                transp_calcMatrixPrimalBurgers2d, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel&
                %Rmatrix(jacobianMatrix), rproblemLevel&
                %Rmatrix(consistentMassMatrix))
          else
            call gfsc_buildJacobianFCT(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                transp_calcMatrixPrimalBurgers2d, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel&
                %Rmatrix(jacobianMatrix))
          end if
          
        case (AFCSTAB_FEMTVD)
          call gfsc_buildJacobianTVD(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
              transp_calcMatrixPrimalBurgers2d, rtimestep%dStep, hstep,&
              .false., rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%Rmatrix(jacobianMatrix),&
              bisExtendedSparsity)
          
        case (AFCSTAB_FEMGP)
          call gfsc_buildJacobianGP(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rproblemLevel&
              %Rmatrix(consistentMassMatrix), rsolution,&
              rsolution0, transp_calcMatrixPrimalBurgers2d,&
              rtimestep%theta, rtimestep%dStep, hstep, .false.,&
              rproblemLevel%Rafcstab(convectionAFC), rproblemLevel&
              %Rmatrix(jacobianMatrix), bisExtendedSparsity)
        end select
        
        
      case(VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D
        
        select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
        case (AFCSTAB_FEMFCT,&
              AFCSTAB_FEMFCT_CLASSICAL)
          
          call parlst_getvalue_int(p_rparlist,&
              rcollection%SquickAccess(1),&
              'imassantidiffusiontype', imassantidiffusiontype)
          
          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call gfsc_buildJacobianFCT(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                transp_calcMatrixPrimalBuckLev1d, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel&
                %Rmatrix(jacobianMatrix), rproblemLevel&
                %Rmatrix(consistentMassMatrix))
          else
            call gfsc_buildJacobianFCT(rproblemLevel&
                %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                transp_calcMatrixPrimalBuckLev1d, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel&
                %Rmatrix(jacobianMatrix))
          end if
          
        case (AFCSTAB_FEMTVD)
          call gfsc_buildJacobianTVD(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
              transp_calcMatrixPrimalBuckLev1d, rtimestep%dStep, hstep,&
              .false., rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%Rmatrix(jacobianMatrix),&
              bisExtendedSparsity)
          
        case (AFCSTAB_FEMGP)
          call gfsc_buildJacobianGP(rproblemLevel&
              %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rproblemLevel&
              %Rmatrix(consistentMassMatrix), rsolution,&
              rsolution0, transp_calcMatrixPrimalBuckLev1d,&
              rtimestep%theta, rtimestep%dStep, hstep, .false.,&
              rproblemLevel%Rafcstab(convectionAFC), rproblemLevel&
              %Rmatrix(jacobianMatrix), bisExtendedSparsity)
        end select

      end select


    case ('dual')

      !-------------------------------------------------------------------------
      ! We are in dual mode
      !-------------------------------------------------------------------------

      ! @FAQ2: What type of velocity are we?
      select case(abs(ivelocitytype))
      case default
        ! The user-defined callback function for matrix coefficients
        ! is used if present; otherwise an error is thrown
        if (present(fcb_calcMatrixDual)) then
          
          select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
          case (AFCSTAB_FEMFCT,&
                AFCSTAB_FEMFCT_CLASSICAL)
            
            call parlst_getvalue_int(p_rparlist,&
                rcollection%SquickAccess(1),&
                'imassantidiffusiontype', imassantidiffusiontype)
            
            ! Should we apply consistent mass antidiffusion?
            if (imassantidiffusiontype .eq. MASS_CONSISTENT) then

              select case(rproblemLevel%rtriangulation%ndim)
              case (NDIM1D)
                call gfsc_buildJacobianFCT(rproblemLevel&
                    %Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                    rsolution, fcb_calcMatrixDual, rtimestep&
                    %theta, rtimestep%dStep, hstep, .false.,&
                    rproblemLevel%Rafcstab(convectionAFC),&
                    rproblemLevel%Rmatrix(jacobianMatrix),&
                    rproblemLevel%Rmatrix(consistentMassMatrix))

              case (NDIM2D)
                call gfsc_buildJacobianFCT(rproblemLevel&
                    %Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                    rsolution, fcb_calcMatrixDual, rtimestep&
                    %theta, rtimestep%dStep, hstep, .false.,&
                    rproblemLevel%Rafcstab(convectionAFC),&
                    rproblemLevel%Rmatrix(jacobianMatrix),&
                    rproblemLevel%Rmatrix(consistentMassMatrix))

              case (NDIM3D)
                call gfsc_buildJacobianFCT(rproblemLevel&
                    %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                    rsolution, fcb_calcMatrixDual, rtimestep&
                    %theta, rtimestep%dStep, hstep, .false.,&
                    rproblemLevel%Rafcstab(convectionAFC),&
                    rproblemLevel%Rmatrix(jacobianMatrix),&
                    rproblemLevel%Rmatrix(consistentMassMatrix))
              end select
              
            else
              
              select case(rproblemLevel%rtriangulation%ndim)
              case (NDIM1D)
                call gfsc_buildJacobianFCT(rproblemLevel&
                    %Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                    rsolution, fcb_calcMatrixDual, rtimestep&
                    %theta, rtimestep%dStep, hstep, .false.,&
                    rproblemLevel%Rafcstab(convectionAFC),&
                    rproblemLevel%Rmatrix(jacobianMatrix))

              case (NDIM2D)
                call gfsc_buildJacobianFCT(rproblemLevel&
                    %Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                    rsolution, fcb_calcMatrixDual, rtimestep&
                    %theta, rtimestep%dStep, hstep, .false.,&
                    rproblemLevel%Rafcstab(convectionAFC),&
                    rproblemLevel%Rmatrix(jacobianMatrix))
                
              case (NDIM3D)
                call gfsc_buildJacobianFCT(rproblemLevel&
                    %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                    rsolution, fcb_calcMatrixDual, rtimestep&
                    %theta, rtimestep%dStep, hstep, .false.,&
                    rproblemLevel%Rafcstab(convectionAFC),&
                    rproblemLevel%Rmatrix(jacobianMatrix))
              end select

            end if
            
          case (AFCSTAB_FEMTVD)
            
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsc_buildJacobianTVD(rproblemLevel&
                  %Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
                  fcb_calcMatrixDual, rtimestep%dStep, hstep,&
                  .false., rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%Rmatrix(jacobianMatrix),&
                  bisExtendedSparsity)

            case (NDIM2D)
              call gfsc_buildJacobianTVD(rproblemLevel&
                  %Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
                  fcb_calcMatrixDual, rtimestep%dStep, hstep,&
                  .false., rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%Rmatrix(jacobianMatrix),&
                  bisExtendedSparsity)
              
            case (NDIM3D)
              call gfsc_buildJacobianTVD(rproblemLevel&
                  %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ), rsolution,&
                  fcb_calcMatrixDual, rtimestep%dStep, hstep,&
                  .false., rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%Rmatrix(jacobianMatrix),&
                  bisExtendedSparsity)
            end select
            
          case (AFCSTAB_FEMGP)
            select case(rproblemLevel%rtriangulation%ndim)
            case (NDIM1D)
              call gfsc_buildJacobianGP(rproblemLevel&
                  %Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                  rsolution, rsolution0, fcb_calcMatrixDual,&
                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%Rmatrix(jacobianMatrix),&
                  bisExtendedSparsity)
              
            case (NDIM2D)
              call gfsc_buildJacobianGP(rproblemLevel&
                  %Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                  rsolution, rsolution0, fcb_calcMatrixDual,&
                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%Rmatrix(jacobianMatrix),&
                  bisExtendedSparsity)
              
            case (NDIM3D)
              call gfsc_buildJacobianGP(rproblemLevel&
                  %Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                  rsolution, rsolution0, fcb_calcMatrixDual,&
                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                  rproblemLevel%Rafcstab(convectionAFC),&
                  rproblemLevel%Rmatrix(jacobianMatrix),&
                  bisExtendedSparsity)
            end select
            
          end select

        else ! callback function not present
          
          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobianThetaScheme')
          call sys_halt()
          
        end if

      case (VELOCITY_ZERO)
        ! zero velocity, do nothing

      
      case(VELOCITY_CONSTANT,&
           VELOCITY_TIMEDEP) 
        ! linear velocity
        
        select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
        case (AFCSTAB_FEMFCT,&
             AFCSTAB_FEMFCT_CLASSICAL)
          
          call parlst_getvalue_int(p_rparlist,&
              rcollection%SquickAccess(1),&
              'imassantidiffusiontype', imassantidiffusiontype)
          
          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call gfsc_buildJacobianFCT(rsolution, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel&
                %Rmatrix(jacobianMatrix), rproblemLevel&
                %Rmatrix(consistentMassMatrix))
          else
            call gfsc_buildJacobianFCT(rsolution, rtimestep%theta,&
                rtimestep%dStep, hstep, .false., rproblemLevel&
                %Rafcstab(convectionAFC), rproblemLevel%Rmatrix(jacobianMatrix))
          end if
          
          
        case (AFCSTAB_FEMTVD)
          call gfsc_buildJacobianTVD(rsolution, rtimestep%dStep, hstep,&
              .false., rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%Rmatrix(jacobianMatrix),&
              bisExtendedSparsity)
          
        case (AFCSTAB_FEMGP)
          call gfsc_buildJacobianGP(rproblemLevel&
              %Rmatrix(consistentMassMatrix), rsolution,&
              rsolution0, rtimestep%theta, rtimestep%dStep, hstep,&
              .false., rproblemLevel%Rafcstab(convectionAFC),&
              rproblemLevel%Rmatrix(jacobianMatrix),&
              bisExtendedSparsity)
        end select

      end select

    case DEFAULT
      call output_line('Invalid mode!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobianThetaScheme')
      call sys_halt()
    end select
    

    ! Impose boundary conditions in strong sence (if any)
    if (rsolver%rboundaryCondition%bStrongBdrCond) then
      call bdrf_filterMatrix(rsolver%rboundaryCondition,&
          rproblemLevel%Rmatrix(jacobianMatrix), 1.0_DP)
    end if

    ! Ok, we updated the Jacobian matrix successfully. Now we still have to
    ! link it to the solver hierarchy. This is done recursively.

    ! What type of flow are we?
    select case(imasstype)
    case (MASS_LUMPED,&
          MASS_CONSISTENT)
      call flagship_updateSolverMatrix(rproblemLevel, rsolver,&
          jacobianMatrix, SYSTEM_INTERLEAVEFORMAT,&
          UPDMAT_JAC_TRANSIENT, rproblemLevel%ilev, rproblemLevel&
          %ilev)

    case DEFAULT
      call flagship_updateSolverMatrix(rproblemLevel, rsolver,&
          jacobianMatrix, SYSTEM_INTERLEAVEFORMAT, UPDMAT_JAC_STEADY,&
          rproblemLevel%ilev, rproblemLevel%ilev)
    end select
    
    ! Finally, we have to update the content of the solver hierarchy
    call solver_updateContent(rsolver)
    
    ! Stop time measurement for matrix evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine transp_calcJacobianThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcRhsRungeKuttaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rrhs, istep, rcollection,&
      rsource, fcb_coeffVecBdrPrimal_sim, fcb_coeffVecBdrDual_sim)

!<description>
    ! This subroutine computes the right-hand side vector
    ! used in the explicit Runge-Kutta scheme.
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(in) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolution0

    ! number of explicit step
    integer, intent(in) :: istep

    ! OPTIONAL: load vector specified by the application
    type(t_vectorBlock), intent(in), optional :: rsource

    ! OPTIONAL: user-defined callback functions
    include 'intf_transpCoeffVecBdr.inc'
    optional :: fcb_coeffVecBdrPrimal_sim
    optional :: fcb_coeffVecBdrDual_sim
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! right-hand side vector
    type(t_vectorBlock), intent(inout) :: rrhs

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    real(DP) :: dweight
    integer :: transportMatrix
    integer :: consistentMassMatrix, lumpedMassMatrix
    integer :: convectionAFC, diffusionAFC
    integer :: imassantidiffusiontype


    print *, "WARNING: This subroutine has not been tested!"
    stop

    ! Start time measurement for residual/rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection, 'rtimerAssemblyVector')
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'transportmatrix', transportMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'consistentmassmatrix', consistentMassMatrix)


    ! Compute the right-hand side
    !
    !   $$ rhs = weight*(1-theta)*dt*L(u)*u $$

    dweight = rtimestep%DmultistepWeights(istep)*&
              rtimestep%dStep*(1.0_DP-rtimestep%theta)
    call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
        rsolution%rvectorBlock(1), rrhs%RvectorBlock(1), dweight,&
        0.0_DP)


    ! Perform algebraic flux correction for the convective term if required
    !
    !   $$ rhs = rhs + f^*(u^n+1,u^n) $$

    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'convectionAFC', convectionAFC)
       
    ! What kind of stabilisation should be applied?
    select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
          
    case (AFCSTAB_FEMFCT,&
          AFCSTAB_FEMFCT_CLASSICAL)

      dweight = rtimestep%DmultistepWeights(istep)*rtimestep%dStep
      call parlst_getvalue_int(p_rparlist, rcollection&
          %SquickAccess(1), 'imassantidiffusiontype',&
          imassantidiffusiontype)

      ! Should we apply consistent mass antidiffusion?
      if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
        call gfsc_buildResidualFCT(rproblemLevel&
            %Rmatrix(lumpedMassMatrix), rsolution, rtimestep%theta,&
            dweight, .true., rrhs, rproblemLevel&
            %Rafcstab(convectionAFC), rproblemLevel&
            %Rmatrix(consistentMassMatrix))
      else
        call gfsc_buildResidualFCT(rproblemLevel&
            %Rmatrix(lumpedMassMatrix), rsolution, rtimestep%theta,&
            dweight, .true., rrhs, rproblemLevel&
            %Rafcstab(convectionAFC))
      end if
          
    case (AFCSTAB_FEMTVD)
      call gfsc_buildResidualTVD(rsolution, dweight, rrhs,&
          rproblemLevel%Rafcstab(convectionAFC))

    case (AFCSTAB_FEMGP)
      call gfsc_buildResidualGP(rproblemLevel&
          %Rmatrix(consistentMassMatrix), rsolution, rsolution0&
          , rtimestep%theta, dweight, rrhs, rproblemLevel&
          %Rafcstab(convectionAFC))
    end select


    ! Perform algebraic flux correction for the diffusive term if required
    !
    !   $$ rhs = rhs + g^*(u^n+1,u^n) $$

    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'diffusionAFC', diffusionAFC)
    
    ! What kind of stabilisation should be applied?
    select case(rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation)
          
    case (AFCSTAB_SYMMETRIC)
      call gfsc_buildResidualSymm(rsolution, 1.0_DP, rrhs,&
          rproblemLevel%Rafcstab(diffusionAFC))
    end select
    
    ! Apply the given load vector to the residual
    if (present(rsource))&
        call lsysbl_vectorLinearComb(rsource, rrhs, 1.0_DP, 1.0_DP)

    ! Stop time measurement for rhs evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine transp_calcRhsRungeKuttaScheme
    
  !*****************************************************************************

!<subroutine>

  subroutine transp_calcRhsThetaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rrhs, rcollection, rsource,&
      fcb_coeffVecBdrPrimal_sim, fcb_coeffVecBdrDual_sim)

!<description>
    ! This subroutine computes the constant right-hand side
    !
    !  $$ rhs = [M + (1-\theta)\Delta t K^n]u^n + s^n$$
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
    optional :: fcb_coeffVecBdrPrimal_sim
    optional :: fcb_coeffVecBdrDual_sim
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
    character(LEN=SYS_STRLEN) :: smode
    real(DP) :: dscale
    integer :: consistentMassMatrix, lumpedMassMatrix
    integer :: transportMatrix, massMatrix
    integer :: imasstype, ivelocitytype


    !###########################################################################
    ! REMARK: If you want to add a new type of velocity/diffusion,
    ! then search for the tag @FAQ2: in this subroutine and create a
    ! new CASE which performs the corresponding task for the new type
    ! of velocity/ diffusion.
    !###########################################################################
    
    ! Check if the preconditioner has to be initialized
    if (iand(rproblemLevel%iproblemSpec,&
             PROBLEV_MSPEC_INITIALIZE) .ne. 0) then
      call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rcollection)
      rproblemLevel%iproblemSpec = iand(rproblemLevel%iproblemSpec,&
                                        not(PROBLEV_MSPEC_INITIALIZE))
    end if

    !---------------------------------------------------------------------------
    ! We calculate the constant right-hand side vector based on the
    ! transport operator. If boundary conditions are imposed in weak
    ! sense, then additional surface integrals are applied to the
    ! right-hand side of the primal problem. In the dual problem,
    ! weakly imposed boundary conditions are built into the target
    ! functional which is supplied via the vector 'rb'.
    ! ---------------------------------------------------------------------------

    ! Start time measurement for rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection, 'rtimerAssemblyVector')
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)
    
    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist, rcollection&
        %SquickAccess(1), 'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection&
        %SquickAccess(1), 'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection &
        %SquickAccess(1), 'transportmatrix', transportMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection&
        %SquickAccess(1), 'imasstype', imasstype)
    call parlst_getvalue_int(p_rparlist, rcollection&
        %SquickAccess(1), 'ivelocitytype', ivelocitytype)
    call parlst_getvalue_string(p_rparlist,&
        rcollection%SquickAccess(1), 'mode', smode)

    ! Do we have some kind of mass matrix?
    select case(imasstype)
    case (MASS_LUMPED, MASS_CONSISTENT)
      
      !-------------------------------------------------------------------------
      ! Compute the constant right-hand side
      !
      !   $$ rhs = M*u^n+(1-theta)*dt*K(u^n)u^n + b.c.'s $$
      !-------------------------------------------------------------------------
      
      ! Do we have some explicit part?
      if (rtimestep%theta .lt. 1.0_DP) then
      
        ! Compute scaling parameter
        dscale = (1.0_DP-rtimestep%theta) * rtimestep%dStep
        
        ! Build transport term $(1-theta)*dt*K(u^n)u^n$, where
        ! $T(u^n)$ denotes the discrete transport operator of high or
        ! low order evaluated at the old solution values
        call lsyssc_scalarMatVec(rproblemLevel&
            %Rmatrix(transportMatrix), rsolution%rvectorBlock(1),&
            rrhs%RvectorBlock(1), dscale, 0.0_DP)
        
        ! Build transient term $M_L*u^n$ or $M_C*u^n$
        massMatrix = merge(lumpedMassMatrix, consistentMassMatrix,&
                           imasstype .eq. MASS_LUMPED)

        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(massMatrix),&
            rsolution%RvectorBlock(1),&
            rrhs%RvectorBlock(1), 1.0_DP, 1.0_DP)
        
        ! Evaluate bilinear form for boundary integral (if any)
        
        ! --- explicit part ---
        call transp_calcLinfBdrCondQuick(rproblemLevel, rsolver,&
            smode, ivelocitytype, rtimestep%dTime-rtimestep%dStep,&
            -dscale, rrhs%RvectorBlock(1), rcollection,&
            fcb_coeffVecBdrPrimal_sim, fcb_coeffVecBdrDual_sim)
        
        dscale = rtimestep%theta*rtimestep%dStep

        ! --- implicit part ---
        call transp_calcLinfBdrCondQuick(rproblemLevel, rsolver,&
            smode, ivelocitytype, rtimestep%dTime,&
            -dscale, rrhs%RvectorBlock(1), rcollection,&
            fcb_coeffVecBdrPrimal_sim, fcb_coeffVecBdrDual_sim)
        
      else ! theta = 1
        
        ! Build transient term $M_L*u^n$ or $M_C*u^n$
        massMatrix = merge(lumpedMassMatrix, consistentMassMatrix,&
                           imasstype .eq. MASS_LUMPED)

        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(massMatrix),&
            rsolution%RvectorBlock(1),&
            rrhs%RvectorBlock(1), 1.0_DP, 0.0_DP)

        ! Evaluate bilinear form for boundary integral (if any)
        
        dscale = rtimestep%theta*rtimestep%dStep

        ! --- implicit part ---
        call transp_calcLinfBdrCondQuick(rproblemLevel, rsolver,&
            smode, ivelocitytype, rtimestep%dTime,&
            -dscale, rrhs%RvectorBlock(1), rcollection,&
            fcb_coeffVecBdrPrimal_sim, fcb_coeffVecBdrDual_sim)
        
      end if ! theta
      
    case DEFAULT
      
      !-------------------------------------------------------------------------
      ! Initialize the constant right-hand side by zeros
      !
      !   $$ rhs = "0" + b.c.'s $$
      !-------------------------------------------------------------------------
      
      ! Clear right-hand side vector
      call lsysbl_clearVector(rrhs)
      
      ! Evaluate bilinear form for boundary integral (if any)
      call transp_calcLinfBdrCondQuick(rproblemLevel, rsolver,&
          smode, ivelocitytype, rtimestep%dTime, -1.0_DP,&
          rrhs%RvectorBlock(1), rcollection,&
            fcb_coeffVecBdrPrimal_sim, fcb_coeffVecBdrDual_sim)
      
    end select
    
    
    ! Apply the source vector to the right-hand side (if any)
    if (present(rsource))&
        call lsysbl_vectorLinearComb(rsource, rrhs, 1.0_DP, 1.0_DP)
    
    ! Stop time measurement for rhs evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine transp_calcRhsThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcResidualThetaScheme(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rrhs, rres, ite, rcollection,&
      rsource)

!<description>
    ! This subroutine computes the nonlinear residual vector
    ! 
    !   $$ res^{(m)} = rhs - [M-\theta\Delta t K^{(m)}]u^{(m)} - s^{(m)} $$
    !
    ! for the standard two-level theta-scheme, whereby the (scaled)
    ! source term is optional. The constant right-hand side
    !
    !  $$ rhs = [M + (1-\theta)\Delta t K^n]u^n + s^n$$
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

    ! right-hand side vector
    type(t_vectorBlock), intent(in) :: rrhs

    ! iteration number
    integer, intent(in) :: ite

    ! OPTIONAL: source vector
    type(t_vectorBlock), intent(in), optional :: rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! residual vector
    type(t_vectorBlock), intent(inout) :: rres

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    real(DP) :: dscale
    integer :: transportMatrix, massMatrix
    integer :: consistentMassMatrix, lumpedMassMatrix
    integer :: convectionAFC, diffusionAFC
    integer :: imasstype, imassantidiffusiontype, ivelocitytype


    !###########################################################################
    ! REMARK: If you want to add a new type of velocity/diffusion,
    ! then search for the tag @FAQ2: in this subroutine and create a
    ! new CASE which performs the corresponding task for the new type
    ! of velocity/ diffusion.
    !###########################################################################
    
    ! Check if the preconditioner has to be initialized
    if (iand(rproblemLevel%iproblemSpec,&
             PROBLEV_MSPEC_INITIALIZE) .ne. 0) then
      call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
          rsolver, rsolution, rcollection)
      rproblemLevel%iproblemSpec = iand(rproblemLevel%iproblemSpec,&
                                        not(PROBLEV_MSPEC_INITIALIZE))
    end if

    !---------------------------------------------------------------------------
    ! We calculate the nonlinear residual vector based on the
    ! transport operator. If boundary conditions are imposed in weak
    ! sense, then they give rise to additional surface integrals in
    ! the bilinear form, i.e., the transport operator so that no
    ! additional linear forms for boundary conditions are required
    !---------------------------------------------------------------------------
    
    ! Start time measurement for residual evaluation
    p_rtimer => collct_getvalue_timer(rcollection, 'rtimerAssemblyVector')
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)
    
    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'transportmatrix', transportMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'imasstype', imasstype)
    
    ! Do we have some kind of mass matrix?
    select case(imasstype)
    case (MASS_LUMPED, MASS_CONSISTENT)

      !-------------------------------------------------------------------------
      ! Compute the residual for transient flows
      !
      !   $$ res^{(m)} = rhs - [M - dt*theta*K(u^{(m)})]*u^{(m)} $$
      !-------------------------------------------------------------------------

      ! Compute scaling parameter
      dscale = rtimestep%theta*rtimestep%dStep
      
      ! Apply constant right-hand side
      call lsysbl_copyVector(rrhs, rres)

      ! Apply transport operator
      call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
          rsolution%rvectorBlock(1), rres%RvectorBlock(1), dscale, 1.0_DP)
      
      massMatrix = merge(lumpedMassMatrix, consistentMassMatrix,&
                         imasstype .eq. MASS_LUMPED)

      ! Apply mass matrix
      call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(massMatrix),&
          rsolution%RvectorBlock(1), rres%RvectorBlock(1), -1.0_DP, 1.0_DP)
      
    case DEFAULT
      
      !-------------------------------------------------------------------------
      ! Compute the residual for stationary flows
      !
      !   $$ res^{(m)} = rhs + K(u^{(m)})*u^{(m)} $$
      !-------------------------------------------------------------------------

      ! Apply constant right-hand side
      call lsysbl_copyVector(rrhs, rres)
      
      ! Apply transport operator
      call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
          rsolution%rvectorBlock(1), rres%RvectorBlock(1), 1.0_DP, 1.0_DP)

    end select
    
    !-------------------------------------------------------------------------
    ! Perform algebraic flux correction for the convective term (if required)
    !
    !   $$ res = res + f^*(u^(m),u^n) $$
    !-------------------------------------------------------------------------
    
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'convectionAFC', convectionAFC)
    
    if (convectionAFC > 0) then
      
      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)
        
        call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
            'imassantidiffusiontype', imassantidiffusiontype)
        
        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsc_buildResidualFCT(rproblemLevel&
              %Rmatrix(lumpedMassMatrix), rsolution, rtimestep&
              %theta, rtimestep%dStep, (ite .eq. 0), rres,&
              rproblemLevel %Rafcstab(convectionAFC), rproblemLevel &
              %Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildResidualFCT(rproblemLevel&
              %Rmatrix(lumpedMassMatrix), rsolution, rtimestep&
              %theta, rtimestep%dStep, (ite .eq. 0), rres,&
              rproblemLevel %Rafcstab(convectionAFC))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildResidualTVD(rsolution, rtimestep%dStep, rres,&
            rproblemLevel%Rafcstab(convectionAFC))
        
      case (AFCSTAB_FEMGP)
        call gfsc_buildResidualGP(rproblemLevel&
            %Rmatrix(consistentMassMatrix), rsolution,&
            rsolution0, rtimestep%theta, rtimestep%dStep,&
            rres, rproblemLevel%Rafcstab(convectionAFC))
      end select
      
    end if   ! convectionAFC > 0

    !-------------------------------------------------------------------------
    ! Perform algebraic flux correction for the diffusive term (if required)
    !
    !   $$ res = res + g^*(u^n+1,u^n) $$
    !-------------------------------------------------------------------------
    
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'diffusionAFC', diffusionAFC)
    
    if (diffusionAFC > 0) then

      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation)
        
      case (AFCSTAB_SYMMETRIC)
        call gfsc_buildResidualSymm(rsolution, 1.0_DP, rres,&
            rproblemLevel%Rafcstab(diffusionAFC))
      end select
      
    end if   ! diffusionAFC > 0
    
    
    ! Apply the source vector to the residual (if any)
    if (present(rsource))&
        call lsysbl_vectorLinearComb(rsource, rres, -1.0_DP, 1.0_DP)
    
    ! Stop time measurement for residual evaluation
    call stat_stopTimer(p_rtimer)
    
  end subroutine transp_calcResidualThetaScheme

  !*****************************************************************************

!<subroutine>

  subroutine transp_setBoundaryConditions(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rres, rcollection)

!<description>
    ! This subroutine imposes the Dirichlet boundary conditions in 
    ! strong sense by filtering the system matrix, the solution
    !  vector and/or the residual vector explicitly.
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! nonlinear solver structure
    type(t_solver), intent(in) :: rsolver

    ! initial solution vector
    type(t_vectorBlock), intent(in) :: rsolution0
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
    integer :: imatrix
    
    ! Get parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

    ! What type of preconditioner are we?
    select case(rsolver%iprecond)
    case (NLSOL_PRECOND_BLOCKD,&
          NLSOL_PRECOND_DEFCOR,&
          NLSOL_PRECOND_NEWTON_FAILED)

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'systemmatrix', imatrix)
      
    case (NLSOL_PRECOND_NEWTON)

      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'jacobianmatrix', imatrix)
      
    case DEFAULT
      call output_line('Invalid nonlinear preconditioner!',&
          OU_CLASS_ERROR, OU_MODE_STD,'transp_setBoundaryConditions')
      call sys_halt()
    end select
    
    ! Impose boundary conditions for the solution vector and impose
    ! zeros in the residual vector and the off-diagonal positions of
    ! the system matrix which depends on the nonlinear solver
    call bdrf_filterSolution(rsolver%rboundaryCondition,&
        rproblemLevel%Rmatrix(imatrix), rsolution, rres,&
        rsolution0, rtimestep%dTime)

  end subroutine transp_setBoundaryConditions

  !*****************************************************************************

!<subroutine>
  
  subroutine transp_calcBilfBoundaryConditions(rproblemLevel, rsolver,&
      dtime, dscale, fcoeff_buildMatrixScBdr2D_sim, rmatrix,&
      rcollection, cconstrType)
    
!<description>
    ! This subroutine computes the bilinear form arising from the weak
    ! imposition of boundary conditions. For this application only
    ! weakly imposed Dirichlet boundary conditions give rise to
    ! additional contributions to the bilinear form.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling factor
    real(DP), intent(in) :: dscale

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

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_boundaryCondition), pointer :: p_rboundaryCondition
    type(t_parlist), pointer :: p_rparlist
    type(t_collection) :: rcollectionTmp
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_bilinearform) :: rform
    integer, dimension(:), pointer :: p_IbdrCondCpIdx, p_IbdrCondType
    integer :: ivelocitytype, velocityfield
    integer :: ibct, isegment

    ! Evaluate bilinear form for boundary integral and return if
    ! there are no weak boundary conditions available
    p_rboundaryCondition => rsolver%rboundaryCondition
    if (.not.p_rboundaryCondition%bWeakBdrCond) return

    ! Initialize temporal collection structure
    call collct_init(rcollectionTmp)

    ! Attach function parser from boundary conditions to collection
    ! structure and specify its name in quick access string array
    call collct_setvalue_pars(rcollectionTmp, 'rfparser',&
        p_rboundaryCondition%rfparser, .true.)
    rcollectionTmp%SquickAccess(1) = 'rfparser'

    ! Get parameters from parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'ivelocitytype', ivelocitytype)

    ! Attach velocity vector to temporal collection structure (if any)
    if (transp_hasVelocityVector(ivelocityType)) then
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'velocityfield', velocityfield)
      rcollectionTmp%p_rvectorQuickAccess1 =>&
          rproblemLevel%RvectorBlock(velocityfield)
    end if
    
    ! Initialize the bilinear form
    rform%itermCount = 1
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    
    ! We have no constant coefficients
    rform%ballCoeffConstant = .false.
    rform%BconstantCoeff    = .false.
    
    ! Prepare quick access arrays of temporal collection structure
    rcollectionTmp%DquickAccess(1) = dtime
    rcollectionTmp%DquickAccess(2) = dscale


    ! How many spatial dimensions are we?
    select case(rproblemLevel%rtriangulation%ndim)
    case (NDIM2D)
      ! Set pointers
      call storage_getbase_int(p_rboundaryCondition%h_IbdrCondCpIdx,&
          p_IbdrCondCpIdx)
      call storage_getbase_int(p_rboundaryCondition%h_IbdrCondType,&
          p_IbdrCondType)

      ! Loop over all boundary components
      do ibct = 1, p_rboundaryCondition%iboundarycount
        
        ! Loop over all boundary segments
        do isegment = p_IbdrCondCpIdx(ibct),&
                      p_IbdrCondCpIdx(ibct+1)-1
      
          ! Prepare quick access array of temporal collection structure
          rcollectionTmp%IquickAccess(1) = p_IbdrCondType(isegment)
          rcollectionTmp%IquickAccess(2) = isegment

          ! What type of boundary conditions are we?
          select case(p_IbdrCondType(isegment))
          case (BDR_DIRICHLET_WEAK)
            ! Create boundary region
            call bdrf_createRegion(p_rboundaryCondition,&
                ibct, isegment-p_IbdrCondCpIdx(ibct)+1,&
                rboundaryRegion)
                  
            ! Assemble the bilinear form
            call bilf_buildMatrixScalarBdr2D(rform, CUB_G3_1D,&
                .false., rmatrix, fcoeff_buildMatrixScBdr2D_sim,&
                rboundaryRegion, rcollectionTmp, cconstrType)
          end select
          
        end do ! isegment
      end do ! ibct
      
    case default
      call output_line('Unsupported spatial dimension!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcBilfBoundaryConditions')
      call sys_halt()
    end select
    
    ! Release temporal collection structure
    call collct_done(rcollectionTmp)
    
  end subroutine transp_calcBilfBoundaryConditions

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcLinfBoundaryConditions(rproblemLevel, rsolver,&
      dtime, dscale, fcoeff_buildVectorScBdr2D_sim, rvector, rcollection)

!<description>
    ! This subroutine computes the linear form arising from the weak
    ! imposition of boundary conditions. For this application only
    ! weakly imposed Dirichlet and non-homogeneous Neumann boundary
    ! conditions give rise to contributions to the linear form.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(in) :: rsolver
    
    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling factor
    real(DP), intent(in) :: dscale

    ! callback routine for nonconstant coefficient vectors.
    include '../../../../../kernel/DOFMaintenance/intf_coefficientVectorScBdr2D.inc'
!</intput>

!<inputoutput>
    ! residual/right-hand side vector
    type(t_vectorScalar), intent(inout) :: rvector

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_boundaryCondition), pointer :: p_rboundaryCondition
    type(t_parlist), pointer :: p_rparlist
    type(t_collection) :: rcollectionTmp
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_linearForm) :: rform
    integer, dimension(:), pointer :: p_IbdrCondCpIdx, p_IbdrCondType
    integer :: ivelocitytype, velocityfield
    integer :: ibct, isegment
    
    ! Evaluate linear form for boundary integral and return if
    ! there are no weak boundary conditions available
    p_rboundaryCondition => rsolver%rboundaryCondition
    if (.not.p_rboundaryCondition%bWeakBdrCond) return

    ! Initialize temporal collection structure
    call collct_init(rcollectionTmp)
    
    ! Attach function parser from boundary conditions to collection
    ! structure and specify its name in quick access string array
    call collct_setvalue_pars(rcollectionTmp, 'rfparser',&
        p_rboundaryCondition%rfparser, .true.)
    rcollectionTmp%SquickAccess(1) = 'rfparser'
    
    ! Get parameters from parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(1), 'ivelocitytype', ivelocitytype)

    ! Attach velocity vector to temporal collection structure (if any)
    if (transp_hasVelocityVector(ivelocityType)) then
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(1), 'velocityfield', velocityfield)
      rcollectionTmp%p_rvectorQuickAccess1 =>&
          rproblemLevel%RvectorBlock(velocityfield)
    end if

    ! Initialize the linear form
    rform%itermCount = 1
    rform%Idescriptors(1) = DER_FUNC

    ! Prepare quick access arrays of temporal collection structure
    rcollectionTmp%DquickAccess(1) = dtime
    rcollectionTmp%DquickAccess(2) = dscale
    
    ! How many spatial dimensions are we?
    select case(rproblemLevel%rtriangulation%ndim)
    case(NDIM2D)
      ! Set pointers
      call storage_getbase_int(p_rboundaryCondition%h_IbdrCondCpIdx,&
          p_IbdrCondCpIdx)
      call storage_getbase_int(p_rboundaryCondition%h_IbdrCondType,&
          p_IbdrCondType)
      
      ! Loop over all boundary components
      do ibct = 1, p_rboundaryCondition%iboundarycount
        
        ! Loop over all boundary segments
        do isegment = p_IbdrCondCpIdx(ibct),&
                      p_IbdrCondCpIdx(ibct+1)-1
            
          ! Prepare quick access array of temporal collection structure
          rcollectionTmp%IquickAccess(1) = p_IbdrCondType(isegment)
          rcollectionTmp%IquickAccess(2) = isegment
          
          ! What type of boundary conditions are we?
          select case(p_IbdrCondType(isegment))
          case(BDR_DIRICHLET_WEAK,&
               BDR_INHOMNEUMANN_WEAK)
            ! Create boundary segment
            call bdrf_createRegion(p_rboundaryCondition, ibct,&
                isegment-p_IbdrCondCpIdx(ibct)+1, rboundaryRegion)
            
            ! Assemble the linear form
            call linf_buildVectorScalarBdr2d(rform, CUB_G3_1D,&
                .false., rvector, fcoeff_buildVectorScBdr2D_sim,&
                rboundaryRegion, rcollectionTmp)
          end select
          
        end do ! isegment
      end do ! ibct
        
    case default
      call output_line('Unsupported spatial dimension !',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_calcLinfBoundaryConditions')
      call sys_halt()
    end select
    
    ! Release temporal collection structure
    call collct_done(rcollectionTmp)
        
  end subroutine transp_calcLinfBoundaryConditions

  !*****************************************************************************

!<subroutine>
  
  subroutine transp_calcLinfBdrCondQuick(rproblemLevel, rsolver,&
      smode, ivelocitytype, dtime, dscale, rvectorScalar,&
      rcollection, fcb_coeffVecBdrPrimal_sim, fcb_coeffVecBdrDual_sim)
    
!<description>
    ! This subroutine is a shortcut for building the linear form
    ! arising from the weak imposition of boundary conditions. It
    ! calls the more general routine transp_calcLinfBoundaryConditions
    ! using the corresponding callback routines.
!</description>

!<input>
    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(in) :: rsolver

    ! problem mode (primal, dual)
    character(LEN=*), intent(in) :: smode

    ! type of velocity
    integer, intent(in) :: ivelocitytype
    
    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! user-defined callback functions
    include 'intf_transpCoeffVecBdr.inc'
    optional :: fcb_coeffVecBdrPrimal_sim
    optional :: fcb_coeffVecBdrDual_sim
!</intput>

!<inputoutput>
    ! scalar vector where to store the linear form
    type(t_vectorScalar), intent(inout) :: rvectorScalar

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! Are we in primal or dual mode?
    select case(trim(smode))
      
    case ('primal')
      
      ! @FAQ2: What type of velocity are we?
      select case(abs(ivelocitytype))
      case default
        ! The user-defined callback function is used if present;
        ! otherwise an error is throws
        if (present(fcb_coeffVecBdrPrimal_sim)) then
          call transp_calcLinfBoundaryConditions(rproblemLevel,&
              rsolver, dtime, dscale, fcb_coeffVecBdrPrimal_sim,&
              rvectorScalar, rcollection)
          
        else ! callback function not present
          
          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcLinfBdrCondQuick')
          call sys_halt()

        end if
        

      case (VELOCITY_ZERO)
        ! zero velocity, do nothing

        
      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP)
        ! linear velocity
        call transp_calcLinfBoundaryConditions(rproblemLevel, rsolver,&
            dtime, dscale, transp_coeffVecBdrPrimalConst2d,&
            rvectorScalar, rcollection)
        

      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers` equation in space-time
        call transp_calcLinfBoundaryConditions(rproblemLevel, rsolver,&
            dtime, dscale, transp_coeffVecBdrPrimalBurgersSpT2d,&
            rvectorScalar, rcollection)

        
      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time
        call transp_calcLinfBoundaryConditions(rproblemLevel, rsolver,&
            dtime, dscale, transp_coeffVecBdrPrimalBuckLevSpT2d,&
            rvectorScalar, rcollection)

        
      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers` equation in 1D
        
        ! @TODO: Implement weak boundary conditions
        print *, "Weak boundary conditions are not available!"
        stop
        
      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers` equation in 2D
        call transp_calcLinfBoundaryConditions(rproblemLevel, rsolver,&
            dtime, dscale, transp_coeffVecBdrPrimalBurgers2d,&
            rvectorScalar, rcollection)

        
      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D
        
        ! @TODO: Implement weak boundary conditions
        print *, "Weak boundary conditions are not available!"
        stop     
      end select
      
      
    case ('dual')
      
      ! @FAQ2: What type of velocity are we?
      select case(abs(ivelocitytype))
        case default
        ! The user-defined callback function is used if present;
        ! otherwise an error is throws
        if (present(fcb_coeffVecBdrDual_sim)) then
          call transp_calcLinfBoundaryConditions(rproblemLevel,&
              rsolver, dtime, dscale, fcb_coeffVecBdrDual_sim,&
              rvectorScalar, rcollection)
          
        else ! callback function not present
          
          call output_line('Missing user-defined callback function!',&
              OU_CLASS_ERROR,OU_MODE_STD,'transp_calcLinfBdrCondQuick')
          call sys_halt()

        end if
        

      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP)
        ! linear velocity
        call transp_calcLinfBoundaryConditions(rproblemLevel, rsolver,&
            dtime, -dscale, transp_coeffVecBdrDualConst2d,&
            rvectorScalar, rcollection)

        
      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers` equation in space-time
        
        ! @TODO: Implement weak boundary conditions
        print *, "Weak boundary conditions are not available!"
        stop

        
      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time
        
        ! @TODO: Implement weak boundary conditions
        print *, "Weak boundary conditions are not available!"
        stop

        
      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers` equation in 1D
        
        ! @TODO: Implement weak boundary conditions
        print *, "Weak boundary conditions are not available!"
        stop
        

      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers` equation in 2D
        
        ! @TODO: Implement weak boundary conditions
        print *, "Weak boundary conditions are not available!"
        stop

        
      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D
        
        ! @TODO: Implement weak boundary conditions
        print *, "Weak boundary conditions are not available!"
        stop
      end select

      
    case DEFAULT
      call output_line('Invalid mode!',&
          OU_CLASS_ERROR,OU_MODE_STD,'calcLinfBdrCond')
      call sys_halt()
    end select
    
  end subroutine transp_calcLinfBdrCondQuick
    
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

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_problemLevel), pointer :: p_rproblemLevel
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(NDIM3D+1) :: Dvalue
    character(LEN=SYS_STRLEN) :: svelocityname
    integer :: ieq, neq, idim, ndim, nlmin, icomp
    integer :: ivelocitytype, velocityfield, discretisation

    
    ! Check if the velocity "vector" needs to be generated explicitly
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'ivelocitytype', ivelocitytype)
    if ((abs(ivelocitytype) .ne. VELOCITY_CONSTANT) .and.&
        (abs(ivelocitytype) .ne. VELOCITY_TIMEDEP)) return

    ! Get parameter from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'velocityfield', velocityfield)
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'discretisation', discretisation)

    ! Get function parser from collection
    p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

    ! Set minimum problem level
    nlmin = rproblemLevel%ilev
    if (present(nlminOpt)) nlmin = nlminOpt

    ! Initialize variable values
    Dvalue           = 0.0_DP
    Dvalue(NDIM3D+1) = dtime

    ! Loop over all problem levels
    p_rproblemLevel => rproblemLevel
    do while(associated(p_rproblemLevel))

      ! Get number of degrees of freedom and spatial dimension
      p_rspatialDiscr => p_rproblemLevel%Rdiscretisation(discretisation)%RspatialDiscr(1)
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
      
      ! Get vertex coordinates of the current problem level
      call storage_getbase_double2d(&
          p_rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)


      ! Loop over all spatial dimensions
      do idim = 1, ndim

        ! Attach discretisation structure
        p_rproblemLevel%RvectorBlock(velocityfield)&
            %RvectorBlock(idim)%p_rspatialDiscr => p_rspatialDiscr

        ! Get scalar subvector
        call lsyssc_getbase_double(&
            p_rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(idim), p_Ddata)

        ! Retrieve function name from parameter list
        call parlst_getvalue_string(rparlist, ssectionName, 'svelocityname',&
                                    svelocityname, isubString=idim)

        ! Determine corresponding component number from the function parser
        icomp = fparser_getFunctionNumber(p_rfparser, svelocityname)

        ! Loop over all equations of scalar subvector
        do ieq = 1, neq
          Dvalue(1:ndim) = p_DvertexCoords(:,ieq)
          call fparser_evalFunction(p_rfparser, icomp, Dvalue, p_Ddata(ieq))
        end do
      end do

      ! Set update notification in problem level structure
      p_rproblemLevel%iproblemSpec = ior(p_rproblemLevel%iproblemSpec,&
                                         PROBLEV_MSPEC_UPDATE)
      
      ! Proceed to coarser problem level if minimum level has not been reached
      if (p_rproblemLevel%ilev .le. nlmin) exit
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      
    end do

  end subroutine transp_calcVelocityField

  !*****************************************************************************

!<subroutine>

  subroutine transp_setVelocityField(rvector)

!<description>
    ! This subroutine sets the global pointer to the velocity vector
    ! on the given problem level structure. Note that this subroutine
    ! will not work of multiple convection-diffusion-reaction problems
    ! are solved in parallel since there is only one global pointer.
!</description>

!<input>
    ! velocity field
    type(t_vectorBlock), intent(in) :: rvector
!</input>
!</subroutine>

    ! What spatial dimension are we?
    select case(rvector%nblocks)
    case (NDIM1D)
      call transp_setVariable1d(rvector%RvectorBlock(1), 1)

    case (NDIM2D)
      call transp_setVariable2d(rvector%RvectorBlock(1), 1)
      call transp_setVariable2d(rvector%RvectorBlock(2), 2)

    case (NDIM3D)
      call transp_setVariable3d(rvector%RvectorBlock(1), 1)
      call transp_setVariable3d(rvector%RvectorBlock(2), 2)
      call transp_setVariable3d(rvector%RvectorBlock(3), 3)

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_setVelocityField')
      call sys_halt()
    end select
    
  end subroutine transp_setVelocityField

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcLinearizedFCT(rbdrCond, rproblemLevel,&
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
    type(t_vectorScalar) :: rflux0, rflux
    type(t_vectorBlock) :: rdata
    real(DP), dimension(:), pointer :: p_MC, p_ML, p_Cx, p_Cy
    real(DP), dimension(:), pointer :: p_u, p_flux0, p_flux, p_data
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal, p_Ksep
    integer :: h_Ksep, templatematrix, lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, nedge

    ! Get parameters from parameter list which are required unconditionally
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'templatematrix', templateMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'coeffMatrix_CX', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'coeffMatrix_CY', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
        'lumpedmassmatrix', lumpedMassMatrix)

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
    call lsyssc_createVector(rflux0, nedge, .true., ST_DOUBLE)
    call lsyssc_createVector(rflux,  nedge, .true., ST_DOUBLE)
    call lsysbl_createVectorBlock(rsolution, rdata, .false.)
    
    ! Set pointers
    call lsysbl_getbase_double(rsolution, p_u)
    call lsysbl_getbase_double(rdata, p_data)
    call lsyssc_getbase_double(rflux, p_flux)
    call lsyssc_getbase_double(rflux0, p_flux0)

    ! Build the flux
    call buildFlux2d(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix&
        %NEQ, nedge, p_u, rtimestep%dStep, p_MC, p_ML, p_Cx, p_Cy,&
        p_data, p_flux, p_flux0)
   
    ! Build the correction and apply it directly
    call buildCorrection(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
        p_rmatrix%NEQ, nedge, p_ML, p_flux, p_flux0, p_data, p_u)
    
    ! Set boundary conditions explicitly
    call bdrf_filterVectorExplicit(rbdrCond, rsolution, rtimestep%dTime)

    ! Release flux vectors
    call storage_free(h_Ksep)
    call lsyssc_releaseVector(rflux0)
    call lsyssc_releaseVector(rflux)
    call lsysbl_releaseVector(rdata)

  contains

    !***************************************************************************

    subroutine buildFlux2d(Kld, Kcol, Kdiagonal, Ksep, NEQ, NEDGE, u,&
        dscale, MC, ML, Cx, Cy, troc, flux0, flux)

      real(DP), dimension(:), intent(in) :: MC,ML,Cx,Cy,u
      real(DP), intent(in) :: dscale
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NEDGE
      
      integer, dimension(:), intent(inout) :: Ksep
      real(DP), dimension(:), intent(inout) :: flux0,flux
      
      real(DP), dimension(:), intent(out) :: troc     

      ! local variables
      real(DP), dimension(NDIM2D) :: C_ii,C_ij, C_ji
      real(DP) :: k_ii,k_ij,k_ji,d_ij,aux,f_ij,f_ji
      integer :: ii,ij,ji,i,j,iedge

      ! Initialize time rate of change
      call lalg_clearVector(troc)

      ! Initialize edge counter
      iedge = 0
      
      ! Loop over all rows
      do i = 1, NEQ
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Compute coefficient
        C_ii(1) = Cx(ii);   C_ii(2) = Cy(ii)

        ! Compute convection coefficients
        call transp_calcMatrixPrimalConst2d(u(i), u(i),&
            C_ii, C_ii, i, i, k_ii, k_ii, d_ij)

        ! Update the time rate of change vector
        troc(i) = troc(i) + dscale*k_ii*u(i)

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

          ! Compute convection coefficients
          call transp_calcMatrixPrimalConst2d(u(i), u(j),&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Artificial diffusion coefficient
          d_ij = max(-k_ij, 0.0_DP, -k_ji)

          ! Compute auxiliary value
          aux = d_ij*(u(j)-u(i))
          
          ! Update the time rate of change vector
          troc(i) = troc(i) + dscale * (k_ij*u(j) + aux)
          troc(j) = troc(j) + dscale * (k_ji*u(i) - aux)

          ! Compute raw antidiffusive flux
          flux0(iedge) = -aux

        end do
      end do


      ! Scale the time rate of change by the lumped mass matrix
      do i = 1, NEQ
        troc(i) = troc(i)/ML(i)
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
          flux(iedge) = flux0(iedge) + MC(ij)*(troc(i)-troc(j))
          
          ! Update edge counter
          iedge = iedge-1
          
        end do
      end do

    end subroutine buildFlux2d

    !***************************************************************************
    
    subroutine buildCorrection(Kld, Kcol, Kdiagonal, Ksep, NEQ,&
        NEDGE, ML, flux, flux0, data, u)

      
      real(DP), dimension(:), intent(in) :: ML,flux0
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: NEQ,NEDGE
      
      real(DP), dimension(:), intent(inout) :: data,u,flux
      integer, dimension(:), intent(inout) :: Ksep
      
      ! local variables
      real(DP), dimension(:), allocatable :: pp,pm,qp,qm,rp,rm
      real(DP) :: f_ij,diff
      integer :: ij,ji,i,j,iedge,ivar

      ! Allocate temporal memory
      allocate(pp(neq), pm(neq), qp(neq), qm(neq), rp(neq), rm(neq))
      
      ! Initialize vectors
      call lalg_clearVector(pp)
      call lalg_clearVector(pm)
      call lalg_clearVector(qp)
      call lalg_clearVector(qm)
      call lalg_clearVector(rp, 1.0_DP)
      call lalg_clearVector(rm, 1.0_DP)
      
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

          ! Apply minmod prelimiter ...
          f_ij = minmod(flux(iedge), 2*flux0(iedge))
          
          ! ... and store prelimited flux
          flux(iedge) = f_ij
          diff        = u(j)-u(i)

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
        qp(i) = qp(i)*ML(i)
        qm(i) = qm(i)*ML(i)
        
        if (pp(i) > qp(i) + SYS_EPSREAL) rp(i) = qp(i)/pp(i)
        if (pm(i) < qm(i) - SYS_EPSREAL) rm(i) = qm(i)/pm(i)
      end do


      ! Initialize correction
      call lalg_clearVector(data)

      ! Loop over all rows (backward)
      do i = NEQ, 1, -1
        
        ! Loop over all off-diagonal matrix entries IJ which are adjacent to
        ! node J such that I < J. That is, explore the upper triangular matrix.
        do ij = Kld(i+1)-1, Ksep(i)+1, -1
          
          ! Get node number J, the corresponding matrix position JI,
          ! and let the separator point to the preceeding entry.
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)-1

          ! Limit conservative fluxes
          f_ij = flux(iedge)
          if (f_ij > 0.0_DP) then
            f_ij = min(rp(i), rm(j))*f_ij
          else
            f_ij = min(rm(i), rp(j))*f_ij
          end if
          
          ! Apply correction
          data(i) = data(i) + f_ij
          data(j) = data(j) - f_ij

          ! Update edge counter
          iedge = iedge-1
          
        end do
      end do


      do i = 1, NEQ
        u(i) = u(i) + data(i)/ML(i)
      end do

      ! Deallocate temporal memory
      deallocate(pp,pm,qp,qm,rp,rm)

    end subroutine buildCorrection

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
    
  end subroutine transp_calcLinearizedFCT

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVectorAnalytic(rdiscretisation, rform,&
      nelements, npointsPerElement, Dpoints, IdofsTest,&
      rdomainIntSubset, Dcoefficients, rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use fparser
    use scalarpde
    use triangulation
    
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

    ! An array accepting the DOF`s on all elements trial in the trial space.
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
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dtime
    integer :: itermCount, ipoint, iel, ndim, icomp
    
    
    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    p_rfparser => collct_getvalue_pars(rcollection,&
        trim(rcollection%SquickAccess(1)))

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
          call fparser_evalFunction(p_rfparser, icomp, 2, Dpoints(:,:&
              ,iel), Dcoefficients(itermCount,:,iel))
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
    
  end subroutine transp_coeffVectorAnalytic

  !*****************************************************************************

!<subroutine>

  subroutine transp_refFuncAnalytic(cderivative, rdiscretisation,&
      nelements, npointsPerElement, Dpoints, IdofsTest,&
      rdomainIntSubset, Dvalues, rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use triangulation
    
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

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
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
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(NDIM3D+1) :: Dvalue
    integer :: ipoint, iel, ndim, icomp


    ! Initialize values
    Dvalue = 0.0_DP
    
    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    p_rfparser => collct_getvalue_pars(rcollection,&
                                       trim(rcollection%SquickAccess(1)))
   
    ! Moreover, this subroutine assumes that the first quick access integer
    ! value holds the number of the function to be evaluated
    icomp = rcollection%IquickAccess(1)
    
    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    Dvalue(NDIM3D+1) = rcollection%DquickAccess(1)
    
    ! Set number of spatial dimensions
    ndim = size(Dpoints, 1)

    do iel = 1, nelements
      do ipoint = 1, npointsPerElement
        
        ! Set values for function parser
        Dvalue(1:ndim) = Dpoints(:, ipoint, iel)
        
        ! Evaluate function parser
        call fparser_evalFunction(p_rfparser, icomp, Dvalue, Dvalues(ipoint,iel))
      end do
    end do
    
  end subroutine transp_refFuncAnalytic

  !*****************************************************************************

!<subroutine>

  subroutine transp_weightFuncAnalytic(rdiscretisation, nelements,&
      npointsPerElement, Dpoints, IdofsTest, rdomainIntSubset,&
      Dvalues, rcollection)

    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use triangulation
    
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

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in trial space,Number of elements)
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
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(out) :: Dvalues
!</output>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(NDIM3D+1) :: Dvalue
    integer :: ipoint, iel, ndim, icomp


    ! Initialize values
    Dvalue = 0.0_DP
    
    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    p_rfparser => collct_getvalue_pars(rcollection,&
                                       trim(rcollection%SquickAccess(1)))
   
    ! Moreover, this subroutine assumes that the second quick access integer
    ! value holds the number of the function to be evaluated
    icomp = rcollection%IquickAccess(2)
    
    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    Dvalue(NDIM3D+1) = rcollection%DquickAccess(1)
    
    ! Set number of spatial dimensions
    ndim = size(Dpoints, 1)

    do iel = 1, nelements
      do ipoint = 1, npointsPerElement
        
        ! Set values for function parser
        Dvalue(1:ndim) = Dpoints(:, ipoint, iel)
        
        ! Evaluate function parser
        call fparser_evalFunction(p_rfparser, icomp, Dvalue, Dvalues(ipoint,iel))
      end do
    end do
    
  end subroutine transp_weightFuncAnalytic

end module transport_callback

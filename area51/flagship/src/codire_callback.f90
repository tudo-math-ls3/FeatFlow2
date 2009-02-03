!##############################################################################
!# ****************************************************************************
!# <name> codire_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve scalar conservation laws in arbitrary spatial dimensions.
!#
!# The following callback functions are available:
!#
!# 1.) codire_coeffFunctionParser
!#     -> Computes coefficients for the linear form by using a function parser
!#
!# 2.) codire_setBoundary
!#     -> Imposes boundary conditions for nonlinear solver
!#
!# 3.) codire_calcPreconditioner
!#     -> Calculates the nonlinear preconditioner
!#
!# 4.) codire_calcResidual
!#     -> Calculates the nonlinear residual vector
!#
!# 5.) codire_calcRHS
!#     -> Calculates the right-hand side vector
!#
!# </purpose>
!##############################################################################

module codire_callback

  use afcstabilisation
  use boundaryfilter
  use codire_basic
  use collection
  use dofmapping
  use fparser
  use fsystem
  use genoutput
  use groupfemscalar
  use linearsystemblock
  use linearsystemscalar
  use problem
  use solver
  use statistics

  implicit none

  private
  public :: codire_coeffFunctionParser
  public :: codire_coeffTargetFunc
  public :: codire_setBoundary
  public :: codire_calcPreconditioner
  public :: codire_calcJacobian
  public :: codire_applyJacobian
  public :: codire_calcResidual
  public :: codire_calcGalerkinResidual
  public :: codire_calcRHS
  public :: codire_calcVelocityField
  public :: codire_setVelocityField

!<globals>

  !*****************************************************************
  ! Pointers to the ACTIVE velocity field.
  !
  ! This global variable is not good programming style but it is the
  ! only way to allow for an efficient access to the velocity data
  ! from within the callback routines which are called repeatedly

  real(DP), dimension(:), pointer, save :: p_DvelocityX => null()
  real(DP), dimension(:), pointer, save :: p_DvelocityY => null()
  real(DP), dimension(:), pointer, save :: p_DvelocityZ => null()

!</globals>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine codire_coeffFunctionParser(rdiscretisation, rform, nelements,&
                                        npointsPerElement, Dpoints, IdofsTest,&
                                        rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use fparser
    use scalarpde
    use triangulation
    
!<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
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
    type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN) :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN) :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional :: rcollection
!</inputoutput>
  
!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT) :: Dcoefficients
!</output>
    
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: rfparser
    integer(I32) :: ielement
    
    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    rfparser => collct_getvalue_pars(rcollection,&
                                     trim(rcollection%SquickAccess(1)))

    ! Evaluate all coefficients using the function parser
    do ielement = 1, nelements
      call fparser_evalFunction(rfparser, 1, 2, Dpoints(:,:,ielement),&
                                Dcoefficients(1,:,ielement))
    end do
  end subroutine codire_coeffFunctionParser

  ! ***************************************************************************

!<subroutine>

  subroutine codire_coeffTargetFunc(rdiscretisation, rform, nelements,&
                                    npointsPerElement, Dpoints, IdofsTest,&
                                    rdomainIntSubset, Dcoefficients, rcollection)
    
    use basicgeometry
    use collection
    use domainintegration
    use fparser
    use scalarpde
    use triangulation
    
!<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
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
    type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN) :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN) :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN) :: rdomainIntSubset
!</input>

!<inputoutput>
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional :: rcollection
!</inputoutput>
  
!<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT) :: Dcoefficients
!</output>
    
!</subroutine>

    ! local variables
    integer :: ipoint, ielement
    real(DP) :: x,y

    do ielement = 1, nelements
      do ipoint = 1, npointsPerElement
        x = Dpoints(1, ipoint, ielement)
        y = Dpoints(2, ipoint, ielement)

        if (x .ge. 0.9 .and. x .le. 1.1) then
          Dcoefficients(:, ipoint, ielement) = 1.0_DP
        else
          Dcoefficients(:, ipoint, ielement) = 0.0_DP
        end if

!!$        if (sqrt(((x-1.625)**2 + (y-0.25)**2)) .le. 0.125) then
!!$          Dcoefficients(:, ipoint, ielement) = 0.06283185_DP
!!$        else
!!$          Dcoefficients(:, ipoint, ielement) = 0.0_DP
!!$        end if

!!$        ! Are we in the rectangle [1,2] x [0,0.1]?
!!$        if (x .gt. 1.0_DP .and. y .le. 0.1) then
!!$          Dcoefficients(:, ipoint, ielement) = 1.0_DP
!!$        else
!!$          Dcoefficients(:, ipoint, ielement) = 0.0_DP
!!$        end if
      end do
    end do
    
  end subroutine codire_coeffTargetFunc
  
  !*****************************************************************************

!<subroutine>

  subroutine codire_setBoundary(rproblemLevel, rtimestep, rsolver,&
                                rsol, rres, rsol0, rcollection)

!<description>
    ! This subroutine imposes the nonlinear boundary conditions.
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! nonlinear solver structure
    type(t_solver), intent(IN) :: rsolver

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: rsol0
!</input>

!<inputoutput>
    ! multigrid level
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsol

    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rres

    ! collection structure to provide additional
    ! information to the boundary setting routine
    type(t_collection), intent(InOUT) :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    integer :: imatrix

    select case(rsolver%iprecond)
    case (NLSOL_PRECOND_BLOCKD,&
          NLSOL_PRECOND_DEFCOR, &
          NLSOL_PRECOND_NEWTON_FAILED)

      ! Impose boundary conditions for the solution vector and impose
      ! zeros in the residual vector and the off-diagonal positions 
      ! of the system matrix which is obtained from the collection
      imatrix = collct_getvalue_int(rcollection, 'SystemMatrix')
      
      call bdrf_filterSolution(rsolver%rboundaryCondition,&
                               rproblemLevel%rtriangulation,&
                               rproblemLevel%Rmatrix(imatrix),&
                               rsol, rres, rsol0, rtimestep%dTime)

      
    case (NLSOL_PRECOND_NEWTON)

      ! Impose boundary conditions for the solution vector and impose
      ! zeros in the residual vector and the off-diagonal positions 
      ! of the system matrix which is obtained from the collection
      imatrix = collct_getvalue_int(rcollection, 'JacobianMatrix')

      call bdrf_filterSolution(rsolver%rboundaryCondition,&
                               rproblemLevel%rtriangulation,&
                               rproblemLevel%Rmatrix(imatrix),&
                               rsol, rres, rsol0, rtimestep%dTime)


    case DEFAULT
      call output_line('Invalid nonlinear preconditioner!',&
                       OU_CLASS_ERROR, OU_MODE_STD,'codire_setBoundary')
      call sys_halt()
    end select
  end subroutine codire_setBoundary

  !*****************************************************************************

!<subroutine>

  subroutine codire_calcPreconditioner(rproblemLevel, rtimestep, rsolver,&
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
    logical :: bStabilize,bbuildAFC,bbcompatible,bisDivergenceFree
    integer :: systemMatrix
    integer :: lumpedMassMatrix, consistentMassMatrix
    integer :: transportMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: coeffMatrix_S
    integer :: imasstype,iconvectiontype,idiffusiontype
    integer :: convectionAFC, diffusionAFC
    integer :: velocityfield
    integer :: primaldual

    print *, ">>> codire_calcPreconditioner"
    
    ! Start time measurement for matrix evaluation
    rtimer => collct_getvalue_timer(rcollection, 'timerAssemblyMatrix')
    call stat_startTimer(rtimer, STAT_TIMERSHORT)

    ! Check if the preconditioner has to be updated
    if (iand(rproblemLevel%iproblemSpec, PROBLEV_MSPEC_UPDATE) .eq. 0) return

    ! Get parameters from collection which are required unconditionally
    transportMatrix = collct_getvalue_int(rcollection, 'transportmatrix')
    coeffMatrix_CX  = collct_getvalue_int(rcollection, 'coeffMatrix_CX')
    coeffMatrix_CY  = collct_getvalue_int(rcollection, 'coeffMatrix_CY')
    coeffMatrix_CZ  = collct_getvalue_int(rcollection, 'coeffMatrix_CZ')
    coeffMatrix_S   = collct_getvalue_int(rcollection, 'coeffMatrix_S')

    
    !---------------------------------------------------------------------------
    ! Assemble diffusion operator
    !---------------------------------------------------------------------------
    
    idiffusiontype = collct_getvalue_int(rcollection, 'idiffusiontype')
    
    select case(idiffusiontype)
    case (DIFFUSION_ZERO)
      ! zero diffusion, clear the transport matrix
      call lsyssc_clearMatrix(rproblemLevel%Rmatrix(transportMatrix))
      
    case (DIFFUSION_ISOTROPIC)
      ! Isotropic diffusion
      call gfsc_buildDiffusionOperator(&
          rproblemLevel%Rmatrix(coeffMatrix_S), .false., .true.,&
          rproblemLevel%Rmatrix(transportMatrix))

    case (DIFFUSION_ANISOTROPIC)
      ! Anisotropic diffusion
      diffusionAFC = collct_getvalue_int(rcollection, 'diffusionAFC')

      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_SYMMETRIC)
        ! Check if stabilisation structure is initialised
        if (rproblemLevel%Rafcstab(diffusionAFC)%iSpec .eq. AFCSTAB_UNDEFINED)&
            call gfsc_initStabilisation(&
                rproblemLevel%Rmatrix(coeffMatrix_S),&
                rproblemLevel%Rafcstab(diffusionAFC))

        call gfsc_buildDiffusionOperator(&
            rproblemLevel%Rmatrix(coeffMatrix_S), .true., .true.,&
            rproblemLevel%Rmatrix(transportMatrix),&
            rproblemLevel%Rafcstab(diffusionAFC))

      case DEFAULT
        call gfsc_buildDiffusionOperator(&
            rproblemLevel%Rmatrix(coeffMatrix_S), .false., .true.,&
            rproblemLevel%Rmatrix(transportMatrix))
      end select

    case (DIFFUSION_VARIABLE)
      print *, "Variable diffusion matrices are yet not implemented!"
      stop
      
    case DEFAULT
      call output_line('Invalid type of diffusion!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcPreconditioner')
      call sys_halt()
    end select
    
    
    !---------------------------------------------------------------------------
    ! Assemble convective operator
    !---------------------------------------------------------------------------

    primaldual    = collct_getvalue_int(rcollection, 'primaldual')
    ivelocitytype = collct_getvalue_int(rcollection, 'ivelocitytype')
    convectionAFC = collct_getvalue_int(rcollection, 'convectionAFC')


    ! Check if velocity is assumed to be discretely divergence free
    bisDivergenceFree = (ivelocitytype .gt. 0)
    
    ! Check if stabilization should be applied
    bStabilize = AFCSTAB_GALERKIN .ne. &
                 rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation

    ! Check if stabilization structure should be built
    bbuildAFC = bStabilize .and. AFCSTAB_UPWIND .ne.&
                rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation

    ! Check if stabilisation structure is initialised
    if (rproblemLevel%Rafcstab(convectionAFC)%iSpec .eq. AFCSTAB_UNDEFINED)&
        call gfsc_initStabilisation(rproblemLevel%Rmatrix(transportMatrix),&
                                    rproblemLevel%Rafcstab(convectionAFC))


    ! Are we in primal or dual mode
    if (primaldual .eq. 1) then
      
      select case(abs(ivelocitytype))
      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        
      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP)
        ! linear velocity
        velocityfield = collct_getvalue_int(rcollection, 'velocityfield')
        call codire_setVelocityField(rproblemLevel%RvectorBlock(velocityfield))
        
        if (bbuildAFC) then
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                rsol, codire_calcPrimalConvConst1d, bisDivergenceFree,&
                .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
                rproblemLevel%Rafcstab(convectionAFC))

          case (NDIM2D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rsol, codire_calcPrimalConvConst2d, bisDivergenceFree,&
                .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
                rproblemLevel%Rafcstab(convectionAFC))

          case (NDIM3D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                rsol, codire_calcPrimalConvConst3d, bisDivergenceFree,&
                .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
                rproblemLevel%Rafcstab(convectionAFC))
          end select
        else
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                rsol, codire_calcPrimalConvConst1d, bisDivergenceFree,&
                bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix))

          case (NDIM2D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rsol, codire_calcPrimalConvConst2d, bisDivergenceFree,&
                bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix))

          case (NDIM3D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                rsol, codire_calcPrimalConvConst3d, bisDivergenceFree,&
                bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix))
          end select
        end if
        
      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers' equation in space-time
        if (bbuildAFC) then
          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsol, codire_calcConvectionBurgersSpT2d, bisDivergenceFree,&
              .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
              rproblemLevel%Rafcstab(convectionAFC))
        else
          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsol, codire_calcConvectionBurgersSpT2d, bisDivergenceFree,&
              bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix))
        end if

      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time
        if (bbuildAFC) then
          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsol, codire_calcConvectionBuckLevSpT2d, bisDivergenceFree,&
              .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
              rproblemLevel%Rafcstab(convectionAFC))
        else
          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsol, codire_calcConvectionBurgersSpT2d, bisDivergenceFree,&
              bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix))
        end if

      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers' equation in 1D
        if (bbuildAFC) then
          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rsol, codire_calcConvectionBurgers1d, bisDivergenceFree,&
              .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
              rproblemLevel%Rafcstab(convectionAFC))
        else
          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rsol, codire_calcConvectionBurgers1d, bisDivergenceFree,&
              bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix))
        end if
        
      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers' equation in 2D
        if (bbuildAFC) then
          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsol, codire_calcConvectionBurgers2d, bisDivergenceFree,&
              .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
              rproblemLevel%Rafcstab(convectionAFC))
        else
          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsol, codire_calcConvectionBurgers2d, bisDivergenceFree,&
              bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix))
        end if

      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D
        if (bbuildAFC) then
          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rsol, codire_calcConvectionBuckLev1d, bisDivergenceFree,&
              .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
              rproblemLevel%Rafcstab(convectionAFC))
        else
          call gfsc_buildConvectionOperator(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rsol, codire_calcConvectionBuckLev1d, bisDivergenceFree,&
              bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix))
        end if

      case DEFAULT
        call output_line('Invalid velocity profile!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_calcPreconditioner')
        call sys_halt()
      end select

    else

      select case(abs(ivelocitytype))
      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        
      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP)
        ! linear velocity
        velocityfield = collct_getvalue_int(rcollection, 'velocityfield')
        call codire_setVelocityField(rproblemLevel%RvectorBlock(velocityfield))
        
        if (bbuildAFC) then
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                rsol, codire_calcDualConvConst1d, bisDivergenceFree,&
                .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
                rproblemLevel%Rafcstab(convectionAFC))

          case (NDIM2D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rsol, codire_calcDualConvConst2d, bisDivergenceFree,&
                .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
                rproblemLevel%Rafcstab(convectionAFC))

          case (NDIM3D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                rsol, codire_calcDualConvConst3d, bisDivergenceFree,&
                .true., .false., rproblemLevel%Rmatrix(transportMatrix),&
                rproblemLevel%Rafcstab(convectionAFC))
          end select
        else
          select case(rproblemLevel%rtriangulation%ndim)
          case (NDIM1D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                rsol, codire_calcDualConvConst1d, bisDivergenceFree,&
                bstabilize, .false., rproblemLevel%Rmatrix(transportMatrix))

          case (NDIM2D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                rsol, codire_calcDualConvConst2d, bisDivergenceFree,&
                bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix))

          case (NDIM3D)
            call gfsc_buildConvectionOperator(&
                rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
                rsol, codire_calcDualConvConst3d, bisDivergenceFree,&
                bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix))
          end select
        end if
        
!!$      case (VELOCITY_BURGERS_SPACETIME)
!!$        ! nonlinear Burgers' equation in space-time
!!$        
!!$      case (VELOCITY_BUCKLEV_SPACETIME)
!!$        ! nonlinear Buckley-Leverett equation in space-time
!!$        
!!$      case (VELOCITY_BURGERS1D)
!!$        ! nonlinear Burgers' equation in 1D
!!$        
!!$      case (VELOCITY_BURGERS2D)
!!$        ! nonlinear Burgers' equation in 2D
!!$        
!!$      case (VELOCITY_BUCKLEV1D)
!!$        ! nonlinear Buckley-Leverett equation in 1D
        
      case DEFAULT
        call output_line('Invalid velocity profile!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_calcPreconditioner')
        call sys_halt()
      end select

    end if

    ! Remove update notifier for constant velocity
    if (abs(ivelocitytype) .eq. VELOCITY_CONSTANT) then
      rproblemLevel%iproblemSpec = iand(rproblemLevel%iproblemSpec,&
                                        not(PROBLEV_MSPEC_UPDATE))
    end if
    
    
    !---------------------------------------------------------------------------
    ! Assemble the global system operator
    !---------------------------------------------------------------------------
    
    systemMatrix  = collct_getvalue_int(rcollection, 'systemmatrix')
    imasstype     = collct_getvalue_int(rcollection, 'imasstype')

    select case(imasstype)
    case (MASS_LUMPED)
      ! Compute the global operator for transient flow
      !
      !   $ A = ML-theta*dt*L $

      lumpedMassMatrix  = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(lumpedMassMatrix), 1._DP,&
                                   rproblemLevel%Rmatrix(transportMatrix),&
                                   -rtimestep%theta*rtimestep%dStep,&
                                   rproblemLevel%Rmatrix(systemMatrix),&
                                   .false., .false., .true., .true.)
    case (MASS_CONSISTENT)
      ! Compute the global operator for transient flow
      !
      !   $ A = MC-theta*dt*L $

      consistentMassMatrix  = collct_getvalue_int(rcollection, 'consistentmassmatrix')
      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(consistentMassMatrix), 1._DP,&
                                   rproblemLevel%Rmatrix(transportMatrix),&
                                   -rtimestep%theta*rtimestep%dStep,&
                                   rproblemLevel%Rmatrix(systemMatrix),&
                                   .false., .false., .true., .true.)

    case DEFAULT
      ! Compute the global operator for steady-state flow
      !
      !   $ A = -L $
      !
      call lsyssc_copyMatrix(rproblemLevel%Rmatrix(transportMatrix),&
                             rproblemLevel%Rmatrix(systemMatrix))
      call lsyssc_scaleMatrix(rproblemLevel%Rmatrix(systemMatrix), -1._DP)
    
    end select
    

    !---------------------------------------------------------------------------
    ! Impose boundary conditions
    !---------------------------------------------------------------------------
      
    call bdrf_filterMatrix(rsolver%rboundaryCondition,&
                           rproblemLevel%rtriangulation,&
                           rproblemLevel%Rmatrix(systemMatrix), 1.0_DP)
        
    ! Ok, we updated the (nonlinear) system operator successfully. Now we still 
    ! have to link it to the solver hierarchy. This is done recursively.
    call codire_updateSolverMatrix(rproblemLevel, rsolver, systemMatrix, UPDMAT_ALL,&
                                   rproblemLevel%ilev, rproblemLevel%ilev)

    ! Finally, we have to update the content of the solver hierarchy
    call solver_updateContent(rsolver)

    ! Stop time measurement for global operator
    call stat_stopTimer(rtimer)
    
    print *, "<<< codire_calcPreconditioner"

  end subroutine codire_calcPreconditioner

  !*****************************************************************************

!<subroutine>

  subroutine codire_calcJacobian(rproblemLevel, rtimestep, rsolver,&
                                 rsol, rsol0, bfailure, rcollection)

!<description>
    ! This callback subroutine computes the Jacobian matrix.
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
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsol

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_timer), pointer :: rtimer
    real(DP) :: hstep
    integer :: templateMatrix
    integer :: transportMatrix, jacobianMatrix
    integer :: consistentMassMatrix, lumpedMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ
    integer :: coeffMatrix_S
    integer :: convectionAFC, diffusionAFC
    integer :: imasstype, imassantidiffusion, ivelocitytype, idiffusiontype
    integer :: primaldual, velocityfield
    logical :: bStabilize, bisExactStructure


    print *, ">>> codire_calcJacobian"

    ! Start time measurement for matrix evaluation
    rtimer => collct_getvalue_timer(rcollection, 'timerAssemblyMatrix')
    call stat_startTimer(rtimer, STAT_TIMERSHORT)

    ! Get parameters from collection which are required unconditionally
    consistentMassMatrix  = collct_getvalue_int(rcollection, 'consistentmassmatrix')
    lumpedMassMatrix      = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
    transportMatrix       = collct_getvalue_int(rcollection, 'transportmatrix')
    jacobianMatrix        = collct_getvalue_int(rcollection, 'jacobianmatrix')
    coeffMatrix_CX        = collct_getvalue_int(rcollection, 'coeffMatrix_CX')
    coeffMatrix_CY        = collct_getvalue_int(rcollection, 'coeffMatrix_CY')
    coeffMatrix_CZ        = collct_getvalue_int(rcollection, 'coeffMatrix_CZ')
    coeffMatrix_S         = collct_getvalue_int(rcollection, 'coeffMatrix_S')

    !---------------------------------------------------------------------------
    ! If the Newton iteration failed completely, adopt the low-order
    ! preconditioner and return from this subroutine
    ! ---------------------------------------------------------------------------
    if (bfailure) then
      call codire_calcPreconditioner(rproblemLevel, rtimestep, rsolver, rsol, rcollection)
      print *, "<<< codire_calcJacobian"
      return
    end if
    
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
      hstep = ( (1+&
          lsysbl_vectorNorm(rsol, LINALG_NORMEUCLID))*SYS_EPSREAL )**(1._DP/3._DP)
      
    case (PERTURB_SQRTEPS)
      hstep= sqrt(SYS_EPSREAL)

    case DEFAULT
      hstep = max(SYS_EPSREAL, rsolver%p_solverNewton%dperturbationStrategy)
    end select
    

    !---------------------------------------------------------------------------
    ! Assemble diffusion operator
    !---------------------------------------------------------------------------

    idiffusiontype = collct_getvalue_int(rcollection, 'idiffusiontype')
    diffusionAFC   = collct_getvalue_int(rcollection, 'diffusionAFC')

    select case(idiffusiontype)
    case (DIFFUSION_ZERO)
      ! zero diffusion, clear the system matrix
      call lsyssc_clearMatrix(rproblemLevel%Rmatrix(transportMatrix))
      
    case (DIFFUSION_ISOTROPIC,&
          DIFFUSION_ANISOTROPIC)
      ! Isotropic diffusion
      call gfsc_buildDiffusionOperator(&
          rproblemLevel%Rmatrix(coeffMatrix_S), .false., .true.,&
          rproblemLevel%Rmatrix(transportMatrix))
      
    case (DIFFUSION_VARIABLE)
      print *, "Variable diffusion matrices are yet not implemented!"
      stop
      
    case DEFAULT
      call output_line('Invalid type of diffusion!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
      call sys_halt()
    end select
    
    
    !---------------------------------------------------------------------------
    ! Assemble convection operator
    !---------------------------------------------------------------------------

    primaldual    = collct_getvalue_int(rcollection, 'primaldual')
    ivelocitytype = collct_getvalue_int(rcollection, 'ivelocitytype')
    convectionAFC = collct_getvalue_int(rcollection, 'convectionAFC')

    ! Check if stabilization should be applied
    bStabilize = (AFCSTAB_GALERKIN .ne.&
                  rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
    
    ! Are we in primal or dual mode
    if (primaldual .eq. 1) then

      select case (abs(ivelocitytype))
      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        
      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP) 
        ! linear velocity
        velocityfield = collct_getvalue_int(rcollection, 'velocityfield')
        call codire_setVelocityField(rproblemLevel%RvectorBlock(velocityfield))

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsc_buildConvectionJacobian(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rsol, codire_calcPrimalConvConst1d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))

        case (NDIM2D)
          call gfsc_buildConvectionJacobian(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsol, codire_calcPrimalConvConst2d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))

        case (NDIM3D)
          call gfsc_buildConvectionJacobian(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rsol, codire_calcPrimalConvConst3d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))
        end select
      
      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers' equation in space-time
        call gfsc_buildConvectionJacobian(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rsol, codire_calcConvectionBurgersSpT2d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))
        
      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time
        call gfsc_buildConvectionJacobian(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsol,&
            codire_calcConvectionBuckLevSpT2d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))
        
      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers' equation in 1D
        call gfsc_buildConvectionJacobian(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsol,&
            codire_calcConvectionBurgers1d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))
        
      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers' equation in 2D
        call gfsc_buildConvectionJacobian(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsol,&
            codire_calcConvectionBurgers2d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))
        
      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D
        call gfsc_buildConvectionJacobian(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsol,&
            codire_calcConvectionBuckLev1d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))
        
      case DEFAULT
        call output_line('Unsupported velocity type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
        call sys_halt()
      end select

    else
      
      select case (abs(ivelocitytype))
      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        
      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP) 
        ! linear velocity
        velocityfield = collct_getvalue_int(rcollection, 'velocityfield')
        call codire_setVelocityField(rproblemLevel%RvectorBlock(velocityfield))

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsc_buildConvectionJacobian(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rsol, codire_calcDualConvConst1d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))

        case (NDIM2D)
          call gfsc_buildConvectionJacobian(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsol, codire_calcDualConvConst2d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))

        case (NDIM3D)
          call gfsc_buildConvectionJacobian(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rsol, codire_calcDualConvConst3d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))
        end select
        
!!$      case (VELOCITY_BURGERS_SPACETIME)
!!$        ! nonlinear Burgers' equation in space-time
!!$        
!!$      case (VELOCITY_BUCKLEV_SPACETIME)
!!$        ! nonlinear Buckley-Leverett equation in space-time
!!$        
!!$      case (VELOCITY_BURGERS1D)
!!$        ! nonlinear Burgers' equation in 1D
!!$        
!!$      case (VELOCITY_BURGERS2D)
!!$        ! nonlinear Burgers' equation in 2D
!!$        
!!$      case (VELOCITY_BUCKLEV1D)
!!$        ! nonlinear Buckley-Leverett equation in 1D

      case DEFAULT
        call output_line('Unsupported velocity type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
        call sys_halt()
      end select

    end if
    
    
    ! Check if the Jacobian operator needs to be generated
    if (.not.lsyssc_hasMatrixStructure(rproblemLevel%Rmatrix(jacobianMatrix))) then
      
      templateMatrix = collct_getvalue_int(rcollection, 'templateMatrix')

      if ((rproblemLevel%Rafcstab(diffusionAFC)%iextendedJacobian .eq. 0) .and.&
          (rproblemLevel%Rafcstab(convectionAFC)%iextendedJacobian .eq. 0)) then
        
        bisExactStructure = .true.

        ! Adopt the standard sparsity pattern from the template matrix
        call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
                                    rproblemLevel%Rmatrix(jacobianMatrix),&
                                    LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      else

        bisExactStructure = .false.

        ! Extend the standard sparsity pattern of the template matrix
        call afcstab_generateExtSparsity(rproblemLevel%Rmatrix(templateMatrix),&
                                         rproblemLevel%Rmatrix(jacobianMatrix))
      end if
    end if
    
    
    !---------------------------------------------------------------------------
    ! Assemble the global system operator for the high-/low-order contribution
    !---------------------------------------------------------------------------

    imasstype = collct_getvalue_int(rcollection, 'imasstype')

    select case(imasstype)
    case (MASS_LUMPED)
      ! Compute the global Jacobian for transient flow
      !
      !   $ J = ML-theta*dt*L $
      
      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(transportMatrix),&
                                   -rtimestep%theta*rtimestep%dStep,&
                                   rproblemLevel%Rmatrix(lumpedMassMatrix), 1._DP,&
                                   rproblemLevel%Rmatrix(jacobianMatrix),&
                                   .false., .false., .true., bisExactStructure)
    case (MASS_CONSISTENT)
      ! Compute the global Jacobian for transient flow
      !
      !   $ J = MC-theta*dt*L $
      
      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(transportMatrix),&
                                   -rtimestep%theta*rtimestep%dStep,&
                                   rproblemLevel%Rmatrix(consistentMassMatrix), 1._DP,&
                                   rproblemLevel%Rmatrix(jacobianMatrix),&
                                   .false., .false., .true., bisExactStructure)
    case DEFAULT
      ! Compute the global Jacobian for steady-state flow
      !
      !   $ J = -L $
      
      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(transportMatrix), -1._DP,&
                                   rproblemLevel%Rmatrix(jacobianMatrix), 0._DP,&
                                   rproblemLevel%Rmatrix(jacobianMatrix),&
                                   .false., .false., .true., bisExactStructure)
    end select


    !---------------------------------------------------------------------------
    ! Assemble the Jacobian matrix for the diffusion operator
    !---------------------------------------------------------------------------
    
    ! What kind of diffusion are we?
    select case(idiffusiontype)
    case (DIFFUSION_ZERO, DIFFUSION_ISOTROPIC)
      ! zero diffusion or isotropic diffusion, do nothing
      
    case (DIFFUSION_ANISOTROPIC)
      ! Anisotropic diffusion
      select case(rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_SYMMETRIC)
        call gfsc_buildJacobianSymm(rsol, 1._DP, hstep, .false.,&
                                    rproblemLevel%Rafcstab(diffusionAFC),&
                                    rproblemLevel%Rmatrix(jacobianMatrix))
      end select

    case (DIFFUSION_VARIABLE)
      print *, "Variable diffusion matrices are yet not implemented!"
      stop
      
    case DEFAULT
      call output_line('Invalid type of diffusion!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Assemble the Jacobian matrix for the convection operator
    !---------------------------------------------------------------------------
    
    ! What kind of velocity are we?
    select case(abs(ivelocitytype))
    case (VELOCITY_ZERO)
      ! zero velocity, do nothing
      
    case(VELOCITY_CONSTANT,&
         VELOCITY_TIMEDEP) 
      ! linear velocity

      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)
        
        imassantidiffusion = collct_getvalue_int(rcollection, 'imassantidiffusion')

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusion .eq. MASS_CONSISTENT) then
          call gfsc_buildJacobianFCT(rsol, rtimestep%theta, rtimestep%dStep, hstep,&
                                     .false., rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildJacobianFCT(rsol, rtimestep%theta, rtimestep%dStep, hstep,&
                                     .false., rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildJacobianTVD(rsol, rtimestep%dStep, hstep, .false.,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(jacobianMatrix))

      case (AFCSTAB_FEMGP)
        call gfsc_buildJacobianGP(rproblemLevel%Rmatrix(consistentMassMatrix), rsol,&
                                  rsol0, rtimestep%theta, rtimestep%dStep, hstep,&
                                  .false., rproblemLevel%Rafcstab(convectionAFC),&
                                  rproblemLevel%Rmatrix(jacobianMatrix))
      end select

    case(VELOCITY_BURGERS_SPACETIME)
      ! nonlinear Burgers' equation in space-time

      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)

        imassantidiffusion = collct_getvalue_int(rcollection, 'imassantidiffusion')

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusion .eq. MASS_CONSISTENT) then
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsol, codire_calcConvectionBurgersSpT2d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsol, codire_calcConvectionBurgersSpT2d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsol, codire_calcConvectionBurgersSpT2d,&
                                   rtimestep%dStep, hstep, .false.,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(jacobianMatrix))

      case (AFCSTAB_FEMGP)
        call gfsc_buildJacobianGP(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                                  rsol, rsol0, codire_calcConvectionBurgersSpT2d,&
                                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                  rproblemLevel%Rafcstab(convectionAFC),&
                                  rproblemLevel%Rmatrix(jacobianMatrix))
      end select

    case(VELOCITY_BUCKLEV_SPACETIME)
      ! nonlinear Buckley-Leverett equation in space-time

      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)

        imassantidiffusion = collct_getvalue_int(rcollection, 'imassantidiffusion')

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusion .eq. MASS_CONSISTENT) then
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsol, codire_calcConvectionBuckLevSpT2d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsol, codire_calcConvectionBuckLevSpT2d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsol, codire_calcConvectionBuckLevSpT2d,&
                                   rtimestep%dStep, hstep, .false.,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(jacobianMatrix))
        
      case (AFCSTAB_FEMGP)
        call gfsc_buildJacobianGP(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                                  rsol, rsol0, codire_calcConvectionBuckLevSpT2d,&
                                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                  rproblemLevel%Rafcstab(convectionAFC),&
                                  rproblemLevel%Rmatrix(jacobianMatrix))
      end select
      
    case(VELOCITY_BURGERS1D)
      ! nonlinear Burgers' equation in 1D
      
      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)

        imassantidiffusion = collct_getvalue_int(rcollection, 'imassantidiffusion')

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusion .eq. MASS_CONSISTENT) then
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsol, codire_calcConvectionBurgers1d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsol, codire_calcConvectionBurgers1d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                   rsol, codire_calcConvectionBurgers1d,&
                                   rtimestep%dStep, hstep, .false.,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(jacobianMatrix))

      case (AFCSTAB_FEMGP)
        call gfsc_buildJacobianGP(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                                  rsol, rsol0, codire_calcConvectionBurgers1d,&
                                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                  rproblemLevel%Rafcstab(convectionAFC),&
                                  rproblemLevel%Rmatrix(jacobianMatrix))
      end select

    case(VELOCITY_BURGERS2D)
      ! nonlinear Burgers' equation in 2D
      
      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)

        imassantidiffusion = collct_getvalue_int(rcollection, 'imassantidiffusion')

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusion .eq. MASS_CONSISTENT) then
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsol, codire_calcConvectionBurgers2d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsol, codire_calcConvectionBurgers2d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsol, codire_calcConvectionBurgers2d,&
                                   rtimestep%dStep, hstep, .false.,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(jacobianMatrix))

      case (AFCSTAB_FEMGP)
        call gfsc_buildJacobianGP(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                                  rsol, rsol0, codire_calcConvectionBurgers2d,&
                                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                  rproblemLevel%Rafcstab(convectionAFC),&
                                  rproblemLevel%Rmatrix(jacobianMatrix))
      end select
      
    case(VELOCITY_BUCKLEV1D)
      ! nonlinear Buckley-Leverett equation in 1D

      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)

        imassantidiffusion = collct_getvalue_int(rcollection, 'imassantidiffusion')

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusion .eq. MASS_CONSISTENT) then
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsol, codire_calcConvectionBuckLev1d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsol, codire_calcConvectionBuckLev1d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                   rsol, codire_calcConvectionBuckLev1d,&
                                   rtimestep%dStep, hstep, .false.,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(jacobianMatrix))

      case (AFCSTAB_FEMGP)
        call gfsc_buildJacobianGP(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                                  rsol, rsol0, codire_calcConvectionBuckLev1d,&
                                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                  rproblemLevel%Rafcstab(convectionAFC),&
                                  rproblemLevel%Rmatrix(jacobianMatrix))
      end select
      
    case DEFAULT
      call output_line('Unsupported velocity type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcJacobian')
      call sys_halt()
    end select
    
    
    !---------------------------------------------------------------------------
    ! Impose boundary conditions
    !---------------------------------------------------------------------------
    
    call bdrf_filterMatrix(rsolver%rboundaryCondition, rproblemLevel%rtriangulation,&
                           rproblemLevel%Rmatrix(jacobianMatrix), 1.0_DP)

    ! Ok, we updated the Jacobian matrix successfully. Now we still have to
    ! link it to the solver hierarchy. This is done recursively.

    ! What type of flow are we?
    select case(imasstype)
    case (MASS_LUMPED,&
          MASS_CONSISTENT)
      call codire_updateSolverMatrix(rproblemLevel, rsolver,&
                                     jacobianMatrix, UPDMAT_JAC_TRANSIENT,&
                                     rproblemLevel%ilev, rproblemLevel%ilev)

    case DEFAULT
      call codire_updateSolverMatrix(rproblemLevel, rsolver,&
                                     jacobianMatrix, UPDMAT_JAC_STEADY,&
                                     rproblemLevel%ilev, rproblemLevel%ilev)
    end select
    
    ! Finally, we have to update the content of the solver hierarchy
    call solver_updateContent(rsolver)
    
    ! Stop time measurement for matrix evaluation
    call stat_stopTimer(rtimer)

    print *, "<<< codire_calcJacobian"

  end subroutine codire_calcJacobian

  !*****************************************************************************

!<subroutine>

  subroutine codire_applyJacobian(rproblemLevel, rx, ry, cx, cy, rcollection)

!<description>
    ! This subroutine applies the (scaled) Jacobian matrix to 
    ! the vector rx and adds the result to the vector ry.
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

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: jacobianMatrix

    jacobianMatrix = collct_getvalue_int(rcollection, 'jacobianMatrix')

    ! We know where the Jacobian matrix is stored and can apply it
    ! by means of matrix vector multiplication
    call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(jacobianMatrix),&
                             rx%RvectorBlock(1), ry%RvectorBlock(1), cx, cy)

  end subroutine codire_applyJacobian

  !*****************************************************************************

!<subroutine>

  subroutine codire_calcResidual(rproblemLevel, rtimestep, rsolver,&
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
    integer :: transportMatrix
    integer :: consistentMassMatrix, lumpedMassMatrix
    integer :: convectionAFC, diffusionAFC
    integer :: imasstype,imassantidiffusion

    
    print *, ">>> codire_calcResidual"

    ! For nonlinear conservation laws, the global system operator
    ! needs to be updated in each nonlinear iteration. The update
    ! routine is written such that it first determines if the problem
    ! is nonlinear and returns without matrix update otherwise.
    call codire_calcPreconditioner(rproblemLevel, rtimestep, rsolver, rsol, rcollection)


    ! Start time measurement for residual/rhs evaluation
    rtimer => collct_getvalue_timer(rcollection, 'timerAssemblyVector')
    call stat_startTimer(rtimer, STAT_TIMERSHORT)

    ! Get parameters from collection which are required unconditionally
    transportMatrix      = collct_getvalue_int(rcollection, 'transportmatrix')
    lumpedMassMatrix     = collct_getvalue_int(rcollection, 'lumpedmassmatrix')
    consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')

    ! Are we in the zeroth iteration?
    if (ite .eq. 0) then

      !-------------------------------------------------------------------------
      ! In the first nonlinear iteration update we compute
      !
      ! - the residual $ r=dt*L*u^n+f $ and
      !
      ! - the right hand side $ b=[M_L+(1-theta)*dt*L]*u^n $
      !
      !-------------------------------------------------------------------------
      
      imasstype = collct_getvalue_int(rcollection, 'imasstype')
      
      select case(imasstype)
      case (MASS_LUMPED)
        
        ! Compute the initial low-order residual
        !
        !   $  res = dt*L(u^n)*u^n $

        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                                 rsol%rvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 rtimestep%dStep, 0.0_DP)
        
        ! Compute the constant right-hand side, whereby the force vector
        ! is already given in the right-hand side vector
        !
        !   $ rhs = M_L*u^n+(1-theta)*res $

        call lsysbl_copyVector(rres, rrhs)
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                 rsol%RvectorBlock(1),&
                                 rrhs%RvectorBlock(1),&
                                 1.0_DP, 1.0_DP-rtimestep%theta)
      case (MASS_CONSISTENT)
        
        ! Compute the initial low-order residual
        !
        !   $  res = dt*L(u^n)*u^n $

        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                                 rsol%rvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 rtimestep%dStep, 0.0_DP)
        
        ! Compute the constant right-hand side, whereby the force vector
        ! is already given in the right-hand side vector
        !
        !   $ rhs = M_C*u^n+(1-theta)*res $

        call lsysbl_copyVector(rres, rrhs)
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                 rsol%RvectorBlock(1),&
                                 rrhs%RvectorBlock(1),&
                                 1.0_DP, 1.0_DP-rtimestep%theta)
      case DEFAULT
        
        ! Compute the initial low-order residual
        !
        !   $ res = L(u^n)*u^n $
        !
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                                 rsol%rvectorBlock(1),&
                                 rres%RvectorBlock(1), 1.0_DP, 0.0_DP)      
      end select
     

      ! Perform algebraic flux correction for the convective term if required
      !
      !   $ res = res + f^*(u^n+1,u^n) $

      convectionAFC = collct_getvalue_int(rcollection, 'convectionAFC')
            
      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
          
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)

        imassantidiffusion = collct_getvalue_int(rcollection, 'imassantidiffusion')

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusion .eq. MASS_CONSISTENT) then
          call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                     rsol, rtimestep%theta, rtimestep%dStep, .true.,&
                                     rres, rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                     rsol, rtimestep%theta, rtimestep%dStep, .true.,&
                                     rres, rproblemLevel%Rafcstab(convectionAFC))
        end if

      case (AFCSTAB_FEMTVD)
        call gfsc_buildResidualTVD(rsol, rtimestep%dStep, rres,&
                                   rproblemLevel%Rafcstab(convectionAFC))

      case (AFCSTAB_FEMGP)
        call gfsc_buildResidualGP(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                  rsol, rsol0, rtimestep%theta, rtimestep%dStep,&
                                  rres, rproblemLevel%Rafcstab(convectionAFC))
      end select

      
      ! Perform algebraic flux correction for the diffusive term if required
      !
      !   $ res = res + g^*(u^n+1,u^n) $

      diffusionAFC = collct_getvalue_int(rcollection, 'diffusionAFC')

      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation)

      case (AFCSTAB_SYMMETRIC)
        call gfsc_buildResidualSymm(rsol, 1.0_DP, rres, rproblemLevel%Rafcstab(diffusionAFC))
      end select
      
    else   ! ite > 0

      !-------------------------------------------------------------------------
      ! In all subsequent nonlinear iterations only the residual vector is
      ! updated, using the right-hand side vector from the first iteration
      !-------------------------------------------------------------------------

      imasstype = collct_getvalue_int(rcollection, 'imasstype')

      select case(imasstype)
      case (MASS_LUMPED)
        
        ! Compute the low-order residual
        !
        !   $ res = rhs+dt*theta*L(u^(m))*u^(m)-M_L*u $

        call lsysbl_copyVector(rrhs, rrhs)
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                                 rsol%rvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 rtimestep%theta*rtimestep%dStep, 1.0_DP)
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                 rsol%RvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 -1.0_DP, 1.0_DP)
      case (MASS_CONSISTENT)
        
        ! Compute the low-order residual
        !
        !   $ res = rhs+dt*theta*L(u^(m))*u^(m)-M_C*u $

        call lsysbl_copyVector(rrhs, rrhs)
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                                 rsol%rvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 rtimestep%theta*rtimestep%dStep, 1.0_DP)
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                 rsol%RvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 -1.0_DP, 1.0_DP)
      case DEFAULT
        
        ! Compute the low-order residual
        !
        !   $ res = L(u^(m))*u^(m) $
        !
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                                 rsol%rvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 1.0_DP, 0.0_DP)
      end select


      ! Perform algebraic flux correction for the convective term if required
      !
      !   $ res = res + f^*(u^(m),u^n) $

      convectionAFC = collct_getvalue_int(rcollection, 'convectionAFC')

      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
          
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)
        
        imassantidiffusion = collct_getvalue_int(rcollection, 'imassantidiffusion')

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusion .eq. MASS_CONSISTENT) then
          call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                     rsol, rtimestep%theta, rtimestep%dStep, .false.,&
                                     rres, rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                     rsol, rtimestep%theta, rtimestep%dStep, .false.,&
                                     rres, rproblemLevel%Rafcstab(convectionAFC))
        end if

      case (AFCSTAB_FEMTVD)
        call gfsc_buildResidualTVD(rsol, rtimestep%dStep, rres,&
                                   rproblemLevel%Rafcstab(convectionAFC))

      case (AFCSTAB_FEMGP)
        call gfsc_buildResidualGP(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                  rsol, rsol0, rtimestep%theta, rtimestep%dStep,&
                                  rres, rproblemLevel%Rafcstab(convectionAFC))
      end select


      ! Perform algebraic flux correction for the diffusive term if required
      !
      !   $ res = res + g^*(u^n+1,u^n) $

      diffusionAFC = collct_getvalue_int(rcollection, 'convectionAFC')
      
      ! What kind of stabilisation should be applied?
      select case(rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation)

      case (AFCSTAB_SYMMETRIC)
        call gfsc_buildResidualSymm(rsol, 1.0_DP, rres, rproblemLevel%Rafcstab(diffusionAFC))
      end select
      
    end if

    ! Stop time measurement for residual/rhs evaluation
    call stat_stopTimer(rtimer)
    
    print *, "<<< codire_calcResidual"

  end subroutine codire_calcResidual

  !*****************************************************************************

!<subroutine>

  subroutine codire_calcGalerkinResidual(rproblemLevel, rsolution,&
                                         rresidual, rcollection, rf)

!<description>
    ! This subroutine calculates the Galerkin residual for a
    ! (converged) steady-state solution to the primal problem
!</description>

!<input>
    ! solution vector
    type(t_vectorBlock), intent(IN) :: rsolution

    ! OPTIONAL: right-hand side vector
    type(t_vectorBlock), intent(IN), optional :: rf
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! residual vector
    type(t_vectorBlock), intent(INOUT) :: rresidual

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixScalar) :: rmatrix
    logical :: bcompatible,bisDivFree

    integer :: imode = 1

    ! Check if vectors are compatible
    call lsysbl_isVectorCompatible(rsolution, rresidual, bcompatible)
    
    ! Re-create residual vector if not compatible
    if (.not.bcompatible) then
      call lsysbl_releaseVector(rresidual)
      call lsysbl_duplicateVector(rsolution, rresidual,&
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
    end if

    ! Create empty matrix for the global system operator
    call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(CDEQ_MATRIX_TEMPLATE),&
                                rmatrix, LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)

    !------------------------------------------------------------
    ! Assemble diffusive operator
    !------------------------------------------------------------
    
    select case(idiffusiontype)
    case (DIFFUSION_ZERO)
      ! zero diffusion, clear the system matrix
      call lsyssc_clearMatrix(rmatrix)

      
    case (DIFFUSION_ISOTROPIC)
      ! Isotropic diffusion
      call gfsc_buildDiffusionOperator(&
          rproblemLevel%Rmatrix(CDEQ_MATRIX_S), .false., .true., rmatrix)
      
    case (DIFFUSION_ANISOTROPIC)
      ! Assemble high-order transport operator on coarser levels
      call gfsc_buildDiffusionOperator(&
          rproblemLevel%Rmatrix(CDEQ_MATRIX_S), .false., .true., rmatrix)
      
    case DEFAULT
      call output_line('Invalid type of diffusion!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcGalerkinResidual')
      call sys_halt()
    end select

    !---------------------------------------------------------------------------
    ! Assemble convective operator
    !---------------------------------------------------------------------------
    
    ! Set velocity vector for current level
!    call codire_setVelocity(rproblemLevel)

    ! Check if velocity is assumed to be discretely divergence free
    bisDivFree = (ivelocitytype .gt. 0)

    ! Assemble high-order transport operator on coarser levels
    select case(abs(ivelocitytype))
    case (VELOCITY_ZERO)
      ! zero velocity, do nothing
      
    case (VELOCITY_CONSTANT,&
          VELOCITY_TIMEDEP)
      ! linear velocity
      select case(rproblemLevel%rtriangulation%ndim)
      case (NDIM1D)
        select case(imode)
        case (1)
          call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
                                            rsolution, codire_calcPrimalConvConst1d,&
                                            bisDivFree, .false., .false., rmatrix)
        case (2)
          call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
                                            rsolution, codire_calcDualConvConst1d,&
                                            bisDivFree, .false., .false., rmatrix)
        end select
          
      case (NDIM2D)
        select case(imode)
          case (1)
            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
                                              rsolution, codire_calcPrimalConvConst2d,&
                                              bisDivFree, .false., .false., rmatrix)
          case (2)
            call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
                                              rsolution, codire_calcDualConvConst2d,&
                                              bisDivFree, .false., .false., rmatrix)
          end select
          
      case (NDIM3D)
        select case(imode)
        case (1)
          call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CZ),&
                                            rsolution, codire_calcPrimalConvConst3d,&
                                            bisDivFree, .false., .false., rmatrix)
        case (2)
          call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CZ),&
                                            rsolution, codire_calcDualConvConst3d,&
                                            bisDivFree, .false., .false., rmatrix)
        end select

      case DEFAULT
        call output_line('Invalid spatial dimension!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'codire_calcGalerkinResidual')
        call sys_halt()
      end select
      
      
    case (VELOCITY_BURGERS_SPACETIME)
      ! nonlinear Burgers' equation in space-time
      call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
                                        rsolution, codire_calcConvectionBurgersSpT2d,&
                                        bisDivFree, .false., .false., rmatrix)
      
    case (VELOCITY_BUCKLEV_SPACETIME)
      ! nonlinear Buckley-Leverett equation in space-time
      call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
                                        rsolution, codire_calcConvectionBuckLevSpT2d,&
                                        bisDivFree, .false., .false., rmatrix)
      
    case (VELOCITY_BURGERS1D)
      ! nonlinear Burgers' equation in 1D
      call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
                                        rsolution, codire_calcConvectionBurgers1d,&
                                        bisDivFree, .false., .false., rmatrix)

    case (VELOCITY_BURGERS2D)
      ! nonlinear Burgers' equation in 2D
      call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CY),&
                                        rsolution, codire_calcConvectionBurgers2d,&
                                        bisDivFree, .false., .false., rmatrix)
      
    case (VELOCITY_BUCKLEV1D)
      ! nonlinear Buckley-Leverett equation in 1D
      call gfsc_buildConvectionOperator(rproblemLevel%Rmatrix(CDEQ_MATRIX_CX:CDEQ_MATRIX_CX),&
                                        rsolution, codire_calcConvectionBuckLev1d,&
                                        bisDivFree, .false., .false., rmatrix)
      
    case DEFAULT
      call output_line('Invalid velocity profile!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_calcGalerkinResidual')
      call sys_halt()
    end select
    
    ! Compute the Galerkin residual
    if (present(rf)) then
      call lsysbl_copyVector(rf, rresidual)
      call lsyssc_scalarMatVec(rmatrix, rsolution%RvectorBlock(1),&
                               rresidual%RvectorBlock(1), -1._DP, 1._DP)
    else
      call lsyssc_scalarMatVec(rmatrix, rsolution%RvectorBlock(1),&
                               rresidual%RvectorBlock(1), -1._DP, 0._DP)
    end if

    ! Release memory
    call lsyssc_releaseMatrix(rmatrix)

  end subroutine codire_calcGalerkinResidual

  !*****************************************************************************
  
!<subroutine>

  subroutine codire_calcRHS(rproblemLevel, rtimestep, rsolver,&
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
    real(DP) :: dweight
    integer :: transportMatrix
    integer :: consistentMassMatrix, lumpedMassMatrix
    integer :: convectionAFC, diffusionAFC
    integer :: imasstype,imassantidiffusion


    ! Start time measurement for residual/rhs evaluation
    rtimer => collct_getvalue_timer(rcollection, 'timerAssemblyVector')
    call stat_startTimer(rtimer, STAT_TIMERSHORT)

    ! For nonlinear conservation laws, the global system operator
    ! needs to be updated in each nonlinear iteration. The update 
    ! routine is written such that it first determines if the problem
    ! is nonlinear and returns without matrix update otherwise.
    call codire_calcPreconditioner(rproblemLevel, rtimestep, rsolver, rsol, rcollection)


    ! Get parameters from collection which are required unconditionally
    transportMatrix      = collct_getvalue_int(rcollection, 'transportmatrix')
    consistentMassMatrix = collct_getvalue_int(rcollection, 'consistentmassmatrix')
    lumpedMassMatrix     = collct_getvalue_int(rcollection, 'lumpedmassmatrix')


    ! Compute the right-hand side
    !
    !   $ rhs = weight*(1-theta)*dt*L(u)*u $

    dweight = rtimestep%DmultistepWeights(istep)*&
              rtimestep%dStep*(1._DP-rtimestep%theta)
    call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                             rsol%rvectorBlock(1),&
                             rrhs%RvectorBlock(1), dweight, 0._DP)


    ! Perform algebraic flux correction for the convective term if required
    !
    !   $ rhs = rhs + f^*(u^n+1,u^n) $

    convectionAFC = collct_getvalue_int(rcollection, 'convectionAFC')
       
    ! What kind of stabilisation should be applied?
    select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
          
    case (AFCSTAB_FEMFCT,&
          AFCSTAB_FEMFCT_CLASSICAL)

      dweight = rtimestep%DmultistepWeights(istep)*rtimestep%dStep
      imassantidiffusion = collct_getvalue_int(rcollection, 'imassantidiffusion')

      ! Should we apply consistent mass antidiffusion?
      if (imassantidiffusion .eq. MASS_CONSISTENT) then
        call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                   rsol, rtimestep%theta, dweight, .true., rrhs,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(consistentMassMatrix))
      else
        call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                   rsol, rtimestep%theta, dweight, .true., rrhs,&
                                   rproblemLevel%Rafcstab(convectionAFC))
      end if
          
    case (AFCSTAB_FEMTVD)
      call gfsc_buildResidualTVD(rsol, dweight, rrhs,&
                                 rproblemLevel%Rafcstab(convectionAFC))

    case (AFCSTAB_FEMGP)
      call gfsc_buildResidualGP(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                rsol, rsol0, rtimestep%theta, dweight, rrhs,&
                                rproblemLevel%Rafcstab(convectionAFC))
    end select


    ! Perform algebraic flux correction for the diffusive term if required
    !
    !   $ rhs = rhs + g^*(u^n+1,u^n) $

    diffusionAFC = collct_getvalue_int(rcollection, 'diffusionAFC')
    
    ! What kind of stabilisation should be applied?
    select case(rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation)
          
    case (AFCSTAB_SYMMETRIC)
      call gfsc_buildResidualSymm(rsol, 1._DP, rrhs, rproblemLevel%Rafcstab(diffusionAFC))
    end select
    
    ! Stop time measurement for residual/rhs evaluation
    call stat_stopTimer(rtimer)
    
  end subroutine codire_calcRHS

  !*****************************************************************************

!<subroutine>

  subroutine codire_calcVelocityField(rappDescriptor, rproblemLevel,&
                                      dtime, rcollection, nlminOpt)

!<description>
    ! This subroutine calculates the velocity fields from the function
    ! parser. The result is stored separately for each problem level.
!</description>

!<input>
    ! application descriptor
    type(t_codire), intent(IN) :: rappDescriptor

    ! simulation time
    real(DP), intent(IN) :: dtime

    ! OPTIONAL: minimum problem level
    integer, intent(IN), optional :: nlminOpt
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT), target :: rproblemLevel

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>

    ! local variables
    type(t_problemLevel), pointer :: p_rproblemLevel
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(NDIM3D+1) :: Dvalue
    integer :: ieq, neq, idim, ndim, nlmin
    integer :: velocityfield


    ! Set minimum problem level
    nlmin = rproblemLevel%ilev
    if (present(nlminOpt)) nlmin = nlminOpt

    ! Get parameter from collection
    velocityfield = collct_getvalue_int(rcollection, 'velocityfield')

    ! Initialize variable values
    Dvalue = 0.0_DP


    ! Loop over all problem levels
    p_rproblemLevel => rproblemLevel
    do while(associated(p_rproblemLevel))

      ! Get number of degrees of freedom and spatial dimension
      neq  = p_rproblemLevel%rtriangulation%NVT
      ndim = p_rproblemLevel%rtriangulation%ndim

      ! Create/resize velocity vector if required
      if (p_rproblemLevel%RvectorBlock(velocityfield)%NEQ .eq. 0) then
        call lsysbl_createVectorBlock(&
            p_rproblemLevel%rvectorBlock(velocityfield), neq, ndim, .true.)
      elseif (p_rproblemLevel%RvectorBlock(velocityfield)%NEQ .ne. neq*ndim) then
        call lsysbl_resizeVectorBlock(&
            p_rproblemLevel%rvectorBlock(velocityfield), neq, .true.)
      end if

      ! Get vertex coordinates of the current problem level
      call storage_getbase_double2D(&
          p_rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)


      ! Loop over all spatial dimensions
      do idim = 1, ndim

        ! Get scalar subvector
        call lsyssc_getbase_double(&
            p_rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(idim), p_Ddata)

        ! Loop over all equations of scalar subvector
        do ieq = 1, neq
          
          Dvalue(1:ndim)   = p_DvertexCoords(:,ieq)
          Dvalue(NDIM3D+1) = dtime
          
          call fparser_evalFunction(rappDescriptor%rfparserVelocityField,&
                                    idim, Dvalue, p_Ddata(ieq))
        end do
      end do

      ! Specify velocity field as updated
      p_rproblemLevel%iproblemSpec = ior(p_rproblemLevel%iproblemSpec,&
                                         PROBLEV_MSPEC_UPDATE)
      
      ! Proceed to coarser problem level if minimum level has not been reached
      if (p_rproblemLevel%ilev .le. nlmin) exit
      p_rproblemLevel => p_rproblemLevel%p_rproblemLevelCoarse
      
    end do

  end subroutine codire_calcVelocityField

  !*****************************************************************************

!<subroutine>

  subroutine codire_setVelocityField(rvector)

!<description>
    ! This subroutine sets the global pointer to the velocity vector
    ! on the given problem level structure. Note that this subroutine
    ! will not work of multiple convection-diffusion-reaction problems
    ! are solved in parallel since there is only one global pointer.
!</description>

!<input>
    ! velocity field
    type(t_vectorBlock), intent(IN) :: rvector
!</input>

    ! What spatial dimension are we?
    select case(rvector%nblocks)
    case (NDIM1D)
      call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvelocityX)

    case (NDIM2D)
      call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvelocityY)

    case (NDIM3D)
      call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvelocityX)
      call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvelocityY)
      call lsyssc_getbase_double(rvector%RvectorBlock(3), p_DvelocityZ)

    case DEFAULT
      call output_line('Invalid spatial dimension!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'codire_setVelocityField')
      call sys_halt()
    end select
    
  end subroutine codire_setVelocityField

  !*****************************************************************************

!<subroutine>

  pure subroutine codire_calcPrimalConvConst1d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 1D
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij,k_ji
!</output>
!</subroutine>

    k_ij = -p_DvelocityX(j)*C_ij(1)
    k_ji = -p_DvelocityX(i)*C_ji(1)

  end subroutine codire_calcPrimalConvConst1d

  !*****************************************************************************

!<subroutine>

  pure subroutine codire_calcDualConvConst1d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 1D
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij,k_ji
!</output>
!</subroutine>

    k_ij = -p_DvelocityX(j)*C_ji(1)
    k_ji = -p_DvelocityX(i)*C_ij(1)

  end subroutine codire_calcDualConvConst1d

  !*****************************************************************************

!<subroutine>

  pure subroutine codire_calcPrimalConvConst2d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 2D
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij,k_ji
!</output>
!</subroutine>

    k_ij = -p_DvelocityX(j)*C_ij(1)-p_DvelocityY(j)*C_ij(2)
    k_ji = -p_DvelocityX(i)*C_ji(1)-p_DvelocityY(i)*C_ji(2)

  end subroutine codire_calcPrimalConvConst2d

  !*****************************************************************************

!<subroutine>
  
  pure subroutine codire_calcDualConvConst2d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 2D
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij,k_ji
!</output>
!</subroutine>

    k_ij = -p_DvelocityX(j)*C_ji(1)-p_DvelocityY(j)*C_ji(2)
    k_ji = -p_DvelocityX(i)*C_ij(1)-p_DvelocityY(i)*C_ij(2)

  end subroutine codire_calcDualConvConst2d

  !*****************************************************************************

!<subroutine>

  pure subroutine codire_calcPrimalConvConst3d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the primal problem in 3D
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij,k_ji
!</output>
!</subroutine>

    k_ij = -p_DvelocityX(j)*C_ij(1)-p_DvelocityY(j)*C_ij(2)-p_DvelocityZ(j)*C_ij(3)
    k_ji = -p_DvelocityX(i)*C_ji(1)-p_DvelocityY(i)*C_ji(2)-p_DvelocityZ(i)*C_ji(3)

  end subroutine codire_calcPrimalConvConst3d

  !*****************************************************************************

!<subroutine>

  pure subroutine codire_calcDualConvConst3d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for a constant velocity vector of the 
    ! form $v=v(x,y)$ or $v=v(x,y,t)$ for the dual problem in 3D
!</description>
    
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij,k_ji
!</output>
!</subroutine>

    k_ij = -p_DvelocityX(j)*C_ji(1)-p_DvelocityY(j)*C_ji(2)-p_DvelocityZ(j)*C_ji(3)
    k_ji = -p_DvelocityX(i)*C_ij(1)-p_DvelocityY(i)*C_ij(2)-p_DvelocityZ(i)*C_ij(3)

  end subroutine codire_calcDualConvConst3d

  !*****************************************************************************
    
!<subroutine>

  pure subroutine codire_calcConvectionBurgersSpT2d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the 
    ! one-dimensional Burgers equation $du/dt+df(u)/dx=0$, whereby
    ! the flux function is given by $f(u)=0.5*u^2$.
!</description>
   
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij,k_ji
!</output>
!</subroutine>

    k_ij = -0.5_DP*(u_i+u_j)*C_ij(1)-C_ij(2)
    k_ji = -0.5_DP*(u_i+u_j)*C_ji(1)-C_ji(2)

  end subroutine codire_calcConvectionBurgersSpT2d

  !*****************************************************************************
  
!<subroutine>

  pure subroutine codire_calcConvectionBuckLevSpT2d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for space-time formulation of the 
    ! Buckley-Leverett equation $du/dt+df(u)/dx=0$, whereby the
    ! flux function is given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
    !
    ! Here, the characteristic velocity $a(u)=f^\prime(u)$ is given
    ! by $a(u)=\frac{4u(1-u)}{(3u^2-2u+1)^2}$.
!</description>
   
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij, k_ji
!</output>
!</subroutine>

    ! local variables
    real(DP) :: v_i,v_j
    
    v_i = 4*u_i*(1-u_i)/(3*u_i*u_i-2*u_i+1)**2
    v_j = 4*u_j*(1-u_j)/(3*u_j*u_j-2*u_j+1)**2

    k_ij = -v_j*C_ij(1)-C_ij(2)
    k_ji = -v_i*C_ji(1)-C_ji(2)
        
  end subroutine codire_calcConvectionBuckLevSpT2d

  !*****************************************************************************
    
!<subroutine>

  pure subroutine codire_calcConvectionBurgers1d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for Burgers' equation in 1D.
!</description>
   
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij,k_ji
!</output>
!</subroutine>

    k_ij = -0.5_DP*(u_i+u_j)*C_ij(1)
    k_ji = -0.5_DP*(u_i+u_j)*C_ji(1)

  end subroutine codire_calcConvectionBurgers1d

  !*****************************************************************************
    
!<subroutine>

  pure subroutine codire_calcConvectionBurgers2d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji)

!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for Burgers' equation in 2D.
!</description>
   
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij,k_ji
!</output>
!</subroutine>

    k_ij = -0.5_DP*(u_i+u_j)*(C_ij(1)+C_ij(2))
    k_ji = -0.5_DP*(u_i+u_j)*(C_ji(1)+C_ji(2))

  end subroutine codire_calcConvectionBurgers2d

  !*****************************************************************************
  
!<subroutine>

  pure subroutine codire_calcConvectionBuckLev1d(u_i, u_j, C_ij, C_ji, i, j, k_ij, k_ji)
                                                
!<description>
    ! This subroutine computes the convective matrix coefficients
    ! $k_{ij}$ and $k_{ji}$ for the Buckley-Leverett equation 
    ! $du/dt+df(u)/dx=0$ in 1D, whereby the flux function is 
    ! given by $f(u)=u^2/(u^2+0.5*(1-u)^2)$
    !
    ! Here, the characteristic velocity $a(u)=f^\prime(u)$ is given
    ! by $a(u)=\frac{4u(1-u)}{(3u^2-2u+1)^2}$.
!</description>
   
!<input>
    ! solution vector
    real(DP), intent(IN) :: u_i, u_j

    ! coefficients from spatial discretization
    real(DP), dimension(:), intent(IN) :: C_ij, C_ji

    ! nodal indices
    integer, intent(IN) :: i, j
!</input>

!<output>
    ! convective coefficients
    real(DP), intent(OUT) :: k_ij,k_ji
!</output>
!</subroutine>

    k_ij = -(4*u_j*(1-u_j)/(3*u_j*u_j-2*u_j+1)**2)*C_ij(1)
    k_ji = -(4*u_i*(1-u_i)/(3*u_i*u_i-2*u_i+1)**2)*C_ji(1)
    
  end subroutine codire_calcConvectionBuckLev1d
  
  
end module codire_callback

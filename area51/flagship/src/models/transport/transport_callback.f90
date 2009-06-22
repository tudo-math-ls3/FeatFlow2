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
!# 1.) transp_coeffVectorAnalytic
!#     -> Callback routine for the evaluation of linear forms
!#        using an analytic expression for the load-vector
!#
!# 2.) transp_refFuncAnalytic
!#     -> Callback routine for the evaluation of the reference
!#        function for error estimation using an analytic expression
!#
!# 3.) transp_nlsolverCallback
!#     -> Callback routine for the nonlinear solver
!#
!# 4.) transp_calcPreconditioner
!#     -> Calculates the nonlinear preconditioner
!#
!# 5.) transp_calcJacobian
!#     -> Calculates the Jacobian matrix
!#
!# 6.) transp_applyJacobian
!#     -> Applies the Jacobian matrix to a given vector
!#
!# 7.) transp_calcResidual
!#     -> Calculates the nonlinear residual vector
!#
!# 8.) transp_calcRHS
!#     -> Calculates the constant right-hand side vector
!#
!# 9.) transp_setBoundary
!#     -> Imposes boundary conditions for nonlinear solver
!#
!# 10.) transp_calcVelocityField
!#      -> Calculates the velocity field
!#
!# 11.) transp_setVelocityField
!#      -> Sets the velocity field internally
!#
!# 12.) transp_calcLinearizedFCT
!#      -> Calculates the linearized FCT correction
!#
!# 13.) transp_buildVectorScalarBdr
!#      -> Builds the linear form defined in terms of a boundary integral
!#
!# </purpose>
!##############################################################################

module transport_callback

  use afcstabilisation
  use basicgeometry
  use boundaryfilter
  use collection
  use derivatives
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use groupfemscalar
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
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
  public :: transp_coeffVectorAnalytic
  public :: transp_refFuncAnalytic
  public :: transp_nlsolverCallback
  public :: transp_setBoundary
  public :: transp_calcPreconditioner
  public :: transp_calcJacobian
  public :: transp_applyJacobian
  public :: transp_calcResidual
  public :: transp_calcRHS
  public :: transp_calcVelocityField
  public :: transp_setVelocityField
  public :: transp_calcLinearizedFCT
  public :: transp_buildVectorScalarBdr

contains

  ! ***************************************************************************

!<subroutine>

  subroutine transp_coeffVectorAnalytic(rdiscretisation, rform, nelements,&
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
    real(DP), dimension(NDIM3D+1) :: Dvalue
    real(DP) :: dtime
    integer :: ipoint, ielement, ndim, icomp
    
    
    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    rfparser => collct_getvalue_pars(rcollection,&
                                     trim(rcollection%SquickAccess(1)))

    ! Moreover, this subroutine assumes tht the first quick access integer
    ! value holds the number of the function to be evaluated
    icomp = rcollection%IquickAccess(1)
   
    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    dtime = rcollection%DquickAccess(1)

    if (dtime < 0.0) then

      ! Evaluate all coefficients using the function parser
      do ielement = 1, nelements
        call fparser_evalFunction(rfparser, icomp, 2, Dpoints(:,:,ielement),&
                                  Dcoefficients(1,:,ielement))
      end do

    else

      ! Initialize values
      Dvalue = 0.0_DP
      Dvalue(NDIM3D+1) = dtime
      
      ! Set number of spatial dimensions
      ndim = size(Dpoints, 1)
      
      do ielement = 1, nelements
        do ipoint = 1, npointsPerElement
          
          ! Set values for function parser
          Dvalue(1:ndim)   = Dpoints(:, ipoint, ielement)
          
          ! Evaluate function parser
          call fparser_evalFunction(rfparser, icomp, Dvalue, Dcoefficients(1,ipoint,ielement))
        end do
      end do

    end if
    
  end subroutine transp_coeffVectorAnalytic
  
  !*****************************************************************************

!<subroutine>

  subroutine transp_refFuncAnalytic(cderivative, rdiscretisation, nelements,&
                                    npointsPerElement, Dpoints, IdofsTest,&
                                    rdomainIntSubset, Dvalues, rcollection)   
    use basicgeometry
    use collection
    use domainintegration
    use scalarpde
    use triangulation
    
!<description>
    ! This subroutine is called during the calculation of errors. It has to compute
    ! the (analytical) values of a function in a couple of points on a couple
    ! of elements. These values are compared to those of a computed FE function
    ! and used to calculate an error.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points.
!</description>
    
!<input>
    ! This is a DER_xxxx derivative identifier (from derivative.f90) that
    ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
    ! The result must be written to the Dvalue-array below.
    integer, intent(IN) :: cderivative
  
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! DIMENSION(NDIM2D,npointsPerElement,nelements)
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(:,:,:), intent(IN) :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in trial space,Number of elements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional :: rcollection   
!</input>
  
!<output>
    ! This array has to receive the values of the (analytical) function
    ! in all the points specified in Dpoints, or the appropriate derivative
    ! of the function, respectively, according to cderivative.
    !   DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), intent(OUT) :: Dvalues
!</output>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: rfparser
    real(DP), dimension(NDIM3D+1) :: Dvalue
    integer :: ipoint, ielement, ndim, icomp


    ! Initialize values
    Dvalue = 0.0_DP
    
    ! This subroutine assumes that the first quick access string
    ! value holds the name of the function parser in the collection.
    rfparser => collct_getvalue_pars(rcollection,&
                                     trim(rcollection%SquickAccess(1)))
   
    ! Moreover, this subroutine assumes that the first quick access integer
    ! value holds the number of the function to be evaluated
    select case(cderivative)
    case(1)
      icomp = rcollection%IquickAccess(1)
    case(0)
      icomp = rcollection%IquickAccess(2)
    case default
      print *, "Not implemented"
      stop
    end select

    ! This subroutine also assumes that the first quick access double
    ! value holds the simulation time
    Dvalue(NDIM3D+1) = rcollection%DquickAccess(1)
    
    ! Set number of spatial dimensions
    ndim = size(Dpoints, 1)

    do ielement = 1, nelements
      do ipoint = 1, npointsPerElement
        
        ! Set values for function parser
        Dvalue(1:ndim)   = Dpoints(:, ipoint, ielement)
        
        ! Evaluate function parser
        call fparser_evalFunction(rfparser, icomp, Dvalue, Dvalues(ipoint,ielement))
      end do
    end do
    
  end subroutine transp_refFuncAnalytic

  !*****************************************************************************

!<subroutine>

  subroutine transp_nlsolverCallback(rproblemLevel, rtimestep, rsolver,&
                                     rsolution, rsolutionInitial,&
                                     rrhs, rres, istep, ioperationSpec,&
                                     rcollection, istatus, rb)

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

    ! OPTIONAL: constant right-hand side vector
    type(t_vectorBlock), intent(IN), optional :: rb
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

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    integer :: jacobianMatrix

    
    ! Do we have to calculate the preconditioner?
    ! --------------------------------------------------------------------------
    if ((iand(ioperationSpec, NLSOL_OPSPEC_CALCPRECOND)  .ne. 0) .or.&
        (iand(ioperationSpec, NLSOL_OPSPEC_CALCRHS)      .ne. 0) .or.&
        (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0)) then
     
      call transp_calcPreconditioner(rproblemLevel, rtimestep, rsolver,&
                                     rsolution, rcollection)
    end if

    
    ! Do we have to calculate the constant right-hand side?
    ! --------------------------------------------------------------------------
    if ((iand(ioperationSpec, NLSOL_OPSPEC_CALCRHS)  .ne. 0)) then

      call transp_calcRhs(rproblemLevel, rtimestep, rsolver,&
                          rsolution, rsolutionInitial,&
                          rrhs, istep, rcollection)
    end if
          
    
    ! Do we have to calculate the residual
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

      call transp_calcResidual(rproblemLevel, rtimestep, rsolver,&
                               rsolution, rsolutionInitial,&
                               rrhs, rres, istep, rcollection, rb)
     end if
    
    
    ! Do we have to calculate the Jacobian operator?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCJACOBIAN) .ne. 0) then
      
      call transp_calcJacobian(rproblemLevel, rtimestep, rsolver,&
                               rsolution, rsolutionInitial, rcollection)
    end if
    
    
    ! Do we have to impose boundary conditions?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then
      
      call transp_setBoundary(rproblemLevel, rtimestep, rsolver,&
                              rsolution, rsolutionInitial, rres, rcollection)
    end if
    

    ! Do we have to apply the Jacobian operator?
    ! --------------------------------------------------------------------------
    if (iand(ioperationSpec, NLSOL_OPSPEC_APPLYJACOBIAN) .ne. 0) then

      p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'jacobianMatrix', jacobianMatrix)
      call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(jacobianMatrix),&
                               rsolution%RvectorBlock(1),&
                               rres%RvectorBlock(1), 1.0_DP, 1.0_DP)
    end if
    
    
    ! Set status flag
    istatus = 0
    
  end subroutine transp_nlsolverCallback
  
  !*****************************************************************************

!<subroutine>

  subroutine transp_calcPreconditioner(rproblemLevel, rtimestep, rsolver,&
                                       rsolution, rcollection)

!<description>
    ! This subroutine calculates the nonlinear preconditioner and
    ! configures the linear solver structure accordingly.
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
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    character(LEN=SYS_STRLEN) :: smode
    logical :: bStabilize, bbuildAFC, bconservative
    integer :: systemMatrix, transportMatrix
    integer :: lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY, coeffMatrix_CZ, coeffMatrix_S
    integer :: imasstype, ivelocitytype, idiffusiontype
    integer :: convectionAFC, diffusionAFC, velocityfield
    integer :: primaldual

    
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
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'transportmatrix', transportMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'coeffMatrix_CX', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'coeffMatrix_CY', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'coeffMatrix_CZ', coeffMatrix_CZ)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'coeffMatrix_S', coeffMatrix_S)
    
    !---------------------------------------------------------------------------
    ! Assemble diffusion operator
    !---------------------------------------------------------------------------
    
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'idiffusiontype', idiffusiontype)
    
    ! Primal and dual mode are equivalent
    select case(idiffusiontype)
    case (DIFFUSION_ZERO)
      ! zero diffusion, clear the transport matrix
      call lsyssc_clearMatrix(rproblemLevel%Rmatrix(transportMatrix))
      
    case (DIFFUSION_ISOTROPIC)
      ! Isotropic diffusion
      call gfsc_buildDiffusionOperator(rproblemLevel%Rmatrix(coeffMatrix_S),&
                                       .false., .true.,&
                                       rproblemLevel%Rmatrix(transportMatrix))

    case (DIFFUSION_ANISOTROPIC)
      ! Anisotropic diffusion
      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'diffusionAFC', diffusionAFC)

      if (diffusionAFC > 0) then
        
        ! What kind of stabilisation should be applied?
        select case(rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation)

        case (AFCSTAB_DMP)
          ! Satisfy discrete maximum principle
          call gfsc_buildDiffusionOperator(rproblemLevel%Rmatrix(coeffMatrix_S),&
                                           .true., .true.,&
                                           rproblemLevel%Rmatrix(transportMatrix))
          
        case (AFCSTAB_SYMMETRIC)
          ! Satisfy discrete maximum principle and assemble stabilization structure
          call gfsc_buildDiffusionOperator(rproblemLevel%Rmatrix(coeffMatrix_S),&
                                           .true., .true.,&
                                           rproblemLevel%Rmatrix(transportMatrix),&
                                           rproblemLevel%Rafcstab(diffusionAFC))

        case DEFAULT
          ! Compute the standard Galerkin approximation
          call gfsc_buildDiffusionOperator(rproblemLevel%Rmatrix(coeffMatrix_S),&
                                           .false., .true.,&
                                           rproblemLevel%Rmatrix(transportMatrix))
        end select

      else   ! diffusionAFC > 0

        ! Compute the standard Galerkin approximation
        call gfsc_buildDiffusionOperator(rproblemLevel%Rmatrix(coeffMatrix_S),&
                                         .false., .true.,&
                                         rproblemLevel%Rmatrix(transportMatrix))

      end if   ! diffusionAFC > 0


    case (DIFFUSION_VARIABLE)
      call output_line('Variable diffusion matrices are yet not implemented!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_calcPreconditioner')
      call sys_halt()

      ! Set update notification in problem level structure
      rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                       PROBLEV_MSPEC_UPDATE)
      
    case DEFAULT
      call output_line('Invalid type of diffusion!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_calcPreconditioner')
      call sys_halt()
    end select
    
    
    !---------------------------------------------------------------------------
    ! Assemble convective operator
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(p_rparlist, rcollection%SquickAccess(1),&
                                'mode', smode)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'convectionAFC', convectionAFC)

    if (convectionAFC > 0) then

      ! Check if stabilization should be applied
      bStabilize = AFCSTAB_GALERKIN .ne. &
                   rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation

      ! Check if stabilization structure should be built
      bbuildAFC = bStabilize .and. AFCSTAB_UPWIND .ne.&
                  rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation
      
    else   ! convectionAFC > 0

      bStabilize = .false.
      bbuildAFC  = .false.

    end if   ! convectionAFC > 0

    ! Check if conservative or non-conservative formulation should be applied
    bconservative = ivelocitytype .gt. 0


    ! Are we in primal or dual mode?
    select case(trim(smode))
    case('primal')
      
      select case(abs(ivelocitytype))
      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP)
        ! linear velocity

        call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                                 'velocityfield', velocityfield)
        call transp_setVelocityField(rproblemLevel%RvectorBlock(velocityfield))
        
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

        if (abs(ivelocitytype) .eq. VELOCITY_TIMEDEP) then
          ! Set update notification in problem level structure
          rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                           PROBLEV_MSPEC_UPDATE)
        end if
        

      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers' equation in space-time

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
              rsolution, transp_calcMatrixPrimalBurgersSpT2d,&
              bStabilize, .false., rproblemLevel%Rmatrix(transportMatrix),&
              bisConservative = bconservative)

        end if

        ! Set update notification in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         PROBLEV_MSPEC_UPDATE)


      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers' equation in 1D

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

        ! Set update notification in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         PROBLEV_MSPEC_UPDATE)
        

      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers' equation in 2D

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

        ! Set update notification in problem level structure
        rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                         PROBLEV_MSPEC_UPDATE)


      case DEFAULT
        call output_line('Invalid velocity profile!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'transp_calcPreconditioner')
        call sys_halt()
      end select

      
    case('dual')
      
      select case(abs(ivelocitytype))
      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        

      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP)
        ! linear velocity

        call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                                 'velocityfield', velocityfield)
        call transp_setVelocityField(rproblemLevel%RvectorBlock(velocityfield))
        
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

        if (abs(ivelocitytype) .eq. VELOCITY_TIMEDEP) then
          ! Set update notification in problem level structure
          rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                       PROBLEV_MSPEC_UPDATE)
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
                         OU_CLASS_ERROR,OU_MODE_STD,'transp_calcPreconditioner')
        call sys_halt()
      end select

    case DEFAULT
      call output_line('Invalid mode!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_calcPreconditioner')
      call sys_halt()
    end select
        
    !---------------------------------------------------------------------------
    ! Assemble the global system operator
    !---------------------------------------------------------------------------
    
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'systemmatrix', systemMatrix)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'imasstype', imasstype)

    select case(imasstype)
    case (MASS_LUMPED)

      ! Compute the global operator for transient flow
      !
      !   $ A = ML-theta*dt*L $

      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'lumpedmassmatrix', lumpedMassMatrix)

      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(lumpedMassMatrix), 1.0_DP,&
                                   rproblemLevel%Rmatrix(transportMatrix),&
                                   -rtimestep%theta*rtimestep%dStep,&
                                   rproblemLevel%Rmatrix(systemMatrix),&
                                   .false., .false., .true., .true.)
    case (MASS_CONSISTENT)

      ! Compute the global operator for transient flow
      !
      !   $ A = MC-theta*dt*L $

      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'consistentmassmatrix', consistentMassMatrix)

      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(consistentMassMatrix), 1.0_DP,&
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
      call lsyssc_scaleMatrix(rproblemLevel%Rmatrix(systemMatrix), -1.0_DP)
    
    end select
    

    !---------------------------------------------------------------------------
    ! Impose boundary conditions
    !---------------------------------------------------------------------------
      
    call bdrf_filterMatrix(rsolver%rboundaryCondition,&
                           rproblemLevel%rtriangulation,&
                           rproblemLevel%Rmatrix(systemMatrix), 1.0_DP)
        
    ! Ok, we updated the (nonlinear) system operator successfully. Now we still 
    ! have to link it to the solver hierarchy. This is done recursively.
    call flagship_updateSolverMatrix(rproblemLevel, rsolver, systemMatrix,&
                                     SYSTEM_INTERLEAVEFORMAT, UPDMAT_ALL,&
                                     rproblemLevel%ilev, rproblemLevel%ilev)

    ! Finally, we have to update the content of the solver hierarchy
    call solver_updateContent(rsolver)

    ! Stop time measurement for global operator
    call stat_stopTimer(p_rtimer)
    
  end subroutine transp_calcPreconditioner

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcJacobian(rproblemLevel, rtimestep, rsolver,&
                                 rsolution, rsolutionInitial, rcollection)

!<description>
    ! This callback subroutine computes the Jacobian matrix.
!</description>

!<input>
    ! time-stepping algorithm
    type(t_timestep), intent(IN) :: rtimestep

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionInitial
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! solution vector
    type(t_vectorBlock), intent(INOUT) :: rsolution

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
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
    integer :: primaldual, velocityfield
    logical :: bStabilize, bisExactStructure, bisExtendedSparsity


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
      hstep = max(SYS_EPSREAL, rsolver%p_solverNewton%dperturbationStrategy)
    end select
    

    !---------------------------------------------------------------------------
    ! Assemble diffusion operator
    !---------------------------------------------------------------------------

    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'idiffusiontype', idiffusiontype)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'diffusionAFC', diffusionAFC)

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
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
      call sys_halt()
    end select
    
    
    !---------------------------------------------------------------------------
    ! Assemble convection operator
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(p_rparlist, rcollection%SquickAccess(1),&
                                'mode', smode)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'ivelocitytype', ivelocitytype)
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'convectionAFC', convectionAFC)

    ! Check if stabilization should be applied
    bStabilize = (AFCSTAB_GALERKIN .ne.&
                  rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
    
    ! Are we in primal or dual mode?
    select case(trim(smode))
    case('primal')

      select case (abs(ivelocitytype))
      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        
      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP) 
        ! linear velocity
        
        call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                                 'velocityfield', velocityfield)
        call transp_setVelocityField(rproblemLevel%RvectorBlock(velocityfield))

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsc_buildConvectionJacobian(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rsolution, transp_calcMatrixPrimalConst1d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))

        case (NDIM2D)
          call gfsc_buildConvectionJacobian(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsolution, transp_calcMatrixPrimalConst2d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))

        case (NDIM3D)
          call gfsc_buildConvectionJacobian(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rsolution, transp_calcMatrixPrimalConst3d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))
        end select
      
      case (VELOCITY_BURGERS_SPACETIME)
        ! nonlinear Burgers' equation in space-time
        call gfsc_buildConvectionJacobian(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
            rsolution, transp_calcMatrixPrimalBurgersSpT2d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))
        
      case (VELOCITY_BUCKLEV_SPACETIME)
        ! nonlinear Buckley-Leverett equation in space-time
        call gfsc_buildConvectionJacobian(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
            transp_calcMatrixPrimalBuckLevSpT2d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))
        
      case (VELOCITY_BURGERS1D)
        ! nonlinear Burgers' equation in 1D
        call gfsc_buildConvectionJacobian(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
            transp_calcMatrixPrimalBurgers1d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))
        
      case (VELOCITY_BURGERS2D)
        ! nonlinear Burgers' equation in 2D
        call gfsc_buildConvectionJacobian(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY), rsolution,&
            transp_calcMatrixPrimalBurgers2d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))
        
      case (VELOCITY_BUCKLEV1D)
        ! nonlinear Buckley-Leverett equation in 1D
        call gfsc_buildConvectionJacobian(&
            rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX), rsolution,&
            transp_calcMatrixPrimalBuckLev1d, hstep, bStabilize,&
            .false., rproblemLevel%Rmatrix(transportMatrix))
        
      case DEFAULT
        call output_line('Unsupported velocity type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
        call sys_halt()
      end select


    case('dual')
      
      select case (abs(ivelocitytype))
      case (VELOCITY_ZERO)
        ! zero velocity, do nothing
        
      case (VELOCITY_CONSTANT,&
            VELOCITY_TIMEDEP) 
        ! linear velocity

        call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                                 'velocityfield', velocityfield)
        call transp_setVelocityField(rproblemLevel%RvectorBlock(velocityfield))

        select case(rproblemLevel%rtriangulation%ndim)
        case (NDIM1D)
          call gfsc_buildConvectionJacobian(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
              rsolution, transp_calcMatrixDualConst1d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))

        case (NDIM2D)
          call gfsc_buildConvectionJacobian(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
              rsolution, transp_calcMatrixDualConst2d, hstep, bStabilize,&
              .false., rproblemLevel%Rmatrix(transportMatrix))

        case (NDIM3D)
          call gfsc_buildConvectionJacobian(&
              rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CZ),&
              rsolution, transp_calcMatrixDualConst3d, hstep, bStabilize,&
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
                         OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
        call sys_halt()
      end select

    case DEFAULT
      call output_line('Invalid mode!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
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

    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'imasstype', imasstype)

    select case(imasstype)
    case (MASS_LUMPED)
      ! Compute the global Jacobian for transient flow
      !
      !   $ J = ML-theta*dt*L $
      
      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(transportMatrix),&
                                   -rtimestep%theta*rtimestep%dStep,&
                                   rproblemLevel%Rmatrix(lumpedMassMatrix), 1.0_DP,&
                                   rproblemLevel%Rmatrix(jacobianMatrix),&
                                   .false., .false., .true., bisExactStructure)
    case (MASS_CONSISTENT)
      ! Compute the global Jacobian for transient flow
      !
      !   $ J = MC-theta*dt*L $
      
      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(transportMatrix),&
                                   -rtimestep%theta*rtimestep%dStep,&
                                   rproblemLevel%Rmatrix(consistentMassMatrix), 1.0_DP,&
                                   rproblemLevel%Rmatrix(jacobianMatrix),&
                                   .false., .false., .true., bisExactStructure)
    case DEFAULT
      ! Compute the global Jacobian for steady-state flow
      !
      !   $ J = -L $
      
      call lsyssc_MatrixLinearComb(rproblemLevel%Rmatrix(transportMatrix), -1.0_DP,&
                                   rproblemLevel%Rmatrix(jacobianMatrix), 0.0_DP,&
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
        call gfsc_buildJacobianSymm(rsolution, 1.0_DP, hstep, .false.,&
                                    rproblemLevel%Rafcstab(diffusionAFC),&
                                    rproblemLevel%Rmatrix(jacobianMatrix),&
                                    bisExtendedSparsity)
      end select

    case (DIFFUSION_VARIABLE)
      print *, "Variable diffusion matrices are yet not implemented!"
      stop
      
    case DEFAULT
      call output_line('Invalid type of diffusion!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
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
        
        call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                                 'imassantidiffusiontype', imassantidiffusiontype)

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsc_buildJacobianFCT(rsolution, rtimestep%theta, rtimestep%dStep, hstep,&
                                     .false., rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildJacobianFCT(rsolution, rtimestep%theta, rtimestep%dStep, hstep,&
                                     .false., rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildJacobianTVD(rsolution, rtimestep%dStep, hstep, .false.,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(jacobianMatrix),&
                                   bisExtendedSparsity)

      case (AFCSTAB_FEMGP)
        call gfsc_buildJacobianGP(rproblemLevel%Rmatrix(consistentMassMatrix), rsolution,&
                                  rsolutionInitial, rtimestep%theta, rtimestep%dStep, hstep,&
                                  .false., rproblemLevel%Rafcstab(convectionAFC),&
                                  rproblemLevel%Rmatrix(jacobianMatrix),&
                                  bisExtendedSparsity)
      end select

    case(VELOCITY_BURGERS_SPACETIME)
      ! nonlinear Burgers' equation in space-time

      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)

        call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                                 'imassantidiffusiontype', imassantidiffusiontype)

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, transp_calcMatrixPrimalBurgersSpT2d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, transp_calcMatrixPrimalBurgersSpT2d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsolution, transp_calcMatrixPrimalBurgersSpT2d,&
                                   rtimestep%dStep, hstep, .false.,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(jacobianMatrix),&
                                   bisExtendedSparsity)

      case (AFCSTAB_FEMGP)
        call gfsc_buildJacobianGP(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                                  rsolution, rsolutionInitial, transp_calcMatrixPrimalBurgersSpT2d,&
                                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                  rproblemLevel%Rafcstab(convectionAFC),&
                                  rproblemLevel%Rmatrix(jacobianMatrix),&
                                  bisExtendedSparsity)
      end select

    case(VELOCITY_BUCKLEV_SPACETIME)
      ! nonlinear Buckley-Leverett equation in space-time

      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)

        call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                                 'imassantidiffusiontype', imassantidiffusiontype)

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, transp_calcMatrixPrimalBuckLevSpT2d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, transp_calcMatrixPrimalBuckLevSpT2d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsolution, transp_calcMatrixPrimalBuckLevSpT2d,&
                                   rtimestep%dStep, hstep, .false.,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(jacobianMatrix),&
                                   bisExtendedSparsity)
        
      case (AFCSTAB_FEMGP)
        call gfsc_buildJacobianGP(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                                  rsolution, rsolutionInitial, transp_calcMatrixPrimalBuckLevSpT2d,&
                                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                  rproblemLevel%Rafcstab(convectionAFC),&
                                  rproblemLevel%Rmatrix(jacobianMatrix),&
                                  bisExtendedSparsity)
      end select
      
    case(VELOCITY_BURGERS1D)
      ! nonlinear Burgers' equation in 1D
      
      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)

        call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                                 'imassantidiffusiontype', imassantidiffusiontype)

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, transp_calcMatrixPrimalBurgers1d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, transp_calcMatrixPrimalBurgers1d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                   rsolution, transp_calcMatrixPrimalBurgers1d,&
                                   rtimestep%dStep, hstep, .false.,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(jacobianMatrix),&
                                   bisExtendedSparsity)

      case (AFCSTAB_FEMGP)
        call gfsc_buildJacobianGP(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                                  rsolution, rsolutionInitial, transp_calcMatrixPrimalBurgers1d,&
                                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                  rproblemLevel%Rafcstab(convectionAFC),&
                                  rproblemLevel%Rmatrix(jacobianMatrix),&
                                  bisExtendedSparsity)
      end select

    case(VELOCITY_BURGERS2D)
      ! nonlinear Burgers' equation in 2D
      
      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)

        call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                                 'imassantidiffusiontype', imassantidiffusiontype)

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, transp_calcMatrixPrimalBurgers2d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                     rsolution, transp_calcMatrixPrimalBurgers2d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                   rsolution, transp_calcMatrixPrimalBurgers2d,&
                                   rtimestep%dStep, hstep, .false.,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(jacobianMatrix),&
                                   bisExtendedSparsity)

      case (AFCSTAB_FEMGP)
        call gfsc_buildJacobianGP(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CY),&
                                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                                  rsolution, rsolutionInitial, transp_calcMatrixPrimalBurgers2d,&
                                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                  rproblemLevel%Rafcstab(convectionAFC),&
                                  rproblemLevel%Rmatrix(jacobianMatrix),&
                                  bisExtendedSparsity)
      end select
      
    case(VELOCITY_BUCKLEV1D)
      ! nonlinear Buckley-Leverett equation in 1D

      select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
      case (AFCSTAB_FEMFCT,&
            AFCSTAB_FEMFCT_CLASSICAL)

        call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                                 'imassantidiffusiontype', imassantidiffusiontype)

        ! Should we apply consistent mass antidiffusion?
        if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, transp_calcMatrixPrimalBuckLev1d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix),&
                                     rproblemLevel%Rmatrix(consistentMassMatrix))
        else
          call gfsc_buildJacobianFCT(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                     rsolution, transp_calcMatrixPrimalBuckLev1d,&
                                     rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                     rproblemLevel%Rafcstab(convectionAFC),&
                                     rproblemLevel%Rmatrix(jacobianMatrix))
        end if
        
      case (AFCSTAB_FEMTVD)
        call gfsc_buildJacobianTVD(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                   rsolution, transp_calcMatrixPrimalBuckLev1d,&
                                   rtimestep%dStep, hstep, .false.,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(jacobianMatrix),&
                                   bisExtendedSparsity)

      case (AFCSTAB_FEMGP)
        call gfsc_buildJacobianGP(rproblemLevel%Rmatrix(coeffMatrix_CX:coeffMatrix_CX),&
                                  rproblemLevel%Rmatrix(consistentMassMatrix),&
                                  rsolution, rsolutionInitial, transp_calcMatrixPrimalBuckLev1d,&
                                  rtimestep%theta, rtimestep%dStep, hstep, .false.,&
                                  rproblemLevel%Rafcstab(convectionAFC),&
                                  rproblemLevel%Rmatrix(jacobianMatrix),&
                                  bisExtendedSparsity)
      end select
      
    case DEFAULT
      call output_line('Unsupported velocity type!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_calcJacobian')
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
      call flagship_updateSolverMatrix(rproblemLevel, rsolver, jacobianMatrix,&
                                       SYSTEM_INTERLEAVEFORMAT, UPDMAT_JAC_TRANSIENT,&
                                       rproblemLevel%ilev, rproblemLevel%ilev)

    case DEFAULT
      call flagship_updateSolverMatrix(rproblemLevel, rsolver, jacobianMatrix,&
                                       SYSTEM_INTERLEAVEFORMAT, UPDMAT_JAC_STEADY,&
                                       rproblemLevel%ilev, rproblemLevel%ilev)
    end select
    
    ! Finally, we have to update the content of the solver hierarchy
    call solver_updateContent(rsolver)
    
    ! Stop time measurement for matrix evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine transp_calcJacobian

  !*****************************************************************************

!<subroutine>

  subroutine transp_applyJacobian(rproblemLevel, rx, ry, cx, cy, rcollection)

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
    type(t_parlist), pointer :: p_rparlist
    integer :: jacobianMatrix
    
    ! Get parameters from parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'jacobianMatrix', jacobianMatrix)
    
    ! We know where the Jacobian matrix is stored and can apply it
    ! by means of matrix vector multiplication
    call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(jacobianMatrix),&
                             rx%RvectorBlock(1), ry%RvectorBlock(1), cx, cy)

  end subroutine transp_applyJacobian

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcRHS(rproblemLevel, rtimestep, rsolver,&
                            rsolution, rsolutionInitial,&
                            rrhs, istep, rcollection)

!<description>
    ! This subroutine computes the right-hand side vector
    ! used in the explicit Lax-Wendroff time-stepping scheme
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(IN) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(IN) :: rsolution

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

    ! right-hand side vector
    type(t_vectorBlock), intent(INOUT) :: rrhs

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
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


    ! Start time measurement for residual/rhs evaluation
    p_rtimer => collct_getvalue_timer(rcollection, 'rtimerAssemblyVector')
    call stat_startTimer(p_rtimer, STAT_TIMERSHORT)

    ! For nonlinear conservation laws, the global system operator
    ! needs to be updated in each nonlinear iteration. The update 
    ! routine is written such that it first determines if the problem
    ! is nonlinear and returns without matrix update otherwise.
    call transp_calcPreconditioner(rproblemLevel, rtimestep, rsolver, rsolution, rcollection)


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
    !   $ rhs = weight*(1-theta)*dt*L(u)*u $

    dweight = rtimestep%DmultistepWeights(istep)*&
              rtimestep%dStep*(1.0_DP-rtimestep%theta)
    call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                             rsolution%rvectorBlock(1),&
                             rrhs%RvectorBlock(1), dweight, 0.0_DP)


    ! Perform algebraic flux correction for the convective term if required
    !
    !   $ rhs = rhs + f^*(u^n+1,u^n) $

    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'convectionAFC', convectionAFC)
       
    ! What kind of stabilisation should be applied?
    select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
          
    case (AFCSTAB_FEMFCT,&
          AFCSTAB_FEMFCT_CLASSICAL)

      dweight = rtimestep%DmultistepWeights(istep)*rtimestep%dStep
      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'imassantidiffusiontype', imassantidiffusiontype)

      ! Should we apply consistent mass antidiffusion?
      if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
        call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                   rsolution, rtimestep%theta, dweight, .true., rrhs,&
                                   rproblemLevel%Rafcstab(convectionAFC),&
                                   rproblemLevel%Rmatrix(consistentMassMatrix))
      else
        call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                   rsolution, rtimestep%theta, dweight, .true., rrhs,&
                                   rproblemLevel%Rafcstab(convectionAFC))
      end if
          
    case (AFCSTAB_FEMTVD)
      call gfsc_buildResidualTVD(rsolution, dweight, rrhs,&
                                 rproblemLevel%Rafcstab(convectionAFC))

    case (AFCSTAB_FEMGP)
      call gfsc_buildResidualGP(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                rsolution, rsolutionInitial, rtimestep%theta, dweight, rrhs,&
                                rproblemLevel%Rafcstab(convectionAFC))
    end select


    ! Perform algebraic flux correction for the diffusive term if required
    !
    !   $ rhs = rhs + g^*(u^n+1,u^n) $

    call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                             'diffusionAFC', diffusionAFC)
    
    ! What kind of stabilisation should be applied?
    select case(rproblemLevel%Rafcstab(diffusionAFC)%ctypeAFCstabilisation)
          
    case (AFCSTAB_SYMMETRIC)
      call gfsc_buildResidualSymm(rsolution, 1.0_DP, rrhs,&
                                  rproblemLevel%Rafcstab(diffusionAFC))
    end select
    
    ! Stop time measurement for residual/rhs evaluation
    call stat_stopTimer(p_rtimer)

  end subroutine transp_calcRHS
    
  !*****************************************************************************

!<subroutine>

  subroutine transp_calcResidual(rproblemLevel, rtimestep, rsolver,&
                                 rsolution, rsolutionInitial, rrhs, rres,&
                                 ite, rcollection, rb)

!<description>
    ! This subroutine computes the nonlinear residual vector
    ! 
    !   $$ res^{(m)} = rhs - [M-\theta\Delta t K^{(m)}]u^{(m)} $$
    !
    !  and the constant right-hand side (only in the zeroth iteration)
    !
    !  $$ rhs = [M + (1-\theta)\Delta t K^n]u^n + \Delta t b$$
!</description>

!<input>
    ! time-stepping structure
    type(t_timestep), intent(IN) :: rtimestep

    ! solution vector
    type(t_vectorBlock), intent(IN) :: rsolution

    ! initial solution vector
    type(t_vectorBlock), intent(IN) :: rsolutionInitial

    ! number of nonlinear iteration
    integer, intent(IN) :: ite

    ! OPTIONAL: load vector specified by the application
    type(t_vectorBlock), intent(IN), optional :: rb
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(INOUT) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(INOUT) :: rsolver

    ! right-hand side vector
    !   ite=0: the right-hand side vector is calculated
    !   ite>0: the right-hand side vector remains unchanged
    type(t_vectorBlock), intent(INOUT) :: rrhs

    ! residual vector
    !   ite=0: the initial residual vector is calculated
    !   ite>0: the residual vector is initialized by the right-hand
    !          sude vector and updated by the implicit contribution
    type(t_vectorBlock), intent(INOUT) :: rres

    ! collection
    type(t_collection), intent(INOUT) :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_timer), pointer :: p_rtimer
    integer :: transportMatrix, consistentMassMatrix, lumpedMassMatrix
    integer :: convectionAFC, diffusionAFC
    integer :: imasstype, imassantidiffusiontype

    
    ! Start time measurement for residual/rhs evaluation
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

    ! Are we in the zero-th iteration?
    if (ite .eq. 0) then

      !-------------------------------------------------------------------------
      ! In the zero-th nonlinear iteration we calculate the initial
      ! residual vector and the constant right-hand side vector
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'imasstype', imasstype)
      
      select case(imasstype)
      case (MASS_LUMPED)
        
        ! Compute the initial low-order residual
        !
        !   $  res^{(0)} = dt*L(u^n)*u^n $
        
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                                 rsolution%rvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 rtimestep%dStep, 0.0_DP)

        ! Compute the constant right-hand side
        !
        !   $ rhs = M_L*u^n+(1-theta)*res^{(0)} $

        if (rtimestep%theta .lt. 1.0_DP) then
          call lsysbl_copyVector(rres, rrhs)

          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                   rsolution%RvectorBlock(1),&
                                   rrhs%RvectorBlock(1), 1.0_DP, 1.0_DP-rtimestep%theta)
        else
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                   rsolution%RvectorBlock(1),&
                                   rrhs%RvectorBlock(1),&
                                   1.0_DP, 0.0_DP)
        end if

      case (MASS_CONSISTENT)

        ! Compute the initial low-order residual
        !
        !   $  res^{(0)} = dt*L(u^n)*u^n $
        
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                                 rsolution%rvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 rtimestep%dStep, 0.0_DP)
        
        ! Compute the constant right-hand side
        !
        !   $ rhs = M_C*u^n+(1-theta)*res^{(0)} $

        if (rtimestep%theta .lt. 1.0_DP) then
          call lsysbl_copyVector(rres, rrhs)

          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                   rsolution%RvectorBlock(1),&
                                   rrhs%RvectorBlock(1), 1.0_DP, 1.0_DP-rtimestep%theta)
        else
          call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                   rsolution%RvectorBlock(1),&
                                   rrhs%RvectorBlock(1),&
                                   1.0_DP, 0.0_DP)
        end if

      case DEFAULT

        ! Compute the initial low-order residual
        !
        !   $  res^{(0)} = L(u^n)*u^n $

        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                                 rsolution%rvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 1.0_DP, 0.0_DP)        
      end select
     

      ! Perform algebraic flux correction for the convective term if required
      !
      !   $ res^{(0)} = res^{(0)} + f(u^n,u^n) $

      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'convectionAFC', convectionAFC)

      if (convectionAFC > 0) then

        ! What kind of stabilisation should be applied?
        select case(rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation)
          
        case (AFCSTAB_FEMFCT,&
              AFCSTAB_FEMFCT_CLASSICAL)

          call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                                   'imassantidiffusiontype', imassantidiffusiontype)
          
          ! Should we apply consistent mass antidiffusion?
          if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
            call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                       rsolution, rtimestep%theta, rtimestep%dStep, .true.,&
                                       rres, rproblemLevel%Rafcstab(convectionAFC),&
                                       rproblemLevel%Rmatrix(consistentMassMatrix))
          else
            call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                       rsolution, rtimestep%theta, rtimestep%dStep, .true.,&
                                       rres, rproblemLevel%Rafcstab(convectionAFC))
          end if
          
        case (AFCSTAB_FEMTVD)
          call gfsc_buildResidualTVD(rsolution, rtimestep%dStep, rres,&
                                     rproblemLevel%Rafcstab(convectionAFC))
          
        case (AFCSTAB_FEMGP)
          call gfsc_buildResidualGP(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                    rsolution, rsolutionInitial, rtimestep%theta, rtimestep%dStep,&
                                    rres, rproblemLevel%Rafcstab(convectionAFC))
        end select

      end if   ! convectionAFC > 0

      
      ! Perform algebraic flux correction for the diffusive term if required
      !
      !   $ res^{(0)} = res^{(0)} + g(u^n,u^n) $
      
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

    else   ! ite > 0

      !-------------------------------------------------------------------------
      ! In all subsequent nonlinear iterations only the residual vector is
      ! updated, using the right-hand side vector from the zero-th iteration
      !-------------------------------------------------------------------------

      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'imasstype', imasstype)

      select case(imasstype)
      case (MASS_LUMPED)
        
        ! Compute the low-order residual
        !
        !   $ res^{(m)} = rhs - [M_L - dt*theta*L(u^{(m)})]*u^{(m)} $

        call lsysbl_copyVector(rrhs, rres)

        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                                 rsolution%rvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 rtimestep%theta*rtimestep%dStep, 1.0_DP)

        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                 rsolution%RvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 -1.0_DP, 1.0_DP)
      case (MASS_CONSISTENT)
        
        ! Compute the low-order residual
        !
        !   $ res^{(m)} = rhs - [M_C - dt*theta*L(u^{(m)})]*u^{(m)} $

        call lsysbl_copyVector(rrhs, rres)

        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                                 rsolution%rvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 rtimestep%theta*rtimestep%dStep, 1.0_DP)

        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                 rsolution%RvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 -1.0_DP, 1.0_DP)
      case DEFAULT
        
        ! Compute the low-order residual
        !
        !   $ res^{(m)} = L(u^{(m)})*u^{(m)} $
        
        call lsyssc_scalarMatVec(rproblemLevel%Rmatrix(transportMatrix),&
                                 rsolution%rvectorBlock(1),&
                                 rres%RvectorBlock(1),&
                                 1.0_DP, 0.0_DP)
      end select


      ! Perform algebraic flux correction for the convective term if required
      !
      !   $ res = res + f^*(u^(m),u^n) $

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
            call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                       rsolution, rtimestep%theta, rtimestep%dStep, .false.,&
                                       rres, rproblemLevel%Rafcstab(convectionAFC),&
                                       rproblemLevel%Rmatrix(consistentMassMatrix))
          else
            call gfsc_buildResidualFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
                                       rsolution, rtimestep%theta, rtimestep%dStep, .false.,&
                                       rres, rproblemLevel%Rafcstab(convectionAFC))
          end if
          
        case (AFCSTAB_FEMTVD)
          call gfsc_buildResidualTVD(rsolution, rtimestep%dStep, rres,&
                                     rproblemLevel%Rafcstab(convectionAFC))
          
        case (AFCSTAB_FEMGP)
          call gfsc_buildResidualGP(rproblemLevel%Rmatrix(consistentMassMatrix),&
                                    rsolution, rsolutionInitial, rtimestep%theta, rtimestep%dStep,&
                                    rres, rproblemLevel%Rafcstab(convectionAFC))
        end select

      end if   ! convectionAFC > 0

      
      ! Perform algebraic flux correction for the diffusive term if required
      !
      !   $ res = res + g^*(u^n+1,u^n) $

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

    end if   ! ite = 0

    ! Apply the given load vector to the residual
    if (present(rb)) call lsysbl_vectorLinearComb(rb, rres, 1.0_DP, 1.0_DP)

    ! Stop time measurement for residual/rhs evaluation
    call stat_stopTimer(p_rtimer)
    
  end subroutine transp_calcResidual

  !*****************************************************************************

!<subroutine>

  subroutine transp_setBoundary(rproblemLevel, rtimestep, rsolver,&
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
    type(t_parlist), pointer :: p_rparlist
    integer :: imatrix

    
    ! Get parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
    
    select case(rsolver%iprecond)
    case (NLSOL_PRECOND_BLOCKD,&
          NLSOL_PRECOND_DEFCOR, &
          NLSOL_PRECOND_NEWTON_FAILED)

      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'systemmatrix', imatrix)
      
    case (NLSOL_PRECOND_NEWTON)

      call parlst_getvalue_int(p_rparlist, rcollection%SquickAccess(1),&
                               'jacobianmatrix', imatrix)
      
    case DEFAULT
      call output_line('Invalid nonlinear preconditioner!',&
                       OU_CLASS_ERROR, OU_MODE_STD,'transp_setBoundary')
      call sys_halt()
    end select
    
    ! Impose boundary conditions for the solution vector and impose
    ! zeros in the residual vector and the off-diagonal positions of
    ! the system matrix which is obtained from the collection
    
    call bdrf_filterSolution(rsolver%rboundaryCondition,& 
                             rproblemLevel%rtriangulation,&
                             rproblemLevel%Rmatrix(imatrix),&
                             rsolution, rres, rsolutionInitial,&
                             rtimestep%dTime)

  end subroutine transp_setBoundary

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcVelocityField(rparlist, ssectionName, rproblemLevel,&
                                      dtime, rcollection, nlminOpt)

!<description>
    ! This subroutine calculates the velocity fields from the function
    ! parser. The result is stored separately for each problem level.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(IN) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(IN) :: ssectionName

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
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_problemLevel), pointer :: p_rproblemLevel
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(NDIM3D+1) :: Dvalue
    character(LEN=SYS_STRLEN) :: svelocityname
    integer :: ieq, neq, idim, ndim, nlmin, icomp
    integer :: ivelocitytype, velocityfield

    
    ! Check if the velocity "vector" needs to be generated explicitly
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'ivelocitytype', ivelocitytype)
    if ((abs(ivelocitytype) .ne. VELOCITY_CONSTANT) .and.&
        (abs(ivelocitytype) .ne. VELOCITY_TIMEDEP)) return

    ! Get parameter from parameter list
    call parlst_getvalue_int(rparlist, ssectionName,&
                             'velocityfield', velocityfield)

    ! Get function parser from collection
    p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

    ! Set minimum problem level
    nlmin = rproblemLevel%ilev
    if (present(nlminOpt)) nlmin = nlminOpt

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
      call storage_getbase_double2d(&
          p_rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)


      ! Loop over all spatial dimensions
      do idim = 1, ndim

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
          Dvalue(1:ndim)   = p_DvertexCoords(:,ieq)
          Dvalue(NDIM3D+1) = dtime
          
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
    type(t_vectorBlock), intent(IN) :: rvector
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

  subroutine transp_calcLinearizedFCT(rbdrCond, rproblemLevel, rtimestep, rsolution, rcollection)

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

    ! local variables
    type(t_matrixScalar), pointer :: p_rmatrix
    type(t_parlist), pointer :: p_rparlist
    type(t_vectorScalar) :: rflux0, rflux
    type(t_vectorBlock) :: rdata
    real(DP), dimension(:), pointer :: p_MC, p_ML, p_Cx, p_Cy, p_u, p_flux0, p_flux, p_data
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
    call buildFlux2d(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ, nedge, p_u, &
                     rtimestep%dStep, p_MC, p_ML, p_Cx, p_Cy, p_data, p_flux, p_flux0)

    ! Build the correction
    call buildCorrection(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_rmatrix%NEQ,&
                         nedge, p_ML, p_flux, p_flux0, p_data, p_u)
    
    ! Set boundary conditions explicitly
    call bdrf_filterVectorExplicit(rbdrCond, rproblemLevel%rtriangulation,&
                                   rsolution, rtimestep%dTime,&
                                   rproblemLevel%p_rproblem%rboundary)

    ! Release flux vectors
    call storage_free(h_Ksep)
    call lsyssc_releaseVector(rflux0)
    call lsyssc_releaseVector(rflux)
    call lsysbl_releaseVector(rdata)

  contains

    !***************************************************************************

    subroutine buildFlux2d(Kld, Kcol, Kdiagonal, Ksep, NEQ, NEDGE, u,&
                           dscale, MC, ML, Cx, Cy, troc, flux0, flux)

      real(DP), dimension(:), intent(IN) :: MC,ML,Cx,Cy,u
      real(DP), intent(IN) :: dscale
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NEDGE
      
      integer, dimension(:), intent(INOUT) :: Ksep
      real(DP), dimension(:), intent(INOUT) :: flux0,flux
      
      real(DP), dimension(:), intent(OUT) :: troc     

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
        call transp_calcMatrixPrimalConst2d(u(i), u(i), C_ii, C_ii, i, i, k_ii, k_ii, d_ij)

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
          call transp_calcMatrixPrimalConst2d(u(i), u(j), C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
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
    
    subroutine buildCorrection(Kld, Kcol, Kdiagonal, Ksep, NEQ, NEDGE,&
                               ML, flux, flux0, data, u)

      
      real(DP), dimension(:), intent(IN) :: ML,flux0
      integer, dimension(:), intent(IN) :: Kld,Kcol,Kdiagonal
      integer, intent(IN) :: NEQ,NEDGE
      
      real(DP), dimension(:), intent(INOUT) :: data,u,flux
      integer, dimension(:), intent(INOUT) :: Ksep
      
      ! local variables
      real(DP), dimension(:), allocatable :: pp,pm,qp,qm,rp,rm
      real(DP) :: f_ij,diff
      integer :: ij,ji,i,j,iedge,ivar

      allocate(pp(neq), pm(neq), qp(neq), qm(neq), rp(neq), rm(neq))
      
      pp=0.0_DP; pm=0.0_DP
      qp=0.0_DP; qm=0.0_DP
      rp=0.0_DP; rm=0.0_DP
      
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
          f_ij = minmod(flux(iedge), flux0(iedge))
          
          ! ... and store prelimited flux
          flux(iedge) = f_ij

          diff = u(j)-u(i)
!!$          f_ij = flux(iedge)
!!$          if (f_ij * diff > 0) then
!!$            f_ij = 0.0_DP;  flux(iedge) = f_ij
!!$          end if

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
        if (pp(i) >  SYS_EPSREAL) rp(i) = min(1.0_DP, ML(i)*qp(i)/pp(i))
        if (pm(i) < -SYS_EPSREAL) rm(i) = min(1.0_DP, ML(i)*qm(i)/pm(i))
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

    end subroutine buildCorrection

    !***************************************************************************

    pure elemental function minmod(a,b)
      real(DP), intent(IN) :: a,b
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

  !*****************************************************************************

!<subroutine>

  subroutine transp_buildVectorScalarBdr(rdiscretisation, bclear,&
                                         rvectorScalar, rboundaryRegion, rcollection)

    use boundary

!<description>

    ! This subroutine builds the discretised vector for the linear
    ! form defined in terms of a boundary integral.

!</description>

!<input>

    ! The underlying discretisation structure which is to be used to
    ! create the vector.
    type(t_spatialDiscretisation), intent(IN), target :: rdiscretisation

    ! Whether to clear the vector before calculating the entries.
    ! If .FALSE., the new entries are added to the existing entries.
    logical, intent(IN) :: bclear

    ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
    ! to calculate. If not specified, the computation is done over
    ! the whole boundary.
    type(t_boundaryRegion), intent(IN), optional :: rboundaryRegion

    ! OPTIONAL: A collection structure. This structure is 
    ! given to the callback function for calculating the function
    ! which should be discretised in the linear form.
    type(t_collection), intent(INOUT), target, optional :: rcollection

!</input>

!<inputoutput>
    ! The FE vector. Calculated entries are imposed to this vector.
    type(t_vectorScalar), intent(INOUT) :: rvectorScalar
!</inputoutput>

!</subroutine>

    
    ! If the vector is not set up as new vector, it has to be unsorted.
    ! If it's a new vector, we switch off the sorting.
    if (bclear) then
      rvectorScalar%isortStrategy = -abs(rvectorScalar%isortStrategy)
    end if
    
    ! The vector must be unsorted, otherwise we can't set up the vector.
    if (rvectorScalar%isortStrategy .gt. 0) then
      call output_line('Vector must be unsorted!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_buildVectorScalarBdr')
      call sys_halt()
    end if

    select case(rvectorScalar%cdataType)
      
    case(ST_DOUBLE)
      call transp_buildVectorBdrDble_conf(rdiscretisation, bclear,&
                                          rvectorScalar, rboundaryRegion, rcollection)

    case DEFAULT
      call output_line('Single precision vectors currently not supported!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_buildVectorScalarBdr')
    end select

  end subroutine transp_buildVectorScalarBdr

!*****************************************************************************

!<subroutine>

  subroutine transp_buildVectorBdrDble_conf(rdiscretisation, bclear,&
                                            rvectorScalar, rboundaryRegion, rcollection)

    use dofmapping
    use boundary
    use triangulation

!<description>

    ! This subroutine builds the discretised vector for the linear
    ! form defined in terms of a boundary integral.
    !
    ! Double-precision version.

!</description>

!<input>

    ! The underlying discretisation structure which is to be used to
    ! create the vector.
    type(t_spatialDiscretisation), intent(IN), target :: rdiscretisation
    
    ! Whether to clear the vector before calculating the entries.
    ! If .FALSE., the new entries are added to the existing entries.
    logical, intent(IN) :: bclear

    ! OPTIONAL: A t_boundaryRegion specifying the boundary region where
    ! to calculate. If not specified, the computation is done over
    ! the whole boundary.
    type(t_boundaryRegion), intent(IN), optional :: rboundaryRegion

    ! OPTIONAL: A collection structure. This structure is 
    ! given to the callback function for calculating the function
    ! which should be discretised in the linear form.
    type(t_collection), intent(INOUT), target, optional :: rcollection

!</input>

!<inputoutput>
    ! The FE vector. Calculated entries are imposed to this vector.
    type(t_vectorScalar), intent(INOUT) :: rvectorScalar
!</inputoutput>

!</subroutine>
    
    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation

    ! The boundary structure - to shorten some things...
    type(t_boundary), pointer :: p_rboundary

    ! The boundary region used for computing the integral
    type(t_boundaryRegion) :: rboundaryReg

    ! local variables
    integer :: ibdc


    if (rvectorScalar%h_Ddata .eq. ST_NOHANDLE) then
      
      call lsyssc_createVecByDiscr(rdiscretisation, rvectorScalar, .true.)
      
    else
      
      if (bclear) call lsyssc_clearVector(rvectorScalar)
      
    end if

    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rdiscretisation%p_rtriangulation

    ! Get a pointer to the boundary - for easier access.
    p_rboundary => rdiscretisation%p_rboundary


    ! If the boundary region is specified, call transp_integral2D_conf
    ! for that boundary region. Otherwise, call transp_integral2D_conf
    ! for all possible boundary regions.
    if (present(rboundaryRegion)) then
      call transp_integral2D_conf(p_rboundary, p_rtriangulation,&
                                  rboundaryRegion, rvectorScalar, rcollection)
    else
      ! Create a boundary region for each boundary component and call
      ! the calculation routine for that.
      do ibdc = 1, boundary_igetNBoundComp(p_rboundary)
        call boundary_createRegion(p_rboundary, ibdc, 0, rboundaryReg)
        call transp_integral2D_conf(p_rboundary, p_rtriangulation,&
                                    rboundaryReg, rvectorScalar, rcollection)
      end do
    end if
    
  end subroutine transp_buildVectorBdrDble_conf

!*****************************************************************************

!<subroutine>

  subroutine transp_integral2D_conf(rboundary, rtriangulation,&
                                    rboundaryRegion, rvectorScalar, rcollection)

    use boundary
    use triangulation

!<description>

    ! This routine calculates the linear form 

!</description>

!<input>
    
    ! A boundary object that specifies the analytical boundary of the domain.
    type(t_boundary), intent(IN) :: rboundary
    
    ! A triangulation object that specifies the mesh in the domain.
    type(t_triangulation), intent(IN) :: rtriangulation
    
    ! A t_boundaryRegion specifying the boundary region where
    ! to calculate. 
    type(t_boundaryRegion), intent(IN) :: rboundaryRegion
    
    ! OPTIONAL: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(INOUT), optional :: rcollection

!</input>

!<inputoutput>
    ! The FE vector. Calculated entries are imposed to this vector.
    type(t_vectorScalar), intent(INOUT) :: rvectorScalar
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_fparser), pointer :: rfparser
    real(DP), dimension(:), pointer :: p_Ddata
    real(DP), dimension(NDIM3D+1, 4) :: DevalData
    real(DP), dimension(NDIM2D, 3) :: DvelocityField
    real(DP), dimension(3) :: DtargetFunc
    real(DP), dimension(NDIM2D) :: Dnormal
    real(DP) :: dpar1, dpar2, dpos1, dpos2, dvalue
    integer :: ibdc, NELbdc, ibdcoffset, iedge, iel, ilocaledge, i1, i2, nve, ipoint
    integer :: icompVelocity_X, icompVelocity_Y, icompTargetFuncBdrInt

    ! The triangulation structure - to shorten some things...
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer, dimension(:), pointer :: p_IedgesAtBoundary
    integer, dimension(:), pointer :: p_IelementsAtBoundary
    real(DP), dimension(:), pointer :: p_DedgeParameterValue
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Get some pointers and arrays for quicker access
    
    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx,&
        p_IboundaryCpIdx)
    call storage_getbase_int (rtriangulation%h_IedgesAtBoundary,&
        p_IedgesAtBoundary)
    call storage_getbase_int (rtriangulation%h_IelementsAtBoundary,&
        p_IelementsAtBoundary)
    call storage_getbase_int2d (rtriangulation%h_IedgesAtElement,&
        p_IedgesAtElement)
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,&
        p_IverticesAtElement)
    call storage_getbase_double (rtriangulation%h_DedgeParameterValue,&
        p_DedgeParameterValue)
    call storage_getbase_double (rtriangulation%h_DvertexParameterValue,&
        p_DvertexParameterValue)

    call lsyssc_getbase_double(rvectorScalar, p_Ddata)


    ! This subroutine assumes that the first two quick access string
    ! values hold the names of the function parsers in the collection.
    rfparser => collct_getvalue_pars(rcollection, trim(rcollection%SquickAccess(1)))

    icompVelocity_X = fparser_getFunctionNumber(rfparser, 'fp_velocity_x')
    icompVelocity_Y = fparser_getFunctionNumber(rfparser, 'fp_velocity_y')
    icompTargetFuncBdrInt = fparser_getFunctionNumber(rfparser, 'fp_targetFuncBdrInt')

    ! Boundary component?
    ibdc = rboundaryRegion%iboundCompIdx

    ! Number of elements on that boundary component?
    NELbdc = p_IboundaryCpIdx(ibdc+1)-p_IboundaryCpIdx(ibdc)
    
    ! Position of the boundary component?
    ibdcoffset = p_IboundaryCpIdx(ibdc)

    ! Clear auxiliary array for velocity field evaluation
    DevalData = 0.0_DP

    ! Set simulation time
    DevalData(NDIM3D+1,:) = rcollection%DquickAccess(1)


    ! Loop through the edges on the boundary component ibdc.
    do iedge = 1, NELbdc
      if (boundary_isInRegion(rboundaryRegion, ibdc,&
          p_DedgeParameterValue(iedge))) then
        
        ! Element orientation; i.e. the local number of the boundary edge 
        do ilocaledge = 1, ubound(p_IedgesAtElement,1)
          if (p_IedgesAtElement(ilocaledge, p_IelementsAtBoundary(iedge)) .eq. &
              p_IedgesAtBoundary(iedge)) exit
        end do
        
        ! Get global element number
        iel = p_IelementsAtBoundary(iedge)

        ! Get number of edges per element
        nve = tria_getNVE(p_IverticesAtElement, iel)

        ! Get global vertex numbers
        i1 = p_IverticesAtElement(ilocaledge, iel)
        i2 = p_IverticesAtElement(mod(ilocaledge, nve)+1, iel)
        
        ! Save the start parameter value of the edge -- in length
        ! parametrisation.
        dpar1 = p_DvertexParameterValue(iedge)
        
        ! Save the end parameter value. Be careful: The last edge
        ! must be treated differently!
        if (iedge .ne. NELbdc) then
          dpar2 = p_DvertexParameterValue(iedge+1)
        else
          dpar2 = boundary_dgetMaxParVal(rboundary,ibdc)
        end if

        ! Compute parameter positions
        dpos1 = boundary_convertParameter(rboundary, &
                ibdc, dpar1, rboundaryRegion%cparType, BDR_PAR_LENGTH)
            
        dpos2 = boundary_convertParameter(rboundary, &
                ibdc, dpar2, rboundaryRegion%cparType, BDR_PAR_LENGTH)

        ! Compute vertex positions
        call boundary_getCoords(rboundary, ibdc, dpar1,&
                                DevalData(1,1), DevalData(2,1))
        call boundary_getCoords(rboundary, ibdc, 0.5*(dpar1+dpar2),&
                                DevalData(1,2), DevalData(2,2))
        call boundary_getCoords(rboundary, ibdc, dpar2,&
                                DevalData(1,3), DevalData(2,3))
        
        ! Compute outward unit normal vector
        call boundary_getNormalVec2D(rboundary, ibdc, 0.5*(dpar1+dpar2),&
                                     Dnormal(1), Dnormal(2))

        ! Compute entries of velocity field at cubature points
        do ipoint = 1, 3
          call fparser_evalFunction(rfparser, icompVelocity_X, DevalData(:,ipoint), DvelocityField(1,ipoint))
          call fparser_evalFunction(rfparser, icompVelocity_Y, DevalData(:,ipoint), DvelocityField(2,ipoint))
        end do
        
        ! Compute entries of target functional at cubature points
        do ipoint = 1, 3
          call fparser_evalFunction(rfparser, icompTargetFuncBdrInt, DevalData(:,ipoint), DtargetFunc(ipoint))
        end do

        ! Compute the boundary integral for the first FE-function
        dvalue =   ( DtargetFunc(1) * sum(DvelocityField(:,1)*Dnormal) +&
                   2*DtargetFunc(2) * sum(DvelocityField(:,2)*Dnormal) ) *&
                   abs(dpos2-dpos1) / 6.0_DP

!!!        dvalue = (DtargetFunc(1) + 2*DtargetFunc(2)) * abs(dpos2-dpos1) / 6.0_DP

        p_Ddata(i1) = p_Ddata(i1) + dvalue

        ! Compute the boundary integral for the first FE-function
        dvalue =   ( DtargetFunc(3) * sum(DvelocityField(:,3)*Dnormal) +&
                   2*DtargetFunc(2) * sum(DvelocityField(:,2)*Dnormal) ) *&
                   abs(dpos2-dpos1) / 6.0_DP

!!!        dvalue = (DtargetFunc(3) + 2*DtargetFunc(2)) * abs(dpos2-dpos1) / 6.0_DP

        p_Ddata(i2) = p_Ddata(i2) + dvalue

      end if
    end do

  end subroutine transp_integral2D_conf

end module transport_callback

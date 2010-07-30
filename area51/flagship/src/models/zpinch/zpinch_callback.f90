!##############################################################################
!# ****************************************************************************
!# <name> zpinch_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all callback functions which are required to
!# solve the simplified MHD equations in arbitrary spatial dimensions.
!#
!# The following callback functions are available:
!#
!# 1.) zpinch_nlsolverCallback
!#     -> Callback routine for the nonlinear solver
!#
!# 2.) zpinch_calcVelocityField
!#     -> Calculates the velocity field for the transport model
!#
!# 3.) zpinch_calcLorentzforceTerm
!#     -> Calculates the Lorentz force term for given solutions
!#
!# 4.) zpinch_calcLinearisedFCT
!#     -> Calculates the linearised FCT correction
!#
!# 5.) zpinch_getVariable
!#     -> Extracts a single variable from the vector of conservative
!#        variables stores in interleave or block format
!#
!# </purpose>
!##############################################################################

module zpinch_callback

  use afcstabilisation
  use basicgeometry
  use boundaryfilter
  use collection
  use derivatives
  use euler_basic
  use euler_callback
  use euler_callback1d
  use euler_callback2d
  use euler_callback3d
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
  use scalarpde
  use solveraux
  use storage
  use timestepaux
  use transport_callback
  use transport_callback1d
  use transport_callback2d
  use transport_callback3d
  use trilinearformevaluation
  use zpinch_callback2d

  implicit none

  private
  public :: zpinch_nlsolverCallback
  public :: zpinch_calcVelocityField
  public :: zpinch_calcLorentzforceTerm
  public :: zpinch_calcLinearisedFCT
  public :: zpinch_checkPressure

contains

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_nlsolverCallback(rproblemLevel, rtimestep,&
      rsolver, rsolution, rsolution0, rrhs, rres, istep,&
      ioperationSpec, rcollection, istatus, rsource)

!<description>
    ! This subroutine is called by the nonlinear solver and it is
    ! responsible to assemble preconditioner, right-hand side vector,
    ! residual vector, etc. for the scalar transport model
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

    ! local parameter which is saved
    real(DP), save :: dtimeEuler = 0.0_DP
    real(DP), save :: dtimeTransport = 0.0_DP

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_vectorBlock), pointer :: p_rsolutionEuler, p_rsolutionTransport
    type(t_vectorBlock) :: rforce
    real(DP) :: dscaleLorentzForceTerm
    integer :: ilorentzForceType
    integer(i32) :: iSpec
    integer :: isystemFormat, jacobianMatrix

    if (trim(rsolver%ssolverName) .eq. 'NonlinearSolverEuler') then

      ! Set the first string quick access array to the section name
      ! of the Euler model which is stored in the third array
      rcollection%SquickAccess(1) = rcollection%SquickAccess(3)

      ! Do we have to calculate the preconditioner?
      ! --------------------------------------------------------------------------
      if (iand(ioperationSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) then

        ! Compute the preconditioner
        call euler_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
            rsolver, rsolution, rcollection)
      end if


      ! Do we have to calculate the residual?
      ! --------------------------------------------------------------------------
      if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

        if ((istep .eq. 0) .and.&
            (dtimeEuler .ne. rtimestep%dTime)) then

          ! Update time variable so that no re-evaluation of the
          ! const right-hand side takes place for current time step
          dtimeEuler = rtimestep%dTime

          ! Compute the constant right-hand side including the
          ! given explicit part of the Lorentz force term
          call euler_calcRhsThetaScheme(rproblemLevel, rtimestep,&
              rsolver, rsolution0, rrhs, rcollection, rsource)
        end if

        ! Set pointers to current solution vectors stored in the
        ! first and second quick access vector
        p_rsolutionEuler     => rcollection%p_rvectorQuickAccess1
        p_rsolutionTransport => rcollection%p_rvectorQuickAccess2

        ! Set pointer to parameter list
        p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

        ! Get parameters from parameter list
        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(2), 'ilorentzforcetype', ilorentzForceType)
        
        ! Calculate scaling for implicit part of the Lorentz force
        dscaleLorentzForceTerm = -rtimestep%theta * rtimestep%dStep

        if ((ilorentzForceType .ne. 0) .and.&
            (dscaleLorentzForceTerm .ne. 0.0_DP)) then

          ! Compute the implicit part of the Lorentz force
          call zpinch_calcLorentzforceTerm(p_rparlist,&
              rcollection%SquickAccess(2), rcollection%SquickAccess(3),&
              rcollection%SquickAccess(4), rproblemLevel,&
              p_rsolutionEuler, p_rsolutionTransport, rtimestep%dTime,&
              dscaleLorentzForceTerm, rforce, rcollection)

          ! Compute the residual including the pre-computed implicit
          ! part of the Lorentz force term
          call euler_calcResidualThetaScheme(rproblemLevel, rtimestep,&
              rsolver, rsolution, rsolution0, rrhs, rres, istep,&
              rcollection, rforce)

          ! Release temporal memory
          call lsysbl_releaseVector(rforce)

        else

          ! Compute the residual without the Lorentz force term
          call euler_calcResidualThetaScheme(rproblemLevel, rtimestep,&
              rsolver, rsolution, rsolution0, rrhs, rres, istep, rcollection)

        end if
      end if


      ! Do we have to impose boundary conditions?
      ! --------------------------------------------------------------------------
      if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

        ! Impose boundary conditions
        call euler_setBoundaryConditions(rproblemLevel, rtimestep,&
            rsolver, rsolution, rsolution0, rres, rcollection)
      end if


    elseif (trim(rsolver%ssolverName) .eq. 'NonlinearSolverTransport') then

      ! Set the first string quick access array to the section name
      ! of the transport model which is stored in the fourth array
      rcollection%SquickAccess(1) = rcollection%SquickAccess(4)

      ! Set pointer to parameter list
      p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')
      
      ! Get configuration from Euler section of parameter list 
      call parlst_getvalue_int(p_rparlist,&
          rcollection%SquickAccess(3), 'isystemformat', isystemFormat)
      
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

      ! Do we have to calculate the constant right-hand side?
      ! --------------------------------------------------------------------------
      if ((iand(iSpec, NLSOL_OPSPEC_CALCRHS)  .ne. 0)) then

        ! What type of system format are we?
        select case(isystemFormat)
          
        case (SYSTEM_INTERLEAVEFORMAT)
          
          ! Compute the preconditioner in interleaved format
          call transp_calcPrecondThetaScheme(rproblemLevel,&
              rtimestep, rsolver, rsolution, rcollection,&
              zpinch_calcMatDiagConvIntlP2d_sim,&
              zpinch_calcMatRusConvIntlP2d_sim,&
              zpinch_calcMatDiagConvIntlD2d_sim,&
              zpinch_calcMatRusConvIntlD2d_sim,&
              transp_coeffMatBdrConvectionP2d,&
              transp_coeffMatBdrConvectionD2d)

        case (SYSTEM_BLOCKFORMAT)

          ! Compute the preconditioner in block format
          call transp_calcPrecondThetaScheme(rproblemLevel,&
              rtimestep, rsolver, rsolution, rcollection,&
              zpinch_calcMatDiagConvBlockP2d_sim,&
              zpinch_calcMatRusConvBlockP2d_sim,&
              zpinch_calcMatDiagConvBlockD2d_sim,&
              zpinch_calcMatRusConvBlockD2d_sim,&
              transp_coeffMatBdrConvectionP2d,&
              transp_coeffMatBdrConvectionD2d)

        case DEFAULT
          call output_line('Invalid system format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'zpinch_nlsolverCallback')
          call sys_halt()
        end select
        
        ! Compute the right-hand side
        call transp_calcRhsRungeKuttaScheme(rproblemLevel,&
            rtimestep, rsolver, rsolution, rsolution0,&
            rrhs, istep, rcollection,&
            fcb_coeffVecBdrPrimal_sim=transp_coeffVecBdrConvectionP2d,&
            fcb_coeffVecBdrDual_sim=transp_coeffVecBdrConvectionD2d)
        
        ! Remove specifier for the preconditioner (if any)
        iSpec = iand(iSpec, not(NLSOL_OPSPEC_CALCPRECOND))
      end if


      ! Do we have to calculate the residual?
      ! --------------------------------------------------------------------------
      if (iand(iSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

        if (istep .eq. 0) then
          
          if (dtimeTransport .ne. rtimestep%dTime) then

            ! Update time variable so that no re-evaluation of the
            ! const right-hand side takes place for current time step
            dtimeTransport = rtimestep%dTime

            ! What type of system format are we?
            select case(isystemFormat)
              
            case (SYSTEM_INTERLEAVEFORMAT)
              
              ! Compute the preconditioner in interleaved format
              call transp_calcPrecondThetaScheme(rproblemLevel,&
                  rtimestep, rsolver, rsolution0, rcollection,&
                  zpinch_calcMatDiagConvIntlP2d_sim,&
                  zpinch_calcMatRusConvIntlP2d_sim,&
                  zpinch_calcMatDiagConvIntlD2d_sim,&
                  zpinch_calcMatRusConvIntlD2d_sim,&
                  transp_coeffMatBdrConvectionP2d,&
                  transp_coeffMatBdrConvectionD2d)

            case (SYSTEM_BLOCKFORMAT)

              ! Compute the preconditioner in block format
              call transp_calcPrecondThetaScheme(rproblemLevel,&
                  rtimestep, rsolver, rsolution0, rcollection,&
                  zpinch_calcMatDiagConvBlockP2d_sim,&
                  zpinch_calcMatRusConvBlockP2d_sim,&
                  zpinch_calcMatDiagConvBlockD2d_sim,&
                  zpinch_calcMatRusConvBlockD2d_sim,&
                  transp_coeffMatBdrConvectionP2d,&
                  transp_coeffMatBdrConvectionD2d)

            case DEFAULT
              call output_line('Invalid system format!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'zpinch_nlsolverCallback')
              call sys_halt()
            end select
            
            ! Assemble the constant right-hand side
            call transp_calcRhsThetaScheme(rproblemLevel, rtimestep,&
                rsolver, rsolution0, rrhs, rcollection,&
                fcb_coeffVecBdrPrimal_sim=transp_coeffVecBdrConvectionP2d,&
                fcb_coeffVecBdrDual_sim=transp_coeffVecBdrConvectionD2d)

          end if

          ! Set pointer to the solution vector of the Euler model
          ! from the current time step
          p_rsolutionEuler => rcollection%p_rvectorQuickAccess1

          ! Calculate the velocity vector using the solution of the
          ! Euler model from the current time step
          call zpinch_calcVelocityField(p_rparlist,&
              rcollection%SquickAccess(1), rproblemLevel,&
              p_rsolutionEuler, rcollection)

          ! What type of system format are we?
          select case(isystemFormat)
            
          case (SYSTEM_INTERLEAVEFORMAT)
            
            ! Compute the preconditioner in interleaved format
            call transp_calcPrecondThetaScheme(rproblemLevel,&
                rtimestep, rsolver, rsolution, rcollection,&
                zpinch_calcMatDiagConvIntlP2d_sim,&
                zpinch_calcMatRusConvIntlP2d_sim,&
                zpinch_calcMatDiagConvIntlD2d_sim,&
                zpinch_calcMatRusConvIntlD2d_sim,&
                transp_coeffMatBdrConvectionP2d,&
                transp_coeffMatBdrConvectionD2d)

          case (SYSTEM_BLOCKFORMAT)
            
            ! Compute the preconditioner in block format
            call transp_calcPrecondThetaScheme(rproblemLevel,&
                rtimestep, rsolver, rsolution, rcollection,&
                zpinch_calcMatDiagConvBlockP2d_sim,&
                zpinch_calcMatRusConvBlockP2d_sim,&
                zpinch_calcMatDiagConvBlockD2d_sim,&
                zpinch_calcMatRusConvBlockD2d_sim,&
                transp_coeffMatBdrConvectionP2d,&
                transp_coeffMatBdrConvectionD2d)

          case DEFAULT
            call output_line('Invalid system format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'zpinch_nlsolverCallback')
            call sys_halt()
          end select
          
          ! Remove specifier for the preconditioner (if any)
          iSpec = iand(iSpec, not(NLSOL_OPSPEC_CALCPRECOND))
        end if
        
        ! Compute the residual
        call transp_calcResidualThetaScheme(rproblemLevel,&
            rtimestep, rsolver, rsolution, rsolution0,&
            rrhs, rres, istep, rcollection)
      end if


      ! Do we have to calculate the preconditioner?
      ! --------------------------------------------------------------------------
      if (iand(iSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) then
        
        ! What type of system format are we?
        select case(isystemFormat)
          
        case (SYSTEM_INTERLEAVEFORMAT)
          
          ! Compute the preconditioner in interleaved format
          call transp_calcPrecondThetaScheme(rproblemLevel,&
              rtimestep, rsolver, rsolution, rcollection,&
              zpinch_calcMatDiagConvIntlP2d_sim,&
              zpinch_calcMatRusConvIntlP2d_sim,&
              zpinch_calcMatDiagConvIntlD2d_sim,&
              zpinch_calcMatRusConvIntlD2d_sim,&
              transp_coeffMatBdrConvectionP2d,&
              transp_coeffMatBdrConvectionD2d)

        case (SYSTEM_BLOCKFORMAT)
          
          ! Compute the preconditioner in block format
          call transp_calcPrecondThetaScheme(rproblemLevel,&
              rtimestep, rsolver, rsolution, rcollection,&
              zpinch_calcMatDiagConvBlockP2d_sim,&
              zpinch_calcMatRusConvBlockP2d_sim,&
              zpinch_calcMatDiagConvBlockD2d_sim,&
              zpinch_calcMatRusConvBlockD2d_sim,&
              transp_coeffMatBdrConvectionP2d,&
              transp_coeffMatBdrConvectionD2d)

        case DEFAULT
          call output_line('Invalid system format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'zpinch_nlsolverCallback')
          call sys_halt()
        end select
      end if


      ! Do we have to calculate the Jacobian operator?
      ! --------------------------------------------------------------------------
      if (iand(iSpec, NLSOL_OPSPEC_CALCJACOBIAN) .ne. 0) then

        ! What type of system format are we?
        select case(isystemFormat)
          
        case (SYSTEM_INTERLEAVEFORMAT)
          
          ! Compute the Jacobian matrix in interleaved format
          call transp_calcJacobianThetaScheme(rproblemLevel,&
              rtimestep, rsolver, rsolution, rsolution0, rcollection,&
              zpinch_calcMatRusConvIntlP2d_sim,&
              zpinch_calcMatRusConvIntlD2d_sim,&
              transp_coeffMatBdrConvectionP2d,&
              transp_coeffMatBdrConvectionD2d)

        case (SYSTEM_BLOCKFORMAT)

          ! Compute the Jacobian matrix in block format
          call transp_calcJacobianThetaScheme(rproblemLevel,&
              rtimestep, rsolver, rsolution, rsolution0, rcollection,&
              zpinch_calcMatRusConvBlockP2d_sim,&
              zpinch_calcMatRusConvBlockD2d_sim,&
              transp_coeffMatBdrConvectionP2d,&
              transp_coeffMatBdrConvectionD2d)

        case DEFAULT
          call output_line('Invalid system format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'zpinch_nlsolverCallback')
          call sys_halt()
        end select
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

        ! Set pointers to parameter list
        p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

        ! Get position of Jacobian matrix
        call parlst_getvalue_int(p_rparlist,&
            rcollection%SquickAccess(1), 'jacobianMatrix', jacobianMatrix)

        ! Apply Jacobian matrix
        call lsyssc_scalarMatVec(&
            rproblemLevel%Rmatrix(jacobianMatrix),&
            rsolution%RvectorBlock(1), rres%RvectorBlock(1),&
            1.0_DP, 1.0_DP)
      end if
      
    else

      call output_line('Invalid nonlinear solver!',&
          OU_CLASS_ERROR,OU_MODE_STD,'zpinch_nlsolverCallback')
      call sys_halt()
      
    end if

    ! Set status flag
    istatus = 0

  end subroutine zpinch_nlsolverCallback

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcVelocityField(rparlist, ssectionName,&
      rproblemLevel, rsolution, rcollection)

!<description>
    ! This subroutine initializes the velocity field from the solution
    ! of the compressible Euler model. The result is stored separately
    ! for each problem level.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! solution vector of compressible Euler model
    type(t_vectorBlock), intent(in) :: rsolution
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout), target :: rproblemLevel

    ! collection
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: velocityfield
    integer :: neq, ndim

    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'velocityfield', velocityfield)

    ! Get number of degrees of freedom and spatial dimension
    neq  = rproblemLevel%rtriangulation%NVT
    ndim = rproblemLevel%rtriangulation%ndim

    ! Create/resize velocity vector if required
    if (rproblemLevel%RvectorBlock(velocityfield)%NEQ .eq. 0) then
      call lsysbl_createVectorBlock(&
          rproblemLevel%rvectorBlock(velocityfield),&
          neq, ndim, .true.)
    elseif (rproblemLevel%RvectorBlock(velocityfield)%NEQ&
            .ne. neq*ndim) then
      call lsysbl_resizeVectorBlock(&
          rproblemLevel%rvectorBlock(velocityfield),&
          neq, .true.)
    end if

    ! Set x-velocity
    call euler_getVariable(rsolution, 'velocity_x',&
        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(1))
    call zpinch_setVariable2d(&
        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(1), 1)

    ! Set y-velocity
    call euler_getVariable(rsolution, 'velocity_y',&
        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(2))
    call zpinch_setVariable2d(&
        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(2), 2)

    ! Set the global solution vector at the current time step
    call lsysbl_duplicateVector(rsolution,&
        rproblemLevel%RvectorBlock(2),&
        LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_COPY)
    call zpinch_setVariable2d(rproblemLevel%RvectorBlock(2), 3)

    ! Set update notification in problem level structure
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     PROBLEV_MSPEC_UPDATE)

  end subroutine zpinch_calcVelocityField

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcLorentzforceTerm(rparlist, ssectionName,&
      ssectionNameEuler, ssectionNameTransport, rproblemLevel,&
      rsolutionEuler, rsolutionTransport, dtime, dscale, rforce,&
      rcollection)

!<description>
    ! This subroutine evaluates the Lorentz force term based on the
    ! given solution vectors and applies stores it in vector rsource.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section names in parameter list
    character(LEN=*), intent(in) :: ssectionName
    character(LEN=*), intent(in) :: ssectionNameEuler
    character(LEN=*), intent(in) :: ssectionNameTransport

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! solution vector for transport model
    type(t_vectorBlock), intent(in) :: rsolutionTransport

    ! solution vector for Euler model
    type(t_vectorBlock), intent(in) :: rsolutionEuler

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling parameter
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! source vector to be assembled
    type(t_vectorBlock), intent(inout) :: rforce

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DdataTransport, p_DdataEuler, p_Ddata
    real(DP), dimension(:), pointer :: p_DlumpedMassMatrix
    character(LEN=SYS_STRLEN) :: slorentzforceName
    real(DP) :: dcurrentDrive, deffectiveRadius
    integer :: isystemFormat, lumpedMassMatrix
    integer :: neq, nvar, icomp, icoordinatesystem
    logical :: bcompatible


    ! Check if solution vector and Lorentz force vector are compatible
    call lsysbl_isVectorCompatible(rsolutionEuler, rforce, bcompatible)
    if (.not.bcompatible)&
        call lsysbl_resizeVectorBlock(rforce, rsolutionEuler, .false.)

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName,&
        'slorentzforcename', slorentzforceName)
    call parlst_getvalue_double(rparlist, ssectionName,&
        'deffectiveradius', deffectiveRadius)
    call parlst_getvalue_int(rparlist, ssectionName,&
        'icoordinatesystem', icoordinatesystem)
    call parlst_getvalue_int(rparlist, ssectionNameEuler,&
        'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(rparlist, ssectionNameEuler,&
        'isystemformat', isystemFormat)
    
    ! Get lumped and consistent mass matrix
    if (lumpedMassMatrix .gt. 0) then
      call lsyssc_getbase_double(&
          rproblemLevel%Rmatrix(lumpedMassMatrix), p_DlumpedMassMatrix)
    else
      call output_line('Lumped mass matrix is not available!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'zpinch_calcLorentzforceTerm')
      call sys_halt()
    end if
    
    ! Set pointer to global solution vectors
    call lsysbl_getbase_double(rforce, p_Ddata)
    call lsysbl_getbase_double(rsolutionEuler, p_DdataEuler)
    call lsysbl_getbase_double(rsolutionTransport, p_DdataTransport)

    ! Set pointer to the vertex coordinates
    call storage_getbase_double2D(&
        rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! Set dimensions
    neq  = rsolutionTransport%NEQ
    nvar = euler_getNVAR(rproblemLevel)

    ! Get function parser from collection structure
    p_rfparser => collct_getvalue_pars(rcollection, 'rfparser')

    ! Get the number of the component used for
    ! evaluating the Lorentz force term
    icomp = fparser_getFunctionNumber(p_rfparser, slorentzforceName)

    ! Evaluate the function parser
    call fparser_evalFunction(p_rfparser, icomp, (/dtime/), dcurrentDrive)

    ! Multiply scaling parameter by the time step
    dcurrentDrive = dscale * dcurrentDrive

    ! What type of system format are we?
    select case(isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)

      ! What type of coordinate system are we?
      select case(icoordinatesystem)
      case(1)
        call calcForceXYInterleaveFormat(dcurrentDrive, deffectiveRadius,&
            neq, nvar, p_DvertexCoords, p_DlumpedMassMatrix,&
            p_DdataTransport, p_DdataEuler, p_Ddata)

      case(2)
        call calcForceRZInterleaveFormat(dcurrentDrive, deffectiveRadius,&
            neq, nvar, p_DvertexCoords, p_DlumpedMassMatrix,&
            p_DdataTransport, p_DdataEuler, p_Ddata)

      case default
        call output_line('Invalid type of coordinate system!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcLorentzforceTerm')
        call sys_halt()
      end select

    case (SYSTEM_BLOCKFORMAT)
      ! What type of coordinate system are we?
      select case(icoordinatesystem)
      case(1)
        call calcForceXYBlockFormat(dcurrentDrive, deffectiveRadius,&
            neq, nvar, p_DvertexCoords, p_DlumpedMassMatrix,&
            p_DdataTransport, p_DdataEuler, p_Ddata)

      case(2)
        call calcForceRZBlockFormat(dcurrentDrive, deffectiveRadius,&
            neq, nvar, p_DvertexCoords, p_DlumpedMassMatrix,&
            p_DdataTransport, p_DdataEuler, p_Ddata)

      case default
        call output_line('Invalid type of coordinate system!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcLorentzforceTerm')
        call sys_halt()
      end select

    case DEFAULT
      call output_line('Invalid system format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcLorentzforceTerm')
      call sys_halt()
    end select


  contains

    ! Here, the real working routines follow

    !**************************************************************
    ! Calculate the Lorentz force term in x-y coordinates.
    ! The system is stored in interleave format.

    subroutine calcForceXYInterleaveFormat(dcurrentDrive,&
        deffectiveRadius, neq, nvar, DvertexCoords,&
        DmassMatrix, DdataTransport, DdataEuler, DdataForce)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(nvar,neq), intent(in) :: DdataEuler
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      real(DP), intent(in) :: dcurrentDrive, deffectiveRadius
      integer, intent(in) :: neq, nvar

      real(DP), dimension(nvar,neq), intent(out) :: DdataForce

      ! local variables
      real(DP) :: drad, daux, x1, x2
      integer :: i

      ! Loop over all rows
      !$omp parallel do private(x1,x2,drad,daux)
      do i = 1, neq

        ! Get coordinates at node i
        x1 = DvertexCoords(1, i)
        x2 = DvertexCoords(2, i)
        
        ! Compute unit vector into origin
        drad = sqrt(x1*x1 + x2*x2)
        if (drad .gt. SYS_EPSREAL) then
          x1 = x1/drad
          x2 = x2/drad
        else
          x1 = 0.0_DP
          x2 = 0.0_DP
        end if
        
        ! Compute Lorentz source term
        daux = dcurrentDrive * DmassMatrix(i)*DdataTransport(i) /&
               max(drad, deffectiveRadius)

        ! Impose source values
        DdataForce(1,i) = 0.0_DP
        DdataForce(2,i) = daux * x1
        DdataForce(3,i) = daux * x2
        DdataForce(4,i) = daux * (DdataEuler(2,i)*x1 +&
                                  DdataEuler(3,i)*x2) / DdataEuler(1,i)
      end do
      !$omp end parallel do

    end subroutine calcForceXYInterleaveFormat


    !**************************************************************
    ! Calculate the Lorentz force term in r-z coordinates.
    ! The system is stored in interleave format.

    subroutine calcForceRZInterleaveFormat(dcurrentDrive,&
        deffectiveRadius, neq, nvar, DvertexCoords, DmassMatrix,&
        DdataTransport, DdataEuler, DdataForce)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(nvar,neq), intent(in) :: DdataEuler
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      real(DP), intent(in) :: dcurrentDrive, deffectiveRadius
      integer, intent(in) :: neq, nvar
      
      real(DP), dimension(nvar,neq), intent(out) :: DdataForce
      
      ! local variables
      real(DP) :: drad, daux, x1
      integer :: i

      ! Loop over all rows
      !$omp parallel do private(x1,drad,daux)
      do i = 1, neq

        ! Get x-coordinate at node i
        x1 = DvertexCoords(1, i); drad = x1
        
        ! Compute unit vector into origin
        if (drad .gt. SYS_EPSREAL) then
          x1 = 1.0_DP
        else
          x1 = 0.0_DP
        end if

        ! Compute Lorentz source term
        daux = dcurrentDrive * DmassMatrix(i)*DdataTransport(i) /&
               max(drad, deffectiveRadius)
        
        ! Impose source values
        DdataForce(1,i) = 0.0_DP
        DdataForce(2,i) = daux * x1
        DdataForce(3,i) = 0.0_DP
        DdataForce(4,i) = daux * DdataEuler(2,i)*x1 / DdataEuler(1,i)
      end do
      !$omp end parallel do

    end subroutine calcForceRZInterleaveFormat


    !**************************************************************
    ! Calculate the Lorentz force term in x-y coordinates.
    ! The system is stored in block format.

    subroutine calcForceXYBlockFormat(dcurrentDrive,&
        deffectiveRadius, neq, nvar, DvertexCoords,&
        DmassMatrix, DdataTransport, DdataEuler, DdataForce)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(neq,nvar), intent(in) :: DdataEuler
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      real(DP), intent(in) :: dcurrentDrive, deffectiveRadius
      integer, intent(in) :: neq, nvar

      real(DP), dimension(neq,nvar), intent(out) :: DdataForce

      ! local variables
      real(DP) :: drad, daux, x1, x2
      integer :: i

      ! Loop over all rows
      !$omp parallel do private(x1,x2,drad,daux)
      do i = 1, neq

        ! Get coordinates at node i
        x1 = DvertexCoords(1, i)
        x2 = DvertexCoords(2, i)
        
        ! Compute unit vector into origin
        drad = sqrt(x1*x1 + x2*x2)
        if (drad .gt. SYS_EPSREAL) then
          x1 = x1/drad
          x2 = x2/drad
        else
          x1 = 0.0_DP
          x2 = 0.0_DP
        end if
        
        ! Compute Lorentz source term
        daux = dcurrentDrive * DmassMatrix(i)*DdataTransport(i) /&
               max(drad, deffectiveRadius)

        ! Impose source values
        DdataForce(i,1) = 0.0_DP
        DdataForce(i,2) = daux * x1
        DdataForce(i,3) = daux * x2
        DdataForce(i,4) = daux * (DdataEuler(i,2)*x1 +&
                                  DdataEuler(i,3)*x2) / DdataEuler(i,1)
      end do
      !$omp end parallel do

    end subroutine calcForceXYBlockFormat


    !**************************************************************
    ! Calculate the Lorentz force term in r-z coordinates.
    ! The system is stored in block format.

    subroutine calcForceRZBlockFormat(dcurrentDrive,&
        deffectiveRadius, neq, nvar, DvertexCoords, DmassMatrix,&
        DdataTransport, DdataEuler, DdataForce)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(neq,nvar), intent(in) :: DdataEuler
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      real(DP), intent(in) :: dcurrentDrive, deffectiveRadius
      integer, intent(in) :: neq, nvar
      
      real(DP), dimension(neq,nvar), intent(out) :: DdataForce
      
      ! local variables
      real(DP) :: drad, daux, x1
      integer :: i

      ! Loop over all rows
      !$omp parallel do private(x1,drad,daux)
      do i = 1, neq

        ! Get x-coordinate at node i
        x1 = DvertexCoords(1, i); drad = x1
        
        ! Compute unit vector into origin
        if (drad .gt. SYS_EPSREAL) then
          x1 = 1.0_DP
        else
          x1 = 0.0_DP
        end if

        ! Compute Lorentz source term
        daux = dcurrentDrive * DmassMatrix(i)*DdataTransport(i) /&
               max(drad, deffectiveRadius)
        
        ! Impose source values
        DdataForce(i,1) = 0.0_DP
        DdataForce(i,2) = daux * x1
        DdataForce(i,3) = 0.0_DP
        DdataForce(i,4) = daux * DdataEuler(i,2)*x1 / DdataEuler(i,1)
      end do
      !$omp end parallel do

    end subroutine calcForceRZBlockFormat

  end subroutine zpinch_calcLorentzforceTerm

  !*****************************************************************************

!<subroutine>
  
  subroutine zpinch_calcLinearisedFCT(RbdrCond, rproblemLevel, rtimestep,&
      rsolverEuler, rsolverTransport, Rsolution, rcollection, Rsource)

!<description>
    ! This subroutine calculates the linearised FCT correction
!</description>

!<input>
    ! array of boundary condition structure
    type(t_boundaryCondition), dimension(:), intent(in) :: RbdrCond

    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! OPTIONAL: source vector
    type(t_vectorBlock), dimension(:), intent(in), optional :: Rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structures
    type(t_solver), intent(inout) :: rsolverEuler
    type(t_solver), intent(inout) :: rsolverTransport

    ! solution vectors
    type(t_vectorBlock), dimension(:), intent(inout) :: Rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_timestep) :: rtimestepAux
    type(t_vectorBlock), dimension(:), pointer :: p_Rpredictor
    type(t_vectorBlock), pointer :: p_rpredictorEuler,p_rpredictorTransport
    type(t_parlist), pointer :: p_rparlist
    character(len=SYS_STRLEN), dimension(:), pointer :: SfailsafeVariables
    real(DP), dimension(:), pointer :: p_DalphaEuler, p_DalphaTransport
    integer, dimension(2) :: IposAFC
    integer :: convectionAFC,inviscidAFC,lumpedMassMatrix,consistentMassMatrix
    integer :: imassantidiffusiontype,nfailsafe,ivariable,nvariable
    integer :: ilimitersynchronisation
    
    ! Set pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(3), 'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(4), 'convectionAFC', convectionAFC)

    ! Do we have to apply linearised FEM-FCT?
    if ((inviscidAFC .le. 0) .or. (convectionAFC .le. 0)) return
    if ((rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation&
         .ne. AFCSTAB_FEMFCT_LINEARISED) .or.&
        (rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation&
         .ne. AFCSTAB_FEMFCT_LINEARISED)) return

    ! Set positions of stabilisation structures
    IposAFC=(/inviscidAFC, convectionAFC/)

    ! Get more parameters from parameter list
     call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(2),&
        'ilimitersynchronisation', ilimitersynchronisation)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(3),&
        'lumpedmassmatrix', lumpedmassmatrix)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(3),&
        'consistentmassmatrix', consistentmassmatrix)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(3),&
        'nfailsafe', nfailsafe)
    call parlst_getvalue_int(p_rparlist,&
        rcollection%SquickAccess(3),&
        'imassantidiffusiontype', imassantidiffusiontype)
    
    !---------------------------------------------------------------------------
    ! Linearised FEM-FCT algorithm
    !---------------------------------------------------------------------------

    ! Initialize dummy timestep
    rtimestepAux%dStep = 1.0_DP
    rtimestepAux%theta = 0.0_DP
    
    !--- compressible Euler model ----------------------------------------------

    ! Set pointer to predictor
    p_rpredictorEuler => rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorPredictor

    ! Set the first string quick access array to the section name
    ! of the Euler model which is stored in the third array
    rcollection%SquickAccess(1) = rcollection%SquickAccess(3)

    ! Compute low-order "right-hand side" without theta parameter
    call euler_calcRhsThetaScheme(rproblemLevel, rtimestepAux,&
        rsolverEuler, Rsolution(1), p_rpredictorEuler, rcollection, Rsource(1))

    ! Compute low-order predictor
    call lsysbl_invertedDiagMatVec(&
        rproblemLevel%Rmatrix(lumpedMassMatrix),&
        p_rpredictorEuler, 1.0_DP, p_rpredictorEuler)

    ! Compute the raw antidiffusive fluxes
    call euler_calcFluxFCT(rproblemLevel, p_rpredictorEuler,&
        Rsolution(1), 0.0_DP, 1.0_DP, 1.0_DP, .true., rcollection)

    !--- transport model -------------------------------------------------------

    ! Set pointer to predictor
    p_rpredictorTransport => rproblemLevel%Rafcstab(convectionAFC)%p_rvectorPredictor

    ! Set the first string quick access array to the section name
    ! of the transport model which is stored in the fourth array
    rcollection%SquickAccess(1) = rcollection%SquickAccess(4)

    ! Compute the preconditioner
    call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
        rsolverTransport, Rsolution(2), rcollection)

    ! Compute low-order "right-hand side" without theta parameter
    call transp_calcRhsThetaScheme(rproblemLevel, rtimestepAux,&
        rsolverTransport, Rsolution(2), p_rpredictorTransport,&
        rcollection,&
        fcb_coeffVecBdrPrimal_sim = transp_coeffVecBdrConvectionP2d,&
        fcb_coeffVecBdrDual_sim = transp_coeffVecBdrConvectionD2d)

    ! Compute low-order predictor
    call lsysbl_invertedDiagMatVec(&
        rproblemLevel%Rmatrix(lumpedMassMatrix),&
        p_rpredictorTransport, 1.0_DP, p_rpredictorTransport)

    ! Should we apply consistent mass antidiffusion?
    if (imassantidiffusiontype .eq. MASS_CONSISTENT) then
      call gfsc_buildFluxFCT(rproblemLevel%Rafcstab(convectionAFC),&
          p_rpredictorTransport, Rsolution(2),&
          rtimestepAux%theta, rtimestepAux%dStep, 1.0_DP, .true.,&
          rproblemLevel%Rmatrix(consistentMassMatrix))
    else
      call gfsc_buildFluxFCT(rproblemLevel%Rafcstab(convectionAFC),&
          p_rpredictorTransport, Rsolution(2),&
          rtimestepAux%theta, rtimestepAux%dStep, 1.0_DP, .true.)
    end if

    !---------------------------------------------------------------------------
    ! Perform failsafe flux correction (if required)
    !---------------------------------------------------------------------------
    
    if (nfailsafe .gt. 0) then
      
      ! Get number of failsafe variables
      nvariable = max(1,&
          parlst_querysubstrings(p_rparlist,&
          rcollection%SquickAccess(3), 'sfailsafevariable'))
      
      ! Allocate character array that stores all failsafe variable names
      allocate(SfailsafeVariables(nvariable))
      
      ! Initialize character array with failsafe variable names
      do ivariable = 1, nvariable
        call parlst_getvalue_string(p_rparlist,&
            rcollection%SquickAccess(3), 'sfailsafevariable',&
            Sfailsafevariables(ivariable), isubstring=ivariable)
      end do

      ! Set up the predictor, that is, make a virtual copy of the
      ! predictor vectors from the Euler system and the scalar tracer
      ! equation which can be passed to the failsafe subroutine.
      allocate(p_rpredictor(2))
      call lsysbl_duplicateVector(p_rpredictorEuler, p_Rpredictor(1),&
          LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_SHARE)
      call lsysbl_duplicateVector(p_rpredictorTransport, p_Rpredictor(2),&
          LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_SHARE)
    end if

    ! Set the first string quick access array to the section name
    ! of the Euler model which is stored in the third array
    rcollection%SquickAccess(1) = rcollection%SquickAccess(3)
    
    ! What type of limiter synchronisation is performed?
    select case(ilimitersynchronisation)
      
    case (0)   ! no synchronisation
      
      if (nfailsafe .gt. 0) then
        
        ! Compute linearised FEM-FCT correction
        call euler_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
            AFCSTAB_FCTALGO_CORRECT, Rsolution(1), rcollection)
        
        ! Apply failsafe flux correction
        call afcstab_failsafeLimiting(rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            SfailsafeVariables, rtimestep%dStep, nfailsafe,&
            euler_getVariable, Rsolution(1), p_rpredictorEuler)

        ! Apply linearised FEM-FCT correction
        call gfsc_buildConvVectorFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))

        !!!
        ! NOTE: We may add failsafe for scalar transport problem here!
        !!!
                
        ! Deallocate temporal memory
        call lsysbl_releaseVector(p_Rpredictor(1))
        call lsysbl_releaseVector(p_Rpredictor(2))
        deallocate(SfailsafeVariables, p_Rpredictor)
        
      else
        
        ! Compute linearised FEM-FCT correction
        call euler_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(1), rcollection)
        
        call gfsc_buildConvVectorFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))
      end if
      

    case (1)   ! minimum synchronisation
      
      ! Compute linearised FEM-FCT correction
      call euler_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
          rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_CORRECT, Rsolution(1), rcollection)
      
      call gfsc_buildConvVectorFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
          rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
          rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_CORRECT, Rsolution(2))
      
      ! Compute minimum correction factor
      call lsyssc_getbase_double(&
          rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorAlpha, p_DalphaEuler)
      call lsyssc_getbase_double(&
          rproblemLevel%Rafcstab(convectionAFC)%p_rvectorAlpha, p_DalphaTransport)
      
      p_DalphaEuler = min(p_DalphaEuler, p_DalphaTransport)
      call lalg_copyVector(p_DalphaEuler, p_DalphaTransport)
      
      if (nfailsafe .gt. 0) then

        ! Apply failsafe flux correction
        call afcstab_failsafeLimiting(rproblemLevel%Rafcstab(IposAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            SfailsafeVariables, rtimestep%dStep, nfailsafe,&
            zpinch_getVariable, Rsolution, p_Rpredictor)
        
        ! Deallocate temporal memory
        call lsysbl_releaseVector(p_Rpredictor(1))
        call lsysbl_releaseVector(p_Rpredictor(2))
        deallocate(SfailsafeVariables, p_Rpredictor)
        
      else
        
        ! Apply linearised FEM-FCT correction
        call euler_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(1), rcollection)
        
        call gfsc_buildConvVectorFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))
      end if


    case (2)   ! Euler first, transport second
      
      ! Compute linearised FEM-FCT correction
      call euler_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
          rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_CORRECT, Rsolution(1), rcollection)
      
      ! Copy correction factor for Euler system to transport model
      call afcstab_duplicateStabilisation(&
          rproblemLevel%Rafcstab(inviscidAFC),&
          rproblemLevel%Rafcstab(convectionAFC), AFCSTAB_DUP_EDGELIMITER)
      
      ! Compute linearised FEM-FCT correction (without initialisation)
      call gfsc_buildConvVectorFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
          rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
          rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_INITALPHA-AFCSTAB_FCTALGO_CORRECT, Rsolution(2))
      
      ! Copy final correction factor back to Euler system
      call afcstab_duplicateStabilisation(&
          rproblemLevel%Rafcstab(convectionAFC),&
          rproblemLevel%Rafcstab(inviscidAFC), AFCSTAB_DUP_EDGELIMITER)
      
      if (nfailsafe .gt. 0) then

        ! Apply failsafe flux correction
        call afcstab_failsafeLimiting(rproblemLevel%Rafcstab(IposAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            SfailsafeVariables, rtimestep%dStep, nfailsafe,&
            zpinch_getVariable, Rsolution, p_Rpredictor)
        
        ! Deallocate temporal memory
        call lsysbl_releaseVector(p_Rpredictor(1))
        call lsysbl_releaseVector(p_Rpredictor(2))
        deallocate(SfailsafeVariables, p_Rpredictor)

      else
        
        ! Apply linearised FEM-FCT correction
        call euler_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(1), rcollection)
        
        call gfsc_buildConvVectorFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))
      end if


    case (3)   ! Transport first, Euler second
      
      ! Compute linearised FEM-FCT correction      
      call gfsc_buildConvVectorFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
          rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
          rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_CORRECT, Rsolution(2))
      
      ! Copy correction factor for transport model to Euler system
      call afcstab_duplicateStabilisation(&
          rproblemLevel%Rafcstab(convectionAFC),&
          rproblemLevel%Rafcstab(inviscidAFC), AFCSTAB_DUP_EDGELIMITER)
      
      ! Compute linearised FEM-FCT correction (without initialisation)
      call euler_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
          rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_INITALPHA-AFCSTAB_FCTALGO_CORRECT,&
          Rsolution(1), rcollection)
      
      ! Copy final correction factor back to transport model
      call afcstab_duplicateStabilisation(&
          rproblemLevel%Rafcstab(inviscidAFC),&
          rproblemLevel%Rafcstab(convectionAFC), AFCSTAB_DUP_EDGELIMITER)

      if (nfailsafe .gt. 0) then

        ! Apply failsafe flux correction
        call afcstab_failsafeLimiting(rproblemLevel%Rafcstab(IposAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            SfailsafeVariables, rtimestep%dStep, nfailsafe,&
            zpinch_getVariable, Rsolution, p_Rpredictor)
        
        ! Deallocate temporal memory
        call lsysbl_releaseVector(p_Rpredictor(1))
        call lsysbl_releaseVector(p_Rpredictor(2))
        deallocate(SfailsafeVariables, p_Rpredictor)

      else
      
        ! Apply linearised FEM-FCT correction
        call euler_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(1), rcollection)
        
        call gfsc_buildConvVectorFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))
      end if
      

    case (4)   ! Euler system only
      
      if (nfailsafe .gt. 0) then

        ! Compute linearised FEM-FCT correction
        call euler_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
            AFCSTAB_FCTALGO_CORRECT, Rsolution(1), rcollection)

        ! Copy correction factor for Euler system to transport model
        call afcstab_duplicateStabilisation(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rafcstab(convectionAFC), AFCSTAB_DUP_EDGELIMITER)

        ! Apply failsafe flux correction
        call afcstab_failsafeLimiting(rproblemLevel%Rafcstab(IposAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            SfailsafeVariables, rtimestep%dStep, nfailsafe,&
            zpinch_getVariable, Rsolution, p_Rpredictor)
        
        ! Deallocate temporal memory
        call lsysbl_releaseVector(p_Rpredictor(1))
        call lsysbl_releaseVector(p_Rpredictor(2))
        deallocate(SfailsafeVariables, p_Rpredictor)
        
      else
        
        ! Apply linearised FEM-FCT correction      
        call euler_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(1), rcollection)
        
        ! Copy correction factor for Euler system to transport model
        call afcstab_duplicateStabilisation(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rafcstab(convectionAFC), AFCSTAB_DUP_EDGELIMITER)
        
        ! Apply linearised FEM-FCT correction
        call gfsc_buildConvVectorFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))
      end if
      
      
    case (5)   ! Transport model only
      
      if (nfailsafe .gt. 0) then

        ! Compute linearised FEM-FCT correction
        call gfsc_buildConvVectorFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
            AFCSTAB_FCTALGO_CORRECT, Rsolution(2))

        ! Copy correction factor for transport model to Euler system
        call afcstab_duplicateStabilisation(&
            rproblemLevel%Rafcstab(convectionAFC),&
            rproblemLevel%Rafcstab(inviscidAFC), AFCSTAB_DUP_EDGELIMITER)

        ! Apply failsafe flux correction
        call afcstab_failsafeLimiting(rproblemLevel%Rafcstab(IposAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrix),&
            SfailsafeVariables, rtimestep%dStep, nfailsafe,&
            zpinch_getVariable, Rsolution, p_Rpredictor)
        
        ! Deallocate temporal memory
        call lsysbl_releaseVector(p_Rpredictor(1))
        call lsysbl_releaseVector(p_Rpredictor(2))
        deallocate(SfailsafeVariables, p_Rpredictor)
        
      else
        
        ! Apply linearised FEM-FCT correction
        call gfsc_buildConvVectorFCT(rproblemLevel%Rmatrix(lumpedMassMatrix),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))
        
        ! Copy correction factor for transport model to Euler system
        call afcstab_duplicateStabilisation(&
            rproblemLevel%Rafcstab(convectionAFC),&
            rproblemLevel%Rafcstab(inviscidAFC), AFCSTAB_DUP_EDGELIMITER)
        
        ! Apply linearised FEM-FCT correction
        call euler_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(1), rcollection)
      end if


    case default
      call output_line('Invalid type of limiter synchronisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcLinearisedFCT')
      call sys_halt()
    end select
    
    
    ! Impose boundary conditions for the solution vector
    select case(rproblemLevel%rtriangulation%ndim)
    case (NDIM1D)
      call bdrf_filterVectorExplicit(rbdrCond(1), Rsolution(1),&
          rtimestep%dTime, euler_calcBoundaryvalues1d)

    case (NDIM2D)
      call bdrf_filterVectorExplicit(rbdrCond(1), Rsolution(1),&
          rtimestep%dTime, euler_calcBoundaryvalues2d)

    case (NDIM3D)
      call bdrf_filterVectorExplicit(rbdrCond(1), Rsolution(1),&
          rtimestep%dTime, euler_calcBoundaryvalues3d)
    end select

    ! Impose boundary conditions for the solution vector
    call bdrf_filterVectorExplicit(rbdrCond(2), Rsolution(2),&
        rtimestep%dTime)
    
  end subroutine zpinch_calcLinearisedFCT

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_getVariable(RvectorBlock, cvariable, rvectorScalar)

!<description>
    ! This subroutine extracts a single variable from the vector of
    ! conservative variables which is stored in interleave of block 
    ! format. In contrast to the corresponding subroutine from the
    ! Euler model, this subroutine accepts arrays of block vectors,
    ! e.g., to pass the solution from the Euler system and the
    ! scalar tracer equation simultaneously.
    !
    ! This subroutine is a wrapper which calls either the getVariable
    ! subroutine for the Euler system or returns the scalar tracer.
!</description>

!<input>
    ! Vector of conservative variables
    type(t_vectorBlock), dimension(:), intent(in) :: RvectorBlock

    ! Identifier for the variable
    character(LEN=*), intent(in) :: cvariable
!</input>

!<inputoutput>
    ! Extracted single variable
    type(t_vectorScalar), intent(inout) :: rvectorScalar
!</inputoutput>
!</subroutine>

    if (trim(cvariable) .eq. 'advect') then
      call lsyssc_copyVector(RvectorBlock(2)%RvectorBlock(1), rvectorScalar)
    else
      call euler_getVariable(RvectorBlock(1), cvariable, rvectorScalar)
    end if

  end subroutine zpinch_getVariable

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcLinearisedFCT_OLD(rbdrCondEuler,&
      rbdrcondTransport, rproblemLevel, rtimestep,&
      rsolutionEuler, rsolutionTransport, rcollection)

!<description>
    ! This subroutine calculates the linearised FCT correction
!</description>

!<input>
    ! boundary condition structure
    type(t_boundaryCondition), intent(in) :: rbdrCondEuler
    type(t_boundaryCondition), intent(in) :: rbdrCondTransport

    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solution vectors
    type(t_vectorBlock), intent(inout) :: rsolutionEuler
    type(t_vectorBlock), intent(inout) :: rsolutionTransport

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixScalar), pointer :: p_rmatrixMass, p_rmatrixTransport
    type(t_parlist), pointer :: p_rparlist
    type(t_afcstab), pointer :: p_rafcstabEuler, p_rafcstabTransport
    type(t_vectorBlock) :: rdataEuler, rdataTransport
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    real(DP), dimension(:,:), pointer :: p_DcoefficientsAtEdge
    real(DP), dimension(:), pointer :: p_ML, p_Cx, p_Cy
    real(DP), dimension(:), pointer :: p_ppEuler, p_pmEuler,&
        p_qpEuler, p_qmEuler, p_rpEuler, p_rmEuler
    real(DP), dimension(:), pointer :: p_ppTransport, p_pmTransport,&
        p_qpTransport, p_qmTransport, p_rpTransport, p_rmTransport
    real(DP), dimension(:), pointer :: p_uEuler, p_uTransport
    real(DP), dimension(:), pointer :: p_fluxEuler, p_dataEuler
    real(DP), dimension(:), pointer :: p_fluxTransport, p_dataTransport, p_alpha
    character(LEN=SYS_STRLEN) :: ssectionNameEuler, ssectionNameTransport
    integer :: transportmatrix, convectionAFC, inviscidAFC
    integer :: lumpedMassMatrix,  consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY


    print *, "HERE"
    stop

    ! Get parameter list
    p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

    ! Get section names
    ssectionNameEuler = rcollection%SquickAccess(1)
    ssectionNameTransport = rcollection%SquickAccess(2)

    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist, ssectionNameEuler,&
        'coeffMatrix_CX', coeffMatrix_CX)
    call parlst_getvalue_int(p_rparlist, ssectionNameEuler,&
        'coeffMatrix_CY', coeffMatrix_CY)
    call parlst_getvalue_int(p_rparlist, ssectionNameEuler,&
        'consistentmassmatrix', consistentMassMatrix)
    call parlst_getvalue_int(p_rparlist, ssectionNameEuler,&
        'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(p_rparlist, ssectionNameEuler,&
        'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(p_rparlist, ssectionNameTransport,&
        'transportmatrix', transportMatrix)
    call parlst_getvalue_int(p_rparlist, ssectionNameTransport,&
        'convectionAFC', convectionAFC)

    ! Set pointers to template matrix and stabilisation structure
    p_rmatrixMass => rproblemLevel%Rmatrix(lumpedMassMatrix)
    p_rmatrixTransport => rproblemLevel%Rmatrix(transportMatrix)
    p_rafcstabTransport => rproblemLevel%Rafcstab(convectionAFC)
    p_rafcstabEuler => rproblemLevel%Rafcstab(inviscidAFC)

    ! Let us check if the linearised FCT algorithm needs to be applied
    if ((p_rafcstabTransport%ctypeAFCstabilisation .ne.&
        AFCSTAB_FEMFCT_LINEARISED) .or.&
        (p_rafcstabEuler%ctypeAFCstabilisation .ne.&
        AFCSTAB_FEMFCT_LINEARISED)) return

    ! Let us check if the edge-based data structure has been generated
    if((iand(p_rafcstabEuler%iSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .and.&
       (iand(p_rafcstabEuler%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(p_rmatrixTransport,&
          p_rafcstabEuler)
    end if

    if((iand(p_rafcstabTransport%iSpec, AFCSTAB_HAS_EDGESTRUCTURE) .eq. 0) .and.&
       (iand(p_rafcstabTransport%iSpec, AFCSTAB_HAS_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(p_rmatrixTransport,&
          p_rafcstabTransport)
    end if

    ! Set pointers
    call lsysbl_getbase_double(rsolutionEuler, p_uEuler)
    call lsysbl_getbase_double(rsolutionTransport, p_uTransport)
    call afcstab_getbase_IverticesAtEdge(p_rafcstabTransport, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(p_rafcstabTransport, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(lumpedMassMatrix), p_ML)
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(coeffMatrix_CX), p_Cx)
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(coeffMatrix_CY), p_Cy)
    call lsyssc_getbase_double(p_rafcstabEuler%RnodalVectors(1), p_ppEuler)
    call lsyssc_getbase_double(p_rafcstabEuler%RnodalVectors(2), p_pmEuler)
    call lsyssc_getbase_double(p_rafcstabEuler%RnodalVectors(3), p_qpEuler)
    call lsyssc_getbase_double(p_rafcstabEuler%RnodalVectors(4), p_qmEuler)
    call lsyssc_getbase_double(p_rafcstabEuler%RnodalVectors(5), p_rpEuler)
    call lsyssc_getbase_double(p_rafcstabEuler%RnodalVectors(6), p_rmEuler)
    call lsyssc_getbase_double(p_rafcstabEuler%RedgeVectors(1), p_fluxEuler)
    call lsyssc_getbase_double(p_rafcstabEuler%RedgeVectors(2), p_alpha)
    call lsyssc_getbase_double(p_rafcstabTransport%RnodalVectors(1), p_ppTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RnodalVectors(2), p_pmTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RnodalVectors(3), p_qpTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RnodalVectors(4), p_qmTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RnodalVectors(5), p_rpTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RnodalVectors(6), p_rmTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RedgeVectors(1), p_fluxTransport)

    ! Create auxiliary vectors
    call lsysbl_createVectorBlock(rsolutionEuler, rdataEuler, .true.)
    call lsysbl_createVectorBlock(rsolutionTransport, rdataTransport, .true.)
    call lsysbl_getbase_double(rdataEuler, p_dataEuler)
    call lsysbl_getbase_double(rdataTransport, p_dataTransport)


    ! Initialize alpha with ones
    call lalg_setVector(p_alpha, 1.0_DP)

    ! Build the raw antidiffusive fluxes
    ! - Euler model -
    call buildFluxEuler2d(p_IverticesAtEdge, p_rmatrixMass%NEQ,&
        NVAR2D, p_rafcstabEuler%NEDGE, p_Cx, p_Cy, p_uEuler, p_fluxEuler)

    ! - transport model -
    call buildFluxTransport2d(p_IverticesAtEdge, p_DcoefficientsAtEdge,&
        p_rmatrixMass%NEQ, p_rafcstabTransport%NEDGE, p_uTransport,&
        p_fluxTransport)


    ! Build the correction coefficients
    ! - density -
    call buildCorrection(p_IverticesAtEdge, p_rmatrixMass%NEQ, NVAR2D,&
        p_rafcstabEuler%NEDGE, 1, rtimestep%dStep, p_ML, p_uEuler,&
        p_fluxEuler, p_alpha, p_ppEuler, p_pmEuler,&
        p_qpEuler, p_qmEuler, p_rpEuler, p_rmEuler)

    ! - pressure -
    call buildCorrection(p_IverticesAtEdge, p_rmatrixMass%NEQ, NVAR2D,&
        p_rafcstabEuler%NEDGE, 4, rtimestep%dStep, p_ML, p_uEuler,&
        p_fluxEuler, p_alpha, p_ppEuler, p_pmEuler,&
        p_qpEuler, p_qmEuler , p_rpEuler, p_rmEuler)

    ! - tracer -
    call buildCorrection(p_IverticesAtEdge, p_rmatrixMass%NEQ, 1,&
        p_rafcstabTransport%NEDGE, 1, rtimestep%dStep, p_ML,&
        p_uTransport, p_fluxTransport, p_alpha, p_ppTransport,&
        p_pmTransport, p_qpTransport, p_qmTransport, p_rpTransport,&
        p_rmTransport)


    ! Apply correction to low-order solution
    ! - Euler model -
    call applyCorrection(p_IverticesAtEdge, p_rmatrixMass%NEQ, NVAR2D,&
        p_rafcstabEuler%NEDGE, rtimestep%dStep, p_ML,&
        p_fluxEuler, p_alpha, p_uEuler, p_dataEuler)

    ! - tranport model -
    call applyCorrection(p_IverticesAtEdge, p_rmatrixMass%NEQ, 1,&
        p_rafcstabTransport%NEDGE, rtimestep%dStep, p_ML,&
        p_fluxTransport, p_alpha, p_uTransport, p_dataTransport)


    ! Set boundary conditions explicitly
    call bdrf_filterVectorExplicit(rbdrCondEuler, rsolutionEuler,&
        rtimestep%dTime, euler_calcBoundaryvalues2d)
    call bdrf_filterVectorExplicit(rbdrCondTransport,&
        rsolutionTransport, rtimestep%dTime)

    ! Release flux vectors
    call lsysbl_releaseVector(rdataEuler)
    call lsysbl_releaseVector(rdataTransport)

  contains

    !***************************************************************************

!<subroutine>

    subroutine buildFluxEuler2d(IverticesAtEdge, NEQ, NVAR, NEDGE, Cx, Cy, u, flux)

!<input>
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEQ,NVAR,NEDGE
      real(DP), dimension(:), intent(in) :: Cx,Cy
      real(DP), dimension(NVAR,NEQ), intent(in) :: u
!</input>

      !<output>
      real(DP), dimension(NVAR,NEDGE), intent(out) :: flux
!</output>
!</subroutine>

      ! local variables
      real(DP), dimension(NVAR) :: K_ij,K_ji,D_ij,Diff
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      integer :: ij,ji,i,j,iedge


      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        ji = IverticesAtEdge(4, iedge)

        ! Compute coefficients
        C_ij(1) = Cx(ij); C_ji(1) = Cx(ji)
        C_ij(2) = Cy(ij); C_ji(2) = Cy(ji)

        ! Calculate diffusion coefficient
!!$        call euler_calcMatRusDissMatD2d(u(:,i), u(:,j),&
!!$            C_ij, C_ji, i, j, 1.0_DP, K_ij, K_ji, D_ij)

        ! Compute solution difference
        Diff = u(:,i)-u(:,j)

        ! Compute the raw antidiffusive flux
        flux(:,iedge) = D_ij*Diff
      end do

    end subroutine buildFluxEuler2d

    !***************************************************************************

!<subroutine>

    subroutine buildFluxTransport2d(IverticesAtEdge,&
        DcoefficientsAtEdge, NEQ, NEDGE, u, flux)

!<input>
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      integer, intent(in) :: NEQ,NEDGE
      real(DP), dimension(:), intent(in) :: u
!</input>

!<output>
      real(DP), dimension(:), intent(out) :: flux
!</output>
!</subroutine>

      ! local variables
      real(DP) :: d_ij
      integer :: i,j,ij,iedge


      !$omp parallel do private(i,j,ij,d_ij)
      do iedge = 1, NEDGE

        ! Get node numbers and matrix position
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)

        ! Get artificial diffusion coefficient
        d_ij = DcoefficientsAtEdge(1, iedge)

        ! Compute raw antidiffusive flux
        flux(iedge) = d_ij * (u(i)-u(j))

      end do
      !$omp end parallel do

    end subroutine buildFluxTransport2d

    !***************************************************************************

!<subroutine>

    subroutine  buildCorrection(IverticesAtEdge, NEQ, NVAR, NEDGE,&
        ivar, dscale, ML, u, flux, alpha, pp, pm, qp, qm, rp, rm)

!<input>
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEQ,NVAR,NEDGE,ivar
      real(DP), intent(in) :: dscale
      real(DP), dimension(:), intent(in) :: ML
      real(DP), dimension(NVAR,NEQ), intent(in) :: u
!</input>

!<inputoutput>
      real(DP), dimension(:), intent(inout) :: alpha
      real(DP), dimension(NVAR,NEDGE), intent(inout) :: flux
!</inputoutput>

!<output>
      real(DP), dimension(:), intent(out) :: pp,pm,qp,qm,rp,rm
!</output>
!</subroutine>


      ! local variables
      real(DP) :: f_ij,f0_ij,diff,aux
      integer :: i,j,iedge

      real(DP) :: u_i,v_i,u_j,v_j,p_ij,p_ji,r_i,r_j,p_i ,p_j

      ! Initialize vectors
      call lalg_clearVector(pp)
      call lalg_clearVector(pm)
      call lalg_clearVector(qp)
      call lalg_clearVector(qm)
      call lalg_setVector(rp, 1.0_DP)
      call lalg_setVector(rm, 1.0_DP)

      select case(ivar)

      case (4)
        ! - pressure -

        ! Loop over all edges
        do iedge = 1, NEDGE

          ! Get node numbers and matrix position
          i  = IverticesAtEdge(1, iedge)
          j  = IverticesAtEdge(2, iedge)

          ! Velocities
          u_i = u(2,i)/u(1,i);   v_i = u(3,i)/u(1,i)
          u_j = u(2,j)/u(1,j);   v_j = u(3,j)/u(1,j)

          ! Solution differences
          p_i = (GAMMA-1) * ( u(4,i) - 0.5 * (u_i*u_i + v_i*v_i) * u(1,i) )
          p_j = (GAMMA-1) * ( u(4,j) - 0.5 * (u_j*u_j + v_j*v_j) * u(1,j) )

          ! Pressure difference
          diff = p_j-p_i

          ! Pressure variables
          p_ij = (GAMMA-1) * ( flux(4,iedge) + 0.5 * (u_i*u_i + v_i*v_i)*flux(1,iedge) -&
                               u_i*flux(2,iedge) - v_i*flux(3,iedge) )
          p_ji =-(GAMMA-1) * ( flux(4,iedge) + 0.5 * (u_j*u_j + v_j*v_j)*flux(1,iedge) -&
                               u_j*flux(2,iedge) - v_j*flux(3,iedge) )

          ! Sums of raw antidiffusive fluxes
          pp(i) = pp(i) + max(0.0_DP, p_ij)
          pp(j) = pp(j) + max(0.0_DP, p_ji)
          pm(i) = pm(i) + min(0.0_DP, p_ij)
          pm(j) = pm(j) + min(0.0_DP, p_ji)

          ! Sums of admissible edge contributions
          qp(i) = max(qp(i),  diff)
          qp(j) = max(qp(j), -diff)
          qm(i) = min(qm(i),  diff)
          qm(j) = min(qm(j), -diff)

        end do

      case default
        ! - conservative variable -

        ! Loop over all edges
        do iedge = 1, NEDGE

          ! Get node numbers and matrix position
          i = IverticesAtEdge(1, iedge)
          j = IverticesAtEdge(2, iedge)

          ! Flux correction in conservative variables
          f_ij = flux(ivar,iedge)
          diff = u(ivar,j)-u(ivar,i)

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

      end select


      ! Compute nodal correction factors
      !$omp parallel do
      do i = 1, NEQ
        qp(i) = qp(i)*ML(i)/dscale
        qm(i) = qm(i)*ML(i)/dscale

        rp(i) = min(1._DP, qp(i)/(pp(i)+SYS_EPSREAL) )
        rm(i) = min(1._DP, qm(i)/(pm(i)-SYS_EPSREAL) )
      end do
      !$omp end parallel do


      select case(ivar)

      case(4)
        ! - pressure -

        ! Loop over all edges
        !$omp parallel do private(i,j,f_ij,p_ij,p_ji,u_i,v_i,u_j,v_j,r_i,r_j)
        do iedge = 1, NEDGE

          ! Get node numbers
          i  = IverticesAtEdge(1, iedge)
          j  = IverticesAtEdge(2, iedge)

          ! Velocities
          u_i = u(2,i)/u(1,i);   v_i = u(3,i)/u(1,i)
          u_j = u(2,j)/u(1,j);   v_j = u(3,j)/u(1,j)

          ! Pressure variables
          p_ij = (GAMMA-1) * ( flux(4,iedge) + 0.5 * (u_i*u_i + v_i*v_i)*flux(1,iedge) -&
                               u_i*flux(2,iedge) - v_i*flux(3,iedge) )
          p_ji =-(GAMMA-1) * ( flux(4,iedge) + 0.5 * (u_j*u_j + v_j*v_j)*flux(1,iedge) -&
                               u_j*flux(2,iedge) - v_j*flux(3,iedge) )


          if (p_ij > SYS_EPSREAL) then
            r_i = rp(i)
          elseif (p_ij < -SYS_EPSREAL) then
            r_i = rm(i)
          else
            r_i = 1.0
          end if

          if (p_ji > SYS_EPSREAL) then
            r_j = rp(j)
          elseif (p_ji < -SYS_EPSREAL) then
            r_j = rm(j)
          else
            r_j = 1.0
          end if

          ! Compute correction factor
          alpha(iedge) = min(alpha(iedge), r_i, r_j)

        end do
        !$omp end parallel do


      case default
        ! - conservative variable -

        ! Loop over all edges
        !$omp parallel do private(i,j,f_ij)
        do iedge = 1, NEDGE

          ! Get node numbers
          i  = IverticesAtEdge(1, iedge)
          j  = IverticesAtEdge(2, iedge)

          ! Flux correction in conservative variables
          f_ij = flux(ivar,iedge)

          ! Limit conservative fluxes
          if (f_ij > SYS_EPSREAL) then
            alpha(iedge) = min(alpha(iedge), rp(i), rm(j))
          elseif (f_ij < SYS_EPSREAL) then
            alpha(iedge) = min(alpha(iedge), rm(i), rp(j))
          end if
        end do
        !$omp end parallel do

      end select

    end subroutine buildCorrection

    !***************************************************************************

!<subroutine>

    subroutine applyCorrection(IverticesAtEdge, NEQ, NVAR, NEDGE,&
        dscale, ML, flux, alpha, u, data)

      !<input>
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEQ,NVAR,NEDGE
      real(DP), intent(in) :: dscale
      real(DP), dimension(:), intent(in) :: ML
      real(DP), dimension(NVAR,NEDGE), intent(in) :: flux
!<input>

!<inputoutput>
      real(DP), dimension(:), intent(inout) :: alpha
      real(DP), dimension(NVAR,NEQ), intent(inout) :: u
!</inputoutput>

!<output>
      real(DP), dimension(NVAR,NEQ), intent(out) :: data
!</output>
!</subroutine>

      ! local variables
      real(DP), dimension(NVAR) :: f_ij
      integer :: i,j,iedge

      ! Initialize correction
      call lalg_clearVector(data)

      ! Loop over all edges
      do iedge = 1, NEDGE

        ! Get node numbers and matrix position
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)

        ! Limit raw antidiffusive flux
        f_ij = alpha(iedge)*flux(:,iedge)

        ! Apply correction
        data(:,i) = data(:,i) + f_ij
        data(:,j) = data(:,j) - f_ij
      end do

      ! Loop over all rows
      !$omp parallel do
      do i = 1, NEQ

        ! Compute flux-corrected solution
        u(:,i) = u(:,i) + dscale * data(:,i)/ML(i)
      end do
      !$omp end parallel do

    end subroutine applyCorrection

  end subroutine zpinch_calcLinearisedFCT_OLD

  !*****************************************************************************

  function zpinch_checkPressure(rvector) result(bpositive)

    type(t_vectorBlock), intent(in) :: rvector

    logical :: bpositive

    ! local variabels
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: neq

    ! Compute number of equations
    neq = rvector%NEQ/NVAR2D

    call lsysbl_getbase_double(rvector, p_Ddata)

    bpositive =  docheck(4, neq, p_Ddata)

  contains

    function docheck(nvar, neq, data) result(bpositive)

      integer, intent(in) :: nvar, neq
      real(DP), dimension(nvar, neq), intent(in) :: data

      logical :: bpositive

      real(DP) :: p
      integer :: ieq

      bpositive = .true.

      !$omp parallel do private(p)
      do ieq = 1, neq

        p = 0.4_DP * (data(4,ieq) - 0.5_DP * &
            (data(2,ieq)**2 + data(3,ieq)**2) / data(1,ieq) )

        if (p .le. 0.0_DP) bpositive = .false.
      end do
      !$omp end parallel do

    end function docheck

  end function zpinch_checkPressure

end module zpinch_callback

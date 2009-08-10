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
!# 2.) zpinch_initVelocityField
!#     -> Initializes the velocity field for the transport model
!#
!# 3.) zpinch_initDensityAveraging
!#     -> Initializes the density averaged mass matrices
!#
!# 4.) zpinch_initLorentzforceTerm
!#     -> Initializes the Lorentz force term for given solutions
!#
!# 5.) zpinch_calcLinearizedFCT
!#     -> Calculates the linearized FCT correction
!#
!# </purpose>
!##############################################################################

module zpinch_callback

  use afcstabilisation
  use boundaryfilter
  use collection
  use derivatives
  use euler_basic
  use euler_callback
  use euler_callback2d
  use flagship_basic
  use fparser
  use fsystem
  use genoutput
  use basicgeometry
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
  use transport_callback2d
  use trilinearformevaluation
  use zpinch_callback2d

  implicit none

  private
  public :: zpinch_nlsolverCallback
  public :: zpinch_initVelocityField
  public :: zpinch_initDensityAveraging
  public :: zpinch_initLorentzforceTerm
  public :: zpinch_calcLinearizedFCT
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
    real(DP) :: dscale
    integer(i32) :: iSpec
    integer :: jacobianMatrix
    
    select case(rsolver%ssolverName)
      
    case('NonlinearSolverEuler')

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
        
        ! Calculate scaling for implicit part of the Lorentz force
        dscale = -rtimestep%theta * rtimestep%dStep
        
        if (dscale .ne. 0.0_DP) then
          
          ! Set pointer to parameter list
          p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

          ! Compute the implicit part of the Lorentz force
          call zpinch_initLorentzforceTerm(p_rparlist,&
              rcollection%SquickAccess(2), rcollection%SquickAccess(3),&
              rcollection%SquickAccess(4), rproblemLevel,&
              p_rsolutionEuler, p_rsolutionTransport, rtimestep%dTime,&
              dscale, rforce, rcollection)
          
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
              rsolver, rsolution, rsolution0, rrhs, rres, istep,&
              rcollection)
          
        end if
      end if
      
      
      ! Do we have to impose boundary conditions?
      ! --------------------------------------------------------------------------
      if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then
        
        ! Impose boundary conditions
        call euler_setBoundaryConditions(rproblemLevel, rtimestep,&
            rsolver, rsolution, rsolution0, rres, rcollection)
      end if
      

    case('NonlinearSolverTransport')
      
      ! Set the first string quick access array to the section name
      ! of the Euler model which is stored in the fourth array
      rcollection%SquickAccess(1) = rcollection%SquickAccess(4)

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
        call transp_calcPrecondThetaScheme(rproblemLevel,&
            rtimestep, rsolver, rsolution, rcollection,&
            zpinch_calcMatRusConvectionP2d,&
            zpinch_calcMatRusConvectionD2d,&
            transp_coeffMatBdrConvectionP2d,&
            transp_coeffMatBdrConvectionD2d)

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

            ! Compute the preconditioner
            call transp_calcPrecondThetaScheme(rproblemLevel,&
                  rtimestep, rsolver, rsolution0, rcollection,&
                  zpinch_calcMatRusConvectionP2d,&
                  zpinch_calcMatRusConvectionD2d,&
                  transp_coeffMatBdrConvectionP2d,&
                  transp_coeffMatBdrConvectionD2d)
              
            ! Assemble the constant right-hand side
            call transp_calcRhsThetaScheme(rproblemLevel, rtimestep,&
                rsolver, rsolution0, rrhs, rcollection,&
                fcb_coeffVecBdrPrimal_sim=transp_coeffVecBdrConvectionP2d,&
                fcb_coeffVecBdrDual_sim=transp_coeffVecBdrConvectionD2d)
            
          end if

          ! Set pointer to parameter list
          p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

          ! Set pointer to the solution vector of the Euler model
          ! from the current time step
          p_rsolutionEuler => rcollection%p_rvectorQuickAccess1
          
          ! Calculate the density averaged mass matrices using the
          ! solution of the Euler model from the current time step
          call zpinch_initDensityAveraging(p_rparlist,&
              rcollection%SquickAccess(3), rcollection%SquickAccess(4),&
              rproblemLevel, p_rsolutionEuler, rcollection)
          
          ! Calculate the velocity vector using the solution of the
          ! Euler model from the current time step
          call zpinch_initVelocityField(p_rparlist,&
              rcollection%SquickAccess(1), rproblemLevel,&
              p_rsolutionEuler, rcollection)

          ! Compute the preconditioner
          call transp_calcPrecondThetaScheme(rproblemLevel,&
              rtimestep, rsolver, rsolution, rcollection,&
              zpinch_calcMatRusConvectionP2d,&
              zpinch_calcMatRusConvectionD2d,&
              transp_coeffMatBdrConvectionP2d,&
              transp_coeffMatBdrConvectionD2d)

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
        
        ! Compute the preconditioner
        call transp_calcPrecondThetaScheme(rproblemLevel,&
            rtimestep, rsolver, rsolution, rcollection,&
            zpinch_calcMatRusConvectionP2d,&
            zpinch_calcMatRusConvectionD2d,&
            transp_coeffMatBdrConvectionP2d,&
            transp_coeffMatBdrConvectionD2d)
      end if
      
      
      ! Do we have to calculate the Jacobian operator?
      ! --------------------------------------------------------------------------
      if (iand(iSpec, NLSOL_OPSPEC_CALCJACOBIAN) .ne. 0) then
        
        ! Compute the Jacobian matrix
        call transp_calcJacobianThetaScheme(rproblemLevel,&
            rtimestep, rsolver, rsolution, rsolution0, rcollection,&
            zpinch_calcMatRusConvectionP2d,&
            zpinch_calcMatRusConvectionD2d,&
            transp_coeffMatBdrConvectionP2d,&
            transp_coeffMatBdrConvectionD2d)
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
        call parlst_getvalue_int(p_rparlist, rcollection&
            %SquickAccess(1), 'jacobianMatrix', jacobianMatrix)
        
        ! Apply Jacobian matrix
        call lsyssc_scalarMatVec(&
            rproblemLevel%Rmatrix(jacobianMatrix),&
            rsolution%RvectorBlock(1), rres%RvectorBlock(1),&
            1.0_DP, 1.0_DP)
      end if
      
    end select
    
    ! Set status flag
    istatus = 0
    
  end subroutine zpinch_nlsolverCallback

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_initVelocityField(rparlist, ssectionName,&
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
    
!!$    ! Set x-velocity, i.e., momentum in x-direction
!!$    call euler_getVariable(rsolution, 'momentum_x',&
!!$        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(1))
!!$    call zpinch_setVariable2d(&
!!$        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(1), 1)
!!$    
!!$    ! Set y-velocity, i.e., momentum in y-direction
!!$    call euler_getVariable(rsolution, 'momentum_y',&
!!$        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(2))
!!$    call zpinch_setVariable2d(&
!!$        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(2), 2)

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
    call lsysbl_duplicateVector(rsolution, rproblemLevel&
        %RvectorBlock(2), LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_COPY)
    call zpinch_setVariable2d(rproblemLevel%RvectorBlock(2), 3)
    
    ! Set update notification in problem level structure
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     PROBLEV_MSPEC_UPDATE)

  end subroutine zpinch_initVelocityField

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_initDensityAveraging(rparlist,&
      ssectionNameEuler, ssectionNameTransport,&
      rproblemlevel, rsolutionEuler, rcollection)

!<description>
    ! This subroutine initializes the density averaged mass matrices
    ! for the transport model based on the solution from the Euler
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section names in parameter list
    character(LEN=*), intent(in) :: ssectionNameEuler
    character(LEN=*), intent(in) :: ssectionNameTransport

    ! solution vector for Euler model
    type(t_vectorBlock), intent(in) :: rsolutionEuler
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_vectorScalar) :: rvector
    real(DP), dimension(:), pointer :: p_ML, p_Density
    integer :: lumpedMassMatrix, lumpedMassMatrixDensity
    integer :: i
    
    ! OK, IN THE CURRENT IMPLEMENTATION WE DO NOT NEED A
    ! DENSITY_AVERAGED MASS MATRIX, HENCE RETURN
 
!!$    ! Get global configuration from parameter list
!!$    call parlst_getvalue_int(rparlist,&
!!$        ssectionNameEuler, 'lumpedmassmatrix', lumpedMassMatrix)
!!$    call parlst_getvalue_int(rparlist,&
!!$        ssectionNameTransport, 'lumpedmassmatrix', lumpedMassMatrixDensity)
!!$    
!!$    ! Get density distribution from the solution of the Euler model
!!$    ! and create block vector which is attached to the collection
!!$    call euler_getVariable(rsolutionEuler, 'density', rvector)
!!$    call lsyssc_getbase_double(rvector, p_Density)
!!$    
!!$    ! Create density averaged lumped mass matrix
!!$    call lsyssc_duplicateMatrix(&
!!$        rproblemLevel%Rmatrix(lumpedMassMatrix),&
!!$        rproblemLevel%Rmatrix(lumpedMassMatrixDensity),&
!!$        LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
!!$    
!!$    call lsyssc_getbase_double(&
!!$        rproblemLevel%Rmatrix(lumpedMassMatrixDensity), p_ML)
!!$    
!!$    !$omp parallel do
!!$    do i = 1, size(p_Density)
!!$      p_ML(i) = p_ML(i)*p_Density(i)
!!$    end do
!!$    !$omp end parallel do
!!$    
!!$    ! Release temporal vector
!!$    call lsyssc_releaseVector(rvector)
    
  end subroutine zpinch_initDensityAveraging

    !*****************************************************************************

!<subroutine>

  subroutine zpinch_initLorentzforceTerm(rparlist, ssectionName,&
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
    real(DP) :: dcurrentDrive
    integer :: isystemFormat, lumpedMassMatrix
    integer :: neq, nvar, icomp, icoords
    logical :: bcompatible


    ! Check if solution vector and Lorentz force vector are compatible
    call lsysbl_isVectorCompatible(rsolutionEuler, rforce, bcompatible)
    if (.not.bcompatible)&
        call lsysbl_resizeVectorBlock(rforce, rsolutionEuler, .false.)
    
    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName,&
        'slorentzforcename', slorentzforceName)
    call parlst_getvalue_int(rparlist, ssectionName,&
        'icoords', icoords)
    call parlst_getvalue_int(rparlist, ssectionNameEuler,&
        'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(rparlist, ssectionNameEuler,&
        'isystemformat', isystemFormat)
    
    ! Get lumped and consistent mass matrix
    call lsyssc_getbase_double(&
        rproblemLevel%Rmatrix(lumpedMassMatrix), p_DlumpedMassMatrix)
    
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
    call fparser_evalFunction(p_rfparser, icomp,&
        (/dtime/), dcurrentDrive)
    
    ! Multiply scaling parameter by the time step
    dcurrentDrive = dscale * dcurrentDrive

    ! What type of system format are we?
    select case(isystemFormat)

    case (SYSTEM_INTERLEAVEFORMAT)

      ! What type of coordinate system are we?
      select case(icoords)
      case(1)
        call calcForceXYInterleaveFormat(dcurrentDrive, neq, nvar,&
            p_DvertexCoords, p_DlumpedMassMatrix, p_DdataTransport,&
            p_DdataEuler, p_Ddata)

      case(2)
        call calcForceRZInterleaveFormat(dcurrentDrive, neq, nvar,&
            p_DvertexCoords, p_DlumpedMassMatrix, p_DdataTransport,&
            p_DdataEuler, p_Ddata)
        
      case default
        call output_line('Invalid type of coordinate system!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_initLorentzforceTerm')
        call sys_halt()
      end select

    case (SYSTEM_BLOCKFORMAT)
      ! What type of coordinate system are we?
      select case(icoords)
      case(1)
        call calcForceXYInterleaveFormat(dcurrentDrive, neq, nvar,&
            p_DvertexCoords, p_DlumpedMassMatrix, p_DdataTransport,&
            p_DdataEuler, p_Ddata)
        
      case(2)
        call calcForceXYInterleaveFormat(dcurrentDrive, neq, nvar,&
            p_DvertexCoords, p_DlumpedMassMatrix, p_DdataTransport,&
            p_DdataEuler, p_Ddata)
        
      case default
        call output_line('Invalid type of coordinate system!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_initLorentzforceTerm')
        call sys_halt()
      end select

    case DEFAULT
      call output_line('Invalid system format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'zpinch_initLorentzforceTerm')
      call sys_halt()
    end select


  contains
    
    ! Here, the real working routines follow
    
    !**************************************************************
    ! Calculate the Lorentz force term in x-y coordinates.
    ! The system is stored in interleave format.
    
    subroutine calcForceXYInterleaveFormat(dconst, neq, nvar,&
        DvertexCoords, DmassMatrix, DdataTransport, DdataEuler,&
        DdataForce)
      
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(nvar,neq), intent(in) :: DdataEuler
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      real(DP), intent(in) :: dconst
      integer, intent(in) :: neq, nvar
      
      real(DP), dimension(nvar,neq), intent(out) :: DdataForce
      
      ! local variables
      real(DP) :: drad, dang, daux1, daux2, x1, x2
      integer :: i
      
      ! Loop over all rows
      !$omp parallel do private(x1,x2,drad,daux1,daux2)
      do i = 1, neq
        
        ! Get coordinates at node i
        x1 = DvertexCoords(1, i)
        x2 = DvertexCoords(2, i)

        ! Compute unit vector
        drad = sqrt(x1*x1 + x2*x2)
        
!!$        ! Compute polar coordinates
!!$        drad = sqrt(x1*x1 + x2*x2)
!!$        dang = atan2(x2, x1)
        
        ! Compute unit vector into origin
        if (drad .gt. 1e-4) then
          x1 = x1/drad
          x2 = x2/drad
        else
          x1 = 0.0; x2 = 0.0
        end if
        
!!$        ! Compute source term
!!$        daux1 = dconst * DmassMatrix(i) * DdataTransport(i) *&
!!$                         DdataEuler(1, i) / max(drad, 1.0e-4_DP)
!!$        daux2 = dconst * DmassMatrix(i) * DdataTransport(i) *&
!!$                         (DdataEuler(2, i) * x1 + DdataEuler(3,i) * x2) /&
!!$                         max(drad, 1.0e-4_DP)

        ! Compute source term
        daux1 = dconst * DmassMatrix(i) * DdataTransport(i) / max(drad, 1.0e-4_DP)
        daux2 = dconst * DmassMatrix(i) * DdataTransport(i) / DdataEuler(1, i) *&
                         (DdataEuler(2, i) * x1 + DdataEuler(3,i) * x2) /&
                         max(drad, 1.0e-4_DP)
               
        ! Impose source values
        DdataForce(1, i) = 0.0_DP
        DdataForce(2, i) = daux1 * x1
        DdataForce(3, i) = daux1 * x2
        DdataForce(4, i) = daux2
      end do
      !$omp end parallel do
      
    end subroutine calcForceXYInterleaveFormat

    !**************************************************************
    ! Calculate the Lorentz force term in r-z coordinates.
    ! The system is stored in interleave format.
    
    subroutine calcForceRZInterleaveFormat(dconst, neq, nvar,&
        DvertexCoords, DmassMatrix, DdataTransport, DdataEuler,&
        DdataForce)
      
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(nvar,neq), intent(in) :: DdataEuler
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      real(DP), intent(in) :: dconst
      integer, intent(in) :: neq, nvar
      
      real(DP), dimension(nvar,neq), intent(out) :: DdataForce
      
      ! local variables
      real(DP) :: drad, dang, daux1, daux2, x1, x2
      integer :: i
      
      ! Loop over all rows
      !$omp parallel do private(x1,x2,drad,daux1,daux2)
      do i = 1, neq
        
        ! Get coordinates at node i
        x1 = DvertexCoords(1, i)
        x2 = DvertexCoords(2, i)
        
        ! Get coordinates at node i
        drad = x1
        
        ! Compute unit vector into origin
        if (drad .gt. 1e-4) then
          x1 = 1.0; x2 = 0.0
        else
          x1 = 0.0; x2 = 0.0
        end if
        
!!$        ! Compute source term
!!$        daux1 = dconst * DmassMatrix(i) * DdataTransport(i) *&
!!$                         DdataEuler(1, i) / max(drad, 1.0e-4_DP)
!!$        daux2 = dconst * DmassMatrix(i) * DdataTransport(i) *&
!!$                         (DdataEuler(2, i) * x1 + DdataEuler(3,i) * x2) /&
!!$                         max(drad, 1.0e-4_DP)

        ! Compute source term
        daux1 = dconst * DmassMatrix(i) * DdataTransport(i) / max(drad, 1.0e-4_DP)
        daux2 = dconst * DmassMatrix(i) * DdataTransport(i) / DdataEuler(1, i) *&
                         (DdataEuler(2, i) * x1 + DdataEuler(3,i) * x2) /&
                         max(drad, 1.0e-4_DP)
        
        ! Impose source values
        DdataForce(1, i) = 0.0_DP
        DdataForce(2, i) = daux1 * x1
        DdataForce(3, i) = daux1 * x2
        DdataForce(4, i) = daux2
      end do
      !$omp end parallel do

    end subroutine calcForceRZInterleaveFormat

  end subroutine zpinch_initLorentzforceTerm

  !*****************************************************************************
  
  !<subroutine>

  subroutine zpinch_calcLinearizedFCT(rbdrCondEuler,&
      rbdrcondTransport, rproblemLevel, rtimestep,&
      rsolutionEuler, rsolutionTransport, rcollection)

!<description>
    ! This subroutine calculates the linearized FCT correction
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
    real(DP), dimension(:), pointer :: p_MC, p_ML, p_Cx, p_Cy
    real(DP), dimension(:), pointer :: p_ppEuler, p_pmEuler,&
        p_qpEuler, p_qmEuler, p_rpEuler, p_rmEuler
    real(DP), dimension(:), pointer :: p_ppTransport, p_pmTransport,&
        p_qpTransport, p_qmTransport, p_rpTransport, p_rmTransport
    real(DP), dimension(:), pointer :: p_uEuler, p_uTransport
    real(DP), dimension(:), pointer :: p_flux0Euler, p_fluxEuler, p_dataEuler
    real(DP), dimension(:), pointer :: p_flux0Transport, p_fluxTransport,&
        p_dataTransport, p_alpha
    integer, dimension(:), pointer :: p_Kdiagonal
    character(LEN=SYS_STRLEN) :: ssectionNameEuler, ssectionNameTransport
    integer :: transportmatrix, convectionAFC, inviscidAFC
    integer :: lumpedMassMatrix,  consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY

    real(DP) :: d1

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
    
    ! Let us check if the edge-based data structure has been generated
    if((iand(p_rafcstabEuler%iSpec, AFCSTAB_EDGESTRUCTURE) .eq. 0) .and.&
       (iand(p_rafcstabEuler%iSpec, AFCSTAB_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(p_rmatrixTransport,&
          p_rafcstabEuler)
    end if

    if((iand(p_rafcstabTransport%iSpec, AFCSTAB_EDGESTRUCTURE) .eq. 0) .and.&
       (iand(p_rafcstabTransport%iSpec, AFCSTAB_EDGEORIENTATION) .eq. 0)) then
      call afcstab_generateVerticesAtEdge(p_rmatrixTransport,&
          p_rafcstabTransport)
    end if

    ! Set pointers
    call lsysbl_getbase_double(rsolutionEuler, p_uEuler)
    call lsysbl_getbase_double(rsolutionTransport, p_uTransport)
    call lsyssc_getbase_Kdiagonal(p_rmatrixTransport, p_Kdiagonal)
    call afcstab_getbase_IverticesAtEdge(p_rafcstabTransport, p_IverticesAtEdge)
    call afcstab_getbase_DcoeffsAtEdge(p_rafcstabTransport, p_DcoefficientsAtEdge)
    call lsyssc_getbase_double(rproblemLevel%Rmatrix(consistentMassMatrix), p_MC)
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
    call lsyssc_getbase_double(p_rafcstabEuler%RedgeVectors(2), p_flux0Euler)
    call lsyssc_getbase_double(p_rafcstabEuler%RedgeVectors(3), p_alpha)
    call lsyssc_getbase_double(p_rafcstabTransport%RnodalVectors(1), p_ppTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RnodalVectors(2), p_pmTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RnodalVectors(3), p_qpTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RnodalVectors(4), p_qmTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RnodalVectors(5), p_rpTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RnodalVectors(6), p_rmTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RedgeVectors(1), p_fluxTransport)
    call lsyssc_getbase_double(p_rafcstabTransport%RedgeVectors(2), p_flux0Transport)

    ! Create auxiliary vectors
    call lsysbl_createVectorBlock(rsolutionEuler, rdataEuler, .true.)
    call lsysbl_createVectorBlock(rsolutionTransport, rdataTransport, .true.)
    call lsysbl_getbase_double(rdataEuler, p_dataEuler)
    call lsysbl_getbase_double(rdataTransport, p_dataTransport)

    ! Compute low-order approximation for transport model
    call lsyssc_scalarMatVec(p_rmatrixTransport,&
        rsolutionTransport%RvectorBlock(1),&
        rdataTransport%RvectorBlock(1), 1.0_DP, 0.0_DP)
    call lsyssc_invertedDiagMatVec(p_rmatrixMass,&
        rdataTransport%RvectorBlock(1), 1.0_DP,&
        rdataTransport%RvectorBlock(1))
    
    ! Compute force term
    call zpinch_initLorentzforceTerm(p_rparlist, 'zpinch',&
        ssectionNameEuler, ssectionNameTransport, rproblemLevel,&
        rsolutionEuler, rsolutionTransport, rtimestep%dTime, 1.0_DP,&
        rdataEuler, rcollection)
   

    ! CONSISTENCY CHECK #1:
    d1 = minval(p_uTransport)
    if (d1 .le. -1e-8) then
      print *, "Tracer has become negative in the low-order method", d1
      pause
    end if

    ! Initialize alpha with ones
    call lalg_setVector(p_alpha, 1.0_DP)

    
    ! Build the flux for the Euler model
    call buildFluxEuler2d(p_Kdiagonal, p_IverticesAtEdge,&
        p_rmatrixMass%NEQ, NVAR2D, p_rafcstabEuler%NEDGE,&
        p_MC, p_ML, p_Cx, p_Cy, p_uEuler,&
        p_fluxEuler, p_flux0Euler, p_dataEuler)

    ! Build the flux for the transport model
    call buildFluxTransport2d(p_Kdiagonal, p_IverticesAtEdge,&
        p_DcoefficientsAtEdge, p_rmatrixMass%NEQ, p_rafcstabTransport&
        %NEDGE, p_MC, p_ML, p_Cx, p_Cy, p_uTransport, p_dataTransport,&
        p_fluxTransport, p_flux0Transport)

    
    ! Build the correction coefficients
    ! - density -
    call buildCorrectionEuler(p_IverticesAtEdge, p_rmatrixMass%NEQ,&
        NVAR2D, p_rafcstabEuler%NEDGE, 1, rtimestep%dStep, p_ML, p_uEuler,&
        p_fluxEuler, p_flux0Euler, p_alpha, p_ppEuler, p_pmEuler,&
        p_qpEuler, p_qmEuler, p_rpEuler, p_rmEuler)

    ! - pressure -
    call buildCorrectionEuler(p_IverticesAtEdge, p_rmatrixMass%NEQ,&
        NVAR2D, p_rafcstabEuler%NEDGE, 4, rtimestep%dStep, p_ML, p_uEuler,&
        p_fluxEuler, p_flux0Euler, p_alpha, p_ppEuler, p_pmEuler,&
        p_qpEuler, p_qmEuler , p_rpEuler, p_rmEuler)
    
    ! - tracer -
    call buildCorrectionEuler(p_IverticesAtEdge, p_rmatrixMass%NEQ,&
        1, p_rafcstabTransport%NEDGE, 1, rtimestep%dStep, p_ML, p_uTransport,&
        p_fluxTransport, p_flux0Transport, p_alpha, p_ppTransport,&
        p_pmTransport, p_qpTransport, p_qmTransport, p_rpTransport,&
        p_rmTransport)
    
    
    ! Apply correction to low-order solution of the Euler model
    call applyCorrectionEuler(p_IverticesAtEdge, p_rmatrixMass%NEQ,&
        NVAR2D, p_rafcstabEuler%NEDGE, rtimestep%dStep, p_ML,&
        p_fluxEuler, p_alpha, p_uEuler, p_ppEuler, p_pmEuler, p_dataEuler)

    ! Set boundary conditions explicitly
    call bdrf_filterVectorExplicit(rbdrCondEuler, rsolutionEuler,&
        rtimestep%dTime, euler_calcBoundaryvalues2d)
    
    ! Release flux vectors
    call lsysbl_releaseVector(rdataEuler)

    
    ! Apply correction to low-order solution of the transport model
    call applyCorrectionTransport(p_IverticesAtEdge, p_rmatrixMass&
        %NEQ, p_rafcstabTransport%NEDGE, rtimestep%dStep, p_ML,&
        p_fluxTransport, p_alpha, p_uTransport, p_dataTransport)
    
    ! Set boundary conditions explicitly
    call bdrf_filterVectorExplicit(rbdrCondTransport,&
        rsolutionTransport, rtimestep%dTime)
    
    ! Release flux vectors
    call lsysbl_releaseVector(rdataTransport)


    ! CONSISTENCY CHECK #2:
    d1 = minval(p_uTransport)
    if (d1 .le. -1e-8) then
      print *, "Tracer has become negative in the flux-correction meth&
          &od", d1
      pause
    end if

  contains
   
    !***************************************************************************
!<subroutine>    

    subroutine buildFluxEuler2d(Kdiagonal, IverticesAtEdge,&
        NEQ, NVAR, NEDGE, MC, ML, Cx, Cy, u, flux, flux0, troc)

!<input>
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEQ,NVAR,NEDGE
      real(DP), dimension(:), intent(in) :: MC,ML,Cx,Cy
      real(DP), dimension(NVAR,NEQ), intent(in) :: u
!</input>

      real(DP), dimension(NVAR,NEQ), intent(inout) :: troc

!<output>
      real(DP), dimension(NVAR,NEDGE), intent(out) :: flux0,flux
      
!</output>
!</subroutine>
      
      ! local variables
      real(DP), dimension(NVAR) :: K_ij,K_ji,D_ij,Diff,F_ij,F_ji
      real(DP), dimension(NDIM2D) :: C_ij,C_ji
      integer :: ij,ji,i,j,iedge

!!$      ! Initialize time rate of change
!!$      call lalg_clearVector(troc)
      
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
        
        ! Calculate low-order flux
        call euler_calcFluxRusanov2d(u(:,i), u(:,j),&
            C_ij, C_ji, i, j, 1.0_DP, F_ij, F_ji)
        
        ! Update the time rate of change vector
        troc(:,i) = troc(:,i) + F_ij
        troc(:,j) = troc(:,j) + F_ji
        
        ! Calculate diffusion coefficient
        call euler_calcMatrixRusanovDiag2d(u(:,i), u(:,j),&
            C_ij, C_ji, i, j, 1.0_DP, K_ij, K_ji, D_ij)
        
        ! Compute solution difference
        Diff = u(:,i)-u(:,j)         
        
        ! Compute the raw antidiffusive flux
        flux0(:,iedge) = D_ij*Diff
      end do

      ! Scale the time rate of change by the lumped mass matrix
      !$omp parallel do
      do i = 1, NEQ
        troc(:,i) = troc(:,i)/ML(i)
      end do
      !$omp end parallel do

      ! Loop over all edges
      !$omp parallel do private(i,j,ij)
      do iedge = 1, NEDGE
      
        ! Get node numbers and matrix positions
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        ij = IverticesAtEdge(3, iedge)
        
        ! Apply mass antidiffusion
        flux(:,iedge) = flux0(:,iedge) + MC(ij)*(troc(:,i)-troc(:,j))
      end do
      !$omp end parallel do

    end subroutine buildFluxEuler2d

    !***************************************************************************
!<subroutine>

    subroutine  buildCorrectionEuler(IverticesAtEdge, NEQ, NVAR, NEDGE,&
        ivar, dscale, ML, u, flux, flux0, alpha, pp, pm, qp, qm, rp, rm)

!<input>
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEQ,NVAR,NEDGE,ivar
      real(DP), intent(in) :: dscale
      real(DP), dimension(:), intent(in) :: ML
      real(DP), dimension(NVAR,NEQ), intent(in) :: u
      real(DP), dimension(NVAR,NEDGE), intent(in) :: flux0
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

      real(DP) :: u_i,v_i,u_j,v_j,p_ij,p_ji,r_i,r_j,p_i,p_j
      
      ! Initialize vectors
      call lalg_clearVector(pp)
      call lalg_clearVector(pm)
      call lalg_clearVector(qp)
      call lalg_clearVector(qm)
      call lalg_setVector(rp, 1.0_DP)
      call lalg_setVector(rm, 1.0_DP)

      if (ivar .eq. 4) then

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
            p_ij = (GAMMA-1) * ( flux(4, iedge) + 0.5 * (u_i*u_i + v_i*v_i)*flux(1, iedge) -&
                                 u_i*flux(2, iedge) - v_i*flux(3, iedge) )
            p_ji =-(GAMMA-1) * ( flux(4, iedge) + 0.5 * (u_j*u_j + v_j*v_j)*flux(1, iedge) -&
                                 u_j*flux(2, iedge) - v_j*flux(3, iedge) )
            
!!$            ! MinMod prelimiting of antidiffusive fluxes
!!$            if (p_ij * diff > 0 .or. -p_ji*diff > 0) then
!!$              alpha(iedge) = 0.0; p_ij = 0; p_ji = 0
!!$            end if

!!$            ! MinMod prelimiting of antidiffusive fluxes
!!$            if (p_ij > SYS_EPSREAL .and. p0_ij > SYS_EPSREAL) then
!!$              aux = min(p_ij, p0_ij)
!!$              alpha(iedge) = min(alpha(iedge), aux/p_ij)
!!$            elseif (p_ij < -SYS_EPSREAL .and. p0_ij < -SYS_EPSREAL) then
!!$              aux = max(p_ij, p0_ij)
!!$              alpha(iedge) = min(alpha(iedge), aux/p_ij)
!!$            else
!!$              aux = 0.0; alpha(iedge) = 0.0
!!$            end if
            
!!$            if (p_ji > SYS_EPSREAL .and. p0_ji > SYS_EPSREAL) then
!!$              aux = min(p_ji, p0_ji)
!!$              alpha(iedge) = min(alpha(iedge), aux/p_ji)
!!$            elseif (p_ji < -SYS_EPSREAL .and. p0_ji < -SYS_EPSREAL) then
!!$              aux = max(p_ji, p0_ji)
!!$              alpha(iedge) = min(alpha(iedge), aux/p_ji)
!!$            else
!!$              aux = 0.0; alpha(iedge) = 0.0
!!$            end if
            
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

      else
      
      ! Loop over all edges
      do iedge = 1, NEDGE
        
        ! Get node numbers and matrix position
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        
        ! Flux correction in conservative variables
        f_ij  = flux(ivar,iedge)
        f0_ij = flux0(ivar,iedge)
        diff  = u(ivar,j)-u(ivar,i)
        
        ! MinMod prelimiting of antidiffusive fluxes
        if (f_ij > SYS_EPSREAL .and. f0_ij > SYS_EPSREAL) then
          aux = min(f_ij, f0_ij)
          alpha(iedge) = min(alpha(iedge), aux/f_ij)
        elseif (f_ij < -SYS_EPSREAL .and. f0_ij < -SYS_EPSREAL) then
          aux = max(f_ij, f0_ij)
          alpha(iedge) = min(alpha(iedge), aux/f_ij)
        else
          aux = 0.0; alpha(iedge) = 0.0
        end if

        f_ij = aux; flux(ivar,iedge) = f_ij
        
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

      end if

      ! Compute nodal correction factors
      !$omp parallel do
      do i = 1, NEQ
        qp(i) = qp(i)*ML(i)/dscale
        qm(i) = qm(i)*ML(i)/dscale

        rp(i) = min(1._DP, qp(i)/(pp(i)+SYS_EPSREAL) )
        rm(i) = min(1._DP, qm(i)/(pm(i)-SYS_EPSREAL) )
      end do
      !$omp end parallel do


      if (ivar .eq. 4) then

        ! Loop over all edges
        !$omp parallel do private(i,j,f_ij)
        do iedge = 1, NEDGE
          
          ! Get node numbers
          i  = IverticesAtEdge(1, iedge)
          j  = IverticesAtEdge(2, iedge)
          
          ! Velocities
          u_i = u(2,i)/u(1,i);   v_i = u(3,i)/u(1,i)
          u_j = u(2,j)/u(1,j);   v_j = u(3,j)/u(1,j)
          
          ! Pressure variables
          p_ij = (GAMMA-1) * ( flux(4, iedge) + 0.5 * (u_i*u_i + v_i*v_i)*flux(1, iedge) -&
                               u_i*flux(2, iedge) - v_i*flux(3, iedge) )
          p_ji =-(GAMMA-1) * ( flux(4, iedge) + 0.5 * (u_j*u_j + v_j*v_j)*flux(1, iedge) -&
                               u_j*flux(2, iedge) - v_j*flux(3, iedge) )


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
        
      else

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

      end if
      
    end subroutine buildCorrectionEuler

    !***************************************************************************

!<subroutine>

    subroutine applyCorrectionEuler(IverticesAtEdge, NEQ, NVAR, NEDGE,&
        dscale, ML, flux, alpha, u, pp, pm, data)
      
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
      real(DP), dimension(:), intent(out) :: pp,pm
      real(DP), dimension(NVAR,NEQ), intent(out) :: data
!</output>
!</subroutine>   
      
      ! local variables
      real(DP), dimension(:,:), pointer :: ufail
      logical, dimension(:), pointer :: bfail
      real(DP), dimension(NVAR) :: f_ij
      real(DP) :: diff
      integer :: i,j,iedge

      ! Initialize correction
      call lalg_clearVector(data)

!!$      allocate(ufail(NVAR, NEQ), bfail(NEQ))
!!$
!!$      ! Initialize vectors
!!$      call lalg_clearVector(pp)
!!$      call lalg_clearVector(pm)
!!$      call lalg_copyVector(u, ufail)

      ! Loop over all edges
      do iedge = 1, NEDGE
      
        ! Get node numbers and matrix position
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        
        ! Limit raw antidiffusive flux
        f_ij = alpha(iedge)*flux(:,iedge)
        
        where (abs(f_ij) .le. SYS_EPSREAL) f_ij = 0.0_DP

        ! Apply correction
        data(:,i) = data(:,i) + f_ij
        data(:,j) = data(:,j) - f_ij
        
!!$        ! Compute pressure difference at nodes I and J
!!$        diff = (GAMMA-1) * ( u(4,j) - u(4,i) &
!!$                         - 0.5_DP*(u(2,j)**2 + u(3,j)**2)/u(1,j)&
!!$                         + 0.5_DP*(u(2,i)**2 + u(3,i)**2)/u(1,i))
!!$        
!!$        ! Compute upper and lower bounds for the pressure
!!$        pp(i) = max(pp(i),  diff)
!!$        pp(j) = max(pp(j), -diff)
!!$        pm(i) = min(pm(i),  diff)
!!$        pm(j) = min(pm(j), -diff)
      end do
      
      ! Loop over all rows
      !$omp parallel do private(diff)
      do i = 1, NEQ

        ! Compute flux-corrected solution
        u(:,i) = u(:,i) + dscale * data(:,i)/ML(i)

!!$        ! Compute pressure difference between the low-order solution
!!$        ! and flux-corrected solution value at node I
!!$        diff = (GAMMA-1) * ( u(4,i) - ufail(4,i) &
!!$                           - 0.5_DP*(u(2,i)**2 + u(3,i)**2)/u(1,i)&
!!$                           + 0.5_DP*(ufail(2,i)**2 + ufail(3,i)**2)/ufail(1,i))
!!$        
!!$        ! Check if pressure violates upper/lower bounds
!!$        bfail(i) = (pm(i) .gt. diff) .or. (diff .gt. pp(i))
      end do
      !$omp end parallel do

      return
      
      ! Do we need to apply failsave limiting?
      if (all(.not.bfail)) then
        deallocate(ufail, bfail)
        return
      end if
     
      ! Loop over all edges
      !$omp parallel do private(i,j)
      do iedge = 1, NEDGE
        
        ! Get node numbers and matrix position
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        
        ! Do we have to cancel the antidiffusive flux for this edge?
        if (bfail(i) .or. bfail(j)) alpha(iedge) = 0.0_DP
      end do
      !$omp end parallel do

      ! Deallocate temporal memory
      deallocate(bfail)


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
        u(:,i) = ufail(:,i) + dscale * data(:,i)/ML(i)
      end do
      !$omp end parallel do

      ! Deallocate temporal memory
      deallocate(ufail)

    end subroutine applyCorrectionEuler

    !***************************************************************************

!<subroutine>

    subroutine buildFluxTransport2d(Kdiagonal, IverticesAtEdge,&
        DcoefficientsAtEdge, NEQ, NEDGE, MC, ML, Cx, Cy,&
        u, troc, flux, flux0)

!<input>
      integer, dimension(:), intent(in) :: Kdiagonal
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      real(DP), dimension(:,:), intent(in) :: DcoefficientsAtEdge
      integer, intent(in) :: NEQ,NEDGE
      real(DP), dimension(:), intent(in) :: MC,ML,Cx,Cy,u,troc
!</input>

!<output>
      real(DP), dimension(:), intent(out) :: flux,flux0
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
        flux0(iedge) = d_ij * (u(i)-u(j))

        ! Apply mass antidiffusion
        flux(iedge) = flux0(iedge) + MC(ij)*(troc(i)-troc(j))
      end do
      !$omp end parallel do

    end subroutine buildFluxTransport2d

!!$    !***************************************************************************
!!$
!!$!<subroutine>
!!$    
!!$    subroutine buildCorrectionTransport(IverticesAtEdge, NEQ, NEDGE,&
!!$        dscale, ML, u, flux, flux0, alpha, pp, pm, qp, qm, rp, rm)
!!$
!!$!<input>
!!$      integer, dimension(:,:) :: IverticesAtEdge
!!$      integer, intent(in) :: NEQ,NEDGE
!!$      real(DP), intent(in) :: dscale
!!$      real(DP), dimension(:), intent(in) :: ML,u,flux0
!!$!</input>
!!$      
!!$!<inputoutput>
!!$      real(DP), dimension(:), intent(inout) :: alpha
!!$      real(DP), dimension(:), intent(inout) :: flux
!!$!</inputoutput>
!!$   
!!$!<output>   
!!$      real(DP), dimension(:), intent(out) :: pp,pm,qp,qm,rp,rm
!!$!</output>
!!$      
!!$      ! local variables
!!$      real(DP) ::  f_ij,f0_ij,diff,aux
!!$      integer :: ij,ji,i,j,iedge
!!$
!!$      
!!$      ! Initialize vectors
!!$      call lalg_clearVector(pp)
!!$      call lalg_clearVector(pm)
!!$      call lalg_clearVector(qp)
!!$      call lalg_clearVector(qm)
!!$      call lalg_setVector(rp, 1.0_DP)
!!$      call lalg_setVector(rm, 1.0_DP)
!!$      
!!$      ! Loop over all edges
!!$      do iedge = 1, NEDGE
!!$      
!!$        ! Get node numbers
!!$        i  = IverticesAtEdge(1, iedge)
!!$        j  = IverticesAtEdge(2, iedge)
!!$        
!!$        ! Flux correction in conservative variables
!!$        f_ij  = flux(iedge)
!!$        f0_ij = flux0(iedge)
!!$        diff  = u(j)-u(i)
!!$
!!$        ! MinMod prelimiting of antidiffusive fluxes
!!$        if (f_ij > SYS_EPSREAL .and. f0_ij > SYS_EPSREAL) then
!!$          aux = min(f_ij, f0_ij)
!!$          alpha(iedge) = min(alpha(iedge), aux/f_ij)
!!$        elseif (f_ij < -SYS_EPSREAL .and. f0_ij < -SYS_EPSREAL) then
!!$          aux = max(f_ij, f0_ij)
!!$          alpha(iedge) = min(alpha(iedge), aux/f_ij)
!!$        else
!!$          f_ij = 0.0; alpha(iedge) = 0.0
!!$        end if

!!$
!!$
!!$        if (f_ij * f0_ij < 0) then
!!$          f_ij = 0.0; alpha(iedge) = 0.0
!!$        elseif (f_ij > 0) then
!!$          f_ij = min(f_ij, f0_ij)
!!$        else
!!$          f_ij = max(f_ij, f0_ij)
!!$        end if

!!$        if (f_ij > SYS_EPSREAL .and. f0_ij > SYS_EPSREAL) then
!!$          f_ij = min(f_ij, f0_ij)
!!$        elseif (f_ij < -SYS_EPSREAL .and. f0_ij < -SYS_EPSREAL) then
!!$          f_ij = max(f_ij, f0_ij)
!!$        else
!!$          f_ij = 0.0
!!$        end if
!!$        
!!$        ! Store flux
!!$        flux(iedge) = f_ij
        
!!$        ! Sums of raw antidiffusive fluxes
!!$        pp(i) = pp(i) + max(0.0_DP,  f_ij)
!!$        pp(j) = pp(j) + max(0.0_DP, -f_ij)
!!$        pm(i) = pm(i) + min(0.0_DP,  f_ij)
!!$        pm(j) = pm(j) + min(0.0_DP, -f_ij)
!!$        
!!$        ! Sums of admissible edge contributions
!!$        qp(i) = max(qp(i),  diff)
!!$        qp(j) = max(qp(j), -diff)
!!$        qm(i) = min(qm(i),  diff)
!!$        qm(j) = min(qm(j), -diff)
!!$      end do
!!$      
!!$      
!!$      ! Compute nodal correction factors
!!$      !$omp parallel do
!!$      do i = 1, NEQ
!!$        qp(i) = qp(i)*ML(i)/dscale
!!$        qm(i) = qm(i)*ML(i)/dscale
!!$        
!!$        if (pp(i) > qp(i) + SYS_EPSREAL) rp(i) = qp(i)/pp(i)
!!$        if (pm(i) < qm(i) - SYS_EPSREAL) rm(i) = qm(i)/pm(i)

!!$
!!$        rp(i) = min(1._DP, qp(i)/(pp(i)+SYS_EPSREAL) )
!!$        rm(i) = min(1._DP, qm(i)/(pm(i)-SYS_EPSREAL) )
!!$        
!!$      end do
!!$      !$omp end parallel do
      
!!$      print *, "Checking |RP|<=|Q|"
!!$      do i = 1, NEQ
!!$        if (abs(rp(i)*pp(i)) .gt. abs(qp(i)))&
!!$            print *, "+", i, abs(rp(i)*pp(i))/abs(qp(i))
!!$        if (abs(rm(i)*pm(i)) .gt. abs(qm(i)))&
!!$            print *, "-", i, abs(rm(i)*pm(i))/abs(qm(i))
!!$      end do
!!$      pause

!!$      ! Loop over all edges
!!$      !$omp parallel do private(i,j,f_ij)
!!$      do iedge = 1, NEDGE
!!$        
!!$        ! Get node numbers
!!$        i  = IverticesAtEdge(1, iedge)
!!$        j  = IverticesAtEdge(2, iedge)
!!$        
!!$        ! Flux correction
!!$        f_ij = flux(iedge)
!!$        
!!$        ! Limit fluxes
!!$        if (f_ij > SYS_EPSREAL) then
!!$          alpha(iedge) = min(alpha(iedge), rp(i), rm(j))
!!$        elseif (f_ij < -SYS_EPSREAL) then
!!$          alpha(iedge) = min(alpha(iedge), rm(i), rp(j))
!!$        end if
!!$      end do
!!$      !$omp end parallel do
!!$
!!$    end subroutine buildCorrectionTransport

    !***************************************************************************
    
!<subroutine>

    subroutine applyCorrectionTransport(IverticesAtEdge, NEQ, NEDGE,&
        dscale, ML, flux, alpha, u, data)

!<input>
      integer, dimension(:,:), intent(in) :: IverticesAtEdge
      integer, intent(in) :: NEQ,NEDGE
      real(DP), intent(in) :: dscale
      real(DP), dimension(:), intent(in) :: ML,flux, alpha
!</input>

!<inputoutput>     
      real(DP), dimension(:), intent(inout) :: u
!</inputoutput>

!<output>
      real(DP), dimension(:), intent(out) :: data
!</output>
!<subroutine>

      ! local variables
      real(DP) :: f_ij,diff
      integer :: ij,ji,i,j,iedge
      
      ! Initialize correction
      call lalg_clearVector(data)
      
      ! Loop over all edges
      do iedge = 1, NEDGE
      
        ! Get node numbers
        i  = IverticesAtEdge(1, iedge)
        j  = IverticesAtEdge(2, iedge)
        
        ! Limit conservative fluxes
        f_ij = alpha(iedge)*flux(iedge)
        
        ! Apply correction
        data(i) = data(i) + f_ij
        data(j) = data(j) - f_ij
      end do

      ! Loop over all rows
      !$omp parallel do
      do i = 1, NEQ
        
        if (dscale * data(i)/ML(i) + u(i) .lt. -1e-3) then
          print *, dscale*data(i)/ML(i), u(i), dscale * data(i)/ML(i) + u(i)
          pause
        end if

        ! Compute flux-corrected solution
        u(i) = u(i) + dscale * data(i)/ML(i)
      end do
      !$omp end parallel do
      
    end subroutine applyCorrectionTransport
    
  end subroutine zpinch_calcLinearizedFCT

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

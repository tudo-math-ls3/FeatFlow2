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
        
        ! Set pointer to parameter list
        p_rparlist => collct_getvalue_parlst(rcollection, 'rparlist')

        ! Set pointers to current solution vectors stored in the
        ! first and second quick access vector
        p_rsolutionEuler     => rcollection%p_rvectorQuickAccess1
        p_rsolutionTransport => rcollection%p_rvectorQuickAccess2
        
        ! Calculate scaling for implicit part of the Lorentz force
        dscale = -rtimestep%theta * rtimestep%dStep

        if (dscale .ne. 0.0_DP) then
          
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
            zpinch_calcMatrixPrimalConst2d,&
            zpinch_calcMatrixDualConst2d,&
            transp_coeffMatBdrPrimalConst2d,&
            transp_coeffMatBdrDualConst2d)

        ! Compute the right-hand side
        call transp_calcRhsRungeKuttaScheme(rproblemLevel,&
            rtimestep, rsolver, rsolution, rsolution0,&
            rrhs, istep, rcollection,&
            fcb_coeffVecBdrPrimal_sim=transp_coeffVecBdrPrimalConst2d,&
            fcb_coeffVecBdrDual_sim=transp_coeffVecBdrDualConst2d)
        
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
                zpinch_calcMatrixPrimalConst2d,&
                zpinch_calcMatrixDualConst2d,&
                transp_coeffMatBdrPrimalConst2d,&
                transp_coeffMatBdrDualConst2d)
            
            ! Assemble the constant right-hand side
            call transp_calcRhsThetaScheme(rproblemLevel, rtimestep,&
                rsolver, rsolution0, rrhs, rcollection,&
                fcb_coeffVecBdrPrimal_sim=transp_coeffVecBdrPrimalConst2d,&
                fcb_coeffVecBdrDual_sim=transp_coeffVecBdrDualConst2d)
            
          end if

          ! Set pointer to the solution vector of the Euler model
          ! from the current time step
          p_rsolutionEuler => rcollection%p_rvectorQuickAccess1
          
          ! Calculate the density averaged mass matrices using the
          ! solution of the Euler model from the current time step
          call zpinch_initDensityAveraging(p_rparlist,&
              rcollection%SquickAccess(1), rproblemLevel,&
              p_rsolutionEuler, rcollection)
          
          ! Calculate the velocity vector using the solution of the
          ! Euler model from the current time step
          call zpinch_initVelocityField(p_rparlist,&
              rcollection%SquickAccess(1), rproblemLevel,&
              p_rsolutionEuler, rcollection)

          ! Compute the preconditioner
          call transp_calcPrecondThetaScheme(rproblemLevel,&
              rtimestep, rsolver, rsolution, rcollection,&
              zpinch_calcMatrixPrimalConst2d,&
              zpinch_calcMatrixDualConst2d,&
              transp_coeffMatBdrPrimalConst2d,&
              transp_coeffMatBdrDualConst2d)

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
            zpinch_calcMatrixPrimalConst2d,&
            zpinch_calcMatrixDualConst2d,&
            transp_coeffMatBdrPrimalConst2d,&
            transp_coeffMatBdrDualConst2d)
      end if
      
      
      ! Do we have to calculate the Jacobian operator?
      ! --------------------------------------------------------------------------
      if (iand(iSpec, NLSOL_OPSPEC_CALCJACOBIAN) .ne. 0) then
        
        ! Compute the Jacobian matrix
        call transp_calcJacobianThetaScheme(rproblemLevel,&
            rtimestep, rsolver, rsolution, rsolution0, rcollection,&
            zpinch_calcMatrixPrimalConst2d,&
            zpinch_calcMatrixDualConst2d,&
            transp_coeffMatBdrPrimalConst2d,&
            transp_coeffMatBdrDualConst2d)
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
    
    ! Set x-velocity, i.e., momentum in x-direction
    call euler_getVariable(rsolution, 'momentum_x',&
        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(1))
    call zpinch_setVariable2d(&
        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(1), 1)
    
    ! Set y-velocity, i.e., momentum in y-direction
    call euler_getVariable(rsolution, 'momentum_y',&
        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(2))
    call zpinch_setVariable2d(&
        rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(2), 2)
    
    ! Set update notification in problem level structure
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     PROBLEV_MSPEC_UPDATE)

  end subroutine zpinch_initVelocityField

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_initDensityAveraging(rparlist,&
      ssectionNameTransport, rproblemlevel, rsolutionEuler,&
      rcollection)

!<description>
    ! This subroutine initializes the density averaged mass matrices
    ! for the transport model based on the solution from the Euler
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section names in parameter list
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
    type(t_trilinearform) :: rform
    type(t_vectorScalar) :: rvector
    integer :: lumpedMassMatrix, consistentMassMatrix


    ! Get global configuration from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionNameTransport, 'lumpedmassmatrix', lumpedMassMatrix)
    call parlst_getvalue_int(rparlist,&
        ssectionNameTransport, 'consistentmassmatrix', consistentMassMatrix)
    
    ! Get density distribution from the solution of the Euler model
    ! and create block vector which is attached to the collection
    call euler_getVariable(rsolutionEuler, 'density', rvector)
    
    ! We have variable coefficients
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff    = .true.

    ! Initialize the bilinear form
    rform%itermCount = 1
    rform%Dcoefficients(1)  = 1.0_DP
    rform%Idescriptors(1,1) = DER_FUNC
    rform%Idescriptors(2,1) = DER_FUNC
    rform%Idescriptors(3,1) = DER_FUNC

    ! Create density averaged consistent mass matrix
    call trilf_buildMatrixScalar(rform, .true.,&
        rproblemLevel%Rmatrix(consistentMassMatrix), rvector)
    
    ! Create density averaged lumped mass matrix
    call lsyssc_duplicateMatrix(rproblemLevel&
        %Rmatrix(consistentMassMatrix), rproblemLevel &
        %Rmatrix(lumpedMassMatrix), LSYSSC_DUP_SHARE, LSYSSC_DUP_COPY)
    
    call lsyssc_lumpMatrixScalar(rproblemLevel&
        %Rmatrix(lumpedMassMatrix), LSYSSC_LUMP_DIAG)
    
    ! Release temporal vector
    call lsyssc_releaseVector(rvector)
    
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
    call lsyssc_getbase_double(rproblemLevel&
        %Rmatrix(lumpedMassMatrix), p_DlumpedMassMatrix)
    
    ! Set pointer to global solution vectors
    call lsysbl_getbase_double(rforce, p_Ddata)
    call lsysbl_getbase_double(rsolutionEuler, p_DdataEuler)
    call lsysbl_getbase_double(rsolutionTransport, p_DdataTransport)
    
    ! Set pointer to the vertex coordinates
    call storage_getbase_double2D(rproblemLevel%rtriangulation&
        %h_DvertexCoords, p_DvertexCoords)
    
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
      do i = 1, neq
        
        ! Get coordinates at node i
        x1 = DvertexCoords(1, i)
        x2 = DvertexCoords(2, i)
        
        ! Compute polar coordinates
        drad = sqrt(x1*x1 + x2*x2)
        dang = atan2(x2, x1)
        
        ! Compute unit vector into origin
        if (drad .gt. 1e-4) then
          x1 = cos(dang)
          x2 = sin(dang)
        else
          x1 = 0.0; x2 = 0.0
        end if
        
        ! Compute source term
        daux1 = dconst * DmassMatrix(i) * DdataTransport(i) *&
                         DdataEuler(1, i) / max(drad, 1.0e-4_DP)
        daux2 = dconst * DmassMatrix(i) * DdataTransport(i) *&
                         (DdataEuler(2, i) * x1 + DdataEuler(3,i) * x2) /&
                         max(drad, 1.0e-4_DP)
               
        ! Impose source values
        DdataForce(1, i) = 0.0_DP
        DdataForce(2, i) = daux1 * x1
        DdataForce(3, i) = daux1 * x2
        DdataForce(4, i) = daux2
      end do
      
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
        
        ! Compute source term
        daux1 = dconst * DmassMatrix(i) * DdataTransport(i) *&
                         DdataEuler(1, i) / max(drad, 1.0e-4_DP)
        daux2 = dconst * DmassMatrix(i) * DdataTransport(i) *&
                         (DdataEuler(2, i) * x1 + DdataEuler(3,i) * x2) /&
                         max(drad, 1.0e-4_DP)
        
        ! Impose source values
        DdataForce(1, i) = 0.0_DP
        DdataForce(2, i) = daux1 * x1
        DdataForce(3, i) = daux1 * x2
        DdataForce(4, i) = daux2
      end do

    end subroutine calcForceRZInterleaveFormat

  end subroutine zpinch_initLorentzforceTerm

  !*****************************************************************************
  
  !<subroutine>

  subroutine zpinch_calcLinearizedFCT(rbdrCond, rproblemLevel,&
      rtimestep, rsolution, rcollection)

!<description>
    ! This subroutine calculates the linearized FCT correction
!</description>

!<input>
    ! boundary condition structure
    type(t_boundaryCondition), intent(in) :: rbdrCond

    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep    
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solution vector
    type(t_vectorBlock), intent(inout) :: rsolution

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixScalar), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rsolutionEuler
    type(t_parlist), pointer :: p_rparlist
    type(t_vectorScalar) :: rflux0, rflux
    type(t_vectorBlock) :: rdata
    real(DP), dimension(:), pointer :: p_MC, p_ML, p_Cx, p_Cy
    real(DP), dimension(:), pointer :: p_u, p_flux0, p_flux, p_data
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal, p_Ksep
    integer :: h_Ksep, templatematrix, lumpedMassMatrix, consistentMassMatrix
    integer :: coeffMatrix_CX, coeffMatrix_CY
    integer :: i, nedge

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

    ! Multiply the low-order solution by the density-averaged mass
    ! matrix evaluated at the low-order density profile
    do i = 1, size(p_u)
      p_u(i) = p_ML(i) * p_u(i)
    end do

    ! Set pointer to solution from Euler model
    p_rsolutionEuler => rcollection%p_rvectorQuickAccess1

    ! Calculate the density averaged mass matrix for the new density
    call zpinch_initDensityAveraging(p_rparlist,&
        rcollection%SquickAccess(1), rproblemLevel,&
        p_rsolutionEuler, rcollection)
    
    ! Calculate the new velocity field based on the new momentum
    call zpinch_initVelocityField(p_rparlist,&
        rcollection%SquickAccess(1), rproblemLevel,&
        p_rsolutionEuler, rcollection)
    
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
        call zpinch_calcMatrixPrimalConst2d(u(i), u(i),&
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
          call zpinch_calcMatrixPrimalConst2d(u(i), u(j),&
              C_ij, C_ji, i, j, k_ij, k_ji, d_ij)
          
          ! Artificial diffusion coefficient
          d_ij = max(abs(k_ij), abs(k_ji))

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
      call lalg_setVector(rp, 1.0_DP)
      call lalg_setVector(rm, 1.0_DP)
      
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
        u(i) = (u(i) + data(i))/ML(i)
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
    
  end subroutine zpinch_calcLinearizedFCT

end module zpinch_callback

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

#include "hydro.h"

  use afcstabilisation
  use basicgeometry
  use boundarycondaux
  use boundaryfilter
  use collection
  use derivatives
  use hydro_basic
  use hydro_callback
  use hydro_callback1d
  use hydro_callback2d
  use hydro_callback3d
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

  !*****************************************************************************

!<constants>

!<constantblock>
  ! Minimum number of equations for OpenMP parallelisation: If the number of
  ! equations is below this value, then no parallelisation is performed.
#ifndef ZPINCH_LORENTZFORCE_NEQMIN_OMP
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: ZPINCH_LORENTZFORCE_NEQMIN_OMP = 10000
#else
  integer, public            :: ZPINCH_LORENTZFORCE_NEQMIN_OMP = 10000
#endif
#endif
  
!</constantblock>

!</constants>

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
    real(DP), save :: dtimeHydro = 0.0_DP
    real(DP), save :: dtimeTransport = 0.0_DP

    ! local variables
    type(t_parlist), pointer :: p_rparlist
    type(t_vectorBlock), pointer :: p_rsolutionHydro
    type(t_vectorBlock), pointer :: p_rsolutionTransport
    type(t_vectorBlock) :: rforce
    character(LEN=SYS_STRLEN) :: ssectionName
    character(LEN=SYS_STRLEN) :: ssectionNameHydro
    character(LEN=SYS_STRLEN) :: ssectionNameTransport
    real(DP) :: dscale
    integer(i32) :: iSpec
    integer :: isystemFormat, jacobianMatrix

    ! Get section names
    call collct_getvalue_string(rcollection,&
        'ssectionname', ssectionName)
    call collct_getvalue_string(rcollection,&
        'ssectionnameHydro', ssectionNameHydro)
    call collct_getvalue_string(rcollection,&
        'ssectionnameTransport', ssectionNameTransport)


    if (trim(rsolver%ssolverName) .eq. 'NonlinearSolverHydro') then
      
      ! Do we have to calculate the preconditioner?
      ! --------------------------------------------------------------------------
      if (iand(ioperationSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) then
        
        ! Compute the preconditioner
        call hydro_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
            rsolver, rsolution, ssectionNameHydro, rcollection)
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

        if ((istep .eq. 0) .and.&
            (dtimeHydro .ne. rtimestep%dTime)) then

          ! Update time variable so that no re-evaluation of the
          ! const right-hand side takes place for current time step
          dtimeHydro = rtimestep%dTime

          ! Compute the constant right-hand side including the
          ! given explicit part of the Lorentz force term
          call hydro_calcRhsThetaScheme(rproblemLevel, rtimestep,&
              rsolver, rsolution0, rrhs, ssectionNameHydro, rcollection,&
              rsource)
        end if

        ! Set pointers to current solution vectors stored in the
        ! first and second quick access vector
        p_rsolutionHydro     => rcollection%p_rvectorQuickAccess1
        p_rsolutionTransport => rcollection%p_rvectorQuickAccess2
        
        ! Set pointer to parameter list (we can always use the
        ! parameter list from the same section since all three
        ! sections in the collection have the same data)
        p_rparlist => collct_getvalue_parlst(rcollection,&
            'rparlist', ssectionName=ssectionName)
        
        ! Calculate scaling for implicit parts
        dscale = -rtimestep%theta * rtimestep%dStep

        if (dscale .ne. 0.0_DP) then

          ! Compute the implicit part of the Lorentz force
          call zpinch_calcLorentzforceTerm(p_rparlist, ssectionName,&
              ssectionNameHydro, ssectionNameTransport, rproblemLevel,&
              p_rsolutionHydro, p_rsolutionTransport, rtimestep%dTime,&
              dscale, .true., rforce, rcollection)
          
          ! Compute the residual including the implicit parts
          call hydro_calcResidualThetaScheme(rproblemLevel, rtimestep,&
              rsolver, rsolution, rsolution0, rrhs, rres, istep,&
              ssectionNameHydro, rcollection, rforce)

          ! Release temporal memory
          call lsysbl_releaseVector(rforce)
          
        else
          
          ! Compute the residual without implicit parts
          call hydro_calcResidualThetaScheme(rproblemLevel, rtimestep,&
              rsolver, rsolution, rsolution0, rrhs, rres, istep,&
              ssectionNameHydro, rcollection)

        end if
      end if


      ! Do we have to impose boundary conditions?
      ! --------------------------------------------------------------------------
      if (iand(ioperationSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

        ! Impose boundary conditions
        call hydro_setBoundaryCondition(rproblemLevel, rtimestep,&
            rsolver, rsolution, rsolution0, rres, ssectionNameHydro, rcollection)
      end if

      
      ! Set status flag
      istatus = 0


    elseif (trim(rsolver%ssolverName) .eq. 'NonlinearSolverTransport') then   

      ! Set pointer to parameter list (we can always use the parameter
      ! list from the same section since all three sections in the
      ! collection have the same data)
      p_rparlist => collct_getvalue_parlst(rcollection,&
          'rparlist', ssectionName=ssectionName)

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
        
        ! Get configuration from hydrodynamic section of parameter list 
        call parlst_getvalue_int(p_rparlist,&
            ssectionNameHydro, 'isystemformat', isystemFormat)
        
        ! What type of system format are we?
        select case(isystemFormat)
          
        case (SYSTEM_INTERLEAVEFORMAT)
          
          ! Compute the preconditioner in interleaved format
          call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
              rsolver, rsolution, ssectionNameTransport, rcollection,&
              zpinch_calcMatDiagConvIntlP2d_sim,&
              zpinch_calcMatRusConvIntlP2d_sim,&
              zpinch_calcMatDiagConvIntlD2d_sim,&
              zpinch_calcMatRusConvIntlD2d_sim,&
              fcb_coeffMatBdrPrimal2d_sim=transp_coeffMatBdrConvP2d_sim,&
              fcb_coeffMatBdrDual2d_sim=transp_coeffMatBdrConvD2d_sim)

        case (SYSTEM_BLOCKFORMAT)

          ! Compute the preconditioner in block format
          call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
              rsolver, rsolution, ssectionNameTransport, rcollection,&
              zpinch_calcMatDiagConvBlockP2d_sim,&
              zpinch_calcMatRusConvBlockP2d_sim,&
              zpinch_calcMatDiagConvBlockD2d_sim,&
              zpinch_calcMatRusConvBlockD2d_sim,&
              fcb_coeffMatBdrPrimal2d_sim=transp_coeffMatBdrConvP2d_sim,&
              fcb_coeffMatBdrDual2d_sim=transp_coeffMatBdrConvD2d_sim)

        case default
          call output_line('Invalid system format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'zpinch_nlsolverCallback')
          call sys_halt()
        end select
        
        ! Compute the right-hand side
        call transp_calcRhsRungeKuttaScheme(rproblemLevel, rtimestep,&
            rsolver, rsolution, rsolution0, rrhs, istep,&
            ssectionNameTransport, rcollection, rsource,&
            fcb_coeffVecBdrPrimal2d_sim=transp_coeffVecBdrConvP2d_sim,&
            fcb_coeffVecBdrDual2d_sim=transp_coeffVecBdrConvD2d_sim)
        
        ! Remove specifier for the preconditioner (if any)
        iSpec = iand(iSpec, not(NLSOL_OPSPEC_CALCPRECOND))
      end if


      ! Do we have to calculate the residual?
      ! --------------------------------------------------------------------------
      if (iand(iSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

        ! Get configuration from hydrodynamic section of parameter list 
        call parlst_getvalue_int(p_rparlist,&
            ssectionNameHydro, 'isystemformat', isystemFormat)

        if (istep .eq. 0) then
          
          if (dtimeTransport .ne. rtimestep%dTime) then

            ! Update time variable so that no re-evaluation of the
            ! const right-hand side takes place for current time step
            dtimeTransport = rtimestep%dTime

            ! What type of system format are we?
            select case(isystemFormat)
              
            case (SYSTEM_INTERLEAVEFORMAT)
              
              ! Compute the preconditioner in interleaved format
              call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
                  rsolver, rsolution0, ssectionNameTransport, rcollection,&
                  zpinch_calcMatDiagConvIntlP2d_sim,&
                  zpinch_calcMatRusConvIntlP2d_sim,&
                  zpinch_calcMatDiagConvIntlD2d_sim,&
                  zpinch_calcMatRusConvIntlD2d_sim,&
                  fcb_coeffMatBdrPrimal2d_sim=transp_coeffMatBdrConvP2d_sim,&
                  fcb_coeffMatBdrDual2d_sim=transp_coeffMatBdrConvD2d_sim)

            case (SYSTEM_BLOCKFORMAT)

              ! Compute the preconditioner in block format
              call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
                  rsolver, rsolution0, ssectionNameTransport, rcollection,&
                  zpinch_calcMatDiagConvBlockP2d_sim,&
                  zpinch_calcMatRusConvBlockP2d_sim,&
                  zpinch_calcMatDiagConvBlockD2d_sim,&
                  zpinch_calcMatRusConvBlockD2d_sim,&
                  fcb_coeffMatBdrPrimal2d_sim=transp_coeffMatBdrConvP2d_sim,&
                  fcb_coeffMatBdrDual2d_sim=transp_coeffMatBdrConvD2d_sim)

            case default
              call output_line('Invalid system format!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'zpinch_nlsolverCallback')
              call sys_halt()
            end select
            
            ! Assemble the constant right-hand side
            call transp_calcRhsThetaScheme(rproblemLevel, rtimestep, rsolver,&
                rsolution0, rrhs, ssectionNameTransport, rcollection, rsource,&
                fcb_coeffVecBdrPrimal2d_sim=transp_coeffVecBdrConvP2d_sim,&
                fcb_coeffVecBdrDual2d_sim=transp_coeffVecBdrConvD2d_sim)

          end if

          ! Set pointers to current solution vector of the
          ! hydrodynamic model stored in the first quick access vector
          p_rsolutionHydro => rcollection%p_rvectorQuickAccess1
          
          ! Calculate the velocity vector using the solution of the
          ! Hydrodynamic model from the current time step
          call zpinch_calcVelocityField(p_rparlist, ssectionNameTransport,&
              rproblemLevel, p_rsolutionHydro, rcollection)

          ! What type of system format are we?
          select case(isystemFormat)
            
          case (SYSTEM_INTERLEAVEFORMAT)
            
            ! Compute the preconditioner in interleaved format
            call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
                rsolver, rsolution, ssectionNameTransport, rcollection,&
                zpinch_calcMatDiagConvIntlP2d_sim,&
                zpinch_calcMatRusConvIntlP2d_sim,&
                zpinch_calcMatDiagConvIntlD2d_sim,&
                zpinch_calcMatRusConvIntlD2d_sim,&
                fcb_coeffMatBdrPrimal2d_sim=transp_coeffMatBdrConvP2d_sim,&
                fcb_coeffMatBdrDual2d_sim=transp_coeffMatBdrConvD2d_sim)

          case (SYSTEM_BLOCKFORMAT)
            
            ! Compute the preconditioner in block format
            call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
                rsolver, rsolution, ssectionNameTransport, rcollection,&
                zpinch_calcMatDiagConvBlockP2d_sim,&
                zpinch_calcMatRusConvBlockP2d_sim,&
                zpinch_calcMatDiagConvBlockD2d_sim,&
                zpinch_calcMatRusConvBlockD2d_sim,&
                fcb_coeffMatBdrPrimal2d_sim=transp_coeffMatBdrConvP2d_sim,&
                fcb_coeffMatBdrDual2d_sim=transp_coeffMatBdrConvD2d_sim)

          case default
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
            rrhs, rres, istep, ssectionNameTransport, rcollection,&
            fcb_coeffVecBdrPrimal2d_sim=transp_coeffVecBdrConvP2d_sim,&
            fcb_coeffVecBdrDual2d_sim=transp_coeffVecBdrConvD2d_sim)
      end if


      ! Do we have to calculate the preconditioner?
      ! --------------------------------------------------------------------------
      if (iand(iSpec, NLSOL_OPSPEC_CALCPRECOND) .ne. 0) then
        
        ! Get configuration from hydrodynamic section of parameter list 
        call parlst_getvalue_int(p_rparlist,&
            ssectionNameHydro, 'isystemformat', isystemFormat)
        
        ! What type of system format are we?
        select case(isystemFormat)
          
        case (SYSTEM_INTERLEAVEFORMAT)
          
          ! Compute the preconditioner in interleaved format
          call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
              rsolver, rsolution, ssectionNameTransport, rcollection,&
              zpinch_calcMatDiagConvIntlP2d_sim,&
              zpinch_calcMatRusConvIntlP2d_sim,&
              zpinch_calcMatDiagConvIntlD2d_sim,&
              zpinch_calcMatRusConvIntlD2d_sim,&
              fcb_coeffMatBdrPrimal2d_sim=transp_coeffMatBdrConvP2d_sim,&
              fcb_coeffMatBdrDual2d_sim=transp_coeffMatBdrConvD2d_sim)

        case (SYSTEM_BLOCKFORMAT)
          
          ! Compute the preconditioner in block format
          call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
              rsolver, rsolution, ssectionNameTransport, rcollection,&
              zpinch_calcMatDiagConvBlockP2d_sim,&
              zpinch_calcMatRusConvBlockP2d_sim,&
              zpinch_calcMatDiagConvBlockD2d_sim,&
              zpinch_calcMatRusConvBlockD2d_sim,&
              fcb_coeffMatBdrPrimal2d_sim=transp_coeffMatBdrConvP2d_sim,&
              fcb_coeffMatBdrDual2d_sim=transp_coeffMatBdrConvD2d_sim)

        case default
          call output_line('Invalid system format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'zpinch_nlsolverCallback')
          call sys_halt()
        end select
      end if


      ! Do we have to calculate the Jacobian operator?
      ! --------------------------------------------------------------------------
      if (iand(iSpec, NLSOL_OPSPEC_CALCJACOBIAN) .ne. 0) then

        ! Get configuration from hydrodynamic section of parameter list 
        call parlst_getvalue_int(p_rparlist,&
            ssectionNameHydro, 'isystemformat', isystemFormat)
        
        ! What type of system format are we?
        select case(isystemFormat)
          
        case (SYSTEM_INTERLEAVEFORMAT)
          
          ! Compute the Jacobian matrix in interleaved format
          call transp_calcJacobianThetaScheme(rproblemLevel, rtimestep,&
              rsolver, rsolution, rsolution0, ssectionNameTransport, rcollection,&
              zpinch_calcMatRusConvIntlP2d_sim,&
              zpinch_calcMatRusConvIntlD2d_sim,&
              fcb_coeffMatBdrPrimal2d_sim=transp_coeffMatBdrConvP2d_sim,&
              fcb_coeffMatBdrDual2d_sim=transp_coeffMatBdrConvD2d_sim)

        case (SYSTEM_BLOCKFORMAT)

          ! Compute the Jacobian matrix in block format
          call transp_calcJacobianThetaScheme(rproblemLevel, rtimestep,&
              rsolver, rsolution, rsolution0, ssectionNameTransport, rcollection,&
              zpinch_calcMatRusConvBlockP2d_sim,&
              zpinch_calcMatRusConvBlockD2d_sim,&
              fcb_coeffMatBdrPrimal2d_sim=transp_coeffMatBdrConvP2d_sim,&
              fcb_coeffMatBdrDual2d_sim=transp_coeffMatBdrConvD2d_sim)

        case default
          call output_line('Invalid system format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'zpinch_nlsolverCallback')
          call sys_halt()
        end select
      end if


      ! Do we have to impose boundary conditions?
      ! --------------------------------------------------------------------------
      if (iand(iSpec, NLSOL_OPSPEC_CALCRESIDUAL) .ne. 0) then

        ! Impose boundary conditions
        call transp_setBoundaryCondition(rproblemLevel, rtimestep,&
            rsolver, rsolution, rsolution0, rres, ssectionNameTransport, rcollection)
      end if


      ! Do we have to apply the Jacobian operator?
      ! --------------------------------------------------------------------------
      if (iand(iSpec, NLSOL_OPSPEC_APPLYJACOBIAN) .ne. 0) then

        ! Get position of Jacobian matrix from transport section of parameter list
        call parlst_getvalue_int(p_rparlist,&
            ssectionNameTransport, 'jacobianMatrix', jacobianMatrix)

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
    ! of the compressible hydrodynamic model. The result is stored separately
    ! for each problem level.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! solution vector of compressible hydrodynamic model
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

    ! What dimension are we?
    select case(ndim)
    case (NDIM1D)
      ! Set x-velocity
      call hydro_getVariable(rsolution, 'velocity_x',&
          rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(1))

    case (NDIM2D)
      ! Set x-velocity
      call hydro_getVariable(rsolution, 'velocity_x',&
          rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(1))

      ! Set y-velocity
      call hydro_getVariable(rsolution, 'velocity_y',&
          rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(2))

    case (NDIM3D)
      ! Set x-velocity
      call hydro_getVariable(rsolution, 'velocity_x',&
          rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(1))
      
      ! Set y-velocity
      call hydro_getVariable(rsolution, 'velocity_y',&
          rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(2))

      ! Set z-velocity
      call hydro_getVariable(rsolution, 'velocity_z',&
          rproblemLevel%RvectorBlock(velocityfield)%RvectorBlock(3))

    case default
      call output_line('Invalid spatial dimension!',&
          OU_CLASS_ERROR, OU_MODE_STD, 'zpinch_calcLorentzforceTerm')
      call sys_halt()
    end select

    ! Set update notification in problem level structure
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     PROBLEV_MSPEC_UPDATE)

  end subroutine zpinch_calcVelocityField

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_calcLorentzforceTerm(rparlist, ssectionName,&
      ssectionNameHydro, ssectionNameTransport, rproblemLevel,&
      rsolutionHydro, rsolutionTransport, dtime, dscale, bclear,&
      rforce, rcollection)

!<description>
    ! This subroutine evaluates the Lorentz force term based on the
    ! given solution vectors and stores it in vector rforce.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section names in parameter list
    character(LEN=*), intent(in) :: ssectionName
    character(LEN=*), intent(in) :: ssectionNameHydro
    character(LEN=*), intent(in) :: ssectionNameTransport

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! solution vector for transport model
    type(t_vectorBlock), intent(in) :: rsolutionTransport

    ! solution vector for hydrodynamic model
    type(t_vectorBlock), intent(in) :: rsolutionHydro

    ! simulation time
    real(DP), intent(in) :: dtime

    ! scaling parameter
    real(DP), intent(in) :: dscale

    ! Clear force vector?
    logical, intent(in) :: bclear
!</input>

!<inputoutput>
    ! source vectors to be assembled
    type(t_vectorBlock), intent(inout) :: rforce

    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    real(DP), dimension(:), pointer :: p_DdataTransport, p_DdataHydro
    real(DP), dimension(:), pointer :: p_DmassMatrix, p_DdataForce
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    character(LEN=SYS_STRLEN) :: slorentzforceName
    real(DP) :: dcurrentDrive, deffectiveRadius
    integer :: isystemFormat, massMatrix
    integer :: neq, nvar, icomp, icoordsystem, ilorentzforcetype
    logical :: bcompatible

    
    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist, ssectionName,&
        'slorentzforcename', slorentzforceName)
    call parlst_getvalue_double(rparlist, ssectionName,&
        'deffectiveradius', deffectiveRadius)
    call parlst_getvalue_int(rparlist, ssectionName,&
        'icoordsystem', icoordsystem)
    call parlst_getvalue_int(rparlist, ssectionName,&
        'ilorentzforcetype', ilorentzforcetype)
    call parlst_getvalue_int(rparlist, ssectionNameHydro,&
        'isystemformat', isystemFormat)
    
    ! Check if solution vector and Lorentz force vector are compatible
    call lsysbl_isVectorCompatible(rsolutionHydro, rforce, bcompatible)
    if (.not.bcompatible)&
        call lsysbl_resizeVectorBlock(rforce, rsolutionHydro, .false.)
        
    ! Get lumped or consistent mass matrix
    select case(ilorentzforcetype)
    case (MASS_LUMPED)
      call parlst_getvalue_int(rparlist,&
          ssectionNameHydro, 'lumpedmassmatrix', massmatrix)
      if (massMatrix .gt. 0) then
        call lsyssc_getbase_double(&
            rproblemLevel%Rmatrix(massMatrix), p_DmassMatrix)
      else
        call output_line('Lumped mass matrix is not available!',&
            OU_CLASS_ERROR, OU_MODE_STD, 'zpinch_calcLorentzforceTerm')
        call sys_halt()
      end if
      
    case (MASS_CONSISTENT)
      call parlst_getvalue_int(rparlist,&
          ssectionNameHydro, 'consistentmassmatrix', massmatrix)
      if (massMatrix .gt. 0) then
        call lsyssc_getbase_double(&
            rproblemLevel%Rmatrix(massMatrix), p_DmassMatrix)
        call lsyssc_getbase_Kld(rproblemLevel%Rmatrix(massmatrix), p_Kld)
        call lsyssc_getbase_Kcol(rproblemLevel%Rmatrix(massmatrix), p_Kcol)
      else
        call output_line('Consistent mass matrix is not available!',&
            OU_CLASS_ERROR, OU_MODE_STD, 'zpinch_calcLorentzforceTerm')
        call sys_halt()
      end if
    end select
    
    ! Set pointer to global solution vectors
    call lsysbl_getbase_double(rforce, p_DdataForce)
    call lsysbl_getbase_double(rsolutionHydro, p_DdataHydro)
    call lsysbl_getbase_double(rsolutionTransport, p_DdataTransport)

    ! Set pointer to the vertex coordinates
    call storage_getbase_double2D(&
        rproblemLevel%rtriangulation%h_DvertexCoords, p_DvertexCoords)

    ! Set dimensions
    neq  = rsolutionTransport%NEQ
    nvar = hydro_getNVAR(rproblemLevel)

    ! Get function parser from collection structure
    p_rfparser => collct_getvalue_pars(rcollection,&
        'rfparser', ssectionName=ssectionName)

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
      select case(icoordsystem)
        
      case(COORDS_CARTESIAN)
        ! Lorentz force term in Cartesian coordinate system
        
        select case(ilorentzforcetype)
        case (MASS_LUMPED)
          call calcLorentzForceXYIntlLumped(dcurrentDrive, deffectiveRadius,&
              neq, nvar, bclear, p_DvertexCoords, p_DmassMatrix,&
              p_DdataTransport, p_DdataHydro, p_DdataForce)
          
        case (MASS_CONSISTENT)
          call calcLorentzForceXYIntlConsistent(dcurrentDrive, deffectiveRadius,&
              neq, nvar, bclear, p_DvertexCoords, p_Kld, p_Kcol, p_DmassMatrix,&
              p_DdataTransport, p_DdataHydro, p_DdataForce)
          
        case default
          call output_line('Unsupported Lorentz force type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcLorentzforceTerm')
          call sys_halt()
        end select
        
        
      case(COORDS_AXIALSYMMETRY)
        ! Lorentz force term in axi-symmetric coordinate system
        
        select case(ilorentzforcetype)
        case (MASS_LUMPED)
          call calcLorentzForceRZIntlLumped(dcurrentDrive, deffectiveRadius,&
              neq, nvar, bclear, p_DvertexCoords, p_DmassMatrix,&
              p_DdataTransport, p_DdataHydro, p_DdataForce)
          
        case (MASS_CONSISTENT)
          call calcLorentzForceRZIntlConsistent(dcurrentDrive, deffectiveRadius,&
              neq, nvar, bclear, p_DvertexCoords, p_Kld, p_Kcol, p_DmassMatrix,&
              p_DdataTransport, p_DdataHydro, p_DdataForce)
          
        case default
          call output_line('Unsupported Lorentz force type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcLorentzforceTerm')
          call sys_halt()
        end select
        
      case default
        call output_line('Invalid type of coordinate system!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcLorentzforceTerm')
        call sys_halt()
      end select
      
      
    case (SYSTEM_BLOCKFORMAT)
      ! What type of coordinate system are we?
      select case(icoordsystem)
        
      case(COORDS_CARTESIAN)
        ! Lorentz force term in Cartesian coordinate system
        
        select case(ilorentzforcetype)
        case (MASS_LUMPED)
          call calcLorentzForceXYBlockLumped(dcurrentDrive, deffectiveRadius,&
              neq, nvar, bclear, p_DvertexCoords, p_DmassMatrix,&
              p_DdataTransport, p_DdataHydro, p_DdataForce)
          
        case (MASS_CONSISTENT)
          call calcLorentzForceXYBlockConsistent(dcurrentDrive, deffectiveRadius,&
              neq, nvar, bclear, p_DvertexCoords, p_Kld, p_Kcol, p_DmassMatrix,&
              p_DdataTransport, p_DdataHydro, p_DdataForce)
          
        case default
          call output_line('Unsupported Lorentz force type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcLorentzforceTerm')
          call sys_halt()
        end select
        
      case(COORDS_AXIALSYMMETRY)
        ! Lorentz force term in axi-symmetric coordinate system
        
        select case(ilorentzforcetype)
        case (MASS_LUMPED)
          call calcLorentzForceRZBlockLumped(dcurrentDrive, deffectiveRadius,&
              neq, nvar, bclear, p_DvertexCoords, p_DmassMatrix,&
              p_DdataTransport, p_DdataHydro, p_DdataForce)
          
        case (MASS_CONSISTENT)
          call calcLorentzForceRZBlockConsistent(dcurrentDrive, deffectiveRadius,&
              neq, nvar, bclear, p_DvertexCoords, p_Kld, p_Kcol, p_DmassMatrix,&
              p_DdataTransport, p_DdataHydro, p_DdataForce)
          
        case default
          call output_line('Unsupported Lorentz force type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcLorentzforceTerm')
          call sys_halt()
        end select
        
      case default
        call output_line('Invalid type of coordinate system!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcLorentzforceTerm')
        call sys_halt()
      end select
      
    case default
      call output_line('Invalid system format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'zpinch_calcLorentzforceTerm')
      call sys_halt()
    end select
    
    
  contains

    ! Here, the real working routines follow

    !**************************************************************
    ! Calculate the Lorentz force term in x-y coordinates.
    ! The system is stored in interleave format.

    subroutine calcLorentzForceXYIntlLumped(dcurrentDrive,&
        deffectiveRadius, neq, nvar, bclear, DvertexCoords,&
        DmassMatrix, DdataTransport, DdataHydro, DdataForce)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(nvar,neq), intent(in) :: DdataHydro
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      real(DP), intent(in) :: dcurrentDrive, deffectiveRadius
      integer, intent(in) :: neq, nvar
      logical, intent(in) :: bclear

      real(DP), dimension(nvar,neq), intent(out) :: DdataForce

      ! local variables
      real(DP) :: drad, daux, x1, x2
      integer :: ieq

      if (bclear) then

        ! Loop over all rows
        !$omp parallel do private(x1,x2,drad,daux)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Get coordinates at node i
          x1 = DvertexCoords(1,ieq)
          x2 = DvertexCoords(2,ieq)
          
          ! Compute unit vector starting at the origin
          drad = sqrt(x1*x1 + x2*x2)
          if (drad .gt. SYS_EPSREAL_DP) then
            x1 = x1/drad
            x2 = x2/drad
          else
            x1 = 0.0_DP
            x2 = 0.0_DP
          end if
          
          ! Compute scaling parameter
          daux = dcurrentDrive * DmassMatrix(ieq)*max(0.0_DP, DdataTransport(ieq)) /&
              max(drad, deffectiveRadius)
          
          ! Compute Lorentz force term
          DdataForce(1,ieq) = 0.0_DP
          DdataForce(2,ieq) = daux * x1
          DdataForce(3,ieq) = daux * x2
          DdataForce(4,ieq) = daux * (X_MOMENTUM_1T_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)*x1+&
                                      Y_MOMENTUM_1T_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)*x2)/&
                                       DENSITY_1T_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)
        end do
        !$omp end parallel do

      else

        ! Loop over all rows
        !$omp parallel do private(x1,x2,drad,daux)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Get coordinates at node i
          x1 = DvertexCoords(1,ieq)
          x2 = DvertexCoords(2,ieq)
          
          ! Compute unit vector starting at the origin
          drad = sqrt(x1*x1 + x2*x2)
          if (drad .gt. SYS_EPSREAL_DP) then
            x1 = x1/drad
            x2 = x2/drad
          else
            x1 = 0.0_DP
            x2 = 0.0_DP
          end if
          
          ! Compute scaling parameter
          daux = dcurrentDrive * DmassMatrix(ieq)*max(0.0_DP, DdataTransport(ieq)) /&
              max(drad, deffectiveRadius)
          
          ! Compute Lorentz force term
          DdataForce(2,ieq) = DdataForce(2,ieq) + daux * x1
          DdataForce(3,ieq) = DdataForce(3,ieq) + daux * x2
          DdataForce(4,ieq) = DdataForce(4,ieq)&
                            + daux * (X_MOMENTUM_1T_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)*x1+&
                                      Y_MOMENTUM_1T_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)*x2)/&
                                         DENSITY_1T_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)
        end do
        !$omp end parallel do

      end if

    end subroutine calcLorentzForceXYIntlLumped

    !**************************************************************
    ! Calculate the Lorentz force term in x-y coordinates.
    ! The system is stored in interleave format.

    subroutine calcLorentzForceXYIntlConsistent(dcurrentDrive,&
        deffectiveRadius, neq, nvar, bclear, DvertexCoords,&
        Kld, Kcol, DmassMatrix, DdataTransport, DdataHydro, DdataForce)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(nvar,neq), intent(in) :: DdataHydro
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      integer, dimension(:), intent(in) :: Kld, Kcol
      real(DP), intent(in) :: dcurrentDrive, deffectiveRadius
      integer, intent(in) :: neq, nvar
      logical, intent(in) :: bclear

      real(DP), dimension(nvar,neq), intent(out) :: DdataForce

      ! local variables
      real(DP), dimension(NVAR2D) :: Ddata
      real(DP) :: drad, daux, x1, x2
      integer :: ieq,ia,jeq

      if (bclear) then

        ! Loop over all rows
        !$omp parallel do private(Ddata,daux,drad,ia,jeq,x1,x2)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)
            
            ! Get coordinates at node jeq
            x1 = DvertexCoords(1,jeq)
            x2 = DvertexCoords(2,jeq)
            
            ! Compute unit vector starting at the origin
            drad = sqrt(x1*x1 + x2*x2)
            if (drad .gt. SYS_EPSREAL_DP) then
              x1 = x1/drad
              x2 = x2/drad
            else
              x1 = 0.0_DP
              x2 = 0.0_DP
            end if
            
            ! Compute scaling parameter
            daux = dcurrentDrive * DmassMatrix(jeq)*max(0.0_DP, DdataTransport(jeq)) /&
                max(drad, deffectiveRadius)
            
            ! Update Lorentz force term
            Ddata(2) = Ddata(2) + daux * x1
            Ddata(3) = Ddata(3) + daux * x2
            Ddata(4) = Ddata(4) + daux * (X_MOMENTUM_1T_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)*x1+&
                                          Y_MOMENTUM_1T_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)*x2)/&
                                             DENSITY_1T_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)
          end do

          ! Overwrite the Lorentz force term
          DdataForce(:,ieq) = Ddata
        end do
        !$omp end parallel do

      else

        ! Loop over all rows
        !$omp parallel do private(Ddata,daux,drad,ia,jeq,x1,x2)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq

          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)
          
            ! Get coordinates at node i
            x1 = DvertexCoords(1,jeq)
            x2 = DvertexCoords(2,jeq)
            
            ! Compute unit vector starting at the origin
            drad = sqrt(x1*x1 + x2*x2)
            if (drad .gt. SYS_EPSREAL_DP) then
              x1 = x1/drad
              x2 = x2/drad
            else
              x1 = 0.0_DP
              x2 = 0.0_DP
            end if
            
            ! Compute scaling parameter
            daux = dcurrentDrive * DmassMatrix(jeq)*max(0.0_DP, DdataTransport(jeq)) /&
                max(drad, deffectiveRadius)
            
            ! Update Lorentz force term
            Ddata(2) = Ddata(2) + daux * x1
            Ddata(3) = Ddata(3) + daux * x2
            Ddata(4) = Ddata(4)&
                     + daux * (X_MOMENTUM_1T_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)*x1+&
                               Y_MOMENTUM_1T_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)*x2)/&
                                  DENSITY_1T_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)
          end do

          ! Update  Lorentz force term
          DdataForce(:,ieq) = DdataForce(:,ieq) + Ddata
        end do
        !$omp end parallel do

      end if

    end subroutine calcLorentzForceXYIntlConsistent
    

    !**************************************************************
    ! Calculate the Lorentz force term in r-z coordinates.
    ! The system is stored in interleave format.

    subroutine calcLorentzForceRZIntlLumped(dcurrentDrive,&
        deffectiveRadius, neq, nvar, bclear, DvertexCoords, DmassMatrix,&
        DdataTransport, DdataHydro, DdataForce)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(nvar,neq), intent(in) :: DdataHydro
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      real(DP), intent(in) :: dcurrentDrive, deffectiveRadius
      integer, intent(in) :: neq, nvar
      logical, intent(in) :: bclear
      
      real(DP), dimension(nvar,neq), intent(out) :: DdataForce
      
      ! local variables
      real(DP) :: drad, daux, x1
      integer :: ieq

      if (bclear) then

        ! Loop over all rows
        !$omp parallel do private(x1,drad,daux)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Get x-coordinate at node i
          x1 = DvertexCoords(1,ieq); drad = abs(x1)
          
          ! Compute unit vector starting at the origin
          if (drad .gt. SYS_EPSREAL_DP) then
            x1 = sign(1.0_DP, x1)
          else
            x1 = 0.0_DP
          end if
          
          ! Compute scaling parameter
          daux = dcurrentDrive * DmassMatrix(ieq)*max(0.0_DP, DdataTransport(ieq)) /&
              max(drad, deffectiveRadius)
          
          ! Compute Lorentz source term        
          DdataForce(1,ieq) = 0.0_DP
          DdataForce(2,ieq) = daux * x1
          DdataForce(3,ieq) = 0.0_DP
          DdataForce(4,ieq) = daux * X_MOMENTUM_1T_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)*x1/&
                                        DENSITY_1T_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)
        end do
        !$omp end parallel do

      else

        ! Loop over all rows
        !$omp parallel do private(x1,drad,daux)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Get x-coordinate at node i
          x1 = DvertexCoords(1,ieq); drad = abs(x1)
          
          ! Compute unit vector starting at the origin
          if (drad .gt. SYS_EPSREAL_DP) then
            x1 = sign(1.0_DP, x1)
          else
            x1 = 0.0_DP
          end if
          
          ! Compute scaling parameter
          daux = dcurrentDrive * DmassMatrix(ieq)*max(0.0_DP, DdataTransport(ieq)) /&
              max(drad, deffectiveRadius)
          
          ! Compute Lorentz source term
          DdataForce(2,ieq) = DdataForce(2,ieq) + daux * x1          
          DdataForce(4,ieq) = DdataForce(4,ieq)&
                            + daux * X_MOMENTUM_1T_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)*x1/&
                                        DENSITY_1T_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)
        end do
        !$omp end parallel do

      end if

    end subroutine calcLorentzForceRZIntlLumped

    
    !**************************************************************
    ! Calculate the Lorentz force term in r-z coordinates.
    ! The system is stored in interleave format.

    subroutine calcLorentzForceRZIntlConsistent(dcurrentDrive,&
        deffectiveRadius, neq, nvar, bclear, DvertexCoords, Kld,&
        Kcol, DmassMatrix, DdataTransport, DdataHydro, DdataForce)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(nvar,neq), intent(in) :: DdataHydro
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      integer, dimension(:), pointer :: Kld, Kcol
      real(DP), intent(in) :: dcurrentDrive, deffectiveRadius
      integer, intent(in) :: neq, nvar
      logical, intent(in) :: bclear
      
      real(DP), dimension(nvar,neq), intent(out) :: DdataForce
      
      ! local variables
      real(DP), dimension(NVAR2D) :: Ddata
      real(DP) :: drad, daux, x1
      integer :: ieq,ia,jeq

      if (bclear) then

        ! Loop over all rows
        !$omp parallel do private(Ddata,daux,drad,ia,jeq,x1)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)

            ! Get x-coordinate at node i
            x1 = DvertexCoords(1,jeq); drad = abs(x1)
            
            ! Compute unit vector starting at the origin
            if (drad .gt. SYS_EPSREAL_DP) then
              x1 = sign(1.0_DP, x1)
            else
              x1 = 0.0_DP
            end if
            
            ! Compute scaling parameter
            daux = dcurrentDrive * DmassMatrix(jeq)*max(0.0_DP, DdataTransport(jeq)) /&
                max(drad, deffectiveRadius)
            
            ! Update Lorentz source term        
            Ddata(2) = Ddata(2) + daux * x1
            Ddata(4) = Ddata(4) + daux * X_MOMENTUM_1T_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)*x1/&
                                            DENSITY_1T_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)
          end do
          
          ! Overwrite the Lorentz force term
          DdataForce(:,ieq) = Ddata
        end do
        !$omp end parallel do

      else

        ! Loop over all rows
        !$omp parallel do private(Ddata,daux,drad,ia,jeq,x1)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Clear temporal data
          Ddata = 0.0_DP
          
          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1
            
            ! Get nodal degree of freedom
            jeq = Kcol(ia)
            
            ! Get x-coordinate at node i
            x1 = DvertexCoords(1,jeq); drad = abs(x1)
            
            ! Compute unit vector starting at the origin
            if (drad .gt. SYS_EPSREAL_DP) then
              x1 = sign(1.0_DP, x1)
            else
              x1 = 0.0_DP
            end if
            
            ! Compute scaling parameter
            daux = dcurrentDrive * DmassMatrix(jeq)*max(0.0_DP, DdataTransport(jeq)) /&
                max(drad, deffectiveRadius)
          
            ! Update Lorentz source term
            Ddata(2) = Ddata(2) + daux * x1          
            Ddata(4) = Ddata(4)&
                     + daux * X_MOMENTUM_1T_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)*x1/&
                                 DENSITY_1T_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)
          end do

          ! Update Lorentz source term
          DdataForce(:,ieq) = DdataForce(:,ieq) + Ddata
        end do
        !$omp end parallel do

      end if

    end subroutine calcLorentzForceRZIntlConsistent


    !**************************************************************
    ! Calculate the Lorentz force term in x-y coordinates.
    ! The system is stored in block format.

    subroutine calcLorentzForceXYBlockLumped(dcurrentDrive,&
        deffectiveRadius, neq, nvar, bclear, DvertexCoords,&
        DmassMatrix, DdataTransport, DdataHydro, DdataForce)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(neq,nvar), intent(in) :: DdataHydro
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      real(DP), intent(in) :: dcurrentDrive, deffectiveRadius
      integer, intent(in) :: neq, nvar
      logical, intent(in) :: bclear

      real(DP), dimension(neq,nvar), intent(out) :: DdataForce

      ! local variables
      real(DP) :: drad, daux, x1, x2
      integer :: ieq

      if (bclear) then

        ! Loop over all rows
        !$omp parallel do private(x1,x2,drad,daux)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Get coordinates at node i
          x1 = DvertexCoords(1,ieq)
          x2 = DvertexCoords(2,ieq)
          
          ! Compute unit vector starting at the origin
          drad = sqrt(x1*x1 + x2*x2)
          if (drad .gt. SYS_EPSREAL_DP) then
            x1 = x1/drad
            x2 = x2/drad
          else
            x1 = 0.0_DP
            x2 = 0.0_DP
          end if
          
          ! Compute scaling parameter
          daux = dcurrentDrive * DmassMatrix(ieq)*max(0.0_DP, DdataTransport(ieq)) /&
              max(drad, deffectiveRadius)
          
          ! Compute Lorentz source term
          DdataForce(ieq,1) = 0.0_DP
          DdataForce(ieq,2) = daux * x1
          DdataForce(ieq,3) = daux * x2
          DdataForce(ieq,4) = daux * (X_MOMENTUM_1L_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)*x1+&
                                      Y_MOMENTUM_1L_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)*x2)/&
                                         DENSITY_1L_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)
        end do
        !$omp end parallel do

      else

        ! Loop over all rows
        !$omp parallel do private(x1,x2,drad,daux)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Get coordinates at node i
          x1 = DvertexCoords(1,ieq)
          x2 = DvertexCoords(2,ieq)
          
          ! Compute unit vector starting at the origin
          drad = sqrt(x1*x1 + x2*x2)
          if (drad .gt. SYS_EPSREAL_DP) then
            x1 = x1/drad
            x2 = x2/drad
          else
            x1 = 0.0_DP
            x2 = 0.0_DP
          end if
          
          ! Compute scaling parameter
          daux = dcurrentDrive * DmassMatrix(ieq)*max(0.0_DP, DdataTransport(ieq)) /&
              max(drad, deffectiveRadius)
          
          ! Compute Lorentz source term
          DdataForce(ieq,2) = DdataForce(ieq,2) + daux * x1
          DdataForce(ieq,3) = DdataForce(ieq,3) + daux * x2
          DdataForce(ieq,4) = DdataForce(ieq,4)&
                            + daux * (X_MOMENTUM_1L_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)*x1+&
                                      Y_MOMENTUM_1L_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)*x2)/&
                                         DENSITY_1L_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)
        end do
        !$omp end parallel do

      end if

    end subroutine calcLorentzForceXYBlockLumped


    !**************************************************************
    ! Calculate the Lorentz force term in x-y coordinates.
    ! The system is stored in block format.

    subroutine calcLorentzForceXYBlockConsistent(dcurrentDrive,&
        deffectiveRadius, neq, nvar, bclear, DvertexCoords, Kld,&
        Kcol, DmassMatrix, DdataTransport, DdataHydro, DdataForce)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(neq,nvar), intent(in) :: DdataHydro
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      integer, dimension(:), intent(in) :: Kld, Kcol
      real(DP), intent(in) :: dcurrentDrive, deffectiveRadius
      integer, intent(in) :: neq, nvar
      logical, intent(in) :: bclear

      real(DP), dimension(neq,nvar), intent(out) :: DdataForce

      ! local variables
      real(DP), dimension(NVAR1D) :: Ddata
      real(DP) :: drad, daux, x1, x2
      integer :: ieq,ia,jeq

      if (bclear) then

        ! Loop over all rows
        !$omp parallel do private(Ddata,daux,drad,ia,jeq,x1,x2)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)

            ! Get coordinates at node jeq
            x1 = DvertexCoords(1,jeq)
            x2 = DvertexCoords(2,jeq)
            
            ! Compute unit vector starting at the origin
            drad = sqrt(x1*x1 + x2*x2)
            if (drad .gt. SYS_EPSREAL_DP) then
              x1 = x1/drad
              x2 = x2/drad
            else
              x1 = 0.0_DP
              x2 = 0.0_DP
            end if
            
            ! Compute scaling parameter
            daux = dcurrentDrive * DmassMatrix(jeq)*max(0.0_DP, DdataTransport(jeq)) /&
                max(drad, deffectiveRadius)
            
            ! Update Lorentz source term
            Ddata(2) = Ddata(2) + daux * x1
            Ddata(3) = Ddata(3) + daux * x2
            Ddata(4) = Ddata(4) + daux * (X_MOMENTUM_1L_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)*x1+&
                                          Y_MOMENTUM_1L_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)*x2)/&
                                             DENSITY_1L_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)
          end do

          ! Overwrite Lorentz source term
          DdataForce(ieq,:) = Ddata
        end do
        !$omp end parallel do

      else

        ! Loop over all rows
        !$omp parallel do private(Ddata,daux,drad,ia,jeq,x1,x2)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)
          
            ! Get coordinates at node jeq
            x1 = DvertexCoords(1,jeq)
            x2 = DvertexCoords(2,jeq)
            
            ! Compute unit vector starting at the origin
            drad = sqrt(x1*x1 + x2*x2)
            if (drad .gt. SYS_EPSREAL_DP) then
              x1 = x1/drad
              x2 = x2/drad
            else
              x1 = 0.0_DP
              x2 = 0.0_DP
            end if
            
            ! Compute scaling parameter
            daux = dcurrentDrive * DmassMatrix(jeq)*max(0.0_DP, DdataTransport(jeq)) /&
                max(drad, deffectiveRadius)
            
            ! Update Lorentz source term
            Ddata(2) = Ddata(2) + daux * x1
            Ddata(3) = Ddata(3) + daux * x2
            Ddata(4) = Ddata(4)&
                     + daux * (X_MOMENTUM_1L_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)*x1+&
                               Y_MOMENTUM_1L_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)*x2)/&
                                  DENSITY_1L_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)
          end do

          ! Update Lorentz source term
          DdataForce(ieq,:) = DdataForce(ieq,:) + Ddata
        end do
        !$omp end parallel do

      end if

    end subroutine calcLorentzForceXYBlockConsistent


    !**************************************************************
    ! Calculate the Lorentz force term in r-z coordinates.
    ! The system is stored in block format.

    subroutine calcLorentzForceRZBlockLumped(dcurrentDrive,&
        deffectiveRadius, neq, nvar, bclear, DvertexCoords, DmassMatrix,&
        DdataTransport, DdataHydro, DdataForce)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(neq,nvar), intent(in) :: DdataHydro
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      real(DP), intent(in) :: dcurrentDrive, deffectiveRadius
      integer, intent(in) :: neq, nvar
      logical, intent(in) :: bclear
      
      real(DP), dimension(neq,nvar), intent(out) :: DdataForce
      
      ! local variables
      real(DP) :: drad, daux, x1
      integer :: ieq

      if (bclear) then

        ! Loop over all rows
        !$omp parallel do private(x1,drad,daux)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Get x-coordinate at node i
          x1 = DvertexCoords(1,ieq); drad = abs(x1)
          
          ! Compute unit vector starting at the origin
          if (drad .gt. SYS_EPSREAL_DP) then
            x1 = sign(1.0_DP,x1)
          else
            x1 = 0.0_DP
          end if
          
          ! Compute scaling parameter
          daux = dcurrentDrive * DmassMatrix(ieq)*max(0.0_DP, DdataTransport(ieq)) /&
              max(drad, deffectiveRadius)
          
          ! Compute Lorentz source term
          DdataForce(ieq,1) = 0.0_DP
          DdataForce(ieq,2) = daux * x1
          DdataForce(ieq,3) = 0.0_DP
          DdataForce(ieq,4) = daux * X_MOMENTUM_1L_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)*x1/&
                                        DENSITY_1L_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)
        end do
        !$omp end parallel do

      else

        ! Loop over all rows
        !$omp parallel do private(x1,drad,daux)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Get x-coordinate at node i
          x1 = DvertexCoords(1,ieq); drad = abs(x1)
          
          ! Compute unit vector starting at the origin
          if (drad .gt. SYS_EPSREAL_DP) then
            x1 = sign(1.0_DP,x1)
          else
            x1 = 0.0_DP
          end if
          
          ! Compute scaling parameter
          daux = dcurrentDrive * DmassMatrix(ieq)*max(0.0_DP, DdataTransport(ieq)) /&
              max(drad, deffectiveRadius)
          
          ! Compute Lorentz source term
          DdataForce(ieq,2) = DdataForce(ieq,2)  + daux * x1
          DdataForce(ieq,4) = DdataForce(ieq,4)&
                            + daux * X_MOMENTUM_1L_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)*x1/&
                                        DENSITY_1L_FROM_CONSVAR(DdataHydro,NVAR2D,ieq)
        end do
        !$omp end parallel do
        
      end if

    end subroutine calcLorentzForceRZBlockLumped


    !**************************************************************
    ! Calculate the Lorentz force term in r-z coordinates.
    ! The system is stored in block format.

    subroutine calcLorentzForceRZBlockConsistent(dcurrentDrive,&
        deffectiveRadius, neq, nvar, bclear, DvertexCoords, Kld,&
        Kcol, DmassMatrix, DdataTransport, DdataHydro, DdataForce)

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      real(DP), dimension(neq,nvar), intent(in) :: DdataHydro
      real(DP), dimension(:), intent(in) :: DmassMatrix, DdataTransport
      integer, dimension(:), intent(in) :: Kld, Kcol
      real(DP), intent(in) :: dcurrentDrive, deffectiveRadius
      integer, intent(in) :: neq, nvar
      logical, intent(in) :: bclear
      
      real(DP), dimension(neq,nvar), intent(out) :: DdataForce
      
      ! local variables
      real(DP), dimension(NVAR2D) :: Ddata
      real(DP) :: drad, daux, x1
      integer :: ieq,ia,jeq

      if (bclear) then

        ! Loop over all rows
        !$omp parallel do private(Ddata,daux,drad,ia,jeq,x1)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)
            
            ! Get x-coordinate at node jeq
            x1 = DvertexCoords(1,jeq); drad = abs(x1)
            
            ! Compute unit vector starting at the origin
            if (drad .gt. SYS_EPSREAL_DP) then
              x1 = sign(1.0_DP,x1)
            else
              x1 = 0.0_DP
            end if
            
            ! Compute scaling parameter
            daux = dcurrentDrive * DmassMatrix(jeq)*max(0.0_DP, DdataTransport(jeq)) /&
                max(drad, deffectiveRadius)
            
            ! Update Lorentz source term
            Ddata(2) = Ddata(2) + daux * x1
            Ddata(4) = Ddata(4) + daux * X_MOMENTUM_1L_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)*x1/&
                                            DENSITY_1L_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)
          end do

          ! Overwrite Lorentz source term
          DdataForce(ieq,:) = Ddata
        end do
        !$omp end parallel do

      else

        ! Loop over all rows
        !$omp parallel do private(Ddata,daux,drad,ia,jeq,x1)&
        !$omp if (neq > ZPINCH_LORENTZFORCE_NEQMIN_OMP)
        do ieq = 1, neq
          
          ! Clear temporal data
          Ddata = 0.0_DP

          ! Loop over all coupled degrees of freedom
          do ia = Kld(ieq), Kld(ieq+1)-1

            ! Get nodal degree of freedom
            jeq = Kcol(ia)

            ! Get x-coordinate at node jeq
            x1 = DvertexCoords(1,jeq); drad = abs(x1)
            
            ! Compute unit vector starting at the origin
            if (drad .gt. SYS_EPSREAL_DP) then
              x1 = sign(1.0_DP,x1)
            else
              x1 = 0.0_DP
            end if
            
            ! Compute scaling parameter
            daux = dcurrentDrive * DmassMatrix(jeq)*max(0.0_DP, DdataTransport(jeq)) /&
                max(drad, deffectiveRadius)
            
            ! Update Lorentz source term
            Ddata(2) = Ddata(2)  + daux * x1
            Ddata(4) = Ddata(4)&
                     + daux * X_MOMENTUM_1L_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)*x1/&
                                 DENSITY_1L_FROM_CONSVAR(DdataHydro,NVAR2D,jeq)
          end do

          ! Update Lorentz source term
          DdataForce(ieq,:) = DdataForce(ieq,:) + Ddata
        end do
        !$omp end parallel do
        
      end if

    end subroutine calcLorentzForceRZBlockConsistent
    
  end subroutine zpinch_calcLorentzforceTerm

  !*****************************************************************************

!<subroutine>
  
  subroutine zpinch_calcLinearisedFCT(RbdrCond, rproblemLevel, rtimestep,&
      rsolverHydro, rsolverTransport, Rsolution, ssectionName,&
      ssectionNameHydro, ssectionNameTransport, rcollection, Rsource)

!<description>
    ! This subroutine calculates the linearised FCT correction
!</description>

!<input>
    ! array of boundary condition structure
    type(t_boundaryCondition), dimension(:), intent(in) :: RbdrCond

    ! time-stepping algorithm
    type(t_timestep), intent(in) :: rtimestep

    ! section names in parameter list
    character(LEN=*), intent(in) :: ssectionName
    character(LEN=*), intent(in) :: ssectionNameHydro
    character(LEN=*), intent(in) :: ssectionNameTransport

    ! OPTIONAL: source vector
    type(t_vectorBlock), dimension(:), intent(in), optional :: Rsource
!</input>

!<inputoutput>
    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structures
    type(t_solver), intent(inout) :: rsolverHydro
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
    type(t_vectorBlock), pointer :: p_rpredictorHydro,p_rpredictorTransport
    type(t_parlist), pointer :: p_rparlist
    character(len=SYS_STRLEN), dimension(:), pointer :: SfailsafeVariables
    real(DP), dimension(:), pointer :: p_DalphaHydro, p_DalphaTransport
    integer, dimension(2) :: IposAFC
    integer :: convectionAFC,inviscidAFC
    integer :: lumpedMassMatrixHydro
    integer :: lumpedMassMatrixTransport,consistentMassMatrixTransport
    integer :: imassantidiffusiontypeHydro,imassantidiffusiontypeTransport
    integer :: ilimitersynchronisation,nfailsafe,ivariable,nvariable


    ! Set pointer to parameter list
    p_rparlist => collct_getvalue_parlst(rcollection,&
        'rparlist', ssectionName=ssectionName)

    ! Get parameters from parameter list
    call parlst_getvalue_int(p_rparlist,&
        ssectionNameHydro, 'inviscidAFC', inviscidAFC)
    call parlst_getvalue_int(p_rparlist,&
        ssectionNameTransport, 'convectionAFC', convectionAFC)

    ! Do we have to apply linearised FEM-FCT?
    if ((inviscidAFC .le. 0) .or. (convectionAFC .le. 0)) return
    if ((rproblemLevel%Rafcstab(inviscidAFC)%ctypeAFCstabilisation&
         .ne. AFCSTAB_FEMFCT_LINEARISED) .or.&
        (rproblemLevel%Rafcstab(convectionAFC)%ctypeAFCstabilisation&
         .ne. AFCSTAB_FEMFCT_LINEARISED)) return

    ! Set positions of stabilisation structures
    IposAFC=(/inviscidAFC, convectionAFC/)

    ! Get more parameters from parameter list.
    call parlst_getvalue_int(p_rparlist, ssectionName,&
        'ilimitersynchronisation', ilimitersynchronisation)
    
    call parlst_getvalue_int(p_rparlist, ssectionNameHydro,&
        'lumpedmassmatrix', lumpedmassmatrixHydro)
    call parlst_getvalue_int(p_rparlist, ssectionNameHydro,&
        'nfailsafe', nfailsafe)
    call parlst_getvalue_int(p_rparlist, ssectionNameHydro,&
        'imassantidiffusiontype', imassantidiffusiontypeHydro)

    call parlst_getvalue_int(p_rparlist, ssectionNameTransport,&
        'lumpedmassmatrix', lumpedmassmatrixTransport)
    call parlst_getvalue_int(p_rparlist, ssectionNameTransport,&
        'consistentmassmatrix', consistentmassmatrixTransport)
    call parlst_getvalue_int(p_rparlist, ssectionNameTransport,&
        'imassantidiffusiontype', imassantidiffusiontypeTransport)
    
    !---------------------------------------------------------------------------
    ! Linearised FEM-FCT algorithm
    !---------------------------------------------------------------------------

    !--- compressible hydrodynamic model ---------------------------------------
    
    ! Should we apply consistent mass antidiffusion?
    if (imassantidiffusiontypeHydro .eq. MASS_CONSISTENT) then

      ! Initialize dummy timestep
      rtimestepAux%dStep = 1.0_DP
      rtimestepAux%theta = 0.0_DP
          
      ! Set pointer to predictor
      p_rpredictorHydro => rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorPredictor

      ! Compute low-order "right-hand side" without theta parameter
      if (present(Rsource)) then
        if (Rsource(1)%NEQ .eq. 0) then
          ! ... without source term
          call hydro_calcRhsThetaScheme(rproblemLevel, rtimestepAux,&
              rsolverHydro, Rsolution(1), p_rpredictorHydro,&
              ssectionNameHydro,rcollection)
        else
          ! ... with source term
          call hydro_calcRhsThetaScheme(rproblemLevel, rtimestepAux,&
              rsolverHydro, Rsolution(1), p_rpredictorHydro,&
              ssectionNameHydro, rcollection, Rsource(1))
        end if
      else
        ! ... without source term
        call hydro_calcRhsThetaScheme(rproblemLevel, rtimestepAux,&
            rsolverHydro, Rsolution(1), p_rpredictorHydro,&
            ssectionNameHydro, rcollection)
      end if

      ! Compute low-order predictor
      call lsysbl_invertedDiagMatVec(&
          rproblemLevel%Rmatrix(lumpedMassMatrixHydro),&
          p_rpredictorHydro, 1.0_DP, p_rpredictorHydro)

      ! Build the raw antidiffusive fluxes with contribution from
      ! consistent mass matrix
      call hydro_calcFluxFCT(rproblemLevel, Rsolution(1), 0.0_DP,&
          1.0_DP, 1.0_DP, .true., p_rpredictorHydro,&
          ssectionNameHydro, rcollection)

    else
      ! Build the raw antidiffusive fluxes without contribution from
      ! consistent mass matrix
      call hydro_calcFluxFCT(rproblemLevel, Rsolution(1), 0.0_DP,&
          1.0_DP, 1.0_DP, .true., Rsolution(1),&
          ssectionNameHydro, rcollection)
    end if

    !--- transport model -------------------------------------------------------

    ! Should we apply consistent mass antidiffusion?
    if (imassantidiffusiontypeTransport .eq. MASS_CONSISTENT) then
      
      ! Initialize dummy timestep
      rtimestepAux%dStep = 1.0_DP
      rtimestepAux%theta = 0.0_DP
      
      ! Set pointer to predictor
      p_rpredictorTransport => rproblemLevel%Rafcstab(convectionAFC)%p_rvectorPredictor
      
      ! Compute the preconditioner
      call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
          rsolverTransport, Rsolution(2), ssectionNameTransport, rcollection)
      
      ! Compute low-order "right-hand side" without theta parameter
      if (present(Rsource)) then
        if (Rsource(2)%NEQ .eq. 0) then
          ! ... without source term
          call transp_calcRhsThetaScheme(rproblemLevel, rtimestepAux,&
              rsolverTransport, Rsolution(2), p_rpredictorTransport,&
              ssectionNameTransport, rcollection,&
              fcb_coeffVecBdrPrimal2d_sim = transp_coeffVecBdrConvP2d_sim,&
              fcb_coeffVecBdrDual2d_sim = transp_coeffVecBdrConvD2d_sim)
        else
          ! ... with source term
          call transp_calcRhsThetaScheme(rproblemLevel, rtimestepAux,&
              rsolverTransport, Rsolution(2), p_rpredictorTransport,&
              ssectionNameTransport, rcollection, Rsource(2),&
              fcb_coeffVecBdrPrimal2d_sim = transp_coeffVecBdrConvP2d_sim,&
              fcb_coeffVecBdrDual2d_sim = transp_coeffVecBdrConvD2d_sim)
        end if
      else
        ! ... without source term
        call transp_calcRhsThetaScheme(rproblemLevel, rtimestepAux,&
            rsolverTransport, Rsolution(2), p_rpredictorTransport,&
            ssectionNameTransport, rcollection,&
            fcb_coeffVecBdrPrimal2d_sim = transp_coeffVecBdrConvP2d_sim,&
            fcb_coeffVecBdrDual2d_sim = transp_coeffVecBdrConvD2d_sim)
      end if
      
      ! Compute low-order predictor
      call lsysbl_invertedDiagMatVec(&
          rproblemLevel%Rmatrix(lumpedMassMatrixTransport),&
          p_rpredictorTransport, 1.0_DP, p_rpredictorTransport)
      
      ! Build the raw antidiffusive fluxes with contribution from
      ! consistent mass matrix
      call gfsc_buildFluxFCT(rproblemLevel%Rafcstab(convectionAFC),&
          Rsolution(2), 0.0_DP, 1.0_DP, 1.0_DP, .true.,&
          rproblemLevel%Rmatrix(consistentMassMatrixTransport),&
          p_rpredictorTransport)
    else
      ! Build the raw antidiffusive fluxes without contribution from
      ! consistent mass matrix
      call gfsc_buildFluxFCT(rproblemLevel%Rafcstab(convectionAFC),&
          Rsolution(2), 0.0_DP, 1.0_DP, 1.0_DP, .true.)
    end if

    !---------------------------------------------------------------------------
    ! Perform failsafe flux correction (if required)
    !---------------------------------------------------------------------------
    
    if (nfailsafe .gt. 0) then
      
      ! Get number of failsafe variables
      nvariable = max(1, parlst_querysubstrings(p_rparlist,&
          ssectionNameHydro, 'sfailsafevariable'))
      
      ! Allocate character array that stores all failsafe variable names
      allocate(SfailsafeVariables(nvariable))
      
      ! Initialize character array with failsafe variable names
      do ivariable = 1, nvariable
        call parlst_getvalue_string(p_rparlist,&
            ssectionNameHydro, 'sfailsafevariable',&
            Sfailsafevariables(ivariable), isubstring=ivariable)
      end do

      ! Set up the predictor, that is, make a virtual copy of the
      ! predictor vectors from the hydrodynamic system and the scalar tracer
      ! equation which can be passed to the failsafe subroutine.
      allocate(p_rpredictor(2))
      call lsysbl_duplicateVector(p_rpredictorHydro, p_Rpredictor(1),&
          LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_SHARE)
      call lsysbl_duplicateVector(p_rpredictorTransport, p_Rpredictor(2),&
          LSYSSC_DUP_TEMPLATE, LSYSSC_DUP_SHARE)
    end if

    ! What type of limiter synchronisation is performed?
    select case(ilimitersynchronisation)
      
    case (0)   ! no synchronisation
      
      if (nfailsafe .gt. 0) then
        
        ! Compute linearised FEM-FCT correction for the hydrodynamic model
        call hydro_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
            AFCSTAB_FCTALGO_CORRECT, Rsolution(1),&
            ssectionNameHydro, rcollection)
        
        ! Apply failsafe flux correction for the hydrodynamic model
        call afcstab_failsafeLimiting(rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrixHydro),&
            SfailsafeVariables, rtimestep%dStep, nfailsafe,&
            hydro_getVariable, Rsolution(1), p_rpredictorHydro)

        ! Compute and apply linearised FEM-FCT correction for the
        ! transport model
        call gfsc_buildConvectionVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrixTransport),&
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
        
        ! Compute and apply linearised FEM-FCT correction for the
        ! hydrodynamic model
        call hydro_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(1),&
            ssectionNameHydro, rcollection)
        
        ! Compute and apply linearised FEM-FCT correction for the
        ! transport model
        call gfsc_buildConvectionVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrixTransport),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))
      end if
      

    case (1)   ! minimum synchronisation
      
      ! Compute linearised FEM-FCT correction for the hydrodynamic model
      call hydro_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
          rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_CORRECT, Rsolution(1),&
          ssectionNameHydro, rcollection)
      
      ! Compute linearised FEM-FCT correction for the transport model
      call gfsc_buildConvectionVectorFCT(&
          rproblemLevel%Rmatrix(lumpedMassMatrixTransport),&
          rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
          rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_CORRECT, Rsolution(2))
      
      ! Set pointers
      call lsyssc_getbase_double(&
          rproblemLevel%Rafcstab(inviscidAFC)%p_rvectorAlpha, p_DalphaHydro)
      call lsyssc_getbase_double(&
          rproblemLevel%Rafcstab(convectionAFC)%p_rvectorAlpha, p_DalphaTransport)

      ! Compute minimum correction factor from both models
      p_DalphaHydro = min(p_DalphaHydro, p_DalphaTransport)
      call lalg_copyVector(p_DalphaHydro, p_DalphaTransport)
      
      if (nfailsafe .gt. 0) then

        ! Apply failsafe flux correction for both models
        call afcstab_failsafeLimiting(rproblemLevel%Rafcstab(IposAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrixHydro),&
            SfailsafeVariables, rtimestep%dStep, nfailsafe,&
            zpinch_getVariable, Rsolution, p_Rpredictor)
        
        ! Deallocate temporal memory
        call lsysbl_releaseVector(p_Rpredictor(1))
        call lsysbl_releaseVector(p_Rpredictor(2))
        deallocate(SfailsafeVariables, p_Rpredictor)
        
      else
        
        ! Apply linearised FEM-FCT correction for the hydrodynamic model
        call hydro_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(1),&
            ssectionNameHydro, rcollection)
        
        ! Apply linearised FEM-FCT correction for the transportmodel
        call gfsc_buildConvectionVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrixTransport),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))
      end if


    case (2)   ! Hydrodynamic model first, transport model second
      
      ! Compute linearised FEM-FCT correction for the hydrodynamic model
      call hydro_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
          rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_CORRECT, Rsolution(1),&
          ssectionNameHydro, rcollection)
      
      ! Copy correction factor for hydrodynamic model to transport model
      call afcstab_copyStabilisation(&
          rproblemLevel%Rafcstab(inviscidAFC),&
          rproblemLevel%Rafcstab(convectionAFC), AFCSTAB_DUP_EDGELIMITER)
      
      ! Compute linearised FEM-FCT correction for the transport model
      ! (without initialisation)
      call gfsc_buildConvectionVectorFCT(&
          rproblemLevel%Rmatrix(lumpedMassMatrixTransport),&
          rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
          rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_INITALPHA-AFCSTAB_FCTALGO_CORRECT, Rsolution(2))
      
      ! Copy final correction factor back to hydrodynamic system
      call afcstab_copyStabilisation(&
          rproblemLevel%Rafcstab(convectionAFC),&
          rproblemLevel%Rafcstab(inviscidAFC), AFCSTAB_DUP_EDGELIMITER)
      
      if (nfailsafe .gt. 0) then

        ! Apply failsafe flux correction for both models
        call afcstab_failsafeLimiting(rproblemLevel%Rafcstab(IposAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrixHydro),&
            SfailsafeVariables, rtimestep%dStep, nfailsafe,&
            zpinch_getVariable, Rsolution, p_Rpredictor)
        
        ! Deallocate temporal memory
        call lsysbl_releaseVector(p_Rpredictor(1))
        call lsysbl_releaseVector(p_Rpredictor(2))
        deallocate(SfailsafeVariables, p_Rpredictor)

      else
        
        ! Apply linearised FEM-FCT correction for the hydrodynamic model
        call hydro_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(1),&
            ssectionNameHydro, rcollection)
        
        ! Apply linearised FEM-FCT correction for the transport model
        call gfsc_buildConvectionVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrixTransport),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))
      end if


    case (3)   ! Transport model first, hydrodynamic second second
      
      ! Compute linearised FEM-FCT correction for the transport model
      call gfsc_buildConvectionVectorFCT(&
          rproblemLevel%Rmatrix(lumpedMassMatrixTransport),&
          rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
          rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_CORRECT, Rsolution(2))
      
      ! Copy correction factor for transport model to hydrodynamic model
      call afcstab_copyStabilisation(&
          rproblemLevel%Rafcstab(convectionAFC),&
          rproblemLevel%Rafcstab(inviscidAFC), AFCSTAB_DUP_EDGELIMITER)
      
      ! Compute linearised FEM-FCT correction for the hydrodynamic model
      ! (without initialisation)
      call hydro_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
          rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
          AFCSTAB_FCTALGO_INITALPHA-AFCSTAB_FCTALGO_CORRECT,&
          Rsolution(1), ssectionNameHydro, rcollection)
      
      ! Copy final correction factor back to transport model
      call afcstab_copyStabilisation(&
          rproblemLevel%Rafcstab(inviscidAFC),&
          rproblemLevel%Rafcstab(convectionAFC), AFCSTAB_DUP_EDGELIMITER)

      if (nfailsafe .gt. 0) then

        ! Apply failsafe flux correction for both models
        call afcstab_failsafeLimiting(rproblemLevel%Rafcstab(IposAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrixHydro),&
            SfailsafeVariables, rtimestep%dStep, nfailsafe,&
            zpinch_getVariable, Rsolution, p_Rpredictor)
        
        ! Deallocate temporal memory
        call lsysbl_releaseVector(p_Rpredictor(1))
        call lsysbl_releaseVector(p_Rpredictor(2))
        deallocate(SfailsafeVariables, p_Rpredictor)

      else
      
        ! Apply linearised FEM-FCT correction for the hydrodynamic model
        call hydro_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(1),&
            ssectionNameHydro, rcollection)
        
        ! Apply linearised FEM-FCT correction for the transport model
        call gfsc_buildConvectionVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrixTransport),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))
      end if
      

    case (4)   ! Hydrodynamic system only
      
      if (nfailsafe .gt. 0) then

        ! Compute and apply linearised FEM-FCT correction for
        ! hydrodynamic model
        call hydro_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
            AFCSTAB_FCTALGO_CORRECT, Rsolution(1),&
            ssectionNameHydro, rcollection)

        ! Copy correction factor for hydrodynamic model to transport model
        call afcstab_copyStabilisation(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rafcstab(convectionAFC), AFCSTAB_DUP_EDGELIMITER)

        ! Apply failsafe flux correction for both models
        call afcstab_failsafeLimiting(rproblemLevel%Rafcstab(IposAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrixHydro),&
            SfailsafeVariables, rtimestep%dStep, nfailsafe,&
            zpinch_getVariable, Rsolution, p_Rpredictor)
        
        ! Deallocate temporal memory
        call lsysbl_releaseVector(p_Rpredictor(1))
        call lsysbl_releaseVector(p_Rpredictor(2))
        deallocate(SfailsafeVariables, p_Rpredictor)
        
      else
        
        ! Compute and apply linearised FEM-FCT correction for
        ! hydrodynamic model
        call hydro_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(1),&
            ssectionNameHydro, rcollection)
        
        ! Copy correction factor for hydrodynamic model to transport model
        call afcstab_copyStabilisation(&
            rproblemLevel%Rafcstab(inviscidAFC),&
            rproblemLevel%Rafcstab(convectionAFC), AFCSTAB_DUP_EDGELIMITER)
        
        ! Apply linearised FEM-FCT correction for transport model
        call gfsc_buildConvectionVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrixTransport),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))
      end if
      
      
    case (5)   ! Transport model only
      
      if (nfailsafe .gt. 0) then

        ! Compute and apply linearised FEM-FCT correction for
        ! transport model
        call gfsc_buildConvectionVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrixTransport),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD-&
            AFCSTAB_FCTALGO_CORRECT, Rsolution(2))

        ! Copy correction factor for transport model to hydrodynamic model
        call afcstab_copyStabilisation(&
            rproblemLevel%Rafcstab(convectionAFC),&
            rproblemLevel%Rafcstab(inviscidAFC), AFCSTAB_DUP_EDGELIMITER)

        ! Apply failsafe flux correction for both models
        call afcstab_failsafeLimiting(rproblemLevel%Rafcstab(IposAFC),&
            rproblemLevel%Rmatrix(lumpedMassMatrixHydro),&
            SfailsafeVariables, rtimestep%dStep, nfailsafe,&
            zpinch_getVariable, Rsolution, p_Rpredictor)
        
        ! Deallocate temporal memory
        call lsysbl_releaseVector(p_Rpredictor(1))
        call lsysbl_releaseVector(p_Rpredictor(2))
        deallocate(SfailsafeVariables, p_Rpredictor)
        
      else
        
        ! Compute and apply linearised FEM-FCT correction for
        ! hydrodynamic model
        call gfsc_buildConvectionVectorFCT(&
            rproblemLevel%Rmatrix(lumpedMassMatrixTransport),&
            rproblemLevel%Rafcstab(convectionAFC), Rsolution(2),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_STANDARD+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(2))
        
        ! Copy correction factor for transport model to hydrodynamic system
        call afcstab_copyStabilisation(&
            rproblemLevel%Rafcstab(convectionAFC),&
            rproblemLevel%Rafcstab(inviscidAFC), AFCSTAB_DUP_EDGELIMITER)
        
        ! Apply linearised FEM-FCT correction  for transport model
        call hydro_calcCorrectionFCT(rproblemLevel, Rsolution(1),&
            rtimestep%dStep, .false., AFCSTAB_FCTALGO_CORRECT+&
            AFCSTAB_FCTALGO_SCALEBYMASS, Rsolution(1),&
            ssectionNameHydro, rcollection)
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
          rtimestep%dTime, hydro_calcBoundaryvalues1d)

    case (NDIM2D)
      call bdrf_filterVectorExplicit(rbdrCond(1), Rsolution(1),&
          rtimestep%dTime, hydro_calcBoundaryvalues2d)

    case (NDIM3D)
      call bdrf_filterVectorExplicit(rbdrCond(1), Rsolution(1),&
          rtimestep%dTime, hydro_calcBoundaryvalues3d)
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
    ! Hydrodynamic model, this subroutine accepts arrays of block vectors,
    ! e.g., to pass the solution from the hydrodynamic system and the
    ! scalar tracer equation simultaneously.
    !
    ! This subroutine is a wrapper which calls either the getVariable
    ! subroutine for the hydrodynamic system or returns the scalar tracer.
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
      call hydro_getVariable(RvectorBlock(1), cvariable, rvectorScalar)
    end if

  end subroutine zpinch_getVariable

end module zpinch_callback

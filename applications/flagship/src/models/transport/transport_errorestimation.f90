!##############################################################################
!# ****************************************************************************
!# <Name> transport_errorestimation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to perform
!# error estimation for scalar conservation laws in arbitrary spatial
!# dimensions.
!#
!# The following routines are available:
!#
!# 1.) transp_errestTargetFunc
!#     -> Estimates the error in the quantity of interest
!#
!# 2.) transp_errestRecovery
!#     -> Estimates the solution error using recovery techniques
!#
!# 3.) transp_errestExact
!#     -> Estimates the solution error using a given analytical solution
!#
!# 4.) transp_errestDispersionGHill
!#     -> Estimates the dispersion error for the rotating Gaussian hill
!#
!#
!# The following internal routines are available:
!#
!# 1.) transp_calcProtectionLayer
!#     -> Calculates a protection layer for grid adaptation
!#
!# </purpose>
!##############################################################################

module transport_errorestimation

#include "../../flagship.h"

!$use omp_lib
  use afcstabbase
  use collection
  use cubature
  use derivatives
  use fparser
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use pprocerror
  use pprocgradients
  use pprocindicator
  use problem
  use solveraux
  use stdoperators
  use storage
  use timestepaux
  use triangulation

  ! Modules from transport model
  use transport_basic
  use transport_callback
  use transport_callback2d

  implicit none

  private

  public :: transp_errestRecovery
  public :: transp_errestTargetFunc
  public :: transp_errestExact
  public :: transp_errestDispersionGHill

contains

  !*****************************************************************************

!<subroutine>

  subroutine transp_errestTargetFunc(rparlist, ssectionName,&
      rproblemLevel, rtimestep, rsolver, rsolutionPrimal,&
      rsolutionDual, rcollection, rtargetError, dtargetError, rrhs)

!<description>
    ! This subroutine estimates the error in the quantity of interest
!</description>

!<input>
    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! time-stepping algorithm
    type(t_timestep), intent(inout) :: rtimestep

    ! primal solution vector
    type(t_vectorBlock), intent(in), target :: rsolutionPrimal

    ! dual solution vector
    type(t_vectorBlock), intent(in) :: rsolutionDual

    ! OPTIONAL: right-hand side vector
    type(t_vectorBlock), intent(in), optional :: rrhs
!</input>

!<inputoutput>
    ! parameter list
    type(t_parlist), intent(inout) :: rparlist

    ! problem level structure
    type(t_problemLevel), intent(inout) :: rproblemLevel

    ! solver structure
    type(t_solver), intent(inout) :: rsolver

    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>

!<output>
    ! element-wise error distribution
    type(t_vectorScalar), intent(out) :: rtargetError

    ! global error in target qunatity
    real(DP), intent(out) :: dtargetError
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: serrorestimatorName
    character(LEN=SYS_STRLEN) :: sexacttargetfuncName
    character(LEN=SYS_STRLEN) :: sexactsolutionName
    character(LEN=SYS_STRLEN) :: stargetfuncName
    character(LEN=SYS_STRLEN) :: svelocityName

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorBlock) :: rvector1, rvector2
    type(t_vectorScalar) :: rvectorScalar
    type(t_matrixScalar) :: rmatrix
    type(t_collection) :: rcollectionTmp
    real(DP), dimension(:), pointer :: p_DlumpedMassMatrix, p_DtargetError
    real(DP), dimension(:), pointer :: p_DsolutionDual, p_Dresidual
    real(DP) :: dexactTargetError, dexactTargetFunc, dprotectLayerTolerance
    real(DP) :: daux, dtargetFunc, dStep, theta
    integer :: ieq, idim, convectionAFC, diffusionAFC, NEQ
    integer :: cconvectionStabilisation, cdiffusionStabilisation
    integer :: lumpedMassMatrix, templateMatrix, velocityfield
    integer :: itargetfunctype, iexactsolutiontype, imasstype
    integer :: nprotectLayers, igridindicator


    !---------------------------------------------------------------------------
    ! Goal-oriented error estimation part 1: Galerkin orthogonality error
    !---------------------------------------------------------------------------

    ! Create vector for Galerkin residual
    call lsysbl_createVectorBlock(rsolutionPrimal, rvector1, .true., ST_DOUBLE)
    call lsysbl_createVectorBlock(rsolutionPrimal, rvector2, .true., ST_DOUBLE)

    ! Ok, this is a little bit tricky. We need to compute the residual
    ! vector for the steady-state Galerkin scheme for the zeroth
    ! iteration. To this end, we switch off all types of
    ! stabilisation, mass matrices, time stepping parameter, etc., and
    ! force the velocity field, and hence, the preconditioner to be
    ! updated. Then the steady-state residual vector is evaluated.

    ! Get parameters from parameter list
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'convectionAFC', convectionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'diffusionAFC', diffusionAFC)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'imasstype', imasstype)

    ! Set mass type to 'no mass matrix'
    call parlst_setvalue(rparlist, ssectionName, 'imasstype', '0')

    ! Set time-stepping parameters
    dStep = rtimestep%dStep; rtimestep%dStep = 1.0_DP
    theta = rtimestep%theta; rtimestep%theta = 1.0_DP

    ! Set stabilisation to standard Galerkin
    cconvectionStabilisation =&
        rproblemLevel%Rafcstab(convectionAFC)%cafcstabType
    rproblemLevel%Rafcstab(convectionAFC)%cafcstabType = AFCSTAB_GALERKIN

    cdiffusionStabilisation =&
        rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType
    rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType = AFCSTAB_GALERKIN

    ! Set update notifiers for the discrete transport operator and the
    ! preconditioner in the problem level structure
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     TRANSP_TROPER_UPDATE)
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     TRANSP_PRECOND_UPDATE)

    ! Calculate the standard Galerkin preconditioner
    ! (required for rhs and residual calculation)
    call transp_calcPrecondThetaScheme(rproblemLevel, rtimestep,&
        rsolver, rsolutionPrimal, ssectionName, rcollection)

    ! Calculate the standard Galerkin right-hand side vector
    call transp_calcRhsThetaScheme(rproblemLevel, rtimestep, rsolver,&
        rsolutionPrimal, rvector1, ssectionName, rcollection, rrhs)

    ! Calculate the standard Galerkin residual
    call transp_calcResidualThetaScheme(rproblemLevel, rtimestep,&
        rsolver, rsolutionPrimal, rsolutionPrimal, rvector1,&
        rvector2, 0, ssectionName, rcollection)

    ! Ok, now we have to switch on all types of stabilisation again
    rproblemLevel%Rafcstab(convectionAFC)%cafcstabType =&
        cconvectionStabilisation
    rproblemLevel%Rafcstab(diffusionAFC)%cafcstabType =&
        cdiffusionStabilisation

    ! ... and we reset the mass type
    call parlst_setvalue(rparlist, ssectionName, 'imasstype',&
        trim(sys_si(imasstype, 3)))

    ! ... and the time-stepping structure
    rtimestep%dStep = dStep
    rtimestep%theta = theta

    ! Again, set update notifiers for the discrete transport operator
    ! and the preconditioner in the problem level structure
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     TRANSP_TROPER_UPDATE)
    rproblemLevel%iproblemSpec = ior(rproblemLevel%iproblemSpec,&
                                     TRANSP_PRECOND_UPDATE)

    ! We need the lumped mass matrix for scaling
    call parlst_getvalue_int(rparlist, ssectionName,&
        'lumpedMassMatrix', lumpedMassMatrix)

    if (lumpedMassMatrix > 0) then
      ! Set pointer to the lumped mass matrix
      call lsyssc_getbase_double(&
          rproblemLevel%Rmatrix(lumpedMassMatrix), p_DlumpedMassMatrix)
      NEQ = rproblemLevel%Rmatrix(lumpedMassMatrix)%NEQ

    else
      ! Compute the lumped mass matrix explicitly
      call parlst_getvalue_int(rparlist, ssectionName,&
          'templatematrix', templateMatrix)
      call lsyssc_duplicateMatrix(rproblemLevel%Rmatrix(templateMatrix),&
          rmatrix, LSYSSC_DUP_SHARE, LSYSSC_DUP_EMPTY)
      call stdop_assembleSimpleMatrix(rmatrix, DER_FUNC, DER_FUNC)
      call lsyssc_lumpMatrixScalar(rmatrix, LSYSSC_LUMP_DIAG)
      call lsyssc_getbase_double(rmatrix, p_DlumpedMassMatrix)
      NEQ = rmatrix%NEQ

    end if   ! lumpedMassMatrix > 0

    ! Set pointers
    call lsysbl_getbase_double(rvector2, p_Dresidual)
    call lsysbl_getbase_double(rsolutionDual, p_DsolutionDual)

    ! Now we compute the global error and its nodal contributions
    dtargetError = 0.0_DP
    do ieq = 1, NEQ
      dtargetError     = dtargetError + p_Dresidual(ieq)*p_DsolutionDual(ieq)
      p_Dresidual(ieq) = abs(p_Dresidual(ieq)*p_DsolutionDual(ieq))/&
                     p_DlumpedMassMatrix(ieq)
    end do
    dtargetError = abs(dtargetError)

    ! Compute the element-wise error.  We create a scalar vector and
    ! compute the local L1-norms of nodal error vector which yields
    ! the local values of the a posteriori error estimator.
    call lsyssc_createVector(rtargetError,&
        rproblemLevel%rtriangulation%NEL, .false.)
    call lsyssc_getbase_double(rtargetError, p_DtargetError)
    call pperr_scalar(PPERR_L1ERROR, daux, rvector2%RvectorBlock(1),&
        DelementError=p_DtargetError)

    ! Release temporal vectors
    call lsysbl_releaseVector(rvector1)
    call lsysbl_releaseVector(rvector2)

    ! Release temporal mass matrix if required
    if (lumpedMassMatrix .le. 0) call lsyssc_releaseMatrix(rmatrix)


    ! Check if an exact solution is available or if there is an exact
    ! expression for the target functional. If neither is available,
    ! then we are unable to compute the effectivity index.
    call parlst_getvalue_int(rparlist, ssectionName,&
        'iexactsolutiontype', iexactsolutiontype)
    call parlst_getvalue_string(rparlist, ssectionName,&
        'sexacttargetfuncname', sexacttargetfuncname, '')

    if ((iexactsolutiontype .ne. SOLUTION_ANALYTIC_POINTVALUE) .and.&
        (trim(sexacttargetfuncname) .eq. '')) then

      call output_lbrk()
      call output_line('Error Analysis')
      call output_line('--------------')
      call output_line('estimated error in quantity of interest: '//trim(sys_sdEP(dtargetError,15,6)))
      call output_line('exact error in quantity of interest:     '//'n.a.')
      call output_line('effectivity index:                       '//'n.a.')
      call output_line('relative effectivity index:              '//'n.a.')
      call output_lbrk()

    else

      !-------------------------------------------------------------------------
      ! Compute the effectivity index
      !-------------------------------------------------------------------------

      ! Get global configuration from parameter list
      call parlst_getvalue_int(rparlist, ssectionName,&
          'itargetfunctype', itargetfuncType)
      call parlst_getvalue_string(rparlist, ssectionName,&
          'sexactsolutionname', sexactsolutionname, '')
      
      select case(itargetfunctype)
      case (TFUNC_ZERO)
        call output_line('Zero target functional is not implemented yet!',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_errestTargetFunc')
        call sys_halt()


      case(TFUNC_VOLINTG)
        ! Get global configuration from parameter list
        call parlst_getvalue_string(rparlist, ssectionName,&
            'stargetfuncname', stargetfuncName)
        
        ! Get function parser from collection structure
        p_rfparser => collct_getvalue_pars(rcollection,&
            'rfparser', ssectionName=ssectionName)

        ! Initialise temporal collection structure
        call collct_init(rcollectionTmp)
        
        ! Prepare quick access arrays of the temporal collection structure
        rcollectionTmp%SquickAccess(1) = ''
        rcollectionTmp%SquickAccess(2) = 'rfparser'
        rcollectionTmp%DquickAccess(1) = rtimestep%dTime
        rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, sexactsolutionName)
        rcollectionTmp%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)

        ! Attach user-defined collection structure to temporal collection
        ! structure (may be required by the callback function)
        rcollectionTmp%p_rnextCollection => rcollection
        
        ! Attach function parser from boundary conditions to collection
        ! structure and specify its name in quick access string array
        call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)
        

        if (trim(sexacttargetfuncname) .ne. '') then
          ! Evaluate exact value of the quantity of interest
          call fparser_evalFunction(p_rfparser, sexacttargetfuncName,&
              (/rtimestep%dTime/), dexactTargetFunc)

          ! Compute the approximate value of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dtargetFunc,&
              rsolutionPrimal%RvectorBlock(1), rcollection=rcollectionTmp,&
              ffunctionWeight=transp_weightFuncAnalytic)
          
          ! Compute exact error in target functional. Note that the
          ! routine pperr_scalar computes the value $J(0-u_h)$ so that we
          ! have to change the sign to minus
          dexactTargetError = dexactTargetFunc+dtargetFunc
          
        else
          
          ! Compute the exact error of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dexactTargetError,&
              rsolutionPrimal%RvectorBlock(1), transp_refFuncAnalytic,&
              rcollectionTmp, transp_weightFuncAnalytic)
          
          ! Create empty vector
          call lsyssc_createVector(&
              rsolutionPrimal%RvectorBlock(1)%p_rspatialDiscr,&
              rvectorScalar, .true.)

          ! Compute the exact value of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dexactTargetFunc,&
              rvectorScalar, transp_refFuncAnalytic,&
              rcollectionTmp, transp_weightFuncAnalytic)

          ! Release temporal vector
          call lsyssc_releaseVector(rvectorScalar)
        end if

        ! Release temporal collection structure
        call collct_done(rcollectionTmp)


      case (TFUNC_SURFINTG)
        ! Get global configuration from parameter list
        call parlst_getvalue_int(rparlist, ssectionName,&
            'velocityfield', velocityfield)
        call parlst_getvalue_string(rparlist, ssectionName,&
            'stargetfuncname', stargetfuncName)

        ! Get function parser from collection structure
        p_rfparser => collct_getvalue_pars(rcollection,&
            'rfparser', ssectionName=ssectionName)
        
        ! Initialise temporal collection structure
        call collct_init(rcollectionTmp)

        ! Prepare quick access arrays of the temporal collection structure
        rcollectionTmp%SquickAccess(1) = ''
        rcollectionTmp%SquickAccess(2) = 'rfparser'
        rcollectionTmp%DquickAccess(1) = rtimestep%dTime
        rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, sexactsolutionName)
        rcollectionTmp%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)

        ! ... and also the numbers of the exact velocity function
        do idim = 1, rproblemLevel%rtriangulation%ndim
          call parlst_getvalue_string(rparlist, ssectionName,&
              'svelocityname', svelocityname, isubString=idim)
          rcollectionTmp%IquickAccess(idim+2) =&
              fparser_getFunctionNumber(p_rfparser, svelocityname)
        end do

        ! Attach primal solution vector and velocity fied to first and
        ! second quick access vector of the temporal collection structure
        rcollectionTmp%p_rvectorQuickAccess1 => rsolutionPrimal
        rcollectionTmp%p_rvectorQuickAccess2 => rproblemLevel%RvectorBlock(velocityfield)
        
        ! Attach user-defined collection structure to temporal collection
        ! structure (may be required by the callback function)
        rcollectionTmp%p_rnextCollection => rcollection
        
        ! Attach function parser from boundary conditions to collection
        ! structure and specify its name in quick access string array
        call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)
        

        if (trim(sexacttargetfuncname) .ne. '') then
          ! Evaluate exact value of the quantity of interest
          call fparser_evalFunction(p_rfparser, sexacttargetfuncName,&
              (/rtimestep%dTime/), dexactTargetFunc)
          
          ! Prepare quick access arrays of the temporal collection structure
          rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, '@null')
          
          ! Compute the approximate value of the quantity of interest
          call pperr_scalarBoundary2D(0, CUB_G3_1D, dtargetFunc,&
              ffunctionReference=transp_errorBdrInt2D_sim,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncBdrInt2D_sim)

          ! Compute exact error in target functional. Note that the
          ! routine pperr_scalar computes the value $J(0-u_h)$ so that we
          ! have to change the sign to minus
          dexactTargetError = dexactTargetFunc+dtargetFunc

        else

          ! Compute the exact error of the quantity of interest
          call pperr_scalarBoundary2D(0, CUB_G3_1D, dexactTargetError,&
              ffunctionReference=transp_errorBdrInt2D_sim,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncBdrInt2D_sim)

          ! Compute the exact value of the quantity of interest
          call pperr_scalarBoundary2D(0, CUB_G3_1D, dexactTargetFunc,&
              ffunctionReference=transp_refFuncBdrInt2D_sim,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncBdrInt2D_sim)
        end if

        ! Release temporal collection structure
        call collct_done(rcollectionTmp)


      case (TFUNC_MIXINTG)
        ! Get global configuration from parameter list
        call parlst_getvalue_int(rparlist, ssectionName,&
            'velocityfield', velocityfield)

        ! Get the name of the function used for evaluating the
        ! volume integral part of the target functional
        if (parlst_querysubstrings(rparlist, ssectionName, 'stargetfuncname') .eq. 0) then
          call parlst_getvalue_string(rparlist, ssectionName,&
              'stargetfuncname', stargetfuncname)
        else
          call parlst_getvalue_string(rparlist, ssectionName,&
              'stargetfuncname', stargetfuncname, isubstring=1)
        end if

        ! Get function parser from collection structure
        p_rfparser => collct_getvalue_pars(rcollection,&
            'rfparser', ssectionName=ssectionName)
        
        ! Initialise temporal collection structure
        call collct_init(rcollectionTmp)

        ! Prepare quick access arrays of the temporal collection structure
        rcollectionTmp%SquickAccess(1) = ''
        rcollectionTmp%SquickAccess(2) = 'rfparser'
        rcollectionTmp%DquickAccess(1) = rtimestep%dTime
        rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, sexactsolutionName)
        rcollectionTmp%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)
        
        ! ... and also the numbers of the exact velocity function
        do idim = 1, rproblemLevel%rtriangulation%ndim
          call parlst_getvalue_string(rparlist, ssectionName,&
              'svelocityname', svelocityname, isubString=idim)
          rcollectionTmp%IquickAccess(idim+2) =&
              fparser_getFunctionNumber(p_rfparser, svelocityname)
        end do

        ! Attach primal solution vector and velocity fied to first and
        ! second quick access vector of the temporal collection structure
        rcollectionTmp%p_rvectorQuickAccess1 => rsolutionPrimal
        rcollectionTmp%p_rvectorQuickAccess2 => rproblemLevel%RvectorBlock(velocityfield)

        ! Attach user-defined collection structure to temporal collection
        ! structure (may be required by the callback function)
        rcollectionTmp%p_rnextCollection => rcollection
        
        ! Attach function parser from boundary conditions to collection
        ! structure and specify its name in quick access string array
        call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)


        if (trim(sexacttargetfuncname) .ne. '') then
          ! Evaluate exact value of the quantity of interest
          call fparser_evalFunction(p_rfparser, sexacttargetfuncName,&
              (/rtimestep%dTime/), dexactTargetFunc)

          ! Prepare quick access arrays of the collection
          rcollectionTmp%IquickAccess(1) = fparser_getFunctionNumber(p_rfparser, '@null')

          ! Compute the approximate value of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dtargetFunc,&
              rsolutionPrimal%RvectorBlock(1), rcollection=rcollectionTmp,&
              ffunctionWeight=transp_weightFuncAnalytic)
          
          ! Get the name of the function used for evaluating the
          ! surface integral part of the target functional
          if (parlst_querysubstrings(rparlist, ssectionName, 'stargetfuncname') .eq. 0) then
            call parlst_getvalue_string(rparlist, ssectionName,&
                'stargetfuncname', stargetfuncname)
          else
            call parlst_getvalue_string(rparlist, ssectionName,&
                'stargetfuncname', stargetfuncname, isubstring=2)
          end if

          ! Prepare quick access arrays of the collection
          rcollectionTmp%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)

          ! Compute the approximate value of the quantity of interest
          call pperr_scalarBoundary2D(0, CUB_G3_1D, daux,&
              ffunctionReference=transp_errorBdrInt2D_sim,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncBdrInt2D_sim)

          ! Add boundary contribution
          dtargetFunc = dtargetFunc + daux

          ! Compute exact error in target functional. Note that the
          ! routine pperr_scalar computes the value $J(0-u_h)$ so that we
          ! have to change the sign to minus
          dexactTargetError = dexactTargetFunc+dtargetFunc

        else

          ! Compute the exact error of the quantity of interest
          call pperr_scalar(PPERR_MEANERROR, dexactTargetError,&
              rsolutionPrimal%RvectorBlock(1), transp_refFuncAnalytic,&
              rcollectionTmp, transp_weightFuncAnalytic)

          ! Create empty vector
          call lsyssc_createVector(&
              rsolutionPrimal%RvectorBlock(1)%p_rspatialDiscr,&
              rvectorScalar, .true.)

          ! Compute the exact value of the quantity of interest.
          call pperr_scalar(PPERR_MEANERROR, dexactTargetFunc,&
              rvectorScalar, transp_refFuncAnalytic,&
              rcollectionTmp, transp_weightFuncAnalytic)
          
          ! Release temporal vector
          call lsyssc_releaseVector(rvectorScalar)

          ! Get the name of the function used for evaluating the
          ! surface integral part of the target functional
          if (parlst_querysubstrings(rparlist, ssectionName, 'stargetfuncname') .eq. 0) then
            call parlst_getvalue_string(rparlist, ssectionName,&
                'stargetfuncname', stargetfuncname)
          else
            call parlst_getvalue_string(rparlist, ssectionName,&
                'stargetfuncname', stargetfuncname, isubstring=2)
          end if

          ! Prepare quick access arrays of the collection
          rcollectionTmp%DquickAccess(1) = rtimestep%dTime
          rcollectionTmp%IquickAccess(2) = fparser_getFunctionNumber(p_rfparser, stargetfuncName)

          ! Compute the exact error of the quantity of interest at the boundary
          call pperr_scalarBoundary2D(0, CUB_G3_1D, daux,&
              ffunctionReference=transp_errorBdrInt2D_sim,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncBdrInt2D_sim)

          ! Add boundary contribution
          dexactTargetError = dexactTargetError + daux

          ! Compute the exact value of the quantity of interest at the boundary
          call pperr_scalarBoundary2D(0, CUB_G3_1D, daux,&
              ffunctionReference=transp_refFuncBdrInt2D_sim,&
              rcollection=rcollectionTmp,&
              rdiscretisation=rsolutionPrimal%RvectorBlock(1)%p_rspatialdiscr,&
              ffunctionWeight=transp_weightFuncBdrInt2D_sim)

          ! Add boundary contribution
          dexactTargetFunc = dexactTargetFunc + daux
        end if

        ! Release temporal collection structure
        call collct_done(rcollectionTmp)


      case default
        call output_line('Invalid type of target functional!',&
            OU_CLASS_ERROR,OU_MODE_STD,'transp_errestTargetFunc')
        call sys_halt()
      end select

      call output_lbrk()
      call output_line('Error Analysis')
      call output_line('--------------')
      call output_line('estimated error in quantity of interest: '//trim(sys_sdEP(dtargetError,15,6)))
      call output_line('exact error in quantity of interest:     '//trim(sys_sdEP(abs(dexactTargetError),15,6)))
      call output_line('effectivity index:                       '//trim(sys_sdEP(dtargetError/abs(dexactTargetError),15,6)))
      call output_line('relative effectivity index:              '//trim(sys_sdEP(abs( (dtargetError-abs(dexactTargetError)) /&
                                                                                        dexactTargetFunc ),15,6)))
      call output_lbrk()

    end if


    !---------------------------------------------------------------------------
    ! Apply the adaptation strategy
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist, ssectionName,&
        'errorestimator', serrorestimatorName)
    call parlst_getvalue_int(rparlist,&
        trim(serrorestimatorName), 'igridindicator', igridindicator)

    ! What type of grid indicator are we?
    select case(igridIndicator)

    case (ERREST_ASIS)
      ! That is simple, do nothing.

    case default
      call output_line('Invalid type of grid indicator!',&
          OU_CLASS_ERROR,OU_MODE_STD,'transp_errestTargetFunc')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Calculate protection layers
    !---------------------------------------------------------------------------

    call parlst_getvalue_string(rparlist,&
        ssectionName, 'errorestimator', serrorestimatorName)
    call parlst_getvalue_int(rparlist,&
        trim(serrorestimatorName), 'nprotectLayers', nprotectLayers, 0)
    call parlst_getvalue_double(rparlist, trim(serrorestimatorName),&
        'dprotectLayerTolerance', dprotectLayerTolerance, 0.0_DP)

    if (nprotectLayers > 0)&
        call transp_calcProtectionLayer(rproblemLevel%rtriangulation,&
        nprotectLayers, dprotectLayerTolerance, rtargetError)

  end subroutine transp_errestTargetFunc

  !*****************************************************************************

!<subroutine>

  subroutine transp_errestRecovery(rparlist, ssectionName,&
      rproblemLevel, rsolution, dtime, rerror, derror, rcollection)

!<description>
    ! This subroutine estimates the error of the discrete solution by
    ! using recovery procedures such as the superconvergent patch
    ! recovery technique or L2-projection. If an exact solution value
    ! is avaialable, it computes the effectivity index of the error
    ! estimator. Moreover it applies the specified strategy to convert
    ! the error estimator into a grid indicator for mesh adaptation.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! simulation time
    real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>

!<output>
    ! element-wise error distribution
    type(t_vectorScalar), intent(out) :: rerror

    ! global error
    real(DP), intent(out) :: derror
!</output>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: serrorestimatorName
    character(LEN=SYS_STRLEN) :: sexactsolutionName

    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_vectorScalar) :: rvectorScalar
    type(t_collection) :: rcollectionTmp
    real(DP), dimension(:), pointer :: p_Ddata, p_Derror, p_DelementError
    real(DP) :: dnoiseFilter, dabsFilter, dsolution, dvalue,&
                dexacterror, dprotectLayerTolerance
    integer :: i, ierrorEstimator, igridindicator, iexactsolutiontype
    integer :: nprotectLayers, iexactsolution, nexactsolution


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'errorestimator', serrorestimatorName)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'iexactsolutiontype', iexactsolutiontype, 0)
    call parlst_getvalue_int(rparlist,&
        trim(serrorestimatorName), 'ierrorestimator', ierrorestimator)
    call parlst_getvalue_int(rparlist,&
        trim(serrorestimatorName), 'igridindicator', igridindicator)
    call parlst_getvalue_int(rparlist,&
        trim(serrorestimatorName), 'nprotectLayers', nprotectLayers, 0)
    call parlst_getvalue_double(rparlist,&
        trim(serrorestimatorName), 'dprotectLayerTolerance',&
        dprotectLayerTolerance, 0.0_DP)
    nexactsolution = parlst_querysubstrings(rparlist,&
        ssectionName, 'sexactsolutionname')


    !---------------------------------------------------------------------------
    ! Perform recovery-based error estimation
    !---------------------------------------------------------------------------

    ! What type of error estimator are we?
    select case(ierrorEstimator)

    case (ERREST_L2PROJECTION)
      call lsyssc_createVector(rerror, rproblemLevel%rtriangulation%NEL, .false.)
      call ppgrd_calcGradientError(rsolution%RvectorBlock(1), derror,&
                                   PPGRD_INTERPOL, 0, rerror)

    case (ERREST_SPR_VERTEX)
      call lsyssc_createVector(rerror, rproblemLevel%rtriangulation%NEL, .false.)
      call ppgrd_calcGradientError(rsolution%RvectorBlock(1), derror,&
                                   PPGRD_ZZTECHNIQUE, PPGRD_NODEPATCH, rerror)

    case (ERREST_SPR_ELEMENT)
      call lsyssc_createVector(rerror, rproblemLevel%rtriangulation%NEL, .false.)
      call ppgrd_calcGradientError(rsolution%RvectorBlock(1), derror,&
                                   PPGRD_ZZTECHNIQUE, PPGRD_ELEMPATCH, rerror)

    case (ERREST_SPR_FACE)
      call lsyssc_createVector(rerror, rproblemLevel%rtriangulation%NEL, .false.)
      call ppgrd_calcGradientError(rsolution%RvectorBlock(1), derror,&
                                   PPGRD_ZZTECHNIQUE, PPGRD_FACEPATCH, rerror)

    case (ERREST_LIMAVR)
      call lsyssc_createVector(rerror, rproblemLevel%rtriangulation%NEL, .false.)
      call ppgrd_calcGradientError(rsolution%RvectorBlock(1), derror,&
                                   PPGRD_LATECHNIQUE, 0, rerror)

    case (ERREST_SECONDDIFF)
      call parlst_getvalue_double(rparlist,&
          trim(serrorestimatorName), 'dnoiseFilter', dnoiseFilter)
      call parlst_getvalue_double(rparlist,&
          trim(serrorestimatorName), 'dabsFilter', dabsFilter)
      call ppind_secondDifference(rsolution%RvectorBlock(1),&
          dnoiseFilter, dabsFilter, rerror)

      derror = 1.0

    case default
      call output_line('Invalid type of error estimator!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_errestRecovery')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Compute the effectivity index
    !---------------------------------------------------------------------------

    select case(iexactsolutiontype)
    case (SOLUTION_ANALYTIC_POINTVALUE)

      ! Get function parser from collection
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)

      ! Initialise temporal collection structure
      call collct_init(rcollectionTmp)

      ! Prepare quick access arrays of the temporal collection structure
      rcollectionTmp%SquickAccess(1) = ''
      rcollectionTmp%SquickAccess(2) = 'rfparser'
      rcollectionTmp%DquickAccess(1) = dtime
      
      do iexactsolution = 1, nexactsolution
        call parlst_getvalue_string(rparlist, ssectionName,&
            'sexactsolutionname', sexactsolutionName, '',&
            isubstring=iexactsolution)
        rcollectionTmp%IquickAccess(iexactsolution) =&
            fparser_getFunctionNumber(p_rfparser, sexactsolutionName)
      end do

      ! Make sure that the callback routine will throw an error if a
      ! derivative is required which is not given by the user.
      do iexactsolution = nexactsolution+1, DER_DERIV3D_Z
        rcollectionTmp%IquickAccess(iexactsolution) = 0
      end do
      
      ! Attach user-defined collection structure to temporal collection
      ! structure (may be required by the callback function)
      rcollectionTmp%p_rnextCollection => rcollection

      ! Attach function parser from boundary conditions to collection
      ! structure and specify its name in quick access string array
      call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)

      ! Calculate the H1-error of the reference solution
      call pperr_scalar(PPERR_H1ERROR, dexacterror, rsolution%RvectorBlock(1),&
                        transp_refFuncAnalytic, rcollectionTmp)

      call output_lbrk()
      call output_line('Error Analysis')
      call output_line('--------------')
      call output_line('estimated H1-error: '//trim(sys_sdEP(derror,15,6)))
      call output_line('exact H1-error:     '//trim(sys_sdEP(dexacterror,15,6)))
      call output_line('effectivity index:  '//trim(sys_sdEP(derror/dexacterror,15,6)))
      call output_lbrk()

      ! Release temporal collection structure
      call collct_done(rcollectionTmp)

    case default
      call output_lbrk()
      call output_line('Error Analysis')
      call output_line('--------------')
      call output_line('estimated H1-error: '//trim(sys_sdEP(derror,15,6)))
      call output_lbrk()

    end select


    !---------------------------------------------------------------------------
    ! Apply the adaptation strategy
    !---------------------------------------------------------------------------

    ! What type of grid indicator are we?
    select case(igridIndicator)

    case (ERREST_ASIS)
      ! That is simple, do nothing.


    case (ERREST_EQUIDIST)
      ! Equidistribute the relative percentage error
      call output_lbrk()
      call output_line('Equidistribution strategy')
      call output_line('-------------------------')

      ! We need the global norm of the scalar error variable
      call pperr_scalar(PPERR_H1ERROR, dsolution,&
          rsolution%RvectorBlock(1))

      call output_line('Total percentage error:       '//trim(sys_sdEP(derror/dsolution,15,6)))

      ! Compute global permissible element error
      dvalue = sqrt(dsolution**2 + derror**2)/sqrt(real(rerror%NEQ, DP))

      call output_line('Permissible percentage error: '//trim(sys_sdEP(dvalue,15,6))//' x TOL')
      call output_lbrk()

      ! Scale element error by global permissible error
      call lsyssc_scaleVector(rerror, 1.0_DP/dvalue)


    case (-ERREST_EQUIDIST)
      ! Equidistribute the relative percentage error
      call output_lbrk()
      call output_line('Equidistribution strategy')
      call output_line('-------------------------')

      ! We need the global norm of the scalar error variable
      call lsyssc_createVector(rvectorScalar, rerror%NEQ, .true.)
      call lsyssc_getbase_double(rvectorScalar, p_DelementError)
      call pperr_scalar(PPERR_H1ERROR, dsolution,&
          rsolution%RvectorBlock(1), DelementError=p_DelementError)

      call output_line('Total percentage error: '//trim(sys_sdEP(derror/dsolution,15,6)))
      call output_lbrk()

      ! Compute local permissible element error for each elements
      call lsyssc_getbase_double(rerror, p_Derror)
      call lsyssc_getbase_double(rvectorScalar, p_Ddata)
      dvalue = dsolution/sqrt(real(rerror%NEQ, DP))

      do i = 1, size(p_Derror)
        p_Derror(i) = p_Derror(i)/(0.5*(sqrt(p_Ddata(i)) + dvalue))
      end do


    case (ERREST_LOGEQUIDIST)
      ! Equidistribute the logarithmic error
      call output_lbrk()
      call output_line('Logarithmic equidistribution strategy')
      call output_line('-------------------------------------')

      ! Set pointer
      call lsyssc_getbase_double(rerror, p_Ddata)

      ! Determine largest error value
      dvalue = -SYS_MAXREAL_DP
      do i = 1, size(p_Ddata)
        dvalue = max(dvalue, p_Ddata(i))
      end do

      call output_line('Maximum error: '//trim(sys_sdEP(derror/dvalue,15,6)))

      ! Normalise error by largest value
      call lsyssc_scaleVector(rerror, 1.0_DP/dvalue)

      ! Initialise mean value
      dvalue = 0.0_DP

      ! Loop over all contributions
      do i = 1, size(p_Ddata)
        p_Ddata(i) = log(max(exp(-20.0_DP), p_Ddata(i)))
        dvalue   = dvalue + p_Ddata(i)
      end do

      ! Calculate mean
      dvalue = dvalue/real(size(p_Ddata), DP)

      call output_line('Mean value:    '//trim(sys_sdEP(derror/dvalue,15,6)))
      call output_lbrk()

      ! Subtract mean value from grid indicator
      do i = 1, size(p_Ddata)
        p_Ddata(i) = p_Ddata(i)-dvalue
      end do


    case (ERREST_AUTORMS)
      ! Automatic treshold for RMS
      call output_lbrk()
      call output_line('Automatic treshold for RMS')
      call output_line('--------------------------')

      ! Set pointer
      call lsyssc_getbase_double(rerror, p_Ddata)

      ! Initialise mean value
      dvalue = 0.0_DP

      ! Loop over all  contributions
      do i = 1, size(p_Ddata)
        dvalue = dvalue + p_Ddata(i)**2
      end do

      ! Calculate root mean value
      dvalue = sqrt(dvalue/real(size(p_Ddata), DP))

      call output_line('RMS value: '//trim(sys_sdEP(derror/dvalue,15,6)))
      call output_lbrk()

      ! Normalise grid indicator by RMS
      call lsyssc_scaleVector(rerror, 1.0_DP/dvalue)


    case default
      call output_line('Invalid type of grid indicator!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'transp_errestRecovery')
      call sys_halt()
    end select


    !---------------------------------------------------------------------------
    ! Calculate protection layers
    !---------------------------------------------------------------------------

    if (nprotectLayers > 0)&
        call transp_calcProtectionLayer(rproblemLevel%rtriangulation,&
        nprotectLayers, dprotectLayerTolerance, rerror)

  end subroutine transp_errestRecovery

  !*****************************************************************************

!<subroutine>

  subroutine transp_errestExact(rparlist, ssectionName,&
      rproblemLevel, rsolution, dtime, rerrorL1, derrorL1,&
      rerrorL2, derrorL2, rerrorH1, derrorH1, rcollection)

!<description>
    ! This subroutine estimates the error of the discrete solution by
    ! comparing it to the analytically given exact solution.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! simulation time
    real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>

!<output>
    ! OPTIONAL: element-wise L1-error distribution
    type(t_vectorScalar), intent(out), optional :: rerrorL1

    ! OPTIONAL: element-wise L2-error distribution
    type(t_vectorScalar), intent(out), optional :: rerrorL2

    ! OPTIONAL: element-wise H1-error distribution
    type(t_vectorScalar), intent(out), optional :: rerrorH1

    ! OPTIONAL: global L1-error
    real(DP), intent(out), optional :: derrorL1

    ! OPTIONAL: global L2-error
    real(DP), intent(out), optional :: derrorL2

    ! OPTIONAL: global H1-error
    real(DP), intent(out), optional :: derrorH1
!</output>
!</subroutine>
    
    ! section names
    character(LEN=SYS_STRLEN) :: serrorestimatorName
    character(LEN=SYS_STRLEN) :: sexactsolutionName
    
    ! local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_collection) :: rcollectionTmp
    real(DP), dimension(:), pointer :: p_Derror
    integer :: iexactsolutiontype, iexactsolution, nexactsolution

    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'errorestimator', serrorestimatorName)
    call parlst_getvalue_int(rparlist,&
        ssectionName, 'iexactsolutiontype', iexactsolutiontype, 0)
    nexactsolution = parlst_querysubstrings(rparlist,&
        ssectionName, 'sexactsolutionname')

    !---------------------------------------------------------------------------
    ! Compute solution error
    !---------------------------------------------------------------------------

    select case(iexactsolutiontype)
    case (SOLUTION_ANALYTIC_POINTVALUE)
      
      ! Get function parser from collection
      p_rfparser => collct_getvalue_pars(rcollection,&
          'rfparser', ssectionName=ssectionName)
      
      ! Initialise temporal collection structure
      call collct_init(rcollectionTmp)
      
      ! Prepare quick access arrays of the temporal collection structure
      rcollectionTmp%SquickAccess(1) = ''
      rcollectionTmp%SquickAccess(2) = 'rfparser'
      rcollectionTmp%DquickAccess(1) = dtime

      do iexactsolution = 1, nexactsolution
        call parlst_getvalue_string(rparlist, ssectionName,&
            'sexactsolutionname', sexactsolutionName, '',&
            isubstring=iexactsolution)
        rcollectionTmp%IquickAccess(iexactsolution) =&
            fparser_getFunctionNumber(p_rfparser, sexactsolutionName)
      end do

      ! Make sure that the callback routine will throw an error if a
      ! derivative is required which is not given by the user.
      do iexactsolution = nexactsolution+1, DER_DERIV3D_Z
        rcollectionTmp%IquickAccess(iexactsolution) = 0
      end do
      
      ! Attach user-defined collection structure to temporal collection
      ! structure (may be required by the callback function)
      rcollectionTmp%p_rnextCollection => rcollection
      
      ! Attach function parser from boundary conditions to collection
      ! structure and specify its name in quick access string array
      call collct_setvalue_pars(rcollectionTmp, 'rfparser', p_rfparser, .true.)
      
      call output_lbrk()
      call output_separator(OU_SEP_DOLLAR, OU_CLASS_MSG, OU_MODE_STD+OU_MODE_BENCHLOG)
      call output_line('Error Analysis', OU_CLASS_MSG, OU_MODE_STD+OU_MODE_BENCHLOG)
      call output_separator(OU_SEP_MINUS, OU_CLASS_MSG, OU_MODE_STD+OU_MODE_BENCHLOG)

      if (present(derrorL1)) then
        ! Calculate the L1-error of the reference solution
        if (present(rerrorL1)) then
          call lsyssc_createVector(rerrorL1,&
              rproblemLevel%rtriangulation%NEL, .false.)
          call lsyssc_getbase_double(rerrorL1, p_Derror)
          call pperr_scalar(PPERR_L1ERROR, derrorL1, rsolution%RvectorBlock(1),&
              transp_refFuncAnalytic, rcollectionTmp, DelementError=p_Derror)
        else
          call pperr_scalar(PPERR_L1ERROR, derrorL1, rsolution%RvectorBlock(1),&
              transp_refFuncAnalytic, rcollectionTmp)
        end if
        call output_line('Exact L1-error: '//trim(sys_sdEP(derrorL1,15,6)),&
            OU_CLASS_MSG, OU_MODE_STD+OU_MODE_BENCHLOG)
      end if

      if (present(derrorL2)) then
        ! Calculate the L2-error of the reference solution
        if (present(rerrorL2)) then
          call lsyssc_createVector(rerrorL2,&
              rproblemLevel%rtriangulation%NEL, .false.)
          call lsyssc_getbase_double(rerrorL2, p_Derror)
          call pperr_scalar(PPERR_L2ERROR, derrorL2, rsolution%RvectorBlock(1),&
              transp_refFuncAnalytic, rcollectionTmp, DelementError=p_Derror)
        else
          call pperr_scalar(PPERR_L2ERROR, derrorL2, rsolution%RvectorBlock(1),&
              transp_refFuncAnalytic, rcollectionTmp)
        end if
        call output_line('Exact L2-error: '//trim(sys_sdEP(derrorL2,15,6)),&
            OU_CLASS_MSG, OU_MODE_STD+OU_MODE_BENCHLOG)
      end if

      if (present(derrorH1)) then
        ! Calculate the H1-error of the reference solution
        if (present(rerrorH1)) then
          call lsyssc_createVector(rerrorH1,&
              rproblemLevel%rtriangulation%NEL, .false.)
          call lsyssc_getbase_double(rerrorH1, p_Derror)
          call pperr_scalar(PPERR_H1ERROR, derrorH1, rsolution%RvectorBlock(1),&
              transp_refFuncAnalytic, rcollectionTmp, DelementError=p_Derror)
        else
          call pperr_scalar(PPERR_H1ERROR, derrorH1, rsolution%RvectorBlock(1),&
              transp_refFuncAnalytic, rcollectionTmp)
        end if
        call output_line('Exact H1-error: '//trim(sys_sdEP(derrorH1,15,6)),&
            OU_CLASS_MSG, OU_MODE_STD+OU_MODE_BENCHLOG)
      end if

      call output_separator(OU_SEP_DOLLAR, OU_CLASS_MSG, OU_MODE_STD+OU_MODE_BENCHLOG)
      call output_lbrk()

    case default
      call output_line('Analytical solution is not available!',&
          OU_CLASS_WARNING,OU_MODE_STD,'transp_errestExact')
      if (present(derrorL1)) derrorL1 = -1.0_DP
      if (present(derrorL2)) derrorL2 = -1.0_DP
      if (present(derrorH1)) derrorH1 = -1.0_DP
    end select
   
  end subroutine transp_errestExact

  !*****************************************************************************

!<subroutine>

  subroutine transp_errestDispersionGHill(rparlist, ssectionName,&
      rproblemLevel, rsolution, dtime, derror, rcollection)

!<description>
    ! This subroutine estimates the dispersion error of the discrete solution
    ! to the rotating Gaussian hill benchmark.
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list and collection structure
    character(LEN=*), intent(in) :: ssectionName

    ! solution vector
    type(t_vectorBlock), intent(in) :: rsolution

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! simulation time
    real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
    ! collection structure
    type(t_collection), intent(inout), target :: rcollection
!</inputoutput>

!<output>
    ! dispersion error
    real(DP), intent(out) :: derror
!</output>

    ! Local variables
    type(t_fparser), pointer :: p_rfparser
    type(t_collection) :: rcollectionTmp
    character(LEN=SYS_STRLEN) :: sdiffusionName
    real(DP), dimension(1) :: Dunity = (/1.0_DP/)
    real(DP) :: dxhat,dyhat,dalpha

    ! Initialise temporal collection structure
    call collct_init(rcollectionTmp)
        
    ! Compute integral quantities
    rcollectionTmp%IquickAccess(1) = 1
    call pperr_scalar(PPERR_MEANERROR, dxhat, rsolution%RvectorBlock(1),&
        rcollection=rcollectionTmp, ffunctionWeight=transp_weightFuncGHill)
    
    rcollectionTmp%IquickAccess(1) = 2
    call pperr_scalar(PPERR_MEANERROR, dyhat, rsolution%RvectorBlock(1),&
        rcollection=rcollectionTmp, ffunctionWeight=transp_weightFuncGHill)

    rcollectionTmp%IquickAccess(1) = 3
    rcollectionTmp%DquickAccess(1) = dxhat
    rcollectionTmp%DquickAccess(2) = dyhat
    call pperr_scalar(PPERR_MEANERROR, derror, rsolution%RvectorBlock(1),&
        rcollection=rcollectionTmp, ffunctionWeight=transp_weightFuncGHill)

    ! Retrieve name/number of expression describing the diffusion coefficient
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'sdiffusionname', sdiffusionName)

    ! Get function parser from collection
    p_rfparser => collct_getvalue_pars(rcollection,&
        'rfparser', ssectionName=ssectionName)
    
    ! Evaluate the constant coefficient from the function parser
    call fparser_evalFunction(p_rfparser, sdiffusionName, Dunity, dalpha)
    
    ! Compute dispersion error
    derror = derror / (4*dtime*dalpha) - 1.0_DP
    
    call output_line('Dispersion-error: '//trim(sys_sdEP(derror,15,6)))

  end subroutine transp_errestDispersionGHill

  !*****************************************************************************

!<subroutine>

  subroutine transp_calcProtectionLayer(rtriangulation, nlayers, dtreshold, rvector)

!<description>
    ! This subroutine takes the given vector rvector as elementwise
    ! indicator and calculates a prescribed number of protection
    ! layers for all markers which exceed the given threshold dtreshold.
!</description>

!<input>
    ! Triangulation
    type(t_triangulation), intent(in) :: rtriangulation

    ! Number of protection layers
    integer, intent(in) :: nlayers

    ! Threshold of protection layer
    real(DP), intent(in) :: dtreshold
!</input>

!<inputoutput>
    ! Vector with elementwise markers
    type(t_vectorScalar), intent(inout) :: rvector
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer, dimension(:,:), pointer :: p_IverticesAtElement,p_IneighboursAtElement
    logical, dimension(:), pointer :: p_BisactiveElement
    integer :: h_BisactiveElement
    integer :: ilayer

    ! Create auxiliary memory
    h_BisactiveElement = ST_NOHANDLE
    call storage_new('transp_calcProtectionLayer',' BisactiveElement',&
        rtriangulation%NEL, ST_LOGICAL, h_BisactiveElement, ST_NEWBLOCK_NOINIT)
    call storage_getbase_logical(h_BisactiveElement, p_BisactiveElement)
    
    ! Set pointers
    call storage_getbase_int2D(rtriangulation%h_IneighboursAtElement,&
                               p_IneighboursAtElement)
    call storage_getbase_int2D(rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
    call lsyssc_getbase_double(rvector, p_Ddata)
    
    ! Compute protection layers
    do ilayer = 1, nlayers
      
      ! Reset activation flag
      p_BisActiveElement = .false.
      
      ! Compute a single-width protection layer
      call doProtectionLayerUniform(p_IverticesAtElement, p_IneighboursAtElement,&
          rtriangulation%NEL, dtreshold, p_Ddata, p_BisActiveElement)
    end do
    
    ! Release memory
    call storage_free(h_BisactiveElement)
    
  contains
    
    !**************************************************************
    ! Compute one uniformly distributed protection layer
    
    subroutine doProtectionLayerUniform(IverticesAtElement,&
        IneighboursAtElement, NEL, dthreshold, Ddata, BisactiveElement)
      
      integer, dimension(:,:), intent(in) :: IverticesAtElement
      integer, dimension(:,:), intent(in) :: IneighboursAtElement
      real(DP), intent(in) :: dthreshold
      integer, intent(in) :: NEL
      
      real(DP), dimension(:), intent(inout) :: Ddata
      logical, dimension(:), intent(inout) :: BisactiveElement
      
      
      ! local variables
      integer :: iel,ive,jel
      
      ! Loop over all elements in triangulation
      do iel = 1, NEL
        
        ! Do nothing if element belongs to active layer
        if (BisactiveElement(iel)) cycle

        ! Do nothing if element indicator does not exceed threshold
        if (Ddata(iel) .le. dthreshold) cycle

        ! Loop over neighbouring elements
        do ive = 1, tria_getNVE(IverticesAtElement, iel)

          ! Get number of neighbouring element
          jel = IneighboursAtElement(ive, iel)

          ! Do nothing at the boundary
          if (jel .eq. 0) cycle

          ! Check if element belongs to active layer
          if (BisactiveElement(jel)) then
            ! If yes, then just update the element indicator
            Ddata(jel) = max(Ddata(jel), Ddata(iel))
          else
            ! Otherwise, we have to check if the neighbouring element
            ! exceeds the prescribed threshold level. If this is the case
            ! it will be processed later or has already been processed
            if (Ddata(jel) .lt. dthreshold) then
              Ddata(jel) = max(Ddata(jel), Ddata(iel))
              BisactiveElement(jel) = .true.
            end if
          end if
        end do
      end do

    end subroutine doProtectionLayerUniform
    
  end subroutine transp_calcProtectionLayer
  
end module transport_errorestimation

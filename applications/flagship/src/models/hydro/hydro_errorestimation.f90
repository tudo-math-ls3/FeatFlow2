!##############################################################################
!# ****************************************************************************
!# <Name> hydro_errorestimation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all routines which are required to perform
!# error estimation for the compressible Euler/Navier-Stokes equations.
!#
!# The following routines are available:
!#
!# 1.) hydro_errestRecovery
!#     -> Estimates the solution error using recovery techniques
!#
!#
!# The following internal routines are available:
!#
!# 1.) hydro_calcProtectionLayer
!#     -> Calculates a protection layer for grid adaptation
!#
!# </purpose>
!##############################################################################

module hydro_errorestimation

#include "../../flagship.h"

!$use omp_lib
  use collection
  use fsystem
  use genoutput
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use pprocgradients
  use pprocindicator
  use problem
  use storage
  use triangulation

  ! Modules from hydrodynamic model
  use hydro_basic

  implicit none

  private

  public :: hydro_errestRecovery

contains

  !*****************************************************************************

!<subroutine>

  subroutine hydro_errestRecovery(rparlist, ssectionName,&
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

    ! section name in parameter list
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
    character(LEN=SYS_STRLEN) :: sindatfileName
    character(LEN=SYS_STRLEN) :: serrorestimatorName
    character(LEN=SYS_STRLEN) :: sexactsolutionName

    ! local variables
    type(t_vectorScalar) :: rvectorScalar, rvectorTmp
    real(DP), dimension(:), pointer :: p_Ddata, p_DdataTmp
    character(LEN=SYS_STRLEN) :: serrorvariable
    real(DP) :: dnoiseFilter, dabsFilter, dvalue,&
                dprotectLayerTolerance, derrorTmp
    integer :: ierrorEstimator, igridindicator, iexactsolutiontype
    integer :: nprotectLayers, ierrorVariable, nerrorVariables


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'indatfile', sindatfileName)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'errorestimator', serrorestimatorName)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'sexactsolutionname', sexactsolutionName, '')
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


    !---------------------------------------------------------------------------
    ! Perform recovery-based error estimation
    !---------------------------------------------------------------------------

    nerrorVariables = parlst_querysubstrings(rparlist,&
        trim(serrorestimatorName), 'serrorvariable')

    ! Loop over all error variables
    do ierrorVariable = 1, nerrorVariables

      ! Get name of error variable
      call parlst_getvalue_string(rparlist, trim(serrorestimatorName),&
          'serrorvariable', serrorVariable, isubString=ierrorVariable)

      ! Extract scalar variable from vector of conservative variables
      call hydro_getVariable(rsolution, serrorVariable, rvectorScalar)

      ! What type of error estimator are we?
      select case(ierrorEstimator)

      case (ERREST_L2PROJECTION)
        call lsyssc_createVector(rvectorTmp,&
            rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
            PPGRD_INTERPOL, 0, rvectorTmp)

      case (ERREST_SPR_VERTEX)
        call lsyssc_createVector(rvectorTmp,&
            rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
            PPGRD_ZZTECHNIQUE, PPGRD_NODEPATCH, rvectorTmp)

      case (ERREST_SPR_ELEMENT)
        call lsyssc_createVector(rvectorTmp,&
            rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
            PPGRD_ZZTECHNIQUE, PPGRD_ELEMPATCH, rvectorTmp)

      case (ERREST_SPR_FACE)
        call lsyssc_createVector(rvectorTmp,&
            rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
            PPGRD_ZZTECHNIQUE, PPGRD_FACEPATCH, rvectorTmp)

      case (ERREST_LIMAVR)
        call lsyssc_createVector(rvectorTmp,&
            rproblemLevel%rtriangulation%NEL, .false.)
        call ppgrd_calcGradientError(rvectorScalar, derrorTmp,&
            PPGRD_LATECHNIQUE, 0, rvectorTmp)

      case (ERREST_SECONDDIFF)
        call parlst_getvalue_double(rparlist,&
            trim(serrorestimatorName), 'dnoiseFilter', dnoiseFilter)
        call parlst_getvalue_double(rparlist,&
            trim(serrorestimatorName), 'dabsFilter', dabsFilter)
        call ppind_secondDifference(rvectorScalar, dnoiseFilter,&
            dabsFilter, rvectorTmp)

        ! This is no error estimator
        derrorTmp = 1.0

      case default
        call output_line('Invalid type of error estimator!',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_errestRecovery')
        call sys_halt()
      end select


      ! Compute the root-mean square value
      call lsyssc_getbase_double(rvectorTmp, p_Ddata)
      dvalue = sqrt(sum(p_Ddata*p_Ddata)/real(rvectorTmp%NEQ, DP))
      if (abs(dvalue) .gt. SYS_EPSREAL_DP) then
        dvalue = 1.0_DP/dvalue
        call lsyssc_scaleVector(rvectorTmp, dvalue)
      end if

      ! Compute the global and element-wise error
      if (ierrorVariable .eq. 1) then

        ! Initialise the global error
        derror = derrorTmp

        ! Initialise the element-wise error
        call lsyssc_copyVector(rvectorTmp, rerror)
      else

        ! Update the global error
        derror = max(derror, derrorTmp)

        ! Update the element-wise error
        call lsyssc_getbase_double(rvectorTmp, p_DdataTmp)
        call lsyssc_getbase_double(rerror, p_Ddata)
        p_Ddata = max(p_Ddata, p_DdataTmp)
!!$        call lsyssc_vectorLinearComb(rvectorTmp, rerror, 1.0_DP, 1.0_DP)
      end if

      ! Release scalar variable and temporal error
      call lsyssc_releaseVector(rvectorScalar)
      call lsyssc_releaseVector(rvectorTmp)
    end do

!!$    ! Scale the global and element-wise error by the number of error variables
!!$    if (nerrorVariables .gt. 1) then
!!$      dvalue = 1.0_DP/real(nerrorVariables, DP)
!!$      derror = derror*dvalue
!!$      call lsyssc_scaleVector(rerror, dvalue)
!!$    end if


    !---------------------------------------------------------------------------
    ! Calculate protection layers
    !---------------------------------------------------------------------------

    if (nprotectLayers > 0)&
        call hydro_calcProtectionLayer(rproblemLevel%rtriangulation,&
        nprotectLayers, dprotectLayerTolerance, rerror)

  end subroutine hydro_errestRecovery

  !*****************************************************************************

!<subroutine>

  subroutine hydro_calcProtectionLayer(rtriangulation, nlayers, dtreshold, rvector)

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
    call storage_new('hydro_calcProtectionLayer',' BisactiveElement',&
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
    
  end subroutine hydro_calcProtectionLayer

end module hydro_errorestimation

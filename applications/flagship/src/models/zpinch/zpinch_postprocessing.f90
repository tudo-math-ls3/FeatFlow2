!##############################################################################
!# ****************************************************************************
!# <Name> zpinch_postprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains all postprocessing routines which are required
!# to solve the time-dependent magnetohydrodynamic equations in the
!# one-, two- or three-dimensional domain $\Omega$.
!#
!# The following routines are available:
!#
!# 1.) zpinch_outputSolution
!#     -> Outputs the solution vector to file in UCD format
!#
!# 2.) zpinch_outputStatistics
!#      -> Outputs the application statitics
!#
!# 3.) zpinch_projectSolution
!#      -> Performs conservative projection of the solution from
!#         a given FE-space to another FE-space
!#
!# </purpose>
!##############################################################################

module zpinch_postprocessing

#include "flagship.h"
#include "models/hydro/hydro.h"

!$ use omp_lib
  use basicgeometry
  use collection
  use flagship_basic
  use fsystem
  use genoutput
  use hydro_basic
  use hydro_basic1d
  use hydro_basic2d
  use hydro_basic3d
  use hydro_postprocessing
  use linearsystemblock
  use linearsystemscalar
  use paramlist
  use problem
  use statistics
  use transport_postprocessing
  use ucd

  implicit none

  private

  public :: zpinch_outputSolution
  public :: zpinch_outputStatistics
  public :: zpinch_projectSolution
  
contains

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_outputSolution(rparlist, ssectionName,&
      rproblemLevel, Rsolution, dtime)

!<description>
    ! This subroutine exports the solution vector to file in UCD format
!</description>

!<input>
    ! parameter list
    type(t_parlist), intent(in) :: rparlist

    ! section name in parameter list
    character(LEN=*), intent(in) :: ssectionName

    ! problem level structure
    type(t_problemLevel), intent(in) :: rproblemLevel

    ! array of solution vectors
    type(t_vectorBlock), dimension(:), intent(in), optional :: Rsolution

    ! OPTIONAL: simulation time
    real(DP), intent(in), optional :: dtime
!</input>
!</subroutine>

    ! section names
    character(LEN=SYS_STRLEN) :: soutputName
    character(LEN=SYS_STRLEN) :: ssectionNameHydro
    character(LEN=SYS_STRLEN) :: sucdsolution

    ! persistent variable
    integer, save :: ifilenumber = 1

    ! local variables
    type(t_ucdExport) :: rexport
    type(t_vectorScalar) :: rvector1, rvector2, rvector3
    real(DP), dimension(:), pointer :: p_Dsolution, p_Ddata1, p_Ddata2, p_Ddata3
    character(len=SYS_NAMELEN) :: cvariable
    integer :: iformatUCD, isystemFormat, isize, ndim, nvar, ivariable, nvariable


    ! Get global configuration from parameter list
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'output', soutputName)
    call parlst_getvalue_string(rparlist,&
        trim(soutputName), 'sucdsolution', sucdsolution)
    call parlst_getvalue_int(rparlist,&
        trim(soutputName), 'iformatucd', iformatUCD, 0)
    call parlst_getvalue_string(rparlist,&
        ssectionName, 'subapplication', ssectionNameHydro, isubstring=1)
    call parlst_getvalue_int(rparlist,&
        ssectionNameHydro, 'isystemformat', isystemformat)
    
    if (iformatUCD .eq. 0) then
      call output_line('No valid output format is specified!',&
          OU_CLASS_WARNING,OU_MODE_STD,'zpinch_outputSolution')
      return
    end if

    ! Initialise the UCD exporter
    call flagship_initUCDexport(rproblemLevel, sucdsolution,&
        iformatUCD, rexport, ifilenumber)

    ! Increase filenumber by one
    ifilenumber = ifilenumber+1

    ! Set simulation time
    if (present(dtime)) call ucd_setSimulationTime(rexport, dtime)

    ! Add solution vectors
    if (present(Rsolution)) then

      ! Add solution for hydrodynamic model
      call lsysbl_getbase_double(Rsolution(1), p_Dsolution)
      nvar  = hydro_getNVAR(rproblemLevel)
      isize = size(p_Dsolution)/nvar
      ndim  = rproblemLevel%rtriangulation%ndim
      
      ! Create auxiliary vectors
      select case(ndim)
      case (NDIM1D)
        call lsyssc_createVector(rvector1, isize, .false.)
        call lsyssc_getbase_double(rvector1, p_Ddata1)

      case (NDIM2D)
        call lsyssc_createVector(rvector1, isize, .false.)
        call lsyssc_createVector(rvector2, isize, .false.)
        call lsyssc_getbase_double(rvector1, p_Ddata1)
        call lsyssc_getbase_double(rvector2, p_Ddata2)

      case (NDIM3D)
        call lsyssc_createVector(rvector1, isize, .false.)
        call lsyssc_createVector(rvector2, isize, .false.)
        call lsyssc_createVector(rvector3, isize, .false.)
        call lsyssc_getbase_double(rvector1, p_Ddata1)
        call lsyssc_getbase_double(rvector2, p_Ddata2)
        call lsyssc_getbase_double(rvector3, p_Ddata3)

      case default
        call output_line('Invalid number of spatial dimensions',&
            OU_CLASS_ERROR,OU_MODE_STD,'hydro_outputSolution')
        call sys_halt()
      end select

      ! Get number of variables to be written
      nvariable = max(1,&
          parlst_querysubstrings(rparlist,&
          trim(soutputName), 'sucdvariable'))

      select case(isystemFormat)
      case(SYSTEM_INTERLEAVEFORMAT)
        
        ! Loop over all variables
        do ivariable = 1, nvariable
          
          ! Get variable name
          call parlst_getvalue_string(rparlist, trim(soutputName),&
              'sucdvariable', cvariable, isubstring=ivariable)
          
          if (trim(cvariable) .eq. 'velocity') then
            
            ! Special treatment of velocity vector
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarInterleaveFormat1d(rvector1%NEQ, NVAR1D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)

            case (NDIM2D)
              call hydro_getVarInterleaveFormat2d(rvector1%NEQ, NVAR2D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call hydro_getVarInterleaveFormat2d(rvector2%NEQ, NVAR2D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call hydro_getVarInterleaveFormat3d(rvector1%NEQ, NVAR3D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call hydro_getVarInterleaveFormat3d(rvector2%NEQ, NVAR3D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call hydro_getVarInterleaveFormat3d(rvector3%NEQ, NVAR3D,&
                  'velocity_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select

          elseif (trim(cvariable) .eq. 'momentum') then
            
            ! Special treatment of momentum vector
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarInterleaveFormat1d(rvector1%NEQ, NVAR1D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'momentum', p_Ddata1)

            case (NDIM2D)
              call hydro_getVarInterleaveFormat2d(rvector1%NEQ, NVAR2D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call hydro_getVarInterleaveFormat2d(rvector2%NEQ, NVAR2D,&
                  'momentum_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'momentum',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call hydro_getVarInterleaveFormat3d(rvector1%NEQ, NVAR3D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call hydro_getVarInterleaveFormat3d(rvector2%NEQ, NVAR3D,&
                  'momentum_y', p_Dsolution, p_Ddata2)
              call hydro_getVarInterleaveFormat3d(rvector3%NEQ, NVAR3D,&
                  'momentum_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'momentum',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select

          elseif(trim(cvariable) .eq. 'advect') then

            ! Special treatment for tracer quantity
            call lsysbl_getbase_double(Rsolution(2), p_Ddata1)
            call ucd_addVariableVertexBased(rexport, 'advect',&
                UCD_VAR_STANDARD, p_Ddata1)

          else

            ! Standard treatment for scalar quantity
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarInterleaveFormat1d(rvector1%NEQ, nvar,&
                  cvariable, p_Dsolution, p_Ddata1)
            case (NDIM2D)
              call hydro_getVarInterleaveFormat2d(rvector1%NEQ, nvar,&
                  cvariable, p_Dsolution, p_Ddata1)
            case (NDIM3D)
              call hydro_getVarInterleaveFormat3d(rvector1%NEQ, nvar,&
                  cvariable, p_Dsolution, p_Ddata1)
            end select
            call ucd_addVariableVertexBased(rexport, cvariable,&
                UCD_VAR_STANDARD, p_Ddata1)
            
          end if
        end do
        
      case (SYSTEM_BLOCKFORMAT)

        ! Loop over all variables
        do ivariable = 1, nvariable
          
          ! Get variable name
          call parlst_getvalue_string(rparlist, trim(soutputName),&
              'sucdvariable', cvariable, isubstring=ivariable)
          
          if (trim(cvariable) .eq. 'velocity') then

            ! Special treatment of velocity vector
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarBlockFormat1d(rvector1%NEQ, NVAR1D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'velocity', p_Ddata1)

            case (NDIM2D)
              call hydro_getVarBlockFormat2d(rvector1%NEQ, NVAR2D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call hydro_getVarBlockFormat2d(rvector2%NEQ, NVAR2D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call hydro_getVarBlockFormat3d(rvector1%NEQ, NVAR3D,&
                  'velocity_x', p_Dsolution, p_Ddata1)
              call hydro_getVarBlockFormat3d(rvector2%NEQ, NVAR3D,&
                  'velocity_y', p_Dsolution, p_Ddata2)
              call hydro_getVarBlockFormat3d(rvector3%NEQ, NVAR3D,&
                  'velocity_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'velocity',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select
        
          elseif (trim(cvariable) .eq. 'momentum') then

            ! Special treatment of momentum vector
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarBlockFormat1d(rvector1%NEQ, NVAR1D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call ucd_addVarVertBasedVec(rexport, 'momentum', p_Ddata1)

            case (NDIM2D)
              call hydro_getVarBlockFormat2d(rvector1%NEQ, NVAR2D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call hydro_getVarBlockFormat2d(rvector2%NEQ, NVAR2D,&
                  'momentum_y', p_Dsolution, p_Ddata2)
              call ucd_addVarVertBasedVec(rexport, 'momentum',&
                  p_Ddata1, p_Ddata2)

            case (NDIM3D)
              call hydro_getVarBlockFormat3d(rvector1%NEQ, NVAR3D,&
                  'momentum_x', p_Dsolution, p_Ddata1)
              call hydro_getVarBlockFormat3d(rvector2%NEQ, NVAR3D,&
                  'momentum_y', p_Dsolution, p_Ddata2)
              call hydro_getVarBlockFormat3d(rvector3%NEQ, NVAR3D,&
                  'momentum_z', p_Dsolution, p_Ddata3)
              call ucd_addVarVertBasedVec(rexport, 'momentum',&
                  p_Ddata1, p_Ddata2, p_Ddata3)
            end select

          elseif(trim(cvariable) .eq. 'advect') then
            
            ! Special treatment for tracer quantity
            call lsysbl_getbase_double(Rsolution(2), p_Ddata1)
            call ucd_addVariableVertexBased(rexport, 'advect',&
                UCD_VAR_STANDARD, p_Ddata1)

          else
            
            ! Standard treatment for scalar quantity
            select case(ndim)
            case (NDIM1D)
              call hydro_getVarBlockFormat1d(rvector1%NEQ, nvar,&
                  cvariable, p_Dsolution, p_Ddata1)
            case (NDIM2D)
              call hydro_getVarBlockFormat2d(rvector1%NEQ, nvar,&
                  cvariable, p_Dsolution, p_Ddata1)
            case (NDIM3D)
              call hydro_getVarBlockFormat3d(rvector1%NEQ, nvar,&
                  cvariable, p_Dsolution, p_Ddata1)
            end select
            call ucd_addVariableVertexBased(rexport, cvariable,&
                UCD_VAR_STANDARD, p_Ddata1)
            
          end if
        end do

      case default
        call output_line('Invalid system format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'zpinch_outputSolution')
        call sys_halt()
      end select

      ! Release temporal memory
      call lsyssc_releaseVector(rvector1)
      call lsyssc_releaseVector(rvector2)
      call lsyssc_releaseVector(rvector3)

    end if

    ! Write UCD file
    call ucd_write(rexport)
    call ucd_release(rexport)

  end subroutine zpinch_outputSolution

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_outputStatistics(rtimerTotal, ssectionName, rcollection)

!<description>
    ! This subroutine output application statistics
!</description>

!<input>
    ! timer for total time measurement
    type(t_timer), intent(in) :: rtimerTotal

    ! section name in collection structure
    character(LEN=*), intent(in) :: ssectionName
!</input>

!<inputoutput>
    ! collection structure
    type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    ! local variable
    type(t_timer), pointer :: p_rtimerSolution
    type(t_timer), pointer :: p_rtimerAdaptation
    type(t_timer), pointer :: p_rtimerErrorEstimation
    type(t_timer), pointer :: p_rtimerTriangulation
    type(t_timer), pointer :: p_rtimerAssemblyCoeff
    type(t_timer), pointer :: p_rtimerAssemblyMatrix
    type(t_timer), pointer :: p_rtimerAssemblyVector
    type(t_timer), pointer :: p_rtimerPrePostprocess
    real(DP) :: dfraction


    ! Get timer objects from collection
    p_rtimerSolution => collct_getvalue_timer(rcollection,&
        'rtimerSolution', ssectionName=ssectionName)
    p_rtimerAdaptation => collct_getvalue_timer(rcollection,&
        'rtimerAdaptation', ssectionName=ssectionName)
    p_rtimerErrorEstimation => collct_getvalue_timer(rcollection,&
        'rtimerErrorEstimation', ssectionName=ssectionName)
    p_rtimerTriangulation => collct_getvalue_timer(rcollection,&
        'rtimerTriangulation', ssectionName=ssectionName)
    p_rtimerAssemblyCoeff => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyCoeff', ssectionName=ssectionName)
    p_rtimerAssemblyMatrix => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyMatrix', ssectionName=ssectionName)
    p_rtimerAssemblyVector => collct_getvalue_timer(rcollection,&
        'rtimerAssemblyVector', ssectionName=ssectionName)
    p_rtimerPrePostprocess => collct_getvalue_timer(rcollection,&
        'rtimerPrePostprocess', ssectionName=ssectionName)

    ! Output statistics
    call output_lbrk()
    call output_line('Time measurement: '//trim(adjustl(ssectionName)))
    call output_line('-----------------')

    call stat_subTimers(p_rtimerAssemblyMatrix, p_rtimerSolution)
    call stat_subTimers(p_rtimerAssemblyVector, p_rtimerSolution)

    dfraction = 100.0_DP/rtimerTotal%delapsedReal

    call output_line('                                Real time (seconds) CPU time (seconds)       Frac')
    call output_line('Time for computing solution   : '//&
                     trim(sys_sdEL(p_rtimerSolution%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerSolution%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerSolution%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerSolution%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerSolution%delapsedCPU/&
                     (p_rtimerSolution%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for mesh adaptivity      : '//&
                     trim(sys_sdEL(p_rtimerAdaptation%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAdaptation%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerAdaptation%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAdaptation%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerAdaptation%delapsedCPU/&
                     (p_rtimerAdaptation%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for error estimation     : '//&
                     trim(sys_sdEL(p_rtimerErrorEstimation%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerErrorEstimation%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerErrorEstimation%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerErrorEstimation%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerErrorEstimation%delapsedCPU/&
                     (p_rtimerErrorEstimation%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for triangulation        : '//&
                     trim(sys_sdEL(p_rtimerTriangulation%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerTriangulation%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerTriangulation%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerTriangulation%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerTriangulation%delapsedCPU/&
                     (p_rtimerTriangulation%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for coefficient assembly : '//&
                     trim(sys_sdEL(p_rtimerAssemblyCoeff%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAssemblyCoeff%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerAssemblyCoeff%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAssemblyCoeff%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerAssemblyCoeff%delapsedCPU/&
                     (p_rtimerAssemblyCoeff%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for matrix assembly      : '//&
                     trim(sys_sdEL(p_rtimerAssemblyMatrix%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAssemblyMatrix%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerAssemblyMatrix%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAssemblyMatrix%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerAssemblyMatrix%delapsedCPU/&
                     (p_rtimerAssemblyMatrix%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for vector assembly      : '//&
                     trim(sys_sdEL(p_rtimerAssemblyVector%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAssemblyVector%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerAssemblyVector%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerAssemblyVector%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerAssemblyVector%delapsedCPU/&
                     (p_rtimerAssemblyVector%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_line('Time for pre-/post-processing : '//&
                     trim(sys_sdEL(p_rtimerPrePostprocess%delapsedReal,3))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerPrePostprocess%delapsedReal,2),5)//' %  '//&
                     trim(sys_sdEL(p_rtimerPrePostprocess%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*p_rtimerPrePostprocess%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(p_rtimerPrePostprocess%delapsedCPU/&
                     (p_rtimerPrePostprocess%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_lbrk()
    call output_line('Time for total simulation     : '//&
                     trim(sys_sdEL(rtimerTotal%delapsedReal,3))//'           '//&
                     trim(sys_sdEL(rtimerTotal%delapsedCPU,5))//'  '//&
                     sys_trimr(sys_sd(dfraction*rtimerTotal%delapsedCPU,2),7)//' %  '//&
                     sys_trimr(sys_sd(rtimerTotal%delapsedCPU/&
                     (rtimerTotal%delapsedReal+SYS_EPSREAL_DP),2),5))
    call output_lbrk()

  end subroutine zpinch_outputStatistics

  !*****************************************************************************

!<subroutine>

  subroutine zpinch_projectSolution(RsourceVector, RdestVector)
    
!<description>
    ! This subroutine performs conservative projection of the given solution
    ! stored in rsourceVector to another FE-space and stores the result in
    ! rdestVector. An FCT algorithm is used ensure monotonicity preservation.
!</description>

!<input>
    ! Source vector
    type(t_vectorBlock), dimension(:), intent(inout) :: rsourceVector
!</input>

!<inputoutput>
    ! Destination vector
    type(t_vectorBlock), dimension(:), intent(inout) :: rdestVector
!</inputoutput>
!</subroutine>

    ! Project solution of hydrodynamic model
    call hydro_projectSolution(RsourceVector(1), RdestVector(1))

    ! Project solution of scalar transport model
    call transp_projectSolution(RsourceVector(2), RdestVector(2))

  end subroutine zpinch_projectSolution

end module zpinch_postprocessing

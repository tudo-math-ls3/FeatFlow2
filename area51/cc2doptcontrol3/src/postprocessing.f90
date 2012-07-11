!##############################################################################
!# ****************************************************************************
!# <name> postprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains various routines for the postprocessing during a
!# space-time optimisation.
!#
!# The following routines can be found here:
!#
!# 1.) optcpp_initpostprocessing
!#     -> Initialise the postprocessing
!#
!# 2.) optcpp_donepostprocessing
!#     -> Clean up the postprocessing
!#
!# 3.) optcpp_postprocessSingleSol
!#     -> Postprocessing of a solution at a specific timestep
!#
!# 4.) optcpp_postprocessSpaceTimeVec
!#     -> Postprocessing of a space-time vector
!#
!# 5.) optcpp_postprocSpaceVisOutput
!#     -> Write out visualisation files of all solutions in a space-time vector
!#
!# 6.) optcpp_calculateControl
!#     -> Calculates the control from a solution vector.
!# </purpose>
!##############################################################################

module postprocessing

  use fsystem
  use genoutput
  use io
  use paramlist
  use mprimitives
  
  use basicgeometry
  use spatialdiscretisation
  use timediscretisation
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  
  use scalarpde
  use linearformevaluation
  use bilinearformevaluation
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use collection
  use feevaluation
  use pprocerror
  use ucd
  use pprocnavierstokes
  use vectorio
  
  use spacetimevectors
  use analyticsolution
  
  use constantsdiscretisation
  use structuresdiscretisation
  use structuresoptcontrol
  use structuresgeneral
  use structuresoperatorasm
  use structuresoptflow
  use structurespostproc
  
  use assemblytemplates
  use spacematvecassembly
  
  use kktsystemspaces
  use kktsystem

  use optcanalysis  
  
  use user_callback

  implicit none
  
  private
  
!!<constants>
!
!!<constantblock description="Type identifiers for creating a space-time vector.">
!
!  ! A zero space-time vector.
!  integer, parameter :: CCSTV_ZERO = 0
!
!  ! The space-time vector is created as copy of a single vector.
!  integer, parameter :: CCSTV_STATIONARY = 1
!
!  ! The space-time vector is created by reading in files.
!  integer, parameter :: CCSTV_READFILES = 2
!
!  ! The space-time vector is created as forward simulation starting from an
!  ! initial solution.
!  integer, parameter :: CCSTV_FORWARDSIMULATION = 3
!
!!</constantblock>
!
!!</constants>

  public :: optcpp_initpostprocessing
  public :: optcpp_donepostprocessing
  public :: optcpp_postprocessing
  public :: optcpp_postprocessSubstep
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_initpostprocessing (&
      rpostproc,rparlist,ssection,rphysics,rsettingsSpaceDiscr)
  
!<description>
  ! Initialises the postprocessing.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters can be found.
  character(len=*), intent(in) :: ssection

  ! Physics of the problem
  type(t_settings_physics), intent(in), target :: rphysics
  
  ! Spatial discretisation settings
  type(t_settings_spacediscr), intent(in), target :: rsettingsSpaceDiscr
!</input>
  
!<inputoutput>
  ! Parameters about the postprocessing.
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>
  
!</subroutine>

    ! Fetch data
    call struc_initPostprocParams (rparlist,ssection,rphysics,rpostproc)

    rpostproc%p_rphysics => rphysics
    rpostproc%p_rsettingsDiscr => rsettingsSpaceDiscr
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_donepostprocessing (rpostproc)
  
!<description>
  ! Clean up the postprocessing structure
!</description>
  
!<inputoutput>
  ! Parameters about the postprocessing.
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>
  
!</subroutine>

    ! Clean up.
    call struc_donePostprocParams (rpostproc)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_unshiftSubvectors (rpostproc,iorder,&
      dtimePrimal1,dtimeDual1,rvector1,&
      dtimePrimal2,dtimeDual2,rvector2,&
      dtimePrimal3,dtimeDual3,rvector3,&
      dtime,rvector)
  
!<description>
  ! 'Unshifts' time shifts in the components of a vector.
  ! Interpolates quadratically a solution from rvectorPrev, rvectorCurr
  ! and rvectorNext into rvectorTarget such that it represents
  ! the solution at time dtimeTarget.
  !
  ! note: there must be dtime1 <= dtime2 <= dtime3 !
!</description>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!<input>
  ! Order of the interpolation. =1: linear, =2: quadratic.
  integer, intent(in) :: iorder

  ! 1st simulation time in the primal and dual equation.
  ! and the vector.
  real(dp), intent(in) :: dtimePrimal1
  real(dp), intent(in) :: dtimeDual1
  type(t_vectorBlock), intent(inout) :: rvector1

  ! 2nd simulation time in the primal and dual equation.
  ! and the vector.
  real(dp), intent(in) :: dtimePrimal2
  real(dp), intent(in) :: dtimeDual2
  type(t_vectorBlock), intent(inout) :: rvector2

  ! 3rd simulation time in the primal and dual equation.
  ! and the vector.
  real(dp), intent(in) :: dtimePrimal3
  real(dp), intent(in) :: dtimeDual3
  type(t_vectorBlock), intent(inout) :: rvector3
!</input>

!<inputoutput>
  ! Destination simulation time in the primal and dual equation.
  ! and the vector.
  real(dp), intent(in) :: dtime
  type(t_vectorBlock), intent(inout) :: rvector
!</inputoutput>

    ! local variables
    real(DP) :: tmidPrimal, tmidDual, tprimal, tdual

    ! General quadratic interpolation:
    !
    !   p(t1) = p1
    !   p(t2) = p2
    !   p(t3) = p3
    !
    ! Substitute: tmid := sigma(t2) := (t1+t3-2*t2)/(t1-t3)
    ! => t1 -> -1,  t3 -> 1
    ! then:
    !
    !   p(-1)   = p1
    !   p(tmid) = p2
    !   p(1)    = p3
    !
    ! The linear interpolation polynom with
    !
    !   p(-1)   = p1
    !   p(1)    = p3
    !
    ! is
    !
    !   p(t) = (-(1/2)*t+1/2)*p1 + ((1/2)*t+1/2)*p3
    !
    ! and through 3 disjuct points
    !
    !  p(t) =  (1/2)*(t2*t^2-t2^2*t+t2^2-t^2+t-t2)*p1/(t2^2-1)
    !         +(1/2)*(2*t^2-2)*p2/(t2^2-1)
    !         +(1/2)*(-t2*t^2+t2^2*t+t2^2-t^2-t+t2)*p3/(t2^2-1)
    
    select case (rpostproc%p_rphysics%cequation)
    case (0,1)
      ! Map the time midpoint(s)
      tmidPrimal = (dtimePrimal1+dtimePrimal3-2.0_DP*dtimePrimal2)/(dtimePrimal1-dtimePrimal3)
      tmidDual = (dtimeDual1+dtimeDual3-2.0_DP*dtimeDual2)/(dtimeDual1-dtimeDual3)

      tprimal = (dtimePrimal1+dtimePrimal3-2.0_DP*dtime)/(dtimePrimal1-dtimePrimal3)
      tdual = (dtimeDual1+dtimeDual3-2.0_DP*dtime)/(dtimeDual1-dtimeDual3)
    
      ! Interpolate the subvectors.
      call interpolateScalar (iorder,&
          rvector1%RvectorBlock(1),rvector2%RvectorBlock(1),rvector3%RvectorBlock(1),&
          tmidPrimal,tprimal,rvector%RvectorBlock(1))
      call interpolateScalar (iorder,&
          rvector1%RvectorBlock(2),rvector2%RvectorBlock(2),rvector3%RvectorBlock(2),&
          tmidPrimal,tprimal,rvector%RvectorBlock(2))
      call interpolateScalar (iorder,&
          rvector1%RvectorBlock(3),rvector2%RvectorBlock(3),rvector3%RvectorBlock(3),&
          tmidDual,tdual,rvector%RvectorBlock(3))
          
      call interpolateScalar (iorder,&
          rvector1%RvectorBlock(4),rvector2%RvectorBlock(4),rvector3%RvectorBlock(4),&
          tmidDual,tdual,rvector%RvectorBlock(4))
      call interpolateScalar (iorder,&
          rvector1%RvectorBlock(5),rvector2%RvectorBlock(5),rvector3%RvectorBlock(5),&
          tmidDual,tdual,rvector%RvectorBlock(5))
      call interpolateScalar (iorder,&
          rvector1%RvectorBlock(6),rvector2%RvectorBlock(6),rvector3%RvectorBlock(6),&
          tmidPrimal,tprimal,rvector%RvectorBlock(6))
    end select
    
  contains
  
    ! linear interpolation routine for a scalar vector.
    ! rvec1 is at time -1, rvec2 at time tmid in [-1,1], dtime3 at time 3.
    ! The result at time t in [-1,1] is written to rvec.
    subroutine interpolateScalar (iorder,rvec1,rvec2,rvec3,tmid,t,rvec)
      type(t_vectorScalar), intent(in) :: rvec1, rvec2, rvec3
      type(t_vectorScalar), intent(inout) :: rvec
      real(DP), intent(in) :: tmid,t
      integer :: iorder
      
      select case (iorder)
      case (1)
        ! Left or right from the time midpoint?
        if (t .le. tmid) then
          ! Linear interpolation at the left interval
          call lsyssc_copyVector (rvec2,rvec)
          call lsyssc_vectorLinearComb (rvec1,rvec,lincoeff1(-1.0_DP,tmid,t),lincoeff2(-1.0_DP,tmid,t))
        else
          ! Linear interpolation at the right interval
          call lsyssc_copyVector (rvec3,rvec)
          call lsyssc_vectorLinearComb (rvec2,rvec,lincoeff1(tmid,1.0_DP,t),lincoeff2(tmid,1.0_DP,t))
        end if
      case (2)
        call output_line ("Not yet implemented!")
        call sys_halt()
      end select
    end subroutine

    ! Linear interpolation formula
  
    real(DP) function lincoeff1 (t1,t2,t)
      real(DP), intent(in) :: t1,t2,t
      if (t1 .eq. t2) then
        lincoeff1 = 1.0_DP
      else
        lincoeff1 = (t-t2)/(t1-t2)
      end if
    end function

    real(DP) function lincoeff2 (t1,t2,t)
      real(DP), intent(in) :: t1,t2,t
      lincoeff2 = 1.0_DP-lincoeff1 (t1,t2,t)
    end function

    ! Coefficient functions for the linear and quadratic interpolation formula.

    real(DP) function quadcoeff1 (t2,t)
      real(DP), intent(in) :: t2,t
      quadcoeff1 = 0.5_DP*(t2*t**2-t2**2*t+t2**2-t**2+t-t2)/(t2**2-1.0_DP)
    end function

    real(DP) function quadcoeff2 (t2,t)
      real(DP), intent(in) :: t2,t
      quadcoeff2 = 0.5_DP*(2.0_DP*t**2-2.0_DP)/(t2**2-1.0_DP)
    end function

    real(DP) function quadcoeff3 (t2,t)
      real(DP), intent(in) :: t2,t
      quadcoeff3 = 0.5_DP*(-t2*t**2+t2**2*t+t2**2-t**2-t+t2)/(t2**2-1.0_DP)
    end function

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_singleVisualisation (rpostproc,rvector,cspace,&
      rtriangulation,rphysics,roptcontrol,ifileid,dtime,itag)
  
!<description>
  ! Writes a visualisation file for the given solution rvector.
!</description>

!<input>
  ! Postprocessing settings
  type(t_optcPostprocessing), intent(in) :: rpostproc
  
  ! Time associated with this solution
  real(DP), intent(in) :: dtime
  
  ! File id to be appended to the filename.
  integer, intent(in) :: ifileid
  
  ! Vector data to write into the file.
  type(t_vectorBlock), intent(in) :: rvector
  
  ! Identifier that defines the type of rvector.
  integer, intent(in) :: cspace
  
  ! Underlying triangulation
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics
  
  ! Underlying optimal control parameters
  type(t_settings_optcontrol), intent(in) :: roptcontrol
  
  ! OPTIONAL: An integer tag which is included into filenames of
  ! output files. If not specified, the tag is not included.
  integer, intent(in), optional :: itag
!</input>

!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sfile,ssuffix
    type(t_ucdExport) :: rexport
    integer :: icomp
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2,p_Ddata3
    integer :: nvt,nmt,nel
    
    ! Determine filename header.
    ! Current space? 
    select case (cspace)
    case (CCSPACE_PRIMAL)
      sfile = trim(rpostproc%sfilenameUCD)//"_primal"
    case (CCSPACE_DUAL)
      sfile = trim(rpostproc%sfilenameUCD)//"_dual"
    case (CCSPACE_CONTROL)
      sfile = trim(rpostproc%sfilenameUCD)//"_control"
    end select
    
    ! Create a filename for the visualisation output
    if (present(itag)) then
      ssuffix = "."//trim(sys_siL(itag,10))//"."//trim(sys_si0L(ifileid,5))
    else
      ssuffix = "."//trim(sys_si0L(ifileid,5))
    end if

    ! Open the visualisation file or cancel
    select case (rpostproc%ioutputUCD)
    case (0)
      return
      
    case (1)
      sfile = trim(sfile)//".gmv"//trim(ssuffix)
      call ucd_startGMV (rexport,&
          UCD_FLAG_STANDARD+UCD_FLAG_IGNOREDEADNODES,rtriangulation,&
          sfile)

    case (2)
      sfile = trim(sfile)//".avs"//trim(ssuffix)
      call ucd_startAVS (rexport,&
          UCD_FLAG_STANDARD+UCD_FLAG_IGNOREDEADNODES,rtriangulation,&
          sfile)
          
    case (3)
      sfile = trim(sfile)//".vtk"//trim(ssuffix)
      call ucd_startVTK (rexport,&
          UCD_FLAG_STANDARD+UCD_FLAG_IGNOREDEADNODES,rtriangulation,&
          sfile)

    case (4)
      sfile = trim(sfile)//".vtk"//trim(ssuffix)
      call ucd_startVTK (rexport,&
          UCD_FLAG_STANDARD+UCD_FLAG_IGNOREDEADNODES+&
          UCD_FLAG_ONCEREFINED+UCD_FLAG_AUTOINTERPOLATE,rtriangulation,&
          sfile)
          
    case default
      call output_line ("Invalid UCD output type.", &
          OU_CLASS_ERROR,OU_MODE_STD,"optcpp_singleVisualisation")
      call sys_halt()
    end select

    call output_line ("Writing vis. file: "//trim(sfile))    
    
    ! Solution time
    call ucd_setSimulationTime(rexport,dtime)
    
    ! Underlying equation
    select case (rphysics%cequation)
    
    ! ********************************************
    ! Stokes / Navier-Stokes
    ! ********************************************
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)
    
      ! Current space?
      select case (cspace)
      
      ! -------------------------------------------
      ! Primal equation
      ! -------------------------------------------
      case (CCSPACE_PRIMAL)
      
        select case (rpostproc%p_rsettingsDiscr%ielementType)
        
        ! --------------------
        ! Q2/QP1.
        ! Manual output.
        ! --------------------
        case (4)
          ! Data arrays
          call lsyssc_getbase_double (Rvector%RvectorBlock(1),p_Ddata1)
          call lsyssc_getbase_double (Rvector%RvectorBlock(2),p_Ddata2)
          
          ! Size of the arrays
          nvt = Rvector%p_rblockDiscr%p_rtriangulation%nvt
          nmt = Rvector%p_rblockDiscr%p_rtriangulation%nmt
          nel = Rvector%p_rblockDiscr%p_rtriangulation%nel
          
          ! Write the velocity field
          call ucd_addVarVertBasedVec(rexport, "velocity_p", UCD_VAR_STANDARD,&
              p_Ddata1(1:nvt), p_Ddata2(1:nvt), &
              DdataMid_X=p_Ddata1(nvt+1:nvt+nmt), &
              DdataMid_Y=p_Ddata2(nvt+1:nvt+nmt), &
              DdataElem_X=p_Ddata1(nvt+nmt+1:nvt+nmt+nel), &
              DdataElem_Y=p_Ddata2(nvt+nmt+1:nvt+nmt+nel))

          ! Write the pressure data
          call ucd_addVectorByVertex (rexport, &
              "pressure_p", UCD_VAR_STANDARD, Rvector%RvectorBlock(3))
        
        ! --------------------
        ! Standard output
        ! --------------------
        case default
      
          ! Write the velocity field
          call ucd_addVectorFieldByVertex (rexport, &
              "velocity_p", UCD_VAR_STANDARD,Rvector%RvectorBlock(1:2))
              
          ! Write the pressure data
          call ucd_addVectorByVertex (rexport, &
              "pressure_p", UCD_VAR_STANDARD, Rvector%RvectorBlock(3))
              
        end select
      
      ! -------------------------------------------
      ! Dual equation
      ! -------------------------------------------
      case (CCSPACE_DUAL)

        select case (rpostproc%p_rsettingsDiscr%ielementType)
        
        ! --------------------
        ! Q2/QP1.
        ! Manual output.
        ! --------------------
        case (4)
          ! Data arrays
          call lsyssc_getbase_double (Rvector%RvectorBlock(1),p_Ddata1)
          call lsyssc_getbase_double (Rvector%RvectorBlock(2),p_Ddata2)
          
          ! Size of the arrays
          nvt = Rvector%p_rblockDiscr%p_rtriangulation%nvt
          nmt = Rvector%p_rblockDiscr%p_rtriangulation%nmt
          nel = Rvector%p_rblockDiscr%p_rtriangulation%nel
          
          ! Write the velocity field
          call ucd_addVarVertBasedVec(rexport, "velocity_d", UCD_VAR_STANDARD,&
              p_Ddata1(1:nvt), p_Ddata2(1:nvt), &
              DdataMid_X=p_Ddata1(nvt+1:nvt+nmt), &
              DdataMid_Y=p_Ddata2(nvt+1:nvt+nmt), &
              DdataElem_X=p_Ddata1(nvt+nmt+1:nvt+nmt+nel), &
              DdataElem_Y=p_Ddata2(nvt+nmt+1:nvt+nmt+nel))

          ! Write the pressure data
          call ucd_addVectorByVertex (rexport, &
              "pressure_d", UCD_VAR_STANDARD, Rvector%RvectorBlock(3))

        ! --------------------
        ! Standard output
        ! --------------------
        case default
          ! Write the velocity field
          call ucd_addVectorFieldByVertex (rexport, &
              "velocity_d", UCD_VAR_STANDARD,Rvector%RvectorBlock(1:2))
              
          ! Write the pressure data
          call ucd_addVectorByVertex (rexport, &
              "pressure_d", UCD_VAR_STANDARD, Rvector%RvectorBlock(3))
        
        end select
      
      ! -------------------------------------------
      ! Control equation
      ! -------------------------------------------
      case (CCSPACE_CONTROL)
      
        select case (rpostproc%p_rsettingsDiscr%ielementType)
        
        ! --------------------
        ! Standard output
        ! --------------------
        case default
          ! The shape of the control variable depends on what is controlled...
          ! Initialise a component counter.
          icomp = 1
          
          ! -------------------------------------------------
          ! Distributed control
          ! -------------------------------------------------
          if (roptcontrol%dalphaC .ge. 0.0_DP) then
            call ucd_addVectorFieldByVertex (rexport, &
                "control", UCD_VAR_STANDARD,Rvector%RvectorBlock(icomp:icomp+1))
            icomp = icomp + 2
          end if
        end select
      
      end select

    ! ********************************************
    ! Heat equation
    ! ********************************************
    case (CCEQ_HEAT2D)
    
      ! Current space?
      select case (cspace)
      
      ! -------------------------------------------
      ! Primal equation
      ! -------------------------------------------
      case (CCSPACE_PRIMAL)
      
        ! Write the solution data
        call ucd_addVectorByVertex (rexport, &
            "solution_p", UCD_VAR_STANDARD, Rvector%RvectorBlock(1))
      
      ! -------------------------------------------
      ! Dual equation
      ! -------------------------------------------
      case (CCSPACE_DUAL)

        ! Write the solution data
        call ucd_addVectorByVertex (rexport, &
            "solution_d", UCD_VAR_STANDARD, Rvector%RvectorBlock(1))
      
      ! -------------------------------------------
      ! Control equation
      ! -------------------------------------------
      case (CCSPACE_CONTROL)
      
        ! The shape of the control variable depends on what is controlled...
        ! Initialise a component counter.
        icomp = 1
        
        ! -------------------------------------------------
        ! Distributed control
        ! -------------------------------------------------
        if (roptcontrol%dalphaC .ge. 0.0_DP) then
          call ucd_addVectorByVertex (rexport, &
              "control", UCD_VAR_STANDARD,Rvector%RvectorBlock(icomp))
          icomp = icomp + 1
        end if

      end select
    
    end select

    ! Write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_singleVectorDump (rpostproc,rvector,cspace,&
      ifileid,itag)
  
!<description>
  ! Dumps a vector into a file.
!</description>

!<input>
  ! Postprocessing settings
  type(t_optcPostprocessing), intent(in) :: rpostproc
  
  ! File id to be appended to the filename.
  integer, intent(in) :: ifileid
  
  ! Vector data to write into the file.
  type(t_vectorBlock), intent(in) :: rvector
  
  ! Identifier that defines the type of rvector.
  integer, intent(in) :: cspace
  
  ! OPTIONAL: An integer tag which is included into filenames of
  ! output files. If not specified, the tag is not included.
  integer, intent(in), optional :: itag
!</input>

!</subroutine>

    ! local variables
    character(len=sys_strlen) :: sfile
    character(len=20) :: sname
    integer :: cwrite
    
    ! Determine filename header.
    ! Current space? 
    select case (cspace)
    case (CCSPACE_PRIMAL)
      sfile = rpostproc%sfilenamePrimalSol
      sname = "primalsol"
      cwrite = rpostproc%cwritePrimalSol
    case (CCSPACE_DUAL)
      sfile = rpostproc%sfilenameDualSol
      sname = "dualsol"
      cwrite = rpostproc%cwriteDualSol
    case (CCSPACE_CONTROL)
      sfile = rpostproc%sfilenameControl
      sname = "control"
      cwrite = rpostproc%cwriteControl
    end select
    
    ! Create a filename for the visualisation output
    if (present(itag)) then
      sfile = trim(sfile)//"."//trim(sys_siL(itag,10))//"."//trim(sys_si0L(ifileid,5))
    else
      sfile = trim(sfile)//"."//trim(sys_si0L(ifileid,5))
    end if
    
    ! Write the vector
    select case (cwrite)
    case (0)
      return
      
    case (1)
      call output_line ("Writing dump file: "//trim(sfile))

      call vecio_writeBlockVectorHR (rvector, sname, .true.,&
          0, sfile)

    case (2)
      call output_line ("Writing dump file: "//trim(sfile))

      call vecio_writeBlockVectorHR (rvector, sname, .true.,&
          0, sfile, "(E20.10)")
          
    case default
      call output_line ("Invalid output type.", &
          OU_CLASS_ERROR,OU_MODE_STD,"optcpp_singleVectorDump")
      call sys_halt()
    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_appendLineToDatafile (sline,sheadline,sfilename,bfirst,&
      itag,iunit,bnewline)
  
!<description>
  ! Dumps a vector into a file.
!</description>

!<input>
  ! Filename
  character(len=*), intent(in) :: sfilename
  
  ! Headlime written to the beginning of the file
  character(len=*), intent(in) :: sheadline
  
  ! The line to append to the file
  character(len=*), intent(in) :: sline
  
  ! Identifies whether or not this is the first line to be written to the file.
  logical, intent(in) :: bfirst
  
  ! OPTIONAL: An integer tag which is included into filenames of
  ! output files. If not specified, the tag is not included.
  integer, intent(in), optional :: itag
  
  ! OPTIONEL: Flag specifying whether or not a line break is appended
  ! to the line. Default value is TRUE.
  logical, intent(in), optional :: bnewline
!</input>

!<inputoutput>
  ! OPTIONAL: An identifier for the file.
  ! If 0 is passed here, a new file is opened and iunit is replaced by the
  ! file handle. Then, the caller is responsible for closing the file.
  ! If a value != 0 is passed here, the data is directly written
  ! to this file without opening it.
  integer, intent(inout), optional :: iunit
!</inputoutput>

!</subroutine>

    ! local variables
    character(len=sys_strlen) :: sfile
    integer :: ifilehandle,cflag
    logical :: bfileExists

    ! If this is the "first" dataset, overwrite the file.
    ! Otherwise append.
    cflag = SYS_APPEND
    if (bfirst) cflag = SYS_REPLACE

    ! Create a filename
    sfile = sfilename
    if (present(itag)) then
      sfile = trim(sfile)//"."//trim(sys_siL(itag,10))
    end if
    
    ifilehandle = 0
    if (present(iunit)) ifilehandle = iunit
    
    if (ifilehandle .eq. 0) then
      ! Open the file, probably write the headline, write the data
      call io_openFileForWriting(sfile, ifilehandle, &
          cflag, bfileExists,.true.)
    else 
      bfileExists = .true.
    end if
        
    if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
      ! Write a headline
      write (ifilehandle,'(A)') trim(sheadline)
    end if
    
    if (present(bnewline)) then
      if (.not. bnewline) then
        write (ifilehandle,"(A)",ADVANCE="NO") trim(sline)
      else
        write (ifilehandle,"(A)") trim(sline)
      end if
    else
      write (ifilehandle,"(A)") trim(sline)
    end if
    
    if (present(iunit)) then
      ! Return the file handle
      iunit = ifilehandle
    else
      ! Close, finish.
      close (ifilehandle)
    end if
    
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine optcpp_getComponentName (sname,rphysics,cspace,icomponent,cderiv)
  
!<description>
  ! Determins the name of a component of a vector.
!</description>

!<input>
  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics

  ! Underlying space. A CCSPACE_xxxx constant.
  integer, intent(in) :: cspace
  
  ! Referred component
  integer, intent(in) :: icomponent
  
  ! Derivative quantifier.
  ! One of the DER_DERIVxx_xxxx constants
  integer, intent(in) :: cderiv
!</input>

!<output>
  ! Name of the component.
  ! ="", if the component is invalid.
  character(len=*), intent(out) :: sname
!</output>

!</subroutine>

    ! local variables
    character(LEN=10), dimension(3,3), parameter :: SfctnamesStokesPrimal = reshape (&
      (/ "y1        ","y1_x      ","y1_y      " , &
         "y2        ","y2_x      ","y2_y      " , &
         "p         ","p_x       ","p_y       "  /) ,&
       (/ 3,3 /) )

    character(LEN=10), dimension(3,3), parameter :: SfctnamesStokesDual = reshape (&
      (/ "lambda1   ","lambda1_x ","lambda1_y " , &
         "lambda2   ","lambda2_x ","lambda2_y " , &
         "xi        ","xi_x      ","xi_y      "  /) ,&
       (/ 3,3 /) )

    character(LEN=10), dimension(3,1), parameter :: SfctnamesHeatPrimal = reshape (&
      (/ "y        ","y_x      ","y_y      " /) ,&
       (/ 3,1 /) )

    character(LEN=10), dimension(3,1), parameter :: SfctnamesHeatDual = reshape (&
      (/ "u        ","u_x      ","u_y      " /) ,&
       (/ 3,1 /) )


    sname = ""

    ! Which equation do we have?
    select case (rphysics%cequation)
    
    ! -------------------------------------------------------------
    ! Stokes/Navier Stokes.
    ! -------------------------------------------------------------
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)

      if ((icomponent .ge. 1) .and. (icomponent .le. 3)) then
      
        ! Underlying space
        select case (cspace)
        
        ! -------------------------------
        ! Primal space
        ! -------------------------------
        case (CCSPACE_PRIMAL)
        
          ! Derivative?
          select case (cderiv)
          case (DER_FUNC2D)
            sname = SfctnamesStokesPrimal(1,icomponent)

          case (DER_DERIV2D_X)
            sname = SfctnamesStokesPrimal(2,icomponent)

          case (DER_DERIV2D_Y)
            sname = SfctnamesStokesPrimal(3,icomponent)
          end select
        
        ! -------------------------------
        ! Dual space
        ! -------------------------------
        case (CCSPACE_DUAL)

          ! Derivative?
          select case (cderiv)
          case (DER_FUNC2D)
            sname = SfctnamesStokesDual(1,icomponent)

          case (DER_DERIV2D_X)
            sname = SfctnamesStokesDual(2,icomponent)

          case (DER_DERIV2D_Y)
            sname = SfctnamesStokesDual(3,icomponent)
          end select
        
        end select
      
      end if

    ! -------------------------------------------------------------
    ! Heat equation.
    ! -------------------------------------------------------------
    case (CCEQ_HEAT2D)

      if (icomponent .eq. 1) then

        ! Underlying space
        select case (cspace)
        
        ! -------------------------------
        ! Primal space
        ! -------------------------------
        case (CCSPACE_PRIMAL)
        
          ! Derivative?
          select case (cderiv)
          case (DER_FUNC2D)
            sname = SfctnamesHeatPrimal(1,icomponent)

          case (DER_DERIV2D_X)
            sname = SfctnamesHeatPrimal(2,icomponent)

          case (DER_DERIV2D_Y)
            sname = SfctnamesHeatPrimal(3,icomponent)
          end select
        
        ! -------------------------------
        ! Dual space
        ! -------------------------------
        case (CCSPACE_DUAL)

          ! Derivative?
          select case (cderiv)
          case (DER_FUNC2D)
            sname = SfctnamesHeatDual(1,icomponent)

          case (DER_DERIV2D_X)
            sname = SfctnamesHeatDual(2,icomponent)

          case (DER_DERIV2D_Y)
            sname = SfctnamesHeatDual(3,icomponent)
          end select
        
        end select
        
      end if

    end select    

  end subroutine  
  

!******************************************************************************

!<subroutine>

  subroutine optcpp_singleEvalPtValues (rpostproc,rvector,cspace,rphysics,&
      bfirst,isolutionid,dtime,itag)

!<description>
  ! Evaluates the solution in a number of points as configured in the DAT file.
!</description>
  
!<input>
  ! Solution vector
  type(t_vectorBlock), intent(in), target :: rvector
  
  ! Space corresponding to rvector
  integer, intent(in) :: cspace
  
  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics

  ! Id of the solution (Number of the timestep).
  integer, intent(in) :: isolutionid
  
  ! Must be set to TRUE for the first value.
  logical, intent(in) :: bfirst
  
  ! Solution time. =0 for stationary simulations.
  real(DP), intent(in) :: dtime

  ! OPTIONAL: An integer tag which is included into filenames of
  ! output files. If not specified, the tag is not included.
  integer, intent(in), optional :: itag
!</input>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: npoints,i,iunit
    integer :: cderiv
    real(dp), dimension(1) :: Dvalues
    character(LEN=SYS_STRLEN) :: sstr,stemp,sheadline

    if (.not. associated(rpostproc%p_DcoordsPointEvalPrimal)) return
    
    ! Get the number of points to evaluate
    npoints = ubound(rpostproc%p_DcoordsPointEvalPrimal,2)
        
    ! Prepare terminal output
    call output_lbrk()
    call output_line ('Point values')
    call output_line ('------------')

    ! File handle for writing to a data file
    iunit = 0

    ! Evaluate the function in these points.
    do i=1,npoints
    
      ! Which derivative to evaluate?
      cderiv = DER_FUNC2D
      select case (rpostproc%p_ItypePointEvalPrimal(2,i))
      case (0)
        cderiv = DER_FUNC2D
      case (1)
        cderiv = DER_DERIV2D_X
      case (2)
        cderiv = DER_DERIV2D_Y
      end select
      
      ! Name of the component?
      call optcpp_getComponentName (&
          stemp,rphysics,cspace,rpostproc%p_ItypePointEvalPrimal(1,i),cderiv)
      
      ! Component must be valid, otherwise we cannot evaluate it.
      if (stemp .ne. "") then

        ! Evaluate the function
        call fevl_evaluate (cderiv, Dvalues(1:1), &
            rvector%RvectorBlock(rpostproc%p_ItypePointEvalPrimal(1,i)), &
            rpostproc%p_DcoordsPointEvalPrimal(:,i:i),cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
      
        ! Print th evalue to the terminal
        write (sstr,"(A10,A,F9.4,A,F9.4,A,E17.10)") stemp,&
            "(",rpostproc%p_DcoordsPointEvalPrimal(1,i),",",&
            rpostproc%p_DcoordsPointEvalPrimal(2,i),") = ",Dvalues(1)
            
        call output_line(trim(sstr),coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
      
        ! Write to DAT file?
        if ((rpostproc%sfilenamePointValuesPrimal .ne. "") .and. &
            (rpostproc%iwritePointValues .ne. 0)) then

          ! Open a new file?
          if (iunit .eq. 0) then
          
            ! Write the result to a text file.
            ! Format: timestep current-time value value value ...
            sheadline = "# timestep time x y type deriv value x y type deriv value ..."

            ! Header of the line
            stemp = &
                trim(sys_si0L(isolutionid,5)) // " " // &
                trim(sys_sdEL(dtime,10)) 

            call optcpp_appendLineToDatafile (stemp,sheadline,&
                rpostproc%sfilenamePointValuesPrimal,bfirst,itag,iunit,bnewline=.false.)
          end if

          ! Data to write
          stemp = " " //&
              trim(sys_sdEL(rpostproc%p_DcoordsPointEvalPrimal(1,i),5)) // " " // &
              trim(sys_sdEL(rpostproc%p_DcoordsPointEvalPrimal(2,i),5)) // " " // &
              trim(sys_si0L(rpostproc%p_ItypePointEvalPrimal(1,i),2)) // " " // &
              trim(sys_si0L(rpostproc%p_ItypePointEvalPrimal(2,i),2)) // " " // &
              trim(sys_sdEL(Dvalues(1),10))
              
          call optcpp_appendLineToDatafile (stemp,sheadline,&
              rpostproc%sfilenamePointValuesPrimal,bfirst,itag,iunit,bnewline=.false.)
              
        end if
      
      end if
      
    end do
    
    if (iunit .ne. 0.0_DP) then
      ! Final newline. Close file.
      write (iunit,"(A)") ""
      close (iunit)
    end if
      
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine optcpp_singleEvalForces (rpostproc,rvector,cspace,rphysics,&
      bfirst,isolutionid,dtime,itag)

!<description>
  ! Evaluates forces and other physical data from a vector.
!</description>
  
!<input>
  ! Solution vector
  type(t_vectorBlock), intent(in), target :: rvector
  
  ! Space corresponding to rvector
  integer, intent(in) :: cspace
  
  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics

  ! Id of the solution (Number of the timestep).
  integer, intent(in) :: isolutionid
  
  ! Must be set to TRUE for the first value.
  logical, intent(in) :: bfirst
  
  ! Solution time. =0 for stationary simulations.
  real(DP), intent(in) :: dtime

  ! OPTIONAL: An integer tag which is included into filenames of
  ! output files. If not specified, the tag is not included.
  integer, intent(in), optional :: itag
!</input>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!</subroutine>

    ! Forces on the object
    real(DP), dimension(NDIM2D) :: Dforces,Derr
    real(DP) :: dflux, denergy
    type(t_boundaryRegion) :: rregion
    character(LEN=SYS_STRLEN) :: stemp

    ! Which equation do we have?
    select case (rphysics%cequation)
    
    ! -------------------------------------------------------------
    ! Stokes/Navier Stokes.
    ! -------------------------------------------------------------
    case (CCEQ_STOKES2D,CCEQ_NAVIERSTOKES2D)

      ! Underlying space
      select case (cspace)
      
      ! -------------------------------
      ! Primal space
      ! -------------------------------
      case (CCSPACE_PRIMAL)
      
        ! -------------------------------
        ! Drag/lift forces
        ! -------------------------------
        if (rpostproc%icalcForces .ne. 0) then

          call output_lbrk()
          call output_line ('Body forces real bd., primal space, bdc/horiz/vert')

          ! Calculate drag-/lift coefficients on the 2nd boundary component.
          ! This is for the benchmark channel!
          call boundary_createRegion (rvector%p_rblockDiscr%p_rboundary, &
              2, 0, rregion)
              
          rregion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
          call ppns2D_bdforces_uniform (rvector,rregion,Dforces,CUB_G1_1D,&
              rpostproc%dbdForcesCoeff1,rpostproc%dbdForcesCoeff2)

          ! Output to the terminal        
          call output_line (' 2 / ' &
              //trim(sys_sdEP(Dforces(1),15,6)) // ' / '&
              //trim(sys_sdEP(Dforces(2),15,6)) )
          
          ! Output to a file?
          if (rpostproc%iwriteBodyForces .ne. 0) then
          
            ! Write the result to a text file.
            stemp = trim(sys_si0L(isolutionid,5)) // ' ' &
                // trim(sys_sdEL(dtime,10)) // ' ' &
                // trim(sys_si0L(rpostproc%ibodyForcesBdComponent,5)) // ' ' &
                // trim(sys_sdEL(Dforces(1),10)) // ' '&
                // trim(sys_sdEL(Dforces(2),10))

            call optcpp_appendLineToDatafile (stemp,&
                "# timestep time bdc horiz vert",rpostproc%sfilenameBodyForces,bfirst,itag)
          
          end if
        end if
      
        ! -------------------------------
        ! Flux
        ! -------------------------------
        if (rpostproc%icalcFlux .ne. 0) then
          
          ! Calculate the flux.
          call ppns2D_calcFluxThroughLine (rvector,&
              rpostproc%Dfluxline(1:2),rpostproc%Dfluxline(3:4),dflux)
              
          call output_lbrk()
          call output_line ('Flux, primal space = '//trim(sys_sdEP(dflux,15,6)))

          ! Write the result to a text file?
          if (rpostproc%iwriteflux .ne. 0) then
            ! Format: timestep current-time value
            stemp = trim(sys_si0L(isolutionid,5)) // ' ' &
                // trim(sys_sdEL(dtime,10)) // ' ' &
                // trim(sys_sdEL(dflux,10))

            call optcpp_appendLineToDatafile (stemp,&
                "# timestep time flux",rpostproc%sfilenameFlux,bfirst,itag)
          end if
          
        end if
      
        ! -------------------------------
        ! Kinetic energy
        ! -------------------------------
        if (rpostproc%icalcKineticEnergy .ne. 0) then
        
          ! Perform error analysis to calculate and add 1/2||u||^2_{L^2}.
          call pperr_scalar (PPERR_L2ERROR,Derr(1),rvector%RvectorBlock(1))
          call pperr_scalar (PPERR_L2ERROR,Derr(2),rvector%RvectorBlock(2))
                             
          denergy = 0.5_DP*(Derr(1)**2+Derr(2)**2)

          call output_lbrk()
          call output_line ('||y_1||_L2         = '//trim(sys_sdEP(Derr(1),15,6)))
          call output_line ('||y_2||_L2         = '//trim(sys_sdEP(Derr(2),15,6)))
          call output_line ('||y||_L2           = '//&
              trim(sys_sdEP(sqrt(Derr(1)**2+Derr(2)**2),15,6)))
          call output_line ('1/2||y||^2_L2      = '//trim(sys_sdEP(denergy,15,6)))
          
          ! Write the result to a text file?
          if (rpostproc%iwriteKineticEnergy .ne. 0) then
            stemp = trim(sys_si0L(isolutionid,5)) // ' ' &
                // trim(sys_sdEL(dtime,10)) // ' ' &
                // trim(sys_sdEL(denergy,10)) // ' ' &
                // trim(sys_sdEL(sqrt(Derr(1)**2+Derr(2)**2),10)) // ' ' &
                // trim(sys_sdEL(Derr(1),10)) // ' ' &
                // trim(sys_sdEL(Derr(2),10))

            call optcpp_appendLineToDatafile (stemp,&
                "# timestep time 1/2||u||^2_L2 ||u||_L2 ||u_1||_L2 ||u_2||_L2",&
                rpostproc%sfilenameKineticEnergy,bfirst,itag)
          end if
          
        end if

      ! -------------------------------
      ! Dual space
      ! -------------------------------
      case (CCSPACE_DUAL)

        ! Nothing to do up to now.
      
      end select
      
    ! -------------------------------------------------------------
    ! Heat equation.
    ! -------------------------------------------------------------
    case (CCEQ_HEAT2D)

      ! Nothing to do up to now.

    end select    

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_spaceTimeVisualisation (rpostproc,rsolution,cspace,&
      rphysics,roptcontrol,itag)
  
!<description>
  ! Writes visualisation files for a solution vector.
!</description>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!<input>
  ! Underlying physics
  type(t_settings_physics), intent(in) :: rphysics
  
  ! Underlying optimal control parameters
  type(t_settings_optcontrol), intent(in) :: roptcontrol

  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_spaceTimeVectorAccess), intent(inout) :: rsolution

  ! Space corresponding to rvector
  integer, intent(in) :: cspace
  
  ! OPTIONAL: An integer tag which is included into filenames of
  ! output files. If not specified, the tag is not included.
  integer, intent(in), optional :: itag
!</input>

!</subroutine>

    ! local variables
    integer :: i,iindex
    real(DP) :: dtime, dtimerel
    type(t_vectorBlock), pointer :: p_rvector
    
    ! Just write out the solutions. Corrections concerning the position
    ! of the solutions in time will be added later.
    do i=1,rsolution%p_rspaceTimeVector%NEQtime
    
      select case (cspace)
      
      case (CCSPACE_PRIMAL)
      
        ! Current time
        call tdiscr_getTimestep(rsolution%p_rspaceTimeVector%p_rtimediscr,i-1,dtimeEnd=dtime)
        
        ! Rescale the time.
        call mprim_linearRescale(dtime,&
            rsolution%p_rspaceTimeVector%p_rtimediscr%dtimeInit,&
            rsolution%p_rspaceTimeVector%p_rtimediscr%dtimeMax,&
            0.0_DP,1.0_DP,dtimerel)

        ! Get the solution
        iindex = -1
        call sptivec_getFreeBufferFromPool (rsolution,iindex,p_rvector)
        call sptivec_getTimestepDataByTime (rsolution%p_rspaceTimeVector, dtimerel, p_rvector)
        
        call optcpp_singleVisualisation (rpostproc,p_rvector,cspace,&
            rsolution%p_rspaceTimeVector%p_rspaceDiscr%p_rtriangulation,&
            rphysics,roptcontrol,i-1,dtime,itag)

      case (CCSPACE_DUAL,CCSPACE_CONTROL)
      
        ! Current time
        call tdiscr_getTimestep(rsolution%p_rspaceTimeVector%p_rtimediscr,i,dtimeStart=dtime)
        
        ! Rescale the time.
        call mprim_linearRescale(dtime,&
            rsolution%p_rspaceTimeVector%p_rtimediscr%dtimeInit,&
            rsolution%p_rspaceTimeVector%p_rtimediscr%dtimeMax,&
            0.0_DP,1.0_DP,dtimerel)

        ! Get the solution
        iindex = -1
        call sptivec_getFreeBufferFromPool (rsolution,iindex,p_rvector)
        call sptivec_getTimestepDataByTime (rsolution%p_rspaceTimeVector, dtimerel, p_rvector)
        
        call optcpp_singleVisualisation (rpostproc,p_rvector,cspace,&
            rsolution%p_rspaceTimeVector%p_rspaceDiscr%p_rtriangulation,&
            rphysics,roptcontrol,i-1,dtime,itag)
            
      end select
    
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_spaceTimeDump (rpostproc,rsolution,cspace,itag)
  
!<description>
  ! Dumps a vector in a file.
!</description>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!<input>
  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_spaceTimeVectorAccess), intent(inout) :: rsolution

  ! Space corresponding to rvector
  integer, intent(in) :: cspace
  
  ! OPTIONAL: An integer tag which is included into filenames of
  ! output files. If not specified, the tag is not included.
  integer, intent(in), optional :: itag
!</input>

!</subroutine>

    ! local variables
    integer :: i
    type(t_vectorBlock), pointer :: p_rvector
    
    ! Cancel if there is nothing to dump
    select case (cspace)
    case (CCSPACE_PRIMAL)
      if (rpostproc%cwritePrimalSol .eq. 0) return
    case (CCSPACE_DUAL)
      if (rpostproc%cwriteDualSol .eq. 0) return
    case (CCSPACE_CONTROL)
      if (rpostproc%cwriteControl .eq. 0) return
    end select
    
    ! Just write out the solutions. Corrections concerning the position
    ! of the solutions in time will be added later.
    do i=1,rsolution%p_rspaceTimeVector%NEQtime
    
      ! Get the solution
      call sptivec_getVectorFromPool (rsolution, i, p_rvector)
      
      ! Dump it
      call optcpp_singleVectorDump (rpostproc,p_rvector,cspace,i-1,itag)
    
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_spaceTimeFunctionals (rpostproc,&
      ispacelevel, itimelevel, roperatorAsmHier, &
      rkktsystem)
  
!<description>
  ! Calculates values of the space-time functionals.
!</description>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!<input>
  ! Space-level corresponding to the solutions
  integer, intent(in) :: ispacelevel

  ! Time-level corresponding to the solutions
  integer, intent(in) :: itimelevel

  ! A hierarchy of operator assembly structures for all levels.
  ! This is a central discretisation structure passed to all assembly routines.
  type(t_spacetimeOpAsmHierarchy), intent(in) :: roperatorAsmHier

  ! Structure defining the KKT system.
  type(t_kktsystem), intent(inout) :: rkktsystem
!</input>

!</subroutine>

    ! local variables
    real(DP), dimension(5) :: Derror

    ! Should we calculate the functional?
    if (rpostproc%icalcFunctionalValues .ne. 0) then
      ! Subdomain is the observation area
      call optcana_nonstatFunctional (Derror,&
          ispacelevel, itimelevel, roperatorAsmHier, &
          rkktsystem)

      call output_line ('||y-z||       = '//trim(sys_sdEL(Derror(2),10)))
      call output_line ('||y(T)-z(T)|| = '//trim(sys_sdEL(Derror(3),10)))
      call output_line ('||u||         = '//trim(sys_sdEL(Derror(4),10)))
      call output_line ('||u_Gamma||   = '//trim(sys_sdEL(Derror(5),10)))
      call output_line ('J(y,u)        = '//trim(sys_sdEL(Derror(1),10)))
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_writeSpaceTimeFct (rpostproc,&
      ispacelevel, itimelevel, roperatorAsmHier, &
      rkktsystem,itag)
  
!<description>
  ! Calculates values of the space-time functionals
  ! and writes them into a text file.
!</description>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!<input>
  ! Space-level corresponding to the solutions
  integer, intent(in) :: ispacelevel

  ! Time-level corresponding to the solutions
  integer, intent(in) :: itimelevel

  ! A hierarchy of operator assembly structures for all levels.
  ! This is a central discretisation structure passed to all assembly routines.
  type(t_spacetimeOpAsmHierarchy), intent(in) :: roperatorAsmHier

  ! Structure defining the KKT system.
  type(t_kktsystem), intent(inout) :: rkktsystem

  ! OPTIONAL: An integer tag which is included into filenames of
  ! output files. If not specified, the tag is not included.
  integer, intent(in), optional :: itag
!</input>

!</subroutine>

    ! local variables
    real(DP), dimension(5) :: Derror
    character(len=SYS_STRLEN) :: sfile
    real(DP) :: dtimeend,dtstep,dtimestart
    integer :: iunit
    integer :: idoftime
    type(t_spacetimeOperatorAsm) :: roperatorAsm

    ! Should we calculate the functional?
    if (rpostproc%cwriteFunctionalStatistics .ne. 0) then
    
      ! Get the corresponding operator assembly structure
      call stoh_getOpAsm_slvtlv (&
          roperatorAsm,roperatorAsmHier,ispacelevel,itimelevel)

      sfile = rpostproc%sfunctionalStatisticsFilename
    
      ! Prepare creation of a new file
      iunit = 0
    
      ! Subdomain is the observation area.
      ! Calculate the value in the endpoint of every time interval
      ! as well as in the midpoint. Write them into the file.
      !
      ! Loop over all timesteps (+ 1 for the final solution)
      do idoftime = 1,rkktsystem%p_rprimalSol%p_rvector%NEQtime+1
      
        ! Get the corresponding time interval
        call tdiscr_getTimestep(roperatorAsm%p_rtimeDiscrPrimal,idofTime-1,&
            dtimeend,dtstep,dtimestart)
            
        ! Calculate the values for the beginning and the timestep midpoint.
        ! If we are at the end, only calculate the values for the time endpoint
        ! (which is the beginning of an interval after the final time).
        !
        ! Get the values
        call optcana_nonstatFunctionalAtTime (Derror,dtimestart,&
            ispacelevel, itimelevel, roperatorAsmHier, rkktsystem)
            
        ! Write
        call optcpp_appendLineToDatafile (&
            sys_adjustl(sys_sdE(dtimestart,10),22)//&
            sys_adjustl(sys_sdE(Derror(1),10),22)//&
            sys_adjustl(sys_sdE(Derror(2),10),22)//&
            sys_adjustl(sys_sdE(Derror(3),10),22)//&
            sys_adjustl(sys_sdE(Derror(4),10),22)//&
            sys_adjustl(sys_sdE(Derror(5),10),22),&
            "# time                "//&
            "J(y,u)                "//&
            "||y-z||_{L^2}         "//&
            "||y(T)-z(T)||_{L^2}   "//&
            "||u||_{L^2}           "//&
            "||u||_{L^2(Gamma_C)}  ",&
            sfile,idoftime .eq. 1,&
            itag,iunit)
        ! Centre
        if (idoftime .ne. rkktsystem%p_rprimalSol%p_rvector%NEQtime+1) then
          ! Get the values
          call optcana_nonstatFunctionalAtTime (Derror,0.5_DP*dtimestart + 0.5_DP*dtimeend,&
              ispacelevel, itimelevel, roperatorAsmHier, rkktsystem)
          call optcpp_appendLineToDatafile (&
              sys_adjustl(sys_sdE(0.5_DP*dtimestart + 0.5_DP*dtimeend,10),22)//&
              sys_adjustl(sys_sdE(Derror(1),10),22)//&
              sys_adjustl(sys_sdE(Derror(2),10),22)//&
              sys_adjustl(sys_sdE(Derror(3),10),22)//&
              sys_adjustl(sys_sdE(Derror(4),10),22)//&
              sys_adjustl(sys_sdE(Derror(5),10),22),&
              "",&
              sfile,.false.,&
              itag,iunit)
        end if
        
      end do
      
      ! Close the file
      close (iunit)
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_postprocessing (rpostproc,&
      ispacelevel,itimelevel,rkktsystem,&
      rsettings,itag)
  
!<description>
  ! Postprocessing of a space-time solution.
  ! Writes the solution into a visualisation file, calculates forces,...
!</description>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!<input>
  ! The global settings structure, passed to callback routines.
  type(t_settings_optflow), intent(inout), target :: rsettings

  ! Space-level corresponding to the solutions
  integer, intent(in) :: ispacelevel

  ! Time-level corresponding to the solutions
  integer, intent(in) :: itimelevel

  ! Structure defining the KKT system.
  type(t_kktsystem), intent(inout) :: rkktsystem

  ! OPTIONAL: An integer tag which is included into filenames of
  ! output files. If not specified, the tag is not included.
  integer, intent(in), optional :: itag
!</input>

!</subroutine>

    ! local variables
    
    ! Global postprocessing: Function values.
    call optcpp_spaceTimeFunctionals (rpostproc,&
        ispacelevel, itimelevel, rsettings%roperatorAsmHier, &
        rkktsystem)
        
    call optcpp_writeSpaceTimeFct (rpostproc,&
        ispacelevel, itimelevel, rsettings%roperatorAsmHier, &
        rkktsystem)
        
    call output_lbrk()
    
    ! ***************************************************************************
    ! Primal solution
    ! ***************************************************************************
    
    ! Visualisation
    call optcpp_spaceTimeVisualisation (rpostproc,&
        rkktsystem%p_rprimalSol%p_rvectorAccess,CCSPACE_PRIMAL,&
        rsettings%rphysics,rsettings%rsettingsOptControl,itag)

    ! Data dump
    call optcpp_spaceTimeDump (rpostproc,&
        rkktsystem%p_rprimalSol%p_rvectorAccess,CCSPACE_PRIMAL,itag)

    ! ***************************************************************************
    ! Dual solution
    ! ***************************************************************************
    
    ! Visualisation
    call optcpp_spaceTimeVisualisation (rpostproc,&
        rkktsystem%p_rdualSol%p_rvectorAccess,CCSPACE_DUAL,&
        rsettings%rphysics,rsettings%rsettingsOptControl,itag)

    ! Data dump
    call optcpp_spaceTimeDump (rpostproc,&
        rkktsystem%p_rdualSol%p_rvectorAccess,CCSPACE_DUAL,itag)

    ! ***************************************************************************
    ! Control
    ! ***************************************************************************
    
    ! Visualisation
    call optcpp_spaceTimeVisualisation (rpostproc,&
        rkktsystem%p_rcontrol%p_rvectorAccess,CCSPACE_CONTROL,&
        rsettings%rphysics,rsettings%rsettingsOptControl,itag)

    ! Data dump
    call optcpp_spaceTimeDump (rpostproc,&
        rkktsystem%p_rcontrol%p_rvectorAccess,CCSPACE_CONTROL,itag)
  

!      ! -------------------------------------------------------------------------
!      ! Error analysis
!      ! -------------------------------------------------------------------------
!
!      ! If error analysis has to be performed, we can calculate
!      ! the real error.
!      if (rpostproc%icalcError .eq. 1) then
!        call output_lbrk()
!        call optcana_analyticalError (rsettings%rglobalData,&
!            rsettings%rsettingsOptControl%rconstraints,&
!            rvector,rpostproc%ranalyticRefFunction,&
!            DerrorU,DerrorP,DerrorLambda,DerrorXi,.true.)
!        call output_lbrk()
!        call output_line ('||y-y0||_[0,T]           = '//trim(sys_sdEL(DerrorU(1),10)))
!        call output_line ('||y-y0||_[0,T)           = '//trim(sys_sdEL(DerrorU(2),10)))
!        call output_line ('||y-y0||_(0,T]           = '//trim(sys_sdEL(DerrorU(3),10)))
!        call output_line ('||y-y0||_(0,T)           = '//trim(sys_sdEL(DerrorU(4),10)))
!        
!        call output_line ('||p-p0||_[0,T]           = '//trim(sys_sdEL(DerrorP(1),10)))
!        call output_line ('||p-p0||_[0,T)           = '//trim(sys_sdEL(DerrorP(2),10)))
!        call output_line ('||p-p0||_(0,T]           = '//trim(sys_sdEL(DerrorP(3),10)))
!        call output_line ('||p-p0||_(0,T)           = '//trim(sys_sdEL(DerrorP(4),10)))
!        
!        call output_line ('||lambda-lambda0||_[0,T] = '//trim(sys_sdEL(DerrorLambda(1),10)))
!        call output_line ('||lambda-lambda0||_[0,T) = '//trim(sys_sdEL(DerrorLambda(2),10)))
!        call output_line ('||lambda-lambda0||_(0,T] = '//trim(sys_sdEL(DerrorLambda(3),10)))
!        call output_line ('||lambda-lambda0||_(0,T) = '//trim(sys_sdEL(DerrorLambda(4),10)))
!
!        call output_line ('||xi-xi0||_[0,T]         = '//trim(sys_sdEL(DerrorXi(1),10)))
!        call output_line ('||xi-xi0||_[0,T)         = '//trim(sys_sdEL(DerrorXi(2),10)))
!        call output_line ('||xi-xi0||_(0,T]         = '//trim(sys_sdEL(DerrorXi(3),10)))
!        call output_line ('||xi-xi0||_(0,T)         = '//trim(sys_sdEL(DerrorXi(4),10)))
!      end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_postprocessSubstep (rpostproc,&
      ispacelevel,itimelevel,rkktsystem,&
      rsettings,cpostprocessFlags,itag)
  
!<description>
  ! Postprocessing of a space-time solution after a nonlinear
  ! iteration.
  ! Writes the solution into a visualisation file, calculates forces,...
  ! This routine is assumed to be called in a nonlinear iteration for
  ! the postprocessing of intermediate solutions and uses therefore
  ! slightly different settings.
!</description>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!<input>
  ! The global settings structure, passed to callback routines.
  type(t_settings_optflow), intent(inout), target :: rsettings

  ! Space-level corresponding to the solutions
  integer, intent(in) :: ispacelevel

  ! Time-level corresponding to the solutions
  integer, intent(in) :: itimelevel

  ! Structure defining the KKT system.
  type(t_kktsystem), intent(inout) :: rkktsystem
  
  ! Defines how to apply the postprocessing.
  ! Whether to postprocess intermediate solutions.
  ! =1: Calculate functional values, errors,...
  ! =2: Write postprocessing files with unique filename.
  ! =3: Calculate functional values, errors,... and
  !     write postprocessing files with unique filename.
  integer, intent(in) :: cpostprocessFlags

  ! OPTIONAL: An integer tag which is included into filenames of
  ! output files. If not specified, the tag is not included.
  integer, intent(in), optional :: itag
!</input>

!</subroutine>

    ! local variables
    
    ! Global postprocessing: Function values.
    if (iand(cpostprocessFlags,1) .ne. 0) then
      call optcpp_spaceTimeFunctionals (rpostproc,&
          ispacelevel, itimelevel, rsettings%roperatorAsmHier, &
          rkktsystem)
      call optcpp_writeSpaceTimeFct (rpostproc,&
          ispacelevel, itimelevel, rsettings%roperatorAsmHier, &
          rkktsystem,itag)
    end if
    
    if (iand(cpostprocessFlags,2) .ne. 0) then
    
      if (iand(cpostprocessFlags,1) .ne. 0) &
          call output_lbrk()
    
      ! ***************************************************************************
      ! Primal solution
      ! ***************************************************************************
      
      ! Visualisation
      call optcpp_spaceTimeVisualisation (rpostproc,&
          rkktsystem%p_rprimalSol%p_rvectorAccess,CCSPACE_PRIMAL,&
          rsettings%rphysics,rsettings%rsettingsOptControl,itag)

      ! Data dump
      call optcpp_spaceTimeDump (rpostproc,&
          rkktsystem%p_rprimalSol%p_rvectorAccess,CCSPACE_PRIMAL,itag)

      ! ***************************************************************************
      ! Dual solution
      ! ***************************************************************************
      
      ! Visualisation
      call optcpp_spaceTimeVisualisation (rpostproc,&
          rkktsystem%p_rdualSol%p_rvectorAccess,CCSPACE_DUAL,&
          rsettings%rphysics,rsettings%rsettingsOptControl,itag)

      ! Data dump
      call optcpp_spaceTimeDump (rpostproc,&
          rkktsystem%p_rdualSol%p_rvectorAccess,CCSPACE_DUAL,itag)

      ! ***************************************************************************
      ! Control
      ! ***************************************************************************
      
      ! Visualisation
      call optcpp_spaceTimeVisualisation (rpostproc,&
          rkktsystem%p_rcontrol%p_rvectorAccess,CCSPACE_CONTROL,&
          rsettings%rphysics,rsettings%rsettingsOptControl,itag)
          
    end if

!    ! Data dump
!    call optcpp_spaceTimeDump (rpostproc,rcontrol%p_rvectorAccess,CCSPACE_CONTROL,itag)
  

!      ! -------------------------------------------------------------------------
!      ! Error analysis
!      ! -------------------------------------------------------------------------
!
!      ! If error analysis has to be performed, we can calculate
!      ! the real error.
!      if (rpostproc%icalcError .eq. 1) then
!        call output_lbrk()
!        call optcana_analyticalError (rsettings%rglobalData,&
!            rsettings%rsettingsOptControl%rconstraints,&
!            rvector,rpostproc%ranalyticRefFunction,&
!            DerrorU,DerrorP,DerrorLambda,DerrorXi,.true.)
!        call output_lbrk()
!        call output_line ('||y-y0||_[0,T]           = '//trim(sys_sdEL(DerrorU(1),10)))
!        call output_line ('||y-y0||_[0,T)           = '//trim(sys_sdEL(DerrorU(2),10)))
!        call output_line ('||y-y0||_(0,T]           = '//trim(sys_sdEL(DerrorU(3),10)))
!        call output_line ('||y-y0||_(0,T)           = '//trim(sys_sdEL(DerrorU(4),10)))
!        
!        call output_line ('||p-p0||_[0,T]           = '//trim(sys_sdEL(DerrorP(1),10)))
!        call output_line ('||p-p0||_[0,T)           = '//trim(sys_sdEL(DerrorP(2),10)))
!        call output_line ('||p-p0||_(0,T]           = '//trim(sys_sdEL(DerrorP(3),10)))
!        call output_line ('||p-p0||_(0,T)           = '//trim(sys_sdEL(DerrorP(4),10)))
!        
!        call output_line ('||lambda-lambda0||_[0,T] = '//trim(sys_sdEL(DerrorLambda(1),10)))
!        call output_line ('||lambda-lambda0||_[0,T) = '//trim(sys_sdEL(DerrorLambda(2),10)))
!        call output_line ('||lambda-lambda0||_(0,T] = '//trim(sys_sdEL(DerrorLambda(3),10)))
!        call output_line ('||lambda-lambda0||_(0,T) = '//trim(sys_sdEL(DerrorLambda(4),10)))
!
!        call output_line ('||xi-xi0||_[0,T]         = '//trim(sys_sdEL(DerrorXi(1),10)))
!        call output_line ('||xi-xi0||_[0,T)         = '//trim(sys_sdEL(DerrorXi(2),10)))
!        call output_line ('||xi-xi0||_(0,T]         = '//trim(sys_sdEL(DerrorXi(3),10)))
!        call output_line ('||xi-xi0||_(0,T)         = '//trim(sys_sdEL(DerrorXi(4),10)))
!      end if
    
  end subroutine

end module

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
  use storage
  use genoutput
  use io
  use basicgeometry
  use linearsolver
  use boundary
  use element
  use cubature
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation
  use linearformevaluation
  use matrixfilters
  use vectorfilters
  use discretebc
  use discretefbc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use pprocerror
  use pprocgradients
  use triasearch
  use mprimitives
  use derivatives
  use feevaluation
  use pprocintegrals
  use spacetimelinearsystem
  use spacematvecassembly
  
  use collection
  use convection
  
  use ucd
  use vectorio
  
  use pprocnavierstokes

  use analyticsolution
  use optcanalysis
  use spatialbcdef
  use spacetimevectors
  use spacematvecassembly
  
  use constantsoptc
  use structuresoptc
  use timediscretisation
  use structuresoptflow
  use user_callback
  
  !use spacetimediscretisation
  
  implicit none
  
  private
  
  public :: t_optcPostprocessing  
  public :: optcpp_initpostprocessing
  public :: optcpp_donepostprocessing
  public :: optcpp_postprocessSingleSol
  public :: optcpp_postprocessSpaceTimeVec
  public :: optcpp_postprocSpaceVisOutput
 
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
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_initpostprocessing (rpostproc,rphysics,cspace,rboundaryConditions,&
      rtimeDiscr,rspaceDiscr,rspaceDiscrPrimal)
  
!<description>
  ! Initialises the postprocessing
!</description>

!<input>
  ! Physics of the problem
  type(t_settings_physics), intent(in), target :: rphysics

  ! Space that is available in rsolution. One of the CCSPACE_xxxx constants.
  integer, intent(in) :: cspace
  
  ! Boundary conditions to use.
  type(t_optcBDC), intent(in), target  :: rboundaryConditions

  ! Underlying space discretisation
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr

  ! Underlying space discretisation in the primal space
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscrPrimal

  ! Underlying time discretisation
  type(t_timeDiscretisation), intent(in), target :: rtimeDiscr
!</input>
  
!<inputoutput>
  ! Parameters about the postprocessing.
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>
  
!</subroutine>

    ! Fetch data
    rpostproc%cspace = cspace
    rpostproc%p_rspaceDiscrPrimal => rspaceDiscrPrimal
    rpostproc%p_rspaceDiscr => rspaceDiscr
    rpostproc%p_rtimeDiscr => rtimeDiscr
    rpostproc%p_rboundaryConditions => rboundaryConditions
    rpostproc%p_rphysics => rphysics
    
    select case (rphysics%cequation)
    case (0,1)
      ! Stokes, Navier-Stokes, 2D

      ! Initialise a P1/Q1/P0/Q0 structure for visualisationoutput.
      call spdiscr_duplicateBlockDiscr(rpostproc%p_rspaceDiscr,rpostproc%rspaceDiscrLinear)
      
      call spdiscr_deriveDiscr_triquad (&
                  rpostproc%p_rspaceDiscr%RspatialDiscr(1), &
                  EL_P1, EL_Q1, CUB_TRZ_T, CUB_G2X2, &
                  rpostproc%rspaceDiscrLinear%RspatialDiscr(1))

      call spdiscr_deriveDiscr_triquad (&
                  rpostproc%p_rspaceDiscr%RspatialDiscr(2), &
                  EL_P1, EL_Q1, CUB_TRZ_T, CUB_G2X2, &
                  rpostproc%rspaceDiscrLinear%RspatialDiscr(2))

      call spdiscr_deriveDiscr_triquad (&
                  rpostproc%p_rspaceDiscr%RspatialDiscr(3), &
                  EL_P0, EL_Q0, CUB_TRZ_T, CUB_G2X2, &
                  rpostproc%rspaceDiscrLinear%RspatialDiscr(3))

      call spdiscr_deriveDiscr_triquad (&
                  rpostproc%p_rspaceDiscr%RspatialDiscr(4), &
                  EL_P1, EL_Q1, CUB_TRZ_T, CUB_G2X2, &
                  rpostproc%rspaceDiscrLinear%RspatialDiscr(4))

      call spdiscr_deriveDiscr_triquad (&
                  rpostproc%p_rspaceDiscr%RspatialDiscr(5), &
                  EL_P1, EL_Q1, CUB_TRZ_T, CUB_G2X2, &
                  rpostproc%rspaceDiscrLinear%RspatialDiscr(5))

      call spdiscr_deriveDiscr_triquad (&
                  rpostproc%p_rspaceDiscr%RspatialDiscr(6), &
                  EL_P0, EL_Q0, CUB_TRZ_T, CUB_G2X2, &
                  rpostproc%rspaceDiscrLinear%RspatialDiscr(6))
    end select

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

    ! Release the discretisations structure
    call spdiscr_releaseBlockDiscr (rpostproc%rspaceDiscrLinear)

    ! Clean up other data.
    rpostproc%cspace = CCSPACE_PRIMAL
    nullify(rpostproc%p_rspaceDiscr)
    nullify(rpostproc%p_rtimeDiscr)
    nullify(rpostproc%p_rboundaryConditions)

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

!******************************************************************************

!<subroutine>

  subroutine optcpp_evaluatePoints (rpostproc,rsolution,bfirst,itimestep,&
      dtimePrimal,dtimeDual)

!<description>
  ! Evaluates the solution in a number of points as configured in the DAT file.
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(in), target :: rsolution
  
  ! Number of the timestep.
  integer, intent(in) :: itimestep
  
  ! Must be set to TRUE for the first value.
  logical, intent(in) :: bfirst
  
  ! Solution time. =0 for stationary simulations. Primal equation
  real(DP), intent(in) :: dtimePrimal

  ! Solution time. =0 for stationary simulations. Primal equation
  real(DP), intent(in) :: dtimeDual
!</input>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: npoints,i,iunit
    integer :: iderType, cflag
    logical :: bfileExists
    real(dp), dimension(:), allocatable :: Dvalues
    character(LEN=SYS_STRLEN) :: sstr,stemp
    character(LEN=10), dimension(3,6), parameter :: Sfctnames = reshape (&
      (/ "       u1 ","     u1_x ","     u1_y " , &
         "       u2 ","     u2_x ","     u2_y " , &
         "        p ","      p_x ","      p_y " , &
         "  lambda1 ","lambda1_x ","lambda1_y " , &
         "  lambda2 ","lambda2_x ","lambda2_y " , &
         "       xi ","     xi_x ","     xi_y " /) ,&
       (/ 3,6 /) )

    if (.not. associated(rpostproc%p_DcoordsPointEval)) return
    
    ! Currently, only Stokes/Navier Stokes supported.
    select case (rpostproc%p_rphysics%cequation)
    
    case (0,1)

      ! Get the number of points to evaluate
      npoints = ubound(rpostproc%p_DcoordsPointEval,2)
          
      ! Allocate memory for the values
      allocate(Dvalues(npoints))
      
      ! Evaluate the function in these points.
      do i=1,npoints
        ! Which derivative to evaluate?
        select case (rpostproc%p_ItypePointEval(2,i))
        case (0)
          iderType = DER_FUNC2D
        case (1)
          iderType = DER_DERIV2D_X
        case (2)
          iderType = DER_DERIV2D_Y
        case default
          iderType = DER_FUNC2D
        end select
        
        call fevl_evaluate (iderType, Dvalues(i:i), &
            rsolution%RvectorBlock(rpostproc%p_ItypePointEval(1,i)), &
            rpostproc%p_DcoordsPointEval(:,i:i),cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
      end do
      
      ! Print the values to the terminal
      call output_lbrk()
      call output_line ('Point values')
      call output_line ('------------')
      do i=1,npoints
        write (sstr,"(A10,A,F9.4,A,F9.4,A,E16.10)") &
            Sfctnames(1+rpostproc%p_ItypePointEval(2,i),&
            rpostproc%p_ItypePointEval(1,i)),&
            "(",rpostproc%p_DcoordsPointEval(1,i),",",rpostproc%p_DcoordsPointEval(2,i),") = ",Dvalues(i)
        call output_line(trim(sstr),coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      end do
      
      ! Write to DAT file?
      if ((rpostproc%sfilenamePointValues .eq. "") .or. (rpostproc%iwritePointValues .eq. 0)) then
        deallocate(Dvalues)
        return
      end if
      
      ! When writing to a file is enabled, delete the file in the first timestep.
      cflag = SYS_APPEND
      if (bfirst) cflag = SYS_REPLACE
      
      ! Write the result to a text file.
      ! Format: timestep current-time value value value ...
      call io_openFileForWriting(rpostproc%sfilenamePointValues, iunit, &
          cflag, bfileExists,.true.)
      if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
        ! Write a headline
        write (iunit,'(A)') &
          '# timestep time(primal) time(dual) x y type deriv value x y type deriv value ...'
      end if
      
      stemp = &
          trim(sys_siL(itimestep,10)) // ' ' // &
          trim(sys_sdEL(dtimePrimal,10)) // ' ' // &
          trim(sys_sdEL(dtimeDual,10))
      write (iunit,ADVANCE='NO',FMT='(A)') trim(stemp)
      do i=1,npoints
        stemp = ' ' //&
            trim(sys_sdEL(rpostproc%p_DcoordsPointEval(1,i),5)) // ' ' // &
            trim(sys_sdEL(rpostproc%p_DcoordsPointEval(2,i),5)) // ' ' // &
            trim(sys_siL(rpostproc%p_ItypePointEval(1,i),2)) // ' ' // &
            trim(sys_siL(rpostproc%p_ItypePointEval(2,i),2)) // ' ' // &
            trim(sys_sdEL(Dvalues(i),10))
        write (iunit,ADVANCE='NO',FMT='(A)') trim(stemp)
      end do
      write (iunit,ADVANCE='YES',FMT='(A)') ""
      close (iunit)
      
      deallocate(Dvalues)
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_postprocessSingleSol (rpostproc,roptcontrol,rsettings,&
      ifileid,bfirstFile,blastFile,&
      dtimePrimal,dtimeDual,rvector,rrhs,&
      dtimePrimalPrev,dtimeDualPrev,rvectorPrev,rrhsPrev,&
      dtimePrimalNext,dtimeDualNext,rvectorNext,rrhsNext)
  
!<description>
  ! Postprocessing of a single solution at a definite time.
  ! Writes the solution into a visualisation file, calculates forces,...
!</description>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!<input>
  ! The global settings structure, passed to callback routines.
  type(t_settings_optflow), intent(inout), target :: rsettings

  ! Id of the file on the hard disc; added to the filename.
  integer, intent(in) :: ifileid
  
  ! Parameters about the optimal control
  type(t_settings_optcontrol), intent(in) :: roptControl
  
  ! Must be set to TRUE for the first file of a sequence of files or
  ! of the file does not belong to a sequence.
  logical, intent(in) :: bfirstFile

  ! Must be set to TRUE for the last file of a sequence of files or
  ! of the file does not belong to a sequence.
  logical, intent(in) :: blastFile

  ! Current simulation time in the primal and dual equation.
  real(dp), intent(in) :: dtimePrimal
  real(dp), intent(in) :: dtimeDual
  
  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_vectorBlock), intent(inout) :: rvector

  ! The RHS of the system.
  type(t_vectorBlock), intent(inout) :: rrhs
  
  ! Solution vector of the previous timestep. May coincide with
  ! rvector if there is none; then, bfirstFile=.true. 
  type(t_vectorBlock), intent(inout) :: rvectorPrev

  ! The RHS of the previous timestep.
  type(t_vectorBlock), intent(inout) :: rrhsPrev

  ! Simulation time in the primal and dual equation of the previous timestep
  real(dp), intent(in) :: dtimePrimalPrev
  real(dp), intent(in) :: dtimeDualPrev

  ! Solution vector of the next timestep. May coincide with
  ! rvector if there is none; then blastFile=.true.
  type(t_vectorBlock), intent(inout) :: rvectorNext

  ! The RHS of the next timestep.
  type(t_vectorBlock), intent(inout) :: rrhsNext

  ! Simulation time in the primal and dual equation of the next timestep
  real(dp), intent(in) :: dtimePrimalNext
  real(dp), intent(in) :: dtimeDualNext
!</input>

!</subroutine>

  ! local variables
  integer(I32) :: celement
  integer :: cflag,iunit
  real(DP) :: dflux,denergy
  logical :: bfileexists

  ! We need some more variables for postprocessing - i.e. writing
  ! a GMV file.
  real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2

  ! A pointer to the triangulation.
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! A vector accepting Q1 data
  type(t_vectorBlock) :: rprjVector,rcontrolVector
  
  ! Output block for UCD output to GMV file
  type(t_ucdExport) :: rexport
  
  ! Forces on the object
  real(DP), dimension(NDIM2D) :: Dforces,Derr
  type(t_boundaryRegion) :: rregion
  
  ! Divergence
  !type(t_vectorScalar), target :: rtempVector
  
  character(SYS_STRLEN) :: sfilename,stemp

  ! Discrete boundary conditions
  type(t_discreteBC) :: rdiscreteBC
  
  ! Discrete fictitious boundary conditions
  type(t_discreteFBC) :: rdiscreteFBC
  
    ! Local variables
    type(t_vectorBlock) :: rvectemp

    select case (rpostproc%p_rphysics%cequation)
    case (0,1)
      ! Stokes, Navier-Stokes, 2D
      
      ! -------------------------------------------------------------------------
      ! Body forces & flux calculation
      ! -------------------------------------------------------------------------

      ! When writing to a file is enabled, delete the file in the first timestep.
      cflag = SYS_APPEND
      if (bfirstFile) cflag = SYS_REPLACE

      ! If we have a uniform discretisation, calculate the body forces on the
      ! 2nd boundary component - if it exists.
      if ((rpostproc%p_rspaceDiscr%RspatialDiscr(1)% &
          ccomplexity .eq. SPDISC_UNIFORM) .and. &
          (boundary_igetNBoundComp(rvector%p_rblockDiscr%p_rboundary) .ge. &
              rpostproc%ibodyForcesBdComponent) .and.&
          (rpostproc%icalcForces .ne. 0)) then

        call output_lbrk()
        call output_line ('Body forces real bd., bdc/horiz/vert')

        ! Create a temp vector that contains the solution at the desired time.
        ! We evaluate in the dual time, so for CN this is in the time midpoint.
        call lsysbl_createVectorBlock (rvector,rvectemp)
        call optcpp_unshiftSubvectors (rpostproc,1,&
            dtimePrimalPrev,dtimeDualPrev,rvectorPrev,&
            dtimePrimal,dtimeDual,rvector,&
            dtimePrimalNext,dtimeDualNext,rvectorNext,&
            dtimeDual,rvectemp)

        ! Calculate drag-/lift coefficients on the 2nd boundary component.
        ! This is for the benchmark channel!
        call boundary_createRegion (rvector%p_rblockDiscr%p_rboundary, &
            2, 0, rregion)
        rregion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
        call ppns2D_bdforces_uniform (rvectemp,rregion,Dforces,CUB_G1_1D,&
            rpostproc%dbdForcesCoeff1,rpostproc%dbdForcesCoeff2)
            
        ! Release the temp vector
        call lsysbl_releaseVector (rvectemp)
        
        call output_line (' 2 / ' &
            //trim(sys_sdEP(Dforces(1),15,6)) // ' / '&
            //trim(sys_sdEP(Dforces(2),15,6)) )
        
        if (rpostproc%iwriteBodyForces .ne. 0) then
          ! Write the result to a text file.
          ! Format: timestep current-time value
          call io_openFileForWriting(rpostproc%sfilenameBodyForces, iunit, &
              cflag, bfileExists,.true.)
          if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
            ! Write a headline
            write (iunit,'(A)') '# timestep time bdc horiz vert'
          end if
          stemp = trim(sys_siL(ifileid,10)) // ' ' &
              // trim(sys_sdEL(dtimeDual,10)) // ' ' &
              // trim(sys_siL(rpostproc%ibodyForcesBdComponent,10)) // ' ' &
              // trim(sys_sdEL(Dforces(1),10)) // ' '&
              // trim(sys_sdEL(Dforces(2),10))
          write (iunit,'(A)') trim(stemp)
          close (iunit)
        end if
        
      endif
      
      if (rpostproc%icalcFlux .ne. 0) then
        
        ! Calculate the flux.
        call ppns2D_calcFluxThroughLine (rvector,&
            rpostproc%Dfluxline(1:2),rpostproc%Dfluxline(3:4),dflux)
            
        call output_lbrk()
        call output_line ('Flux = '//trim(sys_sdEP(dflux,15,6)))

        if (rpostproc%iwriteflux .ne. 0) then
          ! Write the result to a text file.
          ! Format: timestep current-time value
          call io_openFileForWriting(rpostproc%sfilenameFlux, iunit, &
              cflag, bfileExists,.true.)
          if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
            ! Write a headline
            write (iunit,'(A)') '# timestep time flux'
          end if
          stemp = trim(sys_siL(ifileid,10)) // ' ' &
              // trim(sys_sdEL(dtimePrimal,10)) // ' ' &
              // trim(sys_sdEL(dflux,10))
          write (iunit,'(A)') trim(stemp)
          close (iunit)
        end if
        
      end if
      
      if (rpostproc%icalcKineticEnergy .ne. 0) then
      
        ! Perform error analysis to calculate and add 1/2||u||^2_{L^2}.
        call pperr_scalar (rvector%RvectorBlock(1),PPERR_L2ERROR,Derr(1))
        call pperr_scalar (rvector%RvectorBlock(2),PPERR_L2ERROR,Derr(2))
                           
        denergy = 0.5_DP*(Derr(1)**2+Derr(2)**2)

        call output_lbrk()
        call output_line ('||u_1||_L2         = '//trim(sys_sdEP(Derr(1),15,6)))
        call output_line ('||u_2||_L2         = '//trim(sys_sdEP(Derr(2),15,6)))
        call output_line ('||u||_L2           = '//&
            trim(sys_sdEP(sqrt(Derr(1)**2+Derr(2)**2),15,6)))
        call output_line ('1/2||u||^2_L2      = '//trim(sys_sdEP(denergy,15,6)))
        
        if (rpostproc%iwriteKineticEnergy .ne. 0) then
          ! Write the result to a text file.
          ! Format: timestep current-time value
          call io_openFileForWriting(rpostproc%sfilenameKineticEnergy, iunit, &
              cflag, bfileExists,.true.)
          if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
            ! Write a headline
            write (iunit,'(A)') '# timestep time 1/2||u||^2_L2 ||u||_L2 ||u_1||_L2 ||u_2||_L2'
          end if
          stemp = trim(sys_siL(ifileid,10)) // ' ' &
              // trim(sys_sdEL(dtimePrimal,10)) // ' ' &
              // trim(sys_sdEL(denergy,10)) // ' ' &
              // trim(sys_sdEL(sqrt(Derr(1)**2+Derr(2)**2),10)) // ' ' &
              // trim(sys_sdEL(Derr(1),10)) // ' ' &
              // trim(sys_sdEL(Derr(2),10))
          write (iunit,'(A)') trim(stemp)
          close (iunit)
        end if

      end if
      
      ! -------------------------------------------------------------------------
      ! Writing out of the final solution
      ! -------------------------------------------------------------------------
      
      if (rpostproc%sfinalSolutionFileName .ne. "") then
        ! Write the current solution to disc as it is.
        sfilename = trim(rpostproc%sfinalSolutionFileName)//'.'//sys_si0(ifileid,5)
        
        call output_lbrk ()
        call output_line ('Writing solution file: '//trim(sfilename))
        
        call vecio_writeBlockVectorHR (rvector, "vector"//sys_si0(ifileid,5), .false.,&
            0, sfilename,  "(E20.10)")
      end if

      ! -------------------------------------------------------------------------
      ! Writing out of the final coptrol. Note that we have to calculate it
      ! from the final dual solution.
      ! -------------------------------------------------------------------------
      
      if (rpostproc%sfinalControlFileName .ne. "") then
      
        ! Do the projection into a vector in the discretisation
        ! if the primal space.
        call lsysbl_createVecBlockByDiscr (rpostproc%p_rspaceDiscrPrimal,&
            rcontrolVector,.true.)
      
        call optcpp_calcControl (rpostproc%p_rphysics,rvector,roptControl%rconstraints,&
            roptControl%dalphaC,dtimeDual,roptcontrol%ispaceTimeFormulation,rcontrolVector)
      
        ! Write the current solution to disc as it is.
        sfilename = trim(rpostproc%sfinalControlFileName)//'.'//sys_si0(ifileid,5)
        
        call output_lbrk ()
        call output_line ('Writing control file: '//trim(sfilename))
        
        call vecio_writeBlockVectorHR (rcontrolVector, "vector"//sys_si0(ifileid,5), .false.,&
            0, sfilename,  "(E20.10)")
            
        call lsysbl_releaseVector (rcontrolVector)
      end if

      ! -------------------------------------------------------------------------
      ! Visualisation output
      ! -------------------------------------------------------------------------

      if (rpostproc%ioutputUCD .ne. 0) then

        ! If we have a simple Q1~ discretisation, calculate the streamfunction.
  !      IF (rvector%p_rblockDiscretisation%RspatialDiscr(1)% &
  !          ccomplexity .EQ. SPDISC_UNIFORM) THEN
  !          
  !        celement = rvector%p_rblockDiscretisation%RspatialDiscr(1)% &
  !                  RelementDistr(1)%itrialElement
  !                  
  !        IF (elem_getPrimaryElement(celement) .EQ. EL_Q1T) THEN
  !        
  !          ! Create a temporary vector 
  !          CALL lsyssc_createVecByDiscr (rvector%RvectorBlock(3)%p_rspatialDiscretisation,&
  !              rtempVector,.TRUE.)
  !
  !          ! Calculate divergence = B1^T u1 + B2^T u2
  !          CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(rproblem%nlmax)%rmatrixB1,&
  !              rBmatrix,LSYSSC_TR_VIRTUAL)
  !          CALL lsyssc_scalarMatVec (&
  !              rBmatrix, rvector%RvectorBlock(1), &
  !              rtempVector, 1.0_DP, 0.0_DP)
  !          CALL lsyssc_transposeMatrix (rproblem%RlevelInfo(rproblem%nlmax)%rmatrixB2,&
  !              rBmatrix,LSYSSC_TR_VIRTUAL)
  !          CALL lsyssc_scalarMatVec (&
  !              rBmatrix, rvector%RvectorBlock(2), &
  !              rtempVector, 1.0_DP, 1.0_DP)
  !          
  !          CALL output_lbrk()
  !          CALL output_line ('Divergence = ' &
  !              //TRIM(sys_sdEP(lsyssc_vectorNorm(rtempVector,LINALG_NORML2),15,6)) )
  !              
  !          CALL lsyssc_releaseVector (rtempVector)
  !        
  !        END IF
  !        
  !      END IF    

        ! Initialise boundary condition structures    
        call bcasm_initDiscreteBC(rdiscreteBC)
        call bcasm_initDiscreteFBC(rdiscreteFBC)

        ! The pressure discretisation substructure stays the old.
        !
        ! Now set up a new solution vector based on this discretisation,
        ! allocate memory.
        call lsysbl_createVecBlockByDiscr (rpostproc%rspaceDiscrLinear,rprjVector,.false.)
        
        ! Then take our original solution vector and convert it according to the
        ! new discretisation:
        call spdp_projectSolution (rvector,rprjVector)
        
        call sbc_assembleBDconditions (rsettings%roptcBDC,dtimePrimal,dtimeDual,rpostproc%rspaceDiscrLinear,&
            rpostproc%p_rtimeDiscr,CCSPACE_PRIMALDUAL,rdiscreteBC,rsettings%rglobalData)
        call sbc_assembleFBDconditions (dtimePrimal,rpostproc%rspaceDiscrLinear,rpostproc%p_rtimeDiscr,&
            CCSPACE_PRIMALDUAL,rdiscreteFBC,rsettings%rglobalData)
        
        ! Filter the solution vector to implement discrete BC's.
        call vecfil_discreteBCsol (rprjVector,rdiscreteBC)

        ! Filter the solution vector to implement discrete BC's for fictitious 
        ! boundary components.
        call vecfil_discreteFBCsol (rprjVector,rdiscreteFBC)
        
        ! Release boundary condition structures.
        call bcasm_releaseDiscreteBC(rdiscreteBC)
        call bcasm_releaseDiscreteFBC(rdiscreteFBC)
        
        ! Now we have a Q1/Q1/Q0 solution in rprjVector.
        !
        ! From the attached discretisation, get the underlying triangulation
        p_rtriangulation => rpostproc%p_rspaceDiscr%p_rtriangulation
        
        ! Check if we have a filename where to write GMV output to.
        if (rpostproc%sfilenameUCD .ne. "") then
        
          ! Start UCD export:
          sfilename = trim(rpostproc%sfilenameUCD)//'.'//sys_si0(ifileid,5)
          
          call output_lbrk ()
          call output_line ('Writing visualisation file: '//trim(sfilename))
          
          select case (rpostproc%ioutputUCD)
          case (1)
            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfilename)

          case (2)
            call ucd_startAVS (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfilename)
                
          case (3)
            call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfilename)
                
          case default
            call output_line ('Invalid UCD ooutput type.', &
                              OU_CLASS_ERROR,OU_MODE_STD,'fbsim_writeUCD')
            stop
          end select
          
          ! Write the configuration of the application as comment block
          ! to the output file.
          call ucd_addCommentLine (rexport,'Configuration:')
          call ucd_addCommentLine (rexport,'---------------')
          call ucd_addParameterList (rexport,rsettings%p_rparlist)
          call ucd_addCommentLine (rexport,'---------------')

          ! Write velocity field
          call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
          call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
          
          call ucd_addVarVertBasedVec (rexport,'velocity_p',&
              p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
          
          ! Write out cell based or node based pressure.
          call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
          celement = rprjVector%p_rblockDiscr%RspatialDiscr(3)% &
                    RelementDistr(1)%celement
                    
          if ((elem_getPrimaryElement(celement) .eq. EL_Q1) .or. &
              ((elem_getPrimaryElement(celement) .eq. EL_P1))) then
            call ucd_addVariableVertexBased (rexport,'pressure_p',UCD_VAR_STANDARD, &
                p_Ddata(1:p_rtriangulation%NVT))
          else
            call ucd_addVariableElementBased (rexport,'pressure_p',UCD_VAR_STANDARD, &
                p_Ddata(1:p_rtriangulation%NEL))
          end if
          
          ! Dual velocity field
          call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
          call lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata2)

          call ucd_addVarVertBasedVec (rexport,'velocity_d',&
              p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
          
          ! Write out cell based or node based dual pressure.
          call lsyssc_getbase_double (rprjVector%RvectorBlock(6),p_Ddata)
          celement = rprjVector%p_rblockDiscr%RspatialDiscr(6)% &
                    RelementDistr(1)%celement
                    
          if ((elem_getPrimaryElement(celement) .eq. EL_Q1) .or. &
              ((elem_getPrimaryElement(celement) .eq. EL_P1))) then
            call ucd_addVariableVertexBased (rexport,'pressure_d',UCD_VAR_STANDARD, &
                p_Ddata(1:p_rtriangulation%NVT))
          else
            call ucd_addVariableElementBased (rexport,'pressure_d',UCD_VAR_STANDARD, &
                p_Ddata(1:p_rtriangulation%NEL))
          end if

          ! Control u = P[min/max](-1/alpha lambda)
          call optcpp_calcControl (rpostproc%p_rphysics,rprjVector,roptControl%rconstraints,&
              roptControl%dalphaC,dtimeDual,roptcontrol%ispaceTimeFormulation)
          
          call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
          call lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata2)
          call ucd_addVarVertBasedVec (rexport,'control',&
              p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
          
          ! If we have a simple Q1~ discretisation, calculate the streamfunction.
          if (rvector%p_rblockDiscr%RspatialDiscr(1)%ccomplexity .eq. SPDISC_UNIFORM) then
              
            celement = rpostproc%p_rspaceDiscr%RspatialDiscr(1)%RelementDistr(1)%celement
                      
            if (elem_getPrimaryElement(celement) .eq. EL_Q1T) then
                
              call ppns2D_streamfct_uniform (rvector,rprjVector%RvectorBlock(1))
              
              call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
              call ucd_addVariableVertexBased (rexport,'streamfunction',&
                  UCD_VAR_STANDARD, p_Ddata)
                  
            end if
            
          end if

          ! Include the RHS.        
          call spdp_projectSolution (rrhs,rprjVector)
          call lsyssc_getbase_double (rprjVector%RvectorBlock(4),p_Ddata)
          call lsyssc_getbase_double (rprjVector%RvectorBlock(5),p_Ddata2)
          call ucd_addVarVertBasedVec (rexport,'rhs',&
              p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
          
          ! Write the file to disc, that's it.
          call ucd_write (rexport)
          call ucd_release (rexport)
          
        end if
        
        ! Release the auxiliary vector
        call lsysbl_releaseVector (rprjVector)
        
      end if
      
    end select
    
    ! Evaluate point values.
    call optcpp_evaluatePoints (rpostproc,rvector,bfirstFile,ifileid,dtimePrimal,dtimeDual)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_postprocessSpaceTimeVec (rpostproc,rvector,rrhs,roptcontrol,rsettings)
  
!<description>
  ! Postprocessing of a space-time vector.
  ! Writes the solution into a visualisation file, calculates forces,...
!</description>

!<inputoutput>
  ! Postprocessing structure
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>

!<input>
  ! The global settings structure, passed to callback routines.
  type(t_settings_optflow), intent(inout), target :: rsettings

  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_spaceTimeVector), intent(IN) :: rvector

  ! The RHS vector of the system.
  type(t_spaceTimeVector), intent(IN) :: rrhs
  
  ! Parameters about the optimal control
  type(t_settings_optcontrol), intent(inout) :: roptControl
!</input>

!</subroutine>

  ! local variables
  type(t_vectorBlock), pointer :: p_rvecTemp1,p_rrhsTemp1
  type(t_vectorBlock), pointer :: p_rvecTemp2,p_rrhsTemp2
  type(t_vectorBlock), pointer :: p_rvecTemp3,p_rrhsTemp3
  integer :: istep
  real(dp) :: dtimePrimal,dtimeDual,dtstep
  real(dp) :: dtimePrimalNext,dtimeDualNext
  real(dp) :: dtimePrimalPrev,dtimeDualPrev
  real(DP), dimension(4) :: DerrorU,DerrorP,DerrorLambda,DerrorXi,Derror
  type(t_spaceTimeVectorAccess) :: raccessPoolSol,raccessPoolRhs

    select case (rpostproc%p_rphysics%cequation)
    case (0,1)
      ! Stokes, Navier-Stokes, 2D
      ! -------------------------------------------------------------------------
      ! Visualisation output, force calculation
      ! -------------------------------------------------------------------------

      ! Create a vector access pool for the RHS and the solution.
      call sptivec_createAccessPool (rvector,raccessPoolSol,3)
      call sptivec_createAccessPool (rrhs,raccessPoolRhs,3)

      ! Write a file for every timestep
      do istep = 1,rvector%NEQtime

        ! Get the solution and RHS of the current timestep
        call tdiscr_getTimestep(rpostproc%p_rtimediscr,istep-1,dtimePrimal,dtstep)
        dtimeDual = dtimePrimal - (1.0_DP-rpostproc%p_rtimeDiscr%dtheta)*dtstep
        call sptivec_getVectorFromPool (raccessPoolSol,istep,p_rvecTemp2)
        call sptivec_getVectorFromPool (raccessPoolRhs,istep,p_rrhsTemp2)
        
        ! Previous/next timestep.
        if (istep .eq. 1) then
          ! No previous timestep.
          dtimePrimalPrev = dtimePrimal
          dtimeDualPrev = dtimeDual
          p_rvecTemp1 => p_rvecTemp2
          p_rrhsTemp1 => p_rrhsTemp2
        else
          ! Get the solution and RHS of the current timestep
          call tdiscr_getTimestep(rpostproc%p_rtimediscr,istep-2,dtimePrimalPrev,dtstep)
          dtimeDualPrev = dtimePrimalPrev - (1.0_DP-rpostproc%p_rtimeDiscr%dtheta)*dtstep
          call sptivec_getVectorFromPool (raccessPoolSol,istep-1,p_rvecTemp1)
          call sptivec_getVectorFromPool (raccessPoolRhs,istep-1,p_rrhsTemp1)
        end if

        if (istep .eq. rvector%NEQtime) then
          ! No next timestep
          dtimePrimalNext = dtimePrimal
          dtimeDualNext = dtimeDual
          p_rvecTemp3 => p_rvecTemp2
          p_rrhsTemp3 => p_rrhsTemp2
        else
          ! Get the solution and RHS of the current timestep
          call tdiscr_getTimestep(rpostproc%p_rtimediscr,istep,dtimePrimalNext,dtstep)
          dtimeDualNext = dtimePrimalNext - (1.0_DP-rpostproc%p_rtimeDiscr%dtheta)*dtstep
          call sptivec_getVectorFromPool (raccessPoolSol,istep+1,p_rvecTemp3)
          call sptivec_getVectorFromPool (raccessPoolRhs,istep+1,p_rrhsTemp3)
        end if
        
        call optcpp_postprocessSingleSol (rpostproc,roptControl,rsettings,&
            istep-1,istep .eq. 1,istep .eq. rvector%NEQtime,&
            dtimePrimal,dtimeDual,p_rvecTemp2,p_rrhsTemp2,&
            dtimePrimalPrev,dtimeDualPrev,p_rvecTemp1,p_rrhsTemp1,&
            dtimePrimalNext,dtimeDualNext,p_rvecTemp3,p_rrhsTemp3)
      end do
      
      call sptivec_releaseAccessPool (raccessPoolRhs)
      call sptivec_releaseAccessPool (raccessPoolSol)
      
      ! -------------------------------------------------------------------------
      ! Error analysis
      ! -------------------------------------------------------------------------

      ! If error analysis has to be performed, we can calculate
      ! the real error.
      if (rpostproc%icalcError .eq. 1) then
        call output_lbrk()
        call optcana_analyticalError (rsettings%rglobalData,&
            rsettings%rsettingsOptControl%rconstraints,&
            rvector,rpostproc%ranalyticRefFunction,&
            DerrorU,DerrorP,DerrorLambda,DerrorXi,.true.)
        call output_lbrk()
        call output_line ('||y-y0||_[0,T]           = '//trim(sys_sdEL(DerrorU(1),10)))   
        call output_line ('||y-y0||_[0,T)           = '//trim(sys_sdEL(DerrorU(2),10)))   
        call output_line ('||y-y0||_(0,T]           = '//trim(sys_sdEL(DerrorU(3),10)))   
        call output_line ('||y-y0||_(0,T)           = '//trim(sys_sdEL(DerrorU(4),10)))   
        
        call output_line ('||p-p0||_[0,T]           = '//trim(sys_sdEL(DerrorP(1),10)))   
        call output_line ('||p-p0||_[0,T)           = '//trim(sys_sdEL(DerrorP(2),10)))   
        call output_line ('||p-p0||_(0,T]           = '//trim(sys_sdEL(DerrorP(3),10)))   
        call output_line ('||p-p0||_(0,T)           = '//trim(sys_sdEL(DerrorP(4),10)))   
        
        call output_line ('||lambda-lambda0||_[0,T] = '//trim(sys_sdEL(DerrorLambda(1),10)))   
        call output_line ('||lambda-lambda0||_[0,T) = '//trim(sys_sdEL(DerrorLambda(2),10)))   
        call output_line ('||lambda-lambda0||_(0,T] = '//trim(sys_sdEL(DerrorLambda(3),10)))   
        call output_line ('||lambda-lambda0||_(0,T) = '//trim(sys_sdEL(DerrorLambda(4),10)))   

        call output_line ('||xi-xi0||_[0,T]         = '//trim(sys_sdEL(DerrorXi(1),10)))   
        call output_line ('||xi-xi0||_[0,T)         = '//trim(sys_sdEL(DerrorXi(2),10)))   
        call output_line ('||xi-xi0||_(0,T]         = '//trim(sys_sdEL(DerrorXi(3),10)))   
        call output_line ('||xi-xi0||_(0,T)         = '//trim(sys_sdEL(DerrorXi(4),10)))   
      end if

      ! Should we calculate the functional?
      if (rpostproc%icalcFunctionalValues .ne. 0) then
        call output_lbrk()
        call optcana_nonstatFunctional (rsettings%rglobalData,rsettings%rphysicsPrimal,roptControl%rconstraints,&
            rvector,roptcontrol%rtargetFunction,roptControl%dalphaC,roptControl%dgammaC,Derror)
        call output_line ('||y-z||       = '//trim(sys_sdEL(Derror(1),10)))
        call output_line ('||u||         = '//trim(sys_sdEL(Derror(2),10)))
        call output_line ('||y(T)-z(T)|| = '//trim(sys_sdEL(Derror(3),10)))
        call output_line ('J(y,u)        = '//trim(sys_sdEL(Derror(4),10)))
      end if
    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_postprocSpaceVisOutput (rsettings,rspaceDiscr,rspaceDiscrPrimal,&
      rtimeDiscr,rvector,rrhs,ioutputUCD,sfilename)
  
!<description>
  ! For every sub-solution in the global space-time vector rvector,
  ! a visualisation file is written to disc.
  !
  ! Used for debug purposes, as this routine does not use a
  ! postprocessing structure.
!</description>

!<input>
  ! Global settings structure, passed to callback routines.
  type(t_settings_optflow), intent(inout), target :: rsettings

  ! A space-time discretisation structure defining the discretisation of
  ! rvector.
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscr

  ! A space-time discretisation structure defining the discretisation of
  ! the primal space.
  type(t_blockDiscretisation), intent(in), target :: rspaceDiscrPrimal
  
  ! Underlying time discretisation
  type(t_timeDiscretisation), intent(in) :: rtimeDiscr

  ! A space-time vector. For every timestep, a GMV is written.
  type(t_spaceTimeVector), intent(in) :: rvector

  ! A space-time vector representing the RHS.
  type(t_spaceTimeVector), intent(in) :: rrhs
  
  ! Type of output file to generate from solutions.
  ! 0=disabled
  ! 1=GMV
  ! 2=AVS
  ! 3=Paraview (VTK)
  ! 4=Matlab
  integer, intent(in) :: ioutputUCD

  ! A path + basic filename for the GMV-files. A number '.00000','.00001',...
  ! is appended for every timestep.
  character(LEN=*), intent(in) :: sfilename
!</input>

!</subroutine>

    ! local variables
    type(t_optcPostprocessing) :: rpostproc
    type(t_vectorBlock) :: rvecTemp,rrhsTemp
    integer :: istep
    real(dp) :: dtimePrimal,dtimeDual,dtstep
    
    ! Create a default postprocessing structure
    call optcpp_initpostprocessing (rpostproc,rsettings%rphysicsPrimal,CCSPACE_PRIMALDUAL,&
        rsettings%roptcBDC,rtimeDiscr,rspaceDiscr,rspaceDiscrPrimal)
    
    ! Only apply the visualisation output
    rpostproc%ioutputUCD = ioutputUCD
    rpostproc%sfilenameUCD = sfilename
    
    ! Create a temp vector for the output
    call lsysbl_createVectorBlock(rspaceDiscr,rvecTemp)
    call lsysbl_createVectorBlock(rspaceDiscr,rrhsTemp)
    
    ! Write a file for every timestep
    do istep = 1,rvector%NEQtime
      call tdiscr_getTimestep(rtimediscr,istep-1,dtimePrimal,dtstep)
      dtimeDual = dtimePrimal - (1.0_DP-rtimediscr%dtheta)*dtstep
      call sptivec_getTimestepData(rvector,istep,rvecTemp)
      call sptivec_getTimestepData(rrhs,istep,rrhsTemp)
      call optcpp_postprocessSingleSol (rpostproc,rsettings%rsettingsOptControl,rsettings,&
          istep, istep .eq. 1, istep .eq. rvector%NEQtime,&
          dtimePrimal,dtimeDual,rvecTemp,rrhsTemp,&
          dtimePrimal,dtimeDual,rvecTemp,rrhsTemp,&
          dtimePrimal,dtimeDual,rvecTemp,rrhsTemp)
    end do
    
    ! Release all created stuff
    call lsysbl_releaseVector (rrhsTemp)
    call lsysbl_releaseVector (rvecTemp)
    call optcpp_donepostprocessing (rpostproc)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine optcpp_calcControl (rphysics,rvectorSol,rconstraints,dalpha,dtimeDual,&
      ispacetimeFormulation,rvectorControl)
  
!<description>
  ! Uses the dual solution in rvectorSol to calculate the discrete 
  ! control into rvectorControl.
!</description>

!<input>
  ! Physics of the problem
  type(t_settings_physics), intent(in) :: rphysics

  ! Constraints in the optimal control problem.
  type(t_optcconstraintsSpaceTime), intent(in) :: rconstraints

  ! Regularisation parameter $\alpha$.
  real(DP), intent(in) :: dalpha
  
  ! Time corresponding to the dual solution
  real(DP), intent(in) :: dtimeDual
  
  ! Formulation of the Space-time problem.
  ! =0: Formulation for the generation of reference results from papers
  ! =1: usual formulation as specified in the DFG applicance
  ! The two formulations differ in a "-"-sign in front of the dual velocity.
  integer, intent(in) :: ispaceTimeFormulation
!</input>

!<inputoutput>
  ! Solution vector. If rvectorControl is not specified, the dual solution
  ! is replaced by the control.
  type(t_vectorBlock), intent(inout) :: rvectorSol

  ! Receives the discrete control. If not specified, the dual solution
  ! in rvectorSol is replaced by the control.
  type(t_vectorBlock), intent(inout), optional :: rvectorControl
!</inputoutput>

!</subroutine>

    type(t_optcconstraintsSpace) :: rconstrSpace

    select case (rphysics%cequation)
    case (0,1)
      ! Stokes, Navier-Stokes, 2D

      ! Copy, scale and impose constraints if necessary.
      if (present(rvectorControl)) then
        call lsyssc_copyVector (rvectorSol%RvectorBlock(4),rvectorControl%RvectorBlock(1))
        call lsyssc_copyVector (rvectorSol%RvectorBlock(5),rvectorControl%RvectorBlock(2))
        select case (ispaceTimeFormulation)
        case (0)
          call lsyssc_scaleVector (rvectorControl%RvectorBlock(1),1.0_DP/dalpha)
          call lsyssc_scaleVector (rvectorControl%RvectorBlock(2),1.0_DP/dalpha)
        case (1)
          call lsyssc_scaleVector (rvectorControl%RvectorBlock(1),-1.0_DP/dalpha)
          call lsyssc_scaleVector (rvectorControl%RvectorBlock(2),-1.0_DP/dalpha)
        end select
        
        if (rconstraints%ccontrolConstraints .ne. 0) then
          select case (rconstraints%cconstraintsType)
          case (0)
            call smva_projectControlTstepConst (rvectorControl%RvectorBlock(1),&
                rconstraints%dumin1,rconstraints%dumax1)
            call smva_projectControlTstepConst (rvectorControl%RvectorBlock(2),&
                rconstraints%dumin2,rconstraints%dumax2)
          case (1)
            ! Initialise the space constraints.
            call stlin_initSpaceConstraints (rconstraints,dtimeDual,&
                rvectorSol%p_rblockDiscr,rconstrSpace)
            
            ! Implement the constraints
            call smva_projectControlTstepVec (rvectorControl%RvectorBlock(1),&
                rconstrSpace%p_rvectorumin%RvectorBlock(1),&
                rconstrSpace%p_rvectorumax%RvectorBlock(1))
            call smva_projectControlTstepVec (rvectorControl%RvectorBlock(2),&
                rconstrSpace%p_rvectorumin%RvectorBlock(2),&
                rconstrSpace%p_rvectorumax%RvectorBlock(2))
            
            ! Done.    
            call stlin_doneSpaceConstraints (rconstrSpace)
          end select
        end if
      else
        select case (ispaceTimeFormulation)
        case (0)
          call lsyssc_scaleVector (rvectorSol%RvectorBlock(4),1.0_DP/dalpha)
          call lsyssc_scaleVector (rvectorSol%RvectorBlock(5),1.0_DP/dalpha)
        case (1)
          call lsyssc_scaleVector (rvectorSol%RvectorBlock(4),-1.0_DP/dalpha)
          call lsyssc_scaleVector (rvectorSol%RvectorBlock(5),-1.0_DP/dalpha)
        end select
        
        if (rconstraints%ccontrolConstraints .ne. 0) then
          select case (rconstraints%cconstraintsType)
          case (0)
            call smva_projectControlTstepConst (rvectorSol%RvectorBlock(4),&
                rconstraints%dumin1,rconstraints%dumax1)
            call smva_projectControlTstepConst (rvectorSol%RvectorBlock(5),&
                rconstraints%dumin2,rconstraints%dumax2)
          case (1)
            ! Initialise the space constraints.
            call stlin_initSpaceConstraints (rconstraints,dtimeDual,&
                rvectorSol%p_rblockDiscr,rconstrSpace)
            
            ! Implement the constraints
            call smva_projectControlTstepVec (rvectorSol%RvectorBlock(4),&
                rconstrSpace%p_rvectorumin%RvectorBlock(1),&
                rconstrSpace%p_rvectorumax%RvectorBlock(1))
            call smva_projectControlTstepVec (rvectorSol%RvectorBlock(5),&
                rconstrSpace%p_rvectorumin%RvectorBlock(2),&
                rconstrSpace%p_rvectorumax%RvectorBlock(2))
            
            ! Done.    
            call stlin_doneSpaceConstraints (rconstrSpace)
          end select
        end if
      end if
      
    end select

  end subroutine

!  !****************************************************************************
!
!!<subroutine>
!
!  subroutine cc_initpostprocessing (rspaceDiscr,rtimeDiscr,rpostprocessing)
!
!!<description>
!  ! Initialises the given postprocessing structure rpostprocessing
!  ! according to the main problem rproblem. The structure can then be used
!  ! to generate postprocessing data.
!!</description>
!
!!<input>
!  ! A space-time discretisation structure defining the discretisation of
!  ! rvector.
!  type(t_spatialDiscretisation), intent(in), target :: rspaceDiscr
!  
!  ! Underlying time discretisation
!  type(t_timeDiscretisation), intent(in),target :: rtimeDiscr
!
!!</input>
!
!!<output>  
!  type(t_optcPostprocessing), intent(out) :: rpostprocessing
!!</output>
!
!!</subroutine>
!
!    ! Rememnber the discretisation structures.
!    rpostprocessing%p_rspaceDiscr => rspaceDiscr
!    rpostprocessing%p_rtimeDiscr => rtimeDiscr
!
!    ! For postprocessing, we need discretisation structures in the Q0 and Q1 space,
!    ! later perhaps in the Q2 space. For this purpose, derive the corresponding
!    ! discretisation structure using the 'main' discretisation structure on the
!    ! maximum level.
!    !
!    ! For simplicity, we use only the discretisation structure of the X-velocity
!    ! to derive everything.
!    
!    p_rdiscr => rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation
!
!    ! Piecewise constant space:
!    call spdiscr_deriveSimpleDiscrSc (&
!                 rspaceDiscr%RspatialDiscr(1), &
!                 EL_Q0, CUB_G1X1, &
!                 rpostprocessing%rdiscrConstant)
!
!    ! Piecewise linear space:
!    call spdiscr_deriveSimpleDiscrSc (&
!                 rspaceDiscr%RspatialDiscr(1), &
!                 EL_Q1, CUB_G2X2, &
!                 rpostprocessing%rdiscrLinear)
!  
!    ! Piecewise quadratic space:
!    call spdiscr_deriveSimpleDiscrSc (&
!                 rspaceDiscr%RspatialDiscr(1), &
!                 EL_Q2, CUB_G3X3, &
!                 rpostprocessing%rdiscrQuadratic)
!  
!  end subroutine
!
!  !****************************************************************************
!  
!!<subroutine>
!
!  subroutine cc_clearpostprocessing (rpostprocessing)
!
!!<description>
!  ! Releases all calculated data in the given postprocessing structure
!  ! so that it can be allocated again in the calculation routines.
!  ! This routine must be called at the end of the postprocessing routines
!  ! to release the temporary memory that was allocated for the vectors
!  ! in the postprocessing structure.
!!</description>
!
!!<inputoutput>  
!  type(t_optcPostprocessing), intent(INOUT) :: rpostprocessing
!!</inputoutput>
!
!!</subroutine>
!
!    ! Release all vectors which might be allocated.
!    call lsyssc_releaseVector (rpostprocessing%rvectorVelX)
!    call lsyssc_releaseVector (rpostprocessing%rvectorVelY)
!    call lsyssc_releaseVector (rpostprocessing%rvectorPressure)
!    call lsyssc_releaseVector (rpostprocessing%rvectorPressureCells)
!    call lsyssc_releaseVector (rpostprocessing%rvectorStreamfunction)
!    call lsyssc_releaseVector (rpostprocessing%rvectorH1err)
!    call lsyssc_releaseVector (rpostprocessing%rvectorH1errCells)
!
!  end subroutine
!
!  !****************************************************************************
!  
!!<subroutine>
!
!  subroutine cc_donepostprocessing (rpostprocessing)
!
!!<description>
!  ! Releases a given problem structure. All allocated memory of this structure
!  ! is released.
!!</description>
!
!!<inputoutput>  
!  type(t_optcPostprocessing), intent(INOUT) :: rpostprocessing
!!</inputoutput>
!
!!</subroutine>
!
!    ! Release all vectors -- if there are still some allocated
!    call cc_clearpostprocessing (rpostprocessing)
!
!    ! Release the discretisation structures allocated above.
!    call spdiscr_releaseDiscr(rpostprocessing%rdiscrQuadratic)
!    call spdiscr_releaseDiscr(rpostprocessing%rdiscrLinear)
!    call spdiscr_releaseDiscr(rpostprocessing%rdiscrConstant)
!  
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_printControlFunctionalStat (rproblem,dtime,rvector)
!  
!!<description>
!  ! Calculates and prints the value of the optimal control functional J(y,u) 
!  ! in the stationary case.
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT), target :: rproblem
!!</inputoutput>
!
!!<input>
!  ! The solution vector which is to be evaluated by the postprocessing routines.
!  type(t_vectorBlock), intent(IN) :: rvector
!  
!  ! Current simulation time.
!  real(dp), intent(in) :: dtime
!!</input>
!
!!</subroutine>
!    
!    ! local variables
!    real(DP), dimension(3) :: Derror
!    real(DP) :: dalphaC
!
!    ! Initialise the collection for the assembly process with callback routines.
!    ! Basically, this stores the simulation time in the collection if the
!    ! simulation is nonstationary.
!    call cc_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)
!
!    ! Analyse the deviation from the target velocity field.
!    call parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
!                                'dalphaC',dalphaC,0.1_DP)
!    call optcana_stationaryFunctional (rvector,dalphaC,Derror,&
!        rproblem%rcollection)
!
!    call output_line ('||y-z||_L2: '//trim(sys_sdEL(Derror(1),2)))
!    call output_line ('||u||_L2  : '//trim(sys_sdEL(Derror(2),2)))
!    call output_line ('J(y,u)    : '//trim(sys_sdEL(Derror(3),2)))
!    
!    ! Release assembly-stuff from the collection
!    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
!    
!  end subroutine
!
!  !****************************************************************************
!  
!!<subroutine>
!
!  subroutine cc_createVector (rvector,cvectorType,rinitialSol)
!
!!<description>
!  ! Creates a space-time vector based on parameters.
!!</description>
!
!!<input>
!  ! Identification number that defines how to create rvector. A CCSTV_xxxx constant.
!  ! =CCSTV_ZERO: A zero space-time vector.
!  ! =CCSTV_STATIONARY: The space-time vector is created as copy of the single vector rinitialSol.
!  ! =CCSTV_READFILES: The space-time vector is created by reading in files.
!  ! =CCSTV_FORWARDSIMULATION: The space-time vector is created as forward simulation starting from an
!  ! initial solution rinitialSol.
!  integer, intent(in) :: cvectorType
!
!  ! OPTIONAL: Initial solution vector. Only primal solution. Must be at the same
!  ! spatial level and discretised with the same discretisation as the timesteps
!  ! in rvector.
!  type(t_vectorBlock), intent(in), optional :: rinitialSol
!!</input>
!
!!<inputoutput>  
!  ! Space-time vector to be created.
!  type(t_spacetimevector), intent(inout) :: rvector
!  
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables
!    integer :: istep
!    type(t_vectorBlock) :: rvectorTemp
!
!    select case (cvectorType)
!    case (CCSTV_ZERO)
!      ! Initialise by zero.
!      call sptivec_clearVector (rvector)
!      
!    case (CCSTV_STATIONARY)
!      ! Copy the initial solution to all timesteps.
!      
!    end select
!
!  end subroutine



!  !****************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE ppns2D_fbdforces (rvector,rregion,Dforces,cformulation,CcubU,CcubP)
!
!!<description>
!  ! Calculates the drag-/lift-forces acting on a part of the real
!  ! boundary for a vector rvector with the solution of the 2D
!  ! (Navier-)Stokes equation. It's assumed that
!  !  rvector%rvectorBlock(1) = X-velocity,
!  !  rvector%rvectorBlock(2) = Y-velocity,
!  !  rvector%rvectorBlock(3) = pressure,
!  ! and that X- and Y-velocity is discretised with the same 
!  ! finite element.
!  ! rregion specifies a boundary region where to calculate
!  ! the force integral. Dforces(1:NDIM2D) receives the forces in the
!  ! X- and Y-direction.
!  !
!  ! Double precision version
!!</description>
!
!!<input>
!  ! The FE solution vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!  
!  ! Boundary region where to calculate the boundary forces.
!  ! Can be created e.g. by boundary_createRegion.
!  TYPE(t_boundaryRegion), INTENT(OUT) :: rregion
!  
!  ! Type identifier that specifies the formulation of the Navier Stokes
!  ! equation that should be used for performing the integration.
!  ! ... = gradient formulation,
!  ! ... = deformation formulation
!  INTEGER, INTENT(IN) :: cformulation
!  
!  ! OPTIONAL: Array with cubature formula identifiers when integrating
!  ! the velocity. If specified, there must be a cubature formula identifier
!  ! for each of the element distributions specified in the discretisation
!  ! structure that is associated to rvector.
!  ! If not specified, the standard cubature formula in the discretisation
!  ! structure for integrating linear forms is used.
!  INTEGER, DIMENSION(:), INTENT(IN) :: CcubU
!
!  ! OPTIONAL: Array with cubature formula identifiers when integrating
!  ! the pressure. If specified, there must be a cubature formula identifier
!  ! for each of the element distributions specified in the discretisation
!  ! structure that is associated to rvector.
!  ! If not specified, the standard cubature formula in the discretisation
!  ! structure for integrating linear forms is used.
!  INTEGER, DIMENSION(:), INTENT(IN) :: CcubP
!!</input>
!
!!<output>
!  ! Array receiving the forces acting on the boundary specified by rregion.
!  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
!  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dforces
!!</output>
!
!!</subroutine>
!
!  ! Array to tell the element which derivatives to calculate.
!  LOGICAL, DIMENSION(EL_MAXNDER) :: BderU, BderP
!  
!  ! Cubature point coordinates on the reference element.
!  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: DxiU, DxiP
!
!  ! For every cubature point on the reference element,
!  ! the corresponding cubature weight
!  REAL(DP), DIMENSION(CUB_MAXCUBP) :: DomegaU, DomegaP
!  
!  ! number of cubature points on the reference element
!  INTEGER :: ncubpU, ncubpP
!
!  ! Pointer to the vector entries
!  REAL(DP), DIMENSION(:), POINTER :: p_DdataUX, p_DdataUY, p_DdataP
!
!  ! An allocateable array accepting the DOF's of a set of elements.
!  integer, DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrialU, IdofsTrialP
!  
!  ! Allocateable arrays for the values of the basis functions.
!  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTrialU,DbasTrialP
!  
!  ! Number of entries in the vector - for quicker access
!  INTEGER(I32) :: neqU, neqP
!
!  ! Arrays for saving Jacobian determinants and matrices
!  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
!  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Djac
!
!  ! Pointer to KVERT of the triangulation
!  INTEGER(I32), DIMENSION(:,:), POINTER :: p_IverticesAtElement
!  
!  ! Pointer to DCORVG of the triangulation
!  REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
!  
!  ! Current element distribution
!  TYPE(t_elementDistribution), POINTER :: p_elementDistribution
!
!  ! Number of elements in a block. Normally =BILF_NELEMSIM,
!  ! except if there are less elements in the discretisation.
!  INTEGER :: nelementsPerBlock
!  
!
!  BderU = .FALSE.
!  BderP = .FALSE.
!
!  ! Get the vector data
!  neqU = rvector%RvectorBlock(1)%NEQ
!  neqP = rvector%RvectorBlock(3)%NEQ
!  
!  IF (rvector%cdataType .NE. ST_DOUBLE) THEN
!    PRINT *,'ppns2D_fbdforces: Unsupported vector precision.'
!    STOP
!  END IF
!
!  ! We support only uniform and conformal discretisation structures.
!  IF (.NOT. ASSOCIATED(rvector%p_rblockDiscretisation)) THEN
!    PRINT *,'ppns2D_fbdforces: No discretisation structure!'
!    STOP
!  END IF
!
!  IF ((rvector%p_rblockDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) .AND. &
!      (rvector%p_rblockDiscretisation%ccomplexity .NE. SPDISC_CONFORMAL)) THEN
!    PRINT *,'ppns2D_fbdforces: Discretisation too complex!'
!    STOP
!  END IF
!  
!  ! Get pointers to the subvectors from the block vector
!  CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataUX)
!  CALL lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataUY)
!  CALL lsyssc_getbase_double (rvector%RvectorBlock(3),p_DdataP)
!
!  ! Now loop over the different element distributions in the discretisation.
!  ! Each element distribution has a different element combination.
!  ! E.g. distribution 1 may describe P1~/P1~/P0 triangular elements, while
!  ! distrbution 2 may describe Q1~/Q1~/Q0.
!
!  !DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
!
!
!  !END DO ! icurrentElementDistr
!  
!  Dforces = 0.0_DP
!  
!  END SUBROUTINE
!
!
!  !****************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE ppns2D_bdforces (rvector,rregion,Dforces,cformulation,CcubU,CcubP)
!
!!<description>
!  ! Calculates the drag-/lift-forces acting on a part of the real
!  ! boundary for a vector rvector with the solution of the 2D
!  ! (Navier-)Stokes equation. It's assumed that
!  !  rvector%rvectorBlock(1) = X-velocity,
!  !  rvector%rvectorBlock(2) = Y-velocity,
!  !  rvector%rvectorBlock(3) = pressure,
!  ! and that X- and Y-velocity is discretised with the same 
!  ! finite element.
!  ! rregion specifies a boundary region where to calculate
!  ! the force integral. Dforces(1:NDIM2D) receives the forces in the
!  ! X- and Y-direction.
!  !
!  ! Double precision version
!!</description>
!
!!<input>
!  ! The FE solution vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!  
!  ! Boundary region where to calculate the boundary forces.
!  ! Can be created e.g. by boundary_createRegion.
!  TYPE(t_boundaryRegion), INTENT(OUT) :: rregion
!  
!  ! Type identifier that specifies the formulation of the Navier Stokes
!  ! equation that should be used for performing the integration.
!  ! ... = gradient formulation,
!  ! ... = deformation formulation
!  INTEGER, INTENT(IN) :: cformulation
!  
!  ! 1D Cubature formula identifier to use for the line integration.
!  ! One of the CUB_xxxx_1D constants in the cubature.f90.
!  INTEGER, DIMENSION(:), INTENT(IN) :: CcubU
!!</input>
!
!!<output>
!  ! Array receiving the forces acting on the boundary specified by rregion.
!  ! Note: These are the drag-/lift-FORCES, not the coefficients!!!
!  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dforces
!!</output>
!
!!</subroutine>
!
!  ! Spatial discretisation structure of velocity and pressure
!  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscrU, p_rdiscrP
!  
!  ! Element distribution of velocity and pressure
!  TYPE(t_elementDistribution), POINTER :: p_elemDistrU,p_elemDistrP
!
!
!  ! Array to tell the element which derivatives to calculate.
!  LOGICAL, DIMENSION(EL_MAXNDER) :: BderU, BderP
!  
!  ! Cubature point coordinates on the reference element.
!  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi1D, DXi2D
!
!  ! For every cubature point on the reference element,
!  ! the corresponding cubature weight
!  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
!  
!  ! number of cubature points on the reference element
!  INTEGER :: ncubp
!
!  ! Pointer to the vector entries
!  REAL(DP), DIMENSION(:), POINTER :: p_DdataUX, p_DdataUY, p_DdataP
!
!  ! An allocateable array accepting the DOF's of a set of elements.
!  integer, DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrialU, IdofsTrialP
!  
!  ! Allocateable arrays for the values of the basis functions.
!  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTrialU,DbasTrialP
!  
!  ! Number of entries in the vector - for quicker access
!  INTEGER(I32) :: neqU, neqP
!
!  ! Arrays for saving Jacobian determinants and matrices
!  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
!  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Djac
!
!  ! Pointer to KVERT of the triangulation
!  INTEGER(I32), DIMENSION(:,:), POINTER :: p_IverticesAtElement
!  
!  ! Pointer to DCORVG of the triangulation
!  REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
!  
!
!  ! Number of elements in a block. Normally =BILF_NELEMSIM,
!  ! except if there are less elements in the discretisation.
!  INTEGER :: nelementsPerBlock
!  
!  BderU = .FALSE.
!  BderP = .FALSE.
!
!  ! Get the vector data
!  neqU = rvector%RvectorBlock(1)%NEQ
!  neqP = rvector%RvectorBlock(3)%NEQ
!  
!  IF (rvector%cdataType .NE. ST_DOUBLE) THEN
!    PRINT *,'ppns2D_bdforces: Unsupported vector precision.'
!    STOP
!  END IF
!
!  ! We support only uniform and conformal discretisation structures.
!  IF (.NOT. ASSOCIATED(rvector%p_rblockDiscretisation)) THEN
!    PRINT *,'ppns2D_bdforces: No discretisation structure!'
!    STOP
!  END IF
!
!  IF ((rvector%p_rblockDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) .AND. &
!      (rvector%p_rblockDiscretisation%ccomplexity .NE. SPDISC_CONFORMAL)) THEN
!    PRINT *,'ppns2D_bdforces: Discretisation too complex!'
!    STOP
!  END IF
!  
!  ! Get pointers to the subvectors from the block vector
!  CALL lsyssc_getbase_double (rvector%RvectorBlock(1),p_DdataUX)
!  CALL lsyssc_getbase_double (rvector%RvectorBlock(2),p_DdataUY)
!  CALL lsyssc_getbase_double (rvector%RvectorBlock(3),p_DdataP)
!  
!  ! Get pointers to the spatial discretisation structures of the
!  ! velocity and pressure
!  p_rdiscrU => rvector%RvectorBlock(1)%p_spatialDiscretisation
!  p_rdiscrP => rvector%RvectorBlock(3)%p_spatialDiscretisation
!
!  ! Now loop over the different element distributions in the discretisation.
!  ! Each element distribution has a different element combination.
!  ! E.g. distribution 1 may describe P1~/P1~/P0 triangular elements, while
!  ! distrbution 2 may describe Q1~/Q1~/Q0.
!  !
!  ! For every element distribution, we have to loop through all elements
!  ! along the boundary component, calculate their contribution to the
!  ! integral and add everything together.
!
!  DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
!
!    ! Get the element distribution structure of that FE-space comination
!    p_relemDistrU => p_rdiscrU%RelementDistr(RelementDistr)
!    p_relemDistrP => p_rdiscrP%RelementDistr(RelementDistr)
!    
!    ! By that we know
!    ! - what's the element
!    ! - which transformation from the reference element to use
!    
!    
!    
!    
!
!  END DO ! icurrentElementDistr
!  
!  END SUBROUTINE

end module

!##############################################################################
!# ****************************************************************************
!# <name> ccpostprocessing </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains postprocessing routines for the CC2D solver.
!#
!# The following routines can be found here:
!#
!# 1.) cc_initpostprocessing
!#     -> Initialise the postprocessing.
!#
!# 2.) cc_donepostprocessing
!#     -> Clean up the postprocessing.
!#
!# 3.) cc_postprocessingStationary
!#     -> Postprocessing of a solution in a stationary solution process.
!#        Evaluate the solution of the stationary solver, write GMV-file.
!#
!# 4.) cc_postprocessingNonstat
!#     -> Postprocessing of a solution in a nonstationary solution process.
!#        Evaluate the solution of the nonstationary solver, write GMV-file.
!#
!# Auxiliary routines:
!#
!# 1.) cc_errorAnalysis
!#     -> Perform error analysis (comparison to analytic solution;
!#        as configured in the DAT files).
!#
!# 2.) cc_calculateBodyForces
!#     -> Calculate the body forces on an object.
!#
!# 3.) cc_calculateDivergence
!#     -> Calculate the divergence.
!#
!# 4.) cc_writeUCD
!#     -> Write UCD (AVS, GMV,...) output.
!#
!# 5.) cc_writeFilm
!#     -> Write film output (=raw solution vectors).
!#
!# 5.) cc_evaluatePoints
!#     -> Evaluates the solution in a number of points.
!#
!# </purpose>
!##############################################################################

module ccpostprocessing

  use fsystem
  use storage
  use linearsolver
  use boundary
  use linearalgebra
  use cubature
  use matrixfilters
  use vectorfilters
  use discretebc
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use bilinearformevaluation
  use linearformevaluation
  use statistics
  use element
  use multilevelprojection
  use vectorio
  use io
  use fparser

  use collection
  use convection

  use ucd

  use pprocnavierstokes
  use pprocerror

  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop
  use feevaluation2

  use ccboundaryconditionparser
  use ccgeneraldiscretisation
  use ccmatvecassembly
  use ccbasic
  use cccallback

  use analyticprojection

  implicit none

!<types>

!<typeblock>

  ! Postprocessing structure. Whenever the postprocessing calculation routines are
  ! called, they fill this structure using the current solution vector (and other
  ! information if necessary). The information in this structure can then be used
  ! for GMV output e.g.
  type t_c2d2postprocessing

    ! A discretisation structure that describes a piecewise constant discretisation
    ! (usually P0 or Q0).
    type(t_spatialDiscretisation) :: rdiscrConstant

    ! A discretisation structure that describes a piecewise linear discretisation
    ! (usually P1 or Q1).
    type(t_spatialDiscretisation) :: rdiscrLinear

    ! A discretisation structure that describes a piecewise quadratic discretisation
    ! (usually P2 or Q2).
    type(t_spatialDiscretisation) :: rdiscrQuadratic

    ! Whether nonstationary postprocessing should be used or not.
    logical              :: bnonstationaryPostprocessing

    ! Next file extension for UCD output file.
    integer              :: inextFileSuffixUCD = 0

    ! Next file extension for Film output file.
    integer              :: inextFileSuffixFilm = 0

    ! A vector that describes the X-velocity field in the vertices
    type(t_vectorScalar) :: rvectorVelX

    ! A vector that describes the Y-velocity field in the vertices
    type(t_vectorScalar) :: rvectorVelY

    ! A vector that describes the pressure field in the vertices
    type(t_vectorScalar) :: rvectorPressure

    ! A vector that describes the pressure field in the cells
    type(t_vectorScalar) :: rvectorPressureCells

    ! A vector that describes the streamfunction
    type(t_vectorScalar) :: rvectorStreamfunction

    ! A vector that describes the H1-error of the velocity field in the vertices
    type(t_vectorScalar) :: rvectorH1err

    ! A vector that describes the H1-error of the pressure in the cells
    type(t_vectorScalar) :: rvectorH1errCells

    ! (Possibly interpolated) solution from time step t^{n-2}, for time-space
    ! discretisation error analysis
    type(t_vectorBlock) :: roldestSolution

    ! (Possibly interpolated) solution from time step t^{n-1}, for time-space
    ! discretisation error analysis
    type(t_vectorBlock) :: rolderSolution

    ! (Possibly interpolated) solution from time step t^{n}, for time-space
    ! discretisation error analysis
    type(t_vectorBlock) :: roldSolution

    ! Time step scheme identifier. One of the TSCHM_xxxx-constants.
    ! Usually TSCHM_ONESTEP for one step schemes
    integer :: ctimestepType = -1

    ! Whether or not to compute L2 errors and if yes, which callback to evaluate analytic
    ! reference solution
    integer :: icalcL2

    ! Whether or not to compute H1 errors and if yes, which callback to evaluate analytic
    ! reference solution
    integer :: icalcH1

    ! Whether or not to compute time plus space discretisation errors
    integer :: icalcTimeSpaceDiscrErrors

    ! How to measure the time plus space discretisation error
    integer :: itimeSpaceDiscrErrorMethod

    ! Parser object that encapsules error definitions for the L2 error.
    type(t_fparser) :: rrefFunctionL2

    ! Parser object that encapsules error definitions for the H1 error.
    type(t_fparser) :: rrefFunctionH1

  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc_postprocessingStationary (rproblem,rvector,rpostprocessing)

!<description>
  ! Postprocessing of solutions of stationary simulations.
  ! Writes the solution into a GMV file, calculates forces,...
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! Postprocessing structure.
  type(t_c2d2postprocessing), intent(inout) :: rpostprocessing
!</inputoutput>

!<input>
  ! The solution vector which is to be evaluated by the postprocessing routines.
  type(t_vectorBlock), intent(in) :: rvector
!</input>

!</subroutine>

    ! local variables
    type(t_timer) :: rtimer

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    ! Calculate body forces.
    call cc_calculateBodyForces (rvector,0.0_DP,rproblem)

    ! Calculate point values
    call cc_evaluatePoints (rvector,0.0_DP,rproblem)

    ! Calculate flux values
    call cc_evaluateFlux (rvector,0.0_DP,rproblem)

    ! Calculate the divergence
    call cc_calculateDivergence (rvector,rproblem)

    ! Error analysis, comparison to reference function.
    rpostprocessing%ctimestepType = -1
    call cc_errorAnalysis (rvector,0.0_DP,rpostprocessing,rproblem)

    ! Write the UCD export file (GMV, AVS,...) as configured in the DAT file.
    call cc_writeUCD (rpostprocessing, rvector, rproblem)

    ! Gather statistics
    call stat_stopTimer(rtimer)
    rproblem%rstatistics%dtimePostprocessing = &
      rproblem%rstatistics%dtimePostprocessing + rtimer%delapsedReal

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_postprocessingNonstat (rproblem,rvectorPrev,&
      dtimePrev, rvector, dtime, rvectorInt, dtimeInt, istep, rpostprocessing, &
      bisIntermediateStage)

!<description>
  ! Postprocessing of solutions of transient simulations.
  ! Writes the solution into a GMV file, calculates forces,...
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! Postprocessing structure. Defines what to do with solution vectors.
  type(t_c2d2postprocessing), intent(inout) :: rpostprocessing
!</inputoutput>

!<input>
  ! Solution vector of the previous timestep. Coincides with rvector
  ! if there is no previous timestep.
  type(t_vectorBlock), intent(in) :: rvectorPrev

  ! Time of the previous timestep. Coincides with dtime
  ! if there is no previous timestep.
  real(DP), intent(in) :: dtimePrev

  ! Solution vector of the current timestep.
  type(t_vectorBlock), intent(in) :: rvector

  ! Time of the current timestep.
  real(DP), intent(in) :: dtime

  ! Interpolated (in time) solution vector. Velocity and pressure
  ! represent the same point in time.
  type(t_vectorBlock), intent(in) :: rvectorInt

  ! Time of the interpolated solution vector.
  real(DP), intent(in) :: dtimeInt

  ! Number of the timestep. =0: initial solution
  integer, intent(in) :: istep

  ! Whether or not visualisation output should get considered (it should not in
  ! intermediate stages of a multi-stage time stepping scheme (like FS or DIRK)
  logical, intent(in) :: bisIntermediateStage
!</input>

!</subroutine>

    ! local variables
    type(t_timer) :: rtimer
    real(DP) :: dminTime, dmaxTime, dtimeDifferenceUCD
    real(DP) :: dtimeDifferenceFilm, dpptime
    integer :: itime1,itime2,iinterpolateSolutionUCD,iinterpolateSolutionFilm
    integer :: ipostprocTimeInterpSolution
    integer :: iwriteSolDeltaSteps
    type(t_vectorBlock) :: rintVector
    real(DP) :: dweight,dwriteSolDeltaTime

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    ! Whether to apply postprocessing to rvector or rvectorInt
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "ipostprocTimeInterpSolution", ipostprocTimeInterpSolution, 1)

    if (.not. bisIntermediateStage) then
      ! Think about writing out the solution...
      call parlst_getvalue_double (rproblem%rparamList, "CC-DISCRETISATION", &
           "dwriteSolDeltaTime", dwriteSolDeltaTime, 0.0_DP)
      call parlst_getvalue_int (rproblem%rparamList, "CC-DISCRETISATION", &
           "iwriteSolDeltaSteps", iwriteSolDeltaSteps, 1)
      if (iwriteSolDeltaSteps .lt. 1) iwriteSolDeltaSteps = 1

      ! Figure out if we have to write the solution. This is the case if
      ! 1.) Previous and current time is the same (= first call) or
      ! 2.) The solution "crossed the next timestep".

      if ((dwriteSolDeltaTime .gt. 0.0_DP) .and. (dtimePrev .ne. dtime)) then
        itime1 = int((dtimePrev-rproblem%rtimedependence%dtimeInit)/dwriteSolDeltaTime)
        itime2 = int((dtime-rproblem%rtimedependence%dtimeInit)/dwriteSolDeltaTime)
      else
        itime1 = 0
        itime2 = 1
      end if

      if (itime1 .ne. itime2) then
        ! Write the raw solution
        call cc_writeSolution (rproblem,rvector,dtime)
      else
        ! The second option is that the timestep matches.
        if (mod(istep,iwriteSolDeltaSteps) .eq. 0) then
          ! Write the raw solution
          call cc_writeSolution (rproblem,rvector,dtime)
        end if
      end if
    end if

    if (ipostprocTimeInterpSolution .ne. 0) then
      if (.not. bisIntermediateStage) then
        ! Calculate body forces.
        call cc_calculateBodyForces (rvectorInt,dtimeInt,rproblem)

        ! Calculate point values
        call cc_evaluatePoints (rvectorInt,dtimeInt,rproblem)

        ! Calculate flux values
        call cc_evaluateFlux (rvector,dtimeInt,rproblem)
      end if

      ! Calculate the divergence
      call cc_calculateDivergence (rvectorInt,rproblem)

      ! Error analysis, comparison to reference function.
      call cc_errorAnalysis (rvectorInt,dtimeInt,rpostprocessing,rproblem)
    else
      if (.not. bisIntermediateStage) then
        ! Calculate body forces.
        call cc_calculateBodyForces (rvector,dtime,rproblem)

        ! Calculate point values
        call cc_evaluatePoints (rvector,dtime,rproblem)

        ! Calculate flux values
        call cc_evaluateFlux (rvector,dtime,rproblem)

        ! Calculate the divergence
        call cc_calculateDivergence (rvector,rproblem)
      end if

      ! Error analysis, comparison to reference function.
      call cc_errorAnalysis (rvectorInt,dtimeInt,rpostprocessing,rproblem)

    end if

    if (.not. bisIntermediateStage) then
      ! Write the UCD export file (GMV, AVS,...) as configured in the DAT file.
      !
      ! In a nonstationary simulation, first check if we are allowed
      ! to write something.

      call parlst_getvalue_double (rproblem%rparamList, "CC-POSTPROCESSING", &
           "DMINTIMEUCD", dminTime, -1.E100_DP)
      call parlst_getvalue_double (rproblem%rparamList, "CC-POSTPROCESSING", &
           "DMAXTIMEUCD", dmaxTime, 1.E100_DP)
      call parlst_getvalue_double (rproblem%rparamList, "CC-POSTPROCESSING", &
           "DTIMEDIFFERENCEUCD", dtimeDifferenceUCD, 0.0_DP)

      if ((dtime .ge. dminTime-1000.0*SYS_EPSREAL_DP) .and. &
           (dtime .le. dmaxTime+1000.0*SYS_EPSREAL_DP)) then

        ! Figure out if we have to write the solution. This is the case if
        ! 1.) Precious and current time is the same or
        ! 2.) The solution crossed the next ucd timestep.

        if ((dtimeDifferenceUCD .gt. 0.0_DP) .and. (dtimePrev .ne. dtime)) then
          itime1 = int((dtimePrev-rproblem%rtimedependence%dtimeInit)/dtimeDifferenceUCD)
          itime2 = int((dtime-rproblem%rtimedependence%dtimeInit)/dtimeDifferenceUCD)
        else
          itime1 = 0
          itime2 = 1
        end if

        if (itime1 .ne. itime2) then
          ! Probably we have to interpolate the solution to the point dtime in time.
          call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
               "IINTERPOLATESOLUTIONUCD", iinterpolateSolutionUCD,1)
          if ((iinterpolateSolutionUCD .eq. 0) .or. (dtimeDifferenceUCD .eq. 0.0_DP) &
               .or. (dtimePrev .eq. dtime)) then
            ! No interpolation
            call cc_writeUCD (rpostprocessing, rvector, rproblem, dtime)
          else
            ! Interpolate and write out the interpolated solution.
            call lsysbl_copyVector (rvectorPrev,rintVector)
            dpptime = real(itime2, DP) * dtimeDifferenceUCD + &
                      rproblem%rtimedependence%dtimeInit
            dweight = (dpptime-dtimePrev) / (dtime-dtimePrev)
            call lsysbl_vectorLinearComb (rvector,rintVector,dweight,1.0_DP-dweight)
            call cc_writeUCD (rpostprocessing, rintVector, rproblem, dpptime)
            call lsysbl_releaseVector (rintVector)
          end if
        end if
      end if

      ! Write film output (raw data vectors)
      !
      ! First check if we are allowed to write something.
      call parlst_getvalue_double (rproblem%rparamList, "CC-POSTPROCESSING", &
           "DMINTIMEFILM", dminTime, -1.E100_DP)
      call parlst_getvalue_double (rproblem%rparamList, "CC-POSTPROCESSING", &
           "DMAXTIMEFILM", dmaxTime, 1.E100_DP)
      call parlst_getvalue_double (rproblem%rparamList, "CC-POSTPROCESSING", &
           "DTIMEDIFFERENCEFILM", dtimeDifferenceFilm, 0.0_DP)

      if ((dtime .ge. dminTime-1000.0*SYS_EPSREAL_DP) .and. &
           (dtime .le. dmaxTime+1000.0*SYS_EPSREAL_DP)) then

        ! Figure out if we have to write the solution. This is the case if
        ! 1.) Precious and current time is the same ot
        ! 2.) The solution crossed the next ucd timestep.

        if ((dtimeDifferenceFilm .gt. 0.0_DP) .and. (dtimePrev .ne. dtime)) then
          itime1 = int((dtimePrev-rproblem%rtimedependence%dtimeInit)/dtimeDifferenceFilm)
          itime2 = int((dtime-rproblem%rtimedependence%dtimeInit)/dtimeDifferenceFilm)
        else
          itime1 = 0
          itime2 = 1
        end if

        if (itime1 .ne. itime2) then
          ! Probably we have to interpolate the solution to the point dtime in time.
          call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
               "IINTERPOLATESOLUTIONFILM", iinterpolateSolutionFilm,1)
          if ((iinterpolateSolutionFilm .eq. 0) .or. (dtimeDifferenceFilm .eq. 0.0_DP) &
               .or. (dtimePrev .eq. dtime)) then
            ! No interpolation
            call cc_writeFilm (rpostprocessing, rvector, rproblem, dtime)
          else
            ! Interpolate and write out the interpolated solution.
            call lsysbl_copyVector (rvectorPrev,rintVector)
            dpptime = real(itime2,DP) * dtimeDifferenceFilm + &
                      rproblem%rtimedependence%dtimeInit
            dweight = (dpptime-dtimePrev) / (dtime-dtimePrev)
            call lsysbl_vectorLinearComb (rvector,rintVector,dweight,1.0_DP-dweight)
            call cc_writeFilm (rpostprocessing, rintVector, rproblem, dpptime)
            call lsysbl_releaseVector (rintVector)
          end if
        end if
      end if
    end if

    ! Gather statistics
    call stat_stopTimer(rtimer)
    rproblem%rstatistics%dtimePostprocessing = &
      rproblem%rstatistics%dtimePostprocessing + rtimer%delapsedReal

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine fcalc_error (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)

  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration

!<description>
  ! This routine is called during the calculation of L2/H1 errors.
  ! It evaluates errors being given as expressions.
!</description>

!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation

  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements

  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement

  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! A pointer to a collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection

!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>

!</subroutine>

    ! local variables
    integer :: ctype,icomponent
    integer :: ielement,ipoint, iparsercomp
    real(DP), dimension(size(EXPRVARIABLES)) :: Rval

    ! Type. 0=L2 error, 1=H1 error
    ctype = rcollection%IquickAccess(1)

    ! Underlying component of the expression to evaluate
    icomponent = rcollection%IquickAccess(2)

    ! Fill the evaluation structure
    Rval(:) = 0.0_DP
    Rval(7:11) = rcollection%DquickAccess(1:5)

    select case (ctype)

    ! ---------------------------------
    ! L2 error
    ! ---------------------------------
    case (0)
      do ielement = 1,nelements
        do ipoint = 1,npointsPerElement

          ! Get the point coordinates.
          Rval(1) = Dpoints(1,ipoint,ielement)
          Rval(2) = Dpoints(2,ipoint,ielement)

          ! Evaluate
          call fparser_evalFunction (rcollection%p_rfparserQuickAccess1, &
              icomponent, Rval, Dvalues(ipoint,ielement))

        end do
      end do

    ! ---------------------------------
    ! H1 error
    ! ---------------------------------
    case (1)

      ! ---------------------------------
      ! X or Y derivative?
      ! ---------------------------------
      select case (cderivative)

      ! Component in the parser array.
      ! 1,3,5 = X-derivative, 2,4,6 = Y-derivative
      case (DER_DERIV2D_X)
        iparsercomp = NDIM2D*icomponent-1

      case (DER_DERIV2D_Y)
        iparsercomp = NDIM2D*icomponent

      end select

      do ielement = 1,nelements
        do ipoint = 1,npointsPerElement

          ! Get the point coordinates.
          Rval(1) = Dpoints(1,ipoint,ielement)
          Rval(2) = Dpoints(2,ipoint,ielement)

          ! Evaluate
          call fparser_evalFunction (rcollection%p_rfparserQuickAccess1, &
              iparsercomp, Rval, Dvalues(ipoint,ielement))

        end do
      end do

    end select
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_errorAnalysis (rsolution,dtime,rpostprocessing,rproblem)

!<description>
  ! Performs error analysis on a given solution rsolution as specified
  ! in the .DAT file.
  ! The result of the error analysis is written to the standard output.
!</description>

!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(in) :: rsolution

  ! Solution time. =0 for stationary simulations.
  real(DP), intent(in) :: dtime

  ! Postprocessing structure. Defines what to do with solution vectors.
  type(t_c2d2postprocessing), intent(inout), target :: rpostprocessing
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(3) :: Derr
    real(DP) :: derrorVel, derrorP, denergy
    integer :: icalcEnergy
    integer :: iwriteErrorAnalysisL2,iwriteErrorAnalysisH1,iwriteKineticEnergy
    character(len=SYS_STRLEN) :: sfilenameErrorAnalysisL2
    character(len=SYS_STRLEN) :: sfilenameErrorAnalysisH1
    character(len=SYS_STRLEN) :: sfilenameKineticEnergy
    character(len=SYS_STRLEN) :: stemp
    integer :: iunit,cflag,icubError
    logical :: bfileExists
    real(DP) :: dtimebackup
    type(t_scalarCubatureInfo) :: rcubatureInfoUV,rcubatureInfoP
    type(t_collection) :: rlocalCollection
    integer :: iglobaltimestep, cignoreTimeStep
    integer :: ctypeInitialSolution

    ! all subsequent variables for calculation of accumulated time plus space
    ! discretisation error ...

    ! ... with trapezoidal rule
    real(DP) :: derrorL2Vel_old = 0.0_DP, derrorL2P_old = 0.0_DP
    real(DP) :: derrorH1Vel_old = 0.0_DP, derrorH1P_old = 0.0_DP

    ! ... with quadratic/cubic Lagrange time-interpolation and 2- or 3-Point-Gauss
    ! quadrature integration
    real(DP), save :: doldestTime = 0.0_DP, dolderTime = 0.0_DP, doldTime = 0.0_DP
    integer :: ipostprocTimeInterpSolution
    real(DP) :: dtime_GaussPt1, dtime_GaussPt2, dtime_GaussPt3
                                         ! Gauss points for 2- or 3-Point-Gauss quadrature
    real(DP), dimension(3) :: Derr2, Derr3
    type(t_vectorBlock) :: rvector_sol_GaussPt1, rvector_sol_GaussPt2
    type(t_vectorBlock) :: rvector_sol_GaussPt3, rvector_auxSum

    !      \int_0^T || u_h(t) - u_ref(t) ||^2_L2 dt
    ! evaluated with trapezoidal rule or 2- or 3-point Gauss quadrature rule (note that
    ! trailing square!)
    real(DP), save :: derrorL2VelTimeSpace = 0.0_DP
    real(DP), save :: derrorL2PTimeSpace = 0.0_DP     ! analog
    real(DP), save :: derrorH1VelTimeSpace = 0.0_DP   ! analog
    real(DP), save :: derrorH1PTimeSpace = 0.0_DP     ! analog


    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "ICALCKINETICENERGY", icalcEnergy, 1)

    call parlst_getvalue_string (rproblem%rparamList,"CC-POSTPROCESSING",&
                                 "scubError",stemp,"")
    if (stemp .eq. "") then
      call parlst_getvalue_int (rproblem%rparamList,"CC-POSTPROCESSING",&
                                "icubError",icubError,int(CUB_GEN_AUTO))
    else
      icubError = cub_igetID(stemp)
    end if

    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "IWRITEERRORANALYSISL2", iwriteErrorAnalysisL2, 0)
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "IWRITEERRORANALYSISH1", iwriteErrorAnalysisH1, 0)
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "IWRITEKINETICENERGY", iwriteKineticEnergy, 0)

    call parlst_getvalue_string (rproblem%rparamList, "CC-POSTPROCESSING", &
        "SFILENAMEERRORANALYSISL2", sfilenameErrorAnalysisL2, "", bdequote=.true.)
    if (sfilenameErrorAnalysisL2 .eq. "") &
      iwriteErrorAnalysisL2 = 0

    call parlst_getvalue_string (rproblem%rparamList, "CC-POSTPROCESSING", &
        "SFILENAMEERRORANALYSISH1", sfilenameErrorAnalysisH1, "", bdequote=.true.)
    if (sfilenameErrorAnalysisH1 .eq. "") &
      iwriteErrorAnalysisH1 = 0

    call parlst_getvalue_string (rproblem%rparamList, "CC-POSTPROCESSING", &
        "SFILENAMEKINETICENERGY", sfilenameKineticEnergy, "", bdequote=.true.)
    if (sfilenameKineticEnergy .eq. "") &
      iwriteKineticEnergy = 0

    ! Create an cubature info structure which contains our cubature rule
    call spdiscr_createDefCubStructure(&
        rsolution%RvectorBlock(1)%p_rspatialDiscr,rcubatureInfoUV,int(icubError,I32))
    call spdiscr_createDefCubStructure(&
        rsolution%RvectorBlock(3)%p_rspatialDiscr,rcubatureInfoP,int(icubError,I32))

    ! Time plus space discretisation error analysis
    if (rpostprocessing%icalcTimeSpaceDiscrErrors .ne. 0) then
      if (rpostprocessing%icalcL2 .ne. 0 .or. &
          rpostprocessing%icalcH1 .ne. 0) then
        ! Whether to apply postprocessing to rvector or rvectorInt
        call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
             "IPOSTPROCTIMEINTERPSOLUTION", ipostprocTimeInterpSolution, 1)

        if (ipostprocTimeInterpSolution .ne. 1) then
          call output_line ("ERROR: Can not gauge time error if " // &
               "[CC-POSTPROCESSING] ipostprocTimeInterpSolution <> 1")
          stop
        end if
      end if

      ! For one-step schemes holds
      !     global timestep = timestep.
      ! For Fractional-Step-Theta holds
      !     global timestep = timestep/3, local timestep=timestep%3.
      !
      ! Note: The intermediate time steps of the Fractional-Step-Theta method *must* be
      !       taken into account *as well* in order to correctly measure the time plus
      !       space discretisation error. Merely using the solution at the end of macro
      !       time steps (i.e. every third time step`s solution) leads not only for
      !       Fractional-Step-Theta, but also for simpler time stepping schemes like
      !       Backward Euler and Crank-Nicolson to an overestimation of the error
      !       reduction rates! That is no theoretical argument yet, it has shown in
      !       practice for 3 different analytic solution pairs for the transient Stokes
      !       problem.
      if (rproblem%itimedependence .ne. 0) then
        iglobaltimestep = rproblem%rtimedependence%itimestep
        cignoreTimeStep = 0
        if (rpostprocessing%ctimestepType .eq. TSCHM_DIRK34La .or. &
            rpostprocessing%ctimestepType .eq. TSCHM_DIRK34Lb .or. &
            rpostprocessing%ctimestepType .eq. TSCHM_DIRK44L  .or. &
            rpostprocessing%ctimestepType .eq. TSCHM_DIRK54L) then
          ! The velocity and pressure solution in the intermediate stages of the time
          ! stepping schemes DIRK34L (variant a and b) and DIRK44L are typically wrong.
          !
          ! Note: Divisor 3 = nsubsteps - 1, given that all these schemes have 4 stages
          !       with an explicit first stage such that 3 time steps are sufficient for a
          !       complete sweep.
          iglobaltimestep = int(rproblem%rtimedependence%itimestep / 3)
          select case (mod(rproblem%rtimedependence%itimestep,3))
          case (0)
            cignoreTimeStep = 0
          case (1,2)
            cignoreTimeStep = 1
          end select
        end if
        if (rpostprocessing%ctimestepType .eq. TSCHM_SDIRK2) then
          ! The velocity and pressure solution in the intermediate stages of the time
          ! stepping scheme SDIRK2 are typically wrong.
          !
          ! Note: Divisor 4 = nsubsteps, given that this time stepping scheme has 4 stages
          !       with an implicit first stage such that really 4 time steps are needed
          !       for a complete sweep.
          iglobaltimestep = int(rproblem%rtimedependence%itimestep / 4)
          select case (mod(rproblem%rtimedependence%itimestep,4))
          case (0)
            cignoreTimeStep = 0
          case (1,2,3,4)
            cignoreTimeStep = 1
          end select
        end if
        if (rpostprocessing%ctimestepType .eq. TSCHM_SDIRK3PR) then
          ! The velocity and pressure solution in the intermediate stages of the time
          ! stepping scheme SDIRK3PR is typically wrong as does show a manual application
          ! of the schemeto the nonlinear Stokes problem with homogeneous right hand side
          ! and inhomogeneous pure Dirichlet boundary conditions for the velocity set
          ! according to the analytic solution (u1, u2, p) = (t^2, t^2, -2t(x+y-1)):
          ! Without loss of generality let T=[0,2] and time step size t_m = 2. Then, one
          ! sweep of SDIRK3PR suffices to reach T=2, starting from T=0. The velocity
          ! solution is exactly approximated in every intermediate stage, but the pressure
          ! solution is only correct when stage 5 is completed.
          !
          ! Note: Divisor 5 = nsubsteps, given that this time stepping scheme has 5 stages
          !       with an implicit first stage such that really 5 time steps are needed
          !       for a complete sweep.
          iglobaltimestep = int(rproblem%rtimedependence%itimestep / 5)
          select case (mod(rproblem%rtimedependence%itimestep,5))
          case (0)
            cignoreTimeStep = 0
          case (1,2,3,4)
            cignoreTimeStep = 1
          end select
        end if
      end if
    end if


    ! The error analysis might take place at an arbitrary time.
    ! Therefore, modify the "current" time to the time where we want to
    ! do the analysis. Save the "current" time and restore afterwards.
    ! This is necessary, as the user callback routines take the time
    ! for evaluation from here.
    dtimebackup = rproblem%rtimedependence%dtime
    rproblem%rtimedependence%dtime = dtime

    if ((rpostprocessing%icalcL2 .ne. 0) .or. &
        (rpostprocessing%icalcH1 .ne. 0) .or. &
        (icalcEnergy .ne. 0)) then
      call output_lbrk()
      call output_line ("Error Analysis")
      call output_line ("--------------")
    end if

    ! When writing to a file is enabled, delete the file in the first timestep.
    cflag = SYS_APPEND
    if (rproblem%rtimedependence%itimeStep .eq. 0) cflag = SYS_REPLACE

    ! ===============================================================
    ! Error computation. L2 error
    ! ===============================================================

    if (rpostprocessing%icalcL2 .ne. 0) then

      select case (rpostprocessing%icalcL2)

      ! -------------------------------
      ! Compute using callback functions
      ! -------------------------------
      case (1)
        call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

        ! Calculate space discretisation error

        ! Perform error analysis to calculate and add 1/2||u-z||_{L^2}.
        call pperr_scalar (PPERR_L2ERROR,Derr(1),rsolution%RvectorBlock(1),&
                          ffunction_TargetX,rproblem%rcollection,&
                          rcubatureInfo=rcubatureInfoUV)

        call pperr_scalar (PPERR_L2ERROR,Derr(2),rsolution%RvectorBlock(2),&
                          ffunction_TargetY,rproblem%rcollection,&
                          rcubatureInfo=rcubatureInfoUV)

        derrorVel = sqrt(Derr(1)**2+Derr(2)**2)

        call pperr_scalar (PPERR_L2ERROR,Derr(3),rsolution%RvectorBlock(3),&
                          ffunction_TargetP,rproblem%rcollection,&
                          rcubatureInfo=rcubatureInfoP)

        derrorP = Derr(3)


        ! Calculate accumulated time plus space discretisation error
        if (rproblem%itimedependence .ne. 0 .and. &
            rpostprocessing%icalcTimeSpaceDiscrErrors .ne. 0 .and. &
            iglobaltimestep .ge. 0 .and. cignoreTimeStep .eq. 0) then
          call evalTimeSpaceError(&
               ! Evaluate analytic reference solution via hardcoded callback functions
               1, &
               ! Evaluate L2 error
               PPERR_L2ERROR, &
               ! previous L2 errors (used by trapezoidal rule only)
               derrorL2Vel_old, derrorL2P_old, &
               rproblem%rcollection, derrorL2VelTimeSpace, derrorL2PTimeSpace)
        end if


        call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

      ! -------------------------------
      ! Compute via expressions
      ! -------------------------------
      case (2)

        ! Prepare a collection
        rlocalCollection%IquickAccess(1) = 0       ! L2-error
        rlocalCollection%DquickAccess(:) = 0.0_DP
        rlocalCollection%DquickAccess(1) = dtime   ! Time
        rlocalCollection%p_rfparserQuickAccess1 => rpostprocessing%rrefFunctionL2

        ! Calculate space discretisation error

        ! Perform error analysis to calculate and add 1/2||u-z||_{L^2}.
        ! Use the callback function "fcalc_error".
        rlocalCollection%IquickAccess(2) = 1  ! component
        call pperr_scalar (PPERR_L2ERROR,Derr(1),rsolution%RvectorBlock(1),&
                          fcalc_error,rlocalCollection,&
                          rcubatureInfo=rcubatureInfoUV)

        rlocalCollection%IquickAccess(2) = 2  ! component
        call pperr_scalar (PPERR_L2ERROR,Derr(2),rsolution%RvectorBlock(2),&
                          fcalc_error,rlocalCollection,&
                          rcubatureInfo=rcubatureInfoUV)

        derrorVel = sqrt(Derr(1)**2+Derr(2)**2)

        rlocalCollection%IquickAccess(2) = 3  ! component
        call pperr_scalar (PPERR_L2ERROR,Derr(3),rsolution%RvectorBlock(3),&
                          fcalc_error,rlocalCollection,&
                          rcubatureInfo=rcubatureInfoP)

        derrorP = Derr(3)


        ! Calculate accumulated time plus space discretisation error
        if (rproblem%itimedependence .ne. 0 .and. &
            rpostprocessing%icalcTimeSpaceDiscrErrors .ne. 0 .and. &
            iglobaltimestep .ge. 0 .and. cignoreTimeStep .eq. 0) then
          call evalTimeSpaceError(&
               ! Evaluate analytic reference solution via expressions
               2, &
               ! Evaluate L2 error
               PPERR_L2ERROR, &
               ! previous L2 errors (used by trapezoidal rule only)
               derrorL2Vel_old, derrorL2P_old, &
               rlocalCollection, derrorL2VelTimeSpace, derrorL2PTimeSpace)
        end if

      end select

      call output_line ("||u-reference||_L2 = "//trim(sys_sdEP(derrorVel,15,6)) )
      if (cignoreTimeStep .eq. 0) then
        call output_line ("||p-reference||_L2 = "//trim(sys_sdEP(derrorP,15,6)) )
      else
        call output_line ("||p-reference||_L2 = infeasible in intermediate time steps")
      end if

      if (rproblem%itimedependence .ne. 0 .and. &
          rpostprocessing%icalcTimeSpaceDiscrErrors .ne. 0) then
        if (cignoreTimeStep .ne. 0) then
          call output_line ("||u-reference||_{0,T;L2(Omega)} = " // &
               "not evaluated in intermediate time steps")
          call output_line ("||p-reference||_{0,T;L2(Omega)} = " // &
               "infeasible in intermediate time steps")
        else
          ! At T=T_start there is no time error yet. So skip output in that case.
          select case (rpostprocessing%itimeSpaceDiscrErrorMethod)
          case (0) ! l2 error
            call output_line ("||u-reference||_l2(0,T;L2(Omega)) = " // &
                              trim(sys_sdEP(sqrt(derrorL2VelTimeSpace / &
                                                 (iglobaltimestep + 1.0_DP)), 15, 6)) )
            call output_line ("||p-reference||_l2(0,T;L2(Omega)) = " // &
                              trim(sys_sdEP(sqrt(derrorL2PTimeSpace / &
                                                 (iglobaltimestep + 1.0_DP)), 15, 6)) )

          case (1) ! Trapezoidal rule for numeric integration means that only 2 data
                   ! points are needed. Which are available in time step 1: start solution
                   ! and the solution in time step 1
            call output_line ("||u-reference||_{0,T;L2(Omega)} = " // &
                              trim(sys_sdEP(sqrt(derrorL2VelTimeSpace), 15, 6)) )
            call output_line ("||p-reference||_{0,T;L2(Omega)} = " // &
                              trim(sys_sdEP(sqrt(derrorL2PTimeSpace), 15, 6)) )

            ! Save this time step`s space discretisation error for re-use in the next
            ! time step:
            derrorL2Vel_old = derrorVel
            derrorL2P_old = derrorP

          case (2,3) ! Quadratic Lagrange time-interpolation of the solution requires 3
                     ! data points. Before completion of time step 2 is hence no error
                     ! calculation available.
            if (iglobaltimestep .eq. 1) then
              call output_line ("||u-reference||_{0,T;L2(Omega)} = (postponed)  " // &
                                "(time plus space discretisation")
              call output_line ("||p-reference||_{0,T;L2(Omega)} = (postponed)  " // &
                                " errors are only available from")
              call output_line ("                                               " // &
                                " time step 2 onwards.)")
            else if (iglobaltimestep .ge. 2) then
              ! iglobaltimestep >= 2, cignoreTimeStep = 0
              call output_line ("||u-reference||_{0,T;L2(Omega)} = " // &
                                trim(sys_sdEP(sqrt(derrorL2VelTimeSpace), 15, 6)) )
              call output_line ("||p-reference||_{0,T;L2(Omega)} = " // &
                                trim(sys_sdEP(sqrt(derrorL2PTimeSpace), 15, 6)) )
            end if

          case (4,5) ! Cubic Lagrange time-interpolation of the solution requires 4
                     ! data points. Before completion of time step 3 is hence no error
                     ! calculation available.
            if (iglobaltimestep .lt. 3) then
              call output_line ("||u-reference||_{0,T;L2(Omega)} = (postponed)  " // &
                                "(time plus space discretisation")
              call output_line ("||p-reference||_{0,T;L2(Omega)} = (postponed)  " // &
                                " errors are only available from")
              call output_line ("                                               " // &
                                " time step 3 onwards.)")
            else if (iglobaltimestep .ge. 3) then
              ! iglobaltimestep >= 3, cignoreTimeStep = 0
              call output_line ("||u-reference||_{0,T;L2(Omega)} = " // &
                                trim(sys_sdEP(sqrt(derrorL2VelTimeSpace), 15, 6)) )
              call output_line ("||p-reference||_{0,T;L2(Omega)} = " // &
                                trim(sys_sdEP(sqrt(derrorL2PTimeSpace), 15, 6)) )
            end if
          end select
        end if
      end if

      ! -------------------------------
      ! Write to file
      ! -------------------------------
      if (iwriteErrorAnalysisL2 .ne. 0) then
        ! Write the result to a text file.
        ! Format: timestep current-time value
        call io_openFileForWriting(sfilenameErrorAnalysisL2, iunit, &
            cflag, bfileExists,.true.)
        if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
          ! Write a headline
! SB: Not sure whether we would want to extend this, too
          write (iunit,"(A)") "# timestep time ||u-reference||_L2 ||p-reference||_L2"
        end if
        stemp = trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // " " &
            // trim(sys_sdEL(dtime,10)) // " " &
            // trim(sys_sdEL(derrorVel,10)) // " " &
            // trim(sys_sdEL(derrorP,10))
        write (iunit,"(A)") trim (stemp)
        close (iunit)
      end if

    end if

    ! ===============================================================
    ! Error computation. H1 error
    ! ===============================================================

    if (rpostprocessing%icalcH1 .ne. 0) then

      select case (rpostprocessing%icalcH1)

      ! -------------------------------
      ! Compute using callback functions
      ! -------------------------------
      case (1)
        call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

        ! Calculate space discretisation error

        ! Perform error analysis to calculate and add ||u-z||_{H^1}.
        call pperr_scalar (PPERR_H1ERROR,Derr(1),rsolution%RvectorBlock(1),&
                          ffunction_TargetX,rproblem%rcollection,&
                          rcubatureInfo=rcubatureInfoUV)

        call pperr_scalar (PPERR_H1ERROR,Derr(2),rsolution%RvectorBlock(2),&
                          ffunction_TargetY,rproblem%rcollection,&
                          rcubatureInfo=rcubatureInfoUV)

        derrorVel = -1.0_DP
!SB: Note: This if condition is never true as pperr_scalar_conf() runs sqrt() on return value
!          for PPERR_H1ERROR!
        if ((Derr(1) .ne. -1.0_DP) .and. (Derr(2) .ne. -1.0_DP)) then
          derrorVel = sqrt(Derr(1)**2+Derr(2)**2)
        end if

        call pperr_scalar (PPERR_H1ERROR,Derr(3),rsolution%RvectorBlock(3),&
                          ffunction_TargetP,rproblem%rcollection,&
                          rcubatureInfo=rcubatureInfoP)

        derrorP = -1.0_DP
!SB: Note: This if condition is never true as pperr_scalar_conf() runs sqrt() on return value
!          for PPERR_H1ERROR!
        if (Derr(3) .ne. -1.0_DP) then
          derrorP = Derr(3)
        end if


        ! Calculate accumulated time plus space discretisation error
        if (rproblem%itimedependence .ne. 0 .and. &
            rpostprocessing%icalcTimeSpaceDiscrErrors .ne. 0 .and. &
            iglobaltimestep .ge. 0 .and. cignoreTimeStep .eq. 0) then
          call evalTimeSpaceError(&
               ! Evaluate analytic reference solution via hardcoded callback functions
               1, &
               ! Evaluate H1 error
               PPERR_H1ERROR, &
               ! previous H1 errors (used by trapezoidal rule only)
               derrorH1Vel_old, derrorH1P_old, &
               rproblem%rcollection, derrorH1VelTimeSpace, derrorH1PTimeSpace)
        end if

        call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

      ! -------------------------------
      ! Compute via expressions
      ! -------------------------------
      case (2)

        ! Prepare a collection
        rlocalCollection%IquickAccess(1) = 1       ! H1-error
        rlocalCollection%DquickAccess(:) = 0.0_DP
        rlocalCollection%DquickAccess(1) = dtime   ! Time
        rlocalCollection%p_rfparserQuickAccess1 => rpostprocessing%rrefFunctionH1

        ! Calculate space discretisation error

        ! Perform error analysis to calculate and add ||u-z||_{H^1}.
        ! Use the callback function "fcalc_error".
        rlocalCollection%IquickAccess(2) = 1  ! component
        call pperr_scalar (PPERR_H1ERROR,Derr(1),rsolution%RvectorBlock(1),&
                          fcalc_error,rlocalCollection,&
                          rcubatureInfo=rcubatureInfoUV)

        rlocalCollection%IquickAccess(2) = 2  ! component
        call pperr_scalar (PPERR_H1ERROR,Derr(2),rsolution%RvectorBlock(2),&
                          fcalc_error,rlocalCollection,&
                          rcubatureInfo=rcubatureInfoUV)

        derrorVel = -1.0_DP
!SB: Note: This if condition is never true as pperr_scalar_conf() runs sqrt() on return value
!          for PPERR_H1ERROR!
        if ((Derr(1) .ne. -1.0_DP) .and. (Derr(2) .ne. -1.0_DP)) then
          derrorVel = sqrt(Derr(1)**2+Derr(2)**2)
        end if

        rlocalCollection%IquickAccess(2) = 3  ! component
        call pperr_scalar (PPERR_H1ERROR,Derr(3),rsolution%RvectorBlock(3),&
                          fcalc_error,rlocalCollection,&
                          rcubatureInfo=rcubatureInfoP)

        derrorP = -1.0_DP
!SB: Note: This if condition is never true as pperr_scalar_conf() runs sqrt() on return value
!          for PPERR_H1ERROR!
        if (Derr(3) .ne. -1.0_DP) then
          derrorP = Derr(3)
        end if


        ! Calculate accumulated time plus space discretisation error
        if (rproblem%itimedependence .ne. 0 .and. &
            rpostprocessing%icalcTimeSpaceDiscrErrors .ne. 0 .and. &
            iglobaltimestep .ge. 0 .and. cignoreTimeStep .eq. 0) then
          call evalTimeSpaceError(&
               ! Evaluate analytic reference solution via expressions
               2, &
               ! Evaluate H1 error
               PPERR_H1ERROR, &
               ! previous H1 errors (used by trapezoidal rule only)
               derrorH1Vel_old, derrorH1P_old, &
               rlocalCollection, derrorH1VelTimeSpace, derrorH1PTimeSpace)
        end if

      end select

      ! derrorVel/derrorP=-1 indicates that the error is not available.
      if (derrorVel .ne. -1.0_DP) then
        call output_line ("||u-reference||_H1 = "//trim(sys_sdEP(derrorVel,15,6)),&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      else
!SB: Note: This if condition is never true as pperr_scalar_conf() runs sqrt() on return value
!          for PPERR_H1ERROR!
        call output_line ("||u-reference||_H1 = not available",&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      end if

      if (derrorP .ne. -1.0_DP) then
        call output_line ("||p-reference||_H1 = "//trim(sys_sdEP(derrorP,15,6)),&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      else
!SB: Note: This if condition is never true as pperr_scalar_conf() runs sqrt() on return value
!          for PPERR_H1ERROR!
        call output_line ("||p-reference||_H1 = not available",&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      end if


      if (rproblem%itimedependence .ne. 0 .and. &
          rpostprocessing%icalcTimeSpaceDiscrErrors .ne. 0) then
        if (cignoreTimeStep .ne. 0) then
          call output_line ("||u-reference||_{0,T;H1(Omega)} = " // &
               "not evaluated in intermediate time steps")
          call output_line ("||p-reference||_{0,T;H1(Omega)} = " // &
               "infeasible in intermediate time steps")
        else
          ! At T=T_start there is no time error yet. So skip output in that case.
          select case (rpostprocessing%itimeSpaceDiscrErrorMethod)
          case (0) ! l2 error
            call output_line ("||u-reference||_l2(0,T;H1(Omega)) = " // &
                              trim(sys_sdEP(sqrt(derrorH1VelTimeSpace / &
                                                 (iglobaltimestep + 1.0_DP)), 15, 6)) )
            call output_line ("||p-reference||_l2(0,T;H1(Omega)) = " // &
                              trim(sys_sdEP(sqrt(derrorH1PTimeSpace / &
                                                 (iglobaltimestep + 1.0_DP)), 15, 6)) )

          case (1) ! Trapezoidal rule for numeric integration means that only 2 data
                   ! points are needed. Which are available in time step 1: start solution
                   ! and the solution in time step 1
            call output_line ("||u-reference||_{0,T;H1(Omega)} = " // &
                              trim(sys_sdEP(sqrt(derrorH1VelTimeSpace), 15, 6)) )
            call output_line ("||p-reference||_{0,T;H1(Omega)} = " // &
                              trim(sys_sdEP(sqrt(derrorH1PTimeSpace), 15, 6)) )

            ! Save this time step`s space discretisation error for re-use in the next time
            ! step:
            derrorH1Vel_old = derrorVel
            derrorH1P_old = derrorP

          case (2,3) ! Quadratic Lagrange time-interpolation of the solution requires 3
                     ! data points. Before completion of time step 2 is hence no error
                     ! calculation available.
            if (iglobaltimestep .eq. 1) then
              call output_line ("||u-reference||_{0,T;H1(Omega)} = (postponed)  " // &
                                "(time plus space discretisation")
              call output_line ("||p-reference||_{0,T;H1(Omega)} = (postponed)  " // &
                                " errors are only available from")
              call output_line ("                                               " // &
                                " time step 2 onwards.)")
            else if (iglobaltimestep .ge. 2) then
              ! iglobaltimestep >= 2, cignoreTimeStep = 0
              call output_line ("||u-reference||_{0,T;H1(Omega)} = " // &
                                trim(sys_sdEP(sqrt(derrorH1VelTimeSpace), 15, 6)) )
              call output_line ("||p-reference||_{0,T;H1(Omega)} = " // &
                                trim(sys_sdEP(sqrt(derrorH1PTimeSpace), 15, 6)) )
            end if

          case (4,5) ! Cubic Lagrange time-interpolation of the solution requires 4
                     ! data points. Before completion of time step 3 is hence no error
                     ! calculation available.
            if (iglobaltimestep .lt. 3) then
              call output_line ("||u-reference||_{0,T;H1(Omega)} = (postponed)  " // &
                                "(time plus space discretisation")
              call output_line ("||p-reference||_{0,T;H1(Omega)} = (postponed)  " // &
                                " errors are only available from")
              call output_line ("                                               " // &
                                " time step 3 onwards.)")
            else if (iglobaltimestep .ge. 3) then
              ! iglobaltimestep >= 3, cignoreTimeStep = 0
              call output_line ("||u-reference||_{0,T;H1(Omega)} = " // &
                                trim(sys_sdEP(sqrt(derrorH1VelTimeSpace), 15, 6)) )
              call output_line ("||p-reference||_{0,T;H1(Omega)} = " // &
                                trim(sys_sdEP(sqrt(derrorH1PTimeSpace), 15, 6)) )
            end if
          end select
        end if
      end if

      ! -------------------------------
      ! Write to file
      ! -------------------------------
      if (iwriteErrorAnalysisH1 .ne. 0) then
        ! Write the result to a text file.
        ! Format: timestep current-time value
        call io_openFileForWriting(sfilenameErrorAnalysisH1, iunit, &
            cflag, bfileExists,.true.)
        if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
          ! Write a headline
! SB: Not sure whether we would want to extend this, too
          write (iunit,"(A)") "# timestep time ||u-reference||_H1"
        end if
! SB: Not sure whether we would want to extend this, too
        stemp = trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // " " &
            // trim(sys_sdEL(dtime,10)) // " " &
            // trim(sys_sdEL(derrorVel,10))
        write (iunit,"(A)") trim(stemp)
        close (iunit)
      end if

    end if


    if (rproblem%itimedependence .ne. 0 .and. &
        rpostprocessing%icalcTimeSpaceDiscrErrors .ne. 0 .and. &
        (rpostprocessing%icalcL2 .ne. 0 .or. rpostprocessing%icalcH1 .ne. 0) .and. &
        ! Depending on the particular time stepping scheme, it may be necessary to skip
        ! some intermediate steps. Not the case, though, for the Fractional Step Theta
        ! scheme. Omitting the two intermediate steps leads to overestimation of error
        ! reduction rates!
        cignoreTimeStep .eq. 0) then

      call parlst_getvalue_int (rproblem%rparamList,"CC-DISCRETISATION",&
           "ctypeInitialSolution",ctypeInitialSolution,0)
      if (rpostprocessing%itimeSpaceDiscrErrorMethod .ge. 2 .and. &
          ctypeInitialSolution .lt. 3) then
        call output_line ("Code is configured to assume that the analytic", &
             OU_CLASS_WARNING, OU_MODE_STD, "cc_errorAnalysis")
        call output_line (" solution at T=" // &
             trim(sys_sdEL(rproblem%rtimedependence%dtimeInit,10)) // &
             " equals the start", &
             OU_CLASS_WARNING, OU_MODE_STD, "cc_errorAnalysis")
        call output_line ("solution. A safer choice is", &
             OU_CLASS_WARNING, OU_MODE_STD, "cc_errorAnalysis")
        call output_line ("   [CC-DISCRETISATION]", &
             OU_CLASS_WARNING, OU_MODE_STD, "cc_errorAnalysis")
        call output_line ("   ctypeInitialSolution = 3", &
             OU_CLASS_WARNING, OU_MODE_STD, "cc_errorAnalysis")
      end if

      select case (rpostprocessing%itimeSpaceDiscrErrorMethod)
      case (0) ! nothing needs to be stored to calculate l2 error
      case (1) ! Trapezoidal rule for numeric integration means that only
               !    derrorVel, derrorL2Vel_old, derrorH1Vel_old
               !    derrorP, derrorL2P_old, derrorH1P_old
               ! are needed; no explicit storing of older solutions required
        if (iglobaltimestep .eq. 0) then
          ! When applying the trapezoidal rule the first time in time step 1, we will need
          ! this value.
          doldTime = rproblem%rtimedependence%dtimeInit

        else if (iglobaltimestep .gt. 0) then
          doldTime = dtime
        end if

      case (2,3) ! Quadratic Lagrange time-interpolation of the solution requires 3 data
                 ! points: the solution from the current time step one and the two
                 ! preceding time steps.
        if (iglobaltimestep .eq. 0) then

          call lsysbl_createVector (rsolution, rpostprocessing%rolderSolution, .false.)
          call lsysbl_createVector (rsolution, rpostprocessing%roldSolution, .false.)
          ! Store current solution which will be used in time step 2 as penultimate time
          ! step
          call lsysbl_copyVector (rsolution, rpostprocessing%rolderSolution)
          dolderTime = dtime

        else if (iglobaltimestep .eq. 1) then
          ! Store current solution which will be used in time step 2 as the one from the
          ! previous time step
          call lsysbl_copyVector (rsolution, rpostprocessing%roldSolution)
          doldTime = dtime

        ! For (non-intermediate) time steps 1 and 2 we already saved the solution
        ! snapshots in the right order. For later time steps shift them by one such that
        ! the oldest solution snapshot, i.e. from t^{n-1}, is evicted.
        else if (iglobaltimestep .ge. 2) then
          call lsysbl_copyVector(rpostprocessing%roldSolution, &
                                 rpostprocessing%rolderSolution)
          call lsysbl_copyVector(rsolution, rpostprocessing%roldSolution)
          dolderTime = doldTime
          doldTime = dtime
        end if

      case (4,5) ! Cubic Lagrange time-interpolation of the solution requires 4 data
                 ! points: the solution from the current time step one and the three
                 ! preceding time steps.
        if (iglobaltimestep .eq. 0) then
          call lsysbl_createVector (rsolution, rpostprocessing%roldestSolution, .false.)
          call lsysbl_createVector (rsolution, rpostprocessing%rolderSolution, .false.)
          call lsysbl_createVector (rsolution, rpostprocessing%roldSolution, .false.)
          ! Store current solution which will be used in time step 3 as antepenultimate
          ! time step
          call lsysbl_copyVector (rsolution, rpostprocessing%roldestSolution)
          doldestTime = dtime

        else if (iglobaltimestep .eq. 1) then
          ! Store current solution which will be used in time step 3 as the one from the
          ! previous time step
          call lsysbl_copyVector (rsolution, rpostprocessing%rolderSolution)
          dolderTime = dtime

        else if (iglobaltimestep .eq. 2) then
          ! Store current solution which will be used in time step 3 as the one from the
          ! previous time step
          call lsysbl_copyVector (rsolution, rpostprocessing%roldSolution)
          doldTime = dtime

        ! For time steps 0, 1 and 2 we already saved the solution snapshots in the right
        ! order. For later time steps shift them by one such that the oldest solution
        ! snapshot, i.e. from t^{n-2}, is evicted.
        else if (iglobaltimestep .ge. 3) then
          call lsysbl_copyVector(rpostprocessing%rolderSolution, &
                                 rpostprocessing%roldestSolution)
          call lsysbl_copyVector(rpostprocessing%roldSolution, &
                                 rpostprocessing%rolderSolution)
          call lsysbl_copyVector(rsolution, rpostprocessing%roldSolution)
          doldestTime = dolderTime
          dolderTime = doldTime
          doldTime = dtime
        end if

      end select

    end if

    ! ===============================================================
    ! Kinetic energy
    ! ===============================================================

    if (icalcEnergy .ne. 0) then

      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

      ! Perform error analysis to calculate and add 1/2||u||^2_{L^2}.
      call pperr_scalar (PPERR_L2ERROR,Derr(1),rsolution%RvectorBlock(1),&
          rcubatureInfo=rcubatureInfoUV)
      call pperr_scalar (PPERR_L2ERROR,Derr(2),rsolution%RvectorBlock(2),&
          rcubatureInfo=rcubatureInfoUV)

      denergy = 0.5_DP*(Derr(1)**2+Derr(2)**2)

      call output_line ("||u_1||_L2         = "//trim(sys_sdEP(Derr(1),15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      call output_line ("||u_2||_L2         = "//trim(sys_sdEP(Derr(2),15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      call output_line ("||u||_L2           = "//&
          trim(sys_sdEP(sqrt(Derr(1)**2+Derr(2)**2),15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      call output_line ("1/2||u||^2_L2      = "//trim(sys_sdEP(denergy,15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )

      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

      if (iwriteKineticEnergy .ne. 0) then
        ! Write the result to a text file.
        ! Format: timestep current-time value
        call io_openFileForWriting(sfilenameKineticEnergy, iunit, &
            cflag, bfileExists,.true.)
        if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
          ! Write a headline
          write (iunit,"(A)") "# timestep time 1/2||u||^2_L2 ||u||_L2 ||u_1||_L2 ||u_2||_L2"
        end if
        stemp = trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // " " &
            // trim(sys_sdEL(dtime,10)) // " " &
            // trim(sys_sdEL(denergy,10)) // " " &
            // trim(sys_sdEL(sqrt(Derr(1)**2+Derr(2)**2),10)) // " " &
            // trim(sys_sdEL(Derr(1),10)) // " " &
            // trim(sys_sdEL(Derr(2),10))
        write (iunit,"(A)") trim(stemp)
        close (iunit)
      end if

    end if

    ! Release the cubature information structures
    call spdiscr_releaseCubStructure(rcubatureInfoUV)
    call spdiscr_releaseCubStructure(rcubatureInfoP)

    ! Restore the "current" time.
    rproblem%rtimedependence%dtime = dtimebackup


  contains

!<subroutine>

    subroutine evalTimeSpaceError(ccallbackMethod, cerrorType, &
         derrorOldVel, derrorOldP, &
         rcollection, dtimeSpaceErrorVel, dtimeSpaceErrorP)

!<description>
    ! Calculate accumulated time plus space discretisation error for velocity and pressure
    ! solution.
    ! The square of this error accumulated in dtimeSpaceErrorVel and dtimeSpaceErrorP.
!</description>

!<input>
      ! Type of callback to evaluate analytic reference solution
      ! =1: hardcoded in ffunction_TargetX, ffunction_TargetY and ffunction_TargetP
      ! =2: use function parser to evaluate reference function given in configuration file
      integer, intent(in) :: ccallbackMethod

      ! Type of error to compute. Bitfield. This is a combination of the
      ! PPERR_xxxx-constants, which specifies what to compute.
      ! Example: PPERR_L2ERROR computes the $L_2$-error.
      ! (used by all numerical integration schemes except trapezoidal rule)
      integer, intent(in) :: cerrorType

      ! Space discretisation error for velocity, measured in error norm cerrorType, from
      ! previous time step (used by trapezoidal rule only)
      real(DP), intent(in) :: derrorOldVel

      ! Space discretisation error for pressure, measured in error norm cerrorType, from
      ! previous time step (used by trapezoidal rule only)
      real(DP), intent(in) :: derrorOldP
!</input>

!<inputoutput>
      ! A pointer to a collection structure to provide additional
      ! information to the coefficient routine.
      type(t_collection), intent(inout) :: rcollection

      ! Time plus space discretisation error for velocity, measured in error norm
      ! cerrorType
      real(DP), intent(inout) :: dtimeSpaceErrorVel

      ! Time plus space discretisation error for pressure, measured in error norm
      ! cerrorType
      real(DP), intent(inout) :: dtimeSpaceErrorP
!</inputoutput>

!</subroutine>

      real(DP) :: derrorVel_GaussPt1, derrorVel_GaussPt2, derrorVel_GaussPt3
      real(DP) :: derrorP_GaussPt1, derrorP_GaussPt2, derrorP_GaussPt3
      integer :: igaussQuadDegree

      igaussQuadDegree = 0
      select case (rpostprocessing%itimeSpaceDiscrErrorMethod)
      case (0) ! build square of l2 error: sum up square of (L2 or H1, per time step)
               ! space discretisation error, ! divide this sum later by the number of time
               ! steps performed
        if (iglobaltimestep .ge. 0 .and. cignoreTimeStep .eq. 0) then
          dtimeSpaceErrorVel = dtimeSpaceErrorVel + derrorVel**2
          dtimeSpaceErrorP   = dtimeSpaceErrorP   + derrorP**2
        end if


      case (1)
        ! Approximate
        !    || sol(t) - sol_ref(t) ||^2_L2(0,T;\Omega)
        !    = sqrt( \int_{t^{start}}^{t^{end}} || sol(t) - sol_ref(t) ||^2_L2 dt )
        ! by trapezoidal rule
        !    sum_{all time steps} time step * 1/2 * (sol^{n+1}(x) + sol^{n}(x))
        ! taking as sol^{n+1}(x) in every time step t^{n+1} (including t^{1} and t^{end})
        ! the values of
        !       || u(t^{n+1}) - u_ref(t^{n+1}) ||^2_L2
        ! and
        !       || p(t^{n+1}) - p_ref(t^{n+1}) ||^2_L2
        ! that got just calculated in derrorVel and derrorP and as sol^{n}(x) the
        ! corresponding values from the previous time step (for the start solution holds
        ! n=0 and both derrorOldVel and derrorOldP are assumed to be initialised to zero).
        if (iglobaltimestep .gt. 0 .and. cignoreTimeStep .eq. 0) then
          dtimeSpaceErrorVel = dtimeSpaceErrorVel + &
               (dtime - doldTime) * 0.5_DP * (derrorVel**2 + derrorOldVel**2)
          dtimeSpaceErrorP = dtimeSpaceErrorP + &
               (dtime - doldTime) * 0.5_DP * (derrorP**2 + derrorOldP**2)
        end if


      case (2,3)
        ! Time-interpolate the solution with a quadratic Lagrange polynomial, then
        ! evaluate
        !    || sol(t) - sol_ref(t) ||^2_L2(0,T;\Omega)
        !    = sqrt( \int_{t^{start}}^{t^{end}} || sol(t) - sol_ref(t) ||^2_L2 dt )
        ! with 2-point or 3-point Gaussian quadrature rule.
        if (iglobaltimestep .ge. 2 .and. cignoreTimeStep .eq. 0) then

          ! Create temporary vectors
          call lsysbl_createVector (rsolution, rvector_auxSum, .false.)
          call lsysbl_createVector (rsolution, rvector_sol_GaussPt1, .false.)
          call lsysbl_createVector (rsolution, rvector_sol_GaussPt2, .false.)
          if (rpostprocessing%itimeSpaceDiscrErrorMethod .eq. 2) then
            igaussQuadDegree = 2
          else
            igaussQuadDegree = 3
            call lsysbl_createVector (rsolution, rvector_sol_GaussPt3, .false.)
          end if

          if (iglobaltimestep .eq. 2) then
            ! Special case treatment

            ! Calculate
            !   \int_t^{0}^t^{1} || u(t) - u_ref(t) ||^2_L2 dt
            ! and
            !   \int_t^{0}^t^{1} || p(t) - p_ref(t) ||^2_L2 dt
            !
            ! Why is the time interval [t^{0}, t^{1}] a special case?
            ! We need to make up for the fact of not having been able to calculate time
            ! discretisation error in time step 1 due to missing data points for quadratic
            ! Lagrange time-interpolation of the solution.
            !
            ! Time plus space discretisation error for velocity in interval [t^{0}, t^{1}]
            select case (ccallbackMethod)
            case (1)
              call evalTimeSpaceErrorQuadLagrange(cerrorType, 1, &
                   igaussQuadDegree, ffunction_TargetX, rcollection, dolderTime, doldTime)
              call evalTimeSpaceErrorQuadLagrange(cerrorType, 2, &
                   igaussQuadDegree, ffunction_TargetY, rcollection, dolderTime, doldTime)

            case (2)
              rcollection%IquickAccess(2) = 1  ! component
              call evalTimeSpaceErrorQuadLagrange(cerrorType, 1, &
                   igaussQuadDegree, fcalc_error, rcollection, dolderTime, doldTime)
              rcollection%IquickAccess(2) = 2  ! component
              call evalTimeSpaceErrorQuadLagrange(cerrorType, 2, &
                   igaussQuadDegree, fcalc_error, rcollection, dolderTime, doldTime)
            end select

            if (igaussQuadDegree .eq. 2) then
              derrorVel_GaussPt1 = sqrt(Derr(1)**2 + Derr(2)**2)
              derrorVel_GaussPt2 = sqrt(Derr2(1)**2 + Derr2(2)**2)

              ! Apply Gauss quadrature rule of degree n=2
              ! which is defined on the interval [-1,1]. Transforming it to [t^{0}, t^{1}]
              ! leads to a factor of
              !               1/2 (t^{1} - t^{0})
              ! due to the rules for integration by substitution. Weighing factors of
              ! summands for Gauss quadrature rule of degree n=2 are both 1.0.
              dtimeSpaceErrorVel = dtimeSpaceErrorVel + &
                   0.5_DP*(doldTime - dolderTime) * &
                   (1.0_DP * derrorVel_GaussPt1**2 + &
                    1.0_DP * derrorVel_GaussPt2**2)

            else ! 3-Point Gauss quadrature
              derrorVel_GaussPt1 = sqrt(Derr(1)**2 + Derr(2)**2)
              derrorVel_GaussPt2 = sqrt(Derr2(1)**2 + Derr2(2)**2)
              derrorVel_GaussPt3 = sqrt(Derr3(1)**2 + Derr3(2)**2)

              ! Apply Gauss quadrature rule of degree n=3
              ! which is defined on the interval [-1,1]. Transforming it to [t^{0}, t^{1}]
              ! leads to a factor of
              !               1/2 (t^{1} - t^{0})
              ! due to the rules for integration by substitution. Weighing factors of
              ! summands for Gauss quadrature rule of degree n=3 are 5/9, 8/9, 5/9.
              dtimeSpaceErrorVel = dtimeSpaceErrorVel + &
                   0.5_DP*(doldTime - dolderTime) * &
                   (5.0_DP/9.0_DP * derrorVel_GaussPt1**2 + &
                    8.0_DP/9.0_DP * derrorVel_GaussPt2**2 + &
                    5.0_DP/9.0_DP * derrorVel_GaussPt3**2)
            end if ! 2-Point or 3-Point Gauss quadrature
            !   \int_t^{0}^t^{1} || u(t) - u_ref(t) ||^2_L2 dt
            ! calculated.


            ! Time plus space discretisation error for pressure in interval [t^{0}, t^{1}]
            select case (ccallbackMethod)
            case (1)
              call evalTimeSpaceErrorQuadLagrange(cerrorType, 3, &
                   igaussQuadDegree, ffunction_TargetP, rcollection, dolderTime, doldTime)

            case (2)
              rcollection%IquickAccess(2) = 3  ! component
              call evalTimeSpaceErrorQuadLagrange(cerrorType, 3, &
                   igaussQuadDegree, fcalc_error, rcollection, dolderTime, doldTime)
            end select

            if (igaussQuadDegree .eq. 2) then
              derrorP_GaussPt1 = Derr(3)
              derrorP_GaussPt2 = Derr2(3)

              ! Apply Gauss quadrature rule of degree n=2
              ! which is defined on the interval [-1,1]. Transforming it to [t^{0}, t^{1}]
              ! leads to a factor of
              !               1/2 (t^{1} - t^{0})
              ! due to the rules for integration by substitution. Weighing factors of
              ! summands for Gauss quadrature rule of degree n=2 are both 1.0.
              dtimeSpaceErrorP = dtimeSpaceErrorP + &
                   0.5_DP*(doldTime - dolderTime) * &
                   (1.0_DP * derrorP_GaussPt1**2 + &
                    1.0_DP * derrorP_GaussPt2**2)

            else ! 3-Point Gauss quadrature
              derrorP_GaussPt1 = Derr(3)
              derrorP_GaussPt2 = Derr2(3)
              derrorP_GaussPt3 = Derr3(3)

              ! Apply Gauss quadrature rule of degree n=3
              ! which is defined on the interval [-1,1]. Transforming it to [t^{0}, t^{1}]
              ! leads to a factor of
              !               1/2 (t^{1} - t^{0})
              ! due to the rules for integration by substitution. Weighing factors of
              ! summands for Gauss quadrature rule of degree n=3 are 5/9, 8/9, 5/9.
              dtimeSpaceErrorP = dtimeSpaceErrorP + &
                   0.5_DP*(doldTime - dolderTime) * &
                   (5.0_DP/9.0_DP * derrorP_GaussPt1**2 + &
                    8.0_DP/9.0_DP * derrorP_GaussPt2**2 + &
                    5.0_DP/9.0_DP * derrorP_GaussPt3**2)
            end if ! 2-Point or 3-Point Gauss quadrature
            !   \int_t^{0}^t^{1} || p(t) - p_ref(t) ||^2_L2 dt
            ! calculated

          end if ! end special case treatment: iglobaltimestep .eq. 2


          ! Calculate
          !   \int_t^{n}^t^{n+1} || u(t) - u_ref(t) ||^2_L2 dt
          ! and
          !   \int_t^{n}^t^{n+1} || p(t) - p_ref(t) ||^2_L2 dt
          ! (n >= 1)

          ! Time and space discretisation error for velocity in interval [t^{n}, t^{n+1}]
          select case (ccallbackMethod)
          case (1)
            call evalTimeSpaceErrorQuadLagrange(cerrorType, 1, &
                 igaussQuadDegree, ffunction_TargetX, rcollection, doldTime, dtime)
            call evalTimeSpaceErrorQuadLagrange(cerrorType, 2, &
                 igaussQuadDegree, ffunction_TargetY, rcollection, doldTime, dtime)

          case (2)
            rcollection%IquickAccess(2) = 1  ! component
            call evalTimeSpaceErrorQuadLagrange(cerrorType, 1, &
                 igaussQuadDegree, fcalc_error, rcollection, doldTime, dtime)
            rcollection%IquickAccess(2) = 2  ! component
            call evalTimeSpaceErrorQuadLagrange(cerrorType, 2, &
                 igaussQuadDegree, fcalc_error, rcollection, doldTime, dtime)
          end select

          if (igaussQuadDegree .eq. 2) then
            derrorVel_GaussPt1 = sqrt(Derr(1)**2 + Derr(2)**2)
            derrorVel_GaussPt2 = sqrt(Derr2(1)**2 + Derr2(2)**2)

            ! Apply Gauss quadrature rule of degree n=2
            ! which is defined on the interval [-1,1]. Transforming it to [t^{n}, t^{n+1}]
            ! leads to a factor of
            !               1/2 (t^{n+1} - t^{n})
            ! due to the rules for integration by substitution. Weighing factors of
            ! summands for Gauss quadrature rule of degree n=2 are both 1.0.
            dtimeSpaceErrorVel = dtimeSpaceErrorVel + &
                 0.5_DP*(dtime - doldTime) * &
                 (1.0_DP * derrorVel_GaussPt1**2 + &
                  1.0_DP * derrorVel_GaussPt2**2)

          else ! 3-Point Gauss quadrature
            derrorVel_GaussPt1 = sqrt(Derr(1)**2 + Derr(2)**2)
            derrorVel_GaussPt2 = sqrt(Derr2(1)**2 + Derr2(2)**2)
            derrorVel_GaussPt3 = sqrt(Derr3(1)**2 + Derr3(2)**2)

            ! Apply Gauss quadrature rule of degree n=3
            ! which is defined on the interval [-1,1]. Transforming it to [t^{n}, t^{n+1}]
            ! leads to a factor of
            !               1/2 (t^{n+1} - t^{n})
            ! due to the rules for integration by substitution. Weighing factors of
            ! summands for Gauss quadrature rule of degree n=3 are 5/9, 8/9, 5/9.
            dtimeSpaceErrorVel = dtimeSpaceErrorVel + &
                 0.5_DP*(dtime - doldTime) * &
                 (5.0_DP/9.0_DP * derrorVel_GaussPt1**2 + &
                  8.0_DP/9.0_DP * derrorVel_GaussPt2**2 + &
                  5.0_DP/9.0_DP * derrorVel_GaussPt3**2)
          end if ! 2-Point or 3-Point Gauss quadrature
          !   \int_t^{0}^t^{n+1} || u(t) - u_ref(t) ||^2_L2 dt
          ! calculated.


          ! Time plus space discretisation error for pressure in interval [t^{n}, t^{n+1}]
          select case (ccallbackMethod)
          case (1)
            call evalTimeSpaceErrorQuadLagrange(cerrorType, 3, &
                 igaussQuadDegree, ffunction_TargetP, rcollection, doldTime, dtime)

          case (2)
            rcollection%IquickAccess(2) = 3  ! component
            call evalTimeSpaceErrorQuadLagrange(cerrorType, 3, &
                 igaussQuadDegree, fcalc_error, rcollection, doldTime, dtime)
          end select

          if (igaussQuadDegree .eq. 2) then
            derrorP_GaussPt1 = Derr(3)
            derrorP_GaussPt2 = Derr2(3)

            ! Apply Gauss quadrature rule of degree n=2
            ! which is defined on the interval [-1,1]. Transforming it to [t^{n}, t^{n+1}]
            ! leads to a factor of
            !               1/2 (t^{n+1} - t^{n})
            ! due to the rules for integration by substitution. Weighing factors of
            ! summands for Gauss quadrature rule of degree n=2 are both 1.0.
            dtimeSpaceErrorP = dtimeSpaceErrorP + &
                 0.5_DP*(dtime - doldTime) * &
                 (1.0_DP * derrorP_GaussPt1**2 + &
                  1.0_DP * derrorP_GaussPt2**2)

          else ! 3-Point Gauss quadrature
            derrorP_GaussPt1 = Derr(3)
            derrorP_GaussPt2 = Derr2(3)
            derrorP_GaussPt3 = Derr3(3)

            ! Apply Gauss quadrature rule of degree n=3
            ! which is defined on the interval [-1,1]. Transforming it to [t^{n}, t^{n+1}]
            ! leads to a factor of
            !               1/2 (t^{n+1} - t^{n})
            ! due to the rules for integration by substitution. Weighing factors of
            ! summands for Gauss quadrature rule of degree n=3 are 5/9, 8/9, 5/9.
            dtimeSpaceErrorP = dtimeSpaceErrorP + &
                 0.5_DP*(dtime - doldTime) * &
                 (5.0_DP/9.0_DP * derrorP_GaussPt1**2 + &
                  8.0_DP/9.0_DP * derrorP_GaussPt2**2 + &
                  5.0_DP/9.0_DP * derrorP_GaussPt3**2)
          end if ! 2-Point or 3-Point Gauss quadrature
          !   \int_t^{0}^t^{n+1} || p(t) - p_ref(t) ||^2_L2 dt
          ! calculated

          if (igaussQuadDegree .eq. 3) then
            call lsysbl_releaseVector (rvector_sol_GaussPt3)
          end if
          call lsysbl_releaseVector (rvector_sol_GaussPt2)
          call lsysbl_releaseVector (rvector_sol_GaussPt1)
          call lsysbl_releaseVector (rvector_auxSum)

        end if  ! (iglobaltimestep .ge. 2 .and. cignoreTimeStep .eq. 0)


      case (4,5)
        ! Time-interpolate the solution with a cubic Lagrange polynomial, then evaluate
        !    || sol(t) - sol_ref(t) ||^2_L2(0,T;\Omega)
        !    = sqrt( \int_{t^{start}}^{t^{end}} || sol(t) - sol_ref(t) ||^2_L2 dt )
        ! with 2-point or 3-point Gaussian quadrature rule.
        if (iglobaltimestep .ge. 3 .and. cignoreTimeStep .eq. 0) then

          ! Create temporary vectors
          call lsysbl_createVector (rsolution, rvector_auxSum, .false.)
          call lsysbl_createVector (rsolution, rvector_sol_GaussPt1, .false.)
          call lsysbl_createVector (rsolution, rvector_sol_GaussPt2, .false.)
          if (rpostprocessing%itimeSpaceDiscrErrorMethod .eq. 4) then
            igaussQuadDegree = 2
          else
            igaussQuadDegree = 3
            call lsysbl_createVector (rsolution, rvector_sol_GaussPt3, .false.)
          end if

          if (iglobaltimestep .eq. 3) then
            ! Special case treatment

            ! Why are the time intervals [t^{0}, t^{1}], [t^{1}, t^{2}] a special case?
            ! We need to make up for the fact of not having been able to calculate time
            ! discretisation error in time step 1 nor 2 due to missing data points for
            ! cubic Lagrange time-interpolation of the solution.

            ! velocity part 1:
            ! Calculate
            !   \int_t^{0}^t^{1} || u(t) - u_ref(t) ||^2_L2 dt
            ! and
            !   \int_t^{0}^t^{1} || p(t) - p_ref(t) ||^2_L2 dt
            !
            ! Time plus space discretisation error for velocity in interval [t^{0}, t^{1}]
            select case (ccallbackMethod)
            case (1)
              call evalTimeSpaceErrorCubicLagrange(cerrorType, 1, &
                   igaussQuadDegree, ffunction_TargetX, rcollection, &
                   doldestTime, dolderTime)
              call evalTimeSpaceErrorCubicLagrange(cerrorType, 2, &
                   igaussQuadDegree, ffunction_TargetY, rcollection, &
                   doldestTime, dolderTime)

            case (2)
              rcollection%IquickAccess(2) = 1  ! component
              call evalTimeSpaceErrorCubicLagrange(cerrorType, 1, &
                   igaussQuadDegree, fcalc_error, rcollection, &
                   doldestTime, dolderTime)
              rcollection%IquickAccess(2) = 2  ! component
              call evalTimeSpaceErrorCubicLagrange(cerrorType, 2, &
                   igaussQuadDegree, fcalc_error, rcollection, &
                   doldestTime, dolderTime)
            end select

            if (igaussQuadDegree .eq. 2) then
              derrorVel_GaussPt1 = sqrt(Derr(1)**2 + Derr(2)**2)
              derrorVel_GaussPt2 = sqrt(Derr2(1)**2 + Derr2(2)**2)

              ! Apply Gauss quadrature rule of degree n=2
              ! which is defined on the interval [-1,1]. Transforming it to [t^{0}, t^{1}]
              ! leads to a factor of
              !               1/2 (t^{1} - t^{0})
              ! due to the rules for integration by substitution. Weighing factors of
              ! summands for Gauss quadrature rule of degree n=2 are both 1.0.
              dtimeSpaceErrorVel = dtimeSpaceErrorVel + &
                   0.5_DP*(dolderTime - doldestTime) * &
                   (1.0_DP * derrorVel_GaussPt1**2 + &
                    1.0_DP * derrorVel_GaussPt2**2)

            else ! 3-Point Gauss quadrature
              derrorVel_GaussPt1 = sqrt(Derr(1)**2 + Derr(2)**2)
              derrorVel_GaussPt2 = sqrt(Derr2(1)**2 + Derr2(2)**2)
              derrorVel_GaussPt3 = sqrt(Derr3(1)**2 + Derr3(2)**2)

              ! Apply Gauss quadrature rule of degree n=3
              ! which is defined on the interval [-1,1]. Transforming it to [t^{0}, t^{1}]
              ! leads to a factor of
              !               1/2 (t^{1} - t^{0})
              ! due to the rules for integration by substitution. Weighing factors of
              ! summands for Gauss quadrature rule of degree n=3 are 5/9, 8/9, 5/9.
              dtimeSpaceErrorVel = dtimeSpaceErrorVel + &
                   0.5_DP*(dolderTime - doldestTime) * &
                   (5.0_DP/9.0_DP * derrorVel_GaussPt1**2 + &
                    8.0_DP/9.0_DP * derrorVel_GaussPt2**2 + &
                    5.0_DP/9.0_DP * derrorVel_GaussPt3**2)
            end if ! 2-Point or 3-Point Gauss quadrature
            !   \int_t^{0}^t^{1} || u(t) - u_ref(t) ||^2_L2 dt
            ! calculated.


            ! velocity part 2:
            ! Calculate
            !   \int_t^{1}^t^{2} || u(t) - u_ref(t) ||^2_L2 dt
            ! and
            !   \int_t^{1}^t^{2} || p(t) - p_ref(t) ||^2_L2 dt
            !
            ! Time plus space discretisation error for velocity in interval [t^{1}, t^{2}]
            select case (ccallbackMethod)
            case (1)
              call evalTimeSpaceErrorCubicLagrange(cerrorType, 1, &
                   igaussQuadDegree, ffunction_TargetX, rcollection, &
                   dolderTime, doldTime)
              call evalTimeSpaceErrorCubicLagrange(cerrorType, 2, &
                   igaussQuadDegree, ffunction_TargetY, rcollection, &
                   dolderTime, doldTime)

            case (2)
              rcollection%IquickAccess(2) = 1  ! component
              call evalTimeSpaceErrorCubicLagrange(cerrorType, 1, &
                   igaussQuadDegree, fcalc_error, rcollection, &
                   dolderTime, doldTime)
              rcollection%IquickAccess(2) = 2  ! component
              call evalTimeSpaceErrorCubicLagrange(cerrorType, 2, &
                   igaussQuadDegree, fcalc_error, rcollection, &
                   dolderTime, doldTime)
            end select

            if (igaussQuadDegree .eq. 2) then
              derrorVel_GaussPt1 = sqrt(Derr(1)**2 + Derr(2)**2)
              derrorVel_GaussPt2 = sqrt(Derr2(1)**2 + Derr2(2)**2)

              ! Apply Gauss quadrature rule of degree n=2
              ! which is defined on the interval [-1,1]. Transforming it to [t^{1}, t^{2}]
              ! leads to a factor of
              !               1/2 (t^{2} - t^{1})
              ! due to the rules for integration by substitution. Weighing factors of
              ! summands for Gauss quadrature rule of degree n=2 are both 1.0.
              dtimeSpaceErrorVel = dtimeSpaceErrorVel + &
                   0.5_DP*(doldTime - dolderTime) * &
                   (1.0_DP * derrorVel_GaussPt1**2 + &
                    1.0_DP * derrorVel_GaussPt2**2)

            else ! 3-Point Gauss quadrature
              derrorVel_GaussPt1 = sqrt(Derr(1)**2 + Derr(2)**2)
              derrorVel_GaussPt2 = sqrt(Derr2(1)**2 + Derr2(2)**2)
              derrorVel_GaussPt3 = sqrt(Derr3(1)**2 + Derr3(2)**2)

              ! Apply Gauss quadrature rule of degree n=3
              ! which is defined on the interval [-1,1]. Transforming it to [t^{1}, t^{2}]
              ! leads to a factor of
              !               1/2 (t^{2} - t^{1})
              ! due to the rules for integration by substitution. Weighing factors of
              ! summands for Gauss quadrature rule of degree n=3 are 5/9, 8/9, 5/9.
              dtimeSpaceErrorVel = dtimeSpaceErrorVel + &
                   0.5_DP*(doldTime - dolderTime) * &
                   (5.0_DP/9.0_DP * derrorVel_GaussPt1**2 + &
                    8.0_DP/9.0_DP * derrorVel_GaussPt2**2 + &
                    5.0_DP/9.0_DP * derrorVel_GaussPt3**2)
            end if ! 2-Point or 3-Point Gauss quadrature
            !   \int_t^{0}^t^{2} || u(t) - u_ref(t) ||^2_L2 dt
            ! calculated.


            ! pressure part 1:
            ! Time plus space discretisation error for pressure in interval [t^{0}, t^{1}]
            select case (ccallbackMethod)
            case (1)
              call evalTimeSpaceErrorCubicLagrange(cerrorType, 3, &
                   igaussQuadDegree, ffunction_TargetP, rcollection, &
                   doldestTime, dolderTime)

            case (2)
              rcollection%IquickAccess(2) = 3  ! component
              call evalTimeSpaceErrorCubicLagrange(cerrorType, 3, &
                   igaussQuadDegree, fcalc_error, rcollection, &
                   doldestTime, dolderTime)
            end select

            if (igaussQuadDegree .eq. 2) then
              derrorP_GaussPt1 = Derr(3)
              derrorP_GaussPt2 = Derr2(3)

              ! Apply Gauss quadrature rule of degree n=2
              ! which is defined on the interval [-1,1]. Transforming it to [t^{0}, t^{1}]
              ! leads to a factor of
              !               1/2 (t^{1} - t^{0})
              ! due to the rules for integration by substitution. Weighing factors of
              ! summands for Gauss quadrature rule of degree n=2 are both 1.0.
              dtimeSpaceErrorP = dtimeSpaceErrorP + &
                   0.5_DP*(dolderTime - doldestTime) * &
                   (1.0_DP * derrorP_GaussPt1**2 + &
                    1.0_DP * derrorP_GaussPt2**2)

            else ! 3-Point Gauss quadrature
              derrorP_GaussPt1 = Derr(3)
              derrorP_GaussPt2 = Derr2(3)
              derrorP_GaussPt3 = Derr3(3)

              ! Apply Gauss quadrature rule of degree n=3
              ! which is defined on the interval [-1,1]. Transforming it to [t^{0}, t^{1}]
              ! leads to a factor of
              !               1/2 (t^{1} - t^{0})
              ! due to the rules for integration by substitution. Weighing factors of
              ! summands for Gauss quadrature rule of degree n=3 are 5/9, 8/9, 5/9.
              dtimeSpaceErrorP = dtimeSpaceErrorP + &
                   0.5_DP*(dolderTime - doldestTime) * &
                   (5.0_DP/9.0_DP * derrorP_GaussPt1**2 + &
                    8.0_DP/9.0_DP * derrorP_GaussPt2**2 + &
                    5.0_DP/9.0_DP * derrorP_GaussPt3**2)
            end if ! 2-Point or 3-Point Gauss quadrature
            !   \int_t^{0}^t^{1} || p(t) - p_ref(t) ||^2_L2 dt
            ! calculated


            ! pressure part 2:
            ! Time plus space discretisation error for pressure in interval [t^{1}, t^{2}]
            select case (ccallbackMethod)
            case (1)
              call evalTimeSpaceErrorCubicLagrange(cerrorType, 3, &
                   igaussQuadDegree, ffunction_TargetP, rcollection, &
                   dolderTime, doldTime)

            case (2)
              rcollection%IquickAccess(2) = 3  ! component
              call evalTimeSpaceErrorCubicLagrange(cerrorType, 3, &
                   igaussQuadDegree, fcalc_error, rcollection, &
                   dolderTime, doldTime)
            end select

            if (igaussQuadDegree .eq. 2) then
              derrorP_GaussPt1 = Derr(3)
              derrorP_GaussPt2 = Derr2(3)

              ! Apply Gauss quadrature rule of degree n=2
              ! which is defined on the interval [-1,1]. Transforming it to [t^{1}, t^{2}]
              ! leads to a factor of
              !               1/2 (t^{2} - t^{1})
              ! due to the rules for integration by substitution. Weighing factors of
              ! summands for Gauss quadrature rule of degree n=2 are both 1.0.
              dtimeSpaceErrorP = dtimeSpaceErrorP + &
                   0.5_DP*(doldTime - dolderTime) * &
                   (1.0_DP * derrorP_GaussPt1**2 + &
                    1.0_DP * derrorP_GaussPt2**2)

            else ! 3-Point Gauss quadrature
              derrorP_GaussPt1 = Derr(3)
              derrorP_GaussPt2 = Derr2(3)
              derrorP_GaussPt3 = Derr3(3)

              ! Apply Gauss quadrature rule of degree n=3
              ! which is defined on the interval [-1,1]. Transforming it to [t^{1}, t^{2}]
              ! leads to a factor of
              !               1/2 (t^{2} - t^{1})
              ! due to the rules for integration by substitution. Weighing factors of
              ! summands for Gauss quadrature rule of degree n=3 are 5/9, 8/9, 5/9.
              dtimeSpaceErrorP = dtimeSpaceErrorP + &
                   0.5_DP*(doldTime - dolderTime) * &
                   (5.0_DP/9.0_DP * derrorP_GaussPt1**2 + &
                    8.0_DP/9.0_DP * derrorP_GaussPt2**2 + &
                    5.0_DP/9.0_DP * derrorP_GaussPt3**2)
            end if ! 2-Point or 3-Point Gauss quadrature
            !   \int_t^{0}^t^{2} || p(t) - p_ref(t) ||^2_L2 dt
            ! calculated
          end if ! end special case treatment: iglobaltimestep .eq. 3


          ! Calculate
          !   \int_t^{n}^t^{n+1} || u(t) - u_ref(t) ||^2_L2 dt
          ! and
          !   \int_t^{n}^t^{n+1} || p(t) - p_ref(t) ||^2_L2 dt
          ! (n >= 2)

          ! Time and space discretisation error for velocity in interval [t^{n}, t^{n+1}]
          select case (ccallbackMethod)
          case (1)
            call evalTimeSpaceErrorCubicLagrange(cerrorType, 1, &
                 igaussQuadDegree, ffunction_TargetX, rcollection, doldTime, dtime)
            call evalTimeSpaceErrorCubicLagrange(cerrorType, 2, &
                 igaussQuadDegree, ffunction_TargetY, rcollection, doldTime, dtime)

          case (2)
            rcollection%IquickAccess(2) = 1  ! component
            call evalTimeSpaceErrorCubicLagrange(cerrorType, 1, &
                 igaussQuadDegree, fcalc_error, rcollection, doldTime, dtime)
            rcollection%IquickAccess(2) = 2  ! component
            call evalTimeSpaceErrorCubicLagrange(cerrorType, 2, &
                 igaussQuadDegree, fcalc_error, rcollection, doldTime, dtime)
          end select

          if (igaussQuadDegree .eq. 2) then
            derrorVel_GaussPt1 = sqrt(Derr(1)**2 + Derr(2)**2)
            derrorVel_GaussPt2 = sqrt(Derr2(1)**2 + Derr2(2)**2)

            ! Apply Gauss quadrature rule of degree n=2
            ! which is defined on the interval [-1,1]. Transforming it to [t^{n}, t^{n+1}]
            ! leads to a factor of
            !               1/2 (t^{n+1} - t^{n})
            ! due to the rules for integration by substitution. Weighing factors of
            ! summands for Gauss quadrature rule of degree n=2 are both 1.0.
            dtimeSpaceErrorVel = dtimeSpaceErrorVel + &
                 0.5_DP*(dtime - doldTime) * &
                 (1.0_DP * derrorVel_GaussPt1**2 + &
                  1.0_DP * derrorVel_GaussPt2**2)

          else ! 3-Point Gauss quadrature
            derrorVel_GaussPt1 = sqrt(Derr(1)**2 + Derr(2)**2)
            derrorVel_GaussPt2 = sqrt(Derr2(1)**2 + Derr2(2)**2)
            derrorVel_GaussPt3 = sqrt(Derr3(1)**2 + Derr3(2)**2)

            ! Apply Gauss quadrature rule of degree n=3
            ! which is defined on the interval [-1,1]. Transforming it to [t^{n}, t^{n+1}]
            ! leads to a factor of
            !               1/2 (t^{n+1} - t^{n})
            ! due to the rules for integration by substitution. Weighing factors of
            ! summands for Gauss quadrature rule of degree n=3 are 5/9, 8/9, 5/9.
            dtimeSpaceErrorVel = dtimeSpaceErrorVel + &
                 0.5_DP*(dtime - doldTime) * &
                 (5.0_DP/9.0_DP * derrorVel_GaussPt1**2 + &
                  8.0_DP/9.0_DP * derrorVel_GaussPt2**2 + &
                  5.0_DP/9.0_DP * derrorVel_GaussPt3**2)
          end if ! 2-Point or 3-Point Gauss quadrature
          !   \int_t^{0}^t^{n+1} || u(t) - u_ref(t) ||^2_L2 dt
          ! calculated.


          ! Time plus space discretisation error for pressure in interval [t^{n}, t^{n+1}]
          select case (ccallbackMethod)
          case (1)
            call evalTimeSpaceErrorCubicLagrange(cerrorType, 3, &
                 igaussQuadDegree, ffunction_TargetP, rcollection, doldTime, dtime)

          case (2)
            rcollection%IquickAccess(2) = 3  ! component
            call evalTimeSpaceErrorCubicLagrange(cerrorType, 3, &
                 igaussQuadDegree, fcalc_error, rcollection, doldTime, dtime)
          end select

          if (igaussQuadDegree .eq. 2) then
            derrorP_GaussPt1 = Derr(3)
            derrorP_GaussPt2 = Derr2(3)

            ! Apply Gauss quadrature rule of degree n=2
            ! which is defined on the interval [-1,1]. Transforming it to [t^{n}, t^{n+1}]
            ! leads to a factor of
            !               1/2 (t^{n+1} - t^{n})
            ! due to the rules for integration by substitution. Weighing factors of
            ! summands for Gauss quadrature rule of degree n=2 are both 1.0.
            dtimeSpaceErrorP = dtimeSpaceErrorP + &
                 0.5_DP*(dtime - doldTime) * &
                 (1.0_DP * derrorP_GaussPt1**2 + &
                  1.0_DP * derrorP_GaussPt2**2)

          else ! 3-Point Gauss quadrature
            derrorP_GaussPt1 = Derr(3)
            derrorP_GaussPt2 = Derr2(3)
            derrorP_GaussPt3 = Derr3(3)

            ! Apply Gauss quadrature rule of degree n=3
            ! which is defined on the interval [-1,1]. Transforming it to [t^{n}, t^{n+1}]
            ! leads to a factor of
            !               1/2 (t^{n+1} - t^{n})
            ! due to the rules for integration by substitution. Weighing factors of
            ! summands for Gauss quadrature rule of degree n=3 are 5/9, 8/9, 5/9.
            dtimeSpaceErrorP = dtimeSpaceErrorP + &
                 0.5_DP*(dtime - doldTime) * &
                 (5.0_DP/9.0_DP * derrorP_GaussPt1**2 + &
                  8.0_DP/9.0_DP * derrorP_GaussPt2**2 + &
                  5.0_DP/9.0_DP * derrorP_GaussPt3**2)
          end if ! 2-Point or 3-Point Gauss quadrature
          !   \int_t^{0}^t^{n+1} || p(t) - p_ref(t) ||^2_L2 dt
          ! calculated

          if (igaussQuadDegree .eq. 3) then
            call lsysbl_releaseVector (rvector_sol_GaussPt3)
          end if
          call lsysbl_releaseVector (rvector_sol_GaussPt2)
          call lsysbl_releaseVector (rvector_sol_GaussPt1)
          call lsysbl_releaseVector (rvector_auxSum)

        end if  ! (iglobaltimestep .ge. 3 .and. cignoreTimeStep .eq. 0)

      end select  ! rpostprocessing%itimeSpaceDiscrErrorMethod

    end subroutine evalTimeSpaceError

    ! ---------------------------------------------------------------------------------

!<subroutine>

    subroutine evalTimeSpaceErrorQuadLagrange(cerrorType, iblock, igaussQuadDegree, &
         ffunctionReference, rcollection, dtimeStart, dtimeEnd)

!<description>
    ! Function that Lagrange (time-)interpolates the solutions from the last three time
    ! steps (t^{n-1}, t^{n}, t^{n+1} with t^{n+1} being the current one, i.e. the one from
    ! dtimebackup) with a quadratic polynomial, evaluates it in the Gauss points of the
    ! Gauss quadrature rule of degree n=2 or n=3 in time interval [timeStart, timeEnd],
    ! calculates the L2 or H1 error in these Gauss points against the reference solution`s
    ! result in these Gauss points and returns the value.
    ! The result is stored in Derr(iblock), Derr2(iblock) for Gauss quadrature of degree
    ! n=2 and Derr(iblock), Derr2(iblock) and Derr3(iblock) for Gauss quadrature of
    ! degree n=3.
!</description>

      ! Type of error to compute. Bitfield. This is a combination of the
      ! PPERR_xxxx-constants, which specifies what to compute.
      ! Example: PPERR_L2ERROR computes the $L_2$-error.
      integer, intent(in) :: cerrorType

      ! solution block to evaluate
      integer, intent(in) :: iblock

      ! degree of Gauss quadrature rule used for numerical integration
      ! (valid choices: 2 and 3)
      integer, intent(in) :: igaussQuadDegree

      ! A callback function that provides the analytical reference
      ! function to which the error should be computed.
      include '../../../kernel/Postprocessing/intf_refFunctionSc.inc'

      ! A pointer to a collection structure to provide additional
      ! information to the coefficient routine.
      type(t_collection), intent(inout) :: rcollection

      ! start point of time interval to evaluate
      real(DP), intent(in) :: dtimeStart

      ! end point of time interval to evaluate
      real(DP), intent(in) :: dtimeEnd
!</subroutine>


      ! local variable
      real(DP) :: dtimeBackup2


      select case (igaussQuadDegree)
      case (2)
        ! calculate the two Gauss points for Gauss quadrature rule of degree n=2
        ! in [dtimeStart, dtimeEnd]
        dtime_GaussPt1 = 0.5_DP * (dtimeStart + dtimeEnd &
                                - 1.0_DP/sqrt(3.0_DP) * (dtimeEnd - dtimeStart))
        dtime_GaussPt2 = 0.5_DP * (dtimeStart + dtimeEnd &
                                + 1.0_DP/sqrt(3.0_DP) * (dtimeEnd - dtimeStart))

        ! Quadratic Lagrange interpolate FEM solution in first Gauss point,
        !  rvector_sol(1st Gauss Pt,block) = \phi0(1st Gauss Pt) * solutionblock(t^{n-1})+
        !                                    \phi1(1st Gauss Pt) * solutionblock(t^{n}) +
        !                                    \phi2(1st Gauss Pt) * solutionblock(t^{n+1})
        call lsyssc_vectorLinearComb(&
             rpostprocessing%rolderSolution%RvectorBlock(iblock), &
             rpostprocessing%roldSolution%RvectorBlock(iblock), &
             timeQuadLagrangePhi_0(dtime_GaussPt1), &
             timeQuadLagrangePhi_1(dtime_GaussPt1), &
             rvector_auxSum%RvectorBlock(iblock))
        call lsyssc_vectorLinearComb(&
             rsolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeQuadLagrangePhi_2(dtime_GaussPt1), &
             1.0_DP, &
             rvector_sol_GaussPt1%RvectorBlock(iblock))
        ! ... and in second Gauss point
        call lsyssc_vectorLinearComb(&
             rpostprocessing%rolderSolution%RvectorBlock(iblock), &
             rpostprocessing%roldSolution%RvectorBlock(iblock), &
             timeQuadLagrangePhi_0(dtime_GaussPt2), &
             timeQuadLagrangePhi_1(dtime_GaussPt2), &
             rvector_auxSum%RvectorBlock(iblock))
        call lsyssc_vectorLinearComb(&
             rsolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeQuadLagrangePhi_2(dtime_GaussPt2), &
             1.0_DP, &
             rvector_sol_GaussPt2%RvectorBlock(iblock))

        ! Perform error analysis to calculate ||u - u_h||_{L^2} or ||u - u_h||_{H1} over
        ! [timeStart, timeEnd] in Gauss point s1...
        dtimeBackup2 = rcollection%Dquickaccess(1)
        ! manipulate rcollection%Dquickaccess(1) directly, do NOT use
        !     call cc_initCollectForAssembly (rproblem, rcollection)
        ! as that would also reset rcollection%Iquickaccess(1) to the value of
        ! rproblem%itimedependence whereas rcollection%Iquickaccess(1) is already being
        ! used to identify whether L2 or H1 errors are to be calculated in case reference
        ! solution is given in configuration files (evaluation), not hardcoded in
        ! ffunction_Target*.
        rcollection%Dquickaccess(1) = dtime_GaussPt1
        call pperr_scalar (cerrorType, Derr(iblock), &
             rvector_sol_GaussPt1%RvectorBlock(iblock), &
             ffunctionReference, rcollection, rcubatureInfo=rcubatureInfoUV)
        ! ... and s2
        rcollection%Dquickaccess(1) = dtime_GaussPt2
        call pperr_scalar (cerrorType, Derr2(iblock), &
             rvector_sol_GaussPt2%RvectorBlock(iblock), &
             ffunctionReference, rcollection, rcubatureInfo=rcubatureInfoUV)
        rcollection%Dquickaccess(1) = dtimeBackup2   ! restore value


      case (3)
        ! calculate the three Gauss points for Gauss quadrature rule of degree n=3
        ! in [dtimeStart, dtimeEnd]
        dtime_GaussPt1 = 0.5_DP * (dtimeStart + dtimeEnd &
                                - sqrt(3.0_DP/5.0_DP) * (dtimeEnd - dtimeStart))
        dtime_GaussPt2 = 0.5_DP * (dtimeStart + dtimeEnd)
        dtime_GaussPt3 = 0.5_DP * (dtimeStart + dtimeEnd &
                                + sqrt(3.0_DP/5.0_DP) * (dtimeEnd - dtimeStart))

        ! Cubic Lagrange interpolate FEM solution in first Gauss point,
        !  rvector_sol(1st Gauss Pt,block) = \phi0(1st Gauss Pt) * solutionblock(t^{n-1})+
        !                                    \phi1(1st Gauss Pt) * solutionblock(t^{n}) +
        !                                    \phi2(1st Gauss Pt) * solutionblock(t^{n+1})
        call lsyssc_vectorLinearComb(&
             rpostprocessing%rolderSolution%RvectorBlock(iblock), &
             rpostprocessing%roldSolution%RvectorBlock(iblock), &
             timeQuadLagrangePhi_0(dtime_GaussPt1), &
             timeQuadLagrangePhi_1(dtime_GaussPt1), &
             rvector_auxSum%RvectorBlock(iblock))
        call lsyssc_vectorLinearComb(&
             rsolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeQuadLagrangePhi_2(dtime_GaussPt1), &
             1.0_DP, &
             rvector_sol_GaussPt1%RvectorBlock(iblock))
        ! ... and in second Gauss point
        call lsyssc_vectorLinearComb(&
             rpostprocessing%rolderSolution%RvectorBlock(iblock), &
             rpostprocessing%roldSolution%RvectorBlock(iblock), &
             timeQuadLagrangePhi_0(dtime_GaussPt2), &
             timeQuadLagrangePhi_1(dtime_GaussPt2), &
             rvector_auxSum%RvectorBlock(iblock))
        call lsyssc_vectorLinearComb(&
             rsolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeQuadLagrangePhi_2(dtime_GaussPt2), &
             1.0_DP, &
             rvector_sol_GaussPt2%RvectorBlock(iblock))
        ! ... and in third Gauss point
        call lsyssc_vectorLinearComb(&
             rpostprocessing%rolderSolution%RvectorBlock(iblock), &
             rpostprocessing%roldSolution%RvectorBlock(iblock), &
             timeQuadLagrangePhi_0(dtime_GaussPt3), &
             timeQuadLagrangePhi_1(dtime_GaussPt3), &
             rvector_auxSum%RvectorBlock(iblock))
        call lsyssc_vectorLinearComb(&
             rsolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeQuadLagrangePhi_2(dtime_GaussPt3), &
             1.0_DP, &
             rvector_sol_GaussPt3%RvectorBlock(iblock))

        ! Perform error analysis to calculate ||u - u_h||_{L^2} or ||u - u_h||_{H1} over
        ! [timeStart, timeEnd] in Gauss point s1...
        dtimeBackup2 = rcollection%Dquickaccess(1)
        ! manipulate rcollection%Dquickaccess(1) directly, do NOT use
        !     call cc_initCollectForAssembly (rproblem, rcollection)
        ! as that would also reset rcollection%Iquickaccess(1) to the value of
        ! rproblem%itimedependence whereas rcollection%Iquickaccess(1) is already being
        ! used to identify whether L2 or H1 errors are to be calculated in case reference
        ! solution is given in configuration files (evaluation), not hardcoded in
        ! ffunction_Target*.
        rcollection%Dquickaccess(1) = dtime_GaussPt1
        call pperr_scalar (cerrorType, Derr(iblock), &
             rvector_sol_GaussPt1%RvectorBlock(iblock), &
             ffunctionReference, rcollection, rcubatureInfo=rcubatureInfoUV)
        ! ... and s2
        rcollection%Dquickaccess(1) = dtime_GaussPt2
        call pperr_scalar (cerrorType, Derr2(iblock), &
             rvector_sol_GaussPt2%RvectorBlock(iblock),&
             ffunctionReference, rcollection, rcubatureInfo=rcubatureInfoUV)
        ! ... and s3
        rcollection%Dquickaccess(1) = dtime_GaussPt3
        call pperr_scalar (cerrorType, Derr3(iblock), &
             rvector_sol_GaussPt3%RvectorBlock(iblock), &
             ffunctionReference, rcollection, rcubatureInfo=rcubatureInfoUV)
        rcollection%Dquickaccess(1) = dtimeBackup2   ! restore value


      case default
        call output_line ("Invalid degree for Gauss quadrature rule.", &
                          OU_CLASS_ERROR, OU_MODE_STD, "evalTimeSpaceErrorQuadLagrange")
        call output_line ("Valid choices for igaussQuadDegree: 2 or 3.", &
                          OU_CLASS_ERROR, OU_MODE_STD, "evalTimeSpaceErrorQuadLagrange")
        call sys_halt()

      end select

    end subroutine evalTimeSpaceErrorQuadLagrange

    ! ---------------------------------------------------------------------------------

!<subroutine>

    subroutine evalTimeSpaceErrorCubicLagrange(cerrorType, iblock, igaussQuadDegree, &
         ffunctionReference, rcollection, dtimeStart, dtimeEnd)

!<description>
    ! Function that Lagrange (time-)interpolates the solutions from the last four time
    ! steps (t^{n-2}, t^{n-1}, t^{n}, t^{n+1} with t^{n+1} being the current one, i.e. the
    ! one from dtimebackup) with a cubic polynomial, evaluates it in the the Gauss points
    ! of the Gauss quadrature rule of degree n=2 or n=3 in time interval [timeStart,
    ! timeEnd], calculates the L2 or H1 error in these Gauss points against the reference
    ! solution`s result in these Gauss points and returns the value.
    ! The result is stored in Derr(iblock), Derr2(iblock) for Gauss quadrature of degree
    ! n=2 and Derr(iblock), Derr2(iblock) and Derr3(iblock) for Gauss quadrature of
    ! degree n=3.
!</description>

      ! Type of error to compute. Bitfield. This is a combination of the
      ! PPERR_xxxx-constants, which specifies what to compute.
      ! Example: PPERR_L2ERROR computes the $L_2$-error.
      integer, intent(in) :: cerrorType

      ! solution block to evaluate
      integer, intent(in) :: iblock

      ! degree of Gauss quadrature rule used for numerical integration
      ! (valid choices: 2 and 3)
      integer, intent(in) :: igaussQuadDegree

      ! A callback function that provides the analytical reference
      ! function to which the error should be computed.
      include '../../../kernel/Postprocessing/intf_refFunctionSc.inc'

      ! A pointer to a collection structure to provide additional
      ! information to the coefficient routine.
      type(t_collection), intent(inout) :: rcollection

      ! start point of time interval to evaluate
      real(DP), intent(in) :: dtimeStart

      ! end point of time interval to evaluate
      real(DP), intent(in) :: dtimeEnd
!</subroutine>


      ! local variable
      real(DP) :: dtimeBackup2


      select case (igaussQuadDegree)
      case (2)
        ! calculate the two Gauss points for Gauss quadrature rule of degree n=2
        ! in [dtimeStart, dtimeEnd]
        dtime_GaussPt1 = 0.5_DP * (dtimeStart + dtimeEnd &
                                - 1.0_DP/sqrt(3.0_DP) * (dtimeEnd - dtimeStart))
        dtime_GaussPt2 = 0.5_DP * (dtimeStart + dtimeEnd &
                                + 1.0_DP/sqrt(3.0_DP) * (dtimeEnd - dtimeStart))

        ! Cubic Lagrange interpolate FEM solution in first Gauss point,
        !   rvector_sol(1st Gauss Pt,block) = \phi0(1st Gauss Pt) * solutionblock(t^{n-2}) +
        !                                     \phi1(1st Gauss Pt) * solutionblock(t^{n-1}) +
        !                                     \phi2(1st Gauss Pt) * solutionblock(t^{n})   +
        !                                     \phi3(1st Gauss Pt) * solutionblock(t^{n+1})
        call lsyssc_vectorLinearComb(&
             rpostprocessing%roldestSolution%RvectorBlock(iblock), &
             rpostprocessing%rolderSolution%RvectorBlock(iblock), &
             timeCubicLagrangePhi_0(dtime_GaussPt1), &
             timeCubicLagrangePhi_1(dtime_GaussPt1), &
             rvector_auxSum%RvectorBlock(iblock))
        call lsyssc_vectorLinearComb(&
             rpostprocessing%roldSolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeCubicLagrangePhi_2(dtime_GaussPt1), &
             1.0_DP)
        call lsyssc_vectorLinearComb(&
             rsolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeCubicLagrangePhi_3(dtime_GaussPt1), &
             1.0_DP, &
             rvector_sol_GaussPt1%RvectorBlock(iblock))
        ! ... and in second Gauss point
        call lsyssc_vectorLinearComb(&
             rpostprocessing%roldestSolution%RvectorBlock(iblock), &
             rpostprocessing%rolderSolution%RvectorBlock(iblock), &
             timeCubicLagrangePhi_0(dtime_GaussPt2), &
             timeCubicLagrangePhi_1(dtime_GaussPt2), &
             rvector_auxSum%RvectorBlock(iblock))
        call lsyssc_vectorLinearComb(&
             rpostprocessing%roldSolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeCubicLagrangePhi_2(dtime_GaussPt2), &
             1.0_DP)
        call lsyssc_vectorLinearComb(&
             rsolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeCubicLagrangePhi_3(dtime_GaussPt2), &
             1.0_DP, &
             rvector_sol_GaussPt2%RvectorBlock(iblock))

        ! Perform error analysis to calculate ||u - u_h||_{L^2} or ||u - u_h||_{H1} over
        ! [timeStart, timeEnd] in Gauss point s1...
        dtimeBackup2 = rcollection%Dquickaccess(1)
        ! manipulate rcollection%Dquickaccess(1) directly, do NOT use
        !     call cc_initCollectForAssembly (rproblem, rcollection)
        ! as that would also reset rcollection%Iquickaccess(1) to the value of
        ! rproblem%itimedependence whereas rcollection%Iquickaccess(1) is already being used
        ! to identify whether L2 or H1 errors are to be calculated in case reference
        ! solution is given in configuration files (evaluation), not hardcoded in
        ! ffunction_Target*.
        rcollection%Dquickaccess(1) = dtime_GaussPt1
        call pperr_scalar (cerrorType, Derr(iblock), rvector_sol_GaussPt1%RvectorBlock(iblock),&
             ffunctionReference, rcollection, rcubatureInfo=rcubatureInfoUV)
        ! ... and s2
        rcollection%Dquickaccess(1) = dtime_GaussPt2
        call pperr_scalar (cerrorType, Derr2(iblock),rvector_sol_GaussPt2%RvectorBlock(iblock),&
             ffunctionReference, rcollection, rcubatureInfo=rcubatureInfoUV)
        rcollection%Dquickaccess(1) = dtimeBackup2   ! restore value


      case (3)
        ! calculate the three Gauss points for Gauss quadrature rule of degree n=3
        ! in [dtimeStart, dtimeEnd]
        dtime_GaussPt1 = 0.5_DP * (dtimeStart + dtimeEnd &
                                - sqrt(3.0_DP/5.0_DP) * (dtimeEnd - dtimeStart))
        dtime_GaussPt2 = 0.5_DP * (dtimeStart + dtimeEnd)
        dtime_GaussPt3 = 0.5_DP * (dtimeStart + dtimeEnd &
                                + sqrt(3.0_DP/5.0_DP) * (dtimeEnd - dtimeStart))

        ! Cubic Lagrange interpolate FEM solution in first Gauss point,
        !   rvector_sol(1st Gauss Pt,block) = \phi0(1st Gauss Pt) * solutionblock(t^{n-2}) +
        !                                     \phi1(1st Gauss Pt) * solutionblock(t^{n-1}) +
        !                                     \phi2(1st Gauss Pt) * solutionblock(t^{n})   +
        !                                     \phi3(1st Gauss Pt) * solutionblock(t^{n+1})
        call lsyssc_vectorLinearComb(&
             rpostprocessing%roldestSolution%RvectorBlock(iblock), &
             rpostprocessing%rolderSolution%RvectorBlock(iblock), &
             timeCubicLagrangePhi_0(dtime_GaussPt1), &
             timeCubicLagrangePhi_1(dtime_GaussPt1), &
             rvector_auxSum%RvectorBlock(iblock))
        call lsyssc_vectorLinearComb(&
             rpostprocessing%roldSolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeCubicLagrangePhi_2(dtime_GaussPt1), &
             1.0_DP)
        call lsyssc_vectorLinearComb(&
             rsolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeCubicLagrangePhi_3(dtime_GaussPt1), &
             1.0_DP, &
             rvector_sol_GaussPt1%RvectorBlock(iblock))
        ! ... and in second Gauss point
        call lsyssc_vectorLinearComb(&
             rpostprocessing%roldestSolution%RvectorBlock(iblock), &
             rpostprocessing%rolderSolution%RvectorBlock(iblock), &
             timeCubicLagrangePhi_0(dtime_GaussPt2), &
             timeCubicLagrangePhi_1(dtime_GaussPt2), &
             rvector_auxSum%RvectorBlock(iblock))
        call lsyssc_vectorLinearComb(&
             rpostprocessing%roldSolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeCubicLagrangePhi_2(dtime_GaussPt2), &
             1.0_DP)
        call lsyssc_vectorLinearComb(&
             rsolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeCubicLagrangePhi_3(dtime_GaussPt2), &
             1.0_DP, &
             rvector_sol_GaussPt2%RvectorBlock(iblock))
        ! ... and in third Gauss point
        call lsyssc_vectorLinearComb(&
             rpostprocessing%roldestSolution%RvectorBlock(iblock), &
             rpostprocessing%rolderSolution%RvectorBlock(iblock), &
             timeCubicLagrangePhi_0(dtime_GaussPt3), &
             timeCubicLagrangePhi_1(dtime_GaussPt3), &
             rvector_auxSum%RvectorBlock(iblock))
        call lsyssc_vectorLinearComb(&
             rpostprocessing%roldSolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeCubicLagrangePhi_2(dtime_GaussPt3), &
             1.0_DP)
        call lsyssc_vectorLinearComb(&
             rsolution%RvectorBlock(iblock), &
             rvector_auxSum%RvectorBlock(iblock), &
             timeCubicLagrangePhi_3(dtime_GaussPt3), &
             1.0_DP, &
             rvector_sol_GaussPt3%RvectorBlock(iblock))

        ! Perform error analysis to calculate ||u - u_h||_{L^2} or ||u - u_h||_{H1} over
        ! [timeStart, timeEnd] in Gauss point s1...
        dtimeBackup2 = rcollection%Dquickaccess(1)
        ! manipulate rcollection%Dquickaccess(1) directly, do NOT use
        !     call cc_initCollectForAssembly (rproblem, rcollection)
        ! as that would also reset rcollection%Iquickaccess(1) to the value of
        ! rproblem%itimedependence whereas rcollection%Iquickaccess(1) is already being used
        ! to identify whether L2 or H1 errors are to be calculated in case reference
        ! solution is given in configuration files (evaluation), not hardcoded in
        ! ffunction_Target*.
        rcollection%Dquickaccess(1) = dtime_GaussPt1
        call pperr_scalar (cerrorType, Derr(iblock), rvector_sol_GaussPt1%RvectorBlock(iblock),&
             ffunctionReference, rcollection, rcubatureInfo=rcubatureInfoUV)
        ! ... and s2
        rcollection%Dquickaccess(1) = dtime_GaussPt2
        call pperr_scalar (cerrorType, Derr2(iblock),rvector_sol_GaussPt2%RvectorBlock(iblock),&
             ffunctionReference, rcollection, rcubatureInfo=rcubatureInfoUV)
        ! ... and s3
        rcollection%Dquickaccess(1) = dtime_GaussPt3
        call pperr_scalar (cerrorType, Derr3(iblock),rvector_sol_GaussPt3%RvectorBlock(iblock),&
             ffunctionReference, rcollection, rcubatureInfo=rcubatureInfoUV)
        rcollection%Dquickaccess(1) = dtimeBackup2   ! restore value


      case default
        call output_line ("Invalid degree for Gauss quadrature rule.", &
                          OU_CLASS_ERROR, OU_MODE_STD, "evalTimeSpaceErrorCubicLagrange")
        call output_line ("Valid choices for igaussQuadDegree: 2 or 3.", &
                          OU_CLASS_ERROR, OU_MODE_STD, "evalTimeSpaceErrorCubicLagrange")
        call sys_halt()

      end select

    end subroutine evalTimeSpaceErrorCubicLagrange

    ! ---------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------
    ! Function procedures to compute the basis function in time
    ! ---------------------------------------------------------------------------------

    ! These three routines evalute the 3 (time-)interpolation quadratic Lagrange basis
    ! polynomials \phi_0(t), \phi_1(t), \phi_2(t) on [-1,1]. The interpolated data points
    ! are (t^{n-1}, .), (t^{n}, .) and (t^{n+1}, .).
    ! Be sure to pass as argument a time value in [t^{n-1}, t^{n}].
    ! We are assuming time step sizes greater than zero (otherwise we would have divisions
    ! by zero in the following) which seems a safe assumption.
    function timeQuadLagrangePhi_0(t)
      real(DP) :: timeQuadLagrangePhi_0, t

      ! Lagrange basis polynomial l0 := (t-t1)/(t0-t1)*(t-t2)/(t0-t2)
      ! with t0, t1 and t2 being the three given interpolation points t^{n-1}, t^{n}, t^{n+1}
      timeQuadLagrangePhi_0 = &
           (t - doldTime)/(dolderTime - doldTime) * (t - dtime)/(dolderTime - dtime)

    end function timeQuadLagrangePhi_0

    ! ---------------------------------------------------------------------------------

    function timeQuadLagrangePhi_1(t)
      real(DP) :: timeQuadLagrangePhi_1, t

      ! Lagrange basis polynomial l1 := (t-t0)/(t1-t0)*(t-t2)/(t1-t2)
      ! with t0, t1 and t2 being the three given interpolation points t^{n-1}, t^{n}, t^{n+1}
      timeQuadLagrangePhi_1 = &
           (t - dolderTime)/(doldTime - dolderTime) * (t - dtime)/(doldTime - dtime)

    end function timeQuadLagrangePhi_1

    ! ---------------------------------------------------------------------------------

    function timeQuadLagrangePhi_2(t)
      real(DP) :: timeQuadLagrangePhi_2, t

      ! Lagrange basis polynomial l2 := (t-t0)/(t2-t0)*(t-t1)/(t2-t1)
      ! with t0, t1 and t2 being the three given interpolation points t^{n-1}, t^{n}, t^{n+1}
      timeQuadLagrangePhi_2 = &
           (t - dolderTime)/(dtime - dolderTime) * (t - doldTime)/(dtime - doldTime)

    end function timeQuadLagrangePhi_2

    ! ---------------------------------------------------------------------------------

    ! These four routines evalute the 4 (time-)interpolation Lagrange basis polynomials
    ! \phi_0(t), \phi_1(t), \phi_2(t), \phi_3(t) on [-1,1]. The interpolated data points are
    ! (t^{n-2}, .), (t^{n-1}, .), (t^{n}, .) and (t^{n+1}, .).
    ! Be sure to pass as argument a time value in [t^{n-2}, t^{n}].
    ! We are assuming time step sizes greater than zero (otherwise we would have divisions
    ! by zero in the following) which seems a safe assumption.
    function timeCubicLagrangePhi_0(t)
        real(DP) :: timeCubicLagrangePhi_0, t

        ! Lagrange basis polynomial l0 := (t-t1)/(t0-t1) * (t-t2)/(t0-t2) * (t-t3)/(t0-t3)
        ! with t0, t1, t2, t3 being the four given interpolation points t^{n-2}, t^{n-1},
        ! t^{n}, t^{n+1}
        timeCubicLagrangePhi_0 = &
             (t - dolderTime)/(doldestTime - dolderTime) * &
             (t - doldTime)/(doldestTime - doldTime) * &
             (t - dtime)/(doldestTime - dtime)

    end function timeCubicLagrangePhi_0

    ! ---------------------------------------------------------------------------------

    function timeCubicLagrangePhi_1(t)
        real(DP) :: timeCubicLagrangePhi_1, t

        ! Lagrange basis polynomial l1 := (t-t0)/(t1-t0) * (t-t2)/(t1-t2) * (t-t3)/(t1-t3)
        ! with t0, t1, t2, t3 being the four given interpolation points t^{n-2}, t^{n-1},
        ! t^{n}, t^{n+1}
        timeCubicLagrangePhi_1 = &
             (t - doldestTime)/(dolderTime - doldestTime) * &
             (t - doldTime)/(dolderTime - doldTime) * &
             (t - dtime)/(dolderTime - dtime)

    end function timeCubicLagrangePhi_1

    ! ---------------------------------------------------------------------------------

    function timeCubicLagrangePhi_2(t)
        real(DP) :: timeCubicLagrangePhi_2, t

        ! Lagrange basis polynomial l2 := (t-t0)/(t2-t0) * (t-t1)/(t2-t1) * (t-t3)/(t2-t3)
        ! with t0, t1, t2, t3 being the four given interpolation points t^{n-2}, t^{n-1},
        ! t^{n}, t^{n+1}
        timeCubicLagrangePhi_2 = &
             (t - doldestTime)/(doldTime - doldestTime) * &
             (t - dolderTime)/(doldTime - dolderTime) * &
             (t - dtime)/(doldTime - dtime)

    end function timeCubicLagrangePhi_2

    ! ---------------------------------------------------------------------------------

    function timeCubicLagrangePhi_3(t)
        real(DP) :: timeCubicLagrangePhi_3, t

        ! Lagrange basis polynomial l3 := (t-t0)/(t3-t0) * (t-t1)/(t3-t1) * (t-t2)/(t3-t2)
        ! with t0, t1, t2, t3 being the four given interpolation points t^{n-2}, t^{n-1},
        ! t^{n}, t^{n+1}
        timeCubicLagrangePhi_3 = &
             (t - doldestTime)/(dtime - doldestTime) * &
             (t - dolderTime)/(dtime - dolderTime) * &
             (t - doldTime)/(dtime - doldTime)

    end function timeCubicLagrangePhi_3

  end subroutine

! *****************************************************************

!<subroutine>

  subroutine ffunctionBDForcesVisco (cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset, &
                Dvalues,rcollection)

  use basicgeometry
  use triangulation
  use scalarpde
  use domainintegration
  use spatialdiscretisation
  use collection

!<description>
  ! Called in the postprocessing during the calculation of the body forces.
  ! Returns a nonconstant coefficient in the body force integral.
  !
  ! Wrapper to the ffunctionViscoModel callback routine which has
  ! nearly the same interface.
!</description>

!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(in)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation

  ! Number of elements, where the coefficients must be computed.
  integer, intent(in)                                         :: nelements

  ! Number of points per element, where the coefficients must be computed
  integer, intent(in)                                         :: npointsPerElement

  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(in)  :: Dpoints

  ! An array accepting the DOF`s on all elements trial in the trial space.
  ! DIMENSION(\#local DOF`s in trial space,Number of elements)
  integer, dimension(:,:), intent(in) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It is usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional
  ! information to the coefficient routine.
  type(t_collection), intent(inout), optional      :: rcollection

!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(out)                      :: Dvalues
!</output>

!</subroutine>

    ! cderivative will always be DER_FUNC here. Call ffunctionViscoModel
    ! to calculate the viscosity into Dvalues at first.
    call ffunctionViscoModel (1,rdiscretisation, &
              nelements,npointsPerElement,Dpoints, &
              IdofsTest,rdomainIntSubset, &
              Dvalues,rcollection)

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_calculateBodyForces (rsolution,dtime,rproblem)

!<description>
  ! Calculates body forces as configured in the .DAT file.
  ! The result is written to the standard output.
!</description>

!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(in), target :: rsolution

  ! Evaluation time. Must be set to 0.0 for stationary simulations.
  real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_collection) :: rcollection

    ! Forces on the object
    real(DP), dimension(NDIM2D) :: Dforces
    real(DP) :: dbdForcesCoeff1,dbdForcesCoeff2
    type(t_boundaryRegion) :: rregion
    integer :: cformulation
    integer :: ibodyForcesFormulation,icalcBodyForces,ibodyForcesBdComponent

    integer :: iwriteBodyForces
    character(len=SYS_STRLEN) :: sfilenameBodyForces
    character(len=SYS_STRLEN) :: stemp
    integer :: iunit
    integer :: cflag
    logical :: bfileExists
    type(t_vectorScalar) :: rcharfct

    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "icalcBodyForces", icalcBodyForces, 1)
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "ibodyForcesFormulation", ibodyForcesFormulation, -1)
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "ibodyForcesBdComponent", ibodyForcesBdComponent, 2)
    call parlst_getvalue_double (rproblem%rparamList, "CC-POSTPROCESSING", &
        "dbdForcesCoeff1", dbdForcesCoeff1, rproblem%rphysics%dnu)
    call parlst_getvalue_double (rproblem%rparamList, "CC-POSTPROCESSING", &
        "dbdForcesCoeff2", dbdForcesCoeff2, 0.1_DP * 0.2_DP**2)

    ! Probably cancel the calculation
    if (icalcBodyForces .eq. 0) return

    ! Information about writing body forces to a file.
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "IWRITEBODYFORCES", iwriteBodyForces, 0)
    call parlst_getvalue_string (rproblem%rparamList, "CC-POSTPROCESSING", &
        "SFILENAMEBODYFORCES", sfilenameBodyForces, "", bdequote=.true.)
    if (sfilenameBodyForces .eq. "") &
      iwriteBodyForces = 0

    ! When writing to a file is enabled, delete the file in the first timestep.
    cflag = SYS_APPEND
    if (rproblem%rtimedependence%itimeStep .eq. 0) cflag = SYS_REPLACE


    ! If we have a uniform discretisation, calculate the body forces on the
    ! 2nd boundary component - if it exists.
    if ((rsolution%p_rblockDiscr%RspatialDiscr(1)% &
         ccomplexity .eq. SPDISC_UNIFORM) .and. &
        (boundary_igetNBoundComp(rproblem%rboundary) .ge. ibodyForcesBdComponent)) then

      ! Calculate drag-/lift coefficients on the 2nd boundary component.
      ! This is for the benchmark channel!
      call boundary_createRegion (rproblem%rboundary, &
          ibodyForcesBdComponent, 0, rregion)
      rregion%iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND

      select case (icalcBodyForces)
      case (1)
        ! Old implementation:
        call ppns2D_bdforces_uniform (rsolution,rregion,Dforces,CUB_G4_1D,&
            dbdForcesCoeff1,dbdForcesCoeff2)

      case (2)
        ! Extended calculation method.
        !
        ! Select the tensor formulation to use.
        select case (ibodyForcesFormulation)
        case (0)
          cformulation = PPNAVST_GRADIENTTENSOR_SIMPLE
        case (1)
          cformulation = PPNAVST_GRADIENTTENSOR
        case (2)
          cformulation = PPNAVST_DEFORMATIONTENSOR
        case default
          ! Automatic mode.
          ! Use the deformation tensor formulation for the forces if
          ! we are in the deformation tensor formulation
          cformulation = PPNAVST_GRADIENTTENSOR_SIMPLE
          if (rproblem%rphysics%isubequation .eq. 1) &
            cformulation = PPNAVST_DEFORMATIONTENSOR
        end select

        ! Prepare the collection. The "next" collection points to the user defined
        ! collection.
        call ccmva_prepareViscoAssembly (rproblem,rproblem%rphysics,&
            rcollection,rsolution,rproblem%rcollection)

        if (rproblem%rphysics%cviscoModel .eq. 0) then
          call ppns2D_bdforces_line (rsolution,rregion,Dforces,CUB_G4_1D,&
              dbdForcesCoeff1,dbdForcesCoeff2,cformulation)
        else
          call ppns2D_bdforces_line (rsolution,rregion,Dforces,CUB_G4_1D,&
              dbdForcesCoeff1,dbdForcesCoeff2,cformulation,ffunctionBDForcesVisco,&
              rcollection,ntempArrays=5)
        end if

      case (3)
        ! Old implementation:
        call ppns2D_bdforces_uniform (rsolution,rregion,Dforces,CUB_G4_2D,&
            dbdForcesCoeff1,dbdForcesCoeff2)

      case (4)

        ! Create a characteristic function of the boundary segment
        call lsyssc_createVector (&
            rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation%RspatialDiscr(1),&
            rcharfct,.true.)

        call anprj_charFctRealBdComp (rregion,rcharfct)

        ! Select the tensor formulation to use.
        select case (ibodyForcesFormulation)
        case (0)
          cformulation = PPNAVST_GRADIENTTENSOR_SIMPLE
        case (1)
          cformulation = PPNAVST_GRADIENTTENSOR
        case (2)
          cformulation = PPNAVST_DEFORMATIONTENSOR
        case default
          ! Automatic mode.
          ! Use the deformation tensor formulation for the forces if
          ! we are in the deformation tensor formulation
          cformulation = PPNAVST_GRADIENTTENSOR_SIMPLE
          if (rproblem%rphysics%isubequation .eq. 1) &
            cformulation = PPNAVST_DEFORMATIONTENSOR
        end select

        ! Prepare a collection structure in the form necessary for
        ! the computation of a nonconstant viscosity.
        !
        ! Prepare the collection. The "next" collection points to the user defined
        ! collection.
        call ccmva_prepareViscoAssembly (rproblem,rproblem%rphysics,&
            rcollection,rsolution,rproblem%rcollection)

        if (rproblem%rphysics%cviscoModel .eq. 0) then
          call ppns2D_bdforces_vol(rsolution,rcharfct,Dforces,&
              dbdForcesCoeff1,dbdForcesCoeff2,cformulation)
        else
          call ppns2D_bdforces_vol(rsolution,rcharfct,Dforces,&
              dbdForcesCoeff1,dbdForcesCoeff2,cformulation,ffunctionBDForcesVisco,&
              rcollection,ntempArrays=5)
        end if

        call lsyssc_releaseVector(rcharfct)

      end select

      call output_lbrk()
      call output_line ("Body forces")
      call output_line ("-----------")
      call output_line ("Body forces real bd., bdc/horiz/vert",&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
      call output_line (" "//trim(sys_siL(ibodyForcesBdComponent,10)) // " / " &
          //trim(sys_sdEP(Dforces(1),15,6)) // " / "&
          //trim(sys_sdEP(Dforces(2),15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )

      if (iwriteBodyForces .ne. 0) then
        ! Write the result to a text file.
        ! Format: timestep current-time value
        call io_openFileForWriting(sfilenameBodyForces, iunit, &
            cflag, bfileExists,.true.)
        if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
          ! Write a headline
          write (iunit,"(A)") "# timestep time bdc horiz vert"
        end if
        stemp = trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // " " &
            // trim(sys_sdEL(dtime,10)) // " " &
            // trim(sys_siL(ibodyForcesBdComponent,10)) // " " &
            // trim(sys_sdEL(Dforces(1),10)) // " "&
            // trim(sys_sdEL(Dforces(2),10))
        write (iunit,"(A)") trim(stemp)
        close (iunit)
      end if

    end if

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_calculateDivergence (rsolution,rproblem)

!<description>
  ! Calculates the divergence of a solution.
  ! The result is written to the standard output.
!</description>

!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(in) :: rsolution
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: ieltype,icuberror
    character(LEN=SYS_STRLEN) :: stemp
    real(DP) :: ddivergence,ddivergenceL2norm
    type(t_fev2Vectors) :: revalVectors
    type(t_vectorScalar), target :: rtempVector
    type(t_scalarCubatureInfo) :: rcubatureInfo

    call spdiscr_getElemGroupInfo (rsolution%p_rblockDiscr%RspatialDiscr(1),1,ieltype)

    ! -------------------------------------------------------------
    ! Calculate the l2-norm of the divergence vector -- which
    ! is approximately the divergence in the points that represent
    ! the degrees of freedom of the FE space.
    !
    ! Note that this operation only makes sense for Langrangian-
    ! type finite elements (whose DOFs correspond to point values).
    ! -------------------------------------------------------------

    ddivergence = 0.0_DP

    ! Create a temporary vector
    call lsyssc_createVector (rsolution%RvectorBlock(3)%p_rspatialDiscr,&
        rtempVector,.true.)

    ! Calculate divergence = D1 u1 + D2 u2
    call lsyssc_matVec (&
        rproblem%RlevelInfo(rproblem%nlmax)%rasmTempl%rmatrixD1, rsolution%RvectorBlock(1), &
        rtempVector, 1.0_DP, 0.0_DP)
    call lsyssc_matVec (&
        rproblem%RlevelInfo(rproblem%nlmax)%rasmTempl%rmatrixD2, rsolution%RvectorBlock(2), &
        rtempVector, 1.0_DP, 1.0_DP)

    ddivergence = lsyssc_vectorNorm(rtempVector,LINALG_NORML2)

    call lsyssc_releaseVector (rtempVector)

    ! -------------------------------------------------------------
    ! Calculate the L2-norm of the divergence. This is a
    ! pointwise calculation.
    ! -------------------------------------------------------------

    ! Set up the cubature
    call parlst_getvalue_string (rproblem%rparamList,"CC-POSTPROCESSING",&
                                "scubError",stemp,"")
    if (stemp .eq. "") then
      call parlst_getvalue_int (rproblem%rparamList,"CC-POSTPROCESSING",&
                                "icubError",icubError,int(CUB_GEN_AUTO,I32))
    else
      icubError = cub_igetID(stemp)
    end if

    ! Create an cubature info structure which contains our cubature rule
    call spdiscr_createDefCubStructure(&
        rsolution%RvectorBlock(1)%p_rspatialDiscr,rcubatureInfo,int(icubError,I32))

    ! Calculate the divergence via block assembly methods.
    call fev2_addVectorFieldToEvalList (revalVectors,1,&
        rsolution%RvectorBlock(1),rsolution%RvectorBlock(2))

    call bma_buildIntegral (ddivergenceL2norm,BMA_CALC_STANDARD,&
        bma_fcalc_divergenceL2norm, revalVectors=revalVectors,&
        rcubatureInfo=rcubatureInfo)

    call fev2_releaseVectorList (revalVectors)

    ! Release cubature
    call spdiscr_releaseCubStructure(rcubatureInfo)

    ! Taking the square root gives the L2 norm
    ddivergenceL2norm = sqrt(ddivergenceL2norm)

    ! -------------------------------------------------------------
    ! Print
    ! -------------------------------------------------------------

    call output_lbrk()
    call output_line ("Divergence")
    call output_line ("----------")
    call output_line ("||def(div u)||_l2 = " &
        //trim(sys_sdEP(ddivergence,15,6)),coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
    call output_line ("||div u||_L2      = " &
        //trim(sys_sdEP(ddivergenceL2norm,15,6)),coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_evaluatePoints (rsolution,dtime,rproblem)

!<description>
  ! Evaluates the solution in a number of points as configured in the DAT file.
!</description>

!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(in), target :: rsolution

  ! Solution time. =0 for stationary simulations.
  real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: npoints,i,iunit,iwritePointValues
    integer :: iderType, cflag
    logical :: bfileExists
    real(dp), dimension(:), allocatable :: Dvalues
    real(dp), dimension(:,:), allocatable :: Dcoords
    integer, dimension(:), allocatable :: Itypes
    integer, dimension(:), allocatable :: Ider
    character(LEN=SYS_STRLEN) :: sparam
    character(LEN=SYS_STRLEN) :: sstr,sfilenamePointValues,stemp
    character(LEN=10), dimension(3,3), parameter :: Sfctnames = reshape (&
      (/ "       u1 ","     u1_x ","     u1_y " , &
         "       u2 ","     u2_x ","     u2_y " , &
         "        p ","      p_x ","      p_y " /) ,&
       (/ 3,3 /) )

    ! Get the number of points to evaluate
    npoints = parlst_querysubstrings (rproblem%rparamList, "CC-POSTPROCESSING", &
        "CEVALUATEPOINTVALUES")

    if (npoints .eq. 0) return

    ! Allocate memory for the values
    allocate(Dvalues(npoints))
    allocate(Dcoords(NDIM2D,npoints))
    allocate(Itypes(npoints))
    allocate(Ider(npoints))

    ! Read the points
    do i=1,npoints
      call parlst_getvalue_string (rproblem%rparamList, "CC-POSTPROCESSING", &
          "CEVALUATEPOINTVALUES", sparam, "", i)
      read (sparam,*) Dcoords(1,i),Dcoords(2,i),Itypes(i),Ider(i)
    end do

    ! Evaluate the function in these points.
    do i=1,npoints
      select case (Ider(i))
      case (0)
        iderType = DER_FUNC2D
      case (1)
        iderType = DER_DERIV2D_X
      case (2)
        iderType = DER_DERIV2D_Y
      case default
        iderType = DER_FUNC2D
      end select

      call fevl_evaluate (iderType, Dvalues(i:i), rsolution%RvectorBlock(Itypes(i)), &
          Dcoords(:,i:i),cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
    end do

    ! Print the values to the terminal
    call output_lbrk()
    call output_line ("Point values")
    call output_line ("------------")
    do i=1,npoints
      write (sstr,"(A10,A,F9.4,A,F9.4,A,E16.10)") Sfctnames(1+Ider(i),Itypes(i)),&
          "(",Dcoords(1,i),",",Dcoords(2,i),") = ",Dvalues(i)
      call output_line(trim(sstr),coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
    end do

    ! Get information about writing the stuff into a DAT file.
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "IWRITEPOINTVALUES", iwritePointValues, 0)
    call parlst_getvalue_string (rproblem%rparamList, "CC-POSTPROCESSING", &
        "SFILENAMEPOINTVALUES", sfilenamePointValues, """""", bdequote=.true.)
    if (sfilenamePointValues .eq. "") iwritePointValues = 0

    ! When writing to a file is enabled, delete the file in the first timestep.
    cflag = SYS_APPEND
    if (rproblem%rtimedependence%itimeStep .eq. 0) cflag = SYS_REPLACE

    if (iwritePointValues .ne. 0) then
      ! Write the result to a text file.
      ! Format: timestep current-time value value value ...
      call io_openFileForWriting(sfilenamePointValues, iunit, &
          cflag, bfileExists,.true.)
      if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
        ! Write a headline
        write (iunit,"(A)") &
          "# timestep time x y type deriv value x y type deriv value ..."
      end if
      stemp = &
          trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // " " // &
          trim(sys_sdEL(dtime,10))
      write (iunit,ADVANCE="NO",FMT="(A)") trim(stemp)
      do i=1,npoints
        stemp = " " //&
            trim(sys_sdEL(Dcoords(1,i),5)) // " " // &
            trim(sys_sdEL(Dcoords(2,i),5)) // " " // &
            trim(sys_siL(Itypes(i),2)) // " " // &
            trim(sys_siL(Ider(i),2)) // " " // &
            trim(sys_sdEL(Dvalues(i),10))
        write (iunit,ADVANCE="NO",FMT="(A)") trim(stemp)
      end do
      write (iunit,ADVANCE="YES",FMT="(A)") ""
      close (iunit)
    end if

    deallocate(Ider)
    deallocate(Itypes)
    deallocate(Dcoords)
    deallocate(Dvalues)

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_evaluateFlux (rsolution,dtime,rproblem)

!<description>
  ! Evaluates the flux thropugh a set of lines as configured in the DAT file.
!</description>

!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(in), target :: rsolution

  ! Simulation time. =0 for stationary simulations
  real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: nlines,i,iunit,iwriteFluxValues
    integer :: cflag
    logical :: bfileExists
    real(dp), dimension(:), allocatable :: Dvalues
    real(dp), dimension(:,:,:), allocatable :: Dcoords
    character(LEN=SYS_STRLEN) :: sparam
    character(LEN=SYS_STRLEN) :: sstr,sfilenameFluxValues,stemp

    ! Get the number of points to evaluate
    nlines = parlst_querysubstrings (rproblem%rparamList, "CC-POSTPROCESSING", &
        "CEVALUATEFLUXVALUES")

    if (nlines .eq. 0) return

    ! Allocate memory for the values
    allocate(Dvalues(nlines))
    allocate(Dcoords(NDIM2D,nlines,2))

    ! Read the points
    do i=1,nlines
      call parlst_getvalue_string (rproblem%rparamList, "CC-POSTPROCESSING", &
          "CEVALUATEFLUXVALUES", sparam, "", i)
      read (sparam,*) Dcoords(1,i,1),Dcoords(2,i,1),Dcoords(1,i,2),Dcoords(2,i,2)
    end do

    ! Calculate the flux
    do i=1,nlines
      call ppns2D_calcFluxThroughLine (rsolution,Dcoords(1:2,i,1),Dcoords(1:2,i,2),Dvalues(i))
    end do

    ! Print the values to the terminal
    call output_lbrk()
    call output_line ("Flux values")
    call output_line ("-----------")
    do i=1,nlines
      write (sstr,"(A,F9.4,A,F9.4,A,F9.4,A,F9.4,A,E16.10)") "flux (",&
          Dcoords(1,i,1),",",Dcoords(2,i,1),")->(",&
          Dcoords(1,i,2),",",Dcoords(2,i,2),") = ",Dvalues(i)
      call output_line(trim(sstr),coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
    end do

    ! Get information about writing the stuff into a DAT file.
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "IWRITEFLUXVALUES", iwriteFluxValues, 0)
    call parlst_getvalue_string (rproblem%rparamList, "CC-POSTPROCESSING", &
        "SFILENAMEFLUXVALUES", sfilenameFluxValues, """""",bdequote=.true.)
    if (sfilenameFluxValues .eq. "") iwriteFluxValues = 0

    ! When writing to a file is enabled, delete the file in the first timestep.
    cflag = SYS_APPEND
    if (rproblem%rtimedependence%itimeStep .eq. 0) cflag = SYS_REPLACE

    if (iwriteFluxValues .ne. 0) then
      ! Write the result to a text file.
      ! Format: timestep current-time value value value ...
      call io_openFileForWriting(sfilenameFluxValues, iunit, &
          cflag, bfileExists,.true.)
      if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
        ! Write a headline
        write (iunit,"(A)") &
          "# timestep time x1 y1 x2 y2 value1 x1 y1 x2 y2 value2 ..."
      end if
      stemp = &
          trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // " " // &
          trim(sys_sdEL(dtime,10))
      write (iunit,ADVANCE="NO",FMT="(A)") trim(stemp)
      do i=1,nlines
        stemp = " " //&
            trim(sys_sdEL(Dcoords(1,i,1),5)) // " " // &
            trim(sys_sdEL(Dcoords(2,i,1),5)) // " " // &
            trim(sys_sdEL(Dcoords(1,i,2),5)) // " " // &
            trim(sys_sdEL(Dcoords(2,i,2),5)) // " " // &
            trim(sys_sdEL(Dvalues(i),10))
        write (iunit,ADVANCE="NO",FMT="(A)") trim(stemp)
      end do
      write (iunit,ADVANCE="YES",FMT="(A)") ""
      close (iunit)
    end if

    deallocate(Dcoords)
    deallocate(Dvalues)

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_writeUCD (rpostprocessing,rvector,rproblem,dtime)

!<description>
  ! Writes an UCD postprocessing file as configured in the DAT file.
  ! (-> GMV, AVS, Paraview,...)
!</description>

!<input>
  ! Solution vector.
  type(t_vectorBlock), intent(in) :: rvector

  ! OPTIONAL: Simulation time.
  ! Must be ommitted in stationary simulations.
  real(DP), intent(in), optional :: dtime
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(inout), target :: rproblem

  ! Postprocessing structure. Must have been initialised prior
  ! to calling this routine.
  ! The time stamp of the last written out GMV is updated.
  type(t_c2d2postprocessing), intent(inout) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! local variables

    ! We need some more variables for postprocessing - i.e. writing
    ! a GMV file.
    real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2

    ! A pointer to the triangulation.
    type(t_triangulation), pointer :: p_rtriangulation

    ! A vector accepting Q1 data
    type(t_vectorBlock) :: rprjVector

    ! A discretisation structure for Q1
    type(t_blockDiscretisation) :: rprjDiscretisation

    ! A dynamic level information structure containing the BC"s.
    type(t_dynamicLevelInfo), target :: rdynamicInfo

    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport

    ! Backup of current simulation time.
    real(DP) :: dtimebackup

    integer :: ioutputUCD,ilevelUCD
    integer(I32) :: ieltype
    type(t_collection) :: rcollection

    character(SYS_STRLEN) :: sfile,sfilename


    ! Type of output:
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
                              "IOUTPUTUCD", ioutputUCD, 0)
    if (ioutputUCD .eq. 0) return

    ! Level of output:
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
                              "ILEVELUCD", ilevelUCD, 0)
    if (ilevelUCD .le. 0) then
      ilevelUCD = rproblem%NLMAX+ilevelUCD
    end if

    ilevelUCD = min(rproblem%NLMAX,max(rproblem%NLMIN,ilevelUCD))

    ! The solution vector is probably not in the way, GMV likes it!
    ! GMV for example does not understand Q1~ vectors!
    ! Therefore, we first have to convert the vector to a form that
    ! GMV understands.
    ! GMV understands only Q1/P1 solutions! So the task is now to
    ! create a Q1/P1 solution from rvector and write that out.
    !
    ! For this purpose, first create a "derived" simple discretisation
    ! structure based on Q1/P1 by copying the main guiding block
    ! discretisation structure and modifying the discretisation
    ! structures of the two velocity subvectors:

    call spdiscr_duplicateBlockDiscr(rvector%p_rblockDiscr,rprjDiscretisation)

    call spdiscr_deriveDiscr_triquad (rvector%p_rblockDiscr%RspatialDiscr(1), &
                 EL_P1, EL_Q1, rprjDiscretisation%RspatialDiscr(1))

    call spdiscr_deriveDiscr_triquad (rvector%p_rblockDiscr%RspatialDiscr(2), &
                 EL_P1, EL_Q1, rprjDiscretisation%RspatialDiscr(2))

    call spdiscr_deriveDiscr_triquad (rvector%p_rblockDiscr%RspatialDiscr(3), &
                 EL_P0, EL_Q0, rprjDiscretisation%RspatialDiscr(3))

    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVector (rprjDiscretisation,rprjVector,.false.)

    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (rvector,rprjVector)

    if (present(dtime)) then
      ! Only for the postprocessing, switch to time dtime.
      dtimebackup = rproblem%rtimedependence%dtime
      rproblem%rtimedependence%dtime = dtime
    end if

    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

    ! Initialise the dynamic level information structure
    call cc_initDynamicLevelInfo (rdynamicInfo)

    ! Discretise the boundary conditions according to the Q1/Q1/Q0
    ! discretisation for implementing them into a solution vector.
    call cc_assembleBDconditions (rproblem,rprjDiscretisation,&
        rdynamicInfo,rproblem%rcollection,.true.)

    ! The same way, discretise boundary conditions of fictitious boundary components.
    call cc_assembleFBDconditions (rproblem,rprjDiscretisation,&
        rdynamicInfo,rproblem%rcollection)

    ! Filter the solution vector to implement discrete BC`s.
    call vecfil_discreteBCsol (rprjVector,rdynamicInfo%rdiscreteBC)
    call vecfil_discreteFBCsol (rprjVector,rdynamicInfo%rdiscreteFBC)

    ! Basic filename
    call parlst_getvalue_string (rproblem%rparamList, "CC-POSTPROCESSING", &
                                 "SFILENAMEUCD", sfilename, "", bdequote=.true.)

    ! Create the actual filename
    sfile = trim(adjustl(sfilename))//"."//sys_si0(rpostprocessing%inextFileSuffixUCD,5)

    ! Now we have a Q1/Q1/Q0 solution in rprjVector -- on the level NLMAX.
    ! The next step is to project it down to level ilevelUCD.
    ! Due to the fact that solutions are usually 2-level-ordered,
    ! this can be shortened by taking only the first NVT vertices
    ! of the solution vector!

    ! From the attached discretisation, get the underlying triangulation
    ! of that level
    p_rtriangulation => rproblem%RlevelInfo(ilevelUCD)%rtriangulation

    ! Start UCD export to GMV file:
    call output_lbrk ()
    call output_line ("Writing visualisation file: "//sfile)

    select case (ioutputUCD)
    case (1)
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)

    case (2)
      call ucd_startAVS (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)

    case (3)
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)

    case (5)
      call ucd_startBGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)

    case default
      call output_line ("Invalid visualisation output type.", &
                        OU_CLASS_ERROR,OU_MODE_STD,"cc_writeUCD")
      call sys_halt()
    end select

    ! Is there a simulation time?
    if (present(dtime)) &
      call ucd_setSimulationTime (rexport,dtime)

    ! Write the configuration of the application as comment block
    ! to the output file.
    call ucd_addCommentLine (rexport,"Configuration:")
    call ucd_addCommentLine (rexport,"---------------")
    call ucd_addParameterList (rexport,rproblem%rparamList)
    call ucd_addCommentLine (rexport,"---------------")

    ! Get the velocity field
    call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)

    ! Moving frame velocity subtraction deactivated, gives pictures
    ! that can hardly be interpreted.

!    ! Is the moving-frame formulatino active?
!    call parlst_getvalue_int (rproblem%rparamList,"CC-DISCRETISATION",&
!        "imovingFrame",imovingFrame,0)
!
!    ! In this case, the postprocessing data must be modified by the
!    ! moving frame velocity.
!    if (imovingFrame .ne. 0) then
!
!      ! Get the velocity and acceleration from the callback routine.
!      call getMovingFrameVelocity (Dvelocity,Dacceleration,rproblem%rcollection)
!
!      ! Subtract the moving frame velocity from the postprocessing
!      ! data.
!      call lsyssc_addConstant (rprjVector%RvectorBlock(1),-Dvelocity(1))
!      call lsyssc_addConstant (rprjVector%RvectorBlock(2),-Dvelocity(2))
!
!    end if

    ! Write the velocity field

    ! CALL ucd_addVariableVertexBased (rexport,"X-vel",UCD_VAR_XVELOCITY, &
    !     p_Ddata(1:p_rtriangulation%NVT))
    ! CALL ucd_addVariableVertexBased (rexport,"Y-vel",UCD_VAR_YVELOCITY, &
    !     p_Ddata2(1:p_rtriangulation%NVT))
    call ucd_addVarVertBasedVec (rexport,"velocity",&
        p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))

    ! Write pressure
    call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
    call ucd_addVariableElementBased (rexport,"pressure",UCD_VAR_STANDARD, &
        p_Ddata(1:p_rtriangulation%NEL))

    ! If we have a simple Q1 discretisation in the pressure, write it out as it is
    if (rvector%p_rblockDiscr%RspatialDiscr(3)% &
        ccomplexity .eq. SPDISC_UNIFORM) then
      call spdiscr_getElemGroupInfo (rvector%p_rblockDiscr%RspatialDiscr(3),1,ieltype)

      if (elem_getPrimaryElement(ieltype) .eq. EL_Q1) then
        call lsyssc_getbase_double (rvector%RvectorBlock(3),p_Ddata)
        call ucd_addVariableVertexBased (rexport,"pressure",UCD_VAR_STANDARD, &
            p_Ddata(1:p_rtriangulation%NVT))
      else
        ! If this is QP1 or something else, project to Q1.
        call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        call spdp_projectToVertices (rvector%RvectorBlock(3),p_Ddata)
        call ucd_addVariableVertexBased (rexport,"pressure",UCD_VAR_STANDARD, &
            p_Ddata(1:p_rtriangulation%NVT))
      end if
    end if

    ! Write the solution in "raw" format. This does a kernal-internal projection
    ! into the Q1 space and writes out the solution in the vertices.
    ! Boundary conditions are not imposed.
    call ucd_addVectorFieldByVertex (rexport, "velocity_raw", UCD_VAR_STANDARD, &
        (/ rvector%RvectorBlock(1),rvector%RvectorBlock(2) /) )

    call ucd_addVectorByVertex (rexport, "pressure_raw", UCD_VAR_STANDARD, &
        rvector%RvectorBlock(3), DER_FUNC)

    ! If we have a simple Q1~ discretisation, calculate the streamfunction.
    if (rvector%p_rblockDiscr%RspatialDiscr(1)% &
        ccomplexity .eq. SPDISC_UNIFORM) then

      call spdiscr_getElemGroupInfo (rvector%p_rblockDiscr%RspatialDiscr(1),1,ieltype)

      if (elem_getPrimaryElement(ieltype) .eq. EL_Q1T) then

        call ppns2D_streamfct_uniform (rvector,rprjVector%RvectorBlock(1))

        call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        call ucd_addVariableVertexBased (rexport,"streamfunction",&
            UCD_VAR_STANDARD, p_Ddata(1:p_rtriangulation%NVT))

      end if

    end if

    ! Write out the viscosity if nonconstant
    if (rproblem%rphysics%cviscoModel .ne. 0) then

      ! Prepare the collection. The "next" collection points to the user defined
      ! collection.
      call ccmva_prepareViscoAssembly (rproblem,rproblem%rphysics,&
          rcollection,rvector,rproblem%rcollection)

      ! Project the viscosity to the Q1 space.
      call anprj_discrDirect(rprjVector%RvectorBlock(1), ffunctionViscoModel,&
          rcollection,ntempArrays=5)

      ! Write the viscosity
      call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
      call ucd_addVariableVertexBased (rexport,"viscosity",UCD_VAR_STANDARD, p_Ddata)
    end if

    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! Release the auxiliary vector
    call lsysbl_releaseVector (rprjVector)

    ! Release the discretisation structure.
    call spdiscr_releaseBlockDiscr (rprjDiscretisation)

    ! Release the dynamic level information structure including the
    ! discretised BC`s, not used anymore.
    call cc_doneDynamicLevelInfo (rdynamicInfo)

    ! Clean up the collection (as we are done with the assembly, that is it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

    if (present(dtime)) then
      ! Restore the current simulation time.
      rproblem%rtimedependence%dtime = dtimebackup
    end if

    if (present(dtime)) then
      ! Update time stamp of last written out GMV.
      rpostprocessing%inextFileSuffixUCD = rpostprocessing%inextFileSuffixUCD + 1
    end if

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_writeFilm (rpostprocessing,rvector,rproblem,dtime)

!<description>
  ! Writes Film output (raw data vectors) to a file as configured in the
  ! DAT file.
  !
  ! Note: This file is usually only used in a nonstationary simulation.
  ! In a stationary simulation, Film output makes no sense!
!</description>

!<input>
  ! Solution vector.
  type(t_vectorBlock), intent(in) :: rvector

  ! Simulation time.
  real(DP), intent(in) :: dtime
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(inout), target :: rproblem

  ! Postprocessing structure. Must have been initialised prior
  ! to calling this routine.
  ! The time stamp of the last written out Film file is updated.
  type(t_c2d2postprocessing), intent(inout) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ioutputFilm,ilevelFilm

    type(t_vectorBlock) :: rvector1,rvector2
    type(t_vectorScalar) :: rvectorTemp
    character(LEN=SYS_STRLEN) :: sfile,sfilename
    integer :: ilev
    integer :: NEQ
    type(t_interlevelProjectionBlock) :: rprojection
    logical :: bformatted

    ! Type of output:
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
                              "IOUTPUTFILM", ioutputFilm, 0)
    if (ioutputFilm .eq. 0) return

    ! Basic filename
    call parlst_getvalue_string (rproblem%rparamList, "CC-POSTPROCESSING", &
                                 "SFILENAMEFILM", sfilename, "", bdequote=.true.)

    ! Create the actual filename
    sfile = trim(adjustl(sfilename))//"."//sys_si0(rpostprocessing%inextFileSuffixFilm,5)

    ! Level of output:
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
                              "ILEVELFILM", ilevelFilm, 0)
    if (ilevelFilm .le. 0) then
      ilevelFilm = rproblem%NLMAX+ilevelFilm
    end if

    ilevelFilm = min(rproblem%NLMAX,max(rproblem%NLMIN,ilevelFilm))

    if (ilevelFilm .lt. rproblem%NLMIN) then
      call output_line ("Warning: Level for solution vector is < NLMIN! " // &
          "Writing out at level NLMIN!", &
          OU_CLASS_WARNING,OU_MODE_STD,"cc_releasePreconditioner")
      call sys_halt()
      ilevelFilm = rproblem%NLMIN
    end if

    ! Write formatted output?
    bformatted = ioutputFilm .ne. 2

    ! Interpolate the solution down to level istart.
    call lsysbl_copyVector (rvector,rvector1)   ! creates new rvector1!

    do ilev = rproblem%NLMAX,ilevelFilm+1,-1

      ! Initialise a vector for the lower level and a prolongation structure.
      call lsysbl_createVector (&
          rproblem%RlevelInfo(ilev-1)%rdiscretisation,rvector2,.false.)

      call mlprj_initProjectionVec (rprojection,rvector2)

      ! Interpolate to the next higher level.
      ! (Do not "restrict"! Restriction would be for the dual space = RHS vectors!)

      NEQ = mlprj_getTempMemoryVec (rprojection,rvector2,rvector1)
      if (NEQ .ne. 0) call lsyssc_createVector (rvectorTemp,NEQ,.false.)
      call mlprj_performInterpolation (rprojection,rvector2,rvector1, &
                                       rvectorTemp)
      if (NEQ .ne. 0) call lsyssc_releaseVector (rvectorTemp)

      ! Swap rvector1 and rvector2. Release the fine grid vector.
      call lsysbl_swapVectors (rvector1,rvector2)
      call lsysbl_releaseVector (rvector2)

      call mlprj_doneProjection (rprojection)

    end do

    call output_lbrk ()
    call output_line ("Writing Film file: "//sfile)

    ! Write out the solution.
    if (bformatted) then
      call vecio_writeBlockVectorHR (rvector1, "SOLUTION", .true.,&
         0, sfile, "(E22.15)")
    else
      call vecio_writeBlockVectorHR (rvector1, "SOLUTION", .true.,0, sfile)
    end if

    ! Release temp memory.
    call lsysbl_releaseVector (rvector1)

    ! Update time stamp of last written out Film file.
    rpostprocessing%inextFileSuffixFilm = rpostprocessing%inextFileSuffixFilm + 1

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine cc_initpostprocessing (rproblem,rpostprocessing)

!<description>
  ! Initialises the given postprocessing structure rpostprocessing
  ! according to the main problem rproblem. The structure can then be used
  ! to generate postprocessing data.
!</description>

!<input>
  ! A problem structure that describes the main problem to solve.
  type(t_problem), intent(in),target :: rproblem
!</input>

!<output>
  ! Postprocessing structure.
  type(t_c2d2postprocessing), intent(out) :: rpostprocessing
!</output>

!</subroutine>

  ! local variables
  type(t_blockDiscretisation), pointer :: p_rdiscr
  character(len=PARLST_LENLINEBUF) :: sexpression

    ! For postprocessing, we need discretisation structures in the Q0 and Q1 space,
    ! later perhaps in the Q2 space. For this purpose, derive the corresponding
    ! discretisation structure using the "main" discretisation structure on the
    ! maximum level.
    !
    ! For simplicity, we use only the discretisation structure of the X-velocity
    ! to derive everything.

    p_rdiscr => rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation

    ! Piecewise constant space:
    call spdiscr_deriveDiscr_triquad (p_rdiscr%RspatialDiscr(1), &
                 EL_P0, EL_Q0,rpostprocessing%rdiscrConstant)

    ! Piecewise linear space:
    call spdiscr_deriveDiscr_triquad (p_rdiscr%RspatialDiscr(1), &
                 EL_P1, EL_Q1, rpostprocessing%rdiscrLinear)

    ! Piecewise quadratic space:
    call spdiscr_deriveDiscr_triquad (p_rdiscr%RspatialDiscr(1), &
                 EL_P2, EL_Q2, rpostprocessing%rdiscrQuadratic)

    ! Initialise the time/file suffix when the first UCD file is to be written out.
    rpostprocessing%bnonstationaryPostprocessing = (rproblem%itimedependence .ne. 0)
    if (rproblem%itimedependence .ne. 0) then
      call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
         "ISTARTSUFFIXUCD", rpostprocessing%inextFileSuffixUCD, 1)
      call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
         "ISTARTSUFFIXFILM", rpostprocessing%inextFileSuffixFilm, 1)
    end if

    ! Prepare computation of L2/H1 errors
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "IERRORANALYSISL2", rpostprocessing%icalcL2, 0)
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "IERRORANALYSISH1", rpostprocessing%icalcH1, 0)
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "IERRORANALYSISTIMESPACE", rpostprocessing%icalcTimeSpaceDiscrErrors, 0)
    call parlst_getvalue_int (rproblem%rparamList, "CC-POSTPROCESSING", &
        "IERRORANALYSISTIMESPACEMETHOD", rpostprocessing%itimeSpaceDiscrErrorMethod, 0)


    ! Initialise a parser for the expressions.
    call fparser_create (rpostprocessing%rrefFunctionL2,NDIM2D+1)
    call fparser_create (rpostprocessing%rrefFunctionH1,NDIM2D*(NDIM2D+1))

    ! Parse expressions specifying the reference functions for the L2/H1 error.
    !
    ! L2 error
    if (rpostprocessing%icalcL2 .eq. 2) then
      ! Read parameters
      call parlst_getvalue_string (rproblem%rparamList,"CC-POSTPROCESSING",&
          "srefL2ExpressionU1",sexpression,"0",bdequote=.true.)
      call fparser_parseFunction (rpostprocessing%rrefFunctionL2, 1, sexpression, EXPRVARIABLES)

      call parlst_getvalue_string (rproblem%rparamList,"CC-POSTPROCESSING",&
          "srefL2ExpressionU2",sexpression,"0",bdequote=.true.)
      call fparser_parseFunction (rpostprocessing%rrefFunctionL2, 2, sexpression, EXPRVARIABLES)

      call parlst_getvalue_string (rproblem%rparamList,"CC-POSTPROCESSING",&
          "srefL2ExpressionP",sexpression,"0",bdequote=.true.)
      call fparser_parseFunction (rpostprocessing%rrefFunctionL2, 3, sexpression, EXPRVARIABLES)
    end if

    ! H1 error
    if (rpostprocessing%icalcH1 .eq. 2) then
      ! Read parameters
      call parlst_getvalue_string (rproblem%rparamList,"CC-POSTPROCESSING",&
          "srefH1ExpressionU1X",sexpression,"0",bdequote=.true.)
      call fparser_parseFunction (rpostprocessing%rrefFunctionH1, 1, sexpression, EXPRVARIABLES)

      call parlst_getvalue_string (rproblem%rparamList,"CC-POSTPROCESSING",&
          "srefH1ExpressionU1Y",sexpression,"0",bdequote=.true.)
      call fparser_parseFunction (rpostprocessing%rrefFunctionH1, 2, sexpression, EXPRVARIABLES)

      call parlst_getvalue_string (rproblem%rparamList,"CC-POSTPROCESSING",&
          "srefH1ExpressionU2X",sexpression,"0",bdequote=.true.)
      call fparser_parseFunction (rpostprocessing%rrefFunctionH1, 3, sexpression, EXPRVARIABLES)

      call parlst_getvalue_string (rproblem%rparamList,"CC-POSTPROCESSING",&
          "srefH1ExpressionU2Y",sexpression,"0",bdequote=.true.)
      call fparser_parseFunction (rpostprocessing%rrefFunctionH1, 4, sexpression, EXPRVARIABLES)

      call parlst_getvalue_string (rproblem%rparamList,"CC-POSTPROCESSING",&
          "srefH1ExpressionPX",sexpression,"0",bdequote=.true.)
      call fparser_parseFunction (rpostprocessing%rrefFunctionH1, 5, sexpression, EXPRVARIABLES)

      call parlst_getvalue_string (rproblem%rparamList,"CC-POSTPROCESSING",&
          "srefH1ExpressionPY",sexpression,"0",bdequote=.true.)
      call fparser_parseFunction (rpostprocessing%rrefFunctionH1, 6, sexpression, EXPRVARIABLES)
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine cc_copyPostprocessing (rpostprocessingSrc,rpostprocessingDst)

!<description>
  ! Copies the state of the postprocessing from rpostprocessingSrc
  ! to rpostprocessingDst.
  ! This is typically called to back up or to restore the current state
  ! of a running postprocessing structure.
!</description>

!<input>
  ! Source Postprocessing structure.
  type(t_c2d2postprocessing), intent(in) :: rpostprocessingSrc
!</input>

!<inputoutput>
  ! Destination Postprocessing structure.
  type(t_c2d2postprocessing), intent(inout) :: rpostprocessingDst
!</inputoutput>

!</subroutine>

    ! Initialise the time/file suffix when the first UCD file is to be written out.
    if (rpostprocessingSrc%bnonstationaryPostprocessing) then
      rpostprocessingDst%bnonstationaryPostprocessing = &
          rpostprocessingSrc%bnonstationaryPostprocessing

      rpostprocessingDst%inextFileSuffixUCD = rpostprocessingSrc%inextFileSuffixUCD
      rpostprocessingDst%inextFileSuffixFilm = rpostprocessingSrc%inextFileSuffixFilm
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine cc_clearpostprocessing (rpostprocessing)

!<description>
  ! Releases all calculated data in the given postprocessing structure
  ! so that it can be allocated again in the calculation routines.
  ! This routine must be called at the end of the postprocessing routines
  ! to release the temporary memory that was allocated for the vectors
  ! in the postprocessing structure.
!</description>

!<inputoutput>
  ! Postprocessing structure.
  type(t_c2d2postprocessing), intent(inout) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! Release all vectors which might be allocated.
    if (rpostprocessing%rvectorVelX%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorVelX)
    if (rpostprocessing%rvectorVelY%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorVelY)
    if (rpostprocessing%rvectorPressure%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorPressure)
    if (rpostprocessing%rvectorPressureCells%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorPressureCells)
    if (rpostprocessing%rvectorStreamfunction%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorStreamfunction)
    if (rpostprocessing%rvectorH1err%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorH1err)
    if (rpostprocessing%rvectorH1errCells%NEQ .ne. 0) &
      call lsyssc_releaseVector (rpostprocessing%rvectorH1errCells)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine cc_donepostprocessing (rpostprocessing)

!<description>
  ! Releases a given problem structure. All allocated memory of this structure
  ! is released.
!</description>

!<inputoutput>
  type(t_c2d2postprocessing), intent(inout) :: rpostprocessing
!</inputoutput>

!</subroutine>

    ! Release all vectors -- if there are still some allocated
    call cc_clearpostprocessing (rpostprocessing)

    ! Release the discretisation structures allocated above.
    call spdiscr_releaseDiscr(rpostprocessing%rdiscrQuadratic)
    call spdiscr_releaseDiscr(rpostprocessing%rdiscrLinear)
    call spdiscr_releaseDiscr(rpostprocessing%rdiscrConstant)

    call fparser_release (rpostprocessing%rrefFunctionH1)
    call fparser_release (rpostprocessing%rrefFunctionL2)

  end subroutine

end module

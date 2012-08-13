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
  
  use collection
  use convection
  
  use ucd
  
  use pprocnavierstokes
  use pprocerror
  
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
    call cc_errorAnalysis (rvector,0.0_DP,rproblem)
    
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
      dtimePrev,rvector,dtime,rvectorInt,dtimeInt,istep,rpostprocessing)
  
!<description>
  ! Postprocessing of solutions of stationary simulations.
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
  real(dp), intent(in) :: dtimePrev

  ! Solution vector of the current timestep.
  type(t_vectorBlock), intent(in) :: rvector

  ! Time of the current timestep.
  real(dp), intent(in) :: dtime
  
  ! Interpolated (in time) solution vector. Velocity and pressure
  ! represent the same point in time.
  type(t_vectorBlock), intent(in) :: rvectorInt
  
  ! Time of the interpolated solution vector.
  real(dp), intent(in) :: dtimeInt
  
  ! Number of the timestep. =0: initial solution
  integer, intent(in) :: istep
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
    real(dp) :: dweight,dwriteSolDeltaTime

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    ! Whether to apply postprocessing to rvector ot rvectorInt
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'ipostprocTimeInterpSolution', ipostprocTimeInterpSolution, 1)

    ! Think about writing out the solution...
    call parlst_getvalue_double (rproblem%rparamList, 'CC-DISCRETISATION', &
        'dwriteSolDeltaTime', dwriteSolDeltaTime, 0.0_DP)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
        'iwriteSolDeltaSteps', iwriteSolDeltaSteps, 1)
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
    
    if (ipostprocTimeInterpSolution .ne. 0) then
      ! Calculate body forces.
      call cc_calculateBodyForces (rvectorInt,dtimeInt,rproblem)
      
      ! Calculate point values
      call cc_evaluatePoints (rvectorInt,dtimeInt,rproblem)

      ! Calculate flux values
      call cc_evaluateFlux (rvector,dtimeInt,rproblem)

      ! Calculate the divergence
      call cc_calculateDivergence (rvectorInt,rproblem)

      ! Error analysis, comparison to reference function.
      call cc_errorAnalysis (rvectorInt,dtimeInt,rproblem)
    else
      ! Calculate body forces.
      call cc_calculateBodyForces (rvector,dtime,rproblem)
      
      ! Calculate point values
      call cc_evaluatePoints (rvector,dtime,rproblem)

      ! Calculate flux values
      call cc_evaluateFlux (rvector,dtime,rproblem)

      ! Calculate the divergence
      call cc_calculateDivergence (rvector,rproblem)

      ! Error analysis, comparison to reference function.
      call cc_errorAnalysis (rvector,dtime,rproblem)
    end if
    
    ! Write the UCD export file (GMV, AVS,...) as configured in the DAT file.
    !
    ! In a nonstationary simulation, first check if we are allowed
    ! to write something.

    call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'DMINTIMEUCD', dminTime, -1.E100_DP)
    call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'DMAXTIMEUCD', dmaxTime, 1.E100_DP)
    call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'DTIMEDIFFERENCEUCD', dtimeDifferenceUCD, 0.0_DP)
                                    
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
        call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
            'IINTERPOLATESOLUTIONUCD', iinterpolateSolutionUCD,1)
        if ((iinterpolateSolutionUCD .eq. 0) .or. (dtimeDifferenceUCD .eq. 0.0_DP) &
            .or. (dtimePrev .eq. dtime)) then
          ! No interpolation
          call cc_writeUCD (rpostprocessing, rvector, rproblem, dtime)
        else
          ! Interpolate and write out the interpolated solution.
          call lsysbl_copyVector (rvectorPrev,rintVector)
          dpptime = real(itime2,dp)*dtimeDifferenceUCD + rproblem%rtimedependence%dtimeInit
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
    call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'DMINTIMEFILM', dminTime, -1.E100_DP)
    call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'DMAXTIMEFILM', dmaxTime, 1.E100_DP)
    call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'DTIMEDIFFERENCEFILM', dtimeDifferenceFilm, 0.0_DP)
                                    
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
        call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
            'IINTERPOLATESOLUTIONFILM', iinterpolateSolutionFilm,1)
        if ((iinterpolateSolutionFilm .eq. 0) .or. (dtimeDifferenceFilm .eq. 0.0_DP) &
            .or. (dtimePrev .eq. dtime)) then
          ! No interpolation
          call cc_writeFilm (rpostprocessing, rvector, rproblem, dtime)
        else
          ! Interpolate and write out the interpolated solution.
          call lsysbl_copyVector (rvectorPrev,rintVector)
          dpptime = real(itime2,dp)*dtimeDifferenceFilm + rproblem%rtimedependence%dtimeInit
          dweight = (dpptime-dtimePrev) / (dtime-dtimePrev)
          call lsysbl_vectorLinearComb (rvector,rintVector,dweight,1.0_DP-dweight)
          call cc_writeFilm (rpostprocessing, rintVector, rproblem, dpptime)
          call lsysbl_releaseVector (rintVector)
        end if
      end if
    end if
    
    ! Gather statistics
    call stat_stopTimer(rtimer)
    rproblem%rstatistics%dtimePostprocessing = &
      rproblem%rstatistics%dtimePostprocessing + rtimer%delapsedReal

  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_errorAnalysis (rsolution,dtime,rproblem)

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
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>
    
    ! local variables
    real(DP),dimension(3) :: Derr
    real(DP) :: derrorVel, derrorP, denergy
    integer :: icalcL2,icalcH1,icalcEnergy
    integer :: iwriteErrorAnalysisL2,iwriteErrorAnalysisH1,iwriteKineticEnergy
    character(len=SYS_STRLEN) :: sfilenameErrorAnalysisL2
    character(len=SYS_STRLEN) :: sfilenameErrorAnalysisH1
    character(len=SYS_STRLEN) :: sfilenameKineticEnergy
    character(len=SYS_STRLEN) :: stemp
    integer :: iunit
    integer :: cflag
    integer :: icubError
    logical :: bfileExists
    real(DP) :: dtimebackup
    type(t_scalarCubatureInfo) :: rcubatureInfoUV,rcubatureInfoP
    
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IERRORANALYSISL2', icalcL2, 0)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IERRORANALYSISH1', icalcH1, 0)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'ICALCKINETICENERGY', icalcEnergy, 1)

    call parlst_getvalue_string (rproblem%rparamList,'CC-POSTPROCESSING',&
                                 'scubError',stemp,'')
    if (stemp .eq. '') then
      call parlst_getvalue_int (rproblem%rparamList,'CC-POSTPROCESSING',&
                                'icubError',icubError,int(CUB_GEN_AUTO))
    else
      icubError = cub_igetID(stemp)
    end if

    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IWRITEERRORANALYSISL2', iwriteErrorAnalysisL2, 0)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IWRITEERRORANALYSISH1', iwriteErrorAnalysisH1, 0)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IWRITEKINETICENERGY', iwriteKineticEnergy, 0)

    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'SFILENAMEERRORANALYSISL2', sfilenameErrorAnalysisL2, "", bdequote=.true.)
    if (sfilenameErrorAnalysisL2 .eq. '') &
      iwriteErrorAnalysisL2 = 0
    
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'SFILENAMEERRORANALYSISH1', sfilenameErrorAnalysisH1, "", bdequote=.true.)
    if (sfilenameErrorAnalysisH1 .eq. '') &
      iwriteErrorAnalysisH1 = 0
    
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'SFILENAMEKINETICENERGY', sfilenameKineticEnergy, "", bdequote=.true.)
    if (sfilenameKineticEnergy .eq. '') &
      iwriteKineticEnergy = 0
      
    ! Create an cubature info structure which contains our cubature rule
    call spdiscr_createDefCubStructure(&
        rsolution%RvectorBlock(1)%p_rspatialDiscr,rcubatureInfoUV,int(icubError,I32))
    call spdiscr_createDefCubStructure(&
        rsolution%RvectorBlock(3)%p_rspatialDiscr,rcubatureInfoP,int(icubError,I32))
    
    ! The error analysis might take place at an arbitrary time.
    ! Therefore, modify the 'current' time to the time where we want to
    ! do the analysis. Save the 'current' time and restore afterwards.
    ! This is necessary, as the user callback routines take the time
    ! for evaluation from here.
    dtimebackup = rproblem%rtimedependence%dtime
    rproblem%rtimedependence%dtime = dtime
    
    if ((icalcL2 .ne. 0) .or. (icalcH1 .ne. 0) .or. (icalcEnergy .ne. 0)) then
      call output_lbrk()
      call output_line ('Error Analysis')
      call output_line ('--------------')
    end if
    
    ! When writing to a file is enabled, delete the file in the first timestep.
    cflag = SYS_APPEND
    if (rproblem%rtimedependence%itimeStep .eq. 0) cflag = SYS_REPLACE
    
    if (icalcL2 .ne. 0) then
    
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    
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
      
      call output_line ('||u-reference||_L2 = '//trim(sys_sdEP(derrorVel,15,6)) )
      call output_line ('||p-reference||_L2 = '//trim(sys_sdEP(derrorP,15,6)) )
      
      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
      if (iwriteErrorAnalysisL2 .ne. 0) then
        ! Write the result to a text file.
        ! Format: timestep current-time value
        call io_openFileForWriting(sfilenameErrorAnalysisL2, iunit, &
            cflag, bfileExists,.true.)
        if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
          ! Write a headline
          write (iunit,'(A)') '# timestep time ||u-reference||_L2 ||p-reference||_L2'
        end if
        stemp = trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // ' ' &
            // trim(sys_sdEL(dtime,10)) // ' ' &
            // trim(sys_sdEL(derrorVel,10)) // ' ' &
            // trim(sys_sdEL(derrorP,10))
        write (iunit,'(A)') trim (stemp)
        close (iunit)
      end if
      
    end if

    if (icalcH1 .ne. 0) then
    
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    
      ! Perform error analysis to calculate and add ||u-z||_{H^1}.
      call pperr_scalar (PPERR_H1ERROR,Derr(1),rsolution%RvectorBlock(1),&
                         ffunction_TargetX,rproblem%rcollection,&
                         rcubatureInfo=rcubatureInfoUV)

      call pperr_scalar (PPERR_H1ERROR,Derr(2),rsolution%RvectorBlock(2),&
                         ffunction_TargetY,rproblem%rcollection,&
                         rcubatureInfo=rcubatureInfoUV)
                         
      if ((Derr(1) .ne. -1.0_DP) .and. (Derr(2) .ne. -1.0_DP)) then
        derrorVel = sqrt(Derr(1)**2+Derr(2)**2)

        call output_line ('||u-reference||_H1 = '//trim(sys_sdEP(derrorVel,15,6)),&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      else
        call output_line ('||u-reference||_H1 = not available',&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      end if
      
      call pperr_scalar (PPERR_H1ERROR,Derr(3),rsolution%RvectorBlock(3),&
                         ffunction_TargetP,rproblem%rcollection,&
                         rcubatureInfo=rcubatureInfoP)

      if (Derr(3) .ne. -1.0_DP) then
        derrorP = Derr(3)

        call output_line ('||p-reference||_H1 = '//trim(sys_sdEP(derrorP,15,6)),&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      else
        call output_line ('||p-reference||_H1 = not available',&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      end if
      
      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
      if (iwriteErrorAnalysisH1 .ne. 0) then
        ! Write the result to a text file.
        ! Format: timestep current-time value
        call io_openFileForWriting(sfilenameErrorAnalysisH1, iunit, &
            cflag, bfileExists,.true.)
        if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
          ! Write a headline
          write (iunit,'(A)') '# timestep time ||u-reference||_H1'
        end if
        stemp = trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // ' ' &
            // trim(sys_sdEL(dtime,10)) // ' ' &
            // trim(sys_sdEL(derrorVel,10))
        write (iunit,'(A)') trim(stemp)
        close (iunit)
      end if
      
    end if
    
    if (icalcEnergy .ne. 0) then
    
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    
      ! Perform error analysis to calculate and add 1/2||u||^2_{L^2}.
      call pperr_scalar (PPERR_L2ERROR,Derr(1),rsolution%RvectorBlock(1),&
          rcubatureInfo=rcubatureInfoUV)
      call pperr_scalar (PPERR_L2ERROR,Derr(2),rsolution%RvectorBlock(2),&
          rcubatureInfo=rcubatureInfoUV)
                         
      denergy = 0.5_DP*(Derr(1)**2+Derr(2)**2)

      call output_line ('||u_1||_L2         = '//trim(sys_sdEP(Derr(1),15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      call output_line ('||u_2||_L2         = '//trim(sys_sdEP(Derr(2),15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      call output_line ('||u||_L2           = '//&
          trim(sys_sdEP(sqrt(Derr(1)**2+Derr(2)**2),15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      call output_line ('1/2||u||^2_L2      = '//trim(sys_sdEP(denergy,15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      
      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
      if (iwriteKineticEnergy .ne. 0) then
        ! Write the result to a text file.
        ! Format: timestep current-time value
        call io_openFileForWriting(sfilenameKineticEnergy, iunit, &
            cflag, bfileExists,.true.)
        if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
          ! Write a headline
          write (iunit,'(A)') '# timestep time 1/2||u||^2_L2 ||u||_L2 ||u_1||_L2 ||u_2||_L2'
        end if
        stemp = trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // ' ' &
            // trim(sys_sdEL(dtime,10)) // ' ' &
            // trim(sys_sdEL(denergy,10)) // ' ' &
            // trim(sys_sdEL(sqrt(Derr(1)**2+Derr(2)**2),10)) // ' ' &
            // trim(sys_sdEL(Derr(1),10)) // ' ' &
            // trim(sys_sdEL(Derr(2),10))
        write (iunit,'(A)') trim(stemp)
        close (iunit)
      end if

    end if
    
    ! Release the cubature information structures
    call spdiscr_releaseCubStructure(rcubatureInfoUV)
    call spdiscr_releaseCubStructure(rcubatureInfoP)
    
    ! Restore the 'current' time.
    rproblem%rtimedependence%dtime = dtimebackup

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
    
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'icalcBodyForces', icalcBodyForces, 1)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'ibodyForcesFormulation', ibodyForcesFormulation, -1)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'ibodyForcesBdComponent', ibodyForcesBdComponent, 2)
    call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'dbdForcesCoeff1', dbdForcesCoeff1, rproblem%rphysics%dnu)
    call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'dbdForcesCoeff2', dbdForcesCoeff2, 0.1_DP * 0.2_DP**2)
    
    ! Probably cancel the calculation
    if (icalcBodyForces .eq. 0) return
    
    ! Information about writing body forces to a file.
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IWRITEBODYFORCES', iwriteBodyForces, 0)
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'SFILENAMEBODYFORCES', sfilenameBodyForces, "", bdequote=.true.)
    if (sfilenameBodyForces .eq. '') &
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
        call ppns2D_bdforces_uniform (rsolution,rregion,Dforces,CUB_G1_1D,&
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
        rcollection%p_rnextCollection => rproblem%rcollection
        call ccmva_prepareViscoAssembly (rproblem,rproblem%rphysics,&
            rcollection,rsolution)
          
        if (rproblem%rphysics%cviscoModel .eq. 0) then
          call ppns2D_bdforces_line (rsolution,rregion,Dforces,CUB_G1_1D,&
              dbdForcesCoeff1,dbdForcesCoeff2,cformulation)
        else
          call ppns2D_bdforces_line (rsolution,rregion,Dforces,CUB_G1_1D,&
              dbdForcesCoeff1,dbdForcesCoeff2,cformulation,ffunctionBDForcesVisco,&
              rcollection,ntempArrays=5)
        end if
        
      case (3)
        ! Old implementation:
        call ppns2D_bdforces_uniform (rsolution,rregion,Dforces,CUB_G4_2D,&
            dbdForcesCoeff1,dbdForcesCoeff2)

      case (4)
    
        ! Create a characteristic function of the boundary segment
        call lsyssc_createVecByDiscr (&
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
        rcollection%p_rnextCollection => rproblem%rcollection
        call ccmva_prepareViscoAssembly (rproblem,rproblem%rphysics,&
            rcollection,rsolution)
          
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
      call output_line ('Body forces')
      call output_line ('-----------')
      call output_line ('Body forces real bd., bdc/horiz/vert',&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG)
      call output_line (' '//trim(sys_siL(ibodyForcesBdComponent,10)) // ' / ' &
          //trim(sys_sdEP(Dforces(1),15,6)) // ' / '&
          //trim(sys_sdEP(Dforces(2),15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      
      if (iwriteBodyForces .ne. 0) then
        ! Write the result to a text file.
        ! Format: timestep current-time value
        call io_openFileForWriting(sfilenameBodyForces, iunit, &
            cflag, bfileExists,.true.)
        if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
          ! Write a headline
          write (iunit,'(A)') '# timestep time bdc horiz vert'
        end if
        stemp = trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // ' ' &
            // trim(sys_sdEL(dtime,10)) // ' ' &
            // trim(sys_siL(ibodyForcesBdComponent,10)) // ' ' &
            // trim(sys_sdEL(Dforces(1),10)) // ' '&
            // trim(sys_sdEL(Dforces(2),10))
        write (iunit,'(A)') trim(stemp)
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
    integer(I32) :: ieltype
    type(t_vectorScalar), target :: rtempVector
    
    if (rsolution%p_rblockDiscr%RspatialDiscr(1)% &
        ccomplexity .eq. SPDISC_UNIFORM) then
        
      ieltype = rsolution%p_rblockDiscr%RspatialDiscr(1)% &
                RelementDistr(1)%celement

      select case (elem_getPrimaryElement(ieltype))

      case (EL_Q1T, EL_P1T)
      
        ! Create a temporary vector
        call lsyssc_createVecByDiscr (rsolution%RvectorBlock(3)%p_rspatialDiscr,&
            rtempVector,.true.)

        ! Calculate divergence = D1 u1 + D2 u2
        call lsyssc_scalarMatVec (&
            rproblem%RlevelInfo(rproblem%nlmax)%rasmTempl%rmatrixD1, rsolution%RvectorBlock(1), &
            rtempVector, 1.0_DP, 0.0_DP)
        call lsyssc_scalarMatVec (&
            rproblem%RlevelInfo(rproblem%nlmax)%rasmTempl%rmatrixD2, rsolution%RvectorBlock(2), &
            rtempVector, 1.0_DP, 1.0_DP)
        
        call output_lbrk()
        call output_line ('Divergence')
        call output_line ('----------')
        call output_line ('Divergence = ' &
            //trim(sys_sdEP(lsyssc_vectorNorm(rtempVector,LINALG_NORML2),15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
            
        call lsyssc_releaseVector (rtempVector)
      
      end select
      
    end if
    
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
    npoints = parlst_querysubstrings (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'CEVALUATEPOINTVALUES')
        
    if (npoints .eq. 0) return
    
    ! Allocate memory for the values
    allocate(Dvalues(npoints))
    allocate(Dcoords(NDIM2D,npoints))
    allocate(Itypes(npoints))
    allocate(Ider(npoints))
    
    ! Read the points
    do i=1,npoints
      call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
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
    call output_line ('Point values')
    call output_line ('------------')
    do i=1,npoints
      write (sstr,"(A10,A,F9.4,A,F9.4,A,E16.10)") Sfctnames(1+Ider(i),Itypes(i)),&
          "(",Dcoords(1,i),",",Dcoords(2,i),") = ",Dvalues(i)
      call output_line(trim(sstr),coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
    end do
    
    ! Get information about writing the stuff into a DAT file.
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IWRITEPOINTVALUES', iwritePointValues, 0)
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'SFILENAMEPOINTVALUES', sfilenamePointValues, """""", bdequote=.true.)
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
        write (iunit,'(A)') &
          '# timestep time x y type deriv value x y type deriv value ...'
      end if
      stemp = &
          trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // ' ' // &
          trim(sys_sdEL(dtime,10))
      write (iunit,ADVANCE='NO',FMT='(A)') trim(stemp)
      do i=1,npoints
        stemp = ' ' //&
            trim(sys_sdEL(Dcoords(1,i),5)) // ' ' // &
            trim(sys_sdEL(Dcoords(2,i),5)) // ' ' // &
            trim(sys_siL(Itypes(i),2)) // ' ' // &
            trim(sys_siL(Ider(i),2)) // ' ' // &
            trim(sys_sdEL(Dvalues(i),10))
        write (iunit,ADVANCE='NO',FMT='(A)') trim(stemp)
      end do
      write (iunit,ADVANCE='YES',FMT='(A)') ""
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
    integer :: iderType, cflag
    logical :: bfileExists
    real(dp), dimension(:), allocatable :: Dvalues
    real(dp), dimension(:,:,:), allocatable :: Dcoords
    character(LEN=SYS_STRLEN) :: sparam
    character(LEN=SYS_STRLEN) :: sstr,sfilenameFluxValues,stemp

    ! Get the number of points to evaluate
    nlines = parlst_querysubstrings (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'CEVALUATEFLUXVALUES')
        
    if (nlines .eq. 0) return
    
    ! Allocate memory for the values
    allocate(Dvalues(nlines))
    allocate(Dcoords(NDIM2D,nlines,2))
    
    ! Read the points
    do i=1,nlines
      call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
          "CEVALUATEFLUXVALUES", sparam, "", i)
      read (sparam,*) Dcoords(1,i,1),Dcoords(2,i,1),Dcoords(1,i,2),Dcoords(2,i,2)
    end do
    
    ! Calculate the flux
    do i=1,nlines
      call ppns2D_calcFluxThroughLine (rsolution,Dcoords(1:2,i,1),Dcoords(1:2,i,2),Dvalues(i))
    end do
    
    ! Print the values to the terminal
    call output_lbrk()
    call output_line ('Flux values')
    call output_line ('-----------')
    do i=1,nlines
      write (sstr,"(A,F9.4,A,F9.4,A,F9.4,A,F9.4,A,E16.10)") "flux (",&
          Dcoords(1,i,1),",",Dcoords(2,i,1),")->(",&
          Dcoords(1,i,2),",",Dcoords(2,i,2),") = ",Dvalues(i)
      call output_line(trim(sstr),coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
    end do
    
    ! Get information about writing the stuff into a DAT file.
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IWRITEFLUXVALUES', iwriteFluxValues, 0)
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'SFILENAMEFLUXVALUES', sfilenameFluxValues, """""",bdequote=.true.)
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
        write (iunit,'(A)') &
          '# timestep time x1 y1 x2 y2 value1 x1 y1 x2 y2 value2 ...'
      end if
      stemp = &
          trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // ' ' // &
          trim(sys_sdEL(dtime,10))
      write (iunit,ADVANCE='NO',FMT='(A)') trim(stemp)
      do i=1,nlines
        stemp = ' ' //&
            trim(sys_sdEL(Dcoords(1,i,1),5)) // ' ' // &
            trim(sys_sdEL(Dcoords(2,i,1),5)) // ' ' // &
            trim(sys_sdEL(Dcoords(1,i,2),5)) // ' ' // &
            trim(sys_sdEL(Dcoords(2,i,2),5)) // ' ' // &
            trim(sys_sdEL(Dvalues(i),10))
        write (iunit,ADVANCE='NO',FMT='(A)') trim(stemp)
      end do
      write (iunit,ADVANCE='YES',FMT='(A)') ""
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

    ! A dynamic level information structure containing the BC's.
    type(t_dynamicLevelInfo), target :: rdynamicInfo
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    
    ! Backup of current simulation time.
    real(DP) :: dtimebackup
    
    integer :: ioutputUCD,ilevelUCD
    integer(I32) :: ieltype
    type(t_collection) :: rcollection
    
    ! Parameters used for the moving frame formulation
    integer :: imovingFrame
    real(DP), dimension(NDIM2D) :: Dvelocity,Dacceleration
    
    character(SYS_STRLEN) :: sfile,sfilename
    
    ! Type of output:
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                              'IOUTPUTUCD', ioutputUCD, 0)
    if (ioutputUCD .eq. 0) return

    ! Level of output:
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                              'ILEVELUCD', ilevelUCD, 0)
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
    ! For this purpose, first create a 'derived' simple discretisation
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
    call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)
    
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
                            
    ! Connect the vector to the BC`s
    rprjVector%p_rdiscreteBC => rdynamicInfo%rdiscreteBC
    
    ! The same way, discretise boundary conditions of fictitious boundary components.
    call cc_assembleFBDconditions (rproblem,rprjDiscretisation,&
        rdynamicInfo,rproblem%rcollection)
    rprjVector%p_rdiscreteBCfict => rdynamicInfo%rdiscreteFBC
    
    ! Filter the solution vector to implement discrete BC`s.
    call vecfil_discreteBCsol (rprjVector)

    ! Filter the solution vector to implement discrete BC`s for fictitious
    ! boundary components.
    call vecfil_discreteFBCsol (rprjVector)
    
    ! Basic filename
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                 'SFILENAMEUCD', sfilename, "", bdequote=.true.)
                                 
    ! Create the actual filename
    sfile = trim(adjustl(sfilename))//'.'//sys_si0(rpostprocessing%inextFileSuffixUCD,5)
                                 
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
    call output_line ('Writing visualisation file: '//sfile)
    
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
      call output_line ('Invalid visualisation output type.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'cc_writeUCD')
      call sys_halt()
    end select
        
    ! Is there a simulation time?
    if (present(dtime)) &
      call ucd_setSimulationTime (rexport,dtime)
    
    ! Write the configuration of the application as comment block
    ! to the output file.
    call ucd_addCommentLine (rexport,'Configuration:')
    call ucd_addCommentLine (rexport,'---------------')
    call ucd_addParameterList (rexport,rproblem%rparamList)
    call ucd_addCommentLine (rexport,'---------------')

    ! Get the velocity field
    call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
    call lsyssc_getbase_double (rprjVector%RvectorBlock(2),p_Ddata2)
    
    ! Moving frame velocity subtraction deactivated, gives pictures
    ! that can hardly be interpreted.
    
!    ! Is the moving-frame formulatino active?
!    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
!        'imovingFrame',imovingFrame,0)
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

    ! CALL ucd_addVariableVertexBased (rexport,'X-vel',UCD_VAR_XVELOCITY, &
    !     p_Ddata(1:p_rtriangulation%NVT))
    ! CALL ucd_addVariableVertexBased (rexport,'Y-vel',UCD_VAR_YVELOCITY, &
    !     p_Ddata2(1:p_rtriangulation%NVT))
    call ucd_addVarVertBasedVec (rexport,'velocity',&
        p_Ddata(1:p_rtriangulation%NVT),p_Ddata2(1:p_rtriangulation%NVT))
    
    ! Write pressure
    call lsyssc_getbase_double (rprjVector%RvectorBlock(3),p_Ddata)
    call ucd_addVariableElementBased (rexport,'pressure',UCD_VAR_STANDARD, &
        p_Ddata(1:p_rtriangulation%NEL))

    ! If we have a simple Q1 discretisation in the pressure, write it out as it is
    if (rvector%p_rblockDiscr%RspatialDiscr(3)% &
        ccomplexity .eq. SPDISC_UNIFORM) then
      ieltype = rvector%p_rblockDiscr%RspatialDiscr(3)% &
                RelementDistr(1)%celement
                
      if (elem_getPrimaryElement(ieltype) .eq. EL_Q1) then
        call lsyssc_getbase_double (rvector%RvectorBlock(3),p_Ddata)
        call ucd_addVariableVertexBased (rexport,'pressure',UCD_VAR_STANDARD, &
            p_Ddata(1:p_rtriangulation%NVT))
      else
        ! If this is QP1 or something else, project to Q1.
        call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        call spdp_projectToVertices (rvector%RvectorBlock(3),p_Ddata)
        call ucd_addVariableVertexBased (rexport,'pressure',UCD_VAR_STANDARD, &
            p_Ddata(1:p_rtriangulation%NVT))
      end if
    end if
    
    ! If we have a simple Q1~ discretisation, calculate the streamfunction.
    if (rvector%p_rblockDiscr%RspatialDiscr(1)% &
        ccomplexity .eq. SPDISC_UNIFORM) then
        
      ieltype = rvector%p_rblockDiscr%RspatialDiscr(1)% &
                RelementDistr(1)%celement
                
      if (elem_getPrimaryElement(ieltype) .eq. EL_Q1T) then
          
        call ppns2D_streamfct_uniform (rvector,rprjVector%RvectorBlock(1))
        
        call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
        call ucd_addVariableVertexBased (rexport,'streamfunction',&
            UCD_VAR_STANDARD, p_Ddata(1:p_rtriangulation%NVT))
            
      end if
      
    end if
    
    ! Write out the viscosity if nonconstant
    if (rproblem%rphysics%cviscoModel .ne. 0) then
      
      ! Prepare the collection. The "next" collection points to the user defined 
      ! collection.
      call ccmva_prepareViscoAssembly (rproblem,rproblem%rphysics,&
          rcollection,rvector)
      rcollection%p_rnextCollection => rproblem%rcollection

      ! Project the viscosity to the Q1 space.
      call anprj_discrDirect(rprjVector%RvectorBlock(1), ffunctionViscoModel,&
          rcollection,ntempArrays=5)

      ! Write the viscosity
      call lsyssc_getbase_double (rprjVector%RvectorBlock(1),p_Ddata)
      call ucd_addVariableVertexBased (rexport,'viscosity',UCD_VAR_STANDARD, p_Ddata)
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
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                              'IOUTPUTFILM', ioutputFilm, 0)
    if (ioutputFilm .eq. 0) return

    ! Basic filename
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
                                 'SFILENAMEFILM', sfilename, "", bdequote=.true.)
                                 
    ! Create the actual filename
    sfile = trim(adjustl(sfilename))//'.'//sys_si0(rpostprocessing%inextFileSuffixFilm,5)
                                 
    ! Level of output:
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
                              'ILEVELFILM', ilevelFilm, 0)
    if (ilevelFilm .le. 0) then
      ilevelFilm = rproblem%NLMAX+ilevelFilm
    end if
    
    ilevelFilm = min(rproblem%NLMAX,max(rproblem%NLMIN,ilevelFilm))
    
    if (ilevelFilm .lt. rproblem%NLMIN) then
      call output_line ('Warning: Level for solution vector is < NLMIN! ' // &
          'Writing out at level NLMIN!', &
          OU_CLASS_WARNING,OU_MODE_STD,'cc_releasePreconditioner')
      call sys_halt()
      ilevelFilm = rproblem%NLMIN
    end if
    
    ! Write formatted output?
    bformatted = ioutputFilm .ne. 2

    ! Interpolate the solution down to level istart.
    call lsysbl_copyVector (rvector,rvector1)   ! creates new rvector1!

    do ilev = rproblem%NLMAX,ilevelFilm+1,-1
      
      ! Initialise a vector for the lower level and a prolongation structure.
      call lsysbl_createVectorBlock (&
          rproblem%RlevelInfo(ilev-1)%rdiscretisation,rvector2,.false.)
      
      call mlprj_initProjectionVec (rprojection,rvector2)
      
      ! Interpolate to the next higher level.
      ! (Do not 'restrict'! Restriction would be for the dual space = RHS vectors!)

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
    call output_line ('Writing Film file: '//sfile)

    ! Write out the solution.
    if (bformatted) then
      call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,&
         0, sfile, '(E22.15)')
    else
      call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,0, sfile)
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

    ! For postprocessing, we need discretisation structures in the Q0 and Q1 space,
    ! later perhaps in the Q2 space. For this purpose, derive the corresponding
    ! discretisation structure using the 'main' discretisation structure on the
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
      call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
         'ISTARTSUFFIXUCD', rpostprocessing%inextFileSuffixUCD, 1)
      call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
         'ISTARTSUFFIXFILM', rpostprocessing%inextFileSuffixFilm, 1)
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
  
  end subroutine

end module

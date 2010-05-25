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
  use geometry
  use collection
  use convection
  use transformation
  use ucd
  use elementpreprocessing
  use pprocnavierstokes
  use pprocerror
  use dofmapping
  use ccboundaryconditionparser
  use ccgeneraldiscretisation
  use ccmatvecassembly
  use ccbasic
  use cccallback
  use geometryaux
  use analyticprojection
  use timestepping
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
    
    !rpostprocessing%rvectorScalarFB    
    type(t_vectorBlock) :: rvectorBlockProj
  
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

    if(rproblem%iParticles .gt. 0)then    
    ! Drag/Lift Calculation
      call cc_forcesNonStat(rvector,rproblem)
    end if

    ! Calculate body forces.
    call cc_calculateBodyForces (rvector,rproblem)
    
    ! Calculate point values
    call cc_evaluatePoints (rvector,rproblem)
    
    ! Calculate the divergence
    call cc_calculateDivergence (rvector,rproblem)
    
    ! Error analysis, comparison to reference function.
    call cc_errorAnalysis (rvector,rproblem)
    
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
      dtimePrev,rvector,dtime,istep,rpostprocessing,rtimestepping)
  
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
  
  ! Number of the timestep. =0: initial solution
  integer, intent(in) :: istep
  
  type (t_explicitTimeStepping), intent(in) ::rtimestepping  
  
!</input>

!</subroutine>

    ! local variables
    type(t_timer) :: rtimer
    real(DP) :: dminTime, dmaxTime, dtimeDifferenceUCD
    real(DP) :: dtimeDifferenceFilm, dpptime
    integer :: itime1,itime2,iinterpolateSolutionUCD,iinterpolateSolutionFilm
    integer :: iwriteSolDeltaSteps
    type(t_vectorBlock) :: rintVector
    real(dp) :: dweight,dwriteSolDeltaTime

    call stat_clearTimer(rtimer)
    call stat_startTimer(rtimer)

    ! Think about writing out the solution...
    call parlst_getvalue_double (rproblem%rparamList, 'CC-DISCRETISATION', &
        'dwriteSolDeltaTime', dwriteSolDeltaTime, 0.0_DP)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
        'iwriteSolDeltaSteps', iwriteSolDeltaSteps, 1)
    if (iwriteSolDeltaSteps .lt. 1) iwriteSolDeltaSteps = 1
    
    ! Figure out if we have to write the solution. This is the case if
    ! 1.) Precious and current time is the same or
    ! 2.) The solution crossed the next timestep.
  
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
    
    if(rproblem%iParticles .gt. 0)then    
    ! Drag/Lift Calculation
      call cc_forcesNonStat(rvector,rproblem)
    end if
    
    ! Calculate body forces.
    call cc_calculateBodyForces (rvector,rproblem)
    
    ! Calculate point values
    call cc_evaluatePoints (rvector,rproblem)

    ! Calculate the divergence
    call cc_calculateDivergence (rvector,rproblem)

    ! Error analysis, comparison to reference function.
    call cc_errorAnalysis (rvector,rproblem)
    
    ! move the particle
    call cc_updateParticlePosition(rproblem,rtimestepping%dtstep)
    !call cc_moveParticles(rproblem,rtimestepping%dtlaststep)
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
                                    
    if ((dtime .ge. dminTime-1000.0*SYS_EPSREAL) .and. &
        (dtime .le. dmaxTime+1000.0*SYS_EPSREAL)) then
    
      ! Figure out if we have to write the solution. This is the case if
      ! 1.) Precious and current time is the same ot
      ! 2.) The solution crossed the next ucd timestep.
    
      if ((dtimeDifferenceUCD .gt. 0.0_DP) .and. (dtimePrev .ne. dtime)) then
        itime1 = int((dtimePrev-rproblem%rtimedependence%dtimeInit)/dtimeDifferenceUCD)
        itime2 = int((dtime-rproblem%rtimedependence%dtimeInit)/dtimeDifferenceUCD)
      else
        itime1 = 0
        itime2 = 1
      end if
    
      if (itime1 .ne. itime2) then
        ! Proably we have to interpolate the solution to the point dtime in time.
        call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
            'IINTERPOLATESOLUTIONUCD', iinterpolateSolutionUCD,1)
        if ((iinterpolateSolutionUCD .eq. 0) .or. (dtimeDifferenceUCD .eq. 0.0_DP) &
            .or. (dtimePrev .eq. dtime)) then
          ! No interpolation
          call cc_writeUCD (rpostprocessing, rvector, rproblem, &
              rproblem%rtimedependence%dtime)
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
                                    
    if ((dtime .ge. dminTime-1000.0*SYS_EPSREAL) .and. &
        (dtime .le. dmaxTime+1000.0*SYS_EPSREAL)) then
    
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
        ! Proably we have to interpolate the solution to the point dtime in time.
        call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
            'IINTERPOLATESOLUTIONFILM', iinterpolateSolutionFilm,1)
        if ((iinterpolateSolutionFilm .eq. 0) .or. (dtimeDifferenceFilm .eq. 0.0_DP) &
            .or. (dtimePrev .eq. dtime)) then
          ! No interpolation
          call cc_writeFilm (rpostprocessing, rvector, rproblem, &
              rproblem%rtimedependence%dtime)
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

  subroutine cc_errorAnalysis (rsolution,rproblem)

!<description>
  ! Performs error analysis on a given solution rsolution as specified
  ! in the .DAT file.
  ! The result of the error analysis is written to the standard output.
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
    logical :: bfileExists
    
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IERRORANALYSISL2', icalcL2, 0)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IERRORANALYSISH1', icalcH1, 0)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'ICALCKINETICENERGY', icalcEnergy, 1)

    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IWRITEERRORANALYSISL2', iwriteErrorAnalysisL2, 0)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IWRITEERRORANALYSISH1', iwriteErrorAnalysisH1, 0)
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IWRITEKINETICENERGY', iwriteKineticEnergy, 0)

    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'SFILENAMEERRORANALYSISL2', stemp, '''''')
    read(stemp,*) sfilenameErrorAnalysisL2
    if (sfilenameErrorAnalysisL2 .eq. '') &
      iwriteErrorAnalysisL2 = 0
    
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'SFILENAMEERRORANALYSISH1', stemp, '''''')
    read(stemp,*) sfilenameErrorAnalysisH1
    if (sfilenameErrorAnalysisH1 .eq. '') &
      iwriteErrorAnalysisH1 = 0
    
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'SFILENAMEKINETICENERGY', stemp, '''''')
    read(stemp,*) sfilenameKineticEnergy
    if (sfilenameKineticEnergy .eq. '') &
      iwriteKineticEnergy = 0
    
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
      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1),&
                         ffunction_TargetX,rproblem%rcollection)

      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2),&
                         ffunction_TargetY,rproblem%rcollection)
                         
      derrorVel = sqrt(0.5_DP*(Derr(1)**2+Derr(2)**2))

      call pperr_scalar (rsolution%RvectorBlock(3),PPERR_L2ERROR,Derr(3),&
                         ffunction_TargetP,rproblem%rcollection)

      derrorP = sqrt(Derr(3))
      
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
        write (iunit,'(A)') trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // ' ' &
            // trim(sys_sdEL(rproblem%rtimedependence%dtime,10)) // ' ' &
            // trim(sys_sdEL(derrorVel,10)) // ' ' &
            // trim(sys_sdEL(derrorP,10))
        close (iunit)
      end if
      
    end if

    if (icalcH1 .ne. 0) then
    
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    
      ! Perform error analysis to calculate and add ||u-z||_{H^1}.
      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_H1ERROR,Derr(1),&
                         ffunction_TargetX,rproblem%rcollection)

      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_H1ERROR,Derr(2),&
                         ffunction_TargetY,rproblem%rcollection)
                         
      derrorVel = (0.5_DP*(Derr(1)**2+Derr(2)**2))

      call output_line ('||u-reference||_H1 = '//trim(sys_sdEP(derrorVel,15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      
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
        write (iunit,'(A)') trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // ' ' &
            // trim(sys_sdEL(rproblem%rtimedependence%dtime,10)) // ' ' &
            // trim(sys_sdEL(derrorVel,10))
        close (iunit)
      end if
      
    end if
    
    if (icalcEnergy .ne. 0) then
    
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    
      ! Perform error analysis to calculate and add 1/2||u||^2_{L^2}.
      call pperr_scalar (rsolution%RvectorBlock(1),PPERR_L2ERROR,Derr(1))
      call pperr_scalar (rsolution%RvectorBlock(2),PPERR_L2ERROR,Derr(2))
                         
      denergy = 0.5_DP*(Derr(1)**2+Derr(2)**2)

      call output_line ('||u||^2_L2         = '//trim(sys_sdEP(denergy,15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
      
      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
      if (iwriteKineticEnergy .ne. 0) then
        ! Write the result to a text file.
        ! Format: timestep current-time value
        call io_openFileForWriting(sfilenameKineticEnergy, iunit, &
            cflag, bfileExists,.true.)
        if ((cflag .eq. SYS_REPLACE) .or. (.not. bfileexists)) then
          ! Write a headline
          write (iunit,'(A)') '# timestep time ||u||^2_L2'
        end if
        write (iunit,'(A)') trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // ' ' &
            // trim(sys_sdEL(rproblem%rtimedependence%dtime,10)) // ' ' &
            // trim(sys_sdEL(denergy,10))
        close (iunit)
      end if

    end if

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

  subroutine cc_calculateBodyForces (rsolution,rproblem)

!<description>
  ! Calculates body forces as configured in the .DAT file.
  ! The result is written to the standard output.
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(in), target :: rsolution
!</input>

!<inputoutput>
  ! Problem structure.
  type(t_problem), intent(inout) :: rproblem
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
        'dbdForcesCoeff1', dbdForcesCoeff1, rproblem%dnu)
    call parlst_getvalue_double (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'dbdForcesCoeff2', dbdForcesCoeff2, 0.1_DP * 0.2_DP**2)
    
    ! Probably cancel the calculation
    if (icalcBodyForces .eq. 0) return
    
    ! Information about writing body forces to a file.
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IWRITEBODYFORCES', iwriteBodyForces, 0)
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'SFILENAMEBODYFORCES', stemp, '''''')
    read(stemp,*) sfilenameBodyForces
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
      rregion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      
      select case (icalcBodyForces)
      case (1)
        ! Old implementation:
        call ppns2D_bdforces_uniform (rsolution,rregion,Dforces,CUB_G1_1D,&
            dbdForcesCoeff1,dbdForcesCoeff2,cformulation)
        
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
          if (rproblem%isubequation .eq. 1) &
            cformulation = PPNAVST_DEFORMATIONTENSOR
        end select
          
        ! Prepare a collection structure in the form necessary for
        ! the computation of a nonconstant viscosity<
        !
        ! IquickAccess(1) = cviscoModel
        rcollection%IquickAccess(1) = rproblem%cviscoModel
        
        ! IquickAccess(2) = type of the tensor
        rcollection%IquickAccess(2) = rproblem%isubequation
        
        ! DquickAccess(1) = nu
        rcollection%DquickAccess(1) = rproblem%dnu
        
        ! DquickAccess(2) = dviscoexponent
        ! DquickAccess(3) = dviscoEps
        rcollection%DquickAccess(2) = rproblem%dviscoexponent
        rcollection%DquickAccess(3) = rproblem%dviscoEps
        
        ! The first quick access array specifies the evaluation point
        ! of the velocity -- if it exists.
        rcollection%p_rvectorQuickAccess1 => rsolution
          
        if (rproblem%cviscoModel .eq. 0) then
          call ppns2D_bdforces_line (rsolution,rregion,Dforces,CUB_G1_1D,&
              dbdForcesCoeff1,dbdForcesCoeff2,cformulation)
        else
          call ppns2D_bdforces_line (rsolution,rregion,Dforces,CUB_G1_1D,&
              dbdForcesCoeff1,dbdForcesCoeff2,cformulation,ffunctionBDForcesVisco,rcollection)
        end if
        
      case (3)
        ! Old implementation:
        call ppns2D_bdforces_uniform (rsolution,rregion,Dforces,CUB_G4_2D,&
            dbdForcesCoeff1,dbdForcesCoeff2,cformulation)

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
          if (rproblem%isubequation .eq. 1) &
            cformulation = PPNAVST_DEFORMATIONTENSOR
        end select
          
        ! Prepare a collection structure in the form necessary for
        ! the computation of a nonconstant viscosity<
        !
        ! IquickAccess(1) = cviscoModel
        rcollection%IquickAccess(1) = rproblem%cviscoModel
        
        ! IquickAccess(2) = type of the tensor
        rcollection%IquickAccess(2) = rproblem%isubequation
        
        ! DquickAccess(1) = nu
        rcollection%DquickAccess(1) = rproblem%dnu
        
        ! DquickAccess(2) = dviscoexponent
        ! DquickAccess(3) = dviscoEps
        rcollection%DquickAccess(2) = rproblem%dviscoexponent
        rcollection%DquickAccess(3) = rproblem%dviscoEps
        
        ! The first quick access array specifies the evaluation point
        ! of the velocity -- if it exists.
        rcollection%p_rvectorQuickAccess1 => rsolution
          
        if (rproblem%cviscoModel .eq. 0) then
          call ppns2D_bdforces_vol(rsolution,rcharfct,Dforces,&
              dbdForcesCoeff1,dbdForcesCoeff2,cformulation)
        else
          call ppns2D_bdforces_vol(rsolution,rcharfct,Dforces,&
              dbdForcesCoeff1,dbdForcesCoeff2,cformulation,ffunctionBDForcesVisco ,rcollection)
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
        write (iunit,'(A)') trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // ' ' &
            // trim(sys_sdEL(rproblem%rtimedependence%dtime,10)) // ' ' &
            // trim(sys_siL(ibodyForcesBdComponent,10)) // ' ' &
            // trim(sys_sdEL(Dforces(1),10)) // ' '&
            // trim(sys_sdEL(Dforces(2),10))
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
                
      if (elem_getPrimaryElement(ieltype) .eq. EL_Q1T) then
      
        ! Create a temporary vector 
        call lsyssc_createVecByDiscr (rsolution%RvectorBlock(3)%p_rspatialDiscr,&
            rtempVector,.true.)

        ! Calculate divergence = D1 u1 + D2 u2
        call lsyssc_scalarMatVec (&
            rproblem%RlevelInfo(rproblem%nlmax)%rstaticInfo%rmatrixD1, rsolution%RvectorBlock(1), &
            rtempVector, 1.0_DP, 0.0_DP)
        call lsyssc_scalarMatVec (&
            rproblem%RlevelInfo(rproblem%nlmax)%rstaticInfo%rmatrixD2, rsolution%RvectorBlock(2), &
            rtempVector, 1.0_DP, 1.0_DP)
        
        call output_lbrk()
        call output_line ('Divergence')
        call output_line ('----------')
        call output_line ('Divergence = ' &
            //trim(sys_sdEP(lsyssc_vectorNorm(rtempVector,LINALG_NORML2),15,6)),&
          coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
            
        call lsyssc_releaseVector (rtempVector)
      
      end if
      
    end if    
    
  end subroutine

!******************************************************************************

!<subroutine>

  subroutine cc_evaluatePoints (rsolution,rproblem)

!<description>
  ! Evaluates the solution in a number of points as configured in the DAT file.
!</description>
  
!<input>
  ! Solution vector to compute the norm/error from.
  type(t_vectorBlock), intent(in), target :: rsolution
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
    character(LEN=SYS_STRLEN) :: sstr,sfilenamePointValues
    character(LEN=10), dimension(3,3), parameter :: Sfctnames = RESHAPE (&
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
      call output_line(TRIM(sstr),coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
    end do
    
    ! Get information about writing the stuff into a DAT file.
    call parlst_getvalue_int (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'IWRITEPOINTVALUES', iwritePointValues, 0)
    call parlst_getvalue_string (rproblem%rparamList, 'CC-POSTPROCESSING', &
        'SFILENAMEPOINTVALUES', sstr, '''''')
    read(sstr,*) sfilenamePointValues
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
      write (iunit,ADVANCE='NO',FMT='(A)') &
          trim(sys_siL(rproblem%rtimedependence%itimeStep,10)) // ' ' // &
          trim(sys_sdEL(rproblem%rtimedependence%dtime,10))
      do i=1,npoints
        write (iunit,ADVANCE='NO',FMT='(A)') ' ' //&
            trim(sys_sdEL(Dcoords(1,i),5)) // ' ' // &
            trim(sys_sdEL(Dcoords(2,i),5)) // ' ' // &
            trim(sys_siL(Itypes(i),2)) // ' ' // &
            trim(sys_siL(Ider(i),2)) // ' ' // &
            trim(sys_sdEL(Dvalues(i),10))
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
    real(DP), dimension(:,:), pointer :: p_DvertexCoords,p_Ddata2D
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
    
    integer :: ioutputUCD,ilevelUCD, j,iin,iloop
    integer :: ieltype,ipart,ipolyHandle 
    
    ! Parameters used for the moving frame formulation
    integer :: imovingFrame
    real(DP), dimension(NDIM2D) :: Dvelocity,Dacceleration
    
    character(SYS_STRLEN) :: sfile,sfilename,sfilepoly
    
    type(t_geometryObject), pointer :: p_rgeometryObject
    
    type(t_geometryObject) :: rgeometryObject1
    
    real(DP), dimension(:,:), pointer :: p_Dvertices 
    
    integer, dimension(64) :: iconnect
    
    type(t_particleCollection), pointer :: p_rparticleCollection
    
    real(DP), dimension(2) :: Dp1,Dp2,DpointA
    real(dp), dimension(:,:), pointer :: p_DbdyEdg
    real(DP), dimension(:,:), pointer :: p_Dcoord
    real(DP), dimension(:), pointer   :: p_Ddata3   
    real(DP), dimension(:), pointer   :: p_Ddata4   
       
    
    call storage_getbase_double2D(rproblem%h_DedgesAtBoundary,p_DbdyEdg)    
    
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
    ! GMV understands only Q1 solutions! So the task is now to create
    ! a Q1 solution from rvector and write that out.
    !
    ! For this purpose, first create a 'derived' simple discretisation
    ! structure based on Q1 by copying the main guiding block discretisation
    ! structure and modifying the discretisation structures of the
    ! two velocity subvectors:
    
    call spdiscr_duplicateBlockDiscr(rvector%p_rblockDiscr,rprjDiscretisation)
    
    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(1), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(1))

    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(2), &
                 EL_Q1, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(2))

    call spdiscr_deriveSimpleDiscrSc (&
                 rvector%p_rblockDiscr%RspatialDiscr(3), &
                 EL_Q0, CUB_G2X2, &
                 rprjDiscretisation%RspatialDiscr(3))
                 
    ! The pressure discretisation substructure stays the old.
    !
    ! Now set up a new solution vector based on this discretisation,
    ! allocate memory.
    call lsysbl_createVecBlockByDiscr (rprjDiscretisation,rprjVector,.false.)
    
    ! Then take our original solution vector and convert it according to the
    ! new discretisation:
    call spdp_projectSolution (rvector,rprjVector)

    if (present(dtime)) then
      ! Only for the postprocessing, weitch to time dtime.
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
                                 'SFILENAMEUCD', sfile, '')
                                 
    ! Remove possible ''-characters
    read(sfile,*) sfilename
    
    ! Create the actual filename
    sfile = trim(adjustl(sfilename))//'.'//sys_si0(rpostprocessing%inextFileSuffixUCD,5)
    sfilepoly=trim(sfilename)//'.poly'//adjustl('')
    sfilepoly=trim(adjustl(sfilepoly))//'.'//sys_si0(rpostprocessing%inextFileSuffixUCD,5)
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
    call output_line ('Writing GMV file: '//sfile)
    
    select case (ioutputUCD)
    case (1)
      call ucd_startGMV (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)

    case (2)
      call ucd_startAVS (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile)
          
    case (3)
      call ucd_startVTK (rexport,UCD_FLAG_STANDARD,p_rtriangulation,sfile,sfilepoly)
          
    case DEFAULT
      call output_line ('Invalid UCD ooutput type.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'mysubroutine')
      stop
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
    
    ! Get the arrays with information of the source mesh.
    call storage_getbase_double2D (p_rtriangulation%h_DvertexCoords,&
        p_Dcoord)
    
!    allocate(p_Ddata3(p_rtriangulation%NVT))    
!    
!    do iloop=1,p_rtriangulation%NVT
!      Dp1(:)=p_DbdyEdg(:,7)
!      Dp2(:)=p_DbdyEdg(:,8)
!      DpointA(1)=p_Dcoord(1,iloop)
!      DpointA(2)=p_Dcoord(2,iloop)
!      p_Ddata3(iloop)=gaux_calcDistPEdg2D(DpointA,Dp1,Dp2)
!    end do    
!    
!    call ucd_addVariableVertexBased (rexport,'dist',UCD_VAR_STANDARD, &
!        p_Ddata3(1:p_rtriangulation%NVT))
!    
!    deallocate(p_Ddata3)
    
    ! Is the moving-frame formulatino active?
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
        'imovingFrame',imovingFrame,0)
        
    ! In this case, the postprocessing data must be modified by the
    ! moving frame velocity.
    if (imovingFrame .ne. 0) then
    
      ! Get the velocity and acceleration from the callback routine.
      call getMovingFrameVelocity (Dvelocity,Dacceleration,rproblem%rcollection)

      ! Subtract the moving frame velocity from the postprocessing
      ! data.
      call lsyssc_addConstant (rprjVector%RvectorBlock(1),-Dvelocity(1))
      call lsyssc_addConstant (rprjVector%RvectorBlock(2),-Dvelocity(2))
    
    end if

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
            p_Ddata(1:p_rtriangulation%NEL))
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
    
    
    allocate(p_Ddata3(p_rtriangulation%NVT))        
    allocate(p_Ddata4(p_rtriangulation%NVT))            
    if(rproblem%iParticles .gt. 0)then
      p_rparticleCollection => collct_getvalue_particles(rproblem%rcollection,'particles')
      do ipart=1,p_rparticleCollection%nparticles
        p_rgeometryObject => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject    
        
        call geom_polygonise(p_rgeometryObject,ipolyHandle)
        ! Get the vertices
        call storage_getbase_double2D(ipolyHandle, p_Dvertices)
        call ucd_addPolygon(rexport,p_Dvertices,ipart+1)
        call storage_free(ipolyHandle)
        ipolyHandle = ST_NOHANDLE
        ! postprocess forces        
        do iloop=1,p_rtriangulation%NVT
          call geom_isInGeometry (p_rgeometryObject, (/p_Dcoord(1,iloop),p_Dcoord(2,iloop)/), iin)
          if(iin .eq. 1)then
            p_Ddata3(iloop)=p_rparticleCollection%p_rParticles(ipart)%rForceWall(1)          
            p_Ddata4(iloop)=p_rparticleCollection%p_rParticles(ipart)%rForceWall(2)          
          else
            p_Ddata3(iloop)=0.0_dp          
            p_Ddata4(iloop)=0.0_dp          
          end if
        end do
        
      end do
    end if
    
    call ucd_addVarVertBasedVec (rexport,'forcewall',&
        p_Ddata3(1:p_rtriangulation%NVT),p_Ddata4(1:p_rtriangulation%NVT))
    
    deallocate(p_Ddata3)
    deallocate(p_Ddata4)
    
    ! Write the file to disc, that is it.
    call ucd_write (rexport)
    
    call ucd_release (rexport)    
    
!    deallocate(p_Data3)
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
                                 'SFILENAMEFILM', sfile, '')
                                 
    ! Remove possible ''-characters
    read(sfile,*) sfilename
    
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
    call spdiscr_deriveSimpleDiscrSc (&
                 p_rdiscr%RspatialDiscr(1), &
                 EL_Q0, CUB_G1X1, &
                 rpostprocessing%rdiscrConstant)

    ! Piecewise linear space:
    call spdiscr_deriveSimpleDiscrSc (&
                 p_rdiscr%RspatialDiscr(1), &
                 EL_Q1, CUB_G2X2, &
                 rpostprocessing%rdiscrLinear)
  
    ! Piecewise quadratic space:
    call spdiscr_deriveSimpleDiscrSc (&
                 p_rdiscr%RspatialDiscr(1), &
                 EL_Q2, CUB_G3X3, &
                 rpostprocessing%rdiscrQuadratic)
  
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
  
! ***************************************************************************

  !<subroutine>
  subroutine cc_forcesNonStat(rvector,rproblem)
  !<description>
  ! A routine that calculated the forces acting on an object
  ! and moves the object according to these forces
  !</description>

  ! structure for a geometry object
  !<inputoutput>
  type(t_problem), intent(INOUT) :: rproblem
  !</inputoutput>  

  !<input>
  type(t_vectorBlock), intent(IN) :: rvector
  !</input>
  
  !</subroutine>

  ! Local variables
  ! pointer to the entries of the alpha vector  
  real(DP), dimension(:), pointer :: p_Dvector  

  ! save the particle mass centers in a temp variable
  real(DP) :: dcenterx
  real(DP) :: dcentery

  ! pointer to the triangulation structure
  type(t_triangulation), pointer :: p_rtriangulation
  
  integer :: ipart
  
  type(t_particleCollection), pointer :: p_rparticleCollection

  type(t_geometryObject), pointer :: p_rgeometryObject        

  ! get the particle_collection out of the collection
  p_rparticleCollection => collct_getvalue_particles(rproblem%rcollection,'particles')

  ! loop over all particles to calculate the hydrodynamic forces for
  ! each particle  
  do ipart=1,p_rparticleCollection%nparticles
  
    ! if the vector contains data, we release it
    if (p_rparticleCollection%p_rParticles(ipart)%rvectorScalarFB%NEQ .ne. 0) &
      call lsyssc_releaseVector (p_rparticleCollection%p_rParticles(ipart)%rvectorScalarFB)

    ! create the fictitious boundary fem-vector from the discretisation 
    call lsyssc_createVecByDiscr(rvector%RvectorBlock(1)%p_rspatialDiscr, &
    p_rparticleCollection%p_rParticles(ipart)%rvectorScalarFB,.true.)

    ! get a pointer to the entries of the vector
    call lsyssc_getbase_double(p_rparticleCollection%p_rParticles(ipart)%rvectorScalarFB,p_Dvector)
    
    !print *,p_Dvector;
    
    
    ! get the pointer to the current geometry object
    p_rgeometryObject => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject    
    
    ! get a pointer to the triangulation
    p_rtriangulation => &
    p_rparticleCollection%p_rParticles(ipart)%rvectorScalarFB%p_rspatialDiscr%p_rtriangulation
    
    rproblem%rcollection%DQuickaccess(7)=ipart
    ! calculate the fem-function by a l2-projection
    ! here for every particle this has to be evaluated individually
    call anprj_discrDirect (p_rparticleCollection%p_rParticles(ipart)%rvectorScalarFB,cc_Particle,&
                            rproblem%rcollection,iorder=1)  
    
    
    
    !print *,p_Dvector;
    call output_lbrk ()
    call output_separator(OU_SEP_EQUAL)
    call output_line ('Q1 Vector recalculated ')
    call output_separator(OU_SEP_EQUAL)   
    
    ! get the center coordinates of the object
    ! to pass them to the forces function
    dcenterx=p_rgeometryObject%rcoord2D%Dorigin(1)
    dcentery=p_rgeometryObject%rcoord2D%Dorigin(2)
    
    ! we call this function to calculate the hydrodynamic forces
    call cc_forcesIntegrationNonStat(rproblem,rvector,&
         p_rparticleCollection%p_rParticles(ipart)%rvectorScalarFB,&
         dcenterx,dcentery,&
         p_rparticleCollection%p_rParticles(ipart)%rResForceX,&
         p_rparticleCollection%p_rParticles(ipart)%rResForceY,&
         p_rparticleCollection%p_rParticles(ipart)%dTorque,&
         rproblem%dnu,0.1_dp*0.2_dp**2)
    
    call output_lbrk()
    call output_line ('Drag forces')
    call output_line ('-----------')  
    print*, p_rparticleCollection%p_rParticles(ipart)%rResForceX(1)," / ",&
    p_rparticleCollection%p_rParticles(ipart)%rResForceY(1)
    
  end do  ! end particles
  
  end subroutine
  
! ***************************************************************************  

!<subroutine>
  subroutine cc_forcesIntegrationNonStat(rproblem,rvectorSol,&
             rvectorAlpha,dcenterx,dcentery,DforceX,DforceY,Dtor,df1,df2)
!<description>
  ! This routine calculates the error of a given finite element function
  ! in rvector to a given analytical callback function ffunctionReference.
  ! 2D version for double-precision vectors.
!</description>

  ! The body forces are defined as the integrals
  !
  !    Dforces(1) = 2/df2 * int_s [df1 dut/dn n_y - p n_x] ds 
  !    Dforces(2) = 2/df2 * int_s [df1 dut/dn n_x + p n_y] ds 

!<input>
  ! The FE solution vector. Represents a scalar FE function.
  type(t_vectorBlock), intent(in), target :: rvectorSol
  
  type(t_vectorScalar), intent(in), target :: rvectorAlpha
  
  !  parameter
  !integer,intent(IN) :: ipart
  
  ! center of the particle, we need is to
  ! calculate the torque
  real(DP), intent(IN):: dcenterx
  real(DP), intent(IN):: dcentery
  
  ! OPTIONAL: 1st weighting factor for the integral.
  ! If neglected, df1=1.0 is assumed.
  real(DP), intent(IN), optional      :: df1

  ! OPTIONAL: 2nd weighting factor for the integral.
  ! If neglected, df2=2.0 is assumed.
  real(DP), intent(IN), optional      :: df2
  
!</input>
  real(dp), dimension(2), intent(inout) :: DforceX
  real(dp), dimension(2), intent(inout) :: DforceY  
  real(dp), dimension(2), intent(inout) :: Dtor 
  type(t_problem), intent(INOUT) :: rproblem
!</subroutine>

    ! local variables
    integer :: i,k,icurrentElementDistr, ICUBP, NVE
    integer(I32) :: IEL, IELmax, IELset
    real(DP) :: OM, DN1, DN2, DN3,dpp,xtorque,ytorque,atq,atqy
    real(DP) :: ah1,ah2,du1x,du1y,du2x,du2y,Dfx,Dfy,dalx,daly,dTorque
    
    ! Cubature point coordinates on the reference element
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi
    
    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(CUB_MAXCUBP) :: Domega
    
    ! number of cubature points on the reference element
    integer :: ncubp
    
    ! Number of local degees of freedom for test functions
    integer :: indofTrial,indofFunc1,indofFunc2
    
    ! The triangulation structure - to shorten some things...
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! A pointer to an element-number list
    integer(I32), dimension(:), pointer :: p_IelementList
    
    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), allocatable :: p_DcubPtsRef
    
    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj
    
    ! Array for saving Jacobian 
    real(DP), dimension(7) :: Dj
    
    ! Current element distribution
    type(t_elementDistribution), pointer :: p_relementDistributionU
    
    type(t_elementDistribution), pointer :: p_relementDistributionP

    type(t_elementDistribution), pointer :: p_relementDistributionA
    
    ! Number of elements in the current element distribution
    integer(PREC_ELEMENTIDX) :: NEL

    ! Pointer to the values of the function that are computed by the callback routine.
    real(DP), dimension(:,:,:), allocatable :: Dcoefficients
    
    real(DP) :: negGrad
    
    ! Number of elements in a block. Normally =BILF_NELEMSIM,
    ! except if there are less elements in the discretisation.
    integer :: nelementsPerBlock
    
    ! A t_domainIntSubset structure that is used for storing information
    ! and passing it to callback routines.
    type(t_evalElementSet) :: rintSubset
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsTrial
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsFunc1
    ! An allocateable array accepting the DOF's of a set of elements.
    integer, dimension(:,:), allocatable, target :: IdofsFunc2
    
  
    ! Type of transformation from the reference to the real element 
    integer :: ctrafoType
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer :: cevaluationTag
    
    real(dp) :: dpf1,dpf2,robx,roby,length
    
    character(len=SYS_STRLEN) :: sfilenameBodyForces
    integer :: iunit
    integer :: cflag
    logical :: bfileExists    
    
    ! Prepare the weighting coefficients
    dpf1 = 1.0_DP
    dpf2 = 2.0_DP
    if (present(df1)) dpf1 = df1
    if (present(df2)) dpf2 = df2    

    ! make the l2 projection to get the normal vector

    ! Get a pointer to the triangulation - for easier access.
    p_rtriangulation => rvectorAlpha%p_rspatialDiscr%p_rtriangulation
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(PPERR_NELEMSIM,p_rtriangulation%NEL)
    
    ! Now loop over the different element distributions (=combinations
    ! of trial and test functions) in the discretisation.

    do icurrentElementDistr = 1,rvectorSol%p_rblockDiscr%RspatialDiscr(1)%inumFESpaces
    
      ! Activate the current element distribution
      p_relementDistributionU => &
      rvectorSol%p_rblockDiscr%RspatialDiscr(1)%RelementDistr(icurrentElementDistr)
      
      p_relementDistributionA =>&
      rvectorAlpha%p_rspatialDiscr%RelementDistr(icurrentElementDistr)

      p_relementDistributionP => &
      rvectorSol%p_rblockDiscr%RspatialDiscr(3)%RelementDistr(icurrentElementDistr)

    
      ! Cancel if this element distribution is empty.
      if (p_relementDistributionU%NEL .EQ. 0) cycle

      ! Get the number of local DOF's for trial functions
      indofTrial = elem_igetNDofLoc(p_relementDistributionU%celement)
      
      indofFunc1 = elem_igetNDofLoc(p_relementDistributionA%celement) 
      
      indofFunc2 = elem_igetNDofLoc(p_relementDistributionP%celement)
      
      ! Get the number of corner vertices of the element
      NVE = elem_igetNVE(p_relementDistributionU%celement)
      
      if (NVE .NE. elem_igetNVE(p_relementDistributionA%celement)) then
        print *,'cc_forcesIntegration: element spaces incompatible!'
        call sys_halt()
      end if      

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoType = elem_igetTrafoType(p_relementDistributionU%celement)

      ! Initialise the cubature formula,
      ! Get cubature weights and point coordinates on the reference element
      ! Now Dxi stores the point coordinates of the cubature points on the reference element
      call cub_getCubPoints(p_relementDistributionU%ccubTypeEval, ncubp, Dxi, Domega)
      
      ! Allocate some memory to hold the cubature points on the reference element
      allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))
      
      ! Reformat the cubature points; they are in the wrong shape!
      do i=1,ncubp
        do k=1,ubound(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i) = Dxi(i,k)
        end do
      end do
      
      ! Allocate memory for the DOF's of all the elements.
      allocate(IdofsTrial(indofTrial,nelementsPerBlock))
      
      allocate(IdofsFunc1(indofFunc1,nelementsPerBlock))
    
      allocate(IdofsFunc2(indofFunc2,nelementsPerBlock))      

      ! Allocate memory for the coefficients
      allocate(Dcoefficients(ncubp,nelementsPerBlock,13))
    
      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistributionU%celement)
      
      cevaluationTag = ior(cevaluationTag,elem_getEvaluationTag(p_relementDistributionP%celement))
      
      cevaluationTag = ior(cevaluationTag,elem_getEvaluationTag(p_relementDistributionA%celement))     
                      
      ! Make sure that we have determinants.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_DETJ)
      
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REALPOINTS)

      ! p_IelementList must point to our set of elements in the discretisation
      ! with that combination of trial functions
      call storage_getbase_int (p_relementDistributionU%h_IelementList, &
                                p_IelementList)
                     
      ! Get the number of elements there.
      NEL = p_relementDistributionU%NEL
    
      ! Initialize the forces of the element with zero
      Dfx     = 0.0_dp
      Dfy     = 0.0_dp
      dTorque = 0.0_dp
      
          
      ! Prepare the call to the evaluation routine of the analytic function.    
      CALL elprep_init(rintSubset)    
  
      ! Loop over the elements - blockwise.
      do IELset = 1, NEL, PPERR_NELEMSIM
      
        ! We always handle LINF_NELEMSIM elements simultaneously.
        ! How many elements have we actually here?
        ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
        ! elements simultaneously.
        
        IELmax = min(NEL,IELset-1+PPERR_NELEMSIM)
        
        ! Calculate the global DOF's into IdofsTrial.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF's of our LINF_NELEMSIM elements simultaneously.
        
        !--------------------------------------------------------------------------------
        call dof_locGlobMapping_mult(rvectorSol%p_rblockDiscr%RspatialDiscr(1), &
                                     p_IelementList(IELset:IELmax),IdofsTrial)
        !--------------------------------------------------------------------------------                                     
        !--------------------------------------------------------------------------------
        call dof_locGlobMapping_mult(rvectorAlpha%p_rspatialDiscr, &
                                     p_IelementList(IELset:IELmax),IdofsFunc1)
        !--------------------------------------------------------------------------------                                     
        !--------------------------------------------------------------------------------
        call dof_locGlobMapping_mult(rvectorSol%p_rblockDiscr%RspatialDiscr(3), &
                                     p_IelementList(IELset:IELmax),IdofsFunc2)
        !--------------------------------------------------------------------------------                                     
        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (rintSubset,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp))
        p_Ddetj => rintSubset%p_Ddetj

        ! In the next loop, we don't have to evaluate the coordinates
        ! on the reference elements anymore.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

        ! At this point, we must select the correct domain integration and coefficient
        ! calculation routine, depending which type of error we should compute!
        
        !----------------------------------------------------------------------------
        !                         EVALUATION PHASE
        !----------------------------------------------------------------------------
        !
        ! We need to build the following system:
        !    /                                                              \
        !   | |-p(x_i)             |   |du1/dx (x_i) du1/dy (x_i) ...|       |   | -dalpha/dx (x_i) |
        !   | |       -p(x_i)      | + |...  you know this works  ...| + u^t | * | -dalpha/dy (x_i) |
        !    \                                                              /
        !
        ! 
        ! Get the pressure in the cubature points
        ! Save the result to Dcoefficients(:,:,1)

        ! Build the p matrix
        call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(3), rintSubset, &
                p_relementDistributionP%celement, &
                IdofsFunc2, DER_FUNC, Dcoefficients(:,1:IELmax-IELset+1_I32,1))                

        ! Build the jacobi matrix of this (u1,u2,u3)
        ! First Row -------------------------------
        ! Save the result to Dcoefficients(:,:,2:4)
        call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(1), rintSubset, &
                p_relementDistributionU%celement, &
                IdofsTrial, DER_DERIV2D_X, Dcoefficients(:,1:IELmax-IELset+1_I32,2))  

        call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(1), rintSubset, &
                p_relementDistributionU%celement, &
                IdofsTrial, DER_DERIV2D_Y, Dcoefficients(:,1:IELmax-IELset+1_I32,3))  

        ! Second Row -------------------------------
        ! Save the result to Dcoefficients(:,:,4:5)
        call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(2), rintSubset, &
                p_relementDistributionU%celement, &
                IdofsTrial, DER_DERIV2D_X, Dcoefficients(:,1:IELmax-IELset+1_I32,4))  

        call fevl_evaluate_sim3 (rvectorSol%RvectorBlock(2), rintSubset, &
                p_relementDistributionU%celement, &
                IdofsTrial, DER_DERIV2D_Y, Dcoefficients(:,1:IELmax-IELset+1_I32,5))  

        ! Build the alpha vector
        ! Save the result to Dcoefficients(:,:,6:7)
        call fevl_evaluate_sim3 (rvectorAlpha, rintSubset,&
                p_relementDistributionA%celement, IdofsFunc1, DER_DERIV2D_X,&
                Dcoefficients(:,1:IELmax-IELset+1_I32,6))
                
        call fevl_evaluate_sim3 (rvectorAlpha, rintSubset,&
                p_relementDistributionA%celement, IdofsFunc1, DER_DERIV2D_Y,&
                Dcoefficients(:,1:IELmax-IELset+1_I32,7))

        ! Loop through elements in the set and for each element,
        ! loop through the DOF's and cubature points to calculate the
        ! integral: int_Omega (-p * I + Dj(u)) * (-grad(alpha)) dx
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp
            
            OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
            
            ! get the pressure
            dpp  = Dcoefficients(icubp,iel,1)
            
            ! x and y derivative of u1
            du1x = Dcoefficients(icubp,iel,2)
            du1y = Dcoefficients(icubp,iel,3)
            
            ! x and y derivative of u2
            du2x = Dcoefficients(icubp,iel,4)
            du2y = Dcoefficients(icubp,iel,5)
            
            dalx = Dcoefficients(icubp,iel,6)
            daly = Dcoefficients(icubp,iel,7)
            
            dn1  = -dalx
            dn2  = -daly
            
            ! gradient tensor
            ah1 = -dpp*dn1+dpf1*(du1x*dn1+du1y*dn2)
            ah2 = -dpp*dn2+dpf1*(du2x*dn1+du2y*dn2)

            ! deformation tensor
!            ah1 = -dpp*dn1+dpf1*(2.0_dp*du1x*dn1+(du1y+du2x)*dn2)
!            ah2 = -dpp*dn2+dpf1*((du1y+du2x)*dn1+2.0_dp*du2y*dn2)
            
            Dfx = Dfx + ah1 * om         
            Dfy = Dfy + ah2 * om
            
            ! for the torque calculate in 2d:
            ! (x-x_i) .perpdot. (sigma * n) 
            ! calculate the (x-x_i) part
            xtorque = rintSubset%p_DpointsReal(1,icubp,iel) - dcenterx
            ytorque = rintSubset%p_DpointsReal(2,icubp,iel) - dcentery
            
            robx = -1.0_dp*(rintSubset%p_DpointsReal(2,icubp,iel) - dcentery)
            roby = rintSubset%p_DpointsReal(1,icubp,iel) - dcenterx
            
            ! calculate the perpdot product
            atq = xtorque * ah2 - ytorque * ah1
            !atq = robx * ah1 + roby * ah2
            
            ! add up the forces
            dTorque = dTorque + atq * OM

          end do ! ICUBP 

        end do ! IEL
        
      end do ! IELset
      
      ! assign the forces
      DforceX(1) = Dfx
      DforceY(1) = Dfy
      Dtor(1)    = dTorque
      ! Output some values to get some readable stuff on the
      ! the screen
      print *,"--------------------------"
      print *,"torque force1: ",dTorque
      print *,"--------------------------"
      print *,"ForceY: ",Dfy
      print *,"--------------------------"
      print *,"ForceX: ",Dfx
            
      ! calculate also the drag and lift coefficients
      Dfx = Dfx * 2.0_dp/dpf2
      Dfy = Dfy * 2.0_dp/dpf2

      sfilenameBodyForces='ns/DLBFMDT001'
      cflag = SYS_APPEND      
      ! Write the result to a text file.
      ! Format: timestep current-time value
      call io_openFileForWriting(sfilenameBodyForces, iunit, &
          cflag, bfileExists,.true.)

      write (iunit,'(A)')trim(sys_sdEL(rproblem%rtimedependence%dtime,10)) // ' ' &
          // trim(sys_sdEL(Dfx,10)) // ' '&
          // trim(sys_sdEL(Dfy,10))
      close (iunit)
      
      ! save the coefficients
      rproblem%dCoefficientDrag = Dfx 
      rproblem%dCoefficientLift = Dfy
      
      ! Release memory
      call elprep_releaseElementSet(rintSubset)

      deallocate(p_DcubPtsRef)
      deallocate(Dcoefficients)
      deallocate(IdofsTrial)
      deallocate(IdofsFunc1)
      deallocate(IdofsFunc2)

    end do ! icurrentElementDistr
    
  end subroutine
  
! ***************************************************************************  
  
!<subroutine>
  subroutine cc_updateParticlePosition(rproblem,dtimestep)
!<description>
  ! 
!</description>

!<inputoutput>
  type(t_problem), intent(INOUT) :: rproblem
  real(dp), intent(in) :: dtimestep
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ipart,i,iedge1,iedge2,iseg,jpart,iwallforce,ipartforce
    real(DP) :: ah1,ah2,ah3,dvelx,dvely,dvelz,dmasssl,ddmasssl,dvolume,dimomir
    real(DP) :: dfx,dfy,length,dmasssl2,ddmasssl2,dvolume2
    real(DP) :: dCenterX,dCenterY,dCenterXold,dCenterYold,domega,CenterX2,CenterY2
    type(t_geometryObject), pointer :: p_rgeometryObject,p_rgeometryObject2
    
    ! this structure holds all the particles used in the application
    type(t_particleCollection), pointer :: p_rparticleCollection
    
    ! parameters for collision
    real(DP) :: eps1,eps2,pdist, ddist, deltaMesh, FWx,FWy,distout,dovlap
    real(dp) :: ew1,ew2,dwall,t,FPx,FPy,dprad,dprad2,dcx,dcy,dri,drj,dcxi
    real(dp) :: dcyi,dcyj,dcxj,dcuxi,dcuyi,dcuxj,dcuyj
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:,:), pointer :: p_Dcoords
    integer, dimension(:,:), pointer :: p_Iedges    
    real(DP), dimension(2) :: Dp1,Dp2,DpointA,DX,Dnormal,r,DV,DCollNormal
    real(dp), dimension(:,:), pointer :: p_DbdyEdg   
    real(dp),dimension(:),pointer :: Ddistances
    
    call cc_calcDistPart(rproblem,dtimestep)
    
    call storage_getbase_double2D(rproblem%h_DedgesAtBoundary,p_DbdyEdg)    
    
    p_rtriangulation => rproblem%RlevelInfo(rproblem%NLMAX)%rtriangulation
    
    call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords, p_Dcoords)
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtEdge, p_Iedges)
    
    iedge1=p_Iedges(1,1)
    iedge2=p_Iedges(2,1)
    
    Dp1(:)=p_Dcoords(:,iedge1)
    Dp2(:)=p_Dcoords(:,iedge2)

    deltaMesh = 0.015625_dp
    
    deltaMesh=sqrt((Dp1(1)-Dp2(1))**2 + (Dp1(2)-Dp2(2))**2)
    
    ! deltaMesh^2
    eps1 = deltaMesh**2
    
    ! deltaMesh
    eps2 = deltaMesh
    
    pdist = 0.5_dp * deltaMesh
    
    ew1=0.5_dp*(1.0_dp/eps1)
    ew2=0.5_dp*(1.0_dp/eps2)
    
    iwallforce=0
    ipartforce=1

    p_rparticleCollection => collct_getvalue_particles(rproblem%rcollection,'particles')

    p_rparticleCollection%dtime = 0.05_dp

    ! loop over all the particles and move eac of them
    ! individually
    do ipart=1,p_rparticleCollection%nparticles
  
      p_rgeometryObject => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject    
    
      ! position
      dCenterX=p_rgeometryObject%rcoord2D%Dorigin(1)
      dCenterY=p_rgeometryObject%rcoord2D%Dorigin(2)
      
      ! old position
      dCenterXold=p_rgeometryObject%rcoord2D%Dorigin(1)
      dCenterYold=p_rgeometryObject%rcoord2D%Dorigin(2)
      
      ! some physical parameters
      dvolume  = (p_rparticleCollection%p_rParticles(ipart)%drad)**2 * SYS_PI
      dmasssl  = p_rparticleCollection%p_rParticles(ipart)%drho * dvolume 
      ddmasssl = (p_rparticleCollection%p_rParticles(ipart)%drho-1.0_dp) * dvolume 
      
      ! Calculate the moment of inertia by formula,
      ! this formula is valid only for a circle
      ! note: for other shaped an approximation can
      ! be calculated by the routine cc_torque:
      ! there we calculated :int_Particle |r|^2 dParticle
      ! the moment of inertia can then be calculated by
      ! int_Particle [(COG_particle - X) .perpdot. u] dParticle / int_Particle |r|^2 dParticle
      ! The part : "int_Particle [(COG_particle - X) .perpdot. u] dParticle"
      ! we already calculated earlier; so we just need to divide...
      ! moment of inertia 
      dimomir  = &
      dmasssl*((p_rparticleCollection%p_rParticles(ipart)%drad)**2)/4.0_dp 
      
      ! Mean resistance forces for the given time step
      dfx = 0.5_dp * (p_rparticleCollection%p_rParticles(ipart)%rResForceX(2) &
       + p_rparticleCollection%p_rParticles(ipart)%rResForceX(1))
      dfy = 0.5_dp * (p_rparticleCollection%p_rParticles(ipart)%rResForceY(2) &
       + p_rparticleCollection%p_rParticles(ipart)%rResForceY(1))

      allocate(Ddistances(rproblem%isegments))
      FWy=0.0_dp
      FWx=0.0_dp
      !----------------------------------------------------------------------------------------   
      !                        CALCULATE WALL REPULSIVE FORCES
      !---------------------------------------------------------------------------------------- 
      if(iwallforce .eq. 1)then
        do iseg=1,rproblem%isegments
        
          Dp1(:)=p_DbdyEdg(:,2*(iseg-1)+1) 
          Dp2(:)=p_DbdyEdg(:,2*(iseg-1)+2)
          DpointA(1)=dCenterX
          DpointA(2)=dCenterY
          call gaux_calcDistPEdg2D(DpointA,Dp1,Dp2,ddist,t)
          
          DX(1)=Dp1(1)+t*(Dp2(1)-Dp1(1))
          DX(2)=Dp1(2)+t*(Dp2(2)-Dp1(2))
          Ddistances(iseg)=ddist
          call output_line ('SegmentNo: = '//trim(sys_si(iseg,2)) )
          dwall=1.5_dp  
          if(ddist .gt. dwall*p_rparticleCollection%p_rParticles(ipart)%drad + pdist)then
            FWx=FWx+0.0_dp
            FWy=FWy+0.0_dp
            call output_line ('Case 1: WallforcesX: = '//trim(sys_sdEP(FWx,15,6)) )
            call output_line ('Case 1: WallforcesY: = '//trim(sys_sdEP(FWy,15,6)) )
          elseif((ddist .gt. dwall*p_rparticleCollection%p_rParticles(ipart)%drad).and. &
                 (ddist .le. dwall*p_rparticleCollection%p_rParticles(ipart)%drad + pdist))then
            FWx=FWx+&
            ew1*(dCenterX-DX(1))*(dwall*p_rparticleCollection%p_rParticles(ipart)%drad + pdist - ddist)**2             
            FWy=FWy+&
            ew1*(dCenterY-DX(2))*(dwall*p_rparticleCollection%p_rParticles(ipart)%drad + pdist - ddist)**2
            call output_line ('Case 2: WallforcesX: = '//trim(sys_sdEP(FWx,15,6)) )
            call output_line ('Case 2: WallforcesY: = '//trim(sys_sdEP(FWy,15,6)) )
          elseif(ddist .le. dwall*p_rparticleCollection%p_rParticles(ipart)%drad)then
            FWx=FWx+&
            ew2*(dCenterX-DX(1))*(dwall*p_rparticleCollection%p_rParticles(ipart)%drad - ddist)      
            FWy=FWy+&
            ew2*(dCenterY-DX(2))*(dwall*p_rparticleCollection%p_rParticles(ipart)%drad - ddist)
            call output_line ('Case 3: WallforcesX: = '//trim(sys_sdEP(FWx,15,6)) )
            call output_line ('Case 3: WallforcesY: = '//trim(sys_sdEP(FWy,15,6)) )
          end if
        end do
      end if ! end if iwallforce
      p_rparticleCollection%p_rParticles(ipart)%rForceWall(1)=FWx
      p_rparticleCollection%p_rParticles(ipart)%rForceWall(2)=FWy

      !----------------------------------------------------------------------------------------   
      !                        CALCULATE PARTICLE REPULSIVE FORCES
      !---------------------------------------------------------------------------------------- 
      !initialize the forces for this particle by 0!
      FPx = 0.0_dp     
      FPy = 0.0_dp
      if(ipartforce .eq. 0)then     
      do jpart=1,p_rparticleCollection%nparticles
        if(jpart .eq. ipart)cycle
        p_rgeometryObject2=>p_rparticleCollection%p_rParticles(jpart)%rgeometryObject
        ddist=rproblem%dDistMatrix(ipart,jpart)
        dprad=p_rparticleCollection%p_rParticles(ipart)%drad
        dprad2=p_rparticleCollection%p_rParticles(jpart)%drad
        CenterX2=p_rgeometryObject2%rcoord2D%Dorigin(1)
        CenterY2=p_rgeometryObject2%rcoord2D%Dorigin(2)               
        if(ddist .gt. (dprad+dprad2+pdist))then
          FPx = FPx + 0.0_dp
          FPy = FPy + 0.0_dp
          call output_line ('Case 1: CollforcesX: = '//trim(sys_sdEP(0.0_dp,15,6)) )
          call output_line ('Case 1: CollforcesY: = '//trim(sys_sdEP(0.0_dp,15,6)) )          
        elseif((ddist .ge. dprad+dprad2).and.(ddist .le. dprad+dprad2+pdist))then
          FPx = FPx + (1.0_dp/eps1)*(dCenterX-CenterX2)*(dprad+dprad2+pdist-ddist)**2
          FPy = FPy + (1.0_dp/eps1)*(dCenterY-CenterY2)*(dprad+dprad2+pdist-ddist)**2
          call output_line ('Case 2: CollforcesX: = '//trim(sys_sdEP(&
          (1.0_dp/eps1)*(dCenterX-CenterX2)*(dprad+dprad2+pdist-ddist)**2,15,6)) )
          call output_line ('Case 2: CollforcesY: = '//trim(sys_sdEP(&
          (1.0_dp/eps1)*(dCenterY-CenterY2)*(dprad+dprad2+pdist-ddist)**2,15,6)) )          
        elseif(ddist .le. dprad+dprad2)then
          FPx = FPx + (1.0_dp/eps2)*(dCenterX-CenterX2)*(dprad+dprad2-ddist)
          FPy = FPy + (1.0_dp/eps2)*(dCenterY-CenterY2)*(dprad+dprad2-ddist)
          call output_line ('Case 3: CollforcesX: = '//trim(sys_sdEP(&
          (1.0_dp/eps2)*(dCenterX-CenterX2)*(dprad+dprad2-ddist),15,6)) )
          call output_line ('Case 3: CollforcesY: = '//trim(sys_sdEP(&
          (1.0_dp/eps2)*(dCenterY-CenterY2)*(dprad+dprad2-ddist),15,6)) )          
        end if
      
      end do
      end if ! ipartforce
      
      !----------------------------------------------------------------------------------------   
      !               CALCULATE CHANGE IN TRANSLATIONAL AND ANGULAR VELOCITY
      !----------------------------------------------------------------------------------------   
      ! Velocity difference for the given time step
      dvelx   = dtimestep*(FPx+FWx+1.0_dp*dfx+ddmasssl*0.0_dp)/dmasssl
      dvely   = dtimestep*(FPy+FWy+1.0_dp*dfy+ddmasssl*-9.81_dp)/dmasssl
      !dvely   = dtimestep*(FWy+1.0_dp*dfy+ddmasssl*0.0_dp)/dmasssl
   
      ! integrate the torque to get the angular velocity
      ! avel_new=time * 0.5*(torque_old + torque_new)/dimomir
      domega = dtimestep*0.5_dp*(p_rparticleCollection%p_rParticles(ipart)%dTorque(2) &
               + p_rparticleCollection%p_rParticles(ipart)%dTorque(1)) /dimomir

      !----------------------------------------------------------------------------------------   
      !                          UPDATE THE TORQUE
      !----------------------------------------------------------------------------------------   
      ! set the new values for the torque
      ! torque_old=torque_new
      p_rparticleCollection%p_rParticles(ipart)%dTorque(2)= &
      p_rparticleCollection%p_rParticles(ipart)%dTorque(1)
      
      !----------------------------------------------------------------------------------------   
      !                          UPDATE THE ANGLES AND ANGULAR VELOCITY
      !----------------------------------------------------------------------------------------   
      ! 
      ! a_new = a_old + time*(avel_old + 0.5*avel_new)
      ! Overwrite rotation angle
      p_rgeometryObject%rcoord2D%drotation = &
      p_rgeometryObject%rcoord2D%drotation + &
      dtimestep * (p_rparticleCollection%p_rParticles(ipart)%dangVelocity + 0.5_dp*domega)
      
      ! Recalculate SIN and COS values of angle
      p_rgeometryObject%rcoord2D%dsin_rotation = sin(p_rgeometryObject%rcoord2D%drotation)
      p_rgeometryObject%rcoord2D%dcos_rotation = cos(p_rgeometryObject%rcoord2D%drotation)
      
      ! Update the angular velocity
      ! avel_new = avel_old + avel_new
      p_rparticleCollection%p_rParticles(ipart)%dangVelocity= &
      p_rparticleCollection%p_rParticles(ipart)%dangVelocity + domega
    
      !----------------------------------------------------------------------------------------   
      !                          UPDATE THE FORCES
      !----------------------------------------------------------------------------------------   
      ! save the old forces for the next time step
      p_rparticleCollection%p_rParticles(ipart)%rResForceX(2) = &
      p_rparticleCollection%p_rParticles(ipart)%rResForceX(1)      

      p_rparticleCollection%p_rParticles(ipart)%rResForceY(2) = &
      p_rparticleCollection%p_rParticles(ipart)%rResForceY(1)      
      
      !----------------------------------------------------------------------------------------   
      !                          UPDATE THE PARTICLE POSITION
      !----------------------------------------------------------------------------------------   
      ! update the position of the x-center
      dCenterX = dCenterX + dtimestep * &
      (p_rparticleCollection%p_rParticles(ipart)%dtransVelX+0.5_dp*dvelx)

      ! update the position of the geometry object
      p_rgeometryObject%rcoord2D%Dorigin(1)=dCenterX
      
      ! if we are dealing with the reference particle
      ! update the position of the y-center
      dCenterY = dCenterY + dtimestep * &
      (p_rparticleCollection%p_rParticles(ipart)%dtransVelY + &
      0.5_dp*dvely)
      p_rgeometryObject%rcoord2D%Dorigin(2)=dCenterY


      !----------------------------------------------------------------------------------------   
      !                          UPDATE THE PARTICLE VELOCITY
      !----------------------------------------------------------------------------------------   
      ! save old velocity
      p_rparticleCollection%p_rParticles(ipart)%dtransVelXold=p_rparticleCollection%p_rParticles(ipart)%dtransVelX
      p_rparticleCollection%p_rParticles(ipart)%dtransVelYold=p_rparticleCollection%p_rParticles(ipart)%dtransVelY
      ! save the current velocity
      p_rparticleCollection%p_rParticles(ipart)%dtransVelX=&
      p_rparticleCollection%p_rParticles(ipart)%dtransVelX + dvelx
      p_rparticleCollection%p_rParticles(ipart)%dtransVelY=&
      p_rparticleCollection%p_rParticles(ipart)%dtransVelY + dvely
      

      ! normal collision model
      if(ipartforce .eq. 1)then
      do jpart=1,p_rparticleCollection%nparticles
        if(jpart .eq. ipart)cycle
        ! get the data of the 2nd particle      
        p_rgeometryObject2=>p_rparticleCollection%p_rParticles(jpart)%rgeometryObject
        ! get the distance between the centers
        ddist=rproblem%dDistMatrix(ipart,jpart)
        
        dprad=p_rparticleCollection%p_rParticles(ipart)%drad
        dprad2=p_rparticleCollection%p_rParticles(jpart)%drad
        
        if(ddist .le. (2.33_dp*dprad))then
          p_rparticleCollection%dtime=0.01_dp
        else
          p_rparticleCollection%dtime=0.1_dp
        end if
        
        ! get the center of the 2nd particle
        CenterX2=p_rgeometryObject2%rcoord2D%Dorigin(1)
        CenterY2=p_rgeometryObject2%rcoord2D%Dorigin(2)               
        
        ! the distance between the outlines
        distout=ddist-dprad+dprad2
        if(distout .lt. deltaMesh)then
        dovlap=deltaMesh-distout
        Dnormal(1)=CenterX2-dCenterX
        Dnormal(2)=CenterY2-dCenterY
        
        dvolume2  = &
        (p_rparticleCollection%p_rParticles(jpart)%drad)**2 * SYS_PI
        dmasssl2  = &
        p_rparticleCollection%p_rParticles(jpart)%drho * dvolume2 
        !AMASLDIJ=AMASSLDI+AMASSLDJ       
        dri=dmasssl*dovlap/(dmasssl2+dmasssl)
        drj=dmasssl2*dovlap/(dmasssl2+dmasssl)
        
        dcxi=dri*Dnormal(1)/dovlap
        dcyi=dri*Dnormal(2)/dovlap

        dcxj=drj*Dnormal(1)/dovlap
        dcyj=drj*Dnormal(2)/dovlap

        p_rgeometryObject%rcoord2D%Dorigin(1)=&
        p_rgeometryObject%rcoord2D%Dorigin(1)+dcxi
        
        p_rgeometryObject%rcoord2D%Dorigin(2)=&
        p_rgeometryObject%rcoord2D%Dorigin(2)+dcyi
        
        p_rgeometryObject2%rcoord2D%Dorigin(1)=&
        p_rgeometryObject2%rcoord2D%Dorigin(1)+dcxj
        
        p_rgeometryObject2%rcoord2D%Dorigin(2)=&
        p_rgeometryObject2%rcoord2D%Dorigin(2)+dcyj
        
        dcuxi=dcxi/dtimestep
        dcuyi=dcyi/dtimestep

        dcuxj=dcxj/dtimestep
        dcuyj=dcyj/dtimestep
  
        p_rparticleCollection%p_rParticles(ipart)%dtransVelX=&
        p_rparticleCollection%p_rParticles(ipart)%dtransVelX+dcuxi

        p_rparticleCollection%p_rParticles(ipart)%dtransVelY=&
        p_rparticleCollection%p_rParticles(ipart)%dtransVelY+dcuyi

        p_rparticleCollection%p_rParticles(jpart)%dtransVelX=&
        p_rparticleCollection%p_rParticles(jpart)%dtransVelX+dcuxj

        p_rparticleCollection%p_rParticles(jpart)%dtransVelY=&
        p_rparticleCollection%p_rParticles(jpart)%dtransVelY+dcuyj


        
        end if
        
      end do ! end jpart      
      end if ! ipartforce
      
      !----------------------------------------------------------------------------------------   
      !                          Collision with walls
      !----------------------------------------------------------------------------------------
      do iseg=1,rproblem%isegments
      
        Dp1(:)=p_DbdyEdg(:,2*(iseg-1)+1) 
        Dp2(:)=p_DbdyEdg(:,2*(iseg-1)+2)
        DpointA(1)=dCenterX
        DpointA(2)=dCenterY
        call gaux_calcDistPEdg2D(DpointA,Dp1,Dp2,ddist,t)
        
        DX(1)=Dp1(1)+t*(Dp2(1)-Dp1(1))
        DX(2)=Dp1(2)+t*(Dp2(2)-Dp1(2))
        Ddistances(iseg)=ddist
        
        if( Ddistances(iseg) .lt. p_rparticleCollection%p_rParticles(ipart)%drad+deltaMesh )then
          !p_rparticleCollection%dtime=0.01_dp
          DV(1)=p_rparticleCollection%p_rParticles(ipart)%dtransVelX
          DV(2)=p_rparticleCollection%p_rParticles(ipart)%dtransVelY
        
          Dp1(:)=p_DbdyEdg(:,2*(iseg-1)+1) 
          Dp2(:)=p_DbdyEdg(:,2*(iseg-1)+2)        
        
          Dnormal(1)=-1.0_dp*(Dp2(2)-Dp1(2))
          Dnormal(2)= (Dp2(1)-Dp1(1))    
  
          ! normalize the "normal"        
          length=sqrt(Dnormal(1)**2 + Dnormal(2)**2)
          Dnormal(1)=Dnormal(1)/length
          Dnormal(2)=Dnormal(2)/length
          
          ! reflect velocity vector 
          r(:)=DV(:)-2.0_dp*Dnormal(:)*(Dnormal(1)*DV(1)+Dnormal(2)*DV(2))
    
          ! compute the distance between the outline of
          ! the particle and the wall
          dcx=(p_rparticleCollection%p_rParticles(ipart)%drad-Ddistances(iseg)+deltaMesh) * Dnormal(1)
          dcy=(p_rparticleCollection%p_rParticles(ipart)%drad-Ddistances(iseg)+deltaMesh) * Dnormal(2)
          p_rparticleCollection%p_rParticles(ipart)%dtransVelX=&
          p_rparticleCollection%p_rParticles(ipart)%dtransVelX+dcx*dtimestep
          p_rparticleCollection%p_rParticles(ipart)%dtransVelY=&
          p_rparticleCollection%p_rParticles(ipart)%dtransVelY+dcy*dtimestep
          
          p_rgeometryObject%rcoord2D%Dorigin(:)=&
          p_rgeometryObject%rcoord2D%Dorigin(:) + &
          (p_rparticleCollection%p_rParticles(ipart)%drad-Ddistances(iseg)+deltaMesh) * Dnormal(:)

!          ! the position is updated, now update the velocity
!          p_rparticleCollection%p_rParticles(ipart)%dtransVelX=0.984*r(1)
!          p_rparticleCollection%p_rParticles(ipart)%dtransVelY=0.984*r(2)
!          call output_line ('VelX: = '//trim(sys_sdEP(DV(1),15,6)) )                    
!          call output_line ('Vely: = '//trim(sys_sdEP(DV(2),15,6)) )                    
!          call output_line ('RefVelX: = '//trim(sys_sdEP(r(1),15,6)) )                    
!          call output_line ('RefVelX: = '//trim(sys_sdEP(r(2),15,6)) )                    
!          call output_line ('NormalX: = '//trim(sys_sdEP(Dnormal(1),15,6)) )                    
!          call output_line ('NormalY: = '//trim(sys_sdEP(Dnormal(2),15,6)) ) 
!          call output_line ('dcx: = '//trim(sys_sdEP(dcx,15,6)) )                    
!          call output_line ('dcy: = '//trim(sys_sdEP(dcy,15,6)) ) 
!          call output_line ('NewPoX: = '//trim(sys_sdEP( p_rgeometryObject%rcoord2D%Dorigin(1),15,6)) )                    
!          call output_line ('NewPoY: = '//trim(sys_sdEP( p_rgeometryObject%rcoord2D%Dorigin(2),15,6)) ) 
          
        end if
      end do
      
      !----------------------------------------------------------------------------------------   
      !                          SCREEN OUTPUT
      !----------------------------------------------------------------------------------------   
      ! Output some values to get some readable stuff on the
      ! the screen
      call output_line('--------------------------')
      call output_line ('torque force new = '//trim(sys_sdEP(p_rparticleCollection%p_rParticles(ipart)%dTorque(1),15,6)) )
      call output_line('--------------------------')
      call output_line ('torque force old = '//trim(sys_sdEP(p_rparticleCollection%p_rParticles(ipart)%dTorque(2),15,6)) )
      call output_line('--------------------------')
      call output_line ('drag force new = '//trim(sys_sdEP(p_rparticleCollection%p_rParticles(ipart)%rResForceX(1),15,6)) )
      call output_line('--------------------------')
      call output_line ('drag force old = '//trim(sys_sdEP(p_rparticleCollection%p_rParticles(ipart)%rResForceX(2),15,6)) )
      call output_line('--------------------------')
      call output_line ('lift force new = '//trim(sys_sdEP(p_rparticleCollection%p_rParticles(ipart)%rResForceY(1),15,6)) )
      call output_line('--------------------------')
      call output_line ('lift force old = '//trim(sys_sdEP(p_rparticleCollection%p_rParticles(ipart)%rResForceY(2),15,6)) )
      call output_line('--------------------------')
      call output_line ('DAngVel = '//trim(sys_sdEP(p_rparticleCollection%p_rParticles(ipart)%dangVelocity,15,6)) )
      call output_line('--------------------------')
      call output_line ('dCenterX = '//trim(sys_sdEP(dCenterX,15,6)) )
      call output_line('--------------------------')
      call output_line ('dCenterY = '//trim(sys_sdEP(dCenterY,15,6)) )
      call output_line('--------------------------')
      call output_line ('VelocityX new = '//trim(sys_sdEP(p_rparticleCollection%p_rParticles(ipart)%dtransVelX,15,6)) )
      call output_line('--------------------------')
      call output_line ('VelocityY new = '//trim(sys_sdEP(p_rparticleCollection%p_rParticles(ipart)%dtransVelY,15,6)) )
      call output_line('--------------------------')
      call output_line ('VelocityX old = '//trim(sys_sdEP(p_rparticleCollection%p_rParticles(ipart)%dtransVelXold,15,6)) )
      call output_line('--------------------------')
      call output_line ('VelocityY old = '//trim(sys_sdEP(p_rparticleCollection%p_rParticles(ipart)%dtransVelYold,15,6)) )
      call output_line('--------------------------')
      
      deallocate(Ddistances)
      
    end do ! end ipart
    
  end subroutine

! ***************************************************************************  
  
!<subroutine>
  subroutine cc_Particle(cderivative,rdiscretisation, &
                nelements,npointsPerElement,Dpoints, &
                IdofsTest,rdomainIntSubset,&
                Dvalues,rcollection)
  
  use basicgeometry
  use triangulation
  use collection
  use scalarpde
  use domainintegration
  
!<description>
  ! This subroutine is called during the calculation of errors. It has to compute
  ! the (analytical) values of a function in a couple of points on a couple
  ! of elements. These values are compared to those of a computed FE function
  ! and used to calculate an error.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points) in in real coordinates.
  ! According to the terms in the linear form, the routine has to compute
  ! simultaneously for all these points.
!</description>
  
!<input>
  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
  ! The result must be written to the Dvalue-array below.
  integer, intent(IN)                                         :: cderivative

  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
  
  ! Number of elements, where the coefficients must be computed.
  integer, intent(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  integer, intent(IN)                                         :: npointsPerElement
  
  ! This is an array of all points on all the elements where coefficients
  ! are needed.
  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
  real(DP), dimension(:,:,:), intent(IN)                      :: Dpoints

  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(\#local DOF's in trial space,Number of elements)
  integer, dimension(:,:), intent(IN) :: IdofsTest

  ! This is a t_domainIntSubset structure specifying more detailed information
  ! about the element set that is currently being integrated.
  ! It's usually used in more complex situations (e.g. nonlinear matrices).
  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

  ! Optional: A collection structure to provide additional 
  ! information to the coefficient routine. 
  type(t_collection), intent(INOUT), optional      :: rcollection
  
!</input>

!<output>
  ! This array has to receive the values of the (analytical) function
  ! in all the points specified in Dpoints, or the appropriate derivative
  ! of the function, respectively, according to cderivative.
  !   DIMENSION(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
!</output>
  
!</subroutine>
  ! local variables
  integer :: i,j,iin,ipart
  type(t_geometryObject), pointer :: p_rgeometryObject
  
  type(t_particleCollection), pointer :: p_rparticleCollection

  p_rparticleCollection => collct_getvalue_particles(rcollection,'particles')
  !p_rgeometryObject => collct_getvalue_geom(rcollection,'mini')
  
  ipart=rcollection%DQuickaccess(7)
  p_rgeometryObject => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject    
  select case (cderivative)
  case (DER_FUNC)
  
  ! loop over all elements and calculate the
  ! values in the cubature points
  do i=1,nelements
    do j=1,npointsPerElement 
      
      ! Get the distance to the center
      call geom_isInGeometry (p_rgeometryObject, (/Dpoints(1,j,i),Dpoints(2,j,i)/), iin)
      ! check if it is inside      
      if(iin .eq. 1)then 
        Dvalues(j,i) =  1.0_DP 
      else
        Dvalues(j,i) = 0.0_DP
      end if
      
    end do
  end do    
    
  case (DER_DERIV_X)
    ! Not really useful in the case at hand
    Dvalues (:,:) = 0.0_dp
  case (DER_DERIV_Y)
    ! not much better than the above case
    Dvalues (:,:) = 0.0_dp
  case DEFAULT
    ! Unknown. Set the result to 0.0.
    Dvalues = 0.0_DP
  end select
  
  
  end subroutine

! *************************************************************************** 

!<subroutine>
  subroutine cc_calcDistPart(rproblem,dtimestep)
!<description>
  ! 
!</description>

!<inputoutput>
  type(t_problem), intent(INOUT) :: rproblem
  real(dp), intent(in) :: dtimestep
!</inputoutput>

!</subroutine>
    real(DP) :: dCenterX1,dCenterY1,dCenterX2,dCenterY2
    integer  :: ipart,jpart
    type(t_geometryObject), pointer :: p_rgeometryObject1,p_rgeometryObject2
    
    ! this structure holds all the particles used in the application
    type(t_particleCollection), pointer :: p_rparticleCollection

    p_rparticleCollection => collct_getvalue_particles(rproblem%rcollection,'particles')

    ! loop over all the particles and move eac of them
    ! individually
    do ipart=1,p_rparticleCollection%nparticles
    
      p_rgeometryObject1 => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject    
      
      ! position
      dCenterX1=p_rgeometryObject1%rcoord2D%Dorigin(1)
      dCenterY1=p_rgeometryObject1%rcoord2D%Dorigin(2)
      
      rproblem%dDistMatrix(ipart,ipart)=0.0_dp
      
      do jpart=ipart+1,p_rparticleCollection%nparticles
        ! position
        p_rgeometryObject2 => p_rparticleCollection%p_rParticles(jpart)%rgeometryObject
        dCenterX2=p_rgeometryObject2%rcoord2D%Dorigin(1)
        dCenterY2=p_rgeometryObject2%rcoord2D%Dorigin(2)
        
        rproblem%dDistMatrix(ipart,jpart)=sqrt((dCenterX2-dCenterX1)**2+(dCenterY2-dCenterY1)**2)
        rproblem%dDistMatrix(jpart,ipart)=rproblem%dDistMatrix(ipart,jpart)
      
      end do

    end do

  end subroutine
  
! *************************************************************************** 
  
! <subroutine>
  subroutine cc_moveParticles(rproblem,dtimestep)
!<description>
  ! 
!</description>

!<inputoutput>
  type(t_problem), intent(INOUT) :: rproblem
  real(dp), intent(in) :: dtimestep
!</inputoutput>

!</subroutine>

    ! local variables
    integer  :: ipart,isegments,n,model
    real(DP) :: dnu,dvelx,dvely,dmasssl,ddmasssl,dvolume,dimomir
    real(DP) :: dfx,dfy,dfx1,dfy1,dtrq,dtrq1
    real(DP) :: dCenterX,dCenterY,DistX,DistY,domega
    real(dp), dimension(:,:), pointer :: p_DbdyEdg
    type(t_geometryObject), pointer :: p_rgeometryObject
    
    ! this structure holds all the particles used in the application
    type(t_particleCollection), pointer :: p_rparticleCollection

    call storage_getbase_double2D(rproblem%h_DedgesAtBoundary,p_DbdyEdg)

    p_rparticleCollection => collct_getvalue_particles(rproblem%rcollection,'particles')

    n = p_rparticleCollection%nparticles
    isegments = rproblem%isegments
    dnu = rproblem%dnu
    
    !----------------------------------------------------------------------------------------   
    !                           COLLISION OF PARTICLES WITH WALL
    !----------------------------------------------------------------------------------------   
!     call collisionwall(p_rparticleCollection,dtimestep)
    call collisionwall_adv(p_rparticleCollection,isegments,p_DbdyEdg,dtimestep)
    
    !----------------------------------------------------------------------------------------   
    !                         COLLISION OF PARTICLES WITH OTHER PARTICLES
    !----------------------------------------------------------------------------------------   
    ! model name (1,2,3,4 or 5)
    ! since we have five models at the time being
    model = 3
    if (model .eq. 1) then
    call collisionparticle1(p_rparticleCollection,dnu,dtimestep)
    ! end if
    end if
    if (model .eq. 2) then
    call collisionparticle2(p_rparticleCollection,dnu,dtimestep)
    ! end if
    end if
    if (model .eq. 3) then
    call collisionparticle3(p_rparticleCollection,dnu,dtimestep)
    ! end if
    end if
    if (model .eq. 4) then
    call collisionparticle4(p_rparticleCollection,dnu,dtimestep)
    ! end if
    end if
    if (model .eq. 5) then
    call collisionparticle5(p_rparticleCollection,dnu,dtimestep)
    ! end if
    end if
    

    ! loop over all the particles and move each of them
    ! individually
    do ipart=1,n
      
      p_rgeometryObject => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject
      
      ! position
      dCenterX=p_rgeometryObject%rcoord2D%Dorigin(1)
      dCenterY=p_rgeometryObject%rcoord2D%Dorigin(2)
      
      ! some physical parameters
      dvolume  = (p_rparticleCollection%p_rParticles(ipart)%drad)**2 * SYS_PI
      dmasssl  = p_rparticleCollection%p_rParticles(ipart)%drho * dvolume 
      ddmasssl = (p_rparticleCollection%p_rParticles(ipart)%drho-1.0_dp) * dvolume 
      
      ! Calculate the moment of inertia by formula,
      ! this formula is valid only for a circle
      ! note: for other shaped an approximation can
      ! be calculated by the routine cc_torque:
      ! there we calculated :int_Particle |r|^2 dParticle
      ! the moment of inertia can then be calculated by
      ! int_Particle [(COG_particle - X) .perpdot. u] dParticle / int_Particle |r|^2 dParticle
      ! The part : "int_Particle [(COG_particle - X) .perpdot. u] dParticle"
      ! we already calculated earlier; so we just need to divide...
      ! moment of inertia 
      dimomir  = &
      dmasssl*((p_rparticleCollection%p_rParticles(ipart)%drad)**2)/4.0_dp 
      
      ! Mean resistance forces for the given time step

      ! hydrodynamic drag/lift forces
      dfx = 0.5_dp * (p_rparticleCollection%p_rParticles(ipart)%rResForceX(2) &
       + p_rparticleCollection%p_rParticles(ipart)%rResForceX(1))
      dfy = 0.5_dp * (p_rparticleCollection%p_rParticles(ipart)%rResForceY(2) &
       + p_rparticleCollection%p_rParticles(ipart)%rResForceY(1))

      ! particle collision forces
      dfx1 = p_rparticleCollection%p_rParticles(ipart)%pForceX(1)
      dfy1 = p_rparticleCollection%p_rParticles(ipart)%pForceY(1)
      
      ! Velocity difference for the given time step
      dvelx   = dtimestep*(1.0_dp*dfx+1.0_dp*dfx1+ddmasssl*0.0_dp)/dmasssl
      dvely   = dtimestep*(1.0_dp*dfy+1.0_dp*dfy1+ddmasssl*-9.81_dp)/dmasssl
      
      ! integrate the torque to get the angular velocity

      ! torque due to hydrodynamic forces
      dtrq = 0.5_dp*(p_rparticleCollection%p_rParticles(ipart)%dTorque(2) &
       + p_rparticleCollection%p_rParticles(ipart)%dTorque(1))

      ! torque due to particle collision forces
      dtrq1 = p_rparticleCollection%p_rParticles(ipart)%pTorque(1)
      
      ! angular velocity difference for the given time step
      domega = dtimestep*(dtrq + dtrq1)/dimomir

      !----------------------------------------------------------------------------------------   
      !                          UPDATE THE TORQUE
      !----------------------------------------------------------------------------------------   
      ! set the new values for the torque
      ! torque_old=torque_new
      p_rparticleCollection%p_rParticles(ipart)%dTorque(2)= &
      p_rparticleCollection%p_rParticles(ipart)%dTorque(1)

      p_rparticleCollection%p_rParticles(ipart)%pTorque(1) = 0.0_dp
      
      !----------------------------------------------------------------------------------------   
      !                          UPDATE THE ANGLES AND ANGULAR VELOCITY
      !----------------------------------------------------------------------------------------   
      ! 
      ! a_new = a_old + time*(avel_old + 0.5*avel_new)
      ! Overwrite rotation angle
      p_rgeometryObject%rcoord2D%drotation = &
      p_rgeometryObject%rcoord2D%drotation + &
      dtimestep * (p_rparticleCollection%p_rParticles(ipart)%dangVelocity + 0.5_dp*domega)
      
      ! Recalculate SIN and COS values of angle
      p_rgeometryObject%rcoord2D%dsin_rotation = sin(p_rgeometryObject%rcoord2D%drotation)
      p_rgeometryObject%rcoord2D%dcos_rotation = cos(p_rgeometryObject%rcoord2D%drotation)
      
      ! Update the angular velocity
      ! avel_new = avel_old + avel_new
      p_rparticleCollection%p_rParticles(ipart)%dangVelocity= &
      p_rparticleCollection%p_rParticles(ipart)%dangVelocity + domega
      
      !----------------------------------------------------------------------------------------   
      !                          UPDATE THE FORCES
      !----------------------------------------------------------------------------------------   
      ! save the old forces for the next time step
      p_rparticleCollection%p_rParticles(ipart)%rResForceX(2) = &
      p_rparticleCollection%p_rParticles(ipart)%rResForceX(1)

      p_rparticleCollection%p_rParticles(ipart)%rResForceY(2) = &
      p_rparticleCollection%p_rParticles(ipart)%rResForceY(1)

      p_rparticleCollection%p_rParticles(ipart)%pForceX(1) = 0.0_dp

      p_rparticleCollection%p_rParticles(ipart)%pForceY(1) = 0.0_dp
      
      !----------------------------------------------------------------------------------------   
      !                          UPDATE THE PARTICLE POSITION
      !----------------------------------------------------------------------------------------   
      ! update the position of the x-center
      dCenterX = dCenterX + dtimestep * &
      (p_rparticleCollection%p_rParticles(ipart)%dtransVelX+0.5_dp*dvelx)

      ! update the position of the geometry object
      p_rgeometryObject%rcoord2D%Dorigin(1)=dCenterX
      
      ! update the position of the y-center
       dCenterY = dCenterY + dtimestep * &
      (p_rparticleCollection%p_rParticles(ipart)%dtransVelY+0.5_dp*dvely)
      
      ! update the position of the geometry object      
      p_rgeometryObject%rcoord2D%Dorigin(2)=dCenterY

      !----------------------------------------------------------------------------------------   
      !                          UPDATE THE PARTICLE VELOCITY
      !----------------------------------------------------------------------------------------   
      ! save the current velocity
      p_rparticleCollection%p_rParticles(ipart)%dtransVelX=&
      p_rparticleCollection%p_rParticles(ipart)%dtransVelX + dvelx
      p_rparticleCollection%p_rParticles(ipart)%dtransVelY=&
      p_rparticleCollection%p_rParticles(ipart)%dtransVelY + dvely
      
      !----------------------------------------------------------------------------------------   
      !                          SCREEN OUTPUT
      !----------------------------------------------------------------------------------------   
      ! Output some values to get some readable stuff on the
      ! the screen
      print *,"--------------------------"
      print *,"Total torque force: ",p_rparticleCollection%p_rParticles(ipart)%dTorque(1)
      print *,"--------------------------"
      print *,"domega: ",domega
      print *,"--------------------------"
      print *,"DAngVel: ",p_rparticleCollection%p_rParticles(ipart)%dangVelocity
      print *,"--------------------------"
      print *,"New Position X: ",dCenterX
      print *,"--------------------------"
      print *,"New Position Y: ",dCenterY
      print *,"--------------------------"
      print *,"VelocityYfromFlow: ",dvely
      print *,"--------------------------"
      print *,"VelocityYTotal: ",p_rparticleCollection%p_rParticles(ipart)%dtransVelY
      print *,"--------------------------"
      
    end do ! end ipart
    
  end subroutine
  
! ***************************************************************************   
  
!<subroutine>
  subroutine collisionwall(p_rparticleCollection,dtimestep)
!<description>
  ! 
!</description>

!<inputoutput>
  
  real(dp), intent(IN)     :: dtimestep
  type(t_particleCollection), pointer, intent(INOUT)  :: p_rparticleCollection
  
!</inputoutput>

!</subroutine>

  ! local variables
  integer  :: ipart
  real(DP) :: dCenterX,dCenterY,DistX,DistY,VelX,VelY,rad
  real(DP) :: alpha
  type(t_geometryObject), pointer :: p_rgeometryObject
  
  do ipart = 1,p_rparticleCollection%nparticles

    p_rgeometryObject => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject
    
    ! position
    dCenterX=p_rgeometryObject%rcoord2D%Dorigin(1)
    dCenterY=p_rgeometryObject%rcoord2D%Dorigin(2)

    ! velocity
    VelX=p_rparticleCollection%p_rParticles(ipart)%dtransVelX
    VelY=p_rparticleCollection%p_rParticles(ipart)%dtransVelY

    ! radius
    rad=p_rparticleCollection%p_rParticles(ipart)%drad
    
    ! detect boundary or wall
    DistX = min((abs(dCenterX-0.0_dp)-rad), (abs(dCenterX-2.0_dp)-rad))
    DistY = min((abs(dCenterY-0.0_dp)-rad), (abs(dCenterY-6.0_dp)-rad))
    
    alpha = 0.5_dp
    ! resolve collision of the particle with the boundary
    if (DistX-dtimestep*abs(VelX).le.0.0_dp)then
    p_rparticleCollection%p_rParticles(ipart)%dtransVelX = alpha*(-VelX)
    end if

    if (DistY-dtimestep*abs(VelY).le.0.0_dp)then
    p_rparticleCollection%p_rParticles(ipart)%dtransVelY = alpha*(-VelY)
    end if
    
  end do ! end ipart
    
  end subroutine
  
! ***************************************************************************   
  
!<subroutine>
  subroutine collisionwall_adv(p_rparticleCollection,isegments,p_DbdyEdg,dtimestep)
!<description>
  ! 
!</description>

!<inputoutput>
  integer, intent(IN)      :: isegments
  real(dp), intent(IN)     :: dtimestep
  real(dp), dimension(:,:), pointer, intent(IN)       :: p_DbdyEdg
  type(t_particleCollection), pointer, intent(INOUT)  :: p_rparticleCollection
!</inputoutput>
!</subroutine>
  ! local variables
  integer  :: ipart,iseg
  real(DP) :: dCenterX,dCenterY,DistX,DistY,VelX,VelY,rad
  real(DP) :: alpha
  real(DP), dimension(2)          :: Dp1,Dp2
  type(t_geometryObject), pointer :: p_rgeometryObject
  
!     do iseg=1,isegments
!         Dp1(:)=p_DbdyEdg(:,2*(iseg-1)+1)
!         Dp2(:)=p_DbdyEdg(:,2*(iseg-1)+2)
!     end do
  
  do ipart = 1,p_rparticleCollection%nparticles

    p_rgeometryObject => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject
    
    ! position
    dCenterX=p_rgeometryObject%rcoord2D%Dorigin(1)
    dCenterY=p_rgeometryObject%rcoord2D%Dorigin(2)

    ! velocity
    VelX=p_rparticleCollection%p_rParticles(ipart)%dtransVelX
    VelY=p_rparticleCollection%p_rParticles(ipart)%dtransVelY

    ! radius
    rad=p_rparticleCollection%p_rParticles(ipart)%drad
    
    ! detect boundary or wall
    DistX = min((abs(dCenterX-0.0_dp)-rad), (abs(dCenterX-2.0_dp)-rad))
    DistY = min((abs(dCenterY-0.0_dp)-rad), (abs(dCenterY-6.0_dp)-rad))
    
    alpha = 0.5_dp
    ! resolve collision of the particle with the boundary
    if (DistX-dtimestep*abs(VelX).le.0.0_dp)then
    p_rparticleCollection%p_rParticles(ipart)%dtransVelX = alpha*(-VelX)
    end if

    if (DistY-dtimestep*abs(VelY).le.0.0_dp)then
    p_rparticleCollection%p_rParticles(ipart)%dtransVelY = alpha*(-VelY)
    end if
    
  end do ! end ipart
    
  end subroutine
  
! ***************************************************************************
  
!<subroutine>
  subroutine collisionparticle1(p_rparticleCollection,dnu,dtimestep)
!<description>
  ! 
!</description>

!<inputoutput>
  
  real(dp), intent(IN)                               :: dnu,dtimestep
  type(t_particleCollection), pointer, intent(INOUT) :: p_rparticleCollection
  
!</inputoutput>

!</subroutine>

  ! local variables
  integer  :: ipart,jpart,n
  real(DP) :: fx,fy
  real(DP) :: dCenterXi,dCenterYi,dCenterXj,dCenterYj,DistX,DistY
  real(DP) :: Dist,Dist_exact,VelXi,VelYi,VelXj,VelYj,radi,radj
  real(DP) :: ex,ey,rvx,rvy,g_v,g_v1
  real(DP) :: d0,meu,beta,epsln,d1,epsp,eps_p,lembda
  real(DP) :: kDist,kDist1
  type(t_geometryObject), pointer :: p_rgeometryObjecti,p_rgeometryObjectj
  
  n = p_rparticleCollection%nparticles

  do ipart = 1,n
    do jpart = 1,n

      p_rgeometryObjecti => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject
      p_rgeometryObjectj => p_rparticleCollection%p_rParticles(jpart)%rgeometryObject
      
      ! radius
      radi=p_rparticleCollection%p_rParticles(ipart)%drad

      radj=p_rparticleCollection%p_rParticles(jpart)%drad
      
      ! position
      dCenterXi=p_rgeometryObjecti%rcoord2D%Dorigin(1)
      dCenterYi=p_rgeometryObjecti%rcoord2D%Dorigin(2)

      dCenterXj=p_rgeometryObjectj%rcoord2D%Dorigin(1)
      dCenterYj=p_rgeometryObjectj%rcoord2D%Dorigin(2)
      
      ! velocity
      VelXi=p_rparticleCollection%p_rParticles(ipart)%dtransVelX
      VelYi=p_rparticleCollection%p_rParticles(ipart)%dtransVelY

      VelXj=p_rparticleCollection%p_rParticles(jpart)%dtransVelX
      VelYj=p_rparticleCollection%p_rParticles(jpart)%dtransVelY
      
      if (ipart .ne. jpart) then
      
      ! detect particle collision
      Dist = sqrt((dCenterXj-dCenterXi)**2+(dCenterYj-dCenterYi)**2)

      Dist_exact = Dist-(radi+radj)

      ex = (dCenterXj-dCenterXi)/Dist
      ey = (dCenterYj-dCenterYi)/Dist

      rvx = VelXj-VelXi
      rvy = VelYj-VelYi

      g_v = rvx*ex+rvy*ey
      g_v1 = rvx*(-ey)+rvy*ex
      
      ! ok let's start the collision code (b/w particles)
      d0 = 0.08_dp
      meu = 0.001_dp
      meu = dnu
      epsp = 1e-04_dp
      eps_p = 1e-04_dp
      
      ! calculate forces on ith particle
      
!         if ((0.0_dp .le. Dist_exact) .and. (Dist_exact .le. d0)) then
      if ((0.0_dp .le. Dist_exact) .and. (Dist_exact .le. d0) .and. &
      (g_v .le. 0.0_dp)) then
      fx = 1.0_dp/epsp*(dCenterXi-dCenterXj)*(d0-Dist_exact)**2
      fy = 1.0_dp/epsp*(dCenterYi-dCenterYj)*(d0-Dist_exact)**2

      p_rparticleCollection%p_rParticles(ipart)%pForceX(1) = &
      p_rparticleCollection%p_rParticles(ipart)%pForceX(1) + fx
      p_rparticleCollection%p_rParticles(ipart)%pForceY(1) = &
      p_rparticleCollection%p_rParticles(ipart)%pForceY(1) + fy
      ! end if
      end if

      if ((Dist_exact .lt. 0.0_dp)) then
      fx = 1.0_dp/eps_p*(dCenterXi-dCenterXj)*(-Dist_exact)
      fy = 1.0_dp/eps_p*(dCenterYi-dCenterYj)*(-Dist_exact)

      p_rparticleCollection%p_rParticles(ipart)%pForceX(1) = &
      p_rparticleCollection%p_rParticles(ipart)%pForceX(1) + fx
      p_rparticleCollection%p_rParticles(ipart)%pForceY(1) = &
      p_rparticleCollection%p_rParticles(ipart)%pForceY(1) + fy
      ! no torque due to collision forces (acting normally)
      ! end if
      end if
      
      ! end if
      end if

    end do ! end jpart
  end do ! end ipart
    
  end subroutine
  
! ***************************************************************************
  
!<subroutine>
  subroutine collisionparticle2(p_rparticleCollection,dnu,dtimestep)
!<description>
  ! 
!</description>

!<inputoutput>
  
  real(dp), intent(IN)                               :: dnu,dtimestep
  type(t_particleCollection), pointer, intent(INOUT) :: p_rparticleCollection
  
!</inputoutput>

!</subroutine>

  ! local variables
  integer  :: ipart,jpart,n
  real(DP) :: fx,fy
  real(DP) :: dCenterXi,dCenterYi,dCenterXj,dCenterYj,DistX,DistY
  real(DP) :: Dist,Dist_exact,VelXi,VelYi,VelXj,VelYj,radi,radj
  real(DP) :: ex,ey,rvx,rvy,g_v,g_v1
  real(DP) :: d0,meu,beta,epsln,d1,epsp,eps_p,lembda
  real(DP) :: kDist,kDist1
  type(t_geometryObject), pointer :: p_rgeometryObjecti,p_rgeometryObjectj
  
  n = p_rparticleCollection%nparticles

  do ipart = 1,n
    do jpart = 1,n

      p_rgeometryObjecti => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject
      p_rgeometryObjectj => p_rparticleCollection%p_rParticles(jpart)%rgeometryObject
      
      ! radius
      radi=p_rparticleCollection%p_rParticles(ipart)%drad

      radj=p_rparticleCollection%p_rParticles(jpart)%drad
      
      ! position
      dCenterXi=p_rgeometryObjecti%rcoord2D%Dorigin(1)
      dCenterYi=p_rgeometryObjecti%rcoord2D%Dorigin(2)

      dCenterXj=p_rgeometryObjectj%rcoord2D%Dorigin(1)
      dCenterYj=p_rgeometryObjectj%rcoord2D%Dorigin(2)
      
      ! velocity
      VelXi=p_rparticleCollection%p_rParticles(ipart)%dtransVelX
      VelYi=p_rparticleCollection%p_rParticles(ipart)%dtransVelY

      VelXj=p_rparticleCollection%p_rParticles(jpart)%dtransVelX
      VelYj=p_rparticleCollection%p_rParticles(jpart)%dtransVelY
      
      if (ipart .ne. jpart) then
      
      ! detect particle collision
      Dist = sqrt((dCenterXj-dCenterXi)**2+(dCenterYj-dCenterYi)**2)

      Dist_exact = Dist-(radi+radj)

      ex = (dCenterXj-dCenterXi)/Dist
      ey = (dCenterYj-dCenterYi)/Dist

      rvx = VelXj-VelXi
      rvy = VelYj-VelYi

      g_v = rvx*ex+rvy*ey
      g_v1 = rvx*(-ey)+rvy*ex
      
      ! ok let's start the collision code (b/w particles)
      d0 = 0.08_dp
      meu = 0.001_dp
      meu = dnu
      beta = 1.0_dp/1e-03_dp
      epsln = 1e-05_dp
      d1 = 1e-01_dp
      kDist = 0.0_dp
      kDist1 = 0.0_dp
      
      ! calculate forces on ith particle
!         if ((0.0_dp .le. Dist_exact) .and. (Dist_exact .le. d0)) then
      if ((0.0_dp .le. Dist_exact) .and. (Dist_exact .le. d0) .and. &
      (g_v .le. 0.0_dp)) then
      print *, 'The test is ok! This is case:1.0'
      kDist=beta*meu*6.0_dp*SYS_PI*(1.0_dp/Dist_exact)* &
      (radi**2 * radj**2)/(radi + radj)**2

      kDist1=beta*meu*6.0_dp*SYS_PI*log(d0/Dist_exact)* &
      (radi**2 * radj**2)/(radi + radj)**2
      ! end if
      end if
      if (Dist_exact .lt. 0.0_dp) then
      print *, 'The test is ok! This is case:-1.0'
      kDist=-beta*meu*6.0_dp*SYS_PI*(1.0_dp/Dist_exact)* &
      (radi**2 * radj**2)/(radi + radj)**2

      kDist1=beta*meu*6.0_dp*SYS_PI*log(-d0/Dist_exact)* &
      (radi**2 * radj**2)/(radi + radj)**2
      ! end if
      end if
      
!         if ((epsln .le. Dist_exact) .and. (Dist_exact .le. d0)) then
!         print *, 'The test is ok! This is case:1.0'
!         kDist=beta*meu*6.0_dp*SYS_PI*(1.0_dp/Dist_exact)* &
!         (radi**2 * radj**2)/(radi + radj)**2
! 
!         kDist1=beta*meu*6.0_dp*SYS_PI*log(d0/Dist_exact)* &
!         (radi**2 * radj**2)/(radi + radj)**2
!         ! end if
!         end if
! 
!         if ((0.0_dp .le. Dist_exact) .and. (Dist_exact .lt. epsln)) then
!         print *, 'The test is ok! This is case:1.1'
!         kDist=beta*meu*6.0_dp*SYS_PI*(1.0_dp/d1)* &
!         (radi**2 * radj**2)/(radi + radj)**2
! 
!         kDist1=beta*meu*6.0_dp*SYS_PI*log(d0/d1)* &
!         (radi**2 * radj**2)/(radi + radj)**2
!         ! end if
!         end if
! 
!         if ((-epsln .le. Dist_exact) .and. (Dist_exact .lt. 0.0_dp)) then
!         print *, 'The test is ok! This is case:-1.0'
!         kDist=beta*meu*6.0_dp*SYS_PI*(1.0_dp/d1)* &
!         (radi**2 * radj**2)/(radi + radj)**2
! 
!         kDist1=beta*meu*6.0_dp*SYS_PI*log(d0/d1)* &
!         (radi**2 * radj**2)/(radi + radj)**2
!         ! end if
!         end if
! 
!         if (Dist_exact .lt. -epsln) then
!         print *, 'The test is ok! This is case:-1.1'
!         kDist=-beta*meu*6.0_dp*SYS_PI*(1.0_dp/Dist_exact)* &
!         (radi**2 * radj**2)/(radi + radj)**2
! 
!         kDist1=beta*meu*6.0_dp*SYS_PI*log(-d0/Dist_exact)* &
!         (radi**2 * radj**2)/(radi + radj)**2
!         ! end if
!         end if
      
      fx = (kDist*g_v*ex + kDist1*g_v1*(-ey))
      fy = (kDist*g_v*ey + kDist1*g_v1*ex)

      p_rparticleCollection%p_rParticles(ipart)%pForceX(1) = &
      p_rparticleCollection%p_rParticles(ipart)%pForceX(1) + fx
      p_rparticleCollection%p_rParticles(ipart)%pForceY(1) = &
      p_rparticleCollection%p_rParticles(ipart)%pForceY(1) + fy
      ! torque due to collision forces (not acting normally)
      p_rparticleCollection%p_rParticles(ipart)%pTorque(1) = &
      p_rparticleCollection%p_rParticles(ipart)%pTorque(1) + radi*(ex*fy-ey*fx)
      
      ! end if
      end if

    end do ! end jpart
  end do ! end ipart
    
  end subroutine
  
! ***************************************************************************
  
!<subroutine>
  subroutine collisionparticle3(p_rparticleCollection,dnu,dtimestep)
!<description>
  ! 
!</description>

!<inputoutput>
  
  real(dp), intent(IN)                               :: dnu,dtimestep
  type(t_particleCollection), pointer, intent(INOUT) :: p_rparticleCollection
  
!</inputoutput>

!</subroutine>

  ! local variables
  integer  :: ipart,jpart,i,k,n,neq
  real(DP) :: fx,fy
  real(DP) :: dCenterXi,dCenterYi,dCenterXj,dCenterYj,DistX,DistY
  real(DP) :: Dist,Dist_exact,VelXi,VelYi,VelXj,VelYj,radi,radj
  real(DP) :: ex,ey,rvx,rvy,g_v,g_v1
  real(DP) :: d0,meu,beta,epsln,d1,epsp,eps_p,lembda
  real(DP) :: kDist,kDist1
  real(DP) :: dvelx,dvely,dfx,dfy,dfx1,dfy1,dmasssl,ddmasssl,dvolume
  real(DP), dimension(:), allocatable   :: aprioriVelX,aprioriVelY,b,soln
  real(DP), dimension(:,:), allocatable :: A
  type(t_geometryObject), pointer :: p_rgeometryObjecti,p_rgeometryObjectj
  
  n = p_rparticleCollection%nparticles
  allocate (aprioriVelX(n),aprioriVelY(n))
  aprioriVelX = 0.0_dp
  aprioriVelY = 0.0_dp
  k = 0

  do ipart = 1,n
    do jpart = 1,n

      p_rgeometryObjecti => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject
      p_rgeometryObjectj => p_rparticleCollection%p_rParticles(jpart)%rgeometryObject
      
      ! radius
      radi=p_rparticleCollection%p_rParticles(ipart)%drad

      radj=p_rparticleCollection%p_rParticles(jpart)%drad
      
      ! position
      dCenterXi=p_rgeometryObjecti%rcoord2D%Dorigin(1)
      dCenterYi=p_rgeometryObjecti%rcoord2D%Dorigin(2)

      dCenterXj=p_rgeometryObjectj%rcoord2D%Dorigin(1)
      dCenterYj=p_rgeometryObjectj%rcoord2D%Dorigin(2)
      
      ! velocity
      VelXi=p_rparticleCollection%p_rParticles(ipart)%dtransVelX
      VelYi=p_rparticleCollection%p_rParticles(ipart)%dtransVelY

      VelXj=p_rparticleCollection%p_rParticles(jpart)%dtransVelX
      VelYj=p_rparticleCollection%p_rParticles(jpart)%dtransVelY
      
      if (ipart .ne. jpart) then
      
      ! detect particle collision
      Dist = sqrt((dCenterXj-dCenterXi)**2+(dCenterYj-dCenterYi)**2)

      Dist_exact = Dist-(radi+radj)

      ex = (dCenterXj-dCenterXi)/Dist
      ey = (dCenterYj-dCenterYi)/Dist

      rvx = VelXj-VelXi
      rvy = VelYj-VelYi

      g_v = rvx*ex+rvy*ey
      g_v1 = rvx*(-ey)+rvy*ex
      
      ! ok let's start the collision code (b/w particles)
      d0 = 0.05_dp
      meu = 0.001_dp
      meu = dnu
      
      if ((ipart .lt. jpart) .and. (Dist_exact + dtimestep*g_v .le. 0.0_dp)) then

      if (k .eq. 0) then
      allocate (A(n*(n+3)/2,n*(n+3)/2),b(n*(n+3)/2),soln(n*(n+3)/2))
      A = 0.0_dp
      b = 0.0_dp
      soln = 0.0_dp
      do i = 1,2*n
        A(i,i) = 1.0_dp
      end do ! end i
      ! end if
      end if

      k = k + 1

      lembda = -Dist_exact/dtimestep

      A(2*n+k,ipart) = -ex
      A(2*n+k,jpart) = ex
      A(2*n+k,n+ipart) = -ey
      A(2*n+k,n+jpart) = ey

      ! some physical parameters
      dvolume  = radi**2 * SYS_PI
      dmasssl  = p_rparticleCollection%p_rParticles(ipart)%drho * dvolume
      ddmasssl = (p_rparticleCollection%p_rParticles(ipart)%drho-1.0_dp) * dvolume

      A(ipart,2*n+k) = dtimestep*(ex)/dmasssl
      A(jpart,2*n+k) = dtimestep*(-ex)/dmasssl
      A(n+ipart,2*n+k) = dtimestep*(ey)/dmasssl
      A(n+jpart,2*n+k) = dtimestep*(-ey)/dmasssl

      b(2*n+k) = lembda
      ! end if
      end if

      ! end if
      end if
      
    end do ! end jpart
    
    ! a priori velocity of the particle
    
    dvolume  = radi**2 * SYS_PI
    dmasssl  = p_rparticleCollection%p_rParticles(ipart)%drho * dvolume
    ddmasssl = (p_rparticleCollection%p_rParticles(ipart)%drho-1.0_dp) * dvolume

    ! hydrodynamic drag/lift forces
    dfx = 0.5_dp * (p_rparticleCollection%p_rParticles(ipart)%rResForceX(2) &
     + p_rparticleCollection%p_rParticles(ipart)%rResForceX(1))
    dfy = 0.5_dp * (p_rparticleCollection%p_rParticles(ipart)%rResForceY(2) &
     + p_rparticleCollection%p_rParticles(ipart)%rResForceY(1))
    
    ! Velocity difference for the given time step
    dvelx   = dtimestep*(1.0_dp*dfx+ddmasssl*0.0_dp)/dmasssl
    dvely   = dtimestep*(1.0_dp*dfy+ddmasssl*-9.81_dp)/dmasssl
    
    ! save the current velocity
    aprioriVelX(ipart)=p_rparticleCollection%p_rParticles(ipart)%dtransVelX + dvelx
    aprioriVelY(ipart)=p_rparticleCollection%p_rParticles(ipart)%dtransVelY + dvely
    
  end do ! end ipart
  
  if (k .ne. 0) then
  do ipart = 1,n
    b(ipart) = aprioriVelX(ipart)
    b(n+ipart) = aprioriVelY(ipart)
  end do ! end ipart

  neq = 2*n + k

  ! calculate new velocities of the particle
  call gauss_pivot (A,b,soln,neq)

  ! now we need the forces due to particle collision, so
  ! calculate forces on ith particle
  do ipart = 1,n
    
    dvolume  = (p_rparticleCollection%p_rParticles(ipart)%drad)**2 * SYS_PI
    dmasssl  = p_rparticleCollection%p_rParticles(ipart)%drho * dvolume
    ddmasssl = (p_rparticleCollection%p_rParticles(ipart)%drho-1.0_dp) * dvolume

    ! hydrodynamic drag/lift forces
    dfx = 0.5_dp * (p_rparticleCollection%p_rParticles(ipart)%rResForceX(2) &
     + p_rparticleCollection%p_rParticles(ipart)%rResForceX(1))
    dfy = 0.5_dp * (p_rparticleCollection%p_rParticles(ipart)%rResForceY(2) &
     + p_rparticleCollection%p_rParticles(ipart)%rResForceY(1))

    ! force of gravity + hydrodynamic force
    fx = dfx + ddmasssl*0.0_dp
    fy = dfy + ddmasssl*-9.81_dp
    
    ! force = mass*(new_vel - old_vel)/timestep
    ! force due to particle collision = force - (hydrodynamic force + gravity force)
    p_rparticleCollection%p_rParticles(ipart)%pForceX(1) = dmasssl*(soln(ipart)- &
    p_rparticleCollection%p_rParticles(ipart)%dtransVelX)/dtimestep - fx
    p_rparticleCollection%p_rParticles(ipart)%pForceY(1) = dmasssl*(soln(n+ipart)- &
    p_rparticleCollection%p_rParticles(ipart)%dtransVelY)/dtimestep - fy
    ! no torque due to collision forces (acting normally)
    
  end do ! end ipart

  deallocate (A,b,soln)
  ! end if
  end if
  
  deallocate (aprioriVelX,aprioriVelY)
    
  end subroutine
  
! ***************************************************************************
  
!<subroutine>
  subroutine collisionparticle4(p_rparticleCollection,dnu,dtimestep)
!<description>
  ! 
!</description>

!<inputoutput>
  
  real(dp), intent(IN)                               :: dnu,dtimestep
  type(t_particleCollection), pointer, intent(INOUT) :: p_rparticleCollection
  
!</inputoutput>

!</subroutine>

  ! local variables
  integer  :: ipart,jpart,i,k,n,neq,m,s
  real(DP) :: fx,fy
  real(DP) :: dCenterXi,dCenterYi,dCenterXj,dCenterYj,DistX,DistY
  real(DP) :: Dist,Dist_exact,VelXi,VelYi,VelXj,VelYj,radi,radj
  real(DP) :: ex,ey,rvx,rvy,g_v,g_v1
  real(DP) :: d0,meu,beta,epsln,d1,epsp,eps_p,lembda
  real(DP) :: kDist,kDist1
  real(DP) :: dvelx,dvely,dfx,dfy,dfx1,dfy1,dmasssl,ddmasssl,dvolume
  real(DP) :: dmasssli,dvolumei,dmassslj,dvolumej
  real(DP), dimension(:), allocatable   :: aprioriVelX,aprioriVelY,b,soln
  real(DP), dimension(:,:), allocatable :: A,gamma1
  type(t_geometryObject), pointer :: p_rgeometryObjecti,p_rgeometryObjectj
  
  n = p_rparticleCollection%nparticles
  allocate (aprioriVelX(n),aprioriVelY(n))
  aprioriVelX = 0.0_dp
  aprioriVelY = 0.0_dp
  k = 0

  do ipart = 1,n
    do jpart = 1,n

      p_rgeometryObjecti => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject
      p_rgeometryObjectj => p_rparticleCollection%p_rParticles(jpart)%rgeometryObject
      
      ! radius
      radi=p_rparticleCollection%p_rParticles(ipart)%drad

      radj=p_rparticleCollection%p_rParticles(jpart)%drad
      
      ! position
      dCenterXi=p_rgeometryObjecti%rcoord2D%Dorigin(1)
      dCenterYi=p_rgeometryObjecti%rcoord2D%Dorigin(2)

      dCenterXj=p_rgeometryObjectj%rcoord2D%Dorigin(1)
      dCenterYj=p_rgeometryObjectj%rcoord2D%Dorigin(2)
      
      ! velocity
      VelXi=p_rparticleCollection%p_rParticles(ipart)%dtransVelX
      VelYi=p_rparticleCollection%p_rParticles(ipart)%dtransVelY

      VelXj=p_rparticleCollection%p_rParticles(jpart)%dtransVelX
      VelYj=p_rparticleCollection%p_rParticles(jpart)%dtransVelY
      
      if (ipart .ne. jpart) then
      
      ! detect particle collision
      Dist = sqrt((dCenterXj-dCenterXi)**2+(dCenterYj-dCenterYi)**2)

      Dist_exact = Dist-(radi+radj)

      ex = (dCenterXj-dCenterXi)/Dist
      ey = (dCenterYj-dCenterYi)/Dist

      rvx = VelXj-VelXi
      rvy = VelYj-VelYi

      g_v = rvx*ex+rvy*ey
      g_v1 = rvx*(-ey)+rvy*ex
      
      ! ok let's start the collision code (b/w particles)
      d0 = 0.05_dp
      meu = 0.001_dp
      meu = dnu
      
      if ((ipart .lt. jpart) .and. (Dist_exact + dtimestep*g_v .le. 0.0_dp)) then

      if (k .eq. 0) then
      allocate (A(n*(n+3)/2,n*(n+3)/2),b(n*(n+3)/2),soln(n*(n+3)/2),gamma1(n,n))
      A = 0.0_dp
      b = 0.0_dp
      soln = 0.0_dp
      gamma1 = 0.0_dp
      do i = 1,2*n
        A(i,i) = 1.0_dp
      end do ! end i
      ! end if
      end if

      k = k + 1

      lembda = -Dist_exact/dtimestep

      A(2*n+k,ipart) = -ex
      A(2*n+k,jpart) = ex
      A(2*n+k,n+ipart) = -ey
      A(2*n+k,n+jpart) = ey

      ! some physical parameters
      dvolume  = radi**2 * SYS_PI
      dmasssl  = p_rparticleCollection%p_rParticles(ipart)%drho * dvolume
      ddmasssl = (p_rparticleCollection%p_rParticles(ipart)%drho-1.0_dp) * dvolume

      A(ipart,2*n+k) = dtimestep*(ex)/dmasssl
      A(jpart,2*n+k) = dtimestep*(-ex)/dmasssl
      A(n+ipart,2*n+k) = dtimestep*(ey)/dmasssl
      A(n+jpart,2*n+k) = dtimestep*(-ey)/dmasssl

      b(2*n+k) = lembda
      ! end if
      end if

      ! end if
      end if
      
    end do ! end jpart
    
    ! a priori velocity of the particle
    
    dvolume  = radi**2 * SYS_PI
    dmasssl  = p_rparticleCollection%p_rParticles(ipart)%drho * dvolume
    ddmasssl = (p_rparticleCollection%p_rParticles(ipart)%drho-1.0_dp) * dvolume

    ! hydrodynamic drag/lift forces
    dfx = 0.5_dp * (p_rparticleCollection%p_rParticles(ipart)%rResForceX(2) &
     + p_rparticleCollection%p_rParticles(ipart)%rResForceX(1))
    dfy = 0.5_dp * (p_rparticleCollection%p_rParticles(ipart)%rResForceY(2) &
     + p_rparticleCollection%p_rParticles(ipart)%rResForceY(1))
    
    ! Velocity difference for the given time step
    dvelx   = dtimestep*(1.0_dp*dfx+ddmasssl*0.0_dp)/dmasssl
    dvely   = dtimestep*(1.0_dp*dfy+ddmasssl*-9.81_dp)/dmasssl
    
    ! save the current velocity
    aprioriVelX(ipart)=p_rparticleCollection%p_rParticles(ipart)%dtransVelX + dvelx
    aprioriVelY(ipart)=p_rparticleCollection%p_rParticles(ipart)%dtransVelY + dvely
    
  end do ! end ipart
  
  if (k .ne. 0) then
  do ipart = 1,n
    b(ipart) = aprioriVelX(ipart)
    b(n+ipart) = aprioriVelY(ipart)
  end do ! end ipart

  neq = 2*n + k

  ! calculate new velocities of the particle
  call gauss_pivot (A,b,soln,neq)

!     ! case for gluey particle
!     m = 0
!     s = 0
!     beta = 0.5_dp ! parameter to restrict the magnitude of force (beta > 0)
!     do ipart = 1,n
!       do jpart = 1,n
!         if (ipart .lt. jpart) then
!         m = m + 1
!         if (gamma1(ipart,jpart) .eq. 1) then
!         s = s+1
!         gamma(m) = gamma(m) - beta*dt*soln(2*n+s)
! !         if (gamma(m) .gt. 0.0_dp)
!         if (gamma(m) .lt. 0.0_dp) then
!         dvolumei  = (p_rparticleCollection%p_rParticles(ipart)%drad)**2 * SYS_PI
!         dmasssli  = p_rparticleCollection%p_rParticles(ipart)%drho * dvolumei
!         dvolumej  = (p_rparticleCollection%p_rParticles(jpart)%drad)**2 * SYS_PI
!         dmassslj  = p_rparticleCollection%p_rParticles(jpart)%drho * dvolumej
! 
!         soln(ipart) = gamma(m)/dmasssli
!         soln(jpart) = gamma(m)/dmassslj
!         soln(n+ipart) = gamma(m)/dmasssli
!         soln(n+jpart) = gamma(m)/dmassslj
!         gamma(m) = 0.0_dp
! !         gamma(m) = 0.01_dp
!         ! end if
!         end if
!         ! end if
!         end if
!         ! end if
!         end if
!       end do ! end jpart
!     end do ! end ipart

  ! now we need the forces due to particle collision, so
  ! calculate forces on ith particle
  do ipart = 1,n
    
    dvolume  = (p_rparticleCollection%p_rParticles(ipart)%drad)**2 * SYS_PI
    dmasssl  = p_rparticleCollection%p_rParticles(ipart)%drho * dvolume
    ddmasssl = (p_rparticleCollection%p_rParticles(ipart)%drho-1.0_dp) * dvolume

    ! hydrodynamic drag/lift forces
    dfx = 0.5_dp * (p_rparticleCollection%p_rParticles(ipart)%rResForceX(2) &
     + p_rparticleCollection%p_rParticles(ipart)%rResForceX(1))
    dfy = 0.5_dp * (p_rparticleCollection%p_rParticles(ipart)%rResForceY(2) &
     + p_rparticleCollection%p_rParticles(ipart)%rResForceY(1))

    ! force of gravity + hydrodynamic force
    fx = dfx + ddmasssl*0.0_dp
    fy = dfy + ddmasssl*-9.81_dp
    
    ! force = mass*(new_vel - old_vel)/timestep
    ! force due to particle collision = force - (hydrodynamic force + gravity force)
    p_rparticleCollection%p_rParticles(ipart)%pForceX(1) = dmasssl*(soln(ipart)- &
    p_rparticleCollection%p_rParticles(ipart)%dtransVelX)/dtimestep - fx
    p_rparticleCollection%p_rParticles(ipart)%pForceY(1) = dmasssl*(soln(n+ipart)- &
    p_rparticleCollection%p_rParticles(ipart)%dtransVelY)/dtimestep - fy
    ! no torque due to collision forces (acting normally)
    
  end do ! end ipart

  deallocate (A,b,soln,gamma1)
  ! end if
  end if
  
  deallocate (aprioriVelX,aprioriVelY)
    
  end subroutine
  
! ***************************************************************************
  
!<subroutine>
  subroutine collisionparticle5(p_rparticleCollection,dnu,dtimestep)
!<description>
  ! 
!</description>

!<inputoutput>
  real(dp), intent(IN)                               :: dnu,dtimestep
  type(t_particleCollection), pointer, intent(INOUT) :: p_rparticleCollection
!</inputoutput>
!</subroutine>
  ! local variables
  integer  :: ipart,jpart,n
  real(DP) :: fx,fy,fx1,fy1,fx2,fy2
  real(DP) :: dCenterXi,dCenterYi,dCenterXj,dCenterYj,DistX,DistY
  real(DP) :: Dist,Dist_exact,VelXi,VelYi,VelXj,VelYj,radi,radj
  real(DP) :: ex,ey,rvx,rvy,g_v,g_v1
  real(DP) :: d0,meu,beta,epsln,d1,epsp,eps_p,lembda
  real(DP) :: kDist,kDist1
  real(DP) :: dmasssli,dvolumei,dmassslj,dvolumej
  real(DP) :: cr,u1n,u2n,u1,u2,u1x,u1y,u2x,u2y
  type(t_geometryObject), pointer :: p_rgeometryObjecti,p_rgeometryObjectj
  
  n = p_rparticleCollection%nparticles

  do ipart = 1,n
    do jpart = 1,n

      p_rgeometryObjecti => p_rparticleCollection%p_rParticles(ipart)%rgeometryObject
      p_rgeometryObjectj => p_rparticleCollection%p_rParticles(jpart)%rgeometryObject
      
      ! radius
      radi=p_rparticleCollection%p_rParticles(ipart)%drad

      radj=p_rparticleCollection%p_rParticles(jpart)%drad
      
      ! position
      dCenterXi=p_rgeometryObjecti%rcoord2D%Dorigin(1)
      dCenterYi=p_rgeometryObjecti%rcoord2D%Dorigin(2)

      dCenterXj=p_rgeometryObjectj%rcoord2D%Dorigin(1)
      dCenterYj=p_rgeometryObjectj%rcoord2D%Dorigin(2)
      
      ! velocity
      VelXi=p_rparticleCollection%p_rParticles(ipart)%dtransVelX
      VelYi=p_rparticleCollection%p_rParticles(ipart)%dtransVelY

      VelXj=p_rparticleCollection%p_rParticles(jpart)%dtransVelX
      VelYj=p_rparticleCollection%p_rParticles(jpart)%dtransVelY
      
      if (ipart .ne. jpart) then
      
      ! detect particle collision
      Dist = sqrt((dCenterXj-dCenterXi)**2+(dCenterYj-dCenterYi)**2)

      Dist_exact = Dist-(radi+radj)

      ex = (dCenterXj-dCenterXi)/Dist
      ey = (dCenterYj-dCenterYi)/Dist

      rvx = VelXj-VelXi
      rvy = VelYj-VelYi

      g_v = rvx*ex+rvy*ey
      g_v1 = rvx*(-ey)+rvy*ex
      
      ! ok let's start the collision code (b/w particles)
      d0 = 0.05_dp
      meu = 0.001_dp
      meu = dnu
      cr = 0.97_dp
      
      ! calculate forces on ith particle
      
      if (Dist_exact + dtimestep*g_v .le. 0.0_dp) then
!         if ((0.0_dp .le. Dist_exact) .and. (Dist_exact .le. d0)) then
      ! some physical parameters
      dvolumei  = radi**2 * SYS_PI
      dmasssli  = p_rparticleCollection%p_rParticles(ipart)%drho * dvolumei
      dvolumej  = radj**2 * SYS_PI
      dmassslj  = p_rparticleCollection%p_rParticles(jpart)%drho * dvolumej

      ! normal comp. of velocity
      u1n = VelXi*ex+VelYi*ey
      u2n = VelXj*ex+VelYj*ey

      u1 = cr*g_v*dmassslj/(dmasssli+dmassslj) + &
      (dmasssli*u1n + dmassslj*u2n)/(dmasssli+dmassslj)
      u2 = cr*g_v*dmasssli/(dmasssli+dmassslj) + &
      (dmasssli*u1n + dmassslj*u2n)/(dmasssli+dmassslj)

      u1x = u1*ex
      u1y = u1*ey
      u2x = -u2*ex
      u2y = -u2*ey

      fx1 = dmasssli*(u1x-VelXi)/dtimestep
      fy1 = dmasssli*(u1y-VelYi)/dtimestep
      fx2 = dmassslj*(u2x-VelXj)/dtimestep
      fy2 = dmassslj*(u2y-VelYj)/dtimestep

      p_rparticleCollection%p_rParticles(ipart)%pForceX(1) = &
      p_rparticleCollection%p_rParticles(ipart)%pForceX(1) + fx1
      p_rparticleCollection%p_rParticles(ipart)%pForceY(1) = &
      p_rparticleCollection%p_rParticles(ipart)%pForceY(1) + fy1
      p_rparticleCollection%p_rParticles(jpart)%pForceX(1) = &
      p_rparticleCollection%p_rParticles(jpart)%pForceX(1) + fx2
      p_rparticleCollection%p_rParticles(jpart)%pForceY(1) = &
      p_rparticleCollection%p_rParticles(jpart)%pForceY(1) + fy2
      ! no torque due to collision forces (acting normally)
      ! end if
      end if
      
      ! end if
      end if

    end do ! end jpart
  end do ! end ipart
    
  end subroutine
  
! ***************************************************************************
  
!<subroutine>

  ! subroutine gauss_pivot (A,b,soln,neq,error)

  subroutine gauss_pivot (A,b,soln,neq)

  !<description>
  ! This subroutine sloves system of n linear equations
  ! using Gaussian elimination and the maximum pivot technique


    ! <inout>
    integer, intent(IN) :: neq
    real(DP), dimension(:,:), intent(IN)     :: A
    real(DP), dimension(:), intent(IN)       :: b
    real(DP), dimension(:), intent(OUT)      :: soln
    ! integer, intent(OUT)                   :: error

    ! <local>
    real(DP), parameter                      :: epsilon = 1.0e-6
    real(DP), dimension(:,:), allocatable    :: A1
    real(DP), dimension(:), allocatable      :: temp1
    real(DP)                                 :: factor,temp
    integer                                  :: irow,ipeak,jrow

    allocate (A1(neq,neq), temp1(neq))
    A1 = A
    soln = b

    do irow = 1,neq

      ! find peak pivot for column irow in rows irow to neq
      ipeak = irow
      do jrow = irow+1,neq
        if (abs(A1(jrow,irow)) .gt. abs(A1(ipeak,irow))) then
        ipeak = jrow
        end if
      end do ! end jrow

      ! check for singular equations
      if (abs(A1(ipeak,irow)) .lt. epsilon) then
      ! error = 1
      return
      end if

      ! swap equations irow & ipeak
      if (ipeak .ne. irow) then
      temp1(:) = A1(ipeak,:)
      A1(ipeak,:) = A1(irow,:)
      A1(irow,:) = temp1(:)
      temp = soln(ipeak)
      soln(ipeak) = soln(irow)
      soln(irow) = temp
      end if

      ! make zeros for the off-diagonal entries
      do jrow = 1,neq
        if (jrow .ne. irow) then
        factor = -A1(jrow,irow)/A1(irow,irow)
        A1(jrow,:) = A1(irow,:)*factor + A1(jrow,:)
        soln(jrow) = soln(irow)*factor + soln(jrow)
        end if
      end do ! end jrow
    end do ! end irow

    ! divide by the diagonal entry to get final solution
    do irow = 1,neq
      soln(irow) = soln(irow)/A1(irow,irow)
      A1(irow,irow) = 1.0_dp
    end do ! end irow

    ! set error flag to 0 and return
    ! error = 0

    deallocate (A1, temp1)
    
  end subroutine gauss_pivot

! ***************************************************************************  
  
!<subroutine>
  subroutine collisionWall_normal(p_rparticleCollection,dnu,dtimestep)
!<description>
  ! 
!</description>

!<input>
  ! the viscosity
  real(dp), intent(in)                               :: dnu
  ! the timestepsize, that was used to calculate the current solution
  real(dp), intent(in)                               :: dtimestep  
!</input>  

!<inputoutput>  
  type(t_particleCollection), pointer, intent(inout) :: p_rparticleCollection
!</inputoutput>

!</subroutine>
  ! local variables
  integer  :: ipart,jpart,i,k,n,neq,m,s
  real(DP) :: fx,fy
  real(DP) :: dCenterXi,dCenterYi,dCenterXj,dCenterYj,DistX,DistY
  real(DP) :: Dist,Dist_exact,VelXi,VelYi,VelXj,VelYj,radi,radj
  real(DP) :: ex,ey,rvx,rvy,g_v,g_v1
  real(DP) :: d0,meu,beta,epsln,d1,epsp,eps_p,lembda
  real(DP) :: kDist,kDist1
  real(DP) :: dvelx,dvely,dfx,dfy,dfx1,dfy1,dmasssl,ddmasssl,dvolume
  real(DP) :: dmasssli,dvolumei,dmassslj,dvolumej
  type(t_geometryObject), pointer :: p_rgeometryObjecti,p_rgeometryObjectj

  end subroutine ! collisionWall_normal
  
! ***************************************************************************  
  
!<subroutine>
  subroutine collisionParticle_normal(p_rparticleCollection,dnu,dtimestep)
!<description>
  ! 
!</description>

!<input>
  ! the viscosity
  real(dp), intent(in)                               :: dnu
  ! the timestepsize, that was used to calculate the current solution
  real(dp), intent(in)                               :: dtimestep  
!</input>  

!<inputoutput>  
  type(t_particleCollection), pointer, intent(inout) :: p_rparticleCollection
!</inputoutput>

!</subroutine>
  ! local variables
  integer  :: ipart,jpart,i,k,n,neq,m,s
  real(DP) :: fx,fy
  real(DP) :: dCenterXi,dCenterYi,dCenterXj,dCenterYj,DistX,DistY
  real(DP) :: Dist,Dist_exact,VelXi,VelYi,VelXj,VelYj,radi,radj
  real(DP) :: ex,ey,rvx,rvy,g_v,g_v1
  real(DP) :: d0,meu,beta,epsln,d1,epsp,eps_p,lembda
  real(DP) :: kDist,kDist1
  real(DP) :: dvelx,dvely,dfx,dfy,dfx1,dfy1,dmasssl,ddmasssl,dvolume
  real(DP) :: dmasssli,dvolumei,dmassslj,dvolumej
  type(t_geometryObject), pointer :: p_rgeometryObjecti,p_rgeometryObjectj

  end subroutine ! collisionParticle_normal

! ***************************************************************************
end module

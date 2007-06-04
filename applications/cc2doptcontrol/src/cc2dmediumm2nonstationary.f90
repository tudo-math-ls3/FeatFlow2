!##############################################################################
!# ****************************************************************************
!# <name> cc2dmediumm2nonstationary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises a fully coupled time dependent solver for the 
!# coupled Navier-Stokes optimal control problem.
!#
!# The following routines can be found here:
!#
!# 1.) c2d2_initParTimeDependence
!#     -> Initialise the parameters of the time dependent solver from DAT file
!#        parameters.
!#
!# 2.) c2d2_solveNonstationary
!#     -> Solves the space-time coupled system.
!# </purpose>
!##############################################################################

MODULE cc2dmediumm2nonstationary

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE coarsegridcorrection
  USE spdiscprojection
  USE nonlinearsolver
  USE paramlist
  USE linearsolverautoinitialise
  USE matrixrestriction
  USE paramlist
  USE timestepping
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmedium_callback

  USE cc2dmediumm2nonlinearcore
  USE cc2dmediumm2nonlinearcoreinit
  USE cc2dmediumm2stationary
  USE adaptivetimestep
  USE cc2dmediumm2timeanalysis
  USE cc2dmediumm2boundary
  USE cc2dmediumm2discretisation
  USE cc2dmediumm2postprocessing
  USE cc2dmediumm2timesupersystem
    
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE c2d2_initParTimeDependence (rproblem,ssection,rparams)
  
!<description>
  ! Initialises parameters in the problem structure according to whether the
  ! simulation is time dependent or not. Initialises the time stepping scheme,
  ! start proceduce, error bounds, etc.
  !
  ! Note: This will not allocate an memory but only initialise the parameters
  ! in rproblem according to the parameters in rparams from the DAT file.
!</description>
  
!<input>
  ! A parameter list from a DAT file containing parameters that configure the
  ! time dependence.
  TYPE(t_parlist), INTENT(IN) :: rparams
  
  ! The name of the section in the parameter list containing the parameters
  ! for the time dependent simulation.
  CHARACTER(LEN=*), INTENT(IN) :: ssection
!</input>

!<inputoutput>
  ! The problem structure to be initialised.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>
  
!</subroutine>

    ! Fetch the parameters. Initialise with standard settings if they don't 
    ! exist.
    CALL parlst_getvalue_int (rparams,ssection,'itimedependence',  &
        rproblem%itimedependence, 0)
    CALL parlst_getvalue_int (rparams,ssection,'niterations',     &
        rproblem%rtimedependence%niterations, 1000)
    CALL parlst_getvalue_double (rparams,ssection,'dtimestart',   &
        rproblem%rtimedependence%dtimeInit, 0.0_DP)
    CALL parlst_getvalue_double (rparams,ssection,'dtimemax',     &
        rproblem%rtimedependence%dtimeMax, 20.0_DP)
    CALL parlst_getvalue_int (rparams, ssection,'ctimeStepScheme', &
          rproblem%rtimedependence%ctimeStepScheme, 0)
    CALL parlst_getvalue_double (rparams,ssection,'dtimeStepTheta',    &
        rproblem%rtimedependence%dtimeStepTheta, 1.0_DP)
    rproblem%rtimedependence%itimeStep = 0

  END SUBROUTINE  

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_solveNonstationaryDirect (rproblem,rvectorTmp,rrhsTmp)
  
!<description>
  ! Solve the nonstationary optimal control problem. This allocates
  ! memory for the solution vector and solves the space-time coupled
  ! system.
  ! This is the 'single-grid' variant of the coupled solver, i.e.
  ! it calls the solver once without any multigrid schemes in space/time.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
  
  ! Temporary vector which defines the shape of the solution vector
  ! in each time step.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvectorTmp

  ! Temporary vector which defines the shape of the RHS vector
  ! in each time step.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rrhsTmp
!</inputoutput>

!</subroutine>

    TYPE(t_ccoptSpaceTimeDiscretisation) :: rsupermatrix
    TYPE(t_spacetimeVector) :: rx,rd
    TYPE(t_vectorBlock) :: rtempVector
    INTEGER :: i

    ! Initialise the supersystem on the maximum level
    CALL c2d2_initParamsSupersystem (rproblem,1,rproblem%NLMAX,&
        rsupermatrix, rx, rd)
        
    ! Allocate memory for the 3rd temp vector
    CALL lsysbl_createVecBlockIndirect (rrhsTmp,rtempVector,.FALSE.)
    
    ! Call the solver for the space/time coupled system. We only solve on level NLMAX.
    CALL c2d2_solveSupersystem (rproblem, rsupermatrix, rx, rd, &
      rvectorTmp, rrhsTmp, rtempVector)
      
    ! Postprocessing of all solution vectors.
    DO i = 0,rsupermatrix%niterations
    
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + i*rsupermatrix%dtstep
      rproblem%rtimedependence%itimeStep = i
    
      CALL sptivec_getTimestepData (rx, i, rvectorTmp)
    
      CALL c2d2_postprocessingNonstat (rproblem,rvectorTmp)  
      
    END DO
    
    ! Release memory, finish.
    CALL lsysbl_releaseVector (rtempVector)
    
    CALL c2d2_doneParamsSupersystem (rsupermatrix,rx,rd)

  END SUBROUTINE
  
END MODULE

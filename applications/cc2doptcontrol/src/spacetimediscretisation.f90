!##############################################################################
!# ****************************************************************************
!# <name> spacetimediscretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the space time discretisation.
!# The structure t_ccoptSpaceTimeDiscretisation defines the general
!# discretisation of the underlying problem in space and time and can be
!# initialised and cleaned up with the routines in this module. The following
!# routines can be found here:
!#
!# 1.) sptidis_initDiscretisation
!#     -> Initialise a space-time discretisation structure.
!#
!# 2.) sptidis_doneDiscretisation
!#     -> Clean up a space-time discretisation structure.
!##############################################################################

MODULE spacetimediscretisation

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
  USE l2projection
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic

  IMPLICIT NONE

!<types>

!<typeblock>

  ! Defines the basic shape of the supersystem which realises the coupling
  ! between all timesteps.
  TYPE t_ccoptSpaceTimeDiscretisation
  
    ! Spatial refinement level of this matrix.
    INTEGER :: ilevel = 0
  
    ! Number of time steps
    !INTEGER :: niterations         = 0

    ! Absolute start time of the simulation
    !REAL(DP) :: dtimeInit          = 0.0_DP     
    
    ! Maximum time of the simulation
    !REAL(DP) :: dtimeMax           = 1.0_DP
    
    ! Time step length of the time discretisation
    !REAL(DP) :: dtstep

    ! Defines the discretisation in time
    TYPE(t_timeDiscretisation) :: rtimeDiscr
    
    ! Number of time-DOF's. For Theta-schemes, this coincides with
    ! the number of intervals, but for more complicated schemes,
    ! this may differ.
    INTEGER :: NEQtime = 0

    ! Time-stepping scheme;
    ! 0=one step FD scheme (Euler, CN)
    ! 2=dG(0)
    !INTEGER :: ctimeStepScheme     = 0
    
    ! Parameter for one step scheme (THETA) if itimeStepScheme=0;
    ! =0:Forward Euler(instable), =1: Backward Euler, =0.5: Crank-Nicolson
    !REAL(DP) :: dtimeStepTheta     = 1.0_DP
    
    ! Regularisation parameter for the control $\alpha$. Must be <> 0.
    ! A value of 0.0 disables the terminal condition.
    REAL(DP) :: dalphaC = 1.0_DP
    
    ! Regularisation parameter for the terminal condition 
    ! $\gamma/2*||y(T)-z(T)||$.
    ! A value of 0.0 disables the terminal condition.
    REAL(DP) :: dgammaC = 0.0_DP
    
    ! Type of implementation of the terminal condition.
    ! =0: implement terminal condition in a weak sense by filtering.
    ! =1: implement terminal condition in a strong sense by modifying the matrix.
    INTEGER :: itypeTerminalCondition = 0

    ! Problem-related structure that provides the templates for
    ! matrices/vectors on the spatial level of the matrix.
    TYPE(t_problem_lvl), POINTER :: p_rlevelInfo

  END TYPE

!</typeblock>

!</types>

CONTAINS

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE sptidis_initDiscretisation (rproblem,ilevelTime,ilevelSpace,&
      rspaceTimeDiscr)
  
!<description>
  ! Initialises a space time discretisation structure according to the 
  ! parameters in the DAT file and the problem structure.
!</description>

!<input>
  ! The problem structure describing the whole problem.
  TYPE(t_problem), INTENT(IN), TARGET :: rproblem
  
  ! 'Refinement level in time'. =1: Initialise as described in
  ! rproblem. >1: Refine (irefineInTime-1) times regularly in time.
  ! (-> #timesteps * 2**(irefineInTime-1) )
  INTEGER, INTENT(IN) :: ilevelTime

  ! 'Refinement level in space'. =1: Calculate on the coarse mesh
  ! >0: calculate on level ilevelSpace. Must be <= rproblem%NLMAX,
  ! >= rproblem%NLMIN.
  INTEGER, INTENT(IN) :: ilevelSpace
!</input>

!<inputoutput>
  ! Supersystem-structure to be initialised.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(OUT) :: rspaceTimeDiscr
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: niterations

    ! Copy most relevant data from the problem structure.
    rspaceTimeDiscr%ilevel          = ilevelSpace
    rspaceTimeDiscr%p_rlevelInfo    => rproblem%RlevelInfo(ilevelSpace)
    niterations = rproblem%rtimedependence%niterations * 2**MAX(0,ilevelTime-1)
    
    ! Initialise the time discretisation
    SELECT CASE (rproblem%rtimedependence%ctimeStepScheme)
    CASE (0) 
      CALL tdiscr_initTheta(&
          rproblem%rtimedependence%dtimeInit,&
          rproblem%rtimedependence%dtimeMax,&
          niterations,rproblem%rtimedependence%dtimeStepTheta,&
          rspaceTimeDiscr%rtimeDiscr)
    CASE (2)
      CALL tdiscr_initdG0(rproblem%rtimedependence%dtimeInit,&
          rproblem%rtimedependence%dtimeMax,&
          niterations, rspaceTimeDiscr%rtimeDiscr)
    CASE DEFAULT
      PRINT *,'c2d2_initParamsSupersystem: Unsupported time discretisation.'
      CALL sys_halt()
    END SELECT
    
    ! Save the number of time-DOF's
    rspaceTimeDiscr%NEQtime = tdiscr_igetNDofGlob(rspaceTimeDiscr%rtimeDiscr)
    
    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dalphaC',rspaceTimeDiscr%dalphaC,1.0_DP)
    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dgammaC',rspaceTimeDiscr%dgammaC,0.0_DP)
    CALL parlst_getvalue_int (rproblem%rparamList,'OPTIMALCONTROL',&
        'itypeTerminalCondition',rspaceTimeDiscr%itypeTerminalCondition,0)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE sptidis_doneDiscretisation (rspaceTimeDiscr)
  
!<description>
  ! Cleans up a given space time discrtisation structure.
!</description>

!<inputoutput>
  ! Supersystem-structure to be cleaned up.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(OUT) :: rspaceTimeDiscr
!</inputoutput>

!</subroutine>

    ! Currently there is nothing to do!

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE sptidis_infoDiscretisation (rspaceTimeDiscr)
  
!<description>
  ! Prints statistical information about a time discretisation.
!</description>

!<inputoutput>
  ! Supersystem-structure to be cleaned up.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
!</inputoutput>

!</subroutine>
    
    CALL output_line ('NEQ in time   : '//sys_siL(rspaceTimeDiscr%NEQTime,10))
    CALL output_line ('NEQ in space  : '//&
      sys_siL(dof_igetNDofGlobBlock(rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation),10))

  END SUBROUTINE

END MODULE

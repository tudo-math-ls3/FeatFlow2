!##############################################################################
!# ****************************************************************************
!# <name> ccmainproblem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a program how to solve a stationary or nonstationary
!# Navier-Stokes problem 
!#
!#              $$- \nu Laplace(u) + u*grad(u) + \Nabla p = f $$
!#              $$ \Nable \cdot p = 0$$
!#
!# on a 2D domain for a 2D function $u=(u_1,u_2)$ and a pressure $p$.
!#
!# The routine splits up the tasks of reading the domain, creating 
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines. For the nonlinearity, the nonlinear solver is invoked.
!# </purpose>
!##############################################################################

MODULE ccmainproblem

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
  USE statistics
  
  USE collection
  USE convection
    
  USE ccbasic
  USE ccinitgeneralparameters
  USE ccinitparamtriang
  USE ccboundarycondition
  USE ccgeneraldiscretisation
  USE ccpostprocessing
  USE ccstationary
  USE ccnonstationary
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc2dmain
  
!<description>
  ! This is a 'separated' Navier-Stokes solver for solving a Navier-Stokes
  ! problem. The different tasks of the problem are separated into
  ! subroutines. The problem uses a problem-specific structure for the 
  ! communication: All subroutines add their generated information to the
  ! structure, so that the other subroutines can work with them.
  ! (This is somehow a cleaner implementation than using a collection!).
  ! For the communication to callback routines of black-box subroutines
  ! (matrix-assembly), a collection is used.
  !
  ! The following tasks are performed by the subroutines:
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Set up RHS
  ! 4.) Set up matrix
  ! 5.) Create solver structure
  ! 6.) Solve the problem
  ! 7.) Write solution to GMV file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! A problem structure for our problem
    TYPE(t_problem), POINTER :: p_rproblem
    
    ! A structure for the solution vector and the RHS vector of the problem.
    TYPE(t_vectorBlock) :: rvector,rrhs
    
    ! A structure for the postprocessing.
    TYPE(t_c2d2postprocessing) :: rpostprocessing
    
    ! Timer objects for stopping time
    TYPE(t_timer) :: rtimerTotal
    TYPE(t_timer) :: rtimerGridGeneration
    TYPE(t_timer) :: rtimerMatrixGeneration
    TYPE(t_timer) :: rtimerSolver
    
    INTEGER :: i
    
    ! Ok, let's start. 
    
    ! Initialise the timers by zero:
    CALL stat_clearTimer(rtimerTotal)
    CALL stat_clearTimer(rtimerGridGeneration)
    CALL stat_clearTimer(rtimerMatrixGeneration)
    CALL stat_clearTimer(rtimerSolver)

    ! Start the timer
    CALL stat_startTimer(rtimerTotal)
    
    ! Allocate memory fo rthe problem structure -- it's rather large!
    ALLOCATE (p_rproblem)
    
    ! Initialise the collection
    CALL collct_init (p_rproblem%rcollection)
    
    ! Initialise the parameter list object. This creates an empty parameter list.
    CALL parlst_init (p_rproblem%rparamList)
    
    ! Read parameters from the INI/DAT files into the parameter list. 
    CALL cc2d_getDAT (p_rproblem%rparamList)
    
    ! Ok, parameters are read in.
    ! Get the output levels during the initialisation phase and during the program.
    CALL cc_initOutput (p_rproblem)
    
    ! Print the configuration to the terminal
    IF (p_rproblem%MSHOW_Initialisation .GE. 2) THEN
      CALL output_line ('Parameters:')
      CALL parlst_info (p_rproblem%rparamList)
    END IF
    
    ! Evaluate these parameters and initialise global data in the problem
    ! structure for global access.
    CALL cc_initParameters (p_rproblem)
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    !
    ! Parametrisation & Triangulation
    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Initialising parametrisation / triangulation...')
    END IF
    
    CALL stat_startTimer(rtimerGridGeneration)
    
    CALL cc_initParamTriang (p_rproblem)
    
    CALL stat_stopTimer(rtimerGridGeneration)
    CALL output_lbrk ()
    CALL output_line ("Time for mesh generation: "//&
      TRIM(sys_sdL(rtimerGridGeneration%delapsedReal,10)))
    
    ! Print mesh information
    IF (p_rproblem%MSHOW_Initialisation .GE. 2) THEN
      CALL output_lbrk ()
      CALL output_line ('Mesh statistics:')
      CALL output_lbrk ()
      DO i=p_rproblem%NLMIN,p_rproblem%NLMAX
        CALL tria_infoStatistics (p_rproblem%RlevelInfo(i)%rtriangulation,&
            i .EQ. p_rproblem%NLMIN,i)
      END DO
    END IF
    
    ! Discretisation
    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Initialising discretisation...')
    END IF
    CALL cc_initDiscretisation (p_rproblem)    

    IF (p_rproblem%MSHOW_Initialisation .GE. 2) THEN
      CALL output_lbrk ()
      CALL output_line ('Discretisation statistics:')
      DO i=p_rproblem%NLMIN,p_rproblem%NLMAX
        CALL output_lbrk ()
        CALL output_line ('Level '//sys_siL(i,5))
        CALL dof_infoDiscrBlock (p_rproblem%RlevelInfo(i)%rdiscretisation,.FALSE.)
      END DO
    END IF
    
    ! And all the other stuff...
    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Initialising postprocessing...')
    END IF
    CALL cc_initPostprocessing (p_rproblem,rpostprocessing)
    
    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Initialising matrices/vectors...')
    END IF
    
    CALL stat_startTimer(rtimerMatrixGeneration)
    
    CALL cc_allocMatVec (p_rproblem,rvector,rrhs)    
    
    CALL stat_stopTimer(rtimerMatrixGeneration)
    CALL output_lbrk ()
    CALL output_line ("Time for matrix initialisation: "//&
      TRIM(sys_sdL(rtimerGridGeneration%delapsedReal,10)))

    ! On all levels, generate the static matrices used as templates
    ! for the system matrix (Laplace, B, Mass,...)
    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Generating basic matrices...')
    END IF
    CALL cc_generateBasicMatrices (p_rproblem)

    ! Create the solution vector -- zero or read from file.
    IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
      CALL output_separator (OU_SEP_MINUS)
      CALL output_line('Initialising initial solution vector...')
    END IF
    CALL cc_initInitialSolution (p_rproblem,rvector)

    ! Now choose the algorithm. Stationary or time-dependent simulation?
    IF (p_rproblem%itimedependence .EQ. 0) THEN
    
      ! Stationary simulation
      !
      ! Generate the RHS vector.
      IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
        CALL output_separator (OU_SEP_MINUS)
        CALL output_line('Generating RHS vector...')
      END IF
      CALL cc_generateBasicRHS (p_rproblem,rrhs)
      
      ! Generate discrete boundary conditions
      IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
        CALL output_separator (OU_SEP_MINUS)
        CALL output_line('Generating discrete boundary conditions...')
      END IF
      CALL cc_initDiscreteBC (p_rproblem,rvector,rrhs)

      ! Implementation of boundary conditions
      CALL cc_implementBC (p_rproblem,rvector,rrhs,.TRUE.,.TRUE.)
    
      ! Solve the problem
      IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
        CALL output_separator (OU_SEP_MINUS)
        CALL output_line('Invoking stationary solver...')
        CALL output_separator (OU_SEP_MINUS)
      END IF
      
      CALL stat_startTimer(rtimerSolver)
      
      CALL cc_solveStationary (p_rproblem,rvector,rrhs)
      
      CALL stat_stopTimer(rtimerSolver)
    
      ! Postprocessing
      CALL cc_postprocessingStationary (p_rproblem,rvector,rpostprocessing)
      
    ELSE
    
      ! Time dependent simulation with explicit time stepping.
      !
      ! Generate the RHS vector for the first time step.
      IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
        CALL output_separator (OU_SEP_MINUS)
        CALL output_line('Generating RHS vector...')
      END IF
      CALL cc_generateBasicRHS (p_rproblem,rrhs)
      
      ! Initialise the boundary conditions, but 
      ! don't implement any boundary conditions as the nonstationary solver
      ! doesn't like this.
      IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
        CALL output_separator (OU_SEP_MINUS)
        CALL output_line('Generating discrete boundary conditions of first time step...')
      END IF
      CALL cc_initDiscreteBC (p_rproblem,rvector,rrhs)
      
      ! Call the nonstationary solver to solve the problem.
      IF (p_rproblem%MSHOW_Initialisation .GE. 1) THEN
        CALL output_separator (OU_SEP_MINUS)
        CALL output_line('Invoking nonstationary solver...')
        CALL output_separator (OU_SEP_MINUS)
      END IF
      
      CALL stat_startTimer(rtimerSolver)
      
      CALL cc_solveNonstationary (p_rproblem,rvector,rrhs,rpostprocessing)
      
      CALL stat_stopTimer(rtimerSolver)
      
    END IF
    
    ! (Probably) write final solution vector
    CALL cc_writeSolution (p_rproblem,rvector)
    
    ! Cleanup
    CALL cc_doneMatVec (p_rproblem,rvector,rrhs)
    CALL cc_doneBC (p_rproblem)
    CALL cc_doneDiscretisation (p_rproblem)
    CALL cc_donepostprocessing (rpostprocessing)    
    CALL cc_doneParamTriang (p_rproblem)
    
    ! Release parameters from the DAT/INI files from the problem structure.
    CALL cc_doneParameters (p_rproblem)

    ! Release the parameter list
    CALL parlst_done (p_rproblem%rparamList)
    
    ! Print some statistical data about the collection - anything forgotten?
    CALL output_lbrk ()
    CALL output_line ('Remaining collection statistics:')
    CALL output_line ('--------------------------------')
    CALL output_lbrk ()
    CALL collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection and the problem structure.
    CALL collct_done (p_rproblem%rcollection)
    
    DEALLOCATE(p_rproblem)
    
    ! Stop the timer
    CALL stat_stopTimer(rtimerTotal)
    
    ! Print the time for the total computation
    CALL output_lbrk ()
    CALL output_line ("Time for initial mesh generation:       "//&
        TRIM(sys_sdL(rtimerGridGeneration%delapsedReal,10)))
    CALL output_line ("Time for initial matrix initialisation: "//&
      TRIM(sys_sdL(rtimerGridGeneration%delapsedReal,10)))
    CALL output_line ("Total time for solver:                  "//&
      TRIM(sys_sdL(rtimerSolver%delapsedReal,10)))
    CALL output_line ("Total time:                             "//&
        TRIM(sys_sdL(rtimerTotal%delapsedReal,10)))
    
  END SUBROUTINE

END MODULE

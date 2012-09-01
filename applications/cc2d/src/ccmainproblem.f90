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

module ccmainproblem

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
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
  use statistics
  use dofmapping
  
  use collection
  use convection
    
  use ccbasic
  use ccinitgeneralparameters
  use ccinitparamtriang
  use ccgeneraldiscretisation
  use ccpostprocessing
  use ccstationary
  use ccnonstationary
  use ccboundarycondition
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc2dmain
  
!<description>
  ! This is a "separated" Navier-Stokes solver for solving a Navier-Stokes
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
    type(t_problem), pointer :: p_rproblem
    
    ! A structure for the solution vector and the RHS vector of the problem.
    type(t_vectorBlock) :: rvector,rrhs
    
    ! A structure for the postprocessing.
    type(t_c2d2postprocessing) :: rpostprocessing
    
    ! Timer objects for stopping time
    type(t_timer) :: rtimerTotal
    type(t_timer) :: rtimerGridGeneration
    type(t_timer) :: rtimerMatrixGeneration,rtimerRHSgeneration
    type(t_timer) :: rtimerSolver
    
    integer :: i
    
    ! Ok, let us start.
    
    ! Initialise the timers by zero:
    call stat_clearTimer(rtimerTotal)
    call stat_clearTimer(rtimerGridGeneration)
    call stat_clearTimer(rtimerMatrixGeneration)
    call stat_clearTimer(rtimerSolver)

    ! Start the timer
    call stat_startTimer(rtimerTotal)
    
    ! Allocate memory for the problem structure -- it is rather large!
    allocate (p_rproblem)
    
    ! Initialise the collection
    call collct_init (p_rproblem%rcollection)
    
    ! Initialise the parameter list object. This creates an empty parameter list.
    call parlst_init (p_rproblem%rparamList)
    
    ! Read parameters from the INI/DAT files into the parameter list.
    call cc2d_getDAT (p_rproblem%rparamList)
    
    ! Ok, parameters are read in.
    ! Get the output levels during the initialisation phase and during the program.
    call cc_initOutput (p_rproblem)
    
    ! Print the configuration to the terminal
    if (p_rproblem%MSHOW_Initialisation .ge. 2) then
      call output_line ("Parameters:")
      call output_lbrk ()
      call parlst_info (p_rproblem%rparamList)
    end if
    
    ! Evaluate these parameters and initialise global data in the problem
    ! structure for global access.
    call cc_initParameters (p_rproblem)
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    !
    ! Parametrisation & Triangulation
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line("Initialising parametrisation / triangulation...")
    end if
    
    call stat_startTimer(rtimerGridGeneration)
    
    call cc_initParamTriang (p_rproblem)
    
    call stat_stopTimer(rtimerGridGeneration)
    call output_lbrk ()
    call output_line ("Time for mesh generation: "//&
      trim(sys_sdL(rtimerGridGeneration%delapsedReal,10)))
      
    p_rproblem%rstatistics%dtimeGridGeneration = &
      p_rproblem%rstatistics%dtimeGridGeneration + rtimerGridGeneration%delapsedReal
    
    ! Print mesh information
    if (p_rproblem%MSHOW_Initialisation .ge. 2) then
      call output_lbrk ()
      call output_line ("Mesh statistics:")
      call output_lbrk ()
      do i=p_rproblem%NLMIN,p_rproblem%NLMAX
        call tria_infoStatistics (p_rproblem%RlevelInfo(i)%rtriangulation,&
            i .eq. p_rproblem%NLMIN,i)
      end do
    end if
    
    ! Discretisation
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line("Initialising discretisation...")
    end if
    call cc_initDiscretisation (p_rproblem)

    if (p_rproblem%MSHOW_Initialisation .ge. 2) then
      call output_lbrk ()
      call output_line ("Discretisation statistics:")
      do i=p_rproblem%NLMIN,p_rproblem%NLMAX
        call output_lbrk ()
        call output_line ("Level "//sys_siL(i,5))
        call dof_infoDiscrBlock (p_rproblem%RlevelInfo(i)%rdiscretisation,.false.)
      end do
    end if
    
    ! And all the other stuff...
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line("Initialising postprocessing...")
    end if
    call cc_initPostprocessing (p_rproblem,rpostprocessing)
    
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line("Initialising matrices/vectors...")
    end if
    
    call stat_startTimer(rtimerMatrixGeneration)
    
    call cc_allocMatVec (p_rproblem,rvector,rrhs)
    
    call stat_stopTimer(rtimerMatrixGeneration)
    call output_lbrk ()
    call output_line ("Time for matrix initialisation: "//&
      trim(sys_sdL(rtimerGridGeneration%delapsedReal,10)))

    ! On all levels, generate the template matrices used as templates
    ! for the system matrix (Laplace, B, Mass,...)
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line("Generating basic matrices...")
    end if
    call cc_generateBasicMat (p_rproblem)

    ! Create the solution vector -- zero or read from file.
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line("Initialising initial solution vector...")
    end if
    call cc_initInitialSolution (p_rproblem,rvector)

    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line("Initialising RHS vector...")
    end if
    call cc_initRHSAssembly (p_rproblem%rparamList,rrhs%p_rblockDiscr,p_rproblem%rrhsAssembly)

    ! Now choose the algorithm. Stationary or time-dependent simulation?
    if (p_rproblem%itimedependence .eq. 0) then
    
      ! Stationary simulation
      !
      ! Generate the RHS vector.
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line("Generating RHS vector...")
      end if
      call stat_clearTimer(rtimerRHSgeneration)
      call stat_startTimer(rtimerRHSgeneration)
      call cc_generateBasicRHS (p_rproblem,&
          p_rproblem%RlevelInfo(p_rproblem%NLMAX)%rasmTempl,&
          p_rproblem%rrhsAssembly,rrhs)
      call stat_stopTimer(rtimerRHSgeneration)
      p_rproblem%rstatistics%dtimeRHSAssembly = &
          p_rproblem%rstatistics%dtimeRHSAssembly + rtimerRHSgeneration%delapsedReal
      
      ! Generate discrete boundary conditions
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line("Generating discrete boundary conditions...")
      end if
      call cc_initDiscreteBC (p_rproblem,rvector,rrhs)

      ! Implementation of boundary conditions
      call cc_implementBC (p_rproblem,rvector,rrhs,.true.,.true.)
    
      ! Solve the problem
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line("Invoking stationary solver...")
        call output_separator (OU_SEP_MINUS)
      end if
      
      call stat_startTimer(rtimerSolver)
      
      call cc_solveStationary (p_rproblem,rvector,rrhs)
      
      call stat_stopTimer(rtimerSolver)
    
      ! Postprocessing
      call cc_postprocessingStationary (p_rproblem,rvector,rpostprocessing)
      
    else
    
      ! Time dependent simulation with explicit time stepping.
      !
      ! Generate the RHS vector for the first time step.
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line("Generating RHS vector...")
      end if
      call stat_clearTimer(rtimerRHSgeneration)
      call stat_startTimer(rtimerRHSgeneration)
      call cc_generateBasicRHS (p_rproblem,&
          p_rproblem%RlevelInfo(p_rproblem%NLMAX)%rasmTempl,&
          p_rproblem%rrhsAssembly,rrhs)
      call stat_stopTimer(rtimerRHSgeneration)
      p_rproblem%rstatistics%dtimeRHSAssembly = &
          p_rproblem%rstatistics%dtimeRHSAssembly + rtimerRHSgeneration%delapsedReal
      
      ! Initialise the boundary conditions, but
      ! do not implement any boundary conditions as the nonstationary solver
      ! does not like this.
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line("Generating discrete boundary conditions of first time step...")
      end if
      call cc_initDiscreteBC (p_rproblem,rvector,rrhs)
      
      ! Call the nonstationary solver to solve the problem.
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line("Invoking nonstationary solver...")
        call output_separator (OU_SEP_MINUS)
      end if
      
      call stat_startTimer(rtimerSolver)
      
      call cc_solveNonstationary (p_rproblem,rvector,rrhs,rpostprocessing)
      
      call stat_stopTimer(rtimerSolver)
      
    end if
    
    ! Gather statistics
    p_rproblem%rstatistics%dtimeSolver = &
      p_rproblem%rstatistics%dtimeSolver + rtimerSolver%delapsedReal
    
    ! (Probably) write (final) solution vector
    call cc_writeSolution (p_rproblem,rvector)
    
    ! Cleanup
    call cc_doneRHSAssembly (p_rproblem%rrhsAssembly)
    call cc_doneMatVec (p_rproblem,rvector,rrhs)
    call cc_doneBC (p_rproblem)
    call cc_doneDiscretisation (p_rproblem)
    call cc_donepostprocessing (rpostprocessing)
    call cc_doneParamTriang (p_rproblem)
    
    ! Release parameters from the DAT/INI files from the problem structure.
    call cc_doneParameters (p_rproblem)

    ! Release the parameter list
    call parlst_done (p_rproblem%rparamList)
    
    ! Print some statistical data about the collection - anything forgotten?
    call output_lbrk ()
    call output_line ("Remaining collection statistics:")
    call output_line ("--------------------------------")
    call output_lbrk ()
    call collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection and the problem structure.
    call collct_done (p_rproblem%rcollection)
    
    ! Stop the timer
    call stat_stopTimer(rtimerTotal)

    ! Gather statistics
    p_rproblem%rstatistics%dtimeTotal = &
      p_rproblem%rstatistics%dtimeTotal + rtimerTotal%delapsedReal
    
    ! Print the time for the total computation
    call output_lbrk ()
    call output_line ("Total time:                             "//&
        trim(sys_sdL(p_rproblem%rstatistics%dtimeTotal,10)))

    call output_line ("Time for initial mesh generation:       "//&
        trim(sys_sdL(rtimerGridGeneration%delapsedReal,10)))
        
    call output_line ("Time for initial matrix assembly:       "//&
      trim(sys_sdL(rtimerMatrixGeneration%delapsedReal,10)))
      
    call output_line ("Total time for grid generation:         "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeGridGeneration,10)))
      
    call output_line ("Total time for complete solver:         "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeSolver,10)))
      
    call output_line ("Total time for nonlinear solver:        "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeNonlinearSolver,10)))
      
    call output_line ("Total time for defect calculation:      "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeDefectCalculation,10)))

    call output_line ("Total time for optimal damping:         "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeOptimalCorrection,10)))

    call output_line ("Total time for RHS assembly:            "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeRHSAssembly,10)))
      
    call output_line ("Total time for matrix assembly:         "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeMatrixAssembly,10)))
      
    call output_line ("Total time for linear solver:           "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeLinearSolver,10)))
      
    call output_line ("Total time for factorisation:           "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimeLinearSolverFactorisation,10)))

    call output_line ("Total time for postprocessing:          "//&
      trim(sys_sdL(p_rproblem%rstatistics%dtimePostprocessing,10)))
      
    call output_line ("Total #iterations nonlinear solver:     "//&
      trim(sys_siL(p_rproblem%rstatistics%nnonlinearIterations,10)))
      
    call output_line ("Total #iterations linear solver:        "//&
      trim(sys_siL(p_rproblem%rstatistics%nlinearIterations,10)))

    call output_line ("Total number of calculated timesteps:   "//&
      trim(sys_siL(p_rproblem%rstatistics%ntimesteps,10)))

    ! That is it.
    deallocate(p_rproblem)
    
  end subroutine

end module

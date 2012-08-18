!##############################################################################
!# ****************************************************************************
!# <name> cc2dmedium_method2 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstation program how to solve a stationary
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
!# routines.
!#
!# For the nonlinearity, the nonlinear solver is invoked. The
!# defect that is setted up there is preconditioned by a linear Multigrid
!# solver with a simple-VANCA smoother/preconditioner for
!# 2D saddle point problems, Jacobi-Type. As coarse grid solver,
!# UMFPACK is used.
!# </purpose>
!##############################################################################

module cc2dmedium_method2

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use statistics
  
  use collection
  use convection
    
  use cc2dmediumm2basic
  use cc2dmediumm2init
  use cc2dmediumm2boundary
  use cc2dmediumm2discretisation
  use cc2dmediumm2postprocessing
  use cc2dmediumm2stationary
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc2dmedium2_getLogFiles (slogfile,serrorfile)
  
!<description>
  ! Temporarily reads the output DAT file to get the names of the output
  ! files.
!</description>

!<output>
  ! Name of the message log file.
  character(LEN=*), intent(OUT) :: slogfile
  
  ! Name of the error log file.
  character(LEN=*), intent(OUT) :: serrorfile
!</output>

!</subroutine>

    type(t_parlist) :: rparlist
    character(LEN=SYS_STRLEN) :: sstring

    ! Init parameter list that accepts parameters for output files
    call parlst_init (rparlist)

    ! Read parameters that configure the output
    call parlst_readfromfile (rparlist, './data/output.dat')
    
    ! Now the real initialisation of the output including log file stuff!
    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'smsgLog',sstring,'')
    read(sstring,*) slogfile

    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'serrorLog',sstring,'')
    read(sstring,*) serrorfile
    
    ! That temporary parameter list is not needed anymore.
    call parlst_done (rparlist)
    
    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc2dmedium2_getDAT (rparamList)
  
!<description>
  ! Reads in all DAT files into the parameter list rparlist
!</description>

!<inputoutput>
  ! The parameter list where the values of the DAT files should be stored.
  ! The structure must have been initialised, the parameters are just added
  ! to the list.
  type(t_parlist), intent(INOUT) :: rparamList
!</inputoutput>

!</subroutine>

    ! Each 'readfromfile' command adds the parameter of the specified file
    ! to the parameter list.
    call parlst_readfromfile (rparamList, './data/discretisation.dat')
    call parlst_readfromfile (rparamList, './data/linsol_cc2d.dat')
    call parlst_readfromfile (rparamList, './data/nonlinsol_cc2d.dat')
    call parlst_readfromfile (rparamList, './data/output.dat')
    call parlst_readfromfile (rparamList, './data/paramtriang.dat')
    call parlst_readfromfile (rparamList, './data/bdconditions.dat')
    call parlst_readfromfile (rparamList, './data/timediscr.dat')
    call parlst_readfromfile (rparamList, './data/postprocessing.dat')
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc2dmedium2
  
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
    type(t_problem), pointer :: p_rproblem
    
    ! A structure for the solution vector and the RHS vector of the problem.
    type(t_vectorBlock) :: rvector,rrhs
    
    ! A structure for the postprocessing.
    type(t_c2d2postprocessing) :: rpostprocessing
    
    ! Timer objects for stopping time
    type(t_timer) :: rtimerTotal
    type(t_timer) :: rtimerGridGeneration
    type(t_timer) :: rtimerMatrixGeneration
    type(t_timer) :: rtimerSolver
    
    integer :: i
    
    ! Ok, let's start.
    
    ! Initialise the timers by zero:
    call stat_clearTimer(rtimerTotal)
    call stat_clearTimer(rtimerGridGeneration)
    call stat_clearTimer(rtimerMatrixGeneration)
    call stat_clearTimer(rtimerSolver)

    ! Start the timer
    call stat_startTimer(rtimerTotal)
    
    ! Allocate memory fo rthe problem structure -- it's rather large!
    allocate (p_rproblem)
    
    ! Initialise the collection
    call collct_init (p_rproblem%rcollection)
    do i=1,NNLEV
      call collct_addlevel_all (p_rproblem%rcollection)
    end do
    
    ! Initialise the parameter list object. This creates an empty parameter list.
    call parlst_init (p_rproblem%rparamList)
    
    ! Add the parameter list to the collection so that the parameters
    ! from the DAT/INI files are available everywhere where we have the
    ! collection.
    call collct_setvalue_parlst(p_rproblem%rcollection,'INI',&
                                p_rproblem%rparamList,.true.)

    ! Read parameters from the INI/DAT files into the parameter list.
    call cc2dmedium2_getDAT (p_rproblem%rparamList)
    
    ! Ok, parameters are read in.
    ! Get the output levels during the initialisation phase and during the program.
    call c2d2_initOutput (p_rproblem)
    
    ! Print the configuration to the terminal
    if (p_rproblem%MSHOW_Initialisation .ge. 2) then
      call output_line ('Parameters:')
      call parlst_info (p_rproblem%rparamList)
    end if
    
    ! Evaluate these parameters and initialise global data in the problem
    ! structure for global access.
    call c2d2_initParameters (p_rproblem)
    
    ! So now the different steps - one after the other.
    !
    ! Initialisation
    !
    ! Parametrisation & Triangulation
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising parametrisation / triangulation...')
    end if
    
    call stat_startTimer(rtimerGridGeneration)
    
    call c2d2_initParamTriang (p_rproblem)
    
    call stat_stopTimer(rtimerGridGeneration)
    call output_lbrk ()
    call output_line ("Time for mesh generation: "//&
      trim(sys_sdL(rtimerGridGeneration%delapsedReal,10)))
    
    ! Print mesh information
    if (p_rproblem%MSHOW_Initialisation .ge. 2) then
      call output_lbrk ()
      call output_line ('Mesh statistics:')
      call output_lbrk ()
      do i=p_rproblem%NLMIN,p_rproblem%NLMAX
        call tria_infoStatistics (p_rproblem%RlevelInfo(i)%rtriangulation,&
            i .eq. p_rproblem%NLMIN,i)
      end do
    end if
    
    ! Discretisation
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising discretisation...')
    end if
    call c2d2_initDiscretisation (p_rproblem)

    if (p_rproblem%MSHOW_Initialisation .ge. 2) then
      call output_lbrk ()
      call output_line ('Discretisation statistics:')
      do i=p_rproblem%NLMIN,p_rproblem%NLMAX
        call output_lbrk ()
        call output_line ('Level '//sys_siL(i,5))
        call dof_infoDiscrBlock (p_rproblem%RlevelInfo(i)%p_rdiscretisation,.false.)
      end do
    end if
    
    ! And all the other stuff...
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising postprocessing...')
    end if
    call c2d2_initPostprocessing (p_rproblem,rpostprocessing)
    
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising matrices/vectors...')
    end if
    
    call stat_startTimer(rtimerMatrixGeneration)
    
    call c2d2_allocMatVec (p_rproblem,rvector,rrhs)
    
    call stat_stopTimer(rtimerMatrixGeneration)
    call output_lbrk ()
    call output_line ("Time for matrix initialisation: "//&
      trim(sys_sdL(rtimerGridGeneration%delapsedReal,10)))

    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising analytic boundary conditions...')
    end if
    call c2d2_initAnalyticBC (p_rproblem)

    ! On all levels, generate the static matrices used as templates
    ! for the system matrix (Laplace, B, Mass,...)
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Generating basic matrices...')
    end if
    call c2d2_generateBasicMatrices (p_rproblem)

    ! Create the solution vector -- zero or read from file.
    if (p_rproblem%MSHOW_Initialisation .ge. 1) then
      call output_separator (OU_SEP_MINUS)
      call output_line('Initialising initial solution vector...')
    end if
    call c2d2_initInitialSolution (p_rproblem,rvector)

    ! Now choose the algorithm. Stationary or time-dependent simulation?
    if (p_rproblem%itimedependence .eq. 0) then
    
      ! Stationary simulation
      !
      ! Generate the RHS vector.
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line('Generating RHS vector...')
      end if
      call c2d2_generateBasicRHS (p_rproblem,rrhs)
      
      ! Generate discrete boundary conditions
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line('Generating discrete boundary conditions...')
      end if
      call c2d2_initDiscreteBC (p_rproblem,rvector,rrhs)

      ! Implementation of boundary conditions
      call c2d2_implementBC (p_rproblem,rvector,rrhs,.true.,.true.)
    
      ! Solve the problem
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line('Invoking stationary solver...')
        call output_separator (OU_SEP_MINUS)
      end if
      
      call stat_startTimer(rtimerSolver)
      
      call c2d2_solve (p_rproblem,rvector,rrhs)
      
      call stat_stopTimer(rtimerSolver)
    
      ! Postprocessing
      call c2d2_postprocessingStationary (p_rproblem,rvector,rpostprocessing)
      
    else
    
      ! Time dependent simulation with explicit time stepping.
      !
      ! Generate the RHS vector for the first time step.
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line('Generating RHS vector...')
      end if
      call c2d2_generateBasicRHS (p_rproblem,rrhs)
      
      ! Initialise the boundary conditions, but
      ! don't implement any boundary conditions as the nonstationary solver
      ! doesn't like this.
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line('Generating discrete boundary conditions of first time step...')
      end if
      call c2d2_initDiscreteBC (p_rproblem,rvector,rrhs)
      
      ! Call the nonstationary solver to solve the problem.
      if (p_rproblem%MSHOW_Initialisation .ge. 1) then
        call output_separator (OU_SEP_MINUS)
        call output_line('Invoking nonstationary solver...')
        call output_separator (OU_SEP_MINUS)
      end if
      
      call stat_startTimer(rtimerSolver)
      
      call c2d2_solveNonstationary (p_rproblem,rvector,rrhs,rpostprocessing)
      
      call stat_stopTimer(rtimerSolver)
      
    end if
    
    ! (Probably) write final solution vector
    call c2d2_writeSolution (p_rproblem,rvector)
    
    ! Cleanup
    call c2d2_doneMatVec (p_rproblem,rvector,rrhs)
    call c2d2_doneBC (p_rproblem)
    call c2d2_doneDiscretisation (p_rproblem)
    call c2d2_donepostprocessing (rpostprocessing)
    call c2d2_doneParamTriang (p_rproblem)
    
    ! Release parameters from the DAT/INI files from the problem structure.
    call c2d2_doneParameters (p_rproblem)

    ! Release the parameter list
    call collct_deleteValue (p_rproblem%rcollection,'INI')
    call parlst_done (p_rproblem%rparamList)
    
    ! Print some statistical data about the collection - anything forgotten?
    call output_lbrk ()
    call output_line ('Remaining collection statistics:')
    call output_line ('--------------------------------')
    call output_lbrk ()
    call collct_printStatistics (p_rproblem%rcollection)
    
    ! Finally release the collection and the problem structure.
    call collct_done (p_rproblem%rcollection)
    
    deallocate(p_rproblem)
    
    ! Stop the timer
    call stat_stopTimer(rtimerTotal)
    
    ! Print the time for the total computation
    call output_lbrk ()
    call output_line ("Time for initial mesh generation:       "//&
        trim(sys_sdL(rtimerGridGeneration%delapsedReal,10)))
    call output_line ("Time for initial matrix initialisation: "//&
      trim(sys_sdL(rtimerGridGeneration%delapsedReal,10)))
    call output_line ("Total Time for solver:                  "//&
      trim(sys_sdL(rtimerSolver%delapsedReal,10)))
    call output_line ("Total time:                             "//&
        trim(sys_sdL(rtimerTotal%delapsedReal,10)))
    
  end subroutine

end module

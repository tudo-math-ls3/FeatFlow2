!##############################################################################
!# ****************************************************************************
!# <name> main_program </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module solves an optimal control problem for the stationary and
!# nonstationary Navier-Stokes optimal control problem
!#
!#  $$ min J(y,u) = 1/2||y-z||_{L^2} + \gamma/2||y(T)-z(T)||_{L^2} + \alphga/2||u||^2 $$
!#
!#  $$- \nu Delta(y) + y*\Nabla(y) + \Nabla p = f $$
!#  $$ \Nabla \cdot y = 0$$
!#  $$- \nu Delta(\lambda) - y*\Nabla(\lambda) + \lambda\Nabla y + \Nabla \xi = y-z $$
!#  $$ \Nabla \cdot \lambda = 0$$
!#
!#
!# on a 2D domain for a 2D function $y=(y_1,y_2)$, a pressure $p$,
!# a dual velocity $\lambda$ and a dual pressure $\xi$. $u$ is the control
!# and $z$ a desired target function field.
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
!# solver with a simple-VANKA smoother/preconditioner for
!# 2D saddle point problems, Jacobi-Type. As coarse grid solver,
!# UMFPACK is used.
!# </purpose>
!##############################################################################

module main_program

  use fsystem
  use genoutput
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
  use fparser
  use statistics
  
  use collection
  use convection
    
  use externalstorage
  use paramlist
  
  use spacetimevectors
  use analyticsolution
  
  use structurespostproc
  use structuresmain
  
  use structuresoptflow
  use initoptflow
  use newtoniteration
  
  implicit none

!<globals>

  ! Directory containing the data files.
  character(LEN=SYS_STRLEN), save :: DIR_DATA = "./data";

!</globals>
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc2doptc_getLogFiles (slogfile,serrorfile)
  
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

    ! Init parameter list that accepts parameters for output files
    call parlst_init (rparlist)

    ! Read parameters that configure the output
    call parlst_readfromfile (rparlist, trim(DIR_DATA)//'/output.dat')
    
    ! Now the real initialisation of the output including log file stuff!
    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
        'smsgLog',slogfile,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
        'serrorLog',serrorfile,"",bdequote=.true.)
    
    ! That temporary parameter list is not needed anymore.
    call parlst_done (rparlist)
    
    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc2doptc_evalParameters ()
  
!<description>
  ! Evaluates command line parameters.
  ! In the current implementation, command line parameters are passed as
  ! text file. This routine searches in the main directory for a file
  ! "cmdline.dat". If this file is found, it's opened and evaluated.
  ! Every line may contain a command line parameter in the form of
  ! a DAT file (name=value pairs).
  !
  ! Supported command line parameters:
  !   "datdirectory = [Directory, where DAT files can be found]"
!</description>

!</subroutine>

    ! local variables
    type(t_parlist) :: rparamList
    logical :: bexists

    ! Figure out if the file exists.
    inquire(file='./cmdline.dat', exist=bexists)
    
    if (bexists) then
      ! Read the file
      call parlst_init (rparamList)
      call parlst_readfromfile (rparamList, './cmdline.dat')
      
      ! Evaluate parameters
      call parlst_getvalue_string ( &
          rparamList, "","datdirectory", DIR_DATA, DIR_DATA,bdequote=.true.)

      call parlst_done (rparamList)
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine main_getDat (rparamList)
  
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

    logical :: bexists
    character(LEN=SYS_STRLEN) :: smaster
    
    ! Check if a command line parameter specifies the master.dat file.
    call sys_getcommandLineArg(1,smaster,sdefault=trim(DIR_DATA)//'/master.dat')

    ! Read the file 'master.dat'.
    ! If that does not exist, try to manually read files with parameters from a
    ! couple of files.
    inquire(file=smaster, exist=bexists)
    
    if (bexists) then
      ! Read the master file. That either one contains all parameters or
      ! contains references to subfiles with data.
      call parlst_readfromfile (rparamList, trim(smaster),trim(DIR_DATA))
    else
      call output_line("Master file not found: "//trim(smaster),&
          OU_CLASS_WARNING,ssubroutine='main_getDat')
      call output_line("Reading standard parameters.",&
          OU_CLASS_WARNING,ssubroutine='main_getDat')
    
      ! Each 'readfromfile' command adds the parameter of the specified file
      ! to the parameter list.
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/main.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/bdconditions.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/discretisation.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/flows.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/linsol.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/optcontrol.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/output.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/paramtriang.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/postprocessing.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/spacetimesolver.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/timediscr.dat')
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc2doptccalculate (rsettings,rparlist)
  
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

!<input>
  ! Basic program settings
  type(t_settings_main), intent(in) :: rsettings
  
  ! Parameter list containing the DAT file parameters
  type(t_parlist), intent(in) :: rparlist
!</input>

!</subroutine>

    ! Local variables
    !
    ! The following structure configures our solver.
    ! We allocate this dynamically since the structure may be too large
    ! for the stack.
    type(t_settings_optflow), pointer :: p_rsettingsSolver
    
    ! Global solution, temp vector
    type(t_spacetimevector) :: rsolution,rtemp,rrhsDiscrete

    ! Structure for the Newton solver
    type(t_newtonParameters) :: rsolver

    ! Postprocessing data.
    type(t_optcPostprocessing) :: rpostproc
    
    ! Some timers
    type(t_timer) :: rtotalTime,rinitTime,rsolverTime,rtimePostProc,rstartvectime
    type(t_timer) :: rrhsvectime

    ! Ok, let us start.
    !
    ! Initialise the external storage management.
    
    call exstor_init (999,100)
    !CALL exstor_attachDirectory('./ff2storage')
    
    ! Ok, parameters are read in.
    ! Print the parameters to the terminal.
    if (rsettings%routput%ioutputInit .ge. 2) then
      call parlst_info (rparlist)
      call output_separator (OU_SEP_EQUAL)
    end if
    
    ! Measure the total time.
    call stat_clearTimer (rtotalTime)
    call stat_startTimer (rtotalTime)
    
!    ! Initialise the settings of the solver,
!    ! allocate all template matrices etc.
    allocate(p_rsettingsSolver)
!    allocate(p_rnlstsolver)
    
    call stat_clearTimer (rinitTime)
    call stat_startTimer (rinitTime)

    if (rsettings%routput%ioutputInit .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising solver structures.")
    end if
    call init_initOptFlow (rparlist,rsettings,p_rsettingsSolver,rsettings%routput%ioutputInit)
    
    if (rsettings%routput%ioutputInit .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising the descent algorithm loop.")
    end if
    call newtonit_initParams (rsolver,p_rsettingsSolver,rsettings%ssectionSpaceTimeSolver,rparlist)
    call newtonit_initStructure (rsolver)
    call newtonit_initData (rsolver)
    
!    ! Create a temp vector
!    call sptivec_initVector (rtemp,&
!        p_rsettingsSolver%rtimeHierarchy%p_rtimeLevels(p_rsettingsSolver%rtimeHierarchy%nlevels),&
!        p_rsettingsSolver%rfeHierPrimalDual% &
!        p_rfeSpaces(p_rsettingsSolver%rfeHierPrimalDual%nlevels)%p_rdiscretisation)
!
    call stat_stopTimer (rinitTime)

    call output_separator (OU_SEP_EQUAL)
    call output_line ("Time for initialisation            = "//&
        sys_sdL(rinitTime%delapsedReal,10))
    call output_separator (OU_SEP_EQUAL)
!        
!    ! Solve the system
!    call stat_clearTimer (rsolverTime)
!    call stat_startTimer (rsolverTime)
!    call nlstslv_solve (p_rsettingsSolver,p_rnlstsolver,rpostproc,rsolution,rrhsdiscrete,rtemp)
!    call stat_stopTimer (rsolverTime)
!    
!    call output_separator (OU_SEP_EQUAL)
!
!    ! Pipe the solution through our postprocessing routines
!    call output_line ("Postprocessing of the final solution...")
!    call stat_clearTimer (rtimePostProc)
!    call stat_startTimer (rtimePostProc)
!    call optcpp_postprocessSpaceTimeVec (rpostproc,rsolution,rrhsDiscrete,&
!        p_rsettingsSolver%rsettingsOptControl,p_rsettingsSolver)
!    call stat_stopTimer (rtimePostProc)
!    
!    ! Sum up the time for the postprocesing during the simulation
!    call stat_addTimers (p_rnlstsolver%rtimePostprocessing,rtimePostProc)
!
    call output_separator (OU_SEP_EQUAL)
!    
!    ! Print out statistics about our solver.
!    call stnlsinit_printSolverStatistics (p_rnlstsolver)
!    
!    ! Release all data
    call newtonit_doneData (rsolver)
    call newtonit_doneStructure (rsolver)
    call newtonit_done (rsolver)

    call init_doneOptFlow (p_rsettingsSolver)
    
    call stat_stopTimer (rtotalTime)

    call output_separator (OU_SEP_EQUAL)
    call output_line ("Time for initialisation            = "//&
        sys_sdL(rinitTime%delapsedReal,10))
    call output_line ("Time for creation of the start vec = "//&
        sys_sdL(rstartvectime%delapsedReal,10))
    call output_line ("Time for postprocessing            = "//&
        sys_sdL(rtimePostProc%delapsedReal,10))
    call output_line ("Time for solving                   = "//&
        sys_sdL(rsolverTime%delapsedReal,10))
    call output_line ("Total time                         = "//&
        sys_sdL(rtotalTime%delapsedReal,10))
    call output_separator (OU_SEP_EQUAL)
    
    ! Information about external storage usage
    call output_lbrk ()
    call exstor_info (bprintHandles=.true.)
    
    ! Clean up the external storage management
    call exstor_done ()
    
  end subroutine

  ! ***************************************************************************

  subroutine cc2doptcmain
    
    ! Program parameters
    type(t_parlist) :: rparlist
    type(t_settings_main) :: rsettings
    
    ! The very first thing in every application:
    ! Initialise system-wide settings:
    call system_init()
    
    ! Read the program parameters.
    call parlst_init (rparlist)
    call main_getDat (rparlist)
    
    ! Read the basic parameter settings from the "MAIN" section.
    call smain_initMainParams (rsettings,rparlist,"MAIN")
    
    ! Initialise log file for output.
    call output_init (rsettings%routput%smsgLog,rsettings%routput%serrorLog)
    OU_LINE_LENGTH = 132
    cdefaultDateTimeLogPolicy = OU_DTP_ADDDATETIME
    cdatetimeLogFormat = 1
    
    ! Now we can really start!
    !
    ! Initialise the storage management:
    call storage_init(999, 100)
    
    ! Initialise the parser
    call fparser_init ()
    
    ! Call the problem to solve.
    call output_lbrk ()
    call output_line ('Calculating cc2doptc-Problem')
    call output_separator (OU_SEP_MINUS)
    
    call cc2doptccalculate (rsettings,rparlist)

    ! Release the parser
    call fparser_done ()

    ! Print out heap statistics - just to check if everything
    ! is cleaned up.
    ! This should display 'Handles in use=0' and 'Memory in use=0'!
    call output_lbrk ()
    call storage_info(.true.)
    
    ! Clean up the storage management, parameter list, finish
    call storage_done()
    call parlst_done (rparlist)
    
  end subroutine

end module

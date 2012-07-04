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
  use kktsystem
  use kktsystemhierarchy
  use initoptflow
  use newtoniteration
  
  use postprocessing
  
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

  subroutine main_printSolverStatistics (rnewtonStatistics,&
      dtimeInit,dtimeStartVec,dtimePostproc,dtimeTotal)
  
!<description>
  ! Prints out statistics about the solution process.
!</description>

!<input>
  ! Statistics of the Newton solver
  type(t_newtonitSolverStat), intent(in) :: rnewtonStatistics
  
  ! Time for initialisation
  real(DP), intent(in) :: dtimeInit
  
  ! Time for creation of the start vector
  real(DP), intent(in) :: dtimeStartVec

  ! Time for final postprocessing
  real(DP), intent(in) :: dtimePostproc

  ! Total computation time
  real(DP), intent(in) :: dtimeTotal
!</input>

!</subroutine>

    ! Global time
    call output_line ("Total time                                   = "//&
        sys_sdL(dtimeTotal,10))
    call output_lbrk()

    call output_line ("Time for initialisation                      = "//&
        sys_sdL(dtimeInit,10))
    call output_line ("Time for creation of the start vector        = "//&
        sys_sdL(dtimeStartVec,10))
    call output_line ("Time for postprocessing                      = "//&
        sys_sdL(dtimePostproc+rnewtonStatistics%rtimePostprocessing%delapsedReal,10))
    call output_line ("Time for solving                             = "//&
        sys_sdL(rnewtonStatistics%rtotalTime%delapsedReal &
               -rnewtonStatistics%rtimePostprocessing%delapsedReal,10))
    call output_lbrk()

    ! More detailed statistics of the nonlinear solver
    call output_line ("  Nonlinear iterations                       = "//&
        sys_siL(rnewtonStatistics%niterations,1))
    call output_line ("  Time for defect calculations               = "//&
        sys_sdL(rnewtonStatistics%rtimeDefect%delapsedReal,10))
    call output_line ("    Time for forward simulations             = "//&
        sys_sdL(rnewtonStatistics%rtimeForward%delapsedReal,10))
    call output_line ("    Time for backward simulations            = "//&
        sys_sdL(rnewtonStatistics%rtimeBackward%delapsedReal,10))
    call output_line ("    #Iterations nonlinear solver in space    = "//&
        sys_siL(rnewtonStatistics%rspaceslSolverStat%niterations,10))
    call output_line ("    #Iterations linear solver in space       = "//&
        sys_siL(rnewtonStatistics%rspaceslSolverStat%rlssSolverStat%niterations,10))
    call output_line ("  Time for interpolation to lower levels     = "//&
        sys_sdL(rnewtonStatistics%rtimeProlRest%delapsedReal,10))
    call output_line ("  Time linear solver in space                = "//&
        sys_sdL(rnewtonStatistics%rspaceslSolverStat%rtotalTime%delapsedReal,10))
    call output_line ("  #Iterations linear solver in space         = "//&
        sys_siL(rnewtonStatistics%rspaceslSolverStat%niterations,10))
    call output_lbrk()

    ! More detailed statistics about the linear space-time solver    
    call output_line ("  #Iterations linear solver in space         = "//&
        sys_siL(rnewtonStatistics%rnewtonlinSolverStat%niterations,10))
    call output_line ("  Time for linear space-time solver          = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtotalTime%delapsedReal,10))

    call output_line ("    Time for smoothing                       = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeSmoothing%delapsedReal,10))
    call output_line ("      Time for smoothing (fine)              = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeSmoothingFine%delapsedReal,10))
    call output_line ("    Time for prolongation/restriction        = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeProlRest%delapsedReal,10))
    call output_line ("    Time for coarse grid solving             = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeSolverCoarse%delapsedReal,10))
    call output_lbrk()

    call output_line ("    Time for defect calculations             = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeDefect%delapsedReal,10))
    call output_line ("      Time for forward simulations           = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeForwardLin%delapsedReal,10))
    call output_line ("      Time for backward simulations          = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeBackwardLin%delapsedReal,10))
    call output_line ("      #Iterations linear solver in space     = "//&
        sys_siL(rnewtonStatistics%rnewtonlinSolverStat%rspaceslSolverStat%&
                rlssSolverStat%niterations,10))
    call output_line ("      Time linear solver in space            = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rspaceslSolverStat%&
                rtotalTime%delapsedReal,10))
    call output_lbrk()

    call output_line ("    Time for defect calculations (fine)      = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeDefectFine%delapsedReal,10))
    call output_line ("      Time for forward simulations (fine)    = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeForwardLinFine%delapsedReal,10))
    call output_line ("      Time for backward simulations (fine)   = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeBackwardLinFine%delapsedReal,10))
    call output_line ("      #Iterations linear solver (sp., fine)  = "//&
        sys_siL(rnewtonStatistics%rnewtonlinSolverStat%rspaceslSolverStatFine%&
                rlssSolverStat%niterations,10))
    call output_line ("      Time linear solver in space (fine)     = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rspaceslSolverStatFine%&
                rtotalTime%delapsedReal,10))
    call output_lbrk()

    call output_line ("    Time for defect calculations (coarse)    = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeDefectCoarse%delapsedReal,10))
    call output_line ("      Time for forward simulations (coarse)  = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeForwardLinCoarse%delapsedReal,10))
    call output_line ("      Time for backward simulations (coarse) = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rtimeBackwardLinCoarse%delapsedReal,10))
    call output_line ("      #Iterations linear solver (sp, coarse) = "//&
        sys_siL(rnewtonStatistics%rnewtonlinSolverStat%rspaceslSolverStatCoarse%&
                rlssSolverStat%niterations,10))
    call output_line ("      Time linear solver in space (coarse)   = "//&
        sys_sdL(rnewtonStatistics%rnewtonlinSolverStat%rspaceslSolverStatCoarse%rtotalTime%delapsedReal,10))

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
    
    ! Structure for the Newton solver
    type(t_spacetimeNewton) :: rsolver

    ! Structure encapsuling a hierarchy of solutions
    type(t_kktsystemHierarchy) :: rkktsystemHierarchy
    
    ! Structure holding the solution on the maximum level
    type(t_kktSystem) :: rsolution
    
    ! Some timers
    type(t_timer) :: rtotalTime,rinitTime,rsolverTime,rtimePostProc,rstartvectime
    type(t_newtonitSolverStat) :: rnewtonStatistics

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
    
    ! Create a solution
    call kkth_initHierarchy (rkktsystemHierarchy,&
        p_rsettingsSolver%roperatorAsmHier,&
        p_rsettingsSolver%rspaceTimeHierPrimal,&
        p_rsettingsSolver%rspaceTimeHierDual,&
        p_rsettingsSolver%rspaceTimeHierControl,.false.)
        
    ! Allocate a solution
    call kkth_initKKTSystem (rsolution,rkktsystemHierarchy,&
        rkktsystemHierarchy%nlevels,p_rsettingsSolver%roperatorAsmHier)
    
    if (rsettings%routput%ioutputInit .ge. 1) then
      call output_lbrk()
      call output_line ("Initialising the descent algorithm loop.")
    end if
    call newtonit_init (rsolver,p_rsettingsSolver,rsettings%ssectionSpaceTimeSolver,rparlist)
    call newtonit_initStructure (rsolver,rkktsystemHierarchy,&
        p_rsettingsSolver%rprjHierSpaceTimePrimal,p_rsettingsSolver%rprjHierSpaceTimeDual,&
        p_rsettingsSolver%rprjHierSpaceTimeControl)
    
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
        
    ! Solve the system
    call stat_clearTimer (rsolverTime)
    call stat_startTimer (rsolverTime)
    call newtonit_solve (rsolver,rsolution,rnewtonStatistics)
    call stat_stopTimer (rsolverTime)
    
    call output_separator (OU_SEP_EQUAL)

    ! Pipe the solution through our postprocessing routines
    call output_line ("Postprocessing of the final solution...")

    call stat_clearTimer (rtimePostProc)
    call stat_startTimer (rtimePostProc)
    call optcpp_postprocessing (p_rsettingsSolver%rpostproc,&
        p_rsettingsSolver%rspaceTimeHierPrimal%p_rfeHierarchy%nlevels,&
        p_rsettingsSolver%rspaceTimeHierPrimal%p_rtimeHierarchy%nlevels,&
        rsolution%p_rprimalSol,rsolution%p_rdualSol,rsolution%p_rcontrol,&
        p_rsettingsSolver)
!    call optcpp_postprocessSpaceTimeVec (rpostproc,rsolution,rrhsDiscrete,&
!        p_rsettingsSolver%rsettingsOptControl,p_rsettingsSolver)
    call stat_stopTimer (rtimePostProc)
    
    call output_separator (OU_SEP_EQUAL)

!    ! Release all data
    call newtonit_doneStructure (rsolver)
    call newtonit_done (rsolver)
    
    ! Release the solution
    call kkth_doneHierarchy (rkktSystemHierarchy)
    call kkt_doneKKTsystem (rsolution)

    ! Release solver structures
    call init_doneOptFlow (p_rsettingsSolver)
    
    call stat_stopTimer (rtotalTime)

    call output_separator (OU_SEP_EQUAL)
    call main_printSolverStatistics (&
        rnewtonStatistics,rinitTime%delapsedReal,rstartvectime%delapsedReal,&
        rtimePostProc%delapsedReal,rtotalTime%delapsedReal)
    
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

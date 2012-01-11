!##############################################################################
!# ****************************************************************************
!# <name> structuresmain </name>
!# ****************************************************************************
!#
!# <purpose>
!# Structures of the main application
!# </purpose>
!##############################################################################

module structuresmain

  use fsystem
  use storage
  use paramlist
  
  implicit none
  
  private
  
  public :: t_settings_genoutput
  public :: t_settings_main
  public :: smain_initMainParams
  
!<types>

!<typeblock>
  
  ! Structure encapsuling output settings
  type t_settings_genoutput
  
    ! Level of output during the initialisation phase
    ! =0: no output,
    ! =1: basic output
    ! =2: standard output
    integer :: ioutputInit = 2

    ! Log file for messages.
    ! ='': Disable log file.
    character(len=SYS_STRLEN) :: smsgLog = "log/output.log"

    ! Log file for error messages; usually coincides with smsgLog to print
    ! errors into the same log file as standard messages.
    ! ='': Use log file.
    character(len=SYS_STRLEN) :: serrorLog = ""
    
  end type

!</typeblock>

!<typeblock>

  ! Structure of the main program
  type t_settings_main
  
    ! Section defining the output parameters
    character(len=SYS_STRLEN) :: ssectionOutput = "GENERALOUTPUT";

    ! Section defining Debug flags
    character(len=SYS_STRLEN) :: ssectionDebugFlags = "DEBUG";

    ! Section defining the underlying mesh
    character(len=SYS_STRLEN) :: ssectionParamTriang = "PARAMTRIANG";

    ! Name of the section defining the basic space discretisation
    ! (refinement, stabilisation,...)
    character(len=SYS_STRLEN) :: ssectionDiscrSpace = "CC-DISCRETISATION";

    ! Section defining prolongation/restriction in space
    character(len=SYS_STRLEN) :: ssectionProlRestSpace = "CC-PROLREST";

    ! Section defining the basic time discretisation
    character(len=SYS_STRLEN) :: ssectionTimeDiscretisation = "TIME-DISCRETISATION";

    ! Section defining the space-time refinement
    character(len=SYS_STRLEN) :: ssectionSpaceTimeRefinement = "SPACETIME-REFINEMENT";

    ! Section defining the optimal control problem.
    character(len=SYS_STRLEN) :: ssectionOptControl = "OPTIMALCONTROL";

    ! Section defining the nonlinear space-time solver
    character(len=SYS_STRLEN) :: ssectionSpaceTimeSolver = "TIME-SOLVER";

    ! Section defining the space-time preprocessing
    character(len=SYS_STRLEN) :: ssectionSpaceTimePreprocessing = "TIME-PREPROCESSING";

    ! Section defining the space-time postprocessing
    character(len=SYS_STRLEN) :: ssectionSpaceTimePostprocessing = "TIME-POSTPROCESSING";
  
    ! Structure encapsuling output settings
    type(t_settings_genoutput) :: routput
  
  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine smain_initOutputParams (routput,rparlist,ssection)
  
!<description>
  ! Reads basic output settings based on the parameters in the DAT file.
!</description>
  
!<input>
  ! Parameter list containing the output parameters.
  type(t_parlist), intent(in) :: rparlist
  
  ! Section containing the parameters for output
  character(len=*), intent(in) :: ssection
!</input>

!<inputoutput>
  ! Output structure to initialise
  type(t_settings_genoutput), intent(out) :: routput
!</inputoutput>

!</subroutine>

    ! Get the output level for the whole application -- during the
    ! initialisation phase and during the rest of the program.
    call parlst_getvalue_int (rparlist,ssection,&
        'ioutputInit',routput%ioutputInit,2)

    call parlst_getvalue_string (rparlist,ssection,&
        'smsgLog',routput%smsgLog,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssection,&
        'serrorLog',routput%serrorlog,"",bdequote=.true.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine smain_initMainParams (rsettings,rparlist,ssection)
  
!<description>
  ! Reads basic program settings from on the parameters in the DAT file.
!</description>
  
!<input>
  ! Parameter list containing the parameters.
  type(t_parlist), intent(in) :: rparlist
  
  ! Section containing the parameters
  character(len=*), intent(in) :: ssection
!</input>

!<inputoutput>
  ! Structure to initialise
  type(t_settings_main), intent(out) :: rsettings
!</inputoutput>

!</subroutine>

    ! Get parameters
    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionOutput",rsettings%ssectionOutput,&
        rsettings%ssectionOutput,bdequote=.true.)
    
    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionDebugFlags",rsettings%ssectionDebugFlags,&
        rsettings%ssectionDebugFlags,bdequote=.true.)
    
    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionParamTriang",rsettings%ssectionParamTriang,&
        rsettings%ssectionParamTriang,bdequote=.true.)
    
    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionDiscrSpace",rsettings%ssectionDiscrSpace,&
        rsettings%ssectionDiscrSpace,bdequote=.true.)
    
    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionProlRestSpace",rsettings%ssectionProlRestSpace,&
        rsettings%ssectionProlRestSpace,bdequote=.true.)
    
    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionTimeDiscretisation",rsettings%ssectionTimeDiscretisation,&
        rsettings%ssectionTimeDiscretisation,bdequote=.true.)
    
    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionSpaceTimeRefinement",rsettings%ssectionSpaceTimeRefinement,&
        rsettings%ssectionSpaceTimeRefinement,bdequote=.true.)
    
    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionOptControl",rsettings%ssectionOptControl,&
        rsettings%ssectionOptControl,bdequote=.true.)
    
    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionSpaceTimeSolver",rsettings%ssectionSpaceTimeSolver,&
        rsettings%ssectionSpaceTimeSolver,bdequote=.true.)
    
    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionSpaceTimePostprocessing",rsettings%ssectionSpaceTimePostprocessing,&
        rsettings%ssectionSpaceTimePreprocessing,bdequote=.true.)
    
    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionSpaceTimePostprocessing",rsettings%ssectionSpaceTimePostprocessing,&
        rsettings%ssectionSpaceTimePostprocessing,bdequote=.true.)

    ! Read output settings
    call smain_initOutputParams (rsettings%routput,rparlist,rsettings%ssectionOutput)

  end subroutine


end module

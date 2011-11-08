!##############################################################################
!# ****************************************************************************
!# <name> ccinitgeneralparameters </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains basic initialisation routines for CC2D:
!# Initialisation of the main structures with parameters:
!#
!# 1.) cc_getLogFiles
!#     -> Get information for LOG files from DAT files.
!#
!# 2.) cc2d_getDAT
!#     -> Read all DAT files.
!#
!# 3.) cc_initParameters
!#     -> Init the problem structure with data from the INI/DAT files
!#
!# 4.) cc_doneParameters
!#     -> Clean up the problem structure
!#
!# </purpose>
!##############################################################################

module ccinitgeneralparameters

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
  
  use collection
  use convection
    
  use ccbasic
  use ccnonstationary
  
  implicit none
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine cc_getLogFiles (slogfile,serrorfile,sbenchlogfile)
  
!<description>
  ! Temporarily reads the output DAT file to get the names of the output
  ! files.
!</description>

!<output>
  ! Name of the message log file.
  character(LEN=*), intent(OUT) :: slogfile
  
  ! Name of the error log file.
  character(LEN=*), intent(OUT) :: serrorfile

  ! Name of the benchmark log file.
  character(LEN=*), intent(OUT) :: sbenchlogfile
!</output>

!</subroutine>

    type(t_parlist) :: rparlist
    character(LEN=SYS_STRLEN) :: sstring,smaster
    logical :: bexists

    ! Init parameter list that accepts parameters for output files
    call parlst_init (rparlist)
    
    ! Check if a command line parameter specifies the master.dat file.
    call sys_getcommandLineArg(1,smaster,sdefault='./data/master.dat')

    ! Read parameters that configure the output
    inquire(file=smaster, exist=bexists)
    
    if (bexists) then
      ! Read the master file. That either one contains all parameters or
      ! contains references to subfiles with data.
      call parlst_readfromfile (rparlist, smaster)
    else
      call parlst_readfromfile (rparlist, './data/output.dat')
    end if
    
    ! Now the real initialisation of the output including log file stuff!
    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'smsgLog',sstring,'''''')
    read(sstring,*) slogfile

    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'serrorLog',sstring,'''''')
    read(sstring,*) serrorfile

    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'sbenchLog',sstring,'''''')
    read(sstring,*) sbenchlogfile
    
    ! That temporary parameter list is not needed anymore.
    call parlst_done (rparlist)
    
    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc2d_getDAT (rparamList)
  
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
    call sys_getcommandLineArg(1,smaster,sdefault='./data/master.dat')

    ! Read the file 'master.dat'.
    ! If that does not exist, try to manually read files with parameters from a
    ! couple of files.
    inquire(file=smaster, exist=bexists)
    
    if (bexists) then
      ! Read the master file. That either one contains all parameters or
      ! contains references to subfiles with data.
      call parlst_readfromfile (rparamList, smaster)
    else
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
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initOutput (rproblem)
  
!<description>
  ! Initialises basic output settings based on the parameters in the DAT file.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! Get the output level for the whole application -- during the
    ! initialisation phase and during the rest of the program.
    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MSHOW_Initialisation',rproblem%MSHOW_Initialisation,2)

    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MT_OutputLevel',rproblem%MT_OutputLevel,2)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initParameters (rproblem)
  
!<description>
  ! Initialises the structure rproblem with data from the initialisation
  ! files.
  !
  ! The parameters in rproblem\%rparameters are evaluated.
  ! Important parameters are written to the problem structure
  ! rproblem and/or the enclosed collection rproblem\%rcollection.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    real(DP) :: dnu
    integer :: ilvmin,ilvmax,i1

    ! Get the output level for the whole application -- during the
    ! initialisation phase and during the rest of the program.
    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MSHOW_Initialisation',rproblem%MSHOW_Initialisation,2)

    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MT_OutputLevel',rproblem%MT_OutputLevel,2)

    ! Get the viscosity parameter, save it to the problem structure
    ! as well as into the collection.
    ! Note that the parameter in the DAT file is 1/nu !
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'RE',dnu,1000.0_DP)

    dnu = 1.0_DP/dnu
    rproblem%dnu = dnu
    
    ! Get min/max level from the parameter file.
    !
    ! ilvmin receives the minimal level where to discretise for supporting
    ! the solution process.
    ! ilvmax receives the level where we want to solve.
    
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMIN',ilvmin,2)
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMAX',ilvmax,4)

    ! Initialise the level in the problem structure
    if(ilvmin .le. 0) then
      rproblem%NLMIN = max(1,ilvmax+ilvmin)
    else
      rproblem%NLMIN = ilvmin
    end if
    
    rproblem%NLMAX = ilvmax
    
    ! Allocate memory for all the levels.
    allocate(rproblem%RlevelInfo(1:ilvmax))

    ! Which type of problem to discretise? (Stokes, Navier-Stokes,...)
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iEquation',i1,0)
    rproblem%iequation = i1

    ! Type of subproblem (gradient tensor, deformation tensor,...)
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'isubEquation',i1,0)
    rproblem%isubEquation = i1

    ! Type of boundary conditions
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iBoundary',rproblem%iboundary,0)

    ! Time dependence
    call cc_initParTimeDependence (rproblem,'TIME-DISCRETISATION',&
        rproblem%rparamList)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneParameters (rproblem)
  
!<description>
  ! Cleans up parameters read from the DAT files. Removes all references to
  ! parameters from the collection rproblem\%rcollection that were
  ! set up in cc_initParameters.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! Deallocate memory
    deallocate(rproblem%RlevelInfo)

  end subroutine

end module

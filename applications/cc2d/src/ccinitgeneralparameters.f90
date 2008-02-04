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
!# 2.) cc2dmedium2_getDAT
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

MODULE ccinitgeneralparameters

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
  
  USE collection
  USE convection
    
  USE ccbasic
  USE ccnonstationary
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_getLogFiles (slogfile,serrorfile)
  
!<description>
  ! Temporarily reads the output DAT file to get the names of the output
  ! files.
!</description>

!<output>
  ! Name of the message log file.
  CHARACTER(LEN=*), INTENT(OUT) :: slogfile
  
  ! Name of the error log file.
  CHARACTER(LEN=*), INTENT(OUT) :: serrorfile
!</output>

!</subroutine>

    TYPE(t_parlist) :: rparlist
    CHARACTER(LEN=SYS_STRLEN) :: sstring

    ! Init parameter list that accepts parameters for output files
    CALL parlst_init (rparlist)

    ! Read parameters that configure the output
    CALL parlst_readfromfile (rparlist, './data/output.dat')
    
    ! Now the real initialisation of the output including log file stuff!
    CALL parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'smsgLog',sstring,'')
    READ(sstring,*) slogfile

    CALL parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'serrorLog',sstring,'')
    READ(sstring,*) serrorfile
    
    ! That temporary parameter list is not needed anymore.
    CALL parlst_done (rparlist)
    
    END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc2dmedium2_getDAT (rparamList)
  
!<description>
  ! Reads in all DAT files into the parameter list rparlist
!</description>

!<inputoutput>
  ! The parameter list where the values of the DAT files should be stored.
  ! The structure must have been initialised, the parameters are just added
  ! to the list.
  TYPE(t_parlist), INTENT(INOUT) :: rparamList
!</inputoutput>

!</subroutine>

    ! Each 'readfromfile' command adds the parameter of the specified file 
    ! to the parameter list.
    CALL parlst_readfromfile (rparamList, './data/discretisation.dat')
    CALL parlst_readfromfile (rparamList, './data/linsol_cc2d.dat')
    CALL parlst_readfromfile (rparamList, './data/nonlinsol_cc2d.dat')
    CALL parlst_readfromfile (rparamList, './data/output.dat')
    CALL parlst_readfromfile (rparamList, './data/paramtriang.dat')
    CALL parlst_readfromfile (rparamList, './data/bdconditions.dat')
    CALL parlst_readfromfile (rparamList, './data/timediscr.dat')
    CALL parlst_readfromfile (rparamList, './data/postprocessing.dat')
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initOutput (rproblem)
  
!<description>
  ! Initialises basic output settings based on the parameters in the DAT file.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! Get the output level for the whole application -- during the
    ! initialisation phase and during the rest of the program.
    CALL parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MSHOW_Initialisation',rproblem%MSHOW_Initialisation,2)

    CALL parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MT_OutputLevel',rproblem%MT_OutputLevel,2)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_initParameters (rproblem)
  
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
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    REAL(DP) :: dnu
    INTEGER :: ilvmin,ilvmax,i1

    ! Get the output level for the whole application -- during the
    ! initialisation phase and during the rest of the program.
    CALL parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MSHOW_Initialisation',rproblem%MSHOW_Initialisation,2)

    CALL parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MT_OutputLevel',rproblem%MT_OutputLevel,2)

    ! Get the viscosity parameter, save it to the problem structure
    ! as well as into the collection.
    ! Note that the parameter in the DAT file is 1/nu !
    CALL parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'RE',dnu,1000.0_DP)

    dnu = 1E0_DP/dnu
    rproblem%dnu = dnu
    
    ! By default, X- and Y-velocity matrix are coupled.
    rproblem%bdecoupledXY = .FALSE.
    
    ! Get min/max level from the parameter file.
    !
    ! ilvmin receives the minimal level where to discretise for supporting
    ! the solution process.
    ! ilvmax receives the level where we want to solve.
    
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMIN',ilvmin,2)
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'NLMAX',ilvmax,4)

    ! Initialise the level in the problem structure
    rproblem%NLMIN = ilvmin
    rproblem%NLMAX = ilvmax
    
    ! Allocate memory for all the levels.
    ALLOCATE(rproblem%RlevelInfo(1:ilvmax))

    ! Which type of problem to discretise? (Stokes, Navier-Stokes,...)
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iEquation',i1,0)
    rproblem%iequation = i1

    ! Type of subproblem (gradient tensor, deformation tensor,...)
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'isubEquation',i1,0)
    rproblem%isubEquation = i1

    ! Type of boundary conditions
    CALL parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iBoundary',rproblem%iboundary,0)

    ! Time dependence
    CALL cc_initParTimeDependence (rproblem,'TIME-DISCRETISATION',&
        rproblem%rparamList)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE cc_doneParameters (rproblem)
  
!<description>
  ! Cleans up parameters read from the DAT files. Removes all references to
  ! parameters from the collection rproblem\%rcollection that were
  ! set up in cc_initParameters.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  TYPE(t_problem), INTENT(INOUT) :: rproblem
!</inputoutput>

!</subroutine>

    ! Deallocate memory
    DEALLOCATE(rproblem%RlevelInfo)

  END SUBROUTINE

END MODULE

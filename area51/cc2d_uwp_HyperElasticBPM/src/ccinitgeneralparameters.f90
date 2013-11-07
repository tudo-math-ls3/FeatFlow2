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
!# 5.) cc_initRhsAssembly
!#     -> Initialises an assembly structure for creating the RHS.
!#
!# 6.) cc_doneRhsAssembly
!#     -> Releases the RHS assembly structure.
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
  use vectorio
  
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
  character(LEN=*), intent(out) :: slogfile
  
  ! Name of the error log file.
  character(LEN=*), intent(out) :: serrorfile

  ! Name of the benchmark log file.
  character(LEN=*), intent(out) :: sbenchlogfile
!</output>

!</subroutine>

    type(t_parlist) :: rparlist
    character(LEN=SYS_STRLEN) :: smaster
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
                                'smsgLog',slogfile,'',bdequote=.true.)

    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'serrorLog',serrorfile,'',bdequote=.true.)

    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
                                'sbenchLog',sbenchlogfile,'',bdequote=.true.)
    
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
  type(t_parlist), intent(inout) :: rparamList
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
  type(t_problem), intent(inout) :: rproblem
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

  subroutine cc_initPhysics (rproblem)
  
!<description>
  ! Reads parameters about the physics of the problem.
!</description>
  
!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! is the convective terms in fluid balance of momentum to be ignored ?
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iEquation',rproblem%rphysics%iequation,0)

    ! Initial Solidity
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                              'nSo',rproblem%rphysics%dnSo,0.9_DP)

    ! A  parameter that determines the strength of the deformation dependency
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'kappa',rproblem%rphysics%kappa,0)
    ! Initial Porosity
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                              'nFo',rproblem%rphysics%dnFo,0.1_DP)

    ! Initial Permeability
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                              'kFo',rproblem%rphysics%dkFo,0.05_DP)

    ! Pore fluid density
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                              'rhoFR',rproblem%rphysics%drhoFR,1000.0_DP)

    ! Solid skeleton density
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                              'rhoSR',rproblem%rphysics%drhoSR,7800.0_DP)

    ! Lame constants for solid skeleton \lambda and \mu
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                              'lambda',rproblem%rphysics%dlambda,1000.0_DP)

    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                              'mu',rproblem%rphysics%dmu,1000.0_DP)

    ! Body load
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                              'b',rproblem%rphysics%b,9.81_DP)

    ! Shall we include the pore fluid shear stress effect
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                              'IncShear',rproblem%rphysics%IncShear,1.0_DP)



    ! Type of subproblem (gradient tensor, deformation tensor,...) 
    ! see ffunctionViscoModel in ccmatvecassembly.f90
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'isubEquation',rproblem%rphysics%isubEquation,0)

    ! Get the viscosity model
    ! Standard = 0 = constant viscosity
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'cviscoModel',rproblem%rphysics%cviscoModel,0)
                                 
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'dviscoexponent',rproblem%rphysics%dviscoexponent,2.0_DP)
                                 
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'dviscoEps',rproblem%rphysics%dviscoEps,0.01_DP)

    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'dviscoYield',rproblem%rphysics%dviscoYield,1.0_DP)

    ! Get the viscosity parameter, save it to the problem structure
    ! as well as into the collection.
    ! Note that the parameter in the DAT file is 1/nu !
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'RE',rproblem%rphysics%dnu,1000.0_DP)
    
    rproblem%rphysics%dnu = 1.0_DP/rproblem%rphysics%dnu*rproblem%rphysics%IncShear
    
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
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    integer :: ilvmin,ilvmax

    ! Get the output level for the whole application -- during the
    ! initialisation phase and during the rest of the program.
    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MSHOW_Initialisation',rproblem%MSHOW_Initialisation,2)

    call parlst_getvalue_int (rproblem%rparamList,'GENERALOUTPUT',&
                              'MT_OutputLevel',rproblem%MT_OutputLevel,2)

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

    ! Type of boundary conditions
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iBoundary',rproblem%iboundary,0)

    ! Time dependence
    call cc_initParTimeDependence (rproblem,'TIME-DISCRETISATION',&
        rproblem%rparamList)
        
    ! Get the physics of the problem.
    call cc_initPhysics (rproblem)

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
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! Deallocate memory
    deallocate(rproblem%RlevelInfo)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initRhsAssembly (rparlist,rdiscretisation,rrhsAssembly)
  
!<description>
  ! Initialises a RHS assembly structure based on the parameters in the DAT
  ! file.
!</description>
  
!<input>
  ! Parameter list with parameters.
  type(t_parlist), intent(in) :: rparlist
  
  ! Discretisation structure of the topmost level where the RHS exists.
  type(t_blockDiscretisation), intent(inout), target :: rdiscretisation
!</input>

!<output>
  ! RHS assembly structure, initialised by parameters in rparlist.
  type(t_rhsAssembly), intent(out) :: rrhsAssembly
!</output>

!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sarray

    ! Just get the parameters.
    call parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
        'irhs',rrhsAssembly%ctype,0)
    
    call parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
        'irhsFirstIndex',rrhsAssembly%ifirstindex,0)

    call parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
        'irhsFileCount',rrhsAssembly%inumfiles,0)

    call parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
        'irhsFormatted',rrhsAssembly%iformatted,1)

    call parlst_getvalue_double (rparlist,'CC-DISCRETISATION',&
        'drhsTimeInit',rrhsAssembly%dtimeInit,0.0_DP)

    call parlst_getvalue_double (rparlist,'CC-DISCRETISATION',&
        'drhsTimeMax',rrhsAssembly%dtimeMax,0.0_DP)

    call parlst_getvalue_double (rparlist,'CC-DISCRETISATION',&
        'drhsMultiplyX',rrhsAssembly%dmultiplyX,1.0_DP)

    call parlst_getvalue_double (rparlist,'CC-DISCRETISATION',&
        'drhsMultiplyY',rrhsAssembly%dmultiplyY,1.0_DP)
    
    call parlst_getvalue_string (rparlist,'CC-DISCRETISATION',&
        'sfilenameRHS',rrhsAssembly%sfilename,'',bdequote=.true.)

    if ((rrhsAssembly%ctype .eq. 3) .or. (rrhsAssembly%ctype .eq. 4)) then
      if (rrhsAssembly%sfilename .eq. "") then
        call output_line ("No filename for the RHS specified!", &
            OU_CLASS_ERROR,OU_MODE_STD,"cc_initRhs")
        call sys_halt()
      end if
    end if

    ! Is this a stationary solution based on a file?
    if ((rrhsAssembly%ctype .eq. 3) .or. (rrhsAssembly%ctype .eq. 4)) then
    
      ! Create a RHS vector and read it from the file.
      call lsysbl_createVectorBlock (rdiscretisation,rrhsAssembly%rrhsVector)
      
      ! Read in the RHS now?
      if (rrhsAssembly%ctype .eq. 3) then
        call output_line ("Reading RHS file """//trim(rrhsAssembly%sfilename)//""".",&
            coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
        call vecio_readBlockVectorHR (rrhsAssembly%rrhsVector, sarray, .true.,&
            0, rrhsAssembly%sfilename, rrhsAssembly%iformatted .ne. 0)
      else
        ! Create also the 2nd vector. Both vectors are used as temp vectors
        ! when reading in the RHS during the simulation.
        call lsysbl_createVectorBlock (rdiscretisation,rrhsAssembly%rrhsVector2)
      end if
        
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneRhsAssembly (rrhsAssembly)
  
!<description>
  ! Cleans up a RHS structure.
!</description>
  
!<inputoutput>
  ! RHS assembly structure to be cleaned up.
  type(t_rhsAssembly), intent(inout) :: rrhsAssembly
!</inputoutput>

!</subroutine>

    ! Is this a stationary solution based on a file?
    if (rrhsAssembly%ctype .eq. 3) then
      call lsysbl_releaseVector (rrhsAssembly%rrhsVector)
    else if (rrhsAssembly%ctype .eq. 4) then
      call lsysbl_releaseVector (rrhsAssembly%rrhsVector)
      call lsysbl_releaseVector (rrhsAssembly%rrhsVector2)
    end if

    rrhsAssembly%ctype = 0
    rrhsAssembly%icurrentRhs = -1
    rrhsAssembly%icurrentRhs2 = -1

  end subroutine

end module

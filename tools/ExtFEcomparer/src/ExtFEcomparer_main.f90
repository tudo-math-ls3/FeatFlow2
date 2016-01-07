module ExtFEcomparer_main

! Include basic Feat-2 modules
  use genoutput
  use paramlist
  use collection

  use ExtFEcomparer_typedefs

  use ExtFEcomparer_parameters
  use ExtFEcomparer_discretisation
  use ExtFEcomparer_paramtriang
  use ExtFEcomparer_core
  use ExtFEcomparer_vector
  use ExtFEcomparer_postprocessing
  use ExtFEcomparer_benchmark

  use triangulation

  implicit none

contains


subroutine start_ExtFEcomparer_main()

    ! figure out in which modus to start:
    ! normal one - which calculates what we
    ! set in the dat-files - or benchmark-modus -
    ! which calls a series of test to the program
    ! and ignores all user-inputs. This is usually
    ! something run ie in the regression tests.
    ! we outsources this in a specific file and module
    ! How do they work?
    ! Basically, they set up a FEM-Function on their
    ! own and configure the postprocessing structure.
    ! after that, they call the usual routines.
    ! Usually you can ignore this, the usual case
    ! is that you want to enter the normal modus,
    ! and that is the case.
    integer :: modus

    call ExtFE_getModus(modus)

    select case(modus)
        case(1,2)
            call ExtFEcomparer_benchMod(modus)
        case default
            call ExtFEcomparer_normal_mod()
    end select

end subroutine



  subroutine ExtFEcomparer_normal_mod()


    ! We need our problem structures
    type (t_problem) ,pointer :: p_rproblem1, p_rproblem2

    ! a postprocessing structure
    type (t_postprocessing), pointer :: p_rpostprocessing

    ! Print some info on the screen
    call output_line("Starting in normal modus")
    call output_separator(OU_SEP_MINUS)
    call output_lbrk()


    allocate(p_rproblem1)
    allocate(p_rproblem2)
    allocate(p_rpostprocessing)


    ! read in the parameter list.
    ! First step:
    ! Initialise the parameter list object. This creates an empty parameter list.
    call parlst_init (p_rproblem1%rparamlist)
    call parlst_init (p_rproblem2%rparamlist)

    ! Now read in the parameters
    call output_line("read in the parameters")
    call ExtFEcomparer_get_parameters(p_rproblem1%rparamlist)
    call ExtFEcomparer_get_parameters(p_rproblem2%rparamlist)

    ! Print out the parameters as some output - it is always good to have them
    call output_line ("Parameters:")
    call output_lbrk ()
    call parlst_info (p_rproblem1%rparamList)
    call output_lbrk ()

    ! Attention: Now both problem structures contain all information,
    ! we still need to figure out which one we need. That is done
    ! by 2 sections in "function_specific_settings.dat"

    ! Write parameters in the corresponding problem structures
    call output_line("Init the parameters")
    call output_lbrk()
    call ExtFEcomparer_init_parameters(p_rproblem1,"ExtFE-FIRST")
    call ExtFEcomparer_init_parameters(p_rproblem2,"ExtFE-SECOND")

    ! We have the parameters, so we recreate the parametrisation
    ! for every FE-function
    call output_line("Create the parametrisation and triangulation")
    call output_lbrk()
    call ExtFEcomparer_recreate_ParamTriang(p_rproblem1)
    call ExtFEcomparer_recreate_ParamTriang(p_rproblem2)


    ! We have the parametrisation and triangulation,
    ! so we need a discretisation now
    call output_line("Create the discretisation")
    call output_lbrk()
    call ExtFEcomparer_init_discretisation(p_rproblem1)
    call ExtFEcomparer_init_discretisation(p_rproblem2)

    call output_line("Load the vector")
    call output_lbrk()
    call ExtFEcomparer_load_vector(p_rproblem1)
    call ExtFEcomparer_load_vector(p_rproblem2)

    ! Init the postprocessing
    call output_line("Init the postprocessing structure")
    call output_lbrk()
    call ExtFEcomparer_init_postprocessing(p_rpostprocessing, &
                p_rproblem1%rparamlist)

    call output_line("Do the requested calculations")
    call ExtFEcomparer_calculate(p_rproblem1,p_rproblem2, &
                        p_rpostprocessing)

    ! Do the postprocessing
    call output_lbrk()
    call output_line("Do the file i/o: create requested files")
    call ExtFEcomparer_postprocess(p_rpostprocessing)

    ! General clean-up
    call parlst_done(p_rproblem1%rparamlist)
    call parlst_done(p_rproblem2%rparamlist)

    call ExtFEcomparer_doneDiscretisation(p_rproblem1)
    call ExtFEcomparer_doneDiscretisation(p_rproblem2)

    call ExtFEcomparer_doneParamTriang(p_rproblem1)
    call ExtFEcomparer_doneParamTriang(p_rproblem2)

    call ExtFEcomparer_done_postprocessing(p_rpostprocessing)

    call lsysbl_releaseVector(p_rproblem1%coeffVector)
    call lsysbl_releaseVector(p_rproblem2%coeffVector)

    if(associated(p_rproblem1%iElemList)) then
        deallocate(p_rproblem1%iElemList)
    end if
    if(associated(p_rproblem2%iElemList)) then
        deallocate(p_rproblem2%iElemList)
    end if
    deallocate(p_rproblem1)
    deallocate(p_rproblem2)
    deallocate(p_rpostprocessing)

end subroutine


subroutine ExtFE_getModus(modus)
    !<output>
    integer :: modus
    !</output>
		
		! internal variables
		character(LEN=ExtFE_STRLEN) :: masterDat

    type(t_parlist) :: rparlist


    ! Init the parlist
    call parlst_init(rparlist)

    ! We do no search in the dat-files as this modus will
    ! not be set up in the dat-files
		! The only file we search is the master.dat as this is
		! set up in the benchmark
		call ExtFEcomparer_getCmdlMasterDat(masterDat)

		! read the master-dat file
		call parlst_readfromfile(rparlist,masterDat)

		! Now we parse the command line. Thus the
		! commandline overrules the master.dat
    call ExtFEcomparer_parseCmdlArguments(rparlist)

    ! Now we search in the section "ExtFE-Modus" for the variable modus
    call parlst_getvalue_int(rparlist,"ExtFE-Modus","modus",modus,0 )

    ! We no longer need the parlist
    call parlst_done(rparlist)

end subroutine

end module

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

  use triangulation

  implicit none

contains



  subroutine start_ExtFEcomparer_main()


    ! We need our problem structures
    type (t_problem) ,pointer :: p_rproblem1, p_rproblem2

    ! a postprocessing structure
    type (t_postprocessing), pointer :: p_rpostprocessing

    allocate(p_rproblem1)
    allocate(p_rproblem2)
    allocate(p_rpostprocessing)


    ! read in the parameter list.
    ! First step:
    ! Initialise the parameter list object. This creates an empty parameter list.
    call parlst_init (p_rproblem1%rparamlist)
    call parlst_init (p_rproblem2%rparamlist)

    ! Now read in the parameters
    call output_line("read in the parameters from the files")
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

    deallocate(p_rproblem1)
    deallocate(p_rproblem2)
    deallocate(p_rpostprocessing)

end subroutine

end module

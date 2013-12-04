module ExtFEcomparer_main

! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  use fparser
  use paramlist
  use collection

  use triangulation
  use meshgeneration

  use element
  use cubature
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use bilinearformevaluation

  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use blockmatassemblystdop

  use ExtFEcomparer_typedefs

  use ExtFEcomparer_parameters
  use ExtFEcomparer_discretisation
  use ExtFEcomparer_parametrisation
  use ExtFEcomparer_core

  implicit none

contains

!<description>
! The routine start_ExtFEcomparer_main has only
! 1 purpose: to decide if it is a 2d or 3d-problem
! and call the according subroutine
! As there is no 3d-Code yet, this is for future
! purpose
!</description>

  subroutine start_ExtFEcomparer_main()

    call ExtFEcomparer_2D()

  end subroutine


!<description>
! This routine calls the according routines
! for a 2d-problem
!</description>

  subroutine ExtFEcomparer_2D()

    ! We need our problem structures
    type (t_problem) ,pointer :: p_rproblem1, p_rproblem2

    allocate(p_rproblem1)
    allocate(p_rproblem2)

    call collct_init(p_rproblem1%rcollection)
    call collct_init(p_rproblem2%rcollection)

    ! read in the parameter list.
    ! First step:
    ! Initialise the parameter list object. This creates an empty parameter list.
    call parlst_init (p_rproblem1%rparamlist)
    call parlst_init (p_rproblem2%rparamlist)

    ! Now read in the parameters
    call output_line("read in the parameters from the files")
    call ExtFEcomparer_get_parameters(p_rproblem1%rparamlist)
    call ExtFEcomparer_get_parameters(p_rproblem2%rparamlist)

    ! Attention: Now both problem structures contain all information,
    ! we still need to figure out which one we need. That is done
    ! by 2 sections in "configuration.dat"

    ! Write parameters in the corresponding problem structures
    ! and load the vector
    call output_line("Load the vectors and read out domain information")
    call ExtFEcomparer_init_parameters(p_rproblem1,"FIRST")
    call ExtFEcomparer_init_parameters(p_rproblem2,"SECOND")

    ! Recreate everything
    call output_line("Start to recreacte the parametrisations")
    call recreate_parametrisation_2D(p_rproblem1)
    call recreate_parametrisation_2D(p_rproblem2)
    call output_line("Start to recreate the discretisations")
    call reinitDiscretisation(p_rproblem1)
    call reinitDiscretisation(p_rproblem2)

    ! Do the actual computation
    call output_line("We have everything, so we do the actual computation now")
    call output_lbrk()
    call calculate_difference_2D(p_rproblem1,p_rproblem2)

    ! Clean up
    call ExtFEcomparer_doneDiscretisation(p_rproblem1)
    call ExtFEcomparer_doneParamTriang(p_rproblem1)
    call ExtFEcomparer_doneDiscretisation(p_rproblem2)
    call ExtFEcomparer_doneParamTriang(p_rproblem2)
    call parlst_done(p_rproblem1%rparamlist)
    call parlst_done(p_rproblem2%rparamlist)
    call ExtFEcomparer_doneParameters(p_rproblem1)
    call ExtFEcomparer_doneParameters(p_rproblem2)
    call lsysbl_releaseVector(p_rproblem1%coeffVector)
    call lsysbl_releaseVector(p_rproblem2%coeffVector)

    call collct_done(p_rproblem1%rcollection)
    call collct_done(p_rproblem2%rcollection)
    deallocate(p_rproblem1)
    deallocate(p_rproblem2)

  end subroutine


end module

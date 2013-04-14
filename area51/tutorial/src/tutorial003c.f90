!##############################################################################
!# Tutorial 003c: Using the parser
!##############################################################################

module tutorial003c

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage
  use fparser

  implicit none
  private
  
  public :: start_tutorial003c

contains

  ! ***************************************************************************

  subroutine start_tutorial003c

    ! Declare some variables
    integer :: ix,iy
    real(DP) :: dval1, dval2
    
    type(t_fparser) :: rparser
    character(LEN=1), dimension(2), parameter :: Svariables = (/ "X", "Y" /)
    real(DP), dimension(2) :: Dparams

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("Hello world. This is FEAT-2. Tutorial 003c")
    call output_separator (OU_SEP_MINUS)

    ! =================================
    ! Initialise the parser.
    ! =================================

    call fparser_init()

    ! =================================
    ! Evaluate the function in the
    ! points of a 5x5 mesh.
    ! =================================

    ! ----------------
    ! Parse a function
    ! ----------------
    call fparser_create (rparser, 2)
    call fparser_parseFunction (rparser, 1, "SIN(_PI*X) * SIN(_PI*Y)", Svariables)
    call fparser_parseFunction (rparser, 2, "16*X*Y*(1-X)*(1-Y)", Svariables)

    ! ----------------
    ! Evaluate
    ! ----------------
    do ix = 0,4
      do iy = 0,4
      
        ! Prepare the input variables
        Dparams = (/ real(ix,DP) / 4.0_DP, &
                     real(iy,DP) / 4.0_DP /)
             
        ! Evaluate        
        call fparser_evalFunction (rparser, 1, Dparams, dval1)
        call fparser_evalFunction (rparser, 2, Dparams, dval2)
        
        call output_line (" f ("  // trim(sys_sdL(Dparams(1),3)) // &
                          ", "    // trim(sys_sdL(Dparams(2),3)) // &
                          ") = "  // trim(sys_sdL(dval1,3))      // &
                          ",  " // &
                          " g ("  // trim(sys_sdL(Dparams(1),3)) // &
                          ", "    // trim(sys_sdL(Dparams(2),3)) // &
                          ") = "  // trim(sys_sdL(dval2,3)) )
      end do
    end do
    
    ! ----------------
    ! Release
    ! ----------------
    call fparser_release (rparser)
    
    ! =================================
    ! Cleanup
    ! =================================
    
    call fparser_done ()
    
  end subroutine

end module

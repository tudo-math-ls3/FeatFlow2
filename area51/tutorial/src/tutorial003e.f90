!##############################################################################
!# Tutorial 003e: Random number generator
!##############################################################################

module tutorial003e

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use random

  implicit none
  private
  
  public :: start_tutorial003e

contains

  ! ***************************************************************************

  subroutine start_tutorial003e
  
    ! Declare some variables
    type(t_random) :: rrng
    integer :: i
    integer, dimension(4) :: Ival
    real(DP), dimension(4) :: Dval

    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("Hello world. This is FEAT-2. Tutorial 003e")
    call output_separator (OU_SEP_MINUS)
    
    ! =================================
    ! Generate random numbers
    ! =================================
    
    ! Initialise the pseudo-random generator.
    call rng_init(rrng, 7635)
    
    ! Generate some integer numbers.
    do i=1,4
      call rng_get(rrng, Ival(i))
    end do

    ! Generate some double numbers.
    do i=1,4
      call rng_get(rrng, Dval(i))
    end do

    ! =================================
    ! Output
    ! =================================
    
    do i=1,4
      call output_line (" " // trim(sys_siL(Ival(i),10)), bnolinebreak=(i .ne. 4))
    end do

    do i=1,4
      call output_line (" " // sys_sdL(Dval(i),2), bnolinebreak=(i .ne. 4))
    end do

  end subroutine

end module

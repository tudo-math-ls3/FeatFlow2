!##############################################################################
!# Tutorial 001a: Hello world
!##############################################################################

module tutorial001a

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  
  use storage
  use spatialdiscretisation

  implicit none
  private
  
  public :: start_tutorial001a

contains

  ! ***************************************************************************

  subroutine start_tutorial001a

    ! Print a message
    call output_line ("Hello world. This is FEAT-2. Tutorial 001a.")
    call output_line ("Zahlen:" // sys_siL1(10) // sys_siL1(20) )

  end subroutine

!************************************************************************

!<function>

  function sys_siL1(ivalue) result(soutput)

!<description>
    ! This routine converts an integer value to a string of length idigits,
    ! filled up with white spaces.
!</description>

!<input>
    ! value to be converted
    integer, intent(in) :: ivalue
!</input>

!<result>
    ! String representation of the value (left-aligned),
    ! fixed length of idigits characters
    character (len=len_trim(adjustl(sys_si2(ivalue,16)))) :: soutput
!</result>

!</function>

    soutput = adjustl(sys_si2(ivalue,16))
  end function

!************************************************************************

!<function>

  pure character (len=32) function sys_si2(ivalue, idigits) result(soutput)

!<description>
    ! This routine converts an integer value to a string of length idigits.
!</description>

!<result>
    ! String representation of the value, filled with white spaces.
    ! At most 32 characters supported.
!</result>

!<input>

    ! value to be converted
    integer, intent(in) :: ivalue

    !number of decimals
    integer, intent(in) :: idigits
!</input>
!</function>

    character (len=16) :: sformat
    character (len=2)  :: saux
    ! idigits can not be simply adjusted to 16 because some compilers
    ! do not accept that idigits is changed within this function if
    ! the function is called with a hard-coded integer instead of a
    ! variable, i.e.
    !   sys_sli0(foo, 1)
    ! would result in a crash
    if (idigits .gt. 16) then
!      write(6, *) "*** WARNING! Too many decimal places requested in sys_si! ***"
      write(saux, "(i2)") 16
    else if (idigits .lt. 10) then
      write(saux, "(i1)") idigits
    else
      write(saux, "(i2)") idigits
    endif

    sformat = "(i" // trim(saux) // ")"
    write (unit = soutput, fmt = trim(sformat)) ivalue

  end function sys_si2

end module

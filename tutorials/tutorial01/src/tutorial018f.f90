!##############################################################################
!# Tutorial 018f: Collection handling, passing user defined objects ("void*")
!##############################################################################

module tutorial018f

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput

  use collection

  implicit none
  private
  
  public :: start_tutorial018f
  
  ! A small local type with some data
  type t_problem
    integer, dimension(5) :: Iarray
    real(DP), dimension(5) :: Darray
  end type

contains

  ! ************************************************************************
  ! AUXILIARY SUBROUTINES
  ! ************************************************************************

  ! Casts the user defined object to a generic object
  subroutine t_problem_cast(rproblem, rgenericObject)

  ! User defined object
  type(t_problem), intent(in), target :: rproblem

  ! The generic object
  type(t_genericObject), intent(out) :: rgenericObject

#define CASTTYPE t_problem
#define CASTVAR  rproblem
#define CASTGOBJ rgenericObject
#include <casttogenobject.h>
    
  end subroutine

  ! ************************************************************************

  ! Casts a generic object to a pointer to a t_problem object
  subroutine t_problem_uncast(rgenericObject, p_rproblem)

  ! The generic object
  type(t_genericObject), intent(in) :: rgenericObject

  ! User defined object
  type(t_problem), pointer :: p_rproblem

#define CASTTYPE t_problem
#define CASTVAR  p_rproblem
#define CASTGOBJ rgenericObject
#include <uncastfromgenobject.h>

  end subroutine

  ! ************************************************************************
  ! ************************************************************************
  ! ************************************************************************

  subroutine printarray (rcollection)
  type(t_collection), intent(inout) :: rcollection
  
    type(t_problem), pointer :: p_rproblem
    integer :: i
    
    ! Get out type from the quick-access generic object
    call t_problem_uncast(rcollection%RgenObjectQuickAccess(1),p_rproblem)
    
    do i=1,ubound(p_rproblem%Iarray,1)
      call output_line (sys_si(p_rproblem%Iarray(i),5))
    end do

    do i=1,ubound(p_rproblem%Darray,1)
      call output_line (sys_sd(p_rproblem%Darray(i),10))
    end do

  end subroutine

  ! ***************************************************************************

  subroutine start_tutorial018f

    ! Declare some variables.
    type(t_collection) :: rcollection
    type(t_problem), target :: rproblem
    
    ! Print a message
    call output_lbrk()
    call output_separator (OU_SEP_STAR)
    call output_line ("This is FEAT-2. Tutorial 018f")
    call output_separator (OU_SEP_MINUS)

    ! Fill our user defined type with data
    rproblem%Iarray(1:5) = (/ 5,4,3,2,1 /)
    rproblem%Darray(1:5) = (/ 10.0_DP,9.0_DP,8.0_DP,7.0_DP,6.0_DP /)
    
    ! Put it to the quick-access generic object of the collection.
    call t_problem_cast(rproblem, rcollection%RgenObjectQuickAccess(1))
    
    ! Pass the collection to the subroutine
    call printarray (rcollection)
    
  end subroutine

end module

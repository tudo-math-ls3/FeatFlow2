    ! Casts an a "generic object" into an arbitrary structure.
    ! This is the body of an "uncast" subroutine.
    !
    ! Call:
    !
    !  ! Casts the user defined object to a generic object
    !  subroutine tmystruc_uncast(rmyobject, p_rgenericObject)
    !
    !  ! The generic object
    !  type(t_genericObject), intent(in) :: rgenericObject
    !
    !  ! User defined object
    !  type(t_myStructure), pointer :: p_rmyStructure
    !
    !  #define CASTTYPE t_myStructure
    !  #define CASTVAR  p_rmyStructure
    !  #define CASTGOBJ rgenericObject
    !  #include "uncastfromgenobject.h"
    !
    !  end subroutine

    ! Internal data structure
    type t_void_ptr
      type(CASTTYPE), pointer :: p_robj => null()
    end type

    ! Internal variables
    type(t_void_ptr) :: rptr

    ! Transfer the generic object to the void pointer structure
    rptr = transfer(CASTGOBJ, rptr)

    ! Unwrap list from void pointer structure
    CASTVAR => rptr%p_robj

! Remove Defines
#undef CASTTYPE
#undef CASTVAR
#undef CASTGOBJ

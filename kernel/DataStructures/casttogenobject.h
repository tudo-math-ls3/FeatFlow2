    ! Casts an arbitrary structure into a "generic object".
    ! This is the body of a "cast" subroutine.
    !
    ! Call:
    !
    !  ! Casts the user defined object to a generic object
    !  subroutine tmystruc_uncast(rmyobject, rgenericObject)
    !
    !  ! User defined object
    !  type(t_myStructure), intent(in), target :: rmyStructure
    !
    !  ! The generic object
    !  type(t_genericObject), intent(out) :: rgenericObject
    !
    !  #define CASTTYPE t_myStructure
    !  #define CASTVAR  rmyStructure
    !  #define CASTGOBJ rgenericObject
    !  #include "casttogenobject.h"
    !
    !  end subroutine

    ! Internal data structure
    type t_void_ptr
      type(CASTTYPE), pointer :: p_robj => null()
    end type t_void_ptr
    
    ! Internal variables
    type(t_void_ptr) :: rptr

    ! Wrap list by void pointer structure
    rptr%p_robj => CASTVAR
    
    ! Transfer the void pointer structure to the generic object
    CASTGOBJ = transfer(rptr, CASTGOBJ)

! Remove Defines
#undef CASTTYPE
#undef CASTVAR
#undef CASTGOBJ

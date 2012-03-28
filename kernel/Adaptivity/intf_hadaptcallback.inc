  ! This file describes the interfaces of the callback routines
  ! ***************************************************************************

  interface

!<subroutine>

    subroutine fcb_hadaptCallback(iOperation, rcollection)

      use fsystem
      use collection
      use triangulation

!<input>
      ! Identifier for the grid modification operation
      integer, intent(in) :: iOperation
!</input>

!<inputoutput>
      ! Collection
      type(t_collection), intent(inout) :: rcollection
!</inputoutput>
!</subroutine>

    end subroutine fcb_hadaptCallback

  end interface

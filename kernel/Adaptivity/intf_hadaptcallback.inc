  ! This file describes the interfaces of the callback routines
  ! ***************************************************************************

  interface

!<subroutine>

  subroutine fcb_hadaptCallback(rcollection,iOperation,Ivertices,Ielements)

  use fsystem
  use collection
  use triangulation

!<input>
  ! Collection
  type(t_collection), intent(IN) :: rcollection
 
  ! Identifier for the grid modification operation
  integer, intent(IN) :: iOperation

  ! Array of vertices involved in the adaptivity step
  integer, dimension(:), intent(IN) :: Ivertices

  ! Array of elements involved in the adaptivity step
  integer, dimension(:), intent(IN) :: Ielements
!</input>
!</subroutine>

  end subroutine

  end interface

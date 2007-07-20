  ! This file describes the interfaces of the callback routines
  ! ***************************************************************************

  INTERFACE

!<subroutine>

  SUBROUTINE fcb_hadaptCallback(rcollection,iOperation,Ivertices,Ielements)

  USE fsystem
  USE collection
  USE triangulation

!<input>
  ! Collection
  TYPE(t_collection), INTENT(IN) :: rcollection
 
  ! Identifier for the grid modification operation
  INTEGER, INTENT(IN) :: iOperation

  ! Array of vertices involved in the adaptivity step
  INTEGER(PREC_VERTEXIDX), DIMENSION(:), INTENT(IN) :: Ivertices

  ! Array of elements involved in the adaptivity step
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: Ielements
!</input>
!</subroutine>

  END SUBROUTINE

  END INTERFACE

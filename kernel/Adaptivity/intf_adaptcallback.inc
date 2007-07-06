  ! This file describes the interfaces of the callback routines
  ! ***************************************************************************

  INTERFACE

!<subroutine>

  SUBROUTINE fcb_adjustDimension(neq)

  USE fsystem

!<input>
   ! total number of equations
   INTEGER, INTENT(IN) :: neq
!</input>
!</subroutine>

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE fcb_insertVertex(nsvt,ksvt,ivt)

  USE fsystem

!<input>
    ! number of surrounding vertices
    INTEGER, INTENT(IN) :: nsvt

    ! array of numbers of surrounding vertices
    INTEGER, INTENT(IN), DIMENSION(nsvt) :: ksvt

    ! number of new vertex
    INTEGER, INTENT(IN) :: ivt
!</input>
!</subroutine>

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE fcb_removeVertex(ivt)

  USE fsystem

!<input>
    ! number of removed vertex
    INTEGER, INTENT(IN) :: ivt
!</input>
!</subroutine>

  END SUBROUTINE

  END INTERFACE

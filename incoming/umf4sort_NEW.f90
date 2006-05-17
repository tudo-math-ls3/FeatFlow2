!##############################################################################
!# ****************************************************************************
!# <name> umf4sort </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines necessary for making the storage of the 
!# matrices compatible between the UMFPACK4 solver format (0-based) and the 
!# original FEAT format (1-based).
!#
!# </purpose>
!##############################################################################

MODULE umf4sort

  USE fsystem

  IMPLICIT NONE

CONTAINS

!<subroutine>

  SUBROUTINE sortCSR (Da, Kcol, Kld, neq)

!<description>

  ! Sorts the entries in each row of the matrix in ascending order.
  ! This is necessary for the UMFPACK4 solver, which expects a matrix 
  ! in this format.
  ! The input matrix is assumed to be in storage technique 7 with the 
  ! first element in each row to be the diagonal element!

!</description>

!<inputoutput>

  ! On input:  the matrix entries/column numbers to be resorted
  ! On output: the resorted matrix entries/column numbers
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: Da
  INTEGER(PREC_MATRIDX), DIMENSION(:), INTENT(INOUT) :: Kcol
  INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kld

  ! Dimension of the matrix
  INTEGER(I32), INTENT(IN) :: neq

!</inputoutput>

!</subroutine>

  ! local variables
  REAL(DP) :: aux
  INTEGER(I32) :: i, j

  ! loop through each row
  DO i = 1, neq

    ! Take the diagonal element
    aux = Da(Kld(i))

    ! Loop through each column in this row.
    ! Shift every entry until the diagonal is reached.
    DO j = Kld(i)+1, Kld(i+1)-1

      ! Check if we reached the position of the diagonal entry...
      IF (Kcol(J)>i) EXIT

        Kcol(j-1) = KCOL(j)
        Da(j-1) = Da(j)

    END DO

    ! If we have reached the diagonal, we can stop and save our
    ! diagonal entry from the first position there. The rest of the
    ! line is in ascending order according to the specifications of
    ! storage technique 7.

    Kcol(j-1) = i
    Da(j-1) = aux

  END DO
          
  END SUBROUTINE

! ***************************************************************************

!<subroutine>

  SUBROUTINE M7IDSH (Kcol,Kld,neq)

!<description>

  ! Performs an index shift "-1" on all entries on the KCOL and KLD
  ! array. This is due to convert these arrays to 0-based for the
  ! use with UMFPACK4

!</description>

!<inputoutput>

  ! On input:  the matrix column/row numbers to be shifted
  ! On output: the shifted matrix column/row numbers
  INTEGER(I32), DIMENSION(:), INTENT(INOUT) :: Kcol, Kld

  ! Dimension of the matrix
  INTEGER(I32), INTENT(IN) :: neq

!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER(I32) :: i, na
       
    na = Kld(neq+1)-1
      
    DO I=1,na
    Kcol (i) = Kcol(i)-1
    END DO

    DO i=1,neq+1
    Kld(i) = Kld(i)-1
    END DO

  END SUBROUTINE
      

!<subroutine>

  SUBROUTINE M7IDSB (KCOL,KLD,NEQ)

!<description>

  ! Performs an index shift "+1" on all entries on the KCOL and KLD
  ! array. Can be used for debugging purposes to transform a matrix
  ! back after being modified by M7IDSH.

!</description>

!<inputoutput>

  ! On input:  the matrix column/row numbers to be shifted
  ! On output: the shifted matrix column/row numbers
  INTEGER(I32), DIMENSION(:), INTENT(INOUT) :: Kcol, Kld

  ! Dimension of the matrix
  INTEGER(I32), INTENT(IN) :: neq

!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER(I32) :: i, na

  na = Kld(neq+1)

  DO i=1,na
    Kcol (i) = Kcol(i)+1
  END DO

  DO i=1,neq+1
    Kld(i) = Kld(i)+1
  END DO

  END SUBROUTINE
  
  END MODULE
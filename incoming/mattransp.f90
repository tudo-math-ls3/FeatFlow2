!#########################################################################
!# ***********************************************************************
!# <name> mattransp </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains all routines necessary for resorting matrix
!# or vector entries.
!#
!# For storing matrices, we use the CSR format, realised as ''Format-7''
!# matrices. Format-7 is a variant of the CSR-format used in old FEAT
!# programs, where the diagonal element is stored first in each line.
!#
!# The following routines can be found in this module:
!#
!# 1.) mattransp_transposeMatrixStructure
!#     -> Transpose matrix structure
!#
!# 2.) mattransp_transposeMatrix
!#     -> Transpose matrix
!#
!# 3.) mattransp_pivotizeMatrix
!#     -> Pivotize matrix
!#
!# </purpose>
!#########################################################################

MODULE mattransp

  USE fsystem
  
  IMPLICIT NONE
  
  CONTAINS
  
!<subroutine>
  SUBROUTINE mattransp_transposeMatrixStructure (nrow, ncol, Icol, Irow, &
                                                 Itmp, IcolDest, IrowDest)
  
  !<description>
    ! This routine accepts the structure of a structure-7 or structure-9
    ! matrix and creates the structure of the transposed matrix from it.
    !
    ! The resulting structure is not pivoted anymore with the diagonal
    ! entry in front if a structure-7 matrix is to be transposed!
  !</description>
    
  !<input>
    ! Number of rows in the source matrix
    INTEGER(I32), INTENT(IN) :: nrow
    
    ! Number of columns in the source matrix
    INTEGER(I32), INTENT(IN) :: ncol
    
    ! Column structure of the source matrix
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: Icol
    
    ! Row structure of the source matrix
    INTEGER(I32), DIMENSION(nrow+1), INTENT(IN) :: Irow
  !</input>
    
  !<output>
    ! Column structure of the destination matrix
    ! The array must be of the same size as Icol!
    INTEGER(I32), DIMENSION(:), INTENT(OUT) :: IcolDest

    ! Row structure of the destination matrix
    INTEGER(I32), DIMENSION(ncol+1), INTENT(OUT) :: IrowDest
  !</output>
    
  !<inputoutput>
    ! Auxiliary array of size ncol
    INTEGER(I32), DIMENSION(ncol), INTENT(INOUT) :: Itmp
  !</inputoutput>
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, isize, icolumn, ncolumn
    
    ! determin the number of matrix entries
    isize = Irow(nrow+1)-1
    
    ! clear row struxture of the destination matrix
    DO i=1, ncol+1
      IrowDest(i) = 0
    END DO
    
    ! Count how many entries <> 0 are in each column. Note this into
    ! the IrowDest array shifted by 1.
    DO i=1, isize
      IrowDest(Icol(i)+1) = IrowDest(Icol(i)+1)+1
    END DO
    
    ! Now build the final IrowDest by adding up the IrowDest entries.
    IrowDest(1) = 1
    DO i=1, ncol+1
      IrowDest(i) = IrowDest(i)+IrowDest(i-1)
    END DO
    
    ! Now IcolDest must be created. This requires another loop trough
    ! the matrix structure. Itmp receives the index of how many entries
    ! have been written to each row.
    
    ! clear auxiliary vector
    DO i=1, ncol
      Itmp(i) = 0
    END DO

    DO i=1, nrow
      ncolumn = Irow(i+1)-Irow(i)
      DO j=1, ncolumn
        ! Get the column of the item in question -> new row number.
        icolumn = Icol(Irow(i)+j-1)
        ! Rows get columns by transposing, therefore note i as column
        ! number in IcolDest
        IcolDest(IrowDest(icolumn)+Itmp(icolumn)) = i
        ! Increment running index of that row
        Itmp(icolumn) = Itmp(icolumn)+1
      END DO
    END DO
    
  END SUBROUTINE mattransp_transposeMatrixStructure

!<subroutine>
  SUBROUTINE mattransp_transposeMatrix (nrow, ncol, Da, Icol, Irow, &
                                        Itmp, DaDest, IcolDest, IrowDest)
  
  !<description>
    ! This routine accepts a structure-7 or structure-9
    ! matrix and creates the transposed matrix from it.
    !
    ! The resulting matrix is not pivoted anymore with the diagonal
    ! entry in front if a structure-7 matrix is to be transposed!
  !</description>
    
  !<input>
    ! Number of rows in the source matrix
    INTEGER(I32), INTENT(IN) :: nrow
    
    ! Number of columns in the source matrix
    INTEGER(I32), INTENT(IN) :: ncol
    
    ! The entries of the source matrix
    REAL(DP), DIMENSION(:), INTENT(IN) :: Da
    
    ! Column structure of the source matrix
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: Icol
    
    ! Row structure of the source matrix
    INTEGER(I32), DIMENSION(nrow+1), INTENT(IN) :: Irow
  !</input>
    
  !<output>
    ! The entries of the destination matrix
    ! The array must be of the same size as Da or Icol
    REAL(DP), DIMENSION(:), INTENT(OUT) :: DaDest
    
    ! Column structure of the destination matrix
    ! The array must be of the same size as Icol!
    INTEGER(I32), DIMENSION(:), INTENT(OUT) :: IcolDest

    ! Row structure of the destination matrix
    INTEGER(I32), DIMENSION(ncol+1), INTENT(OUT) :: IrowDest
  !</output>
    
  !<inputoutput>
    ! Auxiliary array of size ncol
    INTEGER(I32), DIMENSION(ncol), INTENT(INOUT) :: Itmp
  !</inputoutput>
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, isize, icolumn, ncolumn
    
    ! determin the number of matrix entries
    isize = Irow(nrow+1)-1
    
    ! clear row struxture of the destination matrix
    DO i=1, ncol+1
      IrowDest(i) = 0
    END DO
    
    ! Count how many entries <> 0 are in each column. Note this into
    ! the IrowDest array shifted by 1.
    DO i=1, isize
      IrowDest(Icol(i)+1) = IrowDest(Icol(i)+1)+1
    END DO
    
    ! Now build the final IrowDest by adding up the IrowDest entries.
    IrowDest(1) = 1
    DO i=1, ncol+1
      IrowDest(i) = IrowDest(i)+IrowDest(i-1)
    END DO
    
    ! Now IcolDest must be created. This requires another loop trough
    ! the matrix structure. Itmp receives the index of how many entries
    ! have been written to each row.
    
    ! clear auxiliary vector
    DO i=1, ncol
      Itmp(i) = 0
    END DO

    DO i=1, nrow
      ncolumn = Irow(i+1)-Irow(i)
      DO j=1, ncolumn
        ! Get the column of the item in question -> new row number.
        icolumn = Icol(Irow(i)+j-1)
        ! Rows get columns by transposing, therefore note i as column
        ! number in IcolDest
        IcolDest(IrowDest(icolumn)+Itmp(icolumn)) = i
        ! Copy the matrix entry:
        DaDest(IrowDest(icolumn)+Itmp(icolumn)) = Da(Irow(i)+j-1)
        ! Increment running index of that row
        Itmp(icolumn) = Itmp(icolumn)+1
      END DO
    END DO
    
  END SUBROUTINE mattransp_transposeMatrix

!<subroutine>
  SUBROUTINE mattransp_pivotizeMatrix (neq, Da, Icol, Irow, berror)
  
  !<description>
    ! This routine repivots a matrix in structure 7: The diagonal element
    ! of each row is moved to the front.
  !</description>
    
  !<input>
    ! Number of rows/columns/equations
    INTEGER(I32), INTENT(IN) :: neq
  !</input>
  
  !<output>
    ! If .TRUE. a diagonal element was not found
    LOGICAL, OPTIONAL, INTENT(OUT) :: berror
  !</output>
        
  !<inputoutput>
    ! Matrix entries
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: Da
    
    ! Matrix column structure
    INTEGER(I32), DIMENSION(:), INTENT(INOUT) :: Icol
    
    ! Matrix row structure
    INTEGER(I32), DIMENSION(neq+1), INTENT(INOUT) :: Irow
  !</inputoutput>
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j
    LOGICAL :: bnodiag
    
    bnodiag = .FALSE.
    
    ! Loop through the rows
    DO i=1, neq
      ! Is the line already pivoted? Most unlikely, but we test for sure
      IF (Icol(Irow(i)) .NE. i) THEN

        ! Find the position of the diagonal element
        bnodiag = .TRUE.
        DO j=Irow(i), Irow(i+1)-1
          IF (j .EQ. i) THEN
            bnodiag = .FALSE.
            EXIT
          END IF
        END DO
        ! Oops, diagonal element not found - cancel
        IF (bnodiag) EXIT
        
        ! Ringshift the slice Icol(Irow(i):j) by 1 to the right.
        ! This puts the diagonal element to the front
        ! The same operation is also done for Da(Irow(i):j)
        Icol(Irow(i):j) = CSHIFT(Icol(Irow(i):j),-1)
        Da  (Irow(i):j) = CSHIFT(Da  (Irow(i):j),-1)
      END IF
    END DO

    ! Report error status
    IF (PRESENT(berror)) berror = bnodiag    
  
  END SUBROUTINE mattransp_pivotizeMatrix

END MODULE

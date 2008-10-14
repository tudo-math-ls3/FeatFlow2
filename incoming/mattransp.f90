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

module mattransp

  use fsystem
  
  implicit none
  
  contains
  
!<subroutine>
  subroutine mattransp_transposeMatrixStructure (nrow, ncol, Icol, Irow, &
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
    integer(I32), intent(IN) :: nrow
    
    ! Number of columns in the source matrix
    integer(I32), intent(IN) :: ncol
    
    ! Column structure of the source matrix
    integer(I32), dimension(:), intent(IN) :: Icol
    
    ! Row structure of the source matrix
    integer(I32), dimension(nrow+1), intent(IN) :: Irow
  !</input>
    
  !<output>
    ! Column structure of the destination matrix
    ! The array must be of the same size as Icol!
    integer(I32), dimension(:), intent(OUT) :: IcolDest

    ! Row structure of the destination matrix
    integer(I32), dimension(ncol+1), intent(OUT) :: IrowDest
  !</output>
    
  !<inputoutput>
    ! Auxiliary array of size ncol
    integer(I32), dimension(ncol), intent(INOUT) :: Itmp
  !</inputoutput>
!</subroutine>
    
    !local variables
    integer(I32) :: i, j, isize, icolumn, ncolumn
    
    ! determin the number of matrix entries
    isize = Irow(nrow+1)-1
    
    ! clear row struxture of the destination matrix
    do i=1, ncol+1
      IrowDest(i) = 0
    end do
    
    ! Count how many entries <> 0 are in each column. Note this into
    ! the IrowDest array shifted by 1.
    do i=1, isize
      IrowDest(Icol(i)+1) = IrowDest(Icol(i)+1)+1
    end do
    
    ! Now build the final IrowDest by adding up the IrowDest entries.
    IrowDest(1) = 1
    do i=1, ncol+1
      IrowDest(i) = IrowDest(i)+IrowDest(i-1)
    end do
    
    ! Now IcolDest must be created. This requires another loop trough
    ! the matrix structure. Itmp receives the index of how many entries
    ! have been written to each row.
    
    ! clear auxiliary vector
    do i=1, ncol
      Itmp(i) = 0
    end do

    do i=1, nrow
      ncolumn = Irow(i+1)-Irow(i)
      do j=1, ncolumn
        ! Get the column of the item in question -> new row number.
        icolumn = Icol(Irow(i)+j-1)
        ! Rows get columns by transposing, therefore note i as column
        ! number in IcolDest
        IcolDest(IrowDest(icolumn)+Itmp(icolumn)) = i
        ! Increment running index of that row
        Itmp(icolumn) = Itmp(icolumn)+1
      end do
    end do
    
  end subroutine mattransp_transposeMatrixStructure

!<subroutine>
  subroutine mattransp_transposeMatrix (nrow, ncol, Da, Icol, Irow, &
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
    integer(I32), intent(IN) :: nrow
    
    ! Number of columns in the source matrix
    integer(I32), intent(IN) :: ncol
    
    ! The entries of the source matrix
    real(DP), dimension(:), intent(IN) :: Da
    
    ! Column structure of the source matrix
    integer(I32), dimension(:), intent(IN) :: Icol
    
    ! Row structure of the source matrix
    integer(I32), dimension(nrow+1), intent(IN) :: Irow
  !</input>
    
  !<output>
    ! The entries of the destination matrix
    ! The array must be of the same size as Da or Icol
    real(DP), dimension(:), intent(OUT) :: DaDest
    
    ! Column structure of the destination matrix
    ! The array must be of the same size as Icol!
    integer(I32), dimension(:), intent(OUT) :: IcolDest

    ! Row structure of the destination matrix
    integer(I32), dimension(ncol+1), intent(OUT) :: IrowDest
  !</output>
    
  !<inputoutput>
    ! Auxiliary array of size ncol
    integer(I32), dimension(ncol), intent(INOUT) :: Itmp
  !</inputoutput>
!</subroutine>
    
    !local variables
    integer(I32) :: i, j, isize, icolumn, ncolumn
    
    ! determin the number of matrix entries
    isize = Irow(nrow+1)-1
    
    ! clear row struxture of the destination matrix
    do i=1, ncol+1
      IrowDest(i) = 0
    end do
    
    ! Count how many entries <> 0 are in each column. Note this into
    ! the IrowDest array shifted by 1.
    do i=1, isize
      IrowDest(Icol(i)+1) = IrowDest(Icol(i)+1)+1
    end do
    
    ! Now build the final IrowDest by adding up the IrowDest entries.
    IrowDest(1) = 1
    do i=1, ncol+1
      IrowDest(i) = IrowDest(i)+IrowDest(i-1)
    end do
    
    ! Now IcolDest must be created. This requires another loop trough
    ! the matrix structure. Itmp receives the index of how many entries
    ! have been written to each row.
    
    ! clear auxiliary vector
    do i=1, ncol
      Itmp(i) = 0
    end do

    do i=1, nrow
      ncolumn = Irow(i+1)-Irow(i)
      do j=1, ncolumn
        ! Get the column of the item in question -> new row number.
        icolumn = Icol(Irow(i)+j-1)
        ! Rows get columns by transposing, therefore note i as column
        ! number in IcolDest
        IcolDest(IrowDest(icolumn)+Itmp(icolumn)) = i
        ! Copy the matrix entry:
        DaDest(IrowDest(icolumn)+Itmp(icolumn)) = Da(Irow(i)+j-1)
        ! Increment running index of that row
        Itmp(icolumn) = Itmp(icolumn)+1
      end do
    end do
    
  end subroutine mattransp_transposeMatrix

!<subroutine>
  subroutine mattransp_pivotizeMatrix (neq, Da, Icol, Irow, berror)
  
  !<description>
    ! This routine repivots a matrix in structure 7: The diagonal element
    ! of each row is moved to the front.
  !</description>
    
  !<input>
    ! Number of rows/columns/equations
    integer(I32), intent(IN) :: neq
  !</input>
  
  !<output>
    ! If .TRUE. a diagonal element was not found
    logical, optional, intent(OUT) :: berror
  !</output>
        
  !<inputoutput>
    ! Matrix entries
    real(DP), dimension(:), intent(INOUT) :: Da
    
    ! Matrix column structure
    integer(I32), dimension(:), intent(INOUT) :: Icol
    
    ! Matrix row structure
    integer(I32), dimension(neq+1), intent(INOUT) :: Irow
  !</inputoutput>
!</subroutine>
    
    !local variables
    integer(I32) :: i, j
    logical :: bnodiag
    
    bnodiag = .false.
    
    ! Loop through the rows
    do i=1, neq
      ! Is the line already pivoted? Most unlikely, but we test for sure
      if (Icol(Irow(i)) .ne. i) then

        ! Find the position of the diagonal element
        bnodiag = .true.
        do j=Irow(i), Irow(i+1)-1
          if (j .eq. i) then
            bnodiag = .false.
            exit
          end if
        end do
        ! Oops, diagonal element not found - cancel
        if (bnodiag) exit
        
        ! Ringshift the slice Icol(Irow(i):j) by 1 to the right.
        ! This puts the diagonal element to the front
        ! The same operation is also done for Da(Irow(i):j)
        Icol(Irow(i):j) = cshift(Icol(Irow(i):j),-1)
        Da  (Irow(i):j) = cshift(Da  (Irow(i):j),-1)
      end if
    end do

    ! Report error status
    if (present(berror)) berror = bnodiag    
  
  end subroutine mattransp_pivotizeMatrix

end module

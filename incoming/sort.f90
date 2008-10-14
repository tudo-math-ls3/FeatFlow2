!#########################################################################
!# ***********************************************************************
!# <name> sort </name>
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
!# 1.) sort_vecSort_double
!#     -> Resorting of a vector
!#
!# 2.) sort_vecSort_single
!#     -> Resorting of a vector
!#
!# 3.) sort_matSort_double
!#     -> Resorting of a matrix
!#
!# 4.) sort_matSort_single
!#     -> Resorting of a matrix
!#
!# 5.) sort_matResort
!#     -> Resorting back the matrix structure
!#
!# 6.) sort_crSort
!#     -> Applying bubble sort at column and row entries (only internally used)
!#
!# </purpose>
!#########################################################################

module sort

  use fsystem
  
  implicit none
  
  !<constants>
  !<constantblock description="">
  
  ! Auxiliary column array size
  integer, parameter :: SORT_BUFFERSIZE = 100
  
  !</constantblock>
  !</constants>
  
  contains

!<subroutine>
  subroutine sort_vecSort_double (Dx, Dd, Itr, neq)
  
  !<description>
    ! Resorts the entries in the vector Dx corresponding to Itr.
    ! The result is written to Dd.
  !</description>
    
  !<input>

    ! Number of equations
    integer(I32) :: neq
  
    ! Source vector to be sorted
    real(DP), dimension(neq), intent(IN) :: Dx
    
    ! Array with permutation of 1..neq
    integer(I32), dimension(neq), intent(IN) :: Itr
  !</input>
    
  !<output>
    ! The resorted vector
    real(DP), dimension(neq), intent(OUT) :: Dd
  !</output>
    
!</subroutine>
    
    ! local variable
    integer(I32) :: ieq
    
    do ieq=1, neq
      Dd(ieq) = Dx(Itr(ieq))
    end do
  
  end subroutine sort_vecSort_double

!<subroutine>
  subroutine sort_vecSort_single (Fx, Fd, Itr, neq)
  
  !<description>
    ! Resorts the entries in the vector Fx corresponding to Itr.
    ! The result is written to Fd.
  !</description>
    
  !<input>

    ! Number of equations
    integer(I32) :: neq
  
    ! Source vector to be sorted
    real(SP), dimension(neq), intent(IN) :: Fx
    
    ! Array with permutation of 1..neq
    integer(I32), dimension(neq), intent(IN) :: Itr
  !</input>
    
  !<output>
    ! The resorted vector
    real(SP), dimension(neq), intent(OUT) :: Fd
  !</output>
    
!</subroutine>
    
    ! local variable
    integer(I32) :: ieq
    
    do ieq=1, neq
      Fd(ieq) = Fx(Itr(ieq))
    end do
  
  end subroutine sort_vecSort_single

!<subroutine>
  subroutine sort_matSort_double (Da, DaH, Icol, IcolH, &
                                  Ild, IldH, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! Double precision version, storage technique 7 only.
  !</description>
    
  !<input>
    
    ! Number of equations
    integer(I32), intent(IN) :: neq
    
    ! Source matrix
    real(DP), dimension(:), intent(IN) :: DaH
    
    ! Column structure of source matrix
    integer(I32), dimension(:), intent(IN) :: IcolH
    
    ! Row positions of source matrix
    integer(I32), dimension(neq+1), intent(IN) :: IldH
    
    ! Permutation of 1..neq describing the sorting
    integer(I32), dimension(neq), intent(IN) :: Itr1
    
    ! Permutation of 1..neq describing the sorting back
    integer(I32), dimension(neq), intent(IN) :: Itr2

  !</input>
    
  !<output>

    ! Destination matrix
    real(DP), dimension(:), intent(OUT) :: Da

    ! Column structure of destination matrix
    integer(I32), dimension(:), intent(OUT) :: Icol
    
    ! Row positions of destination matrix
    integer(I32), dimension(neq+1), intent(OUT) :: Ild
    
  !</output>
    
!</subroutine>
    
    !local variables
    integer(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx
    integer(I32), dimension(SORT_BUFFERSIZE) :: Ih1, Ih2
  
    ildIdx = 1
    do i=1, neq
      iidx = Itr1(i)
      Da(ildIdx) = DaH(IldH(iidx))
      Ild(i) = ildIdx
      Icol(ildIdx) = i
      ildIdx = ildIdx + 1
      
      ih1Idx = IldH(iidx)+1
      ih2Idx = IldH(iidx+1)-1
      isize = ih2Idx - ih1Idx + 1
      do j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      end do
      call sort_crSort(Ih1,Ih2,isize)
      do j=1, isize
        Da(ildIdx) = DaH(Ih2(j))
        Icol(ildIdx) = Ih1(j) !Achtung: Ist dies so richtig ? Ih1(j)->Ih2(j) ?
        ildIdx = ildIdx + 1
      end do
    end do
    Ild(neq+1) = IldH(neq+1)

  end subroutine sort_matSort_double

!<subroutine>
  subroutine sort_matSort_single (Fa, FaH, Icol, IcolH, &
                                  Ild, IldH, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! Single precision version, storage technique 7 only.
  !</description>
    
  !<input>
    
    ! Number of equations
    integer(I32), intent(IN) :: neq
    
    ! Source matrix
    real(SP), dimension(:), intent(IN) :: FaH
    
    ! Column structure of source matrix
    integer(I32), dimension(:), intent(IN) :: IcolH
    
    ! Row positions of source matrix
    integer(I32), dimension(neq+1), intent(IN) :: IldH
    
    ! Permutation of 1..neq describing the sorting
    integer(I32), dimension(neq), intent(IN) :: Itr1
    
    ! Permutation of 1..neq describing the sorting back
    integer(I32), dimension(neq), intent(IN) :: Itr2

  !</input>
    
  !<output>

    ! Destination matrix
    real(SP), dimension(:), intent(OUT) :: Fa

    ! Column structure of destination matrix
    integer(I32), dimension(:), intent(OUT) :: Icol
    
    ! Row positions of destination matrix
    integer(I32), dimension(neq+1), intent(OUT) :: Ild
    
  !</output>
    
!</subroutine>
    
    !local variables
    integer(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx
    integer(I32), dimension(SORT_BUFFERSIZE) :: Ih1, Ih2
  
    ildIdx = 1
    do i=1, neq
      iidx = Itr1(i)
      Fa(ildIdx) = FaH(IldH(iidx))
      Ild(i) = ildIdx
      Icol(ildIdx) = i
      ildIdx = ildIdx + 1
      
      ih1Idx = IldH(iidx)+1
      ih2Idx = IldH(iidx+1)-1
      isize = ih2Idx - ih1Idx + 1
      do j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      end do
      call sort_crSort(Ih1,Ih2,isize)
      do j=1, isize
        Fa(ildIdx) = FaH(Ih2(j))
        Icol(ildIdx) = Ih1(j) !Achtung: Ist dies so richtig ? Ih1(j)->Ih2(j) ?
        ildIdx = ildIdx + 1
      end do
    end do
    Ild(neq+1) = IldH(neq+1)

  end subroutine sort_matSort_single

!<subroutine>
  subroutine sort_matResort (Icol, IcolH, Ild, IldH, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2:
    ! Sorts the structure of the given matrix back.
    ! Storage technique 7. Only the structure of a given matrix is
    ! resorted, the routine doesn't handle the entries.
  !</description>
    
  !<input>

    ! Number of equations
    integer(I32) , intent(IN) :: neq

    ! Permutation of 1..neq describing how to resort
    integer(I32), dimension(neq), intent(IN) :: Itr1
    
    ! Permutation of 1..neq describing how to sort
    integer(I32), dimension(neq), intent(IN) :: Itr2

  !</input>
    
  !<inputoutput>

    ! Column description of matrix
    integer(I32), dimension(:), intent(INOUT) :: Icol

    ! Row description of matrix
    integer(I32), dimension(neq+1), intent(INOUT) :: Ild
    
    ! Column structure of source matrix -> resorted matrix
    integer(I32), dimension(:), intent(INOUT) :: IcolH
    
    ! Row positions of source matrix -> resorted matrix
    integer(I32), dimension(neq+1), intent(INOUT) :: IldH
    
  !</inputoutput>
!</subroutine>
    
    !local variables
    integer(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx
    integer(I32), dimension(SORT_BUFFERSIZE) :: Ih1, Ih2
    
    ildIdx = 1
    do i=1, neq
      iidx = Itr2(i)
      Ild(i) = ildIdx
      Icol(ildIdx) = i
      ildIdx = ildIdx+1
      ih1Idx = IldH(iidx)+1
      ih2Idx = IldH(iidx+1)-1
      isize = ih2Idx-ih1Idx+1
      do j=ih1Idx, Ih2Idx
        Ih1(j-ih1Idx+1)=Itr1(IcolH(j))
        Ih2(j-ih1Idx+1)=j
      end do
      call sort_crSort(Ih1,Ih2,isize)
      do j=1, isize
        Icol(ildIdx) = Ih1(j)
        ildIdx = ildIdx+1
      end do
    end do
    Ild(neq+1) = IldH(neq+1) 
  
  end subroutine sort_matResort

!<subroutine>
  subroutine sort_crSort (Ih1, Ih2, nsize)
  
  !<description>
    ! Performs bubble sort for the vector Ih1 and does the
    ! same swaps also on vector Ih2
  !</description>
    
  !<input>
    
    ! number of elements to sort
    integer(I32), intent(IN) :: nsize
    
  !</input>
    
  !<inputoutput>
      
    ! column vector
    integer(I32), dimension(SORT_BUFFERSIZE), intent(INOUT) :: Ih1
      
    ! row vector
    integer(I32), dimension(SORT_BUFFERSIZE), intent(INOUT) :: Ih2
      
  !</inputoutput>
!</subroutine>
    
    !local variables
    integer(I32) :: iaux, icomp
    logical :: bmore
    
    bmore = .true.
    ! While there are swaps necessary...
    do while (bmore)
      bmore = .false.
      ! Test for all elements minus the last one
      do icomp=1, nsize-1
        ! If the order is correct, othewise swap entries
        if (Ih1(icomp) .gt. Ih1(icomp+1)) then
          iaux = Ih1(icomp)
          Ih1(icomp) = Ih1(icomp+1)
          Ih1(icomp+1) = iaux
          iaux = Ih2(icomp)
          Ih2(icomp) = Ih2(icomp+1)
          Ih2(icomp+1) = iaux
          bmore = .true.
        endif
      end do  
    end do
  
  end subroutine sort_crSort

end module

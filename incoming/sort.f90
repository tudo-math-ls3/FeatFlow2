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

MODULE sort

  USE fsystem
  
  IMPLICIT NONE
  
  !<constants>
  !<constantblock description="">
  
  ! Auxiliary column array size
  INTEGER, PARAMETER :: SORT_BUFFERSIZE = 100
  
  !</constantblock>
  !</constants>
  
  CONTAINS

!<subroutine>
  SUBROUTINE sort_vecSort_double (Dx, Dd, Itr, neq)
  
  !<description>
    ! Resorts the entries in the vector Dx corresponding to Itr.
    ! The result is written to Dd.
  !</description>
    
  !<input>

    ! Number of equations
    INTEGER(I32) :: neq
  
    ! Source vector to be sorted
    REAL(DP), DIMENSION(neq), INTENT(IN) :: Dx
    
    ! Array with permutation of 1..neq
    INTEGER(I32), DIMENSION(neq), INTENT(IN) :: Itr
  !</input>
    
  !<output>
    ! The resorted vector
    REAL(DP), DIMENSION(neq), INTENT(OUT) :: Dd
  !</output>
    
!</subroutine>
    
    ! local variable
    INTEGER(I32) :: ieq
    
    DO ieq=1, neq
      Dd(ieq) = Dx(Itr(ieq))
    END DO
  
  END SUBROUTINE sort_vecSort_double

!<subroutine>
  SUBROUTINE sort_vecSort_single (Fx, Fd, Itr, neq)
  
  !<description>
    ! Resorts the entries in the vector Fx corresponding to Itr.
    ! The result is written to Fd.
  !</description>
    
  !<input>

    ! Number of equations
    INTEGER(I32) :: neq
  
    ! Source vector to be sorted
    REAL(SP), DIMENSION(neq), INTENT(IN) :: Fx
    
    ! Array with permutation of 1..neq
    INTEGER(I32), DIMENSION(neq), INTENT(IN) :: Itr
  !</input>
    
  !<output>
    ! The resorted vector
    REAL(SP), DIMENSION(neq), INTENT(OUT) :: Fd
  !</output>
    
!</subroutine>
    
    ! local variable
    INTEGER(I32) :: ieq
    
    DO ieq=1, neq
      Fd(ieq) = Fx(Itr(ieq))
    END DO
  
  END SUBROUTINE sort_vecSort_single

!<subroutine>
  SUBROUTINE sort_matSort_double (Da, DaH, Icol, IcolH, &
                                  Ild, IldH, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! Double precision version, storage technique 7 only.
  !</description>
    
  !<input>
    
    ! Number of equations
    INTEGER(I32), INTENT(IN) :: neq
    
    ! Source matrix
    REAL(DP), DIMENSION(:), INTENT(IN) :: DaH
    
    ! Column structure of source matrix
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: IcolH
    
    ! Row positions of source matrix
    INTEGER(I32), DIMENSION(neq+1), INTENT(IN) :: IldH
    
    ! Permutation of 1..neq describing the sorting
    INTEGER(I32), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing the sorting back
    INTEGER(I32), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<output>

    ! Destination matrix
    REAL(DP), DIMENSION(:), INTENT(OUT) :: Da

    ! Column structure of destination matrix
    INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Icol
    
    ! Row positions of destination matrix
    INTEGER(I32), DIMENSION(neq+1), INTENT(OUT) :: Ild
    
  !</output>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx
    INTEGER(I32), DIMENSION(SORT_BUFFERSIZE) :: Ih1, Ih2
  
    ildIdx = 1
    DO i=1, neq
      iidx = Itr1(i)
      Da(ildIdx) = DaH(IldH(iidx))
      Ild(i) = ildIdx
      Icol(ildIdx) = i
      ildIdx = ildIdx + 1
      
      ih1Idx = IldH(iidx)+1
      ih2Idx = IldH(iidx+1)-1
      isize = ih2Idx - ih1Idx + 1
      DO j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      END DO
      CALL sort_crSort(Ih1,Ih2,isize)
      DO j=1, isize
        Da(ildIdx) = DaH(Ih2(j))
        Icol(ildIdx) = Ih1(j) !Achtung: Ist dies so richtig ? Ih1(j)->Ih2(j) ?
        ildIdx = ildIdx + 1
      END DO
    END DO
    Ild(neq+1) = IldH(neq+1)

  END SUBROUTINE sort_matSort_double

!<subroutine>
  SUBROUTINE sort_matSort_single (Fa, FaH, Icol, IcolH, &
                                  Ild, IldH, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! Single precision version, storage technique 7 only.
  !</description>
    
  !<input>
    
    ! Number of equations
    INTEGER(I32), INTENT(IN) :: neq
    
    ! Source matrix
    REAL(SP), DIMENSION(:), INTENT(IN) :: FaH
    
    ! Column structure of source matrix
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: IcolH
    
    ! Row positions of source matrix
    INTEGER(I32), DIMENSION(neq+1), INTENT(IN) :: IldH
    
    ! Permutation of 1..neq describing the sorting
    INTEGER(I32), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing the sorting back
    INTEGER(I32), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<output>

    ! Destination matrix
    REAL(SP), DIMENSION(:), INTENT(OUT) :: Fa

    ! Column structure of destination matrix
    INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Icol
    
    ! Row positions of destination matrix
    INTEGER(I32), DIMENSION(neq+1), INTENT(OUT) :: Ild
    
  !</output>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx
    INTEGER(I32), DIMENSION(SORT_BUFFERSIZE) :: Ih1, Ih2
  
    ildIdx = 1
    DO i=1, neq
      iidx = Itr1(i)
      Fa(ildIdx) = FaH(IldH(iidx))
      Ild(i) = ildIdx
      Icol(ildIdx) = i
      ildIdx = ildIdx + 1
      
      ih1Idx = IldH(iidx)+1
      ih2Idx = IldH(iidx+1)-1
      isize = ih2Idx - ih1Idx + 1
      DO j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      END DO
      CALL sort_crSort(Ih1,Ih2,isize)
      DO j=1, isize
        Fa(ildIdx) = FaH(Ih2(j))
        Icol(ildIdx) = Ih1(j) !Achtung: Ist dies so richtig ? Ih1(j)->Ih2(j) ?
        ildIdx = ildIdx + 1
      END DO
    END DO
    Ild(neq+1) = IldH(neq+1)

  END SUBROUTINE sort_matSort_single

!<subroutine>
  SUBROUTINE sort_matResort (Icol, IcolH, Ild, IldH, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2:
    ! Sorts the structure of the given matrix back.
    ! Storage technique 7. Only the structure of a given matrix is
    ! resorted, the routine doesn't handle the entries.
  !</description>
    
  !<input>

    ! Number of equations
    INTEGER(I32) , INTENT(IN) :: neq

    ! Permutation of 1..neq describing how to resort
    INTEGER(I32), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing how to sort
    INTEGER(I32), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<inputoutput>

    ! Column description of matrix
    INTEGER(I32), DIMENSION(:), INTENT(INOUT) :: Icol

    ! Row description of matrix
    INTEGER(I32), DIMENSION(neq+1), INTENT(INOUT) :: Ild
    
    ! Column structure of source matrix -> resorted matrix
    INTEGER(I32), DIMENSION(:), INTENT(INOUT) :: IcolH
    
    ! Row positions of source matrix -> resorted matrix
    INTEGER(I32), DIMENSION(neq+1), INTENT(INOUT) :: IldH
    
  !</inputoutput>
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx
    INTEGER(I32), DIMENSION(SORT_BUFFERSIZE) :: Ih1, Ih2
    
    ildIdx = 1
    DO i=1, neq
      iidx = Itr2(i)
      Ild(i) = ildIdx
      Icol(ildIdx) = i
      ildIdx = ildIdx+1
      ih1Idx = IldH(iidx)+1
      ih2Idx = IldH(iidx+1)-1
      isize = ih2Idx-ih1Idx+1
      DO j=ih1Idx, Ih2Idx
        Ih1(j-ih1Idx+1)=Itr1(IcolH(j))
        Ih2(j-ih1Idx+1)=j
      END DO
      CALL sort_crSort(Ih1,Ih2,isize)
      DO j=1, isize
        Icol(ildIdx) = Ih1(j)
        ildIdx = ildIdx+1
      END DO
    END DO
    Ild(neq+1) = IldH(neq+1) 
  
  END SUBROUTINE sort_matResort

!<subroutine>
  SUBROUTINE sort_crSort (Ih1, Ih2, nsize)
  
  !<description>
    ! Performs bubble sort for the vector Ih1 and does the
    ! same swaps also on vector Ih2
  !</description>
    
  !<input>
    
    ! number of elements to sort
    INTEGER(I32), INTENT(IN) :: nsize
    
  !</input>
    
  !<inputoutput>
      
    ! column vector
    INTEGER(I32), DIMENSION(SORT_BUFFERSIZE), INTENT(INOUT) :: Ih1
      
    ! row vector
    INTEGER(I32), DIMENSION(SORT_BUFFERSIZE), INTENT(INOUT) :: Ih2
      
  !</inputoutput>
!</subroutine>
    
    !local variables
    INTEGER(I32) :: iaux, icomp
    LOGICAL :: bmore
    
    bmore = .TRUE.
    ! While there are swaps necessary...
    DO WHILE (bmore)
      bmore = .FALSE.
      ! Test for all elements minus the last one
      DO icomp=1, nsize-1
        ! If the order is correct, othewise swap entries
        IF (Ih1(icomp) .GT. Ih1(icomp+1)) THEN
          iaux = Ih1(icomp)
          Ih1(icomp) = Ih1(icomp+1)
          Ih1(icomp+1) = iaux
          iaux = Ih2(icomp)
          Ih2(icomp) = Ih2(icomp+1)
          Ih2(icomp+1) = iaux
          bmore = .TRUE.
        ENDIF
      END DO  
    END DO
  
  END SUBROUTINE sort_crSort

END MODULE

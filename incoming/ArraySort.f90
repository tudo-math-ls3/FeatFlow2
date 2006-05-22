!#########################################################################
!# ***********************************************************************
!# <name> ArraySort </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains all routines and constant necessary to sort an 
!# 2D array by a given subarray.
!#
!# The following routine can be found in this module:
!#
!# 1.) arraySort_sortByIndex
!#     -> Sorts a 2D array(1..size,1..n) by array(index,1..n)
!#
!# </purpose>
!#########################################################################

MODULE ArraySort

  USE fsystem
  
  IMPLICIT NONE
  
  !<constants>
  !<constantblock description="">
  
  ! Cutoff value for Quicksort
  INTEGER, PARAMETER :: SORT_QUICKSORT_CUTOFF = 20
  
  ! Heapsort (bottom up). Dependable allrounder (default)
  INTEGER, PARAMETER :: SORT_HEAP = 0
  
  ! Quicksort (randomized, with cutoff), then Insertionsort
  INTEGER, PARAMETER :: SORT_QUICK = 1
  
  ! Insertionsort (stable, for small or presorted problems)
  INTEGER, PARAMETER :: SORT_INSERT = 2
  
  ! Mergesort (stable)
  INTEGER, PARAMETER :: SORT_MERGE = 3
  
  !</constantblock>
  !</constants>
  
  CONTAINS
!<subroutine>
  SUBROUTINE arraySort_sortByIndex (Ielem, iindex, nnode, nindex, cmethod)
    
    IMPLICIT NONE
    
  !<description>
    ! Sorting routine for a 2D array Ielem(1..nindex,1..nnode) by
    ! the subarray Ielem(iindex,1..nnode) using the method given in
    ! cmethod.
  !</description>
    
  !<input>
    ! Number of nodes
    INTEGER(I32), INTENT(IN) :: nnode
    
    ! Number of indices per node
    INTEGER(I32), INTENT(IN) :: nindex
    
    ! Index number by which to sort the nodes
    INTEGER(I32), INTENT(IN) :: iindex
    
    ! Method to use for sorting (optional). Defaults to Heapsort
    INTEGER(I32), OPTIONAL, INTENT(IN) :: cmethod
  !</input>
        
  !<inputoutput>
    ! 2D array containing the n nodes Ielem(1..nindex,inode)
    INTEGER(I32), DIMENSION(nindex,nnode), INTENT(INOUT) :: Ielem
  !</inputoutput>
!</subroutine>
        
    IF (PRESENT(cmethod)) THEN
      SELECT CASE (cmethod)
        CASE (SORT_HEAP)
          CALL heapSort
        CASE (SORT_QUICK)
          CALL RANDOM_SEED
          CALL quickSort(1,nnode)
          CALL insertSort
        CASE (SORT_INSERT)
          CALL insertSort
        CASE (SORT_MERGE)
          CALL mergeSort(1,nnode)
        CASE DEFAULT
          STOP 'arraySort_sortByIndex: unknown Method:' // sys_i6(cmethod)
      END SELECT
    ELSE
      CALL heapSort
    END IF
    
    CONTAINS

    SUBROUTINE reheap(istart, istop)
      INTEGER(I32), INTENT(IN) :: istart, istop
      INTEGER(I32)::ielem1, ielem2
      INTEGER(I32)::i,j, k, idx

      if(istop.eq.istart) return
      ! Follow patho of bigger children
      i=istart
      j=ishft(i,1)
      ! While there is a child...
      DO WHILE ((j .LE. istop) .AND. (j .GT. 0))
        ! Go one level up
        i = j
        ! Sole child => exit loop
        IF(i .EQ. istop) EXIT
        ! Path of bigger child
        if(Ielem(iindex,i+1) .GT. Ielem(iindex,i)) i=i+1
        j=ishft(i,1)
      END DO
      ! Search for the correct position along the path
      DO WHILE ((i .GT. istart) .AND. &
                (Ielem(iindex,i) .LT. Ielem(iindex,istart)))
        i=ishft(i,-1)
      END DO
      ! Move the nodes
      k=i
      DO idx=1, nindex
        ielem1=Ielem(idx,istart)
        DO WHILE (i .GE. istart)
          ielem2=Ielem(idx,i)
          Ielem(idx,i)=ielem1
          ielem1=ielem2
          i=ishft(i,-1)
        END DO
        i=k
      END DO
    END SUBROUTINE reheap

    SUBROUTINE heapsort()
    
      INTEGER(I32) :: i, t
      
      ! Heap creation phase (Maxheap)
      DO i=ishft(nnode,-1), 1, -1
         CALL reheap(i,nnode)
      END DO
      ! Selection phase
      DO i=nnode, 2, -1
        CALL swapNode(1,i)
        call reheap(1,i-1)
      END DO
    
    END SUBROUTINE heapsort
    
    RECURSIVE SUBROUTINE quicksort (ilower, iupper)
      INTEGER(I32), INTENT(IN) :: ilower, iupper
      
      INTEGER(I32) :: l, u, i, j, t
      REAL(DP) :: r
      
      l = ilower
      u = iupper
      
      DO WHILE ((u-l)>=SORT_QUICKSORT_CUTOFF)
        ! 1.) Choose pivot
        CALL RANDOM_NUMBER(r)
        i = l+FLOOR(r*(u-l+1))
        t = Ielem(iindex,i)
        CALL swapNode(l,i)
        ! 2.) Partition the field and reposition pivot at index j
        i = l
        j = u
        DO
          i = i+1
          DO WHILE ((i .LT. u) .AND. (Ielem(iindex,i) .LT. t))
            i = i+1
          END DO
          DO WHILE (Ielem(iindex,j) .GT. t)
            j = j-1
          END DO
          IF (i .GE. j) EXIT
          CALL swapNode(i,j)
          j = j-1
        END DO
        CALL swapNode(l,j)
        ! 3.) Rekursion (intelligent)
        IF ((j-l) .GT. (u-j)) THEN
          CALL quicksort(j+1,u)
          u = j-1
        ELSE
          CALL quicksort(l,j-1)
          l = j+1
        END IF
      END DO
      
    END SUBROUTINE quicksort
    
    SUBROUTINE insertSort()
      
      INTEGER(I32) :: i, j, k, idx, t
      
      DO i=2, nnode
        t = Ielem(iindex,i)
        j = i-1
        DO WHILE (Ielem(iindex,j)>t)
          j = j-1
          IF (j .EQ. 0) EXIT
        END DO
        j = j+1
        ! Ringshift of Ielem(:,j:i) to the right
        CALL circShiftRight(j,i)
        
      END DO
    
    END SUBROUTINE insertSort
    
    RECURSIVE SUBROUTINE mergesort(ilow, ihigh)
    
      ! A modification of an algorithm by
      ! Jason Harrison, University of
      ! British Columbia.
      ! Porting to F90 by J-P Moreau, Paris.
    
      INTEGER(I32), INTENT(IN) :: ilow, ihigh
      
      INTEGER(I32) :: imid, iend_lo, istart_hi, ilo, ihi
      
      
      ! Nothing to sort
      IF (ilow .GE. ihigh) RETURN
      
      ilo = ilow
      ihi = ihigh
      
      ! Find midpoint of field
      imid = (ilo+ihi)/2
      
      ! Recursive sorting with mergesort
      CALL mergesort(ilow,      imid)
      CALL mergesort(imid+1, ihigh )
      
      ! Merge the two sorted lists in a stable way
      istart_hi = imid+1
      iend_lo = imid
      DO WHILE ((ilo .le. iend_lo) .and. (istart_hi .le. ihi))
        IF (Ielem(iindex,ilo) .le. Ielem(iindex,istart_hi)) THEN
	  ilo = ilo+1
	ELSE
	  CALL circShiftRight(ilo, istart_hi)
	  ilo = ilo+1
	  iend_lo = iend_lo+1
	  istart_hi = istart_hi+1
	ENDIF
      END DO
    
    END SUBROUTINE mergesort

    SUBROUTINE swapNode(i,j)
      ! Swaps node Ielem(:,i) and Ielem(:,j)
      
      INTEGER(I32), INTENT(IN) :: i, j
      
      INTEGER(I32) :: idx, t
      
      DO idx=1, nindex
        t            = Ielem(idx,i)
        Ielem(idx,i) = Ielem(idx,j)
        Ielem(idx,j) = t
      END DO
      
    END SUBROUTINE swapNode

    SUBROUTINE circShiftRight(j,i)
      ! Circular shift of Ielem(:,j:i) one to the right
      INTEGER(I32), INTENT(IN) :: i, j
      
      INTEGER(I32) :: k, t, idx
      
      DO idx=1, nindex
        t = Ielem(idx,i)
        DO k=i, j+1, -1
          Ielem(idx,k) = Ielem(idx,k-1)
        END DO
        Ielem(idx,j) = t
      END DO 
      
    END SUBROUTINE circShiftRight

  END SUBROUTINE arraySort_sortByIndex

END MODULE

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

module ArraySort

  use fsystem
  
  implicit none
  
  !<constants>
  !<constantblock description="">
  
  ! Cutoff value for Quicksort and Mergesort
  integer, parameter :: SORT_CUTOFF = 20
  
  ! Heapsort (bottom up). Dependable allrounder (default)
  integer, parameter :: SORT_HEAP = 0
  
  ! Quicksort (randomized, with cutoff), then Insertionsort
  integer, parameter :: SORT_QUICK = 1
  
  ! Insertionsort (stable, for small or presorted problems)
  integer, parameter :: SORT_INSERT = 2
  
  ! Mergesort (stable)
  integer, parameter :: SORT_MERGE = 3
  
  ! Stable sorting algorithm
  ! Defaults to mergesort, but this can change in the future
  integer, parameter :: SORT_STABLE = 10
    
  !</constantblock>
  !</constants>
  
  contains
!<subroutine>
  subroutine arraySort_sortByIndex (Ielem, iindex, nnode, nindex, cmethod)
    
    implicit none
    
  !<description>
    ! Sorting routine for a 2D array Ielem(1..nindex,1..nnode) by
    ! the subarray Ielem(iindex,1..nnode) using the method given in
    ! cmethod.
  !</description>
    
  !<input>
    ! Number of nodes
    integer(I32), intent(IN) :: nnode
    
    ! Number of indices per node
    integer(I32), intent(IN) :: nindex
    
    ! Index number by which to sort the nodes
    integer(I32), intent(IN) :: iindex
    
    ! Method to use for sorting (optional). Defaults to Heapsort
    integer(I32), optional, intent(IN) :: cmethod
  !</input>
        
  !<inputoutput>
    ! 2D array containing the n nodes Ielem(1..nindex,inode)
    integer(I32), dimension(nindex,nnode), intent(INOUT) :: Ielem
  !</inputoutput>
!</subroutine>
        
    if (present(cmethod)) then
      select case (cmethod)
        case (SORT_HEAP)
          call heapSort
        case (SORT_QUICK)
          call RANDOM_SEED
          call quickSort(1,nnode)
          call insertSort(1,nnode)
        case (SORT_INSERT)
          call insertSort(1,nnode)
        case (SORT_MERGE,SORT_STABLE)
          call mergeSort(1,nnode)
        case DEFAULT
          stop 'arraySort_sortByIndex: unknown Method:' // sys_i6(cmethod)
      end select
    else
      call heapSort
    end if
    
    contains

    subroutine reheap(istart, istop)
      integer(I32), intent(IN) :: istart, istop
      integer(I32)::ielem1, ielem2
      integer(I32)::i,j, k, idx

      if(istop.eq.istart) return
      ! Follow patho of bigger children
      i=istart
      j=ishft(i,1)
      ! While there is a child...
      do while ((j .le. istop) .and. (j .gt. 0))
        ! Go one level up
        i = j
        ! Sole child => exit loop
        if(i .eq. istop) exit
        ! Path of bigger child
        if(Ielem(iindex,i+1) .gt. Ielem(iindex,i)) i=i+1
        j=ishft(i,1)
      end do
      ! Search for the correct position along the path
      do while ((i .gt. istart) .and. &
                (Ielem(iindex,i) .lt. Ielem(iindex,istart)))
        i=ishft(i,-1)
      end do
      ! Move the nodes
      k=i
      do idx=1, nindex
        ielem1=Ielem(idx,istart)
        do while (i .ge. istart)
          ielem2=Ielem(idx,i)
          Ielem(idx,i)=ielem1
          ielem1=ielem2
          i=ishft(i,-1)
        end do
        i=k
      end do
    end subroutine reheap

    subroutine heapsort()
    
      integer(I32) :: i, t
      
      ! Heap creation phase (Maxheap)
      do i=ishft(nnode,-1), 1, -1
         call reheap(i,nnode)
      end do
      ! Selection phase
      do i=nnode, 2, -1
        call swapNode(1,i)
        call reheap(1,i-1)
      end do
    
    end subroutine heapsort
    
    recursive subroutine quicksort (ilower, iupper)
      integer(I32), intent(IN) :: ilower, iupper
      
      integer(I32) :: l, u, i, j, t
      real(DP) :: r
      
      l = ilower
      u = iupper
      
      do while ((u-l) .ge. SORT_CUTOFF)
        ! 1.) Choose pivot
        call random_number(r)
        i = l+floor(r*(u-l+1))
        t = Ielem(iindex,i)
        call swapNode(l,i)
        ! 2.) Partition the field and reposition pivot at index j
        i = l
        j = u
        do
          i = i+1
          do while ((i .lt. u) .and. (Ielem(iindex,i) .lt. t))
            i = i+1
          end do
          do while (Ielem(iindex,j) .gt. t)
            j = j-1
          end do
          if (i .ge. j) exit
          call swapNode(i,j)
          j = j-1
        end do
        call swapNode(l,j)
        ! 3.) Rekursion (intelligent)
        if ((j-l) .gt. (u-j)) then
          call quicksort(j+1,u)
          u = j-1
        else
          call quicksort(l,j-1)
          l = j+1
        end if
      end do
      
    end subroutine quicksort
    
    subroutine insertSort(ilower, iupper)
      integer(I32), intent(IN) :: ilower, iupper
      
      integer(I32) :: i, j, k, idx, t, istop
      
      istop = ilower-1
      do i=ilower+1, iupper
        t = Ielem(iindex,i)
        j = i-1
        do while (Ielem(iindex,j) .gt. t)
          j = j-1
          if (j .eq. istop) exit
        end do
        j = j+1
        ! Ringshift of Ielem(:,j:i) to the right
        call circShiftRight(j,i)
        
      end do
    
    end subroutine insertSort
    
    recursive subroutine mergesort(ilow, ihigh)
    
      ! A modification of an algorithm by
      ! Jason Harrison, University of
      ! British Columbia.
      ! Porting to F90 by J-P Moreau, Paris.
    
      integer(I32), intent(IN) :: ilow, ihigh
      
      integer(I32) :: imid, iend_lo, istart_hi, ilo, ihi
      
      
      ! Nothing to sort
      if (ilow .ge. ihigh) return
      
      if ((ihigh-ilow) .lt. SORT_CUTOFF) then
        call insertsort(ilow, ihigh)
	return
      endif
      
      ilo = ilow
      ihi = ihigh
      
      ! Find midpoint of field
      imid = (ilo+ihi)/2
      
      ! Recursive sorting with mergesort
      call mergesort(ilow,      imid)
      call mergesort(imid+1, ihigh )
      
      ! Merge the two sorted lists in a stable way
      istart_hi = imid+1
      iend_lo = imid
      do while ((ilo .le. iend_lo) .and. (istart_hi .le. ihi))
        if (Ielem(iindex,ilo) .le. Ielem(iindex,istart_hi)) then
	  ilo = ilo+1
	else
	  call circShiftRight(ilo, istart_hi)
	  ilo = ilo+1
	  iend_lo = iend_lo+1
	  istart_hi = istart_hi+1
	endif
      end do
    
    end subroutine mergesort

    subroutine swapNode(i,j)
      ! Swaps node Ielem(:,i) and Ielem(:,j)
      
      integer(I32), intent(IN) :: i, j
      
      integer(I32) :: idx, t
      
      do idx=1, nindex
        t            = Ielem(idx,i)
        Ielem(idx,i) = Ielem(idx,j)
        Ielem(idx,j) = t
      end do
      
    end subroutine swapNode

    subroutine circShiftRight(j,i)
      ! Circular shift of Ielem(:,j:i) one to the right
      integer(I32), intent(IN) :: i, j
      
      integer(I32) :: k, t, idx
      
      do idx=1, nindex
        t = Ielem(idx,i)
        do k=i, j+1, -1
          Ielem(idx,k) = Ielem(idx,k-1)
        end do
        Ielem(idx,j) = t
      end do 
      
    end subroutine circShiftRight

  end subroutine arraySort_sortByIndex

end module

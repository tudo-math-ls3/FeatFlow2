!##############################################################################
!# ****************************************************************************
!# <name> fparser </name>
!# ****************************************************************************
!#
!# <purpose>
!# 
!# This module contains a couple of standard sorting algorithms that can be
!# used to sort integer, real and double-precision memory blocks.
!# as well as arrays.
!#
!# The module contains the following subroutines:
!#
!# 1.) sort_int
!#     -> Sort an integer array
!#
!# 2.) sort_int32
!#     -> Sort an integer(I32) array
!#
!# 3.) sort_sp
!#     -> Sorts a single precision array
!#
!# 4.) sort_dp
!#     -> Sorts a double precision array
!#
!# 5.) arraySort_sortByIndex
!#     -> Sorts an array of integer arrays. One of the positions in the
!#        arrays os chosen as a key for the sorting process.
!#
!# Auxiliary subroutines:
!#
!# 1.) sort_randomSeedOnce
!#     -> Intialise the random number generator
!#
!# </purpose>
!##############################################################################

MODULE sort

  USE fsystem
  USE error
  USE genoutput

  IMPLICIT NONE

!<constants>

! <constantblock description="sort algorithms">
  
  ! heap sort (default); reliable allrounder
  INTEGER, PARAMETER :: SORT_HEAP       = 0

  ! quicksort, then insertsort
  INTEGER, PARAMETER :: SORT_QUICK      = 1

  ! insertsort; for small arrays (stable)
  INTEGER, PARAMETER :: SORT_INSERT     = 2

  ! bubblesort; for very small arrays
  INTEGER, PARAMETER :: SORT_BUBBLE     = 3

  ! mergesort; stable, n*log(n) complexity, log(n)*c stack memory consumption
  ! A modification of an algorithm by Jason Harrison, University of British Columbia.
  ! Ported to F90 by J-P Moreau, Paris.
  ! Further modified and turned into a stable algorithm by Jens F. Acker, University
  ! of Dortmund.
  ! For small field sizes, insertsort is used
  INTEGER, PARAMETER :: SORT_MERGE      = 4

  ! Stable sorting algorithm
  ! Defaults to mergesort, but this can change in the future
  INTEGER, PARAMETER :: SORT_STABLE = 10
    
!</constantblock>

!<constantblock description="General configuration parameters for sorting algorithms">

  ! cutoff value for hybridized sorting
  ! insert best value for your computer here or
  ! leave the original guesstimated value
  INTEGER, PARAMETER :: SORT_CUTOFF     = 35
  
! </constantblock>

!</constants>

CONTAINS

!<subroutine>
  subroutine sort_int(Iarray, csortMethod, Imapping)

    !<description>
    !Sorting routine for integer data type. If more than one vector must be
    !sorted, or if several vectors may be sorted the same way, the mapping
    !array Imapping can be computed. This array must be created outside the routine
    !and must have the same length as Iarray. Then, the other vectors are computed
    !by Vec(Imapping(i)). The array Imapping is an optional argument. The sorting
    !method used is told by the parameter csortMethod. If this optional argument
    !is not used, the routine performs heapsort.
    !</description>

    !<input>

    ! sort algorithm: SORT_HEAP, SORT_QUICK, SORT_INSERT
    integer, optional :: csortMethod

    !</input>

    !<inoutput>

    ! integer array to be sorted
    integer, dimension(:) :: Iarray

    !optional mapping vector (if more than 1 vector may be sorted)
    integer, dimension(:), optional :: Imapping
    !</inoutput>
!</subroutine>

    !if the optional argument csortMethod is present
    if (present(csortMethod)) then
      select case (csortMethod)
      case(SORT_BUBBLE)
        call sort_bubblesort(size(Iarray), Iarray, Imapping)
      case(SORT_HEAP)
        call heapsort(Iarray, Imapping)
      case(SORT_QUICK)
        call sort_randomSeedOnce
        call quicksort(Iarray, Imapping)
        call insertsort(Iarray, Imapping)
      case(SORT_INSERT)
        call insertsort(Iarray, Imapping)
      case(SORT_MERGE)
        call mergesort(Iarray, Imapping)
      case default
        CALL output_line ('Unknown sorting algorithm: '//TRIM(sys_siL(csortMethod,10)), &
            OU_CLASS_ERROR,OU_MODE_STD,'sort_int')        
        STOP
      end select
    else
      call heapsort(Iarray, Imapping)
    endif

  contains

    subroutine sort_bubblesort(nitems, isort, imap)

      integer :: nitems
      integer, dimension(:) :: isort
      integer, dimension(:) :: imap

      integer :: iswp,i
      logical :: bswap

      bswap=.false.

      do

        bswap=.false.

        do i=1,nitems-1

          if (isort(i).gt.isort(i+1)) then

            iswp=isort(i)
            isort(i)=isort(i+1)
            isort(i+1)=iswp

            iswp=imap(i)
            imap(i)=imap(i+1)
            imap(i+1)=iswp

            bswap=.true.
          endif
        enddo

        if (.not.bswap) exit
      enddo

    end subroutine


    !----------------------------------------------------------------


    subroutine reheap(Iarray, istart, istop, Imapping)
      integer::istart,istop
      integer,dimension(1:istop)::Iarray
      integer,dimension(1:istop), optional:: Imapping
      integer::t1,t2
      integer :: i1, i2
      integer::i,j
      !if(istop.eq.start) return ! nothing to correct

      if (present(Imapping)) then

        !trace the path of the bigger children
        i=istart
        j=ishft(i,1)
        do while((j.le.istop).and.(j>0)) !as long as there is a child
          i=j !go one level higher
          if(i.eq.istop) exit !case of just one child
          if(Iarray(i+1).gt.Iarray(i)) i=i+1 !path of bigger child
          j=ishft(i,1)
        enddo

        !search correct position along the path
        do while ((i.gt.istart).and.(Iarray(i).lt.Iarray(istart)))
          i=ishft(i,-1)
        enddo

        !move data
        t1=Iarray(istart)
        i1=Imapping(istart)

        do while (i>=istart)
          t2=Iarray(i)
          Iarray(i)=t1
          t1=t2

          i2=Imapping(i)
          Imapping(i)=i1

          i1=i2
          i=ishft(i,-1)
        enddo
      else
        !trace the path of the bigger children
        i=istart
        j=ishft(i,1)
        do while((j.le.istop).and.(j>0)) !as long as there is a child
          i=j !go one level higher
          if(i.eq.istop) exit !case of just one child
          if(Iarray(i+1).gt.Iarray(i)) i=i+1 !path of bigger child
          j=ishft(i,1)
        enddo

        !search correct position along the path
        do while ((i.gt.istart).and.(Iarray(i).lt.Iarray(istart)))
          i=ishft(i,-1)
        enddo

        !move data
        t1=Iarray(istart)

        do while (i>=istart)
          t2=Iarray(i)
          Iarray(i)=t1
          t1=t2
          i=ishft(i,-1)
        enddo
      endif

    end subroutine reheap


    !----------------------------------------------------------------


    subroutine heapsort(Iarray, Imapping)
      integer,dimension(:)::Iarray
      integer, dimension(:), optional :: Imapping
      integer::t
      integer :: t2
      integer :: i,n
      n = ubound(Iarray,1)
      ! heap creation phase (Maxheap)

      if (present(Imapping)) then
        do i=ishft(n,-1),1,-1
          call reheap(Iarray,i,n,Imapping)
        enddo
        ! selection phase
        do i=n,2,-1
          t2=Imapping(i)
          Imapping(i)=Imapping(1)
          Imapping(1)=t2

          t=Iarray(i)
          Iarray(i)=Iarray(1)
          Iarray(1)=t
          call reheap(Iarray,1,i-1,Imapping)
        enddo
      else
        do i=ishft(n,-1),1,-1
          call reheap(Iarray,i,n)
        enddo
        ! selection phase
        do i=n,2,-1

          t=Iarray(i)
          Iarray(i)=Iarray(1)
          Iarray(1)=t
          call reheap(Iarray,1,i-1)
        enddo
      endif
    end subroutine heapsort


    !----------------------------------------------------------------


    recursive subroutine quicksort(Iarray, Imapping)
      integer,dimension(:)::Iarray
      integer, dimension(:), optional :: Imapping
      integer::t,temp

      integer :: t2, temp2

      integer::l,u,i,j
      real:: r
      l=1
      u=ubound(Iarray,1)

      if (present(Imapping)) then
        do while ((u-l)>=SORT_CUTOFF)
          ! 1.) choice of pivot
          call random_number(r)
          i=l+floor(r*(u-l+1))
          t=Iarray(i)
          Iarray(i)=Iarray(l)
          Iarray(l)=t

          t2=Imapping(i)
          Imapping(i)=Imapping(l)
          Imapping(l)=t2

          ! 2.) partitioning and positioning of the pivot to position j
          i=l
          j=u
          do
            i=i+1
            do while ((i.lt.u).and.(Iarray(i).lt.t))
              i=i+1
            enddo
            do while (Iarray(j).gt.t)
              j=j-1
            enddo
            if (i.ge.j) exit
            temp=Iarray(i)
            temp2=Imapping(i)

            Iarray(i)=Iarray(j)
            Iarray(j)=temp

            Imapping(i)=Imapping(j)
            Imapping(j)=temp2

            j=j-1
          enddo

          Iarray(l)=Iarray(j)
          Iarray(j)=t

          Imapping(l)=Imapping(j)
          Imapping(j)=t2

          ! 3.) recursion (intelligent)
          if((j-l).gt.(u-j)) then
            call quicksort(Iarray(j+1:u), Imapping(j+1:u))
            u=j-1
          else
            call quicksort(Iarray(l:j-1), Imapping(l:j-1))
            l=j+1
          endif
        enddo
      else

        do while ((u-l)>=SORT_CUTOFF)
          ! 1.) choice of pivot
          call random_number(r)
          i=l+floor(r*(u-l+1))
          t=Iarray(i)
          Iarray(i)=Iarray(l)
          Iarray(l)=t

          ! 2.) partitioning and positioning of the pivot to position j
          i=l
          j=u
          do
            i=i+1
            do while ((i.lt.u).and.(Iarray(i).lt.t))
              i=i+1
            enddo
            do while (Iarray(j).gt.t)
              j=j-1
            enddo
            if (i.ge.j) exit
            temp=Iarray(i)
            Iarray(i)=Iarray(j)
            Iarray(j)=temp
            j=j-1
          enddo
          Iarray(l)=Iarray(j)
          Iarray(j)=t

          ! 3.) recursion (intelligent)
          if((j-l).gt.(u-j)) then
            call quicksort(Iarray(j+1:u))
            u=j-1
          else
            call quicksort(Iarray(l:j-1))
            l=j+1
          endif
        enddo

      endif
    end subroutine quicksort


    !----------------------------------------------------------------


    subroutine insertsort(Iarray, Imapping)
      integer,dimension(:)::Iarray
      integer, dimension(:), optional :: Imapping
      integer::t
      integer :: t2
      integer::i,j,k

      if (present(Imapping)) then
        do i=2,ubound(Iarray,1)
          t=Iarray(i)
          t2 = Imapping(i)
          j=i-1

          do while (Iarray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          enddo

!This statement crashed for the sun compiler on linux
!          Iarray(j+2:i)=Iarray(j+1:i-1)
!          Iarray(j+1)=t
          Iarray(j+1:i)=cshift(Iarray(j+1:i),-1)

!This statement crashed for the sun compiler on linux
!          Imapping(j+2:i)=Imapping(j+1:i-1)
!          Imapping(j+1)=t2
          Imapping(j+1:i)=cshift(Imapping(j+1:i),-1)
        enddo
      else
        do i=2,ubound(Iarray,1)
          t=Iarray(i)
          j=i-1
          do while (Iarray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          enddo
!This statement crashed for the sun compiler on linux
!          Iarray(j+2:i)=Iarray(j+1:i-1)
!          Iarray(j+1)=t
          Iarray(j+1:i)=cshift(Iarray(j+1:i),-1)
        enddo
      endif
    end subroutine insertsort


    !----------------------------------------------------------------


    recursive subroutine mergesort(Iarray, Imapping)
      integer,dimension(:)::Iarray
      integer, dimension(:), optional :: Imapping
      integer::imid, ilen
      integer::ilo,iend_lo,istart_hi

      ilen = ubound(Iarray,1)
      ! Nothing to sort
      if (ilen .lt. 2) return
      ! Sort small fields with insertsort.(stable)
      if (ilen .lt. SORT_CUTOFF) then
        call insertsort(Iarray,Imapping)
        return
      endif

      ! Find midpoint of field
      imid = (ilen+1)/2

      ! Recursive sorting with mergesort
      if (present(Imapping)) then
        call mergesort(Iarray(1:imid),Imapping(1:imid))
        call mergesort(Iarray(imid+1:ilen),Imapping(imid+1:ilen))
      else
        call mergesort(Iarray(1:imid))
        call mergesort(Iarray(imid+1:ilen))
      endif

      ! Merge the two sorted lists in a stable way
      ilo = 1
      iend_lo   = imid
      istart_hi = imid+1
      ! Test of presence of Imapping was moved out of loop
      ! for performance reasons
      if (present(Imapping)) then
        do while ((ilo .le. iend_lo) .and. (istart_hi .le. ilen))
          if (Iarray(ilo) .le. Iarray(istart_hi)) then
            ilo = ilo+1
          else
            Iarray(ilo:istart_hi) = cshift(Iarray(ilo:istart_hi),-1)
            Imapping(ilo:istart_hi) = cshift(Imapping(ilo:istart_hi),-1)
            ilo       = ilo+1
            iend_lo   = iend_lo+1
            istart_hi = istart_hi+1
          endif
        enddo
      else
        do while ((ilo .le. iend_lo) .and. (istart_hi .le. ilen))
          if (Iarray(ilo) .le. Iarray(istart_hi)) then
            ilo = ilo+1
          else
            Iarray(ilo:istart_hi) = cshift(Iarray(ilo:istart_hi),-1)
            ilo       = ilo+1
            iend_lo   = iend_lo+1
            istart_hi = istart_hi+1
          endif
        enddo
      endif
    end subroutine mergesort

  end subroutine sort_int


!************************************************************************


!<subroutine>
  subroutine sort_i32(iarray, csortMethod, Imapping)


    !<description>
    !Sorting routine for i32 data type. If more than one vector must be
    !sorted, or if several vectors may be sorted the same way, the mapping
    !array Imapping can be computed. This array must be created outside the routine
    !and must have the same length as Iarray. Then, the other vectors are computed
    !by Vec(Imapping(i)). The array Imapping is an optional argument. The sorting
    !method used is told by the parameter csortMethod. If this optional argument
    !is not used, the routine performs heapsort.
    !</description>


    !<input>

    ! sort algorithm: SORT_HEAP, SORT_QUICK, SORT_INSERT
    integer, optional :: csortMethod

    !</input>

    !<inoutput>

    ! integer array to be sorted
    integer(i32),dimension(:) :: Iarray

    !optional mapping vector (if more than 1 vector may be sorted)
    integer, dimension(:), optional :: Imapping
    !</inoutput>
!</subroutine>

    !if the optional argument csortMethod is present
    if (present(csortMethod)) then
      select case (csortMethod)
      case(SORT_HEAP)
        call heapsort(Iarray, Imapping)
      case(SORT_QUICK)
        call sort_randomSeedOnce
        call quicksort(Iarray, Imapping)
        call insertsort(Iarray, Imapping)
      case(SORT_INSERT)
        call insertsort(Iarray, Imapping)
      case(SORT_MERGE)
        call mergesort(Iarray, Imapping)
      case default
        CALL output_line ('Unknown sorting algorithm: '//TRIM(sys_siL(csortMethod,10)), &
            OU_CLASS_ERROR,OU_MODE_STD,'sort_int')        
      end select
    else
      call heapsort(Iarray, Imapping)
    endif

  contains

    !----------------------------------------------------------------

    subroutine reheap(Iarray, istart, istop, Imapping)
      integer::istart,istop
      integer(i32),dimension(1:istop)::Iarray
      integer,dimension(1:istop), optional:: Imapping
      integer(i32)::t1,t2
      integer :: i1, i2
      integer::i,j
      !if(istop.eq.start) return ! nothing to correct

      if (present(Imapping)) then

        !trace the path of the bigger children
        i=istart
        j=ishft(i,1)
        do while((j.le.istop).and.(j>0)) !as long as there is a child
          i=j !go one level higher
          if(i.eq.istop) exit !case of just one child
          if(Iarray(i+1).gt.Iarray(i)) i=i+1 !path of bigger child
          j=ishft(i,1)
        enddo

        !search correct position along the path
        do while ((i.gt.istart).and.(Iarray(i).lt.Iarray(istart)))
          i=ishft(i,-1)
        enddo

        !move data
        t1=Iarray(istart)
        i1=Imapping(istart)

        do while (i>=istart)
          t2=Iarray(i)
          Iarray(i)=t1
          t1=t2

          i2=Imapping(i)
          Imapping(i)=i1
          i1=i2
          i=ishft(i,-1)
        enddo
      else
        !trace the path of the bigger children
        i=istart
        j=ishft(i,1)
        do while((j.le.istop).and.(j>0)) !as long as there is a child
          i=j !go one level higher
          if(i.eq.istop) exit !case of just one child
          if(Iarray(i+1).gt.Iarray(i)) i=i+1 !path of bigger child
          j=ishft(i,1)
        enddo

        !search correct position along the path
        do while ((i.gt.istart).and.(Iarray(i).lt.Iarray(istart)))
          i=ishft(i,-1)
        enddo

        !move data
        t1=Iarray(istart)

        do while (i>=istart)
          t2=Iarray(i)
          Iarray(i)=t1
          t1=t2
          i=ishft(i,-1)
        enddo
      endif

    end subroutine reheap


    !----------------------------------------------------------------


    subroutine heapsort(Iarray, Imapping)
      integer(i32),dimension(:)::Iarray
      integer, dimension(:), optional :: Imapping
      integer(i32)::t
      integer :: t2
      integer :: i,n
      n=ubound(Iarray,1)
      ! heap creation phase (Maxheap)

      if (present(Imapping)) then
        do i=ishft(n,-1),1,-1
          call reheap(Iarray,i,n,Imapping)
        enddo
        ! selection phase
        do i=n,2,-1
          t2=Imapping(i)
          Imapping(i)=Imapping(1)
          Imapping(1)=t2

          t=Iarray(i)
          Iarray(i)=Iarray(1)
          Iarray(1)=t
          call reheap(Iarray,1,i-1,Imapping)
        enddo
      else
        do i=ishft(n,-1),1,-1
          call reheap(Iarray,i,n)
        enddo
        ! selection phase
        do i=n,2,-1

          t=Iarray(i)
          Iarray(i)=Iarray(1)
          Iarray(1)=t
          call reheap(Iarray,1,i-1)
        enddo
      endif
    end subroutine heapsort


    !----------------------------------------------------------------


    recursive subroutine quicksort(Iarray, Imapping)
      integer(i32),dimension(:)::Iarray
      integer, dimension(:), optional :: Imapping
      integer(i32)::t,temp

      integer :: t2, temp2

      integer::l,u,i,j
      real:: r
      l=1
      u=ubound(Iarray,1)

      if (present(Imapping)) then
        do while ((u-l)>=SORT_CUTOFF)
          ! 1.) choice of pivot
          call random_number(r)
          i=l+floor(r*(u-l+1))
          t=Iarray(i)
          Iarray(i)=Iarray(l)
          Iarray(l)=t

          t2=Imapping(i)
          Imapping(i)=Imapping(l)
          Imapping(l)=t2

          ! 2.) partitioning and positioning of the pivot to position j
          i=l
          j=u
          do
            i=i+1
            do while ((i.lt.u).and.(Iarray(i).lt.t))
              i=i+1
            enddo
            do while (Iarray(j).gt.t)
              j=j-1
            enddo
            if (i.ge.j) exit
            temp=Iarray(i)
            temp2=Imapping(i)

            Iarray(i)=Iarray(j)
            Iarray(j)=temp

            Imapping(i)=Imapping(j)
            Imapping(j)=temp2

            j=j-1
          enddo

          Iarray(l)=Iarray(j)
          Iarray(j)=t

          Imapping(l)=Imapping(j)
          Imapping(j)=t2

          ! 3.) recursion (intelligent)
          if((j-l).gt.(u-j)) then
            call quicksort(Iarray(j+1:u), Imapping(j+1:u))
            u=j-1
          else
            call quicksort(Iarray(l:j-1), Imapping(j+1:u))
            l=j+1
          endif
        enddo
      else

        do while ((u-l)>=SORT_CUTOFF)
          ! 1.) choice of pivot
          call random_number(r)
          i=l+floor(r*(u-l+1))
          t=Iarray(i)
          Iarray(i)=Iarray(l)
          Iarray(l)=t

          ! 2.) partitioning and positioning of the pivot to position j
          i=l
          j=u
          do
            i=i+1
            do while ((i.lt.u).and.(Iarray(i).lt.t))
              i=i+1
            enddo
            do while (Iarray(j).gt.t)
              j=j-1
            enddo
            if (i.ge.j) exit
            temp=Iarray(i)
            Iarray(i)=Iarray(j)
            Iarray(j)=temp
            j= j-1
          enddo
          Iarray(l)=Iarray(j)
          Iarray(j)=t

          ! 3.) recursion (intelligent)
          if((j-l).gt.(u-j)) then
            call quicksort(Iarray(j+1:u))
            u=j-1
          else
            call quicksort(Iarray(l:j-1))
            l=j+1
          endif
        enddo

      endif
    end subroutine quicksort


    !----------------------------------------------------------------


    subroutine insertsort(Iarray, Imapping)
      integer(i32),dimension(:)::Iarray
      integer, dimension(:), optional :: Imapping
      integer(i32)::t
      integer :: t2
      integer::i,j

      if (present(Imapping)) then
        do i=2,ubound(Iarray,1)
          t=Iarray(i)
          t2 = Imapping(i)
          j=i-1
          do while (Iarray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          enddo
          !Iarray(j+2:i)=Iarray(j+1:i-1)
          !Iarray(j+1)=t
          Iarray(j+1:i)=cshift(Iarray(j+1:i),-1)

          !Imapping(j+2:i)=Imapping(j+1:i-1)
          !Imapping(j+1)=t2
          Imapping(j+1:i)=cshift(Imapping(j+1:i),-1)
        enddo
      else
        do i=2,ubound(Iarray,1)
          t=Iarray(i)
          j=i-1
          do while (Iarray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          enddo
          !Iarray(j+2:i)=Iarray(j+1:i-1)
          !Iarray(j+1)=t
          Iarray(j+1:i)=cshift(Iarray(j+1:i),-1)
        enddo
      endif
    end subroutine insertsort


    !----------------------------------------------------------------


    recursive subroutine mergesort(Iarray, Imapping)
      integer(i32),dimension(:)::Iarray
      integer, dimension(:), optional :: Imapping
      integer::imid, ilen
      integer::ilo,iend_lo,istart_hi

      ilen = ubound(Iarray,1)
      ! Nothing to sort
      if (ilen .lt. 2) return
      ! Sort small fields with insertsort.(stable)
      if (ilen .lt. SORT_CUTOFF) then
        call insertsort(Iarray,Imapping)
        return
      endif

      ! Find midpoint of field
      imid = (ilen+1)/2

      ! Recursive sorting with mergesort
      if (present(Imapping)) then
        call mergesort(Iarray(1:imid),Imapping(1:imid))
        call mergesort(Iarray(imid+1:ilen),Imapping(imid+1:ilen))
      else
        call mergesort(Iarray(1:imid))
        call mergesort(Iarray(imid+1:ilen))
      endif

      ! Merge the two sorted lists in a stable way
      ilo = 1
      iend_lo   = imid
      istart_hi = imid+1
      ! Test of presence of Imapping was moved out of loop
      ! for performance reasons
      if (present(Imapping)) then
        do while ((ilo .le. iend_lo) .and. (istart_hi .le. ilen))
          if (Iarray(ilo) .le. Iarray(istart_hi)) then
            ilo = ilo+1
          else
            Iarray(ilo:istart_hi) = cshift(Iarray(ilo:istart_hi),-1)
            Imapping(ilo:istart_hi) = cshift(Imapping(ilo:istart_hi),-1)
            ilo       = ilo+1
            iend_lo   = iend_lo+1
            istart_hi = istart_hi+1
          endif
        enddo
      else
        do while ((ilo .le. iend_lo) .and. (istart_hi .le. ilen))
          if (Iarray(ilo) .le. Iarray(istart_hi)) then
            ilo = ilo+1
          else
            Iarray(ilo:istart_hi) = cshift(Iarray(ilo:istart_hi),-1)
            ilo       = ilo+1
            iend_lo   = iend_lo+1
            istart_hi = istart_hi+1
          endif
        enddo
      endif
    end subroutine mergesort

  end subroutine sort_i32


!************************************************************************


!<subroutine>
  subroutine sort_sp(Darray, csortMethod)

    !<description>
    ! sort routine for single precision
    !</description>

    !<input>

    !sort algorithm: SORT_HEAP, SORT_QUICK, SORT_INSERT
    integer, optional :: csortMethod

    !</input>

    !<inoutput>

    ! singe precision array to be sorted
    real(sp), dimension(:) :: Darray
    !</inoutput>
!</subroutine>

    !if the optional argument csortMethod is present
    if (present(csortMethod)) then
       select case (csortMethod)
       case(SORT_HEAP)
          call heapsort(Darray)
       case(SORT_QUICK)
          call sort_randomSeedOnce
          call quicksort(Darray)
          call insertsort(Darray)
       case(SORT_INSERT)
          call insertsort(Darray)
       case(SORT_MERGE)
          call mergesort(Darray)
       case default
        CALL output_line ('Unknown sorting algorithm: '//TRIM(sys_siL(csortMethod,10)), &
            OU_CLASS_ERROR,OU_MODE_STD,'sort_int')        
      end select
    else
       call heapsort(Darray)
    endif
  contains

    !----------------------------------------------------------------

    subroutine reheap(Darray, istart, istop)
      integer::istart,istop
      real(sp),dimension(1:istop)::Darray
      real(sp)::t1,t2
      integer::i,j
      !if(istop.eq.istart) return ! nothing to correct

      !trace the path of the bigger children
      i=istart
      j=ishft(i,1)
      do while((j.le.istop).and.(j.gt.0)) !as long as there is a child
        i=j !go one level higher
        if(i.eq.istop) exit !case of just one child
        if(Darray(i+1).gt.Darray(i)) i=i+1 !path of bigger child
        j=ishft(i,1)
      enddo

      !search correct position along the path
      do while ((i.gt.istart).and.(Darray(i).lt.Darray(istart)))
        i=ishft(i,-1)
      enddo

      !move data
      t1=Darray(istart)
      do while (i.ge.istart)
        t2=Darray(i)
        Darray(i)=t1
        t1=t2
        i=ishft(i,-1)
      enddo
    end subroutine reheap


    !----------------------------------------------------------------


    subroutine heapsort(Darray)
      real(sp), dimension(:) :: Darray
      real(sp)::t
      integer :: i,n
      n=ubound(Darray,1)
      ! heap creation phase (maxheap)
      do i=ishft(n,-1),1,-1
        call reheap(Darray,i,n)
      enddo
      ! selection phase
      do i=n,2,-1
        t=Darray(i)
        Darray(i)=Darray(1)
        Darray(1)=t
        call reheap(Darray,1,i-1)
      enddo
    end subroutine heapsort


    !----------------------------------------------------------------


    recursive subroutine quicksort(Darray)
      real(sp), dimension(:) :: Darray
      real(sp)::t,temp
      integer::l,u,i,j
      real:: r
      l=1
      u=ubound(Darray,1)
      do while ((u-l)>=SORT_CUTOFF)
        ! 1.) choice of pivot
        call random_number(r)
        i=l+floor(r*(u-l+1))
        t=Darray(i)
        Darray(i)=Darray(l)
        Darray(l)=t
        ! 2.) partitioning and positioning of the pivot to position j
        i=l
        j=u
        do
          i=i+1
          do while ((i.lt.u).and.(Darray(i).lt.t))
            i=i+1
          enddo
          do while (Darray(j).gt.t)
            j=j-1
          enddo
          if (i.ge.j) exit
          temp=Darray(i)
          Darray(i)=Darray(j)
          Darray(j)=temp
          j=j-1
        enddo
        Darray(l)=Darray(j)
        Darray(j)=t
        ! 3.) recursion (intelligent)
        if((j-l).gt.(u-j)) then
           call quicksort(Darray(j+1:u))
           u=j-1
        else
           call quicksort(Darray(l:j-1))
           l=j+1
        endif
      enddo
    end subroutine quicksort


    !----------------------------------------------------------------


    subroutine insertsort(Darray)
      real(sp), dimension(:) :: Darray
      real(sp)::t
      integer::i,j
      do i=2,ubound(Darray,1)
        t=Darray(i)
        j=i-1
        do while (Darray(j)>t)
          j=j-1
          if (j .eq. 0) exit
        enddo
        !Darray(j+2:i)=Darray(j+1:i-1)
        !Darray(j+1)=t
        Darray(j+1:i)=cshift(Darray(j+1:i),-1)
      enddo
    end subroutine insertsort


    !----------------------------------------------------------------


    recursive subroutine mergesort(Darray)
      real(SP), dimension(:) :: Darray
      integer::imid, ilen
      integer::ilo,iend_lo,istart_hi

      ilen = ubound(Darray,1)
      ! Nothing to sort
      if (ilen .lt. 2) return
      ! Sort small fields with insertsort.(stable)
      if (ilen .lt. SORT_CUTOFF) then
        call insertsort(Darray)
        return
      endif

      ! Find midpoint of field
      imid = (ilen+1)/2

      ! Recursive sorting with mergesort
      call mergesort(Darray(1:imid))
      call mergesort(Darray(imid+1:ilen))

      ! Merge the two sorted lists in a stable way
      ilo = 1
      iend_lo   = imid
      istart_hi = imid+1
      do while ((ilo .le. iend_lo) .and. (istart_hi .le. ilen))
        if (Darray(ilo) .le. Darray(istart_hi)) then
          ilo = ilo+1
        else
    Darray(ilo:istart_hi) = cshift(Darray(ilo:istart_hi),-1)
    ilo       = ilo+1
    iend_lo   = iend_lo+1
    istart_hi = istart_hi+1
  endif
      enddo
    end subroutine mergesort

  end subroutine sort_sp


!************************************************************************


!<subroutine>
  subroutine sort_dp(Darray, csortMethod, Imapping)

    !<description>
    !Sorting routine for double precision arrays. If more than one vector must be
    !sorted, or if several vectors may be sorted the same way, the mapping
    !array Imapping can be computed. This array must be created outside the routine
    !and must have the same length as Iarray. Then, the other vectors are computed
    !by Vec(Imapping(i)). The array Imapping is an optional argument. The sorting
    !method used is told by the parameter csortMethod. If this optional argument
    !is not used, the routine performs heapsort.
    !</description>

    !<input>

    !sort algorithm: SORT_HEAP, SORT_QUICK, SORT_INSERT
    integer, optional :: csortMethod

    !</input>

    !<inoutput>

    ! double precision array to be sorted
    real(dp), dimension(:) :: Darray

    !optional mapping vector (if more than 1 vector may be sorted)
    integer, dimension(:), optional :: Imapping

    !</inoutput>
!</subroutine>

    !if the optional argument csortMethod is present
    if (present(csortMethod)) then
      select case (csortMethod)
      case(SORT_HEAP)
        call heapsort(Darray, Imapping)

      case(SORT_QUICK)
        call sort_randomSeedOnce
        call quicksort(Darray, Imapping)
        call insertsort(Darray, Imapping)

      case(SORT_INSERT)

        call insertsort(Darray, Imapping)

      case(SORT_MERGE)
        call mergesort(Darray, Imapping)

      case default
        CALL output_line ('Unknown sorting algorithm: '//TRIM(sys_siL(csortMethod,10)), &
            OU_CLASS_ERROR,OU_MODE_STD,'sort_int')        
      end select
    else

      call heapsort(Darray, Imapping)

    endif

  contains

    !----------------------------------------------------------------

    subroutine reheap(Darray, istart, istop, Imapping)
      integer::istart,istop
      real(dp),dimension(1:istop)::Darray
      integer,dimension(1:istop), optional:: Imapping
      real(dp)::t1,t2
      integer :: i1, i2
      integer::i,j
      !if(istop.eq.istart) return ! nothing to correct

      if (present(Imapping)) then

        !trace the path of the bigger children
        i=istart
        j=ishft(i,1)
        do while((j.le.istop).and.(j.gt.0)) !as long as there is a child
          i=j !go one level higher
          if(i.eq.istop) exit !case of just one child
          if(Darray(i+1).gt.Darray(i)) i=i+1 !path of bigger child
          j=ishft(i,1)
        enddo

        !search correct position along the path
        do while ((i.gt.istart).and.(Darray(i).lt.Darray(istart)))
          i=ishft(i,-1)
        enddo

        !move data
        t1=Darray(istart)
        i1=Imapping(istart)

        do while (i.ge.istart)
          t2=Darray(i)
          Darray(i)=t1
          t1=t2

          i2=Imapping(i)
          Imapping(i)=i1
          i1=i2
          i=ishft(i,-1)
        enddo

      else

        !trace the path of the bigger children
        i=istart
        j=ishft(i,1)
        do while((j.le.istop).and.(j.gt.0)) !as long as there is a child
          i=j !go one level higher
          if(i.eq.istop) exit !case of just one child
          if(Darray(i+1).gt.Darray(i)) i=i+1 !path of bigger child
          j=ishft(i,1)
        enddo

        !search correct position along the path
        do while ((i.gt.istart).and.(Darray(i).lt.Darray(istart)))
          i=ishft(i,-1)
        enddo

        !move data
        t1=Darray(istart)
        do while (i.ge.istart)
          t2=Darray(i)
          Darray(i)=t1
          t1=t2
          i=ishft(i,-1)
        enddo
      endif
    end subroutine reheap


    !----------------------------------------------------------------


    subroutine heapsort(Darray, Imapping)
      real(dp), dimension(:) :: Darray
      integer, dimension(:), optional :: Imapping
      real(dp)::t
      integer :: t2
      integer :: i,n
      n=ubound(Darray,1)
      ! heap creation phase (maxheap)

      if (present(Imapping)) then
        do i=ishft(n,-1),1,-1
          call reheap(Darray,i,n, Imapping)
        enddo
        ! selection phase
        do i=n,2,-1
          t2=Imapping(i)
          Imapping(i)=Imapping(1)
          Imapping(1)=t2

          t=Darray(i)
          Darray(i)=Darray(1)
          Darray(1)=t
          call reheap(Darray,1,i-1, Imapping)
        enddo
      else

        do i=ishft(n,-1),1,-1
          call reheap(Darray,i,n)
        enddo
        ! selection phase
        do i=n,2,-1
          t=Darray(i)
          Darray(i)=Darray(1)
          Darray(1)=t
          call reheap(Darray,1,i-1)
        enddo
      endif
    end subroutine heapsort


    !----------------------------------------------------------------


    recursive subroutine quicksort(Darray, Imapping)
      real(dp), dimension(:) :: Darray
      integer, dimension(:), optional :: Imapping
      real(dp)::t, temp
      integer :: t2, temp2
      integer::l,u,i,j
      real:: r
      l=1
      u=ubound(Darray,1)

      if (present(Imapping)) then
        do while ((u-l)>=SORT_CUTOFF)
          ! 1.) choice of pivot
          call random_number(r)
          i=l+floor(r*(u-l+1))
          t=Darray(i)
          Darray(i)=Darray(l)
          Darray(l)=t

          t2=Imapping(i)
          Imapping(i)=Imapping(l)
          Imapping(l)=t2

          ! 2.) partitioning and positioning of the pivot to position j
          i=l
          j=u
          do
            i=i+1
            do while ((i.lt.u).and.(Darray(i).lt.t))
              i=i+1
            enddo
            do while (Darray(j).gt.t)
              j=j-1
            enddo
            if (i.ge.j) exit
            temp=Darray(i)
            temp2=Imapping(i)

            Darray(i)=Darray(j)
            Darray(j)=temp

            Imapping(i)=Imapping(j)
            Imapping(j)=temp2

            j=j-1
          enddo

          Darray(l)=Darray(j)
          Darray(j)=t

          Imapping(l)=Imapping(j)
          Imapping(j)=t2

          ! 3.) recursion (intelligent)
          if((j-l).gt.(u-j)) then
            call quicksort(Darray(j+1:u), Imapping(j+1:u))
            u=j-1
          else
            call quicksort(Darray(l:j-1), Imapping(l:j-1))
            l=j+1
          endif
        enddo

      else
        do while ((u-l)>=SORT_CUTOFF)
          ! 1.) choice of pivot
          call random_number(r)
          i=l+floor(r*(u-l+1))
          t=Darray(i)
          Darray(i)=Darray(l)
          Darray(l)=t
          ! 2.) partitioning and positioning of the pivot to position j
          i=l
          j=u
          do
            i=i+1
            do while ((i.lt.u).and.(Darray(i).lt.t))
              i=i+1
            enddo
            do while (Darray(j)>t)
              j=j-1
            enddo
            if (i.ge.j) exit
            temp=Darray(i)
            Darray(i)=Darray(j)
            Darray(j)=temp
            j = j-1
          enddo

          Darray(l)=Darray(j)
          Darray(j)=t
          ! 3.) recursion (intelligent)
          if((j-l).gt.(u-j)) then
            call quicksort(Darray(j+1:u))
            u=j-1
          else
            call quicksort(Darray(l:j-1))
            l=j+1
          endif
        enddo
      endif

    end subroutine quicksort


    !----------------------------------------------------------------


    subroutine insertsort(Darray, Imapping)
      real(dp), dimension(:) :: Darray
      integer, dimension(:), optional :: Imapping

      real(dp)::t
      integer :: t2
      integer::i,j

      if (present(Imapping)) then
        do i=2, ubound(Darray,1)
          t=Darray(i)
          t2 = Imapping(i)
          j=i-1
          do while (Darray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          enddo
          !Darray(j+2:i)=Darray(j+1:i-1)
          !Darray(j+1)=t
          Darray(j+1:i)=cshift(Darray(j+1:i),-1)

          !Imapping(j+2:i)=Imapping(j+1:i-1)
          !Imapping(j+1)=t2
          Imapping(j+1:i)=cshift(Imapping(j+1:i),-1)

        enddo
      else
        do i=2,ubound(Darray,1)
          t=Darray(i)
          j=i-1
          do while (Darray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          enddo
          !Darray(j+2:i)=Darray(j+1:i-1)
          !Darray(j+1)=t
          Darray(j+1:i)=cshift(Darray(j+1:i),-1)
        enddo
      endif
    end subroutine insertsort


    !----------------------------------------------------------------


    recursive subroutine mergesort(Darray, Imapping)
      real(DP), dimension(:) :: Darray
      integer, dimension(:), optional :: Imapping
      integer::imid, ilen
      integer::ilo,iend_lo,istart_hi

      ilen = ubound(Darray,1)
      ! Nothing to sort
      if (ilen .lt. 2) return
      ! Sort small fields with insertsort.(stable)
      if (ilen .lt. SORT_CUTOFF) then
        call insertsort(Darray,Imapping)
        return
      endif

      ! Find midpoint of field
      imid = (ilen+1)/2

      ! Recursive sorting with mergesort
      if (present(Imapping)) then
        call mergesort(Darray(1:imid),Imapping(1:imid))
        call mergesort(Darray(imid+1:ilen),Imapping(imid+1:ilen))
      else
        call mergesort(Darray(1:imid))
        call mergesort(Darray(imid+1:ilen))
      endif

      ! Merge the two sorted lists in a stable way
      ilo = 1
      iend_lo   = imid
      istart_hi = imid+1
      ! Test of presence of Imapping was moved out of loop
      ! for performance reasons
      if (present(Imapping)) then
        do while ((ilo .le. iend_lo) .and. (istart_hi .le. ilen))
          if (Darray(ilo) .le. Darray(istart_hi)) then
            ilo = ilo+1
          else
            Darray(ilo:istart_hi) = cshift(Darray(ilo:istart_hi),-1)
            Imapping(ilo:istart_hi) = cshift(Imapping(ilo:istart_hi),-1)
            ilo       = ilo+1
            iend_lo   = iend_lo+1
            istart_hi = istart_hi+1
          endif
        enddo
      else
        do while ((ilo .le. iend_lo) .and. (istart_hi .le. ilen))
          if (Darray(ilo) .le. Darray(istart_hi)) then
            ilo = ilo+1
          else
            Darray(ilo:istart_hi) = cshift(Darray(ilo:istart_hi),-1)
            ilo       = ilo+1
            iend_lo   = iend_lo+1
            istart_hi = istart_hi+1
          endif
        enddo
      endif
    end subroutine mergesort

  end subroutine sort_dp


!************************************************************************


!<subroutine>
  subroutine sort_randomSeedOnce()
    !<description>
    ! This routine ensures that the random seed is only initialised once.
    !
    ! If randomised algorithms are used, a repeated initialisation can
    ! be causing problems, especially if the system time or a fixed start
    ! value are used. This routine enables those algorithms to ensure that
    ! the random number generator is only initialised when needed.
    !</description>
!</subroutine>
    logical, save :: bnotUsed = .TRUE.
!    print *, "START: bnotUsed=", bnotUsed
    if (bnotUsed) then
      call random_seed
      bnotUsed = .FALSE.
    endif
!    print *, "END: bnotUsed=", bnotUsed
  end subroutine sort_randomSeedOnce

!************************************************************************

!<subroutine>
  SUBROUTINE arraySort_sortByIndex (Ielem, iindex, cmethod)
    
    IMPLICIT NONE
    
  !<description>
    ! Sorting routine for a 2D array Ielem(1..nindex,1..nnode) by
    ! the subarray Ielem(iindex,1..nnode) using the method given in
    ! cmethod.
  !</description>
    
  !<input>
    ! Index number of the entry in Ielem that should be used as
    ! a key for sorting.
    INTEGER, INTENT(IN) :: iindex
    
    ! Method to use for sorting (optional). Defaults to Heapsort
    INTEGER(I32), OPTIONAL, INTENT(IN) :: cmethod
  !</input>
    
  !<inputoutput>
    ! 2D array containing the n nodes Ielem(1..nindex,inode)
    INTEGER(I32), DIMENSION(:,:), INTENT(INOUT) :: Ielem
  !</inputoutput>
!</subroutine>

    INTEGER(I32) :: nindex,nnode
        
    nindex = UBOUND(Ielem,1)
    nnode = UBOUND(Ielem,2)
        
    IF (PRESENT(cmethod)) THEN
      SELECT CASE (cmethod)
        CASE (SORT_HEAP)
          CALL heapSort
        CASE (SORT_QUICK)
          CALL RANDOM_SEED
          CALL quickSort(1,nnode)
          CALL insertSort(1,nnode)
        CASE (SORT_INSERT)
          CALL insertSort(1,nnode)
        CASE (SORT_MERGE,SORT_STABLE)
          CALL mergeSort(1,nnode)
        CASE DEFAULT
          CALL output_line('unknown Method:' // sys_i6(cmethod),&
              OU_CLASS_ERROR,OU_MODE_STD,'arraySort_sortByIndex')
          STOP 
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
      
      DO WHILE ((u-l) .GE. SORT_CUTOFF)
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
    
    SUBROUTINE insertSort(ilower, iupper)
      INTEGER(I32), INTENT(IN) :: ilower, iupper
      
      INTEGER(I32) :: i, j, k, idx, t, istop
      
      istop = ilower-1
      DO i=ilower+1, iupper
        t = Ielem(iindex,i)
        j = i-1
        DO WHILE (Ielem(iindex,j) .GT. t)
          j = j-1
          IF (j .EQ. istop) EXIT
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
      
      IF ((ihigh-ilow) .LT. SORT_CUTOFF) THEN
        call insertsort(ilow, ihigh)
        RETURN
      ENDIF
      
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

END MODULE sort


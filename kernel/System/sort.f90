!##############################################################################
!# ****************************************************************************
!# <name> sort </name>
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
!# 2.) sort_i32
!#     -> Sort an integer array
!#
!# 3.) sort_sp
!#     -> Sorts a single precision array
!#
!# 4.) sort_dp
!#     -> Sorts a double precision array
!#
!# 5.) arraySort_sortByIndex_int
!#     -> Sorts an array of integer arrays. One of the positions in the
!#        arrays is chosen as a key for the sorting process.
!#
!# 6.) arraySort_sortByIndex_dp
!#     -> Sorts an array of double precision arrays. One of the positions in
!#        the arrays is chosen as a key for the sorting process.
!#
!# 7.) arraySort_sortByIndex_sp
!#     -> Sorts an array of single precision arrays. One of the positions in
!#        the arrays is chosen as a key for the sorting process.
!#
!# Auxiliary subroutines:
!#
!# 1.) sort_randomSeedOnce
!#     -> Initialise the random number generator
!#
!# </purpose>
!##############################################################################

module sort

!$use omp_lib
  use fsystem
  use error
  use genoutput
  use storage

  implicit none

  private

!<constants>

! <constantblock description="sort algorithms">

  ! heap sort (default); reliable allrounder
  integer, parameter, public :: SORT_HEAP       = 0

  ! quicksort, then insertsort
  integer, parameter, public :: SORT_QUICK      = 1

  ! insertsort; for small arrays (stable)
  integer, parameter, public :: SORT_INSERT     = 2

  ! bubblesort; for very small arrays
  integer, parameter, public :: SORT_BUBBLE     = 3

  ! mergesort; stable, n*log(n) complexity, log(n)*c stack memory consumption
  ! A modification of an algorithm by Jason Harrison, University of British Columbia.
  ! Ported to F90 by J-P Moreau, Paris.
  ! Further modified and turned into a stable algorithm by Jens F. Acker, University
  ! of Dortmund.
  ! For small field sizes, insertsort is used
  integer, parameter, public :: SORT_MERGE      = 4

  ! Stable sorting algorithm
  ! Defaults to mergesort, but this can change in the future
  integer, parameter, public :: SORT_STABLE = 10

!</constantblock>

!<constantblock description="General configuration parameters for sorting algorithms">

  ! cutoff value for hybridised sorting
  ! insert best value for your computer here or
  ! leave the original guesstimated value
  integer, parameter, public :: SORT_CUTOFF     = 35

! </constantblock>

!</constants>

  public :: sort_int
  public :: sort_i32
  public :: sort_sp
  public :: sort_dp
  public :: arraySort_sortByIndex_int
  public :: arraySort_sortByIndex_dp
  public :: arraySort_sortByIndex_sp
  public :: sort_randomSeedOnce

contains

!<subroutine>
  subroutine sort_int(Iarray, csortMethod, Imapping, Itemp)

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

    ! OPTIONAL: Temporary 2D array containing n nodes
    ! Ielem(1..inode). If not specified, the array is
    ! automatically allocated if necessary.
    integer, dimension(:), intent(inout), target, optional :: Itemp

    !</input>

    !<inputoutput>

    ! integer array to be sorted
    integer, dimension(:) :: Iarray

    !optional mapping vector (if more than 1 vector may be sorted)
    integer, dimension(:), optional :: Imapping

    !</inputoutput>
!</subroutine>

    integer, dimension(:), pointer :: p_Itemp,p_Itemp2

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
      case (SORT_MERGE,SORT_STABLE)

        ! We need temporary storage, that algorithm is not in-situ!
        if (present(Imapping)) then
          allocate(p_Itemp2(size(Iarray)))
        else
          ! Dummy
          nullify(p_Itemp2)
        end if

        if (present(Itemp)) then
          p_Itemp => Itemp
          call mergesort(Iarray, p_Itemp, Imapping, p_Itemp2)
        else
          allocate(p_Itemp(size(Iarray)))
          call mergesort(Iarray, p_Itemp, Imapping, p_Itemp2)
          deallocate(p_Itemp)
        end if

        if (present(Imapping)) deallocate(p_Itemp2)

      case default
        call output_line ('Unknown sorting algorithm: '//trim(sys_siL(csortMethod,10)), &
            OU_CLASS_ERROR,OU_MODE_STD,'sort_int')
        call sys_halt()
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
      integer :: istart,istop
      integer, dimension(1:istop) :: Iarray
      integer, dimension(1:istop), optional :: Imapping
      integer :: t1,t2
      integer :: i1, i2
      integer :: i,j
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
      integer, dimension(:) :: Iarray
      integer, dimension(:), optional :: Imapping
      integer :: t
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
      integer, dimension(:) :: Iarray
      integer, dimension(:), optional :: Imapping
      integer :: t,temp
      integer :: t2, temp2
      integer :: l,u,i,j
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
            if (j+1.lt.u) call quicksort(Iarray(j+1:u), Imapping(j+1:u))
            u=j-1
          else
            if(l.lt.j-1) call quicksort(Iarray(l:j-1), Imapping(l:j-1))
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
            if (j+1.lt.u) call quicksort(Iarray(j+1:u))
            u=j-1
          else
            if (l.lt.j-1) call quicksort(Iarray(l:j-1))
            l=j+1
          endif
        enddo

      endif
    end subroutine quicksort


    !----------------------------------------------------------------


    subroutine insertsort(Iarray, Imapping)
      integer, dimension(:) :: Iarray
      integer, dimension(:), optional :: Imapping
      integer :: t
      integer :: t2
      integer :: i,j

      if (present(Imapping)) then
        do i=2,ubound(Iarray,1)
          t=Iarray(i)
          t2 = Imapping(i)
          j=i-1

          do while (Iarray(j)>t)
            j=j-1
            if (j .eq. 0) exit
          enddo

!This statement crashed for the SunStudio compiler on Linux
!          Iarray(j+2:i)=Iarray(j+1:i-1)
!          Iarray(j+1)=t
          Iarray(j+1:i)=cshift(Iarray(j+1:i),-1)

!This statement crashed for the SunStudio compiler on Linux
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
!This statement crashed for the SunStudio compiler on Linux
!          Iarray(j+2:i)=Iarray(j+1:i-1)
!          Iarray(j+1)=t
          Iarray(j+1:i)=cshift(Iarray(j+1:i),-1)
        enddo
      endif
    end subroutine insertsort


    !----------------------------------------------------------------


    recursive subroutine mergesort(Iarray, Itemp, Imapping, ImappingTemp)
      integer, dimension(:) :: Iarray
      integer, dimension(:) :: Itemp
      integer, dimension(:), optional :: Imapping,ImappingTemp
      integer ::imid, ilen
      integer ::ilo,iend_lo,istart_hi
      integer :: idest

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
        call mergesort(Iarray(1:imid),Itemp(1:imid),&
            Imapping(1:imid),ImappingTemp(1:imid))
        call mergesort(Iarray(imid+1:ilen),Itemp(imid+1:ilen),&
            Imapping(imid+1:ilen),ImappingTemp(imid+1:ilen))
      else
        call mergesort(Iarray(1:imid),Itemp(1:imid))
        call mergesort(Iarray(imid+1:ilen),Itemp(imid+1:ilen))
      endif

      ! Merge the two sorted lists in a stable way.
      ! Put the result in the temp array, copy back later to the
      ! original array.
      ilo = 1
      iend_lo   = imid
      istart_hi = imid+1
      idest = 1
      ! Test of presence of Imapping was moved out of loop
      ! for performance reasons
      if (present(Imapping)) then

        do while ((ilo .le. imid) .and. (istart_hi .le. ilen))
          if (Iarray(ilo) .le. Iarray(istart_hi)) then
            Itemp(idest) = Iarray(ilo)
            ImappingTemp(idest) = Imapping(ilo)
            ilo = ilo+1
          else
            Itemp(idest) = Iarray(istart_hi)
            ImappingTemp(idest) = Imapping(istart_hi)
            istart_hi = istart_hi+1
          end if
          idest = idest+1
        end do

        ! Copy the rest of the array. Only one of them may still
        ! contain data!
        do while (ilo .le. imid)
          Itemp(idest) = Iarray(ilo)
          ImappingTemp(idest) = Imapping(ilo)
          ilo = ilo+1
          idest = idest+1
        end do

        do while (istart_hi .le. ilen)
          Itemp(idest) = Iarray(istart_hi)
          ImappingTemp(idest) = Imapping(istart_hi)
          istart_hi = istart_hi+1
          idest = idest+1
        end do

        ! Copy back from the temp array.
        do idest = 1,ilen
          Iarray(idest) = Itemp(idest)
          Imapping(idest) = ImappingTemp(idest)
        end do

      else

        do while ((ilo .le. imid) .and. (istart_hi .le. ilen))
          if (Iarray(ilo) .le. Iarray(istart_hi)) then
            Itemp(idest) = Iarray(ilo)
            ilo = ilo+1
          else
            Itemp(idest) = Iarray(istart_hi)
            istart_hi = istart_hi+1
          end if
          idest = idest+1
        end do

        ! Copy the rest of the array. Only one of them may still
        ! contain data!
        do while (ilo .le. imid)
          Itemp(idest) = Iarray(ilo)
          ilo = ilo+1
          idest = idest+1
        end do

        do while (istart_hi .le. ilen)
          Itemp(idest) = Iarray(istart_hi)
          istart_hi = istart_hi+1
          idest = idest+1
        end do

        ! Copy back from the temp array.
        do idest = 1,ilen
          Iarray(idest) = Itemp(idest)
        end do

      endif
    end subroutine mergesort

  end subroutine sort_int


!************************************************************************


!<subroutine>
  subroutine sort_i32(iarray, csortMethod, Imapping, Itemp)


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

    ! OPTIONAL: Temporary 2D array containing n nodes
    ! Ielem(1..inode). If not specified, the array is
    ! automatically allocated if necessary.
    integer, dimension(:), intent(inout), target, optional :: Itemp

    !</input>

    !<inputoutput>

    ! integer array to be sorted
    integer, dimension(:) :: Iarray

    !optional mapping vector (if more than 1 vector may be sorted)
    integer, dimension(:), optional :: Imapping
    !</inputoutput>
!</subroutine>

    integer :: hhandle,hhandle2
    integer, dimension(:), pointer :: p_Itemp,p_Itemp2

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
      case (SORT_MERGE,SORT_STABLE)

        ! We need temporary storage, that algorithm is not in-situ!
        if (present(Imapping)) then
          call storage_new ('sort_i32', 'temp', size(Iarray), ST_INT, hhandle2, &
              ST_NEWBLOCK_NOINIT)
          call storage_getbase_int (hhandle2,p_Itemp2)
        else
          ! Dummy
          hhandle2 = ST_NOHANDLE
          nullify(p_Itemp2)
        end if

        if (present(Itemp)) then
          p_Itemp => Itemp
          call mergesort(Iarray, p_Itemp, Imapping,p_Itemp2)
        else
          call storage_new ('sort_i32', 'temp', size(Iarray), ST_INT, hhandle, &
              ST_NEWBLOCK_NOINIT)
          ! Probably allocate memory for shifting the Imapping array
          call storage_getbase_int (hhandle,p_Itemp)
          call mergesort(Iarray, p_Itemp, Imapping, p_Itemp2)
          call storage_free (hhandle)
        end if

        if (present(Imapping)) call storage_free (hhandle2)

      case default
        call output_line ('Unknown sorting algorithm: '//trim(sys_siL(csortMethod,10)), &
            OU_CLASS_ERROR,OU_MODE_STD,'sort_int')
      end select
    else
      call heapsort(Iarray, Imapping)
    endif

  contains

    !----------------------------------------------------------------

    subroutine reheap(Iarray, istart, istop, Imapping)
      integer :: istart,istop
      integer, dimension(1:istop) :: Iarray
      integer, dimension(1:istop), optional :: Imapping
      integer :: t1,t2
      integer :: i1, i2
      integer :: i,j
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
      integer, dimension(:) :: Iarray
      integer, dimension(:), optional :: Imapping
      integer :: t
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
      integer, dimension(:) :: Iarray
      integer, dimension(:), optional :: Imapping
      integer :: t,temp
      integer :: t2, temp2
      integer :: l,u,i,j
      real :: r

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
            if (j+1.lt.u) call quicksort(Iarray(j+1:u), Imapping(j+1:u))
            u=j-1
          else
            if (l.lt.j-1) call quicksort(Iarray(l:j-1), Imapping(j+1:u))
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
            if (j+1.lt.u) call quicksort(Iarray(j+1:u))
            u=j-1
          else
            if (l.lt.j-1) call quicksort(Iarray(l:j-1))
            l=j+1
          endif
        enddo

      endif
    end subroutine quicksort


    !----------------------------------------------------------------


    subroutine insertsort(Iarray, Imapping)
      integer, dimension(:) :: Iarray
      integer, dimension(:), optional :: Imapping
      integer :: t
      integer :: t2
      integer :: i,j

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


    recursive subroutine mergesort(Iarray, Itemp, Imapping, ImappingTemp)
      integer, dimension(:) :: Iarray
      integer, dimension(:) :: Itemp
      integer, dimension(:), optional :: Imapping,ImappingTemp
      integer :: imid, ilen
      integer :: ilo,iend_lo,istart_hi
      integer :: idest

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
        call mergesort(Iarray(1:imid),Itemp(1:imid),&
            Imapping(1:imid),ImappingTemp(1:imid))
        call mergesort(Iarray(imid+1:ilen),Itemp(imid+1:ilen),&
            Imapping(imid+1:ilen),ImappingTemp(imid+1:ilen))
      else
        call mergesort(Iarray(1:imid),Itemp(1:imid))
        call mergesort(Iarray(imid+1:ilen),Itemp(imid+1:ilen))
      endif

      ! Merge the two sorted lists in a stable way.
      ! Put the result in the temp array, copy back later to the
      ! original array.
      ilo = 1
      iend_lo   = imid
      istart_hi = imid+1
      idest = 1
      ! Test of presence of Imapping was moved out of loop
      ! for performance reasons
      if (present(Imapping)) then

        do while ((ilo .le. imid) .and. (istart_hi .le. ilen))
          if (Iarray(ilo) .le. Iarray(istart_hi)) then
            Itemp(idest) = Iarray(ilo)
            ImappingTemp(idest) = Imapping(ilo)
            ilo = ilo+1
          else
            Itemp(idest) = Iarray(istart_hi)
            ImappingTemp(idest) = Imapping(istart_hi)
            istart_hi = istart_hi+1
          end if
          idest = idest+1
        end do

        ! Copy the rest of the array. Only one of them may still
        ! contain data!
        do while (ilo .le. imid)
          Itemp(idest) = Iarray(ilo)
          ImappingTemp(idest) = Imapping(ilo)
          ilo = ilo+1
          idest = idest+1
        end do

        do while (istart_hi .le. ilen)
          Itemp(idest) = Iarray(istart_hi)
          ImappingTemp(idest) = Imapping(istart_hi)
          istart_hi = istart_hi+1
          idest = idest+1
        end do

        ! Copy back from the temp array.
        do idest = 1,ilen
          Iarray(idest) = Itemp(idest)
          Imapping(idest) = ImappingTemp(idest)
        end do

      else

        do while ((ilo .le. imid) .and. (istart_hi .le. ilen))
          if (Iarray(ilo) .le. Iarray(istart_hi)) then
            Itemp(idest) = Iarray(ilo)
            ilo = ilo+1
          else
            Itemp(idest) = Iarray(istart_hi)
            istart_hi = istart_hi+1
          end if
          idest = idest+1
        end do

        ! Copy the rest of the array. Only one of them may still
        ! contain data!
        do while (ilo .le. imid)
          Itemp(idest) = Iarray(ilo)
          ilo = ilo+1
          idest = idest+1
        end do

        do while (istart_hi .le. ilen)
          Itemp(idest) = Iarray(istart_hi)
          istart_hi = istart_hi+1
          idest = idest+1
        end do

        ! Copy back from the temp array.
        do idest = 1,ilen
          Iarray(idest) = Itemp(idest)
        end do

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

    !<inputoutput>

    ! singe precision array to be sorted
    real(sp), dimension(:) :: Darray
    !</inputoutput>
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
        call output_line ('Unknown sorting algorithm: '//trim(sys_siL(csortMethod,10)), &
            OU_CLASS_ERROR,OU_MODE_STD,'sort_int')
      end select
    else
       call heapsort(Darray)
    endif
  contains

    !----------------------------------------------------------------

    subroutine reheap(Darray, istart, istop)
      integer :: istart,istop
      real(sp), dimension(1:istop) :: Darray
      real(sp) :: t1,t2
      integer :: i,j
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
      real(sp) :: t
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
      real(sp) :: t,temp
      integer :: l,u,i,j
      real :: r

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
          if (j+1.lt.u) call quicksort(Darray(j+1:u))
           u=j-1
        else
          if (l.lt.j-1) call quicksort(Darray(l:j-1))
           l=j+1
        endif
      enddo
    end subroutine quicksort


    !----------------------------------------------------------------


    subroutine insertsort(Darray)
      real(sp), dimension(:) :: Darray
      real(sp) :: t
      integer :: i,j
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
      integer :: imid, ilen
      integer :: ilo,iend_lo,istart_hi

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
  subroutine sort_dp(Darray, csortMethod, Imapping, Dtemp)

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

    !<inputoutput>

    ! double precision array to be sorted
    real(dp), dimension(:), intent(inout) :: Darray

    !optional mapping vector (if more than 1 vector may be sorted)
    integer, dimension(:), intent(inout), optional :: Imapping

    ! OPTIONAL: Temporary 2D array containing n nodes
    ! Ielem(1..inode). If not specified, the array is
    ! automatically allocated if necessary.
    real(DP), dimension(:), intent(inout), target, optional :: Dtemp

    !</inputoutput>
!</subroutine>

    integer :: hhandle,hhandle2
    real(DP), dimension(:), pointer :: p_Dtemp
    integer, dimension(:), pointer :: p_Itemp2

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

      case (SORT_MERGE,SORT_STABLE)

        ! We need temporary storage, that algorithm is not in-situ!
        if (present(Imapping)) then
          call storage_new ('sort_dp', 'temp', size(Darray), ST_INT, hhandle2, &
              ST_NEWBLOCK_NOINIT)
          call storage_getbase_int (hhandle2,p_Itemp2)
        else
          ! Dummy
          hhandle2 = ST_NOHANDLE
          nullify(p_Itemp2)
        end if

        if (present(Dtemp)) then
          p_Dtemp => Dtemp
          call mergesort(Darray, p_Dtemp, Imapping,p_Itemp2)
        else
          call storage_new ('sort_dp', 'temp', size(Darray), ST_DOUBLE, hhandle, &
              ST_NEWBLOCK_NOINIT)
          ! Probably allocate memory for shifting the Imapping array
          call storage_getbase_double (hhandle,p_Dtemp)
          call mergesort(Darray, p_Dtemp, Imapping, p_Itemp2)
          call storage_free (hhandle)
        end if
        if (present(Imapping)) call storage_free (hhandle2)

      case default
        call output_line ('Unknown sorting algorithm: '//trim(sys_siL(csortMethod,10)), &
            OU_CLASS_ERROR,OU_MODE_STD,'sort_dp')
      end select
    else

      call heapsort(Darray, Imapping)

    endif

  contains

    !----------------------------------------------------------------

    subroutine reheap(Darray, istart, istop, Imapping)
      integer :: istart,istop
      real(dp), dimension(1:istop) :: Darray
      integer, dimension(1:istop), optional :: Imapping
      real(dp) :: t1,t2
      integer :: i1, i2
      integer :: i,j
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
      real(dp) :: t
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
      real(dp) :: t, temp
      integer :: t2, temp2
      integer :: l,u,i,j
      real :: r

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
            if (j+1.lt.u) call quicksort(Darray(j+1:u), Imapping(j+1:u))
            u=j-1
          else
            if (l.lt.j-1) call quicksort(Darray(l:j-1), Imapping(l:j-1))
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
            if (j+1.lt.u) call quicksort(Darray(j+1:u))
            u=j-1
          else
            if (l.lt.j-1) call quicksort(Darray(l:j-1))
            l=j+1
          endif
        enddo
      endif

    end subroutine quicksort


    !----------------------------------------------------------------


    subroutine insertsort(Darray, Imapping)
      real(dp), dimension(:) :: Darray
      integer, dimension(:), optional :: Imapping

      real(dp) :: t
      integer :: t2
      integer :: i,j

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


    recursive subroutine mergesort(Darray, Dtemp, Imapping, ImappingTemp)
      real(DP), dimension(:) :: Darray
      real(DP), dimension(:) :: Dtemp
      integer, dimension(:), optional :: Imapping,ImappingTemp
      integer :: imid, ilen
      integer :: ilo,iend_lo,istart_hi
      integer :: idest

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
        call mergesort(Darray(1:imid),Dtemp(1:imid),&
            Imapping(1:imid),ImappingTemp(1:imid))
        call mergesort(Darray(imid+1:ilen),Dtemp(imid+1:ilen),&
            Imapping(imid+1:ilen),ImappingTemp(imid+1:ilen))
      else
        call mergesort(Darray(1:imid),Dtemp(1:imid))
        call mergesort(Darray(imid+1:ilen),Dtemp(imid+1:ilen))
      endif

      ! Merge the two sorted lists in a stable way.
      ! Put the result in the temp array, copy back later to the
      ! original array.
      ilo = 1
      iend_lo   = imid
      istart_hi = imid+1
      idest = 1
      ! Test of presence of Imapping was moved out of loop
      ! for performance reasons
      if (present(Imapping)) then

        do while ((ilo .le. imid) .and. (istart_hi .le. ilen))
          if (Darray(ilo) .le. Darray(istart_hi)) then
            Dtemp(idest) = Darray(ilo)
            ImappingTemp(idest) = Imapping(ilo)
            ilo = ilo+1
          else
            Dtemp(idest) = Darray(istart_hi)
            ImappingTemp(idest) = Imapping(istart_hi)
            istart_hi = istart_hi+1
          end if
          idest = idest+1
        end do

        ! Copy the rest of the array. Only one of them may still
        ! contain data!
        do while (ilo .le. imid)
          Dtemp(idest) = Darray(ilo)
          ImappingTemp(idest) = Imapping(ilo)
          ilo = ilo+1
          idest = idest+1
        end do

        do while (istart_hi .le. ilen)
          Dtemp(idest) = Darray(istart_hi)
          ImappingTemp(idest) = Imapping(istart_hi)
          istart_hi = istart_hi+1
          idest = idest+1
        end do

        ! Copy back from the temp array.
        do idest = 1,ilen
          Darray(idest) = Dtemp(idest)
          Imapping(idest) = ImappingTemp(idest)
        end do

      else

        do while ((ilo .le. imid) .and. (istart_hi .le. ilen))
          if (Darray(ilo) .le. Darray(istart_hi)) then
            Dtemp(idest) = Darray(ilo)
            ilo = ilo+1
          else
            Dtemp(idest) = Darray(istart_hi)
            istart_hi = istart_hi+1
          end if
          idest = idest+1
        end do

        ! Copy the rest of the array. Only one of them may still
        ! contain data!
        do while (ilo .le. imid)
          Dtemp(idest) = Darray(ilo)
          ilo = ilo+1
          idest = idest+1
        end do

        do while (istart_hi .le. ilen)
          Dtemp(idest) = Darray(istart_hi)
          istart_hi = istart_hi+1
          idest = idest+1
        end do

        ! Copy back from the temp array.
        do idest = 1,ilen
          Darray(idest) = Dtemp(idest)
        end do

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
    logical, save :: bnotUsed = .true.

    if (bnotUsed) then
      call random_seed
      bnotUsed = .false.
    endif

  end subroutine sort_randomSeedOnce

!************************************************************************

!<subroutine>
  subroutine arraySort_sortByIndex_int (Ielem, iindex, cmethod, Itemp)

    implicit none

  !<description>
    ! Sorting routine for a 2D array Ielem(1..nindex,1..nnode) by
    ! the subarray Ielem(iindex,1..nnode) using the method given in
    ! cmethod.
  !</description>

  !<input>
    ! Index number of the entry in Ielem that should be used as
    ! a key for sorting.
    integer, intent(in) :: iindex

    ! Method to use for sorting (optional). Defaults to Heapsort
    integer, optional, intent(in) :: cmethod

    ! OPTIONAL: Temporary 2D array containing n nodes
    ! Ielem(1..nindex,inode). If not specified, the array is
    ! automatically allocated if necessary.
    integer, dimension(:,:), intent(inout), target, optional :: Itemp
  !</input>

  !<inputoutput>
    ! 2D array containing the n nodes Ielem(1..nindex,inode)
    integer, dimension(:,:), intent(inout) :: Ielem
  !</inputoutput>
!</subroutine>

    integer :: nindex,nnode
    integer, dimension(2) :: Isize
    integer :: hhandle
    integer, dimension(:,:), pointer :: p_Itemp

    nindex = ubound(Ielem,1)
    nnode = ubound(Ielem,2)

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
          ! We need temporary storage, that algorithm is not in-situ!
          if (present(Itemp)) then
            p_Itemp => Itemp
            call mergeSort(1,nnode)
          else
            Isize = ubound(Ielem)
            call storage_new ('arraySort_sortByIndex_int', &
                'Itemp', Isize, ST_INT, hhandle, &
                ST_NEWBLOCK_NOINIT)
            call storage_getbase_int2d (hhandle,p_Itemp)
            call mergeSort(1,nnode)
            call storage_free (hhandle)
          end if
        case DEFAULT
          call output_line('unknown Method:' // sys_i6(cmethod),&
              OU_CLASS_ERROR,OU_MODE_STD,'arraySort_sortByIndex')
          call sys_halt()
      end select
    else
      call heapSort
    end if

    contains

    subroutine reheap(istart, istop)
      integer, intent(in) :: istart, istop
      integer :: ielem1, ielem2
      integer :: i,j, k, idx

      if(istop.eq.istart) return
      ! Follow path of bigger children
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

      integer :: i

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
      integer, intent(in) :: ilower, iupper

      integer :: l, u, i, j, t
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
          if (j+1.lt.u) call quicksort(j+1,u)
          u = j-1
        else
          if (l.lt.j-1) call quicksort(l,j-1)
          l = j+1
        end if
      end do

    end subroutine quicksort

    subroutine insertSort(ilower, iupper)
      integer, intent(in) :: ilower, iupper

      integer :: i, j, t, istop

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

      integer, intent(in) :: ilow, ihigh

      integer :: imid, istart_hi, ilo, ihi
      integer :: idx,idest

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
      call mergesort(ilow,   imid)
      call mergesort(imid+1, ihigh )

      ! Merge the two sorted lists in a stable way.
      ! Put the result in the temp array, copy back later to the
      ! original array.
      istart_hi = imid+1
      idest = ilow
      do while ((ilo .le. imid) .and. (istart_hi .le. ihigh))
        if (Ielem(iindex,ilo) .le. Ielem(iindex,istart_hi)) then
          do idx = 1,ubound(Ielem,1)
            p_Itemp(idx,idest) = Ielem(idx,ilo)
          end do
          ilo = ilo+1
        else
          do idx = 1,ubound(Ielem,1)
            p_Itemp(idx,idest) = Ielem(idx,istart_hi)
          end do
          istart_hi = istart_hi+1
        end if
        idest = idest+1
      end do

      ! Copy the rest of the array. Only one of them may still
      ! contain data!
      do while (ilo .le. imid)
        do idx = 1,ubound(Ielem,1)
          p_Itemp(idx,idest) = Ielem(idx,ilo)
        end do
        ilo = ilo+1
        idest = idest+1
      end do

      do while (istart_hi .le. ihigh)
        do idx = 1,ubound(Ielem,1)
          p_Itemp(idx,idest) = Ielem(idx,istart_hi)
        end do
        istart_hi = istart_hi+1
        idest = idest+1
      end do

      ! Copy back from the temp array.
      do idest = ilow,ihigh
        do idx = 1,ubound(Ielem,1)
          Ielem(idx,idest) = p_Itemp(idx,idest)
        end do
      end do

    end subroutine mergesort

    subroutine swapNode(i,j)
      ! Swaps node Ielem(:,i) and Ielem(:,j)

      integer, intent(in) :: i, j

      integer :: idx, t

      do idx=1, nindex
        t            = Ielem(idx,i)
        Ielem(idx,i) = Ielem(idx,j)
        Ielem(idx,j) = t
      end do

    end subroutine swapNode

    subroutine circShiftRight(j,i)
      ! Circular shift of Ielem(:,j:i) one to the right
      integer, intent(in) :: i, j

      integer :: k, t, idx

      do idx=1, nindex
        t = Ielem(idx,i)
        do k=i, j+1, -1
          Ielem(idx,k) = Ielem(idx,k-1)
        end do
        Ielem(idx,j) = t
      end do

    end subroutine circShiftRight

  end subroutine arraySort_sortByIndex_int

!************************************************************************

!<subroutine>
  subroutine arraySort_sortByIndex_dp (Delem, iindex, cmethod, Dtemp)

    implicit none

  !<description>
    ! Sorting routine for a 2D array Delem(1..nindex,1..nnode) by
    ! the subarray Delem(iindex,1..nnode) using the method given in
    ! cmethod.
  !</description>

  !<input>
    ! Index number of the entry in Ielem that should be used as
    ! a key for sorting.
    integer, intent(in) :: iindex

    ! Method to use for sorting (optional). Defaults to Heapsort
    integer, optional, intent(in) :: cmethod

    ! OPTIONAL: Temporary 2D array containing n nodes
    ! Delem(1..nindex,inode). If not specified, the array is
    ! automatically allocated if necessary.
    real(DP), dimension(:,:), intent(inout), target, optional :: Dtemp
  !</input>

  !<inputoutput>
    ! 2D array containing the n nodes Delem(1..nindex,inode)
    real(DP), dimension(:,:), intent(inout) :: Delem
  !</inputoutput>
!</subroutine>

    integer :: nindex,nnode
    integer, dimension(2) :: Isize
    integer :: hhandle
    real(DP), dimension(:,:), pointer :: p_Dtemp

    nindex = ubound(Delem,1)
    nnode = ubound(Delem,2)

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
          ! We need temporary storage, that algorithm is not in-situ!
          if (present(Dtemp)) then
            p_Dtemp => Dtemp
            call mergeSort(1,nnode)
          else
            Isize = ubound(Delem)
            call storage_new ('arraySort_sortByIndex_dp', &
                'Dtemp', Isize, ST_DOUBLE, hhandle, &
                ST_NEWBLOCK_NOINIT)
            call storage_getbase_double2d (hhandle,p_Dtemp)
            call mergeSort(1,nnode)
            call storage_free (hhandle)
          end if

        case DEFAULT
          call output_line('unknown Method:' // sys_i6(cmethod),&
              OU_CLASS_ERROR,OU_MODE_STD,'arraySort_sortByIndex')
          call sys_halt()
      end select
    else
      call heapSort
    end if

    contains

    subroutine reheap(istart, istop)
      integer, intent(in) :: istart, istop
      real(DP) :: delem1, delem2
      integer :: i,j, k, idx

      if(istop.eq.istart) return
      ! Follow path of bigger children
      i=istart
      j=ishft(i,1)
      ! While there is a child...
      do while ((j .le. istop) .and. (j .gt. 0))
        ! Go one level up
        i = j
        ! Sole child => exit loop
        if(i .eq. istop) exit
        ! Path of bigger child
        if(Delem(iindex,i+1) .gt. Delem(iindex,i)) i=i+1
        j=ishft(i,1)
      end do
      ! Search for the correct position along the path
      do while ((i .gt. istart) .and. &
                (Delem(iindex,i) .lt. Delem(iindex,istart)))
        i=ishft(i,-1)
      end do
      ! Move the nodes
      k=i
      do idx=1, nindex
        delem1=Delem(idx,istart)
        do while (i .ge. istart)
          delem2=Delem(idx,i)
          Delem(idx,i)=delem1
          delem1=delem2
          i=ishft(i,-1)
        end do
        i=k
      end do
    end subroutine reheap

    subroutine heapsort()

      integer :: i

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
      integer, intent(in) :: ilower, iupper

      integer :: l, u, i, j
      real(DP) :: t
      real(DP) :: r

      l = ilower
      u = iupper

      do while ((u-l) .ge. SORT_CUTOFF)
        ! 1.) Choose pivot
        call random_number(r)
        i = l+floor(r*(u-l+1))
        t = Delem(iindex,i)
        call swapNode(l,i)
        ! 2.) Partition the field and reposition pivot at index j
        i = l
        j = u
        do
          i = i+1
          do while ((i .lt. u) .and. (Delem(iindex,i) .lt. t))
            i = i+1
          end do
          do while (Delem(iindex,j) .gt. t)
            j = j-1
          end do
          if (i .ge. j) exit
          call swapNode(i,j)
          j = j-1
        end do
        call swapNode(l,j)
        ! 3.) Rekursion (intelligent)
        if ((j-l) .gt. (u-j)) then
          if (j+1.lt.u) call quicksort(j+1,u)
          u = j-1
        else
          if (l.lt.j-1) call quicksort(l,j-1)
          l = j+1
        end if
      end do

    end subroutine quicksort

    subroutine insertSort(ilower, iupper)
      integer, intent(in) :: ilower, iupper

      integer :: i, j, istop
      real(DP) :: t

      istop = ilower-1
      do i=ilower+1, iupper
        t = Delem(iindex,i)
        j = i-1
        do while (Delem(iindex,j) .gt. t)
          j = j-1
          if (j .eq. istop) exit
        end do
        j = j+1
        ! Ringshift of Delem(:,j:i) to the right
        call circShiftRight(j,i)

      end do

    end subroutine insertSort

    recursive subroutine mergesort(ilow, ihigh)

      ! A modification of an algorithm by
      ! Jason Harrison, University of
      ! British Columbia.
      ! Porting to F90 by J-P Moreau, Paris.

      integer, intent(in) :: ilow, ihigh

      integer :: imid, istart_hi, ilo, ihi
      integer :: idx,idest

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
      call mergesort(ilow,   imid)
      call mergesort(imid+1, ihigh )

      ! Merge the two sorted lists in a stable way.
      ! Put the result in the temp array, copy back later to the
      ! original array.
      istart_hi = imid+1
      idest = ilow
      do while ((ilo .le. imid) .and. (istart_hi .le. ihigh))
        if (Delem(iindex,ilo) .le. Delem(iindex,istart_hi)) then
          do idx = 1,ubound(Delem,1)
            p_Dtemp(idx,idest) = Delem(idx,ilo)
          end do
          ilo = ilo+1
        else
          do idx = 1,ubound(Delem,1)
            p_Dtemp(idx,idest) = Delem(idx,istart_hi)
          end do
          istart_hi = istart_hi+1
        end if
        idest = idest+1
      end do

      ! Copy the rest of the array. Only one of them may still
      ! contain data!
      do while (ilo .le. imid)
        do idx = 1,ubound(Delem,1)
          p_Dtemp(idx,idest) = Delem(idx,ilo)
        end do
        ilo = ilo+1
        idest = idest+1
      end do

      do while (istart_hi .le. ihigh)
        do idx = 1,ubound(Delem,1)
          p_Dtemp(idx,idest) = Delem(idx,istart_hi)
        end do
        istart_hi = istart_hi+1
        idest = idest+1
      end do

      ! Copy back from the temp array.
      do idest = ilow,ihigh
        do idx = 1,ubound(Delem,1)
          Delem(idx,idest) = p_Dtemp(idx,idest)
        end do
      end do

    end subroutine mergesort

    subroutine swapNode(i,j)
      ! Swaps node Delem(:,i) and Delem(:,j)

      integer, intent(in) :: i, j

      integer :: idx
      real(DP) :: t

      do idx=1, nindex
        t            = Delem(idx,i)
        Delem(idx,i) = Delem(idx,j)
        Delem(idx,j) = t
      end do

    end subroutine swapNode

    subroutine circShiftRight(j,i)
      ! Circular shift of Delem(:,j:i) one to the right
      integer, intent(in) :: i, j

      integer :: k, idx
      real(DP) :: t

      do idx=1, nindex
        t = Delem(idx,i)
        do k=i, j+1, -1
          Delem(idx,k) = Delem(idx,k-1)
        end do
        Delem(idx,j) = t
      end do

    end subroutine circShiftRight

  end subroutine arraySort_sortByIndex_dp

!************************************************************************

!<subroutine>
  subroutine arraySort_sortByIndex_sp (Selem, iindex, cmethod, Stemp)

    implicit none

  !<description>
    ! Sorting routine for a 2D array Selem(1..nindex,1..nnode) by
    ! the subarray Selem(iindex,1..nnode) using the method given in
    ! cmethod.
  !</description>

  !<input>
    ! Index number of the entry in Ielem that should be used as
    ! a key for sorting.
    integer, intent(in) :: iindex

    ! Method to use for sorting (optional). Defaults to Heapsort
    integer, optional, intent(in) :: cmethod

    ! OPTIONAL: Temporary 2D array containing n nodes
    ! Delem(1..nindex,inode). If not specified, the array is
    ! automatically allocated if necessary.
    real(SP), dimension(:,:), intent(inout), target, optional :: Stemp
  !</input>

  !<inputoutput>
    ! 2D array containing the n nodes Delem(1..nindex,inode)
    real(SP), dimension(:,:), intent(inout) :: Selem
  !</inputoutput>
!</subroutine>

    integer :: nindex,nnode
    integer, dimension(2) :: Isize
    integer :: hhandle
    real(SP), dimension(:,:), pointer :: p_Stemp

    nindex = ubound(Selem,1)
    nnode = ubound(Selem,2)

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
          ! We need temporary storage, that algorithm is not in-situ!
          if (present(Stemp)) then
            p_Stemp => Stemp
            call mergeSort(1,nnode)
          else
            Isize = ubound(Selem)
            call storage_new ('arraySort_sortByIndex_sp', &
                'Stemp', Isize, ST_SINGLE, hhandle, &
                ST_NEWBLOCK_NOINIT)
            call storage_getbase_single2d (hhandle,p_Stemp)
            call mergeSort(1,nnode)
            call storage_free (hhandle)
          end if
        case DEFAULT
          call output_line('unknown Method:' // sys_i6(cmethod),&
              OU_CLASS_ERROR,OU_MODE_STD,'arraySort_sortByIndex')
          call sys_halt()
      end select
    else
      call heapSort
    end if

    contains

    subroutine reheap(istart, istop)
      integer, intent(in) :: istart, istop
      real(SP) :: selem1, selem2
      integer :: i,j, k, idx

      if(istop.eq.istart) return
      ! Follow path of bigger children
      i=istart
      j=ishft(i,1)
      ! While there is a child...
      do while ((j .le. istop) .and. (j .gt. 0))
        ! Go one level up
        i = j
        ! Sole child => exit loop
        if(i .eq. istop) exit
        ! Path of bigger child
        if(Selem(iindex,i+1) .gt. Selem(iindex,i)) i=i+1
        j=ishft(i,1)
      end do
      ! Search for the correct position along the path
      do while ((i .gt. istart) .and. &
                (Selem(iindex,i) .lt. Selem(iindex,istart)))
        i=ishft(i,-1)
      end do
      ! Move the nodes
      k=i
      do idx=1, nindex
        selem1=Selem(idx,istart)
        do while (i .ge. istart)
          selem2=Selem(idx,i)
          Selem(idx,i)=selem1
          selem1=selem2
          i=ishft(i,-1)
        end do
        i=k
      end do
    end subroutine reheap

    subroutine heapsort()

      integer :: i

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
      integer, intent(in) :: ilower, iupper

      integer :: l, u, i, j
      real(SP) :: t
      real(DP) :: r

      l = ilower
      u = iupper

      do while ((u-l) .ge. SORT_CUTOFF)
        ! 1.) Choose pivot
        call random_number(r)
        i = l+floor(r*(u-l+1))
        t = Selem(iindex,i)
        call swapNode(l,i)
        ! 2.) Partition the field and reposition pivot at index j
        i = l
        j = u
        do
          i = i+1
          do while ((i .lt. u) .and. (Selem(iindex,i) .lt. t))
            i = i+1
          end do
          do while (Selem(iindex,j) .gt. t)
            j = j-1
          end do
          if (i .ge. j) exit
          call swapNode(i,j)
          j = j-1
        end do
        call swapNode(l,j)
        ! 3.) Rekursion (intelligent)
        if ((j-l) .gt. (u-j)) then
          if (j+1.lt.u) call quicksort(j+1,u)
          u = j-1
        else
          if (l.lt.j-1) call quicksort(l,j-1)
          l = j+1
        end if
      end do

    end subroutine quicksort

    subroutine insertSort(ilower, iupper)
      integer, intent(in) :: ilower, iupper

      integer :: i, j, istop
      real(SP) :: t

      istop = ilower-1
      do i=ilower+1, iupper
        t = Selem(iindex,i)
        j = i-1
        do while (Selem(iindex,j) .gt. t)
          j = j-1
          if (j .eq. istop) exit
        end do
        j = j+1
        ! Ringshift of Selem(:,j:i) to the right
        call circShiftRight(j,i)

      end do

    end subroutine insertSort

    recursive subroutine mergesort(ilow, ihigh)

      ! A modification of an algorithm by
      ! Jason Harrison, University of
      ! British Columbia.
      ! Porting to F90 by J-P Moreau, Paris.

      integer, intent(in) :: ilow, ihigh

      integer :: imid, istart_hi, ilo, ihi
      integer :: idx,idest

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

      ! Merge the two sorted lists in a stable way.
      ! Put the result in the temp array, copy back later to the
      ! original array.
      istart_hi = imid+1
      idest = ilow
      do while ((ilo .le. imid) .and. (istart_hi .le. ihigh))
        if (Selem(iindex,ilo) .le. Selem(iindex,istart_hi)) then
          do idx = 1,ubound(Selem,1)
            p_Stemp(idx,idest) = Selem(idx,ilo)
          end do
          ilo = ilo+1
        else
          do idx = 1,ubound(Selem,1)
            p_Stemp(idx,idest) = Selem(idx,istart_hi)
          end do
          istart_hi = istart_hi+1
        end if
        idest = idest+1
      end do

      ! Copy the rest of the array. Only one of them may still
      ! contain data!
      do while (ilo .le. imid)
        do idx = 1,ubound(Selem,1)
          p_Stemp(idx,idest) = Selem(idx,ilo)
        end do
        ilo = ilo+1
        idest = idest+1
      end do

      do while (istart_hi .le. ihigh)
        do idx = 1,ubound(Selem,1)
          p_Stemp(idx,idest) = Selem(idx,istart_hi)
        end do
        istart_hi = istart_hi+1
        idest = idest+1
      end do

      ! Copy back from the temp array.
      do idest = ilow,ihigh
        do idx = 1,ubound(Selem,1)
          Selem(idx,idest) = p_Stemp(idx,idest)
        end do
      end do

    end subroutine mergesort

    subroutine swapNode(i,j)
      ! Swaps node Selem(:,i) and Selem(:,j)

      integer, intent(in) :: i, j

      integer :: idx
      real(SP) :: t

      do idx=1, nindex
        t            = Selem(idx,i)
        Selem(idx,i) = Selem(idx,j)
        Selem(idx,j) = t
      end do

    end subroutine swapNode

    subroutine circShiftRight(j,i)
      ! Circular shift of Selem(:,j:i) one to the right
      integer, intent(in) :: i, j

      integer :: k, idx
      real(SP) :: t

      do idx=1, nindex
        t = Selem(idx,i)
        do k=i, j+1, -1
          Selem(idx,k) = Selem(idx,k-1)
        end do
        Selem(idx,j) = t
      end do

    end subroutine circShiftRight

  end subroutine arraySort_sortByIndex_sp

end module sort


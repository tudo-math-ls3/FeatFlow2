program demo
  implicit none

  type t_vectorScalar
     integer :: neq
     real, dimension(:), pointer :: data
  end type t_vectorScalar


  type t_vectorBlock
     integer :: nblocks
     type(t_vectorScalar), dimension(:), pointer :: RvectorBlock
  end type t_vectorBlock


  integer, parameter :: neq = 100

  type(t_vectorBlock) :: Q

  real :: h_i

  integer :: i
  

  ! Initialisierung
  allocate(Q%RvectorBlock(2))
  Q%nblocks = 2

  allocate(Q%RvectorBlock(1)%data(neq))
  allocate(Q%RvectorBlock(2)%data(neq))
  Q%RvectorBlock(1)%neq = neq
  Q%RvectorBlock(2)%neq = neq
  
  
  ! Variable h, Knoten i
  h_i = Q%Rvector(1)%data(i)


  ! Aufräumen
  deallocate(Q%RvectorBlock(1)%data)
  deallocate(Q%RvectorBlock(2)%data)
  deallocate(Q%RvectorBlock)
  




!!$  real, dimension(:), allocatable, target :: myalloc
!!$  real, dimension(:), pointer             :: myptr => NULL()
!!$
!!$  print *,"Allocatable"
!!$  print *, allocated(myalloc)
!!$
!!$  allocate(myalloc(10))
!!$  print *, allocated(myalloc)
!!$
!!$!  deallocate(myalloc)
!!$!  print *, allocated(myalloc)
!!$
!!$
!!$  print *, "Pointer"
!!$  print *, associated(myptr)
!!$
!!$  allocate(myptr(10))
!!$  print *, associated(myptr)
!!$
!!$  deallocate(myptr)
!!$  print *, associated(myptr)
!!$
!!$
!!$  myptr => myalloc
!!$  print *, associated(myptr)



!!$  WRITE(*, FMT='(A)', ADVANCE='NO') "erste Zeile"
!!$  WRITE(*, FMT='(A)') " mit Erweiterung"

end program demo

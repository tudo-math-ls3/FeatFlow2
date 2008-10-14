
!Questions:
! 1) INTEGER(?), REAL(?)=REAL and DOUBLE
! 2) KVERT - its analog
! 3) create a module => include

module quad2tri

  use triangulation

  implicit none  

  ! include statement
  include 'cbasictria.inc'

  contains 

!<subroutine>    
  subroutine quadToTriang_aux1 (nel, Kvert_quad, Kvert_triang)

    !<description>
    ! Purpose: Convert quad mesh to triangular mesh
  
    !This routine creates a triangular KVERT structure from a 
    !quadrilateral KVERT structure.

    !</description>

  !<input>
    !nel    : Number of quad elements
    integer,intent(IN)                    :: nel	

    !Kvert_quad : array [1..4,1..nel] of integer
    !         KVERT structure of the quad mesh
    integer(I32), dimension(:,:), intent(IN)  :: Kvert_quad
  !</input>

  !<output>
    !Kvert_triang : array [1..4,1..2*nel] of integer
    !         KVERT structure of the tri mesh
    integer(I32), dimension(:,:), intent(OUT)  :: Kvert_triang
  !</output>
!</subroutine>
  
    ! local variables
    integer :: i
      
    ! Copy the old KVERT two times, once to the first half and
    ! once to the second half of Kvert_triang:
    call LCP3(Kvert_quad,Kvert_triang,TRIA_MAXNVE2D*nel)
    call LCP3(Kvert_quad,Kvert_triang(1,nel+1),TRIA_MAXNVE2D*nel)

    ! Correct both halfes of Kvert_triang:
    do i=1,nel
      !Set the 4th entry in the first half of Kvert_triang to 0.
      !So the first triangle in each QUAD element consists of the
      !first three vertices 1..3.
      Kvert_triang(4,i) = 0
        
      !The second triangle in each quad consists of vertices 1,3,4.
      !So in the 2nd half, shift entry 3->2 and 4->3 and set the 4th
      !entry to 0.
      Kvert_triang(2,nel+i) = Kvert_triang(3,nel+i)
      Kvert_triang(3,nel+i) = Kvert_triang(4,nel+i)
      Kvert_triang(4,nel+i) = 0
        
      !That's it.
    end do
      
  end subroutine quadToTriang_aux1 

end module quad2tri

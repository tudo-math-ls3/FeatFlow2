
!Questions:
! 1) INTEGER(?), REAL(?)=REAL and DOUBLE
! 2) KVERT - its analog
! 3) create a module => include

MODULE quad2tri

  USE triangulation

  IMPLICIT NONE  

  ! include statement
  INCLUDE 'cbasictria.inc'

  CONTAINS 

!<subroutine>    
  SUBROUTINE quadToTriang_aux1 (nel, Kvert_quad, Kvert_triang)

    !<description>
    ! Purpose: Convert quad mesh to triangular mesh
  
    !This routine creates a triangular KVERT structure from a 
    !quadrilateral KVERT structure.

    !</description>

  !<input>
    !nel    : Number of quad elements
    INTEGER,INTENT(IN)                    :: nel	

    !Kvert_quad : array [1..4,1..nel] of integer
    !         KVERT structure of the quad mesh
    INTEGER(I32), DIMENSION(:,:), INTENT(IN)  :: Kvert_quad
  !</input>

  !<output>
    !Kvert_triang : array [1..4,1..2*nel] of integer
    !         KVERT structure of the tri mesh
    INTEGER(I32), DIMENSION(:,:), INTENT(OUT)  :: Kvert_triang
  !</output>
!</subroutine>
  
    ! local variables
    INTEGER :: i
      
    ! Copy the old KVERT two times, once to the first half and
    ! once to the second half of Kvert_triang:
    CALL LCP3(Kvert_quad,Kvert_triang,TRIA_MAXNVE2D*nel)
    CALL LCP3(Kvert_quad,Kvert_triang(1,nel+1),TRIA_MAXNVE2D*nel)

    ! Correct both halfes of Kvert_triang:
    DO i=1,nel
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
    END DO
      
  END SUBROUTINE quadToTriang_aux1 

END MODULE quad2tri
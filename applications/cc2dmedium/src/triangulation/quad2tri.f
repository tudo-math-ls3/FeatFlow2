***********************************************************************
* Convert quad mesh to triangular mesh;
* Auxiliary routine 1: Convert KVERT
*
* This routine creates a triangular KVERT structure from a 
* quadrilateral KVERT structure.
*
* In: 
*   NEL    : Number of quad elements
*   KVERTQ : array [1..4,1..NEL] of integer
*            KVERT structure of the quad mesh
* Out:
*   KVERTT : array [1..4,1..2*NEL] of integer
*            KVERT structure of the tri mesh
***********************************************************************
      
      SUBROUTINE Q2TAX1 (NEL,KVERTQ,KVERTT)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      
C     parameters

      INTEGER NEL,KVERTQ(NNVE,*),KVERTT(NNVE,*)
      
C     local variables

      INTEGER I
      
C     Copy the old KVERT two times, once to the first half and
C     once to the second half of KVERTT:

      CALL LCP3(KVERTQ,KVERTT,4*NEL)
      CALL LCP3(KVERTQ,KVERTT(1,NEL+1),4*NEL)

C     Correct both halfes of KVERTT:
      
      DO I=1,NEL

C       Set the 4th entry in the first half of KVERTT to 0.
C       So the first triangle in each QUAD element consists of the
C       first three vertices 1..3.
      
        KVERTT(4,I) = 0
        
C       The second triangle in each quad consists of vertices 1,3,4.
C       So in the 2nd half, shift entry 3->2 and 4->3 and set the 4th
C       entry to 0.

        KVERTT(2,NEL+I) = KVERTT(3,NEL+I)
        KVERTT(3,NEL+I) = KVERTT(4,NEL+I)
        KVERTT(4,NEL+I) = 0
        
C       That's it.
        
      END DO
      
      END 
      
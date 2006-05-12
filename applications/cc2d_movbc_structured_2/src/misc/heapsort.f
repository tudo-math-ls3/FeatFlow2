************************************************************************
* This file implements heapsort strategies for arrays.
************************************************************************

************************************************************************
* Heapsort for integer arrays
*
* This routine accepts a 2D array, i.e. an array of arrays with data
* of the same type. Each sub-array is treated as a node. The routine
* sorts these nodes in ascending order according to a criterium using
* standard heapsort. An index identifies the guiding element in
* each array-node heapsort uses for determining the order.
*
* In:
*   N      - Number of sub-arrays / array-nodes
*   SIZE   - Size of each node = number of elements in each sub-array
*   ELMS   - array [1..SIZE,1..N] of integer
*            The elements that are to be sorted
*   IDX    - integer, 1..SIZE. Index of the key integer inside of a
*            subarray that serves for determining the order
*
* Out:
*   ELMS   - Resorted element arrays.
*
************************************************************************

* Remark: Although "bottom-up" heap-sort would be a little bit cheaper,
* the standard heap-sort is easier to implement, so we use that one...

      SUBROUTINE HSRTIN (N,SIZE,ELMS,IDX)
      
      IMPLICIT NONE
      
      INTEGER N,SIZE,IDX
      INTEGER ELMS(SIZE,N)
     
C local variables

      INTEGER I,J
      INTEGER VTMP

C Heap setup phase.
C Loop through all sub-trees and ensure the inverse heap property:
C largest element on the root. Of course we don't have to take care
C about the leaves...

      DO I=N/2,1,-1
        CALL HSRHIN (N,SIZE,I,ELMS,IDX)
      END DO

C Now we have to move the largest elements to the bottom. Loop through
C all elements:

      DO I=N,2,-1
      
C Exchange the current root I with the topmost element, which is
C the currently largest:
      
        DO J=1,SIZE
          VTMP = ELMS(J,I)
          ELMS(J,I) = ELMS(J,1)
          ELMS(J,1) = VTMP
        END DO
      
C The main root lacks the inversed heap property; reheap to
C reensure it, Root=topmost root, only I-1 elements in the
C tree (as the elements I..N are already sorted):

        CALL HSRHIN (I-1,SIZE,1,ELMS,IDX)
      
      END DO

      END

************************************************************************
* Bottom-up heapsort, reheap subroutine
*
* This subroutine implements the reheap phase for the heapsort
* algorithm. It assures that a subtree of the whole tree fulfills
* the inverse heap property (largest element on the root).
*
* In:
*   N      - Number of sub-arrays / array-nodes
*   SIZE   - Size of each node = number of elements in each sub-array
*   ROOT   - Index of the root inside of the element tree
*   ELMS   - array [1..SIZE,1..N] of integer
*            The elements that are to be sorted
*   IDX    - integer, 1..SIZE. Index of the key integer inside of a
*            subarray that serves for determining the order
*
* Out:
*   ELMS   - Resorted element arrays.
************************************************************************

      SUBROUTINE HSRHIN (N,SIZE,ROOT,ELMS,IDX)
      
      IMPLICIT NONE
      
      INTEGER N,SIZE,ROOT,IDX
      INTEGER ELMS(SIZE,N)
     
C local variables

      INTEGER CNDE,RCH,LCH,XC,I
      INTEGER VAL,VTMP
      
C We use the "inversed heap property" (largest element on the root
C instead of smallest element), since we want to sort in ascending
C order.
C
C We assume that the two subtrees below our root fulfill the inversed 
C heap property, only the root might violate this.
C
C We start with the node on the root:

      CNDE = ROOT
      
C Its value is taken from the sub-array

      VAL = ELMS(IDX,CNDE)
      
C Crippled Do-While-Loop
      
10    CONTINUE
      
C CNDE*2 gives the left child, CNDE*2+1 the right child.
C XC denotes a found child that is larger and has to be exchanged
C with our current node.

        XC = 0
        LCH=2*CNDE
        RCH=2*CNDE+1
      
C One or two children?

        IF (RCH.LE.N) THEN

C Is the right child larger

          IF (ELMS(IDX,RCH).GT.VAL) THEN
            
            XC = RCH
        
C Left child even larger?

            IF (ELMS(IDX,LCH).GT.ELMS(IDX,XC)) XC = LCH
          
          ELSE
          
C Left child larger?

            IF (ELMS(IDX,LCH).GT.VAL) XC = LCH

          END IF
        
        ELSE IF (LCH.LE.N) THEN
          IF (ELMS(IDX,LCH).GT.VAL) THEN
        
C           Only one child - and it is larger
          
            XC = LCH
          END IF
          
        END IF
      
C If a larger child was found (XC<>0), exchange that with our current
C node and continue the search in the lower sub-tree.
C The search ends if no larger element was found or the
C bottom of the tree is reached.

        IF (XC.NE.0) THEN
          DO I=1,SIZE
            VTMP = ELMS(I,XC)
            ELMS(I,XC) = ELMS(I,CNDE)
            ELMS(I,CNDE) = VTMP
          END DO
          CNDE = XC
        END IF
      IF (XC.NE.0) GOTO 10
      
      END
      
!#########################################################################
!# ***********************************************************************
!# <name> cmsort </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains all routines necessary for Cuthill McKee
!# matrix renumbering.
!#
!# For storing matrices, we use the CSR format, realised as ''Format-7''
!# matrices. Format-7 is a variant of the CSR-format used in old FEAT
!# programs, where the diagonal element is stored first in each line.
!#
!# The following routines can be found in this module:
!#
!# 1.) cmsort_calcColumnNumbering
!#     -> Calculate column numbering
!#
!# 2.) cmsort_calcPermutation
!#     -> Calculate permutation
!#
!# </purpose>
!#########################################################################

MODULE cmsort

  USE fsystem

  IMPLICIT NONE
  
  CONTAINS

!<subroutine>
  SUBROUTINE cmsort_calcColumnNumbering (Ild, Icol, Icon, Ideg, neq, ndeg)
  
    !<description>
    ! Purpose: Cuthill McKee matrix renumbering
    !          Calculate column numbering
  
    !The algorithm of Cuthill McKee interprets the system matrix as
    !adjacense matrix. Every Row denotes a note in the corresponing graph.
    !In the first step this function sorts the columns in every row
    !for increasing degree of the nodes in the matrix.
    !The matrix must have symmetric structure!
    !(For FE matrices this is always the case...)

    !Matrix storage technique 7 only!
    !</description>
  
  !<input>
  
    ! Number of equations
    INTEGER(I32),INTENT(IN)                    :: neq

    ! Maximum number of entries != 0 in every row of the matrix
    INTEGER(I32), INTENT(IN)                   :: ndeg    
   
    ! Row description of matrix
    INTEGER(I32), DIMENSION(neq+1), INTENT(IN) :: Ild
  
    ! Column description of matrix
    INTEGER(I32), DIMENSION(:), INTENT(IN)     :: Icol
    
    ! Auxiliary vector; must be at least as long as the
    ! maximum number of entries != 0 in every row of the matrix
    ! INTENT(IN) would create here compiler errors
    INTEGER(I32), DIMENSION(ndeg)              :: Ideg
  
  !</input>
  
  !<inputoutput>
  
    ! Auxiliary vector; the column numbers of KCOL are assigned to this in
    ! the order of increasing degree. When calling the routine the user
    ! must copy the content of KCOL to this! These values are then
    ! resorted.
    INTEGER(I32), DIMENSION(:), INTENT(INOUT)  :: Icon
  !</inputoutput>

!</subroutine>

    ! local variables
    INTEGER(I32) :: ieq, idegIdx, ildIdx, idegIdx1, idegIdx2
    INTEGER(I32) :: idegMin, iidxMin
  
    ! Clear auxiliary vector
    DO idegIdx=1, ndeg
      Ideg(idegIdx) = 0
    END DO
    
    ! Later in the algorithm the column numbers of the non-diagonal
    ! entries are changed. The numbers of the diagonal entries are not
    ! touched. Therefore we copy the column numbers of the diagonal
    ! entries directly.
    DO ieq=1, neq
      Icon(Ild(ieq)) = Icol(Ild(ieq))
    END DO
    
    ! Loop about all rows in the matrix: ieq = current row.
    ! Every row corresponds to an entry in the solution vector = a node
    ! in the graph of the matrix.
    DO ieq=1, neq
    
      ! Copy the column numbers of the non-diagonal entries in the current
      ! column to the auxiliary vector Ideg. Ideg contains always the
      ! following nodes of the current node ieq, which are not processed yet.
      ! The entries of the vector are set to 0 one after the other in the
      ! order of the degree of the following nodes.
      DO ildIdx=Ild(ieq)+1, Ild(ieq+1)-1
        Ideg(ildIdx-Ild(ieq)) = Icol(ildIdx)
      END DO
      
      ! Loop about every column in the current row. The entries in the
      ! row (=column numbers) of the matrix represent the node numbers of
      ! those nodes in the graph of the matrix, which are connected to the
      ! current node ieq.
      
      ! The algorithm now has to sort these nodes for increasing degree.
      ! For this purpose the node with the smallest degree is searched and
      ! deleted from Ideg, then the node with the 2nd smallest degree
      ! and so on.
      ! We don't have many nodes, so we use a simple O(n^2) algorithm here
      ! which simply loops over all nodes.
      ! It's only necessary to search in the non-diagonal entries,
      ! because we want to sort the adjacent nodes of the current one
      ! for increasing degree.
      DO idegIdx1=1, Ild(ieq+1)-Ild(ieq)-1
      
        ! iidxMin will receive the index in the Ideg-array of the node with
	! the smallest degree. We start with node 1 in the current column:
	iidxMin=1
	
	! the variable idegMin always contains the degree of the node, which
	! is described by the index iidxMin
	
	! If Ideg(iidxMin)=0, this node has already been processed. In this
	! case idegMin is set to neq - the largest upper bound which implies
	! that iidxMin is later on surely replaced by a node with a smaller
	! degree.
	
	! Otherwise idegMin receives the degree of the node described by
	! iidxMin, which is calculated by the number of elements in Icol,
	! i.e. the difference of the indices in the Ild array.
	IF (Ideg(iidxMin) .EQ. 0) THEN
	  idegMin = neq
	ELSE
	  idegMin = Ild(Ideg(iidxMin)+1)-Ild(Ideg(iidxMin))-1
	ENDIF
	
	! Compare iidxMin with every node in that line to find that with
	! minimum degree:
	DO idegIdx2=1, Ild(ieq+1)-Ild(ieq)-1

          ! If Ideg(idegIdx2)=0, idegIdx2 has already been processed. Set
	  ! idegIdx=neq ; here a lower bound to prevent the current node
	  ! iidxMin from being replaced by idegIdx2.
	  
	  ! Otherwise set idegIdx to the degree of the node with index
	  ! idegIdx2
	  IF (Ideg(idegIdx2) .EQ. 0) THEN
	    idegIdx = neq
	  ELSE
	    idegIdx = Ild(Ideg(idegIdx2)+1)-Ild(Ideg(idegIdx2))-1
          ENDIF
         
	  ! If now idegIdx=grad(iidxMin) < grad(idegIdx2)=idegMin set
	  ! iidxMin to the new node with smaller degree.
	  IF (idegIdx .LT. idegMin) THEN
	    iidxMin=idegIdx2
	    idegMin=idegIdx
	  ENDIF

	END DO
	
	! At this point, iidxMin contains the index in Ideg of that node with
	! minimal degree, which has not yet been processed.
	! Ideg(iidxMin) contains the column number of that column of the matrix,
	! which corresponds with the node with the next higher degree.
	
	! This value is written to Icon. If this subroutine is finished,
	! Icon therefore contains the column numbers of the Icol array -
	! sorted for increasing degree of the nodes in the graph of the matrix.
	Icon(Ild(ieq)+idegIdx1) = Ideg(iidxMin)
	
	! Set the Ideg-value of the node with minimum degree to 0 to indicate
	! that it has been precessed. This node is then ignored on the next
	! searching process for nodes with minimum degree.
	Ideg(iidxMin) = 0
	
      END DO
      
      ! Clear auxiliary vector; only some entries were used. This is only for
      ! reasons of safetyness, as if the upper loops are processed correctly,
      ! (no nodes were forgotten), all Ideg-arrays should already be 0.
      DO idegIdx=1, Ild(ieq+1)-Ild(ieq)
        Ideg(idegIdx) = 0
      END DO
      
    ! Process next line:  
    END DO 
  
  END SUBROUTINE cmsort_calcColumnNumbering

!<subroutine>
  SUBROUTINE cmsort_calcPermutation (Ild, Icon, neq, Itr1, Itr2)
  
    !<description>
    !Purpose: Cuthill McKee matrix renumbering
    !         Calculate permutation

    ! Uses the Ild/Icon vectors of cmsort_calcColumnNumbering to calc-
    ! ulate permutation vectors Itr1 and Itr2; these describe how a
    ! vector must be sorted/sorted back.

    ! Itr1 must be initialized with 0 on call of this subroutine!
    ! Itr2 should be initialized with 0 on call of this subroutine; ther-
    ! wise only the entries != 0 are calculated!

    ! Matrix storage technique 7 only!
    !</description>
    
    !<input>
    
    ! Number of equations
    INTEGER(I32), INTENT(IN)                   :: neq
    
    ! Row description of matrix
    INTEGER(I32), DIMENSION(neq+1), INTENT(IN) :: Ild
    
    ! Auxiliary vector, calculated with cmsort_calcColumnNumbering
    INTEGER(I32), DIMENSION(:), INTENT(IN)     :: Icon
    !</input>
    
    !<output>
    
    ! The permutation vector that describes how the solution vector
    ! has to be restored
    INTEGER(I32), DIMENSION(neq), INTENT(OUT) :: Itr1
    
    ! Describes how a resorted vector can be sorted back 
    INTEGER(I32), DIMENSION(neq), INTENT(OUT) :: Itr2

    !</output>    

  !</subroutine>
    
    ! local variables
    INTEGER(I32) :: ineqIdx, icount, ildIdx, icolIdx
    
!    DO icount=1, neq
!      IF (Itr1(icount) .EQ. 0) Itr1(icount) = icount
!    END DO
    ! Commented out, old FEAT-Code. The current code is using the fact that
    ! the Itr1-vector is filled with 0 on call to this routine.

    ! Itr1 should receive a permutation of the numbers 1..neq for
    ! resorting the entries, Itr2 should receive the inverse permutation.

    ! The Cuthill-McCee algorithm processes levelwise. One node is picked
    ! as the starting node and receives the first number. The nodes
    ! in the neighbourhood receive the next numbers in increasing order
    ! of the degrees of these nodes. Then the neighbourhood of these
    ! nodes is analyzed for increasing degree to give the next node
    ! numbers, and so on.
    !
    ! For this purpose the permutation vectors Itr1 and Itr2 are build
    ! subsequently, depending on which nodes have already been processed.        

    ! We fix the first node (=first line of the matrix), i.e. we set
    ! Itr1(1)=1. This of course implies that the inverse permutation
    ! does not have to process node 1, i.e. Itr2(1)=1. But this will
    ! be computed by the algorithm later.
    
!    Itr1(1)=1
!    Itr2(1)=1
    ! Old FEAT-code, commented out.
    ! The current code uses the fact, that Itr1 is initialized with 0 on
    ! call to this routine. In in the processing a 0 is found in Itr1, a
    ! new matrix block is found. This is especially the case in the
    ! beginning of the algorithm, as no matrix block has been found yet and
    ! Itr1 is completely 0. But also later on if the matrix contains
    ! independent blocks, this will work correctly in contrast to the
    ! old FEAT code.

    ! The renumbering strategy of Cuthill-McKee orders the nodes for
    ! increasing degree in the graph of the matrix. For the new numbering
    ! we use the index variable ineqIdx, which is increased for every node
    ! and specifies always the new node number for the next found node.
    ! We start with node number ineqIdx=0 and increase ineqIdx before
    ! the node number of the next found node is set.
    ineqIdx = 0

    ! Now loop over all nodes (=lines in the matrix). The following loop
    ! processes every node in a kind of quere exactly one time, so it
    ! has neq passes:
    DO icount=1, neq

      ! Look at line Itr1(icount) (=0, if a new matrix block starts, like
      ! in the beginning). We look at every element in that line and save
      ! the numbers of the following nodes in Itr1. These are processed later.

      ! If Itr1(icount)=0, a new matrix block starts. This is the case if
      ! a new independent block arises in the matrix without connection
      ! to previous elements. This is normally not the case with finite
      ! elements, but may arise in special cases (e.g. compressed matrices,
      ! where in a line there is only a diagonal entry = 1x1-block).
      ! In this case we ignore the number of the first line of that block
      ! (i.e. we give the next free number directly), and continue with
      ! processing the neighbours of the corresponding node in the graph
      ! of the matrix.
      IF (Itr1(icount) .EQ. 0) THEN

        ! New block. all previous blocks are contained in line 1..icount-1.
        ! The next line number to fix is therefore icount.          
        Itr1(icount) = icount

        ! Now the block structure of the matrix (representing the connected 
        ! components in the graph of the matrix) and the structure of the
        ! loop implies:
        !    Itr1(icount) = icount = ineqIdx + 1   !!!
        ! because: If we processed icount-1 elements and start a new block,
        ! the nodes contained in the previous blocks don't contain neighbours
        ! in the current block or in later blocks. Therefore all indices
        ! are <= icount-1, and at the same time ineqIdx=icount-1 holds
        ! (is not set but implicitly the case).

        ! The inverse of the entry is calculated in the DO-loop below,
        ! which also increments ineqIdx appropriately.

      END IF

      ! Now loop about the elements in the line and collect the neighbours
      ! of the corresponding node in the graph:        
      DO ildIdx=Ild(Itr1(icount)), Ild(Itr1(icount)+1)-1

        ! Collect the column numbers in the current line in the order of
        ! increasing degree of the nodes. icolIdx will receive the column
        ! numbers subsequently.
        icolIdx=Icon(ildIdx)

        ! Test if we already calculated the permutation for the current node.
        ! For that purpose analyze if the Itr2-entry of the currently
        ! processed neighbour contains a value:
        IF (Itr2(icolIdx) .EQ. 0) THEN

          ! This is a new node, which follows our current node Itr1(icount).
          ! Give it the next free number.
	  
          ! Obviously for the very first element there is:
          !   Itr1(1)=1, icolIdx=Icon(1)=1, Itr2(1)=0, ineqIdx=0.
          ! That way this routine correctly sets Itr2(1)=1 and correctly
          ! calculates the inverse of the permutation.
          ineqIdx       = ineqIdx+1
          Itr2(icolIdx) = ineqIdx

          ! Now insert the new node in Itr1. That way Itr1 also serves as a kind
          ! of queue of the nodes that have not yet been processed - subsequently
          ! Itr1(1), Itr1(2), Itr1(3),... is build up. Because in all cases
          ! icount <= ineqIdx holds, not calculated elements in Itr1 are never
          ! accessed.
          Itr1(ineqIdx) = icolIdx

        END IF  

      END DO
          
    END DO

  END SUBROUTINE cmsort_calcPermutation	  
        

END MODULE

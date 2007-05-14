!#########################################################################
!# ***********************************************************************
!# <name> sortstrategy </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains different routines to calculate resorting
!# strategies for scalar vectors.
!# A sorting strategy is simply a permutation of the numbers 1..NEQ
!# how to permute a scalar vector and its inverse permutation. 
!# Both permutations are usually assigned as one large array
!# to a (scalar) vector to indicate how it's resorted.
!#
!# Which information are necessary for the calculation of a permutation
!# is completely algorithm dependent. Therefore, this module has full
!# access to the whole (spatial) discretisation including triangulation
!# the domain, matrix data and more. It depends on the algorithm which
!# information is actually used.
!#
!# The following routines can be found in this module:
!#
!# 1.) sstrat_calcCuthillMcKee
!#     -> Calculate column numbering using the Cuthill McKee algorithm
!#
!# 2.) sstrat_calcXYZsorting
!#     -> Calculates rowwise renumbering based on Cartesian coordinates
!#
!# 2.) sstrat_calcFEASTsorting
!#     -> Calculates FEAST rowwise renumbering based cell adjacencies
!#
!# Auxiliary routines:
!#
!# 1.) sstrat_calcColNumberingCM7
!#     -> Calculate Cuthill-Mc-Kee resorting strategy for a structure-7
!#        matrix
!#
!# 2.) sstrat_calcColNumberingCM9
!#     -> Calculate Cuthill-Mc-Kee resorting strategy for a structure-9
!#        matrix
!#
!# 3.) sstrat_calcPermutationCM
!#     -> Auxiliary routine to calculate a permutation from the output
!#        of sstrat_calcColNumberingCM7 or sstrat_calcColNumberingCM9.
!#
!# 4.) sstrat_calcInversePermutation
!#     -> Calculates an inverse permutation
!#
!# </purpose>
!#########################################################################

MODULE sortstrategy

  USE fsystem
  USE storage
  USE genoutput
  USE linearalgebra
  USE spatialdiscretisation
  USE element
  USE triangulation
  USE linearsystemscalar

  IMPLICIT NONE
  
!<constants>

!<constantblock description="Sort strategy identifiers">

  ! No sort strategy; this must be =0!
  INTEGER, PARAMETER :: SSTRAT_UNSORTED     = 0

  ! Cuthill-McKee sort strategy
  INTEGER, PARAMETER :: SSTRAT_CM           = 1

  ! Reverse Cuthill-McKee sort strategy
  INTEGER, PARAMETER :: SSTRAT_RCM          = 2
  
  ! Row-wise sorting for point coordinate. 
  ! (As calculated by sstrat_calcXYZsorting with idirection=0.)
  ! Only for special type of discretisations ($Q_1$, $\tilde Q_1$), where 
  ! the DOF's can be identified with X/Y/Z coordinates.
  ! Coincides with the sorting strategy of FEAST for simple-type domains
  ! like the unit square.
  INTEGER, PARAMETER :: SSTRAT_XYZCOORD     = 3
  
  ! Column-wise sorting for point coordinate. 
  ! (As calculated by sstrat_calcXYZsorting with idirection=1.)
  ! Only for special type of discretisations ($Q_1$, $\tilde Q_1$), where 
  ! the DOF's can be identified with X/Y/Z coordinates.
  INTEGER, PARAMETER :: SSTRAT_ZYXCOORD     = 4
  
  ! General FEAST renumbering.
  ! The DOF's are numbered rowwise, independent of the geometrical
  ! structure of the domain.
  ! Only for special type of discretisations ($Q_1$) and tensor product meshes.
  INTEGER, PARAMETER :: SSTRAT_FEAST        = 5
  
!</constantblock>

!</constants>

CONTAINS

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE sstrat_calcCuthillMcKee (rmatrix,Ipermutation)
  
  !<description>
    ! Computes a column renumbering strategy using the algorithm
    ! of Cuthill-McKee. The algorithm acceps a scalar matrix rmatrix and
    ! uses its structure to calculate the renumbering. The result
    ! Ipermutation then receives the permutation and its inverse.
  !</description>
    
  !<input>
    ! Matrix which should be used to calculate the renumbering strategy
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
  !</input>
    
  !<output>
    ! The permutation vector for sorting and its inverse.
    ! With NEQ=NEQ(matrix):
    !   Ipermutation(1:NEQ)       = permutation,
    !   Ipermutation(NEQ+1:2*NEQ) = inverse permutation.
    INTEGER(PREC_VECIDX), DIMENSION(2*rmatrix%neq), INTENT(OUT) :: Ipermutation
  !</output>    

  !</subroutine>
  
  ! local variables
  INTEGER :: h_Ideg,h_IcolTmp
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Ideg
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_IcolTmp, p_Kcol,p_Kdiag
  INTEGER(PREC_VECIDX) :: NEQ
  
  NEQ = rmatrix%NEQ
  
  ! Currently, only matrix structure 7 and 9 are supported:
  SELECT CASE (rmatrix%cmatrixFormat)
  CASE (LSYSSC_MATRIX9)
    ! At first, duplicate KCOL and also get a temporary Ideg array
    h_IcolTmp = ST_NOHANDLE
    CALL storage_copy(rmatrix%h_Kcol,h_IcolTmp)
    CALL storage_new('sstrat_calcCuthillMcKee', 'KDEG', rmatrix%NEQ, &
                     ST_INT, h_Ideg, ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int(h_IcolTmp, p_IcolTmp)
    CALL lsyssc_getbase_Kcol(rmatrix, p_Kcol)
    CALL storage_getbase_int(h_Ideg, p_Ideg)
    CALL lsyssc_getbase_Kld(rmatrix, p_Kld)
    CALL lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiag)
    
    ! Calculate the strategy, calculate p_IcolTmp
    ! BETA STATUS, ROUTINE NOT TESTED!!!
    CALL sstrat_calcColNumberingCM9 (p_Kld, p_Kcol, p_Kdiag,&
                                     p_IcolTmp, p_Ideg, NEQ, NEQ)
    
    ! Use p_IcolTmp to calculate the actual resorting permutation
    ! and its inverse. Clear the target vector before calling the
    ! calculation routine, as we are creating a new permutation!
    CALL lalg_clearVectorInt(Ipermutation(1:NEQ*2))
    CALL sstrat_calcPermutationCM (p_Kld, p_IcolTmp, NEQ, &
                                   Ipermutation(1:NEQ), Ipermutation(NEQ+1:NEQ*2))

    ! Release temp data.
    CALL storage_free (h_Ideg)
    CALL storage_free (h_IcolTmp)    

  CASE (LSYSSC_MATRIX7)
    ! At first, duplicate KCOL and also get a temporary Ideg array
    h_IcolTmp = ST_NOHANDLE
    CALL storage_copy(rmatrix%h_Kcol,h_IcolTmp)
    CALL storage_new('sstrat_calcCuthillMcKee', 'KDEG', rmatrix%NEQ, &
                     ST_INT, h_Ideg, ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_int(h_IcolTmp, p_IcolTmp)
    CALL lsyssc_getbase_Kcol(rmatrix, p_Kcol)
    CALL storage_getbase_int(h_Ideg, p_Ideg)
    CALL lsyssc_getbase_Kld(rmatrix, p_Kld)
    
    ! Calculate the strategy, calculate p_IcolTmp
    CALL sstrat_calcColNumberingCM7 (p_Kld, p_Kcol, p_IcolTmp, p_Ideg, NEQ, NEQ)
    
    ! Use p_IcolTmp to calculate the actual resorting permutation
    ! and its inverse. Clear the target vector before calling the
    ! calculation routine, as we are creating a new permutation!
    CALL lalg_clearVectorInt(Ipermutation(1:NEQ*2))
    CALL sstrat_calcPermutationCM (p_Kld, p_IcolTmp, NEQ, &
                                   Ipermutation(1:NEQ), Ipermutation(NEQ+1:NEQ*2))

    ! Release temp data.
    CALL storage_free (h_Ideg)
    CALL storage_free (h_IcolTmp)    
  
  CASE DEFAULT
    PRINT *,'sstrat_calcCuthillMcKee: Unsupported matrix format'
    STOP
  END SELECT
    

  END SUBROUTINE     

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE sstrat_calcColNumberingCM7 (Ild, Icol, Icon, Ideg, neq, ndeg)
  
    !<description>
    ! Purpose: Cuthill McKee matrix renumbering
    !          Calculate column numbering
  
    ! The algorithm of Cuthill McKee interprets the system matrix as
    ! adjacense matrix. Every Row denotes a node in the corresponing graph.
    ! In the first step this function sorts the columns in every row
    ! for increasing degree of the nodes in the matrix.
    ! The matrix must have symmetric structure!
    ! (For FE matrices this is always the case...)
    !
    ! Matrix storage technique 7 version.
    !</description>
  
  !<input>
  
    ! Number of equations
    INTEGER(I32),INTENT(IN)                    :: neq

    ! Maximum number of entries != 0 in every row of the matrix
    INTEGER(PREC_VECIDX), INTENT(IN)                   :: ndeg    
   
    ! Row description of matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(IN) :: Ild
  
    ! Column description of matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)     :: Icol
    
  !</input>

  !<inputoutput>
    ! Auxiliary vector; must be at least as long as the
    ! maximum number of entries != 0 in every row of the matrix
    INTEGER(PREC_VECIDX), DIMENSION(ndeg), INTENT(INOUT) :: Ideg

    ! Auxiliary vector; the column numbers of KCOL are assigned to this in
    ! the order of increasing degree. When calling the routine the user
    ! must copy the content of KCOL to this! These values are then
    ! resorted.
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT)  :: Icon
  !</inputoutput>

!</subroutine>

    ! local variables
    INTEGER(I32) :: ieq, idegIdx, ildIdx, idegIdx1, idegIdx2
    INTEGER(I32) :: idegMin, iidxMin
  
    ! Clear auxiliary vector
    Ideg(1:ndeg) = 0
    !DO idegIdx=1, ndeg
    !  Ideg(idegIdx) = 0
    !END DO
    
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

        ! The variable idegMin always contains the degree of the node, which
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
  
  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE sstrat_calcColNumberingCM9 (Ild, Icol, Idiag, Icon, Ideg, &
                                         neq, ndeg)
  
    !<description>
    ! Purpose: Cuthill McKee matrix renumbering
    !          Calculate column numbering
  
    ! The algorithm of Cuthill McKee interprets the system matrix as
    ! adjacense matrix. Every Row denotes a node in the corresponing graph.
    ! In the first step this function sorts the columns in every row
    ! for increasing degree of the nodes in the matrix.
    ! The matrix must have symmetric structure!
    ! (For FE matrices this is always the case...)
    !
    ! Matrix storage technique 9 version.
    !</description>
  
  !<input>
  
    ! Number of equations
    INTEGER(I32),INTENT(IN)                    :: neq

    ! Maximum number of entries != 0 in every row of the matrix
    INTEGER(PREC_VECIDX), INTENT(IN)                   :: ndeg    
   
    ! Row description of matrix
    INTEGER(PREC_VECIDX), DIMENSION(neq+1), INTENT(IN) :: Ild
  
    ! Column description of matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)     :: Icol
    
    ! Incides of diagonal elements in structure 9 matrix
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)     :: Idiag
    
  !</input>

  !<inputoutput>
    ! Auxiliary vector; must be at least as long as the
    ! maximum number of entries != 0 in every row of the matrix
    INTEGER(PREC_VECIDX), DIMENSION(ndeg), INTENT(INOUT) :: Ideg

    ! Auxiliary vector; the column numbers of KCOL are assigned to this in
    ! the order of increasing degree. When calling the routine the user
    ! must copy the content of KCOL to this! These values are then
    ! resorted.
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT)  :: Icon
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
      Icon(Idiag(ieq)) = Icol(Idiag(ieq))
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
      DO ildIdx=Ild(ieq), Ild(ieq+1)-1
        Ideg(ildIdx-Ild(ieq)+1) = Icol(ildIdx)
      END DO

      ! Set Ideg of the diagonal entry to 0 to prevent it from being 
      ! processed later.
      Ideg(Idiag(ieq)-Ild(ieq)+1) = 0
      
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
      DO idegIdx1=1, Ild(ieq+1)-Ild(ieq)
      
        ! It's only necessary to search in the non-diagonal entries,
        ! because we want to sort the adjacent nodes of the current one
        ! for increasing degree.
        IF (Icol(Ild(ieq)+idegIdx1-1) .NE. ieq) THEN
        
          ! iidxMin will receive the index in the Ideg-array of the node with
          ! the smallest degree. We start with node 1 in the current column:
          iidxMin=1

          ! The variable idegMin always contains the degree of the node, which
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
          DO idegIdx2=1, Ild(ieq+1)-Ild(ieq)

            ! If Ideg(idegIdx2)=0, idegIdx2 has already been processed (or it's
            ! the diagonal element which does not have to be processed). Set
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
    
        END IF
      
      END DO
        
      ! Clear auxiliary vector; only some entries were used. This is only for
      ! reasons of safetyness, as if the upper loops are processed correctly,
      ! (no nodes were forgotten), all Ideg-arrays should already be 0.
      DO idegIdx=1, Ild(ieq+1)-Ild(ieq)
        Ideg(idegIdx) = 0
      END DO
      
    ! Process next line:  
    END DO 
  
  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE sstrat_calcPermutationCM (Ild, Icon, neq, Itr1, Itr2)
  
    !<description>
    ! Purpose: Cuthill McKee matrix renumbering
    !          Calculate permutation

    ! Uses the Ild/Icon vectors of sstrat_calcColumnNumbering to calc-
    ! ulate permutation vectors Itr1 and Itr2; these describe how a
    ! vector must be sorted/sorted back.

    ! Itr1 must be initialized with 0 on call of this subroutine!
    ! Itr2 should be initialized with 0 on call of this subroutine,
    ! as only the entries != 0 are calculated!

    ! Matrix storage technique 7 version!
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

  END SUBROUTINE     

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE sstrat_calcXYZsorting (rdiscretisation,Ipermutation,idirection)
  
  !<description>
    ! Computes a column renumbering strategy based on the coordinates of
    ! DOF's. The routine supports 2D and 3D triangulations; in a 2D
    ! triangulation, the Z-coordinate is ignored.
    !
    ! idirection specifies the sorting direction.
    ! The standard sorting direction is idirection=0. In this case, 
    ! the DOF's are sorted rowwise, i.e. first for the X-coordinate,
    ! then for the Y-coordinate. This sorting strategy coincides with the 
    ! FEAST sorting strategy in simple situations like a unit square.
    !
    ! This sorting strategy can only applied for special type discretisations
    ! where the DOF's of the finite elements coincides with some sort of
    ! point coordinates in the domain (like $Q_1$, edge midpoint based
    ! $\tilde Q_1$). If this is not the case, the routine will stop the 
    ! program.
    !
    ! The algorithm acceps a scalar discretisation structure rdiscretisation and
    ! uses its structure to calculate the renumbering. The result
    ! Ipermutation then receives the permutation and its inverse.
  !</description>
    
  !<input>
    ! Spatial discretisation structure that specifies the DOF's and the
    ! triangulation
    TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscretisation
    
    ! OPTIONAL: Specifies the sorting direction.
    ! =0: Sort first for X-, then for Y-, then for Z-coordinate.
    !     In 2D this is rowwise sorting.
    ! =1: Sort first for Z-, then for Y-, then for X-coordinate.
    !     In 2D this is columnwise sorting.
    ! If not specified, idirection=0 is assumed.
    INTEGER, INTENT(IN),OPTIONAL :: idirection
  !</input>
    
  !<output>
    ! The permutation vector for sorting and its inverse.
    ! With NEQ=NEQ(matrix):
    !   Ipermutation(1:NEQ)       = permutation,
    !   Ipermutation(NEQ+1:2*NEQ) = inverse permutation.
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Ipermutation
  !</output>    

  !</subroutine>
  
    ! local variables
    REAL(DP), DIMENSION(:,:), POINTER :: p_Dcoords
    INTEGER(PREC_EDGEIDX), DIMENSION(2) :: Isize
    INTEGER :: hhandle
    INTEGER :: idir
    
    IF (rdiscretisation%ndimension .EQ. 0) THEN
      CALL output_line ('Discretisation not initialised.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcRowwise')
      STOP
    END IF
    
    idir = 0
    IF (PRESENT(idirection)) idir=idirection

    ! Depending on the discrisation, choose the correct implementation.
    SELECT CASE (rdiscretisation%ccomplexity)
    CASE (SPDISC_UNIFORM)
    
      ! FE-space?
      SELECT CASE (elem_getPrimaryElement(&
                       rdiscretisation%RelementDistribution(1)%itrialElement))
      CASE (EL_Q1,EL_P1)
      
        ! $Q_1$-element. Take the vertex coordinates as DOF's and sort for that.
        CALL storage_getbase_double2d (&
            rdiscretisation%p_rtriangulation%h_DvertexCoords,p_Dcoords)
            
        CALL sortCoords (p_Dcoords, Ipermutation(1:UBOUND(p_Dcoords,2)), idir)
      
      CASE (EL_Q1T)
      
        ! $\tilde Q_1$-element. Take the edge midpoint coordinates as DOF's
        ! and sort for that. We have to calculate the midpoints for that...
        Isize(1) = rdiscretisation%p_rtriangulation%ndim
        Isize(2) = rdiscretisation%p_rtriangulation%NMT
        CALL storage_new2D ('rowwiseSorting', 'Dmidpoints', Isize, ST_DOUBLE, &
                            hhandle, ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_double2d (hhandle,p_Dcoords)
        
        ! Call tria_getPointsOnEdge with npointsPerEdge=1; this calculates
        ! the midpoint coordinates.
        CALL tria_getPointsOnEdge (rdiscretisation%p_rtriangulation,p_Dcoords,1)
        
        ! Sort for the midpoint coordinates
        CALL sortCoords (p_Dcoords, Ipermutation(1:UBOUND(p_Dcoords,2)), idir)
        
        ! Release temp array, finish
        CALL storage_free (hhandle)
      
      CASE DEFAULT
        CALL output_line ('Element type not supported.', &
                          OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcRowwise')
        STOP
      END SELECT
    
    CASE DEFAULT
      CALL output_line ('Discretisation too complex.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcRowwise')
      STOP
    END SELECT
  
    ! Calculate the inverse permutation, that's it.
    CALL sstrat_calcInversePermutation (Ipermutation(1:UBOUND(p_Dcoords,2)), &
                                        Ipermutation(UBOUND(p_Dcoords,2)+1:) )
    
  CONTAINS
  
    SUBROUTINE sortCoords (Dcoords, Ipermutation,idirection)

    ! Calculates the rowwise sorting. Dcoords must contain a 2D or 3D
    ! array of point coordinates. In Ipermutation, the routine returns
    ! the permutation of these points.
    ! idirection specifies the sorting direction.
    ! =0: points are first ordered for X-, then for Y- and at 
    !     the end for the Z-coordinate (in 3D).
    ! =1: points are first ordered for Z- (in 3D), then for Y- and at 
    !     the end for the X-coordinate.
    
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dcoords
    INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Ipermutation
    INTEGER, INTENT(IN) :: idirection
    
      ! local variables
      INTEGER(I32), DIMENSION(2) :: Isize
      INTEGER :: h_Dsort
      INTEGER :: i
      REAL(DP), DIMENSION(:,:), POINTER :: p_Dsort
      
      ! Allocate a 2D array with (dim(Dcoords)+1,#coords) elements.
      Isize(1) = UBOUND(Dcoords,1)+1
      Isize(2) = UBOUND(Dcoords,2)
      CALL storage_new2D ('rowwiseSorting', 'Dsort', Isize, ST_DOUBLE, &
                          h_Dsort, ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double2d (h_Dsort,p_Dsort)
      
      ! In the first element of each ndim+1-tupel, store the number of
      ! the point. In the 2nd/3rd,... element, store the coordinate.
      DO i=1,Isize(2)
        p_Dsort(1,i) = REAL(i,DP)
        p_Dsort(2:,i) = Dcoords(:,i)
      END DO
      
      ! Sort the array. First for the last coordinate, then for the
      ! last but one, etc.
      ! Use a stable sorting algorithm to prevent the previous sorting
      ! from getting destroyed.
      SELECT CASE (idirection)
      CASE (0)
        DO i=1,UBOUND(Dcoords,1)
          CALL arraySort_sortByIndex_dp(p_Dsort,1+i,SORT_STABLE)
        END DO
        
      CASE (1)
        DO i=UBOUND(Dcoords,1),1,-1
          CALL arraySort_sortByIndex_dp(p_Dsort,1+i,SORT_STABLE)
        END DO
      
      END SELECT
      
      ! The first element in each ndim+2-tupel is now the permutation.
      ! Do a type conversion to int to get it.
      DO i=1,Isize(2)
        Ipermutation(i) = INT(p_Dsort(1,i),I32)
      END DO
      
      ! Release the temp array, that's it.
      CALL storage_free (h_Dsort)
      
    END SUBROUTINE

  END SUBROUTINE     

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE sstrat_calcFEASTsorting (rdiscretisation,Ipermutation,ifirstVertex)
  
  !<description>
    ! Computes a column renumbering strategy based on the tensor product
    ! structure of the domain. The DOF's are numbered rowwise.
    ! ifirstVertex must specify the first vertex (in the lower left corner)
    ! of the domain. From here, a geometrical search is started to find
    ! all DOF's.
    !
    ! The renumbering strategy is only applicable to special-type 
    ! discretisations (like $Q_1$) and tensor product meshes (FEAST macros).
    ! If this is not the case, the routine will stop the program.
    !
    ! The algorithm acceps a scalar discretisation structure rdiscretisation and
    ! uses its structure to calculate the renumbering. The result
    ! Ipermutation then receives the permutation and its inverse.
  !</description>
    
  !<input>
    ! Spatial discretisation structure that specifies the DOF's and the
    ! triangulation
    TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscretisation
    
    ! OPTIONAL: First vertex (lower left corner) of the domain.
    ! If not specified, ifirstVertex=1 is assumed.
    INTEGER, INTENT(IN), OPTIONAL :: ifirstVertex
  !</input>
    
  !<output>
    ! The permutation vector for sorting and its inverse.
    ! With NEQ=NEQ(matrix):
    !   Ipermutation(1:NEQ)       = permutation,
    !   Ipermutation(NEQ+1:2*NEQ) = inverse permutation.
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Ipermutation
  !</output>    

  !</subroutine>
  
    ! local variables
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IelementsAtEdge
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
    INTEGER(I32), DIMENSION(:), POINTER :: p_IelementsAtVertexIdx
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertex
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_VERTEXIDX) :: ivt
    
    IF (rdiscretisation%ndimension .EQ. 0) THEN
      CALL output_line ('Discretisation not initialised.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcFEASTsorting')
      STOP
    END IF
    
    ! Depending on the discrisation, choose the correct implementation.
    SELECT CASE (rdiscretisation%ccomplexity)
    CASE (SPDISC_UNIFORM)
    
      ! FE-space?
      SELECT CASE (elem_getPrimaryElement(&
                       rdiscretisation%RelementDistribution(1)%itrialElement))
      CASE (EL_Q1)
      
        ! Get geometrical data
        CALL storage_getbase_int2d (rdiscretisation%p_rtriangulation%h_IelementsAtEdge,&
            p_IelementsAtEdge)
        CALL storage_getbase_int2d (&
            rdiscretisation%p_rtriangulation%h_IverticesAtElement,&
            p_IverticesAtElement)
        CALL storage_getbase_int2d (&
            rdiscretisation%p_rtriangulation%h_IedgesAtElement,&
            p_IedgesAtElement)
            
        ! Get the first element on vertex ifirstVertex
        ivt = 1
        IF (PRESENT(ifirstVertex)) ivt = ifirstVertex
        CALL storage_getbase_int (rdiscretisation%p_rtriangulation%h_IelementsAtVertex,&
            p_IelementsAtVertex)
        CALL storage_getbase_int (&
            rdiscretisation%p_rtriangulation%h_IelementsAtVertexIdx,&
            p_IelementsAtVertexIdx)
        iel = p_IelementsAtVertex(p_IelementsAtVertexIdx(ivt))
      
        CALL sortForFeastQ1 (p_IelementsAtEdge, p_IverticesAtElement, p_IedgesAtElement,&
            Ipermutation(1: rdiscretisation%p_rtriangulation%NVT), ivt, iel, &
            rdiscretisation%p_rtriangulation%NVT)
      
      CASE DEFAULT
        CALL output_line ('Element type not supported.', &
                          OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcFEASTsorting')
        STOP
      END SELECT
    
    CASE DEFAULT
      CALL output_line ('Discretisation too complex.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcFEASTsorting')
      STOP
    END SELECT
  
    ! Calculate the inverse permutation, that's it.
    CALL sstrat_calcInversePermutation (Ipermutation(1:rdiscretisation%p_rtriangulation%NVT), &
                                        Ipermutation(rdiscretisation%p_rtriangulation%NVT+1:) )
    
  CONTAINS
  
    SUBROUTINE sortForFeastQ1 (IelementsAtEdge,IverticesAtElement,IedgesAtElement,&
                               Ipermutation,ivt,iel,NVT)

    ! Calculates the rowwise sorting in a FEAST like style. 
    ! IelementsAtEdge, IverticesAtElement and IedgesAtElement specify the 
    ! adjacencies in the macro.
    ! In Ipermutation, the routine returns the permutation of these points.
    ! ivt specifies the lower left corner of the macro. iel specifies
    ! the element in the lower left corner that contains ivt.
    
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), INTENT(IN) :: IelementsAtEdge
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElement
    INTEGER(I32), DIMENSION(:), INTENT(OUT) :: Ipermutation
    INTEGER(PREC_VERTEXIDX), INTENT(IN) :: ivt
    INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: iel
    INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVT

      ! local variables
      INTEGER(PREC_VERTEXIDX) :: icornervertex,icornerelement,ipermidx
      INTEGER(PREC_VERTEXIDX) :: icurrentvertex,icurrentelement
      INTEGER(PREC_EDGEIDX) :: iedgeright,iedgetop
      INTEGER :: ilocalvertex
      INTEGER, PARAMETER :: NVE = 4
      
      ! Current position in the mesh
      icurrentvertex = ivt
      icurrentelement = iel
      
      ! ipermidx counts how many vertices we found.
      ipermidx = 0
      
      ! Loop through the columns until we find the macro border at the top
      DO
        
        ! icornervertex remembers the current lower left corner. icornerelement
        ! the corresponding element.
        icornervertex = icurrentvertex
        icornerelement = icurrentelement
        
        ! We are here:
        !
        !    |   |
        !    +---+---
        !    |IEL|
        !  IVT---+---
        !
        ! Add the first vertex to the permutation
        ipermidx = ipermidx+1
        Ipermutation(ipermidx) = icurrentvertex
        
        ! Get the local number of the vertex in the element
        DO ilocalvertex = 1,NVE
          IF (IverticesAtElement(ilocalvertex,icurrentelement) .EQ. icurrentvertex) EXIT
        END DO
        
        ! Get the edges to the neighbour elements
        !
        !    |   |     
        !    +---+---
        !    |   X iedgeright
        !  IVT---+---
        
        iedgeright = IedgesAtElement(MOD(ilocalvertex,NVE)+1,icurrentElement)-NVT
      
        ! Loop through the macro row until we find the right border of the macro
        DO WHILE (IelementsAtEdge(2,iedgeright) .NE. 0) 
        
          ! Step right to the next element 
          !
          !    |   |     
          !    +---+---+
          !    |   |   |           
          !    +--IVT--+
          
          icurrentvertex = IverticesAtElement(MOD(ilocalvertex,NVE)+1,icurrentElement)
          
          IF (IelementsAtEdge(2,iedgeright) .NE. icurrentElement) THEN
            icurrentElement = IelementsAtEdge(2,iedgeright)
          ELSE
            icurrentElement = IelementsAtEdge(1,iedgeright)
          END IF
          
          ! Add the vertex to the permutation
          ipermidx = ipermidx+1
          Ipermutation(ipermidx) = icurrentvertex
          
          ! Get the local number of the vertex in the element
          DO ilocalvertex = 1,NVE
            IF (IverticesAtElement(ilocalvertex,icurrentelement) &
                .EQ. icurrentvertex) EXIT
          END DO
        
          ! Get the edges to the neighbour elements
          !
          !    |   |   |  
          !    +---+---+--
          !    |   |   X iedgeright
          !    +--IVT--+--
          
          iedgeright = IedgesAtElement(MOD(ilocalvertex,NVE)+1,icurrentElement)-NVT
          
        END DO
        
        ! We have reached the end of the row
        !
        !        |   |   |
        !     ---+---+---+
        !        |   |   X iedgeright
        !     ---+--IVT--IVT2
        !
        ! Remember the last vertex IVT2 in the row
        
        icurrentvertex = IverticesAtElement(MOD(ilocalvertex,NVE)+1,icurrentelement)
        ipermidx = ipermidx+1
        Ipermutation(ipermidx) = icurrentvertex
        
        ! Hop back to the first element and find the local number of the first
        ! vertex in the row again.
        !
        !    |   |
        !    +---+---
        !    |IEL|
        !  IVT---+---

        icurrentvertex = icornervertex
        icurrentelement = icornerelement
        
        DO ilocalvertex = 1,NVE
          IF (IverticesAtElement(ilocalvertex,icurrentelement) .EQ. icurrentvertex) EXIT
        END DO
        
        ! Get the edge that leads to the element row above us.
        !
        !    | iedgetop
        !    +-X-+---
        !    |IEL|
        !  IVT---+---
        
        iedgetop = IedgesAtElement(MOD(ilocalvertex+1,NVE)+1,icurrentElement)-NVT
        
        ! Get the vertex and the element there. Note: If there is no neighbour
        ! element, the current element number gets =0 which is the terminal criterion.
        !
        !    |IEL|     
        !  IVT---+---
        !    |   |
        !    +---+---
      
        icurrentvertex = IverticesAtElement(MOD(ilocalvertex+2,NVE)+1,icurrentElement)
      
        IF (IelementsAtEdge(2,iedgetop) .NE. icurrentElement) THEN
          icurrentElement = IelementsAtEdge(2,iedgetop)
        ELSE
          icurrentElement = IelementsAtEdge(1,iedgetop)
        END IF
      
        IF (icurrentelement .EQ. 0) THEN
          EXIT
        ELSE
          ! There is a neighbour element on top.
          ! Get the edge that leads to the right and continue in the new
          ! element row.
          DO ilocalvertex = 1,NVE
            IF (IverticesAtElement(ilocalvertex,icurrentelement) .EQ. icurrentvertex) EXIT
          END DO

          iedgeright = IedgesAtElement(MOD(ilocalvertex,NVE)+1,icurrentElement)-NVT
        END IF
      END DO
      
      ! Ok, we have done all element rows now.
      !
      !   |   |   |   |   |
      !   +---+---+---+---+
      !   |   |   |   |   |           
      !   +---+---+---+---+
      !
      ! What's still missing is the final edge-row on top! This is another loop.
      ! Switch the element number back to the last remembered one, then we have: 
      !
      !  IVT--+---+---+---+
      !   |IEL|   |   |   |           
      !   +---+---+---+---+
      !   |   |   |   |   |
      
      icurrentelement = icornerelement
      
      ! Add the vertex to the permutation
      ipermidx = ipermidx+1
      Ipermutation(ipermidx) = icurrentvertex
      
      ! Local number and edge to the right? 
      !
      !  IVT--+--           
      !   |   X iedgeright             
      !   +---+-            
      !   |   |            

      DO ilocalvertex = 1,NVE
        IF (IverticesAtElement(ilocalvertex,icurrentelement) .EQ. icurrentvertex) EXIT
      END DO
      
      iedgeright = IedgesAtElement(MOD(ilocalvertex+1,NVE)+1,icurrentElement)-NVT
    
      ! Loop through the macro row until we find the right border of the macro
      DO WHILE (IelementsAtEdge(2,iedgeright) .NE. 0) 
      
        ! Step right to the next element 
        !
        !    +--IVT--+
        !    |   |IEL|           
        !    +---+---+
        !    |   |     
        
        icurrentvertex = IverticesAtElement(MOD(ilocalvertex+2,NVE)+1,icurrentElement)
        
        IF (IelementsAtEdge(2,iedgeright) .NE. icurrentElement) THEN
          icurrentElement = IelementsAtEdge(2,iedgeright)
        ELSE
          icurrentElement = IelementsAtEdge(1,iedgeright)
        END IF
        
        ! Add the vertex to the permutation
        ipermidx = ipermidx+1
        Ipermutation(ipermidx) = icurrentvertex
        
        ! Get the local number of the vertex in the element
        DO ilocalvertex = 1,NVE
          IF (IverticesAtElement(ilocalvertex,icurrentelement) &
              .EQ. icurrentvertex) EXIT
        END DO
      
        ! Get the edges to the neighbour elements
        !
        !    +--IVT--+--
        !    |   |   X iedgeright
        !    +---+---+--
        !    |   |   |  
        
        iedgeright = IedgesAtElement(MOD(ilocalvertex+1,NVE)+1,icurrentElement)-NVT
        
      END DO    
      
      ! Remember the last vertex IVT2 in the row
      !
      !     --IVT--IVT2
      !        |   |           
      !     ---+---+  
      !        |   |  
      
      icurrentvertex = IverticesAtElement(MOD(ilocalvertex+2,NVE)+1,icurrentelement)
      ipermidx = ipermidx+1
      Ipermutation(ipermidx) = icurrentvertex

    END SUBROUTINE

  END SUBROUTINE     

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE sstrat_calcInversePermutation (IpermutationSource, IpermutationDest)
  
  !<description>
    ! Computes the inverse of a permutation. IpermutationSource is a given
    ! permutation. IpermutationDest will receive the inverse permutation.
  !</description>
    
  !<input>
    ! A permutation.
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: IpermutationSource
  !</input>
    
  !<output>
    ! An array of the same size as IpermutationSource. Receives the inverse
    ! permutation.
    INTEGER(I32), DIMENSION(:), INTENT(OUT) :: IpermutationDest
  !</output>    

  !</subroutine>
  
    INTEGER :: i
    DO i=1,SIZE(IpermutationSource)
      IpermutationDest(IpermutationSource(i)) = i
    END DO
  
  END SUBROUTINE

END MODULE

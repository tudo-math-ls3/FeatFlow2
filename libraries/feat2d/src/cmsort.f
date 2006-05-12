************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek, M.Koester          *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* CUTCE0                                                               *
*                                                                      *
* Purpose: Cuthill McKee matrix renumbering                            *
*          Calculate column numbering                                  *
*                                                                      *
* The algorithm of Cuthill McKee interprets the system matrix as       *
* adjacense matrix. Every Row denotes a note in the corresponing graph.*
* In the first step this function sorts the columns in every row       *
* for increasing degree of the nodes in the matrix.                    *
*                                                                      *
* The matrix must have symmetric structure!!!                          *
* (For FE matrices this is always the case...)                         *
*                                                                      *
* Matrix storage technique 7 only!                                     *
*                                                                      *
* Version from  10/15/04                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KCOL     I*4(NA)    Column description of matrix                     *
* KLD      I*4(NEQ+1) Row description of matrix                        *
* NEQ      I*4        Number of equations                              *
* KCON     I*4(NA)    Auxiliary vector; the column numbers of KCOL     *
*                     are assigned to this in the order of increasing  *
*                     degree. When calling the routine the user must   *
*                     copy the content of KCOL to this!!! These values *
*                     are then resorted.                               *
* KDEG     I*4(NDEG)  Auxil. vector; must be at least as long as the   *
*                     maximum number of entries != 0 in every row      *
*                     of the matrix                                    *
* NDEG     I*4        Maximum number of entries != 0 in every row      *
*                     of the matrix                                    *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KCON     I*4(NA)    Resorted KCOL array                              *
*                                                                      *
************************************************************************

      SUBROUTINE CUTCE0(KLD,KCOL,KCON,KDEG,NEQ,NDEG)
      
      IMPLICIT NONE

C parameters
      
      INTEGER NEQ,NDEG,KLD(*),KCOL(*),KCON(*),KDEG(NDEG)

C local variables:

      INTEGER IDEG, IEQ, ILD, INDMIN, MINDEG, IDEG1, IDEG2, INDDEG

C Clear aux. vector

        DO IDEG=1,NDEG
          KDEG(IDEG)=0
        END DO  
        
C Later in the algorithm the column numbers of the non-diagonal entries
C are changed. The numbers of the diagonal entries are not touched.
C Therefore we copy the column numbers of the diagonal entries
C directly.

        DO IEQ=1,NEQ
          KCON(KLD(IEQ)) = KCOL(KLD(IEQ))
        END DO

C Loop about all rows in the matrix:   IEQ = current row.        
C Every row corresponds to an entry in the solution vector = a node
C in the graph of the matrix.

        DO IEQ=1,NEQ

C Copy the column numbers of the non-diagonal entries in the current
C column to the auxiliary vector KDEG. KDEG contains always the
C following nodes of the current node IEQ, which are not processed yet.
C The entries of the vector are set to 0 one after the other in the
C order of the degree of the following nodes.

          DO ILD=KLD(IEQ)+1,KLD(IEQ+1)-1
            KDEG(ILD-KLD(IEQ))=KCOL(ILD)
          END DO  

C Loop about every column in the current row. The entries in the
C row (=column numbers) of the matrix represent the node numbers of
C those nodes in the graph of the matrix, which are connected to the
C current node IEQ.
C
C The algorithm now has to sort these nodes for increasing degree.
C For this purpose the node with the smallest degree is searched and
C deleted from KDEG, then the node with the 2nd smallest degree
C and so on.
C We don't have many nodes, so we use a simple O(n^2) algorithm here
C which simply loops over all nodes.
C It's only necessary to search in the non-diagonal entries,
C because we want to sort the adjacent nodes of the current one
C for increasing degree.

          DO IDEG1=1,KLD(IEQ+1)-KLD(IEQ)-1

C INDMIN will receive the index in the KDEG-array of the node with
C the smallest degree. We start with node 1 in the current column:

            INDMIN=1

C The variable MINDEG always contains the degree of the node, which
C is described by the index INDMIN. 
C
C If KDEG(INDMIN)=0, this node has already been processed. In this
C case MINDEF is set to NEQ - the largest upper bound which implies
C that INDMIN is later on surely replaced by a node with a smaller
C degree.
C
C Otherwise MINDEG receives the degree of the node described by INDMIN,
C which is calculated by the number of elements in KCOL, i.e. the
C difference of the indices in the KLD array.

            IF (KDEG(INDMIN).EQ.0) THEN
              MINDEG=NEQ
            ELSE
              MINDEG=KLD(KDEG(INDMIN)+1)-KLD(KDEG(INDMIN))-1
            END IF

C Compare INDMIN with every node in that line to find that with minimum
C degree:

            DO IDEG2=1,KLD(IEQ+1)-KLD(IEQ)-1

C If KDEG(IDEG2)=0, IDEG2 has already been processed. Set INDDEG=NEQ -
C here a lower bound to prevent the current node INDMIN from being
C replaced by IDEG2.
C
C Otherwise set INDDEG to the degree of the node with index IDEG2.

              IF (KDEG(IDEG2).EQ.0) THEN
                INDDEG=NEQ
              ELSE
                INDDEG=KLD(KDEG(IDEG2)+1)-KLD(KDEG(IDEG2))-1
              END IF

C If now INDDEG=grad(INDMIN) < grad(IDEG2)=MINDEG, set INDMIN to
C the new node with smaller degree.

              IF (INDDEG.LT.MINDEG) THEN
                INDMIN=IDEG2
                MINDEG=INDDEG
              END IF

            END DO

C At this point, INDMIN contains the index in KDEG of that node with
C minimum degree, which has not yet been processed.
C KDEG(INDMIN) contains the column number of that column of the matrix,
C which corresponds with the node with the next higher degree.
C
C This value is written to KCON. If this subroutine is finished,
C KCON therefore contains the column numbers of the KCOL-array -
C sorted for increasing degree of the nodes in the graph of the matrix.
     
            KCON(KLD(IEQ)+IDEG1) = KDEG(INDMIN)

C Set the KDEG-value of the node with minimum degree to 0 to indicate
C that is has been processed. This node is then ignored in the next
C searching processes for nodes with minimum degree.
            
            KDEG(INDMIN) = 0
     
          END DO

C Clear auxiliary vector; only some entries were used. This is only for
C reasons of safetyness, as if the upper loops are processed correctly,
C (no nodes were forgotten), all KDEG-arrays should already be 0.

          DO IDEG=1,KLD(IEQ)+1,KLD(IEQ+1)-1
            KDEG(IDEG)=0
          END DO  

C Process next line:

        END DO

99999 END      

************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek, M.Koester          *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* CUTCE1                                                               *
*                                                                      *
* Purpose: Cuthill McKee matrix renumbering                            *
*          Calculate permutation                                       *
*                                                                      *
* Uses the KLD/KCON vectors of CUTCE0 to calculate permutation vectors *
* KTR1 and KTR2; these describe how a vector must be sorted/sorted     *
* back.                                                                *
*                                                                      *
* KTR1 must be initialized with 0 on call of this subroutine!!!        *
* KTR2 should be initialized with 0 on call of this subroutine; ther-  *
* wise only the entries != 0 are calculated!!!                         *
*                                                                      *
* Matrix storage technique 7 only!                                     *
*                                                                      *
* Version from  10/15/04                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KLD      I*4(NEQ+1) Row description of matrix                        *
* KCON     I*4(NA)    Auxiliary vector, calculated with CUTCE0         *
* NEQ      I*4        Number of equations                              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KTR1     I*4(NEQ)   The permutation vector that describes how the    *
*                     solution vector has to be resorted.              *
* KTR2     I*4(NEQ)   Describes how a resorted vector can be sorted    *
*                     back.                                            *
*                                                                      *
************************************************************************

      SUBROUTINE CUTCE1(KLD,KCON,NEQ,KTR1,KTR2)
      
      IMPLICIT NONE
      
C parameters
      
      INTEGER NEQ,KLD(*),KCON(*),KTR1(NEQ),KTR2(NEQ)
      
C local variables

      INTEGER INDNEQ, ICOUNT, ILD, ICOL

C        DO ICOUNT=1,NEQ
C          IF (KTR1(ICOUNT).EQ.0) KTR1(ICOUNT)=ICOUNT
C        END DO
C Commented out, old FEAT-Code. The current code is using the fact that
C the KTR1-vector is filled with 0 on call to this routine.

C KTR1 should receive a permutation of the numbers 1..NEQ for
C resorting the entries, KTR2 should receive the inverse permutation.

C The Cuthill-McCee algorithm processes levelwise. One node is picked
C as the starting node and receives the first number. The nodes
C in the neighbourhood receive the next numbers in increasing order
C of the degrees of these nodes. Then the neighbourhood of these
C nodes is analyzed for increasing degree to give the next node
C numbers, and so on.
C
C For this purpose the permutation vectors KTR1 and KTR2 are build
C subsequently, depending on which nodes have already been processed.        
C
C We fix the first node (=first line of the matrix), i.e. we set
C KTR1(1)=1. This of course implies that the inverse permutation
C does not have to process node 1, i.e. KTR2(1)=1. But this will
C be computed by the algorithm later.
        
C       KTR1(1)=1
C       KTR2(1)=1
C Old FEAT-code, commented out.
C The current code uses the fact, that KTR1 is initialized with 0 on
C call to this routine. In in the processing a 0 is found in KTR1, a
C new matrix block is found. This is especially the case in the
C beginning of the algorithm, as no matrix block has been found yet and
C KTR1 is completely 0. But also later on if the matrix contains
C independent blocks, this will work correctly in contrast to the
C old FEAT code.
C
C The renumbering strategy of Cuthill-McKee orders the nodes for
C increasing degree in the graph of the matrix. For the new numbering
C we use the index variable INDNEQ, which is increased for every node
C and specifies always the new node number for the next found node.
C We start with node number INDNEQ=0 and increase INDNEQ before
C the node number of the next found node is set.

        INDNEQ=0

C Now loop over all nodes (=lines in the matrix). The following loop
C processes every node in a kind of quere exactly one time, so it
C has NEQ passes:

        DO ICOUNT=1,NEQ

C Look at line KTR1(ICOUNT) (=0, if a new matrix block starts, like
C in the beginning). We look at every element in that line and save
C the numbers of the following nodes in KTR1. These are processed later.
C
C If KTR1(ICOUNT)=0, a new matrix block starts. This is the case if
C a new independent block arises in the matrix without connection
C to previous elements. This is normally not the case with finite
C elements, but may arise in special cases (e.g. compressed matrices,
C where in a line there is only a diagonal entry = 1x1-block).
C In this case we ignore the number of the first line of that block
C (i.e. we give the next free number directly), and continue with
C processing the neighbours of the corresponding node in the graph
C of the matrix.

          IF (KTR1(ICOUNT).EQ.0) THEN

C New block. all previous blocks are contained in line 1..ICOUNT-1.
C The next line number to fix is therefore ICOUNT.          

            KTR1(ICOUNT) = ICOUNT

C Now the block structure of the matrix (representing the connected 
C components in the graph of the matrix) and the structure of the
C loop implies:
C          KTR1(ICOUNT) = ICOUNT = INDNEQ + 1   !!!
C because: If we processed ICOUNT-1 elements and start a new block,
C the nodes contained in the previous blocks don't contain neighbours
C in the current block or in later blocks. Therefore all indices
C are <= ICOUNT-1, and at the same time INDNEQ=ICOUNT-1 holds
C (is not set but implicitly the case).
C
C The inverse of the entry is calculated in the DO-loop below,
C which also increments INDNEQ appropriately.

          END IF

C Now loop about the elements in the line and collect the neighbours
C of the corresponding node in the graph:        
        
          DO ILD=KLD(KTR1(ICOUNT)),KLD(KTR1(ICOUNT)+1)-1

C Collect the column numbers in the current line in the order of
C increasing degree of the nodes. ICOL will receive the column
C numbers subsequently.
          
            ICOL=KCON(ILD)

C Test if we already calculated the permutation for the current node.
C For that purpose analyze if the KTR2-entry of the currently
C processed neighbour contains a value:
            
            IF (KTR2(ICOL).EQ.0) THEN

C This is a new node, which follows our current node KTR1(ICOUNT).
C Give it the next free number.
C
C Obviously for the very first element there is:
C   KTR1(1)=1, ICOL=KCON(1)=1, KTR2(1)=0, INDNEQ=0.
C That way this routine correctly sets KTR2(1)=1 and correctly
C calculates the inverse of the permutation.

              INDNEQ      =INDNEQ+1
              KTR2(ICOL  )=INDNEQ

C Now insert the new node in KTR1. That way KTR1 also serves as a kind
C of queue of the nodes that have not yet been processed - subsequently
C KTR1(1), KTR1(2), KTR1(3),... is build up. Because in all cases
C ICOUNT <= INDNEQ holds, not calculated elements in KTR1 are never
C accessed.
              
              KTR1(INDNEQ)=ICOL
            END IF  
          END DO
          
        END DO
        
99999 END      

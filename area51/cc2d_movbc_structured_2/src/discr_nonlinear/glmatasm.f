************************************************************************
* This file contains routines that allow to assemble a global matrix
* from smaller sub-matrices
************************************************************************

************************************************************************
* Assemble global matrix
*
* This routine accepts a set of structure 7 / structure 9 matrices
* and assembles a global strucure 9 matrix out of them. The resulting
* matrix can in general not be assumed to be symmetric.
*
* The caller has to specify the structure of the global matrix by
* passing handles to sub-matrix blocks. There are NXBLOC matrices
* in horizontal and NYBLOC matrices in vertical direction stacked
* to each other. Starting addresses IA=ICOL=ILD=0 indicate a NULL block
* in the global matrix.
* The global matrix has the structure
*
*     A(1,1)        A(2,1) ... A(NXBLOC,1)
*     A(1,2)        A(2,2) ... A(NXBLOC,2)
*     ...           ...        ...
*     A(NBLOCY,1)   ...        A(NXBLOC,NYBLOC)
*
* At least one matrix block per row must exist!
* All matrices in the same row must have the same number of equations!
* The columns in each matrix block are assumed to be stored in
* ascending order except for the diagonal element in matrix structure 7
* matrix blocks.
*
* In:
*   NXBLOC : Number of matrix blocks in horizontal direction
*   NYBLOC : Number of matrix blocks in vertical direction
*   DA     : array [1..*] of double
*            Array that contains all matrix entries
*   KCOL   : array [1..*] of integer
*            Array with all KCOL entries
*   KLD    : array [1..*] of integer
*            Array with all KLD entries
*   IA     : array [1..NYBLOC,1..NXBLOC] of integer
*            Starting address of the matrix entries of all matrices
*            stored in DA. IA(i,j) is the starting address of the
*            matrix block (i,j). Its matrix entries can be found
*            at DA(IA(i,j)).
*   ICOL   : array [1..NYBLOC,1..NXBLOC] of integer
*            Starting address of the KCOL row structures of all matrices
*            stored in KCOL. ICOL(i,j) is the starting address of the
*            column structure of matrix block (i,j). The column 
*            structure itself is to be found at KCOL(ICOL(i,j)).
*   ILD    : array [1..NYBLOC,1..NXBLOC] of integer
*            Starting address of the KLD row structures of all matrices
*            stored in KLD. ILD(i,j) is the starting address of the
*            row structure of matrix block (i,j). The row 
*            structure itself is to be found at KLD(ICOL(i,j)).
*   KNEQ   : array [1..NYBLOC,1..NXBLOC] of integer
*            Number of equations in each matrix. KNEQ(i,j) gives
*            the number of equations in matrix block (i,j).
*   BSORT  : Sort the entries in ascending order.
*            =true:  Any leading diagonal elements of matrix structure 7
*                    are put to the "right" position in the global
*                    matrix. Can be used to prepare the call to an
*                    UMFPACK solvers e.g..
*            =false: Don't sort anything, simply copy the matrices
*                    as they are.
*
* Out:
*   LA     : Handle to the entries of the global matrix 
*   LCOL   : Handle to the column structure of the global matrix
*   LLD    : Handle to the row structure of the global matrix
*   NEQ    : Number of rows in the global matrix.
*            =-1 indicates an error.
*            =0 indicates an empty global matrix.
*            NEQ=-1 or NEQ=0 implies LA=LCOL=LLD=0!
************************************************************************

      SUBROUTINE MTXA79 (NXBLOC,NYBLOC,DA,KCOL,KLD,IA,ICOL,ILD,KNEQ,
     *                   BSORT,LA,LCOL,LLD,NEQ)
 
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
C     parameters
      
      INTEGER NXBLOC,NYBLOC,KCOL(*),KLD(*)
      INTEGER IA(NYBLOC,NXBLOC),ICOL(NYBLOC,NXBLOC),ILD(NYBLOC,NXBLOC)
      INTEGER KNEQ(NYBLOC,NXBLOC)
      INTEGER LA,LCOL,LLD,NEQ,PLEN
      LOGICAL BSORT
      DOUBLE PRECISION DA(*)
      
C     local variables

      INTEGER I,J,KR,KC,KA,ROW,COL,NEQ1,NEQ2,NA,RLEN,TLEN,P,GROW,CCOL
      INTEGER MA,MCOL,KAUX
      DOUBLE PRECISION DAUX
      INTEGER NROW(NYBLOC),NCOL(NXBLOC)
      
C     Small check

      NEQ  = 0
      LA   = 0
      LCOL = 0
      LLD  = 0
      
      IF ((NXBLOC.EQ.0).OR.(NYBLOC.EQ.0)) RETURN

C     Clear the number of rows/columns in the matrix blocks

      CALL LCL3(NROW,NYBLOC)
      CALL LCL3(NCOL,NXBLOC)

C     First let's calculate the global number of rows in that matrix.
C     We loop through all rows and columns of the matrix and check
C     at the same time, that all matrices in a row are either NULL
C     blocks or have all the same NEQ!

      NEQ  = 0
      NA   = 0
      
      DO I=1,NYBLOC
        
        NEQ1 = 0
        NEQ2 = 0
        
        DO J=1,NXBLOC
          
          IF (KNEQ(I,J).NE.0) THEN

C           In NEQ1 calculate the minimum number of rows, in NEQ2
C           the maximum number of rows in the matrix row I:
            
            NEQ2 = MAX(NEQ2,KNEQ(I,J))
            IF ((NEQ1.EQ.0).AND.(NEQ2.NE.0)) THEN
              NEQ1 = NEQ2
            ELSE
              NEQ1 = MIN(NEQ1,KNEQ(I,J))
            END IF
            
C           Calculate the number of columns into NCOL.
C           Loop through the lines of the matrix block (i,j) and
C           calculate the maximum column number - this is always
C           given by the last column of each row of each matrix block:
            
            IF ((ICOL(I,J).NE.0).AND.(ILD(I,J).NE.0)) THEN

              KR  = ILD(I,J)
              KC  = ICOL(I,J)
              COL = 0
              DO ROW = 1,KNEQ(I,J)
                COL = MAX(COL,KCOL(KC-1+KLD(KR+ROW)-1))
              END DO

              NCOL(J) = MAX(NCOL(J),COL)
              
C             Add the number of entries in that matrix to NA.
C             Then at the end, NA gives the number of elements in the
C             global matrix.

              NA = NA + KLD(KR+KNEQ(I,J)+1-1)-1
              
            END IF
            
          END IF ! KNEQ <> 0
          
        END DO ! I
        
C       NEQ1 holds the minimum NEQ, NEQ2 the maximum; they must not
C       differ!        
        
        IF (NEQ1.NE.NEQ2) THEN
          NEQ = -1
          RETURN
        END IF
        
C       Store the NEQ in NROW and add it to NEQ

        NROW(I) = NEQ2
        NEQ = NEQ + NEQ2
        
      END DO
      
C     Ok, we now know how many rows and columns belong to each
C     matrix block. Allocate memory for the matrices.
C     We need NEQ+1 entries for the new KLD:

      CALL ZNEW (NEQ+1,-3,LLD,'KLDNEW')
      
C     Allocate memory for KCOL and DA. The number of entries was
C     calculated into NA:

      CALL ZNEW (NA,-3,LCOL,'LCLNEW')
      CALL ZNEW (NA,-1,LA,'LANEW ')
      
C     MCOL and MA point into the global matrix where to continue
C     writing:
      
      MCOL = L(LCOL)
      MA   = L(LA)
      
C     Create the KLD-entry of the first line
      
      KWORK(L(LLD)) = 1
      
C     In GROW, save the number of the current global row:

      GROW = 1
      
C     Now loop through all matrix rows:

      DO I=1,NYBLOC
        
C       In each matrix row, loop through all sub-rows:

        DO ROW=1,NROW(I)

C         RLEN counts the number of elements in the current row.
C         COL counts the number of real columns before our current
C         matrix block.
          
          RLEN = 0
          COL  = 0
          
C         In each sub-row, loop through all *existing* matrix blocks:
          
          DO J=1,NXBLOC
          
            IF ((IA(I,J).NE.0).AND.(ICOL(I,J).NE.0).AND.(ILD(I,J).NE.0)
     *          .AND.(KNEQ(I,J).NE.0)) THEN

C             Where is the row ROW of the matrix block (I,J)?
C             How long is it?
      
              KR = ILD(I,J)-1+ROW
              KC = ICOL(I,J)-1+KLD(KR)
              KA = IA(I,J)-1+KLD(KR)
              PLEN = KLD(KR+1)-KLD(KR)
              
C             Append the row to the new matrix. For standard matrices,
C             we can use standard COPY commands. If sorting is activated
C             and we have a structure 7 matrix, we have to corrent by 
C             hand the position of the diagonal element.

              CALL LCP1(DA(KA),DWORK(MA+RLEN),PLEN)
              CALL LCP3(KCOL(KC),KWORK(MCOL+RLEN),PLEN)
              
              IF (BSORT.AND.
     *            ((PLEN.GT.1).AND.(KCOL(KC).GT.KCOL(KC+1)))) THEN
      
C               This is most likely a structure 7 matrix - at least the
C               first element is at the wrong position. Grab the first 
C               element:

                DAUX = DWORK(MA)
                KAUX = KWORK(MCOL)

C               Loop through each column in this row.
C               Shift every entry until the diagonal is reached.

                DO CCOL = 1, PLEN-1
        
C                 Check if we reached the position of the diagonal 
C                 entry...
                
                  IF (KWORK(MCOL+RLEN+CCOL).GT.KAUX) GOTO 10
        
                  KWORK(MCOL+RLEN+CCOL-1) = KWORK(MCOL+RLEN+CCOL)
                  DWORK(MA+RLEN+CCOL-1)   = DWORK(MA+RLEN+CCOL)
        
                END DO
                
10              CONTINUE

C               Insert the diagonal element:

                KWORK(MCOL+RLEN+CCOL-1) = KAUX
                DWORK(MA+RLEN+CCOL-1)   = DAUX

              END IF
              
C             Increase the column numbers according ro the column
C             we are in. That number is counted in COL.

              DO P=MCOL+RLEN,MCOL+RLEN+PLEN-1
                KWORK(P) = KWORK(P)+COL
              END DO
              
C             Increase the number of elements in the global row
C             by the number of elements we just appended:

              RLEN = RLEN + PLEN
              
            END IF

C           Increase COL to point to the next "real" column number:

            COL = COL + NCOL(J)
            
          END DO ! J
          
C         The row is finished. Create the KLD-pointer for the
C         next row or the end of the matrix, resp.

          KWORK(L(LLD)+GROW+1-1) = KWORK(L(LLD)+GROW-1) + RLEN
          
C         Increase MA and MCOL to point where to write next:

          MA   = MA + RLEN
          MCOL = MCOL + RLEN
          
C         Finally increment the "global row counter"

          GROW = GROW + 1
          
        END DO ! ROW
        
      END DO ! I
      
      END
      
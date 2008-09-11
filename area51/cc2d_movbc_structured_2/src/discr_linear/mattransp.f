************************************************************************
* This file contains routines to modify structure-7 matrices
************************************************************************

************************************************************************
* Transpose matrix structure
*
* This routine accepts the structure of a structure-7 or structure-9
* matrix and creates the structure of the transposed matrix from it.
*
* The resulting structure is not pivoted anymore with the diagonal
* entry in front if a structure-7 matrix is to be transposed!
*
* In:
*   NROW   : Number of rows in the matrix
*   NCOL   : Number of columns in the matrix
*   KCOL   : array [1..NA] of integer
*            Column structure of the matrix
*   KLD    : array [1..NROW+1] of integer
*            Row structure of the matrix
*   KTMP   : array [1..NCOL] of integer
*            Auxiliary array
*
* Out:
*   KCOLD  : array [1..NA] of integer
*            Will be filled with the column structure of the transposed
*            matrix. The array must have the same size as KCOL.
*   KLDD   : array [1..NCOL+1] of integer
*            Will be filled with the row structure of the transposed
*            matrix. The array must have the same size as KCOL.
************************************************************************

      SUBROUTINE TRST79 (NROW,NCOL,KCOL,KLD,KTMP,KCOLD,KLDD)
      
      IMPLICIT NONE
      
C     parameter
      
      INTEGER NCOL,NROW
      INTEGER KCOL(*),KLD(*),KCOLD(*),KLDD(*),KTMP(*)
      
C     local variables

      INTEGER NA,I,J,KP
      
C     size of the matrix?      

      NA = KLD(NROW+1) - 1
      
C     At first determine the number of elements in each column.
C     For this purpose, build the KLDD array.
C
C     Count how many entries <> 0 are in each column. Note this
C     into the KLDD array.
C     This requires one loop through the matrix structure. The
C     corresponding number is written into KLDD, shifted by 1
C     (i.e. KLDD(2) is the number of entries of the 1st column).
C     This helps to create the real KLDD more easily later.

      CALL LCL3 (KLDD,NCOL+1)
      DO I = 1,NA
        KLDD(KCOL(I)+1) = KLDD(KCOL(I)+1) + 1
      END DO
      
C     Now build the real KLDD. this consists of indices, where 
C     each row starts. Row 1 starts at position 1:

      KLDD (1) = 1
      
C     Adding the number of entries KLLD(i) in the row to
C     KLLD(i-1) gives the new row pointer of row i of the transposed
C     matrix.

      DO I = 2,NCOL+1
        KLDD(I) = KLDD(I) + KLDD(I-1)
      END DO

C     That's it for KLDD. Now, KCOLD must be created.
C     This requires another loop through the matrix structure.
C     KTMP receives the index how many entries have been written
C     to each column.

      CALL LCL3 (KTMP,NCOL)
      
      DO I = 1,NROW
      
        DO J = 1,KLD (I+1)-KLD(I)
        
C         Get the column of the item in question:
        
          KP = KCOL(KLD(I)+J-1)          
          
C         Note the right column number in KCOLD:
          
          KCOLD (KLDD(KP)+KTMP(KP)) = I    
          
C         Increment running index of that row

          KTMP (KP) = KTMP (KP) + 1      
          
        END DO
        
      END DO
      
C     KCOLD also finished, that's it.

      END

************************************************************************
* Transpose matrix 
*
* This routine accepts a matrix in structure-7 or structure-9 and
* builds the transposed matrix from it.
*
* The resulting matrix is not pivoted anymore with the diagonal
* entry in front if a structure-7 matrix is to be transposed!
*
* In:
*   NROW   : Number of rows in the matrix
*   NCOL   : Number of columns in the matrix
*   DA     : array [1..NA] of double
*            The matrix entries
*   KCOL   : array [1..NA] of integer
*            Column structure of the matrix
*   KLD    : array [1..NROW+1] of integer
*            Row structure of the matrix
*   KTMP   : array [1..NCOL] of integer
*            Auxiliary array
*
* Out:
*   DT     : array [1..NA] of double
*            The entries of the transposed matrix.
*   KCOLD  : array [1..NA] of integer
*            Will be filled with the column structure of the transposed
*            matrix. The array must have the same size as KCOL.
*   KLDD   : array [1..NCOL+1] of integer
*            Will be filled with the row structure of the transposed
*            matrix. The array must have the same size as KCOL.
************************************************************************

      SUBROUTINE TRM79 (NROW,NCOL,DA,KCOL,KLD,KTMP,DT,KCOLD,KLDD)
      
      IMPLICIT NONE
      
C     parameter
      
      INTEGER NCOL,NROW
      INTEGER KCOL(*),KLD(*),KCOLD(*),KLDD(*),KTMP(*)
      DOUBLE PRECISION DA(*),DT(*)
      
C     local variables

      INTEGER NA,I,J,KP
      
C     size of the matrix?      

      NA = KLD(NROW+1) - 1
      
C     At first determine the number of elements in each column.
C     For this purpose, build the KLDD array.
C
C     Count how many entries <> 0 are in each column. Note this
C     into the KLDD array.
C     This requires one loop through the matrix structure. The
C     corresponding number is written into KLDD, shifted by 1
C     (i.e. KLDD(2) is the number of entries of the 1st column).
C     This helps to create the real KLDD more easily later.

      CALL LCL3 (KLDD,NCOL+1)
      DO I = 1,NA
        KLDD(KCOL(I)+1) = KLDD(KCOL(I)+1) + 1
      END DO
      
C     Now build the real KLDD. this consists of indices, where 
C     each row starts. Row 1 starts at position 1:

      KLDD (1) = 1
      
C     Adding the number of entries KLLD(i) in the row to
C     KLLD(i-1) gives the new row pointer of row i of the transposed
C     matrix.

      DO I = 2,NCOL+1
        KLDD(I) = KLDD(I) + KLDD(I-1)
      END DO

C     That's it for KLDD. Now, KCOLD must be created.
C     This requires another loop through the matrix structure.
C     KTMP receives the index how many entries have been written
C     to each column.

      CALL LCL3 (KTMP,NCOL)
      
      DO I = 1,NROW
      
        DO J = 1,KLD (I+1)-KLD(I)
        
C         Get the column of the item in question:
        
          KP = KCOL(KLD(I)+J-1)          
          
C         Note the right column nomber in KCOLD:
          
          KCOLD (KLDD(KP)+KTMP(KP)) = I    
          
C         Copy the matrix entry:
          
          DT (KLDD(KP)+KTMP(KP)) = DA (KLD(I)+J-1) 
          
C         Increment running index of that row

          KTMP (KP) = KTMP (KP) + 1      
          
        END DO
        
      END DO
      
C     KCOLD and DT are also finished, that's it.

      END

************************************************************************
* Pivot matrix
*
* This routine repivots a matrix in structure 7: The diagonal element
* of each row is moved to the front.
*
* In:
*   NEQ    : Number of rows in the matrix
*   DA     : array [1..NA] of double
*            The matrix entries
*   KCOL   : array [1..NA] of integer
*            Column structure of the matrix
*   KLD    : array [1..NROW+1] of integer
*            Row structure of the matrix
*
* Out:
*   DA     : array [1..NA] of double
*            The matrix entries, pivoted
*   KCOL   : array [1..NA] of integer
*            Column structure of the matrix, pivoted
*   ERR    : Error indicator
*            =0: no error
*            =1: diagonal element not found
************************************************************************

      SUBROUTINE MT7PIV (NEQ,DA,KCOL,KLD)
      
      IMPLICIT NONE
      
C     parameters
      
      INTEGER NEQ,KCOL(*),KLD(*),ERR
      DOUBLE PRECISION DA(*)
      
C     local variables
      
      INTEGER I,J,K,KAUX
      DOUBLE PRECISION DAUX
      
      ERR = 0
      
C     Loop through the rows

      DO I=1,NEQ
      
C       Is the line already pivoted? Most unlikely, but we test for sure

        IF (KCOL(KLD(I)).NE.I) THEN
        
C         Find the position of the diagonal element

          DO J = KLD(I),KLD(I+1)-1
            IF (J.EQ.I) GOTO 10
          END DO
          
C         Oops, diagonal element not found - cancel

          ERR = 1
          RETURN

10        CONTINUE

C         Shift the line by 1. Put the diagonal element to the front.

          DAUX = DA(J)
          KAUX = KCOL(J)
          DO K = J-1,KLD(I),-1
            DA(K+1)   = DA(K)
            KCOL(K+1) = KCOL(K)
          END DO
          DA(KLD(I))   = DAUX
          KCOL(KLD(I)) = KAUX

        END IF
      
      END DO ! I
      
      END
      
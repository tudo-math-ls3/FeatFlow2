************************************************************************
* CSRSRT
*
* Purpose: Sorts the entries in each row of the matrix in ascending
*          order. This is necessary for the UMFPACK4 solver, which
*          expects a matrix in this format.
*          The input matrix is assumed to be in storage technique 7
*          with the first element in each row to be the diagonal
*          element!
*
* Parameters:
*  DA, KCOL, KLD -    Input matrix to be resorted
*  NEQ           -    Dimension of the matrix
*
* Output
*  DA, KCOL      -    resorted matrix entries/column numbers
************************************************************************

      SUBROUTINE CSRSRT (DA, KCOL, KLD, NEQ)

        IMPLICIT NONE

C Parameters

        DOUBLE PRECISION DA(*)
        INTEGER KCOL(*), KLD(*),NEQ
        
C local variables

        DOUBLE PRECISION AUX
        INTEGER I,J

C loop through each row

        DO I = 1, NEQ
        
C Take the diagonal element
          AUX = DA(KLD(I))
          
C Loop through each column in this row.
C Shift every entry until the diagonal is reached.

          DO J = KLD(i)+1, KLD(i+1)-1

C Check if we reached the position of the diagonal entry...
          
            IF (KCOL(J).GT.I) GOTO 10

            KCOL(J-1) = KCOL(J)
            DA(J-1) = DA(J)

          END DO

C If we have reached the diagonal, we can stop and save our diagonal
C entry from the first position there. The rest of the line
C is in ascending order according to the specifications of storage
C technique 7.

10        KCOL(J-1) = I
          DA(J-1) = AUX

       END DO
          
       END
       
       
************************************************************************
* M7IDSH
*
* Purpose: Performs an index shift "-1" on all entries on the KCOL and
*          KLD array. This is due to convert these arrays to 0-based
*          for the use with UMFPACK4
************************************************************************

      SUBROUTINE M7IDSH (KCOL,KLD,NEQ)
      
      IMPLICIT NONE
      
      INTEGER KCOL(*), KLD(*), NEQ
      
      INTEGER I,NA
      
      NA = KLD(NEQ+1)-1
      
C      print *,'kcol'
      DO I=1,NA
        KCOL (I) = KCOL(I)-1
C        print *,kcol(i)
      END DO
      
C      print *,'kld'
      DO I=1,NEQ+1
        KLD(I) = KLD(I)-1
C        print *,kld(i)
      END DO

      END
      
************************************************************************
* M7IDSB
*
* Purpose: Performs an index shift "+1" on all entries on the KCOL and
*          KLD array. Can be used for debugging purposes to transform
*          a matrix back after being modified by M7IDSH.
************************************************************************

      SUBROUTINE M7IDSB (KCOL,KLD,NEQ)
      
      IMPLICIT NONE
      
      INTEGER KCOL(*), KLD(*), NEQ
      
      INTEGER I,NA
      
      NA = KLD(NEQ+1)
      
C      print *,'kcol'
      DO I=1,NA
        KCOL (I) = KCOL(I)+1
C        print *,kcol(i)
      END DO
      
C      print *,'kld'
      DO I=1,NEQ+1
        KLD(I) = KLD(I)+1
C        print *,kld(i)
      END DO

      END
      
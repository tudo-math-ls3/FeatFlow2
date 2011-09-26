************************************************************************
* ZLCHCK
*
* Purpose: L-Array check procedure for debugging.
*
* Parameters:
*  MODE   - 0=Store current content of L-array.
*           1=Check if content of L-array has changed.
*           2=Print current status of L-array.
*  INFO   - Debug information parameter for output
* 
* When checking the content of the L-array, if the content is different
* from the current status, the program is immediately halted!
* The content of the L-array and the desired/saved status are
* printed on screen, together with the user specified INFO parameter, 
* which allows the user to find the code fragment where the error has 
* occurred.
************************************************************************

      SUBROUTINE ZLCHCK (MODE,INFO)
      IMPLICIT NONE
      
      INTEGER NNARR
      PARAMETER (NNARR=299)
      
      COMMON NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      INTEGER NWORK, IWORK, IWMAX, L
      DOUBLE PRECISION DWORK

      INTEGER MODE, INFO, I, J

      INTEGER LSAV (NNARR)
      SAVE LSAV
      
      IF (MODE.EQ.0) THEN
        DO I=1,NNARR
          LSAV (I) = L(I)
        END DO
      ELSE IF (MODE.EQ.1) THEN 
        DO I=1,NNARR
          IF (LSAV(I).NE.L(I)) THEN 
            PRINT *,'L-Array accidently changed!!! Program halted !!!'
            PRINT *,'Debug-position: ',INFO
            PRINT *,'Saved/desired L-status:'
            PRINT *,(LSAV(J),J=1,NNARR)
            PRINT *,'Current L-status:'
            PRINT *,(L(J),J=1,NNARR)
            CALL EXIT(0)
          END IF
        END DO
      ELSE IF (MODE.EQ.2) THEN
        PRINT *,'Debug-position: ',INFO
        PRINT *,'Current L-status:'
        PRINT *,(L(J),J=1,NNARR)
      END IF
      
999   END

************************************************************************
* ZCHKSM
*
* Purpose: Calculate a simple checksum of an array.
*
* Parameters:
*  DX     - An arbitrary array - is treated as INTEGER array
*  LEN    - Length of the array, if treated as INTEGER's.
*           (If the array is REAL or INTEGER, this is the actual
*            length of the array. If the array is DOUBLE, LEN must
*            be 2*length of the actual array)
*  NAME   - Name of the array - for debugging purposes
*  ID     - ID of the array - for debugging purposes
* 
* The name+ID of the array is printed on the screen, together with
* its checksum (all values added). This routine is designed to help
* the user to debug a program, so it's easier to check if the content
* of an array has accidently changed.
************************************************************************

      SUBROUTINE ZCHKSM (DX,LEN,NAME,ID)
      IMPLICIT NONE
      INTEGER LEN, DX
      DIMENSION DX(LEN)
      INTEGER CS,I,ID
      CHARACTER NAME*(10)
      
      CS = 0
      
      DO I=1,LEN
        CS = CS+DX(I)
      END DO
      
      WRITE (*,*)
      WRITE (*,'(A,A10,A,I3,A,I12)') 'Checksum for array ',NAME,' (ID:',
     *                               ID,'): ',CS 
      
      END

      
************************************************************************
* Read output parameters
*
* Read the parameters from the .DAT file concerning the OUTPUT
* parameters. This is necessary on program start to prepare all the
* output. It will open the file CFNAME, read the parameters and
* close the file. 
*
* In:
*   MDATA  - Unit number to use for reading process
*   CFNAME - Name of the file
*   MFILE1 - A unit filename for output to a file
*
* Out:
*   MSHOW  - Message level for initialisation routines (from file).
*            Controls if the initialisation routines print out the
*            readed parameters to screen
*
* Modifies the parameters in the /OUTPUT/ and /FILOUT/ COMMON blocks
* according to the .DAT file. Opens the file unit MFILE for user output
* to a file and stores it globally into the /FILOUT/ COMMON block.
* The filename that should be used for file output is taken from the
* DAT file. 
************************************************************************
      
      SUBROUTINE RDOUT (MDATA,CFNAME,MFILE1,MSHOW)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      INCLUDE 'cerr.inc'
      
C     parameters
      
      CHARACTER CFNAME*(*)
      INTEGER MSHOW,MDATA,MFILE1
      
C     local variables

      INTEGER I,IFMTS

C     Open DAT file
      
      IFMTS = 1
      CALL OF0 (MDATA,CFNAME,IFMTS)
      
C     Ignore the first 6 lines

      DO I=1,6
        READ (MDATA,*)
      END DO
      
C     Check version number

      READ (MDATA,*) I
      IF (I.NE.100) THEN
        WRITE (*,'(A)') 'Error in reading output parameters'
        WRITE (*,'(A)') 'Version number incorrect'
        STOP
      END IF

C     Ignore separator
      
      READ (MDATA,*)
      
C     Read the name of protocol file

      READ (MDATA,*) CFLPRT

C     This data must be read in directly without GETINT.
      
      READ(MDATA,*) M
      M=ABS(M)

      READ(MDATA,*) MT
      MT=ABS(MT)

      READ(MDATA,*) ICHECK
      ICHECK=ABS(ICHECK)

      READ(MDATA,*) MSHOW
      MSHOW=ABS(MSHOW)
      
      IF (MSHOW.GE.2) THEN
        WRITE(MTERM,1) 'Message level for file output:      M      = ',M
        WRITE(MTERM,1) 'Message level for terminal output:  MT     = ',
     *                  MT
        WRITE(MTERM,1) 'Level for tracing:                  ICHECK = ',
     *                  ICHECK
        WRITE(MTERM,1) 'Output level in initialisation:     MSHOW  = ',
     *                  MSHOW
        WRITE(MTERM,*)
      END IF
      
      CLOSE (MDATA)

C     Open file for user output; globally available to all routines

      MFILE = MFILE1
      CALL OF0 (MFILE1,CFLPRT,IFMTS)
      
1     FORMAT (A,I4)
9000  FORMAT(79('-'))

      END


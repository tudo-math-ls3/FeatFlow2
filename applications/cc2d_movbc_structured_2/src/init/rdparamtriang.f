************************************************************************
* Read parametrization / triangulation parameters
*
* This routine reads the parameters about the parametrization and
* triangulation from a given file and stores them in a COMMON block.
*
* In:
*   MDATA  - Unit number to use for reading process
*   MSHOW  - Level of output in the initialization phase
*   CFNAME - Name of the file
************************************************************************

      SUBROUTINE RDPARD (MDATA,MSHOW,CFNAME)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cparametrization.inc'
      
C     parameters
      
      CHARACTER CFNAME*(*)
      INTEGER MSHOW,MDATA
      
C     local variables

      INTEGER I,IFMTS
      CHARACTER CSTR*(255)

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
        WRITE (*,'(A)') 'Error in reading parametr. parameters'
        WRITE (*,'(A)') 'Version number incorrect'
        STOP
      END IF
      
      WRITE (CSTR,'(A)') ' Parametrization/triangulation parameters'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      WRITE (CSTR,'(A)') '------------------------------------------'
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Ignore separator
      
      READ (MDATA,*)

C     Read the data

      READ(MDATA,*) IMESH
      IF (IMESH.NE.1) IMESH=0

      READ(MDATA,*) IRMESH
      IRMESH=ABS(IRMESH)
      
      WRITE (CSTR,1000) 
     *      'Type of parametrization       : IMESH  = ',IMESH
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

      WRITE (CSTR,1000) 
     *      'Type of parametrization       : IMESH  = ',IMESH
      
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      READ(MDATA,*) CPARM
      READ(MDATA,*) CMESH

      WRITE (CSTR,1003) 
     *      'Parametrization file          : CPARM  = ',CPARM
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      WRITE (CSTR,1003) 
     *      'Coarse grid file              : CMESH  = ',CMESH
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     That's it.
      
      WRITE (MTERM,*)
      
      CLOSE (MDATA)
      
1000  FORMAT (A,I4)
1001  FORMAT (A,E16.8)
1002  FORMAT (A,2I4)
1003  FORMAT (A,A)

9000  FORMAT(79('-'))
9001  FORMAT(60('-'))

      END
      
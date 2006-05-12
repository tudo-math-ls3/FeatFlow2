C ********************************************************************
C This file introduces subroutines for easy reading of parameters from
C a DAT file.
C
C Every parameter read from the file is written to screen and/or
C to the output file, depending on MSHOW.
C ********************************************************************
      
C Read INTEGER-value from file IUNIT to INT. Write it to the
C screen/file after the output string NAME.

      SUBROUTINE GETINT(IUNIT,MSHOW,NAME,INTG)
      IMPLICIT NONE
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'

      CHARACTER NAME*(*)
      INTEGER IUNIT,MSHOW
      INTEGER INTG
      CHARACTER CSTR*255
      
      READ(IUNIT,*) INTG
      WRITE(CSTR,'(A,I10)') NAME,INTG
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      END

C Read multiple INTEGER-value from file IUNIT to INT. Write it to the
C screen/file after the output string NAME.
C CNT integers are read as a sequence

      SUBROUTINE GETINM(IUNIT,MSHOW,NAME,INTG,CNT)
      IMPLICIT NONE
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'

      CHARACTER NAME*(*)
      INTEGER IUNIT,CNT,MSHOW
      INTEGER INTG(CNT)
      
      INTEGER J
      
      READ(IUNIT,*) (INTG(J),J=1,CNT)
      IF (MSHOW.GE.0) THEN
        WRITE(MFILE,'(A$)') NAME
        WRITE(MFILE,'(I10$)') (INTG(J),J=1,CNT)
        WRITE(MFILE,*)
      END IF
      IF (MSHOW.GE.2) THEN
        WRITE(MTERM,'(A$)') NAME
        WRITE(MTERM,'(I10$)') (INTG(J),J=1,CNT)
        WRITE(MTERM,*)
      END IF
      
      END

C Read Integer-value from file IUNIT. Interpret it as
C boolean: 0=false, <> 0 = true.
C (So boolean values are treated as integers in the file!)
C
C Write the value to the screen/file after the output string NAME.

      SUBROUTINE GETBOL(IUNIT,MSHOW,NAME,BOL)
      IMPLICIT NONE
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'

      CHARACTER NAME*(*)
      INTEGER IUNIT,MSHOW
      LOGICAL BOL
      
      INTEGER INT
      CHARACTER CSTR*255
      
      READ(IUNIT,*) INT
      BOL = INT.NE.0
      WRITE(CSTR,'(A,I10)') NAME,INT
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      END

C Read DOUBLE-value from file IUNIT to DBL. Write it to the
C screen/file after the output string NAME.

      SUBROUTINE GETDBL(IUNIT,MSHOW,NAME,DBL)
      IMPLICIT NONE
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'

      CHARACTER NAME*(*)
      INTEGER IUNIT,MSHOW
      DOUBLE PRECISION DBL
      CHARACTER CSTR*255
      
      READ(IUNIT,*) DBL
      WRITE(CSTR,'(A,D25.10)') NAME,DBL
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      END

C Read STRING-value from file IUNIT to STR. Write it to the
C screen/file after the output string NAME.

      SUBROUTINE GETSTR(IUNIT,MSHOW,NAME,STR)
      IMPLICIT NONE
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'

      CHARACTER NAME*(*),STR*(*)
      INTEGER IUNIT,MSHOW
      CHARACTER*255 CSTR
      
      READ(IUNIT,*) STR
      WRITE(CSTR,'(A,A)') NAME,STR
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)
      
      END

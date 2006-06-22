************************************************************************
* This file provides additional file handling routines to simplify
* reading/writing of arrays from/to disc.
************************************************************************

      SUBROUTINE ORA1X(MFILE,DX,ILEN,IER,IFMT,ARR)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION DX(*)
      INTEGER MFILE,ILEN,IFMT
      LOGICAL BFMT
      CHARACTER ARR*15
      
      INTEGER ILENBK
      
C     Different handling whether we use formatted or non-formatted
C     input.
      
      BFMT = IFMT.NE.0
      
      ILENBK = ILEN

      IF (BFMT) THEN
        READ (MFILE,'(2A15,2I15)',ERR=99997,END=99997)
     *        ARR,CFORM,ITYPE,ILEN
      ELSE
        READ (MFILE,ERR=99997,END=99997) ARR,ITYPE,ILEN
      ENDIF

C     If the array is longer than allowed, cancel!

      IF (ILEN.GT.ILENBK) THEN
        IER = -1
        RETURN
      END IF
      
C     If the type does not match, cancel

      IF (ITYPE.NE.1) THEN
        IER = -2
        RETURN
      END IF
      
C     Read the data

      IF (BFMT) THEN
        READ (MFILE,CFORM,ERR=99997,END=99997) (DX(I),I=1,ILEN)
      ELSE
        READ (MFILE,ERR=99997,END=99997) (DX(I),I=1,ILEN)
      ENDIF

C     That's it

      IER = 0
      RETURN
      
C     Error

99997 CONTINUE
      IER = -3

      END


      SUBROUTINE OWA1X(MFILE,DX,ILEN,IER,IFMT,CFORM,ARR)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION DX(*)
      INTEGER MFILE,ILEN,IFMT
      LOGICAL BFMT
      CHARACTER CFORM*15,ARR*15
      
      INTEGER ILENBK
      
C     Different handling whether we use formatted or non-formatted
C     input.
      
      BFMT = IFMT.NE.0
      
      ILENBK = ILEN

      IF (BFMT) THEN
        WRITE (MFILE,'(2A15,2I15)',ERR=99997,END=99997)
     *         ARR,CFORM,1,ILEN
      ELSE
        WRITE (MFILE,ERR=99997,END=99997) ARR,1,ILEN
      ENDIF

C     If the array is longer than allowed, cancel!

      IF (ILEN.GT.ILENBK) THEN
        IER = -1
        RETURN
      END IF
      
C     Read the data

      IF (BFMT) THEN
        READ (MFILE,CFORM,ERR=99997,END=99997) (DX(I),I=1,ILEN)
      ELSE
        READ (MFILE,ERR=99997,END=99997) (DX(I),I=1,ILEN)
      ENDIF

C     That's it

      IER = 0
      RETURN
      
C     Error

99997 CONTINUE
      IER = -2

      END

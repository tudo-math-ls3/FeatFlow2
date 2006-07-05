***********************************************************************
* This file contains additional routines for the dynamic
* string handling.
***********************************************************************

***********************************************************************
* Create a DString from a Fortran string
*
* Creates a new DString from a standard Fortran string. Trims
* leading/trailing space characters if desired.
*        string = CSTR   or   string = TRIM(CSTR)
*
* In:
*  CSTR    - The source Fortran string
*  BTRIM   - whether to trim leading+trailing space characters
*
* Out:
*  Return value = handle to the new DString
***********************************************************************

      INTEGER FUNCTION STNEWC (BTRIM,CSTR)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      CHARACTER*(*) CSTR
      LOGICAL BTRIM
      
      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL
      
      INTEGER LTMP
      
C Create an empty DString

      LTMP = STNEWL(LEN(CSTR))
      
C Copy the Fortran string

      CALL ST_REP (1,1,STLEN(LTMP),KWORK(L(LTMP)+1),CSTR)
      
C Trim spaces

      IF (BTRIM) THEN
        CALL STTRM (LTMP)
      END IF
      
      STNEWC = LTMP
      
      END
      
***********************************************************************
* Concatenate with a Fortran string
*
* Appends a Fortran string to a DString with or without trimming:
*        CSTR1 = CSTR1+CSTR2  or   CSTR1 = CSTR1+TRIM(CSTR2)
*
* In:
*  LSTR1   - Handle to the source/destination string where the Fortran
*            string should be appended to
*  CSTR2   - The source Fortran string
*  BTRIM   - whether to trim leading+trailing space characters
*            from the Fortran string
*
* Out:
*  LSTR1   - The modified DString
***********************************************************************

      SUBROUTINE STCATC (LSTR1, BTRIM2, CSTR2)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR1
      LOGICAL BTRIM2
      CHARACTER*(*) CSTR2
      
      INTEGER LTMP
      
      INTEGER STNEWC
      EXTERNAL STNEWC
      
C Create a DString from the Fortran string

      LTMP = STNEWC(BTRIM2,CSTR2)
      
C and concatenate it

      CALL STCAT (LSTR1,LTMP)
      
      CALL STDIS(LTMP)
      
      END

***********************************************************************
* Insert a Fortran string
*
* Inserts a Fortran string to a DString with or without trimming
* at a given position.
*
* In:
*  LSTR1   - Handle to the source/destination DString where the Fortran
*            string should be inserted to
*  CSTR2   - The source Fortran string
*  POS     - The position where the string should be inserted
*  BTRIM   - whether to trim leading+trailing space characters
*            from the Fortran string
*
* Out:
*  LSTR1   - The modified DString
***********************************************************************

      SUBROUTINE STINSC (LSTR, POS, BTRIM2, CSTR2)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR,POS
      LOGICAL BTRIM2
      CHARACTER*(*) CSTR2
      
      INTEGER LTMP
      
      INTEGER STNEWC
      EXTERNAL STNEWC
      
C Create a DString from the Fortran string

      LTMP = STNEWC(BTRIM2,CSTR2)
      
C and insert it

      CALL STINS (LSTR,LTMP,POS)
      
      CALL STDIS(LTMP)
      
      END

***********************************************************************
* Concatenate with an integer
*
* Formats and appends an integer value to a DString:
*        CSTR = CSTR+dstring(I)
*
* In:
*  LSTR   - Handle to the source/destination string where the integer
*           should be appended to
*  I      - The integer value
*  NMINLN - minimum length of the integer value; if the string 
*           representation of I is shorter than NMINLN it will be
*           filled with leading space characters. May be 0.
*  BFIL0  - if .TRUE. the string representation of I is not filled
*           with leading space characters but with leading "0"-
*           characters.
*
* Out:
*  LSTR   - The modified DString
*
* Remark: NMINLN<=32 !
***********************************************************************

      SUBROUTINE STCATI (LSTR, I, NMINLN, BFIL0)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR,NMINLN,I
      LOGICAL BFIL0
      
      INTEGER STNEWC,STNEWL,STLEN
      EXTERNAL STNEWC,STNEWL,STLEN
      CHARACTER STGCH
      EXTERNAL  STGCH
      
      INTEGER LTMP1,LTMP2,NMINL
      LOGICAL BISNEG
      CHARACTER*(32) CTMP
      
C Determine the maximum length of the integer

      NMINL = MIN(NMINLN,32)
      
C Use WRITE to format the integer and create a temporary DString
C of it
      CALL ST_SET (1,32,' ',CTMP)
      WRITE (CTMP,'(I32)') I
      LTMP1 = STNEWC(.TRUE.,CTMP)
      
C Do we have a '-'-sign in front? Delete it, but remember!

      BISNEG = STGCH(LTMP1,1).EQ.'-'
      
C Create an empty DString as long as necessary to fill the string
C representation to NMINL characters

      LTMP2 = STNEWL(MAX(NMINL-STLEN(LTMP1),0))
      
C Fill it with "0" if desired

      IF (BFIL0) THEN

        CALL STSET(LTMP2,1,STLEN(LTMP2),'0')

C The "0" will be put in front of the number. delete it from the
C string and add in in front of the zeroes.

        IF (BISNEG) THEN
          CALL STDEL(LTMP1,1,1)
          CALL STINCH (LTMP2,'-',1)
        END IF

      END IF
      
C Concatenate everything - the leading spaces/zeroes and the
C string representation

      CALL STCAT (LSTR,LTMP2)
      CALL STCAT (LSTR,LTMP1)
      
      CALL STDIS(LTMP2)
      CALL STDIS(LTMP1)
      
      END

***********************************************************************
* Concatenate with an floating point number
*
* Formats and appends a floating point value to a DString:
*        CSTR = CSTR+dstring(FMT,D)
* The floating point value will be formatted with the standard
* Fortran format string FMT, like e.g.:
*        CALL STCATD (LSTR,D,'(D24.16)')
* Leading/trailing space characters are trimmed.
*
* Remark: The format string is allowed only to format the given
* floating point value up to 256 characters! If the caller does not
* take care of that, the stack may be destroyed which might lead
* to a program crash!!!
*
* In:
*  LSTR   - Handle to the source/destination string where the integer
*           should be appended to
*  D      - The floating point value
*  FMT    - The Fortran format string how D should be formatted.
*
* Out:
*  LSTR   - The modified DString
***********************************************************************

      SUBROUTINE STCATD (LSTR, D, FMT)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR
      DOUBLE PRECISION D
      CHARACTER*(*) FMT
      
      INTEGER STNEWC,STLEN
      EXTERNAL STNEWC,STLEN
      
      INTEGER LTMP
      CHARACTER*(256) CTMP
      
C Format the floating point value as desired

      CALL ST_SET (1,256,' ',CTMP)
      WRITE (CTMP,FMT) D
      LTMP = STNEWC(.TRUE.,CTMP)
      
C and concatenate

      CALL STCAT (LSTR,LTMP)
      CALL STDIS(LTMP)
      
      END

***********************************************************************
* Output a string to an output channel
*
* Writes the given string to the output channel MFILE. 
*
* In:
*  LSTR   - Handle to the source string
*  MFILE  - Number of the output channel where to write the string to.
*           If =0, the standard output channel ("*") will be used.
*  BNL    - Whether to break the line / write a newline-character after
*           writing the string.
***********************************************************************

      SUBROUTINE STOUT (LSTR, MFILE, BNL)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR,MFILE
      LOGICAL BNL
      
      INTEGER STLEN,STNEWC,STMID
      EXTERNAL STLEN,STNEWC,STMID
      
      CHARACTER*32 FMT
      CHARACTER*128 STR
      INTEGER LTMP,POS,I,LEN
      
C For a proper output with the WRITE command we have to use
C Fortran strings - which only are allowed to be fixed length.
C We use a buffer of 128 characters and write the string in
C portions of this block length.

C Write the string in blocks a 128 characters:

      POS = 1
      LEN = STLEN(LSTR)
      WRITE (FMT,'(A)') '(A128$)'
      DO I=1,STLEN(LSTR)/128
        LTMP=STMID(LSTR,POS,128)
        CALL STPUT (LTMP,STR)
        CALL STDIS (LTMP)
        IF (MFILE.EQ.0) THEN
          WRITE (*,FMT) STR
        ELSE
          WRITE (MFILE,FMT) STR
        END IF
        POS = POS+128
        LEN = LEN-128
      END DO
      
C Modify the format string for the last < 128 characters
      
      LTMP=STNEWC(.FALSE.,'(A')
      CALL STCATI(LTMP,LEN,0,.FALSE.)
      IF (.NOT.BNL) THEN
        CALL STCATC(LTMP,.FALSE.,'$')
      END IF
      CALL STCATC(LTMP,.FALSE.,')')
      
      CALL STPUT (LTMP,FMT)
      CALL STDIS (LTMP)
      
C and use this to output the rest
      
      LTMP=STMID(LSTR,POS,LEN)
      CALL STPUT (LTMP,STR)
      CALL STDIS (LTMP)
      IF (MFILE.EQ.0) THEN
        WRITE (*,FMT) STR
      ELSE
        WRITE (MFILE,FMT) STR
      END IF
      
      END

***********************************************************************
* Read a line from an input channel
*
* Reads a full line of text from the input channel MFILE. Trailing 
* blank characters are trimmed.
*
* In:
*  MFILE  - Number of the output channel where to write the string to.
*           If =0, the standard input channel ("*") will be used.
*
* Out:
*  LSTR   - Handle to a new string containing the line
*           =0: critical error while reading, no characters could
*               be read
*  IERR   - Error indicator.
*           =0: no error
*           <0: end of file reached; LSTR contains all characters
*               that could be read, or 0 if reading was not possible
*           >0: other (Fortran-specific) error; LSTR contains all 
*               characters that could be read, or 0 if reading was 
*               not possible
***********************************************************************

      SUBROUTINE STINL (LSTR, MFILE, IERR)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
C Size of the temp-array

      INTEGER BSIZE
      PARAMETER (BSIZE=256)
      
      INTEGER LSTR,MFILE,IERR
      
      INTEGER STLEN,STNEWC,STLEF
      EXTERNAL STLEN,STNEWC,STLEF
      
      INTEGER LTMP,LTMP2,LTMP3,P,LEN,SLEN,E
      CHARACTER CTMP*(BSIZE),FCFMT*(30)
      
      INTEGER STNEWL
      EXTERNAL STNEWL
      
C FORTRAN IS ONE OF THE MOST STUPID PROGRAMMING LANGUAGES WHEN
C READING LINES :(((
C
C Ugly tricks must be used. We read in a buffer with BSIZE characters.
C This causes Fortran to jump to the next line. To read the
C next BSIZE characters, we have to backspace to the previous line,
C ignore the first BSIZE characters and read the next block, and 
C so on - until the line is completed.
C
C We assume a line to be completely read if more than the last half
C of the Fortran array containing the line is filled with whitespaces.
C
C In the variable FCFMT we build a Fortran READ format string that
C says how many characters to skip in the front...

C At first create an empty target string and an empty TMP-string.

      LSTR = STNEWL(0)
      SLEN = 0

C Now begin reading the line:

      E = 0

10    CONTINUE

C Build a format string that skips the first SLEN characters:

      LTMP3 = STNEWC (.FALSE.,'(TR')
      CALL STCATI (LTMP3,SLEN,0,.FALSE.)
      CALL STCATC (LTMP3,.FALSE.,',A')
      CALL STCATI (LTMP3,BSIZE,0,.FALSE.)
      CALL STCATC (LTMP3,.FALSE.,')')
      CALL STPUT (LTMP3, FCFMT)
      CALL STDIS (LTMP3)
      
C Read the line

      READ (UNIT=MFILE, FMT=FCFMT, IOSTAT=IERR) CTMP
      
C Has there been an error?
      
      IF (IERR.LE.0) THEN

C No. Create a string from the information readed:
        
        LTMP = STNEWC (.FALSE.,CTMP)

C Trim the trailing spaces 

        CALL ST_TPS (STLEN(LTMP),P,LEN,KWORK(L(LTMP)+1))
        LTMP2 = STLEF (LTMP,LEN+P-1)
        
C If the string contains less than BSIZE/2 characters, we assume that 
C we have reached the end of the line. If not, there may be more
C information following. Then we have to backspace and repeat the read
C starting at the position we stopped.
C Otherwise we can stop reading.

        IF (STLEN(LTMP2).LT.(BSIZE/2)) THEN
          E = 1
        ELSE
          BACKSPACE (MFILE)
        END IF

C Append the full string that was readed

        CALL STCAT (LSTR,LTMP)
        
C Advance BSIZE characters in the source line in the file

        SLEN = SLEN + BSIZE

C Delete the temporary space

        CALL STDIS (LTMP2)
        CALL STDIS (LTMP)

C Stop reading if we have reached the end of the file

        IF (IERR.LT.0) E=1

      ELSE
      
C Some kind of error; stop reading

        E = 1
      
      END IF
 
C go on reading

      IF (E.EQ.0) GOTO 10 
 
C Trim the trailing spaces because it's not sure whether they come from
C the Fortran READ command or from the source file!

      CALL ST_TPS (STLEN(LSTR),P,LEN,KWORK(L(LSTR)+1))
      CALL STSLN (LSTR,LEN+P-1)
 
C If nothing could be read, dispose the LSTR, too, and set it to 0
C to indicate an error

      IF (STLEN(LSTR).EQ.0) THEN
        CALL STDIS (LSTR)
        LSTR = 0
      END IF
      
      END

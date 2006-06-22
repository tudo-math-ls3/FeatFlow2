***********************************************************************
* Dynamic string handling - auxiliary routines
*
* The routines in this file are auxiliary routines for the
* string handling implementation. They work on the low-level
* character arrays without any knowledge about the FEAT memory
* management.
***********************************************************************

***********************************************************************
* Character to upper case
*
* Transforms a character into its uppercase representation
*
* In:
*  C       - The character
*
* Out:
*  Return value = uppercase (C)
***********************************************************************

      CHARACTER FUNCTION TOUPRC (C)
      
      IMPLICIT NONE

      CHARACTER C
      
      IF ((ICHAR(C).GE.97).AND.(ICHAR(C).LE.122)) THEN
        TOUPRC = CHAR(ICHAR(C)-32)
      ELSE
        TOUPRC = C
      END IF
      
      END


***********************************************************************
* Character to lower case
*
* Transforms a character into its lowercase representation
*
* In:
*  C       - The character
*
* Out:
*  Return value = lowercase (C)
***********************************************************************

      CHARACTER FUNCTION TOLWRC (C)
      
      IMPLICIT NONE

      CHARACTER C
      
      IF ((ICHAR(C).GE.65).AND.(ICHAR(C).LE.90)) THEN
        TOLWRC = CHAR(ICHAR(C)+32)
      ELSE
        TOLWRC = C
      END IF
      
      END


***********************************************************************
* General string replacement
*
* Copies characters from string CSTR2 into CSTR1 while overwriting
* the old content:
*        CSTR1(POS1:POS1+LEN) = CSTR2(POS2:POS2+LEN)
*
* In:
*  CSTR1   - The destination string
*  CSTR2   - The source string
*  POS1    - Start position in the destination string
*  POS2    - Start position in the source string
*  LEN     - Length of the partial string that sould be copied
***********************************************************************

      SUBROUTINE ST_REP (POS1,POS2,LEN,CSTR1,CSTR2)
      
      IMPLICIT NONE
      
      CHARACTER CSTR1(*),CSTR2(*)
      INTEGER POS1,POS2,LEN
      
      INTEGER I
      
      DO I=1,LEN
        CSTR1(POS1-1+I)=CSTR2(POS2-1+I)
      END DO
      
      END


***********************************************************************
* String-filling
*
* Fills portions of a string with characters
*        CSTR(POS:POS+LEN) = C
*
* In:
*  CSTR   - The destination string
*  POS    - Start position in the destination string
*  LEN    - Length of the partial string that will be overwritten
*  C      - Character that should be copied into the string
***********************************************************************

      SUBROUTINE ST_SET (POS,LEN,C,CSTR)
      
      IMPLICIT NONE

      CHARACTER CSTR(*)
      INTEGER POS,LEN
      CHARACTER C
      
      INTEGER I
      
      DO I=1,LEN
        CSTR(POS-1+I) = C
      END DO
      
      END


***********************************************************************
* Determine string trim-positions
*
* Searches for the first and last character in a string that is not
* a space character.
*
* In:
*  CSTR   - The source string
*  SLEN   - Length of the source string
* 
* Out:
*  POS    - Position of the first character
*  LEN    - Length of the partial string with last character <> ' '
***********************************************************************

      SUBROUTINE ST_TPS (SLEN,POS,LEN,CSTR)
      
      IMPLICIT NONE

      CHARACTER CSTR(*)
      INTEGER POS,LEN,SLEN
      
      INTEGER I,ENDE
      LOGICAL SFOUND
      
      SFOUND = .FALSE.
      POS = 1
      ENDE = 0
      DO I=1,SLEN
        IF (CSTR(I).NE.' ') THEN
          IF (.NOT.SFOUND) THEN
            POS = I
            ENDE = I
            SFOUND = .TRUE.
          ELSE
            ENDE = I
          END IF
        END IF
      END DO
      LEN = ENDE-POS+1
      
      END

***********************************************************************
* Set character
*
* Replaced a character in the string:
*        CSTR(POS) = C
*
* In:
*  CSTR   - The destination string
*  POS    - Character position to be replaced
*  C      - The character
***********************************************************************

      SUBROUTINE ST_SCH (POS,C,CSTR)
      
      IMPLICIT NONE

      CHARACTER CSTR(*)
      INTEGER POS
      CHARACTER C
      
      CSTR(POS) = C
      
      END

***********************************************************************
* Get character
*
* Returns the character at the desired position:
*
* In:
*  CSTR   - The source string
*  POS    - Character position 
*
* Out:
*  Return value = CSTR(POS)
***********************************************************************

      CHARACTER FUNCTION ST_GCH (POS,CSTR)
      
      IMPLICIT NONE

      CHARACTER CSTR(*)
      INTEGER POS
      
      ST_GCH = CSTR(POS)
      
      END

***********************************************************************
* Case-sensitive string-compare
*
* Compares two strings, case-sensitive
*
* In:
*  CSTR1  - The first string
*  CSTR2  - The second string
*  LEN1   - Length of CSTR1
*  LEN2   - Length of CSTR2
*
* Out:
*  Return value = -1, if CSTR1 < CSTR2
*                  0, if CSTR1 = CSTR2
*                  1, if CSTR1 > CSTR2
***********************************************************************

      INTEGER FUNCTION ST_CMP (LEN1,LEN2,CSTR1,CSTR2)
      
      IMPLICIT NONE

      CHARACTER CSTR1(*),CSTR2(*)
      INTEGER LEN1,LEN2
      
      INTEGER I
      
      DO I=1,MIN(LEN1,LEN2)
        IF (CSTR1(I).NE.CSTR2(I)) THEN
          IF (CSTR1(I).LT.CSTR2(I)) THEN
            ST_CMP = -1
            RETURN
          ELSE
            ST_CMP = 1
            RETURN
          END IF
        END IF
      END DO

      IF (LEN1.LT.LEN2) THEN
        ST_CMP = -1
      ELSE IF (LEN1.GT.LEN2) THEN
        ST_CMP = 1
      ELSE 
        ST_CMP = 0
      END IF
      
      END

***********************************************************************
* Search for first occurance of a character
*
* Searches for a character in a string, from left to right
*
* In:
*  CSTR  - The string
*  LEN   - Length of CSTR
*  C     - The character to search for
*
* Out:
*  Return value = First position in CSTR containing C
*                 0, if CSTR doesn't contain C
***********************************************************************

      INTEGER FUNCTION ST_CHR (LEN,C,CSTR)
      
      IMPLICIT NONE

      CHARACTER CSTR(*)
      INTEGER LEN
      CHARACTER C
      
      INTEGER I
      
      DO I=1,LEN
        IF (CSTR(I).EQ.C) THEN
          ST_CHR = I
          RETURN
        END IF
      END DO
      
      ST_CHR = 0
      
      END

***********************************************************************
* Search for last occurance of a character
*
* Searches for a character in a string, from right to left
*
* In:
*  CSTR  - The string
*  LEN   - Length of CSTR
*  C     - The character to search for
*
* Out:
*  Return value = Last position in CSTR containing C
*                 0, if CSTR doesn't contain C
***********************************************************************

      INTEGER FUNCTION ST_RCH (LEN,C,CSTR)
      
      IMPLICIT NONE

      CHARACTER CSTR(*)
      CHARACTER C
      INTEGER LEN
      
      INTEGER I
      
      DO I=LEN,1,-1
        IF (CSTR(I).EQ.C) THEN
          ST_RCH = I
          RETURN
        END IF
      END DO
      
      ST_RCH = 0
      
      END

***********************************************************************
* Search for first occurance of a substring in a string
*
* Tests if CSTR1 contains the substring CSTR2.
*
* In:
*  CSTR1  - The string
*  CSTR2  - The substring to search for
*  LEN1   - Length of CSTR1
*  LEN2   - Length of CSTR2
*  POS    - Position where to start the search
*
* Out:
*  Return value = First position in CSTR where CSTR2 starts
*                 0, if CSTR1 doesn't contain CSTR2
***********************************************************************

      INTEGER FUNCTION ST_STR (LEN1,LEN2,POS,CSTR1,CSTR2)
      
      IMPLICIT NONE

      CHARACTER CSTR1(*),CSTR2(*)
      INTEGER LEN1,LEN2,POS
      
      INTEGER I,CPOS
      
      ST_STR = 0
      
      DO CPOS = POS,LEN1-LEN2+1
        DO I=1,LEN2
          IF (CSTR1(CPOS+I-1).NE.CSTR2(I)) GOTO 1
        END DO
        ST_STR = CPOS
        RETURN
1       CONTINUE
      END DO
      
      END

***********************************************************************
* Search for a word in a string
* 
* Search for a word in a string. The search starts at position SPOS.
*
* In:
*  CSTR  - The string where to search in
*  SLEN  - The length of the string CSTR
*  SPOS  - First position where to start the search
*  C     - separation character for words; should be set to " "
*  BENC  - whether to take care of strings enclosed by '"'-characters
* 
* Out:
*  POS   - Position of the first character <> C or 
*          0, if no character was found
*  LEN   - Length of the word
***********************************************************************

      SUBROUTINE ST_PRS (SLEN,C,BENC,SPOS,POS,LEN,CSTR)
      
      IMPLICIT NONE

      CHARACTER CSTR(*)
      CHARACTER C
      INTEGER SLEN,SPOS,POS,LEN
      LOGICAL BENC
      
      INTEGER I,J
      LOGICAL BESC
      
      BESC = .FALSE.
      
      DO I=SPOS,SLEN
        IF (CSTR(I).NE.C) THEN

C A character was found. Search for the end of the word

          POS = I
          DO J=I,SLEN
          
C As soon we hit a '"'-character, we switch the "escaped"-mode.
C In "escaped"-mode we ignore separation characters!

            IF (BENC) THEN
              IF (CSTR(J).EQ.'"') BESC = .NOT.BESC
            END IF
          
            IF (.NOT.BESC) THEN
              IF (CSTR(J).EQ.C) THEN
C End of the word
                LEN = J-I+1
                RETURN
              END IF
            END IF
            
          END DO

C End of the string

          LEN = SLEN-I+1
          RETURN
          
        END IF
      END DO
      
C Nothing found
      
      POS = 0
      LEN = 0

      END
      
***********************************************************************
* Convert a string to uppercase
*
* Converts the string CSTR to its uppercase representation
*
* In:
*  CSTR   - The string
*  LEN    - Length of CSTR
*
* Out:
*  CSTR   - The string in uppercase letters
***********************************************************************

      SUBROUTINE ST_UPR (LEN,CSTR)
      
      IMPLICIT NONE

      CHARACTER CSTR(*)
      INTEGER LEN
      
      CHARACTER TOUPRC
      EXTERNAL TOUPRC
      
      INTEGER I
      
      DO I=1,LEN
        CSTR(I) = TOUPRC (CSTR(I))
      END DO
      
      END

***********************************************************************
* Convert a string to lowercase
*
* Converts the string CSTR to its lowercase representation
*
* In:
*  CSTR   - The string
*  LEN    - Length of CSTR
*
* Out:
*  CSTR   - The string in lowercase letters
***********************************************************************

      SUBROUTINE ST_LWR (LEN,CSTR)
      
      IMPLICIT NONE

      CHARACTER CSTR(*)
      INTEGER LEN

      CHARACTER TOLWRC
      EXTERNAL TOLWRC
      
      INTEGER I
      
      DO I=1,LEN
        CSTR(I) = TOLWRC (CSTR(I))
      END DO
      
      END

***********************************************************************
* Case-insensitive string-compare
*
* Compares two strings, case-insensitive
*
* In:
*  CSTR1  - The first string
*  CSTR2  - The second string
*  LEN1   - Length of CSTR1
*  LEN2   - Length of CSTR2
*
* Out:
*  Return value = -1, if CSTR1 < CSTR2
*                  0, if CSTR1 = CSTR2
*                  1, if CSTR1 > CSTR2
***********************************************************************

      INTEGER FUNCTION ST_ICM (LEN1,LEN2,CSTR1,CSTR2)
      
      IMPLICIT NONE

      CHARACTER CSTR1(*),CSTR2(*)
      INTEGER LEN1,LEN2
      
      CHARACTER TOUPRC
      EXTERNAL TOUPRC
      
      INTEGER I
      CHARACTER C1,C2
      
      DO I=1,MIN(LEN1,LEN2)
        C1 = TOUPRC (CSTR1(I))
        C2 = TOUPRC (CSTR2(I))
        IF (C1.NE.C2) THEN
          IF (C1.LT.C2) THEN
            ST_ICM = -1
            RETURN
          ELSE
            ST_ICM = 1
            RETURN
          END IF
        END IF
      END DO

      IF (LEN1.LT.LEN2) THEN
        ST_ICM = -1
      ELSE IF (LEN1.GT.LEN2) THEN
        ST_ICM = 1
      ELSE 
        ST_ICM = 0
      END IF
      
      END


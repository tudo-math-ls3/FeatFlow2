***********************************************************************
* This file implements dynamic string handling
*
* The string handling that is used here combines aspects of the string
* handling that is used by Turbo Pascal and C/C++. The names of the
* commands here and the syntax is comparable of C/C++: the "target"
* argument is generally the first one while the other arguments
* are "source" or modification parameters - except where mentioned.
*
* General implementation aspects:
* -------------------------------
* The dynamic strings - we denote them as DSTRINGS in the following -
* are identified by handles and dynamically allocated on the KWORK-
* array of the FEAT-Heap. The implementation follows that of PASCAL-
* strings: Every string consists of a length-information stored in 
* the first integer, while the rest of the allocated memory block 
* stores characters and is accessed as character array. 
*
* Whenever not mentioned, the handles to the "target" string may be
* modified in every string manipulation routine due to re-allocation
* of memory on the FEAT-heap. Therefore, never store a copy of a 
* string-handle anywhere. When you need a backup of a string, use 
* STDUP to create a copy!
* 
* Remark on speed-efficiency
* --------------------------
* One warning: AVOID AT ALL COSTS TO CREATE A STRING AT THE VERY 
* BEGINNING OF THE PROGRAM WHILE YOU NEED TO MODIFY IT AT THE VERY 
* END!!!
* Imagine the following: You create a string at the very beginning.
* Then you allocate more memory on the FEAT-heap to store e.g. matrices
* of 1-2 GB size. After this you trim the string to delete leading/
* trailing space characters. What happens is that the routines will 
* create a new string, copy all non-space characters into that, delete 
* your old string and return the new one. As the FEAT memory management
* immediately fill up the memory hole that arises in that situation by 
* deleting the old string, the whole work-arrays of 1-2 GB size will 
* be moved which might take a while! Therefore it might be very much 
* faster first to create a copy of the string at the end of the heap 
* and then modify that!
*
* Remark on 64-bit machines
* -------------------------
* It's assumed that one INTEGER contains space for 4 characters. 
* Although this is not true on 64-bit machines, there the effect of this 
* treatment is only that the allocated memory blocks are double the 
* size that is necessary for the storing of the actual string, but the 
* implementation should still work.
*
* Structure of the implementation
* -------------------------------
* The string implementation consists roghly of 4 files plus optional
* header files:
*  DSTRINGSAUX.F - Contains the actual low-level character-array 
*    modification routines
*  DSTRINGS1.F   - Contains string modification routines for one or
*    more DStrings like concat, insert, delete,...
*  DSTRINGS2.F   - Contains routines for creating/modification of
*    DStrings with standard Fortran strings and output routines
*  DSTRINGS3.F   - Contains routines for the modification of standard
*    Fortran strings with the help of DStrings
*
* This decomposition was basically introduced to avoid nasty compiler
* warnings because of type-conflicts. THE TYPE CONFLICTS THAT MAY
* ARISE IN THESE ROUTINES ARE INTENTIONAL AND CAME UP DUE TO THE NATURE
* OF THE IMPLEMENTATION! As we use the KWORK-integer-array for
* storing character strings mixed with a length information, we
* frequently have to call low-level modification routines with the
* address of the first character of the string. Obviously this
* might produce a false compiler warning, which can be ignored!
*
* Example of a typical implementation
* -----------------------------------
* Let's assume we want to create a filename of a GMV-file with the
* help of DStrings. We first declare some string-handles and a standard
* Fortran string:
*
*   INTEGER LSTR
*   CHARACTER FNAME*64
*
* IMPORTANT: It's crucial that here FNAME is declared as a Fortrag-
* string with the modifier "*64". Declaring it as character array 
* "FNAME(64)" will not work, then only the first character will
* be affected when trying to obtain the Fortrag-string from the 
* DString!!! (-> Character-Array=List of Fortran-strings with length 1)
*
* At first we create the basic filename from a Fortran string constant:
*
*   LSTR = STNEWC (.FALSE.,'SOL')
*
* Let's assume the Integer variable NUM contains the number of the GMV-
* file. We add this number to the string and fill it with leading "0"
* if the number has less than 4 characters in length:
*
*   CALL STCATI (LSTR,NUM,4,.TRUE.)
*
* Then we add the GMV-extension to the filename:
*
*   CALL STCATC (LSTR,.FALSE.,'.gmv')
*
* Finally we copy the filename to the original Fortran string and
* delete it from the heap:
*
*  CALL STPUT (LSTR, FNAME)
*  CALL STDIS (LSTR)
*
* That's all. If e.g. NUM=20, the Fortran string FNAME now contains
* the string "SOL0020.gmv", filled up with spaces to the end.
*
* List of routines
* ----------------
* Here a brief overview about the routines in this implementation:
*
*  STNEWL - Create new empty string
*  STNEWC - Create a DString from a Fortran string
*  STPUT  - Convert a DString into a Fortran string
*  STOUT  - Output a string to an output channel
*  STINL  - Read a line from a file and write it into a DString
*  STDIS  - Release string
*  
*  STLEN  - Determine the length of the string
*  
*  STCPY  - Copy a string
*  STNCP  - Partially copy a string
*  STDUP  - Duplicate a string
*  STCMP  - String-comparison, case-sensitive
*  STICM  - String-comparison, case-insensitive
*  STNCM  - Partial string-comparison
*  STREP  - General string replacement without length-modification
*  STOWR  - Overwrite a string with automatic length-modification
*  
*  STINS  - Insert one string into another string
*  STINCH - Insert a character
*  STINSC - Insert a Fortran string
*  
*  STDEL  - Deletes portions of a string
*  
*  STCAT  - Append one string to another
*  STCATC - Concatenate with a Fortran string
*  STCATI - Concatenate with an integer
*  STCATD - Concatenate with an floating point number
*  STNCT  - Appends portions of a string to a string
*  
*  STSET  - Fill string with characters
*  STLEF  - Extract the first LEN characters
*  STRIG  - Extract the last LEN characters
*  STMID  - Extract a substring
*  
*  STTPS  - Determine string trim-positions
*  STTRM  - Trim a string
*  STTRP  - Trim a string partially
*  STLJS  - Left-justify string
*  STRJS  - Right-justify string
*  STCJS  - Center-justify string
*  STPRS  - String-parsing
*  
*  STSLN  - Set string length
*  
*  TOUPRC - Convert a character to uppercase
*  TOLWRC - Convert a character to lowercase
*  STSCH  - Set a character
*  STGCH  - Get a character
*  STUPR  - Convert to uppercase
*  STLWR  - Convert to lowercase
*  
*  STCHR  - Search for first occurance of a character
*  STRCH  - Search for last occurance of a character
*  STSTR  - Search for first occurance of a substring in a string
*
* (C) by FeatFlow group, http://www.featflow.de
* Author: Michael Koester, 06.12.2005
***********************************************************************


***********************************************************************
* General string replacement
*
* Copies characters from string CSTR2 into CSTR1 while overwriting
* the old content:
*        CSTR1(POS1:POS1+LEN) = CSTR2(POS2:POS2+LEN)
*
* In:
*  LSTR1   - Handle to the destination string
*  LSTR2   - Handle to the source string
*  POS1    - Start position in the destination string
*  POS2    - Start position in the source string
*  LEN     - Length of the partial string that sould be copied
*
* Out:
*  LSTR1   - Handle to the modified string
*
* If the string CSTR1 is not long enough, the replacement will stop
* at the end of CSTR1 without concerning the rest of CSTR2.
* For string replacement with augmenting CSTR1 use STOWR!
***********************************************************************

      SUBROUTINE STREP (LSTR1,LSTR2,POS1,POS2,LEN)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR1, LSTR2, POS1, POS2, LEN

      INTEGER STLEN
      EXTERNAL STLEN

      IF ((LSTR1.EQ.0).OR.(LSTR2.EQ.0)) RETURN

      CALL ST_REP (MAX(POS1,1+STLEN(LSTR1)),
     *             MAX(POS2,1+STLEN(LSTR2)),
     *             MAX(0,MIN(MIN(LEN,STLEN(LSTR2)-POS2+1),
     *                       STLEN(LSTR1)-POS1+1)),
     *             KWORK(L(LSTR1)+1),KWORK(L(LSTR2)+1) )
      
      END
      
***********************************************************************
* Fill string with characters
*
* Fills portions of a string with characters
*        CSTR(POS:POS+LEN) = C
*
* In:
*  LSTR   - Handle to the destination string
*  POS    - Start position in the destination string
*  LEN    - Length of the partial string that will be overwritten
*  C      - Character that should be copied into the string
*
* If POS=LEN=0, the whole string will be overwritten.
***********************************************************************

      SUBROUTINE STSET (LSTR,POS,LEN,C)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR,POS,LEN
      CHARACTER C
      
      INTEGER STLEN
      EXTERNAL STLEN
      
      INTEGER CPOS,CLEN
      
C Parameter corrections
      
      CPOS = POS
      CLEN = LEN
      
      IF ((CPOS.EQ.0).OR.(CLEN.EQ.0)) THEN
        CPOS = 1
        CLEN = STLEN(LSTR)
      END IF
      
      IF ((LSTR.EQ.0).OR.(CPOS.EQ.0).OR.
     *    (POS.GT.STLEN(LSTR).OR.(CLEN.LE.0))) RETURN

C Determine where to overwrite what

      IF (POS.GE.1) THEN
        CALL ST_SET (CPOS,MIN(STLEN(LSTR)-CPOS+1,MAX(0,LEN)),
     *               C,KWORK(L(LSTR)+1))
      ELSE
        CALL ST_SET (1,MIN(STLEN(LSTR),MAX(0,LEN+POS-1)),
     *               C,KWORK(L(LSTR)+1))
      END IF
      
      END

***********************************************************************
* Create new empty string
*
* Creates a new string on the heap containing LENGTH space characters.
*
* In:
*  LENGTH - The desired length of the string. 
*           =0 is allowed, will cerate an empty string.
*
* Out:
*  Return value = handle to new string
***********************************************************************

      INTEGER FUNCTION STNEWL (LENGTH)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LENGTH
      INTEGER LTMP
      
C Allocate memory: 4 characters per integer + 1 integer for the length
      
      LTMP = 0
      CALL ZNEW (1+(MAX(LENGTH,0)+3)/4,-3,LTMP,'STRING')
      STNEWL = LTMP
      
C Store the length in the first integer
      
      KWORK(L(LTMP)) = MAX(LENGTH,0)
      
C Fill the string with spaces
      
      CALL STSET (LTMP,1,LENGTH,' ')
      
      END

***********************************************************************
* Release string
*
* Deletes a string, release the memory.
*
* In:
*  LSTR   - Handle to the string to be deleted
*
* Out:
*  LSTR = 0
***********************************************************************

      SUBROUTINE STDIS (LSTR)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR
      
      CALL ZDISP (0,LSTR,'STRING')
      
      END
      
***********************************************************************
* Determine the length of the string
*
* Deletes a string, release the memory.
*
* In:
*  LSTR   - Handle to the string 
*
* Out:
*  Return value = Length of the string
***********************************************************************
      
      INTEGER FUNCTION STLEN (LSTR)

      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR
      
C The length can be found directly in the first integer
      
      IF (LSTR.EQ.0) THEN
        STLEN = 0
      ELSE
        STLEN = KWORK(L(LSTR))
      END IF
      
      END
      
***********************************************************************
* Append one string to another
*
* Appends CSTR2 to CSTR1 overwriting the old CSTR1.
*        CSTR1 = CSTR1+CSTR2
*
* In:
*  LSTR1   - Handle to the first string
*            if =0, the string is assumed to be empty;
*            the result will be a copy of CSTR2.
*  LSTR2   - Handle to the second string
*
* Out:
*  LSTR1   - Handle to the concatenated string.
***********************************************************************
      
      SUBROUTINE STCAT (LSTR1, LSTR2)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR1, LSTR2
      INTEGER LTMP
      
      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL

C no second string -> nothing to concatenate

      IF (LSTR2.EQ.0) RETURN
      
C If there's no first string, create an empty new one

      IF (LSTR1.EQ.0) THEN
        LSTR1 = STNEWL(0)
      END IF
      
C Create a string large enough to hold the concatenated result
      
      LTMP = STNEWL(STLEN(LSTR1)+STLEN(LSTR2))
      
C Copy both source strings into it
      
      CALL ST_REP (1,1,STLEN(LSTR1),KWORK(L(LTMP)+1),KWORK(L(LSTR1)+1))
      CALL ST_REP (1+STLEN(LSTR1),1,STLEN(LSTR2),
     *             KWORK(L(LTMP)+1),KWORK(L(LSTR2)+1))
     
C and replace LSTR1 with the result
     
      CALL STDIS (LSTR1)
      LSTR1 = LTMP
      
      END

***********************************************************************
* Appends portions of a string to a string
*
* Appends the dirst LEN characters of CSTR2 to CSTR1 overwriting the 
* old CSTR1.
*        CSTR1 = CSTR1+CSTR2(1:LEN)
*
* In:
*  LSTR1   - Handle to the first string
*            if =0, the string is assumed to be empty;
*            the result will be a copy of CSTR2.
*  LSTR2   - Handle to the second string
*  LEN     - Length of the substring to be appended to CSTR1
*
* Out:
*  LSTR1   - Handle to the concatenated string.
***********************************************************************
      
      SUBROUTINE STNCT (LSTR1, LSTR2, LEN)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR1, LSTR2, LEN
      INTEGER LTMP
      
      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL

      INTEGER CLEN

C no second string -> nothing to concatenate

      IF (LSTR2.EQ.0) RETURN
      
C If there's no first string, create an empty new one

      IF (LSTR1.EQ.0) THEN
        LSTR1 = STNEWL(0)
      END IF
      
C Determine the allowed length

      CLEN = MAX(0,MIN(LEN,STLEN(LSTR2)))
      
C Create a string large enough to hold the concatenated result
      
      LTMP = STNEWL(STLEN(LSTR1)+CLEN)
      
C Copy both source strings into it
      
      CALL ST_REP (1,1,STLEN(LSTR1),KWORK(L(LTMP)+1),KWORK(L(LSTR1)+1))
      CALL ST_REP (1+STLEN(LSTR1),1,CLEN,
     *             KWORK(L(LTMP)+1),KWORK(L(LSTR2)+1))
     
C and replace LSTR1 with the result
     
      CALL STDIS (LSTR1)
      LSTR1 = LTMP
      
      END

***********************************************************************
* Copy a string
*
* Copies CSTR2 to CSTR1 overwriting the old CSTR1 completely:
*        CSTR1 = CSTR2
*
* In:
*  LSTR1   - Handle to the first string;
*            if =0, the string is assumed to be empty;
*            the result will be a copy of CSTR2.
*  LSTR2   - Handle to the second string. 
*
* Out:
*  LSTR1   - Handle to the copy of CSTR2.
***********************************************************************
      
      SUBROUTINE STCPY (LSTR1,LSTR2)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR1, LSTR2

      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL

C no second string -> nothing to copy

      IF (LSTR2.EQ.0) RETURN
      
      IF (LSTR1.EQ.0) THEN
C If there's no first string, create an empty new one
        LSTR1 = STNEWL (STLEN(LSTR2))
      ELSE
C otherwise set the correct length
        CALL STSLN (LSTR1,STLEN(LSTR2))
      END IF
      
C Copy the content

      CALL ST_REP (1,1,STLEN(LSTR2),KWORK(L(LSTR1)+1),KWORK(L(LSTR2)+1))
      
      END

***********************************************************************
* Partially copy a string
*
* Copies the first LEN characters from CSTR2 to CSTR1 overwriting 
* the old CSTR1. If CSTR2 is longer than CSTR1, CSTR1 will be replaced
* completely. If CSTR2 is shorter than CSTR1, the first LEN characters
* are replaced, the rest of CSTR1 is left as is.
*        CSTR1(1:LEN) = CSTR2(1:LEN)
*
* In:
*  LSTR1   - Handle to the first string;
*            if =0, the string is assumed to be empty;
*            the result will be a copy of CSTR2.
*  LSTR2   - Handle to the second string. 
*
* Out:
*  LSTR1   - Handle to the copy of CSTR2(1:LEN)
***********************************************************************
      
      SUBROUTINE STNCP (LSTR1,LSTR2,LEN)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR1, LSTR2, LEN

      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL
      
      INTEGER CLEN

C no second string -> nothing to copy

      IF (LSTR2.EQ.0) RETURN
      
C determine the allowed length

      CLEN = MAX(0,MIN(LEN,STLEN(LSTR2)))
      
      IF (LSTR1.EQ.0) THEN
C If there's no first string, create an empty new one
        LSTR1 = STNEWL (CLEN)
      ELSE
C Otherwise set the correct length. Don't shorten the old string!
        CALL STSLN (LSTR1,MAX(CLEN,STLEN(LSTR1)))
      END IF
      
C Copy the content

      CALL ST_REP (1,1,CLEN,KWORK(L(LSTR1)+1),KWORK(L(LSTR2)+1))
      
      END

***********************************************************************
* Duplicate a string
*
* Create a copy of the string CSTR and return a handle to it
*
* In:
*  LSTR   - Handle to the source string.
*
* Out:
*  Return value = handle to copy of CSTR.
***********************************************************************

      INTEGER FUNCTION STDUP (LSTR)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR

      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL
      
      INTEGER LTMP
      
C Works also if LSTR because STLEN returns 0, so a new
C empty string will be allocated
      
      LTMP = STNEWL (STLEN(LSTR))
      CALL ST_REP (1,1,STLEN(LSTR),
     *             KWORK(L(LTMP)+1),KWORK(L(LSTR)+1))
      STDUP = LTMP
      
      END

***********************************************************************
* Deletes portions of a string
*
* Deletes the substring CSTR(POS:POS+LEN) from the string CSTR.
*
* In:
*  LSTR   - Handle to the source string
*  POS    - Position of the substring
*  LEN    - Length of the substring to be deleted
*
* Out:
*  LSTR   - Handle to the modified string
***********************************************************************

      SUBROUTINE STDEL (LSTR,POS,LEN)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR,POS,LEN

      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL
      
      INTEGER CLEN
      
      IF (LSTR.EQ.0) THEN
        LSTR = STNEWL(0)
        RETURN
      END IF
      
      IF ((LEN.LE.0).OR.((POS-1).GE.STLEN(LSTR))) RETURN

C Copy the characters from the end of the string to fill the gap.
C Cancel if there's nothing to do.
C CLEN receives the number of characters at the end of CSTR that
C have to be moved.

      CLEN = MAX(0,STLEN(LSTR)-(POS+LEN-1))
      IF (CLEN.GE.STLEN(LSTR)) RETURN
      
      CALL ST_REP (MAX(POS,1),STLEN(LSTR)-CLEN+1,CLEN,
     *             KWORK(L(LSTR)+1),KWORK(L(LSTR)+1))
     
C Set the new length of the string.
     
      CALL STSLN (LSTR,MAX(POS,1)-1+CLEN)

      END
      
***********************************************************************
* Insert one string into another string
*
* Inserts the string CSTR2 into CSTR1 at position POS.
*
* In:
*  LSTR1  - Handle to the source/destination string
*  LSTR2  - Handle to the substring that should be inserted into CSTR1
*  POS    - Position where to insert the substring
*           <=0: insert CSTR2 in the very beginning
*           >length(CSTR2): concatenate CSTR1 and CSTR2
*
* Out:
*  LSTR1  - Handle to the modified string
***********************************************************************
      
      SUBROUTINE STINS (LSTR1,LSTR2,POS)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR1,LSTR2,POS

      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL
      
      INTEGER LTMP,CLEN,CPOS
      
      INTEGER STDUP
      EXTERNAL STDUP
      
      IF (LSTR2.EQ.0) RETURN
      IF (LSTR1.EQ.0) THEN
        LSTR1 = STDUP(LSTR2)
        RETURN
      END IF
      
C What's the real position where to insert?
      
      CPOS = MIN(MAX (1,POS), STLEN(LSTR1)+1)
      
C How many leading characters of the first string have to be copied?
      
      CLEN = MIN(CPOS-1,STLEN(LSTR1))
      
C Create an empty string to hold everything
      
      LTMP = STNEWL (STLEN(LSTR1)+STLEN(LSTR2))
      
C Copy the first CLEN strings of CSTR1
      
      CALL ST_REP (1,1,CLEN,KWORK(L(LTMP)+1),KWORK(L(LSTR1)+1))
      
C Copy CSTR2
      
      CALL ST_REP (CLEN+1,1,STLEN(LSTR2),
     *             KWORK(L(LTMP)+1),KWORK(L(LSTR2)+1))
     
C and add the rest of CSTR1

      CALL ST_REP (CLEN+1+STLEN(LSTR2),CLEN+1,STLEN(LSTR1)-CLEN,
     *             KWORK(L(LTMP)+1),KWORK(L(LSTR1)+1))
     
      CALL STDIS(LSTR1)
      LSTR1 = LTMP
      
      END
      
***********************************************************************
* Insert a character 
*
* Inserts the character c into the string CSTR at position POS.
*
* In:
*  LSTR  - Handle to the source/destination string
*  C     - The character to insert
*  POS   - Position where to insert the character
*          <=0: insert CSTR2 in the very beginning
*          >length(CSTR2): concatenate CSTR1 and CSTR2
*
* Out:
*  LSTR1  - Handle to the modified string
***********************************************************************
      
      SUBROUTINE STINCH (LSTR,C,POS)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR,POS
      CHARACTER C

      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL
      
      INTEGER CPOS,LTMP,CLEN
      
      IF (LSTR.EQ.0) THEN
        LSTR = STNEWL(1)
        CALL STSCH (LSTR,1,C)
        RETURN
      END IF
      
C What's the real position where to insert?
      
      CPOS = MIN(MAX (1,POS), STLEN(LSTR)+1)
      
C How many leading characters of the first string have to be copied?
      
      CLEN = MIN(CPOS-1,STLEN(LSTR))
      
C Create an empty string to hold everything
      
      LTMP = STNEWL (STLEN(LSTR)+1)
      
C Copy the first CLEN strings of CSTR1
      
      CALL ST_REP (1,1,CLEN,KWORK(L(LTMP)+1),KWORK(L(LSTR)+1))
      
C insert the character

      CALL ST_SCH (CLEN+1,C,KWORK(L(LTMP)+1))
      
C and add the rest of CSTR1

      CALL ST_REP (CLEN+2,CLEN+1,STLEN(LSTR)-CLEN,
     *             KWORK(L(LTMP)+1),KWORK(L(LSTR)+1))
     
      CALL STDIS(LSTR)
      LSTR = LTMP
      
      END
      
***********************************************************************
* Extract the first LEN characters
*
* Creates a substring of CSTR containing the first LEN characters
* of CSTR.
*
* In:
*  LSTR  - Handle to the source string
*  LEN   - Length of the substring to be extracted
*
* Out:
*  Return value = Handle to a new string containing the first LEN
*                 characters of CSTR.
***********************************************************************
      
      INTEGER FUNCTION STLEF (LSTR,LEN)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR,LEN

      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL
      
      INTEGER LTMP,CLEN
      
C Create a string large enough to hold everything
      
      CLEN = MIN(LEN,STLEN(LSTR))
      LTMP = STNEWL (CLEN)
      STLEF = LTMP
      
C If LSTR was empty, that's all
      
      IF ((LEN.LE.0).OR.(LSTR.EQ.0)) RETURN
      
C Copy the content

      CALL ST_REP (1,1,CLEN,KWORK(L(LTMP)+1),KWORK(L(LSTR)+1))
      
      END
      
***********************************************************************
* Extract the last LEN characters
*
* Creates a substring of CSTR containing the last LEN characters
* of CSTR.
*
* In:
*  LSTR  - Handle to the source string
*  LEN   - Length of the substring to be extracted
*
* Out:
*  Return value = Handle to a new string containing the last LEN
*                 characters of CSTR.
***********************************************************************
      
      INTEGER FUNCTION STRIG (LSTR,LEN)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR,LEN

      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL
      
      INTEGER LTMP,CLEN
      
C Create a string large enough to hold everything
      
      CLEN = MIN(LEN,STLEN(LSTR))
      LTMP = STNEWL (CLEN)
      STRIG = LTMP
      
C If LSTR was empty, that's all
      
      IF ((LEN.LE.0).OR.(LSTR.EQ.0)) RETURN
      
C Copy the content
      
      CALL ST_REP (1,STLEN(LSTR)-CLEN+1,CLEN,
     *             KWORK(L(LTMP)+1),KWORK(L(LSTR)+1))
      
      END
      
***********************************************************************
* Extract a substring
*
* Extracts a substring from a given string:
*        string = CSTR(POS:POS+LEN)
*
* In:
*  LSTR  - Handle to the source string
*  POS   - Position of the substring to be extracted
*  LEN   - Length of the substring to be extracted
*
* Out:
*  Return value = Handle to a new string containing the substring.
***********************************************************************

      INTEGER FUNCTION STMID (LSTR,POS,LEN)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR,POS,LEN

      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL
      
      INTEGER LTMP,CLEN
      
C Create a string large enough to hold everything

      CLEN = MIN(LEN+POS-1,MIN(LEN,STLEN(LSTR)-POS+1))
      LTMP = STNEWL (CLEN)
      STMID = LTMP
      
C If LSTR was empty, that's all
      
      IF ((CLEN.LE.0).OR.(LSTR.EQ.0)) RETURN
      
C Copy the content
      
      CALL ST_REP (1,MAX(1,POS),CLEN,
     *             KWORK(L(LTMP)+1),KWORK(L(LSTR)+1))
      
      END
      
***********************************************************************
* Determine string trim-positions
*
* Searches for the first and last character in a string that is not
* a space character.
*
* In:
*  CSTR   - The source string
* 
* Out:
*  POS    - Position of the first character
*  LEN    - Length of the partial string with last character <> ' '
***********************************************************************

      SUBROUTINE STTPS (LSTR,POS,LEN)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR,POS,LEN

      INTEGER STLEN
      EXTERNAL STLEN
      
      IF (LSTR.EQ.0) THEN
        POS = 1
        LEN = 0
      END IF
      
      CALL ST_TPS (STLEN(LSTR),POS,LEN,KWORK(L(LSTR)+1))

      END

***********************************************************************
* Trim a string
*
* Deletes leading and trailing space-characters:
*        CSTR = TRIM(CSTR)
*
* In:
*  LSTR   - Handle to the source/destination string
* 
* Out:
*  LSTR   - Handle to the modified string
***********************************************************************

      SUBROUTINE STTRM (LSTR)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR

      INTEGER STLEN,STMID
      EXTERNAL STLEN,STMID
      
      INTEGER POS,LEN,LTMP
      
      IF (LSTR.EQ.0) RETURN
      
      CALL ST_TPS (STLEN(LSTR),POS,LEN,KWORK(L(LSTR)+1))
      LTMP = STMID (LSTR,POS,LEN)
      CALL STCPY(LSTR,LTMP)
      CALL STDIS(LTMP)

      END

***********************************************************************
* Trim a string partially
*
* Deletes leading and/or trailing space-characters, depending on a
* parameter
*        CSTR = TRIMPART(CSTR,ITRIM)
*
* In:
*  LSTR   - Handle to the source/destination string
*  ITRIM  - =0: trim nothing
*           =1: trim leading spaces
*           =2: trim trailing spaces
*           =3: trim both, leading and trailing spaces
* 
* Out:
*  LSTR   - Handle to the modified string
***********************************************************************

      SUBROUTINE STTRP (LSTR,ITRIM)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR,ITRIM

      INTEGER STLEN,STMID
      EXTERNAL STLEN,STMID
      
      INTEGER POS,LEN,LTMP
      
      IF (LSTR.EQ.0) RETURN
      
      CALL ST_TPS (STLEN(LSTR),POS,LEN,KWORK(L(LSTR)+1))
      
      IF (ITRIM.EQ.1) THEN
C       Only leading spaces
        LEN = STLEN(LSTR)-POS+1
      ELSE IF (ITRIM.EQ.2) THEN
C       Only trailing spaces
        LEN = LEN+POS-1
        POS = 1
      END IF
      
      LTMP = STMID (LSTR,POS,LEN)
      CALL STCPY(LSTR,LTMP)
      CALL STDIS(LTMP)

      END

***********************************************************************
* Overwrite a string
*
* Overwrites a string with another string starting at position POS:
*        CSTR1(POS:POS+LEN) = CSTR2
*
* If the destination string is not long enough, it will be extended.
*
* In:
*  LSTR1   - Handle to the source/destination string
*            if =0, a new string will be allocated.
*  LSTR2   - Handle to the source string with that CSTR1 should be
*            overwritten
*  POS1    - First position in CSTR1 that will be overwritten
*            if =0, CSTR1 is overwritten from the beginning
*  POS2    - First character in CSTR2 that is copied to CSTR1.
*            if =0, the first character of CSTR2 will be taken
*  LEN     - Number of characters that will be overwritten.
*            if =0, copy the whole CSTR2 into CSTR1.
* 
* Out:
*  LSTR1   - Handle to the modified string
***********************************************************************

      SUBROUTINE STOWR (LSTR1,LSTR2,POS1,POS2,LEN)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR1,LSTR2
      INTEGER POS1,POS2,LEN

      INTEGER STLEN,STDUP,STNEWL
      EXTERNAL STLEN,STDUP,STNEWL
      
      INTEGER CPOS1,CPOS2,CLEN

      IF (LSTR2.EQ.0) RETURN
      IF (LSTR1.EQ.0) THEN
        LSTR1 = STNEWL(0)
        RETURN
      END IF
      
C Determine position and length of what to copy from/to
      
      CPOS1 = MIN(MAX(1,POS1),STLEN(LSTR1)+1)
      CPOS2 = MIN(MAX(1,POS2),STLEN(LSTR2)+1)
      CLEN = MAX(0,MIN(LEN,STLEN(LSTR2)-CPOS2+1))

C Augment CSTR1 if necessary

      IF ((CPOS1+LEN).GT.STLEN(LSTR1)) THEN
        CALL STSLN (LSTR1,CPOS1+LEN)
      END IF

C Copy the content

      CALL ST_REP (CPOS1,CPOS2,LEN,KWORK(L(LSTR1)+1),KWORK(L(LSTR2)+1))

      END
      
***********************************************************************
* Set string length
*
* Modifies the length of a string. If necessary, the memory will be
* reallocated.
*
* In:
*  LSTR   - Handle to the source/destination string
*  LEN    - Desired length for the string
* 
* Out:
*  LSTR   - Handle to the modified string
***********************************************************************

      SUBROUTINE STSLN (LSTR,LEN)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INTEGER LEN
      
      INTEGER STLEN,STNEWL
      EXTERNAL STLEN,STNEWL

      INTEGER LSTR
      
      INTEGER LTMP,CLEN
      
C How long do we have to make the string?
      
      CLEN = MAX(0,LEN)
      
C Theoretically we could use ZDISP for extending/shrinking
C the memory. But if we do that often with strings that might be
C on the beginning of the heap, this will lead to very time consuming
C string manipulations as the whole heap must be shifted to make space/
C fill gaps. So it's better to reallocate the string at the end of the
C heap and copy the content, so subsequent calls will not affect
C the whole heap anymore!

      LTMP = STNEWL (CLEN)
      CALL ST_REP (1,1,CLEN,KWORK(L(LTMP)+1),KWORK(L(LSTR)+1))
      
      CALL STDIS (LSTR)
      LSTR = LTMP
      
      END
      
***********************************************************************
* Set a character
*
* Replaced a character in the string:
*        CSTR(POS) = C
*
* In:
*  LSTR   - Handle to the destination string
*  POS    - Character position to be replaced
*  C      - The character
***********************************************************************
 
      SUBROUTINE STSCH (LSTR,POS,C)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR
      INTEGER POS
      CHARACTER C
      
      INTEGER STDUP,STLEN      
      EXTERNAL STDUP,STLEN
      
      IF (LSTR.EQ.0) RETURN
      IF ((POS.LE.0).OR.(POS.GT.(STLEN(LSTR)+1))) RETURN
      
      CALL ST_SCH (POS,C,KWORK(L(LSTR)+1))
      
      END

***********************************************************************
* Get a character
*
* Returns the character at the desired position:
*
* In:
*  LSTR   - Handle to the source string
*  POS    - Character position 
*
* Out:
*  Return value = CSTR(POS)
*
* If there was an error (no string, POS out of bounds), this routine
* returns the \0-character!
***********************************************************************

      CHARACTER FUNCTION STGCH (LSTR,POS)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR
      INTEGER POS

      CHARACTER ST_GCH
      EXTERNAL ST_GCH
      
      INTEGER STDUP,STLEN
      EXTERNAL STDUP,STLEN
      
      STGCH=CHAR(0)
      
      IF (LSTR.EQ.0) RETURN
      IF ((POS.LE.0).OR.(POS.GT.(STLEN(LSTR)+1))) RETURN

      STGCH = ST_GCH(POS,KWORK(L(LSTR)+1))
      
      END

***********************************************************************
* Case-sensitive string-comparison
*
* Compares two strings, case-sensitive
*
* In:
*  LSTR1  - Handle to the first string
*  LSTR2  - Handle to the second string
*
* Out:
*  Return value = -1, if CSTR1 < CSTR2
*                  0, if CSTR1 = CSTR2
*                  1, if CSTR1 > CSTR2
***********************************************************************

      INTEGER FUNCTION STCMP(LSTR1,LSTR2)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR1,LSTR2

      INTEGER STLEN,ST_CMP
      EXTERNAL STLEN,ST_CMP
      
      IF ((LSTR1.EQ.0).AND.(LSTR1.EQ.0)) THEN
        STCMP = 0
        RETURN
      ELSE IF ((LSTR1.EQ.0).AND.(LSTR1.NE.0)) THEN
        STCMP = -1
        RETURN
      ELSE IF ((LSTR1.NE.0).AND.(LSTR1.EQ.0)) THEN
        STCMP = 1
        RETURN
      END IF
      
      STCMP = ST_CMP(STLEN(LSTR1),STLEN(LSTR2),
     *               KWORK(L(LSTR1)+1),KWORK(L(LSTR1)+1))
      
      END

***********************************************************************
* Case-insensitive string-comparison
*
* Compares two strings, case-insensitive
*
* In:
*  LSTR1  - Handle to the first string
*  LSTR2  - Handle to the second string
*
* Out:
*  Return value = -1, if CSTR1 < CSTR2
*                  0, if CSTR1 = CSTR2
*                  1, if CSTR1 > CSTR2
***********************************************************************

      INTEGER FUNCTION STICM(LSTR1,LSTR2)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR1,LSTR2

      INTEGER STLEN,ST_ICM
      EXTERNAL STLEN,ST_ICM
      
      IF ((LSTR1.EQ.0).AND.(LSTR1.EQ.0)) THEN
        STICM = 0
        RETURN
      ELSE IF ((LSTR1.EQ.0).AND.(LSTR1.NE.0)) THEN
        STICM = -1
        RETURN
      ELSE IF ((LSTR1.NE.0).AND.(LSTR1.EQ.0)) THEN
        STICM = 1
        RETURN
      END IF
      
      STICM = ST_ICM(STLEN(LSTR1),STLEN(LSTR2),
     *               KWORK(L(LSTR1)+1),KWORK(L(LSTR1)+1))
      
      END

***********************************************************************
* Partial string-comparison, case-sensitive
*
* Compares two strings, case-sensitive
* Compares at most the first LEN characters
*
* In:
*  LSTR1  - Handle to the first string
*  LSTR2  - Handle to the second string
*  LEN    - Maximum length that should be compared.
*           If at least one string is shorter, the minimum string
*           length determines how many characters are compared.
*
* Out:
*  Return value = -1, if CSTR1(1:LEN) < CSTR2(1:LEN)
*                  0, if CSTR1(1:LEN) = CSTR2(1:LEN)
*                  1, if CSTR1(1:LEN) > CSTR2(1:LEN)
***********************************************************************

      INTEGER FUNCTION STNCM(LSTR1,LSTR2,LEN)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR1,LSTR2,LEN

      INTEGER STLEN,ST_CMP
      EXTERNAL STLEN,ST_CMP
      
      INTEGER CLEN
      
      IF ((LSTR1.EQ.0).AND.(LSTR1.EQ.0)) THEN
        STNCM = 0
        RETURN
      ELSE IF ((LSTR1.EQ.0).AND.(LSTR1.NE.0)) THEN
        STNCM = -1
        RETURN
      ELSE IF ((LSTR1.NE.0).AND.(LSTR1.EQ.0)) THEN
        STNCM = 1
        RETURN
      END IF
      
C Determine the maximum length
      
      CLEN = MAX(0,MIN(MIN(LEN,STLEN(LSTR1)),STLEN(LSTR2)))
      
C Compare the first CLEN characters only
      
      STNCM = ST_CMP(CLEN,CLEN,
     *               KWORK(L(LSTR1)+1),KWORK(L(LSTR1)+1))
      
      END

***********************************************************************
* Search for first occurance of a character
*
* Searches for a character in a string, from left to right
*
* In:
*  LSTR  - Handle to the string where to search
*  C     - The character to search for
*
* Out:
*  Return value = First position in CSTR containing C
*                 0, if CSTR doesn't contain C
***********************************************************************

      INTEGER FUNCTION STCHR (LSTR,C)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR
      CHARACTER C
      
      INTEGER ST_CHR,STLEN
      EXTERNAL ST_CHR,STLEN
      
      IF (LSTR.EQ.0) THEN
        STCHR = 0
        RETURN
      END IF
      
      STCHR = ST_CHR (STLEN(LSTR),C,KWORK(L(LSTR)+1))
      
      END

***********************************************************************
* Search for last occurance of a character
*
* Searches for a character in a string, from right to left
*
* In:
*  LSTR  - Handle to the string where to search
*  C     - The character to search for
*
* Out:
*  Return value = Last position in CSTR containing C
*                 0, if CSTR doesn't contain C
***********************************************************************

      INTEGER FUNCTION STRCH (LSTR,C)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR
      CHARACTER C
      
      INTEGER ST_RCH,STLEN
      EXTERNAL ST_RCH,STLEN

      IF (LSTR.EQ.0) THEN
        STRCH = 0
        RETURN
      END IF
      
      STRCH = ST_RCH (STLEN(LSTR),C,KWORK(L(LSTR)+1))
      
      END

***********************************************************************
* Search for first occurance of a substring in a string
*
* Tests if CSTR1 contains the substring CSTR2.
*
* In:
*  LSTR1  - Handle to the string where to search
*  LSTR2  - The substring to search for
*  POS    - Position where to start the search
*
* Out:
*  Return value = First position in CSTR where CSTR2 starts
*                 0, if CSTR1 doesn't contain CSTR2
***********************************************************************

      INTEGER FUNCTION STSTR (LSTR1,LSTR2,POS)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR1,LSTR2,POS
      
      INTEGER ST_STR,STLEN
      EXTERNAL ST_STR,STLEN
      
      IF (LSTR1.EQ.0) THEN
        STSTR = 0
        RETURN
      END IF
      IF (LSTR2.EQ.0) THEN
        STSTR = 1
        RETURN
      END IF
      
      STSTR = ST_STR (STLEN(LSTR1),STLEN(LSTR2),POS,
     *                KWORK(L(LSTR1)+1),KWORK(L(LSTR2)+1))
      
      END

***********************************************************************
* Left-justify string
*
* Left-justify the string CSTR without changing its length.
*
* In:
*  LSTR   - Handle to the string
*
* Out:
*  LSTR   - The justified string
***********************************************************************

      SUBROUTINE STLJS (LSTR)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR
      
      INTEGER STLEN
      EXTERNAL STLEN
      
      INTEGER POS,LEN

C Obtain the position/length of the text

      CALL STTPS (LSTR,POS,LEN)
      
C Return if there's nothing to do.
      
      IF ((POS.LE.1).OR.(LEN.EQ.0)) RETURN
      
C Copy the content. Recognize that in this situation LSTR is <> 0
C because of the result of STTPS and the above IF-clause!

      CALL ST_REP (1,POS,LEN,KWORK(L(LSTR)+1),KWORK(L(LSTR)+1))
      
C Fill the rest with spaces. We only need to fill that amount that we
C have moved, the rest is (because of STTPS) already filled with 
C spaces!

      CALL ST_SET (1+LEN,POS-1,' ',KWORK(L(LSTR)+1))

      END
      
***********************************************************************
* Right-justify string
*
* Right-justify the string CSTR without changing its length.
*
* In:
*  LSTR   - Handle to the string
*
* Out:
*  LSTR   - The justified string
***********************************************************************

      SUBROUTINE STRJS (LSTR)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR
      
      INTEGER STLEN,STMID
      EXTERNAL STLEN,STMID
      
      INTEGER POS,LEN,LTMP

C Obtain the position/length of the text

      CALL STTPS (LSTR,POS,LEN)
      
C Return if there's nothing to do.
      
      IF ((POS.LE.0).OR.(LEN.EQ.0).OR.((POS+LEN-1).GE.STLEN(LSTR))) 
     *  RETURN
      
C Create a new string with only that content.
C This is a little bit more work that in STLJS because the later 
C copy process copies "from left to right" - and we don't want to 
C overwrite what we have to copy :)

      LTMP = STMID (LSTR,POS,LEN)
      
C Fill the string with spaces - we only have to fill what contains
C non-space characters...

      CALL ST_SET (POS,LEN,' ',KWORK(L(LSTR)+1))

C Copy the content back to the correct position
      
      CALL ST_REP (1+STLEN(LSTR)-LEN,POS,LEN,
     *             KWORK(L(LSTR)+1),KWORK(L(LSTR)+1))

      CALL STDIS (LTMP)      

      END

***********************************************************************
* Center-justify string
*
* Center-justify the string CSTR without changing its length.
*
* In:
*  LSTR   - Handle to the string
*
* Out:
*  LSTR   - The justified string
***********************************************************************

      SUBROUTINE STCJS (LSTR)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR
      
      INTEGER STLEN,STMID
      EXTERNAL STLEN,STMID
      
      INTEGER POS,LEN,LTMP,M

C Obtain the position/length of the text

      CALL STTPS (LSTR,POS,LEN)
      
C Calculate the position where the string should be copied to

      M = MAX(1,STLEN(LSTR)/2 - LEN/2)
      
C Return if there's nothing to do.
      
      IF ((POS.LE.0).OR.(LEN.EQ.0).OR.(POS.EQ.M)) RETURN
      
C Create a new string with only that content.

      LTMP = STMID (LSTR,POS,LEN)
      
C Fill the string with spaces - we only have to fill what contains
C non-space characters...

      CALL ST_SET (POS,LEN,' ',KWORK(L(LSTR)+1))

C Copy the content back
      
      CALL ST_REP (M,POS,LEN,KWORK(L(LSTR)+1),KWORK(L(LSTR)+1))

      CALL STDIS (LTMP)      

      END

***********************************************************************
* String-parsing
*
* Searches for the next word in a string that is enclosed by spaces.
*
* SPOS tracks the last found space character.
* For the search for the first word, SPOS has to be set to 0.
* The routine always updates this variable to find the next word.
* If the end of the string is reache, i.e. no word was found,
* SPOS is set to -1.
*
* In:
*  LSTR   - The string to parse
*  SPOS   - The position of the last found space character or
*           -1 if the last search failed.
*           Must be set to 0 for first search.
*  C      - separation character for words; should be set to " ".
*  BENC   - Whether to take care of "highlighting"-characters that
*           enclose strings. If set to .TRUE., separation characters
*           that happen in string sequences enclosed by '"'-signs
*           are ignored. This way it's possible to parse strings like
*               abc "def ghi" jkl
*           Returns:
*            1.) abc
*            2.) "def ghi"
*            3.) jkl
*           Can be set to .FALSE. for standard word-parsing.
*
* Out:
*  SPOS   - Updated space character position
*  POS    - Position of the first character <> C of the next found
*           word or 0 if nothing was found
*  LEN    - Length of the word that was found.
***********************************************************************

      SUBROUTINE STPRS (LSTR,C,BENC,SPOS,POS,LEN)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INTEGER LSTR,SPOS,POS,LEN
      CHARACTER C
      LOGICAL BENC
      
      INTEGER STLEN,STMID
      EXTERNAL STLEN,STMID

C Determine if we have anything to search at all

      POS = 0
      LEN = 0

      IF (SPOS.LT.0) RETURN
      
      IF ((SPOS-1).GE.STLEN(LSTR)) THEN
C We have already reached the end
        SPOS = -1
        RETURN
      END IF

C Search for the next word

      CALL ST_PRS (STLEN(LSTR),C,BENC,SPOS,POS,LEN,KWORK(L(LSTR)+1))
      
C Set SPOS to the character after the found word - if anything 
C was found

      IF (POS.EQ.0) THEN
        SPOS = -1
      ELSE
        SPOS = POS+LEN
      END IF

      END
      
***********************************************************************
* Convert to uppercase
*
* Converts a string to its uppercase representation
*
* In:
*  LSTR   - Handle to the string
*
* Out:
*  LSTR   - Handle to the string in uppercase representation
***********************************************************************

      SUBROUTINE STUPR (LSTR)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR
      
      INTEGER STLEN
      EXTERNAL STLEN
      
      CALL ST_UPR (STLEN(LSTR),KWORK(L(LSTR)+1))
      
      END
      
***********************************************************************
* Convert to lowercase
*
* Converts a string to its lowercase representation
*
* In:
*  LSTR   - Handle to the string
*
* Out:
*  LSTR   - Handle to the string in lowercase representation
***********************************************************************

      SUBROUTINE STLWR (LSTR)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR

      INTEGER STLEN
      EXTERNAL STLEN
      
      CALL ST_LWR (STLEN(LSTR),KWORK(L(LSTR)+1))
      
      END
      

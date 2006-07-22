***********************************************************************
* This file offeres a routine for conditional output to terminal /
* into file. It can be used for a quick-and-notsodirty direction
* of string output.
***********************************************************************

***********************************************************************
* Conditional output of standard strings
*
* This routine writes a string CSTR of arbitrary length to the
* output channels MTERM and MFILE as long as these are <> 0.
* the variable MT decides on whether the output is made or not.
*
* In:
*   MT     - Output flag.
*            <= 0: No output
*             = 1: Write CSTR to output channel MFILE
*            >= 2: Write CSTR to putput channels MFILE and MTERM
*   MTERM  - Output channel of terminal output.
*            =0: Prevent output to terminal
*   MFILE  - Output channel of file output.
*            =0: Prevent output to file
*            The file must be open for proper output
*   BNL    - Whether to break the line / write a newline-character after
*            writing the string.
*   CSTR   - string to be written to MTERM / MFILE
*
* The routine automatically recognizes the length of the string
* CSTR, thus removing and trailing space characters (which makes
* FORTRAN output ugly otherwise...)
***********************************************************************

      SUBROUTINE CNOUTS (MT,MTERM,MFILE,BNL,CSTR)
      
      IMPLICIT NONE
      
      INCLUDE 'dstrings.inc'
      
      INTEGER MT,MTERM,MFILE
      CHARACTER CSTR*(*)
      LOGICAL BNL
      
      INTEGER LSTR
      
C     MT < 0 prevents output

      IF (MT.LT.0) RETURN
      
C     Create a truncated DString from the Fortran string.
C     Only trim trailing spaces.

      LSTR = STNEWC (.FALSE.,CSTR)
      
      IF (LSTR.EQ.0) THEN
        
C       Fatal error, no handles. Must be printed anywhere!!!        
        
        IF (MTERM.NE.0) THEN
          WRITE (MTERM,'(A)') 
     *      'CNOUTS Error: Could not generate message! No free handles!'
          WRITE (MTERM,'(A,A)') 'Last message: ',CSTR
        END IF
        IF (MTERM.NE.0) THEN
          WRITE (MFILE,'(A)') 
     *      'CNOUTS Error: Could not generate message! No free handles!'
          WRITE (MTERM,'(A,A)') 'Last message: ',CSTR
        END IF
        IF ((MTERM.EQ.0).AND.(MFILE.EQ.0)) THEN
          WRITE (*,'(A)') 
     *      'CNOUTS Error: Could not generate message! No free handles!'
          WRITE (MTERM,'(A,A)') 'Last message: ',CSTR
        END IF
        
        RETURN
        
      END IF
      
      CALL STTRP (LSTR,2)
      
C     Write to terminal?

      IF ((MT.GE.2).AND.(MTERM.GT.0))
     *  CALL STOUT (LSTR, MTERM, BNL)
     
C     Write to file?
      
      IF ((MT.GE.1).AND.(MFILE.GT.0))
     *  CALL STOUT (LSTR, MFILE, BNL)

C     Delete the string

      CALL STDIS (LSTR)

      END
      
***********************************************************************
* Conditional output of a DString
*
* This routine writes a dstring LSTR of arbitrary length to the
* output channels MTERM and MFILE as long as these are <> 0.
* the variable MT decides on whether the output is made or not.
*
* In:
*   MT     - Output flag.
*            <= 0: No output
*             = 1: Write CSTR to output channel MFILE
*            >= 2: Write CSTR to putput channels MFILE and MTERM
*   MTERM  - Output channel of terminal output.
*            =0: Prevent output to terminal
*   MFILE  - Output channel of file output.
*            =0: Prevent output to file
*            The file must be open for proper output
*   BNL    - Whether to break the line / write a newline-character after
*            writing the string.
*   LSTR   - DString to be written to MTERM / MFILE
***********************************************************************

      SUBROUTINE CNOUTD (MT,MTERM,MFILE,BNL,LSTR)
      
      IMPLICIT NONE
      
      INCLUDE 'dstrings.inc'
      
      INTEGER MT,MTERM,MFILE,LSTR
      LOGICAL BNL
      
C     MT < 0 prevents output

      IF (MT.LT.0) RETURN
      
C     Write to terminal?

      IF ((MT.GE.2).AND.(MTERM.GT.0))
     *  CALL STOUT (LSTR, MTERM, BNL)
     
C     Write to file?
      
      IF ((MT.GE.1).AND.(MFILE.GT.0))
     *  CALL STOUT (LSTR, MFILE, BNL)

      END
      
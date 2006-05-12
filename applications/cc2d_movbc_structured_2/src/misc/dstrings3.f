***********************************************************************
* Convert a DString into a Fortran string
*
* Writes the given string into the Fortrag string CSTR.
*
* Remark: In this routine the "source" string is in the first
* parameter while the second parameter is the "destination"!
*
* In:
*  LSTR   - Handle to the source DString
*  CSTR   - Fortran string where to write the DString to
*
* Out:
*  CSTR   - Fortran string overwritten with the DString,
*           filled up with trailing spaces.
***********************************************************************
     
      SUBROUTINE STPUT (LSTR, CSTR)
      
      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      
      INTEGER LSTR
      CHARACTER*(*) CSTR
      
      INTEGER STLEN
      EXTERNAL STLEN
      
      INTEGER CLEN
      
C Fill the Fortran string with spaces
      
      CALL ST_SET(1,LEN(CSTR),' ',CSTR)
      
C And overwrite it with out DString
      
      CLEN = MIN(STLEN(LSTR),LEN(CSTR))
      CALL ST_REP(1,1,CLEN,CSTR,KWORK(L(LSTR)+1))
      
      END

      

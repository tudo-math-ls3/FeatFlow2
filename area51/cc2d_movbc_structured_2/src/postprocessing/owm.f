************************************************************************
* This file contains routines for writing out matrices in plain
* text format to a file.
************************************************************************

************************************************************************
* Auxiliary routine: Create output format strings
*
* This routine creates different output format strings for the output
* of numbers to a file.
*
* In: 
*   NS      : Number format
*             =1: Double precision format
*             =2: Real format
*             =3: Integer format
*   FMT     : General format string of the output, e.g. '(D20.10)' for
*             double precision output (with NS=1)
*
* Out:
*   CFORM1  : FMT without line break
*   CFORM2  : FMT with line break
*   CFORM3  : Format for character array with length given by FMT,
*             without line break
*   CFORM4  : Format for character array with length given by FMT,
*             with line break
*   CLEN    : Length of the format string
************************************************************************

      SUBROUTINE OWPARM (NS,CFMT,CFORM1,CFORM2,CFORM3,CFORM4,CLEN)
      
      IMPLICIT NONE
      
      INCLUDE 'dstrings.inc'
      
      INTEGER NS,CLEN
      CHARACTER CFMT*(*),CFORM1*(*),CFORM2*(*),CFORM3*(*),CFORM4*(*)
      
C     local variables

      INTEGER LHND,IPOST
      CHARACTER C*(128)
      
C     Cancel on wrong number format:
      
      IF ((NS.LT.1).OR.(NS.GT.3)) GOTO 99999
      
C     Determine the length of the output string
      
      IF (NS.EQ.1) THEN
        WRITE (C,CFMT) 0D0
      ELSE IF (NS.EQ.2) THEN
        WRITE (C,CFMT) 0E0
      ELSE
        WRITE (C,CFMT) 0
      END IF
      LHND = STNEWC (.FALSE.,C)
      CALL STTRP (LHND,2)
      CLEN = STLEN (LHND)
      CALL STDIS (LHND)
      
C     Create CFORM1/CFORM2 from that - with and without the "$" before
C     the bracket for the line break.

      LHND = STNEWC (.TRUE.,CFMT)
      CALL STPUT (LHND,CFORM2)
      
      IPOST = MAX(STLEN (LHND),1)
      CALL STINSC (LHND, IPOST, .FALSE., '$')
      CALL STPUT (LHND, CFORM1)
      
      CALL STDIS (LHND)
      
C     Create empty output string for the case that zeroes are not to
C     be written

      LHND = STNEWC (.FALSE.,'(A')
      CALL STCATI (LHND,CLEN,0,.FALSE.)
      CALL STCATC (LHND,.FALSE.,'$)')
      CALL STPUT (LHND, CFORM3)
      CALL STDIS (LHND)

      LHND = STNEWC (.FALSE.,'(A')
      CALL STCATI (LHND,CLEN,0,.FALSE.)
      CALL STCATC (LHND,.FALSE.,')')
      CALL STPUT (LHND, CFORM4)
      CALL STDIS (LHND)

99999 END

************************************************************************
* Write full double precision matrix into a text file.
*
* This writes an array DA with NROW rows and NCOL columns as a matrix
* into a text file.
*
* In:
*   DA     - array [1..NROW,1..NCOL] of double
*            The matrix.
*   NROW   - Number of rows
*   NCOL   - Number of columns
*   ARR    - Name of the matrix
*   BCPRSS - Suppress zeroes in output.
*            =true : Zeroes are not written
*            =false: Zeroes are written
*   MFILE  - Output channel to use for the matrix output.
*            = 0: Use default channel 69. Open the file 'CFILE', write
*                 the matrix and close the file afterwards.
*            <>0: Write to channel MFILE. Don't close the cahnnel.
*                 CFILE is ignored.
*   CFILE  - If MFILE=0  : Name of the file where to write to.
*            If MFILE<>0 : Ignored.
*   CFMT   - Format string to use for the output; e.g. '(D20.10)'
************************************************************************

      SUBROUTINE OWM11(DA,NROW,NCOL,BCPRSS,MFILE,ARR,CFILE,CFMT)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'dstrings.inc'
      
C     parameters
      
      INTEGER NROW,NCOL,MFILE
      LOGICAL BCPRSS
      DOUBLE PRECISION DA(NROW,NCOL)
      CHARACTER ARR*(*),CFILE*(*),CFMT*(*)

C     local variables

      INTEGER MFL1,IFMT,IND,I,J
      DOUBLE PRECISION DW
      CHARACTER CFORM1*(64),CFORM2*(64),CFORM3*(64),CFORM4*(64)
      
C     Get the format strings from CFMT
      
      CALL OWPARM (1,CFMT,CFORM1,CFORM2,CFORM3,CFORM4,IND)

      IF (NROW.LE.0) RETURN
      IF (NCOL.LE.0) RETURN

C     Open the file if necessary

      MFL1 = MFILE
      IF (MFL1.EQ.0) THEN
        IFMT = 1
        MFL1 = 69
        CALL OF0(MFL1,CFILE,IFMT)
        IF (IER.NE.0) RETURN
      END IF

C     Write all format strings into the file

      WRITE (MFL1,'(5A15,L1,3I15)') ARR,CFORM1,CFORM2,CFORM3,CFORM4,
     *                              BCPRSS,IND,NROW,NCOL

C     Save the matrix

      DO I=0,NROW-1
        DO J=0,NCOL-2
          DW = DA(I,J)
          IF ((.NOT.BCPRSS).OR.(DW.NE.0D0)) THEN
            WRITE (MFL1,CFORM1) DW
          ELSE
            WRITE (MFL1,CFORM3) '.'
          END IF
        END DO
        DW = DA(I,J)
        IF ((.NOT.BCPRSS).OR.(DW.NE.0D0)) THEN
          WRITE (MFL1,CFORM2) DW
        ELSE
          WRITE (MFL1,CFORM4) '.'
        END IF
      END DO

C     Close the file if necessary

      IF (MFILE.EQ.0) CLOSE(MFL1)

99999 END

************************************************************************
* Expand matrix row
*
* This routine expands the row of a matrix-7/matrix-9 matrix structure
* to a full vector.
*
* In:
*   DROW   - Compressed row in the matrix-7 / matrix-9 matrix
*   KCOL   - Column numbers corresponding to the entries in DROW
*   NCOL   - Number of columns in the row
*   NEQ    - Length of the uncompressed row
*
* Out:
*   DROWO  - array [1..NEQ] of double
*            Uncompressed vector. All matrix entries in DROW which
*            are exactly =0 are written as 1D-99 into DROWO,
*            so the caller can later find out which entries in the
*            original matrix have been 0.0.
************************************************************************

      SUBROUTINE SRRW71 (DROW,KCOL,NCOL,NEQ,DROWO)
      
      IMPLICIT NONE

C     parameters
      
      INTEGER NCOL,KCOL(NCOL),NEQ
      DOUBLE PRECISION DROW(NCOL),DROWO(NEQ)
      
C     local variables

      INTEGER I
      
C     Clear the vector

      CALL LCL1(DROWO,NEQ)
      
C     Build the vector

      DO I=1,NCOL
      
        DROWO(KCOL(I)) = DROW(I)
        
C       Replace 0.0 by 1D-99

        IF (DROW(I).EQ.0D0) DROWO(KCOL(I)) = 1D-99
        
      END DO
      
      END

************************************************************************
* Write matrix-7/matrix-9 double precision matrix into a text file.
*
* This writes an a matrix DA/KCOL/KLD in matrix structure 7 or 9
* into a text file.
*
* In:
*   DA,
*   KCOL,
*   KLD    - The matrix in structure 7 or structure 9.
*   NROW   - Number of rows
*   NCOL   - Number of columns; must be =NROW for structure-7 matrices.
*   ARR    - Name of the matrix
*   BCPRSS - Suppress zeroes in output.
*            =true : Zeroes are not written
*            =false: Zeroes are written
*   MFILE  - Output channel to use for the matrix output.
*            = 0: Use default channel 69. Open the file 'CFILE', write
*                 the matrix and close the file afterwards.
*            <>0: Write to channel MFILE. Don't close the cahnnel.
*                 CFILE is ignored.
*   CFILE  - If MFILE=0  : Name of the file where to write to.
*            If MFILE<>0 : Ignored.
*   CFMT   - Format string to use for the output; e.g. '(D20.10)'
************************************************************************

      SUBROUTINE OWM17(DA,KCOL,KLD,NROW,NCOL,BCPRSS,MFILE,
     *                 ARR,CFILE,CFMT)
      
      IMPLICIT NONE
      
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'
      
      INCLUDE 'dstrings.inc'
      
C     parameters
      
      INTEGER NROW,NCOL,MFILE,KCOL(*),KLD(*)
      LOGICAL BCPRSS
      DOUBLE PRECISION DA(*)
      CHARACTER ARR*(*),CFILE*(*),CFMT*(*)

C     local variables

      INTEGER MFL1,IFMT,LVEC,IND,I,J
      DOUBLE PRECISION DW
      CHARACTER CFORM1*(64),CFORM2*(64),CFORM3*(64),CFORM4*(64)
      
C     Get the format strings from CFMT
      
      CALL OWPARM (1,CFMT,CFORM1,CFORM2,CFORM3,CFORM4,IND)

C     Open the file if necessary

      MFL1 = MFILE
      IF (MFL1.EQ.0) THEN
        IFMT = 1
        MFL1 = 69
        CALL OF0(MFL1,CFILE,IFMT)
        IF (IER.NE.0) RETURN
      END IF

C     Write all format strings into the file

      WRITE (MFL1,'(5A15,L1,3I15)') ARR,CFORM1,CFORM2,CFORM3,CFORM4,
     *                              BCPRSS,IND,NROW,NCOL

      IF (NROW.LE.0) RETURN
      IF (NCOL.LE.0) RETURN

C     Allocate a temporary vector where to build each row

      CALL ZNEW (NCOL,1,LVEC,'TMPARR ')  

C     Save the matrix

      DO I=1,NROW
      
C       Extract row I  
      
        CALL SRRW71 (DA(KLD(I)),KCOL(KLD(I)),KLD(I+1)-KLD(I),NCOL,
     *               DWORK(L(LVEC)))
     
C       Write the row.
C       All entries = 0 are only written out if BCPRESS=false.
C       All entries < 1D-90 are written out as 0.
     
        DO J=1,NCOL-1
          DW = DWORK(L(LVEC)+J-1)
          IF ((.NOT.BCPRSS).OR.(DW.NE.0D0)) THEN
            IF (ABS(DW).LE.1D-90) THEN
              WRITE (MFL1,CFORM3) '0'
            ELSE
              WRITE (MFL1,CFORM1) DW
            END IF  
          ELSE
            WRITE (MFL1,CFORM3) '.'
          END IF
        END DO
        DW = DWORK(L(LVEC)+NCOL-1)
        IF ((.NOT.BCPRSS).OR.(DW.NE.0D0)) THEN
          IF (ABS(DW).LE.1D-90) THEN
            WRITE (MFL1,CFORM4) '0'
          ELSE   
            WRITE (MFL1,CFORM2) DW
          END IF  
        ELSE
          WRITE (MFL1,CFORM4) '.'
        END IF
      END DO
      
C     Release temporary vector
      
      CALL ZDISP (0,LVEC,'LVEC')

C     Close the file if necessary

      IF (MFILE.EQ.0) CLOSE(MFL1)

99999 END

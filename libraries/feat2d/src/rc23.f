************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.0)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XRC23                                                                *
*                                                                      *
* Purpose  Call of RC23                                                *
*          Compress DWORK                                              *
*                                                                      *
* Subroutines/functions called   RC23, ZDISP, ZNEW                     *
*                                                                      *
* Version from  11/19/90                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Numbers of the corresponding arrays                  *
* LDIA     I*4    set by ZNEW                                          *
* LDIAS    I*4                                                         *
* IDISP    I*4    1  means call of ZDISP after deletion of zero entries*
* ARR      C*6    Names of blocks (for error messages only)            *
* For the description of the remaining parameters see RC23             *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 -114 Type of at least one array is not double prec.  *
*                 -115 Length of at least one array is < NA            *
*                                                                      *
************************************************************************
*                                                                      *
* RC23                                                                 *
*                                                                      *
* Purpose  Removing entries of small modulus from an array             *
*          Storage technique 3/4                                       *
*          Single precision version                                    *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  11/19/90                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Field of matrices                                    *
* KDIA     I*4    Pointer vectors                                      *
* KDIAS    I*4                                                         *
* NA       I*4    length of VA                                         *
* NDIA     I*4    number of diagonal rows                              *
* NEQ      I*4    number of rows                                       *
* NBLOC    I*4    Number of matrices stored on VA                      *
* KOFF     I*4    Matrix IBLOC starts at position KOFF(IBLOC)+1 on VA  *
* TOL      R*8    entries of  modulus <= TOL  are deleted              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VA       R*4    Compressed matrices                                  *
* KDIA     I*4    described in the same storage technique              *
* KDIAS    I*4                                                         *
* NA       I*4                                                         *
*                                                                      *
************************************************************************
C
      SUBROUTINE XRC23(LA,LDIA,LDIAS,NDIA,NA,NEQ,NBLOC,TOL,IDISP,ARR)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      DIMENSION ARR(*),LA(*)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='XRC23'
      IF (ICHECK.GE.997) CALL OTRC('XRC23 ','11/19/90')
      IER=0
C
      DO 1 IBLOC=1,NBLOC
C *** Check input parameter
      CALL ZTYPE(LA(IBLOC),ITYPE)
      IF (ITYPE.NE.2) THEN
       WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
       CALL WERR(-114,'XRC23 ')
       GOTO 99999
      ENDIF
      CALL ZLEN(LA(IBLOC),ILEN)
      IF (ILEN.LT.NA) THEN
       WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
       CALL WERR(-115,'XRC23 ')
       GOTO 99999
      ENDIF
1     CONTINUE
C
      CALL ZNEW(NBLOC,-3,LOFF,'KOFF  ')
      IF (IER.NE.0) GOTO 99999
C
      DO 2 IBLOC=1,NBLOC
      KWORK(L(LOFF)+IBLOC-1)=L(LA(IBLOC))-1
2     CONTINUE
C
      CALL RC23(VWORK(1),KWORK(L(LDIA)),KWORK(L(LDIAS)),NDIA,
     *          NA,NEQ,NBLOC,KWORK(L(LOFF)),TOL)
C
      CALL ZDISP(0,LOFF,'KOFF  ')
      IF (IDISP.EQ.1) THEN
       DO 3 IBLOC=1,NBLOC
3      CALL ZDISP(NA,LA(IBLOC),ARR(IBLOC))
       CALL ZDISP(NDIA  ,LDIA ,'KDIA  ')
       CALL ZDISP(NDIA+1,LDIAS,'KDIAS ')
      ENDIF
C
99999 END
C
C
C
      SUBROUTINE RC23(VA,KDIA,KDIAS,NDIA,NA,NEQ,NBLOC,KOFF,TOL)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KDIA(*),KDIAS(*),KOFF(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      SUB='RC23'
      IF (ICHECK.GE.997) CALL OTRC('RC23  ','11/19/90')
C
      DO 1 IBLOC=1,NBLOC
      IOFF=KOFF(IBLOC)
      DO 2 IA=1,NA
      IF (ABS(VA(IOFF+IA)).LT.TOL) VA(IOFF+IA)=0.
2     CONTINUE
1     CONTINUE
C
      IDIA0 =2
      IDIA1 =2
      IDIAS1=NEQ+1
C
      DO 10 IDIA=2,NDIA
C
      IF (IDIA0.GT.NDIA) GOTO 50
C
      ILEN=NEQ-ABS(KDIA(IDIA0))
      DO 20 IBLOC=1,NBLOC
      IOFF=KOFF(IBLOC)+KDIAS(IDIA0)
      DO 21 J=0,ILEN-1
      IF (ABS(VA(IOFF+J)).GE.TOL) GOTO 30
21    CONTINUE
20    CONTINUE
      IF (IDIA0.EQ.NDIA) GOTO 50
C
      IDIA0=IDIA0+1
      ILEN=NEQ-ABS(KDIA(IDIA0))
C
30    DO 40 IBLOC=1,NBLOC
      IOFF0=KOFF(IBLOC)+KDIAS(IDIA0)
      IOFF1=KOFF(IBLOC)+IDIAS1
      DO 41 J=0,ILEN-1
      VA(IOFF1+J)=VA(IOFF0+J)
41    CONTINUE
40    CONTINUE
      KDIA (IDIA1)=KDIA(IDIA0)
      KDIAS(IDIA1)=IDIAS1
C
      IDIA0 =IDIA0 +1
      IDIA1 =IDIA1 +1
      IDIAS1=IDIAS1+ILEN
C
10    CONTINUE
C
50    NDIA=IDIA1-1
      KDIAS(IDIA1)=IDIAS1
      NA=IDIAS1-1
C
      END

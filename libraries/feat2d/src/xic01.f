************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XIC01n                                                               *
*                                                                      *
* Purpose  Call IC01n                                                  *
*                                                                      *
* Subroutines/functions called   ZTYPE, IC01n                          *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Numbers of the arrays describing the matrix in       *
* LCOL     I*4    Storage technique n                                  *
* LLD      I*4                                                         *
* LOP      I*4                                                         *
* LX,LB    I*4    Numbers of the solution and the right hand side      *
* OMEGA    R*8    Relaxation parameter                                 *
* For the description of the remaining input parameters see IC01n      *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 -170  X or B are not double precision                *
*                                                                      *
************************************************************************
C
      SUBROUTINE XIC017(LA,LCOL,LLD,LX,LB,NEQ,NIT,ITE,EPS,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /ERRCTL/,/CHAR/
C
      SUB='XIC017'
      IF (ICHECK.GE.997) CALL OTRC('XIC017','01/02/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.1.OR.ITYPE2.NE.1.OR.ITYPE3.NE.1) THEN
       CALL WERR(-170,'XIC017')
       GOTO 99999
      ENDIF
C
      CALL IC017(DWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),
     *           DWORK(L(LX)),DWORK(L(LB)),NEQ,NIT,ITE,EPS,OMEGA)
C
99999 END
C
C
C
      SUBROUTINE XIC01A(LA,LCOL,LLD,LOP,LX,LB,NEQ,NIT,ITE,EPS,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /ERRCTL/,/CHAR/
C
      SUB='XIC01A'
      IF (ICHECK.GE.997) CALL OTRC('XIC01A','01/02/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.1.OR.ITYPE2.NE.1.OR.ITYPE3.NE.1) THEN
       CALL WERR(-170,'XIC01A')
       GOTO 99999
      ENDIF
C
      CALL IC01A(DWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),
     *           KWORK(L(LOP)),DWORK(L(LX)),DWORK(L(LB)),NEQ,
     *           NIT,ITE,EPS,OMEGA)
C
99999 END

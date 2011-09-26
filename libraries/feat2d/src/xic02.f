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
* XIC02n                                                               *
*                                                                      *
* Purpose  Call IC02n                                                  *
*                                                                      *
* Subroutines/functions called   ZTYPE, IC02n                          *
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
* For the description of the remaining input parameters see IC02n      *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 -170  X or B are not single precision                *
*                                                                      *
************************************************************************
C
      SUBROUTINE XIC027(LA,LCOL,LLD,LX,LB,NEQ,NIT,ITE,EPS,OMEGA)
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
      SUB='XIC027'
      IF (ICHECK.GE.997) CALL OTRC('XIC027','01/02/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.2.OR.ITYPE2.NE.2.OR.ITYPE3.NE.2) THEN
       CALL WERR(-170,'XIC027')
       GOTO 99999
      ENDIF
C
      CALL IC027(VWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),
     *           VWORK(L(LX)),VWORK(L(LB)),NEQ,NIT,ITE,EPS,OMEGA)
C
99999 END
C
C
C
      SUBROUTINE XIC02A(LA,LCOL,LLD,LOP,LX,LB,NEQ,NIT,ITE,EPS,OMEGA)
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
      SUB='XIC02A'
      IF (ICHECK.GE.997) CALL OTRC('XIC02A','01/02/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.2.OR.ITYPE2.NE.2.OR.ITYPE3.NE.2) THEN
       CALL WERR(-170,'XIC02A')
       GOTO 99999
      ENDIF
C
      CALL IC02A(VWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),
     *           KWORK(L(LOP)),VWORK(L(LX)),VWORK(L(LB)),NEQ,
     *           NIT,ITE,EPS,OMEGA)
C
99999 END

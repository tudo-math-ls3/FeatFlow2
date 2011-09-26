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
* XIE01n                                                               *
*                                                                      *
* Purpose  Allocate Workspace vector on DWORK                          *
*          Double/double precision version                             *
*          Call IE010                                                  *
*                                                                      *
* Subroutines/functions called   ZTYPE, ZNEW, ZDISP, IE010             *
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
* OMEGA    R*8    Determines preconditioning technique                 *
*                 0 < OMEGA      No Preconditioning                    *
*                 0 = OMEGA      Scaling using diagonal entries        *
*                 0 < OMEGA < 2  SSOR-Preconditioning                  *
*                 2 < OMEGA      Use external subroutine DCG0C         *
* DCG0C    SUBR   EXTERNAL Subroutine DCG0C(DG,NEQ)                    *
*                 Results  DG := C**(-1) * DG  for the precondioning   *
*                 matrix C                                             *
* For the description of the remaining input parameters see IE010      *
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
      SUBROUTINE XIE013(LA,LDIA,LDIAS,NDIA,LX,LB,NEQ,NIT,
     *                  ITE,EPS,OMEGA,DCG0C)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL DCG0C,YLAX13,YIA113
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE013'
      IF (ICHECK.GE.997) CALL OTRC('XIE013','02/01/91')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.1.OR.ITYPE2.NE.1.OR.ITYPE3.NE.1) THEN
       CALL WERR(-170,'XIE013')
       GOTO 99999
      ENDIF
C
      IREQ=4*NEQ
      BNOCON=OMEGA.LT.0D0
      IF (BNOCON) IREQ=3*NEQ
      IREQ=MAX(IREQ,4)
      CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
      L1=L(LWORK)
      L2=L1+NEQ
      L3=L2+NEQ
      L4=L3+NEQ
      IF (BNOCON) L4=L1
C
      KXYPAR(1)=L(LA)
      KXYPAR(2)=L(LDIA)
      KXYPAR(3)=L(LDIAS)
      KXYPAR(4)=NDIA
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX13,YIA113,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX13,YIA113,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ELSE
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX13,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ENDIF
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
99999 END
C
C
C
      SUBROUTINE XIE014(LA,LDIA,LDIAS,NDIA,LX,LB,NEQ,NIT,
     *                  ITE,EPS,OMEGA,DCG0C)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL DCG0C,YLAX14,YIA113
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE014'
      IF (ICHECK.GE.997) CALL OTRC('XIE014','02/01/91')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.1.OR.ITYPE2.NE.1.OR.ITYPE3.NE.1) THEN
       CALL WERR(-170,'XIE014')
       GOTO 99999
      ENDIF
C
      IREQ=4*NEQ
      BNOCON=OMEGA.LT.0D0
      IF (BNOCON) IREQ=3*NEQ
      IREQ=MAX(IREQ,4)
      CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
      L1=L(LWORK)
      L2=L1+NEQ
      L3=L2+NEQ
      L4=L3+NEQ
      IF (BNOCON) L4=L1
C
      KXYPAR(1)=L(LA)
      KXYPAR(2)=L(LDIA)
      KXYPAR(3)=L(LDIAS)
      KXYPAR(4)=NDIA
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX14,YIA113,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX14,YIA113,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ELSE
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX14,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ENDIF
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
99999 END
C
C
C
      SUBROUTINE XIE017(LA,LCOL,LLD,LX,LB,NEQ,NIT,ITE,EPS,OMEGA,DCG0C)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL DCG0C,YLAX17,YIA117,YID117
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE017'
      IF (ICHECK.GE.997) CALL OTRC('XIE017','01/02/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.1.OR.ITYPE2.NE.1.OR.ITYPE3.NE.1) THEN
       CALL WERR(-170,'XIE017')
       GOTO 99999
      ENDIF
C
      IREQ=4*NEQ
      BNOCON=OMEGA.LT.0D0
      IF (BNOCON) IREQ=3*NEQ
      IREQ=MAX(IREQ,4)
      CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
      L1=L(LWORK)
      L2=L1+NEQ
      L3=L2+NEQ
      L4=L3+NEQ
      IF (BNOCON) L4=L1
C
      KXYPAR(1)=L(LA)
      KXYPAR(2)=L(LCOL)
      KXYPAR(3)=L(LLD)
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX17,YIA117,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX17,YID117,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ELSE
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX17,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ENDIF
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
99999 END
C
C
C
      SUBROUTINE XIE018(LA,LCOL,LLD,LX,LB,NEQ,NIT,ITE,EPS,OMEGA,DCG0C)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL DCG0C,YLAX18,YIA117,YID118
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE018'
      IF (ICHECK.GE.997) CALL OTRC('XIE018','01/02/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.1.OR.ITYPE2.NE.1.OR.ITYPE3.NE.1) THEN
       CALL WERR(-170,'XIE018')
       GOTO 99999
      ENDIF
C
      IREQ=4*NEQ
      BNOCON=OMEGA.LT.0D0
      IF (BNOCON) IREQ=3*NEQ
      IREQ=MAX(IREQ,4)
      CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
      L1=L(LWORK)
      L2=L1+NEQ
      L3=L2+NEQ
      L4=L3+NEQ
      IF (BNOCON) L4=L1
C
      KXYPAR(1)=L(LA)
      KXYPAR(2)=L(LCOL)
      KXYPAR(3)=L(LLD)
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX18,YIA117,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX18,YID118,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ELSE
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX18,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ENDIF
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
99999 END
C
C
C
      SUBROUTINE XIE01A(LA,LCOL,LLD,LOP,LX,LB,NEQ,NIT,ITE,EPS,OMEGA,
     *                  DCG0C)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL DCG0C,YLAX1A,YIA11A,YID11A
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE01A'
      IF (ICHECK.GE.997) CALL OTRC('XIE01A','01/02/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.1.OR.ITYPE2.NE.1.OR.ITYPE3.NE.1) THEN
       CALL WERR(-170,'XIE01A')
       GOTO 99999
      ENDIF
C
      IREQ=4*NEQ
      BNOCON=OMEGA.LT.0D0
      IF (BNOCON) IREQ=3*NEQ
      IREQ=MAX(IREQ,4)
      CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
      L1=L(LWORK)
      L2=L1+NEQ
      L3=L2+NEQ
      L4=L3+NEQ
      IF (BNOCON) L4=L1
C
      KXYPAR(1)=L(LA)
      KXYPAR(2)=L(LCOL)
      KXYPAR(3)=L(LLD)
      KXYPAR(4)=L(LOP)
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX1A,YIA11A,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX1A,YID11A,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ELSE
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX1A,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
C
      ENDIF
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
99999 END

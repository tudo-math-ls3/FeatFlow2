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
* XIE33n                                                               *
*                                                                      *
* Purpose  Allocate Workspace vector on DWORK                          *
*          Single/double precision version                             *
*          Call IE010                                                  *
*                                                                      *
* Subroutines/functions called   ZTYPE, ZNEW, ZDISP, IE010             *
*                                                                      *
* Version from  11/12/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Numbers of the arrays describing the matrix in       *
* LCOL     I*4    Storage technique n                                  *
* LLD      I*4    LA  refers to a single precision matrix              *
* LOP      I*4                                                         *
* LX,LB    I*4    Numbers of the solution and the right hand side      *
* OMEGA    R*8    Determines preconditioning technique                 *
*                 0 < OMEGA      No Preconditioning                    *
*                 0 = OMEGA      Scaling using diagonal entries        *
*                 0 < OMEGA < 2  SSOR-Preconditioning                  *
*                 2 < OMEGA      Use external subroutine DCG0C         *
* DCG0C    SUBR   EXTERNAL Subroutine DCG0C(DG,NEQ,OMEGA)              *
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
      SUBROUTINE XIE333(LA,LDIA,LDIAS,NDIA,LX,LB,NEQ,NIT,
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
      EXTERNAL DCG0C,YLAX33,YIA133
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE333'
      IF (ICHECK.GE.997) CALL OTRC('XIE333','02/01/91')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.1.OR.ITYPE2.NE.1.OR.ITYPE3.NE.2) THEN
       CALL WERR(-170,'XIE333')
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
     *            NEQ,NIT,ITE,EPS,YLAX33,YIA133,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX33,YIA133,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
C
      ELSE
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX33,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
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
      SUBROUTINE XIE334(LA,LDIA,LDIAS,NDIA,LX,LB,NEQ,NIT,
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
      EXTERNAL DCG0C,YLAX34,YIA133
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE334'
      IF (ICHECK.GE.997) CALL OTRC('XIE334','02/01/91')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.1.OR.ITYPE2.NE.1.OR.ITYPE3.NE.2) THEN
       CALL WERR(-170,'XIE334')
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
     *            NEQ,NIT,ITE,EPS,YLAX34,YIA133,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX34,YIA133,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
C
      ELSE
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX34,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
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
      SUBROUTINE XIE337(LA,LCOL,LLD,LX,LB,NEQ,NIT,ITE,EPS,OMEGA,DCG0C)
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
      EXTERNAL DCG0C,YLAX37,YIA137,YID137
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE337'
      IF (ICHECK.GE.997) CALL OTRC('XIE337','11/12/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.1.OR.ITYPE2.NE.1.OR.ITYPE3.NE.2) THEN
       CALL WERR(-170,'XIE337')
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
     *            NEQ,NIT,ITE,EPS,YLAX37,YIA137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX37,YID137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
C
      ELSE
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX37,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
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
      SUBROUTINE XIE338(LA,LCOL,LLD,LX,LB,NEQ,NIT,ITE,EPS,OMEGA,DCG0C)
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
      EXTERNAL DCG0C,YLAX38,YIA137,YID138
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE338'
      IF (ICHECK.GE.997) CALL OTRC('XIE338','11/12/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.1.OR.ITYPE2.NE.1.OR.ITYPE3.NE.2) THEN
       CALL WERR(-170,'XIE338')
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
     *            NEQ,NIT,ITE,EPS,YLAX38,YIA137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX38,YID138,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
C
      ELSE
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX38,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
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
      SUBROUTINE XIE33A(LA,LCOL,LLD,LOP,LX,LB,NEQ,NIT,ITE,EPS,OMEGA,
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
      EXTERNAL DCG0C,YLAX3A,YIA13A,YID13A
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE33A'
      IF (ICHECK.GE.997) CALL OTRC('XIE33A','11/12/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE1.NE.1.OR.ITYPE2.NE.1.OR.ITYPE3.NE.2) THEN
       CALL WERR(-170,'XIE33A')
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
     *            NEQ,NIT,ITE,EPS,YLAX3A,YIA13A,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX3A,YID13A,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
C
      ELSE
C
       CALL IE010(DWORK(L(LX)),DWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX3A,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1)
C
      ENDIF
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
99999 END

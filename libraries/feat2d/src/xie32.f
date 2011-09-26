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
* XIE32n                                                               *
*                                                                      *
* Purpose  Allocate Workspace vector on DWORK                          *
*          Call IE020                                                  *
*                                                                      *
* Subroutines/functions called   ZTYPE, ZNEW, ZDISP, IE020             *
*                                                                      *
* Version from  11/12/89                                               *
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
*                 2 < OMEGA      Use external subroutine  CG0C         *
* CG0C     SUBR   EXTERNAL Subroutine  CG0C(VG,NEQ,OMEGA)              *
*                 Results  VG := C**(-1) * VG  for the precondioning   *
*                 matrix C                                             *
* For the description of the remaining input parameters see IE020      *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 -170 X or B are not single precision                 *
*                                                                      *
************************************************************************
C
      SUBROUTINE XIE323(LA,LDIA,LDIAS,NDIA,LX,LB,NEQ,NIT,
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
      EXTERNAL DCG0C,YLAX23,YIA123
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE323'
      IF (ICHECK.GE.997) CALL OTRC('XIE323','02/01/91')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      IF (ITYPE1.NE.ITYPE2.OR.ITYPE1.NE.2) THEN
       CALL WERR(-170,'XIE323')
       GOTO 99999
      ENDIF
C
      IREQ=4*NEQ
      BNOCON=OMEGA.LT.0D0
      IF (BNOCON) IREQ=3*NEQ
      IREQ=MAX(IREQ,4)
      CALL ZNEW(IREQ,-2,LWORK,'WORKCG')
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
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX23,YIA123,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX23,YIA123,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
C
      ELSE
C
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX23,DCG0C,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
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
      SUBROUTINE XIE324(LA,LDIA,LDIAS,NDIA,LX,LB,NEQ,NIT,
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
      EXTERNAL DCG0C,YLAX24,YIA123
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE324'
      IF (ICHECK.GE.997) CALL OTRC('XIE324','02/01/91')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      IF (ITYPE1.NE.ITYPE2.OR.ITYPE1.NE.2) THEN
       CALL WERR(-170,'XIE324')
       GOTO 99999
      ENDIF
C
      IREQ=4*NEQ
      BNOCON=OMEGA.LT.0D0
      IF (BNOCON) IREQ=3*NEQ
      IREQ=MAX(IREQ,4)
      CALL ZNEW(IREQ,-2,LWORK,'WORKCG')
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
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX24,YIA123,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX24,YIA123,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
C
      ELSE
C
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX24,DCG0C,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
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
      SUBROUTINE XIE327(LA,LCOL,LLD,LX,LB,NEQ,NIT,ITE,EPS,OMEGA,CG0C)
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
      EXTERNAL CG0C,YLAX27,YIA127,YID127
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE327'
      IF (ICHECK.GE.997) CALL OTRC('XIE327','11/12/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      IF (ITYPE1.NE.ITYPE2.OR.ITYPE1.NE.2) THEN
       CALL WERR(-170,'XIE327')
       GOTO 99999
      ENDIF
C
      IREQ=4*NEQ
      BNOCON=OMEGA.LT.0D0
      IF (BNOCON) IREQ=3*NEQ
      IREQ=MAX(IREQ,4)
      CALL ZNEW(IREQ,-2,LWORK,'WORKCG')
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
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX27,YIA127,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX27,YID127,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
C
      ELSE
C
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX27,CG0C,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
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
      SUBROUTINE XIE328(LA,LCOL,LLD,LX,LB,NEQ,NIT,ITE,EPS,OMEGA,CG0C)
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
      EXTERNAL CG0C,YLAX28,YIA127,YID128
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE328'
      IF (ICHECK.GE.997) CALL OTRC('XIE328','11/12/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      IF (ITYPE1.NE.ITYPE2.OR.ITYPE1.NE.2) THEN
       CALL WERR(-170,'XIE328')
       GOTO 99999
      ENDIF
C
      IREQ=4*NEQ
      BNOCON=OMEGA.LT.0D0
      IF (BNOCON) IREQ=3*NEQ
      IREQ=MAX(IREQ,4)
      CALL ZNEW(IREQ,-2,LWORK,'WORKCG')
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
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX28,YIA127,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX28,YID128,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
C
      ELSE
C
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX28,CG0C,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
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
      SUBROUTINE XIE32A(LA,LCOL,LLD,LOP,LX,LB,NEQ,NIT,ITE,EPS,OMEGA,
     *                  CG0C)
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
      EXTERNAL CG0C,YLAX2A,YIA12A,YID12A
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='XIE32A'
      IF (ICHECK.GE.997) CALL OTRC('XIE32A','11/12/89')
C
      CALL ZTYPE(LX,ITYPE1)
      CALL ZTYPE(LB,ITYPE2)
      IF (ITYPE1.NE.ITYPE2.OR.ITYPE1.NE.2) THEN
       CALL WERR(-170,'XIE32A')
       GOTO 99999
      ENDIF
C
      IREQ=4*NEQ
      BNOCON=OMEGA.LT.0D0
      IF (BNOCON) IREQ=3*NEQ
      IREQ=MAX(IREQ,4)
      CALL ZNEW(IREQ,-2,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
      L1=L(LWORK)
      L2=L1+NEQ
      L3=L2+NEQ
      L4=L3+NEQ
      IF (OMEGA.LT.0D0) L4=L1
C
      KXYPAR(1)=L(LA)
      KXYPAR(2)=L(LCOL)
      KXYPAR(3)=L(LLD)
      KXYPAR(4)=L(LOP)
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
C
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX2A,YIA12A,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX2A,YID12A,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
C
      ELSE
C
       CALL IE020(VWORK(L(LX)),VWORK(L(LB)),
     *            NEQ,NIT,ITE,EPS,YLAX2A,CG0C,BNOCON,
     *            VWORK(L1),VWORK(L2),VWORK(L3),VWORK(L4),1)
C
      ENDIF
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
99999 END

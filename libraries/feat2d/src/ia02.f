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
* IA02n                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          Jacobi-method                                               *
*          Single/single precision version                             *
*                                                                      *
* Subroutines/functions called  LLI2, LCL2, LCP2, LLC2, LVM2           *
*                                                                      *
* Version from  02/06/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix                                               *
* KCOL     I*4    Pointer vectors for the matrix VA corresponding to   *
* KLD      I*4    the storage technique n                              *
* KOP      I*4                                                         *
* VX       R*4    Starting vector                                      *
* VB       R*4    Right hand side                                      *
* VD       R*4    Help vector                                          *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Maximum number of iterations                         *
* EPS      R*8    Desired precision                                    *
* OMEGA    R*8    Relaxation parameter                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VX       R*4    Solution vector                                      *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 +1  Precision EPS not achieved after NIT iterations  *
*                                                                      *
************************************************************************
C
	SUBROUTINE IA023(VA,KDIA,KDIAS,NDIA,VX,VB,VD,
     *                 NEQ,NIT,ITE,EPS,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KDIA(*),KDIAS(*),VX(*),VB(*),VD(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='IA023'
      IF (ICHECK.GE.997) CALL OTRC('IA023 ','02/06/91')
C
      BMSG2=M.GE.2.OR.MT.GE.2
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IA023 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      DO 1 IEQ=1,NEQ
1     VA(IEQ)=1./VA(IEQ)
C
      DO 2 ITE=1,NIT
      DM1=0D0
      DM2=0D0
      CALL LCP2(VB,VD,NEQ)
      DO 3 IDIA=2,NDIA
      J1=KDIA(IDIA)
      IF (J1.GT.0) THEN
       I1=1
       NEQ1=NEQ-J1
      ELSE
       I1=1-J1
       NEQ1=NEQ+J1
      ENDIF
      J0=KDIAS(IDIA)-I1
      CALL LVM2(VX(I1+J1),VA(I1+J0),VD(I1),NEQ1,-1D0,1D0)
3     CONTINUE
      CALL LVM2(VD,VA,VD,NEQ,1D0,0D0)
      CALL LLC2(VX,VD,NEQ,-OMEGA,OMEGA)
      CALL LLI2(VD,NEQ,DM1,IND)
      CALL LLC2(VX,VD,NEQ,1D0,1D0)
      CALL LLI2(VD,NEQ,DM2,IND)
      CALL LCP2(VD,VX,NEQ)
      IF (DM1.LE.EPS*DM2) GOTO 99998
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
       CALL OMSG(74,'IA023 ')
      ENDIF
2     CONTINUE
C
      IER=1
      DO 10 IEQ=1,NEQ
10    VA(IEQ)=1./VA(IEQ)
      WRITE (CPARAM,'(I15,D25.16)') NIT,DM1/DM2
      CALL OMSG(71,'IA023 ')
      CALL OMSG(75,'IA023 ')
      GOTO 99999
C
99998 IER=0
      DO 11 IEQ=1,NEQ
11    VA(IEQ)=1./VA(IEQ)
      WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
      CALL OMSG(75,'IA023 ')
C
99999 END
C
C
C
      SUBROUTINE IA027(VA,KCOL,KLD,VX,VB,VD,
     *                 NEQ,NIT,ITE,EPS,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),VB(*),VD(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='IA027'
      IF (ICHECK.GE.997) CALL OTRC('IA027 ','02/06/91')
C
      BMSG2=M.GE.2.OR.MT.GE.2
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IA027 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      DO 2 ITE=1,NIT
      DM1=0D0
      DM2=0D0
      CALL LCP2(VB,VD,NEQ)
      DO 3 IEQ=1,NEQ
      DO 4 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
4     VD(IEQ)=VD(IEQ)-VA(ICOL)*VX(KCOL(ICOL))
3     CONTINUE
      DO 5 IEQ=1,NEQ
      VD(IEQ)=(1D0-OMEGA)*VX(IEQ)+OMEGA*VD(IEQ)/VA(KLD(IEQ))
      DM1=MAX(DM1,DBLE(ABS(VD(IEQ)-VX(IEQ))))
      DM2=MAX(DM2,DBLE(ABS(VD(IEQ))))
5	VX(IEQ)=VD(IEQ)
      IF (DM1.LE.EPS*DM2) GOTO 99998
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
       CALL OMSG(74,'IA027 ')
      ENDIF
2     CONTINUE
C
      IER=1
      WRITE (CPARAM,'(I15,D25.16)') NIT,DM1/DM2
      CALL OMSG(71,'IA027 ')
      CALL OMSG(75,'IA027 ')
      GOTO 99999
C
99998 IER=0
      WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
      CALL OMSG(75,'IA027 ')
C
99999 END
C
C
C
      SUBROUTINE IA02A(VA,KCOL,KLD,KOP,VX,VB,VD,
     *                 NEQ,NIT,ITE,EPS,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),VX(*),VB(*),VD(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='IA02A'
      IF (ICHECK.GE.997) CALL OTRC('IA02A ','02/06/91')
C
      BMSG2=M.GE.2.OR.MT.GE.2
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IA02A ')
        GOTO 99999
       ENDIF
      ENDIF
C
      DO 2 ITE=1,NIT
      DM1=0D0
      DM2=0D0
      CALL LCP2(VB,VD,NEQ)
      DO 3 IEQ=1,NEQ
      JOP=KOP(IEQ)
      DO 4 ICOL=KLD(JOP)+1,KLD(JOP+1)-1
4     VD(IEQ)=VD(IEQ)-VA(ICOL)*VX(KCOL(ICOL)+IEQ)
3     CONTINUE
      DO 5 IEQ=1,NEQ
      VD(IEQ)=(1D0-OMEGA)*VX(IEQ)+OMEGA*VD(IEQ)/VA(KLD(KOP(IEQ)))
      DM1=MAX(DM1,DBLE(ABS(VD(IEQ)-VX(IEQ))))
      DM2=MAX(DM2,DBLE(ABS(VD(IEQ))))
5	VX(IEQ)=VD(IEQ)
      IF (DM1.LE.EPS*DM2) GOTO 99998
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
       CALL OMSG(74,'IA02A ')
      ENDIF
2     CONTINUE
C
      IER=1
      WRITE (CPARAM,'(I15,D25.16)') NIT,DM1/DM2
      CALL OMSG(71,'IA02A ')
      CALL OMSG(75,'IA02A ')
      GOTO 99999
C
99998 IER=0
      WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
      CALL OMSG(75,'IA02A ')
C
99999 END

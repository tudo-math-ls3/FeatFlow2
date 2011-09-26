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
* IA03n                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          Jacobi-method                                               *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called  LLI1, LCL1, LCP1, LLC1, LVM3           *
*                                                                      *
* Version from  02/06/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix                                               *
* KCOL     I*4    Pointer vectors for the matrix VA corresponding to   *
* KLD      I*4    the storage technique n                              *
* KOP      I*4                                                         *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* DD       R*8    Help vector                                          *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Maximum number of iterations                         *
* EPS      R*8    Desired precision                                    *
* OMEGA    R*8    Relaxation parameter                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Solution vector                                      *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 +1  Precision EPS not achieved after NIT iterations  *
*                                                                      *
************************************************************************
C
C
	SUBROUTINE IA033(VA,KDIA,KDIAS,NDIA,DX,DB,DD,
     *                 NEQ,NIT,ITE,EPS,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KDIA(*),KDIAS(*),DX(*),DB(*),DD(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='IA033'
      IF (ICHECK.GE.997) CALL OTRC('IA033 ','02/06/91')
C
      BMSG2=M.GE.2.OR.MT.GE.2
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IA033 ')
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
      CALL LCP1(DB,DD,NEQ)
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
      CALL LVM3(DX(I1+J1),VA(I1+J0),DD(I1),NEQ1,-1D0,1D0)
3     CONTINUE
      CALL LVM3(DD,VA,DD,NEQ,1D0,0D0)
      CALL LLC1(DX,DD,NEQ,-OMEGA,OMEGA)
      CALL LLI1(DD,NEQ,DM1,IND)
      CALL LLC1(DX,DD,NEQ,1D0,1D0)
      CALL LLI1(DD,NEQ,DM2,IND)
      CALL LCP1(DD,DX,NEQ)
      IF (DM1.LE.EPS*DM2) GOTO 99998
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
       CALL OMSG(74,'IA033 ')
      ENDIF
2     CONTINUE
C
      IER=1
      DO 10 IEQ=1,NEQ
10    VA(IEQ)=1./VA(IEQ)
      WRITE (CPARAM,'(I15,D25.16)') NIT,DM1/DM2
      CALL OMSG(71,'IA033 ')
      CALL OMSG(75,'IA033 ')
      GOTO 99999
C
99998 IER=0
      DO 11 IEQ=1,NEQ
11    VA(IEQ)=1./VA(IEQ)
      WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
      CALL OMSG(75,'IA033 ')
C
99999 END
C
C
C
      SUBROUTINE IA037(VA,KCOL,KLD,DX,DB,DD,
     *                 NEQ,NIT,ITE,EPS,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='IA037'
      IF (ICHECK.GE.997) CALL OTRC('IA037 ','02/06/91')
C
      BMSG2=M.GE.2.OR.MT.GE.2
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IA037 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      DO 2 ITE=1,NIT
      DM1=0D0
      DM2=0D0
      CALL LCP1(DB,DD,NEQ)
      DO 3 IEQ=1,NEQ
      DO 4 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
4     DD(IEQ)=DD(IEQ)-DBLE(VA(ICOL))*DX(KCOL(ICOL))
3     CONTINUE
      DO 5 IEQ=1,NEQ
      DD(IEQ)=(1D0-OMEGA)*DX(IEQ)+OMEGA*DD(IEQ)/DBLE(VA(KLD(IEQ)))
      DM1=MAX(DM1,ABS(DD(IEQ)-DX(IEQ)))
      DM2=MAX(DM2,ABS(DD(IEQ)))
5	DX(IEQ)=DD(IEQ)
      IF (DM1.LE.EPS*DM2) GOTO 99998
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
       CALL OMSG(74,'IA037 ')
      ENDIF
2     CONTINUE
C
      IER=1
      WRITE (CPARAM,'(I15,D25.16)') NIT,DM1/DM2
      CALL OMSG(71,'IA037 ')
      CALL OMSG(75,'IA037 ')
      GOTO 99999
C
99998 IER=0
      WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
      CALL OMSG(75,'IA037 ')
C
99999 END
C
C
C
      SUBROUTINE IA03A(VA,KCOL,KLD,KOP,DX,DB,DD,
     *                 NEQ,NIT,ITE,EPS,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),DX(*),DB(*),DD(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='IA03A'
      IF (ICHECK.GE.997) CALL OTRC('IA03A ','02/06/91')
C
      BMSG2=M.GE.2.OR.MT.GE.2
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IA03A ')
        GOTO 99999
       ENDIF
      ENDIF
C
      DO 2 ITE=1,NIT
      DM1=0D0
      DM2=0D0
      CALL LCP1(DB,DD,NEQ)
      DO 3 IEQ=1,NEQ
      JOP=KOP(IEQ)
      DO 4 ICOL=KLD(JOP)+1,KLD(JOP+1)-1
4     DD(IEQ)=DD(IEQ)-DBLE(VA(ICOL))*DX(KCOL(ICOL)+IEQ)
3     CONTINUE
      DO 5 IEQ=1,NEQ
      DD(IEQ)=(1D0-OMEGA)*DX(IEQ)+OMEGA*DD(IEQ)/DBLE(VA(KLD(KOP(IEQ))))
      DM1=MAX(DM1,ABS(DD(IEQ)-DX(IEQ)))
      DM2=MAX(DM2,ABS(DD(IEQ)))
5	DX(IEQ)=DD(IEQ)
      IF (DM1.LE.EPS*DM2) GOTO 99998
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
       CALL OMSG(74,'IA03A ')
      ENDIF
2     CONTINUE
C
      IER=1
      WRITE (CPARAM,'(I15,D25.16)') NIT,DM1/DM2
      CALL OMSG(71,'IA03A ')
      CALL OMSG(75,'IA03A ')
      GOTO 99999
C
99998 IER=0
      WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
      CALL OMSG(75,'IA03A ')
C
99999 END

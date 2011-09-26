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
* IC02n                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using SOR-method      *
*          Single/single precision version                             *
*                                                                      *
* Subroutines/functions called  LLI2, LCL2                             *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix                                               *
* KCOL     I*4    Pointer vectors for the matrix VA corresponding to   *
* KLD      I*4    the storage technique n                              *
* KOP      I*4                                                         *
* VX       R*4    Starting vector                                      *
* VB       R*4    Right hand side                                      *
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
      SUBROUTINE IC027(VA,KCOL,KLD,VX,VB,NEQ,NIT,ITE,EPS,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*),VX(*),VB(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='IC027'
      IF (ICHECK.GE.997) CALL OTRC('IC027 ','01/02/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      DM1=0D0
      DM2=1D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IC027 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      DO 2 ITE=1,NIT
      DM1=0D0
      DM2=0D0
      DO 3 IEQ=1,NEQ
      AUX=VB(IEQ)
      DO 4 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
4     AUX=AUX-VA(ICOL)*VX(KCOL(ICOL))
      AUX=OMEGA*(AUX/VA(KLD(IEQ))-VX(IEQ))+VX(IEQ)
      DM1=MAX(DM1,ABS(AUX-VX(IEQ)))
      DM2=MAX(DM2,ABS(AUX))
3     VX(IEQ)=AUX
      IF (DM1.LE.EPS*DM2) GOTO 99998
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
       CALL OMSG(74,'IC027 ')
      ENDIF
2     CONTINUE
C
      IER=1
      WRITE (CPARAM,'(I15,D25.16)') NIT,DM1/DM2
      CALL OMSG(71,'IC027 ')
      CALL OMSG(75,'IC027 ')
      GOTO 99999
C
99998 IER=0
      WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
      CALL OMSG(75,'IC027 ')
C
99999 END
C
C
C
      SUBROUTINE IC02A(VA,KCOL,KLD,KOP,VX,VB,NEQ,NIT,ITE,EPS,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),VX(*),VB(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='IC02A'
      IF (ICHECK.GE.997) CALL OTRC('IC02A ','01/02/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      DM1=0D0
      DM2=1D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI2(VB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL2(VX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IC02A ')
        GOTO 99999
       ENDIF
      ENDIF
C
      DO 2 ITE=1,NIT
      DM1=0D0
      DM2=0D0
      DO 3 IEQ=1,NEQ
      AUX=VB(IEQ)
      JOP=KOP(IEQ)
      DO 4 ICOL=KLD(JOP)+1,KLD(JOP+1)-1
4     AUX=AUX-VA(ICOL)*VX(KCOL(ICOL)+IEQ)
      AUX=OMEGA*(AUX/VA(KLD(JOP))-VX(IEQ))+VX(IEQ)
      DM1=MAX(DM1,ABS(AUX-VX(IEQ)))
      DM2=MAX(DM2,ABS(AUX))
3     VX(IEQ)=AUX
      IF (DM1.LE.EPS*DM2) GOTO 99998
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
       CALL OMSG(74,'IC02A ')
      ENDIF
2     CONTINUE
C
      IER=1
      WRITE (CPARAM,'(I15,D25.16)') NIT,DM1/DM2
      CALL OMSG(71,'IC02A ')
      CALL OMSG(75,'IC02A ')
      GOTO 99999
C
99998 IER=0
      WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
      CALL OMSG(75,'IC02A ')
C
99999 END

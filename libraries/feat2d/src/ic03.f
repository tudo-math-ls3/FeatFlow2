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
* IC03n                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using SOR-method      *
*          Single/double precision version                             *
*                                                                      *
* Subroutines/functions called  LLI1, LCL1                             *
*                                                                      *
* Version from  01/09/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VA       R*4    Matrix                                               *
* KCOL     I*4    Pointer vectors for the matrix VA corresponding to   *
* KLD      I*4    the storage technique n                              *
* KOP      I*4                                                         *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
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
      SUBROUTINE IC037(VA,KCOL,KLD,DX,DB,NEQ,NIT,ITE,EPS,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*),DX(*),DB(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='IC037'
      IF (ICHECK.GE.997) CALL OTRC('IC037 ','01/09/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      DM1=0D0
      DM2=1D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IC037 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      DO 2 ITE=1,NIT
      DM1=0D0
      DM2=0D0
      DO 3 IEQ=1,NEQ
      AUX=DB(IEQ)
      DO 4 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
4     AUX=AUX-DBLE(VA(ICOL))*DX(KCOL(ICOL))
      AUX=OMEGA*(AUX/DBLE(VA(KLD(IEQ)))-DX(IEQ))+DX(IEQ)
      DM1=MAX(DM1,ABS(AUX-DX(IEQ)))
      DM2=MAX(DM2,ABS(AUX))
3     DX(IEQ)=AUX
      IF (DM1.LE.EPS*DM2) GOTO 99998
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
       CALL OMSG(74,'IC037 ')
      ENDIF
2     CONTINUE
C
      IER=1
      WRITE (CPARAM,'(I15,D25.16)') NIT,DM1/DM2
      CALL OMSG(71,'IC037 ')
      CALL OMSG(75,'IC037 ')
      GOTO 99999
C
99998 IER=0
      WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
      CALL OMSG(75,'IC037 ')
C
99999 END
C
C
C
      SUBROUTINE IC03A(VA,KCOL,KLD,KOP,DX,DB,NEQ,NIT,ITE,EPS,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VA(*),KCOL(*),KLD(*),KOP(*),DX(*),DB(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='IC03A'
      IF (ICHECK.GE.997) CALL OTRC('IC03A ','01/09/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      DM1=0D0
      DM2=1D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IC03A ')
        GOTO 99999
       ENDIF
      ENDIF
C
      DO 2 ITE=1,NIT
      DM1=0D0
      DM2=0D0
      DO 3 IEQ=1,NEQ
      AUX=DB(IEQ)
      JOP=KOP(IEQ)
      DO 4 ICOL=KLD(JOP)+1,KLD(JOP+1)-1
4     AUX=AUX-DBLE(VA(ICOL))*DX(KCOL(ICOL)+IEQ)
      AUX=OMEGA*(AUX/DBLE(VA(KLD(JOP)))-DX(IEQ))+DX(IEQ)
      DM1=MAX(DM1,ABS(AUX-DX(IEQ)))
      DM2=MAX(DM2,ABS(AUX))
3     DX(IEQ)=AUX
      IF (DM1.LE.EPS*DM2) GOTO 99998
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
       CALL OMSG(74,'IC03A ')
      ENDIF
2     CONTINUE
C
      IER=1
      WRITE (CPARAM,'(I15,D25.16)') NIT,DM1/DM2
      CALL OMSG(71,'IC03A ')
      CALL OMSG(75,'IC03A ')
      GOTO 99999
C
99998 IER=0
      WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
      CALL OMSG(75,'IC03A ')
C
99999 END

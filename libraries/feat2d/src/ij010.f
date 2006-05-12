************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek, M. Koester         *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* IJ010 - BiCGStab as Coarse grid solver for MG                        *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          a preconditioned BI-CGSTAB method                           *
*          Double precision version                                    *
*                                                                      *
* This routine is exactly the same as II010. The background is the     *
* missing recursion in Multigrid: If MG is used as a preconditioner    *
* in a BiCGStab solver, there mustn't be the same BiCGStab solver      *
* routine be used in the coarse grid solver!                           *
* Therefore this routine should be used as a coarse grid solver while  *
* II010 should be used as a general BiCGStab solver.                   *
*                                                                      *
* Subroutines/functions called  LSP1 , LLC1 , LL21 , LLI1, LCL1        *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* ITE      I*4    Minimum number of iterations (used for BREL=1 only)  *
* NIT      I*4    Maximum number of iterations                         *
* EPS      R*8    Desired precision                                    *
*                 BREL=0: Stop if !!RES!! < EPS                        *
*                 BREL=1: Stop if !!RES!!/!!RES0!! < EPS               *
*                         and a minimum of ITE iterations is performed *
* DAX0     SUBR   EXTERNAL Subroutine DAX0(DX,DAX,NEQ,A1,A2)           *
*                 Results  DAX := A1 * A * DX + A2 * DAX               *
* DCG0C    SUBR   EXTERNAL Subroutine DCG0C(DG,NEQ)                    *
*                 Results  DG := C**(-1) * DG  for the precondioning   *
*                 matrix C                                             *
* BNOCON   L*4    .TRUE.   No Preconditioning                          *
* DR,DD    R*8    Workspace vectors of length NEQ                      *
* DD1,DG   R*8    For BNOCON , DG must be replaced by DR               *
* BREL     I*4    See above                                            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Solution vector                                      *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 +1  Precision EPS not achieved after NIT iterations  *
* RHO      R*8    Convergence rate                                     *
* DABSRS   R*8    Norm of last absolute residuum                       *
* DRELRD   R*8    Norm of last relative residuum                       *
* DEFINI   R*8    Norm of initial residuum                             *
* RHOASM   R*8    asymptotic convergence rate                          *
*                 calculated from the last three residuals             *
*                                                                      *
* 06.08.2003: Slight modification by Michael Koester:                  *
* Now BNOCON works correctly !                                         *
* 21.08.2003: Conv. rate and norm of residuals are now returned.       *
* 12.10.2004: Removed an unnecessary parameter DTMP                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE IJ010(DX,DB,NEQ,NIT,ITE,EPS,DAX0,DCG0C,BNOCON,
     *                 DR,DR0,DP,DPA,DSA,BREL,
     *                 RHO,DABSRS,DRELRS,DEFINI,RHOASM)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*),DB(*),DR(*),DR0(*),DP(*),DPA(*),DSA(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/

C Length of the queue of last residuals for the computation of
C the asymptotic convergence rate

      PARAMETER (IASRLN=3)

C The queue saves the current residual and the two previous
C and the two previous residuals

      DIMENSION RESQUE(IASRLN)
      DOUBLE PRECISION RESQUE

      SUB='II010'
      IF (ICHECK.GE.997) CALL OTRC('II010 ','04/14/93')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
C
C Initialize used vectors with zero
      CALL LCL1(DP,NEQ)
      CALL LCL1(DPA,NEQ)
C
      CALL LCP1(DB,DR,NEQ)
      IF (.NOT.BNOCON) CALL DCG0C(DR,NEQ)
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DR,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'II010 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      RHO0  =1D0
      DALPHA=1D0
      OMEGA0=1D0
C
      CALL LCP1(DB,DR,NEQ)
      CALL DAX0(DX,DR,NEQ,-1D0,1D0)
      IF (.NOT.BNOCON) CALL DCG0C(DR,NEQ)
      CALL LL21(DR,NEQ,RES)

C Michael Koester: Skalierung, damit der Vektor (1111...) die Norm 1 hat

      RES = RES / SQRT (DBLE(NEQ))
      IF (.NOT.((RES.GE.1D-99).AND.(RES.LE.1D99))) RES = 0D0

C Initialize starting residuum

      DEFINI = RES

C initialize the queue of the last residuals with RES
       DO I=1,IASRLN
         RESQUE(I)=RES
       END DO  

      IF (RES.LE.EPS.AND..NOT.BREL) THEN
       ITE=0
       FR=RES
       GOTO 200
      ENDIF
C
      CALL LCP1(DR,DR0,NEQ)

C      call lcp1 (db,dtmp,neq)
C      CALL DAX0(DX,Dtmp,NEQ,-1D0,1D0)
C      CALL LL21 (Dtmp,NEQ,FR2)
C      WRITE (*,'(A,I15,3D25.16)') 'Reference: ',ITE2,FR,FR2,RES
      

C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LSP1(DR0,DR,NEQ,RHO1)
      DBETA=(RHO1*DALPHA)/(RHO0*OMEGA0)
      RHO0 =RHO1
C
      CALL LLC1(DR ,DP,NEQ,1D0,DBETA)
      CALL LLC1(DPA,DP,NEQ,-DBETA*OMEGA0,1D0)
C
      CALL DAX0(DP,DPA,NEQ,1D0,0D0)
      IF (.NOT.BNOCON) CALL DCG0C(DPA,NEQ)
C
      CALL LSP1(DR0,DPA,NEQ,DALPHA)
      DALPHA=RHO1/DALPHA
C
      CALL LLC1(DPA,DR,NEQ,-DALPHA,1D0)
C
      CALL DAX0(DR,DSA,NEQ,1D0,0D0)
      IF (.NOT.BNOCON) CALL DCG0C(DSA,NEQ)
C
      CALL LSP1(DSA,DR ,NEQ,OMEGA1)
      CALL LSP1(DSA,DSA,NEQ,OMEGA2)
      IF (OMEGA1.EQ.0) THEN
        OMEGA0 = 0D0
      ELSE
        IF (OMEGA2.EQ.0) GOTO 150
        OMEGA0=OMEGA1/OMEGA2
      END IF
C
      CALL LLC1(DP ,DX ,NEQ,DALPHA,1D0)
      CALL LLC1(DR ,DX ,NEQ,OMEGA0,1D0)
C
      CALL LLC1(DSA,DR,NEQ,-OMEGA0,1D0)
C
C
      CALL LL21(DR,NEQ,FR)

C Michael Koester: Scaling for the vektor (1111...) to have norm = 1
      FR = FR / SQRT (DBLE(NEQ))
      IF (.NOT.((FR.GE.1D-99).AND.(FR.LE.1D99))) FR = 0D0
     
C shift the queue with the last residuals and add the new
C residual to it
       DO I=1,IASRLN-1
         RESQUE(I)=RESQUE(I+1)
       END DO  
       RESQUE(IASRLN)=FR

C      call lcp1 (db,dtmp,neq)
C      CALL DAX0(DX,Dtmp,NEQ,-1D0,1D0)
C      CALL LL21 (Dtmp,NEQ,FR2)
C      WRITE (*,'(A,I15,3D25.16)') 'Reference: ',ITE2,FR,FR2,RES
      
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'II010 ')
      ENDIF
      IF (BREL) THEN
       IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
      ELSE
       IF (FR.LE.EPS) GOTO 200
      ENDIF
C
100   CONTINUE
C
150   WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'II010 ')
      CALL OMSG(72,'II010 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/DBLE(NIT))
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'II010 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'II010 ')
C
C      call lcp1 (db,dp,neq)
C      CALL DAX0(DX,DP,NEQ,-1D0,1D0)
C      CALL LL21 (DP,NEQ,RES2)
C      CALL LL21 (DB,NEQ,FR2)
C      res3 = res2/fr2
C      droc = RES3**(1D0/DBLE(ITE))
C      print *,'res/fr/res2/roc:',res,fr,res,res2,droc

      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/DBLE(ITE))
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'II010 ')
C
99999 RHO = CAPPA
      DABSRS = FR
      DRELRS = RES
      
      RHOASM = RHO
      IF (RESQUE(1).GE.1D-70) THEN
        NITREL = MIN(ITE,IASRLN-1)
        RHOASM=(FR/RESQUE(1))**(1D0/DBLE(NITREL))
      END IF
      
      END


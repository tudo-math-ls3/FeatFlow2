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
* II210                                                                *
*                                                                      *
* Purpose  Smoothing using the BiCGStab method.                        *
*          Double precision version, Storage technique 7               *
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
* NIT      I*4    Number of smoothing steps                            *
* DAX0     SUBR   EXTERNAL Subroutine DAX0(DX,DAX,NEQ,A1,A2)           *
*                 Results  DAX := A1 * A * DX + A2 * DAX               *
* DCG0C    SUBR   EXTERNAL Subroutine DCG0C(DG,NEQ)                    *
*                 Results  DG := C**(-1) * DG  for the precondioning   *
*                 matrix C                                             *
* BNOCON   L*4    .TRUE.   No Preconditioning                          *
*                 CARE: in contrast to older IIxxx-routines            *
*                 this is working correctly here !!!                   *
* DR,DD    R*8    Workspace vectors of length NEQ                      *
* DD1,DG   R*8    For BNOCON , DG must be replaced by DR               *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8(NEQ) Smoothed vector                                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE II210(DX,DB,NEQ,NIT,DAX0,DCG0C,BNOCON,
     *                 DR,DR0,DP,DPA,DSA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*),DB(*),DR(*),DR0(*),DP(*),DPA(*),DSA(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='II010'
      IF (ICHECK.GE.997) CALL OTRC('II217 ','08/05/03')

C Initialize used vectors with zero
      CALL LCL1(DP,NEQ)
      CALL LCL1(DPA,NEQ)
C
      CALL LCP1(DB,DR,NEQ)
      IF (.NOT.BNOCON) CALL DCG0C(DR,NEQ)
C
C *** Initialization
      RHO0  =1D0
      DALPHA=1D0
      OMEGA0=1D0

      CALL LCP1(DB,DR,NEQ)
      CALL DAX0(DX,DR,NEQ,-1D0,1D0)
      IF (.NOT.BNOCON) CALL DCG0C(DR,NEQ)

      CALL LCP1(DR,DR0,NEQ)

C *** Iterative correction
      DO ITE=1,NIT

        CALL LSP1(DR0,DR,NEQ,RHO1)
        DBETA=(RHO1*DALPHA)/(RHO0*OMEGA0)
        RHO0 =RHO1

        CALL LLC1(DR ,DP,NEQ,1D0,DBETA)
        CALL LLC1(DPA,DP,NEQ,-DBETA*OMEGA0,1D0)

        CALL DAX0(DP,DPA,NEQ,1D0,0D0)
        IF (.NOT.BNOCON) CALL DCG0C(DPA,NEQ)

        CALL LSP1(DR0,DPA,NEQ,DALPHA)
C Cancel here if the norm is nearly zero
        IF (DALPHA.LT.1D-99) GOTO 99999
        DALPHA=RHO1/DALPHA

        CALL LLC1(DPA,DR,NEQ,-DALPHA,1D0)

        CALL DAX0(DR,DSA,NEQ,1D0,0D0)
        IF (.NOT.BNOCON) CALL DCG0C(DSA,NEQ)

        CALL LSP1(DSA,DR ,NEQ,OMEGA1)
        CALL LSP1(DSA,DSA,NEQ,OMEGA2)
C Cancel here if the norm is nearly zero
        IF (OMEGA2.LT.1D-99) GOTO 99999

        OMEGA0=OMEGA1/OMEGA2

        CALL LLC1(DP ,DX ,NEQ,DALPHA,1D0)
        CALL LLC1(DR ,DX ,NEQ,OMEGA0,1D0)

        CALL LLC1(DSA,DR,NEQ,-OMEGA0,1D0)
        
C        call ll21 (dr,neq,dnorm2)
C        print *,'Residuum norm: ',dnorm2
        
      END DO

99999 END

************************************************************************
*                                                                      *
* II217M                                                               *
*                                                                      *
* Purpose  Smoothing using the BiCGStab method.                        *
*          Double precision version, Storage technique 7               *
*                                                                      *
C This subroutine is neary the same like II017 except for it is used   *
C in multigrid as a smoother with calling conventon like XII217.       *
C It has no stopping tolerace parameters but uses only a maximum number*
C of iteration steps for smoothing.                                    *
C The smoothing is only aborted if the residuum is too small.          *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Handles of the arrays describing the matrix in       *
* LCOL     I*4    Storage technique n                                  *
* LLD      I*4                                                         *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Number of smoothing steps                            *
* OMEGA    R*8    Determines preconditioning technique                 *
*                 0 < OMEGA      No Preconditioning                    *
*                 0 = OMEGA      Scaling using diagonal entries        *
*                 0 < OMEGA < 2  SSOR-Preconditioning                  *
*                 2 < OMEGA      Use external subroutine DCG0C         *
* DCG0C    SUBR   EXTERNAL Subroutine DCG0C(DG,NEQ)                    *
*                 Results  DG := C**(-1) * DG  for the precondioning   *
*                 matrix C                                             *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8(NEQ) Smoothed vector                                    *
*                                                                      *
************************************************************************

      SUBROUTINE II217M(LA,LCOL,LLD,DX,DB,NEQ,NIT,OMEGA,DCG0C)
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
      DIMENSION DX(*),DB(*)
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='II217M'
      IF (ICHECK.GE.997) CALL OTRC('II217M','08/01/03')
C
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE3.NE.1) THEN
       CALL WERR(-170,'II217M')
       GOTO 99999
      ENDIF
C
      IREQ=5*NEQ
      BNOCON=OMEGA.LT.0D0
      IREQ=MAX(IREQ,5)
      CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
      L1=L(LWORK)
      L2=L1+NEQ
      L3=L2+NEQ
      L4=L3+NEQ
      L5=L4+NEQ
C
      KXYPAR(1)=L(LA)
      KXYPAR(2)=L(LCOL)
      KXYPAR(3)=L(LLD)
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
C
       CALL II210(DX,DB,
     *            NEQ,NIT,YLAX17,YIA117,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5))
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL II210(DX,DB,
     *            NEQ,NIT,YLAX17,YID117,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5))
C
      ELSE
C
       CALL II210(DX,DB,
     *            NEQ,NIT,YLAX17,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5))
C
      ENDIF
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
99999 END

************************************************************************
*                                                                      *
* II227M                                                               *
*                                                                      *
* Purpose  Smoothing using the BiCGStab method.                        *
*          Single precision version, Storage technique 7               *
*                                                                      *
C This subroutine is neary the same like II017 except for it is used   *
C in multigrid as a smoother with calling conventon like XII217.       *
C It has no stopping tolerace parameters but uses only a maximum number*
C of iteration steps for smoothing.                                    *
C The smoothing is only aborted if the residuum is too small.          *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Handles of the arrays describing the matrix in       *
* LCOL     I*4    Storage technique n (single precision)               *
* LLD      I*4                                                         *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Number of smoothing steps                            *
* OMEGA    R*8    Determines preconditioning technique                 *
*                 0 < OMEGA      No Preconditioning                    *
*                 0 = OMEGA      Scaling using diagonal entries        *
*                 0 < OMEGA < 2  SSOR-Preconditioning                  *
*                 2 < OMEGA      Use external subroutine DCG0C         *
* DCG0C    SUBR   EXTERNAL Subroutine DCG0C(DG,NEQ)                    *
*                 Results  DG := C**(-1) * DG  for the precondioning   *
*                 matrix C                                             *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8(NEQ) Smoothed vector                                    *
*                                                                      *
************************************************************************

      SUBROUTINE II227M(LA,LCOL,LLD,DX,DB,NEQ,NIT,OMEGA,DCG0C)
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
      DIMENSION DX(*),DB(*)
      SAVE /ERRCTL/,/CHAR/,/XYPAR/
C
      SUB='II227M'
      IF (ICHECK.GE.997) CALL OTRC('II227M','08/05/03')
C
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE3.NE.2) THEN
       CALL WERR(-170,'II227M')
       GOTO 99999
      ENDIF
C
      IREQ=5*NEQ
      BNOCON=OMEGA.LT.0D0
      IREQ=MAX(IREQ,5)
      CALL ZNEW(IREQ,-1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
      L1=L(LWORK)
      L2=L1+NEQ
      L3=L2+NEQ
      L4=L3+NEQ
      L5=L4+NEQ
C
      KXYPAR(1)=L(LA)
      KXYPAR(2)=L(LCOL)
      KXYPAR(3)=L(LLD)
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
C
       CALL II210(DX,DB,
     *            NEQ,NIT,YLAX37,YIA137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5))
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL II210(DX,DB,
     *            NEQ,NIT,YLAX37,YID137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5))
C
      ELSE
C
       CALL II210(DX,DB,
     *            NEQ,NIT,YLAX37,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5))
C
      ENDIF
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
99999 END

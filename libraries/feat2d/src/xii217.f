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
* XII217                                                               *
*                                                                      *
* Purpose  Smoothing using the BiCGStab method.                        *
*          Double precision version, Storage technique 7               *
*                                                                      *
* This subroutine is neary the same like XII017 except for it is used  *
* in multigrid as a smoother with calling convention like XII017.      *
* It has no stopping tolerace parameters but uses only a maximum       *
* number of iteration steps for smoothing.                             *
* The smoothing is only aborted if the residuum is too small.          *
*                                                                      *
* Subroutines/functions called  LSP1 , LLC1 , LL21 , LLI1, LCL1,       *
*                               II217                                  *
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
*                 CARE: in contrast to the other IIxxx-routines        *
*                 this is working correctly here !!!                   *
* DR,DD    R*8    Workspace vectors of length NEQ                      *
* DD1,DG   R*8    For BNOCON , DG must be replaced by DR               *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Solution vector                                      *
*                                                                      *
************************************************************************

      SUBROUTINE XII217(LA,LCOL,LLD,DX,DB,NEQ,NIT,OMEGA,DCG0C)
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
      IF (ICHECK.GE.997) CALL OTRC('XII217','08/01/03')
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

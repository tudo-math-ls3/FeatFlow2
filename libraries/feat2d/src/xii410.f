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
* XII41n                                                               *
*                                                                      *
* Purpose  Allocate Workspace vector on DWORK                          *
*          Double/double precision version                             *
*          Call XI010M                                                 *
*                                                                      *
* This subroutine is neary the same like XII017 except for it is used  *
* in multigrid as a coarse grid solver. For that purpose it has some   *
* slight modifications in their parameters.                            *
*                                                                      *
* Subroutines/functions called   ZTYPE, ZNEW, ZDISP, II417             *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Numbers of the arrays describing the matrix in       *
* LCOL     I*4    Storage technique n                                  *
* LLD      I*4                                                         *
* LX,LB    I*4    Numbers of the solution and the right hand side      *
* OMEGA    R*8    Determines preconditioning technique                 *
*                 0 < OMEGA      No Preconditioning                    *
*                 0 = OMEGA      Scaling using diagonal entries        *
*                 0 < OMEGA < 2  SSOR-Preconditioning                  *
*                 2 < OMEGA      Use external subroutine DCG0C         *
* DCG0C    SUBR   EXTERNAL Subroutine DCG0C(DG,NEQ)                    *
*                 Results  DG := C**(-1) * DG  for the precondioning   *
*                 matrix C                                             *
* For the description of the remaining input parameters see II010M     *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 -170  X or B are not double precision                *
*                                                                      *
************************************************************************

      SUBROUTINE XII417(LA,LCOL,LLD,DX,DB,NEQ,NIT,ITE,EPS,OMEGA,DCG0C,
     *                  RHO,DEF,RES,DEFINI,RHOASM,BREL)
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
      SUB='II017M'
      IF (ICHECK.GE.997) CALL OTRC('XI017M','08/01/03')
C
      CALL ZTYPE(LA,ITYPE3)
      IF (ITYPE3.NE.1) THEN
       CALL WERR(-170,'II017M')
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
C      IF (BNOCON) L5=L1
C
      KXYPAR(1)=L(LA)
      KXYPAR(2)=L(LCOL)
      KXYPAR(3)=L(LLD)
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
C
       CALL II410(DX,DB,
     *            NEQ,NIT,ITE,EPS,YLAX17,YIA117,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),BREL,RHO,DEF,RES,DEFINI,RHOASM)
C
      ELSE IF (OMEGA.LT.2D0) THEN
C
       CALL II410(DX,DB,
     *            NEQ,NIT,ITE,EPS,YLAX17,YID117,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),BREL,RHO,DEF,RES,DEFINI,RHOASM)
C
      ELSE
C
       CALL II410(DX,DB,
     *            NEQ,NIT,ITE,EPS,YLAX17,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),BREL,RHO,DEF,RES,DEFINI,RHOASM)
C
      ENDIF
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
99999 END


************************************************************************
* XII41N
*
* Solve system using BiCGStab-solver.
*
* This routine allocates workspace vectors and calls the BiCGStab solver
* for a user defined matrix-vector multiplication routine. The caller 
* has to take care about the initialisation of the matrix vector
* multiplication (i.e. initialisation of the parameters in KXYPAR,
* DXYPAR,....) The Y-routine that performs this operation must be
* provided to this routine by the parameter YMVMUL. 
*
* This subroutine is neary the same like XII01N except for it is used  
* in multigrid as a coarse grid solver. For that purpose it has some   
* slight modifications in their parameters.                            
*
* INPUT    TYPE                                                         
* -----    ----                                                         
* DX,DB    R*8    array [1..NEQ] of double;
*                 Solution and the right hand side       
* YMVMUL   SUBR   EXTERNAL SUBROUTINE YMVMUL (DX,DAX,NEQ,A1,A2)
*                 Performs matrix vector multiplication
*                 DAX = A1*DX + A2*DAX
* OMEGA    R*8    Determines preconditioning technique                  
*                 0 < OMEGA      No Preconditioning                     
*                 0 = OMEGA      Scaling using diagonal entries         
*                 0 < OMEGA < 2  SSOR-Preconditioning                   
*                 2 <= OMEGA      Use external subroutine DCG0C          
* DCG0C    SUBR   EXTERNAL Subroutine DCG0C(DG,NEQ)                     
*                 Results  DG := C**(-1) * DG  for the precondioning    
*                 matrix C                                              
* For the description of the remaining input parameters see II010       
*                                                                       
* OUTPUT   TYPE                                                         
* ------   ----                                                         
* ITE      I*4    Number of iterations                                  
* IER      I*4    Error indicator                                       
*                 -170  X or B are not double precision                 
*                                                                       
************************************************************************

      SUBROUTINE XII41N(DX,DB,YMVMUL,NEQ,NIT,ITE,EPS,OMEGA,DCG0C,
     *                  RHO,DABSRS,DRELRS,DEFINI,RHOASM,BREL)

      IMPLICIT NONE
      
      INTEGER NNARR
      PARAMETER (NNARR=299)
      REAL VWORK(1)
      INTEGER KWORK(1)
      INTEGER IER,ICHECK,NWORK,IWORK,IWMAX,L
      DOUBLE PRECISION DWORK,DXYPAR
      INTEGER KXYPAR
      CHARACTER SUB*6,FMT*15,CPARAM*120
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      
C parameters
      
      INTEGER NEQ,NIT,ITE
      DOUBLE PRECISION DX(NEQ),DB(NEQ)
      DOUBLE PRECISION EPS,OMEGA,RHO,DABSRS,DRELRS,DEFINI,RHOASM
      LOGICAL BREL
      EXTERNAL DCG0C,YLAX17,YIA117,YID117,YMVMUL
      
C local variables

      INTEGER IREQ,L1,L2,L3,L4,L5,IER1,LWORK
      LOGICAL BNOCON
      
      SUB='XII017'
      IF (ICHECK.GE.997) CALL OTRC('XII01N','04/14/93')

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
  
      IF (OMEGA.EQ.0D0) THEN

       CALL II010(DX,DB,
     *            NEQ,NIT,ITE,EPS,YMVMUL,YIA117,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),BREL,
     *            RHO,DABSRS,DRELRS,DEFINI,RHOASM)

      ELSE IF (OMEGA.LT.2D0) THEN

       CALL II010(DX,DB,
     *            NEQ,NIT,ITE,EPS,YMVMUL,YID117,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),BREL,
     *            RHO,DABSRS,DRELRS,DEFINI,RHOASM)

      ELSE

       CALL II010(DX,DB,
     *            NEQ,NIT,ITE,EPS,YMVMUL,DCG0C,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),BREL,
     *            RHO,DABSRS,DRELRS,DEFINI,RHOASM)

      ENDIF

      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
99999 END

************************************************************************
* Matrix-vector-multiplication
*
* This performs the following matrix vector multiplication:
* 
*          DD:= C1*(A*DU) + C2*DD
*
* for the vectors  D=(D1,D2,DP)  and  U=(U1,U2,P)  with the 
* given scalar variables C1,C2. The matrix A is the linearised
* nonlinear system matrix for velocity and pressure. Its structure
* is described in "smat2dns.inc"!
*
* In:
*   DU     - array [1..*] of double
*            Input vector
*   DD     - array [1..*] of double
*            Input vector
*   C1,
*   C2     - The constant multipliers
*   MATDAT - array [1..SZN2MI] of integer 
*            TNS2DMatrixParams matrix structure specifying
*            the system matrix that is to multiply with.
* Out:
*   DD     = C1*(A*DU) + C2*DD
************************************************************************

      SUBROUTINE MTMUL2 (DU,DD,C1,C2,MATDAT)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      
      INCLUDE 'smat2dns.inc'

C     parameters

      DOUBLE PRECISION DD(*),DU(*),C1,C2
      INTEGER MATDAT(SZN2MI)
      
C     local variables

      INTEGER KA,KCOLA,KLDA,NEQA,NEQB
      INTEGER KB1,KCOLB1,KLDB1
      INTEGER KB2,KCOLB2,KLDB2
      
C     Calculation of  D1=C2*D1+C1*(A*U1+B1*P)

      NEQA  = MATDAT(ONEQA)
      NEQB  = MATDAT(ONEQB)

      KA    = L(MATDAT(OLA1))
      KCOLA = L(MATDAT(OLCLA1))
      KLDA  = L(MATDAT(OLLDA1))
      KB1    = L(MATDAT(OLB1))
      KCOLB1 = L(MATDAT(OLCLB1))
      KLDB1  = L(MATDAT(OLLDB1))
      
      CALL LAX17 (DWORK(KA),KWORK(KCOLA),KWORK(KLDA),NEQA,
     *            DU(1),DD(1),C1,C2)
      CALL LAX19 (DWORK(KB1),KWORK(KCOLB1),KWORK(KLDB1),NEQA,
     *            DU(1+2*NEQA),DD(1),C1,1D0)

C     Calculation of  D2=C2*D2+C1*(A*U2+B2*P)

      KA    = L(MATDAT(OLA2))
      KCOLA = L(MATDAT(OLCLA2))
      KLDA  = L(MATDAT(OLLDA2))
      KB2    = L(MATDAT(OLB2))
      KCOLB2 = L(MATDAT(OLCLB2))
      KLDB2  = L(MATDAT(OLLDB2))

      CALL LAX17 (DWORK(KA),KWORK(KCOLA),KWORK(KLDA),NEQA,
     *            DU(1+NEQA),DD(1+NEQA),C1,C2)
      CALL LAX19 (DWORK(KB2),KWORK(KCOLB2),KWORK(KLDB2),NEQA,
     *            DU(1+2*NEQA),DD(1+NEQA),C1,1D0)

C     Calculation of  DP=C2*DP+C1*(B1T*U1+B2T*U2)

      CALL LTX19 (DWORK(KB1),KWORK(KCOLB1),KWORK(KLDB1),
     *            NEQA,NEQB,DU(1),DD(1+2*NEQA),C1,C2)
      CALL LTX19 (DWORK(KB2),KWORK(KCOLB2),KWORK(KLDB2),
     *            NEQA,NEQB,DU(1+NEQA),DD(1+2*NEQA),C1,1D0)

      END
      

************************************************************************
* Basic standard-routines for support of matrix structure 7
*
* These routines implement basic matrix/vector and preconditioning
* routines if matrix structure 7 is used as matrix or preconditioner.
* This is only basic support. For extended support of MV/
* preconditioning, special versions have to be programmed.
************************************************************************

************************************************************************
* Matrix-vector multiplication with matrix in matrix structure 7
*
* This standard callback-routine calls matrix-vector multiplication
* routine, based on a matrix in memory structure 7. The caller of the
* solver must initialize the standard-preconditioner block of
* IPARAM/DPARAM the following way:
*
*  IPARAM(OMATI)   = Handle to DA
*  IPARAM(OMATI+1) = Handle to KCOL
*  IPARAM(OMATI+2) = Handle to KLD
************************************************************************

      SUBROUTINE YJAX17(DX,DAX,NEQ,A1,A2, IPARAM,DPARAM)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'      
      INCLUDE 'ssolvers.inc'
      
      INTEGER NEQ
      INTEGER IPARAM(SZSLVI)
      DOUBLE PRECISION DPARAM(SZSLVD)
      DOUBLE PRECISION DX(NEQ), DAX(NEQ), A1, A2
      
C local variables

      INTEGER KA, KCOL, KLD, START
      
      START = OMATI
      KA = L(IPARAM(START))
      KCOL = L(IPARAM(START+1))
      KLD = L(IPARAM(START+2))

      CALL LAX17(DWORK(KA),KWORK(KCOL),
     *           KWORK(KLD),NEQ,DX,DAX,A1,A2)
     
      END

************************************************************************
* Preconditioning with Jacobi (scaling by diagonal), 
* extended calling style.
*
* This standard callback-routine calls the Jacobi preconditioning 
* routine, based on a matrix in memory structure 7. The preconditioning
* is based on the system matrix directly. The caller of the solver
* must initialize the standard-preconditioner block of
* IPARAM/DPARAM the following way:
*
*  IPARAM(OPRECI)   = Handle to DA
*  IPARAM(OPRECI+1) = Handle to KCOL
*  IPARAM(OPRECI+2) = Handle to KLD
************************************************************************

      SUBROUTINE YJA117(DX,NEQ, IPARAM,DPARAM)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'      
      INCLUDE 'ssolvers.inc'
      
      INTEGER NEQ
      INTEGER IPARAM(SZSLVI)
      DOUBLE PRECISION DPARAM(SZSLVD), DX(NEQ)
      
C local variables

      INTEGER KA, KLD
      
      KA = L(IPARAM(OPRECI))
      KLD = L(IPARAM(OPRECI+2))

      CALL IA117(DWORK(KA),KWORK(KLD),DX,NEQ)
     
      END

************************************************************************
* Preconditioning with SSOR, extended calling style.
*
* This standard callback-routine calls the SSOR preconditioning routine,
* based on a matrix in memory structure 7. The preconditioning is
* based on the system matrix directly. The caller of the solver
* must initialize the standard-preconditioner block of
* IPARAM/DPARAM the following way:
*
*  IPARAM(OPRECI)   = Handle to DA
*  IPARAM(OPRECI+1) = Handle to KCOL
*  IPARAM(OPRECI+2) = Handle to KLD
*
*  DPARAM(OPRECD)  = preconditioning parameter
************************************************************************

      SUBROUTINE YJD117(DX,NEQ, IPARAM,DPARAM)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'      
      INCLUDE 'ssolvers.inc'
      
      INTEGER NEQ
      INTEGER IPARAM(SZSLVI)
      DOUBLE PRECISION DPARAM(SZSLVD), DX(NEQ)
      
C local variables

      INTEGER KA, KCOL, KLD
      
      KA = L(IPARAM(OPRECI))
      KCOL = L(IPARAM(OPRECI+1))
      KLD = L(IPARAM(OPRECI+2))

      CALL ID117(DWORK(KA),KWORK(KCOL),
     *           KWORK(KLD),DX,NEQ,DPARAM(OPRECD))
     
      END

************************************************************************
* Preconditioning with ILU0, extended calling style.
*
* This standard callback-routine calls the ILU0 preconditioning routine,
* based on a matrix in memory structure 7. The preconditioning is
* based on a ILU0-decomposition of a matrix, given in matrix 
* structure 7. The caller of the solver must initialize the standard-
* preconditioner block of IPARAM/DPARAM the following way:
*
*  IPARAM(OPRECI)   = Handle to DA of the ILU0-decomposition
*  IPARAM(OPRECI+1) = Handle to KCOL of the ILU0-decomposition
*  IPARAM(OPRECI+2) = Handle to KLD of the ILU0-decomposition
*
*  DPARAM(OPRECD)   = preconditioning parameter
************************************************************************

      SUBROUTINE YJF117(DX,NEQ, IPARAM,DPARAM)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'      
      INCLUDE 'ssolvers.inc'
      
      INTEGER NEQ
      INTEGER IPARAM(SZSLVI)
      DOUBLE PRECISION DPARAM(SZSLVD), DX(NEQ)
      
C local variables

      INTEGER KA, KCOL, KLD
      
      KA = L(IPARAM(OPRECI))
      KCOL = L(IPARAM(OPRECI+1))
      KLD = L(IPARAM(OPRECI+2))

      CALL IF117(DWORK(KA),KWORK(KCOL),
     *           KWORK(KLD),DX,NEQ,DPARAM(OPRECD))
     
      END
      

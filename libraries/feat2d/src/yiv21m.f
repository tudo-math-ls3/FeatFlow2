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
* YIV21M                                                               *
*                                                                      *
* Purpose: Preconditioning with (M)ILU(s) using SPLIB/LUSOLT           *
*  - Double precision vector                                           *
*  - Double precision (M)ILU(s) matrix, stored in KWORK                *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DG       R*8(NEQ)   Vector to be preconditioned                      *
* NEQ      I*4        Number of equations/Entries in DG                *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DG       R*8(NEQ)   preconditioned vector                            *
*                                                                      *
* External information provided by COMMON blocks:                      *
*                                                                      *
*  KXYPAR(5) = KLU   - Pointer to (M)ILU(s) data structure in KWORK.   *
*                      Let LU denote the Handle of the (M)ILU(s) matrix*
*                      in KWORK that results by a call to IVD4x7.      *
*                      Then KLU = L(LU).                               *
*  KXYPAR(6) = LLU   - Pointer in the (M)ILU(s) data structure to the  *
*                      LU matrix. Result of call to IVD4x7.            *
*  KXYPAR(7) = LLUJ  - Pointer in the (M)ILU(s) data structure to JLU. *
*                      Result of call to IVD4x7.                       *
*  KXYPAR(8) = LUPTR - Pointer in the (M)ILU(s) data structure to the  *
*                      array describing the starting offsets of U in   *
*                      each row. Result of call to IVD4x7.             *
*                                                                      *
************************************************************************

      SUBROUTINE YIV21M(DG,NEQ)      

      IMPLICIT NONE
      
      INTEGER NNARR
      PARAMETER (NNARR=299)

C *** Standard COMMON blocks

      INTEGER NWORK, IWORK, IWMAX, L, KWORK, IER, ICHECK, KXYPAR
      DOUBLE PRECISION DWORK, DXYPAR
      REAL VWORK

      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      DIMENSION VWORK(1),KWORK(1),DWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))

      SAVE /ERRCTL/,/XYPAR/

C Local variables

      INTEGER NEQ, KLU, LLU, LLUJ, LUPTR
      DOUBLE PRECISION DG(*)

      IF (ICHECK.EQ.999) CALL OTRC('YIV21M','08/10/04')

      KLU = KXYPAR(5)
      LLU = KXYPAR(6)
      LLUJ = KXYPAR(7)
      LUPTR = KXYPAR(8)
      
C Solve the system with SPLIB's LUSOLT command
      
      CALL LUSOLT (NEQ,DG,KWORK(KLU+LLU-1),
     *             KWORK(KLU+LLUJ-1),KWORK(KLU+LUPTR-1))

99999 END

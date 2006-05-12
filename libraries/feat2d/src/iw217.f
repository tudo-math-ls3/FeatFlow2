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
* IW217                                                                *
*                                                                      *
* Purpose: Smoothing with standard Richardson(Omega) smoother          *
*          x_(n+1) = x_n + Omega*(b-Ax)                                *
*          Double precision version                                    *
*          Matrix structure 7                                          *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8(NA)    Reference to system matrix                       *
* KCOL     I*4(NA)    Column description of matrix                     *
* KLD      I*4(NEQ+1) Row description of matrix                        *
* NEQ      I*4        Number of equations                              *
* DX       R*8(NEQ)   Starting vector                                  *
* DB       R*8(NEQ)   Right hand side                                  *
* DD       R*8(NEQ)   Auxiliary vector of length NEQ                   *
* NIT      I*4        Number of iterations                             *
* OMEGA    R*8        Relaxation vector                                *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8(NEQ)   Vector after NIT smoothing steps with R.-smoother*
*                                                                      *
************************************************************************

      SUBROUTINE IW217 (DA,KCOL,KLD,DX,DB,DD,NEQ,NIT,OMEGA)

      IMPLICIT NONE

C parameters

      DIMENSION DA(*),KCOL(*),KLD(*)
      
      INTEGER KCOL, KLD, NIT, NEQ
      DOUBLE PRECISION DA, DX(NEQ), DD(NEQ), DB(NEQ), OMEGA
      
C local variables

      INTEGER I

C perform NIT smoothing steps to DX
C
C           dx = dx + Omega*(b-Ax)

      DO I = 1,NIT
C       *** dd = b-Ax = -DA*dx+db
        CALL LCP1 (DB,DD,NEQ)
        CALL LAX17 (DA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)
C       *** dx = dx + Omega*dd
        CALL LLC1 (DD,DX,NEQ,OMEGA,1D0)
      END DO

      END


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
* IVD01M                                                               *
*                                                                      *
* Purpose: Calculate ILU(s) matrix.                                    *
*          Resulting matrix will use SPLIB matrix format.              *
*          Double precision version                                    *
*                                                                      *
* The structure of the matrix will be determined automatically. The    *
* memory structure that is returned in the parameters is the origional *
* MSR structure of the SPLIB library.                                  *
*                                                                      *
* For the computation of the structure the routine allocates the whole *
* memory available by the pseudodynamic memory management. After the   *
* computation the unused memory is freed again.                        *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8(NA)    Reference to system matrix                       *
* KCOL     I*4(NA)    Column description of matrix                     *
* KLD      I*4(NEQ+1) Row description of matrix                        *
* NEQ      I*4        Number of equations                              *
* NS       I*4        Fill-in parameters for ILU(s); 0=standard=ILU(0) *
* OMEGA    R*8        MILU relaxation parameter. If <> 0, compute      *
*                     MILU(Omega,s) matrix instead of ILU(s)           *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LU       I*8        Handle to ILU(s) matrix structure (MSR-format)   *
* LLU      I*8        MSR-internal parameter; Pointer to LU matrix     *
* LLUJ     I*8        MSR internal parameter; Pointer to JLU           *
* LUPTR    I*8        MSR internal parameter; Pointer to array of      *
*                     starting columns of U matrix                     *
* MNEED    I*8        Number of integer values allocated temporarily   *
*                     for computation of the ILU(s) matrix in SPLIB.   *
*                                                                      *
* When calling a Routine that is working with the MSR structure,       *
* the internal parameters LLU, LLUJ, LUPTR must be provided together   *
* with the handle to the matrix. As these are internal parameters,     *
* they mustn't be changed!                                             *
*                                                                      *
************************************************************************

      SUBROUTINE IVD01M (DA,KCOL,KLD,NEQ,NS,OMEGA,LU,LLU,LLUJ,
     *                   LUPTR,MNEED)

      IMPLICIT NONE

C Parameters

      DOUBLE PRECISION DA(*)
      DOUBLE PRECISION OMEGA
      INTEGER KCOL(*), KLD(*)
      INTEGER NEQ, NS, LU, LLU, LLUJ, LUPTR, MNEED

C Common blocks

      INTEGER IER,ICHECK,NWORK,IWORK,IWMAX,L
      INTEGER M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      INTEGER NNARR,KWORK
      REAL VWORK
      DOUBLE PRECISION DWORK
      
      PARAMETER (NNARR=299)
      
      DIMENSION VWORK(1),KWORK(1)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /ERRCTL/,/OUTPUT/

C internal variables

      INTEGER MNEED2      

C Temporarily reserve the whole memory available. Use the INTEGER type,
C because the ILU decomposition expects an integer array which is used
C for double dumbers, too.

      MNEED2 = 0
      CALL ZNEW (MNEED2,-3,LU,'LU    ')

C ZNEW returns in MNEED2, how much memory is reserved

C Calculate ILU 

      IF (IER.EQ.0) THEN
        CALL ILUS (NEQ,DA,KCOL,KLD,NS,OMEGA,LLU,LLUJ,LUPTR,
     *             KWORK(L(LU)),MNEED2,IER,MNEED)
      END IF
      IF (IER.NE.0) THEN 
        CALL ZDISP (0,LU,'LU    ')
        IER = 1
        WRITE (MTERM,'(A)') ' Not enough memory to compute ILU !'
        GOTO 99999
      END IF

C Dispose unused memory; use MNEED therefore, as it returns the actual size 
C of the matrix

      CALL ZDISP (MNEED,LU,'LU    ')

99999 END

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
* IV21M                                                                *
*                                                                      *
* Purpose: Smoothing with ILU(s), SPLIB's MSR matrix format            *
*          Precomputed matrix                                          *
*          x_(n+1) = x_n + Omega*M*(b-Ax)                              *
*          Double precision version                                    *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8(NA)    Reference to system matrix                       *
* KCOL     I*4(NA)    Column description of matrix                     *
* KLD      I*4(NEQ+1) Row description of matrix                        *
* KLU      I*8        Handle to ILU(s) matrix structure (MSR-format)   *
* KLLU     I*8        MSR-internal parameter; Pointer to LU matrix     *
* KLLUJ    I*8        MSR internal parameter; Pointer to JLU           *
* KLUP     I*8        MSR internal parameter; Pointer to array of      *
*                     starting columns of U matrix                     *
* DX       R*8(NEQ)   Starting vector                                  *
* DB       R*8(NEQ)   Right hand side                                  *
* DD       R*8(NEQ)   Auxiliary vector of length NEQ                   *
* NEQ      I*4        Number of equations                              *
* NIT      I*4        Number of iterations                             *
* RLXSM    I*4        Relaxation parameter for smoothing               *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8(NEQ)   Vector after NIT smoothing steps                 *
*                                                                      *
************************************************************************

      SUBROUTINE IV21M (DA,KCOL,KLD,KLU,KLLU,KLUJ,KLUP,DX,DB,DD,
     *                   NEQ,NIT,RLXSM)

      IMPLICIT NONE

C Parameters

      INTEGER KCOL(*),KLD(*),KLU(*),KLLU,KLUJ,KLUP,NEQ,NIT
      DOUBLE PRECISION DA(*),DX(NEQ),DB(NEQ),DD(NEQ),RLXSM
      
C Common blocks

      INTEGER IER,ICHECK,NWORK,IWORK,IWMAX,L
      INTEGER M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      INTEGER NNARR,KWORK
      REAL VWORK
      DOUBLE PRECISION DWORK
      
      PARAMETER (NNARR=299)
      
      DIMENSION VWORK(1),KWORK(1)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /ERRCTL/,/OUTPUT/

C local variables

      INTEGER I

C perform NIT smoothing steps to DX.
C           dx = dx + Omega*M*(b-Ax)

      DO I = 1,NIT
      
C       *** dd = b-Ax = -DA*dx+db

        CALL LCP1 (DB,DD,NEQ)
        CALL LAX17 (DA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)

C Use ILU preconditioning of defect vector for smoothing. This is done
C by solving the linear systems Uy=dd, Ldb=y with the SPLIB-Routine LUSOLT:

         CALL LUSOLT (NEQ,DD,KLU(KLLU),KLU(KLUJ),KLU(KLUP))
         
C Use the relaxation parameter before adding the defect vector
C to the solution vector: dx = dx + omega*dd

        CALL LLC1 (DD,DX,NEQ,RLXSM,1D0)

      END DO

      END

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
* IV21MD                                                               *
*                                                                      *
* Purpose: Smoothing with ILU(s), SPLIB's MSR matrix format            *
*          Direct smoothing, smoothing matrix is computed before and   *
*          destroyed after the smoothing automatically.                *
*          x_(n+1) = x_n + Omega*M*(b-Ax)                              *
*          Double precision version                                    *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8(NA)    Reference to system matrix                       *
* KCOL     I*4(NA)    Column description of matrix                     *
* KLD      I*4(NEQ+1) Row description of matrix                        *
* DX       R*8(NEQ)   Starting vector                                  *
* DB       R*8(NEQ)   Right hand side                                  *
* DD       R*8(NEQ)   Auxiliary vector of length NEQ                   *
* NS       I*4        Fill-in parameters for ILU(s); 0=standard=ILU(0) *
* OMEGA    R*8        MILU relaxation parameter. If <> 0, use          *
*                     MILU(Omega,s) matrix instead of ILU(s)           *
* NEQ      I*4        Number of equations                              *
* NIT      I*4        Number of iterations                             *
* RLXSM    I*4        Relaxation parameter for smoothing               *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8(NEQ)   Vector after NIT smoothing steps                 *
*                                                                      *
************************************************************************

      SUBROUTINE IV21MD (DA,KCOL,KLD,NS,OMEGA,DX,DB,DD,
     *                   NEQ,NIT,RLXSM)

C Parameters

      DOUBLE PRECISION DA(*),DX(NEQ),DB(NEQ),DD(NEQ),RLXSM,OMEGA
      INTEGER KCOL(*),KLD(*),NEQ,NIT,NS

C Common blocks

      INTEGER IER,ICHECK,NWORK,IWORK,IWMAX,L
      INTEGER M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      INTEGER NNARR,KWORK
      REAL VWORK
      DOUBLE PRECISION DWORK
      
      PARAMETER (NNARR=299)
      
      DIMENSION VWORK(1),KWORK(1)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /ERRCTL/,/OUTPUT/

C local variables

      INTEGER LU

C Calculate ILU(s) matrix

C      print *,'Calculating ILU for ',NEQ,' equations...'
      CALL IVD01M (DA,KCOL,KLD,NEQ,NS,OMEGA,LU,LLU,LLUJ,
     *             LUPTR,MNEED)

C Perform smoothing

C      print *,'Smoothing...'
      CALL IV21M (DA,KCOL,KLD,KWORK(L(LU)),LLU,LLUJ,LUPTR,DX,DB,DD,
     *            NEQ,NIT,RLXSM)


C finally destroy ILU matrix again

C      print *,'Destroying ILU...'
      CALL ZDISP (0,LU,'ILkTMP')

      END

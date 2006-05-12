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
* IVD017                                                               *
*                                                                      *
* Purpose: Calculate ILU(s) matrix.                                    *
*          Resulting matrix will be saved in memory structure 7.       *
*          Double precision version                                    *
*                                                                      *
* The structure of the matrix will be determined automatically. The    *
* subroutine returns handles to DA/KCOL/KLD in the variables           *
* LM/LMCOL/LMLD.                                                       *
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
* LM       I*8        Handle to array DA of matrix                     *
* LMCOL    I*8        Handle to KCOL of matrix                         *
* LMLD     I*8        Handle to KLD of matrix                          *
* MNEED    I*8        Number of integer values allocated temporarily   *
*                     for computation of the ILU(s) matrix in SPLIB.   *
*                                                                      *
************************************************************************

      SUBROUTINE IVD017 (DA,KCOL,KLD,NEQ,NS,OMEGA,LM,LMCOL,LMLD,
     *                   MNEED)
     
      IMPLICIT NONE
     
C Parameters

      DOUBLE PRECISION DA(*)
      DOUBLE PRECISION OMEGA
      INTEGER KCOL(*), KLD(*)
      INTEGER NEQ, NS, LM, LMCOL, LMLD, MNEED

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

      INTEGER MNEED2, LTMP, LU, JLU, ILUP, JLU1, LU1I, MSRNA
      INTEGER MNA
      EXTERNAL MSRNA

C Temporarily reserve the whole memory available. Use the INTEGER type,
C because the ILU decomposition expects an integer array which is used
C for double dumbers, too.

      MNEED2 = 0
      CALL ZNEW (MNEED2,-3,LTMP,'LTMP  ')

C ZNEW returns in MNEED2, how much memory is reserved

C Calculate ILU 

      IF (IER.EQ.0) THEN
        CALL ILUS (NEQ,DA,KCOL,KLD,NS,OMEGA,LU,JLU,ILUP,
     *             KWORK(L(LTMP)),MNEED2,IER,MNEED)
      END IF
      IF (IER.NE.0) THEN 
        CALL ZDISP (0,LTMP,'LTMP  ')
        IER = 1
        WRITE (MTERM,'(A)') ' Not enough memory to compute ILU !'
        GOTO 99999
      END IF

C Dispose unused memory; use MNEED therefore, as it returns the actual size 
C of the matrix

      CALL ZDISP (MNEED,LTMP,'LTMP  ')

C Info: The array L(LTMP) contains also the content of the matrix, which is
C saved as DOUBLE's in this integer verctor!

C Determine the size of the matrix; therefore compute the starting addresses
C of the structores in the L(LTMP) array

      JLU1 = L(LTMP)+JLU-1
      LU1I = L(LTMP)+LU-1

C Determine the size of the structure of the storage technique 7; therefore
C compute the number of entries in the ILU(s) matrix

      MNA = MSRNA (KWORK(JLU1),NEQ)

C Reserve memory for matrix in storage technique 7

      CALL ZNEW (NEQ+1,3,LMLD,'LMLD  ')
      CALL ZNEW (MNA,3,LMCOL,'LMCOL ')
      CALL ZNEW (MNA,1,LM,'LM    ')
      
C Copy MSR structure to just reserved memory. There is only one problem:
C SPLIB saves DOUBLE values in the INTEGER array!

C Therefore use the usual trick of FORTRAN77: Provide the address of
C the integer array in the parameter list where a pointer to a DOUBLE
C array is expected...

C There might be a compiler warning here - this is OK !!!

      CALL MSR2T7 (KWORK(LU1I),KWORK(JLU1),
     *             DWORK(L(LM)),KWORK(L(LMCOL)),KWORK(L(LMLD)),
     *             NEQ,1) 

C Deallocate unused memory, finish.

      CALL ZDISP (0,LTMP,'LTMP  ')

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
* IV217                                                                *
*                                                                      *
* Purpose: Smoothing with ILU(s)                                       *
*          Precomputed matrix, storage structure 7                     *
*          x_(n+1) = x_n + Omega*M*(b-Ax)                              *
*          Double precision version                                    *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8(NA)    Reference to system matrix                       *
* KCOL     I*4(NA)    Column description of matrix                     *
* KLD      I*4(NEQ+1) Row description of matrix                        *
* DM       R*8(NA)    Reference to ILU matrix                          *
* KMCOL    I*4(NA)    Column description of ILU matrix                 *
* KMLD     I*4(NEQ+1) Row description of ILU matrix                    *
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

      SUBROUTINE IV217 (DA,DM,KCOL,KLD,KMCOL,KMLD,DX,DB,DD,
     *                  NEQ,NIT,RLXSM)
     
      IMPLICIT NONE
     
C Parameters

      INTEGER KCOL(*),KLD(*),KMCOL(*),KMLD(*),NEQ,NIT
      DOUBLE PRECISION DA(*),DX(NEQ),DB(NEQ),DD(NEQ),RLXSM
      DOUBLE PRECISION DM(*)
      
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

C perform NIT smoothing steps to dx
C           dx = dx + M*(b-Ax)

      DO I = 1,NIT

C       *** dd = b-Ax = -DA*dx+db

        CALL LCP1 (DB,DD,NEQ)
        CALL LAX17 (DA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)
        
C Use ILU preconditioning of the defect vector for smoothing;
C therefore solve the linear systems Uy=dd, Ldb=y with the help
C of the subroutine LUSLV7:

        CALL LUSLV7 (NEQ,DD,DM,KMCOL,KMLD)

C Now multiply the preconditioned defect vector by the relaxation 
C parameter RLXSM and add it to the solution vector:
C   dx = dx + omega*dd     (omega=RLXSM here)

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
* LUSLV7                                                               *
*                                                                      *
* Purpose: Perform a forward then backward solve for a CSR matrix      *
*          containing a unit lower triangular and an upper triangular  *
*          matrix, both stored in a single CSR data structure.         *
*          Double precision version                                    *
*                                                                      *
* This routine is origionally a routine of the SPLIB library. It has   *
* been slightly modified do work with DOUBLE PRECISION arrays in       *
* storage technique 7, expecially as these matrices don't have an      *
* inverted diagonal.                                                   *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* N        I*4        Number of equations NEQ                          *
* LU       R*8(NA)    Pointer to matrix entries DA (not a handle!!!)   *
* KCOL     I*4(NA)    Column description of matrix                     *
* KLD      I*4(NEQ+1) Row description of matrix                        *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* X        R*8(NEQ)   Solution vector                                  *
*                                                                      *
************************************************************************

      subroutine luslv7 (n, x, lu, kcol, kld)

        implicit none
        integer kcol(*),kld(*),n,i,k,icol
        double precision x(n), lu(*)
        double precision x_i

        if (n .le. 0)  return

c       ------------
c       Solve Lx = x : x_i = 1/1 * (x_i - a_i1 x1 - a_i2 x_2 - ...)
c       ------------
        do i = 1, n
        
C Take the right hand side
           x_i =  x(i)

C Subtract the sum of the left hand side
           do k = kld(i)+1,kld(i+1)-1
             icol = kcol (k)
             if (icol.ge.i) goto 10
             x_i = x_i - lu(k) * x(icol)
           end do
           
C Store the value
10         x(i) = x_i

        end do


c       ------------
c       Solve Ux = x: x_i = 1/aii * (x_i - a_in x_n - a_in-1 x_n-1 - ...)
c       ------------
        do i = n, 1, -1
        
C Take the right hand side
           x_i =  x(i)
           
           do k=kld(i+1)-1,kld(i)+1,-1
             icol = kcol (k)
             if (icol.le.i) goto 20
             x_i = x_i - lu(k) * x(icol)
           end do
C Divide by the diagonal and store the value
20         x(i) = x_i / lu(kld(i))

        end do

        return
      end


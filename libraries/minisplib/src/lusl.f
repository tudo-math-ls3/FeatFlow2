c========================================================================
      subroutine lusl (n, y, x, lu, jlu, uptr, pcmeth, a, colind, 
     *                 rwptr ) 

c**********************************************************************
c                   SPLIB Copyright (C) 1995                         **
c                     Indiana University                             **
c**********************************************************************

c-----------------------------------------------------------------------*
c                     *** LUSL: Solve Mx = y ***			*
c									*
c   Preconditioner application for SPLIB.  This routine computes the	*
c   action of the preconditioner on a vector; usually this is a lower	*
c   and upper triangular solve, explaining the name.  In particular,	*
c   (lu, jlu) is usually the MSR data structure for L and U combined	*
c   with the integer vector uptr giving the pointers to the start	*
c   position of the upper triangular part of each row.  However, for	*
c   SSOR preconditioning lu is the scalar relaxation parameter omega	*
c   and uptr is a vector of pointers to the diagonal entries in the	*
c   data structure for A (not LU). In this case the names lu and jlur	*
c   are somewhat deceptive, but have been kept since ILU, MILU, ILUT,	*
c   and ECIMGS all use the MSR data structure for LU.			*
c									*
c   Although x and y can be different, the methods implemented so far	*
c   allow overwriting, so there is an inefficiency in copying y over	*
c   to x unnecessarily.							*
c									*
c-----------------------------------------------------------------------*
c
c  Arguments:
c  ==========
c
c  x       : double precision vector of length n (output).
c  y       : double precision vector of length n (input).
c
c  n       : Order of matrix A, and system Mx = y where M is the
c            preconditioning matrix.
c  A       : double precision array for nonzero entries of matrix A; of
c            length nnz = rwptr(n+1)-1.
c  colind  : Integer array of column indices for entries in A; of
c            length nnz = rwptr(n+1)-1.
c  rwptr   : Integer array of pointers to beginning position of rows
c            in arrays A and colind.  Must be of length n+1.
c  uptr, 
c  lu,jlu  : Data structure for the preconditioner M, with different meanings
c            for different preconditioners.
c
c            ILU(s), MILU(s), ILUt, or ECIMGS --->
c                 (lu,jlu) gives the combined L/U factors for the preconditioner
c                 in a modified sparse row (MSR) format.  L is unit lower 
c                 triangular and U is upper triangular.  The diagonal of U is
c                 stored in lu(1:n) inverted, to replace divisions in the algorithm
c                 with multiplications.  Each i-th row of the combined L/U matrix in
c                 the (lu,jlu) data structure contains the i-th row of L (excluding
c                 the diagonal entry = 1), followed by the i-th row of u.
c                 uptr is a integer vector of length n containing pointers to the
c                 start position of each row of U in (lu, jlu).
c
c            SSOR(omega)   --->
c                 lu is the scalar acceleration parameter omega.
c                 jlu is not used, and uptr is a vector of pointers to the
c                 diagonal entries in the CSR data structure for A.
c
c            TRID(s)   --->
c                 This implements block diagonal preconditioning, where
c                 each diagonal block is from the corresponding tridiagonal part
c                 of the matrix A.  The size of the blocks is given implicitly in
c                 the integer iskip, stored in jlu(n).  Factors are set up for 
c                 contiguous blocks of iskip rows and so by using iskip = n the 
c                 factor can be for the entire matrix's tridiagonal part.
c                 Using iskip = 1 is an inefficient way of performing diagonal
c                 preconditioning, while 1 < iskip < n is for a block tridiagonal
c                 preconditioner.  Gaussian elimination with pivoting is used for 
c                 factoring each block, so it is necessary to store the pivot
c                 sequence in jlu(1:n-1).  lu has the four vectors giving the nonzeros
c                 of the tridiagonal preconditioner.  Vector uptr is not used.
c
c            ILU0   --->
c                 This implements the version of ILU(0) that does not change any
c                 off-diagonal entries of A, and so only needs to store the n-vector
c                 of diagonal pivots.  lu(1:n) stores those pivots, and jlu(1:n) 
c                 contains pointers to the diagonal entry in each row of A.
c
c  pcmeth  : preconditioning method to use.
c		pcmeth = 0  ---> no preconditioner, i.e. m = identity matrix.
c		pcmeth = 1  ---> ILU(c).
c		pcmeth = 2  ---> MILU(c).
c		pcmeth = 3  ---> ILUT(c).
c		pcmeth = 4  ---> SSOR(omega).
c		pcmeth = 5  ---> TRID.
c		pcmeth = 6  ---> ILU0.
c		pcmeth = 7  ---> ECIMGS.
c
c
c  Routines
c  ========
c
c  lusolt  : performs a forward and a backward solve for matrices L and U
c            stored in combined MSR format, with L unit upper triangular.
c            the diagonal entries of U are assumed to be inverted already.
c  lussor  : Performs forward and backward solve for SSOR preconditioner.
c  lusol0  : Performs forward and backward solve for ILU0 preconditioner.
c  tridsl  : Applies the tridiagonal block preconditioner.
c  dcopy   : BLAS1 routine to copy a vector.
c-------------------------------------------------------------------------
c
c     Randall Bramley
c     Department of Computer Science
c     Indiana University, Bloomington
c     bramley@cs.indiana.edu
c     Wed Jul 12 10:03:18 EST 1995
c
c-------------------------------------------------------------------------

      implicit none

c     Declarations of arguments
c     =========================
      double precision y(*),x(*),lu(*)
      double precision a(*)
      integer n
      integer jlu(*),uptr(*),pcmeth
      integer colind(*),rwptr(n+1)

c     Local variables
c     ===============
      integer k, length, ierr

c     Instrumentation variables for timing.
c     =====================================
      include '../include/heads.inc'
      double precision mytime, tstart
      external mytime

c========================================================================
c      Beginning of Executable Statements
c========================================================================

c      tstart = mytime()

c     ----------------
c     Copy y over to x
c     ----------------
      call dcopy(n, y, 1, x, 1)

c     ------------------------------------------
c     Unknown pcmeth (should not be encountered)
c     ------------------------------------------
      if (pcmeth.lt.0 .or. pcmeth .gt. 7) then
         write(*,*) 'Unknown method in lusl'
         write(*,*) 'pcmeth = ', pcmeth
         stop
      end if

c     ------------------
c     No preconditioning
c     ------------------
      if (pcmeth.eq.0) then
         return
      endif

c     --------------------------
c     ILU, MILU, ILUT, or ECIMGS
c     --------------------------
      if ( ((pcmeth.ge.1) .and. (pcmeth.le.3)) .or. pcmeth .eq. 7)
     +       call lusolt (n, x, lu, jlu, uptr) 

c     --------------------
c     SSOR preconditioning
c     --------------------
      if (pcmeth.eq.4)
     +       call lussor(n, x, lu, jlu, a,colind,rwptr ) 

c     --------------------
c     ILU0 preconditioning
c     --------------------
      if (pcmeth.eq.6)
     +       call lusol0 (n, x, lu, jlu, a,colind,rwptr ) 

c     --------------------
c     TRID preconditioning
c     --------------------
      if (pcmeth.eq.5) then
C         do 100 k = 1, n, jlu(n)
C            length = min(jlu(n),n-k+1)
C            call tridsl('N',length,lu(1+k),
C     &                             lu(k+n),
C     &                             lu(k+2*n),
C     &                             lu(k+3*n),
C     &                             jlu(k),x(k),ierr)
C 100     continue
      end if

c      times(5) = times(5) + mytime() - tstart

      return

c===================== End of LUSL =====================================

      end
c=======================================================================
      subroutine lusolt (n, x, lu, jlu, uptr)
c-----------------------------------------------------------------------
c
c Perform a forward then backward solve for a MSR matrix containing
c a unit lower triangular and an upper triangular matrix with inverted
c diagonal, both stored in a single MSR data structure.  The vector of
c pointers uptr() contains the row pointers to the beginning of U in
c each row.
c
c-----------------------------------------------------------------------
        implicit none
        integer jlu(*),uptr(*),n,i,k
        double precision x(n), lu(*)
        double precision x_i

        if (n .le. 0)  return

c       ------------
c       Solve Lx = x
c       ------------
        do i = 1, n-1
           x_i =  x(i)

c*$* prefetch_ref=lu(jlu(i+1)),stride=1,kind=rd
c*$* prefetch_ref=lu(jlu(i+1)+16),stride=1,kind=rd
c*$* prefetch_ref=jlu(jlu(i+1)+32),stride=1,kind=rd

           do k = jlu(i), uptr(i)-1
              x_i = x_i - lu(k)* x(jlu(k))
           end do
           x(i) = x_i
        end do

        x_i =  x(n)
        do k = jlu(n), uptr(n)-1
           x_i = x_i - lu(k)* x(jlu(k))
        end do
        x(n) = x_i


c       ------------
c       Solve Ux = x
c       ------------
        do i = n, 1, -1
           x_i =  x(i)

c*$* prefetch_ref=lu(jlu(i)-1),stride=1,kind=rd
c*$* prefetch_ref=lu(jlu(i)-17),stride=1,kind=rd
c*$* prefetch_ref=jlu(jlu(i)-33),stride=1,kind=rd

           do k = jlu(i+1)-1, uptr(i), -1
              x_i = x_i - lu(k)*x(jlu(k))
           end do
           x(i) = lu(i)*x_i
        end do

        return
      end
c-------------------End of lussor---------------------------------------

c----------------------------------------------------------------------- 
      subroutine lusol0(n, x, pivots, diag, a,colind,rwptr) 
c----------------------------------------------------------------------- 
c
c  Solve the system LUx = x, where L and U were created using
c  ILU(-1) preconditioner, the space saver version of ILU(0).
c  The implicit splitting is
c          LU = (D + LA) (I + (D)^-1 UA)
c
c----------------------------------------------------------------------- 
      implicit none
      integer n, colind(*), rwptr(n+1), diag(n)
      double precision x(n), pivots(n)
      double precision a(*), sum
      integer i, j

      do i = 1, n
         sum = 0.0d0
         do j = rwptr(i), diag(i)-1
            sum = sum + a(j)*x(colind(j))
         end do
         x(i) = pivots(i)*(x(i) - sum)
      end do

      do i = n, 1, -1
         sum =  0.0d0
         do j = diag(i)+1, rwptr(i+1)-1
            sum = sum + a(j)*x(colind(j))
         end do
         x(i) = x(i) - pivots(i)*sum
      end do

      return
      end

c-------------------End of lusol0--------------------------------------- 
c-----------------------------------------------------------------------
      subroutine lussor(n, x, omega, diag, a,colind,rwptr ) 
c-----------------------------------------------------------------------
c
c  Solve the two triangular systems given for the SSOR preconditioner
c  with acceleration parameter omega.
c
c  System solved is 
c      Mx = y, with M = (D/omega + L) * ( I + omega * D^-1 * U).
c
c  diag : vector of pointers to the diagonal entry in each row k of A.
c  Warning: the rows of A must be sorted so that within each row,
c           the subdiagonal entries come first, then the diagonal entry,
c           and then the superdiagonal entries.  Furthermore, the diagonal
c           of A must be nonzero.
c
c-----------------------------------------------------------------------

      implicit none
      integer n, colind(*), rwptr(n+1), diag(n)
      double precision x(n)
      double precision a(*), omega, sum
      integer j, k

      do k = 1, n
         sum = 0.0d0
         do j = rwptr(k), diag(k)-1
            sum = sum + a(j)*x(colind(j))
         end do
         x(k) = (x(k)-sum) * omega / a(diag(k))
      end do

      do k = n, 1, -1
         sum = 0.0d0
         do j = diag(k)+1, rwptr(k+1)-1
            sum = sum + a(j)*x(colind(j))
         end do
         x(k) = x(k) - sum * omega / a(diag(k))
      end do

      return
      end
c-------------------End of lussor---------------------------------------

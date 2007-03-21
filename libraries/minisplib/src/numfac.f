c======================================================================== 
      subroutine  numfac(n, colind,rwptr,  jlu, uptr, a, lu,
     *				relax, ierr, colptrs, milu)

c**********************************************************************
c                   SPLIB Copyright (C) 1995                         **
c                     Indiana University                             **
c**********************************************************************

c======================================================================== 
c
c Numerical factorization with given sparsity pattern identified in
c the MSR data structure  jlu.
c
c Because of the MSR storage format, the method used is a deferred update
c version.  At the k-th step, the updates from the previous rows m,
c m < k and m corresponding to a nonzero in row k of L/U, are applied
c to row k.  A relaxation parameter is used.  When it is 0.0, ILU
c factorization is performed.  Otherwise the parameter is multiplied
c times the discarded fillin and added to the diagonal entry.  So when 
c the parameter is equal to 1.0, classical MILU factorization is performed.
c
c This routine implicitly assumes that the sparsity pattern specified
c by the MSR data structure contains that of A.  Note the comments in
c the code if you want to remove this restriction.
c
c Matrix A is input in CSR format, and the output MSR array contains
c both L and U, with L assumed to be unit lower triangular.  The diagonal
c entries of U are stored in inverted form to avoid divisions when
c applying the preconditioner.
c
c This version differs from earlier SPLIB versions in performing the
c computations in-place, without using a software accumulator register
c row(:).  The reasoning is that for very large systems, the memory
c access patterns when using such a register will cause many cache
c misses.  The technique (suggested by Xiaoge Wang) is to use an integer*4
c vector colptrs which simultaneously identifies which columns correspond
c to allowed fill by being nonzero, and stores the index within the data
c structure for L/U of that allowed fill.
c
c Finally, if a small pivot element is found in the factorization, it
c is replaced with a small quantity based on machine epsilon, the sign
c of the pivot element, and the norm of the current row.
c
c======================================================================== 
c   Matlab code for computing this factorization is:
c
c      U = zeros(n,n);
c      L = eye(n,n);
c   
c      U(1,1:n) = A(1,1:n);
c   
c      for k = 2:n
c         L(k,1:k-1) = (U(1:k-1,1:k-1)'\(A(k,1:k-1)'))';
c         U(k,k:n) = A(k,k:n) - L(k,1:k-1)*U(1:k-1,k:n);
c      end ;
c   
c   
c======================================================================== 
c
c  The process can be viewed best in the dense case, shown here
c  for k = 4 and n = 6.  Rows 1 through 3 have been computed and
c  are complete at the beginning of this stage (they need no
c  further updates).  The  corresponding dense factorization is
c  usually viewed as consisting of a triangular solve to find the
c  entries L(k,1:k-1), followed by a vector*matrix update to get
c  the entries U(k,k:n).  When a saxpy form is used for both the
c  triangular solve and the vector*matrix update, both phases have
c  the same form and can be computed with the same loop.
c  
c  Stage 0: Copy row 1 of A over to row 1 of LU.
c
c  Stage 1: Compute LU(4,1), as LU(4,1)/U(1,1).
c  
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           | slv |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c  
c  Stage 2: Update LU(4,2:n) using multiplier in LU(4,1) and row 1 of LU:
c  
c           -------------------------------------
c           |     | U12 | U13 | U14 | U15 | U16 |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           | M   | upd | upd | upd | upd | upd |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c  
c  Stage 3: Compute LU(4,2), as LU(4,2)/U(2,2).
c  
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     | slv |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c  
c  Stage 4: Update LU(4,3:n) using multiplier in LU(4,2) and row 2 of LU:
c  
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     | U23 | U24 | U25 | U26 |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     | Ml  | upd | upd | upd | upd |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c  
c  Stage 5: Compute LU(4,3), as LU(4,3)/U(3,3).
c  
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     | slv |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c  
c  Stage 6: Update LU(4,4:n) using multiplier in LU(4,3) and row 3 of LU:
c  
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     | U34 | U35 | U36 |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     | M   | upd | upd | upd |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c  
c  Stage 7: Invert diagonal entry U(4,4), which is now complete:
c  
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     | Inv |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c           |     |     |     |     |     |     |
c           |     |     |     |     |     |     |
c           -------------------------------------
c  
c  Stages 8 onwards: Go on to the next row, repeating above process.
c
c  Note:  When applying updates from previous rows to the
c  current one, we only need apply those for which there is a
c  corresponding nonzero in the current row of LU; otherwise the
c  multiplier would be zero.
c
c     Randall Bramley and Xiaoge Wang
c     Department of Computer Science
c     Indiana University
c     email: bramley@cs.indiana.edu
c     Sun Jun 26 09:51:56 EST 1994
c  
c======================================================================== 
c     a, colind, rwptr    :  Matrix A in CSR format
c     lu, jlu, uptr: Data structure for L/U matrix in MSR format.
c                    lu contains the nonzero values of L/U
c                       lu(1:n) diagonal entries listed in order
c                       lu(n+1) is unused,
c                       lu(n+2:nzlu+1) are the off-diagonal entries in
c                               column increasing order.
c                    jlu contains pointers and column indices.
c                       jlu(1:n) pointers to the beginning of each row
c                                in arrays lu and jlu
c                       jlu(n+1) is a pointer to one beyond the last
c                                entry in lu and jlu
c                       jlu(n+2:nnz+1) are the column indices
c                    uptr(1:n) contains pointers the first entry of U in
c                                each row.
c     n             : order of the matrices.
c     relax         : relaxation parameter for MILU methods.
c     mult          : temporary scalar for holding multiplier L(k,j).
c     ierr          : return error code.
c                        ierr = 0  -> all's OK.
c                        ierr < 0  -> row -ierr had a small pivot
c     k             : loop index for current row number
c     j             : loop index for current row number
c     indj          : index for traversing row k; index of j-th entry in
c                        the lu and jlu arrays.
c     indja         : index for traversing row k of A.
c     inds          : index for traversing updating row of L/U when 
c                        processing row k by applying sparse saxpy ops from
c                        previous (updating) rows to row k (updated row).
c     jluj          : index j of current column number in row k
c     jlus          : index s of current column number in row j
c     ijaj          : column index j for entries of A.
c
c  colptrs(1:n)     : used as an indicator array, to perform updates only
c                     on columns that corr to allowed nonzeros.   If
c                     column j is an allowed nonzero entry of the
c                     current row, then colptrs(j) is the
c                     index in arrays lu(:) and jlu(:) of
c                     the corresponding entry (k,j) of L/U
c  milu             : logical indicating whether to relax diagonal entries
c                     or not
c  rwnrm            : row norm of current row; used for determining small
c		      diagonal element replacement.
c  nrw              : number of nonzeros in current row of A.
c
c======================================================================== 

      implicit none

c     ------------------------
c     Declaration of arguments
c     ------------------------
      integer*4  n, colind(*), rwptr(*), jlu(*), uptr(*)
      real*8  a(*), lu(*), relax
      integer*4  colptrs(n), ierr, nrw

c     ---------------
c     Functions used
c     ---------------
      real*8 dlamch
      external dlamch

c     ---------------
c     Local variables
c     ---------------
      integer*4  k, indj, inds, indja
      integer*4  jluj, jlus, ijaj
      logical  milu
      real*8  SMALL
      real*8  rwnrm, mult

c========================================================================
c      Beginning of Executable Statements
c========================================================================

      SMALL = sqrt(dlamch('E'))

c     ------------------------------------------------------------
c     colptrs is used to hold the indices of entries in LU of 
c     row k.  It is initialized to zero here, and then reset after 
c     each row's work.
c     ------------------------------------------------------------
      do k =  1, n
         colptrs( k ) = 0
      end do

c     ------------------
c     Proceed row by row
c     ------------------
      do 200 k = 1, n

c        --------------------------------------------------------------
c        Set up colptrs with indices in lu of allowed nonzeros of row k
c        --------------------------------------------------------------
         do indj = jlu(k), jlu(k+1)-1
            colptrs(jlu(indj)) = indj
            lu(indj) = 0.0d0
         end do

c        ---------------------------------------------------------
c        Set the diagonal entry (not needed for CSR format of ILU)
c        ---------------------------------------------------------
         colptrs(k) =  k

c        ----------------------------------------------------------------
c        Copy row k of A over to LU.  Note that the "if" test in the loop
c        can be removed if it is known that the sparsity pattern of A
c        is contained in that of LU, which is the case for (M)ILU(s).
c        ----------------------------------------------------------------
         rwnrm = 0.0d0
         nrw = 0
         do indja = rwptr(k), rwptr(k+1)-1
            ijaj = colind(indja)
c            if (colptrs(ijaj) .ne. 0) then
                nrw = nrw + 1
                rwnrm =  rwnrm + abs(a(indja))
                lu(colptrs(ijaj)) = a(indja)
c            end if
         end do

c        -------------------------------------------------------------------
c         The first segment of the next loop on indj effectively solves
c         the transposed upper triangular system
c                 U(1:k-1, 1:k-1)'L(k,1:k-1)' = A(k,1:k-1)'
c         via sparse saxpy operations, throwing away disallowed fill.
c         When the loop index indj reaches the k-th column (i.e., the
c         diagonal entry), then the innermost sparse saxpy operation 
c         effectively is applying the previous updates to the corresponding 
c         part of U via sparse vector*matrix, discarding disallowed fill-in
c         entries.  That operation is 
c            U(k,k:n) = A(k,k:n) - U(1:k-1,k:n)*L(k,1:k-1)
c        -------------------------------------------------------------------

         do 100 indj = jlu(k), uptr(k)-1


c           -----------------------------------------------------------
c           jluj is the col number of current entry in row k of L,
c           and index of diagonal entry of U in same column.  For LU in
c           CSR format, that diag entry will require fancier indexing.
c           -----------------------------------------------------------
            jluj = jlu(indj)

c           -----------------------------------------------------------
c           Solve for next unknown in row k of L: L_kj = L_kj/U_jj ...
c           -----------------------------------------------------------
            lu(indj) = lu(indj)*lu(jluj)
            mult = lu(indj)

c           -------------------------------------------------------------
c           ... and use it as a multiplier to update the entries s in row
c           k of L, s = j+1, ... , k-1, and the entries s in row k of U,
c           s = k, ..., n.
c           -------------------------------------------------------------
            if (milu) then
              do inds = uptr(jluj), jlu(jluj+1)-1
                 jlus = jlu(inds)
                 if (colptrs(jlus) .ne. 0) then
                    lu(colptrs(jlus)) = lu(colptrs(jlus))
     &                              - mult*lu(inds)
                 else
                    lu(k) = lu(k) - relax*mult*lu(inds)
                 end if
              end do
            else
              do inds = uptr(jluj), jlu(jluj+1)-1
                 jlus = jlu(inds)
                 if (colptrs(jlus) .ne. 0) then
                    lu(colptrs(jlus)) = lu(colptrs(jlus))
     &                              - mult*lu(inds)
                 end if
              end do
            end if

 100     continue


c        ----------------------------------------------------------
c        Finished with row k of LU; reset colptrs indices to zero 
c        for next row, and invert diagonal entry of U for this row.
c        ----------------------------------------------------------
         do indj = jlu(k), jlu(k+1)-1
            colptrs(jlu(indj)) = 0
         end do

         colptrs(k) =  0

         if (abs(lu(k)) .le. SMALL*rwnrm/nrw) then
            lu(k) = sign(SMALL*rwnrm/nrw, lu(k))
c            ierr = -k
c            return
         end if
         lu(k) = 1.0d0/lu(k)

 200  continue

      ierr  = 0
      return

c================== End of numfac =====================================

      end


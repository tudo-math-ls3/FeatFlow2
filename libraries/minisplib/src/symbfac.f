c======================================================================== 
      subroutine symbfac(n     ,colind,rwptr  ,levfill,nzmax  ,nzlu   ,
     1                   ijlu  ,uptr   ,lastcol,levels ,rowll  ,
     2                   ierr  ,mneed)
      implicit none

c**********************************************************************
c                   SPLIB Copyright (C) 1995                         **
c                     Indiana University                             **
c**********************************************************************

c======================================================================== 
c======================================================================== 
c                                                                      	*
c  Symbolic factorization of a matrix in compressed sparse row format, 	*
c    with resulting factors stored in a single MSR data structure.     	*
c                                                                      	*
c  This routine uses the CSR data structure of A in two integer vectors	*
c    colind, rwptr to set up the data structure for the ILU(levfill) 	*
c    factorization of A in the integer vectors ijlu and uptr.  Both L	*
c    and U are stored in the same structure, and uptr(i) is the pointer	*
c    to the beginning of the i-th row of U in ijlu.			*
c                                                                      	*
c  The algorithm was originally part of pcgpack, and has been		*
c    modified and adapted for SPLIB.  Most of the changes have been	*
c    documentation and explanation of what is going on.			*
c                                                                      	*
c======================================================================== 
c                                                                      	*
c    Method Used                                                       	*
c    ===========                                                      	*
c                                                                      	*
c  The implementation assumes that the diagonal entries are		*
c  nonzero, and remain nonzero throughout the elimination		*
c  process.  The algorithm proceeds row by row.  When computing		*
c  the sparsity pattern of the i-th row, the effect of row		*
c  operations from previous rows is considered.  Only those		*
c  preceding rows j for which (i,j) is nonzero need be considered,	*
c  since otherwise we would not have formed a linear combination	*
c  of rows i and j.							*
c                                                                      	*
c  The method used has some variations possible.  The definition	*
c  of ILU(s) is not well specified enough to get a factorization	*
c  that is uniquely defined, even in the sparsity pattern that		*
c  results.  For s = 0 or 1, there is not much variation, but for	*
c  higher levels of fill the problem is as follows:  Suppose		*
c  during the decomposition while computing the nonzero pattern		*
c  for row i the following principal submatrix is obtained:		*
c       _______________________						*
c       |          |           |					*
c       |          |           |					*
c       |  j,j     |    j,k    |					*
c       |          |           |					*
c       |__________|___________|					*
c       |          |           |					*
c       |          |           |					*
c       |  i,j     |    i,k    |					*
c       |          |           |					*
c       |__________|___________|					*
c  									*
c  Furthermore, suppose that entry (i,j) resulted from an earlier	*
c  fill-in and has level s1, and (j,k) resulted from an earlier		*
c  fill-in and has level s2:						*
c       _______________________						*
c       |          |           |					*
c       |          |           |					*
c       | level 0  | level s2  |					*
c       |          |           |					*
c       |__________|___________|					*
c       |          |           |					*
c       |          |           |					*
c       | level s1 |           |					*
c       |          |           |					*
c       |__________|___________|					*
c  									*
c  When using A(j,j) to annihilate A(i,j), fill-in will be incurred	*
c  in A(i,k).  How should its level be defined?  It would not be	*
c  operated on if A(i,j) or A(j,m) had not been filled in.  The 	*
c  version used here is to define its level as s1 + s2 + 1.  However,	*
c  other reasonable choices would have been min(s1,s2) or max(s1,s2).	*
c  Using the sum gives a more conservative strategy in terms of the	*
c  growth of the number of nonzeros as s increases.			*
c  									*
c  levels(n+2:nzlu    ) stores the levels from previous rows,		*
c  that is, the s2's above.  levels(1:n) stores the fill-levels		*
c  of the current row (row i), which are the s1's above.		*
c  levels(n+1) is not used, so levels is conformant with MSR format.	*
c  									*
c  Vectors used:							*
c  =============							*
c  									*
c  lastcol(n):								*
c  	The integer lastcol(k) is the row index of the last row		*
c  	to have a nonzero in column k, including the current		*
c  	row, and fill-in up to this point.  So for the matrix		*
c  									*
c             |--------------------------|				*
c             | 11   12           15     |				*
c             | 21   22                26|				*
c             |      32  33   34         |				*
c             | 41       43   44         |				*
c             |      52       54  55   56|				*
c             |      62                66|				*
c             ---------------------------				*
c  									*
c             after step 1, lastcol() = [1  0  0  0  1  0]		*
c             after step 2, lastcol() = [2  2  0  0  2  2]		*
c             after step 3, lastcol() = [2  3  3  3  2  3]		*
c             after step 4, lastcol() = [4  3  4  4  4  3]		*
c             after step 5, lastcol() = [4  5  4  5  5  5]		*
c             after step 6, lastcol() = [4  6  4  5  5  6]		*
c									*  
c          Note that on step 2, lastcol(5) = 2 because there is a	*
c          fillin position (2,5) in the matrix.  lastcol() is used	*
c   	to determine if a nonzero occurs in column j because		*
c   	it is a nonzero in the original matrix, or was a fill.		*
c									*  
c  rowll(n):								*
c  	The integer vector rowll is used to keep a linked list of	*
c  	the nonzeros in the current row, allowing fill-in to be		*
c   	introduced sensibly.  rowll is initialized with the		*
c  	original nonzeros of the current row, and then sorted		*
c  	using a shell sort.  A pointer called head         		*
c  	(what ingenuity) is  initialized.  Note that at any		*
c  	point rowll may contain garbage left over from previous		*
c  	rows, which the linked list structure skips over.		*
c  	For row 4 of the matrix above, first rowll is set to		*
c   	rowll() = [3  1  2  5  -  -], where - indicates any integer.	*
c   	Then the vector is sorted, which yields				*
c   	rowll() = [1  2  3  5  -  -].  The vector is then expanded	*
c  	to linked list form by setting head = 1  and         		*
c   	rowll() = [2  3  5  -  7  -], where 7 indicates termination.	*
c									*  
c  ijlu(nzlu):								*
c  	The returned nonzero structure for the LU factors.		*
c  	This is built up row by row in MSR format, with both L		*
c  	and U stored in the data structure.  Another vector, uptr(n),	*
c  	is used to give pointers to the beginning of the upper		*
c  	triangular part of the LU factors in ijlu.			*
c									*  
c  levels(n+2:nzlu):							*
c  	This vector stores the fill level for each entry from		*
c  	all the previous rows, used to compute if the current entry	*
c  	will exceed the allowed levels of fill.  The value in		*
c  	levels(m) is added to the level of fill for the element in	*
c   	the current row that is being reduced, to figure if 		*
c  	a column entry is to be accepted as fill, or rejected.		*
c  	See the method explanation above.				*
c									*  
c  levels(1:n):								*
c  	This vector stores the fill level number for the current	*
c  	row's entries.  If they were created as fill elements		*
c  	themselves, this number is added to the corresponding		*
c  	entry in levels(n+2:nzlu) to see if a particular column		*
c       entry will							*
c  	be created as new fill or not.  NOTE: in practice, the		*
c  	value in levels(1:n) is one larger than the "fill" level of	*
c  	the corresponding row entry, except for the diagonal		*
c  	entry.  That is why the accept/reject test in the code		*
c  	is "if (levels(j) + levels(m) .le. levfill + 1)".		*
c									*  
c======================================================================== 
c									*
c     31 December 1993							*
c     Randall Bramley							*
c     Department of Computer Science					*
c     Indiana University						*
c     email: bramley@cs.indiana.edu					*
c     Wed Jul 12 15:50:08 EST 1995					*
c									*
c======================================================================== 
c======================================================================== 
c
c on entry:
c========== 
c  n       = The order of the matrix A.
c  ija     = Integer array. Matrix A stored in modified sparse row format.
c  levfill = Integer. Level of fill-in allowed.
c  nzmax   = Integer. The maximum number of nonzero entries in the
c           approximate factorization of a.  This is the amount of storage
c           allocated for ijlu.
c
c on return:
c=========== 
c
c nzlu   = The actual number of entries in the approximate factors, plus one.
c ijlu   = Integer array of length nzlu containing pointers to 
c           delimit rows and specify column number for stored 
c           elements of the approximate factors of a.  the l 
c           and u factors are stored as one matrix.
c uptr   = Integer array of length n containing the pointers to        
c
c ierr is an error flag:
c        ierr  = -i --> near zero pivot in step i
c        ierr  = 0  --> all's OK
c        ierr  = 1  --> not enough storage; check mneed.
c        ierr  = 2  --> illegal parameter
c
c mneed   = contains the actual number of elements in ldu, or the amount
c		of additional storage needed for ldu
c
c work arrays:
c=============
c lastcol    = integer array of length n containing last update of the
c              corresponding column.
c levels     = integer array of length n containing the level of
c              fill-in in current row in its first n entries, and
c              level of fill of previous rows of U in remaining part.
c rowll      = integer array of length n containing pointers to implement a
c              linked list for the fill-in elements.
c
c
c external functions:
c====================
c ifix, float, min0, srtr
c
c======================================================================== 


      integer n,colind(*),rwptr(*),ijlu(*),uptr(*),rowll(*), lastcol(*),
     1        levels(*), levfill,nzmax,nzlu
      integer  ierr,   mneed
      integer icolindj,ijlum,i,j,k,m,ibegin,iend,Ujbeg,Ujend
      integer head,prev,lm,actlev,lowct,k1,k2,levp1,lmk,nzi,rowct

c======================================================================== 
c       Beginning of Executable Statements
c======================================================================== 


c     --------------------------------------------------------------
c     Because the first row of the factor contains no strictly lower
c     triangular parts (parts of L), uptr(1) = ijlu(1) = n+2:
c     --------------------------------------------------------------
      ijlu(1)  =  n+2
      uptr(1)  =  n+2

c     --------------------------------------------------------
c     The storage for the nonzeros of LU must be at least n+1, 
c     for a diagonal matrix:
c     --------------------------------------------------------
      nzlu     =  n+1

c     --------------------------------------------------------------------
c     Number of allowed levels plus 1; used for the test of accept/reject.
c     See the notes about the methodology above.
c     --------------------------------------------------------------------
      levp1    =  levfill + 1

c     -------------------------------------------------------------
c     Initially, for all columns there were no nonzeros in the rows
c     above, because there are no rows above the first one.
c     -------------------------------------------------------------
      do i = 1,n
      	lastcol(i) = 0
      end do

c     -------------------
c     Proceed row by row:
c     -------------------

      do 100 i = 1,n

c       ----------------------------------------------------------
c       Because the matrix diagonal entry is nonzero, the level of
c       fill for that diagonal entry is zero:
c       ----------------------------------------------------------
      	levels(i) = 0

c       ----------------------------------------------------------
c       ibegin and iend are the beginning of rows i and i+1, resp.
c       ----------------------------------------------------------
      	ibegin    =  rwptr(i)
      	iend    =  rwptr(i+1)

c       -------------------------------------------------------------
c       Number of offdiagonal nonzeros in the original matrix's row i
c       -------------------------------------------------------------
      	nzi   =  iend - ibegin

c       --------------------------------------------------------
c       If only the diagonal entry in row i is nonzero, skip the
c       fancy stuff; nothing need be done:
c       --------------------------------------------------------
      	if (nzi .gt. 1) then

c           ----------------------------------------------------------
c           Decrement iend, so that it can be used as the ending index
c           in icolind of row i:
c           ----------------------------------------------------------
            iend          =  iend - 1

c           ---------------------------------------------------------
c           rowct keeps count of the number of nondiagonal entries in
c           the current row:
c           ---------------------------------------------------------
            rowct          =  0

c           ------------------------------------------------------------
c           For nonzeros in the current row from the original matrix A,
c           set lastcol to be the current row number, and the levels of
c           the entry to be 1.  Note that this is really the true level
c           of the element, plus 1.  At the same time, load up the work
c           array rowll with the column numbers for the original entries
c           from row i:
c           ------------------------------------------------------------
            do j = ibegin, iend
               icolindj           =  colind(j)
               lastcol(icolindj)  =  i
               if (icolindj .ne. i) then
                  levels(icolindj)   =  1
                  rowct          =  rowct + 1
                  rowll(rowct)   =  icolindj
               end if
            end do

c           ---------------------------------------------------------
c           Sort the entries in rowll, so that the row has its column
c           entries in increasing order.
c           ---------------------------------------------------------
            call srtr(nzi-1,rowll)

c           ---------------------------------------------------------
c           Now set up rowll as a linked list containing the original
c           nonzero column numbers, as described in the methods section:
c           ---------------------------------------------------------
            head  =  rowll(1)
            k1    =  n+1
            do j = nzi-1, 1, -1
               k2        =  rowll(j)
               rowll(k2) =  k1
               k1        = k2
            end do

c           ------------------------------------------------------------
c           Increment count of nonzeros in the LU factors by the number
c           of nonzeros in the original matrix's row i.  Further
c           incrementing will be necessary if any fill-in actually occurs
c           ------------------------------------------------------------
            nzlu  =  nzlu + nzi - 1

c           ------------------------------------------------------------
c           The integer j will be used as a pointer to track through the
c           linked list rowll:
c           ------------------------------------------------------------
            j  =  head
c
c           ------------------------------------------------------------
c           The integer lowct is used to keep count of the number of
c           nonzeros in the current row's strictly lower triangular part,
c           for setting uptr pointers to indicate where in ijlu the upperc
c           triangular part starts. 
c           ------------------------------------------------------------
            lowct =  0

c           ------------------------------------------------------------
c           Fill-in could only have resulted from rows preceding row i,
c           so we only need check those rows with index j < i.
c           Furthermore, if the current row has a zero in column j,
c           there is no need to check the preceding rows; there clearly
c           could not be any fill-in from those rows to this entry.
c           ------------------------------------------------------------
            do 80 while (j .lt. i)

c              ------------------------------------------------------------
c              Increment lower triangular part count, since in this case
c              (j<i) we got another entry in L:
c              ------------------------------------------------------------
               lowct = lowct  + 1

c              ---------------------------------------------------------
c              If the fill level is zero, there is no way to get fill in
c              occuring.  
c              ---------------------------------------------------------
               if (levfill .ne. 0) then

c                 -----------------------------------------------------
c                 Ujbeg is beginning index of strictly upper triangular
c                 part of U's j-th row, and Ujend is the ending index
c                 of it, in ijlu().
c                 -----------------------------------------------------
                  Ujbeg = uptr(j)
                  Ujend = ijlu(j+1) - 1

c                 -----------------------------------------------------
c                 Need to set pointer to previous entry before working
c                 segment of rowll, because if fill occurs that will be
c                 a moving segment.
c                 -----------------------------------------------------
                  prev  =  j

c                 -----------------------------------------------------
c                 lm is the next nonzero pointer in linked list rowll:
c                 -----------------------------------------------------
                  lm    =  rowll(j)

c                 -------------------------------------------------------
c                 lmk is the fill level in this row, caused by
c                 eliminating column entry j.  That is, level s1 from the
c                 methodology explanation above.
c                 -------------------------------------------------------
                  lmk   =  levels(j)

c                 -------------------------------------------------------
c                 Now proceed through the j-th row of U, because in the
c                 elimination we add a multiple of it to row i to zero
c                 out entry (i,j).  If a column entry in row j of U is
c                 zero, there is no need to worry about fill, because it
c                 cannot cause a fill in the corresponding entry of row i
c                 -------------------------------------------------------
                  do 60  m = Ujbeg, Ujend

c                    ----------------------------------------------------
c                    ijlum is the column number of the current nonzero in
c                    row j of U:
c                    ----------------------------------------------------
                     ijlum =  ijlu(m)
c
c                    ---------------------------------------------------
c                    actlev is the actual level (plus 1) of column entry
c                    j in row i, from summing the level contributions
c                    s1 and s2 as explained in the methods section.
c                    Note that the next line could reasonably be
c                    replaced by, e.g., actlev = max(lmk, levels(m)),
c                    but this would cause greater fill-in:
c                    ---------------------------------------------------
                     actlev = lmk + levels(m)

c                    ---------------------------------------------------
c                    If lastcol of the current column entry in U is not
c                    equal to the current row number i, then the current
c                    row has a zero in column j, and the earlier row j
c                    in U has a nonzero, so possible fill can occur.
c                    ---------------------------------------------------
                     if (lastcol(ijlum) .ne. i) then

c                    --------------------------------------------------
c                    If actlev < levfill + 1, then the new entry has an
c                    acceptable fill level and needs to be added to the
c                    data structure.
c                    --------------------------------------------------
                        if (actlev .le. levp1) then

c                          -------------------------------------------
c                          Since the column entry ijlum in the current
c                          row i is to be filled, we need to update
c                          lastcol for that column number.  Also, the
c                          level number of the current entry needs to be
c                          set to actlev.  Note that when we finish 
c                          processing this row, the n-vector levels(1:n)
c                          will be copied over to the corresponding 
c                          trailing part of levels, so that it can be
c                          used in subsequent rows:
c                          -------------------------------------------

                           lastcol(ijlum) = i
                           levels(ijlum) = actlev

c                          -------------------------------------------
c                          Now find location in the linked list rowll
c                          where the fillin entry should be placed.
c                          Chase through the linked list until the next
c                          nonzero column is to the right of the fill
c                          column number.
c                          -------------------------------------------
                           do 50 while (lm .le. ijlum)
                              prev = lm
                              lm   = rowll(lm)
 50                        continue

c                          -------------------------------------------
c                          Insert new entry into the linked list for
c                          row i, and increase the nonzero count for LU
c                          -------------------------------------------
                           rowll(prev)  = ijlum
                           rowll(ijlum) = lm
                           prev       = ijlum
                           nzlu  =  nzlu  + 1
                        endif

c                    -------------------------------------------------
c                    Else clause is for when lastcol(ijlum) = i.  In
c                    this case, the current column has a nonzero, but
c                    it resulted from an earlier fill-in or from an
c                    original matrix entry.  In this case, need to
c                    update the level number for this column to be the
c                    smaller of the two possible fill contributors,
c                    the current fill number or the computed one from
c                    updating this entry from a previous row.
c                    -------------------------------------------------
                     else
                        levels(ijlum) = min0(levels(ijlum),actlev)
                     endif

c                  -------------------------------------------------
c                  Now go and pick up the next column entry from row
c                  j of U:
c                  -------------------------------------------------
 60                continue

c              -------------------------------------------
c              End if clause for levfill not equal to zero
c              -------------------------------------------
               endif

c              ------------------------------------------------------
c              Pick up next nonzero column index from the linked
c              list, and continue processing the i-th row's nonzeros.
c              This ends the first while loop (j < i).
c              ------------------------------------------------------
               j = rowll(j)

 80         continue

c           ---------------------------------------------------------
c           Check to see if we have exceeded the allowed memory
c           storage before storing the results of computing row i's
c           sparsity pattern into the ijlu and uptr data structures.
c           ---------------------------------------------------------
            if (nzlu .gt. nzmax) then
               mneed = idint((dble(n-i)/dble(2*i))*3*nzlu)
               ierr  = 1
               return
            endif

c           ---------------------------------------------------------
c           Storage is adequate, so update ijlu data structure.
c           Row i ends at nzlu + 1:
c           ---------------------------------------------------------
            ijlu(i+1)   =  nzlu + 1

c           ---------------------------------------------------------
c           ... and the upper triangular part of LU begins at
c           lowct entries to right of where row i begins.
c           ---------------------------------------------------------
            uptr(i)     =  ijlu(i)  + lowct

c           -----------------------------------------------------
c           Now chase through linked list for row i, recording
c           information into ijlu.  At same time, put level data
c           into the levels array for use on later rows:
c           -----------------------------------------------------
            j  =  head
            k1 =  ijlu(i)
            do k  =  k1, nzlu
               ijlu(k)    =  j
               levels(k)  =  levels(j)
               j          =  rowll(j)
            end do

         else

c           ---------------------------------------------------------
c           This else clause ends the (nzi > 1) if.  If nzi = 1, then
c           the update of ijlu and uptr is trivial:
c           ---------------------------------------------------------
            ijlu(i+1)   =  nzlu + 1
            uptr(i)     =  ijlu(i)

         endif

c        ----------------------------------------------
c        And you thought we would never get through....
c        ----------------------------------------------

 100  continue

      ierr = 0
      return

c======================== End of symbfac ==============================

      end
c
c
c======================================================================
c
      subroutine srtr(num,q)
      implicit none
c======================================================================== 
c======================================================================== 
c
c  Implement shell sort, with hardwired increments.  The algorithm for
c  sorting entries in A(0:n-1) is as follows:
c----------------------------------------------------------------
c  inc = initialinc(n)
c  while inc >= 1
c     for i = inc to n-1
c         j = i
c         x = A(i)
c         while j >= inc and A(j-inc) > x
c            A(j) = A(j-inc)
c            j    = j-inc
c         end while 
c         A(j) = x
c     end for
c     inc = nextinc(inc,n)
c  end while
c----------------------------------------------------------------
c
c  The increments here are 1, 4, 13, 40, 121, ..., (3**i - 1)/2, ...
c  In this case, nextinc(inc,n) = (inc-1)/3.  Usually shellsort
c  would have the largest increment the largest integer of the form
c  (3**i - 1)/2 that is less than n, but here it is fixed at 121
c  because most sparse matrices have 121 or fewer nonzero entries
c  per row.  If this routine is expanded for a complete sparse
c  factorization routine, or if a large number of levels of fill is
c  allowed, then possibly it should be replaced with more efficient
c  sorting.
c
c  Any set of increments with 1 as the first one will result in a
c  true sorting algorithm.
c
c     Randall Bramley
c     Department of Computer Science
c     Indiana University
c     email: bramley@cs.indiana.edu
c     Mon Jan 17 20:47:45 EST 1994
c
c======================================================================== 
c
      integer num, q(num)
c
      integer iinc(5), key, icn, ih, ii, i, j, jj
      data iinc/1,4,13,40,121/
c
      if   (num .eq. 0) then
            icn   =  0
      elseif   (num .lt. 14) then
                  icn   =  1
            elseif   (num .lt. 41) then
                  icn   =  2
                  elseif   (num .lt. 122) then
                              icn   =  3
                        elseif   (num .lt. 365) then
                                    icn   =  4
                              else
                                    icn   =  5
      end if
      do 40 ii = 1, icn
            ih = iinc(icn + 1 - ii)
            do 30  j = ih+1, num
                  i = j-ih
                  key = q(j)
                  do 10 jj = 1, j-ih, ih
                        if (key .ge. q(i)) then
                              go to 20
                        else
                              q(i + ih) = q(i)
                              i = i - ih
                        end if
 10               continue
 20               continue
                  q(i + ih) = key
 30         continue
 40   continue
c
      return
      end

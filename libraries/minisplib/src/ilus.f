c======================================================================
      subroutine ilus(n,a,colind,rwptr,
     1                s,relax,
     2                lu,jlu,ilup,
     3                iwork,maxstr,
     4                ierr,mneed)
c**********************************************************************
c                   SPLIB Copyright (C) 1995                         **
c                     Indiana University                             **
c**********************************************************************

c======================================================================== 
c                                                                      	*
c      Incomplete LU factorization for s levels of fill-in.	       	*
c                                                                      	*
c======================================================================== 
c                                                                      	*
c  This routine performs incomplete LU factorization, with s levels	*
c  of fill-in.  It is passed the matrix A in CSR format, and a real*8	*
c  parameter relax, which will indicate whether or not to perform	*
c  MILU factorization.  If the parameter is zero, regular ILU is	*
c  performed instead.  A single integer*4 work array is passed in which	*
c  to create the LU factors data structures, and for intermediate	*
c  storage in the computations.  What is returned are three pointers	*
c  to the start locations in the work array, for the combined MSR	*
c  data structure containing the LU factors.  These pointers will be	*
c  aligned on double word boundaries, provided that the initial word	*
c  in the work array is aligned on a double word boundary.		*
c                                                                      	*
c  If not enough work space is provided, an error code is returned	*
c  and the pointers are not set.  An estimate of the additional		*
c  work space needed is returned in mneed, but this is only a guess.	*
c                                                                      	*
c     Randall Bramley and Xiaoge Wang				       	*
c     Department of Computer Science				       	*
c     Indiana University					       	*
c     email: bramley@cs.indiana.edu				       	*
c     Sun Jan  2 11:09:18 EST 1994				       	*
c									*
c======================================================================== 
c
c parameters
c-----------
c
c on entry:
c========== 
c n       = integer*4. the dimension of the matrix a.
c
c a,colind,rwptr = matrix stored in compressed sparse row format.
c
c s       = integer*4; the fill-in parameter.  Fill-in generated during the
c           factorization is kept only if it is caused by a nonzero element
c           in the matrix that was created on the s level or less.  Fill
c           generated by the original nonzeros are of level 1, those created
c           by such entries are level 2, etc.  The original matrix entries
c           are defined as being of level 0.  See symbfac.f for details.
c
c relax   = relaxation parameter.  relax = 0 gives ilu(s), relax = 1
c           gives milu(s), and values between give the multiplier to use
c           before adding discarded fill to the diagonal.
c
c maxstr  = integer*4 giving amount of integer*4 word space available.
c
c intermediate variables:
c========================
c
c rowll   = pointer to temp vector of levels in a row
c lastcol = pointer to temp vector of last column to update a row
c levels  = pointer to temp vector keeping all levels during symbfac.
c colptrs = pointer to temp vector keeping column pointers in numfac.
c
c on return:
c=========== 
c
c lu, jlu = (pointers) matrix stored in modified sparse row (msr) format containing
c           the l and u factors together. the diagonal (stored in
c           lu(1:n) ) is inverted. each i-th row of the lu, jlu matrix 
c           contains the i-th row of l (excluding the diagonal entry=1) 
c           followed by the i-th row of u.  
c                                                                        
c ilup    = (pointer) integer*4 array of length n containing the pointers to        
c           the beginning of each row of u in the matrix lu, jlu. 
c                                                                       
c  ierr is an error flag:
c        ierr  = -i --> near zero pivot in step i
c        ierr  = 0  --> all's OK
c        ierr  = 1  --> not enough storage; check mneed.
c        ierr  = 2  --> illegal parameter
c
c mneed   = contains number of additional storage locations needed for lu
c           in case of failure, or total integer*4 words used from work array
c           in case of success.
c
c Storage requirements:
c======================
c
c   jlu:	nzlu  		integer*4s.
c   ilup:	n 		integer*4s.
c   lu:		nzlu 		reals
c
c   temporary arrays
c   ----------------
c   rowll       n               integer*4s
c   lastcol     n               integer*4s
c   levels      nzlu            integer*4s
c   colptrs     n               integer*4s
c
c  Warnings:
c===========
c
c 1) A must have nonzero diagonal entries, at least.  This does not assure
c    completion of the factorization, however.
c 2) The pointers ilup, jlu, lu are given as indices to the work array
c    that is passed in.  So if the integer*4 work array is passed in as
c    work(ptr), you will need to offset those pointers by the operations
c        jlu  = jlu  + ptr -1
c        ilup = ilup + ptr -1
c        lu   = lu   + ptr -1
c    on return.
c
c======================================================================== 

      implicit none
      include '../include/intdbl.inc'

      integer*4   n, iwork(*), s,  ierr, rwptr(*), colind(*)
      real*8 a(*), relax
      integer*4  mneed, maxstr, nzlu, remain
      logical milu

c     ----------------------
c     Pointers in work array
c     ----------------------
      integer*4 lu, jlu, ilup
      integer*4 lastcol, rowll, levels, colptrs

      ierr  =  0

c     ------------------------------------------------
c     If requested level of fill is negative, set to 0
c     ------------------------------------------------
      s = max(s,0)

c     ----------------------------------------------
c     If relaxation parameter is zero, method is ILU
c     ----------------------------------------------
      milu = (relax .ne. 0.0d0)

c     ----------------------------------------------------------------
c     Compute pointers into work array.  This version uses two arrays
c     dependent on the (unknown in advance) size nzlu, jlu(1:nzlu) and
c     levels(1:nzlu).  To handle this, first the data space for the
c     work arrays in symbfac is allocated, then half the remainder is
c     assigned to jlu, the rest to levels.  The partitioning of the 
c     integer*4 work array is as:
c
c         [ilup   jlu      ...      levels    rowll lastcol],
c         [-n-----nzlu-----...-------nzlu-------n------n---]
c
c     The index of the last entry in the work array is maxstr, and
c     rowll and lastcol start at the far end of the work array.
c     Note that these are aligned on double word boundaries.
c     ----------------------------------------------------------------
      lastcol = maxstr - n + 1
      rowll = lastcol - n

c     ----------------------------------------------------------
c     Set ilup to the far left of the avaiable work array space.
c     ----------------------------------------------------------
      ilup = 1
      jlu  = ilup + n + mod(INT(n),intdbl)

c     ----------------------------------------------------------------
c     Of space remaining, can allocate half to jlu and half to levels.
c     Compute the halfway dividing mark "remain".  Note that
c     growth available = total - lastcol  - rowll  - ilup 
c     ----------------------------------------------------------------
      remain = (maxstr - 3*n) / 2

c     -----------------------------------------------------
c     Now levels = storage for ilup + half remaining words.
c     -----------------------------------------------------
      levels = n                + remain

c     ------------------------------------------------------------------
c     If remain is nonpositive, estimate 4n as needed storage, 2n for LU
c     and 2n for jlu: this will handle a diagonal preconditioner.
c     ------------------------------------------------------------------
      if (remain .le. 0) then
         ierr = 1
         mneed = 4*n + 1
         print *,'e0: ',ierr,lastcol,rowll,jlu,remain,mneed
         return
      end if

c     -------------------------------
c     Perform symbolic factorization:
c     -------------------------------
      print *,'symb'
      call symbfac(n  , colind    ,rwptr     ,s      ,remain ,nzlu ,
     1     iwork(jlu)   ,iwork(ilup), iwork(lastcol),iwork(levels),
     2     iwork(rowll), ierr  ,mneed)
      print *,'symb ok'

c     -----------------------------------------------------------------
c     If error encountered (probably because of lack of storage space),
c     return.
c     -----------------------------------------------------------------
      if (ierr .ne. 0) then
      	print *,'e1: ',ierr
      	return
      endif

c     -----------------------------------------------------------------
c     Check to see if enough space is available for the lu factorization.
c     The partitioning of the integer*4 work array is as:
c
c        [ilup   jlu           LU        colptrs],
c        [-n-----nzlu-----intdbl*nzlu-------n---]
c
c     -----------------------------------------------------------------
       colptrs = maxstr - n + 1
       remain = maxstr - n - nzlu - n
       if (remain .lt. intdbl*nzlu) then
          ierr = 1
          mneed = intdbl*nzlu - remain
          print *,'e2: ',ierr,mneed
          return
       end if

c     ----------------------------------------
c     Set pointer to LU in integer*4 work array:
c     ----------------------------------------
      lu = jlu + nzlu
      lu = lu + mod(lu-1,intdbl)

c     --------------------------------
c     Perform numerical factorization.
c     --------------------------------
      print *,'numb'
      call numfac(n, colind, rwptr, iwork(jlu),iwork(ilup), a, 
     &              iwork(lu),relax, ierr, iwork(colptrs), milu)
      print *,'numb ok'

c     ----------------------------------------------------
c     Return actual integer*4 words used in LU factors.
c     The last two terms are from forcing alignment of the
c     arrays jlu and alu on double word boundaries.
c     ----------------------------------------------------
      mneed = (1 + intdbl)*nzlu + n + mod(INT(n),intdbl)
     &                          + mod(INT(nzlu),intdbl)
c
      print *,'ok: ',mneed
      return

c==================== End of ILUS ======================================= 

      end

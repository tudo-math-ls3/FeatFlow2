!##############################################################################
!# ****************************************************************************
!# <name> iluk </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# In this module contains routines for preconditioners and solving linear
!# systems Ax=b or LUx=b.
!#
!# The following routines can be found here
!#
!#  1.) iluk_ilu
!#      -> Perfoms an (M)ILU-decomposition of a given matrix A
!#
!#  2.) iluk_symbfac
!#      -> The F90 version of the old minisplib symbolic factorization
!#         performs a symbolic factorization on the matrix A, which
!#         determines the sparsity pattern of the matrix A ,the number
!#         of non-zero entries and the desired fill-in
!#
!#  3.) iluk_numfac
!#      -> Performs the actual decomposition of the matrix A with
!#         the pattern determined by the symbolic factorization and the
!#         desired fill-in. The result are the matrices L and U which are
!#         are stored in a kind of modified spare row format, where we get an
!#         additional pointer to the rows of the matrix U
!#
!#  4.) iluk_lusolt
!#      -> performs a forward and a backward solve for matrices L and U
!#         stored in combined MSR format, with L unit upper triangular
!#         computed by iluk_ilu or iluk_symbfac / iluk_numfac, respectively.
!#
!#  5.) iluk_freeDecomp
!#      -> Release memory
!#
!#</purpose>
!##########################################################################

module iluk

!$use omp_lib
  use fsystem
  use storage
  use linearalgebra

  implicit none

  private

!<types>
!<typeblock>

    ! We get an (M)ILU-Decomposition in the modified sparse row format
    ! so we need :
    !
    ! The entries of the matrices L and U. The library routines
    ! return these entries in such a way that we both the
    ! entries of L and U in one array (p_lu with handle h_lu)
    !
    ! The column numbers of the LU entries, we find that in
    ! jlu. The start indices of the rows of L are also in jlu.
    !
    ! The start indices of the rows of U we find in ilup.

  type t_MILUdecomp


    ! this is an array of length nzlu, where nzlu is the number of
    ! nonzero elements in the LU decomposition

    ! This is a handle to an array that stores the non-zero
    ! entries of the lu-decomposition, the entries of both L and U
    ! are stored in one array
    integer  :: h_lu = ST_NOHANDLE

    ! This is a handle to an index array
    ! that tells us the beginning of the rows of L and
    ! the column indices of the entries
    integer  :: h_jlu = ST_NOHANDLE

    ! This is a handle to an index array where we store pointers to the
    ! starting index of each row of U in the array jlu array
    integer  :: h_ilup = ST_NOHANDLE

    ! The number of nonzero elements in the lu decomposition
    integer  :: nzlu = 0

    ! total number of entries in all arrays
    integer  :: isize = 0

  end type

  public :: t_MILUdecomp

!</typeblock>
!</types>

  public :: iluk_ilu
  public :: iluk_symbfac
  public :: iluk_numfac
  public :: iluk_lusolt
  public :: iluk_freeDecomp

contains

  ! ***************************************************************************

!<subroutine>
    subroutine iluk_freeDecomp(rMILUdecomp)
!<description>
  ! Call this routine to release the memory
  ! allocated in a rMILUdecomp structure
!/description>

!<inputoutput>
    type(t_MILUdecomp),intent(inout) :: rMILUdecomp
!</inputoutput>

!</subroutine>
      ! check if the structure is initialised
      if(rMILUdecomp%h_lu .ne. ST_NOHANDLE) then
        call storage_free(rMILUdecomp%h_lu)
        call storage_free(rMILUdecomp%h_jlu)
        call storage_free(rMILUdecomp%h_ilup)
      end if

    end subroutine ! end iluk_freeDecomp

  ! ***************************************************************************

!<subroutine>

    subroutine iluk_lusolt (n, x, lu, jlu, uptr)

  !<description>
  !  Solve Mx = y.
  !  Perform a forward then backward solve for a MSR matrix containing
  !  a unit lower triangular and an upper triangular matrix with inverted
  !  diagonal, both stored in a single MSR data structure.  The vector of
  !  pointers uptr() contains the row pointers to the beginning of U in
  !  each row.
  !</description>

  !<input>
      ! Column numbers of the entries in the matrix.
      integer, dimension(:), intent(in) :: jlu

      ! Pointer to the entries of U of the LU decomposition
      integer, dimension(:), intent(in) :: uptr

      ! Matrix entries.
      real(DP), dimension(:), intent(in) :: lu
  !</input>

  !<inputoutput>
      ! Result vector
      real(DP), dimension(:), intent(inout) :: x
  !</inputoutput>

      ! local variables
      integer :: n,i,k
      real(DP) :: x_i

      if (n .le. 0)  return

!     ------------
!     Solve Lx = x
!     ------------
      do i = 1, n-1
        x_i =  x(i)

!*$* prefetch_ref=lu(jlu(i+1)),stride=1,kind=rd
!*$* prefetch_ref=lu(jlu(i+1)+16),stride=1,kind=rd
!*$* prefetch_ref=jlu(jlu(i+1)+32),stride=1,kind=rd

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


!     ------------
!     Solve Ux = x
!     ------------
      do i = n, 1, -1
        x_i =  x(i)

!*$* prefetch_ref=lu(jlu(i)-1),stride=1,kind=rd
!*$* prefetch_ref=lu(jlu(i)-17),stride=1,kind=rd
!*$* prefetch_ref=jlu(jlu(i)-33),stride=1,kind=rd

        do k = jlu(i+1)-1, uptr(i), -1
            x_i = x_i - lu(k)*x(jlu(k))
        end do
        x(i) = lu(i)*x_i
      end do

    end subroutine

  ! ***************************************************************************

!**********************************************************************
!                   SPLIB Copyright (C) 1995                         **
!                     Indiana University                             **
!**********************************************************************

!<subroutine>
    subroutine iluk_ilu(n,a,colind,rwptr, ifill,relax,&
                        ierr, rILUDecomp)

!<description>
!  Incomplete LU factorization for s levels of fill-in.
!
!  This routine performs incomplete LU factorization, with s levels
!  of fill-in.  It is passed the matrix A in CSR format, and a real*8
!  parameter relax, which will indicate whether or not to perform
!  MILU factorization.  If the parameter is zero, regular ILU is
!  performed instead.  A single integer*4 work array is passed in which
!  to create the LU factors data structures, and for intermediate
!  storage in the computations.  What is returned are three pointers
!  to the start locations in the work array, for the combined MSR
!  data structure containing the LU factors.  These pointers will be
!  aligned on double word boundaries, provided that the initial word
!  in the work array is aligned on a double word boundary.
!
!  If not enough work space is provided, an error code is returned
!  and the pointers are not set.  An estimate of the additional
!  work space needed is returned in mneed, but this is only a guess.
!
!  <verb>
!     Randall Bramley and Xiaoge Wang
!     Department of Computer Science
!     Indiana University
!     email: bramley@cs.indiana.edu
!     Sun Jan  2 11:09:18 EST 1994
!  </verb>
!
!
!  Warnings:
!
! 1) A must have nonzero diagonal entries, at least.  This does not assure
!    completion of the factorization, however.
!
!========================================================================

!</description>

!<input>
      ! the dimension of the rectangular matrix
      integer, intent(in) ::  n

      ! The fill-in parameter.  Fill-in generated during the
      ! factorization is kept only if it is caused by a nonzero element
      ! in the matrix that was created on the s level or less.  Fill
      ! generated by the original nonzeros are of level 1, those created
      ! by such entries are level 2, etc.  The original matrix entries
      ! are defined as being of level 0.  See symbfac.f for details.
      integer, intent(in) ::  ifill

      ! Matrix entries of the matrix to be factorised.
      real(dp), dimension(:), intent(in) :: a

      ! colind is another csr component and points to the
      ! column indices of A
      integer, dimension(:), intent(in) :: rwptr, colind

      ! relaxation parameter.  relax = 0 gives ilu(s), relax = 1
      ! gives milu(s), and values between give the multiplier to use
      ! before adding discarded fill to the diagonal.
      real(dp), intent(in) :: relax

!</input>

!</output>
      ! A return value that describes the result with
      ! that the routine finished (success, failure, whatever...).
      !   ierr  = -i --> near zero pivot in step i,
      !   ierr  = 0  --> all is OK.
      !   ierr  = 1  --> not enough storage.
      !   ierr  = 2  --> illegal parameter.
      integer, intent(out) :: ierr

      ! A structure that stores the return values of this routine,
      ! which is an ILU-decomposition in a modified sparse row format(MSR).
      type(t_MILUdecomp), intent(inout) :: rILUDecomp
!</output>

!</subroutine>

      ! local variables
      integer :: nzlu
      integer :: mneed,s

!     ----------------------
!     current size of arrays
!     ----------------------
      integer ::  isize

      ! if the relax input parameter is non-zero, we
      ! will perform a milu decomposition
      logical :: milu

      ! add integer handles for the temporary arrays
      integer :: h_rowll, h_lastcol, h_levels, h_colptrs

      integer, dimension(:), pointer :: ph_colptrs
      !real(dp), dimension(:), pointer :: p_lu

      ! pointers to allocated handles
      ! two arrays of length n with pointers to the rows and last column
      integer,  dimension(:), pointer :: ph_rowll, ph_lastcol

      real(dp), dimension(:), pointer :: ph_lu

      integer,  dimension(:), pointer :: ph_jlu, ph_ilup, ph_levels

! intermediate variables:
!========================
!
! rowll   = pointer to temp vector of levels in a row
! lastcol = pointer to temp vector of last column to update a row
! levels  = pointer to temp vector keeping all levels during symbfac.
! colptrs = pointer to temp vector keeping column pointers in numfac.
!
!  ierr is an error flag:
!        ierr  = -i --> near zero pivot in step i
!        ierr  = 0  --> everything is OK
!        ierr  = 1  --> not enough storage; check mneed.
!        ierr  = 2  --> illegal parameter
!
! mneed   = contains number of additional storage locations needed for lu
!           in case of failure, or total integer*4 words used from work array
!           in case of success.
!
! Storage requirements:
!======================
!
!   jlu :	nzlu  integer
!   ilup:	n 		integer
!   lu  :	nzlu	reals
!
!   temporary arrays
!   ----------------
!   rowll       n               integer
!   lastcol     n               integer
!   levels      nzlu            integer
!   colptrs     n               integer

      ierr  =  0

!     ------------------------------------------------
!     If requested level of fill is negative, set to 0
!     ------------------------------------------------
      s = max(ifill,0)

!     ----------------------------------------------
!     If relaxation parameter is zero, method is ILU
!     ----------------------------------------------
      milu = (relax .ne. 0.0_dp)

!    ----------------------------------------------------------------

!     ------------------------------------------------------------------
!     estimate 4n as needed storage, 2n for LU
!     and 2n for jlu: this will handle a diagonal preconditioner.
!     ------------------------------------------------------------------
       mneed = 4*n + 1

!    ---------------OOOOOOKKKKKKKKKKKKAAAAAAAAYYYYYYY--------------------
!    for symbfac I need the rowll and lastcol arrays and sufficient storage
!    for it. so set it up before

      isize = 2*n

      ! get memory for the rowll array length n as in the drawing
      call storage_new('ILU', 'h_rowll', n, ST_INT, h_rowll,ST_NEWBLOCK_ZERO)

      ! get memory for the lastcol array length n as in the drawing
      call storage_new('ILU', 'h_lastcol', n, ST_INT, h_lastcol,ST_NEWBLOCK_ZERO)

      ! get memory for the levels array, estimate 2n as length
      call storage_new('ILU', 'h_levels', isize, ST_INT, h_levels,ST_NEWBLOCK_ZERO)

      ! memory for the result arrays in the structure length n as in the drawing
      call storage_new('ILU', 'h_ilup', n, ST_INT, rILUDecomp%h_ilup, ST_NEWBLOCK_ZERO)

      ! In the beginning we dont know what the number on non-zero elements is, so
      ! we estimate for the array jlu 2*n memory
      call storage_new('ILU', 'h_jlu', isize, ST_INT, rILUDecomp%h_jlu, ST_NEWBLOCK_ZERO)

      ! now get the pointers
      call storage_getbase_int(h_rowll, ph_rowll)
      call storage_getbase_int(h_lastcol, ph_lastcol)
      call storage_getbase_int(h_levels, ph_levels)

      call storage_getbase_int(rILUDecomp%h_jlu, ph_jlu)
      call storage_getbase_int(rILUDecomp%h_ilup, ph_ilup)

      do

        ! call symbolic factorization with the pointers we just created
        ! new function calling style without the work array
        ! here we still need mneed to reallocate if there is not enough memory
        call iluk_symbfac(n  , colind    ,rwptr     ,s      ,isize ,nzlu ,&
             ph_jlu   ,ph_ilup, ph_lastcol, ph_levels, ph_rowll, &
             ierr  ,mneed)

  !     -----------------------------------------------------------------
  !     If error encountered (probably because of lack of storage space),
  !     return.
  !     -----------------------------------------------------------------
        if (ierr .ne. 0) then
          ! increase the amount of memory in the arrays by mneed
          isize = isize + max(mneed,n)
          call storage_realloc('ILU', isize, rILUDecomp%h_jlu, ST_NEWBLOCK_ZERO, .false.)

          call storage_realloc('ILU', isize, h_levels, ST_NEWBLOCK_ZERO, .false.)

          ! the pointers are lost after reallocation so get them back
          call storage_getbase_int(h_levels, ph_levels)

          call storage_getbase_int(rILUDecomp%h_jlu, ph_jlu)
        else
          exit
        end if

      end do ! end do while

      ! here we can throw away unneccessary arrays
      call storage_free(h_levels)
      call storage_free(h_rowll)
      call storage_free(h_lastcol)

      ! we now know the exact amount of memory we need
      ! and we reallocate to use less if time is not critical
      call storage_realloc('ilu2', nzlu, rILUDecomp%h_jlu, ST_NEWBLOCK_ZERO, .true.)
      call storage_getbase_int( rILUDecomp%h_jlu, ph_jlu)
!     -----------------------------------------------------------------
!     Check to see if enough space is available for the lu factorization.
!     The partitioning of the work arrays is as:
!
!     ilup    = integer, dimension(n)
!     jlu     = integer, dimension(nzlu)
!     lu      = double precision, dimension(nzlu)
!     colptrs = integer, dimension(n)
!
!     Note: The new implementation does not write
!     -----------------------------------------------------------------

      ! Symbolic factorization is done
      ! We now know exactely the value of nzlu, now we can allocate memory
      ! to store the entries in the LU decomposition
      call storage_new('ILU', 'h_lu', nzlu, ST_DOUBLE, rILUDecomp%h_lu, ST_NEWBLOCK_ZERO)

      ! get memory for another temporary array
      call storage_new('ILU', 'h_colptrs', n, ST_INT, h_colptrs, ST_NEWBLOCK_ZERO)

      ! get the pointers
      call storage_getbase_int(h_colptrs, ph_colptrs)

      call storage_getbase_double(rILUDecomp%h_lu, ph_lu)

!     --------------------------------
!     Perform numerical factorization.
!     --------------------------------
      call iluk_numfac(n, colind, rwptr, ph_jlu,ph_ilup, a, &
                   ph_lu,relax, ierr, ph_colptrs, milu)

      ! store these parameters to handle the memory mangement
      rILUDecomp%nzlu = nzlu
      rILUDecomp%isize = isize

      ! free temporary arrays
      call storage_free(h_colptrs)

      end subroutine    ! free handles

!========================================================================

!<subroutine>
      subroutine  iluk_numfac(n, colind,rwptr,  jlu, uptr, a, lu,&
                  relax, ierr, colptrs, milu)

!<description>
!<verb>
!**********************************************************************
!                   SPLIB Copyright (C) 1995                         **
!                     Indiana University                             **
!**********************************************************************
!</verb>
!
!
! Numerical factorization with given sparsity pattern identified in
! the MSR data structure  jlu.
!
! Because of the MSR storage format, the method used is a deferred updaet
! version.  At the k-th step, the updates from the previous rows m,
! m < k and m corresponding to a nonzero in row k of L/U, are applied
! to row k.  A relaxation parameter is used.  When it is 0.0, ILU
! factorization is performed.  Otherwise the parameter is multiplied
! times the discarded fillin and added to the diagonal entry.  So when
! the parameter is equal to 1.0, classical MILU factorization is performed.
!
! This routine implicitly assumes that the sparsity pattern specified
! by the MSR data structure contains that of A.  Note the comments in
! the code if you want to remove this restriction.
!
! Matrix A is input in CSR format, and the output MSR array contains
! both L and U, with L assumed to be unit lower triangular.  The diagonal
! entries of U are stored in inverted form to avoid divisions when
! applying the preconditioner.
!
! This version differs from earlier SPLIB versions in performing the
! computations in-place, without using a software accumulator register
! row(:).  The reasoning is that for very large systems, the memory
! access patterns when using such a register will cause many cache
! misses.  The technique (suggested by Xiaoge Wang) is to use an integer*4
! vector colptrs which simultaneously identifies which columns correspond
! to allowed fill by being nonzero, and stores the index within the data
! structure for L/U of that allowed fill.
!
! Finally, if a small pivot element is found in the factorization, it
! is replaced with a small quantity based on machine epsilon, the sign
! of the pivot element, and the norm of the current row.
!
! <verb>
!========================================================================
!   Matlab code for computing this factorization is:
!
!      U = zeros(n,n);
!      L = eye(n,n);
!
!      U(1,1:n) = A(1,1:n);
!
!      for k = 2:n
!         L(k,1:k-1) = (U(1:k-1,1:k-1)'\(A(k,1:k-1)'))';        #(cpp fix: ')
!         U(k,k:n) = A(k,k:n) - L(k,1:k-1)*U(1:k-1,k:n);
!      end ;
!
!
!========================================================================
! </verb>
!
!  The process can be viewed best in the dense case, shown here
!  for k = 4 and n = 6.  Rows 1 through 3 have been computed and
!  are complete at the beginning of this stage (they need no
!  further updates).  The  corresponding dense factorization is
!  usually viewed as consisting of a triangular solve to find the
!  entries L(k,1:k-1), followed by a vector*matrix update to get
!  the entries U(k,k:n).  When a saxpy form is used for both the
!  triangular solve and the vector*matrix update, both phases have
!  the same form and can be computed with the same loop.
!
! <verb>
!  Stage 0: Copy row 1 of A over to row 1 of LU.
!
!  Stage 1: Compute LU(4,1), as LU(4,1)/U(1,1).
!
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           | slv |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!
!  Stage 2: Update LU(4,2:n) using multiplier in LU(4,1) and row 1 of LU:
!
!           -------------------------------------
!           |     | U12 | U13 | U14 | U15 | U16 |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           | M   | upd | upd | upd | upd | upd |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!
!  Stage 3: Compute LU(4,2), as LU(4,2)/U(2,2).
!
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     | slv |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!
!  Stage 4: Update LU(4,3:n) using multiplier in LU(4,2) and row 2 of LU:
!
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     | U23 | U24 | U25 | U26 |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     | Ml  | upd | upd | upd | upd |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!
!  Stage 5: Compute LU(4,3), as LU(4,3)/U(3,3).
!
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     | slv |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!
!  Stage 6: Update LU(4,4:n) using multiplier in LU(4,3) and row 3 of LU:
!
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     | U34 | U35 | U36 |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     | M   | upd | upd | upd |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!
!  Stage 7: Invert diagonal entry U(4,4), which is now complete:
!
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     | Inv |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!           |     |     |     |     |     |     |
!           |     |     |     |     |     |     |
!           -------------------------------------
!
!  Stages 8 onwards: Go on to the next row, repeating above process.
! </verb>
!
!  Note:  When applying updates from previous rows to the
!  current one, we only need apply those for which there is a
!  corresponding nonzero in the current row of LU; otherwise the
!  multiplier would be zero.
!
! <verb>
!     Randall Bramley and Xiaoge Wang
!     Department of Computer Science
!     Indiana University
!     email: bramley@cs.indiana.edu
!     Sun Jun 26 09:51:56 EST 1994
! </verb>
!
! <verb>
!========================================================================
!     a, colind, rwptr    :  Matrix A in CSR format
!     lu, jlu, uptr: Data structure for L/U matrix in MSR format.
!                    lu contains the nonzero values of L/U
!                       lu(1:n) diagonal entries listed in order
!                       lu(n+1) is unused,
!                       lu(n+2:nzlu+1) are the off-diagonal entries in
!                               column increasing order.
!                    jlu contains pointers and column indices.
!                       jlu(1:n) pointers to the beginning of each row
!                                in arrays lu and jlu
!                       jlu(n+1) is a pointer to one beyond the last
!                                entry in lu and jlu
!                       jlu(n+2:nnz+1) are the column indices
!                    uptr(1:n) contains pointers the first entry of U in
!                                each row.
!     n             : order of the matrices.
!     relax         : relaxation parameter for MILU methods.
!     mult          : temporary scalar for holding multiplier L(k,j).
!     ierr          : return error code.
!                        ierr = 0  -> everything is OK.
!                        ierr < 0  -> row -ierr had a small pivot
!     k             : loop index for current row number
!     j             : loop index for current row number
!     indj          : index for traversing row k; index of j-th entry in
!                        the lu and jlu arrays.
!     indja         : index for traversing row k of A.
!     inds          : index for traversing updating row of L/U when
!                        processing row k by applying sparse saxpy ops from
!                        previous (updating) rows to row k (updated row).
!     jluj          : index j of current column number in row k
!     jlus          : index s of current column number in row j
!     ijaj          : column index j for entries of A.
!
!  colptrs(1:n)     : used as an indicator array, to perform updates only
!                     on columns that corr to allowed nonzeros.   If
!                     column j is an allowed nonzero entry of the
!                     current row, then colptrs(j) is the
!                     index in arrays lu(:) and jlu(:) of
!                     the corresponding entry (k,j) of L/U
!  milu             : logical indicating whether to relax diagonal entries
!                     or not
!  rwnrm            : row norm of current row; used for determining small
!		      diagonal element replacement.
!  nrw              : number of nonzeros in current row of A.
! </verb>
!
!========================================================================

!</description>

!     ------------------------
!     Declaration of arguments
!     ------------------------
      integer :: n
      integer, dimension(:), intent(in) :: colind,rwptr, jlu, uptr
      real(dp), dimension(:), intent(in) :: a
      real(dp), dimension(:), intent(out) :: lu
      real(dp), intent(in) :: relax
      integer, dimension(:), intent(out) ::  colptrs
      integer, intent(inout) :: ierr

      ! local variables
      integer :: nrw

!     ---------------
!     Functions used
!     ---------------
      !real(dp) dlamch
      !external dlamch

!</subroutine>

!     ---------------
!     Local variables
!     ---------------
      integer ::  k, indj, inds, indja
      integer ::  jluj, jlus, ijaj
      logical  milu
      real(dp) ::  SMALL
      real(dp) :: rwnrm, mult

!========================================================================
!      Beginning of Executable Statements
!========================================================================

      !SMALL = sqrt(dlamch('E'))
      SMALL = SYS_EPSREAL_DP

!     ------------------------------------------------------------
!     colptrs is used to hold the indices of entries in LU of
!     row k.  It is initialised to zero here, and then reset after
!     each row`s work.
!     ------------------------------------------------------------
      do k =  1, n
         colptrs( k ) = 0
      end do

!     ------------------
!     Proceed row by row
!     ------------------
      do k = 1, n

!        --------------------------------------------------------------
!        Set up colptrs with indices in lu of allowed nonzeros of row k
!        --------------------------------------------------------------
         do indj = jlu(k), jlu(k+1)-1
            colptrs(jlu(indj)) = indj
            lu(indj) = 0.0_dp
         end do

!        ---------------------------------------------------------
!        Set the diagonal entry (not needed for CSR format of ILU)
!        ---------------------------------------------------------
         colptrs(k) =  k

!        ----------------------------------------------------------------
!        Copy row k of A over to LU.  Note that the "if" test in the loop
!        can be removed if it is known that the sparsity pattern of A
!        is contained in that of LU, which is the case for (M)ILU(s).
!        ----------------------------------------------------------------
         rwnrm = 0.0_dp
         nrw = 0
         do indja = rwptr(k), rwptr(k+1)-1
            ijaj = colind(indja)
!            if (colptrs(ijaj) .ne. 0) then
                nrw = nrw + 1
                rwnrm =  rwnrm + abs(a(indja))
                lu(colptrs(ijaj)) = a(indja)
!            end if
         end do

!        -------------------------------------------------------------------
!         The first segment of the next loop on indj effectively solves
!         the transposed upper triangular system
!                 U(1:k-1, 1:k-1)'L(k,1:k-1)' = A(k,1:k-1)'             (cpp fix: ')
!         via sparse saxpy operations, throwing away disallowed fill.
!         When the loop index indj reaches the k-th column (i.e., the
!         diagonal entry), then the innermost sparse saxpy operation
!         effectively is applying the previous updates to the corresponding
!         part of U via sparse vector*matrix, discarding disallowed fill-in
!         entries.  That operation is
!            U(k,k:n) = A(k,k:n) - U(1:k-1,k:n)*L(k,1:k-1)
!        -------------------------------------------------------------------

         do indj = jlu(k), uptr(k)-1


!           -----------------------------------------------------------
!           jluj is the col number of current entry in row k of L,
!           and index of diagonal entry of U in same column.  For LU in
!           CSR format, that diag entry will require fancier indexing.
!           -----------------------------------------------------------
            jluj = jlu(indj)

!           -----------------------------------------------------------
!           Solve for next unknown in row k of L: L_kj = L_kj/U_jj ...
!           -----------------------------------------------------------
            lu(indj) = lu(indj)*lu(jluj)
            mult = lu(indj)

!           -------------------------------------------------------------
!           ... and use it as a multiplier to update the entries s in row
!           k of L, s = j+1, ... , k-1, and the entries s in row k of U,
!           s = k, ..., n.
!           -------------------------------------------------------------
            if (milu) then
              do inds = uptr(jluj), jlu(jluj+1)-1
                 jlus = jlu(inds)
                 if (colptrs(jlus) .ne. 0) then
                    lu(colptrs(jlus)) = lu(colptrs(jlus)) &
                                      - mult*lu(inds)
                 else
                    lu(k) = lu(k) - relax*mult*lu(inds)
                 end if
              end do
            else
              do inds = uptr(jluj), jlu(jluj+1)-1
                 jlus = jlu(inds)
                 if (colptrs(jlus) .ne. 0) then
                    lu(colptrs(jlus)) = lu(colptrs(jlus)) &
                                      - mult*lu(inds)
                 end if
              end do
            end if

         end do ! indj


!        ----------------------------------------------------------
!        Finished with row k of LU; reset colptrs indices to zero
!        for next row, and invert diagonal entry of U for this row.
!         ----------------------------------------------------------
         do indj = jlu(k), jlu(k+1)-1
            colptrs(jlu(indj)) = 0
         end do

         colptrs(k) =  0

         if (abs(lu(k)) .le. SMALL*rwnrm/nrw) then
            lu(k) = sign(SMALL*rwnrm/nrw, lu(k))
!            ierr = -k
!            return
         end if
         lu(k) = 1.0_dp/lu(k)

      end do ! k

      ierr  = 0
      return

      end subroutine
!================== End of numfac =====================================


!======================================================================

!<subroutine>
      subroutine iluk_srtr(num,q)
 !<description>
!  Implement shell sort, with hardwired increments.  The algorithm for
!  sorting entries in A(0:n-1) is as follows:
! <verb>
!----------------------------------------------------------------
!  inc = initialinc(n)
!  while inc >= 1
!     for i = inc to n-1
!         j = i
!         x = A(i)
!         while j >= inc and A(j-inc) > x
!            A(j) = A(j-inc)
!            j    = j-inc
!         end while
!         A(j) = x
!     end for
!     inc = nextinc(inc,n)
!  end while
!----------------------------------------------------------------
! </verb>
!
!  The increments here are 1, 4, 13, 40, 121, ..., (3**i - 1)/2, ...
!  In this case, nextinc(inc,n) = (inc-1)/3.  Usually shellsort
!  would have the largest increment the largest integer*4 of the form
!  (3**i - 1)/2 that is less than n, but here it is fixed at 121
!  because most sparse matrices have 121 or fewer nonzero entries
!  per row.  If this routine is expanded for a complete sparse
!  factorization routine, or if a large number of levels of fill is
!  allowed, then possibly it should be replaced with more efficient
!  sorting.
!
!  Any set of increments with 1 as the first one will result in a
!  true sorting algorithm.
!
! <verb>
!     Randall Bramley
!     Department of Computer Science
!     Indiana University
!     email: bramley@cs.indiana.edu
!     Mon Jan 17 20:47:45 EST 1994
! </verb>
!
!========================================================================
!</description>

     integer :: num
     integer, dimension(num) :: q
!

!</subroutine>

     integer :: key, icn, ih, ii, i, j, jj

     integer, dimension(5) :: iinc = (/1,4,13,40,121/)

!
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
      do ii = 1, icn
        ih = iinc(icn + 1 - ii)
        do j = ih+1, num
          i = j-ih
          key = q(j)
          do jj = 1, j-ih, ih
            if (key .ge. q(i)) then
                  exit
            else
                  q(i + ih) = q(i)
                  i = i - ih
            end if
          end do
          q(i + ih) = key
        end do
      end do

      end subroutine


!========================================================================

!**********************************************************************
!                   SPLIB Copyright (C) 1995                         **
!                     Indiana University                             **
!**********************************************************************

!<subroutine>
  subroutine iluk_symbfac(n     ,colind,rwptr  ,levfill,nzmax  ,nzlu, &
                        ijlu  ,uptr   ,lastcol,levels ,rowll  , &
                        ierr  ,mneed)
  !<description>
  !  Symbolic factorization of a matrix in compressed sparse row format,
  !    with resulting factors stored in a single MSR data structure.
  !
  !  This routine uses the CSR data structure of A in two integer*4 vectors
  !    colind, rwptr to set up the data structure for the ILU(levfill)
  !    factorization of A in the integer*4 vectors ijlu and uptr.  Both L
  !    and U are stored in the same structure, and uptr(i) is the pointer
  !    to the beginning of the i-th row of U in ijlu.
  !
  !  The algorithm was originally part of pcgpack, and has been
  !    modified and adapted for SPLIB.  Most of the changes have been
  !    documentation and explanation of what is going on.
!</description>

  !<input>
      ! The order of the matrix A.
      integer, intent(in) :: n

      ! Column pattern of the matrix to be factorised, in CSR format.
      integer, dimension(:), intent(in) :: colind

      ! Row pattern of the matrix to be factorised, in CSR format.
      integer, dimension(:), intent(in) :: rwptr

      ! Level of fill-in allowed.
      integer, intent(in) :: levfill

      ! The maximum number of nonzero entries in the
      ! approximate factorization of a.  This is the amount of storage
      ! allocated for ijlu.
      integer, intent(in) :: nzmax

    !</input>

    !<output>
      ! The actual number of entries in the approximate factors, plus one.
      integer, intent(out) :: nzlu

      ! MSR structure arrays for the factorised matrix.
      integer, dimension(:), intent(out) :: ijlu,uptr

      ! ierr is an error flag:
      !        ierr  = -i --> near zero pivot in step i.
      !        ierr  = 0  --> everything is OK.
      !        ierr  = 1  --> not enough storage; check mneed.
      !        ierr  = 2  --> illegal parameter.
      integer, intent(out) ::  ierr

      ! Contains the actual number of elements in ldu, or the amount
      ! of additional storage needed for ldu
      integer, intent(out) ::  mneed
    !</output>

    !<inputoutput>
      ! Temporary array. Array of length n containing pointers to implement a
      ! linked list for the fill-in elements.
      integer, dimension(:), intent(inout) :: rowll

      ! Temporary array. Array of length n containing last update of the
      ! corresponding column.
      integer, dimension(:), intent(inout) :: lastcol

      ! Temporary array. Array of length n containing the level of
      ! fill-in in current row in its first n entries, and
      ! level of fill of previous rows of U in remaining part.
      integer, dimension(:), intent(inout) :: levels
    !</inputoutput>

!</subroutine>

      ! local variables
      integer :: icolindj,ijlum,i,j,k,m,ibegin,iend,Ujbeg,Ujend

      integer :: head,prev,lm,actlev,lowct,k1,k2,levp1,lmk,nzi,rowct

  !    Method Used
  !    ===========
  !
  !  The implementation assumes that the diagonal entries are
  !  nonzero, and remain nonzero throughout the elimination
  !  process.  The algorithm proceeds row by row.  When computing
  !  the sparsity pattern of the i-th row, the effect of row
  !  operations from previous rows is considered.  Only those
  !  preceding rows j for which (i,j) is nonzero need be considered,
  !  since otherwise we would not have formed a linear combination
  !  of rows i and j.
  !
  !  The method used has some variations possible.  The definition
  !  of ILU(s) is not well specified enough to get a factorization
  !  that is uniquely defined, even in the sparsity pattern that
  !  results.  For s = 0 or 1, there is not much variation, but for
  !  higher levels of fill the problem is as follows:  Suppose
  !  during the decomposition while computing the nonzero pattern
  !  for row i the following principal submatrix is obtained:
  !       _______________________
  !       |          |           |
  !       |          |           |
  !       |  j,j     |    j,k    |
  !       |          |           |
  !       |__________|___________|
  !       |          |           |
  !       |          |           |
  !       |  i,j     |    i,k    |
  !       |          |           |
  !       |__________|___________|
  !
  !  Furthermore, suppose that entry (i,j) resulted from an earlier
  !  fill-in and has level s1, and (j,k) resulted from an earlier
  !  fill-in and has level s2:
  !       _______________________
  !       |          |           |
  !       |          |           |
  !       | level 0  | level s2  |
  !       |          |           |
  !       |__________|___________|
  !       |          |           |
  !       |          |           |
  !       | level s1 |           |
  !       |          |           |
  !       |__________|___________|
  !
  !  When using A(j,j) to annihilate A(i,j), fill-in will be incurred
  !  in A(i,k).  How should its level be defined?  It would not be
  !  operated on if A(i,j) or A(j,m) had not been filled in.  The
  !  version used here is to define its level as s1 + s2 + 1.  However,
  !  other reasonable choices would have been min(s1,s2) or max(s1,s2).
  !  Using the sum gives a more conservative strategy in terms of the
  !  growth of the number of nonzeros as s increases.
  !
  !  levels(n+2:nzlu    ) stores the levels from previous rows,
  !  that is, the s2`s above.  levels(1:n) stores the fill-levels
  !  of the current row (row i), which are the s1`s above.
  !  levels(n+1) is not used, so levels is conformant with MSR format.
  !
  !  Vectors used:
  !  =============
  !
  !  lastcol(n):
  !  	The integer*4 lastcol(k) is the row index of the last row
  !  	to have a nonzero in column k, including the current
  !  	row, and fill-in up to this point.  So for the matrix
  !
  !             |--------------------------|
  !             | 11   12           15     |
  !             | 21   22                26|
  !             |      32  33   34         |
  !             | 41       43   44         |
  !             |      52       54  55   56|
  !             |      62                66|
  !             ---------------------------
  !
  !             after step 1, lastcol() = [1  0  0  0  1  0]
  !             after step 2, lastcol() = [2  2  0  0  2  2]
  !             after step 3, lastcol() = [2  3  3  3  2  3]
  !             after step 4, lastcol() = [4  3  4  4  4  3]
  !             after step 5, lastcol() = [4  5  4  5  5  5]
  !             after step 6, lastcol() = [4  6  4  5  5  6]
  !
  !          Note that on step 2, lastcol(5) = 2 because there is a
  !          fillin position (2,5) in the matrix.  lastcol() is used
  !   	to determine if a nonzero occurs in column j because
  !   	it is a nonzero in the original matrix, or was a fill.
  !
  !  rowll(n):
  !  	The integer*4 vector rowll is used to keep a linked list of
  !  	the nonzeros in the current row, allowing fill-in to be
  !   	introduced sensibly.  rowll is initialised with the
  !  	original nonzeros of the current row, and then sorted
  !  	using a shell sort.  A pointer called head
  !  	(what ingenuity) is  initialised.  Note that at any
  !  	point rowll may contain garbage left over from previous
  !  	rows, which the linked list structure skips over.
  !  	For row 4 of the matrix above, first rowll is set to
  !   	rowll() = [3  1  2  5  -  -], where - indicates any integer*4.
  !   	Then the vector is sorted, which yields
  !   	rowll() = [1  2  3  5  -  -].  The vector is then expanded
  !  	to linked list form by setting head = 1  and
  !   	rowll() = [2  3  5  -  7  -], where 7 indicates termination.
  !
  !  ijlu(nzlu):
  !  	The returned nonzero structure for the LU factors.
  !  	This is built up row by row in MSR format, with both L
  !  	and U stored in the data structure.  Another vector, uptr(n),
  !  	is used to give pointers to the beginning of the upper
  !  	triangular part of the LU factors in ijlu.
  !
  !  levels(n+2:nzlu):
  !  	This vector stores the fill level for each entry from
  !  	all the previous rows, used to compute if the current entry
  !  	will exceed the allowed levels of fill.  The value in
  !  	levels(m) is added to the level of fill for the element in
  !   	the current row that is being reduced, to figure if
  !  	a column entry is to be accepted as fill, or rejected.
  !  	See the method explanation above.
  !
  !  levels(1:n):
  !  	This vector stores the fill level number for the current
  !  	row`s entries.  If they were created as fill elements
  !  	themselves, this number is added to the corresponding
  !  	entry in levels(n+2:nzlu) to see if a particular column
  !       entry will
  !  	be created as new fill or not.  NOTE: in practice, the
  !  	value in levels(1:n) is one larger than the "fill" level of
  !  	the corresponding row entry, except for the diagonal
  !  	entry.  That is why the accept/reject test in the code
  !  	is "if (levels(j) + levels(m) .le. levfill + 1)".
  !
  !========================================================================
  !
  !     31 December 1993
  !     Randall Bramley
  !     Department of Computer Science
  !     Indiana University
  !     email: bramley@cs.indiana.edu
  !     Wed Jul 12 15:50:08 EST 1995
  !
  !========================================================================

  !========================================================================
  !       Beginning of Executable Statements
  !========================================================================


    ! --------------------------------------------------------------
    ! Because the first row of the factor contains no strictly lower
    ! triangular parts (parts of L), uptr(1) = ijlu(1) = n+2:
    ! --------------------------------------------------------------
    ijlu(1)  =  n+2
    uptr(1)  =  n+2

    ! --------------------------------------------------------
    ! The storage for the nonzeros of LU must be at least n+1,
    ! for a diagonal matrix:
    ! --------------------------------------------------------
    nzlu     =  n+1

    ! --------------------------------------------------------------------
    ! Number of allowed levels plus 1; used for the test of accept/reject.
    ! See the notes about the methodology above.
    ! --------------------------------------------------------------------
    levp1    =  levfill + 1

    ! -------------------------------------------------------------
    ! Initially, for all columns there were no nonzeros in the rows
    ! above, because there are no rows above the first one.
    ! -------------------------------------------------------------
    ! do i = 1,n
      ! lastcol(i) = 0
    ! end do

    ! -------------------
    ! Proceed row by row:
    ! -------------------

     do i = 1,n

      ! ----------------------------------------------------------
      ! Because the matrix diagonal entry is nonzero, the level of
      ! fill for that diagonal entry is zero:
      ! ----------------------------------------------------------
      levels(i) = 0

      ! ----------------------------------------------------------
      ! ibegin and iend are the beginning of rows i and i+1, resp.
      ! ----------------------------------------------------------
      ibegin    =  rwptr(i)
      iend    =  rwptr(i+1)

      ! -------------------------------------------------------------
      ! Number of offdiagonal nonzeros in the original matrix` row i
      ! -------------------------------------------------------------
      nzi   =  iend - ibegin

      ! --------------------------------------------------------
      ! If only the diagonal entry in row i is nonzero, skip the
      ! fancy stuff; nothing need be done:
      ! --------------------------------------------------------
      if (nzi .gt. 1) then

        ! ----------------------------------------------------------
        ! Decrement iend, so that it can be used as the ending index
        ! in icolind of row i:
        ! ----------------------------------------------------------
        iend          =  iend - 1

        ! ---------------------------------------------------------
        ! rowct keeps count of the number of nondiagonal entries in
        ! the current row:
        ! ---------------------------------------------------------
        rowct          =  0

        ! ------------------------------------------------------------
        ! For nonzeros in the current row from the original matrix A,
        ! set lastcol to be the current row number, and the levels of
        ! the entry to be 1.  Note that this is really the true level
        ! of the element, plus 1.  At the same time, load up the work
        ! array rowll with the column numbers for the original entries
        ! from row i:
        ! ------------------------------------------------------------
        do j = ibegin, iend
          icolindj           =  colind(j)
          lastcol(icolindj)  =  i
          if (icolindj .ne. i) then
            levels(icolindj)   =  1
            rowct          =  rowct + 1
            rowll(rowct)   =  icolindj
          end if
        end do

        ! ---------------------------------------------------------
        ! Sort the entries in rowll, so that the row has its column
        ! entries in increasing order.
        ! ---------------------------------------------------------
        call iluk_srtr(nzi-1,rowll)

        ! ---------------------------------------------------------
        ! Now set up rowll as a linked list containing the original
        ! nonzero column numbers, as described in the methods section:
        ! ---------------------------------------------------------
        head  =  rowll(1)
        k1    =  n+1
        do j = nzi-1, 1, -1
          k2        =  rowll(j)
          rowll(k2) =  k1
          k1        = k2
        end do

        ! ------------------------------------------------------------
        ! Increment count of nonzeros in the LU factors by the number
        ! of nonzeros in the original matrix` row i.  Further
        ! incrementing will be necessary if any fill-in actually occurs
        ! ------------------------------------------------------------
        nzlu  =  nzlu + nzi - 1

        ! ------------------------------------------------------------
        ! The integer*4 j will be used as a pointer to track through the
        ! linked list rowll:
        ! ------------------------------------------------------------
        j  =  head

        ! ------------------------------------------------------------
        ! The integer*4 lowct is used to keep count of the number of
        ! nonzeros in the current row`s strictly lower triangular part,
        ! for setting uptr pointers to indicate where in ijlu the upperc
        ! triangular part starts.
        ! ------------------------------------------------------------
        lowct =  0

        ! ------------------------------------------------------------
        ! Fill-in could only have resulted from rows preceding row i,
        ! so we only need check those rows with index j < i.
        ! Furthermore, if the current row has a zero in column j,
        ! there is no need to check the preceding rows; there clearly
        ! could not be any fill-in from those rows to this entry.
        ! ------------------------------------------------------------
        do while (j .lt. i)

        !    ------------------------------------------------------------
        !    Increment lower triangular part count, since in this case
        !    (j<i) we got another entry in L:
        !    ------------------------------------------------------------
             lowct = lowct  + 1

        !    ---------------------------------------------------------
        !    If the fill level is zero, there is no way to get fill in
        !    occuring.
        !    ---------------------------------------------------------
             if (levfill .ne. 0) then

        !       -----------------------------------------------------
        !       Ujbeg is beginning index of strictly upper triangular
        !       part of U`s j-th row, and Ujend is the ending index
        !       of it, in ijlu().
        !       -----------------------------------------------------
                Ujbeg = uptr(j)
                Ujend = ijlu(j+1) - 1

        !       -----------------------------------------------------
        !       Need to set pointer to previous entry before working
        !       segment of rowll, because if fill occurs that will be
        !       a moving segment.
        !       -----------------------------------------------------
                prev  =  j

        !       -----------------------------------------------------
        !       lm is the next nonzero pointer in linked list rowll:
        !       -----------------------------------------------------
                lm    =  rowll(j)

        !       -------------------------------------------------------
        !       lmk is the fill level in this row, caused by
        !       eliminating column entry j.  That is, level s1 from the
        !       methodology explanation above.
        !       -------------------------------------------------------
                lmk   =  levels(j)

        !       -------------------------------------------------------
        !       Now proceed through the j-th row of U, because in the
        !       elimination we add a multiple of it to row i to zero
        !       out entry (i,j).  If a column entry in row j of U is
        !       zero, there is no need to worry about fill, because it
        !       cannot cause a fill in the corresponding entry of row i
        !       -------------------------------------------------------
                do m = Ujbeg, Ujend

        !          ----------------------------------------------------
        !          ijlum is the column number of the current nonzero in
        !          row j of U:
        !          ----------------------------------------------------
                   ijlum =  ijlu(m)

        !          ---------------------------------------------------
        !          actlev is the actual level (plus 1) of column entry
        !          j in row i, from summing the level contributions
        !          s1 and s2 as explained in the methods section.
        !          Note that the next line could reasonably be
        !          replaced by, e.g., actlev = max(lmk, levels(m)),
        !          but this would cause greater fill-in:
        !          ---------------------------------------------------
                   actlev = lmk + levels(m)

        !          ---------------------------------------------------
        !          If lastcol of the current column entry in U is not
        !          equal to the current row number i, then the current
        !          row has a zero in column j, and the earlier row j
        !          in U has a nonzero, so possible fill can occur.
        !          ---------------------------------------------------
                   if (lastcol(ijlum) .ne. i) then

        !          --------------------------------------------------
        !          If actlev < levfill + 1, then the new entry has an
        !          acceptable fill level and needs to be added to the
        !          data structure.
        !          --------------------------------------------------
                      if (actlev .le. levp1) then

        !                -------------------------------------------
        !                Since the column entry ijlum in the current
        !                row i is to be filled, we need to update
        !                lastcol for that column number.  Also, the
        !                level number of the current entry needs to be
        !                set to actlev.  Note that when we finish
        !                processing this row, the n-vector levels(1:n)
        !                will be copied over to the corresponding
        !                trailing part of levels, so that it can be
        !                used in subsequent rows:
        !                -------------------------------------------

                        lastcol(ijlum) = i
                        levels(ijlum) = actlev

        !                -------------------------------------------
        !                Now find location in the linked list rowll
        !                where the fillin entry should be placed.
        !                Chase through the linked list until the next
        !                nonzero column is to the right of the fill
        !                column number.
        !                -------------------------------------------
                         do while (lm .le. ijlum)
                            prev = lm
                            lm   = rowll(lm)
                         end do

        !                -------------------------------------------
        !                Insert new entry into the linked list for
        !                row i, and increase the nonzero count for LU
        !                -------------------------------------------
                         rowll(prev)  = ijlum
                         rowll(ijlum) = lm
                         prev       = ijlum
                         nzlu  =  nzlu  + 1
                      endif

        !          -------------------------------------------------
        !          Else clause is for when lastcol(ijlum) = i.  In
        !          this case, the current column has a nonzero, but
        !          it resulted from an earlier fill-in or from an
        !          original matrix entry.  In this case, need to
        !          update the level number for this column to be the
        !          smaller of the two possible fill contributors,
        !          the current fill number or the computed one from
        !          updating this entry from a previous row.
        !          -------------------------------------------------
                   else
                      levels(ijlum) = min0(levels(ijlum),actlev)
                   endif

        !        -------------------------------------------------
        !        Now go and pick up the next column entry from row
        !        j of U:
        !        -------------------------------------------------
                 end do

        !    -------------------------------------------
        !    End if clause for levfill not equal to zero
        !    -------------------------------------------
             endif

        !    ------------------------------------------------------
        !    Pick up next nonzero column index from the linked
        !    list, and continue processing the i-th row`s nonzeros.
        !    This ends the first while loop (j < i).
        !    ------------------------------------------------------
             j = rowll(j)

        end do

        ! ---------------------------------------------------------
        ! Check to see if we have exceeded the allowed memory
        ! storage before storing the results of computing row i`s
        ! sparsity pattern into the ijlu and uptr data structures.
        ! ---------------------------------------------------------
        if (nzlu .gt. nzmax) then
           mneed =  0.5 * idint((dble(n-i)/dble(2*i))*3*nzlu)
           ierr  = 1
           return
        endif

        ! ---------------------------------------------------------
        ! Storage is adequate, so update ijlu data structure.
        ! Row i ends at nzlu + 1:
        ! ---------------------------------------------------------
        ijlu(i+1)   =  nzlu + 1

        ! ---------------------------------------------------------
        ! ... and the upper triangular part of LU begins at
        ! lowct entries to right of where row i begins.
        ! ---------------------------------------------------------
        uptr(i)     =  ijlu(i)  + lowct

        ! -----------------------------------------------------
        ! Now chase through linked list for row i, recording
        ! information into ijlu.  At same time, put level data
        ! into the levels array for use on later rows:
        ! -----------------------------------------------------
        j  =  head
        k1 =  ijlu(i)
        do k  =  k1, nzlu
           ijlu(k)    =  j
           levels(k)  =  levels(j)
           j          =  rowll(j)
        end do

      else

!       ---------------------------------------------------------
!       This else clause ends the (nzi > 1) if.  If nzi = 1, then
!       the update of ijlu and uptr is trivial:
!       ---------------------------------------------------------
        ijlu(i+1)   =  nzlu + 1
        uptr(i)     =  ijlu(i)

      endif

!     ----------------------------------------------
!     And you thought we would never get through....
!     ----------------------------------------------

    end do

    ierr = 0

  end subroutine

end module

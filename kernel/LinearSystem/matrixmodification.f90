!##############################################################################
!# ****************************************************************************
!# <name> matrixmodification </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of basic matrix modification routines.
!# The routines here work directly on the matrix structure/entries and
!# have no relationship with discretisation routines/information or similar.
!#
!# The following routines can be found in this module:
!#
!# 1.) mmod_replaceLinesByUnit
!#     -> Replaces some rows in a scalar matrix by unit vectors
!#
!# 2.) mmod_replaceLinesByZero
!#     -> Replaces some rows in a scalar matrix by zero vectors
!#
!# 3.) mmod_clearOffdiags
!#     -> Replace all off-diagonal entries in some rows of a matrix
!#        by zero.
!#
!# 4.) mmod_mergeLines
!#        -> Merges some rows in a scalar matrix.
!#
!# 5.) mmod_replaceLinesByUnitBlk
!#     -> Replaces some rows in a block matrix by unit vectors
!#
!# 6.) mmod_expandToFullRow
!#     -> Expands one row of a matrix to a full row
!#
!# 7.) mmod_replaceLineByLumpedMass
!#     -> Replaces one row of a matrix by the entries of the lumped mass
!#        matrix.
!# </purpose>
!##############################################################################

module matrixmodification

!$use omp_lib
  use fsystem
  use storage
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use genoutput

  implicit none

  private

  public :: mmod_replaceLinesByUnit
  public :: mmod_replaceLinesByZero
  public :: mmod_clearOffdiags
  public :: mmod_mergeLines
  public :: mmod_replaceLinesByUnitBlk
  public :: mmod_expandToFullRow
  public :: mmod_replaceLineByLumpedMass
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine mmod_replaceLinesByUnit (rmatrix,Irows)

!<description>
    ! This routine replaces some lines of a given scalar matrix by unit vectors.
!</description>

!<input>
    ! A list of row numbers of all the rows which are to be replaced
    ! by unit vectors.
    integer, intent(in), dimension(:) :: Irows
!</input>

!<inputoutput>
    ! The matrix which is to be modified.
    type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
      call output_line('Virtually transposed matrices not supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mmod_replaceLinesByUnit')
      call sys_halt()
    end if

    ! At first we must take care of the matrix type.
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9)
      call replaceLines_format9 (rmatrix,Irows)
    case (LSYSSC_MATRIX7)
      call replaceLines_format7 (rmatrix,Irows)
    end select

  contains

    ! ****************************************
    ! The replacement routine for format 9

    subroutine replaceLines_format9 (rmatrix,Irows)

      integer, intent(in), dimension(:) :: Irows
      type(t_matrixScalar), intent(inout) :: rmatrix

      ! local variables
      integer :: irow
      real(DP), dimension(:), pointer :: p_DA
      real(SP), dimension(:), pointer :: p_FA
      integer, dimension(:), pointer :: p_Kld,p_Kdiagonal

      ! Get Kld and Kdiagonal
      call lsyssc_getbase_Kld(rmatrix,p_Kld)
      call lsyssc_getbase_Kdiagonal(rmatrix,p_Kdiagonal)

      ! Take care of the format of the entries
      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        ! Get the data array
        call lsyssc_getbase_double(rmatrix,p_DA)
        if (.not. associated(p_DA)) then
          call output_line('Matrix has no data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'replaceLines_format9')
          call sys_halt()
        end if

        ! loop through the rows
        do irow = 1,size(Irows)

          ! Clear the row
          p_DA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_DP

          ! And put a unit vector there
          p_DA(p_Kdiagonal(Irows(irow))) = 1.0_DP

        end do

      case (ST_SINGLE)
        ! Get the data array
        call lsyssc_getbase_single(rmatrix,p_FA)
        if (.not. associated(p_FA)) then
          call output_line('Matrix has no data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'replaceLines_format9')
          call sys_halt()
        end if

        ! loop through the rows
        do irow = 1,size(Irows)

          ! Clear the row
          p_FA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_SP

          ! And put a unit vector there
          p_FA(p_Kdiagonal(Irows(irow))) = 1.0_SP

        end do

      case DEFAULT
        call output_line('Unsupported data format',&
            OU_CLASS_ERROR,OU_MODE_STD,'replaceLines_format9')
        call sys_halt()
      end select

    end subroutine replaceLines_format9

    ! ****************************************
    ! The replacement routine for format 7

    subroutine replaceLines_format7 (rmatrix,Irows)

      integer, intent(in), dimension(:) :: Irows
      type(t_matrixScalar), intent(inout) :: rmatrix

      ! local variables
      integer :: irow
      integer, dimension(:), pointer :: p_Kld
      real(DP), dimension(:), pointer :: p_DA
      real(SP), dimension(:), pointer :: p_FA

      ! Get Kld:
      call lsyssc_getbase_Kld(rmatrix,p_Kld)

      ! Take care of the format of the entries
      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        ! Get the data array
        call lsyssc_getbase_double(rmatrix,p_DA)
        if (.not. associated(p_DA)) then
          call output_line('Matrix has no data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'replaceLines_format7')
          call sys_halt()
        end if

        ! loop through the rows
        do irow = 1,size(Irows)

          ! Put a unit vector there
          p_DA(p_Kld(Irows(irow))) = 1.0_DP

          ! and clear the row
          p_DA(p_Kld(Irows(irow))+1:p_Kld(Irows(irow)+1)-1) = 0.0_DP

        end do

      case (ST_SINGLE)
        ! Get the data array
        call lsyssc_getbase_single(rmatrix,p_FA)
        if (.not. associated(p_FA)) then
          call output_line('Matrix has no data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'replaceLines_format7')
          call sys_halt()
        end if

        ! loop through the rows
        do irow = 1,size(Irows)

          ! Put a unit vector there
          p_FA(p_Kld(Irows(irow))) = 1.0_SP

          ! and clear the row
          p_FA(p_Kld(Irows(irow))+1:p_Kld(Irows(irow)+1)-1) = 0.0_SP

        end do

      case DEFAULT
        call output_line('Unsupported data format',&
            OU_CLASS_ERROR,OU_MODE_STD,'replaceLines_format7')
        call sys_halt()
      end select

    end subroutine replaceLines_format7

  end subroutine mmod_replaceLinesByUnit

  ! ***************************************************************************

!<subroutine>

  subroutine mmod_clearOffdiags (rmatrix,Irows)

!<description>
    ! This routine replaces the offdiagonal entries in some lines of
    ! a given scalar matrix by zero. The diagonal elements are not
    ! changed.
!</description>

!<input>
    ! A list of row numbers of all the rows which are to be replaced
    ! by unit vectors.
    integer, intent(in), dimension(:) :: Irows
!</input>

!<inputoutput>
    ! The matrix which is to be modified.
    type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
      call output_line('Virtually transposed matrices not supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mmod_clearOffdiags')
      call sys_halt()
    end if

    ! At first we must take care of the matrix type.
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9)
      call removeOffdiags_format9 (rmatrix,Irows)
    case (LSYSSC_MATRIX7)
      call removeOffdiags_format7 (rmatrix,Irows)
    end select

  contains

    ! ****************************************
    ! The replacement routine for format 9

    subroutine removeOffdiags_format9 (rmatrix,Irows)

      integer, intent(in), dimension(:) :: Irows
      type(t_matrixScalar), intent(inout) :: rmatrix

      ! local variables
      integer :: irow
      real(DP), dimension(:), pointer :: p_DA
      real(SP), dimension(:), pointer :: p_FA
      integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
      real(DP) :: ddiag,fdiag

      ! Get Kld and Kdiagonal
      call lsyssc_getbase_Kld(rmatrix,p_Kld)
      call lsyssc_getbase_Kdiagonal(rmatrix,p_Kdiagonal)

      ! Take care of the format of the entries
      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        ! Get the data array
        call lsyssc_getbase_double(rmatrix,p_DA)

        ! loop through the rows
        do irow = 1,size(Irows)

          ! Get the diagonal
          ddiag = p_DA(p_Kdiagonal(Irows(irow)))

          ! Clear the row
          p_DA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_DP

          ! restore the diagonal
          p_DA(p_Kdiagonal(Irows(irow))) = ddiag

        end do

      case (ST_SINGLE)
        ! Get the data array
        call lsyssc_getbase_single(rmatrix,p_FA)

        ! loop through the rows
        do irow = 1,size(Irows)

          ! Get the diagonal
          fdiag = p_FA(p_Kdiagonal(Irows(irow)))

          ! Clear the row
          p_FA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_SP

          ! restore the diagonal
          p_FA(p_Kdiagonal(Irows(irow))) = fdiag

        end do

      case DEFAULT
        call output_line('Unsupported data format',&
            OU_CLASS_ERROR,OU_MODE_STD,'removeOffdiags_format9')
        call sys_halt()
      end select

    end subroutine removeOffdiags_format9

    ! ****************************************
    ! The replacement routine for format 7

    subroutine removeOffdiags_format7 (rmatrix,Irows)

      integer, intent(in), dimension(:) :: Irows
      type(t_matrixScalar), intent(inout) :: rmatrix

      ! local variables
      integer :: irow
      integer, dimension(:), pointer :: p_Kld
      real(DP), dimension(:), pointer :: p_DA
      real(SP), dimension(:), pointer :: p_FA

      ! Get Kld:
      call lsyssc_getbase_Kld(rmatrix,p_Kld)

      ! Take care of the format of the entries
      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        ! Get the data array
        call lsyssc_getbase_double(rmatrix,p_DA)

        ! loop through the rows
        do irow = 1,size(Irows)

          ! Clear the row except for the diagonal
          p_DA(p_Kld(Irows(irow))+1:p_Kld(Irows(irow)+1)-1) = 0.0_DP

        end do

      case (ST_SINGLE)
        ! Get the data array
        call lsyssc_getbase_single(rmatrix,p_FA)

        ! loop through the rows
        do irow = 1,size(Irows)

          ! Clear the row except for the diagonal
          p_FA(p_Kld(Irows(irow))+1:p_Kld(Irows(irow)+1)-1) = 0.0_SP

        end do

      case DEFAULT
        call output_line('Unsupported data format',&
            OU_CLASS_ERROR,OU_MODE_STD,'removeOffdiags_format7')
        call sys_halt()
      end select

    end subroutine removeOffdiags_format7

  end subroutine mmod_clearOffdiags

  ! ***************************************************************************

!<subroutine>

  subroutine mmod_replaceLinesByZero (rmatrix,Irows)

!<description>
    ! This routine replaces some lines of a given scalar matrix by unit vectors.
!</description>

!<input>
    ! A list of row numbers of all the rows which are to be replaced
    ! by unit vectors.
    integer, intent(in), dimension(:) :: Irows
!</input>

!<inputoutput>
    ! The matrix which is to be modified.
    type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
      call output_line('Virtually transposed matrices not supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mmod_replaceLinesByZero')
      call sys_halt()
    end if

    ! At first we must take care of the matrix type.
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
      call replaceLinesZero_format97 (rmatrix,Irows)
    end select

  contains

    ! ****************************************
    ! The replacement routine for format 9 and 7

    subroutine replaceLinesZero_format97 (rmatrix,Irows)

      integer, intent(in), dimension(:) :: Irows
      type(t_matrixScalar), intent(inout) :: rmatrix

      ! local variables
      integer :: irow
      real(DP), dimension(:), pointer :: p_DA
      real(SP), dimension(:), pointer :: p_FA
      integer, dimension(:), pointer :: p_Kld,p_Kdiagonal

      ! Get Kld and Kdiagonal
      call lsyssc_getbase_Kld(rmatrix,p_Kld)
      call lsyssc_getbase_Kdiagonal(rmatrix,p_Kdiagonal)

      ! Take care of the format of the entries
      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        ! Get the data array
        call lsyssc_getbase_double(rmatrix,p_DA)

        ! loop through the rows
        do irow = 1,size(Irows)
          ! Clear the row
          p_DA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_DP
        end do

      case (ST_SINGLE)
        ! Get the data array
        call lsyssc_getbase_single(rmatrix,p_FA)

        ! loop through the rows
        do irow = 1,size(Irows)
          ! Clear the row
          p_FA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_SP
        end do

      case DEFAULT
        call output_line('Unsupported data format',&
            OU_CLASS_ERROR,OU_MODE_STD,'replaceLinesZero_format97')
        call sys_halt()
      end select

    end subroutine replaceLinesZero_format97

  end subroutine mmod_replaceLinesByZero

  ! ***************************************************************************

!<subroutine>

  subroutine mmod_mergeLines (rmatrix,Irows,bsymmetric)

!<description>
    ! This routine merges some pairs of lines of a given scalar matrix.
    ! Note that the data (if any) will be cleared. If the optional parameter
    ! bsymmetric=.TRUE., then the sparsity pattern will be symmetric.
!</description>

!<input>
    ! A list of row numbers of all the rows which are to be merged.
    integer, intent(in), dimension(:,:) :: Irows

    ! OPTIONAL: If bsymmetric=.TRUE. the sparsity pattern will be symmetric
    logical, intent(in), optional                    :: bsymmetric
!</input>

!<inputoutput>
    ! The matrix which is to be modified.
    type(t_matrixScalar), intent(inout)              :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
     integer, dimension(:), pointer :: p_Kld,p_ImergeWithRow
     integer, dimension(:), pointer :: p_Kcol
     integer :: ieq,ild,jld,irow,jrow,icol,jcol,naIncr
     integer              :: h_ImergeWithRow

    ! Check, if matrix is not a copy of another matrix or if resize is to be enforced
    if (iand(rmatrix%imatrixSpec, LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0 .or.&
        iand(rmatrix%imatrixSpec, LSYSSC_MSPEC_CONTENTISCOPY)   .ne. 0) then
      call output_line('A copied matrix cannot be modified!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mmod_mergeLines')
        call sys_halt()
    end if

    if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
      call output_line('Virtually transposed matrices not supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mmod_mergeLines')
      call sys_halt()
    end if

    ! Allocate temporal memory
    call storage_new('mmod_mergeLines', 'p_ImergeWithRow', rmatrix%NEQ, ST_INT,&
        h_ImergeWithRow, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(h_ImergeWithRow, p_ImergeWithRow)

    ! Get Kld and Kcol
    call lsyssc_getbase_Kld(rmatrix, p_Kld)
    call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

    ! Loop over the list of rows to be merged
    naIncr = 0
    do ieq = 1, size(Irows,2)

      ! Get row numbers to be merged
      irow = Irows(1, ieq)
      jrow = Irows(2, ieq)

      ! If both numbers are identical, then nothing needs to be done
      if (irow .eq. jrow) cycle

      ! Mark rows for merging
      p_ImergeWithRow(irow) = jrow
      p_ImergeWithRow(jrow) = irow

      ! Loop over row number irow
      loop1: do ild = p_Kld(irow), p_Kld(irow+1)-1
        icol = p_Kcol(ild)
        loop2: do jld = p_Kld(jrow), p_Kld(jrow+1)-1
          if (p_Kcol(jld) .eq. icol) cycle loop1
        end do loop2
        naIncr = naIncr+1
      end do loop1

      ! Loop over row number jrow
      loop3: do jld = p_Kld(jrow), p_Kld(jrow+1)-1
        jcol = p_Kcol(jld)
        loop4: do ild = p_Kld(irow), p_Kld(irow+1)-1
          if (p_Kcol(ild) .eq. jcol) cycle loop3
        end do loop4
        naIncr = naIncr+1
      end do loop3
    end do

    ! Ok, we now know the number of nonzero entries that need to be added to the
    ! global matrix, hence, set new matix dimension and resize Kcol
    rmatrix%NA = rmatrix%NA+naIncr
    call storage_realloc('mmod_mergeLines', rmatrix%NA, rmatrix%h_Kcol,&
        ST_NEWBLOCK_NOINIT, .true.)

    ! At first we must take care of the matrix type and modify the structure
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9, LSYSSC_MATRIX9INTL)
      call mergeLines_format9 (rmatrix, p_ImergeWithRow)

      ! Do we have to generate a symmetric sparsity graph?
      if (present(bsymmetric)) then
        if (bsymmetric) call mergeColumns_format9(rmatrix, p_ImergeWithRow)
      end if

    case (LSYSSC_MATRIX7, LSYSSC_MATRIX7INTL)
      call mergeLines_format7 (rmatrix, p_ImergeWithRow)

      ! Do we have to generate a symmetric sparsity graph?
      if (present(bsymmetric)) then
        if (bsymmetric) call mergeColumns_format7(rmatrix, p_ImergeWithRow)
      end if

    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mmod_mergeLines')
    end select

    ! Free temporal storage
    call storage_free(h_ImergeWithRow)

    ! Next, resize the data accordingly (if present)
    if (rmatrix%h_DA .ne. ST_NOHANDLE) then
      select case (rmatrix%cinterleavematrixFormat)
      case (LSYSSC_MATRIXUNDEFINED)
        call storage_realloc('mmod_mergeLines', rmatrix%NA, rmatrix%h_DA,&
            ST_NEWBLOCK_ZERO, .false.)

      case (LSYSSC_MATRIXD)
        call storage_realloc('mmod_mergeLines', rmatrix%NA*rmatrix%NVAR,&
            rmatrix%h_DA, ST_NEWBLOCK_ZERO, .false.)

      case (LSYSSC_MATRIX1)
        call storage_realloc('mmod_mergeLines', rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR,&
            rmatrix%h_DA, ST_NEWBLOCK_ZERO, .false.)

      case DEFAULT
        call output_line('Unsupported interleave matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mmod_mergeLines')
      end select
    end if

  contains

    ! ****************************************
    ! The row merging routine for format 7

    subroutine mergeLines_format7(rmatrix, ImergeWithRow)
      type(t_matrixScalar), intent(inout)            :: rmatrix
      integer, intent(in), dimension(:) :: ImergeWithRow

      ! local variables
      integer, dimension(:), pointer :: p_Kld
      integer, dimension(:), pointer :: p_Kcol
      integer :: ieq,jeq,ild,jld,icol,jcol,na,naIncr

      ! Get Kld and Kcol
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      ! Initialise current position
      na = rmatrix%NA

      ! Loop over all equations of the matrix but in reverse order (!!!)
      do ieq = rmatrix%NEQ, 1, -1

        ! Ok, so what do we have to do with this row?
        if (ImergeWithRow(ieq) .eq. 0) then

          ! Set current position
          naIncr = na

          ! Just copy the row without modifications
          do ild = p_Kld(ieq+1)-1, p_Kld(ieq), -1
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(ild)

            ! Decrease position counter
            naIncr = naIncr-1
          end do

          ! Update Kld for next column and reduce number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr

        elseif (ImergeWithRow(ieq) .gt. ieq) then

          ! The two rows have already been merged, hence
          ! we can adopt all entries from the indicated row
          jeq = ImergeWithRow(ieq)

          ! Set current position
          naIncr = na

          ! Just copy the row without modifications
          do jld = p_Kld(jeq+1)-1, p_Kld(jeq), -1
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(jld)

            ! Decrease position counter
            naIncr = naIncr-1
          end do

          ! Sort row according to format 7 format. First we need to find the position of the diagonal entry
          ! and swap the diagonal entry from the copied row with the entry of the current row. Then we need
          ! to move the diagonal entry from the copied row to its correct position.
          sort: do ild = naIncr+1, na
            if (p_Kcol(ild) .eq. ieq) then
              ! Swap the first entry which corresponds to the diagonal entry of the copied row with current position
              p_Kcol(ild)      = jeq
              p_Kcol(naIncr+1) = ieq

              ! Move the swapped diagonal entry from the copied row to its correct position
              do jld = ild+1, na
                if (p_Kcol(jld) .gt. jeq) exit sort
                p_Kcol(jld-1) = p_Kcol(jld)
                p_Kcol(jld)   = jeq
              end do
              exit sort
            end if
          end do sort

          ! Update Kld for next column and reduce number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr

        else

          ! We actually have to merge the two rows
          jeq = ImergeWithRow(ieq)

          ! First, sort the row that should be merged into the current row in ascending order.
          ! This is mandatory since the matrix is stored in format 7 so that the diagonal
          ! entry is at the first position and not stored in-between.
          do jld = p_Kld(jeq)+1, p_Kld(jeq+1)-1
            if (p_Kcol(jld) .gt. jeq) exit
            p_Kcol(jld-1) = p_Kcol(jld)
            p_Kcol(jld)   = jeq
          end do

          ! Ok, now the row is in ascending order
          ild = p_Kld(ieq+1)-1
          jld = p_Kld(jeq+1)-1

          ! Set current position
          naIncr = na

          ! Loop in reverse order until both rows have been processed
          ! Note that the diagonal entry of row irow is treated below
          do while((ild .gt. p_Kld(ieq)) .and. (jld .ge. p_Kld(jeq)))

            ! Get column indices for both rows and recude position counter by one
            icol = p_Kcol(ild)
            jcol = p_Kcol(jld)

            ! Which is the next column that should be processed
            if (jcol .gt. icol) then
              p_Kcol(naIncr) = jcol
              jld = jld-1
            elseif (jcol .eq. icol) then
              p_Kcol(naIncr) = jcol
              ild = ild-1
              jld = jld-1
            else
              p_Kcol(naIncr) = icol
              ild = ild-1
            end if

            ! Decrease position counter
            naIncr = naIncr-1
          end do

          ! Copy leftover from the row that corresponds to JEQ?
          do while (jld .ge. p_Kld(jeq))
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(jld)

            ! Decrease position counter
            naIncr = naIncr-1
            jld    = jld-1
          end do

          ! Copy leftover from the row that corresponds to IEQ?
          ! Here, we explicitly copy the diagonal entry to the first position.
          do while (ild .ge. p_Kld(ieq))
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(ild)

            ! Decrease position counter
            naIncr = naIncr-1
            ild    = ild-1
          end do

          ! That is it! Update Kld for next column and adjust number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr
        end if

        ! Ok, if we process two consecutive rows, then we must also adjust the starting position
        if ((ImergeWithRow(ieq) .eq. ieq-1) .and. (ieq .ne. 1)) p_Kld(ieq) = na+1
      end do

      ! Consistency check. If na is not zero then something went wrong
      if (na .ne. 0) then
        call output_line('An internal error occured; please contact the Featflow team!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mergeLines_format7')
        call sys_halt()
      end if
    end subroutine mergeLines_format7


    ! ****************************************
    ! The row merging routine for format 9

    subroutine mergeLines_format9(rmatrix, ImergeWithRow)
      type(t_matrixScalar), intent(inout)            :: rmatrix
      integer, intent(in), dimension(:) :: ImergeWithRow

      ! local variables
      integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
      integer, dimension(:), pointer :: p_Kcol
      integer :: ieq,jeq,ild,jld,icol,jcol,na,naIncr

      ! Get Kld, Kcol and Kdiagonal
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

      ! Initialise current position
      na = rmatrix%NA

      ! Loop over all equations of the matrix but in reverse order (!!!)
      do ieq = rmatrix%NEQ, 1, -1

        ! Ok, so what do we have to do with this row?
        if (ImergeWithRow(ieq) .eq. 0) then

          ! Set current position
          naIncr = na

          ! Just copy the row without modifications
          do ild = p_Kld(ieq+1)-1, p_Kld(ieq), -1
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(ild)

            ! Check for diagonal entry
            if (p_Kcol(naIncr) .eq. ieq) p_Kdiagonal(ieq) = naIncr

            ! Decrease position counter
            naIncr = naIncr-1
          end do

          ! Update Kld for next column and adjust number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr

        elseif (ImergeWithRow(ieq) .gt. ieq) then

          ! The two rows have already been merged, hence
          ! we can adopt all entries from the indicated row
          jeq = ImergeWithRow(ieq)

          ! Set current position
          naIncr = na

          ! Just copy the row without modifications
          do jld = p_Kld(jeq+1)-1, p_Kld(jeq), -1
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(jld)

            ! Check for diagonal entry
            if (p_Kcol(naIncr) .eq. ieq) p_Kdiagonal(ieq) = naIncr

            ! Decrease position counter
            naIncr = naIncr-1
          end do

          ! Update Kld for next column and adjust number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr

        else

          ! We actually have to merge the two rows
          jeq = ImergeWithRow(ieq)

          ! The rows are in ascending order
          ild = p_Kld(ieq+1)-1
          jld = p_Kld(jeq+1)-1

          ! Set current position
          naIncr = na

          ! Loop in reverse order until both rows have been processed
          do while((ild .ge. p_Kld(ieq)) .and. (jld .ge. p_Kld(jeq)))

            ! Get column indices for both rows and recude position counter by one
            icol = p_Kcol(ild)
            jcol = p_Kcol(jld)

            ! Which is the next column that should be processed?
            if (jcol .gt. icol) then
              p_Kcol(naIncr) = jcol
              jld = jld-1
            elseif (jcol .eq. icol) then
              p_Kcol(naIncr) = jcol
              ild = ild-1
              jld = jld-1
            else
              p_Kcol(naIncr) = icol
              ild = ild-1
            end if

            ! Check for diagonal entry
            if (p_Kcol(naIncr) .eq. ieq) p_Kdiagonal(ieq) = naIncr

            ! Decrease position counter
            naIncr = naIncr-1
          end do

          ! Copy leftover from the row that corresponds to JEQ?
          do while (jld .ge. p_Kld(jeq))
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(jld)

            ! Check for diagonal entry
            if (p_Kcol(naIncr) .eq. ieq) p_Kdiagonal(ieq) = naIncr

            ! Decrease position counter
            naIncr = naIncr-1
            jld    = jld-1
          end do

          ! Copy leftover from the row that corresponds to IEQ?
          do while (ild .ge. p_Kld(ieq))
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(ild)

            ! Check for diagonal entry
            if (p_Kcol(naIncr) .eq. ieq) p_Kdiagonal(ieq) = naIncr

            ! Decrease position counter
            naIncr = naIncr-1
            ild    = ild-1
          end do

          ! That is it! Update Kld for next column and adjust number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr
        end if

        ! Ok, if we process two consecutive rows, then we must also adjust the starting position
        if ((ImergeWithRow(ieq) .eq. ieq-1) .and. (ieq .ne. 1)) p_Kld(ieq) = na+1
      end do

      ! Consistency check. If na is not zero then something went wrong
      if (na .ne. 0) then
        call output_line('An internal error occured; please contact the Featflow team!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mergeLines_format9')
        call sys_halt()
      end if
    end subroutine mergeLines_format9


    ! ****************************************
    ! The column merging routine for format 7

    subroutine mergeColumns_format7(rmatrix, ImergeWithRow)
      type(t_matrixScalar), intent(inout)            :: rmatrix
      integer, intent(in), dimension(:) :: ImergeWithRow

      ! local variables
      integer, dimension(:), pointer :: p_Kld
      integer, dimension(:), pointer :: p_KldAux
      integer, dimension(:), pointer :: p_Kcol,p_KcolAux
      integer :: ieq,ild,jld,jjld,icol,jcol,naIncr,idxIncr,iidxIncr
      integer              :: h_KldAux,h_KcolAux

      ! Allocate temporal memory
      call storage_new('mergeColumns_format7', 'p_KldAux', rmatrix%NEQ+1, ST_INT,&
          h_KldAux, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(h_KldAux, p_KldAux)

      ! Get Kld and Kcol
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      ! Compute number of nonzero entries that need to be inserted into the matrix.
      ! First, compute the number of nonzero entries present in each row.
      do ieq = 1, rmatrix%NEQ
        p_KldAux(ieq) = p_Kld(ieq+1)-p_Kld(ieq)
      end do

      ! Next, subtract the number of nonzero entries present in each column.
      do ild = 1, rmatrix%NA
        icol = p_Kcol(ild)
        p_KldAux(icol) = p_KldAux(icol)-1
      end do

      ! If an entry is zero, then the number of nonzero entries in the corresponding
      ! row equals the number of nonzero entries in the corresponding column and
      ! nothing needs to be done. A negative entry indicates that the number of
      ! nonzero entries in the corresponding column exceeds the number of nonzero
      ! entries in the corresponding column. Hence, we must fill the missing entries.
      naIncr = 0
      do ieq = 1, rmatrix%NEQ
        p_KldAux(ieq) = -min(p_KldAux(ieq), 0)
        naIncr        =  naIncr+p_KldAux(ieq)
      end do

      ! Return of no new entries need to be considered
      if (naIncr .eq. 0) then
        call storage_free(h_KldAux)
        return
      end if

      ! Make a copy of the original column indices
      h_KcolAux = ST_NOHANDLE
      call storage_copy(rmatrix%h_Kcol, h_KcolAux)
      call storage_getbase_int(h_KcolAux, p_KcolAux)

      ! Ok, we now know the number of nonzero entries that need to be added to the
      ! global matrix, hence, set new matix dimension and resize Kcol
      rmatrix%NA = rmatrix%NA+naIncr
      call storage_realloc('mergeLines_format7', rmatrix%NA, rmatrix%h_Kcol,&
          ST_NEWBLOCK_NOINIT, .true.)
      call lsyssc_getbase_Kcol(rmatrix,p_Kcol)

      ! Adjust the row separator and copy its original content to p_KldAux
      idxIncr     = p_KldAux(1)+p_Kld(2)-p_Kld(1)
      p_KldAux(1) = p_Kld(1)

      do ieq = 2, rmatrix%NEQ
        ! Compute number of nonzero entries in subsequent row
        iidxIncr = p_KldAux(ieq)+p_Kld(ieq+1)-p_Kld(ieq)

        ! Copy original entry
        p_KldAux(ieq) = p_Kld(ieq)

        ! Adjust absolut position for current row
        p_Kld(ieq) = p_Kld(ieq-1)+idxIncr

        ! Swap increment
        idxIncr = iidxIncr
      end do
      p_KldAux(rmatrix%NEQ+1) = p_Kld(rmatrix%NEQ+1)
      p_Kld(rmatrix%NEQ+1)    = p_Kld(rmatrix%NEQ)+idxIncr

      ! In the first step, move the column data to their new positions
      do ieq = rmatrix%NEQ, 1, -1
        ! Compute offset to new position
        idxIncr = p_Kld(ieq)-p_KldAux(ieq)

        ! Copy data
        if (idxIncr .gt. 0) then
          do ild = p_KldAux(ieq+1)-1, p_KldAux(ieq), -1
            p_Kcol(ild+idxIncr) = p_Kcol(ild)
          end do
        end if
      end do

      ! Loop over all matrix rows
      loop1: do ieq = 1, rmatrix%NEQ

        ! If the row has not been merged then skip it
        if (ImergeWithRow(ieq) .eq. 0) cycle

        ! Otherwise, we process all entries of the current row
        ! and check if the corresponding column entries exist
        loop2: do ild = p_KldAux(ieq)+1, p_KldAux(ieq+1)-1

          ! Get column number
          icol = p_KcolAux(ild)

          ! Skip diagonal entries
          if (icol .eq. ieq) cycle loop2

          ! Loop over row number icol and check if entry ieq exists
          loop3: do jld = p_Kld(icol)+1, p_Kld(icol+1)-1

            ! Get column number
            jcol = p_Kcol(jld)

            ! Did we find the entry?
            if (jcol .eq. ieq) exit

            ! Skip left off-diagonal entries
            if (jcol .lt. ieq) cycle

            ! No, the entry was not found
            do jjld = p_Kld(icol+1)-1, jld+1, -1
              p_Kcol(jjld) = p_Kcol(jjld-1)
            end do
            p_Kcol(jld) = ieq
            exit
          end do loop3
        end do loop2
      end do loop1

      ! Free temporal storage
      call storage_free(h_KldAux)
      call storage_free(h_KcolAux)
    end subroutine mergeColumns_format7


    ! ****************************************
    ! The column merging routine for format 9

    subroutine mergeColumns_format9(rmatrix, ImergeWithRow)
      type(t_matrixScalar), intent(inout)            :: rmatrix
      integer, intent(in), dimension(:) :: ImergeWithRow

      ! local variables
      integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
      integer, dimension(:), pointer :: p_KldAux
      integer, dimension(:), pointer :: p_Kcol,p_KcolAux
      integer :: ieq,ild,jld,jjld,icol,jcol,idxIncr,iidxIncr
      integer              :: h_KldAux,h_KcolAux

      ! Allocate temporal memory
      call storage_new('mergeColumns_format9', 'p_KldAux', rmatrix%NEQ+1, ST_INT,&
          h_KldAux, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(h_KldAux, p_KldAux)

      ! Get Kld and Kcol
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      ! Compute number of nonzero entries that need to be inserted into the matrix.
      ! First, compute the number of nonzero entries present in each row.
      do ieq = 1, rmatrix%NEQ
        p_KldAux(ieq) = p_Kld(ieq+1)-p_Kld(ieq)
      end do

      ! Next, subtract the number of nonzero entries present in each column.
      do ild = 1, rmatrix%NA
        icol = p_Kcol(ild)
        p_KldAux(icol) = p_KldAux(icol)-1
      end do

      ! If an entry is zero, then the number of nonzero entries in the corresponding
      ! row equals the number of nonzero entries in the corresponding column and
      ! nothing needs to be done. A negative entry indicates that the number of
      ! nonzero entries in the corresponding column exceeds the number of nonzero
      ! entries in the corresponding column. Hence, we must fill the missing entries.
      naIncr = 0
      do ieq = 1, rmatrix%NEQ
        p_KldAux(ieq) = -min(p_KldAux(ieq), 0)
        naIncr        =  naIncr+p_KldAux(ieq)
      end do

      ! Return of no new entries need to be considered
      if (naIncr .eq. 0) then
        call storage_free(h_KldAux)
        return
      end if

      ! Make a copy of the original column indices
      h_KcolAux = ST_NOHANDLE
      call storage_copy(rmatrix%h_Kcol, h_KcolAux)
      call storage_getbase_int(h_KcolAux, p_KcolAux)

      ! Ok, we now know the number of nonzero entries that need to be added to the
      ! global matrix, hence, set new matix dimension and resize Kcol
      rmatrix%NA = rmatrix%NA+naIncr
      call storage_realloc('mergeLines_format9', rmatrix%NA, rmatrix%h_Kcol,&
          ST_NEWBLOCK_NOINIT, .true.)
      call lsyssc_getbase_Kcol(rmatrix,p_Kcol)

      ! Adjust the row separator and copy its original content to p_KldAux
      idxIncr     = p_KldAux(1)+p_Kld(2)-p_Kld(1)
      p_KldAux(1) = p_Kld(1)

      do ieq = 2, rmatrix%NEQ
        ! Compute number of nonzero entries in subsequent row
        iidxIncr = p_KldAux(ieq)+p_Kld(ieq+1)-p_Kld(ieq)

        ! Copy original entry
        p_KldAux(ieq) = p_Kld(ieq)

        ! Adjust absolut position for current row
        p_Kld(ieq) = p_Kld(ieq-1)+idxIncr

        ! Swap increment
        idxIncr = iidxIncr
      end do
      p_KldAux(rmatrix%NEQ+1) = p_Kld(rmatrix%NEQ+1)
      p_Kld(rmatrix%NEQ+1)    = p_Kld(rmatrix%NEQ)+idxIncr

      ! In the first step, move the column data to their new positions
      do ieq = rmatrix%NEQ, 1, -1
        ! Compute offset to new position
        idxIncr = p_Kld(ieq)-p_KldAux(ieq)

        ! Copy data
        if (idxIncr .gt. 0) then
          do ild = p_KldAux(ieq+1)-1, p_KldAux(ieq), -1
            p_Kcol(ild+idxIncr) = p_Kcol(ild)
          end do
        end if
      end do

      ! Loop over all matrix rows
      loop1: do ieq = 1, rmatrix%NEQ

        ! If the row has not been merged then skip it
        if (ImergeWithRow(ieq) .eq. 0) cycle

        ! Otherwise, we process all entries of the current row
        ! and check if the corresponding column entries exist
        loop2: do ild = p_KldAux(ieq), p_KldAux(ieq+1)-1

          ! Get column number
          icol = p_KcolAux(ild)

          ! Skip diagonal entries
          if (icol .eq. ieq) cycle loop2

          ! Loop over row number icol and check if entry ieq exists
          loop3: do jld = p_Kld(icol), p_Kld(icol+1)-1

            ! Get column number
            jcol = p_Kcol(jld)

            ! Did we find the entry?
            if (jcol .eq. ieq) exit

            ! Skip left off-diagonal entries
            if (jcol .lt. ieq) cycle

            ! No, the entry was not found
            do jjld = p_Kld(icol+1)-1, jld+1, -1
              p_Kcol(jjld) = p_Kcol(jjld-1)
            end do
            p_Kcol(jld) = ieq
            exit
          end do loop3
        end do loop2
      end do loop1

      ! Rebuild the diagonal array
      call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
      call lsyssc_rebuildKdiagonal(p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)

      ! Free temporal storage
      call storage_free(h_KldAux)
      call storage_free(h_KcolAux)
    end subroutine mergeColumns_format9
  end subroutine mmod_mergeLines

  ! ***************************************************************************

!<subroutine>

  subroutine mmod_replaceLinesByUnitBlk (rmatrix,iblockRow,Irows)

  !<description>
    ! This routine replaces some lines in a given block matrix by unit vectors.
  !</description>

  !<input>
    ! Number of the block row where some lines should be replaced by unit
    ! vectors
    integer, intent(in) :: iblockRow

    ! A list of row numbers of all the rows in block row iblockRow which are
    ! to be replaced by unit vectors. These numbers are not global DOF`s
    ! but the starting indices relative to the block row iblockRow
    ! (e.g. "1" identifies the first row in the block row iblockrow).
    integer, intent(in), dimension(:) :: Irows
  !</input>

  !<inputoutput>
    ! The matrix which is to be modified.
    type(t_matrixBlock), intent(inout) :: rmatrix
  !</inputoutput>

  !</subroutine>

    integer :: icol

    ! Loop through all column matrices
    do icol = 1,rmatrix%nblocksPerCol
      if (lsysbl_isSubmatrixPresent(rmatrix,iblockRow,icol)) then
        ! Replace by unit or zero vectors
        if (icol .eq. iblockRow) then
          call mmod_replaceLinesByUnit (rmatrix%RmatrixBlock(iblockRow,icol),Irows)
        else
          call mmod_replaceLinesByZero (rmatrix%RmatrixBlock(iblockRow,icol),Irows)
        end if
      end if
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mmod_expandToFullRow (rmatrix,irow)

!<description>
    ! Expands one row of a matrix to a 'full' row.
    ! For Full matrices, this has no effect.
    ! For sparse matrices, row irow is expanded to a full row.
    ! If data is present in the matrix, it is appropriately re-organised.
!</description>

!<input>
    ! Row number of a row to be expanded to a full row.-
    integer, intent(in) :: irow
!</input>

!<inputoutput>
    ! The matrix which is to be modified.
    type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Kdiagonal
    real(DP), dimension(:), pointer :: p_Da
    integer :: ncol,icoldiff,i,j,k
    
    real(DP), dimension(:), allocatable :: Da
    integer, dimension(:), allocatable :: col
    
    if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
      call output_line("Virtually transposed matrices not supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"mmod_expandToFullRow")
      call sys_halt()
    end if

    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX1)
      ! Full matrix, nothing to be done.

    case (LSYSSC_MATRIX9)
      ! CSR format.

      ! Get the row pointer. Figure out the number of elements
      ! in row irow.
      call lsyssc_getbase_Kld (rmatrix,p_Kld)
      ncol = p_Kld(irow+1)-p_Kld(irow)
      icoldiff = rmatrix%NCOLS - ncol
      
      if (icoldiff .gt. 0) then

        ! Check, if matrix is not a copy of another matrix or if resize is to be enforced
        if (iand(rmatrix%imatrixSpec, LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0) then
          call output_line("A copied matrix cannot be modified!",&
              OU_CLASS_ERROR,OU_MODE_STD,"mmod_expandToFullRow")
            call sys_halt()
        end if

        if (iand(rmatrix%imatrixSpec, LSYSSC_MSPEC_CONTENTISCOPY) .ne. 0) then
          call output_line("A copied matrix cannot be modified!",&
              OU_CLASS_ERROR,OU_MODE_STD,"mmod_expandToFullRow")
            call sys_halt()
        end if

        ! We have to expand; reallocate memory.
        
        allocate(Da(ncol))
        allocate(Col(ncol))
      
        call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
        call lalg_copyVector (p_Kcol(p_Kld(irow):p_Kld(irow+1)-1),Col)
        
        call storage_realloc ("mmod_expandToFullRow",&
            rmatrix%NA+icoldiff,rmatrix%h_Kcol,ST_NEWBLOCK_NOINIT)

        ! Restructure the Kcol array, move the data.
        ! We *on-purpose* use storage_getbase here since lsyssc_getbase
        ! does not work -- it would take rmatrix%NA into account!
        call storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
        call copyVectorInt (p_Kcol,&
            p_Kld(irow+1),p_Kld(rmatrix%NEQ+1)-1,p_Kld(irow+1)+icoldiff)

        ! Restructure the data if there is any
        if (lsyssc_hasMatrixContent(rmatrix)) then

          if (rmatrix%cdataType .ne. ST_DOUBLE) then
            call output_line("Unsupported data type",&
                OU_CLASS_ERROR,OU_MODE_STD,"mmod_expandToFullRow")
            call sys_halt()
          end if

          call lsyssc_getbase_double (rmatrix,p_Da)
          call lalg_copyVector (p_Da(p_Kld(irow):p_Kld(irow+1)-1),Da)

          call storage_realloc ("mmod_expandToFullRow",&
              rmatrix%NA+icoldiff,rmatrix%h_Da,ST_NEWBLOCK_NOINIT)

          ! Move the data.
          ! We *on-purpose* use storage_getbase here since lsyssc_getbase
          ! does not work -- it would take rmatrix%NA into account!
          call storage_getbase_double (rmatrix%h_Da,p_Da)
          call copyVectorDP (p_Da,&
              p_Kld(irow+1),p_Kld(rmatrix%NEQ+1)-1,p_Kld(irow+1)+icoldiff)

        end if

        ! Restructure the Kld array.
        call lalg_vectorAddScalar(p_Kld(irow+1:rmatrix%NEQ+1),icoldiff)
        
        ! Restructure the Kdiagonal array.
        ! icoldiff must be added to the entry of row irow+1,... end,
        ! and the entry of row irow has to be calculated appropriately.
        call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
        call lalg_vectorAddScalar(p_Kdiagonal(irow+1:rmatrix%NEQ),icoldiff)
        p_Kdiagonal(irow) = p_Kld(irow) + irow - 1

        ! Take care of the column structure (p_Kcol) of the row number 'irow'
        ! Since we have a full row in 'irow' then, the p_Kcol array
        ! for this row, simply starts from 1 to rmatrix%NCOLS
        i = p_Kld(irow)-1
        do j=1,p_Kld(irow+1)-p_Kld(irow)
          p_Kcol(j+i) = j
        end do
        
        ! Initialise the new line with zero,
        ! Copy the entries from the old line into the new line.
        if (lsyssc_hasMatrixContent(rmatrix)) then
        
          ! Initialise with zero
          do j=p_Kld(irow), p_Kld(irow+1)-1
            p_Da(j) = 0.0_DP
          end do

          ! Take care of the matrix data (p_Da) of the row number 'irow'
          ! Keep the original values      
          k = p_Kld(irow)-1
          do i=1,ncol
            do j=k+1,p_Kld(irow+1)-1
              if (Col(i) .eq. p_Kcol(j)) then
                p_Da(j) = Da(i)
                k = j
                exit
              end if
            end do
          end do
          
        end if
                
        deallocate(Col)
        deallocate(Da)
        
        ! Fix the matrix size.
        ! After this, lsyssc_getbase_XXXX works again.
        rmatrix%NA = rmatrix%NA+icoldiff

      end if

    case default
      call output_line("Unsupported matrix format!",&
          OU_CLASS_ERROR,OU_MODE_STD,"mmod_expandToFullRow")
      call sys_halt()
    end select

  contains

  ! ***************************************************************************

    subroutine copyVectorInt (Ix,istartpos,iendpos,inewpos)

    ! Copies an integer vector Ix(inewpos:) = Ix(istartpos,iendpos).
    ! Copies backwards; used to copy parts in a vector
    ! to a higher address.

    ! Source vector
    integer, dimension(:), intent(inout) :: Ix

    ! Start position
    integer, intent(in) :: istartpos

    ! End position
    integer, intent(in) :: iendpos
    
    ! New start position
    integer, intent(in) :: inewpos

    integer :: i,ioffset

      ioffset = inewpos-istartpos
      do i = iendpos,istartpos,-1
        Ix(i+ioffset) = Ix(i)
      end do

    end subroutine


  ! ***************************************************************************

    subroutine copyVectorDP (Dx,istartpos,iendpos,inewpos)

    ! Copies an real vector Dx(inewpos:) = Dx(istartpos,iendpos).
    ! Copies backwards; used to copy parts in a vector
    ! to a higher address.

    ! Source vector
    real(DP), dimension(:), intent(inout) :: Dx

    ! Start position
    integer, intent(in) :: istartpos

    ! End position
    integer, intent(in) :: iendpos
    
    ! New start position
    integer, intent(in) :: inewpos

    integer :: i,ioffset

      ioffset = inewpos-istartpos
      do i = iendpos,istartpos,-1
        Dx(i+ioffset) = Dx(i)
      end do

    end subroutine

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine mmod_replaceLineByLumpedMass (rmatrix,irow,rmassMatrix)

!<description>
  ! Replace all the elements of the row number 'irow' of 'rmatrix' 
  ! by the values of the lumped mass matrix.
  ! The routine could be used to implement the
  ! zero mean value constraint.
  !
  ! REMARK: If necessary, the structure of the matrix is automatically
  ! extended to be able to hold the complete irow-th row!
!</description>

!<input>
  ! Row number of a row to be expanded to a full row and
  ! replaced by the diagonal values of rMassMatrix
  integer, intent(in) :: irow

  ! The lumped mass matrix.
  type(t_matrixScalar), intent(in) :: rmassMatrix
!</input>

!<inputoutput>
  ! The matrix which is to be modified.
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_KdiagMass,p_Kld
    real(DP), dimension(:), pointer :: p_Da,p_DaMass
    integer :: i,j
    
    ! Matrices must have the same size
    if (rmatrix%NEQ .ne. rmassMatrix%NEQ) then
      call output_line('Matrices not compatible, different size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mmod_replaceLineByLumpedMass')
      call sys_halt()
    end if    
    
    ! First expand row number 'irow' of the 'rmatrix' to a full row
    call mmod_expandToFullRow (rmatrix,irow)
    
    select case (rmatrix%cmatrixFormat)
    
    case (LSYSSC_MATRIX9)
      ! Check the format of our lumped mass matrix 'rMassMatrix'
      select case (rmassMatrix%cmatrixFormat)

      case (LSYSSC_MATRIX9)
        ! CSR format.

        ! Get the diagonal structure and values of the lumped mass matrix
        call lsyssc_getbase_double (rmassMatrix,p_DaMass)
        call lsyssc_getbase_Kdiagonal (rmassMatrix,p_KdiagMass)
        
        ! Get the row structure and data of the rmatrix
        call lsyssc_getbase_Kld (rmatrix,p_Kld)      
        call lsyssc_getbase_double (rmatrix,p_Da)      
        
        j = p_Kld(irow)-1
        do i = p_Kld(irow), p_Kld(irow+1)-1
          p_Da(i) = p_DaMass(p_KdiagMass(i-j))
        end do

      case (LSYSSC_MATRIXD)
        ! Diagonal format
        
        ! Get the diagonal values of the lumped mass matrix
        call lsyssc_getbase_double (rmassMatrix,p_DaMass)
        
        ! Get the row structure and data of the rmatrix
        call lsyssc_getbase_Kld (rmatrix,p_Kld)      
        call lsyssc_getbase_double (rmatrix,p_Da)
        call lalg_copyVector (p_DaMass(:),p_Da(p_Kld(irow):p_Kld(irow+1)-1))
        
      case default
        call output_line("Unsupported mass matrix format!",&
            OU_CLASS_ERROR,OU_MODE_STD,"mmod_replaceLineByLumpedMass")
        call sys_halt()
        
      end select
      
    case default
      call output_line("Unsupported destination matrix format!",&
          OU_CLASS_ERROR,OU_MODE_STD,"mmod_replaceLineByLumpedMass")
      call sys_halt()

    end select
    
  end subroutine

end module matrixmodification

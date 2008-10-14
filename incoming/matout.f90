!#########################################################################
!# ***********************************************************************
!# <name> matout </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains all routines and constant necessary to output 
!# matrices of different formats.
!#
!# The following routines can be found in this module:
!#
!# 1.) matout_writeFullMatrix
!#     -> writes a full double precision matrix into a text file
!#
!# 2.) matout_writeSparseMatrix
!#     -> writes a sparse double precision matrix into a text file
!#
!# </purpose>
!#########################################################################

module matout

  use fsystem
  use storage
  use io
  
  implicit none

  contains

!<subroutine>
  subroutine matout_writeFullMatrix (Da, nrow, ncol, sarray, &
                                     bnoZero, cfile, sfile, sformat)
  
  !<description>
    ! Write full double precision matrix into a text file.
    !
    ! This writes an array Da with nrow rows and ncol columns as a matrix
    ! into a text file.
  !</description>
    
  !<input>
    ! matrix: array [1..nrow,1..ncol] of double
    real(DP), dimension(nrow,ncol), intent(IN) :: Da
    
    ! number of rows
    integer(I32), intent(IN) :: nrow
    
    ! number of columns
    integer(I32), intent(IN) :: ncol
    
    ! name of the matrix
    character(len=*), intent(IN) :: sarray
    
    ! suppress zeroes in output: yes/no
    logical, intent(IN) :: bnoZero
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel cfile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer(I32), intent(IN) :: cfile
    
    ! name of the file where to write to. Only relevant for cfile=0!
    character(len=*), intent(IN) :: sfile
    
    ! format string to use for the output; e.g. '(D20.10)'
    character(len=*), intent(IN) :: sformat
    
  !</input>
    
!</subroutine>
    
    !local variables
    integer(I32) :: i, j, cf, nchar
    real(DP) :: dval
    character(len=128) :: S
    character(len=6) :: sformatChar
    
    if (cfile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE)
      if (cf .eq. -1) then
        print *, 'matout_writeFullMatrix: Could not open file '// &
                 trim(sfile)
        stop
      end if
    else
      cf = cfile
    end if
    
    ! Get length of output strings
    S(:) = ' '
    write (S,sformat) 0.0_DP
    nchar = len(trim(S))
    
    ! Build array format string
    sformatChar = '(A'//sys_i3(nchar)//')'
    
    ! Write all format strings into the file
    write (cf,'(3A15,L2,3I15)') sarray, sformat, sformatChar, &
                                bnoZero, nchar, nrow, ncol
    
    if (nrow .le. 0) return
    if (ncol .le. 0) return

    ! Write the matrix
    do i=1, nrow
      do j=1, ncol-1
        dval = Da(i,j)
        if ((.not. bnoZero) .or. (dval .ne. 0.0_DP)) then
          write (cf,sformat,ADVANCE='NO') dval
        else
          write (cf,sformatChar, ADVANCE='NO') '.'
        end if
      end do
      dval = Da(i,j)
      if ((.not. bnoZero) .or. (dval .ne. 0.0_DP)) then
        write (cf,sformat) dval
      else
        write (cf,sformatChar) '.'
      end if
    end do
    
    ! Close the file if necessary
    if (cfile .eq. 0) close(cf)
  
  end subroutine matout_writeFullMatrix

!<subroutine>
  subroutine matout_writeSparseMatrix (Da, Icol, Irow, &
                                       nrow, ncol, sarray, &
                                       bnoZero, cfile, sfile, sformat)
  
  !<description>
    ! Write sparse double precision matrix into a text file.
    !
    ! This writes an array Da with nrow rows and ncol columns as a matrix
    ! into a text file.
  !</description>
    
  !<input>
    ! matrix: array [1..na] of double
    real(DP), dimension(:), intent(IN) :: Da
    
    ! Column structure of the matrix
    integer(I32), dimension(:), intent(IN) :: Icol
    
    ! Row structure of the matrix
    integer(I32), dimension(nrow+1), intent(IN) :: Irow
    
    ! number of rows
    integer(I32), intent(IN) :: nrow
    
    ! number of columns; must be =nrow for structure-7 matrices
    integer(I32), intent(IN) :: ncol
    
    ! name of the matrix
    character(len=*), intent(IN) :: sarray
    
    ! suppress zeroes in output: yes/no
    logical, intent(IN) :: bnoZero
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel cfile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer(I32), intent(IN) :: cfile
    
    ! name of the file where to write to. Only relevant for cfile=0!
    character(len=*), intent(IN) :: sfile
    
    ! format string to use for the output; e.g. '(D20.10)'
    character(len=*), intent(IN) :: sformat
    
  !</input>
    
!</subroutine>
    
    !local variables
    integer(I32) :: i, j, k, cf, nchar
    real(DP) :: dval
    character(len=128) :: S
    character(len=6) :: sformatChar
    integer :: h_DrowVec
    real(DP), dimension(:), pointer :: p_DrowVec
    
    if (cfile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE)
      if (cf .eq. -1) then
        print *, 'matout_writeFullMatrix: Could not open file '// &
                 trim(sfile)
        stop
      end if
    else
      cf = cfile
    end if
    
    ! Get length of output strings
    S(:) = ' '
    write (S,sformat) 0.0_DP
    nchar = len(trim(S))
    
    ! Build array format string
    sformatChar = '(A'//sys_i3(nchar)//')'
    
    ! Write all format strings into the file
    write (cf,'(3A15,L2,3I15)') sarray, sformat, sformatChar, &
                                bnoZero, nchar, nrow, ncol

    if (nrow .le. 0) return
    if (ncol .le. 0) return
    
    ! Write the matrix
    call storage_new('matout_writeSparseMatrix','DrowVec',ncol, &
                     ST_DOUBLE, h_DrowVec, ST_NEWBLOCK_NOINIT)
    call storage_getbase_double(h_DrowVec, p_DrowVec)
    do i=1, nrow
      ! Extract row i
      do j=1, ncol
        p_DrowVec(j) = 0.0_DP
      end do
      do j=0, Irow(i+1)-Irow(i)-1
        k = Irow(i)+j
        p_DrowVec(Icol(k)) = Da(k)
        if (Da(k) .eq. 0.0_DP) p_DrowVec(Icol(k)) = 1D-99
      end do
      ! Write row i
      do j=1, ncol-1
        dval = p_DrowVec(j)
        if ((.not. bnoZero) .or. (dval .ne. 0.0_DP)) then
          if (abs(dval) .le. 1.0D-90) then
            write (cf,sformatChar, ADVANCE='NO') '0'
          else
            write (cf,sformat,ADVANCE='NO') dval
          end if
        else
          write (cf,sformatChar, ADVANCE='NO') '.'
        end if
      end do
      dval = p_DrowVec(ncol)
      if ((.not. bnoZero) .or. (dval .ne. 0.0_DP)) then
        if (abs(dval) .le. 1.0D-90) then
          write (cf,sformatChar, ADVANCE='NO') '0'
        else
          write (cf,sformat,ADVANCE='NO') dval
        end if
      else
        write (cf,sformatChar, ADVANCE='NO') '.'
      end if
    end do
    call storage_free(h_DrowVec)
    
    ! Close the file if necessary
    if (cfile .eq. 0) close(cf)
  
  end subroutine matout_writeSparseMatrix

end module

!#########################################################################
!# ***********************************************************************
!# <name> matrixio </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains all routines and constant necessary to output
!# matrices of different formats.
!#
!# The following routines can be found in this module:
!#
!# 1.) matio_writeMatrixHR
!#     -> Writes a matrix in human readable form into a text file
!#
!# 2.) matio_writeBlockMatrixHR
!#     -> Writes a block matrix in human readable form into a text file
!#
!# 3.) matio_spyMatrix
!#     -> Writes a scalar matrix in CSR format into a file which can be
!#        visualised by means of the MATLAB command SPY
!#
!# 4.) matio_spyBlockMatrix
!#     -> Writes a block matrix in CSR format into a file which can be
!#        visualised by means of the MATLAB command SPY
!#
!# 5.) matio_writeMatrixMaple
!#     -> Writes a matrix into a text file using the MAPLE syntax
!#
!# 2.) matio_writeBlockMatrixMaple
!#     -> Writes a block matrix in MAPLE format into a text file
!#
!# The following auxiliary functions can be found here:
!#
!# 1.) matio_writeMatrix1_Dble
!#     -> writes a full double precision matrix into a text file
!#
!# 2.) matio_writeMatrix79_Dble
!#     -> writes a sparse double precision matrix (format 7/9)
!#        into a text file
!#
!# 3.) matio_writeMatrix1_Sngl
!#     -> writes a full single precision matrix into a text file
!#
!# 4.) matio_writeMatrix79_Sngl
!#     -> writes a sparse single precision matrix (format 7/9)
!#        into a text file
!#
!# </purpose>
!#########################################################################

module matrixio

!$use omp_lib
  use fsystem
  use genoutput
  use storage
  use io
  use linearsystemscalar
  use linearsystemblock
  use globalsystem
  
  implicit none
  
  private
  
  public :: matio_writeMatrixHR
  public :: matio_writeBlockMatrixHR
  public :: matio_spyMatrix
  public :: matio_spyBlockMatrix
  public :: matio_writeMatrixMaple
  public :: matio_writeBlockMatrixMaple
  public :: matio_writeMatrix1_Dble
  public :: matio_writeMatrix79_Dble
  public :: matio_writeMatrix1_Sngl
  public :: matio_writeMatrix79_Sngl
  
contains

  ! ***************************************************************************

!<subroutine>
  subroutine matio_writeBlockMatrixHR (rmatrix, sarray,&
                                       bnoZero, ifile, sfile, sformat, dthreshold)
  
  !<description>
    ! This routine writes a block matrix into a text file.
    ! The matrix is written in human readable form.
    ! Note that for this purpose, a new matrix is temporarily created in memory!
  !</description>
    
  !<input>
    ! The matrix to be written out
    type(t_matrixBlock), intent(in) :: rmatrix
    
    ! Name of the matrix
    character(len=*), intent(in) :: sarray
    
    ! Suppress zeroes in output: yes/no
    logical, intent(in) :: bnoZero
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile
    
    ! Format string to use for the output; e.g. '(E20.10)'
    character(len=*), intent(in) :: sformat
    
    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    ! If not present, values are not replaced (i.e. the default value is 0.0).
    real(DP), intent(in), optional :: dthreshold
  !</input>
    
!</subroutine>

    ! local variables
    type(t_matrixBlock) :: rtempMatrix

    ! We have to create a global matrix first!
    call glsys_assembleGlobal (rmatrix,rtempMatrix,.true.,.true.)
                              
    ! Write matrix to the file
    call matio_writeMatrixHR (rtempMatrix%RmatrixBlock(1,1), sarray,&
                              bnoZero, ifile, sfile, sformat,dthreshold)

    ! Release the temporary matrix
    call lsysbl_releaseMatrix (rtempMatrix)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine matio_writeMatrixHR (rmatrix, sarray,&
                                  bnoZero, ifile, sfile, sformat, dthreshold)
  
  !<description>
    ! This routine writes a scalar matrix into a text file.
    ! The matrix is written in human readable form.
  !</description>
    
  !<input>
    ! The matrix to be written out
    type(t_matrixScalar), intent(in) :: rmatrix
    
    ! Name of the matrix
    character(len=*), intent(in) :: sarray
    
    ! Suppress zeroes in output: yes/no
    logical, intent(in) :: bnoZero
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile
    
    ! Format string to use for the output; e.g. '(E20.10)'
    character(len=*), intent(in) :: sformat
    
    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    ! If not present, values are not replaced (i.e. the default value is 0.0).
    real(DP), intent(in), optional :: dthreshold

  !</input>
    
!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_DA
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kld
  real(dp) :: dthres

  ! Replace small values by zero
  dthres = 0.0_DP
  if (present(dthreshold)) dthres = dthreshold

  ! Depending on the matrix format, choose the right routine for writing
  select case (rmatrix%cmatrixFormat)
  case (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
    ! Matrix precision?
    select case (rmatrix%cdataType)
    case (ST_DOUBLE)
      ! Get the data arrays and write the matrix
      if (.not.lsyssc_hasMatrixContent(rmatrix)) then
        call output_line('Matrix has no data',&
            OU_CLASS_ERROR,OU_MODE_STD,'matio_writeMatrixHR')
        call sys_halt()
      end if
      call lsyssc_getbase_double (rmatrix,p_Da)
      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)

      call matio_writeMatrix79_Dble (p_Da, p_Kcol, p_Kld, &
                                      rmatrix%NEQ, rmatrix%NCOLS, sarray, &
                                      bnoZero, ifile, sfile, sformat,dthres)
    case DEFAULT
      call output_line ('Unsupported matrix precision!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'matio_writeFullMatrix')
      call sys_halt()
    end select
  case DEFAULT
    call output_line ('Unknown matrix format!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'matio_writeFullMatrix')
    call sys_halt()
  end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine matio_writeMatrix1_Dble (Da, nrow, ncol, sarray, &
                                      bnoZero, ifile, sfile, sformat,dthreshold)
  
  !<description>
    ! Write full double precision matrix into a text file.
    !
    ! This writes an array Da with nrow rows and ncol columns as a matrix
    ! into a text file.
  !</description>
    
  !<input>
    ! number of rows
    integer, intent(in) :: nrow
    
    ! number of columns
    integer, intent(in) :: ncol
    
    ! matrix: array [1..nrow,1..ncol] of double
    real(DP), dimension(nrow,ncol), intent(in) :: Da
    
    ! name of the matrix
    character(len=*), intent(in) :: sarray
    
    ! suppress zeroes in output: yes/no
    logical, intent(in) :: bnoZero
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile
    
    ! format string to use for the output; e.g. '(D20.10)'
    character(len=*), intent(in) :: sformat
    
    ! Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    real(DP), intent(in) :: dthreshold

  !</input>
    
!</subroutine>
    
    !local variables
    integer :: i, j, cf, nchar
    real(DP) :: dval
    character(len=128) :: S
    character(len=6) :: sformatChar
    
    if (ifile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE)
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'matio_writeMatrix1_Dble')
        call sys_halt()
      end if
    else
      cf = ifile
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
        if ((.not. bnoZero) .or. (abs(dval) .gt. dthreshold)) then
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
    if (ifile .eq. 0) close(cf)
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  subroutine matio_writeMatrix79_Dble (Da, Icol, Irow, &
                                       nrow, ncol, sarray, &
                                       bnoZero, ifile, sfile, sformat,dthreshold)
  
  !<description>
    ! Write sparse double precision matrix in matrix format 9 or
    ! matrix format 9 into a text file.
    !
    ! Double-precision version
  !</description>
    
  !<input>
    ! number of rows
    integer, intent(in) :: nrow
    
    ! number of columns; must be =nrow for structure-7 matrices
    integer, intent(in) :: ncol
    
    ! matrix: array [1..na] of double
    real(DP), dimension(:), intent(in) :: Da
    
    ! Column structure of the matrix
    integer, dimension(:), intent(in) :: Icol
    
    ! Row structure of the matrix
    integer, dimension(nrow+1), intent(in) :: Irow
    
    ! name of the matrix
    character(len=*), intent(in) :: sarray
    
    ! suppress zeroes in output: yes/no
    logical, intent(in) :: bnoZero
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile
    
    ! format string to use for the output; e.g. '(D20.10)'
    character(len=*), intent(in) :: sformat
    
    ! Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    real(DP), intent(in) :: dthreshold

  !</input>
    
!</subroutine>
    
    !local variables
    integer :: i, j, k, cf, nchar
    real(DP) :: dval
    character(len=128) :: S
    character(len=6) :: sformatChar
    integer :: h_DrowVec
    real(DP), dimension(:), pointer :: p_DrowVec
    
    if (ifile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE)
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'matio_writeMatrix79_Dble')
        call sys_halt()
      end if
    else
      cf = ifile
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
    call storage_new('matio_writeMatrix79_Dble','DrowVec',ncol, &
                     ST_DOUBLE, h_DrowVec, ST_NEWBLOCK_NOINIT)
    call storage_getbase_double(h_DrowVec, p_DrowVec)
    do i=1, nrow
      
      ! Extract row i
      if (bnoZero) then
        ! SYS_MAXREAL_DP is written out as '.'
        do j=1, ncol
          p_DrowVec(j) = SYS_MAXREAL_DP
        end do
      else
        do j=1, ncol
          p_DrowVec(j) = 0.0_DP
        end do
      end if
      
      do j=0, Irow(i+1)-Irow(i)-1
        k = Irow(i)+j
        p_DrowVec(Icol(k)) = Da(k)
      end do

      ! Write row i
      do j=1, ncol-1
        dval = p_DrowVec(j)
        if (abs(dval) .lt. dthreshold) dval = 0.0_DP
        if (bnoZero .and. (dval .eq. SYS_MAXREAL_DP)) then
          write (cf,sformatChar, ADVANCE='NO') '.'
        else
          write (cf,sformat,ADVANCE='NO') dval
        end if
      end do
      
      dval = p_DrowVec(ncol)
      if (abs(dval) .lt. dthreshold) dval = 0.0_DP
      if (bnoZero .and. (dval .eq. SYS_MAXREAL_DP)) then
        write (cf,sformatChar, ADVANCE='YES') '.'
      else
        write (cf,sformat,ADVANCE='YES') dval
      end if
    end do
    
    call storage_free(h_DrowVec)
    
    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine matio_writeMatrix1_Sngl (Fa, nrow, ncol, sarray, &
                                      bnoZero, ifile, sfile, sformat, fthreshold)
  
  !<description>
    ! Write full single precision matrix into a text file.
    !
    ! This writes an array Da with nrow rows and ncol columns as a matrix
    ! into a text file.
  !</description>
    
  !<input>
    ! number of rows
    integer, intent(in) :: nrow
    
    ! number of columns
    integer, intent(in) :: ncol
    
    ! matrix: array [1..nrow,1..ncol] of single
    real(SP), dimension(nrow,ncol), intent(in) :: Fa
    
    ! name of the matrix
    character(len=*), intent(in) :: sarray
    
    ! suppress zeroes in output: yes/no
    logical, intent(in) :: bnoZero
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile
    
    ! format string to use for the output; e.g. '(F20.10)'
    character(len=*), intent(in) :: sformat
    
    ! Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    real(SP), intent(in) :: fthreshold

  !</input>
    
!</subroutine>
    
    !local variables
    integer :: i, j, cf, nchar
    real(SP) :: fval
    character(len=128) :: S
    character(len=6) :: sformatChar
    
    if (ifile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE)
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'matio_writeMatrix1_Sngl')
        call sys_halt()
      end if
    else
      cf = ifile
    end if
    
    ! Get length of output strings
    S(:) = ' '
    write (S,sformat) 0.0_SP
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
        fval = Fa(i,j)
        if ((.not. bnoZero) .or. (abs(fval) .gt. fthreshold)) then
          write (cf,sformat,ADVANCE='NO') fval
        else
          write (cf,sformatChar, ADVANCE='NO') '.'
        end if
      end do
      fval = Fa(i,j)
      if ((.not. bnoZero) .or. (fval .ne. 0.0_SP)) then
        write (cf,sformat) fval
      else
        write (cf,sformatChar) '.'
      end if
    end do
    
    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  subroutine matio_writeMatrix79_Sngl (Fa, Icol, Irow, &
                                       nrow, ncol, sarray, &
                                       bnoZero, ifile, sfile, sformat, fthreshold)
  
  !<description>
    ! Write sparse single precision matrix in matrix format 9 or
    ! matrix format 9 into a text file.
    !
    ! Double-precision version
  !</description>
    
  !<input>
    ! number of rows
    integer, intent(in) :: nrow
    
    ! number of columns; must be =nrow for structure-7 matrices
    integer, intent(in) :: ncol
    
    ! matrix: array [1..na] of double
    real(SP), dimension(:), intent(in) :: Fa
    
    ! Column structure of the matrix
    integer, dimension(:), intent(in) :: Icol
    
    ! Row structure of the matrix
    integer, dimension(nrow+1), intent(in) :: Irow
    
    ! name of the matrix
    character(len=*), intent(in) :: sarray
    
    ! suppress zeroes in output: yes/no
    logical, intent(in) :: bnoZero
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile
    
    ! format string to use for the output; e.g. '(D20.10)'
    character(len=*), intent(in) :: sformat
    
    ! Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    real(SP), intent(in) :: fthreshold

  !</input>
    
!</subroutine>
    
    !local variables
    integer :: i, j, k, cf, nchar
    real(SP) :: fval
    character(len=128) :: S
    character(len=6) :: sformatChar
    integer :: h_FrowVec
    real(SP), dimension(:), pointer :: p_FrowVec
    
    if (ifile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE)
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'matio_writeMatrix79_Sngl')
        call sys_halt()
      end if
    else
      cf = ifile
    end if
    
    ! Get length of output strings
    S(:) = ' '
    write (S,sformat) 0.0_SP
    nchar = len(trim(S))
    
    ! Build array format string
    sformatChar = '(A'//sys_i3(nchar)//')'
    
    ! Write all format strings into the file
    write (cf,'(3A15,L2,3I15)') sarray, sformat, sformatChar, &
                                bnoZero, nchar, nrow, ncol

    if (nrow .le. 0) return
    if (ncol .le. 0) return
    
    ! Write the matrix
    call storage_new('matio_writeMatrix79_Sngl','FrowVec',ncol, &
                     ST_SINGLE, h_FrowVec, ST_NEWBLOCK_NOINIT)
    call storage_getbase_single(h_FrowVec, p_FrowVec)
    do i=1, nrow
      
      ! Extract row i
      if (bnoZero) then
        ! SYS_MAXREAL_SP is written out as '.'
        do j=1, ncol
          p_FrowVec(j) = SYS_MAXREAL_SP
        end do
      else
        do j=1, ncol
          p_FrowVec(j) = 0.0_SP
        end do
      end if
      
      do j=0, Irow(i+1)-Irow(i)-1
        k = Irow(i)+j
        p_FrowVec(Icol(k)) = Fa(k)
      end do

      ! Write row i
      do j=1, ncol-1
        fval = p_FrowVec(j)
        if (abs(fval) .lt. fthreshold) fval = 0.0_SP
        if (bnoZero .and. (fval .eq. SYS_MAXREAL_SP)) then
          write (cf,sformatChar, ADVANCE='NO') '.'
        else
          write (cf,sformat,ADVANCE='NO') fval
        end if
      end do
      
      fval = p_FrowVec(ncol)
      if (abs(fval) .lt. fthreshold) fval = 0.0_SP
      if (bnoZero .and. (fval .eq. SYS_MAXREAL_SP)) then
        write (cf,sformatChar, ADVANCE='YES') '.'
      else
        write (cf,sformat,ADVANCE='YES') fval
      end if
    end do
    
    call storage_free(h_FrowVec)
    
    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)
  
  end subroutine


  ! ***************************************************************************

!<subroutine>
  subroutine matio_writeMatrixMaple (rmatrix, sarray,&
                                     ifile, sfile, sformat)
  
  !<description>
    ! This routine writes a scalar matrix into a text file using the MAPLE
    ! syntax.
  !</description>
    
  !<input>
    ! The matrix to be written out
    type(t_matrixScalar), intent(in) :: rmatrix
    
    ! Name of the matrix
    character(len=*), intent(in) :: sarray
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile
    
    ! Format string to use for the output; e.g. '(E20.10)'
    character(len=*), intent(in) :: sformat
    
  !</input>
    
!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_DA
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kld

  ! Depending on the matrix format, choose the right routine for writing
  select case (rmatrix%cmatrixFormat)
  case (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
    ! Matrix precision?
    select case (rmatrix%cdataType)
    case (ST_DOUBLE)
      ! Get the data arrays and write the matrix
      if (.not.lsyssc_hasMatrixContent(rmatrix)) then
        call output_line('Matrix has no data',&
            OU_CLASS_ERROR,OU_MODE_STD,'matio_writeMapleMatrix')
        call sys_halt()
      end if
      call lsyssc_getbase_double (rmatrix,p_Da)
      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)
      call matio_writeMapleMatrix79_D (p_Da, p_Kcol, p_Kld, &
                                       rmatrix%NEQ, rmatrix%NEQ, sarray, &
                                       ifile, sfile, sformat)
    case DEFAULT
      call output_line ('Unsupported matrix precision!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'matio_writeMapleMatrix')
      call sys_halt()
    end select
  case DEFAULT
    call output_line ('Unknown matrix format!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'matio_writeMapleMatrix')
    call sys_halt()
  end select
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
  subroutine matio_writeMapleMatrix79_D (Da, Icol, Irow, &
                                        nrow, ncol, sarray, &
                                        ifile, sfile, sformat)
  
  !<description>
    ! Write sparse double precision matrix in matrix format 9 or
    ! matrix format 9 into a text file using the Maple syntax.
    !
    ! Double-precision version
  !</description>
    
  !<input>
    ! number of rows
    integer, intent(in) :: nrow
    
    ! number of columns; must be =nrow for structure-7 matrices
    integer, intent(in) :: ncol
    
    ! matrix: array [1..na] of double
    real(DP), dimension(:), intent(in) :: Da
    
    ! Column structure of the matrix
    integer, dimension(:), intent(in) :: Icol
    
    ! Row structure of the matrix
    integer, dimension(nrow+1), intent(in) :: Irow
    
    ! name of the matrix
    character(len=*), intent(in) :: sarray
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile
    
    ! format string to use for the output; e.g. '(D20.10)'
    character(len=*), intent(in) :: sformat
  !</input>
    
!</subroutine>
    
    !local variables
    integer :: i, j, cf, nchar
    character(len=32) :: S
    character(len=6) :: sformatChar
    
    if (ifile .eq. 0) then
      call io_openFileForWriting(sfile, cf, SYS_REPLACE)
      if (cf .eq. -1) then
        call output_line ('Could not open file '//trim(sfile), &
                          OU_CLASS_ERROR,OU_MODE_STD,'matio_writeMapleMatrix79_D')
        call sys_halt()
      end if
    else
      cf = ifile
    end if
    
    ! Get length of output strings
    S(:) = ' '
    write (S,sformat) 0.0_DP
    nchar = len(trim(S))
    
    ! Build array format string
    sformatChar = '(A'//sys_i3(nchar)//')'
    
    if (nrow .le. 0) return
    if (ncol .le. 0) return
    
    ! Write a header to the file that declares the matrix.
    write (cf,'(6A)') sarray,' := matrix(',&
        trim(sys_siL(nrow,10)),',',trim(sys_siL(ncol,10)),',0):'
        
    ! Now the entries. This is a sparse matrix, so we insert commands
    ! only for the entries.
    do i=1, nrow
    
      do j=Irow(i),Irow(i+1)-1
        write (s,sformat) Da(j)
        write (cf,'(A)') &
            sarray//'['//trim(sys_siL(i,10))//','//trim(sys_siL(Icol(j),10))//']:='//&
            trim(adjustl(s))//':'
      end do
      
    end do
    
    ! Close the file if necessary
    if (ifile .eq. 0) close(cf)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine matio_writeBlockMatrixMaple (rmatrix, sarray,&
                                          ifile, sfile, sformat,dthreshold)
  
  !<description>
    ! This routine writes a block matrix into a text file using the MAPLE
    ! syntax.
    ! Note that for this purpose, a new matrix is temporarily created in memory!
  !</description>
    
  !<input>
    ! The matrix to be written out
    type(t_matrixBlock), intent(in) :: rmatrix
    
    ! Name of the matrix
    character(len=*), intent(in) :: sarray
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Do not close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(in) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(in) :: sfile
    
    ! Format string to use for the output; e.g. '(E20.10)'
    character(len=*), intent(in) :: sformat
    
    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for beter visualisation.
    ! If not present, a default of 1E-12 is assumed.
    real(DP), intent(in), optional :: dthreshold
  !</input>
    
!</subroutine>

    ! local variables
    type(t_matrixBlock) :: rtempMatrix
    real(DP), dimension(:), pointer :: p_DA
    real(DP) :: dthres

    ! We have to create a global matrix first!
    call glsys_assembleGlobal (rmatrix,rtempMatrix,.true.,.true.)
                              
    ! Replace small values by zero
    dthres = 1E-12_DP
    if (present(dthreshold)) dthres = dthreshold
    if (abs(dthres) .gt. 0.0_DP) then
      call lsyssc_getbase_double (rtempMatrix%RmatrixBlock(1,1),p_DA)
      where (abs(p_Da) .lt. dthres) p_Da = 0.0_DP
    end if
    
    ! Write matrix to the file
    call matio_writeMatrixMaple (rtempMatrix%RmatrixBlock(1,1), sarray,&
                                     ifile, sfile, sformat)

    ! Release the temporary matrix
    call lsysbl_releaseMatrix (rtempMatrix)

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine matio_spyBlockMatrix(sfilename,smatrixName,rmatrix,bdata,cstatus,dthreshold)

!<description>
    ! Writes a block matrix in CSR format to a file so that its sparsity structure
    ! can be visualised in MATLAB by means of the spy command.
    ! If the matrix contains entries, the entries are written to the file
    ! as well, so the matrix can be used as sparse matrix for arbitrary purposes.
    !
    ! To load a matrix written out by this routine, one has simply to type the name
    ! of the ".m"-file that is written out. MATLAB will read the file and
    ! create a sparse matrix with the name smatrixName in memory.
!</description>

!<input>
    ! File name of the MATLAB file without fileextension. A ".m" is appended.
    character(LEN=*), intent(in) :: sfileName
    
    ! Name of the matrix in MATLAB file. This will be the name of the
    ! variable containing the matrix data when reading the file into matlab.
    character(LEN=*), intent(in) :: smatrixName

    ! Source matrix
    type(t_matrixBlock), intent(in) :: rmatrix
    
    ! Whether to spy the real data of the matrix or only its sparsity
    ! pattern
    logical, intent(in) :: bdata

    ! OPTIONAL: status of file
    integer, intent(in), optional :: cstatus

    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for beter visualisation.
    ! If not present, a default of 1E-12 is assumed.
    real(DP), intent(in), optional :: dthreshold
!</input>
!</subroutine>

    ! local variables
    type(t_matrixBlock) :: rtempMatrix

    ! We have to create a global matrix first!
    call glsys_assembleGlobal (rmatrix,rtempMatrix,.true.,.true.)
                              
    ! Write matrix to the file
    call matio_spyMatrix(sfilename,smatrixName,rtempMatrix%RmatrixBlock(1,1),&
        bdata,cstatus,dthreshold)

    ! Release the temporary matrix
    call lsysbl_releaseMatrix (rtempMatrix)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine matio_spyMatrix(sfilename,smatrixName,rmatrix,bdata,cstatus,dthreshold)

!<description>
    ! Writes a scalar matrix in CSR format to a file so that its sparsity structure
    ! can be visualised in MATLAB by means of the spy command.
    ! If the matrix contains entries, the entries are written to the file
    ! as well, so the matrix can be used as sparse matrix for arbitrary purposes.
    !
    ! To load a matrix written out by this routine, one has simply to type the name
    ! of the ".m"-file that is written out. MATLAB will read the file and
    ! create a sparse matrix with the name smatrixName in memory.
!</description>

!<input>
    ! File name of the MATLAB file without fileextension. A ".m" is appended.
    character(LEN=*), intent(in) :: sfileName
    
    ! Name of the matrix in MATLAB file. This will be the name of the
    ! variable containing the matrix data when reading the file into matlab.
    character(LEN=*), intent(in) :: smatrixName

    ! Source matrix
    type(t_matrixScalar), intent(in) :: rmatrix
    
    ! Whether to spy the real data of the matrix or only its sparsity
    ! pattern
    logical, intent(in) :: bdata

    ! OPTIONAL: status of file
    integer, intent(in), optional :: cstatus

    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    ! If not present, a default of 1E-12 is assumed.
    real(DP), intent(in), optional :: dthreshold
!</input>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), pointer :: p_Da
    real(SP), dimension(:), pointer :: p_Fa
    integer, dimension(:), pointer :: p_Kld,p_Kcol
    integer :: iunit,ieq
    real(DP) :: dthres
    character(LEN=10) :: cstat,cpos

    ! Replace small values by zero
    dthres = 1E-12_DP
    if (present(dthreshold)) dthres = dthreshold

    ! Set file status if required
    if (present(cstatus)) then
      select case(cstatus)
      case (IO_NEW)
        cstat="NEW"; cpos="ASIS"
      case (IO_REPLACE)
        cstat="REPLACE"; cpos="ASIS"
      case (IO_OLD)
        cstat="OLD"; cpos="APPEND"
      case DEFAULT
        cstat="UNKNOWN"; cpos ="ASIS"
      end select
    else
      cstat="UNKNOWN"; cpos="ASIS"
    end if
    
    ! Open output file
    iunit=sys_getFreeUnit()
    open (UNIT=iunit,STATUS=trim(cstat),POSITION=trim(cpos),FILE=trim(adjustl(sfilename))//'.m')
    
    ! Which matrix format are we?
    select case(rmatrix%cmatrixFormat)
      
    case (LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9INTL)
      
      call lsyssc_getbase_Kld(rmatrix,p_Kld)
      call lsyssc_getbase_Kcol(rmatrix,p_Kcol)
      
      if (bdata) then
        ! Which matrix type are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_double(rmatrix,p_Da)
          select case(rmatrix%cinterleavematrixFormat)
          case (LSYSSC_MATRIXD)
            call do_spy_mat79matD_double(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,&
                p_Kld,p_Kcol,p_Da,dthres)
          case (LSYSSC_MATRIX1)
            call do_spy_mat79mat1_double(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,&
                rmatrix%NVAR,p_Kld,p_Kcol,p_Da,dthres)
          case DEFAULT
            call output_line ('Unsupported interleave matrix type!', &
                              OU_CLASS_ERROR,OU_MODE_STD,'matio_spyMatrix')
            call sys_halt()
          end select
          write(UNIT=iunit,FMT=30)
          
        case (ST_SINGLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_single(rmatrix,p_Fa)
          select case(rmatrix%cinterleavematrixFormat)
          case (LSYSSC_MATRIXD)
            call do_spy_mat79matD_single(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,&
                p_Kld,p_Kcol,p_Fa,dthres)
          case (LSYSSC_MATRIX1)
            call do_spy_mat79mat1_single(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,&
                rmatrix%NVAR,p_Kld,p_Kcol,p_Fa,dthres)
          case DEFAULT
            call output_line ('Unsupported interleave matrix type!', &
                              OU_CLASS_ERROR,OU_MODE_STD,'matio_spyMatrix')
            call sys_halt()
          end select
          write(UNIT=iunit,FMT=30)

        case DEFAULT
          call output_line ('Unsupported matrix type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'matio_spyMatrix')
          call sys_halt()
        end select
        
      else
        
        ! Output only matrix structure
        write(UNIT=iunit,FMT=10)
        select case(rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIXD)
          call do_spy_mat79matD_double(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,&
              p_Kld,p_Kcol)
        case (LSYSSC_MATRIX1)
          call do_spy_mat79mat1_double(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,&
              rmatrix%NVAR,p_Kld,p_Kcol)
        case DEFAULT
          call output_line ('Unsupported interleave matrix type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'matio_spyMatrix')
          call sys_halt()
        end select
        write(UNIT=iunit,FMT=30)
      end if

    case (LSYSSC_MATRIX7,LSYSSC_MATRIX9)
      
      call lsyssc_getbase_Kld(rmatrix,p_Kld)
      call lsyssc_getbase_Kcol(rmatrix,p_Kcol)
      
      if (bdata) then
        ! Which matrix type are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_double(rmatrix,p_Da)
          call do_spy_mat79matD_double(rmatrix%NEQ,rmatrix%NCOLS,1,p_Kld,p_Kcol,p_Da,dthres)
          write(UNIT=iunit,FMT=30)

        case (ST_SINGLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_single(rmatrix,p_Fa)
          call do_spy_mat79matD_single(rmatrix%NEQ,rmatrix%NCOLS,1,p_Kld,p_Kcol,p_Fa,dthres)
          write(UNIT=iunit,FMT=30)
          
        case DEFAULT
          call output_line ('Unsupported matrix type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'matio_spyMatrix')
          call sys_halt()
        end select

      else
        
        ! Output only matrix structure
        write(UNIT=iunit,FMT=10)
        call do_spy_mat79matD_double(rmatrix%NEQ,rmatrix%NCOLS,1,p_Kld,p_Kcol)
        write(UNIT=iunit,FMT=30)
      end if
      
    case(LSYSSC_MATRIXD)
      
      if (bdata) then
        ! Which matrix type are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_double(rmatrix,p_Da)
          do ieq=1,rmatrix%NEQ
            if (abs(p_Da(ieq)) .ge. dthres) then
              write(UNIT=iunit,FMT=20) ieq,ieq,p_Da(ieq)
            end if
          end do
          write(UNIT=iunit,FMT=30)

        case (ST_SINGLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_single(rmatrix,p_Fa)
          do ieq=1,rmatrix%NEQ
            if (abs(p_Fa(ieq)) .ge. dthres) then
              write(UNIT=iunit,FMT=20) ieq,ieq,p_Fa(ieq)
            end if
          end do
          write(UNIT=iunit,FMT=30)

        case DEFAULT
          call output_line ('Unsupported matrix type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'matio_spyMatrix')
          call sys_halt()
        end select
        
      else
        
        ! Output only matrix structure
        write(UNIT=iunit,FMT=10)
        do ieq=1,rmatrix%NEQ
          write(UNIT=iunit,FMT=20) ieq,ieq,1
        end do
        write(UNIT=iunit,FMT=30)
        
      end if
      
    case(LSYSSC_MATRIX1)

      if (bdata) then
        ! Which matrix type are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_double(rmatrix,p_Da)
          call do_spy_mat1_double(rmatrix%NEQ,rmatrix%NCOLS,p_Da,dthres)
          write(UNIT=iunit,FMT=30)

        case (ST_SINGLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_single(rmatrix,p_Fa)
          call do_spy_mat1_single(rmatrix%NEQ,rmatrix%NCOLS,p_Fa,dthres)
          write(UNIT=iunit,FMT=30)

        case DEFAULT
          call output_line ('Unsupported matrix type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'matio_spyMatrix')
          call sys_halt()
        end select
        
      else
        
        ! Output only matrix structure
        write(UNIT=iunit,FMT=10)
        call do_spy_mat1_double(rmatrix%NEQ,rmatrix%NCOLS)
        write(UNIT=iunit,FMT=30)

      end if
      
    case DEFAULT
      call output_line ('Unsupported matrix format!', &
                         OU_CLASS_ERROR,OU_MODE_STD,'matio_spyMatrix')
      call sys_halt()
      
    end select
    
    ! Close file
    write(UNIT=iunit,FMT=40) smatrixName, rmatrix%NEQ, rmatrix%NCOLS
    close(UNIT=iunit)

10  format("data=[...")
20  format(I10,1X,I10,1X,E15.8,";")
30  format("];")
40  format(A,"=sparse(data(:,1),data(:,2),data(:,3),",I10,",",I10,"); clear data;")

  contains

    ! Here, the real SPY routines follow.

    !**************************************************************
    ! SPY CSR matrix in double precision
    
    subroutine do_spy_mat79matD_double(neq,ncols,nvar,Kld,Kcol,Da,dthres)
      integer, dimension(:), intent(in)  :: Kld
      integer, dimension(:), intent(in)  :: Kcol
      integer, intent(in)                :: neq,ncols
      integer, intent(in)                             :: nvar
      real(DP), dimension(nvar,*), intent(in), optional :: Da
      real(DP), intent(in), optional :: dthres
      real(DP) :: ddata
      integer :: ieq,ild,ivar
      
      if (present(Da)) then
        do ieq=1,neq
          do ild=Kld(ieq),Kld(ieq+1)-1
            do ivar=1,nvar
              ddata = Da(ivar,ild)
              if (abs(ddata) .ge. dthres) then
                write(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                    (ivar-1)*neq+ieq,(ivar-1)*ncols+Kcol(ild),ddata
              end if
            end do
          end do
        end do
      else
        do ieq=1,neq
          do ild=Kld(ieq),Kld(ieq+1)-1
            do ivar=1,nvar
              write(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                  (ivar-1)*neq+ieq,(ivar-1)*ncols+Kcol(ild)
            end do
          end do
        end do
      end if
    end subroutine do_spy_mat79matD_double

    !**************************************************************
    ! SPY CSR matrix in single precision
    
    subroutine do_spy_mat79matD_single(neq,ncols,nvar,Kld,Kcol,Fa,dthres)
      integer, dimension(:), intent(in)  :: Kld
      integer, dimension(:), intent(in)  :: Kcol
      integer, intent(in)                :: neq,ncols
      integer, intent(in)                             :: nvar
      real(SP), dimension(nvar,*), intent(in), optional :: Fa
      real(DP), intent(in), optional :: dthres
      real(SP) :: fdata
      integer :: ieq,ild,ivar

      if (present(Fa)) then
        do ieq=1,neq
          do ild=Kld(ieq),Kld(ieq+1)-1
            do ivar=1,nvar
              fdata = Fa(ivar,ild)
              if (abs(fdata) .ge. dthres) then
                write(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                    (ivar-1)*neq+ieq,(ivar-1)*ncols+Kcol(ild),fdata
              end if
            end do
          end do
        end do
      else
        do ieq=1,neq
          do ild=Kld(ieq),Kld(ieq+1)-1
            do ivar=1,nvar
              write(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                  (ivar-1)*neq+ieq,(ivar-1)*ncols+Kcol(ild)
            end do
          end do
        end do
      end if
    end subroutine do_spy_mat79matD_single

    !**************************************************************
    ! SPY CSR matrix in double precision
    
    subroutine do_spy_mat79mat1_double(neq,ncols,nvar,mvar,Kld,Kcol,Da,dthres)
      integer, dimension(:), intent(in)  :: Kld
      integer, dimension(:), intent(in)  :: Kcol
      integer, intent(in)                :: neq,ncols
      integer, intent(in)                             :: nvar,mvar
      real(DP), dimension(nvar,mvar,*), intent(in), optional :: Da
      real(DP), intent(in), optional :: dthres
      real(DP) :: ddata
      integer :: ieq,ild,ivar,jvar

      if (present(Da)) then
        do ieq=1,neq
          do ild=Kld(ieq),Kld(ieq+1)-1
            do jvar=1,nvar ! local column index
              do ivar=1,mvar ! local row index
                ddata = Da(ivar,jvar,ild)
                if (abs(ddata) .ge. dthres) then
                  write(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                      (ivar-1)*neq+ieq,(jvar-1)*ncols+Kcol(ild),ddata
                end if
              end do
            end do
          end do
        end do
      else
        do ieq=1,neq
          do ild=Kld(ieq),Kld(ieq+1)-1
            do jvar=1,nvar ! local column index
              do ivar=1,mvar ! local row index
                write(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                    (ivar-1)*neq+ieq,(jvar-1)*ncols+Kcol(ild)
              end do
            end do
          end do
        end do
      end if
    end subroutine do_spy_mat79mat1_double

    !**************************************************************
    ! SPY CSR matrix in double precision
    
    subroutine do_spy_mat79mat1_single(neq,ncols,nvar,mvar,Kld,Kcol,Fa,dthres)
      integer, dimension(:), intent(in)  :: Kld
      integer, dimension(:), intent(in)  :: Kcol
      integer, intent(in)                :: neq,ncols
      integer, intent(in)                             :: nvar,mvar
      real(SP), dimension(nvar,mvar,*), intent(in), optional :: Fa
      real(DP), intent(in), optional :: dthres
      real(SP) :: fdata
      integer :: ieq,ild,ivar,jvar

      if (present(Fa)) then
        do ieq=1,neq
          do ild=Kld(ieq),Kld(ieq+1)-1
            do jvar=1,nvar ! local column index
              do ivar=1,mvar ! local row index
                fdata = Fa(ivar,jvar,ild)
                if (abs(fdata) .ge. dthres) then
                  write(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                      (ivar-1)*neq+ieq,(jvar-1)*ncols+Kcol(ild),fdata
                end if
              end do
            end do
          end do
        end do
      else
        do ieq=1,neq
          do ild=Kld(ieq),Kld(ieq+1)-1
            do jvar=1,nvar ! local column index
              do ivar=1,mvar ! local row index
                write(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                    (ivar-1)*neq+ieq,(jvar-1)*ncols+Kcol(ild)
              end do
            end do
          end do
        end do
      end if
    end subroutine do_spy_mat79mat1_single

    !**************************************************************
    ! SPY full matrix in double precision
    
    subroutine do_spy_mat1_double(neq,ncols,Da,dthres)
      integer, intent(in) :: neq,ncols
      real(DP), dimension(:), intent(in), optional :: Da
      real(DP), intent(in), optional :: dthres
      real(DP) :: ddata
      integer :: ieq,icol

      if (present(Da)) then
        do ieq=1,neq
          do icol=1,ncols
            ddata = Da((icol-1)*neq+ieq)
            if (abs(ddata) .ge. dthres) then
              write(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                  ieq,icol,ddata
            end if
          end do
        end do
      else
        do ieq=1,neq
          do icol=1,ncols
            write(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                ieq,icol,Da((icol-1)*neq+ieq)
          end do
        end do
      end if
    end subroutine do_spy_mat1_double

    !**************************************************************
    ! SPY full matrix in single precision
    
    subroutine do_spy_mat1_single(neq,ncols,Fa,dthres)
      integer, intent(in) :: neq,ncols
      real(SP), dimension(:), intent(in), optional :: Fa
      real(DP), intent(in), optional :: dthres
      real(SP) :: fdata
      integer :: ieq,icol

      if (present(Fa)) then
        do ieq=1,neq
          do icol=1,ncols
            fdata = Fa((icol-1)*neq+ieq)
            if (abs(fdata) .ge. dthres) then
              write(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                  ieq,icol,fdata
            end if
          end do
        end do
      else
        do ieq=1,neq
          do icol=1,ncols
            write(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                ieq,icol,Fa((icol-1)*neq+ieq)
          end do
        end do
      end if
    end subroutine do_spy_mat1_single
  end subroutine matio_spyMatrix

end module

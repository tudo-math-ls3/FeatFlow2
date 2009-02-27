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
!#        visualized by means of the MATLAB command SPY
!#
!# 4.) matio_spyBlockMatrix
!#     -> Writes a block matrix in CSR format into a file which can be 
!#        visualized by means of the MATLAB command SPY
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
!# </purpose>
!#########################################################################

module matrixio

  use fsystem
  use storage
  use io
  use linearsystemscalar
  use linearsystemblock
  use globalsystem
  
  implicit none

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
    type(t_matrixBlock), intent(IN) :: rmatrix
    
    ! Name of the matrix
    character(len=*), intent(IN) :: sarray
    
    ! Suppress zeroes in output: yes/no
    logical, intent(IN) :: bnoZero
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! Format string to use for the output; e.g. '(E20.10)'
    character(len=*), intent(IN) :: sformat
    
    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    ! If not present, values are not replaced (i.e. the default value is 0.0).
    real(DP), intent(IN), optional :: dthreshold
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
    type(t_matrixScalar), intent(IN) :: rmatrix
    
    ! Name of the matrix
    character(len=*), intent(IN) :: sarray
    
    ! Suppress zeroes in output: yes/no
    logical, intent(IN) :: bnoZero
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! Format string to use for the output; e.g. '(E20.10)'
    character(len=*), intent(IN) :: sformat
    
    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    ! If not present, values are not replaced (i.e. the default value is 0.0).
    real(DP), intent(IN), optional :: dthreshold

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
    integer, intent(IN) :: nrow
    
    ! number of columns
    integer, intent(IN) :: ncol
    
    ! matrix: array [1..nrow,1..ncol] of double
    real(DP), dimension(nrow,ncol), intent(IN) :: Da
    
    ! name of the matrix
    character(len=*), intent(IN) :: sarray
    
    ! suppress zeroes in output: yes/no
    logical, intent(IN) :: bnoZero
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! format string to use for the output; e.g. '(D20.10)'
    character(len=*), intent(IN) :: sformat
    
    ! Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    real(DP), intent(IN) :: dthreshold

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
                          OU_CLASS_ERROR,OU_MODE_STD,'matio_writeFullMatrix')
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
    integer, intent(IN) :: nrow
    
    ! number of columns; must be =nrow for structure-7 matrices
    integer, intent(IN) :: ncol
    
    ! matrix: array [1..na] of double
    real(DP), dimension(:), intent(IN) :: Da
    
    ! Column structure of the matrix
    integer, dimension(:), intent(IN) :: Icol
    
    ! Row structure of the matrix
    integer, dimension(nrow+1), intent(IN) :: Irow
    
    ! name of the matrix
    character(len=*), intent(IN) :: sarray
    
    ! suppress zeroes in output: yes/no
    logical, intent(IN) :: bnoZero
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! format string to use for the output; e.g. '(D20.10)'
    character(len=*), intent(IN) :: sformat
    
    ! Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for better visualisation.
    real(DP), intent(IN) :: dthreshold

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
                          OU_CLASS_ERROR,OU_MODE_STD,'matio_writeFullMatrix')
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
    call storage_new('matio_writeSparseMatrix','DrowVec',ncol, &
                     ST_DOUBLE, h_DrowVec, ST_NEWBLOCK_NOINIT)
    call storage_getbase_double(h_DrowVec, p_DrowVec)
    do i=1, nrow
      
      ! Extract row i
      if (bnoZero) then
        ! SYS_MAXREAL is written out as '.'
        do j=1, ncol
          p_DrowVec(j) = SYS_MAXREAL
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
        if (bnoZero .and. (dval .eq. SYS_MAXREAL)) then
          write (cf,sformatChar, ADVANCE='NO') '.'
        else
          write (cf,sformat,ADVANCE='NO') dval
        end if
      end do
      
      dval = p_DrowVec(ncol)
      if (abs(dval) .lt. dthreshold) dval = 0.0_DP
      if (bnoZero .and. (dval .eq. SYS_MAXREAL)) then
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
  subroutine matio_writeMatrixMaple (rmatrix, sarray,&
                                     ifile, sfile, sformat)
  
  !<description>
    ! This routine writes a scalar matrix into a text file using the MAPLE
    ! syntax.
  !</description>
    
  !<input>
    ! The matrix to be written out
    type(t_matrixScalar), intent(IN) :: rmatrix
    
    ! Name of the matrix
    character(len=*), intent(IN) :: sarray
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! Format string to use for the output; e.g. '(E20.10)'
    character(len=*), intent(IN) :: sformat
    
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
    integer, intent(IN) :: nrow
    
    ! number of columns; must be =nrow for structure-7 matrices
    integer, intent(IN) :: ncol
    
    ! matrix: array [1..na] of double
    real(DP), dimension(:), intent(IN) :: Da
    
    ! Column structure of the matrix
    integer, dimension(:), intent(IN) :: Icol
    
    ! Row structure of the matrix
    integer, dimension(nrow+1), intent(IN) :: Irow
    
    ! name of the matrix
    character(len=*), intent(IN) :: sarray
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! format string to use for the output; e.g. '(D20.10)'
    character(len=*), intent(IN) :: sformat
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
    type(t_matrixBlock), intent(IN) :: rmatrix
    
    ! Name of the matrix
    character(len=*), intent(IN) :: sarray
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    integer, intent(IN) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    character(len=*), intent(IN) :: sfile
    
    ! Format string to use for the output; e.g. '(E20.10)'
    character(len=*), intent(IN) :: sformat
    
    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for beter visualisation.
    ! If not present, a default of 1E-12 is assumed.
    real(DP), intent(IN), optional :: dthreshold
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
    ! can be visualized in MATLAB by means of the spy command.
    ! If the matrix contains entries, the entries are written to the file
    ! as well, so the matrix can be used as sparse matrix for arbitrary purposes.
    !
    ! To load a matrix written out by this routine, one has simply to type the name
    ! of the ".m"-file that is written out. MATLAB will read the file and
    ! create a sparse matrix with the name smatrixName in memory.
!</description>

!<input>
    ! File name of the MATLAB file without fileextension. A ".m" is appended.
    character(LEN=*), intent(IN) :: sfileName
    
    ! Name of the matrix in MATLAB file. This will be the name of the
    ! variable containing the matrix data when reading the file into matlab.
    character(LEN=*), intent(IN) :: smatrixName

    ! Source matrix
    type(t_matrixBlock), intent(IN) :: rmatrix
    
    ! Whether to spy the real data of the matrix or only its sparsity
    ! pattern
    logical, intent(IN) :: bdata

    ! OPTIONAL: status of file
    integer, intent(IN), optional :: cstatus

    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for beter visualisation.
    ! If not present, a default of 1E-12 is assumed.
    real(DP), intent(IN), optional :: dthreshold
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
    ! can be visualized in MATLAB by means of the spy command.
    ! If the matrix contains entries, the entries are written to the file
    ! as well, so the matrix can be used as sparse matrix for arbitrary purposes.
    !
    ! To load a matrix written out by this routine, one has simply to type the name
    ! of the ".m"-file that is written out. MATLAB will read the file and
    ! create a sparse matrix with the name smatrixName in memory.
!</description>

!<input>
    ! File name of the MATLAB file without fileextension. A ".m" is appended.
    character(LEN=*), intent(IN) :: sfileName
    
    ! Name of the matrix in MATLAB file. This will be the name of the
    ! variable containing the matrix data when reading the file into matlab.
    character(LEN=*), intent(IN) :: smatrixName

    ! Source matrix
    type(t_matrixScalar), intent(IN) :: rmatrix
    
    ! Whether to spy the real data of the matrix or only its sparsity
    ! pattern
    logical, intent(IN) :: bdata

    ! OPTIONAL: status of file
    integer, intent(IN), optional :: cstatus

    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for beter visualisation.
    ! If not present, a default of 1E-12 is assumed.
    real(DP), intent(IN), optional :: dthreshold
!</input>
!</subroutine>
    
    ! local variables
    real(DP), dimension(:), pointer :: Da
    real(SP), dimension(:), pointer :: Fa
    integer, dimension(:), pointer :: Kld,Kcol
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
      
      call lsyssc_getbase_Kld(rmatrix,Kld)
      call lsyssc_getbase_Kcol(rmatrix,Kcol)
      
      if (bdata) then
        ! Which matrix type are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_double(rmatrix,Da)
          select case(rmatrix%cinterleavematrixFormat)
          case (LSYSSC_MATRIXD)
            call do_spy_mat79matD_double(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,&
                Kld,Kcol,Da,dthres)
          case (LSYSSC_MATRIX1)
            call do_spy_mat79mat1_double(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,&
                rmatrix%NVAR,Kld,Kcol,Da,dthres)
          case DEFAULT
            call output_line ('Unsupported interleave matrix type!', &
                              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spyMatrix')
            call sys_halt()
          end select
          write(UNIT=iunit,FMT=30)
          
        case (ST_SINGLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_single(rmatrix,Fa)
          select case(rmatrix%cinterleavematrixFormat)
          case (LSYSSC_MATRIXD)
            call do_spy_mat79matD_single(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,&
                Kld,Kcol,Fa,dthres)
          case (LSYSSC_MATRIX1)
            call do_spy_mat79mat1_single(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,&
                rmatrix%NVAR,Kld,Kcol,Fa,dthres)
          case DEFAULT
            call output_line ('Unsupported interleave matrix type!', &
                              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spyMatrix')
            call sys_halt()
          end select
          write(UNIT=iunit,FMT=30)

        case DEFAULT
          call output_line ('Unsupported matrix type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spyMatrix')
          call sys_halt()
        end select
        
      else
        
        ! Output only matrix structure
        write(UNIT=iunit,FMT=10)
        select case(rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIXD)
          call do_spy_mat79matD_double(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,&
              Kld,Kcol)
        case (LSYSSC_MATRIX1)
          call do_spy_mat79mat1_double(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,&
              rmatrix%NVAR,Kld,Kcol)
        case DEFAULT
          call output_line ('Unsupported interleave matrix type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spyMatrix')
          call sys_halt()
        end select
        write(UNIT=iunit,FMT=30)
      end if

    case (LSYSSC_MATRIX7,LSYSSC_MATRIX9)
      
      call lsyssc_getbase_Kld(rmatrix,Kld)
      call lsyssc_getbase_Kcol(rmatrix,Kcol)
      
      if (bdata) then
        ! Which matrix type are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_double(rmatrix,Da)
          call do_spy_mat79matD_double(rmatrix%NEQ,rmatrix%NCOLS,1,Kld,Kcol,Da,dthres)
          write(UNIT=iunit,FMT=30)

        case (ST_SINGLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_single(rmatrix,Fa)
          call do_spy_mat79matD_single(rmatrix%NEQ,rmatrix%NCOLS,1,Kld,Kcol,Fa,dthres)
          write(UNIT=iunit,FMT=30)
          
        case DEFAULT
          call output_line ('Unsupported matrix type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spyMatrix')
          call sys_halt()
        end select

      else
        
        ! Output only matrix structure
        write(UNIT=iunit,FMT=10)
        call do_spy_mat79matD_double(rmatrix%NEQ,rmatrix%NCOLS,1,Kld,Kcol)
        write(UNIT=iunit,FMT=30)
      end if
      
    case(LSYSSC_MATRIXD)
      
      if (bdata) then
        ! Which matrix type are we?
        select case(rmatrix%cdataType)
        case (ST_DOUBLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_double(rmatrix,Da)
          do ieq=1,rmatrix%NEQ
            if (abs(Da(ieq)) .ge. dthres) then
              write(UNIT=iunit,FMT=20) ieq,ieq,Da(ieq)
            end if
          end do
          write(UNIT=iunit,FMT=30)

        case (ST_SINGLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_single(rmatrix,Fa)
          do ieq=1,rmatrix%NEQ
            if (abs(Fa(ieq)) .ge. dthres) then
              write(UNIT=iunit,FMT=20) ieq,ieq,Fa(ieq)
            end if
          end do
          write(UNIT=iunit,FMT=30)

        case DEFAULT
          call output_line ('Unsupported matrix type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spyMatrix')
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
          call lsyssc_getbase_double(rmatrix,Da)
          call do_spy_mat1_double(rmatrix%NEQ,rmatrix%NCOLS,Da,dthres)
          write(UNIT=iunit,FMT=30)

        case (ST_SINGLE)
          write(UNIT=iunit,FMT=10)
          call lsyssc_getbase_single(rmatrix,Fa)
          call do_spy_mat1_single(rmatrix%NEQ,rmatrix%NCOLS,Fa,dthres)
          write(UNIT=iunit,FMT=30)

        case DEFAULT
          call output_line ('Unsupported matrix type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spyMatrix')
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
                         OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spyMatrix')
      call sys_halt()
      
    end select
    
    ! Close file
    write(UNIT=iunit,FMT=40) smatrixName
    close(UNIT=iunit)

10  format("data=[...")
20  format(I10,1X,I10,1X,E15.8,";")
30  format("];")
40  format(A,"=sparse(data(:,1),data(:,2),data(:,3));")

  contains

    ! Here, the real SPY routines follow.

    !**************************************************************
    ! SPY CSR matrix in double precision
    
    subroutine do_spy_mat79matD_double(neq,ncols,nvar,Kld,Kcol,Da,dthres)
      integer, dimension(:), intent(IN)  :: Kld
      integer, dimension(:), intent(IN)  :: Kcol
      integer, intent(IN)                :: neq,ncols
      integer, intent(IN)                             :: nvar
      real(DP), dimension(nvar,*), intent(IN), optional :: Da
      real(DP), intent(IN), optional :: dthres
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
      integer, dimension(:), intent(IN)  :: Kld
      integer, dimension(:), intent(IN)  :: Kcol
      integer, intent(IN)                :: neq,ncols
      integer, intent(IN)                             :: nvar
      real(SP), dimension(nvar,*), intent(IN), optional :: Fa
      real(DP), intent(IN), optional :: dthres
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
      integer, dimension(:), intent(IN)  :: Kld
      integer, dimension(:), intent(IN)  :: Kcol
      integer, intent(IN)                :: neq,ncols
      integer, intent(IN)                             :: nvar,mvar
      real(DP), dimension(nvar,mvar,*), intent(IN), optional :: Da
      real(DP), intent(IN), optional :: dthres
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
      integer, dimension(:), intent(IN)  :: Kld
      integer, dimension(:), intent(IN)  :: Kcol
      integer, intent(IN)                :: neq,ncols
      integer, intent(IN)                             :: nvar,mvar
      real(SP), dimension(nvar,mvar,*), intent(IN), optional :: Fa
      real(DP), intent(IN), optional :: dthres
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
      integer, intent(IN) :: neq,ncols
      real(DP), dimension(:), intent(IN), optional :: Da
      real(DP), intent(IN), optional :: dthres
      real(DP) :: ddata
      integer :: ieq,icol

      if (present(Da)) then
        do ieq=1,neq
          do icol=1,ncols
            ddata = Da(ncols*(ieq-1)+icol)
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
                ieq,icol,Da(ncols*(ieq-1)+icol)
          end do
        end do
      end if
    end subroutine do_spy_mat1_double

    !**************************************************************
    ! SPY full matrix in single precision
    
    subroutine do_spy_mat1_single(neq,ncols,Fa,dthres)
      integer, intent(IN) :: neq,ncols
      real(SP), dimension(:), intent(IN), optional :: Fa
      real(DP), intent(IN), optional :: dthres
      real(SP) :: fdata
      integer :: ieq,icol

      if (present(Fa)) then
        do ieq=1,neq
          do icol=1,ncols
            fdata = Fa(ncols*(ieq-1)+icol)
            if (abs(fdata) .ge. dthres) then
              write(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                  ieq,icol,Fa(ncols*(ieq-1)+icol)
            end if
          end do
        end do
      else
        do ieq=1,neq
          do icol=1,ncols
            write(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                ieq,icol,Da(ncols*(ieq-1)+icol)
          end do
        end do
      end if
    end subroutine do_spy_mat1_single
  end subroutine matio_spyMatrix

end module

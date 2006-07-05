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

MODULE matrixio

  USE fsystem
  USE storage
  USE io
  USE linearsystemscalar
  
  IMPLICIT NONE

  CONTAINS

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE matio_writeMatrixHR (rmatrix, sarray,&
                                  bnoZero, ifile, sfile, sformat)
  
  !<description>
    ! This routine writes a scalar matrix into a text file.
    ! The matrix is written in human readable form.
  !</description>
    
  !<input>
    ! The matrix to be written out
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
    
    ! Name of the matrix
    CHARACTER(len=*), INTENT(IN) :: sarray
    
    ! Suppress zeroes in output: yes/no
    LOGICAL, INTENT(IN) :: bnoZero
    
    ! Output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN) :: ifile
    
    ! Name of the file where to write to. Only relevant for ifile=0!
    CHARACTER(len=*), INTENT(IN) :: sfile
    
    ! Format string to use for the output; e.g. '(E20.10)'
    CHARACTER(len=*), INTENT(IN) :: sformat
    
  !</input>
    
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_DA
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld

  ! Depending on the matrix format, choose the right routine for writing
  SELECT CASE (rmatrix%cmatrixFormat)
  CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
    ! Matrix precision?
    SELECT CASE (rmatrix%cdataType)
    CASE (ST_DOUBLE)
      ! Get the data arrays and write the matrix
      CALL storage_getbase_double (rmatrix%h_Da,p_Da)
      CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
      CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)
      CALL matio_writeMatrix79_Dble (p_Da, p_Kcol, p_Kld, &
                                      rmatrix%NEQ, rmatrix%NEQ, sarray, &
                                      bnoZero, ifile, sfile, sformat)
    CASE DEFAULT
      PRINT *,'matio_writeFullMatrix: Unsupported matrix precision.'
      STOP
    END SELECT
  CASE DEFAULT
    PRINT *,'matio_writeFullMatrix: Unknown matrix format.'
    STOP
  END SELECT
    
  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE matio_writeMatrix1_Dble (Da, nrow, ncol, sarray, &
                                       bnoZero, ifile, sfile, sformat)
  
  !<description>
    ! Write full double precision matrix into a text file.
    !
    ! This writes an array Da with nrow rows and ncol columns as a matrix
    ! into a text file.
  !</description>
    
  !<input>
    ! number of rows
    INTEGER(I32), INTENT(IN) :: nrow
    
    ! number of columns
    INTEGER(I32), INTENT(IN) :: ncol
    
    ! matrix: array [1..nrow,1..ncol] of double
    REAL(DP), DIMENSION(nrow,ncol), INTENT(IN) :: Da
    
    ! name of the matrix
    CHARACTER(len=*), INTENT(IN) :: sarray
    
    ! suppress zeroes in output: yes/no
    LOGICAL, INTENT(IN) :: bnoZero
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    CHARACTER(len=*), INTENT(IN) :: sfile
    
    ! format string to use for the output; e.g. '(D20.10)'
    CHARACTER(len=*), INTENT(IN) :: sformat
    
  !</input>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, cf, nchar
    REAL(DP) :: dval
    CHARACTER(len=128) :: S
    CHARACTER(len=6) :: sformatChar
    
    IF (ifile .EQ. 0) THEN
      CALL io_openFileForWriting(sfile, cf, SYS_REPLACE)
      IF (cf .EQ. -1) THEN
        PRINT *, 'matio_writeFullMatrix: Could not open file '// &
                 trim(sfile)
        STOP
      END IF
    ELSE
      cf = ifile
    END IF
    
    ! Get length of output strings
    S(:) = ' '
    WRITE (S,sformat) 0.0_DP
    nchar = LEN(trim(S))
    
    ! Build array format string
    sformatChar = '(A'//sys_i3(nchar)//')'
    
    ! Write all format strings into the file
    WRITE (cf,'(3A15,L2,3I15)') sarray, sformat, sformatChar, &
                                bnoZero, nchar, nrow, ncol
    
    IF (nrow .LE. 0) RETURN
    IF (ncol .LE. 0) RETURN

    ! Write the matrix
    DO i=1, nrow
      DO j=1, ncol-1
        dval = Da(i,j)
        IF ((.NOT. bnoZero) .OR. (dval .NE. 0.0_DP)) THEN
          WRITE (cf,sformat,ADVANCE='NO') dval
        ELSE
          WRITE (cf,sformatChar, ADVANCE='NO') '.'
        END IF
      END DO
      dval = Da(i,j)
      IF ((.NOT. bnoZero) .OR. (dval .NE. 0.0_DP)) THEN
        WRITE (cf,sformat) dval
      ELSE
        WRITE (cf,sformatChar) '.'
      END IF
    END DO
    
    ! Close the file if necessary
    IF (ifile .EQ. 0) CLOSE(cf)
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE matio_writeMatrix79_Dble (Da, Icol, Irow, &
                                        nrow, ncol, sarray, &
                                        bnoZero, ifile, sfile, sformat)
  
  !<description>
    ! Write sparse double precision matrix in matrix format 9 or
    ! matrix format 9 into a text file.
    !
    ! Double-precision version
  !</description>
    
  !<input>
    ! number of rows
    INTEGER(I32), INTENT(IN) :: nrow
    
    ! number of columns; must be =nrow for structure-7 matrices
    INTEGER(I32), INTENT(IN) :: ncol
    
    ! matrix: array [1..na] of double
    REAL(DP), DIMENSION(:), INTENT(IN) :: Da
    
    ! Column structure of the matrix
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: Icol
    
    ! Row structure of the matrix
    INTEGER(I32), DIMENSION(nrow+1), INTENT(IN) :: Irow
    
    ! name of the matrix
    CHARACTER(len=*), INTENT(IN) :: sarray
    
    ! suppress zeroes in output: yes/no
    LOGICAL, INTENT(IN) :: bnoZero
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel ifile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN) :: ifile
    
    ! name of the file where to write to. Only relevant for ifile=0!
    CHARACTER(len=*), INTENT(IN) :: sfile
    
    ! format string to use for the output; e.g. '(D20.10)'
    CHARACTER(len=*), INTENT(IN) :: sformat
    
  !</input>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, k, cf, nchar
    REAL(DP) :: dval
    CHARACTER(len=128) :: S
    CHARACTER(len=6) :: sformatChar
    INTEGER :: h_DrowVec
    REAL(DP), DIMENSION(:), POINTER :: p_DrowVec
    
    IF (ifile .EQ. 0) THEN
      CALL io_openFileForWriting(sfile, cf, SYS_REPLACE)
      IF (cf .EQ. -1) THEN
        PRINT *, 'matio_writeFullMatrix: Could not open file '// &
                 trim(sfile)
        STOP
      END IF
    ELSE
      cf = ifile
    END IF
    
    ! Get length of output strings
    S(:) = ' '
    WRITE (S,sformat) 0.0_DP
    nchar = LEN(trim(S))
    
    ! Build array format string
    sformatChar = '(A'//sys_i3(nchar)//')'
    
    ! Write all format strings into the file
    WRITE (cf,'(3A15,L2,3I15)') sarray, sformat, sformatChar, &
                                bnoZero, nchar, nrow, ncol

    IF (nrow .LE. 0) RETURN
    IF (ncol .LE. 0) RETURN
    
    ! Write the matrix
    CALL storage_new('matio_writeSparseMatrix','DrowVec',ncol, &
                     ST_DOUBLE, h_DrowVec, ST_NEWBLOCK_NOINIT)
    CALL storage_getbase_double(h_DrowVec, p_DrowVec)
    DO i=1, nrow
      ! Extract row i
      DO j=1, ncol
        p_DrowVec(j) = 0.0_DP
      END DO
      DO j=0, Irow(i+1)-Irow(i)-1
        k = Irow(i)+j
        p_DrowVec(Icol(k)) = Da(k)
        IF (Da(k) .EQ. 0.0_DP) p_DrowVec(Icol(k)) = SYS_MAXREAL
      END DO
      ! Write row i
      DO j=1, ncol-1
        dval = p_DrowVec(j)
        IF ((.NOT. bnoZero) .OR. (dval .NE. 0.0_DP)) THEN
          IF (dval .EQ. SYS_MAXREAL) THEN
            WRITE (cf,sformatChar, ADVANCE='NO') '0'
          ELSE
            WRITE (cf,sformat,ADVANCE='NO') dval
          END IF
        ELSE
          WRITE (cf,sformatChar, ADVANCE='NO') '.'
        END IF
      END DO
      dval = p_DrowVec(ncol)
      IF ((.NOT. bnoZero) .OR. (dval .NE. 0.0_DP)) THEN
        IF (dval .EQ. SYS_MAXREAL) THEN
          WRITE (cf,sformatChar) '0'
        ELSE
          WRITE (cf,sformat) dval
        END IF
      ELSE
        WRITE (cf,sformatChar) '.'
      END IF
    END DO
    CALL storage_free(h_DrowVec)
    
    ! Close the file if necessary
    IF (ifile .EQ. 0) CLOSE(cf)
  
  END SUBROUTINE

END MODULE

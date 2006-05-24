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

MODULE matout

  USE fsystem
  USE storage
  USE io
  
  IMPLICIT NONE

  CONTAINS

!<subroutine>
  SUBROUTINE matout_writeFullMatrix (Da, nrow, ncol, sarray, &
                                     bnoZero, cfile, sfile, sformat)
  
  !<description>
    ! Write full double precision matrix into a text file.
    !
    ! This writes an array Da with nrow rows and ncol columns as a matrix
    ! into a text file.
  !</description>
    
  !<input>
    ! matrix: array [1..nrow,1..ncol] of double
    REAL(DP), DIMENSION(nrow,ncol), INTENT(IN) :: Da
    
    ! number of rows
    INTEGER(I32), INTENT(IN) :: nrow
    
    ! number of columns
    INTEGER(I32), INTENT(IN) :: ncol
    
    ! name of the matrix
    CHARACTER(len=*), INTENT(IN) :: sarray
    
    ! suppress zeroes in output: yes/no
    LOGICAL, INTENT(IN) :: bnoZero
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel cfile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN) :: cfile
    
    ! name of the file where to write to. Only relevant for cfile=0!
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
    
    IF (cfile .EQ. 0) THEN
      CALL io_openFileForWriting(sfile, cf, SYS_REPLACE)
      IF (cf .EQ. -1) THEN
        PRINT *, 'matout_writeFullMatrix: Could not open file '// &
                 trim(sfile)
        STOP
      END IF
    ELSE
      cf = cfile
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
    IF (cfile .EQ. 0) CLOSE(cf)
  
  END SUBROUTINE matout_writeFullMatrix

!<subroutine>
  SUBROUTINE matout_writeSparseMatrix (Da, Icol, Irow, &
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
    REAL(DP), DIMENSION(:), INTENT(IN) :: Da
    
    ! Column structure of the matrix
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: Icol
    
    ! Row structure of the matrix
    INTEGER(I32), DIMENSION(nrow+1), INTENT(IN) :: Irow
    
    ! number of rows
    INTEGER(I32), INTENT(IN) :: nrow
    
    ! number of columns; must be =nrow for structure-7 matrices
    INTEGER(I32), INTENT(IN) :: ncol
    
    ! name of the matrix
    CHARACTER(len=*), INTENT(IN) :: sarray
    
    ! suppress zeroes in output: yes/no
    LOGICAL, INTENT(IN) :: bnoZero
    
    ! output channel to use for output
    !  = 0: Get temporary channel for file 'sfile'
    ! <> 0: Write to channel cfile. Don't close the channel afterwards.
    !       'sfile' is ignored.
    INTEGER(I32), INTENT(IN) :: cfile
    
    ! name of the file where to write to. Only relevant for cfile=0!
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
    
    IF (cfile .EQ. 0) THEN
      CALL io_openFileForWriting(sfile, cf, SYS_REPLACE)
      IF (cf .EQ. -1) THEN
        PRINT *, 'matout_writeFullMatrix: Could not open file '// &
                 trim(sfile)
        STOP
      END IF
    ELSE
      cf = cfile
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
    CALL storage_new('matout_writeSparseMatrix','DrowVec',ncol, &
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
        IF (Da(k) .EQ. 0.0_DP) p_DrowVec(Icol(k)) = 1D-99
      END DO
      ! Write row i
      DO j=1, ncol-1
        dval = p_DrowVec(j)
        IF ((.NOT. bnoZero) .OR. (dval .NE. 0.0_DP)) THEN
          IF (ABS(dval) .LE. 1.0D-90) THEN
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
        IF (ABS(dval) .LE. 1.0D-90) THEN
          WRITE (cf,sformatChar, ADVANCE='NO') '0'
        ELSE
          WRITE (cf,sformat,ADVANCE='NO') dval
        END IF
      ELSE
        WRITE (cf,sformatChar, ADVANCE='NO') '.'
      END IF
    END DO
    CALL storage_free(h_DrowVec)
    
    ! Close the file if necessary
    IF (cfile .EQ. 0) CLOSE(cf)
  
  END SUBROUTINE matout_writeSparseMatrix

END MODULE

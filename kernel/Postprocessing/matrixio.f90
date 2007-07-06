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
!# 1.) matio_writeBlockMatrixHR
!#     -> Writes a block matrix in human readable form into a text file
!#
!# 3.) matio_spyMatrix
!#     -> Writes a scalar matrix into a file which can be visualized
!#        by means of the MATLAB command SPY
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
  USE linearsystemblock
  USE globalsystem
  
  IMPLICIT NONE

  CONTAINS

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE matio_writeBlockMatrixHR (rmatrix, sarray,&
                                       bnoZero, ifile, sfile, sformat, dthreshold)
  
  !<description>
    ! This routine writes a block matrix into a text file.
    ! The matrix is written in human readable form.
    ! Note that for this purpose, a new matrix is temporarily created in memory!
  !</description>
    
  !<input>
    ! The matrix to be written out
    TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
    
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
    
    ! OPTIONAL: Threshold parameter for the entries. Entries whose absolute
    ! value is below this threshold are replaced by 0.0 for beter visualisation.
    ! If not present, a default of 1E-12 is assumed.
    REAL(DP), INTENT(IN), OPTIONAL :: dthreshold
  !</input>
    
!</subroutine>

    ! local variables
    TYPE(t_matrixBlock) :: rtempMatrix
    REAL(DP), DIMENSION(:), POINTER :: p_DA
    REAL(DP) :: dthres 

    ! We have to create a global matrix first!
    CALL glsys_assembleGlobal (rmatrix,rtempMatrix,.TRUE.,.TRUE.)
                              
    ! Replace small values by zero
    dthres = 1E-12_DP
    IF (PRESENT(dthreshold)) dthres = dthreshold
    IF (ABS(dthres) .GT. 0.0_DP) THEN
      CALL lsyssc_getbase_double (rtempMatrix%RmatrixBlock(1,1),p_DA)
      WHERE (abs(p_Da) .LT. dthres) p_Da = 0.0_DP
    END IF
    
    ! Write matrix to the file
    CALL matio_writeMatrixHR (rtempMatrix%RmatrixBlock(1,1), sarray,&
                              bnoZero, ifile, sfile, sformat)

    ! Release the temporary matrix
    CALL lsysbl_releaseMatrix (rtempMatrix)

  END SUBROUTINE 

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
      CALL lsyssc_getbase_double (rmatrix,p_Da)
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      CALL matio_writeMatrix79_Dble (p_Da, p_Kcol, p_Kld, &
                                      rmatrix%NEQ, rmatrix%NEQ, sarray, &
                                      bnoZero, ifile, sfile, sformat)
    CASE DEFAULT
      PRINT *,'matio_writeFullMatrix: Unsupported matrix precision.'
      CALL sys_halt()
    END SELECT
  CASE DEFAULT
    PRINT *,'matio_writeFullMatrix: Unknown matrix format.'
    CALL sys_halt()
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
    INTEGER :: i, j, cf, nchar
    REAL(DP) :: dval
    CHARACTER(len=128) :: S
    CHARACTER(len=6) :: sformatChar
    
    IF (ifile .EQ. 0) THEN
      CALL io_openFileForWriting(sfile, cf, SYS_REPLACE)
      IF (cf .EQ. -1) THEN
        PRINT *, 'matio_writeFullMatrix: Could not open file '// &
                 trim(sfile)
        CALL sys_halt()
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
    INTEGER :: i, j, k, cf, nchar
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
        CALL sys_halt()
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
      IF (bnoZero) THEN
        ! SYS_MAXREAL is written out as '.'
        DO j=1, ncol
          p_DrowVec(j) = SYS_MAXREAL
        END DO
      ELSE
        DO j=1, ncol
          p_DrowVec(j) = 0.0_DP
        END DO
      END IF
      
      DO j=0, Irow(i+1)-Irow(i)-1
        k = Irow(i)+j
        p_DrowVec(Icol(k)) = Da(k)
      END DO
      
      ! Write row i
      DO j=1, ncol-1
        dval = p_DrowVec(j)
        IF (bnoZero .AND. (dval .EQ. SYS_MAXREAL)) THEN
          WRITE (cf,sformatChar, ADVANCE='NO') '.'
        ELSE
          WRITE (cf,sformat,ADVANCE='NO') dval
        END IF
      END DO
      
      dval = p_DrowVec(ncol)
      IF (bnoZero .AND. (dval .EQ. SYS_MAXREAL)) THEN
        WRITE (cf,sformatChar, ADVANCE='YES') '.'
      ELSE
        WRITE (cf,sformat,ADVANCE='YES') dval
      END IF
    END DO
    
    CALL storage_free(h_DrowVec)
    
    ! Close the file if necessary
    IF (ifile .EQ. 0) CLOSE(cf)
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE matio_spyMatrix(sfilename,smatrixName,rmatrix,bdata,cstatus)

!<description>
    ! Writes a scalar matrix to file so that its sparsity structure
    ! can be visualized in MATLAB by means of the spy command
!</description>

!<input>
    ! file name of the MATLAB file without fileextension
    CHARACTER(LEN=*), INTENT(IN) :: sfileName
    
    ! name of the matrix in MATLAB file
    CHARACTER(LEN=*), INTENT(IN) :: smatrixName

    ! Source matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
    
    ! Whether to spy the real data of the matrix or only its sparsity
    ! pattern
    LOGICAL, INTENT(IN) :: bdata

    ! OPTIONAL: status of file
    INTEGER, INTENT(IN), OPTIONAL :: cstatus
!</input>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: Da
    REAL(SP), DIMENSION(:), POINTER :: Fa
    INTEGER(I32), DIMENSION(:), POINTER :: Kld,Kcol
    INTEGER :: iunit,ieq
    CHARACTER(LEN=10) :: cstat,cpos

    ! Set file status if required
    IF (PRESENT(cstatus)) THEN
      SELECT CASE(cstatus)
      CASE (IO_NEW)
        cstat="NEW"; cpos="ASIS"
      CASE (IO_REPLACE)
        cstat="REPLACE"; cpos="ASIS"
      CASE (IO_OLD)
        cstat="OLD"; cpos="APPEND"
      CASE DEFAULT
        cstat="UNKNOWN"; cpos ="ASIS"
      END SELECT
    ELSE
      cstat="UNKNOWN"; cpos="ASIS"
    END IF
    
    ! Open output file
    iunit=sys_getFreeUnit()
    OPEN (UNIT=iunit,STATUS=TRIM(cstat),POSITION=TRIM(cpos),FILE=TRIM(ADJUSTL(sfilename))//'.m')
    
    ! Which matrix format are we?
    SELECT CASE(rmatrix%cmatrixFormat)
      
    CASE (LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9INTL)
      
      CALL lsyssc_getbase_Kld(rmatrix,Kld)
      CALL lsyssc_getbase_Kcol(rmatrix,Kcol)
      
      IF (bdata) THEN
        ! Which matrix type are we?
        SELECT CASE(rmatrix%cdataType)
        CASE (ST_DOUBLE)
          WRITE(UNIT=iunit,FMT=10)
          CALL lsyssc_getbase_double(rmatrix,Da)
          SELECT CASE(rmatrix%cinterleavematrixFormat)
          CASE (LSYSSC_MATRIXD)
            CALL do_spy_mat79matD_double(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,Kld,Kcol,Da)
          CASE (LSYSSC_MATRIX1)
            CALL do_spy_mat79mat1_double(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,rmatrix%NVAR,Kld,Kcol,Da)
          CASE DEFAULT
            PRINT *, 'lsyssc_spyMatrix: Unsupported interleave matrix type!'
            CALL sys_halt()
          END SELECT
          WRITE(UNIT=iunit,FMT=30)
          
        CASE (ST_SINGLE)
          WRITE(UNIT=iunit,FMT=10)
          CALL lsyssc_getbase_single(rmatrix,Fa)
          SELECT CASE(rmatrix%cinterleavematrixFormat)
          CASE (LSYSSC_MATRIXD)
            CALL do_spy_mat79matD_single(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,Kld,Kcol,Fa)
          CASE (LSYSSC_MATRIX1)
            CALL do_spy_mat79mat1_single(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,rmatrix%NVAR,Kld,Kcol,Fa)
          CASE DEFAULT
            PRINT *, 'lsyssc_spyMatrix: Unsupported interleave matrix type!'
            CALL sys_halt()
          END SELECT
          WRITE(UNIT=iunit,FMT=30)

        CASE DEFAULT
          PRINT *, 'lsyssc_spyMatrix: Unsupported matrix type!'
          CALL sys_halt()
        END SELECT
        
      ELSE
        
        ! Output only matrix structure
        WRITE(UNIT=iunit,FMT=10)
        SELECT CASE(rmatrix%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIXD)
          CALL do_spy_mat79matD_double(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,Kld,Kcol)
        CASE (LSYSSC_MATRIX1)
          CALL do_spy_mat79mat1_double(rmatrix%NEQ,rmatrix%NCOLS,rmatrix%NVAR,rmatrix%NVAR,Kld,Kcol)
        CASE DEFAULT
          PRINT *, 'lsyssc_spyMatrix: Unsupported interleave matrix type!'
          CALL sys_halt()
        END SELECT
        WRITE(UNIT=iunit,FMT=30)
      END IF

    CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX9)
      
      CALL lsyssc_getbase_Kld(rmatrix,Kld)
      CALL lsyssc_getbase_Kcol(rmatrix,Kcol)
      
      IF (bdata) THEN
        ! Which matrix type are we?
        SELECT CASE(rmatrix%cdataType)
        CASE (ST_DOUBLE)
          WRITE(UNIT=iunit,FMT=10)
          CALL lsyssc_getbase_double(rmatrix,Da)
          CALL do_spy_mat79matD_double(rmatrix%NEQ,rmatrix%NCOLS,1,Kld,Kcol,Da)
          WRITE(UNIT=iunit,FMT=30)

        CASE (ST_SINGLE)
          WRITE(UNIT=iunit,FMT=10)
          CALL lsyssc_getbase_single(rmatrix,Fa)
          CALL do_spy_mat79matD_single(rmatrix%NEQ,rmatrix%NCOLS,1,Kld,Kcol,Fa)
          WRITE(UNIT=iunit,FMT=30)
          
        CASE DEFAULT
          PRINT *, 'lsyssc_spyMatrix: Unsupported matrix type!'
          CALL sys_halt()
        END SELECT

      ELSE
        
        ! Output only matrix structure
        WRITE(UNIT=iunit,FMT=10)
        CALL do_spy_mat79matD_double(rmatrix%NEQ,rmatrix%NCOLS,1,Kld,Kcol)
        WRITE(UNIT=iunit,FMT=30)
      END IF
      
    CASE(LSYSSC_MATRIXD)
      
      IF (bdata) THEN
        ! Which matrix type are we?
        SELECT CASE(rmatrix%cdataType)
        CASE (ST_DOUBLE)
          WRITE(UNIT=iunit,FMT=10)
          CALL lsyssc_getbase_double(rmatrix,Da)
          DO ieq=1,rmatrix%NEQ
            WRITE(UNIT=iunit,FMT=20) ieq,ieq,Da(ieq)
          END DO
          WRITE(UNIT=iunit,FMT=30)

        CASE (ST_SINGLE)
          WRITE(UNIT=iunit,FMT=10)
          CALL lsyssc_getbase_single(rmatrix,Fa)
          DO ieq=1,rmatrix%NEQ
            WRITE(UNIT=iunit,FMT=20) ieq,ieq,Fa(ieq)
          END DO
          WRITE(UNIT=iunit,FMT=30)

        CASE DEFAULT
          PRINT *, 'lsyssc_spyMatrix: Unsupported matrix type!'
          CALL sys_halt()
        END SELECT
        
      ELSE
        
        ! Output only matrix structure
        WRITE(UNIT=iunit,FMT=10)
        DO ieq=1,rmatrix%NEQ
          WRITE(UNIT=iunit,FMT=20) ieq,ieq,1
        END DO
        WRITE(UNIT=iunit,FMT=30)
        
      END IF
      
    CASE(LSYSSC_MATRIX1)

      IF (bdata) THEN
        ! Which matrix type are we?
        SELECT CASE(rmatrix%cdataType)
        CASE (ST_DOUBLE)
          WRITE(UNIT=iunit,FMT=10)
          CALL lsyssc_getbase_double(rmatrix,Da)
          CALL do_spy_mat1_double(rmatrix%NEQ,rmatrix%NCOLS,Da)
          WRITE(UNIT=iunit,FMT=30)

        CASE (ST_SINGLE)
          WRITE(UNIT=iunit,FMT=10)
          CALL lsyssc_getbase_single(rmatrix,Fa)
          CALL do_spy_mat1_single(rmatrix%NEQ,rmatrix%NCOLS,Fa)
          WRITE(UNIT=iunit,FMT=30)

        CASE DEFAULT
          PRINT *, 'lsyssc_spyMatrix: Unsupported matrix type!'
          CALL sys_halt()
        END SELECT
        
      ELSE
        
        ! Output only matrix structure
        WRITE(UNIT=iunit,FMT=10)
        CALL do_spy_mat1_double(rmatrix%NEQ,rmatrix%NCOLS)
        WRITE(UNIT=iunit,FMT=30)

      END IF
      
    CASE DEFAULT
      PRINT *, 'lsyssc_spyMatrix: Unsoppurted matrix format!'
      CALL sys_halt()
      
    END SELECT
    
    ! Close file
    WRITE(UNIT=iunit,FMT=40) smatrixName
    CLOSE(UNIT=iunit)

10  FORMAT("data=[...")
20  FORMAT(I10,1X,I10,1X,E15.8,";")
30  FORMAT("];")
40  FORMAT(A,"=sparse(data(:,1),data(:,2),data(:,3));")

  CONTAINS

    ! Here, the real SPY routines follow.

    !**************************************************************
    ! SPY CSR matrix in double precision
    
    SUBROUTINE do_spy_mat79matD_double(neq,ncols,nvar,Kld,Kcol,Da)
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)  :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)  :: Kcol
      INTEGER(PREC_VECIDX), INTENT(IN)                :: neq,ncols
      INTEGER, INTENT(IN)                             :: nvar
      REAL(DP), DIMENSION(nvar,*), INTENT(IN), OPTIONAL :: Da
      INTEGER :: ieq,ild,ivar

      
      IF (PRESENT(Da)) THEN
        DO ieq=1,neq
          DO ild=Kld(ieq),Kld(ieq+1)-1
            DO ivar=1,nvar
              WRITE(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                  (ivar-1)*neq+ieq,(ivar-1)*ncols+Kcol(ild),Da(ivar,ild)
            END DO
          END DO
        END DO
      ELSE
        DO ieq=1,neq
          DO ild=Kld(ieq),Kld(ieq+1)-1
            DO ivar=1,nvar
              WRITE(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                  (ivar-1)*neq+ieq,(ivar-1)*ncols+Kcol(ild)
            END DO
          END DO
        END DO
      END IF
    END SUBROUTINE do_spy_mat79matD_double

    !**************************************************************
    ! SPY CSR matrix in single precision
    
    SUBROUTINE do_spy_mat79matD_single(neq,ncols,nvar,Kld,Kcol,Fa)
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)  :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)  :: Kcol
      INTEGER(PREC_VECIDX), INTENT(IN)                :: neq,ncols
      INTEGER, INTENT(IN)                             :: nvar
      REAL(SP), DIMENSION(nvar,*), INTENT(IN), OPTIONAL :: Fa
      INTEGER :: ieq,ild,ivar

      IF (PRESENT(Fa)) THEN
        DO ieq=1,neq
          DO ild=Kld(ieq),Kld(ieq+1)-1
            DO ivar=1,nvar
              WRITE(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                  (ivar-1)*neq+ieq,(ivar-1)*ncols+Kcol(ild),Fa(ivar,ild)
            END DO
          END DO
        END DO
      ELSE
        DO ieq=1,neq
          DO ild=Kld(ieq),Kld(ieq+1)-1
            DO ivar=1,nvar
              WRITE(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                  (ivar-1)*neq+ieq,(ivar-1)*ncols+Kcol(ild)
            END DO
          END DO
        END DO
      END IF
    END SUBROUTINE do_spy_mat79matD_single

    !**************************************************************
    ! SPY CSR matrix in double precision
    
    SUBROUTINE do_spy_mat79mat1_double(neq,ncols,nvar,mvar,Kld,Kcol,Da)
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)  :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)  :: Kcol
      INTEGER(PREC_VECIDX), INTENT(IN)                :: neq,ncols
      INTEGER, INTENT(IN)                             :: nvar,mvar
      REAL(DP), DIMENSION(nvar,mvar,*), INTENT(IN), OPTIONAL :: Da
      INTEGER :: ieq,ild,ivar,jvar

      IF (PRESENT(Da)) THEN
        DO ieq=1,neq
          DO ild=Kld(ieq),Kld(ieq+1)-1
            DO jvar=1,nvar ! local column index
              DO ivar=1,mvar ! local row index
                WRITE(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                    (ivar-1)*neq+ieq,(jvar-1)*ncols+Kcol(ild),Da(ivar,jvar,ild)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO ieq=1,neq
          DO ild=Kld(ieq),Kld(ieq+1)-1
            DO jvar=1,nvar ! local column index
              DO ivar=1,mvar ! local row index
                WRITE(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                    (ivar-1)*neq+ieq,(jvar-1)*ncols+Kcol(ild)
              END DO
            END DO
          END DO
        END DO
      END IF
    END SUBROUTINE do_spy_mat79mat1_double

    !**************************************************************
    ! SPY CSR matrix in double precision
    
    SUBROUTINE do_spy_mat79mat1_single(neq,ncols,nvar,mvar,Kld,Kcol,Fa)
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)  :: Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)  :: Kcol
      INTEGER(PREC_VECIDX), INTENT(IN)                :: neq,ncols
      INTEGER, INTENT(IN)                             :: nvar,mvar
      REAL(SP), DIMENSION(nvar,mvar,*), INTENT(IN), OPTIONAL :: Fa
      INTEGER :: ieq,ild,ivar,jvar

      IF (PRESENT(Fa)) THEN
        DO ieq=1,neq
          DO ild=Kld(ieq),Kld(ieq+1)-1
            DO jvar=1,nvar ! local column index
              DO ivar=1,mvar ! local row index
                WRITE(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                    (ivar-1)*neq+ieq,(jvar-1)*ncols+Kcol(ild),Fa(ivar,jvar,ild)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO ieq=1,neq
          DO ild=Kld(ieq),Kld(ieq+1)-1
            DO jvar=1,nvar ! local column index
              DO ivar=1,mvar ! local row index
                WRITE(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                    (ivar-1)*neq+ieq,(jvar-1)*ncols+Kcol(ild)
              END DO
            END DO
          END DO
        END DO
      END IF
    END SUBROUTINE do_spy_mat79mat1_single

    !**************************************************************
    ! SPY full matrix in double precision
    
    SUBROUTINE do_spy_mat1_double(neq,ncols,Da)
      INTEGER(PREC_VECIDX), INTENT(IN) :: neq,ncols
      REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: Da
      INTEGER :: ieq,icol

      IF (PRESENT(Da)) THEN
        DO ieq=1,neq
          DO icol=1,ncols
            WRITE(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                ieq,icol,Da(ncols*(ieq-1)+icol)
          END DO
        END DO
      ELSE
        DO ieq=1,neq
          DO icol=1,ncols
            WRITE(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                ieq,icol,Da(ncols*(ieq-1)+icol)
          END DO
        END DO
      END IF
    END SUBROUTINE do_spy_mat1_double

    !**************************************************************
    ! SPY full matrix in single precision
    
    SUBROUTINE do_spy_mat1_single(neq,ncols,Fa)
      INTEGER(PREC_VECIDX), INTENT(IN) :: neq,ncols
      REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: Fa
      INTEGER :: ieq,icol

      IF (PRESENT(Fa)) THEN
        DO ieq=1,neq
          DO icol=1,ncols
            WRITE(UNIT=iunit,FMT='(I10,1X,I10,1X,E15.8,";")') &
                ieq,icol,Fa(ncols*(ieq-1)+icol)
          END DO
        END DO
      ELSE
        DO ieq=1,neq
          DO icol=1,ncols
            WRITE(UNIT=iunit,FMT='(I10,1X,I10,1X,"1.0;")') &
                ieq,icol,Da(ncols*(ieq-1)+icol)
          END DO
        END DO
      END IF
    END SUBROUTINE do_spy_mat1_single
  END SUBROUTINE matio_spyMatrix

END MODULE

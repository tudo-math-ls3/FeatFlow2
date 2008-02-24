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
!# </purpose>
!##############################################################################

MODULE matrixmodification

  USE fsystem
  USE storage
  USE linearsystemscalar
  USE genoutput
  USE sort
  
  IMPLICIT NONE

CONTAINS
 
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mmod_replaceLinesByUnit (rmatrix,Irows)
  
!<description>
    ! This routine replaces some lines of a given scalar matrix by unit vectors.
!</description>

!<input>
    ! A list of row numbers of all the rows which are to be replaced
    ! by unit vectors.
    INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: Irows
!</input>

!<inputoutput>
    ! The matrix which is to be modified.
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>

!</subroutine>

    ! At first we must take care of the matrix type.
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9)
      CALL replaceLines_format9 (rmatrix,Irows)
    CASE (LSYSSC_MATRIX7)
      CALL replaceLines_format7 (rmatrix,Irows)
    END SELECT
    
  CONTAINS
    
    ! ****************************************
    ! The replacement routine for format 9
    
    SUBROUTINE replaceLines_format9 (rmatrix,Irows)
      
      INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: Irows
      TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
      
      ! local variables
      INTEGER(PREC_MATIDX) :: irow
      REAL(DP), DIMENSION(:), POINTER :: p_DA
      REAL(SP), DIMENSION(:), POINTER :: p_FA
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_Kdiagonal
      
      ! Get Kld and Kdiagonal
      CALL lsyssc_getbase_Kld(rmatrix,p_Kld)
      CALL lsyssc_getbase_Kdiagonal(rmatrix,p_Kdiagonal)
      
      ! Take care of the format of the entries
      SELECT CASE (rmatrix%cdataType)
      CASE (ST_DOUBLE)
        ! Get the data array
        CALL lsyssc_getbase_double(rmatrix,p_DA)
        IF (.NOT. ASSOCIATED(p_DA)) THEN
          CALL output_line('Matrix has no data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'replaceLines_format9')
          CALL sys_halt()
        END IF
        
        ! loop through the rows
        DO irow = 1,SIZE(Irows)
          
          ! Clear the row
          p_DA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_DP
          
          ! And put a unit vector there
          p_DA(p_Kdiagonal(Irows(irow))) = 1.0_DP
          
        END DO
        
      CASE (ST_SINGLE)
        ! Get the data array
        CALL lsyssc_getbase_single(rmatrix,p_FA)
        IF (.NOT. ASSOCIATED(p_FA)) THEN
          CALL output_line('Matrix has no data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'replaceLines_format9')
          CALL sys_halt()
        END IF
        
        ! loop through the rows
        DO irow = 1,SIZE(Irows)
          
          ! Clear the row
          p_FA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_SP
          
          ! And put a unit vector there
          p_FA(p_Kdiagonal(Irows(irow))) = 1.0_SP
          
        END DO
        
      CASE DEFAULT
        CALL output_line('Unsupported data format',&
            OU_CLASS_ERROR,OU_MODE_STD,'replaceLines_format9')
        CALL sys_halt()
      END SELECT
      
    END SUBROUTINE replaceLines_format9
    
    ! ****************************************
    ! The replacement routine for format 7
    
    SUBROUTINE replaceLines_format7 (rmatrix,Irows)
      
      INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: Irows
      TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
      
      ! local variables
      INTEGER(PREC_MATIDX) :: irow
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
      REAL(DP), DIMENSION(:), POINTER :: p_DA
      REAL(SP), DIMENSION(:), POINTER :: p_FA
      
      ! Get Kld:
      CALL lsyssc_getbase_Kld(rmatrix,p_Kld)
      
      ! Take care of the format of the entries
      SELECT CASE (rmatrix%cdataType)
      CASE (ST_DOUBLE)
        ! Get the data array
        CALL lsyssc_getbase_double(rmatrix,p_DA)
        IF (.NOT. ASSOCIATED(p_DA)) THEN
          CALL output_line('Matrix has no data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'replaceLines_format7')
          CALL sys_halt()
        END IF
        
        ! loop through the rows
        DO irow = 1,SIZE(Irows)
          
          ! Put a unit vector there
          p_DA(p_Kld(Irows(irow))) = 1.0_DP
          
          ! and clear the row
          p_DA(p_Kld(Irows(irow))+1:p_Kld(Irows(irow)+1)-1) = 0.0_DP
          
        END DO

      CASE (ST_SINGLE)
        ! Get the data array
        CALL lsyssc_getbase_single(rmatrix,p_FA)
        IF (.NOT. ASSOCIATED(p_FA)) THEN
          CALL output_line('Matrix has no data!',&
              OU_CLASS_ERROR,OU_MODE_STD,'replaceLines_format7')
          CALL sys_halt()
        END IF
        
        ! loop through the rows
        DO irow = 1,SIZE(Irows)
          
          ! Put a unit vector there
          p_FA(p_Kld(Irows(irow))) = 1.0_SP
          
          ! and clear the row
          p_FA(p_Kld(Irows(irow))+1:p_Kld(Irows(irow)+1)-1) = 0.0_SP
          
        END DO
        
      CASE DEFAULT
        CALL output_line('Unsupported data format',&
            OU_CLASS_ERROR,OU_MODE_STD,'replaceLines_format7')
        CALL sys_halt()
      END SELECT
      
    END SUBROUTINE replaceLines_format7
    
  END SUBROUTINE mmod_replaceLinesByUnit

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mmod_clearOffdiags (rmatrix,Irows)
  
!<description>
    ! This routine replaces the offdiagonal entries in some lines of 
    ! a given scalar matrix by zero. The diagonal elements are not
    ! changed.
!</description>

!<input>
    ! A list of row numbers of all the rows which are to be replaced
    ! by unit vectors.
    INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: Irows
!</input>

!<inputoutput>
    ! The matrix which is to be modified.
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>

!</subroutine>

    ! At first we must take care of the matrix type.
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9)
      CALL removeOffdiags_format9 (rmatrix,Irows)
    CASE (LSYSSC_MATRIX7)
      CALL removeOffdiags_format7 (rmatrix,Irows)
    END SELECT
    
  CONTAINS
    
    ! ****************************************
    ! The replacement routine for format 9
    
    SUBROUTINE removeOffdiags_format9 (rmatrix,Irows)
      
      INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: Irows
      TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
      
      ! local variables
      INTEGER(PREC_MATIDX) :: irow
      REAL(DP), DIMENSION(:), POINTER :: p_DA
      REAL(SP), DIMENSION(:), POINTER :: p_FA
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_Kdiagonal
      REAL(DP) :: ddiag,fdiag
      
      ! Get Kld and Kdiagonal
      CALL lsyssc_getbase_Kld(rmatrix,p_Kld)
      CALL lsyssc_getbase_Kdiagonal(rmatrix,p_Kdiagonal)
      
      ! Take care of the format of the entries
      SELECT CASE (rmatrix%cdataType)
      CASE (ST_DOUBLE)
        ! Get the data array
        CALL lsyssc_getbase_double(rmatrix,p_DA)
        
        ! loop through the rows
        DO irow = 1,SIZE(Irows)
          
          ! Get the diagonal
          ddiag = p_DA(p_Kdiagonal(Irows(irow)))
          
          ! Clear the row
          p_DA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_DP
          
          ! restore the diagonal
          p_DA(p_Kdiagonal(Irows(irow))) = ddiag
          
        END DO

      CASE (ST_SINGLE)
        ! Get the data array
        CALL lsyssc_getbase_single(rmatrix,p_FA)
        
        ! loop through the rows
        DO irow = 1,SIZE(Irows)
          
          ! Get the diagonal
          fdiag = p_FA(p_Kdiagonal(Irows(irow)))
          
          ! Clear the row
          p_FA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_SP
          
          ! restore the diagonal
          p_FA(p_Kdiagonal(Irows(irow))) = fdiag
          
        END DO
        
      CASE DEFAULT
        CALL output_line('Unsupported data format',&
            OU_CLASS_ERROR,OU_MODE_STD,'removeOffdiags_format9')
        CALL sys_halt()
      END SELECT
      
    END SUBROUTINE removeOffdiags_format9
    
    ! ****************************************
    ! The replacement routine for format 7
    
    SUBROUTINE removeOffdiags_format7 (rmatrix,Irows)
      
      INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: Irows
      TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
      
      ! local variables
      INTEGER(PREC_MATIDX) :: irow
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
      REAL(DP), DIMENSION(:), POINTER :: p_DA
      REAL(SP), DIMENSION(:), POINTER :: p_FA
      
      ! Get Kld:
      CALL lsyssc_getbase_Kld(rmatrix,p_Kld)
      
      ! Take care of the format of the entries
      SELECT CASE (rmatrix%cdataType)
      CASE (ST_DOUBLE)
        ! Get the data array
        CALL lsyssc_getbase_double(rmatrix,p_DA)
        
        ! loop through the rows
        DO irow = 1,SIZE(Irows)
          
          ! Clear the row except for the diagonal
          p_DA(p_Kld(Irows(irow))+1:p_Kld(Irows(irow)+1)-1) = 0.0_DP
          
        END DO

      CASE (ST_SINGLE)
        ! Get the data array
        CALL lsyssc_getbase_single(rmatrix,p_FA)
        
        ! loop through the rows
        DO irow = 1,SIZE(Irows)
          
          ! Clear the row except for the diagonal
          p_FA(p_Kld(Irows(irow))+1:p_Kld(Irows(irow)+1)-1) = 0.0_SP
          
        END DO
        
      CASE DEFAULT
        CALL output_line('Unsupported data format',&
            OU_CLASS_ERROR,OU_MODE_STD,'removeOffdiags_format7')
        CALL sys_halt()
      END SELECT
      
    END SUBROUTINE removeOffdiags_format7
    
  END SUBROUTINE mmod_clearOffdiags
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mmod_replaceLinesByZero (rmatrix,Irows)
  
!<description>
    ! This routine replaces some lines of a given scalar matrix by unit vectors.
!</description>

!<input>
    ! A list of row numbers of all the rows which are to be replaced
    ! by unit vectors.
    INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: Irows
!</input>

!<inputoutput>
    ! The matrix which is to be modified.
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>

!</subroutine>

    ! At first we must take care of the matrix type.
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
      CALL replaceLinesZero_format97 (rmatrix,Irows)
    END SELECT
    
  CONTAINS
    
    ! ****************************************
    ! The replacement routine for format 9 and 7
    
    SUBROUTINE replaceLinesZero_format97 (rmatrix,Irows)
      
      INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: Irows
      TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
    
      ! local variables
      INTEGER(PREC_MATIDX) :: irow
      REAL(DP), DIMENSION(:), POINTER :: p_DA
      REAL(SP), DIMENSION(:), POINTER :: p_FA
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_Kdiagonal
      
      ! Get Kld and Kdiagonal
      CALL lsyssc_getbase_Kld(rmatrix,p_Kld)
      CALL lsyssc_getbase_Kdiagonal(rmatrix,p_Kdiagonal)
      
      ! Take care of the format of the entries
      SELECT CASE (rmatrix%cdataType)
      CASE (ST_DOUBLE)
        ! Get the data array
        CALL lsyssc_getbase_double(rmatrix,p_DA)
        
        ! loop through the rows
        DO irow = 1,SIZE(Irows)
          ! Clear the row
          p_DA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_DP
        END DO

      CASE (ST_SINGLE)
        ! Get the data array
        CALL lsyssc_getbase_single(rmatrix,p_FA)
        
        ! loop through the rows
        DO irow = 1,SIZE(Irows)
          ! Clear the row
          p_FA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_SP
        END DO
        
      CASE DEFAULT
        CALL output_line('Unsupported data format',&
            OU_CLASS_ERROR,OU_MODE_STD,'replaceLinesZero_format97')
        CALL sys_halt()
      END SELECT
      
    END SUBROUTINE replaceLinesZero_format97
    
  END SUBROUTINE mmod_replaceLinesByZero

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE mmod_mergeLines (rmatrix,Irows)
  
!<description>
    ! This routine merges some pairs of lines of a given scalar matrix.
!</description>

!<input>
    ! A list of row numbers of all the rows which are to be merged.
    INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:,:) :: Irows
!</input>

!<inputoutput>
    ! The matrix which is to be modified.
    TYPE(t_matrixScalar), INTENT(INOUT)              :: rmatrix
!</inputoutput>

!</subroutine>

    ! Check, if matrix is not a copy of another matrix or if resize is to be enforced
    IF (IAND(rmatrix%imatrixSpec, LSYSSC_MSPEC_STRUCTUREISCOPY) .NE. 0 .OR.&
        IAND(rmatrix%imatrixSpec, LSYSSC_MSPEC_CONTENTISCOPY)   .NE. 0) THEN
      CALL output_line('A copied matrix cannot be modified!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mmod_mergeLines')
        CALL sys_halt()
    END IF

    ! At first we must take care of the matrix type.
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9, LSYSSC_MATRIX9INTL)
      CALL mergeLines_format9 (rmatrix, Irows)
    CASE (LSYSSC_MATRIX7, LSYSSC_MATRIX7INTL)
      CALL mergeLines_format7 (rmatrix, Irows)
    END SELECT

    print *, "Lines merged"
    STOP

  CONTAINS
    
    ! ****************************************
    ! The merging routine for format 7

    SUBROUTINE mergeLines_format7(rmatrix, Irows)
      TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
      INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:,:) :: Irows
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_ImergeWithRow
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol
      INTEGER(PREC_MATIDX) :: ieq,jeq,ild,jld,irow,jrow,icol,jcol,na,naIncr
      INTEGER              :: h_ImergeWithRow

      print *, Irows
      
      ! Allocate temporal memory
      CALL storage_new('mergeLines_format7', 'p_ImergeWithRow', rmatrix%NEQ, ST_INT,&
          h_ImergeWithRow, ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(h_ImergeWithRow, p_ImergeWithRow)

      ! Get Kld and Kcol
      CALL lsyssc_getbase_Kld(rmatrix, p_Kld)
      CALL lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      
      ! Loop over the list of rows to be merged
      naIncr = 0
      DO ieq = 1, SIZE(Irows,2)
        
        ! Get row numbers to be merged
        irow = Irows(1, ieq)
        jrow = Irows(2, ieq)

        ! If both numbers are identical, then nothing needs to be done
        IF (irow .EQ. jrow) CYCLE

        ! Mark rows for merging
        p_ImergeWithRow(irow) = jrow
        p_ImergeWithRow(jrow) = irow

        ! Treat diagonal entry of row irow separately
        DO jld = p_Kld(jrow)+1, p_Kld(jrow+1)-1
          IF (p_Kcol(jld) .EQ. irow) EXIT
          IF (p_Kcol(jld) .LT. irow) CYCLE
          naIncr = naIncr+1
          EXIT
        END DO
        
        ! Loop over row number irow
        DO ild = p_Kld(irow)+1, p_Kld(irow+1)-1
          icol = p_Kcol(ild)
          DO jld = p_Kld(jrow), p_Kld(jrow+1)-1
            IF (p_Kcol(jld) .EQ. icol) EXIT
            IF (p_Kcol(jld) .LT. icol) EXIT
            naIncr = naIncr+1
            EXIT
          END DO
        END DO

        ! Treat diagonal entry of row jrow separately
        DO ild = p_Kld(irow)+1, p_Kld(irow+1)-1
          IF (p_Kcol(ild) .EQ. jrow) EXIT
          IF (p_Kcol(ild) .LT. jrow) CYCLE
          naIncr = naIncr+1
          EXIT
        END DO
        
        ! Loop over row number jrow
        DO jld = p_Kld(jrow)+1, p_Kld(jrow+1)-1
          jcol = p_Kcol(jld)
          DO ild = p_Kld(irow), p_Kld(irow+1)-1
            IF (p_Kcol(ild) .EQ. jcol) EXIT
            IF (p_Kcol(ild) .LT. jcol) EXIT
            naIncr = naIncr+1
            EXIT
          END DO
        END DO
      END DO

      DO ieq = 1, rmatrix%NEQ
        PRINT *, p_Kcol(p_Kld(ieq):p_Kld(ieq+1)-1)
      END DO
      pause
      print *, "---"

      ! Ok, we now know the number of nonzero entries that need to be added to the
      ! global matrix, hence, set new matix dimension and resize Kcol
      rmatrix%NA = rmatrix%NA+naIncr
      CALL storage_realloc('mergeLines_format7', rmatrix%NA, rmatrix%h_Kcol,&
          ST_NEWBLOCK_NOINIT, .TRUE.)
      CALL lsyssc_getbase_Kcol(rmatrix,p_Kcol)

      ! Initialize current position
      na = rmatrix%NA

      ! Loop over all equations of the matrix but in reverse order (!!!)
      DO ieq = rmatrix%NEQ, 1, -1

        print *, "Processing row",ieq

        ! Ok, so what do we have to do with this row
        IF (p_ImergeWithRow(ieq) .EQ. 0) THEN
          
          ! Set current position
          naIncr = na

          ! Just copy the row without modifications
          DO ild = p_Kld(ieq+1)-1, p_Kld(ieq), -1
            p_Kcol(naIncr) = p_Kcol(ild)
            naIncr         = naIncr-1
          END DO
          
          ! Update Kld for next column and reduce number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr
          
        ELSEIF (p_ImergeWithRow(ieq) .GT. ieq) THEN

          ! The two rows have already been merged, hence 
          ! we can adopt all entries from the indicated row
          jeq    = p_ImergeWithRow(ieq)
          naIncr = na

          ! Just copy the row without modifications
          DO jld = p_Kld(jeq+1)-1, p_Kld(jeq), -1
            IF (p_Kcol(jld) .EQ. ieq) CYCLE
            p_Kcol(naIncr) = p_Kcol(jld)
            naIncr         = naIncr-1
          END DO

          ! Sort row according to format 7 format
          DO ild = naIncr, na
            IF (p_Kcol(ild) .GT. jeq) EXIT
            p_Kcol(ild-1) = p_Kcol(ild)
            p_Kcol(ild)   = jeq
          END DO
          p_Kcol(naIncr) = ieq
          naIncr = naIncr-1
          
          ! Update Kld for next column and reduce number of nonzero entries
          p_Kld(ieq+1) = na+1
          na = naIncr

        ELSE
          
          ! We actually have to merge the two rows
          jeq = p_ImergeWithRow(ieq)

          ! First, sort the row that should b merged into the current row in ascending order.
          ! This is mandatory since the matrix is stored in format 7 so that the diagonal
          ! entry is at the first position and not stored in-between.
          DO jld = p_Kld(jeq)+1, p_Kld(jeq+1)-1
            IF (p_Kcol(jld) .GT. jeq) EXIT
            p_Kcol(jld-1) = p_Kcol(jld)
            p_Kcol(jld)   = jeq
          END DO
          
          ! Ok, now the row is in ascending order
          naIncr = na
          ild = p_Kld(ieq+1)-1
          jld = p_Kld(jeq+1)-1

          ! Loop in reverse order until both rows have been processed
          DO WHILE((ild .GT. p_Kld(ieq)) .AND. (jld .GE. p_Kld(jeq)))

            ! Get column indices for both rows and recude position counter by one
            icol = p_Kcol(ild)
            jcol = p_Kcol(jld)

            ! Which is the next column that should be processed
            IF (jcol .GT. icol) THEN
              p_Kcol(naIncr) = jcol
              jld = jld-1
            ELSEIF (jcol .EQ. icol) THEN
              p_Kcol(naIncr) = jcol
              ild = ild-1
              jld = jld-1
            ELSE
              p_Kcol(naIncr) = icol
              ild = ild-1
            END IF

            ! Decrease position counter
            naIncr = naIncr-1
          END DO

          ! Copy leftover from the row that corresponds to JEQ?
          DO WHILE (jld .GE. p_Kld(jeq))
            p_Kcol(naIncr) = p_Kcol(jld)
            naIncr = naIncr-1
            jld    = jld-1
          END DO

          ! Copy leftover from the row that corresponds to IEQ?
          DO WHILE (ild .GE. p_Kld(ieq))
            p_Kcol(naIncr) = p_Kcol(ild)
            naIncr = naIncr-1
            ild    = ild-1
          END DO

          ! That's it. ! Update Kld for next column and adjust number of nonzero entries
          p_Kld(ieq+1) = na+1
          na = naIncr
        END IF
      END DO

      p_Kld(1) = na

      DO ieq = 1, rmatrix%NEQ
        PRINT *, p_Kcol(p_Kld(ieq):p_Kld(ieq+1)-1)
      END DO
      
      stop

    END SUBROUTINE mergeLines_format7

    
    ! ****************************************
    ! The merging routine for format 9

    SUBROUTINE mergeLines_format9(rmatrix, Irows)
      TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
      INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:,:) :: Irows
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_ImergeWithRow
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol
      INTEGER(PREC_MATIDX) :: ieq,jeq,ild,jld,irow,jrow,icol,jcol,na,naIncr
      INTEGER              :: h_ImergeWithRow

      print *, Irows
      
      ! Allocate temporal memory
      CALL storage_new('mergeLines_format7', 'p_ImergeWithRow', rmatrix%NEQ, ST_INT,&
          h_ImergeWithRow, ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(h_ImergeWithRow, p_ImergeWithRow)

      ! Get Kld and Kcol
      CALL lsyssc_getbase_Kld(rmatrix, p_Kld)
      CALL lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      CALL lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
      
      ! Loop over the list of rows to be merged
      naIncr = 0
      DO ieq = 1, SIZE(Irows,2)
        
        ! Get row numbers to be merged
        irow = Irows(1, ieq)
        jrow = Irows(2, ieq)

        ! If both numbers are identical, then nothing needs to be done
        IF (irow .EQ. jrow) CYCLE

        ! Mark rows for merging
        p_ImergeWithRow(irow) = jrow
        p_ImergeWithRow(jrow) = irow

        ! Loop over row number irow
        DO ild = p_Kld(irow), p_Kld(irow+1)-1
          icol = p_Kcol(ild)
          DO jld = p_Kld(jrow), p_Kld(jrow+1)-1
            IF (p_Kcol(jld) .EQ. icol) EXIT
            IF (p_Kcol(jld) .LT. icol) EXIT
            naIncr = naIncr+1
            EXIT
          END DO
        END DO
        
        ! Loop over row number jrow
        DO jld = p_Kld(jrow), p_Kld(jrow+1)-1
          jcol = p_Kcol(jld)
          DO ild = p_Kld(irow), p_Kld(irow+1)-1
            IF (p_Kcol(ild) .EQ. jcol) EXIT
            IF (p_Kcol(ild) .LT. jcol) EXIT
            naIncr = naIncr+1
            EXIT
          END DO
        END DO
      END DO

      ! Ok, we now know the number of nonzero entries that need to be added to the
      ! global matrix, hence, set new matix dimension and resize Kcol
      rmatrix%NA = rmatrix%NA+naIncr
      CALL storage_realloc('mergeLines_format9', rmatrix%NA, rmatrix%h_Kcol,&
          ST_NEWBLOCK_NOINIT, .TRUE.)
      CALL lsyssc_getbase_Kcol(rmatrix,p_Kcol)

      ! Initialize current position
      na = rmatrix%NA

      ! Loop over all equations of the matrix but in reverse order (!!!)
      DO ieq = rmatrix%NEQ, 1, -1

        ! Ok, so what do we have to do with this row
        IF (p_ImergeWithRow(ieq) .EQ. 0) THEN
          
          ! Set current position
          naIncr = na

          ! Just copy the row without modifications
          DO ild = p_Kld(ieq+1)-1, p_Kld(ieq), -1
            p_Kcol(naIncr) = p_Kcol(ild)
            naIncr         = naIncr-1
          END DO
          
          ! Update Kld for next column and adjust number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr
          
        ELSEIF (p_ImergeWithRow(ieq) .GT. ieq) THEN

          ! The two rows have already been merged, hence 
          ! we can adopt all entries from the indicated row
          jeq    = p_ImergeWithRow(ieq)
          naIncr = na

          ! Just copy the row without modifications
          DO jld = p_Kld(jeq+1)-1, p_Kld(jeq), -1
            p_Kcol(naIncr) = p_Kcol(jld)
            naIncr         = naIncr-1
          END DO
          
          ! Update Kld for next column and adjust number of nonzero entries
          p_Kld(ieq+1) = na+1
          na = naIncr

        ELSE
          
          ! We actually have to merge the two rows
          jeq = p_ImergeWithRow(ieq)

          ! The rows are in ascendin order
          naIncr = na
          ild = p_Kld(ieq+1)-1
          jld = p_Kld(jeq+1)-1

          ! Loop in reverse order until both rows have been processed
          DO WHILE((ild .GE. p_Kld(ieq)) .AND. (jld .GE. p_Kld(jeq)))

            ! Get column indices for both rows and recude position counter by one
            icol = p_Kcol(ild)
            jcol = p_Kcol(jld)

            ! Which is the next column that should be processed
            IF (jcol .GT. icol) THEN
              p_Kcol(naIncr) = jcol
              jld = jld-1
            ELSEIF (jcol .EQ. icol) THEN
              p_Kcol(naIncr) = jcol
              ild = ild-1
              jld = jld-1
            ELSE
              p_Kcol(naIncr) = icol
              ild = ild-1
            END IF

            ! Decrease position counter
            naIncr = naIncr-1
          END DO

          ! Copy leftover from the row that corresponds to JEQ?
          DO WHILE (jld .GE. p_Kld(jeq))
            p_Kcol(naIncr) = p_Kcol(jld)
            naIncr = naIncr-1
            jld    = jld-1
          END DO

          ! Copy leftover from the row that corresponds to IEQ?
          DO WHILE (ild .GE. p_Kld(ieq))
            p_Kcol(naIncr) = p_Kcol(ild)
            naIncr = naIncr-1
            ild    = ild-1
          END DO

          ! That's it. ! Update Kld for next column and adjust number of nonzero entries
          p_Kld(ieq+1) = na+1
          na = naIncr
        END IF
      END DO

      p_Kld(1) = na

      DO ieq = 1, rmatrix%NEQ
        PRINT *, p_Kcol(p_Kld(ieq):p_Kld(ieq+1)-1)
      END DO
      
      stop

    END SUBROUTINE mergeLines_format9
  END SUBROUTINE mmod_mergeLines

END MODULE matrixmodification

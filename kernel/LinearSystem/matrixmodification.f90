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

  SUBROUTINE mmod_mergeLines (rmatrix,Irows,bsymmetric)
  
!<description>
    ! This routine merges some pairs of lines of a given scalar matrix.
    ! Note that the data (if any) will be cleared. If the optional parameter
    ! bsymmetric=.TRUE., then the sparsity pattern will be symmetric.
!</description>

!<input>
    ! A list of row numbers of all the rows which are to be merged.
    INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:,:) :: Irows
    
    ! OPTIONAL: If bsymmetric=.TRUE. the sparsity pattern will be symmetric
    LOGICAL, INTENT(IN), OPTIONAL                    :: bsymmetric
!</input>

!<inputoutput>
    ! The matrix which is to be modified.
    TYPE(t_matrixScalar), INTENT(INOUT)              :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables
     INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_ImergeWithRow
     INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol
     INTEGER(PREC_MATIDX) :: ieq,ild,jld,irow,jrow,icol,jcol,naIncr
     INTEGER              :: h_ImergeWithRow

    ! Check, if matrix is not a copy of another matrix or if resize is to be enforced
    IF (IAND(rmatrix%imatrixSpec, LSYSSC_MSPEC_STRUCTUREISCOPY) .NE. 0 .OR.&
        IAND(rmatrix%imatrixSpec, LSYSSC_MSPEC_CONTENTISCOPY)   .NE. 0) THEN
      CALL output_line('A copied matrix cannot be modified!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mmod_mergeLines')
        CALL sys_halt()
    END IF
    
    ! Allocate temporal memory
    CALL storage_new('mmod_mergeLines', 'p_ImergeWithRow', rmatrix%NEQ, ST_INT,&
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
      
      ! Loop over row number irow
      loop1: DO ild = p_Kld(irow), p_Kld(irow+1)-1
        icol = p_Kcol(ild)
        loop2: DO jld = p_Kld(jrow), p_Kld(jrow+1)-1
          IF (p_Kcol(jld) .EQ. icol) CYCLE loop1
        END DO loop2
        naIncr = naIncr+1
      END DO loop1
      
      ! Loop over row number jrow
      loop3: DO jld = p_Kld(jrow), p_Kld(jrow+1)-1
        jcol = p_Kcol(jld)
        loop4: DO ild = p_Kld(irow), p_Kld(irow+1)-1
          IF (p_Kcol(ild) .EQ. jcol) CYCLE loop3
        END DO loop4
        naIncr = naIncr+1
      END DO loop3
    END DO
    
    ! Ok, we now know the number of nonzero entries that need to be added to the
    ! global matrix, hence, set new matix dimension and resize Kcol
    rmatrix%NA = rmatrix%NA+naIncr
    CALL storage_realloc('mmod_mergeLines', rmatrix%NA, rmatrix%h_Kcol,&
        ST_NEWBLOCK_NOINIT, .TRUE.)
    
    ! At first we must take care of the matrix type and modify the structure
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9, LSYSSC_MATRIX9INTL)
      CALL mergeLines_format9 (rmatrix, p_ImergeWithRow)
      
      ! Do we have to generate a symmetric sparsity graph?
      IF (PRESENT(bsymmetric)) THEN
        IF (bsymmetric) CALL mergeColumns_format9(rmatrix, p_ImergeWithRow)
      END IF
      
    CASE (LSYSSC_MATRIX7, LSYSSC_MATRIX7INTL)
      CALL mergeLines_format7 (rmatrix, p_ImergeWithRow)
      
      ! Do we have to generate a symmetric sparsity graph?
      IF (PRESENT(bsymmetric)) THEN
        IF (bsymmetric) CALL mergeColumns_format7(rmatrix, p_ImergeWithRow)
      END IF

    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'mmod_mergeLines')
    END SELECT

    ! Free temporal storage
    CALL storage_free(h_ImergeWithRow)

    ! Next, resize the data accordingly (if present)
    IF (rmatrix%h_DA .NE. ST_NOHANDLE) THEN
      SELECT CASE (rmatrix%cinterleavematrixFormat)
      CASE (LSYSSC_MATRIXUNDEFINED)
        CALL storage_realloc('mmod_mergeLines', rmatrix%NA, rmatrix%h_DA,&
            ST_NEWBLOCK_ZERO, .FALSE.)
        
      CASE (LSYSSC_MATRIXD)
        CALL storage_realloc('mmod_mergeLines', rmatrix%NA*rmatrix%NVAR,&
            rmatrix%h_DA, ST_NEWBLOCK_ZERO, .FALSE.)
        
      CASE (LSYSSC_MATRIX1)
        CALL storage_realloc('mmod_mergeLines', rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR,&
            rmatrix%h_DA, ST_NEWBLOCK_ZERO, .FALSE.)
        
      CASE DEFAULT
        CALL output_line('Unsupported interleave matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mmod_mergeLines')
      END SELECT
    END IF

  CONTAINS
    
    ! ****************************************
    ! The row merging routine for format 7

    SUBROUTINE mergeLines_format7(rmatrix, ImergeWithRow)
      TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrix
      INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: ImergeWithRow
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol
      INTEGER(PREC_MATIDX) :: ieq,jeq,ild,jld,irow,jrow,icol,jcol,na,naIncr

      ! Get Kld and Kcol
      CALL lsyssc_getbase_Kld(rmatrix, p_Kld)
      CALL lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      
      ! Initialize current position
      na = rmatrix%NA

      ! Loop over all equations of the matrix but in reverse order (!!!)
      DO ieq = rmatrix%NEQ, 1, -1

        ! Ok, so what do we have to do with this row?
        IF (ImergeWithRow(ieq) .EQ. 0) THEN
          
          ! Set current position
          naIncr = na

          ! Just copy the row without modifications
          DO ild = p_Kld(ieq+1)-1, p_Kld(ieq), -1
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(ild)

            ! Decrease position counter
            naIncr = naIncr-1
          END DO
          
          ! Update Kld for next column and reduce number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr
          
        ELSEIF (ImergeWithRow(ieq) .GT. ieq) THEN

          ! The two rows have already been merged, hence 
          ! we can adopt all entries from the indicated row
          jeq = ImergeWithRow(ieq)

          ! Set current position
          naIncr = na

          ! Just copy the row without modifications
          DO jld = p_Kld(jeq+1)-1, p_Kld(jeq), -1
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(jld)
            
            ! Decrease position counter
            naIncr = naIncr-1
          END DO

          ! Sort row according to format 7 format. First we need to find the position of the diagonal entry
          ! and swap the diagonal entry from the copied row with the entry of the current row. Then we need
          ! to move the diagonal entry from the copied row to its correct position.
          sort: DO ild = naIncr+1, na
            IF (p_Kcol(ild) .EQ. ieq) THEN
              ! Swap the first entry which corresponds to the diagonal entry of the copied row with current position
              p_Kcol(ild)      = jeq
              p_Kcol(naIncr+1) = ieq
              
              ! Move the swapped diagonal entry from the copied row to its correct position
              DO jld = ild+1, na
                IF (p_Kcol(jld) .GT. jeq) EXIT sort
                p_Kcol(jld-1) = p_Kcol(jld)
                p_Kcol(jld)   = jeq
              END DO
              EXIT sort
            END IF
          END DO sort
          
          ! Update Kld for next column and reduce number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr

        ELSE
          
          ! We actually have to merge the two rows
          jeq = ImergeWithRow(ieq)

          ! First, sort the row that should be merged into the current row in ascending order.
          ! This is mandatory since the matrix is stored in format 7 so that the diagonal
          ! entry is at the first position and not stored in-between.
          DO jld = p_Kld(jeq)+1, p_Kld(jeq+1)-1
            IF (p_Kcol(jld) .GT. jeq) EXIT
            p_Kcol(jld-1) = p_Kcol(jld)
            p_Kcol(jld)   = jeq
          END DO
          
          ! Ok, now the row is in ascending order
          ild = p_Kld(ieq+1)-1
          jld = p_Kld(jeq+1)-1

          ! Set current position
          naIncr = na

          ! Loop in reverse order until both rows have been processed
          ! Note that the diagonal entry of row irow is treated below
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
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(jld)

            ! Decrease position counter
            naIncr = naIncr-1
            jld    = jld-1
          END DO
          
          ! Copy leftover from the row that corresponds to IEQ?
          ! Here, we explicitly copy the diagonal entry to the first position.
          DO WHILE (ild .GE. p_Kld(ieq))
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(ild)
            
            ! Decrease position counter
            naIncr = naIncr-1
            ild    = ild-1
          END DO

          ! That's it! Update Kld for next column and adjust number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr
        END IF

        ! Ok, if we process two consecutive rows, then we must also adjust the starting position
        IF ((ImergeWithRow(ieq) .EQ. ieq-1) .AND. (ieq .NE. 1)) p_Kld(ieq) = na+1
      END DO

      ! Consistency check. If na is not zero then something went wrong
      IF (na .NE. 0) THEN
        CALL output_line('An internal error occured; please contact the Featflow team!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mergeLines_format7')
        CALL sys_halt()
      END IF
    END SUBROUTINE mergeLines_format7

    
    ! ****************************************
    ! The row merging routine for format 9

    SUBROUTINE mergeLines_format9(rmatrix, ImergeWithRow)
      TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrix
      INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: ImergeWithRow
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_Kdiagonal
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol
      INTEGER(PREC_MATIDX) :: ieq,jeq,ild,jld,irow,jrow,icol,jcol,na,naIncr
      
      ! Get Kld, Kcol and Kdiagonal
      CALL lsyssc_getbase_Kld(rmatrix, p_Kld)
      CALL lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      CALL lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

      ! Initialize current position
      na = rmatrix%NA

      ! Loop over all equations of the matrix but in reverse order (!!!)
      DO ieq = rmatrix%NEQ, 1, -1
        
        ! Ok, so what do we have to do with this row?
        IF (ImergeWithRow(ieq) .EQ. 0) THEN
          
          ! Set current position
          naIncr = na
          
          ! Just copy the row without modifications
          DO ild = p_Kld(ieq+1)-1, p_Kld(ieq), -1
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(ild)

            ! Check for diagonal entry
            IF (p_Kcol(naIncr) .EQ. ieq) p_Kdiagonal(ieq) = naIncr
            
            ! Decrease position counter
            naIncr = naIncr-1
          END DO
          
          ! Update Kld for next column and adjust number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr
          
        ELSEIF (ImergeWithRow(ieq) .GT. ieq) THEN

          ! The two rows have already been merged, hence 
          ! we can adopt all entries from the indicated row
          jeq = ImergeWithRow(ieq)

          ! Set current position
          naIncr = na

          ! Just copy the row without modifications
          DO jld = p_Kld(jeq+1)-1, p_Kld(jeq), -1
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(jld)
            
            ! Check for diagonal entry
            IF (p_Kcol(naIncr) .EQ. ieq) p_Kdiagonal(ieq) = naIncr

            ! Decrease position counter
            naIncr = naIncr-1
          END DO
          
          ! Update Kld for next column and adjust number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr

        ELSE
          
          ! We actually have to merge the two rows
          jeq = ImergeWithRow(ieq)

          ! The rows are in ascending order
          ild = p_Kld(ieq+1)-1
          jld = p_Kld(jeq+1)-1

          ! Set current position
          naIncr = na

          ! Loop in reverse order until both rows have been processed
          DO WHILE((ild .GE. p_Kld(ieq)) .AND. (jld .GE. p_Kld(jeq)))

            ! Get column indices for both rows and recude position counter by one
            icol = p_Kcol(ild)
            jcol = p_Kcol(jld)

            ! Which is the next column that should be processed?
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

            ! Check for diagonal entry
            IF (p_Kcol(naIncr) .EQ. ieq) p_Kdiagonal(ieq) = naIncr

            ! Decrease position counter
            naIncr = naIncr-1
          END DO

          ! Copy leftover from the row that corresponds to JEQ?
          DO WHILE (jld .GE. p_Kld(jeq))
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(jld)

            ! Check for diagonal entry
            IF (p_Kcol(naIncr) .EQ. ieq) p_Kdiagonal(ieq) = naIncr
            
            ! Decrease position counter
            naIncr = naIncr-1
            jld    = jld-1
          END DO

          ! Copy leftover from the row that corresponds to IEQ?
          DO WHILE (ild .GE. p_Kld(ieq))
            ! Copy column index
            p_Kcol(naIncr) = p_Kcol(ild)

            ! Check for diagonal entry
            IF (p_Kcol(naIncr) .EQ. ieq) p_Kdiagonal(ieq) = naIncr

            ! Decrease position counter
            naIncr = naIncr-1
            ild    = ild-1
          END DO

          ! That's it! Update Kld for next column and adjust number of nonzero entries
          p_Kld(ieq+1) = na+1
          na           = naIncr
        END IF

        ! Ok, if we process two consecutive rows, then we must also adjust the starting position
        IF ((ImergeWithRow(ieq) .EQ. ieq-1) .AND. (ieq .NE. 1)) p_Kld(ieq) = na+1
      END DO

      ! Consistency check. If na is not zero then something went wrong
      IF (na .NE. 0) THEN
        CALL output_line('An internal error occured; please contact the Featflow team!',&
            OU_CLASS_ERROR,OU_MODE_STD,'mergeLines_format9')
        CALL sys_halt()
      END IF
    END SUBROUTINE mergeLines_format9


    ! ****************************************
    ! The column merging routine for format 7

    SUBROUTINE mergeColumns_format7(rmatrix, ImergeWithRow)
      TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrix
      INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: ImergeWithRow
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KldAux
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol,p_KcolAux
      INTEGER(PREC_MATIDX) :: ieq,ild,jld,jjld,icol,jcol,na,naIncr,idxIncr,iidxIncr
      INTEGER              :: h_KldAux,h_KcolAux

      ! Allocate temporal memory
      CALL storage_new('mergeColumns_format7', 'p_KldAux', rmatrix%NEQ+1, ST_INT,&
          h_KldAux, ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_int(h_KldAux, p_KldAux)

      ! Get Kld and Kcol
      CALL lsyssc_getbase_Kld(rmatrix, p_Kld)
      CALL lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      ! Compute number of nonzero entries that need to be inserted into the matrix.
      ! First, compute the number of nonzero entries present in each row.
      DO ieq = 1, rmatrix%NEQ
        p_KldAux(ieq) = p_Kld(ieq+1)-p_Kld(ieq)
      END DO

      ! Next, subtract the number of nonzero entries present in each column. 
      DO ild = 1, rmatrix%NA
        icol = p_Kcol(ild)
        p_KldAux(icol) = p_KldAux(icol)-1
      END DO

      ! If an entry is zero, then the number of nonzero entries in the corresponding
      ! row equals the number of nonzero entries in the corresponding column and
      ! nothing needs to be done. A negative entry indicates that the number of 
      ! nonzero entries in the corresponding column exceeds the number of nonzero
      ! entries in the corresponding column. Hence, we must fill the missing entries.
      naIncr = 0
      DO ieq = 1, rmatrix%NEQ
        p_KldAux(ieq) = -MIN(p_KldAux(ieq), 0)
        naIncr        =  naIncr+p_KldAux(ieq)
      END DO

      ! Return of no new entries need to be considered
      IF (naIncr .EQ. 0) THEN
        CALL storage_free(h_KldAux)
        RETURN
      END IF

      ! Make a copy of the original column indices
      h_KcolAux = ST_NOHANDLE
      CALL storage_copy(rmatrix%h_Kcol, h_KcolAux)
      CALL storage_getbase_int(h_KcolAux, p_KcolAux)

      ! Ok, we now know the number of nonzero entries that need to be added to the
      ! global matrix, hence, set new matix dimension and resize Kcol
      rmatrix%NA = rmatrix%NA+naIncr
      CALL storage_realloc('mergeLines_format7', rmatrix%NA, rmatrix%h_Kcol,&
          ST_NEWBLOCK_NOINIT, .TRUE.)
      CALL lsyssc_getbase_Kcol(rmatrix,p_Kcol)

      ! Adjust the row separator and copy its original content to p_KldAux
      idxIncr     = p_KldAux(1)+p_Kld(2)-p_Kld(1)
      p_KldAux(1) = p_Kld(1)

      DO ieq = 2, rmatrix%NEQ
        ! Compute number of nonzero entries in subsequent row
        iidxIncr = p_KldAux(ieq)+p_Kld(ieq+1)-p_Kld(ieq)

        ! Copy original entry
        p_KldAux(ieq) = p_Kld(ieq)

        ! Adjust absolut position for current row
        p_Kld(ieq) = p_Kld(ieq-1)+idxIncr

        ! Swap increment
        idxIncr = iidxIncr
      END DO
      p_KldAux(rmatrix%NEQ+1) = p_Kld(rmatrix%NEQ+1)
      p_Kld(rmatrix%NEQ+1)    = p_Kld(rmatrix%NEQ)+idxIncr
      
      ! In the first step, move the column data to their new positions
      DO ieq = rmatrix%NEQ, 1, -1
        ! Compute offset to new position
        idxIncr = p_Kld(ieq)-p_KldAux(ieq)

        ! Copy data
        IF (idxIncr .GT. 0) THEN
          DO ild = p_KldAux(ieq+1)-1, p_KldAux(ieq), -1
            p_Kcol(ild+idxIncr) = p_Kcol(ild)
          END DO
        END IF
      END DO

      ! Loop over all matrix rows
      loop1: DO ieq = 1, rmatrix%NEQ
        
        ! If the row has not been merged then skip it
        IF (ImergeWithRow(ieq) .EQ. 0) CYCLE
        
        ! Otherwise, we process all entries of the current row
        ! and check if the corresponding column entries exist
        loop2: DO ild = p_KldAux(ieq)+1, p_KldAux(ieq+1)-1

          ! Get column number
          icol = p_KcolAux(ild)

          ! Skip diagonal entries
          IF (icol .EQ. ieq) CYCLE loop2
          
          ! Loop over row number icol and check if entry ieq exists
          loop3: DO jld = p_Kld(icol)+1, p_Kld(icol+1)-1

            ! Get column number
            jcol = p_Kcol(jld)

            ! Did we find the entry?
            IF (jcol .EQ. ieq) EXIT

            ! Skip left off-diagonal entries
            IF (jcol .LT. ieq) CYCLE
            
            ! No, the entry was not found
            DO jjld = p_Kld(icol+1)-1, jld+1, -1
              p_Kcol(jjld) = p_Kcol(jjld-1)
            END DO
            p_Kcol(jld) = ieq
            EXIT
          END DO loop3
        END DO loop2
      END DO loop1
      
      ! Free temporal storage
      CALL storage_free(h_KldAux)
      CALL storage_free(h_KcolAux)      
    END SUBROUTINE mergeColumns_format7

    
    ! ****************************************
    ! The column merging routine for format 9

    SUBROUTINE mergeColumns_format9(rmatrix, ImergeWithRow)
      TYPE(t_matrixScalar), INTENT(INOUT)            :: rmatrix
      INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: ImergeWithRow
      
      ! local variables
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_Kdiagonal
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KldAux
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol,p_KcolAux
      INTEGER(PREC_MATIDX) :: ieq,ild,jld,jjld,icol,jcol,na,naIncr,idxIncr,iidxIncr
      INTEGER              :: h_KldAux,h_KcolAux

      ! Allocate temporal memory
      CALL storage_new('mergeColumns_format9', 'p_KldAux', rmatrix%NEQ+1, ST_INT,&
          h_KldAux, ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_int(h_KldAux, p_KldAux)

      ! Get Kld and Kcol
      CALL lsyssc_getbase_Kld(rmatrix, p_Kld)
      CALL lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      ! Compute number of nonzero entries that need to be inserted into the matrix.
      ! First, compute the number of nonzero entries present in each row.
      DO ieq = 1, rmatrix%NEQ
        p_KldAux(ieq) = p_Kld(ieq+1)-p_Kld(ieq)
      END DO

      ! Next, subtract the number of nonzero entries present in each column. 
      DO ild = 1, rmatrix%NA
        icol = p_Kcol(ild)
        p_KldAux(icol) = p_KldAux(icol)-1
      END DO

      ! If an entry is zero, then the number of nonzero entries in the corresponding
      ! row equals the number of nonzero entries in the corresponding column and
      ! nothing needs to be done. A negative entry indicates that the number of 
      ! nonzero entries in the corresponding column exceeds the number of nonzero
      ! entries in the corresponding column. Hence, we must fill the missing entries.
      naIncr = 0
      DO ieq = 1, rmatrix%NEQ
        p_KldAux(ieq) = -MIN(p_KldAux(ieq), 0)
        naIncr        =  naIncr+p_KldAux(ieq)
      END DO

      ! Return of no new entries need to be considered
      IF (naIncr .EQ. 0) THEN
        CALL storage_free(h_KldAux)
        RETURN
      END IF

      ! Make a copy of the original column indices
      h_KcolAux = ST_NOHANDLE
      CALL storage_copy(rmatrix%h_Kcol, h_KcolAux)
      CALL storage_getbase_int(h_KcolAux, p_KcolAux)

      ! Ok, we now know the number of nonzero entries that need to be added to the
      ! global matrix, hence, set new matix dimension and resize Kcol
      rmatrix%NA = rmatrix%NA+naIncr
      CALL storage_realloc('mergeLines_format9', rmatrix%NA, rmatrix%h_Kcol,&
          ST_NEWBLOCK_NOINIT, .TRUE.)
      CALL lsyssc_getbase_Kcol(rmatrix,p_Kcol)

      ! Adjust the row separator and copy its original content to p_KldAux
      idxIncr     = p_KldAux(1)+p_Kld(2)-p_Kld(1)
      p_KldAux(1) = p_Kld(1)

      DO ieq = 2, rmatrix%NEQ
        ! Compute number of nonzero entries in subsequent row
        iidxIncr = p_KldAux(ieq)+p_Kld(ieq+1)-p_Kld(ieq)

        ! Copy original entry
        p_KldAux(ieq) = p_Kld(ieq)

        ! Adjust absolut position for current row
        p_Kld(ieq) = p_Kld(ieq-1)+idxIncr

        ! Swap increment
        idxIncr = iidxIncr
      END DO
      p_KldAux(rmatrix%NEQ+1) = p_Kld(rmatrix%NEQ+1)
      p_Kld(rmatrix%NEQ+1)    = p_Kld(rmatrix%NEQ)+idxIncr
      
      ! In the first step, move the column data to their new positions
      DO ieq = rmatrix%NEQ, 1, -1
        ! Compute offset to new position
        idxIncr = p_Kld(ieq)-p_KldAux(ieq)

        ! Copy data
        IF (idxIncr .GT. 0) THEN
          DO ild = p_KldAux(ieq+1)-1, p_KldAux(ieq), -1
            p_Kcol(ild+idxIncr) = p_Kcol(ild)
          END DO
        END IF
      END DO

      ! Loop over all matrix rows
      loop1: DO ieq = 1, rmatrix%NEQ
        
        ! If the row has not been merged then skip it
        IF (ImergeWithRow(ieq) .EQ. 0) CYCLE
        
        ! Otherwise, we process all entries of the current row
        ! and check if the corresponding column entries exist
        loop2: DO ild = p_KldAux(ieq), p_KldAux(ieq+1)-1

          ! Get column number
          icol = p_KcolAux(ild)

          ! Skip diagonal entries
          IF (icol .EQ. ieq) CYCLE loop2
          
          ! Loop over row number icol and check if entry ieq exists
          loop3: DO jld = p_Kld(icol), p_Kld(icol+1)-1

            ! Get column number
            jcol = p_Kcol(jld)

            ! Did we find the entry?
            IF (jcol .EQ. ieq) EXIT

            ! Skip left off-diagonal entries
            IF (jcol .LT. ieq) CYCLE
            
            ! No, the entry was not found
            DO jjld = p_Kld(icol+1)-1, jld+1, -1
              p_Kcol(jjld) = p_Kcol(jjld-1)
            END DO
            p_Kcol(jld) = ieq
            EXIT
          END DO loop3
        END DO loop2
      END DO loop1
      
      ! Rebuild the diagonal array
      CALL lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)
      CALL lsyssc_rebuildKdiagonal(p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)

      ! Free temporal storage
      CALL storage_free(h_KldAux)
      CALL storage_free(h_KcolAux)      
    END SUBROUTINE mergeColumns_format9
  END SUBROUTINE mmod_mergeLines

END MODULE matrixmodification

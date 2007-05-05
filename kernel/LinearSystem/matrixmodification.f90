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
!# </purpose>
!##############################################################################

MODULE matrixmodification

  USE fsystem
  USE storage
  USE linearsystemscalar
  
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
        PRINT *,'replaceLines_format9: Matrix has no data!'
        STOP
      END IF
      
      ! loop through the rows
      DO irow = 1,SIZE(Irows)
      
        ! Clear the row
        p_DA(p_Kld(Irows(irow)):p_Kld(Irows(irow)+1)-1) = 0.0_DP
        
        ! And put a unit vector there
        p_DA(p_Kdiagonal(Irows(irow))) = 1.0_DP
      
      END DO
      
    CASE DEFAULT
      PRINT *,'mmod_replaceLinesByUnit: Only double prec. matices supported!'
      STOP
    END SELECT
    
    END SUBROUTINE

    ! ****************************************
    ! The replacement routine for format 7
    
    SUBROUTINE replaceLines_format7 (rmatrix,Irows)
    
    INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: Irows
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
    
    ! local variables
    INTEGER(PREC_MATIDX) :: irow
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
    REAL(DP), DIMENSION(:), POINTER :: p_DA

    ! Get Kld:
    CALL lsyssc_getbase_Kld(rmatrix,p_Kld)

    ! Take care of the format of the entries
    SELECT CASE (rmatrix%cdataType)
    CASE (ST_DOUBLE)
      ! Get the data array
      CALL lsyssc_getbase_double(rmatrix,p_DA)
      IF (.NOT. ASSOCIATED(p_DA)) THEN
        PRINT *,'replaceLines_format7: Matrix has no data!'
        STOP
      END IF
      
      ! loop through the rows
      DO irow = 1,SIZE(Irows)
      
        ! Put a unit vector there
        p_DA(p_Kld(Irows(irow))) = 1.0_DP

        ! and clear the row
        p_DA(p_Kld(Irows(irow))+1:p_Kld(Irows(irow)+1)-1) = 0.0_DP
      
      END DO
      
    CASE DEFAULT
      PRINT *,'mmod_replaceLinesByUnit: Only double prec. matices supported!'
      STOP
    END SELECT

    END SUBROUTINE
  
  END SUBROUTINE

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
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld,p_Kdiagonal
    REAL(DP) :: ddiag
    
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
      
    CASE DEFAULT
      PRINT *,'mmod_clearOffdiags: Only double prec. matices supported!'
      STOP
    END SELECT
    
    END SUBROUTINE

    ! ****************************************
    ! The replacement routine for format 7
    
    SUBROUTINE removeOffdiags_format7 (rmatrix,Irows)
    
    INTEGER(PREC_MATIDX), INTENT(IN), DIMENSION(:) :: Irows
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
    
    ! local variables
    INTEGER(PREC_MATIDX) :: irow
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
    REAL(DP), DIMENSION(:), POINTER :: p_DA

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
      
    CASE DEFAULT
      PRINT *,'mmod_clearOffdiags: Only double prec. matices supported!'
      STOP
    END SELECT

    END SUBROUTINE
  
  END SUBROUTINE

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
      
    CASE DEFAULT
      PRINT *,'mmod_replaceLinesByZero: Only double prec. matices supported!'
      STOP
    END SELECT
    
    END SUBROUTINE

  END SUBROUTINE

END MODULE

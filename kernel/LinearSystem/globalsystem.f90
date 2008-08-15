!##############################################################################
!# ****************************************************************************
!# <name> globalsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to assembla a global 1x1 system matrix from
!# a block matrix containing multiple submatrices.
!#
!# The following routines can be found here:
!#
!# 1.) glsys_assembleGlobal
!#     -> Assemble a global matrix.
!# </purpose>
!##############################################################################

MODULE globalsystem

  USE fsystem
  USE linearsystemscalar
  USE linearsystemblock
  
  IMPLICIT NONE

CONTAINS
 
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE glsys_assembleGlobal (rsourceMatrix,rdestMatrix, &
                                   bstructure, bcontent, &
                                   cmatrixFormat, cdataType)
  
!<description>
  ! This routine assembles a 1x1 block matrix rdestMatrix from a
  ! NxM block matrix rsourceMatrix.
!</description>
                                   
!<input>
  ! The source matrix.
  TYPE(t_matrixBlock), INTENT(IN) :: rsourceMatrix

  ! Whether to assemble a new structure in rdestMatrix.
  ! =TRUE (Standard): Release any old data/structure from rdestMatrix
  !   and create a new one.
  ! =FALSE          : rdestMatrix is assumed to be an existing matrix
  !   of the correct dimension/size which does not have to be rebuild.
  LOGICAL, INTENT(IN) :: bstructure
  
  ! Whether to assemble the matrix content.
  ! =TRUE (Standard): The content of rsourceMatrix is build in rdestMatrix.
  ! =FALSE          : No matrix content is build up in rdestMatrix.
  LOGICAL, INTENT(IN) :: bcontent
  
  ! OPTIONAL: Target format of the matrix rdestMatrix. Standard is Format 9.
  INTEGER, INTENT(IN), OPTIONAL :: cmatrixFormat
  
  ! OPTIONAL: Data type for the entries of rdestMatrix. 
  ! Standard is double precision.
  INTEGER, INTENT(IN), OPTIONAL :: cdataType

!</input>

!<output>
  ! The destination matrix. If this matrix contains valid handles to
  ! allocated memory blocks of the correct size for structure/entries,
  ! this data is overwritten. If not, arrays are (re-)allocated in
  ! the correct size.
  TYPE(t_matrixBlock), INTENT(INOUT) :: rdestMatrix
!</output>

!</subroutine>
    
    ! local variables
    INTEGER :: cdataTypeLocal, cmatrixFormatLocal,i,j
    LOGICAL :: balloc
    TYPE(t_matrixBlock) :: rlocalMatrix
    TYPE(t_matrixScalar) :: rlocalMatrixScalar
    INTEGER(PREC_VECIDX), DIMENSION(MAX(rsourceMatrix%ndiagBlocks,1)+1) :: Icolumns,Irows
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiagonal,p_Kld
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_MATIDX) :: isize
    
    ! Initialise values for data type and matrix format.
    cdataTypeLocal = ST_DOUBLE
    cmatrixFormatLocal = LSYSSC_MATRIX9
    
    IF (.NOT. bstructure) cmatrixFormatLocal = &
      rdestMatrix%RmatrixBlock(1,1)%cmatrixFormat
    IF (.NOT. bcontent) cdataTypeLocal = &
      rdestMatrix%RmatrixBlock(1,1)%cdataType
    
    IF (PRESENT(cdataType)) cdataTypeLocal = cdataType
    IF (PRESENT(cmatrixFormat)) cmatrixFormatLocal = cmatrixFormat
    
    ! Up to now, we don't support everything!
    ! Cancel if the input parameters want too much of us...
    IF (rsourceMatrix%ndiagBlocks .EQ. 0) RETURN
    
    IF (cdataTypeLocal .NE. ST_DOUBLE) THEN
      PRINT *,'glsys_assembleGlobal: Only double precision dest. matrix supported!'
      CALL sys_halt()
    END IF

    IF (cmatrixFormatLocal .NE. LSYSSC_MATRIX9) THEN
      PRINT *,'glsys_assembleGlobal: Only format 9 dest. matrix supported!'
      CALL sys_halt()
    END IF
    
    DO j=1,rsourceMatrix%ndiagBlocks
      DO i=1,rsourceMatrix%ndiagBlocks
      
        IF (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j)) THEN
        
          IF (rsourceMatrix%RmatrixBlock(i,j)%cdataType .NE. ST_DOUBLE) THEN
            PRINT *,'glsys_assembleGlobal: Only double precision source matrices &
                    &supported!'
            CALL sys_halt()
          END IF

          IF (rsourceMatrix%RmatrixBlock(i,j)%cmatrixFormat .NE. LSYSSC_MATRIX9) THEN
            PRINT *,'glsys_assembleGlobal: Only format 9 source matrices supported!'
            CALL sys_halt()
          END IF
          
        END IF
        
      END DO
    END DO

  ! Now start the actual assembly. What and how to assemble depend on
  ! the boolean input variables.
  IF ((.NOT. bstructure) .AND. (.NOT. bcontent)) RETURN ! nothing to do.
  
  ! Problem: Some of the submatrices may be virtually transposed,
  ! but the algorithms here can only handle un-transposed matrices.
  ! So at first, we create an un-transposed global matrix.
  
  CALL lsysbl_duplicateMatrix (rsourceMatrix,rlocalMatrix, &
      LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
  DO j=1,rlocalMatrix%ndiagBlocks
    DO i=1,rlocalMatrix%ndiagBlocks
        
      IF (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j)) THEN
        ! Transpose the submatrix if necessary
        IF (IAND(rsourceMatrix%RmatrixBlock(i,j)%imatrixSpec, &
                LSYSSC_MSPEC_TRANSPOSED) .NE. 0) THEN
          ! Untranspose the source-submatrix to a local matrix
          CALL lsyssc_transposeMatrix (rsourceMatrix%RmatrixBlock(i,j),&
                      rlocalMatrixScalar,LSYSSC_TR_VIRTUAL)
          
          ! Retranspose it - not vitually, but create a real transposed copy.
          IF (bcontent) THEN
            ! Transpose everything
            CALL lsyssc_releaseMatrix(rlocalMatrix%RmatrixBlock(i,j))
            CALL lsyssc_transposeMatrix (rlocalMatrixScalar,&
                        rlocalMatrix%RmatrixBlock(i,j),LSYSSC_TR_ALL)
          ELSE
            ! Transpose only the structure, ignore the content
            CALL lsyssc_releaseMatrix(rlocalMatrix%RmatrixBlock(i,j))
            CALL lsyssc_transposeMatrix (rlocalMatrixScalar,&
                        rlocalMatrix%RmatrixBlock(i,j),LSYSSC_TR_STRUCTURE)
          END IF
                
        END IF
      END IF
    END DO
  END DO
  
  ! Ok, rlocalMatrix is now a transposed-free source matrix.
  
  IF (bstructure .AND. bcontent) THEN
    
    ! Initialise general data of the destination matrix.
    CALL glmatasm_initDestination (rlocalMatrix,rdestMatrix,&
                                   cmatrixFormatLocal,cdataTypeLocal)
    
    ! Assemble the global matrix - structure and entries.
    !
    ! What's the destination matrix structure?
    SELECT CASE (cmatrixFormatLocal)
    CASE (LSYSSC_MATRIX9)
    
      ! Get the basic row/column indices of the destination matrix.
      CALL glmatasm_getOffsets (rlocalMatrix,Icolumns,Irows)

      ! Allocate a KLD in the destination matrix if we don't have a previous
      ! array in the correct size.
      balloc = .TRUE.
      IF ((.NOT. lsyssc_isMatrixStructureShared (rdestMatrix%RmatrixBlock(1,1))) .AND. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Kld .NE. ST_NOHANDLE)) THEN
        CALL storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Kld,isize)
        ! Matrix is for sure not transposed!
        IF (isize .EQ. rlocalMatrix%NEQ+1) THEN
          balloc = .FALSE.
        ELSE
          ! Release the previous memory before allocating a new block.
          CALL storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Kld)
        END IF
      END IF
      
      IF (balloc) CALL storage_new ('glsys_assembleGlobal', 'Kld',  &
                                    rlocalMatrix%NEQ+1_PREC_VECIDX, &
                                    ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kld,&
                                    ST_NEWBLOCK_ZERO)
      
      ! Set up KLD and NA of the destination matrix
      CALL glmatasm_KLD (rlocalMatrix,rdestMatrix%RmatrixBlock(1,1))

      ! Allocate a KCOL in the destination matrix if we don't have a previous
      ! array in the correct size.
      balloc = .TRUE.
      IF ((.NOT. lsyssc_isMatrixStructureShared (rdestMatrix%RmatrixBlock(1,1))) .AND. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Kcol .NE. ST_NOHANDLE)) THEN
        CALL storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Kcol,isize)
        ! Matrix is for sure not transposed!
        IF (isize .EQ. rdestMatrix%RmatrixBlock(1,1)%NA) THEN
          balloc = .FALSE.
        ELSE
          ! Release the previous memory before allocating a new block.
          CALL storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Kcol)
        END IF
      END IF
      IF (balloc) CALL storage_new ('glsys_assembleGlobal', 'Kcol', &
                                    rdestMatrix%RmatrixBlock(1,1)%NA, &
                                    ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kcol,&
                                    ST_NEWBLOCK_NOINIT)
                        
      ! Allocate the data array if we don't have a previous
      ! array in the correct size.
      balloc = .TRUE.
      IF ((.NOT. lsyssc_isMatrixContentShared (rdestMatrix%RmatrixBlock(1,1))) .AND. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Da .NE. ST_NOHANDLE)) THEN
        CALL storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Da,isize)
        IF (isize .EQ. rdestMatrix%RmatrixBlock(1,1)%NA) THEN
          balloc = .FALSE.
        ELSE
          ! Release the previous memory before allocating a new block.
          CALL storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Da)
        END IF
      END IF
      IF (balloc) CALL storage_new ('glsys_assembleGlobal', 'Da', & 
                                    rdestMatrix%RmatrixBlock(1,1)%NA, &
                                    cdataTypeLocal, &
                                    rdestMatrix%RmatrixBlock(1,1)%h_Da,&
                                    ST_NEWBLOCK_NOINIT)
            
      ! Set up KCOL and the matrix entries
      CALL glmatasm_KcolDa99dble (rlocalMatrix,rdestMatrix%RmatrixBlock(1,1),&
                                  Icolumns, Irows)      
      
      ! Allocate a Kdiagonal in the destination matrix if we don't have a previous
      ! array in the correct size.
      balloc = .TRUE.
      IF ((.NOT. lsyssc_isMatrixStructureShared (rdestMatrix%RmatrixBlock(1,1))) .AND. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal .NE. ST_NOHANDLE)) THEN
        CALL storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal,isize)
        ! Matrix is for sure not transposed!
        IF (isize .EQ. rlocalMatrix%NEQ) THEN
          balloc = .FALSE.
        ELSE
          ! Release the previous memory before allocating a new block.
          CALL storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal)
        END IF
      END IF

      ! Allocate a Kdiagonal in the destination matrix      
      IF (balloc) CALL storage_new (&
                        'glsys_assembleGlobal', 'Kdiagonal', rlocalMatrix%NEQ, &
                        ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal,&
                        ST_NEWBLOCK_NOINIT)

      ! Rebuild Kdiagonal
      CALL lsyssc_getbase_Kdiagonal(rdestMatrix%RmatrixBlock(1,1),p_Kdiagonal)
      CALL lsyssc_getbase_Kcol(rdestMatrix%RmatrixBlock(1,1),p_Kcol)
      CALL lsyssc_getbase_Kld(rdestMatrix%RmatrixBlock(1,1),p_Kld)
      CALL lsyssc_rebuildKdiagonal (p_Kcol, p_Kld, p_Kdiagonal, rdestMatrix%NEQ)
      
    END SELECT
    
  ELSE IF (bstructure) THEN
  
    ! Assemble only the structure of the global matrix
    
    ! Initialise general data of the destination matrix.
    CALL glmatasm_initDestination (rlocalMatrix,rdestMatrix,&
                                   cmatrixFormatLocal,cdataTypeLocal)
    
    ! What's the destination matrix structure?
    SELECT CASE (cmatrixFormatLocal)
    CASE (LSYSSC_MATRIX9)
    
      ! Get the basic row/column indices of the destination matrix.
      CALL glmatasm_getOffsets (rlocalMatrix,Icolumns,Irows)

      ! Allocate a KLD in the destination matrix if we don't have a previous
      ! array in the correct size.
      balloc = .TRUE.
      IF ((.NOT. lsyssc_isMatrixStructureShared (rdestMatrix%RmatrixBlock(1,1))) .AND. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Kld .NE. ST_NOHANDLE)) THEN
        CALL storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Kld,isize)
        ! Matrix is for sure not transposed!
        IF (isize .EQ. rlocalMatrix%NEQ+1) THEN
          balloc = .FALSE.
        ELSE
          ! Release the previous memory before allocating a new block.
          CALL storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Kld)
        END IF
      END IF
      
      IF (balloc) CALL storage_new ('glsys_assembleGlobal', 'Kld',  &
                                    rlocalMatrix%NEQ+1_PREC_VECIDX, &
                                    ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kld,&
                                    ST_NEWBLOCK_ZERO)
      
      ! Set up KLD and NA of the destination matrix
      CALL glmatasm_KLD (rlocalMatrix,rdestMatrix%RmatrixBlock(1,1))

      ! Allocate a KCOL in the destination matrix if we don't have a previous
      ! array in the correct size.
      balloc = .TRUE.
      IF ((.NOT. lsyssc_isMatrixStructureShared (rdestMatrix%RmatrixBlock(1,1))) .AND. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Kcol .NE. ST_NOHANDLE)) THEN
        CALL storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Kcol,isize)
        ! Matrix is for sure not transposed!
        IF (isize .EQ. rdestMatrix%RmatrixBlock(1,1)%NA) THEN
          balloc = .FALSE.
        ELSE
          ! Release the previous memory before allocating a new block.
          CALL storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Kcol)
        END IF
      END IF
      IF (balloc) CALL storage_new ('glsys_assembleGlobal', 'Kcol', &
                                    rdestMatrix%RmatrixBlock(1,1)%NA, &
                                    ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kcol,&
                                    ST_NEWBLOCK_NOINIT)
                        
      ! Set up KCOL and the matrix entries
      CALL glmatasm_Kcol99dble (rlocalMatrix,rdestMatrix%RmatrixBlock(1,1),&
                                Icolumns, Irows)      
      
      ! Allocate a Kdiagonal in the destination matrix      
      CALL storage_new ('glsys_assembleGlobal', 'Kdiagonal', rlocalMatrix%NEQ, &
                        ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal,&
                        ST_NEWBLOCK_NOINIT)

      ! Rebuild Kdiagonal
      CALL lsyssc_getbase_Kdiagonal(rdestMatrix%RmatrixBlock(1,1),p_Kdiagonal)
      CALL lsyssc_getbase_Kcol(rdestMatrix%RmatrixBlock(1,1),p_Kcol)
      CALL lsyssc_getbase_Kld(rdestMatrix%RmatrixBlock(1,1),p_Kld)
      CALL lsyssc_rebuildKdiagonal (p_Kcol, p_Kld, p_Kdiagonal, rdestMatrix%NEQ)
      
    END SELECT
    
  ELSE 
  
    ! Assemble only the entries of the global matrix; the structure
    ! is already there.
    ! Initialise general data of the destination matrix.
    CALL glmatasm_initDestination (rlocalMatrix,rdestMatrix,&
                                   cmatrixFormatLocal,cdataTypeLocal)
    
    ! What's the destination matrix structure?
    SELECT CASE (cmatrixFormatLocal)
    CASE (LSYSSC_MATRIX9)
    
      ! Get the basic row/column indices of the destination matrix.
      CALL glmatasm_getOffsets (rlocalMatrix,Icolumns,Irows)

      ! KCol/Kld are assumed to be ok.
      !
      ! Allocate the data array if we don't have a previous
      ! array in the correct size.
      balloc = .TRUE.
      IF ((.NOT. lsyssc_isMatrixContentShared (rdestMatrix%RmatrixBlock(1,1))) .AND. &
          (rdestMatrix%RmatrixBlock(1,1)%h_Da .NE. ST_NOHANDLE)) THEN
        CALL storage_getsize(rdestMatrix%RmatrixBlock(1,1)%h_Da,isize)
        IF (isize .EQ. rdestMatrix%RmatrixBlock(1,1)%NA) THEN
          balloc = .FALSE.
        ELSE
          ! Release the previous memory before allocating a new block.
          CALL storage_free (rdestMatrix%RmatrixBlock(1,1)%h_Da)
        END IF
      END IF
      IF (balloc) CALL storage_new ('glsys_assembleGlobal', 'Da', & 
                                    rdestMatrix%RmatrixBlock(1,1)%NA, &
                                    cdataTypeLocal, &
                                    rdestMatrix%RmatrixBlock(1,1)%h_Da,&
                                    ST_NEWBLOCK_NOINIT)
            
      ! Set up KCOL and the matrix entries
      CALL glmatasm_Da99dble (rlocalMatrix,rdestMatrix%RmatrixBlock(1,1),&
                              Icolumns, Irows)      
      
      ! Allocate a Kdiagonal in the destination matrix      
      CALL storage_new ('glsys_assembleGlobal', 'Kdiagonal', rlocalMatrix%NEQ, &
                        ST_INT, rdestMatrix%RmatrixBlock(1,1)%h_Kdiagonal,&
                        ST_NEWBLOCK_NOINIT)

      ! Rebuild Kdiagonal
      CALL lsyssc_getbase_Kdiagonal(rdestMatrix%RmatrixBlock(1,1),p_Kdiagonal)
      CALL lsyssc_getbase_Kcol(rdestMatrix%RmatrixBlock(1,1),p_Kcol)
      CALL lsyssc_getbase_Kld(rdestMatrix%RmatrixBlock(1,1),p_Kld)
      CALL lsyssc_rebuildKdiagonal (p_Kcol, p_Kld, p_Kdiagonal, rdestMatrix%NEQ)
      
    END SELECT
    
  END IF

  ! Release the local matrix.
  ! Note that only those submatrices are released from the heap that we
  ! created by transposing them above!
  CALL lsysbl_releaseMatrix(rlocalMatrix)
  
  CONTAINS
  
    !------------------------------------------------------------------
    ! Set up general data of the destination matrix
    
    SUBROUTINE glmatasm_initDestination (rsourceMatrix,rdestMatrix, &
                                         cmatrixFormat,cdataType)

    ! The source block matrix 
    TYPE(t_matrixBlock), INTENT(IN) :: rsourceMatrix
    
    ! The matrix to be initialised
    TYPE(t_matrixBlock), INTENT(INOUT) :: rdestMatrix
  
    ! Target format of the matrix rdestMatrix. Standard is Format 9.
    INTEGER, INTENT(IN) :: cmatrixFormat
    
    ! Data type for the entries of rdestMatrix. 
    ! Standard is double precision.
    INTEGER, INTENT(IN) :: cdataType

      IF (rdestMatrix%ndiagBlocks .LT. 1) THEN
        ! Create a new 1x1 matrix if necessary.
        CALL lsysbl_releaseMatrix (rdestMatrix)
        CALL lsysbl_createEmptyMatrix (rdestMatrix,1)
      END IF
      rdestMatrix%NEQ = rsourceMatrix%NEQ
      rdestMatrix%NCOLS = rsourceMatrix%NCOLS
      rdestMatrix%imatrixSpec = rsourceMatrix%imatrixSpec
      
      ! There's no appropriate discretisation or boundary condition
      ! structure for a global matrix! These only apply to
      ! scalar discretisations.
      NULLIFY(rdestMatrix%p_rblockDiscrTest)
      NULLIFY(rdestMatrix%p_rblockDiscrTrial)
      rdestMatrix%bidenticalTrialAndTest = .true.
      NULLIFY(rdestMatrix%p_rdiscreteBC)
      NULLIFY(rdestMatrix%p_rdiscreteBCfict)

      ! Initialise structural information of the submatrix      
      rdestMatrix%RmatrixBlock(1,1)%NEQ = rsourceMatrix%NEQ
      rdestMatrix%RmatrixBlock(1,1)%NCOLS = rsourceMatrix%NCOLS
      rdestMatrix%RmatrixBlock(1,1)%cmatrixFormat = cmatrixFormat
      rdestMatrix%RmatrixBlock(1,1)%cdataType = cdataType
      rdestMatrix%RmatrixBlock(1,1)%imatrixSpec = 0
      rdestMatrix%RmatrixBlock(1,1)%isortStrategy = 0
      rdestMatrix%RmatrixBlock(1,1)%h_IsortPermutation = ST_NOHANDLE
      rdestMatrix%RmatrixBlock(1,1)%dscaleFactor = 1.0_DP

    END SUBROUTINE
  
    !------------------------------------------------------------------
    ! Calculates the row-block offsets in the destination matrix
    SUBROUTINE glmatasm_getOffsets (rsourceMatrix,Icolumns,Irows)

    ! The source block matrix 
    TYPE(t_matrixBlock), INTENT(IN), TARGET :: rsourceMatrix
    
    ! Real column numbers in each block-column of the block matrix.
    ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Icolumns

    ! Real row numbers in each block-column of the block matrix.
    ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Irows

      ! local variables
      INTEGER :: i,j
      INTEGER(PREC_VECIDX) :: irow
      INTEGER(PREC_MATIDX) :: irowoffset
      
      Icolumns(:) = 0
      Irows(:) = 0
      
      ! Loop through all matrix blocks.
      DO i=1,rsourceMatrix%ndiagBlocks
      
        DO j=1,rsourceMatrix%ndiagBlocks
          ! When checking for the presence of the matrix, don't respect
          ! the scaling factor; we only want to get the size of the matrix 
          ! columns/rows!
          IF (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j,.TRUE.)) THEN
            
            Icolumns(j+1) = rsourceMatrix%RmatrixBlock(i,j)%NCOLS
            Irows(i+1) = rsourceMatrix%RmatrixBlock(i,j)%NEQ
          
          END IF
        END DO
        
      END DO
      
      ! Icolumns gives the number of real columns/rows. Add them together
      ! to get the real column/row numbers, each block column/row in the
      ! global matrix starts with.
      Icolumns(1) = 1
      Irows(1) = 1
      DO i=2,rsourceMatrix%ndiagBlocks
        Icolumns(i) = Icolumns(i) + Icolumns(i-1)
        Irows(i) = Irows(i) + Irows(i-1)
      END DO
      
    END SUBROUTINE

    !------------------------------------------------------------------
    ! Calculates NA and the KLD row structure of the global matrix
    SUBROUTINE glmatasm_KLD (rsourceMatrix,rdestMatrix)

    ! The source block matrix 
    TYPE(t_matrixBlock), INTENT(IN), TARGET :: rsourceMatrix
    
    ! The scalar submatrix which KLD is to be initialised.
    ! KLD must exist and be filled with 0.
    TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
    
      ! local variables
      INTEGER :: i,j
      INTEGER(PREC_VECIDX) :: irow
      INTEGER(PREC_MATIDX) :: irowoffset,narow
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld,p_KldDest
      TYPE(t_matrixScalar), POINTER :: p_rmatrix
      
      ! Create a new, empty KLD.
      CALL lsyssc_getbase_Kld (rdestMatrix,p_KldDest)
                        
      ! Loop through all matrix blocks.
      DO i=1,rsourceMatrix%ndiagBlocks
      
        ! Get the starting position of this block-row in the global matrix
        irowoffset = Irows(i)-1
        
        narow = 0
        
        DO j=1,rsourceMatrix%ndiagBlocks
          IF (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j)) THEN
            
            p_rmatrix => rsourceMatrix%RmatrixBlock(i,j)
            CALL lsyssc_getbase_Kld (p_rmatrix,p_Kld)
            
            ! Loop through all lines in the matrix and add the number of
            ! entries per row to get KLD:
            DO irow = 1,p_rmatrix%NEQ
              ! The first entriy in KLD is always one; add the length
              ! of the line to KLD one element shifted, so calculation
              ! ok KLD is easier later.
              p_KldDest(1+irow+irowoffset) = &
                p_KldDest(1+irow+irowoffset) + (p_Kld(irow+1)-p_Kld(irow))
            END DO
            
            ! Calculate the number of entries in this matrix-row-block
            narow = narow + rsourceMatrix%RmatrixBlock(i,j)%NA
          
          END IF
        END DO
        
      END DO
      
      ! Now we have:
      ! KLD(1) = 0, 
      ! KLD(2) = Number of entries in row 1,
      ! KLD(3) = Number of entries in row 2, etc.
      ! Sum up the values to get the actual KLD.
      p_KldDest(1) = 1
      DO irow = 1,rsourceMatrix%NEQ
        p_KldDest(irow+1) = p_KldDest(irow+1) + p_KldDest(irow)
      END DO
      
      ! and so we have NA. 
      rdestMatrix%NA = p_KldDest(rsourceMatrix%NEQ+1)-1
    
    END SUBROUTINE

    !------------------------------------------------------------------
    ! Calculates the KCOL column structure of the global matrix
    ! The KLD/NA structure must already be present in rdestMatrix!
    !
    ! Matrix-structure-9 source -> Matrix-structure-9 destination,
    SUBROUTINE glmatasm_Kcol99dble (rsourceMatrix,rdestMatrix,&
                                    Icolumns, Irows)

    ! The source block matrix 
    TYPE(t_matrixBlock), INTENT(IN), TARGET :: rsourceMatrix
    
    ! The scalar submatrix which KLD is to be initialised
    TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
    
    ! Real column numbers in each block-column of the block matrix.
    ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Icolumns

    ! Real row numbers in each block-column of the block matrix.
    ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Irows

      ! local variables
      INTEGER :: i,j,h_KldTmp,icol
      INTEGER(PREC_VECIDX) :: irow,ncols,irowGlobal
      INTEGER(PREC_MATIDX) :: ioffsetGlobal,ioffsetLocal,narow,isize
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld,p_KldDest,p_KldTmp
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol,p_KcolDest
      TYPE(t_matrixScalar), POINTER :: p_rmatrix
      
      ! Create a copy of KLD which we use for storing the index
      ! how many entries are already allocated in each row.
      h_KldTmp = ST_NOHANDLE
      CALL storage_copy (rdestMatrix%h_Kld,h_KldTmp)
      CALL storage_getbase_int (h_KldTmp,p_KldTmp)
                        
      ! Get KCol,Kld
      CALL lsyssc_getbase_Kcol (rdestMatrix,p_KcolDest)
      CALL lsyssc_getbase_Kld (rdestMatrix,p_KldDest)
      
      ! Loop through all matrix subblocks
      DO i=1,rsourceMatrix%ndiagBlocks
      
        ! Global row index?
        irowGlobal = Irows(i)-1
        
        DO j=1,rsourceMatrix%ndiagBlocks
        
          IF (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j)) THEN
    
            ! Get the submatrix
            p_rmatrix => rsourceMatrix%RmatrixBlock(i,j)
            
            ! Get the local matrix pointers / column structure
            CALL lsyssc_getbase_Kcol (p_rmatrix,p_Kcol)
            CALL lsyssc_getbase_Kld (p_rmatrix,p_Kld)
            
            ! Loop through all rows to append them to the current rows.
            DO irow = 1,p_rmatrix%NEQ
            
              ! How many elements to append?
              ncols = p_Kld(irow+1)-p_Kld(irow)
              
              ! Position of the data?
              ioffsetGlobal = p_KldTmp(irowGlobal+irow)
              ioffsetLocal  = p_Kld(irow)
              
              ! Copy matrix data and the column numbers
              p_KcolDest(ioffsetGlobal:ioffsetGlobal+ncols-1) = &
                p_Kcol(ioffsetLocal:ioffsetLocal+ncols-1)
              
              ! Increase the column numbers by the global column number
              ! of that matrix-block-column
              p_KcolDest(ioffsetGlobal:ioffsetGlobal+ncols-1) = &
                p_KcolDest(ioffsetGlobal:ioffsetGlobal+ncols-1) + Icolumns(j)-1
              
              ! Increase the counter/index position array for how 
              ! many elements are added to that row.
              p_KldTmp(irowGlobal+irow) = p_KldTmp(irowGlobal+irow) + ncols
            
            END DO ! irow
          
          END IF ! neq != 0
    
        END DO ! j
      END DO ! i
      
      ! Release the temp array
      CALL storage_free (h_KldTmp)
                        
    END SUBROUTINE

    !------------------------------------------------------------------
    ! Transfers the entries of the local matrices into the
    ! global matrix.
    ! The KLD/NA structure must already be present in rdestMatrix!
    !
    ! Matrix-structure-9 source -> Matrix-structure-9 destination,
    ! double precision vection
    SUBROUTINE glmatasm_Da99dble (rsourceMatrix,rdestMatrix,&
                                  Icolumns, Irows)

    ! The source block matrix 
    TYPE(t_matrixBlock), INTENT(IN), TARGET :: rsourceMatrix
    
    ! The scalar submatrix which KLD is to be initialised
    TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
    
    ! Real column numbers in each block-column of the block matrix.
    ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Icolumns

    ! Real row numbers in each block-column of the block matrix.
    ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Irows

      ! local variables
      INTEGER :: i,j,h_KldTmp,icol
      REAL(DP) :: dscale
      INTEGER(PREC_VECIDX) :: irow,ncols,irowGlobal
      INTEGER(PREC_MATIDX) :: ioffsetGlobal,ioffsetLocal,narow,isize
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld,p_KldDest,p_KldTmp
      REAL(DP), DIMENSION(:), POINTER :: p_Da,p_DaDest
      TYPE(t_matrixScalar), POINTER :: p_rmatrix
      
      ! Create a copy of KLD which we use for storing the index
      ! how many entries are already allocated in each row.
      h_KldTmp = ST_NOHANDLE
      CALL storage_copy (rdestMatrix%h_Kld,h_KldTmp)
      CALL storage_getbase_int (h_KldTmp,p_KldTmp)
                        
      ! Get Kld
      CALL lsyssc_getbase_Kld (rdestMatrix,p_KldDest)
      
      ! Get destination data arrays
      CALL lsyssc_getbase_double (rdestMatrix,p_DaDest)
    
      ! Loop through all matrix subblocks
      DO i=1,rsourceMatrix%ndiagBlocks
      
        ! Global row index?
        irowGlobal = Irows(i)-1
        
        DO j=1,rsourceMatrix%ndiagBlocks
        
          IF (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j)) THEN
    
            ! Get the submatrix
            p_rmatrix => rsourceMatrix%RmatrixBlock(i,j)
            
            ! Get the local matrix pointers structure
            CALL lsyssc_getbase_Kld (p_rmatrix,p_Kld)
            CALL lsyssc_getbase_double (p_rmatrix,p_Da)
            dscale = p_rmatrix%dscaleFactor
            
            ! Loop through all rows to append them to the current rows.
            DO irow = 1,p_rmatrix%NEQ
            
              ! How many elements to append?
              ncols = p_Kld(irow+1)-p_Kld(irow)
              
              ! Position of the data?
              ioffsetGlobal = p_KldTmp(irowGlobal+irow)
              ioffsetLocal  = p_Kld(irow)
              
              ! Copy matrix data 
              p_DaDest(ioffsetGlobal:ioffsetGlobal+ncols-1) = &
                dscale * p_Da(ioffsetLocal:ioffsetLocal+ncols-1)
              
              ! Increase the counter/index position array for how 
              ! many elements are added to that row.
              p_KldTmp(irowGlobal+irow) = p_KldTmp(irowGlobal+irow) + ncols
            
            END DO ! irow
          
          END IF ! neq != 0
    
        END DO ! j
      END DO ! i
      
      ! Release the temp array
      CALL storage_free (h_KldTmp)
                        
    END SUBROUTINE

    !------------------------------------------------------------------
    ! Calculates the KCOL column structure of the global matrix
    ! and transfers the entries of the local matrices into the
    ! global matrix.
    ! The KLD/NA structure must already be present in rdestMatrix!
    !
    ! Matrix-structure-9 source -> Matrix-structure-9 destination,
    ! double precision vection
    SUBROUTINE glmatasm_KcolDa99dble (rsourceMatrix,rdestMatrix,&
                                      Icolumns, Irows)

    ! The source block matrix 
    TYPE(t_matrixBlock), INTENT(IN), TARGET :: rsourceMatrix
    
    ! The scalar submatrix which KLD is to be initialised
    TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
    
    ! Real column numbers in each block-column of the block matrix.
    ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Icolumns

    ! Real row numbers in each block-column of the block matrix.
    ! DIMENSION(rsourceMatrix%ndiagBlocks+1)
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Irows

      ! local variables
      INTEGER :: i,j,h_KldTmp,icol
      REAL(DP) :: dscale
      INTEGER(PREC_VECIDX) :: irow,ncols,irowGlobal
      INTEGER(PREC_MATIDX) :: ioffsetGlobal,ioffsetLocal,narow,isize
      INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld,p_KldDest,p_KldTmp
      INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol,p_KcolDest
      REAL(DP), DIMENSION(:), POINTER :: p_Da,p_DaDest
      TYPE(t_matrixScalar), POINTER :: p_rmatrix
      
      ! Create a copy of KLD which we use for storing the index
      ! how many entries are already allocated in each row.
      h_KldTmp = ST_NOHANDLE
      CALL storage_copy (rdestMatrix%h_Kld,h_KldTmp)
      CALL storage_getbase_int (h_KldTmp,p_KldTmp)
                        
      ! Get KCol,Kld
      CALL lsyssc_getbase_Kcol (rdestMatrix,p_KcolDest)
      CALL storage_getbase_int (rdestMatrix%h_Kld,p_KldDest)
      
      ! Get destination data arrays
      CALL lsyssc_getbase_double (rdestMatrix,p_DaDest)
    
      ! Loop through all matrix subblocks
      DO i=1,rsourceMatrix%ndiagBlocks
      
        ! Global row index?
        irowGlobal = Irows(i)-1
        
        DO j=1,rsourceMatrix%ndiagBlocks
        
          IF (lsysbl_isSubmatrixPresent (rsourceMatrix,i,j)) THEN
    
            ! Get the submatrix
            p_rmatrix => rsourceMatrix%RmatrixBlock(i,j)
            
            ! Get the local matrix pointers / column structure
            CALL lsyssc_getbase_Kcol (p_rmatrix,p_Kcol)
            CALL lsyssc_getbase_Kld (p_rmatrix,p_Kld)
            CALL lsyssc_getbase_double (p_rmatrix,p_Da)
            dscale = p_rmatrix%dscaleFactor
            
            ! Loop through all rows to append them to the current rows.
            DO irow = 1,p_rmatrix%NEQ
            
              ! How many elements to append?
              ncols = p_Kld(irow+1)-p_Kld(irow)
              
              ! Position of the data?
              ioffsetGlobal = p_KldTmp(irowGlobal+irow)
              ioffsetLocal  = p_Kld(irow)
              
              ! Copy matrix data and the column numbers
              p_DaDest(ioffsetGlobal:ioffsetGlobal+ncols-1) = &
                dscale * p_Da(ioffsetLocal:ioffsetLocal+ncols-1)
              p_KcolDest(ioffsetGlobal:ioffsetGlobal+ncols-1) = &
                p_Kcol(ioffsetLocal:ioffsetLocal+ncols-1)
              
              ! Increase the column numbers by the global column number
              ! of that matrix-block-column
              p_KcolDest(ioffsetGlobal:ioffsetGlobal+ncols-1) = &
                p_KcolDest(ioffsetGlobal:ioffsetGlobal+ncols-1) + Icolumns(j)-1
              
              ! Increase the counter/index position array for how 
              ! many elements are added to that row.
              p_KldTmp(irowGlobal+irow) = p_KldTmp(irowGlobal+irow) + ncols
            
            END DO ! irow
          
          END IF ! neq != 0
    
        END DO ! j
      END DO ! i
      
      ! Release the temp array
      CALL storage_free (h_KldTmp)
                        
    END SUBROUTINE

  END SUBROUTINE


END MODULE
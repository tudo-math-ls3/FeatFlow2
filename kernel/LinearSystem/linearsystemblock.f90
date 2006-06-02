!##############################################################################
!# ****************************************************************************
!# <name> linearsystemblock </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic definitions and routines to maintain block
!# matrices and block vectors. A block matrix is realised as 2D-array of
!# scalar matrices, while a block vector is realised as 1D-array of
!# scalar vectors.
!#
!# The following routines can be found here:
!#
!#  1.) lsysbl_createVecBlockDirect
!#      -> Create a block vector by specifying the size of the subblocks
!# 
!#  2.) lsysbl_createVecBlockIndirect
!#      -> Create a block vector by copying the structure of another block
!#         vector
!#
!#  3.) lsysbl_createVecBlockIndMat
!#      -> Create a block vector according to a block matrix.
!#
!#  4.) lsysbl_createMatFromScalar
!#      -> Creates a 1x1 block matrix from a scalar matrix
!#
!#  5.) lsysbl_createVecFromScalar
!#      -> Creates a 1-block vector from a scalar vector
!#
!#  6.) lsysbl_enforceStructure
!#      -> Enforces the structure of a given block vector in another
!#         block vector
!#
!#  7.) lsysbl_assignDiscretIndirect
!#      -> Assign discretisation related information of one vector
!#         to another
!#
!#  8.) lsysbl_assignDiscretIndirectMat
!#      -> Assign discretisation related information of a matrix
!#         to a vector tp make it compatible.
!#
!#  9.) lsysbl_updateMatStrucInfo
!#      -> Recalculate structural data of a block matrix from
!#         the submatrices
!#
!# 10.) lsysbl_releaseVector
!#      -> Release a block vector from memory
!#
!# 11.) lsysbl_blockMatVec
!#      -> Multiply a block matrix with a block vector
!#
!# 12.) lsysbl_vectorCopy
!#       -> Copy a block vector over to another one
!#
!# 13.) lsysbl_vectorScale
!#      -> Scale a block vector by a constant
!#
!# 14.) lsysbl_vectorClear
!#      -> Clear a block vector
!#
!# 15.) lsysbl_vectorLinearComb
!#      -> Linear combination of two block vectors
!#
!# 16.) lsysbl_scalarProduct
!#      -> Calculate a scalar product of two vectors
!#
!# 17.) lsysbl_setSortStrategy
!#      -> Assigns a sorting strategy/permutation to every subvector
!#
!# 18.) lsysbl_sortVectorInSitu
!#      -> Resort the entries of all subvectors according to an assigned
!#         sorting strategy
!#
!# 19.) lsysbl_isVectorCompatible
!#      -> Checks whether two vectors are compatible to each other
!#
!# 20.) lsysbl_isMatrixCompatible
!#      -> Checks whether a matrix and a vector are compatible to each other
!#
!# 21.) lsysbl_isMatrixSorted
!#      -> Checks if a block matrix is sorted
!#
!# 22.) lsysbl_isVectorSorted
!#      -> Checks if a block vector is sorted
!# </purpose>
!##############################################################################

MODULE linearsystemblock

  USE fsystem
  USE storage
  USE linearsystemscalar
  USE linearalgebra
  
  IMPLICIT NONE

!<constants>

!<constantblock description="Global constants for block matrices">

  ! Maximum number of blocks that is allowed in 'usual' block matrices,
  ! supported by the routines in this file.
  INTEGER, PARAMETER :: LSYSBL_MAXBLOCKS = 16
  
!</constantblock>

!<constantblock description="Flags for the matrix specification bitfield">

  ! Standard matrix
  INTEGER, PARAMETER :: LSYSBS_MSPEC_GENERAL           =        0
  
  ! Block matrix is a scalar (i.e. 1x1) matrix.
  INTEGER, PARAMETER :: LSYSBS_MSPEC_SCALAR            =     2**0

  ! Block matrix is of saddle-point type:
  !  (A  B1)
  !  (B2 0 )
  INTEGER, PARAMETER :: LSYSBS_MSPEC_SADDLEPOINT       =     2**1

  ! Block matrix is nearly of saddle-point type 
  !  (A  B1)
  !  (B2 C )
  ! with C~0 being a stabilisation matrix
  INTEGER, PARAMETER :: LSYSBS_MSPEC_NEARLYSADDLEPOINT =     2**2

!</constantblock>

!</constants>

!<types>
  
!<typeblock>
  
  ! A block vector that can be used by block linear algebra routines.
  
  TYPE t_vectorBlock
    
    ! Total number of equations in the vector
    INTEGER                    :: NEQ = 0
    
    ! Handle identifying the vector entries or = ST_NOHANDLE if not
    ! allocated.
    INTEGER                    :: h_Ddata = ST_NOHANDLE
    
    ! Data type of the entries in the vector. Either ST_SINGLE or
    ! ST_DOUBLE. The subvectors are all of the same type.
    INTEGER                    :: cdataType = ST_DOUBLE

    ! Is set to true, if the handle h_Ddata belongs to another vector,
    ! i.e. when this vector shares data with another vector.
    LOGICAL                    :: bisCopy   = .FALSE.

    ! Number of blocks allocated in RvectorBlock
    INTEGER                    :: nblocks = 0
    
    ! A 1D array with scalar vectors for all the blocks.
    ! The handle identifier inside of these blocks are set to h_Ddata.
    ! The iidxFirstEntry index pointer in each subblock is set according
    ! to the position inside of the large vector.
    TYPE(t_vectorScalar), DIMENSION(LSYSBL_MAXBLOCKS) :: RvectorBlock
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! A block matrix that can be used by scalar linear algebra routines.
  
  TYPE t_matrixBlock
    
    ! Total number of equations in the vector
    INTEGER                    :: NEQ         = 0
    
    ! Number of diagonal blocks in the matrix
    INTEGER                    :: ndiagBlocks = 0

    ! Matrix specification tag. This is a bitfield coming from an OR 
    ! combination of different LSYSBS_MSPEC_xxxx constants and specifies
    ! various details of the matrix. If it is =LSYSBS_MSPEC_GENERAL,
    ! the matrix is a usual matrix that needs no special handling.
    INTEGER(I32) :: imatrixSpec = LSYSBS_MSPEC_GENERAL

    ! A 2D array with scalar matrices for all the blocks
    TYPE(t_matrixScalar), &
      DIMENSION(LSYSBL_MAXBLOCKS,LSYSBL_MAXBLOCKS) :: RmatrixBlock
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! This structure encapsules a pointer to a block matrix.
  
  TYPE t_matrixBlockPointer
    
    ! A pointer to a block matrix
    TYPE(t_matrixBlock), POINTER        :: ptr
    
  END TYPE
  
!</typeblock>

!</types>

  INTERFACE lsysbl_createVectorBlock
    MODULE PROCEDURE lsysbl_createVecBlockDirect 
    MODULE PROCEDURE lsysbl_createVecBlockIndirect 
  END INTERFACE

CONTAINS

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_isVectorCompatible (rvector1,rvector2,bcompatible)
  
!<description>
  ! Checks whether two vectors are compatible to each other.
  ! Two vectors are compatible if
  ! - they have the same size,
  ! - they have the same structure of subvectors,
  ! - they have the same sorting strategy OR they are both unsorted
!</description>

!<input>
  ! The first vector
  TYPE(t_vectorBlock), INTENT(IN) :: rvector1
  
  ! The second vector
  TYPE(t_vectorBlock), INTENT(IN) :: rvector2
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether the vectors are compatible or not.
  ! If not given, an error will inform the user if the two vectors are
  ! not compatible and the program will halt.
  LOGICAL, INTENT(OUT), OPTIONAL :: bcompatible
!</output>

!</subroutine>

  ! local variables
  INTEGER :: i

  ! We assume that we are not compatible
  IF (PRESENT(bcompatible)) bcompatible = .FALSE.
  
  ! Vectors must have the same size and number of blocks
  IF (rvector1%NEQ .NE. rvector2%NEQ) THEN
    IF (PRESENT(bcompatible)) THEN
      bcompatible = .FALSE.
      RETURN
    ELSE
      PRINT *,'Vectors not compatible, different size!'
      STOP
    END IF
  END IF

  IF (rvector1%nblocks .NE. rvector2%nblocks) THEN
    IF (PRESENT(bcompatible)) THEN
      bcompatible = .FALSE.
      RETURN
    ELSE
      PRINT *,'Vectors not compatible, different block structure!'
      STOP
    END IF
  END IF
  
  ! All the subblocks must be the same size and must be sorted the same way.
  DO i=1,rvector1%nblocks
    IF (rvector1%RvectorBlock(i)%NEQ .NE. rvector2%RvectorBlock(i)%NEQ) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Vectors not compatible, different block structure!'
        STOP
      END IF
    END IF
    
    ! isortStrategy < 0 means unsorted. Both unsorted is ok.
    
    IF ((rvector1%RvectorBlock(i)%isortStrategy .GT. 0) .OR. &
        (rvector2%RvectorBlock(i)%isortStrategy .GT. 0)) THEN

      IF (rvector1%RvectorBlock(i)%isortStrategy .NE. &
          rvector2%RvectorBlock(i)%isortStrategy) THEN
        IF (PRESENT(bcompatible)) THEN
          bcompatible = .FALSE.
          RETURN
        ELSE
          PRINT *,'Vectors not compatible, differently sorted!'
          STOP
        END IF
      END IF

      IF (rvector1%RvectorBlock(i)%h_isortPermutation .NE. &
          rvector2%RvectorBlock(i)%h_isortPermutation) THEN
        IF (PRESENT(bcompatible)) THEN
          bcompatible = .FALSE.
          RETURN
        ELSE
          PRINT *,'Vectors not compatible, differently sorted!'
          STOP
        END IF
      END IF
    END IF
  END DO

  ! Ok, they are compatible
  IF (PRESENT(bcompatible)) bcompatible = .TRUE.

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_isMatrixCompatible (rvector,rmatrix,bcompatible)
  
!<description>
  ! Checks whether a vector and a matrix are compatible to each other.
  ! A vector and a matrix are compatible if
  ! - they have the same NEQ,
  ! - they have the same structure of subvectors and diagonal blocks,
  ! - they have the same sorting strategy OR they are both unsorted
  ! - they have the same spatial discretisation,
  ! - they have the same boundary conditions.
!</description>

!<input>
  ! The vector
  TYPE(t_vectorBlock), INTENT(IN) :: rvector
  
  ! The matrix
  TYPE(t_matrixBlock), INTENT(IN) :: rmatrix
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether vector and matrix are compatible or not.
  ! If not given, an error will inform the user if the vector/matrix are
  ! not compatible and the program will halt.
  LOGICAL, INTENT(OUT), OPTIONAL :: bcompatible
!</output>

!</subroutine>

  ! local variables
  INTEGER :: i
  LOGICAL :: b1,b2,b3

  ! We assume that we are not compatible
  IF (PRESENT(bcompatible)) bcompatible = .FALSE.
  
  ! Vector/Matrix must have the same size and number of blocks
  IF (rvector%NEQ .NE. rmatrix%NEQ) THEN
    IF (PRESENT(bcompatible)) THEN
      bcompatible = .FALSE.
      RETURN
    ELSE
      PRINT *,'Vector/Matrix not compatible, different size!'
      STOP
    END IF
  END IF

  IF (rvector%nblocks .NE. rmatrix%ndiagBlocks) THEN
    IF (PRESENT(bcompatible)) THEN
      bcompatible = .FALSE.
      RETURN
    ELSE
      PRINT *,'Vector/Matrix not compatible, different block structure!'
      STOP
    END IF
  END IF
  
  ! All the subblocks must be the same size and must be sorted the same way.
  DO i=1,rvector%nblocks

    b1 = ASSOCIATED(rvector%RvectorBlock(i)%p_rdiscreteBC)
    b2 = ASSOCIATED(rmatrix%RmatrixBlock(i,i)%p_rdiscreteBC)
    b3 = ASSOCIATED(rvector%RvectorBlock(i)%p_rdiscreteBC, &
                    rmatrix%RmatrixBlock(i,i)%p_rdiscreteBC)
    IF ((b1 .OR. b2) .AND. .NOT. (b1 .AND. b2 .AND. b3)) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Vector/Matrix not compatible, different boundary conditions!'
        STOP
      END IF
    END IF

    b1 = ASSOCIATED(rvector%RvectorBlock(i)%p_rdiscreteBCfict)
    b2 = ASSOCIATED(rmatrix%RmatrixBlock(i,i)%p_rdiscreteBCfict)
    b3 = ASSOCIATED(rvector%RvectorBlock(i)%p_rdiscreteBCfict, &
                    rmatrix%RmatrixBlock(i,i)%p_rdiscreteBCfict)
    IF ((b1 .OR. b2) .AND. .NOT. (b1 .AND. b2 .AND. b3)) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Vector/Matrix not compatible, different fict. boundary conditions!'
        STOP
      END IF
    END IF

    b1 = ASSOCIATED(rvector%RvectorBlock(i)%p_rspatialDiscretisation)
    b2 = ASSOCIATED(rmatrix%RmatrixBlock(i,i)%p_rspatialDiscretisation)
    b3 = ASSOCIATED(rvector%RvectorBlock(i)%p_rspatialDiscretisation, &
                    rmatrix%RmatrixBlock(i,i)%p_rspatialDiscretisation)
    IF ((b1 .OR. b2) .AND. .NOT. (b1 .AND. b2 .AND. b3)) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Vector/Matrix not compatible, different discretisation!'
        STOP
      END IF
    END IF

    IF (rvector%RvectorBlock(i)%NEQ .NE. rmatrix%RmatrixBlock(i,i)%NEQ) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Vector/Matrix not compatible, different block structure!'
        STOP
      END IF
    END IF

    ! isortStrategy < 0 means unsorted. Both unsorted is ok.

    IF ((rvector%RvectorBlock(i)%isortStrategy .GT. 0) .OR. &
        (rmatrix%RmatrixBlock(i,i)%isortStrategy .GT. 0)) THEN

      IF (rvector%RvectorBlock(i)%isortStrategy .NE. &
          rmatrix%RmatrixBlock(i,i)%isortStrategy) THEN
        IF (PRESENT(bcompatible)) THEN
          bcompatible = .FALSE.
          RETURN
        ELSE
          PRINT *,'Vector/Matrix not compatible, differently sorted!'
          STOP
        END IF
      END IF

      IF (rvector%RvectorBlock(i)%h_isortPermutation .NE. &
          rmatrix%RmatrixBlock(i,i)%h_isortPermutation) THEN
        IF (PRESENT(bcompatible)) THEN
          bcompatible = .FALSE.
          RETURN
        ELSE
          PRINT *,'Vector/Matrix not compatible, differently sorted!'
          STOP
        END IF
      END IF
    END IF
  END DO

  ! Ok, they are compatible
  IF (PRESENT(bcompatible)) bcompatible = .TRUE.

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_createVecBlockDirect (rx, Isize, bclear, cdataType)
  
!<description>
  ! Initialises the vector block structure rx. Isize is an array
  ! of integers containing the length of the individual blocks.
  ! Memory is allocated on the heap to hold vectors of the size
  ! according to Isize.
!</description>

!<input>
  ! An array with length-tags for the different blocks
  INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Isize
  
  ! Optional: If set to YES, the vector will be filled with zero initially.
  LOGICAL, INTENT(IN), OPTIONAL             :: bclear
  
  ! OPTIONAL: Data type identifier for the entries in the vector. 
  ! Either ST_SINGLE or ST_DOUBLE. If not present, ST_DOUBLE is assumed.
  INTEGER, INTENT(IN),OPTIONAL              :: cdataType
!</input>

!<output>
  
  ! Destination structure. Memory is allocated for each of the blocks.
  TYPE(t_vectorBlock),INTENT(OUT) :: rx
  
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i,n
  INTEGER :: cdata
  
  cdata = ST_DOUBLE
  IF (PRESENT(cdataType)) cdata = cdataType
  
  ! Allocate one large vector holding all data.
  CALL storage_new1D ('lsysbl_createVecBlockDirect', 'Vector', SUM(Isize), cdata, &
                      rx%h_Ddata, ST_NEWBLOCK_NOINIT)
  rx%cdataType = cdata
  
  ! Initialise the sub-blocks. Save a pointer to the starting address of
  ! each sub-block.
  ! Denote in the subvector that the handle belongs to us - not to
  ! the subvector.
  
  n=1
  DO i = 1,SIZE(Isize)
    IF (Isize(i) .GT. 0) THEN
      rx%RvectorBlock(i)%NEQ = Isize(i)
      rx%RvectorBlock(i)%iidxFirstEntry = n
      rx%RvectorBlock(i)%h_Ddata = rx%h_Ddata
      rx%RvectorBlock(i)%cdataType = rx%cdataType
      rx%RvectorBlock(i)%bisCopy = .TRUE.
      rx%RvectorBlock(i)%cdataType = cdata
      n = n+Isize(i)
    ELSE
      rx%RvectorBlock(i)%NEQ = 0
      rx%RvectorBlock(i)%iidxFirstEntry = 0
      rx%RvectorBlock(i)%h_Ddata = ST_NOHANDLE
    END IF
  END DO
  
  rx%NEQ = n
  rx%nblocks = SIZE(Isize)
  
  ! The data of the vector belongs to us (we created the handle), 
  ! not to somebody else.
  rx%bisCopy = .FALSE.
  
  ! Warning: don't reformulate the following check into one IF command
  ! as this might give problems with some compilers!
  IF (PRESENT(bclear)) THEN
    IF (bclear) THEN
      CALL lsysbl_vectorClear (rx)
    END IF
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_createVecBlockIndirect (rtemplate,rx,bclear,cdataType)
  
!<description>
  ! Initialises the vector block structure rx. rtemplate is an
  ! existing block vector structure.
  ! Memory is allocated on the heap for rx, the different components
  ! of rx will have the same size as the corresponding ones in rtemplate.
!</description>
  
!<input>
  ! A template vector structure
  TYPE(t_vectorBlock),INTENT(IN) :: rtemplate
  
  ! Optional: If set to YES, the vector will be filled with zero initially.
  ! Otherwise the content of rx is undefined.
  LOGICAL, INTENT(IN), OPTIONAL             :: bclear
  
  ! OPTIONAL: Data type identifier for the entries in the vector. 
  ! Either ST_SINGLE or ST_DOUBLE. If not present, the new vector will
  ! receive the same data type as rtemplate.
  INTEGER, INTENT(IN),OPTIONAL              :: cdataType
!</input>

!<output>
  
  ! Destination structure. Memory is allocated for each of the blocks.
  TYPE(t_vectorBlock),INTENT(OUT) :: rx
  
!</output>
  
!</subroutine>

  INTEGER :: cdata
  
  cdata = rtemplate%cdataType
  IF (PRESENT(cdataType)) cdata = cdataType

  ! Copy the vector structure with all pointers.
  ! The handle identifying the vector data on the heap will be overwritten later.
  rx = rtemplate
  
  ! Allocate one large new vector holding all data.
  CALL storage_new1D ('lsysbl_createVecBlockDirect', 'Vector', rtemplate%NEQ, cdata, &
                      rx%h_Ddata, ST_NEWBLOCK_NOINIT)
  rx%cdataType = cdata
  
  WHERE (rtemplate%RvectorBlock%h_Ddata .NE. ST_NOHANDLE)
    ! Initialise the data type
    rx%RvectorBlock%cdataType = cdata
                     
    ! Put the new handle to all subvectors
    rx%RvectorBlock(:)%h_Ddata = rx%h_Ddata
  
    ! Note in the subvectors, that the handle belongs to somewhere else... to us!
    rx%RvectorBlock%bisCopy = .TRUE.
  END WHERE
  
  ! Our handle belongs to us, so it's not a copy of another vector.
  rx%bisCopy = .FALSE.

  ! Warning: don't reformulate the following check into one IF command
  ! as this might give problems with some compilers!
  IF (PRESENT(bclear)) THEN
    IF (bclear) THEN
      CALL lsysbl_vectorClear (rx)
    END IF
  END IF
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_createVecBlockIndMat (rtemplateMat,rx, bclear, cdataType)
  
!<description>
  ! Initialises the vector block structure rx. rtemplateMat is an
  ! existing block matrix structure. The vector rx will be created
  ! according to the size of the blocks of the submatrices.
  !
  ! Memory is allocated on the heap for rx. The different components
  ! of rx will have the same size as the corresponding diagonal
  ! blocks in rtemplateMat.
  ! The sorting strategies of the subvectors are initialised
  ! with the sorting strategies of the diagonal blocks of rtemplateMat.
!</description>
  
!<input>
  ! A template vector structure
  TYPE(t_matrixBlock), INTENT(IN) :: rtemplateMat
  
  ! OPTIONAL: If set to TRUE, the vector will be filled with zero initially.
  ! Otherwise the content of rx is undefined.
  LOGICAL, INTENT(IN), OPTIONAL   :: bclear
  
  ! OPTIONAL: Data type identifier for the entries in the vector. 
  ! Either ST_SINGLE or ST_DOUBLE. If not present, ST_DOUBLE is assumed.
  INTEGER, INTENT(IN),OPTIONAL              :: cdataType
!</input>

!<output>
  ! Destination structure. Memory is allocated for each of the blocks.
  TYPE(t_vectorBlock),INTENT(OUT) :: rx
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i,n
  INTEGER :: cdata
  
  cdata = ST_DOUBLE
  IF (PRESENT(cdataType)) cdata = cdataType
  
  ! Allocate one large vector holding all data.
  CALL storage_new1D ('lsysbl_createVecBlockDirect', 'Vector', rtemplateMat%NEQ, &
                      cdata, rx%h_Ddata, ST_NEWBLOCK_NOINIT)
  
  ! Initialise the sub-blocks. Save a pointer to the starting address of
  ! each sub-block.
  
  n=1
  DO i = 1,rtemplateMat%ndiagBlocks
    IF (rtemplateMat%RmatrixBlock(i,i)%NEQ .GT. 0) THEN
      rx%RvectorBlock(i)%NEQ = rtemplateMat%RmatrixBlock(i,i)%NEQ
      
      ! Take the handle of the complete-solution vector, but set the index of
      ! the first entry to a value >= 1 - so that it points to the first
      ! entry in the global solution vector!
      rx%RvectorBlock(i)%h_Ddata = rx%h_Ddata
      rx%RvectorBlock(i)%iidxFirstEntry = n
      
      ! Give the vector the discretisation of the matrix
      rx%RvectorBlock(i)%p_rspatialDiscretisation => &
        rtemplateMat%RmatrixBlock(i,i)%p_rspatialDiscretisation
        
      ! Give the vector the discrete boundary conditions
      rx%RvectorBlock(i)%p_rdiscreteBC => rtemplateMat%RmatrixBlock(i,i)%p_rdiscreteBC
      rx%RvectorBlock(i)%p_rdiscreteBCfict => &
        rtemplateMat%RmatrixBlock(i,i)%p_rdiscreteBCfict
        
      ! Give the vector the same sorting strategy as the matrix, so that
      ! the matrix and vector get compatible. Otherwise, things
      ! like matrix vector multiplication won't work...
      rx%RvectorBlock(i)%isortStrategy = &
        rtemplateMat%RmatrixBlock(i,i)%isortStrategy
      rx%RvectorBlock(i)%h_IsortPermutation = &
        rtemplateMat%RmatrixBlock(i,i)%h_IsortPermutation
      
      ! Denote in the subvector that the handle belongs to us - not to
      ! the subvector.
      rx%RvectorBlock(i)%bisCopy = .TRUE.
      
      ! Set the data type
      rx%RvectorBlock(i)%cdataType = cdata
      
      n = n+rtemplateMat%RmatrixBlock(i,i)%NEQ
    ELSE
      ! Let's hope this situation (an empty equation) never occurs - 
      ! might produce some errors elsewhere :)
      rx%RvectorBlock(i)%NEQ = 0
      rx%RvectorBlock(i)%iidxFirstEntry = 0
    END IF
  END DO
  
  rx%NEQ = rtemplateMat%NEQ
  rx%nblocks = rtemplateMat%ndiagBlocks
  rx%cdataType = cdata
  
  ! Warning: don't reformulate the following check into one IF command
  ! as this might give problems with some compilers!
  IF (PRESENT(bclear)) THEN
    IF (bclear) THEN
      CALL lsysbl_vectorClear (rx)
    END IF
  END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_assignDiscretIndirect (rtemplate,rx)
  
!<description>
  ! Assigns discretisation-related information (spatial discretisation, 
  ! boundary conditions,...) to rx. The vector rtemplate is used as a
  ! template, so at the end of the routine, rx and rtemplate share the
  ! same discretisation.
!</description>
  
!<input>
  ! A template vector structure
  TYPE(t_vectorBlock),INTENT(IN) :: rtemplate
!</input>

!<inputoutput>
  ! Destination structure. Discretisation-related information of rtemplate
  ! is assigned to rx.
  TYPE(t_vectorBlock),INTENT(INOUT) :: rx
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i

  ! Simply modify all pointers of all subvectors, that's it.
  DO i=1,rx%nblocks
    ! Spatial discretisation
    rx%RvectorBlock(i)%p_rspatialDiscretisation => &
      rtemplate%RvectorBlock(i)%p_rspatialDiscretisation 

    ! Discrete boundary conditions
    rx%RvectorBlock(i)%p_rdiscreteBC => &
      rtemplate%RvectorBlock(i)%p_rdiscreteBC

    ! Discrete fictitious boundary conditions
    rx%RvectorBlock(i)%p_rdiscreteBCfict => &
      rtemplate%RvectorBlock(i)%p_rdiscreteBCfict
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_assignDiscretIndirectMat (rtemplateMat,rx)
  
!<description>
  ! Assigns discretisation-related information (spatial discretisation, 
  ! boundary conditions,...) of a matrix to rx. The matrix rtemplateMat is 
  ! used as a template, so at the end of the routine, rx and rtemplate 
  ! share the same discretisation and boundary conditions.
  ! (More precisely, the blocks in rx and the diagonal blocks of
  !  rtemplateMat!)
  !
  ! In case the vector is completely incompatiblew to the matrix
  ! (different number of blocks, different NEQ), an error is thrown.
!</description>
  
!<input>
  ! A template vector structure
  TYPE(t_matrixBlock),INTENT(IN) :: rtemplateMat
!</input>

!<inputoutput>
  ! Destination structure. Discretisation-related information of rtemplate
  ! is assigned to rx.
  TYPE(t_vectorBlock),INTENT(INOUT) :: rx
!</output>
  
!</subroutine>

  ! local variables
  INTEGER :: i

  ! There must be at least some compatibility between the matrix
  ! and the vector!
  IF ((rtemplateMat%NEQ .NE. rx%NEQ) .OR. &
      (rtemplateMat%ndiagBlocks .NE. rx%nblocks)) THEN
    PRINT *,'lsysbl_assignDiscretIndirectMat error: Matrix/Vector incompatible!'
    STOP
  END IF

  ! Simply modify all pointers of all subvectors, that's it.
  DO i=1,rx%nblocks
    ! Spatial discretisation
    rx%RvectorBlock(i)%p_rspatialDiscretisation => &
      rtemplateMat%RmatrixBlock(i,i)%p_rspatialDiscretisation 

    ! Discrete boundary conditions
    rx%RvectorBlock(i)%p_rdiscreteBC => &
      rtemplateMat%RmatrixBlock(i,i)%p_rdiscreteBC

    ! Discrete fictitious boundary conditions
    rx%RvectorBlock(i)%p_rdiscreteBCfict => &
      rtemplateMat%RmatrixBlock(i,i)%p_rdiscreteBCfict
  END DO
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_releaseVector (rx)
  
!<description>
  ! Releases the memory that is reserved by rx.
  ! Remark: The memory associated to the vector is only released if this
  !  vector structure is the owner of the data.
!</description>
  
!<inputoutput>
  
  ! Vector structure that is to be released.
  TYPE(t_vectorBlock),INTENT(INOUT) :: rx
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i
  
  IF (rx%h_Ddata .EQ. ST_NOHANDLE) THEN
    PRINT *,'lsysbl_releaseVector warning: releasing unused vector.'
  END IF
  
  ! Release the data - if the handle is not a copy of another vector!
  IF ((.NOT. rx%bisCopy) .AND. (rx%h_Ddata .NE. ST_NOHANDLE)) THEN
    CALL storage_free (rx%h_Ddata)
  ELSE
    rx%h_Ddata = ST_NOHANDLE
  END IF
  
  ! Clean up the structure
  DO i = 1,rx%nblocks
    CALL lsyssc_releaseVector (rx%RvectorBlock(i))
  END DO
  rx%NEQ = 0
  rx%nblocks = 0
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE lsysbl_releaseMatrix (rmatrix)
 
!<description>
  ! Releases the memory that is reserved by rmatrix.
!</description>
  
!<inputoutput>
  ! Block matrix to be released
  TYPE(t_matrixBlock), INTENT(INOUT)                :: rMatrix
!</inputoutput>

!</subroutine>
  
  ! local variables
  INTEGER :: x,y
  
  ! loop through all the submatrices and release them
  DO x=1,rmatrix%ndiagBlocks
    DO y=1,rmatrix%ndiagBlocks
      
      ! Only release the matrix if there is one.
      IF (rmatrix%RmatrixBlock(y,x)%NA .NE. 0) THEN
        CALL lsyssc_releaseMatrix (rmatrix%RmatrixBlock(y,x))
      END IF
      
    END DO
  END DO
  
  ! Clean up the other variables, finish
  rmatrix%NEQ = 0
  rmatrix%ndiagBlocks = 0
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE lsysbl_blockMatVec (rmatrix, rx, ry, cx, cy)
 
!<description>
  ! Performs a matrix vector multiplicationwith a given scalar matrix:
  !    $$ Dy   =   cx * rMatrix * rx   +   cy * ry $$
  ! Vector and matrix must be compatible to each other (same size, sorting 
  ! strategy,...).
!</description>
  
!<input>
  
  ! Block matrix
  TYPE(t_matrixBlock), INTENT(IN)                   :: rmatrix

  ! Block vector to multiply with the matrix.
  TYPE(t_vectorBlock), INTENT(IN)                   :: rx
  
  ! Multiplicative factor for rx
  REAL(DP), INTENT(IN)                              :: cx

  ! Multiplicative factor for ry
  REAL(DP), INTENT(IN)                              :: cy
  
!</input>

!<inputoutput>
  
  ! Additive vector. Receives the result of the matrix-vector multiplication
  TYPE(t_vectorBlock), INTENT(INOUT)                :: ry
  
!</inputoutput>

!</subroutine>
  
  ! local variables
  INTEGER :: x,y,mvok
  REAL(DP)     :: cyact
  
  ! The vectors must be compatible to each other.
  CALL lsysbl_isVectorCompatible (rx,ry)

  ! and compatible to the matrix
  CALL lsysbl_isMatrixCompatible (rx,rmatrix)

  ! loop through all the sub matrices and multiply.
  DO y=1,rmatrix%ndiagBlocks
  
    ! The original vector ry must only be respected once in the first
    ! matrix-vector multiplication. The additive contributions of
    ! the other blocks must simply be added to ry without ry being multiplied
    ! by cy again!
  
    cyact = cy
    MVOK = NO
    DO x=1,rMatrix%ndiagBlocks
      
      ! Only call the MV when there is a scalar matrix that we can use!
      IF (rMatrix%RmatrixBlock(y,x)%NA .NE. 0) THEN
        CALL lsyssc_scalarMatVec (rMatrix%RmatrixBlock(y,x), rx%RvectorBlock(x), &
                                  ry%RvectorBlock(x), cx, cyact)
        cyact = 1.0_DP
        mvok = YES
      END IF
      
    END DO
    
    ! If mvok=NO here, there was no matrix in the current line!
    ! A very special case, but nevertheless may happen. In this case,
    ! simply scale the vector ry by cyact!
    
    IF (mvok .EQ.NO) THEN
      CALL lsysbl_vectorScale (ry,cy)
    END IF
    
  END DO
   
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_vectorCopy (rx,ry)
  
!<description>
  ! Copies vector data: ry = rx.
  ! Both vectors must have the same size. All structural data of rx is
  ! transferred to ry, so rx and ry are compatible to each other afterwards.
!<input>
  
  ! Source vector
  TYPE(t_vectorBlock),INTENT(IN) :: rx
  
!</input>

!<inputoutput>
  
  ! Destination vector
  TYPE(t_vectorBlock),INTENT(INOUT) :: ry
  
  !</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: h_Ddata, cdataType
  LOGICAL :: bisCopy
  
  ! First, make a backup of some crucial data so that it does not
  ! get destroyed.
  h_Ddata = ry%h_Ddata
  cdataType = ry%cdataType
  bisCopy = ry%bisCopy
  
  ! Then transfer all structural information of rx to ry.
  ! This automatically makes both vectors compatible to each other.
  ry = rx
  
  ! Restore crucial data
  ry%h_Ddata = h_Ddata
  ry%bisCopy = bisCopy
  
  ! Restore the handle in all subvectors
  ry%RvectorBlock(1:ry%nblocks)%h_Ddata = h_Ddata
  ry%RvectorBlock(1:ry%nblocks)%cdataType = cdataType
  
  ! And finally copy the data. 
  ! Use storage_copy to copy data by handles.
  CALL storage_copy(rx%h_Ddata, ry%h_Ddata)
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_vectorScale (rx,c)
  
!<description>
  ! Scales a vector vector rx: rx = c * rx
!</description>
  
!<inputoutput>
  
  ! Source and destination vector
  TYPE(t_vectorBlock), INTENT(INOUT) :: rx
  
!</inputoutput>

!<input>
  
  ! Multiplication factor
  REAL(DP), INTENT(IN) :: c
  
!</input>
  
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata
  
  ! Taje care of the data type!
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the pointer and scale the whole data array.
    CALL storage_getbase_double(rx%h_Ddata,p_Ddata)
    CALL lalg_vectorScaleDble (p_Ddata,c)  

  CASE (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    CALL storage_getbase_single(rx%h_Ddata,p_Fdata)
    CALL lalg_vectorScaleSngl (p_Fdata,REAL(c,SP))  

  CASE DEFAULT
    PRINT *,'lsysbl_vectorScale: Unsupported data type!'
    STOP
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_vectorClear (rx)
  
!<description>
  ! Clears the block vector dx: Dx = 0
!<inputoutput>
  
  ! Destination vector to be cleared
  TYPE(t_vectorBlock), INTENT(INOUT) :: rx
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource
  REAL(SP), DIMENSION(:), POINTER :: p_Ssource
  
  ! Take care of the data type
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the pointer and scale the whole data array.
    CALL storage_getbase_double(rx%h_Ddata,p_Dsource)
    CALL lalg_vectorClearDble (p_Dsource)
  
  CASE (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    CALL storage_getbase_single(rx%h_Ddata,p_Ssource)
    CALL lalg_vectorClearSngl (p_Ssource)

  CASE DEFAULT
    PRINT *,'lsysbl_vectorClear: Unsupported data type!'
    STOP
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsysbl_vectorLinearComb (rx,ry,cx,cy)
  
!<description>
  ! Performs a linear combination: ry = cx * rx  +  cy * ry
  ! Both vectors must be compatible to each other (same size, sorting 
  ! strategy,...).
!</description>  
  
!<input>
  
  ! First source vector
  TYPE(t_vectorBlock), INTENT(IN)    :: rx
  
  ! Scaling factor for Dx
  REAL(DP), INTENT(IN)               :: cx

  ! Scaling factor for Dy
  REAL(DP), INTENT(IN)               :: cy
  
!</input>

!<inputoutput>
  
  ! Second source vector; also receives the result
  TYPE(t_vectorBlock), INTENT(INOUT) :: ry
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource, p_Ddest
  REAL(SP), DIMENSION(:), POINTER :: p_Ssource, p_Sdest
  
  ! The vectors must be compatible to each other.
  CALL lsysbl_isVectorCompatible (rx,ry)

  IF (rx%cdataType .NE. ry%cdataType) THEN
    PRINT *,'lsysbl_vectorLinearComb: different data types not supported!'
    STOP
  END IF
  
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the pointers and copy the whole data array.
    CALL storage_getbase_double(rx%h_Ddata,p_Dsource)
    CALL storage_getbase_double(ry%h_Ddata,p_Ddest)
    
    CALL lalg_vectorLinearCombDble (p_Dsource,p_Ddest,cx,cy)

  CASE (ST_SINGLE)
    ! Get the pointers and copy the whole data array.
    CALL storage_getbase_single(rx%h_Ddata,p_Ssource)
    CALL storage_getbase_single(ry%h_Ddata,p_Sdest)
    
    CALL lalg_vectorLinearCombSngl (p_Ssource,p_Sdest,cx,cy)
  
  CASE DEFAULT
    PRINT *,'lsysbl_vectorLinearComb: Unsupported data type!'
    STOP
  END SELECT

  END SUBROUTINE
  
  !****************************************************************************
!<function>
  
  REAL(DP) FUNCTION lsysbl_scalarProduct (rx, ry)
  
!<description>
  ! Calculates a scalar product of two block vectors.
  ! Both vectors must be compatible to each other (same size, sorting 
  ! strategy,...).
!</description>
  
!<input>
  ! First vector
  TYPE(t_vectorBlock), INTENT(IN)                  :: rx

  ! Second vector
  TYPE(t_vectorBlock), INTENT(IN)                  :: ry

!</input>

!<result>
  ! The scalar product <rx,ry> of the two block vectors.
!</result>

!</function>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: h_Ddata1dp
  REAL(DP), DIMENSION(:), POINTER :: h_Ddata2dp
  INTEGER(PREC_VECIDX) :: i
  REAL(DP) :: res
  
  ! The vectors must be compatible to each other.
  CALL lsysbl_isVectorCompatible (rx,ry)

  ! Is there data at all?
  res = 0.0_DP
  
  IF ( (rx%NEQ .EQ. 0) .OR. (ry%NEQ .EQ. 0) .OR. (rx%NEQ .NE. rx%NEQ)) THEN
    PRINT *,'Error in lsyssc_scalarProduct: Vector dimensions wrong!'
    STOP
  END IF
  
  ! Take care of the data type before doing a scalar product!
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    
    CALL storage_getbase_double (rx%h_Ddata,h_Ddata1dp)
    
    SELECT CASE (ry%cdataType)
    CASE (ST_DOUBLE)
      ! Get the data arrays
      CALL storage_getbase_double (ry%h_Ddata,h_Ddata2dp)
      
      ! Perform the scalar product
      res = 0.0_DP
      DO i=1,rx%NEQ
        res = res + h_Ddata1dp(i)*h_Ddata2dp(i)
      END DO
      
    CASE DEFAULT
      PRINT *,'lsyssc_scalarProduct: Not supported precision combination'
      STOP
    END SELECT
    
  CASE DEFAULT
    PRINT *,'lsyssc_scalarProduct: Not supported precision combination'
    STOP
  END SELECT
  
  ! Return the scalar product, finish
  lsysbl_scalarProduct = res

  END FUNCTION

  !****************************************************************************
!<subroutine>
  
  REAL(DP) FUNCTION lsysbl_vectorNorm (rx,cnorm,iposMax)
  
!<description>
  ! Calculates the norm of a vector. cnorm specifies the type of norm to
  ! calculate.
!</description>
  
!<input>
  ! Vector to calculate the norm of.
  TYPE(t_vectorBlock), INTENT(IN)                  :: rx

  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  INTEGER, INTENT(IN) :: cnorm
!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  INTEGER(I32), INTENT(OUT), OPTIONAL :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! < 0, if an error occurred (unknown norm).
!</result>

!<result>
  ! The scalar product <rx,ry> of the two block vectors.
!</result>

!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata

  ! Is there data at all?
  IF (rx%h_Ddata .EQ. ST_NOHANDLE) THEN
    PRINT *,'Error in lsysbl_vectorNorm: Vector empty!'
    STOP
  END IF
  
  ! Take care of the data type before doing a scalar product!
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the array and calculate the norm
    CALL storage_getbase_double (rx%h_Ddata,p_Ddata)
    lsysbl_vectorNorm = lalg_normDble (p_Ddata,cnorm,iposMax) 
    
  CASE (ST_SINGLE)
    ! Get the array and calculate the norm
    CALL storage_getbase_single (rx%h_Ddata,p_Fdata)
    lsysbl_vectorNorm = lalg_normSngl (p_Fdata,cnorm,iposMax) 
    
  CASE DEFAULT
    PRINT *,'lsysbl_vectorNorm: Unsupported data type!'
    STOP
  END SELECT
  
  END FUNCTION

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_createMatFromScalar (rscalarMat,rmatrix)
  
!<description>
  ! This routine creates a 1x1 block matrix rmatrix from a scalar matrix
  ! rscalarMat. Both, rscalarMat and rmatrix will share the same handles,
  ! so changing the content of rmatrix will change rscalarMat, too.
  ! Therefore, the imatrixSpec flag of the submatrix in the block matrix
  ! will be set to LSYSSC_MSPEC_ISCOPY to indicate that the matrix is a
  ! copy of another one.
!</description>
  
!<input>
  ! The scalar matrix which should provide the data
  TYPE(t_matrixScalar), INTENT(IN) :: rscalarMat
!</input>

!<output>
  ! The 1x1 block matrix, created from rscalarMat.
  TYPE(t_matrixBlock), INTENT(OUT) :: rmatrix
!</output>
  
!</subroutine>

    ! Fill the rmatrix structure with data.
    rmatrix%NEQ         = rscalarMat%NEQ
    rmatrix%ndiagBlocks = 1
    rmatrix%imatrixSpec = LSYSBS_MSPEC_SCALAR
    
    ! Copy the content of the scalar matrix structure into the
    ! first block of the block matrix
    rmatrix%RmatrixBlock(1,1)             = rscalarMat

    ! The matrix is a copy of another one. Note this!
    rmatrix%RmatrixBlock(1,1)%imatrixSpec = LSYSSC_MSPEC_ISCOPY
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_createVecFromScalar (rscalarVec,rvector)
  
!<description>
  ! This routine creates a 1 block vector rvector from a scalar vector
  ! rscalarVec. Both, rscalarVec and rvector will share the same handles,
  ! so changing the content of rvector will change rscalarVec, too.
  ! Therefore, the bisCopy flag of the subvector in the block vector
  ! will be set to TRUE.
!</description>
  
!<input>
  ! The scalar vector which should provide the data
  TYPE(t_vectorScalar), INTENT(IN) :: rscalarVec
!</input>

!<output>
  ! The 1x1 block matrix, created from rscalarMat.
  TYPE(t_vectorBlock), INTENT(OUT) :: rvector
!</output>
  
!</subroutine>

    ! Fill the rvector structure with data.
    rvector%NEQ         = rscalarVec%NEQ
    rvector%cdataType   = rscalarVec%cdataType
    rvector%h_Ddata     = rscalarVec%h_Ddata
    rvector%nblocks     = 1
    
    ! Copy the content of the scalar matrix structure into the
    ! first block of the block vector
    rvector%RvectorBlock(1)             = rscalarVec

    ! The handle belongs to another vector - note this in the structure.
    rvector%RvectorBlock(1)%bisCopy     = .TRUE.
    
    ! The data of the vector actually belongs to another one - to the
    ! scalar one. Note this in the structure.
    rvector%bisCopy = .TRUE.
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_enforceStructure (rtemplateVec,rvector)
  
!<description>
  ! This routine enforces the structure of the vector rtemplale
  ! in the vector rvector. 
  !
  ! WARNING: This routine should be used with care if you knwo
  !          what you are doing !!!
  ! All structural data of rtemplate (boundary conditions, spatial 
  ! discretisation, size,sorting,...) are copied from rtemplate to 
  ! rvector without checking the compatibility!!!
  !
  ! The only check in this routine is that the memory rvector provides
  ! must be at least as large as rtemplate%NEQ; otherwise an error
  ! is thrown. The data type of rvector is also not changed.
!</description>
  
!<input>
  ! A template vector specifying structural data
  TYPE(t_vectorBlock), INTENT(IN) :: rtemplateVec
!</input>

!input<output>
  ! The vector which structural data should be overwritten
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: cdata,h_Ddata
  INTEGER(PREC_VECIDX) :: length
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata

    ! Only basic check: there must be enough memory
    SELECT CASE (rvector%cdataType)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double (rvector%h_Ddata,p_Ddata)
      length = SIZE(p_Ddata)
    CASE (ST_SINGLE)
      CALL storage_getbase_double (rvector%h_Ddata,p_Ddata)
      length = SIZE(p_Ddata)
    CASE DEFAULT
      PRINT *,'lsysbl_enforceStructure: Unsupported data type'
      STOP
    END SELECT
    
    IF (length .LT. rtemplateVec%NEQ) THEN
      PRINT *,'lsysbl_enforceStructure: Destination vector too small!'
      STOP
    END IF
  
    ! Get data type and handle from rvector
    cdata = rvector%cdataType
    h_Ddata = rvector%h_Ddata
    
    ! Overwrite rvector
    rvector = rtemplateVec
    
    ! Restore the data array
    rvector%cdataType = cdata
    rvector%h_Ddata = h_Ddata
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_updateMatStrucInfo (rmatrix)
  
!<description>
  ! Updates structural information of a block matrix by recalculation
  ! from submatrices.
  ! This routine must be called, if the application changes structural
  ! information in one or more of the matrix blocks of a submatrix.
  ! In this case, this routine recalculates crucial information
  ! of the global matrix (number of blocks, global NEQ,...) by checking 
  ! the submatrices.
!</description>
  
!<inputoutput>
  ! The block matrix which is to be updated.
  TYPE(t_matrixBlock), INTENT(INOUT) :: rmatrix
!</inputoutput>

!</subroutine>

  INTEGER :: i,j,NEQ
  
  ! At first, recalculate the number of diagonal blocks in the
  ! block matrix
  DO i=LSYSBL_MAXBLOCKS,1,-1
    IF (rmatrix%RmatrixBlock(i,i)%NA .NE. 0) EXIT
  END DO
  
  ! i is the new number of diagonal blocks
  rmatrix%ndiagBlocks = i
  
  ! Calculate the new NEQ
  NEQ = 0
  DO j=1,i
    NEQ = NEQ + rmatrix%RmatrixBlock(i,i)%NEQ
  END DO
  rmatrix%NEQ = NEQ

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_sortVectorInSitu (rvector,rtemp,bsort)
  
!<description>
  ! This routine sorts a block vector or unsorts it.
  ! If bsort=TRUE, the vector is sorted, otherwise it's unsorted.
  !
  ! The sorting uses the associated permutation of every subvector,
  ! so before calling this routine, a permutation should be assigned
  ! to every subvector, either manually or using lsysbl_setSortStrategy.
  ! The associated sorting strategy tag will change to
  !  + |subvector%isortStrategy|  - if bsort = TRUE
  !  - |subvector%isortStrategy|  - if bsort = false
  ! so the absolute value of rvector%isortStrategy indicates the sorting
  ! strategy and the sign determines whether the (sub)vector is
  ! actually sorted for the associated sorting strategy.
!</description>

!<input>
  ! Whether to sort the vector (TRUE) or sort it back to unsorted state
  ! (FALSE).
  LOGICAL, INTENT(IN) :: bsort
!</input>
  
!<inputoutput>
  ! The vector which is to be resorted
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector

  ! A scalar temporary vector. Must be of the same data type as rvector.
  ! Must be at least as large as the longest subvector in rvector.
  TYPE(t_vectorScalar), INTENT(INOUT) :: rtemp
!</inputoutput>

!</subroutine>

  INTEGER :: iblock
  
  ! Loop over the blocks
  IF (bsort) THEN
    DO iblock = 1,rvector%nblocks
      CALL lsyssc_sortVectorInSitu (&
          rvector%RvectorBlock(iblock), rtemp,&
          ABS(rvector%RvectorBlock(iblock)%isortStrategy))
    END DO
  ELSE
    DO iblock = 1,rvector%nblocks
      CALL lsyssc_sortVectorInSitu (&
          rvector%RvectorBlock(iblock), rtemp,&
          -ABS(rvector%RvectorBlock(iblock)%isortStrategy))
    END DO
  END IF

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsysbl_setSortStrategy (rvector,IsortStrategy,Hpermutations)
  
!<description>
  ! This routine simultaneously connects all subvectors of a block
  ! vector with a sorting strategy. IsortStrategy is an array of tags
  ! that identify the sorting strategy of each subvector. Hpermutation
  ! is an array of handles. Each handle identifies the permutation
  ! that is to assign to the corresponding subvector.
  !
  ! The routine does not resort the vector. Only the identifier in
  ! IsortStrategy and the handle of the permutation are saved to the
  ! vector. To indicate that the subvector is not sorted, the negative
  ! value of IsortStrategy is saved to subvector%isortStrategy.
!</description>

!<inputoutput>
  ! The vector which is to be resorted
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!</inputoutput>

!<input>
  ! An array of sorting strategy identifiers (SSTRAT_xxxx), each for one
  ! subvector of the global vector. 
  ! The negative value of this identifier is saved to the corresponding
  ! subvector.
  INTEGER, DIMENSION(rvector%nblocks), INTENT(IN) :: IsortStrategy
  
  ! An array of handles. Each handle corresponds to a subvector and
  ! defines a permutation how to resort the subvector.
  ! Each permutation associated to a handle must be of the form
  !    array [1..2*NEQ] of integer  (NEQ=NEQ(subvector))
  ! with entries (1..NEQ)       = permutation
  ! and  entries (NEQ+1..2*NEQ) = inverse permutation.
  INTEGER, DIMENSION(rvector%nblocks), INTENT(IN) :: Hpermutations

!</input>
  
!</subroutine>

  ! Install the sorting strategy in every block. 
  rvector%RvectorBlock%isortStrategy = -ABS(IsortStrategy)
  rvector%RvectorBlock%h_IsortPermutation = Hpermutations

  END SUBROUTINE

  !****************************************************************************
!<function>
  
  LOGICAL FUNCTION lsysbl_isMatrixSorted (rmatrix)
  
!<description>
  ! Loops through all submatrices in a block matrix and checks if any of
  ! them is sorted. If yes, the block matrix is given the state 'sorted'.
  ! If not, the matrix is stated 'unsorted'
!</description>
  
!<input>
  ! Block matrix to check
  TYPE(t_matrixBlock), INTENT(IN)                  :: rmatrix
!</input>

!<result>
  ! Whether the matrix is sorted or not.
!</result>

!</function>

  INTEGER(PREC_VECIDX) :: i,j,k
  
  ! Loop through all matrices and get the largest sort-identifier
  k = 0
  DO j=1,rmatrix%ndiagBlocks
    DO i=1,rmatrix%ndiagBlocks
      IF (rmatrix%RmatrixBlock(i,j)%NEQ .NE. 0) THEN
        k = MAX(k,rmatrix%RmatrixBlock(i,j)%isortStrategy)
      END IF
    END DO
  END DO
  
  ! If k is greater than 0, at least one submatrix is sorted
  lsysbl_isMatrixSorted = k .GT. 0

  END FUNCTION

  !****************************************************************************
!<function>
  
  LOGICAL FUNCTION lsysbl_isVectorSorted (rvector)
  
!<description>
  ! Loops through all subvectors in a block vector and checks if any of
  ! them is sorted. If yes, the block vector is given the state 'sorted'.
  ! If not, the vector is stated 'unsorted'
!</description>
  
!<input>
  ! Block matrix to check
  TYPE(t_vectorBlock), INTENT(IN)                  :: rvector
!</input>

!<result>
  ! Whether the vector is sorted or not.
!</result>

!</function>

  INTEGER(PREC_VECIDX) :: i,k
  
  ! Loop through all matrices and get the largest sort-identifier
  k = 0
  DO i=1,rvector%nblocks
    IF (rvector%RvectorBlock(i)%NEQ .NE. 0) THEN
      k = MAX(k,rvector%RvectorBlock(i)%isortStrategy)
    END IF
  END DO
  
  ! If k is greater than 0, at least one submatrix is sorted
  lsysbl_isVectorSorted = k .GT. 0

  END FUNCTION

END MODULE

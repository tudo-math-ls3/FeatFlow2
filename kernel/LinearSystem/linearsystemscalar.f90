!##############################################################################
!# ****************************************************************************
!# <name> linearsystemscalar </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic definitions and routines to maintain scalar
!# linear systems. A linear system is realised by matrices and vectors, where
!# the matrix does not decompose into smaller 'block-matrices'. Therefore,
!# we denote the matrices/vectors here as 'scalar', thus containing no 
!# ''blocks''.
!#
!# For storing matrices, we use the CSR format, realised as ''Format-7''
!# and ''Format-9'' matrices. Format-9 matrices describe the usual CSR-format
!# with the extension, that an index is stored for each line to the first
!# element on/above the diagonal. Format-7 is a variant of the CSR-format
!# used in old FEAT programs, where the diagonal element is stored first in
!# each line.
!#
!# The following routines can be found in this module:
!#
!#  1.) lsyssc_createVector
!#      -> Create a simple scalar vector of length NEQ
!#
!#  2.) lsyssc_createVecByDiscr
!#      -> Create a vector based on a scalar discretisation structure
!#
!#  3.) lsyssc_scalarProduct
!#      -> Calculate the scalar product of two vectors
!#
!#  4.) lsyssc_scalarMatVec
!#      -> Multiply a scalar matrix with a scalar vector
!#
!#  5.) lsyssc_releaseMatrix
!#      -> Release a scalar matrix from memory.
!#
!#  6.) lsyssc_releaseVector
!#      -> Release a scalar vector from memory.
!#
!#  7.) lsyssc_duplicateMatrix
!#      -> Create a duplicate of a given matrix or matrix-structure
!#
!#  8.) lsyssc_duplicateVector
!#      -> Create a duplicate of a given vector
!#
!#  9.) lsyssc_sortVectorInSitu
!#      -> Resort the entries of a vector or unsort them
!#
!# 10.) lsyssc_vectorActivateSorting
!#      -> Resort the entries of a vector or unsort them according to
!#         a previously attached sorting strategy
!#
!# 11.) lsyssc_sortMatrix
!#      -> Resort the entries of a matrix or unsort them
!#
!# 12.) lsyssc_isVectorCompatible
!#      -> Checks whether two vectors are compatible to each other
!#
!# 13.) lsyssc_isMatrixCompatible
!#      -> Checks whether a matrix and a vector are compatible to each other
!#
!# 14.) lsyssc_getbase_double
!#      -> Get a pointer to the double precision data array of the vector
!#
!# 15.) lsyssc_getbase_single
!#      -> Get a pointer to the single precision data array of the vector
!#
!# 16.) lsyssc_addIndex
!#      -> Auxiliary routine. Adds an integer to each elememt of an integer 
!#         array.
!#
!# 17.) lsyssc_vectorNorm
!#      -> Calculate the norm of a vector.
!#
!# 18.) lsyssc_invertedDiagMatVec
!#      -> Multiply a vector with the inverse of the diagonal of a scalar
!#         matrix
!#
!# 19.) lsyssc_clearMatrix
!#      -> Clears a matrix, i.e. overwrites all entries with 0.0
!#
!# 20.) lsyssc_convertMatrix
!#      -> Allows to convert a matrix to another matrix structure.
!#
!# 21.) lsyssc_copyVector
!#       -> Copy a vector over to another one
!#
!# 22.) lsyssc_scaleVector
!#      -> Scale a vector by a constant
!#
!# 23.) lsyssc_clearVector
!#      -> Clear a vector
!#
!# 24.) lsyssc_vectorLinearComb
!#      -> Linear combination of two vectors
!#
!# 25.) lsyssc_copyMatrix
!#      -> Copies a matrix to another one provided that they have the same 
!#         structure.
!#
!# 26.) lsyssc_transposeMatrix
!#      -> Transposes a scalar matrix.
!#
!# Sometimes useful auxiliary routines:
!#
!# 1.) lsyssc_rebuildKdiagonal (Kcol, Kld, Kdiagonal, neq)
!#     -> Rebuild the Kdiagonal array in a matrix of format 9
!# </purpose>
!##############################################################################

MODULE linearsystemscalar

  USE fsystem
  USE storage
  USE spatialdiscretisation
  USE discretebc
  
  IMPLICIT NONE

!<constants>

!<constantblock description="Global format flags for matrices">

  ! Unidentified matrix format
  INTEGER, PARAMETER :: LSYSSC_MATRIXUNDEFINED = 0
  
  ! Identifier for matrix format 1 - full matrix
  INTEGER, PARAMETER :: LSYSSC_MATRIX1 = 1
  
  ! Identifier for matrix format 9 - CSR
  INTEGER, PARAMETER :: LSYSSC_MATRIX9 = 9

  ! Identifier for matrix format 7 - CSR with diagonal element in front
  INTEGER, PARAMETER :: LSYSSC_MATRIX7 = 7

!</constantblock>

!<constantblock description="Flags for the matrix specification bitfield">

  ! Standard matrix
  INTEGER, PARAMETER :: LSYSSC_MSPEC_STANDARD =        0
  
  ! Matrix structure is a copy of another matrix, shared via the same
  ! handles. 
  INTEGER, PARAMETER :: LSYSSC_MSPEC_STRUCTUREISCOPY = 2**0

  ! Matrix content is a copy of another matrix, shared via the same
  ! handles. 
  INTEGER, PARAMETER :: LSYSSC_MSPEC_CONTENTISCOPY   = 2**1
  
  ! Complete matrix is duplicate of another matrix and shares structure
  ! and entries via the same pointers
  INTEGER, PARAMETER :: LSYSSC_MSPEC_ISCOPY = LSYSSC_MSPEC_STRUCTUREISCOPY +&
                                              LSYSSC_MSPEC_CONTENTISCOPY

  ! Matrix is saved transposed.
  ! To use a matrix in a transposed way, the application has to
  ! 1.) set this flag in the imatrixSpec bitfield
  ! 2.) exchange the values in t_matrixScalar%NEQ and t_matrixScalar%NCOLS 
  !     of the matrix structure.
  INTEGER, PARAMETER :: LSYSSC_MSPEC_TRANSPOSED =      2**2

  ! Matrix not present in memory
  INTEGER, PARAMETER :: LSYSSC_MSPEC_NOTINMEMORY =     2**3
  
!</constantblock>

!<constantblock description="KIND values for matrix/vector data">
  
  ! kind value for indices in matrices
  INTEGER, PARAMETER :: PREC_MATIDX = I32

  ! kind value for indices in vectors
  INTEGER, PARAMETER :: PREC_VECIDX = I32

  ! kind value for precision that should be used in matrices
  INTEGER, PARAMETER :: PREC_MATRIX = DP

!</constantblock>

!<constantblock description="Constants for duplicating a matrix">
  
  ! Duplicate by ownership. What belongs to the matrix is duplicated,
  ! what belongs to another matrix is left as-is.
  !INTEGER, PARAMETER :: LSYSSC_DUP_ASIS      = 0

  ! Duplicate nothing, simply copy the structure and mark the handles
  ! as belonging to another matrix.
  INTEGER, PARAMETER :: LSYSSC_DUP_NONE      = 1

  ! Copy the handles of the structure and mark them as belonging to
  ! another matrix. No content is created.
  INTEGER, PARAMETER :: LSYSSC_DUP_STRNOCONT = 2

  ! Duplicate the content, share the structure
  INTEGER, PARAMETER :: LSYSSC_DUP_CONTENT   = 3

  ! Duplicate the structure, share the content
  INTEGER, PARAMETER :: LSYSSC_DUP_STRUCTURE = 4

  ! Duplicate both, structure and content
  INTEGER, PARAMETER :: LSYSSC_DUP_ALL       = 10
  
  ! Allocate memory for the matrix structure in the same size as the original
  ! matrix. No content is created.
  INTEGER, PARAMETER :: LSYSSC_DUP_EMPTYSTRUC = 6

  ! Allocate memory for the matrix structure and content in the same size 
  ! as the original matrix.
  INTEGER, PARAMETER :: LSYSSC_DUP_EMPTYALL   = 9
  
  ! ------------
  
  ! Don't set up the content/structure of the destination matrix, ignore
  ! any previous structure/content
  INTEGER, PARAMETER :: LSYSSC_DUP_IGNORE = 0
  
  ! Removes any existing matrix content/structure from the destination matrix.
  ! Releases memory if necessary.
  INTEGER, PARAMETER :: LSYSSC_DUP_REMOVE = 1
  
  ! Removes any existing matrix content from the destination matrix.
  ! No memory is released, handles are simply dismissed.
  INTEGER, PARAMETER :: LSYSSC_DUP_DISMISS = 2
  
  ! The destination matrix recveives the same handles for matrix content/structure
  ! as the source matrix  and therefore shares the same content/structure.
  INTEGER, PARAMETER :: LSYSSC_DUP_SHARE = 3
  
  ! The destination matrix gets a copy of the content of rsourceMatrix.
  ! If necessary, new memory is allocated.
  ! If the destination matrix  already contains allocated memory, content/structure
  ! data is simply copied from rsourceMatrix into that.
  INTEGER, PARAMETER :: LSYSSC_DUP_COPY = 4
  
  ! Duplicate by ownership. What belongs to the source matrix is copied 
  ! (the same as LSYSSC_DUP_COPY). What belongs even to another matrix than
  ! the source matrix is shared (the same as LSYSSC_DUP_SHARE, .
  INTEGER, PARAMETER :: LSYSSC_DUP_ASIS = 5
  
  ! New memory is allocated for the structure/content in the same size as 
  ! in tzhe source matrix but no data is copied; the arrays are left uninitialised.
  INTEGER, PARAMETER :: LSYSSC_DUP_EMPTY = 6 
                                               
  ! Copy the basic matrix information but don't copy handles.
  ! Set all handles of dynamic information to ST_NOHANDLE, so the matrix
  ! 'looks like' the old but has no dynamic data associated.
  INTEGER, PARAMETER :: LSYSSC_DUP_TEMPLATE   = 7
                 
!</constantblock>

!<constantblock description="Constants for transposing a matrix">
  
  ! Transpose the full matrix
  INTEGER, PARAMETER :: LSYSSC_TR_ALL       = 0

  ! Transpose only the matrix structure
  INTEGER, PARAMETER :: LSYSSC_TR_STRUCTURE = 1   

  ! Don't transpose the matrix, simply mark the matrix as transposed
  ! by changing the flag in imatrixSpec
  INTEGER, PARAMETER :: LSYSSC_TR_VIRTUAL    = 2

  ! Don't transpose the matrix. Copy the matrix in memory and mark the 
  ! matrix as transposed by changing the flag in imatrixSpec
  INTEGER, PARAMETER :: LSYSSC_TR_VIRTUALCOPY = 3

!</constantblock>


!</constants>

!<types>
!<typeblock>
  
  ! A scalar vector that can be used by scalar linear algebra routines.
  
  TYPE t_vectorScalar
  
    ! Length of the vector; not necessarily = SIZE(p_Ddata), as this
    ! scalar vector might be a subvector of a larger vector allocated
    ! on the heap!
    INTEGER(PREC_VECIDX) :: NEQ       = 0
    
    ! Data type of the entries in the vector. Either ST_SINGLE or
    ! ST_DOUBLE.
    INTEGER         :: cdataType = ST_DOUBLE
    
    ! Handle identifying the vector entries.
    ! = ST_NOHANDLE if not allocated on the global heap. 
    INTEGER         :: h_Ddata   = ST_NOHANDLE
    
    ! Flag whether or not the vector is resorted.
    !  <0: Vector is unsorted, sorting strategy is prepared in 
    !      h_IsortPermutation for a possible resorting of the entries.
    !  =0: Vector is unsorted, no sorting strategy attached.
    !  >0: Vector is sorted according to a sorting strategy.
    ! The value identifies the sorting strategy used; this is
    ! usually one of the SSTRAT_xxxx constants from the module
    ! 'sortstrategy'.
    ! If <> 0, the absolute value 
    !               |isortStrategy| > 0
    ! indicates the sorting strategy to use, while the sign
    ! indicates whether the sorting strategy is active on the
    ! vector (+) or not (-).
    ! The value is usually one of the SSTRAT_xxxx constants from
    ! the module 'sortstrategy'.
    INTEGER         :: isortStrategy   = 0
    
    ! Handle to renumbering strategy for resorting the vector.
    ! The renumbering strategy is a vector
    !   array [1..2*NEQ] of integer
    ! The first NEQ entries (1..NEQ) represent the permutation how to
    ! sort an unsorted vector. The second NEQ entries (NEQ+1..2*NEQ)
    ! represent the inverse permutation.
    ! Whether or not the vector is actually sorted depends on the
    ! flag isortStrategy!
    INTEGER         :: h_IsortPermutation = ST_NOHANDLE
    
    ! Start position of the vector data in the array identified by
    ! h_Ddata. Normally = 1. Can be set to > 1 if the vector is a subvector
    ! in a larger memory block allocated on the heap.
    INTEGER(PREC_VECIDX) :: iidxFirstEntry = 1
    
    ! Is set to true, if the handle h_Ddata belongs to another vector,
    ! i.e. when this vector shares data with another vector.
    ! This is usually used for block vectors containing a couple of
    ! scalar subvectors; in this case, the block vector is the actual
    ! 'owner' of the handle!
    LOGICAL              :: bisCopy        = .FALSE.
    
    ! A pointer to the spatial discretisation or NULL(), if the vector is
    ! just an array, not belonging to any discretisation.
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscretisation => NULL()
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! A scalar matrix that can be used by scalar linear algebra routines.
  
  TYPE t_matrixScalar
    
    ! Format-tag. Identifies the format of the matrix. 
    ! Can tage one of the LSYSSC_MATRIXx format flags.
    INTEGER      :: cmatrixFormat = LSYSSC_MATRIXUNDEFINED
    
    ! Matrix specification tag. This is a bitfield coming from an OR 
    ! combination of different LSYSSC_MSPEC_xxxx constants and specifies
    ! various details of the matrix. If it is =LSYSSC_MSPEC_STANDARD,
    ! the matrix is a usual matrix that needs no special handling.
    INTEGER(I32) :: imatrixSpec = LSYSSC_MSPEC_STANDARD
    
    ! Number of elements in the matrix
    INTEGER(PREC_MATIDX) :: NA  = 0
    
    ! Number of equations = rows in the matrix.
    ! Remark: When the application sets the LSYSSC_MSPEC_TRANSPOSED flag
    !  in imatrixSpec to indicate a transposed matrix, this has to be 
    !  exchanged with NCOLS to indicate the number of equations in the
    !  transposed matrix!
    INTEGER(PREC_VECIDX) :: NEQ = 0

    ! Number of columns in the matrix.
    ! Remark: When the application sets the LSYSSC_MSPEC_TRANSPOSED flag
    !  in imatrixSpec to indicate a transposed matrix, this has to be 
    !  exchanged with NROWS to indicate the number of columns in the
    !  transposed matrix!
    INTEGER(PREC_VECIDX) :: NCOLS = 0
    
    ! Multiplier for matrix entries. All entries in the matrix are
    ! scaled by this multiplier when doing Matrix-vector multiplication.
    REAL(DP)         :: dscaleFactor = 1.0_DP
    
    ! Flag whether or not the matrix is resorted.
    !  <0: Matrix is unsorted, sorting strategy is prepared in 
    !      h_IsortPermutation for a possible resorting of the entries.
    !  =0: Matrix is unsorted, no sorting strategy attached.
    !  >0: Matrix is sorted according to a sorting strategy.
    ! If <> 0, the absolute value 
    !               |isortStrategy| > 0
    ! indicates the sorting strategy to use, while the sign
    ! indicates whether the sorting strategy is active on the
    ! matrix (+) or not (-).
    ! The value is usually one of the SSTRAT_xxxx constants from
    ! the module 'sortstrategy'.
    INTEGER         :: isortStrategy   = 0
    
    ! Handle to renumbering strategy for resorting the matrix.
    ! The renumbering strategy is a vector
    !   array [1..2*NEQ] of integer
    ! The first NEQ entries (1..NEQ) represent the permutation how to
    ! sort an unsorted matrix. The second NEQ entries (NEQ+1..2*NEQ)
    ! represent the inverse permutation.
    ! Whether or not the matrix is actually sorted depends on the
    ! flag isortStrategy!
    INTEGER         :: h_IsortPermutation = ST_NOHANDLE
    
    
    ! Data type of the entries in the vector. Either ST_SINGLE or
    ! ST_DOUBLE.
    INTEGER         :: cdataType = ST_DOUBLE
    
    ! Format-7 and Format-9: Handle identifying the elements in the matrix
    !REAL(PREC_MATRIX), DIMENSION(:), POINTER       :: DA         => NULL()
    INTEGER                    :: h_DA = ST_NOHANDLE
    
    ! Format-7 and Format-9: Handle identifying the column structure
    !INTEGER(PREC_MATIDX), DIMENSION(:), POINTER    :: KCOL       => NULL()
    INTEGER                    :: h_Kcol = ST_NOHANDLE
    
    ! Format-7 and Format-9: Handle identifying the row structure
    !INTEGER, DIMENSION(:), POINTER            :: KLD        => NULL()
    INTEGER                    :: h_Kld = ST_NOHANDLE
    
    ! Format-9: Similar to row structure. Handle top array of length NEQ.
    ! For each row, pointer to the first element on the upper
    ! triangular part of the matrix. For matrices with elements on the
    ! diagonal, this points to the diagonal element in each row.
    ! For matrices with no elements in the upper right part of the
    ! matrix, this points to the first element on the next line.
    !INTEGER, DIMENSION(:), POINTER            :: Kdiagonal  => NULL()
    INTEGER                    :: h_Kdiagonal = ST_NOHANDLE
    
    ! A pointer to the spatial discretisation
    TYPE(t_spatialDiscretisation), POINTER     :: p_rspatialDiscretisation => NULL()
    
  END TYPE
  
!</typeblock>

!</types>

CONTAINS

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_isVectorCompatible (rvector1,rvector2,bcompatible)
  
!<description>
  ! Checks whether two vectors are compatible to each other, i.e. share
  ! the same structure, size and sorting strategy.
!</description>

!<input>
  ! The first vector
  TYPE(t_vectorScalar), INTENT(IN) :: rvector1
  
  ! The second vector
  TYPE(t_vectorScalar), INTENT(IN) :: rvector2
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether the vectors are compatible or not.
  ! If not given, an error will inform the user if the two vectors are
  ! not compatible and the program will halt.
  LOGICAL, INTENT(OUT), OPTIONAL :: bcompatible
!</output>

!</subroutine>

  ! We assume that we are not compatible
  IF (PRESENT(bcompatible)) bcompatible = .FALSE.
  
  ! Vectors must have the same size
  IF (rvector1%NEQ .NE. rvector2%NEQ) THEN
    IF (PRESENT(bcompatible)) THEN
      bcompatible = .FALSE.
      RETURN
    ELSE
      PRINT *,'Vectors not compatible, different size!'
      STOP
    END IF
  END IF
    
  ! isortStrategy < 0 means unsorted. Both unsorted is ok.
  
  IF ((rvector1%isortStrategy .GT. 0) .OR. &
      (rvector2%isortStrategy .GT. 0)) THEN

    IF (rvector1%isortStrategy .NE. &
        rvector2%isortStrategy) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Vectors not compatible, differently sorted!'
        STOP
      END IF
    END IF

    IF (rvector1%h_isortPermutation .NE. &
        rvector2%h_isortPermutation) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Vectors not compatible, differently sorted!'
        STOP
      END IF
    END IF
  END IF

  ! Ok, they are compatible
  IF (PRESENT(bcompatible)) bcompatible = .TRUE.

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_isMatrixCompatible (rvector,rmatrix,btransposed,bcompatible)
  
!<description>
  ! Checks whether a vector and a matrix are compatible to each other, i.e. 
  ! share the same structure, size and sorting strategy.
!</description>

!<input>
  ! The vector
  TYPE(t_vectorScalar), INTENT(IN) :: rvector
  
  ! The matrix
  TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
  
  ! Check for rvector being compatible from the left or from the right
  ! to the matrix rmatrix.
  ! =FALSE: Check whether matrix-vector product A*x is possible
  ! =TRUE : Check whether matrix-vector product x^T*A = A^T*x is possible
  LOGICAL, INTENT(IN)              :: btransposed
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether vector and matrix are compatible or not.
  ! If not given, an error will inform the user if the vector/matrix are
  ! not compatible and the program will halt.
  LOGICAL, INTENT(OUT), OPTIONAL :: bcompatible
!</output>

!</subroutine>

  INTEGER(PREC_VECIDX) :: NCOLS

  ! We assume that we are not compatible
  IF (PRESENT(bcompatible)) bcompatible = .FALSE.
  
  IF (.NOT. btransposed) THEN
    NCOLS = rmatrix%NCOLS
  ELSE
    NCOLS = rmatrix%NEQ
  END IF
  
  ! Vector/Matrix must have the same size 
  IF (rvector%NEQ .NE. NCOLS) THEN
    IF (PRESENT(bcompatible)) THEN
      bcompatible = .FALSE.
      RETURN
    ELSE
      PRINT *,'Vector/Matrix not compatible, different block structure!'
      STOP
    END IF
  END IF

  ! isortStrategy < 0 means unsorted. Both unsorted is ok.

  IF ((rvector%isortStrategy .GT. 0) .OR. &
      (rmatrix%isortStrategy .GT. 0)) THEN

    IF (rvector%isortStrategy .NE. &
        rmatrix%isortStrategy) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Vector/Matrix not compatible, differently sorted!'
        STOP
      END IF
    END IF

    IF (rvector%h_isortPermutation .NE. &
        rmatrix%h_isortPermutation) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Vector/Matrix not compatible, differently sorted!'
        STOP
      END IF
    END IF
  END IF

  ! Ok, they are compatible
  IF (PRESENT(bcompatible)) bcompatible = .TRUE.

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_getbase_double (rvector,p_Ddata)
  
!<description>
  ! Returns a pointer to the double precision data array of the vector.
  ! An error is thrown if the vector is not double precision.
!</description>

!<input>
  ! The vector
  TYPE(t_vectorScalar), INTENT(IN) :: rvector
!</input>

!<output>
  ! Pointer to the double precision data array of the vector.
  ! NULL() if the vector has no data array.
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
!</output>

!</subroutine>

  ! Do we have data at all?
 IF ((rvector%NEQ .EQ. 0) .OR. (rvector%h_Ddata .EQ. ST_NOHANDLE)) THEN
   NULLIFY(p_Ddata)
   RETURN
 END IF

  ! Check that the vector is really double precision
  IF (rvector%cdataType .NE. ST_DOUBLE) THEN
    PRINT *,'lsyssc_getbase_double: Vector is of wrong precision!'
    STOP
  END IF

  ! Get the data array
  CALL storage_getbase_double (rvector%h_Ddata,p_Ddata)
  
  ! Modify the starting address/length to get the real array.
  p_Ddata => p_Ddata(rvector%iidxFirstEntry:rvector%iidxFirstEntry+rvector%NEQ-1)
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_getbase_single (rvector,p_Fdata)
  
!<description>
  ! Returns a pointer to the single precision data array of the vector.
  ! An error is thrown if the vector is not single precision.
!</description>

!<input>
  ! The vector
  TYPE(t_vectorScalar), INTENT(IN) :: rvector
!</input>

!<output>
  ! Pointer to the double precision data array of the vector.
  ! NULL() if the vector has no data array.
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata
!</output>

!</subroutine>

  ! Do we have data at all?
 IF ((rvector%NEQ .EQ. 0) .OR. (rvector%h_Ddata .EQ. ST_NOHANDLE)) THEN
   NULLIFY(p_Fdata)
   RETURN
 END IF

  ! Check that the vector is really double precision
  IF (rvector%cdataType .NE. ST_SINGLE) THEN
    PRINT *,'lsyssc_getbase_single: Vector is of wrong precision!'
    STOP
  END IF

  ! Get the data array
  CALL storage_getbase_single (rvector%h_Ddata,p_Fdata)
  
  ! Modify the starting address/length to get the real array.
  p_Fdata => p_Fdata(rvector%iidxFirstEntry:rvector%iidxFirstEntry+rvector%NEQ-1)
  
  END SUBROUTINE

  !****************************************************************************

!<subroutine>
  
  SUBROUTINE lsyssc_createVector (rvector,NEQ,bclear,cdataType)
  
!<description>
  ! This creates a simple scalar vector of length NEQ. Memory is 
  ! allocated on the heap and can be released by lsyssc_releaseVector.
!</description>
  
!<input>
  
  ! Desired length of the vector
  INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ

  ! Whether to fill the vector with zero initially
  LOGICAL, INTENT(IN)                               :: bclear

  ! OPTIONAL: Data type of the vector.
  ! If not specified, ST_DOUBLE is assumed.
  INTEGER, INTENT(IN), OPTIONAL                     :: cdataType  
  
!</input>

!<output>
  ! Scalar vector structure
  TYPE(t_vectorScalar), INTENT(OUT)                 :: rvector
!</output>

!</subroutine>

  INTEGER ::cdata

    cdata = ST_DOUBLE
    IF (PRESENT(cdataType)) cdata=cdataType

    ! The INTENT(OUT) already initialises rvector with the most important
    ! information. The rest comes now:
    !
    ! Size:
    rvector%NEQ = MAX(0,NEQ)
    
    ! Handle - if NEQ > 0
    IF (rvector%NEQ .GT. 0) THEN
      IF (bclear) THEN
        CALL storage_new1D ('lsyssc_createVector', 'ScalarVector', rvector%NEQ, &
                            cdata, rvector%h_Ddata,ST_NEWBLOCK_ZERO)
      ELSE
        CALL storage_new1D ('lsyssc_createVector', 'ScalarVector', rvector%NEQ, &
                            cdata, rvector%h_Ddata,ST_NEWBLOCK_NOINIT)
      END IF
    END IF
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_createVecByDiscr (rdiscretisation,rx,bclear,cdataType)
  
!<description>
  ! Initialises the vector structure rx based on a discretisation
  ! structure rDiscretisation. 
  !
  ! Memory is allocated on the heap for rx accordint to the number of 
  ! DOF's indicated by the spatial discretisation structures in 
  ! rdiscretisation. 
!</description>
  
!<input>
  ! A block discretisation structure specifying the spatial discretisations
  ! for all the subblocks in rx.
  TYPE(t_spatialDiscretisation),INTENT(IN), TARGET :: rdiscretisation
  
  ! Optional: If set to YES, the vector will be filled with zero initially.
  ! Otherwise the content of rx is undefined.
  LOGICAL, INTENT(IN), OPTIONAL             :: bclear
  
  ! OPTIONAL: Data type identifier for the entries in the vector. 
  ! Either ST_SINGLE or ST_DOUBLE. If not present, ST_DOUBLE is used.
  INTEGER, INTENT(IN),OPTIONAL              :: cdataType
!</input>

!<output>
  ! Destination structure. Memory is allocated appropriately.
  ! A pointer to rdiscretisation is saved to r.
  TYPE(t_vectorScalar),INTENT(OUT) :: rx
!</output>
  
!</subroutine>

  INTEGER :: cdata
  INTEGER(PREC_VECIDX) :: NEQ
  
  cdata = ST_DOUBLE
  IF (PRESENT(cdataType)) cdata = cdataType
  
  ! Get NEQ:
  NEQ = dof_igetNDofGlob(rdiscretisation)
  
  ! Create a new vector with that block structure
  CALL lsyssc_createVector (rx, NEQ, bclear, cdataType)
  
  ! Initialise further data of the block vector
  rx%p_rspatialDiscretisation => rdiscretisation
  
  END SUBROUTINE
  
  !****************************************************************************
!<subroutine>
  
  REAL(DP) FUNCTION lsyssc_scalarProduct (rx, ry)
  
!<description>
  ! Calculates a scalar product of two vectors.
!</description>
  
!<input>
  ! First vector
  TYPE(t_vectorScalar), INTENT(IN)                  :: rx

  ! Second vector
  TYPE(t_vectorScalar), INTENT(IN)                  :: ry

!</input>

!<result>
  ! The scalar product (rx,ry) of the two vectors.
!</result>

!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata1dp
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata2dp
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata1dp
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata2dp
  REAL(DP) :: res
  INTEGER(PREC_VECIDX) i
  
  ! Vectors must be compatible!
  CALL lsyssc_isVectorCompatible (rx,ry)
  
  ! Is there data at all?
  res = 0.0_DP
  
  IF ( (rx%NEQ .EQ. 0) .OR. (ry%NEQ .EQ. 0) .OR. (rx%NEQ .NE. rx%NEQ)) THEN
    PRINT *,'Error in lsyssc_scalarProduct: Vector dimensions wrong!'
    STOP
  END IF
  
  IF (rx%cdataType .NE. ry%cdataType) THEN
    PRINT *,'lsyssc_scalarProduct: Data types different!'
    STOP
  END IF
  
  ! Take care of the data type before doing a scalar product!
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
     
    ! Get the data arrays
    CALL lsyssc_getbase_double (rx,p_Ddata1dp)
    CALL lsyssc_getbase_double (ry,p_Ddata2dp)
    
    ! Perform the scalar product
    res = p_Ddata1dp(1)*p_Ddata2dp(1)
    DO i=2,rx%NEQ
      res = res + p_Ddata1dp(i)*p_Ddata2dp(i)
    END DO
    
  CASE (ST_SINGLE)

    ! Get the data arrays
    CALL lsyssc_getbase_single (rx,p_Fdata1dp)
    CALL lsyssc_getbase_single (ry,p_Fdata2dp)
    
    ! Perform the scalar product
    res = p_Fdata1dp(1)*p_Fdata2dp(1)
    DO i=2,rx%NEQ
      res = res + p_Fdata1dp(i)*p_Fdata2dp(i)
    END DO
    
  CASE DEFAULT
    PRINT *,'lsyssc_scalarProduct: Not supported precision combination'
    STOP
  END SELECT
  
  ! Return the scalar product, finish
  lsyssc_scalarProduct = res

  END FUNCTION

  !****************************************************************************

!<subroutine>
  
  SUBROUTINE lsyssc_scalarMatVec (rmatrix, rx, ry, cx, cy)
  
!<description>
  ! Performs a matrix vector multiplicationwith a given scalar matrix:
  !    $$ Dy   =   cx * rMatrix * rx   +   cy * ry $$
!</description>
  
!<input>
  
  ! Scalar matrix
  TYPE(t_matrixScalar), INTENT(IN)                  :: rmatrix

  ! Vector to multiply with the matrix.
  TYPE(t_vectorScalar), INTENT(IN)                  :: rx
  
  ! Multiplicative factor for rx
  REAL(DP), INTENT(IN)                              :: cx

  ! Multiplicative factor for ry
  REAL(DP), INTENT(IN)                              :: cy
  
!</input>

!<inputoutput>
  ! Additive vector. Receives the result of the matrix-vector multiplication
  TYPE(t_vectorScalar), INTENT(INOUT)               :: ry
!</inputoutput>

!</subroutine>
  
    ! Vectors must be compatible to the matrix.
    CALL lsyssc_isMatrixCompatible (rx,rmatrix,.FALSE.)
    CALL lsyssc_isMatrixCompatible (ry,rmatrix,.TRUE.)
    
    ! rx and ry must have at least the same data type!
    IF (rx%cdataType .NE. ry%cdataType) THEN
      PRINT *,'MV with different data types for rx and ry not supported!'
      STOP
    END IF
    
    IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .EQ. 0) THEN
      ! Select the right MV multiplication routine from the matrix format
      SELECT CASE (rmatrix%cmatrixFormat)
      
      CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX9)
      
        ! Take care of the precision of the entries
        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
          ! Format 7 and Format 9 multiplication
          SELECT CASE (rx%cdataType)
          
          CASE (ST_DOUBLE)
            ! double precision matrix, double precision vectors
            CALL lsyssc_LAX79doubledouble (rmatrix,rx,ry,cx,cy)
          
          CASE DEFAULT
            PRINT *,'Only double precision vectors supported for now in MV!'
            STOP
            
          END SELECT
          
        CASE DEFAULT
          PRINT *,'Only double precision matrices supported for now in MV!'
          STOP
        END SELECT
        
      CASE DEFAULT
        PRINT *,'Unknown matrix format in MV-multiplication!'
        STOP
      END SELECT
      
    ELSE
      ! Transposed matrix
      ! Select the right MV multiplication routine from the matrix format
      SELECT CASE (rmatrix%cmatrixFormat)
      
      CASE (LSYSSC_MATRIX9)
      
        ! Take care of the precision of the entries
        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
          ! Format 9 multiplication
          SELECT CASE (rx%cdataType)
          
          CASE (ST_DOUBLE)
            ! double precision matrix, double precision vectors
            CALL lsyssc_LTX9doubledouble (rmatrix,rx,ry,cx,cy)
          
          CASE DEFAULT
            PRINT *,'Only double precision vectors supported for now in MV!'
            STOP
            
          END SELECT
          
        CASE DEFAULT
          PRINT *,'Only double precision matrices supported for now in MV!'
          STOP
        END SELECT
        
      CASE (LSYSSC_MATRIX7)
      
        ! Take care of the precision of the entries
        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
          ! Format 7 multiplication
          SELECT CASE (rx%cdataType)
          
          CASE (ST_DOUBLE)
            ! double precision matrix, double precision vectors
            CALL lsyssc_LTX7doubledouble (rmatrix,rx,ry,cx,cy)
          
          CASE DEFAULT
            PRINT *,'Only double precision vectors supported for now in MV!'
            STOP
            
          END SELECT
          
        CASE DEFAULT
          PRINT *,'Only double precision matrices supported for now in MV!'
          STOP
        END SELECT
        
      CASE DEFAULT
        PRINT *,'Unknown matrix format in MV-multiplication!'
        STOP
      END SELECT
      
    END IF
  
  CONTAINS
  
    ! Now the real MV multiplication routines follow.
    ! We create them in the scoping unit of the procedure to prevent
    ! direct calls from outside.
    
    !**************************************************************
    ! Format 7 and Format 9 multiplication
    ! double precision matrix,
    ! double precision vectors
    
    SUBROUTINE lsyssc_LAX79doubledouble (rmatrix,rx,ry,cx,cy)

    ! Save arguments as above - given as parameters as some compilers
    ! might have problems with scoping units...
    TYPE(t_matrixScalar), INTENT(IN)                  :: rmatrix
    TYPE(t_vectorScalar), INTENT(IN)                  :: rx
    REAL(DP), INTENT(IN)                              :: cx
    REAL(DP), INTENT(IN)                              :: cy
    TYPE(t_vectorScalar), INTENT(INOUT)               :: ry

    REAL(DP), DIMENSION(:), POINTER :: p_DA, p_Dx, p_Dy
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_VECIDX) :: irow,icol
    REAL(DP) :: dtmp
    INTEGER(PREC_VECIDX) :: NEQ

      ! Get the matrix
      CALL storage_getbase_double (rmatrix%h_DA,p_DA)
      CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
      CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)
      
      ! Get NEQ - from the matrix, not from the vector!
      NEQ = rmatrix%NEQ

      ! Get the vectors
      CALL lsyssc_getbase_double (rx,p_Dx)
      CALL lsyssc_getbase_double (ry,p_Dy)
      
      ! Perform the multiplication
      IF (cx .NE. 0.0_DP) THEN
      
        IF (cy .EQ. 0.0_DP) THEN
        
          ! cy = 0. We have simply to make matrix*vector without adding ry.
          ! Multiply the first entry in each line of the matrix with the
          ! corresponding entry in rx and add it to ry.
          ! Don't multiply with cy, this comes later.
          !
          ! What is this complicated IF-THEN structure for?
          ! Well, to prevent an initialisation of rx with zero in case cy=0!
          
          DO irow=1,NEQ
            icol = p_Kcol(p_Kld(irow))
            p_Dy(irow) = p_Dx(icol) * p_DA(p_Kld(irow))
          END DO
          
          ! Now we have an initial ry where we can do a usual MV
          ! with the rest of the matrix...
          
        ELSE 
        
          ! cy <> 0. We have to perform matrix*vector + vector.
          ! What we actually calculate here is:
          !    ry  =  cx * A * x  +  cy * y
          !        =  cx * ( A * x  +  cy/cx * y).
          !
          ! Scale down y:
        
          dtmp = cy/cx
          IF (dtmp .NE. 1.0_DP) THEN
            CALL lalg_scaleVectorDble(p_Dy,dtmp)
          END IF
          
          ! Multiply the first entry in each line of the matrix with the
          ! corresponding entry in rx and add it to the (scaled) ry.
          
          DO irow=1,NEQ
            ICOL = p_Kcol(p_Kld(irow))
            p_Dy(irow) = p_Dx(icol)*p_DA(p_Kld(irow)) + p_Dy(irow) 
          END DO
          
        ENDIF
        
        ! Multiply the rest of rx with the matrix and add it to ry:
        
        DO irow=1,NEQ
          DO icol = p_Kld(irow)+1,p_Kld(irow+1)-1
            p_Dy(irow) = p_Dy(irow) + p_DA(icol)*p_Dx(p_Kcol(icol))
          END DO
        END DO
        
        ! Scale by cx, finish.
        
        IF (cx .NE. 1.0_DP) THEN
          CALL lalg_scaleVectorDble (p_Dy,cx)
        END IF
        
      ELSE 
        ! cx = 0. The formula is just a scaling of the vector ry!
        CALL lalg_scaleVectorDble(p_Dy,cy)
      ENDIF
      
    END SUBROUTINE
   
    !**************************************************************
    ! Format 7 multiplication, transposed matrix
    ! double precision matrix,
    ! double precision vectors
    
    SUBROUTINE lsyssc_LTX7doubledouble (rmatrix,rx,ry,cx,cy)

    ! Save arguments as above - given as parameters as some compilers
    ! might have problems with scoping units...
    TYPE(t_matrixScalar), INTENT(IN)                  :: rmatrix
    TYPE(t_vectorScalar), INTENT(IN)                  :: rx
    REAL(DP), INTENT(IN)                              :: cx
    REAL(DP), INTENT(IN)                              :: cy
    TYPE(t_vectorScalar), INTENT(INOUT)               :: ry

    REAL(DP), DIMENSION(:), POINTER :: p_DA, p_Dx, p_Dy
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_VECIDX) :: irow,icol
    REAL(DP) :: dtmp
    INTEGER(PREC_VECIDX) :: NEQ

      ! Get the matrix
      CALL storage_getbase_double (rmatrix%h_DA,p_DA)
      CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
      CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)
      
      ! NCOLS(real matrix) = NEQ(saved matrix structure) !
      NEQ = rmatrix%NCOLS

      ! Get the vectors
      CALL lsyssc_getbase_double (rx,p_Dx)
      CALL lsyssc_getbase_double (ry,p_Dy)
      
      ! Perform the multiplication.
      IF (cx .NE. 0.0_DP) THEN
      
        IF (cy .EQ. 0.0_DP) THEN
        
          ! cy = 0. We have simply to make matrix*vector without adding ry.
          ! Multiply the first entry in each line of the matrix with the
          ! corresponding entry in rx and add it to ry.
          ! Don't multiply with cy, this comes later.
          !
          ! What is this complicated IF-THEN structure for?
          ! Well, to prevent an initialisation of rx with zero in case cy=0!
          
          DO irow=1,NEQ
            p_Dy(irow) = p_Dx(irow)*p_DA(p_Kld(irow)) 
          END DO
          
          ! Now we have an initial ry where we can do a usual MV
          ! with the rest of the matrix...
          
        ELSE 
        
          ! cy <> 0. We have to perform matrix*vector + vector.
          ! What we actually calculate here is:
          !    ry  =  cx * A^t * x  +  cy * y
          !        =  cx * ( A^t * x  +  cy/cx * y).
          !        =  cx * ( (x^t * A)  +  cy/cx * y^t)^t.
          !
          ! Scale down y:
        
          dtmp = cy/cx
          IF (dtmp .NE. 1.0_DP) THEN
            CALL lalg_scaleVectorDble(p_Dy,dtmp)
          END IF
          
          ! Multiply the first entry in each line of the matrix with the
          ! corresponding entry in rx and add it to the (scaled) ry.
          
          DO irow=1,NEQ
            p_Dy(irow) = p_Dy(irow) + p_Dx(irow)*p_DA(p_Kld(irow)) 
          END DO
          
        ENDIF
        
        ! Multiply the rest of rx with the matrix and add it to ry:
        
        DO irow=1,NEQ
          DO icol = p_Kld(irow)+1,p_Kld(irow+1)-1
            p_Dy(p_Kcol(icol)) = p_Dy(p_Kcol(icol)) + p_Dx(irow)*p_DA(icol)
          END DO
        END DO
        
        ! Scale by cx, finish.
        
        IF (cx .NE. 1.0_DP) THEN
          CALL lalg_scaleVectorDble (p_Dy,cx)
        END IF
        
      ELSE 
        ! cx = 0. The formula is just a scaling of the vector ry!
        CALL lalg_scaleVectorDble(p_Dy,cy)
      ENDIF
      
    END SUBROUTINE

    !**************************************************************
    ! Format 9 multiplication, transposed matrix
    ! double precision matrix,
    ! double precision vectors
    
    SUBROUTINE lsyssc_LTX9doubledouble (rmatrix,rx,ry,cx,cy)

    ! Save arguments as above - given as parameters as some compilers
    ! might have problems with scoping units...
    TYPE(t_matrixScalar), INTENT(IN)                  :: rmatrix
    TYPE(t_vectorScalar), INTENT(IN)                  :: rx
    REAL(DP), INTENT(IN)                              :: cx
    REAL(DP), INTENT(IN)                              :: cy
    TYPE(t_vectorScalar), INTENT(INOUT)               :: ry

    REAL(DP), DIMENSION(:), POINTER :: p_DA, p_Dx, p_Dy
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_VECIDX) :: irow,icol
    REAL(DP) :: dtmp
    INTEGER(PREC_VECIDX) :: NEQ

      ! Get the matrix
      CALL storage_getbase_double (rmatrix%h_DA,p_DA)
      CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
      CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)
      
      ! NCOLS(real matrix) = NEQ(saved matrix structure) !
      NEQ = rmatrix%NCOLS

      ! Get the vectors
      CALL lsyssc_getbase_double (rx,p_Dx)
      CALL lsyssc_getbase_double (ry,p_Dy)
      
      ! Perform the multiplication.
      IF (cx .NE. 0.0_DP) THEN
      
        IF (cy .EQ. 0.0_DP) THEN
        
          ! cy = 0. Clear the output vector at first
          CALL lalg_clearVectorDble (p_Dy)
          
          ! Now we have an empty ry where we can do a usual MV
          ! with the rest of the matrix...
          
        ELSE 
        
          ! cy <> 0. We have to perform matrix*vector + vector.
          ! What we actually calculate here is:
          !    ry  =  cx * A^t * x  +  cy * y
          !        =  cx * ( A^t * x  +  cy/cx * y).
          !        =  cx * ( (x^t * A)  +  cy/cx * y^t)^t.
          !
          ! Scale down y:
        
          dtmp = cy/cx
          IF (dtmp .NE. 1.0_DP) THEN
            CALL lalg_scaleVectorDble(p_Dy,dtmp)
          END IF
          
        ENDIF
        
        ! Multiply rx with the matrix and add it to ry:
        
        DO irow=1,NEQ
          DO icol = p_Kld(irow),p_Kld(irow+1)-1
            p_Dy(p_Kcol(icol)) = p_Dy(p_Kcol(icol)) + p_Dx(irow)*p_DA(icol)
          END DO
        END DO
        
        ! Scale by cx, finish.
        
        IF (cx .NE. 1.0_DP) THEN
          CALL lalg_scaleVectorDble (p_Dy,cx)
        END IF
        
      ELSE 
        ! cx = 0. The formula is just a scaling of the vector ry!
        CALL lalg_scaleVectorDble(p_Dy,cy)
      ENDIF
      
    END SUBROUTINE

  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_duplicateVector (roldVector,rnewVector)
  
!<description>
  ! Duplicates an existing vector. All structural data from roldVector
  ! is assigned to rnewVector. A new handle for vector data will be allocated
  ! in rnewVector and the data of the vector in roldVector is copied
  ! to the new data array.
!</description>
  
!<input>
  ! Vector to copy
  TYPE(t_vectorScalar), INTENT(IN)                :: roldVector
!</input>

!<output>
  ! The new vector which will be a copy of roldvector
  TYPE(t_vectorScalar), INTENT(OUT)               :: rnewVector
!</output>

!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource,p_Ddest
  REAL(SP), DIMENSION(:), POINTER :: p_Fsource,p_Fdest
  INTEGER(PREC_VECIDX) :: NEQ

  ! At first, copy all 'local' data.
  rnewVector = roldVector
  
  ! Then allocate a new array for the content in the same data
  ! format as the vector and copy the data.
  ! We can't use storage_copy here, as the content of roldVector
  ! maybe > NEQ!
  
  NEQ = roldVector%NEQ
  CALL storage_new ('lsyssc_duplicateVector','vec-copy',NEQ,&
                    roldVector%cdataType, rnewVector%h_Ddata, &
                    ST_NEWBLOCK_NOINIT)

  ! The new vector starts at index position 1
  rnewVector%iidxFirstEntry = 1

  ! Take care of the index of the first entry when copying the data!

  SELECT CASE (roldVector%cdataType)
  CASE (ST_DOUBLE)
    CALL lsyssc_getbase_double (roldVector,p_Dsource)
    CALL lsyssc_getbase_double (rnewVector,p_Ddest)
    CALL lalg_copyVectorDble (p_Dsource,p_Ddest)
    
  CASE (ST_SINGLE)
    CALL lsyssc_getbase_single (roldVector,p_Fsource)
    CALL lsyssc_getbase_single (rnewVector,p_Fdest)
    CALL lalg_copyVectorSngl (p_Fsource,p_Fdest)

  CASE DEFAULT
    PRINT *,'lsyssc_duplicateVector: Unsupported data type!'
    STOP
  END SELECT
   
  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_duplicateMatrix (rsourceMatrix,rdestMatrix,&
                                                idupStructure, idupContent)
  
!<description>
  ! Duplicates an existing matrix, creates a new matrix rdestMatrix based
  ! on a template matrix rsourceMatrix.
  ! Duplicating a matrix does not necessarily mean that new memory is
  ! allocated and the matrix entries are copied to that. The two flags
  ! iduipStructure and idupContent decide on how to set up rdestMatrix.
  ! Depending on their setting, it's possible to copy only the handles
  ! of such dynamic information, so that both matrices share the same
  ! information.
  !
  ! We follow the following convention:
  !  Structure = Column structure, Row structure, Sorting permutation,
  !              Discretisation-related information
  !  Content   = Enties in a matrix.
  !
  ! Remark: There is never memory allocated on the heap for the sorting
  !  permutation. A matrix is never the 'owner' of a permutation, i.e.
  !  does not maintain it. Therefore, copying a permutation means
  !  copying the corresponding handle. The application must keep track
  !  of the permutations.
!</description>
  
!<input>
  ! Source matrix.
  TYPE(t_matrixScalar), INTENT(IN)               :: rsourceMatrix
  
  ! Duplication flag that decides on how to set up the structure
  ! of rdestMatrix. One of the LSYSSC_DUP_xxxx flags:
  ! LSYSSC_DUP_IGNORE     : Don't set up the structure of rdestMatrix. Any
  !   matrix structure is ignored and therefore preserved.
  ! LSYSSC_DUP_REMOVE     : Removes any existing matrix structure from 
  !   rdestMatrix if there is any. Releases memory if necessary.
  !   Does not delete 'static' information like NEQ,NCOLS,NA,...
  ! LSYSSC_DUP_DISMISS    : Removes any existing matrix structure from 
  !   rdestMatrix if there is any. No memory is released, handles are simply
  !   dismissed. Does not delete 'static' information like NEQ,NCOLS,NA,...
  ! LSYSSC_DUP_SHARE      : rdestMatrix receives the same handles for
  !   structural data as rsourceMatrix and therefore shares the same structure.
  ! LSYSSC_DUP_COPY       : rdestMatrix gets a copy of the structure of 
  !   rsourceMatrix. If necessary, new memory is allocated for the structure. 
  !   If rdestMatrix already contains allocated memory, structural data
  !   is simply copied from rsourceMatrix into that.
  ! LSYSSC_DUP_ASIS       : Duplicate by ownership. If the structure of 
  !   rsourceMatrix belongs to rsourceMatrix, rdestMatrix gets a copy
  !   of the structure; new memory is allocated if necessary (the same as 
  !   LSYSSC_DUP_COPY). If the structure of rsourceMatrix belongs to another
  !   matrix than rsourceMatrix, rdestMatrix receives the same handles as
  !   rsourceMatrix and is therefore a third matrix sharing the same structure
  !   (the same as LSYSSC_DUP_SHARE, so rsourceMatrix, rdestMatrix and the 
  !   other matrix have the same structure).
  ! LSYSSC_DUP_EMPTY      : New memory is allocated for the structure in the
  !   same size as the structure in rsourceMatrix but no data is copied;
  !   the arrays are left uninitialised.
  ! LSYSSC_DUP_TEMPLATE   : Copies statis structural information about the
  !   structure (NEQ, NCOLS,...) to the destination matrix. Dynamic information
  !   is removed from the destination matrix, all handles are reset.
  INTEGER, INTENT(IN)                            :: idupStructure
  
  ! Duplication flag that decides on how to set up the content
  ! of rdestMatrix. One of the LSYSSC_DUP_xxxx flags:
  ! LSYSSC_DUP_IGNORE     : Don't set up the content of rdestMatrix. Any
  !   matrix content is ignored and therefore preserved.
  ! LSYSSC_DUP_REMOVE     : Removes any existing matrix content from 
  !   rdestMatrix if there is any. Releases memory if necessary.
  ! LSYSSC_DUP_DISMISS    : Removes any existing matrix content from 
  !   rdestMatrix if there is any. No memory is released, handles are simply
  !   dismissed.
  ! LSYSSC_DUP_SHARE      : rdestMatrix receives the same handles for
  !   matrix content data as rsourceMatrix and therefore shares the same content.
  ! LSYSSC_DUP_COPY       : rdestMatrix gets a copy of the content of rsourceMatrix.
  !   If necessary, new memory is allocated for the content.
  !   If rdestMatrix already contains allocated memory, content data
  !   is simply copied from rsourceMatrix into that.
  ! LSYSSC_DUP_ASIS       : Duplicate by ownership. If the content of 
  !   rsourceMatrix belongs to rsourceMatrix, rdestMatrix gets a copy
  !   of the content; new memory is allocated if necessary (the same as 
  !   LSYSSC_DUP_COPY). If the content of rsourceMatrix belongs to another
  !   matrix than rsourceMatrix, rdestMatrix receives the same handles as
  !   rsourceMatrix and is therefore a third matrix sharing the same content
  !   (the same as LSYSSC_DUP_SHARE, so rsourceMatrix, rdestMatrix and the 
  !   other matrix have the same content).
  ! LSYSSC_DUP_EMPTY      : New memory is allocated for the content in the
  !   same size as the structure in rsourceMatrix but no data is copied;
  !   the arrays are left uninitialised.
  ! LSYSSC_DUP_TEMPLATE   : Copies statis structural information about the
  !   structure (NEQ, NCOLS,...) to the destination matrix. Dynamic information
  !   is removed from the destination matrix, all handles are reset.
  INTEGER, INTENT(IN)                            :: idupContent
!</input>

!<output>
  ! Destination matrix.
  TYPE(t_matrixScalar), INTENT(INOUT)            :: rdestMatrix
!</output>  

!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX) :: isize
  
    ! What do we have to do for the structure?
    SELECT CASE (idupStructure)
    CASE (LSYSSC_DUP_IGNORE)
      ! Nothing
    CASE (LSYSSC_DUP_REMOVE)
      ! Remove the structure - if there is any.
      CALL removeStructure(rdestMatrix, .TRUE.)
      
    CASE (LSYSSC_DUP_DISMISS)
      ! Dismiss the structure - if there is any.
      CALL removeStructure(rdestMatrix, .FALSE.)
      
    CASE (LSYSSC_DUP_TEMPLATE)
      ! Remove the structure - if there is any.
      CALL removeStructure(rdestMatrix, .TRUE.)
      
      ! Copy static structural information
      CALL copyStaticStructure(rsourceMatrix,rdestMatrix)

    CASE (LSYSSC_DUP_SHARE)
      ! Share the structure between rsourceMatrix and rdestMatrix
      CALL shareStructure(rsourceMatrix, rdestMatrix)   
      
    CASE (LSYSSC_DUP_COPY)
      ! Copy the structure of rsourceMatrix to rdestMatrix
      CALL copyStructure(rsourceMatrix, rdestMatrix)  
      
    CASE (LSYSSC_DUP_ASIS)
      ! What's with the source matrix. Does the structure belong to it?
      IF (IAND(rsourceMatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .NE. 0) THEN
        ! Copy the structure
        CALL copyStructure(rsourceMatrix, rdestMatrix)
      ELSE
        ! Share the structure
        CALL shareStructure(rsourceMatrix, rdestMatrix)  
      END IF
      
    CASE (LSYSSC_DUP_EMPTY)
      ! Start with sharing the structure
      CALL shareStructure(rsourceMatrix, rdestMatrix)
      ! Reset ownership
      rdestMatrix%imatrixSpec = IAND(rdestMatrix%imatrixSpec,&
                                     NOT(LSYSSC_MSPEC_STRUCTUREISCOPY))

      ! And at last recreate the arrays.
      ! Which source matrix do we have?  
      SELECT CASE (rsourceMatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9)
        CALL storage_new ('lsyssc_duplicateMatrix', 'KCOL', &
            rsourceMatrix%NA, ST_INT, &
            rdestMatrix%h_Kcol, ST_NEWBLOCK_NOINIT)

        CALL storage_new ('lsyssc_duplicateMatrix', 'KLD', &
            rdestMatrix%NEQ+1, ST_INT, &
            rdestMatrix%h_Kld, ST_NEWBLOCK_NOINIT)

        CALL storage_new ('lsyssc_duplicateMatrix', 'Kdiagonal', &
            rdestMatrix%NEQ, ST_INT, &
            rdestMatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT)
        
      CASE (LSYSSC_MATRIX7)
        CALL storage_new ('lsyssc_duplicateMatrix', 'KCOL', &
            rdestMatrix%NA, ST_INT, &
            rdestMatrix%h_Kcol, ST_NEWBLOCK_NOINIT)

        CALL storage_new ('lsyssc_duplicateMatrix', 'KLD', &
            rdestMatrix%NEQ+1, ST_INT, &
            rdestMatrix%h_Kld, ST_NEWBLOCK_NOINIT)

      END SELECT
    END SELECT
    
    ! -----
    ! Ok, handling of the structure is finished. The data is handled similar.
    ! -----
    
    ! What do we have to do for the content?
    SELECT CASE (idupContent)
    CASE (LSYSSC_DUP_IGNORE)
      ! Nothing
    CASE (LSYSSC_DUP_REMOVE)
      ! Remove the structure - if there is any.
      CALL removeContent(rdestMatrix, .TRUE.)
      
    CASE (LSYSSC_DUP_DISMISS)
      ! Dismiss the structure - if there is any.
      CALL removeContent(rdestMatrix, .FALSE.)
      
    CASE (LSYSSC_DUP_SHARE)
      ! Share the structure between rsourceMatrix and rdestMatrix
      CALL shareContent(rsourceMatrix, rdestMatrix)   
      
    CASE (LSYSSC_DUP_COPY)
      ! Copy the structure of rsourceMatrix to rdestMatrix
      CALL copyContent(rsourceMatrix, rdestMatrix)  
      
    CASE (LSYSSC_DUP_TEMPLATE)
      ! Remove the structure - if there is any.
      CALL removeContent(rdestMatrix, .TRUE.)
      
      ! Copy static content information
      CALL copyStaticContent(rsourceMatrix,rdestMatrix)

    CASE (LSYSSC_DUP_ASIS)
      ! What's with the source matrix. Does the structure belong to it?
      IF (IAND(rsourceMatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .NE. 0) THEN
        ! Copy the structure
        CALL copyContent(rsourceMatrix, rdestMatrix)
      ELSE
        ! Share the structure
        CALL shareContent(rsourceMatrix, rdestMatrix)  
      END IF
      
    CASE (LSYSSC_DUP_EMPTY)
      ! Start with sharing the structure
      CALL shareContent(rsourceMatrix, rdestMatrix)
      ! Reset ownership
      rdestMatrix%imatrixSpec = IAND(rdestMatrix%imatrixSpec,&
                                     NOT(LSYSSC_MSPEC_CONTENTISCOPY))

      ! And at last recreate the arrays.
      ! Which source matrix do we have?  
      SELECT CASE (rsourceMatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
        ! Create a new content array in the same data type as the original matrix
        CALL storage_new('lsyssc_duplicateMatrix', 'DA', &
            rdestMatrix%NA, rsourceMatrix%cdataType, &
            rdestMatrix%h_DA, ST_NEWBLOCK_NOINIT)
      END SELECT
    END SELECT
    
    ! -----
    ! Final check. Check if we destroyed the matrix. May only happen if the
    ! user on-purpose calls this routine two times with different source
    ! but the same destination matrix.
    ! -----
    
    SELECT CASE (rdestMatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9)
    
      ! Check length of DA
      IF (rdestMatrix%h_DA .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_DA,isize)
        IF (isize .NE. rdestMatrix%NA) THEN
          PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA != length(DA)!'
          STOP
        END IF
      END IF
      
      ! Check length of KCOL
      IF (rdestMatrix%h_Kcol .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_Kcol,isize)
        IF (isize .NE. rdestMatrix%NA) THEN
          PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA != length(KCOL)!'
          STOP
        END IF
      END IF

      ! Check length of KLD
      IF (rdestMatrix%h_Kld .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_Kld,isize)
        
        ! Be careful, matrix may be transposed.
        IF (IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .EQ. 0) THEN
          IF (isize .NE. rdestMatrix%NEQ+1) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NEQ+1 != length(KLD)!'
            STOP
          END IF
        ELSE
          IF (isize .NE. rdestMatrix%NCOLS+1) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NEQ+1 != length(KLD)!'
            STOP
          END IF
        END IF
      END IF
      
      ! Check length of Kdiagonal
      IF (rdestMatrix%h_Kdiagonal .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_Kdiagonal,isize)
        
        ! Be careful, matrix may be transposed.
        IF (IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .EQ. 0) THEN
          IF (isize .NE. rdestMatrix%NEQ) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NEQ+1 != length(Kdiag)!'
            STOP
          END IF
        ELSE
          IF (isize .NE. rdestMatrix%NCOLS) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NEQ+1 != length(Kdiag)!'
            STOP
          END IF
        END IF
      END IF
      
    CASE (LSYSSC_MATRIX7)
    
      ! Check length of DA
      IF (rdestMatrix%h_DA .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_DA,isize)
        IF (isize .NE. rdestMatrix%NA) THEN
          PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA != length(DA)!'
          STOP
        END IF
      END IF
      
      ! Check length of KCOL
      IF (rdestMatrix%h_Kcol .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_Kcol,isize)
        IF (isize .NE. rdestMatrix%NA) THEN
          PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA != length(KCOL)!'
          STOP
        END IF
      END IF

      ! Check length of KLD
      IF (rdestMatrix%h_Kld .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_Kld,isize)
        
        ! Be careful, matrix may be transposed.
        IF (IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .EQ. 0) THEN
          IF (isize .NE. rdestMatrix%NEQ) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NEQ+1 != length(KLD)!'
            STOP
          END IF
        ELSE
          IF (isize .NE. rdestMatrix%NCOLS) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NEQ+1 != length(KLD)!'
            STOP
          END IF
        END IF
      END IF
      
    END SELECT

  CONTAINS
    
    !--------------------------------------------------------
    ! Auxiliary routine: Remove the structure from a matrix.
    
    SUBROUTINE removeStructure(rmatrix, brelease)
    
    ! The matrix to to be processed.
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
    
    ! Whether to release data from the heap (if it belongs to the matrix)
    ! or simply to overwrite the handles with ST_NOHANDLE
    LOGICAL, INTENT(IN) :: brelease
    
      ! This is a matrix-dependent task
      SELECT CASE (rmatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9)
        ! Release the handles from the heap?
        ! Only release it if the data belongs to this matrix.
        IF (brelease .AND. &
            (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .EQ. 0)) THEN
          IF (rmatrix%h_Kcol .NE. ST_NOHANDLE) CALL storage_free (rmatrix%h_Kcol)
          IF (rmatrix%h_Kld .NE. ST_NOHANDLE) CALL storage_free (rmatrix%h_Kld)
          IF (rmatrix%h_Kdiagonal .NE. ST_NOHANDLE) CALL storage_free (rmatrix%h_Kdiagonal)
        END IF
        
        ! Reset the handles
        rmatrix%h_Kld = ST_NOHANDLE
        rmatrix%h_Kcol = ST_NOHANDLE
        rmatrix%h_Kdiagonal = ST_NOHANDLE
        
        ! Reset the ownership-status
        rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,&
                                   NOT(LSYSSC_MSPEC_STRUCTUREISCOPY))
      CASE (LSYSSC_MATRIX7)
        ! Release the handles from the heap?
        ! Only release it if the data belongs to this matrix.
        IF (brelease .AND. &
            (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .EQ. 0)) THEN
          IF (rmatrix%h_Kcol .NE. ST_NOHANDLE) CALL storage_free (rmatrix%h_Kcol)
          IF (rmatrix%h_Kld .NE. ST_NOHANDLE) CALL storage_free (rmatrix%h_Kld)
        END IF
        
        ! Reset the handles
        rmatrix%h_Kld = ST_NOHANDLE
        rmatrix%h_Kcol = ST_NOHANDLE
        
      END SELECT

      ! Reset the ownership-status
      rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,&
                                  NOT(LSYSSC_MSPEC_STRUCTUREISCOPY))
    
    END SUBROUTINE

    !--------------------------------------------------------
    ! Auxiliary routine: Shares the structure between rsourceMatrix
    ! and rdestMatrix
    
    SUBROUTINE shareStructure(rsourceMatrix, rdestMatrix)
    
    ! The source matrix 
    TYPE(t_matrixScalar), INTENT(IN) :: rsourceMatrix

    ! The destination matrix 
    TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
    
      ! local variables
      INTEGER(I32) :: iflag,iflag2
    
      ! Overwrite structural data
      rdestMatrix%NA     = rsourceMatrix%NA
      rdestMatrix%NEQ    = rsourceMatrix%NEQ
      rdestMatrix%NCOLS  = rsourceMatrix%NCOLS
      rdestMatrix%cmatrixFormat       = rsourceMatrix%cmatrixFormat
      rdestMatrix%isortStrategy       = rsourceMatrix%isortStrategy
      rdestMatrix%h_IsortPermutation  = rsourceMatrix%h_IsortPermutation
      
      ! Transfer all flags except the 'dup' flags
      iflag = IAND(rsourceMatrix%imatrixSpec,NOT(LSYSSC_MSPEC_ISCOPY))
      iflag2 = IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
      rdestMatrix%imatrixSpec = IOR(iflag,iflag2)

      ! Transfer discretisation-related information,
      rdestMatrix%p_rspatialDiscretisation => rsourceMatrix%p_rspatialDiscretisation

      ! Which source matrix do we have?  
      SELECT CASE (rsourceMatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9)
        rdestMatrix%h_Kcol = rsourceMatrix%h_Kcol
        rdestMatrix%h_Kld  = rsourceMatrix%h_Kld
        rdestMatrix%h_Kdiagonal = rsourceMatrix%h_Kdiagonal
        
      CASE (LSYSSC_MATRIX7)
        ! Overwrite structural data
        rdestMatrix%h_Kcol = rsourceMatrix%h_Kcol
        rdestMatrix%h_Kld  = rsourceMatrix%h_Kld
      
      END SELECT
      
      ! Indicate via the matrixSpec-flag that we are not
      ! the owner of the structure.
      rdestMatrix%imatrixSpec = IOR(rdestMatrix%imatrixSpec,&
                                    LSYSSC_MSPEC_STRUCTUREISCOPY)
    
    END SUBROUTINE

    !--------------------------------------------------------
    ! Auxiliary routine: Copy static structural information
    ! (NEQ, NCOLS, NA,...) from rsourceMatrix to rdestMatrix
    
    SUBROUTINE copyStaticStructure(rsourceMatrix, rdestMatrix)
    
    ! The source matrix 
    TYPE(t_matrixScalar), INTENT(IN) :: rsourceMatrix

    ! The destination matrix 
    TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
    
      ! local variables
      INTEGER(I32) :: iflag,iflag2
    
      ! Overwrite structural data
      rdestMatrix%NA     = rsourceMatrix%NA
      rdestMatrix%NEQ    = rsourceMatrix%NEQ
      rdestMatrix%NCOLS  = rsourceMatrix%NCOLS
      rdestMatrix%cmatrixFormat       = rsourceMatrix%cmatrixFormat
      rdestMatrix%isortStrategy       = rsourceMatrix%isortStrategy
      rdestMatrix%h_IsortPermutation  = rsourceMatrix%h_IsortPermutation

      ! Transfer all flags except the 'dup' flags
      iflag = IAND(rsourceMatrix%imatrixSpec,NOT(LSYSSC_MSPEC_ISCOPY))
      iflag2 = IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
      rdestMatrix%imatrixSpec = IOR(iflag,iflag2)
      
      ! Transfer discretisation-related information,
      rdestMatrix%p_rspatialDiscretisation => rsourceMatrix%p_rspatialDiscretisation

    END SUBROUTINE

    !--------------------------------------------------------
    ! Auxiliary routine: Copy the structure of rsourceMatrix
    ! to rdestMatrix
    
    SUBROUTINE copyStructure(rsourceMatrix, rdestMatrix)
    
    ! The source matrix 
    TYPE(t_matrixScalar), INTENT(IN) :: rsourceMatrix

    ! The destination matrix 
    TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
    
      ! local variables
      INTEGER(I32) :: iflag,iflag2
      INTEGER(PREC_MATIDX) :: isize
      INTEGER(PREC_VECIDX) :: NEQ
      LOGICAL :: bremove
    
      ! Overwrite structural data
      CALL copyStaticStructure(rsourceMatrix, rdestMatrix)

      ! Check the structure if it exists and if it has the right
      ! size - then we can overwrite!
      bRemove = .FALSE.
      
      ! But at first, check if rdestMatrix is the owner of the matrix
      ! structure:
      IF (IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .NE. 0) THEN
        ! No, the structure belongs to someone else, but we want to have
        ! our own. Detach the foreign matrix structure.
        bremove = .TRUE.
        
      ELSE
        
        ! Ok, rdestMatrix owns its own matrix structure or does not have any.
        ! Check the structure if it's large enough so we can overwrite
        ! it. If the size of the arrays is not equal to those
        ! in the source matrix, release the structure and allocate a new one.
        
        IF (IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .NE. 0) THEN
          NEQ = rdestMatrix%NCOLS
        ELSE
          NEQ = rdestMatrix%NEQ
        END IF
        
        ! Which source matrix do we have?  
        SELECT CASE (rsourceMatrix%cmatrixFormat)
        CASE (LSYSSC_MATRIX9)
          IF ((rdestMatrix%h_Kcol .NE. ST_NOHANDLE) .AND. &
              (rdestMatrix%h_Kld .NE. ST_NOHANDLE) .AND. &
              (rdestMatrix%h_Kdiagonal .NE. ST_NOHANDLE)) THEN
        
            CALL storage_getsize (rdestMatrix%h_Kcol,isize)
            bremove = bremove .OR. (isize .NE. rdestMatrix%NA)
            
            CALL storage_getsize (rsourceMatrix%h_Kld,isize)
            bremove = bremove .OR. (isize .NE. NEQ+1)
            
            CALL storage_getsize (rsourceMatrix%h_Kdiagonal,isize)
            bremove = bremove .OR. (isize .NE. NEQ)
          
          ELSE

            ! Remove any partial information if there is any.
            bremove = .TRUE.
            
          END IF
          
        CASE (LSYSSC_MATRIX7)
        
          IF ((rdestMatrix%h_Kcol .NE. ST_NOHANDLE) .AND. &
              (rdestMatrix%h_Kld .NE. ST_NOHANDLE)) THEN
 
            CALL storage_getsize (rdestMatrix%h_Kcol,isize)
            bremove = bremove .OR. (isize .NE. rdestMatrix%NA)
            
            CALL storage_getsize (rsourceMatrix%h_Kld,isize)
            bremove = bremove .OR. (isize .NE. NEQ+1)

          ELSE

            ! Remove any partial information if there is any.
            bremove = .TRUE.
            
          END IF
        
        END SELECT
      
      END IF
    
      ! Remove the old matrix structure if we should do so.
      IF (bremove) CALL removeStructure(rdestMatrix, .TRUE.)

      ! Duplicate structural data from the source matrix.
      ! Storage_copy allocates new memory if necessary, as the handles are all
      ! set to ST_NOHANDLE with the above removeStructure!
      ! If the handles exist, they have the correct size, so we can overwrite
      ! the entries.
      !
      ! Which source matrix do we have?  
      SELECT CASE (rsourceMatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9)
        CALL storage_copy (rsourceMatrix%h_Kcol,rdestMatrix%h_Kcol)
        CALL storage_copy (rsourceMatrix%h_Kld,rdestMatrix%h_Kld)
        CALL storage_copy (rsourceMatrix%h_Kdiagonal,rdestMatrix%h_Kdiagonal)
        
      CASE (LSYSSC_MATRIX7)
        CALL storage_copy (rsourceMatrix%h_Kcol,rdestMatrix%h_Kcol)
        CALL storage_copy (rsourceMatrix%h_Kld,rdestMatrix%h_Kld)
      
      END SELECT
      
      ! Indicate via the matrixSpec-flag that we are the owner of the structure.
      rdestMatrix%imatrixSpec = IAND(rdestMatrix%imatrixSpec,&
                                    NOT(LSYSSC_MSPEC_STRUCTUREISCOPY))
    
    END SUBROUTINE
  
    !--------------------------------------------------------
    ! Auxiliary routine: Remove the content from a matrix.
    
    SUBROUTINE removeContent(rmatrix, brelease)
    
    ! The matrix to to be processed.
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
    
    ! Whether to release data from the heap (if it belongs to the matrix)
    ! or simply to overwrite the handles with ST_NOHANDLE
    LOGICAL, INTENT(IN) :: brelease
    
      ! This is a matrix-dependent task
      SELECT CASE (rmatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
        ! Release the handles from the heap?
        ! Only release it if the data belongs to this matrix.
        IF (brelease .AND. &
            (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .EQ. 0)) THEN
          IF (rmatrix%h_Da .NE. ST_NOHANDLE) CALL storage_free (rmatrix%h_Da)
        END IF
        
        ! Reset the handles
        rmatrix%h_Da = ST_NOHANDLE
        
      END SELECT

      ! Reset the ownership-status
      rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,&
                                  NOT(LSYSSC_MSPEC_CONTENTISCOPY))
    
    END SUBROUTINE

    !--------------------------------------------------------
    ! Auxiliary routine: Shares the content between rsourceMatrix
    ! and rdestMatrix
    
    SUBROUTINE shareContent(rsourceMatrix, rdestMatrix)
    
    ! The source matrix 
    TYPE(t_matrixScalar), INTENT(IN) :: rsourceMatrix

    ! The destination matrix 
    TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
    
      ! Remove the old content - if there is any.
      CALL removeContent(rdestMatrix, .TRUE.)

      ! Overwrite structural data
      rdestMatrix%dscaleFactor = rsourceMatrix%dscaleFactor
      rdestMatrix%cdataType    = rsourceMatrix%cdataType

      ! Which source matrix do we have?  
      SELECT CASE (rsourceMatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
        rdestMatrix%h_Da = rsourceMatrix%h_Da
      END SELECT
      
      ! Indicate via the matrixSpec-flag that we are not
      ! the owner of the structure.
      rdestMatrix%imatrixSpec = IOR(rdestMatrix%imatrixSpec,&
                                    LSYSSC_MSPEC_CONTENTISCOPY)
    
    END SUBROUTINE

    !--------------------------------------------------------
    ! Auxiliary routine: Copy static content-related
    ! information (scale factor,...) rsourceMatrix
    ! to rdestMatrix
    
    SUBROUTINE copyStaticContent(rsourceMatrix, rdestMatrix)
    
    ! The source matrix 
    TYPE(t_matrixScalar), INTENT(IN) :: rsourceMatrix

    ! The destination matrix 
    TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
    
      ! Overwrite structural data
      rdestMatrix%dscaleFactor = rsourceMatrix%dscaleFactor

    END SUBROUTINE

    !--------------------------------------------------------
    ! Auxiliary routine: Copy the content of rsourceMatrix
    ! to rdestMatrix
    
    SUBROUTINE copyContent(rsourceMatrix, rdestMatrix)
    
    ! The source matrix 
    TYPE(t_matrixScalar), INTENT(IN) :: rsourceMatrix

    ! The destination matrix 
    TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
    
      ! local variables
      LOGICAL :: bremove
    
      ! Overwrite structural data
      CALL copyStaticContent(rsourceMatrix, rdestMatrix)
      
      ! Check the content if it exists and if it has the right
      ! size - then we can overwrite!
      bRemove = .FALSE.
      
      ! But at first, check if rdestMatrix is the owner of the matrix
      ! structure:
      IF (IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .NE. 0) THEN
        ! No, the content belongs to someone else, but we want to have
        ! our own. Detach the foreign matrix content.
        bremove = .TRUE.
        
      ELSE
        
        ! Ok, rdestMatrix owns some matrix content.
        ! Check the structure if it's large enough so we can overwrite
        ! it. If the size of the arrays is not equal to those
        ! in the source matrix, release the content and allocate a new one.
        
        ! Which source matrix do we have?  
        SELECT CASE (rsourceMatrix%cmatrixFormat)
        CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
        
          IF (rdestMatrix%h_Da .NE. ST_NOHANDLE) THEN
        
            CALL storage_getsize (rdestMatrix%h_Da,isize)
            bremove = bremove .OR. (isize .NE. rdestMatrix%NA)
          
            ! Check the data type
            bremove = bremove .OR. (rdestMatrix%cdataType .NE. rsourceMatrix%cdataType)
          
          ELSE
          
            ! Remove any partial information if there is any
            bremove = .TRUE.
            
          END IF
          
        END SELECT
      
      END IF
    
      ! Remove the old content - if we should do so
      IF (bremove) CALL removeContent(rdestMatrix, .TRUE.)

      rdestMatrix%cdataType = rsourceMatrix%cdataType

      ! Duplicate content data from the source matrix.
      ! Storage_copy allocates new memory as the handles are all
      ! set to ST_NOHANDLE with the above removeStructure!
      !
      ! Which source matrix do we have?  
      SELECT CASE (rsourceMatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
        CALL storage_copy (rsourceMatrix%h_Da,rdestMatrix%h_Da)
      END SELECT
      
      ! Indicate via the matrixSpec-flag that we are the owner of the structure.
      rdestMatrix%imatrixSpec = IAND(rdestMatrix%imatrixSpec,&
                                    NOT(LSYSSC_MSPEC_CONTENTISCOPY))
    
    END SUBROUTINE
  
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_releaseVector (rvector)
  
!<description>
  ! Releases a vector from memory. The vector structure is cleaned up.
  ! Remark: The memory associated to the vector is only released if this
  !  vector structure is the owner of the data.
!</description>
  
!<inputoutput>
  
  ! Vector to release.
  TYPE(t_vectorScalar), INTENT(INOUT)               :: rvector
  
!</inputoutput>

!</subroutine>

  IF (rvector%h_Ddata .EQ. ST_NOHANDLE) THEN
    PRINT *,'lsyssc_releaseVector warning: releasing unused vector.'
  END IF
  
  ! Clean up the data structure.
  ! Don't release the vector data if the handle belongs to another
  ! vector!
  IF ((.NOT. rvector%bisCopy) .AND. (rvector%h_Ddata .NE. ST_NOHANDLE)) THEN
    CALL storage_free(rvector%h_Ddata)
  ELSE
    rvector%h_Ddata = ST_NOHANDLE
  END IF
  rvector%NEQ = 0
  rvector%cdataType = ST_NOHANDLE
  rvector%iidxFirstEntry = 1
  rvector%bisCopy = .FALSE.
  rvector%isortStrategy = 0
  rvector%h_IsortPermutation = ST_NOHANDLE
  rvector%p_rspatialDiscretisation => NULL()
   
  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_releaseMatrix (rmatrix)
  
!<description>
  ! Releases a matrix from memory. The structure is cleaned up.
  ! All data vectors belonging to the structure are released;
  ! data structures belonging to other matrices are not.
!</description>
  
!<inputoutput>
  
  ! Matrix to release.
  TYPE(t_matrixScalar), INTENT(INOUT)               :: rmatrix
  
!</inputoutput>

!</subroutine>

  ! Which handles do we have to release?
  !
  ! Release the matrix data if the handle is not a copy of another matrix
  IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .EQ. 0) THEN
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
      ! Release matrix data, structure 9,7
      IF (rmatrix%h_DA .NE. ST_NOHANDLE) THEN
        CALL storage_free(rmatrix%h_DA)
      END IF
    END SELECT
  END IF
  
  ! Release the structure if it doesn't belong to another vector
  IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .EQ. 0) THEN
    ! Release matrix structure
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9)
      IF (rmatrix%h_Kcol .NE. ST_NOHANDLE)      CALL storage_free(rmatrix%h_Kcol)
      IF (rmatrix%h_Kld .NE. ST_NOHANDLE)       CALL storage_free(rmatrix%h_Kld)
      IF (rmatrix%h_Kdiagonal .NE. ST_NOHANDLE) CALL storage_free(rmatrix%h_Kdiagonal)
    CASE (LSYSSC_MATRIX7)
      IF (rmatrix%h_Kcol .NE. ST_NOHANDLE) CALL storage_free(rmatrix%h_Kcol)
      IF (rmatrix%h_Kld .NE. ST_NOHANDLE)  CALL storage_free(rmatrix%h_Kld)
    END SELECT
  END IF
  
  ! Clean up the rest
  rmatrix%h_DA        = ST_NOHANDLE
  rmatrix%h_Kcol      = ST_NOHANDLE
  rmatrix%h_Kld       = ST_NOHANDLE
  rmatrix%h_Kdiagonal = ST_NOHANDLE
  rmatrix%cdataType   = ST_DOUBLE
  rmatrix%NA  = 0
  rmatrix%NEQ = 0
  rmatrix%NCOLS = 0
  rmatrix%dScaleFactor = 1.0_DP
  rmatrix%isortStrategy = 0
  rmatrix%h_IsortPermutation = ST_NOHANDLE
  rmatrix%p_rspatialDiscretisation => NULL()

  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_clearMatrix (rmatrix)
  
!<description>
  ! Clears the entries in a matrix. All entries are overwritten with 0.0.
!</description>
  
!<inputoutput>
  
  ! Matrix to release.
  TYPE(t_matrixScalar), INTENT(INOUT)               :: rmatrix
  
!</inputoutput>

!</subroutine>

  IF (rmatrix%NEQ .LE. 0) RETURN ! Empty matrix

  ! Which matrix type do we have?
  SELECT CASE (rmatrix%cmatrixFormat)
  CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
    ! Get the handle, the associated memory and clear that.
    IF (rmatrix%h_DA .NE. ST_NOHANDLE) THEN
      CALL storage_clear(rmatrix%h_DA)
    END IF
  END SELECT
  
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_convertMatrix (rmatrix,cmatrixFormat,bresortEntries)
  
!<description>
  ! Tries to convert a matrix rmatrix into a different matrix format.
  ! If the matrix cannot be converted (due to format incompatibility),
  ! an error is thrown.
!</description>
  
!<input>
  ! Destination format of the matrix. One of the LSYSSC_MATRIXx constants.
  INTEGER, INTENT(IN)                 :: cmatrixFormat
  
  ! OPTIONAL: Whether to resort the entries of the matrix. Standard = TRUE.
  ! If set to FALSE, the entries are not resorted; helpful to set up
  ! convert only the structure, not the entries of a matrix.
  LOGICAL, INTENT(IN), OPTIONAL       :: bresortEntries
!</input>
  
!<inputoutput>
  ! Matrix to convert.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER(PREC_VECIDX) :: i
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld,p_Kdiagonal
  LOGICAL :: bentries

  ! Matrix is already in that format.
  IF (rmatrix%cmatrixFormat .EQ. cmatrixFormat) RETURN
  
  ! Empty matrix
  IF (rmatrix%NEQ .LE. 0) RETURN 
  
  bentries = .TRUE.
  IF (PRESENT(bresortEntries)) bentries = bresortEntries

  ! Which matrix type do we have?
  SELECT CASE (rmatrix%cmatrixFormat)
  CASE (LSYSSC_MATRIX9)
  
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX7)
    
      CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
      CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)
      CALL storage_getbase_int (rmatrix%h_Kdiagonal,p_Kdiagonal)

      ! Check that the matrix can be converted. There's a format error
      ! if there's no diagonal element.
      DO i=1,rmatrix%NEQ
        IF (p_Kdiagonal(i) .NE. i) THEN
          PRINT *,'lsyssc_convertMatrix: incompatible Format-9 matrix!'
          STOP
        END IF
      END DO
      
      IF ((.NOT. bentries) .OR. (rmatrix%h_DA .EQ. ST_NOHANDLE)) THEN
      
        ! No matrix entries, only resort the structure
        CALL lsyssc_unsortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)
        
        ! Release diagonal
        CALL storage_free (rmatrix%h_Kdiagonal)

        rmatrix%cmatrixFormat = LSYSSC_MATRIX7
      
      ELSE
    
        ! Convert from structure 7 to structure 9. Use the sortCSRxxxx 
        ! routine below.
        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
        
          CALL storage_getbase_double (rmatrix%h_DA,p_Ddata)
          CALL lsyssc_unsortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ, p_Ddata)
          ! Release diagonal
          CALL storage_free (rmatrix%h_Kdiagonal)

          rmatrix%cmatrixFormat = LSYSSC_MATRIX7
          
        CASE DEFAULT
          PRINT *,'lsyssc_convertMatrix: Unsupported data type!'
          STOP
        END SELECT

      END IF

    CASE DEFAULT
      PRINT *,'lsyssc_convertMatrix: Cannot convert matrix!'
      STOP
    END SELECT

  CASE (LSYSSC_MATRIX7)
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9)

      ! Convert from structure 7 to structure 9. Use the sortCSRxxxx 
      ! routine below.
      CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
      CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)

      ! Create a pointer to the diagonal
      CALL storage_new ('lsyssc_convertMatrix', 'Kdiagonal', &
            rmatrix%NEQ, ST_INT, rmatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_int (rmatrix%h_Kdiagonal,p_Kdiagonal)

      IF ((.NOT. bentries) .OR. (rmatrix%h_DA .EQ. ST_NOHANDLE)) THEN
      
        ! No matrix entries, only resort the structure
        CALL lsyssc_sortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)
        
        rmatrix%cmatrixFormat = LSYSSC_MATRIX9
      
      ELSE

        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
          CALL storage_getbase_double (rmatrix%h_DA,p_Ddata)
          CALL lsyssc_sortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ, p_Ddata)
          
          rmatrix%cmatrixFormat = LSYSSC_MATRIX9
          
        CASE DEFAULT
          PRINT *,'lsyssc_convertMatrix: Unsupported data type!'
          STOP
        END SELECT
      
      END IF
      
    CASE DEFAULT
      PRINT *,'lsyssc_convertMatrix: Cannot convert matrix!'
      STOP
    END SELECT

  CASE DEFAULT
    PRINT *,'lsyssc_convertMatrix: Cannot convert matrix!'
    STOP
  END SELECT
  
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_rebuildKdiagonal (Kcol, Kld, Kdiagonal, neq)

!<description>
  ! Internal auxiliary routine.
  ! Rebuilds the Kdiagonal-array of a structure-9 matrix.
!</description>

!<input>
  ! Row structure in the matrix
  INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kld

  ! Column structure of the matrix
  INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Kcol

  ! Dimension of the matrix
  INTEGER(I32), INTENT(IN) :: neq
!</input>

!<output>
  ! Pointers to the diagonal entries of the matrix.
  INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(OUT) :: Kdiagonal
!</output>

!</subroutine>

  ! local variables
  INTEGER(I32) :: i, j

  ! loop through each row
  DO i = 1, neq

    ! Loop through each column in this row.
    ! Shift every entry until the diagonal is reached.
    DO j = Kld(i)+1, Kld(i+1)-1

      ! Check if we reached the position of the diagonal entry...
      IF (Kcol(j)>i) EXIT

    END DO

    ! Save the position of the diagonal entry
    Kdiagonal(i) = j

  END DO
          
  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_sortCSRdouble (Kcol, Kld, Kdiagonal, neq, Da)

!<description>
  ! Internal auxiliary routine.
  ! Sorts the entries in each row of the matrix in ascending order.
  ! The input matrix is assumed to be in storage technique 7 with the 
  ! first element in each row to be the diagonal element!
  ! Creates the Kdiagonal array for a structure-9 matrix.
  !
  ! Double precision version
!</description>

!<input>
  ! Row pointer in the matrix
  INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kld

  ! Dimension of the matrix
  INTEGER(I32), INTENT(IN) :: neq
!</input>

!<inputoutput>
  ! OPTIONAL:
  ! On input:  the matrix entries to be resorted,
  ! On output: the resorted matrix entries
  REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Da
  
  ! On input:  the column numbers to be resorted,
  ! On output: the resorted column numbers
  INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Kcol
!</inputoutput>

!<output>
  ! Pointers to the diagonal entries of the matrix.
  INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(OUT) :: Kdiagonal
!</output>

!</subroutine>

  ! local variables
  REAL(DP) :: aux
  INTEGER(I32) :: i, j

  IF (PRESENT(Da)) THEN

    ! loop through each row
    DO i = 1, neq

      ! Take the diagonal element
      aux = Da(Kld(i))

      ! Loop through each column in this row.
      ! Shift every entry until the diagonal is reached.
      DO j = Kld(i)+1, Kld(i+1)-1

        ! Check if we reached the position of the diagonal entry...
        IF (Kcol(J)>i) EXIT

        Kcol(j-1) = KCOL(j)
        Da(j-1) = Da(j)

      END DO

      ! If we have reached the diagonal, we can stop and save our
      ! diagonal entry from the first position there. The rest of the
      ! line is in ascending order according to the specifications of
      ! storage technique 7.

      Kcol(j-1) = i
      Da(j-1) = aux
      
      ! Save the position of the diagonal entry
      Kdiagonal(i) = j

    END DO
          
  ELSE

    ! loop through each row
    DO i = 1, neq

      ! Loop through each column in this row.
      ! Shift every entry until the diagonal is reached.
      DO j = Kld(i)+1, Kld(i+1)-1

        ! Check if we reached the position of the diagonal entry...
        IF (Kcol(J)>i) EXIT

        Kcol(j-1) = KCOL(j)

      END DO

      ! If we have reached the diagonal, we can stop and save our
      ! diagonal entry from the first position there. The rest of the
      ! line is in ascending order according to the specifications of
      ! storage technique 7.

      Kcol(j-1) = i
      
      ! Save the position of the diagonal entry
      Kdiagonal(i) = j

    END DO
          
  END IF
          
  END SUBROUTINE
  
  !****************************************************************************

!<subroutine>  
  
  SUBROUTINE lsyssc_unsortCSRdouble (Kcol, Kld, Kdiagonal, neq, Da)

!<description>
  ! Internal auxiliary routine.
  ! Unorts the entries in each row of the matrix in ascending order.
  ! This searches in each row of a matrix for the diagonal element
  ! and shifts it to the front.
  !
  ! Double precision version
!</description>

!<input>
  ! Row pointer in the matrix
  INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kld

  ! Dimension of the matrix
  INTEGER(I32), INTENT(IN) :: neq

  ! Pointers to the diagonal entries of the matrix.
  INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kdiagonal
!</input>

!<inputoutput>
  ! OPTIONAL:
  ! On input:  the matrix entries to be resorted,
  ! On output: the resorted matrix entries
  REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Da
  
  ! On input:  the column numbers to be resorted,
  ! On output: the resorted column numbers
  INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: Kcol
!</inputoutput>

!</subroutine>

  ! local variables
  REAL(DP) :: aux
  INTEGER(I32) :: i, j

  IF (PRESENT(Da)) THEN

    ! loop through each row
    DO i = 1, neq

      ! Take the diagonal element
      aux = Da(Kdiagonal(i))

      ! Loop through each column in this row.
      ! Shift every entry one element to the right.
      DO j = Kdiagonal(i)-1,Kld(i)+1,-1

        Kcol(j) = Kcol(j-1)
        Da(j) = Da(j-1)

      END DO
      
      ! Put the diagonal to the front.
      Kcol(Kdiagonal(i)) = i
      Da(Kdiagonal(i)) = aux
      
    END DO
    
  ELSE
    ! loop through each row
    DO i = 1, neq

      ! Loop through each column in this row.
      ! Shift every entry one element to the right.
      DO j = Kdiagonal(i)-1,Kld(i)+1,-1
        Kcol(j) = Kcol(j-1)
      END DO
      
      ! Put the diagonal to the front.
      Kcol(Kdiagonal(i)) = i
      
    END DO

  END IF
          
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_sortVectorInSitu (rvector,rtemp,isortStrategy,h_IsortPermutation)
  
!<description>
  ! Resorts the entries of the given vector rvector or unsorts it.
  ! rvector is the vector to be resorted, rtemp is a temporary vector and
  ! isortStrategy is a type flag identifying the sorting algorithm.
  !
  ! This routine can also be used to assign an unsorted vector a sorting 
  ! strategy without actually resorting it. For this purpose, isortStrategy
  ! should be '- SSTRAT_xxxx' and h_IsortPermutation a handle to a
  ! sorting strategy permutation.
!</description>
  
!<inputoutput>
  ! Vector to resort
  TYPE(t_vectorScalar), INTENT(INOUT)               :: rvector

  ! A temporary vector. Must be of the same data type as rvector.
  ! Must be at least as large as rvector.
  TYPE(t_vectorScalar), INTENT(INOUT)               :: rtemp
!</inputoutput>

!<input>
  ! Identifier for the sorting strategy to apply to the vector.
  ! This is usually one of the SSTRAT_xxxx constants from the module
  ! 'sortstrategy', although it's actually used here as follows:
  ! <=0: Calculate the unsorted vector
  !  >0: Resort the vector according to a permutation;
  !      this is either the permutation specified in the vector
  !      or that one identified by h_IsortPermutation
 INTEGER, INTENT(IN)                                :: isortStrategy
  
  ! OPTIONAL: Handle to permutation to use for resorting the vector.
  ! 
  ! The array must be of the form
  !    array [1..2*NEQ] of integer
  ! with entries (1..NEQ)       = permutation
  ! and  entries (NEQ+1..2*NEQ) = inverse permutation.
  !
  ! If not specified, the associated permutation rvector%h_IsortPermutation
  ! is used.
  ! If specified and isortStrategy>0, the vector is resorted according to 
  ! the permutation h_IsortPermutation (probably unsorted before if necessary).
  !
  ! In any case, the associated permutation of the vector
  ! rvector%h_IsortPermutation is changed to h_IsortPermutation.
  ! Remark: The memory associated to the previous sorting strategy
  !  in this case is not released automatically; this has to be done by the
  !  application!
  INTEGER, INTENT(IN), OPTIONAL                     :: h_IsortPermutation
!</input> 

!</subroutine>

  ! local variables
  INTEGER :: h_Iperm
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata,p_Fdata2
  INTEGER(PREC_VECIDX) :: NEQ
  
    ! Desired sorting strategy and currently active sorting strategy identical?
    IF (.NOT. PRESENT(h_IsortPermutation)) THEN
    
      IF (isortStrategy .EQ. rvector%isortStrategy) RETURN
      IF ((isortStrategy .LE. 0) .AND. (rvector%isortStrategy .LE. 0)) RETURN
    
    ELSE
    
      IF ((isortStrategy .LE. 0) .AND. (rvector%isortStrategy .LE. 0)) THEN
        IF (h_IsortPermutation .NE. rvector%h_IsortPermutation) THEN
          ! Vector is unsorted and should stay unsorted, but
          ! permutation should change.
          rvector%isortStrategy = isortStrategy
          rvector%h_IsortPermutation = h_IsortPermutation
        END IF
        RETURN
      END IF
      IF ((isortStrategy .GT. 0) .AND. (rvector%isortStrategy .GT. 0) .AND.&
          (h_IsortPermutation .EQ. rvector%h_IsortPermutation)) RETURN
    
    END IF
    
    NEQ = rvector%NEQ
    
    ! Get pointers to the vector data
    SELECT CASE (rvector%cdataType)
    CASE (ST_DOUBLE)
      CALL lsyssc_getbase_double(rvector,p_Ddata)
      CALL lsyssc_getbase_double(rtemp,p_Ddata2)
    CASE (ST_SINGLE)   
      CALL lsyssc_getbase_single(rvector,p_Fdata)
      CALL lsyssc_getbase_single(rtemp,p_Fdata2)
    CASE DEFAULT
      PRINT *,'lsyssc_sortVectorInSitu: unsuppported data type'
      STOP
    END SELECT

    
    ! Sort the vector back?
    IF (isortStrategy .LE. 0) THEN
      ! Do it - with the associated permutation.
      CALL storage_getbase_int(rvector%h_IsortPermutation,p_Iperm)
      
      SELECT CASE (rvector%cdataType)
      CASE (ST_DOUBLE)
        ! Copy the entries to the temp vector
        CALL lalg_copyVectorDble (p_Ddata,p_Ddata2)
        ! Then sort back. Use the inverse permutation.
        CALL lalg_vectorSortDble (p_Ddata2,p_Ddata,p_Iperm(NEQ+1:NEQ*2))
        
      CASE (ST_SINGLE)   
        ! Copy the entries to the temp vector
        CALL lalg_copyVectorSngl (p_Fdata,p_Fdata2)
        ! Then sort back. Use the inverse permutation.
        CALL lalg_vectorSortSngl (p_Fdata2,p_Fdata,p_Iperm(NEQ+1:NEQ*2))
        
      CASE DEFAULT
        PRINT *,'lsyssc_sortVectorInSitu: unsuppported data type'
        STOP
        
      END SELECT
        
      ! Inform the vector about which sorting strategy we now use.
      rvector%isortStrategy = isortStrategy
      RETURN
    END IF
    
    ! Get the actual sorting strategy.
    h_Iperm = rvector%h_IsortPermutation
    IF (PRESENT(h_IsortPermutation)) h_Iperm = h_IsortPermutation 
    
    ! Do we have to sort back before resorting?
    IF ((h_Iperm .NE. rvector%h_IsortPermutation) .AND. &
        (rvector%h_IsortPermutation .NE. ST_NOHANDLE)) THEN
        
      ! Sort back at first - with the associated permutation
      CALL storage_getbase_int(rvector%h_IsortPermutation,p_Iperm)
      
      SELECT CASE (rvector%cdataType)
      CASE (ST_DOUBLE)
        ! Copy the entries to the temp vector
        CALL lalg_copyVectorDble (p_Ddata,p_Ddata2)
        ! Then sort back. Use the inverse permutation.
        CALL lalg_vectorSortDble (p_Ddata2,p_Ddata,p_Iperm(NEQ+1:NEQ*2))
        
      CASE (ST_SINGLE)   
        ! Copy the entries to the temp vector
        CALL lalg_copyVectorSngl (p_Fdata,p_Fdata2)
        ! Then sort back. Use the inverse permutation.
        CALL lalg_vectorSortSngl (p_Fdata2,p_Fdata,p_Iperm(NEQ+1:NEQ*2))
        
      CASE DEFAULT
        PRINT *,'lsyssc_sortVectorInSitu: unsuppported data type'
        STOP
        
      END SELECT
      
      ! Change the sorting strategy in the vector to the one
      ! we are now going to use. Throw away the old handle.
      rvector%h_IsortPermutation = h_Iperm
      
    END IF
    
    ! Now sort the vector according to h_Iperm
    CALL storage_getbase_int(h_Iperm,p_Iperm)
    
    SELECT CASE (rvector%cdataType)
    CASE (ST_DOUBLE)
      ! Copy the entries to the temp vector
      CALL lalg_copyVectorDble (p_Ddata,p_Ddata2)
      ! Then do the sorting with the given permutation.
      CALL lalg_vectorSortDble (p_Ddata2,p_Ddata,p_Iperm(1:NEQ))
      
    CASE (ST_SINGLE)   
      ! Copy the entries to the temp vector
      CALL lalg_copyVectorSngl (p_Fdata,p_Fdata2)
      ! Then do the sorting with the given permutation.
      CALL lalg_vectorSortSngl (p_Fdata2,p_Fdata,p_Iperm(1:NEQ))
      
    CASE DEFAULT
      PRINT *,'lsyssc_sortVectorInSitu: unsuppported data type'
      STOP
      
    END SELECT
    
    ! Inform the vector about which sorting strategy we now use.
    rvector%isortStrategy = isortStrategy

    ! If h_IsortPermutation was given, change the permutation
    IF (PRESENT(h_IsortPermutation)) rvector%h_IsortPermutation = h_IsortPermutation
  
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_vectorActivateSorting (rvector,bsort,rtemp)
  
!<description>
  ! Resorts the entries of the given vector rvector or unsorts it
  ! according to the resorting strategy associated to rvector.
!</description>
  
!<inputoutput>
  ! Vector to resort. The sorting strategy must have been attached to
  ! rvector before with lsyssc_sortVectorInSitu, otherwise nothing happens.
  TYPE(t_vectorScalar), INTENT(INOUT)               :: rvector

  ! OPTIONAL: A temporary vector. 
  ! Must be of the same data type as rvector. Must be at least as 
  ! large as rvector. If not specified, a temporary vector is created
  ! and released automatically on the heap.
  TYPE(t_vectorScalar), INTENT(INOUT), TARGET, OPTIONAL :: rtemp
!</inputoutput>

!<input>
  ! Whether to sort or unsort.
  ! =TRUE : Activate sorting (if not activated)
  ! =FALSE: Unsort vector (if sorted)
 LOGICAL, INTENT(IN) :: bsort
!</input> 

!</subroutine>

  ! local variables
  TYPE(t_vectorScalar), POINTER :: p_rtemp
  TYPE(t_vectorScalar), TARGET :: rtempLocal
  
  ! Cancel if there's nothing to do.
  IF (rvector%isortStrategy .EQ. 0) RETURN
  IF ((.NOT. bsort) .AND. (rvector%isortStrategy .LE. 0)) RETURN
  IF (bsort .AND. (rvector%isortStrategy .GT. 0)) RETURN
  
  ! Temporary vector available? If not, create a new one based on rvector.
  IF (PRESENT(rtemp)) THEN
    p_rtemp => rtemp
  ELSE
    p_rtemp => rtempLocal
    CALL lsyssc_duplicateVector (rvector,rtempLocal)
  END IF

  ! Perform the sorting or unsorting
  IF (bsort) THEN
    CALL lsyssc_sortVectorInSitu (rvector,p_rtemp,ABS(rvector%isortStrategy))
  ELSE
    CALL lsyssc_sortVectorInSitu (rvector,p_rtemp,-ABS(rvector%isortStrategy))
  END IF
  
  ! Remove the temp vector if it's ours.
  IF (.NOT. PRESENT(rtemp)) THEN
    CALL lsyssc_releaseVector (rtempLocal)
  END IF
  
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_sortMatrix (rmatrix,bsortEntries,&
                                isortStrategy,h_IsortPermutation)
  
!<description>
  ! Matrix sorting. Sorts the entries and/or structure of a given
  ! matrix or unsorts them according to a permutation.
  !
  ! This routine can also be used to assign an unsorted vector a sorting 
  ! strategy without actually resorting it. For this purpose, isortStrategy
  ! should be '- SSTRAT_xxxx' and h_IsortPermutation a handle to a
  ! sorting strategy permutation.
!</description>
  
!<inputoutput>
  ! Vector to resort
  TYPE(t_matrixScalar), INTENT(INOUT)               :: rmatrix
!</inputoutput>

!<input>
  ! Sort the entries or only the structure of the matrix.
  ! = FALSE: Only sort the structure of the matrix,
  ! = TRUE : Sort both, entries and structure of the matrix.
  LOGICAL bsortEntries

  ! Identifier for the sorting strategy to apply to the matrix.
  ! This is usually one of the SSTRAT_xxxx constants from the module
  ! 'sortstrategy', although it's actually used here as follows:
  ! <=0: Calculate the unsorted matrix
  !  >0: Resort the vector according to a permutation;
  !      this is either the permutation specified in the vector
  !      or that one identified by h_IsortPermutation
 INTEGER, INTENT(IN)                                :: isortStrategy
  
  ! OPTIONAL: Handle to permutation to use for resorting the matrix.
  ! 
  ! The array must be of the form
  !    array [1..2*NEQ] of integer
  ! with entries (1..NEQ)       = permutation
  ! and  entries (NEQ+1..2*NEQ) = inverse permutation.
  !
  ! If not specified, the associated permutation rmatrix%h_IsortPermutation
  ! is used.
  ! If specified and isortStrategy>0, the matrix is resorted according to 
  ! the permutation h_IsortPermutation (probably unsorted before if necessary).
  !
  ! In any case, the associated permutation of the matrix
  ! rmatrix%h_IsortPermutation is changed to h_IsortPermutation.
  ! Remark: The memory associated to the previous sorting strategy
  !  in this case is not released automatically; this has to be done by the
  !  application!
  INTEGER, INTENT(IN), OPTIONAL                     :: h_IsortPermutation
!</input> 

!</subroutine>

  ! local variables
  INTEGER :: h_Iperm
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm
  INTEGER(PREC_VECIDX) :: NEQ
  
    ! Desired sorting strategy and currently active sorting strategy identical?
    IF (.NOT. PRESENT(h_IsortPermutation)) THEN
    
      IF (isortStrategy .EQ. rmatrix%isortStrategy) RETURN
      IF ((isortStrategy .LE. 0) .AND. (rmatrix%isortStrategy .LE. 0)) RETURN
    
    ELSE
    
      IF ((isortStrategy .LE. 0) .AND. (rmatrix%isortStrategy .LE. 0)) THEN
        IF (h_IsortPermutation .NE. rmatrix%h_IsortPermutation) THEN
          ! Matrix is unsorted and should stay unsorted, but
          ! permutation should change.
          rmatrix%isortStrategy = isortStrategy
          rmatrix%h_IsortPermutation = h_IsortPermutation
        END IF
        RETURN
      END IF
      IF ((isortStrategy .GT. 0) .AND. (rmatrix%isortStrategy .GT. 0) .AND.&
          (h_IsortPermutation .EQ. rmatrix%h_IsortPermutation)) RETURN
    
    END IF
    
    NEQ = rmatrix%NEQ
    
    ! Sort the matrix back?
    IF (isortStrategy .LE. 0) THEN
      ! Get the permutation that describes how to resort the matrix:
      CALL storage_getbase_int(rmatrix%h_IsortPermutation,p_Iperm)
      
      ! Exchange the roles of the first and second part of p_Iperm.
      ! This makes the permutation to the inverse permutation and vice versa.
      ! Call the resort-subroutine with that to do the actual sorting.
      CALL do_matsort (rmatrix,p_Iperm(NEQ+1:NEQ*2),p_Iperm(1:NEQ),bsortEntries)
      
      ! Inform the vector about which sorting strategy we now use.
      rmatrix%isortStrategy = isortStrategy
      RETURN
    END IF
    
    ! Get the actual sorting strategy.
    h_Iperm = rmatrix%h_IsortPermutation
    IF (PRESENT(h_IsortPermutation)) h_Iperm = h_IsortPermutation 
    
    ! Do we have to sort back before resorting?
    IF ((h_Iperm .NE. rmatrix%h_IsortPermutation) .AND. &
        (rmatrix%h_IsortPermutation .NE. ST_NOHANDLE)) THEN

      ! Sort back at first - with the associated permutation
      CALL storage_getbase_int(rmatrix%h_IsortPermutation,p_Iperm)
      
      ! Exchange the roles of the first and second part of p_Iperm.
      ! This makes the permutation to the inverse permutation and vice versa.
      ! Call the resort-subroutine with that to do the actual sorting.
      CALL do_matsort (rmatrix,p_Iperm(NEQ+1:NEQ*2),p_Iperm(1:NEQ),bsortEntries)
    END IF
    
    ! Now sort the vector according to h_Iperm
    CALL storage_getbase_int(h_Iperm,p_Iperm)

    ! This time, we don't exchange the roles of the permutation and
    ! its inverse :-)
    CALL do_matsort (rmatrix,p_Iperm(1:NEQ),p_Iperm(NEQ+1:NEQ*2),bsortEntries)
    
    ! Inform the vector about which sorting strategy we now use.
    rmatrix%isortStrategy = isortStrategy
    
    ! If h_IsortPermutation was given, change the permutation
    IF (PRESENT(h_IsortPermutation)) rmatrix%h_IsortPermutation = h_IsortPermutation
  
  CONTAINS
    
    !----------------------------------------------------------------
    ! Sort matrix or matrix entries.
    ! This calls the actual resorting routine, depending
    ! on the information tags in the matrix.
    
    SUBROUTINE do_matsort (rmatrix,Itr1,Itr2,bsortEntries)
    
    ! The matrix to be resorted
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
    
    ! The transformation to use
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Itr1
    
    ! The inverse transformation
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Itr2
    
    ! TRUE  = sort matrix structure + entries
    ! FALSE = sort only matrix structure
    LOGICAL, INTENT(IN) :: bsortEntries
    
    ! local variables
    
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld,p_KldTmp,p_Kdiag
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol,p_KcolTmp
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_DdataTmp
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata,p_FdataTmp
    TYPE(t_matrixScalar) :: rtempMatrix
    INTEGER(PREC_VECIDX) :: NEQ
    
      NEQ = rmatrix%NEQ
    
      ! Which matrix configuration do we have?
      SELECT CASE (rmatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9)
        ! Duplicate the matrix before resorting.
        ! We need either a copy only of the structure or of the full matrix.
        IF (.NOT. bsortEntries) THEN
          CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                       LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)
        ELSE
          CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                       LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
        END IF
        
        ! Get the structure of the original and the temporary matrix
        CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
        CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)
        CALL storage_getbase_int (rmatrix%h_Kdiagonal,p_Kdiag)
        CALL storage_getbase_int (rtempMatrix%h_Kcol,p_KcolTmp)
        CALL storage_getbase_int (rtempMatrix%h_Kld,p_KldTmp)
        
        IF (.NOT. bsortEntries) THEN
        
          ! Sort only the structure of the matrix, keep the entries
          ! unchanged.
          CALL lsyssc_sortMat9Struc (p_Kcol, p_KcolTmp, p_Kld, p_KldTmp, &
                                     p_Kdiag, Itr1, Itr2, NEQ)        
        
        ELSE
        
          ! Sort structure + entries
          SELECT CASE (rmatrix%cdataType)
          CASE (ST_DOUBLE)
            ! Double precision version
            CALL storage_getbase_double (rmatrix%h_Da,p_Ddata)
            CALL storage_getbase_double (rtempMatrix%h_Da,p_DdataTmp)
            CALL lsyssc_sortMat9_double (p_Ddata,p_DdataTmp,p_Kcol, p_KcolTmp, &
                                         p_Kld, p_KldTmp, p_Kdiag, &
                                         Itr1, Itr2, NEQ)        
          CASE (ST_SINGLE)
            ! Single precision version
            CALL storage_getbase_single (rmatrix%h_Da,p_Fdata)
            CALL storage_getbase_single (rtempMatrix%h_Da,p_FdataTmp)
            CALL lsyssc_sortMat9_single (p_Fdata,p_FdataTmp,p_Kcol, p_KcolTmp, &
                                         p_Kld, p_KldTmp, p_Kdiag, &
                                         Itr1, Itr2, NEQ)        
          CASE DEFAULT
            PRINT *,'lsyssc_sortMatrix: Unsupported data type.'
            STOP
          END SELECT
        END IF

        ! Remove temp matrix
        CALL lsyssc_releaseMatrix(rtempMatrix)
        
      CASE (LSYSSC_MATRIX7)
        ! Duplicate the matrix before resorting.
        ! We need either a copy only of the structure or of the full matrix.
        IF (.NOT. bsortEntries) THEN
          CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                       LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)
        ELSE
          CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                       LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
        END IF
        
        ! Get the structure of the original and the temporary matrix
        CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
        CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)
        CALL storage_getbase_int (rtempMatrix%h_Kcol,p_KcolTmp)
        CALL storage_getbase_int (rtempMatrix%h_Kld,p_KldTmp)
        
        IF (.NOT. bsortEntries) THEN
        
          ! Sort only the structure of the matrix, keep the entries
          ! unchanged.
          CALL lsyssc_sortMat7Struc (p_Kcol, p_KcolTmp, p_KldTmp, p_KldTmp, &
                                     Itr1, Itr2, NEQ)        
        
        ELSE
        
          ! Sort structure + entries
          SELECT CASE (rmatrix%cdataType)
          CASE (ST_DOUBLE)
            ! Double precision version
            CALL storage_getbase_double (rmatrix%h_Da,p_Ddata)
            CALL storage_getbase_double (rtempMatrix%h_Da,p_DdataTmp)
            CALL lsyssc_sortMat7_double (p_Ddata,p_DdataTmp,p_Kcol, p_KcolTmp, &
                                         p_Kld, p_KldTmp, &
                                         Itr1, Itr2, NEQ)        
          CASE (ST_SINGLE)
            ! Single precision version
            CALL storage_getbase_single (rmatrix%h_Da,p_Fdata)
            CALL storage_getbase_single (rtempMatrix%h_Da,p_FdataTmp)
            CALL lsyssc_sortMat7_single (p_Fdata,p_FdataTmp,p_Kcol, p_KcolTmp, &
                                         p_Kld, p_KldTmp, &
                                         Itr1, Itr2, NEQ)        
          CASE DEFAULT
            PRINT *,'lsyssc_sortMatrix: Unsupported data type.'
            STOP
          END SELECT
        END IF

        ! Remove temp matrix
        CALL lsyssc_releaseMatrix(rtempMatrix)
    
      CASE DEFAULT
        PRINT *,'lsyssc_sortMatrix: Unsupported matrix format!'
        STOP
        
      END SELECT

    END SUBROUTINE
    
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat7_double (Da, DaH, Icol, IcolH, &
                                     Ild, IldH, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! Double precision version.
    !
    ! Storage technique 7 version.
  !</description>
    
  !<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Source matrix
    REAL(DP), DIMENSION(:), INTENT(IN) :: DaH
    
    ! Column structure of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: IcolH
    
    ! Row positions of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(neq+1), INTENT(IN) :: IldH
    
    ! Permutation of 1..neq describing the sorting
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing the sorting back
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<output>

    ! Destination matrix
    REAL(DP), DIMENSION(:), INTENT(OUT) :: Da

    ! Column structure of destination matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Icol
    
    ! Row positions of destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Ild
    
  !</output>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx

    ! Temporary variables for saving data of each line
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Ih1, Ih2
  
    ! Get memory for Ih1 and Ih2:
    ALLOCATE(Ih1(neq),Ih2(neq))
  
    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    DO i=1, neq

      ! Which row should be moved to here?
      iidx = Itr1(i)

      ! Copy the diagonal element.
      ! The new row starts at position ildIdx. Save this to Ild.
      Da(ildIdx) = DaH(IldH(iidx))
      Ild(i) = ildIdx
      Icol(ildIdx) = i
      ildIdx = ildIdx + 1
      
      ! Get the start- and end-index of the row that should be moved
      ! to here.
      ih1Idx = IldH(iidx)+1
      ih2Idx = IldH(iidx+1)-1

      ! The row has length...
      isize = ih2Idx - ih1Idx + 1

      ! Loop through the source row.
      ! IcolH defines the column numbers in the source row.
      ! Itr2 is the inverse permutation. Therefore, IcolH(Itr2) maps the
      ! column numbers in the source matrix to the column numbers in the
      ! destination matrix.
      ! Keeping that in mind, we set up:
      !  Ih1 := column numbers in the destination matrix.
      !  Ih2 := positions of the entries in the current row in the source matrix
      DO j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      END DO

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      CALL lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j=1, isize
        ! Get the matrix entry of the current column and write it to DA
        Da(ildIdx) = DaH(Ih2(j))

        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j) !Achtung: Ist dies so richtig ? Ih1(j)->Ih2(j) ?

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO

    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat7_single (Fa, FaH, Icol, IcolH, &
                                     Ild, IldH, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! Single precision version.
    !
    ! Storage technique 7 version.
  !</description>
    
  !<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Source matrix
    REAL(SP), DIMENSION(:), INTENT(IN) :: FaH
    
    ! Column structure of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: IcolH
    
    ! Row positions of source matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(IN) :: IldH
    
    ! Permutation of 1..neq describing the sorting
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing the sorting back
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<output>

    ! Destination matrix
    REAL(SP), DIMENSION(:), INTENT(OUT) :: Fa

    ! Column structure of destination matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Icol
    
    ! Row positions of destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Ild
    
  !</output>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx
    
    ! Temporary variables for saving data of each line
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Ih1, Ih2
  
    ! Get memory for Ih1 and Ih2:
    ALLOCATE(Ih1(neq),Ih2(neq))
  
    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    DO i=1, neq
      ! Which row should be moved to here?
      iidx = Itr1(i)

      ! Copy the diagonal element.
      ! The new row starts at position ildIdx. Save this to Ild.
      Fa(ildIdx) = FaH(IldH(iidx))
      Ild(i) = ildIdx
      Icol(ildIdx) = i
      ildIdx = ildIdx + 1
      
      ! Get the start- and end-index of the row that should be moved
      ! to here.
      ih1Idx = IldH(iidx)+1
      ih2Idx = IldH(iidx+1)-1

      ! The row has length...
      isize = ih2Idx - ih1Idx + 1

      ! Loop through the source row.
      ! IcolH defines the column numbers in the source row.
      ! Itr2 is the inverse permutation. Therefore, IcolH(Itr2) maps the
      ! column numbers in the source matrix to the column numbers in the
      ! destination matrix.
      ! Keeping that in mind, we set up:
      !  Ih1 := column numbers in the destination matrix.
      !  Ih2 := positions of the entries in the current row in the source matrix
      DO j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      END DO

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      CALL lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j=1, isize
        ! Get the matrix entry of the current column and write it to DA
        Fa(ildIdx) = FaH(Ih2(j))

        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j) 

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO

    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)
    
    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat7Struc (Icol, IcolH, Ild, IldH, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the structure of the given matrix, corresponding to Itr1/Itr2.
    ! Only the structure of a given matrix is resorted, the routine 
    ! doesn't handle the entries.
    !
    ! Storage technique 7 version. 
  !</description>
    
  !<input>

    ! Number of equations
    INTEGER(PREC_VECIDX) , INTENT(IN) :: neq

    ! Permutation of 1..neq describing how to resort
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing how to sort
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<inputoutput>

    ! Column description of matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Icol

    ! Row description of matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(INOUT) :: Ild
    
    ! Column structure of source matrix -> resorted matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: IcolH
    
    ! Row positions of source matrix -> resorted matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(INOUT) :: IldH
    
  !</inputoutput>
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx
    ! Temporary variables for saving data of each line
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Ih1, Ih2
  
    ! Get memory for Ih1 and Ih2:
    ALLOCATE(Ih1(neq),Ih2(neq))
    
    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    DO i=1, neq
      ! Which row should be moved to here?
      iidx = Itr1(i)

      ! Copy the diagonal element.
      ! The new row starts at position ildIdx. Save this to Ild.
      Ild(i) = ildIdx
      Icol(ildIdx) = i
      ildIdx = ildIdx+1

      ! Get the start- and end-index of the row that should be moved
      ! to here.
      ih1Idx = IldH(iidx)+1
      ih2Idx = IldH(iidx+1)-1

      ! The row has length...
      isize = ih2Idx-ih1Idx+1

      ! Loop through the source row.
      ! IcolH defines the column numbers in the source row.
      ! Itr2 is the inverse permutation. Therefore, IcolH(Itr2) maps the
      ! column numbers in the source matrix to the column numbers in the
      ! destination matrix.
      ! Keeping that in mind, we set up:
      !  Ih1 := column numbers in the destination matrix.
      !  Ih2 := positions of the entries in the current row in the source matrix
      DO j=ih1Idx, Ih2Idx
        Ih1(j-ih1Idx+1)=Itr2(IcolH(j))
        Ih2(j-ih1Idx+1)=j
      END DO

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      CALL lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j=1, isize
        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j)

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx+1
      END DO
    END DO
    Ild(neq+1) = IldH(neq+1) 
    
    DEALLOCATE(Ih1,Ih2)
  
  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat9_double (Da, DaH, Icol, IcolH, &
                                     Ild, IldH, Idiag, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! Double precision version.
    !
    ! Storage technique 9 version.
  !</description>
    
  !<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Source matrix
    REAL(DP), DIMENSION(:), INTENT(IN) :: DaH
    
    ! Column structure of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: IcolH
    
    ! Row positions of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(neq+1), INTENT(IN) :: IldH

    ! Permutation of 1..neq describing the sorting
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing the inverse permutation
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<output>

    ! Destination matrix
    REAL(DP), DIMENSION(:), INTENT(OUT) :: Da

    ! Column structure of destination matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Icol
    
    ! Row positions of destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Ild
    
    ! Positions of diagonal elements in the destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Idiag
    
  !</output>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Ih1, Ih2
  
    ! Get memory for Ih1 and Ih2:
    ALLOCATE(Ih1(neq),Ih2(neq))
  
    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1
    
    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    DO i=1, neq
    
      ! Which row should be moved to here?
      iidx = Itr1(i)
      
      ! The new row starts at position ildIdx. Save this to Ild.
      Ild(i) = ildIdx
      
      ! Get the start- and end-index of the row that should be moved
      ! to here.
      ih1Idx = IldH(iidx)
      ih2Idx = IldH(iidx+1)-1
      
      ! The row has length...
      isize = ih2Idx - ih1Idx + 1
      
      ! Loop through the source row.
      ! IcolH defines the column numbers in the source row.
      ! Itr2 is the inverse permutation. Therefore, IcolH(Itr2) maps the
      ! column numbers in the source matrix to the column numbers in the
      ! destination matrix.
      ! Keeping that in mind, we set up:
      !  Ih1 := column numbers in the destination matrix.
      !  Ih2 := positions of the entries in the current row in the source matrix
      DO j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      END DO
      
      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      CALL lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      DO idiagidx=1,isize
        IF (Ih1(idiagidx) .GE. i) EXIT
      END DO
      
      ! Write the position of the diagonal element to Idiag
      Idiag(i) = ildIdx + idiagidx-1
      
      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j = 1,isize
        
        ! Get the matrix entry of the current column and write it to DA
        Da(ildIdx) = DaH(Ih2(j))
        
        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j) 
        
        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO
    
    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat9_single (Fa, FaH, Icol, IcolH, &
                                     Ild, IldH, Idiag, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! Single precision version.
    !
    ! Storage technique 9 version.
  !</description>
    
  !<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Source matrix
    REAL(SP), DIMENSION(:), INTENT(IN) :: FaH
    
    ! Column structure of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: IcolH
    
    ! Row positions of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(neq+1), INTENT(IN) :: IldH

    ! Permutation of 1..neq describing the sorting
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing the inverse permutation
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<output>

    ! Destination matrix
    REAL(SP), DIMENSION(:), INTENT(OUT) :: Fa

    ! Column structure of destination matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Icol
    
    ! Row positions of destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Ild
    
    ! Positions of diagonal elements in the destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Idiag
    
  !</output>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Ih1, Ih2
  
    ! Get memory for Ih1 and Ih2:
    ALLOCATE(Ih1(neq),Ih2(neq))
  
    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1
    
    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    DO i=1, neq
    
      ! Which row should be moved to here?
      iidx = Itr1(i)
      
      ! The new row starts at position ildIdx. Save this to Ild.
      Ild(i) = ildIdx
      
      ! Get the start- and end-index of the row that should be moved
      ! to here.
      ih1Idx = IldH(iidx)
      ih2Idx = IldH(iidx+1)-1
      
      ! The row has length...
      isize = ih2Idx - ih1Idx + 1
      
      ! Loop through the source row.
      ! IcolH defines the column numbers in the source row.
      ! Itr2 is the inverse permutation. Therefore, IcolH(Itr2) maps the
      ! column numbers in the source matrix to the column numbers in the
      ! destination matrix.
      ! Keeping that in mind, we set up:
      !  Ih1 := column numbers in the destination matrix.
      !  Ih2 := positions of the entries in the current row in the source matrix
      DO j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      END DO
      
      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      CALL lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      DO idiagidx=1,isize
        IF (Ih1(idiagidx) .GE. i) EXIT
      END DO
      
      ! Write the position of the diagonal element to Idiag
      Idiag(i) = ildIdx + idiagidx-1
      
      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j = 1,isize
        
        ! Get the matrix entry of the current column and write it to DA
        Fa(ildIdx) = FaH(Ih2(j))
        
        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j) 
        
        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO
    
    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat9Struc (Icol, IcolH, &
                                   Ild, IldH, Idiag, Itr1, Itr2, neq)
  
  !<description>
    ! Resorts the structure of the given matrix, corresponding to Itr1/Itr2.
    ! Only the structure of a given matrix is resorted, the routine 
    ! doesn't handle the entries.
    !
    ! Storage technique 9 version. 
  !</description>
    
  !<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Column structure of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: IcolH
    
    ! Row positions of source matrix
    INTEGER(PREC_VECIDX), DIMENSION(neq+1), INTENT(IN) :: IldH

    ! Permutation of 1..neq describing the sorting
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr1
    
    ! Permutation of 1..neq describing the inverse permutation
    INTEGER(PREC_VECIDX), DIMENSION(neq), INTENT(IN) :: Itr2

  !</input>
    
  !<output>

    ! Column structure of destination matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: Icol
    
    ! Row positions of destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Ild
    
    ! Positions of diagonal elements in the destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(OUT) :: Idiag
    
  !</output>
    
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    INTEGER(I32), DIMENSION(:), ALLOCATABLE :: Ih1, Ih2
  
    ! Get memory for Ih1 and Ih2:
    ALLOCATE(Ih1(neq),Ih2(neq))
  
    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1
    
    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    DO i=1, neq
    
      ! Which row should be moved to here?
      iidx = Itr1(i)
      
      ! The new row starts at position ildIdx. Save this to Ild.
      Ild(i) = ildIdx
      
      ! Get the start- and end-index of the row that should be moved
      ! to here.
      ih1Idx = IldH(iidx)
      ih2Idx = IldH(iidx+1)-1
      
      ! The row has length...
      isize = ih2Idx - ih1Idx + 1
      
      ! Loop through the source row.
      ! IcolH defines the column numbers in the source row.
      ! Itr2 is the inverse permutation. Therefore, IcolH(Itr2) maps the
      ! column numbers in the source matrix to the column numbers in the
      ! destination matrix.
      ! Keeping that in mind, we set up:
      !  Ih1 := column numbers in the destination matrix.
      !  Ih2 := positions of the entries in the current row in the source matrix
      DO j=ih1Idx, ih2Idx
        Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
        Ih2(j-ih1Idx+1) = j
      END DO
      
      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      CALL lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      DO idiagidx=1,isize
        IF (Ih1(idiagidx) .GE. i) EXIT
      END DO
      
      ! Write the position of the diagonal element to Idiag
      Idiag(i) = ildIdx + idiagidx-1
      
      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j = 1,isize
        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j) 
        
        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO
    
    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortCR (Ih1, Ih2)
  
  !<description>
    ! Performs bubble sort for the vector Ih1 and does the
    ! same swaps also on vector Ih2.
  !</description>
    
  !<inputoutput>
      
    ! column vector
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Ih1
      
    ! row vector; same size as Ih1
    INTEGER(PREC_VECIDX), DIMENSION(SIZE(Ih1)), INTENT(INOUT) :: Ih2
      
  !</inputoutput>
!</subroutine>
    
    !local variables
    INTEGER(PREC_VECIDX) :: iaux, icomp
    LOGICAL :: bmore
    
    bmore = .TRUE.
    ! While there are swaps necessary...
    DO WHILE (bmore)
      bmore = .FALSE.
      ! Test for all elements minus the last one
      DO icomp=1, SIZE(Ih1)-1
        ! If the order is correct, next entry; othewise swap entries
        IF (Ih1(icomp) .GT. Ih1(icomp+1)) THEN
          iaux = Ih1(icomp)
          Ih1(icomp) = Ih1(icomp+1)
          Ih1(icomp+1) = iaux
          iaux = Ih2(icomp)
          Ih2(icomp) = Ih2(icomp+1)
          Ih2(icomp+1) = iaux
          bmore = .TRUE.
        ENDIF
      END DO  
    END DO
  
  END SUBROUTINE 
  
  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_addIndex (h_Ix,ivalue,istartpos,ilength)
  
!<description>

  ! This is an auxiliary routine. It accepts a handle to an integer array
  ! and adds the value ivalue to each entry of this vector.
  ! This can be used e.g. to convert the index vector of a matrix
  ! von 1-based to 0-based for an external library (UMFPACK).

!</description>
  
!<input>
  ! Handle to an integer data array.
  INTEGER, INTENT(IN) :: h_Ix

  ! The value to add to every entry.
  INTEGER(I32) :: ivalue

  ! OPTIONAL: Starting position of a part of the vector to modify; usually = 1
  INTEGER(I32), OPTIONAL :: istartpos

  ! OPTIONAL: Length of the part of the vector to modify; <= SIZE(Ix)
  INTEGER(I32), OPTIONAL :: ilength

!</input>

!</subroutine>
    
    ! Actual length
    INTEGER(I32) :: iactlength,istart,i
    
    ! Vector
    INTEGER(I32), DIMENSION(:), POINTER :: p_Ix
    
    ! Get the array
    CALL storage_getbase_int (h_Ix,p_Ix)
    IF (PRESENT(istartpos)) THEN
      istart = MIN(istartpos,SIZE(p_Ix)+1)
    ELSE
      istart = 1
    END IF
    IF (PRESENT(ilength)) THEN
      iactlength = MIN(ilength,SIZE(p_Ix)-istart+1)
    ELSE
      iactlength = SIZE(p_Ix)-istart+1
    END IF
    p_Ix => p_Ix(istart:istart + iactlength - 1)
    
    ! Increase by ivalue
    DO i=1,iactlength
      p_Ix(i) = p_Ix(i) + ivalue
    END DO
    
  END SUBROUTINE
  
  !****************************************************************************
!<subroutine>
  
  REAL(DP) FUNCTION lsyssc_vectorNorm (rx,cnorm,iposMax)
  
!<description>
  ! Calculates the norm of a vector. cnorm specifies the type of norm to
  ! calculate.
!</description>
  
!<input>
  ! Vector to calculate the norm of.
  TYPE(t_vectorScalar), INTENT(IN)                  :: rx

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
  ! The scalar product (rx,ry) of the two block vectors.
!</result>

!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata

  ! Is there data at all?
  IF (rx%h_Ddata .EQ. ST_NOHANDLE) THEN
    PRINT *,'Error in lsyssc_vectorNorm: Vector empty!'
    STOP
  END IF
  
  ! Take care of the data type before doing a scalar product!
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the array and calculate the norm
    CALL lsyssc_getbase_double (rx,p_Ddata)
    lsyssc_vectorNorm = lalg_normDble (p_Ddata,cnorm,iposMax) 
    
  CASE (ST_SINGLE)
    ! Get the array and calculate the norm
    CALL lsyssc_getbase_single (rx,p_Fdata)
    lsyssc_vectorNorm = lalg_normSngl (p_Fdata,cnorm,iposMax) 
    
  CASE DEFAULT
    PRINT *,'lsyssc_vectorNorm: Unsupported data type!'
    STOP
  END SELECT
  
  END FUNCTION

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_invertedDiagMatVec (rmatrix,rvectorSrc,dscale,rvectorDst)
  
!<description>
  ! This routine multiplies the weighted inverted diagonal domega*D^{-1}
  ! of the matrix rmatrix with the vector rvectorSrc and stores the result 
  ! into the vector rvectorDst:
  !   rvectorDst = dscale * D^{-1} * rvectorSrc
  ! Both, rvectorSrc and rvectorDst may coincide.
!</description>
  
!<input>
  ! The matrix. 
  TYPE(t_matrixScalar), INTENT(IN) :: rmatrix

  ! The source vector.
  TYPE(t_vectorScalar), INTENT(IN) :: rvectorSrc

  ! A multiplication factor. Standard value is 1.0_DP
  REAL(DP), INTENT(IN) :: dscale
!</input>

!<inputoutput>
  ! The destination vector which receives the result.
  TYPE(t_vectorScalar), INTENT(INOUT) :: rvectorDst
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: h_Idiag
  INTEGER(PREC_VECIDX) :: i
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiag
  REAL(DP), DIMENSION(:), POINTER :: p_Da, p_Dvec, p_Dvec2
  REAL(SP), DIMENSION(:), POINTER :: p_Fa, p_Fvec, p_Fvec2
  REAL(DP) :: dmyscale
  REAL(SP) :: fmyscale

  ! Let's hope, the matrix and the vectors are compatible.
  ! As we multiply with D^{-1}, this forces the matrix to be quadratic 
  ! and both vectors to be of the same length!
  CALL lsyssc_isVectorCompatible (rvectorSrc,rvectorDst)
  CALL lsyssc_isMatrixCompatible (rvectorSrc,rmatrix,.FALSE.)

  IF (rvectorSrc%cdataType .NE. rvectorDst%cdataType) THEN
    PRINT *,'lsyssc_invertedDiagMatVec: Vectors have different precisions!'
    STOP
  END IF

  ! Which matrix structure do we have?
  SELECT CASE (rmatrix%cmatrixFormat)
  CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
    ! Get a pointer to the diagonal - either h_KLD or h_Kdiagonal
    IF (rmatrix%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
      h_Idiag = rmatrix%h_Kdiagonal
    ELSE
      h_Idiag = rmatrix%h_Kld
    END IF
    CALL storage_getbase_int(h_Idiag,p_Kdiag)
    
    ! Data type?
    SELECT CASE (rmatrix%cdataType)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double (rmatrix%h_Da,p_Da)
      dmyscale = dscale * rmatrix%dscaleFactor

      ! And the vector(s)?
      SELECT CASE (rvectorSrc%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rvectorSrc,p_Dvec)
        CALL lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let's go...
        DO i=1,rvectorSrc%NEQ
          p_Dvec2(i) = p_Dvec(i)*dmyscale/p_Da(p_Kdiag(i))
        END DO
        
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rvectorSrc,p_Fvec)
        CALL lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let's go...
        DO i=1,rvectorSrc%NEQ
          p_Fvec2(i) = p_Fvec(i)*dmyscale/p_Da(p_Kdiag(i))
        END DO
        
      CASE DEFAULT
        PRINT *,'lsyssc_invertedDiagMatVec: unsupported vector precision!'
        STOP
      END SELECT
    CASE (ST_SINGLE)
      CALL storage_getbase_single (rmatrix%h_Da,p_Fa)
      fmyscale = REAL(dscale * rmatrix%dscaleFactor,SP)

      ! And the vector(s)?
      SELECT CASE (rvectorSrc%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rvectorSrc,p_Dvec)
        CALL lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let's go...
        DO i=1,rvectorSrc%NEQ
          p_Dvec2(i) = p_Dvec(i)*fmyscale/p_Fa(p_Kdiag(i))
        END DO
        
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rvectorSrc,p_Fvec)
        CALL lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let's go...
        DO i=1,rvectorSrc%NEQ
          p_Fvec2(i) = p_Fvec(i)*fmyscale/p_Fa(p_Kdiag(i))
        END DO
        
      CASE DEFAULT
        PRINT *,'lsyssc_invertedDiagMatVec: unsupported vector precision!'
        STOP
      END SELECT
    CASE DEFAULT
      PRINT *,'lsyssc_invertedDiagMatVec: unsupported matrix precision!'
      STOP
    END SELECT
    
  CASE DEFAULT
    PRINT *,'lsyssc_invertedDiagMatVec: unsupported matrix format!'
    STOP
  END SELECT
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_copyVector (rx,ry,bforceOnlyData)
  
!<description>
  ! Copies vector data: ry = rx.
  ! Both vectors must have the same size. All structural data of rx is
  ! transferred to ry, so rx and ry are compatible to each other afterwards.
!</description>

!<input>
  ! Source vector
  TYPE(t_vectorScalar),INTENT(IN) :: rx
  
  ! OPTIONAL: Copy only data.
  ! If set to TRUE, this forces the routine to copy only the vector data of rx to ry.
  ! Structural data is not copied.
  LOGICAL, INTENT(IN), OPTIONAL :: bforceOnlyData
!</input>

!<inputoutput>
  ! Destination vector
  TYPE(t_vectorScalar),INTENT(INOUT) :: ry
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: h_Ddata, cdataType
  LOGICAL :: bisCopy,bonlyData
  INTEGER(PREC_VECIDX) :: ioffset
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource,p_Ddest
  REAL(SP), DIMENSION(:), POINTER :: p_Fsource,p_Fdest
  
  IF (PRESENT(bforceOnlyData)) THEN 
    bonlyData = bforceOnlyData
  ELSE
    bonlyData = .FALSE.
  END IF
  
  IF (.NOT. bonlyData) THEN
    ! Standard case: Copy the whole vector information
    !
    ! First, make a backup of some crucial data so that it does not
    ! get destroyed.
    h_Ddata = ry%h_Ddata
    cdataType = ry%cdataType
    bisCopy = ry%bisCopy
    ioffset = rx%iidxFirstEntry
    
    ! Then transfer all structural information of rx to ry.
    ! This automatically makes both vectors compatible to each other.
    ry = rx
    
    ! Restore crucial data
    ry%h_Ddata = h_Ddata
    ry%bisCopy = bisCopy
    ry%iidxFirstEntry = ioffset 
  END IF
  
  IF (rx%cdataType .NE. ry%cdataType) THEN
    PRINT *,'lsyssc_copyVector: different data types not supported!'
    STOP
  END IF
  
  ! And finally copy the data. 
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    CALL lsyssc_getbase_double (rx,p_Dsource)
    CALL lsyssc_getbase_double (ry,p_Ddest)
    CALL lalg_copyVectorDble (p_Dsource,p_Ddest)
    
  CASE (ST_SINGLE)
    CALL lsyssc_getbase_single (rx,p_Fsource)
    CALL lsyssc_getbase_single (ry,p_Fdest)
    CALL lalg_copyVectorSngl (p_Fsource,p_Fdest)

  CASE DEFAULT
    PRINT *,'lsyssc_copyVector: Unsupported data type!'
    STOP
  END SELECT
   
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_scaleVector (rx,c)
  
!<description>
  ! Scales a vector vector rx: rx = c * rx
!</description>
  
!<inputoutput>
  ! Source and destination vector
  TYPE(t_vectorScalar), INTENT(INOUT) :: rx
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
    CALL lsyssc_getbase_double(rx,p_Ddata)
    CALL lalg_scaleVectorDble (p_Ddata,c)  

  CASE (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    CALL lsyssc_getbase_single(rx,p_Fdata)
    CALL lalg_scaleVectorSngl (p_Fdata,REAL(c,SP))  

  CASE DEFAULT
    PRINT *,'lsyssc_scaleVector: Unsupported data type!'
    STOP
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_clearVector (rx)
  
!<description>
  ! Clears the block vector dx: Dx = 0
!</description>
  
!<inputoutput>
  ! Destination vector to be cleared
  TYPE(t_vectorScalar), INTENT(INOUT) :: rx
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource
  REAL(SP), DIMENSION(:), POINTER :: p_Ssource
  
  ! Take care of the data type
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the pointer and scale the whole data array.
    CALL lsyssc_getbase_double(rx,p_Dsource)
    CALL lalg_clearVectorDble (p_Dsource)
  
  CASE (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    CALL lsyssc_getbase_single(rx,p_Ssource)
    CALL lalg_clearVectorSngl (p_Ssource)

  CASE DEFAULT
    PRINT *,'lsyssc_clearVector: Unsupported data type!'
    STOP
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_vectorLinearComb (rx,ry,cx,cy)
  
!<description>
  ! Performs a linear combination: ry = cx * rx  +  cy * ry
  ! Both vectors must be compatible to each other (same size, sorting 
  ! strategy,...).
!</description>  
  
!<input>
  ! First source vector
  TYPE(t_vectorScalar), INTENT(IN)   :: rx
  
  ! Scaling factor for Dx
  REAL(DP), INTENT(IN)               :: cx

  ! Scaling factor for Dy
  REAL(DP), INTENT(IN)               :: cy
!</input>

!<inputoutput>
  ! Second source vector; also receives the result
  TYPE(t_vectorScalar), INTENT(INOUT) :: ry
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource, p_Ddest
  REAL(SP), DIMENSION(:), POINTER :: p_Ssource, p_Sdest
  
  ! The vectors must be compatible to each other.
  CALL lsyssc_isVectorCompatible (rx,ry)

  IF (rx%cdataType .NE. ry%cdataType) THEN
    PRINT *,'lsyssc_vectorLinearComb: different data types not supported!'
    STOP
  END IF
  
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the pointers and copy the whole data array.
    CALL lsyssc_getbase_double(rx,p_Dsource)
    CALL lsyssc_getbase_double(ry,p_Ddest)
    
    CALL lalg_vectorLinearCombDble (p_Dsource,p_Ddest,cx,cy)

  CASE (ST_SINGLE)
    ! Get the pointers and copy the whole data array.
    CALL lsyssc_getbase_single(rx,p_Ssource)
    CALL lsyssc_getbase_single(ry,p_Sdest)
    
    CALL lalg_vectorLinearCombSngl (p_Ssource,p_Sdest,cx,cy)
  
  CASE DEFAULT
    PRINT *,'lsyssc_vectorLinearComb: Unsupported data type!'
    STOP
  END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_copyMatrix (rsourceMatrix,rdestMatrix)
  
!<description>
  ! Copies a matrix data: rsourceMatrix = rdestMatrix.
  ! Both matrices must have the same size and the same structure.
  ! All structural (discretisation related) data of rsourceMatrix is
  ! transferred to rrdestMatrix.
!</description>
  
!<input>
  ! Source vector
  TYPE(t_matrixScalar),INTENT(IN) :: rsourceMatrix
!</input>

!<inputoutput>
  ! Destination vector
  TYPE(t_matrixScalar),INTENT(INOUT) :: rdestMatrix
!</inputoutput>
  
!</subroutine>

  ! local variables
  TYPE(t_matrixScalar) :: rtempMatrix
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource,p_Ddest
  REAL(SP), DIMENSION(:), POINTER :: p_Fsource,p_Fdest
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolSource,p_KcolDest
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldSource,p_KldDest
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalSource,p_KdiagonalDest
  
  ! Is it possible at all to copy the matrix? Both matrices must have
  ! the same size, otherwise the memory does not fit.
  IF (rsourceMatrix%cmatrixFormat .NE. rdestMatrix%cmatrixFormat) THEN
    PRINT *,'lsyssc_copyMatrix: Different matrix formats not allowed!'
    STOP
  END IF
  
  IF (rsourceMatrix%NA .NE. rdestMatrix%NA) THEN
    PRINT *,'lsyssc_copyMatrix: Matrices have different size!'
    STOP
  END IF
  
  ! NEQ/NCOLS is irrelevant. It's only important that we have enough memory
  ! and the same structure!
  
  IF (rsourceMatrix%cdataType .NE. rdestMatrix%cdataType) THEN
    PRINT *,'lsyssc_copyMatrix: Matrices have different data types!'
    STOP
  END IF
  
  ! First, make a backup of the matrix for restoring some cricial data.
  rtempMatrix = rdestMatrix
  
  ! Copy the matrix structure
  rdestMatrix = rsourceMatrix
  
  ! Restore crucial data
  rdestMatrix%imatrixSpec = rtempMatrix%imatrixSpec

  ! Which structure do we actually have?
  SELECT CASE (rsourceMatrix%cmatrixFormat)
  CASE (LSYSSC_MATRIX9)
  
    ! Restore crucial format-specific data
    rdestMatrix%h_DA        = rtempMatrix%h_DA
    rdestMatrix%h_Kcol      = rtempMatrix%h_Kcol
    rdestMatrix%h_Kld       = rtempMatrix%h_Kld
    rdestMatrix%h_Kdiagonal = rtempMatrix%h_Kdiagonal
    
    ! And finally copy the data. 
    SELECT CASE (rsourceMatrix%cdataType)
    
    CASE (ST_DOUBLE)
      CALL storage_getbase_double (rsourceMatrix%h_DA,p_Dsource)
      CALL storage_getbase_double (rdestMatrix%h_DA,p_Ddest)
      CALL lalg_copyVectorDble (p_Dsource,p_Ddest)

    CASE (ST_SINGLE)
      CALL storage_getbase_single (rsourceMatrix%h_DA,p_Fsource)
      CALL storage_getbase_single (rdestMatrix%h_DA,p_Fdest)
      CALL lalg_copyVectorSngl (p_Fsource,p_Fdest)

    CASE DEFAULT
      PRINT *,'lsyssc_copyMatrix: Unsupported data type!'
      STOP
    END SELECT
    
    CALL storage_getbase_int (rsourceMatrix%h_Kcol,p_KcolSource)
    CALL storage_getbase_int (rdestMatrix%h_Kcol,p_KcolDest)
    CALL lalg_copyVectorInt (p_KcolSource,p_KcolDest)

    CALL storage_getbase_int (rsourceMatrix%h_Kld,p_KldSource)
    CALL storage_getbase_int (rdestMatrix%h_Kld,p_KldDest)
    CALL lalg_copyVectorInt (p_KldSource,p_KldDest)
    
    CALL storage_getbase_int (rsourceMatrix%h_Kdiagonal,p_KdiagonalSource)
    CALL storage_getbase_int (rdestMatrix%h_Kdiagonal,p_KdiagonalDest)
    CALL lalg_copyVectorInt (p_KdiagonalSource,p_KdiagonalDest)
      
  CASE (LSYSSC_MATRIX7)

    ! Restore crucial format-specific data
    rdestMatrix%h_DA        = rtempMatrix%h_DA
    rdestMatrix%h_Kcol      = rtempMatrix%h_Kcol
    rdestMatrix%h_Kld       = rtempMatrix%h_Kld
    
    ! And finally copy the data. 
    SELECT CASE (rsourceMatrix%cdataType)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double (rsourceMatrix%h_DA,p_Dsource)
      CALL storage_getbase_double (rdestMatrix%h_DA,p_Ddest)
      CALL lalg_copyVectorDble (p_Dsource,p_Ddest)

    CASE (ST_SINGLE)
      CALL storage_getbase_single (rsourceMatrix%h_DA,p_Fsource)
      CALL storage_getbase_single (rdestMatrix%h_DA,p_Fdest)
      CALL lalg_copyVectorSngl (p_Fsource,p_Fdest)

    CASE DEFAULT
      PRINT *,'storage_copyMatrix: Unsupported data type!'
      STOP
    END SELECT

    CALL storage_getbase_int (rsourceMatrix%h_Kcol,p_KcolSource)
    CALL storage_getbase_int (rdestMatrix%h_Kcol,p_KcolDest)
    CALL lalg_copyVectorInt (p_KcolSource,p_KcolDest)

    CALL storage_getbase_int (rsourceMatrix%h_Kld,p_KldSource)
    CALL storage_getbase_int (rdestMatrix%h_Kld,p_KldDest)
    CALL lalg_copyVectorInt (p_KldSource,p_KldDest)
      
  CASE DEFAULT
    PRINT *,'lsyssc_copyMatrix: Unsupported matrix format!'
    STOP
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>
  SUBROUTINE lsyssc_transpMatStruct79double (nrow, ncol, Icol, Irow, &
                                             Itmp, IcolDest, IrowDest)
  
  !<description>
    ! This routine accepts the structure of a structure-7 or structure-9
    ! matrix and creates the structure of the transposed matrix from it.
    !
    ! The resulting structure is of structure 9, i.e. not pivoted anymore 
    ! with the diagonal entry in front if a structure-7 matrix is to be 
    ! transposed!
  !</description>
    
  !<input>
    ! Number of rows in the source matrix
    INTEGER(PREC_VECIDX), INTENT(IN) :: nrow
    
    ! Number of columns in the source matrix
    INTEGER(PREC_VECIDX), INTENT(IN) :: ncol
    
    ! Column structure of the source matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Icol
    
    ! Row structure of the source matrix
    INTEGER(PREC_MATIDX), DIMENSION(nrow+1), INTENT(IN) :: Irow
  !</input>
    
  !<output>
    ! Column structure of the destination matrix
    ! The array must be of the same size as Icol!
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: IcolDest

    ! Row structure of the destination matrix
    INTEGER(PREC_MATIDX), DIMENSION(ncol+1), INTENT(OUT) :: IrowDest
  !</output>
    
  !<inputoutput>
    ! Auxiliary array of size ncol
    INTEGER(PREC_VECIDX), DIMENSION(ncol), INTENT(INOUT) :: Itmp
  !</inputoutput>
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, isize, icolumn, ncolumn
    
    ! determin the number of matrix entries
    isize = Irow(nrow+1)-1
    
    ! clear row struxture of the destination matrix
    DO i=1, ncol+1
      IrowDest(i) = 0
    END DO
    
    ! Count how many entries <> 0 are in each column. Note this into
    ! the IrowDest array shifted by 1.
    DO i=1, isize
      IrowDest(Icol(i)+1) = IrowDest(Icol(i)+1)+1
    END DO
    
    ! Now build the final IrowDest by adding up the IrowDest entries.
    IrowDest(1) = 1
    DO i=2, ncol+1
      IrowDest(i) = IrowDest(i)+IrowDest(i-1)
    END DO
    
    ! Now IcolDest must be created. This requires another loop trough
    ! the matrix structure. Itmp receives the index of how many entries
    ! have been written to each row.
    
    ! clear auxiliary vector
    DO i=1, ncol
      Itmp(i) = 0
    END DO

    DO i=1, nrow
      ncolumn = Irow(i+1)-Irow(i)
      DO j=1, ncolumn
        ! Get the column of the item in question -> new row number.
        icolumn = Icol(Irow(i)+j-1)
        ! Rows get columns by transposing, therefore note i as column
        ! number in IcolDest
        IcolDest(IrowDest(icolumn)+Itmp(icolumn)) = i
        ! Increment running index of that row
        Itmp(icolumn) = Itmp(icolumn)+1
      END DO
    END DO
    
  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE lsyssc_transpMat79double (nrow, ncol, Da, Icol, Irow, &
                                       Itmp, DaDest, IcolDest, IrowDest)
  
  !<description>
    ! Auxiliary routine.
    ! This routine accepts a structure-7 or structure-9
    ! matrix and creates the transposed matrix from it.
    !
    ! The resulting matrix is of format 9, i.e. not pivoted anymore with 
    ! the diagonal entry in front if a structure-7 matrix is to be transposed!
  !</description>
    
  !<input>
    ! Number of rows in the source matrix
    INTEGER(PREC_VECIDX), INTENT(IN) :: nrow
    
    ! Number of columns in the source matrix
    INTEGER(PREC_VECIDX), INTENT(IN) :: ncol
    
    ! The entries of the source matrix
    REAL(DP), DIMENSION(:), INTENT(IN) :: Da
    
    ! Column structure of the source matrix
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Icol
    
    ! Row structure of the source matrix.
    ! DIMENSION(nrow+1)
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Irow
  !</input>
    
  !<output>
    ! The entries of the destination matrix
    ! The array must be of the same size as Da or Icol
    REAL(DP), DIMENSION(:), INTENT(OUT) :: DaDest
    
    ! Column structure of the destination matrix
    ! The array must be of the same size as Icol!
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(OUT) :: IcolDest

    ! Row structure of the destination matrix.
    ! DIMENSION(ncol+1)
    INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(OUT) :: IrowDest
  !</output>
    
  !<inputoutput>
    ! Auxiliary array.
    ! DIMENSION(ncol)
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Itmp
  !</inputoutput>
  
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j, isize, icolumn, ncolumn
    
    IF ((SIZE(Itmp) .NE. ncol) .OR. (SIZE(IrowDest) .NE. ncol+1) .OR. &
        (SIZE(Irow) .NE. nrow+1)) THEN
      PRINT *,'lsyssc_transpMat79double: array parameters have wrong size!'
      STOP
    END IF
    
    ! determine the number of matrix entries
    isize = Irow(nrow+1)-1
    
    ! clear row structure of the destination matrix
    DO i=1, ncol+1
      IrowDest(i) = 0
    END DO
    
    ! Count how many entries <> 0 are in each column. Note this
    ! into the IrowDest array.
    ! This requires one loop through the matrix structure. The
    ! corresponding number is written into KLDD, shifted by 1
    ! (i.e. IrowDest(2) is the number of entries of the 1st column).
    ! This helps to create the real IrowDest more easily later.
    DO i=1, isize
      IrowDest(Icol(i)+1) = IrowDest(Icol(i)+1)+1
    END DO
    
    ! Now build the real IrowDest. This consists of indices, where 
    ! each row starts. Row 1 starts at position 1:
    IrowDest(1) = 1
    
    ! Adding the number of entries IrowDest(i) in the row to
    ! IrowDest(i-1) gives the new row pointer of row i of the transposed
    ! matrix.
    DO i=2, ncol+1
      IrowDest(i) = IrowDest(i)+IrowDest(i-1)
    END DO
    
    ! That's it for IrowDest. 
    ! Now IcolDest must be created. This requires another loop trough
    ! the matrix structure. Itmp receives the index of how many entries
    ! have been written to each row.
    
    ! clear auxiliary vector
    DO i=1, ncol
      Itmp(i) = 0
    END DO

    DO i=1, nrow
      ! Loop through the row
      ncolumn = Irow(i+1)-Irow(i)
      DO j=1, ncolumn
        ! Get the column of the item in question -> new row number.
        icolumn = Icol(Irow(i)+j-1)
        ! Rows get columns by transposing, therefore note i as column
        ! number in IcolDest
        IcolDest(IrowDest(icolumn)+Itmp(icolumn)) = i
        ! Copy the matrix entry:
        DaDest(IrowDest(icolumn)+Itmp(icolumn)) = Da(Irow(i)+j-1)
        ! Increment running index of that row
        Itmp(icolumn) = Itmp(icolumn)+1
      END DO
    END DO
    
  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_pivotiseMatrix7double (neq, Icol, Irow, Da, berror)
  
  !<description>
    ! Auxiliary routine.
    ! This routine repivots a matrix in structure 7: The diagonal element
    ! of each row is moved to the front.
  !</description>
    
  !<input>
    ! Number of rows/columns/equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
  !</input>
  
  !<output>
    ! If .TRUE. a diagonal element was not found
    LOGICAL, OPTIONAL, INTENT(OUT) :: berror
  !</output>
        
  !<inputoutput>
    ! OPTIONAL: Matrix entries.
    ! If not specified, only the matrix structure is pivotised.
    REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Da
    
    ! Matrix column structure
    INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Icol
    
    ! Matrix row structure
    INTEGER(PREC_MATIDX), DIMENSION(neq+1), INTENT(INOUT) :: Irow
  !</inputoutput>
!</subroutine>
    
    !local variables
    INTEGER(I32) :: i, j
    LOGICAL :: bnodiag
    
    bnodiag = .FALSE.
    
    IF (PRESENT(Da)) THEN
      ! Structure + entries
      !
      ! Loop through the rows
      DO i=1, neq
        ! Is the line already pivoted? Most unlikely, but we test for sure
        IF (Icol(Irow(i)) .NE. i) THEN

          ! Find the position of the diagonal element
          bnodiag = .TRUE.
          DO j=Irow(i), Irow(i+1)-1
            IF (j .EQ. i) THEN
              bnodiag = .FALSE.
              EXIT
            END IF
          END DO
          ! Oops, diagonal element not found - cancel
          IF (bnodiag) EXIT
          
          ! Ringshift the slice Icol(Irow(i):j) by 1 to the right.
          ! This puts the diagonal element to the front
          ! The same operation is also done for Da(Irow(i):j)
          Icol(Irow(i):j) = CSHIFT(Icol(Irow(i):j),-1)
          Da  (Irow(i):j) = CSHIFT(Da  (Irow(i):j),-1)
        END IF
      END DO

    ELSE
      ! Only structure
      !
      ! Loop through the rows
      DO i=1, neq
        ! Is the line already pivoted? Most unlikely, but we test for sure
        IF (Icol(Irow(i)) .NE. i) THEN

          ! Find the position of the diagonal element
          bnodiag = .TRUE.
          DO j=Irow(i), Irow(i+1)-1
            IF (j .EQ. i) THEN
              bnodiag = .FALSE.
              EXIT
            END IF
          END DO
          ! Oops, diagonal element not found - cancel
          IF (bnodiag) EXIT
          
          ! Ringshift the slice Icol(Irow(i):j) by 1 to the right.
          ! This puts the diagonal element to the front
          ! The same operation is also done for Da(Irow(i):j)
          Icol(Irow(i):j) = CSHIFT(Icol(Irow(i):j),-1)
          Da  (Irow(i):j) = CSHIFT(Da  (Irow(i):j),-1)
        END IF
      END DO
      
    END IF

    ! Report error status
    IF (PRESENT(berror)) berror = bnodiag    
  
  END SUBROUTINE 

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_transposeMatrix (rmatrix,rtransposedMatrix,itransFlag)
  
!<description>
  ! This routine transposes a matrix rmatrix and creates the transposed
  ! matrix in rtransposedMatrix.
  ! itransFlag decides (if specified) how the creation of the
  ! transposed matrix is to be performed.
  !
  ! Remark: Building a transposed matrix involves switching the trial- and
  !  test-functions of the underlying discretisation. This routine does
  !  *not* perform this switching! If the matrices base on a 
  !  discretisation, the application has to set up the corresponding
  !  discretisation structure manually and to attach it to the transposed
  !  matrix!
!</description>

!<input>
  ! OPTIONAL: transposed-flag; one of the LSYSSC_TR_xxxx constants:
  ! =LSYSSC_TR_ALL or not specified: transpose the matrix completely.
  !     If there is a matrix in rtransposedMatrix, it's overwritten 
  !       or an error is thrown, depending on whether the
  !       structural data (NA,...) of rtransposedMatrix matches rmatrix.
  !     If there's no matrix in rtransposedMatrix, a new matrix is
  !       created.
  ! =LSYSSC_TR_STRUCTURE           : Transpose only the matrix structure; 
  !     the content of the transposed matrix is invalid afterwards.
  !     If there is a matrix in rtransposedMatrix, the structure is 
  !       overwritten or an error is thrown, depending on whether the
  !       structural data (NA,...) of rtransposedMatrix matches rmatrix.
  !     If there's no matrix in rtransposedMatrix, a new matrix structure
  !       without entries is created.
  ! =LSYSSC_TR_VIRTUAL             : Actually don't touch the matrix 
  !     structure, but invert the 'transposed' flag in imatrixSpec. 
  !
  !     rtransposedMatrix is created by copying the data handles from
  !     rmatrix - any previous data in rtransposedMatrix is overwritten. 
  !     Afterwards, rtransposedMatrix is marked as transposed 
  !     without modifying the structures, and all matrix-vector operations
  !     will be performed with the transposed matrix. 
  !     But caution: there may be some algorithms that don't work with such 
  !       'virtually' transposed matrices!
  ! =LSYSSC_TR_VIRTUALCOPY         : The same as LSYSSC_TR_VIRTUAL, but
  !     creates a duplicate of the source matrix in memory, thus resulting
  !     in rtransposedMatrix being a totally independent matrix.
  INTEGER, INTENT(IN), OPTIONAL :: itransFlag
  
  ! The matrix to be transposed.
  TYPE(t_matrixScalar),INTENT(IN) :: rmatrix
!</input>

!<inputoutput>
  ! The transposed matrix. 
  ! If the structure is empty, a new matrix is created.
  ! If the structure exists and structural data matches rmatrix, 
  !   it's overwritten.
  ! If the structure exists and structural data does not match rmatrix, 
  !   an error is thrown.
  TYPE(t_matrixScalar),INTENT(INOUT) :: rtransposedMatrix
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER :: itrans, h_Itemp
    INTEGER(PREC_VECIDX) :: ntemp
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolSource,p_KcolDest
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldSource,p_KldDest
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiagonal
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Itemp
    REAL(DP), DIMENSION(:), POINTER :: p_DaSource,p_DaDest
    
    itrans = LSYSSC_TR_ALL
    IF (PRESENT(itransFlag)) itrans = itransFlag
    
    ! How to perform the creation of the transposed matrix?
    IF (itrans .EQ. LSYSSC_TR_VIRTUAL) THEN
    
      ! Duplicate the matrix structure, copy all handles.
      ! Don't allocate any memory.
      CALL lsyssc_duplicateMatrix (rmatrix,rtransposedMatrix,&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    
      ! Only change the 'transposed' flag in imatrixSpec
      rtransposedMatrix%imatrixSpec = &
        IEOR(rtransposedMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED)
      
      ! Exchange NEQ and NCOLS as these always describe the actual matrix.
      ntemp = rtransposedMatrix%NEQ
      rtransposedMatrix%NEQ = rtransposedMatrix%NCOLS
      rtransposedMatrix%NCOLS = ntemp
   
      ! That's it.
      RETURN
    END IF

    IF (itrans .EQ. LSYSSC_TR_VIRTUALCOPY) THEN
    
      ! Duplicate the complete matrix structure in memory
      CALL lsyssc_duplicateMatrix (rmatrix,rtransposedMatrix,&
                                    LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
    
      ! Only change the 'transposed' flag in imatrixSpec
      rtransposedMatrix%imatrixSpec = &
        IEOR(rtransposedMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED)
      
      ! Exchange NEQ and NCOLS as these always describe the actual matrix.
      ntemp = rtransposedMatrix%NEQ
      rtransposedMatrix%NEQ = rtransposedMatrix%NCOLS
      rtransposedMatrix%NCOLS = ntemp
   
      ! That's it.
      RETURN
    END IF

    ! Does the destination matrix exist?
    IF (rtransposedMatrix%NEQ .NE. 0) THEN
    
      ! Make sure the destination matrix can accept all data of the source
      ! matrix. Matrix format and size must match.
      IF (rMatrix%cmatrixFormat .NE. rtransposedMatrix%cmatrixFormat) THEN
        PRINT *,'lsyssc_transposeMatrix: Different matrix formats not allowed!'
        STOP
      END IF

      IF (rmatrix%NA .NE. rtransposedMatrix%NA) THEN
        PRINT *,'lsyssc_transposeMatrix: Matrices have different size!'
        STOP
      END IF

      ! NEQ/NCOLS is irrelevant. It's only important that we have enough memory
      ! and the same structure!

      IF (itrans .EQ. LSYSSC_TR_ALL) THEN
        IF (rmatrix%cdataType .NE. rtransposedMatrix%cdataType) THEN
          PRINT *,'lsyssc_transposeMatrix: Matrices have different data types!'
          STOP
        END IF
      END IF

    ELSE
      ! Otherwise, we must create a new matrix that accepts the transposed
      ! data. 
      ! Copy the matrix structure, reset all handles, so we have a 'template'
      ! matrix with the same data as the original one but without dynamic
      ! information attached.
      CALL lsyssc_duplicateMatrix (rmatrix,rtransposedMatrix,&
                                    LSYSSC_DUP_TEMPLATE,LSYSSC_DUP_TEMPLATE)
    END IF

    ! Otherwise, check that we are able to create the transposed matrix at all.
    IF (itrans .EQ. LSYSSC_TR_ALL) THEN
      IF (rmatrix%cdataType .NE. ST_DOUBLE) THEN
        PRINT *,'lsyssc_transposeMatrix: Unsupported data type!'
        STOP
      END IF
    END IF
    
    ! Now what should we do...
    SELECT CASE (itrans)
    CASE (LSYSSC_TR_STRUCTURE)
      ! Only the structure is to be transposed. Which structure do we have?
      SELECT CASE (rmatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9)
        ! Is the source matrix saved tranposed? In that case simply copy
        ! the data arrays, this results in the transposed matrix :)
        IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .NE. 0) THEN
        
          CALL storage_copy (rmatrix%h_Kcol,rtransposedMatrix%h_Kcol)
          CALL storage_copy (rmatrix%h_Kld,rtransposedMatrix%h_Kld)
          CALL storage_copy (rmatrix%h_Kdiagonal,rtransposedMatrix%h_Kdiagonal)
          rtransposedMatrix%imatrixSpec = &
            IAND(rtransposedMatrix%imatrixSpec,NOT(LSYSSC_MSPEC_TRANSPOSED))
          
        ELSE
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          CALL storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                            ST_INT, rtransposedMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
          CALL storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1, &
                            ST_INT, rtransposedMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
          
          ! Get Kcol/Kld
          CALL storage_getbase_int(rmatrix%h_Kcol,p_KcolSource)
          CALL storage_getbase_int(rmatrix%h_Kld,p_KldSource)

          CALL storage_getbase_int(rtransposedMatrix%h_Kcol,p_KcolDest)
          CALL storage_getbase_int(rtransposedMatrix%h_Kld,p_KldDest)
          
          ! We need a temporary array
          CALL storage_new ('lsyssc_transposeMatrix','Itemp',rmatrix%NCOLS,&
                            ST_INT, h_Itemp, ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int (h_Itemp,p_Itemp)
          
          ! Calculate the transposed matrix structure:
          CALL lsyssc_transpMatStruct79double (rmatrix%NEQ, rmatrix%NCOLS, &
                     p_KcolSource, p_KldSource, &
                     p_Itemp, p_KcolDest, p_KldDest)
                     
          ! Release the temp array
          CALL storage_free (h_Itemp)

          ! (Re-)calculate Kdiagonal
          CALL storage_new ('lsyssc_transposeMatrix', 'Kdiagonal', rmatrix%NCOLS, &
                            ST_INT, rtransposedMatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int(rtransposedMatrix%h_Kdiagonal,p_Kdiagonal)
          CALL lsyssc_rebuildKdiagonal (p_KcolDest, p_KldDest, p_Kdiagonal, &
                                        rmatrix%NCOLS)

        END IF

      CASE (LSYSSC_MATRIX7)
        ! Is the source matrix saved tranposed? In that case simply copy
        ! the data arrays, this results in the transposed matrix :)
        IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .NE. 0) THEN
          CALL storage_copy (rmatrix%h_Kcol,rtransposedMatrix%h_Kcol)
          CALL storage_copy (rmatrix%h_Kld,rtransposedMatrix%h_Kld)
          rtransposedMatrix%imatrixSpec = &
            IAND(rtransposedMatrix%imatrixSpec,NOT(LSYSSC_MSPEC_TRANSPOSED))
          
        ELSE
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          CALL storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                            ST_INT, rtransposedMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
          CALL storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1, &
                            ST_INT, rtransposedMatrix%h_Kld,ST_NEWBLOCK_NOINIT)

          ! Get Kcol/Kld
          CALL storage_getbase_int(rmatrix%h_Kcol,p_KcolSource)
          CALL storage_getbase_int(rmatrix%h_Kld,p_KldSource)

          CALL storage_getbase_int(rtransposedMatrix%h_Kcol,p_KcolDest)
          CALL storage_getbase_int(rtransposedMatrix%h_Kld,p_KldDest)
          
          ! We need a temporary array
          CALL storage_new ('lsyssc_transposeMatrix','Itemp',rmatrix%NA,&
                            ST_INT, h_Itemp, ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int (h_Itemp,p_Itemp)
          
          ! Calculate the transposed matrix structure:
          CALL lsyssc_transpMatStruct79double (rmatrix%NEQ, rmatrix%NCOLS, &
                     p_KcolSource, p_KldSource, &
                     p_Itemp, p_KcolDest, p_KldDest)
                     
          ! Release the temp array
          CALL storage_free (h_Itemp)

          ! Pivotise the matrix to move the diagonal element to the front.
          CALL lsyssc_pivotiseMatrix7double (rmatrix%NEQ, p_KcolDest, p_KldDest)
        END IF
        
      CASE DEFAULT
        PRINT *,'lsyssc_transposeMatrix: Unsupported matrix format.'
        STOP
      END SELECT
      
    CASE (LSYSSC_TR_ALL)
      ! Transpose the full matrix - structure+content
      SELECT CASE (rmatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9)
        ! Is the source matrix saved tranposed? In that case simply copy
        ! the data arrays, this results in the transposed matrix :)
        IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .NE. 0) THEN
        
          CALL storage_copy (rmatrix%h_Kcol,rtransposedMatrix%h_Kcol)
          CALL storage_copy (rmatrix%h_Kld,rtransposedMatrix%h_Kld)
          CALL storage_copy (rmatrix%h_Kdiagonal,rtransposedMatrix%h_Kdiagonal)
          rtransposedMatrix%imatrixSpec = &
            IAND(rtransposedMatrix%imatrixSpec,NOT(LSYSSC_MSPEC_TRANSPOSED))
          
        ELSE
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          CALL storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                            ST_INT, rtransposedMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
          CALL storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1, &
                            ST_INT, rtransposedMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
          CALL storage_new ('lsyssc_transposeMatrix', 'Da', rmatrix%NA, &
                            rtransposedMatrix%cdataType, rtransposedMatrix%h_Da,&
                            ST_NEWBLOCK_NOINIT)

          ! Get Kcol/Kld
          CALL storage_getbase_int(rmatrix%h_Kcol,p_KcolSource)
          CALL storage_getbase_int(rmatrix%h_Kld,p_KldSource)

          CALL storage_getbase_int(rtransposedMatrix%h_Kcol,p_KcolDest)
          CALL storage_getbase_int(rtransposedMatrix%h_Kld,p_KldDest)
          
          ! Get the data array(s)
          CALL storage_getbase_double(rMatrix%h_Da,p_DaSource)
          CALL storage_getbase_double(rtransposedMatrix%h_Da,p_DaDest)
          
          ! We need a temporary array
          CALL storage_new ('lsyssc_transposeMatrix','Itemp',rmatrix%NCOLS,&
                            ST_INT, h_Itemp, ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int (h_Itemp,p_Itemp)
          
          ! Calculate the transposed matrix structure:
          CALL lsyssc_transpMat79double (rmatrix%NEQ, rmatrix%NCOLS, &
                p_DaSource, p_KcolSource, p_KldSource, &
                p_Itemp, p_DaDest, p_KcolDest, p_KldDest)
                     
          ! Release the temp array
          CALL storage_free (h_Itemp)
          
          ! (Re-)calculate Kdiagonal
          CALL storage_new ('lsyssc_transposeMatrix', 'Kdiagonal', rmatrix%NCOLS, &
                            ST_INT, rtransposedMatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int(rtransposedMatrix%h_Kdiagonal,p_Kdiagonal)
          CALL lsyssc_rebuildKdiagonal (p_KcolDest, p_KldDest, p_Kdiagonal, &
                                        rmatrix%NCOLS)
          
        END IF

      CASE (LSYSSC_MATRIX7)
        ! Is the source matrix saved tranposed? In that case simply copy
        ! the data arrays, this results in the transposed matrix :)
        IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .NE. 0) THEN
        
          CALL storage_copy (rmatrix%h_Kcol,rtransposedMatrix%h_Kcol)
          CALL storage_copy (rmatrix%h_Kld,rtransposedMatrix%h_Kld)
          rtransposedMatrix%imatrixSpec = &
            IAND(rtransposedMatrix%imatrixSpec,NOT(LSYSSC_MSPEC_TRANSPOSED))
          
        ELSE
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          CALL storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                            ST_INT, rtransposedMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
          CALL storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1, &
                            ST_INT, rtransposedMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
          CALL storage_new ('lsyssc_transposeMatrix', 'Da', rmatrix%NA, &
                            rtransposedMatrix%cdataType, rtransposedMatrix%h_Da,&
                            ST_NEWBLOCK_NOINIT)

          ! Get Kcol/Kld
          CALL storage_getbase_int(rmatrix%h_Kcol,p_KcolSource)
          CALL storage_getbase_int(rmatrix%h_Kld,p_KldSource)

          CALL storage_getbase_int(rtransposedMatrix%h_Kcol,p_KcolDest)
          CALL storage_getbase_int(rtransposedMatrix%h_Kld,p_KldDest)
          
          ! Get the data array(s)
          CALL storage_getbase_double(rMatrix%h_Da,p_DaSource)
          CALL storage_getbase_double(rtransposedMatrix%h_Da,p_DaDest)

          ! We need a temporary array
          CALL storage_new ('lsyssc_transposeMatrix','Itemp',rmatrix%NA,&
                            ST_INT, h_Itemp, ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int (h_Itemp,p_Itemp)
          
          ! Calculate the transposed matrix structure:
          CALL lsyssc_transpMat79double (rmatrix%NEQ, rmatrix%NCOLS, &
                p_DaSource, p_KcolSource, p_KldSource, &
                p_Itemp, p_DaDest, p_KcolDest, p_KldDest)
                     
          ! Release the temp array
          CALL storage_free (h_Itemp)
          
          ! Pivotise the matrix to move the diagonal element to the front.
          CALL lsyssc_pivotiseMatrix7double (rmatrix%NEQ, p_KcolDest, p_KldDest, p_DaDest)
        END IF
      
      CASE DEFAULT
        PRINT *,'lsyssc_transposeMatrix: Unsupported matrix format.'
        STOP
      END SELECT
      
    CASE DEFAULT ! = LSYSSC_TR_ALL 
    END SELECT
  
    ! Exchange NEQ and NCOLS as these always describe the actual matrix.
    ntemp = rtransposedMatrix%NEQ
    rtransposedMatrix%NEQ = rtransposedMatrix%NCOLS
    rtransposedMatrix%NCOLS = ntemp
    
    ! That's it.
  
  END SUBROUTINE

END MODULE

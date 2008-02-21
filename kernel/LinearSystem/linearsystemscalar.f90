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
!#  3.) lsyssc_createVecIndMat
!#      -> Create a vector according to the matrix
!#
!#  4.) lsyssc_scalarProduct
!#      -> Calculate the scalar product of two vectors
!#
!#  5.) lsyssc_scalarMatVec
!#      -> Multiply a scalar matrix with a scalar vector
!#
!#  6.) lsyssc_releaseMatrix
!#      -> Release a scalar matrix from memory.
!#
!#  7.) lsyssc_releaseVector
!#      -> Release a scalar vector from memory.
!#
!#  8.) lsyssc_duplicateMatrix
!#      -> Create a duplicate of a given matrix or matrix-structure
!#
!#  9.) lsyssc_duplicateVector
!#      -> Create a duplicate of a given vector
!#
!# 10.) lsyssc_sortVectorInSitu
!#      -> Resort the entries of a vector or unsort them
!#
!# 11.) lsyssc_vectorActivateSorting
!#      -> Resort the entries of a vector or unsort them according to
!#         a previously attached sorting strategy
!#
!# 12.) lsyssc_synchroniseSortVecVec
!#      -> Synchronises the sorting strategy of a vector according to another
!#         vector.
!#
!# 13.) lsyssc_synchroniseSortMatVec
!#      -> Synchronises the sorting strategy of a vector according to a matrix.
!#
!# 14.) lsyssc_sortMatrix
!#      -> Resort the entries of a matrix or unsort them.
!#
!# 15.) lsyssc_unsortMatrix
!#      -> Unsorts a matrix.
!#
!# 16.) lsyssc_isVectorCompatible
!#      -> Checks whether two vectors are compatible to each other
!#
!# 17.) lsyssc_isMatrixCompatible = lsyssc_isMatrixVectorCompatible /
!#                                  lsyssc_isMatrixMatrixCompatible
!#      -> Checks whether a matrix and a vector are compatible to each other
!#
!# 18.) lsyssc_getbase_double
!#      -> Get a pointer to the double precision data array of a vector or a matrix
!#
!# 19.) lsyssc_getbase_single
!#      -> Get a pointer to the single precision data array of a vector or a matrix
!#
!# 20.) lsyssc_getbase_int
!#      -> Get a pointer to the integer data array of a vector
!#
!# 21.) lsyssc_getbase_Kcol
!#      -> Get a pointer to the integer data array Kcol of a matrix
!#         (if the matrix has one)
!#
!# 22.) lsyssc_getbase_Kld
!#      -> Get a pointer to the integer data array Kld of a matrix
!#         (if the matrix has one)
!#
!# 23.) lsyssc_getbase_Kdiagonal
!#      -> Get a pointer to the integer data array Kdiagonal of a matrix
!#         (if the matrix has one)
!#
!# 24.) lsyssc_addIndex
!#      -> Auxiliary routine. Adds an integer to each elememt of an integer 
!#         array.
!#
!# 25.) lsyssc_vectorNorm
!#      -> Calculate the norm of a vector.
!#
!# 26.) lsyssc_invertedDiagMatVec
!#      -> Multiply a vector with the inverse of the diagonal of a scalar
!#         matrix
!#
!# 27.) lsyssc_clearMatrix
!#      -> Clears a matrix, i.e. overwrites all entries with 0.0 or 
!#         with a defined value
!#
!# 28.) lsyssc_initialiseIdentityMatrix
!#      -> Initialises the content of a matrix to an identity matrix
!#
!# 29.) lsyssc_convertMatrix
!#      -> Allows to convert a matrix to another matrix structure.
!#
!# 30.) lsyssc_copyVector
!#       -> Copy a vector over to another one
!#
!# 31.) lsyssc_scaleVector
!#      -> Scale a vector by a constant
!#
!# 32.) lsyssc_clearVector
!#      -> Clear a vector, i.e. overwrites all entries with 0.0 or 
!#         with a defined value
!#
!# 33.) lsyssc_vectorLinearComb
!#      -> Linear combination of two vectors
!#
!# 34.) lsyssc_copyMatrix
!#      -> Copies a matrix to another one provided that they have the same 
!#         structure.
!#
!# 35.) lsyssc_transposeMatrix,
!#      lsyssc_transposeMatrixInSitu
!#      -> Transposes a scalar matrix.
!#
!# 36.) lsyssc_allocEmptyMatrix
!#      -> Allocates memory for the entries of a matrix.
!#
!# 37.) lsyssc_lumpMatrixScalar
!#      -> Performs lumping of a given matrix
!#
!# 38.) lsyssc_scaleMatrix
!#      -> Scale a matrix by a constant
!#
!# 39.) lsyssc_multMatMat
!#      -> Multiplies two matrices
!#
!# 40.) lsyssc_matrixLinearComb
!#      -> Adds two matrices
!#
!# 41.) lsyssc_swapVectors
!#      -> Swap two vectors
!#
!# 42.) lsyssc_isMatrixStructureShared
!#      -> Tests if the structure of a matrix is shared with another matrix
!#
!# 43.) lsyssc_isMatrixContentShared
!#      -> Tests if the content of a matrix is shared with another matrix
!#
!# 44.) lsyssc_resizeVector
!#      -> Resize the vector and reallocate memory if necessary.
!#
!# 45.) lsyssc_resizeMatrix
!#      -> Resize the matrix and reallocate memory if necessary.
!#
!# 46.) lsyssc_createDiagMatrixStruc
!#      -> Creates a diagonal matrix, does not allocate memory for the entries.
!#
!# 47.) lsyssc_clearOffdiags
!#      -> Clear all offdiagonal entries in a matrix.
!#
!# 48.) lsyssc_isMatrixSorted
!#      -> Checks if a matrix is currently sorted.
!#
!# 49.) lsyssc_isVectorSorted
!#      -> Checks if a vector is currently sorted.
!#
!# 50.) lsyssc_hasMatrixStructure
!#      -> Check if a matrix has a structure in memory or not.
!#
!# 51.) lsyssc_hasMatrixContent
!#      -> Check if a matrix has a content in memory or not.
!#
!# 52.) lsyssc_releaseMatrixContent
!#      -> Releases the content of the matrix, the structure will stay unchanged.
!#
!# 53.) lsyssc_spreadVector
!#      -> Spreads a scalar vector into another scalar vector
!#
!# 54.) lsyssc_spreadMatrix
!#      -> Spreads a scalar matrix into another scalar matrix
!#
!# 55.) lsyssc_packVector
!#      -> Packs a scalar vector into another scalar vector
!#
!# Sometimes useful auxiliary routines:
!#
!# 1.) lsyssc_rebuildKdiagonal (Kcol, Kld, Kdiagonal, neq)
!#     -> Rebuild the Kdiagonal array in a matrix of format 9
!#
!# 2.) lsyssc_infoMatrix
!#     -> Outputs information about the matrix (mostly used for debugging)
!#
!# 3.) lsyssc_infoVector
!#     -> Outputs information about the vector (mostly used for debugging)
!#
!# 4.) lsyssc_auxcopy_da
!#     -> Copies the data vector DA of a matrix to another without checks.
!#        If the destination array does not exist, it's created.
!#
!# 5.) lsyssc_auxcopy_Kcol
!#     -> Copies the data vector Kcol of a matrix to another without checks
!#        If the destination array does not exist, it's created.
!#
!# 6.) lsyssc_auxcopy_Kld
!#     -> Copies the data vector Kld of a matrix to another without checks
!#        If the destination array does not exist, it's created.
!#
!# 7.) lsyssc_auxcopy_Kdiagonal
!#     -> Copies the data vector Kdiagonal of a matrix to another without checks
!#        If the destination array does not exist, it's created.
!# </purpose>
!##############################################################################

MODULE linearsystemscalar

  USE fsystem
  USE storage
  USE spatialdiscretisation
  USE dofmapping
  USE genoutput

  IMPLICIT NONE

!<constants>

!<constantblock description="Global constants for scalar vectors/matrices">

  ! Maximum number of tags that can be assigned to a scalar vector or matrix.
  INTEGER, PARAMETER :: LSYSSC_MAXTAGS = 16
  
!</constantblock>

!<constantblock description="Global format flags for matrices">

  ! Unidentified matrix format
  INTEGER, PARAMETER :: LSYSSC_MATRIXUNDEFINED = 0
  
  ! Identifier for matrix format 1 - full matrix.
  ! Important matrix properties defining the matrix:
  ! NEQ = Number of rows, NCOLS = Number of columns, NA = Number of entries,
  ! h_Da = handle to matrix entries
  INTEGER, PARAMETER :: LSYSSC_MATRIX1 = 1
  
  ! Identifier for matrix format 9 - CSR
  ! Important matrix properties defining the matrix:
  ! NEQ = Number of rows, NCOLS = Number of columns, NA = Number of entries,
  ! h_Da        = handle to matrix entries,
  ! h_Kcol      = handle to column structure,
  ! h_Kld       = handle to row structure,
  ! h_Kdiagonal = handle to diagonal pointer
  INTEGER, PARAMETER :: LSYSSC_MATRIX9 = 9

  ! Identifier for matrix format 7 - CSR with diagonal element in front
  ! Important matrix properties defining the matrix:
  ! NEQ = Number of rows, NCOLS = Number of columns, NA = Number of entries,
  ! h_Da        = handle to matrix entries,
  ! h_Kcol      = handle to column structure,
  ! h_Kld       = handle to row structure
  INTEGER, PARAMETER :: LSYSSC_MATRIX7 = 7

  ! Identifier for matrix format D - Diagonal matrix, only entries on the
  ! main diagonal.
  ! Important matrix properties defining the matrix:
  ! NEQ = NCOLS = NA = Number of rows = Number of columns = Number of entries,
  ! h_Da        = handle to matrix entries
  INTEGER, PARAMETER :: LSYSSC_MATRIXD = 20

  ! Identifier for matrix format 7intl - CSR interleaved with diagonal element in front
  ! Important matrix properties defining the matrix:
  ! NEQ = Number of rows, NCOLS = Number of columns, 
  ! NA = Number of entries, NVAR = Number of variables 
  ! h_Da        = handle to matrix entries,
  ! h_Kcol      = handle to column structure,
  ! h_Kld       = handle to row structure
  INTEGER, PARAMETER :: LSYSSC_MATRIX7INTL = 70

  ! Identifier for matrix format 9intl - CSR interleaved
  ! Important matrix properties defining the matrix:
  ! NEQ = Number of rows, NCOLS = Number of columns, 
  ! NA = Number of entries, NVAR = Number of variables 
  ! h_Da        = handle to matrix entries,
  ! h_Kcol      = handle to column structure,
  ! h_Kld       = handle to row structure,
  ! h_Kdiagonal = handle to diagonal pointer
  INTEGER, PARAMETER :: LSYSSC_MATRIX9INTL = 90

!</constantblock>

!<constantblock description="Flags for the matrix specification bitfield">

  ! Standard matrix
  INTEGER(I32), PARAMETER :: LSYSSC_MSPEC_STANDARD =        0
  
  ! Matrix structure is a copy of another matrix, shared via the same
  ! handles. 
  INTEGER(I32), PARAMETER :: LSYSSC_MSPEC_STRUCTUREISCOPY = 2**0

  ! Matrix content is a copy of another matrix, shared via the same
  ! handles. 
  INTEGER(I32), PARAMETER :: LSYSSC_MSPEC_CONTENTISCOPY   = 2**1
  
  ! Complete matrix is duplicate of another matrix and shares structure
  ! and entries via the same pointers
  INTEGER(I32), PARAMETER :: LSYSSC_MSPEC_ISCOPY = LSYSSC_MSPEC_STRUCTUREISCOPY +&
                                                   LSYSSC_MSPEC_CONTENTISCOPY

  ! Matrix is saved transposed.
  ! To use a matrix in a transposed way, the application has to
  ! 1.) set this flag in the imatrixSpec bitfield
  ! 2.) exchange the values in t_matrixScalar\%NEQ and t_matrixScalar\%NCOLS 
  !     of the matrix structure.
  INTEGER(I32), PARAMETER :: LSYSSC_MSPEC_TRANSPOSED =      2**2

  ! Matrix not present in memory
  INTEGER(I32), PARAMETER :: LSYSSC_MSPEC_NOTINMEMORY =     2**3
  
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
  ! If the destination matrix  already contains allocated memory that belongs
  ! to that matrix, content/structure data is simply copied from rsourceMatrix into that.
  ! Note that this respects the ownership! I.e. if the destination matrix is not
  ! the owner of the content/structure data arrays, new memory is allocated to
  ! prevent the actual owner from getting destroyed!
  INTEGER, PARAMETER :: LSYSSC_DUP_COPY = 4
  
  ! The destination matrix gets a copy of the content of rsourceMatrix.
  ! If necessary, new memory is allocated.
  ! If the destination matrix  already contains allocated memory, content/structure 
  ! data is simply copied from rsourceMatrix into that.
  ! The ownership of the content/data arrays is not respected, i.e. if the
  ! destination matrix is not the owner, the actual owner of the data arrays is
  ! modified, too!
  INTEGER, PARAMETER :: LSYSSC_DUP_COPYOVERWRITE = 5
  
  ! Duplicate by ownership. What belongs to the source matrix is copied 
  ! (the same as LSYSSC_DUP_COPY). What belongs even to another matrix than
  ! the source matrix is shared (the same as LSYSSC_DUP_SHARE, .
  INTEGER, PARAMETER :: LSYSSC_DUP_ASIS = 6
  
  ! New memory is allocated for the structure/content in the same size as 
  ! in the source matrix but no data is copied; the arrays are left uninitialised.
  INTEGER, PARAMETER :: LSYSSC_DUP_EMPTY = 7 
                                               
  ! Copy the basic matrix information but don't copy handles.
  ! Set all handles of dynamic information to ST_NOHANDLE, so the matrix
  ! 'looks like' the old but has no dynamic data associated.
  INTEGER, PARAMETER :: LSYSSC_DUP_TEMPLATE   = 8
                 
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

!<constantblock description="Constants for initialising matrix entries when allocating">
  
  ! Let the entries of a matrix undefined
  INTEGER, PARAMETER :: LSYSSC_SETM_UNDEFINED = -1

  ! Clear the entries of a matrix when allocating
  INTEGER, PARAMETER :: LSYSSC_SETM_ZERO      = 0

  ! Set the entries of a matrix to 1 when allocating
  INTEGER, PARAMETER :: LSYSSC_SETM_ONE       = 1

!</constantblock>

!<constantblock description="Constants for lumping of matrices">
  
  ! Standard lumping; extract diagonal from a given matrix.
  INTEGER, PARAMETER :: LSYSSC_LUMP_STD       = 0

  ! Diagonal lumping; add all offdiagonal entries to the diagonal and take the diagonial
  INTEGER, PARAMETER :: LSYSSC_LUMP_DIAG      = 1

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
    
    ! Number of local variables; in general, scalar vectors of size
    ! NEQ posses NEQ entries. However, scalar vectors can be
    ! interleaved, that is, each of the NEQ entries stores NVAR local
    ! variables. In this case, NEQ remains unmodified but NVAR>1 such
    ! that the physical length of the vector is NEQ*NVAR.
    INTEGER         :: NVAR      = 1

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
    ! Looking from another viewpoint with the background of how a vector
    ! is renumbered, one can say:
    !  p_IsortPermutation (position in sorted vector) = position in unsorted vector.
    !  p_IsortPermutation (NEQ+position in unsorted vector) = position in sorted vector.
    ! Whether or not the vector is actually sorted depends on the
    ! flag isortStrategy!
    INTEGER         :: h_IsortPermutation = ST_NOHANDLE
    
    ! Start position of the vector data in the array identified by
    ! h_Ddata. Normally = 1. Can be set to > 1 if the vector is a subvector
    ! in a larger memory block allocated on the heap.
    INTEGER(PREC_VECIDX) :: iidxFirstEntry = 1
    
    ! Integer tags. This array of integer values can be used to store
    ! auxiliary tags which depend on the application.
    INTEGER, DIMENSION(LSYSSC_MAXTAGS) :: ITags = 0

    ! Real tags. This array of real values can be used to store
    ! auxiliary
    ! tags which depend on the application.
    REAL(DP), DIMENSION(LSYSSC_MAXTAGS) :: DTags = 0._DP

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
    ! Can take one of the LSYSSC_MATRIXx format flags.
    INTEGER      :: cmatrixFormat = LSYSSC_MATRIXUNDEFINED

    ! Format-tag. Identifies the format of the interleaved matrix.
    ! Can take one of the LSYSSC_MATRIX1/D format flags.
    INTEGER      :: cinterleavematrixFormat = LSYSSC_MATRIXUNDEFINED
    
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

    ! Number of variables for interleaved matrix.
    ! Remark: When an interleaved matrix is defined, then each
    ! position of the symbolic matrix structure can store a quadratic
    ! submatrix of dimension NVAR x NVAR. This is not to be confused
    ! with layers which means that several matrices share the same
    ! structure but their values can be treated individually. For
    ! interleaved matrices, the matrix-vector multiplication, etc. is
    ! adapted to the local dense blocks.
    INTEGER :: NVAR = 1
    
    ! Multiplier for matrix entries. All entries in the matrix are
    ! scaled by this multiplier when doing Matrix-vector multiplication.
    ! Note: This parameter is not supported by all algorithms, many
    ! algorithms will simply stop when this factor is <> 1, or they
    ! might ignore it. Therefore, use this factor with care!
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
    ! Looking from another viewpoint with the background of how a matrix
    ! is renumbered, one can say:
    !  p_IsortPermutation (column in sorted matrix) = column in unsorted matrix.
    !  p_IsortPermutation (NEQ+column in unsorted matrix) = column in sorted matrix.
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
    ! diagonal, this points to the diagonal element in each row
    ! (or to the first off-diagonal element above the diagonal if there
    ! is no diagonal element, respectively).
    ! For matrices with no elements in the upper right part of the
    ! matrix, this points to the first element on the next line.
    !INTEGER, DIMENSION(:), POINTER            :: Kdiagonal  => NULL()
    INTEGER                    :: h_Kdiagonal = ST_NOHANDLE

    ! Integer tags. This array of integer values can be used to store
    ! auxiliary tags which depend on the application.
    INTEGER, DIMENSION(LSYSSC_MAXTAGS) :: ITags = 0

    ! Real tags. This array of real values can be used to store
    ! auxiliary
    ! tags which depend on the application.
    REAL(DP), DIMENSION(LSYSSC_MAXTAGS) :: DTags = 0._DP
    
    ! A pointer to the spatial discretisation
    TYPE(t_spatialDiscretisation), POINTER     :: p_rspatialDiscretisation => NULL()
    
  END TYPE
  
!</typeblock>

!</types>

  INTERFACE lsyssc_getbase_double
    MODULE PROCEDURE lsyssc_getbaseVector_double
    MODULE PROCEDURE lsyssc_getbaseMatrixDA_double
  END INTERFACE

  INTERFACE lsyssc_getbase_single
    MODULE PROCEDURE lsyssc_getbaseVector_single
    MODULE PROCEDURE lsyssc_getbaseMatrixDA_single
  END INTERFACE

  INTERFACE lsyssc_getbase_int
    MODULE PROCEDURE lsyssc_getbaseVector_int
    MODULE PROCEDURE lsyssc_getbaseMatrixDA_int
  END INTERFACE

  INTERFACE lsyssc_isMatrixCompatible
    MODULE PROCEDURE lsyssc_isMatrixVectorCompatible
    MODULE PROCEDURE lsyssc_isMatrixMatrixCompatible
  END INTERFACE

  INTERFACE lsyssc_createVector
    MODULE PROCEDURE lsyssc_createVector
    MODULE PROCEDURE lsyssc_createVectorIntl
  END INTERFACE

  INTERFACE lsyssc_resizeVector
    MODULE PROCEDURE lsyssc_resizeVectorDirect
    MODULE PROCEDURE lsyssc_resizeVectorIndirect
    MODULE PROCEDURE lsyssc_resizeVectorIndMat
  END INTERFACE

  INTERFACE lsyssc_resizeMatrix
    MODULE PROCEDURE lsyssc_resizeMatrixDirect
    MODULE PROCEDURE lsyssc_resizeMatrixIndirect
  END INTERFACE

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
  IF (rvector1%NEQ .NE. rvector2%NEQ .OR. &
      rvector1%NVAR .NE. rvector2%NVAR) THEN
    IF (PRESENT(bcompatible)) THEN
      bcompatible = .FALSE.
      RETURN
    ELSE
      PRINT *,'Vectors not compatible, different size!'
      CALL sys_halt()
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
        CALL sys_halt()
      END IF
    END IF

    IF (rvector1%h_isortPermutation .NE. &
        rvector2%h_isortPermutation) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Vectors not compatible, differently sorted!'
        CALL sys_halt()
      END IF
    END IF
  END IF

  ! Ok, they are compatible
  IF (PRESENT(bcompatible)) bcompatible = .TRUE.

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_isMatrixVectorCompatible (rvector,rmatrix,btransposed,bcompatible)
  
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
  ! =FALSE: Check whether matrix-vector product $A*x$ is possible
  ! =TRUE : Check whether matrix-vector product $x^T*A = A^T*x$ is possible
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
      CALL sys_throwFPE()
      CALL sys_halt()
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
        CALL sys_halt()
      END IF
    END IF

    IF (rvector%h_isortPermutation .NE. &
        rmatrix%h_isortPermutation) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Vector/Matrix not compatible, differently sorted!'
        CALL sys_halt()
      END IF
    END IF
  END IF

  ! Ok, they are compatible
  IF (PRESENT(bcompatible)) bcompatible = .TRUE.

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_isMatrixMatrixCompatible (rmatrix1,rmatrix2,bcompatible)
  
!<description>
  ! Checks whether two matrices are compatible to each other, i.e. 
  ! share the same structure, size and sorting strategy.
!</description>

!<input>
  ! The vector
  TYPE(t_matrixScalar), INTENT(IN) :: rmatrix1
  
  ! The matrix
  TYPE(t_matrixScalar), INTENT(IN) :: rmatrix2
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether vector and matrix are compatible or not.
  ! If not given, an error will inform the user if the vector/matrix are
  ! not compatible and the program will halt.
  LOGICAL, INTENT(OUT), OPTIONAL :: bcompatible
!</output>

!</subroutine>

  ! We assume that we are not compatible
  IF (PRESENT(bcompatible)) bcompatible = .FALSE.
  
  ! Matrices must have the same size 
  IF ((rmatrix1%NEQ .NE. rmatrix2%NEQ) .OR. &
      (rmatrix1%NCOLS .NE. rmatrix2%NCOLS) .OR. &
      (rmatrix1%NA .NE. rmatrix2%NA)) THEN
    IF (PRESENT(bcompatible)) THEN
      bcompatible = .FALSE.
      RETURN
    ELSE
      PRINT *,'Matrices not compatible, different block structure!'
      CALL sys_halt()
    END IF
  END IF

  ! Matrices must have the same sparsity pattern
  IF (rmatrix1%NVAR .EQ. rmatrix2%NVAR) THEN
    
    ! We can perform restrictive checks, since both matrices have the
    ! same number of internal variables
    IF ((rmatrix1%cmatrixFormat .NE. rmatrix2%cmatrixFormat) .OR. &
        (rmatrix1%cinterleavematrixFormat .NE. rmatrix2%cinterleavematrixFormat)) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Matrices not compatible, different sparsity pattern!'
        CALL sys_halt()
      END IF
    END IF

  ELSE

    ! We have to be more careful.
    IF ((rmatrix1%cmatrixFormat .NE. rmatrix2%cmatrixFormat) .AND.&
        (rmatrix1%cmatrixFormat .NE. 10*rmatrix2%cmatrixFormat) .AND.&
        (10*rmatrix1%cmatrixFormat .NE. rmatrix2%cmatrixFormat)) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Matrices not compatible, different sparsity pattern!'
        CALL sys_halt()
      END IF
    END IF

  END IF

  ! isortStrategy < 0 means unsorted. Both unsorted is ok.

  IF ((rmatrix1%isortStrategy .GT. 0) .OR. &
      (rmatrix2%isortStrategy .GT. 0)) THEN

    IF (rmatrix1%isortStrategy .NE. &
        rmatrix2%isortStrategy) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Matrices not compatible, differently sorted!'
        CALL sys_halt()
      END IF
    END IF

    IF (rmatrix1%h_isortPermutation .NE. &
        rmatrix2%h_isortPermutation) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        PRINT *,'Matrices not compatible, differently sorted!'
        CALL sys_halt()
      END IF
    END IF
  END IF

  ! Ok, they are compatible
  IF (PRESENT(bcompatible)) bcompatible = .TRUE.

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_checkDiscretisation (rvector,rdiscretisation,bcompatible)
  
!<description>
  ! Checks whether a given discretisation structure rdiscretisation is
  ! basically compatible to a given vector (concerning NEQ, cubature formulas
  ! on element distribution, element compatibility,...)
!</description>

!<input>
  ! The vector which is to be checked.
  TYPE(t_vectorScalar), INTENT(IN) :: rvector

  ! A scalar discretisation structure to check against the vector.
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscretisation
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether vector and matrix are compatible or not.
  ! If not given, an error will inform the user if the vector/matrix are
  ! not compatible and the program will halt.
  LOGICAL, INTENT(OUT), OPTIONAL :: bcompatible
!</output>

!</subroutine>

    INTEGER(PREC_DOFIDX) :: NEQ

    ! We assume that we are compatible
    IF (PRESENT(bcompatible)) bcompatible = .TRUE.

    ! NEQ must be correct.
    NEQ = dof_igetNDofGlob(rdiscretisation)
    IF (NEQ .NE. rvector%NEQ) THEN
      IF (PRESENT(bcompatible)) THEN
        bcompatible = .FALSE.
        RETURN
      ELSE
        CALL output_line('Vector not compatible to discretisation, different NEQ!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_checkDiscretisation')
        CALL sys_halt()
      END IF
    END IF
    
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_getbaseVector_double (rvector,p_Ddata)
  
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
    CALL output_line('Trying to access empty vector!', &
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbase_double')
    CALL sys_halt()
  END IF

  ! Check that the vector is really double precision
  IF (rvector%cdataType .NE. ST_DOUBLE) THEN
    CALL output_line('Vector is of wrong precision!', &
       OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbase_double')
    CALL sys_halt()
  END IF

  ! Get the data array
  CALL storage_getbase_double (rvector%h_Ddata,p_Ddata)
  
  ! Modify the starting address/length to get the real array.
  p_Ddata => p_Ddata(rvector%iidxFirstEntry:&
      rvector%iidxFirstEntry+rvector%NEQ*rvector%NVAR-1)
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_getbaseVector_single (rvector,p_Fdata)
  
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
    CALL sys_halt()
  END IF

  ! Get the data array
  CALL storage_getbase_single (rvector%h_Ddata,p_Fdata)
  
  ! Modify the starting address/length to get the real array.
  p_Fdata => p_Fdata(rvector%iidxFirstEntry:&
      rvector%iidxFirstEntry+rvector%NEQ*rvector%NVAR-1)
  
  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_getbaseVector_int (rvector,p_Idata)
  
!<description>
  ! Returns a pointer to the integer data array of the vector.
  ! An error is thrown if the vector is not integer.
!</description>

!<input>
  ! The vector
  TYPE(t_vectorScalar), INTENT(IN) :: rvector
!</input>

!<output>
  ! Pointer to the integer data array of the vector.
  ! NULL() if the vector has no data array.
  INTEGER, DIMENSION(:), POINTER :: p_Idata
!</output>

!</subroutine>

  ! Do we have data at all?
  IF ((rvector%NEQ .EQ. 0) .OR. (rvector%h_Ddata .EQ. ST_NOHANDLE)) THEN
    NULLIFY(p_Idata)
    RETURN
  END IF

  ! Check that the vector is really integer
  IF (rvector%cdataType .NE. ST_INT) THEN
    PRINT *,'lsyssc_getbase_int: Vector is of wrong precision!'
    CALL sys_halt()
  END IF

  ! Get the data array
  CALL storage_getbase_int (rvector%h_Ddata,p_Idata)
  
  ! Modify the starting address/length to get the real array.
  p_Idata => p_Idata(rvector%iidxFirstEntry:&
      rvector%iidxFirstEntry+rvector%NEQ*rvector%NVAR-1)
  
  END SUBROUTINE

  ! ***************************************************************************

!<function>

  LOGICAL FUNCTION lsyssc_isExplicitMatrix1D (rmatrix)

!<description>
  ! Checks whether a given matrix explicitly exists in memory as 1D array,
  ! which can be accessed via lsyssc_getbase_double. This is the usual case
  ! for most of the matrices (format-7, format-9,...) but might be different
  ! for special-type matrices that do not exist in memory, i.e. where it's
  ! only known how to apply them to a vector.
!</description>

!<input>
    ! The matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
!</input>

!<result>
    ! TRUE if the matrix is an explicit matrix that can be accessed with
    ! lsyssc_getbase_double.
    ! FALSE if the matrix is not accessible in that way.
!</result>

!</function>

    ! Check the matrix type. If we have a matrix that is known to be in memory,
    ! we can return TRUE.
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX1, LSYSSC_MATRIX9,     LSYSSC_MATRIX7, &
          LSYSSC_MATRIXD, LSYSSC_MATRIX7INTL, LSYSSC_MATRIX9INTL)
      lsyssc_isExplicitMatrix1D = .TRUE.
      
    CASE DEFAULT
      lsyssc_isExplicitMatrix1D = .FALSE.
      
    END SELECT

  END FUNCTION

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_getbaseMatrixDA_double (rmatrix, p_Ddata)

!<description>
    ! Returns a pointer to the double precision data array of the matrix.
    ! An error is thrown if the matrix is not double precision.
!</description>

!<input>
    ! The matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
!</input>

!<output>
    ! Pointer to the double precision data array of the matrix.
    ! NULL() if the matrix has no data array.
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
!</output>

!</subroutine>

    ! Check if the matrix format allows to access DA
    IF (.NOT. lsyssc_isExplicitMatrix1D(rmatrix)) THEN
      PRINT *,'lsyssc_getbaseMatrixDA_double: Matrix does not exist explicitely!'
      CALL sys_halt()
    END IF

    ! Do we have data at all?
    IF ((rmatrix%NA .EQ. 0) .OR. (rmatrix%h_DA .EQ. ST_NOHANDLE)) THEN
      NULLIFY(p_Ddata)
      RETURN
    END IF

    ! Check that the matrix is really double precision
    IF (rmatrix%cdataType .NE. ST_DOUBLE) THEN
      PRINT *,'lsyssc_getbaseMatrixDA_double: Matrix is of wrong precision!'
      CALL sys_halt()
    END IF
    
    ! Get the data array
    SELECT CASE(rmatrix%cinterleavematrixFormat)
    CASE (LSYSSC_MATRIX1)
      CALL storage_getbase_double (rmatrix%h_DA,p_Ddata,rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR)
    CASE (LSYSSC_MATRIXD)
      CALL storage_getbase_double (rmatrix%h_DA,p_Ddata,rmatrix%NA*rmatrix%NVAR)
    CASE DEFAULT
      CALL storage_getbase_double (rmatrix%h_DA,p_Ddata,rmatrix%NA)
    END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_getbaseMatrixDA_single (rmatrix, p_Fdata)

!<description>
    ! Returns a pointer to the single precision data array of the matrix.
    ! An error is thrown if the matrix is not single precision.
!</description>

!<input>
    ! The matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
!</input>

!<output>
    ! Pointer to the single precision data array of the matrix.
    ! NULL() if the matrix has no data array.
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
!</output>

!</subroutine>

    ! Check if the matrix format allows to access DA
    IF (.NOT. lsyssc_isExplicitMatrix1D(rmatrix)) THEN
      PRINT *,'lsyssc_getbaseMatrixDA_single: Matrix does not exist explicitely!'
      CALL sys_halt()
    END IF

    ! Do we have data at all?
    IF ((rmatrix%NA .EQ. 0) .OR. (rmatrix%h_DA .EQ. ST_NOHANDLE)) THEN
      NULLIFY(p_Fdata)
      RETURN
    END IF

    ! Check that the matrix is really single precision
    IF (rmatrix%cdataType .NE. ST_SINGLE) THEN
      PRINT *,'lsyssc_getbaseMatrix_single: Matrix is of wrong precision!'
      CALL sys_halt()
    END IF

    ! Get the data array
    SELECT CASE(rmatrix%cinterleavematrixFormat)
    CASE (LSYSSC_MATRIX1)
      CALL storage_getbase_single (rmatrix%h_DA,p_Fdata,rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR)
    CASE(LSYSSC_MATRIXD)
      CALL storage_getbase_single (rmatrix%h_DA,p_Fdata,rmatrix%NA*rmatrix%NVAR)
    CASE DEFAULT
      CALL storage_getbase_single (rmatrix%h_DA,p_Fdata,rmatrix%NA)
    END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_getbaseMatrixDA_int (rmatrix, p_Idata)

!<description>
    ! Returns a pointer to the integer data array of the matrix.
    ! An error is thrown if the matrix is not integer.
!</description>

!<input>
    ! The matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
!</input>

!<output>
    ! Pointer to the integer data array of the matrix.
    ! NULL() if the matrix has no data array.
    INTEGER, DIMENSION(:), POINTER :: p_Idata
!</output>

!</subroutine>

    ! Check if the matrix format allows to access DA
    IF (.NOT. lsyssc_isExplicitMatrix1D(rmatrix)) THEN
      PRINT *,'lsyssc_getbaseMatrixDA_int: Matrix does not exist explicitely!'
      CALL sys_halt()
    END IF

    ! Do we have data at all?
    IF ((rmatrix%NA .EQ. 0) .OR. (rmatrix%h_DA .EQ. ST_NOHANDLE)) THEN
      NULLIFY(p_Idata)
      RETURN
    END IF

    ! Check that the matrix is really integer
    IF (rmatrix%cdataType .NE. ST_INT) THEN
      PRINT *,'lsyssc_getbaseMatrix_int: Matrix is of wrong precision!'
      CALL sys_halt()
    END IF

    ! Get the data array
    SELECT CASE(rmatrix%cinterleavematrixFormat)
    CASE(LSYSSC_MATRIX1)
      CALL storage_getbase_int (rmatrix%h_DA,p_Idata,rmatrix%NA*rmatrix%NVAR)
    CASE(LSYSSC_MATRIXD)
      CALL storage_getbase_int (rmatrix%h_DA,p_Idata,rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR)
    CASE DEFAULT
      CALL storage_getbase_int (rmatrix%h_DA,p_Idata,rmatrix%NA)
    END SELECT

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_getbase_Kcol (rmatrix,p_Kcol)

!<description>
    ! Returns a pointer to the column data array of the matrix.
    ! An error is thrown if the matrix does not provide column array.
!</description>

!<input>
    ! The matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
!</input>

!<output>
    ! Pointer to the column array of the matrix.
    ! NULL() if the matrix has no column array.
    INTEGER, DIMENSION(:), POINTER :: p_Kcol
!</output>

!</subroutine>

    ! Is matrix in correct format?
    IF ((rmatrix%cmatrixFormat /= LSYSSC_MATRIX7) .AND. &
        (rmatrix%cmatrixFormat /= LSYSSC_MATRIX7INTL) .AND.&
        (rmatrix%cmatrixFormat /= LSYSSC_MATRIX9) .AND. &
        (rmatrix%cmatrixFormat /= LSYSSC_MATRIX9INTL)) THEN
      PRINT *,'lsyssc_getbase_Kcol: matrix format does not provide KCOL!'
      CALL sys_halt()
    END IF

    ! Do we have a column array at all?
    IF ((rmatrix%NA .EQ. 0) .OR. (rmatrix%h_Kcol .EQ. ST_NOHANDLE)) THEN
      NULLIFY(p_Kcol)
      RETURN
    END IF

    ! Get the column array
    CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol,rmatrix%NA)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_getbase_Kld (rmatrix,p_Kld)

!<description>
    ! Returns a pointer to the column offset data array of the matrix.
    ! An error is thrown if the matrix does not provide column offset array.
!</description>

!<input>
    ! The matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
!</input>

!<output>
    ! Pointer to the column offset array of the matrix.
    ! NULL() if the matrix has no column offset array.
    INTEGER, DIMENSION(:), POINTER :: p_Kld
!</output>

!</subroutine>

    ! Is matrix in correct format?
    IF ((rmatrix%cmatrixFormat /= LSYSSC_MATRIX7) .AND. &
        (rmatrix%cmatrixFormat /= LSYSSC_MATRIX7INTL) .AND.&
        (rmatrix%cmatrixFormat /= LSYSSC_MATRIX9) .AND. &
        (rmatrix%cmatrixFormat /= LSYSSC_MATRIX9INTL)) THEN
      PRINT *,'lsyssc_getbase_Kld: matrix format does not provide KLD!'
      CALL sys_halt()
    END IF

    ! Do we have a column offset array at all?
    IF ((rmatrix%NEQ .EQ. 0) .OR. (rmatrix%h_Kld .EQ. ST_NOHANDLE)) THEN
      NULLIFY(p_Kld)
      RETURN
    END IF

    ! Get the column offset array.
    ! Take care: If the matrix is virtually transposed, NCOLS and NEQ are
    ! exchanged!
    IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .EQ. 0) THEN
      CALL storage_getbase_int (rmatrix%h_Kld,p_Kld,rmatrix%NEQ+1)
    ELSE
      CALL storage_getbase_int (rmatrix%h_Kld,p_Kld,rmatrix%NCOLS+1)
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)

!<description>
    ! Returns a pointer to the diagonal data array of the matrix.
    ! An error is thrown if the matrix does not provide diagonal array.
!</description>

!<input>
    ! The matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
!</input>

!<output>
    ! Pointer to the diagonal array of the matrix.
    ! NULL() if the matrix has no diagonal array.
    INTEGER, DIMENSION(:), POINTER :: p_Kdiagonal
!</output>

!</subroutine>

    ! Is matrix in correct format?
    IF ((rmatrix%cmatrixFormat /= LSYSSC_MATRIX9) .AND.&
        (rmatrix%cmatrixFormat /= LSYSSC_MATRIX9INTL)) THEN
      PRINT *,'lsyssc_getbase_Kdiagonal: matrix format does not provide KDIAGONAL!'
      CALL sys_halt()
    END IF

    ! Do we have a column offset array at all?
    IF ((rmatrix%NEQ .EQ. 0) .OR. (rmatrix%h_Kdiagonal .EQ. ST_NOHANDLE)) THEN
      NULLIFY(p_Kdiagonal)
      RETURN
    END IF

    ! Get the column offset array
    CALL storage_getbase_int (rmatrix%h_Kdiagonal,p_Kdiagonal,rmatrix%NEQ)

  END SUBROUTINE

  !****************************************************************************

!<subroutine>
  
  SUBROUTINE lsyssc_createVector (rvector,NEQ,bclear,cdataType,NEQMAX)
  
!<description>
  ! This creates a simple scalar vector of length NEQ. Memory is 
  ! allocated on the heap and can be released by lsyssc_releaseVector.
  !
  ! Note, if the optional parameter NEQMAX is given, then memory is
  ! allocated for a vector of length NEQMAX but only length NEQ is
  ! assigned to the vector. The vector can be resized arbitrarily.
  ! Note that no memory reallocation is required if NEQ < NEQMAX.
  ! In order to keep the actual size of the memory transparent from
  ! the user, NEQMAX is not stored directly. It can only be obtained,
  ! by getting the size of the associated storage block.
!</description>
  
!<input>
  
  ! Desired length of the vector
  INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ

  ! Whether to fill the vector with zero initially
  LOGICAL, INTENT(IN)                               :: bclear

  ! OPTIONAL: Data type of the vector.
  ! If not specified, ST_DOUBLE is assumed.
  INTEGER, INTENT(IN), OPTIONAL                     :: cdataType  

  ! OPTIONAL: Maximum length of the vector
  INTEGER(PREC_VECIDX), INTENT(IN), OPTIONAL        :: NEQMAX
  
!</input>

!<output>
  ! Scalar vector structure
  TYPE(t_vectorScalar), INTENT(OUT)                 :: rvector
!</output>

!</subroutine>

  INTEGER ::cdata
  INTEGER(PREC_VECIDX) :: isize

    cdata = ST_DOUBLE
    IF (PRESENT(cdataType)) cdata=cdataType
    isize = MAX(0,NEQ)
    IF (PRESENT(NEQMAX)) isize=MAX(isize,NEQMAX)

    ! The INTENT(OUT) already initialises rvector with the most important
    ! information. The rest comes now:
    !
    ! Size:
    rvector%NEQ = MAX(0,NEQ)
    
    ! Handle - if NEQ > 0
    IF (rvector%NEQ .GT. 0) THEN
      IF (bclear) THEN
        CALL storage_new1D ('lsyssc_createVector', 'ScalarVector', isize, &
                            cdata, rvector%h_Ddata,ST_NEWBLOCK_ZERO)
      ELSE
        CALL storage_new1D ('lsyssc_createVector', 'ScalarVector', isize, &
                            cdata, rvector%h_Ddata,ST_NEWBLOCK_NOINIT)
      END IF
    END IF
  
  END SUBROUTINE

!****************************************************************************

!<subroutine>
  
  SUBROUTINE lsyssc_createVectorIntl (rvector,NEQ,NVAR,bclear,cdataType,NEQMAX)
  
!<description>
  ! This creates an interleaved scalar vector of length NEQ with NVAR
  ! local variables. Memory is allocated on the heap and can be
  ! released by lsyssc_releaseVector.
  !
  ! Note, if the optional parameter NEQMAX is given, then memory is
  ! allocated for a vector of length NEQMAX but only length NEQ is
  ! assigned to the vector. The vector can be resized arbitrarily.
  ! Note that no memory reallocation is required if NEQ < NEQMAX.
  ! In order to keep the actual size of the memory transparent from
  ! the user, NEQMAX is not stored directly. It can only be obtained,
  ! by getting the size of the associated storage block.
!</description>
  
!<input>
  
  ! Desired length of the vector
  INTEGER(PREC_VECIDX), INTENT(IN)                  :: NEQ

  ! Desired number of local variables
  INTEGER, INTENT(IN)                               :: NVAR

  ! Whether to fill the vector with zero initially
  LOGICAL, INTENT(IN)                               :: bclear

  ! OPTIONAL: Data type of the vector.
  ! If not specified, ST_DOUBLE is assumed.
  INTEGER, INTENT(IN), OPTIONAL                     :: cdataType  

  ! OPTIONAL: Maximum length of the vector
  INTEGER(PREC_VECIDX), INTENT(IN), OPTIONAL        :: NEQMAX
  
!</input>

!<output>
  ! Scalar vector structure
  TYPE(t_vectorScalar), INTENT(OUT)                 :: rvector
!</output>

!</subroutine>

  INTEGER ::cdata
  INTEGER(PREC_VECIDX) :: isize

    cdata = ST_DOUBLE
    IF (PRESENT(cdataType)) cdata=cdataType
    isize = MAX(0,NEQ)
    IF (PRESENT(NEQMAX)) isize=MAX(isize,NEQMAX)
    isize = isize*MAX(0,NVAR)

    ! The INTENT(OUT) already initialises rvector with the most important
    ! information. The rest comes now:
    !
    ! Size:
    rvector%NEQ = MAX(0,NEQ)
    rvector%NVAR= MAX(0,NVAR)
    
    ! Handle - if NEQ > 0
    IF (rvector%NEQ*rvector%NVAR .GT. 0) THEN
      IF (bclear) THEN
        CALL storage_new1D ('lsyssc_createVector', 'ScalarVector', isize, &
                            cdata, rvector%h_Ddata,ST_NEWBLOCK_ZERO)
      ELSE
        CALL storage_new1D ('lsyssc_createVector', 'ScalarVector', isize, &
                            cdata, rvector%h_Ddata,ST_NEWBLOCK_NOINIT)
      END IF
    END IF
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_createVecByDiscr (rdiscretisation,rx,bclear,cdataType,NEQMAX)
  
!<description>
  ! Initialises the vector structure rx based on a discretisation
  ! structure rDiscretisation. 
  !
  ! Memory is allocated on the heap for rx accordint to the number of 
  ! DOF's indicated by the spatial discretisation structures in 
  ! rdiscretisation.
  !
  ! Note, if the optional parameter NEQMAX is given, then memory is
  ! allocated for a vector of length NEQMAX but only length NEQ is
  ! assigned to the vector. The vector can be resized arbitrarily.
  ! Note that no memory reallocation is required if NEQ < NEQMAX.
  ! In order to keep the actual size of the memory transparent from
  ! the user, NEQMAX is not stored directly. It can only be obtained,
  ! by getting the size of the associated storage block.
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

  ! OPTIONAL: Maximum length of the vector
  INTEGER(PREC_VECIDX), INTENT(IN), OPTIONAL :: NEQMAX
!</input>

!<output>
  ! Destination structure. Memory is allocated appropriately.
  ! A pointer to rdiscretisation is saved to r.
  TYPE(t_vectorScalar),INTENT(OUT) :: rx
!</output>
  
!</subroutine>

  INTEGER :: cdata
  INTEGER(PREC_VECIDX) :: NEQ
  LOGICAL :: bcl
  
  cdata = ST_DOUBLE
  IF (PRESENT(cdataType)) cdata = cdataType
  
  bcl = .FALSE.
  IF (PRESENT(bclear)) bcl = bclear
  
  ! Get NEQ:
  NEQ = dof_igetNDofGlob(rdiscretisation)
  
  ! Create a new vector with that block structure
  CALL lsyssc_createVector (rx, NEQ, bcl, cdataType, NEQMAX)
  
  ! Initialise further data of the block vector
  rx%p_rspatialDiscretisation => rdiscretisation
  
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_createVecIndMat (rtemplateMat,rx,bclear,cdataType)

!<description>
    ! Initializes the scalar vector structure rx. rtemplateMat is an 
    ! existing scalar matrix structure. The vector rx will be created
    ! according to the size of the template matrix.
    !
    ! Memory is allocated on the heap for rx. The vector rx will have
    ! the same size as the number of matrix columns.
    ! The sorting strategy of the vector is initialized with the 
    ! sorting strategy of the template matrix rtemplateMat.
!</description>

!<input>
    ! A template matrix structure
    TYPE(t_matrixScalar), INTENT(IN) :: rtemplateMat

    ! OPTIONAL: If set to TRUE, the vector will be filled with zero initially.
    ! Otherwise the content of rx is undefined.
    LOGICAL, INTENT(IN), OPTIONAL   :: bclear
    
    ! OPTIONAL: Data type identifier for the entries in the vector. 
    ! Either ST_SINGLE or ST_DOUBLE. If not present, ST_DOUBLE is assumed.
    INTEGER, INTENT(IN),OPTIONAL              :: cdataType
!</input>
    
!<output>
    ! Destination structure. Memory is allocated for each of the blocks.
    TYPE(t_vectorScalar),INTENT(OUT) :: rx
!</output>

!</subroutine>

    ! local variables
    INTEGER :: cdata

    cdata = ST_DOUBLE
    IF (PRESENT(cdataType)) cdata = cdataType
    
    ! Allocate memory for vector
    CALL storage_new1D ('lsyssc_createVecIndMat', 'Vector', &
        rtemplateMat%NCOLS*rtemplateMat%NVAR, &
        cdata, rx%h_Ddata, ST_NEWBLOCK_NOINIT)

    ! Set structure
    rx%NEQ       = rtemplateMat%NCOLS
    rx%NVAR      = rtemplateMat%NVAR
    rx%cdataType = cdata

    ! Transfer discretization pointers from the matrix to the vector
    rx%p_rspatialDiscretisation => rtemplateMat%p_rspatialDiscretisation

    ! Transfer sorting strategy from the matrix to the vector
    rx%isortStrategy      = rtemplateMat%isortStrategy
    rx%h_IsortPermutation = rtemplateMat%h_IsortPermutation

    IF (PRESENT(bclear)) THEN
      IF (bclear) THEN
        CALL lsyssc_clearVector (rx)
      END IF
    END IF

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_resizeVectorDirect (rvector, NEQ, bclear, bcopy, NEQMAX)

!<description>
  ! Resizes the vector structure to the new size NEQ which is given explicitely.
  !
  ! If NEQ is smaller than the real memory allocated for the vector, then
  ! only NEQ is reset. Otherwise, the memory is physically reallocated.
  !
  ! If the parameter bclear=.TRUE. the complete vector is clear.
  ! This is done both in case the vector is reallocated physically or not.
  !
  ! If the parameter bcopy=.TRUE. the content of the existing vector is 
  ! copied. This is only required, if the vector needs to be reallocated
  ! physically. Otherwise, the data still exists in the memory.
  !
  ! Remark: The parameter bclear has higher priority than the parameter bcopy.
  ! That is, if both parameters are given, then no data is copied if the
  ! vector shoud be cleared afterwards ;-)
  !
  ! If the optional parameter NEQMAX is specified, then memory of size
  ! NEQMAX is allocated if reallocate is required. Otherwise, no reallocation
  ! is performed an NEQMAX is just neglected.
!</description>

!<input>

    ! Desired length of the vector
    INTEGER(PREC_VECIDX), INTENT(IN)           :: NEQ  

    ! Whether to fill the vector with zero initially
    LOGICAL, INTENT(IN)                        :: bclear

    ! OPTIONAL: Whether to copy the content of the vector to the resized one
    LOGICAL, INTENT(IN), OPTIONAL              :: bcopy

    ! OPTIONAL: Maximum length of the vector
    INTEGER(PREC_VECIDX), INTENT(IN), OPTIONAL :: NEQMAX

!</input>

!<inputoutput>

    ! Scalar vector structure
    TYPE(t_vectorScalar), INTENT(INOUT)         :: rvector

!</inputoutput>

!</subroutine>

    INTEGER(PREC_VECIDX) :: iNEQ,isize
    LOGICAL :: bdocopy

    ! Check, that vector is not a copy of another (possibly larger) vector
    IF (rvector%bisCopy) THEN
      PRINT *, "lsyssc_resizeVectorDirect: A copied vector cannot be resized!"
      CALL sys_halt()
    END IF
    
    ! Check, if vector has been initialized before.
    IF (rvector%NEQ == 0 .OR. rvector%h_Ddata == ST_NOHANDLE) THEN
      PRINT *, "lsyssc_resizeVectorDirect: A vector can only be resized &
          & if it has been created correctly!"
      CALL sys_halt()
    END IF
    
    ! Set working dimensions
    iNEQ = MAX(0,NEQ)
    IF (PRESENT(NEQMAX)) iNEQ = MAX(iNEQ,NEQMAX)

    ! Set copy/clear attributes
    bdocopy = (.NOT.bclear)
    IF (PRESENT(bcopy)) bdocopy = (bdocopy .AND. bcopy)
    
    ! If the vector should be cleared, then the sorting strategy (if any)
    ! can be ignored and reset. Otherwise, the vector needs to be unsorted
    ! prior to copying some part of it. Afterwards, no sorting strategy is
    ! available in any case.
    IF (bdocopy .AND. rvector%isortStrategy > 0) THEN
      CALL lsyssc_vectorActivateSorting(rvector,.FALSE.)
    END IF
    
    ! Reset sorting strategy, there is none
    rvector%isortStrategy      = 0
    rvector%h_iSortPermutation = ST_NOHANDLE

    ! Get current size of vector memory
    CALL storage_getsize(rvector%h_Ddata, isize)

    ! Update NEQ. Note that NVAR cannot be modified.
    rvector%NEQ = NEQ

    ! Do we really have to reallocate the vector physically?
    IF (rvector%NVAR*rvector%NEQ > isize) THEN
      
      ! Yes, so adopt the new size. Note that some extra memory is
      ! allocated if the optional argument NEQMAX is specified
      isize = iNEQ*rvector%NVAR
      
      ! Reallocate the memory for handle h_Ddata
      CALL storage_realloc('lsyssc_resizeVectorDirect', isize, rvector%h_Ddata, &
          ST_NEWBLOCK_NOINIT, bdocopy)

    ELSEIF (PRESENT(NEQMAX)) THEN
      
      ! The available memory suffices for the vector, i.e. isize <= NVAR*NEQ.
      ! Let's check if the user supplied a new upper limit which makes it 
      ! mandatory to "shrink" the allocated memory. Note that memory for
      ! at least NEQ vector intries is allocated.
      IF (isize > rvector%NVAR*iNEQ) THEN
        
        ! Compute new size, i.e. MAX(0,NEQ,NEQMAX)
        isize = rvector%NVAR*iNEQ

        IF (isize == 0) THEN
          ! If nothing is left, then the vector can also be released.
          CALL lsyssc_releaseVector(rvector)
          RETURN
        ELSE
          CALL storage_realloc('lsyssc_resizeVectorDirect', isize, rvector%h_Ddata, &
              ST_NEWBLOCK_NOINIT, bdocopy)
        END IF
      END IF
    END IF
    
    ! Should the vector be cleared?
    IF (bclear) CALL lsyssc_clearVector(rvector)

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_resizeVectorIndirect (rvector, rvectorTemplate, bclear, bcopy)

!<description>
  ! Resizes the vector structure so that it exhibits the same memory layout 
  ! as the template vector. Note that this subroutine can only be used
  ! if the template vector has the same internal structure, e.g., data type, NVAR
  ! as the vector to be resized.
  !
  ! If NEQ is smaller than the real memory allocated for the vector, then
  ! only NEQ is reset. Otherwise, the memory is physically reallocated.
  !
  ! If the parameter bclear=.TRUE. the complete vector is cleared.
  ! This is done both in case the vector is reallocated physically or not.
  !
  ! If the parameter bcopy=.TRUE. the content of the existing vector is 
  ! copied. This is only required, if the vector needs to be reallocated
  ! physically. Otherwise, the data still exists in the memory.
  !
  ! Remark: The parameter bclear has higher priority than the parameter bcopy.
  ! That is, if both parameters are given, then no data is copied if the
  ! vector shoud be cleared afterwards ;-)
!</description>

!<input>

    ! Scalar template vector
    TYPE(t_vectorScalar), INTENT(IN)           :: rvectorTemplate

    ! Whether to fill the vector with zero initially
    LOGICAL, INTENT(IN)                        :: bclear

    ! OPTIONAL: Whether to copy the content of the vector to the resized one
    LOGICAL, INTENT(IN), OPTIONAL              :: bcopy

!</input>

!<inputoutput>

    ! Scalar vector structure
    TYPE(t_vectorScalar), INTENT(INOUT)         :: rvector

!</inputoutput>

!</subroutine>

    INTEGER(PREC_VECIDX) :: isize,iNEQMAX

    ! Check, if vector is a copy of another (possibly larger) vector
    IF (rvector%bisCopy) THEN
      PRINT *, "lsyssc_resizeVectorIndirect: A copied vector cannot be resized!"
      CALL sys_halt()
    END IF

    ! Check, if vector has been initialized before. If this is not
    ! the case, then create a new vector by duplicating the template vector.
    ! Moreover, no data is copied and the vector is cleared if bclear=.TRUE.
    IF ((rvector%NEQ == 0) .OR.&
        (rvector%h_Ddata == ST_NOHANDLE)) THEN
      
      ! At first, copy all 'local' data.
      rvector = rvectorTemplate

      ! Then allocate a new array for the content in the same data
      ! format as the template vector.
      CALL storage_getsize(rvectorTemplate%h_Ddata, isize)
            
      IF (bclear) THEN
        CALL storage_new ('lsyssc_resizeVectorIndirect','ScalarVector',isize,&
            rvectorTemplate%cdataType, rvector%h_Ddata, ST_NEWBLOCK_ZERO)
      ELSE
        CALL storage_new ('lsyssc_resizeVectorIndirect','ScalarVector',isize,&
            rvectorTemplate%cdataType, rvector%h_Ddata, ST_NEWBLOCK_NOINIT)
      END IF

    ELSE
      
      ! Check if vectors are compatible
      IF ((rvector%cdataType /= rvectorTemplate%cdataType) .OR.&
          (rvector%NVAR /= rvectorTemplate%NVAR)) THEN
        PRINT *, "lsyssc_resizeVectorIndirect: Vectors are incompatible!"
        CALL sys_halt()
      END IF

      ! Get NEQMAX from template vector
      CALL storage_getsize(rvectorTemplate%h_Ddata, iNEQMAX)

      ! Resize vector directly
      CALL lsyssc_resizeVectorDirect(rvector, rvectorTemplate%NEQ, bclear, bcopy, iNEQMAX)

    END IF
    
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_resizeVectorIndMat (rmatrix, rvector, bclear, bcopy)

!<description>
    ! Resizes the vector structure so that it is compatible with the 
    ! scalar template matrix. If the vector does not exist, then it
    ! is created according to the matrix structure.
    !
    ! If it already exists, then it is resized accordingly so that
    ! its size coincides with the number of columns of the matrix.
!</description>

!<input>

    ! Scalar template matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix

    ! Whether to fill the vector with zero initially
    LOGICAL, INTENT(IN)                        :: bclear

    ! OPTIONAL: Whether to copy the content of the vector to the resized one
    LOGICAL, INTENT(IN), OPTIONAL              :: bcopy

!</input>

!<inputoutput>

    ! Scalar vector structure
    TYPE(t_vectorScalar), INTENT(INOUT)         :: rvector

!</inputoutput>

!</subroutine>

    ! Check, if vector is a copy of another (possibly larger) vector
    IF (rvector%bisCopy) THEN
      PRINT *, "lsyssc_resizeVectorIndirect: A copied vector cannot be resized!"
      CALL sys_halt()
    END IF
    
    ! Does the vector exist?
    IF (rvector%NEQ == 0 .OR.&
        rvector%h_Ddata == ST_NOHANDLE) THEN

      ! Create new vector
      CALL lsyssc_createVecIndMat(rmatrix, rvector, bclear)

    ELSE
      
      ! Check if vector/matrix are compatible
      IF (rvector%NVAR /= rmatrix%NVAR) THEN
        PRINT *, "lsyssc_resizeVectorIndMat: Vector/Matrix incompatible!"
        CALL sys_halt()
      END IF

      ! Resize vector directly
      CALL lsyssc_resizeVectorDirect(rvector, rmatrix%NCOLS, bclear, bcopy)
      
    END IF
    
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_resizeMatrixDirect (rmatrix, NEQ, NCOLS, NA, bclear, bcopy, &
      NEQMAX, NCOLSMAX, NAMAX, bforce)

!<description>
  ! Resizes the matrix structure to the new values NEQ, NCOLS and NA which are
  ! given explicitely.
  !
  ! IF NEQ, NCOLS, NA is smaller than the real memory allocated for the 
  ! specific handle, then only the values are adjusted and no reallocation
  ! takes place. Otherwise, the handles are reallocated accordingly. 
  !
  ! If the parameter bclear=.TRUE. the complete matrix is cleared.
  ! This is done both in case the matrix is reallocated physically or not.
  !
  ! If the parameter bcopy=.TRUE. the content of the existing matrix is 
  ! copied. This is only required, if the matrix needs to be reallocated
  ! physically. Otherwise, the data still exists in the memory.
  !
  ! Remark: The parameter bclear has higher priority than the parameter bcopy.
  ! That is, if both parameters are given, then no data is copied if the
  ! matrix shoud be cleared afterwards ;-)
  !
  ! If the optional parameters NEQMAX, NCOLSMAX or NAMAX are specified, then 
  ! memory is allocated for these values rather than NEQ, NCOLS or NA.
  !
  ! By default, a matrix cannot be resized, if it is copied (even partially)
  ! from another matrix. This is to prevent inconsistent data, i.e., the
  ! copied matrix is resized correctly but the original matrix still has
  ! old values NEQ, NCOLS, and NA. However, if the optional parameter
  ! bforce=.TRUE. then a copied matrix is resized virtually. That means,
  ! the atomic data NEQ, NCOLS, and NA are set to the new values if the
  ! corresponding handles for the structure and data are associated
  ! to memory blocks of sufficient size.
  ! Usability: Suppose you have one template matrix and a copy thereof.
  ! If you resize the template matrix, then the copied matrix still has its
  ! old dimensions whereas all of its pointers already point to larger
  ! memory blocks. Hence, it suffices to adjust, e.g., NEQ, NCOLS, and NA.
!</description>

!<input>

    ! Desired number of equations
    INTEGER(PREC_VECIDX), INTENT(IN)           :: NEQ

    ! Desired number of columns
    INTEGER(PREC_VECIDX), INTENT(IN)           :: NCOLS

    ! Desired number of elements
    INTEGER(PREC_MATIDX), INTENT(IN)           :: NA

    ! Whether to fill the matrix with zero initially
    LOGICAL, INTENT(IN)                        :: bclear

    ! OPTIONAL: Whether to copy the content of the matrix to the resized one
    LOGICAL, INTENT(IN), OPTIONAL              :: bcopy

    ! OPTIONAL: Maximum number of equations
    INTEGER(PREC_VECIDX), INTENT(IN), OPTIONAL :: NEQMAX

    ! OPTIONAL: Maximum number of columns
    INTEGER(PREC_VECIDX), INTENT(IN), OPTIONAL :: NCOLSMAX

    ! OPTIONAL: Maximum number of elements
    INTEGER(PREC_MATIDX), INTENT(IN), OPTIONAL :: NAMAX

    ! OPTIONAL: Wether to enforce resize even if matrix is copied from another matrix
    LOGICAL, INTENT(IN), OPTIONAL              :: bforce

!</input>

!<inputoutput>

    ! Scalar matrix structure
    TYPE(t_matrixScalar), INTENT(INOUT)        :: rmatrix

!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX) :: iNA,isize,isizeNew
    INTEGER(PREC_VECIDX) :: iNEQ,iNCOLS
    LOGICAL :: bdocopy,bdoresize,btransposed

    ! Check if resize should be forced
    bdoresize = .FALSE.
    IF (PRESENT(bforce)) bdoresize=bforce

    ! Check, if matrix is not a copy of another matrix or if resize is to be enforced
    IF ((IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) /= 0 .OR.&
         IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY)   /= 0) .AND.&
         .NOT.bdoresize) THEN
      PRINT *, "lsyssc_resizeMatrixDirect: A copied matrix can only be resized if&
          & this is forced explicitely!"
      CALL sys_halt()
    END IF

    ! Check, if matrix has been initialized before.
    IF (rmatrix%NEQ == 0 .OR. rmatrix%NCOLS == 0 .OR. rmatrix%NA == 0) THEN
      PRINT *, "lsyssc_resizeMatrixDirect: A matrix can only be resized &
          & if it has been created correctly!"
      CALL sys_halt()
    END IF

    ! Set working dimensions
    iNA    = MAX(0,NA)
    IF (PRESENT(NAMAX))    iNA    = MAX(iNA,NAMAX)
    iNEQ   = MAX(0,NEQ)
    IF (PRESENT(NEQMAX))   iNEQ   = MAX(iNEQ,NEQMAX)
    iNCOLS = MAX(0,NCOLS)
    IF (PRESENT(NCOLSMAX)) iNCOLS = MAX(iNCOLS,NCOLSMAX)  

    ! Set copy/clear attributes
    bdocopy = (.NOT.bclear)
    IF (PRESENT(bcopy)) bdocopy = (bdocopy .AND. bcopy)

    ! Set transposed indicator
    btransposed = (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) == 0)

    ! If the matrix should be cleared, then the sorting strategy (if any)
    ! can be ignored and reset. Otherwise, the matrix needs to be unsorted
    ! prior to copying some part of it. Afterwards, no sorting strategy is
    ! available in any case.
    IF (bdocopy .AND. rmatrix%isortStrategy > 0) THEN
      CALL lsyssc_unsortMatrix(rmatrix,.TRUE.)
    END IF

    ! Reset sorting strategy, there is none
    rmatrix%isortStrategy      = 0
    rmatrix%h_isortPermutation = ST_NOHANDLE

    ! Update NEQ, NCOLS and NA.
    rmatrix%NA    = NA
    rmatrix%NEQ   = NEQ
    rmatrix%NCOLS = NCOLS

    ! Now, the matrix has been virtually resized, that is, it already states
    ! the new dimensions. In order to resize the matrix physically, its structure
    ! and/or data need to be (re-)allocated).
    ! Remark: If the matrix was a copy of another matrix and if resize was
    ! not forced, then this routine would have stopped earlier. Hence, we can
    ! be sure, that resize is admissible in general. In other words, the structure/
    ! data is resized unconditionally if it is not copied from another matrix.
    ! If it is a copy of another matrix, then we check if corresponding handles
    ! are associated to large enough memory blocks. If not, then we stop with 
    ! an error. In any case, copied parts of the matrix are not cleared.
    
    ! What kind of matrix are we?
    SELECT CASE(rmatrix%cmatrixFormat)
      
    CASE (LSYSSC_MATRIX1,LSYSSC_MATRIXD)
      
      ! For these matrices, only the data array needs to be resized since
      ! there is no actual matrix structure.
      
      ! Are we copy or not?
      IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) == 1) THEN

        ! Check if handle coincides with matrix dimensions
        CALL storage_getsize(rmatrix%h_DA,isize)
        IF (rmatrix%NA > isize) THEN
          PRINT *, "lsyssc_resizeMatrixDirect: Dimensions of copied matrix mismatch!"
          CALL sys_halt()
        END IF

      ELSE
      
        ! Do we really have to reallocate the matrix data physically?
        IF (rmatrix%h_DA /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrix%h_DA,isize)
          IF (rmatrix%NA > isize) THEN
            
            ! Yes, so adopt the new size or reserve some extra memory if required.
            isize = iNA
            
            ! Reallocate the memory
            CALL storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_DA,&
                ST_NEWBLOCK_NOINIT, bdocopy)
            
          ELSEIF (PRESENT(NAMAX)) THEN
            
            ! The available memory suffices for the matrix. Let's check if the user
            ! suplied a new upper limit which makes it mandatory to "shrink" the
            ! allocated memory. Note that memory for at least NA matrix entries
            ! is allocated.
            IF (isize > iNA) THEN
              
              ! Set new size
              isize = iNA
              
              IF (isize == 0) THEN
                ! If nothing is left, then the matrix can also be released.
                CALL lsyssc_releaseMatrix(rmatrix)
                RETURN
              ELSE
                CALL storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_DA,&
                    ST_NEWBLOCK_NOINIT, bdocopy)
              END IF
            END IF
          END IF
          
          ! Should the matrix be cleared?
          IF (bclear) CALL storage_clear(rmatrix%h_DA)
        END IF
      END IF
      
    CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
      
      ! First, resize the structure of the matrix, i.e., KLD, KCOL and 
      ! possibly KDIAGONAL for matrices stored in matrix format 9. Here,
      ! it is necessary to distinguish between virtually transposed matrices
      ! and non-transposed ones.
      
      ! Are we copy or not?
      IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) == 1) THEN
        
        ! Check if handle coincides with matrix dimensions
        CALL storage_getsize(rmatrix%h_Kld,isize)
        IF ((rmatrix%NEQ+1 > isize) .OR. btransposed .AND. (rmatrix%NCOLS+1 > isize)) THEN
          PRINT *, "lsyssc_resizeMatrixDirect: Dimensions of copied matrix mismatch!"
          CALL sys_halt()
        END IF
        
        CALL storage_getsize(rmatrix%h_Kcol,isize)
        IF (rmatrix%NA > isize) THEN
          PRINT *, "lsyssc_resizeMatrixDirect: Dimensions of copied matrix mismatch!"
          CALL sys_halt()
        END IF
        
        IF ((rmatrix%cmatrixFormat == LSYSSC_MATRIX9) .OR.&
            (rmatrix%cmatrixFormat == LSYSSC_MATRIX9INTL)) THEN
          CALL storage_getsize(rmatrix%h_Kdiagonal,isize)
          IF ((rmatrix%NEQ > isize) .OR. btransposed .AND. (rmatrix%NCOLS > isize)) THEN
            PRINT *, "lsyssc_resizeMatrixDirect: Dimensions of copied matrix mismatch!"
            CALL sys_halt()
          END IF
        END IF
        
      ELSE   ! The structure of the matrix is not a copy of another matrix
      
        ! Do we really have to reallocate the array KLD?
        IF (rmatrix%h_Kld /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrix%h_Kld,isize)
          IF ((rmatrix%NEQ+1 > isize) .OR. btransposed .AND. (rmatrix%NCOLS+1 > isize)) THEN
            
            ! Yes, so adopt the new size or reserve some extra memory if required.
            isize = MERGE(iNCOLS+1,iNEQ+1,btransposed)
            
            ! Reallocate the memory
            CALL storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Kld,&
                ST_NEWBLOCK_NOINIT, bdocopy)
            
          ELSEIF (PRESENT(NEQMAX) .OR. btransposed .AND. PRESENT(NCOLSMAX)) THEN
            
            isizeNew = MERGE(iNCOLS+1,iNEQ+1,btransposed)

            ! The available memory suffices for the matrix. Let's check if the user
            ! suplied a new upper limit which makes it mandatory to "shrink" the
            ! allocated memory.
            IF (isize > isizeNew) THEN
              
              IF (isizeNew == 0) THEN
                ! If nothing is left, then the matrix can also be released.
                CALL lsyssc_releaseMatrix(rmatrix)
                RETURN
              ELSE
                CALL storage_realloc('lsyssc_resizeMatrixDirect', isizeNew, rmatrix%h_Kld,&
                    ST_NEWBLOCK_NOINIT, bdocopy)
              END IF
            END IF
          END IF
          
          ! Should the matrix be cleared?
          IF (bclear) CALL storage_clear(rmatrix%h_Kld)
        END IF
        
        ! Do we really have to reallocate the array KCOL?
        IF (rmatrix%h_Kcol /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrix%h_Kcol,isize)
          IF (rmatrix%NA > isize) THEN
            
            ! Yes, so adopt the new size or reserve some extra memory if required.
            isize = iNA
            
            ! Reallocate the memory
            CALL storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Kcol,&
                ST_NEWBLOCK_NOINIT, bdocopy)
            
          ELSEIF (PRESENT(NAMAX)) THEN
            
            ! The available memory suffices for the matrix. Let's check if the user
            ! suplied a new upper limit which makes it mandatory to "shrink" the
            ! allocated memory.
            IF (isize > iNA) THEN
              
              ! Set new size
              isize = iNA
              
              IF (isize == 0) THEN
                ! If nothing is left, then the matrix can also be released.
                CALL lsyssc_releaseMatrix(rmatrix)
                RETURN
              ELSE
                CALL storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Kcol,&
                    ST_NEWBLOCK_NOINIT, bdocopy)
              END IF
            END IF
          END IF
          
          ! Should the matrix be cleared?
          IF (bclear) CALL storage_clear(rmatrix%h_Kcol)
        END IF
        
        ! If the matrix is stored in format 9, then the diagonal array must be resized.
        IF ((rmatrix%cmatrixFormat == LSYSSC_MATRIX9) .OR.&
            (rmatrix%cmatrixFormat == LSYSSC_MATRIX9INTL)) THEN
          
          ! Do we really have to reallocate the array KCOL?
          IF (rmatrix%h_Kdiagonal /= ST_NOHANDLE) THEN
            CALL storage_getsize(rmatrix%h_Kdiagonal,isize)
            IF ((rmatrix%NEQ > isize) .OR. btransposed .AND. (rmatrix%NCOLS > isize)) THEN
              
              ! Yes, so adopt the new size or reserve some extra memory if required.
              isize = MERGE(iNCOLS,iNEQ,btransposed)
              
              ! Reallocate the memory
              CALL storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Kdiagonal,&
                  ST_NEWBLOCK_NOINIT, bdocopy)
              
            ELSEIF (PRESENT(NEQMAX) .OR. btransposed .AND. PRESENT(NCOLSMAX)) THEN
              
              isizeNew = MERGE(iNCOLS,iNEQ,btransposed)

              ! The available memory suffices for the matrix. Let's check if the user
              ! suplied a new upper limit which makes it mandatory to "shrink" the
              ! allocated memory.
              IF (isize > isizeNew) THEN
                
                IF (isizeNew == 0) THEN
                  ! If nothing is left, then the matrix can also be released.
                  CALL lsyssc_releaseMatrix(rmatrix)
                  RETURN
                ELSE
                  CALL storage_realloc('lsyssc_resizeMatrixDirect', isizeNew, rmatrix%h_Kdiagonal,&
                      ST_NEWBLOCK_NOINIT, bdocopy)
                END IF
              END IF
            END IF
            
            ! Should the matrix be cleared?
            IF (bclear) CALL storage_clear(rmatrix%h_Kdiagonal)
          END IF
        END IF
        
      END IF
      
      ! Ok, the matrix structure has been resized according to the prescribed dimensions.
      ! Now, let us resize the data array of the matrix.
      
      ! Are we copy or not?
      IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) == 1) THEN

        ! Check if handle coincides with matrix dimensions
        CALL storage_getsize(rmatrix%h_DA,isize)

        ! What kind of interleave matrix are we (if any)?
        SELECT CASE(rmatrix%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIXUNDEFINED)
          IF (rmatrix%NA > isize) THEN
            PRINT *, "lsyssc_resizeMatrixDirect: Dimensions of copied matrix mismatch!"
            CALL sys_halt()
          END IF
          
        CASE (LSYSSC_MATRIX1)
          IF (rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR > isize) THEN
            PRINT *, "lsyssc_resizeMatrixDirect: Dimensions of copied matrix mismatch!"
            CALL sys_halt()
          END IF

        CASE (LSYSSC_MATRIXD)
          IF (rmatrix%NA*rmatrix%NVAR > isize) THEN
            PRINT *, "lsyssc_resizeMatrixDirect: Dimensions of copied matrix mismatch!"
            CALL sys_halt()
          END IF

        CASE DEFAULT
          PRINT *, "lsyssc_resizeMatrixDirect: Unsupported interleave matrix format!"
          CALL sys_halt()
        END SELECT

      ELSE   ! The content of the matrix is not a copy of another matrix

        ! What kind of interleave matrix are we (if any)?
        SELECT CASE(rmatrix%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIXUNDEFINED)
          
          ! Do we really have to reallocate the matrix physically?
          IF (rmatrix%h_DA /= ST_NOHANDLE) THEN
            CALL storage_getsize(rmatrix%h_DA,isize)
            IF (rmatrix%NA > isize) THEN
              
              ! Yes, so adopt the new size or reserve some extra memory if required.
              isize = iNA
              
              ! Reallocate the memory
              CALL storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_DA,&
                  ST_NEWBLOCK_NOINIT, bdocopy)
              
            ELSEIF (PRESENT(NAMAX)) THEN
              
              ! The available memory suffices for the matrix. Let's check if the user
              ! suplied a new upper limit which makes it mandatory to "shrink" the
              ! allocated memory. Note that memory for at least NA matrix entries
              ! is allocated.
              IF (isize > iNA) THEN
                
                ! Set new size
                isize = iNA
                
                IF (isize == 0) THEN
                  ! If nothing is left, then the matrix can also be released.
                  CALL lsyssc_releaseMatrix(rmatrix)
                  RETURN
                ELSE
                  CALL storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_DA,&
                      ST_NEWBLOCK_NOINIT, bdocopy)
                END IF
              END IF
            END IF
            
            ! Should the matrix be cleared?
            IF (bclear) CALL storage_clear(rmatrix%h_DA)
          END IF
          
        CASE (LSYSSC_MATRIX1)
          
          ! Do we really have to reallocate the matrix physically?
          IF (rmatrix%h_DA /= ST_NOHANDLE) THEN
            CALL storage_getsize(rmatrix%h_DA,isize)
            IF (rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR > isize) THEN
              
              ! Yes, so adopt the new size or reserve some extra memory if required.
              isize = iNA*rmatrix%NVAR*rmatrix%NVAR
              
              ! Reallocate the memory
              CALL storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_DA,&
                  ST_NEWBLOCK_NOINIT, bdocopy)
              
            ELSEIF (PRESENT(NAMAX)) THEN
              
              ! The available memory suffices for the matrix. Let's check if the user
              ! suplied a new upper limit which makes it mandatory to "shrink" the
              ! allocated memory. Note that memory for at least NA matrix entries
              ! is allocated.
              IF (isize > iNA*rmatrix%NVAR*rmatrix%NVAR) THEN
                
                ! Set new size
                isize = iNA*rmatrix%NVAR*rmatrix%NVAR
                
                IF (isize == 0) THEN
                  ! If nothing is left, then the matrix can also be released.
                  CALL lsyssc_releaseMatrix(rmatrix)
                  RETURN
                ELSE
                  CALL storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_DA,&
                      ST_NEWBLOCK_NOINIT, bdocopy)
                END IF
              END IF
            END IF
            
            ! Should the matrix be cleared?
            IF (bclear) CALL storage_clear(rmatrix%h_DA)
          END IF
          
        CASE (LSYSSC_MATRIXD)
          
          ! Do we really have to reallocate the matrix physically?
          IF (rmatrix%h_DA /= ST_NOHANDLE) THEN
            CALL storage_getsize(rmatrix%h_DA,isize)
            IF (rmatrix%NA*rmatrix%NVAR > isize) THEN
              
              ! Yes, so adopt the new size or reserve some extra memory if required.
              isize = iNA*rmatrix%NVAR
              
              ! Reallocate the memory
              CALL storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_DA,&
                  ST_NEWBLOCK_NOINIT, bdocopy)
              
            ELSEIF (PRESENT(NAMAX)) THEN
              
              ! The available memory suffices for the matrix. Let's check if the user
              ! suplied a new upper limit which makes it mandatory to "shrink" the
              ! allocated memory. Note that memory for at least NA matrix entries
              ! is allocated.
              IF (isize > iNA*rmatrix%NVAR) THEN
                
                ! Set new size
                isize = iNA*rmatrix%NVAR
                
                IF (isize == 0) THEN
                  ! If nothing is left, then the matrix can also be released.
                  CALL lsyssc_releaseMatrix(rmatrix)
                  RETURN
                ELSE
                  CALL storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_DA,&
                      ST_NEWBLOCK_NOINIT, bdocopy)
                END IF
              END IF
            END IF
            
            ! Should the matrix be cleared?
            IF (bclear) CALL storage_clear(rmatrix%h_DA)
          END IF
          
        CASE DEFAULT
          PRINT *, "lsyssc_resizeMatrixDirect: Unsupported interleave matrix format!"
          CALL sys_halt()
        END SELECT
        
      END IF
      
    CASE DEFAULT
      PRINT *, "lsyssc_resizeMatrixDirect: Unsupported matrix format!"
      CALL sys_halt()
    END SELECT
  END SUBROUTINE lsyssc_resizeMatrixDirect

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_resizeMatrixIndirect(rmatrix, rmatrixTemplate, bclear, bcopy, bforce)

!<description>
  ! Resizes the matrix so that it exhibits the same memory layout as the
  ! template matrix. Note that this subroutine can only be used if the
  ! template matrix has the same internal structure, e.g., data type, NVAR
  ! as the matrix to be resized.
  !
  ! If the parameter bclear=.TRUE. the complete matrix is cleared.
  ! This is done both in case the matrix is reallocated physically or not.
  !
  ! If the parameter bcopy=.TRUE. the content of the existing matrix is 
  ! copied. This is only required, if the matrix needs to be reallocated
  ! physically. Otherwise, the data still exists in the memory.
  !
  ! Remark: The parameter bclear has higher priority than the parameter bcopy.
  ! That is, if both parameters are given, then no data is copied if the
  ! matrix shoud be cleared afterwards ;-)
  !
  ! By default, a matrix cannot be resized, if it is copied (even partially)
  ! from another matrix. This is to prevent inconsistent data, i.e., the
  ! copied matrix is resized correctly but the original matrix still has
  ! old values NEQ, NCOLS, and NA. However, if the optional parameter
  ! bforce=.TRUE. then a copied matrix is resized virtually. That means,
  ! the atomic data NEQ, NCOLS, and NA are set to the new values if the
  ! corresponding handles for the structure and data are associated
  ! to memory blocks of sufficient size.
  ! Usability: Suppose you have one template matrix and a copy thereof.
  ! If you resize the template matrix, then the copied matrix still has its
  ! old dimensions whereas all of its pointers already point to larger
  ! memory blocks. Hence, it suffices to adjust, e.g., NEQ, NCOLS, and NA.
!</description>

!<input>

    ! Scalar template matrix
    TYPE(t_matrixScalar), INTENT(IN)           :: rmatrixTemplate

    ! Whether to fill the vector with zero initially
    LOGICAL, INTENT(IN)                        :: bclear

    ! OPTIONAL: Whether to copy the content of the vector to the resized one
    LOGICAL, INTENT(IN), OPTIONAL              :: bcopy

    ! OPTIONAL: Whether to enforce resize even if matrix is copied from another matrix
    LOGICAL, INTENT(IN), OPTIONAL              :: bforce

!</input>

!<inputoutput>

    ! Scalar matrix structure
    TYPE(t_matrixScalar), INTENT(INOUT)        :: rmatrix

!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX) :: isize,isizeTmp
    LOGICAL :: bdocopy,bdoresize

    ! Check if resize should be forced
    bdoresize = .FALSE.
    IF (PRESENT(bforce)) bdoresize=bforce

    ! Check, if matrix is not a copy of another matrix
    IF ((IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) /= 0 .OR.&
         IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY)   /= 0) .AND.&
         .NOT.bdoresize) THEN
      PRINT *, "lsyssc_resizeMatrixDirect: A copied matrix can only be resized if&
          & this is forced explicitely!"
      CALL sys_halt()
    END IF

    ! Check, if matrix has been initialized before.
    IF (rmatrix%cmatrixFormat == LSYSSC_MATRIXUNDEFINED) THEN

      ! At first, copy all 'local' data.
      rmatrix = rmatrixTemplate

      ! Then allocate memory for the matrix structure and matrix data.
      IF (bclear) THEN
        IF (rmatrixTemplate%h_Kld /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrixTemplate%h_Kld, isize)
          CALL storage_new ('lsyssc_resizeMatrixIndirect', 'h_Kld', isize,&
              ST_INT, rmatrix%h_Kld, ST_NEWBLOCK_ZERO)
        END IF
        IF (rmatrixTemplate%h_Kcol /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrixTemplate%h_Kcol, isize)
          CALL storage_new ('lsyssc_resizeMatrixIndirect', 'h_Kcol', isize,&
              ST_INT, rmatrix%h_Kcol, ST_NEWBLOCK_ZERO)
        END IF
        IF (rmatrixTemplate%h_Kdiagonal /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrixTemplate%h_Kdiagonal, isize)
          CALL storage_new ('lsyssc_resizeMatrixIndirect', 'h_Kdiagonal', isize,&
              ST_INT, rmatrix%h_Kdiagonal, ST_NEWBLOCK_ZERO)
        END IF
        IF (rmatrixTemplate%h_DA /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrixTemplate%h_DA, isize)
          CALL storage_new ('lsyssc_resizeMatrixIndirect', 'h_DA', isize,&
              rmatrixTemplate%cdataType, rmatrix%h_DA, ST_NEWBLOCK_ZERO)
        END IF

      ELSE   ! Matrix should not be cleared

        IF (rmatrixTemplate%h_Kld /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrixTemplate%h_Kld, isize)
          CALL storage_new ('lsyssc_resizeMatrixIndirect', 'h_Kld', isize,&
              ST_INT, rmatrix%h_Kld, ST_NEWBLOCK_NOINIT)
        END IF
        IF (rmatrixTemplate%h_Kcol /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrixTemplate%h_Kcol, isize)
          CALL storage_new ('lsyssc_resizeMatrixIndirect', 'h_Kcol', isize,&
              ST_INT, rmatrix%h_Kcol, ST_NEWBLOCK_NOINIT)
        END IF
        IF (rmatrixTemplate%h_Kdiagonal /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrixTemplate%h_Kdiagonal, isize)
          CALL storage_new ('lsyssc_resizeMatrixIndirect', 'h_Kdiagonal', isize,&
              ST_INT, rmatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT)
        END IF
        IF (rmatrixTemplate%h_DA /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrixTemplate%h_DA, isize)
          CALL storage_new ('lsyssc_resizeMatrixIndirect', 'h_DA', isize,&
              rmatrixTemplate%cdataType, rmatrix%h_DA, ST_NEWBLOCK_NOINIT)
        END IF
      END IF
      
    ELSE

      ! The matrix has been initialized before.

      ! Check if matrices are compatible except for the dimensions
      IF ((rmatrix%cmatrixFormat           /= rmatrixTemplate%cmatrixFormat) .OR.&
          (rmatrix%cinterleavematrixFormat /= rmatrixTemplate%cinterleavematrixFormat) .OR.&
          (rmatrix%NVAR                    /= rmatrixTemplate%NVAR) .OR.&
          (rmatrix%cdataType               /= rmatrixTemplate%cdataType) .OR.&
          IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) /=&
          IAND(rmatrixTemplate%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED)) THEN
        PRINT *, "lsyssc_resizeMatrixDirect: Matrices are incompatible!"
        CALL sys_halt()
      END IF
      
      ! Set copy/clear attributes
      bdocopy = (.NOT.bclear)
      IF (PRESENT(bcopy)) bdocopy = (bdocopy .AND. bcopy)

      ! If the matrix should be cleared, then the sorting strategy (if any)
      ! can be ignored and reset. Otherwise, the matrix needs to be unsorted
      ! prior to copying some part of it. Afterwards, no sorting strategy is
      ! available in any case.
      IF (bdocopy .AND. rmatrix%isortStrategy > 0) THEN
        CALL lsyssc_unsortMatrix(rmatrix,.TRUE.)
      END IF
      
      ! Reset sorting strategy, there is none
      rmatrix%isortStrategy      = 0
      rmatrix%h_isortPermutation = ST_NOHANDLE
      
      ! Update NA
      rmatrix%NA = rmatrixTemplate%NA

      ! Are we copy or not?
      IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) == 1) THEN

        ! Check if handle coincides with matrix dimensions
        CALL storage_getsize(rmatrix%h_DA,isize)
        IF (rmatrix%NA > isize) THEN
          PRINT *, "lsyssc_resizeMatrixIndirect: Dimensions of copied matrix mismatch!"
          CALL sys_halt()
        END IF

      ELSE
        
        ! Do we have to reallocate the handle?
        IF (rmatrix%h_DA /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrix%h_DA, isize)
          IF (rmatrix%NA > isize) THEN
            
            ! Yes, we have to reallocate the handle. 
            isize = rmatrix%NA

            ! Also consider the size of the template matrix.
            IF (rmatrixTemplate%h_DA /= ST_NOHANDLE) THEN
              CALL storage_getsize(rmatrixTemplate%h_DA, isizeTmp)
              isize = MAX(0,isize,isizeTmp)
            END IF
            
            ! Reallocate the memory
            CALL storage_realloc('lsyssc_resizeMatrixIndirect', isize, rmatrix%h_DA,&
                ST_NEWBLOCK_NOINIT, bdocopy)
          END IF
        END IF
      END IF
      
      ! Update NEQ and NCOLS.
      rmatrix%NEQ   = rmatrixTemplate%NEQ
      rmatrix%NCOLS = rmatrixTemplate%NCOLS

      ! Are we copy or not?
      IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) == 1) THEN
        
        ! Check if handle coincides with matrix simensions
        IF (rmatrix%h_Kld /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrix%h_Kld,isize)
          IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) == 1 .AND.&
              rmatrix%NCOLS+1 > isize .OR. rmatrix%NEQ+1 > isize) THEN
            PRINT *, "lsyssc_resizeMatrixIndirect: Dimensions of copied matrix mismatch!"
            CALL sys_halt()
          END IF
        END IF

        ! Check if handle coincides with matrix simensions
        IF (rmatrix%h_Kcol /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrix%h_Kcol, isize)
          IF (rmatrix%NA > isize) THEN
            PRINT *, "lsyssc_resizeMatrixIndirect: Dimensions of copied matrix mismatch!"
            CALL sys_halt()
          END IF
        END IF

        ! Check if handle coincides with matrix simensions
        IF (rmatrix%h_Kdiagonal /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrix%h_Kdiagonal, isize)
          IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) == 1 .AND.&
              rmatrix%NCOLS > isize .OR. rmatrix%NEQ > isize) THEN
            PRINT *, "lsyssc_resizeMatrixIndirect: Dimensions of copied matrix mismatch!"
            CALL sys_halt()
          END IF
        END IF
        
      ELSE  ! Matrix is no copy
        
        ! Do we have to reallocate the handle?
        IF (rmatrix%h_Kld /= ST_NOHANDLE) THEN
          
          ! Do we process a virtually transposed matrix?
          IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) == 1) THEN
            CALL storage_getsize(rmatrix%h_Kld, isize)
            IF (rmatrix%NCOLS+1 > isize) THEN
              
              ! Yes, we have to reallocate the handle. 
              isize = rmatrix%NCOLS+1

              ! Also consider the size of the template matrix.
              IF (rmatrixTemplate%h_Kld /= ST_NOHANDLE) THEN
                CALL storage_getsize(rmatrixTemplate%h_Kld, isizeTmp)
                isize = MAX(0,isize,isizeTmp)
              END IF
              
              ! Reallocate the memory
              CALL storage_realloc('lsyssc_resizeMatrixIndirect', isize, rmatrix%h_Kld,&
                  ST_NEWBLOCK_NOINIT, bdocopy)
            END IF
          ELSE
            CALL storage_getsize(rmatrix%h_Kld, isize)
            IF (rmatrix%NEQ+1 > isize) THEN
              
              ! Yes, we have to reallocate the handle. 
              isize = rmatrix%NEQ+1

              ! Also consider the size of the template matrix.
              IF (rmatrixTemplate%h_Kld /= ST_NOHANDLE) THEN
                CALL storage_getsize(rmatrixTemplate%h_Kld, isizeTmp)
                isize = MAX(0,isize,isizeTmp)
              END IF

              ! Reallocate the memory
              CALL storage_realloc('lsyssc_resizeMatrixIndirect', isize, rmatrix%h_Kld,&
                  ST_NEWBLOCK_NOINIT, bdocopy)
            END IF
          END IF
        END IF

        ! Do we have to reallocate the handle?
        IF (rmatrix%h_Kcol /= ST_NOHANDLE) THEN
          CALL storage_getsize(rmatrix%h_Kcol, isize)
          IF (rmatrix%NA > isize) THEN
            
            ! Yes, we have to reallocate the handle. 
            isize = rmatrix%NA

            ! Also consider the size of the template matrix.
            IF (rmatrixTemplate%h_Kcol /= ST_NOHANDLE) THEN
              CALL storage_getsize(rmatrixTemplate%h_Kcol, isizeTmp)
              isize = MAX(0,isize,isizeTmp)
            END IF
            
            ! Reallocate the memory
            CALL storage_realloc('lsyssc_resizeMatrixIndirect', isize, rmatrix%h_Kcol,&
                ST_NEWBLOCK_NOINIT, bdocopy)
          END IF
        END IF

        ! Do we have to reallocate the handle?
        IF (rmatrix%h_Kdiagonal /= ST_NOHANDLE) THEN

          ! Do we process a virtually transposed matrix?
          IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) == 1) THEN
            CALL storage_getsize(rmatrix%h_Kdiagonal, isize)
            IF (rmatrix%NCOLS > isize) THEN
              
              ! Yes, we have to reallocate the handle. 
              isize = rmatrix%NCOLS

              ! Also consider the size of the template matrix.
              IF (rmatrixTemplate%h_Kdiagonal /= ST_NOHANDLE) THEN
                CALL storage_getsize(rmatrixTemplate%h_Kdiagonal, isizeTmp)
                isize = MAX(0,isize,isizeTmp)
              END IF
              
              ! Reallocate the memory
              CALL storage_realloc('lsyssc_resizeMatrixIndirect', isize, rmatrix%h_Kdiagonal,&
                  ST_NEWBLOCK_NOINIT, bdocopy)
            END IF
          ELSE
            CALL storage_getsize(rmatrix%h_Kdiagonal, isize)
            IF (rmatrix%NEQ > isize) THEN
              
              ! Yes, we have to reallocate the handle. 
              isize = rmatrix%NEQ

              ! Also consider the size of the template matrix.
              IF (rmatrixTemplate%h_Kdiagonal /= ST_NOHANDLE) THEN
                CALL storage_getsize(rmatrixTemplate%h_Kdiagonal, isizeTmp)
                isize = MAX(0,isize,isizeTmp)
              END IF
              
              ! Reallocate the memory
              CALL storage_realloc('lsyssc_resizeMatrixIndirect', isize, rmatrix%h_Kdiagonal,&
                  ST_NEWBLOCK_NOINIT, bdocopy)
            END IF
          END IF
          
        END IF
      END IF
    END IF

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
  
  ! Vectors must be compatible!
  CALL lsyssc_isVectorCompatible (rx,ry)
  
  ! Is there data at all?
  res = 0.0_DP
  
  IF ((rx%NEQ*rx%NVAR .EQ. 0) .OR. &
      (ry%NEQ*ry%NVAR .EQ. 0) .OR. &
      (rx%NEQ .NE. rx%NEQ) .OR. &
      (rx%NVAR .NE. ry%NVAR)) THEN
    PRINT *,'Error in lsyssc_scalarProduct: Vector dimensions wrong!'
    CALL sys_halt()
  END IF
  
  IF (rx%cdataType .NE. ry%cdataType) THEN
    PRINT *,'lsyssc_scalarProduct: Data types different!'
    CALL sys_halt()
  END IF
  
  ! Take care of the data type before doing a scalar product!
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
     
    ! Get the data arrays
    CALL lsyssc_getbase_double (rx,p_Ddata1dp)
    CALL lsyssc_getbase_double (ry,p_Ddata2dp)
    
    ! Perform the scalar product
    res=lalg_scalarProductDble(p_Ddata1dp,p_Ddata2dp)
    
  CASE (ST_SINGLE)

    ! Get the data arrays
    CALL lsyssc_getbase_single (rx,p_Fdata1dp)
    CALL lsyssc_getbase_single (ry,p_Fdata2dp)
    
    ! Perform the scalar product
    res=lalg_scalarProductSngl(p_Fdata1dp,p_Fdata2dp)
    
  CASE DEFAULT
    PRINT *,'lsyssc_scalarProduct: Not supported precision combination'
    CALL sys_halt()
  END SELECT
  
  ! Return the scalar product, finish
  lsyssc_scalarProduct = res

  END FUNCTION

  !****************************************************************************
!<subroutine>
  
  REAL(DP) FUNCTION lsyssc_scalarProductMatVec (rx, ry)
  
!<description>
  ! Calculates a scalar product of two vectors, whereby the first
    ! vector is a diagonal matrix.
!</description>
  
!<input>
  ! First vector given as diagonal matrix
  TYPE(t_MatrixScalar), INTENT(IN)                  :: rx

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
  
  ! Matrix must be diagonal matrix
  IF (rx%cmatrixFormat .NE. LSYSSC_MATRIXD) THEN
    PRINT *, 'lsyssc_scalarProductMatVec: Matrix must be diagonal matrix!'
    CALL sys_halt()
  END IF

  ! Vectors must be compatible!
  CALL lsyssc_isMatrixVectorCompatible (ry,rx,.FALSE.)
  
  ! Is there data at all?
  res = 0.0_DP
  
  IF ((rx%NEQ*rx%NVAR .EQ. 0) .OR. &
      (ry%NEQ*ry%NVAR .EQ. 0) .OR. &
      (rx%NEQ .NE. rx%NEQ) .OR. &
      (rx%NVAR .NE. ry%NVAR)) THEN
    PRINT *,'Error in lsyssc_scalarProductMatVec: Vector dimensions wrong!'
    CALL sys_halt()
  END IF
  
  IF (rx%cdataType .NE. ry%cdataType) THEN
    PRINT *,'lsyssc_scalarProductMatVec: Data types different!'
    CALL sys_halt()
  END IF

  ! Take care of the data type before doing a scalar product!
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
     
    ! Get the data arrays
    CALL lsyssc_getbase_double (rx,p_Ddata1dp)
    CALL lsyssc_getbase_double (ry,p_Ddata2dp)
    
    ! Perform the scalar product
    res=lalg_scalarProductDble(p_Ddata1dp,p_Ddata2dp)
    
  CASE (ST_SINGLE)

    ! Get the data arrays
    CALL lsyssc_getbase_single (rx,p_Fdata1dp)
    CALL lsyssc_getbase_single (ry,p_Fdata2dp)
    
    ! Perform the scalar product
    res=lalg_scalarProductSngl(p_Fdata1dp,p_Fdata2dp)
    
  CASE DEFAULT
    PRINT *,'lsyssc_scalarProduct: Not supported precision combination'
    CALL sys_halt()
  END SELECT
  
  ! Return the scalar product, finish
  lsyssc_scalarProductMatVec = res

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

    ! If the scale factor is =0, we have nothing to do.
    IF (rmatrix%dscaleFactor .EQ. 0.0_DP) RETURN
    
    ! Vectors must be compatible to the matrix.
    CALL lsyssc_isMatrixCompatible (rx,rmatrix,.FALSE.)
    CALL lsyssc_isMatrixCompatible (ry,rmatrix,.TRUE.)
    
    ! rx and ry must have at least the same data type!
    IF (rx%cdataType .NE. ry%cdataType) THEN
      PRINT *,'MV with different data types for rx and ry not supported!'
      CALL sys_halt()
    END IF
    
    ! Up to now, all matrix types use h_Da. So if that's not associated,
    ! there is for sure an error!
    IF (rmatrix%h_Da .EQ. ST_NOHANDLE) THEN
      PRINT *,'lsyssc_scalarMatVec: Matrix has no data!'
      CALL sys_halt()
    END IF
    
    ! Handle the scaling factor by multiplication of cx with dscaleFactor.
    !
    ! Now, which matrix format do we have?
    
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
            CALL lsyssc_LAX79doubledouble (&
                rmatrix,rx,ry,cx*rmatrix%dscaleFactor,cy)
          
          CASE DEFAULT
            PRINT *,'Only double precision vectors supported for now in MV!'
            CALL sys_halt()
            
          END SELECT
          
        CASE DEFAULT
          PRINT *,'Only double precision matrices supported for now in MV!'
          CALL sys_halt()
        END SELECT
        
      CASE (LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9INTL)
        
        ! Take care of the precision of the entries
        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
          ! Format 7 and Format 9 multiplication
          SELECT CASE (rx%cdataType)
          
          CASE (ST_DOUBLE)
            ! double precision matrix, double precision vectors
            SELECT CASE(rmatrix%cinterleavematrixFormat)
            CASE (LSYSSC_MATRIX1)
              CALL lsyssc_LAX79INTL1doubledouble (&
                  rmatrix,rx,ry,cx*rmatrix%dscaleFactor,cy)

            CASE (LSYSSC_MATRIXD)
              CALL lsyssc_LAX79INTLDdoubledouble (&
                  rmatrix,rx,ry,cx*rmatrix%dscaleFactor,cy)

            CASE DEFAULT
              PRINT *, 'Invalid interleave matrix format!'
              CALL sys_halt()
            END SELECT
          
          CASE DEFAULT
            PRINT *,'Only double precision vectors supported for now in MV!'
            CALL sys_halt()
            
          END SELECT
          
        CASE DEFAULT
          PRINT *,'Only double precision matrices supported for now in MV!'
          CALL sys_halt()
        END SELECT

      CASE (LSYSSC_MATRIXD)
      
        ! Take care of the precision of the entries
        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
          ! Format D multiplication
          SELECT CASE (rx%cdataType)
          
          CASE (ST_DOUBLE)
            ! double precision matrix, double precision vectors
            CALL lsyssc_LATXDdoubledouble (&
                rmatrix,rx,ry,cx*rmatrix%dscaleFactor,cy)
          
          CASE DEFAULT
            PRINT *,'Only double precision vectors supported for now in MV!'
            CALL sys_halt()
            
          END SELECT
          
        CASE DEFAULT
          PRINT *,'Only double precision matrices supported for now in MV!'
          CALL sys_halt()
        END SELECT

      CASE DEFAULT
        PRINT *,'Unknown matrix format in MV-multiplication!'
        CALL sys_halt()
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
            CALL lsyssc_LTX9doubledouble (&
                rmatrix,rx,ry,cx*rmatrix%dscaleFactor,cy)
          
          CASE DEFAULT
            PRINT *,'Only double precision vectors supported for now in MV!'
            CALL sys_halt()
            
          END SELECT
          
        CASE DEFAULT
          PRINT *,'Only double precision matrices supported for now in MV!'
          CALL sys_halt()
        END SELECT
        
      CASE (LSYSSC_MATRIX7)
      
        ! Take care of the precision of the entries
        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
          ! Format 7 multiplication
          SELECT CASE (rx%cdataType)
          
          CASE (ST_DOUBLE)
            ! double precision matrix, double precision vectors
            CALL lsyssc_LTX7doubledouble (&
                rmatrix,rx,ry,cx*rmatrix%dscaleFactor,cy)
          
          CASE DEFAULT
            PRINT *,'Only double precision vectors supported for now in MV!'
            CALL sys_halt()
            
          END SELECT
          
        CASE DEFAULT
          PRINT *,'Only double precision matrices supported for now in MV!'
          CALL sys_halt()
        END SELECT
        
      CASE (LSYSSC_MATRIXD)
      
        ! Take care of the precision of the entries
        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
          ! Format D multiplication
          SELECT CASE (rx%cdataType)
          
          CASE (ST_DOUBLE)
            ! double precision matrix, double precision vectors
            CALL lsyssc_LATXDdoubledouble (&
                rmatrix,rx,ry,cx*rmatrix%dscaleFactor,cy)
          
          CASE DEFAULT
            PRINT *,'Only double precision vectors supported for now in MV!'
            CALL sys_halt()
            
          END SELECT
          
        CASE DEFAULT
          PRINT *,'Only double precision matrices supported for now in MV!'
          CALL sys_halt()
        END SELECT

      CASE DEFAULT
        PRINT *,'Unknown matrix format in MV-multiplication!'
        CALL sys_halt()
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
    INTEGER(PREC_MATIDX) :: ia
    INTEGER(PREC_VECIDX) :: irow,icol
    REAL(DP) :: dtmp
    INTEGER(PREC_VECIDX) :: NEQ

      ! Get the matrix
      CALL lsyssc_getbase_double (rmatrix,p_DA)
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      
      ! Get NEQ - from the matrix, not from the vector!
      NEQ = rmatrix%NEQ

      ! Get the vectors
      CALL lsyssc_getbase_double (rx,p_Dx)
      CALL lsyssc_getbase_double (ry,p_Dy)
      
      ! By commenting in the following two lines,
      ! one can gain a slight speedup in the matrix vector
      ! multiplication by avoiding some checks and by using
      ! arrays with undefined length.
      ! Makes a difference of 57 to 40 MFLOP/s, but should only
      ! be used in RELEASE-mode!
      !
      ! #ifdef RELEASE
      ! CALL lsyssc_qLAX79doubledouble (p_DA,p_Kcol,p_Kld,p_Dx,p_Dy,cx,cy,NEQ)
      ! RETURN
      ! #endif
      
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
       
!$omp parallel do default(shared) private(irow,icol,ia)
          DO irow=1,NEQ
            ia   = p_Kld(irow)
            icol = p_Kcol(ia)
            p_Dy(irow) = p_Dx(icol) * p_DA(ia)
          END DO
!$omp end parallel do

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
          
!$omp parallel do default(shared) private(irow,icol,ia)
          DO irow=1,NEQ
            ia   = p_Kld(irow)
            icol = p_Kcol(ia)
            p_Dy(irow) = p_Dx(icol)*p_DA(ia) + p_Dy(irow) 
          END DO
!$omp end parallel do
          
        ENDIF
        
        ! Multiply the rest of rx with the matrix and add it to ry:
        
!$omp parallel do default(shared) private(irow,icol,ia)
        DO irow=1,NEQ
          DO ia = p_Kld(irow)+1,p_Kld(irow+1)-1
            icol = p_Kcol(ia)
            p_Dy(irow) = p_Dy(irow) + p_DA(ia)*p_Dx(icol)
          END DO
        END DO
!$omp end parallel do
        
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
    ! Format 7 and Format 9 multiplication
    ! double precision matrix,
    ! double precision vectors
    ! 'Quick' method avoiding some checks, thus being faster.
    
    SUBROUTINE lsyssc_qLAX79doubledouble (DA,Kcol,Kld,Dx,Dy,cx,cy,NEQ)

    ! The matrix
    INTEGER(PREC_VECIDX), DIMENSION(*), INTENT(IN) :: KCOL
    INTEGER(PREC_MATIDX), DIMENSION(*), INTENT(IN) :: KLD
    REAL(DP), DIMENSION(*), INTENT(IN) :: DA
    
    ! Size of the vectors
    INTEGER(PREC_VECIDX) :: NEQ

    ! The vectors
    REAL(DP), DIMENSION(*), INTENT(IN) :: Dx
    REAL(DP), DIMENSION(*), INTENT(INOUT) :: Dy
    
    ! Multiplication factors for the vectors.
    REAL(DP) :: cx,cy
    
    INTEGER(PREC_MATIDX) :: ia
    INTEGER(PREC_VECIDX) :: irow,icol
    REAL(DP) :: dtmp

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
       
!$omp parallel do default(shared) private(irow,icol,ia)
          DO irow=1,NEQ
            ia   = Kld(irow)
            icol = Kcol(ia)
            Dy(irow) = Dx(icol) * DA(ia)
          END DO
!$omp end parallel do

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
            CALL lalg_scaleVectorDble(Dy(1:NEQ),dtmp)
          END IF
          
          ! Multiply the first entry in each line of the matrix with the
          ! corresponding entry in rx and add it to the (scaled) ry.
          
!$omp parallel do default(shared) private(irow,icol,ia)
          DO irow=1,NEQ
            ia   = Kld(irow)
            icol = Kcol(ia)
            Dy(irow) = Dx(icol)*DA(ia) + Dy(irow) 
          END DO
!$omp end parallel do
          
        ENDIF
        
        ! Multiply the rest of rx with the matrix and add it to ry:
        
!$omp parallel do default(shared) private(irow,icol,ia)
        DO irow=1,NEQ
          DO ia = Kld(irow)+1,Kld(irow+1)-1
            icol = Kcol(ia)
            Dy(irow) = Dy(irow) + DA(ia)*Dx(icol)
          END DO
        END DO
!$omp end parallel do
        
        ! Scale by cx, finish.
        
        IF (cx .NE. 1.0_DP) THEN
          CALL lalg_scaleVectorDble (Dy(1:NEQ),cx)
        END IF
        
      ELSE 
        ! cx = 0. The formula is just a scaling of the vector ry!
        CALL lalg_scaleVectorDble(Dy(1:NEQ),cy)
      ENDIF
   
    END SUBROUTINE
    
    !**************************************************************
    ! Format 7 and Format 9 full interleaved multiplication
    ! double precision matrix,
    ! double precision vectors
    
    SUBROUTINE lsyssc_LAX79INTL1doubledouble (rmatrix,rx,ry,cx,cy)

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
    INTEGER(PREC_VECIDX) :: irow,icol,ia
    REAL(DP) :: dtmp
    INTEGER(PREC_VECIDX) :: NEQ
    INTEGER :: ivar,jvar
    INTEGER :: NVAR

      ! Get the matrix
      CALL lsyssc_getbase_double (rmatrix,p_DA)
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      
      ! Get NEQ - from the matrix, not from the vector!
      NEQ = rmatrix%NEQ

      ! Get NVAR - from the matrix, not from the vector!
      NVAR = rmatrix%NVAR

      ! Check if vectors are compatible
      IF (rx%NVAR /= NVAR .OR. ry%NVAR /= NVAR) THEN
        PRINT *, "lsyssc_LAX79INTL1doubledouble: Matrix/Vector is incompatible!"
        CALL sys_halt()
      END IF

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
!$omp parallel do default(shared) private(irow,icol,ia,ivar,jvar,dtmp)
          DO irow=1,NEQ
            ia   = p_Kld(irow)
            icol = p_Kcol(ia)
            
            ! Here, we compute
            !   y(ivar,irow) = SUM_jvar ( A(ivar,jvar,ia)*x(jvar,icol) )
            DO ivar=1,NVAR
              dtmp = 0
              DO jvar=1,NVAR
                dtmp = dtmp + p_Dx(NVAR*(icol-1)+jvar)&
                    * p_DA(NVAR*NVAR*(ia-1)+NVAR*(jvar-1)+ivar)
              END DO
              p_Dy(NVAR*(irow-1)+ivar) = dtmp
            END DO
          END DO
!$omp end parallel do

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
!$omp parallel do default(shared) private(irow,icol,ia,ivar,jvar,dtmp)
          DO irow=1,NEQ
            ia   = p_Kld(irow)
            icol = p_Kcol(ia)

            ! Here, we compute
            !   y(ivar,irow) = y(ivar,irow) + SUM_jvar ( A(ivar,jvar,ia)*x(jvar,icol) )
            DO ivar=1,NVAR
              dtmp = 0
              DO jvar=1,NVAR
                dtmp = dtmp + p_Dx(NVAR*(icol-1)+jvar)&
                    * p_DA(NVAR*NVAR*(ia-1)+NVAR*(jvar-1)+ivar)
              END DO
              p_Dy(NVAR*(irow-1)+ivar) = dtmp + p_Dy(NVAR*(irow-1)+ivar)
            END DO
          END DO
!$omp end parallel do
          
        ENDIF
        
        ! Multiply the rest of rx with the matrix and add it to ry:
!$omp parallel do default(shared) private(irow,icol,ia,ivar,jvar,dtmp)
          DO irow=1,NEQ
            DO ia = p_Kld(irow)+1,p_Kld(irow+1)-1
              icol = p_Kcol(ia)
              
              ! Here, we compute
              !   y(ivar,irow) = y(ivar,irow) + SUM_jvar ( A(ivar,jvar,ia)*x(jvar,icol) )
              DO ivar=1,NVAR
                dtmp = 0
                DO jvar=1,NVAR
                  dtmp = dtmp + p_Dx(NVAR*(icol-1)+jvar)&
                      * p_DA(NVAR*NVAR*(ia-1)+NVAR*(jvar-1)+ivar)
                END DO
                p_Dy(NVAR*(irow-1)+ivar) = dtmp + p_Dy(NVAR*(irow-1)+ivar)
              END DO
            END DO
          END DO
!$omp end parallel do
        
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
    ! Format 7 and Format 9 diagonal interleaved multiplication
    ! double precision matrix,
    ! double precision vectors
    
    SUBROUTINE lsyssc_LAX79INTLDdoubledouble (rmatrix,rx,ry,cx,cy)

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
    INTEGER(PREC_VECIDX) :: irow,icol,ia
    REAL(DP) :: dtmp
    INTEGER(PREC_VECIDX) :: NEQ
    INTEGER :: ivar,jvar
    INTEGER :: NVAR

      ! Get the matrix
      CALL lsyssc_getbase_double (rmatrix,p_DA)
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      
      ! Get NEQ - from the matrix, not from the vector!
      NEQ = rmatrix%NEQ

      ! Get NVAR - from the matrix, not from the vector!
      NVAR = rmatrix%NVAR

      ! Check if vectors are compatible
      IF (rx%NVAR /= NVAR .OR. ry%NVAR /= NVAR) THEN
        PRINT *, "lsyssc_LAX79INTLDdoubledouble: Matrix/Vector is incompatible!"
        CALL sys_halt()
      END IF

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
!$omp parallel do default(shared) private(irow,icol,ia,ivar)
          DO irow=1,NEQ
            ia   = p_Kld(irow)
            icol = p_Kcol(ia)

            DO ivar=1,NVAR
              p_Dy(NVAR*(irow-1)+ivar) = p_Dx(NVAR*(icol-1)+ivar)&
                  * p_DA(NVAR*(ia-1)+ivar)
            END DO
          END DO
!$omp end parallel do

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
!$omp parallel do default(shared) private(irow,icol,ia,ivar)
          DO irow=1,NEQ
            ia   = p_Kld(irow)
            icol = p_Kcol(ia)
            
            DO ivar=1,NVAR
              p_Dy(NVAR*(irow-1)+ivar) = p_Dx(NVAR*(icol-1)+ivar)&
                  * p_DA(NVAR*(ia-1)+ivar)&
                  + p_Dy(NVAR*(irow-1)+ivar)
            END DO
          END DO
!$omp end parallel do
          
        ENDIF
        
        ! Multiply the rest of rx with the matrix and add it to ry:
!$omp parallel do default(shared) private(irow,icol,ia,ivar)
          DO irow=1,NEQ
            DO ia = p_Kld(irow)+1,p_Kld(irow+1)-1
              icol = p_Kcol(ia)
              
              DO ivar=1,NVAR
                p_Dy(NVAR*(irow-1)+ivar) = p_Dx(NVAR*(icol-1)+ivar)&
                    * p_DA(NVAR*(ia-1)+ivar)&
                    + p_Dy(NVAR*(irow-1)+ivar)
              END DO
            END DO
          END DO
!$omp end parallel do
           
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
    ! Format D (diagonal matrix) multiplication
    ! double precision matrix,
    ! double precision vectors
    ! As we have  diagonal matrix, this is used for both, MV and
    ! tranposed MV.
    
    SUBROUTINE lsyssc_LATXDdoubledouble (rmatrix,rx,ry,cx,cy)

    ! Save arguments as above - given as parameters as some compilers
    ! might have problems with scoping units...
    TYPE(t_matrixScalar), INTENT(IN)                  :: rmatrix
    TYPE(t_vectorScalar), INTENT(IN)                  :: rx
    REAL(DP), INTENT(IN)                              :: cx
    REAL(DP), INTENT(IN)                              :: cy
    TYPE(t_vectorScalar), INTENT(INOUT)               :: ry

    REAL(DP), DIMENSION(:), POINTER :: p_DA, p_Dx, p_Dy
    REAL(DP) :: dtmp
    INTEGER(PREC_VECIDX) :: irow,NEQ
    INTEGER :: ivar,NVAR

      ! Get the matrix - it's an 1D array
      CALL lsyssc_getbase_double (rmatrix,p_DA)
      
      ! Get NEQ - from the matrix, not from the vector!
      NEQ = rmatrix%NEQ

      ! Get NVAR - from the vector not from the matrix!
      NVAR = rx%NVAR

      IF (NVAR /= ry%NVAR) THEN
        PRINT *, "Internal structure of vectors is not compatible!"
        CALL sys_halt()
      END IF

      ! Get the vectors
      CALL lsyssc_getbase_double (rx,p_Dx)
      CALL lsyssc_getbase_double (ry,p_Dy)

      IF (NVAR == 1) THEN
      
        ! Perform the multiplication
        IF (cx .NE. 0.0_DP) THEN
          
          IF (cy .EQ. 0.0_DP) THEN
            
            ! cy = 0. Multiply cx*A with X and write to Y.
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(irow)
            DO irow = 1,NEQ
              p_Dy(irow) = cx*p_Da(irow)*p_Dx(irow)
            END DO
!$omp end parallel do
          
          ELSE
        
            ! Full multiplication: cx*A*X + cy*Y
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(irow)
            DO irow = 1,NEQ
              p_Dy(irow) = cy*p_Dy(irow) + cx*p_Da(irow)*p_Dx(irow) 
            END DO
!$omp end parallel do
        
          END IF
        
        ELSE 
          
          ! cx = 0. The formula is just a scaling of the vector ry!
          CALL lalg_scaleVectorDble(p_Dy,cy)
          
        ENDIF
      
      ELSE

        ! Perform the multiplication
        IF (cx .NE. 0.0_DP) THEN
          
          IF (cy .EQ. 0.0_DP) THEN
            
            ! cy = 0. Multiply cx*A with X and write to Y.
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(irow)
            DO irow = 1,NEQ
              DO ivar=1,NVAR
                p_Dy(NVAR*(irow-1)+ivar) = cx*p_Da(irow)*p_Dx(NVAR*(irow-1)+ivar)
              END DO
            END DO
!$omp end parallel do
          
          ELSE
        
            ! Full multiplication: cx*A*X + cy*Y
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(irow)
            DO irow = 1,NEQ
              DO ivar=1,NVAR
                p_Dy(NVAR*(irow-1)+ivar) = cy*p_Dy(NVAR*(irow-1)+ivar) + cx*p_Da(irow)*p_Dx(NVAR*(irow-1)+ivar) 
              END DO
            END DO
!$omp end parallel do
        
          END IF
        
        ELSE 
          
          ! cx = 0. The formula is just a scaling of the vector ry!
          CALL lalg_scaleVectorDble(p_Dy,cy)
          
        ENDIF

      END IF

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
      CALL lsyssc_getbase_double (rmatrix,p_DA)
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      
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
          
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(irow)
          DO irow=1,NEQ
            p_Dy(irow) = p_Dx(irow)*p_DA(p_Kld(irow)) 
          END DO
!$omp end parallel do
          
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
          
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(irow)
          DO irow=1,NEQ
            p_Dy(irow) = p_Dy(irow) + p_Dx(irow)*p_DA(p_Kld(irow)) 
          END DO
!$omp end parallel do
          
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
    INTEGER(PREC_MATIDX) :: ia
    INTEGER(PREC_VECIDX) :: irow,icol
    REAL(DP) :: dtmp
    INTEGER(PREC_VECIDX) :: NEQ

      ! Get the matrix
      CALL lsyssc_getbase_double (rmatrix,p_DA)
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      
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
          DO ia = p_Kld(irow),p_Kld(irow+1)-1
            icol = p_Kcol(ia)
            p_Dy(icol) = p_Dy(icol) + p_Dx(irow)*p_DA(ia)
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
  
  SUBROUTINE lsyssc_duplicateVector (rx,ry,cdupStructure,cdupContent)
  
!<description>
  ! Duplicates an existing vector: ry := rx.
  ! Creates a new vector ry based on a template vector rx.
  ! Duplicating a vector does not necessarily mean that new memory is
  ! allocated and the vector entries are copied to that. The two flags
  ! cdupStructure and cdupContent decide on how to set up ry.
  ! Depending on their setting, it's possible to copy only the handle
  ! of such dynamic information, so that both vectors share the same
  ! information.
  !
  ! We follow the following convention:
  !  Structure = NEQ, sorting permutation, discretisation-related information.
  !  Content   = Enties in the vector.
  !
  ! Remark 1: There is never memory allocated on the heap for the sorting
  !  permutation. A vector is never the 'owner' of a permutation, i.e.
  !  does not maintain it. Therefore, copying a permutation means
  !  copying the corresponding handle. The application must keep track
  !  of the permutations.
  ! Remark 2: The vector is never resorted! If rx is sorted while
  !  ry is unsorted and cdupStructure=LSYSSC_DUP_IGNORE,
  !  the state of the new vector is undefined as no reference to the
  !  sorting strategy is transferred! The data is simply copied.
!</description>
  
!<input>
  ! Vector to copy
  TYPE(t_vectorScalar), INTENT(IN)                :: rx

  ! Duplication flag that decides on how to set up the structure
  ! of ry. Not all flags are possible!
  ! One of the LSYSSC_DUP_xxxx flags:
  ! LSYSSC_DUP_IGNORE     : Ignore the structure of rx
  ! LSYSSC_DUP_COPY or
  ! LSYSSC_DUP_COPYOVERWRITE or
  ! LSYSSC_DUP_TEMPLATE   : Structural data is copied from rx
  !   to ry (NEQ, sorting strategy, pointer to discretisation structure).
  INTEGER, INTENT(IN)                            :: cdupStructure
  
  ! Duplication flag that decides on how to set up the content
  ! of ry. Not all flags are possible!
  ! One of the LSYSSC_DUP_xxxx flags:
  ! LSYSSC_DUP_IGNORE     : Ignore the content of rx.
  ! LSYSSC_DUP_REMOVE     : Removes any existing content from 
  !   ry if there is any. Releases memory if necessary.
  ! LSYSSC_DUP_DISMISS    : Removes any existing content from 
  !   ry if there is any. No memory is released, handles are simply
  !   dismissed.
  ! LSYSSC_DUP_SHARE      : ry receives the same handles for
  !   content data as rx and therefore shares the same content.
  !   If necessary, the old data in ry is released.
  ! LSYSSC_DUP_COPY       : ry gets a copy of the content of rx.
  !   If necessary, new memory is allocated for the content.
  !   If ry already contains allocated memory that belongs
  !   to that vector, content/structure data is simply copied from rx
  !   into that.
  !   Note that this respects the ownership! I.e. if the destination vector is not
  !   the owner of the content/structure data arrays, new memory is allocated to
  !   prevent the actual owner from getting destroyed!
  ! LSYSSC_DUP_COPYOVERWRITE:   The destination vector gets a copy of the content 
  !   of rx. If necessary, new memory is allocated.
  !   If the destination vector already contains allocated memory, content
  !   data is simply copied from rx into that.
  !   The ownership of the content/data arrays is not respected, i.e. if the
  !   destination vector is not the owner, the actual owner of the data arrays is
  !   modified, too!
  ! LSYSSC_DUP_ASIS       : Duplicate by ownership. If the content of 
  !   rx belongs to rx, ry gets a copy
  !   of the content; new memory is allocated if necessary (the same as 
  !   LSYSSC_DUP_COPY). If the content of rx belongs to another
  !   vector than rx, ry receives the same handles as
  !   rx and is therefore a third vector sharing the same content
  !   (the same as LSYSSC_DUP_SHARE, so rx, ry and the 
  !   other vector have the same content).
  ! LSYSSC_DUP_EMPTY      : New memory is allocated for the content in the
  !   same size in rsourceMatrix but no data is copied;
  !   the arrays are left uninitialised.
  INTEGER, INTENT(IN)                            :: cdupContent
!</input>

!<inputoutput>
  ! The new vector which will be a copy of roldvector
  TYPE(t_vectorScalar), INTENT(INOUT)               :: ry
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER(PREC_VECIDX) :: ioffset,cdataType
    LOGICAL :: bisCopy
    INTEGER :: h_Ddata
    
    ! First the structure
    
    IF (cdupStructure .NE. LSYSSC_DUP_IGNORE) THEN
      ! Standard case: Copy the whole vector information
      !
      ! First, make a backup of some crucial data so that it does not
      ! get destroyed.
      h_Ddata = ry%h_Ddata
      cdataType = ry%cdataType
      bisCopy = ry%bisCopy
      ioffset = ry%iidxFirstEntry
      
      ! Then transfer all structural information of rx to ry.
      ! This automatically makes both vectors compatible to each other.
      ry = rx
      
      ! Restore crucial data
      ry%h_Ddata = h_Ddata
      ry%bisCopy = bisCopy
      ry%cdataType = cdataType
      ry%iidxFirstEntry = ioffset 
    END IF

    ! Now the content.

    SELECT CASE (cdupContent)
    CASE (LSYSSC_DUP_IGNORE) 
      ! Nothing to do
    
    CASE (LSYSSC_DUP_REMOVE)
      ! Release vector data
      IF ((.NOT. ry%bisCopy) .AND. (ry%h_Ddata .NE. ST_NOHANDLE)) &
        CALL storage_free (ry%h_Ddata)
        
      ry%h_Ddata = ST_NOHANDLE
      ry%bisCopy = .FALSE.
      
    CASE (LSYSSC_DUP_DISMISS)
      ! Dismiss data
      ry%h_Ddata = ST_NOHANDLE
      ry%bisCopy = .FALSE.
    
    CASE (LSYSSC_DUP_SHARE)
      ! Share information. Release memory if necessary.
      IF ((.NOT. ry%bisCopy) .AND. (ry%h_Ddata .NE. ST_NOHANDLE)) &
        CALL storage_free (ry%h_Ddata)
        
      ry%h_Ddata = rx%h_Ddata
      ry%bisCopy = .TRUE.
      
    CASE (LSYSSC_DUP_COPY)
      ! Copy information. If necessary, allocate new data -- by setting
      ! h_Ddata to ST_NOHANDLE before calling storage_copy.
      IF (ry%bisCopy) THEN
        ry%h_Ddata = ST_NOHANDLE
        ry%bisCopy = .FALSE.
      END IF

      IF (ry%h_Ddata .EQ. ST_NOHANDLE) &
        CALL allocDestinationVector (rx,ry)
      
      CALL copy_data (rx,ry)
    
    CASE (LSYSSC_DUP_COPYOVERWRITE)
      ! Copy information, regardless of whether ry is the owner or not.
      ! If no memory is allocated, allocate new memory.
      
      IF (ry%h_Ddata .EQ. ST_NOHANDLE) &
        CALL allocDestinationVector (rx,ry)
      
      CALL copy_data (rx,ry)
    
    CASE (LSYSSC_DUP_ASIS)
    
      ! Copy by ownership. This is either LSYSSC_COPY or LSYSSC_SHARE,
      ! depending on whether rx is the owner of the data or not.
      IF (rx%bisCopy) THEN
        ! rx shares it's data and thus ry will also.
        ry%h_Ddata = rx%h_Ddata
        ry%bisCopy = .TRUE.
      ELSE
        ! The data belongs to rx and thus it must also belong to ry --
        ! so copy it.
        ry%h_Ddata = ST_NOHANDLE
        ry%bisCopy = .FALSE.
        CALL allocDestinationVector (rx,ry)
        CALL copy_data (rx,ry)
      END IF
      
    CASE (LSYSSC_DUP_EMPTY)
    
      ! Allocate new memory if ry is empty. Don't initialise.
      ! If ry contains data, we don't have to do anything.
      IF (ry%h_Ddata .EQ. ST_NOHANDLE) THEN
        CALL allocDestinationVector (rx,ry)
      END IF
    
    CASE DEFAULT
    
      PRINT *,'lsyssc_duplicateVector: cdupContent unknown!'
      CALL sys_halt()
    
    END SELECT
    
  CONTAINS
  
    SUBROUTINE allocDestinationVector (rx,ry)
    
    ! Allocates new memory ry in the size of rx. The memory is not initialised.
    
    TYPE(t_vectorScalar), INTENT(IN) :: rx
    TYPE(t_vectorScalar), INTENT(INOUT) :: ry

      ! local variables
      INTEGER(PREC_VECIDX) :: NEQ,NVAR
    
      NEQ = rx%NEQ
      NVAR= rx%NVAR
      ry%cdataType = rx%cdataType
      CALL storage_new ('lsyssc_duplicateVector','ry',NEQ*NVAR,&
                        rx%cdataType, ry%h_Ddata, &
                        ST_NEWBLOCK_NOINIT)
    
    END SUBROUTINE
    
    ! ---------------------------------------------------------------
    
    SUBROUTINE copy_data (rx,ry)
    
    ! Copies the content of rx to ry. Takes care of the data type
    
    TYPE(t_vectorScalar), INTENT(IN) :: rx
    TYPE(t_vectorScalar), INTENT(INOUT) :: ry
    
      ! local variables
      REAL(DP), DIMENSION(:), POINTER :: p_Dsource,p_Ddest
      REAL(SP), DIMENSION(:), POINTER :: p_Fsource,p_Fdest    
    
      ! Take care of the data type!
      SELECT CASE (rx%cdataType)
      CASE (ST_DOUBLE)
        ! Get the pointer and scale the whole data array.
        CALL lsyssc_getbase_double(rx,p_Dsource)
        CALL lsyssc_getbase_double(ry,p_Ddest)
        CALL lalg_copyVectorDble (p_Dsource,p_Ddest)  

      CASE (ST_SINGLE)
        ! Get the pointer and scale the whole data array.
        CALL lsyssc_getbase_single(rx,p_Fsource)
        CALL lsyssc_getbase_single(ry,p_Fdest)
        CALL lalg_copyVectorSngl (p_Fsource,p_Fdest)  

      CASE DEFAULT
        CALL output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateVector')
        CALL sys_halt()
      END SELECT

    END SUBROUTINE    
    
  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_duplicateMatrix (rsourceMatrix,rdestMatrix,&
                                     cdupStructure, cdupContent)
  
!<description>
  ! Duplicates an existing matrix, creates a new matrix rdestMatrix based
  ! on a template matrix rsourceMatrix.
  ! Duplicating a matrix does not necessarily mean that new memory is
  ! allocated and the matrix entries are copied to that. The two flags
  ! cdupStructure and cdupContent decide on how to set up rdestMatrix.
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
  !   If rdestMatrix already contains allocated memory that belongs
  !   to that matrix, content/structure data is simply copied from rsourceMatrix
  !   into that.
  !   Note that this respects the ownership! I.e. if the destination matrix is not
  !   the owner of the content/structure data arrays, new memory is allocated to
  !   prevent the actual owner from getting destroyed!
  ! LSYSSC_DUP_COPYOVERWRITE:   The destination matrix gets a copy of the content 
  !   of rsourceMatrix. If necessary, new memory is allocated.
  !   If the destination matrix  already contains allocated memory, content/structure 
  !   data is simply copied from rsourceMatrix into that.
  !   The ownership of the content/data arrays is not respected, i.e. if the
  !   destination matrix is not the owner, the actual owner of the data arrays is
  !   modified, too!
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
  ! LSYSSC_DUP_TEMPLATE   : Copies static structural information about the
  !   structure (NEQ, NCOLS,...) to the destination matrix. Dynamic information
  !   is removed from the destination matrix, all handles are reset.
  INTEGER, INTENT(IN)                            :: cdupStructure
  
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
  !   If rdestMatrix already contains allocated memory that belongs
  !   to that matrix, content/structure data is simply copied from rsourceMatrix
  !   into that.
  !   Note that this respects the ownership! I.e. if the destination matrix is not
  !   the owner of the content/structure data arrays, new memory is allocated to
  !   prevent the actual owner from getting destroyed!
  ! LSYSSC_DUP_COPYOVERWRITE:   The destination matrix gets a copy of the content 
  !   of rsourceMatrix. If necessary, new memory is allocated.
  !   If the destination matrix  already contains allocated memory, content/structure 
  !   data is simply copied from rsourceMatrix into that.
  !   The ownership of the content/data arrays is not respected, i.e. if the
  !   destination matrix is not the owner, the actual owner of the data arrays is
  !   modified, too!
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
  ! LSYSSC_DUP_TEMPLATE   : Copies static structural information about the
  !   structure (NEQ, NCOLS,...) to the destination matrix. Dynamic information
  !   is removed from the destination matrix, all handles are reset.
  INTEGER, INTENT(IN)                            :: cdupContent
!</input>

!<output>
  ! Destination matrix.
  TYPE(t_matrixScalar), INTENT(INOUT)            :: rdestMatrix
!</output>  

!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX) :: isize
  
    ! What do we have to do for the structure?
    SELECT CASE (cdupStructure)
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
      ! Remove the structure - if there is any.
      CALL removeStructure(rdestMatrix, .TRUE.)

      ! Share the structure between rsourceMatrix and rdestMatrix
      CALL shareStructure(rsourceMatrix, rdestMatrix)   
      
    CASE (LSYSSC_DUP_COPY)
      ! Copy the structure of rsourceMatrix to rdestMatrix
      CALL copyStructure(rsourceMatrix, rdestMatrix, .FALSE.)  
      
    CASE (LSYSSC_DUP_COPYOVERWRITE)
      ! Copy the structure of rsourceMatrix to rdestMatrix, don't respect ownership
      CALL copyStructure(rsourceMatrix, rdestMatrix, .TRUE.)  
      
    CASE (LSYSSC_DUP_ASIS)
      ! What's with the source matrix. Does the structure belong to it?
      IF (IAND(rsourceMatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .NE. 0) THEN
        ! Copy the structure
        CALL copyStructure(rsourceMatrix, rdestMatrix, .FALSE.)
      ELSE
        ! Remove the structure - if there is any.
        CALL removeStructure(rdestMatrix, .TRUE.)

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
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
        CALL storage_new ('lsyssc_duplicateMatrix', 'KCOL', &
            rsourceMatrix%NA, ST_INT, &
            rdestMatrix%h_Kcol, ST_NEWBLOCK_NOINIT)

        CALL storage_new ('lsyssc_duplicateMatrix', 'KLD', &
            rdestMatrix%NEQ+1_I32, ST_INT, &
            rdestMatrix%h_Kld, ST_NEWBLOCK_NOINIT)

        CALL storage_new ('lsyssc_duplicateMatrix', 'Kdiagonal', &
            rdestMatrix%NEQ, ST_INT, &
            rdestMatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT)
        
      CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
        CALL storage_new ('lsyssc_duplicateMatrix', 'KCOL', &
            rdestMatrix%NA, ST_INT, &
            rdestMatrix%h_Kcol, ST_NEWBLOCK_NOINIT)

        CALL storage_new ('lsyssc_duplicateMatrix', 'KLD', &
            rdestMatrix%NEQ+1_I32, ST_INT, &
            rdestMatrix%h_Kld, ST_NEWBLOCK_NOINIT)

      CASE (LSYSSC_MATRIXD)
        ! Nothing to do
      END SELECT
    END SELECT
    
    ! -----
    ! Ok, handling of the structure is finished. The data is handled similar.
    ! -----
    
    ! What do we have to do for the content?
    SELECT CASE (cdupContent)
    CASE (LSYSSC_DUP_IGNORE)
      ! Nothing
    CASE (LSYSSC_DUP_REMOVE)
      ! Remove the content - if there is any.
      CALL removeContent(rdestMatrix, .TRUE.)
      
    CASE (LSYSSC_DUP_DISMISS)
      ! Dismiss the structure - if there is any.
      CALL removeContent(rdestMatrix, .FALSE.)
      
    CASE (LSYSSC_DUP_SHARE)
      ! Remove the content - if there is any.
      CALL removeContent(rdestMatrix, .TRUE.)

      ! Share the structure between rsourceMatrix and rdestMatrix
      CALL shareContent(rsourceMatrix, rdestMatrix)   
      
    CASE (LSYSSC_DUP_COPY)
      ! Copy the structure of rsourceMatrix to rdestMatrix
      CALL copyContent(rsourceMatrix, rdestMatrix, .FALSE.)  
      
    CASE (LSYSSC_DUP_COPYOVERWRITE)
      ! Copy the structure of rsourceMatrix to rdestMatrix, don't respect ownership
      CALL copyContent(rsourceMatrix, rdestMatrix, .TRUE.)  
      
    CASE (LSYSSC_DUP_TEMPLATE)
      ! Remove the structure - if there is any.
      CALL removeContent(rdestMatrix, .TRUE.)
      
      ! Copy static content information
      CALL copyStaticContent(rsourceMatrix,rdestMatrix)

    CASE (LSYSSC_DUP_ASIS)
      ! What's with the source matrix. Does the structure belong to it?
      IF (IAND(rsourceMatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .NE. 0) THEN
        ! Copy the structure
        CALL copyContent(rsourceMatrix, rdestMatrix, .FALSE.)
      ELSE
        ! Remove the content - if there is any.
        CALL removeContent(rdestMatrix, .TRUE.)

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
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7,LSYSSC_MATRIXD)
        ! Create a new content array in the same data type as the original matrix
        CALL storage_new('lsyssc_duplicateMatrix', 'DA', rdestMatrix%NA, &
            rsourceMatrix%cdataType, rdestMatrix%h_DA, ST_NEWBLOCK_NOINIT)
      CASE (LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7INTL)
        ! Create a new content array in the same data type as the original matrix
        SELECT CASE(rsourceMatrix%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIX1)
          CALL storage_new('lsyssc_duplicateMatrix', 'DA', &
              INT(rdestMatrix%NA*rdestMatrix%NVAR*rdestMatrix%NVAR,I32), &
              rsourceMatrix%cdataType, rdestMatrix%h_DA, ST_NEWBLOCK_NOINIT)
        CASE (LSYSSC_MATRIXD)
          CALL storage_new('lsyssc_duplicateMatrix', 'DA', &
              INT(rdestMatrix%NA*rdestMatrix%NVAR,I32),rsourceMatrix%cdataType, &
              rdestMatrix%h_DA, ST_NEWBLOCK_NOINIT)
        CASE DEFAULT
          PRINT *, 'lsyssc_duplicateMatrix: wrong matrix format of interleaved matrix'
          CALL sys_halt()
        END SELECT
      END SELECT
    END SELECT
    
    ! -----
    ! Final check. Check if we destroyed the matrix. May only happen if the
    ! user on-purpose calls this routine two times with different source
    ! but the same destination matrix.
    ! -----
    
    SELECT CASE (rdestMatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
    
      ! Check length of DA
      IF (rdestMatrix%h_DA .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_DA,isize)
        SELECT CASE(rdestMatrix%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIX1)
          IF (isize .LT. rdestMatrix%NA * rdestMatrix%NVAR * rdestMatrix%NVAR) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA < length(DA)!'
            CALL sys_halt()
          END IF
        CASE (LSYSSC_MATRIXD)
          IF (isize .LT. rdestMatrix%NA * rdestMatrix%NVAR) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA < length(DA)!'
            CALL sys_halt()
          END IF
        CASE DEFAULT
          IF (isize .LT. rdestMatrix%NA) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA < length(DA)!'
            CALL sys_halt()
          END IF
        END SELECT
      END IF
      
      ! Check length of KCOL
      IF (rdestMatrix%h_Kcol .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_Kcol,isize)
        IF (isize .LT. rdestMatrix%NA) THEN
          PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA < length(KCOL)!'
          CALL sys_halt()
        END IF
      END IF

      ! Check length of KLD
      IF (rdestMatrix%h_Kld .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_Kld,isize)
        
        ! Be careful, matrix may be transposed.
        IF (IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .EQ. 0) THEN
          IF (isize .LT. rdestMatrix%NEQ+1) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NEQ+1 < length(KLD)!'
            CALL sys_halt()
          END IF
        ELSE
          IF (isize .LT. rdestMatrix%NCOLS+1) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NEQ+1 < length(KLD)!'
            CALL sys_halt()
          END IF
        END IF
      END IF
      
      ! Check length of Kdiagonal
      IF (rdestMatrix%h_Kdiagonal .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_Kdiagonal,isize)
        
        ! Be careful, matrix may be transposed.
        IF (IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .EQ. 0) THEN
          IF (isize .LT. rdestMatrix%NEQ) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NEQ+1 < length(Kdiag)!'
            CALL sys_halt()
          END IF
        ELSE
          IF (isize .LT. rdestMatrix%NCOLS) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NEQ+1 < length(Kdiag)!'
            CALL sys_halt()
          END IF
        END IF
      END IF
      
    CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
    
      ! Check length of DA
      IF (rdestMatrix%h_DA .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_DA,isize)
        SELECT CASE(rdestMatrix%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIX1)
          IF (isize .LT. rdestMatrix%NA * rdestMatrix%NVAR*rdestMatrix%NVAR) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA < length(DA)!'
            CALL sys_halt()
          END IF
        CASE (LSYSSC_MATRIXD)
          IF (isize .LT. rdestMatrix%NA * rdestMatrix%NVAR) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA < length(DA)!'
            CALL sys_halt()
          END IF
        CASE DEFAULT
          IF (isize .LT. rdestMatrix%NA) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA < length(DA)!'
            CALL sys_halt()
          END IF
        END SELECT
      END IF
      
      ! Check length of KCOL
      IF (rdestMatrix%h_Kcol .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_Kcol,isize)
        IF (isize .LT. rdestMatrix%NA) THEN
          PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA < length(KCOL)!'
          CALL sys_halt()
        END IF
      END IF

      ! Check length of KLD
      IF (rdestMatrix%h_Kld .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_Kld,isize)
        
        ! Be careful, matrix may be transposed.
        IF (IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .EQ. 0) THEN
          IF (isize .LT. rdestMatrix%NEQ+1) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NEQ+1 < length(KLD)!'
            CALL sys_halt()
          END IF
        ELSE
          IF (isize .LT. rdestMatrix%NCOLS+1) THEN
            PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NEQ+1 < length(KLD)!'
            CALL sys_halt()
          END IF
        END IF
      END IF
      
    CASE (LSYSSC_MATRIXD)
    
      ! Check length of DA
      IF (rdestMatrix%h_DA .NE. ST_NOHANDLE) THEN
        CALL storage_getsize (rdestMatrix%h_DA,isize)
        IF (isize .NE. rdestMatrix%NA) THEN
          PRINT *,'lsyssc_duplicateMatrix: Matrix destroyed; NA != length(DA)!'
          CALL sys_halt()
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
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
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
      CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
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
        
      CASE (LSYSSC_MATRIXD,LSYSSC_MATRIX1)
        ! Nothing to do

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
      rdestMatrix%NVAR   = rsourceMatrix%NVAR
      rdestMatrix%cmatrixFormat           = rsourceMatrix%cmatrixFormat
      rdestMatrix%cinterleavematrixFormat = rsourceMatrix%cinterleavematrixFormat
      rdestMatrix%isortStrategy           = rsourceMatrix%isortStrategy
      rdestMatrix%h_IsortPermutation      = rsourceMatrix%h_IsortPermutation
      
      ! Transfer all flags except the 'dup' flags
      iflag = IAND(rsourceMatrix%imatrixSpec,NOT(LSYSSC_MSPEC_ISCOPY))
      iflag2 = IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
      rdestMatrix%imatrixSpec = IOR(iflag,iflag2)

      ! Transfer discretisation-related information,
      rdestMatrix%p_rspatialDiscretisation => rsourceMatrix%p_rspatialDiscretisation

      ! Which source matrix do we have?  
      SELECT CASE (rsourceMatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
        rdestMatrix%h_Kcol = rsourceMatrix%h_Kcol
        rdestMatrix%h_Kld  = rsourceMatrix%h_Kld
        rdestMatrix%h_Kdiagonal = rsourceMatrix%h_Kdiagonal
        
      CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
        ! Overwrite structural data
        rdestMatrix%h_Kcol = rsourceMatrix%h_Kcol
        rdestMatrix%h_Kld  = rsourceMatrix%h_Kld
      
      CASE (LSYSSC_MATRIXD,LSYSSC_MATRIX1)
        ! Nothing to do

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
      rdestMatrix%NVAR   = rsourceMatrix%NVAR
      rdestMatrix%cmatrixFormat           = rsourceMatrix%cmatrixFormat
      rdestMatrix%cinterleavematrixFormat = rsourceMatrix%cinterleavematrixFormat
      rdestMatrix%isortStrategy           = rsourceMatrix%isortStrategy
      rdestMatrix%h_IsortPermutation      = rsourceMatrix%h_IsortPermutation

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
    
    SUBROUTINE copyStructure(rsourceMatrix, rdestMatrix, bignoreOwner)
    
    ! The source matrix 
    TYPE(t_matrixScalar), INTENT(IN) :: rsourceMatrix

    ! The destination matrix 
    TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
    
    ! Whether to respect the ownership of the arrays and allocate 
    ! memory automatically if a matrix is not the owner of the arrays.
    LOGICAL, INTENT(IN) :: bignoreOwner

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
      IF ((.NOT. bignoreOwner) .AND. &
          (IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .NE. 0)) THEN
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
        CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
          IF ((rdestMatrix%h_Kcol .NE. ST_NOHANDLE) .AND. &
              (rdestMatrix%h_Kld .NE. ST_NOHANDLE) .AND. &
              (rdestMatrix%h_Kdiagonal .NE. ST_NOHANDLE)) THEN
        
            CALL storage_getsize (rdestMatrix%h_Kcol,isize)
            bremove = bremove .OR. (isize .LT. rdestMatrix%NA)
            
            CALL storage_getsize (rsourceMatrix%h_Kld,isize)
            bremove = bremove .OR. (isize .LT. NEQ+1)
            
            CALL storage_getsize (rsourceMatrix%h_Kdiagonal,isize)
            bremove = bremove .OR. (isize .LT. NEQ)
          
          ELSE

            ! Remove any partial information if there is any.
            bremove = .TRUE.
            
          END IF
          
        CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
        
          IF ((rdestMatrix%h_Kcol .NE. ST_NOHANDLE) .AND. &
              (rdestMatrix%h_Kld .NE. ST_NOHANDLE)) THEN
 
            CALL storage_getsize (rdestMatrix%h_Kcol,isize)
            bremove = bremove .OR. (isize .LT. rdestMatrix%NA)
            
            CALL storage_getsize (rsourceMatrix%h_Kld,isize)
            bremove = bremove .OR. (isize .LT. NEQ+1)

          ELSE

            ! Remove any partial information if there is any.
            bremove = .TRUE.
            
          END IF
        
        CASE (LSYSSC_MATRIXD,LSYSSC_MATRIX1)
          ! Nothing to do

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
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
        CALL lsyssc_auxcopy_Kcol (rsourceMatrix,rdestMatrix)
        CALL lsyssc_auxcopy_Kld (rsourceMatrix,rdestMatrix)
        CALL lsyssc_auxcopy_Kdiagonal (rsourceMatrix,rdestMatrix)
        
      CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
        CALL lsyssc_auxcopy_Kcol (rsourceMatrix,rdestMatrix)
        CALL lsyssc_auxcopy_Kld (rsourceMatrix,rdestMatrix)
      
      CASE (LSYSSC_MATRIXD,LSYSSC_MATRIX1)
        ! Nothing to do

      END SELECT
      
      IF (.NOT. bignoreOwner) THEN
        ! Indicate via the matrixSpec-flag that we are the owner of the structure.
        rdestMatrix%imatrixSpec = IAND(rdestMatrix%imatrixSpec,&
                                      NOT(LSYSSC_MSPEC_STRUCTUREISCOPY))
      END IF
    
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
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7,&
          LSYSSC_MATRIX7INTL,LSYSSC_MATRIXD)
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
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7,&
          LSYSSC_MATRIX7INTL,LSYSSC_MATRIXD)
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
    
    SUBROUTINE copyContent(rsourceMatrix, rdestMatrix, bignoreOwner)
    
    ! The source matrix 
    TYPE(t_matrixScalar), INTENT(IN) :: rsourceMatrix

    ! The destination matrix 
    TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
    
    ! Whether to respect the ownership of the arrays and allocate 
    ! memory automatically if a matrix is not the owner of the arrays.
    LOGICAL, INTENT(IN) :: bignoreOwner
    
      ! local variables
      LOGICAL :: bremove
      REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2
      REAL(SP), DIMENSION(:), POINTER :: p_Fdata,p_Fdata2
      
      ! Overwrite structural data
      CALL copyStaticContent(rsourceMatrix, rdestMatrix)
      
      ! Check the content if it exists and if it has the right
      ! size - then we can overwrite!
      bRemove = .FALSE.
      
      ! But at first, check if rdestMatrix is the owner of the matrix
      ! structure:
      IF ((.NOT. bignoreOwner) .AND. &
          (IAND(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .NE. 0)) THEN
          
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
        CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7,LSYSSC_MATRIXD)
        
          IF (rdestMatrix%h_Da .NE. ST_NOHANDLE) THEN
        
            CALL storage_getsize (rdestMatrix%h_Da,isize)
            bremove = bremove .OR. (isize .LT. rdestMatrix%NA)
          
            ! Check the data type
            bremove = bremove .OR. (rdestMatrix%cdataType .NE. rsourceMatrix%cdataType)
          
          ELSE
          
            ! Remove any partial information if there is any
            bremove = .TRUE.
            
          END IF

          CASE (LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7INTL)
        
          IF (rdestMatrix%h_Da .NE. ST_NOHANDLE) THEN
        
            CALL storage_getsize (rdestMatrix%h_Da,isize)
            SELECT CASE(rdestMatrix%cinterleavematrixFormat)
            CASE (LSYSSC_MATRIX1)
              bremove = bremove .OR. &
                        (isize .LT. rdestMatrix%NA*rdestMatrix%NVAR*rdestMatrix%NVAR)
            CASE (LSYSSC_MATRIXD)
              bremove = bremove .OR. (isize .LT. rdestMatrix%NA*rdestMatrix%NVAR)
            CASE DEFAULT
              PRINT *, 'copyContent: wrong interleave matrix format'
              CALL sys_halt()
            END SELECT
          
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
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7,&
          LSYSSC_MATRIX7INTL,LSYSSC_MATRIXD,LSYSSC_MATRIX1)
        
        CALL lsyssc_auxcopy_da (rsourceMatrix,rdestMatrix)
      
      END SELECT
      
      IF (.NOT. bignoreOwner) THEN
        ! Indicate via the matrixSpec-flag that we are the owner of the structure.
        rdestMatrix%imatrixSpec = IAND(rdestMatrix%imatrixSpec,&
                                      NOT(LSYSSC_MSPEC_CONTENTISCOPY))
      END IF
    
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
  rvector%NVAR = 1
  rvector%cdataType = ST_DOUBLE
  rvector%iidxFirstEntry = 1
  rvector%bisCopy = .FALSE.
  rvector%isortStrategy = 0
  rvector%h_IsortPermutation = ST_NOHANDLE
  rvector%p_rspatialDiscretisation => NULL()
   
  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_releaseMatrixContent (rmatrix)
  
!<description>
  ! Releases the content of a matrix from memory. The structure stays unchanged.
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
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7&
          &,LSYSSC_MATRIX7INTL,LSYSSC_MATRIXD,LSYSSC_MATRIX1)
        ! Release matrix data, structure 9,7
        IF (rmatrix%h_DA .NE. ST_NOHANDLE) THEN
          CALL storage_free(rmatrix%h_DA)
        END IF
      END SELECT
    END IF
    
    rmatrix%h_DA        = ST_NOHANDLE
    
    ! Set the copy-flag appropriately
    rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,NOT(LSYSSC_MSPEC_CONTENTISCOPY))

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

  ! Matrix template; initialised by the default initialisation strategy
  ! of Fortran 90.
  TYPE(t_matrixScalar) :: rmatrixTemplate

  ! Release the matrix content.
  CALL lsyssc_releaseMatrixContent(rmatrix)
  
  ! Release the structure if it doesn't belong to another vector
  IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .EQ. 0) THEN
    ! Release matrix structure
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
      IF (rmatrix%h_Kcol .NE. ST_NOHANDLE)      CALL storage_free(rmatrix%h_Kcol)
      IF (rmatrix%h_Kld .NE. ST_NOHANDLE)       CALL storage_free(rmatrix%h_Kld)
      IF (rmatrix%h_Kdiagonal .NE. ST_NOHANDLE) CALL storage_free(rmatrix%h_Kdiagonal)
      
    CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
      IF (rmatrix%h_Kcol .NE. ST_NOHANDLE) CALL storage_free(rmatrix%h_Kcol)
      IF (rmatrix%h_Kld .NE. ST_NOHANDLE)  CALL storage_free(rmatrix%h_Kld)
      
    CASE (LSYSSC_MATRIXD)
      ! Nothing to do
    END SELECT
  END IF
  
  ! Clean up the rest. For this purpose, overwrite the matrix structure
  ! by a matrix structure which is initialised by the default Fortran
  ! initialisation routines.
  rmatrix = rmatrixTemplate
  !rmatrix%h_DA        = ST_NOHANDLE
  !rmatrix%h_Kcol      = ST_NOHANDLE
  !rmatrix%h_Kld       = ST_NOHANDLE
  !rmatrix%h_Kdiagonal = ST_NOHANDLE
  !rmatrix%cdataType   = ST_DOUBLE
  !rmatrix%cmatrixFormat = LSYSSC_MATRIXUNDEFINED
  !rmatrix%imatrixSpec = 0
  !rmatrix%NA  = 0
  !rmatrix%NEQ = 0
  !rmatrix%NCOLS = 0
  !rmatrix%NVAR = 1
  !rmatrix%dscaleFactor = 1.0_DP
  !rmatrix%isortStrategy = 0
  !rmatrix%h_IsortPermutation = ST_NOHANDLE
  !rmatrix%p_rspatialDiscretisation => NULL()

  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_clearMatrix (rmatrix,dvalue)
  
!<description>
  ! Clears the entries in a matrix. All entries are overwritten with 0.0 or
  ! with dvalue (if specified).
!</description>
  
!<inputoutput>
  
  ! Matrix to clear.
  TYPE(t_matrixScalar), INTENT(INOUT)               :: rmatrix
  
  ! OPTIONAL: Value to write into the matrix.
  ! If not specified, all matrix entries are set to 0.0.
  ! If specified, all matrix entries are set to dvalue.
  REAL(DP), INTENT(IN), OPTIONAL :: dvalue
  
!</inputoutput>

!</subroutine>

  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata

  IF (rmatrix%NEQ .LE. 0) RETURN ! Empty matrix

  IF (lsyssc_isExplicitMatrix1D(rmatrix)) THEN
    
    ! Get the data array and clear it -- depending on the data type.
    SELECT CASE (rmatrix%cdataType)
    CASE (ST_DOUBLE)
      CALL lsyssc_getbase_double (rmatrix,p_Ddata)
      IF (.NOT. PRESENT(dvalue)) THEN
        CALL lalg_clearVectorDble (p_Ddata)
      ELSE
        CALL lalg_setVectorDble (p_Ddata,dvalue)
      END IF
    CASE (ST_SINGLE)
      CALL lsyssc_getbase_single (rmatrix,p_Fdata)
      IF (.NOT. PRESENT(dvalue)) THEN
        CALL lalg_clearVectorSngl (p_Fdata)
      ELSE
        CALL lalg_setVectorSngl (p_Fdata,REAL(dvalue,SP))
      END IF
    CASE DEFAULT
      PRINT *,'lsyssc_clearMatrix: Unsupported Data type!'
      CALL sys_halt()
    END SELECT
  
  END IF

  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_initialiseIdentityMatrix (rmatrix)
  
!<description>
  ! Initialises the matrix rmatrix to an identity matrix.
  ! The matrix structure must already have been set up. If necessary, new data
  ! for the matrix is allocated on the heap. 
  ! The scaling factor of the matrix remains unchanged!
!</description>
  
!<inputoutput>
  
  ! Matrix to release.
  TYPE(t_matrixScalar), INTENT(INOUT)               :: rmatrix
  
!</inputoutput>

!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiagonal
    INTEGER(PREC_VECIDX) :: i

    IF (rmatrix%NEQ .LE. 0) RETURN ! Empty matrix

    ! Which matrix type do we have?
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7,&
        LSYSSC_MATRIXD,LSYSSC_MATRIX1)
      ! If necessary, allocate memory.
      IF (rmatrix%h_DA .EQ. ST_NOHANDLE) THEN
        CALL lsyssc_allocEmptyMatrix (rmatrix,LSYSSC_SETM_UNDEFINED)
      END IF
    END SELECT

    ! Now choose -- depending on the matrix format -- how to initialise the
    ! matrix data.
    SELECT CASE (rmatrix%cmatrixFormat)
    CASE (LSYSSC_MATRIX9)
      ! Clear the old content
      CALL lsyssc_clearMatrix (rmatrix)
      
      ! Get the structure and the data.
      ! Put the diagonal elements to 1.
      CALL lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
      
      SELECT CASE (rmatrix%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rmatrix,p_Ddata)
        DO i=1,rmatrix%NEQ
          p_Ddata(p_Kdiagonal(i)) = 1.0_DP
        END DO
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rmatrix,p_Fdata)
        DO i=1,rmatrix%NEQ
          p_Fdata(p_Kdiagonal(i)) = 1.0_SP
        END DO
      CASE DEFAULT
        PRINT *,'lsyssc_initialiseIdentityMatrix: Unsupported data type!'
        CALL sys_halt()
      END SELECT

    CASE (LSYSSC_MATRIX7)
      ! Clear the old content
      CALL lsyssc_clearMatrix (rmatrix)
      
      ! Get the structure and the data.
      ! Put the diagonal elements to 1.
      CALL lsyssc_getbase_Kld (rmatrix,p_Kdiagonal)
      
      SELECT CASE (rmatrix%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rmatrix,p_Ddata)
        DO i=1,rmatrix%NEQ
          p_Ddata(p_Kdiagonal(i)) = 1.0_DP
        END DO
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rmatrix,p_Fdata)
        DO i=1,rmatrix%NEQ
          p_Fdata(p_Kdiagonal(i)) = 1.0_SP
        END DO
      CASE DEFAULT
        PRINT *,'lsyssc_initialiseIdentityMatrix: Unsupported data type!'
        CALL sys_halt()
      END SELECT

    CASE (LSYSSC_MATRIX1)
      ! Clear the old content
      CALL lsyssc_clearMatrix (rmatrix)
      
      ! Get the structure and the data.
      ! Put the diagonal elements to 1.
      CALL lsyssc_getbase_Kld (rmatrix,p_Kdiagonal)
      
      SELECT CASE (rmatrix%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rmatrix,p_Ddata)
        DO i=1,rmatrix%NEQ
          p_Ddata(i*(rmatrix%NEQ-1)+i) = 1.0_DP
        END DO
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rmatrix,p_Fdata)
        DO i=1,rmatrix%NEQ
          p_Fdata(i*(rmatrix%NEQ-1)+i) = 1.0_SP
        END DO
      CASE DEFAULT
        PRINT *,'lsyssc_initialiseIdentityMatrix: Unsupported data type!'
        CALL sys_halt()
      END SELECT

    CASE (LSYSSC_MATRIXD)
      ! Put all elements to 1.0
      SELECT CASE (rmatrix%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rmatrix,p_Ddata)
        p_Ddata(:) = 1.0_DP
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rmatrix,p_Fdata)
        p_Fdata(:) = 1.0_SP
      CASE DEFAULT
        PRINT *,'lsyssc_initialiseIdentityMatrix: Unsupported data type!'
        CALL sys_halt()
      END SELECT
    
    CASE DEFAULT
      PRINT *,'lsyssc_initialiseIdentityMatrix: Unsupported matrix format!'
      CALL sys_halt()
    
    END SELECT
    
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_convertMatrix (rmatrix,cmatrixFormat,bconvertEntries)
  
!<description>
  ! Tries to convert a matrix rmatrix in-situ into a different matrix format.
  ! If necessary, memory is allocated for the new structure. (May be the case
  ! if the structure belongs to another matrix).
  ! If the matrix cannot be converted (due to format incompatibility),
  ! an error is thrown.
!</description>
  
!<input>
  ! Destination format of the matrix. One of the LSYSSC_MATRIXx constants.
  INTEGER, INTENT(IN)                 :: cmatrixFormat
  
  ! OPTIONAL: Whether to convert the entries of the matrix. Standard = TRUE.
  ! If set to FALSE, the entries are not converted, only the structure
  ! is converted (thus leaving the entries to an undefined state).
  LOGICAL, INTENT(IN), OPTIONAL       :: bconvertEntries
!</input>
  
!<inputoutput>
  ! Matrix to convert.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER(PREC_VECIDX) :: i
  INTEGER :: ihandle
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata,p_Ddata2
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld,p_Kdiagonal
  LOGICAL :: bentries,bstrucOwner,bentryOwner

  ! Matrix is already in that format.
  IF (rmatrix%cmatrixFormat .EQ. cmatrixFormat) RETURN
  
  ! Empty matrix
  IF (rmatrix%NEQ .LE. 0) RETURN 
  
  bentries = .TRUE.
  IF (PRESENT(bconvertEntries)) bentries = bconvertEntries
  
  ! Be careful with the structure during conversion: If the structure belongs
  ! to another matrix, we must not deallocate it!
  bstrucOwner = IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .EQ. 0
  
  ! The same for the entries!
  bentryOwner = IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .EQ. 0

  ! Which matrix type do we have?
  SELECT CASE (rmatrix%cmatrixFormat)
  CASE (LSYSSC_MATRIX9)
  
    SELECT CASE (cmatrixFormat)
    CASE (LSYSSC_MATRIX7)
    
      IF (.NOT. bstrucOwner) THEN
      
        ! Copy the structure so we are allowed to modify it.
        ihandle = ST_NOHANDLE
        CALL storage_copy(rmatrix%h_Kcol,ihandle)
        rmatrix%h_Kcol = ihandle

        ihandle = ST_NOHANDLE
        CALL storage_copy(rmatrix%h_Kld,ihandle)
        rmatrix%h_Kld = ihandle

        ! Kdiagonal is thrown away later, we don't need to copy it.
        
        rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,NOT(LSYSSC_MSPEC_STRUCTUREISCOPY))
        
      END IF
    
      IF ((.NOT. bentryOwner) .AND. bentries .AND. (rmatrix%h_DA .NE. ST_NOHANDLE)) THEN
      
        ! Copy the entries so we are allowed to modify them.
        ihandle = ST_NOHANDLE
        CALL storage_copy(rmatrix%h_Da,ihandle)
        rmatrix%h_Da = ihandle
        
        rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,NOT(LSYSSC_MSPEC_CONTENTISCOPY))
        
      END IF
    
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      CALL lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)

      ! Check that the matrix can be converted. There's a format error
      ! if there's no diagonal element.
      DO i=1,rmatrix%NEQ
        IF (p_Kcol(p_Kdiagonal(i)) .NE. i) THEN
          PRINT *,'lsyssc_convertMatrix: incompatible Format-9 matrix!'
          CALL sys_halt()
        END IF
      END DO
      
      IF ((.NOT. bentries) .OR. (rmatrix%h_DA .EQ. ST_NOHANDLE)) THEN
      
        ! No matrix entries, only resort the structure
        CALL lsyssc_unsortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)
        
        ! Release diagonal pointer if it belongs to us
        IF (bstrucOwner) THEN
          CALL storage_free (rmatrix%h_Kdiagonal)
        ELSE
          rmatrix%h_Kdiagonal = ST_NOHANDLE
        END IF

        rmatrix%cmatrixFormat = LSYSSC_MATRIX7
      
      ELSE
    
        ! Convert from structure 9 to structure 7. Use the sortCSRxxxx 
        ! routine below.
        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
        
          CALL lsyssc_getbase_double (rmatrix,p_Ddata)
          CALL lsyssc_unsortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ, p_Ddata)
          
          ! Release diagonal pointer if it belongs to us
          IF (bstrucOwner) THEN
            CALL storage_free (rmatrix%h_Kdiagonal)
          ELSE
            rmatrix%h_Kdiagonal = ST_NOHANDLE
          END IF

          rmatrix%cmatrixFormat = LSYSSC_MATRIX7
          
        CASE DEFAULT
          PRINT *,'lsyssc_convertMatrix: Unsupported data type!'
          CALL sys_halt()
        END SELECT

      END IF

    CASE (LSYSSC_MATRIXD)

      IF (bentries .AND. (rmatrix%h_DA .NE. ST_NOHANDLE)) THEN
      
        ! Convert from structure 9 to structure D by copying the
        ! diagonal to the front and reallocation of the memory
        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
        
          CALL lsyssc_getbase_double (rmatrix,p_Ddata)

          IF (.NOT. bentryOwner) THEN
            ! Allocate new memory for the entries
            rmatrix%h_Da = ST_NOHANDLE
            CALL storage_new ('lsyssc_convertMatrix', 'Da', &
                  rmatrix%NEQ, ST_DOUBLE, rmatrix%h_Da, ST_NEWBLOCK_NOINIT)
            CALL storage_getbase_double (rmatrix%h_Da,p_Ddata2)
            rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,NOT(LSYSSC_MSPEC_CONTENTISCOPY))
          ELSE
            ! Destinatinon pointer points to source matrix.
            ! This overwrites the original entries which is no problem as
            ! extracting the diagonal is always a 'compression' overwriting
            ! information that is not used anymore.
            p_Ddata2 => p_Ddata
          END IF
          
          CALL lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
          DO i=1,rmatrix%NEQ
            p_Ddata2(i) = p_Ddata(p_Kdiagonal(i))
          END DO

          rmatrix%cmatrixFormat = LSYSSC_MATRIXD
          
        CASE DEFAULT
          PRINT *,'lsyssc_convertMatrix: Unsupported data type!'
          CALL sys_halt()
        END SELECT

        ! Reallocate the entries-array to have only the diagonal entries.
        ! In case we are not the owner, we have allocated new memory and thus
        ! don't need to resize it again.
        IF (bentryOwner) THEN
          CALL storage_realloc ('lsyssc_convertMatrix', rmatrix%NEQ, &
                                rmatrix%h_Da, ST_NEWBLOCK_NOINIT, bentries)
        END IF
      
        rmatrix%NA = rmatrix%NEQ

      END IF
    
      ! Release unused information
      IF (bstrucOwner) THEN
        CALL storage_free (rmatrix%h_Kdiagonal)
        CALL storage_free (rmatrix%h_Kcol)
        CALL storage_free (rmatrix%h_Kld)
      ELSE
        rmatrix%h_Kdiagonal = ST_NOHANDLE
        rmatrix%h_Kcol = ST_NOHANDLE
        rmatrix%h_Kld = ST_NOHANDLE
      END IF
      
    CASE DEFAULT
      PRINT *,'lsyssc_convertMatrix: Cannot convert matrix!'
      CALL sys_halt()
    END SELECT

  CASE (LSYSSC_MATRIX7)
  
    SELECT CASE (cmatrixFormat)
    CASE (LSYSSC_MATRIX9)
    
      IF (.NOT. bstrucOwner) THEN
      
        ! Duplicate the structure in memory so we are allowed to modify it.
        ihandle = ST_NOHANDLE
        CALL storage_copy (rmatrix%h_Kcol,ihandle)
        rmatrix%h_Kcol = ihandle
        
        ihandle = ST_NOHANDLE
        CALL storage_copy (rmatrix%h_Kld,ihandle)
        rmatrix%h_Kld = ihandle
        
        rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,NOT(LSYSSC_MSPEC_STRUCTUREISCOPY))
        
      END IF

      IF ((.NOT. bentryOwner) .AND. bentries .AND. (rmatrix%h_DA .NE. ST_NOHANDLE)) THEN
      
        ! Copy the entries so we are allowed to modify them.
        ihandle = ST_NOHANDLE
        CALL storage_copy(rmatrix%h_Da,ihandle)
        rmatrix%h_Da = ihandle
        
        rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,NOT(LSYSSC_MSPEC_CONTENTISCOPY))
        
      END IF
    
      ! Convert from structure 7 to structure 9. Use the sortCSRxxxx 
      ! routine below.
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)

      ! Create a pointer to the diagonal
      CALL storage_new ('lsyssc_convertMatrix', 'Kdiagonal', &
            rmatrix%NEQ, ST_INT, rmatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT)

      ! WARNING: lsyssc_getbase_Kdiagonal does not(!) work here, because
      ! rmatrix is still in format LSYSSC_MATRIX7 and hence does not provide
      ! the array Kdiagonal ;-)
      CALL storage_getbase_int (rmatrix%h_Kdiagonal,p_Kdiagonal,rmatrix%NEQ)

      IF ((.NOT. bentries) .OR. (rmatrix%h_DA .EQ. ST_NOHANDLE)) THEN
      
        ! No matrix entries, only resort the structure
        CALL lsyssc_sortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)
        
        rmatrix%cmatrixFormat = LSYSSC_MATRIX9
      
      ELSE

        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
          CALL lsyssc_getbase_double (rmatrix,p_Ddata)
          CALL lsyssc_sortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ, p_Ddata)
          
          rmatrix%cmatrixFormat = LSYSSC_MATRIX9
          
        CASE DEFAULT
          PRINT *,'lsyssc_convertMatrix: Unsupported data type!'
          CALL sys_halt()
        END SELECT
      
      END IF
      
    CASE (LSYSSC_MATRIXD)

      IF (bentries .AND. (rmatrix%h_DA .NE. ST_NOHANDLE)) THEN
      
        ! Convert from structure 7 to structure D by copying the
        ! diagonal to the front and reallocation of the memory
        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
        
          CALL lsyssc_getbase_double (rmatrix,p_Ddata)

          IF (.NOT. bentryOwner) THEN
            ! Allocate new memory for the entries
            rmatrix%h_Da = ST_NOHANDLE
            CALL storage_new ('lsyssc_convertMatrix', 'Da', &
                  rmatrix%NEQ, ST_DOUBLE, rmatrix%h_Da, ST_NEWBLOCK_NOINIT)
            CALL storage_getbase_double (rmatrix%h_Da,p_Ddata2)
            rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,NOT(LSYSSC_MSPEC_CONTENTISCOPY))
          ELSE
            ! Destinatinon pointer points to source matrix.
            ! This overwrites the original entries which is no problem as
            ! extracting the diagonal is always a 'compression' overwriting
            ! information that is not used anymore.
            p_Ddata2 => p_Ddata
          END IF
          
          CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
          DO i=1,rmatrix%NEQ
            p_Ddata2(i) = p_Ddata(p_Kld(i))
          END DO

          rmatrix%cmatrixFormat = LSYSSC_MATRIXD
          
        CASE DEFAULT
          PRINT *,'lsyssc_convertMatrix: Unsupported data type!'
          CALL sys_halt()
        END SELECT

        ! Reallocate the entries-array to have only the diagonal entries.
        ! In case we are not the owner, we have allocated new memory and thus
        ! don't need to resize it again.
        IF (bentryOwner) THEN
          CALL storage_realloc ('lsyssc_convertMatrix', rmatrix%NEQ, &
                                rmatrix%h_Da, ST_NEWBLOCK_NOINIT, bentries)
        END IF
      
        rmatrix%NA = rmatrix%NEQ

      END IF
    
      ! Release unused information
      IF (bstrucOwner) THEN
        CALL storage_free (rmatrix%h_Kcol)
        CALL storage_free (rmatrix%h_Kld)
      ELSE
        rmatrix%h_Kcol = ST_NOHANDLE
        rmatrix%h_Kld = ST_NOHANDLE
      END IF
      
    CASE DEFAULT
      PRINT *,'lsyssc_convertMatrix: Cannot convert matrix!'
      CALL sys_halt()
    END SELECT

  CASE (LSYSSC_MATRIX9INTL)
    
    SELECT CASE (cmatrixFormat)
    CASE (LSYSSC_MATRIX7INTL)

      IF (.NOT. bstrucOwner) THEN
      
        ! Duplicate the structure in memory so we are allowed to modify it.
        ihandle = ST_NOHANDLE
        CALL storage_copy (rmatrix%h_Kcol,ihandle)
        rmatrix%h_Kcol = ihandle
        
        ihandle = ST_NOHANDLE
        CALL storage_copy (rmatrix%h_Kld,ihandle)
        rmatrix%h_Kld = ihandle
        
        ! Kdiagonal is thrown away later, we don't need to copy it.
        
        rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,NOT(LSYSSC_MSPEC_STRUCTUREISCOPY))
        
      END IF

      IF ((.NOT. bentryOwner) .AND. bentries .AND. (rmatrix%h_DA .NE. ST_NOHANDLE)) THEN
      
        ! Copy the entries so we are allowed to modify them.
        ihandle = ST_NOHANDLE
        CALL storage_copy(rmatrix%h_Da,ihandle)
        rmatrix%h_Da = ihandle
        
        rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,NOT(LSYSSC_MSPEC_CONTENTISCOPY))
        
      END IF
    
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      CALL lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)

      ! Check that the matrix can be converted. There's a format error
      ! if there's no diagonal element.
      DO i=1,rmatrix%NEQ
        IF (p_Kcol(p_Kdiagonal(i)) .NE. i) THEN
          PRINT *,'lsyssc_convertMatrix: incompatible Format-9 matrix!'
          CALL sys_halt()
        END IF
      END DO

      IF ((.NOT. bentries) .OR. (rmatrix%h_DA .EQ. ST_NOHANDLE)) THEN
        
        ! No matrix entries, only resort the structure
        CALL lsyssc_unsortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)
        
        ! Release diagonal pointer if it belongs to us
        IF (bstrucOwner) THEN
          CALL storage_free (rmatrix%h_Kdiagonal)
        ELSE
          rmatrix%h_Kdiagonal = ST_NOHANDLE
        END IF
        
        rmatrix%cmatrixFormat = LSYSSC_MATRIX7INTL
        
      ELSE

        ! Convert from structure 9 to structure 7. Use the sortCSRxxxx 
        ! routine below.
        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
        
          CALL lsyssc_getbase_double (rmatrix,p_Ddata)
          SELECT CASE(rmatrix%cinterleavematrixFormat)
          CASE (LSYSSC_MATRIX1)
            CALL lsyssc_unsortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal,&
                & rmatrix%NEQ, p_Ddata, rmatrix%NVAR*rmatrix%NVAR)
                
          CASE (LSYSSC_MATRIXD)
            CALL lsyssc_unsortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal,&
                & rmatrix%NEQ, p_Ddata, rmatrix%NVAR)
                
          CASE DEFAULT
            PRINT *, 'lsyssc_convertMatrix: Unsupported interleave matrix!'
            CALL sys_halt()
          END SELECT

          ! Release diagonal pointer if it belongs to us
          IF (bstrucOwner) THEN
            CALL storage_free (rmatrix%h_Kdiagonal)
          ELSE
            rmatrix%h_Kdiagonal = ST_NOHANDLE
          END IF

          rmatrix%cmatrixFormat = LSYSSC_MATRIX7INTL
          
        CASE DEFAULT
          PRINT *,'lsyssc_convertMatrix: Unsupported data type!'
          CALL sys_halt()
        END SELECT

      END IF

    CASE DEFAULT
      PRINT *,'lsyssc_convertMatrix: Cannot convert matrix!'
      CALL sys_halt()
    END SELECT
    
  CASE (LSYSSC_MATRIX7INTL)

    SELECT CASE (cmatrixFormat)
    CASE (LSYSSC_MATRIX9INTL)
      
      IF (.NOT. bstrucOwner) THEN
      
        ! Duplicate the structure in memory so we are allowed to modify it.
        ihandle = ST_NOHANDLE
        CALL storage_copy (rmatrix%h_Kcol,ihandle)
        rmatrix%h_Kcol = ihandle
        
        ihandle = ST_NOHANDLE
        CALL storage_copy (rmatrix%h_Kld,ihandle)
        rmatrix%h_Kld = ihandle
        
        rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,NOT(LSYSSC_MSPEC_STRUCTUREISCOPY))
        
      END IF

      IF ((.NOT. bentryOwner) .AND. bentries .AND. (rmatrix%h_DA .NE. ST_NOHANDLE)) THEN
      
        ! Copy the entries so we are allowed to modify them.
        ihandle = ST_NOHANDLE
        CALL storage_copy(rmatrix%h_Da,ihandle)
        rmatrix%h_Da = ihandle
        
        rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,NOT(LSYSSC_MSPEC_CONTENTISCOPY))
        
      END IF
    
      ! Convert from structure 7 to structure 9. Use the sortCSRxxxx 
      ! routine below.
      CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
      
      ! Create a pointer to the diagonal
      CALL storage_new ('lsyssc_convertMatrix', 'Kdiagonal', &
          rmatrix%NEQ, ST_INT, rmatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT)
      CALL lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)

      IF ((.NOT. bentries) .OR. (rmatrix%h_DA .EQ. ST_NOHANDLE)) THEN
        
        ! No matrix entries, only resort the structure
        CALL lsyssc_sortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)
        
        rmatrix%cmatrixFormat = LSYSSC_MATRIX9INTL
      
      ELSE

        SELECT CASE (rmatrix%cdataType)
        CASE (ST_DOUBLE)
          CALL lsyssc_getbase_double (rmatrix,p_Ddata)
          SELECT CASE(rmatrix%cinterleavematrixFormat)
          CASE (LSYSSC_MATRIX1)
            CALL lsyssc_sortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal,&
                & rmatrix%NEQ, p_Ddata, rmatrix%NVAR*rmatrix%NVAR)
                
          CASE (LSYSSC_MATRIXD)
            CALL lsyssc_sortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal,&
                & rmatrix%NEQ, p_Ddata, rmatrix%NVAR)
                
          CASE DEFAULT
            PRINT *, 'lsyssc_convertMatrix: Unsupported interleave matrix format!'
            CALL sys_halt()
          END SELECT
          rmatrix%cmatrixFormat = LSYSSC_MATRIX9INTL
          
        CASE DEFAULT
          PRINT *,'lsyssc_convertMatrix: Unsupported data type!'
          CALL sys_halt()
        END SELECT
        
      END IF

    CASE DEFAULT
      PRINT *,'lsyssc_convertMatrix: Cannot convert matrix!'
      CALL sys_halt()
    END SELECT

  CASE DEFAULT
    PRINT *,'lsyssc_convertMatrix: Cannot convert matrix!'
    CALL sys_halt()
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
  INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld

  ! Column structure of the matrix
  INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Kcol

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
    ! Search for the first element on or above the diagonal.
    DO j = Kld(i), Kld(i+1)-1

      ! Check if we reached the position of the diagonal entry...
      IF (Kcol(j) .GE. i) EXIT

    END DO

    ! Save the position of the diagonal entry
    Kdiagonal(i) = j

  END DO
          
  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_sortCSRdouble (Kcol, Kld, Kdiagonal, neq, Da, nintl)

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

  ! Dimension of the interleaved submatrices (if any)
  INTEGER, INTENT(IN), OPTIONAL :: nintl
!</input>

!<inputoutput>
  ! OPTIONAL:
  ! On input:  the matrix entries to be resorted,
  ! On output: the resorted matrix entries
  REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Da
  
  ! On input:  the column numbers to be resorted,
  ! On output: the resorted column numbers
  INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Kcol
!</inputoutput>

!<output>
  ! Pointers to the diagonal entries of the matrix.
  INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(OUT) :: Kdiagonal
!</output>

!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: Daux
  REAL(DP) :: aux
  INTEGER :: h_Daux
  INTEGER(I32) :: i,j

  IF (PRESENT(Da)) THEN

    ! Check if size of interleave matrix is specified. If this is the
    ! case then each "move" is performed for a local vector
    IF (PRESENT(nintl)) THEN

      ! Allocate memory for auxiliary vector
      CALL storage_new1D ('lsyssc_sortCSRdouble', 'Daux', INT(nintl,I32), &
          ST_DOUBLE, h_Daux,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double(h_Daux,Daux)
      
      ! loop through each row
      DO i = 1, neq
        
        ! Take the diagonal element
        Daux = Da(nintl*(Kld(i)-1)+1:nintl*Kld(i))
        
        ! Loop through each column in this row.
        ! Shift every entry until the diagonal is reached.
        DO j = Kld(i)+1, Kld(i+1)-1
          
          ! Check if we reached the position of the diagonal entry...
          IF (Kcol(j)>i) EXIT
          
          Kcol(j-1) = Kcol(j)
          Da(nintl*(j-2)+1:nintl*(j-1)) = Da(nintl*(j-1)+1:nintl*j)
          
        END DO
        
        ! If we have reached the diagonal, we can stop and save our
        ! diagonal entry from the first position there. The rest of the
        ! line is in ascending order according to the specifications of
        ! storage technique 7.
        
        Kcol(j-1) = i
        Da(nintl*(j-2)+1:nintl*(j-1)) = Daux
        
        ! Save the position of the diagonal entry
        Kdiagonal(i) = j-1
        
      END DO

      ! Free auxiliary memory
      CALL storage_free(h_Daux)
      
    ELSE

      ! loop through each row
      DO i = 1, neq
        
        ! Take the diagonal element
        aux = Da(Kld(i))
        
        ! Loop through each column in this row.
        ! Shift every entry until the diagonal is reached.
        DO j = Kld(i)+1, Kld(i+1)-1
          
          ! Check if we reached the position of the diagonal entry...
          IF (Kcol(j)>i) EXIT
          
          Kcol(j-1) = Kcol(j)
          Da(j-1) = Da(j)
          
        END DO
        
        ! If we have reached the diagonal, we can stop and save our
        ! diagonal entry from the first position there. The rest of the
        ! line is in ascending order according to the specifications of
        ! storage technique 7.
        
        Kcol(j-1) = i
        Da(j-1) = aux
        
        ! Save the position of the diagonal entry
        Kdiagonal(i) = j-1
        
      END DO

    END IF

  ELSE

    ! loop through each row
    DO i = 1, neq

      ! Loop through each column in this row.
      ! Shift every entry until the diagonal is reached.
      DO j = Kld(i)+1, Kld(i+1)-1

        ! Check if we reached the position of the diagonal entry...
        IF (Kcol(j)>i) EXIT

        Kcol(j-1) = Kcol(j)

      END DO

      ! If we have reached the diagonal, we can stop and save our
      ! diagonal entry from the first position there. The rest of the
      ! line is in ascending order according to the specifications of
      ! storage technique 7.

      Kcol(j-1) = i
      
      ! Save the position of the diagonal entry
      Kdiagonal(i) = j-1

    END DO
          
  END IF
  
  END SUBROUTINE
  
  !****************************************************************************

!<subroutine>  
  
  SUBROUTINE lsyssc_unsortCSRdouble (Kcol, Kld, Kdiagonal, neq, Da, nintl)

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

  ! Dimension of the interleave submatrices (if any)
  INTEGER, INTENT(IN), OPTIONAL :: nintl
!</input>

!<inputoutput>
  ! OPTIONAL:
  ! On input:  the matrix entries to be resorted,
  ! On output: the resorted matrix entries
  REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Da
  
  ! On input:  the column numbers to be resorted,
  ! On output: the resorted column numbers
  INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Kcol
!</inputoutput>

!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: Daux
  REAL(DP) :: aux
  INTEGER :: h_Daux
  INTEGER(I32) :: i, j

  IF (PRESENT(Da)) THEN

    ! Check if size of interleave matrix is specified. If this is the
    ! case then each "move" is performed for a local vector
    IF (PRESENT(nintl)) THEN

      ! Allocate memory for auxiliary vector
      CALL storage_new1D ('lsyssc_sortCSRdouble', 'Daux', INT(nintl,I32), &
          ST_DOUBLE, h_Daux,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double(h_Daux,Daux)

      ! loop through each row
      DO i = 1, neq
        
        ! Take the diagonal element
        Daux = Da(nintl*(Kdiagonal(i)-1)+1:nintl*Kdiagonal(i))
        
        ! Loop through each column in this row.
        ! Shift every entry one element to the right.
        DO j = Kdiagonal(i),Kld(i)+1,-1
          
          Kcol(j) = Kcol(j-1)
          Da(nintl*(j-1)+1:nintl*j) = Da(nintl*(j-2)+1:nintl*(j-1))
          
        END DO
        
        ! Put the diagonal to the front.
        Kcol(Kld(i)) = i
        Da(nintl*(Kld(i)-1)+1:nintl*Kld(i)) = Daux

      END DO

    ELSE
      
      ! loop through each row
      DO i = 1, neq
        
        ! Take the diagonal element
        aux = Da(Kdiagonal(i))
        
        ! Loop through each column in this row.
        ! Shift every entry one element to the right.
        DO j = Kdiagonal(i),Kld(i)+1,-1
          
          Kcol(j) = Kcol(j-1)
          Da(j) = Da(j-1)
          
        END DO
        
        ! Put the diagonal to the front.
        Kcol(Kld(i)) = i
        Da(Kld(i)) = aux
        
      END DO
      
    END IF
      
  ELSE

    ! loop through each row
    DO i = 1, neq

      ! Loop through each column in this row.
      ! Shift every entry one element to the right.
      DO j = Kdiagonal(i),Kld(i)+1,-1
        Kcol(j) = Kcol(j-1)
      END DO
      
      ! Put the diagonal to the front.
      Kcol(Kld(i)) = i
      
    END DO

  END IF
          
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_synchroniseSortVecVec (rvectorSrc,rvectorDst,rtemp)
  
!<description>
  ! Synchronises the sorting strategy of rvectorDest according to rvectorSrc:
  ! If rvectorSrc is unsorted, rvectorDest will be unsorted (without changing
  ! the attached sorting permutation).
  ! If rvectorSrc is sorted differently to rvectorDest, rvectorDest is
  ! sorted according to rvectorSrc.
  ! Therefore if the routine is finished, rvectorDest 'looks like'
  ! rvectorSrc according to the sorting strategy.
!</description>

!<input>
  ! Source vector defining the sorting strategy.
  TYPE(t_vectorScalar), INTENT(IN)               :: rvectorSrc
!</input>

!<inputoutput>
  ! Destination vector; is resorted according to the sort strategy in rvectorSrc
  ! or is unsorted, if rvectorSrc is unsorted.
  ! Must have the same size as rvectorSrc.
  TYPE(t_vectorScalar), INTENT(INOUT)            :: rvectorDst

  ! A temporary vector. Must be of the same data type as rvector.
  ! Must be at least as large as rvectorDst.
  TYPE(t_vectorScalar), INTENT(INOUT)            :: rtemp
!</inputoutput>

!</subroutine>

    IF (rvectorSrc%NEQ .NE. rvectorDst%NEQ) THEN
      PRINT *,'lsyssc_synchroniseSortVecVec: Vectors have different size!'
      CALL sys_halt()
    END IF

    IF (rtemp%NEQ .LT. rvectorDst%NEQ) THEN
      PRINT *,'lsyssc_synchroniseSortVecVec: Auxiliary vector too small!'
      CALL sys_halt()
    END IF

    ! If both are unsorted or both sorted in the same way, there's nothing to do.
    IF ((rvectorSrc%isortStrategy .EQ. rvectorDst%isortStrategy) .OR. &
        ((rvectorSrc%isortStrategy .LT. 0) .AND. (rvectorDst%isortStrategy .LT. 0))) &
      RETURN

    ! Should rvectorDst be unsorted?
    IF (rvectorSrc%isortStrategy .LT. 0) THEN
      CALL lsyssc_vectorActivateSorting (rvectorDst,.FALSE.,rtemp)
      RETURN
    END IF

    ! rvectorDst is differently sorted than rvectorDst; synchronise them!
    CALL lsyssc_sortVectorInSitu (rvectorDst,rtemp,&
         rvectorSrc%isortStrategy,rvectorSrc%h_IsortPermutation)

  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
  
  SUBROUTINE lsyssc_synchroniseSortMatVec (rmatrixSrc,rvectorDst,rtemp)
  
!<description>
  ! Synchronises the sorting strategy of rvectorDest according to rmatrixSrc:
  ! If rmatrixSrc is unsorted, rvectorDest will be unsorted (without changing
  ! the attached sorting permutation).
  ! If rmatrixSrc is sorted differently to rvectorDest, rvectorDest is
  ! sorted according to rmatrixSrc.
  ! Therefore if the routine is finished, rvectorDest is compatible to
  ! rmatrixSrc according to the sorting strategy.
!</description>

!<input>
  ! Source matrix defining the sorting strategy.
  TYPE(t_matrixScalar), INTENT(IN)               :: rmatrixSrc
!</input>

!<inputoutput>
  ! Destination vector; is resorted according to the sort strategy in rmatrixSrc
  ! or is unsorted, if rmatrixSrc is unsorted.
  ! Must have the same size (NEQ) as rmatrixSrc.
  TYPE(t_vectorScalar), INTENT(INOUT)            :: rvectorDst

  ! A temporary vector. Must be of the same data type as rvector.
  ! Must be at least as large as rvectorDst.
  TYPE(t_vectorScalar), INTENT(INOUT)            :: rtemp
!</inputoutput>

!</subroutine>

    IF (rmatrixSrc%NEQ .NE. rvectorDst%NEQ) THEN
      PRINT *,'lsyssc_synchroniseSortMatVec: Matrix and vector have different size!'
      CALL sys_halt()
    END IF

    IF (rtemp%NEQ .LT. rvectorDst%NEQ) THEN
      PRINT *,'lsyssc_synchroniseSortMatVec: Auxiliary vector too small!'
      CALL sys_halt()
    END IF


    ! If both are unsorted or both sorted in the same way, there's nothing to do.
    IF ((rmatrixSrc%isortStrategy .EQ. rvectorDst%isortStrategy) .OR. &
        ((rmatrixSrc%isortStrategy .LT. 0) .AND. (rvectorDst%isortStrategy .LT. 0))) &
      RETURN

    ! Should rvectorDst be unsorted?
    IF (rmatrixSrc%isortStrategy .LT. 0) THEN
      CALL lsyssc_vectorActivateSorting (rvectorDst,.FALSE.,rtemp)
      RETURN
    END IF

    ! rvectorDst is differently sorted than rvectorDst; synchronise them!
    CALL lsyssc_sortVectorInSitu (rvectorDst,rtemp,&
         rmatrixSrc%isortStrategy,rmatrixSrc%h_IsortPermutation)

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
      CALL sys_halt()
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
        CALL sys_halt()
        
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
        (rvector%h_IsortPermutation .NE. ST_NOHANDLE) .AND. &
        (rvector%isortStrategy .GT. 0)) THEN
        
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
        CALL sys_halt()
        
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
      CALL sys_halt()
      
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
    CALL lsyssc_copyVector (rvector,rtempLocal)
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
  !
  ! WARNING: This routine does NOT change any information (structure or content)
  !  which is marked as 'shared' with another matrix! Therefore, if the structure
  !  (column/row structure) of a matrix A belongs to another matrix B and 
  !  bsortEntries=TRUE, only the entries of matrix A are sorted, not the structure!
  !  So for sorting multiple matrices which all share the same structure, first
  !  all 'child' matrices have to be sorted before the 'parent' matrix is sorted!
!</description>
  
!<inputoutput>
  ! Vector to resort
  TYPE(t_matrixScalar), INTENT(INOUT), TARGET       :: rmatrix
!</inputoutput>

!<input>
  ! Sort the entries or only the structure of the matrix.
  ! = FALSE: Only sort the structure of the matrix,
  ! = TRUE : Sort both, entries and structure of the matrix.
  LOGICAL, INTENT(IN) :: bsortEntries

  ! OPTIONAL: Identifier for the sorting strategy to apply to the matrix.
  ! This is usually one of the SSTRAT_xxxx constants from the module
  ! 'sortstrategy', although it's actually used here as follows:
  ! <=0: Calculate the unsorted matrix
  !  >0: Resort the vector according to a permutation;
  !      this is either the permutation specified in the vector
  !      or that one identified by h_IsortPermutation
  ! If not specified, the sorting of the matrix is activated using
  ! the isortStrategy specifier in the matrix -- i.e. this 'activates'
  ! a previously attached sorting.
 INTEGER, INTENT(IN), OPTIONAL                      :: isortStrategy
  
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
  INTEGER :: h_Iperm, isortStrat
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Iperm
  INTEGER(PREC_VECIDX) :: NEQ
  TYPE(t_matrixScalar), POINTER :: p_rmatrix
  LOGICAL :: bsortEntriesTmp
  
  isortStrat = ABS(rmatrix%isortStrategy)
  IF (PRESENT(isortStrategy)) isortStrat = isortStrategy
  
    ! Desired sorting strategy and currently active sorting strategy identical?
    IF (.NOT. PRESENT(h_IsortPermutation)) THEN
    
      IF (isortStrat .EQ. rmatrix%isortStrategy) RETURN
      IF ((isortStrat .LE. 0) .AND. (rmatrix%isortStrategy .LE. 0)) RETURN
    
    ELSE
    
      IF ((isortStrat .LE. 0) .AND. (rmatrix%isortStrategy .LE. 0)) THEN
        IF (h_IsortPermutation .NE. rmatrix%h_IsortPermutation) THEN
          ! Matrix is unsorted and should stay unsorted, but
          ! permutation should change.
          rmatrix%isortStrategy = isortStrat
          rmatrix%h_IsortPermutation = h_IsortPermutation
        END IF
        RETURN
      END IF
      IF ((isortStrat .GT. 0) .AND. (rmatrix%isortStrategy .GT. 0) .AND.&
          (h_IsortPermutation .EQ. rmatrix%h_IsortPermutation)) RETURN
    
    END IF
    
    NEQ = rmatrix%NEQ
    
    ! Sort the matrix back?
    IF (isortStrat .LE. 0) THEN
      ! Get the permutation that describes how to resort the matrix:
      CALL storage_getbase_int(rmatrix%h_IsortPermutation,p_Iperm)
      
      ! Exchange the roles of the first and second part of p_Iperm.
      ! This makes the permutation to the inverse permutation and vice versa.
      ! Call the resort-subroutine with that to do the actual sorting.
      CALL do_matsort (rmatrix,p_Iperm(NEQ+1:NEQ*2),p_Iperm(1:NEQ),bsortEntries)
      
      ! Inform the vector about which sorting strategy we now use.
      rmatrix%isortStrategy = isortStrat
      RETURN
    END IF
    
    p_rmatrix => rmatrix
    bsortEntriesTmp = bsortEntries
    
    ! Get the actual sorting strategy.
    h_Iperm = rmatrix%h_IsortPermutation
    IF (PRESENT(h_IsortPermutation)) h_Iperm = h_IsortPermutation 
    
    ! Do we have to sort back before resorting?
    IF ((h_Iperm .NE. rmatrix%h_IsortPermutation) .AND. &
        (rmatrix%h_IsortPermutation .NE. ST_NOHANDLE) .AND. &
        (rmatrix%isortStrategy .GT. 0)) THEN

      ! That's a little bit tricky now if structure and/or content of rmatrix
      ! is shared with another matrix!
      ! If that's the case, we make a copy of our matrix and work with that.
      IF (.NOT. bsortEntries) THEN
        IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .NE. 0) THEN
          ! The user wants to sort the structure, but the structure belongs to 
          ! another matrix! Sorting would destroy that matrix, so we don't allow
          ! that. Therefore, there's nothing to do here.
          RETURN
        END IF
      ELSE
        IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY) .EQ. LSYSSC_MSPEC_ISCOPY) THEN
          ! Structure and content belongs to another matrix! Sorting would 
          ! destroy that matrix, so we don't allow that. Therefore, there'
          ! s nothing to do here.
          RETURN
        END IF
        IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .EQ. &
            LSYSSC_MSPEC_CONTENTISCOPY) THEN
          ! The user wants to sort structure and content, but the content belongs to 
          ! another matrix! Sorting would destroy that matrix, so we don't allow
          ! that. Therefore, we only sort the structure.
          bsortEntriesTmp = .FALSE.  
        END IF
        IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .EQ. &
            LSYSSC_MSPEC_STRUCTUREISCOPY) THEN
          ! The user wants to sort structure and content, but the structure belongs to 
          ! another matrix! Sorting would destroy that matrix, so we don't allow
          ! that. This situation is a little but tricky, as we want to sort
          ! the entries without sorting the structure AND we must sort back
          ! at first. We have no chance, we allocate another matrix as a copy
          ! of rmatrix and work with that. The first unsorting below will therefore
          ! unsort structure AND content, so we can sort it later.
          ALLOCATE (p_rmatrix)
          CALL lsyssc_duplicateMatrix (rmatrix,p_rmatrix,LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
        END IF
      END IF

      ! Sort back at first - with the associated permutation
      CALL storage_getbase_int(p_rmatrix%h_IsortPermutation,p_Iperm)
      
      ! Exchange the roles of the first and second part of p_Iperm.
      ! This makes the permutation to the inverse permutation and vice versa.
      ! Call the resort-subroutine with that to do the actual sorting.
      CALL do_matsort (p_rmatrix,p_Iperm(NEQ+1:NEQ*2),p_Iperm(1:NEQ),bsortEntries)
      
      ! If p_rmatrix is our local copy, manually change the ownership of the
      ! structure. This will prevent the sorting routine below from sorting
      ! the structure, which saves some time. But we have to remember to
      ! switch the ownership back before releasing our local matrix!
      IF (.NOT. ASSOCIATED(p_rmatrix,rmatrix)) THEN
        p_rmatrix%imatrixSpec = IOR(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY)
      END IF
      
    END IF
    
    ! Now sort the vector according to h_Iperm
    CALL storage_getbase_int(h_Iperm,p_Iperm)

    ! This time, we don't exchange the roles of the permutation and
    ! its inverse :-)
    CALL do_matsort (p_rmatrix,p_Iperm(1:NEQ),p_Iperm(NEQ+1:NEQ*2),bsortEntries)

    ! If p_rmatrix is our local copy, copy the content to rmatrix and release the
    ! local copy. Because of the small 'hack' above, we first restore the 
    ! ownership status and then release the matrix.    
    IF (.NOT. ASSOCIATED(p_rmatrix,rmatrix)) THEN
      CALL lsyssc_duplicateMatrix (p_rmatrix,rmatrix,LSYSSC_DUP_COPY,LSYSSC_DUP_IGNORE)
      p_rmatrix%imatrixSpec = IAND(rmatrix%imatrixSpec,NOT(LSYSSC_MSPEC_STRUCTUREISCOPY))
      CALL lsyssc_releaseMatrix (p_rmatrix)
    END IF
    
    ! Inform the vector about which sorting strategy we now use.
    rmatrix%isortStrategy = isortStrat
    
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
      
        IF (.NOT. bsortEntries) THEN
        
          ! Sort only the structure of the matrix, keep the entries
          ! unchanged. Don't do anything if the structure does not belong
          ! to the matrix!
          IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .EQ. 0) THEN
          
            ! Duplicate the matrix before resorting.
            CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                        LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)
                                        
            ! Get the structure of the original and the temporary matrix
            CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
            CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
            CALL lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiag)
            CALL lsyssc_getbase_Kcol (rtempMatrix,p_KcolTmp)
            CALL lsyssc_getbase_Kld (rtempMatrix,p_KldTmp)
            
            ! Sort
            CALL lsyssc_sortMat9Struc (p_Kcol, p_KcolTmp, p_Kld, p_KldTmp, &
                                       p_Kdiag, Itr1, Itr2, NEQ)        
                                       
            ! Remove temp matrix
            CALL lsyssc_releaseMatrix(rtempMatrix)
                                       
          END IF
        
        ELSE
        
          ! Don't do anything if structure & entries don't belong to our matrix.
          IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY) .EQ. LSYSSC_MSPEC_ISCOPY) &
            RETURN
        
          ! If the structure does not belong to the matrix, only sort the!
          ! entries. Otherwise, sort structure + entries
          
          IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .NE. 0) THEN

            ! Get the structure of the matrix
            CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
            CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
            CALL lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiag)

            ! Create a copy of the matrix entries
            CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                        LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

            SELECT CASE (rmatrix%cdataType)
            CASE (ST_DOUBLE)
              ! Double precision version
              CALL lsyssc_getbase_double (rmatrix,p_Ddata)
              CALL lsyssc_getbase_double (rtempMatrix,p_DdataTmp)
              CALL lsyssc_sortMat9Ent_double (p_Ddata,p_DdataTmp,p_Kcol, &
                                              p_Kld, Itr1, Itr2, NEQ)        
            CASE (ST_SINGLE)
              ! Single precision version
              CALL lsyssc_getbase_single (rmatrix,p_Fdata)
              CALL lsyssc_getbase_single (rtempMatrix,p_FdataTmp)
              CALL lsyssc_sortMat9Ent_single (p_Fdata,p_FdataTmp,p_Kcol, &
                                              p_Kld, Itr1, Itr2, NEQ)        
            CASE DEFAULT
              PRINT *,'lsyssc_sortMatrix: Unsupported data type.'
              CALL sys_halt()
            END SELECT
          
            ! Remove temp matrix
            CALL lsyssc_releaseMatrix(rtempMatrix)
          
          ELSE
          
            ! Duplicate the matrix before resorting.
            ! We need either a copy only of the structure or of the full matrix.
            CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                         LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
            
            ! Get the structure of the original and the temporary matrix
            CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
            CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
            CALL lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiag)
            CALL lsyssc_getbase_Kcol (rtempMatrix,p_KcolTmp)
            CALL lsyssc_getbase_Kld (rtempMatrix,p_KldTmp)
          
            SELECT CASE (rmatrix%cdataType)
            CASE (ST_DOUBLE)
              ! Double precision version
              CALL lsyssc_getbase_double (rmatrix,p_Ddata)
              CALL lsyssc_getbase_double (rtempMatrix,p_DdataTmp)
              CALL lsyssc_sortMat9_double (p_Ddata,p_DdataTmp,p_Kcol, p_KcolTmp, &
                                           p_Kld, p_KldTmp, p_Kdiag, &
                                           Itr1, Itr2, NEQ)        
            CASE (ST_SINGLE)
              ! Single precision version
              CALL lsyssc_getbase_single (rmatrix,p_Fdata)
              CALL lsyssc_getbase_single (rtempMatrix,p_FdataTmp)
              CALL lsyssc_sortMat9_single (p_Fdata,p_FdataTmp,p_Kcol, p_KcolTmp, &
                                           p_Kld, p_KldTmp, p_Kdiag, &
                                           Itr1, Itr2, NEQ)        
            CASE DEFAULT
              PRINT *,'lsyssc_sortMatrix: Unsupported data type.'
              CALL sys_halt()
            END SELECT
            
            ! Remove temp matrix
            CALL lsyssc_releaseMatrix(rtempMatrix)
         
          END IF   
         
        END IF

      CASE (LSYSSC_MATRIX7)

        IF (.NOT. bsortEntries) THEN
        
          ! Duplicate the matrix before resorting.
          CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                       LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)
        
          ! Get the structure of the original and the temporary matrix
          CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
          CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
          CALL lsyssc_getbase_Kcol (rtempMatrix,p_KcolTmp)
          CALL lsyssc_getbase_Kld (rtempMatrix,p_KldTmp)
        
          ! Sort only the structure of the matrix, keep the entries
          ! unchanged.
          CALL lsyssc_sortMat7Struc (p_Kcol, p_KcolTmp, p_KldTmp, p_KldTmp, &
                                     Itr1, Itr2, NEQ)        
        
          ! Remove temp matrix
          CALL lsyssc_releaseMatrix(rtempMatrix)
        
        ELSE
        
          ! Don't do anything if structure & entries don't belong to our matrix.
          IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY) .EQ. LSYSSC_MSPEC_ISCOPY) &
            RETURN
        
          ! If the structure does not belong to the matrix, only sort the!
          ! entries. Otherwise, sort structure + entries
          
          IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .NE. 0) THEN

            ! Get the structure of the matrix
            CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
            CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
        
            ! Create a copy of the matrix entries
            CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                        LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        
            ! Sort only the entries
            SELECT CASE (rmatrix%cdataType)
            CASE (ST_DOUBLE)
              ! Double precision version
              CALL lsyssc_getbase_double (rmatrix,p_Ddata)
              CALL lsyssc_getbase_double (rtempMatrix,p_DdataTmp)
              CALL lsyssc_sortMat7Ent_double (p_Ddata,p_DdataTmp,p_Kcol, &
                                              p_Kld, Itr1, Itr2, NEQ)        
            CASE (ST_SINGLE)
              ! Single precision version
              CALL lsyssc_getbase_single (rmatrix,p_Fdata)
              CALL lsyssc_getbase_single (rtempMatrix,p_FdataTmp)
              CALL lsyssc_sortMat7Ent_single (p_Fdata,p_FdataTmp,p_Kcol, &
                                              p_Kld, Itr1, Itr2, NEQ)        
            CASE DEFAULT
              PRINT *,'lsyssc_sortMatrix: Unsupported data type.'
              CALL sys_halt()
            END SELECT

            ! Remove temp matrix
            CALL lsyssc_releaseMatrix(rtempMatrix)

          ELSE        
        
            ! Duplicate the matrix before resorting.
            CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                         LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
        
            ! Get the structure of the original and the temporary matrix
            CALL lsyssc_getbase_Kcol (rmatrix,p_Kcol)
            CALL lsyssc_getbase_Kld (rmatrix,p_Kld)
            CALL lsyssc_getbase_Kcol (rtempMatrix,p_KcolTmp)
            CALL lsyssc_getbase_Kld (rtempMatrix,p_KldTmp)
        
            ! Sort structure + entries
            SELECT CASE (rmatrix%cdataType)
            CASE (ST_DOUBLE)
              ! Double precision version
              CALL lsyssc_getbase_double (rmatrix,p_Ddata)
              CALL lsyssc_getbase_double (rtempMatrix,p_DdataTmp)
              CALL lsyssc_sortMat7_double (p_Ddata,p_DdataTmp,p_Kcol, p_KcolTmp, &
                                          p_Kld, p_KldTmp, &
                                          Itr1, Itr2, NEQ)        
            CASE (ST_SINGLE)
              ! Single precision version
              CALL lsyssc_getbase_single (rmatrix,p_Fdata)
              CALL lsyssc_getbase_single (rtempMatrix,p_FdataTmp)
              CALL lsyssc_sortMat7_single (p_Fdata,p_FdataTmp,p_Kcol, p_KcolTmp, &
                                          p_Kld, p_KldTmp, &
                                          Itr1, Itr2, NEQ)        
            CASE DEFAULT
              PRINT *,'lsyssc_sortMatrix: Unsupported data type.'
              CALL sys_halt()
            END SELECT
            
            ! Remove temp matrix
            CALL lsyssc_releaseMatrix(rtempMatrix)
            
          END IF
            
        END IF

      CASE (LSYSSC_MATRIXD)
        
        ! D-matrices can only be sorted if we are allowed to sort the content --
        ! as they have no explicit structure!
        
        IF (bsortEntries .AND. &
            (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .EQ. 0)) THEN

          ! Duplicate the matrix before resorting.
          CALL lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                       LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
        
          ! Sort entries; there is no structure.
          SELECT CASE (rmatrix%cdataType)
          CASE (ST_DOUBLE)
            ! Double precision version
            CALL lsyssc_getbase_double (rmatrix,p_Ddata)
            CALL lsyssc_getbase_double (rtempMatrix,p_DdataTmp)
            CALL lalg_vectorSortDble (p_DdataTmp, p_Ddata, Itr1)
          CASE (ST_SINGLE)
            ! Single precision version
            CALL lsyssc_getbase_single (rmatrix,p_Fdata)
            CALL lsyssc_getbase_single (rtempMatrix,p_FdataTmp)
            CALL lalg_vectorSortSngl (p_FdataTmp, p_Fdata, Itr1)
          CASE DEFAULT
            PRINT *,'lsyssc_sortMatrix: Unsupported data type.'
            CALL sys_halt()
          END SELECT
          
          ! Remove temp matrix
          CALL lsyssc_releaseMatrix(rtempMatrix)

        END IF

      CASE DEFAULT
        PRINT *,'lsyssc_sortMatrix: Unsupported matrix format!'
        CALL sys_halt()
        
      END SELECT

    END SUBROUTINE
    
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat7_double (Da, DaH, Icol, IcolH, &
                                     Ild, IldH, Itr1, Itr2, neq)
  
!<description>
    ! Resorts structure and entries of the given matrix, corresponding to Itr1/Itr2.
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

  SUBROUTINE lsyssc_sortMat7_single (Fa, FaH, Icol, IcolH, &
                                     Ild, IldH, Itr1, Itr2, neq)
  
!<description>
    ! Resorts structure and entries of the given matrix, corresponding to Itr1/Itr2.
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

  SUBROUTINE lsyssc_sortMat7_int (Ia, IaH, Icol, IcolH, &
                                     Ild, IldH, Itr1, Itr2, neq)
  
!<description>
    ! Resorts structure and entries of the given matrix, corresponding to Itr1/Itr2.
    ! Integer version.
    !
    ! Storage technique 7 version.
!</description>
    
!<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Source matrix
    INTEGER, DIMENSION(:), INTENT(IN) :: IaH
    
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
    INTEGER, DIMENSION(:), INTENT(OUT) :: Ia

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
      Ia(ildIdx) = IaH(IldH(iidx))
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
        Ia(ildIdx) = IaH(Ih2(j))

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

  SUBROUTINE lsyssc_sortMat7Ent_double (Da, DaH, IcolH, &
                                        IldH, Itr1, Itr2, neq)
  
!<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! The structure is not changed.
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
      ! The new row starts at position ildIdx. 
      Da(ildIdx) = DaH(IldH(iidx))
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

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO

    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat7Ent_single (Fa, FaH, IcolH, &
                                        IldH, Itr1, Itr2, neq)
  
!<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! The structure is not changed.
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
      ! The new row starts at position ildIdx. 
      Fa(ildIdx) = FaH(IldH(iidx))
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

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO

    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat7Ent_int (Ia, IaH, IcolH, &
                                     IldH, Itr1, Itr2, neq)
  
!<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! The structure is not changed.
    ! Integer version.
    !
    ! Storage technique 7 version.
!</description>
    
!<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Source matrix
    INTEGER, DIMENSION(:), INTENT(IN) :: IaH
    
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
    INTEGER, DIMENSION(:), INTENT(OUT) :: Ia

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
      ! The new row starts at position ildIdx. 
      Ia(ildIdx) = IaH(IldH(iidx))
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
        Ia(ildIdx) = IaH(Ih2(j))

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO

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
    ! Resorts structure and entries of the given matrix, corresponding to Itr1/Itr2.
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
    ! Resorts structure and entries of the given matrix, corresponding to Itr1/Itr2.
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

  SUBROUTINE lsyssc_sortMat9_int (Ia, IaH, Icol, IcolH, &
                                     Ild, IldH, Idiag, Itr1, Itr2, neq)
  
!<description>
    ! Resorts structure and entries of the given matrix, corresponding to Itr1/Itr2.
    ! Double precision version.
    !
    ! Storage technique 9 version.
!</description>
    
!<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Source matrix
    INTEGER, DIMENSION(:), INTENT(IN) :: IaH
    
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
    INTEGER, DIMENSION(:), INTENT(OUT) :: Ia

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
        Ia(ildIdx) = IaH(Ih2(j))
        
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

  SUBROUTINE lsyssc_sortMat9Ent_double (Da, DaH, IcolH, IldH, &
                                        Itr1, Itr2, neq)
  
!<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! The structure is not changed.
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
      
      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j = 1,isize
        
        ! Get the matrix entry of the current column and write it to DA
        Da(ildIdx) = DaH(Ih2(j))
        
        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO
    
    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat9Ent_single (Fa, FaH, IcolH, IldH, &
                                        Itr1, Itr2, neq)
  
!<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! The structure is not changed.
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
      
      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j = 1,isize
        
        ! Get the matrix entry of the current column and write it to DA
        Fa(ildIdx) = FaH(Ih2(j))
        
        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO
    
    ! Release temp variables
    DEALLOCATE(Ih1,Ih2)

  END SUBROUTINE 

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_sortMat9Ent_int (Ia, IaH, IcolH, IldH, &
                                     Itr1, Itr2, neq)
  
!<description>
    ! Resorts the entries of the given matrix, corresponding to Itr1/Itr2.
    ! The structure is not changed.
    ! Double precision version.
    !
    ! Storage technique 9 version.
!</description>
    
!<input>
    
    ! Number of equations
    INTEGER(PREC_VECIDX), INTENT(IN) :: neq
    
    ! Source matrix
    INTEGER, DIMENSION(:), INTENT(IN) :: IaH
    
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
    INTEGER, DIMENSION(:), INTENT(OUT) :: Ia

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
      
      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      DO j = 1,isize
        
        ! Get the matrix entry of the current column and write it to DA
        Ia(ildIdx) = IaH(Ih2(j))
        
        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      END DO
    END DO
    
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
  
  SUBROUTINE lsyssc_unsortMatrix (rmatrix,bsortEntries)
  
!<description>
  ! This routine deactivates the sorting of matrix rmatrix, the matrix
  ! is unsorted. If bsortEntries=TRUE, the entries of the matrix are unsorted
  ! together with the structure; otherwise, only the structure is unsorted,
  ! thus leaving the entries in an undefined state.
  ! rmatrix%isortStrategy is set to negative to indicate that the sorting
  ! is deactivated.
!</description>
  
!<inputoutput>
  ! Vector to resort
  TYPE(t_matrixScalar), INTENT(INOUT), TARGET       :: rmatrix
!</inputoutput>

!<input>
  ! Sort the entries or only the structure of the matrix.
  ! = FALSE: Only sort the structure of the matrix,
  ! = TRUE : Sort both, entries and structure of the matrix.
  LOGICAL, INTENT(IN) :: bsortEntries

!</input>

!</subroutine>

    ! Call the sort-routine, deactivate the sorting -- if activated.
    CALL lsyssc_sortMatrix (rmatrix,bsortEntries,-ABS(rmatrix%isortStrategy))

  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_addIndex (h_Ix,ivalue,istartpos,ilength)
  
!<description>

  ! This is an auxiliary routine. It accepts a handle to an integer array
  ! and adds the value ivalue to each entry of this vector.
  ! This can be used e.g. to convert the index vector of a matrix
  ! from 1-based to 0-based for an external library (UMFPACK).

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
    CALL sys_halt()
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
    CALL sys_halt()
  END SELECT
  
  END FUNCTION

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_invertedDiagMatVec (rmatrix,rvectorSrc,dscale,rvectorDst)
  
!<description>
  ! This routine multiplies the weighted inverted diagonal $domega*D^{-1}$
  ! of the matrix rmatrix with the vector rvectorSrc and stores the result 
  ! into the vector rvectorDst:
  !   $$rvectorDst = dscale * D^{-1} * rvectorSrc$$
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
  INTEGER :: h_Idiag,ivar,NVAR
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
    CALL sys_halt()
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
      CALL lsyssc_getbase_double (rmatrix,p_Da)
      dmyscale = dscale * rmatrix%dscaleFactor

      ! And the vector(s)?
      SELECT CASE (rvectorSrc%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rvectorSrc,p_Dvec)
        CALL lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let's go...
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i)
        DO i=1,rvectorSrc%NEQ
          p_Dvec2(i) = p_Dvec(i)*dmyscale/p_Da(p_Kdiag(i))
        END DO
!$omp end parallel do
        
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rvectorSrc,p_Fvec)
        CALL lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let's go...
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i)
        DO i=1,rvectorSrc%NEQ
          p_Fvec2(i) = p_Fvec(i)*dmyscale/p_Da(p_Kdiag(i))
        END DO
!$omp end parallel do
        
      CASE DEFAULT
        PRINT *,'lsyssc_invertedDiagMatVec: unsupported vector precision!'
        CALL sys_halt()
      END SELECT
      
    CASE (ST_SINGLE)
      CALL lsyssc_getbase_single (rmatrix,p_Fa)
      fmyscale = REAL(dscale * rmatrix%dscaleFactor,SP)

      ! And the vector(s)?
      SELECT CASE (rvectorSrc%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rvectorSrc,p_Dvec)
        CALL lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let's go...
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i)
        DO i=1,rvectorSrc%NEQ
          p_Dvec2(i) = p_Dvec(i)*fmyscale/p_Fa(p_Kdiag(i))
        END DO
!$omp end parallel do
        
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rvectorSrc,p_Fvec)
        CALL lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let's go...
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i)
        DO i=1,rvectorSrc%NEQ
          p_Fvec2(i) = p_Fvec(i)*fmyscale/p_Fa(p_Kdiag(i))
        END DO
!$omp end parallel do
        
      CASE DEFAULT
        PRINT *,'lsyssc_invertedDiagMatVec: unsupported vector precision!'
        CALL sys_halt()
      END SELECT

    CASE DEFAULT
      PRINT *,'lsyssc_invertedDiagMatVec: unsupported matrix precision!'
      CALL sys_halt()
    END SELECT

  CASE (LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7INTL)
    ! Get a pointer to the diagonal - either h_KLD or h_Kdiagonal
    IF (rmatrix%cmatrixFormat .EQ. LSYSSC_MATRIX9INTL) THEN
      h_Idiag = rmatrix%h_Kdiagonal
    ELSE
      h_Idiag = rmatrix%h_Kld
    END IF
    CALL storage_getbase_int(h_Idiag,p_Kdiag)
    
    NVAR=rmatrix%NVAR

    ! Data type?
    SELECT CASE (rmatrix%cdataType)
    CASE (ST_DOUBLE)
      CALL lsyssc_getbase_double (rmatrix,p_Da)
      dmyscale = dscale * rmatrix%dscaleFactor

      ! And the vector(s)?
      SELECT CASE (rvectorSrc%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rvectorSrc,p_Dvec)
        CALL lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let's go...
        SELECT CASE(rmatrix%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIX1)
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i,ivar)
          DO i=1,rvectorSrc%NEQ
            DO ivar=1,NVAR
              p_Dvec2(NVAR*(i-1)+ivar) = p_Dvec(NVAR*(i-1)+ivar)&
                  *dmyscale/p_Da(NVAR*NVAR*(p_Kdiag(i)-1)+NVAR*ivar)
            END DO
          END DO
!$omp end parallel do

        CASE (LSYSSC_MATRIXD)
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i,ivar)
          DO i=1,rvectorSrc%NEQ
            DO ivar=1,NVAR
              p_Dvec2(NVAR*(i-1)+ivar) = p_Dvec(NVAR*(i-1)+ivar)&
                  *dmyscale/p_Da(NVAR*(p_Kdiag(i)-1)+ivar)
            END DO
          END DO
!$omp end parallel do

        CASE DEFAULT
          PRINT *, 'lsyssc_invertedDiagMatVec: unsupported interleave &
              &matrix format!'
          CALL sys_halt()
        END SELECT
        
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rvectorSrc,p_Fvec)
        CALL lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let's go...
        SELECT CASE(rmatrix%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIX1)
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i,ivar)
          DO i=1,rvectorSrc%NEQ
            DO ivar=1,NVAR
              p_Fvec2(NVAR*(i-1)+ivar) = p_Fvec(NVAR*(i-1)+ivar)&
                  *dmyscale/p_Da(NVAR*NVAR*(p_Kdiag(i)-1)+NVAR*ivar)
            END DO
          END DO
!$omp end parallel do

        CASE (LSYSSC_MATRIXD)
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i,ivar)
          DO i=1,rvectorSrc%NEQ
            DO ivar=1,NVAR
              p_Fvec2(NVAR*(i-1)+ivar) = p_Fvec(NVAR*(i-1)+ivar)&
                  *dmyscale/p_Da(NVAR*(p_Kdiag(i)-1)+ivar)
            END DO
          END DO
!$omp end parallel do

        CASE DEFAULT
          PRINT *, 'lsyssc_invertedDiagMatVec: unsupported interleave &
              &matrix format!'
          CALL sys_halt()
        END SELECT
        
      CASE DEFAULT
        PRINT *,'lsyssc_invertedDiagMatVec: unsupported vector precision!'
        CALL sys_halt()
      END SELECT
      
    CASE (ST_SINGLE)
      CALL lsyssc_getbase_single (rmatrix,p_Fa)
      fmyscale = REAL(dscale * rmatrix%dscaleFactor,SP)

      ! And the vector(s)?
      SELECT CASE (rvectorSrc%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rvectorSrc,p_Dvec)
        CALL lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let's go...
        SELECT CASE(rmatrix%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIX1)
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i,ivar)
          DO i=1,rvectorSrc%NEQ
            DO ivar=1,NVAR
              p_Dvec2(NVAR*(i-1)+ivar) = p_Dvec(NVAR*(i-1)+ivar)&
                  *fmyscale/p_Fa(NVAR*NVAR*(p_Kdiag(i)-1)+NVAR*ivar)
            END DO
          END DO
!$omp end parallel do

        CASE (LSYSSC_MATRIXD)
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i,ivar)
          DO i=1,rvectorSrc%NEQ
            DO ivar=1,NVAR
              p_Dvec2(NVAR*(i-1)+ivar) = p_Dvec(NVAR*(i-1)+ivar)&
                  *fmyscale/p_Fa(NVAR*(p_Kdiag(i)-1)+ivar)
            END DO
          END DO
!$omp end parallel do
        
        CASE DEFAULT
          PRINT *, 'lsyssc_invertedDiagMatVec: unsupported interleave &
              &matrix format!'
          CALL sys_halt()
        END SELECT
        
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rvectorSrc,p_Fvec)
        CALL lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let's go...
        SELECT CASE(rmatrix%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIX1)
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i,ivar)
          DO i=1,rvectorSrc%NEQ
            DO ivar=1,NVAR
              p_Fvec2(NVAR*(i-1)+ivar) = p_Fvec(NVAR*(i-1)+ivar)&
                  *fmyscale/p_Fa(NVAR*NVAR*(p_Kdiag(i)-1)+NVAR*ivar)
            END DO
          END DO
!$omp end parallel do

        CASE (LSYSSC_MATRIXD)
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i,ivar)
          DO i=1,rvectorSrc%NEQ
            DO ivar=1,NVAR
              p_Fvec2(NVAR*(i-1)+ivar) = p_Fvec(NVAR*(i-1)+ivar)&
                  *fmyscale/p_Fa(NVAR*(p_Kdiag(i)-1)+ivar)
            END DO
          END DO
!$omp end parallel do
        
        CASE DEFAULT
          PRINT *, 'lsyssc_invertedDiagMatVec: unsupported interleave &
              &matrix format!'
          CALL sys_halt()
        END SELECT

      CASE DEFAULT
        PRINT *,'lsyssc_invertedDiagMatVec: unsupported vector precision!'
        CALL sys_halt()
      END SELECT

    CASE DEFAULT
      PRINT *,'lsyssc_invertedDiagMatVec: unsupported matrix precision!'
      CALL sys_halt()
    END SELECT
    
  CASE (LSYSSC_MATRIXD)
    
    ! Data type?
    SELECT CASE (rmatrix%cdataType)
    CASE (ST_DOUBLE)
      CALL lsyssc_getbase_double (rmatrix,p_Da)
      dmyscale = dscale * rmatrix%dscaleFactor

      ! And the vector(s)?
      SELECT CASE (rvectorSrc%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rvectorSrc,p_Dvec)
        CALL lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let's go...
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i)
        DO i=1,rvectorSrc%NEQ
          p_Dvec2(i) = p_Dvec(i)*dmyscale/p_Da(i)
        END DO
!$omp end parallel do
        
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rvectorSrc,p_Fvec)
        CALL lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let's go...
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i)
        DO i=1,rvectorSrc%NEQ
          p_Fvec2(i) = p_Fvec(i)*dmyscale/p_Da(i)
        END DO
!$omp end parallel do
        
      CASE DEFAULT
        PRINT *,'lsyssc_invertedDiagMatVec: unsupported vector precision!'
        CALL sys_halt()
      END SELECT
      
    CASE (ST_SINGLE)
      CALL lsyssc_getbase_single (rmatrix,p_Fa)
      fmyscale = REAL(dscale * rmatrix%dscaleFactor,SP)

      ! And the vector(s)?
      SELECT CASE (rvectorSrc%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rvectorSrc,p_Dvec)
        CALL lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let's go...
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i)
        DO i=1,rvectorSrc%NEQ
          p_Dvec2(i) = p_Dvec(i)*fmyscale/p_Fa(i)
        END DO
!$omp end parallel do
        
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rvectorSrc,p_Fvec)
        CALL lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let's go...
!$omp parallel do&
!$omp&default(shared) &
!$omp&private(i)
        DO i=1,rvectorSrc%NEQ
          p_Fvec2(i) = p_Fvec(i)*fmyscale/p_Fa(i)
        END DO
!$omp end parallel do
        
      CASE DEFAULT
        PRINT *,'lsyssc_invertedDiagMatVec: unsupported vector precision!'
        CALL sys_halt()
      END SELECT

    CASE DEFAULT
      PRINT *,'lsyssc_invertedDiagMatVec: unsupported matrix precision!'
      CALL sys_halt()
    END SELECT

  CASE DEFAULT
    PRINT *,'lsyssc_invertedDiagMatVec: unsupported matrix format!'
    CALL sys_halt()
  END SELECT
    
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_copyVector (rx,ry)
  
!<description>
  ! Copies vector data: ry = rx.
  ! Both vectors must have the same size. All structural data of rx is
  ! transferred to ry, so rx and ry are compatible to each other afterwards.
  ! If ry is empty, new memory is allocated automatically. Otherwise,
  ! ry is overwritten.
!</description>

!<input>
  ! Source vector
  TYPE(t_vectorScalar),INTENT(IN) :: rx
!</input>

!<inputoutput>
  ! Destination vector
  TYPE(t_vectorScalar),INTENT(INOUT) :: ry
!</inputoutput>
  
!</subroutine>

    CALL lsyssc_duplicateVector (rx,ry,&
        LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE)

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
    CALL sys_halt()
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_clearVector (rx,dvalue)
  
!<description>
  ! Clears the block vector dx: Dx = 0 (or Dx = dvalue if dvalue is specified)
!</description>
  
!<inputoutput>
  ! Destination vector to be cleared
  TYPE(t_vectorScalar), INTENT(INOUT) :: rx

  ! OPTIONAL: Value to write into the matrix.
  ! If not specified, all matrix entries are set to 0.0.
  ! If specified, all matrix entries are set to dvalue.
  REAL(DP), INTENT(IN), OPTIONAL :: dvalue
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
    IF (.NOT. PRESENT(dvalue)) THEN
      CALL lalg_clearVectorDble (p_Dsource)
    ELSE
      CALL lalg_setVectorDble (p_Dsource,dvalue)
    END IF
  
  CASE (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    CALL lsyssc_getbase_single(rx,p_Ssource)
    IF (.NOT. PRESENT(dvalue)) THEN
      CALL lalg_clearVectorSngl (p_Ssource)
    ELSE
      CALL lalg_setVectorSngl (p_Ssource,REAL(dvalue,SP))
    END IF

  CASE DEFAULT
    PRINT *,'lsyssc_clearVector: Unsupported data type!'
    CALL sys_halt()
  END SELECT
  
  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_vectorLinearComb (rx,ry,cx,cy,rdest)
  
!<description>
  ! Performs a linear combination: rdest or ry = cx * rx  +  cy * ry
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
  ! Second source vector; also receives the result if rdest is not specified.
  TYPE(t_vectorScalar), INTENT(INOUT), TARGET :: ry

  ! OPTIONAL: Destination vector. If not specified, ry will be overwritten.
  TYPE(t_vectorScalar), INTENT(INOUT), OPTIONAL, TARGET :: rdest
!</inputoutput>
  
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Dsource, p_Ddest
  REAL(SP), DIMENSION(:), POINTER :: p_Ssource, p_Sdest
  TYPE(t_vectorScalar), POINTER :: p_rdest
  
  ! The vectors must be compatible to each other.
  CALL lsyssc_isVectorCompatible (rx,ry)

  IF (rx%cdataType .NE. ry%cdataType) THEN
    PRINT *,'lsyssc_vectorLinearComb: different data types not supported!'
    CALL sys_halt()
  END IF
  
  p_rdest => ry
  IF (PRESENT(rdest)) THEN
    
    p_rdest => rdest
    
    ! The vectors must be compatible to each other or rdest must be a clean
    ! vector.
    IF (rdest%NEQ .NE. 0) THEN
      CALL lsyssc_isVectorCompatible (rx,rdest)

      IF (rx%cdataType .NE. rdest%cdataType) THEN
        PRINT *,'lsyssc_vectorLinearComb: different data types not supported!'
        CALL sys_halt()
      END IF
    END IF
    
    ! Copy ry to rdest. Change rdest afterwards.
    CALL lsyssc_copyVector (ry,rdest)
    
  END IF
  
  SELECT CASE (rx%cdataType)
  CASE (ST_DOUBLE)
    ! Get the pointers and copy the whole data array.
    CALL lsyssc_getbase_double(rx,p_Dsource)
    CALL lsyssc_getbase_double(p_rdest,p_Ddest)
    
    CALL lalg_vectorLinearCombDble (p_Dsource,p_Ddest,cx,cy)

  CASE (ST_SINGLE)
    ! Get the pointers and copy the whole data array.
    CALL lsyssc_getbase_single(rx,p_Ssource)
    CALL lsyssc_getbase_single(p_rdest,p_Sdest)
    
    CALL lalg_vectorLinearCombSngl (p_Ssource,p_Sdest,REAL(cx,SP),REAL(cy,SP))
  
  CASE DEFAULT
    PRINT *,'lsyssc_vectorLinearComb: Unsupported data type!'
    CALL sys_halt()
  END SELECT

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_copyMatrix (rsourceMatrix,rdestMatrix)
  
!<description>
  ! Copies a matrix data: rsourceMatrix = rdestMatrix.
  ! All structural (discretisation related) data of rsourceMatrix is
  ! transferred to rdestMatrix.
  !
  ! If the destination matrix contains data, the data is overwritten,
  ! regardless of whether the data belongs to rdestmatrix or not.
  ! If the destination matrix is empty, new memory is allocated.
  !
  ! Therefore, lsyssc_copyMatrix coincides with a call to
  ! lsyssc_duplicateMatrix (.,.,LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE).
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

    CALL lsyssc_duplicateMatrix (rsourceMatrix,rdestMatrix,&
        LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE)

!  ! local variables
!  TYPE(t_matrixScalar) :: rtempMatrix
!  REAL(DP), DIMENSION(:), POINTER :: p_Dsource,p_Ddest
!  REAL(SP), DIMENSION(:), POINTER :: p_Fsource,p_Fdest
!  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_KcolSource,p_KcolDest
!  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KldSource,p_KldDest
!  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_KdiagonalSource,p_KdiagonalDest
!  
!  ! Is it possible at all to copy the matrix? Both matrices must have
!  ! the same size, otherwise the memory does not fit.
!  IF (rsourceMatrix%cmatrixFormat .NE. rdestMatrix%cmatrixFormat) THEN
!    PRINT *,'lsyssc_copyMatrix: Different matrix formats not allowed!'
!    CALL sys_halt()
!  END IF
!  
!  IF (rsourceMatrix%NA .NE. rdestMatrix%NA) THEN
!    PRINT *,'lsyssc_copyMatrix: Matrices have different size!'
!    CALL sys_halt()
!  END IF
!  
!  ! NEQ/NCOLS is irrelevant. It's only important that we have enough memory
!  ! and the same structure!
!  
!  IF (rsourceMatrix%cdataType .NE. rdestMatrix%cdataType) THEN
!    PRINT *,'lsyssc_copyMatrix: Matrices have different data types!'
!    CALL sys_halt()
!  END IF
!  
!  ! First, make a backup of the matrix for restoring some critical data.
!  rtempMatrix = rdestMatrix
!  
!  ! Copy the matrix structure
!  rdestMatrix = rsourceMatrix
!  
!  ! Restore crucial data
!  rdestMatrix%imatrixSpec = rtempMatrix%imatrixSpec
!
!  ! Which structure do we actually have?
!  SELECT CASE (rsourceMatrix%cmatrixFormat)
!  CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
!  
!    ! Restore crucial format-specific data
!    rdestMatrix%h_DA        = rtempMatrix%h_DA
!    rdestMatrix%h_Kcol      = rtempMatrix%h_Kcol
!    rdestMatrix%h_Kld       = rtempMatrix%h_Kld
!    rdestMatrix%h_Kdiagonal = rtempMatrix%h_Kdiagonal
!    
!    ! And finally copy the data. 
!    SELECT CASE (rsourceMatrix%cdataType)
!    
!    CASE (ST_DOUBLE)
!      CALL lsyssc_getbase_double (rsourceMatrix,p_Dsource)
!      CALL lsyssc_getbase_double (rdestMatrix,p_Ddest)
!      CALL lalg_copyVectorDble (p_Dsource,p_Ddest)
!
!    CASE (ST_SINGLE)
!      CALL lsyssc_getbase_single (rsourceMatrix,p_Fsource)
!      CALL lsyssc_getbase_single (rdestMatrix,p_Fdest)
!      CALL lalg_copyVectorSngl (p_Fsource,p_Fdest)
!
!    CASE DEFAULT
!      PRINT *,'lsyssc_copyMatrix: Unsupported data type!'
!      CALL sys_halt()
!    END SELECT
!    
!    CALL lsyssc_getbase_Kcol (rsourceMatrix,p_KcolSource)
!    CALL lsyssc_getbase_Kcol (rdestMatrix,p_KcolDest)
!    CALL lalg_copyVectorInt (p_KcolSource,p_KcolDest)
!
!    CALL lsyssc_getbase_Kld (rsourceMatrix,p_KldSource)
!    CALL lsyssc_getbase_Kld (rdestMatrix,p_KldDest)
!    CALL lalg_copyVectorInt (p_KldSource,p_KldDest)
!    
!    CALL lsyssc_getbase_Kdiagonal (rsourceMatrix,p_KdiagonalSource)
!    CALL lsyssc_getbase_Kdiagonal (rdestMatrix,p_KdiagonalDest)
!    CALL lalg_copyVectorInt (p_KdiagonalSource,p_KdiagonalDest)
!      
!  CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
!
!    ! Restore crucial format-specific data
!    rdestMatrix%h_DA        = rtempMatrix%h_DA
!    rdestMatrix%h_Kcol      = rtempMatrix%h_Kcol
!    rdestMatrix%h_Kld       = rtempMatrix%h_Kld
!    
!    ! And finally copy the data. 
!    SELECT CASE (rsourceMatrix%cdataType)
!    CASE (ST_DOUBLE)
!      CALL lsyssc_getbase_double (rsourceMatrix,p_Dsource)
!      CALL lsyssc_getbase_double (rdestMatrix,p_Ddest)
!      CALL lalg_copyVectorDble (p_Dsource,p_Ddest)
!
!    CASE (ST_SINGLE)
!      CALL lsyssc_getbase_single (rsourceMatrix,p_Fsource)
!      CALL lsyssc_getbase_single (rdestMatrix,p_Fdest)
!      CALL lalg_copyVectorSngl (p_Fsource,p_Fdest)
!
!    CASE DEFAULT
!      PRINT *,'storage_copyMatrix: Unsupported data type!'
!      CALL sys_halt()
!    END SELECT
!
!    CALL lsyssc_getbase_Kcol (rsourceMatrix,p_KcolSource)
!    CALL lsyssc_getbase_Kcol (rdestMatrix,p_KcolDest)
!    CALL lalg_copyVectorInt (p_KcolSource,p_KcolDest)
!
!    CALL lsyssc_getbase_Kld (rsourceMatrix,p_KldSource)
!    CALL lsyssc_getbase_Kld (rdestMatrix,p_KldDest)
!    CALL lalg_copyVectorInt (p_KldSource,p_KldDest)
!
!  CASE (LSYSSC_MATRIXD)
!
!    ! Restore crucial format-specific data
!    rdestMatrix%h_DA        = rtempMatrix%h_DA
!    
!    ! And finally copy the data. 
!    SELECT CASE (rsourceMatrix%cdataType)
!    CASE (ST_DOUBLE)
!      CALL lsyssc_getbase_double (rsourceMatrix,p_Dsource)
!      CALL lsyssc_getbase_double (rdestMatrix,p_Ddest)
!      CALL lalg_copyVectorDble (p_Dsource,p_Ddest)
!
!    CASE (ST_SINGLE)
!      CALL lsyssc_getbase_single (rsourceMatrix,p_Fsource)
!      CALL lsyssc_getbase_single (rdestMatrix,p_Fdest)
!      CALL lalg_copyVectorSngl (p_Fsource,p_Fdest)
!
!    CASE DEFAULT
!      PRINT *,'storage_copyMatrix: Unsupported data type!'
!      CALL sys_halt()
!    END SELECT
!
!  CASE DEFAULT
!    PRINT *,'lsyssc_copyMatrix: Unsupported matrix format!'
!    CALL sys_halt()
!  END SELECT
  
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
      CALL sys_halt()
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

  SUBROUTINE lsyssc_transposeMatrixInSitu (rmatrix,itransFlag)
  
!<description>
  ! This routine in-situ transposes a matrix rmatrix and creates the transposed
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
!</input>

!<inputoutput>
  ! The matrix to be transposed.
  TYPE(t_matrixScalar),INTENT(INOUT) :: rmatrix
!</inputoutput>
!</subroutine>

    TYPE(t_matrixScalar) :: rmatrix2
    
    CALL lsyssc_transposeMatrix (rmatrix,rmatrix2,itransFlag)

    ! Transfer the copy-flags
    rmatrix2%imatrixSpec = IOR(IAND(rmatrix2%imatrixSpec,NOT(LSYSSC_MSPEC_ISCOPY)),&
                               IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY))

    ! Probably remove old data
    IF (.NOT. lsyssc_isMatrixContentShared(rmatrix,rmatrix2)) &
      CALL lsyssc_releaseMatrixContent(rmatrix)
    
    ! Prevent lsyssc_releaseMatrix from releasing the content if it's still there
    rmatrix%imatrixSpec = IOR(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY)
    
    IF (.NOT. lsyssc_isMatrixStructureShared(rmatrix,rmatrix2)) &
      CALL lsyssc_releaseMatrix(rmatrix)
      
    ! Create a hard-copy of rmatrix.
    ! This is one of the very rare cases where we on purpose assign a matrix
    ! structure to another...
    rmatrix = rmatrix2
    
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
        CALL sys_halt()
      END IF

      IF (rmatrix%NA .NE. rtransposedMatrix%NA) THEN
        PRINT *,'lsyssc_transposeMatrix: Matrices have different size!'
        CALL sys_halt()
      END IF

      ! NEQ/NCOLS is irrelevant. It's only important that we have enough memory
      ! and the same structure!

      IF (itrans .EQ. LSYSSC_TR_ALL) THEN
        IF (rmatrix%cdataType .NE. rtransposedMatrix%cdataType) THEN
          PRINT *,'lsyssc_transposeMatrix: Matrices have different data types!'
          CALL sys_halt()
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
        CALL sys_halt()
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
        
          CALL lsyssc_auxcopy_Kcol (rmatrix,rtransposedMatrix)
          CALL lsyssc_auxcopy_Kld (rmatrix,rtransposedMatrix)
          CALL lsyssc_auxcopy_Kdiagonal (rmatrix,rtransposedMatrix)
          
          rtransposedMatrix%imatrixSpec = &
            IAND(rtransposedMatrix%imatrixSpec,NOT(LSYSSC_MSPEC_TRANSPOSED))
          
        ELSE
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          CALL storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                            ST_INT, rtransposedMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
          CALL storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1_I32, &
                            ST_INT, rtransposedMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
          
          ! Get Kcol/Kld. Don't use the lsyssc_getbase routines as we just created
          ! the arrays for a matrix which is not yet valid!
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
          CALL lsyssc_auxcopy_Kcol (rmatrix,rtransposedMatrix)
          CALL lsyssc_auxcopy_Kld (rmatrix,rtransposedMatrix)
          rtransposedMatrix%imatrixSpec = &
            IAND(rtransposedMatrix%imatrixSpec,NOT(LSYSSC_MSPEC_TRANSPOSED))
          
        ELSE
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          CALL storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                            ST_INT, rtransposedMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
          CALL storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1_I32, &
                            ST_INT, rtransposedMatrix%h_Kld,ST_NEWBLOCK_NOINIT)

          ! Get Kcol/Kld. Don't use the lsyssc_getbase routines as we just created
          ! the arrays for a matrix which is not yet valid!
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
        
      CASE (LSYSSC_MATRIXD)
        ! Nothing to do

      CASE DEFAULT
        PRINT *,'lsyssc_transposeMatrix: Unsupported matrix format.'
        CALL sys_halt()
      END SELECT
      
    CASE (LSYSSC_TR_ALL)
      ! Transpose the full matrix - structure+content
      SELECT CASE (rmatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX9)
        ! Is the source matrix saved tranposed? In that case simply copy
        ! the data arrays, this results in the transposed matrix :)
        IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .NE. 0) THEN
        
          IF (rmatrix%h_Da .NE. ST_NOHANDLE) THEN
            CALL lsyssc_auxcopy_da (rmatrix,rtransposedMatrix)
          END IF
          CALL lsyssc_auxcopy_Kcol (rmatrix,rtransposedMatrix)
          CALL lsyssc_auxcopy_Kld (rmatrix,rtransposedMatrix)
          CALL lsyssc_auxcopy_Kdiagonal (rmatrix,rtransposedMatrix)
          rtransposedMatrix%imatrixSpec = &
            IAND(rtransposedMatrix%imatrixSpec,NOT(LSYSSC_MSPEC_TRANSPOSED))
          
        ELSE
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          CALL storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                            ST_INT, rtransposedMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
          CALL storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1_I32, &
                            ST_INT, rtransposedMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
          CALL storage_new ('lsyssc_transposeMatrix', 'Da', rmatrix%NA, &
                            rtransposedMatrix%cdataType, rtransposedMatrix%h_Da,&
                            ST_NEWBLOCK_NOINIT)

          ! Get Kcol/Kld. Don't use the lsyssc_getbase routines as we just created
          ! the arrays for a matrix which is not yet valid!
          CALL storage_getbase_int(rmatrix%h_Kcol,p_KcolSource)
          CALL storage_getbase_int(rmatrix%h_Kld,p_KldSource)

          CALL storage_getbase_int(rtransposedMatrix%h_Kcol,p_KcolDest)
          CALL storage_getbase_int(rtransposedMatrix%h_Kld,p_KldDest)
          
          ! Get the data array(s)
          CALL lsyssc_getbase_double(rMatrix,p_DaSource)
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
        
          IF (rmatrix%h_Da .NE. ST_NOHANDLE) THEN
            CALL lsyssc_auxcopy_da (rmatrix,rtransposedMatrix)
          END IF

          CALL lsyssc_auxcopy_Kcol (rmatrix,rtransposedMatrix)
          CALL lsyssc_auxcopy_Kld (rmatrix,rtransposedMatrix)
          rtransposedMatrix%imatrixSpec = &
            IAND(rtransposedMatrix%imatrixSpec,NOT(LSYSSC_MSPEC_TRANSPOSED))
          
        ELSE
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          CALL storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                            ST_INT, rtransposedMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
          CALL storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1_I32, &
                            ST_INT, rtransposedMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
          CALL storage_new ('lsyssc_transposeMatrix', 'Da', rmatrix%NA, &
                            rtransposedMatrix%cdataType, rtransposedMatrix%h_Da,&
                            ST_NEWBLOCK_NOINIT)

          ! Get Kcol/Kld. Don't use the lsyssc_getbase routines as we just created
          ! the arrays for a matrix which is not yet valid!
          CALL storage_getbase_int(rmatrix%h_Kcol,p_KcolSource)
          CALL storage_getbase_int(rmatrix%h_Kld,p_KldSource)

          CALL storage_getbase_int(rtransposedMatrix%h_Kcol,p_KcolDest)
          CALL storage_getbase_int(rtransposedMatrix%h_Kld,p_KldDest)
          
          ! Get the data array(s)
          CALL lsyssc_getbase_double(rMatrix,p_DaSource)
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
          
          ! Pivotise the matrix to move the diagonal element to the front.
          CALL lsyssc_pivotiseMatrix7double (rmatrix%NEQ, &
              p_KcolDest, p_KldDest, p_DaDest)
        END IF
      
      CASE (LSYSSC_MATRIXD)
        ! Nothing to do

      CASE DEFAULT
        PRINT *,'lsyssc_transposeMatrix: Unsupported matrix format.'
        CALL sys_halt()
      END SELECT
      
    CASE DEFAULT ! = LSYSSC_TR_ALL 
    
    END SELECT
  
    ! Exchange NEQ and NCOLS as these always describe the actual matrix.
    ntemp = rtransposedMatrix%NEQ
    rtransposedMatrix%NEQ = rtransposedMatrix%NCOLS
    rtransposedMatrix%NCOLS = ntemp
    
    ! That's it.
  
  END SUBROUTINE

  !****************************************************************************
  
!<subroutine>
      
  SUBROUTINE lsyssc_auxcopy_DA (rsourceMatrix,rdestMatrix)
  
!<description>
  ! Auxiliary routine: 
  ! Copies the content of rsourceMatrix%DA to rdestMatrix%DA
  ! without additional checks.
  ! If the destination array does not exist, it's created.
!</description>

!<input>
  ! Source matrix
  TYPE(t_matrixScalar), INTENT(IN) :: rsourceMatrix
!</input>
 
!<inputoutput>
  ! Destination matrix
  TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
!</inputoutput>

!</subroutine>
  
    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: p_Da1,p_Da2
    REAL(SP), DIMENSION(:), POINTER :: p_Fa1,p_Fa2
  
    IF (rsourceMatrix%h_Da .EQ. ST_NOHANDLE) THEN
      PRINT *,'lsyssc_auxcopy_DA: Source matrix undefined!'
      CALL sys_halt()
    END IF

    IF (rdestMatrix%h_Da .EQ. ST_NOHANDLE) THEN
      CALL storage_new ('lsyssc_auxcopy_DA', 'DA', rsourceMatrix%NA, &
          rsourceMatrix%cdataType, rdestMatrix%h_Da,ST_NEWBLOCK_NOINIT)
    END IF

    SELECT CASE (rsourceMatrix%cdataType)
    CASE (ST_DOUBLE)
      CALL lsyssc_getbase_double(rsourceMatrix,p_Da1)
      CALL lsyssc_getbase_double(rdestMatrix,p_Da2)
      CALL lalg_copyVectorDble (p_Da1,p_Da2)
      
    CASE (ST_SINGLE)
      CALL lsyssc_getbase_single(rsourceMatrix,p_Fa1)
      CALL lsyssc_getbase_single(rdestMatrix,p_Fa2)
      CALL lalg_copyVectorSngl (p_Fa1,p_Fa2)

    CASE DEFAULT
      PRINT *,'lsyssc_transposeMatrix: Unsupported data type!'
      CALL sys_halt()
      
    END SELECT

  END SUBROUTINE
  
  !****************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_auxcopy_Kcol (rsourceMatrix,rdestMatrix)
  
!<description>
  ! Auxiliary routine: 
  ! Copies the content of rsourceMatrix%Kcol to rdestMatrix%Kcol
  ! without additional checks.
  ! If the destination array does not exist, it's created.
!</description>
  
!<input>
  ! Source matrix
  TYPE(t_matrixScalar), INTENT(IN) :: rsourceMatrix
!</input>
 
!<inputoutput>
  ! Destination matrix
  TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kcol1,p_Kcol2
  
    CALL lsyssc_getbase_Kcol (rsourceMatrix,p_Kcol1)
    IF (rdestMatrix%h_Kcol .EQ. ST_NOHANDLE) THEN
      CALL storage_new ('lsyssc_auxcopy_Kcol', 'Kcol', rsourceMatrix%NA, &
          ST_INT, rdestMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
    END IF
    CALL lsyssc_getbase_Kcol (rdestMatrix,p_Kcol2)
    
    CALL lalg_copyVectorInt (p_Kcol1,p_Kcol2)
    
  END SUBROUTINE
    
  !****************************************************************************
    
!<subroutine>

  SUBROUTINE lsyssc_auxcopy_Kld (rsourceMatrix,rdestMatrix)
  
!<description>
  ! Auxiliary routine: 
  ! Copies the content of rsourceMatrix%Kld to rdestMatrix%Kld
  ! without additional checks.
  ! If the destination array does not exist, it's created.
!</description>
  
!<input>
  ! Source matrix
  TYPE(t_matrixScalar), INTENT(IN) :: rsourceMatrix
!</input>
 
!<inputoutput>
  ! Destination matrix
  TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kld1,p_Kld2
  
    CALL lsyssc_getbase_Kld (rsourceMatrix,p_Kld1)
    IF (rdestMatrix%h_Kld .EQ. ST_NOHANDLE) THEN
      CALL storage_new ('lsyssc_auxcopy_Kld', 'Kld', &
          rsourceMatrix%NEQ+1_PREC_VECIDX, &
          ST_INT, rdestMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
    END IF
    CALL lsyssc_getbase_Kld (rdestMatrix,p_Kld2)
    
    CALL lalg_copyVectorInt (p_Kld1,p_Kld2)
    
  END SUBROUTINE
    
  !****************************************************************************
    
!<subroutine>

  SUBROUTINE lsyssc_auxcopy_Kdiagonal (rsourceMatrix,rdestMatrix)
  
!<description>
  ! Auxiliary routine: 
  ! Copies the content of rsourceMatrix%Kdiagonal to rdestMatrix%Kdiagonal
  ! without additional checks.
  ! If the destination array does not exist, it's created.
!</description>
  
!<input>
  ! Source matrix
  TYPE(t_matrixScalar), INTENT(IN) :: rsourceMatrix
!</input>
 
!<inputoutput>
  ! Destination matrix
  TYPE(t_matrixScalar), INTENT(INOUT) :: rdestMatrix
!</inputoutput>
  
!</subroutine>
  
    ! local variables
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kdiag1,p_Kdiag2
  
    CALL lsyssc_getbase_Kdiagonal (rsourceMatrix,p_Kdiag1)
    IF (rdestMatrix%h_Kdiagonal .EQ. ST_NOHANDLE) THEN
      CALL storage_new ('lsyssc_auxcopy_Kdiagonal', 'Kdiagonal', &
          rsourceMatrix%NEQ, &
          ST_INT, rdestMatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
    END IF
    CALL lsyssc_getbase_Kdiagonal (rdestMatrix,p_Kdiag2)
    
    CALL lalg_copyVectorInt (p_Kdiag1,p_Kdiag2)
    
  END SUBROUTINE
      
  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_allocEmptyMatrix (rmatrixScalar,iclear,cdataType)
  
!<description>
  ! This routine allocates memory for the matrix entries without computing
  ! the entries. This can be used to attach an 'empty' matrix to a matrix
  ! structure. The number of entries NA as well as the type of the matrix
  ! cmatrixFormat must be initialised in rmatrixScalar.
!</description>

!<input>
  ! Whether and how to fill the matrix with initial values.
  ! One of the LSYSSC_SETM_xxxx constants:
  ! LSYSSC_SETM_UNDEFINED : Don't initialise the matrix,
  ! LSYSSC_SETM_ZERO      : Clear the matrix / fill it with 0.0,
  ! LSYSSC_SETM_ONE       : Fill the matrix with 1.0. (Used e.g.
  !                         for UMFPACK who needs a non-zero
  !                         matrix for symbolic factorisation.)
  INTEGER, INTENT(IN) :: iclear
  
  ! OPTIONAL: Data type of the matrix (ST_SINGLE, ST_DOUBLE)
  ! If not present, the standard data type cdataType-variable in the matrix
  ! is used (usually ST_DOUBLE).
  INTEGER, INTENT(IN), OPTIONAL :: cdataType
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER(PREC_MATIDX) :: NA
  INTEGER :: cdType,NVAR
  REAL(SP), DIMENSION(:), POINTER :: p_Fa
  REAL(DP), DIMENSION(:), POINTER :: p_Da
  
  IF (rmatrixScalar%h_DA .NE. ST_NOHANDLE) THEN
    PRINT *,'lsyssc_allocEmptyMatrix: Cannot create empty matrix; exists already!'
    CALL sys_halt()
  END IF
  
  NA = rmatrixScalar%NA
  NVAR = rmatrixScalar%NVAR

  IF (PRESENT(cdataType)) THEN
    cdType = cdataType
    ! Change the data type in the matrix according to the new we will use.
    rmatrixScalar%cdataType = cdType
  ELSE
    cdType = rmatrixScalar%cdataType
  END IF

  ! Which matrix structure do we have?
  SELECT CASE (rmatrixScalar%cmatrixFormat) 
  CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX7,LSYSSC_MATRIXD,LSYSSC_MATRIX1)
    
    ! Check if the matrix entries exist. If not, allocate the matrix.
    IF (rmatrixScalar%h_DA .EQ. ST_NOHANDLE) THEN
    
      IF (iclear .GE. LSYSSC_SETM_ZERO) THEN
        CALL storage_new1D ('lsyssc_allocEmptyMatrix', 'DA', &
                            NA, cdType, rmatrixScalar%h_DA, &
                            ST_NEWBLOCK_ZERO)
      ELSE
        CALL storage_new1D ('lsyssc_allocEmptyMatrix', 'DA', &
                            NA, cdType, rmatrixScalar%h_DA, &
                            ST_NEWBLOCK_NOINIT)
      END IF
      
      IF (iclear .GE. LSYSSC_SETM_ONE) THEN
        SELECT CASE (cdType)
        CASE (ST_DOUBLE)
          CALL lsyssc_getbase_double (rmatrixScalar,p_Da)
          CALL lalg_setVectorDble (p_Da,1.0_DP)
        CASE (ST_SINGLE)
          CALL lsyssc_getbase_single (rmatrixScalar,p_Fa)
          CALL lalg_setVectorSngl (p_Fa,1.0_SP)
        CASE DEFAULT
          PRINT *,'lsyssc_allocEmptyMatrix: Unknown data type!'
          CALL sys_halt()
        END SELECT
      END IF
      
    END IF
    
  CASE (LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7INTL)
    
    ! Check if the matrix entries exist. If not, allocate the matrix.
    IF (rmatrixScalar%h_DA .EQ. ST_NOHANDLE) THEN
    
      IF (iclear .GE. LSYSSC_SETM_ZERO) THEN
        SELECT CASE (rmatrixScalar%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIX1)
          CALL storage_new1D ('lsyssc_allocEmptyMatrix', 'DA', &
              INT(NA*NVAR*NVAR,I32), cdType, rmatrixScalar%h_DA, ST_NEWBLOCK_ZERO)
        CASE (LSYSSC_MATRIXD)
          CALL storage_new1D ('lsyssc_allocEmptyMatrix', 'DA', &
              INT(NA*NVAR,I32), cdType, rmatrixScalar%h_DA, ST_NEWBLOCK_ZERO)
        CASE DEFAULT
          PRINT *, 'lsyssc_allocEmptyMatrix: Unsupported interl&
              &eave matrix format'
          CALL sys_halt()
        END SELECT

      ELSE
        SELECT CASE (rmatrixScalar%cinterleavematrixFormat)
        CASE (LSYSSC_MATRIX1)
          CALL storage_new1D ('lsyssc_allocEmptyMatrix', 'DA', &
              INT(NA*NVAR*NVAR,I32), cdType, rmatrixScalar%h_DA, ST_NEWBLOCK_NOINIT)
        CASE (LSYSSC_MATRIXD)
          CALL storage_new1D ('lsyssc_allocEmptyMatrix', 'DA', &
              INT(NA*NVAR,I32), cdType, rmatrixScalar%h_DA, ST_NEWBLOCK_NOINIT)
        CASE DEFAULT
          PRINT *, 'lsyssc_allocEmptyMatrix: Unsupported interl&
              &eave matrix format'
          CALL sys_halt()
        END SELECT
          
      END IF
      
      IF (iclear .GE. LSYSSC_SETM_ONE) THEN
        SELECT CASE (cdType)
        CASE (ST_DOUBLE)
          CALL lsyssc_getbase_double (rmatrixScalar,p_Da)
          CALL lalg_setVectorDble (p_Da,1.0_DP)
        CASE (ST_SINGLE)
          CALL lsyssc_getbase_single (rmatrixScalar,p_Fa)
          CALL lalg_setVectorSngl (p_Fa,1.0_SP)
        CASE DEFAULT
          PRINT *,'lsyssc_allocEmptyMatrix: Unknown data type!'
          CALL sys_halt()
        END SELECT
      END IF
      
    END IF
    
  CASE DEFAULT
    PRINT *,'lsyssc_allocEmptyMatrix: Not supported matrix structure!'
    CALL sys_halt()
  END SELECT
    
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_createDiagMatrixStruc (rmatrix,NEQ,cmatrixFormat)
  
!<description>
  ! Creates a diagonal matrix in matrix format cmatrixType. Initialises
  ! the structure but does not allocate any memory for the entries (if the
  ! matrix has any).
  !
  ! The matrix is created 'directly' without any discretisation structure
  ! in the background. NEQ specifies the size of the matrix.
!</description>

!<input>
  ! Matrix type, e.g. LSYSSC_MATRIXD or LSYSSC_MATRIX9,...
  INTEGER, INTENT(IN) :: cmatrixFormat
  
  ! Number of equations, the matrix should hold.
  INTEGER(PREC_MATIDX), INTENT(IN) :: NEQ
!</input>

!<inputoutput>
  ! The FE matrix to be initialised.
  TYPE(t_matrixScalar), INTENT(OUT) :: rmatrix
!</inputoutput>

!</subroutine>

    rmatrix%NEQ = NEQ
    rmatrix%NCOLS = NEQ
    rmatrix%NA = NEQ
    rmatrix%cmatrixFormat = cmatrixFormat
    
    ! Depending on the format, choose the initialisation routine.
    SELECT CASE (cmatrixFormat)
    CASE (LSYSSC_MATRIXD,LSYSSC_MATRIX1)
      ! Nothing to do, there is no structure.
      
    CASE (LSYSSC_MATRIX9)
      ! Create KCOL, KLD, KDiagonal
      CALL storage_new1D ('lsyssc_createDiagMatrix', 'KCOL', NEQ, &
          ST_INT, rmatrix%h_Kcol, ST_NEWBLOCK_ORDERED)
      CALL storage_new1D ('lsyssc_createDiagMatrix', 'KLD', NEQ+1_PREC_MATIDX, &
          ST_INT, rmatrix%h_Kld, ST_NEWBLOCK_ORDERED)
      CALL storage_new1D ('lsyssc_createDiagMatrix', 'KDiagonal', NEQ, &
          ST_INT, rmatrix%h_Kdiagonal, ST_NEWBLOCK_ORDERED)
          
    CASE (LSYSSC_MATRIX7)
      ! Create KCOL, KLD.
      CALL storage_new1D ('lsyssc_createDiagMatrix', 'KCOL', NEQ, &
          ST_INT, rmatrix%h_Kcol, ST_NEWBLOCK_ORDERED)
      CALL storage_new1D ('lsyssc_createDiagMatrix', 'KLD', NEQ+1_PREC_MATIDX, &
          ST_INT, rmatrix%h_Kld, ST_NEWBLOCK_ORDERED)
          
    CASE DEFAULT
    
      PRINT *,'lsyssc_createDiagMatrix: unsupported matrix format!'
      CALL sys_halt()
      
    END SELECT
  
  END SUBROUTINE

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_lumpMatrixScalar (rmatrixScalar,clumpType,bconvertToDiag)
  
!<description>
  ! This routine performs (mass) lumping with a scalar matrix. Lumping
  ! is in this sense a (more or less) algebraic operation to convert a given
  ! matrix into a diagonal matrix. There are different types of lumping
  ! supported:\\
  ! 'Simple lumping' simply takes the diagonal from a given matrix, forcing
  !   all off-diagonal entries to be zero.\\
  ! 'Diagonal lumping' adds all off-diagonal entries to the main diagonal and
  !   extracts the diagonal entries of the resulting matrix to form a
  !   diagonal matrix.\\
  !
  ! If bconvertToDiag=TRUE, the matrix structure will be changed into a
  ! diagonal matrix, which is the standard case. When setting this parameter
  ! to FALSE, the matrix structure stays unchanged.
!</description>

!<input>
  ! Type of lumping to perform. One of the LSYSSC_LUMP_xxxx constants.
  ! LSYSSC_LUMP_STD:  Perform standard lumping.
  ! LSYSSC_LUMP_DIAG: Perform diagonal lumping.
  INTEGER, INTENT(IN) :: clumpType
  
  ! OPTIONAL: Whether to convert the matrix structure to a diagonal matrix.
  ! TRUE:  Convert the matrix to a diagonal matrix. This is the default case
  !        if the parameter is not present.
  ! FALSE: Matrix structure stays unchanged.
  LOGICAL, INTENT(IN), OPTIONAL :: bconvertToDiag
!</input>

!<inputoutput>
  ! The matrix to be lumped. Will be overwritten by a diagonal matrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>
 
  ! local variables
  INTEGER(PREC_MATIDX) :: irow,icol,j
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kdiagonal, p_Kld
  REAL(DP), DIMENSION(:), POINTER :: p_Da
  REAL(SP), DIMENSION(:), POINTER :: p_Fa
  LOGICAL :: bconvert

  IF (rmatrixScalar%NEQ .EQ. 0) RETURN
  
  bconvert = .TRUE.
  IF (PRESENT(bconvertToDiag)) bconvert = bconvertToDiag

  ! Which type of lumping do we have to do?
  
  SELECT CASE (clumpType)
  CASE (LSYSSC_LUMP_STD)
    ! Standard lumping. Nothing to do here.
  
  CASE (LSYSSC_LUMP_DIAG)
    ! Diagonal lumping. We have to add all off-diagonal entries to the main
    ! diagonal.
    !
    ! What's the matrix format?
    SELECT CASE (rmatrixScalar%cmatrixFormat)
    
    CASE (LSYSSC_MATRIX7)
      ! Get the structure
      CALL lsyssc_getbase_Kld (rmatrixScalar,p_Kld)
      CALL lsyssc_getbase_Kcol (rmatrixScalar,p_Kcol)
      
      ! Get the data and perform the lumping
      SELECT CASE (rmatrixScalar%cdataType)
      
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rmatrixScalar,p_Da)
        
        ! Add all off-diagonals to the diagonal
        DO irow=1,rmatrixScalar%NEQ
          j = p_Kld(irow)
          DO icol = p_Kld(irow)+1,p_Kld(irow+1)-1
            p_Da(j) = p_Da(j) + p_Da(icol)
          END DO
        END DO
      
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rmatrixScalar,p_Fa)

        ! Add all off-diagonals to the diagonal
        DO irow=1,rmatrixScalar%NEQ
          j = p_Kld(irow)
          DO icol = p_Kld(irow)+1,p_Kld(irow+1)-1
            p_Fa(j) = p_Fa(j) + p_Fa(icol)
          END DO
        END DO
      
      CASE DEFAULT
        PRINT *,'lsyssc_lumpMatrixScalar: Unsupported matrix precision'
        CALL sys_halt()
      END SELECT
      
    CASE (LSYSSC_MATRIX9)
      ! Get the structure
      CALL lsyssc_getbase_Kld (rmatrixScalar,p_Kld)
      CALL lsyssc_getbase_Kdiagonal (rmatrixScalar,p_Kdiagonal)
      CALL lsyssc_getbase_Kcol (rmatrixScalar,p_Kcol)
      
      ! Get the data and perform the lumping
      SELECT CASE (rmatrixScalar%cdataType)
      
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double (rmatrixScalar,p_Da)
        
        ! Add all off-diagonals to the diagonal
        DO irow=1,rmatrixScalar%NEQ
          j = p_Kdiagonal(irow)
          DO icol = p_Kld(irow),j-1
            p_Da(j) = p_Da(j) + p_Da(icol)
          END DO
          DO icol = j+1,p_Kld(irow+1)-1
            p_Da(j) = p_Da(j) + p_Da(icol)
          END DO
        END DO
      
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single (rmatrixScalar,p_Fa)

        ! Add all off-diagonals to the diagonal
        DO irow=1,rmatrixScalar%NEQ
          j = p_Kdiagonal(irow)
          DO icol = p_Kld(irow),j-1
            p_Fa(j) = p_Fa(j) + p_Fa(icol)
          END DO
          DO icol = j+1,p_Kld(irow+1)-1
            p_Fa(j) = p_Fa(j) + p_Fa(icol)
          END DO
        END DO
      
      CASE DEFAULT
        PRINT *,'lsyssc_lumpMatrixScalar: Unsupported matrix precision'
        CALL sys_halt()
      END SELECT
      
    CASE DEFAULT
      PRINT *,'lsyssc_lumpMatrixScalar: Unsupported matrix format'
      CALL sys_halt()
    END SELECT
  
  END SELECT

  IF (bconvert) THEN
    ! Convert the given matrix into a diagonal matrix by extracting 
    ! the diagonal entries.
    CALL lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIXD)
  ELSE
    ! Clear all off-diagonal entries, so the matrix will be a diagonal
    ! matrix in the original structure. 
    CALL lsyssc_clearOffdiags (rmatrixScalar)
  END IF

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_clearOffdiags (rmatrix)
  
!<description>
  ! Sets all off-diagonal entries in the matrix rmatrix to zero.
!</description>

!<inputoutput>
  ! The matrix which is to be modified.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>

!</subroutine>

  ! At first we must take care of the matrix type.
  SELECT CASE (rmatrix%cmatrixFormat)
  CASE (LSYSSC_MATRIX9)
    CALL removeOffdiags_format9 (rmatrix)
  CASE (LSYSSC_MATRIX7)
    CALL removeOffdiags_format7 (rmatrix)
  CASE (LSYSSC_MATRIXD)
    ! Nothing to do
  CASE (LSYSSC_MATRIX1)
    CALL removeOffdiags_format1 (rmatrix) 
  CASE DEFAULT
    PRINT *,'lsyssc_clearOffdiags: Unsupported matrix format'
    CALL sys_halt()
  END SELECT
  
  CONTAINS
   
    ! ****************************************
    ! The replacement routine for format 9
    
    SUBROUTINE removeOffdiags_format9 (rmatrix)
    
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
      DO irow = 1,rmatrix%NEQ
      
        ! Get the diagonal
        ddiag = p_DA(p_Kdiagonal(irow))
      
        ! Clear the row
        p_DA(p_Kld(irow):p_Kld(irow+1)-1) = 0.0_DP
        
        ! restore the diagonal
        p_DA(p_Kdiagonal(irow)) = ddiag
      
      END DO
      
    CASE (ST_SINGLE)
      ! Get the data array
      CALL lsyssc_getbase_single(rmatrix,p_FA)
      
      ! loop through the rows
      DO irow = 1,rmatrix%NEQ
      
        ! Get the diagonal
        fdiag = p_FA(p_Kdiagonal(irow))
      
        ! Clear the row
        p_FA(p_Kld(irow):p_Kld(irow+1)-1) = 0.0_SP
        
        ! restore the diagonal
        p_FA(p_Kdiagonal(irow)) = Fdiag
      
      END DO
      
    CASE DEFAULT
      PRINT *,'removeOffdiags_format9: Unsupported matrix precision!'
      CALL sys_halt()
      
    END SELECT
    
    END SUBROUTINE

    ! ****************************************
    ! The replacement routine for format 7
    
    SUBROUTINE removeOffdiags_format7 (rmatrix)
    
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
      DO irow = 1,rmatrix%NEQ
      
        ! Clear the row except for the diagonal
        p_DA(p_Kld(irow)+1:p_Kld(irow+1)-1) = 0.0_DP
      
      END DO
      
    CASE (ST_SINGLE)
      ! Get the data array
      CALL lsyssc_getbase_single(rmatrix,p_FA)
      
      ! loop through the rows
      DO irow = 1,rmatrix%NEQ
      
        ! Clear the row except for the diagonal
        p_FA(p_Kld(irow)+1:p_Kld(irow+1)-1) = 0.0_SP
      
      END DO
      
    CASE DEFAULT
      PRINT *,'removeOffdiags_format7: Unsupported matrix precision!'
      CALL sys_halt()
    END SELECT

    END SUBROUTINE
  
    ! ****************************************
    ! The replacement routine for format 1
    
    SUBROUTINE removeOffdiags_format1 (rmatrix)
    
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
    
    ! local variables
    INTEGER(PREC_MATIDX) :: irow,icol
    REAL(DP), DIMENSION(:), POINTER :: p_DA
    REAL(SP), DIMENSION(:), POINTER :: p_FA
    REAL(DP) :: ddata
    REAL(SP) :: fdata

    ! Take care of the format of the entries
    SELECT CASE (rmatrix%cdataType)
    CASE (ST_DOUBLE)
      ! Get the data array
      CALL lsyssc_getbase_double(rmatrix,p_DA)
      
      ! loop through the rows ald colums, set off-diagonal entries to zero.
      DO icol = 1,rmatrix%NCOLS
        IF (icol .LE. rmatrix%NEQ) THEN
          ! Remember the diagonal
          ddata = p_Da((icol-1)*rmatrix%NEQ+icol)
          
          ! Fill the column with zero
          p_Da ((icol-1)*rmatrix%NEQ+1 : icol*rmatrix%NEQ) = 0.0_DP
          
          ! Restore the diagonal
          p_Da((icol-1)*rmatrix%NEQ+icol) = ddata
        ELSE
          ! Fill the column with zero
          p_Da ((icol-1)*rmatrix%NEQ+1 : icol*rmatrix%NEQ) = 0.0_DP
        END IF
      END DO
      
    CASE (ST_SINGLE)
      ! Get the data array
      CALL lsyssc_getbase_single(rmatrix,p_FA)
      
      ! loop through the rows ald colums, set off-diagonal entries to zero.
      DO icol = 1,rmatrix%NCOLS
        IF (icol .LE. rmatrix%NEQ) THEN
          ! Remember the diagonal
          fdata = p_Fa((icol-1)*rmatrix%NEQ+icol)
          
          ! Fill the column with zero
          p_Fa ((icol-1)*rmatrix%NEQ+1 : icol*rmatrix%NEQ) = 0.0_DP
          
          ! Restore the diagonal
          p_Fa((icol-1)*rmatrix%NEQ+icol) = fdata
        ELSE
          ! Fill the column with zero
          p_Fa ((icol-1)*rmatrix%NEQ+1 : icol*rmatrix%NEQ) = 0.0_DP
        END IF
      END DO
      
    CASE DEFAULT
      PRINT *,'removeOffdiags_format7: Unsupported matrix precision!'
      CALL sys_halt()
    END SELECT
    
    END SUBROUTINE

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_scaleMatrix (rmatrix,c)
  
!<description>
  ! Scales a matrix rmatrix: rmatrix = c * rmatrix
!</description>
  
!<inputoutput>
  ! Source and destination matrix
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
!</inputoutput>

!<input>
  ! Multiplication factor
  REAL(DP), INTENT(IN) :: c
!</input>
  
!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata
  
  ! Take care of the data type!
  SELECT CASE (rmatrix%cdataType)
  CASE (ST_DOUBLE)
    ! Get the pointer and scale the whole data array.
    CALL lsyssc_getbase_double(rmatrix,p_Ddata)
    CALL lalg_scaleVectorDble (p_Ddata,c)  

  CASE (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    CALL lsyssc_getbase_single(rmatrix,p_Fdata)
    CALL lalg_scaleVectorSngl (p_Fdata,REAL(c,SP))  

  CASE DEFAULT
    PRINT *,'lsyssc_scaleMatrix: Unsupported data type!'
    CALL sys_halt()
  END SELECT
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_multMatMat (rmatrixA,rmatrixB,rmatrixC,bmemory&
      &,bsymb,bnumb,bisExactStructure)

!<description>
    ! Computes the matrix-matrix-product 
    !   rmatrixC := rmatrixA * rmatrixB
    !
    ! All matrix formats and most combinations of matrix formats are
    ! supported for matrices A and B. The resulting matrix C=A*B is
    ! generated in the "encompassing" format. That is, if A and B are
    ! both diagonal matrices, then C is also a diagonal matrix. If
    ! for example, A is a full matrix and B is a diagonal matrix (or
    ! vice versa) then C is also a full matrix so that all matrix
    ! entries can be stored. The only combination which is not
    ! supported is full and sparse matrices A and B. In general,
    ! sparse matrices are used for large data. Hence, it does not
    ! make sense to multiply a sparse matrix and a full one and store
    ! the result in a full matrix which is quite likely to run out of
    ! memory. Both CSR formats 7 and 9 can be combined. Note that if
    ! at least one matrix A and/or B is stored in format 9, then the
    ! resulting matrix C will also be stored in format 9.
!</description>

!<input>
    ! source matrices
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrixA,rmatrixB

    ! BMEMORY = FALSE: Do not allocate required memory for C=A*B.
    ! BMEMORY = TRUE:  Generate all required structures for C=A*B 
    LOGICAL :: bmemory

    ! Compute symbolic matrix-matrix-product
    ! BSYMB = FALSE: Do not generate the required matrix structures.
    !                This may be useful, if the same matrices A and B
    !                need to be multiplied several times, but the
    !                sparsity patterns do not change
    ! BSYMN = TRUE:  Generate all required matrix structures for C=A*B
    LOGICAL :: bsymb

    ! Compute numerical matrix-matrix-product
    ! BNUMB = FALSE: Do not perform the numerical multiplication
    ! BNUMB = TRUE:  Perform numerical multiplication
    LOGICAL :: bnumb

    ! OPTIONAL: Indicates whether the resulting matrix C has the
    ! required symbolic structure or is a superset of the symbolic
    ! matrix-matrix product. In some cases, this may allow for much
    ! more efficient implementation.
    LOGICAL, OPTIONAL :: bisExactStructure
!</input>

!<inputoutput>
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixC
!</inputoutput>

!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: DaA,DaB,DaC,Daux
    REAL(SP), DIMENSION(:), POINTER :: FaA,FaB,FaC,Faux
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: KcolA,KcolB,KcolC,Kaux
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: KldA,KldB,KldC,KdiagonalC
    INTEGER :: h_Kaux,h_Daux,h_Faux
    LOGICAL :: bfast

    h_Kaux=ST_NOHANDLE
    h_Daux=ST_NOHANDLE
    h_Faux=ST_NOHANDLE

    ! Check if fast implementation can be used
    bfast=.FALSE.
    IF (PRESENT(bisExactStructure)) bfast = bisExactStructure
    bfast = bfast .OR. (bmemory .AND. bsymb)

    ! Check if both matrices are compatible
    IF (rmatrixA%NCOLS /= rmatrixB%NEQ) THEN
      PRINT *, 'lsyssc_multMatMat: number of columns of matrix A is not&
          & compatible with number of rows of matrix B'
      CALL sys_halt()
    END IF

    ! Check if both matrices have the same sorting
    IF (rmatrixA%isortStrategy /= rmatrixB%isortStrategy) THEN
      PRINT *, 'lsyssc_multMatMat: incompatible sorting strategies'
      CALL sys_halt()
    END IF

    ! Release matrix if required and set common variables
    IF (bmemory) THEN
      CALL lsyssc_releaseMatrix(rmatrixC)
      IF (rmatrixA%cdataType == ST_DOUBLE .OR. rmatrixB&
          &%cdataType == ST_DOUBLE) THEN
        rmatrixC%cdataType = ST_DOUBLE
      ELSE
        rmatrixC%cdataType = ST_SINGLE
      END IF
    END IF
    
    ! Set sorting strategy for matrix C
    rmatrixC%isortStrategy = rmatrixA%isortStrategy
    rmatrixC%h_IsortPermutation = rmatrixA%h_IsortPermutation
    
    ! Perform matrix-matrix multiplication
    SELECT CASE(rmatrixA%cmatrixFormat)

    CASE (LSYSSC_MATRIX1)   ! A is full matrix ----------------------
      
      SELECT CASE(rmatrixB%cmatrixFormat)
        
      CASE (LSYSSC_MATRIX1) ! B is full matrix - - - - - - - - - - -

        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = LSYSSC_MATRIX1
          rmatrixC%NA = rmatrixA%NEQ*rmatrixB%NCOLS
          CALL storage_new('lsyssc_multMatMat','h_Da',rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= LSYSSC_MATRIX1) THEN
          PRINT *, 'lsyssc_multMatMat: destination matrix has incompati&
              &ble format'
          CALL sys_halt()
        END IF

        ! symbolical multiplication?
        IF (bsymb) THEN
          rmatrixC%NEQ = rmatrixA%NEQ
          rmatrixC%NCOLS = rmatrixB%NCOLS
        END IF
        
        ! numerical multiplication?
        IF (bnumb) THEN
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double          
          SELECT CASE(rmatrixA%cdataType)

          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat1mat1mul_doubledouble(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,rmatrixB%NCOLS,DaA,DaB,DaC)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat1mat1mul_doublesingle(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,rmatrixB%NCOLS,DaA,FaB,DaC)
              
            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat1mat1mul_singledouble(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,rmatrixB%NCOLS,FaA,DaB,DaC)
              
            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL do_mat1mat1mul_singlesingle(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,rmatrixB%NCOLS,FaA,FaB,FaC)
              
            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE DEFAULT
            PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - - 

        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = LSYSSC_MATRIX1
          rmatrixC%NA = rmatrixA%NA
          CALL storage_new('lsyssc_multMatMat','h_Da',rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= LSYSSC_MATRIX1) THEN
          PRINT *, 'lsyssc_multMatMat: destination matrix has incompati&
              &ble format'
          CALL sys_halt()
        END IF

        ! symbolical multiplication?
        IF (bsymb) THEN
          rmatrixC%NEQ = rmatrixA%NEQ
          rmatrixC%NCOLS = rmatrixA%NCOLS
        END IF

        ! numerical multiplication?
        IF (bnumb) THEN
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double          
          SELECT CASE(rmatrixA%cdataType)

          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat1matDmul_doubledouble(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,DaA,DaB,DaC)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat1matDmul_doublesingle(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,DaA,FaB,DaC)

            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)

            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat1matDmul_singledouble(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,FaA,DaB,DaC)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL do_mat1matDmul_singlesingle(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,FaA,FaB,FaC)
              
            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT
        
          CASE DEFAULT
            PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE DEFAULT
        PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
        CALL sys_halt()
      END SELECT

      !--------------------------------------------------------------
    
    CASE (LSYSSC_MATRIXD)   ! A is diagonal matrix ------------------
      
      SELECT CASE(rmatrixB%cmatrixFormat)

      CASE (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - -
       
        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = LSYSSC_MATRIXD
          rmatrixC%NA = rmatrixA%NA
          CALL storage_new('lsyssc_multMatMat','h_Da',rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= LSYSSC_MATRIXD) THEN
          PRINT *, 'lsyssc_multMatMat: destination matrix has incompati&
              &ble format'
          CALL sys_halt()
        END IF

        ! symbolic multiplication?
        IF (bsymb) THEN
          rmatrixC%NEQ = rmatrixA%NEQ
          rmatrixC%NCOLS = rmatrixA%NCOLS
        END IF
        
        ! numerical multiplication?
        IF (bnumb) THEN
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double 
          SELECT CASE(rmatrixA%cdataType)

          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              DaC=DaA*DaB

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              DaC=DaA*FaB

            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE (ST_SINGLE)

            SELECT CASE(rmatrixC%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              DaC=FaA*DaB

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              FaC=FaA*FaB

            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE DEFAULT
            PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE (LSYSSC_MATRIX1) ! B is full matrix - - - - - - - - - - - 

        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = LSYSSC_MATRIX1
          rmatrixC%NA = rmatrixB%NA
          CALL storage_new('lsyssc_multMatMat','h_Da',rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= LSYSSC_MATRIX1) THEN
          PRINT *, 'lsyssc_multMatMat: destination matrix has incompati&
              &ble format'
          CALL sys_halt()
        END IF

        ! symbolical multiplication?
        IF (bsymb) THEN
          rmatrixC%NEQ = rmatrixB%NEQ
          rmatrixC%NCOLS = rmatrixB%NCOLS
        END IF

        ! numerical multiplication?
        IF (bnumb) THEN
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double          
          SELECT CASE(rmatrixA%cdataType)

          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_matDmat1mul_doubledouble(rmatrixB%NEQ,rmatrixB&
                  &%NCOLS,DaA,DaB,DaC)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_matDmat1mul_doublesingle(rmatrixB%NEQ,rmatrixB&
                  &%NCOLS,DaA,FaB,DaC)

            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)

            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_matDmat1mul_singledouble(rmatrixB%NEQ,rmatrixB&
                  &%NCOLS,FaA,DaB,DaC)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL do_matDmat1mul_singlesingle(rmatrixB%NEQ,rmatrixB&
                  &%NCOLS,FaA,FaB,FaC)
              
            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT
        
          CASE DEFAULT
            PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! B is CSR matrix - - - - 

        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = rmatrixB%cmatrixFormat
          rmatrixC%NA = rmatrixB%NA
          CALL storage_new('lsyssc_multMatMat','h_Da',rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= rmatrixB%cmatrixFormat) THEN
          PRINT *, 'lsyssc_multMatMat: destination matrix has incompati&
              &ble format'
          CALL sys_halt()
        END IF

        ! symbolical multiplication?
        IF (bsymb) THEN
          rmatrixC%NEQ = rmatrixB%NEQ
          rmatrixC%NCOLS = rmatrixB%NCOLS
          rmatrixC%h_Kld = rmatrixB%h_Kld
          rmatrixC%h_Kcol = rmatrixB%h_Kcol
          rmatrixC%h_Kdiagonal = rmatrixB%h_Kdiagonal
          rmatrixC%imatrixSpec = IOR(rmatrixC%imatrixSpec,&
              &LSYSSC_MSPEC_STRUCTUREISCOPY)
        END IF

        ! numerical multiplication?
        IF (bnumb) THEN
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          SELECT CASE(rmatrixA%cdataType)

          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              IF (bfast) THEN
                CALL do_matDmat79mul_doubledouble(DaA,DaB,KldB,KcolB&
                    &,rmatrixB%NEQ,DaC)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_matDmat79mul_doubledouble(DaA,DaB,KldB,KcolB&
                    &,rmatrixB%NEQ,DaC,KldC,KcolC)
              END IF

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kld(rmatrixB,KcolB)
              IF (bfast) THEN
                CALL do_matDmat79mul_doublesingle(DaA,FaB,KldB,KcolB&
                    &,rmatrixB%NEQ,DaC)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_matDmat79mul_doublesingle(DaA,FaB,KldB,KcolB&
                    &,rmatrixB%NEQ,DaC,KldC,KcolC)
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)
            
            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              IF (bfast) THEN
                CALL do_matDmat79mul_singledouble(FaA,DaB,KldB,KcolB&
                    &,rmatrixB%NEQ,DaC)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_matDmat79mul_singledouble(FaA,DaB,KldB,KcolB&
                    &,rmatrixB%NEQ,DaC,KldC,KcolC)
              END IF

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              IF (bfast) THEN
                CALL do_matDmat79mul_singlesingle(FaA,FaB,KldB,KcolB&
                    &,rmatrixB%NEQ,FaC)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_matDmat79mul_singlesingle(FaA,FaB,KldB,KcolB&
                    &,rmatrixB%NEQ,FaC,KldC,KcolC)
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE DEFAULT
            PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE DEFAULT
        PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
        CALL sys_halt()
      END SELECT
      
      !--------------------------------------------------------------

    CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! A is CSR matrix ----------
      
      SELECT CASE(rmatrixB%cmatrixFormat)
        
      CASE (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - -
        
        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = rmatrixA%cmatrixFormat
          rmatrixC%NA = rmatrixA%NA
          CALL storage_new('lsyssc_multMatMat','h_Da',rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= rmatrixA%cmatrixFormat) THEN
          PRINT *, 'lsyssc_multMatMat: destination matrix has incompati&
              &ble format'
          CALL sys_halt()
        END IF

        ! symbolical multiplication?
        IF (bsymb) THEN
          rmatrixC%NEQ = rmatrixA%NEQ
          rmatrixC%NCOLS = rmatrixA%NCOLS
          rmatrixC%h_Kld = rmatrixA%h_Kld
          rmatrixC%h_Kcol = rmatrixA%h_Kcol
          rmatrixC%h_Kdiagonal = rmatrixA%h_Kdiagonal
          rmatrixC%imatrixSpec = IOR(rmatrixC%imatrixSpec,&
              &LSYSSC_MSPEC_STRUCTUREISCOPY)
        END IF

        ! numerical multiplication?
        IF (bnumb) THEN
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          SELECT CASE(rmatrixA%cdataType)

          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              IF (bfast) THEN
                CALL do_mat79matDmul_doubledouble(DaA,KldA,KcolA&
                    &,rmatrixA%NEQ,DaB,DaC)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_mat79matDmul_doubledouble(DaA,KldA,KcolA&
                    &,rmatrixA%NEQ,DaB,DaC,KldC,KcolC)
              END IF

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              IF (bfast) THEN
                CALL do_mat79matDmul_doublesingle(DaA,KldA,KcolA&
                    &,rmatrixA%NEQ,FaB,DaC)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_mat79matDmul_doublesingle(DaA,KldA,KcolA&
                    &,rmatrixA%NEQ,FaB,DaC,KldC,KcolC)
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)
            
            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              IF (bfast) THEN
                CALL do_mat79matDmul_singledouble(FaA,KldA,KcolA&
                    &,rmatrixA%NEQ,DaB,DaC)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_mat79matDmul_singledouble(FaA,KldA,KcolA&
                    &,rmatrixA%NEQ,DaB,DaC,KldC,KcolC)
              END IF

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              IF (bfast) THEN
                CALL do_mat79matDmul_singlesingle(FaA,KldA,KcolA&
                    &,rmatrixA%NEQ,FaB,FaC)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_mat79matDmul_singlesingle(FaA,KldA,KcolA&
                    &,rmatrixA%NEQ,FaB,FaC,KldC,KcolC)
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE DEFAULT
            PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! B is CSR matrix - - - - 
        
        ! Set pointers
        CALL lsyssc_getbase_Kld(rmatrixA,KldA)
        CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
        CALL lsyssc_getbase_Kld(rmatrixB,KldB)
        CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)

        ! memory allocation?
        IF (bmemory) THEN
          ! Duplicate structure of matrix A or B depending on which
          ! matrix is given in the encompassing matrix format
          IF (rmatrixA%cmatrixFormat >= rmatrixB%cmatrixFormat) THEN
            CALL lsyssc_duplicateMatrix(rmatrixA,rmatrixC&
                &,LSYSSC_DUP_EMPTY,LSYSSC_DUP_REMOVE)
          ELSE
            CALL lsyssc_duplicateMatrix(rmatrixB,rmatrixC&
                &,LSYSSC_DUP_EMPTY,LSYSSC_DUP_REMOVE)
          END IF
          
          ! Set auxiliary pointers
          CALL storage_new('lsyssc_multMatMat','Kaux',rmatrixB%NCOLS,&
              &ST_INT,h_Kaux,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int(h_Kaux,Kaux)
          
          ! Compute number of nonzero matrix entries: NA
          rmatrixC%NA = do_mat79mat79mul_computeNA(KldA,KcolA&
              &,rmatrixA%NEQ,KldB,KcolB,Kaux)
          CALL storage_new('lsyssc_multMatMat','h_Da',rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
          CALL storage_realloc('lsyssc_multMatMat',rmatrixC%NA&
              &,rmatrixC%h_Kcol,ST_NEWBLOCK_NOINIT,.FALSE.)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= MAX(rmatrixA%cmatrixFormat,&
            &rmatrixB%cmatrixFormat)) THEN
          PRINT *, 'lsyssc_multMatMat: destination matrix has incompati&
              &ble format'
          CALL sys_halt()
        END IF

        ! symbolical multiplication?
        IF (bsymb) THEN

          ! Set pointers
          CALL lsyssc_getbase_Kld(rmatrixC,KldC)
          CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
          IF (h_Kaux == ST_NOHANDLE) THEN
            CALL storage_new('lsyssc_multMatMat','Kaux',MAX(rmatrixA&
                &%NEQ,MAX(rmatrixA%NCOLS,rmatrixB%NCOLS)),ST_INT&
                &,h_Kaux,ST_NEWBLOCK_NOINIT)
          ELSE
            CALL storage_realloc('lsyssc_multMatMat',MAX(rmatrixA%NEQ&
                &,MAX(rmatrixA%NCOLS,rmatrixB%NCOLS)),h_Kaux&
                &,ST_NEWBLOCK_NOINIT,.FALSE.)
          END IF
          CALL storage_getbase_int(h_Kaux,Kaux)

          IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX9) THEN
            CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
            CALL do_mat79mat79mul_symb(rmatrixA%NEQ,rmatrixA%NCOLS&
                &,rmatrixB%NCOLS,KldA,KcolA,KldB,KcolB,KldC,KcolC&
                &,Kaux,KdiagonalC)
          ELSE
            CALL do_mat79mat79mul_symb(rmatrixA%NEQ,rmatrixA%NCOLS&
                &,rmatrixB%NCOLS,KldA,KcolA,KldB,KcolB,KldC,KcolC&
                &,Kaux)
          END IF
        END IF

        ! numerical multiplication?
        IF (bnumb) THEN

          ! Set pointers
          CALL lsyssc_getbase_Kld(rmatrixC,KldC)
          CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
          IF (h_Daux == ST_NOHANDLE) THEN
            CALL storage_new('lsyssc_multMatMat','Daux',MAX(rmatrixA&
                &%NEQ,MAX(rmatrixA%NCOLS,rmatrixB%NCOLS)),ST_DOUBLE&
                &,h_Daux,ST_NEWBLOCK_NOINIT)
          ELSE
            CALL storage_realloc('lsyssc_multMatMat',MAX(rmatrixA%NEQ&
                &,MAX(rmatrixA%NCOLS,rmatrixB%NCOLS)),h_Daux&
                &,ST_NEWBLOCK_NOINIT,.FALSE.)
          END IF
          CALL storage_getbase_double(h_Daux,Daux)

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double   
          SELECT CASE(rmatrixA%cdataType)
            
          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE) 
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat79mat79mul_numb_dbledble(rmatrixA%NEQ&
                  &,rmatrixA%NCOLS,rmatrixB%NCOLS,KldA,KcolA,DaA,KldB&
                  &,KcolB,DaB,KldC,KcolC,DaC,Daux)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat79mat79mul_numb_dblesngl(rmatrixA%NEQ&
                  &,rmatrixA%NCOLS,rmatrixB%NCOLS,KldA,KcolA,DaA,KldB&
                  &,KcolB,FaB,KldC,KcolC,DaC,Daux)
              
            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)
            
            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE) 
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat79mat79mul_numb_sngldble(rmatrixA%NEQ&
                  &,rmatrixA%NCOLS,rmatrixB%NCOLS,KldA,KcolA,FaA,KldB&
                  &,KcolB,DaB,KldC,KcolC,DaC,Daux)

            CASE (ST_SINGLE)
              ! We need an additional vector FAUX
              CALL storage_new('lsyssc_multMatMat','Faux',MAX(rmatrixA&
                  &%NEQ,MAX(rmatrixA%NCOLS,rmatrixB%NCOLS)),ST_SINGLE&
                  &,h_Faux,ST_NEWBLOCK_NOINIT)
              CALL storage_getbase_single(h_Faux,Faux)
              
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL do_mat79mat79mul_numb_snglsngl(rmatrixA%NEQ&
                  &,rmatrixA%NCOLS,rmatrixB%NCOLS,KldA,KcolA,FaA,KldB&
                  &,KcolB,FaB,KldC,KcolC,FaC,Faux)

            CASE DEFAULT
              PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE DEFAULT
            PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF
        
      CASE DEFAULT
        PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
        CALL sys_halt()
      END SELECT
            
      !--------------------------------------------------------------
      
    CASE DEFAULT
      PRINT *, 'lsyssc_multMatMat: Unsupported data type!'
      CALL sys_halt()
    END SELECT

    ! Clear auxiliary vectors
    IF (h_Kaux /= ST_NOHANDLE) CALL storage_free(h_Kaux)
    IF (h_Daux /= ST_NOHANDLE) CALL storage_free(h_Daux)
    IF (h_Faux /= ST_NOHANDLE) CALL storage_free(h_Faux)

  CONTAINS

    ! Here, the real MM multiplication routines follow.

    !**************************************************************
    ! Format 1 multiplication
    ! double precision matrix A (format 1)
    ! double precision matrix B (format 1)
    ! double precision matrix C (format 1)

    SUBROUTINE do_mat1mat1mul_doubledouble(n,m,k,Da1,Da2,Da3)
      
      ! Remark: MATRIX1 is stored row-wise. The intrinsic function
      !         MATMUL requires the matrix to be stored columnwise
      !         Hence, compute C = A*B = (B'*A')'
      
      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m,k
      REAL(DP), DIMENSION(m,n), INTENT(IN)    :: Da1
      REAL(DP), DIMENSION(k,m), INTENT(IN)    :: Da2
      REAL(DP), DIMENSION(k,n), INTENT(INOUT) :: Da3

      ! Remark: The tuned BLAS3 routine performs much better than the
      ! intrinsic MATMUL Fortran90 routine. Hence, the BLAS routine
      ! is used whenever possible, that is, when all matrices have
      ! the same precision.
      CALL DGEMM('N','N',k,n,m,1.0_DP,Da2,k,Da1,m,0.0_DP,Da3,k)
    END SUBROUTINE do_mat1mat1mul_doubledouble

    !**************************************************************
    ! Format 1 multiplication
    ! double precision matrix A (format 1)
    ! single precision matrix B (format 1)
    ! double precision matrix C (format 1)

    SUBROUTINE do_mat1mat1mul_doublesingle(n,m,k,Da1,Fa2,Da3)
      
      ! Remark: MATRIX1 is stored row-wise. The intrinsic function
      !         MATMUL requires the matrix to be stored columnwise
      !         Hence, compute C = A*B = (B'*A')'
      
      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m,k
      REAL(DP), DIMENSION(m,n), INTENT(IN)    :: Da1
      REAL(SP), DIMENSION(k,m), INTENT(IN)    :: Fa2
      REAL(DP), DIMENSION(k,n), INTENT(INOUT) :: Da3

      CALL DGEMM('N','N',k,n,m,1.0_DP,REAL(Fa2,DP),k,Da1,m,0.0_DP,Da3,k)
!!$      Da3=MATMUL(Fa2,Da1)
    END SUBROUTINE do_mat1mat1mul_doublesingle

    !**************************************************************
    ! Format 1 multiplication
    ! single precision matrix A (format 1)
    ! double precision matrix B (format 1)
    ! double precision matrix C (format 1)

    SUBROUTINE do_mat1mat1mul_singledouble(n,m,k,Fa1,Da2,Da3)
      
      ! Remark: MATRIX1 is stored row-wise. The intrinsic function
      !         MATMUL requires the matrix to be stored columnwise
      !         Hence, compute C = A*B = (B'*A')'
      
      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m,k
      REAL(SP), DIMENSION(m,n), INTENT(IN)    :: Fa1
      REAL(DP), DIMENSION(k,m), INTENT(IN)    :: Da2
      REAL(DP), DIMENSION(k,n), INTENT(INOUT) :: Da3
 
      CALL DGEMM('N','N',k,n,m,1.0_DP,Da2,k,REAL(Fa1,DP),m,0.0_DP,Da3,k)
!!$      Da3=MATMUL(Da2,Fa1)
    END SUBROUTINE do_mat1mat1mul_singledouble

    !**************************************************************
    ! Format 1 multiplication
    ! single precision matrix A (format 1)
    ! single precision matrix B (format 1)
    ! single precision matrix C (format 1)

    SUBROUTINE do_mat1mat1mul_singlesingle(n,m,k,Fa1,Fa2,Fa3)
      
      ! Remark: MATRIX1 is stored row-wise. The intrinsic function
      !         MATMUL requires the matrix to be stored columnwise
      !         Hence, compute C = A*B = (B'*A')'
      
      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m,k
      REAL(SP), DIMENSION(m,n), INTENT(IN)    :: Fa1
      REAL(SP), DIMENSION(k,m), INTENT(IN)    :: Fa2
      REAL(SP), DIMENSION(k,n), INTENT(INOUT) :: Fa3
      
      ! Remark: The tuned BLAS3 routine performs much better than the
      ! intrinsic MATMUL Fortran90 routine. Hence, the BLAS routine
      ! is used whenever possible, that is, when all matrices have
      ! the same precision.
      CALL SGEMM('N','N',k,n,m,1E0,Fa2,k,Fa1,m,0E0,Fa3,k)
    END SUBROUTINE do_mat1mat1mul_singlesingle

    !**************************************************************
    ! Format D-1 multiplication
    ! double precision matrix A (format D)
    ! double precision matrix B (format 1)
    ! double precision matrix C (format 1)

    SUBROUTINE do_matDmat1mul_doubledouble(n,m,Da1,Da2,Da3)
      
      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m
      REAL(DP), DIMENSION(n), INTENT(IN)    :: Da1
      REAL(DP), DIMENSION(m,n), INTENT(IN)    :: Da2
      REAL(DP), DIMENSION(m,n), INTENT(INOUT) :: Da3

      INTEGER :: i
      
!$omp parallel do&
!$omp&default(shared)&
!$omp&private(i)
      DO i=1,n
        Da3(:,i)=Da1(i)*Da2(:,i)
      END DO
!$omp end parallel do
    END SUBROUTINE do_matDmat1mul_doubledouble

    !**************************************************************
    ! Format D-1 multiplication
    ! single precision matrix A (format D)
    ! double precision matrix B (format 1)
    ! double precision matrix C (format 1)

    SUBROUTINE do_matDmat1mul_singledouble(n,m,Fa1,Da2,Da3)
      
      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m
      REAL(SP), DIMENSION(n), INTENT(IN)    :: Fa1
      REAL(DP), DIMENSION(m,n), INTENT(IN)    :: Da2
      REAL(DP), DIMENSION(m,n), INTENT(INOUT) :: Da3

      INTEGER :: i
      
!$omp parallel do&
!$omp&default(shared)&
!$omp&private(i)
      DO i=1,n
        Da3(:,i)=Fa1(i)*Da2(:,i)
      END DO
!$omp end parallel do
    END SUBROUTINE do_matDmat1mul_singledouble

    !**************************************************************
    ! Format D-1 multiplication
    ! double precision matrix A (format D)
    ! single precision matrix B (format 1)
    ! double precision matrix C (format 1)

    SUBROUTINE do_matDmat1mul_doublesingle(n,m,Da1,Fa2,Da3)
      
      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m
      REAL(DP), DIMENSION(n), INTENT(IN)    :: Da1
      REAL(SP), DIMENSION(m,n), INTENT(IN)    :: Fa2
      REAL(DP), DIMENSION(m,n), INTENT(INOUT) :: Da3

      INTEGER :: i

!$omp parallel do&
!$omp&default(shared)&
!$omp&private(i)
      DO i=1,n
        Da3(:,i)=Da1(i)*Fa2(:,i)
      END DO
!$omp  end parallel do
    END SUBROUTINE do_matDmat1mul_doublesingle

    !**************************************************************
    ! Format D-1 multiplication
    ! single precision matrix A (format D)
    ! single precision matrix B (format 1)
    ! single precision matrix C (format 1)

    SUBROUTINE do_matDmat1mul_singlesingle(n,m,Fa1,Fa2,Fa3)
      
      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m
      REAL(SP), DIMENSION(n), INTENT(IN)    :: Fa1
      REAL(SP), DIMENSION(m,n), INTENT(IN)    :: Fa2
      REAL(SP), DIMENSION(m,n), INTENT(INOUT) :: Fa3

      INTEGER :: i

!$omp parallel do&
!$omp&default(shared)&
!$omp&private(i)
      DO i=1,n
        Fa3(:,i)=Fa1(i)*Fa2(:,i)
      END DO
!$omp end parallel do
    END SUBROUTINE do_matDmat1mul_singlesingle

    !**************************************************************
    ! Format 1-D multiplication
    ! double precision matrix A (format 1)
    ! double precision matrix B (format D)
    ! double precision matrix C (format 1)

    SUBROUTINE do_mat1matDmul_doubledouble(n,m,Da1,Da2,Da3)
      
      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m
      REAL(DP), DIMENSION(m,n), INTENT(IN)    :: Da1
      REAL(DP), DIMENSION(m),   INTENT(IN)    :: Da2
      REAL(DP), DIMENSION(m,n), INTENT(INOUT) :: Da3

      INTEGER :: i

!$omp parallel do&
!$omp&default(shared)&
!$omp&private(i)
      DO i=1,m
        Da3(i,:)=Da1(i,:)*Da2(i)
      END DO
!$omp end parallel do
    END SUBROUTINE do_mat1matDmul_doubledouble

    !**************************************************************
    ! Format 1-D multiplication
    ! single precision matrix A (format 1)
    ! double precision matrix B (format D)
    ! double precision matrix C (format 1)

    SUBROUTINE do_mat1matDmul_singledouble(n,m,Fa1,Da2,Da3)
      
      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m
      REAL(SP), DIMENSION(m,n), INTENT(IN)    :: Fa1
      REAL(DP), DIMENSION(m),   INTENT(IN)    :: Da2
      REAL(DP), DIMENSION(m,n), INTENT(INOUT) :: Da3

      INTEGER :: i

!$omp parallel do&
!$omp&default(shared)&
!$omp&private(i)
      DO i=1,m
        Da3(i,:)=Fa1(i,:)*Da2(i)
      END DO
!$omp end parallel do
    END SUBROUTINE do_mat1matDmul_singledouble

    !**************************************************************
    ! Format 1-D multiplication
    ! double precision matrix A (format 1)
    ! single precision matrix B (format D)
    ! double precision matrix C (format 1)

    SUBROUTINE do_mat1matDmul_doublesingle(n,m,Da1,Fa2,Da3)
      
      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m
      REAL(DP), DIMENSION(m,n), INTENT(IN)    :: Da1
      REAL(SP), DIMENSION(m),   INTENT(IN)    :: Fa2
      REAL(DP), DIMENSION(m,n), INTENT(INOUT) :: Da3

      INTEGER :: i

!$omp parallel do&
!$omp&default(shared)&
!$omp&private(i)  
      DO i=1,m
        Da3(i,:)=Da1(i,:)*Fa2(i)
      END DO
!$omp end parallel do
    END SUBROUTINE do_mat1matDmul_doublesingle

    !**************************************************************
    ! Format 1-D multiplication
    ! single precision matrix A (format 1)
    ! single precision matrix B (format D)
    ! single precision matrix C (format 1)

    SUBROUTINE do_mat1matDmul_singlesingle(n,m,Fa1,Fa2,Fa3)
      
      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m
      REAL(SP), DIMENSION(m,n), INTENT(IN)    :: Fa1
      REAL(SP), DIMENSION(m),   INTENT(IN)    :: Fa2
      REAL(SP), DIMENSION(m,n), INTENT(INOUT) :: Fa3

      INTEGER :: i

!$omp parallel do&
!$omp&default(shared)&
!$omp&private(i)
      DO i=1,m
        Fa3(i,:)=Fa1(i,:)*Fa2(i)
      END DO
!$omp end parallel do
    END SUBROUTINE do_mat1matDmul_singlesingle

    !**************************************************************
    ! Format D-7/9 multiplication
    ! double precision matrix A (format D)
    ! double precision matrix B (format 7 or format 9)
    ! double precision matrix C (format 7 or format 9)

    SUBROUTINE do_matDmat79mul_doubledouble(Da1,Da2,Kld2,Kcol2&
        &,neq,Da3,Kld3,Kcol3)

      REAL(DP), DIMENSION(:), INTENT(IN)    :: Da1
      REAL(DP), DIMENSION(:), INTENT(IN)    :: Da2
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: Da3
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld2
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol2
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kld3
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kcol3
      INTEGER(PREC_DOFIDX), INTENT(IN) :: neq

      INTEGER :: ieq,ild2,ild3,ildend3
      
      IF (PRESENT(Kld3) .AND. PRESENT(Kcol3)) THEN
        
        DO ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          DO ild2=Kld2(ieq),Kld2(ieq+1)-1
            DO ild3=ild3,ildend3
              IF (Kcol3(ild3) == Kcol2(ild2)) EXIT
              Da3(ild3)=0
            END DO
            Da3(ild3)=Da1(ieq)*Da2(ild2); ild3=ild3+1
          END DO
          Da3(ild3:ildend3)=0
        END DO

      ELSE

        DO ieq=1,neq
          DO ild2=Kld2(ieq),Kld2(ieq+1)-1
            Da3(ild2)=Da1(ieq)*Da2(ild2)
          END DO
        END DO

      END IF
    END SUBROUTINE do_matDmat79mul_doubledouble

    !**************************************************************
    ! Format D-7/9 multiplication
    ! single precision matrix A (format D)
    ! double precision matrix B (format 7 or format 9)
    ! double precision matrix C (format 7 or format 9)

    SUBROUTINE do_matDmat79mul_singledouble(Fa1,Da2,Kld2,Kcol2&
        &,neq,Da3,Kld3,Kcol3)

      REAL(SP), DIMENSION(:), INTENT(IN)    :: Fa1
      REAL(DP), DIMENSION(:), INTENT(IN)    :: Da2
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: Da3
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld2
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol2
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kld3
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kcol3
      INTEGER(PREC_DOFIDX), INTENT(IN) :: neq

      INTEGER :: ieq,ild2,ild3,ildend3
      
      IF (PRESENT(Kld3) .AND. PRESENT(Kcol3)) THEN
        
        DO ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          DO ild2=Kld2(ieq),Kld2(ieq+1)-1
            DO ild3=ild3,ildend3
              IF (Kcol3(ild3) == Kcol2(ild2)) EXIT
              Da3(ild3)=0
            END DO
            Da3(ild3)=Fa1(ieq)*Da2(ild2); ild3=ild3+1
          END DO
          Da3(ild3:ildend3)=0
        END DO

      ELSE
        
        DO ieq=1,neq
          DO ild2=Kld2(ieq),Kld2(ieq+1)-1
            Da3(ild2)=Fa1(ieq)*Da2(ild2)
          END DO
        END DO
        
      END IF
    END SUBROUTINE do_matDmat79mul_singledouble

    !**************************************************************
    ! Format D-7/9 multiplication
    ! double precision matrix A (format D)
    ! single precision matrix B (format 7 or format 9)
    ! double precision matrix C (format 7 or format 9)

    SUBROUTINE do_matDmat79mul_doublesingle(Da1,Fa2,Kld2,Kcol2&
        &,neq,Da3,Kld3,Kcol3)

      REAL(DP), DIMENSION(:), INTENT(IN)    :: Da1
      REAL(SP), DIMENSION(:), INTENT(IN)    :: Fa2
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: Da3
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld2
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol2
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kld3
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kcol3
      INTEGER(PREC_DOFIDX), INTENT(IN) :: neq

      INTEGER :: ieq,ild2,ild3,ildend3
      
      IF (PRESENT(Kld3) .AND. PRESENT(Kcol3)) THEN
        
        DO ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          DO ild2=Kld2(ieq),Kld2(ieq+1)-1
            DO ild3=ild3,ildend3
              IF (Kcol3(ild3) == Kcol2(ild2)) EXIT
              Da3(ild3)=0
            END DO
            Da3(ild3)=Da1(ieq)*Fa2(ild2); ild3=ild3+1
          END DO
          Da3(ild3:ildend3)=0
        END DO

      ELSE
        
        DO ieq=1,neq
          DO ild2=Kld2(ieq),Kld2(ieq+1)-1
            Da3(ild2)=Da1(ieq)*Fa2(ild2)
          END DO
        END DO
        
      END IF
    END SUBROUTINE do_matDmat79mul_doublesingle

    !**************************************************************
    ! Format D-7/9 multiplication
    ! single precision matrix A (format D)
    ! single precision matrix B (format 7 or format 9)
    ! single precision matrix C (format 7 or format 9)

    SUBROUTINE do_matDmat79mul_singlesingle(Fa1,Fa2,Kld2,Kcol2&
        &,neq,Fa3,Kld3,Kcol3)

      REAL(SP), DIMENSION(:), INTENT(IN)    :: Fa1
      REAL(SP), DIMENSION(:), INTENT(IN)    :: Fa2
      REAL(SP), DIMENSION(:), INTENT(INOUT) :: Fa3
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld2
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol2
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kld3
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kcol3
      INTEGER(PREC_DOFIDX), INTENT(IN) :: neq

      INTEGER :: ieq,ild2,ild3,ildend3
      
      IF (PRESENT(Kld3) .AND. PRESENT(Kcol3)) THEN
        
        DO ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          DO ild2=Kld2(ieq),Kld2(ieq+1)-1
            DO ild3=ild3,ildend3
              IF (Kcol3(ild3) == Kcol2(ild2)) EXIT
              Fa3(ild3)=0
            END DO
            Fa3(ild3)=Fa1(ieq)*Fa2(ild2); ild3=ild3+1
          END DO
          Fa3(ild3:ildend3)=0
        END DO

      ELSE
        
        DO ieq=1,neq
          DO ild2=Kld2(ieq),Kld2(ieq+1)-1
            Fa3(ild2)=Fa1(ieq)*Fa2(ild2)
          END DO
        END DO
        
      END IF
    END SUBROUTINE do_matDmat79mul_singlesingle

    !**************************************************************
    ! Format 7/9-D multiplication
    ! double precision matrix A (format 7 or format 9)
    ! double precision matrix B (format D)
    ! double precision matrix C (format 7 or format 9)
    
    SUBROUTINE do_mat79matDmul_doubledouble(Da1,Kld1,Kcol1,neq,Da2&
        &,Da3,Kld3,Kcol3)

      REAL(DP), DIMENSION(:), INTENT(IN)    :: Da1
      REAL(DP), DIMENSION(:), INTENT(IN)    :: Da2
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: Da3
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld1
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol1
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kld3
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kcol3
      INTEGER(PREC_DOFIDX), INTENT(IN) :: neq
      
      INTEGER :: ieq,ild1,ild3,ildend3
      
      IF (PRESENT(Kld3) .AND. PRESENT(Kcol3)) THEN
        
        DO ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          DO ild1=Kld1(ieq),Kld1(ieq+1)-1
            DO ild3=ild3,ildend3
              IF (Kcol3(ild3) == Kcol1(ild1)) EXIT
              Da3(ild3)=0
            END DO
            Da3(ild3)=Da1(ild1)*Da2(Kcol1(ild1)); ild3=ild3+1
          END DO
          Da3(ild3:ildend3)=0
        END DO
        
      ELSE
        
        DO ieq=1,neq
          DO ild1=Kld1(ieq),Kld1(ieq+1)-1
            Da3(ild1)=Da1(ild1)*Da2(Kcol1(ild1))
          END DO
        END DO

      END IF
    END SUBROUTINE do_mat79matDmul_doubledouble
      
    !**************************************************************
    ! Format 7/9-D multiplication
    ! single precision matrix A (format 7 or format 9)
    ! double precision matrix B (format D)
    ! double precision matrix C (format 7 or format 9)
    
    SUBROUTINE do_mat79matDmul_singledouble(Fa1,Kld1,Kcol1,neq,Da2&
        &,Da3,Kld3,Kcol3)

      REAL(SP), DIMENSION(:), INTENT(IN)    :: Fa1
      REAL(DP), DIMENSION(:), INTENT(IN)    :: Da2
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: Da3
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld1
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol1
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kld3
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kcol3
      INTEGER(PREC_DOFIDX), INTENT(IN) :: neq
      
      INTEGER :: ieq,ild1,ild3,ildend3
      
      IF (PRESENT(Kld3) .AND. PRESENT(Kcol3)) THEN
        
        DO ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          DO ild1=Kld1(ieq),Kld1(ieq+1)-1
            DO ild3=ild3,ildend3
              IF (Kcol3(ild3) == Kcol1(ild1)) EXIT
              Da3(ild3)=0
            END DO
            Da3(ild3)=Fa1(ild1)*Da2(Kcol1(ild1)); ild3=ild3+1
          END DO
          Da3(ild3:ildend3)=0
        END DO

      ELSE
        
        DO ieq=1,neq
          DO ild1=Kld1(ieq),Kld1(ieq+1)-1
            Da3(ild1)=Fa1(ild1)*Da2(Kcol1(ild1))
          END DO
        END DO

      END IF
    END SUBROUTINE do_mat79matDmul_singledouble

    !**************************************************************
    ! Format 7/9-D multiplication
    ! double precision matrix A (format 7 or format 9)
    ! single precision matrix B (format D)
    ! double precision matrix C (format 7 or format 9)
    
    SUBROUTINE do_mat79matDmul_doublesingle(Da1,Kld1,Kcol1,neq,Fa2&
        &,Da3,Kld3,Kcol3)

      REAL(DP), DIMENSION(:), INTENT(IN)    :: Da1
      REAL(SP), DIMENSION(:), INTENT(IN)    :: Fa2
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: Da3
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld1
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol1
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kld3
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kcol3
      INTEGER(PREC_DOFIDX), INTENT(IN) :: neq

      INTEGER :: ieq,ild1,ild3,ildend3
      
      IF (PRESENT(Kld3) .AND. PRESENT(Kcol3)) THEN
        
        DO ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          DO ild1=Kld1(ieq),Kld1(ieq+1)-1
            DO ild3=ild3,ildend3
              IF (Kcol3(ild3) == Kcol1(ild1)) EXIT
              Da3(ild3)=0
            END DO
            Da3(ild3)=Da1(ild1)*Fa2(Kcol1(ild1)); ild3=ild3+1
          END DO
          Da3(ild3:ildend3)=0
        END DO
        
      ELSE
        
        DO ieq=1,neq
          DO ild1=Kld1(ieq),Kld1(ieq+1)-1
            Da3(ild1)=Da1(ild1)*Fa2(Kcol1(ild1))
          END DO
        END DO

      END IF
    END SUBROUTINE do_mat79matDmul_doublesingle

    !**************************************************************
    ! Format 7/9-D multiplication
    ! single precision matrix A (format 7 or format 9)
    ! single precision matrix B (format D)
    ! single precision matrix C (format 7 or format 9)
    
    SUBROUTINE do_mat79matDmul_singlesingle(Fa1,Kld1,Kcol1,neq,Fa2&
        &,Fa3,Kld3,Kcol3)

      REAL(SP), DIMENSION(:), INTENT(IN)    :: Fa1
      REAL(SP), DIMENSION(:), INTENT(IN)    :: Fa2
      REAL(SP), DIMENSION(:), INTENT(INOUT) :: Fa3
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld1
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol1
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kld3
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kcol3
      INTEGER(PREC_DOFIDX), INTENT(IN) :: neq

      INTEGER :: ieq,ild1,ild3,ildend3

      IF (PRESENT(Kld3) .AND. PRESENT(Kcol3)) THEN
        
        DO ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          DO ild1=Kld1(ieq),Kld1(ieq+1)-1
            DO ild3=ild3,ildend3
              IF (Kcol3(ild3) == Kcol1(ild1)) EXIT
              Fa3(ild3)=0
            END DO
            Fa3(ild3)=Fa1(ild1)*Fa2(Kcol1(ild1)); ild3=ild3+1
          END DO
          Fa3(ild3:ildend3)=0
        END DO

      ELSE
        
        DO ieq=1,neq
          DO ild1=Kld1(ieq),Kld1(ieq+1)-1
            Fa3(ild1)=Fa1(ild1)*Fa2(Kcol1(ild1))
          END DO
        END DO

      END IF
    END SUBROUTINE do_mat79matDmul_singlesingle

    !**************************************************************
    ! Format 7/9-7/9 multiplication
    ! Compute the number of nonzero matrix entries of C:=A * B
    ! 
    ! Method: A' * A = sum [over i=1, nrow] a(i)^T a(i)
    !         where a(i) = i-th row of A. We must be careful not
    !         to add the elements already accounted for.
    !
    ! Remark: This subroutine is a modified version of the AMUBDG
    !         subroutine taken from the SPARSEKIT library written
    !         by Youcef Saad.

    FUNCTION do_mat79mat79mul_computeNA(KldA,KcolA,neq,KldB,KcolB&
        &,Kaux) RESULT(NA)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA,KcolA,KldB,KcolB
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Kaux
      INTEGER(PREC_DOFIDX), INTENT(IN) :: neq
      INTEGER :: NA

      INTEGER :: ieq,jeq,ild,irow,icol,idg,ndg,last

      ! Initialization
      Kaux=0; NA=0

      DO ieq=1,neq
        
        ! For each row of matrix A
        ndg=0
        
        ! End-of-linked list
        last=-1

        DO ild=KldA(ieq),KldA(ieq+1)-1

          ! Row number to be added
          irow=KcolA(ild)
          
          DO jeq=KldB(irow),KldB(irow+1)-1
            icol=KcolB(jeq)

            ! Add one element to the linked list
            IF (Kaux(icol) == 0) THEN
              ndg = ndg+1
              Kaux(icol) = last
              last = icol
            END IF
            
          END DO
        END DO

        ! We are done with row IEQ
        NA = NA+ndg
        
        ! Reset Kaux to zero
        DO idg=1,ndg
          jeq = Kaux(last)
          Kaux(last) = 0
          last = jeq
        END DO
        
      END DO
    END FUNCTION do_mat79mat79mul_computeNA

    !**************************************************************
    ! Format 7/9-7/9 multiplication
    ! Perform symbolical matrix-matrix-multiplication C := A * B
    ! 
    ! This subroutine is based on the SYMBMM routine from the
    ! Sparse Matrix Multiplication Package (SMMP) written by 
    ! R.E. Bank and C.C. Douglas which is freely available at:
    ! http://cs-www.cs.yale.edu/homes/douglas-craig/Codes/smmp.tgz

    SUBROUTINE do_mat79mat79mul_symb(n,m,l,KldA,KcolA,KldB,KcolB,KldC&
        &,KcolC,Kindex,KdiagonalC)

      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m,l
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA,KldB
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA,KcolB
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: KldC,Kindex
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: KcolC
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT), OPTIONAL :: KdiagonalC

      INTEGER :: i,j,k,istart,ilength,jj,jk,kaux,ioff

      ! KINDEX is used to store links. If an entry in KINDEX os
      ! nonzero, it is a pointer to the next column with a nonzero.
      ! The links are determined as they are found, and are unordered.
      Kindex = 0
      KldC(1)= 1

      ! The main loop consists of three components: initialization, a
      ! long loop that merges row lists, and code to copy the links
      ! into the KCOLC vector.
      DO i=1,n
        
        ! Initialization
        istart = -1
        ilength = 0

        ! The start column ISTART is reset and the number of column
        ! entries for the I-th row is assumed empty. Now, merge the
        ! row lists as follows
        DO jj=KldA(i),KldA(i+1)-1
          j=KcolA(jj)
          
          ! Determine the intersection of the current row I in matrix
          ! A with the nonzeros in column J of matrix B
          DO k=KldB(j),KldB(j+1)-1
            jk=KcolB(k)
            IF (Kindex(jk) == 0) THEN
              Kindex(jk) = istart
              istart = jk
              ilength = ilength+1
            END IF
          END DO
        END DO

        ! Copy the links into the KCOLC vector as column indices
        ! Since the column vectors are given in unordered way, they
        ! must be ordered before inserted inco KCOLC
        KldC(i+1)=KldC(i)+ilength
        
        IF (PRESENT(KdiagonalC)) THEN
          ! If KDIAGONALC is present, then the matrix C is stored in
          ! CSR9 format, that is, all entries in KCOLC are numbered
          ! continuously but the position of the diagonal entry is
          ! kept in KDIAGONALC
          DO j=KldC(i),KldC(i+1)-1
            KcolC(j) = istart
            istart = Kindex(istart)
            Kindex(KcolC(j)) = 0
            
            ! Perform local insertionsort to find the correct
            ! position of the provisional entry KCOLC(J)
            kaux = KcolC(j)
            DO jj=j-1,KldC(i),-1
              IF (KcolC(jj) <= kaux) GOTO 10
              KcolC(jj+1) = KcolC(jj)
              IF (KcolC(jj+1) == i) KdiagonalC(i) = jj+1
            END DO
            jj=KldC(i)-1
10          KcolC(jj+1) = kaux
            IF (KcolC(jj+1) == i) KdiagonalC(i) = jj+1
          END DO

        ELSE
          ! If KDIAGONALC is not present, then the matrix C is stored
          ! in CSR7 format, that is, all entries in KCOLC are
          ! numbered continuously except for the diagonal entry which
          ! is stored in the first position KLDC(I) of row I.
          
          ! The variable IOFF is used to indicate, whether the
          ! diagonal entry has already been inserted or not.
          ioff=0
          
          DO j=KldC(i),KldC(i+1)-1
            KcolC(j) = istart
            istart = Kindex(istart)
            Kindex(KcolC(j)) = 0
            
            IF (KcolC(j) == i) THEN
              ! The current entry is the diagonal entry, hence, shift
              ! all entries to the right by one and move the current
              ! entry to the first position
              KcolC(KldC(i):j)=CSHIFT(KcolC(KldC(i):j),-1)
              
              ! Set IOFF to one so that the first entry of the I-th
              ! row is neglected in subsequent iterations
              ioff=1
            ELSE
              ! Perform local insertsiosort to find the correct
              ! position of the provisional entry KCOLC(J)
              kaux = KcolC(j); jj=j
              DO jj=j-1,KldC(i),-1
                IF (KcolC(jj) <= kaux) GOTO 20
                KcolC(jj+1) = KcolC(jj)
              END DO
              jj=KldC(i)-1
20            KcolC(jj+1) = kaux
            END IF
          END DO
        END IF

        Kindex(i) = 0
      END DO
    END SUBROUTINE do_mat79mat79mul_symb

    !**************************************************************
    ! Format 7/9-7/9 multiplication
    ! Perform numerical matrix-matrix-multiplication C := A * B
    !
    ! double precision matrix A (format 7 or format 9)
    ! double precision matrix B (format 7 or format 9)
    ! double precision matrix C (format 7 or format 9)
    !
    ! This subroutine is based on the NUMBMM routine from the
    ! Sparse Matrix Multiplication Package (SMMP) written by 
    ! R.E. Bank and C.C. Douglas which is freely available at:
    ! http://cs-www.cs.yale.edu/homes/douglas-craig/Codes/smmp.tgz
    
    SUBROUTINE do_mat79mat79mul_numb_dbledble(n,m,l,KldA,KcolA&
        &,DaA,KldB,KcolB,DaB,KldC,KcolC,DaC,Dtemp)

      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m,l
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA,KldB
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA,KcolB
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: KldC
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: KcolC
      REAL(DP), DIMENSION(:), INTENT(IN)    :: DaA
      REAL(DP), DIMENSION(:), INTENT(IN)    :: DaB
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: DaC,Dtemp

      INTEGER :: i,j,k,jj,jk
      REAL(DP) :: daux

      ! DTEMP is used to store partial sums.
      Dtemp = 0

      ! The main loop forms the partial sums and then copies the
      ! completed sums into the correct locations in the sparse matrix
      DO i=1,n
        
        DO jj=KldA(i),KldA(i+1)-1
          j = KcolA(jj)
          daux = DaA(jj)
          
          ! Accumulate the product for row J
          DO k=KldB(j),KldB(j+1)-1
            jk=KcolB(k)
            Dtemp(jk)=Dtemp(jk)+daux*DaB(k)
          END DO
        END DO
        
        ! Store product for row J
        DO j=KldC(i),KldC(i+1)-1
          jj=KcolC(j)
          DaC(j) = Dtemp(jj)
          Dtemp(jj) = 0
        END DO
      END DO
    END SUBROUTINE do_mat79mat79mul_numb_dbledble

    !**************************************************************
    ! Format 7/9-7/9 multiplication
    ! Perform numerical matrix-matrix-multiplication C := A * B
    !
    ! single precision matrix A (format 7 or format 9)
    ! double precision matrix B (format 7 or format 9)
    ! double precision matrix C (format 7 or format 9)
    !
    ! This subroutine is based on the NUMBMM routine from the
    ! Sparse Matrix Multiplication Package (SMMP) written by 
    ! R.E. Bank and C.C. Douglas which is freely available at:
    ! http://cs-www.cs.yale.edu/homes/douglas-craig/Codes/smmp.tgz
    
    SUBROUTINE do_mat79mat79mul_numb_sngldble(n,m,l,KldA,KcolA&
        &,FaA,KldB,KcolB,DaB,KldC,KcolC,DaC,Dtemp)

      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m,l
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA,KldB
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA,KcolB
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: KldC
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: KcolC
      REAL(SP), DIMENSION(:), INTENT(IN)    :: FaA
      REAL(DP), DIMENSION(:), INTENT(IN)    :: DaB
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: DaC,Dtemp

      INTEGER :: i,j,k,jj,jk
      REAL(DP) :: daux

      ! DTEMP is used to store partial sums.
      Dtemp = 0

      ! The main loop forms the partial sums and then copies the
      ! completed sums into the correct locations in the sparse matrix
      DO i=1,n
        
        DO jj=KldA(i),KldA(i+1)-1
          j = KcolA(jj)
          daux = FaA(jj)
          
          ! Accumulate the product for row J
          DO k=KldB(j),KldB(j+1)-1
            jk=KcolB(k)
            Dtemp(jk)=Dtemp(jk)+daux*DaB(k)
          END DO
        END DO
        
        ! Store product for row J
        DO j=KldC(i),KldC(i+1)-1
          jj=KcolC(j)
          DaC(j) = Dtemp(jj)
          Dtemp(jj) = 0
        END DO
      END DO
    END SUBROUTINE do_mat79mat79mul_numb_sngldble

    !**************************************************************
    ! Format 7/9-7/9 multiplication
    ! Perform numerical matrix-matrix-multiplication C := A * B
    !
    ! double precision matrix A (format 7 or format 9)
    ! single precision matrix B (format 7 or format 9)
    ! double precision matrix C (format 7 or format 9)
    !
    ! This subroutine is based on the NUMBMM routine from the
    ! Sparse Matrix Multiplication Package (SMMP) written by 
    ! R.E. Bank and C.C. Douglas which is freely available at:
    ! http://cs-www.cs.yale.edu/homes/douglas-craig/Codes/smmp.tgz
    
    SUBROUTINE do_mat79mat79mul_numb_dblesngl(n,m,l,KldA,KcolA&
        &,DaA,KldB,KcolB,FaB,KldC,KcolC,DaC,Dtemp)

      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m,l
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA,KldB
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA,KcolB
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: KldC
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: KcolC
      REAL(DP), DIMENSION(:), INTENT(IN)    :: DaA
      REAL(SP), DIMENSION(:), INTENT(IN)    :: FaB
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: DaC,Dtemp

      INTEGER :: i,j,k,jj,jk
      REAL(DP) :: daux

      ! DTEMP is used to store partial sums.
      Dtemp = 0

      ! The main loop forms the partial sums and then copies the
      ! completed sums into the correct locations in the sparse matrix
      DO i=1,n
        
        DO jj=KldA(i),KldA(i+1)-1
          j = KcolA(jj)
          daux = DaA(jj)
          
          ! Accumulate the product for row J
          DO k=KldB(j),KldB(j+1)-1
            jk=KcolB(k)
            Dtemp(jk)=Dtemp(jk)+daux*FaB(k)
          END DO
        END DO
        
        ! Store product for row J
        DO j=KldC(i),KldC(i+1)-1
          jj=KcolC(j)
          DaC(j) = Dtemp(jj)
          Dtemp(jj) = 0
        END DO
      END DO
    END SUBROUTINE do_mat79mat79mul_numb_dblesngl

    !**************************************************************
    ! Format 7/9-7/9 multiplication
    ! Perform numerical matrix-matrix-multiplication C := A * B
    !
    ! single precision matrix A (format 7 or format 9)
    ! single precision matrix B (format 7 or format 9)
    ! single precision matrix C (format 7 or format 9)
    !
    ! This subroutine is based on the NUMBMM routine from the
    ! Sparse Matrix Multiplication Package (SMMP) written by 
    ! R.E. Bank and C.C. Douglas which is freely available at:
    ! http://cs-www.cs.yale.edu/homes/douglas-craig/Codes/smmp.tgz
    
    SUBROUTINE do_mat79mat79mul_numb_snglsngl(n,m,l,KldA,KcolA&
        &,FaA,KldB,KcolB,FaB,KldC,KcolC,FaC,Ftemp)

      INTEGER(PREC_DOFIDX), INTENT(IN) :: n,m,l
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA,KldB
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA,KcolB
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: KldC
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: KcolC
      REAL(SP), DIMENSION(:), INTENT(IN)    :: FaA
      REAL(SP), DIMENSION(:), INTENT(IN)    :: FaB
      REAL(SP), DIMENSION(:), INTENT(INOUT) :: FaC,Ftemp

      INTEGER :: i,j,k,jj,jk
      REAL(SP) :: faux

      ! FTEMP is used to store partial sums.
      Ftemp = 0

      ! The main loop forms the partial sums and then copies the
      ! completed sums into the correct locations in the sparse matrix
      DO i=1,n
        
        DO jj=KldA(i),KldA(i+1)-1
          j = KcolA(jj)
          faux = FaA(jj)
          
          ! Accumulate the product for row J
          DO k=KldB(j),KldB(j+1)-1
            jk=KcolB(k)
            Ftemp(jk)=Ftemp(jk)+faux*FaB(k)
          END DO
        END DO
        
        ! Store product for row J
        DO j=KldC(i),KldC(i+1)-1
          jj=KcolC(j)
          FaC(j) = Ftemp(jj)
          Ftemp(jj) = 0
        END DO
      END DO
    END SUBROUTINE do_mat79mat79mul_numb_snglsngl
  END SUBROUTINE lsyssc_multMatMat

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE lsyssc_matrixLinearComb (rmatrixA,cA,rmatrixB,cB,rmatrixC,&
      bmemory,bsymb,bnumb,bisExactStructure)

    !<description>
    ! Adds constant times a matrix to another matrix
    !   rmatrixC = cA*rmatrixA + cB*rmatrixB
    !
    ! All matrix formats and most combinations of matrix formats are
    ! supported for matrices A and B. The resulting matrix C=ca*A+cb*B
    ! is generated in the "encompassing" format. That is, if A and B are
    ! both diagonal matrices, then C is also a diagonal matrix. If
    ! for example, A is a full matrix and B is a diagonal matrix (or
    ! vice versa) then C is also a full matrix so that all matrix
    ! entries can be stored.
    ! Moreover, rmatrixC can be the same as rmatrixA or rmatrixB, but
    ! then, BMEMORY and BSYMB must both be FALSE and the user must
    ! take care, the matrix C corresponds to the "encompassing" matrix.
    !</description>

!<input>
    ! source matrices
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrixA,rmatrixB

    ! scaling factors
    REAL(DP), INTENT(IN) :: cA,cB

    ! BMEMORY = FALSE: Do not allocate required memory for C=A+B.
    ! BMEMORY = TRUE:  Generate all required structures for C=A+B 
    LOGICAL, INTENT(IN) :: bmemory

    ! Compute symbolic matrix-matrix-addition
    ! BSYMB = FALSE: Do not generate the required matrix structures.
    !                This may be useful, if the same matrices A and B
    !                need to be added several times, but the
    !                sparsity patterns do not change
    ! BSYMN = TRUE:  Generate all required matrix structures for C=A+B
    LOGICAL, INTENT(IN) :: bsymb

    ! Compute numerical matrix-matrix-addition
    ! BNUMB = FALSE: Do not perform the numerical addition
    ! BNUMB = TRUE:  Perform numerical addition
    LOGICAL, INTENT(IN) :: bnumb
    
    ! OPTIONAL: Indicates whether the resulting matrix C has the
    ! required symbolic structure or is a superset of the symbolic
    ! matrix-matrix product. In some cases, this may allow for much
    ! more efficient implementation.
    ! Standard: FALSE
    LOGICAL, INTENT(IN), OPTIONAL :: bisExactStructure
!</input>

!<inputoutput>
    ! Output matrix. May coincode with rmatrixB.
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixC
!</inputoutput>

!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER :: DaA,DaB,DaC
    REAL(SP), DIMENSION(:), POINTER :: FaA,FaB,FaC
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: KldA,KldB,KldC,KdiagonalC,Kaux
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: KcolA,KcolB,KcolC
    INTEGER :: h_Kaux,isizeIntl
    LOGICAL :: bfast

    h_Kaux=ST_NOHANDLE
    
    ! Check if fast implementation can be used
    bfast=.FALSE.
    IF (PRESENT(bisExactStructure)) bfast = bisExactStructure
    bfast = bfast .OR. (bmemory .AND. bsymb)

    ! Check if both matrices are compatible
    IF (rmatrixA%NEQ /= rmatrixB%NEQ .OR. &
        & rmatrixA%NCOLS /= rmatrixB%NCOLS) THEN
      PRINT *, 'lsyssc_matrixLinearComb: number of rows/columns of matrix A &
          &is not compatible with number of rows/columns of matrix B'
      CALL sys_halt()
    END IF

    ! Check if both matrices have the same sorting
    IF (rmatrixA%isortStrategy /= rmatrixB%isortStrategy) THEN
      PRINT *, 'lsyssc_matrixLinearComb: incompatible sorting strategies'
      CALL sys_halt()
    END IF

    ! Release matrix if required and set common variables
    IF (bmemory) THEN
      CALL lsyssc_releaseMatrix(rmatrixC)
      IF (rmatrixA%cdataType == ST_DOUBLE .OR. rmatrixB%cdataType == ST_DOUBLE) THEN
        rmatrixC%cdataType = ST_DOUBLE
      ELSE
        rmatrixC%cdataType = ST_SINGLE
      END IF
    END IF
    
    ! Set sorting strategy for matrix C
    rmatrixC%isortStrategy = rmatrixA%isortStrategy
    rmatrixC%h_IsortPermutation = rmatrixA%h_IsortPermutation

    SELECT CASE(rmatrixA%cmatrixFormat)

    CASE (LSYSSC_MATRIX1)   ! A is full matrix ---------------------------------

      SELECT CASE(rmatrixB%cmatrixFormat)

      CASE (LSYSSC_MATRIX1) ! B is full matrix  - - - - - - - - - - - - - - - -
        
        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = LSYSSC_MATRIX1
          rmatrixC%NA = rmatrixA%NEQ*rmatrixA%NCOLS
          CALL storage_new('lsyssc_matrixLinearComb','h_Da',rmatrixC%NA,&
              rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= LSYSSC_MATRIX1) THEN
          PRINT *, 'lsyssc_matrixLinearComb: destination matrix has incompatible format'
          CALL sys_halt()
        END IF
        
        ! symbolic addition?
        IF (bsymb) THEN
          rmatrixC%NEQ = rmatrixA%NEQ
          rmatrixC%NCOLS = rmatrixA%NCOLS
        END IF
        
        ! numerical addition?
        IF (bnumb) THEN
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          SELECT CASE(rmatrixA%cdataType)

          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lalg_copyVectorDble(DaA,DaC)
              CALL lalg_vectorLinearCombDble(DaB,DaC,cB,cA)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lalg_copyVectorDble(DaA,DaC)
              CALL lalg_vectorLinearCombDble(REAL(FaB,DP),DaC,cB,cA)

            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lalg_copyVectorDble(REAL(FaA,DP),DaC)
              CALL lalg_vectorLinearCombDble(DaB,DaC,cB,cA)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL lalg_copyVectorSngl(FaA,FaC)
              CALL lalg_vectorLinearCombSngl(FaB,FaC,REAL(cB,SP),REAL(cA,SP))

            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE DEFAULT
            PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF
        
      CASE (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - - - - - - - -
        
        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = LSYSSC_MATRIX1
          rmatrixC%NA = rmatrixA%NEQ*rmatrixA%NCOLS
          CALL storage_new('lsyssc_matrixLinearComb','h_Da',rmatrixC%NA,&
              rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= LSYSSC_MATRIX1) THEN
          PRINT *, 'lsyssc_matrixLinearComb: destination matrix has incompatible format'
          CALL sys_halt()
        END IF
        
        ! symbolical addition?
        IF (bsymb) THEN
          rmatrixC%NEQ = rmatrixA%NEQ
          rmatrixC%NCOLS = rmatrixA%NCOLS
        END IF
        
        ! numerical addition?
        IF (bnumb) THEN
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double          
          SELECT CASE(rmatrixA%cdataType)

          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat1matDadd_doubledouble(rmatrixA%NEQ,rmatrixA%NCOLS,DaA,cA,DaB,cB,DaC)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat1matDadd_doublesingle(rmatrixA%NEQ,rmatrixA%NCOLS,DaA,cA,FaB,cB,DaC)

            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat1matDadd_singledouble(rmatrixA%NEQ,rmatrixA%NCOLS,FaA,cA,DaB,cB,DaC)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL do_mat1matDadd_singlesingle(rmatrixA%NEQ,rmatrixA%NCOLS,FaA,cA,FaB,cB,FaC)

            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE DEFAULT
            PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! B is CSR matrix - - - - - - - - - -
        
        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = LSYSSC_MATRIX1
          rmatrixC%NA = rmatrixA%NEQ*rmatrixA%NCOLS
          CALL storage_new('lsyssc_matrixLinearComb','h_Da',rmatrixC%NA,&
              rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= LSYSSC_MATRIX1) THEN
          PRINT *, 'lsyssc_matrixLinearComb: destination matrix has incompatible format'
          CALL sys_halt()
        END IF
        
        ! symbolical addition?
        IF (bsymb) THEN
          rmatrixC%NEQ = rmatrixA%NEQ
          rmatrixC%NCOLS = rmatrixA%NCOLS
        END IF
        
        ! numerical addition?
        IF (bnumb) THEN

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double          
          SELECT CASE(rmatrixA%cdataType)
            
          CASE (ST_DOUBLE)
            
            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              CALL do_mat1mat79add_doubledouble(rmatrixA%NEQ,rmatrixA%NCOLS,&
                  DaA,cA,KldB,KcolB,DaB,cB,DaC)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              CALL do_mat1mat79add_doublesingle(rmatrixA%NEQ,rmatrixA%NCOLS,&
                  DaA,cA,KldB,KcolB,FaB,cB,DaC)

            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              CALL do_mat1mat79add_singledouble(rmatrixA%NEQ,rmatrixA%NCOLS,&
                  FaA,cA,KldB,KcolB,DaB,cB,DaC)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              CALL do_mat1mat79add_singlesingle(rmatrixA%NEQ,rmatrixA%NCOLS,&
                  FaA,cA,KldB,KcolB,FaB,cB,FaC)

            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE DEFAULT
            PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE DEFAULT
        PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
        CALL sys_halt()
      END SELECT

      !-------------------------------------------------------------------------

    CASE (LSYSSC_MATRIXD)   ! A is diagonal matrix -----------------------------

      SELECT CASE(rmatrixB%cmatrixFormat)

      CASE (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - - - - - - - -
        
        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = LSYSSC_MATRIXD
          rmatrixC%NA = rmatrixA%NA
          CALL storage_new('lsyssc_matrixLinearComb','h_Da',rmatrixC%NA,&
              rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= LSYSSC_MATRIXD) THEN
          PRINT *, 'lsyssc_matrixLinearComb: destination matrix has incompatible format'
          CALL sys_halt()
        END IF
        
        ! symbolical addition?
        IF (bsymb) THEN
          rmatrixC%NEQ = rmatrixA%NEQ
          rmatrixC%NCOLS = rmatrixA%NCOLS
        END IF
        
        ! numerical addition?
        IF (bnumb) THEN
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double          
          SELECT CASE(rmatrixA%cdataType)
            
          CASE (ST_DOUBLE)
            
            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              DaC=cA*DaA+cb*DaB

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              DaC=cA*DaA+cb*FaB
              
            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              DaC=cA*FaA+cb*DaB

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              FaC=cA*FaA+cb*FaB

            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE DEFAULT
            PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE (LSYSSC_MATRIX1) ! B is full matrix - - - - - - - - - - -

        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = LSYSSC_MATRIX1
          rmatrixC%NA = rmatrixB%NEQ*rmatrixB%NCOLS
          CALL storage_new('lsyssc_matrixLinearComb','h_Da',rmatrixC%NA,&
              rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= LSYSSC_MATRIX1) THEN
          PRINT *, 'lsyssc_matrixLinearComb: destination matrix has incompatible format'
          CALL sys_halt()
        END IF
        
        ! symbolical addition?
        IF (bsymb) THEN
          rmatrixC%NEQ = rmatrixB%NEQ
          rmatrixC%NCOLS = rmatrixB%NCOLS
        END IF
        
        ! numerical addition?
        IF (bnumb) THEN

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double          
          SELECT CASE(rmatrixA%cdataType)

          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat1matDadd_doubledouble(rmatrixB%NEQ,rmatrixB%NCOLS,DaB,cB,DaA,cA,DaC)
              
            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat1matDadd_singledouble(rmatrixB%NEQ,rmatrixB%NCOLS,FaB,cB,DaA,cA,DaC)

            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL do_mat1matDadd_doublesingle(rmatrixA%NEQ,rmatrixA%NCOLS,DaB,cB,FaA,cA,DaC)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL do_mat1matDadd_singlesingle(rmatrixB%NEQ,rmatrixB%NCOLS,FaB,cB,FaA,cA,FaC)

            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE DEFAULT
            PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF
        
      CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! B is CSR matrix - - - - - - - - - -
        
        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = rmatrixB%cmatrixFormat
          rmatrixC%NA = rmatrixB%NA
          CALL storage_new('lsyssc_matrixLinearComb','h_Da',rmatrixC%NA,&
              rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= rmatrixB%cmatrixFormat) THEN
          PRINT *, 'lsyssc_matrixLinearComb: destination matrix has incompatible format'
          CALL sys_halt()
        END IF
        
        ! symbolic addition?
        IF (bsymb) THEN
          CALL lsyssc_duplicateMatrix(rmatrixB,rmatrixC,LSYSSC_DUP_COPY,LSYSSC_DUP_IGNORE)
        END IF
        
        ! numerical addition?
        IF (bnumb) THEN
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          SELECT CASE(rmatrixA%cdataType)

          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                CALL do_mat79matDadd_doubledouble(KldB,KcolB,1,1,rmatrixB%NEQ,&
                    DaB,cB,DaA,cA,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_mat79matDadd_doubledouble(KldB,KcolB,1,1,rmatrixB%NEQ,&
                    DaB,cB,DaA,cA,DaC,KldC,KcolC)
              END IF
              
            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                CALL do_mat79matDadd_singledouble(KldB,KcolB,1,1,rmatrixB%NEQ,&
                    FaB,cB,DaA,cA,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_mat79matDadd_singledouble(KldB,KcolB,1,1,rmatrixB%NEQ,&
                    FaB,cB,DaA,cA,DaC,KldC,KcolC)
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE (ST_SINGLE)
            
            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                CALL do_mat79matDadd_doublesingle(KldB,KcolB,1,1,rmatrixB%NEQ,&
                    DaB,cB,FaA,cA,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_mat79matDadd_doublesingle(KldB,KcolB,1,1,rmatrixB%NEQ,&
                    DaB,cB,FaA,cA,DaC,KldC,KcolC)
              END IF
              
            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                CALL do_mat79matDadd_singlesingle(KldB,KcolB,1,1,rmatrixB%NEQ,&
                    FaB,cB,FaA,cA,FaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_mat79matDadd_singlesingle(KldB,KcolB,1,1,rmatrixB%NEQ,&
                    FaB,cB,FaA,cA,FaC,KldC,KcolC)
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE DEFAULT
            PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE (LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9INTL) ! B is interleave CSR matrix
        
        ! Set size of interleave block
        isizeIntl=rmatrixB%NVAR
        IF (rmatrixB%cinterleavematrixFormat == LSYSSC_MATRIX1) isizeIntl=isizeIntl*isizeIntl
        
        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = rmatrixB%cmatrixFormat
          rmatrixC%cinterleavematrixFormat = rmatrixB%cinterleavematrixFormat
          rmatrixC%NA = rmatrixB%NA
          rmatrixC%NVAR = rmatrixB%NVAR
          CALL storage_new('lsyssc_matrixLinearComb','h_Da',INT(isizeIntl*rmatrixC%NA,I32),&
              rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF

        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= rmatrixB%cmatrixFormat .OR. &
            rmatrixC%cinterleavematrixFormat /= rmatrixB%cinterleavematrixFormat) THEN
          PRINT *, 'lsyssc_matrixLinearComb: destination matrix has incompatible format'
          CALL sys_halt()
        END IF

        ! symbolic addition?
        IF (bsymb) THEN
          CALL lsyssc_duplicateMatrix(rmatrixB,rmatrixC,LSYSSC_DUP_COPY,LSYSSC_DUP_IGNORE)
        END IF

        ! numerical addition?
        IF (bnumb) THEN

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          SELECT CASE(rmatrixA%cdataType)
            
          CASE (ST_DOUBLE)
            
            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7INTL) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_doubledouble(KldB,KcolB,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixB%NEQ,DaB,cB,DaA,cA,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                ELSE
                  CALL do_mat79matDadd_doubledouble(KldB,KcolB,rmatrixC%NVAR,1,&
                      rmatrixB%NEQ,DaB,cB,DaA,cA,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                END IF

              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_doubledouble(KldB,KcolB,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixB%NEQ,DaB,cB,DaA,cA,DaC,KldC,KcolC)
                ELSE
                  CALL do_mat79matDadd_doubledouble(KldB,KcolB,rmatrixC%NVAR,1,&
                      rmatrixB%NEQ,DaB,cB,DaA,cA,DaC,KldC,KcolC)
                END IF
              END IF
              
            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7INTL) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_singledouble(KldB,KcolB,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixB%NEQ,FaB,cB,DaA,cA,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                ELSE
                  CALL do_mat79matDadd_singledouble(KldB,KcolB,rmatrixC%NVAR,1,&
                      rmatrixB%NEQ,FaB,cB,DaA,cA,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                END IF

              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_singledouble(KldB,KcolB,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixB%NEQ,FaB,cB,DaA,cA,DaC,KldC,KcolC)
                ELSE
                  CALL do_mat79matDadd_singledouble(KldB,KcolB,rmatrixC%NVAR,1,&
                      rmatrixB%NEQ,FaB,cB,DaA,cA,DaC,KldC,KcolC)
                END IF
              END IF

            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)
            
            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7INTL) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_doublesingle(KldB,KcolB,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixB%NEQ,DaB,cB,FaA,cA,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                ELSE
                  CALL do_mat79matDadd_doublesingle(KldB,KcolB,rmatrixC%NVAR,1,&
                      rmatrixB%NEQ,DaB,cB,FaA,cA,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                END IF

              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_doublesingle(KldB,KcolB,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixB%NEQ,DaB,cB,FaA,cA,DaC,KldC,KcolC)
                ELSE
                  CALL do_mat79matDadd_doublesingle(KldB,KcolB,rmatrixC%NVAR,1,&
                      rmatrixB%NEQ,DaB,cB,FaA,cA,DaC,KldC,KcolC)
                END IF
              END IF
              
            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL lsyssc_getbase_Kld(rmatrixB,KldB)
              CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7INTL) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_singlesingle(KldB,KcolB,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixB%NEQ,FaB,cB,FaA,cA,FaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                ELSE
                  CALL do_mat79matDadd_singlesingle(KldB,KcolB,rmatrixC%NVAR,1,&
                      rmatrixB%NEQ,FaB,cB,FaA,cA,FaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                END IF

              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_singlesingle(KldB,KcolB,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixB%NEQ,FaB,cB,FaA,cA,FaC,KldC,KcolC)
                ELSE
                  CALL do_mat79matDadd_singlesingle(KldB,KcolB,rmatrixC%NVAR,1,&
                      rmatrixB%NEQ,FaB,cB,FaA,cA,FaC,KldC,KcolC)
                END IF
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE DEFAULT
            PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE DEFAULT
        PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
        CALL sys_halt()
      END SELECT

      !-------------------------------------------------------------------------

    CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! A is CSR matrix ---------------------
      
      SELECT CASE(rmatrixB%cmatrixFormat)
        
      CASE (LSYSSC_MATRIX1) ! B is full matrix - - - - - - - - - - - - - - - - -
        
        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = LSYSSC_MATRIX1
          rmatrixC%NA = rmatrixB%NEQ*rmatrixB%NCOLS
          CALL storage_new('lsyssc_matrixLinearComb','h_Da',rmatrixC%NA,&
              rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= LSYSSC_MATRIX1) THEN
          PRINT *, 'lsyssc_matrixLinearComb: destination matrix has incompatible format'
          CALL sys_halt()
        END IF
        
        ! symbolical addition?
        IF (bsymb) THEN
          rmatrixC%NEQ = rmatrixB%NEQ
          rmatrixC%NCOLS = rmatrixB%NCOLS
        END IF
        
        ! numerical addition?
        IF (bnumb) THEN
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double          
          SELECT CASE(rmatrixA%cdataType)
            
          CASE (ST_DOUBLE)
            
            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              CALL do_mat1mat79add_doubledouble(rmatrixB%NEQ,rmatrixB%NCOLS,&
                  DaB,cB,KldA,KcolA,DaA,cA,DaC)

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              CALL do_mat1mat79add_singledouble(rmatrixB%NEQ,rmatrixB%NCOLS,&
                  FaB,cB,KldA,KcolA,DaA,cA,DaC)
              
            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE (ST_SINGLE)
            
            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              CALL do_mat1mat79add_doublesingle(rmatrixB%NEQ,rmatrixB%NCOLS,&
                  DaB,cB,KldA,KcolA,FaA,cA,DaC)
              
            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              CALL do_mat1mat79add_singlesingle(rmatrixB%NEQ,rmatrixB%NCOLS,&
                  FaB,cB,KldA,KcolA,FaA,cA,FaC)
              
            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE DEFAULT
            PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - - - - - - - -
      
        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = rmatrixA%cmatrixFormat
          rmatrixC%NA = rmatrixA%NA
          CALL storage_new('lsyssc_matrixLinearComb','h_Da',rmatrixC%NA,&
              rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= rmatrixA%cmatrixFormat) THEN
          PRINT *, 'lsyssc_matrixLinearComb: destination matrix has incompatible format'
          CALL sys_halt()
        END IF
        
        ! symbolic addition?
        IF (bsymb) THEN
          CALL lsyssc_duplicateMatrix(rmatrixA,rmatrixC,LSYSSC_DUP_COPY,LSYSSC_DUP_IGNORE)
        END IF
        
        ! numerical addition?
        IF (bnumb) THEN

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          SELECT CASE(rmatrixA%cdataType)

          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                CALL do_mat79matDadd_doubledouble(KldA,KcolA,1,1,rmatrixA%NEQ,&
                    DaA,cA,DaB,cB,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_mat79matDadd_doubledouble(KldA,KcolA,1,1,rmatrixA%NEQ,&
                    DaA,cA,DaB,cB,DaC,KldC,KcolC)
              END IF
              
            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                CALL do_mat79matDadd_doublesingle(KldA,KcolA,1,1,rmatrixA%NEQ,&
                    DaA,cA,FaB,cB,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_mat79matDadd_doublesingle(KldA,KcolA,1,1,rmatrixA%NEQ,&
                    DaA,cA,FaB,cB,DaC,KldC,KcolC)
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE (ST_SINGLE)
            
            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                CALL do_mat79matDadd_singledouble(KldA,KcolA,1,1,rmatrixA%NEQ,&
                    FaA,cA,DaB,cB,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_mat79matDadd_singledouble(KldA,KcolA,1,1,rmatrixA%NEQ,&
                    FaA,cA,DaB,cB,DaC,KldC,KcolC)
              END IF

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                CALL do_mat79matDadd_singlesingle(KldA,KcolA,1,1,rmatrixA%NEQ,&
                    FaA,cA,FaB,cB,FaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                CALL do_mat79matDadd_singlesingle(KldA,KcolA,1,1,rmatrixA%NEQ,&
                    FaA,cA,FaB,cB,FaC,KldC,KcolC)
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE DEFAULT
            PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! B is CSR matrix - - - - - - - - - -
                
        ! Set pointers
        CALL lsyssc_getbase_Kld(rmatrixA,KldA)
        CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
        CALL lsyssc_getbase_Kld(rmatrixB,KldB)
        CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
        
        ! memory allocation?
        IF (bmemory) THEN
          ! Duplicate structure of matrix A or B depending on which
          ! matrix is given in the encompassing matrix format
          IF (rmatrixA%cmatrixFormat >= rmatrixB%cmatrixFormat) THEN
            CALL lsyssc_duplicateMatrix(rmatrixA,rmatrixC,LSYSSC_DUP_EMPTY,LSYSSC_DUP_REMOVE)
          ELSE
            CALL lsyssc_duplicateMatrix(rmatrixB,rmatrixC,LSYSSC_DUP_EMPTY,LSYSSC_DUP_REMOVE)
          END IF
          
          ! Set auxiliary pointers
          CALL storage_new('lsyssc_matrixLinearComb','Kaux',rmatrixB%NCOLS,&
              &ST_INT,h_Kaux,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int(h_Kaux,Kaux)
          
          ! Compute number of nonzero matrix entries: NA
          rmatrixC%NA=do_mat79mat79add_computeNA(rmatrixA%NEQ,rmatrixA%NCOLS,KldA,KcolA,KldB,KcolB,Kaux)
          CALL storage_new('lsyssc_matrixLinearComb','h_Da',rmatrixC%NA,&
              rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
          CALL storage_realloc('lsyssc_matrixLinearComb',rmatrixC%NA,&
              rmatrixC%h_Kcol,ST_NEWBLOCK_NOINIT,.FALSE.)
        END IF
        
        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= MAX(rmatrixA%cmatrixFormat,rmatrixB%cmatrixFormat)) THEN
          PRINT *, 'lsyssc_matrixLinearComb: destination matrix has incompatible format'
          CALL sys_halt()
        END IF
        
        ! symbolic addition?
        IF (bsymb) THEN
          
          ! Set pointers
          CALL lsyssc_getbase_Kld(rmatrixC,KldC)
          CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                    
          IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX9) THEN
            CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
            CALL do_mat79mat79add_symb(rmatrixC%NEQ,rmatrixC%NCOLS,KldA,KcolA,&
                rmatrixA%cmatrixFormat,KldB,KcolB,rmatrixB%cmatrixFormat,KldC,KcolC,KdiagonalC)
          ELSE
            CALL do_mat79mat79add_symb(rmatrixC%NEQ,rmatrixC%NCOLS,KldA,KcolA,&
                rmatrixA%cmatrixFormat,KldB,KcolB,rmatrixB%cmatrixFormat,KldC,KcolC)
          END IF
        END IF
        
        ! numerical addition?
        IF (bnumb) THEN
          
          ! Set pointers
          CALL lsyssc_getbase_Kld(rmatrixC,KldC)
          CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double   
          SELECT CASE(rmatrixA%cdataType)
            
          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              
              IF ((rmatrixA%cmatrixFormat == rmatrixB%cmatrixFormat) .AND. &
                  (rmatrixA%cmatrixFormat == rmatrixC%cmatrixFormat) .AND. (bfast)) THEN
                
                ! We rely on the user who tells us that we should assume the
                ! same structure for the matrices. In this case, we can
                ! directly call a BLAS routine to do the matrix combination.
                ! That's MUCH faster!
                ! Note that matrix B might coincide with matrix C, so we
                ! first check this to prevent unnecessary copy operations!

                IF (.NOT. ASSOCIATED (DaB,DaC)) THEN
                  CALL lalg_copyVectorDble (DaB,DaC)
                END IF
                
                CALL lalg_vectorLinearCombDble (DaA,DaC,cA,cB)                
                
              ELSE IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX9) THEN
                
                CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                CALL do_mat79mat79add_numb_dbledble(1,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,DaA,cA,KldB,KcolB,DaB,cB,KldC,KcolC,KdiagonalC,DaC)
                
              ELSE
                
                CALL do_mat79mat79add_numb_dbledble(1,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,DaA,cA,KldB,KcolB,DaB,cB,KldC,KcolC,KldC,DaC)
                    
              END IF
              
            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              
              IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX9) THEN
                
                CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                CALL do_mat79mat79add_numb_dblesngl(1,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,DaA,cA,KldB,KcolB,FaB,cB,KldC,KcolC,KdiagonalC,DaC)
                
              ELSE
                
                CALL do_mat79mat79add_numb_dblesngl(1,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,DaA,cA,KldB,KcolB,FaB,cB,KldC,KcolC,KldC,DaC)
                
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE (ST_SINGLE)
            
            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE) 
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)

              IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX9) THEN
                
                CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                CALL do_mat79mat79add_numb_sngldble(1,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,FaA,cA,KldB,KcolB,DaB,cB,KldC,KcolC,KdiagonalC,DaC)

              ELSE
                
                CALL do_mat79mat79add_numb_sngldble(1,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,FaA,cA,KldB,KcolB,DaB,cB,KldC,KcolC,KldC,DaC)
                
              END IF
              
            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              
              IF ((rmatrixA%cmatrixFormat == rmatrixB%cmatrixFormat) .AND. &
                  (rmatrixA%cmatrixFormat == rmatrixC%cmatrixFormat) .AND. (bfast)) THEN
                  
                ! We rely on the user who tells us that we should assume the
                ! same structure for the matrices. In this case, we can
                ! directly call a BLAS routine to do the matrix combination.
                ! That's MUCH faster!
                ! Note that matrix B might coincide with matrix C, so we
                ! first check this to prevent unnecessary copy operations!
                
                IF (.NOT. ASSOCIATED (DaB,DaC)) THEN
                  CALL lalg_copyVectorSngl (FaB,FaC)
                END IF
                
                CALL lalg_vectorLinearCombSngl (FaA,FaC,REAL(cA,SP),REAL(cB,SP))                

              ELSE IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX9) THEN
                
                CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                CALL do_mat79mat79add_numb_snglsngl(1,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,FaA,cA,KldB,KcolB,FaB,cB,KldC,KcolC,KdiagonalC,FaC)

              ELSE
                
                CALL do_mat79mat79add_numb_snglsngl(1,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,FaA,cA,KldB,KcolB,FaB,cB,KldC,KcolC,KldC,FaC)
                
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE DEFAULT
            PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE DEFAULT
        PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
        CALL sys_halt()
      END SELECT
      
      !-------------------------------------------------------------------------
      
    CASE (LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9INTL) ! A is interleave CSR matrix -
      
      ! Set size of interleave block
      isizeIntl=rmatrixB%NVAR
      IF (rmatrixB%cinterleavematrixFormat == LSYSSC_MATRIX1) isizeIntl=isizeIntl*isizeIntl

      SELECT CASE(rmatrixB%cmatrixFormat)

      CASE (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - - - - - - - -

        ! memory allocation?
        IF (bmemory) THEN
          rmatrixC%cmatrixFormat = rmatrixA%cmatrixFormat
          rmatrixC%cinterleavematrixFormat = rmatrixA%cinterleavematrixFormat
          rmatrixC%NA = rmatrixA%NA
          rmatrixC%NVAR = rmatrixA%NVAR
          CALL storage_new('lsyssc_matrixLinearComb','h_Da',INT(isizeIntl*rmatrixC%NA,I32),&
              rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        END IF

        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= rmatrixA%cmatrixFormat .OR. &
            rmatrixC%cinterleavematrixFormat /= rmatrixA%cinterleavematrixFormat) THEN
          PRINT *, 'lsyssc_matrixLinearComb: destination matrix has incompatible format'
          CALL sys_halt()
        END IF

        ! symbolic addition?
        IF (bsymb) THEN
          CALL lsyssc_duplicateMatrix(rmatrixA,rmatrixC,LSYSSC_DUP_COPY,LSYSSC_DUP_IGNORE)
        END IF

        ! numerical addition?
        IF (bnumb) THEN
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          SELECT CASE(rmatrixA%cdataType)
            
          CASE (ST_DOUBLE)
            
            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7INTL) THEN
                  CALL lsyssc_getbase_Kcol(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_doubledouble(KldA,KcolA,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixA%NEQ,DaA,cA,DaB,cB,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                ELSE
                  CALL do_mat79matDadd_doubledouble(KldA,KcolA,rmatrixC%NVAR,1,&
                      rmatrixA%NEQ,DaA,cA,DaB,cB,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                END IF

              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_doubledouble(KldA,KcolA,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixA%NEQ,DaA,cA,DaB,cB,DaC,KldC,KcolC)
                ELSE
                  CALL do_mat79matDadd_doubledouble(KldA,KcolA,rmatrixC%NVAR,1,&
                      rmatrixA%NEQ,DaA,cA,DaB,cB,DaC,KldC,KcolC)
                END IF
              END IF

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7INTL) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_doublesingle(KldA,KcolA,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixA%NEQ,DaA,cA,FaB,cB,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                ELSE
                  CALL do_mat79matDadd_doublesingle(KldA,KcolA,rmatrixC%NVAR,1,&
                      rmatrixA%NEQ,DaA,cA,FaB,cB,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                END IF
                
              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_doublesingle(KldA,KcolA,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixA%NEQ,DaA,cA,FaB,cB,DaC,KldC,KcolC)
                ELSE
                  CALL do_mat79matDadd_doublesingle(KldA,KcolA,rmatrixC%NVAR,1,&
                      rmatrixA%NEQ,DaA,cA,FaB,cB,DaC,KldC,KcolC)
                END IF
              END IF

            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)
            
            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7INTL) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_singledouble(KldA,KcolA,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixA%NEQ,FaA,cA,DaB,cB,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                ELSE
                  CALL do_mat79matDadd_singledouble(KldA,KcolA,rmatrixC%NVAR,1,&
                      rmatrixA%NEQ,FaA,cA,DaB,cB,DaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                END IF

              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_singledouble(KldA,KcolA,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixA%NEQ,FaA,cA,DaB,cB,DaC,KldC,KcolC)
                ELSE
                  CALL do_mat79matDadd_singledouble(KldA,KcolA,rmatrixC%NVAR,1,&
                      rmatrixA%NEQ,FaA,cA,DaB,cB,DaC,KldC,KcolC)
                END IF
              END IF

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)
              CALL lsyssc_getbase_Kld(rmatrixA,KldA)
              CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
              IF (bfast) THEN
                IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX7INTL) THEN
                  CALL lsyssc_getbase_Kld(rmatrixC,KdiagonalC)
                ELSE
                  CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                END IF
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_singlesingle(KldA,KcolA,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixA%NEQ,FaA,cA,FaB,cB,FaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                ELSE
                  CALL do_mat79matDadd_singlesingle(KldA,KcolA,rmatrixC%NVAR,1,&
                      rmatrixA%NEQ,FaA,cA,FaB,cB,FaC,Kdiag3=KdiagonalC,na=rmatrixC%NA)
                END IF

              ELSE
                CALL lsyssc_getbase_Kld(rmatrixC,KldC)
                CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
                IF (rmatrixC%cinterleavematrixFormat == LSYSSC_MATRIX1) THEN
                  CALL do_mat79matDadd_singlesingle(KldA,KcolA,rmatrixC%NVAR,rmatrixC%NVAR,&
                      rmatrixA%NEQ,FaA,cA,FaB,cB,FaC,KldC,KcolC)
                ELSE
                  CALL do_mat79matDadd_singlesingle(KldA,KcolA,rmatrixC%NVAR,1,&
                      rmatrixA%NEQ,FaA,cA,FaB,cB,FaC,KldC,KcolC)
                END IF
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE DEFAULT
            PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF
                
      CASE (LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9INTL) ! B is interleave CSR matrix-

        ! Check if interleave matrices are compatible
        IF (rmatrixA%NVAR /= rmatrixB% NVAR .OR. &
            rmatrixA%cinterleavematrixFormat /= rmatrixB%cinterleavematrixFormat) THEN
          PRINT *, 'lsyssc_matrixLinearComb: incompatible interleave matrices!'
          CALL sys_halt()
        END IF
        
        ! Set size of interleave block
        isizeIntl=rmatrixA%NVAR
        IF (rmatrixA%cinterleavematrixFormat == LSYSSC_MATRIX1) isizeIntl=isizeIntl*isizeIntl
        
        ! Set pointers
        CALL lsyssc_getbase_Kld(rmatrixA,KldA)
        CALL lsyssc_getbase_Kcol(rmatrixA,KcolA)
        CALL lsyssc_getbase_Kld(rmatrixB,KldB)
        CALL lsyssc_getbase_Kcol(rmatrixB,KcolB)
        
        ! memory allocation?
        IF (bmemory) THEN
          ! Duplicate structure of matrix A or B depending on which
          ! matrix is given in the encompassing matrix format
          IF (rmatrixA%cmatrixFormat >= rmatrixB%cmatrixFormat) THEN
            CALL lsyssc_duplicateMatrix(rmatrixA,rmatrixC,LSYSSC_DUP_EMPTY,LSYSSC_DUP_REMOVE)
          ELSE
            CALL lsyssc_duplicateMatrix(rmatrixB,rmatrixC,LSYSSC_DUP_EMPTY,LSYSSC_DUP_REMOVE)
          END IF
          
          ! Set auxiliary pointers
          CALL storage_new('lsyssc_matrixLinearComb','Kaux',rmatrixB%NCOLS,ST_INT,h_Kaux,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int(h_Kaux,Kaux)
          
          ! Compute number of nonzero matrix entries: NA
          rmatrixC%NA=do_mat79mat79add_computeNA(rmatrixA%NEQ,rmatrixA%NCOLS,KldA,KcolA,KldB,KcolB,Kaux)
          CALL storage_new('lsyssc_matrixLinearComb','h_Da',INT(isizeIntl*rmatrixC%NA,I32),&
              rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
          CALL storage_realloc('lsyssc_matrixLinearComb',rmatrixC%NA,&
              rmatrixC%h_Kcol,ST_NEWBLOCK_NOINIT,.FALSE.)
        END IF

        ! Check if matrix is given in the correct format
        IF (rmatrixC%cmatrixFormat /= MAX(rmatrixA%cmatrixFormat,rmatrixB%cmatrixFormat) .OR. &
          rmatrixC%cinterleavematrixFormat /= rmatrixA%cinterleavematrixFormat .OR. &
          rmatrixC%NA /= rmatrixA%NVAR) THEN
          PRINT *, 'lsyssc_matrixLinearComb: destination matrix has incompatible format'
          CALL sys_halt()
        END IF

        ! symbolic addition?
        IF (bsymb) THEN
          
          ! Set pointers
          CALL lsyssc_getbase_Kld(rmatrixC,KldC)
          CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
          
          IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX9INTL) THEN
            CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
            CALL do_mat79mat79add_symb(rmatrixC%NEQ,rmatrixC%NCOLS,KldA,KcolA,&
                rmatrixA%cmatrixFormat,KldB,KcolB,rmatrixB%cmatrixFormat,KldC,KcolC,KdiagonalC)
          ELSE
            CALL do_mat79mat79add_symb(rmatrixC%NEQ,rmatrixC%NCOLS,KldA,KcolA,&
                rmatrixA%cmatrixFormat,KldB,KcolB,rmatrixB%cmatrixFormat,KldC,KcolC)
          END IF
        END IF

        ! numerical addition?
        IF (bnumb) THEN
          
          ! Set pointers
          CALL lsyssc_getbase_Kld(rmatrixC,KldC)
          CALL lsyssc_getbase_Kcol(rmatrixC,KcolC)
          
          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double   
          SELECT CASE(rmatrixA%cdataType)
            
          CASE (ST_DOUBLE)

            SELECT CASE(rmatrixB%cdataType)

            CASE (ST_DOUBLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)

              IF ((rmatrixA%cmatrixFormat == rmatrixB%cmatrixFormat) .AND. &
                  (rmatrixA%cmatrixFormat == rmatrixC%cmatrixFormat) .AND. (bfast)) THEN
                
                ! We rely on the user who tells us that we should assume the
                ! same structure for the matrices. In this case, we can
                ! directly call a BLAS routine to do the matrix combination.
                ! That's MUCH faster!
                ! Note that matrix B might coincide with matrix C, so we
                ! first check this to prevent unnecessary copy operations!
                
                IF (.NOT. ASSOCIATED (DaB,DaC)) THEN
                  CALL lalg_copyVectorDble (DaB,DaC)
                END IF
                
                CALL lalg_vectorLinearCombDble (DaA,DaC,cA,cB)                
                
              ELSE IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX9INTL) THEN
                
                CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                CALL do_mat79mat79add_numb_dbledble(isizeIntl,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,DaA,cA,KldB,KcolB,DaB,cB,KldC,KcolC,KdiagonalC,DaC)
                
              ELSE
                
                CALL do_mat79mat79add_numb_dbledble(isizeIntl,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,DaA,cA,KldB,KcolB,DaB,cB,KldC,KcolC,KldC,DaC)
                
              END IF

            CASE (ST_SINGLE)
              CALL lsyssc_getbase_double(rmatrixA,DaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              
              IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX9INTL) THEN
                
                CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                CALL do_mat79mat79add_numb_dblesngl(isizeIntl,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,DaA,cA,KldB,KcolB,FaB,cB,KldC,KcolC,KdiagonalC,DaC)
                
              ELSE
                
                CALL do_mat79mat79add_numb_dblesngl(isizeIntl,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,DaA,cA,KldB,KcolB,FaB,cB,KldC,KcolC,KldC,DaC)
                
              END IF
              
            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT

          CASE (ST_SINGLE)
            
            SELECT CASE(rmatrixB%cdataType)
              
            CASE (ST_DOUBLE) 
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_double(rmatrixB,DaB)
              CALL lsyssc_getbase_double(rmatrixC,DaC)
              
              IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX9INTL) THEN
                
                CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                CALL do_mat79mat79add_numb_sngldble(isizeIntl,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,FaA,cA,KldB,KcolB,DaB,cB,KldC,KcolC,KdiagonalC,DaC)

              ELSE
                
                CALL do_mat79mat79add_numb_sngldble(isizeIntl,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,FaA,cA,KldB,KcolB,DaB,cB,KldC,KcolC,KldC,DaC)
                
              END IF
              
            CASE (ST_SINGLE)
              CALL lsyssc_getbase_single(rmatrixA,FaA)
              CALL lsyssc_getbase_single(rmatrixB,FaB)
              CALL lsyssc_getbase_single(rmatrixC,FaC)

              IF ((rmatrixA%cmatrixFormat == rmatrixB%cmatrixFormat) .AND. &
                  (rmatrixA%cmatrixFormat == rmatrixC%cmatrixFormat) .AND. (bfast)) THEN
                
                ! We rely on the user who tells us that we should assume the
                ! same structure for the matrices. In this case, we can
                ! directly call a BLAS routine to do the matrix combination.
                ! That's MUCH faster!
                ! Note that matrix B might coincide with matrix C, so we
                ! first check this to prevent unnecessary copy operations!
                
                IF (.NOT. ASSOCIATED (DaB,DaC)) THEN
                  CALL lalg_copyVectorSngl (FaB,FaC)
                END IF
                
                CALL lalg_vectorLinearCombSngl (FaA,FaC,REAL(cA,SP),REAL(cB,SP))                
                
              ELSE IF (rmatrixC%cmatrixFormat == LSYSSC_MATRIX9INTL) THEN
                
                CALL lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
                CALL do_mat79mat79add_numb_snglsngl(isizeIntl,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,FaA,cA,KldB,KcolB,FaB,cB,KldC,KcolC,KdiagonalC,FaC)

              ELSE
                
                CALL do_mat79mat79add_numb_snglsngl(isizeIntl,rmatrixC%NEQ,rmatrixC%NCOLS,&
                    KldA,KcolA,FaA,cA,KldB,KcolB,FaB,cB,KldC,KcolC,KldC,FaC)
                
              END IF

            CASE DEFAULT
              PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
              CALL sys_halt()
            END SELECT
            
          CASE DEFAULT
            PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
            CALL sys_halt()
          END SELECT
        END IF

      CASE DEFAULT
        PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
        CALL sys_halt()
      END SELECT

      !-------------------------------------------------------------------------

    CASE DEFAULT
      PRINT *, 'lsyssc_matrixLinearComb: Unsupported data type!'
      CALL sys_halt()
    END SELECT

    ! Clear auxiliary vectors
    IF (h_Kaux /= ST_NOHANDLE) CALL storage_free(h_Kaux)

  CONTAINS

    ! Here, the real MM addition routines follow.

    !**************************************************************
    ! Format 1-D addition
    ! double precision matrix A (format 1)
    ! double precision matrix B (format D)
    ! double precision matrix C (format 1)

    SUBROUTINE do_mat1matDadd_doubledouble(neq,ncols,Da1,c1,Da2,c2,Da3)

      INTEGER(PREC_DOFIDX), INTENT(IN)              :: neq,ncols
      REAL(DP), INTENT(IN)                          :: c1,c2
      REAL(DP), DIMENSION(ncols,neq), INTENT(IN)    :: Da1
      REAL(DP), DIMENSION(neq), INTENT(IN)          :: Da2
      REAL(DP), DIMENSION(ncols,neq), INTENT(INOUT) :: Da3
      INTEGER(PREC_DOFIDX)                          :: ieq

      CALL DCOPY(INT(neq*ncols),Da1,1,Da3,1)
      CALL DSCAL(INT(neq*ncols),c1,Da3,1)
!$omp  parallel do&
!$omp& default(shared)&
!$omp& private(ieq)
      DO ieq=1,neq
        Da3(ieq,ieq)=Da3(ieq,ieq)+c2*Da2(ieq)
      END DO
!$omp  end parallel do

!!$      REAL(DP), DIMENSION(m*n), INTENT(IN)    :: Da1
!!$      REAL(DP), DIMENSION(m*n), INTENT(INOUT) :: Da3
!!$      CALL lalg_copyVectorDble(Da1,Da3)
!!$      CALL lalg_scaleVectorDble(Da3,c1)
!!$      DO i=0,n-1
!!$        Da3(i*m+i+1)=Da3(i*m+i+1)+c2*Da2(i+1)
!!$      END DO
    END SUBROUTINE do_mat1matDadd_doubledouble

    !**************************************************************
    ! Format 1-D addition
    ! single precision matrix A (format 1)
    ! double precision matrix B (format D)
    ! double precision matrix C (format 1)

    SUBROUTINE do_mat1matDadd_singledouble(neq,ncols,Fa1,c1,Da2,c2,Da3)

      INTEGER(PREC_DOFIDX), INTENT(IN)              :: neq,ncols
      REAL(DP), INTENT(IN)                          :: c1,c2
      REAL(SP), DIMENSION(ncols,neq), INTENT(IN)    :: Fa1
      REAL(DP), DIMENSION(neq), INTENT(IN)          :: Da2
      REAL(DP), DIMENSION(ncols,neq), INTENT(INOUT) :: Da3
      INTEGER(PREC_DOFIDX)                          :: ieq,icols

!$omp  parallel do &
!$omp& default(shared) &
!$omp& private(ieq,icols)
      DO ieq=1,neq
        DO icols=1,ncols
          Da3(icols,ieq) = c1*Fa1(icols,ieq)
        END DO
      END DO
!$omp  end parallel do

!$omp  parallel do&
!$omp& default(shared)&
!$omp& private(ieq)
      DO ieq=1,neq
        Da3(ieq,ieq)=Da3(ieq,ieq)+c2*Da2(ieq)
      END DO
!$omp  end parallel do

!!$      REAL(SP), DIMENSION(m*n), INTENT(IN)    :: Fa1
!!$      REAL(DP), DIMENSION(m*n), INTENT(INOUT) :: Da3
!!$      CALL lalg_scaleVectorDble(Da3,c1)
!!$      DO i=0,n-1
!!$        Da3(i*m+i+1)=Da3(i*m+i+1)+c2*Da2(i+1)
!!$      END DO
    END SUBROUTINE do_mat1matDadd_singledouble

    !**************************************************************
    ! Format 1-D addition
    ! double precision matrix A (format 1)
    ! single precision matrix B (format D)
    ! double precision matrix C (format 1)

    SUBROUTINE do_mat1matDadd_doublesingle(neq,ncols,Da1,c1,Fa2,c2,Da3)

      INTEGER(PREC_DOFIDX), INTENT(IN)              :: neq,ncols
      REAL(DP), INTENT(IN)                          :: c1,c2
      REAL(DP), DIMENSION(ncols,neq), INTENT(IN)    :: Da1
      REAL(SP), DIMENSION(neq), INTENT(IN)          :: Fa2
      REAL(DP), DIMENSION(ncols,neq), INTENT(INOUT) :: Da3
      INTEGER(PREC_DOFIDX) :: ieq

      CALL DCOPY(INT(neq*ncols),Da1,1,Da3,1)
      CALL DSCAL(INT(neq*ncols),c1,Da3,1)
!$omp  parallel do&
!$omp& default(shared)&
!$omp& private(ieq)
      DO ieq=1,neq
        Da3(ieq,ieq)=Da3(ieq,ieq)+c2*Fa2(ieq)
      END DO
!$omp  end parallel do

!!$      REAL(DP), DIMENSION(m*n), INTENT(IN)    :: Da1
!!$      REAL(DP), DIMENSION(m*n), INTENT(INOUT) :: Da3
!!$      CALL lalg_copyVectorDble(Da1,Da3)
!!$      CALL lalg_scaleVectorDble(Da3,c1)
!!$      DO i=0,n-1
!!$        Da3(i*m+i+1)=Da3(i*m+i+1)+c2*Fa2(i+1)
!!$      END DO
    END SUBROUTINE do_mat1matDadd_doublesingle

    !**************************************************************
    ! Format 1-D addition
    ! single precision matrix A (format 1)
    ! single precision matrix B (format D)
    ! single precision matrix C (format 1)

    SUBROUTINE do_mat1matDadd_singlesingle(neq,ncols,Fa1,c1,Fa2,c2,Fa3)

      INTEGER(PREC_DOFIDX), INTENT(IN) :: neq,ncols
      REAL(DP), INTENT(IN) :: c1,c2
      REAL(SP), DIMENSION(ncols,neq), INTENT(IN)    :: Fa1
      REAL(SP), DIMENSION(neq), INTENT(IN)      :: Fa2
      REAL(SP), DIMENSION(ncols,neq), INTENT(INOUT) :: Fa3
      INTEGER(PREC_DOFIDX) :: ieq

      CALL SCOPY(INT(neq*ncols),Fa1,1,Fa3,1)
      CALL SSCAL(INT(neq*ncols),REAL(c1,SP),Fa3,1)
!$omp  parallel do&
!$omp& default(shared)&
!$omp& private(ieq)
      DO ieq=1,neq
        Fa3(ieq,ieq)=Fa3(ieq,ieq)+c2*Fa2(ieq)
      END DO
!$omp  end parallel do

!!$      REAL(SP), DIMENSION(m*n), INTENT(IN)    :: Fa1
!!$      REAL(SP), DIMENSION(m*n), INTENT(INOUT) :: Fa3
!!$      CALL lalg_copyVectorSngl(Fa1,Fa3)
!!$      CALL lalg_scaleVectorSngl(Fa3,REAL(c1,SP))
!!$      DO i=0,n-1
!!$        Fa3(i*m+i+1)=Fa3(i*m+i+1)+c2*Fa2(i+1)
!!$      END DO
    END SUBROUTINE do_mat1matDadd_singlesingle

    !**************************************************************
    ! Format 1-7/9 addition
    ! double precision matrix A (format 1)
    ! double precision matrix B (format 7 or format 9)
    ! double precision matrix C (format 1)

    SUBROUTINE do_mat1mat79add_doubledouble(neq,ncols,Da1,c1,Kld2,Kcol2,Da2,c2,Da3)

      ! REMARK: The matrix A (format 1) is stored row-wise. Hence,
      ! this subroutine handles the transposed of matrix A. Therefore,
      ! the row and column indices IEQ and ICOL are swapped when the
      ! contribution of matrix B is applied to matrix C.

      INTEGER(PREC_DOFIDX), INTENT(IN)               :: neq,ncols
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld2
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol2
      REAL(DP), INTENT(IN)                           :: c1,c2
      REAL(DP), DIMENSION(ncols,neq), INTENT(IN)     :: Da1
      REAL(DP), DIMENSION(:), INTENT(IN)             :: Da2
      REAL(DP), DIMENSION(ncols,neq), INTENT(INOUT)  :: Da3
      INTEGER(PREC_DOFIDX) :: ild,ieq,icol

      CALL DCOPY(INT(neq*ncols),Da1,1,Da3,1)
      CALL DSCAL(INT(neq*ncols),c1,Da3,1)
      
      DO ieq=1,neq
        DO ild=Kld2(ieq),Kld2(ieq+1)-1
          icol=Kcol2(ild)
          Da3(icol,ieq)=Da3(icol,ieq)+c2*Da2(ild)
        END DO
      END DO

!!$      REAL(DP), DIMENSION(m*n), INTENT(IN)    :: Da1
!!$      REAL(DP), DIMENSION(m*n), INTENT(INOUT) :: Da3
!!$      CALL lalg_copyVectorDble(Da1,Da3)
!!$      CALL lalg_scaleVectorDble(Da3,c1)
!!$      
!!$      DO i=1,n
!!$        DO ild=Kld2(i),Kld2(i+1)-1
!!$          j=Kcol2(ild)-1
!!$          Da3(j*m+i)=Da3(j*m+i)+c2*Da2(ild)
!!$        END DO
!!$      END DO
    END SUBROUTINE do_mat1mat79add_doubledouble

    !**************************************************************
    ! Format 1-7/9 addition
    ! single precision matrix A (format 1)
    ! double precision matrix B (format 7 or format 9)
    ! double precision matrix C (format 1)

    SUBROUTINE do_mat1mat79add_singledouble(neq,ncols,Fa1,c1,Kld2,Kcol2,Da2,c2,Da3)

      ! REMARK: The matrix A (format 1) is stored row-wise. Hence,
      ! this subroutine handles the transposed of matrix A. Therefore,
      ! the row and column indices IEQ and ICOL are swapped when the
      ! contribution of matrix B is applied to matrix C.

      INTEGER(PREC_DOFIDX), INTENT(IN)               :: neq,ncols
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld2
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol2
      REAL(DP), INTENT(IN)                           :: c1,c2
      REAL(SP), DIMENSION(ncols,neq), INTENT(IN)     :: Fa1
      REAL(DP), DIMENSION(:), INTENT(IN)             :: Da2
      REAL(DP), DIMENSION(ncols,neq), INTENT(INOUT)  :: Da3
      INTEGER(PREC_DOFIDX) :: ild,ieq,icol

!$omp  parallel do &
!$omp& default(shared) &
!$omp& private(ieq,icol)
      DO ieq=1,neq
        DO icol=1,ncols
          Da3(icol,ieq)=c1*Fa1(icol,ieq)
        END DO
      END DO
!$omp  end parallel do
      
      DO ieq=1,neq
        DO ild=Kld2(ieq),Kld2(ieq+1)-1
          icol=Kcol2(ild)
          Da3(icol,ieq)=Da3(icol,ieq)+c2*Da2(ild)
        END DO
      END DO
    END SUBROUTINE do_mat1mat79add_singledouble

    !**************************************************************
    ! Format 1-7/9 addition
    ! double precision matrix A (format 1)
    ! single precision matrix B (format 7 or format 9)
    ! double precision matrix C (format 1)

    SUBROUTINE do_mat1mat79add_doublesingle(neq,ncols,Da1,c1,Kld2,Kcol2,Fa2,c2,Da3)

      ! REMARK: The matrix A (format 1) is stored row-wise. Hence,
      ! this subroutine handles the transposed of matrix A. Therefore,
      ! the row and column indices IEQ and ICOL are swapped when the
      ! contribution of matrix B is applied to matrix C.

      INTEGER(PREC_DOFIDX), INTENT(IN)               :: neq,ncols
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld2
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol2
      REAL(DP), INTENT(IN)                           :: c1,c2
      REAL(DP), DIMENSION(ncols,neq), INTENT(IN)     :: Da1
      REAL(SP), DIMENSION(:), INTENT(IN)             :: Fa2
      REAL(DP), DIMENSION(ncols,neq), INTENT(INOUT)  :: Da3
      INTEGER(PREC_DOFIDX) :: ild,ieq,icol

      CALL DCOPY(INT(neq*ncols),Da1,1,Da3,1)
      CALL DSCAL(INT(neq*ncols),c1,Da3,1)

      DO ieq=1,neq
        DO ild=Kld2(ieq),Kld2(ieq+1)-1
          icol=Kcol2(ild)
          Da3(icol,ieq)=Da3(icol,ieq)+c2*Fa2(ild)
        END DO
      END DO
      
!!$      REAL(DP), DIMENSION(m*n), INTENT(IN)    :: Da1
!!$      REAL(DP), DIMENSION(m*n), INTENT(INOUT) :: Da3
!!$      CALL lalg_copyVectorDble(Da1,Da3)
!!$      CALL lalg_scaleVectorDble(Da3,c1)
!!$
!!$      DO i=1,n
!!$        DO ild=Kld2(i),Kld2(i+1)-1
!!$          j=Kcol2(ild)
!!$          Da3(j*m+i)=Da3(j*m+i)+c2*Fa2(ild)
!!$        END DO
!!$      END DO
    END SUBROUTINE do_mat1mat79add_doublesingle

    !**************************************************************
    ! Format 1-7/9 addition
    ! single precision matrix A (format 1)
    ! single precision matrix B (format 7 or format 9)
    ! single precision matrix C (format 1)

    SUBROUTINE do_mat1mat79add_singlesingle(neq,ncols,Fa1,c1,Kld2,Kcol2,Fa2,c2,Fa3)

      ! REMARK: The matrix A (format 1) is stored row-wise. Hence,
      ! this subroutine handles the transposed of matrix A. Therefore,
      ! the row and column indices IEQ and ICOL are swapped when the
      ! contribution of matrix B is applied to matrix C.

      INTEGER(PREC_DOFIDX), INTENT(IN)               :: neq,ncols
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld2
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol2
      REAL(DP), INTENT(IN)                           :: c1,c2
      REAL(SP), DIMENSION(ncols,neq), INTENT(IN)     :: Fa1
      REAL(SP), DIMENSION(:), INTENT(IN)             :: Fa2
      REAL(SP), DIMENSION(ncols,neq), INTENT(INOUT)  :: Fa3
      INTEGER(PREC_DOFIDX) :: ild,ieq,icol
      
      CALL SCOPY(INT(neq*ncols),Fa1,1,Fa3,1)
      CALL SSCAL(INT(neq*ncols),REAL(c1,SP),Fa3,1)
      
      DO ieq=1,neq
        DO ild=Kld2(ieq),Kld2(ieq+1)-1
          icol=Kcol2(ild)
          Fa3(icol,ieq)=Fa3(icol,ieq)+c2*Fa2(ild)
        END DO
      END DO
      
!!$      REAL(SP), DIMENSION(m*n), INTENT(IN)    :: Fa1
!!$      REAL(SP), DIMENSION(m*n), INTENT(INOUT) :: Fa3
!!$      CALL lalg_copyVectorSngl(Fa1,Fa3)
!!$      CALL lalg_scaleVectorSngl(Fa3,REAL(c1,SP))
!!$      
!!$      DO i=1,n
!!$        DO ild=Kld2(i),Kld2(i+1)-1
!!$          j=Kcol2(ild)
!!$          Fa3(j*m+i)=Fa3(j*m+i)+c2*Fa2(ild)
!!$        END DO
!!$      END DO
    END SUBROUTINE do_mat1mat79add_singlesingle

    !**************************************************************
    ! Format 7/9-D addition
    ! double precision matrix A (format 7 or format 9, interleave possible)
    ! double precision matrix B (format D)
    ! double precision matrix C (format 7 or format 9, interleave possible)

    SUBROUTINE do_mat79matDadd_doubledouble(Kld1,Kcol1,nvar,mvar,neq,&
        Da1,c1,Da2,c2,Da3,Kld3,Kcol3,Kdiag3,na)
      
      INTEGER(PREC_DOFIDX), INTENT(IN)                         :: neq
      INTEGER, INTENT(IN)                                      :: nvar,mvar
      INTEGER(PREC_MATIDX), INTENT(IN), OPTIONAL               :: na
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)           :: Kld1
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)           :: Kcol1
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kld3,Kdiag3
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kcol3
      REAL(DP), INTENT(IN)                                     :: c1,c2
      REAL(DP), DIMENSION(nvar,mvar,*), INTENT(IN)             :: Da1
      REAL(DP), DIMENSION(:), INTENT(IN)                       :: Da2
      REAL(DP), DIMENSION(nvar,mvar,*), INTENT(INOUT)          :: Da3
      
      INTEGER(PREC_DOFIDX) :: ieq
      INTEGER(PREC_VECIDX) :: ild1,ild3,ildend3
      INTEGER(PREC_MATIDX) :: icol1,icol3
      INTEGER :: ivar

      IF (PRESENT(Kld3) .AND. PRESENT(Kcol3)) THEN
      
        ! The structure of matrix A may be only a subset of C, so that
        ! each entry must be visited separately
        
        IF (nvar /= mvar) THEN

          ! Loop over all rows
          DO ieq=1,neq
            ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
            
            ! Loop over all columns of current row
            DO ild1=Kld1(ieq),Kld1(ieq+1)-1
              icol1=Kcol1(ild1)
              
              ! Skip positions in resulting matrix if required
              DO ild3=ild3,ildend3
                icol3=Kcol3(ild3)
                IF (icol3 == icol1) EXIT
                Da3(:,1,ild3)=0._DP
              END DO
              
              Da3(:,1,ild3)=c1*Da1(:,1,ild1)
              ! Diagonal entry?
              IF (icol3 == ieq) Da3(:,1,ild3)=Da3(:,1,ild3)+c2*Da2(ieq)
              ild3=ild3+1
            END DO
            Da3(:,1,ild3:ildend3)=0._DP
          END DO
        
        ELSE
          
          ! Loop over all rows
          DO ieq=1,neq
            ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
            
            ! Loop over all columns of current row
            DO ild1=Kld1(ieq),Kld1(ieq+1)-1
              icol1=Kcol1(ild1)
              
              ! Skip positions in resulting matrix if required
              DO ild3=ild3,ildend3
                icol3=Kcol3(ild3)
                IF (icol3 == icol1) EXIT
                Da3(:,:,ild3)=0._DP
              END DO
              
              Da3(:,:,ild3)=c1*Da1(:,:,ild1)
              ! Diagonal entry?
              IF (icol3 == ieq) THEN
                DO ivar=1,nvar
                  Da3(ivar,ivar,ild3)=Da3(ivar,ivar,ild3)+c2*Da2(ieq)
                END DO
              END IF
              ild3=ild3+1
            END DO
            Da3(:,:,ild3:ildend3)=0._DP
          END DO

        END IF

      ELSEIF (PRESENT(Kdiag3) .AND. PRESENT(na)) THEN
        
        ! Structure of matrices A and C is identical

        CALL DCOPY(INT(nvar*mvar*na),Da1,1,Da3,1)
        CALL DSCAL(INT(nvar*mvar*na),c1,Da3,1)

        IF (nvar /= mvar) THEN
!$omp  parallel do &
!$omp& default(shared) &
!$omp& private(ieq,ild3)        
          DO ieq=1,neq
            ild3=Kdiag3(ieq)
            Da3(:,1,ild3)=Da3(:,1,ild3)+c2*Da2(ieq)
          END DO
!$omp  end parallel do
        ELSE
!$omp  parallel do &
!$omp& default(shared) &
!$omp& private(ivar,ieq,ild3)        
          DO ieq=1,neq
            ild3=Kdiag3(ieq)
            DO ivar=1,nvar
              Da3(ivar,ivar,ild3)=Da3(ivar,ivar,ild3)+c2*Da2(ieq)
            END DO
          END DO
!$omp  end parallel do
        END IF

      ELSE
        PRINT *, "do_mat79matDadd_doubledouble: either Kld,Kcol or Kdiag must be present."
        CALL sys_halt()        
      END IF
    END SUBROUTINE do_mat79matDadd_doubledouble

    !**************************************************************
    ! Format 7/9-D addition
    ! single precision matrix A (format 7 or format 9, interleave possible)
    ! double precision matrix B (format D)
    ! double precision matrix C (format 7 or format 9, interleave possible)

    SUBROUTINE do_mat79matDadd_singledouble(Kld1,Kcol1,nvar,mvar,neq,&
        Fa1,c1,Da2,c2,Da3,Kld3,Kcol3,Kdiag3,na)
      
      INTEGER(PREC_DOFIDX), INTENT(IN)                         :: neq
      INTEGER, INTENT(IN)                                      :: nvar,mvar
      INTEGER(PREC_MATIDX), INTENT(IN), OPTIONAL               :: na
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)           :: Kld1
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)           :: Kcol1
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kld3,Kdiag3
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kcol3
      REAL(DP), INTENT(IN)                                     :: c1,c2
      REAL(SP), DIMENSION(nvar,mvar,*), INTENT(IN)             :: Fa1
      REAL(DP), DIMENSION(:), INTENT(IN)                       :: Da2
      REAL(DP), DIMENSION(nvar,mvar,*), INTENT(INOUT)          :: Da3

      INTEGER(PREC_DOFIDX) :: ieq
      INTEGER(PREC_VECIDX) :: ild1,ild3,ildend3
      INTEGER(PREC_MATIDX) :: icol1,icol3,ia
      INTEGER :: ivar
      
      IF (PRESENT(Kld3) .AND. PRESENT(Kcol3)) THEN
        
        ! The structure of matrix A may be only a subset of C, so that
        ! each entry must be visited separately
        
        IF (nvar /= mvar) THEN
          
          ! Loop over all rows
          DO ieq=1,neq
            ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
            
            ! Loop over all columns of current row
            DO ild1=Kld1(ieq),Kld1(ieq+1)-1
              icol1=Kcol1(ild1)
              
              ! Skip positions in resulting matrix if required
              DO ild3=ild3,ildend3
                icol3=Kcol3(ild3)
                IF (icol3 == icol1) EXIT
                Da3(:,1,ild3)=0._DP
              END DO
              
              Da3(:,1,ild3)=c1*Fa1(:,1,ild1)
              ! Diagonal entry?
              IF (icol3 == ieq) Da3(:,1,ild3)=Da3(:,1,ild3)+c2*Da2(ieq)
              ild3=ild3+1
            END DO
            Da3(:,1,ild3:ildend3)=0._DP
          END DO
          
        ELSE

          ! Loop over all rows
          DO ieq=1,neq
            ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
            
            ! Loop over all columns of current row
            DO ild1=Kld1(ieq),Kld1(ieq+1)-1
              icol1=Kcol1(ild1)
              
              ! Skip positions in resulting matrix if required
              DO ild3=ild3,ildend3
                icol3=Kcol3(ild3)
                IF (icol3 == icol1) EXIT
                Da3(:,:,ild3)=0._DP
              END DO
              
              Da3(:,:,ild3)=c1*Fa1(:,:,ild1)
              ! Diagonal entry?
              IF (icol3 == ieq) THEN
                DO ivar=1,nvar
                  Da3(ivar,ivar,ild3)=Da3(ivar,ivar,ild3)+c2*Da2(ieq)
                END DO
              END IF
              ild3=ild3+1
            END DO
            Da3(:,:,ild3:ildend3)=0._DP
          END DO

        END IF

      ELSEIF (PRESENT(Kdiag3) .AND. PRESENT(na)) THEN
        
        ! Structure of matrices A and C is identical

!$omp  parallel do &
!$omp& default(shared) &
!$omp& private(ia)
        DO ia=1,na
          Da3(:,:,ia)=c1*Fa1(:,:,ia)
        END DO
!$omp  end parallel do
        
        IF (nvar /= mvar) THEN
!$omp  parallel do &
!$omp& default(shared) &
!$omp& private(ieq,ild3)
          DO ieq=1,neq
            ild3=Kdiag3(ieq)
            Da3(:,1,ild3)=Da3(:,1,ild3)+c2*Da2(ieq)
          END DO
!$omp  end parallel do
        ELSE
!$omp  parallel do &
!$omp& default(shared) &
!$omp& private(ivar,ieq,ild3)
          DO ieq=1,neq
            ild3=Kdiag3(ieq)
            DO ivar=1,nvar
              Da3(ivar,ivar,ild3)=Da3(ivar,ivar,ild3)+c2*Da2(ieq)
            END DO
          END DO
!$omp  end parallel do
        END IF
        
      ELSE
        PRINT *, "do_mat79matDadd_singledouble: either Kld,Kcol or Kdiag must be present."
        CALL sys_halt()        
      END IF
    END SUBROUTINE do_mat79matDadd_singledouble

    !**************************************************************
    ! Format 7/9-D addition
    ! double precision matrix A (format 7 or format 9, interleave possible)
    ! single precision matrix B (format D)
    ! double precision matrix C (format 7 or format 9, interleave possible)

    SUBROUTINE do_mat79matDadd_doublesingle(Kld1,Kcol1,nvar,mvar,neq,&
        Da1,c1,Fa2,c2,Da3,Kld3,Kcol3,Kdiag3,na)

      INTEGER(PREC_DOFIDX), INTENT(IN)                         :: neq
      INTEGER, INTENT(IN)                                      :: nvar,mvar
      INTEGER(PREC_MATIDX), INTENT(IN), OPTIONAL               :: na
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)           :: Kld1
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)           :: Kcol1
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kld3,Kdiag3
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kcol3
      REAL(DP), INTENT(IN)                                     :: c1,c2
      REAL(DP), DIMENSION(nvar,mvar,*), INTENT(IN)             :: Da1
      REAL(SP), DIMENSION(:), INTENT(IN)                       :: Fa2
      REAL(DP), DIMENSION(nvar,mvar,*), INTENT(INOUT)          :: Da3

      INTEGER(PREC_DOFIDX) :: ieq
      INTEGER(PREC_VECIDX) :: ild1,ild3,ildend3
      INTEGER(PREC_MATIDX) :: icol1,icol3
      INTEGER :: ivar
      
      IF (PRESENT(Kld3) .AND. PRESENT(Kcol3)) THEN
        
        ! The structure of matrix A may be only a subset of C, so that
        ! each entry must be visited separately
        
        IF (nvar /= mvar) THEN
          
          ! Loop over all rows
          DO ieq=1,neq
            ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
            
            ! Loop over all columns of current row
            DO ild1=Kld1(ieq),Kld1(ieq+1)-1
              icol1=Kcol1(ild1)
              
              ! Skip positions in resulting matrix if required
              DO ild3=ild3,ildend3
                icol3=Kcol3(ild3)
                IF (icol3 == icol1) EXIT
                Da3(:,1,ild3)=0._DP
              END DO
              
              Da3(:,1,ild3)=c1*Da1(:,1,ild1)
              ! Diagonal entry?
              IF (icol3 == ieq) Da3(:,1,ild3)=Da3(:,1,ild3)+c2*Fa2(ieq)
              ild3=ild3+1
            END DO
            Da3(:,1,ild3:ildend3)=0._DP
          END DO
        
        ELSE

          ! Loop over all rows
          DO ieq=1,neq
            ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
            
            ! Loop over all columns of current row
            DO ild1=Kld1(ieq),Kld1(ieq+1)-1
              icol1=Kcol1(ild1)
              
              ! Skip positions in resulting matrix if required
              DO ild3=ild3,ildend3
                icol3=Kcol3(ild3)
                IF (icol3 == icol1) EXIT
                Da3(:,:,ild3)=0._DP
              END DO
              
              Da3(:,:,ild3)=c1*Da1(:,:,ild1)
              ! Diagonal entry?
              IF (icol3 == ieq) THEN
                DO ivar=1,nvar
                  Da3(ivar,ivar,ild3)=Da3(ivar,ivar,ild3)+c2*Fa2(ieq)
                END DO
              END IF
              ild3=ild3+1
            END DO
            Da3(:,:,ild3:ildend3)=0._DP
          END DO
          
        END IF
        
      ELSEIF (PRESENT(Kdiag3) .AND. PRESENT(na)) THEN

        ! Structure of matrices A and C is identical
        
        CALL DCOPY(INT(nvar*mvar*na),Da1,1,Da3,1)
        CALL DSCAL(INT(nvar*mvar*na),c1,Da3,1)

        IF (nvar /= mvar) THEN
!$omp  parallel do &
!$omp& default(shared) &
!$omp& private(ivar,ieq,ild3)
          DO ieq=1,neq
            ild3=Kdiag3(ieq)
            DO ivar=1,nvar
              Da3(ivar,1,ild3)=Da3(ivar,1,ild3)+c2*Fa2(ieq)
            END DO
          END DO
!$omp  end parallel do
        ELSE
!$omp  parallel do &
!$omp& default(shared) &
!$omp& private(ivar,ieq,ild3)
          DO ieq=1,neq
            ild3=Kdiag3(ieq)
            DO ivar=1,nvar
              Da3(ivar,ivar,ild3)=Da3(ivar,ivar,ild3)+c2*Fa2(ieq)
            END DO
          END DO
!$omp  end parallel do
        END IF
     
      ELSE
        PRINT *, "do_mat79matDadd_doublesingle: either Kld,Kcol or Kdiag must be present."
        CALL sys_halt()        
      END IF
    END SUBROUTINE do_mat79matDadd_doublesingle

    !**************************************************************
    ! Format 7/9-D addition
    ! single precision matrix A (format 7 or format 9, interleave possible)
    ! single precision matrix B (format D)
    ! single precision matrix C (format 7 or format 9, interleave possible)

    SUBROUTINE do_mat79matDadd_singlesingle(Kld1,Kcol1,nvar,mvar,neq,&
        Fa1,c1,Fa2,c2,Fa3,Kld3,Kcol3,Kdiag3,na)

      INTEGER(PREC_DOFIDX), INTENT(IN)                         :: neq
      INTEGER, INTENT(IN)                                      :: nvar,mvar
      INTEGER(PREC_MATIDX), INTENT(IN), OPTIONAL               :: na
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN)           :: Kld1
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN)           :: Kcol1
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kld3,Kdiag3
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Kcol3
      REAL(DP), INTENT(IN)                                     :: c1,c2
      REAL(SP), DIMENSION(nvar,mvar,*), INTENT(IN)             :: Fa1
      REAL(SP), DIMENSION(:), INTENT(IN)                       :: Fa2
      REAL(SP), DIMENSION(nvar,mvar,*), INTENT(INOUT)          :: Fa3

      INTEGER(PREC_DOFIDX) :: ieq
      INTEGER(PREC_VECIDX) :: ild1,ild3,ildend3
      INTEGER(PREC_MATIDX) :: icol1,icol3
      INTEGER :: ivar
      
      IF (PRESENT(Kld3) .AND. PRESENT(Kcol3)) THEN

        ! The structure of matrix A may be only a subset of C, so that
        ! each entry must be visited separately
        
        IF (nvar /= mvar) THEN

          ! Loop over all rows
          DO ieq=1,neq
            ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
            
            ! Loop over all columns of current row
            DO ild1=Kld1(ieq),Kld1(ieq+1)-1
              icol1=Kcol1(ild1)
              
              ! Skip positions in resulting matrix if required
              DO ild3=ild3,ildend3
                icol3=Kcol3(ild3)
                IF (icol3 == icol1) EXIT
                Fa3(:,1,ild3)=0._SP
              END DO
              
              Fa3(:,1,ild3)=c1*Fa1(:,1,ild1)
              ! Diagonal entry?
              IF (icol3 == ieq) Fa3(:,1,ild3)=Fa3(:,1,ild3)+c2*Fa2(ieq)
              ild3=ild3+1
            END DO
            Fa3(:,1,ild3:ildend3)=0._SP
          END DO

        ELSE

          ! Loop over all rows
          DO ieq=1,neq
            ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
            
            ! Loop over all columns of current row
            DO ild1=Kld1(ieq),Kld1(ieq+1)-1
              icol1=Kcol1(ild1)
              
              ! Skip positions in resulting matrix if required
              DO ild3=ild3,ildend3
                icol3=Kcol3(ild3)
                IF (icol3 == icol1) EXIT
                Fa3(:,:,ild3)=0._SP
              END DO
              
              Fa3(:,:,ild3)=c1*Fa1(:,:,ild1)
              ! Diagonal entry?
              IF (icol3 == ieq) THEN
                DO ivar=1,nvar
                  Fa3(ivar,ivar,ild3)=Fa3(ivar,ivar,ild3)+c2*Fa2(ieq)
                END DO
              END IF
              ild3=ild3+1
            END DO
            Fa3(:,1,ild3:ildend3)=0._SP
          END DO

        END IF
        
      ELSEIF (PRESENT(Kdiag3) .AND. PRESENT(na)) THEN

        ! Structure of matrices A and C is identical
        
        CALL SCOPY(INT(nvar*mvar*na),Fa1,1,Fa3,1)
        CALL SSCAL(INT(nvar*mvar*na),REAL(c1,SP),Fa3,1)

        IF (nvar /= mvar) THEN
!$omp  parallel do &
!$omp& default(shared) &
!$omp& private(ieq,ild3)      
          DO ieq=1,neq
            ild3=Kdiag3(ieq)
            Fa3(:,1,ild3)=Fa3(:,1,ild3)+c2*Fa2(ieq)
          END DO
!$omp  end parallel do
        ELSE
!$omp  parallel do &
!$omp& default(shared) &
!$omp& private(ivar,ieq,ild3)      
          DO ieq=1,neq
            ild3=Kdiag3(ieq)
            DO ivar=1,nvar
              Fa3(ivar,ivar,ild3)=Fa3(ivar,ivar,ild3)+c2*Fa2(ieq)
            END DO
          END DO
!$omp  end parallel do
        END IF

      ELSE
        PRINT *, "do_mat79matDadd_singlesingle: either Kld,Kcol or Kdiag must be present."
        CALL sys_halt()        
      END IF
    END SUBROUTINE do_mat79matDadd_singlesingle
    
    !**************************************************************
    ! Format 7/9-7/9 addition
    ! Compute the number of nonzero matrix entries of C:=A + B
    ! 
    ! Remark: This subroutine is a modified version of the APLBDG
    !         subroutine taken from the SPARSEKIT library written
    !         by Youcef Saad.

    FUNCTION do_mat79mat79add_computeNA(neq,ncols,KldA,KcolA,KldB&
        &,KcolB,Kaux) RESULT(NA)

      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA,KldB
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolB,KcolA
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: Kaux
      INTEGER(PREC_DOFIDX), INTENT(IN) :: neq,ncols
      INTEGER :: NA
      
      INTEGER :: ieq,jeq,ild,icol,idg,ndg,last

      ! Initialization
      Kaux=0; NA=0
      
      DO ieq=1,neq

        ! For each row of matrix A
        ndg=0

        ! End-of-linked list
        last=-1

        ! Row of matrix A
        DO ild=KldA(ieq),KldA(ieq+1)-1
          
          ! Column number to be added
          icol=KcolA(ild)

          ! Add element to the linked list
          ndg = ndg+1
          Kaux(icol) = last
          last = icol
        END DO
        
        ! Row of matrix B
        DO ild=KldB(ieq),KldB(ieq+1)-1
          
          ! Column number to be added
          icol=KcolB(ild)

          ! Add element to the linked list
          IF (Kaux(icol) == 0) THEN
            ndg = ndg+1
            Kaux(icol) = last
            last = icol
          END IF
        END DO

        ! We are done with row IEQ
        NA = NA+ndg

        ! Reset KAUX to zero
        DO idg=1,ndg
          jeq = Kaux(last)
          Kaux(last) = 0
          last = jeq
        END DO

      END DO
    END FUNCTION do_mat79mat79add_computeNA

    !**************************************************************
    ! Format 7/9-7/9 addition
    ! Perform symbolical matrix-matrix-addition C := A + B
    ! 
    ! Remark: This subroutine is a modified version of the APLB1
    !         subroutine taken from the SPARSEKIT library written
    !         by Youcef Saad.

    SUBROUTINE do_mat79mat79add_symb(neq,ncols,KldA,KcolA,cmatrixFormatA,&
        &KldB,KcolB,cmatrixFormatB,KldC,KcolC,Kdiagonal)

      INTEGER(PREC_DOFIDX), INTENT(IN) :: neq,ncols
      INTEGER, INTENT(IN) :: cmatrixFormatA,cmatrixFormatB
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA,KldB
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA,KcolB
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT) :: KldC
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(INOUT) :: KcolC
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Kdiagonal
      
      INTEGER :: ieq,ildA,ildB,ildC,ildendA,ildendB,icolA,icolB,icolC

      ! Initialization
      KldC(1)=1; ildC=1

      ! Loop over all rows
      DO ieq=1,neq

        ! Initialize column pointers for matrix A and B
        ildA=KldA(ieq); ildendA=KldA(ieq+1)-1
        ildB=KldB(ieq); ildendB=KldB(ieq+1)-1

        ! Check if diagonal entry needs to be stored at leading
        ! position for storage format CSR7. Then, both matrices A and
        ! B must be stored in storage format CSR7 so that the pointer
        ! ILDA and ILDB need to be increased by one (see below).
        IF (.NOT.PRESENT(Kdiagonal)) THEN
          KcolC(ildC) = ieq
          ildC = ildC+1
        END IF

        ! In any case, skip diagonal entry for matrix A and/or B if
        ! they are storage format 7
        IF (cmatrixFormatA == LSYSSC_MATRIX7) ildA = ildA+1
        IF (cmatrixFormatB == LSYSSC_MATRIX7) ildB = ildB+1

        ! For each row IEQ loop over the columns of matrix A and B
        ! simultaneously and collect the corresponding matrix entries
        DO
        
          ! Find next column number for matrices A and B
          IF (ildA <= ildendA) THEN
            icolA = KcolA(ildA)
          ELSE
            icolA = ncols+1
          END IF
          
          IF (ildB <= ildendB) THEN
            icolB = KcolB(ildB)
          ELSE
            icolB = ncols+1
          END IF

          ! We need to consider three different cases. Since the
          ! diagonal element needs to be treated separately, the
          ! auxiliary variable ICOLC is used to store the provisional
          ! position of the next matrix entry
          IF (icolA == icolB) THEN
            ! Processing same column in matrix A and B 
            icolC = icolA
            ildA = ildA+1
            ildB = ildB+1
            
          ELSEIF (icolA < icolB) THEN
            ! Processing column in matrix A only
            icolC = icolA
            ildA = ildA+1

          ELSE
            ! Processing column in matrix B only
            icolC = icolB
            ildB = ildB+1

          END IF

          ! Now, we have the provisional position ICOLC of the next
          ! matrix entry. If it corresponds to the diagonal entry,
          ! then special care must be taken.
          IF (icolC == ieq) THEN
            
            ! If KDIAGONAL is present, then all entries in row IEQ
            ! are stored continuously but the diagonal entry is
            ! additionally marked in KDIAGONAL. IF KDIAGONAL is not
            ! present then the diagonal entry is already stored in
            ! the first position of each row (see above).
            IF (PRESENT(Kdiagonal)) THEN
              Kdiagonal(ieq)=ildC
              KcolC(ildC) = icolC
              ildC = ildC+1
            END IF

          ELSE
            ! Off-diagonal entries are handled as usual
            KcolC(ildC)=icolC
            ildC = ildC+1
          END IF
          
          ! Check if column IEQ is completed for both matrices A and B
          IF (ildA > ildendA .AND. ildB > ildendB) EXIT
        END DO

        KldC(ieq+1)=ildC
      END DO
    END SUBROUTINE do_mat79mat79add_symb

    !**************************************************************
    ! Format 7/9-7/9 addition
    ! Perform numerical matrix-matrix-addition C := ca * A + cb * B
    ! 
    ! Remark: This subroutine is a modified version of the APLB1
    !         subroutine taken from the SPARSEKIT library written
    !         by Youcef Saad.
    !
    ! The vectors KLDC and KCOLC are not necessary at first glance.
    ! However, if the matrix C is allowed to have even a larger
    ! sparsity pattern as the "sum" of A and B, then they are
    ! required to find the correct positions in the final matrix C.
    !
    ! In this subroutine, KDIAGC is the vector which points to the
    ! position of the diagonal entries. If matrix C is stored in
    ! format CSR7 then KDIAGC corresponds to KLDC(1:NEQ). If matrix C
    ! is stored in format CSR9 then KDIAGC corresponds to KDIAGONALC.

    SUBROUTINE do_mat79mat79add_numb_dbledble(isizeIntl,neq,ncols,&
        KldA,KcolA,DaA,cA,KldB,KcolB,DaB,cB,KldC,KcolC,KdiagC,DaC)
      
      INTEGER(PREC_DOFIDX), INTENT(IN)               :: neq,ncols
      INTEGER, INTENT(IN)                            :: isizeIntl
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA,KldB,KldC,KdiagC
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA,KcolB,KcolC
      REAL(DP), INTENT(IN)                           :: cA,cB
      REAL(DP), DIMENSION(isizeIntl,*), INTENT(IN)   :: DaA
      REAL(DP), DIMENSION(isizeIntl,*), INTENT(IN)   :: DaB
      REAL(DP), DIMENSION(isizeIntl,*), INTENT(INOUT):: DaC
      
      INTEGER(PREC_DOFIDX) :: ieq
      INTEGER(PREC_VECIDX) :: ildA,ildB,ildC,ildendA,ildendB,ildendC
      INTEGER(PREC_MATIDX) :: icolA,icolB,icolC,idiagC

      ! Loop over all ROWS
      DO ieq=1,neq
        
        ! Initialize column pointers for matrix A, B and C
        ildA=KldA(ieq); ildendA=KldA(ieq+1)-1
        ildB=KldB(ieq); ildendB=KldB(ieq+1)-1
        ildC=KldC(ieq); ildendC=KldC(ieq+1)-1

        ! Initialize pointer to diagonal entry
        idiagC = KdiagC(ieq)
        
        ! Since the final value of the diagonal entry
        !    c_i,i = ca * a_i,i + cb * b_i,i
        ! is updated step-by-step, set the diagonal entry to zero
        DaC(:,idiagC) = 0._DP
        
        ! For each row IEQ loop over the columns of matrix A and B
        ! simultaneously and collect the corresponding matrix entries
        DO
          
          ! Find next column number for matrices A and B
          IF (ildA <= ildendA) THEN
            icolA = KcolA(ildA)
          ELSE
            icolA = ncols+1
          END IF
          
          IF (ildB <= ildendB) THEN
            icolB = KcolB(ildB)
          ELSE
            icolB = ncols+1
          END IF
          
          ! First, check if at least for one of the two matrices A
          ! and/or B the diagonal entry which requires special
          ! treatment has been reached. In this case, update (!!!)
          ! the diagonal entry of the resulting matrix C immediately
          ! and proceed to the next iteration
          IF (icolA == ieq) THEN
            IF (icolB == ieq) THEN
              
              ! 1. Case: For both matrices A and B the diagonal entry
              ! has been reached
              DaC(:,idiagC)=cA*DaA(:,ildA)+cB*DaB(:,ildB)
              ildA = ildA+1
              ildB = ildB+1
            ELSE

              ! 2. Case: For matrix A the diagonal entry has been
              ! reached. Hence, skip matrix B.
              DaC(:,idiagC)=DaC(:,idiagC)+cA*DaA(:,ildA)
              ildA = ildA+1
            END IF
          ELSEIF (icolB == ieq) THEN

            !   3. Case: For matrix B the diagonal entry has been
            !      reached. Hence, skip matrix A.
            DaC(:,idiagC)=DaC(:,idiagC)+cB*DaB(:,ildB)
            ildB = ildB+1

          ELSE
            
            ! For both matrices A and B we have to process off-
            ! -diagonal entries. Consider three different cases.
            ! 1.) The next column is the same for both matrices A and B
            ! 2.) The next column number is only present in matrix A
            ! 3.) The next column number is only present in matrix B
            !
            ! Since matrix C is implicitly allowed to posses matrix
            ! entries which are present neither in matrix A nor B, the
            ! position ILDC needs to be updated of the current column
            ! is not the diagonal entry whose position is known a
            ! priori
            
            IF (icolA == icolB) THEN
              
              ! 1. Case: Processing same column in matrix A and B

              ! Update column number for matrix C
              DO ildC=ildc,ildendC
                icolC=KcolC(ildC)
                IF (icolC == icolA) EXIT
                IF (icolC /= ieq) DaC(:,ildC) = 0._DP
              END DO

              DaC(:,ildC)=cA*DaA(:,ildA)+cB*DaB(:,ildB)
              ildA = ildA+1
              ildB = ildB+1
              ildC = ildC+1
              
            ELSEIF (icolA < icolB) THEN
              
              ! 2. Case: Processing column in matrix A only

              ! Update column number for matrix C
              DO ildC=ildC,ildendC
                icolC=KcolC(ildC)
                IF (icolC == icolA) EXIT
                IF (icolC /= ieq) DaC(:,ildC) = 0._DP
              END DO

              DaC(:,ildC)=cA*DaA(:,ildA)
              ildA = ildA+1
              ildC = ildC+1
              
            ELSE
              
              ! 3. Case: Processing column in matrix B only
              
              ! Update column number for matrix C
              DO ildC=ildC,ildendC
                icolC=KcolC(ildC)
                IF (icolC == icolB) EXIT
                IF (icolC /= ieq) DaC(:,ildC) = 0._DP
              END DO

              DaC(:,ildC)=cB*DaB(:,ildB)
              ildB = ildB+1
              ildC = ildC+1

            END IF
          END IF
          
          ! Check if column IEQ is completed for both matrices A and B
          IF (ildA > ildendA .AND. ildB > ildendB) EXIT
        END DO

        ! Since matrix C is allowed to have additional column entries
        ! which are not present in the "sum" of A and B, the
        ! remainder of C needs to be nullified by hand
        IF (ildC <= ildendC) THEN
          IF (KcolC(ildC) == ieq) THEN
            DaC(:,ildC+1:1:ildendC) = 0._DP
          ELSE
            DaC(:,ildC:1:ildendC)   = 0._DP
          END IF
        END IF
      END DO
    END SUBROUTINE do_mat79mat79add_numb_dbledble

    !**************************************************************
    ! Format 7/9-7/9 addition
    ! Perform numerical matrix-matrix-addition C := ca * A + cb * B
    ! 
    ! Remark: This subroutine is a modified version of the APLB1
    !         subroutine taken from the SPARSEKIT library written
    !         by Youcef Saad.
    !
    ! The vectors KLDC and KCOLC are not necessary at first glance.
    ! However, if the matrix C is allowed to have even a larger
    ! sparsity pattern as the "sum" of A and B, then they are
    ! required to find the correct positions in the final matrix C.
    !
    ! In this subroutine, KDIAGC is the vector which points to the
    ! position of the diagonal entries. If matrix C is stored in
    ! format CSR7 then KDIAGC corresponds to KLDC(1:NEQ). If matrix C
    ! is stored in format CSR9 then KDIAGC corresponds to KDIAGONALC.

    SUBROUTINE do_mat79mat79add_numb_dblesngl(isizeIntl,neq,ncols,&
        KldA,KcolA,DaA,cA,KldB,KcolB,FaB,cB,KldC,KcolC,KdiagC,DaC)
      
      INTEGER(PREC_DOFIDX), INTENT(IN)               :: neq,ncols
      INTEGER, INTENT(IN)                            :: isizeIntl
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA,KldB,KldC,KdiagC
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA,KcolB,KcolC
      REAL(DP), INTENT(IN)                           :: cA,cB
      REAL(DP), DIMENSION(isizeIntl,*), INTENT(IN)   :: DaA
      REAL(SP), DIMENSION(isizeIntl,*), INTENT(IN)   :: FaB
      REAL(DP), DIMENSION(isizeIntl,*), INTENT(INOUT):: DaC
      
      INTEGER(PREC_DOFIDX) :: ieq
      INTEGER(PREC_VECIDX) :: ildA,ildB,ildC,ildendA,ildendB,ildendC
      INTEGER(PREC_MATIDX) :: icolA,icolB,icolC,idiagC
      
      ! Loop over all ROWS
      DO ieq=1,neq
        
        ! Initialize column pointers for matrix A, B and C
        ildA=KldA(ieq); ildendA=KldA(ieq+1)-1
        ildB=KldB(ieq); ildendB=KldB(ieq+1)-1
        ildC=KldC(ieq); ildendC=KldC(ieq+1)-1

        ! Initialize pointer to diagonal entry
        idiagC = KdiagC(ieq)
        
        ! Since the final value of the diagonal entry
        !    c_i,i = ca * a_i,i + cb * b_i,i
        ! is updated step-by-step, set the diagonal entry to zero
        DaC(:,idiagC) = 0._DP
        
        ! For each row IEQ loop over the columns of matrix A and B
        ! simultaneously and collect the corresponding matrix entries
        DO
          
          ! Find next column number for matrices A and B
          IF (ildA <= ildendA) THEN
            icolA = KcolA(ildA)
          ELSE
            icolA = ncols+1
          END IF
          
          IF (ildB <= ildendB) THEN
            icolB = KcolB(ildB)
          ELSE
            icolB = ncols+1
          END IF
          
          ! First, check if at least for one of the two matrices A
          ! and/or B the diagonal entry which requires special
          ! treatment has been reached. In this case, update (!!!)
          ! the diagonal entry of the resulting matrix C immediately
          ! and proceed to the next iteration
          IF (icolA == ieq) THEN
            IF (icolB == ieq) THEN
          
              ! 1. Case: For both matrices A and B the diagonal entry
              ! has been reached
              DaC(:,idiagC)=cA*DaA(:,ildA)+cB*FaB(:,ildB)
              ildA = ildA+1
              ildB = ildB+1
            ELSE

              ! 2. Case: For matrix A the diagonal entry has been
              ! reached. Hence, skip matrix B.
              DaC(:,idiagC)=DaC(:,idiagC)+cA*DaA(:,ildA)
              ildA = ildA+1
            END IF
          ELSEIF (icolB == ieq) THEN

            !   3. Case: For matrix B the diagonal entry has been
            !      reached. Hence, skip matrix A.
            DaC(:,idiagC)=DaC(:,idiagC)+cB*FaB(:,ildB)
            ildB = ildB+1

          ELSE
            
            ! For both matrices A and B we have to process off-
            ! -diagonal entries. Consider three different cases.
            ! 1.) The next column is the same for both matrices A and B
            ! 2.) The next column number is only present in matrix A
            ! 3.) The next column number is only present in matrix B
            !
            ! Since matrix C is implicitly allowed to posses matrix
            ! entries which are present neither in matrix A nor B, the
            ! position ILDC needs to be updated of the current column
            ! is not the diagonal entry whose position is known a
            ! priori
            
            IF (icolA == icolB) THEN
              
              ! 1. Case: Processing same column in matrix A and B

              ! Update column number for matrix C
              DO ildC=ildc,ildendC
                icolC=KcolC(ildC)
                IF (icolC == icolA) EXIT
                IF (icolC /= ieq) DaC(:,ildC) = 0._DP
              END DO

              DaC(:,ildC)=cA*DaA(:,ildA)+cB*FaB(:,ildB)
              ildA = ildA+1
              ildB = ildB+1
              ildC = ildC+1
              
            ELSEIF (icolA < icolB) THEN
              
              ! 2. Case: Processing column in matrix A only

              ! Update column number for matrix C
              DO ildC=ildC,ildendC
                icolC=KcolC(ildC)
                IF (icolC == icolA) EXIT
                IF (icolC /= ieq) DaC(:,ildC) = 0._DP
              END DO

              DaC(:,ildC)=cA*DaA(:,ildA)
              ildA = ildA+1
              ildC = ildC+1
              
            ELSE
              
              ! 3. Case: Processing column in matrix B only
              
              ! Update column number for matrix C
              DO ildC=ildC,ildendC
                icolC=KcolC(ildC)
                IF (icolC == icolB) EXIT
                IF (icolC /= ieq) DaC(:,ildC) = 0._DP
              END DO

              DaC(:,ildC)=cB*FaB(:,ildB)
              ildB = ildB+1
              ildC = ildC+1

            END IF
          END IF
          
          ! Check if column IEQ is completed for both matrices A and B
          IF (ildA > ildendA .AND. ildB > ildendB) EXIT
        END DO

        ! Since matrix C is allowed to have additional column entries
        ! which are not present in the "sum" of A and B, the
        ! remainder of C needs to be nullified by hand
        IF (ildC <= ildendC) THEN
          IF (KcolC(ildC) == ieq) THEN
            DaC(:,ildC+1:1:ildendC) = 0._DP
          ELSE
            DaC(:,ildC:1:ildendC)   = 0._DP
          END IF
        END IF
      END DO
    END SUBROUTINE do_mat79mat79add_numb_dblesngl

    !**************************************************************
    ! Format 7/9-7/9 addition
    ! Perform numerical matrix-matrix-addition C := ca * A + cb * B
    ! 
    ! Remark: This subroutine is a modified version of the APLB1
    !         subroutine taken from the SPARSEKIT library written
    !         by Youcef Saad.
    !
    ! The vectors KLDC and KCOLC are not necessary at first glance.
    ! However, if the matrix C is allowed to have even a larger
    ! sparsity pattern as the "sum" of A and B, then they are
    ! required to find the correct positions in the final matrix C.
    !
    ! In this subroutine, KDIAGC is the vector which points to the
    ! position of the diagonal entries. If matrix C is stored in
    ! format CSR7 then KDIAGC corresponds to KLDC(1:NEQ). If matrix C
    ! is stored in format CSR9 then KDIAGC corresponds to KDIAGONALC.

    SUBROUTINE do_mat79mat79add_numb_sngldble(isizeIntl,neq,ncols,&
        KldA,KcolA,FaA,cA,KldB,KcolB,DaB,cB,KldC,KcolC,KdiagC,DaC)
      
      INTEGER(PREC_DOFIDX), INTENT(IN)               :: neq,ncols
      INTEGER, INTENT(IN)                            :: isizeIntl
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA,KldB,KldC,KdiagC
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA,KcolB,KcolC
      REAL(DP), INTENT(IN)                           :: cA,cB
      REAL(SP), DIMENSION(isizeIntl,*), INTENT(IN)   :: FaA
      REAL(DP), DIMENSION(isizeIntl,*), INTENT(IN)   :: DaB
      REAL(DP), DIMENSION(isizeIntl,*), INTENT(INOUT):: DaC

      INTEGER(PREC_DOFIDX) :: ieq
      INTEGER(PREC_VECIDX) :: ildA,ildB,ildC,ildendA,ildendB,ildendC
      INTEGER(PREC_MATIDX) :: icolA,icolB,icolC,idiagC
      
      ! Loop over all ROWS
      DO ieq=1,neq
        
        ! Initialize column pointers for matrix A, B and C
        ildA=KldA(ieq); ildendA=KldA(ieq+1)-1
        ildB=KldB(ieq); ildendB=KldB(ieq+1)-1
        ildC=KldC(ieq); ildendC=KldC(ieq+1)-1

        ! Initialize pointer to diagonal entry
        idiagC = KdiagC(ieq)
        
        ! Since the final value of the diagonal entry
        !    c_i,i = ca * a_i,i + cb * b_i,i
        ! is updated step-by-step, set the diagonal entry to zero
        DaC(:,idiagC) = 0._DP
        
        ! For each row IEQ loop over the columns of matrix A and B
        ! simultaneously and collect the corresponding matrix entries
        DO
          
          ! Find next column number for matrices A and B
          IF (ildA <= ildendA) THEN
            icolA = KcolA(ildA)
          ELSE
            icolA = ncols+1
          END IF
          
          IF (ildB <= ildendB) THEN
            icolB = KcolB(ildB)
          ELSE
            icolB = ncols+1
          END IF
          
          ! First, check if at least for one of the two matrices A
          ! and/or B the diagonal entry which requires special
          ! treatment has been reached. In this case, update (!!!)
          ! the diagonal entry of the resulting matrix C immediately
          ! and proceed to the next iteration
          IF (icolA == ieq) THEN
            IF (icolB == ieq) THEN
          
              ! 1. Case: For both matrices A and B the diagonal entry
              ! has been reached
              DaC(:,idiagC)=cA*FaA(:,ildA)+cB*DaB(:,ildB)
              ildA = ildA+1
              ildB = ildB+1
            ELSE

              ! 2. Case: For matrix A the diagonal entry has been
              ! reached. Hence, skip matrix B.
              DaC(:,idiagC)=DaC(:,idiagC)+cA*FaA(:,ildA)
              ildA = ildA+1
            END IF
          ELSEIF (icolB == ieq) THEN

            !   3. Case: For matrix B the diagonal entry has been
            !      reached. Hence, skip matrix A.
            DaC(:,idiagC)=DaC(:,idiagC)+cB*DaB(:,ildB)
            ildB = ildB+1

          ELSE
            
            ! For both matrices A and B we have to process off-
            ! -diagonal entries. Consider three different cases.
            ! 1.) The next column is the same for both matrices A and B
            ! 2.) The next column number is only present in matrix A
            ! 3.) The next column number is only present in matrix B
            !
            ! Since matrix C is implicitly allowed to posses matrix
            ! entries which are present neither in matrix A nor B, the
            ! position ILDC needs to be updated of the current column
            ! is not the diagonal entry whose position is known a
            ! priori
            
            IF (icolA == icolB) THEN
              
              ! 1. Case: Processing same column in matrix A and B

              ! Update column number for matrix C
              DO ildC=ildc,ildendC
                icolC=KcolC(ildC)
                IF (icolC == icolA) EXIT
                IF (icolC /= ieq) DaC(:,ildC) = 0._DP
              END DO

              DaC(:,ildC)=cA*FaA(:,ildA)+cB*DaB(:,ildB)
              ildA = ildA+1
              ildB = ildB+1
              ildC = ildC+1
              
            ELSEIF (icolA < icolB) THEN
              
              ! 2. Case: Processing column in matrix A only

              ! Update column number for matrix C
              DO ildC=ildC,ildendC
                icolC=KcolC(ildC)
                IF (icolC == icolA) EXIT
                IF (icolC /= ieq) DaC(:,ildC) = 0._DP
              END DO

              DaC(:,ildC)=cA*FaA(:,ildA)
              ildA = ildA+1
              ildC = ildC+1
              
            ELSE
              
              ! 3. Case: Processing column in matrix B only
              
              ! Update column number for matrix C
              DO ildC=ildC,ildendC
                icolC=KcolC(ildC)
                IF (icolC == icolB) EXIT
                IF (icolC /= ieq) DaC(:,ildC) = 0._DP
              END DO

              DaC(:,ildC)=cB*DaB(:,ildB)
              ildB = ildB+1
              ildC = ildC+1

            END IF
          END IF
          
          ! Check if column IEQ is completed for both matrices A and B
          IF (ildA > ildendA .AND. ildB > ildendB) EXIT
        END DO

        ! Since matrix C is allowed to have additional column entries
        ! which are not present in the "sum" of A and B, the
        ! remainder of C needs to be nullified by hand
        IF (ildC <= ildendC) THEN
          IF (KcolC(ildC) == ieq) THEN
            DaC(:,ildC+1:1:ildendC) = 0._DP
          ELSE
            DaC(:,ildC:1:ildendC)   = 0._DP
          END IF
        END IF
      END DO
    END SUBROUTINE do_mat79mat79add_numb_sngldble

    !**************************************************************
    ! Format 7/9-7/9 addition
    ! Perform numerical matrix-matrix-addition C := ca * A + cb * B
    ! 
    ! Remark: This subroutine is a modified version of the APLB1
    !         subroutine taken from the SPARSEKIT library written
    !         by Youcef Saad.
    !
    ! The vectors KLDC and KCOLC are not necessary at first glance.
    ! However, if the matrix C is allowed to have even a larger
    ! sparsity pattern as the "sum" of A and B, then they are
    ! required to find the correct positions in the final matrix C.
    !
    ! In this subroutine, KDIAGC is the vector which points to the
    ! position of the diagonal entries. If matrix C is stored in
    ! format CSR7 then KDIAGC corresponds to KLDC(1:NEQ). If matrix C
    ! is stored in format CSR9 then KDIAGC corresponds to KDIAGONALC.

    SUBROUTINE do_mat79mat79add_numb_snglsngl(isizeIntl,neq,ncols,&
        KldA,KcolA,FaA,cA,KldB,KcolB,FaB,cB,KldC,KcolC,KdiagC,FaC)
      
      INTEGER(PREC_DOFIDX), INTENT(IN)               :: neq,ncols
      INTEGER, INTENT(IN)                            :: isizeIntl
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: KldA,KldB,KldC,KdiagC
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: KcolA,KcolB,KcolC
      REAL(DP), INTENT(IN)                           :: cA,cB
      REAL(SP), DIMENSION(isizeIntl,*), INTENT(IN)   :: FaA
      REAL(SP), DIMENSION(isizeIntl,*), INTENT(IN)   :: FaB
      REAL(SP), DIMENSION(isizeIntl,*), INTENT(INOUT):: FaC
      
      INTEGER(PREC_DOFIDX) :: ieq
      INTEGER(PREC_VECIDX) :: ildA,ildB,ildC,ildendA,ildendB,ildendC
      INTEGER(PREC_MATIDX) :: icolA,icolB,icolC,idiagC
      
      ! Loop over all ROWS
      DO ieq=1,neq
        
        ! Initialize column pointers for matrix A, B and C
        ildA=KldA(ieq); ildendA=KldA(ieq+1)-1
        ildB=KldB(ieq); ildendB=KldB(ieq+1)-1
        ildC=KldC(ieq); ildendC=KldC(ieq+1)-1

        ! Initialize pointer to diagonal entry
        idiagC = KdiagC(ieq)
        
        ! Since the final value of the diagonal entry
        !    c_i,i = ca * a_i,i + cb * b_i,i
        ! is updated step-by-step, set the diagonal entry to zero
        FaC(:,idiagC) = 0._SP
        
        ! For each row IEQ loop over the columns of matrix A and B
        ! simultaneously and collect the corresponding matrix entries
        DO
          
          ! Find next column number for matrices A and B
          IF (ildA <= ildendA) THEN
            icolA = KcolA(ildA)
          ELSE
            icolA = ncols+1
          END IF
          
          IF (ildB <= ildendB) THEN
            icolB = KcolB(ildB)
          ELSE
            icolB = ncols+1
          END IF
          
          ! First, check if at least for one of the two matrices A
          ! and/or B the diagonal entry which requires special
          ! treatment has been reached. In this case, update (!!!)
          ! the diagonal entry of the resulting matrix C immediately
          ! and proceed to the next iteration
          IF (icolA == ieq) THEN
            IF (icolB == ieq) THEN
          
              ! 1. Case: For both matrices A and B the diagonal entry
              ! has been reached
              FaC(:,idiagC)=cA*FaA(:,ildA)+cB*FaB(:,ildB)
              ildA = ildA+1
              ildB = ildB+1
            ELSE

              ! 2. Case: For matrix A the diagonal entry has been
              ! reached. Hence, skip matrix B.
              FaC(:,idiagC)=FaC(:,idiagC)+cA*FaA(:,ildA)
              ildA = ildA+1
            END IF
          ELSEIF (icolB == ieq) THEN

            !   3. Case: For matrix B the diagonal entry has been
            !      reached. Hence, skip matrix A.
            FaC(:,idiagC)=FaC(:,idiagC)+cB*FaB(:,ildB)
            ildB = ildB+1

          ELSE
            
            ! For both matrices A and B we have to process off-
            ! -diagonal entries. Consider three different cases.
            ! 1.) The next column is the same for both matrices A and B
            ! 2.) The next column number is only present in matrix A
            ! 3.) The next column number is only present in matrix B
            !
            ! Since matrix C is implicitly allowed to posses matrix
            ! entries which are present neither in matrix A nor B, the
            ! position ILDC needs to be updated of the current column
            ! is not the diagonal entry whose position is known a
            ! priori
            
            IF (icolA == icolB) THEN
              
              ! 1. Case: Processing same column in matrix A and B

              ! Update column number for matrix C
              DO ildC=ildc,ildendC
                icolC=KcolC(ildC)
                IF (icolC == icolA) EXIT
                IF (icolC /= ieq) FaC(:,ildC) = 0._SP
              END DO

              FaC(:,ildC)=cA*FaA(:,ildA)+cB*FaB(:,ildB)
              ildA = ildA+1
              ildB = ildB+1
              ildC = ildC+1
              
            ELSEIF (icolA < icolB) THEN
              
              ! 2. Case: Processing column in matrix A only

              ! Update column number for matrix C
              DO ildC=ildC,ildendC
                icolC=KcolC(ildC)
                IF (icolC == icolA) EXIT
                IF (icolC /= ieq) FaC(:,ildC) = 0._SP
              END DO

              FaC(:,ildC)=cA*FaA(:,ildA)
              ildA = ildA+1
              ildC = ildC+1
              
            ELSE
              
              ! 3. Case: Processing column in matrix B only
              
              ! Update column number for matrix C
              DO ildC=ildC,ildendC
                icolC=KcolC(ildC)
                IF (icolC == icolB) EXIT
                IF (icolC /= ieq) FaC(:,ildC) = 0._SP
              END DO

              FaC(:,ildC)=cB*FaB(:,ildB)
              ildB = ildB+1
              ildC = ildC+1

            END IF
          END IF
          
          ! Check if column IEQ is completed for both matrices A and B
          IF (ildA > ildendA .AND. ildB > ildendB) EXIT
        END DO

        ! Since matrix C is allowed to have additional column entries
        ! which are not present in the "sum" of A and B, the
        ! remainder of C needs to be nullified by hand
        IF (ildC <= ildendC) THEN
          IF (KcolC(ildC) == ieq) THEN
            FaC(:,ildC+1:1:ildendC) = 0._SP
          ELSE
            FaC(:,ildC:1:ildendC)   = 0._SP
          END IF
        END IF
      END DO
    END SUBROUTINE do_mat79mat79add_numb_snglsngl
  END SUBROUTINE lsyssc_matrixLinearComb

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_swapVectors(rvector1,rvector2)

!<description>
    ! This subroutine swaps the content of two different vectors
!</description>

!<inputoutput>
    ! first scalar vector
    TYPE(t_vectorScalar), INTENT(INOUT) :: rvector1

    ! second scalar vector
    TYPE(t_vectorScalar), INTENT(INOUT) :: rvector2
!</inputoutput>
!</subroutine>

    ! local variables
    TYPE(t_vectorScalar) :: rvector

    rvector  = rvector1
    rvector%p_rspatialDiscretisation => rvector1%p_rspatialDiscretisation
    
    rvector1 = rvector2
    rvector1%p_rspatialDiscretisation => rvector2%p_rspatialDiscretisation
    
    rvector2 = rvector
    rvector2%p_rspatialDiscretisation => rvector%p_rspatialDiscretisation
  END SUBROUTINE lsyssc_swapVectors

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_infoMatrix(rmatrix)

!<description>
    ! This subroutine outputs information about the matrix
!</description>

!<input>
    ! scalar matrix
    TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
!</input>
!</subroutine>

    ! local variables
    INTEGER(PREC_MATIDX) :: isize

    CALL output_line ('ScalarMatrix:')
    CALL output_line ('-------------')
    CALL output_line ('cmatrixFormat:           '//TRIM(sys_siL(rmatrix%cmatrixFormat,15)))
    CALL output_line ('cinterleavematrixFormat: '//TRIM(sys_siL(rmatrix%cinterleaveMatrixFormat,15)))
    CALL output_line ('cdataType:               '//TRIM(sys_siL(rmatrix%cdataType,15)))
    CALL output_line ('imatrixSpec:             '//TRIM(sys_siL(rmatrix%imatrixSpec,15)))
    CALL output_line ('NA:                      '//TRIM(sys_siL(rmatrix%NA,15)))
    CALL output_line ('NEQ:                     '//TRIM(sys_siL(rmatrix%NEQ,15)))
    CALL output_line ('NCOLS:                   '//TRIM(sys_siL(rmatrix%NCOLS,15)))
    CALL output_line ('NVAR:                    '//TRIM(sys_siL(rmatrix%NVAR,15)))
    CALL output_line ('dscaleFactor:            '//TRIM(sys_sdL(rmatrix%dscaleFactor,2)))
    CALL output_line ('isortStrategy:           '//TRIM(sys_siL(rmatrix%isortStrategy,15)))
    CALL output_line ('h_IsortPermutation:      '//TRIM(sys_siL(rmatrix%h_IsortPermutation,15)))
    CALL output_line ('h_DA:                    '//TRIM(sys_siL(rmatrix%h_DA,15)))
    IF (rmatrix%h_DA /= ST_NOHANDLE) THEN
      CALL storage_getsize(rmatrix%h_DA,isize)
      SELECT CASE(rmatrix%cinterleaveMatrixFormat)
        CASE (LSYSSC_MATRIX1)
          CALL output_line ('DA memory usage:         '//&
              TRIM(sys_sdL(100/REAL(isize,DP)*rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR,2))//'%')
        CASE (LSYSSC_MATRIXD)
          CALL output_line ('DA memory usage:         '//&
              TRIM(sys_sdL(100/REAL(isize,DP)*rmatrix%NA*rmatrix%NVAR,2))//'%')
        CASE DEFAULT
          CALL output_line ('DA memory usage:         '//TRIM(sys_sdL(100/REAL(isize,DP)*rmatrix%NA,2))//'%')
        END SELECT
    END IF
    CALL output_line ('h_Kcol:                  '//TRIM(sys_siL(rmatrix%h_Kcol,15)))
    IF (rmatrix%h_Kcol /= ST_NOHANDLE) THEN
      CALL storage_getsize(rmatrix%h_Kcol,isize)
      CALL output_line ('Kcol memory usage:       '//TRIM(sys_sdL(100/REAL(isize,DP)*rmatrix%NA,2))//'%')
    END IF
    CALL output_line ('h_Kld:                   '//TRIM(sys_siL(rmatrix%h_Kld,15)))
    IF (rmatrix%h_Kld /= ST_NOHANDLE) THEN
      CALL storage_getsize(rmatrix%h_Kld,isize)
      IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) == 0) THEN
        CALL output_line ('Kld memory usage:        '//TRIM(sys_sdL(100/REAL(isize,DP)*(rmatrix%NEQ+1),2))//'%')
      ELSE
        CALL output_line ('Kld memory usage:        '//TRIM(sys_sdL(100/REAL(isize,DP)*(rmatrix%NCOLS+1),2))//'%')
      END IF
    END IF
    CALL output_line ('h_Kdiagonal:             '//TRIM(sys_siL(rmatrix%h_Kdiagonal,15)))
    IF (rmatrix%h_Kdiagonal /= ST_NOHANDLE) THEN
      CALL storage_getsize(rmatrix%h_Kdiagonal,isize)
      IF (IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) == 0) THEN
        CALL output_line ('Kdiagonal memory usage:  '//TRIM(sys_sdL(100/REAL(isize,DP)*(rmatrix%NEQ),2))//'%')
      ELSE
        CALL output_line ('Kdiagonl memory usage:   '//TRIM(sys_sdL(100/REAL(isize,DP)*(rmatrix%NCOLS),2))//'%')
      END IF
    END IF
  END SUBROUTINE lsyssc_infoMatrix

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_infoVector(rvector)

!<description>
    ! This subroutine outputs information about the vector
!</description>

!<input>
    ! scalar vector
    TYPE(t_vectorScalar), INTENT(IN) :: rvector
!</input>
!</subroutine>

    ! local variables
    INTEGER(PREC_VECIDX) :: isize

    CALL output_line ('ScalarVector:')
    CALL output_line ('-------------')
    CALL output_line ('cdataType:              '//TRIM(sys_siL(rvector%cdataType,15)))
    CALL output_line ('NEQ:                    '//TRIM(sys_siL(rvector%NEQ,15)))
    CALL output_line ('NVAR:                   '//TRIM(sys_siL(rvector%NVAR,15)))
    CALL output_line ('isortStrategy:          '//TRIM(sys_siL(rvector%isortStrategy,15)))
    CALL output_line ('h_IsortPermutation:     '//TRIM(sys_siL(rvector%h_IsortPermutation,15)))
    CALL output_line ('h_Ddata:                '//TRIM(sys_siL(rvector%h_Ddata,15)))
    IF (rvector%h_Ddata /= ST_NOHANDLE) THEN
      CALL storage_getsize(rvector%h_Ddata,isize)
      CALL output_line ('Ddata memory usage:     '//TRIM(sys_sdL(100/REAL(isize,DP)*rvector%NEQ*rvector%NVAR,2))//'%')
    END IF
    CALL output_line ('iidxFirstEntry:         '//TRIM(sys_siL(rvector%iidxFirstEntry,15)))
    WRITE(*,FMT='(A)')       '-------------------------'
  END SUBROUTINE lsyssc_infoVector

  ! ***************************************************************************
  
!<function>

  ELEMENTAL LOGICAL FUNCTION lsyssc_isMatrixStructureShared (rmatrix,&
      rmatrix2) RESULT (bresult)
  
!<description>
  ! Checks whether the structure of rmatrix belongs to that matrix or if it's
  ! shared among rmatrix and another matrix.
  !
  ! Note: If both matrices rmatrix and rmatrix2 are specified and both matrices
  !   don't contain structure data, the return value is of course FALSE.
!</description>

!<input>
  ! A scalar matrix to be checked.
  TYPE(t_matrixScalar), INTENT(IN) :: rmatrix
  
  ! OPTIONAL: A second matrix to be compared with rmatrix.
  ! If not specified, lsyssc_isMatrixStructureShared checks if rmatrix
  !   shares its structure with any other matrix.
  ! If specified, lsyssc_isMatrixStructureShared checks if rmatrix
  !   shares its entries with rmatrix2.
  TYPE(t_matrixScalar), INTENT(IN), OPTIONAL :: rmatrix2
!</input>

!<result>
  ! TRUE, if the structure arrays in rmatrix are shared among rmatrix and
  !   another matrix or with the matrix rmatrix2, respectively.
  ! FALSE, if the structure arrays in rmatrix belong to rmatrix.
!</result>

!</function>

    IF (.NOT. PRESENT(rmatrix2)) THEN
      ! General check
      bresult = IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .NE. 0
    ELSE
      ! Check if the structure is shared among rmatrix and rmatrix2.
      ! This depends on the matrix format...
      ! The matrix is declared as 'structure is shared' if at least one
      ! of the handles that define the structure is shared among the matrices.
      bresult = .FALSE.
      IF (rmatrix%cmatrixFormat .NE. rmatrix2%cmatrixFormat) RETURN
      SELECT CASE (rmatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX1,LSYSSC_MATRIXD)
        ! No structure, full matrix
        bresult = .FALSE.
      CASE (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
        bresult = ((rmatrix%h_Kcol .NE. ST_NOHANDLE) .AND. &
                   (rmatrix%h_Kcol .EQ. rmatrix2%h_Kcol)) .OR. &
                  ((rmatrix%h_Kld .NE. ST_NOHANDLE) .AND. &
                   (rmatrix%h_Kld .EQ. rmatrix2%h_Kld))
      CASE (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
        bresult = ((rmatrix%h_Kcol .NE. ST_NOHANDLE) .AND. &
                   (rmatrix%h_Kcol .EQ. rmatrix2%h_Kcol)) .OR. &
                  ((rmatrix%h_Kld .NE. ST_NOHANDLE) .AND. &
                   (rmatrix%h_Kld .EQ. rmatrix2%h_Kld)) .OR. &
                  ((rmatrix%h_Kdiagonal .NE. ST_NOHANDLE) .AND. &
                   (rmatrix%h_Kdiagonal .EQ. rmatrix2%h_Kdiagonal))
      END SELECT
    END IF

  END FUNCTION

  ! ***************************************************************************
  
!<function>

  ELEMENTAL LOGICAL FUNCTION lsyssc_isMatrixContentShared (rmatrix, &
      rmatrix2) RESULT (bresult)
  
!<description>
  ! Checks whether the content of rmatrix belongs to that matrix or if it's
  ! shared among rmatrix and another matrix.
  !
  ! Note: If both matrices rmatrix and rmatrix2 are specified and both matrices
  !   don''t contain content data, the return value is of course FALSE.
!</description>

!<input>
  ! A scalar matrix
  TYPE(t_matrixScalar), INTENT(IN) :: rmatrix

  ! OPTIONAL: A second matrix to be compared with rmatrix.
  ! If not specified, lsyssc_isMatrixContentShared checks if rmatrix
  !   shares its content with any other matrix.
  ! If specified, lsyssc_isMatrixStructureShared checks if rmatrix
  !   shares its content with rmatrix2.
  TYPE(t_matrixScalar), INTENT(IN), OPTIONAL :: rmatrix2
!</input>

!<result>
  ! TRUE, if the content arrays in rmatrix are shared among rmatrix and
  !   another matrix or with the matrix rmatrix2, respectively.
  ! FALSE, if the content arrays in rmatrix belong to rmatrix or if
  !   rmatrix does not have content data.
!</result>

!</function>

    IF (.NOT. PRESENT(rmatrix2)) THEN
      bresult = IAND(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .NE. 0
    ELSE
      ! Check if the content is shared among rmatrix and rmatrix2.
      ! This depends on the matrix format...
      ! The matrix is declared as 'content is shared' if at least one
      ! of the handles that define the structure is shared among the matrices.
      bresult = .FALSE.
      IF (rmatrix%cmatrixFormat .NE. rmatrix2%cmatrixFormat) RETURN
      SELECT CASE (rmatrix%cmatrixFormat)
      CASE (LSYSSC_MATRIX1,LSYSSC_MATRIX7,LSYSSC_MATRIX9,LSYSSC_MATRIXD, &
            LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9INTL)
        bresult = ((rmatrix%h_Da .NE. ST_NOHANDLE) .AND. &
                   (rmatrix%h_Da .EQ. rmatrix2%h_Da))
      END SELECT
    END IF

  END FUNCTION

  !****************************************************************************
  
!<function>
  
  ELEMENTAL LOGICAL FUNCTION lsyssc_isMatrixPresent (rmatrix, &
      bignoreScaleFactor) RESULT (bispresent)
  
!<description>
  ! This routine checks if the matrix is defined or empty.
!</description>
  
!<input>
  ! The matrix to be checked
  TYPE(t_matrixScalar), INTENT(IN)               :: rmatrix
  
  ! OPTIONAL: Whether to check the scaling factor.
  ! FALSE: A scaling factor of 0.0 disables a submatrix. 
  !        This is the standard setting.
  ! TRUE: The scaling factor is ignored.
  LOGICAL, INTENT(IN), OPTIONAL :: bignoreScaleFactor
!</input>

!<output>
  ! Whether the matrix structure realises an existing matrix exists or not.
!</output>  

!</function>
    LOGICAL :: bscale

    IF (PRESENT(bignoreScaleFactor)) THEN
      bscale = bignoreScaleFactor
    ELSE
      bscale = .FALSE.
    END IF

    bispresent = &
      (rmatrix%cmatrixFormat .NE. LSYSSC_MATRIXUNDEFINED) &
      .AND. ((.NOT. bscale) .OR. &
             (rmatrix%dscaleFactor .NE. 0.0_DP))

  END FUNCTION
    
  !****************************************************************************

!<function>
  
  PURE LOGICAL FUNCTION lsyssc_isMatrixSorted (rmatrix)
  
!<description>
  ! Returns whether a matrix is sorted or not.
!</description>
  
!<input>
  ! Matrix to check
  TYPE(t_matrixScalar), INTENT(IN)                  :: rmatrix
!</input>

!<result>
  ! Whether the matrix is sorted or not.
!</result>

!</function>

    lsyssc_isMatrixSorted = rmatrix%isortStrategy .GT. 0

  END FUNCTION

  !****************************************************************************

!<function>
  
  PURE LOGICAL FUNCTION lsyssc_isVectorSorted (rvector)
  
!<description>
  ! Returns whether a vector is sorted or not.
!</description>
  
!<input>
  ! Vector to check
  TYPE(t_vectorScalar), INTENT(IN)                  :: rvector
!</input>

!<result>
  ! Whether the vector is sorted or not.
!</result>

!</function>

    lsyssc_isVectorSorted = rvector%isortStrategy .GT. 0

  END FUNCTION

  !****************************************************************************

!<function>
  
  PURE LOGICAL FUNCTION lsyssc_hasMatrixStricture (rmatrix)
  
!<description>
  ! Returns whether a matrix has a structure or not.
  !
  ! Note that some matrix types (e.g. matrix-type 1 = full matrix)
  ! don't have a structure at all, so the routine always returns
  ! FALSE in such a case.
!</description>
  
!<input>
  ! Matrix to check
  TYPE(t_matrixScalar), INTENT(IN)                  :: rmatrix
!</input>

!<result>
  ! Whether the matrix has strucure arrays in memory or not.
!</result>

!</function>

    ! All up to now implemented matrix types use Kcol if they have a 
    ! structure...
    lsyssc_hasMatrixStricture = rmatrix%h_Kcol .NE. ST_NOHANDLE

  END FUNCTION

  !****************************************************************************

!<function>
  
  PURE LOGICAL FUNCTION lsyssc_hasMatrixContent (rmatrix)
  
!<description>
  ! Returns whether a matrix has a content or not.
!</description>
  
!<input>
  ! Matrix to check
  TYPE(t_matrixScalar), INTENT(IN)                  :: rmatrix
!</input>

!<result>
  ! Whether the matrix has a content array in memory or not.
!</result>

!</function>

    ! All up to now implemented matrix types save their data in h_Da.
    lsyssc_hasMatrixContent = rmatrix%h_Da .NE. ST_NOHANDLE

  END FUNCTION

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_spreadVector(rvector1, rvector2)

!<description>
    ! This subroutine spreads a scalar vector which is not stored in interleave
    ! format into another vector which is stored in interleave format.
!</description>

!<input>
    ! Scalar source vector
    TYPE(t_vectorScalar), INTENT(IN)    :: rvector1
!</input>

!<inputoutput>
    ! Scalar destination vector in interleave format
    TYPE(t_vectorScalar), INTENT(INOUT) :: rvector2
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER     :: p_Ddata1,p_Ddata2
    REAL(SP), DIMENSION(:), POINTER     :: p_Fdata1,p_Fdata2
    INTEGER(I32), DIMENSION(:), POINTER :: p_Idata1,p_Idata2
    
    ! Source vector must not be stored in interleave format
    IF (rvector1%NVAR .NE. 1) THEN
      CALL output_line('Source vector must not be stored in interleave format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadVector')
      CALL sys_halt()
    END IF

    ! Vectors must have the same size
    IF (rvector1%NEQ .NE. rvector2%NEQ) THEN
      CALL output_line('Vectors not compatible, different size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadVector')
      CALL sys_halt()
    END IF

    ! Vectors must have the data type
    IF (rvector1%cdataType .NE. rvector2%cdataType) THEN
      CALL output_line('Vectors not compatible, different data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadVector')
      CALL sys_halt()
    END IF

    ! isortStrategy < 0 means unsorted. Both unsorted is ok.
    
    IF ((rvector1%isortStrategy .GT. 0) .OR. &
        (rvector2%isortStrategy .GT. 0)) THEN
      
      IF (rvector1%isortStrategy .NE. &
          rvector2%isortStrategy) THEN
        CALL output_line('Vectors not compatible, differently sorted!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadVector')
        CALL sys_halt()
      END IF
    END IF

    IF (rvector1%h_isortPermutation .NE. &
        rvector2%h_isortPermutation) THEN
      CALL output_line('Vectors not compatible, differently sorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadVector')
      CALL sys_halt()
    END IF

    SELECT CASE (rvector1%cdataType)
    CASE (ST_DOUBLE)
      CALL lsyssc_getbase_double(rvector1, p_Ddata1)
      CALL lsyssc_getbase_double(rvector2, p_Ddata2)
      CALL do_spreadDble(p_Ddata1, rvector2%NVAR, rvector2%NEQ, p_Ddata2)
      
    CASE (ST_SINGLE)
      CALL lsyssc_getbase_single(rvector1, p_Fdata1)
      CALL lsyssc_getbase_single(rvector2, p_Fdata2)
      CALL do_spreadSngl(p_Fdata1, rvector2%NVAR, rvector2%NEQ, p_Fdata2)

    CASE (ST_INT)
      CALL lsyssc_getbase_int(rvector1, p_Idata1)
      CALL lsyssc_getbase_int(rvector2, p_Idata2)
      CALL do_spreadInt(p_Idata1, rvector2%NVAR, rvector2%NEQ, p_Idata2)
      
    CASE DEFAULT
      CALL output_line('Unsupported data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadVector')
      CALL sys_halt()
    END SELECT

  CONTAINS

    ! Here, the real working routines follow.

    !**************************************************************

    SUBROUTINE do_spreadDble(Ddata1, NVAR, NEQ, Ddata2)
      REAL(DP), DIMENSION(:), INTENT(IN)         :: Ddata1
      INTEGER, INTENT(IN)                        :: NVAR
      INTEGER(PREC_DOFIDX), INTENT(IN)           :: NEQ
      REAL(DP), DIMENSION(NVAR,NEQ), INTENT(OUT) :: Ddata2

      INTEGER(PREC_DOFIDX) :: ieq

      DO ieq = 1, NEQ
        Ddata2(:,ieq) = Ddata1(ieq)
      END DO
    END SUBROUTINE do_spreadDble

    !**************************************************************

    SUBROUTINE do_spreadSngl(Fdata1, NVAR, NEQ, Fdata2)
      REAL(SP), DIMENSION(:), INTENT(IN)         :: Fdata1
      INTEGER, INTENT(IN)                        :: NVAR
      INTEGER(PREC_DOFIDX), INTENT(IN)           :: NEQ
      REAL(SP), DIMENSION(NVAR,NEQ), INTENT(OUT) :: Fdata2

      INTEGER(PREC_DOFIDX) :: ieq

      DO ieq = 1, NEQ
        Fdata2(:,ieq) = Fdata1(ieq)
      END DO
    END SUBROUTINE do_spreadSngl

    !**************************************************************

    SUBROUTINE do_spreadInt(Idata1, NVAR, NEQ, Idata2)
      INTEGER(I32), DIMENSION(:), INTENT(IN)         :: Idata1
      INTEGER, INTENT(IN)                            :: NVAR
      INTEGER(PREC_DOFIDX), INTENT(IN)               :: NEQ
      INTEGER(I32), DIMENSION(NVAR,NEQ), INTENT(OUT) :: Idata2

      INTEGER(PREC_DOFIDX) :: ieq

      DO ieq = 1, NEQ
        Idata2(:,ieq) = Idata1(ieq)
      END DO
    END SUBROUTINE do_spreadInt
  END SUBROUTINE lsyssc_spreadVector

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_spreadMatrix(rmatrix1, rmatrix2)

!<description>
    ! This subroutine spreads a scalar matrix which is not stored in interleave
    ! format into another matrix which is stored in interleave format.
!</description>

!<input>
    ! Scalar source matrix
    TYPE(t_matrixScalar), INTENT(IN)    :: rmatrix1
!</input>

!<inputoutput>
    ! Scalar destination matrix in interleave format
    TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix2
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER     :: p_Ddata1,p_Ddata2
    REAL(SP), DIMENSION(:), POINTER     :: p_Fdata1,p_Fdata2
    INTEGER(I32), DIMENSION(:), POINTER :: p_Idata1,p_Idata2

    ! Source matrices must not be stored in interleave format
    IF (rmatrix1%NVAR .NE. 1) THEN
      CALL output_line('Source matrix must not be stored in interleave format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadMatrix')
      CALL sys_halt()
    END IF

    ! Check if matrices are compatible
    CALL lsyssc_isMatrixMatrixCompatible(rmatrix1,rmatrix2)
    
    ! Ok, now we can copy the matrices
    SELECT CASE(rmatrix2%cinterleavematrixFormat)
      
    CASE (LSYSSC_MATRIXUNDEFINED)
      ! Destination matrix is identical to source matrix
      CALL lsyssc_duplicateMatrix (rmatrix1,rmatrix2,&
          LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE)
      
    CASE (LSYSSC_MATRIX1)
      SELECT CASE (rmatrix1%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double(rmatrix1, p_Ddata1)
        CALL lsyssc_getbase_double(rmatrix2, p_Ddata2)
        CALL do_spreadDble(p_Ddata1, rmatrix2%NVAR, rmatrix2%NVAR,&
                           rmatrix2%NA, p_Ddata2)
        
      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single(rmatrix1, p_Fdata1)
        CALL lsyssc_getbase_single(rmatrix2, p_Fdata2)
        CALL do_spreadSngl(p_Fdata1, rmatrix2%NVAR, rmatrix2%NVAR,&
                           rmatrix2%NA, p_Fdata2)

      CASE (ST_INT)
        CALL lsyssc_getbase_int(rmatrix1, p_Idata1)
        CALL lsyssc_getbase_int(rmatrix2, p_Idata2)
        CALL do_spreadInt(p_Idata1, rmatrix2%NVAR, rmatrix2%NVAR,&
                          rmatrix2%NA, p_Idata2)

      CASE DEFAULT
        CALL output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadMatrix')
        CALL sys_halt()
      END SELECT
      
    CASE (LSYSSC_MATRIXD)
      SELECT CASE (rmatrix1%cdataType)
      CASE (ST_DOUBLE)
        CALL lsyssc_getbase_double(rmatrix1, p_Ddata1)
        CALL lsyssc_getbase_double(rmatrix2, p_Ddata2)
        CALL do_spreadDble(p_Ddata1, rmatrix2%NVAR, 1,&
                           rmatrix2%NA, p_Ddata2)

      CASE (ST_SINGLE)
        CALL lsyssc_getbase_single(rmatrix1, p_Fdata1)
        CALL lsyssc_getbase_single(rmatrix2, p_Fdata2)
        CALL do_spreadSngl(p_Fdata1, rmatrix2%NVAR, 1,&
                           rmatrix2%NA, p_Fdata2)

      CASE (ST_INT)
        CALL lsyssc_getbase_int(rmatrix1, p_Idata1)
        CALL lsyssc_getbase_int(rmatrix2, p_Idata2)
        CALL do_spreadInt(p_Idata1, rmatrix2%NVAR, 1,&
                          rmatrix2%NA, p_Idata2)

      CASE DEFAULT
        CALL output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadMatrix')
        CALL sys_halt()
      END SELECT
      
    CASE DEFAULT
      CALL output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadMatrix')
      CALL sys_halt()
    END SELECT
    
  CONTAINS

    ! Here, the real working routines follow.

    !**************************************************************

    SUBROUTINE do_spreadDble(Ddata1, NVAR, MVAR, NA, Ddata2)
      REAL(DP), DIMENSION(:), INTENT(IN)             :: Ddata1
      INTEGER, INTENT(IN)                            :: NVAR,MVAR
      INTEGER(PREC_MATIDX), INTENT(IN)               :: NA
      REAL(DP), DIMENSION(NVAR,MVAR,NA), INTENT(OUT) :: Ddata2

      INTEGER(PREC_MATIDX) :: ia

      DO ia = 1, NA
        Ddata2(:,:,ia) = Ddata1(ia)
      END DO
    END SUBROUTINE do_spreadDble

    !**************************************************************

    SUBROUTINE do_spreadSngl(Fdata1, NVAR, MVAR, NA, Fdata2)
      REAL(SP), DIMENSION(:), INTENT(IN)             :: Fdata1
      INTEGER, INTENT(IN)                            :: NVAR,MVAR
      INTEGER(PREC_MATIDX), INTENT(IN)               :: NA
      REAL(SP), DIMENSION(NVAR,MVAR,NA), INTENT(OUT) :: Fdata2

      INTEGER(PREC_MATIDX) :: ia

      DO ia = 1, NA
        Fdata2(:,:,ia) = Fdata1(ia)
      END DO
    END SUBROUTINE do_spreadSngl

    !**************************************************************

    SUBROUTINE do_spreadInt(Idata1, NVAR, MVAR, NA, Idata2)
      INTEGER(I32), DIMENSION(:), INTENT(IN)             :: Idata1
      INTEGER, INTENT(IN)                                :: NVAR,MVAR
      INTEGER(PREC_MATIDX), INTENT(IN)                   :: NA
      INTEGER(I32), DIMENSION(NVAR,MVAR,NA), INTENT(OUT) :: Idata2

      INTEGER(PREC_MATIDX) :: ia

      DO ia = 1, NA
        Idata2(:,:,ia) = Idata1(ia)
      END DO
    END SUBROUTINE do_spreadInt   
  END SUBROUTINE lsyssc_spreadMatrix

  !****************************************************************************

!<subroutine>

  SUBROUTINE lsyssc_packVector(rvector1, rvector2, ivar)

!<description>
    ! This subroutine packs a scalar vector which is stored in interleave
    ! format into another vector which is not stored in interleave format.
!</description>

!<input>
    ! Scalar source vector in interleave format
    TYPE(t_vectorScalar), INTENT(IN)    :: rvector1

    ! Number of the variable to pack
    INTEGER, INTENT(IN)                 :: ivar
!</input>

!<inputoutput>
    ! Scalar destination vector
    TYPE(t_vectorScalar), INTENT(INOUT) :: rvector2
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER     :: p_Ddata1,p_Ddata2
    REAL(SP), DIMENSION(:), POINTER     :: p_Fdata1,p_Fdata2
    INTEGER(I32), DIMENSION(:), POINTER :: p_Idata1,p_Idata2

    ! Source vector must be stored in interleave format
    IF (rvector1%NVAR .LT. ivar) THEN
      CALL output_line('Source vector does not provide variable IVAR!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_packVector')
      CALL sys_halt()
    END IF

    ! Destination vector must not be stored in interleave format
    IF (rvector2%NVAR .NE. 1) THEN
      CALL output_line('Destination vector must not be stored in interleave format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_packVector')
      CALL sys_halt()
    END IF

    ! Vectors must have the same size
    IF (rvector1%NEQ .NE. rvector2%NEQ) THEN
      CALL output_line('Vectors not compatible, different size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_packVector')
      CALL sys_halt()
    END IF

    ! Vectors must have the data type
    IF (rvector1%cdataType .NE. rvector2%cdataType) THEN
      CALL output_line('Vectors not compatible, different data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_packVector')
      CALL sys_halt()
    END IF

    ! isortStrategy < 0 means unsorted. Both unsorted is ok.
    
    IF ((rvector1%isortStrategy .GT. 0) .OR. &
        (rvector2%isortStrategy .GT. 0)) THEN
      
      IF (rvector1%isortStrategy .NE. &
          rvector2%isortStrategy) THEN
        CALL output_line('Vectors not compatible, differently sorted!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_packVector')
        CALL sys_halt()
      END IF
    END IF

    IF (rvector1%h_isortPermutation .NE. &
        rvector2%h_isortPermutation) THEN
      CALL output_line('Vectors not compatible, differently sorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_packVector')
      CALL sys_halt()
    END IF

    SELECT CASE (rvector1%cdataType)
    CASE (ST_DOUBLE)
      CALL lsyssc_getbase_double(rvector1, p_Ddata1)
      CALL lsyssc_getbase_double(rvector2, p_Ddata2)
      CALL do_packDble(p_Ddata1, rvector1%NVAR, rvector1%NEQ, ivar, p_Ddata2)
      
    CASE (ST_SINGLE)
      CALL lsyssc_getbase_single(rvector1, p_Fdata1)
      CALL lsyssc_getbase_single(rvector2, p_Fdata2)
      CALL do_packSngl(p_Fdata1, rvector1%NVAR, rvector1%NEQ, ivar, p_Fdata2)

    CASE (ST_INT)
      CALL lsyssc_getbase_int(rvector1, p_Idata1)
      CALL lsyssc_getbase_int(rvector2, p_Idata2)
      CALL do_packInt(p_Idata1, rvector1%NVAR, rvector1%NEQ, ivar, p_Idata2)
      
    CASE DEFAULT
      CALL output_line('Unsupported data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_packVector')
      CALL sys_halt()
    END SELECT

  CONTAINS

    ! Here, the real working routines follow.

    !**************************************************************

    SUBROUTINE do_packDble(Ddata1, NVAR, NEQ, ivar, Ddata2)
      REAL(DP), DIMENSION(NVAR,NEQ), INTENT(IN) :: Ddata1
      INTEGER, INTENT(IN)                       :: NVAR
      INTEGER(PREC_DOFIDX), INTENT(IN)          :: NEQ
      INTEGER, INTENT(IN)                       :: ivar
      REAL(DP), DIMENSION(:), INTENT(OUT)       :: Ddata2

      INTEGER(PREC_DOFIDX) :: ieq

      DO ieq = 1, NEQ
        Ddata2(ieq) = Ddata1(ivar,ieq)
      END DO
    END SUBROUTINE do_packDble

    !**************************************************************

    SUBROUTINE do_packSngl(Fdata1, NVAR, NEQ, ivar, Fdata2)
      REAL(SP), DIMENSION(NVAR,NEQ), INTENT(IN) :: Fdata1
      INTEGER, INTENT(IN)                       :: NVAR
      INTEGER(PREC_DOFIDX), INTENT(IN)          :: NEQ
      INTEGER, INTENT(IN)                       :: ivar
      REAL(SP), DIMENSION(:), INTENT(OUT)       :: Fdata2

      INTEGER(PREC_DOFIDX) :: ieq

      DO ieq = 1, NEQ
        Fdata2(ieq) = Fdata1(ivar,ieq)
      END DO
    END SUBROUTINE do_packSngl

    !**************************************************************

    SUBROUTINE do_packInt(Idata1, NVAR, NEQ, ivar, Idata2)
      INTEGER(I32), DIMENSION(NVAR,NEQ), INTENT(IN) :: Idata1
      INTEGER, INTENT(IN)                           :: NVAR
      INTEGER(PREC_DOFIDX), INTENT(IN)              :: NEQ
      INTEGER, INTENT(IN)                           :: ivar
      INTEGER(I32), DIMENSION(:), INTENT(OUT)       :: Idata2

      INTEGER(PREC_DOFIDX) :: ieq

      DO ieq = 1, NEQ
        Idata2(ieq) = Idata1(ivar,ieq)
      END DO
    END SUBROUTINE do_packInt
  END SUBROUTINE lsyssc_packVector
END MODULE

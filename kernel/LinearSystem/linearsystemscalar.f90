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
!#  1.) lsyssc_createVecDirect / lsyssc_createVecDirectIntl
!#      -> Create a simple scalar vector of length NEQ
!#
!#  2.) lsyssc_createVecByDiscr / lsyssc_createVecByDiscrIntl
!#      -> Create a vector based on a scalar discretisation structure
!#
!#  3.) lsyssc_createVecIndMat
!#      -> Create a vector according to the matrix
!#
!#  4.) lsyssc_scalarProduct
!#      -> Calculate the scalar product of two vectors
!#
!#  5.) lsyssc_scalarMatVec
!#      -> Multiply a scalar matrix (or its transpose) with a scalar vector
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
!# 10.) lsyssc_setSortStrategy
!#      -> Attaches a sorting strategy to a vector/matrix.
!#         The sorting can be activated/deactivated at any time.
!#
!# 11.) lsyssc_sortVector
!#      -> Resort the entries of a vector or unsort them according to a
!#         previously attached sorting strategy
!#
!# 12.) lsyssc_sortMatrix
!#      -> Resort the entries of a matrix or unsort them according to a
!#         previously attached sorting strategy.
!#
!# 13.) lsyssc_synchroniseSort = lsyssc_synchroniseSortVecVec
!#      -> Synchronises the sorting strategy of a vector according to another
!#         vector.
!#
!# 14.) lsyssc_synchroniseSort = lsyssc_synchroniseSortMatVec
!#      -> Synchronises the sorting strategy of a vector according to a matrix.
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
!#      lsyssc_transposeMatrixInSitu,
!#      lsyssc_transposeMatrixDirect
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
!# 40.) lsyssc_matrixLinearComb / lsyssc_matrixLinearCombIndexed
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
!# 56.) lsyssc_createFullMatrix
!#      -> Create a full square or rectangular matrix in matrix format 1
!#
!# 57.) lsyssc_assignDiscrDirectMat
!#      -> Assign a discretisation to a matrix
!#
!# 58.) lsyssc_createFpdbObjectVec
!#      -> Creates an ObjectItem representing a scalar vector
!#
!# 59.) lsyssc_createFpdbObjectMat
!#      -> Creates an ObjectItem representing a scalar matrix
!#
!# 60.) lsyssc_restoreFpdbObjectVec
!#      -> Restores a scalar vector from an ObjectItem
!#
!# 61.) lsyssc_restoreFpdbObjectMat
!#      -> Restores a scalar matrix from an ObjectItem
!#
!# 62.) lsyssc_unshareMatrix
!#      -> Renders a matrix independent, resets the sharing state
!#
!# 63.) lsyssc_unshareVector
!#      -> Renders a vector independent, resets the sharing state
!#
!# 64.) lsyssc_isExplicitMatrix1D
!#      -> Checks whether a given matrix explicitly exists in memory as
!#         1D array
!#
!# 65.) lsyssc_checkDiscretisation
!#      -> Checks whether a given discretisation structure rdiscretisation is
!#         basically compatible to a given vector
!#
!# 66.) lsyssc_setDataTypeMatrix
!#      -> Converts the data type of a scalar matrix
!#
!# 67.) lsyssc_setDataTypeVector
!#      -> Converts the data type of a scalar vector
!#
!# 68.) lsyssc_moveMatrix
!#      -> Moves a matrix to another matrix.
!#
!# 69.) lsyssc_addConstant
!#      -> Adds a constant to a vector.
!#
!# 70.) lsyssc_swapMatrices
!#      -> Swaps two matrices
!#
!# 71.) lsyssc_sortMatrixCols
!#      -> Sorts the columns of a format 9 matrix to be in ascending order.
!#
!# 72.) lsyssc_createRowCMatrix
!#      -> Creates a row-compressed matrix (containing only nonzero rows).
!#
!# 73.) lsyssc_getbase_KrowIdx
!#      -> Returns the pointer of the KrowIdx-array
!#
!# 74.) lsyssc_copyH2D_Vector
!#      -> Copies the data of a vector from the host memory
!#         to the memory of the coprocessor device.
!#
!# 75.) lsyssc_copyD2H_Vector
!#      -> Copies the data of a vector from the memory of the
!#         coprocessor device to the host memory
!#
!# 76.) lsyssc_initPerfConfig
!#      -> Initialises the global performance configuration
!#
!# 77.) lsyssc_createMatrixSymmPart
!#      -> Create the symmetric part of the matrix
!#
!# 78.) lsyssc_createMatrixAsymmPart
!#      -> Create the antisymmetric part of the matrix
!#
!# 79.) lsyssc_copyH2D_Matrix
!#      -> Copies the data of a matrix from the host memory
!#         to the memory of the coprocessor device.
!#
!# 80.) lsyssc_copyD2H_Matrix
!#      -> Copies the data of a matrix from the memory of the
!#         coprocessor device to the host memory
!#
!# 81.) lsyssc_calcDeterminant
!#      -> Calculates the determinant of a square matrix
!#
!# 82.) lsyssc_calcGerschgorin
!#      -> Apply the theorem of Gerschgorin to get approximations
!#         for the min/max eigenvalues of a matrix
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
!#     -> Copies the data vector Da of a matrix to another without checks.
!#        If the destination array does not exist, it is created.
!#
!# 5.) lsyssc_auxcopy_Kcol
!#     -> Copies the data vector Kcol of a matrix to another without checks
!#        If the destination array does not exist, it is created.
!#
!# 6.) lsyssc_auxcopy_Kld
!#     -> Copies the data vector Kld of a matrix to another without checks
!#        If the destination array does not exist, it is created.
!#
!# 7.) lsyssc_auxcopy_Kdiagonal
!#     -> Copies the data vector Kdiagonal of a matrix to another without checks
!#        If the destination array does not exist, it is created.
!#
!# 8.) lsyssc_createEmptyMatrixStub
!#     -> Creates an empty matrix without allocating data. Basic tags are
!#        initialised.
!#
!# 9.) lsyssc_createEmptyMatrix9
!#     -> Manually creates an empty matrix in format 9 for being
!#        build with lsyssc_setRowMatrix9.
!#
!# 10.) lsyssc_setRowMatrix9
!#      -> Adds or replaces a row in a format-9 matrix.
!#
!# 11.) lsyssc_genEdgeList
!#      -> Generates a list of edges of a given matrix
!#
!# 12.) lsyssc_regroupEdgeList
!#      -> Regroups a list of edges into independent groups which can be
!#         processed individually, i.e. in parallel
!#
!# 13.) lsyssc_calcDimsFromMatrix
!#      -> Calculates dimensions NA, NEQ, and NEDGE from matrix
!#
!# 14.) lsyssc_createGraphFromMatrix
!#      -> Creates a graph data structure from a given matrix
!#
!# 15.) lsyssc_createMatrixFromGraph
!#      -> Creates a matrix from a graph data structure
!# </purpose>
!##############################################################################

module linearsystemscalar

!$use omp_lib
  use arraylistInt
  use arraylistInt_Int
  use dofmapping
  use fpersistence
  use fsystem
  use genoutput
  use graph
  use linearalgebra
  use mapInt_Int
  use perfconfig
  use spatialdiscretisation
  use storage
  use sortstrategybase
  use uuid

  implicit none

  private

!<constants>

!<constantblock description="Global constants for scalar vectors/matrices">

  ! Maximum number of tags that can be assigned to a scalar vector or matrix.
  integer, parameter, public :: LSYSSC_MAXTAGS = 16

!</constantblock>

!<constantblock description="Global format flags for matrices">

  ! Unidentified matrix format
  integer, parameter, public :: LSYSSC_MATRIXUNDEFINED = 0

  ! Identifier for matrix format 1 - full matrix.
  ! Important matrix properties defining the matrix:
  ! NEQ = Number of rows, NCOLS = Number of columns, NA = Number of entries,
  ! h_Da = handle to matrix entries
  integer, parameter, public :: LSYSSC_MATRIX1 = 1

  ! Identifier for matrix format 9 - CSR
  ! Important matrix properties defining the matrix:
  ! NEQ = Number of rows, NCOLS = Number of columns, NA = Number of entries,
  ! h_Da        = handle to matrix entries,
  ! h_Kcol      = handle to column structure,
  ! h_Kld       = handle to row structure,
  ! h_Kdiagonal = handle to diagonal pointer
  integer, parameter, public :: LSYSSC_MATRIX9 = 9

  ! Identifier for matrix format 9ROWC - CSR row compressed.
  ! The same as format 9, but all zero-rows removed.
  ! Used to provide a compressed difference to a 'parent' matrix.
  ! All non-zero rows must have exactly the same length as in the 'parent'
  ! matrix. KrowIdx gives an index for every row in the compressed matrix
  ! which is the coresponding row in the 'parent' matrix.
  ! Important matrix properties defining the matrix:
  ! NEQ = Number of rows, NCOLS = Number of columns, NA = Number of entries,
  ! h_Da        = handle to matrix entries,
  ! h_Kcol      = handle to column structure,
  ! h_Kld       = handle to row structure,
  ! h_KrowIdx   = handle to the row index structure
  integer, parameter, public :: LSYSSC_MATRIX9ROWC = 900

  ! Identifier for matrix format 7 - CSR with diagonal element in front
  ! Important matrix properties defining the matrix:
  ! NEQ = Number of rows, NCOLS = Number of columns, NA = Number of entries,
  ! h_Da        = handle to matrix entries,
  ! h_Kcol      = handle to column structure,
  ! h_Kld       = handle to row structure
  integer, parameter, public :: LSYSSC_MATRIX7 = 7

  ! Identifier for matrix format D - Diagonal matrix, only entries on the
  ! main diagonal.
  ! Important matrix properties defining the matrix:
  ! NEQ = NCOLS = NA = Number of rows = Number of columns = Number of entries,
  ! h_Da        = handle to matrix entries
  integer, parameter, public :: LSYSSC_MATRIXD = 20

  ! Identifier for matrix format 7intl - CSR interleaved with diagonal element in front
  ! Important matrix properties defining the matrix:
  ! NEQ = Number of rows, NCOLS = Number of columns,
  ! NA = Number of entries, NVAR = Number of variables
  ! h_Da        = handle to matrix entries,
  ! h_Kcol      = handle to column structure,
  ! h_Kld       = handle to row structure
  integer, parameter, public :: LSYSSC_MATRIX7INTL = 70

  ! Identifier for matrix format 9intl - CSR interleaved
  ! Important matrix properties defining the matrix:
  ! NEQ = Number of rows, NCOLS = Number of columns,
  ! NA = Number of entries, NVAR = Number of variables
  ! h_Da        = handle to matrix entries,
  ! h_Kcol      = handle to column structure,
  ! h_Kld       = handle to row structure,
  ! h_Kdiagonal = handle to diagonal pointer
  integer, parameter, public :: LSYSSC_MATRIX9INTL = 90

!</constantblock>

!<constantblock description="Flags for the matrix specification bitfield">

  ! Standard matrix
  integer(I32), parameter, public :: LSYSSC_MSPEC_STANDARD =        0

  ! Matrix structure is a copy of another matrix, shared via the same
  ! handles.
  integer(I32), parameter, public :: LSYSSC_MSPEC_STRUCTUREISCOPY = 2**0

  ! Matrix content is a copy of another matrix, shared via the same
  ! handles.
  integer(I32), parameter, public :: LSYSSC_MSPEC_CONTENTISCOPY   = 2**1

  ! Complete matrix is duplicate of another matrix and shares structure
  ! and entries via the same pointers
  integer(I32), parameter, public :: LSYSSC_MSPEC_ISCOPY = LSYSSC_MSPEC_STRUCTUREISCOPY +&
                                                   LSYSSC_MSPEC_CONTENTISCOPY

  ! Matrix is saved transposed.
  ! To use a matrix in a transposed way, the application has to
  ! 1.) set this flag in the imatrixSpec bitfield
  ! 2.) exchange the values in t_matrixScalar\%NEQ and t_matrixScalar\%NCOLS
  !     of the matrix structure.
  integer(I32), parameter, public :: LSYSSC_MSPEC_TRANSPOSED =      2**2

  ! Matrix not present in memory
  integer(I32), parameter, public :: LSYSSC_MSPEC_NOTINMEMORY =     2**3

!</constantblock>

!<constantblock description="KIND values for matrix/vector data">

  ! kind value for indices in matrices
  integer, parameter, public :: PREC_MATIDX = I32

  ! kind value for indices in vectors
  integer, parameter, public :: PREC_VECIDX = I32

  ! kind value for precision that should be used in matrices
  integer, parameter, public :: PREC_MATRIX = DP

!</constantblock>

!<constantblock description="Constants for duplicating a matrix">

  ! Do not set up the content/structure of the destination matrix, ignore
  ! any previous structure/content
  integer, parameter, public :: LSYSSC_DUP_IGNORE = 0

  ! Removes any existing matrix content/structure from the destination matrix.
  ! Releases memory if necessary.
  integer, parameter, public :: LSYSSC_DUP_REMOVE = 1

  ! Removes any existing matrix content from the destination matrix.
  ! No memory is released, handles are simply dismissed.
  integer, parameter, public :: LSYSSC_DUP_DISMISS = 2

  ! The destination matrix recveives the same handles for matrix content/structure
  ! as the source matrix  and therefore shares the same content/structure.
  integer, parameter, public :: LSYSSC_DUP_SHARE = 3

  ! The destination matrix gets a copy of the content of rsourceMatrix.
  ! If necessary, new memory is allocated.
  ! If the destination matrix  already contains allocated memory that belongs
  ! to that matrix, content/structure data is simply copied from rsourceMatrix into that.
  ! Note that this respects the ownership! I.e. if the destination matrix is not
  ! the owner of the content/structure data arrays, new memory is allocated to
  ! prevent the actual owner from getting destroyed!
  integer, parameter, public :: LSYSSC_DUP_COPY = 4

  ! The destination matrix gets a copy of the content of rsourceMatrix.
  ! If necessary, new memory is allocated.
  ! If the destination matrix  already contains allocated memory, content/structure
  ! data is simply copied from rsourceMatrix into that.
  ! The ownership of the content/data arrays is not respected, i.e. if the
  ! destination matrix is not the owner, the actual owner of the data arrays is
  ! modified, too!
  integer, parameter, public :: LSYSSC_DUP_COPYOVERWRITE = 5

  ! Duplicate by ownership. What belongs to the source matrix is copied
  ! (the same as LSYSSC_DUP_COPY). What belongs even to another matrix than
  ! the source matrix is shared (the same as LSYSSC_DUP_SHARE, .
  integer, parameter, public :: LSYSSC_DUP_ASIS = 6

  ! New memory is allocated for the structure/content in the same size as
  ! in the source matrix but no data is copied; the arrays are left uninitialised.
  integer, parameter, public :: LSYSSC_DUP_EMPTY = 7

  ! Copy the basic matrix information but do not copy handles.
  ! Set all handles of dynamic information to ST_NOHANDLE, so the matrix
  ! 'looks like' the old but has no dynamic data associated.
  integer, parameter, public :: LSYSSC_DUP_TEMPLATE   = 8

!</constantblock>

!<constantblock description="Constants for transposing a matrix">

  ! Do not transpose the matrix, simply mark the matrix as transposed
  ! by changing the flag in imatrixSpec
  integer, parameter, public :: LSYSSC_TR_VIRTUAL    = 0

  ! Transpose only the matrix structure
  integer, parameter, public :: LSYSSC_TR_STRUCTURE = 1

  ! Transpose only the matrix entries
  integer, parameter, public :: LSYSSC_TR_CONTENT   = 2

  ! Transpose the full matrix
  integer, parameter, public :: LSYSSC_TR_ALL       = 3

  ! Do not transpose the matrix. Copy the matrix in memory and mark the
  ! matrix as transposed by changing the flag in imatrixSpec
  integer, parameter, public :: LSYSSC_TR_VIRTUALCOPY = 4

!</constantblock>

!<constantblock description="Constants for initialising matrix entries when allocating">

  ! Let the entries of a matrix undefined
  integer, parameter, public :: LSYSSC_SETM_UNDEFINED = -1

  ! Clear the entries of a matrix when allocating
  integer, parameter, public :: LSYSSC_SETM_ZERO      = 0

  ! Set the entries of a matrix to 1 when allocating
  integer, parameter, public :: LSYSSC_SETM_ONE       = 1

!</constantblock>

!<constantblock description="Constants for lumping of matrices">

  ! Standard lumping; extract diagonal from a given matrix.
  integer, parameter, public :: LSYSSC_LUMP_STD       = 0

  ! Diagonal lumping; add all offdiagonal entries to the diagonal and take the diagonial
  integer, parameter, public :: LSYSSC_LUMP_DIAG      = 1

!</constantblock>

!<constantblock description="Constants for content stored at each edge">

  ! Store only the node numbers
  integer, parameter, public :: LSYSSC_EDGELIST_NODESONLY   = 0

  ! Store only the absolute matrix positions
  integer, parameter, public :: LSYSSC_EDGELIST_POSONLY     = 1

  ! Store the node numbers and the matrix positions
  integer, parameter, public :: LSYSSC_EDGELIST_NODESANDPOS = 2

!</constantblock>

!<constantblock description="Constants for defining the algorithm for edge coloring">

  ! Greedy edge-coloring algorithm
  integer, parameter, public :: LSYSSC_EDGECOLORING_GREEDY  = 0

  ! Edge-coloring algorithm due to Nishizeki, Terada, Leven
  integer, parameter, public :: LSYSSC_EDGECOLORING_NTL     = 1

!</constantblock>

!</constants>

!<types>
!<typeblock>

  ! A scalar vector that can be used by scalar linear algebra routines.

  type t_vectorScalar

    ! Universally unique identifier
    type(t_uuid) :: ruuid

    ! Length of the vector; not necessarily = SIZE(p_Ddata), as this
    ! scalar vector might be a subvector of a larger vector allocated
    ! on the heap!
    integer :: NEQ = 0

    ! Number of local variables; in general, scalar vectors of size
    ! NEQ posses NEQ entries. However, scalar vectors can be
    ! interleaved, that is, each of the NEQ entries stores NVAR local
    ! variables. In this case, NEQ remains unmodified but NVAR>1 such
    ! that the physical length of the vector is NEQ*NVAR.
    integer :: NVAR = 1

    ! Data type of the entries in the vector. Either ST_SINGLE or
    ! ST_DOUBLE.
    integer :: cdataType = ST_DOUBLE

    ! Handle identifying the vector entries.
    ! = ST_NOHANDLE if not allocated on the global heap.
    integer :: h_Ddata = ST_NOHANDLE
    
    ! Specifies if the vector is resorted.
    logical :: bisSorted = .false.
    
    ! Pointer to the sorting strategy. May point to NULL()
    ! if no sorting strategy is attached.
    type(t_sortStrategy), pointer :: p_rsortStrategy => null()

    ! Start position of the vector data in the array identified by
    ! h_Ddata. Normally = 1. Can be set to > 1 if the vector is a subvector
    ! in a larger memory block allocated on the heap.
    integer :: iidxFirstEntry = 1

    ! Integer tags. This array of integer values can be used to store
    ! auxiliary tags which depend on the application.
    integer, dimension(LSYSSC_MAXTAGS) :: ITags = 0

    ! Real tags. This array of real values can be used to store
    ! auxiliary tags which depend on the application.
    real(DP), dimension(LSYSSC_MAXTAGS) :: DTags = 0._DP

    ! Is set to true, if the handle h_Ddata belongs to another vector,
    ! i.e. when this vector shares data with another vector.
    ! This is usually used for block vectors containing a couple of
    ! scalar subvectors; in this case, the block vector is the actual
    ! 'owner' of the handle!
    logical :: bisCopy = .false.

    ! A pointer to the spatial discretisation or NULL(), if the vector is
    ! just an array, not belonging to any discretisation.
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr => null()

  end type

  public :: t_vectorScalar

!</typeblock>

!<typeblock>

  ! A scalar matrix that can be used by scalar linear algebra routines.

  type t_matrixScalar

    ! Universally unique identifier
    type(t_uuid) :: ruuid

    ! Format-tag. Identifies the format of the matrix.
    ! Can take one of the LSYSSC_MATRIXx format flags.
    integer :: cmatrixFormat = LSYSSC_MATRIXUNDEFINED

    ! Format-tag. Identifies the format of the interleaved matrix.
    ! Can take one of the LSYSSC_MATRIX1/D format flags.
    integer :: cinterleavematrixFormat = LSYSSC_MATRIXUNDEFINED

    ! Matrix specification tag. This is a bitfield coming from an OR
    ! combination of different LSYSSC_MSPEC_xxxx constants and specifies
    ! various details of the matrix. If it is =LSYSSC_MSPEC_STANDARD,
    ! the matrix is a usual matrix that needs no special handling.
    integer(I32) :: imatrixSpec = LSYSSC_MSPEC_STANDARD

    ! Number of elements in the matrix
    integer :: NA  = 0

    ! Number of equations = rows in the matrix.
    ! Remark: When the application sets the LSYSSC_MSPEC_TRANSPOSED flag
    !  in imatrixSpec to indicate a transposed matrix, this has to be
    !  exchanged with NCOLS to indicate the number of equations in the
    !  transposed matrix!
    integer :: NEQ = 0

    ! Number of columns in the matrix.
    ! Remark: When the application sets the LSYSSC_MSPEC_TRANSPOSED flag
    !  in imatrixSpec to indicate a transposed matrix, this has to be
    !  exchanged with NROWS to indicate the number of columns in the
    !  transposed matrix!
    integer :: NCOLS = 0

    ! Number of variables for interleaved matrix.
    ! Remark: When an interleaved matrix is defined, then each
    ! position of the symbolic matrix structure can store a quadratic
    ! submatrix of dimension NVAR x NVAR. This is not to be confused
    ! with layers which means that several matrices share the same
    ! structure but their values can be treated individually. For
    ! interleaved matrices, the matrix-vector multiplication, etc. is
    ! adapted to the local dense blocks.
    integer :: NVAR = 1

    ! Number of nonzero rows. Only Format 9ROWC.
    integer :: NNZROWS = 0

    ! Multiplier for matrix entries. All entries in the matrix are
    ! scaled by this multiplier when doing Matrix-vector multiplication.
    ! Note: This parameter is not supported by all algorithms, many
    ! algorithms will simply stop when this factor is <> 1, or they
    ! might ignore it. Therefore, use this factor with care!
    real(DP) :: dscaleFactor = 1.0_DP

    ! Specifies if the columns in the matrix are resorted.
    logical :: bcolumnsSorted = .false.

    ! Pointer to the sorting strategy for the columns. May point to NULL()
    ! if no sorting strategy is attached.
    type(t_sortStrategy), pointer :: p_rsortStrategyColumns => null()

    ! Specifies if the rows in the matrix are resorted.
    logical :: browsSorted = .false.

    ! Pointer to the sorting strategy for the rows. May point to NULL()
    ! if no sorting strategy is attached.
    ! May point to p_rsortStrategyColumns if the sorting strategy of
    ! rows and columns match.
    type(t_sortStrategy), pointer :: p_rsortStrategyRows => null()

    ! Data type of the entries in the vector. Either ST_SINGLE or
    ! ST_DOUBLE.
    integer :: cdataType = ST_DOUBLE

    ! Format-7 and Format-9: Handle identifying the elements in the matrix
    !REAL(PREC_MATRIX), DIMENSION(:), POINTER :: Da         => NULL()
    integer :: h_Da = ST_NOHANDLE

    ! Format-7 and Format-9: Handle identifying the column structure
    !integer, DIMENSION(:), POINTER :: KCOL       => NULL()
    integer :: h_Kcol = ST_NOHANDLE

    ! Format-7 and Format-9: Handle identifying the row structure
    ! Formal-9ROWC: Handle identifyinf the row structure of the parent matrix.
    !INTEGER, DIMENSION(:), POINTER :: KLD        => NULL()
    integer :: h_Kld = ST_NOHANDLE

    ! Format-9ROWC: Handle to an array saving the nonzero row indices.
    ! KrowIdx(1:NNZROWS) saves for every row in the compressed
    ! matrix the corresponding row in the 'parent' matrix.
    ! KrowIdx(NNZROWS+1:NEQ+NNZROWS) saves for every row in the 'parent' 
    ! matrix the corresponding row number in the compressed matrix -- 
    ! or 0 if the row does not exist in the compressed matrix.
    integer :: h_KrowIdx = ST_NOHANDLE

    ! Format-9: Similar to row structure. Handle top array of length NEQ.
    ! For each row, pointer to the first element on the upper
    ! triangular part of the matrix. For matrices with elements on the
    ! diagonal, this points to the diagonal element in each row
    ! (or to the first off-diagonal element above the diagonal if there
    ! is no diagonal element, respectively).
    ! For matrices with no elements in the upper right part of the
    ! matrix, this points to the first element on the next line.
    !INTEGER, DIMENSION(:), POINTER :: Kdiagonal  => NULL()
    integer :: h_Kdiagonal = ST_NOHANDLE

    ! Integer tags. This array of integer values can be used to store
    ! auxiliary tags which depend on the application.
    integer, dimension(LSYSSC_MAXTAGS) :: ITags = 0

    ! Real tags. This array of real values can be used to store
    ! auxiliary
    ! tags which depend on the application.
    real(DP), dimension(LSYSSC_MAXTAGS) :: DTags = 0._DP

    ! A pointer to the spatial discretisation for the trial functions.
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrTrial => null()

    ! A pointer to the spatial discretisation for the test functions.
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrTest => null()

    ! Flag: Trial and Test functions in all element distributions
    ! are the same. If FALSE, there is at least one element distribution
    ! with different trial and test functions.
    logical :: bidenticalTrialAndTest = .true.

  end type

  public :: t_matrixScalar

!</typeblock>

!</types>

  !*****************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: lsyssc_perfconfig

  !*****************************************************************************

  interface lsyssc_createVector
    module procedure lsyssc_createVecDirect
    module procedure lsyssc_createVecDirectIntl
    module procedure lsyssc_createVecIndMat
    module procedure lsyssc_createVecByDiscr
    module procedure lsyssc_createVecByDiscrIntl
  end interface

  public :: lsyssc_createVector
  public :: lsyssc_createVecDirect
  public :: lsyssc_createVecDirectIntl
  public :: lsyssc_createVecIndMat
  public :: lsyssc_createVecByDiscr
  public :: lsyssc_createVecByDiscrIntl

  interface lsyssc_resizeVector
    module procedure lsyssc_resizeVectorDirect
    module procedure lsyssc_resizeVectorIndMat
    module procedure lsyssc_resizeVectorIndirect
    module procedure lsyssc_resizeVectorByDiscr
  end interface

  public :: lsyssc_resizeVector
  public :: lsyssc_resizeVectorDirect
  public :: lsyssc_resizeVectorIndMat
  public :: lsyssc_resizeVectorIndirect
  public :: lsyssc_resizeVectorByDiscr

  interface lsyssc_resizeMatrix
    module procedure lsyssc_resizeMatrixDirect
    module procedure lsyssc_resizeMatrixIndirect
  end interface

  public :: lsyssc_resizeMatrix
  public :: lsyssc_resizeMatrixDirect
  public :: lsyssc_resizeMatrixIndirect

  interface lsyssc_getbase_double
    module procedure lsyssc_getbaseVector_double
    module procedure lsyssc_getbaseMatrixDa_double
  end interface

  public :: lsyssc_getbase_double

  interface lsyssc_getbase_single
    module procedure lsyssc_getbaseVector_single
    module procedure lsyssc_getbaseMatrixDa_single
  end interface

  public :: lsyssc_getbase_single

  interface lsyssc_getbase_int
    module procedure lsyssc_getbaseVector_int
    module procedure lsyssc_getbaseMatrixDa_int
  end interface

  public :: lsyssc_getbase_int

  interface lsyssc_isMatrixCompatible
    module procedure lsyssc_isMatrixVectorCompatible
    module procedure lsyssc_isMatrixMatrixCompatible
  end interface

  public :: lsyssc_isMatrixCompatible
  public :: lsyssc_isMatrixVectorCompatible
  public :: lsyssc_isMatrixMatrixCompatible

  interface lsyssc_matrixLinearComb
    module procedure lsyssc_matrixLinearComb
    module procedure lsyssc_matrixLinearComb2
  end interface lsyssc_matrixLinearComb
  
  public :: lsyssc_matrixLinearComb

  interface lsyssc_setSortStrategy
    module procedure lsyssc_setSortStrategyVec
    module procedure lsyssc_setSortStrategyMat
  end interface

  public :: lsyssc_setSortStrategy

  public :: lsyssc_matrixLinearCombIndexed
  public :: lsyssc_initPerfConfig
  public :: lsyssc_scalarProduct
  public :: lsyssc_scalarMatVec
  public :: lsyssc_releaseMatrix
  public :: lsyssc_releaseVector
  public :: lsyssc_duplicateMatrix
  public :: lsyssc_duplicateVector
  public :: lsyssc_synchroniseSortVecVec
  public :: lsyssc_synchroniseSortMatVec
  public :: lsyssc_isVectorCompatible
  public :: lsyssc_getbase_Kcol
  public :: lsyssc_getbase_Kld
  public :: lsyssc_getbase_Kdiagonal
  public :: lsyssc_addIndex
  public :: lsyssc_vectorNorm
  public :: lsyssc_invertedDiagMatVec
  public :: lsyssc_clearMatrix
  public :: lsyssc_initialiseIdentityMatrix
  public :: lsyssc_convertMatrix
  public :: lsyssc_copyVector
  public :: lsyssc_scaleVector
  public :: lsyssc_clearVector
  public :: lsyssc_vectorLinearComb
  public :: lsyssc_copyMatrix
  public :: lsyssc_transposeMatrix
  public :: lsyssc_transposeMatrixInSitu
  public :: lsyssc_transposeMatrixDirect
  public :: lsyssc_allocEmptyMatrix
  public :: lsyssc_lumpMatrixScalar
  public :: lsyssc_scaleMatrix
  public :: lsyssc_multMatMat
  public :: lsyssc_swapVectors
  public :: lsyssc_swapMatrices
  public :: lsyssc_isMatrixStructureShared
  public :: lsyssc_isMatrixContentShared
  public :: lsyssc_createDiagMatrixStruc
  public :: lsyssc_clearOffdiags
  public :: lsyssc_hasMatrixStructure
  public :: lsyssc_hasMatrixContent
  public :: lsyssc_releaseMatrixContent
  public :: lsyssc_spreadVector
  public :: lsyssc_spreadMatrix
  public :: lsyssc_packVector
  public :: lsyssc_createFullMatrix
  public :: lsyssc_assignDiscrDirectMat
  public :: lsyssc_createFpdbObjectVec
  public :: lsyssc_createFpdbObjectMat
  public :: lsyssc_restoreFpdbObjectVec
  public :: lsyssc_restoreFpdbObjectMat
  public :: lsyssc_unshareMatrix
  public :: lsyssc_unshareVector
  public :: lsyssc_isExplicitMatrix1D
  public :: lsyssc_checkDiscretisation
  public :: lsyssc_setDataTypeMatrix
  public :: lsyssc_setDataTypeVector
  public :: lsyssc_moveMatrix
  public :: lsyssc_addConstant
  public :: lsyssc_sortMatrixCols
  public :: lsyssc_createRowCMatrix
  public :: lsyssc_getbase_KrowIdx
  public :: lsyssc_copyH2D_Vector
  public :: lsyssc_copyD2H_Vector
  public :: lsyssc_copyH2D_Matrix
  public :: lsyssc_copyD2H_Matrix
  public :: lsyssc_createMatrixSymmPart
  public :: lsyssc_createMatrixAsymmPart
  public :: lsyssc_calcDeterminant
  public :: lsyssc_calcGerschgorin
  public :: lsyssc_sortVector
  public :: lsyssc_sortMatrix

  public :: lsyssc_rebuildKdiagonal
  public :: lsyssc_infoMatrix
  public :: lsyssc_infoVector
  public :: lsyssc_auxcopy_da
  public :: lsyssc_auxcopy_Kcol
  public :: lsyssc_auxcopy_Kld
  public :: lsyssc_auxcopy_Kdiagonal
  public :: lsyssc_createEmptyMatrixStub
  public :: lsyssc_setRowMatrix9
  public :: lsyssc_createEmptyMatrix9
  public :: lsyssc_genEdgeList
  public :: lsyssc_regroupEdgeList
  public :: lsyssc_calcDimsFromMatrix
  public :: lsyssc_createGraphFromMatrix
  public :: lsyssc_createMatrixFromGraph

contains

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_initPerfConfig(rperfconfig)

!<description>
  ! This routine initialises the global performance configuration
!</description>

!<input>
  ! OPTIONAL: performance configuration that should be used to initialise
  ! the global performance configuration. If not present, the values of
  ! the legacy constants is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
!</subroutine>

    if (present(rperfconfig)) then
      lsyssc_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(lsyssc_perfconfig)
    end if

  end subroutine lsyssc_initPerfConfig

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_isVectorCompatible (rvector1,rvector2,bcompatible)

!<description>
  ! Checks whether two vectors are compatible to each other, i.e. share
  ! the same structure, size and sorting strategy.
!</description>

!<input>
  ! The first vector
  type(t_vectorScalar), intent(in) :: rvector1

  ! The second vector
  type(t_vectorScalar), intent(in) :: rvector2
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether the vectors are compatible or not.
  ! If not given, an error will inform the user if the two vectors are
  ! not compatible and the program will halt.
  logical, intent(out), optional :: bcompatible
!</output>

!</subroutine>

  ! We assume that we are not compatible
  if (present(bcompatible)) bcompatible = .false.

  ! Vectors must have the same size
  if (rvector1%NEQ  .ne. rvector2%NEQ .or. &
      rvector1%NVAR .ne. rvector2%NVAR) then
    if (present(bcompatible)) then
      bcompatible = .false.
      return
    else
      call output_line("Vectors not compatible, different size!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isVectorCompatible")
      call sys_halt()
    end if
  end if

  ! Check the sorting.
  if (rvector1%bisSorted .neqv. rvector2%bisSorted) then
    if (present(bcompatible)) then
      bcompatible = .false.
      return
    else
      call output_line("Matrix/vector not compatible, differently sorted!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixVectorCompatible")
      call sys_halt()
    end if
  end if
  
  if (rvector1%bisSorted) then
    if (rvector1%p_rsortStrategy%ctype .ne. rvector2%p_rsortStrategy%ctype) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line("Matrix/vector not compatible, differently sorted!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixVectorCompatible")
        call sys_halt()
      end if
    end if
  end if

  ! Ok, they are compatible
  if (present(bcompatible)) bcompatible = .true.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_isMatrixVectorCompatible (rvector,rmatrix,btransposed,bcompatible)

!<description>
  ! Checks whether a vector and a matrix are compatible to each other, i.e.
  ! share the same structure, size and sorting strategy.
!</description>

!<input>
  ! The vector
  type(t_vectorScalar), intent(in) :: rvector

  ! The matrix
  type(t_matrixScalar), intent(in) :: rmatrix

  ! Check for rvector being compatible from the left or from the right
  ! to the matrix rmatrix.
  ! =FALSE: Check whether matrix-vector product $A*x$ is possible
  ! =TRUE : Check whether matrix-vector product $x^T*A = A^T*x$ is possible
  logical, intent(in) :: btransposed
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether vector and matrix are compatible or not.
  ! If not given, an error will inform the user if the vector/matrix are
  ! not compatible and the program will halt.
  logical, intent(out), optional :: bcompatible
!</output>

!</subroutine>

  integer :: NCOLS

  ! We assume that we are not compatible
  if (present(bcompatible)) bcompatible = .false.

  if (.not. btransposed) then
    NCOLS = rmatrix%NCOLS
  else
    NCOLS = rmatrix%NEQ
  end if

  ! Vector/Matrix must have the same size
  if (rvector%NEQ .ne. NCOLS) then
    if (present(bcompatible)) then
      bcompatible = .false.
      return
    else
      call output_line("Vector/Matrix not compatible, different block structure!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixVectorCompatible")
      call sys_halt()
    end if
  end if

  ! Check the sorting.
  if (.not. btransposed) then

    ! Sorting of the columns must match the vector

    if (rmatrix%bcolumnsSorted .neqv. rvector%bisSorted) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line("Matrix/vector not compatible, differently sorted!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixVectorCompatible")
        call sys_halt()
      end if
    end if
    
    if (rmatrix%bcolumnsSorted) then
      if (rmatrix%p_rsortStrategyColumns%ctype .ne. rvector%p_rsortStrategy%ctype) then
        if (present(bcompatible)) then
          bcompatible = .false.
          return
        else
          call output_line("Matrix/vector not compatible, differently sorted!",&
              OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixVectorCompatible")
          call sys_halt()
        end if
      end if
    end if
    
  else

    ! Sorting of the rows must match the vector

    if (rmatrix%browsSorted .neqv. rvector%bisSorted) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line("Matrix/vector not compatible, differently sorted!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixVectorCompatible")
        call sys_halt()
      end if
    end if
    
    if (rmatrix%browsSorted) then
      if (rmatrix%p_rsortStrategyRows%ctype .ne. rvector%p_rsortStrategy%ctype) then
        if (present(bcompatible)) then
          bcompatible = .false.
          return
        else
          call output_line("Matrix/vector not compatible, differently sorted!",&
              OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixVectorCompatible")
          call sys_halt()
        end if
      end if
    end if

  end if
    
  ! Ok, they are compatible
  if (present(bcompatible)) bcompatible = .true.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_isMatrixMatrixCompatible (rmatrix1,rmatrix2,bcompatible)

!<description>
  ! Checks whether two matrices are compatible to each other, i.e.
  ! share the same structure, size and sorting strategy.
!</description>

!<input>
  ! The vector
  type(t_matrixScalar), intent(in) :: rmatrix1

  ! The matrix
  type(t_matrixScalar), intent(in) :: rmatrix2
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether vector and matrix are compatible or not.
  ! If not given, an error will inform the user if the vector/matrix are
  ! not compatible and the program will halt.
  logical, intent(out), optional :: bcompatible
!</output>

!</subroutine>

  ! We assume that we are not compatible
  if (present(bcompatible)) bcompatible = .false.

  ! Matrices must have the same size
  if ((rmatrix1%NEQ .ne. rmatrix2%NEQ) .or. &
      (rmatrix1%NCOLS .ne. rmatrix2%NCOLS) .or. &
      (rmatrix1%NA .ne. rmatrix2%NA)) then
    if (present(bcompatible)) then
      bcompatible = .false.
      return
    else
      call output_line("Matrices not compatible, different block structure!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixMatrixCompatible")
      call sys_halt()
    end if
  end if

  ! Matrices must have the same sparsity pattern
  if (rmatrix1%NVAR .eq. rmatrix2%NVAR) then

    ! We can perform restrictive checks, since both matrices have the
    ! same number of internal variables
    if ((rmatrix1%cmatrixFormat .ne. rmatrix2%cmatrixFormat) .or. &
        (rmatrix1%cinterleavematrixFormat .ne. rmatrix2%cinterleavematrixFormat)) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line("Matrices not compatible, different sparsity pattern!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixMatrixCompatible")
        call sys_halt()
      end if
    end if

  else

    ! We have to be more careful.
    if ((rmatrix1%cmatrixFormat .ne. rmatrix2%cmatrixFormat) .and.&
        (rmatrix1%cmatrixFormat .ne. 10*rmatrix2%cmatrixFormat) .and.&
        (10*rmatrix1%cmatrixFormat .ne. rmatrix2%cmatrixFormat)) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line("Matrices not compatible, different sparsity pattern!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixMatrixCompatible")
        call sys_halt()
      end if
    end if

  end if

  ! Check the sorting.
  if ((rmatrix1%bcolumnsSorted .neqv. rmatrix2%bcolumnsSorted) .or. &
      (rmatrix1%browsSorted .neqv. rmatrix2%browsSorted)) then
    if (present(bcompatible)) then
      bcompatible = .false.
      return
    else
      call output_line("Matrices not compatible, differently sorted!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixMatrixCompatible")
      call sys_halt()
    end if
  end if
  
  if (rmatrix1%bcolumnsSorted) then
    if (rmatrix1%p_rsortStrategyColumns%ctype .ne. rmatrix2%p_rsortStrategyColumns%ctype) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line("Matrices not compatible, differently sorted!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixMatrixCompatible")
        call sys_halt()
      end if
    end if
  end if

  if (rmatrix1%browsSorted) then
    if (rmatrix1%p_rsortStrategyRows%ctype .ne. rmatrix2%p_rsortStrategyRows%ctype) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line("Matrices not compatible, differently sorted!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_isMatrixMatrixCompatible")
        call sys_halt()
      end if
    end if
  end if

  ! Ok, they are compatible
  if (present(bcompatible)) bcompatible = .true.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_checkDiscretisation (rvector,rdiscretisation,bcompatible)

!<description>
  ! Checks whether a given discretisation structure rdiscretisation is
  ! basically compatible to a given vector (concerning NEQ, cubature formulas
  ! on element distribution, element compatibility,...)
!</description>

!<input>
  ! The vector which is to be checked.
  type(t_vectorScalar), intent(in) :: rvector

  ! A scalar discretisation structure to check against the vector.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether vector and matrix are compatible or not.
  ! If not given, an error will inform the user if the vector/matrix are
  ! not compatible and the program will halt.
  logical, intent(out), optional :: bcompatible
!</output>

!</subroutine>

    integer :: NEQ

    ! We assume that we are compatible
    if (present(bcompatible)) bcompatible = .true.

    ! NEQ must be correct.
    NEQ = dof_igetNDofGlob(rdiscretisation)
    if (NEQ .ne. rvector%NEQ) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        call output_line("Vector not compatible to discretisation, different NEQ!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_checkDiscretisation")
        call sys_halt()
      end if
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_setDataTypeMatrix (rmatrix, cdataType)

!<description>
  ! Converts the matrix to another data type.
!</description>

!<input>
  ! Target datatype of the matrix
  integer, intent(in) :: cdataType
!</input>

!<inputoutput>
  ! Matrix that should be converted
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>
!</subroutine>

    ! Check if matrix needs conversion
  if (rmatrix%cdataType .eq. cdataType) return

    ! Set data type
    rmatrix%cdataType = cdataType

    ! Check if matrix has data
    if (rmatrix%h_Da .ne. ST_NOHANDLE)&
        call storage_setdatatype (rmatrix%h_Da, cdataType)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_setDataTypeVector (rvector, cdataType)

!<description>
  ! Converts the vector to another data type.
!</description>

!<input>
  ! Target datatype of the vector
  integer, intent(in) :: cdataType
!</input>

!<inputoutput>
  ! Vector that should be converted
  type(t_vectorScalar), intent(inout) :: rvector
!</inputoutput>
!</subroutine>

    ! Check if vector needs conversion
    if (rvector%cdataType .eq. cdataType) return

    ! Set data type
    rvector%cdataType = cdataType

    ! Check if matrix has data
    if (rvector%h_Ddata .ne. ST_NOHANDLE)&
        call storage_setdatatype (rvector%h_Ddata, cdataType)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_getbaseVector_double (rvector,p_Ddata)

!<description>
  ! Returns a pointer to the double precision data array of the vector.
  ! An error is thrown if the vector is not double precision.
!</description>

!<input>
  ! The vector
  type(t_vectorScalar), intent(in) :: rvector
!</input>

!<output>
  ! Pointer to the double precision data array of the vector.
  ! NULL() if the vector has no data array.
  real(DP), dimension(:), pointer :: p_Ddata
!</output>

!</subroutine>

  ! Do we have data at all?
  if ((rvector%NEQ .eq. 0) .or. (rvector%h_Ddata .eq. ST_NOHANDLE)) then
    call output_line('Trying to access empty vector!', &
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbase_double')
    call sys_halt()
  end if

  ! Check that the vector is really double precision
  if (rvector%cdataType .ne. ST_DOUBLE) then
    call output_line('Vector is of wrong precision!', &
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbase_double')
    call sys_halt()
  end if

  ! Get the data array
  call storage_getbase_double (rvector%h_Ddata,p_Ddata)

  ! Modify the starting address/length to get the real array.
  p_Ddata => p_Ddata(rvector%iidxFirstEntry:&
      rvector%iidxFirstEntry+rvector%NEQ*rvector%NVAR-1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_getbaseVector_single (rvector,p_Fdata)

!<description>
  ! Returns a pointer to the single precision data array of the vector.
  ! An error is thrown if the vector is not single precision.
!</description>

!<input>
  ! The vector
  type(t_vectorScalar), intent(in) :: rvector
!</input>

!<output>
  ! Pointer to the double precision data array of the vector.
  ! NULL() if the vector has no data array.
  real(SP), dimension(:), pointer :: p_Fdata
!</output>

!</subroutine>

  ! Do we have data at all?
  if ((rvector%NEQ .eq. 0) .or. (rvector%h_Ddata .eq. ST_NOHANDLE)) then
    nullify(p_Fdata)
    return
  end if

  ! Check that the vector is really double precision
  if (rvector%cdataType .ne. ST_SINGLE) then
    call output_line('Vector is of wrong precision!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbaseVector_single')
    call sys_halt()
  end if

  ! Get the data array
  call storage_getbase_single (rvector%h_Ddata,p_Fdata)

  ! Modify the starting address/length to get the real array.
  p_Fdata => p_Fdata(rvector%iidxFirstEntry:&
      rvector%iidxFirstEntry+rvector%NEQ*rvector%NVAR-1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_getbaseVector_int (rvector,p_Idata)

!<description>
  ! Returns a pointer to the integer data array of the vector.
  ! An error is thrown if the vector is not integer.
!</description>

!<input>
  ! The vector
  type(t_vectorScalar), intent(in) :: rvector
!</input>

!<output>
  ! Pointer to the integer data array of the vector.
  ! NULL() if the vector has no data array.
  integer, dimension(:), pointer :: p_Idata
!</output>

!</subroutine>

  ! Do we have data at all?
  if ((rvector%NEQ .eq. 0) .or. (rvector%h_Ddata .eq. ST_NOHANDLE)) then
    nullify(p_Idata)
    return
  end if

  ! Check that the vector is really integer
  if (rvector%cdataType .ne. ST_INT) then
    call output_line('Vector is of wrong precision!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbaseVector_int')
    call sys_halt()
  end if

  ! Get the data array
  call storage_getbase_int (rvector%h_Ddata,p_Idata)

  ! Modify the starting address/length to get the real array.
  p_Idata => p_Idata(rvector%iidxFirstEntry:&
      rvector%iidxFirstEntry+rvector%NEQ*rvector%NVAR-1)

  end subroutine

  ! ***************************************************************************

!<function>

  logical function lsyssc_isExplicitMatrix1D (rmatrix)

!<description>
  ! Checks whether a given matrix explicitly exists in memory as 1D array,
  ! which can be accessed via lsyssc_getbase_double. This is the usual case
  ! for most of the matrices (format-7, format-9,...) but might be different
  ! for special-type matrices that do not exist in memory, i.e. where it is
  ! only known how to apply them to a vector.
!</description>

!<input>
    ! The matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<result>
    ! TRUE if the matrix is an explicit matrix that can be accessed with
    ! lsyssc_getbase_double.
    ! FALSE if the matrix is not accessible in that way.
!</result>

!</function>

    ! Check the matrix type. If we have a matrix that is known to be in memory,
    ! we can return TRUE.
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX1, LSYSSC_MATRIX9,     LSYSSC_MATRIX7, &
          LSYSSC_MATRIXD, LSYSSC_MATRIX7INTL, LSYSSC_MATRIX9INTL,&
          LSYSSC_MATRIX9ROWC)
      lsyssc_isExplicitMatrix1D = .true.

    case default
      lsyssc_isExplicitMatrix1D = .false.

    end select

  end function

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_getbaseMatrixDa_double (rmatrix, p_Ddata)

!<description>
    ! Returns a pointer to the double precision data array of the matrix.
    ! An error is thrown if the matrix is not double precision.
!</description>

!<input>
    ! The matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<output>
    ! Pointer to the double precision data array of the matrix.
    ! NULL() if the matrix has no data array.
    real(DP), dimension(:), pointer :: p_Ddata
!</output>

!</subroutine>

    ! Check if the matrix format allows to access Da
    if (.not. lsyssc_isExplicitMatrix1D(rmatrix)) then
      call output_line('Matrix does not exist explicitely!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbaseMatrixDa_double')
      call sys_halt()
    end if

    ! Do we have data at all?
    if ((rmatrix%NA .eq. 0) .or. (rmatrix%h_Da .eq. ST_NOHANDLE)) then
      nullify(p_Ddata)
      return
    end if

    ! Check that the matrix is really double precision
    if (rmatrix%cdataType .ne. ST_DOUBLE) then
      call output_line('Matrix is of wrong precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbaseMatrixDa_double')
      call sys_halt()
    end if

    ! Get the data array
    select case (rmatrix%cinterleavematrixFormat)
    case (LSYSSC_MATRIX1)
      call storage_getbase_double (rmatrix%h_Da,p_Ddata,rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR)
    case (LSYSSC_MATRIXD)
      call storage_getbase_double (rmatrix%h_Da,p_Ddata,rmatrix%NA*rmatrix%NVAR)
    case default
      call storage_getbase_double (rmatrix%h_Da,p_Ddata,rmatrix%NA)
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_getbaseMatrixDa_single (rmatrix, p_Fdata)

!<description>
    ! Returns a pointer to the single precision data array of the matrix.
    ! An error is thrown if the matrix is not single precision.
!</description>

!<input>
    ! The matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<output>
    ! Pointer to the single precision data array of the matrix.
    ! NULL() if the matrix has no data array.
    real(SP), dimension(:), pointer :: p_Fdata
!</output>

!</subroutine>

    ! Check if the matrix format allows to access Da
    if (.not. lsyssc_isExplicitMatrix1D(rmatrix)) then
      call output_line('Matrix does not exist explicitely!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbaseMatrixDa_single')
      call sys_halt()
    end if

    ! Do we have data at all?
    if ((rmatrix%NA .eq. 0) .or. (rmatrix%h_Da .eq. ST_NOHANDLE)) then
      nullify(p_Fdata)
      return
    end if

    ! Check that the matrix is really single precision
    if (rmatrix%cdataType .ne. ST_SINGLE) then
      call output_line('Matrix is of wrong precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbaseMatrixDa_single')
      call sys_halt()
    end if

    ! Get the data array
    select case (rmatrix%cinterleavematrixFormat)
    case (LSYSSC_MATRIX1)
      call storage_getbase_single (rmatrix%h_Da,p_Fdata,rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR)
    case (LSYSSC_MATRIXD)
      call storage_getbase_single (rmatrix%h_Da,p_Fdata,rmatrix%NA*rmatrix%NVAR)
    case default
      call storage_getbase_single (rmatrix%h_Da,p_Fdata,rmatrix%NA)
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_getbaseMatrixDa_int (rmatrix, p_Idata)

!<description>
    ! Returns a pointer to the integer data array of the matrix.
    ! An error is thrown if the matrix is not integer.
!</description>

!<input>
    ! The matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<output>
    ! Pointer to the integer data array of the matrix.
    ! NULL() if the matrix has no data array.
    integer, dimension(:), pointer :: p_Idata
!</output>

!</subroutine>

    ! Check if the matrix format allows to access Da
    if (.not. lsyssc_isExplicitMatrix1D(rmatrix)) then
      call output_line('Matrix does not exist explicitely!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbaseMatrixDa_int')
      call sys_halt()
    end if

    ! Do we have data at all?
    if ((rmatrix%NA .eq. 0) .or. (rmatrix%h_Da .eq. ST_NOHANDLE)) then
      nullify(p_Idata)
      return
    end if

    ! Check that the matrix is really integer
    if (rmatrix%cdataType .ne. ST_INT) then
      call output_line('Matrix is of wrong precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbaseMatrixDa_int')
      call sys_halt()
    end if

    ! Get the data array
    select case (rmatrix%cinterleavematrixFormat)
    case (LSYSSC_MATRIX1)
      call storage_getbase_int (rmatrix%h_Da,p_Idata,rmatrix%NA*rmatrix%NVAR)
    case (LSYSSC_MATRIXD)
      call storage_getbase_int (rmatrix%h_Da,p_Idata,rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR)
    case default
      call storage_getbase_int (rmatrix%h_Da,p_Idata,rmatrix%NA)
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_getbase_Kcol (rmatrix,p_Kcol)

!<description>
    ! Returns a pointer to the column data array of the matrix.
    ! An error is thrown if the matrix does not provide column array.
!</description>

!<input>
    ! The matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<output>
    ! Pointer to the column array of the matrix.
    ! NULL() if the matrix has no column array.
    integer, dimension(:), pointer :: p_Kcol
!</output>

!</subroutine>

    ! Is matrix in correct format?
    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7INTL) .and.&
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9ROWC) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9INTL)) then
      call output_line('Matrix format does not provide KCOL!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbase_Kcol')
      call sys_halt()
    end if

    ! Do we have a column array at all?
    if ((rmatrix%NA .eq. 0) .or. (rmatrix%h_Kcol .eq. ST_NOHANDLE)) then
      nullify(p_Kcol)
      return
    end if

    ! Get the column array
    call storage_getbase_int (rmatrix%h_Kcol,p_Kcol,rmatrix%NA)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_getbase_Kld (rmatrix,p_Kld)

!<description>
    ! Returns a pointer to the column offset data array of the matrix.
    ! An error is thrown if the matrix does not provide column offset array.
!</description>

!<input>
    ! The matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<output>
    ! Pointer to the column offset array of the matrix.
    ! NULL() if the matrix has no column offset array.
    integer, dimension(:), pointer :: p_Kld
!</output>

!</subroutine>

    ! Is matrix in correct format?
    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7INTL) .and.&
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9ROWC) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9INTL)) then
      call output_line('Matrix format does not provide KLD!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbase_Kld')
      call sys_halt()
    end if

    ! Do we have a column offset array at all?
    if ((rmatrix%NEQ .eq. 0) .or. (rmatrix%h_Kld .eq. ST_NOHANDLE)) then
      nullify(p_Kld)
      return
    end if

    ! Get the column offset array.
    ! Take care: If the matrix is virtually transposed, NCOLS and NEQ are
    ! exchanged!
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9ROWC)
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .eq. 0) then
        call storage_getbase_int (rmatrix%h_Kld,p_Kld,rmatrix%nnzrows+1)
      else
        call storage_getbase_int (rmatrix%h_Kld,p_Kld,rmatrix%NCOLS+1)
      end if
    case default
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .eq. 0) then
        call storage_getbase_int (rmatrix%h_Kld,p_Kld,rmatrix%NEQ+1)
      else
        call storage_getbase_int (rmatrix%h_Kld,p_Kld,rmatrix%NCOLS+1)
      end if
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_getbase_KrowIdx (rmatrix,p_KrowIdx)

!<description>
    ! Returns a pointer to the data array specifying the nonzero rows in
    ! the matrix. An error is thrown if the matrix does not provide column
    ! offset array.
!</description>

!<input>
    ! The matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<output>
    ! Pointer to the nonzero-rows array of the matrix.
    ! NULL() if the matrix has no column offset array.
    integer, dimension(:), pointer :: p_KrowIdx
!</output>

!</subroutine>

    ! Is matrix in correct format?
    if (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9ROWC)  then
      call output_line('Matrix format does not provide KNZROWS!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbase_KrowIdx')
      call sys_halt()
    end if

    ! Do we have a column offset array at all?
    if ((rmatrix%NEQ .eq. 0) .or. (rmatrix%h_KrowIdx .eq. ST_NOHANDLE)) then
      nullify(p_KrowIdx)
      return
    end if

    ! Get the column offset array.
    call storage_getbase_int (rmatrix%h_KrowIdx,p_KrowIdx)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)

!<description>
    ! Returns a pointer to the diagonal data array of the matrix.
    ! An error is thrown if the matrix does not provide diagonal array.
!</description>

!<input>
    ! The matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<output>
    ! Pointer to the diagonal array of the matrix.
    ! NULL() if the matrix has no diagonal array.
    integer, dimension(:), pointer :: p_Kdiagonal
!</output>

!</subroutine>

    ! Is matrix in correct format?
    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and.&
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9INTL)) then
      call output_line('Matrix format does not provide KDIAGONAL!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_getbase_Kdiagonal')
      call sys_halt()
    end if

    ! Do we have a column offset array at all?
    if ((rmatrix%NEQ .eq. 0) .or. (rmatrix%h_Kdiagonal .eq. ST_NOHANDLE)) then
      nullify(p_Kdiagonal)
      return
    end if

    ! Get the column offset array
    call storage_getbase_int (rmatrix%h_Kdiagonal,p_Kdiagonal,rmatrix%NEQ)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_createVecDirect (rvector,NEQ,bclear,cdataType,NEQMAX)

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
  integer, intent(in) :: NEQ

  ! OPTIONAL: Whether to fill the vector with zero initially.
  ! If not specified, the vector is left uninitialised.
  logical, intent(in), optional :: bclear

  ! OPTIONAL: Data type of the vector.
  ! If not specified, ST_DOUBLE is assumed.
  integer, intent(in), optional :: cdataType

  ! OPTIONAL: Maximum length of the vector
  integer, intent(in), optional :: NEQMAX

!</input>

!<output>
  ! Scalar vector structure
  type(t_vectorScalar), intent(out) :: rvector
!</output>

!</subroutine>

  integer :: cdata
  integer :: isize
  logical :: bclr

    cdata = ST_DOUBLE
    if (present(cdataType)) cdata = cdataType
    isize = max(0,NEQ)
    if (present(NEQMAX)) isize=max(isize,NEQMAX)

    ! The INTENT(out) already initialises rvector with the most important
    ! information. The rest comes now:

    ! Datatype
    rvector%cdataType = cdata

    ! Size:
    rvector%NEQ = max(0,NEQ)

    ! Handle - if NEQ > 0
    if (rvector%NEQ .gt. 0) then
    
      bclr = .false.
      if (present(bclear)) bclr = bclear
    
      if (bclr) then
        call storage_new ('lsyssc_createVecDirect', 'ScalarVector', isize, &
                            cdata, rvector%h_Ddata, ST_NEWBLOCK_ZERO)
      else
        call storage_new ('lsyssc_createVecDirect', 'ScalarVector', isize, &
                            cdata, rvector%h_Ddata, ST_NEWBLOCK_NOINIT)
      end if
    end if

  end subroutine

!****************************************************************************

!<subroutine>

  subroutine lsyssc_createVecDirectIntl (rvector,NEQ,NVAR,bclear,cdataType,NEQMAX)

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
  integer, intent(in) :: NEQ

  ! Desired number of local variables
  integer, intent(in) :: NVAR

  ! Whether to fill the vector with zero initially
  logical, intent(in) :: bclear

  ! OPTIONAL: Data type of the vector.
  ! If not specified, ST_DOUBLE is assumed.
  integer, intent(in), optional :: cdataType

  ! OPTIONAL: Maximum length of the vector
  integer, intent(in), optional :: NEQMAX

!</input>

!<output>
  ! Scalar vector structure
  type(t_vectorScalar), intent(out) :: rvector
!</output>

!</subroutine>

  integer :: cdata
  integer :: isize

    cdata = ST_DOUBLE
    if (present(cdataType)) cdata = cdataType
    isize = max(0,NEQ)
    if (present(NEQMAX)) isize=max(isize,NEQMAX)
    isize = isize*max(0,NVAR)

    ! The INTENT(out) already initialises rvector with the most important
    ! information. The rest comes now:

    ! Datatype
    rvector%cdataType = cdata

    ! Size:
    rvector%NEQ = max(0,NEQ)
    rvector%NVAR= max(0,NVAR)

    ! Handle - if NEQ > 0
    if (rvector%NEQ*rvector%NVAR .gt. 0) then
      if (bclear) then
        call storage_new ('lsyssc_createVector', 'ScalarVector', isize, &
                            cdata, rvector%h_Ddata,ST_NEWBLOCK_ZERO)
      else
        call storage_new ('lsyssc_createVector', 'ScalarVector', isize, &
                            cdata, rvector%h_Ddata,ST_NEWBLOCK_NOINIT)
      end if
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_createVecByDiscr (rdiscretisation,rx,bclear,cdataType,NEQMAX)

!<description>
  ! Initialises the vector structure rx based on a discretisation
  ! structure rDiscretisation.
  !
  ! Memory is allocated on the heap for rx accordint to the number of
  ! DOF`s indicated by the spatial discretisation structures in
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
  ! A discretisation structure specifying the spatial discretisations
  ! for all the subblocks in rx.
  type(t_spatialDiscretisation),intent(in), target :: rdiscretisation

  ! OPTIONAL: If set to YES, the vector will be filled with zero initially.
  ! Otherwise the content of rx is undefined.
  logical, intent(in), optional :: bclear

  ! OPTIONAL: Data type identifier for the entries in the vector.
  ! Either ST_SINGLE or ST_DOUBLE. If not present, ST_DOUBLE is used.
  integer, intent(in), optional :: cdataType

  ! OPTIONAL: Maximum length of the vector
  integer, intent(in), optional :: NEQMAX
!</input>

!<output>
  ! Destination structure. Memory is allocated appropriately.
  ! A pointer to rdiscretisation is saved to r.
  type(t_vectorScalar),intent(out) :: rx
!</output>

!</subroutine>

  integer :: NEQ
  logical :: bcl

  bcl = .false.
  if (present(bclear)) bcl = bclear

  ! Get NEQ:
  NEQ = dof_igetNDofGlob(rdiscretisation)

  ! Create a new vector with that block structure
  call lsyssc_createVector (rx, NEQ, bcl, cdataType, NEQMAX)

  ! Initialise further data of the block vector
  rx%p_rspatialDiscr => rdiscretisation

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_createVecByDiscrIntl (rdiscretisation,rx,NVAR,bclear,&
                                          cdataType,NEQMAX)

!<description>
  ! Initialises the vector structure rx based on a discretisation
  ! structure rDiscretisation.
  !
  ! Memory is allocated on the heap for rx accordint to the number of
  ! DOF`s indicated by the spatial discretisation structures in
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
  ! A discretisation structure specifying the spatial discretisations
  ! for all the subblocks in rx.
  type(t_spatialDiscretisation),intent(in), target :: rdiscretisation

  ! Desired number of local variables
  integer, intent(in) :: NVAR

  ! OPTIONAL: If set to YES, the vector will be filled with zero initially.
  ! Otherwise the content of rx is undefined.
  logical, intent(in), optional :: bclear

  ! OPTIONAL: Data type identifier for the entries in the vector.
  ! Either ST_SINGLE or ST_DOUBLE. If not present, ST_DOUBLE is used.
  integer, intent(in), optional :: cdataType

  ! OPTIONAL: Maximum length of the vector
  integer, intent(in), optional :: NEQMAX
!</input>

!<output>
  ! Destination structure. Memory is allocated appropriately.
  ! A pointer to rdiscretisation is saved to r.
  type(t_vectorScalar),intent(out) :: rx
!</output>

!</subroutine>

  integer :: NEQ
  logical :: bcl

  bcl = .false.
  if (present(bclear)) bcl = bclear

  ! Get NEQ:
  NEQ = dof_igetNDofGlob(rdiscretisation)

  ! Create a new vector with that block structure
  call lsyssc_createVector (rx, NEQ, NVAR, bcl, cdataType, NEQMAX)

  ! Initialise further data of the block vector
  rx%p_rspatialDiscr => rdiscretisation

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_createVecIndMat (rtemplateMat,rx,bclear,btransposed,&
                                     cdataType,NEQMAX)

!<description>
    ! Initialises the scalar vector structure rx. rtemplateMat is an
    ! existing scalar matrix structure. The vector rx will be created
    ! according to the size of the template matrix.
    !
    ! Memory is allocated on the heap for rx. The vector rx will have
    ! the same size as the number of matrix columns.
    ! The sorting strategy of the vector is initialised with the
    ! sorting strategy of the template matrix rtemplateMat.
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
    ! A template matrix structure
    type(t_matrixScalar), intent(in) :: rtemplateMat

    ! OPTIONAL: If set to TRUE, the vector will be filled with zero initially.
    ! Otherwise the content of rx is undefined.
    logical, intent(in), optional :: bclear

    ! OPTIONAL: If not specified or set to FALSE, the vector will be
    ! created as 'right' vector of the matrix (so matrix vector multiplication
    ! (A x) is possible).
    ! If set to TRUE, the vector will be a 'left' vector, i.e. matrix
    ! vector multiplication (A^T x) is possible.
    logical, intent(in), optional :: btransposed

    ! OPTIONAL: Data type identifier for the entries in the vector.
    ! Either ST_SINGLE or ST_DOUBLE. If not present, ST_DOUBLE is assumed.
    integer, intent(in), optional :: cdataType

    ! OPTIONAL: Maximum length of the vector
    integer, intent(in), optional :: NEQMAX
!</input>

!<output>
    ! Destination structure. Memory is allocated for each of the blocks.
    type(t_vectorScalar),intent(out) :: rx
!</output>

!</subroutine>

    ! local variables
    integer :: cdata
    integer :: NEQ,NCOLS,isize

    cdata = ST_DOUBLE
    if (present(cdataType)) cdata = cdataType

    ! Transfer discretization pointers from the matrix to the vector
    rx%p_rspatialDiscr => rtemplateMat%p_rspatialDiscrTrial

    NEQ   = rtemplateMat%NEQ
    NCOLS = rtemplateMat%NCOLS

    if (present(btransposed)) then
      if (btransposed) then
        NEQ   = rtemplateMat%NCOLS
        NCOLS = rtemplateMat%NEQ

        ! Transfer discretization pointers from the matrix to the vector
        rx%p_rspatialDiscr => rtemplateMat%p_rspatialDiscrTest
      end if
    end if

    isize = max(0,NCOLS)
    if (present(NEQMAX)) isize=max(isize,NEQMAX)
    isize = isize*max(0,rtemplateMat%NVAR)

    ! Allocate memory for vector
    call storage_new ('lsyssc_createVecIndMat', 'Vector', &
                        isize, cdata, rx%h_Ddata, ST_NEWBLOCK_NOINIT)

    ! Set structure
    rx%NEQ       = NCOLS
    rx%NVAR      = rtemplateMat%NVAR
    rx%cdataType = cdata

    ! Transfer sorting strategy from the matrix to the vector
    rx%bisSorted = rtemplateMat%bcolumnsSorted
    rx%p_rsortStrategy => rtemplateMat%p_rsortStrategyColumns

    if (present(bclear)) then
      if (bclear) then
        call lsyssc_clearVector (rx)
      end if
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_resizeVectorDirect (rvector, NEQ, bclear, bcopy, NEQMAX)

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
  ! vector should be cleared afterwards ;-)
  !
  ! If the optional parameter NEQMAX is specified, then memory of size
  ! NEQMAX is allocated if reallocate is required. Otherwise, no reallocation
  ! is performed and NEQMAX is just neglected.
!</description>

!<input>

    ! Desired length of the vector
    integer, intent(in) :: NEQ

    ! Whether to fill the vector with zero initially
    logical, intent(in) :: bclear

    ! OPTIONAL: Whether to copy the content of the vector to the resized one
    logical, intent(in), optional :: bcopy

    ! OPTIONAL: Maximum length of the vector
    integer, intent(in), optional :: NEQMAX

!</input>

!<inputoutput>

    ! Scalar vector structure
    type(t_vectorScalar), intent(inout) :: rvector

!</inputoutput>

!</subroutine>

    integer :: iNEQ,isize
    logical :: bdocopy

    ! Check, that vector is not a copy of another (possibly larger) vector
    if (rvector%bisCopy) then
      call output_line('A copied vector cannot be resized!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeVectorDirect')
      call sys_halt()
    end if

    ! Check, if vector has been initialised before.
    if (rvector%NEQ .eq. 0 .or. rvector%h_Ddata .eq. ST_NOHANDLE) then
      call output_line('A vector can only be resized if it has been created correctly!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeVectorDirect')
      call sys_halt()
    end if

    ! Set working dimensions
    iNEQ = max(0,NEQ)
    if (present(NEQMAX)) iNEQ = max(iNEQ,NEQMAX)

    ! Set copy/clear attributes
    bdocopy = (.not.bclear)
    if (present(bcopy)) bdocopy = (bdocopy .and. bcopy)

    ! If the vector should be cleared, then the sorting strategy (if any)
    ! can be ignored and reset. Otherwise, the vector needs to be unsorted
    ! prior to copying some part of it. Afterwards, no sorting strategy is
    ! available in any case.
    if (bdocopy .and. rvector%bisSorted) then
      call lsyssc_sortVector(rvector,.false.)
    end if

    ! Reset sorting strategy, there is none
    rvector%bisSorted = .false.
    nullify(rvector%p_rsortStrategy)

    ! Get current size of vector memory
    call storage_getsize(rvector%h_Ddata, isize)

    ! Update NEQ. Note that NVAR cannot be modified.
    rvector%NEQ = max(0,NEQ)

    ! Do we really have to reallocate the vector physically?
    if (rvector%NEQ*rvector%NVAR > isize) then

      ! Yes, so adopt the new size. Note that some extra memory is
      ! allocated if the optional argument NEQMAX is specified
      isize = iNEQ*rvector%NVAR

      ! Reallocate the memory for handle h_Ddata
      call storage_realloc('lsyssc_resizeVectorDirect',&
          isize, rvector%h_Ddata, ST_NEWBLOCK_NOINIT, bdocopy)

    elseif (present(NEQMAX)) then

      ! The available memory suffices for the vector, i.e. isize <= NVAR*NEQ.
      ! Let us check if the user supplied a new upper limit which makes it
      ! mandatory to "shrink" the allocated memory. Note that memory for
      ! at least NEQ vector intries is allocated.
      if (isize > iNEQ*rvector%NVAR) then

        ! Compute new size, i.e. MAX(0,NEQ,NEQMAX)
        isize = iNEQ*rvector%NVAR

        if (isize .eq. 0) then
          ! If nothing is left, then the vector can also be released.
          call lsyssc_releaseVector(rvector)
          return
        else
          call storage_realloc('lsyssc_resizeVectorDirect',&
              isize, rvector%h_Ddata, ST_NEWBLOCK_NOINIT, bdocopy)
        end if
      end if
    end if

    ! Should the vector be cleared?
    if (bclear) call lsyssc_clearVector(rvector)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_resizeVectorByDiscr(rdiscretisation, rx,  bclear, bcopy, NEQMAX)

!<description>
  ! Resize the vector structure rx based on a discretisation
  ! structure rDiscretisation.
  !
  ! Memory is reallocated on the heap for rx accordint to the number of
  ! DOF`s indicated by the spatial discretisation structures in
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
  ! A discretisation structure specifying the spatial discretisations
  ! for all the subblocks in rx.
  type(t_spatialDiscretisation),intent(in), target :: rdiscretisation

  ! Whether to fill the vector with zero initially
  logical, intent(in) :: bclear

  ! OPTIONAL: Whether to copy the content of the vector to the resized one
  logical, intent(in), optional :: bcopy

  ! OPTIONAL: Maximum length of the vector
  integer, intent(in), optional :: NEQMAX

!</input>

!<inputoutput>

  ! Scalar vector structure
  type(t_vectorScalar), intent(inout) :: rx

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: NEQ

    ! Get NEQ:
    NEQ = dof_igetNDofGlob(rdiscretisation)

    ! Resize scalar vector
    call lsyssc_resizeVectorDirect (rx, NEQ, bclear, bcopy, NEQMAX)

    ! Attach spatial discretisation structure
    rx%p_rspatialDiscr => rdiscretisation


  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_resizeVectorIndirect (rvector, rvectorTemplate, bclear,&
                                          bcopy, NEQMAX)

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
    type(t_vectorScalar), intent(in) :: rvectorTemplate

    ! Whether to fill the vector with zero initially
    logical, intent(in) :: bclear

    ! OPTIONAL: Whether to copy the content of the vector to the resized one
    logical, intent(in), optional :: bcopy

    ! OPTIONAL: Maximum length of the vector
    integer, intent(in), optional :: NEQMAX

!</input>

!<inputoutput>

    ! Scalar vector structure
    type(t_vectorScalar), intent(inout) :: rvector

!</inputoutput>

!</subroutine>

    integer :: isize,iNEQMAX

    ! Check, if vector is a copy of another (possibly larger) vector
    if (rvector%bisCopy) then
      call output_line('A copied vector cannot be resized!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeVectorIndirect')
      call sys_halt()
    end if

    ! Check, if vector has been initialised before. If this is not
    ! the case, then create a new vector by duplicating the template vector.
    ! Moreover, no data is copied and the vector is cleared if bclear=.TRUE.
    if ((rvector%NEQ     .eq. 0) .or.&
        (rvector%h_Ddata .eq. ST_NOHANDLE)) then

      ! At first, copy all 'local' data.
      rvector = rvectorTemplate

      ! Then allocate a new array for the content in the same data
      ! format as the template vector.
      call storage_getsize(rvectorTemplate%h_Ddata, isize)
      if (present(NEQMAX)) isize = max(isize,NEQMAX)

      if (bclear) then
        call storage_new ('lsyssc_resizeVectorIndirect','ScalarVector',isize,&
            rvectorTemplate%cdataType, rvector%h_Ddata, ST_NEWBLOCK_ZERO)
      else
        call storage_new ('lsyssc_resizeVectorIndirect','ScalarVector',isize,&
            rvectorTemplate%cdataType, rvector%h_Ddata, ST_NEWBLOCK_NOINIT)
      end if

    else

      ! Check if vectors are compatible
      if ((rvector%cdataType .ne. rvectorTemplate%cdataType) .or.&
          (rvector%NVAR      .ne. rvectorTemplate%NVAR)) then
        call output_line('Vectors are incompatible!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeVectorIndirect')
        call sys_halt()
      end if

      ! Get NEQMAX from template vector
      call storage_getsize(rvectorTemplate%h_Ddata, iNEQMAX)
      if (present(NEQMAX)) iNEQMAX = max(iNEQMAX,NEQMAX)

      ! Resize vector directly
      call lsyssc_resizeVectorDirect(rvector, rvectorTemplate%NEQ, bclear, bcopy, iNEQMAX)

    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_resizeVectorIndMat (rmatrix, rvector, bclear, bcopy, NEQMAX)

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
    type(t_matrixScalar), intent(in) :: rmatrix

    ! Whether to fill the vector with zero initially
    logical, intent(in) :: bclear

    ! OPTIONAL: Whether to copy the content of the vector to the resized one
    logical, intent(in), optional :: bcopy

    ! OPTIONAL: Maximum length of the vector
    integer, intent(in), optional :: NEQMAX

!</input>

!<inputoutput>

    ! Scalar vector structure
    type(t_vectorScalar), intent(inout) :: rvector

!</inputoutput>

!</subroutine>

    ! Check, if vector is a copy of another (possibly larger) vector
    if (rvector%bisCopy) then
      call output_line('A copied vector cannot be resized!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeVectorIndMat')
      call sys_halt()
    end if

    ! Does the vector exist?
    if (rvector%NEQ     .eq. 0 .or.&
        rvector%h_Ddata .eq. ST_NOHANDLE) then

      ! Create new vector
      call lsyssc_createVecIndMat(rmatrix, rvector, bclear)

    else

      ! Check if vector/matrix are compatible
      if (rvector%NVAR .ne. rmatrix%NVAR) then
        call output_line('Vector/Matrix incompatible!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeVectorIndMat')
        call sys_halt()
      end if

      ! Resize vector directly
      call lsyssc_resizeVectorDirect(rvector, rmatrix%NCOLS, bclear, bcopy, NEQMAX)

    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_resizeMatrixDirect (rmatrix, NEQ, NCOLS, NA, bclear, bcopy, &
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
    integer, intent(in) :: NEQ

    ! Desired number of columns
    integer, intent(in) :: NCOLS

    ! Desired number of elements
    integer, intent(in) :: NA

    ! Whether to fill the matrix with zero initially
    logical, intent(in) :: bclear

    ! OPTIONAL: Whether to copy the content of the matrix to the resized one
    logical, intent(in), optional :: bcopy

    ! OPTIONAL: Maximum number of equations
    integer, intent(in), optional :: NEQMAX

    ! OPTIONAL: Maximum number of columns
    integer, intent(in), optional :: NCOLSMAX

    ! OPTIONAL: Maximum number of elements
    integer, intent(in), optional :: NAMAX

    ! OPTIONAL: Wether to enforce resize even if matrix is copied from another matrix
    logical, intent(in), optional :: bforce

!</input>

!<inputoutput>

    ! Scalar matrix structure
    type(t_matrixScalar), intent(inout) :: rmatrix

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iNA,isize,isizeNew
    integer :: iNEQ,iNCOLS
    logical :: bdocopy,bdoresize,btransposed

    ! Check if resize should be forced
    bdoresize = .false.
    if (present(bforce)) bdoresize=bforce

    ! Check, if matrix is not a copy of another matrix or if resize is to be enforced
    if ((iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0 .or.&
         iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY)   .ne. 0) .and.&
         .not.bdoresize) then
      call output_line('A copied matrix can only be resized if forced explicitely!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixDirect')
      call sys_halt()
    end if

    ! Check, if matrix has been initialised before.
    if (rmatrix%NEQ .eq. 0 .or. rmatrix%NCOLS .eq. 0 .or. rmatrix%NA .eq. 0) then
      call output_line('A matrix can only be resized if it has been created correctly!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixDirect')
      call sys_halt()
    end if

    ! Set working dimensions
    iNA    = max(0,NA)
    if (present(NAMAX))    iNA    = max(iNA,NAMAX)
    iNEQ   = max(0,NEQ)
    if (present(NEQMAX))   iNEQ   = max(iNEQ,NEQMAX)
    iNCOLS = max(0,NCOLS)
    if (present(NCOLSMAX)) iNCOLS = max(iNCOLS,NCOLSMAX)

    ! Set copy/clear attributes
    bdocopy = (.not.bclear)
    if (present(bcopy)) bdocopy = (bdocopy .and. bcopy)

    ! Set transposed indicator
    btransposed = (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0)

    ! If the matrix should be cleared, then the sorting strategy (if any)
    ! can be ignored and reset. Otherwise, the matrix needs to be unsorted
    ! prior to copying some part of it. Afterwards, no sorting strategy is
    ! available in any case.
    if (bdocopy .and. (rmatrix%bcolumnsSorted .or. rmatrix%browsSorted)) then
      call lsyssc_sortMatrix(rmatrix,.false.)
    end if

    ! Reset sorting strategy, there is none
    rmatrix%bcolumnsSorted = .false.
    rmatrix%browsSorted = .false.
    nullify(rmatrix%p_rsortStrategyColumns)
    nullify(rmatrix%p_rsortStrategyRows)

    ! Update NEQ, NCOLS and NA.
    rmatrix%NA    = max(0,NA)
    rmatrix%NEQ   = max(0,NEQ)
    rmatrix%NCOLS = max(0,NCOLS)

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
    select case (rmatrix%cmatrixFormat)

    case (LSYSSC_MATRIX1,LSYSSC_MATRIXD)

      ! For these matrices, only the data array needs to be resized since
      ! there is no actual matrix structure.

      ! Are we copy or not?
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .ne. 0) then

        ! Check if handle coincides with matrix dimensions
        call storage_getsize(rmatrix%h_Da,isize)
        if (rmatrix%NA*rmatrix%NVAR > isize) then
          call output_line('Dimensions of copied matrix mismatch!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixDirect')
          call sys_halt()
        end if

      else

        ! Do we really have to reallocate the matrix data physically?
        if (rmatrix%h_Da .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrix%h_Da,isize)
          if (rmatrix%NA*rmatrix%NVAR > isize) then

            ! Yes, so adopt the new size or reserve some extra memory if required.
            isize = iNA*rmatrix%NVAR

            ! Reallocate the memory
            call storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Da,&
                ST_NEWBLOCK_NOINIT, bdocopy)

          elseif (present(NAMAX)) then

            ! The available memory suffices for the matrix. Let us check if the user
            ! suplied a new upper limit which makes it mandatory to "shrink" the
            ! allocated memory. Note that memory for at least NA matrix entries
            ! is allocated.
            if (isize > iNA*rmatrix%NVAR) then

              ! Set new size
              isize = iNA*rmatrix%NVAR

              if (isize .eq. 0) then
                ! If nothing is left, then the matrix can also be released.
                call lsyssc_releaseMatrix(rmatrix)
                return
              else
                call storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Da,&
                    ST_NEWBLOCK_NOINIT, bdocopy)
              end if
            end if
          end if

          ! Should the matrix be cleared?
          if (bclear) call storage_clear(rmatrix%h_Da)
        end if
      end if

    case (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)

      ! First, resize the structure of the matrix, i.e., KLD, KCOL and
      ! possibly KDIAGONAL for matrices stored in matrix format 9. Here,
      ! it is necessary to distinguish between virtually transposed matrices
      ! and non-transposed ones.

      ! Are we copy or not?
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0) then

        ! Check if handle coincides with matrix dimensions
        call storage_getsize(rmatrix%h_Kld,isize)
        if ((rmatrix%NEQ+1 > isize) .or. btransposed .and. (rmatrix%NCOLS+1 > isize)) then
          call output_line('Dimensions of copied matrix mismatch!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixDirect')
          call sys_halt()
        end if

        call storage_getsize(rmatrix%h_Kcol,isize)
        if (rmatrix%NA > isize) then
          call output_line('Dimensions of copied matrix mismatch!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixDirect')
          call sys_halt()
        end if

        if ((rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX9) .or.&
            (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX9INTL)) then
          call storage_getsize(rmatrix%h_Kdiagonal,isize)
          if ((rmatrix%NEQ > isize) .or. btransposed .and. (rmatrix%NCOLS > isize)) then
            call output_line('Dimensions of copied matrix mismatch!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixDirect')
            call sys_halt()
          end if
        end if

      else   ! The structure of the matrix is not a copy of another matrix

        ! Do we really have to reallocate the array KLD?
        if (rmatrix%h_Kld .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrix%h_Kld,isize)
          if ((rmatrix%NEQ+1 > isize) .or. btransposed .and. (rmatrix%NCOLS+1 > isize)) then

            ! Yes, so adopt the new size or reserve some extra memory if required.
            isize = merge(iNCOLS+1,iNEQ+1,btransposed)

            ! Reallocate the memory
            call storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Kld,&
                ST_NEWBLOCK_NOINIT, bdocopy)

          elseif (present(NEQMAX) .or. btransposed .and. present(NCOLSMAX)) then

            isizeNew = merge(iNCOLS+1,iNEQ+1,btransposed)

            ! The available memory suffices for the matrix. Let us check if the user
            ! suplied a new upper limit which makes it mandatory to "shrink" the
            ! allocated memory.
            if (isize > isizeNew) then

              if (isizeNew .eq. 0) then
                ! If nothing is left, then the matrix can also be released.
                call lsyssc_releaseMatrix(rmatrix)
                return
              else
                call storage_realloc('lsyssc_resizeMatrixDirect', isizeNew, rmatrix%h_Kld,&
                    ST_NEWBLOCK_NOINIT, bdocopy)
              end if
            end if
          end if

          ! Should the matrix be cleared?
          if (bclear) call storage_clear(rmatrix%h_Kld)
        end if

        ! Do we really have to reallocate the array KCOL?
        if (rmatrix%h_Kcol .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrix%h_Kcol,isize)
          if (rmatrix%NA > isize) then

            ! Yes, so adopt the new size or reserve some extra memory if required.
            isize = iNA

            ! Reallocate the memory
            call storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Kcol,&
                ST_NEWBLOCK_NOINIT, bdocopy)

          elseif (present(NAMAX)) then

            ! The available memory suffices for the matrix. Let us check if the user
            ! suplied a new upper limit which makes it mandatory to "shrink" the
            ! allocated memory.
            if (isize > iNA) then

              ! Set new size
              isize = iNA

              if (isize .eq. 0) then
                ! If nothing is left, then the matrix can also be released.
                call lsyssc_releaseMatrix(rmatrix)
                return
              else
                call storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Kcol,&
                    ST_NEWBLOCK_NOINIT, bdocopy)
              end if
            end if
          end if

          ! Should the matrix be cleared?
          if (bclear) call storage_clear(rmatrix%h_Kcol)
        end if

        ! If the matrix is stored in format 9, then the diagonal array must be resized.
        if ((rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX9) .or.&
            (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX9INTL)) then

          ! Do we really have to reallocate the array KCOL?
          if (rmatrix%h_Kdiagonal .ne. ST_NOHANDLE) then
            call storage_getsize(rmatrix%h_Kdiagonal,isize)
            if ((rmatrix%NEQ > isize) .or. btransposed .and. (rmatrix%NCOLS > isize)) then

              ! Yes, so adopt the new size or reserve some extra memory if required.
              isize = merge(iNCOLS,iNEQ,btransposed)

              ! Reallocate the memory
              call storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Kdiagonal,&
                  ST_NEWBLOCK_NOINIT, bdocopy)

            elseif (present(NEQMAX) .or. btransposed .and. present(NCOLSMAX)) then

              isizeNew = merge(iNCOLS,iNEQ,btransposed)

              ! The available memory suffices for the matrix. Let us check if the user
              ! suplied a new upper limit which makes it mandatory to "shrink" the
              ! allocated memory.
              if (isize > isizeNew) then

                if (isizeNew .eq. 0) then
                  ! If nothing is left, then the matrix can also be released.
                  call lsyssc_releaseMatrix(rmatrix)
                  return
                else
                  call storage_realloc('lsyssc_resizeMatrixDirect', isizeNew, rmatrix%h_Kdiagonal,&
                      ST_NEWBLOCK_NOINIT, bdocopy)
                end if
              end if
            end if

            ! Should the matrix be cleared?
            if (bclear) call storage_clear(rmatrix%h_Kdiagonal)
          end if
        end if

      end if

      ! Ok, the matrix structure has been resized according to the prescribed dimensions.
      ! Now, let us resize the data array of the matrix.

      ! Are we copy or not?
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .ne. 0) then

        ! Check if handle coincides with matrix dimensions
        call storage_getsize(rmatrix%h_Da,isize)

        ! What kind of interleave matrix are we (if any)?
        select case (rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIXUNDEFINED)
          if (rmatrix%NA > isize) then
            call output_line('Dimensions of copied matrix mismatch!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixDirect')
            call sys_halt()
          end if

        case (LSYSSC_MATRIX1)
          if (rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR > isize) then
            call output_line('Dimensions of copied matrix mismatch!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixDirect')
            call sys_halt()
          end if

        case (LSYSSC_MATRIXD)
          if (rmatrix%NA*rmatrix%NVAR > isize) then
            call output_line('Dimensions of copied matrix mismatch!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixDirect')
            call sys_halt()
          end if

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixDirect')
          call sys_halt()
        end select

      else   ! The content of the matrix is not a copy of another matrix

        ! What kind of interleave matrix are we (if any)?
        select case (rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIXUNDEFINED)

          ! Do we really have to reallocate the matrix physically?
          if (rmatrix%h_Da .ne. ST_NOHANDLE) then
            call storage_getsize(rmatrix%h_Da,isize)
            if (rmatrix%NA > isize) then

              ! Yes, so adopt the new size or reserve some extra memory if required.
              isize = iNA

              ! Reallocate the memory
              call storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Da,&
                  ST_NEWBLOCK_NOINIT, bdocopy)

            elseif (present(NAMAX)) then

              ! The available memory suffices for the matrix. Let us check if the user
              ! suplied a new upper limit which makes it mandatory to "shrink" the
              ! allocated memory. Note that memory for at least NA matrix entries
              ! is allocated.
              if (isize > iNA) then

                ! Set new size
                isize = iNA

                if (isize .eq. 0) then
                  ! If nothing is left, then the matrix can also be released.
                  call lsyssc_releaseMatrix(rmatrix)
                  return
                else
                  call storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Da,&
                      ST_NEWBLOCK_NOINIT, bdocopy)
                end if
              end if
            end if

            ! Should the matrix be cleared?
            if (bclear) call storage_clear(rmatrix%h_Da)
          end if

        case (LSYSSC_MATRIX1)

          ! Do we really have to reallocate the matrix physically?
          if (rmatrix%h_Da .ne. ST_NOHANDLE) then
            call storage_getsize(rmatrix%h_Da,isize)
            if (rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR > isize) then

              ! Yes, so adopt the new size or reserve some extra memory if required.
              isize = iNA*rmatrix%NVAR*rmatrix%NVAR

              ! Reallocate the memory
              call storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Da,&
                  ST_NEWBLOCK_NOINIT, bdocopy)

            elseif (present(NAMAX)) then

              ! The available memory suffices for the matrix. Let us check if the user
              ! suplied a new upper limit which makes it mandatory to "shrink" the
              ! allocated memory. Note that memory for at least NA matrix entries
              ! is allocated.
              if (isize > iNA*rmatrix%NVAR*rmatrix%NVAR) then

                ! Set new size
                isize = iNA*rmatrix%NVAR*rmatrix%NVAR

                if (isize .eq. 0) then
                  ! If nothing is left, then the matrix can also be released.
                  call lsyssc_releaseMatrix(rmatrix)
                  return
                else
                  call storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Da,&
                      ST_NEWBLOCK_NOINIT, bdocopy)
                end if
              end if
            end if

            ! Should the matrix be cleared?
            if (bclear) call storage_clear(rmatrix%h_Da)
          end if

        case (LSYSSC_MATRIXD)

          ! Do we really have to reallocate the matrix physically?
          if (rmatrix%h_Da .ne. ST_NOHANDLE) then
            call storage_getsize(rmatrix%h_Da,isize)
            if (rmatrix%NA*rmatrix%NVAR > isize) then

              ! Yes, so adopt the new size or reserve some extra memory if required.
              isize = iNA*rmatrix%NVAR

              ! Reallocate the memory
              call storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Da,&
                  ST_NEWBLOCK_NOINIT, bdocopy)

            elseif (present(NAMAX)) then

              ! The available memory suffices for the matrix. Let us check if the user
              ! suplied a new upper limit which makes it mandatory to "shrink" the
              ! allocated memory. Note that memory for at least NA matrix entries
              ! is allocated.
              if (isize > iNA*rmatrix%NVAR) then

                ! Set new size
                isize = iNA*rmatrix%NVAR

                if (isize .eq. 0) then
                  ! If nothing is left, then the matrix can also be released.
                  call lsyssc_releaseMatrix(rmatrix)
                  return
                else
                  call storage_realloc('lsyssc_resizeMatrixDirect', isize, rmatrix%h_Da,&
                      ST_NEWBLOCK_NOINIT, bdocopy)
                end if
              end if
            end if

            ! Should the matrix be cleared?
            if (bclear) call storage_clear(rmatrix%h_Da)
          end if

        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixDirect')
          call sys_halt()
        end select

      end if

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixDirect')
      call sys_halt()
    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_resizeMatrixIndirect(rmatrix, rmatrixTemplate, bclear, bcopy, &
      NEQMAX, NCOLSMAX, NAMAX, bforce)

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

    ! Scalar template matrix
    type(t_matrixScalar), intent(in) :: rmatrixTemplate

    ! Whether to fill the vector with zero initially
    logical, intent(in) :: bclear

    ! OPTIONAL: Whether to copy the content of the vector to the resized one
    logical, intent(in), optional :: bcopy

    ! OPTIONAL: Maximum number of equations
    integer, intent(in), optional :: NEQMAX

    ! OPTIONAL: Maximum number of columns
    integer, intent(in), optional :: NCOLSMAX

    ! OPTIONAL: Maximum number of elements
    integer, intent(in), optional :: NAMAX

    ! OPTIONAL: Whether to enforce resize even if matrix is copied from another matrix
    logical, intent(in), optional :: bforce

!</input>

!<inputoutput>

    ! Scalar matrix structure
    type(t_matrixScalar), intent(inout) :: rmatrix

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: isize,isizeTmp
    integer :: iNEQ,iNCOLS,iNA
    logical :: bdocopy,bdoresize,btransposed

    ! Check if resize should be forced
    bdoresize = .false.
    if (present(bforce)) bdoresize=bforce

    ! Check, if matrix is not a copy of another matrix
    if ((iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0 .or.&
         iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY)   .ne. 0) .and.&
         .not.bdoresize) then
      call output_line('A copied matrix can only be resized if this is forced explicitely!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeVectorIndMat')
      call sys_halt()
    end if

    ! Set transposed indicator
    btransposed = (iand(rmatrixTemplate%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0)

    ! Check, if matrix has been initialised before.
    if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIXUNDEFINED) then

      ! At first, copy all 'local' data.
      rmatrix = rmatrixTemplate

      ! Then allocate memory for the matrix structure and matrix data.
      if (bclear) then
        if (rmatrixTemplate%h_Kld .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrixTemplate%h_Kld, isize)
          iNEQ = isize
          if (present(NEQMAX)) iNEQ=max(iNEQ, NEQMAX+1)
          iNCOLS = isize
          if (present(NCOLSMAX)) iNCOLS=max(iNCOLS, NCOLSMAX+1)
          isize = merge(iNCOLS,iNEQ,btransposed)
          call storage_new ('lsyssc_resizeMatrixIndirect', 'h_Kld', isize,&
              ST_INT, rmatrix%h_Kld, ST_NEWBLOCK_ZERO)
        end if
        if (rmatrixTemplate%h_Kcol .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrixTemplate%h_Kcol, isize)
          if (present(NAMAX)) isize = max(isize,NAMAX)
          call storage_new ('lsyssc_resizeMatrixIndirect', 'h_Kcol', isize,&
              ST_INT, rmatrix%h_Kcol, ST_NEWBLOCK_ZERO)
        end if
        if (rmatrixTemplate%h_Kdiagonal .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrixTemplate%h_Kdiagonal, isize)
          iNEQ = isize
          if (present(NEQMAX)) iNEQ=max(iNEQ, NEQMAX+1)
          iNCOLS = isize
          if (present(NCOLSMAX)) iNCOLS=max(iNCOLS, NCOLSMAX+1)
          isize = merge(iNCOLS,iNEQ,btransposed)
          call storage_new ('lsyssc_resizeMatrixIndirect', 'h_Kdiagonal', isize,&
              ST_INT, rmatrix%h_Kdiagonal, ST_NEWBLOCK_ZERO)
        end if
        if (rmatrixTemplate%h_Da .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrixTemplate%h_Da, isize)
          if (present(NAMAX)) isize = max(isize,NAMAX)
          call storage_new ('lsyssc_resizeMatrixIndirect', 'h_Da', isize,&
              rmatrixTemplate%cdataType, rmatrix%h_Da, ST_NEWBLOCK_ZERO)
        end if

      else   ! Matrix should not be cleared

        if (rmatrixTemplate%h_Kld .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrixTemplate%h_Kld, isize)
          iNEQ = isize
          if (present(NEQMAX)) iNEQ=max(iNEQ, NEQMAX+1)
          iNCOLS = isize
          if (present(NCOLSMAX)) iNCOLS=max(iNCOLS, NCOLSMAX+1)
          isize = merge(iNCOLS,iNEQ,btransposed)
          call storage_new ('lsyssc_resizeMatrixIndirect', 'h_Kld', isize,&
              ST_INT, rmatrix%h_Kld, ST_NEWBLOCK_NOINIT)
        end if
        if (rmatrixTemplate%h_Kcol .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrixTemplate%h_Kcol, isize)
          if (present(NAMAX)) isize = max(isize,NAMAX)
          call storage_new ('lsyssc_resizeMatrixIndirect', 'h_Kcol', isize,&
              ST_INT, rmatrix%h_Kcol, ST_NEWBLOCK_NOINIT)
        end if
        if (rmatrixTemplate%h_Kdiagonal .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrixTemplate%h_Kdiagonal, isize)
          iNEQ = isize
          if (present(NEQMAX)) iNEQ=max(iNEQ, NEQMAX+1)
          iNCOLS = isize
          if (present(NCOLSMAX)) iNCOLS=max(iNCOLS, NCOLSMAX+1)
          isize = merge(iNCOLS,iNEQ,btransposed)
          call storage_new ('lsyssc_resizeMatrixIndirect', 'h_Kdiagonal', isize,&
              ST_INT, rmatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT)
        end if
        if (rmatrixTemplate%h_Da .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrixTemplate%h_Da, isize)
          if (present(NAMAX)) isize = max(isize,NAMAX)
          call storage_new ('lsyssc_resizeMatrixIndirect', 'h_Da', isize,&
              rmatrixTemplate%cdataType, rmatrix%h_Da, ST_NEWBLOCK_NOINIT)
        end if
      end if

    else

      ! The matrix has been initialised before.

      ! Check if matrices are compatible except for the dimensions
      if ((rmatrix%cmatrixFormat           .ne. rmatrixTemplate%cmatrixFormat) .or.&
          (rmatrix%cinterleavematrixFormat .ne. rmatrixTemplate%cinterleavematrixFormat) .or.&
          (rmatrix%NVAR                    .ne. rmatrixTemplate%NVAR) .or.&
          (rmatrix%cdataType               .ne. rmatrixTemplate%cdataType) .or.&
          iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne.&
          iand(rmatrixTemplate%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED)) then
        call output_line('Matrices are incompatible!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixIndirect')
        call sys_halt()
      end if

      ! Set working dimensions
      iNA    = rmatrixTemplate%NA
      if (present(NAMAX)) iNA = max(iNA,NAMAX)

      iNEQ   = rmatrixTemplate%NEQ
      iNCOLS = rmatrixTemplate%NCOLS

      if (btransposed) then
        if (present(NEQMAX)) iNCOLS = max(iNCOLS,NEQMAX)
        if (present(NCOLSMAX)) iNEQ = max(iNEQ,NCOLSMAX)
      else
        if (present(NEQMAX)) iNEQ = max(iNEQ,NEQMAX)
        if (present(NCOLSMAX)) iNCOLS = max(iNCOLS,NCOLSMAX)
      end if

      ! Set copy/clear attributes
      bdocopy = (.not.bclear)
      if (present(bcopy)) bdocopy = (bdocopy .and. bcopy)

      ! If the matrix should be cleared, then the sorting strategy (if any)
      ! can be ignored and reset. Otherwise, the matrix needs to be unsorted
      ! prior to copying some part of it. Afterwards, no sorting strategy is
      ! available in any case.
      if (bdocopy .and. (rmatrix%bcolumnsSorted .or. rmatrix%browsSorted)) then
        call lsyssc_sortMatrix(rmatrix,.false.)
      end if

      ! Reset sorting strategy, there is none
      rmatrix%bcolumnsSorted = .false.
      rmatrix%browsSorted = .false.
      nullify(rmatrix%p_rsortStrategyColumns)
      nullify(rmatrix%p_rsortStrategyRows)

      ! Update NA
      rmatrix%NA = rmatrixTemplate%NA

      ! Are we copy or not?
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .ne. 0) then

        ! Check if handle coincides with matrix dimensions
        call storage_getsize(rmatrix%h_Da,isize)
        if (rmatrix%NA*rmatrix%NVAR > isize) then
          call output_line('Dimensions of copied matrix mismatch!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixIndirect')
          call sys_halt()
        end if

      else

        ! Do we have to reallocate the handle?
        if (rmatrix%h_Da .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrix%h_Da, isize)
          if (rmatrix%NA*rmatrix%NVAR > isize) then

            ! Yes, we have to reallocate the handle.
            isize = iNA*rmatrix%NVAR

            ! Also consider the size of the template matrix.
            if (rmatrixTemplate%h_Da .ne. ST_NOHANDLE) then
              call storage_getsize(rmatrixTemplate%h_Da, isizeTmp)
              isize = max(0,isize,isizeTmp)
            end if

            ! Reallocate the memory
            call storage_realloc('lsyssc_resizeMatrixIndirect',&
                isize, rmatrix%h_Da, ST_NEWBLOCK_NOINIT, bdocopy)

          elseif (present(NAMAX)) then

            ! The available memory suffices for the matrix, i.e. isize <= NVAR*NA.
            ! Let us check if the user supplied a new upper limit which makes it
            ! mandatory to "shrink" the allocated memory. Note that memory for
            ! at least NA matrix intries is allocated.
            if (isize > iNA*rmatrix%NVAR) then

              ! Compute new size, i.e. MAX(0,NA,NAMAX)
              isize = iNA*rmatrix%NVAR

              if (isize .eq. 0) then
                ! If nothing is left, then the matrix can also be released.
                call lsyssc_releaseMatrix(rmatrix)
                return
              else
                call storage_realloc('lsyssc_resizeMatrixIndirect',&
                    isize, rmatrix%h_Da, ST_NEWBLOCK_NOINIT, bdocopy)
              end if
            end if
          end if
        end if
      end if

      ! Update NEQ and NCOLS.
      rmatrix%NEQ   = rmatrixTemplate%NEQ
      rmatrix%NCOLS = rmatrixTemplate%NCOLS

      ! Are we copy or not?
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0) then

        ! Check if handle coincides with matrix dimensions
        if (rmatrix%h_Kld .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrix%h_Kld,isize)
          if (btransposed .and. rmatrix%NCOLS+1 > isize .or. rmatrix%NEQ+1 > isize) then
            call output_line('Dimensions of copied matrix mismatch!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixIndirect')
            call sys_halt()
          end if
        end if

        ! Check if handle coincides with matrix dimensions
        if (rmatrix%h_Kcol .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrix%h_Kcol, isize)
          if (rmatrix%NA > isize) then
            call output_line('Dimensions of copied matrix mismatch!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixIndirect')
            call sys_halt()
          end if
        end if

        ! Check if handle coincides with matrix dimensions
        if (rmatrix%h_Kdiagonal .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrix%h_Kdiagonal, isize)
          if (btransposed .and. rmatrix%NCOLS > isize .or. rmatrix%NEQ > isize) then
            call output_line('Dimensions of copied matrix mismatch!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_resizeMatrixIndirect')
            call sys_halt()
          end if
        end if

      else  ! Matrix is no copy

        ! Do we have to reallocate the handle?
        if (rmatrix%h_Kld .ne. ST_NOHANDLE) then

          ! Do we process a virtually transposed matrix?
          if (btransposed) then
            call storage_getsize(rmatrix%h_Kld, isize)
            if (rmatrix%NCOLS+1 > isize) then

              ! Yes, we have to reallocate the handle.
              isize = iNCOLS+1

              ! Also consider the size of the template matrix.
              if (rmatrixTemplate%h_Kld .ne. ST_NOHANDLE) then
                call storage_getsize(rmatrixTemplate%h_Kld, isizeTmp)
                isize = max(0,isize,isizeTmp)
              end if

              ! Reallocate the memory
              call storage_realloc('lsyssc_resizeMatrixIndirect',&
                  isize, rmatrix%h_Kld, ST_NEWBLOCK_NOINIT, bdocopy)

            elseif (present(NCOLSMAX)) then

              ! The available memory suffices for the matrix, i.e. isize <= NCOLS.
              ! Let us check if the user supplied a new upper limit which makes it
              ! mandatory to "shrink" the allocated memory. Note that memory for
              ! at least NCOLS matrix intries is allocated.
              if (isize > iNCOLS+1) then

                ! Compute new size, i.e. MAX(0,NCOLS,NCOLSMAX)
                isize = iNCOLS

                if (isize .eq. 0) then
                  ! If nothing is left, then the matrix can also be released.
                  call lsyssc_releaseMatrix(rmatrix)
                  return
                else
                  call storage_realloc('lsyssc_resizeMatrixIndirect',&
                      isize+1, rmatrix%h_Kld, ST_NEWBLOCK_NOINIT, bdocopy)
                end if
              end if
            end if
          else
            call storage_getsize(rmatrix%h_Kld, isize)
            if (rmatrix%NEQ+1 > isize) then

              ! Yes, we have to reallocate the handle.
              isize = iNEQ+1

              ! Also consider the size of the template matrix.
              if (rmatrixTemplate%h_Kld .ne. ST_NOHANDLE) then
                call storage_getsize(rmatrixTemplate%h_Kld, isizeTmp)
                isize = max(0,isize,isizeTmp)
              end if

              ! Reallocate the memory
              call storage_realloc('lsyssc_resizeMatrixIndirect',&
                  isize, rmatrix%h_Kld, ST_NEWBLOCK_NOINIT, bdocopy)

            elseif (present(NEQMAX)) then

              ! The available memory suffices for the matrix, i.e. isize <= NEQ.
              ! Let us check if the user supplied a new upper limit which makes it
              ! mandatory to "shrink" the allocated memory. Note that memory for
              ! at least NEQ matrix intries is allocated.
              if (isize > iNEQ+1) then

                ! Compute new size, i.e. MAX(0,NEQ,NEQMAX)
                isize = iNEQ

                if (isize .eq. 0) then
                  ! If nothing is left, then the matrix can also be released.
                  call lsyssc_releaseMatrix(rmatrix)
                  return
                else
                  call storage_realloc('lsyssc_resizeMatrixIndirect',&
                      isize+1, rmatrix%h_Kld, ST_NEWBLOCK_NOINIT, bdocopy)
                end if
              end if
            end if
          end if
        end if

        ! Do we have to reallocate the handle?
        if (rmatrix%h_Kcol .ne. ST_NOHANDLE) then
          call storage_getsize(rmatrix%h_Kcol, isize)
          if (rmatrix%NA > isize) then

            ! Yes, we have to reallocate the handle.
            isize = iNA

            ! Also consider the size of the template matrix.
            if (rmatrixTemplate%h_Kcol .ne. ST_NOHANDLE) then
              call storage_getsize(rmatrixTemplate%h_Kcol, isizeTmp)
              isize = max(0,isize,isizeTmp)
            end if

            ! Reallocate the memory
            call storage_realloc('lsyssc_resizeMatrixIndirect',&
                isize, rmatrix%h_Kcol, ST_NEWBLOCK_NOINIT, bdocopy)

          elseif (present(NAMAX)) then

            ! The available memory suffices for the matrix, i.e. isize <= NA.
            ! Let us check if the user supplied a new upper limit which makes it
            ! mandatory to "shrink" the allocated memory. Note that memory for
            ! at least NA matrix intries is allocated.
            if (isize > iNA) then

              ! Compute new size, i.e. MAX(0,NA,NAMAX)
              isize = iNA

              if (isize .eq. 0) then
                ! If nothing is left, then the matrix can also be released.
                call lsyssc_releaseMatrix(rmatrix)
                return
              else
                call storage_realloc('lsyssc_resizeMatrixIndirect',&
                    isize, rmatrix%h_Kcol, ST_NEWBLOCK_NOINIT, bdocopy)
              end if
            end if
          end if
        end if

        ! Do we have to reallocate the handle?
        if (rmatrix%h_Kdiagonal .ne. ST_NOHANDLE) then

          ! Do we process a virtually transposed matrix?
          if (btransposed) then
            call storage_getsize(rmatrix%h_Kdiagonal, isize)
            if (rmatrix%NCOLS > isize) then

              ! Yes, we have to reallocate the handle.
              isize = iNCOLS

              ! Also consider the size of the template matrix.
              if (rmatrixTemplate%h_Kdiagonal .ne. ST_NOHANDLE) then
                call storage_getsize(rmatrixTemplate%h_Kdiagonal, isizeTmp)
                isize = max(0,isize,isizeTmp)
              end if

              ! Reallocate the memory
              call storage_realloc('lsyssc_resizeMatrixIndirect',&
                  isize, rmatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT, bdocopy)

            elseif (present(NCOLSMAX)) then

              ! The available memory suffices for the matrix, i.e. isize <= NCOLS.
              ! Let us check if the user supplied a new upper limit which makes it
              ! mandatory to "shrink" the allocated memory. Note that memory for
              ! at least NCOLS matrix intries is allocated.
              if (isize > iNCOLS) then

                ! Compute new size, i.e. MAX(0,NCOLS,NCOLSMAX)
                isize = iNCOLS

                if (isize .eq. 0) then
                  ! If nothing is left, then the matrix can also be released.
                  call lsyssc_releaseMatrix(rmatrix)
                  return
                else
                  call storage_realloc('lsyssc_resizeMatrixIndirect',&
                      isize, rmatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT, bdocopy)
                end if
              end if
            end if
          else
            call storage_getsize(rmatrix%h_Kdiagonal, isize)
            if (rmatrix%NEQ > isize) then

              ! Yes, we have to reallocate the handle.
              isize = iNEQ

              ! Also consider the size of the template matrix.
              if (rmatrixTemplate%h_Kdiagonal .ne. ST_NOHANDLE) then
                call storage_getsize(rmatrixTemplate%h_Kdiagonal, isizeTmp)
                isize = max(0,isize,isizeTmp)
              end if

              ! Reallocate the memory
              call storage_realloc('lsyssc_resizeMatrixIndirect',&
                  isize, rmatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT, bdocopy)

            elseif (present(NEQMAX)) then

              ! The available memory suffices for the matrix, i.e. isize <= NEQ.
              ! Let us check if the user supplied a new upper limit which makes it
              ! mandatory to "shrink" the allocated memory. Note that memory for
              ! at least NEQ matrix intries is allocated.
              if (isize > iNEQ) then

                ! Compute new size, i.e. MAX(0,NEQ,NEQMAX)
                isize = iNEQ

                if (isize .eq. 0) then
                  ! If nothing is left, then the matrix can also be released.
                  call lsyssc_releaseMatrix(rmatrix)
                  return
                else
                  call storage_realloc('lsyssc_resizeMatrixIndirect',&
                      isize, rmatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT, bdocopy)
                end if
              end if
            end if
          end if

        end if
      end if
    end if

  end subroutine

  !****************************************************************************

!<function>

  real(DP) function lsyssc_scalarProduct (rx, ry)

!<description>
  ! Calculates a scalar product of two vectors.
!</description>

!<input>
  ! First vector
  type(t_vectorScalar), intent(in) :: rx

  ! Second vector
  type(t_vectorScalar), intent(in) :: ry

!</input>

!<result>
  ! The scalar product (rx,ry) of the two vectors.
!</result>

!</function>

  ! local variables
  real(DP), dimension(:), pointer :: p_Ddata1dp
  real(DP), dimension(:), pointer :: p_Ddata2dp
  real(SP), dimension(:), pointer :: p_Fdata1dp
  real(SP), dimension(:), pointer :: p_Fdata2dp
  real(DP) :: res

  ! Vectors must be compatible!
  call lsyssc_isVectorCompatible (rx,ry)

  ! Is there data at all?
  res = 0.0_DP

  if ((rx%NEQ*rx%NVAR .eq. 0) .or. &
      (ry%NEQ*ry%NVAR .eq. 0) .or. &
      (rx%NEQ .ne. rx%NEQ) .or. &
      (rx%NVAR .ne. ry%NVAR)) then
    call output_line('Vector dimensions wrong!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarProduct')
    call sys_halt()
  end if

  if (rx%cdataType .ne. ry%cdataType) then
    call output_line('Data types different!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarProduct')
    call sys_halt()
  end if

  ! Take care of the data type before doing a scalar product!
  select case (rx%cdataType)
  case (ST_DOUBLE)

    ! Get the data arrays
    call lsyssc_getbase_double (rx,p_Ddata1dp)
    call lsyssc_getbase_double (ry,p_Ddata2dp)

    ! Perform the scalar product
    res = lalg_scalarProduct(p_Ddata1dp,p_Ddata2dp)

  case (ST_SINGLE)

    ! Get the data arrays
    call lsyssc_getbase_single (rx,p_Fdata1dp)
    call lsyssc_getbase_single (ry,p_Fdata2dp)

    ! Perform the scalar product
    res=lalg_scalarProduct(p_Fdata1dp,p_Fdata2dp)

  case default
    call output_line('Not supported precision combination!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarProduct')
    call sys_halt()
  end select

  ! Return the scalar product, finish
  lsyssc_scalarProduct = res

  end function

  !****************************************************************************
!<function>

  real(DP) function lsyssc_scalarProductMatVec (rx, ry)

!<description>
  ! Calculates a scalar product of two vectors, whereby the first
  ! vector is a diagonal matrix.
!</description>

!<input>
  ! First vector given as diagonal matrix
  type(t_MatrixScalar), intent(in) :: rx

  ! Second vector
  type(t_vectorScalar), intent(in) :: ry

!</input>

!<result>
  ! The scalar product (rx,ry) of the two vectors.
!</result>

!</function>

  ! local variables
  real(DP), dimension(:), pointer :: p_Ddata1dp
  real(DP), dimension(:), pointer :: p_Ddata2dp
  real(SP), dimension(:), pointer :: p_Fdata1dp
  real(SP), dimension(:), pointer :: p_Fdata2dp
  real(DP) :: res

  ! Matrix must be diagonal matrix
  if (rx%cmatrixFormat .ne. LSYSSC_MATRIXD) then
    call output_line('Matrix must be diagonal matrix!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarProductMatVec')
    call sys_halt()
  end if

  ! Vectors must be compatible!
  call lsyssc_isMatrixVectorCompatible (ry,rx,.false.)

  ! Is there data at all?
  res = 0.0_DP

  if ((rx%NEQ*rx%NVAR .eq. 0) .or. &
      (ry%NEQ*ry%NVAR .eq. 0) .or. &
      (rx%NEQ .ne. rx%NEQ) .or. &
      (rx%NVAR .ne. ry%NVAR)) then
    call output_line('Vector dimensions wrong!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarProductMatVec')
    call sys_halt()
  end if

  if (rx%cdataType .ne. ry%cdataType) then
    call output_line('Data types different!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarProductMatVec')
    call sys_halt()
  end if

  ! Take care of the data type before doing a scalar product!
  select case (rx%cdataType)
  case (ST_DOUBLE)

    ! Get the data arrays
    call lsyssc_getbase_double (rx,p_Ddata1dp)
    call lsyssc_getbase_double (ry,p_Ddata2dp)

    ! Perform the scalar product
    res=lalg_scalarProduct(p_Ddata1dp,p_Ddata2dp)

  case (ST_SINGLE)

    ! Get the data arrays
    call lsyssc_getbase_single (rx,p_Fdata1dp)
    call lsyssc_getbase_single (ry,p_Fdata2dp)

    ! Perform the scalar product
    res=lalg_scalarProduct(p_Fdata1dp,p_Fdata2dp)

  case default
    call output_line('Not supported precision combination!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarProductMatVec')
    call sys_halt()
  end select

  ! Return the scalar product, finish
  lsyssc_scalarProductMatVec = res

  end function

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_scalarMatVec (rmatrix, rx, ry, cx, cy, btranspose, rperfconfig)

!<description>
  ! Performs a matrix vector multiplication with a given scalar matrix:
  !   <tex> $$  Dy   =   cx * rMatrix * rx   +   cy * ry  $$ </tex>
!</description>

    use lsyssc_scalarMatVecAux

!<input>

  ! Scalar matrix
  type(t_matrixScalar), intent(in) :: rmatrix

  ! Vector to multiply with the matrix.
  type(t_vectorScalar), intent(in) :: rx

  ! Multiplicative factor for rx
  real(DP), intent(in) :: cx

  ! Multiplicative factor for ry
  real(DP), intent(in) :: cy

  ! OPTIONAL: Specifies whether a multiplication with the matrix itself
  ! (.false.) or with its transpose (.true.) should be performed. If not
  ! given, .false. is assumed.
  logical, optional, intent(in) :: btranspose

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! Additive vector. Receives the result of the matrix-vector multiplication
  type(t_vectorScalar), intent(inout) :: ry
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Da, p_Dx, p_Dy
    real(SP), dimension(:), pointer :: p_Fa, p_Fx, p_Fy
    integer, dimension(:), pointer :: p_Kld
    integer, dimension(:), pointer :: p_Kcol
    integer, dimension(:), pointer :: p_KrowIdx
    integer :: NEQ
    logical :: bvirt_trans
    logical :: btrans

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => lsyssc_perfconfig
    end if

    bvirt_trans = .false.
    btrans = .false.

    ! Should we multiply by the matrix transpose?
    if(present(btranspose)) btrans = btranspose

    ! If the scale factor is = 0 or cx = 0, we do not want to perform a MV
    ! multiplication. We just need to take care of the cy then...
    if (cx*rmatrix%dscaleFactor .eq. 0.0_DP) then

      if(cy .eq. 0.0_DP) then
        ! Simply clear ry
        call lsyssc_clearVector(ry)

      else if(cy .ne. 1.0_DP) then
        ! Scale ry by cy
        call lsyssc_scaleVector(ry,cy)

      end if

      ! And exit this routine
      return

    end if

    ! Vectors must be compatible to the matrix.
    call lsyssc_isMatrixCompatible (rx,rmatrix, btrans)
    call lsyssc_isMatrixCompatible (ry,rmatrix, .not. btrans)

    ! rx and ry must have at least the same data type!
    if (rx%cdataType .ne. ry%cdataType) then
      call output_line('Different data types for rx and ry not supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
      call sys_halt()
    end if

    ! Up to now, all matrix types use h_Da. So if that is not associated,
    ! there is for sure an error!
    if (rmatrix%h_Da .eq. ST_NOHANDLE) then
      call output_line('Matrix has no data!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
      call sys_halt()
    end if

    ! Is the matrix 'virtually transposed' ?
    bvirt_trans = (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0)

    ! By what do we multiply now? The matrix itself or its transposed?
    if (bvirt_trans .eqv. btrans) then

      ! We are multiplying by the matrix itself.
      ! Now is the matrix virtually transposed?
      if(bvirt_trans) then
        ! Yes, so exchange NEQ with NCOLS
        NEQ = rmatrix%NCOLS
      else
        ! No, so simply copy NEQ
        NEQ = rmatrix%NEQ
      end if

      ! Select the right MV multiplication routine from the matrix format
      select case (rmatrix%cmatrixFormat)

      case (LSYSSC_MATRIX7,LSYSSC_MATRIX9)
        ! Set pointers
        call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
        call lsyssc_getbase_Kld (rmatrix,p_Kld)

        ! Take care of the precision of the matrix
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Format 7 and Format 9 multiplication
          call lsyssc_getbase_double (rmatrix,p_Da)

          ! Take care of the precision of the vectors
          select case (rx%cdataType)
          case (ST_DOUBLE)
            ! Set pointers
            call lsyssc_getbase_double (rx,p_Dx)
            call lsyssc_getbase_double (ry,p_Dy)

            ! double precision matrix, double precision vectors
            call lsyssc_LAX79DbleDble (p_Kld,p_Kcol,p_Da,p_Dx,p_Dy,&
                cx*rmatrix%dscaleFactor,cy,NEQ,rx%NVAR,p_rperfconfig)

          case (ST_SINGLE)
            ! Set pointers
            call lsyssc_getbase_single (rx,p_Fx)
            call lsyssc_getbase_single (ry,p_Fy)

            ! double precision matrix, single precision vectors
            call lsyssc_LAX79DbleSngl (p_Kld,p_Kcol,p_Da,p_Fx,p_Fy,&
                cx*rmatrix%dscaleFactor,real(cy,SP),NEQ,rx%NVAR,p_rperfconfig)

          case default
            call output_line('Invalid combination of matrix/vector precision!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
            call sys_halt()
          end select

        case (ST_SINGLE)
          ! Format 7 and Format 9 multiplication
          call lsyssc_getbase_single (rmatrix,p_Fa)

          ! Take care of the precision of the vectors
          select case (rx%cdataType)

          case (ST_DOUBLE)
            ! Set pointers
            call lsyssc_getbase_double (rx,p_Dx)
            call lsyssc_getbase_double (ry,p_Dy)

            ! single precision matrix, double precision vectors
            call lsyssc_LAX79SnglDble (p_Kld,p_Kcol,p_Fa,p_Dx,p_Dy,&
                real(cx*rmatrix%dscaleFactor,SP),cy,NEQ,rx%NVAR,p_rperfconfig)

          case (ST_SINGLE)
            ! Set pointers
            call lsyssc_getbase_single (rx,p_Fx)
            call lsyssc_getbase_single (ry,p_Fy)

            ! single precision matrix, single precision vectors
            call lsyssc_LAX79SnglSngl (p_Kld,p_Kcol,p_Fa,p_Fx,p_Fy,&
                real(cx*rmatrix%dscaleFactor,SP),real(cy,SP),NEQ,rx%NVAR,p_rperfconfig)

          case default
            call output_line('Invalid combination of matrix/vector precision!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
            call sys_halt()
          end select

        case default
          call output_line('Invalid matrix precision!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
          call sys_halt()
        end select

      case (LSYSSC_MATRIX9ROWC)
        ! Set pointers
        call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
        call lsyssc_getbase_Kld (rmatrix,p_Kld)
        call lsyssc_getbase_KrowIdx (rmatrix,p_KrowIdx)

        ! Take care of the precision of the matrix
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Format 9 multiplication
          call lsyssc_getbase_double (rmatrix,p_Da)

          ! Take care of the precision of the vectors
          select case (rx%cdataType)
          case (ST_DOUBLE)
            ! Set pointers
            call lsyssc_getbase_double (rx,p_Dx)
            call lsyssc_getbase_double (ry,p_Dy)

            ! double precision matrix, double precision vectors
            call lsyssc_LAX9rowcDbleDble (p_Kld,p_Kcol,p_KrowIdx,p_Da,p_Dx,p_Dy,&
                cx*rmatrix%dscaleFactor,cy,rmatrix%nnzrows,rx%NVAR,p_rperfconfig)

          case (ST_SINGLE)
            ! Set pointers
            call lsyssc_getbase_single (rx,p_Fx)
            call lsyssc_getbase_single (ry,p_Fy)

            ! double precision matrix, double precision vectors
            call lsyssc_LAX9rowcDbleSngl (p_Kld,p_Kcol,p_KrowIdx,p_Da,p_Fx,p_Fy,&
                cx*rmatrix%dscaleFactor,real(cy,SP),rmatrix%nnzrows,rx%NVAR,p_rperfconfig)

          case default
            call output_line('Invalid combination of matrix/vector precision!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
            call sys_halt()
          end select

        case (ST_SINGLE)
          ! Format 9 multiplication
          call lsyssc_getbase_single (rmatrix,p_Fa)

          ! Take care of the precision of the vectors
          select case (rx%cdataType)
          case (ST_DOUBLE)
            ! Set pointers
            call lsyssc_getbase_double (rx,p_Dx)
            call lsyssc_getbase_double (ry,p_Dy)

            ! double precision matrix, double precision vectors
            call lsyssc_LAX9rowcSnglDble (p_Kld,p_Kcol,p_KrowIdx,p_Fa,p_Dx,p_Dy,&
                real(cx*rmatrix%dscaleFactor,SP),cy,rmatrix%nnzrows,rx%NVAR,p_rperfconfig)

          case (ST_SINGLE)
            ! Set pointers
            call lsyssc_getbase_single (rx,p_Fx)
            call lsyssc_getbase_single (ry,p_Fy)

            ! double precision matrix, double precision vectors
            call lsyssc_LAX9rowcSnglSngl (p_Kld,p_Kcol,p_KrowIdx,p_Fa,p_Fx,p_Fy,&
                real(cx*rmatrix%dscaleFactor,SP),real(cy,SP),rmatrix%nnzrows,rx%NVAR,&
                p_rperfconfig)

          case default
            call output_line('Invalid combination of matrix/vector precision!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
            call sys_halt()
          end select
        end select

      case (LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9INTL)
        ! Set pointers
        call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
        call lsyssc_getbase_Kld (rmatrix,p_Kld)

        ! Take care of the precision of the matrix
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Format 7 and Format 9 multiplication
          call lsyssc_getbase_double (rmatrix,p_Da)

          ! Take care of the precision of the vectors
          select case (rx%cdataType)
          case (ST_DOUBLE)
            ! Set pointers
            call lsyssc_getbase_double (rx,p_Dx)
            call lsyssc_getbase_double (ry,p_Dy)

            ! double precision matrix, double precision vectors
            select case (rmatrix%cinterleavematrixFormat)
            case (LSYSSC_MATRIX1)
              call lsyssc_LAX79INTL1DbleDble (p_Kld,p_Kcol,p_Da,p_Dx,p_Dy,&
                  cx*rmatrix%dscaleFactor,cy,NEQ,rx%NVAR,p_rperfconfig)

            case (LSYSSC_MATRIXD)
              call lsyssc_LAX79INTLDDbleDble (p_Kld,p_Kcol,p_Da,p_Dx,p_Dy,&
                  cx*rmatrix%dscaleFactor,cy,NEQ,rx%NVAR,p_rperfconfig)

            case default
              call output_line('Invalid interleave matrix format!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
              call sys_halt()
            end select

          case (ST_SINGLE)
            ! Set pointers
            call lsyssc_getbase_single (rx,p_Fx)
            call lsyssc_getbase_single (ry,p_Fy)

            ! double precision matrix, single precision vectors
            select case (rmatrix%cinterleavematrixFormat)
            case (LSYSSC_MATRIX1)
              call lsyssc_LAX79INTL1DbleSngl (p_Kld,p_Kcol,p_Da,p_Fx,p_Fy,&
                  cx*rmatrix%dscaleFactor,real(cy,SP),NEQ,rx%NVAR,p_rperfconfig)

            case (LSYSSC_MATRIXD)
              call lsyssc_LAX79INTLDDbleSngl (p_Kld,p_Kcol,p_Da,p_Fx,p_Fy,&
                  cx*rmatrix%dscaleFactor,real(cy,SP),NEQ,rx%NVAR,p_rperfconfig)

            case default
              call output_line('Invalid interleave matrix format!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
              call sys_halt()
            end select

          case default
            call output_line('Invalid combination of matrix/vector precision!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
            call sys_halt()
          end select

        case (ST_SINGLE)
          ! Format 7 and Format 9 multiplication
          call lsyssc_getbase_single (rmatrix,p_Fa)

          ! Take care of the precision of the vectors
          select case (rx%cdataType)
          case (ST_DOUBLE)
            ! Set pointers
            call lsyssc_getbase_double (rx,p_Dx)
            call lsyssc_getbase_double (ry,p_Dy)

            ! single precision matrix, double precision vectors
            select case (rmatrix%cinterleavematrixFormat)
            case (LSYSSC_MATRIX1)
              call lsyssc_LAX79INTL1SnglDble (p_Kld,p_Kcol,p_Fa,p_Dx,p_Dy,&
                  real(cx*rmatrix%dscaleFactor,SP),cy,NEQ,rx%NVAR,p_rperfconfig)

            case (LSYSSC_MATRIXD)
              call lsyssc_LAX79INTLDSnglDble (p_Kld,p_Kcol,p_Fa,p_Dx,p_Dy,&
                  real(cx*rmatrix%dscaleFactor,SP),cy,NEQ,rx%NVAR,p_rperfconfig)

            case default
              call output_line('Invalid interleave matrix format!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
              call sys_halt()
            end select

          case (ST_SINGLE)
            ! Set pointers
            call lsyssc_getbase_single (rx,p_Fx)
            call lsyssc_getbase_single (ry,p_Fy)

            ! single precision matrix, single precision vectors
            select case (rmatrix%cinterleavematrixFormat)
            case (LSYSSC_MATRIX1)
              call lsyssc_LAX79INTL1SnglSngl (p_Kld,p_Kcol,p_Fa,p_Fx,p_Fy,&
                  real(cx*rmatrix%dscaleFactor,SP),real(cy,SP),NEQ,rx%NVAR,&
                  p_rperfconfig)

            case (LSYSSC_MATRIXD)
              call lsyssc_LAX79INTLDSnglSngl (p_Kld,p_Kcol,p_Fa,p_Fx,p_Fy,&
                  real(cx*rmatrix%dscaleFactor,SP),real(cy,SP),NEQ,rx%NVAR,&
                  p_rperfconfig)

            case default
              call output_line('Invalid interleave matrix format!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
              call sys_halt()
            end select

          case default
            call output_line('Invalid combination of matrix/vector precision!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
            call sys_halt()
          end select

        case default
          call output_line('Invalid matrix precision!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
          call sys_halt()
        end select

      case (LSYSSC_MATRIXD)

        ! Take care of the precision of the matrix
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Format D multiplication
          call lsyssc_getbase_double (rmatrix,p_Da)

          ! Take care of the precision of the vectors
          select case (rx%cdataType)
          case (ST_DOUBLE)
            ! Set pointer
            call lsyssc_getbase_double (rx,p_Dx)
            call lsyssc_getbase_double (ry,p_Dy)

            ! double precision matrix, double precision vectors
            call lsyssc_LATXDDbleDble (p_Da,p_Dx,p_Dy,&
                cx*rmatrix%dscaleFactor,cy,NEQ,rx%NVAR,p_rperfconfig)

          case (ST_SINGLE)
            ! Set pointer
            call lsyssc_getbase_single (rx,p_Fx)
            call lsyssc_getbase_single (ry,p_Fy)

            ! double precision matrix, single precision vectors
            call lsyssc_LATXDDbleSngl (p_Da,p_Fx,p_Fy,&
                cx*rmatrix%dscaleFactor,real(cy,SP),NEQ,rx%NVAR,p_rperfconfig)

          case default
            call output_line('Invalid combination of matrix/vector precision!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
            call sys_halt()
          end select

        case (ST_SINGLE)
          ! Format D multiplication
          call lsyssc_getbase_single (rmatrix,p_Fa)

          ! Take care of the precision of the vectors
          select case (rx%cdataType)
          case (ST_DOUBLE)
            ! Set pointer
            call lsyssc_getbase_double (rx,p_Dx)
            call lsyssc_getbase_double (ry,p_Dy)

            ! single precision matrix, double precision vectors
            call lsyssc_LATXDSnglDble (p_Fa,p_Dx,p_Dy,&
                real(cx*rmatrix%dscaleFactor,SP),cy,NEQ,rx%NVAR,p_rperfconfig)

          case (ST_SINGLE)
            ! Set pointer
            call lsyssc_getbase_single (rx,p_Fx)
            call lsyssc_getbase_single (ry,p_Fy)

            ! single precision matrix, single precision vectors
            call lsyssc_LATXDSnglSngl (p_Fa,p_Fx,p_Fy,&
                real(cx*rmatrix%dscaleFactor,SP),real(cy,SP),NEQ,rx%NVAR,&
                p_rperfconfig)

          case default
            call output_line('Invalid combination of matrix/vector precision!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
            call sys_halt()
          end select

        case default
          call output_line('Invalid matrix precision!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
          call sys_halt()
        end select

      case default
        call output_line('Invalid matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
        call sys_halt()
      end select

    else
      ! We are multiplying by the transposed matrix.
      ! Now is the matrix virtually transposed?
      if(bvirt_trans) then
        ! Yes, so exchange NEQ with NCOLS
        NEQ = rmatrix%NCOLS
      else
        ! No, so simply copy NEQ
        NEQ = rmatrix%NEQ
      end if

      ! Select the right MV multiplication routine from the matrix format
      select case (rmatrix%cmatrixFormat)

      case (LSYSSC_MATRIX7,LSYSSC_MATRIX9)
        ! Set pointers
        call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
        call lsyssc_getbase_Kld (rmatrix,p_Kld)

        ! Take care of the precision of the matrix
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Format 7 and Format 9 multiplication
          call lsyssc_getbase_double (rmatrix,p_Da)

          ! Take care of the precision of the vectors
          select case (rx%cdataType)
          case (ST_DOUBLE)
            ! Set pointers
            call lsyssc_getbase_double (rx,p_Dx)
            call lsyssc_getbase_double (ry,p_Dy)

            ! double precision matrix, double precision vectors
            call lsyssc_LTX79DbleDble (p_Kld,p_Kcol,p_Da,p_Dx,p_Dy,&
                cx*rmatrix%dscaleFactor,cy,NEQ,rx%NVAR,p_rperfconfig)

          case (ST_SINGLE)
            ! Set pointers
            call lsyssc_getbase_single (rx,p_Fx)
            call lsyssc_getbase_single (ry,p_Fy)

            ! double precision matrix, single precision vectors
            call lsyssc_LTX79DbleSngl (p_Kld,p_Kcol,p_Da,p_Fx,p_Fy,&
                cx*rmatrix%dscaleFactor,real(cy,SP),NEQ,rx%NVAR,p_rperfconfig)

          case default
            call output_line('Invalid combination of matrix/vector precision!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
            call sys_halt()
          end select

        case (ST_SINGLE)
          ! Format 7 and Format 9 multiplication
          call lsyssc_getbase_single (rmatrix,p_Fa)

          ! Take care of the precision of the vectors
          select case (rx%cdataType)
          case (ST_DOUBLE)
            ! Set pointers
            call lsyssc_getbase_double (rx,p_Dx)
            call lsyssc_getbase_double (ry,p_Dy)

            ! single precision matrix, double precision vectors
            call lsyssc_LTX79SnglDble (p_Kld,p_Kcol,p_Fa,p_Dx,p_Dy,&
                real(cx*rmatrix%dscaleFactor,SP),cy,NEQ,rx%NVAR,p_rperfconfig)

          case (ST_SINGLE)
            ! Set pointers
            call lsyssc_getbase_single (rx,p_Fx)
            call lsyssc_getbase_single (ry,p_Fy)

            ! single precision matrix, single precision vectors
            call lsyssc_LTX79SnglSngl (p_Kld,p_Kcol,p_Fa,p_Fx,p_Fy,&
                real(cx*rmatrix%dscaleFactor,SP),real(cy,SP),NEQ,rx%NVAR,&
                p_rperfconfig)

          case default
            call output_line('Invalid combination of matrix/vector precision!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
            call sys_halt()
          end select

        case default
          call output_line('Invalid matrix precision!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
          call sys_halt()
        end select

      case (LSYSSC_MATRIXD)

        ! Take care of the precision of the matrix
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          ! Format D multiplication
          call lsyssc_getbase_double (rmatrix,p_Da)

          ! Take care of the precision of the vectors
          select case (rx%cdataType)
          case (ST_DOUBLE)
            ! Set pointer
            call lsyssc_getbase_double (rx,p_Dx)
            call lsyssc_getbase_double (ry,p_Dy)

            ! double precision matrix, double precision vectors
            call lsyssc_LATXDDbleDble (p_Da,p_Dx,p_Dy,&
                cx*rmatrix%dscaleFactor,cy,NEQ,rx%NVAR,p_rperfconfig)

          case (ST_SINGLE)
            ! Set pointer
            call lsyssc_getbase_single (rx,p_Fx)
            call lsyssc_getbase_single (ry,p_Fy)

            ! double precision matrix, single precision vectors
            call lsyssc_LATXDDbleSngl (p_Da,p_Fx,p_Fy,&
                cx*rmatrix%dscaleFactor,real(cy,SP),NEQ,rx%NVAR,p_rperfconfig)

          case default
            call output_line('Invalid combination of matrix/vector precision!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
            call sys_halt()
          end select

        case (ST_SINGLE)
          ! Format D multiplication
          call lsyssc_getbase_single (rmatrix,p_Fa)

          ! Take care of the precision of the vectors
          select case (rx%cdataType)
          case (ST_DOUBLE)
            ! Set pointer
            call lsyssc_getbase_double (rx,p_Dx)
            call lsyssc_getbase_double (ry,p_Dy)

            ! single precision matrix, double precision vectors
            call lsyssc_LATXDSnglDble (p_Fa,p_Dx,p_Dy,&
                real(cx*rmatrix%dscaleFactor,SP),cy,NEQ,rx%NVAR,p_rperfconfig)

          case (ST_SINGLE)
            ! Set pointer
            call lsyssc_getbase_single (rx,p_Fx)
            call lsyssc_getbase_single (ry,p_Fy)

            ! single precision matrix, single precision vectors
            call lsyssc_LATXDSnglSngl (p_Fa,p_Fx,p_Fy,&
                real(cx*rmatrix%dscaleFactor,SP),real(cy,SP),NEQ,rx%NVAR,&
                p_rperfconfig)

          case default
            call output_line('Invalid combination of matrix/vector precision!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
            call sys_halt()
          end select

        case default
          call output_line('Invalid matrix precision!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
          call sys_halt()
        end select

      case default
        call output_line('Invalid matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scalarMatVec')
        call sys_halt()
      end select

    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_duplicateVector (rx,ry,cdupStructure,cdupContent)

!<description>
  ! Duplicates an existing vector: ry := rx.
  ! Creates a new vector ry based on a template vector rx.
  ! Duplicating a vector does not necessarily mean that new memory is
  ! allocated and the vector entries are copied to that. The two flags
  ! cdupStructure and cdupContent decide on how to set up ry.
  ! Depending on their setting, it is possible to copy only the handle
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
  type(t_vectorScalar), intent(in) :: rx

  ! Duplication flag that decides on how to set up the structure
  ! of ry. Not all flags are possible!
  ! One of the LSYSSC_DUP_xxxx flags:
  ! LSYSSC_DUP_IGNORE     : Ignore the structure of rx
  ! LSYSSC_DUP_COPY or
  ! LSYSSC_DUP_COPYOVERWRITE or
  ! LSYSSC_DUP_TEMPLATE   : Structural data is copied from rx
  !   to ry (NEQ, sorting strategy, pointer to discretisation structure).
  integer, intent(in) :: cdupStructure

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
  integer, intent(in) :: cdupContent
!</input>

!<inputoutput>
  ! The new vector which will be a copy of roldvector
  type(t_vectorScalar), intent(inout) :: ry
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ioffset,cdataType
    logical :: bisCopy
    integer :: h_Ddata

    ! First the structure

    if (cdupStructure .ne. LSYSSC_DUP_IGNORE) then
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
    end if

    ! Now the content.

    select case (cdupContent)
    case (LSYSSC_DUP_IGNORE)
      ! Nothing to do

    case (LSYSSC_DUP_REMOVE)
      ! Release vector data
      if ((.not. ry%bisCopy) .and. (ry%h_Ddata .ne. ST_NOHANDLE)) &
        call storage_free (ry%h_Ddata)

      ry%h_Ddata = ST_NOHANDLE
      ry%bisCopy = .false.

    case (LSYSSC_DUP_DISMISS)
      ! Dismiss data
      ry%h_Ddata = ST_NOHANDLE
      ry%bisCopy = .false.

    case (LSYSSC_DUP_SHARE)
      ! Share information. Release memory if necessary.
      if ((.not. ry%bisCopy) .and. (ry%h_Ddata .ne. ST_NOHANDLE)) &
        call storage_free (ry%h_Ddata)

      ry%h_Ddata = rx%h_Ddata
      ry%iidxFirstEntry = rx%iidxFirstEntry
      ry%bisCopy = .true.

    case (LSYSSC_DUP_COPY)
      ! Copy information. If necessary, allocate new data -- by setting
      ! h_Ddata to ST_NOHANDLE before calling storage_copy.
      if (ry%bisCopy) then
        ry%h_Ddata = ST_NOHANDLE
        ry%bisCopy = .false.
      end if

      if (ry%h_Ddata .eq. ST_NOHANDLE) &
        call allocDestinationVector (rx,ry)

      call copy_data (rx,ry)

    case (LSYSSC_DUP_COPYOVERWRITE)
      ! Copy information, regardless of whether ry is the owner or not.
      ! If no memory is allocated, allocate new memory.

      if (ry%h_Ddata .eq. ST_NOHANDLE) &
        call allocDestinationVector (rx,ry)

      call copy_data (rx,ry)

    case (LSYSSC_DUP_ASIS)

      ! Copy by ownership. This is either LSYSSC_COPY or LSYSSC_SHARE,
      ! depending on whether rx is the owner of the data or not.
      if (rx%bisCopy) then
        ! rx shares it is data and thus ry will also.
        ry%h_Ddata = rx%h_Ddata
        ry%iidxFirstEntry = rx%iidxFirstEntry
        ry%bisCopy = .true.
      else
        ! The data belongs to rx and thus it must also belong to ry --
        ! so copy it.
        ry%h_Ddata = ST_NOHANDLE
        ry%bisCopy = .false.
        call allocDestinationVector (rx,ry)
        call copy_data (rx,ry)
      end if

    case (LSYSSC_DUP_EMPTY)

      ! Allocate new memory if ry is empty. Do not initialise.
      ! If ry contains data, we do not have to do anything.
      if (ry%h_Ddata .eq. ST_NOHANDLE) then
        call allocDestinationVector (rx,ry)
      end if

    case default

      call output_line('cdupContent unknown!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateVector')
      call sys_halt()

    end select

  contains

    subroutine allocDestinationVector (rx,ry)

    ! Allocates new memory ry in the size of rx. The memory is not initialised.

    type(t_vectorScalar), intent(in) :: rx
    type(t_vectorScalar), intent(inout) :: ry

      ! local variables
      integer :: NEQ,NVAR

      NEQ = rx%NEQ
      NVAR= rx%NVAR
      ry%cdataType = rx%cdataType
      call storage_new ('lsyssc_duplicateVector','ry',NEQ*NVAR,&
                        rx%cdataType, ry%h_Ddata, &
                        ST_NEWBLOCK_NOINIT)

    end subroutine

    ! ---------------------------------------------------------------

    subroutine copy_data (rx,ry)

    ! Copies the content of rx to ry. Takes care of the data type

    type(t_vectorScalar), intent(in) :: rx
    type(t_vectorScalar), intent(inout) :: ry

      ! local variables
      real(DP), dimension(:), pointer :: p_Dsource,p_Ddest
      real(SP), dimension(:), pointer :: p_Fsource,p_Fdest

      if ((rx%h_Ddata .eq. ry%h_Ddata) .and.&
          (rx%iidxFirstEntry .eq. ry%iidxFirstEntry)) then
        ! Eehm... forget it.
        return
      end if

      ! Take care of the data type!
      select case (rx%cdataType)
      case (ST_DOUBLE)
        ! Get the pointer and scale the whole data array.
        call lsyssc_getbase_double(rx,p_Dsource)
        call lsyssc_getbase_double(ry,p_Ddest)
        call lalg_copyVectorDble (p_Dsource,p_Ddest)

      case (ST_SINGLE)
        ! Get the pointer and scale the whole data array.
        call lsyssc_getbase_single(rx,p_Fsource)
        call lsyssc_getbase_single(ry,p_Fdest)
        call lalg_copyVectorSngl (p_Fsource,p_Fdest)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateVector')
        call sys_halt()
      end select

    end subroutine

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_duplicateMatrix (rsourceMatrix,rdestMatrix,&
                                     cdupStructure, cdupContent)

!<description>
  ! Duplicates an existing matrix, creates a new matrix rdestMatrix based
  ! on a template matrix rsourceMatrix.
  ! Duplicating a matrix does not necessarily mean that new memory is
  ! allocated and the matrix entries are copied to that. The two flags
  ! cdupStructure and cdupContent decide on how to set up rdestMatrix.
  ! Depending on their setting, it is possible to copy only the handles
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
  type(t_matrixScalar), intent(in) :: rsourceMatrix

  ! Duplication flag that decides on how to set up the structure
  ! of rdestMatrix. One of the LSYSSC_DUP_xxxx flags:
  ! LSYSSC_DUP_IGNORE     : Do not set up the structure of rdestMatrix. Any
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
  integer, intent(in) :: cdupStructure

  ! Duplication flag that decides on how to set up the content
  ! of rdestMatrix. One of the LSYSSC_DUP_xxxx flags:
  ! LSYSSC_DUP_IGNORE     : Do not set up the content of rdestMatrix. Any
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
  integer, intent(in) :: cdupContent
!</input>

!<output>
  ! Destination matrix.
  type(t_matrixScalar), intent(inout) :: rdestMatrix
!</output>

!</subroutine>

    ! local variables
    integer :: isize

    ! What do we have to do for the structure?
    select case (cdupStructure)
    case (LSYSSC_DUP_IGNORE)
      ! Nothing
    case (LSYSSC_DUP_REMOVE)
      ! Remove the structure - if there is any.
      call removeStructure(rdestMatrix, .true.)

    case (LSYSSC_DUP_DISMISS)
      ! Dismiss the structure - if there is any.
      call removeStructure(rdestMatrix, .false.)

    case (LSYSSC_DUP_TEMPLATE)
      ! Remove the structure - if there is any.
      call removeStructure(rdestMatrix, .true.)

      ! Copy static structural information
      call copyStaticStructure(rsourceMatrix,rdestMatrix)

    case (LSYSSC_DUP_SHARE)
      ! Remove the structure - if there is any.
      call removeStructure(rdestMatrix, .true.)

      ! Share the structure between rsourceMatrix and rdestMatrix
      call shareStructure(rsourceMatrix, rdestMatrix)

    case (LSYSSC_DUP_COPY)
      ! Copy the structure of rsourceMatrix to rdestMatrix
      call copyStructure(rsourceMatrix, rdestMatrix, .false.)

    case (LSYSSC_DUP_COPYOVERWRITE)
      ! Copy the structure of rsourceMatrix to rdestMatrix, do not respect ownership
      call copyStructure(rsourceMatrix, rdestMatrix, .true.)

    case (LSYSSC_DUP_ASIS)
      ! What is with the source matrix. Does the structure belong to it?
      if (iand(rsourceMatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0) then
        ! Copy the structure
        call copyStructure(rsourceMatrix, rdestMatrix, .false.)
      else
        ! Remove the structure - if there is any.
        call removeStructure(rdestMatrix, .true.)

        ! Share the structure
        call shareStructure(rsourceMatrix, rdestMatrix)
      end if

    case (LSYSSC_DUP_EMPTY)
      ! Start with sharing the structure
      call shareStructure(rsourceMatrix, rdestMatrix)
      ! Reset ownership
      rdestMatrix%imatrixSpec = iand(rdestMatrix%imatrixSpec,&
                                     not(LSYSSC_MSPEC_STRUCTUREISCOPY))

      ! And at last recreate the arrays.
      ! Which source matrix do we have?
      select case (rsourceMatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
        call storage_new ('lsyssc_duplicateMatrix', 'KCOL', &
            rsourceMatrix%NA, ST_INT, &
            rdestMatrix%h_Kcol, ST_NEWBLOCK_NOINIT)

        call storage_new ('lsyssc_duplicateMatrix', 'KLD', &
            rdestMatrix%NEQ+1, ST_INT, &
            rdestMatrix%h_Kld, ST_NEWBLOCK_NOINIT)

        call storage_new ('lsyssc_duplicateMatrix', 'Kdiagonal', &
            rdestMatrix%NEQ, ST_INT, &
            rdestMatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT)

      case (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
        call storage_new ('lsyssc_duplicateMatrix', 'KCOL', &
            rdestMatrix%NA, ST_INT, &
            rdestMatrix%h_Kcol, ST_NEWBLOCK_NOINIT)

        call storage_new ('lsyssc_duplicateMatrix', 'KLD', &
            rdestMatrix%NEQ+1, ST_INT, &
            rdestMatrix%h_Kld, ST_NEWBLOCK_NOINIT)

      case (LSYSSC_MATRIXD)
        ! Nothing to do
      end select
    end select

    ! -----
    ! Ok, handling of the structure is finished. The data is handled similar.
    ! -----

    ! What do we have to do for the content?
    select case (cdupContent)
    case (LSYSSC_DUP_IGNORE)
      ! Nothing
    case (LSYSSC_DUP_REMOVE)
      ! Remove the content - if there is any.
      call removeContent(rdestMatrix, .true.)

    case (LSYSSC_DUP_DISMISS)
      ! Dismiss the structure - if there is any.
      call removeContent(rdestMatrix, .false.)

    case (LSYSSC_DUP_SHARE)
      ! Remove the content - if there is any.
      call removeContent(rdestMatrix, .true.)

      ! Share the structure between rsourceMatrix and rdestMatrix
      call shareContent(rsourceMatrix, rdestMatrix)

    case (LSYSSC_DUP_COPY)
      ! Copy the structure of rsourceMatrix to rdestMatrix
      call copyContent(rsourceMatrix, rdestMatrix, .false.)

    case (LSYSSC_DUP_COPYOVERWRITE)
      ! Copy the structure of rsourceMatrix to rdestMatrix, do not respect ownership
      call copyContent(rsourceMatrix, rdestMatrix, .true.)

    case (LSYSSC_DUP_TEMPLATE)
      ! Remove the structure - if there is any.
      call removeContent(rdestMatrix, .true.)

      ! Copy static content information
      call copyStaticContent(rsourceMatrix,rdestMatrix)

    case (LSYSSC_DUP_ASIS)
      ! What is with the source matrix. Does the structure belong to it?
      if (iand(rsourceMatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .ne. 0) then
        ! Copy the structure
        call copyContent(rsourceMatrix, rdestMatrix, .false.)
      else
        ! Remove the content - if there is any.
        call removeContent(rdestMatrix, .true.)

        ! Share the structure
        call shareContent(rsourceMatrix, rdestMatrix)
      end if

    case (LSYSSC_DUP_EMPTY)
      ! Start with sharing the structure
      call shareContent(rsourceMatrix, rdestMatrix)
      ! Reset ownership
      rdestMatrix%imatrixSpec = iand(rdestMatrix%imatrixSpec,&
                                     not(LSYSSC_MSPEC_CONTENTISCOPY))

      ! And at last recreate the arrays.
      ! Which source matrix do we have?
      select case (rsourceMatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9,LSYSSC_MATRIX9ROWC,LSYSSC_MATRIX7,LSYSSC_MATRIXD)
        ! Create a new content array in the same data type as the original matrix
        call storage_new('lsyssc_duplicateMatrix', 'Da', rdestMatrix%NA, &
            rsourceMatrix%cdataType, rdestMatrix%h_Da, ST_NEWBLOCK_NOINIT)
      case (LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7INTL)
        ! Create a new content array in the same data type as the original matrix
        select case (rsourceMatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIX1)
          call storage_new('lsyssc_duplicateMatrix', 'Da', &
              rdestMatrix%NA*rdestMatrix%NVAR*rdestMatrix%NVAR, &
              rsourceMatrix%cdataType, rdestMatrix%h_Da, ST_NEWBLOCK_NOINIT)
        case (LSYSSC_MATRIXD)
          call storage_new('lsyssc_duplicateMatrix', 'Da', &
              rdestMatrix%NA*rdestMatrix%NVAR,rsourceMatrix%cdataType, &
              rdestMatrix%h_Da, ST_NEWBLOCK_NOINIT)
        case default
          call output_line('Wrong matrix format of interleaved matrix!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
          call sys_halt()
        end select
      end select
    end select

    ! -----
    ! Final check. Check if we destroyed the matrix. May only happen if the
    ! user on-purpose calls this routine two times with different source
    ! but the same destination matrix.
    ! -----

    select case (rdestMatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)

      ! Check length of Da
      if (rdestMatrix%h_Da .ne. ST_NOHANDLE) then
        call storage_getsize (rdestMatrix%h_Da,isize)
        select case (rdestMatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIX1)
          if (isize .lt. rdestMatrix%NA * rdestMatrix%NVAR * rdestMatrix%NVAR) then
            call output_line('Matrix destroyed; NA < length(Da)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        case (LSYSSC_MATRIXD)
          if (isize .lt. rdestMatrix%NA * rdestMatrix%NVAR) then
            call output_line('Matrix destroyed; NA < length(Da)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        case default
          if (isize .lt. rdestMatrix%NA) then
            call output_line('Matrix destroyed; NA < length(Da)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        end select
      end if

      ! Check length of KCOL
      if (rdestMatrix%h_Kcol .ne. ST_NOHANDLE) then
        call storage_getsize (rdestMatrix%h_Kcol,isize)
        if (isize .lt. rdestMatrix%NA) then
          call output_line('Matrix destroyed; NA < length(KCOL)!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
          call sys_halt()
        end if
      end if

      ! Check length of KLD
      if (rdestMatrix%h_Kld .ne. ST_NOHANDLE) then
        call storage_getsize (rdestMatrix%h_Kld,isize)

        ! Be careful, matrix may be transposed.
        if (iand(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .eq. 0) then
          if (isize .lt. rdestMatrix%NEQ+1) then
            call output_line('Matrix destroyed; NEQ+1 < length(KLD)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        else
          if (isize .lt. rdestMatrix%NCOLS+1) then
            call output_line('Matrix destroyed; NEQ+1 < length(KLD)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        end if
      end if

      ! Check length of Kdiagonal
      if (rdestMatrix%h_Kdiagonal .ne. ST_NOHANDLE) then
        call storage_getsize (rdestMatrix%h_Kdiagonal,isize)

        ! Be careful, matrix may be transposed.
        if (iand(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .eq. 0) then
          if (isize .lt. rdestMatrix%NEQ) then
            call output_line('Matrix destroyed; NEQ+1 < length(Kdiag)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        else
          if (isize .lt. rdestMatrix%NCOLS) then
            call output_line('Matrix destroyed; NEQ+1 < length(Kdiag)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        end if
      end if

    case (LSYSSC_MATRIX9ROWC)

      ! Check length of Da
      if (rdestMatrix%h_Da .ne. ST_NOHANDLE) then
        call storage_getsize (rdestMatrix%h_Da,isize)
        select case (rdestMatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIX1)
          if (isize .lt. rdestMatrix%NA * rdestMatrix%NVAR * rdestMatrix%NVAR) then
            call output_line('Matrix destroyed; NA < length(Da)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        case (LSYSSC_MATRIXD)
          if (isize .lt. rdestMatrix%NA * rdestMatrix%NVAR) then
            call output_line('Matrix destroyed; NA < length(Da)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        case default
          if (isize .lt. rdestMatrix%NA) then
            call output_line('Matrix destroyed; NA < length(Da)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        end select
      end if

      ! Check length of KCOL
      if (rdestMatrix%h_Kcol .ne. ST_NOHANDLE) then
        call storage_getsize (rdestMatrix%h_Kcol,isize)
        if (isize .lt. rdestMatrix%NA) then
          call output_line('Matrix destroyed; NA < length(KCOL)!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
          call sys_halt()
        end if
      end if

      ! Check length of KLD
      if (rdestMatrix%h_Kld .ne. ST_NOHANDLE) then
        call storage_getsize (rdestMatrix%h_Kld,isize)

        ! Be careful, matrix may be transposed.
        if (iand(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .eq. 0) then
          if (isize .lt. rdestMatrix%nnzrows+1) then
            call output_line('Matrix destroyed; NEQ+1 < length(KLD)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        else
          if (isize .lt. rdestMatrix%NCOLS+1) then
            call output_line('Matrix destroyed; NEQ+1 < length(KLD)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        end if
      end if

    case (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)

      ! Check length of Da
      if (rdestMatrix%h_Da .ne. ST_NOHANDLE) then
        call storage_getsize (rdestMatrix%h_Da,isize)
        select case (rdestMatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIX1)
          if (isize .lt. rdestMatrix%NA * rdestMatrix%NVAR*rdestMatrix%NVAR) then
            call output_line('Matrix destroyed; NA < length(Da)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        case (LSYSSC_MATRIXD)
          if (isize .lt. rdestMatrix%NA * rdestMatrix%NVAR) then
            call output_line('Matrix destroyed; NA < length(Da)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        case default
          if (isize .lt. rdestMatrix%NA) then
            call output_line('Matrix destroyed; NA < length(Da)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        end select
      end if

      ! Check length of KCOL
      if (rdestMatrix%h_Kcol .ne. ST_NOHANDLE) then
        call storage_getsize (rdestMatrix%h_Kcol,isize)
        if (isize .lt. rdestMatrix%NA) then
          call output_line('Matrix destroyed; NA < length(KCOL)!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
          call sys_halt()
        end if
      end if

      ! Check length of KLD
      if (rdestMatrix%h_Kld .ne. ST_NOHANDLE) then
        call storage_getsize (rdestMatrix%h_Kld,isize)

        ! Be careful, matrix may be transposed.
        if (iand(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .eq. 0) then
          if (isize .lt. rdestMatrix%NEQ+1) then
            call output_line('Matrix destroyed; NEQ+1 < length(KLD)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        else
          if (isize .lt. rdestMatrix%NCOLS+1) then
            call output_line('Matrix destroyed; NEQ+1 < length(KLD)!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
            call sys_halt()
          end if
        end if
      end if

    case (LSYSSC_MATRIXD)

      ! Check length of Da
      if (rdestMatrix%h_Da .ne. ST_NOHANDLE) then
        call storage_getsize (rdestMatrix%h_Da,isize)
        if (isize .ne. rdestMatrix%NA) then
          call output_line('Matrix destroyed; NA != length(Da)!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
          call sys_halt()
        end if
      end if

    end select

  contains

    !--------------------------------------------------------
    ! Auxiliary routine: Remove the structure from a matrix.

    subroutine removeStructure(rmatrix, brelease)

    ! The matrix to to be processed.
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! Whether to release data from the heap (if it belongs to the matrix)
    ! or simply to overwrite the handles with ST_NOHANDLE
    logical, intent(in) :: brelease

      ! This is a matrix-dependent task
      select case (rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
        ! Release the handles from the heap?
        ! Only release it if the data belongs to this matrix.
        if (brelease .and. &
            (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .eq. 0)) then
          if (rmatrix%h_Kcol .ne. ST_NOHANDLE) call storage_free (rmatrix%h_Kcol)
          if (rmatrix%h_Kld .ne. ST_NOHANDLE) call storage_free (rmatrix%h_Kld)
          if (rmatrix%h_Kdiagonal .ne. ST_NOHANDLE) call storage_free (rmatrix%h_Kdiagonal)
        end if

        ! Reset the handles
        rmatrix%h_Kld = ST_NOHANDLE
        rmatrix%h_Kcol = ST_NOHANDLE
        rmatrix%h_Kdiagonal = ST_NOHANDLE

        ! Reset the ownership-status
        rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,&
                                   not(LSYSSC_MSPEC_STRUCTUREISCOPY))
      case (LSYSSC_MATRIX9ROWC)
        ! Release the handles from the heap?
        ! Only release it if the data belongs to this matrix.
        if (brelease .and. &
            (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .eq. 0)) then
          if (rmatrix%h_Kcol .ne. ST_NOHANDLE) call storage_free (rmatrix%h_Kcol)
          if (rmatrix%h_Kld .ne. ST_NOHANDLE) call storage_free (rmatrix%h_Kld)
          if (rmatrix%h_Kdiagonal .ne. ST_NOHANDLE) call storage_free (rmatrix%h_Kdiagonal)
          if (rmatrix%h_KrowIdx .ne. ST_NOHANDLE) call storage_free (rmatrix%h_KrowIdx)
        end if

        ! Reset the handles
        rmatrix%h_Kld = ST_NOHANDLE
        rmatrix%h_Kcol = ST_NOHANDLE
        rmatrix%h_Kdiagonal = ST_NOHANDLE
        rmatrix%h_KrowIdx = ST_NOHANDLE

        ! Reset the ownership-status
        rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,&
                                   not(LSYSSC_MSPEC_STRUCTUREISCOPY))
      case (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
        ! Release the handles from the heap?
        ! Only release it if the data belongs to this matrix.
        if (brelease .and. &
            (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .eq. 0)) then
          if (rmatrix%h_Kcol .ne. ST_NOHANDLE) call storage_free (rmatrix%h_Kcol)
          if (rmatrix%h_Kld .ne. ST_NOHANDLE) call storage_free (rmatrix%h_Kld)
        end if

        ! Reset the handles
        rmatrix%h_Kld = ST_NOHANDLE
        rmatrix%h_Kcol = ST_NOHANDLE

      case (LSYSSC_MATRIXD,LSYSSC_MATRIX1)
        ! Nothing to do

      end select

      ! Reset the ownership-status
      rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,&
                                  not(LSYSSC_MSPEC_STRUCTUREISCOPY))

    end subroutine

    !--------------------------------------------------------
    ! Auxiliary routine: Shares the structure between rsourceMatrix
    ! and rdestMatrix

    subroutine shareStructure(rsourceMatrix, rdestMatrix)

    ! The source matrix
    type(t_matrixScalar), intent(in) :: rsourceMatrix

    ! The destination matrix
    type(t_matrixScalar), intent(inout) :: rdestMatrix

      ! local variables
      integer(I32) :: iflag,iflag2

      ! Overwrite structural data
      rdestMatrix%NA     = rsourceMatrix%NA
      rdestMatrix%NEQ    = rsourceMatrix%NEQ
      rdestMatrix%NCOLS  = rsourceMatrix%NCOLS
      rdestMatrix%NVAR   = rsourceMatrix%NVAR
      rdestMatrix%cmatrixFormat           = rsourceMatrix%cmatrixFormat
      rdestMatrix%cinterleavematrixFormat = rsourceMatrix%cinterleavematrixFormat
      rdestMatrix%bcolumnsSorted          = rsourceMatrix%bcolumnsSorted
      rdestMatrix%browsSorted             = rsourceMatrix%browsSorted
      rdestMatrix%p_rsortStrategyColumns => rsourceMatrix%p_rsortStrategyColumns
      rdestMatrix%p_rsortStrategyRows => rsourceMatrix%p_rsortStrategyRows

      ! Transfer all flags except the 'dup' flags
      iflag = iand(rsourceMatrix%imatrixSpec,not(LSYSSC_MSPEC_ISCOPY))
      iflag2 = iand(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
      rdestMatrix%imatrixSpec = ior(iflag,iflag2)

      ! Transfer discretisation-related information,
      rdestMatrix%p_rspatialDiscrTest => rsourceMatrix%p_rspatialDiscrTest
      rdestMatrix%p_rspatialDiscrTrial => rsourceMatrix%p_rspatialDiscrTrial
      rdestMatrix%bidenticalTrialAndTest = rsourceMatrix%bidenticalTrialAndTest

      ! Which source matrix do we have?
      select case (rsourceMatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
        rdestMatrix%h_Kcol = rsourceMatrix%h_Kcol
        rdestMatrix%h_Kld  = rsourceMatrix%h_Kld
        rdestMatrix%h_Kdiagonal = rsourceMatrix%h_Kdiagonal

        ! Indicate via the matrixSpec-flag that we are not
        ! the owner of the structure. Exception: If there is no
        ! structure, there is no ownership!
        if ((rsourceMatrix%h_Kcol .ne. 0) .and.&
            (rsourceMatrix%h_Kld .ne. 0) .and.&
            (rsourceMatrix%h_Kdiagonal .ne. 0)) then
          rdestMatrix%imatrixSpec = ior(rdestMatrix%imatrixSpec,&
                                        LSYSSC_MSPEC_STRUCTUREISCOPY)
        else
          rdestMatrix%imatrixSpec = iand(rdestMatrix%imatrixSpec,&
                                        not(LSYSSC_MSPEC_STRUCTUREISCOPY))
        end if

      case (LSYSSC_MATRIX9ROWC)
        rdestMatrix%h_Kcol = rsourceMatrix%h_Kcol
        rdestMatrix%h_Kld  = rsourceMatrix%h_Kld
        rdestMatrix%h_KrowIdx = rsourceMatrix%h_KrowIdx
        rdestMatrix%NNZROWS = rsourceMatrix%NNZROWS

        ! Indicate via the matrixSpec-flag that we are not
        ! the owner of the structure. Exception: If there is no
        ! structure, there is no ownership!
        if ((rsourceMatrix%h_Kcol .ne. 0) .and.&
            (rsourceMatrix%h_Kld .ne. 0) .and.&
            (rsourceMatrix%h_KrowIdx .ne. 0)) then
          rdestMatrix%imatrixSpec = ior(rdestMatrix%imatrixSpec,&
                                        LSYSSC_MSPEC_STRUCTUREISCOPY)
        else
          rdestMatrix%imatrixSpec = iand(rdestMatrix%imatrixSpec,&
                                        not(LSYSSC_MSPEC_STRUCTUREISCOPY))
        end if

      case (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
        ! Overwrite structural data
        rdestMatrix%h_Kcol = rsourceMatrix%h_Kcol
        rdestMatrix%h_Kld  = rsourceMatrix%h_Kld

        ! Indicate via the matrixSpec-flag that we are not
        ! the owner of the structure. Exception: If there is no
        ! structure, there is no ownership!
        if ((rsourceMatrix%h_Kcol .ne. 0) .and.&
            (rsourceMatrix%h_Kld .ne. 0)) then
          rdestMatrix%imatrixSpec = ior(rdestMatrix%imatrixSpec,&
                                        LSYSSC_MSPEC_STRUCTUREISCOPY)
        else
          rdestMatrix%imatrixSpec = iand(rdestMatrix%imatrixSpec,&
                                        not(LSYSSC_MSPEC_STRUCTUREISCOPY))
        end if

      case (LSYSSC_MATRIXD,LSYSSC_MATRIX1)
        ! Nothing to do. No ownership.
        rdestMatrix%imatrixSpec = iand(rdestMatrix%imatrixSpec,&
                                      not(LSYSSC_MSPEC_STRUCTUREISCOPY))

      end select

    end subroutine

    !--------------------------------------------------------
    ! Auxiliary routine: Copy static structural information
    ! (NEQ, NCOLS, NA,...) from rsourceMatrix to rdestMatrix

    subroutine copyStaticStructure(rsourceMatrix, rdestMatrix)

    ! The source matrix
    type(t_matrixScalar), intent(in) :: rsourceMatrix

    ! The destination matrix
    type(t_matrixScalar), intent(inout) :: rdestMatrix

      ! local variables
      integer(I32) :: iflag,iflag2

      ! Overwrite structural data
      rdestMatrix%NA     = rsourceMatrix%NA
      rdestMatrix%NEQ    = rsourceMatrix%NEQ
      rdestMatrix%NCOLS  = rsourceMatrix%NCOLS
      rdestMatrix%NVAR   = rsourceMatrix%NVAR
      rdestMatrix%NNZROWS = rsourceMatrix%NNZROWS
      rdestMatrix%cmatrixFormat           = rsourceMatrix%cmatrixFormat
      rdestMatrix%cinterleavematrixFormat = rsourceMatrix%cinterleavematrixFormat
      rdestMatrix%bcolumnsSorted          = rsourceMatrix%bcolumnsSorted
      rdestMatrix%browsSorted             = rsourceMatrix%browsSorted
      rdestMatrix%p_rsortStrategyColumns => rsourceMatrix%p_rsortStrategyColumns
      rdestMatrix%p_rsortStrategyRows    => rsourceMatrix%p_rsortStrategyRows

      ! Transfer all flags except the 'dup' flags
      iflag = iand(rsourceMatrix%imatrixSpec,not(LSYSSC_MSPEC_ISCOPY))
      iflag2 = iand(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
      rdestMatrix%imatrixSpec = ior(iflag,iflag2)

      ! Transfer discretisation-related information,
      rdestMatrix%p_rspatialDiscrTrial => rsourceMatrix%p_rspatialDiscrTrial
      rdestMatrix%p_rspatialDiscrTest => rsourceMatrix%p_rspatialDiscrTest
      rdestMatrix%bidenticalTrialAndTest = rsourceMatrix%bidenticalTrialAndTest

    end subroutine

    !--------------------------------------------------------
    ! Auxiliary routine: Copy the structure of rsourceMatrix
    ! to rdestMatrix

    subroutine copyStructure(rsourceMatrix, rdestMatrix, bignoreOwner)

    ! The source matrix
    type(t_matrixScalar), intent(in) :: rsourceMatrix

    ! The destination matrix
    type(t_matrixScalar), intent(inout) :: rdestMatrix

    ! Whether to respect the ownership of the arrays and allocate
    ! memory automatically if a matrix is not the owner of the arrays.
    logical, intent(in) :: bignoreOwner

      ! local variables
      integer :: isize
      integer :: NEQ,NNZROWS
      logical :: bremove

      ! Overwrite structural data
      call copyStaticStructure(rsourceMatrix, rdestMatrix)

      ! Check the structure if it exists and if it has the right
      ! size - then we can overwrite!
      bRemove = .false.

      ! But at first, check if rdestMatrix is the owner of the matrix
      ! structure:
      if ((.not. bignoreOwner) .and. &
          (iand(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0)) then
        ! No, the structure belongs to someone else, but we want to have
        ! our own. Detach the foreign matrix structure.
        bremove = .true.

      else

        ! Ok, rdestMatrix owns its own matrix structure or does not have any.
        ! Check the structure if it is large enough so we can overwrite
        ! it. If the size of the arrays is not equal to those
        ! in the source matrix, release the structure and allocate a new one.

        if (iand(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
          NEQ = rdestMatrix%NCOLS
        else
          NEQ = rdestMatrix%NEQ
        end if

        ! Which source matrix do we have?
        select case (rsourceMatrix%cmatrixFormat)
        case (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
          if ((rdestMatrix%h_Kcol .ne. ST_NOHANDLE) .and. &
              (rdestMatrix%h_Kld .ne. ST_NOHANDLE) .and. &
              (rdestMatrix%h_Kdiagonal .ne. ST_NOHANDLE)) then

            call storage_getsize (rdestMatrix%h_Kcol,isize)
            bremove = bremove .or. (isize .lt. rdestMatrix%NA)

            call storage_getsize (rsourceMatrix%h_Kld,isize)
            bremove = bremove .or. (isize .lt. NEQ+1)

            call storage_getsize (rsourceMatrix%h_Kdiagonal,isize)
            bremove = bremove .or. (isize .lt. NEQ)

          else

            ! Remove any partial information if there is any.
            bremove = .true.

          end if

        case (LSYSSC_MATRIX9ROWC)
          NNZROWS = rdestMatrix%NNZROWS

          if ((rdestMatrix%h_Kcol .ne. ST_NOHANDLE) .and. &
              (rdestMatrix%h_Kld .ne. ST_NOHANDLE) .and. &
              (rdestMatrix%h_KrowIdx .ne. ST_NOHANDLE) .and. &
              (rdestMatrix%h_Kdiagonal .ne. ST_NOHANDLE)) then

            call storage_getsize (rdestMatrix%h_Kcol,isize)
            bremove = bremove .or. (isize .lt. rdestMatrix%NA)

            call storage_getsize (rsourceMatrix%h_Kld,isize)
            bremove = bremove .or. (isize .lt. NEQ+1)

            call storage_getsize (rsourceMatrix%h_Kdiagonal,isize)
            bremove = bremove .or. (isize .lt. NEQ)

            call storage_getsize (rsourceMatrix%h_KrowIdx,isize)
            bremove = bremove .or. (isize .lt. NNZROWS)

          else

            ! Remove any partial information if there is any.
            bremove = .true.

          end if

        case (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)

          if ((rdestMatrix%h_Kcol .ne. ST_NOHANDLE) .and. &
              (rdestMatrix%h_Kld .ne. ST_NOHANDLE)) then

            call storage_getsize (rdestMatrix%h_Kcol,isize)
            bremove = bremove .or. (isize .lt. rdestMatrix%NA)

            call storage_getsize (rsourceMatrix%h_Kld,isize)
            bremove = bremove .or. (isize .lt. NEQ+1)

          else

            ! Remove any partial information if there is any.
            bremove = .true.

          end if

        case (LSYSSC_MATRIXD,LSYSSC_MATRIX1)
          ! Nothing to do

        end select

      end if

      ! Remove the old matrix structure if we should do so.
      if (bremove) call removeStructure(rdestMatrix, .true.)

      ! Duplicate structural data from the source matrix.
      ! Storage_copy allocates new memory if necessary, as the handles are all
      ! set to ST_NOHANDLE with the above removeStructure!
      ! If the handles exist, they have the correct size, so we can overwrite
      ! the entries.
      !
      ! Which source matrix do we have?
      select case (rsourceMatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
        call lsyssc_auxcopy_Kcol (rsourceMatrix,rdestMatrix)
        call lsyssc_auxcopy_Kld (rsourceMatrix,rdestMatrix)
        call lsyssc_auxcopy_Kdiagonal (rsourceMatrix,rdestMatrix)

      case (LSYSSC_MATRIX9ROWC)
        call lsyssc_auxcopy_Kcol (rsourceMatrix,rdestMatrix)
        call lsyssc_auxcopy_Kld (rsourceMatrix,rdestMatrix)
        call lsyssc_auxcopy_Kdiagonal (rsourceMatrix,rdestMatrix)
        call lsyssc_auxcopy_KrowIdx (rsourceMatrix,rdestMatrix)

      case (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
        call lsyssc_auxcopy_Kcol (rsourceMatrix,rdestMatrix)
        call lsyssc_auxcopy_Kld (rsourceMatrix,rdestMatrix)

      case (LSYSSC_MATRIXD,LSYSSC_MATRIX1)
        ! Nothing to do

      end select

      if (.not. bignoreOwner) then
        ! Indicate via the matrixSpec-flag that we are the owner of the structure.
        rdestMatrix%imatrixSpec = iand(rdestMatrix%imatrixSpec,&
                                      not(LSYSSC_MSPEC_STRUCTUREISCOPY))
      end if

    end subroutine

    !--------------------------------------------------------
    ! Auxiliary routine: Remove the content from a matrix.

    subroutine removeContent(rmatrix, brelease)

    ! The matrix to to be processed.
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! Whether to release data from the heap (if it belongs to the matrix)
    ! or simply to overwrite the handles with ST_NOHANDLE
    logical, intent(in) :: brelease

      ! This is a matrix-dependent task
      select case (rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9,LSYSSC_MATRIX9ROWC,LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7,&
          LSYSSC_MATRIX7INTL,LSYSSC_MATRIXD)
        ! Release the handles from the heap?
        ! Only release it if the data belongs to this matrix.
        if (brelease .and. &
            (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .eq. 0)) then
          if (rmatrix%h_Da .ne. ST_NOHANDLE) call storage_free (rmatrix%h_Da)
        end if

        ! Reset the handles
        rmatrix%h_Da = ST_NOHANDLE

      end select

      ! Reset the ownership-status
      rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,&
                                  not(LSYSSC_MSPEC_CONTENTISCOPY))

    end subroutine

    !--------------------------------------------------------
    ! Auxiliary routine: Shares the content between rsourceMatrix
    ! and rdestMatrix

    subroutine shareContent(rsourceMatrix, rdestMatrix)

    ! The source matrix
    type(t_matrixScalar), intent(in) :: rsourceMatrix

    ! The destination matrix
    type(t_matrixScalar), intent(inout) :: rdestMatrix

      ! Remove the old content - if there is any.
      call removeContent(rdestMatrix, .true.)

      ! Overwrite structural data
      rdestMatrix%dscaleFactor = rsourceMatrix%dscaleFactor
      rdestMatrix%cdataType    = rsourceMatrix%cdataType

      ! Which source matrix do we have?
      select case (rsourceMatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9,LSYSSC_MATRIX9ROWC,LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7,&
          LSYSSC_MATRIX7INTL,LSYSSC_MATRIXD)

        rdestMatrix%h_Da = rsourceMatrix%h_Da

        ! Indicate via the matrixSpec-flag that we are not
        ! the owner of the structure. Only exception: If the handle
        ! is zero, there is nothing that can be owned!
        if (rsourceMatrix%h_Da .ne. 0) then
          rdestMatrix%imatrixSpec = ior(rdestMatrix%imatrixSpec,&
                                        LSYSSC_MSPEC_CONTENTISCOPY)
        else
          rdestMatrix%imatrixSpec = iand(rdestMatrix%imatrixSpec,&
                                        not(LSYSSC_MSPEC_CONTENTISCOPY))
        end if

      end select

    end subroutine

    !--------------------------------------------------------
    ! Auxiliary routine: Copy static content-related
    ! information (scale factor,...) rsourceMatrix
    ! to rdestMatrix

    subroutine copyStaticContent(rsourceMatrix, rdestMatrix)

    ! The source matrix
    type(t_matrixScalar), intent(in) :: rsourceMatrix

    ! The destination matrix
    type(t_matrixScalar), intent(inout) :: rdestMatrix

      ! Overwrite structural data
      rdestMatrix%dscaleFactor = rsourceMatrix%dscaleFactor

    end subroutine

    !--------------------------------------------------------
    ! Auxiliary routine: Copy the content of rsourceMatrix
    ! to rdestMatrix

    subroutine copyContent(rsourceMatrix, rdestMatrix, bignoreOwner)

    ! The source matrix
    type(t_matrixScalar), intent(in) :: rsourceMatrix

    ! The destination matrix
    type(t_matrixScalar), intent(inout) :: rdestMatrix

    ! Whether to respect the ownership of the arrays and allocate
    ! memory automatically if a matrix is not the owner of the arrays.
    logical, intent(in) :: bignoreOwner

      ! local variables
      logical :: bremove

      ! Overwrite structural data
      call copyStaticContent(rsourceMatrix, rdestMatrix)

      ! Check the content if it exists and if it has the right
      ! size - then we can overwrite!
      bRemove = .false.

      ! But at first, check if rdestMatrix is the owner of the matrix
      ! structure:
      if ((.not. bignoreOwner) .and. &
          (iand(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .ne. 0)) then

        ! No, the content belongs to someone else, but we want to have
        ! our own. Detach the foreign matrix content.
        bremove = .true.

      else

        ! Ok, rdestMatrix owns some matrix content.
        ! Check the structure if it is large enough so we can overwrite
        ! it. If the size of the arrays is not equal to those
        ! in the source matrix, release the content and allocate a new one.

        ! Which source matrix do we have?
        select case (rsourceMatrix%cmatrixFormat)
        case (LSYSSC_MATRIX9,LSYSSC_MATRIX9ROWC,LSYSSC_MATRIX7,LSYSSC_MATRIXD)

          if (rdestMatrix%h_Da .ne. ST_NOHANDLE) then

            call storage_getsize (rdestMatrix%h_Da,isize)
            bremove = bremove .or. (isize .lt. rdestMatrix%NA)

            ! Check the data type
            bremove = bremove .or. (rdestMatrix%cdataType .ne. rsourceMatrix%cdataType)

          else

            ! Remove any partial information if there is any
            bremove = .true.

          end if

          case (LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7INTL)

          if (rdestMatrix%h_Da .ne. ST_NOHANDLE) then

            call storage_getsize (rdestMatrix%h_Da,isize)
            select case (rdestMatrix%cinterleavematrixFormat)
            case (LSYSSC_MATRIX1)
              bremove = bremove .or. &
                        (isize .lt. rdestMatrix%NA*rdestMatrix%NVAR*rdestMatrix%NVAR)
            case (LSYSSC_MATRIXD)
              bremove = bremove .or. (isize .lt. rdestMatrix%NA*rdestMatrix%NVAR)
            case default
              call output_line('Wrong interleave matrix format!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_duplicateMatrix')
              call sys_halt()
            end select

            ! Check the data type
            bremove = bremove .or. (rdestMatrix%cdataType .ne. rsourceMatrix%cdataType)

          else

            ! Remove any partial information if there is any
            bremove = .true.

          end if

        end select

      end if

      ! Remove the old content - if we should do so
      if (bremove) call removeContent(rdestMatrix, .true.)

      rdestMatrix%cdataType = rsourceMatrix%cdataType

      ! Duplicate content data from the source matrix.
      ! Storage_copy allocates new memory as the handles are all
      ! set to ST_NOHANDLE with the above removeStructure!
      !
      ! Which source matrix do we have?
      select case (rsourceMatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7,&
          LSYSSC_MATRIX7INTL,LSYSSC_MATRIXD,LSYSSC_MATRIX1)

        call lsyssc_auxcopy_da (rsourceMatrix,rdestMatrix)

      end select

      if (.not. bignoreOwner) then
        ! Indicate via the matrixSpec-flag that we are the owner of the structure.
        rdestMatrix%imatrixSpec = iand(rdestMatrix%imatrixSpec,&
                                      not(LSYSSC_MSPEC_CONTENTISCOPY))
      end if

    end subroutine

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_releaseVector (rvector, bquiet)

!<description>
  ! Releases a vector from memory. The vector structure is cleaned up.
  ! Remark: The memory associated to the vector is only released if this
  !  vector structure is the owner of the data.
!</description>

!<inputoutput>

  ! Vector to release.
  type(t_vectorScalar), intent(inout) :: rvector

!</inputoutput>

!<input>
  ! OPTIONAL: Specifies whether a warning should be printed when released an
  ! empty vector (bquiet = .false.) or whether to remain silent in this case.
  ! If not specified, bquiet = .false. is used.
  logical, optional, intent(in) :: bquiet
!</input>

!</subroutine>

  ! local variables
  logical :: bwarn

  ! Check whether the vector is empty or not
  if (rvector%h_Ddata .eq. ST_NOHANDLE) then

    ! Shout or shut up?
    bwarn = .true.
    if(present(bquiet)) bwarn = .not. bquiet

    if(bwarn) then
      call output_line('Releasing unused vector!',&
          OU_CLASS_WARNING,OU_MODE_STD,'lsyssc_releaseVector')
    end if

  end if

  ! Clean up the data structure.
  ! Do not release the vector data if the handle belongs to another
  ! vector!
  if ((.not. rvector%bisCopy) .and. (rvector%h_Ddata .ne. ST_NOHANDLE)) then
    call storage_free(rvector%h_Ddata)
  else
    rvector%h_Ddata = ST_NOHANDLE
  end if
  rvector%NEQ = 0
  rvector%NVAR = 1
  rvector%cdataType = ST_DOUBLE
  rvector%iidxFirstEntry = 1
  rvector%bisCopy = .false.
  rvector%bisSorted = .false.
  nullify(rvector%p_rsortStrategy)
  rvector%p_rspatialDiscr => null()

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_releaseMatrixContent (rmatrix)

!<description>
  ! Releases the content of a matrix from memory. The structure stays unchanged.
!</description>

!<inputoutput>

  ! Matrix to release.
  type(t_matrixScalar), intent(inout) :: rmatrix

!</inputoutput>

!</subroutine>

    ! Which handles do we have to release?
    !
    ! Release the matrix data if the handle is not a copy of another matrix
    if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .eq. 0) then
      select case (rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9,LSYSSC_MATRIX9ROWC,LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7&
          &,LSYSSC_MATRIX7INTL,LSYSSC_MATRIXD,LSYSSC_MATRIX1)
        ! Release matrix data, structure 9,7
        if (rmatrix%h_Da .ne. ST_NOHANDLE) then
          call storage_free(rmatrix%h_Da)
        end if
      end select
    end if

    rmatrix%h_Da        = ST_NOHANDLE

    ! Set the copy-flag appropriately
    rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_CONTENTISCOPY))

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_releaseMatrix (rmatrix)

!<description>
  ! Releases a matrix from memory. The structure is cleaned up.
  ! All data vectors belonging to the structure are released;
  ! data structures belonging to other matrices are not.
!</description>

!<inputoutput>

  ! Matrix to release.
  type(t_matrixScalar), intent(inout) :: rmatrix

!</inputoutput>

!</subroutine>

  ! Matrix template; initialised by the default initialisation strategy
  ! of Fortran 90.
  type(t_matrixScalar) :: rmatrixTemplate

  ! Release the matrix content.
  call lsyssc_releaseMatrixContent(rmatrix)

  ! Release the structure if it does not belong to another vector
  if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .eq. 0) then
    ! Release matrix structure
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
      if (rmatrix%h_Kcol .ne. ST_NOHANDLE)      call storage_free(rmatrix%h_Kcol)
      if (rmatrix%h_Kld .ne. ST_NOHANDLE)       call storage_free(rmatrix%h_Kld)
      if (rmatrix%h_Kdiagonal .ne. ST_NOHANDLE) call storage_free(rmatrix%h_Kdiagonal)

    case (LSYSSC_MATRIX9ROWC)
      if (rmatrix%h_Kcol .ne. ST_NOHANDLE)      call storage_free(rmatrix%h_Kcol)
      if (rmatrix%h_Kld .ne. ST_NOHANDLE)       call storage_free(rmatrix%h_Kld)
      if (rmatrix%h_Kdiagonal .ne. ST_NOHANDLE) call storage_free(rmatrix%h_Kdiagonal)
      if (rmatrix%h_KrowIdx .ne. ST_NOHANDLE)   call storage_free(rmatrix%h_KrowIdx)

    case (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
      if (rmatrix%h_Kcol .ne. ST_NOHANDLE) call storage_free(rmatrix%h_Kcol)
      if (rmatrix%h_Kld .ne. ST_NOHANDLE)  call storage_free(rmatrix%h_Kld)

    case (LSYSSC_MATRIXD)
      ! Nothing to do
    end select
  end if

  ! Clean up the rest. For this purpose, overwrite the matrix structure
  ! by a matrix structure which is initialised by the default Fortran
  ! initialisation routines.
  rmatrix = rmatrixTemplate
  !rmatrix%h_Da        = ST_NOHANDLE
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
  !rmatrix%bcolumnsSorted = .false.
  !rmatrix%browsSorted = .false.
  !nullify(rmatrix%p_rsortStrategyColumns)
  !nullify(rmatrix%p_rsortStrategyRows)
  !rmatrix%p_rspatialDiscretisation => NULL()

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_clearMatrix (rmatrix,dvalue)

!<description>
  ! Clears the entries in a matrix. All entries are overwritten with 0.0 or
  ! with dvalue (if specified).
!</description>

!<inputoutput>

  ! Matrix to clear.
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! OPTIONAL: Value to write into the matrix.
  ! If not specified, all matrix entries are set to 0.0.
  ! If specified, all matrix entries are set to dvalue.
  real(DP), intent(in), optional :: dvalue

!</inputoutput>

!</subroutine>

  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata

  if (rmatrix%NEQ .le. 0) return ! Empty matrix

  if (lsyssc_isExplicitMatrix1D(rmatrix)) then

    ! Get the data array and clear it -- depending on the data type.
    select case (rmatrix%cdataType)
    case (ST_DOUBLE)
      call lsyssc_getbase_double (rmatrix,p_Ddata)
      if (.not. present(dvalue)) then
        call lalg_clearVectorDble (p_Ddata)
      else
        call lalg_setVectorDble (p_Ddata,dvalue)
      end if
    case (ST_SINGLE)
      call lsyssc_getbase_single (rmatrix,p_Fdata)
      if (.not. present(dvalue)) then
        call lalg_clearVectorSngl (p_Fdata)
      else
        call lalg_setVectorSngl (p_Fdata,real(dvalue,SP))
      end if
    case default
      call output_line('Unsupported Data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_clearMatrix')
      call sys_halt()
    end select

  end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_initialiseIdentityMatrix (rmatrix)

!<description>
  ! Initialises the matrix rmatrix to an identity matrix.
  ! The matrix structure must already have been set up. If necessary, new data
  ! for the matrix is allocated on the heap.
  ! The scaling factor of the matrix remains unchanged!
!</description>

!<inputoutput>

  ! Matrix to release.
  type(t_matrixScalar), intent(inout) :: rmatrix

!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    real(SP), dimension(:), pointer :: p_Fdata
    integer, dimension(:), pointer :: p_Kdiagonal
    integer :: i

    if (rmatrix%NEQ .le. 0) return ! Empty matrix

    ! Which matrix type do we have?
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7,&
        LSYSSC_MATRIXD,LSYSSC_MATRIX1)
      ! If necessary, allocate memory.
      if (rmatrix%h_Da .eq. ST_NOHANDLE) then
        call lsyssc_allocEmptyMatrix (rmatrix,LSYSSC_SETM_UNDEFINED)
      end if
    end select

    ! Now choose -- depending on the matrix format -- how to initialise the
    ! matrix data.
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9)
      ! Clear the old content
      call lsyssc_clearMatrix (rmatrix)

      ! Get the structure and the data.
      ! Put the diagonal elements to 1.
      call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)

      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrix,p_Ddata)
        do i=1,rmatrix%NEQ
          p_Ddata(p_Kdiagonal(i)) = 1.0_DP
        end do
      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrix,p_Fdata)
        do i=1,rmatrix%NEQ
          p_Fdata(p_Kdiagonal(i)) = 1.0_SP
        end do
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_initialiseIdentityMatrix')
        call sys_halt()
      end select

    case (LSYSSC_MATRIX7)
      ! Clear the old content
      call lsyssc_clearMatrix (rmatrix)

      ! Get the structure and the data.
      ! Put the diagonal elements to 1.
      call lsyssc_getbase_Kld (rmatrix,p_Kdiagonal)

      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrix,p_Ddata)
        do i=1,rmatrix%NEQ
          p_Ddata(p_Kdiagonal(i)) = 1.0_DP
        end do
      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrix,p_Fdata)
        do i=1,rmatrix%NEQ
          p_Fdata(p_Kdiagonal(i)) = 1.0_SP
        end do
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_initialiseIdentityMatrix')
        call sys_halt()
      end select

    case (LSYSSC_MATRIX1)
      ! Clear the old content
      call lsyssc_clearMatrix (rmatrix)

      ! Get the structure and the data.
      ! Put the diagonal elements to 1.
      call lsyssc_getbase_Kld (rmatrix,p_Kdiagonal)

      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrix,p_Ddata)
        do i=1,rmatrix%NEQ
          p_Ddata(i*(rmatrix%NEQ-1)+i) = 1.0_DP
        end do
      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrix,p_Fdata)
        do i=1,rmatrix%NEQ
          p_Fdata(i*(rmatrix%NEQ-1)+i) = 1.0_SP
        end do
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_initialiseIdentityMatrix')
        call sys_halt()
      end select

    case (LSYSSC_MATRIXD)
      ! Put all elements to 1.0
      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrix,p_Ddata)
        p_Ddata(:) = 1.0_DP
      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrix,p_Fdata)
        p_Fdata(:) = 1.0_SP
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_initialiseIdentityMatrix')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_initialiseIdentityMatrix')
      call sys_halt()

    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_convertMatrix (rmatrix,cmatrixFormat,bconvertEntries)

!<description>
  ! Tries to convert a matrix rmatrix in-situ into a different matrix format.
  ! If necessary, memory is allocated for the new structure. (May be the case
  ! if the structure belongs to another matrix).
  ! If the matrix cannot be converted (due to format incompatibility),
  ! an error is thrown.
!</description>

!<input>
  ! Destination format of the matrix. One of the LSYSSC_MATRIXx constants.
  integer, intent(in) :: cmatrixFormat

  ! OPTIONAL: Whether to convert the entries of the matrix. Standard = TRUE.
  ! If set to FALSE, the entries are not converted, only the structure
  ! is converted (thus leaving the entries to an undefined state).
  logical, intent(in), optional :: bconvertEntries
!</input>

!<inputoutput>
  ! Matrix to convert.
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  integer :: ihandle
  real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2
  real(SP), dimension(:), pointer :: p_Fdata,p_Fdata2
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
  logical :: bentries,bstrucOwner,bentryOwner
  integer :: NA,iA,iEQ
  integer :: j,nrows,ncols
  integer :: h_Da,h_Kcol,h_Kld,h_Kdiagonal

  ! Matrix is already in that format.
  if (rmatrix%cmatrixFormat .eq. cmatrixFormat) return

  ! Empty matrix
  if (rmatrix%NEQ .le. 0) return

  bentries = .true.
  if (present(bconvertEntries)) bentries = bconvertEntries

  ! Be careful with the structure during conversion: If the structure belongs
  ! to another matrix, we must not deallocate it!
  bstrucOwner = iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .eq. 0

  ! The same for the entries!
  bentryOwner = iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .eq. 0

  ! Which matrix type do we have?
  select case (rmatrix%cmatrixFormat)
  case (LSYSSC_MATRIX9)

    select case (cmatrixFormat)
    case (LSYSSC_MATRIX7)

      if (.not. bstrucOwner) then

        ! Copy the structure so we are allowed to modify it.
        ihandle = ST_NOHANDLE
        call storage_copy(rmatrix%h_Kcol,ihandle)
        rmatrix%h_Kcol = ihandle

        ihandle = ST_NOHANDLE
        call storage_copy(rmatrix%h_Kld,ihandle)
        rmatrix%h_Kld = ihandle

        ! Kdiagonal is thrown away later, we do not need to copy it.

        rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_STRUCTUREISCOPY))

      end if

      if ((.not. bentryOwner) .and. bentries .and. (rmatrix%h_Da .ne. ST_NOHANDLE)) then

        ! Copy the entries so we are allowed to modify them.
        ihandle = ST_NOHANDLE
        call storage_copy(rmatrix%h_Da,ihandle)
        rmatrix%h_Da = ihandle

        rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_CONTENTISCOPY))

      end if

      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)
      call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)

      ! Check that the matrix can be converted. There is a format error
      ! if there is no diagonal element.
      do i=1,rmatrix%NEQ
        if (p_Kcol(p_Kdiagonal(i)) .ne. i) then
          call output_line('Incompatible Format-9 matrix!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
          call sys_halt()
        end if
      end do

      if ((.not. bentries) .or. (rmatrix%h_Da .eq. ST_NOHANDLE)) then

        ! No matrix entries, only resort the structure
        call lsyssc_unsortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)

        ! Release diagonal pointer if it belongs to us
        if (bstrucOwner) then
          call storage_free (rmatrix%h_Kdiagonal)
        else
          rmatrix%h_Kdiagonal = ST_NOHANDLE
        end if

        rmatrix%cmatrixFormat = LSYSSC_MATRIX7

      else

        ! Convert from structure 9 to structure 7. Use the sortCSRxxxx
        ! routine below.
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double (rmatrix,p_Ddata)
          call lsyssc_unsortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ, p_Ddata)

        case (ST_SINGLE)
          call lsyssc_getbase_single (rmatrix,p_Fdata)
          call lsyssc_unsortCSRsingle (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ, p_Fdata)

        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
          call sys_halt()
        end select

        ! Release diagonal pointer if it belongs to us
        if (bstrucOwner) then
          call storage_free (rmatrix%h_Kdiagonal)
        else
          rmatrix%h_Kdiagonal = ST_NOHANDLE
        end if

        rmatrix%cmatrixFormat = LSYSSC_MATRIX7
      end if

    case (LSYSSC_MATRIXD)

      if (bentries .and. (rmatrix%h_Da .ne. ST_NOHANDLE)) then

        ! Convert from structure 9 to structure D by copying the
        ! diagonal to the front and reallocation of the memory
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)

          call lsyssc_getbase_double (rmatrix,p_Ddata)

          if (.not. bentryOwner) then
            ! Allocate new memory for the entries
            rmatrix%h_Da = ST_NOHANDLE
            call storage_new ('lsyssc_convertMatrix', 'Da', &
                  rmatrix%NEQ, ST_DOUBLE, rmatrix%h_Da, ST_NEWBLOCK_NOINIT)
            call storage_getbase_double (rmatrix%h_Da,p_Ddata2)
            rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_CONTENTISCOPY))
          else
            ! Destinatinon pointer points to source matrix.
            ! This overwrites the original entries which is no problem as
            ! extracting the diagonal is always a 'compression' overwriting
            ! information that is not used anymore.
            p_Ddata2 => p_Ddata
          end if

          call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
          do i=1,rmatrix%NEQ
            p_Ddata2(i) = p_Ddata(p_Kdiagonal(i))
          end do

          rmatrix%cmatrixFormat = LSYSSC_MATRIXD

        case (ST_SINGLE)

          call lsyssc_getbase_single (rmatrix,p_Fdata)

          if (.not. bentryOwner) then
            ! Allocate new memory for the entries
            rmatrix%h_Da = ST_NOHANDLE
            call storage_new ('lsyssc_convertMatrix', 'Da', &
                  rmatrix%NEQ, ST_SINGLE, rmatrix%h_Da, ST_NEWBLOCK_NOINIT)
            call storage_getbase_single (rmatrix%h_Da,p_Fdata2)
            rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_CONTENTISCOPY))
          else
            ! Destinatinon pointer points to source matrix.
            ! This overwrites the original entries which is no problem as
            ! extracting the diagonal is always a 'compression' overwriting
            ! information that is not used anymore.
            p_Fdata2 => p_Fdata
          end if

          call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
          do i=1,rmatrix%NEQ
            p_Fdata2(i) = p_Fdata(p_Kdiagonal(i))
          end do

          rmatrix%cmatrixFormat = LSYSSC_MATRIXD

        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
          call sys_halt()
        end select

        ! Reallocate the entries-array to have only the diagonal entries.
        ! In case we are not the owner, we have allocated new memory and thus
        ! do not need to resize it again.
        if (bentryOwner) then
          call storage_realloc ('lsyssc_convertMatrix', rmatrix%NEQ, &
                                rmatrix%h_Da, ST_NEWBLOCK_NOINIT, bentries)
        end if

        rmatrix%NA = rmatrix%NEQ

      end if

      ! Release unused information
      if (bstrucOwner) then
        call storage_free (rmatrix%h_Kdiagonal)
        call storage_free (rmatrix%h_Kcol)
        call storage_free (rmatrix%h_Kld)
      else
        rmatrix%h_Kdiagonal = ST_NOHANDLE
        rmatrix%h_Kcol = ST_NOHANDLE
        rmatrix%h_Kld = ST_NOHANDLE
      end if

    case (LSYSSC_MATRIX1)

      if (bentries) then

        ! The structure is thrown away later so we do not care if the
        ! matrix is the owner of the structure or not at the moment
        call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
        call lsyssc_getbase_Kld (rmatrix,p_Kld)

        ! Allocate new memory for full matrix
        ihandle = ST_NOHANDLE
        call storage_new ('lsyssc_convertMatrix', 'Da', rmatrix%NEQ*rmatrix%NCOLS,&
            rmatrix%cdataType, ihandle, ST_NEWBLOCK_ZERO)

        ! Convert from format 9 to format 1 by looping over all rows and
        ! columns and copying matrix entries one by one.
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double(rmatrix, p_Ddata)
          call storage_getbase_double(ihandle, p_Ddata2)

          do iEQ = 1,rmatrix%NEQ
            do iA = p_Kld(ieq), p_Kld(ieq+1)-1
              p_Ddata2(rmatrix%NEQ*(p_Kcol(iA)-1)+iEQ) = p_Ddata(iA)
            end do
          end do

        case (ST_SINGLE)
          call lsyssc_getbase_single(rmatrix, p_Fdata)
          call storage_getbase_single(ihandle, p_Fdata2)

          do iEQ = 1,rmatrix%NEQ
            do iA = p_Kld(ieq), p_Kld(ieq+1)
              p_Fdata2(rmatrix%NEQ*(p_Kcol(iA)-1)+iEQ) = p_Fdata(iA)
            end do
          end do

        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
          call sys_halt()
        end select

        if (bentryOwner .and. (rmatrix%h_Da .ne. ST_NOHANDLE)) then
          call storage_free(rmatrix%h_Da)
        end if

        rmatrix%h_Da = ihandle
        rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_CONTENTISCOPY))

      end if

      rmatrix%cmatrixFormat = LSYSSC_MATRIX1
      rmatrix%NA = rmatrix%NEQ*rmatrix%NCOLS
      rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_STRUCTUREISCOPY))

      ! Release unused information
      if (bstrucOwner) then
        call storage_free (rmatrix%h_Kdiagonal)
        call storage_free (rmatrix%h_Kcol)
        call storage_free (rmatrix%h_Kld)
      else
        rmatrix%h_Kdiagonal = ST_NOHANDLE
        rmatrix%h_Kcol = ST_NOHANDLE
        rmatrix%h_Kld = ST_NOHANDLE
      end if

    case default
      call output_line('Cannot convert matrix!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
      call sys_halt()
    end select

  case (LSYSSC_MATRIX7)

    select case (cmatrixFormat)
    case (LSYSSC_MATRIX9)

      if (.not. bstrucOwner) then

        ! Duplicate the structure in memory so we are allowed to modify it.
        ihandle = ST_NOHANDLE
        call storage_copy (rmatrix%h_Kcol,ihandle)
        rmatrix%h_Kcol = ihandle

        ihandle = ST_NOHANDLE
        call storage_copy (rmatrix%h_Kld,ihandle)
        rmatrix%h_Kld = ihandle

        rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_STRUCTUREISCOPY))

      end if

      if ((.not. bentryOwner) .and. bentries .and. (rmatrix%h_Da .ne. ST_NOHANDLE)) then

        ! Copy the entries so we are allowed to modify them.
        ihandle = ST_NOHANDLE
        call storage_copy(rmatrix%h_Da,ihandle)
        rmatrix%h_Da = ihandle

        rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_CONTENTISCOPY))

      end if

      ! Convert from structure 7 to structure 9. Use the sortCSRxxxx
      ! routine below.
      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)

      ! Create a pointer to the diagonal
      call storage_new ('lsyssc_convertMatrix', 'Kdiagonal', &
            rmatrix%NEQ, ST_INT, rmatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT)

      ! WARNING: lsyssc_getbase_Kdiagonal does not(!) work here, because
      ! rmatrix is still in format LSYSSC_MATRIX7 and hence does not provide
      ! the array Kdiagonal ;-)
      call storage_getbase_int (rmatrix%h_Kdiagonal,p_Kdiagonal,rmatrix%NEQ)

      if ((.not. bentries) .or. (rmatrix%h_Da .eq. ST_NOHANDLE)) then

        ! No matrix entries, only resort the structure
        call lsyssc_sortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)

        rmatrix%cmatrixFormat = LSYSSC_MATRIX9

      else

        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double (rmatrix,p_Ddata)
          call lsyssc_sortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ, p_Ddata)

          rmatrix%cmatrixFormat = LSYSSC_MATRIX9

        case (ST_SINGLE)
          call lsyssc_getbase_single (rmatrix,p_Fdata)
          call lsyssc_sortCSRsingle (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ, p_Fdata)

          rmatrix%cmatrixFormat = LSYSSC_MATRIX9

        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
          call sys_halt()
        end select

      end if

    case (LSYSSC_MATRIXD)

      if (bentries .and. (rmatrix%h_Da .ne. ST_NOHANDLE)) then

        ! Convert from structure 7 to structure D by copying the
        ! diagonal to the front and reallocation of the memory
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)

          call lsyssc_getbase_double (rmatrix,p_Ddata)

          if (.not. bentryOwner) then
            ! Allocate new memory for the entries
            rmatrix%h_Da = ST_NOHANDLE
            call storage_new ('lsyssc_convertMatrix', 'Da', &
                  rmatrix%NEQ, ST_DOUBLE, rmatrix%h_Da, ST_NEWBLOCK_NOINIT)
            call storage_getbase_double (rmatrix%h_Da,p_Ddata2)
            rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_CONTENTISCOPY))
          else
            ! Destinatinon pointer points to source matrix.
            ! This overwrites the original entries which is no problem as
            ! extracting the diagonal is always a 'compression' overwriting
            ! information that is not used anymore.
            p_Ddata2 => p_Ddata
          end if

          call lsyssc_getbase_Kld (rmatrix,p_Kld)
          do i=1,rmatrix%NEQ
            p_Ddata2(i) = p_Ddata(p_Kld(i))
          end do

          rmatrix%cmatrixFormat = LSYSSC_MATRIXD

        case (ST_SINGLE)

          call lsyssc_getbase_single (rmatrix,p_Fdata)

          if (.not. bentryOwner) then
            ! Allocate new memory for the entries
            rmatrix%h_Da = ST_NOHANDLE
            call storage_new ('lsyssc_convertMatrix', 'Da', &
                  rmatrix%NEQ, ST_SINGLE, rmatrix%h_Da, ST_NEWBLOCK_NOINIT)
            call storage_getbase_single (rmatrix%h_Da,p_Fdata2)
            rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_CONTENTISCOPY))
          else
            ! Destinatinon pointer points to source matrix.
            ! This overwrites the original entries which is no problem as
            ! extracting the diagonal is always a 'compression' overwriting
            ! information that is not used anymore.
            p_Fdata2 => p_Fdata
          end if

          call lsyssc_getbase_Kld (rmatrix,p_Kld)
          do i=1,rmatrix%NEQ
            p_Fdata2(i) = p_Fdata(p_Kld(i))
          end do

          rmatrix%cmatrixFormat = LSYSSC_MATRIXD

        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
          call sys_halt()
        end select

        ! Reallocate the entries-array to have only the diagonal entries.
        ! In case we are not the owner, we have allocated new memory and thus
        ! do not need to resize it again.
        if (bentryOwner) then
          call storage_realloc ('lsyssc_convertMatrix', rmatrix%NEQ, &
                                rmatrix%h_Da, ST_NEWBLOCK_NOINIT, bentries)
        end if

        rmatrix%NA = rmatrix%NEQ

      end if

      ! Release unused information
      if (bstrucOwner) then
        call storage_free (rmatrix%h_Kcol)
        call storage_free (rmatrix%h_Kld)
      else
        rmatrix%h_Kcol = ST_NOHANDLE
        rmatrix%h_Kld = ST_NOHANDLE
      end if

    case (LSYSSC_MATRIX1)

      if (bentries) then

        ! The structure is thrown away later so we do not care if the
        ! matrix is the owner of the structure or not at the moment
        call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
        call lsyssc_getbase_Kld (rmatrix,p_Kld)

        ! Allocate new memory for full matrix
        ihandle = ST_NOHANDLE
        call storage_new ('lsyssc_convertMatrix', 'Da', rmatrix%NEQ*rmatrix%NCOLS,&
            rmatrix%cdataType, ihandle, ST_NEWBLOCK_ZERO)

        ! Convert from format 9 to format 1 by looping over all rows and
        ! columns and copying matrix entries one by one.
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double(rmatrix, p_Ddata)
          call storage_getbase_double(ihandle, p_Ddata2)

          do iEQ = 1,rmatrix%NEQ
            do iA = p_Kld(ieq), p_Kld(ieq+1)-1
              p_Ddata2(rmatrix%NCOLS*(iEQ-1)+p_Kcol(iA)) = p_Ddata(iA)
            end do
          end do

        case (ST_SINGLE)
          call lsyssc_getbase_single(rmatrix, p_Fdata)
          call storage_getbase_single(ihandle, p_Fdata2)

          do iEQ = 1,rmatrix%NEQ
            do iA = p_Kld(ieq), p_Kld(ieq+1)
              p_Fdata2(rmatrix%NCOLS*(iEQ-1)+p_Kcol(iA)) = p_Fdata(iA)
            end do
          end do

        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
          call sys_halt()
        end select

        if (bentryOwner .and. (rmatrix%h_Da .ne. ST_NOHANDLE)) then
          call storage_free(rmatrix%h_Da)
        end if

        rmatrix%h_Da = ihandle

      end if

      rmatrix%cmatrixFormat = LSYSSC_MATRIX1
      rmatrix%NA = rmatrix%NEQ*rmatrix%NCOLS

      ! Release unused information
      if (bstrucOwner) then
        call storage_free (rmatrix%h_Kdiagonal)
        call storage_free (rmatrix%h_Kcol)
        call storage_free (rmatrix%h_Kld)
      else
        rmatrix%h_Kdiagonal = ST_NOHANDLE
        rmatrix%h_Kcol = ST_NOHANDLE
        rmatrix%h_Kld = ST_NOHANDLE
      end if

    case default
      call output_line('Cannot convert matrix!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
      call sys_halt()
    end select

  case (LSYSSC_MATRIX1)

    select case (cmatrixFormat)
    case (LSYSSC_MATRIX9)

      ! Get the matrix data
      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrix,p_Ddata2)
      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrix,p_Fdata2)
      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
        call sys_halt()
      end select

      ! If the matrix is virtually transposed, NCOLS and NEQ are exchanged.
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
        nrows = rmatrix%NCOLS
        ncols = rmatrix%NEQ
      else
        nrows = rmatrix%NEQ
        ncols = rmatrix%NCOLS
      end if

      ! Count the number of nonzeroes
      NA = 0
      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        do j = 1,ncols
          do i = 1,nrows
            if (p_Ddata2((j-1)*nrows+i) .ne. 0.0_DP) NA = NA+1
          end do
        end do

      case (ST_SINGLE)
        do j = 1,ncols
          do i = 1,nrows
            if (p_Fdata2((j-1)*nrows+i) .ne. 0.0_SP) NA = NA+1
          end do
        end do
      end select

      ! Allocate Kcol, Kld, Da,Kdiagonal
      call storage_new ('lsyssc_convertMatrix', 'Kcol', &
            NA, ST_INT, h_Kcol, ST_NEWBLOCK_NOINIT)
      call storage_new ('lsyssc_convertMatrix', 'Kld', &
            nrows+1, ST_INT, h_Kld, ST_NEWBLOCK_NOINIT)
      call storage_new ('lsyssc_convertMatrix', 'Kdiagonal', &
            nrows, ST_INT, h_Kdiagonal, ST_NEWBLOCK_NOINIT)

      call storage_getbase_int (h_Kcol,p_Kcol)
      call storage_getbase_int (h_Kld,p_Kld)
      call storage_getbase_int (h_Kdiagonal,p_Kdiagonal)

      ! Transfer the entries -- rowwise.
      ! If we do not have to transfer the entries, simply set up the structure.
      if (bentries) then

        ! Allocate memory for the entries
        call storage_new ('lsyssc_convertMatrix', 'Da', &
              NA, rmatrix%cdataType, h_Da, ST_NEWBLOCK_NOINIT)
        call storage_getbase_double (h_Da,p_Ddata)

        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          iA = 0
          do i = 1,nrows
            p_Kld(i) = iA+1
            do j = 1,ncols
              if (p_Ddata2((j-1)*nrows+i) .ne. 0.0_DP) then
                iA = iA + 1
                p_Ddata(iA) = p_Ddata2((j-1)*nrows+i)
                p_Kcol(iA) = j
              end if

              ! p_Kdiagonal is set to the first entry on or above the diagonal.
              if (i .eq. j) p_Kdiagonal(i) = iA
            end do
            if (i .gt. ncols) p_Kdiagonal(i) = iA  ! Kdiagonal was not set
          end do
          p_Kld(nrows+1) = NA+1

        case (ST_SINGLE)
          iA = 0
          do i = 1,nrows
            p_Kld(i) = iA+1
            do j = 1,ncols
              if (p_Fdata2((j-1)*nrows+i) .ne. 0.0_SP) then
                iA = iA + 1
                p_Fdata(iA) = p_Fdata2((j-1)*nrows+i)
                p_Kcol(iA) = j
              end if

              ! p_Kdiagonal is set to the first entry on or above the diagonal.
              if (i .eq. j) p_Kdiagonal(i) = iA
            end do
            if (i .gt. ncols) p_Kdiagonal(i) = iA  ! Kdiagonal was not set
          end do
          p_Kld(nrows+1) = NA+1
        end select

        ! Release the content if it belongs to this matrix.
        ! Replace by the new content.
        if (bentryOwner) then
          call storage_free (rmatrix%h_Da)
        end if

        rmatrix%h_Da = h_Da

      else

        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          iA = 0
          do i = 1,nrows
            p_Kld(i) = iA+1
            do j = 1,ncols
              if (p_Ddata2((j-1)*nrows+i) .ne. 0.0_DP) then
                iA = iA + 1
                p_Kcol(iA) = j
              end if

              ! p_Kdiagonal is set to the first entry on or above the diagonal.
              if (i .eq. j) p_Kdiagonal(i) = iA
            end do
            if (i .gt. ncols) p_Kdiagonal(i) = iA  ! Kdiagonal was not set
          end do
          p_Kld(nrows+1) = NA+1

        case (ST_SINGLE)
          iA = 0
          do i = 1,nrows
            p_Kld(i) = iA+1
            do j = 1,ncols
              if (p_Fdata2((j-1)*nrows+i) .ne. 0.0_SP) then
                iA = iA + 1
                p_Kcol(iA) = j
              end if

              ! p_Kdiagonal is set to the first entry on or above the diagonal.
              if (i .eq. j) p_Kdiagonal(i) = iA
            end do
            if (i .gt. ncols) p_Kdiagonal(i) = iA  ! Kdiagonal was not set
          end do
          p_Kld(nrows+1) = NA+1
        end select

        ! Do not touch the content. The result is left in an undefined state!

      end if

      ! Switch the matrix format
      rmatrix%cmatrixFormat = LSYSSC_MATRIX9

      rmatrix%h_Kcol = h_Kcol
      rmatrix%h_Kld = h_Kld
      rmatrix%h_Kdiagonal = h_Kdiagonal
      rmatrix%NA = NA

      ! Switch off any duplication flag; the content now belongs to this matrix.
      rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_ISCOPY))

    end select

  case (LSYSSC_MATRIX9INTL)

    select case (cmatrixFormat)
    case (LSYSSC_MATRIX7INTL)

      if (.not. bstrucOwner) then

        ! Duplicate the structure in memory so we are allowed to modify it.
        ihandle = ST_NOHANDLE
        call storage_copy (rmatrix%h_Kcol,ihandle)
        rmatrix%h_Kcol = ihandle

        ihandle = ST_NOHANDLE
        call storage_copy (rmatrix%h_Kld,ihandle)
        rmatrix%h_Kld = ihandle

        ! Kdiagonal is thrown away later, we do not need to copy it.

        rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_STRUCTUREISCOPY))

      end if

      if ((.not. bentryOwner) .and. bentries .and. (rmatrix%h_Da .ne. ST_NOHANDLE)) then

        ! Copy the entries so we are allowed to modify them.
        ihandle = ST_NOHANDLE
        call storage_copy(rmatrix%h_Da,ihandle)
        rmatrix%h_Da = ihandle

        rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_CONTENTISCOPY))

      end if

      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)
      call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)

      ! Check that the matrix can be converted. There is a format error
      ! if there is no diagonal element.
      do i=1,rmatrix%NEQ
        if (p_Kcol(p_Kdiagonal(i)) .ne. i) then
          call output_line('Incompatible Format-9 matrix!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
          call sys_halt()
        end if
      end do

      if ((.not. bentries) .or. (rmatrix%h_Da .eq. ST_NOHANDLE)) then

        ! No matrix entries, only resort the structure
        call lsyssc_unsortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)

        ! Release diagonal pointer if it belongs to us
        if (bstrucOwner) then
          call storage_free (rmatrix%h_Kdiagonal)
        else
          rmatrix%h_Kdiagonal = ST_NOHANDLE
        end if

        rmatrix%cmatrixFormat = LSYSSC_MATRIX7INTL

      else

        ! Convert from structure 9 to structure 7. Use the sortCSRxxxx
        ! routine below.
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)

          call lsyssc_getbase_double (rmatrix,p_Ddata)
          select case (rmatrix%cinterleavematrixFormat)
          case (LSYSSC_MATRIX1)
            call lsyssc_unsortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal,&
                & rmatrix%NEQ, p_Ddata, rmatrix%NVAR*rmatrix%NVAR)

          case (LSYSSC_MATRIXD)
            call lsyssc_unsortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal,&
                & rmatrix%NEQ, p_Ddata, rmatrix%NVAR)

          case default
            call output_line('Unsupported interleave matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
            call sys_halt()
          end select

        case (ST_SINGLE)

          call lsyssc_getbase_single (rmatrix,p_Fdata)
          select case (rmatrix%cinterleavematrixFormat)
          case (LSYSSC_MATRIX1)
            call lsyssc_unsortCSRsingle (p_Kcol, p_Kld, p_Kdiagonal,&
                & rmatrix%NEQ, p_Fdata, rmatrix%NVAR*rmatrix%NVAR)

          case (LSYSSC_MATRIXD)
            call lsyssc_unsortCSRsingle (p_Kcol, p_Kld, p_Kdiagonal,&
                & rmatrix%NEQ, p_Fdata, rmatrix%NVAR)

          case default
            call output_line('Unsupported interleave matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
            call sys_halt()
          end select

        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
          call sys_halt()
        end select

        ! Release diagonal pointer if it belongs to us
        if (bstrucOwner) then
          call storage_free (rmatrix%h_Kdiagonal)
        else
          rmatrix%h_Kdiagonal = ST_NOHANDLE
        end if

        rmatrix%cmatrixFormat = LSYSSC_MATRIX7INTL
      end if

    case default
      call output_line('Cannot convert matrix!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
      call sys_halt()
    end select

  case (LSYSSC_MATRIX7INTL)

    select case (cmatrixFormat)
    case (LSYSSC_MATRIX9INTL)

      if (.not. bstrucOwner) then

        ! Duplicate the structure in memory so we are allowed to modify it.
        ihandle = ST_NOHANDLE
        call storage_copy (rmatrix%h_Kcol,ihandle)
        rmatrix%h_Kcol = ihandle

        ihandle = ST_NOHANDLE
        call storage_copy (rmatrix%h_Kld,ihandle)
        rmatrix%h_Kld = ihandle

        rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_STRUCTUREISCOPY))

      end if

      if ((.not. bentryOwner) .and. bentries .and. (rmatrix%h_Da .ne. ST_NOHANDLE)) then

        ! Copy the entries so we are allowed to modify them.
        ihandle = ST_NOHANDLE
        call storage_copy(rmatrix%h_Da,ihandle)
        rmatrix%h_Da = ihandle

        rmatrix%imatrixSpec = iand(rmatrix%imatrixSpec,not(LSYSSC_MSPEC_CONTENTISCOPY))

      end if

      ! Convert from structure 7 to structure 9. Use the sortCSRxxxx
      ! routine below.
      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)

      ! Create a pointer to the diagonal
      call storage_new ('lsyssc_convertMatrix', 'Kdiagonal', &
          rmatrix%NEQ, ST_INT, rmatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT)
      call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)

      if ((.not. bentries) .or. (rmatrix%h_Da .eq. ST_NOHANDLE)) then

        ! No matrix entries, only resort the structure
        call lsyssc_sortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)

        rmatrix%cmatrixFormat = LSYSSC_MATRIX9INTL

      else

        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double (rmatrix,p_Ddata)
          select case (rmatrix%cinterleavematrixFormat)
          case (LSYSSC_MATRIX1)
            call lsyssc_sortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal,&
                & rmatrix%NEQ, p_Ddata, rmatrix%NVAR*rmatrix%NVAR)

          case (LSYSSC_MATRIXD)
            call lsyssc_sortCSRdouble (p_Kcol, p_Kld, p_Kdiagonal,&
                & rmatrix%NEQ, p_Ddata, rmatrix%NVAR)

          case default
            call output_line('Unsupported interleave matrix format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
            call sys_halt()
          end select

          rmatrix%cmatrixFormat = LSYSSC_MATRIX9INTL

        case (ST_SINGLE)
          call lsyssc_getbase_single (rmatrix,p_Fdata)
          select case (rmatrix%cinterleavematrixFormat)
          case (LSYSSC_MATRIX1)
            call lsyssc_sortCSRsingle (p_Kcol, p_Kld, p_Kdiagonal,&
                & rmatrix%NEQ, p_Fdata, rmatrix%NVAR*rmatrix%NVAR)

          case (LSYSSC_MATRIXD)
            call lsyssc_sortCSRsingle (p_Kcol, p_Kld, p_Kdiagonal,&
                & rmatrix%NEQ, p_Fdata, rmatrix%NVAR)

          case default
            call output_line('Unsupported interleave matrix format!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
            call sys_halt()
          end select

          rmatrix%cmatrixFormat = LSYSSC_MATRIX9INTL

        case default
          call output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
          call sys_halt()
        end select

      end if

    case default
      call output_line('Cannot convert matrix!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
      call sys_halt()
    end select

  case default
    call output_line('Cannot convert matrix!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_convertMatrix')
    call sys_halt()
  end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_copyH2D_Vector(rvector, bclear, btranspose, istream)

!<description>
    ! This subroutine copies the data of a vector from the host memory
    ! to the memory of the coprocessor device. If no device is
    ! available, then an error is thrown.
!</description>

!<input>
    ! Scalar vector
    type(t_vectorScalar), intent(in) :: rvector

    ! Switch: If TRUE, then memory is allocated on the device and
    ! cleared. If FALSE, then memory is allocated on the device and
    ! the content from the host memory is transferred.
    logical :: bclear

    ! OPTIONAL: if TRUE then the memory is transposed.
    logical, intent(in), optional :: btranspose

    ! OPTIONAL: stream for asynchronious transfer.
    ! If istream is present and if asynchroneous transfer is supported
    ! then all memory transfers are carried out asynchroneously
    integer(I64), intent(in), optional :: istream
!</input>

!</subroutine>

    ! Clear the vector?
    if (bclear) then
      call storage_clearMemoryOnDevice(rvector%h_Ddata)
    else
      ! No, we have to copy the content from host to device memory
      call storage_syncMemoryHostDevice(rvector%h_Ddata,&
          ST_SYNCBLOCK_COPY_H2D, btranspose, istream)
    end if

  end subroutine lsyssc_copyH2D_Vector

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_copyD2H_Vector(rvector, bclear, btranspose, istream)

!<description>
    ! This subroutine copies the data of a vector from the memory of
    ! the coprocessor device to the host memory. If no device is
    ! available, then an error is thrown.
!</description>

!<input>
    ! Scalar vector
    type(t_vectorScalar), intent(in) :: rvector

    ! Switch: If TRUE, then memory is allocated on the device and
    ! cleared. If FALSE, then memory is allocated on the device and
    ! the content from the host memory is transferred.
    logical :: bclear

    ! OPTIONAL: if TRUE then the memory is transposed.
    logical, intent(in), optional :: btranspose

    ! OPTIONAL: stream for asynchronious transfer.
    ! If istream is present and if asynchroneous transfer is supported
    ! then all memory transfers are carried out asynchroneously
    integer(I64), intent(in), optional :: istream
!</input>

!</subroutine>

    ! Clear the vector?
    if (bclear) then
      ! Yes, we can simply overwrite the content in host memory
      call storage_syncMemoryHostDevice(rvector%h_Ddata,&
          ST_SYNCBLOCK_COPY_D2H, btranspose, istream)
    else
      ! No, we have to copy-add the content from device to host memory
      call storage_syncMemoryHostDevice(rvector%h_Ddata,&
          ST_SYNCBLOCK_ACCUMULATE_D2H, btranspose, istream)
    end if

  end subroutine lsyssc_copyD2H_Vector

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_copyH2D_Matrix(rmatrix, bclear, btranspose, istream)

!<description>
    ! This subroutine copies the data of a matrix from the host memory
    ! to the memory of the coprocessor device. If no device is
    ! available, then an error is thrown.
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! Switch: If TRUE, then memory is allocated on the device and
    ! cleared. If FALSE, then memory is allocated on the device and
    ! the content from the host memory is transferred.
    logical :: bclear

    ! OPTIONAL: if TRUE then the memory is transposed.
    logical, intent(in), optional :: btranspose

    ! OPTIONAL: stream for asynchronious transfer.
    ! If istream is present and if asynchroneous transfer is supported
    ! then all memory transfers are carried out asynchroneously
    integer(I64), intent(in), optional :: istream
!</input>

!</subroutine>

    ! Clear the matrix?
    if (bclear) then
      call storage_clearMemoryOnDevice(rmatrix%h_Da)
    else
      ! No, we have to copy the content from host to device memory
      call storage_syncMemoryHostDevice(rmatrix%h_Da,&
          ST_SYNCBLOCK_COPY_H2D, btranspose, istream)
    end if

  end subroutine lsyssc_copyH2D_Matrix

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_copyD2H_Matrix(rmatrix, bclear, btranspose, istream)

!<description>
    ! This subroutine copies the data of a matrix from the memory of
    ! the coprocessor device to the host memory. If no device is
    ! available, then an error is thrown.
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! Switch: If TRUE, then memory is allocated on the device and
    ! cleared. If FALSE, then memory is allocated on the device and
    ! the content from the host memory is transferred.
    logical :: bclear

    ! OPTIONAL: if TRUE then the memory is transposed.
    logical, intent(in), optional :: btranspose

    ! OPTIONAL: stream for asynchronious transfer.
    ! If istream is present and if asynchroneous transfer is supported
    ! then all memory transfers are carried out asynchroneously
    integer(I64), intent(in), optional :: istream
!</input>

!</subroutine>

    ! Clear the vector?
    if (bclear) then
      ! Yes, we can simply overwrite the content in host memory
      call storage_syncMemoryHostDevice(rmatrix%h_Da,&
          ST_SYNCBLOCK_COPY_D2H, btranspose, istream)
    else
      ! No, we have to copy-add the content from device to host memory
      call storage_syncMemoryHostDevice(rmatrix%h_Da,&
          ST_SYNCBLOCK_ACCUMULATE_D2H, btranspose, istream)
    end if

  end subroutine lsyssc_copyD2H_Matrix

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_createMatrixSymmPart(rmatrix, dscale, rmatrixDest)

!<description>
    ! This subroutine computes the symmetric parts dscale*(A+A^T) of the
    ! matrix rmatrix. If rmatrixDest is not present then the resulting
    ! matrix is stored in rmatrix.
!</description>

!<input>
    ! Scaling parameter
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! Source matrix
    type(t_matrixScalar), intent(inout), target :: rmatrix

    ! OPTIONAL: destination matrix
    type(t_matrixScalar), intent(inout), target, optional :: rmatrixDest
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixScalar), pointer :: p_rmatrix
    real(DP), dimension(:), pointer :: p_Da,p_Db
    real(SP), dimension(:), pointer :: p_Fa,p_Fb
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Kdiagonal,p_Ksep
    integer :: h_Ksep

    ! Check if matrix is a square matrix
    if (rmatrix%NEQ .ne. rmatrix%NCOLS) then
      call output_line('Matrix is not a square matrix!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_createMatrixSymmPart')
      call sys_halt()
    end if

    ! Check if destination matrix is present
    if (present(rmatrixDest)) then
      ! Check if source and destination matrix are compatible
      call lsyssc_isMatrixCompatible(rmatrix, rmatrixDest)
      p_rmatrix => rmatrixDest
    else
      p_rmatrix => rmatrix
    end if

    ! What matrix format are we?
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIXD)
      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call lalg_scaleVector(p_Da, 2._DP*dscale)
      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call lalg_scaleVector(p_Fa, real(2._SP*dscale,SP))
      end select

    case (LSYSSC_MATRIX1)
      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call lsyssc_getbase_double(p_rmatrix, p_Db)
        call calcSymmPartMat1Dble(rmatrix%NEQ, rmatrix%NCOLS, dscale, p_Da, p_Db)
      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call lsyssc_getbase_single(p_rmatrix, p_Fb)
        call calcSymmPartMat1Sngl(rmatrix%NEQ, rmatrix%NCOLS, real(dscale,SP), p_Fa, p_Fb)
      end select

    case (LSYSSC_MATRIX7)
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrix%NEQ+1)
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call lsyssc_getbase_double(p_rmatrix, p_Db)
        call calcSymmPartMat7Dble(p_Kld, p_Kcol, p_Ksep, 1, dscale, p_Da, p_Db)
      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call lsyssc_getbase_single(p_rmatrix, p_Fb)
        call calcSymmPartMat7Sngl(p_Kld, p_Kcol, p_Ksep, 1, real(dscale,SP), p_Fa, p_Fb)
      end select

      ! Release diagonal separator
      call storage_free(h_Ksep)

    case (LSYSSC_MATRIX7INTL)
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrix%NEQ+1)
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      ! What interleave matrix format are we?
      select case (rmatrix%cinterleavematrixFormat)

      case (LSYSSC_MATRIXD)
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double(rmatrix, p_Da)
          call lsyssc_getbase_double(p_rmatrix, p_Db)
          call calcSymmPartMat7Dble(p_Kld, p_Kcol, p_Ksep,&
                                    rmatrix%NVAR, dscale, p_Da, p_Db)
        case (ST_SINGLE)
          call lsyssc_getbase_single(rmatrix, p_Fa)
          call lsyssc_getbase_single(p_rmatrix, p_Fb)
          call calcSymmPartMat7Sngl(p_Kld, p_Kcol, p_Ksep,&
                                    rmatrix%NVAR, real(dscale,SP), p_Fa, p_Fb)
        end select

      case (LSYSSC_MATRIX1)
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double(rmatrix, p_Da)
          call lsyssc_getbase_double(p_rmatrix, p_Db)
          call calcSymmPartMat7Dble(p_Kld, p_Kcol, p_Ksep,&
                                    rmatrix%NVAR*rmatrix%NVAR, dscale, p_Da, p_Db)
        case (ST_SINGLE)
          call lsyssc_getbase_single(rmatrix, p_Fa)
          call lsyssc_getbase_single(p_rmatrix, p_Fb)
          call calcSymmPartMat7Sngl(p_Kld, p_Kcol, p_Ksep,&
                                    rmatrix%NVAR*rmatrix%NVAR, real(dscale,SP), p_Fa, p_Fb)
        end select

      case default
        call output_line('Invalid interleave matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_createMatrixSymmPart')
        call sys_halt()
      end select

      ! Release diagonal separator
      call storage_free(h_Ksep)

    case (LSYSSC_MATRIX9)
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrix%NEQ+1)
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call lsyssc_getbase_double(p_rmatrix, p_Db)
        call calcSymmPartMat9Dble(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                  1, dscale, p_Da, p_Db)
      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call lsyssc_getbase_single(p_rmatrix, p_Fb)
        call calcSymmPartMat9Sngl(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                  1, real(dscale,SP), p_Fa, p_Fb)
      end select

      ! Release diagonal separator
      call storage_free(h_Ksep)

    case (LSYSSC_MATRIX9INTL)
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrix%NEQ+1)
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

      ! What interleave matrix format are we?
      select case (rmatrix%cinterleavematrixFormat)

      case (LSYSSC_MATRIXD)
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double(rmatrix, p_Da)
          call lsyssc_getbase_double(p_rmatrix, p_Db)
          call calcSymmPartMat9Dble(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rmatrix%NVAR, dscale, p_Da, p_Db)
        case (ST_SINGLE)
          call lsyssc_getbase_single(rmatrix, p_Fa)
          call lsyssc_getbase_single(p_rmatrix, p_Fb)
          call calcSymmPartMat9Sngl(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rmatrix%NVAR, real(dscale,SP), p_Fa, p_Fb)
        end select

      case (LSYSSC_MATRIX1)
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double(rmatrix, p_Da)
          call lsyssc_getbase_double(p_rmatrix, p_Db)
          call calcSymmPartMat9Dble(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rmatrix%NVAR*rmatrix%NVAR, dscale, p_Da, p_Db)
        case (ST_SINGLE)
          call lsyssc_getbase_single(rmatrix, p_Fa)
          call lsyssc_getbase_single(p_rmatrix, p_Fb)
          call calcSymmPartMat9Sngl(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                    rmatrix%NVAR*rmatrix%NVAR, real(dscale,SP), p_Fa, p_Fb)
        end select

      case default
        call output_line('Invalid interleave matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_createMatrixSymmPart')
        call sys_halt()
      end select

      ! Release diagonal separator
      call storage_free(h_Ksep)

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_WARNING,OU_MODE_STD,'lsyssc_createMatrixSymmPart')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Compute symmetric part for dense matrix stored in matrix format 1

    subroutine calcSymmPartMat1Dble(NEQ, NCOLS, dscale, Da, Db)

      ! input parameters
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NCOLS

      ! input/output parameters
      real(DP), dimension(NEQ,NCOLS), intent(inout) :: Da,Db

      ! local variables
      integer :: i,j

      ! Loop over all rows
      do i = 1, NEQ
        Db(i,i) = dscale * 2._DP * Da(i,i)
        do j = i+1, NCOLS
          Db(i,j) = dscale * (Da(i,j) + Da(j,i))
          Db(j,i) = Db(i,j)
        end do
      end do

    end subroutine calcSymmPartMat1Dble

    !**************************************************************
    ! Compute symmetric part for dense matrix stored in matrix format 1

    subroutine calcSymmPartMat1Sngl(NEQ, NCOLS, fscale, Fa, Fb)

      ! input parameters
      real(SP), intent(in) :: fscale
      integer, intent(in) :: NEQ,NCOLS

      ! input/output parameters
      real(SP), dimension(NEQ,NCOLS), intent(inout) :: Fa,Fb

      ! local variables
      integer :: i,j

      ! Loop over all rows
      do i = 1, NEQ
        Fb(i,i) = fscale * 2._SP * Fa(i,i)
        do j = i+1, NCOLS
          Fb(i,j) = fscale * (Fa(i,j) + Fa(j,i))
          Fb(j,i) = Fb(i,j)
        end do
      end do

    end subroutine calcSymmPartMat1Sngl

    !**************************************************************
    ! Compute symmetric part for CRS matrix stored in matrix format 7

    subroutine calcSymmPartMat7Dble(Kld, Kcol, Ksep, NVAR, dscale, Da, Db)

      ! input parameters
      integer, dimension(:), intent(in) :: Kld,Kcol
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,*), intent(inout) :: Da,Db
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      integer :: i,ii,ij,j,ji

      ! Loop over all rows
      do i = 1, size(Kld)-1

        ! Diagonal entry
        ii = Kld(i)
        Db(:,ii) = dscale * 2._DP * Da(:,ii)

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J, that is, the upper
        ! right triangular part of the matrix is traversed
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          Db(:,ij) = dscale * (Da(:,ij) + Da(:,ji))
          Db(:,ji) = Db(:,ij)
        end do
      end do

    end subroutine calcSymmPartMat7Dble

    !**************************************************************
    ! Compute symmetric part for CRS matrix stored in matrix format 7

    subroutine calcSymmPartMat7Sngl(Kld, Kcol, Ksep, NVAR, fscale, Fa, Fb)

      ! input parameters
      integer, dimension(:), intent(in) :: Kld,Kcol
      real(SP), intent(in) :: fscale
      integer, intent(in) :: NVAR

      ! input/output parameters
      real(SP), dimension(NVAR,*), intent(inout) :: Fa,Fb
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      integer :: i,ii,ij,j,ji

      ! Loop over all rows
      do i = 1, size(Kld)-1

        ! Diagonal entry
        ii = Kld(i)
        Fb(:,ii) = fscale * 2._SP * Fa(:,ii)

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J, that is, the upper
        ! right triangular part of the matrix is traversed
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          Fb(:,ij) = fscale * (Fa(:,ij) + Fa(:,ji))
          Fb(:,ji) = Fb(:,ij)
        end do
      end do

    end subroutine calcSymmPartMat7Sngl

    !**************************************************************
    ! Compute symmetric part for CRS matrix stored in matrix format 9

    subroutine calcSymmPartMat9Dble(Kld, Kcol, Kdiagonal, Ksep, NVAR, dscale, Da, Db)

      ! input parameters
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,*), intent(inout) :: Da,Db
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      integer :: i,ii,ij,j,ji

      ! Loop over all rows
      do i = 1, size(Kld)-1

        ! Diagonal entry
        ii = Kdiagonal(i)
        Db(:,ii) = dscale * 2._DP * Da(:,ii)

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J, that is, the upper
        ! right triangular part of the matrix is traversed
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          Db(:,ij) = dscale * (Da(:,ij) + Da(:,ji))
          Db(:,ji) = Db(:,ij)
        end do
      end do

    end subroutine calcSymmPartMat9Dble

    !**************************************************************
    ! Compute symmetric part for CRS matrix stored in matrix format 9

    subroutine calcSymmPartMat9Sngl(Kld, Kcol, Kdiagonal, Ksep, NVAR, fscale, Fa, Fb)

      ! input parameters
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      real(SP), intent(in) :: fscale
      integer, intent(in) :: NVAR

      ! input/output parameters
      real(SP), dimension(NVAR,*), intent(inout) :: Fa,Fb
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      integer :: i,ii,ij,j,ji

      ! Loop over all rows
      do i = 1, size(Kld)-1

        ! Diagonal entry
        ii = Kdiagonal(i)
        Fb(:,ii) = fscale * 2._SP * Fa(:,ii)

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J, that is, the upper
        ! right triangular part of the matrix is traversed
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          Fb(:,ij) = fscale * (Fa(:,ij) + Fa(:,ji))
          Fb(:,ji) = Fb(:,ij)
        end do
      end do

    end subroutine calcSymmPartMat9Sngl

  end subroutine lsyssc_createMatrixSymmPart

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_createMatrixAsymmPart(rmatrix, dscale, rmatrixDest)

!<description>
    ! This subroutine computes the antisymmetric parts dscale*(A-A^T)
    ! of the matrix rmatrix. If rmatrixDest is not present then the
    ! resulting matrix is stored in rmatrix.
!</description>

!<input>
    ! Scaling parameter
    real(DP), intent(in) :: dscale
!</input>

!<inputoutput>
    ! Source matrix
    type(t_matrixScalar), intent(inout), target :: rmatrix

    ! OPTIONAL: destination matrix
    type(t_matrixScalar), intent(inout), target, optional :: rmatrixDest
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixScalar), pointer :: p_rmatrix
    real(DP), dimension(:), pointer :: p_Da,p_Db
    real(SP), dimension(:), pointer :: p_Fa,p_Fb
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Kdiagonal,p_Ksep
    integer :: h_Ksep

    ! Check if matrix is a square matrix
    if (rmatrix%NEQ .ne. rmatrix%NCOLS) then
      call output_line('Matrix is not a square matrix!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_createMatrixAsymmPart')
      call sys_halt()
    end if

    ! Check if destination matrix is present
    if (present(rmatrixDest)) then
      ! Check if source and destination matrix are compatible
      call lsyssc_isMatrixCompatible(rmatrix, rmatrixDest)
      p_rmatrix => rmatrixDest
    else
      p_rmatrix => rmatrix
    end if

    ! What matrix format are we?
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIXD)
      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call lalg_scaleVector(p_Da, 2._DP*dscale)
      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call lalg_scaleVector(p_Fa, real(2._SP*dscale,SP))
      end select

    case (LSYSSC_MATRIX1)
      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call lsyssc_getbase_double(p_rmatrix, p_Db)
        call calcAsymmPartMat1Dble(rmatrix%NEQ, rmatrix%NCOLS, dscale, p_Da, p_Db)
      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call lsyssc_getbase_single(p_rmatrix, p_Fb)
        call calcAsymmPartMat1Sngl(rmatrix%NEQ, rmatrix%NCOLS, real(dscale,SP), p_Fa, p_Fb)
      end select

    case (LSYSSC_MATRIX7)
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrix%NEQ+1)
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call lsyssc_getbase_double(p_rmatrix, p_Db)
        call calcAsymmPartMat7Dble(p_Kld, p_Kcol, p_Ksep, 1, dscale, p_Da, p_Db)
      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call lsyssc_getbase_single(p_rmatrix, p_Fb)
        call calcAsymmPartMat7Sngl(p_Kld, p_Kcol, p_Ksep, 1, real(dscale,SP), p_Fa, p_Fb)
      end select

      ! Release diagonal separator
      call storage_free(h_Ksep)

    case (LSYSSC_MATRIX7INTL)
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrix%NEQ+1)
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      ! What interleave matrix format are we?
      select case (rmatrix%cinterleavematrixFormat)

      case (LSYSSC_MATRIXD)
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double(rmatrix, p_Da)
          call lsyssc_getbase_double(p_rmatrix, p_Db)
          call calcAsymmPartMat7Dble(p_Kld, p_Kcol, p_Ksep,&
                                     rmatrix%NVAR, dscale, p_Da, p_Db)
        case (ST_SINGLE)
          call lsyssc_getbase_single(rmatrix, p_Fa)
          call lsyssc_getbase_single(p_rmatrix, p_Fb)
          call calcAsymmPartMat7Sngl(p_Kld, p_Kcol, p_Ksep,&
                                     rmatrix%NVAR, real(dscale,SP), p_Fa, p_Fb)
        end select

      case (LSYSSC_MATRIX1)
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double(rmatrix, p_Da)
          call lsyssc_getbase_double(p_rmatrix, p_Db)
          call calcAsymmPartMat7Dble(p_Kld, p_Kcol, p_Ksep,&
                                     rmatrix%NVAR*rmatrix%NVAR, dscale, p_Da, p_Db)
        case (ST_SINGLE)
          call lsyssc_getbase_single(rmatrix, p_Fa)
          call lsyssc_getbase_single(p_rmatrix, p_Fb)
          call calcAsymmPartMat7Sngl(p_Kld, p_Kcol, p_Ksep,&
                                     rmatrix%NVAR*rmatrix%NVAR, real(dscale,SP), p_Fa, p_Fb)
        end select

      case default
        call output_line('Invalid interleave matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_createMatrixAsymmPart')
        call sys_halt()
      end select

      ! Release diagonal separator
      call storage_free(h_Ksep)

    case (LSYSSC_MATRIX9)
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrix%NEQ+1)
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        call lsyssc_getbase_double(p_rmatrix, p_Db)
        call calcAsymmPartMat9Dble(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                   1, dscale, p_Da, p_Db)
      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        call lsyssc_getbase_single(p_rmatrix, p_Fb)
        call calcAsymmPartMat9Sngl(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                   1, real(dscale,SP), p_Fa, p_Fb)
      end select

      ! Release diagonal separator
      call storage_free(h_Ksep)

    case (LSYSSC_MATRIX9INTL)
      ! Create diagonal separator
      h_Ksep = ST_NOHANDLE
      call storage_copy(rmatrix%h_Kld, h_Ksep)
      call storage_getbase_int(h_Ksep, p_Ksep, rmatrix%NEQ+1)
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

      ! What interleave matrix format are we?
      select case (rmatrix%cinterleavematrixFormat)

      case (LSYSSC_MATRIXD)
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double(rmatrix, p_Da)
          call lsyssc_getbase_double(p_rmatrix, p_Db)
          call calcAsymmPartMat9Dble(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                     rmatrix%NVAR, dscale, p_Da, p_Db)
        case (ST_SINGLE)
          call lsyssc_getbase_single(rmatrix, p_Fa)
          call lsyssc_getbase_single(p_rmatrix, p_Fb)
          call calcAsymmPartMat9Sngl(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                     rmatrix%NVAR, real(dscale,SP), p_Fa, p_Fb)
        end select

      case (LSYSSC_MATRIX1)
        select case (rmatrix%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double(rmatrix, p_Da)
          call lsyssc_getbase_double(p_rmatrix, p_Db)
          call calcAsymmPartMat9Dble(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                     rmatrix%NVAR*rmatrix%NVAR, dscale, p_Da, p_Db)
        case (ST_SINGLE)
          call lsyssc_getbase_single(rmatrix, p_Fa)
          call lsyssc_getbase_single(p_rmatrix, p_Fb)
          call calcAsymmPartMat9Sngl(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
                                     rmatrix%NVAR*rmatrix%NVAR, real(dscale,SP), p_Fa, p_Fb)
        end select

      case default
        call output_line('Invalid interleave matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_createMatrixAsymmPart')
        call sys_halt()
      end select

      ! Release diagonal separator
      call storage_free(h_Ksep)

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_WARNING,OU_MODE_STD,'lsyssc_createMatrixAsymmPart')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Compute symmetric part for dense matrix stored in matrix format 1

    subroutine calcAsymmPartMat1Dble(NEQ, NCOLS, dscale, Da, Db)

      ! input parameters
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NEQ,NCOLS

      ! input/output parameters
      real(DP), dimension(NEQ,NCOLS), intent(inout) :: Da,Db

      ! local variables
      integer :: i,j

      ! Loop over all rows
      do i = 1, NEQ
        Db(i,i) = 0._DP
        do j = i+1, NCOLS
          Db(i,j) = dscale * (Da(i,j) - Da(j,i))
          Db(j,i) = Db(i,j)
        end do
      end do

    end subroutine calcAsymmPartMat1Dble

    !**************************************************************
    ! Compute symmetric part for dense matrix stored in matrix format 1

    subroutine calcAsymmPartMat1Sngl(NEQ, NCOLS, fscale, Fa, Fb)

      ! input parameters
      real(SP), intent(in) :: fscale
      integer, intent(in) :: NEQ,NCOLS

      ! input/output parameters
      real(SP), dimension(NEQ,NCOLS), intent(inout) :: Fa,Fb

      ! local variables
      integer :: i,j

      ! Loop over all rows
      do i = 1, NEQ
        Fb(i,i) = 0._SP
        do j = i+1, NCOLS
          Fb(i,j) = fscale * (Fa(i,j) - Fa(j,i))
          Fb(j,i) = Fb(i,j)
        end do
      end do

    end subroutine calcAsymmPartMat1Sngl

    !**************************************************************
    ! Compute symmetric part for CRS matrix stored in matrix format 7

    subroutine calcAsymmPartMat7Dble(Kld, Kcol, Ksep, NVAR, dscale, Da, Db)

      ! input parameters
      integer, dimension(:), intent(in) :: Kld,Kcol
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,*), intent(inout) :: Da,Db
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      integer :: i,ii,ij,j,ji

      ! Loop over all rows
      do i = 1, size(Kld)-1

        ! Diagonal entry
        ii = Kld(i)
        Db(:,ii) = 0._DP

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J, that is, the upper
        ! right triangular part of the matrix is traversed
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          Db(:,ij) = dscale * (Da(:,ij) - Da(:,ji))
          Db(:,ji) = Db(:,ij)
        end do
      end do

    end subroutine calcAsymmPartMat7Dble

    !**************************************************************
    ! Compute symmetric part for CRS matrix stored in matrix format 7

    subroutine calcAsymmPartMat7Sngl(Kld, Kcol, Ksep, NVAR, fscale, Fa, Fb)

      ! input parameters
      integer, dimension(:), intent(in) :: Kld,Kcol
      real(SP), intent(in) :: fscale
      integer, intent(in) :: NVAR

      ! input/output parameters
      real(SP), dimension(NVAR,*), intent(inout) :: Fa,Fb
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      integer :: i,ii,ij,j,ji

      ! Loop over all rows
      do i = 1, size(Kld)-1

        ! Diagonal entry
        ii = Kld(i)
        Fb(:,ii) = 0._SP

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J, that is, the upper
        ! right triangular part of the matrix is traversed
        do ij = Ksep(i)+1, Kld(i+1)-1

          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

          Fb(:,ij) = fscale * (Fa(:,ij) - Fa(:,ji))
          Fb(:,ji) = Fb(:,ij)
        end do
      end do

    end subroutine calcAsymmPartMat7Sngl

    !**************************************************************
    ! Compute symmetric part for CRS matrix stored in matrix format 9

    subroutine calcAsymmPartMat9Dble(Kld, Kcol, Kdiagonal, Ksep, NVAR, dscale, Da, Db)

      ! input parameters
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      real(DP), intent(in) :: dscale
      integer, intent(in) :: NVAR

      ! input/output parameters
      real(DP), dimension(NVAR,*), intent(inout) :: Da,Db
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      integer :: i,ii,ij,j,ji

      ! Loop over all rows
      do i = 1, size(Kld)-1

        ! Diagonal entry
        ii = Kdiagonal(i)
        Db(:,ii) = 0._DP

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J, that is, the upper
        ! right triangular part of the matrix is traversed
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          Db(:,ij) = dscale * (Da(:,ij) - Da(:,ji))
          Db(:,ji) = Db(:,ij)
        end do
      end do

    end subroutine calcAsymmPartMat9Dble

    !**************************************************************
    ! Compute symmetric part for CRS matrix stored in matrix format 9

    subroutine calcAsymmPartMat9Sngl(Kld, Kcol, Kdiagonal, Ksep, NVAR, fscale, Fa, Fb)

      ! input parameters
      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      real(SP), intent(in) :: fscale
      integer, intent(in) :: NVAR

      ! input/output parameters
      real(SP), dimension(NVAR,*), intent(inout) :: Fa,Fb
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      integer :: i,ii,ij,j,ji

      ! Loop over all rows
      do i = 1, size(Kld)-1

        ! Diagonal entry
        ii = Kdiagonal(i)
        Fb(:,ii) = 0._SP

        ! Loop over all off-diagonal matrix entries IJ which are
        ! adjacent to node J such that I < J, that is, the upper
        ! right triangular part of the matrix is traversed
        do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get node number J, the corresponding matrix positions JI,
          ! and let the separator point to the next entry
          j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

          Fb(:,ij) = fscale * (Fa(:,ij) - Fa(:,ji))
          Fb(:,ji) = Fb(:,ij)
        end do
      end do

    end subroutine calcAsymmPartMat9Sngl

  end subroutine lsyssc_createMatrixAsymmPart

  !****************************************************************************

!<function>

  function lsyssc_calcDeterminant(rmatrix, Icol, Irow) result(ddeterminant)

!<description>
    ! This subroutine calculates the determinant of a (sub-)matrix. If
    ! the optional arrays Irow and/or Icol are given, then the
    ! determinant of the submatrix is calculated.
!</description>

!</function>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! OPTIONAL: column indices to which submatrix is restricted
    integer, dimension(:), intent(in), optional :: Icol

    ! OPTIONAL: row indices to which submatrix is restricted
    integer, dimension(:), intent(in), optional :: Irow
!</input>

!<result>
    ! Value of the determinant
    real(DP) :: ddeterminant
!</result>

    ! local variables
    real(DP), dimension(:), pointer :: p_Da,p_DaSub
    real(SP), dimension(:), pointer :: p_Fa,p_FaSub
    integer, dimension(:), pointer :: p_Kld,p_Kcol
    integer, dimension(:), pointer :: p_IrowIdx,p_IcolIdx
    integer :: h_IrowIdx,h_IcolIdx,NEQ,NCOLS

    ! Check if matrix is square matrix
    if (rmatrix%NEQ .ne. rmatrix%NCOLS) then
      call output_line('Matrix is not square matrix!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_calcDeterminant')
      call sys_halt()
    end if

    ! Initialisation
    h_IcolIdx = ST_NOHANDLE; h_IrowIdx = ST_NOHANDLE

    ! What type of matrix are we?
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX1)
      !-------------------------------------------------------------------------
      ! Compute determinant of full matrix
      !-------------------------------------------------------------------------

      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        nullify(p_DaSub)

        ! Collect submatrix?
        if (present(Icol)) then
          if (present(Irow)) then
            ! Both row and column restrictions
            call collectSubmatrixMat1Dble(rmatrix%NEQ, rmatrix%NCOLS, p_Da,&
                                          Icol, Irow, NEQ, NCOLS, p_DaSub)
            ddeterminant = calcDetMat1Dble(NEQ, NCOLS, p_DaSub)
            deallocate(p_DaSub)
          else
            ! Only column restrictions
            call storage_new('lsyssc_calcDeterminant', 'IrowIdx', rmatrix%NEQ,&
                ST_INT, h_IrowIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IrowIdx, p_IrowIdx)
            call collectSubmatrixMat1Dble(rmatrix%NEQ, rmatrix%NCOLS, p_Da,&
                                          Icol, p_IrowIdx, NEQ, NCOLS, p_DaSub)
            ddeterminant = calcDetMat1Dble(NEQ, NCOLS, p_DaSub)
            deallocate(p_DaSub)
            call storage_free(h_IrowIdx)
          end if
        else
          if (present(Irow)) then
            ! Only row restrictions
            call storage_new('lsyssc_calcDeterminant', 'IcolIdx', rmatrix%NEQ,&
                ST_INT, h_IcolIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IcolIdx, p_IcolIdx)
            call collectSubmatrixMat1Dble(rmatrix%NEQ, rmatrix%NCOLS, p_Da,&
                                          p_IcolIdx, Irow, NEQ, NCOLS, p_DaSub)
            ddeterminant = calcDetMat1Dble(NEQ, NCOLS, p_DaSub)
            deallocate(p_DaSub)
            call storage_free(h_IcolIdx)
          else
            ! No restrictions at all
            ddeterminant = calcDetMat1Dble(rmatrix%NEQ, rmatrix%NCOLS, p_Da)
          end if
        end if

      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        nullify(p_FaSub)

        ! Collect submatrix?
        if (present(Icol)) then
          if (present(Irow)) then
            ! Both row and column restrictions
            call collectSubmatrixMat1Sngl(rmatrix%NEQ, rmatrix%NCOLS, p_Fa,&
                                          Icol, Irow, NEQ, NCOLS, p_FaSub)
            ddeterminant = calcDetMat1Sngl(NEQ, NCOLS, p_FaSub)
            deallocate(p_FaSub)
          else
            ! Only column restrictions
            call storage_new('lsyssc_calcDeterminant', 'IrowIdx', rmatrix%NEQ,&
                ST_INT, h_IrowIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IrowIdx, p_IrowIdx)
            call collectSubmatrixMat1Sngl(rmatrix%NEQ, rmatrix%NCOLS, p_Fa,&
                                          Icol, p_IrowIdx, NEQ, NCOLS, p_FaSub)
            ddeterminant = calcDetMat1Sngl(NEQ, NCOLS, p_FaSub)
            deallocate(p_FaSub)
            call storage_free(h_IrowIdx)
          end if
        else
          if (present(Irow)) then
            ! Only row restrictions
            call storage_new('lsyssc_calcDeterminant', 'IcolIdx', rmatrix%NEQ,&
                ST_INT, h_IcolIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IcolIdx, p_IcolIdx)
            call collectSubmatrixMat1Sngl(rmatrix%NEQ, rmatrix%NCOLS, p_Fa,&
                                          p_IcolIdx, Irow, NEQ, NCOLS, p_FaSub)
            ddeterminant = calcDetMat1Sngl(NEQ, NCOLS, p_FaSub)
            deallocate(p_FaSub)
            call storage_free(h_IcolIdx)
          else
            ! No restrictions at all
            ddeterminant = calcDetMat1Sngl(rmatrix%NEQ, rmatrix%NCOLS, p_Fa)
          end if
        end if

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_calcDeterminant')
        call sys_halt()
      end select

    case (LSYSSC_MATRIXD)
      !-------------------------------------------------------------------------
      ! Compute determinant of diagonal matrix
      !-------------------------------------------------------------------------

      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)

        ! Collect submatrix?
        if (present(Icol)) then
          if (present(Irow)) then
            ! Both row and column restrictions
            ddeterminant = calcDetMatDDble(p_Da, Icol, Irow)
          else
            ! Only column restrictions
            call storage_new('lsyssc_calcDeterminant', 'IrowIdx', rmatrix%NEQ,&
                ST_INT, h_IrowIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IrowIdx, p_IrowIdx)
            ddeterminant = calcDetMatDDble(p_Da, Icol, p_IrowIdx)
            call storage_free(h_IrowIdx)
          end if
        else
          if (present(Irow)) then
            ! Only row restrictions
            call storage_new('lsyssc_calcDeterminant', 'IcolIdx', rmatrix%NEQ,&
                ST_INT, h_IcolIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IcolIdx, p_IcolIdx)
            ddeterminant = calcDetMatDDble(p_Da, p_IcolIdx, Irow)
            call storage_free(h_IcolIdx)
          else
            ! No restrictions at all
            ddeterminant = calcDetQuickMatDDble(p_Da)
          end if
        end if

      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)

        ! Collect submatrix?
        if (present(Icol)) then
          if (present(Irow)) then
            ! Both row and column restrictions
            ddeterminant = calcDetMatDSngl(p_Fa, Icol, Irow)
          else
            ! Only column restrictions
            call storage_new('lsyssc_calcDeterminant', 'IrowIdx', rmatrix%NEQ,&
                ST_INT, h_IrowIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IrowIdx, p_IrowIdx)
            ddeterminant = calcDetMatDSngl(p_Fa, Icol, p_IrowIdx)
            call storage_free(h_IrowIdx)
          end if
        else
          if (present(Irow)) then
            ! Only row restrictions
            call storage_new('lsyssc_calcDeterminant', 'IcolIdx', rmatrix%NEQ,&
                ST_INT, h_IcolIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IcolIdx, p_IcolIdx)
            ddeterminant = calcDetMatDSngl(p_Fa, p_IcolIdx, Irow)
            call storage_free(h_IcolIdx)
          else
            ! No restrictions at all
            ddeterminant = calcDetQuickMatDSngl(p_Fa)
          end if
        end if

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_calcDeterminant')
        call sys_halt()
      end select

    case (LSYSSC_MATRIX7,LSYSSC_MATRIX9)
      !-------------------------------------------------------------------------
      ! Compute determinant of sparse matrix
      !-------------------------------------------------------------------------

      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

      select case (rmatrix%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix, p_Da)
        nullify(p_DaSub)

        ! Collect submatrix?
        if (present(Icol)) then
          if (present(Irow)) then
            ! Both row and column restrictions
            call collectSubmatrixMat79Dble(p_Kld, p_Kcol, p_Da,&
                                           Icol, Irow, NEQ, NCOLS, p_DaSub)
            ddeterminant = calcDetMat1Dble(NEQ, NCOLS, p_DaSub)
            deallocate(p_DaSub)
          else
            ! Only column restrictions
            call storage_new('lsyssc_calcDeterminant', 'IrowIdx', rmatrix%NEQ,&
                ST_INT, h_IrowIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IrowIdx, p_IrowIdx)
            call collectSubmatrixMat79Dble(p_Kld, p_Kcol, p_Da,&
                                           Icol, p_IrowIdx, NEQ, NCOLS, p_DaSub)
            ddeterminant = calcDetMat1Dble(NEQ, NCOLS, p_DaSub)
            deallocate(p_DaSub)
            call storage_free(h_IrowIdx)
          end if
        else
          if (present(Irow)) then
            ! Only row restrictions
            call storage_new('lsyssc_calcDeterminant', 'IcolIdx', rmatrix%NEQ,&
                ST_INT, h_IcolIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IcolIdx, p_IcolIdx)
            call collectSubmatrixMat79Dble(p_Kld, p_Kcol, p_Da,&
                                           p_IcolIdx, Irow, NEQ, NCOLS, p_DaSub)
            ddeterminant = calcDetMat1Dble(NEQ, NCOLS, p_DaSub)
            deallocate(p_DaSub)
            call storage_free(h_IcolIdx)
          else
            ! No restrictions at all
            call storage_new('lsyssc_calcDeterminant', 'IcolIdx', rmatrix%NEQ,&
                ST_INT, h_IcolIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IcolIdx, p_IcolIdx)
            call collectSubmatrixMat79Dble(p_Kld, p_Kcol, p_Da,&
                                           p_IcolIdx, p_IcolIdx, NEQ, NCOLS, p_DaSub)
            ddeterminant = calcDetMat1Dble(rmatrix%NEQ, rmatrix%NCOLS, p_DaSub)
            deallocate(p_DaSub)
            call storage_free(h_IcolIdx)
          end if
        end if

      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix, p_Fa)
        nullify(p_FaSub)

        ! Collect submatrix?
        if (present(Icol)) then
          if (present(Irow)) then
            ! Both row and column restrictions
            call collectSubmatrixMat79Sngl(p_Kld, p_Kcol, p_Fa,&
                                           Icol, Irow, NEQ, NCOLS, p_FaSub)
            ddeterminant = calcDetMat1Sngl(NEQ, NCOLS, p_FaSub)
            deallocate(p_FaSub)
          else
            ! Only column restrictions
            call storage_new('lsyssc_calcDeterminant', 'IrowIdx', rmatrix%NEQ,&
                ST_INT, h_IrowIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IrowIdx, p_IrowIdx)
            call collectSubmatrixMat79Sngl(p_Kld, p_Kcol, p_Fa,&
                                           Icol, p_IrowIdx, NEQ, NCOLS, p_FaSub)
            ddeterminant = calcDetMat1Sngl(NEQ, NCOLS, p_FaSub)
            deallocate(p_FaSub)
            call storage_free(h_IrowIdx)
          end if
        else
          if (present(Irow)) then
            ! Only row restrictions
            call storage_new('lsyssc_calcDeterminant', 'IcolIdx', rmatrix%NEQ,&
                ST_INT, h_IcolIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IcolIdx, p_IcolIdx)
            call collectSubmatrixMat79Sngl(p_Kld, p_Kcol, p_Fa,&
                                           p_IcolIdx, Irow, NEQ, NCOLS, p_FaSub)
            ddeterminant = calcDetMat1Sngl(NEQ, NCOLS, p_FaSub)
            deallocate(p_FaSub)
            call storage_free(h_IcolIdx)
          else
            ! No restrictions at all
            call storage_new('lsyssc_calcDeterminant', 'IcolIdx', rmatrix%NEQ,&
                ST_INT, h_IcolIdx, ST_NEWBLOCK_ORDERED)
            call storage_getbase_int(h_IcolIdx, p_IcolIdx)
            call collectSubmatrixMat79Sngl(p_Kld, p_Kcol, p_Fa,&
                                           p_IcolIdx, p_IcolIdx, NEQ, NCOLS, p_FaSub)
            ddeterminant = calcDetMat1Sngl(rmatrix%NEQ, rmatrix%NCOLS, p_FaSub)
            deallocate(p_FaSub)
            call storage_free(h_IcolIdx)
          end if
        end if

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_calcDeterminant')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_WARNING,OU_MODE_STD,'lsyssc_calcDeterminant')
      call sys_halt()
    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Collect submatrix from full matrix with double-valued data

    subroutine collectSubmatrixMat1Dble(NEQ, NCOLS, Da, IcolIdx, IrowIdx,&
                                        NEQSub, NCOLSSub, p_DaSub)

      real(DP), dimension(:), intent(in) :: Da
      integer, dimension(:), intent(in) :: IcolIdx,IrowIdx
      integer, intent(in) :: NEQ,NCOLS

      real(DP), dimension(:), pointer :: p_DaSub
      integer, intent(out) :: NEQSub,NCOLSSub

      ! local variables
      integer :: idx1,idx2,icol,irow

      ! Initialisation
      NEQSub = 0
      NCOLSSub = 0

      ! Count number of rows
      do idx1 = 1, size(IrowIdx)
        if (IrowIdx(idx1) .ne. 0) NEQSub = NEQSub+1
      end do

      ! Count number of columns
      do idx2 = 1, size(IcolIdx)
        if (IcolIdx(idx2) .ne. 0) NCOLSSub = NCOLSSub+1
      end do

      ! Allocate storage for submatrix (if required)
      if (associated(p_DaSub) .and. size(p_DaSub) .lt. NEQSub*NCOLSSub) then
        deallocate(p_DaSub)
        allocate(p_DaSub(NEQSub*NCOLSSub))
      else
        allocate(p_DaSub(NEQSub*NCOLSSub))
      end if

      ! Initialise row counter
      irow = 0

      ! Loop over all rows
      row: do idx1 = 1, size(IrowIdx)

        ! Skip row?
        if (IrowIdx(idx1) .eq. 0) cycle row

        ! Increase row counter
        irow = irow+1

        ! Initialise column counter
        icol = 0

        ! Loop over all columns
        col: do idx2 = 1, size(IcolIdx)

          ! Skip column?
          if (IcolIdx(idx2) .eq. 0) cycle col

          ! Increase column counter
          icol = icol+1

          ! Copy entry A(i,j) to A_sub(i_loc,j_loc)
          p_DaSub(NCOLSSub*(irow-1)+icol) = Da(NCOLS*(IrowIdx(idx1)-1)+IcolIdx(idx2))
        end do col
      end do row

    end subroutine collectSubmatrixMat1Dble

    !**************************************************************
    ! Collect submatrix from full matrix with single-valued data

    subroutine collectSubmatrixMat1Sngl(NEQ, NCOLS, Fa, IcolIdx, IrowIdx,&
                                        NEQSub, NCOLSSub, p_FaSub)

      real(SP), dimension(:), intent(in) :: Fa
      integer, dimension(:), intent(in) :: IcolIdx,IrowIdx
      integer, intent(in) :: NEQ,NCOLS

      real(SP), dimension(:), pointer :: p_FaSub
      integer, intent(out) :: NEQSub,NCOLSSub

      ! local variables
      integer :: idx1,idx2,icol,irow

      ! Initialisation
      NEQSub = 0
      NCOLSSub = 0

      ! Count number of rows
      do idx1 = 1, size(IrowIdx)
        if (IrowIdx(idx1) .ne. 0) NEQSub = NEQSub+1
      end do

      ! Count number of columns
      do idx2 = 1, size(IcolIdx)
        if (IcolIdx(idx2) .ne. 0) NCOLSSub = NCOLSSub+1
      end do

      ! Allocate storage for submatrix (if required)
      if (associated(p_FaSub) .and. size(p_FaSub) .lt. NEQSub*NCOLSSub) then
        deallocate(p_FaSub)
        allocate(p_FaSub(NEQSub*NCOLSSub))
      else
        allocate(p_FaSub(NEQSub*NCOLSSub))
      end if

      ! Initialise row counter
      irow = 0

      ! Loop over all rows
      row: do idx1 = 1, size(IrowIdx)

        ! Skip row?
        if (IrowIdx(idx1) .eq. 0) cycle row

        ! Increase row counter
        irow = irow+1

        ! Initialise column counter
        icol = 0

        ! Loop over all columns
        col: do idx2 = 1, size(IcolIdx)

          ! Skip column?
          if (IcolIdx(idx2) .eq. 0) cycle col

          ! Increase column counter
          icol = icol+1

          ! Copy entry A(i,j) to A_sub(i_loc,j_loc)
          p_FaSub(NCOLSSub*(irow-1)+icol) = Fa(NCOLS*(IrowIdx(idx1)-1)+IcolIdx(idx2))
        end do col
      end do row

    end subroutine collectSubmatrixMat1Sngl

    !**************************************************************
    ! Collect submatrix from CSR matrix with double-valued data

    subroutine collectSubmatrixMat79Dble(Kld, Kcol, Da, IcolIdx, IrowIdx,&
                                         NEQSub, NCOLSSub, p_DaSub)

      real(DP), dimension(:), intent(in) :: Da
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, dimension(:), intent(in) :: IcolIdx,IrowIdx

      real(DP), dimension(:), pointer :: p_DaSub
      integer, intent(out) :: NEQSub,NCOLSSub

      ! local variables
      integer :: idx1,idx2,icol,irow,ia

      ! Initialisation
      NEQSub = 0
      NCOLSSub = 0

      ! Count number of rows
      do idx1 = 1, size(IrowIdx)
        if (IrowIdx(idx1) .ne. 0) NEQSub = NEQSub+1
      end do

      ! Count number of columns
      do idx2 = 1, size(IcolIdx)
        if (IcolIdx(idx2) .ne. 0) NCOLSSub = NCOLSSub+1
      end do

      ! Allocate storage for submatrix (if required)
      if (associated(p_DaSub) .and. size(p_DaSub) .lt. NEQSub*NCOLSSub) then
        deallocate(p_DaSub)
        allocate(p_DaSub(NEQSub*NCOLSSub))
      else
        allocate(p_DaSub(NEQSub*NCOLSSub))
      end if

      ! Initialise row counter
      irow = 0

      ! Loop over all rows
      row: do idx1 = 1, size(IrowIdx)

        ! Skip row?
        if (IrowIdx(idx1) .eq. 0) cycle row

        ! Increase row counter
        irow = irow+1

        ! Initialise column counter
        icol = 0

        ! Loop over all columns
        col: do idx2 = 1, size(IcolIdx)

          ! Skip column?
          if (IcolIdx(idx2) .eq. 0) cycle col

          ! Increase column counter
          icol = icol+1

          ! Look for entry A(i,j) in sparse matrix
          do ia = Kld(IrowIdx(idx1)),Kld(IrowIdx(idx1)+1)-1
            if (Kcol(ia) .eq. IcolIdx(idx2)) then
              p_DaSub(NEQSub*(icol-1)+irow) = Da(ia)
              cycle col
            end if
          end do

          ! Fill with zeris
          p_DaSub(NEQSub*(icol-1)+irow) = 0.0_DP
        end do col
      end do row

    end subroutine collectSubmatrixMat79Dble

    !**************************************************************
    ! Collect submatrix from CSR matrix with single-valued data

    subroutine collectSubmatrixMat79Sngl(Kld, Kcol, Fa, IcolIdx, IrowIdx,&
                                         NEQSub, NCOLSSub, p_FaSub)

      real(SP), dimension(:), intent(in) :: Fa
      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, dimension(:), intent(in) :: IcolIdx,IrowIdx

      real(SP), dimension(:), pointer :: p_FaSub
      integer, intent(out) :: NEQSub,NCOLSSub

      ! local variables
      integer :: idx1,idx2,icol,irow,ia

      ! Initialisation
      NEQSub = 0
      NCOLSSub = 0

      ! Count number of rows
      do idx1 = 1, size(IrowIdx)
        if (IrowIdx(idx1) .ne. 0) NEQSub = NEQSub+1
      end do

      ! Count number of columns
      do idx2 = 1, size(IcolIdx)
        if (IcolIdx(idx2) .ne. 0) NCOLSSub = NCOLSSub+1
      end do

      ! Allocate storage for submatrix (if required)
      if (associated(p_FaSub) .and. size(p_FaSub) .lt. NEQSub*NCOLSSub) then
        deallocate(p_FaSub)
        allocate(p_FaSub(NEQSub*NCOLSSub))
      else
        allocate(p_FaSub(NEQSub*NCOLSSub))
      end if

      ! Initialise row counter
      irow = 0

      ! Loop over all rows
      row: do idx1 = 1, size(IrowIdx)

        ! Skip row?
        if (IrowIdx(idx1) .eq. 0) cycle row

        ! Increase row counter
        irow = irow+1

        ! Initialise column counter
        icol = 0

        ! Loop over all columns
        col: do idx2 = 1, size(IcolIdx)

          ! Skip column?
          if (IcolIdx(idx2) .eq. 0) cycle col

          ! Increase column counter
          icol = icol+1

          ! Look for entry A(i,j) in sparse matrix
          do ia = Kld(IrowIdx(idx1)),Kld(IrowIdx(idx1)+1)-1
            if (Kcol(ia) .eq. IcolIdx(idx2)) then
              p_FaSub(NEQSub*(icol-1)+irow) = Fa(ia)
              cycle col
            end if
          end do

          ! Fill with zeris
          p_FaSub(NEQSub*(icol-1)+irow) = 0.0_SP
        end do col
      end do row

    end subroutine collectSubmatrixMat79Sngl

    !**************************************************************
    ! Compute determinant of full matrix with double-valued data

    function calcDetMat1Dble(NEQ, NCOLS, Da) result(ddeterminant)

      integer, intent(in) :: NEQ,NCOLS

      real(DP), dimension(NEQ,NCOLS), intent(inout) :: Da

      real(DP) :: ddeterminant

      ! local variable
      integer, dimension(:), allocatable :: Iperm
      integer :: iaux,icol,ieq,ipermCount,ipiv,jeq
      real(DP) :: daux

      ! Allocate and initialize permutation vector
      allocate(Iperm(NEQ)); ipermCount = 0
      do ieq = 1,NEQ
        Iperm(ieq) = ieq
      end do

      ! Perform Gaussian elimination with partial pivoting, that is,
      ! swap rows so that the diagonal entry with the larges modules
      ! is used as pivoting element in each step

      ! Initialise determinant
      ddeterminant = 1.0_DP

      ! Perform Gaussian elimination
      row: do ieq = 1, NEQ

        ! Row-wise pivoting
        ipiv = ieq

        do jeq = ieq,NEQ
          if (abs(Da(Iperm(jeq),ieq)) .gt. abs(Da(Iperm(ipiv),ieq))) then
            ipiv = jeq
          end if
        end do

        ! Swap rows IPIV <-> IEQ
        if (Iperm(ipiv) .ne. Iperm(ieq)) then
          iaux        = Iperm(ieq)
          Iperm(ieq)  = Iperm(ipiv)
          Iperm(ipiv) = iaux
          ipermCount  = ipermCount+1
        end if

        ! Get pivot element
        daux = Da(Iperm(ieq),ieq)

        ! Check for zero determinant
        if (abs(daux) .le. SYS_EPSREAL_DP) then
          ! Return zero determinant
          ddeterminant = 0.0_DP
          return
        end if

        ! Elimination step
        col: do jeq = ieq+1,NEQ
          Da(Iperm(jeq),ieq) = Da(Iperm(jeq),ieq) / daux
          do icol = ieq+1,NCOLS
            Da(Iperm(jeq),icol) = Da(Iperm(jeq),icol)&
                                - Da(Iperm(ieq),icol)*Da(Iperm(jeq),ieq)
          end do
        end do col

        ! Update determinant
        ddeterminant = ddeterminant*daux
      end do row

      ! Adjust sign of determinant
      if (mod(ipermCount,2) .ne. 0)&
          ddeterminant = -ddeterminant

      ! Deallocate temporal memory
      deallocate(Iperm)

    end function calcDetMat1Dble

    !**************************************************************
    ! Compute determinant of full matrix with single-valued data

    function calcDetMat1Sngl(NEQ, NCOLS, Fa) result(ddeterminant)

      integer, intent(in) :: NEQ,NCOLS

      real(SP), dimension(NEQ,NCOLS), intent(inout) :: Fa

      real(DP) :: ddeterminant

      ! local variable
      integer, dimension(:), allocatable :: Iperm
      integer :: iaux,icol,ieq,ipermCount,ipiv,jeq
      real(SP) :: faux

      ! Allocate and initialize permutation vector
      allocate(Iperm(NEQ)); ipermCount = 0
      do ieq = 1,NEQ
        Iperm(ieq) = ieq
      end do

      ! Perform Gaussian elimination with partial pivoting, that is,
      ! swap rows so that the diagonal entry with the larges modules
      ! is used as pivoting element in each step

      ! Initialise determinant
      ddeterminant = 1.0_DP

      ! Perform Gaussian elimination
      row: do ieq = 1, NEQ

        ! Row-wise pivoting
        ipiv = ieq

        do jeq = ieq,NEQ
          if (abs(Fa(Iperm(jeq),ieq)) .gt. abs(Fa(Iperm(ipiv),ieq))) then
            ipiv = jeq
          end if
        end do

        ! Swap rows IPIV <-> IEQ
        if (ipiv .ne. ieq) then
          iaux        = Iperm(ieq)
          Iperm(ieq)  = Iperm(ipiv)
          Iperm(ipiv) = iaux
          ipermCount  = ipermCount+1
        end if

        ! Get pivot element
        faux = Fa(Iperm(ieq),ieq)

        ! Check for zero determinant
        if (abs(faux) .le. SYS_EPSREAL_SP) then
          ! Return zero determinant
          ddeterminant = 0.0_DP
          return
        end if

        ! Elimination
        col: do jeq = ieq+1,NEQ
          Fa(Iperm(jeq),ieq) = Fa(Iperm(jeq),ieq) / faux
          do icol = ieq+1,NCOLS
            Fa(Iperm(jeq),icol) = Fa(Iperm(jeq),icol)&
                                - Fa(Iperm(ieq),icol)*Fa(Iperm(jeq),ieq)
          end do
        end do col

        ! Update determinant
        ddeterminant = ddeterminant*real(faux,DP)
      end do row

      ! Adjust sign of determinant
      if (mod(ipermCount,2) .ne. 0)&
          ddeterminant = -ddeterminant

      ! Deallocate temporal memory
      deallocate(Iperm)

    end function calcDetMat1Sngl

    !**************************************************************
    ! Compute determinant of diagonal matrix with double-valued data
    ! without any restriction on the rows and/or columns

    function calcDetQuickMatDDble(Da) result(ddeterminant)

      real(DP), dimension(:), intent(in) :: Da

      real(DP) :: ddeterminant

      ! local variable
      integer :: idx

      ! Initialisation
      ddeterminant = 1.0_DP

      ! Compute product of diagonal entries
      do idx = 1,size(Da)
        ddeterminant = ddeterminant*Da(idx)
      end do

    end function calcDetQuickMatDDble

    !**************************************************************
    ! Compute determinant of diagonal matrix with single-valued data
    ! without any restriction on the rows and/or columns

    function calcDetQuickMatDSngl(Fa) result(ddeterminant)

      real(SP), dimension(:), intent(in) :: Fa

      real(DP) :: ddeterminant

      ! local variable
      integer :: idx

      ! Initialisation
      ddeterminant = 1.0_DP

      ! Compute product of diagonal entries
      do idx = 1,size(Fa)
        ddeterminant = ddeterminant*real(Fa(idx),DP)
      end do

    end function calcDetQuickMatDSngl

    !**************************************************************
    ! Compute determinant of diagonal matrix with double-valued data

    function calcDetMatDDble(Da, IcolIdx, IrowIdx) result(ddeterminant)

      real(DP), dimension(:), intent(in) :: Da
      integer, dimension(:), intent(in) :: IcolIdx,IrowIdx

      real(DP) :: ddeterminant

      ! local variable
      integer :: idx1,idx2,ncol,nrow

      ! Initialisation
      idx1 = 1; idx2 = 1
      nrow = 0; ncol = 0
      ddeterminant = 1.0_DP

      ! Compute product of diagonal entries
      do while ((idx1 .le. size(IrowIdx)) .and.&
                (idx2 .le. size(IcolIdx)))

        ! Skip row?
        if (IrowIdx(idx1) .eq. 0) then
          idx1 = idx1+1; cycle
        end if

        ! Skip column?
        if (IcolIdx(idx2) .eq. 0) then
          idx2 = idx2+1; cycle
        end if

        ! Update number of processed rows and columns
        nrow = nrow + 1
        ncol = ncol + 1

        ! Are we at a diagonal entry of the original matrix?
        if (IrowIdx(idx1) .eq. IcolIdx(idx2)) then
          ! Update determinant
          ddeterminant = ddeterminant*Da(IrowIdx(idx1))
        else
          ! Return zero determinant
          ddeterminant = 0.0_DP
          return
        end if

        ! Increase row and column index
        idx1 = idx1 + 1
        idx2 = idx2 + 1
      end do

      ! Check if matrix is square matrix
      if (nrow .ne. ncol) then
        call output_line('Submatrix is not square matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_calcDeterminant')
        call sys_halt()
      end if

    end function calcDetMatDDble

    !**************************************************************
    ! Compute determinant of diagonal matrix with single-valued data

    function calcDetMatDSngl(Fa, IcolIdx, IrowIdx) result(ddeterminant)

      real(SP), dimension(:), intent(in) :: Fa
      integer, dimension(:), intent(in) :: IcolIdx,IrowIdx

      real(DP) :: ddeterminant

      ! local variable
      integer :: idx1,idx2,ncol,nrow

      ! Initialisation
      idx1 = 1; idx2 = 1
      nrow = 0; ncol = 0
      ddeterminant = 1.0_DP

      ! Compute product of diagonal entries
      do while ((idx1 .le. size(IrowIdx)) .and.&
                (idx2 .le. size(IcolIdx)))

        ! Skip row?
        if (IrowIdx(idx1) .eq. 0) then
          idx1 = idx1+1; cycle
        end if

        ! Skip column?
        if (IcolIdx(idx2) .eq. 0) then
          idx2 = idx2+1; cycle
        end if

        ! Update number of processed rows and columns
        nrow = nrow + 1
        ncol = ncol + 1

        ! Are we at a diagonal entry of the original matrix?
        if (IrowIdx(idx1) .eq. IcolIdx(idx2)) then

          ! Update determinant
          ddeterminant = ddeterminant*real(Fa(IrowIdx(idx1)),DP)
        else
          ! Return zero determinant
          ddeterminant = 0.0_DP
          return
        end if

        ! Increase row and column index
        idx1 = idx1 + 1
        idx2 = idx2 + 1
      end do

      ! Check if matrix is square matrix
      if (nrow .ne. ncol) then
        call output_line('Submatrix is not square matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_calcDeterminant')
        call sys_halt()
      end if

    end function calcDetMatDSngl

  end function lsyssc_calcDeterminant

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_rebuildKdiagonal (Kcol, Kld, Kdiagonal, neq)

!<description>
  ! Internal auxiliary routine.
  ! Rebuilds the Kdiagonal-array of a structure-9 matrix.
!</description>

!<input>
  ! Row structure in the matrix
  integer, dimension(:), intent(in) :: Kld

  ! Column structure of the matrix
  integer, dimension(:), intent(inout) :: Kcol

  ! Dimension of the matrix
  integer, intent(in) :: neq
!</input>

!<output>
  ! Pointers to the diagonal entries of the matrix.
  integer, dimension(:), intent(out) :: Kdiagonal
!</output>

!</subroutine>

  ! local variables
  integer(I32) :: i, j

  ! loop through each row
  do i = 1, neq

    ! Loop through each column in this row.
    ! Search for the first element on or above the diagonal.
    do j = Kld(i), Kld(i+1)-1

      ! Check if we reached the position of the diagonal entry...
      if (Kcol(j) .ge. i) exit

    end do

    ! Save the position of the diagonal entry
    Kdiagonal(i) = j

  end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortCSRdouble (Kcol, Kld, Kdiagonal, neq, Da, nintl)

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
  integer, dimension(:), intent(in) :: Kld

  ! Dimension of the matrix
  integer, intent(in) :: neq

  ! Dimension of the interleaved submatrices (if any)
  integer, intent(in), optional :: nintl
!</input>

!<inputoutput>
  ! OPTIONAL:
  ! On input:  the matrix entries to be resorted,
  ! On output: the resorted matrix entries
  real(DP), dimension(:), intent(inout), optional :: Da

  ! On input:  the column numbers to be resorted,
  ! On output: the resorted column numbers
  integer, dimension(:), intent(inout) :: Kcol
!</inputoutput>

!<output>
  ! Pointers to the diagonal entries of the matrix.
  integer, dimension(:), intent(out) :: Kdiagonal
!</output>

!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: Daux
  real(DP) :: aux
  integer :: h_Daux
  integer(I32) :: i,j

  if (present(Da)) then

    ! Check if size of interleave matrix is specified. If this is the
    ! case then each "move" is performed for a local vector
    if (present(nintl)) then

      ! Allocate memory for auxiliary vector
      call storage_new ('lsyssc_sortCSRdouble', 'Daux', nintl, &
          ST_DOUBLE, h_Daux,ST_NEWBLOCK_NOINIT)
      call storage_getbase_double(h_Daux,Daux)

      ! loop through each row
      do i = 1, neq

        ! Take the diagonal element
        Daux = Da(nintl*(Kld(i)-1)+1:nintl*Kld(i))

        ! Loop through each column in this row.
        ! Shift every entry until the diagonal is reached.
        do j = Kld(i)+1, Kld(i+1)-1

          ! Check if we reached the position of the diagonal entry...
          if (Kcol(j)>i) exit

          Kcol(j-1) = Kcol(j)
          Da(nintl*(j-2)+1:nintl*(j-1)) = Da(nintl*(j-1)+1:nintl*j)

        end do

        ! If we have reached the diagonal, we can stop and save our
        ! diagonal entry from the first position there. The rest of the
        ! line is in ascending order according to the specifications of
        ! storage technique 7.

        Kcol(j-1) = i
        Da(nintl*(j-2)+1:nintl*(j-1)) = Daux

        ! Save the position of the diagonal entry
        Kdiagonal(i) = j-1

      end do

      ! Free auxiliary memory
      call storage_free(h_Daux)

    else

      ! loop through each row
      do i = 1, neq

        ! Take the diagonal element
        aux = Da(Kld(i))

        ! Loop through each column in this row.
        ! Shift every entry until the diagonal is reached.
        do j = Kld(i)+1, Kld(i+1)-1

          ! Check if we reached the position of the diagonal entry...
          if (Kcol(j)>i) exit

          Kcol(j-1) = Kcol(j)
          Da(j-1) = Da(j)

        end do

        ! If we have reached the diagonal, we can stop and save our
        ! diagonal entry from the first position there. The rest of the
        ! line is in ascending order according to the specifications of
        ! storage technique 7.

        Kcol(j-1) = i
        Da(j-1) = aux

        ! Save the position of the diagonal entry
        Kdiagonal(i) = j-1

      end do

    end if

  else

    ! loop through each row
    do i = 1, neq

      ! Loop through each column in this row.
      ! Shift every entry until the diagonal is reached.
      do j = Kld(i)+1, Kld(i+1)-1

        ! Check if we reached the position of the diagonal entry...
        if (Kcol(j)>i) exit

        Kcol(j-1) = Kcol(j)

      end do

      ! If we have reached the diagonal, we can stop and save our
      ! diagonal entry from the first position there. The rest of the
      ! line is in ascending order according to the specifications of
      ! storage technique 7.

      Kcol(j-1) = i

      ! Save the position of the diagonal entry
      Kdiagonal(i) = j-1

    end do

  end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortCSRsingle (Kcol, Kld, Kdiagonal, neq, Fa, nintl)

!<description>
  ! Internal auxiliary routine.
  ! Sorts the entries in each row of the matrix in ascending order.
  ! The input matrix is assumed to be in storage technique 7 with the
  ! first element in each row to be the diagonal element!
  ! Creates the Kdiagonal array for a structure-9 matrix.
  !
  ! Single precision version
!</description>

!<input>
  ! Row pointer in the matrix
  integer, dimension(:), intent(in) :: Kld

  ! Dimension of the matrix
  integer, intent(in) :: neq

  ! Dimension of the interleaved submatrices (if any)
  integer, intent(in), optional :: nintl
!</input>

!<inputoutput>
  ! OPTIONAL:
  ! On input:  the matrix entries to be resorted,
  ! On output: the resorted matrix entries
  real(SP), dimension(:), intent(inout), optional :: Fa

  ! On input:  the column numbers to be resorted,
  ! On output: the resorted column numbers
  integer, dimension(:), intent(inout) :: Kcol
!</inputoutput>

!<output>
  ! Pointers to the diagonal entries of the matrix.
  integer, dimension(:), intent(out) :: Kdiagonal
!</output>

!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: Faux
  real(DP) :: aux
  integer :: h_Faux
  integer(I32) :: i,j

  if (present(Fa)) then

    ! Check if size of interleave matrix is specified. If this is the
    ! case then each "move" is performed for a local vector
    if (present(nintl)) then

      ! Allocate memory for auxiliary vector
      call storage_new ('lsyssc_sortCSRsingle', 'Faux', nintl, &
          ST_DOUBLE, h_Faux,ST_NEWBLOCK_NOINIT)
      call storage_getbase_double(h_Faux,Faux)

      ! loop through each row
      do i = 1, neq

        ! Take the diagonal element
        Faux = Fa(nintl*(Kld(i)-1)+1:nintl*Kld(i))

        ! Loop through each column in this row.
        ! Shift every entry until the diagonal is reached.
        do j = Kld(i)+1, Kld(i+1)-1

          ! Check if we reached the position of the diagonal entry...
          if (Kcol(j)>i) exit

          Kcol(j-1) = Kcol(j)
          Fa(nintl*(j-2)+1:nintl*(j-1)) = Fa(nintl*(j-1)+1:nintl*j)

        end do

        ! If we have reached the diagonal, we can stop and save our
        ! diagonal entry from the first position there. The rest of the
        ! line is in ascending order according to the specifications of
        ! storage technique 7.

        Kcol(j-1) = i
        Fa(nintl*(j-2)+1:nintl*(j-1)) = Faux

        ! Save the position of the diagonal entry
        Kdiagonal(i) = j-1

      end do

      ! Free auxiliary memory
      call storage_free(h_Faux)

    else

      ! loop through each row
      do i = 1, neq

        ! Take the diagonal element
        aux = Fa(Kld(i))

        ! Loop through each column in this row.
        ! Shift every entry until the diagonal is reached.
        do j = Kld(i)+1, Kld(i+1)-1

          ! Check if we reached the position of the diagonal entry...
          if (Kcol(j)>i) exit

          Kcol(j-1) = Kcol(j)
          Fa(j-1) = Fa(j)

        end do

        ! If we have reached the diagonal, we can stop and save our
        ! diagonal entry from the first position there. The rest of the
        ! line is in ascending order according to the specifications of
        ! storage technique 7.

        Kcol(j-1) = i
        Fa(j-1) = aux

        ! Save the position of the diagonal entry
        Kdiagonal(i) = j-1

      end do

    end if

  else

    ! loop through each row
    do i = 1, neq

      ! Loop through each column in this row.
      ! Shift every entry until the diagonal is reached.
      do j = Kld(i)+1, Kld(i+1)-1

        ! Check if we reached the position of the diagonal entry...
        if (Kcol(j)>i) exit

        Kcol(j-1) = Kcol(j)

      end do

      ! If we have reached the diagonal, we can stop and save our
      ! diagonal entry from the first position there. The rest of the
      ! line is in ascending order according to the specifications of
      ! storage technique 7.

      Kcol(j-1) = i

      ! Save the position of the diagonal entry
      Kdiagonal(i) = j-1

    end do

  end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_unsortCSRdouble (Kcol, Kld, Kdiagonal, neq, Da, nintl)

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
  integer, dimension(:), intent(in) :: Kld

  ! Dimension of the matrix
  integer, intent(in) :: neq

  ! Pointers to the diagonal entries of the matrix.
  integer, dimension(:), intent(in) :: Kdiagonal

  ! Dimension of the interleave submatrices (if any)
  integer, intent(in), optional :: nintl
!</input>

!<inputoutput>
  ! OPTIONAL:
  ! On input:  the matrix entries to be resorted,
  ! On output: the resorted matrix entries
  real(DP), dimension(:), intent(inout), optional :: Da

  ! On input:  the column numbers to be resorted,
  ! On output: the resorted column numbers
  integer, dimension(:), intent(inout) :: Kcol
!</inputoutput>

!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: Daux
  real(DP) :: aux
  integer :: h_Daux
  integer(I32) :: i, j

  if (present(Da)) then

    ! Check if size of interleave matrix is specified. If this is the
    ! case then each "move" is performed for a local vector
    if (present(nintl)) then

      ! Allocate memory for auxiliary vector
      call storage_new ('lsyssc_sortCSRdouble', 'Daux', nintl, &
          ST_DOUBLE, h_Daux,ST_NEWBLOCK_NOINIT)
      call storage_getbase_double(h_Daux,Daux)

      ! loop through each row
      do i = 1, neq

        ! Take the diagonal element
        Daux = Da(nintl*(Kdiagonal(i)-1)+1:nintl*Kdiagonal(i))

        ! Loop through each column in this row.
        ! Shift every entry one element to the right.
        do j = Kdiagonal(i),Kld(i)+1,-1

          Kcol(j) = Kcol(j-1)
          Da(nintl*(j-1)+1:nintl*j) = Da(nintl*(j-2)+1:nintl*(j-1))

        end do

        ! Put the diagonal to the front.
        Kcol(Kld(i)) = i
        Da(nintl*(Kld(i)-1)+1:nintl*Kld(i)) = Daux

      end do

    else

      ! loop through each row
      do i = 1, neq

        ! Take the diagonal element
        aux = Da(Kdiagonal(i))

        ! Loop through each column in this row.
        ! Shift every entry one element to the right.
        do j = Kdiagonal(i),Kld(i)+1,-1

          Kcol(j) = Kcol(j-1)
          Da(j) = Da(j-1)

        end do

        ! Put the diagonal to the front.
        Kcol(Kld(i)) = i
        Da(Kld(i)) = aux

      end do

    end if

  else

    ! loop through each row
    do i = 1, neq

      ! Loop through each column in this row.
      ! Shift every entry one element to the right.
      do j = Kdiagonal(i),Kld(i)+1,-1
        Kcol(j) = Kcol(j-1)
      end do

      ! Put the diagonal to the front.
      Kcol(Kld(i)) = i

    end do

  end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_unsortCSRsingle (Kcol, Kld, Kdiagonal, neq, Fa, nintl)

!<description>
  ! Internal auxiliary routine.
  ! Unorts the entries in each row of the matrix in ascending order.
  ! This searches in each row of a matrix for the diagonal element
  ! and shifts it to the front.
  !
  ! Single precision version
!</description>

!<input>
  ! Row pointer in the matrix
  integer, dimension(:), intent(in) :: Kld

  ! Dimension of the matrix
  integer, intent(in) :: neq

  ! Pointers to the diagonal entries of the matrix.
  integer, dimension(:), intent(in) :: Kdiagonal

  ! Dimension of the interleave submatrices (if any)
  integer, intent(in), optional :: nintl
!</input>

!<inputoutput>
  ! OPTIONAL:
  ! On input:  the matrix entries to be resorted,
  ! On output: the resorted matrix entries
  real(SP), dimension(:), intent(inout), optional :: Fa

  ! On input:  the column numbers to be resorted,
  ! On output: the resorted column numbers
  integer, dimension(:), intent(inout) :: Kcol
!</inputoutput>

!</subroutine>

  ! local variables
  real(SP), dimension(:), pointer :: Faux
  real(SP) :: aux
  integer :: h_Faux
  integer(I32) :: i, j

  if (present(Fa)) then

    ! Check if size of interleave matrix is specified. If this is the
    ! case then each "move" is performed for a local vector
    if (present(nintl)) then

      ! Allocate memory for auxiliary vector
      call storage_new ('lsyssc_sortCSRsingle', 'Faux', nintl, &
          ST_SINGLE, h_Faux,ST_NEWBLOCK_NOINIT)
      call storage_getbase_single(h_Faux,Faux)

      ! loop through each row
      do i = 1, neq

        ! Take the diagonal element
        Faux = Fa(nintl*(Kdiagonal(i)-1)+1:nintl*Kdiagonal(i))

        ! Loop through each column in this row.
        ! Shift every entry one element to the right.
        do j = Kdiagonal(i),Kld(i)+1,-1

          Kcol(j) = Kcol(j-1)
          Fa(nintl*(j-1)+1:nintl*j) = Fa(nintl*(j-2)+1:nintl*(j-1))

        end do

        ! Put the diagonal to the front.
        Kcol(Kld(i)) = i
        Fa(nintl*(Kld(i)-1)+1:nintl*Kld(i)) = Faux

      end do

    else

      ! loop through each row
      do i = 1, neq

        ! Take the diagonal element
        aux = Fa(Kdiagonal(i))

        ! Loop through each column in this row.
        ! Shift every entry one element to the right.
        do j = Kdiagonal(i),Kld(i)+1,-1

          Kcol(j) = Kcol(j-1)
          Fa(j) = Fa(j-1)

        end do

        ! Put the diagonal to the front.
        Kcol(Kld(i)) = i
        Fa(Kld(i)) = aux

      end do

    end if

  else

    ! loop through each row
    do i = 1, neq

      ! Loop through each column in this row.
      ! Shift every entry one element to the right.
      do j = Kdiagonal(i),Kld(i)+1,-1
        Kcol(j) = Kcol(j-1)
      end do

      ! Put the diagonal to the front.
      Kcol(Kld(i)) = i

    end do

  end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_synchroniseSortVecVec (rvectorSrc,rvectorDst,rtemp,bautoUnsort)

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
  type(t_vectorScalar), intent(in) :: rvectorSrc

  ! OPTIONAL: Whether or not to check the target vector if it is already
  ! sorted. Default is TRUE.
  ! If set to TRUE and the target vector has already a sorting structure
  ! attached, the target is sorted back.
  ! If set to FALSE, the sorting of the target vector is deactivated and
  ! rsortStrategy is installed as sorting strategy. The data of the vector
  ! gets invalid. (Usually used for temp vectors.)
  logical, intent(in), optional :: bautoUnsort
!</input>

!<inputoutput>
  ! Destination vector; is resorted according to the sort strategy in rvectorSrc
  ! or is unsorted, if rvectorSrc is unsorted.
  ! Must have the same size as rvectorSrc.
  type(t_vectorScalar), intent(inout) :: rvectorDst

  ! OPTIONAL: A temporary vector. Must be of the same data type as rvector.
  ! Must be at least as large as rvectorDst.
  type(t_vectorScalar), intent(inout) :: rtemp
!</inputoutput>

!</subroutine>

    if (rvectorSrc%NEQ .ne. rvectorDst%NEQ) then
      call output_line("Vectors have different size!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_synchroniseSortVecVec")
      call sys_halt()
    end if

    ! If both are unsorted or both sorted in the same way, there is nothing to do.
    if ((.not. rvectorSrc%bisSorted) .and. (.not. rvectorDst%bisSorted)) return
    
    if ((rvectorSrc%bisSorted) .and. (rvectorDst%bisSorted)) then
      ! Exit if sorting strategy equal, nothing to do.
      if (rvectorSrc%p_rsortStrategy%ctype .eq. rvectorDst%p_rsortStrategy%ctype) return
    end if

    ! Sort the target vector back if sorted.
    if (.not. present(bautoUnsort)) then
      call lsyssc_sortVector (rvectorDst,.false.,rtemp)
    else if (bautoUnsort) then
      call lsyssc_sortVector (rvectorDst,.false.,rtemp)
    end if
    
    ! Synchronise the strategy
    call lsyssc_setSortStrategy (rvectorDst,rvectorSrc%p_rsortStrategy,bautoUnsort)
    
    ! Probably activate the sorting
    call lsyssc_sortVector (rvectorDst,rvectorSrc%bisSorted,rtemp)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_synchroniseSortMatVec (rmatrixSrc,rvectorDst,rtemp,bautoUnsort)

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
  type(t_matrixScalar), intent(in) :: rmatrixSrc

  ! OPTIONAL: Whether or not to check the target vector if it is already
  ! sorted. Default is TRUE.
  ! If set to TRUE and the target vector has already a sorting structure
  ! attached, the target is sorted back.
  ! If set to FALSE, the sorting of the target vector is deactivated and
  ! rsortStrategy is installed as sorting strategy. The data of the vector
  ! gets invalid. (Usually used for temp vectors.)
  logical, intent(in), optional :: bautoUnsort
!</input>

!<inputoutput>
  ! Destination vector; is resorted according to the sort strategy in rmatrixSrc
  ! or is unsorted, if rmatrixSrc is unsorted.
  ! Must have the same size (NEQ) as rmatrixSrc.
  type(t_vectorScalar), intent(inout) :: rvectorDst

  ! OPTIOAL: A temporary vector. Must be of the same data type as rvector.
  ! Must be at least as large as rvectorDst.
  type(t_vectorScalar), intent(inout), optional :: rtemp
!</inputoutput>

!</subroutine>

    if (rmatrixSrc%NEQ .ne. rvectorDst%NEQ) then
      call output_line("Matrix and vector have different size!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_synchroniseSortMatVec")
      call sys_halt()
    end if

    if (rtemp%NEQ .lt. rvectorDst%NEQ) then
      call output_line("Auxiliary vector too small!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_synchroniseSortMatVec")
      call sys_halt()
    end if


    ! If both are unsorted or both sorted in the same way, there is nothing to do.
    if ((.not. rmatrixSrc%bcolumnsSorted) .and. (.not. rvectorDst%bisSorted)) return
    
    if ((rmatrixSrc%bcolumnsSorted) .and. (rvectorDst%bisSorted)) then
      ! Exit if sorting strategy equal, nothing to do.
      if (rmatrixSrc%p_rsortStrategyColumns%ctype .eq. rvectorDst%p_rsortStrategy%ctype) return
    end if
    
    ! Sort the target vector back if sorted.
    if (.not. present(bautoUnsort)) then
      call lsyssc_sortVector (rvectorDst,.false.,rtemp)
    else if (bautoUnsort) then
      call lsyssc_sortVector (rvectorDst,.false.,rtemp)
    end if
    
    ! Synchronise the strategy
    call lsyssc_setSortStrategy (rvectorDst,rmatrixSrc%p_rsortStrategyColumns,bautoUnsort)
    
    ! Probably activate the sorting
    call lsyssc_sortVector (rvectorDst,rmatrixSrc%bcolumnsSorted,rtemp)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_setSortStrategyVec (rvector,rsortStrategy,bautoUnsort)

!<description>
  ! Attaches a sorting strategy to a vector.
  ! the sorting can be activated at aty time by lsyssc_sortVector.
!</description>

!<input>
  ! Sorting strategy.
  type(t_sortStrategy), intent(in), target :: rsortStrategy
  
  ! OPTIONAL: Whether or not to check the target vector if it is already
  ! sorted. Default is TRUE.
  ! If set to TRUE and the target vector has already a sorting structure
  ! attached, the target is sorted back.
  ! If set to FALSE, the sorting of the target vector is deactivated and
  ! rsortStrategy is installed as sorting strategy. The data of the vector
  ! gets invalid. (Usually used for temp vectors.)
  logical, intent(in), optional :: bautoUnsort
!</input>

!<inputoutput>
  ! Vector to attach the sorting to
  type(t_vectorScalar), intent(inout) :: rvector
!</inputoutput>

!</subroutine>

    ! Revert any old sorting strategy.
    if (.not. present(bautoUnsort)) then
      
      if (rvector%bisSorted) then
        call lsyssc_sortVector(rvector,.false.)
      end if
      
    else if (bautoUnsort) then
      
      if (rvector%bisSorted) then
        call lsyssc_sortVector(rvector,.false.)
      end if
    
    end if

    rvector%bisSorted = .false.
    rvector%p_rsortStrategy => rsortStrategy

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_setSortStrategyMat (&
      rmatrix,rsortStrategyColumns,rsortStrategyRows,bautoUnsort)

!<description>
  ! Attaches a sorting strategy to a matrix.
  ! the sorting can be activated at aty time by lsyssc_sortMatrix.
!</description>

!<input>
  ! OPTIONAL: Sorting strategy for the columns.
  ! If not specified, the columns remain unsorted.
  type(t_sortStrategy), intent(in), target, optional :: rsortStrategyColumns

  ! OPTIONAL: Sorting strategy for the rows.
  ! If not specified, the rows remain unsorted.
  type(t_sortStrategy), intent(in), target, optional :: rsortStrategyRows

  ! OPTIONAL: Whether or not to check the target vector if it is already
  ! sorted. Default is TRUE.
  ! If set to TRUE and the target vector has already a sorting structure
  ! attached, the target is sorted back.
  ! If set to FALSE, the sorting of the target vector is deactivated and
  ! rsortStrategy is installed as sorting strategy. The data of the vector
  ! gets invalid. (Usually used for temp vectors.)
  logical, intent(in), optional :: bautoUnsort
!</input>

!<inputoutput>
  ! Matrix to attach the sorting to
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! Revert any old sorting strategy.
    if (.not. present(bautoUnsort)) then

      if (rmatrix%bcolumnsSorted .or. rmatrix%browsSorted) then
        call lsyssc_sortMatrix(rmatrix,.false.)
      end if

    else if (bautoUnsort) then

      if (rmatrix%bcolumnsSorted .or. rmatrix%browsSorted) then
        call lsyssc_sortMatrix(rmatrix,.false.)
      end if

    end if

    rmatrix%bcolumnsSorted = .false.
    rmatrix%browsSorted = .false.
    nullify(rmatrix%p_rsortStrategyColumns)
    nullify(rmatrix%p_rsortStrategyRows)

    if (present(rsortStrategyColumns)) then
      rmatrix%p_rsortStrategyColumns => rsortStrategyColumns
    end if
    
    if (present(rsortStrategyRows)) then
      rmatrix%p_rsortStrategyRows => rsortStrategyRows
    end if
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortVector (rvector,bsort,rtemp)

!<description>
  ! Resorts the entries of the given vector rvector or unsorts it.
  !
  ! The sorting strategy must have been attached to the vector in advance
  ! by lsyssc_setSortStrategy.
!</description>

!<input>
  ! TRUE if the sorting should be activated.
  ! FALSE if the sorting should be deactivated.
  logical, intent(in) :: bsort
!</input>

!<inputoutput>
  ! Vector to resort
  type(t_vectorScalar), intent(inout) :: rvector

  ! OPTIONAL: A temporary vector. Must be of the same data type as rvector.
  ! Must be at least as large as rvector.
  type(t_vectorScalar), intent(inout), target, optional :: rtemp
!</inputoutput>

!</subroutine>

  ! local variables
  integer, dimension(:), pointer :: p_Iperm
  real(DP), dimension(:), pointer :: p_Ddata,p_Ddata2
  real(SP), dimension(:), pointer :: p_Fdata,p_Fdata2
  type(t_vectorScalar), pointer :: p_rtemp
  type(t_vectorScalar), target :: rtempLocal
  integer :: NEQ
  
    ! nothing to do?
    if (.not. associated(rvector%p_rsortStrategy)) return
    if (bsort .eqv. rvector%bisSorted) return
    
    ! If the sorting strategy is "unsorted", there is nothing to do,
    ! we can leave the vector marked as "unsorted".
    if (bsort .and. rvector%p_rsortStrategy%ctype .eq. SSTRAT_UNSORTED) then
      return
    end if

    ! Temporary vector available? If not, create a new one based on rvector.
    if (present(rtemp)) then
      p_rtemp => rtemp
    else
      p_rtemp => rtempLocal
      
      ! Create the temp vector as a pure data array.
      call lsyssc_createVector (rtempLocal,rvector%NEQ)
    end if

    if (p_rtemp%NEQ .lt. rvector%NEQ) then
      call output_line("Temp vector too small!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_sortVector")
      call sys_halt()
    end if

    NEQ = rvector%NEQ
    
    ! Get pointers to the vector data.
    ! Copy the vector to the temp vector
    select case (rvector%cdataType)
    case (ST_DOUBLE)
      call lsyssc_getbase_double(rvector,p_Ddata)
      call lsyssc_getbase_double(p_rtemp,p_Ddata2)
      call lalg_copyVectorDble(p_Ddata,p_Ddata2,size(p_Ddata))
    case (ST_SINGLE)
      call lsyssc_getbase_single(rvector,p_Fdata)
      call lsyssc_getbase_single(p_rtemp,p_Fdata2)
      call lalg_copyVectorSngl(p_Fdata,p_Fdata2,size(p_Fdata))
    case default
      call output_line("Unsuppported data type!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_sortVector")
      call sys_halt()
    end select

    ! Sort the vector?
    if (bsort) then
      call sstrat_getUnsortedPosInfo(rvector%p_rsortStrategy,p_Iperm)
    else
      call sstrat_getSortedPosInfo(rvector%p_rsortStrategy,p_Iperm)
    end if

    select case (rvector%cdataType)
    case (ST_DOUBLE)
      ! Then do the sorting with the given permutation.
      call lalg_vectorSortDble (p_Ddata2,p_Ddata,p_Iperm)

    case (ST_SINGLE)
      ! Then do the sorting with the given permutation.
      call lalg_vectorSortSngl (p_Fdata2,p_Fdata,p_Iperm)

    case default
      call output_line("Unsuppported data type!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_sortVector")
      call sys_halt()

    end select

    ! Mark the vector as (un-)sorted
    rvector%bisSorted = bsort
    
    ! Remove the temp vector if it is ours.
    if (.not. present(rtemp)) then
      call lsyssc_releaseVector (rtempLocal)
    end if
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMatrix (rmatrix,bsort,bincludeEntries)

!<description>
  ! Matrix sorting. Sorts the entries and/or structure of a given
  ! matrix or unsorts them according to a permutation.
  !
  ! The sorting strategy must have been attached to the vector in advance
  ! by lsyssc_setSortStrategy.
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
  type(t_matrixScalar), intent(inout), target :: rmatrix
!</inputoutput>

!<input>
  ! Defines whether or not to sort the matrix.
  ! TRUE activates sorting, FALSE deactivates the sorting.
  logical, intent(in) :: bsort

  ! OPTIONAL: Defines whether or not to include the entries in the sorting.
  ! TRUE sorts structure and entries (default).
  ! FALSE only sorts the structure and ignores the entries.
  ! Note that in this case, the entries are left in a most likely undefined state.
  logical, intent(in), optional :: bincludeEntries
!</input>

!</subroutine>

  ! local variables
  integer, dimension(:), pointer :: p_Iperm,p_IpermInverse
  integer :: NEQ
  type(t_matrixScalar), target :: rmatrixLocal
  type(t_matrixScalar), pointer :: p_rmatrix
  integer :: cstructureDupFlag,ccontentDupFlag
  logical :: bwithEntries,bcols,brows

    ! nothing to do?
    if ((.not. associated(rmatrix%p_rsortStrategyColumns)) .and. &
        (.not. associated(rmatrix%p_rsortStrategyRows))) return
    if (bsort .eqv. (rmatrix%bcolumnsSorted .or. rmatrix%browsSorted)) return
    
    ! The optional parameters
    bwithEntries = .true.
    if (present(bincludeEntries)) bwithEntries = bincludeEntries

    ! Virtually transposed matrices cannot be sorted.
    if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
      call output_line("Sorting of virtually transposed matrices not supported.!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_sortMatrix")
      call sys_halt()
    end if

    ! It is a little bit tricky now if structure and/or content of rmatrix
    ! is shared with another matrix!
    ! If that is the case, we make a copy of our matrix and work with that.
    cstructureDupFlag = LSYSSC_DUP_SHARE
    ccontentDupFlag = LSYSSC_DUP_SHARE

    ! Check if the structure belongs to another matrix. If yes, we
    ! create a copy of the structure and apply the sorting to a copy
    ! of the matrix to prevent it from being destroyed.
    if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .eq. &
        LSYSSC_MSPEC_STRUCTUREISCOPY) then
      cstructureDupFlag = LSYSSC_DUP_COPY
    end if

    if (bwithEntries) then
    
      ! Check if the content belongs to another matrix. If yes, we
      ! create a copy of the structure and apply the sorting to a copy
      ! of the matrix to prevent it from being destroyed.
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .eq. &
          LSYSSC_MSPEC_CONTENTISCOPY) then
        ccontentDupFlag = LSYSSC_DUP_COPY
      end if
    
      if ((cstructureDupFlag .eq. LSYSSC_DUP_COPY) .and. &
          (ccontentDupFlag .eq. LSYSSC_DUP_COPY)) then
        ! The user wants to sort content and/or structure, but the everything belongs to
        ! another matrix! Sorting would destroy that matrix, so we do not allow
        ! that. Therefore, there is nothing to do here.
        return
      end if
      
    else
    
      if ((cstructureDupFlag .eq. LSYSSC_DUP_COPY) .and. &
          .not. lsyssc_hasMatrixContent(rmatrix)) then
        ! Also nothing to do.
        return
      end if
      
    end if
    
    ! Big question: Do we have to do a copy or not?
    ! The routine do_matsort always sorts the structure, but not
    ! necessarily the content...

    p_rmatrix => rmatrix
    if ((cstructureDupFlag .ne. LSYSSC_DUP_SHARE) .or. &
        (ccontentDupFlag .ne. LSYSSC_DUP_SHARE)) then
      ! Copy necessary.
      p_rmatrix => rmatrixLocal
      call lsyssc_duplicateMatrix (rmatrix,rmatrixLocal,cstructureDupFlag,ccontentDupFlag)
    end if
    
    ! Some preparations...
    NEQ = rmatrix%NEQ
    
    ! What to sort?
    brows = .false.
    bcols = .false.
    if (associated(p_rmatrix%p_rsortStrategyRows)) then
      brows = p_rmatrix%p_rsortStrategyRows%ctype .ne. SSTRAT_UNSORTED
    end if
    if (associated(p_rmatrix%p_rsortStrategyColumns)) then
      bcols = p_rmatrix%p_rsortStrategyColumns%ctype .ne. SSTRAT_UNSORTED
    end if

    ! Sort the matrix?
    if (bsort) then
    
      ! Call the sorting routine
      if (bcols .and. brows) then
        call sstrat_getUnsortedPosInfo(p_rmatrix%p_rsortStrategyRows,p_Iperm)
        call sstrat_getSortedPosInfo(p_rmatrix%p_rsortStrategyColumns,p_IpermInverse)

        call do_matsort (p_rmatrix,p_Iperm,p_IpermInverse,bwithEntries)

      else if (bcols) then
        call sstrat_getSortedPosInfo(p_rmatrix%p_rsortStrategyColumns,p_IpermInverse)

        call do_matsort (p_rmatrix,Itr2=p_IpermInverse,bsortEntries=bwithEntries)

      else if (brows) then
        call sstrat_getUnsortedPosInfo(p_rmatrix%p_rsortStrategyRows,p_Iperm)

        call do_matsort (p_rmatrix,Itr1=p_Iperm,bsortEntries=bwithEntries)
      end if
    
    ! Sort back.
    else
    
      ! Exchange the roles of the first and second part of p_Iperm.
      ! This makes the permutation to the inverse permutation and vice versa.
      ! Call the resort-subroutine with that to do the actual sorting.
      if (bcols .and. brows) then
        call sstrat_getUnsortedPosInfo(p_rmatrix%p_rsortStrategyColumns,p_Iperm)
        call sstrat_getSortedPosInfo(p_rmatrix%p_rsortStrategyRows,p_IpermInverse)

        call do_matsort (p_rmatrix,p_IpermInverse,p_Iperm,bwithEntries)
      
      else if (bcols) then
        call sstrat_getUnsortedPosInfo(p_rmatrix%p_rsortStrategyColumns,p_Iperm)

        call do_matsort (p_rmatrix,Itr2=p_Iperm,bsortEntries=bwithEntries)
      
      else if (brows) then
        call sstrat_getSortedPosInfo(p_rmatrix%p_rsortStrategyRows,p_IpermInverse)

        call do_matsort (p_rmatrix,Itr1=p_IpermInverse,bsortEntries=bwithEntries)
      end if
          
    end if
    
    ! Mark the matrix as (un-)sorted
    rmatrix%bcolumnsSorted = .false.
    rmatrix%browsSorted = .false.
    if (associated(rmatrix%p_rsortStrategyColumns)) then
      rmatrix%bcolumnsSorted = bsort .and. (rmatrix%p_rsortStrategyColumns%ctype .ne. SSTRAT_UNSORTED)
    end if
    if (associated(rmatrix%p_rsortStrategyRows)) then
      rmatrix%browsSorted = bsort .and. (rmatrix%p_rsortStrategyRows%ctype .ne. SSTRAT_UNSORTED)
    end if
    
    if (.not. associated(p_rmatrix,rmatrix)) then
      ! Release the temp matrix
      call lsyssc_releaseMatrix (rmatrixLocal)
    end if
    
  contains

    !----------------------------------------------------------------
    ! Sort matrix or matrix entries.
    ! This calls the actual resorting routine, depending
    ! on the information tags in the matrix.
    ! The routine does not replace handles in rmatrix.

    subroutine do_matsort (rmatrix,Itr1,Itr2,bsortEntries)

    ! The matrix to be resorted
    type(t_matrixScalar), intent(inout) :: rmatrix

    ! OPTIONAL: The transformation to use for the rows.
    ! If not specified, the identity is used.
    integer, dimension(:), intent(in), optional :: Itr1

    ! OPTIONAL: The inverse transformation to use for the columns.
    ! If not specified, the identity is used.
    integer, dimension(:), intent(in), optional :: Itr2

    ! TRUE  = sort matrix structure + entries
    ! FALSE = sort only matrix structure
    logical, intent(in) :: bsortEntries

    ! local variables

    integer, dimension(:), pointer :: p_Kld,p_KldTmp,p_Kdiag
    integer, dimension(:), pointer :: p_Kcol,p_KcolTmp
    real(DP), dimension(:), pointer :: p_Ddata,p_DdataTmp
    real(SP), dimension(:), pointer :: p_Fdata,p_FdataTmp
    type(t_matrixScalar) :: rtempMatrix
    integer :: NEQ

      NEQ = rmatrix%NEQ

      ! Which matrix configuration do we have?
      select case (rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9)

        if (.not. bsortEntries) then

          ! Sort only the structure of the matrix, keep the entries
          ! unchanged. Do not do anything if the structure does not belong
          ! to the matrix!
          if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .eq. 0) then

            ! Duplicate the matrix before resorting.
            call lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                        LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)

            ! Get the structure of the original and the temporary matrix
            call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
            call lsyssc_getbase_Kld (rmatrix,p_Kld)
            call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiag)
            call lsyssc_getbase_Kcol (rtempMatrix,p_KcolTmp)
            call lsyssc_getbase_Kld (rtempMatrix,p_KldTmp)

            ! Sort
            call lsyssc_sortMat9Struc (p_Kcol, p_KcolTmp, p_Kld, p_KldTmp, &
                                       p_Kdiag, Itr1, Itr2, NEQ)

            ! Remove temp matrix
            call lsyssc_releaseMatrix(rtempMatrix)

          end if

        else

          ! Do not do anything if structure & entries do not belong to our matrix.
          if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY) .eq. LSYSSC_MSPEC_ISCOPY) &
            return

          ! If the structure does not belong to the matrix, only sort the!
          ! entries. Otherwise, sort structure + entries

          if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0) then

            ! Get the structure of the matrix
            call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
            call lsyssc_getbase_Kld (rmatrix,p_Kld)
            call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiag)

            ! Create a copy of the matrix entries
            call lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                        LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

            select case (rmatrix%cdataType)
            case (ST_DOUBLE)
              ! Double precision version
              call lsyssc_getbase_double (rmatrix,p_Ddata)
              call lsyssc_getbase_double (rtempMatrix,p_DdataTmp)
              call lsyssc_sortMat9Ent_double (p_Ddata,p_DdataTmp,p_Kcol, &
                                              p_Kld, Itr1, Itr2, NEQ)
            case (ST_SINGLE)
              ! Single precision version
              call lsyssc_getbase_single (rmatrix,p_Fdata)
              call lsyssc_getbase_single (rtempMatrix,p_FdataTmp)
              call lsyssc_sortMat9Ent_single (p_Fdata,p_FdataTmp,p_Kcol, &
                                              p_Kld, Itr1, Itr2, NEQ)
            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_sortVector")
              call sys_halt()
            end select

            ! Remove temp matrix
            call lsyssc_releaseMatrix(rtempMatrix)

          else

            ! Duplicate the matrix before resorting.
            ! We need either a copy only of the structure or of the full matrix.
            call lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                         LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)

            ! Get the structure of the original and the temporary matrix
            call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
            call lsyssc_getbase_Kld (rmatrix,p_Kld)
            call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiag)
            call lsyssc_getbase_Kcol (rtempMatrix,p_KcolTmp)
            call lsyssc_getbase_Kld (rtempMatrix,p_KldTmp)

            select case (rmatrix%cdataType)
            case (ST_DOUBLE)
              ! Double precision version
              call lsyssc_getbase_double (rmatrix,p_Ddata)
              call lsyssc_getbase_double (rtempMatrix,p_DdataTmp)
              call lsyssc_sortMat9_double (p_Ddata,p_DdataTmp,p_Kcol, p_KcolTmp, &
                                           p_Kld, p_KldTmp, p_Kdiag, &
                                           Itr1, Itr2, NEQ)
            case (ST_SINGLE)
              ! Single precision version
              call lsyssc_getbase_single (rmatrix,p_Fdata)
              call lsyssc_getbase_single (rtempMatrix,p_FdataTmp)
              call lsyssc_sortMat9_single (p_Fdata,p_FdataTmp,p_Kcol, p_KcolTmp, &
                                           p_Kld, p_KldTmp, p_Kdiag, &
                                           Itr1, Itr2, NEQ)
            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_sortVector")
              call sys_halt()
            end select

            ! Remove temp matrix
            call lsyssc_releaseMatrix(rtempMatrix)

          end if

        end if

      case (LSYSSC_MATRIX7)

        if (.not. bsortEntries) then

          ! Duplicate the matrix before resorting.
          call lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                       LSYSSC_DUP_COPY,LSYSSC_DUP_SHARE)

          ! Get the structure of the original and the temporary matrix
          call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
          call lsyssc_getbase_Kld (rmatrix,p_Kld)
          call lsyssc_getbase_Kcol (rtempMatrix,p_KcolTmp)
          call lsyssc_getbase_Kld (rtempMatrix,p_KldTmp)

          ! Sort only the structure of the matrix, keep the entries
          ! unchanged.
          call lsyssc_sortMat7Struc (p_Kcol, p_KcolTmp, p_KldTmp, p_KldTmp, &
                                     Itr1, Itr2, NEQ)

          ! Remove temp matrix
          call lsyssc_releaseMatrix(rtempMatrix)

        else

          ! Do not do anything if structure & entries do not belong to our matrix.
          if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY) .eq. LSYSSC_MSPEC_ISCOPY) &
            return

          ! If the structure does not belong to the matrix, only sort the!
          ! entries. Otherwise, sort structure + entries

          if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0) then

            ! Get the structure of the matrix
            call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
            call lsyssc_getbase_Kld (rmatrix,p_Kld)

            ! Create a copy of the matrix entries
            call lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                        LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

            ! Sort only the entries
            select case (rmatrix%cdataType)
            case (ST_DOUBLE)
              ! Double precision version
              call lsyssc_getbase_double (rmatrix,p_Ddata)
              call lsyssc_getbase_double (rtempMatrix,p_DdataTmp)
              call lsyssc_sortMat7Ent_double (p_Ddata,p_DdataTmp,p_Kcol, &
                                              p_Kld, Itr1, Itr2, NEQ)
            case (ST_SINGLE)
              ! Single precision version
              call lsyssc_getbase_single (rmatrix,p_Fdata)
              call lsyssc_getbase_single (rtempMatrix,p_FdataTmp)
              call lsyssc_sortMat7Ent_single (p_Fdata,p_FdataTmp,p_Kcol, &
                                              p_Kld, Itr1, Itr2, NEQ)
            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_sortVector")
              call sys_halt()
            end select

            ! Remove temp matrix
            call lsyssc_releaseMatrix(rtempMatrix)

          else

            ! Duplicate the matrix before resorting.
            call lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                         LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)

            ! Get the structure of the original and the temporary matrix
            call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
            call lsyssc_getbase_Kld (rmatrix,p_Kld)
            call lsyssc_getbase_Kcol (rtempMatrix,p_KcolTmp)
            call lsyssc_getbase_Kld (rtempMatrix,p_KldTmp)

            ! Sort structure + entries
            select case (rmatrix%cdataType)
            case (ST_DOUBLE)
              ! Double precision version
              call lsyssc_getbase_double (rmatrix,p_Ddata)
              call lsyssc_getbase_double (rtempMatrix,p_DdataTmp)
              call lsyssc_sortMat7_double (p_Ddata,p_DdataTmp,p_Kcol, p_KcolTmp, &
                                          p_Kld, p_KldTmp, &
                                          Itr1, Itr2, NEQ)
            case (ST_SINGLE)
              ! Single precision version
              call lsyssc_getbase_single (rmatrix,p_Fdata)
              call lsyssc_getbase_single (rtempMatrix,p_FdataTmp)
              call lsyssc_sortMat7_single (p_Fdata,p_FdataTmp,p_Kcol, p_KcolTmp, &
                                          p_Kld, p_KldTmp, &
                                          Itr1, Itr2, NEQ)
            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_sortVector")
              call sys_halt()
            end select

            ! Remove temp matrix
            call lsyssc_releaseMatrix(rtempMatrix)

          end if

        end if

      case (LSYSSC_MATRIXD)

        ! D-matrices can only be sorted if we are allowed to sort the content --
        ! as they have no explicit structure!

        if (bsortEntries .and. &
            (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .eq. 0)) then

          ! Duplicate the matrix before resorting.
          call lsyssc_duplicateMatrix (rmatrix,rtempMatrix,&
                                       LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)

          ! Sort entries; there is no structure.
          select case (rmatrix%cdataType)
          case (ST_DOUBLE)
            ! Double precision version
            call lsyssc_getbase_double (rmatrix,p_Ddata)
            call lsyssc_getbase_double (rtempMatrix,p_DdataTmp)
            call lalg_vectorSortDble (p_DdataTmp, p_Ddata, Itr1)
          case (ST_SINGLE)
            ! Single precision version
            call lsyssc_getbase_single (rmatrix,p_Fdata)
            call lsyssc_getbase_single (rtempMatrix,p_FdataTmp)
            call lalg_vectorSortSngl (p_FdataTmp, p_Fdata, Itr1)
          case default
            call output_line("Unsupported data type!",&
                OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_sortVector")
            call sys_halt()
          end select

          ! Remove temp matrix
          call lsyssc_releaseMatrix(rtempMatrix)

        end if

      case default
        call output_line("Unsupported matrix format!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_sortVector")
        call sys_halt()

      end select

    end subroutine

  end subroutine

  !****************************************************************************
  ! ACTUAL SORTING ROUTINES
  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat7_double (Da, DaH, Icol, IcolH, &
                                     Ild, IldH, Itr1, Itr2, neq)

!<description>
    ! Resorts structure and entries of the given matrix, corresponding to Itr1/Itr2.
    ! Double precision version.
    !
    ! Storage technique 7 version.
!</description>

!<input>

    ! Number of equations
    integer, intent(in) :: neq

    ! Source matrix
    real(DP), dimension(:), intent(in) :: DaH

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting of the rows.
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse sorting of
    ! the columns.
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Destination matrix
    real(DP), dimension(:), intent(out) :: Da

    ! Column structure of destination matrix
    integer, dimension(:), intent(out) :: Icol

    ! Row positions of destination matrix
    integer, dimension(:), intent(out) :: Ild

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq

      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j=1, isize
        ! Get the matrix entry of the current column and write it to Da
        Da(ildIdx) = DaH(Ih2(j))

        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j)

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat7_single (Fa, FaH, Icol, IcolH, &
                                     Ild, IldH, Itr1, Itr2, neq)

!<description>
    ! Resorts structure and entries of the given matrix, corresponding to Itr1/Itr2.
    ! Single precision version.
    !
    ! Storage technique 7 version.
!</description>

!<input>

    ! Number of equations
    integer, intent(in) :: neq

    ! Source matrix
    real(SP), dimension(:), intent(in) :: FaH

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting of the rows.
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse sorting of
    ! the columns.
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Destination matrix
    real(SP), dimension(:), intent(out) :: Fa

    ! Column structure of destination matrix
    integer, dimension(:), intent(out) :: Icol

    ! Row positions of destination matrix
    integer, dimension(:), intent(out) :: Ild

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq
      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j=1, isize
        ! Get the matrix entry of the current column and write it to Da
        Fa(ildIdx) = FaH(Ih2(j))

        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j)

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat7_int (Ia, IaH, Icol, IcolH, &
                                     Ild, IldH, Itr1, Itr2, neq)

!<description>
    ! Resorts structure and entries of the given matrix, corresponding to Itr1/Itr2.
    ! Integer version.
    !
    ! Storage technique 7 version.
!</description>

!<input>

    ! Number of equations
    integer, intent(in) :: neq

    ! Source matrix
    integer, dimension(:), intent(in) :: IaH

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting of the rows.
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse sorting of
    ! the columns.
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Destination matrix
    integer, dimension(:), intent(out) :: Ia

    ! Column structure of destination matrix
    integer, dimension(:), intent(out) :: Icol

    ! Row positions of destination matrix
    integer, dimension(:), intent(out) :: Ild

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq
      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j=1, isize
        ! Get the matrix entry of the current column and write it to Da
        Ia(ildIdx) = IaH(Ih2(j))

        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j)

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat7Ent_double (Da, DaH, IcolH, &
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
    integer, intent(in) :: neq

    ! Source matrix
    real(DP), dimension(:), intent(in) :: DaH

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting of the rows.
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse sorting of
    ! the columns.
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Destination matrix
    real(DP), dimension(:), intent(out) :: Da

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq

      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j=1, isize
        ! Get the matrix entry of the current column and write it to Da
        Da(ildIdx) = DaH(Ih2(j))

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat7Ent_single (Fa, FaH, IcolH, &
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
    integer, intent(in) :: neq

    ! Source matrix
    real(SP), dimension(:), intent(in) :: FaH

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting of the rows.
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse sorting
    ! of the columns.
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Destination matrix
    real(SP), dimension(:), intent(out) :: Fa

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq
      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j=1, isize
        ! Get the matrix entry of the current column and write it to Da
        Fa(ildIdx) = FaH(Ih2(j))

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat7Ent_int (Ia, IaH, IcolH, &
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
    integer, intent(in) :: neq

    ! Source matrix
    integer, dimension(:), intent(in) :: IaH

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting of the rows
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse sorting
    ! of the columns.
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Destination matrix
    integer, dimension(:), intent(out) :: Ia

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq
      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j=1, isize
        ! Get the matrix entry of the current column and write it to Da
        Ia(ildIdx) = IaH(Ih2(j))

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat7Struc (Icol, IcolH, Ild, IldH, Itr1, Itr2, neq)

!<description>
    ! Resorts the structure of the given matrix, corresponding to Itr1/Itr2.
    ! Only the structure of a given matrix is resorted, the routine
    ! does not handle the entries.
    !
    ! Storage technique 7 version.
!</description>

!<input>

    ! Number of equations
    integer , intent(in) :: neq

    ! Permutation of 1..neq describing how to resort the rows.
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing how to sort the columns.
    ! Inverse column permutation.
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<inputoutput>

    ! Column description of matrix
    integer, dimension(:), intent(inout) :: Icol

    ! Row description of matrix
    integer, dimension(:), intent(inout) :: Ild

    ! Column structure of source matrix -> resorted matrix
    integer, dimension(:), intent(inout) :: IcolH

    ! Row positions of source matrix -> resorted matrix
    integer, dimension(:), intent(inout) :: IldH

!</inputoutput>
!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx
    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq
      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j=1, isize
        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j)

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx+1
      end do
    end do
    Ild(neq+1) = IldH(neq+1)

    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat9_double (Da, DaH, Icol, IcolH, &
                                     Ild, IldH, Idiag, Itr1, Itr2, neq)

!<description>
    ! Resorts structure and entries of the given matrix, corresponding to Itr1/Itr2.
    ! Double precision version.
    !
    ! Storage technique 9 version.
!</description>

!<input>

    ! Number of equations
    integer, intent(in) :: neq

    ! Source matrix
    real(DP), dimension(:), intent(in) :: DaH

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse permutation
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Destination matrix
    real(DP), dimension(:), intent(out) :: Da

    ! Column structure of destination matrix
    integer, dimension(:), intent(out) :: Icol

    ! Row positions of destination matrix
    integer, dimension(:), intent(out) :: Ild

    ! Positions of diagonal elements in the destination matrix
    integer, dimension(:), intent(out) :: Idiag

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq

      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      do idiagidx=1,isize
        if (Ih1(idiagidx) .ge. i) exit
      end do

      ! Write the position of the diagonal element to Idiag
      Idiag(i) = ildIdx + idiagidx-1

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j = 1,isize

        ! Get the matrix entry of the current column and write it to Da
        Da(ildIdx) = DaH(Ih2(j))

        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j)

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat9_single (Fa, FaH, Icol, IcolH, &
                                     Ild, IldH, Idiag, Itr1, Itr2, neq)

!<description>
    ! Resorts structure and entries of the given matrix, corresponding to Itr1/Itr2.
    ! Single precision version.
    !
    ! Storage technique 9 version.
!</description>

!<input>

    ! Number of equations
    integer, intent(in) :: neq

    ! Source matrix
    real(SP), dimension(:), intent(in) :: FaH

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse permutation
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Destination matrix
    real(SP), dimension(:), intent(out) :: Fa

    ! Column structure of destination matrix
    integer, dimension(:), intent(out) :: Icol

    ! Row positions of destination matrix
    integer, dimension(:), intent(out) :: Ild

    ! Positions of diagonal elements in the destination matrix
    integer, dimension(:), intent(out) :: Idiag

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq

      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      do idiagidx=1,isize
        if (Ih1(idiagidx) .ge. i) exit
      end do

      ! Write the position of the diagonal element to Idiag
      Idiag(i) = ildIdx + idiagidx-1

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j = 1,isize

        ! Get the matrix entry of the current column and write it to Da
        Fa(ildIdx) = FaH(Ih2(j))

        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j)

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat9_int (Ia, IaH, Icol, IcolH, &
                                     Ild, IldH, Idiag, Itr1, Itr2, neq)

!<description>
    ! Resorts structure and entries of the given matrix, corresponding to Itr1/Itr2.
    ! Double precision version.
    !
    ! Storage technique 9 version.
!</description>

!<input>

    ! Number of equations
    integer, intent(in) :: neq

    ! Source matrix
    integer, dimension(:), intent(in) :: IaH

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse permutation
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Destination matrix
    integer, dimension(:), intent(out) :: Ia

    ! Column structure of destination matrix
    integer, dimension(:), intent(out) :: Icol

    ! Row positions of destination matrix
    integer, dimension(:), intent(out) :: Ild

    ! Positions of diagonal elements in the destination matrix
    integer, dimension(:), intent(out) :: Idiag

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq

      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      do idiagidx=1,isize
        if (Ih1(idiagidx) .ge. i) exit
      end do

      ! Write the position of the diagonal element to Idiag
      Idiag(i) = ildIdx + idiagidx-1

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j = 1,isize

        ! Get the matrix entry of the current column and write it to Da
        Ia(ildIdx) = IaH(Ih2(j))

        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j)

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat9Ent_double (Da, DaH, IcolH, IldH, &
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
    integer, intent(in) :: neq

    ! Source matrix
    real(DP), dimension(:), intent(in) :: DaH

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting of the rows
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse permutation
    ! of the columns
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Destination matrix
    real(DP), dimension(:), intent(out) :: Da

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq

      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      do idiagidx=1,isize
        if (Ih1(idiagidx) .ge. i) exit
      end do

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j = 1,isize

        ! Get the matrix entry of the current column and write it to Da
        Da(ildIdx) = DaH(Ih2(j))

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat9Ent_single (Fa, FaH, IcolH, IldH, &
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
    integer, intent(in) :: neq

    ! Source matrix
    real(SP), dimension(:), intent(in) :: FaH

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting of the rows
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse permutation
    ! of the columns
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Destination matrix
    real(SP), dimension(:), intent(out) :: Fa

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq

      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      do idiagidx=1,isize
        if (Ih1(idiagidx) .ge. i) exit
      end do

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j = 1,isize

        ! Get the matrix entry of the current column and write it to Da
        Fa(ildIdx) = FaH(Ih2(j))

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat9Ent_int (Ia, IaH, IcolH, IldH, &
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
    integer, intent(in) :: neq

    ! Source matrix
    integer, dimension(:), intent(in) :: IaH

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting of the rows
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse permutation
    ! of the columns
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Destination matrix
    integer, dimension(:), intent(out) :: Ia

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq

      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      do idiagidx=1,isize
        if (Ih1(idiagidx) .ge. i) exit
      end do

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j = 1,isize

        ! Get the matrix entry of the current column and write it to Da
        Ia(ildIdx) = IaH(Ih2(j))

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortMat9Struc (Icol, IcolH, &
                                   Ild, IldH, Idiag, Itr1, Itr2, neq)

!<description>
    ! Resorts the structure of the given matrix, corresponding to Itr1/Itr2.
    ! Only the structure of a given matrix is resorted, the routine
    ! does not handle the entries.
    !
    ! Storage technique 9 version.
!</description>

!<input>

    ! Number of equations
    integer, intent(in) :: neq

    ! Column structure of source matrix
    integer, dimension(:), intent(in) :: IcolH

    ! Row positions of source matrix
    integer, dimension(:), intent(in) :: IldH

    ! Permutation of 1..neq describing the sorting of the rows
    integer, dimension(:), intent(in), optional :: Itr1

    ! Permutation of 1..neq describing the inverse permutation
    ! of the columns
    integer, dimension(:), intent(in), optional :: Itr2

!</input>

!<output>

    ! Column structure of destination matrix
    integer, dimension(:), intent(out) :: Icol

    ! Row positions of destination matrix
    integer, dimension(:), intent(out) :: Ild

    ! Positions of diagonal elements in the destination matrix
    integer, dimension(:), intent(out) :: Idiag

!</output>

!</subroutine>

    !local variables
    integer :: i, j, iidx, ildIdx, isize, ih1Idx, ih2Idx, idiagidx

    ! Temporary variables for saving data of each line
    integer, dimension(:), allocatable :: Ih1, Ih2

    ! Get memory for Ih1 and Ih2:
    allocate(Ih1(neq),Ih2(neq))

    ! We build the sorted matrix from the scratch. ildIdx points
    ! to the current 'output' position in the data array of
    ! the matrix.
    ildIdx = 1

    ! Loop through all lines. i is the 'destination' row, i.e.
    ! we fetch the row Itr1(i) and write it to row i.
    do i=1, neq

      ! Which row should be moved to here?
      iidx = i
      if (present(Itr1)) iidx = Itr1(i)

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
      if (.not. present(Itr2)) then
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = IcolH(j)
          Ih2(j-ih1Idx+1) = j
        end do
      else
        do j=ih1Idx, ih2Idx
          Ih1(j-ih1Idx+1) = Itr2(IcolH(j))
          Ih2(j-ih1Idx+1) = j
        end do
      end if

      ! Call a small sorting algorithm to sort both, Ih1 and Ih2, for Ih1.
      ! Afterwards:
      !  Ih1 = column numbers in the destination matrix, sorted
      !  Ih2 = positions of the entries in Ih1 in the source matrix
      call lsyssc_sortCR(Ih1(1:isize),Ih2(1:isize))

      ! Again loop through the row - this time to find the index of
      ! the first element in the upper triangular part of the matrix.
      do idiagidx=1,isize
        if (Ih1(idiagidx) .ge. i) exit
      end do

      ! Write the position of the diagonal element to Idiag
      Idiag(i) = ildIdx + idiagidx-1

      ! Now copy the source row to the current row.
      ! Loop through the columns of the current row:
      do j = 1,isize
        ! Get the column number and write it to Icol
        Icol(ildIdx) = Ih1(j)

        ! Increase the destination pointer in the matrix structure
        ildIdx = ildIdx + 1
      end do
    end do

    ! Close the matrix structure by writing out the last element of Kld.
    Ild(neq+1) = IldH(neq+1)

    ! Release temp variables
    deallocate(Ih1,Ih2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_sortCR (Ih1, Ih2)

!<description>
    ! Performs bubble sort for the vector Ih1 and does the
    ! same swaps also on vector Ih2.
!</description>

!<inputoutput>

    ! column vector
    integer, dimension(:), intent(inout) :: Ih1

    ! row vector; same size as Ih1
    integer, dimension(size(Ih1)), intent(inout) :: Ih2

!</inputoutput>
!</subroutine>

    !local variables
    integer :: iaux, icomp
    logical :: bmore

    bmore = .true.
    ! While there are swaps necessary...
    do while (bmore)
      bmore = .false.
      ! Test for all elements minus the last one
      do icomp=1, size(Ih1)-1
        ! If the order is correct, next entry; othewise swap entries
        if (Ih1(icomp) .gt. Ih1(icomp+1)) then
          iaux = Ih1(icomp)
          Ih1(icomp) = Ih1(icomp+1)
          Ih1(icomp+1) = iaux
          iaux = Ih2(icomp)
          Ih2(icomp) = Ih2(icomp+1)
          Ih2(icomp+1) = iaux
          bmore = .true.
        endif
      end do
    end do

  end subroutine

  !****************************************************************************
  ! SORTING ROUTINES FINISHED
  !****************************************************************************

!<subroutine>

  subroutine lsyssc_addIndex (h_Ix,ivalue,istartpos,ilength)

!<description>

  ! This is an auxiliary routine. It accepts a handle to an integer array
  ! and adds the value ivalue to each entry of this vector.
  ! This can be used e.g. to convert the index vector of a matrix
  ! from 1-based to 0-based for an external library (UMFPACK).

!</description>

!<input>
  ! Handle to an integer data array.
  integer, intent(in) :: h_Ix

  ! The value to add to every entry.
  integer :: ivalue

  ! OPTIONAL: Starting position of a part of the vector to modify; usually = 1
  integer, optional :: istartpos

  ! OPTIONAL: Length of the part of the vector to modify; <= SIZE(Ix)
  integer, optional :: ilength

!</input>

!</subroutine>

    ! Actual length
    integer :: iactlength,istart,i

    ! Vector
    integer, dimension(:), pointer :: p_Ix

    ! Get the array
    call storage_getbase_int (h_Ix,p_Ix)
    if (present(istartpos)) then
      istart = min(istartpos,size(p_Ix)+1)
    else
      istart = 1
    end if
    if (present(ilength)) then
      iactlength = min(ilength,size(p_Ix)-istart+1)
    else
      iactlength = size(p_Ix)-istart+1
    end if
    p_Ix => p_Ix(istart:istart + iactlength - 1)

    ! Increase by ivalue
    do i=1,iactlength
      p_Ix(i) = p_Ix(i) + ivalue
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_addConstant (rvector,dvalue)

!<description>
  ! Adds a constant to a vector.
!</description>

!<inputoutput>
  ! Vector to be modified.
  type(t_vectorScalar), intent(in) :: rvector
!</inputoutput>

!<input>
  ! Value to be added to the vector.
  real(DP), intent(in) :: dvalue
!</input>

!</subroutine>

    ! Vector
    real(DP), dimension(:), pointer :: p_Dx

    select case (rvector%cdataType)

    case (ST_DOUBLE)
      ! Get the data array
      call lsyssc_getbase_double (rvector,p_Dx)

      ! Increase by dvalue
      call lalg_vectorAddScalarDble (p_Dx,dvalue,rvector%NEQ)

    case default
      call output_line('Data type not implemented!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_addConstant')
      call sys_halt()

    end select

  end subroutine

  !****************************************************************************
!<function>

  real(DP) function lsyssc_vectorNorm (rx,cnorm,iposMax)

!<description>
  ! Calculates the norm of a vector. cnorm specifies the type of norm to
  ! calculate.
!</description>

!<input>
  ! Vector to calculate the norm of.
  type(t_vectorScalar), intent(in) :: rx

  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(in) :: cnorm
!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(out), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! < 0, if an error occurred (unknown norm).
!</result>

!<result>
  ! The scalar product (rx,ry) of the two block vectors.
!</result>

!</function>

  ! local variables
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata

  ! Is there data at all?
  if (rx%h_Ddata .eq. ST_NOHANDLE) then
    call output_line('Vector empty!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_vectorNorm')
    call sys_halt()
  end if

  ! Take care of the data type before doing a scalar product!
  select case (rx%cdataType)
  case (ST_DOUBLE)
    ! Get the array and calculate the norm
    call lsyssc_getbase_double (rx,p_Ddata)
    lsyssc_vectorNorm = lalg_norm(p_Ddata,cnorm,iposMax)

  case (ST_SINGLE)
    ! Get the array and calculate the norm
    call lsyssc_getbase_single (rx,p_Fdata)
    lsyssc_vectorNorm = lalg_norm(p_Fdata,cnorm,iposMax)

  case default
    call output_line('Unsupported data type!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_vectorNorm')
    call sys_halt()
  end select

  end function

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_invertedDiagMatVec (rmatrix,rvectorSrc,dscale,rvectorDst,&
                                        rperfconfig)

!<description>
  ! This routine multiplies the weighted inverted diagonal <tex>$domega*D^{-1}$</tex>
  ! of the matrix rmatrix with the vector rvectorSrc and stores the result
  ! into the vector rvectorDst:
  !   <tex>$$rvectorDst = dscale * D^{-1} * rvectorSrc$$ </tex>
  ! Both, rvectorSrc and rvectorDst may coincide.
!</description>

!<input>
  ! The matrix.
  type(t_matrixScalar), intent(in) :: rmatrix

  ! The source vector.
  type(t_vectorScalar), intent(in) :: rvectorSrc

  ! A multiplication factor. Standard value is 1.0_DP
  real(DP), intent(in) :: dscale

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The destination vector which receives the result.
  type(t_vectorScalar), intent(inout) :: rvectorDst
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: h_Idiag,ivar,NVAR
  integer :: i
  integer, dimension(:), pointer :: p_Kdiag
  real(DP), dimension(:), pointer :: p_Da, p_Dvec, p_Dvec2
  real(SP), dimension(:), pointer :: p_Fa, p_Fvec, p_Fvec2
  real(DP) :: dmyscale
  real(SP) :: fmyscale

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

  if (present(rperfconfig)) then
    p_rperfconfig => rperfconfig
  else
    p_rperfconfig => lsyssc_perfconfig
  end if

  ! Let us hope, the matrix and the vectors are compatible.
  ! As we multiply with D^{-1}, this forces the matrix to be quadratic
  ! and both vectors to be of the same length!
  call lsyssc_isVectorCompatible (rvectorSrc,rvectorDst)
  call lsyssc_isMatrixCompatible (rvectorSrc,rmatrix,.false.)

  if (rvectorSrc%cdataType .ne. rvectorDst%cdataType) then
    call output_line('Vectors have different precisions!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
    call sys_halt()
  end if

  ! Which matrix structure do we have?
  select case (rmatrix%cmatrixFormat)
  case (LSYSSC_MATRIX9,LSYSSC_MATRIX7)
    ! Get a pointer to the diagonal - either h_KLD or h_Kdiagonal
    if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX9) then
      h_Idiag = rmatrix%h_Kdiagonal
    else
      h_Idiag = rmatrix%h_Kld
    end if
    call storage_getbase_int(h_Idiag,p_Kdiag)

    ! Data type?
    select case (rmatrix%cdataType)
    case (ST_DOUBLE)
      call lsyssc_getbase_double (rmatrix,p_Da)
      dmyscale = dscale * rmatrix%dscaleFactor

      ! And the vector(s)?
      select case (rvectorSrc%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rvectorSrc,p_Dvec)
        call lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let us go...
        if (rvectorSrc%NVAR .eq. 1) then
          !$omp parallel do default(shared) if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            p_Dvec2(i) = p_Dvec(i)*dmyscale/p_Da(p_Kdiag(i))
          end do
          !$omp end parallel do
        else
          NVAR = rvectorSrc%NVAR
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Dvec2(NVAR*(i-1)+ivar) = p_Dvec(NVAR*(i-1)+ivar)&
                  *dmyscale/p_Da(p_Kdiag(i))
            end do
          end do
          !$omp end parallel do
        end if

      case (ST_SINGLE)
        call lsyssc_getbase_single (rvectorSrc,p_Fvec)
        call lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let us go...
        if (rvectorSrc%NVAR .eq. 1) then
          !$omp parallel do default(shared) if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            p_Fvec2(i) = p_Fvec(i)*dmyscale/p_Da(p_Kdiag(i))
          end do
          !$omp end parallel do
        else
          NVAR = rvectorSrc%NVAR
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Fvec2(NVAR*(i-1)+ivar) = p_Fvec(NVAR*(i-1)+ivar)&
                  *dmyscale/p_Da(p_Kdiag(i))
            end do
          end do
          !$omp end parallel do
        end if

      case default
        call output_line('Unsupported vector precision!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
        call sys_halt()
      end select

    case (ST_SINGLE)
      call lsyssc_getbase_single (rmatrix,p_Fa)
      fmyscale = real(dscale * rmatrix%dscaleFactor,SP)

      ! And the vector(s)?
      select case (rvectorSrc%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rvectorSrc,p_Dvec)
        call lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let us go...
        if (rvectorSrc%NVAR .eq. 1) then
          !$omp parallel do default(shared)
          do i=1,rvectorSrc%NEQ
            p_Dvec2(i) = p_Dvec(i)*fmyscale/p_Fa(p_Kdiag(i))
          end do
          !$omp end parallel do
        else
          NVAR = rvectorSrc%NVAR
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Dvec2(NVAR*(i-1)+ivar) = p_Dvec(NVAR*(i-1)+ivar)*&
                  fmyscale/p_Fa(p_Kdiag(i))
            end do
          end do
          !$omp end parallel do
        end if

      case (ST_SINGLE)
        call lsyssc_getbase_single (rvectorSrc,p_Fvec)
        call lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let us go...
        if (rvectorSrc%NVAR .eq. 1) then
          !$omp parallel do default(shared) if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            p_Fvec2(i) = p_Fvec(i)*fmyscale/p_Fa(p_Kdiag(i))
          end do
          !$omp end parallel do
        else
          NVAR = rvectorSrc%NVAR
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Fvec2(NVAR*(i-1)+ivar) = p_Fvec(NVAR*(i-1)+ivar)*&
                  fmyscale/p_Fa(p_Kdiag(i))
            end do
          end do
          !$omp end parallel do
        end if

      case default
        call output_line('Unsupported vector precision!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported matrix precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
      call sys_halt()
    end select

  case (LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7INTL)

    if (rmatrix%NVAR .ne. rvectorSrc%NVAR) then
      call output_line('Matrix/vectors have different interleaved format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
      call sys_halt()
    end if

    ! Number of interleaved variables
    NVAR = rmatrix%NVAR

    ! Get a pointer to the diagonal - either h_KLD or h_Kdiagonal
    if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX9INTL) then
      h_Idiag = rmatrix%h_Kdiagonal
    else
      h_Idiag = rmatrix%h_Kld
    end if
    call storage_getbase_int(h_Idiag,p_Kdiag)

    ! Data type?
    select case (rmatrix%cdataType)
    case (ST_DOUBLE)
      call lsyssc_getbase_double (rmatrix,p_Da)
      dmyscale = dscale * rmatrix%dscaleFactor

      ! And the vector(s)?
      select case (rvectorSrc%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rvectorSrc,p_Dvec)
        call lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let us go...
        select case (rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIX1)
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Dvec2(NVAR*(i-1)+ivar) = p_Dvec(NVAR*(i-1)+ivar)&
                  *dmyscale/p_Da(NVAR*NVAR*(p_Kdiag(i)-1)+NVAR*ivar)
            end do
          end do
          !$omp end parallel do

        case (LSYSSC_MATRIXD)
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Dvec2(NVAR*(i-1)+ivar) = p_Dvec(NVAR*(i-1)+ivar)&
                  *dmyscale/p_Da(NVAR*(p_Kdiag(i)-1)+ivar)
            end do
          end do
          !$omp end parallel do

        case default
          call output_line('Unsupported interleaved matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
          call sys_halt()
        end select

      case (ST_SINGLE)
        call lsyssc_getbase_single (rvectorSrc,p_Fvec)
        call lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let us go...
        select case (rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIX1)
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Fvec2(NVAR*(i-1)+ivar) = p_Fvec(NVAR*(i-1)+ivar)&
                  *dmyscale/p_Da(NVAR*NVAR*(p_Kdiag(i)-1)+NVAR*ivar)
            end do
          end do
          !$omp end parallel do

        case (LSYSSC_MATRIXD)
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Fvec2(NVAR*(i-1)+ivar) = p_Fvec(NVAR*(i-1)+ivar)&
                  *dmyscale/p_Da(NVAR*(p_Kdiag(i)-1)+ivar)
            end do
          end do
          !$omp end parallel do

        case default
          call output_line('Unsupported interleaved matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported vector precision!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
        call sys_halt()
      end select

    case (ST_SINGLE)
      call lsyssc_getbase_single (rmatrix,p_Fa)
      fmyscale = real(dscale * rmatrix%dscaleFactor,SP)

      ! And the vector(s)?
      select case (rvectorSrc%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rvectorSrc,p_Dvec)
        call lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let us go...
        select case (rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIX1)
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Dvec2(NVAR*(i-1)+ivar) = p_Dvec(NVAR*(i-1)+ivar)&
                  *fmyscale/p_Fa(NVAR*NVAR*(p_Kdiag(i)-1)+NVAR*ivar)
            end do
          end do
          !$omp end parallel do

        case (LSYSSC_MATRIXD)
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Dvec2(NVAR*(i-1)+ivar) = p_Dvec(NVAR*(i-1)+ivar)&
                  *fmyscale/p_Fa(NVAR*(p_Kdiag(i)-1)+ivar)
            end do
          end do
          !$omp end parallel do

        case default
          call output_line('Unsupported interleaved matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
          call sys_halt()
        end select

      case (ST_SINGLE)
        call lsyssc_getbase_single (rvectorSrc,p_Fvec)
        call lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let us go...
        select case (rmatrix%cinterleavematrixFormat)
        case (LSYSSC_MATRIX1)
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Fvec2(NVAR*(i-1)+ivar) = p_Fvec(NVAR*(i-1)+ivar)&
                  *fmyscale/p_Fa(NVAR*NVAR*(p_Kdiag(i)-1)+NVAR*ivar)
            end do
          end do
          !$omp end parallel do

        case (LSYSSC_MATRIXD)
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Fvec2(NVAR*(i-1)+ivar) = p_Fvec(NVAR*(i-1)+ivar)&
                  *fmyscale/p_Fa(NVAR*(p_Kdiag(i)-1)+ivar)
            end do
          end do
          !$omp end parallel do

        case default
          call output_line('Unsupported interleaved matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
          call sys_halt()
        end select

      case default
        call output_line('Unsupported vector precision!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported matrix precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
      call sys_halt()
    end select

  case (LSYSSC_MATRIXD)

    ! Data type?
    select case (rmatrix%cdataType)
    case (ST_DOUBLE)
      call lsyssc_getbase_double (rmatrix,p_Da)
      dmyscale = dscale * rmatrix%dscaleFactor

      ! And the vector(s)?
      select case (rvectorSrc%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rvectorSrc,p_Dvec)
        call lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let us go...
        if (rvectorSrc%NVAR .eq. 1) then
          !$omp parallel do default(shared) if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            p_Dvec2(i) = p_Dvec(i)*dmyscale/p_Da(i)
          end do
          !$omp end parallel do
        else
          NVAR = rvectorSrc%NVAR
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Dvec2(NVAR*(i-1)+ivar) = p_Dvec(NVAR*(i-1)+ivar)&
                  *dmyscale/p_Da(i)
            end do
          end do
          !$omp end parallel do
        end if

      case (ST_SINGLE)
        call lsyssc_getbase_single (rvectorSrc,p_Fvec)
        call lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let us go...
        if (rvectorSrc%NVAR .eq. 1) then
          !$omp parallel do default(shared) if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            p_Fvec2(i) = p_Fvec(i)*dmyscale/p_Da(i)
          end do
          !$omp end parallel do
        else
          NVAR = rvectorSrc%NVAR
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Fvec2(NVAR*(i-1)+ivar) = p_Fvec(NVAR*(i-1)+ivar)&
                  *dmyscale/p_Da(i)
            end do
          end do
          !$omp end parallel do
        end if

      case default
        call output_line('Unsupported vector precision!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
        call sys_halt()
      end select

    case (ST_SINGLE)
      call lsyssc_getbase_single (rmatrix,p_Fa)
      fmyscale = real(dscale * rmatrix%dscaleFactor,SP)

      ! And the vector(s)?
      select case (rvectorSrc%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double (rvectorSrc,p_Dvec)
        call lsyssc_getbase_double (rvectorDst,p_Dvec2)
        ! Let us go...
        if (rvectorSrc%NVAR .eq. 1) then
          !$omp parallel do default(shared) if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            p_Dvec2(i) = p_Dvec(i)*fmyscale/p_Fa(i)
          end do
          !$omp end parallel do
        else
          NVAR = rvectorSrc%NVAR
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Dvec2(NVAR*(i-1)+ivar) = p_Dvec(NVAR*(i-1)+ivar)&
                  *fmyscale/p_Fa(i)
            end do
          end do
          !$omp end parallel do
        end if

      case (ST_SINGLE)
        call lsyssc_getbase_single (rvectorSrc,p_Fvec)
        call lsyssc_getbase_single (rvectorDst,p_Fvec2)
        ! Let us go...
        if (rvectorSrc%NVAR .eq. 1) then
          !$omp parallel do default(shared) if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            p_Fvec2(i) = p_Fvec(i)*fmyscale/p_Fa(i)
          end do
          !$omp end parallel do
        else
          NVAR = rvectorSrc%NVAR
          !$omp parallel do default(shared) private(ivar) &
          !$omp if(rvectorSrc%NEQ > p_rperfconfig%NEQMIN_OMP)
          do i=1,rvectorSrc%NEQ
            do ivar = 1, NVAR
              p_Fvec2(NVAR*(i-1)+ivar) = p_Fvec(NVAR*(i-1)+ivar)&
                  *fmyscale/p_Fa(i)
            end do
          end do
          !$omp end parallel do
        end if

      case default
        call output_line('Unsupported vector precision!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported matrix precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
      call sys_halt()
    end select

  case default
    call output_line('Unsupported matrix format!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_invertedDiagMatVec')
    call sys_halt()
  end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_copyVector (rx,ry)

!<description>
  ! Copies vector data: ry = rx.
  ! Both vectors must have the same size. All structural data of rx is
  ! transferred to ry, so rx and ry are compatible to each other afterwards.
  ! If ry is empty, new memory is allocated automatically. Otherwise,
  ! ry is overwritten.
!</description>

!<input>
  ! Source vector
  type(t_vectorScalar),intent(in) :: rx
!</input>

!<inputoutput>
  ! Destination vector
  type(t_vectorScalar),intent(inout) :: ry
!</inputoutput>

!</subroutine>

    call lsyssc_duplicateVector (rx,ry,&
        LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_scaleVector (rx,c)

!<description>
  ! Scales a vector vector rx: rx = c * rx
!</description>

!<inputoutput>
  ! Source and destination vector
  type(t_vectorScalar), intent(inout) :: rx
!</inputoutput>

!<input>
  ! Multiplication factor
  real(DP), intent(in) :: c
!</input>

!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata

  ! Taje care of the data type!
  select case (rx%cdataType)
  case (ST_DOUBLE)
    ! Get the pointer and scale the whole data array.
    call lsyssc_getbase_double(rx,p_Ddata)
    call lalg_scaleVector(p_Ddata,c)

  case (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    call lsyssc_getbase_single(rx,p_Fdata)
    call lalg_scaleVector(p_Fdata,real(c,SP))

  case default
    call output_line('Unsupported data type!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scaleVector')
    call sys_halt()
  end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_clearVector (rx,dvalue)

!<description>
  ! Clears the block vector dx: Dx = 0 (or Dx = dvalue if dvalue is specified)
!</description>

!<inputoutput>
  ! Destination vector to be cleared
  type(t_vectorScalar), intent(inout) :: rx

  ! OPTIONAL: Value to write into the matrix.
  ! If not specified, all matrix entries are set to 0.0.
  ! If specified, all matrix entries are set to dvalue.
  real(DP), intent(in), optional :: dvalue
!</inputoutput>

!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_Dsource
  real(SP), dimension(:), pointer :: p_Ssource

  ! Take care of the data type
  select case (rx%cdataType)
  case (ST_DOUBLE)
    ! Get the pointer and scale the whole data array.
    call lsyssc_getbase_double(rx,p_Dsource)
    if (.not. present(dvalue)) then
      call lalg_clearVectorDble (p_Dsource)
    else
      call lalg_setVectorDble (p_Dsource,dvalue)
    end if

  case (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    call lsyssc_getbase_single(rx,p_Ssource)
    if (.not. present(dvalue)) then
      call lalg_clearVectorSngl (p_Ssource)
    else
      call lalg_setVectorSngl (p_Ssource,real(dvalue,SP))
    end if

  case default
    call output_line('Unsupported data type!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_clearVector')
    call sys_halt()
  end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_vectorLinearComb (rx,ry,cx,cy,rdest)

!<description>
  ! Performs a linear combination: rdest or ry = cx * rx  +  cy * ry
  ! Both vectors must be compatible to each other (same size, sorting
  ! strategy,...).
!</description>

!<input>
  ! First source vector
  type(t_vectorScalar), intent(in) :: rx

  ! Scaling factor for Dx
  real(DP), intent(in) :: cx

  ! Scaling factor for Dy
  real(DP), intent(in) :: cy
!</input>

!<inputoutput>
  ! Second source vector; also receives the result if rdest is not specified.
  type(t_vectorScalar), intent(inout), target :: ry

  ! OPTIONAL: Destination vector. If not specified, ry will be overwritten.
  type(t_vectorScalar), intent(inout), optional, target :: rdest
!</inputoutput>

!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_Dsource, p_Ddest
  real(SP), dimension(:), pointer :: p_Ssource, p_Sdest
  type(t_vectorScalar), pointer :: p_rdest

  ! The vectors must be compatible to each other.
  call lsyssc_isVectorCompatible (rx,ry)

  if (rx%cdataType .ne. ry%cdataType) then
    call output_line('Different data types not supported!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_vectorLinearComb')
    call sys_halt()
  end if

  p_rdest => ry
  if (present(rdest)) then

    p_rdest => rdest

    ! The vectors must be compatible to each other or rdest must be a clean
    ! vector.
    if (rdest%NEQ .ne. 0) then
      call lsyssc_isVectorCompatible (rx,rdest)

      if (rx%cdataType .ne. rdest%cdataType) then
        call output_line('Different data types not supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_vectorLinearComb')
        call sys_halt()
      end if
    end if

    ! Copy ry to rdest. Change rdest afterwards.
    call lsyssc_copyVector (ry,rdest)

  end if

  select case (rx%cdataType)
  case (ST_DOUBLE)
    ! Get the pointers and copy the whole data array.
    call lsyssc_getbase_double(rx,p_Dsource)
    call lsyssc_getbase_double(p_rdest,p_Ddest)

    call lalg_vectorLinearCombDble (p_Dsource,p_Ddest,cx,cy)

  case (ST_SINGLE)
    ! Get the pointers and copy the whole data array.
    call lsyssc_getbase_single(rx,p_Ssource)
    call lsyssc_getbase_single(p_rdest,p_Sdest)

    call lalg_vectorLinearCombSngl (p_Ssource,p_Sdest,real(cx,SP),real(cy,SP))

  case default
    call output_line('Unsupported data type!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_vectorLinearComb')
    call sys_halt()
  end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_copyMatrix (rsourceMatrix,rdestMatrix)

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
  type(t_matrixScalar),intent(in) :: rsourceMatrix
!</input>

!<inputoutput>
  ! Destination vector
  type(t_matrixScalar),intent(inout) :: rdestMatrix
!</inputoutput>

!</subroutine>

    call lsyssc_duplicateMatrix (rsourceMatrix,rdestMatrix,&
        LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE)

!  ! local variables
!  TYPE(t_matrixScalar) :: rtempMatrix
!  REAL(DP), DIMENSION(:), POINTER :: p_Dsource,p_Ddest
!  REAL(SP), DIMENSION(:), POINTER :: p_Fsource,p_Fdest
!  integer, DIMENSION(:), POINTER :: p_KcolSource,p_KcolDest
!  integer, DIMENSION(:), POINTER :: p_KldSource,p_KldDest
!  integer, DIMENSION(:), POINTER :: p_KdiagonalSource,p_KdiagonalDest
!
!  ! Is it possible at all to copy the matrix? Both matrices must have
!  ! the same size, otherwise the memory does not fit.
!  IF (rsourceMatrix%cmatrixFormat .NE. rdestMatrix%cmatrixFormat) THEN
!    call output_line('Different matrix formats not allowed!',&
!        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_copyMatrix')
!    CALL sys_halt()
!  END IF
!
!  IF (rsourceMatrix%NA .NE. rdestMatrix%NA) THEN
!    call output_line('Matrices have different size!',&
!        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_copyMatrix')
!    CALL sys_halt()
!  END IF
!
!  ! NEQ/NCOLS is irrelevant. It is only important that we have enough memory
!  ! and the same structure!
!
!  IF (rsourceMatrix%cdataType .NE. rdestMatrix%cdataType) THEN
!    call output_line('Matrices have different data types!',&
!        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_copyMatrix')
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
!    rdestMatrix%h_Da        = rtempMatrix%h_Da
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
!    CASE default
!      call output_line('Unsupported data type!',&
!        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_copyMatrix')
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
!    rdestMatrix%h_Da        = rtempMatrix%h_Da
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
!    CASE default
!      call output_line('Unsupported data type!',&
!        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_copyMatrix')
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
!    rdestMatrix%h_Da        = rtempMatrix%h_Da
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
!    CASE default
!      call output_line('Unsupported data type!',&
!        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_copyMatrix')
!      CALL sys_halt()
!    END SELECT
!
!  CASE default
!    call output_line('Unsupported matrix format!',&
!        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_copyMatrix')
!    CALL sys_halt()
!  END SELECT

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine lsyssc_transpMatStruct79double (nrow, ncol, Icol, Irow, &
                                             Itmp, IcolDest, IrowDest, &
                                             Ipermutation)

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
    integer, intent(in) :: nrow

    ! Number of columns in the source matrix
    integer, intent(in) :: ncol

    ! Column structure of the source matrix
    integer, dimension(:), intent(in) :: Icol

    ! Row structure of the source matrix
    integer, dimension(nrow+1), intent(in) :: Irow
!</input>

!<output>
    ! Column structure of the destination matrix
    ! The array must be of the same size as Icol!
    integer, dimension(:), intent(out) :: IcolDest

    ! Row structure of the destination matrix
    integer, dimension(ncol+1), intent(out) :: IrowDest

    ! OPTIONAL: Permutation array that specifies how to get the transposed
    ! matrix from the untransposed matrix. This is an array of length NA.
    integer, dimension(:), intent(out), optional :: Ipermutation
!</output>

!<inputoutput>
    ! Auxiliary array of size ncol
    integer, dimension(ncol), intent(inout) :: Itmp
!</inputoutput>
!</subroutine>

    !local variables
    integer(I32) :: i, j, isize, icolumn, ncolumn

    ! determin the number of matrix entries
    isize = Irow(nrow+1)-1

    ! clear row struxture of the destination matrix
    do i=1, ncol+1
      IrowDest(i) = 0
    end do

    ! Count how many entries <> 0 are in each column. Note this into
    ! the IrowDest array shifted by 1.
    do i=1, isize
      IrowDest(Icol(i)+1) = IrowDest(Icol(i)+1)+1
    end do

    ! Now build the final IrowDest by adding up the IrowDest entries.
    IrowDest(1) = 1
    do i=2, ncol+1
      IrowDest(i) = IrowDest(i)+IrowDest(i-1)
    end do

    ! Now IcolDest must be created. This requires another loop trough
    ! the matrix structure. Itmp receives the index of how many entries
    ! have been written to each row.

    ! clear auxiliary vector
    do i=1, ncol
      Itmp(i) = 0
    end do

    if (.not. present(Ipermutation)) then
      do i=1, nrow
        ncolumn = Irow(i+1)-Irow(i)
        do j=1, ncolumn
          ! Get the column of the item in question -> new row number.
          icolumn = Icol(Irow(i)+j-1)

          ! Rows get columns by transposing, therefore note i as column
          ! number in IcolDest
          IcolDest(IrowDest(icolumn)+Itmp(icolumn)) = i

          ! Increment running index of that row
          Itmp(icolumn) = Itmp(icolumn)+1
        end do
      end do
    else
      do i=1, nrow
        ncolumn = Irow(i+1)-Irow(i)
        do j=1, ncolumn
          ! Get the column of the item in question -> new row number.
          icolumn = Icol(Irow(i)+j-1)

          ! Rows get columns by transposing, therefore note i as column
          ! number in IcolDest
          IcolDest(IrowDest(icolumn)+Itmp(icolumn)) = i

          ! Save where this entry is moved to.
          Ipermutation(IrowDest(icolumn)+Itmp(icolumn)) = Irow(i)+j-1

          ! Increment running index of that row
          Itmp(icolumn) = Itmp(icolumn)+1
        end do
      end do
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine lsyssc_transpMatEntries79double (nrow, ncol, Da, Icol, Irow, &
      Itmp, DaDest, IrowDest, Ipermutation)

!<description>
    ! Auxiliary routine.
    ! This routine accepts a structure-7 or structure-9
    ! matrix and creates the transposed matrix from it.
    ! The 'destination matrix' is assumed to be already in transposed form.
    ! The routine only copies the entries in transposed form from the source
    ! to the destination matrix, ignoring the structure of the destination.
    !
    ! The resulting matrix is of format 9, i.e. not pivoted anymore with
    ! the diagonal entry in front if a structure-7 matrix is to be transposed!
!</description>

!<input>
    ! Number of rows in the source matrix
    integer, intent(in) :: nrow

    ! Number of columns in the source matrix
    integer, intent(in) :: ncol

    ! The entries of the source matrix
    real(DP), dimension(:), intent(in) :: Da

    ! Column structure of the source matrix
    integer, dimension(:), intent(in) :: Icol

    ! Row structure of the source matrix.
    ! DIMENSION(nrow+1)
    integer, dimension(:), intent(in) :: Irow
!</input>

!<output>
    ! The entries of the destination matrix
    ! The array must be of the same size as Da or Icol
    real(DP), dimension(:), intent(out) :: DaDest

    ! Temporary array. Row structure of the destination matrix.
    ! DIMENSION(ncol+1)
    integer, dimension(:), intent(out) :: IrowDest

    ! OPTIONAL: Permutation array that specifies how to get the transposed
    ! matrix from the untransposed matrix. This is an array of length NA.
    integer, dimension(:), intent(out), optional :: Ipermutation
!</output>

!<inputoutput>
    ! Auxiliary array.
    ! DIMENSION(ncol)
    integer, dimension(:), intent(inout) :: Itmp
!</inputoutput>

!</subroutine>

    !local variables
    integer(I32) :: i, j, isize, icolumn, ncolumn

    if ((size(Itmp) .ne. ncol) .or. (size(IrowDest) .ne. ncol+1) .or. &
        (size(Irow) .ne. nrow+1)) then
      call output_line('Array parameters have wrong size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transpMatEntries79double')
      call sys_halt()
    end if

    ! determine the number of matrix entries
    isize = Irow(nrow+1)-1

    ! clear row structure of the destination matrix
    do i=1, ncol+1
      IrowDest(i) = 0
    end do

    ! Count how many entries <> 0 are in each column. Note this
    ! into the IrowDest array.
    ! This requires one loop through the matrix structure. The
    ! corresponding number is written into KLDD, shifted by 1
    ! (i.e. IrowDest(2) is the number of entries of the 1st column).
    ! This helps to create the real IrowDest more easily later.
    do i=1, isize
      IrowDest(Icol(i)+1) = IrowDest(Icol(i)+1)+1
    end do

    ! Now build the real IrowDest. This consists of indices, where
    ! each row starts. Row 1 starts at position 1:
    IrowDest(1) = 1

    ! Adding the number of entries IrowDest(i) in the row to
    ! IrowDest(i-1) gives the new row pointer of row i of the transposed
    ! matrix.
    do i=2, ncol+1
      IrowDest(i) = IrowDest(i)+IrowDest(i-1)
    end do

    ! That is it for IrowDest.
    ! Now, the matrix must be copied to DaDest in transposed form.
    ! This requires another loop trough the matrix structure.
    ! Itmp receives the index of how many entries have been written
    ! to each row.

    ! clear auxiliary vector
    do i=1, ncol
      Itmp(i) = 0
    end do

    if (.not. present(Ipermutation)) then
      do i=1, nrow
        ! Loop through the row
        ncolumn = Irow(i+1)-Irow(i)
        do j=1, ncolumn
          ! Get the column of the item in question -> new row number.
          icolumn = Icol(Irow(i)+j-1)

          ! Copy the matrix entry:
          DaDest(IrowDest(icolumn)+Itmp(icolumn)) = Da(Irow(i)+j-1)

          ! Increment running index of that row
          Itmp(icolumn) = Itmp(icolumn)+1
        end do
      end do
    else
      do i=1, nrow
        ! Loop through the row
        ncolumn = Irow(i+1)-Irow(i)
        do j=1, ncolumn
          ! Get the column of the item in question -> new row number.
          icolumn = Icol(Irow(i)+j-1)

          ! Copy the matrix entry:
          DaDest(IrowDest(icolumn)+Itmp(icolumn)) = Da(Irow(i)+j-1)

          ! Save where this entry is moved to.
          Ipermutation(IrowDest(icolumn)+Itmp(icolumn)) = Irow(i)+j-1

          ! Increment running index of that row
          Itmp(icolumn) = Itmp(icolumn)+1
        end do
      end do
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine lsyssc_transpMat79double (nrow, ncol, Da, Icol, Irow, &
      Itmp, DaDest, IcolDest, IrowDest, Ipermutation)

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
    integer, intent(in) :: nrow

    ! Number of columns in the source matrix
    integer, intent(in) :: ncol

    ! The entries of the source matrix
    real(DP), dimension(:), intent(in) :: Da

    ! Column structure of the source matrix
    integer, dimension(:), intent(in) :: Icol

    ! Row structure of the source matrix.
    ! DIMENSION(nrow+1)
    integer, dimension(:), intent(in) :: Irow
!</input>

!<output>
    ! The entries of the destination matrix
    ! The array must be of the same size as Da or Icol
    real(DP), dimension(:), intent(out) :: DaDest

    ! Column structure of the destination matrix
    ! The array must be of the same size as Icol!
    integer, dimension(:), intent(out) :: IcolDest

    ! Row structure of the destination matrix.
    ! DIMENSION(ncol+1)
    integer, dimension(:), intent(out) :: IrowDest

    ! OPTIONAL: Permutation array that specifies how to get the transposed
    ! matrix from the untransposed matrix. This is an array of length NA.
    integer, dimension(:), intent(out), optional :: Ipermutation
!</output>

!<inputoutput>
    ! Auxiliary array.
    ! DIMENSION(ncol)
    integer, dimension(:), intent(inout) :: Itmp
!</inputoutput>

!</subroutine>

    !local variables
    integer(I32) :: i, j, isize, icolumn, ncolumn

    if ((size(Itmp) .ne. ncol) .or. (size(IrowDest) .ne. ncol+1) .or. &
        (size(Irow) .ne. nrow+1)) then
      call output_line('Array parameters have wrong size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transpMat79double')
      call sys_halt()
    end if

    ! determine the number of matrix entries
    isize = Irow(nrow+1)-1

    ! clear row structure of the destination matrix
    do i=1, ncol+1
      IrowDest(i) = 0
    end do

    ! Count how many entries <> 0 are in each column. Note this
    ! into the IrowDest array.
    ! This requires one loop through the matrix structure. The
    ! corresponding number is written into KLDD, shifted by 1
    ! (i.e. IrowDest(2) is the number of entries of the 1st column).
    ! This helps to create the real IrowDest more easily later.
    do i=1, isize
      IrowDest(Icol(i)+1) = IrowDest(Icol(i)+1)+1
    end do

    ! Now build the real IrowDest. This consists of indices, where
    ! each row starts. Row 1 starts at position 1:
    IrowDest(1) = 1

    ! Adding the number of entries IrowDest(i) in the row to
    ! IrowDest(i-1) gives the new row pointer of row i of the transposed
    ! matrix.
    do i=2, ncol+1
      IrowDest(i) = IrowDest(i)+IrowDest(i-1)
    end do

    ! That is it for IrowDest.
    ! Now IcolDest must be created. This requires another loop trough
    ! the matrix structure. Itmp receives the index of how many entries
    ! have been written to each row.

    ! clear auxiliary vector
    do i=1, ncol
      Itmp(i) = 0
    end do

    if (.not. present(Ipermutation)) then
      do i=1, nrow
        ! Loop through the row
        ncolumn = Irow(i+1)-Irow(i)
        do j=1, ncolumn
          ! Get the column of the item in question -> new row number.
          icolumn = Icol(Irow(i)+j-1)

          ! Rows get columns by transposing, therefore note i as column
          ! number in IcolDest
          IcolDest(IrowDest(icolumn)+Itmp(icolumn)) = i

          ! Copy the matrix entry:
          DaDest(IrowDest(icolumn)+Itmp(icolumn)) = Da(Irow(i)+j-1)

          ! Increment running index of that row
          Itmp(icolumn) = Itmp(icolumn)+1
        end do
      end do
    else
      do i=1, nrow
        ! Loop through the row
        ncolumn = Irow(i+1)-Irow(i)
        do j=1, ncolumn
          ! Get the column of the item in question -> new row number.
          icolumn = Icol(Irow(i)+j-1)

          ! Rows get columns by transposing, therefore note i as column
          ! number in IcolDest
          IcolDest(IrowDest(icolumn)+Itmp(icolumn)) = i

          ! Copy the matrix entry:
          DaDest(IrowDest(icolumn)+Itmp(icolumn)) = Da(Irow(i)+j-1)

          ! Save where this entry is moved to.
          Ipermutation(IrowDest(icolumn)+Itmp(icolumn)) = Irow(i)+j-1

          ! Increment running index of that row
          Itmp(icolumn) = Itmp(icolumn)+1
        end do
      end do
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_pivotiseMatrix7double (neq, Icol, Irow, Da, berror)

!<description>
    ! Auxiliary routine.
    ! This routine repivots a matrix in structure 7: The diagonal element
    ! of each row is moved to the front.
!</description>

!<input>
    ! Number of rows/columns/equations
    integer, intent(in) :: neq
!</input>

!<output>
    ! If .TRUE. a diagonal element was not found
    logical, optional, intent(out) :: berror
!</output>

!<inputoutput>
    ! OPTIONAL: Matrix entries.
    ! If not specified, only the matrix structure is pivotised.
    real(DP), dimension(:), intent(inout), optional :: Da

    ! Matrix column structure
    integer, dimension(:), intent(inout) :: Icol

    ! Matrix row structure
    integer, dimension(:), intent(inout) :: Irow
!</inputoutput>
!</subroutine>

    !local variables
    integer(I32) :: i, j
    logical :: bnodiag

    bnodiag = .false.

    if (present(Da)) then
      ! Structure + entries
      !
      ! Loop through the rows
      do i=1, neq
        ! Is the line already pivoted? Most unlikely, but we test for sure
        if (Icol(Irow(i)) .ne. i) then

          ! Find the position of the diagonal element
          bnodiag = .true.
          do j=Irow(i), Irow(i+1)-1
            if (j .eq. i) then
              bnodiag = .false.
              exit
            end if
          end do
          ! Oops, diagonal element not found - cancel
          if (bnodiag) exit

          ! Ringshift the slice Icol(Irow(i):j) by 1 to the right.
          ! This puts the diagonal element to the front
          ! The same operation is also done for Da(Irow(i):j)
          Icol(Irow(i):j) = cshift(Icol(Irow(i):j),-1)
          Da  (Irow(i):j) = cshift(Da  (Irow(i):j),-1)
        end if
      end do

    else
      ! Only structure
      !
      ! Loop through the rows
      do i=1, neq
        ! Is the line already pivoted? Most unlikely, but we test for sure
        if (Icol(Irow(i)) .ne. i) then

          ! Find the position of the diagonal element
          bnodiag = .true.
          do j=Irow(i), Irow(i+1)-1
            if (j .eq. i) then
              bnodiag = .false.
              exit
            end if
          end do
          ! Oops, diagonal element not found - cancel
          if (bnodiag) exit

          ! Ringshift the slice Icol(Irow(i):j) by 1 to the right.
          ! This puts the diagonal element to the front
          ! The same operation is also done for Da(Irow(i):j)
          Icol(Irow(i):j) = cshift(Icol(Irow(i):j),-1)
          Da  (Irow(i):j) = cshift(Da  (Irow(i):j),-1)
        end if
      end do

    end if

    ! Report error status
    if (present(berror)) berror = bnodiag

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_transposeMatrixInSitu (rmatrix,itransFlag)

!<description>
  ! This routine in-situ transposes a matrix rmatrix and creates the transposed
  ! matrix in rtransposedMatrix.
  ! itransFlag decides (if specified) how the creation of the
  ! transposed matrix is to be performed.
!</description>

!<input>
  ! OPTIONAL: transposed-flag; one of the LSYSSC_TR_xxxx constants:
  ! =LSYSSC_TR_ALL or not specified: transpose the matrix completely.
  ! =LSYSSC_TR_STRUCTURE           : Transpose only the matrix structure;
  !     the content of the transposed matrix is invalid afterwards,
  !     but is not released/destroyed.
  ! =LSYSSC_TR_VIRTUAL             : Actually do not touch the matrix
  !     structure, but invert the "transposed" flag in imatrixSpec.
  !     rmatrix is marked as transposed without modifying the structures,
  !     and all matrix-vector operations will be performed with the transposed
  !     matrix.
  !     But caution: there may be some algorithms that do not work with such
  !       "virtually" transposed matrices!
  ! =LSYSSC_TR_VIRTUALCOPY         : The same as LSYSSC_TR_VIRTUAL, but
  !     creates a duplicate of the source matrix in memory, thus resulting
  !     in rtransposedMatrix being a totally independent matrix.
  !     So afterwards, rmatrix is an independent, virtually transposed matrix.
  integer, intent(in), optional :: itransFlag
!</input>

!<inputoutput>
  ! The matrix to be transposed.
  type(t_matrixScalar),intent(inout) :: rmatrix
!</inputoutput>
!</subroutine>

    type(t_matrixScalar) :: rmatrix2
    integer :: itrans

    ! Transpose the matrix to a temporary copy.
    call lsyssc_transposeMatrix (rmatrix,rmatrix2,itransFlag)

    ! Next thing is to release the old matrix; but that is not as easy
    ! as calling lsyssc_releaseMatrix! The point is that the source
    ! matrix may or may not share its structure/content with another
    ! one -- and the new matrix may do this as well!
    !
    ! We have to separate some cases...

    itrans = LSYSSC_TR_ALL
    if (present(itransFlag)) itrans = itransFlag

    select case (itrans)
    case (LSYSSC_TR_ALL,LSYSSC_TR_VIRTUALCOPY)
      ! These are the simple cases: The new matrix is independent.
      ! Release the old one.
      call lsyssc_releaseMatrix(rmatrix)

    case (LSYSSC_TR_VIRTUAL)
      ! This case slightly harder. The new matrix shares both,
      ! structure and content with the old one. The old one may share
      ! its structure or content with another one. We have to transfer
      ! the copy flags from the old matrix here.
      rmatrix2%imatrixSpec = ior(iand(rmatrix2%imatrixSpec,not(LSYSSC_MSPEC_ISCOPY)),&
                                 iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY))

      ! rmatrix2 has the handles of rmatrix but we want to release
      ! rmatrix now. In this case we have to prevent the data arrays
      ! from being released in case rmatrix is the old owner of the
      ! arrays.
      rmatrix%imatrixSpec = ior(rmatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY)

      ! Now we can savely release the source matrix without destroying
      ! anything.
      call lsyssc_releaseMatrix(rmatrix)

    case (LSYSSC_TR_STRUCTURE)
      ! Ok, this case is a bit strange, somehow, and I am not sure if it
      ! really works this way, but it should :-)
      !
      ! The above lsyssc_transposeMatrix transposed the structure,
      ! so the new matrix has a complete new structure -- we do not
      ! have to take care about.
      !
      ! The content of rmatrix may belong to rmatrix or not.
      ! The rule is that lsyssc_transposeMatrixInSitu must keeps the
      ! content -- it gets invalid, but it is not destroyed. So we
      ! must make sure that releasing rmatrix must not destroy
      ! the content. On the other hand, rmatrix2 must get the
      ! copy state of the content from rmatrix.
      !
      ! So at first make sure, rmatrix2 is the owner of the content
      ! (if rmatrix was the original owner) or takes the copy state
      ! from rmatrix.
      rmatrix2%imatrixSpec = ior(iand(rmatrix2%imatrixSpec,&
                                      not(LSYSSC_MSPEC_CONTENTISCOPY)),&
                                 iand(rmatrix%imatrixSpec,&
                                      LSYSSC_MSPEC_CONTENTISCOPY))

      ! Prevent lsyssc_releaseMatrix from releasing the content if it is still there
      rmatrix%imatrixSpec = ior(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY)

      ! Now we can release the old matrix.
      call lsyssc_releaseMatrix(rmatrix)

    case (LSYSSC_TR_CONTENT)
      ! Sorry, that is not possible. The transpose routine needs the structure of
      ! the original matrix, otherwise it cannot compute how to transpose
      ! the entries!
      call output_line("TRansposing matrix entries without structure is not possible!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_transposeMatrixInSitu")
      call sys_halt()

    end select

    ! Create a hard-copy of rmatrix.
    ! This is one of the very rare cases where we on purpose assign a matrix
    ! structure to another...
    rmatrix = rmatrix2

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_transposeMatrixDirect (rsourceMatrix,rdestMatrix,Ipermutation)

!<description>
  ! Transposes the matrix entries of a matrix according to a matrix permutation
  ! array. Ipermutation must be a permutation array computed by
  ! lsyssc_transposeMatrix. lsyssc_transposeMatrixQuick will then permute
  ! the entries in the matrix rsourceMatrix according to this permutation and
  ! will save the result to rdestMatrix.
  ! The routine assumes that rdestMatrix describes the structure of the
  ! transposed of rsourceMatrix, generated by lsyssc_transposeMatrix.
!</description>

!<input>
  ! The matrix to be transposed.
  type(t_matrixScalar),intent(in) :: rsourceMatrix

  ! Permutation array that specifies how to get the transposed matrix
  ! from the untransposed matrix. This is an array of length NA.
  integer, dimension(:), intent(in) :: Ipermutation
!</input>

!<inputoutput>
  ! The matrix structure that receives the transposed matrix.
  type(t_matrixScalar),intent(inout) :: rdestMatrix
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_DaSource,p_DaDest

    select case (rsourceMatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9,LSYSSC_MATRIX7,LSYSSC_MATRIXD)

      if (rsourceMatrix%NA .eq. 0) then
        call output_line("Source matrix empry!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_transposeMatrixDirect")
        call sys_halt()
      end if

      if (rsourceMatrix%cdataType .ne. ST_DOUBLE) then
        call output_line("Invalid data type!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_transposeMatrixDirect")
        call sys_halt()
      end if

      if ((rsourceMatrix%NA .ne. rdestMatrix%NA) .or. &
          (rsourceMatrix%cdataType .ne. rdestMatrix%cdataType)) then
        call output_line("Source and destination matrix incompatible!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_transposeMatrixDirect")
        call sys_halt()
      end if

      select case (rsourceMatrix%cdataType)
      case (ST_DOUBLE)

        ! Get the data arrays
        call lsyssc_getbase_double(rsourceMatrix,p_DaSource)
        if (rdestMatrix%h_Da .eq. ST_NOHANDLE) then
          call lsyssc_allocEmptyMatrix (rdestMatrix,LSYSSC_SETM_UNDEFINED)
        end if
        call lsyssc_getbase_double(rdestMatrix,p_DaDest)

        ! Permute
        call lalg_vectorSortDble (p_DaSource, p_DaDest, Ipermutation)

      case default
        call output_line("Invalid data type!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_transposeMatrixDirect")
        call sys_halt()
      end select

    case default
      call output_line("Unsupported matrix format!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_transposeMatrixDirect")
      call sys_halt()
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_transposeMatrix (rmatrix,rtransposedMatrix,itransFlag,&
      Ipermutation)

!<description>
  ! This routine transposes a matrix rmatrix and creates the transposed
  ! matrix in rtransposedMatrix.
  ! itransFlag decides (if specified) how the creation of the
  ! transposed matrix is to be performed.
!</description>

!<input>
  ! OPTIONAL: transposed-flag; one of the LSYSSC_TR_xxxx constants:
  ! =LSYSSC_TR_ALL or not specified: transpose the matrix completely.
  !     If there is a matrix in rtransposedMatrix, it is overwritten
  !       or an error is thrown, depending on whether the
  !       structural data (NA,...) of rtransposedMatrix matches rmatrix.
  !     If there is no matrix in rtransposedMatrix, a new matrix is
  !       created.
  ! =LSYSSC_TR_STRUCTURE           : Transpose only the matrix structure;
  !     the content of the transposed matrix is invalid afterwards.
  !     If there is a matrix in rtransposedMatrix, the structure is
  !       overwritten or an error is thrown, depending on whether the
  !       structural data (NA,...) of rtransposedMatrix matches rmatrix.
  !     If there is no matrix in rtransposedMatrix, a new matrix structure
  !       without entries is created.
  ! =LSYSSC_TR_CONTENT             : Transpose only the matrix content/entries;
  !     the structure of the transposed matrix is invalid afterwards.
  !     If there is already a matrix in rtransposedMatrix, the content is
  !       overwritten or an error is thrown, depending on whether the
  !       structural data (NA,...) of rtransposedMatrix matches rmatrix.
  !       Note that rtransposedMatrix may contain a valid structure
  !       of the transposed matrix. In this case, rtransposedMatrix will
  !       receive the transposed matrix entries, thus resulting in
  !       a proper transposed matrix.
  !     If there is no matrix in rtransposedMatrix, a new matrix content
  !       without structure is created.
  ! =LSYSSC_TR_VIRTUAL             : Actually do not touch the matrix
  !     structure, but invert the 'transposed' flag in imatrixSpec.
  !
  !     rtransposedMatrix is created by copying the data handles from
  !     rmatrix - any previous data in rtransposedMatrix is overwritten.
  !     Afterwards, rtransposedMatrix is marked as transposed
  !     without modifying the structures, and all matrix-vector operations
  !     will be performed with the transposed matrix.
  !     But caution: there may be some algorithms that do not work with such
  !       'virtually' transposed matrices!
  ! =LSYSSC_TR_VIRTUALCOPY         : The same as LSYSSC_TR_VIRTUAL, but
  !     creates a duplicate of the source matrix in memory, thus resulting
  !     in rtransposedMatrix being a totally independent matrix.
  integer, intent(in), optional :: itransFlag

  ! The matrix to be transposed.
  type(t_matrixScalar),intent(in) :: rmatrix
!</input>

!<inputoutput>
  ! The transposed matrix.
  ! If the structure is empty, a new matrix is created.
  ! If the structure exists and structural data matches rmatrix,
  !   it is overwritten.
  ! If the structure exists and structural data does not match rmatrix,
  !   an error is thrown.
  type(t_matrixScalar),intent(inout) :: rtransposedMatrix
!</inputoutput>

!<output>
    ! OPTIONAL: Permutation array that specifies how to get the transposed matrix
    ! from the untransposed matrix. This is an array of length NA.
    integer, dimension(:), intent(out), optional :: Ipermutation
!</output>

!</subroutine>

    ! local variables
    integer :: itrans, h_Itemp, h_Kld, h_Kcol
    integer :: ntemp
    integer, dimension(:), pointer :: p_KcolSource,p_KcolDest
    integer, dimension(:), pointer :: p_KldSource,p_KldDest
    integer, dimension(:), pointer :: p_Kdiagonal
    integer, dimension(:), pointer :: p_Itemp
    real(DP), dimension(:), pointer :: p_DaSource,p_DaDest

    itrans = LSYSSC_TR_ALL
    if (present(itransFlag)) itrans = itransFlag

    ! How to perform the creation of the transposed matrix?
    if (itrans .eq. LSYSSC_TR_VIRTUAL) then

      ! Duplicate the matrix structure, copy all handles.
      ! Do not allocate any memory.
      call lsyssc_duplicateMatrix (rmatrix,rtransposedMatrix,&
                                    LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      ! Only change the 'transposed' flag in imatrixSpec
      rtransposedMatrix%imatrixSpec = &
        ieor(rtransposedMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED)

      ! Exchange NEQ and NCOLS as these always describe the actual matrix.
      ntemp = rtransposedMatrix%NEQ
      rtransposedMatrix%NEQ = rtransposedMatrix%NCOLS
      rtransposedMatrix%NCOLS = ntemp

      ! Switch the discretisation structures
      rtransposedMatrix%p_rspatialDiscrTrial => rmatrix%p_rspatialDiscrTest
      rtransposedMatrix%p_rspatialDiscrTest => rmatrix%p_rspatialDiscrTrial

      ! That is it.
      return
    end if

    if (itrans .eq. LSYSSC_TR_VIRTUALCOPY) then

      ! Duplicate the complete matrix structure in memory
      call lsyssc_duplicateMatrix (rmatrix,rtransposedMatrix,&
                                    LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)

      ! Only change the 'transposed' flag in imatrixSpec
      rtransposedMatrix%imatrixSpec = &
        ieor(rtransposedMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED)

      ! Exchange NEQ and NCOLS as these always describe the actual matrix.
      rtransposedMatrix%NEQ = rmatrix%NCOLS
      rtransposedMatrix%NCOLS = rmatrix%NEQ

      ! Switch the discretisation structures
      rtransposedMatrix%p_rspatialDiscrTrial => rmatrix%p_rspatialDiscrTest
      rtransposedMatrix%p_rspatialDiscrTest => rmatrix%p_rspatialDiscrTrial

      ! That is it.
      return
    end if

    ! Does the destination matrix exist?
    if (rtransposedMatrix%NEQ .ne. 0) then

      ! Make sure the destination matrix can accept all data of the source
      ! matrix. Matrix format and size must match.
      if (rMatrix%cmatrixFormat .ne. rtransposedMatrix%cmatrixFormat) then
        call output_line('Different matrix formats not allowed!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
        call sys_halt()
      end if

      if (rmatrix%NA .ne. rtransposedMatrix%NA) then
        call output_line('Matrices have different size!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
        call sys_halt()
      end if

      ! NEQ/NCOLS is irrelevant. It is only important that we have enough memory
      ! and the same structure!

      if (itrans .eq. LSYSSC_TR_ALL) then
        if (rmatrix%cdataType .ne. rtransposedMatrix%cdataType) then
          call output_line('Matrices have different data types!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
          call sys_halt()
        end if
      end if

    else
      ! Otherwise, we must create a new matrix that accepts the transposed
      ! data.
      ! Copy the matrix structure, reset all handles, so we have a 'template'
      ! matrix with the same data as the original one but without dynamic
      ! information attached.
      call lsyssc_duplicateMatrix (rmatrix,rtransposedMatrix,&
                                    LSYSSC_DUP_TEMPLATE,LSYSSC_DUP_TEMPLATE)
    end if

    ! Otherwise, check that we are able to create the transposed matrix at all.
    if (itrans .eq. LSYSSC_TR_ALL) then
      if (rmatrix%cdataType .ne. ST_DOUBLE) then
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
        call sys_halt()
      end if
    end if

    ! Now what should we do...
    select case (itrans)

    case (LSYSSC_TR_STRUCTURE)

      ! Only the structure is to be transposed. Which structure do we have?
      select case (rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9)
        ! Is the source matrix saved tranposed? In that case simply copy
        ! the data arrays, this results in the transposed matrix :)
        if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then

          call lsyssc_auxcopy_Kcol (rmatrix,rtransposedMatrix)
          call lsyssc_auxcopy_Kld (rmatrix,rtransposedMatrix)
          call lsyssc_auxcopy_Kdiagonal (rmatrix,rtransposedMatrix)

          rtransposedMatrix%imatrixSpec = &
            iand(rtransposedMatrix%imatrixSpec,not(LSYSSC_MSPEC_TRANSPOSED))

          if (present(Ipermutation)) then
            call output_line('Calculation of the permutation not possible!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
            call sys_halt()
          end if

        else
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          if (rtransposedMatrix%h_Kcol .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                ST_INT, rtransposedMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
          end if

          if (rtransposedMatrix%h_Kld .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1, &
                ST_INT, rtransposedMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
          end if

          ! Get Kcol/Kld. Do not use the lsyssc_getbase routines as we just created
          ! the arrays for a matrix which is not yet valid!
          call storage_getbase_int(rmatrix%h_Kcol,p_KcolSource)
          call storage_getbase_int(rmatrix%h_Kld,p_KldSource)

          call storage_getbase_int(rtransposedMatrix%h_Kcol,p_KcolDest)
          call storage_getbase_int(rtransposedMatrix%h_Kld,p_KldDest)

          ! We need a temporary array
          call storage_new ('lsyssc_transposeMatrix','Itemp',rmatrix%NCOLS,&
                            ST_INT, h_Itemp, ST_NEWBLOCK_NOINIT)
          call storage_getbase_int (h_Itemp,p_Itemp)

          ! Calculate the transposed matrix structure:
          call lsyssc_transpMatStruct79double (rmatrix%NEQ, rmatrix%NCOLS, &
                     p_KcolSource, p_KldSource, &
                     p_Itemp, p_KcolDest, p_KldDest,Ipermutation)

          ! Release the temp array
          call storage_free (h_Itemp)

          ! (Re-)calculate Kdiagonal
          if (rtransposedMatrix%h_Kdiagonal .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Kdiagonal', rmatrix%NCOLS, &
                ST_INT, rtransposedMatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
          end if
          call storage_getbase_int(rtransposedMatrix%h_Kdiagonal,p_Kdiagonal)
          call lsyssc_rebuildKdiagonal (p_KcolDest, p_KldDest, p_Kdiagonal, &
                                        rmatrix%NCOLS)

        end if

      case (LSYSSC_MATRIX7)
        ! Is the source matrix saved tranposed? In that case simply copy
        ! the data arrays, this results in the transposed matrix :)
        if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
          call lsyssc_auxcopy_Kcol (rmatrix,rtransposedMatrix)
          call lsyssc_auxcopy_Kld (rmatrix,rtransposedMatrix)
          rtransposedMatrix%imatrixSpec = &
            iand(rtransposedMatrix%imatrixSpec,not(LSYSSC_MSPEC_TRANSPOSED))

          if (present(Ipermutation)) then
            call output_line('Calculation of the permutation not possible!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
            call sys_halt()
          end if

        else
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          if (rtransposedMatrix%h_Kcol .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                ST_INT, rtransposedMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
          end if

          if (rtransposedMatrix%h_Kld .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1, &
                ST_INT, rtransposedMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
          end if

          ! Get Kcol/Kld. Do not use the lsyssc_getbase routines as we just created
          ! the arrays for a matrix which is not yet valid!
          call storage_getbase_int(rmatrix%h_Kcol,p_KcolSource)
          call storage_getbase_int(rmatrix%h_Kld,p_KldSource)

          call storage_getbase_int(rtransposedMatrix%h_Kcol,p_KcolDest)
          call storage_getbase_int(rtransposedMatrix%h_Kld,p_KldDest)

          ! We need a temporary array
          call storage_new ('lsyssc_transposeMatrix','Itemp',rmatrix%NA,&
                            ST_INT, h_Itemp, ST_NEWBLOCK_NOINIT)
          call storage_getbase_int (h_Itemp,p_Itemp)

          ! Calculate the transposed matrix structure:
          call lsyssc_transpMatStruct79double (rmatrix%NEQ, rmatrix%NCOLS, &
                     p_KcolSource, p_KldSource, &
                     p_Itemp, p_KcolDest, p_KldDest, Ipermutation)

          ! Release the temp array
          call storage_free (h_Itemp)

          ! Pivotise the matrix to move the diagonal element to the front.
          call lsyssc_pivotiseMatrix7double (rmatrix%NEQ, p_KcolDest, p_KldDest)
        end if

      case (LSYSSC_MATRIXD)
        ! Nothing to do

        if (present(Ipermutation)) then
          call output_line('Calculation of the permutation not possible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
          call sys_halt()
        end if

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
        call sys_halt()
      end select

      ! Exchange NEQ and NCOLS as these always describe the actual matrix.
      rtransposedMatrix%NEQ = rmatrix%NCOLS
      rtransposedMatrix%NCOLS = rmatrix%NEQ

      ! Switch the discretisation structures
      rtransposedMatrix%p_rspatialDiscrTrial => rmatrix%p_rspatialDiscrTest
      rtransposedMatrix%p_rspatialDiscrTest => rmatrix%p_rspatialDiscrTrial

    case (LSYSSC_TR_CONTENT)
      ! Transpose the full matrix - structure+content
      select case (rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9)
        ! Is the source matrix saved tranposed? In that case simply copy
        ! the data arrays, this results in the transposed matrix :)
        if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then

          call lsyssc_auxcopy_da (rmatrix,rtransposedMatrix)
          rtransposedMatrix%imatrixSpec = &
            iand(rtransposedMatrix%imatrixSpec,not(LSYSSC_MSPEC_TRANSPOSED))

          if (present(Ipermutation)) then
            call output_line('Calculation of the permutation not possible!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
            call sys_halt()
          end if

        else
          ! We really have to do some work now :)
          ! Allocate a new Kld and Da.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          call storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1, &
              ST_INT, h_Kld,ST_NEWBLOCK_NOINIT)

          if (rtransposedMatrix%h_Da .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Da', rmatrix%NA, &
                rtransposedMatrix%cdataType, rtransposedMatrix%h_Da,&
                ST_NEWBLOCK_NOINIT)
          end if

          ! Get Kcol/Kld. Do not use the lsyssc_getbase routines as we just created
          ! the arrays for a matrix which is not yet valid!
          call storage_getbase_int(rmatrix%h_Kcol,p_KcolSource)
          call storage_getbase_int(rmatrix%h_Kld,p_KldSource)

          call storage_getbase_int(h_Kld,p_KldDest)

          ! Get the data array(s)
          call lsyssc_getbase_double(rMatrix,p_DaSource)
          call storage_getbase_double(rtransposedMatrix%h_Da,p_DaDest)

          ! We need a temporary array
          call storage_new ('lsyssc_transposeMatrix','Itemp',rmatrix%NCOLS,&
                            ST_INT, h_Itemp, ST_NEWBLOCK_NOINIT)
          call storage_getbase_int (h_Itemp,p_Itemp)

          ! Calculate the transposed matrix structure:
          call lsyssc_transpMatEntries79double (rmatrix%NEQ, rmatrix%NCOLS, &
                p_DaSource, p_KcolSource, p_KldSource, &
                p_Itemp, p_DaDest, p_KldDest, Ipermutation)

          ! Release the temp array
          call storage_free (h_Itemp)
          call storage_free (h_Kld)

        end if

      case (LSYSSC_MATRIX7)
        ! For matrix 7 we have a slight problem. Here, the diagonal
        ! entry is in the front, while the transpose routine will produce
        ! us a matrix with the diagonal entry somewhere in the middle.
        ! Thus, we have (with that routine) no other chance than to
        ! produce the full transposed matrix (including structure), figure
        ! our where the diagonal is and move it to the front.

        ! Is the source matrix saved tranposed? In that case simply copy
        ! the data arrays, this results in the transposed matrix :)
        if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then

          if (rmatrix%h_Da .ne. ST_NOHANDLE) then
            call lsyssc_auxcopy_da (rmatrix,rtransposedMatrix)
          end if

          rtransposedMatrix%imatrixSpec = &
            iand(rtransposedMatrix%imatrixSpec,not(LSYSSC_MSPEC_TRANSPOSED))

          if (present(Ipermutation)) then
            call output_line('Calculation of the permutation not possible!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
            call sys_halt()
          end if

        else
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          call storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                            ST_INT, h_Kcol,ST_NEWBLOCK_NOINIT)
          call storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1, &
                            ST_INT, h_Kld,ST_NEWBLOCK_NOINIT)

          if (rtransposedMatrix%h_Da .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Da', rmatrix%NA, &
                rtransposedMatrix%cdataType, rtransposedMatrix%h_Da,&
                ST_NEWBLOCK_NOINIT)
          end if

          ! Get Kcol/Kld. Do not use the lsyssc_getbase routines as we just created
          ! the arrays for a matrix which is not yet valid!
          call storage_getbase_int(rmatrix%h_Kcol,p_KcolSource)
          call storage_getbase_int(rmatrix%h_Kld,p_KldSource)

          call storage_getbase_int(h_Kcol,p_KcolDest)
          call storage_getbase_int(h_Kld,p_KldDest)

          ! Get the data array(s)
          call lsyssc_getbase_double(rMatrix,p_DaSource)
          call storage_getbase_double(rtransposedMatrix%h_Da,p_DaDest)

          ! We need a temporary array
          call storage_new ('lsyssc_transposeMatrix','Itemp',rmatrix%NCOLS,&
                            ST_INT, h_Itemp, ST_NEWBLOCK_NOINIT)
          call storage_getbase_int (h_Itemp,p_Itemp)

          ! Calculate the transposed matrix structure:
          call lsyssc_transpMat79double (rmatrix%NEQ, rmatrix%NCOLS, &
                p_DaSource, p_KcolSource, p_KldSource, &
                p_Itemp, p_DaDest, p_KcolDest, p_KldDest, Ipermutation)

          ! Release the temp array
          call storage_free (h_Itemp)

          ! Pivotise the matrix to move the diagonal element to the front.
          call lsyssc_pivotiseMatrix7double (rmatrix%NEQ, &
              p_KcolDest, p_KldDest, p_DaDest)

          ! Release the rest.
          call storage_free (h_Kcol)
          call storage_free (h_Kld)
        end if

      case (LSYSSC_MATRIXD)
        call output_line('Not implemented!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
        call sys_halt()

      end select

    case (LSYSSC_TR_ALL)
      ! Transpose the full matrix - structure+content
      select case (rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX9)
        ! Is the source matrix saved tranposed? In that case simply copy
        ! the data arrays, this results in the transposed matrix :)
        if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then

          if (rmatrix%h_Da .ne. ST_NOHANDLE) then
            call lsyssc_auxcopy_da (rmatrix,rtransposedMatrix)
          end if
          call lsyssc_auxcopy_Kcol (rmatrix,rtransposedMatrix)
          call lsyssc_auxcopy_Kld (rmatrix,rtransposedMatrix)
          call lsyssc_auxcopy_Kdiagonal (rmatrix,rtransposedMatrix)
          rtransposedMatrix%imatrixSpec = &
            iand(rtransposedMatrix%imatrixSpec,not(LSYSSC_MSPEC_TRANSPOSED))

          if (present(Ipermutation)) then
            call output_line('Calculation of the permutation not possible!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
            call sys_halt()
          end if

        else
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          if (rtransposedMatrix%h_Kcol .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                ST_INT, rtransposedMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
          end if

          if (rtransposedMatrix%h_Kld .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1, &
                ST_INT, rtransposedMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
          end if

          if (rtransposedMatrix%h_Da .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Da', rmatrix%NA, &
                rtransposedMatrix%cdataType, rtransposedMatrix%h_Da,&
                ST_NEWBLOCK_NOINIT)
          end if

          ! Get Kcol/Kld. Do not use the lsyssc_getbase routines as we just created
          ! the arrays for a matrix which is not yet valid!
          call storage_getbase_int(rmatrix%h_Kcol,p_KcolSource)
          call storage_getbase_int(rmatrix%h_Kld,p_KldSource)

          call storage_getbase_int(rtransposedMatrix%h_Kcol,p_KcolDest)
          call storage_getbase_int(rtransposedMatrix%h_Kld,p_KldDest)

          ! Get the data array(s)
          call lsyssc_getbase_double(rMatrix,p_DaSource)
          call storage_getbase_double(rtransposedMatrix%h_Da,p_DaDest)

          ! We need a temporary array
          call storage_new ('lsyssc_transposeMatrix','Itemp',rmatrix%NCOLS,&
                            ST_INT, h_Itemp, ST_NEWBLOCK_NOINIT)
          call storage_getbase_int (h_Itemp,p_Itemp)

          ! Calculate the transposed matrix structure:
          call lsyssc_transpMat79double (rmatrix%NEQ, rmatrix%NCOLS, &
                p_DaSource, p_KcolSource, p_KldSource, &
                p_Itemp, p_DaDest, p_KcolDest, p_KldDest, Ipermutation)

          ! Release the temp array
          call storage_free (h_Itemp)

          ! (Re-)calculate Kdiagonal
          call storage_new ('lsyssc_transposeMatrix', 'Kdiagonal', rmatrix%NCOLS, &
                            ST_INT, rtransposedMatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
          call storage_getbase_int(rtransposedMatrix%h_Kdiagonal,p_Kdiagonal)
          call lsyssc_rebuildKdiagonal (p_KcolDest, p_KldDest, p_Kdiagonal, &
                                        rmatrix%NCOLS)

        end if

      case (LSYSSC_MATRIX7)
        ! Is the source matrix saved tranposed? In that case simply copy
        ! the data arrays, this results in the transposed matrix :)
        if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then

          if (rmatrix%h_Da .ne. ST_NOHANDLE) then
            call lsyssc_auxcopy_da (rmatrix,rtransposedMatrix)
          end if

          call lsyssc_auxcopy_Kcol (rmatrix,rtransposedMatrix)
          call lsyssc_auxcopy_Kld (rmatrix,rtransposedMatrix)
          rtransposedMatrix%imatrixSpec = &
            iand(rtransposedMatrix%imatrixSpec,not(LSYSSC_MSPEC_TRANSPOSED))

          if (present(Ipermutation)) then
            call output_line('Calculation of the permutation not possible!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
            call sys_halt()
          end if

        else
          ! We really have to do some work now :)
          ! Allocate a new KCol and a new Kld.
          ! The Kld must be of size NCOLS+1, not NEQ+1, since we create a transposed
          ! matrix!
          if (rtransposedMatrix%h_Kcol .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Kcol', rmatrix%NA, &
                ST_INT, rtransposedMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
          end if

          if (rtransposedMatrix%h_Kld .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Kld', rmatrix%NCOLS+1, &
                ST_INT, rtransposedMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
          end if

          if (rtransposedMatrix%h_Da .eq. ST_NOHANDLE) then
            call storage_new ('lsyssc_transposeMatrix', 'Da', rmatrix%NA, &
                rtransposedMatrix%cdataType, rtransposedMatrix%h_Da,&
                ST_NEWBLOCK_NOINIT)
          end if

          ! Get Kcol/Kld. Do not use the lsyssc_getbase routines as we just created
          ! the arrays for a matrix which is not yet valid!
          call storage_getbase_int(rmatrix%h_Kcol,p_KcolSource)
          call storage_getbase_int(rmatrix%h_Kld,p_KldSource)

          call storage_getbase_int(rtransposedMatrix%h_Kcol,p_KcolDest)
          call storage_getbase_int(rtransposedMatrix%h_Kld,p_KldDest)

          ! Get the data array(s)
          call lsyssc_getbase_double(rMatrix,p_DaSource)
          call storage_getbase_double(rtransposedMatrix%h_Da,p_DaDest)

          ! We need a temporary array
          call storage_new ('lsyssc_transposeMatrix','Itemp',rmatrix%NCOLS,&
                            ST_INT, h_Itemp, ST_NEWBLOCK_NOINIT)
          call storage_getbase_int (h_Itemp,p_Itemp)

          ! Calculate the transposed matrix structure:
          call lsyssc_transpMat79double (rmatrix%NEQ, rmatrix%NCOLS, &
                p_DaSource, p_KcolSource, p_KldSource, &
                p_Itemp, p_DaDest, p_KcolDest, p_KldDest, Ipermutation)

          ! Release the temp array
          call storage_free (h_Itemp)

          ! Pivotise the matrix to move the diagonal element to the front.
          call lsyssc_pivotiseMatrix7double (rmatrix%NEQ, &
              p_KcolDest, p_KldDest, p_DaDest)
        end if

      case (LSYSSC_MATRIXD)
        call output_line('Not implemented!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
        call sys_halt()

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_transposeMatrix')
        call sys_halt()
      end select

      ! Exchange NEQ and NCOLS as these always describe the actual matrix.
      rtransposedMatrix%NEQ = rmatrix%NCOLS
      rtransposedMatrix%NCOLS = rmatrix%NEQ

      ! Switch the discretisation structures
      rtransposedMatrix%p_rspatialDiscrTrial => rmatrix%p_rspatialDiscrTest
      rtransposedMatrix%p_rspatialDiscrTest => rmatrix%p_rspatialDiscrTrial

    case default ! = LSYSSC_TR_ALL

    end select

    ! That is it.

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_auxcopy_Da (rsourceMatrix,rdestMatrix)

!<description>
  ! Auxiliary routine:
  ! Copies the content of rsourceMatrix%Da to rdestMatrix%Da
  ! without additional checks.
  ! If the destination array does not exist, it is created.
!</description>

!<input>
  ! Source matrix
  type(t_matrixScalar), intent(in) :: rsourceMatrix
!</input>

!<inputoutput>
  ! Destination matrix
  type(t_matrixScalar), intent(inout) :: rdestMatrix
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Da1,p_Da2
    real(SP), dimension(:), pointer :: p_Fa1,p_Fa2

    if (rsourceMatrix%h_Da .eq. rdestMatrix%h_Da) then
      ! Eehm... forget it.
      return
    end if

    if (rsourceMatrix%h_Da .eq. ST_NOHANDLE) then
      call output_line('Source matrix undefined!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_auxcopy_Da')
      call sys_halt()
    end if

    if (rdestMatrix%h_Da .eq. ST_NOHANDLE) then
      call storage_new ('lsyssc_auxcopy_Da', 'Da', rsourceMatrix%NA, &
          rsourceMatrix%cdataType, rdestMatrix%h_Da,ST_NEWBLOCK_NOINIT)
    end if

    select case (rsourceMatrix%cdataType)
    case (ST_DOUBLE)
      call lsyssc_getbase_double(rsourceMatrix,p_Da1)
      call lsyssc_getbase_double(rdestMatrix,p_Da2)
      call lalg_copyVectorDble (p_Da1,p_Da2)

    case (ST_SINGLE)
      call lsyssc_getbase_single(rsourceMatrix,p_Fa1)
      call lsyssc_getbase_single(rdestMatrix,p_Fa2)
      call lalg_copyVectorSngl (p_Fa1,p_Fa2)

    case default
      call output_line('Unsupported data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_auxcopy_Da')
      call sys_halt()

    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_auxcopy_Kcol (rsourceMatrix,rdestMatrix)

!<description>
  ! Auxiliary routine:
  ! Copies the content of rsourceMatrix%Kcol to rdestMatrix%Kcol
  ! without additional checks.
  ! If the destination array does not exist, it is created.
!</description>

!<input>
  ! Source matrix
  type(t_matrixScalar), intent(in) :: rsourceMatrix
!</input>

!<inputoutput>
  ! Destination matrix
  type(t_matrixScalar), intent(inout) :: rdestMatrix
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Kcol1,p_Kcol2

    if (rsourceMatrix%h_Kcol .eq. rdestMatrix%h_Kcol) then
      ! Eehm... forget it.
      return
    end if

    call lsyssc_getbase_Kcol (rsourceMatrix,p_Kcol1)
    if (rdestMatrix%h_Kcol .eq. ST_NOHANDLE) then
      call storage_new ('lsyssc_auxcopy_Kcol', 'Kcol', rsourceMatrix%NA, &
          ST_INT, rdestMatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
    end if
    call lsyssc_getbase_Kcol (rdestMatrix,p_Kcol2)

    call lalg_copyVectorInt (p_Kcol1,p_Kcol2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_auxcopy_Kld (rsourceMatrix,rdestMatrix)

!<description>
  ! Auxiliary routine:
  ! Copies the content of rsourceMatrix%Kld to rdestMatrix%Kld
  ! without additional checks.
  ! If the destination array does not exist, it is created.
!</description>

!<input>
  ! Source matrix
  type(t_matrixScalar), intent(in) :: rsourceMatrix
!</input>

!<inputoutput>
  ! Destination matrix
  type(t_matrixScalar), intent(inout) :: rdestMatrix
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Kld1,p_Kld2
    integer :: NEQ

    if (rsourceMatrix%h_Kld .eq. rdestMatrix%h_Kld) then
      ! Eehm... forget it.
      return
    end if

    call lsyssc_getbase_Kld (rsourceMatrix,p_Kld1)
    if (rdestMatrix%h_Kld .eq. ST_NOHANDLE) then

      ! Get the correct array length
      if (iand(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
        NEQ = rsourceMatrix%NCOLS
      else
        NEQ = rsourceMatrix%NEQ
      end if

      call storage_new ('lsyssc_auxcopy_Kld', 'Kld', &
          NEQ+1, &
          ST_INT, rdestMatrix%h_Kld,ST_NEWBLOCK_NOINIT)
    end if
    call lsyssc_getbase_Kld (rdestMatrix,p_Kld2)

    call lalg_copyVectorInt (p_Kld1,p_Kld2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_auxcopy_KrowIdx (rsourceMatrix,rdestMatrix)

!<description>
  ! Auxiliary routine:
  ! Copies the content of rsourceMatrix%KnzRows to rdestMatrix%KnzRows
  ! without additional checks.
  ! If the destination array does not exist, it is created.
!</description>

!<input>
  ! Source matrix
  type(t_matrixScalar), intent(in) :: rsourceMatrix
!</input>

!<inputoutput>
  ! Destination matrix
  type(t_matrixScalar), intent(inout) :: rdestMatrix
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Knz1,p_Knz2
    integer :: NEQ

    if (rsourceMatrix%h_KrowIdx .eq. rdestMatrix%h_KrowIdx) then
      ! Eehm... forget it.
      return
    end if

    call lsyssc_getbase_KrowIdx (rsourceMatrix,p_Knz1)
    if (rdestMatrix%h_Kld .eq. ST_NOHANDLE) then

      NEQ = rsourceMatrix%NNZROWS + rsourceMatrix%NEQ

      call storage_new ('lsyssc_auxcopy_KrowIdx', 'KrowIdx', &
          NEQ, &
          ST_INT, rdestMatrix%h_KrowIdx,ST_NEWBLOCK_NOINIT)
    end if
    call lsyssc_getbase_KrowIdx (rdestMatrix,p_Knz2)

    call lalg_copyVectorInt (p_Knz1,p_Knz2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_auxcopy_Kdiagonal (rsourceMatrix,rdestMatrix)

!<description>
  ! Auxiliary routine:
  ! Copies the content of rsourceMatrix%Kdiagonal to rdestMatrix%Kdiagonal
  ! without additional checks.
  ! If the destination array does not exist, it is created.
!</description>

!<input>
  ! Source matrix
  type(t_matrixScalar), intent(in) :: rsourceMatrix
!</input>

!<inputoutput>
  ! Destination matrix
  type(t_matrixScalar), intent(inout) :: rdestMatrix
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Kdiag1,p_Kdiag2
    integer :: NEQ

    if (rsourceMatrix%h_Kdiagonal .eq. rdestMatrix%h_Kdiagonal) then
      ! Eehm... forget it.
      return
    end if

    call lsyssc_getbase_Kdiagonal (rsourceMatrix,p_Kdiag1)
    if (rdestMatrix%h_Kdiagonal .eq. ST_NOHANDLE) then
      ! Get the correct array length
      if (iand(rdestMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
        NEQ = rsourceMatrix%NCOLS
      else
        NEQ = rsourceMatrix%NEQ
      end if

      call storage_new ('lsyssc_auxcopy_Kdiagonal', 'Kdiagonal', &
          NEQ, ST_INT, rdestMatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
    end if
    call lsyssc_getbase_Kdiagonal (rdestMatrix,p_Kdiag2)

    call lalg_copyVectorInt (p_Kdiag1,p_Kdiag2)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_allocEmptyMatrix (rmatrixScalar,iclear,bignoreExisting,cdataType)

!<description>
  ! This routine allocates memory for the matrix entries without computing
  ! the entries. This can be used to attach an 'empty' matrix to a matrix
  ! structure. The number of entries NA as well as the type of the matrix
  ! cmatrixFormat must be initialised in rmatrixScalar.
  !
  ! Memory is allocated if either the matrix has no content attached or
  ! if the content belongs to another matrix. In both case, a new matrix
  ! content array is created and the matrix will be the owner of that
  ! memory block.
!</description>

!<input>
  ! Whether and how to fill the matrix with initial values.
  ! One of the LSYSSC_SETM_xxxx constants:
  ! LSYSSC_SETM_UNDEFINED : Do not initialise the matrix,
  ! LSYSSC_SETM_ZERO      : Clear the matrix / fill it with 0.0,
  ! LSYSSC_SETM_ONE       : Fill the matrix with 1.0. (Used e.g.
  !                         for UMFPACK who needs a non-zero
  !                         matrix for symbolic factorisation.)
  integer, intent(in) :: iclear

  ! OPTIONAL: If set to TRUE, existing matrices are ignored.
  ! Standard value is FALSE which stops the application with an error.
  logical, intent(in), optional :: bignoreExisting

  ! OPTIONAL: Data type of the matrix (ST_SINGLE, ST_DOUBLE)
  ! If not present, the standard data type cdataType-variable in the matrix
  ! is used (usually ST_DOUBLE).
  integer, intent(in), optional :: cdataType
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(inout) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: NA
  integer :: cdType,NVAR
  real(SP), dimension(:), pointer :: p_Fa
  real(DP), dimension(:), pointer :: p_Da

  if ((rmatrixScalar%h_Da .ne. ST_NOHANDLE) .and. &
      (iand(rmatrixScalar%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .eq. 0)) then
    if (present(bignoreExisting)) then
      if (bignoreExisting) return
    end if
    call output_line('Cannot create empty matrix; exists already!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_allocEmptyMatrix')
    call sys_halt()
  end if

  NA = rmatrixScalar%NA
  NVAR = rmatrixScalar%NVAR

  if (present(cdataType)) then
    cdType = cdataType
    ! Change the data type in the matrix according to the new we will use.
    rmatrixScalar%cdataType = cdType
  else
    cdType = rmatrixScalar%cdataType
  end if

  ! Which matrix structure do we have?
  select case (rmatrixScalar%cmatrixFormat)
  case (LSYSSC_MATRIX9,LSYSSC_MATRIX7,LSYSSC_MATRIXD,LSYSSC_MATRIX1,LSYSSC_MATRIX9ROWC)

    ! Check if the matrix entries exist and belongs to us.
    ! If not, allocate the matrix.
    if ((rmatrixScalar%h_Da .eq. ST_NOHANDLE) .or. &
        (iand(rmatrixScalar%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .ne. 0)) then

      if (iclear .ge. LSYSSC_SETM_ZERO) then
        call storage_new ('lsyssc_allocEmptyMatrix', 'Da', &
                            NA, cdType, rmatrixScalar%h_Da, &
                            ST_NEWBLOCK_ZERO)
      else
        call storage_new ('lsyssc_allocEmptyMatrix', 'Da', &
                            NA, cdType, rmatrixScalar%h_Da, &
                            ST_NEWBLOCK_NOINIT)
      end if

      if (iclear .ge. LSYSSC_SETM_ONE) then
        select case (cdType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double (rmatrixScalar,p_Da)
          call lalg_setVectorDble (p_Da,1.0_DP)
        case (ST_SINGLE)
          call lsyssc_getbase_single (rmatrixScalar,p_Fa)
          call lalg_setVectorSngl (p_Fa,1.0_SP)
        case default
          call output_line('Unknown data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_allocEmptyMatrix')
          call sys_halt()
        end select
      end if

    end if

  case (LSYSSC_MATRIX9INTL,LSYSSC_MATRIX7INTL)

    ! Check if the matrix entries exist. If not, allocate the matrix.
    if ((rmatrixScalar%h_Da .eq. ST_NOHANDLE) .or. &
        (iand(rmatrixScalar%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .ne. 0)) then

      if (iclear .ge. LSYSSC_SETM_ZERO) then
        select case (rmatrixScalar%cinterleavematrixFormat)
        case (LSYSSC_MATRIX1)
          call storage_new ('lsyssc_allocEmptyMatrix', 'Da', &
              NA*NVAR*NVAR, cdType, rmatrixScalar%h_Da, ST_NEWBLOCK_ZERO)
        case (LSYSSC_MATRIXD)
          call storage_new ('lsyssc_allocEmptyMatrix', 'Da', &
              NA*NVAR, cdType, rmatrixScalar%h_Da, ST_NEWBLOCK_ZERO)
        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_allocEmptyMatrix')
          call sys_halt()
        end select

      else
        select case (rmatrixScalar%cinterleavematrixFormat)
        case (LSYSSC_MATRIX1)
          call storage_new ('lsyssc_allocEmptyMatrix', 'Da', &
              NA*NVAR*NVAR, cdType, rmatrixScalar%h_Da, ST_NEWBLOCK_NOINIT)
        case (LSYSSC_MATRIXD)
          call storage_new ('lsyssc_allocEmptyMatrix', 'Da', &
              NA*NVAR, cdType, rmatrixScalar%h_Da, ST_NEWBLOCK_NOINIT)
        case default
          call output_line('Unsupported interleave matrix format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_allocEmptyMatrix')
          call sys_halt()
        end select

      end if

      if (iclear .ge. LSYSSC_SETM_ONE) then
        select case (cdType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double (rmatrixScalar,p_Da)
          call lalg_setVectorDble (p_Da,1.0_DP)
        case (ST_SINGLE)
          call lsyssc_getbase_single (rmatrixScalar,p_Fa)
          call lalg_setVectorSngl (p_Fa,1.0_SP)
        case default
          call output_line('Unknown data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_allocEmptyMatrix')
          call sys_halt()
        end select
      end if

    end if

  case default
    call output_line('Not supported matrix structure!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_allocEmptyMatrix')
    call sys_halt()
  end select

  ! The matrix is now the owner of the memory.
  rmatrixScalar%imatrixSpec = &
      iand(rmatrixScalar%imatrixSpec,not(LSYSSC_MSPEC_CONTENTISCOPY))

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_createDiagMatrixStruc (rmatrix,NEQ,cmatrixFormat)

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
  integer, intent(in) :: cmatrixFormat

  ! Number of equations, the matrix should hold.
  integer, intent(in) :: NEQ
!</input>

!<inputoutput>
  ! The FE matrix to be initialised.
  type(t_matrixScalar), intent(out) :: rmatrix
!</inputoutput>

!</subroutine>

    rmatrix%NEQ = NEQ
    rmatrix%NCOLS = NEQ
    rmatrix%NA = NEQ
    rmatrix%cmatrixFormat = cmatrixFormat

    ! Depending on the format, choose the initialisation routine.
    select case (cmatrixFormat)
    case (LSYSSC_MATRIXD,LSYSSC_MATRIX1)
      ! Nothing to do, there is no structure.

    case (LSYSSC_MATRIX9)
      ! Create KCOL, KLD, KDiagonal
      call storage_new ('lsyssc_createDiagMatrixStuc', 'KCOL', NEQ, &
          ST_INT, rmatrix%h_Kcol, ST_NEWBLOCK_ORDERED)
      call storage_new ('lsyssc_createDiagMatrixStruc', 'KLD', NEQ+1, &
          ST_INT, rmatrix%h_Kld, ST_NEWBLOCK_ORDERED)
      call storage_new ('lsyssc_createDiagMatrixStruc', 'KDiagonal', NEQ, &
          ST_INT, rmatrix%h_Kdiagonal, ST_NEWBLOCK_ORDERED)

    case (LSYSSC_MATRIX7)
      ! Create KCOL, KLD.
      call storage_new ('lsyssc_createDiagMatrixStruc', 'KCOL', NEQ, &
          ST_INT, rmatrix%h_Kcol, ST_NEWBLOCK_ORDERED)
      call storage_new ('lsyssc_createDiagMatrixStruc', 'KLD', NEQ+1, &
          ST_INT, rmatrix%h_Kld, ST_NEWBLOCK_ORDERED)

    case default

      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_createDiagMatrixStruc')
      call sys_halt()

    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_createFullMatrix (rmatrix,bclear,NEQ,ncols,cdataType)

!<description>
  ! Creates a simple NEQ*ncols matrix in matrix structure 1.
  ! No discretisation structure or similar will be attached to that matrix.
!</description>

!<input>
  ! Number of rows in the matrix
  integer, intent(in) :: NEQ

  ! Whether to initialise the matrix with zero or not.
  logical, intent(in) :: bclear

  ! OPTIONAL: Number of columns in the matrix. If not specified, NEQ=nrows
  ! is assumed.
  integer, intent(in), optional :: ncols

  ! OPTIONAL: Data type of the entries. A ST_xxxx constant.
  ! Default is ST_DOUBLE.
  integer, intent(in), optional :: cdataType
!</input>

!<inputoutput>
  ! The matrix to be initialised.
  type(t_matrixScalar), intent(out) :: rmatrix
!</inputoutput>

!</subroutine>

    rmatrix%NEQ = NEQ
    rmatrix%NCOLS = NEQ
    if (present(ncols)) rmatrix%NCOLS = ncols
    rmatrix%NA = rmatrix%NEQ * rmatrix%NCOLS
    rmatrix%cmatrixFormat = LSYSSC_MATRIX1
    if (present(cdataType)) rmatrix%cdataType = cdataType

    if (bclear) then
      call storage_new ('lsyssc_createFullMatrix', 'h_Da', &
          rmatrix%NA , rmatrix%cdataType, rmatrix%h_Da,ST_NEWBLOCK_ZERO)
    else
      call storage_new ('lsyssc_createFullMatrix', 'h_Da', &
          rmatrix%NA , rmatrix%cdataType, rmatrix%h_Da,ST_NEWBLOCK_NOINIT)
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_lumpMatrixScalar (rmatrixScalar,clumpType,bconvertToDiag)

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
  integer, intent(in) :: clumpType

  ! OPTIONAL: Whether to convert the matrix structure to a diagonal matrix.
  ! TRUE:  Convert the matrix to a diagonal matrix. This is the default case
  !        if the parameter is not present.
  ! FALSE: Matrix structure stays unchanged.
  logical, intent(in), optional :: bconvertToDiag
!</input>

!<inputoutput>
  ! The matrix to be lumped. Will be overwritten by a diagonal matrix.
  type(t_matrixScalar), intent(inout) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: irow,icol,j
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kdiagonal, p_Kld
  real(DP), dimension(:), pointer :: p_Da
  real(SP), dimension(:), pointer :: p_Fa
  logical :: bconvert

  if (rmatrixScalar%NEQ .eq. 0) return

  bconvert = .true.
  if (present(bconvertToDiag)) bconvert = bconvertToDiag

  ! Which type of lumping do we have to do?

  select case (clumpType)
  case (LSYSSC_LUMP_STD)
    ! Standard lumping. Nothing to do here.

  case (LSYSSC_LUMP_DIAG)
    ! Diagonal lumping. We have to add all off-diagonal entries to the main
    ! diagonal.
    !
    ! What is the matrix format?
    select case (rmatrixScalar%cmatrixFormat)

    case (LSYSSC_MATRIX7)
      ! Get the structure
      call lsyssc_getbase_Kld (rmatrixScalar,p_Kld)
      call lsyssc_getbase_Kcol (rmatrixScalar,p_Kcol)

      ! Get the data and perform the lumping
      select case (rmatrixScalar%cdataType)

      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrixScalar,p_Da)

        ! Add all off-diagonals to the diagonal
        do irow = 1,rmatrixScalar%NEQ
          j = p_Kld(irow)
          do icol = p_Kld(irow)+1,p_Kld(irow+1)-1
            p_Da(j) = p_Da(j) + p_Da(icol)
          end do
        end do

      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixScalar,p_Fa)

        ! Add all off-diagonals to the diagonal
        do irow = 1,rmatrixScalar%NEQ
          j = p_Kld(irow)
          do icol = p_Kld(irow)+1,p_Kld(irow+1)-1
            p_Fa(j) = p_Fa(j) + p_Fa(icol)
          end do
        end do

      case default
        call output_line('Unsupported matrix precision!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_lumpMatrixScalar')
        call sys_halt()
      end select

    case (LSYSSC_MATRIX9)
      ! Get the structure
      call lsyssc_getbase_Kld (rmatrixScalar,p_Kld)
      call lsyssc_getbase_Kdiagonal (rmatrixScalar,p_Kdiagonal)
      call lsyssc_getbase_Kcol (rmatrixScalar,p_Kcol)

      ! Get the data and perform the lumping
      select case (rmatrixScalar%cdataType)

      case (ST_DOUBLE)
        call lsyssc_getbase_double (rmatrixScalar,p_Da)

        ! Add all off-diagonals to the diagonal
        do irow = 1,rmatrixScalar%NEQ
          j = p_Kdiagonal(irow)
          do icol = p_Kld(irow),j-1
            p_Da(j) = p_Da(j) + p_Da(icol)
          end do
          do icol = j+1,p_Kld(irow+1)-1
            p_Da(j) = p_Da(j) + p_Da(icol)
          end do
        end do

      case (ST_SINGLE)
        call lsyssc_getbase_single (rmatrixScalar,p_Fa)

        ! Add all off-diagonals to the diagonal
        do irow = 1,rmatrixScalar%NEQ
          j = p_Kdiagonal(irow)
          do icol = p_Kld(irow),j-1
            p_Fa(j) = p_Fa(j) + p_Fa(icol)
          end do
          do icol = j+1,p_Kld(irow+1)-1
            p_Fa(j) = p_Fa(j) + p_Fa(icol)
          end do
        end do

      case default
        call output_line('Unsupported matrix precision!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_lumpMatrixScalar')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_lumpMatrixScalar')
      call sys_halt()
    end select

  end select

  if (bconvert) then
    ! Convert the given matrix into a diagonal matrix by extracting
    ! the diagonal entries.
    call lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIXD)
  else
    ! Clear all off-diagonal entries, so the matrix will be a diagonal
    ! matrix in the original structure.
    call lsyssc_clearOffdiags (rmatrixScalar)
  end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_clearOffdiags (rmatrix)

!<description>
  ! Sets all off-diagonal entries in the matrix rmatrix to zero.
!</description>

!<inputoutput>
  ! The matrix which is to be modified.
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

  ! At first we must take care of the matrix type.
  select case (rmatrix%cmatrixFormat)
  case (LSYSSC_MATRIX9)
    call removeOffdiags_format9 (rmatrix)
  case (LSYSSC_MATRIX7)
    call removeOffdiags_format7 (rmatrix)
  case (LSYSSC_MATRIXD)
    ! Nothing to do
  case (LSYSSC_MATRIX1)
    call removeOffdiags_format1 (rmatrix)
  case default
    call output_line('Unsupported matrix format!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_clearOffdiags')
    call sys_halt()
  end select

  contains

    ! ****************************************
    ! The replacement routine for format 9

    subroutine removeOffdiags_format9 (rmatrix)

    type(t_matrixScalar), intent(inout) :: rmatrix

    ! local variables
    integer :: irow
    real(DP), dimension(:), pointer :: p_Da
    real(SP), dimension(:), pointer :: p_FA
    integer, dimension(:), pointer :: p_Kld,p_Kdiagonal
    real(DP) :: ddiag,fdiag

    ! Get Kld and Kdiagonal
    call lsyssc_getbase_Kld(rmatrix,p_Kld)
    call lsyssc_getbase_Kdiagonal(rmatrix,p_Kdiagonal)

    ! Take care of the format of the entries
    select case (rmatrix%cdataType)
    case (ST_DOUBLE)
      ! Get the data array
      call lsyssc_getbase_double(rmatrix,p_Da)

      ! loop through the rows
      do irow = 1,rmatrix%NEQ

        ! Get the diagonal
        ddiag = p_Da(p_Kdiagonal(irow))

        ! Clear the row
        p_Da(p_Kld(irow):p_Kld(irow+1)-1) = 0.0_DP

        ! restore the diagonal
        p_Da(p_Kdiagonal(irow)) = ddiag

      end do

    case (ST_SINGLE)
      ! Get the data array
      call lsyssc_getbase_single(rmatrix,p_FA)

      ! loop through the rows
      do irow = 1,rmatrix%NEQ

        ! Get the diagonal
        fdiag = p_FA(p_Kdiagonal(irow))

        ! Clear the row
        p_FA(p_Kld(irow):p_Kld(irow+1)-1) = 0.0_SP

        ! restore the diagonal
        p_FA(p_Kdiagonal(irow)) = Fdiag

      end do

    case default
      call output_line('Unsupported matrix precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_clearOffdiags')
      call sys_halt()

    end select

    end subroutine

    ! ****************************************
    ! The replacement routine for format 7

    subroutine removeOffdiags_format7 (rmatrix)

    type(t_matrixScalar), intent(inout) :: rmatrix

    ! local variables
    integer :: irow
    integer, dimension(:), pointer :: p_Kld
    real(DP), dimension(:), pointer :: p_Da
    real(SP), dimension(:), pointer :: p_FA

    ! Get Kld:
    call lsyssc_getbase_Kld(rmatrix,p_Kld)

    ! Take care of the format of the entries
    select case (rmatrix%cdataType)
    case (ST_DOUBLE)
      ! Get the data array
      call lsyssc_getbase_double(rmatrix,p_Da)

      ! loop through the rows
      do irow = 1,rmatrix%NEQ

        ! Clear the row except for the diagonal
        p_Da(p_Kld(irow)+1:p_Kld(irow+1)-1) = 0.0_DP

      end do

    case (ST_SINGLE)
      ! Get the data array
      call lsyssc_getbase_single(rmatrix,p_FA)

      ! loop through the rows
      do irow = 1,rmatrix%NEQ

        ! Clear the row except for the diagonal
        p_FA(p_Kld(irow)+1:p_Kld(irow+1)-1) = 0.0_SP

      end do

    case default
      call output_line('Unsupported matrix precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_clearOffdiags')
      call sys_halt()
    end select

    end subroutine

    ! ****************************************
    ! The replacement routine for format 1

    subroutine removeOffdiags_format1 (rmatrix)

    type(t_matrixScalar), intent(inout) :: rmatrix

    ! local variables
    integer :: icol
    real(DP), dimension(:), pointer :: p_Da
    real(SP), dimension(:), pointer :: p_FA
    real(DP) :: ddata
    real(SP) :: fdata

    ! Take care of the format of the entries
    select case (rmatrix%cdataType)
    case (ST_DOUBLE)
      ! Get the data array
      call lsyssc_getbase_double(rmatrix,p_Da)

      ! loop through the rows ald colums, set off-diagonal entries to zero.
      do icol = 1,rmatrix%NCOLS
        if (icol .le. rmatrix%NEQ) then
          ! Remember the diagonal
          ddata = p_Da((icol-1)*rmatrix%NEQ+icol)

          ! Fill the column with zero
          p_Da ((icol-1)*rmatrix%NEQ+1 : icol*rmatrix%NEQ) = 0.0_DP

          ! Restore the diagonal
          p_Da((icol-1)*rmatrix%NEQ+icol) = ddata
        else
          ! Fill the column with zero
          p_Da ((icol-1)*rmatrix%NEQ+1 : icol*rmatrix%NEQ) = 0.0_DP
        end if
      end do

    case (ST_SINGLE)
      ! Get the data array
      call lsyssc_getbase_single(rmatrix,p_FA)

      ! loop through the rows ald colums, set off-diagonal entries to zero.
      do icol = 1,rmatrix%NCOLS
        if (icol .le. rmatrix%NEQ) then
          ! Remember the diagonal
          fdata = p_Fa((icol-1)*rmatrix%NEQ+icol)

          ! Fill the column with zero
          p_Fa ((icol-1)*rmatrix%NEQ+1 : icol*rmatrix%NEQ) = 0.0_DP

          ! Restore the diagonal
          p_Fa((icol-1)*rmatrix%NEQ+icol) = fdata
        else
          ! Fill the column with zero
          p_Fa ((icol-1)*rmatrix%NEQ+1 : icol*rmatrix%NEQ) = 0.0_DP
        end if
      end do

    case default
      call output_line('Unsupported matrix precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_clearOffdiags')
      call sys_halt()
    end select

    end subroutine

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_scaleMatrix (rmatrix,c)

!<description>
  ! Scales a matrix rmatrix: rmatrix = c * rmatrix
!</description>

!<inputoutput>
  ! Source and destination matrix
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!<input>
  ! Multiplication factor
  real(DP), intent(in) :: c
!</input>

!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata

  ! Take care of the data type!
  select case (rmatrix%cdataType)
  case (ST_DOUBLE)
    ! Get the pointer and scale the whole data array.
    call lsyssc_getbase_double(rmatrix,p_Ddata)
    call lalg_scaleVector(p_Ddata,c)

  case (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    call lsyssc_getbase_single(rmatrix,p_Fdata)
    call lalg_scaleVector(p_Fdata,real(c,SP))

  case default
    call output_line('Unsupported data type!',&
        OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_scaleMatrix')
    call sys_halt()
  end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_multMatMat (rmatrixA,rmatrixB,rmatrixC,bmemory,&
      bsymb,bnumb,bisExactStructure)

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
    type(t_matrixScalar), intent(in) :: rmatrixA,rmatrixB

    ! BMEMORY = FALSE: Do not allocate required memory for C=A*B.
    ! BMEMORY = TRUE:  Generate all required structures for C=A*B
    logical :: bmemory

    ! Compute symbolic matrix-matrix-product
    ! BSYMB = FALSE: Do not generate the required matrix structures.
    !                This may be useful, if the same matrices A and B
    !                need to be multiplied several times, but the
    !                sparsity patterns do not change
    ! BSYMN = TRUE:  Generate all required matrix structures for C=A*B
    logical :: bsymb

    ! Compute numerical matrix-matrix-product
    ! BNUMB = FALSE: Do not perform the numerical multiplication
    ! BNUMB = TRUE:  Perform numerical multiplication
    logical :: bnumb

    ! OPTIONAL: Indicates whether the resulting matrix C has the
    ! required symbolic structure or is a superset of the symbolic
    ! matrix-matrix product. In some cases, this may allow for much
    ! more efficient implementation.
    logical, optional :: bisExactStructure
!</input>

!<inputoutput>
    type(t_matrixScalar), intent(inout) :: rmatrixC
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: DaA,DaB,DaC,Daux
    real(SP), dimension(:), pointer :: FaA,FaB,FaC,Faux
    integer, dimension(:), pointer :: KcolA,KcolB,KcolC,Kaux
    integer, dimension(:), pointer :: KldA,KldB,KldC,KdiagonalC
    integer :: h_Kaux,h_Daux,h_Faux
    logical :: bfast

    h_Kaux=ST_NOHANDLE
    h_Daux=ST_NOHANDLE
    h_Faux=ST_NOHANDLE

    ! Check if fast implementation can be used
    bfast=.false.
    if (present(bisExactStructure)) bfast = bisExactStructure
    bfast = bfast .or. (bmemory .and. bsymb)

    ! Check if both matrices are compatible
    if (rmatrixA%NCOLS .ne. rmatrixB%NEQ) then
      call output_line("Number of columns of matrix A is not compatible with "//&
          "number of rows of matrix B!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
      call sys_halt()
    end if

    ! Check if both matrices have the same sorting
    if (rmatrixA%bcolumnsSorted .or. rmatrixB%bcolumnsSorted) then
      if (rmatrixA%p_rsortStrategyColumns%ctype .ne. rmatrixB%p_rsortStrategyColumns%ctype) then
        call output_line("Incompatible sorting strategies!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
        call sys_halt()
      end if
    end if

    if (rmatrixA%browsSorted .or. rmatrixB%browsSorted) then
      if (rmatrixA%p_rsortStrategyRows%ctype .ne. rmatrixB%p_rsortStrategyRows%ctype) then
        call output_line("Incompatible sorting strategies!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
        call sys_halt()
      end if
    end if

    ! Release matrix if required and set common variables
    if (bmemory) then
      call lsyssc_releaseMatrix(rmatrixC)
      if (rmatrixA%cdataType .eq. ST_DOUBLE .or. rmatrixB&
          &%cdataType .eq. ST_DOUBLE) then
        rmatrixC%cdataType = ST_DOUBLE
      else
        rmatrixC%cdataType = ST_SINGLE
      end if
    end if

    ! Set sorting strategy for matrix C
    rmatrixC%bcolumnsSorted = rmatrixA%bcolumnsSorted
    rmatrixC%browsSorted = rmatrixA%browsSorted
    rmatrixC%p_rsortStrategyColumns => rmatrixA%p_rsortStrategyColumns
    rmatrixC%p_rsortStrategyRows => rmatrixA%p_rsortStrategyRows

    ! Perform matrix-matrix multiplication
    select case (rmatrixA%cmatrixFormat)

    case (LSYSSC_MATRIX1)   ! A is full matrix ----------------------

      select case (rmatrixB%cmatrixFormat)

      case (LSYSSC_MATRIX1) ! B is full matrix - - - - - - - - - - -

        ! memory allocation?
        if (bmemory) then
          rmatrixC%cmatrixFormat = LSYSSC_MATRIX1
          rmatrixC%NA = rmatrixA%NEQ*rmatrixB%NCOLS
          call storage_new("lsyssc_multMatMat","h_Da",rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (rmatrixC%cmatrixFormat .ne. LSYSSC_MATRIX1) then
          call output_line("Destination matrix has incompatible format!",&
              OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
          call sys_halt()
        end if

        ! symbolical multiplication?
        if (bsymb) then
          rmatrixC%NEQ = rmatrixA%NEQ
          rmatrixC%NCOLS = rmatrixB%NCOLS
        end if

        ! numerical multiplication?
        if (bnumb) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call do_mat1mat1mulDbleDble(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,rmatrixB%NCOLS,DaA,DaB,DaC)

            case (ST_SINGLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call do_mat1mat1mulDbleSngl(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,rmatrixB%NCOLS,DaA,FaB,DaC)

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case (ST_SINGLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call do_mat1mat1mulSnglDble(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,rmatrixB%NCOLS,FaA,DaB,DaC)

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_single(rmatrixC,FaC)
              call do_mat1mat1mulSnglSngl(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,rmatrixB%NCOLS,FaA,FaB,FaC)

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case default
            call output_line("Unsupported data type!",&
                OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
            call sys_halt()
          end select
        end if

      case (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - -

        ! memory allocation?
        if (bmemory) then
          rmatrixC%cmatrixFormat = LSYSSC_MATRIX1
          rmatrixC%NA = rmatrixA%NA
          call storage_new("lsyssc_multMatMat","h_Da",rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (rmatrixC%cmatrixFormat .ne. LSYSSC_MATRIX1) then
          call output_line("Destination matrix has incompatible format!",&
              OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
          call sys_halt()
        end if

        ! symbolical multiplication?
        if (bsymb) then
          rmatrixC%NEQ = rmatrixA%NEQ
          rmatrixC%NCOLS = rmatrixA%NCOLS
        end if

        ! numerical multiplication?
        if (bnumb) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call do_mat1matDmulDbleDble(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,DaA,DaB,DaC)

            case (ST_SINGLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call do_mat1matDmulDbleSngl(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,DaA,FaB,DaC)

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case (ST_SINGLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call do_mat1matDmulSnglDble(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,FaA,DaB,DaC)

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_single(rmatrixC,FaC)
              call do_mat1matDmulSnglSngl(rmatrixA%NEQ,rmatrixA&
                  &%NCOLS,FaA,FaB,FaC)

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case default
            call output_line("Unsupported data type!",&
                OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
            call sys_halt()
          end select
        end if

      case default
        call output_line("Unsupported data type!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
        call sys_halt()
      end select

      !--------------------------------------------------------------

    case (LSYSSC_MATRIXD)   ! A is diagonal matrix ------------------

      select case (rmatrixB%cmatrixFormat)

      case (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - -

        ! memory allocation?
        if (bmemory) then
          rmatrixC%cmatrixFormat = LSYSSC_MATRIXD
          rmatrixC%NA = rmatrixA%NA
          call storage_new("lsyssc_multMatMat","h_Da",rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (rmatrixC%cmatrixFormat .ne. LSYSSC_MATRIXD) then
          call output_line("Destination matrix has incompatible format!",&
              OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
          call sys_halt()
        end if

        ! symbolic multiplication?
        if (bsymb) then
          rmatrixC%NEQ = rmatrixA%NEQ
          rmatrixC%NCOLS = rmatrixA%NCOLS
        end if

        ! numerical multiplication?
        if (bnumb) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              DaC=DaA*DaB

            case (ST_SINGLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              DaC=DaA*FaB

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case (ST_SINGLE)

            select case (rmatrixC%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              DaC=FaA*DaB

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_single(rmatrixC,FaC)
              FaC=FaA*FaB

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case default
            call output_line("Unsupported data type!",&
                OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
            call sys_halt()
          end select
        end if

      case (LSYSSC_MATRIX1) ! B is full matrix - - - - - - - - - - -

        ! memory allocation?
        if (bmemory) then
          rmatrixC%cmatrixFormat = LSYSSC_MATRIX1
          rmatrixC%NA = rmatrixB%NA
          call storage_new("lsyssc_multMatMat","h_Da",rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (rmatrixC%cmatrixFormat .ne. LSYSSC_MATRIX1) then
          call output_line("Destination matrix has incompatible format!",&
              OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
          call sys_halt()
        end if

        ! symbolical multiplication?
        if (bsymb) then
          rmatrixC%NEQ = rmatrixB%NEQ
          rmatrixC%NCOLS = rmatrixB%NCOLS
        end if

        ! numerical multiplication?
        if (bnumb) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call do_matDmat1mulDbleDble(rmatrixB%NEQ,rmatrixB&
                  &%NCOLS,DaA,DaB,DaC)

            case (ST_SINGLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call do_matDmat1mulDbleSngl(rmatrixB%NEQ,rmatrixB&
                  &%NCOLS,DaA,FaB,DaC)

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case (ST_SINGLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call do_matDmat1mulSnglDble(rmatrixB%NEQ,rmatrixB&
                  &%NCOLS,FaA,DaB,DaC)

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_single(rmatrixC,FaC)
              call do_matDmat1mulSnglSngl(rmatrixB%NEQ,rmatrixB&
                  &%NCOLS,FaA,FaB,FaC)

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case default
            call output_line("Unsupported data type!",&
                OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
            call sys_halt()
          end select
        end if

      case (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! B is CSR matrix - - - -

        ! memory allocation?
        if (bmemory) then
          rmatrixC%cmatrixFormat = rmatrixB%cmatrixFormat
          rmatrixC%NA = rmatrixB%NA
          call storage_new("lsyssc_multMatMat","h_Da",rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (rmatrixC%cmatrixFormat .ne. rmatrixB%cmatrixFormat) then
          call output_line("Destination matrix has incompatible format!",&
              OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
          call sys_halt()
        end if

        ! symbolical multiplication?
        if (bsymb) then
          rmatrixC%NEQ = rmatrixB%NEQ
          rmatrixC%NCOLS = rmatrixB%NCOLS
          rmatrixC%h_Kld = rmatrixB%h_Kld
          rmatrixC%h_Kcol = rmatrixB%h_Kcol
          rmatrixC%h_Kdiagonal = rmatrixB%h_Kdiagonal
          rmatrixC%imatrixSpec = ior(rmatrixC%imatrixSpec,&
              &LSYSSC_MSPEC_STRUCTUREISCOPY)
        end if

        ! numerical multiplication?
        if (bnumb) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call lsyssc_getbase_Kld(rmatrixB,KldB)
              call lsyssc_getbase_Kcol(rmatrixB,KcolB)
              if (bfast) then
                call do_matDmat79mulDbleDble(DaA,DaB,KldB,KcolB&
                    &,rmatrixB%NEQ,DaC)
              else
                call lsyssc_getbase_Kld(rmatrixC,KldC)
                call lsyssc_getbase_Kcol(rmatrixC,KcolC)
                call do_matDmat79mulDbleDble(DaA,DaB,KldB,KcolB&
                    &,rmatrixB%NEQ,DaC,KldC,KcolC)
              end if

            case (ST_SINGLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call lsyssc_getbase_Kld(rmatrixB,KldB)
              call lsyssc_getbase_Kld(rmatrixB,KcolB)
              if (bfast) then
                call do_matDmat79mulDbleSngl(DaA,FaB,KldB,KcolB&
                    &,rmatrixB%NEQ,DaC)
              else
                call lsyssc_getbase_Kld(rmatrixC,KldC)
                call lsyssc_getbase_Kcol(rmatrixC,KcolC)
                call do_matDmat79mulDbleSngl(DaA,FaB,KldB,KcolB&
                    &,rmatrixB%NEQ,DaC,KldC,KcolC)
              end if

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case (ST_SINGLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call lsyssc_getbase_Kld(rmatrixB,KldB)
              call lsyssc_getbase_Kcol(rmatrixB,KcolB)
              if (bfast) then
                call do_matDmat79mulSnglDble(FaA,DaB,KldB,KcolB&
                    &,rmatrixB%NEQ,DaC)
              else
                call lsyssc_getbase_Kld(rmatrixC,KldC)
                call lsyssc_getbase_Kcol(rmatrixC,KcolC)
                call do_matDmat79mulSnglDble(FaA,DaB,KldB,KcolB&
                    &,rmatrixB%NEQ,DaC,KldC,KcolC)
              end if

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_single(rmatrixC,FaC)
              call lsyssc_getbase_Kld(rmatrixB,KldB)
              call lsyssc_getbase_Kcol(rmatrixB,KcolB)
              if (bfast) then
                call do_matDmat79mulSnglSngl(FaA,FaB,KldB,KcolB&
                    &,rmatrixB%NEQ,FaC)
              else
                call lsyssc_getbase_Kld(rmatrixC,KldC)
                call lsyssc_getbase_Kcol(rmatrixC,KcolC)
                call do_matDmat79mulSnglSngl(FaA,FaB,KldB,KcolB&
                    &,rmatrixB%NEQ,FaC,KldC,KcolC)
              end if

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case default
            call output_line("Unsupported data type!",&
                OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
            call sys_halt()
          end select
        end if

      case default
        call output_line("Unsupported data type!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
        call sys_halt()
      end select

      !--------------------------------------------------------------

    case (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! A is CSR matrix ----------

      select case (rmatrixB%cmatrixFormat)

      case (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - -

        ! memory allocation?
        if (bmemory) then
          rmatrixC%cmatrixFormat = rmatrixA%cmatrixFormat
          rmatrixC%NA = rmatrixA%NA
          call storage_new("lsyssc_multMatMat","h_Da",rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (rmatrixC%cmatrixFormat .ne. rmatrixA%cmatrixFormat) then
          call output_line("Destination matrix has incompatible format!",&
              OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
          call sys_halt()
        end if

        ! symbolical multiplication?
        if (bsymb) then
          rmatrixC%NEQ = rmatrixA%NEQ
          rmatrixC%NCOLS = rmatrixA%NCOLS
          rmatrixC%h_Kld = rmatrixA%h_Kld
          rmatrixC%h_Kcol = rmatrixA%h_Kcol
          rmatrixC%h_Kdiagonal = rmatrixA%h_Kdiagonal
          rmatrixC%imatrixSpec = ior(rmatrixC%imatrixSpec,&
              &LSYSSC_MSPEC_STRUCTUREISCOPY)
        end if

        ! numerical multiplication?
        if (bnumb) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call lsyssc_getbase_Kld(rmatrixA,KldA)
              call lsyssc_getbase_Kcol(rmatrixA,KcolA)
              if (bfast) then
                call do_mat79matDmulDbleDble(DaA,KldA,KcolA&
                    &,rmatrixA%NEQ,DaB,DaC)
              else
                call lsyssc_getbase_Kld(rmatrixC,KldC)
                call lsyssc_getbase_Kcol(rmatrixC,KcolC)
                call do_mat79matDmulDbleDble(DaA,KldA,KcolA&
                    &,rmatrixA%NEQ,DaB,DaC,KldC,KcolC)
              end if

            case (ST_SINGLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call lsyssc_getbase_Kld(rmatrixA,KldA)
              call lsyssc_getbase_Kcol(rmatrixA,KcolA)
              if (bfast) then
                call do_mat79matDmulDbleSngl(DaA,KldA,KcolA&
                    &,rmatrixA%NEQ,FaB,DaC)
              else
                call lsyssc_getbase_Kld(rmatrixC,KldC)
                call lsyssc_getbase_Kcol(rmatrixC,KcolC)
                call do_mat79matDmulDbleSngl(DaA,KldA,KcolA&
                    &,rmatrixA%NEQ,FaB,DaC,KldC,KcolC)
              end if

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case (ST_SINGLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call lsyssc_getbase_Kld(rmatrixA,KldA)
              call lsyssc_getbase_Kcol(rmatrixA,KcolA)
              if (bfast) then
                call do_mat79matDmulSnglDble(FaA,KldA,KcolA&
                    &,rmatrixA%NEQ,DaB,DaC)
              else
                call lsyssc_getbase_Kld(rmatrixC,KldC)
                call lsyssc_getbase_Kcol(rmatrixC,KcolC)
                call do_mat79matDmulSnglDble(FaA,KldA,KcolA&
                    &,rmatrixA%NEQ,DaB,DaC,KldC,KcolC)
              end if

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_single(rmatrixC,FaC)
              call lsyssc_getbase_Kld(rmatrixA,KldA)
              call lsyssc_getbase_Kcol(rmatrixA,KcolA)
              if (bfast) then
                call do_mat79matDmulSnglSngl(FaA,KldA,KcolA&
                    &,rmatrixA%NEQ,FaB,FaC)
              else
                call lsyssc_getbase_Kld(rmatrixC,KldC)
                call lsyssc_getbase_Kcol(rmatrixC,KcolC)
                call do_mat79matDmulSnglSngl(FaA,KldA,KcolA&
                    &,rmatrixA%NEQ,FaB,FaC,KldC,KcolC)
              end if

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case default
            call output_line("Unsupported data type!",&
                OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
            call sys_halt()
          end select
        end if

      case (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! B is CSR matrix - - - -

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixA,KldA)
        call lsyssc_getbase_Kcol(rmatrixA,KcolA)
        call lsyssc_getbase_Kld(rmatrixB,KldB)
        call lsyssc_getbase_Kcol(rmatrixB,KcolB)

        ! memory allocation?
        if (bmemory) then
          ! Duplicate structure of matrix A or B depending on which
          ! matrix is given in the encompassing matrix format
          if (rmatrixA%cmatrixFormat .ge. rmatrixB%cmatrixFormat) then
            call lsyssc_duplicateMatrix(rmatrixA,rmatrixC&
                &,LSYSSC_DUP_EMPTY,LSYSSC_DUP_REMOVE)
          else
            call lsyssc_duplicateMatrix(rmatrixB,rmatrixC&
                &,LSYSSC_DUP_EMPTY,LSYSSC_DUP_REMOVE)
          end if

          ! Set auxiliary pointers
          call storage_new("lsyssc_multMatMat","Kaux",rmatrixB%NCOLS,&
              &ST_INT,h_Kaux,ST_NEWBLOCK_NOINIT)
          call storage_getbase_int(h_Kaux,Kaux)

          ! Compute number of nonzero matrix entries: NA
          rmatrixC%NA = do_mat79mat79mul_computeNA(KldA,KcolA&
              &,rmatrixA%NEQ,KldB,KcolB,Kaux)
          call storage_new("lsyssc_multMatMat","h_Da",rmatrixC&
              &%NA,rmatrixC%cdataType,rmatrixC%h_Da,ST_NEWBLOCK_NOINIT)
          call storage_realloc("lsyssc_multMatMat",rmatrixC%NA&
              &,rmatrixC%h_Kcol,ST_NEWBLOCK_NOINIT,.false.)
        end if

        ! Check if matrix is given in the correct format
        if (rmatrixC%cmatrixFormat .ne. max(rmatrixA%cmatrixFormat,&
            &rmatrixB%cmatrixFormat)) then
          call output_line("Destination matrix has incompatible format",&
              OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
          call sys_halt()
        end if

        ! symbolical multiplication?
        if (bsymb) then

          ! Set pointers
          call lsyssc_getbase_Kld(rmatrixC,KldC)
          call lsyssc_getbase_Kcol(rmatrixC,KcolC)
          if (h_Kaux .eq. ST_NOHANDLE) then
            call storage_new("lsyssc_multMatMat","Kaux",max(rmatrixA&
                &%NEQ,max(rmatrixA%NCOLS,rmatrixB%NCOLS)),ST_INT&
                &,h_Kaux,ST_NEWBLOCK_NOINIT)
          else
            call storage_realloc("lsyssc_multMatMat",max(rmatrixA%NEQ&
                &,max(rmatrixA%NCOLS,rmatrixB%NCOLS)),h_Kaux&
                &,ST_NEWBLOCK_NOINIT,.false.)
          end if
          call storage_getbase_int(h_Kaux,Kaux)

          if (rmatrixC%cmatrixFormat .eq. LSYSSC_MATRIX9) then
            call lsyssc_getbase_Kdiagonal(rmatrixC,KdiagonalC)
            call do_mat79mat79mul_symb(rmatrixA%NEQ,rmatrixA%NCOLS&
                &,rmatrixB%NCOLS,KldA,KcolA,KldB,KcolB,KldC,KcolC&
                &,Kaux,KdiagonalC)
          else
            call do_mat79mat79mul_symb(rmatrixA%NEQ,rmatrixA%NCOLS&
                &,rmatrixB%NCOLS,KldA,KcolA,KldB,KcolB,KldC,KcolC&
                &,Kaux)
          end if
        end if

        ! numerical multiplication?
        if (bnumb) then

          ! Set pointers
          call lsyssc_getbase_Kld(rmatrixC,KldC)
          call lsyssc_getbase_Kcol(rmatrixC,KcolC)
          if (h_Daux .eq. ST_NOHANDLE) then
            call storage_new("lsyssc_multMatMat","Daux",max(rmatrixA&
                &%NEQ,max(rmatrixA%NCOLS,rmatrixB%NCOLS)),ST_DOUBLE&
                &,h_Daux,ST_NEWBLOCK_NOINIT)
          else
            call storage_realloc("lsyssc_multMatMat",max(rmatrixA%NEQ&
                &,max(rmatrixA%NCOLS,rmatrixB%NCOLS)),h_Daux&
                &,ST_NEWBLOCK_NOINIT,.false.)
          end if
          call storage_getbase_double(h_Daux,Daux)

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call do_mat79mat79mul_numb_dbledble(rmatrixA%NEQ&
                  &,rmatrixA%NCOLS,rmatrixB%NCOLS,KldA,KcolA,DaA,KldB&
                  &,KcolB,DaB,KldC,KcolC,DaC,Daux)

            case (ST_SINGLE)
              call lsyssc_getbase_double(rmatrixA,DaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call do_mat79mat79mul_numb_dblesngl(rmatrixA%NEQ&
                  &,rmatrixA%NCOLS,rmatrixB%NCOLS,KldA,KcolA,DaA,KldB&
                  &,KcolB,FaB,KldC,KcolC,DaC,Daux)

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case (ST_SINGLE)

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_double(rmatrixB,DaB)
              call lsyssc_getbase_double(rmatrixC,DaC)
              call do_mat79mat79mul_numb_sngldble(rmatrixA%NEQ&
                  &,rmatrixA%NCOLS,rmatrixB%NCOLS,KldA,KcolA,FaA,KldB&
                  &,KcolB,DaB,KldC,KcolC,DaC,Daux)

            case (ST_SINGLE)
              ! We need an additional vector FAUX
              call storage_new("lsyssc_multMatMat","Faux",max(rmatrixA&
                  &%NEQ,max(rmatrixA%NCOLS,rmatrixB%NCOLS)),ST_SINGLE&
                  &,h_Faux,ST_NEWBLOCK_NOINIT)
              call storage_getbase_single(h_Faux,Faux)

              call lsyssc_getbase_single(rmatrixA,FaA)
              call lsyssc_getbase_single(rmatrixB,FaB)
              call lsyssc_getbase_single(rmatrixC,FaC)
              call do_mat79mat79mul_numb_snglsngl(rmatrixA%NEQ&
                  &,rmatrixA%NCOLS,rmatrixB%NCOLS,KldA,KcolA,FaA,KldB&
                  &,KcolB,FaB,KldC,KcolC,FaC,Faux)

            case default
              call output_line("Unsupported data type!",&
                  OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
              call sys_halt()
            end select

          case default
            call output_line("Unsupported data type!",&
                OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
            call sys_halt()
          end select
        end if

      case default
        call output_line("Unsupported data type!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
        call sys_halt()
      end select

      !--------------------------------------------------------------

    case default
      call output_line("Unsupported data type!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_multMatMat")
      call sys_halt()
    end select

    ! Clear auxiliary vectors
    if (h_Kaux .ne. ST_NOHANDLE) call storage_free(h_Kaux)
    if (h_Daux .ne. ST_NOHANDLE) call storage_free(h_Daux)
    if (h_Faux .ne. ST_NOHANDLE) call storage_free(h_Faux)

  contains

    ! Here, the real MM multiplication routines follow.

    !**************************************************************
    ! Format 1 multiplication
    ! double precision matrix A (format 1)
    ! double precision matrix B (format 1)
    ! double precision matrix C (format 1)

    subroutine do_mat1mat1mulDbleDble(n,m,k,Da1,Da2,Da3)

      ! Remark: MATRIX1 is stored row-wise. The intrinsic function
      !         MATMUL requires the matrix to be stored columnwise
      !         Hence, compute C = A*B = (B"*A")"          (cpp fix: .")

      integer, intent(in) :: n,m,k
      real(DP), dimension(m,n), intent(in) :: Da1
      real(DP), dimension(k,m), intent(in) :: Da2
      real(DP), dimension(k,n), intent(inout) :: Da3

      ! Remark: The tuned BLAS3 routine performs much better than the
      ! intrinsic MATMUL Fortran90 routine. Hence, the BLAS routine
      ! is used whenever possible, that is, when all matrices have
      ! the same precision.
      call DGEMM("N","N",k,n,m,1.0_DP,Da2,k,Da1,m,0.0_DP,Da3,k)
    end subroutine

    !**************************************************************
    ! Format 1 multiplication
    ! double precision matrix A (format 1)
    ! single precision matrix B (format 1)
    ! double precision matrix C (format 1)

    subroutine do_mat1mat1mulDbleSngl(n,m,k,Da1,Fa2,Da3)

      ! Remark: MATRIX1 is stored row-wise. The intrinsic function
      !         MATMUL requires the matrix to be stored columnwise
      !         Hence, compute C = A*B = (B"*A")"        (cpp fix: .")

      integer, intent(in) :: n,m,k
      real(DP), dimension(m,n), intent(in) :: Da1
      real(SP), dimension(k,m), intent(in) :: Fa2
      real(DP), dimension(k,n), intent(inout) :: Da3

      call DGEMM("N","N",k,n,m,1.0_DP,real(Fa2,DP),k,Da1,m,0.0_DP,Da3,k)
!!$      Da3=MATMUL(Fa2,Da1)
    end subroutine

    !**************************************************************
    ! Format 1 multiplication
    ! single precision matrix A (format 1)
    ! double precision matrix B (format 1)
    ! double precision matrix C (format 1)

    subroutine do_mat1mat1mulSnglDble(n,m,k,Fa1,Da2,Da3)

      ! Remark: MATRIX1 is stored row-wise. The intrinsic function
      !         MATMUL requires the matrix to be stored columnwise
      !         Hence, compute C = A*B = (B"*A")"        (cpp fix: .")

      integer, intent(in) :: n,m,k
      real(SP), dimension(m,n), intent(in) :: Fa1
      real(DP), dimension(k,m), intent(in) :: Da2
      real(DP), dimension(k,n), intent(inout) :: Da3

      call DGEMM("N","N",k,n,m,1.0_DP,Da2,k,real(Fa1,DP),m,0.0_DP,Da3,k)
!!$      Da3=MATMUL(Da2,Fa1)
    end subroutine

    !**************************************************************
    ! Format 1 multiplication
    ! single precision matrix A (format 1)
    ! single precision matrix B (format 1)
    ! single precision matrix C (format 1)

    subroutine do_mat1mat1mulSnglSngl(n,m,k,Fa1,Fa2,Fa3)

      ! Remark: MATRIX1 is stored row-wise. The intrinsic function
      !         MATMUL requires the matrix to be stored columnwise
      !         Hence, compute C = A*B = (B"*A")"        (cpp fix: .")

      integer, intent(in) :: n,m,k
      real(SP), dimension(m,n), intent(in) :: Fa1
      real(SP), dimension(k,m), intent(in) :: Fa2
      real(SP), dimension(k,n), intent(inout) :: Fa3

      ! Remark: The tuned BLAS3 routine performs much better than the
      ! intrinsic MATMUL Fortran90 routine. Hence, the BLAS routine
      ! is used whenever possible, that is, when all matrices have
      ! the same precision.
      call SGEMM("N","N",k,n,m,1E0,Fa2,k,Fa1,m,0E0,Fa3,k)
    end subroutine

    !**************************************************************
    ! Format D-1 multiplication
    ! double precision matrix A (format D)
    ! double precision matrix B (format 1)
    ! double precision matrix C (format 1)

    subroutine do_matDmat1mulDbleDble(n,m,Da1,Da2,Da3)

      integer, intent(in) :: n,m
      real(DP), dimension(n), intent(in) :: Da1
      real(DP), dimension(m,n), intent(in) :: Da2
      real(DP), dimension(m,n), intent(inout) :: Da3

      integer :: i

      do i=1,n
        Da3(:,i)=Da1(i)*Da2(:,i)
      end do
    end subroutine

    !**************************************************************
    ! Format D-1 multiplication
    ! single precision matrix A (format D)
    ! double precision matrix B (format 1)
    ! double precision matrix C (format 1)

    subroutine do_matDmat1mulSnglDble(n,m,Fa1,Da2,Da3)

      integer, intent(in) :: n,m
      real(SP), dimension(n), intent(in) :: Fa1
      real(DP), dimension(m,n), intent(in) :: Da2
      real(DP), dimension(m,n), intent(inout) :: Da3

      integer :: i

      do i=1,n
        Da3(:,i)=Fa1(i)*Da2(:,i)
      end do
    end subroutine

    !**************************************************************
    ! Format D-1 multiplication
    ! double precision matrix A (format D)
    ! single precision matrix B (format 1)
    ! double precision matrix C (format 1)

    subroutine do_matDmat1mulDbleSngl(n,m,Da1,Fa2,Da3)

      integer, intent(in) :: n,m
      real(DP), dimension(n), intent(in) :: Da1
      real(SP), dimension(m,n), intent(in) :: Fa2
      real(DP), dimension(m,n), intent(inout) :: Da3

      integer :: i

      do i=1,n
        Da3(:,i)=Da1(i)*Fa2(:,i)
      end do
    end subroutine

    !**************************************************************
    ! Format D-1 multiplication
    ! single precision matrix A (format D)
    ! single precision matrix B (format 1)
    ! single precision matrix C (format 1)

    subroutine do_matDmat1mulSnglSngl(n,m,Fa1,Fa2,Fa3)

      integer, intent(in) :: n,m
      real(SP), dimension(n), intent(in) :: Fa1
      real(SP), dimension(m,n), intent(in) :: Fa2
      real(SP), dimension(m,n), intent(inout) :: Fa3

      integer :: i

      do i=1,n
        Fa3(:,i)=Fa1(i)*Fa2(:,i)
      end do
    end subroutine

    !**************************************************************
    ! Format 1-D multiplication
    ! double precision matrix A (format 1)
    ! double precision matrix B (format D)
    ! double precision matrix C (format 1)

    subroutine do_mat1matDmulDbleDble(n,m,Da1,Da2,Da3)

      integer, intent(in) :: n,m
      real(DP), dimension(m,n), intent(in) :: Da1
      real(DP), dimension(m),   intent(in) :: Da2
      real(DP), dimension(m,n), intent(inout) :: Da3

      integer :: i

      do i=1,m
        Da3(i,:)=Da1(i,:)*Da2(i)
      end do
    end subroutine

    !**************************************************************
    ! Format 1-D multiplication
    ! single precision matrix A (format 1)
    ! double precision matrix B (format D)
    ! double precision matrix C (format 1)

    subroutine do_mat1matDmulSnglDble(n,m,Fa1,Da2,Da3)

      integer, intent(in) :: n,m
      real(SP), dimension(m,n), intent(in) :: Fa1
      real(DP), dimension(m),   intent(in) :: Da2
      real(DP), dimension(m,n), intent(inout) :: Da3

      integer :: i

      do i=1,m
        Da3(i,:)=Fa1(i,:)*Da2(i)
      end do
    end subroutine

    !**************************************************************
    ! Format 1-D multiplication
    ! double precision matrix A (format 1)
    ! single precision matrix B (format D)
    ! double precision matrix C (format 1)

    subroutine do_mat1matDmulDbleSngl(n,m,Da1,Fa2,Da3)

      integer, intent(in) :: n,m
      real(DP), dimension(m,n), intent(in) :: Da1
      real(SP), dimension(m),   intent(in) :: Fa2
      real(DP), dimension(m,n), intent(inout) :: Da3

      integer :: i

      do i=1,m
        Da3(i,:)=Da1(i,:)*Fa2(i)
      end do
    end subroutine

    !**************************************************************
    ! Format 1-D multiplication
    ! single precision matrix A (format 1)
    ! single precision matrix B (format D)
    ! single precision matrix C (format 1)

    subroutine do_mat1matDmulSnglSngl(n,m,Fa1,Fa2,Fa3)

      integer, intent(in) :: n,m
      real(SP), dimension(m,n), intent(in) :: Fa1
      real(SP), dimension(m),   intent(in) :: Fa2
      real(SP), dimension(m,n), intent(inout) :: Fa3

      integer :: i

      do i=1,m
        Fa3(i,:)=Fa1(i,:)*Fa2(i)
      end do
    end subroutine

    !**************************************************************
    ! Format D-7/9 multiplication
    ! double precision matrix A (format D)
    ! double precision matrix B (format 7 or format 9)
    ! double precision matrix C (format 7 or format 9)

    subroutine do_matDmat79mulDbleDble(Da1,Da2,Kld2,Kcol2&
        &,neq,Da3,Kld3,Kcol3)

      real(DP), dimension(:), intent(in) :: Da1
      real(DP), dimension(:), intent(in) :: Da2
      real(DP), dimension(:), intent(inout) :: Da3
      integer, dimension(:), intent(in) :: Kld2
      integer, dimension(:), intent(in) :: Kcol2
      integer, dimension(:), intent(in), optional :: Kld3
      integer, dimension(:), intent(in), optional :: Kcol3
      integer, intent(in) :: neq

      integer :: ieq,ild2,ild3,ildend3

      if (present(Kld3) .and. present(Kcol3)) then

        do ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          do ild2=Kld2(ieq),Kld2(ieq+1)-1
            do ild3=ild3,ildend3
              if (Kcol3(ild3) .eq. Kcol2(ild2)) exit
              Da3(ild3)=0
            end do
            Da3(ild3)=Da1(ieq)*Da2(ild2); ild3=ild3+1
          end do
          Da3(ild3:ildend3)=0
        end do

      else

        do ieq=1,neq
          do ild2=Kld2(ieq),Kld2(ieq+1)-1
            Da3(ild2)=Da1(ieq)*Da2(ild2)
          end do
        end do

      end if
    end subroutine

    !**************************************************************
    ! Format D-7/9 multiplication
    ! single precision matrix A (format D)
    ! double precision matrix B (format 7 or format 9)
    ! double precision matrix C (format 7 or format 9)

    subroutine do_matDmat79mulSnglDble(Fa1,Da2,Kld2,Kcol2&
        &,neq,Da3,Kld3,Kcol3)

      real(SP), dimension(:), intent(in) :: Fa1
      real(DP), dimension(:), intent(in) :: Da2
      real(DP), dimension(:), intent(inout) :: Da3
      integer, dimension(:), intent(in) :: Kld2
      integer, dimension(:), intent(in) :: Kcol2
      integer, dimension(:), intent(in), optional :: Kld3
      integer, dimension(:), intent(in), optional :: Kcol3
      integer, intent(in) :: neq

      integer :: ieq,ild2,ild3,ildend3

      if (present(Kld3) .and. present(Kcol3)) then

        do ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          do ild2=Kld2(ieq),Kld2(ieq+1)-1
            do ild3=ild3,ildend3
              if (Kcol3(ild3) .eq. Kcol2(ild2)) exit
              Da3(ild3)=0
            end do
            Da3(ild3)=Fa1(ieq)*Da2(ild2); ild3=ild3+1
          end do
          Da3(ild3:ildend3)=0
        end do

      else

        do ieq=1,neq
          do ild2=Kld2(ieq),Kld2(ieq+1)-1
            Da3(ild2)=Fa1(ieq)*Da2(ild2)
          end do
        end do

      end if
    end subroutine

    !**************************************************************
    ! Format D-7/9 multiplication
    ! double precision matrix A (format D)
    ! single precision matrix B (format 7 or format 9)
    ! double precision matrix C (format 7 or format 9)

    subroutine do_matDmat79mulDbleSngl(Da1,Fa2,Kld2,Kcol2&
        &,neq,Da3,Kld3,Kcol3)

      real(DP), dimension(:), intent(in) :: Da1
      real(SP), dimension(:), intent(in) :: Fa2
      real(DP), dimension(:), intent(inout) :: Da3
      integer, dimension(:), intent(in) :: Kld2
      integer, dimension(:), intent(in) :: Kcol2
      integer, dimension(:), intent(in), optional :: Kld3
      integer, dimension(:), intent(in), optional :: Kcol3
      integer, intent(in) :: neq

      integer :: ieq,ild2,ild3,ildend3

      if (present(Kld3) .and. present(Kcol3)) then

        do ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          do ild2=Kld2(ieq),Kld2(ieq+1)-1
            do ild3=ild3,ildend3
              if (Kcol3(ild3) .eq. Kcol2(ild2)) exit
              Da3(ild3)=0
            end do
            Da3(ild3)=Da1(ieq)*Fa2(ild2); ild3=ild3+1
          end do
          Da3(ild3:ildend3)=0
        end do

      else

        do ieq=1,neq
          do ild2=Kld2(ieq),Kld2(ieq+1)-1
            Da3(ild2)=Da1(ieq)*Fa2(ild2)
          end do
        end do

      end if
    end subroutine

    !**************************************************************
    ! Format D-7/9 multiplication
    ! single precision matrix A (format D)
    ! single precision matrix B (format 7 or format 9)
    ! single precision matrix C (format 7 or format 9)

    subroutine do_matDmat79mulSnglSngl(Fa1,Fa2,Kld2,Kcol2&
        &,neq,Fa3,Kld3,Kcol3)

      real(SP), dimension(:), intent(in) :: Fa1
      real(SP), dimension(:), intent(in) :: Fa2
      real(SP), dimension(:), intent(inout) :: Fa3
      integer, dimension(:), intent(in) :: Kld2
      integer, dimension(:), intent(in) :: Kcol2
      integer, dimension(:), intent(in), optional :: Kld3
      integer, dimension(:), intent(in), optional :: Kcol3
      integer, intent(in) :: neq

      integer :: ieq,ild2,ild3,ildend3

      if (present(Kld3) .and. present(Kcol3)) then

        do ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          do ild2=Kld2(ieq),Kld2(ieq+1)-1
            do ild3=ild3,ildend3
              if (Kcol3(ild3) .eq. Kcol2(ild2)) exit
              Fa3(ild3)=0
            end do
            Fa3(ild3)=Fa1(ieq)*Fa2(ild2); ild3=ild3+1
          end do
          Fa3(ild3:ildend3)=0
        end do

      else

        do ieq=1,neq
          do ild2=Kld2(ieq),Kld2(ieq+1)-1
            Fa3(ild2)=Fa1(ieq)*Fa2(ild2)
          end do
        end do

      end if
    end subroutine

    !**************************************************************
    ! Format 7/9-D multiplication
    ! double precision matrix A (format 7 or format 9)
    ! double precision matrix B (format D)
    ! double precision matrix C (format 7 or format 9)

    subroutine do_mat79matDmulDbleDble(Da1,Kld1,Kcol1,neq,Da2&
        &,Da3,Kld3,Kcol3)

      real(DP), dimension(:), intent(in) :: Da1
      real(DP), dimension(:), intent(in) :: Da2
      real(DP), dimension(:), intent(inout) :: Da3
      integer, dimension(:), intent(in) :: Kld1
      integer, dimension(:), intent(in) :: Kcol1
      integer, dimension(:), intent(in), optional :: Kld3
      integer, dimension(:), intent(in), optional :: Kcol3
      integer, intent(in) :: neq

      integer :: ieq,ild1,ild3,ildend3

      if (present(Kld3) .and. present(Kcol3)) then

        do ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          do ild1=Kld1(ieq),Kld1(ieq+1)-1
            do ild3=ild3,ildend3
              if (Kcol3(ild3) .eq. Kcol1(ild1)) exit
              Da3(ild3)=0
            end do
            Da3(ild3)=Da1(ild1)*Da2(Kcol1(ild1)); ild3=ild3+1
          end do
          Da3(ild3:ildend3)=0
        end do

      else

        do ieq=1,neq
          do ild1=Kld1(ieq),Kld1(ieq+1)-1
            Da3(ild1)=Da1(ild1)*Da2(Kcol1(ild1))
          end do
        end do

      end if
    end subroutine

    !**************************************************************
    ! Format 7/9-D multiplication
    ! single precision matrix A (format 7 or format 9)
    ! double precision matrix B (format D)
    ! double precision matrix C (format 7 or format 9)

    subroutine do_mat79matDmulSnglDble(Fa1,Kld1,Kcol1,neq,Da2&
        &,Da3,Kld3,Kcol3)

      real(SP), dimension(:), intent(in) :: Fa1
      real(DP), dimension(:), intent(in) :: Da2
      real(DP), dimension(:), intent(inout) :: Da3
      integer, dimension(:), intent(in) :: Kld1
      integer, dimension(:), intent(in) :: Kcol1
      integer, dimension(:), intent(in), optional :: Kld3
      integer, dimension(:), intent(in), optional :: Kcol3
      integer, intent(in) :: neq

      integer :: ieq,ild1,ild3,ildend3

      if (present(Kld3) .and. present(Kcol3)) then

        do ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          do ild1=Kld1(ieq),Kld1(ieq+1)-1
            do ild3=ild3,ildend3
              if (Kcol3(ild3) .eq. Kcol1(ild1)) exit
              Da3(ild3)=0
            end do
            Da3(ild3)=Fa1(ild1)*Da2(Kcol1(ild1)); ild3=ild3+1
          end do
          Da3(ild3:ildend3)=0
        end do

      else

        do ieq=1,neq
          do ild1=Kld1(ieq),Kld1(ieq+1)-1
            Da3(ild1)=Fa1(ild1)*Da2(Kcol1(ild1))
          end do
        end do

      end if
    end subroutine

    !**************************************************************
    ! Format 7/9-D multiplication
    ! double precision matrix A (format 7 or format 9)
    ! single precision matrix B (format D)
    ! double precision matrix C (format 7 or format 9)

    subroutine do_mat79matDmulDbleSngl(Da1,Kld1,Kcol1,neq,Fa2&
        &,Da3,Kld3,Kcol3)

      real(DP), dimension(:), intent(in) :: Da1
      real(SP), dimension(:), intent(in) :: Fa2
      real(DP), dimension(:), intent(inout) :: Da3
      integer, dimension(:), intent(in) :: Kld1
      integer, dimension(:), intent(in) :: Kcol1
      integer, dimension(:), intent(in), optional :: Kld3
      integer, dimension(:), intent(in), optional :: Kcol3
      integer, intent(in) :: neq

      integer :: ieq,ild1,ild3,ildend3

      if (present(Kld3) .and. present(Kcol3)) then

        do ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          do ild1=Kld1(ieq),Kld1(ieq+1)-1
            do ild3=ild3,ildend3
              if (Kcol3(ild3) .eq. Kcol1(ild1)) exit
              Da3(ild3)=0
            end do
            Da3(ild3)=Da1(ild1)*Fa2(Kcol1(ild1)); ild3=ild3+1
          end do
          Da3(ild3:ildend3)=0
        end do

      else

        do ieq=1,neq
          do ild1=Kld1(ieq),Kld1(ieq+1)-1
            Da3(ild1)=Da1(ild1)*Fa2(Kcol1(ild1))
          end do
        end do

      end if
    end subroutine

    !**************************************************************
    ! Format 7/9-D multiplication
    ! single precision matrix A (format 7 or format 9)
    ! single precision matrix B (format D)
    ! single precision matrix C (format 7 or format 9)

    subroutine do_mat79matDmulSnglSngl(Fa1,Kld1,Kcol1,neq,Fa2&
        &,Fa3,Kld3,Kcol3)

      real(SP), dimension(:), intent(in) :: Fa1
      real(SP), dimension(:), intent(in) :: Fa2
      real(SP), dimension(:), intent(inout) :: Fa3
      integer, dimension(:), intent(in) :: Kld1
      integer, dimension(:), intent(in) :: Kcol1
      integer, dimension(:), intent(in), optional :: Kld3
      integer, dimension(:), intent(in), optional :: Kcol3
      integer, intent(in) :: neq

      integer :: ieq,ild1,ild3,ildend3

      if (present(Kld3) .and. present(Kcol3)) then

        do ieq=1,neq
          ild3=Kld3(ieq); ildend3=Kld3(ieq+1)-1
          do ild1=Kld1(ieq),Kld1(ieq+1)-1
            do ild3=ild3,ildend3
              if (Kcol3(ild3) .eq. Kcol1(ild1)) exit
              Fa3(ild3)=0
            end do
            Fa3(ild3)=Fa1(ild1)*Fa2(Kcol1(ild1)); ild3=ild3+1
          end do
          Fa3(ild3:ildend3)=0
        end do

      else

        do ieq=1,neq
          do ild1=Kld1(ieq),Kld1(ieq+1)-1
            Fa3(ild1)=Fa1(ild1)*Fa2(Kcol1(ild1))
          end do
        end do

      end if
    end subroutine

    !**************************************************************
    ! Format 7/9-7/9 multiplication
    ! Compute the number of nonzero matrix entries of C := A * B
    !
    ! Method: A' * A = sum [over i=1, nrow] a(i)^T a(i)         (cpp fix: .')
    !         where a(i) = i-th row of A. We must be careful not
    !         to add the elements already accounted for.
    !
    ! Remark: This subroutine is a modified version of the AMUBDG
    !         subroutine taken from the SPARSEKIT library written
    !         by Youcef Saad.

    function do_mat79mat79mul_computeNA(KldA,KcolA,neq,KldB,KcolB&
        &,Kaux) result(NA)

      integer, dimension(:), intent(in) :: KldA,KcolA,KldB,KcolB
      integer, dimension(:), intent(inout) :: Kaux
      integer, intent(in) :: neq
      integer :: NA

      integer :: ieq,jeq,ild,irow,icol,idg,ndg,last

      ! Initialization
      Kaux=0; NA=0

      do ieq=1,neq

        ! For each row of matrix A
        ndg=0

        ! End-of-linked list
        last=-1

        do ild=KldA(ieq),KldA(ieq+1)-1

          ! Row number to be added
          irow=KcolA(ild)

          do jeq=KldB(irow),KldB(irow+1)-1
            icol=KcolB(jeq)

            ! Add one element to the linked list
            if (Kaux(icol) .eq. 0) then
              ndg = ndg+1
              Kaux(icol) = last
              last = icol
            end if

          end do
        end do

        ! We are done with row IEQ
        NA = NA+ndg

        ! Reset Kaux to zero
        do idg=1,ndg
          jeq = Kaux(last)
          Kaux(last) = 0
          last = jeq
        end do

      end do
    end function do_mat79mat79mul_computeNA

    !**************************************************************
    ! Format 7/9-7/9 multiplication
    ! Perform symbolical matrix-matrix-multiplication C := A * B
    !
    ! This subroutine is based on the SYMBMM routine from the
    ! Sparse Matrix Multiplication Package (SMMP) written by
    ! R.E. Bank and C.C. Douglas which is freely available at:
    ! http://cs-www.cs.yale.edu/homes/douglas-craig/Codes/smmp.tgz

    subroutine do_mat79mat79mul_symb(n,m,l,KldA,KcolA,KldB,KcolB,KldC&
        &,KcolC,Kindex,KdiagonalC)

      integer, intent(in) :: n,m,l
      integer, dimension(:), intent(in) :: KldA,KldB
      integer, dimension(:), intent(in) :: KcolA,KcolB
      integer, dimension(:), intent(inout) :: KldC,Kindex
      integer, dimension(:), intent(inout) :: KcolC
      integer, dimension(:), intent(inout), optional :: KdiagonalC

      integer :: i,j,k,istart,ilength,jj,jk,kaux,ioff

      ! KINDEX is used to store links. If an entry in KINDEX os
      ! nonzero, it is a pointer to the next column with a nonzero.
      ! The links are determined as they are found, and are unordered.
      Kindex = 0
      KldC(1)= 1

      ! The main loop consists of three components: initialization, a
      ! long loop that merges row lists, and code to copy the links
      ! into the KCOLC vector.
      do i=1,n

        ! Initialization
        istart = -1
        ilength = 0

        ! The start column ISTART is reset and the number of column
        ! entries for the I-th row is assumed empty. Now, merge the
        ! row lists as follows
        do jj=KldA(i),KldA(i+1)-1
          j=KcolA(jj)

          ! Determine the intersection of the current row I in matrix
          ! A with the nonzeros in column J of matrix B
          do k=KldB(j),KldB(j+1)-1
            jk=KcolB(k)
            if (Kindex(jk) .eq. 0) then
              Kindex(jk) = istart
              istart = jk
              ilength = ilength+1
            end if
          end do
        end do

        ! Copy the links into the KCOLC vector as column indices
        ! Since the column vectors are given in unordered way, they
        ! must be ordered before inserted inco KCOLC
        KldC(i+1)=KldC(i)+ilength

        if (present(KdiagonalC)) then
          ! If KDIAGONALC is present, then the matrix C is stored in
          ! CSR9 format, that is, all entries in KCOLC are numbered
          ! continuously but the position of the diagonal entry is
          ! kept in KDIAGONALC
          do j=KldC(i),KldC(i+1)-1
            KcolC(j) = istart
            istart = Kindex(istart)
            Kindex(KcolC(j)) = 0

            ! Perform local insertionsort to find the correct
            ! position of the provisional entry KCOLC(J)
            kaux = KcolC(j)
            do jj=j-1,KldC(i),-1
              if (KcolC(jj) .le. kaux) goto 10
              KcolC(jj+1) = KcolC(jj)
              if (KcolC(jj+1) .eq. i) KdiagonalC(i) = jj+1
            end do
            jj=KldC(i)-1
10          KcolC(jj+1) = kaux
            if (KcolC(jj+1) .eq. i) KdiagonalC(i) = jj+1
          end do

        else
          ! If KDIAGONALC is not present, then the matrix C is stored
          ! in CSR7 format, that is, all entries in KCOLC are
          ! numbered continuously except for the diagonal entry which
          ! is stored in the first position KLDC(I) of row I.

          ! The variable IOFF is used to indicate, whether the
          ! diagonal entry has already been inserted or not.
          ioff=0

          do j=KldC(i),KldC(i+1)-1
            KcolC(j) = istart
            istart = Kindex(istart)
            Kindex(KcolC(j)) = 0

            if (KcolC(j) .eq. i) then
              ! The current entry is the diagonal entry, hence, shift
              ! all entries to the right by one and move the current
              ! entry to the first position
              KcolC(KldC(i):j)=cshift(KcolC(KldC(i):j),-1)

              ! Set IOFF to one so that the first entry of the I-th
              ! row is neglected in subsequent iterations
              ioff=1
            else
              ! Perform local insertsiosort to find the correct
              ! position of the provisional entry KCOLC(J)
              kaux = KcolC(j); jj=j
              do jj=j-1,KldC(i),-1
                if (KcolC(jj) .le. kaux) goto 20
                KcolC(jj+1) = KcolC(jj)
              end do
              jj=KldC(i)-1
20            KcolC(jj+1) = kaux
            end if
          end do
        end if

        Kindex(i) = 0
      end do
    end subroutine

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

    subroutine do_mat79mat79mul_numb_dbledble(n,m,l,KldA,KcolA&
        &,DaA,KldB,KcolB,DaB,KldC,KcolC,DaC,Dtemp)

      integer, intent(in) :: n,m,l
      integer, dimension(:), intent(in) :: KldA,KldB
      integer, dimension(:), intent(in) :: KcolA,KcolB
      integer, dimension(:), intent(inout) :: KldC
      integer, dimension(:), intent(inout) :: KcolC
      real(DP), dimension(:), intent(in) :: DaA
      real(DP), dimension(:), intent(in) :: DaB
      real(DP), dimension(:), intent(inout) :: DaC,Dtemp

      integer :: i,j,k,jj,jk
      real(DP) :: daux

      ! DTEMP is used to store partial sums.
      Dtemp = 0

      ! The main loop forms the partial sums and then copies the
      ! completed sums into the correct locations in the sparse matrix
      do i=1,n

        do jj=KldA(i),KldA(i+1)-1
          j = KcolA(jj)
          daux = DaA(jj)

          ! Accumulate the product for row J
          do k=KldB(j),KldB(j+1)-1
            jk=KcolB(k)
            Dtemp(jk)=Dtemp(jk)+daux*DaB(k)
          end do
        end do

        ! Store product for row J
        do j=KldC(i),KldC(i+1)-1
          jj=KcolC(j)
          DaC(j) = Dtemp(jj)
          Dtemp(jj) = 0
        end do
      end do
    end subroutine

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

    subroutine do_mat79mat79mul_numb_sngldble(n,m,l,KldA,KcolA&
        &,FaA,KldB,KcolB,DaB,KldC,KcolC,DaC,Dtemp)

      integer, intent(in) :: n,m,l
      integer, dimension(:), intent(in) :: KldA,KldB
      integer, dimension(:), intent(in) :: KcolA,KcolB
      integer, dimension(:), intent(inout) :: KldC
      integer, dimension(:), intent(inout) :: KcolC
      real(SP), dimension(:), intent(in) :: FaA
      real(DP), dimension(:), intent(in) :: DaB
      real(DP), dimension(:), intent(inout) :: DaC,Dtemp

      integer :: i,j,k,jj,jk
      real(DP) :: daux

      ! DTEMP is used to store partial sums.
      Dtemp = 0

      ! The main loop forms the partial sums and then copies the
      ! completed sums into the correct locations in the sparse matrix
      do i=1,n

        do jj=KldA(i),KldA(i+1)-1
          j = KcolA(jj)
          daux = FaA(jj)

          ! Accumulate the product for row J
          do k=KldB(j),KldB(j+1)-1
            jk=KcolB(k)
            Dtemp(jk)=Dtemp(jk)+daux*DaB(k)
          end do
        end do

        ! Store product for row J
        do j=KldC(i),KldC(i+1)-1
          jj=KcolC(j)
          DaC(j) = Dtemp(jj)
          Dtemp(jj) = 0
        end do
      end do
    end subroutine

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

    subroutine do_mat79mat79mul_numb_dblesngl(n,m,l,KldA,KcolA&
        &,DaA,KldB,KcolB,FaB,KldC,KcolC,DaC,Dtemp)

      integer, intent(in) :: n,m,l
      integer, dimension(:), intent(in) :: KldA,KldB
      integer, dimension(:), intent(in) :: KcolA,KcolB
      integer, dimension(:), intent(inout) :: KldC
      integer, dimension(:), intent(inout) :: KcolC
      real(DP), dimension(:), intent(in) :: DaA
      real(SP), dimension(:), intent(in) :: FaB
      real(DP), dimension(:), intent(inout) :: DaC,Dtemp

      integer :: i,j,k,jj,jk
      real(DP) :: daux

      ! DTEMP is used to store partial sums.
      Dtemp = 0

      ! The main loop forms the partial sums and then copies the
      ! completed sums into the correct locations in the sparse matrix
      do i=1,n

        do jj=KldA(i),KldA(i+1)-1
          j = KcolA(jj)
          daux = DaA(jj)

          ! Accumulate the product for row J
          do k=KldB(j),KldB(j+1)-1
            jk=KcolB(k)
            Dtemp(jk)=Dtemp(jk)+daux*FaB(k)
          end do
        end do

        ! Store product for row J
        do j=KldC(i),KldC(i+1)-1
          jj=KcolC(j)
          DaC(j) = Dtemp(jj)
          Dtemp(jj) = 0
        end do
      end do
    end subroutine

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

    subroutine do_mat79mat79mul_numb_snglsngl(n,m,l,KldA,KcolA&
        &,FaA,KldB,KcolB,FaB,KldC,KcolC,FaC,Ftemp)

      integer, intent(in) :: n,m,l
      integer, dimension(:), intent(in) :: KldA,KldB
      integer, dimension(:), intent(in) :: KcolA,KcolB
      integer, dimension(:), intent(inout) :: KldC
      integer, dimension(:), intent(inout) :: KcolC
      real(SP), dimension(:), intent(in) :: FaA
      real(SP), dimension(:), intent(in) :: FaB
      real(SP), dimension(:), intent(inout) :: FaC,Ftemp

      integer :: i,j,k,jj,jk
      real(SP) :: faux

      ! FTEMP is used to store partial sums.
      Ftemp = 0

      ! The main loop forms the partial sums and then copies the
      ! completed sums into the correct locations in the sparse matrix
      do i=1,n

        do jj=KldA(i),KldA(i+1)-1
          j = KcolA(jj)
          faux = FaA(jj)

          ! Accumulate the product for row J
          do k=KldB(j),KldB(j+1)-1
            jk=KcolB(k)
            Ftemp(jk)=Ftemp(jk)+faux*FaB(k)
          end do
        end do

        ! Store product for row J
        do j=KldC(i),KldC(i+1)-1
          jj=KcolC(j)
          FaC(j) = Ftemp(jj)
          Ftemp(jj) = 0
        end do
      end do
    end subroutine
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_matrixLinearComb2 (rmatrixA,cA,rmatrixB,cB,rmatrixC,&
      bmemory,bsymb,bnumb,bisExactStructure,rperfconfig)

    !<description>
    ! Performs a linear combination:
    !   rmatrixC = ca*rmatrixA + cB*rmatrixB
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
    !
    ! WARNING: This routine is deprecated and will be removed in future.
    !          Use the (in some cases) much more efficient subroutine
    !          lsyssc_matrixLinearComb2 as replacement.
    !</description>

!<input>
    ! First source matrix
    type(t_matrixScalar), intent(in) :: rmatrixA

    ! scaling factors
    real(DP), intent(in) :: ca,cb

    ! BMEMORY = FALSE: Do not allocate required memory for C=A+B.
    ! BMEMORY = TRUE:  Generate all required structures for C=A+B
    logical, intent(in) :: bmemory

    ! Compute symbolic matrix-matrix-addition
    ! BSYMB = FALSE: Do not generate the required matrix structures.
    !                This may be useful, if the same matrices A and B
    !                need to be added several times, but the
    !                sparsity patterns do not change
    ! BSYMN = TRUE:  Generate all required matrix structures for C=A+B
    logical, intent(in) :: bsymb

    ! Compute numerical matrix-matrix-addition
    ! BNUMB = FALSE: Do not perform the numerical addition
    ! BNUMB = TRUE:  Perform numerical addition
    logical, intent(in) :: bnumb

    ! OPTIONAL: Indicates whether the resulting matrix C has the
    ! required symbolic structure or is a superset of the symbolic
    ! matrix-matrix product. In some cases, this may allow for much
    ! more efficient implementation.
    ! Standard: FALSE
    logical, intent(in), optional :: bisExactStructure

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Second source matrix
    type(t_matrixScalar), intent(inout) :: rmatrixB

    ! Output matrix. May coincide with rmatrixB.
    type(t_matrixScalar), intent(inout) :: rmatrixC
!</inputoutput>

!</subroutine>

#if WARN_DEPREC
    call output_line ("Using deprecated feature. Please update your code.", &
        OU_CLASS_WARNING,OU_MODE_STD,"lsyssc_matrixLinearComb2")
#endif

    call lsyssc_matrixLinearComb(rmatrixA,rmatrixB,ca,cb,&
        bmemory,bsymb,bnumb,bisExactStructure,rmatrixC,rperfconfig)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_matrixLinearComb (rmatrixA,rmatrixB,ca,cb,&
      bmemory,bsymb,bnumb,bisExactStructure,rdest,rperfconfig)

    !<description>
    ! Performs a linear combination:
    !   rdest or rmatrixB = ca*rmatrixA + cB*rmatrixB
    !
    ! Most combinations of matrix formats are supported for matrices A
    ! and B. If the destination matrix rmatrixDest is not given, then
    ! matrix rmatrrixB must be combatible with matrix rmatrixA, that
    ! is, all non-zero entries of matrix rmatrixA must be present in
    ! rmatrixB. If the destination matrix rmatrixC is given and
    ! BMEMORY=FALSE, then matrix rmatrixC must be compatible to both
    ! matrices rmatrixA and rmatrixB. If BMEMORY=TRUE, then the
    ! destination matrix rmatrixC is generated on-the-fly with the
    ! correct sparsity pattern to comprise the non-zero matrix entries
    ! of rmatrixA and rmatrixB.
    ! </description>

!<input>
    ! First source matrix
    type(t_matrixScalar), intent(in) :: rmatrixA

    ! scaling factors
    real(DP), intent(in) :: ca,cb

    ! OPTIONAL Allocate memory for matrix structure and/or data
    ! BMEMORY = FALSE: Do not allocate required memory for C=A+B.
    ! BMEMORY = TRUE:  Generate all required structures for C=A+B
    ! Default: FALSE
    logical, intent(in), optional :: bmemory

    ! OPTIONAL: Compute symbolic matrix-matrix-addition
    ! BSYMB = FALSE: Do not generate the required matrix structures.
    !                This may be useful, if the same matrices A and B
    !                need to be added several times, but the
    !                sparsity patterns do not change
    ! BSYMN = TRUE:  Generate all required matrix structures for C=A+B
    ! Default: FALSE
    logical, intent(in), optional :: bsymb

    ! OPTIONAL: Compute numerical matrix-matrix-addition
    ! BNUMB = FALSE: Do not perform the numerical addition
    ! BNUMB = TRUE:  Perform numerical addition
    ! Default: TRUE
    logical, intent(in), optional :: bnumb

    ! OPTIONAL: Indicates whether the resulting matrix has the
    ! required symbolic structure or is a superset of the symbolic
    ! matrix-matrix product. In some cases, this may allow for much
    ! more efficient implementation.
    ! Standard: FALSE
    logical, intent(in), optional :: bisExactStructure

    ! OPTIONAL: local performance configuration. If not given, the
    ! global performance configuration is used.
    type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
    ! Second source vector; also receives the result if rdest is not specified.
    type(t_matrixScalar), intent(inout), target :: rmatrixB

    ! OPTIONAL: Destination matrix. If not specified, rmatrixB will be overwritten.
    type(t_matrixScalar), intent(inout), target, optional :: rdest
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_matrixScalar), pointer :: p_rdest
    real(DP), dimension(:), pointer :: p_DaA,p_DaB,p_DaC
    real(SP), dimension(:), pointer :: p_FaA,p_FaB,p_FaC
    integer, dimension(:), pointer :: p_KldA,p_KldB,p_KldC
    integer, dimension(:), pointer :: p_KrowIdxA
    integer, dimension(:), pointer :: p_KcolA,p_KcolB,p_KcolC
    integer, dimension(:), pointer :: p_KdiagonalB,p_KdiagonalC,p_Kaux
    integer :: h_Kaux,isizeIntl
    logical :: bfast,bmem,bsym,bnum

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => lsyssc_perfconfig
    end if

    ! Set pointer either to the explicitly given destination matrix or
    ! to the second source matrix which serves as destination matrix
    if (present(rdest)) then
      p_rdest => rdest
    else
      p_rdest => rmatrixB
    end if

    ! Set optional arguments
    bmem = .false.; if(present(bmemory)) bmem=bmemory
    bsym = .false.; if(present(bsymb))   bsym=bsymb
    bnum = .false.; if(present(bnumb))   bnum=bnumb

    ! Check if fast implementation can be used
    bfast=.false.; if (present(bisExactStructure)) bfast = bisExactStructure
    bfast = bfast .or. (bmem .and. bsym)

    ! Do we need to generate the destination matrix?
    if (bmem) then
      if (present(rdest)) then
        ! Prepare destination matrix
        call lsyssc_releaseMatrix(rdest)
        if ((rmatrixA%cdataType .eq. ST_DOUBLE) .or.&
            (rmatrixB%cdataType .eq. ST_DOUBLE)) then
          rdest%cdataType = ST_DOUBLE
        else
          rdest%cdataType = ST_SINGLE
        end if
        rdest%NEQ                = rmatrixA%NEQ
        rdest%NCOLS              = rmatrixA%NCOLS
        rdest%bcolumnsSorted     = rmatrixA%bcolumnsSorted
        rdest%browsSorted        = rmatrixA%browsSorted
        rdest%p_rsortStrategyColumns => rmatrixA%p_rsortStrategyColumns
        rdest%p_rsortStrategyRows => rmatrixA%p_rsortStrategyRows
      else
        call output_line('Destination matrix must be given if BMEMORY=TRUE!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
        call sys_halt()
      end if
    end if

    ! Check if all matrices are compatible
    if ((rmatrixA%NEQ   .ne. rmatrixB%NEQ) .or.&
        (rmatrixA%NEQ   .ne. p_rdest%NEQ)  .or.&
        (rmatrixA%NCOLS .ne. rmatrixB%NCOLS) .or.&
        (rmatrixA%NCOLS .ne. p_rdest%NCOLS)) then
      call output_line('Number of rows/columns is not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
      call sys_halt()
    end if

    ! Check if all matrices have the same sorting
    if (rmatrixA%bcolumnsSorted .or. rmatrixB%bcolumnsSorted) then
      if (rmatrixA%p_rsortStrategyColumns%ctype .ne. rmatrixB%p_rsortStrategyColumns%ctype) then
        call output_line("Incompatible sorting strategies!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearComb")
        call sys_halt()
      end if
    end if

    if (rmatrixA%browsSorted .or. rmatrixB%browsSorted) then
      if (rmatrixA%p_rsortStrategyRows%ctype .ne. rmatrixB%p_rsortStrategyRows%ctype) then
        call output_line("Incompatible sorting strategies!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearComb")
        call sys_halt()
      end if
    end if

    if (rmatrixA%bcolumnsSorted .or. p_rdest%bcolumnsSorted) then
      if (rmatrixA%p_rsortStrategyColumns%ctype .ne. p_rdest%p_rsortStrategyColumns%ctype) then
        call output_line("Incompatible sorting strategies!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearComb")
        call sys_halt()
      end if
    end if

    if (rmatrixA%browsSorted .or. p_rdest%browsSorted) then
      if (rmatrixA%p_rsortStrategyRows%ctype .ne. p_rdest%p_rsortStrategyRows%ctype) then
        call output_line("Incompatible sorting strategies!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearComb")
        call sys_halt()
      end if
    end if

    if (rmatrixB%bcolumnsSorted .or. p_rdest%bcolumnsSorted) then
      if (rmatrixB%p_rsortStrategyColumns%ctype .ne. p_rdest%p_rsortStrategyColumns%ctype) then
        call output_line("Incompatible sorting strategies!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearComb")
        call sys_halt()
      end if
    end if

    if (rmatrixB%browsSorted .or. p_rdest%browsSorted) then
      if (rmatrixB%p_rsortStrategyRows%ctype .ne. p_rdest%p_rsortStrategyRows%ctype) then
        call output_line("Incompatible sorting strategies!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearComb")
        call sys_halt()
      end if
    end if

    ! Check if destination matrix has compatible data type
    if (((rmatrixA%cdatatype .eq. ST_DOUBLE) .or.&
         (rmatrixB%cdatatype .eq. ST_DOUBLE)) .and.&
         (p_rdest%cdatatype .ne. ST_DOUBLE)) then
      call output_line('Data type of destination matrix is not compatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
      call sys_halt()
    end if

    ! What is the matrix format of matrix A?
    select case (rmatrixA%cmatrixFormat)

    case (LSYSSC_MATRIX1)   ! A is full matrix ---------------------------------

      select case (rmatrixB%cmatrixFormat)

      case (LSYSSC_MATRIX1) ! B is full matrix  - - - - - - - - - - - - - - - -

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          rdest%cmatrixFormat = LSYSSC_MATRIX1
          rdest%NA            = rmatrixA%NEQ*rmatrixA%NCOLS
          call storage_new('lsyssc_matrixLinearComb','h_Da',&
              rdest%NA,rdest%cdataType,rdest%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (p_rdest%cmatrixFormat .ne. LSYSSC_MATRIX1) then
          call output_line('Destination matrix has incompatible format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolic addition?
        if (bsym .and. present(rdest)) then
          rdest%NEQ   = rmatrixA%NEQ
          rdest%NCOLS = rmatrixA%NCOLS
        end if

        ! numerical addition?
        if (bnum) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,p_DaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_DaB,p_DaC)
                call lalg_vectorLinearComb(p_DaA,p_DaC,ca,cb)
              else
                call lalg_vectorLinearComb(p_DaA,p_DaB,ca,cb)
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_FaB,p_DaC)
                call lalg_vectorLinearComb(p_DaA,p_DaC,ca,cb)
              else
                call output_line('Destination matrix must not be single precision if'//&
                    'one of the source matrices is in double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case (ST_SINGLE) ! A is single precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_DaB,p_DaC)
                call lalg_vectorLinearComb(p_FaA,p_DaC,real(ca,SP),cb)
              else
                call lalg_vectorLinearComb(p_FaA,p_DaB,real(ca,SP),cb)
              end if

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_single(rmatrixB,p_FaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rdest,p_FaC)
                call lalg_copyVector(p_FaB,p_FaC)
                call lalg_vectorLinearComb(p_FaA,p_FaC,real(ca,SP),real(cb,SP))
              else
                call lalg_vectorLinearComb(p_FaA,p_FaB,real(ca,SP),real(cb,SP))
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - - - - - - - -

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          rdest%cmatrixFormat = LSYSSC_MATRIX1
          rdest%NA            = rmatrixA%NEQ*rmatrixA%NCOLS
          call storage_new('lsyssc_matrixLinearComb','h_Da',&
              rdest%NA,rdest%cdataType,rdest%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (p_rdest%cmatrixFormat .ne. LSYSSC_MATRIX1) then
          call output_line('Destination matrix has incompatible format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolical addition?
        if (bsym .and. present(rdest)) then
          rdest%NEQ   = rmatrixA%NEQ
          rdest%NCOLS = rmatrixA%NCOLS
        end if

        ! numerical addition?
        if (bnum) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_double(rmatrixB,p_DaB)
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_DaA,p_DaC)
                call do_matDmat1addDbleDble(rmatrixB%NEQ,rmatrixB%NCOLS,&
                                             p_DaB,p_DaC,cb,ca)
              else
                call output_line('Destination matrix must not be diagonal matrix if'//&
                    'one of the source matrices is a full matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_DaA,p_DaC)
                call do_matDmat1addSnglDble(rmatrixB%NEQ,rmatrixB%NCOLS,&
                                             p_FaB,p_DaC,real(cb,SP),ca)
              else
                call output_line('Destination matrix must not be diagonal matrix if'//&
                    'one of the source matrices is a full matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case (ST_SINGLE) ! A is single precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rmatrixA,p_FaA)
                call lsyssc_getbase_double(rmatrixB,p_DaB)
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_FaA,p_DaC)
                call do_matDmat1addDbleDble(rmatrixB%NEQ,rmatrixB%NCOLS,&
                                             p_DaB,p_DaC,cb,ca)
              else
                call output_line('Destination matrix must not be diagonal matrix if'//&
                    'one of the source matrices is a full matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rmatrixA,p_FaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_single(rdest,p_FaC)
                call lalg_copyVector(p_FaA,p_FaC)
                call do_matDmat1addSnglSngl(rmatrixB%NEQ,rmatrixB%NCOLS,&
                                             p_FaB,p_FaC,real(cb,SP),real(ca,SP))
              else
                call output_line('Destination matrix must not be diagonal matrix if'//&
                    'one of the source matrices is a full matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! B is CSR matrix - - - - - - - - - -

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          rdest%cmatrixFormat = LSYSSC_MATRIX1
          rdest%NA            = rmatrixA%NEQ*rmatrixA%NCOLS
          call storage_new('lsyssc_matrixLinearComb','h_Da',rdest%NA,&
              rdest%cdataType,rdest%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (p_rdest%cmatrixFormat .ne. LSYSSC_MATRIX1) then
          call output_line('Destination matrix has incompatible format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolical addition?
        if (bsym .and. present(rdest)) then
          rdest%NEQ   = rmatrixA%NEQ
          rdest%NCOLS = rmatrixA%NCOLS
        end if

        ! numerical addition?
        if (bnum) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_double(rmatrixB,p_DaB)
                call lsyssc_getbase_double(rdest,p_DaC)
                call lsyssc_getbase_Kld(rmatrixB,p_KldB)
                call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)
                call lalg_copyVector(p_DaA,p_DaC)
                call do_mat79mat1addDbleDble(rmatrixB%NEQ,rmatrixB%NCOLS,&
                                              p_KldB,p_KcolB,p_DaB,p_DaC,cb,ca)
              else
                call output_line('Destination matrix must not be CSR matrix if'//&
                    'one of the source matrices is a full matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_double(rdest,p_DaC)
                call lsyssc_getbase_Kld(rmatrixB,p_KldB)
                call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)
                call lalg_copyVector(p_DaA,p_DaC)
                call do_mat79mat1addSnglDble(rmatrixB%NEQ,rmatrixB%NCOLS,&
                                              p_KldB,p_KcolB,p_FaB,p_DaC,real(cb,SP),ca)
              else
                call output_line('Destination matrix must not be CSR matrix if'//&
                    'one of the source matrices is a full matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case (ST_SINGLE) ! A is single precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rmatrixA,p_FaA)
                call lsyssc_getbase_double(rmatrixB,p_DaB)
                call lsyssc_getbase_double(rdest,p_DaC)
                call lsyssc_getbase_Kld(rmatrixB,p_KldB)
                call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)
                call lalg_copyVector(p_FaA,p_DaC)
                call do_mat79mat1addDbleDble(rmatrixB%NEQ,rmatrixB%NCOLS,&
                                              p_KldB,p_KcolB,p_DaB,p_DaC,cb,ca)
              else
                call output_line('Destination matrix must not be CSR matrix if'//&
                    'one of the source matrices is a full matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rmatrixA,p_FaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_single(rdest,p_FaC)
                call lsyssc_getbase_Kld(rmatrixB,p_KldB)
                call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)
                call lalg_copyVector(p_FaA,p_FaC)
                call do_mat79mat1addSnglSngl(rmatrixB%NEQ,rmatrixB%NCOLS,&
                                              p_KldB,p_KcolB,p_FaB,p_FaC,real(cb,SP),real(ca,SP))
              else
                call output_line('Destination matrix must not be CSR matrix if'//&
                    'one of the source matrices is a full matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case default
        call output_line('Unsupported combination of matrix formats!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
        call sys_halt()
      end select

      !-------------------------------------------------------------------------

    case (LSYSSC_MATRIXD)   ! A is diagonal matrix -----------------------------

      select case (rmatrixB%cmatrixFormat)

      case (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - - - - - - - -

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          rdest%cmatrixFormat = LSYSSC_MATRIXD
          rdest%NA            = rmatrixA%NA
          call storage_new('lsyssc_matrixLinearComb','h_Da',&
              rdest%NA,rdest%cdataType,rdest%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (p_rdest%cmatrixFormat .ne. LSYSSC_MATRIXD) then
          call output_line('Destination matrix has incompatible format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolical addition?
        if (bsym .and. present(rdest)) then
          rdest%NEQ   = rmatrixA%NEQ
          rdest%NCOLS = rmatrixA%NCOLS
        end if

        ! numerical addition?
        if (bnum) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,p_DaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_DaB,p_DaC)
                call lalg_vectorLinearComb(p_DaA,p_DaC,ca,cb)
              else
                call lalg_vectorLinearComb(p_DaA,p_DaB,ca,cb)
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_FaB,p_DaC)
                call lalg_vectorLinearComb(p_DaA,p_DaC,ca,cb)
              else
                call output_line('Destination matrix must not be single precision if'//&
                    'one of the source matrices is in double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case (ST_SINGLE) ! A is single precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_DaB,p_DaC)
                call lalg_vectorLinearComb(p_FaA,p_DaC,real(ca,SP),cb)
              else
                call lalg_vectorLinearComb(p_FaA,p_DaB,real(ca,SP),cb)
              end if

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_single(rmatrixB,p_FaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rdest,p_FaC)
                call lalg_copyVector(p_FaB,p_FaC)
                call lalg_vectorLinearComb(p_FaA,p_FaC,real(ca,SP),real(cb,SP))
              else
                call lalg_vectorLinearComb(p_FaA,p_FaB,real(ca,SP),real(cb,SP))
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case (LSYSSC_MATRIX1) ! B is full matrix - - - - - - - - - - -

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          rdest%cmatrixFormat = LSYSSC_MATRIX1
          rdest%NA            = rmatrixB%NEQ*rmatrixB%NCOLS
          call storage_new('lsyssc_matrixLinearComb','h_Da',&
              rdest%NA,rdest%cdataType,rdest%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (p_rdest%cmatrixFormat .ne. LSYSSC_MATRIX1) then
          call output_line('Destination matrix has incompatible format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolical addition?
        if (bsym .and. present(rdest)) then
          rdest%NEQ   = rmatrixB%NEQ
          rdest%NCOLS = rmatrixB%NCOLS
        end if

        ! numerical addition?
        if (bnum) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,p_DaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_DaB,p_DaC)
                call do_matDmat1addDbleDble(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                             p_DaA,p_DaC,ca,cb)
              else
                call do_matDmat1addDbleDble(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                             p_DaA,p_DaB,ca,cb)
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_FaB,p_DaC)
                call do_matDmat1addDbleDble(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                             p_DaA,p_DaC,ca,cb)
              else
                call output_line('Destination matrix must not be single precision if'//&
                    'one of the source matrices is in double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case (ST_SINGLE) ! A is single precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_DaB,p_DaC)
                call do_matDmat1addSnglDble(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                             p_FaA,p_DaC,real(ca,SP),cb)
              else
                call do_matDmat1addSnglDble(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                             p_FaA,p_DaB,real(ca,SP),cb)
              end if

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_single(rmatrixB,p_FaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rdest,p_FaC)
                call lalg_copyVector(p_FaB,p_FaC)
                call do_matDmat1addSnglSngl(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                             p_FaA,p_FaC,real(ca,SP),real(cb,SP))
              else
                call do_matDmat1addSnglSngl(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                             p_FaA,p_FaB,real(ca,SP),real(cb,SP))
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! B is CSR matrix - - - - - - - - - -

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          rdest%cmatrixFormat = rmatrixB%cmatrixFormat
          rdest%NA            = rmatrixB%NA
          call storage_new('lsyssc_matrixLinearComb','h_Da',&
              rdest%NA,rdest%cdataType,rdest%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (p_rdest%cmatrixFormat .ne. rmatrixB%cmatrixFormat) then
          call output_line('Destination matrix has incompatible format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolic addition?
        if (bsym .and. present(rdest)) then
          call lsyssc_duplicateMatrix(rmatrixB,rdest,LSYSSC_DUP_COPY,LSYSSC_DUP_IGNORE)
        end if

        ! numerical addition?
        if (bnum) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,p_DaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_DaB,p_DaC)
                else
                  call lsyssc_getbase_Kld(rmatrixB,p_KldB)
                  call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addDbleDble(1,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaB,p_DaC,1.0_DP,0.0_DP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if
                call do_matDmat79addDbleDble(1,1,rdest%NEQ,rmatrixB%NA,&
                                              p_KdiagonalC,p_DaA,p_DaC,ca,cb)
              else
                ! Set diagonal pointer
                if (rmatrixB%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                  call lsyssc_getbase_Kld(rmatrixB,p_KdiagonalB)
                else
                  call lsyssc_getbase_Kdiagonal(rmatrixB,p_KdiagonalB)
                end if

                ! We have to apply only the diagonal of matrix A
                call do_matDmat79addDbleDble(1,1,rmatrixB%NEQ,rmatrixB%NA,&
                                              p_KdiagonalB,p_DaA,p_DaB,ca,cb)
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_FaB,p_DaC)
                else
                  call lsyssc_getbase_Kld(rmatrixB,p_KldB)
                  call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addSnglDble(1,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaB,p_DaC,1.0_SP,0.0_DP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix A
                call do_matDmat79addDbleDble(1,1,rdest%NEQ,rdest%NA,&
                                              p_KdiagonalC,p_DaA,p_DaC,ca,cb)
              else
                call output_line('Destination matrix must not be single precision if'//&
                    'one of the source matrices is in double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case (ST_SINGLE) ! A is single precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_DaB,p_DaC)
                else
                  call lsyssc_getbase_Kld(rmatrixB,p_KldB)
                  call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addDbleDble(1,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaB,p_DaC,1.0_DP,0.0_DP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix A
                call do_matDmat79addSnglDble(1,1,rdest%NEQ,rdest%NA,&
                                              p_KdiagonalC,p_FaA,p_DaC,real(ca,SP),cb)
              else
                ! Set diagonal pointer
                if (rmatrixB%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                  call lsyssc_getbase_Kld(rmatrixB,p_KdiagonalB)
                else
                  call lsyssc_getbase_Kdiagonal(rmatrixB,p_KdiagonalB)
                end if

                ! We have to apply only the diagonal of matrix A
                call do_matDmat79addSnglDble(1,1,rmatrixB%NEQ,rmatrixB%NA,&
                                              p_KdiagonalB,p_FaA,p_DaB,real(ca,SP),cb)
              end if

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_single(rmatrixB,p_FaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rdest,p_FaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_FaB,p_FaC)
                else
                  call lsyssc_getbase_Kld(rmatrixB,p_KldB)
                  call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addSnglSngl(1,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaB,p_FaC,1.0_SP,0.0_SP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix A
                call do_matDmat79addSnglSngl(1,1,rdest%NEQ,rdest%NA,&
                                              p_KdiagonalC,p_FaA,p_FaC,real(ca,SP),real(cb,SP))
              else
                ! Set diagonal pointer
                if (rmatrixB%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                  call lsyssc_getbase_Kld(rmatrixB,p_KdiagonalB)
                else
                  call lsyssc_getbase_Kdiagonal(rmatrixB,p_KdiagonalB)
                end if

                ! We have to apply only the diagonal of matrix A
                call do_matDmat79addSnglSngl(1,1,rmatrixB%NEQ,rmatrixB%NA,&
                                              p_KdiagonalB,p_FaA,p_FaB,real(ca,SP),real(cb,SP))
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case (LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9INTL) ! B is interleave CSR matrix

        ! Set size of interleave block
        isizeIntl=rmatrixB%NVAR
        if (rmatrixB%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) isizeIntl=isizeIntl*isizeIntl

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          rdest%cmatrixFormat           = rmatrixB%cmatrixFormat
          rdest%cinterleavematrixFormat = rmatrixB%cinterleavematrixFormat
          rdest%NA                      = rmatrixB%NA
          rdest%NVAR                    = rmatrixB%NVAR
          call storage_new('lsyssc_matrixLinearComb','h_Da',isizeIntl*rdest%NA,&
              rdest%cdataType,rdest%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (p_rdest%cmatrixFormat           .ne. rmatrixB%cmatrixFormat .or. &
            p_rdest%cinterleavematrixFormat .ne. rmatrixB%cinterleavematrixFormat) then
          call output_line('Destination matrix has incompatible format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolic addition?
        if (bsym .and. present(rdest)) then
          call lsyssc_duplicateMatrix(rmatrixB,rdest,LSYSSC_DUP_COPY,LSYSSC_DUP_IGNORE)
        end if

        ! numerical addition?
        if (bnum) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,p_DaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_DaB,p_DaC)
                else
                  call lsyssc_getbase_Kld(rmatrixB,p_KldB)
                  call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addDbleDble(isizeIntl,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaB,p_DaC,1.0_DP,0.0_DP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix A
                if (rdest%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) then
                  call do_matDmat79addDbleDble(rdest%NVAR,rdest%NVAR,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_DaA,p_DaC,ca,cb)
                else
                  call do_matDmat79addDbleDble(rdest%NVAR,1,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_DaA,p_DaC,ca,cb)
                end if
              else
                ! Set diagonal pointer
                if (rmatrixB%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                  call lsyssc_getbase_Kld(rmatrixB,p_KdiagonalB)
                else
                  call lsyssc_getbase_Kdiagonal(rmatrixB,p_KdiagonalB)
                end if

                ! We have to apply only the diagonal of matrix A
                if (rmatrixB%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) then
                  call do_matDmat79addDbleDble(rmatrixB%NVAR,rmatrixB%NVAR,rmatrixB%NEQ,rmatrixB%NA,&
                                                p_KdiagonalB,p_DaA,p_DaB,ca,cb)
                else
                  call do_matDmat79addDbleDble(rmatrixB%NVAR,1,rmatrixB%NEQ,rmatrixB%NA,&
                                                p_KdiagonalB,p_DaA,p_DaB,ca,cb)
                end if
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_FaB,p_DaC)
                else
                  call lsyssc_getbase_Kld(rmatrixB,p_KldB)
                  call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addSnglDble(isizeIntl,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaB,p_DaC,1.0_SP,0.0_DP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix A
                if (rdest%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) then
                  call do_matDmat79addDbleDble(rdest%NVAR,rdest%NVAR,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_DaA,p_DaC,ca,cb)
                else
                  call do_matDmat79addDbleDble(rdest%NVAR,1,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_DaA,p_DaC,ca,cb)
                end if
              else
                call output_line('Destination matrix must not be single precision if'//&
                    'one of the source matrices is in double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case (ST_SINGLE) ! A is single precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_DaB,p_DaC)
                else
                  call lsyssc_getbase_Kld(rmatrixB,p_KldB)
                  call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addDbleDble(isizeIntl,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaB,p_DaC,1.0_DP,0.0_DP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix A
                if (rdest%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) then
                  call do_matDmat79addSnglDble(rdest%NVAR,rdest%NVAR,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_FaA,p_DaC,real(ca,SP),cb)
                else
                  call do_matDmat79addSnglDble(rdest%NVAR,1,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_FaA,p_DaC,real(ca,SP),cb)
                end if
              else
                ! Set diagonal pointer
                if (rmatrixB%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                  call lsyssc_getbase_Kld(rmatrixB,p_KdiagonalB)
                else
                  call lsyssc_getbase_Kdiagonal(rmatrixB,p_KdiagonalB)
                end if

                ! We have to apply only the diagonal of matrix A
                if (rmatrixB%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) then
                  call do_matDmat79addSnglDble(rmatrixB%NVAR,rmatrixB%NVAR,rmatrixB%NEQ,rmatrixB%NA,&
                                                p_KdiagonalB,p_FaA,p_DaB,real(ca,SP),cb)
                else
                  call do_matDmat79addSnglDble(rmatrixB%NVAR,1,rmatrixB%NEQ,rmatrixB%NA,&
                                                p_KdiagonalB,p_FaA,p_DaB,real(ca,SP),cb)
                end if
              end if

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_single(rmatrixB,p_FaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rdest,p_FaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_FaB,p_FaC)
                else
                  call lsyssc_getbase_Kld(rmatrixB,p_KldB)
                  call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addSnglSngl(isizeIntl,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaB,p_FaC,1.0_SP,0.0_SP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix A
                if (rdest%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) then
                  call do_matDmat79addSnglSngl(rdest%NVAR,rdest%NVAR,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_FaA,p_FaC,real(ca,SP),real(cb,SP))
                else
                  call do_matDmat79addSnglSngl(rdest%NVAR,1,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_FaA,p_FaC,real(ca,SP),real(cb,SP))
                end if
              else
                ! Set diagonal pointer
                if (rmatrixB%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                  call lsyssc_getbase_Kld(rmatrixB,p_KdiagonalB)
                else
                  call lsyssc_getbase_Kdiagonal(rmatrixB,p_KdiagonalB)
                end if

                ! We have to apply only the diagonal of matrix A
                if (rmatrixB%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) then
                  call do_matDmat79addSnglSngl(rmatrixB%NVAR,rmatrixB%NVAR,rmatrixB%NEQ,rmatrixB%NA,&
                                                p_KdiagonalB,p_FaA,p_FaB,real(ca,SP),real(cb,SP))
                else
                  call do_matDmat79addSnglSngl(rmatrixB%NVAR,1,rmatrixB%NEQ,rmatrixB%NA,&
                                                p_KdiagonalB,p_FaA,p_FaB,real(ca,SP),real(cb,SP))
                end if
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case default
        call output_line('Unsupported combination of matrix formats!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
        call sys_halt()
      end select

      !-------------------------------------------------------------------------

    case (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! A is CSR matrix ---------------------

      select case (rmatrixB%cmatrixFormat)

      case (LSYSSC_MATRIX1) ! B is full matrix - - - - - - - - - - - - - - - - -

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          rdest%cmatrixFormat = LSYSSC_MATRIX1
          rdest%NA            = rmatrixB%NEQ*rmatrixB%NCOLS
          call storage_new('lsyssc_matrixLinearComb','h_Da',&
              rdest%NA,rdest%cdataType,rdest%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (p_rdest%cmatrixFormat .ne. LSYSSC_MATRIX1) then
          call output_line('Destination matrix has incompatible format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolical addition?
        if (bsym .and. present(rdest)) then
          rdest%NEQ   = rmatrixB%NEQ
          rdest%NCOLS = rmatrixB%NCOLS
        end if

        ! numerical addition?
        if (bnum) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,p_DaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              call lsyssc_getbase_Kld(rmatrixA,p_KldA)
              call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_DaB,p_DaC)
                call do_mat79mat1addDbleDble(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                              p_KldA,p_KcolA,p_DaA,p_DaC,ca,cb)
              else
                call do_mat79mat1addDbleDble(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                              p_KldA,p_KcolA,p_DaA,p_DaB,ca,cb)
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_double(rdest,p_DaC)
                call lsyssc_getbase_Kld(rmatrixA,p_KldA)
                call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
                call lalg_copyVector(p_FaB,p_DaC)
                call do_mat79mat1addDbleDble(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                              p_KldA,p_KcolA,p_DaA,p_DaC,ca,cb)
              else
                call output_line('Destination matrix must not be single precision if'//&
                    'one of the source matrices is in double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case (ST_SINGLE) ! A is single precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              call lsyssc_getbase_Kld(rmatrixA,p_KldA)
              call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)
                call lalg_copyVector(p_DaB,p_DaC)
                call do_mat79mat1addSnglDble(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                              p_KldA,p_KcolA,p_FaA,p_DaC,real(ca,SP),cb)
              else
                call do_mat79mat1addSnglDble(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                              p_KldA,p_KcolA,p_FaA,p_DaB,real(ca,SP),cb)
              end if

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_single(rmatrixB,p_FaB)
              call lsyssc_getbase_Kld(rmatrixA,p_KldA)
              call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rdest,p_FaC)
                call lalg_copyVector(p_FaB,p_FaC)
                call do_mat79mat1addSnglSngl(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                              p_KldA,p_KcolA,p_FaA,p_FaC,real(ca,SP),real(cb,SP))
              else
                call do_mat79mat1addSnglSngl(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                              p_KldA,p_KcolA,p_FaA,p_FaB,real(ca,SP),real(cb,SP))
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - - - - - - - -

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          rdest%cmatrixFormat = rmatrixA%cmatrixFormat
          rdest%NA            = rmatrixA%NA
          call storage_new('lsyssc_matrixLinearComb','h_Da',&
              rdest%NA,rdest%cdataType,rdest%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if (p_rdest%cmatrixFormat .ne. rmatrixA%cmatrixFormat) then
          call output_line('Destination matrix has incompatible format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolic addition?
        if (bsym .and. present(rdest)) then
          call lsyssc_duplicateMatrix(rmatrixA,rdest,LSYSSC_DUP_COPY,LSYSSC_DUP_IGNORE)
        end if

        ! numerical addition?
        if (bnum) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_double(rmatrixB,p_DaB)
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_DaA,p_DaC)
                else
                  call lsyssc_getbase_Kld(rmatrixA,p_KldA)
                  call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := A + 0*C
                  call do_mat79mat79addDbleDble(1,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaA,p_DaC,1.0_DP,0.0_DP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix B
                call do_matDmat79addDbleDble(1,1,rdest%NEQ,rdest%NA,&
                                              p_KdiagonalC,p_DaB,p_DaC,cb,ca)
              else
                call output_line('Destination matrix must not be diagonal matrix if'//&
                    'one of the source matrices is a CSR matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_DaA,p_DaC)
                else
                  call lsyssc_getbase_Kld(rmatrixA,p_KldA)
                  call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagoanl pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := A + 0*C
                  call do_mat79mat79addDbleDble(1,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaA,p_DaC,1.0_DP,0.0_DP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix B
                call do_matDmat79addSnglDble(1,1,rdest%NEQ,rdest%NA,&
                                              p_KdiagonalC,p_FaB,p_DaC,real(cb,SP),ca)
              else
                call output_line('Destination matrix must not be diagonal matrix if'//&
                    'one of the source matrices is a CSR matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case (ST_SINGLE) ! A is single precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rmatrixA,p_FaA)
                call lsyssc_getbase_double(rmatrixB,p_DaB)
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_FaA,p_DaC)
                else
                  call lsyssc_getbase_Kld(rmatrixA,p_KldA)
                  call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := A + 0*C
                  call do_mat79mat79addSnglDble(1,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaA,p_DaC,1.0_SP,0.0_DP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix B
                call do_matDmat79addDbleDble(1,1,rdest%NEQ,rdest%NA,&
                                              p_KdiagonalC,p_DaB,p_DaC,cb,ca)
              else
                call output_line('Destination matrix must not be diagonal matrix if'//&
                    'one of the source matrices is a CSR matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rmatrixA,p_FaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_single(rdest,p_FaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_FaA,p_FaC)
                else
                  call lsyssc_getbase_Kld(rmatrixA,p_KldA)
                  call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := A + 0*C
                  call do_mat79mat79addSnglSngl(1,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaA,p_FaC,1.0_SP,0.0_SP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix B
                call do_matDmat79addSnglSngl(1,1,rdest%NEQ,rdest%NA,&
                                              p_KdiagonalC,p_FaB,p_FaC,real(cb,SP),real(ca,SP))
              else
                call output_line('Destination matrix must not be diagonal matrix if'//&
                    'one of the source matrices is CSR matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! B is CSR matrix - - - - - - - - - -

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixA,p_KldA)
        call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
        call lsyssc_getbase_Kld(rmatrixB,p_KldB)
        call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          ! Duplicate structure of matrix A or B depending on which
          ! matrix is given in the encompassing matrix format
          if (rmatrixA%cmatrixFormat .ge. rmatrixB%cmatrixFormat) then
            call lsyssc_duplicateMatrix(rmatrixA,rdest,LSYSSC_DUP_EMPTY,LSYSSC_DUP_REMOVE)
          else
            call lsyssc_duplicateMatrix(rmatrixB,rdest,LSYSSC_DUP_EMPTY,LSYSSC_DUP_REMOVE)
          end if

          ! Set auxiliary pointers
          h_Kaux = ST_NOHANDLE
          call storage_new('lsyssc_matrixLinearComb','h_Kaux',&
              rmatrixB%NCOLS,ST_INT,h_Kaux,ST_NEWBLOCK_NOINIT)
          call storage_getbase_int(h_Kaux,p_Kaux)

          ! Compute number of nonzero matrix entries: NA
          rdest%NA=do_mat79mat79add_computeNA(rmatrixA%NEQ,p_KldA,p_KcolA,p_KldB,p_KcolB,p_Kaux)
          call storage_new('lsyssc_matrixLinearComb','h_Da',&
              rdest%NA,rdest%cdataType,rdest%h_Da,ST_NEWBLOCK_NOINIT)
          call storage_new('lsyssc_matrixLinearComb','h_Kcol',&
              rdest%NA,ST_INT,rdest%h_Kcol,ST_NEWBLOCK_NOINIT)

          ! Free temporal memory
          call storage_free(h_Kaux)
        end if

        ! Check if matrix is given in the correct format
        if (p_rdest%cmatrixFormat .ne.&
            max(rmatrixA%cmatrixFormat,rmatrixB%cmatrixFormat)) then
          call output_line('Destination matrix has incompatible format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolic addition?
        if (bsym .and. present(rdest)) then

          ! Set pointers
          call lsyssc_getbase_Kld(rdest,p_KldC)
          call lsyssc_getbase_Kcol(rdest,p_KcolC)

          if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX9) then
            call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
            call do_mat79mat79add_symb(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                       p_KldA,p_KcolA,rmatrixA%cmatrixFormat,&
                                       p_KldB,p_KcolB,rmatrixB%cmatrixFormat,&
                                       p_KldC,p_KcolC,p_KdiagonalC)
          else
            call do_mat79mat79add_symb(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                       p_KldA,p_KcolA,rmatrixA%cmatrixFormat,&
                                       p_KldB,p_KcolB,rmatrixB%cmatrixFormat,&
                                       p_KldC,p_KcolC)
          end if
        end if

        ! numerical addition?
        if (bnum) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,p_DaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. &
                    (rmatrixA%cmatrixFormat .eq. rdest%cmatrixFormat) .and. (bfast)) then

                  call lalg_copyVector(p_DaB, p_DaC)
                  call lalg_vectorLinearComb(p_DaA,p_DaC,ca,cb)
                else
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagoanl pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addDbleDble(1,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaB,p_DaC,1.0_DP,0.0_DP)

                  ! Compute C := ca*A + cb*C
                  call do_mat79mat79addDbleDble(1,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaA,p_DaC,ca,cb)
                end if
              else
                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. (bfast)) then

                  call lalg_vectorLinearComb(p_DaA,p_DaB,ca,cb)
                else
                  ! Set diagonal pointer
                  if (rmatrixB%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                  call lsyssc_getbase_Kld(rmatrixB,p_KdiagonalB)
                else
                  call lsyssc_getbase_Kdiagonal(rmatrixB,p_KdiagonalB)
                end if

                  ! Compute B:= ca*A + cb*B
                  call do_mat79mat79addDbleDble(1,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,p_KdiagonalB,&
                                                 p_DaA,p_DaB,ca,cb)
                end if
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. &
                    (rmatrixA%cmatrixFormat .eq. rdest%cmatrixFormat) .and. (bfast)) then

                  call lalg_copyVector(p_FaB, p_DaC)
                  call lalg_vectorLinearComb(p_DaA,p_DaC,ca,cb)
                else
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addSnglDble(1,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaB,p_DaC,1.0_SP,0.0_DP)

                  ! Compute C := ca*A + cb*C
                  call do_mat79mat79addDbleDble(1,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaA,p_DaC,ca,cb)
                end if
              else
                call output_line('Destination matrix must not be single precision if'//&
                    'one of the source matrices is in double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case (ST_SINGLE) ! A is single precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. &
                    (rmatrixA%cmatrixFormat .eq. rdest%cmatrixFormat) .and. (bfast)) then

                  call lalg_copyVector(p_DaB, p_DaC)
                  call lalg_vectorLinearComb(p_FaA,p_DaC,real(ca,SP),cb)
                else
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addDbleDble(1,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaB,p_DaC,1.0_DP,0.0_DP)

                  ! Compute C := ca*A + cb*C
                  call do_mat79mat79addSnglDble(1,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaA,p_DaC,real(ca,SP),cb)
                end if
              else
                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. (bfast)) then

                  call lalg_vectorLinearComb(p_FaA,p_DaB,real(ca,SP),cb)
                else
                  ! Set diagonal pointer
                  if (rmatrixB%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rmatrixB,p_KdiagonalB)
                  else
                    call lsyssc_getbase_Kdiagonal(rmatrixB,p_KdiagonalB)
                  end if

                  ! Compute B:= ca*A + cb*B
                  call do_mat79mat79addSnglDble(1,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,p_KdiagonalB,&
                                                 p_FaA,p_DaB,real(ca,SP),cb)
                end if
              end if

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_single(rmatrixB,p_FaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rdest,p_FaC)

                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. &
                    (rmatrixA%cmatrixFormat .eq. rdest%cmatrixFormat) .and. (bfast)) then

                  call lalg_copyVector(p_FaB, p_FaC)
                  call lalg_vectorLinearComb(p_FaA,p_FaC,real(ca,SP),real(cb,SP))
                else
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addSnglSngl(1,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaB,p_FaC,1.0_SP,0.0_SP)

                  ! Compute C := ca*A + cb*C
                  call do_mat79mat79addSnglSngl(1,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaA,p_FaC,real(ca,SP),real(cb,SP))
                end if
              else
                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. (bfast)) then

                  call lalg_vectorLinearComb(p_FaA,p_FaB,real(ca,SP),real(cb,SP))
                else
                  ! Set diagonal pointer
                  if (rmatrixB%cmatrixFormat .eq. LSYSSC_MATRIX7) then
                    call lsyssc_getbase_Kld(rmatrixB,p_KdiagonalB)
                  else
                    call lsyssc_getbase_Kdiagonal(rmatrixB,p_KdiagonalB)
                  end if

                  ! Compute B:= ca*A + cb*B
                  call do_mat79mat79addSnglSngl(1,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,p_KdiagonalB,&
                                                 p_FaA,p_FaB,real(ca,SP),real(cb,SP))
                end if
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case default
        call output_line('Unsupported combination of matrix formats!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
        call sys_halt()
      end select

      !-------------------------------------------------------------------------

    case (LSYSSC_MATRIX9ROWC) ! A is differential CSR matrix -------------------

      select case (rmatrixB%cmatrixFormat)

      case (LSYSSC_MATRIX7,LSYSSC_MATRIX9) ! B is CSR matrix - - - - - - - - - -

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixA,p_KldA)
        call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
        call lsyssc_getbase_Kld(rmatrixB,p_KldB)
        call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          call output_line('Not supported!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolic addition?
        if (bsym .and. present(rdest)) then

          call output_line('Not supported!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()

        end if

        ! numerical addition?
        if (bnum) then

          ! Get row indices
          call lsyssc_getbase_KrowIdx(rdest,p_KrowIdxA)

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,p_DaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)

              ! Only 'fast' mode possible at the moment!
              if (.not. bfast) then
                call output_line('Not supported!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)

                call lsyssc_getbase_Kld(rdest,p_KldC)
                call lsyssc_getbase_Kcol(rdest,p_KcolC)

                ! Compute C := B + 0*C
                call lalg_copyVector(p_DaB, p_DaC)

                ! Compute C := ca*A + cb*C
                call do_mat9rowcMat9addDbleDble(rmatrixA%NNZROWS,&
                    rmatrixA%NA,p_KldA,p_KcolA,p_KrowIdxA,&
                    rdest%NA,p_KldC,p_KcolC,&
                    p_DaA,p_DaC,ca,cb)
              else

                ! Compute B:= ca*A + cb*B
                call do_mat9rowcMat9addDbleDble(rmatrixA%NNZROWS,&
                    rmatrixA%NA,p_KldA,p_KcolA,p_KrowIdxA,&
                    rmatrixB%NA,p_KldB,p_KcolB,&
                    p_DaA,p_DaB,ca,cb)

              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case default
        call output_line('Unsupported combination of matrix formats!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
        call sys_halt()
      end select

      !-------------------------------------------------------------------------

    case (LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9INTL) ! A is interleave CSR matrix -

      ! Set size of interleave block
      isizeIntl=rmatrixB%NVAR
      if (rmatrixB%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) isizeIntl=isizeIntl*isizeIntl

      select case (rmatrixB%cmatrixFormat)

      case (LSYSSC_MATRIXD) ! B is diagonal matrix - - - - - - - - - - - - - - -

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          rdest%cmatrixFormat           = rmatrixA%cmatrixFormat
          rdest%cinterleavematrixFormat = rmatrixA%cinterleavematrixFormat
          rdest%NA                      = rmatrixA%NA
          rdest%NVAR                    = rmatrixA%NVAR
          call storage_new('lsyssc_matrixLinearComb','h_Da',isizeIntl*rdest%NA,&
              rdest%cdataType,rdest%h_Da,ST_NEWBLOCK_NOINIT)
        end if

        ! Check if matrix is given in the correct format
        if ((p_rdest%cmatrixFormat .ne. rmatrixA%cmatrixFormat) .or. &
            (p_rdest%cinterleavematrixFormat .ne. rmatrixA%cinterleavematrixFormat)) then
          call output_line('Destination matrix has incompatible format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolic addition?
        if (bsym .and. present(rdest)) then
          call lsyssc_duplicateMatrix(rmatrixA,rdest,LSYSSC_DUP_COPY,LSYSSC_DUP_IGNORE)
        end if

        ! numerical addition?
        if (bnum) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_double(rmatrixB,p_DaB)
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_DaA,p_DaC)
                else
                  call lsyssc_getbase_Kld(rmatrixA,p_KldA)
                  call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := A + 0*C
                  call do_mat79mat79addDbleDble(isizeIntl,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaA,p_DaC,1.0_DP,0.0_DP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix B
                if (rdest%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) then
                  call do_matDmat79addDbleDble(rdest%NVAR,rdest%NVAR,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_DaB,p_DaC,cb,ca)
                else
                  call do_matDmat79addDbleDble(rdest%NVAR,1,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_DaB,p_DaC,cb,ca)
                end if
              else
                call output_line('Destination matrix must not be diagonal matrix if'//&
                    'one of the source matrices is a CSR matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_DaA,p_DaC)
                else
                  call lsyssc_getbase_Kld(rmatrixA,p_KldA)
                  call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := A + 0*C
                  call do_mat79mat79addDbleDble(isizeIntl,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaA,p_DaC,1.0_DP,0.0_DP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix B
                if (rdest%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) then
                  call do_matDmat79addSnglDble(rdest%NVAR,rdest%NVAR,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_FaB,p_DaC,real(cb,SP),ca)
                else
                  call do_matDmat79addSnglDble(rdest%NVAR,1,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_FaB,p_DaC,real(cb,SP),ca)
                end if
              else
                call output_line('Destination matrix must not be diagonal matrix if'//&
                    'one of the source matrices is a CSR matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case (ST_SINGLE) ! A is single precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rmatrixA,p_FaA)
                call lsyssc_getbase_double(rmatrixB,p_DaB)
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_FaA,p_DaC)
                else
                  call lsyssc_getbase_Kld(rmatrixA,p_KldA)
                  call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := A + 0*C
                  call do_mat79mat79addSnglDble(isizeIntl,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaA,p_DaC,1.0_SP,0.0_DP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix B
                if (rdest%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) then
                  call do_matDmat79addDbleDble(rdest%NVAR,rdest%NVAR,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_DaB,p_DaC,cb,ca)
                else
                  call do_matDmat79addDbleDble(rdest%NVAR,1,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_DaB,p_DaC,cb,ca)
                end if
              else
                call output_line('Destination matrix must not be diagonal matrix if'//&
                    'one of the source matrices is a CSR matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rmatrixA,p_FaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_single(rdest,p_FaC)

                ! Can we make use of fast copy operation?
                if (bfast) then
                  call lalg_copyVector(p_FaA,p_FaC)
                else
                  call lsyssc_getbase_Kld(rmatrixA,p_KldA)
                  call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := A + 0*C
                  call do_mat79mat79addSnglSngl(isizeIntl,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaA,p_FaC,1.0_SP,0.0_SP)
                end if

                ! Set diagonal pointer
                if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                  call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                else
                  call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                end if

                ! We have to apply only the diagonal of matrix B
                if (rdest%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) then
                  call do_matDmat79addSnglSngl(rdest%NVAR,rdest%NVAR,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_FaB,p_FaC,real(cb,SP),real(ca,SP))
                else
                  call do_matDmat79addSnglSngl(rdest%NVAR,1,rdest%NEQ,rdest%NA,&
                                                p_KdiagonalC,p_FaB,p_FaC,real(cb,SP),real(ca,SP))
                end if
              else
                call output_line('Destination matrix must not be diagonal matrix if'//&
                    'one of the source matrices is CSR matrix!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case (LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9INTL) ! B is interleave CSR matrix-

        ! Check if interleave matrices are compatible
        if (rmatrixA%NVAR .ne. rmatrixB% NVAR .or. &
            rmatrixA%cinterleavematrixFormat .ne. rmatrixB%cinterleavematrixFormat) then
          call output_line('Source matrices have incompatible interleave format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! Set size of interleave block
        isizeIntl=rmatrixA%NVAR
        if (rmatrixA%cinterleavematrixFormat .eq. LSYSSC_MATRIX1) isizeIntl=isizeIntl*isizeIntl

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrixA,p_KldA)
        call lsyssc_getbase_Kcol(rmatrixA,p_KcolA)
        call lsyssc_getbase_Kld(rmatrixB,p_KldB)
        call lsyssc_getbase_Kcol(rmatrixB,p_KcolB)

        ! memory allocation?
        if (bmem .and. present(rdest)) then
          ! Duplicate structure of matrix A or B depending on which
          ! matrix is given in the encompassing matrix format
          if (rmatrixA%cmatrixFormat .ge. rmatrixB%cmatrixFormat) then
            call lsyssc_duplicateMatrix(rmatrixA,rdest,LSYSSC_DUP_EMPTY,LSYSSC_DUP_REMOVE)
          else
            call lsyssc_duplicateMatrix(rmatrixB,rdest,LSYSSC_DUP_EMPTY,LSYSSC_DUP_REMOVE)
          end if

          ! Set auxiliary pointers
          h_Kaux = ST_NOHANDLE
          call storage_new('lsyssc_matrixLinearComb','h_Kaux',&
              rmatrixB%NCOLS,ST_INT,h_Kaux,ST_NEWBLOCK_NOINIT)
          call storage_getbase_int(h_Kaux,p_Kaux)

          ! Compute number of nonzero matrix entries: NA
          rdest%NA=do_mat79mat79add_computeNA(rmatrixA%NEQ,p_KldA,p_KcolA,p_KldB,p_KcolB,p_Kaux)
          call storage_new('lsyssc_matrixLinearComb','h_Da',isizeIntl*rdest%NA,&
              rdest%cdataType,rdest%h_Da,ST_NEWBLOCK_NOINIT)
          call storage_realloc('lsyssc_matrixLinearComb',rdest%NA,&
              rdest%h_Kcol,ST_NEWBLOCK_NOINIT,.false.)

          ! Free temporal memory
          call storage_free(h_Kaux)
        end if

        ! Check if matrix is given in the correct format
        if (p_rdest%cmatrixFormat .ne.&
              max(rmatrixA%cmatrixFormat,rmatrixB%cmatrixFormat) .or. &
            p_rdest%cinterleavematrixFormat .ne.&
              rmatrixA%cinterleavematrixFormat .or. &
            p_rdest%NA .ne. rmatrixA%NVAR) then
          call output_line('Destination matrix has incompatible format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
          call sys_halt()
        end if

        ! symbolic addition?
        if (bsym .and. present(rdest)) then

          ! Set pointers
          call lsyssc_getbase_Kld(rdest,p_KldC)
          call lsyssc_getbase_Kcol(rdest,p_KcolC)

          if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX9INTL) then
            call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
            call do_mat79mat79add_symb(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                       p_KldA,p_KcolA,rmatrixA%cmatrixFormat,&
                                       p_KldB,p_KcolB,rmatrixB%cmatrixFormat,&
                                       p_KldC,p_KcolC,p_KdiagonalC)
          else
            call do_mat79mat79add_symb(rmatrixA%NEQ,rmatrixA%NCOLS,&
                                       p_KldA,p_KcolA,rmatrixA%cmatrixFormat,&
                                       p_KldB,p_KcolB,rmatrixB%cmatrixFormat,&
                                       p_KldC,p_KcolC)
          end if
        end if

        ! numerical addition?
        if (bnum) then

          ! Find the correct internal subroutine for the specified
          ! data types. Note that the resulting matrix C will be
          ! double if at least one of the source matrices is of type
          ! double
          select case (rmatrixA%cdataType)

          case (ST_DOUBLE) ! A is double precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_double(rmatrixA,p_DaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. &
                    (rmatrixA%cmatrixFormat .eq. rdest%cmatrixFormat) .and. (bfast)) then

                  call lalg_copyVector(p_DaB, p_DaC)
                  call lalg_vectorLinearComb(p_DaA,p_DaC,ca,cb)
                else
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addDbleDble(isizeIntl,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaB,p_DaC,1.0_DP,0.0_DP)

                  ! Compute C := ca*A + cb*C
                  call do_mat79mat79addDbleDble(isizeIntl,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaA,p_DaC,ca,cb)
                end if
              else
                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. (bfast)) then

                  call lalg_vectorLinearComb(p_DaA,p_DaB,ca,cb)
                else
                  ! Set diagonal pointer
                  if (rmatrixB%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rmatrixB,p_KdiagonalB)
                  else
                    call lsyssc_getbase_Kdiagonal(rmatrixB,p_KdiagonalB)
                  end if

                  ! Compute B:= ca*A + cb*B
                  call do_mat79mat79addDbleDble(isizeIntl,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,p_KdiagonalB,&
                                                 p_DaA,p_DaB,ca,cb)
                end if
              end if

            case (ST_SINGLE)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rmatrixA,p_DaA)
                call lsyssc_getbase_single(rmatrixB,p_FaB)
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. &
                    (rmatrixA%cmatrixFormat .eq. rdest%cmatrixFormat) .and. (bfast)) then

                  call lalg_copyVector(p_FaB, p_DaC)
                  call lalg_vectorLinearComb(p_DaA,p_DaC,ca,cb)
                else
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addSnglDble(isizeIntl,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaB,p_DaC,1.0_SP,0.0_DP)

                  ! Compute C := ca*A + cb*C
                  call do_mat79mat79addDbleDble(1,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaA,p_DaC,ca,cb)
                end if
              else
                call output_line('Destination matrix must not be single precision if'//&
                    'one of the source matrices is in double precision!',&
                    OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
                call sys_halt()
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case (ST_SINGLE) ! A is single precision matrix ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            select case (rmatrixB%cdataType)

            case (ST_DOUBLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_double(rmatrixB,p_DaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_double(rdest,p_DaC)

                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. &
                    (rmatrixA%cmatrixFormat .eq. rdest%cmatrixFormat) .and. (bfast)) then

                  call lalg_copyVector(p_DaB, p_DaC)
                  call lalg_vectorLinearComb(p_FaA,p_DaC,real(ca,SP),cb)
                else
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addDbleDble(isizeIntl,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_DaB,p_DaC,1.0_DP,0.0_DP)

                  ! Compute C := ca*A + cb*C
                  call do_mat79mat79addSnglDble(isizeIntl,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaA,p_DaC,real(ca,SP),cb)
                end if
              else
                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. (bfast)) then

                  call lalg_vectorLinearComb(p_FaA,p_DaB,real(ca,SP),cb)
                else
                  ! Set diagonal pointer
                  if (rmatrixB%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rmatrixB,p_KdiagonalB)
                  else
                    call lsyssc_getbase_Kdiagonal(rmatrixB,p_KdiagonalB)
                  end if

                  ! Compute B:= ca*A + cb*B
                  call do_mat79mat79addSnglDble(isizeIntl,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,p_KdiagonalB,&
                                                 p_FaA,p_DaB,real(ca,SP),cb)
                end if
              end if

            case (ST_SINGLE)
              call lsyssc_getbase_single(rmatrixA,p_FaA)
              call lsyssc_getbase_single(rmatrixB,p_FaB)
              ! Do we have an explicit destination matrix?
              if (present(rdest)) then
                call lsyssc_getbase_single(rdest,p_FaC)

                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. &
                    (rmatrixA%cmatrixFormat .eq. rdest%cmatrixFormat) .and. (bfast)) then

                  call lalg_copyVector(p_FaB, p_FaC)
                  call lalg_vectorLinearComb(p_FaA,p_FaC,real(ca,SP),real(cb,SP))
                else
                  call lsyssc_getbase_Kld(rdest,p_KldC)
                  call lsyssc_getbase_Kcol(rdest,p_KcolC)

                  ! Set diagonal pointer
                  if (rdest%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rdest,p_KdiagonalC)
                  else
                    call lsyssc_getbase_Kdiagonal(rdest,p_KdiagonalC)
                  end if

                  ! Compute C := B + 0*C
                  call do_mat79mat79addSnglSngl(isizeIntl,rmatrixB%NEQ,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaB,p_FaC,1.0_SP,0.0_SP)

                  ! Compute C := ca*A + cb*C
                  call do_mat79mat79addSnglSngl(isizeIntl,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rdest%NA,p_KldC,p_KcolC,p_KdiagonalC,&
                                                 p_FaA,p_FaC,real(ca,SP),real(cb,SP))
                end if
              else
                ! Can we make use of fast copy operation?
                if ((rmatrixA%cmatrixFormat .eq. rmatrixB%cmatrixFormat) .and. (bfast)) then

                  call lalg_vectorLinearComb(p_FaA,p_FaB,real(ca,SP),real(cb,SP))
                else
                  ! Set diagonal pointer
                  if (rmatrixB%cmatrixFormat .eq. LSYSSC_MATRIX7INTL) then
                    call lsyssc_getbase_Kld(rmatrixB,p_KdiagonalB)
                  else
                    call lsyssc_getbase_Kdiagonal(rmatrixB,p_KdiagonalB)
                  end if

                  ! Compute B:= ca*A + cb*B
                  call do_mat79mat79addSnglSngl(isizeIntl,rmatrixA%NEQ,&
                                                 rmatrixA%NA,p_KldA,p_KcolA,&
                                                 rmatrixB%NA,p_KldB,p_KcolB,p_KdiagonalB,&
                                                 p_FaA,p_FaB,real(ca,SP),real(cb,SP))
                end if
              end if

            case default
              call output_line('Unsupported combination of data types!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
              call sys_halt()
            end select

          case default
            call output_line('Unsupported data type of source matrix!',&
                OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
            call sys_halt()
          end select
        end if

      case default
        call output_line('Unsupported combination of matrix formats!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
        call sys_halt()
      end select

      !-------------------------------------------------------------------------

    case default
      call output_line('Unsupported matrix format of source matrix!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_matrixLinearComb')
      call sys_halt()
    end select

  contains

    ! Here, the real MM addition routines follow.

    !**************************************************************
    ! Format D-1 addition:    B := ca*A + cb*B
    ! double precision matrix A (format D)
    ! double precision matrix B (format 1)

    subroutine do_matDmat1addDbleDble(neq,ncols,DaA,DaB,da,db)

      integer, intent(in) :: neq,ncols
      real(DP), intent(in) :: da,db
      real(DP), dimension(neq), intent(in) :: DaA
      real(DP), dimension(ncols,neq), intent(inout) :: DaB

      integer :: ieq

      ! Compute B := cb*B
      call lalg_scaleVector(DaB,db)

      ! Compute B := B + diag(A)
      do ieq = 1,neq
        DaB(ieq,ieq) = DaB(ieq,ieq)+da*DaA(ieq)
      end do

    end subroutine

    !**************************************************************
    ! Format D-1 addition:    B := ca*A + cb*B
    ! single precision matrix A (format D)
    ! single precision matrix B (format 1)

    subroutine do_matDmat1addSnglSngl(neq,ncols,FaA,FaB,fa,fb)

      integer, intent(in) :: neq,ncols
      real(SP), intent(in) :: fa,fb
      real(SP), dimension(neq), intent(in) :: FaA
      real(SP), dimension(ncols,neq), intent(inout) :: FaB

      integer :: ieq

      ! Compute B := cb*B
      call lalg_scaleVector(FaB,fb)

      ! Compute B := B + diag(A)
      do ieq = 1,neq
        FaB(ieq,ieq) = FaB(ieq,ieq)+fa*FaA(ieq)
      end do

    end subroutine

    !**************************************************************
    ! Format D-1 addition:    B := ca*A + cb*B
    ! single precision matrix A (format D)
    ! double precision matrix B (format 1)

    subroutine do_matDmat1addSnglDble(neq,ncols,FaA,DaB,fa,db)

      integer, intent(in) :: neq,ncols
      real(SP), intent(in) :: fa
      real(DP), intent(in) :: db
      real(SP), dimension(neq), intent(in) :: FaA
      real(DP), dimension(ncols,neq), intent(inout) :: DaB

      integer :: ieq

      ! Compute B := cb*B
      call lalg_scaleVector(DaB,db)

      ! Compute B := B + diag(A)
      do ieq = 1,neq
        DaB(ieq,ieq) = DaB(ieq,ieq)+fa*FaA(ieq)
      end do

    end subroutine

    !**************************************************************
    ! Format 7/9-1 addition:  B := ca*A + cb*B
    ! double precision matrix A (format 7 or format 9)
    ! double precision matrix B (format 1)
    subroutine do_mat79mat1addDbleDble(neq,ncols,Kld,Kcol,DaA,DaB,da,db)

      ! REMARK: The matrix B (format 1) is stored row-wise. Hence,
      ! this subroutine handles the transposed of matrix B. Therefore,
      ! the row and column indices IEQ and ICOL are swapped when the
      ! contribution of matrix A is applied.

      integer, intent(in) :: neq,ncols
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(DP), intent(in) :: da,db
      real(DP), dimension(:), intent(in) :: DaA
      real(DP), dimension(ncols,neq), intent(inout) :: DaB

      integer :: ild,ieq,icol

      ! Compute B := cb*B
      call lalg_scaleVector(DaB,db)

      ! Compute B := B + A
      !$omp parallel do private(ild,icol) default(shared) &
      !$omp if(NEQ > p_rperfconfig%NEQMIN_OMP)
      do ieq = 1,neq
        do ild = Kld(ieq),Kld(ieq+1)-1
          icol = Kcol(ild)
          DaB(icol,ieq) = DaB(icol,ieq)+da*DaA(ild)
        end do
      end do
      !$omp end parallel do

    end subroutine

    !**************************************************************
    ! Format 7/9-1 addition:  B := ca*A + cb*B
    ! single precision matrix A (format 7 or format 9)
    ! single precision matrix B (format 1)
    subroutine do_mat79mat1addSnglSngl(neq,ncols,Kld,Kcol,FaA,FaB,fa,fb)

      ! REMARK: The matrix B (format 1) is stored row-wise. Hence,
      ! this subroutine handles the transposed of matrix B. Therefore,
      ! the row and column indices IEQ and ICOL are swapped when the
      ! contribution of matrix A is applied.

      integer, intent(in) :: neq,ncols
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(SP), intent(in) :: fa,fb
      real(SP), dimension(:), intent(in) :: FaA
      real(SP), dimension(ncols,neq), intent(inout) :: FaB

      integer :: ild,ieq,icol

      ! Compute B := cb*B
      call lalg_scaleVector(FaB,fb)

      ! Compute B := B + A
      !$omp parallel do private(ild,icol) default(shared) &
      !$omp if(NEQ > p_rperfconfig%NEQMIN_OMP)
      do ieq = 1,neq
        do ild = Kld(ieq),Kld(ieq+1)-1
          icol = Kcol(ild)
          FaB(icol,ieq) = FaB(icol,ieq)+fa*FaA(ild)
        end do
      end do
      !$omp end parallel do

    end subroutine

    !**************************************************************
    ! Format 7/9-1 addition:  B := ca*A + cb*B
    ! single precision matrix A (format 7 or format 9)
    ! double precision matrix B (format 1)
    subroutine do_mat79mat1addSnglDble(neq,ncols,Kld,Kcol,FaA,DaB,fa,db)

      ! REMARK: The matrix B (format 1) is stored row-wise. Hence,
      ! this subroutine handles the transposed of matrix B. Therefore,
      ! the row and column indices IEQ and ICOL are swapped when the
      ! contribution of matrix A is applied.

      integer, intent(in) :: neq,ncols
      integer, dimension(:), intent(in) :: Kld
      integer, dimension(:), intent(in) :: Kcol
      real(SP), intent(in) :: fa
      real(DP), intent(in) :: db
      real(SP), dimension(:), intent(in) :: FaA
      real(DP), dimension(ncols,neq), intent(inout) :: DaB

      integer :: ild,ieq,icol

      ! Compute B := cb*B
      call lalg_scaleVector(DaB,db)

      ! Compute B := B + A
      !$omp parallel do private(ild,icol) default(shared) &
      !$omp if(NEQ > p_rperfconfig%NEQMIN_OMP)
      do ieq = 1,neq
        do ild = Kld(ieq),Kld(ieq+1)-1
          icol = Kcol(ild)
          DaB(icol,ieq) = DaB(icol,ieq)+fa*FaA(ild)
        end do
      end do
      !$omp end parallel do

    end subroutine

    !**************************************************************
    ! Format 7/9-D addition:  B := ca*A + cb*B
    ! double precision matrix A (format D)
    ! double precision matrix B (format 7 or format 9, interleave possible)

    subroutine do_matDmat79addDbleDble(nvar,mvar,neq,na,Kdiagonal,DaA,DaB,da,db)

      integer, intent(in) :: neq,na
      integer, intent(in) :: nvar,mvar
      integer, dimension(:), intent(in) :: Kdiagonal
      real(DP), intent(in) :: da,db
      real(DP), dimension(:), intent(in) :: DaA
      real(DP), dimension(nvar,mvar,na), intent(inout) :: DaB

      integer :: ieq,ild,ivar

      ! Compute B := cb*B
      call lalg_scaleVector(DaB,db)

      if (nvar .eq. mvar) then

        ! Compute B := B + A
        !$omp parallel do private(ild,ivar) default(shared) &
        !$omp if(NEQ > p_rperfconfig%NEQMIN_OMP)
        do ieq = 1,neq
          ild = Kdiagonal(ieq)
          do ivar = 1,nvar
            DaB(ivar,ivar,ild) = DaB(ivar,ivar,ild)+da*DaA(ieq)
          end do
        end do
        !$omp end parallel do

      elseif (nvar .eq. 1) then

        ! Compute B := B + A
        !$omp parallel do private(ild) default(shared) &
        !$omp if(NEQ > p_rperfconfig%NEQMIN_OMP)
        do ieq = 1,neq
          ild = Kdiagonal(ieq)
          DaB(1,:,ild) = DaB(1,:,ild)+da*DaA(ieq)
        end do
        !$omp end parallel do

      elseif (mvar .eq. 1) then

        ! Compute B := B + A
        !$omp parallel do private(ild) default(shared) &
        !$omp if(NEQ > p_rperfconfig%NEQMIN_OMP)
        do ieq = 1,neq
          ild = Kdiagonal(ieq)
          DaB(:,1,ild) = DaB(:,1,ild)+da*DaA(ieq)
        end do
        !$omp end parallel do

      else   ! nvar /= mvar

        call output_line('Destination matrix must satisfy NVAR = MVAR!',&
            OU_CLASS_ERROR,OU_MODE_STD,'do_matDmat79DbleDble')
        call sys_halt()

      end if

    end subroutine

    !**************************************************************
    ! Format 7/9-D addition:  B := ca*A + cb*B
    ! single precision matrix A (format D)
    ! single precision matrix B (format 7 or format 9, interleave possible)

    subroutine do_matDmat79addSnglSngl(nvar,mvar,neq,na,Kdiagonal,FaA,FaB,fa,fb)

      integer, intent(in) :: neq,na
      integer, intent(in) :: nvar,mvar
      integer, dimension(:), intent(in) :: Kdiagonal
      real(SP), intent(in) :: fa,fb
      real(SP), dimension(:), intent(in) :: FaA
      real(SP), dimension(nvar,mvar,na), intent(inout) :: FaB

      integer :: ieq,ild,ivar

      ! Compute B := cb*B
      call lalg_scaleVector(FaB,fb)

      if (nvar .eq. mvar) then

        ! Compute B := B + A
        !$omp parallel do private(ild,ivar) default(shared) &
        !$omp if(NEQ > p_rperfconfig%NEQMIN_OMP)
        do ieq = 1,neq
          ild = Kdiagonal(ieq)
          do ivar = 1,nvar
            FaB(ivar,ivar,ild) = FaB(ivar,ivar,ild)+fa*FaA(ieq)
          end do
        end do
        !$omp end parallel do

      elseif (nvar .eq. 1) then

        ! Compute B := B + A
        !$omp parallel do private(ild) default(shared) &
        !$omp if(NEQ > p_rperfconfig%NEQMIN_OMP)
        do ieq = 1,neq
          ild = Kdiagonal(ieq)
          FaB(1,:,ild) = FaB(1,:,ild)+fa*FaA(ieq)
        end do
        !$omp end parallel do

      elseif (mvar .eq. 1) then

        ! Compute B := B + A
        !$omp parallel do private(ild) default(shared) &
        !$omp if(NEQ > p_rperfconfig%NEQMIN_OMP)
        do ieq = 1,neq
          ild = Kdiagonal(ieq)
          FaB(:,1,ild) = FaB(:,1,ild)+fa*FaA(ieq)
        end do
        !$omp end parallel do

      else   ! nvar /= mvar

        call output_line('Destination matrix must satisfy NVAR = MVAR!',&
            OU_CLASS_ERROR,OU_MODE_STD,'do_matDmat79SnglSngl')
        call sys_halt()

      end if

    end subroutine

    !**************************************************************
    ! Format 7/9-D addition:  B := ca*A + cb*B
    ! single precision matrix A (format D)
    ! double precision matrix B (format 7 or format 9, interleave possible)

    subroutine do_matDmat79addSnglDble(nvar,mvar,neq,na,Kdiagonal,FaA,DaB,fa,db)

      integer, intent(in) :: neq,na
      integer, intent(in) :: nvar,mvar
      integer, dimension(:), intent(in) :: Kdiagonal
      real(SP), intent(in) :: fa
      real(DP), intent(in) :: db
      real(SP), dimension(:), intent(in) :: FaA
      real(DP), dimension(nvar,mvar,na), intent(inout) :: DaB

      integer :: ieq,ild,ivar

      ! Compute B := cb*B
      call lalg_scaleVector(DaB,db)

      if (nvar .eq. mvar) then

        ! Compute B := B + A
        !$omp parallel do private(ild,ivar) default(shared) &
        !$omp if(NEQ > p_rperfconfig%NEQMIN_OMP)
        do ieq = 1,neq
          ild = Kdiagonal(ieq)
          do ivar = 1,nvar
            DaB(ivar,ivar,ild) = DaB(ivar,ivar,ild)+fa*FaA(ieq)
          end do
        end do
        !$omp end parallel do

      elseif (nvar .eq. 1) then

        ! Compute B := B + A
        !$omp parallel do private(ild) default(shared) &
        !$omp if(NEQ > p_rperfconfig%NEQMIN_OMP)
        do ieq = 1,neq
          ild = Kdiagonal(ieq)
          DaB(1,:,ild) = DaB(1,:,ild)+fa*FaA(ieq)
        end do
        !$omp end parallel do

      elseif (mvar .eq. 1) then

        ! Compute B := B + A
        !$omp parallel do private(ild) default(shared) &
        !$omp if(NEQ > p_rperfconfig%NEQMIN_OMP)
        do ieq = 1,neq
          ild = Kdiagonal(ieq)
          DaB(:,1,ild) = DaB(:,1,ild)+fa*FaA(ieq)
        end do
        !$omp end parallel do

      else   ! nvar /= mvar

        call output_line('Destination matrix must satisfy NVAR = MVAR!',&
            OU_CLASS_ERROR,OU_MODE_STD,'do_matDmat79SnglDble')
        call sys_halt()

      end if

    end subroutine

    !**************************************************************
    ! Format 7/9-7/9 addition:  C :=A + B
    ! This routinecomputes the number of nonzero matrix entries.
    !
    ! Remark: This subroutine is a modified version of the APLBDG
    !         subroutine taken from the SPARSEKIT library written
    !         by Youcef Saad.

    function do_mat79mat79add_computeNA(neq,KldA,KcolA,KldB,KcolB,Kaux) result(NA)

      integer, intent(in) :: neq
      integer, dimension(:), intent(in) :: KldA,KldB
      integer, dimension(:), intent(in) :: KcolB,KcolA
      integer, dimension(:), intent(inout) :: Kaux

      integer :: NA

      integer :: ieq,jeq,ild,icol,idg,ndg,last

      ! Initialization
      Kaux = 0; NA = 0

      do ieq = 1,neq

        ! For each row of matrix A
        ndg = 0

        ! End-of-linked list
        last = -1

        ! Row of matrix A
        do ild = KldA(ieq),KldA(ieq+1)-1

          ! Column number to be added
          icol = KcolA(ild)

          ! Add element to the linked list
          ndg = ndg+1
          Kaux(icol) = last
          last = icol
        end do

        ! Row of matrix B
        do ild = KldB(ieq),KldB(ieq+1)-1

          ! Column number to be added
          icol = KcolB(ild)

          ! Add element to the linked list
          if (Kaux(icol) .eq. 0) then
            ndg = ndg+1
            Kaux(icol) = last
            last = icol
          end if
        end do

        ! We are done with row IEQ
        NA = NA+ndg

        ! Reset KAUX to zero
        do idg = 1,ndg
          jeq = Kaux(last)
          Kaux(last) = 0
          last = jeq
        end do

      end do

    end function do_mat79mat79add_computeNA

    !**************************************************************
    ! Format 7/9-7/9 addition:  C := A + B
    ! Ths routine performs symbolical matrix-matrix-addition.
    !
    ! Remark: This subroutine is a modified version of the APLB1
    !         subroutine taken from the SPARSEKIT library written
    !         by Youcef Saad.

    subroutine do_mat79mat79add_symb(neq,ncols,KldA,KcolA,cmatrixFormatA,&
        KldB,KcolB,cmatrixFormatB,KldC,KcolC,Kdiagonal)

      integer, intent(in) :: neq,ncols
      integer, intent(in) :: cmatrixFormatA,cmatrixFormatB
      integer, dimension(:), intent(in) :: KldA,KldB
      integer, dimension(:), intent(in) :: KcolA,KcolB
      integer, dimension(:), intent(inout) :: KldC
      integer, dimension(:), intent(inout) :: KcolC
      integer, dimension(:), intent(inout), optional :: Kdiagonal

      integer :: ieq,ildA,ildB,ildC,ildendA,ildendB,icolA,icolB,icolC

      ! Initialization
      KldC(1) = 1; ildC = 1

      ! Loop over all rows
      do ieq = 1,neq

        ! Initialise column pointers for matrix A and B
        ildA = KldA(ieq); ildendA = KldA(ieq+1)-1
        ildB = KldB(ieq); ildendB = KldB(ieq+1)-1

        ! Check if diagonal entry needs to be stored at leading
        ! position for storage format CSR7. Then, both matrices A and
        ! B must be stored in storage format CSR7 so that the pointer
        ! ILDa and ILDB need to be increased by one (see below).
        if (.not.present(Kdiagonal)) then
          KcolC(ildC) = ieq
          ildC = ildC+1
        end if

        ! In any case, skip diagonal entry for matrix A and/or B if
        ! they are storage format 7
        if (cmatrixFormatA .eq. LSYSSC_MATRIX7) ildA = ildA+1
        if (cmatrixFormatB .eq. LSYSSC_MATRIX7) ildB = ildB+1

        ! For each row IEQ loop over the columns of matrix A and B
        ! simultaneously and collect the corresponding matrix entries
        do

          ! Find next column number for matrices A and B
          if (ildA .le. ildendA) then
            icolA = KcolA(ildA)
          else
            icolA = ncols+1
          end if

          if (ildB .le. ildendB) then
            icolB = KcolB(ildB)
          else
            icolB = ncols+1
          end if

          ! We need to consider three different cases. Since the
          ! diagonal element needs to be treated separately, the
          ! auxiliary variable ICOLC is used to store the provisional
          ! position of the next matrix entry
          if (icolA .eq. icolB) then
            ! Processing same column in matrix A and B
            icolC = icolA
            ildA = ildA+1
            ildB = ildB+1

          elseif (icolA < icolB) then
            ! Processing column in matrix A only
            icolC = icolA
            ildA = ildA+1

          else
            ! Processing column in matrix B only
            icolC = icolB
            ildB = ildB+1

          end if

          ! Now, we have the provisional position ICOLC of the next
          ! matrix entry. If it corresponds to the diagonal entry,
          ! then special care must be taken.
          if (icolC .eq. ieq) then

            ! If KDIAGONAL is present, then all entries in row IEQ
            ! are stored continuously but the diagonal entry is
            ! additionally marked in KDIAGONAL. IF KDIAGONAL is not
            ! present then the diagonal entry is already stored in
            ! the first position of each row (see above).
            if (present(Kdiagonal)) then
              Kdiagonal(ieq) = ildC
              KcolC(ildC) = icolC
              ildC = ildC+1
            end if

          else
            ! Off-diagonal entries are handled as usual
            KcolC(ildC) = icolC
            ildC = ildC+1
          end if

          ! Check if column IEQ is completed for both matrices A and B
          if (ildA > ildendA .and. ildB > ildendB) exit
        end do

        KldC(ieq+1) = ildC
      end do

    end subroutine

    !**************************************************************
    ! Format 7/9-7/9 addition:  B := ca*A + cb*B
    ! This routine perform numerical matrix-matrix-addition.

    subroutine do_mat79mat79addDbleDble(isizeIntl,neq,&
        naA,KldA,KcolA,naB,KldB,KcolB,KdiagonalB,DaA,DaB,da,db)

      integer, intent(in) :: isizeIntl,neq
      integer, intent(in) :: naA,naB
      integer, dimension(:), intent(in) :: KldA,KldB
      integer, dimension(:), intent(in) :: KcolA,KcolB
      integer, dimension(:), intent(in) :: KdiagonalB
      real(DP), intent(in) :: da,db
      real(DP), dimension(isizeIntl,naA), intent(in) :: DaA
      real(DP), dimension(isizeIntl,naB), intent(inout):: DaB

      integer :: ieq,ildA,ildB,icolA,icolB,ild

      ! Compute B := cb*B
      call lalg_scaleVector(DaB,db)

      ! Loop over all rows
      row: do ieq = 1, neq

        ! Initialise column pointer for matrix B
        ildB = KldB(ieq)

        ! Loop over all colums
        col: do ildA = KldA(ieq),KldA(ieq+1)-1

          ! Get column number
          icolA = KcolA(ildA)

          ! Are we on the diagonal?
          if (icolA .eq. ieq) then
            DaB(:,KdiagonalB(ieq)) = DaB(:,KdiagonalB(ieq)) + da*DaA(:,ildA)
          else
            ! Search for row/column number in matrix B
            do ild = ildB, KldB(ieq+1)-1
              icolB = KcolB(ild)
              if (icolB .eq. icolA) then
                ildB = ild
                DaB(:,ildB) = DaB(:,ildB) + da*DaA(:,ildA)
                cycle col
              end if
            end do

            do ild = KldB(ieq), ildB-1
              icolB = KcolB(ild)
              if (icolB .eq. icolA) then
                ildB = ild
                DaB(:,ildB) = DaB(:,ildB) + da*DaA(:,ildA)
                cycle col
              end if
            end do

            ! If we end up here, then the row/column number is not
            ! present in matrix B and we stop with an error
            call output_line('Destination matrix does not provide row/column number!',&
                OU_CLASS_ERROR,OU_MODE_STD,'do_mat79mat79addDbleDble')
            call sys_halt()
          end if

        end do col
      end do row

    end subroutine


    !**************************************************************
    ! Format 7/9-7/9 addition:  B := ca*A + cb*B
    ! This routine perform numerical matrix-matrix-addition.

    subroutine do_mat79mat79addSnglSngl(isizeIntl,neq,&
        naA,KldA,KcolA,naB,KldB,KcolB,KdiagonalB,FaA,FaB,fa,fb)

      integer, intent(in) :: isizeIntl,neq
      integer, intent(in) :: naA,naB
      integer, dimension(:), intent(in) :: KldA,KldB
      integer, dimension(:), intent(in) :: KcolA,KcolB
      integer, dimension(:), intent(in) :: KdiagonalB
      real(SP), intent(in) :: fa,fb
      real(SP), dimension(isizeIntl,naA), intent(in) :: FaA
      real(SP), dimension(isizeIntl,naB), intent(inout):: FaB

      integer :: ieq,ildA,ildB,icolA,icolB,ild

      ! Compute B := cb*B
      call lalg_scaleVector(FaB,fb)

      ! Loop over all rows
      row: do ieq = 1, neq

        ! Initialise column pointer for matrix B
        ildB = KldB(ieq)

        ! Loop over all colums
        col: do ildA = KldA(ieq),KldA(ieq+1)-1

          ! Get column number
          icolA = KcolA(ildA)

          ! Are we on the diagonal?
          if (icolA .eq. ieq) then
            FaB(:,KdiagonalB(ieq)) = FaB(:,KdiagonalB(ieq)) + fa*FaA(:,ildA)
          else
            ! Search for row/column number in matrix B
            do ild = ildB, KldB(ieq+1)-1
              icolB = KcolB(ild)
              if (icolB .eq. icolA) then
                ildB = ild
                FaB(:,ildB) = FaB(:,ildB) + fa*FaA(:,ildA)
                cycle col
              end if
            end do

            do ild = KldB(ieq), ildB-1
              icolB = KcolB(ild)
              if (icolB .eq. icolA) then
                ildB = ild
                FaB(:,ildB) = FaB(:,ildB) + fa*FaA(:,ildA)
                cycle col
              end if
            end do

            ! If we end up here, then the row/column number is not
            ! present in matrix B and we stop with an error
            call output_line('Destination matrix does not provide row/column number!',&
                OU_CLASS_ERROR,OU_MODE_STD,'do_mat79mat79addSnglSngl')
            call sys_halt()
          end if

        end do col
      end do row

    end subroutine

    !**************************************************************
    ! Format 7/9-7/9 addition:  B := ca*A + cb*B
    ! This routine perform numerical matrix-matrix-addition.

    subroutine do_mat79mat79addSnglDble(isizeIntl,neq,&
        naA,KldA,KcolA,naB,KldB,KcolB,KdiagonalB,FaA,DaB,fa,db)

      integer, intent(in) :: isizeIntl,neq
      integer, intent(in) :: naA,naB
      integer, dimension(:), intent(in) :: KldA,KldB
      integer, dimension(:), intent(in) :: KcolA,KcolB
      integer, dimension(:), intent(in) :: KdiagonalB
      real(SP), intent(in) :: fa
      real(DP), intent(in) :: db
      real(SP), dimension(isizeIntl,naA), intent(in) :: FaA
      real(DP), dimension(isizeIntl,naB), intent(inout):: DaB

      integer :: ieq,ildA,ildB,icolA,icolB,ild

      ! Compute B := cb*B
      call lalg_scaleVector(DaB,db)

      ! Loop over all rows
      row: do ieq = 1, neq

        ! Initialise column pointer for matrix B
        ildB = KldB(ieq)

        ! Loop over all colums
        col: do ildA = KldA(ieq),KldA(ieq+1)-1

          ! Get column number
          icolA = KcolA(ildA)

          ! Are we on the diagonal?
          if (icolA .eq. ieq) then
            DaB(:,KdiagonalB(ieq)) = DaB(:,KdiagonalB(ieq)) + fa*FaA(:,ildA)
          else
            ! Search for row/column number in matrix B
            do ild = ildB, KldB(ieq+1)-1
              icolB = KcolB(ild)
              if (icolB .eq. icolA) then
                ildB = ild
                DaB(:,ildB) = DaB(:,ildB) + fa*FaA(:,ildA)
                cycle col
              end if
            end do

            do ild = KldB(ieq), ildB-1
              icolB = KcolB(ild)
              if (icolB .eq. icolA) then
                ildB = ild
                DaB(:,ildB) = DaB(:,ildB) + fa*FaA(:,ildA)
                cycle col
              end if
            end do

            ! If we end up here, then the row/column number is not
            ! present in matrix B and we stop with an error
            call output_line('Destination matrix does not provide row/column number!',&
                OU_CLASS_ERROR,OU_MODE_STD,'do_mat79mat79addSnglDble')
            call sys_halt()
          end if

        end do col
      end do row

    end subroutine

    !**************************************************************
    ! Format 9rowc-9 addition:  B := ca*A + cb*B
    ! This routine perform numerical matrix-matrix-addition.

    subroutine do_mat9rowcMat9addDbleDble(nnzrows,&
        naA,KldA,KcolA,KrowIdxA,naB,KldB,KcolB,DaA,DaB,da,db)

      integer, intent(in) :: nnzrows
      integer, intent(in) :: naA,naB
      integer, dimension(:), intent(in) :: KldA,KldB
      integer, dimension(:), intent(in) :: KcolA,KcolB
      ! Last NNZROWS entries in KrowIdx
      integer, dimension(:), intent(in) :: KrowIdxA
      real(DP), intent(in) :: da,db
      real(DP), dimension(:), intent(in) :: DaA
      real(DP), dimension(:), intent(inout) :: DaB

      integer :: ieq,ild,ioffset

      ! Compute B := cb*B
      call lalg_scaleVector(DaB,db)

      ! Loop over all rows
      row: do ieq = 1, nnzrows

        ! Initialise column pointer for matrix B
        ioffset = KldB(KrowIdxA(ieq)) - KldA(ieq)

        ! Add the row
        do ild = KldA(ieq),KldA(ieq+1)-1
          DaB(ioffset+ild) = ca*DaA(ild)+cb*DaB(ioffset+ild)
        end do

      end do row

    end subroutine

!!$    !**************************************************************
!!$    ! Format 7/9-7/9 addition
!!$    ! Perform numerical matrix-matrix-addition C := ca * A + cb * B
!!$    !
!!$    ! Remark: This subroutine is a modified version of the APLB1
!!$    !         subroutine taken from the SPARSEKIT library written
!!$    !         by Youcef Saad.
!!$    !
!!$    ! The vectors KLDC and KCOLC are not necessary at first glance.
!!$    ! However, if the matrix C is allowed to have even a larger
!!$    ! sparsity pattern as the "sum" of A and B, then they are
!!$    ! required to find the correct positions in the final matrix C.
!!$    !
!!$    ! In this subroutine, KDIAGC is the vector which points to the
!!$    ! position of the diagonal entries. If matrix C is stored in
!!$    ! format CSR7 then KDIAGC corresponds to KLDC(1:NEQ). If matrix C
!!$    ! is stored in format CSR9 then KDIAGC corresponds to KDIAGONALC.
!!$
!!$    subroutine do_mat79mat79addDbleDble(isizeIntl,neq,ncols,&
!!$        KldA,KcolA,KldB,KcolB,DaA,DaB,ca,cb,KldC,KcolC,KdiagC,DaC)
!!$
!!$      integer, intent(in) :: neq,ncols
!!$      integer, intent(in) :: isizeIntl
!!$      integer, dimension(:), intent(in) :: KldA,KldB,KldC,KdiagC
!!$      integer, dimension(:), intent(in) :: KcolA,KcolB,KcolC
!!$      real(DP), intent(in) :: ca,cb
!!$      real(DP), dimension(isizeIntl,*), intent(in) :: DaA
!!$      real(DP), dimension(isizeIntl,*), intent(in) :: DaB
!!$      real(DP), dimension(isizeIntl,*), intent(inout):: DaC
!!$
!!$      integer :: ieq
!!$      integer :: ildA,ildB,ildC,ildendA,ildendB,ildendC
!!$      integer :: icolA,icolB,icolC,idiagC
!!$
!!$      ! Loop over all ROWS
!!$      do ieq=1,neq
!!$
!!$        ! Initialise column pointers for matrix A, B and C
!!$        ildA=KldA(ieq); ildendA=KldA(ieq+1)-1
!!$        ildB=KldB(ieq); ildendB=KldB(ieq+1)-1
!!$        ildC=KldC(ieq); ildendC=KldC(ieq+1)-1
!!$
!!$        ! Initialise pointer to diagonal entry
!!$        idiagC = KdiagC(ieq)
!!$
!!$        ! Since the final value of the diagonal entry
!!$        !    c_i,i = ca * a_i,i + cb * b_i,i
!!$        ! is updated step-by-step, set the diagonal entry to zero
!!$        DaC(:,idiagC) = 0._DP
!!$
!!$        ! For each row IEQ loop over the columns of matrix A and B
!!$        ! simultaneously and collect the corresponding matrix entries
!!$        do
!!$
!!$          ! Find next column number for matrices A and B
!!$          if (ildA .le. ildendA) then
!!$            icolA = KcolA(ildA)
!!$          else
!!$            icolA = ncols+1
!!$          end if
!!$
!!$          if (ildB .le. ildendB) then
!!$            icolB = KcolB(ildB)
!!$          else
!!$            icolB = ncols+1
!!$          end if
!!$
!!$          ! First, check if at least for one of the two matrices A
!!$          ! and/or B the diagonal entry which requires special
!!$          ! treatment has been reached. In this case, update (!!!)
!!$          ! the diagonal entry of the resulting matrix C immediately
!!$          ! and proceed to the next iteration
!!$          if (icolA .eq. ieq) then
!!$            if (icolB .eq. ieq) then
!!$
!!$              ! 1. Case: For both matrices A and B the diagonal entry
!!$              ! has been reached
!!$              DaC(:,idiagC)=ca*DaA(:,ildA)+cB*DaB(:,ildB)
!!$              ildA = ildA+1
!!$              ildB = ildB+1
!!$            else
!!$
!!$              ! 2. Case: For matrix A the diagonal entry has been
!!$              ! reached. Hence, skip matrix B.
!!$              DaC(:,idiagC)=DaC(:,idiagC)+ca*DaA(:,ildA)
!!$              ildA = ildA+1
!!$            end if
!!$          elseif (icolB .eq. ieq) then
!!$
!!$            !   3. Case: For matrix B the diagonal entry has been
!!$            !      reached. Hence, skip matrix A.
!!$            DaC(:,idiagC)=DaC(:,idiagC)+cB*DaB(:,ildB)
!!$            ildB = ildB+1
!!$
!!$          else
!!$
!!$            ! For both matrices A and B we have to process off-
!!$            ! -diagonal entries. Consider three different cases.
!!$            ! 1.) The next column is the same for both matrices A and B
!!$            ! 2.) The next column number is only present in matrix A
!!$            ! 3.) The next column number is only present in matrix B
!!$            !
!!$            ! Since matrix C is implicitly allowed to posses matrix
!!$            ! entries which are present neither in matrix A nor B, the
!!$            ! position ILDC needs to be updated of the current column
!!$            ! is not the diagonal entry whose position is known a
!!$            ! priori
!!$
!!$            if (icolA .eq. icolB) then
!!$
!!$              ! 1. Case: Processing same column in matrix A and B
!!$
!!$              ! Update column number for matrix C
!!$              do ildC=ildc,ildendC
!!$                icolC=KcolC(ildC)
!!$                if (icolC .eq. icolA) exit
!!$                if (icolC .ne. ieq) DaC(:,ildC) = 0._DP
!!$              end do
!!$
!!$              DaC(:,ildC)=ca*DaA(:,ildA)+cB*DaB(:,ildB)
!!$              ildA = ildA+1
!!$              ildB = ildB+1
!!$              ildC = ildC+1
!!$
!!$            elseif (icolA < icolB) then
!!$
!!$              ! 2. Case: Processing column in matrix A only
!!$
!!$              ! Update column number for matrix C
!!$              do ildC=ildC,ildendC
!!$                icolC=KcolC(ildC)
!!$                if (icolC .eq. icolA) exit
!!$                if (icolC .ne. ieq) DaC(:,ildC) = 0._DP
!!$              end do
!!$
!!$              DaC(:,ildC)=ca*DaA(:,ildA)
!!$              ildA = ildA+1
!!$              ildC = ildC+1
!!$
!!$            else
!!$
!!$              ! 3. Case: Processing column in matrix B only
!!$
!!$              ! Update column number for matrix C
!!$              do ildC=ildC,ildendC
!!$                icolC=KcolC(ildC)
!!$                if (icolC .eq. icolB) exit
!!$                if (icolC .ne. ieq) DaC(:,ildC) = 0._DP
!!$              end do
!!$
!!$              DaC(:,ildC)=cB*DaB(:,ildB)
!!$              ildB = ildB+1
!!$              ildC = ildC+1
!!$
!!$            end if
!!$          end if
!!$
!!$          ! Check if column IEQ is completed for both matrices A and B
!!$          if (ildA > ildendA .and. ildB > ildendB) exit
!!$        end do
!!$
!!$        ! Since matrix C is allowed to have additional column entries
!!$        ! which are not present in the "sum" of A and B, the
!!$        ! remainder of C needs to be nullified by hand
!!$        if (ildC .le. ildendC) then
!!$          if (KcolC(ildC) .eq. ieq) then
!!$            DaC(:,ildC+1:1:ildendC) = 0._DP
!!$          else
!!$            DaC(:,ildC:1:ildendC)   = 0._DP
!!$          end if
!!$        end if
!!$      end do
!!$    end subroutine

!!$    !**************************************************************
!!$    ! Format 7/9-7/9 addition
!!$    ! Perform numerical matrix-matrix-addition C := ca * A + cb * B
!!$    !
!!$    ! Remark: This subroutine is a modified version of the APLB1
!!$    !         subroutine taken from the SPARSEKIT library written
!!$    !         by Youcef Saad.
!!$    !
!!$    ! The vectors KLDC and KCOLC are not necessary at first glance.
!!$    ! However, if the matrix C is allowed to have even a larger
!!$    ! sparsity pattern as the "sum" of A and B, then they are
!!$    ! required to find the correct positions in the final matrix C.
!!$    !
!!$    ! In this subroutine, KDIAGC is the vector which points to the
!!$    ! position of the diagonal entries. If matrix C is stored in
!!$    ! format CSR7 then KDIAGC corresponds to KLDC(1:NEQ). If matrix C
!!$    ! is stored in format CSR9 then KDIAGC corresponds to KDIAGONALC.
!!$
!!$    subroutine do_mat79mat79addDbleSngl(isizeIntl,neq,ncols,&
!!$        KldA,KcolA,DaA,ca,KldB,KcolB,FaB,cb,KldC,KcolC,KdiagC,DaC)
!!$
!!$      integer, intent(in) :: neq,ncols
!!$      integer, intent(in) :: isizeIntl
!!$      integer, dimension(:), intent(in) :: KldA,KldB,KldC,KdiagC
!!$      integer, dimension(:), intent(in) :: KcolA,KcolB,KcolC
!!$      real(DP), intent(in) :: ca,cb
!!$      real(DP), dimension(isizeIntl,*), intent(in) :: DaA
!!$      real(SP), dimension(isizeIntl,*), intent(in) :: FaB
!!$      real(DP), dimension(isizeIntl,*), intent(inout):: DaC
!!$
!!$      integer :: ieq
!!$      integer :: ildA,ildB,ildC,ildendA,ildendB,ildendC
!!$      integer :: icolA,icolB,icolC,idiagC
!!$
!!$      ! Loop over all ROWS
!!$      do ieq=1,neq
!!$
!!$        ! Initialise column pointers for matrix A, B and C
!!$        ildA=KldA(ieq); ildendA=KldA(ieq+1)-1
!!$        ildB=KldB(ieq); ildendB=KldB(ieq+1)-1
!!$        ildC=KldC(ieq); ildendC=KldC(ieq+1)-1
!!$
!!$        ! Initialise pointer to diagonal entry
!!$        idiagC = KdiagC(ieq)
!!$
!!$        ! Since the final value of the diagonal entry
!!$        !    c_i,i = ca * a_i,i + cb * b_i,i
!!$        ! is updated step-by-step, set the diagonal entry to zero
!!$        DaC(:,idiagC) = 0._DP
!!$
!!$        ! For each row IEQ loop over the columns of matrix A and B
!!$        ! simultaneously and collect the corresponding matrix entries
!!$        do
!!$
!!$          ! Find next column number for matrices A and B
!!$          if (ildA .le. ildendA) then
!!$            icolA = KcolA(ildA)
!!$          else
!!$            icolA = ncols+1
!!$          end if
!!$
!!$          if (ildB .le. ildendB) then
!!$            icolB = KcolB(ildB)
!!$          else
!!$            icolB = ncols+1
!!$          end if
!!$
!!$          ! First, check if at least for one of the two matrices A
!!$          ! and/or B the diagonal entry which requires special
!!$          ! treatment has been reached. In this case, update (!!!)
!!$          ! the diagonal entry of the resulting matrix C immediately
!!$          ! and proceed to the next iteration
!!$          if (icolA .eq. ieq) then
!!$            if (icolB .eq. ieq) then
!!$
!!$              ! 1. Case: For both matrices A and B the diagonal entry
!!$              ! has been reached
!!$              DaC(:,idiagC)=ca*DaA(:,ildA)+cB*FaB(:,ildB)
!!$              ildA = ildA+1
!!$              ildB = ildB+1
!!$            else
!!$
!!$              ! 2. Case: For matrix A the diagonal entry has been
!!$              ! reached. Hence, skip matrix B.
!!$              DaC(:,idiagC)=DaC(:,idiagC)+ca*DaA(:,ildA)
!!$              ildA = ildA+1
!!$            end if
!!$          elseif (icolB .eq. ieq) then
!!$
!!$            !   3. Case: For matrix B the diagonal entry has been
!!$            !      reached. Hence, skip matrix A.
!!$            DaC(:,idiagC)=DaC(:,idiagC)+cB*FaB(:,ildB)
!!$            ildB = ildB+1
!!$
!!$          else
!!$
!!$            ! For both matrices A and B we have to process off-
!!$            ! -diagonal entries. Consider three different cases.
!!$            ! 1.) The next column is the same for both matrices A and B
!!$            ! 2.) The next column number is only present in matrix A
!!$            ! 3.) The next column number is only present in matrix B
!!$            !
!!$            ! Since matrix C is implicitly allowed to posses matrix
!!$            ! entries which are present neither in matrix A nor B, the
!!$            ! position ILDC needs to be updated of the current column
!!$            ! is not the diagonal entry whose position is known a
!!$            ! priori
!!$
!!$            if (icolA .eq. icolB) then
!!$
!!$              ! 1. Case: Processing same column in matrix A and B
!!$
!!$              ! Update column number for matrix C
!!$              do ildC=ildc,ildendC
!!$                icolC=KcolC(ildC)
!!$                if (icolC .eq. icolA) exit
!!$                if (icolC .ne. ieq) DaC(:,ildC) = 0._DP
!!$              end do
!!$
!!$              DaC(:,ildC)=ca*DaA(:,ildA)+cB*FaB(:,ildB)
!!$              ildA = ildA+1
!!$              ildB = ildB+1
!!$              ildC = ildC+1
!!$
!!$            elseif (icolA < icolB) then
!!$
!!$              ! 2. Case: Processing column in matrix A only
!!$
!!$              ! Update column number for matrix C
!!$              do ildC=ildC,ildendC
!!$                icolC=KcolC(ildC)
!!$                if (icolC .eq. icolA) exit
!!$                if (icolC .ne. ieq) DaC(:,ildC) = 0._DP
!!$              end do
!!$
!!$              DaC(:,ildC)=ca*DaA(:,ildA)
!!$              ildA = ildA+1
!!$              ildC = ildC+1
!!$
!!$            else
!!$
!!$              ! 3. Case: Processing column in matrix B only
!!$
!!$              ! Update column number for matrix C
!!$              do ildC=ildC,ildendC
!!$                icolC=KcolC(ildC)
!!$                if (icolC .eq. icolB) exit
!!$                if (icolC .ne. ieq) DaC(:,ildC) = 0._DP
!!$              end do
!!$
!!$              DaC(:,ildC)=cB*FaB(:,ildB)
!!$              ildB = ildB+1
!!$              ildC = ildC+1
!!$
!!$            end if
!!$          end if
!!$
!!$          ! Check if column IEQ is completed for both matrices A and B
!!$          if (ildA > ildendA .and. ildB > ildendB) exit
!!$        end do
!!$
!!$        ! Since matrix C is allowed to have additional column entries
!!$        ! which are not present in the "sum" of A and B, the
!!$        ! remainder of C needs to be nullified by hand
!!$        if (ildC .le. ildendC) then
!!$          if (KcolC(ildC) .eq. ieq) then
!!$            DaC(:,ildC+1:1:ildendC) = 0._DP
!!$          else
!!$            DaC(:,ildC:1:ildendC)   = 0._DP
!!$          end if
!!$        end if
!!$      end do
!!$    end subroutine

!!$    !**************************************************************
!!$    ! Format 7/9-7/9 addition
!!$    ! Perform numerical matrix-matrix-addition C := ca * A + cb * B
!!$    !
!!$    ! Remark: This subroutine is a modified version of the APLB1
!!$    !         subroutine taken from the SPARSEKIT library written
!!$    !         by Youcef Saad.
!!$    !
!!$    ! The vectors KLDC and KCOLC are not necessary at first glance.
!!$    ! However, if the matrix C is allowed to have even a larger
!!$    ! sparsity pattern as the "sum" of A and B, then they are
!!$    ! required to find the correct positions in the final matrix C.
!!$    !
!!$    ! In this subroutine, KDIAGC is the vector which points to the
!!$    ! position of the diagonal entries. If matrix C is stored in
!!$    ! format CSR7 then KDIAGC corresponds to KLDC(1:NEQ). If matrix C
!!$    ! is stored in format CSR9 then KDIAGC corresponds to KDIAGONALC.
!!$
!!$    subroutine do_mat79mat79addSnglSngl(isizeIntl,neq,ncols,&
!!$        KldA,KcolA,FaA,ca,KldB,KcolB,FaB,cb,KldC,KcolC,KdiagC,FaC)
!!$
!!$      integer, intent(in) :: neq,ncols
!!$      integer, intent(in) :: isizeIntl
!!$      integer, dimension(:), intent(in) :: KldA,KldB,KldC,KdiagC
!!$      integer, dimension(:), intent(in) :: KcolA,KcolB,KcolC
!!$      real(DP), intent(in) :: ca,cb
!!$      real(SP), dimension(isizeIntl,*), intent(in) :: FaA
!!$      real(SP), dimension(isizeIntl,*), intent(in) :: FaB
!!$      real(SP), dimension(isizeIntl,*), intent(inout):: FaC
!!$
!!$      integer :: ieq
!!$      integer :: ildA,ildB,ildC,ildendA,ildendB,ildendC
!!$      integer :: icolA,icolB,icolC,idiagC
!!$
!!$      ! Loop over all ROWS
!!$      do ieq=1,neq
!!$
!!$        ! Initialise column pointers for matrix A, B and C
!!$        ildA=KldA(ieq); ildendA=KldA(ieq+1)-1
!!$        ildB=KldB(ieq); ildendB=KldB(ieq+1)-1
!!$        ildC=KldC(ieq); ildendC=KldC(ieq+1)-1
!!$
!!$        ! Initialise pointer to diagonal entry
!!$        idiagC = KdiagC(ieq)
!!$
!!$        ! Since the final value of the diagonal entry
!!$        !    c_i,i = ca * a_i,i + cb * b_i,i
!!$        ! is updated step-by-step, set the diagonal entry to zero
!!$        FaC(:,idiagC) = 0._SP
!!$
!!$        ! For each row IEQ loop over the columns of matrix A and B
!!$        ! simultaneously and collect the corresponding matrix entries
!!$        do
!!$
!!$          ! Find next column number for matrices A and B
!!$          if (ildA .le. ildendA) then
!!$            icolA = KcolA(ildA)
!!$          else
!!$            icolA = ncols+1
!!$          end if
!!$
!!$          if (ildB .le. ildendB) then
!!$            icolB = KcolB(ildB)
!!$          else
!!$            icolB = ncols+1
!!$          end if
!!$
!!$          ! First, check if at least for one of the two matrices A
!!$          ! and/or B the diagonal entry which requires special
!!$          ! treatment has been reached. In this case, update (!!!)
!!$          ! the diagonal entry of the resulting matrix C immediately
!!$          ! and proceed to the next iteration
!!$          if (icolA .eq. ieq) then
!!$            if (icolB .eq. ieq) then
!!$
!!$              ! 1. Case: For both matrices A and B the diagonal entry
!!$              ! has been reached
!!$              FaC(:,idiagC)=ca*FaA(:,ildA)+cB*FaB(:,ildB)
!!$              ildA = ildA+1
!!$              ildB = ildB+1
!!$            else
!!$
!!$              ! 2. Case: For matrix A the diagonal entry has been
!!$              ! reached. Hence, skip matrix B.
!!$              FaC(:,idiagC)=FaC(:,idiagC)+ca*FaA(:,ildA)
!!$              ildA = ildA+1
!!$            end if
!!$          elseif (icolB .eq. ieq) then
!!$
!!$            !   3. Case: For matrix B the diagonal entry has been
!!$            !      reached. Hence, skip matrix A.
!!$            FaC(:,idiagC)=FaC(:,idiagC)+cB*FaB(:,ildB)
!!$            ildB = ildB+1
!!$
!!$          else
!!$
!!$            ! For both matrices A and B we have to process off-
!!$            ! -diagonal entries. Consider three different cases.
!!$            ! 1.) The next column is the same for both matrices A and B
!!$            ! 2.) The next column number is only present in matrix A
!!$            ! 3.) The next column number is only present in matrix B
!!$            !
!!$            ! Since matrix C is implicitly allowed to posses matrix
!!$            ! entries which are present neither in matrix A nor B, the
!!$            ! position ILDC needs to be updated of the current column
!!$            ! is not the diagonal entry whose position is known a
!!$            ! priori
!!$
!!$            if (icolA .eq. icolB) then
!!$
!!$              ! 1. Case: Processing same column in matrix A and B
!!$
!!$              ! Update column number for matrix C
!!$              do ildC=ildc,ildendC
!!$                icolC=KcolC(ildC)
!!$                if (icolC .eq. icolA) exit
!!$                if (icolC .ne. ieq) FaC(:,ildC) = 0._SP
!!$              end do
!!$
!!$              FaC(:,ildC)=ca*FaA(:,ildA)+cB*FaB(:,ildB)
!!$              ildA = ildA+1
!!$              ildB = ildB+1
!!$              ildC = ildC+1
!!$
!!$            elseif (icolA < icolB) then
!!$
!!$              ! 2. Case: Processing column in matrix A only
!!$
!!$              ! Update column number for matrix C
!!$              do ildC=ildC,ildendC
!!$                icolC=KcolC(ildC)
!!$                if (icolC .eq. icolA) exit
!!$                if (icolC .ne. ieq) FaC(:,ildC) = 0._SP
!!$              end do
!!$
!!$              FaC(:,ildC)=ca*FaA(:,ildA)
!!$              ildA = ildA+1
!!$              ildC = ildC+1
!!$
!!$            else
!!$
!!$              ! 3. Case: Processing column in matrix B only
!!$
!!$              ! Update column number for matrix C
!!$              do ildC=ildC,ildendC
!!$                icolC=KcolC(ildC)
!!$                if (icolC .eq. icolB) exit
!!$                if (icolC .ne. ieq) FaC(:,ildC) = 0._SP
!!$              end do
!!$
!!$              FaC(:,ildC)=cB*FaB(:,ildB)
!!$              ildB = ildB+1
!!$              ildC = ildC+1
!!$
!!$            end if
!!$          end if
!!$
!!$          ! Check if column IEQ is completed for both matrices A and B
!!$          if (ildA > ildendA .and. ildB > ildendB) exit
!!$        end do
!!$
!!$        ! Since matrix C is allowed to have additional column entries
!!$        ! which are not present in the "sum" of A and B, the
!!$        ! remainder of C needs to be nullified by hand
!!$        if (ildC .le. ildendC) then
!!$          if (KcolC(ildC) .eq. ieq) then
!!$            FaC(:,ildC+1:1:ildendC) = 0._SP
!!$          else
!!$            FaC(:,ildC:1:ildendC)   = 0._SP
!!$          end if
!!$        end if
!!$      end do
!!$    end subroutine

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_swapVectors(rvector1,rvector2)

!<description>
    ! This subroutine swaps the content of two different vectors
!</description>

!<inputoutput>
    ! first scalar vector
    type(t_vectorScalar), intent(inout) :: rvector1

    ! second scalar vector
    type(t_vectorScalar), intent(inout) :: rvector2
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorScalar) :: rvector

    ! Vector1 -> Vector
    rvector  = rvector1

    ! Vector2 -> Vector1
    rvector1 = rvector2

    ! Vector -> Vector2
    rvector2 = rvector

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_swapMatrices(rmatrix1,rmatrix2)

!<description>
    ! This subroutine swaps the content of two different matrices
!</description>

!<inputoutput>
    ! first scalar matrix
    type(t_matrixScalar), intent(inout) :: rmatrix1

    ! second scalar matrix
    type(t_matrixScalar), intent(inout) :: rmatrix2
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_matrixScalar) :: rmatrix

    ! Matrix1 -> Matrix
    rmatrix = rmatrix1

    ! Matrix2 -> Matrix1
    rmatrix1 = rmatrix2

    ! Matrix -> Matrix2
    rmatrix2 = rmatrix

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_assignDiscrDirectMat (rmatrix,rdiscrTrial,rdiscrTest)

!<description>
  ! Assigns given discretisation structures for trial/test spaces to a
  ! matrix.
!</description>

!<input>
  ! Discretisation structure for trial functions.
  type(t_spatialDiscretisation), intent(in), target :: rdiscrTrial

  ! OPTIONAL: Discretisation structure for test functions.
  ! If not specified, trial and test functions coincide.
  type(t_spatialDiscretisation), intent(in), target, optional :: rdiscrTest
!</input>

!<inputoutput>
  ! Destination matrix.
  type(t_matrixScalar),intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    if (rmatrix%NCOLS .ne. dof_igetNDofGlob(rdiscrTrial)) then
      call output_line ('Discretisation invalid for the matrix!', &
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_assignDiscrDirectMat')
      call sys_halt()
    end if

    rmatrix%p_rspatialDiscrTrial => rdiscrTrial

    ! Depending on whether rdiscrTest is given, set the pointers for
    ! the test functions.
    rmatrix%bidenticalTrialAndTest = present(rdiscrTest)

    if (present(rdiscrTest)) then

      if (rmatrix%NEQ .ne. dof_igetNDofGlob(rdiscrTest)) then
        call output_line ('Discretisation invalid for the matrix!', &
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_assignDiscrDirectMat')
        call sys_halt()
      end if

      ! Set the block discretisation of the block matrix
      rmatrix%p_rspatialDiscrTest => rdiscrTest
    else
      ! Trial and test functions coincide
      if (rmatrix%NEQ .ne. dof_igetNDofGlob(rdiscrTrial)) then
        call output_line ('Discretisation invalid for the matrix!', &
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_assignDiscrDirectMat')
        call sys_halt()
      end if

      ! Set the block discretisation of the block matrix
      rmatrix%p_rspatialDiscrTest => rdiscrTrial
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_infoMatrix(rmatrix)

!<description>
    ! This subroutine outputs information about the matrix
!</description>

!<input>
    ! scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>
!</subroutine>

    ! local variables
    integer :: isize,isortcols,isortrows
    
    isortcols = 0
    isortrows = 0
    if (rmatrix%bcolumnsSorted) isortcols = 1
    if (rmatrix%browsSorted) isortrows = 1

    call output_line ('ScalarMatrix:')
    call output_line ('-------------')
    call output_line ('UUID:                    '//uuid_conv2String(rmatrix%ruuid))
    call output_line ('cmatrixFormat:           '//trim(sys_siL(rmatrix%cmatrixFormat,15)))
    call output_line ('cinterleavematrixFormat: '//trim(sys_siL(rmatrix%cinterleaveMatrixFormat,15)))
    call output_line ('cdataType:               '//trim(sys_siL(rmatrix%cdataType,15)))
    call output_line ('imatrixSpec:             '//trim(sys_siL(int(rmatrix%imatrixSpec),15)))
    call output_line ('NA:                      '//trim(sys_siL(rmatrix%NA,15)))
    call output_line ('NEQ:                     '//trim(sys_siL(rmatrix%NEQ,15)))
    call output_line ('NCOLS:                   '//trim(sys_siL(rmatrix%NCOLS,15)))
    call output_line ('NVAR:                    '//trim(sys_siL(rmatrix%NVAR,15)))
    call output_line ('dscaleFactor:            '//trim(sys_sdL(rmatrix%dscaleFactor,2)))
    call output_line ('bcolumnsSorted:          '//trim(sys_siL(isortcols,2)))
    if (associated(rmatrix%p_rsortStrategyColumns)) then
      call output_line ('isortStrategy cols:      '//trim(sys_siL(rmatrix%p_rsortStrategyColumns%ctype,15)))
      call output_line ('h_Ipermutation cols:     '//trim(sys_siL(rmatrix%p_rsortStrategyColumns%h_Ipermutation,15)))
    else
      call output_line ('isortStrategy:           '//trim(sys_siL(0,15)))
      call output_line ('h_Ipermutation:          '//trim(sys_siL(0,15)))
    end if
    call output_line ('browsSorted:             '//trim(sys_siL(isortrows,2)))
    if (associated(rmatrix%p_rsortStrategyRows)) then
      call output_line ('isortStrategy rows:      '//trim(sys_siL(rmatrix%p_rsortStrategyRows%ctype,15)))
      call output_line ('h_Ipermutation rows:     '//trim(sys_siL(rmatrix%p_rsortStrategyRows%h_Ipermutation,15)))
    else
      call output_line ('isortStrategy:           '//trim(sys_siL(0,15)))
      call output_line ('h_Ipermutation:          '//trim(sys_siL(0,15)))
    end if
    call output_line ('h_Da:                    '//trim(sys_siL(rmatrix%h_Da,15)))
    if (rmatrix%h_Da .ne. ST_NOHANDLE) then
      call storage_getsize(rmatrix%h_Da,isize)
      select case (rmatrix%cinterleaveMatrixFormat)
        case (LSYSSC_MATRIX1)
          call output_line ('Da memory usage:         '//&
              trim(sys_sdL(100/real(isize,DP)*rmatrix%NA*rmatrix%NVAR*rmatrix%NVAR,2))//'%')
        case (LSYSSC_MATRIXD)
          call output_line ('Da memory usage:         '//&
              trim(sys_sdL(100/real(isize,DP)*rmatrix%NA*rmatrix%NVAR,2))//'%')
        case default
          call output_line ('Da memory usage:         '//trim(sys_sdL(100/real(isize,DP)*rmatrix%NA,2))//'%')
        end select
    end if
    call output_line ('h_Kcol:                  '//trim(sys_siL(rmatrix%h_Kcol,15)))
    if (rmatrix%h_Kcol .ne. ST_NOHANDLE) then
      call storage_getsize(rmatrix%h_Kcol,isize)
      call output_line ('Kcol memory usage:       '//trim(sys_sdL(100/real(isize,DP)*rmatrix%NA,2))//'%')
    end if
    call output_line ('h_Kld:                   '//trim(sys_siL(rmatrix%h_Kld,15)))
    if (rmatrix%h_Kld .ne. ST_NOHANDLE) then
      call storage_getsize(rmatrix%h_Kld,isize)
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .eq. 0) then
        call output_line ('Kld memory usage:        '//trim(sys_sdL(100/real(isize,DP)*(rmatrix%NEQ+1),2))//'%')
      else
        call output_line ('Kld memory usage:        '//trim(sys_sdL(100/real(isize,DP)*(rmatrix%NCOLS+1),2))//'%')
      end if
    end if
    call output_line ('h_Kdiagonal:             '//trim(sys_siL(rmatrix%h_Kdiagonal,15)))
    if (rmatrix%h_Kdiagonal .ne. ST_NOHANDLE) then
      call storage_getsize(rmatrix%h_Kdiagonal,isize)
      if (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .eq. 0) then
        call output_line ('Kdiagonal memory usage:  '//trim(sys_sdL(100/real(isize,DP)*(rmatrix%NEQ),2))//'%')
      else
        call output_line ('Kdiagonl memory usage:   '//trim(sys_sdL(100/real(isize,DP)*(rmatrix%NCOLS),2))//'%')
      end if
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_infoVector(rvector)

!<description>
    ! This subroutine outputs information about the vector
!</description>

!<input>
    ! scalar vector
    type(t_vectorScalar), intent(in) :: rvector
!</input>
!</subroutine>

    ! local variables
    integer :: isize,isort
    
    isort = 0
    if (rvector%bisSorted) isort = 1

    call output_line ('ScalarVector:')
    call output_line ('-------------')
    call output_line ('UUID:                   '//uuid_conv2String(rvector%ruuid))
    call output_line ('cdataType:              '//trim(sys_siL(rvector%cdataType,15)))
    call output_line ('NEQ:                    '//trim(sys_siL(rvector%NEQ,15)))
    call output_line ('NVAR:                   '//trim(sys_siL(rvector%NVAR,15)))
    call output_line ('bisSorted:              '//trim(sys_siL(isort,2)))
    if (associated(rvector%p_rsortStrategy)) then
      call output_line ('isortStrategy:          '//trim(sys_siL(rvector%p_rsortStrategy%ctype,15)))
      call output_line ('h_Ipermutation:         '//trim(sys_siL(rvector%p_rsortStrategy%h_Ipermutation,15)))
    else
      call output_line ('isortStrategy:          '//trim(sys_siL(0,15)))
      call output_line ('h_Ipermutation:         '//trim(sys_siL(0,15)))
    end if
    call output_line ('h_Ddata:                '//trim(sys_siL(rvector%h_Ddata,15)))
    if (rvector%h_Ddata .ne. ST_NOHANDLE) then
      call storage_getsize(rvector%h_Ddata,isize)
      call output_line ('Ddata memory usage:     '//trim(sys_sdL(100/real(isize,DP)*rvector%NEQ*rvector%NVAR,2))//'%')
    end if
    call output_line ('iidxFirstEntry:         '//trim(sys_siL(rvector%iidxFirstEntry,15)))
    write(*,FMT='(A)')       '-------------------------'

  end subroutine

  ! ***************************************************************************

!<function>

  elemental logical function lsyssc_isMatrixStructureShared (rmatrix,&
      rmatrix2) result (bresult)

!<description>
  ! Checks whether the structure of rmatrix belongs to that matrix or if it is
  ! shared among rmatrix and another matrix.
  !
  ! Note: If both matrices rmatrix and rmatrix2 are specified and both matrices
  !   do not contain structure data, the return value is of course FALSE.
!</description>

!<input>
  ! A scalar matrix to be checked.
  type(t_matrixScalar), intent(in) :: rmatrix

  ! OPTIONAL: A second matrix to be compared with rmatrix.
  ! If not specified, lsyssc_isMatrixStructureShared checks if rmatrix
  !   shares its structure with any other matrix.
  ! If specified, lsyssc_isMatrixStructureShared checks if rmatrix
  !   shares its entries with rmatrix2.
  type(t_matrixScalar), intent(in), optional :: rmatrix2
!</input>

!<result>
  ! TRUE, if the structure arrays in rmatrix are shared among rmatrix and
  !   another matrix or with the matrix rmatrix2, respectively.
  ! FALSE, if the structure arrays in rmatrix belong to rmatrix.
!</result>

!</function>

    if (.not. present(rmatrix2)) then
      ! General check
      bresult = iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0
    else
      ! Check if the structure is shared among rmatrix and rmatrix2.
      ! This depends on the matrix format...
      ! The matrix is declared as 'structure is shared' if at least one
      ! of the handles that define the structure is shared among the matrices.
      bresult = .false.
      if (rmatrix%cmatrixFormat .ne. rmatrix2%cmatrixFormat) return
      select case (rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX1,LSYSSC_MATRIXD)
        ! No structure, full matrix
        bresult = .false.
      case (LSYSSC_MATRIX7,LSYSSC_MATRIX7INTL)
        bresult = ((rmatrix%h_Kcol .ne. ST_NOHANDLE) .and. &
                   (rmatrix%h_Kcol .eq. rmatrix2%h_Kcol)) .or. &
                  ((rmatrix%h_Kld .ne. ST_NOHANDLE) .and. &
                   (rmatrix%h_Kld .eq. rmatrix2%h_Kld))
      case (LSYSSC_MATRIX9,LSYSSC_MATRIX9INTL)
        bresult = ((rmatrix%h_Kcol .ne. ST_NOHANDLE) .and. &
                   (rmatrix%h_Kcol .eq. rmatrix2%h_Kcol)) .or. &
                  ((rmatrix%h_Kld .ne. ST_NOHANDLE) .and. &
                   (rmatrix%h_Kld .eq. rmatrix2%h_Kld)) .or. &
                  ((rmatrix%h_Kdiagonal .ne. ST_NOHANDLE) .and. &
                   (rmatrix%h_Kdiagonal .eq. rmatrix2%h_Kdiagonal))
      end select
    end if

  end function

  ! ***************************************************************************

!<function>

  elemental logical function lsyssc_isMatrixContentShared (rmatrix, &
      rmatrix2) result (bresult)

!<description>
  ! Checks whether the content of rmatrix belongs to that matrix or if it is
  ! shared among rmatrix and another matrix.
  !
  ! Note: If both matrices rmatrix and rmatrix2 are specified and both matrices
  !   don''t contain content data, the return value is of course FALSE.
!</description>

!<input>
  ! A scalar matrix
  type(t_matrixScalar), intent(in) :: rmatrix

  ! OPTIONAL: A second matrix to be compared with rmatrix.
  ! If not specified, lsyssc_isMatrixContentShared checks if rmatrix
  !   shares its content with any other matrix.
  ! If specified, lsyssc_isMatrixStructureShared checks if rmatrix
  !   shares its content with rmatrix2.
  type(t_matrixScalar), intent(in), optional :: rmatrix2
!</input>

!<result>
  ! TRUE, if the content arrays in rmatrix are shared among rmatrix and
  !   another matrix or with the matrix rmatrix2, respectively.
  ! FALSE, if the content arrays in rmatrix belong to rmatrix or if
  !   rmatrix does not have content data.
!</result>

!</function>

    if (.not. present(rmatrix2)) then
      bresult = iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .ne. 0
    else
      ! Check if the content is shared among rmatrix and rmatrix2.
      ! This depends on the matrix format...
      ! The matrix is declared as 'content is shared' if at least one
      ! of the handles that define the structure is shared among the matrices.
      bresult = .false.
      if (rmatrix%cmatrixFormat .ne. rmatrix2%cmatrixFormat) return
      select case (rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX1,LSYSSC_MATRIX7,LSYSSC_MATRIX9,LSYSSC_MATRIXD, &
            LSYSSC_MATRIX7INTL,LSYSSC_MATRIX9INTL)
        bresult = ((rmatrix%h_Da .ne. ST_NOHANDLE) .and. &
                   (rmatrix%h_Da .eq. rmatrix2%h_Da))
      end select
    end if

  end function

  !****************************************************************************

!<function>

  elemental logical function lsyssc_isMatrixPresent (rmatrix, &
      bignoreScaleFactor) result (bispresent)

!<description>
  ! This routine checks if the matrix is defined or empty.
!</description>

!<input>
  ! The matrix to be checked
  type(t_matrixScalar), intent(in) :: rmatrix

  ! OPTIONAL: Whether to check the scaling factor.
  ! FALSE: A scaling factor of 0.0 disables a submatrix.
  !        This is the standard setting.
  ! TRUE: The scaling factor is ignored.
  logical, intent(in), optional :: bignoreScaleFactor
!</input>

!<output>
  ! Whether the matrix structure realises an existing matrix exists or not.
!</output>

!</function>
    logical :: bscale

    if (present(bignoreScaleFactor)) then
      bscale = bignoreScaleFactor
    else
      bscale = .false.
    end if

    bispresent = &
      (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIXUNDEFINED) &
      .and. ((.not. bscale) .or. &
             (rmatrix%dscaleFactor .ne. 0.0_DP))

  end function

  !****************************************************************************

!<function>

  pure logical function lsyssc_hasMatrixStructure (rmatrix)

!<description>
  ! Returns whether a matrix has a structure or not.
  !
  ! Note that some matrix types (e.g. matrix-type 1 = full matrix)
  ! do not have a structure at all, so the routine always returns
  ! FALSE in such a case.
!</description>

!<input>
  ! Matrix to check
  type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<result>
  ! Whether the matrix has strucure arrays in memory or not.
!</result>

!</function>

    ! All up to now implemented matrix types use Kcol if they have a
    ! structure...
    lsyssc_hasMatrixStructure = rmatrix%h_Kcol .ne. ST_NOHANDLE

  end function

  !****************************************************************************

!<function>

  pure logical function lsyssc_hasMatrixContent (rmatrix)

!<description>
  ! Returns whether a matrix has a content or not.
!</description>

!<input>
  ! Matrix to check
  type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<result>
  ! Whether the matrix has a content array in memory or not.
!</result>

!</function>

    ! All up to now implemented matrix types save their data in h_Da.
    lsyssc_hasMatrixContent = rmatrix%h_Da .ne. ST_NOHANDLE

  end function

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_spreadVector(rvector1, rvector2)

!<description>
    ! This subroutine spreads a scalar vector which is not stored in interleave
    ! format into another vector which is stored in interleave format.
!</description>

!<input>
    ! Scalar source vector
    type(t_vectorScalar), intent(in) :: rvector1
!</input>

!<inputoutput>
    ! Scalar destination vector in interleave format
    type(t_vectorScalar), intent(inout) :: rvector2
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2
    real(SP), dimension(:), pointer :: p_Fdata1,p_Fdata2
    integer, dimension(:), pointer :: p_Idata1,p_Idata2

    ! Source vector must not be stored in interleave format
    if (rvector1%NVAR .ne. 1) then
      call output_line("Source vector must not be stored in interleave format!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_spreadVector")
      call sys_halt()
    end if

    ! Vectors must have the same size
    if (rvector1%NEQ .ne. rvector2%NEQ) then
      call output_line("Vectors not compatible, different size!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_spreadVector")
      call sys_halt()
    end if

    ! Vectors must have the data type
    if (rvector1%cdataType .ne. rvector2%cdataType) then
      call output_line("Vectors not compatible, different data type!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_spreadVector")
      call sys_halt()
    end if

    ! Check the sorting
    if (rvector1%bisSorted .neqv. rvector2%bisSorted) then
      call output_line("Vectors not compatible, differently sorted!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_spreadVector")
      call sys_halt()
    end if
    
    if (rvector1%bisSorted) then
      if (rvector1%p_rsortStrategy%ctype .ne. rvector2%p_rsortStrategy%ctype) then
        call output_line("Vectors not compatible, differently sorted!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_spreadVector")
        call sys_halt()
      end if
    end if

    select case (rvector1%cdataType)
    case (ST_DOUBLE)
      call lsyssc_getbase_double(rvector1, p_Ddata1)
      call lsyssc_getbase_double(rvector2, p_Ddata2)
      call do_spreadDble(p_Ddata1, rvector2%NVAR, rvector2%NEQ, p_Ddata2)

    case (ST_SINGLE)
      call lsyssc_getbase_single(rvector1, p_Fdata1)
      call lsyssc_getbase_single(rvector2, p_Fdata2)
      call do_spreadSngl(p_Fdata1, rvector2%NVAR, rvector2%NEQ, p_Fdata2)

    case (ST_INT)
      call lsyssc_getbase_int(rvector1, p_Idata1)
      call lsyssc_getbase_int(rvector2, p_Idata2)
      call do_spreadInt(p_Idata1, rvector2%NVAR, rvector2%NEQ, p_Idata2)

    case default
      call output_line("Unsupported data type!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_spreadVector")
      call sys_halt()
    end select

  contains

    ! Here, the real working routines follow.

    !**************************************************************

    subroutine do_spreadDble(Ddata1, NVAR, NEQ, Ddata2)
      real(DP), dimension(:), intent(in) :: Ddata1
      integer, intent(in) :: NVAR
      integer, intent(in) :: NEQ
      real(DP), dimension(NVAR,NEQ), intent(out) :: Ddata2

      integer :: ieq

      do ieq = 1, NEQ
        Ddata2(:,ieq) = Ddata1(ieq)
      end do

    end subroutine

    !**************************************************************

    subroutine do_spreadSngl(Fdata1, NVAR, NEQ, Fdata2)
      real(SP), dimension(:), intent(in) :: Fdata1
      integer, intent(in) :: NVAR
      integer, intent(in) :: NEQ
      real(SP), dimension(NVAR,NEQ), intent(out) :: Fdata2

      integer :: ieq

      do ieq = 1, NEQ
        Fdata2(:,ieq) = Fdata1(ieq)
      end do

    end subroutine

    !**************************************************************

    subroutine do_spreadInt(Idata1, NVAR, NEQ, Idata2)
      integer, dimension(:), intent(in) :: Idata1
      integer, intent(in) :: NVAR
      integer, intent(in) :: NEQ
      integer, dimension(NVAR,NEQ), intent(out) :: Idata2

      integer :: ieq

      do ieq = 1, NEQ
        Idata2(:,ieq) = Idata1(ieq)
      end do

    end subroutine

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_spreadMatrix(rmatrix1, rmatrix2)

!<description>
    ! This subroutine spreads a scalar matrix which is not stored in interleave
    ! format into another matrix which is stored in interleave format.
!</description>

!<input>
    ! Scalar source matrix
    type(t_matrixScalar), intent(in) :: rmatrix1
!</input>

!<inputoutput>
    ! Scalar destination matrix in interleave format
    type(t_matrixScalar), intent(inout) :: rmatrix2
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2
    real(SP), dimension(:), pointer :: p_Fdata1,p_Fdata2
    integer, dimension(:), pointer :: p_Idata1,p_Idata2

    ! Source matrices must not be stored in interleave format
    if (rmatrix1%NVAR .ne. 1) then
      call output_line('Source matrix must not be stored in interleave format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadMatrix')
      call sys_halt()
    end if

    ! Check if matrices are compatible
    call lsyssc_isMatrixMatrixCompatible(rmatrix1,rmatrix2)

    ! Ok, now we can copy the matrices
    select case (rmatrix2%cinterleavematrixFormat)

    case (LSYSSC_MATRIXUNDEFINED)
      ! Destination matrix is identical to source matrix
      call lsyssc_duplicateMatrix (rmatrix1,rmatrix2,&
          LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE)

    case (LSYSSC_MATRIX1)
      select case (rmatrix1%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix1, p_Ddata1)
        call lsyssc_getbase_double(rmatrix2, p_Ddata2)
        call do_spreadDble(p_Ddata1, rmatrix2%NVAR, rmatrix2%NVAR,&
                           rmatrix2%NA, p_Ddata2)

      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix1, p_Fdata1)
        call lsyssc_getbase_single(rmatrix2, p_Fdata2)
        call do_spreadSngl(p_Fdata1, rmatrix2%NVAR, rmatrix2%NVAR,&
                           rmatrix2%NA, p_Fdata2)

      case (ST_INT)
        call lsyssc_getbase_int(rmatrix1, p_Idata1)
        call lsyssc_getbase_int(rmatrix2, p_Idata2)
        call do_spreadInt(p_Idata1, rmatrix2%NVAR, rmatrix2%NVAR,&
                          rmatrix2%NA, p_Idata2)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadMatrix')
        call sys_halt()
      end select

    case (LSYSSC_MATRIXD)
      select case (rmatrix1%cdataType)
      case (ST_DOUBLE)
        call lsyssc_getbase_double(rmatrix1, p_Ddata1)
        call lsyssc_getbase_double(rmatrix2, p_Ddata2)
        call do_spreadDble(p_Ddata1, rmatrix2%NVAR, 1,&
                           rmatrix2%NA, p_Ddata2)

      case (ST_SINGLE)
        call lsyssc_getbase_single(rmatrix1, p_Fdata1)
        call lsyssc_getbase_single(rmatrix2, p_Fdata2)
        call do_spreadSngl(p_Fdata1, rmatrix2%NVAR, 1,&
                           rmatrix2%NA, p_Fdata2)

      case (ST_INT)
        call lsyssc_getbase_int(rmatrix1, p_Idata1)
        call lsyssc_getbase_int(rmatrix2, p_Idata2)
        call do_spreadInt(p_Idata1, rmatrix2%NVAR, 1,&
                          rmatrix2%NA, p_Idata2)

      case default
        call output_line('Unsupported data type!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadMatrix')
        call sys_halt()
      end select

    case default
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_spreadMatrix')
      call sys_halt()
    end select

  contains

    ! Here, the real working routines follow.

    !**************************************************************

    subroutine do_spreadDble(Ddata1, NVAR, MVAR, NA, Ddata2)
      real(DP), dimension(:), intent(in) :: Ddata1
      integer, intent(in) :: NVAR,MVAR
      integer, intent(in) :: NA
      real(DP), dimension(NVAR,MVAR,NA), intent(out) :: Ddata2

      integer :: ia

      do ia = 1, NA
        Ddata2(:,:,ia) = Ddata1(ia)
      end do

    end subroutine

    !**************************************************************

    subroutine do_spreadSngl(Fdata1, NVAR, MVAR, NA, Fdata2)
      real(SP), dimension(:), intent(in) :: Fdata1
      integer, intent(in) :: NVAR,MVAR
      integer, intent(in) :: NA
      real(SP), dimension(NVAR,MVAR,NA), intent(out) :: Fdata2

      integer :: ia

      do ia = 1, NA
        Fdata2(:,:,ia) = Fdata1(ia)
      end do

    end subroutine

    !**************************************************************

    subroutine do_spreadInt(Idata1, NVAR, MVAR, NA, Idata2)
      integer, dimension(:), intent(in) :: Idata1
      integer, intent(in) :: NVAR,MVAR
      integer, intent(in) :: NA
      integer, dimension(NVAR,MVAR,NA), intent(out) :: Idata2

      integer :: ia

      do ia = 1, NA
        Idata2(:,:,ia) = Idata1(ia)
      end do

    end subroutine

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_packVector(rvector1, rvector2, ivar)

!<description>
    ! This subroutine packs a scalar vector which is stored in interleave
    ! format into another vector which is not stored in interleave format.
!</description>

!<input>
    ! Scalar source vector in interleave format
    type(t_vectorScalar), intent(in) :: rvector1

    ! Number of the variable to pack
    integer, intent(in) :: ivar
!</input>

!<inputoutput>
    ! Scalar destination vector
    type(t_vectorScalar), intent(inout) :: rvector2
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2
    real(SP), dimension(:), pointer :: p_Fdata1,p_Fdata2
    integer, dimension(:), pointer :: p_Idata1,p_Idata2

    ! Source vector must be stored in interleave format
    if (rvector1%NVAR .lt. ivar) then
      call output_line('Source vector does not provide variable IVAR!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_packVector')
      call sys_halt()
    end if

    ! Destination vector must not be stored in interleave format
    if (rvector2%NVAR .ne. 1) then
      call output_line('Destination vector must not be stored in interleave format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_packVector')
      call sys_halt()
    end if

    ! Vectors must have the same size
    if (rvector1%NEQ .ne. rvector2%NEQ) then
      call output_line('Vectors not compatible, different size!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_packVector')
      call sys_halt()
    end if

    ! Vectors must have the data type
    if (rvector1%cdataType .ne. rvector2%cdataType) then
      call output_line('Vectors not compatible, different data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_packVector')
      call sys_halt()
    end if

    ! isortStrategy < 0 means unsorted. Both unsorted is ok.

    ! Check the sorting
    if (rvector1%bisSorted .neqv. rvector2%bisSorted) then
      call output_line("Vectors not compatible, differently sorted!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_packVector")
      call sys_halt()
    end if
    
    if (rvector1%bisSorted) then
      if (rvector1%p_rsortStrategy%ctype .ne. rvector2%p_rsortStrategy%ctype) then
        call output_line("Vectors not compatible, differently sorted!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_packVector")
        call sys_halt()
      end if
    end if

    select case (rvector1%cdataType)
    case (ST_DOUBLE)
      call lsyssc_getbase_double(rvector1, p_Ddata1)
      call lsyssc_getbase_double(rvector2, p_Ddata2)
      call do_packDble(p_Ddata1, rvector1%NVAR, rvector1%NEQ, ivar, p_Ddata2)

    case (ST_SINGLE)
      call lsyssc_getbase_single(rvector1, p_Fdata1)
      call lsyssc_getbase_single(rvector2, p_Fdata2)
      call do_packSngl(p_Fdata1, rvector1%NVAR, rvector1%NEQ, ivar, p_Fdata2)

    case (ST_INT)
      call lsyssc_getbase_int(rvector1, p_Idata1)
      call lsyssc_getbase_int(rvector2, p_Idata2)
      call do_packInt(p_Idata1, rvector1%NVAR, rvector1%NEQ, ivar, p_Idata2)

    case default
      call output_line('Unsupported data type!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_packVector')
      call sys_halt()
    end select

  contains

    ! Here, the real working routines follow.

    !**************************************************************

    subroutine do_packDble(Ddata1, NVAR, NEQ, ivar, Ddata2)
      real(DP), dimension(NVAR,NEQ), intent(in) :: Ddata1
      integer, intent(in) :: NVAR
      integer, intent(in) :: NEQ
      integer, intent(in) :: ivar
      real(DP), dimension(:), intent(out) :: Ddata2

      integer :: ieq

      do ieq = 1, NEQ
        Ddata2(ieq) = Ddata1(ivar,ieq)
      end do

    end subroutine

    !**************************************************************

    subroutine do_packSngl(Fdata1, NVAR, NEQ, ivar, Fdata2)
      real(SP), dimension(NVAR,NEQ), intent(in) :: Fdata1
      integer, intent(in) :: NVAR
      integer, intent(in) :: NEQ
      integer, intent(in) :: ivar
      real(SP), dimension(:), intent(out) :: Fdata2

      integer :: ieq

      do ieq = 1, NEQ
        Fdata2(ieq) = Fdata1(ivar,ieq)
      end do

    end subroutine

    !**************************************************************

    subroutine do_packInt(Idata1, NVAR, NEQ, ivar, Idata2)
      integer, dimension(NVAR,NEQ), intent(in) :: Idata1
      integer, intent(in) :: NVAR
      integer, intent(in) :: NEQ
      integer, intent(in) :: ivar
      integer, dimension(:), intent(out) :: Idata2

      integer :: ieq

      do ieq = 1, NEQ
        Idata2(ieq) = Idata1(ivar,ieq)
      end do

    end subroutine

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine lsyssc_createFpdbObjectVec (rfpdbObjectItem, sname, rvector)

!<description>
    ! This subroutine creates an abstract ObjectItem from the scalar
    ! vector that can be stored in the persistence database
!</description>

!<input>
    ! The full qualified name of the object
    character(LEN=*), intent(in) :: sname
!</input>

!<inputoutput>
    ! The ObjectItem  that is created
    type(t_fpdbObjectItem), intent(inout) :: rfpdbObjectItem

    ! The scalar vector
    type(t_vectorScalar), intent(inout), target :: rvector
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fpdbDataItem), pointer :: p_fpdbDataItem

    ! Check if vector has UUID; otherwise create one
    if (uuid_isNil(rvector%ruuid)) then
      call uuid_createUUID(4, rvector%ruuid)
    end if
    rfpdbObjectItem%ruuid = rvector%ruuid

    ! Set the name and type of the object
    rfpdbObjectItem%sname = sname
    rfpdbObjectItem%stype = 't_vectorScalar'

    ! Allocate the array of data items: 11
    allocate(rfpdbObjectItem%p_RfpdbDataItem(11))

    ! Fill the array of data items: NEQ
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(1)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'NEQ'
    p_fpdbDataItem%iinteger = rvector%NEQ

    ! Fill the array of data items: NVAR
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(2)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'NVAR'
    p_fpdbDataItem%iinteger = rvector%NVAR

    ! Fill the array of data items: cdataType
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(3)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'cdataType'
    p_fpdbDataItem%iinteger = rvector%cdataType

    ! Fill the array of data items: h_Ddata
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(4)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'h_Ddata'
    p_fpdbDataItem%iinteger = rvector%h_Ddata

    ! Fill the array of data items: isortStrategy
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(5)
    p_fpdbDataItem%ctype    = FPDB_LOGICAL
    p_fpdbDataItem%sname    = 'bisSorted'
    p_fpdbDataItem%blogical = rvector%bisSorted

    ! Fill the array of data items: h_Ipermutation
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(6)
    p_fpdbDataItem%sname    = 'p_rsortStrategy'
    if (associated(rvector%p_rsortStrategy)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if

    ! Fill the array of data items: iidxFirstEntry
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(7)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'iidxFirstEntry'
    p_fpdbDataItem%iinteger = rvector%iidxFirstEntry

    ! Fill the array of data items: ITags
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(8)
    p_fpdbDataItem%ctype        =  FPDB_INT1D
    p_fpdbDataItem%sname        =  'ITags'
    p_fpdbDataItem%p_Iinteger1D => rvector%ITags

    ! Fill the array of data items: DTags
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(9)
    p_fpdbDataItem%ctype       =  FPDB_DOUBLE1D
    p_fpdbDataItem%sname       =  'DTags'
    p_fpdbDataItem%p_Ddouble1D => rvector%DTags

    ! Fill the array of data items: bisCopy
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(10)
    p_fpdbDataItem%ctype    = FPDB_LOGICAL
    p_fpdbDataItem%sname    = 'bisCopy'
    p_fpdbDataItem%blogical = rvector%bisCopy

    ! Fill the array of data items: p_rspatialDiscr
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(11)
    p_fpdbDataItem%sname = 'p_rspatialDiscr'
    if (associated(rvector%p_rspatialDiscr)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine lsyssc_createFpdbObjectMat (rfpdbObjectItem, sname, rmatrix)

!<description>
    ! This subroutine creates an abstract ObjectItem from the scalar
    ! matrix that can be stored in the persistence database
!</description>

!<input>
    ! The full qualified name of the object
    character(LEN=*), intent(in) :: sname
!</input>

!<inputoutput>
    ! The ObjectItem  that is created
    type(t_fpdbObjectItem), intent(inout) :: rfpdbObjectItem

    ! The scalar matrix
    type(t_matrixScalar), intent(inout), target :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fpdbDataItem), pointer :: p_fpdbDataItem

    ! Check if matrix has UUID; otherwise create one
    if (uuid_isNil(rmatrix%ruuid)) then
      call uuid_createUUID(4, rmatrix%ruuid)
    end if
    rfpdbObjectItem%ruuid = rmatrix%ruuid

    ! Set the name and type of the object
    rfpdbObjectItem%sname = sname
    rfpdbObjectItem%stype = 't_matrixScalar'

    ! Allocate the array of data items: 22
    allocate(rfpdbObjectItem%p_RfpdbDataItem(22))

    ! Fill the array of data items: cmatrixFormat
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(1)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'cmatrixFormat'
    p_fpdbDataItem%iinteger = rmatrix%cmatrixFormat

    ! Fill the array of data items: cinterleavematrixFormat
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(2)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'cinterleavematrixFormat'
    p_fpdbDataItem%iinteger = rmatrix%cinterleavematrixFormat

    ! Fill the array of data items: imatrixSpec
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(3)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'imatrixSpec'
    p_fpdbDataItem%iinteger = rmatrix%imatrixSpec

    ! Fill the array of data items: NA
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(4)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'NA'
    p_fpdbDataItem%iinteger = rmatrix%NA

    ! Fill the array of data items: NEQ
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(5)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'NEQ'
    p_fpdbDataItem%iinteger = rmatrix%NEQ

    ! Fill the array of data items: NCOLS
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(6)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'NCOLS'
    p_fpdbDataItem%iinteger = rmatrix%NCOLS

    ! Fill the array of data items: NVAR
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(7)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'NVAR'
    p_fpdbDataItem%iinteger = rmatrix%NVAR

    ! Fill the array of data items: dscaleFactor
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(8)
    p_fpdbDataItem%ctype   = FPDB_DOUBLE
    p_fpdbDataItem%sname   = 'dscaleFactor'
    p_fpdbDataItem%ddouble = rmatrix%dscaleFactor

    ! Fill the array of data items: bcolumnsSorted
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(9)
    p_fpdbDataItem%ctype    = FPDB_LOGICAL
    p_fpdbDataItem%sname    = 'bcolumnsSorted'
    p_fpdbDataItem%blogical = rmatrix%bcolumnsSorted

    ! Fill the array of data items: browsSorted
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(10)
    p_fpdbDataItem%ctype    = FPDB_LOGICAL
    p_fpdbDataItem%sname    = 'browsSorted'
    p_fpdbDataItem%blogical = rmatrix%browsSorted

    ! Fill the array of data items: p_rsortStrategyColumns
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(11)
    p_fpdbDataItem%sname    = 'p_rsortStrategyColumns'
    if (associated(rmatrix%p_rsortStrategyColumns)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if

    ! Fill the array of data items: p_rsortStrategyRows
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(12)
    p_fpdbDataItem%sname    = 'p_rsortStrategyRows'
    if (associated(rmatrix%p_rsortStrategyRows)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if

    ! Fill the array of data items: cdataType
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(13)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'cdataType'
    p_fpdbDataItem%iinteger = rmatrix%cdataType

    ! Fill the array of data items: h_Da
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(14)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'h_Da'
    p_fpdbDataItem%iinteger = rmatrix%h_Da

    ! Fill the array of data items: h_Kcol
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(15)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'h_Kcol'
    p_fpdbDataItem%iinteger = rmatrix%h_Kcol

    ! Fill the array of data items: h_Kld
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(16)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'h_Kld'
    p_fpdbDataItem%iinteger = rmatrix%h_Kld

    ! Fill the array of data items: h_Kdiagonal
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(17)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'h_Kdiagonal'
    p_fpdbDataItem%iinteger = rmatrix%h_Kdiagonal

    ! Fill the array of data items: ITags
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(18)
    p_fpdbDataItem%ctype        =  FPDB_INT1D
    p_fpdbDataItem%sname        =  'ITags'
    p_fpdbDataItem%p_Iinteger1D => rmatrix%ITags

    ! Fill the array of data items: DTags
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(19)
    p_fpdbDataItem%ctype       =  FPDB_DOUBLE1D
    p_fpdbDataItem%sname       =  'DTags'
    p_fpdbDataItem%p_Ddouble1D => rmatrix%DTags

    ! Fill the array of data items: bidenticalTrialAndTest
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(20)
    p_fpdbDataItem%ctype    = FPDB_LOGICAL
    p_fpdbDataItem%sname    = 'bidenticalTrialAndTest'
    p_fpdbDataItem%blogical = rmatrix%bidenticalTrialAndTest

    ! Fill the array of data items: p_rspatialDiscrTrial
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(21)
    p_fpdbDataItem%sname = 'p_rspatialDiscrTrial'
    if (associated(rmatrix%p_rspatialDiscrTrial)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if


    ! Fill the array of data items: p_rspatialDiscrTest
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(22)
    p_fpdbDataItem%sname = 'p_rspatialDiscrTest'
    if (associated(rmatrix%p_rspatialDiscrTest)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine lsyssc_restoreFpdbObjectVec (rfpdbObjectItem, rvector)

!<description>
    ! This subroutine restores the scalar vector from the abstract ObjectItem
!</description>

!<input>
    ! The object item that is created
    type(t_fpdbObjectItem), intent(in) :: rfpdbObjectItem
!</input>

!<output>
    ! The scalar vector
    type(t_vectorScalar), intent(out), target :: rvector
!</output>
!</subroutine>

    ! local variables
    type(t_fpdbDataItem), pointer :: p_fpdbDataItem

    ! Check if ObjectItem has correct type
    if (trim(rfpdbObjectItem%stype) .ne. "t_vectorScalar") then
      call output_line ("Invalid object type!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    end if

    ! Check if DataItems are associated
    if (.not.associated(rfpdbObjectItem%p_RfpdbDataItem)) then
      call output_line ("Missing data!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    end if

    ! Check if DataItems have correct size
    if (size(rfpdbObjectItem%p_RfpdbDataItem) .ne. 11) then
      call output_line ("Invalid data!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    end if

    ! Restore the UUID
    rvector%ruuid = rfpdbObjectItem%ruuid

    ! Restore the data from the DataItem: NEQ
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(1)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "NEQ")) then
      call output_line ("Invalid data: NEQ!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    else
      rvector%NEQ = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: NVAR
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(2)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "NVAR")) then
      call output_line ("Invalid data: NVAR!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    else
      rvector%NVAR = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: cdataType
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(3)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "cdataType")) then
      call output_line ("Invalid data: cdataType!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    else
      rvector%cdataType = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: h_Ddata
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(4)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "h_Ddata")) then
      call output_line ("Invalid data: h_Ddata!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    else
      rvector%h_Ddata = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: isortStrategy
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(5)
    if ((p_fpdbDataItem%ctype .ne. FPDB_LOGICAL) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "bisSorted")) then
      call output_line ("Invalid data: bisSorted!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    else
      rvector%bisSorted = p_fpdbDataItem%blogical
    end if

    ! Restore the data from the DataItem: h_Ipermutation
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(6)
    if (trim(p_fpdbDataItem%sname) .ne. "p_rsortStrategy") then
      call output_line ("Invalid data: p_rsortStrategy!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rvector%p_rsortStrategy)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ("Invalid data: p_rsortingStrategy!", &
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
        call sys_halt()
      end if
    end if

    ! Restore the data from the DataItem: iidxFirstEntry
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(7)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "iidxFirstEntry")) then
      call output_line ("Invalid data: iidxFirstEntry!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    else
      rvector%iidxFirstEntry = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: ITags
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(8)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT1D) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "ITags")) then
      call output_line ("Invalid data: ITags!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    else
      call fpdb_getdata_int1d(p_fpdbDataItem, rvector%ITags)
    end if

    ! Restore the data from the DataItem: DTags
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(9)
    if ((p_fpdbDataItem%ctype .ne. FPDB_DOUBLE1D) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "DTags")) then
      call output_line ("Invalid data: DTags!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    else
      call fpdb_getdata_double1d(p_fpdbDataItem, rvector%DTags)
    end if

    ! Restore the data from the DataItem: bisCopy
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(10)
    if ((p_fpdbDataItem%ctype .ne. FPDB_LOGICAL) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "bisCopy")) then
      call output_line ("Invalid data: bisCopy!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    else
      rvector%bisCopy = p_fpdbDataItem%blogical
    end if

    ! Restore the data from the DataItem: p_rspatialDiscr
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(11)
    if (trim(p_fpdbDataItem%sname) .ne. "p_rspatialDiscr") then
      call output_line ("Invalid data: p_rspatialDiscr!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rvector%p_rspatialDiscr)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ("Invalid data: p_rspatialDiscr!", &
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectVec")
        call sys_halt()
      end if
    end if

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine lsyssc_restoreFpdbObjectMat (rfpdbObjectItem, rmatrix)

!<description>
    ! This subroutine restores the scalar matrix from the abstract ObjectItem
!</description>

!<input>
    ! The object item that is created
    type(t_fpdbObjectItem), intent(in) :: rfpdbObjectItem
!</input>

!<output>
    ! The scalar matrix
    type(t_matrixScalar), intent(out), target :: rmatrix
!</output>
!</subroutine>

    ! local variables
    type(t_fpdbDataItem), pointer :: p_fpdbDataItem

    ! Check if ObjectItem has correct type
    if (trim(rfpdbObjectItem%stype) .ne. "t_matrixScalar") then
      call output_line ("Invalid object type!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    end if

    ! Check if DataItems are associated
    if (.not.associated(rfpdbObjectItem%p_RfpdbDataItem)) then
      call output_line ("Missing data!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    end if

    ! Check if DataItems have correct size
    if (size(rfpdbObjectItem%p_RfpdbDataItem) .ne. 20) then
      call output_line ("Invalid data!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    end if

    ! Restore the UUID
    rmatrix%ruuid = rfpdbObjectItem%ruuid

    ! Restore the data from the DataItem: cmatrixFormat
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(1)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "cmatrixFormat")) then
      call output_line ("Invalid data: cmatrixFormat!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%cmatrixFormat = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: cinterleavematrixFormat
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(2)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "cinterleavematrixFormat")) then
      call output_line ("Invalid data: cinterleavematrixFormat!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%cinterleavematrixFormat = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: imatrixSpec
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(3)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "imatrixSpec")) then
      call output_line ("Invalid data: imatrixSpec!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%imatrixSpec = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: NA
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(4)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "NA")) then
      call output_line ("Invalid data: NA!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%NA = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: NEQ
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(5)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "NEQ")) then
      call output_line ("Invalid data: NEQ!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%NEQ = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: NCOLS
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(6)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "NCOLS")) then
      call output_line ("Invalid data: NCOLS!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%NCOLS = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: NVAR
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(7)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "NVAR")) then
      call output_line ("Invalid data: NVAR!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%NVAR = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: dscaleFactor
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(8)
    if ((p_fpdbDataItem%ctype .ne. FPDB_DOUBLE) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "dscaleFactor")) then
      call output_line ("Invalid data: dscaleFactor!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%dscaleFactor = p_fpdbDataItem%ddouble
    end if

    ! Restore the data from the DataItem: bcolumnsSorted
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(9)
    if ((p_fpdbDataItem%ctype .ne. FPDB_LOGICAL) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "bcolumnsSorted")) then
      call output_line ("Invalid data: bcolumnsSorted!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%bcolumnsSorted = p_fpdbDataItem%blogical
    end if

    ! Restore the data from the DataItem: browsSorted
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(10)
    if ((p_fpdbDataItem%ctype .ne. FPDB_LOGICAL) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "browsSorted")) then
      call output_line ("Invalid data: browsSorted!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%browsSorted = p_fpdbDataItem%blogical
    end if

    ! Restore the data from the DataItem: p_rsortStrategy
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(11)
    if (trim(p_fpdbDataItem%sname) .ne. "p_rsortStrategyColumns") then
      call output_line ("Invalid data: p_rsortStrategyColumns!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rmatrix%p_rsortStrategyColumns)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ("Invalid data: p_rsortingStrategy!", &
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
        call sys_halt()
      end if
    end if

    ! Restore the data from the DataItem: p_rsortStrategy
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(12)
    if (trim(p_fpdbDataItem%sname) .ne. "p_rsortStrategyRows") then
      call output_line ("Invalid data: p_rsortStrategyRows!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rmatrix%p_rsortStrategyRows)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ("Invalid data: p_rsortingStrategyRows!", &
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
        call sys_halt()
      end if
    end if

    ! Restore the data from the DataItem: cdataType
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(13)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "cdataType")) then
      call output_line ("Invalid data: cdataType!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%cdataType = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: h_Da
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(14)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "h_Da")) then
      call output_line ("Invalid data: h_Da!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%h_Da = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: h_Kcol
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(15)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "h_Kcol")) then
      call output_line ("Invalid data: h_Kcol!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%h_Kcol = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: h_Kld
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(16)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "h_Kld")) then
      call output_line ("Invalid data: h_Kld!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%h_Kld = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: h_Kdiagonal
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(17)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "h_Kdiagonal")) then
      call output_line ("Invalid data: h_Kdiagonal!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%h_Kdiagonal = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: ITags
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(18)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT1D) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "ITags")) then
      call output_line ("Invalid data: ITags!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      call fpdb_getdata_int1d(p_fpdbDataItem, rmatrix%ITags)
    end if

    ! Restore the data from the DataItem: DTags
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(19)
    if ((p_fpdbDataItem%ctype .ne. FPDB_DOUBLE1D) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "DTags")) then
      call output_line ("Invalid data: DTags!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      call fpdb_getdata_double1d(p_fpdbDataItem, rmatrix%DTags)
    end if

    ! Restore the data from the DataItem: bidenticalTrialAndTest
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(20)
    if ((p_fpdbDataItem%ctype .ne. FPDB_LOGICAL) .or.&
        (trim(p_fpdbDataItem%sname) .ne. "bidenticalTrialAndTest")) then
      call output_line ("Invalid data: bidenticalTrialAndTest!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      rmatrix%bidenticalTrialAndTest = p_fpdbDataItem%blogical
    end if

    ! Restore the data from the DataItem: p_rspatialDiscrTrial
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(21)
    if (trim(p_fpdbDataItem%sname) .eq. "p_rspatialDiscrTrial") then
      call output_line ("Invalid data: p_rspatialDiscrTrial!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rmatrix%p_rspatialDiscrTrial)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ("Invalid data: p_rspatialDiscrTrial!", &
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
        call sys_halt()
      end if
    end if

    ! Restore the data from the DataItem: p_rspatialDiscrTest
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(22)
    if (trim(p_fpdbDataItem%sname) .ne. "p_rspatialDiscrTest") then
      call output_line ("Invalid data: p_rspatialDiscrTest!", &
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rmatrix%p_rspatialDiscrTest)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ("Invalid data: p_rspatialDiscrTest!", &
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_restoreFpdbObjectMat")
        call sys_halt()
      end if
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_unshareMatrix (rmatrix,bstructure,bdata)

!<description>
  ! Resets the sharing state of the structure/data arrays in a matrix.
  ! bstructure/bdata decides on if the structure and/or the data is
  ! associated to rmatrix. If the matrix is already the owner of the
  ! structure/data, nothing happens. If the matrix shares its structure/data
  ! with another matrix, the new memory is allocated, the old data is
  ! copied and the new independent structure/data arrays are associated
  ! to the matrix. Therefore, if bstructure=bdata=.true, the matrix
  ! will be made completely independent.
!</description>

!<input>
  ! Revoke the sharing state of the matrix structure. If the matrix
  ! is not the owner of the structure, a new structure array is allocated
  ! in memory and the data of the old structure is copied.
  logical, intent(in) :: bstructure

  ! Revoke the sharing state of the matrix data. If the matrix
  ! is not the owner of the data, a new data array is allocated
  ! in memory and the data of the old structure is copied.
  logical, intent(in) :: bdata
!</input>

!<inputoutput>
  ! Matrix to be changed matrix.
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    type(t_matrixScalar) :: rmatrix2
    integer :: idup1,idup2

    idup1 = LSYSSC_DUP_SHARE
    idup2 = LSYSSC_DUP_SHARE

    if (bstructure .and. &
        (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY) .ne. 0)) then
      idup1 = LSYSSC_DUP_COPY
    end if

    if (bdata .and. &
        (iand(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY) .ne. 0)) then
      idup2 = LSYSSC_DUP_COPY
    end if

    ! Probably nothing to do.
    if ((idup1 .eq. LSYSSC_DUP_COPY) .and. (idup2 .eq. LSYSSC_DUP_COPY)) return

    ! Duplicate the matrix to eventually allocate memory
    call lsyssc_duplicateMatrix (rmatrix,rmatrix2,idup1,idup2)

    ! Next thing is to release the old matrix; but that is not as easy
    ! as calling lsyssc_releaseMatrix! The point is that the source
    ! matrix may or may not share its structure/content with another
    ! one -- and the new matrix may do this as well!
    !
    ! We have to handle structure and content separately.
    ! If the structure/content is to be made independent, we have to reset
    ! the copy flag. Otherwise, we have to take the old copy flag.
    rmatrix2%imatrixSpec = rmatrix%imatrixSpec

    if (idup1 .eq. LSYSSC_DUP_COPY) then
      rmatrix2%imatrixSpec = iand(rmatrix2%imatrixSpec,not(LSYSSC_MSPEC_STRUCTUREISCOPY))
    else
      ! As rmatrix2 is now the owner, prevent the data in rmatrix from being released.
      rmatrix%imatrixSpec = ior(rmatrix%imatrixSpec,LSYSSC_MSPEC_STRUCTUREISCOPY)
    end if

    if (idup2 .eq. LSYSSC_DUP_COPY) then
      rmatrix2%imatrixSpec = iand(rmatrix2%imatrixSpec,not(LSYSSC_MSPEC_CONTENTISCOPY))
    else
      ! As rmatrix2 is now the owner, prevent the data in rmatrix from being released.
      rmatrix%imatrixSpec = ior(rmatrix%imatrixSpec,LSYSSC_MSPEC_CONTENTISCOPY)
    end if

    ! Release the old matrix
    call lsyssc_releaseMatrix(rmatrix)

    ! Create a hard-copy to rmatrix.
    ! This is one of the very rare cases where we on purpose assign a matrix
    ! structure to another...
    rmatrix = rmatrix2

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_unshareVector (rvector)

!<description>
  ! Resets the sharing state of the data arrays in a vector.
  ! If the vector is already the owner of the
  ! structure/data, nothing happens. If the vector shares its structure/data
  ! with another vector, the new memory is allocated, the old data is
  ! copied and the new independent data arrays are associated
  ! to the vector.
!</description>

!<output>
  ! Vector to be changed matrix.
  type(t_vectorScalar), intent(inout) :: rvector
!</output>

!</subroutine>

    type(t_vectorScalar) :: rvector2

    ! If the vector is not the owner, make it the owner.
    if (rvector%bisCopy) then
      ! Copy & unshare the data
      call lsyssc_duplicateVector(rvector,rvector2,LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      ! Release the original vector
      call lsyssc_releaseVector (rvector)

      ! Create a hard-copy of rvector.
      ! This is one of the very rare cases where we on purpose assign a matrix
      ! structure to another...
      rvector = rvector2
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_moveMatrix (rsourceMatrix,rdestMatrix)

!<description>
  ! Moves a matrix from a source matrix to a destination matrix.
  ! The source matrix will remain as 'shared copy' of the destination matrix
  ! and can safely be released.
!</description>

!<inputoutput>
  ! Source matrix
  type(t_matrixScalar),intent(inout) :: rsourceMatrix

  ! Destination matrix
  type(t_matrixScalar),intent(inout) :: rdestMatrix
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: cflags1,cflags2

    ! Duplicate the matrix, share the data.
    call lsyssc_duplicateMatrix (rsourceMatrix, rdestMatrix,&
        LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

    ! Change the ownership from the source to the destination matrix
    cflags1 = iand (rsourceMatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
    cflags2 = iand (rdestMatrix%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
    rsourceMatrix%imatrixSpec = &
      ior(iand(rsourceMatrix%imatrixSpec,not(LSYSSC_MSPEC_ISCOPY)),cflags2)
    rdestMatrix%imatrixSpec = &
      ior(iand(rdestMatrix%imatrixSpec,not(LSYSSC_MSPEC_ISCOPY)),cflags1)

    ! Now, the destination matrix is the owner and the source matrix
    ! the copy...

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_createEmptyMatrixStub (rmatrix,cmatrixFormat,neq,ncol)

!<description>
  ! Creates an empty matrix in format cmatrixFormat without allocating memory.
  ! The caller is responsible for filling in missing information.
!</description>

!<input>
  ! Matrix format
  integer, intent(in) :: cmatrixFormat

  ! Number of rows
  integer, intent(in) :: neq

  ! OPTIONAL: Number of columns. If not present, defaults to neq.
  integer, intent(in), optional :: ncol
!</input>

!<output>
  ! Empty matrix. No data has been allocated.
  type(t_matrixScalar), intent(out) :: rmatrix
!</output>

!</subroutine>

    rmatrix%cmatrixFormat = cmatrixFormat
    rmatrix%neq = neq
    if (present(ncol)) then
      rmatrix%NCOLS = ncol
    else
      rmatrix%NCOLS = neq
    end if

  end subroutine

  ! ***************************************************************************
!<subroutine>

  subroutine lsyssc_createEmptyMatrix9 (rmatrix,neq,na,ncols)

!<description>
  ! Creates an empty matrix in matrix format 9 (CSR). The matrix is designed
  ! to have na entries, neq rows and ncols columns.
!</description>

!<input>
  ! Number of rows.
  integer, intent(in) :: neq

  ! Number of entries in the matrix.
  integer, intent(in) :: na

  ! OPTIONAL: Number of columns. If not present, NEQ is assumed.
  integer, intent(in), optional :: ncols
!</input>

!<output>
  ! Matrix to create.
  type(t_matrixScalar), intent(out), target :: rmatrix
!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Kld

    ! Create basic matrix stub.
    call lsyssc_createEmptyMatrixStub (rmatrix,LSYSSC_MATRIX9,neq,ncols)

    ! Allocate memory
    rmatrix%na = na

    call storage_new ('lsyssc_createMatrixFormat9', 'Da', &
        rmatrix%na, ST_DOUBLE, rmatrix%h_Da, ST_NEWBLOCK_ZERO)
    call storage_new ('lsyssc_createMatrixFormat9', 'Kcol', &
        rmatrix%na, ST_INT, rmatrix%h_Kcol, ST_NEWBLOCK_ZERO)
    call storage_new ('lsyssc_createMatrixFormat9', 'Kld', &
        neq+1, ST_INT, rmatrix%h_Kld, ST_NEWBLOCK_ZERO)
    call storage_new ('lsyssc_createMatrixFormat9', 'Kdiagonal', &
        neq, ST_INT, rmatrix%h_Kdiagonal, ST_NEWBLOCK_ZERO)

    ! Element na+1 in Kld must be = 1 since there are no elements
    ! in the matrix.
    call lsyssc_getbase_Kld (rmatrix,p_Kld)
    p_Kld(neq+1) = 1

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_setRowMatrix9 (rmatrix,irow,nentries,Kcol,Da)

!<description>
  ! Replaces row irow by the column and matrix data specified in
  ! Kcol and Da. If the row does not exist, it is created.
!</description>

!<input>
  ! Row to modify
  integer, intent(in) :: irow

  ! Number of entries, the row should have.
  ! If the row already exists, a value "-1" replaces the existing
  ! entries in the row, the length is automatically calculated.
  integer, intent(in) :: nentries

  ! Column numbers that should be added.
  integer, dimension(:), intent(in) :: Kcol

  ! OPTIONAL: Matrix data to be added.
  real(DP), dimension(:), intent(in), optional :: Da
!</input>

!<output>
  ! Matrix to create.
  type(t_matrixScalar), intent(inout), target :: rmatrix
!</output>

!</subroutine>

    ! local variables
    integer :: i,n
    integer, dimension(:), pointer :: p_Kcol, p_Kld, P_Kdiagonal
    real(DP), dimension(:), pointer :: p_Da

    ! Get pointers to the matrix data.
    call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_Kld (rmatrix,p_Kld)
    call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)

    ! Does the row already exist?
    if ((p_Kld(irow) .gt. 0) .and. (p_Kld(irow+1) .ge. p_Kld(irow))) then
      ! Does the number of entries match?
      n = nentries
      if (n .lt. 0) then
        n = p_Kld(irow+1)-p_Kld(irow+1)
      else
        if (n .ne. (p_Kld(irow+1)-p_Kld(irow))) then
          call output_line('Number of entries in the existing row does not match nentries!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_setRowMatrix9')
          call sys_halt()
        end if
      end if
    else
      n = nentries
      if (n .lt. 0) then
        call output_line('Number of entried that should be put into the row not specified.',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_setRowMatrix9')
        call sys_halt()
      end if

      if ((irow .gt. 1) .and. (p_Kld(irow) .eq. 0)) then
        call output_line('Insertion at arbitrary position currently not supported!'//&
                         ' Only appending allowed.',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_setRowMatrix9')
        call sys_halt()
      end if

      ! Create the row.
      p_Kld(irow) = p_Kld(rmatrix%neq+1)
      p_Kld(irow+1) = p_Kld(rmatrix%neq+1) + n
      p_Kld(rmatrix%neq+1) = p_Kld(rmatrix%neq+1) + n

      if (irow .eq. rmatrix%neq-1) then
        ! Last but one row. set up the last entry in KLD appropriately.
        p_Kld(rmatrix%neq+1) = rmatrix%na+1
      end if

    end if

    ! Copy column data into the row
    p_Kcol(p_Kld(irow):p_Kld(irow+1)-1) = Kcol(1:n)

    ! Find the diagonal
    p_Kdiagonal(irow) = p_Kld(irow+1)
    do i=1,n
      if (Kcol(i) .ge. irow) then
        p_Kdiagonal(irow) = p_Kld(irow)+i-1
        exit
      end if
    end do

    ! Probably copy matrix data.
    if (present(Da)) then

      if (rmatrix%cdataType .ne. ST_DOUBLE) then
        call output_line('Only double precision matrices supported!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_setRowMatrix9')
        call sys_halt()
      end if

      call lsyssc_getbase_double (rmatrix,p_Da)
      p_Da(p_Kld(irow):p_Kld(irow+1)-1) = Da(1:n)

    end if

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_sortMatrixCols(rmatrix)

!<description>
  ! Sorts the columns on a format 9 matrix to be in ascending order.
  ! This routine also recalculates the diagonal pointer array of the matrix.
!</description>

!<inputoutput>
  ! The matrix that is to be sorted.
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

  integer, dimension(:), pointer :: p_IrowPtr, p_IcolIdx, p_Idiag
  real(SP), dimension(:), pointer :: p_Fdata
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP) :: ft
  real(DP) :: dt
  integer :: i,it,j,j1,j2,k

    ! First of all, ensure that the matrix is of format 9
    if(rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) then
      call output_line('Matrix must be of format 9!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'lsyssc_sortMatrixCols')
      call sys_halt()
    end if

    ! Fetch the structure arrays of the matrix
    call storage_getbase_int(rmatrix%h_Kld,  p_IrowPtr)
    call storage_getbase_int(rmatrix%h_Kcol, p_IcolIdx)
    call storage_getbase_int(rmatrix%h_Kdiagonal, p_Idiag)

    ! Note: The sorting algorithm used below is linear insertion sort.
    !       Since we assume the number of non-zeroes per row to be
    !       "small", there shouldn't be any runtime issues here.

    ! Check whether the data array is allocated
    if(rmatrix%h_Da .eq. ST_NOHANDLE) then

      ! matrix data not present
      do i = 1, rmatrix%NEQ
        j1 = p_IrowPtr(i)
        j2 = p_IrowPtr(i+1)-1
        do j = j1, j2
          it = p_IcolIdx(j)
          do k = j-1, j1, -1
            if(p_IcolIdx(k) .le. it) exit
            p_IcolIdx(k+1) = p_IcolIdx(k)
          end do ! k
          p_IcolIdx(k+1) = it
        end do ! j
      end do ! i

    else

      select case (rmatrix%cdataType)
      case (ST_SINGLE)
        ! fetch the data array
        call storage_getbase_single(rmatrix%h_Da, p_Fdata)
        do i = 1, rmatrix%NEQ
          j1 = p_IrowPtr(i)
          j2 = p_IrowPtr(i+1)-1
          do j = j1, j2
            it = p_IcolIdx(j)
            ft = p_Fdata(j)
            do k = j-1, j1, -1
              if(p_IcolIdx(k) .le. it) exit
              p_IcolIdx(k+1) = p_IcolIdx(k)
              p_Fdata(k+1) = p_Fdata(k)
            end do ! k
            p_IcolIdx(k+1) = it
            p_Fdata(k+1) = ft
          end do ! j
        end do ! i

      case (ST_DOUBLE)
        ! fetch the data array
        call storage_getbase_double(rmatrix%h_Da, p_Ddata)
        do i = 1, rmatrix%NEQ
          j1 = p_IrowPtr(i)
          j2 = p_IrowPtr(i+1)-1
          do j = j1, j2
            it = p_IcolIdx(j)
            dt = p_Ddata(j)
            do k = j-1, j1, -1
              if(p_IcolIdx(k) .le. it) exit
              p_IcolIdx(k+1) = p_IcolIdx(k)
              p_Ddata(k+1) = p_Ddata(k)
            end do ! k
            p_IcolIdx(k+1) = it
            p_Ddata(k+1) = dt
          end do ! j
        end do ! i

      case default
        ! Datatype not supported
        call output_line('Only single and double precision matrices are supported!', &
            OU_CLASS_ERROR, OU_MODE_STD, 'lsyssc_sortMatrixCols')
        call sys_halt()
      end select

    end if

    ! Finally, recalculate the diagonal pointer array.
    call lsyssc_rebuildKdiagonal(p_IcolIdx, p_IrowPtr, p_Idiag, rmatrix%NEQ)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsyssc_createRowCMatrix(rmatrix,rdestMatrix)

!<description>
    ! Creates a 'differential' matrix. A 'differential' matrix contains only
    ! the nonzero rows of the original matrix and is therefore faster to handle.
!</description>

!<input>
    ! Source matrix to be compressed.
    type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<output>
    ! Destination matrix that receives the compressed data.
    ! If the matrix does not contain any nonzero data, an empty matrix
    ! is returned.
    type(t_matrixScalar), intent(out) :: rdestMatrix
!</output>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Kld,p_Kld2,p_Kcol,p_Kcol2
    integer, dimension(:), pointer :: p_KrowIdx
    integer :: nnzrows,irow,icol,ioffset,na,inewrow
    integer :: h_KrowIdx,h_Kld2,h_Kcol2,h_Da2
    real(DP) :: dtmp
    real(DP), dimension(:), pointer :: p_Da,p_Da2

    if (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) then
      call output_line('Only format 9 supported!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'lsyssc_createDiffMatrix')
      call sys_halt()
    end if

    if (rmatrix%cdataType .ne. ST_DOUBLE) then
      call output_line('Only double precision supported!', &
          OU_CLASS_ERROR, OU_MODE_STD, 'lsyssc_createDiffMatrix')
      call sys_halt()
    end if

    ! Loop through all rows and calculate the number of nonzero rows.
    ! That is the most simple way to avoid too much reallocation.
    ! Also compute the new size of the matrix.
    nnzrows = 0
    na = 0

    call lsyssc_getbase_Kld (rmatrix,p_Kld)
    call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_double (rmatrix,p_Da)

    do irow = 1,rmatrix%neq
      ! Calculate the norm of the to to check if it nonzero.
      ! If yes, remember that row.
      dtmp = 0.0_DP
      do icol = p_Kld(irow), p_Kld(irow+1)-1
        dtmp = dtmp + p_Da(icol)*p_Da(icol)
      end do

      if (dtmp .ne. 0.0_DP) then
        ! Remember that row
        nnzrows = nnzrows+1
        na = na + p_Kld(irow+1) - p_Kld(irow)
      end if
    end do

    ! Cancel if the matrix is empty.
    if (nnzrows .eq. 0) return

    ! Allocate Kld, Kcol, KrowIdx.
    call storage_new ('lsyssc_createDiffMatrix', 'Kld2', nnzrows+1, &
                      ST_INT, h_Kld2, ST_NEWBLOCK_NOINIT)
    call storage_new ('lsyssc_createDiffMatrix', 'KrowPtr', nnzrows+rmatrix%NEQ, &
                      ST_INT, h_KrowIdx, ST_NEWBLOCK_ZERO)
    call storage_new ('lsyssc_createDiffMatrix', 'Kcol2', na, &
                      ST_INT, h_Kcol2, ST_NEWBLOCK_NOINIT)
    call storage_new ('lsyssc_createDiffMatrix', 'Da2', na, &
                      ST_DOUBLE, h_Da2, ST_NEWBLOCK_NOINIT)

    ! Get the arrays
    call storage_getbase_int (h_Kld2,p_Kld2)
    call storage_getbase_int (h_Kcol2,p_Kcol2)
    call storage_getbase_int (h_KrowIdx,p_KrowIdx)
    call storage_getbase_double (h_Da2,p_Da2)

    ! Loop through the rows of the original matrix, collect the data.
    inewrow = 0
    p_Kld2(1) = 1
    do irow = 1,rmatrix%neq
      ! Calculate the norm of the to to check if it nonzero.
      ! If yes, remember that row.
      dtmp = 0.0_DP
      do icol = p_Kld(irow), p_Kld(irow+1)-1
        dtmp = dtmp + p_Da(icol)*p_Da(icol)
      end do

      if (dtmp .ne. 0.0_DP) then
        ! Remember that row
        inewrow = inewrow+1
        p_KrowIdx(inewrow) = irow

        ! Fill up KrowIdx(nnzrows+1:) with the 'inverse permutation'.
        p_KrowIdx(nnzrows+irow) = inewrow

        ! Sum up the row length of the nonzero rows.
        p_Kld2(inewrow+1) = p_Kld2(inewrow) + (p_Kld(irow+1) - p_Kld(irow))

        ! Copy the matrix data.
        ioffset = p_Kld(irow) - p_Kld2(inewrow)
        do icol = p_Kld2(inewrow), p_Kld2(inewrow+1)-1
          p_Da2(icol) = p_Da(ioffset+icol)
          p_Kcol2(icol) = p_Kcol(ioffset+icol)
        end do

      end if
    end do

    ! Create the destination matrix.
    call lsyssc_createEmptyMatrixStub (rdestMatrix,LSYSSC_MATRIX9ROWC,&
        rmatrix%neq,rmatrix%ncols)

    ! Copy basic properties
    rdestMatrix%bcolumnsSorted          = rmatrix%bcolumnsSorted
    rdestMatrix%browsSorted             = rmatrix%browsSorted
    rdestMatrix%p_rsortStrategyColumns => rmatrix%p_rsortStrategyColumns
    rdestMatrix%p_rsortStrategyRows => rmatrix%p_rsortStrategyRows
    rdestMatrix%p_rspatialDiscrTrial    => rmatrix%p_rspatialDiscrTrial
    rdestMatrix%p_rspatialDiscrTest     => rmatrix%p_rspatialDiscrTest
    rdestMatrix%bidenticalTrialAndTest  = rmatrix%bidenticalTrialAndTest

    ! insert the handles.
    rdestMatrix%h_Da = h_Da2
    rdestMatrix%h_Kcol = h_Kcol2
    rdestMatrix%h_Kld = h_Kld2
    rdestMatrix%h_KrowIdx = h_KrowIdx
    rdestMatrix%nnzrows = nnzrows
    rdestMatrix%na = na

  end subroutine

  !*****************************************************************************

!<subroutine>

  subroutine lsyssc_genEdgeList(rmatrix, h_IedgeList, ccontentType,&
      bisSymmetric, bignoreDiagonal, nedge, Iselection)

!<description>
    ! This subroutine generates a list of edges (i,j) based on all or
    ! a subset of non-zero matrix entries ij by transversing the
    ! matrix. Depending on the flags bisSymmetric and bignoreDiagonal
    ! it is assumed that the matrix is symmetric so that only the
    ! upper-right triangular part is traversed and/or the diagonal
    ! "edges" (i,i) are ignored or included.
!</description>

!<input>
    ! Source matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! Specifier which defines the content to be stored at each edge.
    ! Must be one of the LSYSSC_EDGELIST_xxx constants
    integer, intent(in) :: ccontentType

    ! Flag which defines whether the matrix is structurally symmetric,
    ! that is, an edge (i,j) exists if an only if an edge (j,i) exists
    logical, intent(in) :: bisSymmetric

    ! Flag which defines whether the diagonal edges (i,i) should be
    ! ignored or not.
    logical, intent(in) :: bignoreDiagonal

    ! OPTIONAL: a list of degress of freedoms to which the edge
    ! structure should be restricted.
    integer, dimension(:), intent(in), optional :: Iselection
!</input>

!<inputoutput>
    ! Handle for the edge list
    integer, intent(inout) :: h_IedgeList

    ! On input: if nedge > 0 then this value is adopted as the number
    !           of edges precomputed externally; otherwise the
    !           number of edges is computed internally and returned
    ! On output: the number of edges
    integer, intent(inout) :: nedge
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IedgeList
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Ksep,p_Kdiagonal
    integer, dimension(2) :: Isize
    integer :: h_Ksep,na,ncols,ncontent,neq

    ! Determine the number of edges
    if (nedge .le. 0)&
        call lsyssc_calcDimsFromMatrix(rmatrix, na, neq, ncols, nedge, Iselection)

    ! Determine the type of information stored at each edge
    select case (ccontentType)
    case (LSYSSC_EDGELIST_NODESONLY)
      ncontent = 2
    case (LSYSSC_EDGELIST_POSONLY)
      ncontent = 4
    case (LSYSSC_EDGELIST_NODESANDPOS)
      ncontent = 6
    case default
      call output_line('Unsupported content specifier!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_genEdgeList')
      call sys_halt()
    end select

    ! (Re-)allocate memory
    if (h_IedgeList .eq. ST_NOHANDLE) then
      Isize = (/ncontent, nedge/)
      call storage_new('lsyssc_genEdgeList', 'IedgeList',&
          Isize, ST_INT, h_IedgeList, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_IedgeList, Isize)
      if ((Isize(1) .ne. ncontent) .or. (Isize(2) .ne. nedge)) then
        call storage_free(h_IedgeList)
        Isize = (/ncontent, nedge/)
        call storage_new('lsyssc_genEdgeList', 'IedgeList',&
            Isize, ST_INT, h_IedgeList, ST_NEWBLOCK_NOINIT)
      end if
    end if

    ! Set pointer
    call storage_getbase_int2d(h_IedgeList, p_IedgeList)

    ! Generate edge data structure
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX1)
      if (bisSymmetric) then
        if (present(Iselection)) then
          call genEdgesMat1SymmSel(rmatrix%NEQ, ccontentType,&
              bignoreDiagonal, Iselection, p_IedgeList)
        else
          call genEdgesMat1Symm(rmatrix%NEQ, ccontentType,&
              bignoreDiagonal, p_IedgeList)
        end if
      else
        if (present(Iselection)) then
          call genEdgesMat1Sel(rmatrix%NEQ, rmatrix%NCOLS, ccontentType,&
              bignoreDiagonal, Iselection, p_IedgeList)
        else
          call genEdgesMat1(rmatrix%NEQ, rmatrix%NCOLS, ccontentType,&
              bignoreDiagonal, p_IedgeList)
        end if
      end if

    case (LSYSSC_MATRIX7, LSYSSC_MATRIX7INTL)
      if (bisSymmetric) then
        ! Set pointers
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rmatrix%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rmatrix%NEQ+1)

        ! Generate edge structure
        if (present(Iselection)) then
          call genEdgesMat7SymmSel(p_Kld, p_Kcol, p_Ksep, ccontentType,&
              bignoreDiagonal, Iselection, p_IedgeList)
        else
          call genEdgesMat7Symm(p_Kld, p_Kcol, p_Ksep, ccontentType,&
              bignoreDiagonal, p_IedgeList)
        end if

        ! Release diagonal separator
        call storage_free(h_Ksep)
      else
        if (present(Iselection)) then
          call genEdgesMat7Sel(p_Kld, p_Kcol, ccontentType,&
              bignoreDiagonal, max(rmatrix%NEQ,rmatrix%NCOLS),&
              Iselection, p_IedgeList)
        else
          call genEdgesMat7(p_Kld, p_Kcol, ccontentType,&
              bignoreDiagonal, p_IedgeList)
        end if
      end if

    case (LSYSSC_MATRIX9, LSYSSC_MATRIX9INTL)
      if (bisSymmetric) then
        ! Set pointers
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
        call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

        ! Create diagonal separator
        h_Ksep = ST_NOHANDLE
        call storage_copy(rmatrix%h_Kld, h_Ksep)
        call storage_getbase_int(h_Ksep, p_Ksep, rmatrix%NEQ+1)

        ! Generate edge structure
        if (present(Iselection)) then
          call genEdgesMat9SymmSel(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              ccontentType, bignoreDiagonal, Iselection, p_IedgeList)
        else
          call genEdgesMat9Symm(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep,&
              ccontentType, bignoreDiagonal, p_IedgeList)
        end if

        ! Release diagonal separator
        call storage_free(h_Ksep)
      else
        if (present(Iselection)) then
          call genEdgesMat9Sel(p_Kld, p_Kcol, ccontentType,&
              bignoreDiagonal, max(rmatrix%NEQ,rmatrix%NCOLS),&
              Iselection, p_IedgeList)
        else
          call genEdgesMat9(p_Kld, p_Kcol, ccontentType,&
              bignoreDiagonal, p_IedgeList)
        end if
      end if

    end select

  contains

    ! Here, the working routines follow

    !**************************************************************
    ! Generate edge data structure for symmetric dense matrices
    ! stored in matrix format 1

    subroutine genEdgesMat1Symm(neq, ccontentType, bignoreDiagonal, IedgeList)

      integer, intent(in) :: neq,ccontentType
      logical, intent(in) :: bignoreDiagonal
      integer, dimension(:,:), intent(inout) :: IedgeList

      ! local variables
      integer :: i,j,iedge,ioffset

      ! Initialise edge counter
      iedge = 0

      ! Include diagonal entries?
      ioffset = merge(1,0, bignoreDiagonal)

      select case (ccontentType)

      case (LSYSSC_EDGELIST_NODESONLY)

        ! Loop over all equations
        do i = 1, neq

          ! Loop over all (off-diagonal) matrix entries IJ which are
          ! adjacent to node J such that I <= J or I < J, that  is,
          ! the upper right triangular part of the matrix is traversed
          do j = i+ioffset, neq

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j
          end do
        end do

      case (LSYSSC_EDGELIST_POSONLY)

        ! Loop over all equations
        do i = 1, neq

          ! Loop over all (off-diagonal) matrix entries IJ which are
          ! adjacent to node J such that I <= J or I < J, that  is,
          ! the upper right triangular part of the matrix is traversed
          do j = i+ioffset, neq

            ! Increase edge counter
            iedge = iedge+1

            ! Set matrix positions ij, ji, ii, and jj
            IedgeList(1,iedge) = neq*(i-1)+j
            IedgeList(2,iedge) = neq*(j-1)+i
            IedgeList(3,iedge) = neq*(i-1)+i
            IedgeList(4,iedge) = neq*(j-1)+j
          end do
        end do

      case (LSYSSC_EDGELIST_NODESANDPOS)

        ! Loop over all equations
        do i = 1, neq

          ! Loop over all (off-diagonal) matrix entries IJ which are
          ! adjacent to node J such that I <= J or I < J, that  is,
          ! the upper right triangular part of the matrix is traversed
          do j = i+ioffset, neq

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j

            ! Set matrix positions ij, ji, ii, and jj
            IedgeList(3,iedge) = neq*(i-1)+j
            IedgeList(4,iedge) = neq*(j-1)+i
            IedgeList(5,iedge) = neq*(i-1)+i
            IedgeList(6,iedge) = neq*(j-1)+j
          end do
        end do

      end select

    end subroutine genEdgesMat1Symm

    !**************************************************************
    ! Generate edge data structure for symmetric dense matrices stored
    ! in matrix format 1 restricted to a list of selected entries

    subroutine genEdgesMat1SymmSel(neq, ccontentType, bignoreDiagonal,&
        Iselection, IedgeList)

      integer, dimension(:), intent(in) :: Iselection
      integer, intent(in) :: neq,ccontentType
      logical, intent(in) :: bignoreDiagonal
      integer, dimension(:,:), intent(inout) :: IedgeList

      ! local variables
      logical, dimension(:), allocatable :: BisActive
      integer :: i,j,iedge,ioffset

      ! Generate set of active entries
      allocate(BisActive(neq)); BisActive=.false.
      do i = 1, size(Iselection)
        BisActive(Iselection(i))=.true.
      end do

      ! Initialise edge counter
      iedge = 0

      ! Include diagonal entries?
      ioffset = merge(1,0, bignoreDiagonal)

      select case (ccontentType)

      case (LSYSSC_EDGELIST_NODESONLY)

        ! Loop over all equations
        do i = 1, neq

          ! Check if this row belongs to an active DOF
          if (.not.BisActive(i)) cycle

          ! Loop over all (off-diagonal) matrix entries IJ which are
          ! adjacent to node J such that I <= J or I < J, that  is,
          ! the upper right triangular part of the matrix is traversed
          do j = i+ioffset, neq

            ! Check if this edge connects two DOFs which are marked active
            if (.not.BisActive(j)) cycle

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j
          end do
        end do

      case (LSYSSC_EDGELIST_POSONLY)

        ! Loop over all equations
        do i = 1, neq

          ! Check if this row belongs to an active DOF
          if (.not.BisActive(i)) cycle

          ! Loop over all (off-diagonal) matrix entries IJ which are
          ! adjacent to node J such that I <= J or I < J, that  is,
          ! the upper right triangular part of the matrix is traversed
          do j = i+ioffset, neq

            ! Check if this edge connects two DOFs which are marked active
            if (.not.BisActive(j)) cycle

            ! Increase edge counter
            iedge = iedge+1

            ! Set matrix positions ij, ji, ii, and jj
            IedgeList(1,iedge) = neq*(i-1)+j
            IedgeList(2,iedge) = neq*(j-1)+i
            IedgeList(3,iedge) = neq*(i-1)+i
            IedgeList(4,iedge) = neq*(j-1)+j
          end do
        end do

      case (LSYSSC_EDGELIST_NODESANDPOS)

        ! Loop over all equations
        do i = 1, neq

          ! Check if this row belongs to an active DOF
          if (.not.BisActive(i)) cycle

          ! Loop over all (off-diagonal) matrix entries IJ which are
          ! adjacent to node J such that I <= J or I < J, that  is,
          ! the upper right triangular part of the matrix is traversed
          do j = i+ioffset, neq

            ! Check if this edge connects two DOFs which are marked active
            if (.not.BisActive(j)) cycle

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j

            ! Set matrix positions ij, ji, ii, and jj
            IedgeList(3,iedge) = neq*(i-1)+j
            IedgeList(4,iedge) = neq*(j-1)+i
            IedgeList(5,iedge) = neq*(i-1)+i
            IedgeList(6,iedge) = neq*(j-1)+j
          end do
        end do

      end select

      ! Deallocate temporal memory
      deallocate(BisActive)

    end subroutine genEdgesMat1SymmSel

    !**************************************************************
    ! Generate edge data structure for arbitrary dense matrices
    ! stored in matrix format 1

    subroutine genEdgesMat1(neq, ncols, ccontentType, bignoreDiagonal, IedgeList)

      integer, intent(in) :: neq,ncols,ccontentType
      logical, intent(in) :: bignoreDiagonal
      integer, dimension(:,:), intent(inout) :: IedgeList

      ! local variables
      integer :: i,j,iedge

      ! Initialise edge counter
      iedge = 0

      select case (ccontentType)

      case (LSYSSC_EDGELIST_NODESONLY)

        ! Loop over all equations
        do i = 1, neq

          ! Loop over all (off-diagonal) matrix entries IJ which are
          ! adjacent to node J
          do j = 1, ncols

            ! Ignore the diagonal entries?
            if (bignoreDiagonal .and. i.eq.j) cycle

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j
          end do
        end do

      case (LSYSSC_EDGELIST_POSONLY)

        ! Loop over all equations
        do i = 1, neq

          ! Loop over all (off-diagonal) matrix entries IJ which are
          ! adjacent to node J
          do j = 1, ncols

            ! Ignore the diagonal entries?
            if (bignoreDiagonal .and. i.eq.j) cycle

            ! Increase edge counter
            iedge = iedge+1

            ! Set matrix positions ij, ji, ii, and jj
            IedgeList(1,iedge) = neq*(i-1)+j
            IedgeList(2,iedge) = neq*(j-1)+i
            IedgeList(3,iedge) = neq*(i-1)+i
            IedgeList(4,iedge) = neq*(j-1)+j
          end do
        end do

      case (LSYSSC_EDGELIST_NODESANDPOS)

        ! Loop over all equations
        do i = 1, neq

          ! Loop over all (off-diagonal) matrix entries IJ which are
          ! adjacent to node J
          do j = 1, ncols

            ! Ignore the diagonal entries?
            if (bignoreDiagonal .and. i.eq.j) cycle

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j

            ! Set matrix positions ij, ji, ii, and jj
            IedgeList(3,iedge) = neq*(i-1)+j
            IedgeList(4,iedge) = neq*(j-1)+i
            IedgeList(5,iedge) = neq*(i-1)+i
            IedgeList(6,iedge) = neq*(j-1)+j
          end do
        end do

      end select

    end subroutine genEdgesMat1

    !**************************************************************
    ! Generate edge data structure for arbitrary dense matrices stored
    ! in matrix format 1 restricted to a list of selected entries

    subroutine genEdgesMat1Sel(neq, ncols, ccontentType, bignoreDiagonal,&
        Iselection, IedgeList)

      integer, dimension(:), intent(in) :: Iselection
      integer, intent(in) :: neq,ncols,ccontentType
      logical, intent(in) :: bignoreDiagonal
      integer, dimension(:,:), intent(inout) :: IedgeList

      ! local variables
      logical, dimension(:), allocatable :: BisActive
      integer :: i,j,iedge

      ! Generate set of active entries
      allocate(BisActive(max(neq,ncols))); BisActive=.false.
      do i = 1, size(Iselection)
        BisActive(Iselection(i))=.true.
      end do

      ! Initialise edge counter
      iedge = 0

      select case (ccontentType)

      case (LSYSSC_EDGELIST_NODESONLY)

        ! Loop over all equations
        do i = 1, neq

          ! Check if this row belongs to an active DOF
          if (.not.BisActive(i)) cycle

          ! Loop over all (off-diagonal) matrix entries IJ which are
          ! adjacent to node J
          do j = 1, ncols

            ! Check if this edge connects two DOFs which are marked active
            if (.not.BisActive(j)) cycle

            ! Ignore the diagonal entries?
            if (bignoreDiagonal .and. i.eq.j) cycle

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j
          end do
        end do

      case (LSYSSC_EDGELIST_POSONLY)

        ! Loop over all equations
        do i = 1, neq

          ! Check if this row belongs to an active DOF
          if (.not.BisActive(i)) cycle

          ! Loop over all (off-diagonal) matrix entries IJ which are
          ! adjacent to node J
          do j = 1, ncols

            ! Check if this edge connects two DOFs which are marked active
            if (.not.BisActive(j)) cycle

            ! Ignore the diagonal entries?
            if (bignoreDiagonal .and. i.eq.j) cycle

            ! Increase edge counter
            iedge = iedge+1

            ! Set matrix positions ij, ji, ii, and jj
            IedgeList(1,iedge) = neq*(i-1)+j
            IedgeList(2,iedge) = neq*(j-1)+i
            IedgeList(3,iedge) = neq*(i-1)+i
            IedgeList(4,iedge) = neq*(j-1)+j
          end do
        end do

      case (LSYSSC_EDGELIST_NODESANDPOS)

        ! Loop over all equations
        do i = 1, neq

          ! Check if this row belongs to an active DOF
          if (.not.BisActive(i)) cycle

          ! Loop over all (off-diagonal) matrix entries IJ which are
          ! adjacent to node J
          do j = 1, ncols

            ! Check if this edge connects two DOFs which are marked active
            if (.not.BisActive(j)) cycle

            ! Ignore the diagonal entries?
            if (bignoreDiagonal .and. i.eq.j) cycle

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j

            ! Set matrix positions ij, ji, ii, and jj
            IedgeList(3,iedge) = neq*(i-1)+j
            IedgeList(4,iedge) = neq*(j-1)+i
            IedgeList(5,iedge) = neq*(i-1)+i
            IedgeList(6,iedge) = neq*(j-1)+j
          end do
        end do

      end select

      ! Deallocate temporal memory
      deallocate(BisActive)

    end subroutine genEdgesMat1Sel

    !**************************************************************
    ! Generate edge data structure for symmetric matrix in format 7

    subroutine genEdgesMat7Symm(Kld, Kcol, Ksep, ccontentType,&
        bignoreDiagonal, IedgeList)

      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: ccontentType
      logical, intent(in) :: bignoreDiagonal
      integer, dimension(:,:), intent(inout) :: IedgeList
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      integer :: i,j,ij,ji,iedge

      ! Initialise edge counter
      iedge = 0

      select case (ccontentType)

      case (LSYSSC_EDGELIST_NODESONLY)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          if (.not.bignoreDiagonal) then
            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and i
            IedgeList(1:2,iedge) = i
          end if

          ! Loop over all off-diagonal matrix entries IJ which are
          ! adjacent to node J such that I < J, that is, the upper
          ! right triangular part of the matrix is traversed
          do ij = Ksep(i)+1, Kld(i+1)-1

            ! Get node number J, the corresponding matrix positions JI,
            ! and let the separator point to the next entry
            j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j
          end do
        end do

      case (LSYSSC_EDGELIST_POSONLY)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          if (.not.bignoreDiagonal) then
            ! Increase edge counter
            iedge = iedge+1

            ! Set matrix position ii
            IedgeList(1:4,iedge) = Kld(i)
          end if

          ! Loop over all off-diagonal matrix entries IJ which are
          ! adjacent to node J such that I < J, that is, the upper
          ! right triangular part of the matrix is traversed
          do ij = Ksep(i)+1, Kld(i+1)-1

            ! Get node number J, the corresponding matrix positions JI,
            ! and let the separator point to the next entry
            j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

            ! Increase edge counter
            iedge = iedge+1

            ! Set matrix positions ij, ji, ii, and ,jj
            IedgeList(1,iedge) = ij
            IedgeList(2,iedge) = ji
            IedgeList(3,iedge) = Kld(i)
            IedgeList(4,iedge) = Kld(j)
          end do
        end do

      case (LSYSSC_EDGELIST_NODESANDPOS)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          if (.not.bignoreDiagonal) then
            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and i
            IedgeList(1:2,iedge) = i

            ! Set matrix position ii
            IedgeList(3:6,iedge) = Kld(i)
          end if

          ! Loop over all off-diagonal matrix entries IJ which are
          ! adjacent to node J such that I < J, that is, the upper
          ! right triangular part of the matrix is traversed
          do ij = Ksep(i)+1, Kld(i+1)-1

            ! Get node number J, the corresponding matrix positions JI,
            ! and let the separator point to the next entry
            j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j

            ! Set matrix positions ij, ji, ii, and jj
            IedgeList(3,iedge) = ij
            IedgeList(4,iedge) = ji
            IedgeList(5,iedge) = Kld(i)
            IedgeList(6,iedge) = Kld(j)
          end do
        end do

      end select

    end subroutine genEdgesMat7Symm

    !**************************************************************
    ! Generate edge data structure for symmetric matrix in format 7
    ! restricted to a list of selected entries

    subroutine genEdgesMat7SymmSel(Kld, Kcol, Ksep, ccontentType,&
        bignoreDiagonal, Iselection, IedgeList)

      integer, dimension(:), intent(in) :: Kld,Kcol,Iselection
      integer, intent(in) :: ccontentType
      logical, intent(in) :: bignoreDiagonal
      integer, dimension(:,:), intent(inout) :: IedgeList
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      logical, dimension(:), allocatable :: BisActive
      integer :: i,j,ij,ji,iedge

      ! Generate set of active entries
      allocate(BisActive(size(Kld)-1)); BisActive=.false.
      do i = 1, size(Iselection)
        BisActive(Iselection(i))=.true.
      end do

      ! Initialise edge counter
      iedge = 0

      select case (ccontentType)

      case (LSYSSC_EDGELIST_NODESONLY)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          ! Check if this row belongs to an active DOF
          if (BisActive(i)) then

            if (.not.bignoreDiagonal) then
              ! Increase edge counter
              iedge = iedge+1

              ! Set node numbers i and i
              IedgeList(1:2,iedge) = i
            end if

            ! Loop over all off-diagonal matrix entries IJ which are
            ! adjacent to node J such that I < J, that is, the upper
            ! right triangular part of the matrix is traversed
            do ij = Ksep(i)+1, Kld(i+1)-1

              ! Get node number J, the corresponding matrix positions JI,
              ! and let the separator point to the next entry
              j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

              ! Check if this edge connects two DOFs which are marked active
              if (.not.BisActive(j)) cycle

              ! Increase edge counter
              iedge = iedge+1

              ! Set node numbers i and j
              IedgeList(1,iedge) = i
              IedgeList(2,iedge) = j
            end do

          else

            ! Row is not active, so just update Ksep
            do ij = Ksep(i)+1, Kld(i+1)-1
              j = Kcol(ij); Ksep(j) = Ksep(j)+1
            end do

          end if
        end do

      case (LSYSSC_EDGELIST_POSONLY)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          ! Check if this row belongs to an active DOF
          if (BisActive(i)) then

            if (.not.bignoreDiagonal) then
              ! Increase edge counter
              iedge = iedge+1

              ! Set matrix position ii
              IedgeList(1:4,iedge) = Kld(i)
            end if

            ! Loop over all off-diagonal matrix entries IJ which are
            ! adjacent to node J such that I < J, that is, the upper
            ! right triangular part of the matrix is traversed
            do ij = Ksep(i)+1, Kld(i+1)-1

              ! Get node number J, the corresponding matrix positions JI,
              ! and let the separator point to the next entry
              j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

              ! Check if this edge connects two DOFs which are marked active
              if (.not.BisActive(j)) cycle

              ! Increase edge counter
              iedge = iedge+1

              ! Set matrix positions ij, ji, ii, and jj
              IedgeList(1,iedge) = ij
              IedgeList(2,iedge) = ji
              IedgeList(3,iedge) = Kld(i)
              IedgeList(4,iedge) = Kld(j)
            end do

          else

            ! Row is not active, so just update Ksep
            do ij = Ksep(i)+1, Kld(i+1)-1
              j = Kcol(ij); Ksep(j) = Ksep(j)+1
            end do

          end if
        end do

      case (LSYSSC_EDGELIST_NODESANDPOS)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          ! Check if this row belongs to an active DOF
          if (BisActive(i)) then

            if (.not.bignoreDiagonal) then
              ! Increase edge counter
              iedge = iedge+1

              ! Set node numbers i and i
              IedgeList(1:2,iedge) = i

              ! Set matrix position ii
              IedgeList(3:6,iedge) = Kld(i)
            end if

            ! Loop over all off-diagonal matrix entries IJ which are
            ! adjacent to node J such that I < J, that is, the upper
            ! right triangular part of the matrix is traversed
            do ij = Ksep(i)+1, Kld(i+1)-1

              ! Get node number J, the corresponding matrix positions JI,
              ! and let the separator point to the next entry
              j = Kcol(ij); Ksep(j) = Ksep(j)+1; ji = Ksep(j)

              ! Check if this edge connects two DOFs which are marked active
              if (.not.BisActive(j)) cycle

              ! Increase edge counter
              iedge = iedge+1

              ! Set node numbers i and j
              IedgeList(1,iedge) = i
              IedgeList(2,iedge) = j

              ! Set matrix positions ij, ji, ii, and jj
              IedgeList(3,iedge) = ij
              IedgeList(4,iedge) = ji
              IedgeList(5,iedge) = Kld(i)
              IedgeList(6,iedge) = Kld(j)
            end do

          else

            ! Row is not active, so just update Ksep
            do ij = Ksep(i)+1, Kld(i+1)-1
              j = Kcol(ij); Ksep(j) = Ksep(j)+1
            end do

          end if
        end do

      end select

      ! Deallocate temporal memory
      deallocate(BisActive)

    end subroutine genEdgesMat7SymmSel

    !**************************************************************
    ! Generate edge data structure for arbitrary matrix in format 7

    subroutine genEdgesMat7(Kld, Kcol, ccontentType,&
        bignoreDiagonal, IedgeList)

      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: ccontentType
      logical, intent(in) :: bignoreDiagonal
      integer, dimension(:,:), intent(inout) :: IedgeList

      ! local variables
      integer :: i,j,ij,iedge,ioffset

      ! Initialise edge counter
      iedge = 0

      ! Include diagonal entries?
      ioffset = merge(1,0, bignoreDiagonal)

      select case (ccontentType)

      case (LSYSSC_EDGELIST_NODESONLY)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          ! Loop over all (off-diagonal) matrix entries IJ
          do ij = Kld(i)+ioffset, Kld(i+1)-1

            ! Get node number j
            j = Kcol(ij)

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j
          end do
        end do

      case (LSYSSC_EDGELIST_POSONLY,LSYSSC_EDGELIST_NODESANDPOS)
        call output_line('Positions are not available for unsymmetric matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'genEdgesMat7')
        call sys_halt()
      end select

    end subroutine genEdgesMat7

    !**************************************************************
    ! Generate edge data structure for arbitrary matrix in format 7
    ! restricted to a list of selected entries

    subroutine genEdgesMat7Sel(Kld, Kcol, ccontentType,&
        bignoreDiagonal, ndof, Iselection, IedgeList)

      integer, dimension(:), intent(in) :: Kld,Kcol,Iselection
      integer, intent(in) :: ndof,ccontentType
      logical, intent(in) :: bignoreDiagonal
      integer, dimension(:,:), intent(inout) :: IedgeList

      ! local variables
      logical, dimension(:), allocatable :: BisActive
      integer :: i,j,ij,iedge,ioffset

      ! Generate set of active entries
      allocate(BisActive(ndof)); BisActive=.false.
      do i = 1, size(Iselection)
        BisActive(Iselection(i))=.true.
      end do

      ! Initialise edge counter
      iedge = 0

      ! Include diagonal entries?
      ioffset = merge(1,0, bignoreDiagonal)

      select case (ccontentType)

      case (LSYSSC_EDGELIST_NODESONLY)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          ! Check if this row belongs to an active DOF
          if (.not.BisActive(i)) cycle

          ! Loop over all (off-diagonal) matrix entries IJ
          do ij = Kld(i)+ioffset, Kld(i+1)-1

            ! Get node number j
            j = Kcol(ij)

            ! Check if this edge connects two DOFs which are marked active
            if (.not.BisActive(j)) cycle

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j
          end do
        end do

      case (LSYSSC_EDGELIST_POSONLY,LSYSSC_EDGELIST_NODESANDPOS)
        call output_line('Positions are not available for unsymmetric matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'genEdgesMat7Sel')
        call sys_halt()
      end select

      ! Deallocate temporal memory
      deallocate(BisActive)

    end subroutine genEdgesMat7Sel

    !**************************************************************
    ! Generate edge data structure for symmetric matrix in format 9

    subroutine genEdgesMat9Symm(Kld, Kcol, Kdiagonal, Ksep, ccontentType,&
        bignoreDiagonal, IedgeList)

      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal
      integer, intent(in) :: ccontentType
      logical, intent(in) :: bignoreDiagonal
      integer, dimension(:,:), intent(inout) :: IedgeList
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      integer :: i,j,ij,ji,iedge

      ! Initialise edge counter
      iedge = 0

      select case (ccontentType)

      case (LSYSSC_EDGELIST_NODESONLY)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          if (.not.bignoreDiagonal) then
            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and i
            IedgeList(1:2,iedge) = i
          end if

          ! Loop over all off-diagonal matrix entries IJ which are
          ! adjacent to node J such that I < J, that is, the upper
          ! right triangular part of the matrix is traversed
          do ij = Kdiagonal(i)+1, Kld(i+1)-1

            ! Get node number J, the corresponding matrix positions JI,
            ! and let the separator point to the next entry
            j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j
          end do
        end do

      case (LSYSSC_EDGELIST_POSONLY)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          if (.not.bignoreDiagonal) then
            ! Increase edge counter
            iedge = iedge+1

            ! Set matrix position ii
            IedgeList(1:4,iedge) = Kdiagonal(i)
          end if

          ! Loop over all off-diagonal matrix entries IJ which are
          ! adjacent to node J such that I < J. That is, explore the
          ! upper triangular matrix
          do ij = Kdiagonal(i)+1, Kld(i+1)-1

            ! Get node number J, the corresponding matrix positions JI,
            ! and let the separator point to the next entry
            j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

            ! Increase edge counter
            iedge = iedge+1

            ! Set matrix positions ij, ji, ii, and jj
            IedgeList(1,iedge) = ij
            IedgeList(2,iedge) = ji
            IedgeList(3,iedge) = Kdiagonal(i)
            IedgeList(4,iedge) = Kdiagonal(j)
          end do
        end do

      case (LSYSSC_EDGELIST_NODESANDPOS)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          if (.not.bignoreDiagonal) then
            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and i
            IedgeList(1:2,iedge) = i

            ! Set node numbers i and i
            IedgeList(3:6,iedge) = Kdiagonal(i)
          end if

          ! Loop over all off-diagonal matrix entries IJ which are
          ! adjacent to node J such that I < J. That is, explore the
          ! upper triangular matrix
          do ij = Kdiagonal(i)+1, Kld(i+1)-1

            ! Get node number J, the corresponding matrix positions JI,
            ! and let the separator point to the next entry
            j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j

            ! Set matrix positions ij and ji
            IedgeList(3,iedge) = ij
            IedgeList(4,iedge) = ji
            IedgeList(5,iedge) = Kdiagonal(i)
            IedgeList(6,iedge) = Kdiagonal(j)
          end do
        end do

      end select

    end subroutine genEdgesMat9Symm

    !**************************************************************
    ! Generate edge data structure for symmetric matrix in format 9
    ! restricted to a list of selected entries

    subroutine genEdgesMat9SymmSel(Kld, Kcol, Kdiagonal, Ksep, ccontentType,&
        bignoreDiagonal, Iselection, IedgeList)

      integer, dimension(:), intent(in) :: Kld,Kcol,Kdiagonal,Iselection
      integer, intent(in) :: ccontentType
      logical, intent(in) :: bignoreDiagonal
      integer, dimension(:,:), intent(inout) :: IedgeList
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      logical, dimension(:), allocatable :: BisActive
      integer :: i,j,ij,ji,iedge

      ! Generate set of active entries
      allocate(BisActive(size(Kld)-1)); BisActive=.false.
      do i = 1, size(Iselection)
        BisActive(Iselection(i))=.true.
      end do

      ! Initialise edge counter
      iedge = 0

      select case (ccontentType)

      case (LSYSSC_EDGELIST_NODESONLY)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          ! Check if this row belongs to an active DOF
          if (BisActive(i)) then

            if (.not.bignoreDiagonal) then
              ! Increase edge counter
              iedge = iedge+1

              ! Set node numbers i and i
              IedgeList(1:2,iedge) = i
            end if

            ! Loop over all off-diagonal matrix entries IJ which are
            ! adjacent to node J such that I < J, that is, the upper
            ! right triangular part of the matrix is traversed
            do ij = Kdiagonal(i)+1, Kld(i+1)-1

              ! Get node number J, the corresponding matrix positions JI,
              ! and let the separator point to the next entry
              j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

              ! Check if this edge connects two DOFs which are marked active
              if (.not.BisActive(j)) cycle

              ! Increase edge counter
              iedge = iedge+1

              ! Set node numbers i and j
              IedgeList(1,iedge) = i
              IedgeList(2,iedge) = j
            end do

          else

            ! Row is not active, so just update Ksep
            do ij = Ksep(i)+1, Kld(i+1)-1
              j = Kcol(ij); Ksep(j) = Ksep(j)+1
            end do

          end if
        end do

      case (LSYSSC_EDGELIST_POSONLY)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          ! Check if this row belongs to an active DOF
          if (BisActive(i)) then

            if (.not.bignoreDiagonal) then
              ! Increase edge counter
              iedge = iedge+1

              ! Set node numbers i and i
              IedgeList(1:4,iedge) = Kdiagonal(i)
            end if

            ! Loop over all off-diagonal matrix entries IJ which are
            ! adjacent to node J such that I < J. That is, explore the
            ! upper triangular matrix
            do ij = Kdiagonal(i)+1, Kld(i+1)-1

              ! Get node number J, the corresponding matrix positions JI,
              ! and let the separator point to the next entry
              j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

              ! Check if this edge connects two DOFs which are marked active
              if (.not.BisActive(j)) cycle

              ! Increase edge counter
              iedge = iedge+1

              ! Set matrix positions ij and ji
              IedgeList(1,iedge) = ij
              IedgeList(2,iedge) = ji
              IedgeList(3,iedge) = Kdiagonal(i)
              IedgeList(4,iedge) = Kdiagonal(j)
            end do

          else

            ! Row is not active, so just update Ksep
            do ij = Ksep(i)+1, Kld(i+1)-1
              j = Kcol(ij); Ksep(j) = Ksep(j)+1
            end do

          end if

        end do

      case (LSYSSC_EDGELIST_NODESANDPOS)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          ! Check if this row belongs to an active DOF
          if (BisActive(i)) then

            if (.not.bignoreDiagonal) then
              ! Increase edge counter
              iedge = iedge+1

              ! Set node numbers i and i
              IedgeList(1:2,iedge) = i

              ! Set node numbers i and i
              IedgeList(3:6,iedge) = Kdiagonal(i)
            end if

            ! Loop over all off-diagonal matrix entries IJ which are
            ! adjacent to node J such that I < J. That is, explore the
            ! upper triangular matrix
            do ij = Kdiagonal(i)+1, Kld(i+1)-1

              ! Get node number J, the corresponding matrix positions JI,
              ! and let the separator point to the next entry
              j = Kcol(ij); ji = Ksep(j); Ksep(j) = Ksep(j)+1

              ! Check if this edge connects two DOFs which are marked active
              if (.not.BisActive(j)) cycle

              ! Increase edge counter
              iedge = iedge+1

              ! Set node numbers i and j
              IedgeList(1,iedge) = i
              IedgeList(2,iedge) = j

              ! Set matrix positions ij and ji
              IedgeList(3,iedge) = ij
              IedgeList(4,iedge) = ji
              IedgeList(5,iedge) = Kdiagonal(i)
              IedgeList(6,iedge) = Kdiagonal(j)
            end do

          else

            ! Row is not active, so just update Ksep
            do ij = Ksep(i)+1, Kld(i+1)-1
              j = Kcol(ij); Ksep(j) = Ksep(j)+1
            end do

          end if
        end do

      end select

      ! Deallocate temporal memory
      deallocate(BisActive)

    end subroutine genEdgesMat9SymmSel

    !**************************************************************
    ! Generate edge data structure for arbitrary matrix in format 9

    subroutine genEdgesMat9(Kld, Kcol, ccontentType,&
        bignoreDiagonal, IedgeList)

      integer, dimension(:), intent(in) :: Kld,Kcol
      integer, intent(in) :: ccontentType
      logical, intent(in) :: bignoreDiagonal
      integer, dimension(:,:), intent(inout) :: IedgeList

      ! local variables
      integer :: i,j,ij,iedge

      ! Initialise edge counter
      iedge = 0

      select case (ccontentType)

      case (LSYSSC_EDGELIST_NODESONLY)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          ! Loop over all (off-diagonal) matrix entries IJ
          do ij = Kld(i), Kld(i+1)-1

            ! Get node number j
            j = Kcol(ij)

            ! Ignore the diagonal entries?
            if (bignoreDiagonal .and. i.eq.j) cycle

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j
          end do
        end do

      case (LSYSSC_EDGELIST_POSONLY,LSYSSC_EDGELIST_NODESANDPOS)
        call output_line('Positions are not available for unsymmetric matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'genEdgesMat9')
        call sys_halt()
      end select

    end subroutine genEdgesMat9

    !**************************************************************
    ! Generate edge data structure for arbitrary matrix in format 9
    ! restricted to a list of selected entries

    subroutine genEdgesMat9Sel(Kld, Kcol, ccontentType,&
        bignoreDiagonal, ndof, Iselection, IedgeList)

      integer, dimension(:), intent(in) :: Kld,Kcol,Iselection
      integer, intent(in) :: ndof,ccontentType
      logical, intent(in) :: bignoreDiagonal
      integer, dimension(:,:), intent(inout) :: IedgeList

      ! local variables
      logical, dimension(:), allocatable :: BisActive
      integer :: i,j,ij,iedge

      ! Generate set of active entries
      allocate(BisActive(ndof)); BisActive=.false.
      do i = 1, size(Iselection)
        BisActive(Iselection(i))=.true.
      end do

      ! Initialise edge counter
      iedge = 0

      select case (ccontentType)

      case (LSYSSC_EDGELIST_NODESONLY)

        ! Loop over all rows
        do i = 1, size(Kld)-1

          ! Check if this row belongs to an active DOF
          if (.not.BisActive(i)) cycle

          ! Loop over all (off-diagonal) matrix entries IJ
          do ij = Kld(i), Kld(i+1)-1

            ! Get node number j
            j = Kcol(ij)

            ! Check if this edge connects two DOFs which are marked active
            if (.not.BisActive(j)) cycle

            ! Ignore the diagonal entries?
            if (bignoreDiagonal .and. i.eq.j) cycle

            ! Increase edge counter
            iedge = iedge+1

            ! Set node numbers i and j
            IedgeList(1,iedge) = i
            IedgeList(2,iedge) = j
          end do
        end do

      case (LSYSSC_EDGELIST_POSONLY,LSYSSC_EDGELIST_NODESANDPOS)
        call output_line('Positions are not available for unsymmetric matrix!',&
            OU_CLASS_ERROR,OU_MODE_STD,'genEdgesMat9')
        call sys_halt()
      end select

    end subroutine genEdgesMat9Sel

  end subroutine lsyssc_genEdgeList

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_regroupEdgeList(ndof, IedgeList, h_IedgeListIdx,&
                                    ncolor, ccoloringType)

!<description>
    ! This subroutine regroups the entries of a given list of edges
    ! into independent groups. That is, all edges in one group can be
    ! processed simultaneously without the need for synchronization
!</description>

!<input>
    ! Number of degrees of freedom
    integer, intent(in) :: ndof

    ! OPTIONAL: number of colours required to separate the given edge
    ! list into independent
    integer, intent(in), optional :: ncolor

    ! OPTIONAL: constant defining the type of edge-coloring algorithm
    ! Must be one of the LSYSSC_EDGECOLORING_xxx constants
    integer, intent(in), optional :: ccoloringType
!</input>

!<inputoutput>
    ! On input: ungrouped list of edges
    ! On output: list of edges regrouped into independent groups
    integer, dimension(:,:), intent(inout) :: IedgeList

    ! Handle for the edge list separator
    integer, intent(inout) :: h_IedgeListIdx
!</inputoutput>
!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IedgeListIdx
    integer :: nmaxColor,isize,ccoloring

    ccoloring = LSYSSC_EDGECOLORING_GREEDY
    if (present(ccoloringType)) ccoloring = ccoloringType

    if (present(ncolor)) then
      nmaxColor = ncolor
    else
      ! Determine the maximum number of colours required to separate
      ! the given edge list into independent groups
      nmaxColor = compMaxEdgeColors(ndof, IedgeList, ccoloring)
    end if
    
    ! Concistency check
    if (nmaxColor .le. 0) then
      call output_line('Number of colors for edge cloring must be greater than zero!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_regroupEdgeList')
      call sys_halt()
    end if
    
    ! (Re-)allocate memory
    if (h_IedgeListIdx .eq. ST_NOHANDLE) then
      call storage_new('lsyssc_regroupEdgeList','IedgeListIdx',&
          nmaxColor+1, ST_INT, h_IedgeListIdx, ST_NEWBLOCK_NOINIT)
    else
      call storage_getsize(h_IedgeListIdx, isize)
      if (isize .lt. nmaxColor) then
        call storage_free(h_IedgeListIdx)
        call storage_new('lsyssc_regroupEdgeList','IedgeListIdx',&
            nmaxColor+1, ST_INT, h_IedgeListIdx, ST_NEWBLOCK_NOINIT)
      end if
    end if
    
    ! Set pointer
    call storage_getbase_int(h_IedgeListIdx, p_IedgeListIdx)
    
    ! What type of coloring algorithm are we?
    select case(ccoloring)
    case (LSYSSC_EDGECOLORING_GREEDY)
      if (nmaxColor .le. 64) then
        call genEdgeListIdx_greedy64(ndof, nmaxColor, p_IedgeListIdx, IedgeList)
      else
        call output_line('Greedy edge-coloring algorithm for more than '//&
            '64 colors is not implemented yet!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_regroupEdgeList')
        call sys_halt()
      end if

    case (LSYSSC_EDGECOLORING_NTL)
      if (nmaxColor .le. 64) then
        call genEdgeListIdx_ntl64(ndof, nmaxColor, p_IedgeListIdx, IedgeList)
      else
        call output_line('NLT edge-coloring algorithm for more than '//&
            '64 colors is not implemented yet!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_regroupEdgeList')
        call sys_halt()
      end if
      
    case default
      call output_line('Unsupported type of coloring algorithm!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_regroupEdgeList')
      call sys_halt()
    end select
    
  contains
    
    ! Here, the working routines follow

    !**************************************************************
    ! Compute the maximum number of colors required to color the
    ! edges of the matrix such that no two edges with the same color
    ! share a common node (This is the classical edge-coloring)

    function compMaxEdgeColors(ndof, IedgeList, ccoloringType) result(ncolor)

      integer, intent(in) :: ndof, ccoloringType
      integer, dimension(:,:), intent(in) :: IedgeList

      integer :: ncolor

      integer, dimension(:), pointer :: p_IedgesAtDof
      integer :: h_IedgesAtDof
      integer :: i,j,iedge

      ! Allocate temporal memory
      h_IedgesAtDof = ST_NOHANDLE
      call storage_new('compMaxEdgeColors', 'IedgesAtDof', ndof,&
          ST_INT, h_IedgesAtDof, ST_NEWBLOCK_ZERO)
      call storage_getbase_int(h_IedgesAtDof, p_IedgesAtDof)

      ! First shuffle pass:
      ! Compute the number of edges at each DOF
      do iedge = 1, size(IedgeList,2)
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)
        p_IedgesAtDof(i) = p_IedgesAtDof(i)+1
        p_IedgesAtDof(j) = p_IedgesAtDof(j)+1
      end do

      ! Second shuffle pass:
      ! Compute the maximum number of edges at DOF
      ncolor = 0
      do i = 1, ndof
        ncolor = max(ncolor, p_IedgesAtDof(i))
      end do

      ! Deallocate temporal memory
      call storage_free(h_IedgesAtDof)

      ! What type of coloring algorithm are we?
      select case(ccoloringType)
      case (LSYSSC_EDGECOLORING_GREEDY)
        ncolor = 2*ncolor
        
      case (LSYSSC_EDGECOLORING_NTL)
        ncolor = ncolor+1
        
      case default
        ncolor = 0
      end select

    end function compMaxEdgeColors

    !**************************************************************
    ! Generate index for edge data structure based on greedy
    ! edge-coloring algorithm which may require 2*ncolor colors
    ! Note that this implementation works for 64 colors at most

    subroutine genEdgeListIdx_greedy64(neq, nmaxColor, IedgeListIdx, IedgeList)

      integer, intent(in) :: neq,nmaxColor

      integer, dimension(:), intent(inout) :: IedgeListIdx
      integer, dimension(:,:), intent(inout) :: IedgeList

      integer(I64), dimension(:), pointer :: p_IdofColor
      integer, dimension(:,:), pointer :: p_IedgeListTemp
      integer, dimension(2) :: Isize
      integer :: h_IdofColor,h_IedgeListTemp
      integer :: i,j,iedge,jedge,icolor

      if (nmaxColor .gt. 64) then
        call output_line('Current implementation does not support more than 64 colors!',&
          OU_CLASS_ERROR,OU_MODE_STD,'genEdgeListIdx_greedy')
      call sys_halt()
      end if

      ! Allocate temporal memory
      h_IdofColor = ST_NOHANDLE
      call storage_new('genEdgeListIdx', 'IdofColor', neq,&
          ST_INT64, h_IdofColor, ST_NEWBLOCK_ZERO)
      call storage_getbase_int64(h_IdofColor, p_IdofColor)

      ! Make a temporal backup if the edge list
      Isize = shape(IedgeList)
      h_IedgeListTemp = ST_NOHANDLE
      call storage_new('genEdgeListIdx','IedgeListTemp',&
          Isize, ST_INT, h_IedgeListTemp, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int2d(h_IedgeListTemp, p_IedgeListTemp)
      call lalg_copyVector(IedgeList, p_IedgeListTemp)

      ! First shuffle pass: Loop over all edges and determine the
      ! next free color to be used for that edge
      call lalg_clearVector(IedgeListIdx)
      call lalg_clearVector(p_IdofColor)

      edge1: do iedge = 1, size(p_IedgeListTemp,2)
        ! Get node numbers
        i = p_IedgeListTemp(1,iedge)
        j = p_IedgeListTemp(2,iedge)

        ! Find next unused color for nodes I and J
        color1: do icolor = 0, nmaxColor

          if (.not.btest(p_IdofColor(i), icolor) .and.&
              .not.btest(p_IdofColor(j), icolor)) then

            ! Set color marker for nodes I and J
            p_IdofColor(i) = ibset(p_IdofColor(i), icolor)
            p_IdofColor(j) = ibset(p_IdofColor(j), icolor)

            ! Increase number of edges in this color group by one.
            IedgeListIdx(icolor+2) = IedgeListIdx(icolor+2)+1

            ! That is it
            exit color1
          end if
        end do color1
      end do edge1

      ! Second shuffle pass: Compute starting position of each group
      IedgeListIdx(1) = 1
      do icolor = 2, size(IedgeListIdx)
        IedgeListIdx(icolor) = IedgeListIdx(icolor)&
                             + IedgeListIdx(icolor-1)
      end do

      ! Third shuffle pass: Reorder the edges physically
      call lalg_clearVector(p_IdofColor)

       edge2: do iedge = 1, size(p_IedgeListTemp,2)
        ! Get node numbers
        i = p_IedgeListTemp(1,iedge)
        j = p_IedgeListTemp(2,iedge)

        ! Find next unused color for nodes I and J
        color2: do icolor = 0, nmaxColor

          if (.not.btest(p_IdofColor(i), icolor) .and.&
              .not.btest(p_IdofColor(j), icolor)) then

            ! Set color marker for nodes I and J
            p_IdofColor(i) = ibset(p_IdofColor(i), icolor)
            p_IdofColor(j) = ibset(p_IdofColor(j), icolor)

            ! Get next free position
            jedge = IedgeListIdx(icolor+1)

            ! Increase number of edges in this color group by one
            IedgeListIdx(icolor+1) = IedgeListIdx(icolor+1)+1

            ! Copy edge data
            IedgeList(:,jedge) = p_IedgeListTemp(:,iedge)

            ! That is it
            exit color2
          end if
        end do color2
      end do edge2

      ! Fourth shuffle pass:
      ! Adjust the starting index which was tainted in the third pass
      do icolor = size(IedgeListIdx), 2, -1
        IedgeListIdx(icolor) = IedgeListIdx(icolor-1)
      end do
      IedgeListIdx(1) = 1

      ! Deallocate temporal memory
      call storage_free(h_IdofColor)
      call storage_free(h_IedgeListTemp)

    end subroutine genEdgeListIdx_greedy64


    !**************************************************************
    ! Generate index for edge data structure based on recoloring
    ! algorithm due to Nishizeki, Leven and Terada which makes use
    ! of the constructive proof by Vizing to edge-color any graph
    ! with ncolor or ncolor+1 different colours
    ! Note that this implementation works for 64 colors at most
    !
    ! Implementation details: T. Nishizeki, N. Chiba, Planar
    !                         Graphs: Theory and Algorithms.
    !                         Annals of Discrete Mathematic 32

    subroutine genEdgeListIdx_ntl64(neq, nmaxColor, IedgeListIdx, IedgeList)

      integer, intent(in) :: neq,nmaxColor

      integer, dimension(:), intent(inout) :: IedgeListIdx
      integer, dimension(:,:), intent(inout) :: IedgeList

      ! local constants
      integer, parameter, dimension(0:67) :: lookup67 = &
      (/ 64,  0,  1, 39,  2, 15, 40, 23,&
          3, 12, 16, 59, 41, 19, 24, 54,&
          4, -1, 13, 10, 17, 62, 60, 28,&
         42, 30, 20, 51, 25, 44, 55, 47,&
          5, 32, -1, 38, 14, 22, 11, 58,&
         18, 53, 63,  9, 61, 27, 29, 50,&
         43, 46, 31, 37, 21, 57, 52,  8,&
         26, 49, 45, 36, 56,  7, 48, 35,&
          6, 34, 33, -1 /)

      ! local variables
      type(t_arraylistInt_Int) :: rgraph
      type(it_arraylistInt_Int) :: rgraph_iter
      integer(I64), dimension(:), pointer :: p_IdofColor
      integer, dimension(:,:), allocatable :: Ifan
      integer, dimension(:), allocatable :: Iwork
      integer, dimension(:,:), pointer :: p_IApath,p_IedgeListTemp
      integer, dimension(:), pointer :: p_IedgeColor,p_Iedge
      integer, dimension(2) :: Isize
      integer :: h_IedgeColor,h_IdofColor,h_IApath,h_IedgeListTemp
      integer :: iapathLength,ifanLength,i,j,k,k_mis
      integer :: icol1,icol2,icolor,icolor_center,icolor_fan,icolor_last
      integer :: iedge,iedge1,iedge2,iedge_fan,iedge_last,iedge_mis
      integer :: ivert1,ivertex_fan,ivertex_last,ivertex_mis
      logical :: bendAtCenter


      ! There is no need to store a unique missing colour for each
      ! vertex. For each vertex, a 64-bit integer is used to mark if
      ! colour 0..63 is missing or not. The unique missing colour is
      ! simply the least significant bit the is not set.
#define MISSINGCOLOUR(v) lookup67(mod(iand( not(p_IdofColor(v)), -not(p_IdofColor(v))),67_I64))

      ! Create array list for storing adjacency information
      call alst_create(rgraph, neq, 2*size(IedgeList,2), 1)
      call alst_createTbl(rgraph, neq)

      ! Allocate temporal array for storing edge colors
      h_IedgeColor = ST_NOHANDLE
      call storage_new('genEdgeListIdx', 'IedgeColor', size(IedgeList,2),&
          ST_INT, h_IedgeColor, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(h_IedgeColor, p_IedgeColor)
      call lalg_setVector(p_IedgeColor, -1)

      ! Allocate temporal array for storing alternating path
      h_IApath = ST_NOHANDLE
      call storage_new('genEdgeListIdx', 'IApath', (/2,neq/),&
          ST_INT, h_IApath, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int2d(h_IApath, p_IApath)
      
      ! Allocate temporal array for used colors per vertex
      h_IdofColor = ST_NOHANDLE
      call storage_new('genEdgeListIdx', 'IdofColor', neq,&
          ST_INT64, h_IdofColor, ST_NEWBLOCK_ZERO)
      call storage_getbase_int64(h_IdofColor, p_IdofColor)

      ! Allocate working array for constructing maximal fan sequences
      allocate(Iwork(0:nmaxColor))

      ! Allocate array to store vertex and edge numbers of the maximal
      ! fan sequences to be created for each newly introduced edge
      allocate(Ifan(3,0:nmaxColor))

      ! COLOR: Insert each edge of the graph separately
      ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
      do iedge = 1, size(IedgeList,2)
        
        ! Get vertex numbers
        i = IedgeList(1,iedge)
        j = IedgeList(2,iedge)

        ! Perform RECOLORing of new graph G + (i,j)
        ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
               
        ! Loop over all edges adjacent to vertex j and generate a
        ! colour-to-edgenumber lookup table in the working array Iwork
        Iwork = 0; rgraph_iter = alst_begin(rgraph, j)
        do while (rgraph_iter /= alst_end(rgraph, j))

          ! Get edge number
          call alst_getbase_data(rgraph, rgraph_iter, p_Iedge)
          
          ! Get edge color
          icolor = p_IedgeColor(p_Iedge(1))

          ! Mark edge/color combination in working array
          Iwork(icolor) = p_Iedge(1)

          ! Proceed to next edge
          call alst_next(rgraph_iter)          
        end do

        ! Insert edge (i,j) into graph structure
        call alst_push_back(rgraph, i, j, (/iedge/))
        call alst_push_back(rgraph, j, i, (/iedge/))
        
        ! Create maximal fan sequence at vertex j starting with
        ! uncoloured edge (j,i) in time O(nmaxColor), that is
        !
        ! (j,i_0), (j,i_1), ..., (j,i_S)
        !
        ! such that
        !
        ! edge (j,i_k) is coloured with `the' missing color of vertex  i_{k-1}

        ivertex_fan = i
        icolor_fan  = MISSINGCOLOUR(i)

        Ifan(1,0) = ivertex_fan
        Ifan(2,0) = iedge
        Ifan(3,0) = icolor_fan

        fan_loop: do ifanLength = 1, nmaxColor

          ! Get the edge that is coloured with the missing colour of
          ! the previous vertex; zero if no such edge exists. To
          ! ensure that each colour is processed only once we flip the
          ! sign of the working array after having processed
          ! it. Otherwise, we may run into a cycle.
          iedge_fan         = Iwork(icolor_fan)
          Iwork(icolor_fan) = -Iwork(icolor_fan)
          
          ! Can we expand the maximal fan further?
          if (iedge_fan .le. 0) then
            ! Edge has been visited before or no edge with colour
            ! icolor_fan exist; thus, we are done here
            exit fan_loop
          else
            ! Get new vertex number and its missing colour
            ivertex_fan = IedgeList(1,iedge_fan) + IedgeList(2,iedge_fan) - j
            icolor_fan  = MISSINGCOLOUR(ivertex_fan)

            ! Store data in fan sequence
            Ifan(1,ifanLength) = ivertex_fan
            Ifan(2,ifanLength) = iedge_fan
            Ifan(3,ifanLength) = icolor_fan
          end if
        end do fan_loop
        
        ! Adjust length of maximal fan sequence
        ifanLength = ifanLength-1

        ! Now, the array Ifan(0:ifanLength) contains the vertex
        ! numbers which represent the maximal fan sequence centered
        ! about vertex j starting at vertex i.

        ivertex_last = Ifan(1,ifanLength)
        iedge_last   = Ifan(2,ifanLength)
        icolor_last  = Ifan(3,ifanLength)

        ! Variable ivertex_last represents the number of the last
        ! vertex in the maximal fan sequence and icolor_last is `the'
        ! missing colour of the last vertex in the fan sequence

        if (.not.btest(p_IdofColor(j), icolor_last)) then
          ! Case (1): `the' missing colour of the last vertex with
          ! number ivertex_last in the fan is a missing colour of vertex j
          ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          
          ! Simply shift all colours of fan sequence by one
          do k = 0, ifanLength

            ! Get vertex and edge number
            ivert1 = Ifan(1,k)
            iedge1 = Ifan(2,k)

            ! Make backup of the previous edge colour
            icol1 = p_IedgeColor(iedge1)

            ! Color edge (j,i_k) with the missing colour of i_k
            p_IedgeColor(iedge1) = Ifan(3,k)
            
            ! Update (=clear old,set new) missing colours of vertex i_k
            if (k .ne. 0)  p_IdofColor(ivert1) = ibclr(p_IdofColor(ivert1), icol1)
            p_IdofColor(ivert1) = ibset(p_IdofColor(ivert1), p_IedgeColor(iedge1))
          end do

          ! Update missing colour of fan centering vertex j only for
          ! last edge (i,ivert1) since colour icolor_last has not been
          ! used for colouring the edges within the maximal fan sequence
          p_IdofColor(j) = ibset(p_IdofColor(j), icolor_last)
          
        else
          ! Case (2): an edge of the maximal fan sequence is coloured
          ! with `the' missing colour of the last vertex numbered ivertex
          ! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          
          do k = 0, ifanLength-2
            if (Ifan(3,k) .eq. icolor_last) exit
          end do
          ivertex_mis = Ifan(1,k)
          iedge_mis   = Ifan(2,k)
          k_mis       = k

          ! Now, ivertex_mis is the number of the vertex such that
          ! `the' missing colours of vertices ivertex_mis and the last
          ! vertex in the fan sequence ivertex_last coincide.
          
          ! Construct alternating path with colours icolor_last and
          ! the missing colour icolor_center of the fan center
          ! vertex j starting at the intermediate vertex ivertex_mis
          icolor_center = MISSINGCOLOUR(j)
          iapathLength  = 0
          ivert1        = ivertex_mis
          icol1         = icolor_center
          iedge1        = 0
          apath_loop1: do
            
            rgraph_iter = alst_begin(rgraph, ivert1)
            do while (rgraph_iter /= alst_end(rgraph, ivert1))
              
              ! Get edge number
              call alst_getbase_data(rgraph, rgraph_iter, p_Iedge)
              iedge2 = p_Iedge(1)
              
              ! Check if this edge belongs to the alternating path and
              ! has not been traversed by the alternating path before
              if ((iedge1 .ne. iedge2) .and.&
                  (icol1  .eq. p_IedgeColor(iedge2))) then
                
                ! Update new edge and starting vertex
                iedge1 = iedge2
                ivert1 = IedgeList(1,iedge1) + IedgeList(2,iedge1) - ivert1
                
                ! Update path length
                iapathLength = iapathLength+1
                
                ! Store edge number and colour
                p_IApath(1,iapathLength) = iedge1
                p_IApath(2,iapathLength) = icol1
                
                ! Alternate colour
                icol1 = icolor_center + icolor_last - icol1
                
                cycle apath_loop1
              end if
              
              ! Proceed to next edge
              call alst_next(rgraph_iter)
            end do
            
            ! No further edges could be found
            exit apath_loop1
          end do apath_loop1

          ! Check if path ends at fan centering vertex j
          if (iapathLength .gt. 0) then
            iedge1 = p_IApath(1,iapathLength)
            bendAtCenter = (IedgeList(1,iedge1) .eq. j)&
                      .or. (IedgeList(2,iedge1) .eq. j)
          else
            ! There is no path that could end at vertex j
            bendAtCenter = .false.
          end if

          if (bendAtCenter) then
            ! The alternating path starting at vertex ivertex_last
            ! ends at the fan centering vertex j; thus we have to
            ! construct another alternating path
            
            ! Construct alternating path with colours icolor_last and
            ! `the' missing colour icolor_center of the fan center
            ! vertex j starting at the last vertex ivertex_last
            icolor_center = MISSINGCOLOUR(j)
            iapathLength  = 0
            ivert1        = ivertex_last
            icol1         = icolor_center
            iedge1        = 0
            apath_loop2: do
              
              rgraph_iter = alst_begin(rgraph, ivert1)
              do while (rgraph_iter /= alst_end(rgraph, ivert1))
                
                ! Get edge number of the next incident edge
                call alst_getbase_data(rgraph, rgraph_iter, p_Iedge)
                iedge2 = p_Iedge(1)
                
                ! Check if this edge belongs to the alternating path,
                ! that is, it is coloured with colour icol1 and has not
                ! been traversed by the alternating path before
                if ((iedge1 .ne. iedge2) .and.&
                    (icol1  .eq. p_IedgeColor(iedge2))) then
                  
                  ! Update new edge and  starting vertex
                  iedge1 = iedge2
                  ivert1 = IedgeList(1,iedge1) + IedgeList(2,iedge1) - ivert1
                  
                  ! Update path length
                  iapathLength = iapathLength+1
                  
                  ! Store edge number and colour
                  p_IApath(1,iapathLength) = iedge1
                  p_IApath(2,iapathLength) = icol1
                  
                  ! Alternate colour
                  icol1 = icolor_center + icolor_last - icol1
                  
                  cycle apath_loop2
                end if
                
                ! Proceed to next edge
                call alst_next(rgraph_iter)
              end do
              
              ! No further edges could be found
              exit apath_loop2
            end do apath_loop2

            ! The alternating path starting at vertex ivertex_last
            ! does not end at the fan centering vertex j

            ! Interchange the colors of the alternating path
            do k = 1, iapathLength
              ! Get edge number
              iedge1 = p_IApath(1,k)
              
              ! Determine current and new color of the ede
              icol1 = icolor_center + icolor_last - p_IApath(2,k)
              
              ! Update edge color
              p_IedgeColor(iedge1) = icol1
            end do
          
            if (iapathLength .gt. 0) then
              ! The alternating path does not form one large circle so
              ! that we need to update the list of missing colours at
              ! the start and end vertices of the path accordingly
              p_IdofColor(ivertex_last) = ibclr(p_IdofColor(ivertex_last), icolor_center)
              p_IdofColor(ivertex_last) = ibset(p_IdofColor(ivertex_last), icolor_last)
              
              icol1 = p_IApath(2,iapathLength)
              icol2 = icolor_center + icolor_last - icol1
              
              p_IdofColor(ivert1) = ibclr(p_IdofColor(ivert1), icol1)
              p_IdofColor(ivert1) = ibset(p_IdofColor(ivert1), icol2)
            end if

            ! The missing colour of vertex j is now missing both at j
            ! and ivertex_last; Shift all colours of the fan sequence
            ! by one; note that the zeroth edge is treated separately
            do k = 0, ifanLength-1
              ! Get vertex and edge number
              ivert1 = Ifan(1,k)
              iedge1 = Ifan(2,k)
              
              ! Make backup of the previous edge colour
              icol1 = p_IedgeColor(iedge1)
              
              ! Color edge (j,i_k) with the missing colour of i_k
              p_IedgeColor(iedge1) = Ifan(3,k)
              
              ! Update (=clear old,set new) missing colours of vertex i_k
              if (k .ne. 0) p_IdofColor(ivert1) = ibclr(p_IdofColor(ivert1), icol1)
              p_IdofColor(ivert1) = ibset(p_IdofColor(ivert1), p_IedgeColor(iedge1))

              if (k .eq. ifanLength-1) then
                ! Update missing colour of fan centering vertex j only for
                ! last edge (i,ivert1) since colour icolor_last has not been
                ! used for colouring the edges within the maximal fan sequence
                p_IdofColor(j) = ibset(p_IdofColor(j), p_IedgeColor(iedge1))
              end if
            end do

            ! Color edge (j,ivertex_last) with the missing colour
            ! icolour_center of the fan centering vertex j
            p_IdofColor(ivertex_last) = ibclr(p_IdofColor(ivertex_last), p_IedgeColor(iedge_last))
            p_IedgeColor(iedge_last)  = icolor_center
            p_IdofColor(ivertex_last) = ibset(p_IdofColor(ivertex_last), icolor_center)
            p_IdofColor(j)            = ibset(p_IdofColor(j),           icolor_center)

          else
            ! The alternating path starting at vertex ivertex_last
            ! does not end at the fan centering vertex j

            ! Interchange the colors of the alternating path
            do k = 1, iapathLength
              ! Get edge number
              iedge1 = p_IApath(1,k)
              
              ! Determine current and new color of the ede
              icol1 = icolor_center + icolor_last - p_IApath(2,k)
              
              ! Update edge color
              p_IedgeColor(iedge1) = icol1
            end do

            if (iapathLength .gt. 0) then
              ! The alternating path does not form one large circle so
              ! that we need to update the list of missing colours at
              ! the start and end vertices of the path accordingly
              p_IdofColor(ivertex_mis) = ibclr(p_IdofColor(ivertex_mis), icolor_center)
              p_IdofColor(ivertex_mis) = ibset(p_IdofColor(ivertex_mis), icolor_last)
              
              icol1 = p_IApath(2,iapathLength)
              icol2 = icolor_center + icolor_last - icol1
                            
              p_IdofColor(ivert1) = ibclr(p_IdofColor(ivert1), icol1)
              p_IdofColor(ivert1) = ibset(p_IdofColor(ivert1), icol2)
            end if

            ! The missing colour of vertex j is now missing both at j
            ! and ivertex_mis; Shift all colours of the fan sequence
            ! by one; note that the zeroth edge is treated separately
            do k = 0, k_mis-1
              ! Get vertex and edge number
              ivert1 = Ifan(1,k)
              iedge1 = Ifan(2,k)
              
              ! Make backup of the previous edge colour
              icol1 = p_IedgeColor(iedge1)
              
              ! Color edge (j,i_k) with the missing colour of i_k
              p_IedgeColor(iedge1) = Ifan(3,k)
                           
              ! Update (=clear old,set new) missing colours of vertex i_k
              if (k .ne. 0) p_IdofColor(ivert1) = ibclr(p_IdofColor(ivert1), icol1)
              p_IdofColor(ivert1) = ibset(p_IdofColor(ivert1), p_IedgeColor(iedge1))

              if (k .eq. k_mis-1) then
                ! Update missing colour of fan centering vertex j only for
                ! last edge (i,ivert1) since colour icolor_last has not been
                ! used for colouring the edges within the maximal fan sequence
                p_IdofColor(j) = ibset(p_IdofColor(j), p_IedgeColor(iedge1))
              end if
            end do
            
            ! Color edge (j,ivertex_mis) with the missing colour
            ! icolour_center of the fan centering vertex j
            p_IdofColor(ivertex_mis) = ibclr(p_IdofColor(ivertex_mis), p_IedgeColor(iedge_mis))
            p_IedgeColor(iedge_mis)  = icolor_center
            p_IdofColor(ivertex_mis) = ibset(p_IdofColor(ivertex_mis), icolor_center)
            p_IdofColor(j)           = ibset(p_IdofColor(j),           icolor_center)
          end if
        end if
        
      end do
      
      ! Release temporal array lists
      call alst_release(rgraph)

      ! Deallocate temporal arrays
      deallocate(Iwork)
      deallocate(Ifan)

      ! First shuffle pass: Loop over all edges and count the number
      ! of edges per colour group
      call lalg_clearVector(IedgeListIdx)
      edge1: do iedge = 1, size(IedgeList,2)
        icolor = p_IedgeColor(iedge)
        IedgeListIdx(icolor+2) = IedgeListIdx(icolor+2)+1
      end do edge1

      ! Second shuffle pass: Compute starting position of each group
      IedgeListIdx(1) = 1
      do icolor = 2, size(IedgeListIdx)
        IedgeListIdx(icolor) = IedgeListIdx(icolor)&
                             + IedgeListIdx(icolor-1)
      end do

      ! Make a temporal backup if the edge list
      Isize = shape(IedgeList)
      h_IedgeListTemp = ST_NOHANDLE
      call storage_new('genEdgeListIdx','IedgeListTemp',&
          Isize, ST_INT, h_IedgeListTemp, ST_NEWBLOCK_NOINIT)
      call storage_getbase_int2d(h_IedgeListTemp, p_IedgeListTemp)
      call lalg_copyVector(IedgeList, p_IedgeListTemp)

      ! Third shuffle pass: Reorder the edges physically
      edge2: do iedge = 1, size(IedgeList,2)
        ! Get edge colour and first free position
        icolor = p_IedgeColor(iedge)+1
        iedge1 = IedgeListIdx(icolor)

        ! Copy edge data into its correct position
        IedgeList(:,iedge1) = p_IedgeListTemp(:,iedge)
        IedgeListIdx(icolor) = IedgeListIdx(icolor)+1
      end do edge2

      ! Fourth shuffle pass:
      ! Adjust the starting index which was tainted in the third pass
      do icolor = size(IedgeListIdx), 2, -1
        IedgeListIdx(icolor) = IedgeListIdx(icolor-1)
      end do
      IedgeListIdx(1) = 1
      
      ! Release temporal memory
      call storage_free(h_IApath)
      call storage_free(h_IdofColor)
      call storage_free(h_IedgeColor)
      call storage_free(h_IedgeListTemp)
      
    end subroutine genEdgeListIdx_ntl64

  end subroutine lsyssc_regroupEdgeList

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_calcDimsFromMatrix(rmatrix, na, neq, ncols, nedge,&
      IdofList, IdofsTest, IdofsTrial)

!<description>
    ! This subroutine calculates the dimensions of the matrix:
    ! - number of non-zero matrix entries NA
    ! - number of equations NEQ
    ! - number of columns NCOLS
    ! - number of edges NEDGE (=0 if NEQ /= NCOLS)
    !
    ! If the optional parameters IdofList, IdofsTest or IdofsTrial are
    ! given, then their entries are used as restriction set.
    !
    ! If IdofList is given, then its entries are used both for rows
    ! and columns. If IdofsTest AND IdofsTrial are given, then different
    ! restriction sets are considered for the rows and the columns.
    ! If either IdofsTest OR IdofsTrial are given, then it is assumed
    ! that only the rows or columns are restricted.
!</description>

!<input>
    ! Scalar matrix
    type(t_matrixScalar), intent(in) :: rmatrix

    ! OPTIONAL: a list of degress of freedom to which the edge
    ! structure should be restricted.
    integer, dimension(:), intent(in), optional :: IdofList
    integer, dimension(:), intent(in), optional :: IdofsTest
    integer, dimension(:), intent(in), optional :: IdofsTrial
!</intput>

!<output>
    ! Number of non-zero matrix entries
    integer, intent(out) :: na

    ! Number of equations
    integer, intent(out) :: neq

    ! Number of columns
    integer, intent(out) :: ncols

    ! Number of edges
    integer, intent(out) :: nedge
!</output>
!</subroutine>

    ! local variables
    logical, dimension(:), pointer :: BisActiveRow,BisActiveColumn
    integer, dimension(:), pointer :: p_Kld, p_Kcol
    integer :: i,ij,j

    if (rmatrix%NEQ .ne. rmatrix%NCOLS) then
      call output_line('Matrix must be square matrix!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_calcDimsFromMatrix')
      call sys_halt()
    end if

    ! Do we have any restrictions?
    if (.not.(present(IdofList))  .and.&
        .not.(present(IdofsTest)) .and.&
        .not.(present(IdofsTrial))) then

      na    = rmatrix%NA
      neq   = rmatrix%NEQ
      ncols = rmatrix%NCOLS

      if (neq .eq. ncols) then
        nedge = (na-neq)/2
      else
        nedge = 0
      end if

    else

      ! Do we have a common list of DOFs for rows and columns?
      if (present(IdofList)) then

        ! Use the same set of restrictions for rows and columns
        allocate(BisActiveRow(max(rmatrix%NEQ,rmatrix%NCOLS)))
        BisActiveRow=.false.; BisActiveColumn => BisActiveRow

        ! Set number of rows and columns
        neq = size(IdofList); ncols = neq

        ! Generate set of active rows
        do i = 1, size(IdofList)
          BisActiveRow(IdofList(i)) = .true.
        end do

      else ! use individual lists of DOFs for rows and columns

        allocate(BisActiveRow(rmatrix%NEQ))
        allocate(BisActiveColumn(rmatrix%NCOLS))

        if (present(IdofsTest)) then
          ! Set number of rows
          neq = size(IdofsTest)

          ! Generate set of active rows
          BisActiveRow=.false.
          do i = 1, size(IdofsTest)
            BisActiveRow(IdofsTest(i)) = .true.
          end do
        else
          ! Set number of rows
          neq = rmatrix%NEQ

          ! No restriction for rows
          BisActiveRow=.true.
        end if

        if (present(IdofsTrial)) then
          ! Set number of columns
          ncols = size(IdofsTrial)

          ! Generate set of active columns
          BisActiveColumn=.false.
          do i = 1, size(IdofsTrial)
            BisActiveColumn(IdofsTrial(i)) = .true.
          end do
        else
          ! Set number of columns
          ncols = rmatrix%NCOLS

          ! No restriction for columns
          BisActiveColumn=.true.
        end if
      end if

      !-------------------------------------------------------------------------

      ! Initialisation
      na = 0

      ! What type of matrix are we?
      select case (rmatrix%cmatrixFormat)
      case (LSYSSC_MATRIX1)

        ! Determine number of non-zero matrix entries, equations and edges
        do i = 1, rmatrix%NEQ
          if (.not.BisActiveRow(i)) cycle
          do j = 1, rmatrix%NCOLS
            if (.not.BisActiveColumn(j)) cycle
            na = na+1
          end do
        end do

        ! Compute number of edges
        if (neq .eq. ncols) then
          nedge = (na-neq)/2
        else
          nedge = 0
        end if

      case (LSYSSC_MATRIX7, LSYSSC_MATRIX7INTL,&
            LSYSSC_MATRIX9, LSYSSC_MATRIX9INTL)

        ! Set pointers
        call lsyssc_getbase_Kld(rmatrix, p_Kld)
        call lsyssc_getbase_Kcol(rmatrix, p_Kcol)

        ! Determine number of non-zero matrix entries, equations and edges
        do i=1, rmatrix%NEQ
          if (.not.BisActiveRow(i)) cycle
          do ij = p_Kld(i), p_Kld(i+1)-1
            j = p_Kcol(ij)
            if (.not.BisActiveColumn(j)) cycle
            na = na+1
          end do
        end do

        ! Compute number of edges
        if (neq .eq. ncols) then
          nedge = (na-neq)/2
        else
          nedge = 0
        end if

      case default
        call output_line('Unsupported matrix format!',&
            OU_CLASS_WARNING,OU_MODE_STD,'lsyssc_calcDimsFromMatrix1')
        call sys_halt()
      end select

      ! Deallocate temporal memory
      if (present(IdofList)) then
        deallocate(BisActiveRow)
        nullify(BisActiveColumn)
      else
        deallocate(BisActiveRow,BisActiveColumn)
      end if

    end if

  end subroutine lsyssc_calcDimsFromMatrix

  ! ***************************************************************************

!<subroutine>
  
  subroutine lsyssc_calcGerschgorin (rmatrix,dminEigenval,dmaxEigenval)
  
!<description>
  ! Applies the theorem of Gerschgorin to compute an approximation
  ! of the minimum and maximum eigenvalue of a matrix.
!</description>
  
!<input>
  ! A matrix
  type(t_matrixScalar), intent(in) :: rmatrix
!</input>

!<output>
  ! Approximation for the real part of the minimum eigenvalue
  real(DP), intent(out) :: dminEigenval

  ! Approximation for the real part of the maximum eigenvalue
  real(DP), intent(out) :: dmaxEigenval
!</output>
  
!</subroutine>
    ! local variables
    real(DP), dimension(:), pointer :: p_Da
    integer, dimension(:), pointer :: p_Kld
    integer, dimension(:), pointer :: p_Kdiagonal
    integer :: i,j
    real(DP) :: drad,da

    ! At the moment, we only support format-9 matrices, double
    ! precision.
    
    dminEigenval = SYS_MAXREAL_DP
    dmaxEigenval = 0.0_DP

    if (rmatrix%cdataType .ne. ST_DOUBLE) then
      call output_line("Data type not supported",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_calcGerschgorin")
      call sys_halt()
    end if
    
    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX9)
      
      ! Get the data arrays
      call lsyssc_getbase_double (rmatrix,p_Da)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)
      call lsyssc_getbase_Kdiagonal (rmatrix,p_Kdiagonal)
      
      ! Loop through the rows
      do i=1,rmatrix%neq
      
        ! Calculate the Gerschgorin radius
        drad = 0
        do j=p_Kld(i),p_Kld(i+1)-1
          drad = drad + abs(p_Da(j))
        end do
        ! Subtract the diagonal, it does not belong to the radius
        da = abs(p_Da(p_Kdiagonal(i)))
        drad = drad - da
        
        ! Calculate the approximation for the eigenvalues
        dminEigenval = min(dminEigenval,da-drad)
        dmaxEigenval = max(dmaxEigenval,da-drad)
      
      end do
      
    case default
      call output_line("Unsupported matrix format!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_calcGerschgorin")
      call sys_halt()
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_createGraphFromMatrix(rmatrix, rgraph, nvtMax, nedgeMax)

!<description>
    ! This routine creates a an adjacency graph rgraph from the scalar
    ! matrix rmatrix. If the optional parameters nvtMax or nedgeMax
    ! are present then their values will be adopted for the maximum
    ! number of vertices/edges that can be stored in the graph
    ! initially without reallocation.
!</description>

!<input>
    ! The scalar matrix which should be used to create the sparsity graph
    type(t_matrixScalar), intent(in) :: rmatrix

    ! OPTIONAL: maximum number of vertices
    integer, intent(in), optional :: nvtMax

    ! OPTIONAL: maximum number of edges
    integer, intent(in), optional :: nedgeMax
!</input>

!<output>
    ! The adjacency graph
    type(t_graph), intent(out) :: rgraph
!</output>

!</subroutine>

    ! local variables
    integer :: h_Key
    integer :: nvt,nedge
    integer, dimension(:), pointer :: p_Kld,p_Kcol,p_Key

    ! Set dimensions
    rgraph%NVT   = rmatrix%NEQ
    rgraph%NEDGE = rmatrix%NA

    ! Estimate maximum number of vertices
    if (present(nvtMax)) then
      nvt = max(nvtMax,rgraph%NVT)
    else
      nvt = int(1.5_DP*rgraph%NVT)
    end if

    ! Estimate maximum number of edges
    if (present(nedgeMax)) then
      nedge = max(nedgeMax,rgraph%NEDGE)
    else
      nedge = int(0.5_DP*(rgraph%NEDGE+rgraph%NVT))
    end if

    ! What matrix format are we?
    select case(rmatrix%cmatrixFormat)

    case(LSYSSC_MATRIX7,&
         LSYSSC_MATRIX7INTL)
      rgraph%cgraphFormat = GRPH_GRAPH7

      ! Create map for the list of vertices
      call map_create(rgraph%rVertices,rgraph%NVT+1,0)

      ! Create array of lists for edges
      call alst_create(rgraph%rEdges,nvt,nedge)

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrix,p_Kld)
      call lsyssc_getbase_Kcol(rmatrix,p_Kcol)

      ! Fill list of edges
      call alst_copy(p_Kld,p_Kcol,rgraph%rEdges)

      ! Generate p_Key = array [1,2,3,...,NVT]
      call storage_new('lsyssc_createGraphFromMatrix','p_Key',&
          rgraph%NVT,ST_INT, h_Key,ST_NEWBLOCK_ORDERED)
      call storage_getbase_int(h_Key,p_Key)
      call map_copy(p_Key,rgraph%rVertices)
      call storage_free(h_Key)


    case(LSYSSC_MATRIX9,&
         LSYSSC_MATRIX9INTL)
      rgraph%cgraphFormat = GRPH_GRAPH9

      ! Create map for the list of vertices
      call map_create(rgraph%rVertices,rgraph%NVT+1,0)

      ! Create array of lists for edges
      call alst_create(rgraph%rEdges,nvt,nedge)

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrix,p_Kld)
      call lsyssc_getbase_Kcol(rmatrix,p_Kcol)

      ! Fill list of edges
      call alst_copy(p_Kld,p_Kcol,rgraph%rEdges)

      ! Generate p_Key = array [1,2,3,...,NVT]
      call storage_new('lsyssc_createGraphFromMatrix','p_Key',&
          rgraph%NVT,ST_INT, h_Key,ST_NEWBLOCK_ORDERED)
      call storage_getbase_int(h_Key,p_Key)
      call map_copy(p_Key,rgraph%rVertices)
      call storage_free(h_Key)


    case DEFAULT
      call output_line('Invalid matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_createGraphFromMatrix')
      call sys_halt()
    end select

  end subroutine lsyssc_createGraphFromMatrix

  ! ***************************************************************************

!<subroutine>

  subroutine lsyssc_createMatrixFromGraph(rgraph, rmatrix)

!<description>
    ! This subroutine generates a sparse matrix rmatrix from the
    ! sparsity pattern given by the adjacency graph rgraph.
!</description>

!<input>
    ! The adjacency graph
    type(t_graph), intent(inout) :: rgraph
!</input>

!<inputoutput>
    ! The matrix
    type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    type(it_mapInt_Int) :: rmapIter
    integer, dimension(:), pointer :: p_Idata,p_Kcol,p_Kdiagonal,p_Kld
    integer :: ia,ieq,isize,itable,ncols

    ! Check that matrix and graph have the same format
    select case(rmatrix%cmatrixFormat)

    case(LSYSSC_MATRIX7,&
         LSYSSC_MATRIX7INTL)
      if (rgraph%cgraphFormat .ne. GRPH_GRAPH7) then
        call output_line('Matrix/graph have incompatible format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_createMatrixFromGraph')
        call sys_halt()
      end if

    case(LSYSSC_MATRIX9,&
         LSYSSC_MATRIX9INTL)
      if (rgraph%cgraphFormat .ne. GRPH_GRAPH9) then
        call output_line('Matrix/graph have incompatible format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_createMatrixFromGraph')
        call sys_halt()
      end if

    case DEFAULT
      call output_line('Unsupported matrix format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsyssc_createMatrixFromGraph')
      call sys_halt()
    end select

    ! Set number of edges = number of nonzero matrix entries
    rmatrix%NA = rgraph%NEDGE

    ! Set number of columns/rows. This is a little bit ugly because
    ! the number of vertices (NVT) may be different from the number of
    ! tables.  If vertices are inserted as follows 1,2,3,...,NVT, then
    ! the number of tables is equal to the number of vertices.
    ! However, if only vertices 1,5,9,10 are present in the graph,
    ! then NVT=4 but the largest number in the graph is 10. Hence, the
    ! matrix has NEQ=NCOLS=10 and some rows, e.g.  2,3,4 have no
    ! entries. For a finite element matrix, this does not make sense
    ! but this graph module should be as general as possible.
    rmatrix%NEQ   = max(rgraph%NVT,alst_ntable(rgraph%rEdges))
    rmatrix%NCOLS = rmatrix%NEQ

    if (rgraph%bisDense) then

      ! Convert array list to matrix
      call alst_copy(rgraph%rEdges,rmatrix%h_Kld,rmatrix%h_Kcol)

    else

      ! Check if matrix is empty
      if ((rmatrix%NEQ .eq. 0) .or. rmatrix%NA .eq. 0) return

      ! Convert array list step-by-step
      if (rmatrix%h_Kld .eq. ST_NOHANDLE) then
        call storage_new('lsyssc_createMatrixFromGraph','p_Kld',&
            rmatrix%NEQ+1,ST_INT,rmatrix%h_Kld,ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(rmatrix%h_Kld,isize)
        if (isize < rmatrix%NEQ+1) then
          call storage_realloc('lsyssc_createMatrixFromGraph',&
              rmatrix%NEQ+1,rmatrix%h_Kld,ST_NEWBLOCK_NOINIT,.false.)
        end if
      end if

      if (rmatrix%h_Kcol .eq. ST_NOHANDLE) then
        call storage_new('lsyssc_createMatrixFromGraph','p_Kcol',&
            rmatrix%NA,ST_INT,rmatrix%h_Kcol,ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(rmatrix%h_Kcol,isize)
        if (isize < rmatrix%NA) then
          call storage_realloc('lsyssc_createMatrixFromGraph',&
              rmatrix%NA,rmatrix%h_Kcol,ST_NEWBLOCK_NOINIT,.false.)
        end if
      end if

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrix,p_Kld)
      call lsyssc_getbase_Kcol(rmatrix,p_Kcol)

      ! Initialization
      ia=1

      ! Loop over all equations
      do ieq=1,rmatrix%NEQ

        p_Kld(ieq) = ia

        ! Check if vertex with number IEQ exists
        rmapIter = map_find(rgraph%rVertices,ieq)
        if (.not.map_isNull(rmapIter)) then

          ! Get auxiliary data
          call map_getbase_data(rgraph%rVertices,rmapIter,p_Idata)
          itable = p_Idata(1)

          ! Restore row from table
          call alst_copyTbl(rgraph%rEdges,itable,p_Kcol(ia:),ncols)
        end if

        ia = ia+ncols
      end do
      p_Kld(rmatrix%NEQ+1)=ia
    end if

    ! Do we have to rebuild the diagonal?
    if (rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX9 .or.&
        rmatrix%cmatrixFormat .eq. LSYSSC_MATRIX9INTL) then

      ! Create new memory or resize existing memory
      if (rmatrix%h_Kdiagonal .eq. ST_NOHANDLE) then
        call storage_new('lsyssc_createMatrixFromGraph','p_Kdiagonal',&
            rmatrix%NEQ, ST_INT, rmatrix%h_Kdiagonal,ST_NEWBLOCK_NOINIT)
      else
        call storage_getsize(rmatrix%h_Kdiagonal,isize)
        if (isize < rmatrix%NEQ) then
          call storage_realloc('lsyssc_createMatrixFromGraph',&
              rmatrix%NEQ, rmatrix%h_Kdiagonal, ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if

      ! Set pointers
      call lsyssc_getbase_Kld(rmatrix, p_Kld)
      call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
      call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiagonal)

      ! Rebuild array
      call lsyssc_rebuildKdiagonal(p_Kcol, p_Kld, p_Kdiagonal, rmatrix%NEQ)
    end if

  end subroutine lsyssc_createMatrixFromGraph

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsyssc_matrixLinearCombIndexed (rmatrixA,rmatrixB,ca,cb,Irows)

    !<description>
    ! Performs a linear combination:
    !   rmatrixB = ca*rmatrixA + cB*rmatrixB
    !
    ! Only those rows indexed by Irows are processed. 
    ! Both matrices must have the same matrix format and structure.
    ! </description>

!<input>
    ! Source matrix
    type(t_matrixScalar), intent(in) :: rmatrixA

    ! scaling factors
    real(DP), intent(in) :: ca,cb

    ! OPTIONAL: List of rows to be processed.
    ! If not specified, all rows are processed.
    integer, dimension(:), intent(in), optional :: Irows
!</input>

!<inputoutput>
    ! OPTIONAL: Destination matrix.
    type(t_matrixScalar), intent(inout), target, optional :: rmatrixB
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata1, p_Ddata2
    integer, dimension(:), pointer :: p_Kld
    integer :: i,j

    ! Check if all matrices are compatible
    if ((rmatrixA%NEQ   .ne. rmatrixB%NEQ) .or.&
        (rmatrixA%NCOLS .ne. rmatrixB%NCOLS)) then
      call output_line("Number of rows/columns is not compatible!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
      call sys_halt()
    end if

    ! Matrices must have the same formatcmatrixFormat
    if (rmatrixA%cmatrixFormat   .ne. rmatrixB%cmatrixFormat) then
      call output_line("Matrices have different format!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
      call sys_halt()
    end if

    ! Check if all matrices have the same sorting
    if (rmatrixA%bcolumnsSorted .or. rmatrixB%bcolumnsSorted) then
      if (rmatrixA%p_rsortStrategyColumns%ctype .ne. rmatrixB%p_rsortStrategyColumns%ctype) then
        call output_line("Incompatible sorting strategies!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
        call sys_halt()
      end if
    end if

    if (rmatrixA%browsSorted .or. rmatrixB%browsSorted) then
      if (rmatrixA%p_rsortStrategyRows%ctype .ne. rmatrixB%p_rsortStrategyRows%ctype) then
        call output_line("Incompatible sorting strategies!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
        call sys_halt()
      end if
    end if

    ! Check if destination matrix has compatible data type
    if ((rmatrixA%cdatatype .ne. ST_DOUBLE) .or.&
        (rmatrixB%cdatatype .ne. ST_DOUBLE)) then
      call output_line("Data type of destination matrix not supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
      call sys_halt()
    end if
    
    select case (rmatrixA%cmatrixFormat)
    case (LSYSSC_MATRIX9)

      select case (rmatrixA%cmatrixFormat)
      case (LSYSSC_MATRIX9)
      
        ! Get the data pointers
        call lsyssc_getbase_double (rmatrixA,p_Ddata1)
        call lsyssc_getbase_double (rmatrixB,p_Ddata2)
        call lsyssc_getbase_Kld (rmatrixA,p_Kld)
        
        ! Linear combination of the rows
        if (present(Irows)) then
          ! Subset of the rows
          if ((ca .eq. 0.0_DP) .and. (cb .eq. 0.0_DP)) then
            do i=1,size(Irows)
              do j=p_Kld(Irows(i)), p_Kld(Irows(i+1)-1)
                p_Ddata2(j) = 0.0_DP
              end do
            end do
          else if (ca .eq. 0.0_DP) then
            do i=1,size(Irows)
              do j=p_Kld(Irows(i)), p_Kld(Irows(i)+1)-1
                p_Ddata2(j) = cb * p_Ddata2(j)
              end do
            end do
          else if (cb .eq. 0.0_DP) then
            do i=1,size(Irows)
              do j=p_Kld(Irows(i)), p_Kld(Irows(i)+1)-1
                p_Ddata2(j) = ca * p_Ddata1(j)
              end do
            end do
          else
            do i=1,size(Irows)
              do j=p_Kld(Irows(i)), p_Kld(Irows(i)+1)-1
                p_Ddata2(j) = ca * p_Ddata1(j) + cb * p_Ddata2(j)
              end do
            end do
          end if
        
        else
        
          ! All rows
          if ((ca .eq. 0.0_DP) .and. (cb .eq. 0.0_DP)) then
            do j=1,rmatrixA%NA
              p_Ddata2(j) = 0.0_DP
            end do
          else if (ca .eq. 0.0_DP) then
            do j=1,rmatrixA%NA
              p_Ddata2(j) = cb * p_Ddata2(j)
            end do
          else if (cb .eq. 0.0_DP) then
            do j=1,rmatrixA%NA
              p_Ddata2(j) = ca * p_Ddata1(j)
            end do
          else
            do j=1,rmatrixA%NA
              p_Ddata2(j) = ca * p_Ddata1(j) + cb * p_Ddata2(j)
            end do
          end if
        end if        
      
      case default
        call output_line(&
            "Combination of source and destination matrix format not supported!",&
            OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
        call sys_halt()
      end select

    case default
      call output_line(&
          "Combination of source and destination matrix format not supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"lsyssc_matrixLinearCombIndexed")
      call sys_halt()
    end select      

  end subroutine

end module

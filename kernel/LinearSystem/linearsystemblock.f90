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
!#  3.) lsysbl_createVecBlockByDiscr
!#      -> Create a block vector using a block discretisation structure
!#
!#  4.) lsysbl_createVecBlockIndMat
!#      -> Create a block vector according to a block matrix.
!#
!#  5.) lsysbl_createMatFromScalar
!#      -> Create a 1x1 block matrix from a scalar matrix
!#
!#  6.) lsysbl_createVecFromScalar
!#      -> Create a 1-block vector from a scalar vector
!#
!#  7.) lsysbl_createMatBlockByDiscr
!#      -> Create an empty block matrix using a block discretisation structure
!#
!#  8.) lsysbl_createEmptyMatrix
!#      -> Creates an empty block matrix
!#
!#  9.) lsysbl_duplicateMatrix
!#      -> Duplicates a block matrix by duplicating all sub-matrices
!#      -> Extendend version of lsysbl_copyMatrix
!#
!# 10.) lsysbl_duplicateVector
!#      -> Duplicates a block vector by duplicating all sub-vectors
!#      -> Extended version of lsysbl_copyVector
!#
!# 11.) lsysbl_enforceStructure
!#      -> Enforces the structure of a given block vector in another
!#         block vector
!#
!# 12.) lsysbl_enforceStructureDirect
!#      -> Enforces a subvector structure according to a length definition
!#         to a block vector
!#
!# 13.) lsysbl_enforceStructureDiscr
!#      -> Enforces a subvector structure according to a given discretisation
!#
!# 14.) lsysbl_assignDiscrIndirect
!#      -> Assign discretisation related information of one vector
!#         to another
!#
!# 15.) lsysbl_assignDiscrIndirectMat
!#      -> Assign discretisation related information of a matrix
!#         to a vector to make it compatible.
!#
!# 16.) lsysbl_updateMatStrucInfo
!#      -> Recalculate structural data of a block matrix from
!#         the submatrices
!#
!# 17.) lsysbl_releaseVector
!#      -> Release a block vector from memory
!#
!# 18.) lsysbl_releaseMatrix
!#      -> Releases a block matrix and all submatrices
!#
!# 19.) lsysbl_releaseMatrixRow
!#      -> Releases a row of submatrices from a block matrix
!#
!# 20.) lsysbl_releaseMatrixColumn
!#      -> Releases a column of submatrices from a block matrix
!#
!# 21.) lsysbl_blockMatVec
!#      -> Multiply a block matrix with a block vector
!#
!# 22.) lsysbl_copyVector
!#       -> Copy a block vector to another one
!#
!# 23.) lsysbl_copyMatrix
!#       -> Copy a block matrix to another one
!#
!# 24.) lsysbl_scaleVector
!#      -> Scale a block vector by a constant
!#
!# 25.) lsysbl_clearVector
!#      -> Clear a block vector, i.e. overwrites all entries with 0.0 or 
!#         with a defined value
!#
!# 26.) lsysbl_vectorLinearComb
!#      -> Linear combination of two block vectors
!#
!# 27.) lsysbl_scalarProduct
!#      -> Calculate a scalar product of two vectors
!#
!# 28.) lsysbl_setSortStrategy
!#      -> Assigns a sorting strategy/permutation to every subvector
!#
!# 29.) lsysbl_sortVectorInSitu
!#      -> Resort the entries of all subvectors according to an assigned
!#         sorting strategy
!#
!# 30.) lsysbl_isVectorCompatible
!#      -> Checks whether two vectors are compatible to each other
!#
!# 31.) lsysbl_isMatrixCompatible
!#      -> Checks whether a matrix and a vector are compatible to each other
!#
!# 32.) lsysbl_isMatrixSorted
!#      -> Checks if a block matrix is sorted
!#
!# 33.) lsysbl_isVectorSorted
!#      -> Checks if a block vector is sorted
!#
!# 34.) lsysbl_getbase_double
!#      -> Get a pointer to the double precision data array of the vector
!#
!# 35.) lsysbl_getbase_single
!#      -> Get a pointer to the single precision data array of the vector
!#
!# 36.) lsysbl_vectorNorm
!#      -> Calculates the norm of a vector. the vector is treated as one
!#         long data array.
!#
!# 37.) lsysbl_vectorNormBlock
!#      -> Calculates the norm of all subvectors in a given block vector.
!#
!# 38.) lsysbl_invertedDiagMatVec
!#      -> Multiply a vector with the inverse of the diagonal of a matrix
!#
!# 39.) lsysbl_swapVectors
!#      -> Swap two vectors
!#
!# 40.) lsysbl_deriveSubvector
!#      -> Derives a blockvector as a subset of another blockvector
!#
!# 41.) lsysbl_deriveSubmatrix / lsysbl_extractSubmatrix
!#      -> Extracts a submatrix from a block matrix
!#
!# 42.) lsysbl_isSubmatrixPresent
!#      -> Checks if a submatrix of a blockmatrix is present
!#
!# 43.) lsysbl_resizeVectorBlock
!#      -> Resize a block vector
!#
!# 44.) lsysbl_resizeVecBlockIndMat
!#      -> Resize a block vector according to a block matrix
!#
!# 45.) lsysbl_infoVector
!#      -> Outputs information about the vector (mostly used for debugging)
!#
!# 46.) lsysbl_infoMatrix
!#      -> Outputs information about the matrix (mostly used for debugging)
!#
!# 47.) lsysbl_clearMatrix
!#      -> Clears a matrix, i.e. overwrites all entries with 0.0 or 
!#         with a defined value
!#
!# 48.) lsysbl_convertVecFromScalar
!#      -> Converts a scalar vector to a 1-block vector, whereby the
!#         block vector takes over leadership of the data.
!#
!# 49.) lsysbl_convertMatFromScalar
!#      -> Converts a scalar matrix to a 1x1 block vector, whereby the
!#         block matrix takes over leadership of the data.
!#
!# 50.) lsysbl_insertSubmatrix
!#      -> Insert a block matrix into another block matrix
!#
!# 51.) lsysbl_assignDiscrDirectMat
!#      -> Assign a block discretisation to a matrix
!#
!# 52.) lsysbl_assignDiscrDirectVec
!#      -> Assign a block discretisation to a vector
!#
!# 53.) lsysbl_createFpdbObjectVec
!#      -> Creates an ObjectItem representing a block vector
!#
!# 54.) lsysbl_createFpdbObjectMat
!#      -> Creates an ObjectItem representing a block matrix
!#
!# 55.) lsysbl_restoreFpdbObjectVec
!#      -> Restores a block vector from an ObjectItem
!#
!# 56.) lsysbl_restoreFpdbObjectMat
!#      -> Restores a block matrix from an ObjectItem
!#
!# 57.) lsyssc_unshareMatrix
!#      -> Renders a matrix independent, resets the sharing state
!# 
!# 58.) lsyssc_unshareVector
!#      -> Renders a vector independent, resets the sharing state
!#
!# 59.) lsysbl_moveToSubmatrix
!#      -> Moves a matrix to a submatrix of a larger matrix
!#
!# 60.) lsysbl_allocEmptyMatrix
!#      -> Allocates memory for the entries of a matrix.
!#
!# 61.) lsysbl_reintegrateSubmatrix
!#      -> Reintegrates a submatrix into a larger block matrix which was
!#         extracted from that matrix by lsysbl_deriveSubmatrix.
!#
!# 62.) lsysbl_getVectorMagnitude
!#      -> Compute the vector magnitude.
!#
!# 63.) lsysbl_synchroniseSortVecVec
!#      -> Synchronises the sorting between a vector and another vector
!#
!# 64.) lsysbl_synchroniseSortMatVec
!#      -> Synchrionises the sorting of a vector according to the sorting
!#         of a matrix
!#
!# 65.) lsysbl_createscalarfromvec
!#      -> Create a scalar vector from a block vector
!# </purpose>
!##############################################################################

module linearsystemblock

  use fpersistence
  use fsystem
  use genoutput
  use storage
  use basicgeometry
  use spatialdiscretisation
  use linearsystemscalar
  use linearalgebra
  use dofmapping
  use discretebc
  use discretefbc
  use uuid
  
  implicit none
  
  private

!<constants>

!<constantblock description="Flags for the matrix specification bitfield">

  ! Standard matrix
  integer(I32), parameter, public :: LSYSBS_MSPEC_GENERAL           =        0
  
  ! Block matrix is a scalar (i.e. 1x1) matrix.
  integer(I32), parameter, public :: LSYSBS_MSPEC_SCALAR            =        1

  ! Block matrix is of saddle-point type:
  !  (A  B1)
  !  (B2 0 )
  integer(I32), parameter, public :: LSYSBS_MSPEC_SADDLEPOINT       =        2

  ! Block matrix is nearly of saddle-point type 
  !  (A  B1)
  !  (B2 C )
  ! with C~0 being a stabilisation matrix
  integer(I32), parameter, public :: LSYSBS_MSPEC_NEARLYSADDLEPOINT =        3

  ! The block matrix is a submatrix of another block matrix and not
  ! located at the diagonal of its parent. 
  ! This e.g. modifies the way, boundary conditions are implemented
  ! into a matrix; see in the boundary condition implementation
  ! routines for details.
  integer(I32), parameter, public :: LSYSBS_MSPEC_OFFDIAGSUBMATRIX  =        4

  ! The submatrices in the block matrix all share the same structure.
  ! Submatrices are allowed to be empty. 
  integer(I32), parameter, public :: LSYSBS_MSPEC_GROUPMATRIX       =        5

!</constantblock>

!</constants>

!<types>
  
!<typeblock>
  
  ! A block vector that can be used by block linear algebra routines.
  
  type t_vectorBlock

    ! Universally unique identifier
    type(t_uuid) :: ruuid
    
    ! Total number of equations in the vector
    integer                    :: NEQ = 0
    
    ! Handle identifying the vector entries or = ST_NOHANDLE if not
    ! allocated.
    integer                    :: h_Ddata = ST_NOHANDLE
    
    ! Start position of the vector data in the array identified by
    ! h_Ddata. Normally = 1. Can be set to > 1 if the vector is a subvector
    ! in a larger memory block allocated on the heap.
    integer                    :: iidxFirstEntry = 1

    ! Data type of the entries in the vector. Either ST_SINGLE or
    ! ST_DOUBLE. The subvectors are all of the same type.
    integer                    :: cdataType = ST_DOUBLE

    ! Is set to true, if the handle h_Ddata belongs to another vector,
    ! i.e. when this vector shares data with another vector.
    logical                    :: bisCopy   = .false.

    ! Number of blocks allocated in RvectorBlock
    integer                    :: nblocks = 0
    
    ! Pointer to a block discretisation structure that specifies
    ! the discretisation of all the subblocks in the vector.
    ! Points to NULL(), if there is no underlying block discretisation
    ! structure.
    type(t_blockDiscretisation), pointer :: p_rblockDiscr => null()

    ! A pointer to discretised boundary conditions for real boundary components.
    ! These boundary conditions allow to couple multiple equations in the 
    ! system and don't belong only to one single scalar component of the
    ! solution.
    ! If no system-wide boundary conditions are specified, p_rdiscreteBC
    ! can be set to NULL().
    type(t_discreteBC), pointer  :: p_rdiscreteBC => null()
    
    ! A pointer to discretised boundary conditions for fictitious boundary
    ! components.
    ! These boundary conditions allow to couple multiple equations in the 
    ! system and don't belong only to one single scalar component of the
    ! solution.
    ! If no system-wide boundary conditions are specified, p_rdiscreteBCfict
    ! can be set to NULL().
    type(t_discreteFBC), pointer  :: p_rdiscreteBCfict => null()
    
    ! A 1D array with scalar vectors for all the blocks.
    ! The handle identifier inside of these blocks are set to h_Ddata.
    ! The iidxFirstEntry index pointer in each subblock is set according
    ! to the position inside of the large vector.
    type(t_vectorScalar), dimension(:), pointer :: RvectorBlock => null()
    
  end type
  
  public :: t_vectorBlock
  
!</typeblock>

!<typeblock>
  
  ! A block matrix that can be used by scalar linear algebra routines.
  ! Note that a block matrix is always quadratic with ndiagBlocks
  ! block rows and ndiagBlock block columns!
  ! If necessary, the corresponding block rows/columns are filled
  ! up with zero matrices!
  
  type t_matrixBlock
    
    ! Universally unique identifier
    type(t_uuid) :: ruuid

    ! Total number of equations = rows in the whole matrix.
    ! This usually coincides with NCOLS. NCOLS > NEQ indicates, that 
    ! at least one block row in the block matrix is completely zero.
    integer                    :: NEQ         = 0

    ! Total number of columns in the whole matrix.
    ! This usually coincides with NEQ. NCOLS < NEQ indicates, that
    ! at least one block column in the block matrix is completely zero.
    integer                    :: NCOLS       = 0
    
    ! Number of blocks per row.
    integer                    :: nblocksPerRow = 0
    
    ! Number of blocks per column.
    integer                    :: nblocksPerCol = 0

    ! Matrix specification tag. This is a bitfield coming from an OR 
    ! combination of different LSYSBS_MSPEC_xxxx constants and specifies
    ! various details of the matrix. If it is =LSYSBS_MSPEC_GENERAL,
    ! the matrix is a usual matrix that needs no special handling.
    integer(I32) :: imatrixSpec = LSYSBS_MSPEC_GENERAL
    
    ! Pointer to a block discretisation structure that specifies
    ! the discretisation of the test functions in the subblocks in the matrix.
    ! Points to NULL(), if there is no underlying block discretisation
    ! structure.
    type(t_blockDiscretisation), pointer :: p_rblockDiscrTest => null()

    ! Pointer to a block discretisation structure that specifies
    ! the discretisation of the trial functions in the subblocks in the matrix.
    ! Points to NULL(), if there is no underlying block discretisation
    ! structure.
    type(t_blockDiscretisation), pointer :: p_rblockDiscrTrial => null()

    ! Flag: Trial and Test functions in all element distributions
    ! are the same. If FALSE, there's at least one element distribution
    ! with different trial and test functions.
    logical                          :: bidenticalTrialAndTest = .true.

    ! A pointer to discretised boundary conditions for real boundary 
    ! components. If no boundary conditions are specified, p_rdiscreteBC
    ! can be set to NULL().
    type(t_discreteBC), pointer  :: p_rdiscreteBC => null()
    
    ! A pointer to discretised boundary conditions for fictitious boundary
    ! components. If no fictitious boundary conditions are specified, 
    ! p_rdiscreteBCfict can be set to NULL().
    type(t_discreteFBC), pointer  :: p_rdiscreteBCfict => null()
    
    ! A 2D array with scalar matrices for all the blocks.
    ! A submatrix is assumed to be empty (zero-block) if the corresponding
    ! NEQ=NCOLS=0.
    type(t_matrixScalar), dimension(:,:), pointer :: RmatrixBlock => null()
    
  end type
  
  public :: t_matrixBlock
  
!</typeblock>

!</types>

  interface lsysbl_createVectorBlock
    module procedure lsysbl_createVecBlockDirect 
    module procedure lsysbl_createVecBlockDirectDims
    module procedure lsysbl_createVecBlockIndirect 
    module procedure lsysbl_createVecBlockByDiscr
  end interface

  public :: lsysbl_createVectorBlock
  public :: lsysbl_createVecBlockDirect
  public :: lsysbl_createVecBlockIndirect
  public :: lsysbl_createVecBlockByDiscr

  interface lsysbl_resizeVectorBlock
    module procedure lsysbl_resizeVecBlockDirect
    module procedure lsysbl_resizeVecBlockDirectDims
    module procedure lsysbl_resizeVecBlockIndirect
  end interface

  public :: lsysbl_resizeVectorBlock
  public :: lsysbl_resizeVecBlockIndMat
  public :: lsysbl_resizeVecBlockDirectDims
  public :: lsysbl_resizeVecBlockIndirect
  
  public :: lsysbl_createVecBlockIndMat
  public :: lsysbl_createMatFromScalar
  public :: lsysbl_createVecFromScalar
  public :: lsysbl_createMatBlockByDiscr
  public :: lsysbl_createEmptyMatrix
  public :: lsysbl_duplicateMatrix
  public :: lsysbl_duplicateVector
  public :: lsysbl_enforceStructure
  public :: lsysbl_enforceStructureDirect
  public :: lsysbl_enforceStructureDiscr
  public :: lsysbl_assignDiscrIndirect
  public :: lsysbl_assignDiscrIndirectMat
  public :: lsysbl_updateMatStrucInfo
  public :: lsysbl_releaseVector
  public :: lsysbl_releaseMatrix
  public :: lsysbl_releaseMatrixRow
  public :: lsysbl_releaseMatrixColumn
  public :: lsysbl_blockMatVec
  public :: lsysbl_copyVector
  public :: lsysbl_copyMatrix
  public :: lsysbl_scaleVector
  public :: lsysbl_clearVector
  public :: lsysbl_vectorLinearComb
  public :: lsysbl_scalarProduct
  public :: lsysbl_setSortStrategy
  public :: lsysbl_sortVectorInSitu
  public :: lsysbl_isVectorCompatible
  public :: lsysbl_isMatrixCompatible
  public :: lsysbl_isMatrixSorted
  public :: lsysbl_isVectorSorted
  public :: lsysbl_getbase_double
  public :: lsysbl_getbase_single
  public :: lsysbl_vectorNorm
  public :: lsysbl_vectorNormBlock
  public :: lsysbl_invertedDiagMatVec
  public :: lsysbl_swapVectors
  public :: lsysbl_deriveSubvector
  public :: lsysbl_deriveSubmatrix
  public :: lsysbl_extractSubmatrix
  public :: lsysbl_isSubmatrixPresent
  public :: lsysbl_infoVector
  public :: lsysbl_infoMatrix
  public :: lsysbl_clearMatrix
  public :: lsysbl_convertVecFromScalar
  public :: lsysbl_convertMatFromScalar
  public :: lsysbl_insertSubmatrix
  public :: lsysbl_assignDiscrDirectMat
  public :: lsysbl_assignDiscrDirectVec
  public :: lsysbl_createFpdbObjectVec
  public :: lsysbl_createFpdbObjectMat
  public :: lsysbl_restoreFpdbObjectVec
  public :: lsysbl_restoreFpdbObjectMat
  public :: lsyssc_unshareMatrix
  public :: lsyssc_unshareVector
  public :: lsysbl_moveToSubmatrix
  public :: lsysbl_allocEmptyMatrix
  public :: lsysbl_reintegrateSubmatrix
  public :: lsysbl_getVectorMagnitude
  public :: lsysbl_synchroniseSortVecVec
  public :: lsysbl_synchroniseSortMatVec
  public :: lsysbl_createscalarfromvec
  
contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_isVectorCompatible (rvector1,rvector2,bcompatible)
  
!<description>
  ! Checks whether two vectors are compatible to each other.
  ! Two vectors are compatible if
  ! - they have the same size,
  ! - they have the same structure of subvectors,
  ! - they have the same sorting strategy OR they are both unsorted
!</description>

!<input>
  ! The first vector
  type(t_vectorBlock), intent(IN) :: rvector1
  
  ! The second vector
  type(t_vectorBlock), intent(IN) :: rvector2
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether the vectors are compatible or not.
  ! If not given, an error will inform the user if the two vectors are
  ! not compatible and the program will halt.
  logical, intent(OUT), optional :: bcompatible
!</output>

!</subroutine>

  ! local variables
  integer :: i

  ! We assume that we are not compatible
  if (present(bcompatible)) bcompatible = .false.
  
  ! Vectors must have the same size and number of blocks
  if (rvector1%NEQ .ne. rvector2%NEQ) then
    if (present(bcompatible)) then
      bcompatible = .false.
      return
    else
      print *,'Vectors not compatible, different size!'
      call sys_halt()
    end if
  end if

  if (rvector1%nblocks .ne. rvector2%nblocks) then
    if (present(bcompatible)) then
      bcompatible = .false.
      return
    else
      print *,'Vectors not compatible, different block structure!'
      call sys_halt()
    end if
  end if
  
  ! All the subblocks must be the same size and must be sorted the same way.
  do i=1,rvector1%nblocks
    if (rvector1%RvectorBlock(i)%NEQ .ne. rvector2%RvectorBlock(i)%NEQ) then
      if (present(bcompatible)) then
        bcompatible = .false.
        return
      else
        print *,'Vectors not compatible, different block structure!'
        call sys_halt()
      end if
    end if
    
    ! isortStrategy < 0 means unsorted. Both unsorted is ok.
    
    if ((rvector1%RvectorBlock(i)%isortStrategy .gt. 0) .or. &
        (rvector2%RvectorBlock(i)%isortStrategy .gt. 0)) then

      if (rvector1%RvectorBlock(i)%isortStrategy .ne. &
          rvector2%RvectorBlock(i)%isortStrategy) then
        if (present(bcompatible)) then
          bcompatible = .false.
          return
        else
          print *,'Vectors not compatible, differently sorted!'
          call sys_halt()
        end if
      end if

      if (rvector1%RvectorBlock(i)%h_isortPermutation .ne. &
          rvector2%RvectorBlock(i)%h_isortPermutation) then
        if (present(bcompatible)) then
          bcompatible = .false.
          return
        else
          print *,'Vectors not compatible, differently sorted!'
          call sys_halt()
        end if
      end if
    end if
  end do

  ! Ok, they are compatible
  if (present(bcompatible)) bcompatible = .true.

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_isMatrixCompatible (rvector,rmatrix,btransposed,bcompatible)
  
!<description>
  ! Checks whether a vector and a matrix are compatible to each other.
  ! A vector and a matrix are compatible if
  ! - they have the same NEQ,
  ! - they have the same structure of subvectors and diagonal blocks,
  ! - they have the same sorting strategy OR they are both unsorted
  ! Currently disabled checks:
  ! - they have the same spatial discretisation,
  ! - they have the same boundary conditions.
!</description>

!<input>
  ! The vector
  type(t_vectorBlock), intent(IN) :: rvector
  
  ! The matrix
  type(t_matrixBlock), intent(IN) :: rmatrix

  ! Check for rvector being compatible from the left or from the right
  ! to the matrix rmatrix.
  ! =FALSE: Check whether matrix-vector product $A*x$ is possible
  ! =TRUE : Check whether matrix-vector product $x^T*A = A^T*x$ is possible
  logical, intent(IN)              :: btransposed
!</input>

!<output>
  ! OPTIONAL: If given, the flag will be set to TRUE or FALSE depending on
  ! whether vector and matrix are compatible or not.
  ! If not given, an error will inform the user if the vector/matrix are
  ! not compatible and the program will halt.
  logical, intent(OUT), optional :: bcompatible
!</output>

!</subroutine>

  ! local variables
  integer :: i,j,itmp
  !LOGICAL :: b1,b2,b3

  ! We assume that we are not compatible
  if (present(bcompatible)) bcompatible = .false.
  
  ! Vector/Matrix must have the same size and number of blocks and equations
  if(btransposed) then
    itmp = rmatrix%NEQ
  else
    itmp = rmatrix%NCOLS
  end if
  if (rvector%NEQ .ne. itmp) then
    if (present(bcompatible)) then
      bcompatible = .false.
      return
    else
      print *,'Vector/Matrix not compatible, different size!'
      call sys_halt()
    end if
  end if

  if(btransposed) then
    itmp = rmatrix%nblocksPerCol
  else
    itmp = rmatrix%nblocksPerRow
  end if

  if (rvector%nblocks .ne. itmp) then
    if (present(bcompatible)) then
      bcompatible = .false.
      return
    else
      print *,'Vector/Matrix not compatible, different block structure!'
      call sys_halt()
    end if
  end if
  
!  ! Vector and matrix must share the same BC's.
!  b1 = ASSOCIATED(rvector%p_rdiscreteBC)
!  b2 = ASSOCIATED(rmatrix%p_rdiscreteBC)
!  b3 = ASSOCIATED(rvector%p_rdiscreteBC,rmatrix%p_rdiscreteBC)
!  IF ((b1 .OR. b2) .AND. .NOT. (b1 .AND. b2 .AND. b3)) THEN
!    IF (PRESENT(bcompatible)) THEN
!      bcompatible = .FALSE.
!      RETURN
!    ELSE
!      PRINT *,'Vector/Matrix not compatible, different boundary conditions!'
!      CALL sys_halt()
!    END IF
!  END IF
!
!  b1 = ASSOCIATED(rvector%p_rdiscreteBCfict)
!  b2 = ASSOCIATED(rmatrix%p_rdiscreteBCfict)
!  b3 = ASSOCIATED(rvector%p_rdiscreteBCfict,rmatrix%p_rdiscreteBCfict)
!  IF ((b1 .OR. b2) .AND. .NOT. (b1 .AND. b2 .AND. b3)) THEN
!    IF (PRESENT(bcompatible)) THEN
!      bcompatible = .FALSE.
!      RETURN
!    ELSE
!      PRINT *,'Vector/Matrix not compatible, different fict. boundary conditions!'
!      CALL sys_halt()
!    END IF
!  END IF
  
  ! All the subblocks must be the same size and must be sorted the same way.
  ! Each subvector corresponds to one 'column' in the block matrix.
  !
  ! Loop through all columns (!) of the matrix
  !do j=1,rvector%nblocks
  do j = 1, rmatrix%nblocksPerRow
  
    ! Loop through all rows of the matrix
    !do i=1,rvector%nblocks
    do i = 1, rmatrix%nblocksPerCol
    
      if (rmatrix%RmatrixBlock(i,j)%NEQ .eq. 0) cycle

!        b1 = ASSOCIATED(rvector%RvectorBlock(i)%p_rspatialDiscretisation)
!        b2 = ASSOCIATED(rmatrix%RmatrixBlock(i,j)%p_rspatialDiscretisation)
!        b3 = ASSOCIATED(rvector%RvectorBlock(i)%p_rspatialDiscretisation, &
!                        rmatrix%RmatrixBlock(i,j)%p_rspatialDiscretisation)
!        IF ((b1 .OR. b2) .AND. .NOT. (b1 .AND. b2 .AND. b3)) THEN
!          IF (PRESENT(bcompatible)) THEN
!            bcompatible = .FALSE.
!            RETURN
!          ELSE
!            PRINT *,'Vector/Matrix not compatible, different discretisation!'
!            CALL sys_halt()
!          END IF
!        END IF

      if(btransposed) then
        itmp = i
        if (rvector%RvectorBlock(i)%NEQ .ne. rmatrix%RmatrixBlock(i,j)%NEQ) then
          if (present(bcompatible)) then
            bcompatible = .false.
            return
          else
            print *,'Vector/Matrix not compatible, different block structure!'
            call sys_halt()
          end if
        end if
      else ! not transposed
        itmp = j
        if (rvector%RvectorBlock(j)%NEQ .ne. rmatrix%RmatrixBlock(i,j)%NCOLS) then
          if (present(bcompatible)) then
            bcompatible = .false.
            return
          else
            print *,'Vector/Matrix not compatible, different block structure!'
            call sys_halt()
          end if
        end if
      end if

      ! isortStrategy < 0 means unsorted. Both unsorted is ok.
      

      if ((rvector%RvectorBlock(itmp)%isortStrategy .gt. 0) .or. &
          (rmatrix%RmatrixBlock(i,j)%isortStrategy .gt. 0)) then

        if (rvector%RvectorBlock(itmp)%isortStrategy .ne. &
            rmatrix%RmatrixBlock(i,j)%isortStrategy) then
          if (present(bcompatible)) then
            bcompatible = .false.
            return
          else
            print *,'Vector/Matrix not compatible, differently sorted!'
            call sys_halt()
          end if
        end if

        if (rvector%RvectorBlock(itmp)%h_isortPermutation .ne. &
            rmatrix%RmatrixBlock(i,j)%h_isortPermutation) then
          if (present(bcompatible)) then
            bcompatible = .false.
            return
          else
            print *,'Vector/Matrix not compatible, differently sorted!'
            call sys_halt()
          end if
        end if
      end if
    
    end do
  end do

  ! Ok, they are compatible
  if (present(bcompatible)) bcompatible = .true.

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_getbase_double (rvector,p_Ddata)
  
!<description>
  ! Returns a pointer to the double precision data array of the vector.
  ! An error is thrown if the vector is not double precision.
!</description>

!<input>
  ! The vector
  type(t_vectorBlock), intent(IN) :: rvector
!</input>

!<output>
  ! Pointer to the double precision data array of the vector.
  ! NULL() if the vector has no data array.
  real(DP), dimension(:), pointer :: p_Ddata
!</output>

!</subroutine>

  ! Do we have data at all?
 if ((rvector%NEQ .eq. 0) .or. (rvector%h_Ddata .eq. ST_NOHANDLE)) then
   nullify(p_Ddata)
   return
 end if

  ! Check that the vector is really double precision
  if (rvector%cdataType .ne. ST_DOUBLE) then
    print *,'lsysbl_getbase_double: Vector is of wrong precision!'
    call sys_halt()
  end if

  ! Get the data array
  call storage_getbase_double (rvector%h_Ddata,p_Ddata)

  ! Modify the starting address/length to get the real array.
  p_Ddata => p_Ddata(rvector%iidxFirstEntry:rvector%iidxFirstEntry+rvector%NEQ-1)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_getbase_single (rvector,p_Fdata)
  
!<description>
  ! Returns a pointer to the single precision data array of the vector.
  ! An error is thrown if the vector is not single precision.
!</description>

!<input>
  ! The vector
  type(t_vectorBlock), intent(IN) :: rvector
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
    print *,'lsysbl_getbase_single: Vector is of wrong precision!'
    call sys_halt()
  end if

  ! Get the data array
  call storage_getbase_single (rvector%h_Ddata,p_Fdata)
  
  ! Modify the starting address/length to get the real array.
  p_Fdata => p_Fdata(rvector%iidxFirstEntry:rvector%iidxFirstEntry+rvector%NEQ-1)
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_getbase_int (rvector,p_Idata)
  
!<description>
  ! Returns a pointer to the integer data array of the vector.
  ! An error is thrown if the vector is not integer.
!</description>

!<input>
  ! The vector
  type(t_vectorBlock), intent(IN) :: rvector
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
    print *,'lsysbl_getbase_int: Vector is of wrong precision!'
    call sys_halt()
  end if

  ! Get the data array
  call storage_getbase_int (rvector%h_Ddata,p_Idata)
  
  ! Modify the starting address/length to get the real array.
  p_Idata => p_Idata(rvector%iidxFirstEntry:rvector%iidxFirstEntry+rvector%NEQ-1)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_createVecBlockDirect (rx, Isize, bclear, cdataType)
  
!<description>
  ! Initialises the vector block structure rx. Isize is an array
  ! of integers containing the length of the individual blocks.
  ! Memory is allocated on the heap to hold vectors of the size
  ! according to Isize.
  !
  ! Remark: There is no block discretisation structure attached to the vector!
!</description>

!<input>
  ! An array with length-tags for the different blocks
  integer, dimension(:), intent(IN)         :: Isize
  
  ! Optional: If set to YES, the vector will be filled with zero initially.
  logical, intent(IN), optional             :: bclear
  
  ! OPTIONAL: Data type identifier for the entries in the vector. 
  ! Either ST_SINGLE or ST_DOUBLE. If not present, ST_DOUBLE is assumed.
  integer, intent(IN),optional              :: cdataType
!</input>

!<output>
  
  ! Destination structure. Memory is allocated for each of the blocks.
  type(t_vectorBlock),intent(OUT) :: rx
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,n
  integer :: cdata
  
  cdata = ST_DOUBLE
  if (present(cdataType)) cdata = cdataType
  
  ! rx is initialised by INTENT(OUT) with the most common data.
  ! What is missing is the data array.
  !
  ! Allocate one large vector holding all data.
  call storage_new ('lsysbl_createVecBlockDirect', 'Vector', sum(Isize), cdata, &
                      rx%h_Ddata, ST_NEWBLOCK_NOINIT)
  rx%cdataType = cdata
  
  ! Initialise the sub-blocks. Save a pointer to the starting address of
  ! each sub-block.
  ! Denote in the subvector that the handle belongs to us - not to
  ! the subvector.
  allocate(rx%RvectorBlock(size(Isize)))
  
  n=1
  do i = 1,size(Isize)
    if (Isize(i) .gt. 0) then
      rx%RvectorBlock(i)%NEQ = Isize(i)
      rx%RvectorBlock(i)%iidxFirstEntry = n
      rx%RvectorBlock(i)%h_Ddata = rx%h_Ddata
      rx%RvectorBlock(i)%cdataType = rx%cdataType
      rx%RvectorBlock(i)%bisCopy = .true.
      n = n+Isize(i)
    else
      rx%RvectorBlock(i)%NEQ = 0
      rx%RvectorBlock(i)%iidxFirstEntry = 0
      rx%RvectorBlock(i)%h_Ddata = ST_NOHANDLE
    end if
  end do
  
  rx%NEQ = n-1
  rx%nblocks = size(Isize)
  
  ! The data of the vector belongs to us (we created the handle), 
  ! not to somebody else.
  rx%bisCopy = .false.
  
  ! Warning: don't reformulate the following check into one IF command
  ! as this might give problems with some compilers!
  if (present(bclear)) then
    if (bclear) then
      call lsysbl_clearVector (rx)
    end if
  end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_createVecBlockDirectDims (rx, isize, iblocks, bclear, cdataType)
  
!<description>
  ! Initialises the vector block structure rx. Isize is an integer
  ! which describes the length of each block. Iblocks denotes the
  ! number of similar vector blocks.
  ! Memory is allocated on the heap to hold vectors of the size
  ! according to Isize.
  !
  ! Remark: There is no block discretisation structure attached to the vector!
!</description>

!<input>
  ! An integer describing the length of each block
  integer, intent(IN) :: isize

  ! An integer describing the number of vector blocks
  integer, intent(IN) :: iblocks
  
  ! Optional: If set to YES, the vector will be filled with zero initially.
  logical, intent(IN), optional             :: bclear
  
  ! OPTIONAL: Data type identifier for the entries in the vector. 
  ! Either ST_SINGLE or ST_DOUBLE. If not present, ST_DOUBLE is assumed.
  integer, intent(IN),optional              :: cdataType
!</input>

!<output>
  
  ! Destination structure. Memory is allocated for each of the blocks.
  type(t_vectorBlock),intent(OUT) :: rx
  
!</output>
  
!</subroutine>

  ! local variables
  integer :: i,n
  integer :: cdata
  
  cdata = ST_DOUBLE
  if (present(cdataType)) cdata = cdataType
  
  ! rx is initialised by INTENT(OUT) with the most common data.
  ! What is missing is the data array.
  !
  ! Allocate one large vector holding all data.
  call storage_new ('lsysbl_createVecBlockDirect', 'Vector', &
                      isize*iblocks, cdata, &
                      rx%h_Ddata, ST_NEWBLOCK_NOINIT)
  rx%cdataType = cdata
  
  ! Initialise the sub-blocks. Save a pointer to the starting address of
  ! each sub-block.
  ! Denote in the subvector that the handle belongs to us - not to
  ! the subvector.
  allocate(rx%RvectorBlock(iblocks))
  
  n=1
  do i = 1,iblocks
    if (isize .gt. 0) then
      rx%RvectorBlock(i)%NEQ = isize
      rx%RvectorBlock(i)%iidxFirstEntry = n
      rx%RvectorBlock(i)%h_Ddata = rx%h_Ddata
      rx%RvectorBlock(i)%cdataType = rx%cdataType
      rx%RvectorBlock(i)%bisCopy = .true.
      rx%RvectorBlock(i)%cdataType = cdata
      n = n+isize
    else
      rx%RvectorBlock(i)%NEQ = 0
      rx%RvectorBlock(i)%iidxFirstEntry = 0
      rx%RvectorBlock(i)%h_Ddata = ST_NOHANDLE
    end if
  end do
  
  rx%NEQ = n-1
  rx%nblocks = iblocks
  
  ! The data of the vector belongs to us (we created the handle), 
  ! not to somebody else.
  rx%bisCopy = .false.
  
  ! Warning: don't reformulate the following check into one IF command
  ! as this might give problems with some compilers!
  if (present(bclear)) then
    if (bclear) then
      call lsysbl_clearVector (rx)
    end if
  end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_createVecBlockIndirect (rtemplate,rx,bclear,cdataType)
  
!<description>
  ! Initialises the vector block structure rx. rtemplate is an
  ! existing block vector structure.
  ! Memory is allocated on the heap for rx, the different components
  ! of rx will have the same size as the corresponding ones in rtemplate.
!</description>
  
!<input>
  ! A template vector structure
  type(t_vectorBlock),intent(IN) :: rtemplate
  
  ! Optional: If set to YES, the vector will be filled with zero initially.
  ! Otherwise the content of rx is undefined.
  logical, intent(IN), optional             :: bclear
  
  ! OPTIONAL: Data type identifier for the entries in the vector. 
  ! Either ST_SINGLE or ST_DOUBLE. If not present, the new vector will
  ! receive the same data type as rtemplate.
  integer, intent(IN),optional              :: cdataType
!</input>

!<output>
  
  ! Destination structure. Memory is allocated for each of the blocks.
  type(t_vectorBlock),intent(OUT) :: rx
  
!</output>
  
!</subroutine>

  integer :: cdata,i
  integer :: n
  
  cdata = rtemplate%cdataType
  if (present(cdataType)) cdata = cdataType

  ! Copy the vector structure with all pointers.
  ! The handle identifying the vector data on the heap will be overwritten later.
  rx = rtemplate
  nullify(rx%RvectorBlock)
  
  ! Check if number of diagonal blocks is nonzero, otherwise exit
  if (rx%nblocks <= 0) return
  allocate(rx%RvectorBlock(rx%nblocks))
  rx%RvectorBlock = rtemplate%RvectorBlock
  
  ! Allocate one large new vector holding all data.
  call storage_new ('lsysbl_createVecBlockDirect', 'Vector', rtemplate%NEQ, cdata, &
                      rx%h_Ddata, ST_NEWBLOCK_NOINIT)
  rx%cdataType = cdata
  
  ! As we have a new handle, reset the starting index of the vector
  ! in the array back to 1. This is necessary, as "rx = rtemplate" copied
  ! the of index from rtemplate to rx.
  rx%iidxFirstEntry = 1
  
  where (rtemplate%RvectorBlock%h_Ddata .ne. ST_NOHANDLE)
    ! Initialise the data type
    rx%RvectorBlock%cdataType = cdata
                     
    ! Put the new handle to all subvectors
    rx%RvectorBlock(:)%h_Ddata = rx%h_Ddata
  
    ! Note in the subvectors, that the handle belongs to somewhere else... to us!
    rx%RvectorBlock%bisCopy = .true.
  end where
  
  ! Relocate the starting indices of the subvector.
  n = 1
  do i=1,rx%nblocks
    rx%RvectorBlock(i)%iidxFirstEntry = n
    n = n + rx%RvectorBlock(i)%NEQ
  end do
  
  ! Our handle belongs to us, so it's not a copy of another vector.
  rx%bisCopy = .false.

  ! Warning: don't reformulate the following check into one IF command
  ! as this might give problems with some compilers!
  if (present(bclear)) then
    if (bclear) then
      call lsysbl_clearVector (rx)
    end if
  end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_createVecBlockByDiscr (rblockDiscretisation,rx,bclear,&
                                           cdataType)
  
!<description>
  ! Initialises the vector block structure rx based on a block discretisation
  ! structure rblockDiscretisation. 
  !
  ! Memory is allocated on the heap for rx. The size of the subvectors in rx
  ! is calculated according to the number of DOF's indicated by the
  ! spatial discretisation structures in rblockDiscretisation.
!</description>
  
!<input>
  ! A block discretisation structure specifying the spatial discretisations
  ! for all the subblocks in rx.
  type(t_blockDiscretisation),intent(IN), target :: rblockDiscretisation
  
  ! Optional: If set to YES, the vector will be filled with zero initially.
  ! Otherwise the content of rx is undefined.
  logical, intent(IN), optional             :: bclear
  
  ! OPTIONAL: Data type identifier for the entries in the vector. 
  ! Either ST_SINGLE or ST_DOUBLE. If not present, ST_DOUBLE is used.
  integer, intent(IN),optional              :: cdataType
!</input>

!<output>
  ! Destination structure. Memory is allocated for each of the blocks.
  ! A pointer to rblockDiscretisation is saved to rx.
  type(t_vectorBlock),intent(OUT) :: rx
!</output>
  
!</subroutine>

  integer :: cdata,i
  integer, dimension(:), allocatable :: Isize
  
  cdata = ST_DOUBLE
  if (present(cdataType)) cdata = cdataType
  
  allocate(Isize(max(1,rblockDiscretisation%ncomponents)))
  
  ! Loop to the blocks in the block discretisation. Calculate size (#DOF's)
  ! of all the subblocks.
  Isize(1) = 0             ! Initialisation in case ncomponents=0
  do i=1,rblockDiscretisation%ncomponents
    Isize(i) = dof_igetNDofGlob(rblockDiscretisation%RspatialDiscr(i))
  end do
  
  ! Create a new vector with that block structure
  call lsysbl_createVecBlockDirect (rx, Isize(:), bclear, cdataType)
  
  ! Initialise further data of the block vector
  rx%p_rblockDiscr => rblockDiscretisation
  
  ! Initialise further data in the subblocks
  do i=1,rblockDiscretisation%ncomponents
    rx%RvectorBlock(i)%p_rspatialDiscr=> &
      rblockDiscretisation%RspatialDiscr(i)
  end do
  
  deallocate(Isize)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_createMatBlockByDiscr (rblockDiscretisationTrial,rmatrix,&
      rblockDiscretisationTest)
  
!<description>
  ! Initialises the matrix block structure rmatrix based on a block 
  ! discretisation structure rblockDiscretisation. 
  !
  ! The basic variables of the structure are initialised (number of blocks,
  ! pointer to the block discretisation). The submatrices of the matrix
  ! are not initialised; this must be done by the application afterwards.
!</description>
  
!<input>
  ! A block discretisation structure specifying the spatial discretisations
  ! for the trial functions in all the subblocks in rmatrix.
  ! If rblockDiscretisationTest is not specified, this also specifies
  ! the discretisation of the trial functions, this test and trial
  ! functions coincide.
  type(t_blockDiscretisation),intent(IN), target :: rblockDiscretisationTrial

  ! A block discretisation structure specifying the spatial discretisations
  ! for the test functions in all the subblocks in rmatrix.
  type(t_blockDiscretisation),intent(IN), target, optional :: rblockDiscretisationTest
!</input>

!<output>
  ! Destination structure. Memory is allocated for each of the blocks.
  ! A pointer to rblockDiscretisation is saved to rx.
  type(t_matrixBlock),intent(OUT) :: rmatrix
!</output>
  
!</subroutine>

  integer :: i
  integer :: NEQ,NCOLS
  integer, dimension(:), pointer :: p_Isize,p_Isize2
  
  allocate(p_Isize(max(1,rblockDiscretisationTrial%ncomponents))) 
  allocate(p_Isize2(max(1,rblockDiscretisationTrial%ncomponents))) 
  
  ! Loop to the blocks in the block discretisation. Calculate size (#DOF's)
  ! of all the subblocks.
  p_Isize(1) = 0             ! Initialisation in case ncomponents=0
  do i=1,rblockDiscretisationTrial%ncomponents
    p_Isize(i) = dof_igetNDofGlob(rblockDiscretisationTrial%RspatialDiscr(i))
  end do
  NCOLS = sum(p_Isize(:))
  if(present(rblockDiscretisationTest)) then
    p_Isize2(1) = 0
    do i = 1, rblockDiscretisationTest%ncomponents
      p_Isize(i) = dof_igetNDofGlob(rblockDiscretisationTest%RspatialDiscr(i))
    end do
    NEQ = sum(p_Isize(:))
  else
    NEQ = NCOLS
  end if
  
  ! The 'INTENT(OUT)' already initialised the structure with the most common
  ! values. What is still missing, we now initialise:
  
  ! Initialise a pointer to the block discretisation
  rmatrix%p_rblockDiscrTrial => rblockDiscretisationTrial
  if (present(rblockDiscretisationTest)) then
    rmatrix%p_rblockDiscrTest => rblockDiscretisationTest
    rmatrix%bidenticalTrialAndTest = .false.
  else
    rmatrix%p_rblockDiscrTest => rblockDiscretisationTrial
    rmatrix%bidenticalTrialAndTest = .true.
  end if
  
  ! Initialise NEQ/NCOLS by default for a quadrativ matrix
  rmatrix%NEQ = NEQ
  rmatrix%NCOLS = NCOLS
  rmatrix%nblocksPerRow = rblockDiscretisationTrial%ncomponents
  if(present(rblockDiscretisationTest)) then
    rmatrix%nblocksPerCol = rblockDiscretisationTest%ncomponents
  else
    rmatrix%nblocksPerCol = rmatrix%nblocksPerRow
  end if
  
  ! Check if number of diagonal blocks is nonzeri, otherwise exit
  if ((rmatrix%nblocksPerCol <= 0) .or. (rmatrix%nblocksPerCol <= 0)) return
  
  allocate(rmatrix%RmatrixBlock(rmatrix%nblocksPerCol,rmatrix%nblocksPerRow))
  
  if ((rmatrix%nblocksPerCol .eq. 1) .and. &
      (rmatrix%nblocksPerRow .eq. 1)) then
    rmatrix%imatrixSpec = LSYSBS_MSPEC_SCALAR
  end if
  
  deallocate(p_Isize2)
  deallocate(p_Isize)
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_createEmptyMatrix (rmatrix,nblocksPerCol,nblocksPerRow)
  
!<description>
  ! Creates a basic (nblocksPerCol x nblocksPerRow) block matrix. Reserves
  ! memory for all the submatrices but does not initialise any submatrix.
  !
  ! This routine can be used to create a totally empty matrix without any
  ! discretisation structure in the background.
  !
  ! The caller can manually fill in scalar matrices in rmatrix%RmatrixBlock
  ! as necessary. Afterwards, lsysbl_updateMatStrucInfo can be used to
  ! calculate the actual matrix dimensions (NEQ,NCOLS,...).
!</description>

!<input>
  ! Number of blocks per column in the matrix
  integer, intent(IN) :: nblocksPerCol
  
  ! OPTIONAL: Number of blocks per row in the matrix. If not given, 
  ! nblocksPerCol is used.
  integer, optional, intent(IN) :: nblocksPerRow
!</input>

!<output>
  ! Block matrix structure to be initialised.
  type(t_matrixBlock), intent(OUT) :: rmatrix
!</output>

!</subroutine>

  integer :: nbpr
  
    ! Check if number of blocks is nonzero, otherwise exit
    if (nblocksPerCol .le. 0) return
    if (present(nblocksPerRow)) then
      if(nblocksPerRow .le. 0) return
      nbpr = nblocksPerRow
    else
      nbpr = nblocksPerCol
    end if
    

    ! Allocate memory for the blocks, that's it.
    allocate(rmatrix%RmatrixBlock(nblocksPerCol,nbpr))
    rmatrix%nblocksPerCol = nblocksPerCol
    rmatrix%nblocksPerRow = nbpr

    if ((rmatrix%nblocksPerCol .eq. 1) .and. &
        (rmatrix%nblocksPerRow .eq. 1)) then
      rmatrix%imatrixSpec = LSYSBS_MSPEC_SCALAR
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_createVecBlockIndMat (rtemplateMat,rx, bclear, btransposed,cdataType)
  
!<description>
  ! Initialises the vector block structure rx. rtemplateMat is an
  ! existing block matrix structure. The vector rx will be created
  ! according to the size of the blocks of the submatrices.
  !
  ! Memory is allocated on the heap for rx. The different components
  ! of rx will have the same size as the corresponding column
  ! blocks in rtemplateMat.
  ! The sorting strategies of the subvectors are initialised
  ! with the sorting strategies of the column blocks of rtemplateMat.
!</description>
  
!<input>
  ! A template vector structure
  type(t_matrixBlock), intent(IN) :: rtemplateMat
  
  ! OPTIONAL: If set to TRUE, the vector will be filled with zero initially.
  ! Otherwise the content of rx is undefined.
  logical, intent(IN), optional   :: bclear
  
  ! OPTIONAL: If not specified or set to FALSE, the vector will be
  ! created as 'right' vector of the matrix (so matrix vector multiplication
  ! (A x) is possible).
  ! If set to TRUE, the vector will be a 'left' vector, i.e. matrix
  ! vector multiplication (A^T x) is possible.
  logical, intent(IN), optional :: btransposed
  
  ! OPTIONAL: Data type identifier for the entries in the vector. 
  ! Either ST_SINGLE or ST_DOUBLE. If not present, ST_DOUBLE is assumed.
  integer, intent(IN),optional :: cdataType
!</input>

!<output>
  ! Destination structure. Memory is allocated for each of the blocks.
  type(t_vectorBlock),intent(OUT) :: rx
!</output>
  
!</subroutine>

    ! local variables
    integer :: i,n,j,NVAR,nbpc,nbpr
    integer :: cdata
    logical :: btrans
    integer :: NEQ, NCOLS

    cdata = ST_DOUBLE
    if (present(cdataType)) cdata = cdataType
  
    btrans = .false.
    if (present(btransposed)) then
      btrans = btransposed
    end if
    
    if (btrans) then
      nbpc = rtemplateMat%nblocksPerRow
      nbpr = rtemplateMat%nblocksPerCol
      NEQ = rtemplateMat%NCOLS
      NCOLS = rtemplateMat%NEQ
    else
      nbpc = rtemplateMat%nblocksPerCol
      nbpr = rtemplateMat%nblocksPerRow
      NEQ = rtemplateMat%NEQ
      NCOLS = rtemplateMat%NCOLS
    end if
    
    rx%NEQ = NCOLS

    ! Allocate one large vector holding all data.
    call storage_new ('lsysbl_createVecBlockIndMat', 'Vector', &
                        NCOLS, cdata, rx%h_Ddata, ST_NEWBLOCK_NOINIT)
    
    ! Check if number of bocks is nonzero, otherwise exit
    if (nbpr <= 0) return

    ! Initialise the sub-blocks. Save a pointer to the starting address of
    ! each sub-block.
    allocate(rx%RvectorBlock(nbpr))
    
    n=1
    do i = 1, nbpr
      ! Search for the first matrix in column i of the block matrix -
      ! this will give us information about the vector. Note that the
      ! diagonal blocks does not necessarily have to exist!
      do j = 1, nbpc
      
        if (btrans) then
          NEQ = rtemplateMat%RmatrixBlock(i,j)%NCOLS
          NCOLS = rtemplateMat%RmatrixBlock(i,j)%NEQ
          NVAR = rtemplateMat%RmatrixBlock(i,j)%NVAR
        else
          NEQ = rtemplateMat%RmatrixBlock(j,i)%NEQ
          NCOLS = rtemplateMat%RmatrixBlock(j,i)%NCOLS
          NVAR = rtemplateMat%RmatrixBlock(j,i)%NVAR
        end if
      
        ! Check if the matrix is not empty
        if (NCOLS .gt. 0) then
          
          ! Found a template matrix we can use :-)
          rx%RvectorBlock(i)%NEQ  = NCOLS
          rx%RvectorBlock(i)%NVAR = NVAR

          ! Take the handle of the complete-solution vector, but set the index of
          ! the first entry to a value >= 1 - so that it points to the first
          ! entry in the global solution vector!
          rx%RvectorBlock(i)%h_Ddata = rx%h_Ddata
          rx%RvectorBlock(i)%iidxFirstEntry = n
          
          ! Give the vector the discretisation of the matrix
          if (.not. btrans) then
            rx%RvectorBlock(i)%p_rspatialDiscr => &
              rtemplateMat%RmatrixBlock(j,i)%p_rspatialDiscrTrial

            ! Give the vector the same sorting strategy as the matrix, so that
            ! the matrix and vector get compatible. Otherwise, things
            ! like matrix vector multiplication won't work...
            rx%RvectorBlock(i)%isortStrategy = &
              rtemplateMat%RmatrixBlock(j,i)%isortStrategy
            rx%RvectorBlock(i)%h_IsortPermutation = &
              rtemplateMat%RmatrixBlock(j,i)%h_IsortPermutation
          else
            rx%RvectorBlock(i)%p_rspatialDiscr => &
              rtemplateMat%RmatrixBlock(i,j)%p_rspatialDiscrTest

            ! Give the vector the same sorting strategy as the matrix, so that
            ! the matrix and vector get compatible. Otherwise, things
            ! like matrix vector multiplication won't work...
            rx%RvectorBlock(i)%isortStrategy = &
              rtemplateMat%RmatrixBlock(i,j)%isortStrategy
            rx%RvectorBlock(i)%h_IsortPermutation = &
              rtemplateMat%RmatrixBlock(i,j)%h_IsortPermutation
              
            ! DOES THIS WORK?!?!? I don't know =) MK
          end if
            
          ! Denote in the subvector that the handle belongs to us - not to
          ! the subvector.
          rx%RvectorBlock(i)%bisCopy = .true.
          
          ! Set the data type
          rx%RvectorBlock(i)%cdataType = cdata
          
          n = n+NCOLS
          
          ! Finish this loop, continue with the next column
          exit
          
        end if
        
      end do
      
      if (j .gt. nbpc) then
        ! Let's hope this situation (an empty equation) never occurs - 
        ! might produce some errors elsewhere :)
        rx%RvectorBlock(i)%NEQ  = 0
        rx%RvectorBlock(i)%NVAR = 1
        rx%RvectorBlock(i)%iidxFirstEntry = 0
      end if
      
    end do
    
    rx%nblocks = nbpr
    rx%cdataType = cdata

    ! Transfer the boundary conditions and block discretisation pointers
    ! from the matrix to the vector.
    if (.not. btrans) then
      rx%p_rblockDiscr => rtemplateMat%p_rblockDiscrTrial
    else
      rx%p_rblockDiscr => rtemplateMat%p_rblockDiscrTest
    end if
    rx%p_rdiscreteBC     => rtemplateMat%p_rdiscreteBC
    rx%p_rdiscreteBCfict => rtemplateMat%p_rdiscreteBCfict
    
    ! Warning: don't reformulate the following check into one IF command
    ! as this might give problems with some compilers!
    if (present(bclear)) then
      if (bclear) then
        call lsysbl_clearVector (rx)
      end if
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_assignDiscrIndirect (rtemplate,rx)
  
!<description>
  ! Assigns discretisation-related information (spatial discretisation, 
  ! boundary conditions,...) to rx. The vector rtemplate is used as a
  ! template, so at the end of the routine, rx and rtemplate share the
  ! same discretisation.
!</description>
  
!<input>
  ! A template vector structure
  type(t_vectorBlock),intent(IN) :: rtemplate
!</input>

!<inputoutput>
  ! Destination structure. Discretisation-related information of rtemplate
  ! is assigned to rx.
  type(t_vectorBlock),intent(INOUT) :: rx
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: i

    ! Simply modify all pointers of all subvectors, that's it.
    do i=1,rx%nblocks
      ! Spatial discretisation
      rx%RvectorBlock(i)%p_rspatialDiscr => &
        rtemplate%RvectorBlock(i)%p_rspatialDiscr
    end do
    
    ! Transfer the boundary conditions and block discretisation pointers
    ! from the template to the vector.
    rx%p_rblockDiscr => rtemplate%p_rblockDiscr
    rx%p_rdiscreteBC     => rtemplate%p_rdiscreteBC
    rx%p_rdiscreteBCfict => rtemplate%p_rdiscreteBCfict
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_assignDiscrIndirectMat (rtemplateMat,rx,btransposed)
  
!<description>
  ! Assigns discretisation-related information (spatial discretisation, 
  ! boundary conditions,...) of a matrix to rx. The matrix rtemplateMat is 
  ! used as a template, so at the end of the routine, rx and rtemplate 
  ! share the same discretisation and boundary conditions.
  ! (More precisely, the blocks in rx and the diagonal blocks of
  !  rtemplateMat!)
  !
  ! In case the vector is completely incompatible to the matrix
  ! (different number of blocks, different NEQ), an error is thrown.
!</description>
  
!<input>
  ! A template vector structure
  type(t_matrixBlock),intent(IN) :: rtemplateMat

  ! OPTIONAL: If not specified or set to FALSE, the vector will be
  ! created as 'right' vector of the matrix (so matrix vector multiplication
  ! (A x) is possible).
  ! If set to TRUE, the vector will be a 'left' vector, i.e. matrix
  ! vector multiplication (A^T x) is possible.
  logical, intent(IN), optional :: btransposed
  
!</input>

!<inputoutput>
  ! Destination structure. Discretisation-related information of rtemplate
  ! is assigned to rx.
  type(t_vectorBlock),intent(INOUT) :: rx
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: i,j
    integer :: NEQ, NCOLS,nbpc,nbpr
    logical :: btrans

    btrans = .false.
    if (present(btransposed)) then
      btrans = btransposed
    end if
    
    if (btrans) then
      nbpc = rtemplateMat%nblocksPerRow
      nbpr = rtemplateMat%nblocksPerCol
      NEQ = rtemplateMat%NCOLS
      NCOLS = rtemplateMat%NEQ
    else
      nbpc = rtemplateMat%nblocksPerCol
      nbpr = rtemplateMat%nblocksPerRow
      NEQ = rtemplateMat%NEQ
      NCOLS = rtemplateMat%NCOLS
    end if

    ! There must be at least some compatibility between the matrix
    ! and the vector!
    if ((NCOLS .ne. rx%NEQ) .or. (nbpr .ne. rx%nblocks)) then
      print *,'lsysbl_assignDiscrIndirectMat error: Matrix/Vector incompatible!'
      call sys_halt()
    end if

    ! Simply modify all pointers of all subvectors, that's it.
    if (.not. btrans) then
      do i = 1, nbpr
        do j = 1, nbpc
          ! Search for the first matrix in the 'column' of the block matrix and use
          ! its properties to initialise
          if (rtemplateMat%RmatrixBlock(j,i)%NEQ .ne. 0) then
            ! Spatial discretisation
            rx%RvectorBlock(i)%p_rspatialDiscr => &
              rtemplateMat%RmatrixBlock(j,i)%p_rspatialDiscrTrial
            exit
          end if
        end do
      end do
    else
      do i = 1, nbpr
        do j = 1, nbpc
          ! Search for the first matrix in the 'column' of the block matrix and use
          ! its properties to initialise
          if (rtemplateMat%RmatrixBlock(i,j)%NCOLS .ne. 0) then
            ! Spatial discretisation
            rx%RvectorBlock(i)%p_rspatialDiscr => &
              rtemplateMat%RmatrixBlock(i,j)%p_rspatialDiscrTest
            exit
          end if
        end do
      end do
    end if
    
    ! Transfer the boundary conditions and block discretisation pointers
    ! from the matrix to the vector.
    if (.not. btrans) then
      rx%p_rblockDiscr => rtemplateMat%p_rblockDiscrTrial
    else
      rx%p_rblockDiscr => rtemplateMat%p_rblockDiscrTrial
    end if
    rx%p_rdiscreteBC     => rtemplateMat%p_rdiscreteBC
    rx%p_rdiscreteBCfict => rtemplateMat%p_rdiscreteBCfict

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_assignDiscrDirectMat (rmatrix,rdiscrTrial,rdiscrTest)
  
!<description>
  ! Assigns given discretisation structures for trial/test spaces to a
  ! matrix. If necessary, corrects NEQ/NCOLS in the matrix according to
  ! the discretisation structure(s).
  !
  ! Usually used after a lsysbl_deriveSubmatrix to set the correct
  ! block discretisation structure.
!</description>
  
!<input>
  ! Discretisation structure for trial functions.
  type(t_blockDiscretisation), intent(IN), target :: rdiscrTrial

  ! OPTIONAL: Discretisation structure for test functions.
  ! If not specified, trial and test functions coincide.
  type(t_blockDiscretisation), intent(IN), target, optional :: rdiscrTest
!</input>

!<inputoutput>
  ! Destination matrix.
  type(t_matrixBlock),intent(INOUT) :: rmatrix
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: i,j,nactcols,nactrows

    ! Modify all pointers of all submatrices.
    do j=1,rmatrix%nblocksPerRow
      do i=1,rmatrix%nblocksPerCol
        if (lsysbl_isSubmatrixPresent(rmatrix,i,j,.true.)) then

          if (.not. present(rdiscrTest)) then
            call lsyssc_assignDiscrDirectMat (rmatrix%RmatrixBlock(i,j),&
                rdiscrTrial%RspatialDiscr(j),rdiscrTrial%RspatialDiscr(i))
          else
            call lsyssc_assignDiscrDirectMat (rmatrix%RmatrixBlock(i,j),&
                rdiscrTrial%RspatialDiscr(j),rdiscrTest%RspatialDiscr(i))
          end if

        end if
      end do
    end do
    
    ! Check the number of columns/rows.
    ! If the number of columns/rows in the matrix are smaller than indicated
    ! by dof_igetNDofGlobBlock, we have to update that information as it may
    ! stem from zero block rows/columns.
    ! If it's larger, there's something wrong.
    
    nactcols = dof_igetNDofGlobBlock(rdiscrTrial)
    
    if (rmatrix%NCOLS .gt. nactcols) then
      call output_line ('Discretisation invalid for the matrix!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_assignDiscrDirectMat')
      call sys_halt()
    else
      if (rmatrix%NCOLS .lt. nactcols) then
        rmatrix%NCOLS = nactcols
      end if
    end if

    ! Set the block discretisation of the block matrix
    rmatrix%p_rblockDiscrTrial => rdiscrTrial
    
    ! Depending on whether rdiscrTest is given set the discretisation structure
    ! of the test functions.
    rmatrix%bidenticalTrialAndTest = present(rdiscrTest)
    
    if (present(rdiscrTest)) then
      
      nactrows = dof_igetNDofGlobBlock(rdiscrTest)
      
      if (rmatrix%NEQ .gt. nactrows) then
        call output_line ('Discretisation invalid for the matrix!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_assignDiscrDirectMat')
        call sys_halt()
      else
        if (rmatrix%NEQ .lt. nactrows) then
          rmatrix%NEQ = nactrows
        end if
      end if

      ! Set the block discretisation of the block matrix
      rmatrix%p_rblockDiscrTest => rdiscrTest

    else
      ! Trial and test functions coincide
      
      nactrows = dof_igetNDofGlobBlock(rdiscrTrial)
      
      if (rmatrix%NEQ .gt. nactrows) then
        call output_line ('Discretisation invalid for the matrix!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_assignDiscrDirectMat')
        call sys_halt()
      else
        if (rmatrix%NEQ .lt. nactrows) then
          rmatrix%NEQ = nactrows
        end if
      end if

      ! Set the block discretisation of the block matrix
      rmatrix%p_rblockDiscrTest => rdiscrTrial
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_assignDiscrDirectVec (rvector,rblockdiscr)
  
!<description>
  ! Assigns given discretisation structures for trial/test spaces to a
  ! vector.
  !
  ! Usually used after a lsysbl_deriveSubmatrix to set the correct
  ! block discretisation structure.
!</description>
  
!<input>
  ! Discretisation structure for trial functions.
  type(t_blockDiscretisation), intent(IN), target :: rblockDiscr
!</input>

!<inputoutput>
  ! Destination vector.
  type(t_vectorBlock),intent(INOUT) :: rvector
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: j

    ! Modify all pointers of all submatrices.
    do j=1,rvector%nblocks

      if (rvector%RvectorBlock(j)%NEQ .ne. &
          dof_igetNDofGlob(rblockDiscr%RspatialDiscr(j))) then
        call output_line ('Discretisation invalid for the vector!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_assignDiscrDirectVec')
        call sys_halt()
      end if

      rvector%RvectorBlock(j)%p_rspatialDiscr => rblockDiscr%RspatialDiscr(j)
    end do
    
    if (rvector%NEQ .ne. dof_igetNDofGlobBlock(rblockDiscr)) then
      call output_line ('Discretisation invalid for the matrix!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_assignDiscrDirectMat')
      call sys_halt()
    end if

    ! Set the block discretisation of the block vector
    rvector%p_rblockDiscr => rblockDiscr
    
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_releaseVector (rx)
  
!<description>
  ! Releases the memory that is reserved by rx.
  ! Remark: The memory associated to the vector is only released if this
  !  vector structure is the owner of the data.
!</description>
  
!<inputoutput>
  
  ! Vector structure that is to be released.
  type(t_vectorBlock),intent(INOUT) :: rx
  
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i
  
  if (rx%h_Ddata .eq. ST_NOHANDLE) then
    print *,'lsysbl_releaseVector warning: releasing unused vector.'
  end if
  
  ! Release the data - if the handle is not a copy of another vector!
  if ((.not. rx%bisCopy) .and. (rx%h_Ddata .ne. ST_NOHANDLE)) then
    call storage_free (rx%h_Ddata)
  else
    rx%h_Ddata = ST_NOHANDLE
  end if
  
  ! Clean up the structure
  do i = 1,rx%nblocks
    call lsyssc_releaseVector (rx%RvectorBlock(i))
  end do
  rx%NEQ = 0
  rx%nblocks = 0
  rx%iidxFirstEntry = 1
  
  nullify(rx%p_rblockDiscr)
  nullify(rx%p_rdiscreteBC)
  nullify(rx%p_rdiscreteBCfict)
  if (associated(rx%RvectorBlock)) deallocate(rx%RvectorBlock)
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine lsysbl_releaseMatrix (rmatrix)
 
!<description>
  ! Releases the memory that is reserved by rmatrix (including all submatrices).
!</description>
  
!<inputoutput>
  ! Block matrix to be released
  type(t_matrixBlock), intent(INOUT)                :: rMatrix
!</inputoutput>

!</subroutine>
  
  ! local variables
  integer :: i,j
  
  ! loop through all the submatrices and release them
  do j = 1, rmatrix%nblocksPerRow
    do i = 1, rmatrix%nblocksPerCol
      
      ! Only release the matrix if there is one.
      if (rmatrix%RmatrixBlock(i,j)%NA .ne. 0) then
        call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(i,j))
      end if
      
    end do
  end do
  
  ! Clean up the other variables, finish
  rmatrix%NEQ = 0
  rmatrix%NCOLS = 0
  rmatrix%nblocksPerCol = 0
  rmatrix%nblocksPerRow = 0
  
  nullify(rmatrix%p_rblockDiscrTrial)
  nullify(rmatrix%p_rblockDiscrTest)
  rmatrix%bidenticalTrialAndTest = .true.
  nullify(rmatrix%p_rdiscreteBC)
  nullify(rmatrix%p_rdiscreteBCfict)
  if (associated(rmatrix%RmatrixBlock)) deallocate(rmatrix%RmatrixBlock)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine lsysbl_releaseMatrixRow (rmatrix,irow)
 
!<description>
  ! Releases all submatrices in row irow from the block matrix rmatrix.
!</description>
  
!<input>
  ! Matrix row which should be released.
  integer, intent(IN)                               :: irow
!</input>
  
!<inputoutput>
  ! Block matrix to be released (partially)
  type(t_matrixBlock), intent(INOUT)                :: rMatrix
!</inputoutput>

!</subroutine>
  
  ! local variables
  integer :: j
  
    ! loop through all the submatrices in row irow and release them
    do j = 1, rmatrix%nblocksPerRow
      ! Only release the matrix if there is one.
      if (rmatrix%RmatrixBlock(irow,j)%NA .ne. 0) then
        call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(irow,j))
      end if
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine lsysbl_releaseMatrixColumn (rmatrix,icolumn)
 
!<description>
  ! Releases all submatrices in column icolumn from the block matrix rmatrix.
!</description>
  
!<input>
  ! Matrix column which should be released.
  integer, intent(IN)                               :: icolumn
!</input>
  
!<inputoutput>
  ! Block matrix to be released (partially)
  type(t_matrixBlock), intent(INOUT)                :: rMatrix
!</inputoutput>

!</subroutine>
  
  ! local variables
  integer :: i
  
    ! loop through all the submatrices in row irow and release them
    do i = 1, rmatrix%nblocksPerCol
      ! Only release the matrix if there is one.
      if (rmatrix%RmatrixBlock(i,icolumn)%NA .ne. 0) then
        call lsyssc_releaseMatrix (rmatrix%RmatrixBlock(i,icolumn))
      end if
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine lsysbl_blockMatVec (rmatrix, rx, ry, cx, cy, btransposed)
 
!<description>
  ! Performs a matrix vector multiplicationwith a given scalar matrix:
  !    $$ Dy   =   cx * rMatrix * rx   +   cy * ry $$
  ! Vector and matrix must be compatible to each other (same size, sorting 
  ! strategy,...).
!</description>
  
!<input>
  
  ! Block matrix
  type(t_matrixBlock), intent(IN)                   :: rmatrix

  ! Block vector to multiply with the matrix.
  type(t_vectorBlock), intent(IN)                   :: rx
  
  ! Multiplicative factor for rx
  real(DP), intent(IN)                              :: cx

  ! Multiplicative factor for ry
  real(DP), intent(IN)                              :: cy
  
  ! OPTIONAL: Specifies whether a multiplication with the matrix itself
  ! (.false.) or with its transpose (.true.) should be performed. If not
  ! given, .false. is assumed.
  logical, optional, intent(IN)                     :: btransposed
!</input>

!<inputoutput>
  
  ! Additive vector. Receives the result of the matrix-vector multiplication
  type(t_vectorBlock), intent(INOUT)                :: ry
  
!</inputoutput>

!</subroutine>
  
  ! local variables
  integer :: i,j,mvok,nbpc,nbpr
  real(DP) :: cyact
  logical :: btrans = .false.
  
  ! DEBUG!!!
  real(DP), dimension(:), pointer :: p_Dx, p_Dy
  call lsysbl_getbase_double (rx,p_Dx)
  call lsysbl_getbase_double (ry,p_Dy)
  
    if(present(btransposed)) btrans = btransposed
    
    ! The vectors must be compatible to each other.
    !call lsysbl_isVectorCompatible (rx,ry)

    ! and compatible to the matrix
    call lsysbl_isMatrixCompatible (ry,rmatrix,.not. btrans)
    call lsysbl_isMatrixCompatible (rx,rmatrix,btrans)
    
    if(btrans) then
      nbpc = rmatrix%nblocksPerRow
      nbpr = rmatrix%nblocksPerCol
    else
      nbpc = rmatrix%nblocksPerCol
      nbpr = rmatrix%nblocksPerRow
    end if
    
    ! loop through all the sub matrices and multiply.
    do i = 1, nbpc
    
      ! The original vector ry must only be respected once in the first
      ! matrix-vector multiplication. The additive contributions of
      ! the other blocks must simply be added to ry without ry being multiplied
      ! by cy again!
    
      cyact = cy
      MVOK = NO
      do j = 1, nbpr
        
        ! Only call the MV when there is a scalar matrix that we can use!
        if(.not. btrans) then
          if (lsysbl_isSubmatrixPresent (rmatrix,i,j)) then
            call lsyssc_scalarMatVec (rMatrix%RmatrixBlock(i,j), rx%RvectorBlock(j), &
                                      ry%RvectorBlock(i), cx, cyact)
            cyact = 1.0_DP
            mvok = YES
          end if
        else
          if (lsysbl_isSubmatrixPresent (rmatrix,j,i)) then
            call lsyssc_scalarMatVec (rMatrix%RmatrixBlock(j,i), rx%RvectorBlock(j), &
                                      ry%RvectorBlock(i), cx, cyact)
            cyact = 1.0_DP
            mvok = YES
          end if
        end if
        
      end do
      
      ! If mvok=NO here, there was no matrix in the current line!
      ! A very special case, but nevertheless may happen. In this case,
      ! simply scale the vector ry by cyact!
      
      if (mvok .eq. NO) then
        !call lsysbl_scaleVector (ry,cy)
        call lsyssc_scaleVector(ry%RvectorBlock(i),cy)
      end if
      
    end do
   
  end subroutine
  
  !****************************************************************************

!<subroutine>

  subroutine lsysbl_invertedDiagMatVec (rmatrix,rvectorSrc,dscale,rvectorDst)
  
!<description>
  ! This routine multiplies the weighted inverted diagonal $domega*D^{-1}$
  ! of the diagonal blocks in the matrix rmatrix with the vector rvectorSrc and 
  ! stores the result into the vector rvectorDst:
  !   $$rvectorDst_i = dscale * D_i^{-1} * rvectorSrc_i  , i=1..nblocks$$
  ! Both, rvectorSrc and rvectorDst may coincide.
!</description>
  
!<input>
  ! The matrix. 
  type(t_matrixBlock), intent(IN) :: rmatrix

  ! The source vector.
  type(t_vectorBlock), intent(IN) :: rvectorSrc

  ! A multiplication factor. Standard value is 1.0_DP
  real(DP), intent(IN) :: dscale
!</input>

!<inputoutput>
  ! The destination vector which receives the result.
  type(t_vectorBlock), intent(INOUT) :: rvectorDst
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: iblock
  
    ! Vectors and matrix must be compatible
    call lsysbl_isVectorCompatible (rvectorSrc,rvectorDst)
    call lsysbl_isMatrixCompatible (rvectorSrc,rmatrix,.false.)
    call lsysbl_isMatrixCompatible (rvectorDst,rmatrix,.true.)
    
    ! As both vectors are compatible to each other and compatible to the
    ! matrix, the matrix must be a square matrix!
  
    ! Loop over the blocks
    do iblock = 1,rvectorSrc%nblocks
      ! Multiply with the inverted diagonal of the submatrix.
      call lsyssc_invertedDiagMatVec (rmatrix%RmatrixBlock(iblock,iblock),&
            rvectorSrc%RvectorBlock(iblock),dscale,&
            rvectorDst%RvectorBlock(iblock))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_copyVector (rx,ry)
  
!<description>
  ! Copies vector data: ry = rx.
  ! If the destination vector is empty, a new vector is created and all
  ! data of rx is copied into that.
  ! If the destination vector ry exists in memory, it must be at least
  ! as large as rx. All structural data as well as the content of rx is  
  ! transferred to ry, so rx and ry are compatible to each other afterwards.
!</description>

!<input>
  ! Source vector
  type(t_vectorBlock),intent(IN) :: rx
!</input>

!<inputoutput>
  ! Destination vector
  type(t_vectorBlock),intent(INOUT) :: ry
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: h_Ddata, cdataType
  integer :: isize,NEQ,i
  integer :: ioffset
  logical :: bisCopy
  type(t_vectorScalar), dimension(:), pointer :: p_rblocks
  
  ! If the destination vector does not exist, create a new one
  ! based on rx.
  if (ry%h_Ddata .eq. ST_NOHANDLE) then
    call lsysbl_createVecBlockIndirect (rx,ry,.false.)
  end if
  
  call storage_getsize (ry%h_Ddata, isize)
  NEQ = rx%NEQ
  if (isize .lt. NEQ) then
    print *,'lsysbl_copyVector: Destination vector too small!'
    call sys_halt()
  end if
  
  if (rx%cdataType .ne. ry%cdataType) then
    print *,'lsysbl_copyVector: Destination vector has different type &
            &than source vector!'
    call sys_halt()
  end if
  
  ! First, make a backup of some crucial data so that it does not
  ! get destroyed.
  h_Ddata = ry%h_Ddata
  cdataType = ry%cdataType
  bisCopy = ry%bisCopy
  ioffset = ry%iidxFirstEntry
  p_rblocks => ry%RvectorBlock
  
  ! Then transfer all structural information of rx to ry.
  ! This automatically makes both vectors compatible to each other.
  ry = rx
  
  ! Restore crucial data
  ry%h_Ddata = h_Ddata
  ry%bisCopy = bisCopy
  ry%cdataType = cdataType
  ry%iidxFirstEntry = ioffset
  ry%RvectorBlock => p_rblocks
  
  ! If necessary, allocate new memory for the blocks.
  if (.not. associated(ry%RvectorBlock)) then
    allocate(ry%RvectorBlock(ry%nblocks))
  end if
  
  ! Copy the block structure. Don't destroy crucial data.
  do i=1,size(ry%RvectorBlock)
    h_Ddata = ry%RvectorBlock(i)%h_Ddata
    cdataType = ry%RvectorBlock(i)%cdataType
    bisCopy = ry%RvectorBlock(i)%bisCopy
    ioffset = ry%RvectorBlock(i)%iidxFirstEntry
    
    ! Transfer all structural information of rx to ry.
    ! This automatically makes both vectors compatible to each other.
    ry%RvectorBlock(i) = rx%RvectorBlock(i)
    
    ! Restore crucial data
    ry%RvectorBlock(i)%h_Ddata = h_Ddata
    ry%RvectorBlock(i)%bisCopy = bisCopy
    ry%RvectorBlock(i)%cdataType = cdataType
    ry%RvectorBlock(i)%iidxFirstEntry = ioffset 
  end do
  
  ! And finally copy the data. 
  call lsysbl_copyVectorDirect (rx,ry)
   
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_copyVectorDirect (rx,ry)
  
!<description>
  ! Copies vector data: ry = rx.
  ! Bypasses all structure checks. Directly copies data from rx to ry.
  !
  ! THIS ROUTINE IS USUALLY ONLY USED INTERNALLY! USE THIS ROUTINE WITH CARE!
  ! USE IT ONLY IF YOU NEED SPEED, KNOWING WHAT YOU ARE DOING!
  !
  ! Under normal circumstances, use lsysbl_copyVector or lsysbl_duplicateVector!
!</description>

!<input>
  ! Source vector
  type(t_vectorBlock),intent(IN) :: rx
!</input>

!<inputoutput>
  ! Destination vector
  type(t_vectorBlock),intent(INOUT) :: ry
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: NEQ
  real(DP), dimension(:), pointer :: p_Dsource,p_Ddest
  real(SP), dimension(:), pointer :: p_Fsource,p_Fdest
  
  ! Copy the data
  NEQ = ry%NEQ
  select case (rx%cdataType)
  case (ST_DOUBLE)
    select case (ry%cdataType)
    case (ST_DOUBLE)
      call lsysbl_getbase_double (rx,p_Dsource)
      call lsysbl_getbase_double (ry,p_Ddest)
      call lalg_copyVectorDble (p_Dsource(1:NEQ),p_Ddest(1:NEQ))
      
    case (ST_SINGLE)
      call lsysbl_getbase_double (rx,p_Dsource)
      call lsysbl_getbase_single (ry,p_Fdest)
      call lalg_copyVectorDblSngl (p_Dsource(1:NEQ),p_Fdest(1:NEQ))
      
    case DEFAULT
      call output_line ('Unsupported data type!', &
                         OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_copyVectorDirect')
      call sys_halt()
    end select
    
  case (ST_SINGLE)
    select case (ry%cdataType)
    case (ST_DOUBLE)
      call lsysbl_getbase_single (rx,p_Fsource)
      call lsysbl_getbase_double (ry,p_Ddest)
      call lalg_copyVectorSnglDbl (p_Fsource(1:NEQ),p_Ddest(1:NEQ))
      
    case (ST_SINGLE)
      call lsysbl_getbase_single (rx,p_Fsource)
      call lsysbl_getbase_single (ry,p_Fdest)
      call lalg_copyVectorSngl (p_Fsource(1:NEQ),p_Fdest(1:NEQ))

    case DEFAULT
      call output_line ('Unsupported data type!', &
                         OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_copyVectorDirect')
      call sys_halt()
    end select

  case DEFAULT
    call output_line ('Unsupported data type!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_copyVectorDirect')
    call sys_halt()
  end select
   
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_copyMatrix (rsourceMatrix,rdestMatrix)
  
!<description>
  ! Copies a matrix data: rsourceMatrix = rdestMatrix.
  ! All structural (discretisation related) data of rsourceMatrix is
  ! transferred to rdestMatrix.
  !
  ! If the destination matrix contains data, the data is overwritten,
  ! regardless of whether the data belongs to rdestmatrix or not.
  ! If the destination matrix is empty, new memory is allocated.
  !
  ! Therefore, lsysbl_copyMatrix coincides with a call to
  ! lsysbl_duplicateMatrix (.,.,LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE).
!</description>
  
!<input>
  ! Source vector
  type(t_matrixBlock),intent(IN) :: rsourceMatrix
!</input>

!<inputoutput>
  ! Destination vector
  type(t_matrixBlock),intent(INOUT) :: rdestMatrix
!</inputoutput>
  
!</subroutine>

    call lsysbl_duplicateMatrix (rsourceMatrix,rdestMatrix,&
        LSYSSC_DUP_COPYOVERWRITE,LSYSSC_DUP_COPYOVERWRITE)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_scaleVector (rx,c)
  
!<description>
  ! Scales a vector vector rx: rx = c * rx
!</description>
  
!<inputoutput>
  ! Source and destination vector
  type(t_vectorBlock), intent(INOUT) :: rx
!</inputoutput>

!<input>
  ! Multiplication factor
  real(DP), intent(IN) :: c
!</input>
  
!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata
  
  ! Taje care of the data type!
  select case (rx%cdataType)
  case (ST_DOUBLE)
    ! Get the pointer and scale the whole data array.
    call lsysbl_getbase_double(rx,p_Ddata)
    call lalg_scaleVectorDble (p_Ddata,c)  

  case (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    call lsysbl_getbase_single(rx,p_Fdata)
    call lalg_scaleVectorSngl (p_Fdata,real(c,SP))  

  case DEFAULT
    print *,'lsysbl_scaleVector: Unsupported data type!'
    call sys_halt()
  end select
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_clearVector (rx,dvalue)
  
!<description>
  ! Clears the block vector dx: Dx = 0 (or Dx = dvalue if dvalue is specified)
!</description>

!<inputoutput>
  ! Destination vector to be cleared
  type(t_vectorBlock), intent(INOUT) :: rx

  ! OPTIONAL: Value to write into the matrix.
  ! If not specified, all matrix entries are set to 0.0.
  ! If specified, all matrix entries are set to dvalue.
  real(DP), intent(IN), optional :: dvalue
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_Dsource
  real(SP), dimension(:), pointer :: p_Ssource
  
  ! Take care of the data type
  select case (rx%cdataType)
  case (ST_DOUBLE)
    ! Get the pointer and scale the whole data array.
    call lsysbl_getbase_double(rx,p_Dsource)
    if (.not. present(dvalue)) then
      call lalg_clearVectorDble (p_Dsource)
    else
      call lalg_setVectorDble (p_Dsource,dvalue)
    end if

  case (ST_SINGLE)
    ! Get the pointer and scale the whole data array.
    call lsysbl_getbase_single(rx,p_Ssource)
    if (.not. present(dvalue)) then
      call lalg_clearVectorSngl (p_Ssource)
    else
      call lalg_setVectorSngl (p_Ssource,real(dvalue,SP))
    end if

  case DEFAULT
    print *,'lsysbl_clearVector: Unsupported data type!'
    call sys_halt()
  end select
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_vectorLinearComb (rx,ry,cx,cy,rdest)
  
!<description>
  ! Performs a linear combination: ry = cx * rx  +  cy * ry
  ! If rdest is given, the routine calculates rdest = cx * rx  +  cy * ry
  ! without overwriting ry.
  ! All vectors must be compatible to each other (same size, sorting 
  ! strategy,...).
!</description>  
  
!<input>
  ! First source vector
  type(t_vectorBlock), intent(IN)    :: rx
  
  ! Scaling factor for Dx
  real(DP), intent(IN)               :: cx

  ! Scaling factor for Dy
  real(DP), intent(IN)               :: cy
!</input>

!<inputoutput>
  ! Second source vector; also receives the result if rdest is not given.
  type(t_vectorBlock), intent(INOUT) :: ry
  
  ! OPTIONAL: Destination vector. If not given, ry is overwritten.
  type(t_vectorBlock), intent(INOUT), optional :: rdest
!</inputoutput>
  
!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_Dsource, p_Ddest, p_Ddest2
  real(SP), dimension(:), pointer :: p_Ssource, p_Sdest, p_Sdest2
  
  ! The vectors must be compatible to each other.
  call lsysbl_isVectorCompatible (rx,ry)
  
  if (present(rdest)) then
    call lsysbl_isVectorCompatible (rx,rdest)
  end if

  if (rx%cdataType .ne. ry%cdataType) then
    print *,'lsysbl_vectorLinearComb: different data types not supported!'
    call sys_halt()
  end if
  
  select case (rx%cdataType)
  case (ST_DOUBLE)
    ! Get the pointers and copy the whole data array.
    call lsysbl_getbase_double(rx,p_Dsource)
    if (.not. present(rdest)) then
      call lsysbl_getbase_double(ry,p_Ddest)
    else
      call lsysbl_getbase_double(rdest,p_Ddest)
      call lsysbl_getbase_double(ry,p_Ddest2)
      call lalg_copyVectorDble (p_Ddest2,p_Ddest)
    end if
    
    call lalg_vectorLinearCombDble (p_Dsource,p_Ddest,cx,cy)

  case (ST_SINGLE)
    ! Get the pointers and copy the whole data array.
    call lsysbl_getbase_single(rx,p_Ssource)
    if (.not. present(rdest)) then
      call lsysbl_getbase_single(ry,p_Sdest)
    else
      call lsysbl_getbase_single(rdest,p_Sdest)
      call lsysbl_getbase_single(ry,p_Sdest2)
      call lalg_copyVectorSngl (p_Sdest2,p_Sdest)
    end if
    
    call lalg_vectorLinearCombSngl (p_Ssource,p_Sdest,real(cx,SP),real(cy,SP))
  
  case DEFAULT
    print *,'lsysbl_vectorLinearComb: Unsupported data type!'
    call sys_halt()
  end select

  end subroutine
  
  !****************************************************************************
!<function>
  
  real(DP) function lsysbl_scalarProduct (rx, ry)
  
!<description>
  ! Calculates a scalar product of two block vectors.
  ! Both vectors must be compatible to each other (same size, sorting 
  ! strategy,...).
!</description>
  
!<input>
  ! First vector
  type(t_vectorBlock), intent(IN)                  :: rx

  ! Second vector
  type(t_vectorBlock), intent(IN)                  :: ry

!</input>

!<result>
  ! The scalar product (rx,ry) of the two block vectors.
!</result>

!</function>

  ! local variables
  real(DP), dimension(:), pointer :: p_Ddata1dp
  real(DP), dimension(:), pointer :: p_Ddata2dp
  real(SP), dimension(:), pointer :: p_Fdata1dp
  real(SP), dimension(:), pointer :: p_Fdata2dp
  ! integer :: i
  real(DP) :: res

  ! The vectors must be compatible to each other.
  call lsysbl_isVectorCompatible (rx,ry)

  ! Is there data at all?
  res = 0.0_DP
  
  if ( (rx%NEQ .eq. 0) .or. (ry%NEQ .eq. 0) .or. (rx%NEQ .ne. rx%NEQ)) then
    print *,'Error in lsysbl_scalarProduct: Vector dimensions wrong!'
    call sys_halt()
  end if

  if (rx%cdataType .ne. ry%cdataType) then
    print *,'lsysbl_scalarProduct: Data types different!'
    call sys_halt()
  end if

  ! Take care of the data type before doing a scalar product!
  select case (rx%cdataType)
  case (ST_DOUBLE)

    ! Get the data arrays
    call lsysbl_getbase_double (rx,p_Ddata1dp)
    call lsysbl_getbase_double (ry,p_Ddata2dp)
    
    ! Perform the scalar product
    res=lalg_scalarProduct(p_Ddata1dp,p_Ddata2dp)
!!$      res = 0.0_DP
!!$      DO i=1,rx%NEQ
!!$        res = res + p_Ddata1dp(i)*p_Ddata2dp(i)
!!$      END DO
    
  case (ST_SINGLE)

    ! Get the data arrays
    call lsysbl_getbase_single (rx,p_Fdata1dp)
    call lsysbl_getbase_single (ry,p_Fdata2dp)

    ! Perform the scalar product
    res=lalg_scalarProduct(p_Fdata1dp,p_Fdata2dp)

  case default
    print *,'lsysbl_scalarProduct: Not supported precision combination'
    call sys_halt()
  end select
  
  ! Return the scalar product, finish
  lsysbl_scalarProduct = res

  end function

  !****************************************************************************
!<subroutine>
  
  real(DP) function lsysbl_vectorNorm (rx,cnorm,iposMax)
  
!<description>
  ! Calculates the norm of a vector. cnorm specifies the type of norm to
  ! calculate.
!</description>
  
!<input>
  ! Vector to calculate the norm of.
  type(t_vectorBlock), intent(IN)                  :: rx

  ! Identifier for the norm to calculate. One of the LINALG_NORMxxxx constants.
  integer, intent(IN) :: cnorm
!</input>

!<output>
  ! OPTIONAL: If the MAX norm is to calculate, this returns the
  ! position of the largest element. If another norm is to be
  ! calculated, the result is undefined.
  integer, intent(OUT), optional :: iposMax
!</output>

!<result>
  ! The norm of the given array.
  ! < 0, if an error occurred (unknown norm).
!</result>

!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_Ddata
  real(SP), dimension(:), pointer :: p_Fdata

  ! Is there data at all?
  if (rx%h_Ddata .eq. ST_NOHANDLE) then
    print *,'Error in lsysbl_vectorNorm: Vector empty!'
    call sys_halt()
  end if
  
  ! Take care of the data type before doing a scalar product!
  select case (rx%cdataType)
  case (ST_DOUBLE)
    ! Get the array and calculate the norm
    call lsysbl_getbase_double (rx,p_Ddata)
    lsysbl_vectorNorm = lalg_norm (p_Ddata,cnorm,iposMax) 
    
  case (ST_SINGLE)
    ! Get the array and calculate the norm
    call lsysbl_getbase_single (rx,p_Fdata)
    lsysbl_vectorNorm = lalg_norm (p_Fdata,cnorm,iposMax) 
    
  case DEFAULT
    print *,'lsysbl_vectorNorm: Unsupported data type!'
    call sys_halt()
  end select
  
  end function

  !****************************************************************************
!<subroutine>
  
  subroutine lsysbl_vectorNormBlock (rx,Cnorms,Dnorms,IposMax) 
  
!<description>
  ! Calculates the norms of all subvectors in a given block vector.
  ! Cnorms is an array containing the type of norm to compute for each
  ! subvector.
  !
  ! If the arrays Cnorms or IposMax are smaller than rx%nblocks,
  ! only the first couple of subvectors are processed until the
  ! first array (Cnorms or IposMax) is completely filled.
!</description>
  
!<input>
  ! Vector to calculate the norm of.
  type(t_vectorBlock), intent(IN)                   :: rx

  ! Identifier list. For every subvector in rx, this identifies the norm 
  ! to calculate. Each entry is a LINALG_NORMxxxx constants.
  integer, dimension(:), intent(IN)                 :: Cnorms
!</input>

!<output>
  ! OPTIONAL: For each subvector: if the MAX norm is to calculate, 
  ! this returns the position of the largest element in that subvector. 
  ! If another norm is to be calculated, the result is undefined.
  integer, dimension(:), intent(OUT), optional :: IposMax
!</output>

!<result>
  ! An array of norms for each subvector.
  ! An entry might be < 0, if an error occurred (unknown norm).
  real(DP), dimension(:), intent(out) :: Dnorms
!</result>

!</subroutine>

  ! local variables
  integer :: i,nblk

  nblk = min(rx%nblocks,size(Cnorms))
  if (present(IposMax)) nblk = min(nblk,size(IposMax))

  ! Loop over the subvectors. Don't calculate more subvectors as we 
  ! are allowed to.
  do i=1,nblk
    ! Calculate the norm of that subvector.
    if (present(IposMax)) then
      Dnorms(i) = lsyssc_vectorNorm (rx%RvectorBlock(i),Cnorms(i),IposMax(i))
    else
      Dnorms(i) = lsyssc_vectorNorm (rx%RvectorBlock(i),Cnorms(i))
    end if
  end do
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_createMatFromScalar (rscalarMat,rmatrix,&
                                         rblockDiscrTrial,rblockDiscrTest)
  
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
  type(t_matrixScalar), intent(IN) :: rscalarMat
  
  ! OPTIONAL: A block discretisation structure specifying the trial functions.
  ! A pointer to this will be saved to the matrix.
  type(t_blockDiscretisation), intent(IN), optional, target :: rblockDiscrTrial

  ! OPTIONAL: A block discretisation structure specifying the test functions.
  ! A pointer to this will be saved to the matrix.
  ! If not specified, the trial functions coincide with the test functions.
  type(t_blockDiscretisation), intent(IN), optional, target :: rblockDiscrTest
!</input>

!<output>
  ! The 1x1 block matrix, created from rscalarMat.
  type(t_matrixBlock), intent(OUT) :: rmatrix
!</output>
  
!</subroutine>

    ! Fill the rmatrix structure with data.
    rmatrix%NEQ           = rscalarMat%NEQ * rscalarMat%NVAR
    rmatrix%NCOLS         = rscalarMat%NCOLS * rscalarMat%NVAR
    rmatrix%nblocksPerCol = 1
    rmatrix%nblocksPerRow = 1
    rmatrix%imatrixSpec   = LSYSBS_MSPEC_SCALAR
    
    ! Copy the content of the scalar matrix structure into the
    ! first block of the block matrix
    allocate(rmatrix%RmatrixBlock(1,1))
    rmatrix%RmatrixBlock(1,1) = rscalarMat

    ! The matrix is a copy of another one. Note this!
    rmatrix%RmatrixBlock(1,1)%imatrixSpec = &
      ior(rmatrix%RmatrixBlock(1,1)%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
      
    ! Save a pointer to the discretisation structure if that one
    ! is specified.
    if (present(rblockDiscrTrial)) then
      rmatrix%p_rblockDiscrTrial => rblockDiscrTrial
      if (present(rblockDiscrTest)) then
        rmatrix%p_rblockDiscrTest => rblockDiscrTest
        rmatrix%bidenticalTrialAndTest = .false.
      else
        rmatrix%p_rblockDiscrTest => rblockDiscrTrial
        rmatrix%bidenticalTrialAndTest = .true.
      end if
    else
      nullify(rmatrix%p_rblockDiscrTest)
      nullify(rmatrix%p_rblockDiscrTrial)
      rmatrix%bidenticalTrialAndTest = .true.
    end if
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_createVecFromScalar (rscalarVec,rvector,&
                                         rblockDiscretisation,nblocks)
  
!<description>
  ! This routine creates a 1-block vector rvector from a scalar vector
  ! rscalarVec. Both, rscalarVec and rvector will share the same handles,
  ! so changing the content of rvector will change rscalarVec, too.
  ! Therefore, the bisCopy flag of the subvector in the block vector
  ! will be set to TRUE.
!</description>
  
!<input>
  ! The scalar vector which should provide the data
  type(t_vectorScalar), intent(IN) :: rscalarVec

  ! OPTIONAL: A block discretisation structure specifying the trial functions.
  ! A pointer to this will be saved to the vector.
  type(t_blockDiscretisation), intent(IN), optional, target :: rblockDiscretisation
  
  ! OPTIONAL: Number of blocks to reserve.
  ! Normally, the created scalar vector has only one block. If nblocks
  ! is specified, the resulting vector will have more blocks while only
  ! the first block vector is used.
  integer, intent(IN), optional :: nblocks
!</input>

!<output>
  ! The block vector, created from rscalarVec.
  type(t_vectorBlock), intent(OUT) :: rvector
!</output>
  
!</subroutine>

    integer :: nactblocks
    
    nactblocks = 1
    if (present(nblocks)) nactblocks = max(nblocks,nactblocks)
    if (present(rblockDiscretisation)) &
      nactblocks = max(rblockDiscretisation%ncomponents,nactblocks)

    ! Fill the rvector structure with data.
    rvector%NEQ         = rscalarVec%NEQ * rscalarVec%NVAR
    rvector%cdataType   = rscalarVec%cdataType
    rvector%h_Ddata     = rscalarVec%h_Ddata
    rvector%nblocks     = nactblocks
    
    ! Copy the content of the scalar matrix structure into the
    ! first block of the block vector
    allocate(rvector%RvectorBlock(nactblocks))
    rvector%RvectorBlock(1) = rscalarVec
    
    ! Copy the starting address of the scalar vector to our block vector.
    rvector%iidxFirstEntry  = rscalarVec%iidxFirstEntry

    ! The handle belongs to another vector - note this in the structure.
    rvector%RvectorBlock(1)%bisCopy     = .true.
    
    ! The data of the vector actually belongs to another one - to the
    ! scalar one. Note this in the structure.
    rvector%bisCopy = .true.

    ! Save a pointer to the discretisation structure if that one
    ! is specified.
    if (present(rblockDiscretisation)) then
      rvector%p_rblockDiscr => rblockDiscretisation
    else
      nullify(rvector%p_rblockDiscr)
    end if
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_createScalarFromVec (rvector,rscalarVec,bshare)
  
!<description>
  ! This routine creates a scalar vector rscalarVec from a block vector rvector.
  ! rscalarVec will contain the information from all subvectors in rvector
  ! in the order specified by rvector.
  ! If bshare=TRUE, the routine shares the information of rvector
  ! with rscalarVec without allocating new memory, so changing the content of 
  ! rscalarVec will also change rvector. The bisCopy flag of the subvector 
  ! in the block vector will be set to TRUE.
  !
  ! If the source vector contains exactly one subvector, the routine copies
  ! additional information (pointer to the discretisation structure) to 
  ! rscalarVec. Otherwise, rscalarVec will contain only the basic data.
!</description>
  
!<input>
  ! The block vector which should provide the data
  type(t_vectorBlock), intent(IN) :: rvector
  
  ! OPTIONAL: Whether to share data.
  ! = FALSE: Create a new vector, copy all data.
  ! = TRUE: Create a vector that shares information with rvector.
  !         If that's not possible, create a new vector. (Standard)
  logical, optional :: bshare
!</input>

!<output>
  ! Scalar vector created from rvector.
  type(t_vectorScalar), intent(OUT) :: rscalarVec
!</output>
  
!</subroutine>

    ! local variables
    logical :: bshareData
    integer :: istart
    real(DP), dimension(:), pointer :: p_Dsource,p_Ddest
    real(SP), dimension(:), pointer :: p_Fsource,p_Fdest
    
    bshareData = .true.
    if (present(bshare)) bshareData = bshare
    
    if (bshareData) then
    
      ! Share the data
      rscalarVec%cdataType = rvector%cdataType
      rscalarVec%NEQ = rvector%NEQ
      rscalarVec%h_Ddata = rvector%h_Ddata
      rscalarVec%bisCopy = .true.
      
    else
    
      ! Is there only one subvector? If yes, copy that one!
      ! Allocate new memory if necessary.
      if (rvector%nblocks .eq. 1) then
      
        call lsyssc_duplicateVector (rvector%RvectorBlock(1),rscalarVec,&
            LSYSSC_DUP_COPY,LSYSSC_DUP_COPY)
            
      else
    
        ! Create a new vector.
        call lsyssc_createVector (rscalarVec,rvector%NEQ,.false.,rvector%cdataType)
        
        ! Copy the data

        istart = 1
        select case (rvector%cdataType)
        case (ST_DOUBLE)
          call lsyssc_getbase_double (rscalarVec,p_Ddest)
          call lsysbl_getbase_double (rvector,p_Dsource)
          call lalg_copyVectorDble (p_Dsource(:),p_Ddest(:))
          
        case (ST_SINGLE)
          call lsyssc_getbase_single (rscalarVec,p_Fdest)
          call lsysbl_getbase_single (rvector,p_Fsource)
          call lalg_copyVectorSngl (p_Fsource(:),p_Fdest(:))
          
        case DEFAULT
        
          print *,'lsysbl_createScalarFromVec: Invalid data type!'
          call sys_halt()
          
        end select
        
      end if
        
    end if
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_enforceStructure (rtemplateVec,rvector)
  
!<description>
  ! This routine enforces the structure of the vector rtemplate
  ! in the vector rvector. 
  !
  ! WARNING: This routine should be used with care if you know
  !          what you are doing !!!
  ! All structural data of rtemplate (boundary conditions, spatial 
  ! discretisation, size,sorting,...) are copied from rtemplate to 
  ! rvector without checking the compatibility!!!
  !
  ! The only check in this routine is that rvector%NEQ is
  ! at least as large as rtemplate%NEQ; otherwise an error
  ! is thrown. The data type of rvector is also not changed.
!</description>
  
!<input>
  ! A template vector specifying structural data
  type(t_vectorBlock), intent(IN) :: rtemplateVec
!</input>

!<inputoutput>
  ! The vector which structural data should be overwritten
  type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: cdata,h_Ddata,i
  logical :: biscopy
  integer :: istart,n !,length
  type(t_vectorScalar), dimension(:), pointer :: p_Rblocks

!  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
!  REAL(SP), DIMENSION(:), POINTER :: p_Fdata
!
! Here, we check the real size of the array on the heap, not
! simply NEQ. Allows to enlarge a vector if there's enough space.
!    SELECT CASE (rvector%cdataType)
!    CASE (ST_DOUBLE)
!      CALL storage_getbase_double (rvector%h_Ddata,p_Ddata)
!      length = SIZE(p_Ddata)-rvector%iidxFirstEntry+1
!    CASE (ST_SINGLE)
!      CALL storage_getbase_double (rvector%h_Ddata,p_Fdata)
!      length = SIZE(p_Fdata)-rvector%iidxFirstEntry+1
!    CASE DEFAULT
!      PRINT *,'lsysbl_enforceStructure: Unsupported data type'
!      CALL sys_halt()
!    END SELECT
    
    ! Only basic check: there must be enough memory.
    if (rvector%NEQ .lt. rtemplateVec%NEQ) then
      print *,'lsysbl_enforceStructure: Destination vector too small!'
      call sys_halt()
    end if
  
    ! Get data type and handle from rvector
    cdata = rvector%cdataType
    h_Ddata = rvector%h_Ddata
    istart = rvector%iidxFirstEntry
    biscopy = rvector%bisCopy
    p_Rblocks => rvector%RvectorBlock
    
    ! Overwrite rvector
    rvector = rtemplateVec
    
    ! Restore the data array
    rvector%cdataType = cdata
    rvector%h_Ddata = h_Ddata
    rvector%iidxFirstEntry = istart
    rvector%bisCopy = biscopy
    rvector%RvectorBlock => p_Rblocks
    
    ! If RvectorBlock is too small, reallocate
    if (size(rvector%RvectorBlock) .lt. rtemplateVec%nblocks) then
      deallocate(rvector%RvectorBlock)
      allocate(rvector%RvectorBlock(rtemplateVec%nblocks))
    end if
    
    ! Copy the subvector data
    rvector%RvectorBlock = rtemplateVec%RvectorBlock
    
    ! Relocate the starting indices of the subvectors.
    ! Correct the handles of the subvectors since they point to
    ! the data of rtemplateVec because of the previous copy-command.
    n = istart
    do i=1,rvector%nblocks
      rvector%RvectorBlock(i)%h_Ddata = h_Ddata
      rvector%RvectorBlock(i)%iidxFirstEntry = n
      n = n + rvector%RvectorBlock(i)%NEQ
    end do
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_enforceStructureDirect (Isize,rvector)
  
!<description>
  ! This routine enforces the structure of the vector rtemplate
  ! in the vector rvector. 
  !
  ! WARNING: This routine should be used with care if you know
  !          what you are doing !!!
  ! 
  ! The routine enforces SIZE(Isize) subvectors in rvector and
  ! sets the size of the i'thg subvector to Isize(i).
  !
  ! The only check in this routine is that rvector%NEQ is
  ! at least as large as SUM(Isize); otherwise an error
  ! is thrown. The data type of rvector is also not changed.
!</description>
  
!<input>
  ! An array with size definitions. Isize(j) specifies the length of the j'th
  ! subvector.
  integer, dimension(:), intent(IN) :: Isize
!</input>

!<inputoutput>
  ! The vector which structural data should be overwritten
  type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>
  
!</subroutine>

  ! local variables
  integer :: i
  integer :: n !,length

!  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
!  REAL(SP), DIMENSION(:), POINTER :: p_Fdata
!
    ! Only basic check: there must be enough memory.
    if (rvector%NEQ .lt. sum(Isize)) then
      print *,'lsysbl_enforceStructureDirect: Destination vector too small!'
      call sys_halt()
    end if
  
    ! Set the attributes of the vector
    rvector%nblocks = size(Isize)
   
    ! If RvectorBlock is wrong, reallocate
    if (size(rvector%RvectorBlock) .ne. rvector%nblocks) then
      deallocate(rvector%RvectorBlock)
      allocate(rvector%RvectorBlock(rvector%nblocks))
    end if

    rvector%RvectorBlock(1:rvector%nblocks)%h_Ddata = rvector%h_Ddata
    rvector%RvectorBlock(1:rvector%nblocks)%bisCopy = .true.
    rvector%RvectorBlock(1:rvector%nblocks)%cdataType = rvector%cdataType
    
    ! Relocate the starting indices of the subvectors.
    n = rvector%iidxFirstEntry
    do i=1,rvector%nblocks
      rvector%RvectorBlock(i)%iidxFirstEntry = n
      rvector%RvectorBlock(i)%NEQ = Isize(i)
      n = n + Isize(i)
    end do
    
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_enforceStructureDiscr (rdiscretisation,rx)
  
!<description>
  ! This routine enforces the structure of the discretisation rdiscretisaation
  ! in the vector rvector. 
  !
  ! WARNING: This routine should be used with care if you know
  !          what you are doing !!!
  ! 
  ! The routine can be used to save space if multiple vectors with a different
  ! structure share the same memory array. Nevertheless, the caller must
  ! take care of which information is aactually in the arrays.
  !
  ! The only check in this routine is that the free space in the vector is
  ! at least as large as NEQ specified by the discretisation structure; 
  ! otherwise an error is thrown. The data type of rvector is also not 
  ! changed.
!</description>
  
!<input>
  ! Destination structure. Discretisation-related information of rdiscretisation
  ! is assigned to rx.
  type(t_blockDiscretisation),intent(IN), target :: rdiscretisation
!</input>

!<inputoutput>
  ! Destination vector which structure should be changed.
  type(t_vectorBlock),intent(INOUT) :: rx
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: i
    integer :: iidx

    ! Get current size of vector memory
    call storage_getsize(rx%h_Ddata, iidx)

    ! Check that the destination vector can handle this structure
    if (dof_igetNDofGlobBlock(rdiscretisation) .gt. iidx) then
      print *,'lsysbl_enforceStructureDiscr: Destination vector too small!'
      call sys_halt()
    end if

    iidx = 1

    do i=1,rx%nblocks
      ! Modify all pointers of the spatial discretisation in all subvectors.
      rx%RvectorBlock(i)%p_rspatialDiscr => &
        rdiscretisation%RspatialDiscr(i)
        
      ! Recalculate the size of the subvectors
      rx%RvectorBlock(i)%NEQ = &
        dof_igetNDofGlob(rx%RvectorBlock(i)%p_rspatialDiscr)
      rx%RvectorBlock(i)%iidxFirstEntry = iidx
      
      iidx = iidx + rx%RvectorBlock(i)%NEQ
    end do
    
    ! Total size of the vector
    rx%NEQ = iidx-1
    
    ! Set a pointer to the discretisation structure in the vector, finish
    rx%p_rblockDiscr => rdiscretisation
  
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_updateMatStrucInfo (rmatrix)
  
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
  type(t_matrixBlock), intent(INOUT) :: rmatrix
!</inputoutput>

!</subroutine>

  integer :: j,k,NEQ,NCOLS
  
  ! At first, recalculate the number of blocks in the block matrix
  rmatrix%nblocksPerCol = ubound(rmatrix%RmatrixBlock,1)
  rmatrix%nblocksPerRow = ubound(rmatrix%RmatrixBlock,2)
  
  ! Calculate the new NEQ and NCOLS. Go through all 'columns' and 'rows
  ! of the block matrix and sum up their dimensions.
  NCOLS = 0
  do j = 1, rmatrix%nblocksPerRow
    do k = 1, rmatrix%nblocksPerCol
      if (rmatrix%RmatrixBlock(k,j)%NCOLS .ne. 0) then
        NCOLS = NCOLS + rmatrix%RmatrixBlock(k,j)%NCOLS *&
                        rmatrix%RmatrixBlock(k,j)%NVAR
        exit
      end if
    end do
  end do

  ! Loop in the transposed way to calculate NEQ.
  NEQ = 0
  do j = 1, rmatrix%nblocksPerCol
    do k = 1, rmatrix%nblocksPerRow
      if (rmatrix%RmatrixBlock(j,k)%NEQ .ne. 0) then
        NEQ = NEQ + rmatrix%RmatrixBlock(j,k)%NEQ *&
                    rmatrix%RmatrixBlock(j,k)%NVAR
        exit
      end if
    end do
  end do
  
  rmatrix%NEQ = NEQ
  rmatrix%NCOLS = NCOLS

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_clearMatrix (rmatrix,dvalue)
  
!<description>
  ! Clears the entries in all submatrices of a block matrix. 
  ! All entries are overwritten with 0.0 or with dvalue (if specified).
!</description>
  
!<inputoutput>
  ! The block matrix which is to be updated.
  type(t_matrixBlock), intent(INOUT) :: rmatrix

  ! OPTIONAL: Value to write into the matrix.
  ! If not specified, all matrix entries are set to 0.0.
  ! If specified, all matrix entries are set to dvalue.
  real(DP), intent(IN), optional :: dvalue
!</inputoutput>

!</subroutine>

  integer :: i,j
  
  ! Loop through all blocks and clear the matrices
  ! block matrix
  do i=1,rmatrix%nblocksPerCol
    do j=1,rmatrix%nblocksPerRow
      if (lsysbl_isSubmatrixPresent (rmatrix,i,j)) &
        call lsyssc_clearMatrix (rmatrix%RmatrixBlock(i,j),dvalue)
    end do
  end do

  end subroutine

  !****************************************************************************
  
!<subroutine>
  
  subroutine lsysbl_duplicateMatrix (rsourceMatrix,rdestMatrix,&
                                     cdupStructure, cdupContent)
  
!<description>
  ! Duplicates an existing matrix, creates a new matrix rdestMatrix based
  ! on a template matrix rsourceMatrix. To duplicate a block matrix means here
  ! to copy all the existing submatrices in the block matrix in the same way.
  !
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
  !  does not maintain it. Therefore, copying a permutation in one of
  !  the submatrices means copying the corresponding handle. 
  !  The application must keep track of the permutations.
!</description>
  
!<input>
  ! Source matrix.
  type(t_matrixBlock), intent(IN)               :: rsourceMatrix
  
  ! Duplication flag that decides on how to set up the structure
  ! of rdestMatrix. This duplication flag is applied to all submatrices
  ! of rsourceMatrix.
  !
  ! One of the LSYSSC_DUP_xxxx flags:
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
  integer, intent(IN)                            :: cdupStructure
  
  ! Duplication flag that decides on how to set up the content
  ! of rdestMatrix. This duplication flag is applied to all submatrices
  ! of rsourceMatrix.
  !
  ! One of the LSYSSC_DUP_xxxx flags:
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
  integer, intent(IN)                            :: cdupContent
!</input>

!<output>
  ! Destination matrix.
  type(t_matrixBlock), intent(INOUT)            :: rdestMatrix
!</output>  

!</subroutine>

    ! local variables
    integer :: i,j
    type(t_matrixScalar), dimension(:,:), pointer :: p_rblocks
    
    ! Copy the matrix structure of rsourceMatrix to rdestMatrix. This will also
    ! copy the structures of the submatrices, what we have to correct later.
    p_rblocks => rdestMatrix%RmatrixBlock
    rdestMatrix = rsourceMatrix
    rdestMatrix%RmatrixBlock => p_rblocks
    if (.not. associated(rdestMatrix%RmatrixBlock)) then
      ! Check if number of diagonal blocks is nonzero, otherwise exit
      if ((rdestMatrix%nblocksPerCol <= 0) .or. &
          (rdestMatrix%nblocksPerRow <= 0)) return
      allocate(rdestMatrix%RmatrixBlock(rdestMatrix%nblocksPerCol,rdestMatrix%nblocksPerRow))
    end if
    rdestMatrix%RmatrixBlock = rsourceMatrix%RmatrixBlock
    
    ! For every submatrix in the source matrix, call the 'scalar' variant
    ! of duplicateMatrix. Dismiss all old information and replace it by
    ! the new one.
    do j=1,rsourceMatrix%nblocksPerRow
      do i=1,rsourceMatrix%nblocksPerCol
        if (rdestMatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIXUNDEFINED) then
        
          ! Set the specification flags to "matrix is copy", so the releaseMatrix
          ! routine will not release any memory but only dismiss all information.
          rdestMatrix%RmatrixBlock(i,j)%imatrixSpec = &
              ior(rdestMatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
          call lsyssc_releaseMatrix (rdestMatrix%RmatrixBlock(i,j))
              
          call lsyssc_duplicateMatrix ( &
              rsourceMatrix%RmatrixBlock(i,j), rdestMatrix%RmatrixBlock(i,j),&
              cdupStructure, cdupContent)
                                       
        end if
      end do
    end do
    
    ! Update the structural information of the block matrix for completeness.
    call lsysbl_updateMatStrucInfo(rdestMatrix)
    
  end subroutine
    
  !****************************************************************************
  
!<subroutine>
  
  subroutine lsysbl_duplicateVector (rx,ry,&
                                     cdupStructure, cdupContent)
  
!<description>
  ! Duplicates an existing vector: ry := rx.
  ! Creates a new vector rdestVector based
  ! on a template vector rsourceVector. To duplicate a block vector means here
  ! to copy all the existing subvectors in the block vector in the same way.
  !
  ! Duplicating a vector does not necessarily mean that new memory is
  ! allocated and the vector entries are copied to that. The two flags
  ! cdupStructure and cdupContent decide on how to set up rdestVector.
  ! Depending on their setting, it's possible to copy only then handles
  ! of such dynamic information, so that both vectors share the same
  ! information.
  !
  ! We follow the following convention:
  !  Structure = NEQ, sorting permutation(s), discretisation-related 
  !              information, boundary conditions.
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
  ! Source vector.
  type(t_vectorBlock), intent(IN) :: rx

  ! Duplication flag that decides on how to set up the structure
  ! of ry. Not all flags are possible!
  ! One of the LSYSSC_DUP_xxxx flags:
  ! LSYSSC_DUP_IGNORE     : Ignore the structure of rx
  ! LSYSSC_DUP_COPY or
  ! LSYSSC_DUP_COPYOVERWRITE or
  ! LSYSSC_DUP_TEMPLATE   : Structural data is copied from rx
  !   to ry (NEQ, sorting strategy, pointer to discretisation structure).
  integer, intent(IN)                            :: cdupStructure
  
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
  integer, intent(IN)                            :: cdupContent
!</input>

!<inputoutput>
  ! Destination vector.
  type(t_vectorBlock), intent(INOUT) :: ry
!</inputoutput>

!</subroutine>

    integer :: h_Ddata,i
    integer :: cdataType
    logical :: bisCopy
    type(t_vectorScalar), dimension(:), pointer :: p_Rblocks
    integer ::  isize,ioffset

    ! First the structure

    if (cdupStructure .ne. LSYSSC_DUP_IGNORE) then
      ! First, make a backup of some crucial data so that it does not
      ! get destroyed.
      h_Ddata = ry%h_Ddata
      cdataType = ry%cdataType
      bisCopy = ry%bisCopy
      p_Rblocks => ry%RvectorBlock
      
      ! Then transfer all structural information of rx to ry.
      ! This automatically makes both vectors compatible to each other.
      ry = rx
      
      ! Restore crucial data
      ry%h_Ddata = h_Ddata
      ry%bisCopy = bisCopy
      ry%RvectorBlock => p_Rblocks
      
      ! If necessary, allocate new memory for the blocks.
      ! Note that ry%nblocks is now =rx%nblocks!
      if (.not. associated(ry%RvectorBlock)) then
        allocate(ry%RvectorBlock(ry%nblocks))
      else if (size(ry%RvectorBlock) .lt. rx%nblocks) then
        deallocate(ry%RvectorBlock)
        allocate(ry%RvectorBlock(ry%nblocks))
      end if
      
      ! Copy the block structure. Don't destroy crucial data.
      do i=1,size(ry%RvectorBlock)
        h_Ddata = ry%RvectorBlock(i)%h_Ddata
        cdataType = ry%RvectorBlock(i)%cdataType
        bisCopy = ry%RvectorBlock(i)%bisCopy
        ioffset = ry%RvectorBlock(i)%iidxFirstEntry
        
        ! Transfer all structural information of rx to ry.
        ! This automatically makes both vectors compatible to each other.
        ry%RvectorBlock(i) = rx%RvectorBlock(i)
        
        ! Restore crucial data
        ry%RvectorBlock(i)%h_Ddata = h_Ddata
        ry%RvectorBlock(i)%bisCopy = bisCopy
        ry%RvectorBlock(i)%cdataType = cdataType
        ry%RvectorBlock(i)%iidxFirstEntry = ioffset 
      end do

    end if
    
    ! Now the content

    select case (cdupContent)
    case (LSYSSC_DUP_IGNORE) 
      ! Nothing to do
    
    case (LSYSSC_DUP_REMOVE)
      ! Release vector data
      if ((.not. ry%bisCopy) .and. (ry%h_Ddata .ne. ST_NOHANDLE)) then
        call storage_free (ry%h_Ddata)
      end if
      ry%bisCopy = .false.
      ry%iidxFirstEntry = 1
      
      ry%RvectorBlock(1:ry%nblocks)%h_Ddata = ST_NOHANDLE
      ry%RvectorBlock(1:ry%nblocks)%bisCopy = .false.
      ry%RvectorBlock(1:ry%nblocks)%iidxFirstEntry = 1
      
    case (LSYSSC_DUP_DISMISS)
      ! Dismiss data
      ry%h_Ddata = ST_NOHANDLE
      ry%bisCopy = .false.
      ry%iidxFirstEntry = 1
      
      ry%RvectorBlock(1:ry%nblocks)%h_Ddata = ST_NOHANDLE
      ry%RvectorBlock(1:ry%nblocks)%bisCopy = .false.
      ry%RvectorBlock(1:ry%nblocks)%iidxFirstEntry = 1
    
    case (LSYSSC_DUP_SHARE)
      ! Share information. Release memory if necessary.
      if ((.not. ry%bisCopy) .and. (ry%h_Ddata .ne. ST_NOHANDLE)) then
        call storage_free (ry%h_Ddata)
      end if
      ry%h_Ddata = rx%h_Ddata
      ry%cdataType = rx%cdataType
      ry%bisCopy = .true.
      
      ! Restore the handle in all subvectors
      ry%RvectorBlock(1:ry%nblocks)%h_Ddata = ry%h_Ddata
      ry%RvectorBlock(1:ry%nblocks)%cdataType = ry%cdataType
      ry%RvectorBlock(1:ry%nblocks)%bisCopy = .true.
      
      ! Set the starting positions of all subvectors as well as the starting position
      ! of ry to that of rx -- we now have the same handle as rx.
      call updateIndex (ry,rx%iidxFirstEntry)
      
    case (LSYSSC_DUP_COPY)
      ! Copy information. If necessary, allocate new data -- by setting
      ! h_Ddata to ST_NOHANDLE before calling storage_copy.
      if (ry%bisCopy) then
        ry%h_Ddata = ST_NOHANDLE
        ry%bisCopy = .false.
      end if
      
      if (ry%h_Ddata .ne. ST_NOHANDLE) then    
        call storage_getsize (ry%h_Ddata, isize)
        if (isize .lt. rx%NEQ) then
          ! Reallocate by first releasing the data
          call storage_free (ry%h_Ddata)
        end if
      end if
      
      ! Allocate memory if necessary.
      if (ry%h_Ddata .eq. ST_NOHANDLE) &
        call newMemory (rx,ry)
      
      ! Copy the data directly.
      call lsysbl_copyVectorDirect (rx,ry)
    
    case (LSYSSC_DUP_COPYOVERWRITE)
    
      ! Some basic checks
      if (ry%h_Ddata .ne. ST_NOHANDLE) then    
        call storage_getsize (ry%h_Ddata, isize)
        if (isize .lt. rx%NEQ) then
          ! Reallocate by first releasing the data
          call storage_free (ry%h_Ddata)
          ry%bisCopy = .false.
        end if
      end if
      
      ! Copy information, regardless of whether ry is the owner or not.
      ! Allocate memory if necessary.
      if (ry%h_Ddata .eq. ST_NOHANDLE) &
        call newMemory (rx,ry)
      
      ! Copy the data directly.
      call lsysbl_copyVectorDirect (rx,ry)
    
    case (LSYSSC_DUP_ASIS)
    
      ! Copy by ownership. This is either LSYSSC_COPY or LSYSSC_SHARE,
      ! depending on whether rx is the owner of the data or not.
      if (rx%bisCopy) then
        
        ! rx shares it's data and thus ry will also.
        ry%h_Ddata = rx%h_Ddata
        ry%cdataType = rx%cdataType
        ry%bisCopy = .true.

        ! Set up the sub-blocks
        ry%RvectorBlock(1:ry%nblocks)%h_Ddata = ry%h_Ddata
        ry%RvectorBlock(1:ry%nblocks)%cdataType = ry%cdataType
        ry%RvectorBlock(1:ry%nblocks)%bisCopy = .true.

        ! Set the starting positions of all subvectors as well as the starting position
        ! of ry to that of rx -- we now have the same handle as rx.
        call updateIndex (ry,rx%iidxFirstEntry)
        
      else
      
        ! The data belongs to rx and thus it must also belong to ry --
        ! So copy information. If necessary, allocate new data -- by setting
        ! h_Ddata to ST_NOHANDLE before calling storage_copy.
        if (ry%bisCopy) then
          ry%h_Ddata = ST_NOHANDLE
          ry%bisCopy = .false.
        end if
        
        if (ry%h_Ddata .ne. ST_NOHANDLE) then    
          call storage_getsize (ry%h_Ddata, isize)
          if (isize .lt. rx%NEQ) then
            ! Reallocate by first releasing the data
            call storage_free (ry%h_Ddata)
          end if
        end if
        
        ! Allocate memory if necessary.
        if (ry%h_Ddata .eq. ST_NOHANDLE) &
          call newMemory (rx,ry)
        
        ! Copy the data directly.
        call lsysbl_copyVectorDirect (rx,ry)

      end if
      
    case (LSYSSC_DUP_EMPTY)
    
      ! Allocate new memory if ry is empty. Don't initialise.
      ! If ry contains data, we don't have to do anything.
      if (ry%h_Ddata .eq. ST_NOHANDLE) then
        call newMemory (rx,ry)
      end if

    case DEFAULT
    
      print *,'lsysbl_duplicateVector: cdupContent unknown!'
      call sys_halt()
    
    end select
  
  contains
  
    subroutine newMemory (rx,ry)
    
    ! Creates memory in ry for a vector in the same size and data type as rx.
    
    ! Source vector that specifies the structure
    type(t_vectorBlock), intent(IN) :: rx
    
    ! Destination vector. Memory is allocated here. Block structure is
    ! not changed.
    type(t_vectorBlock), intent(INOUT) :: ry

      integer :: ioffset,i

      ! Allocate memory   
      ry%cdataType = rx%cdataType 
      ry%bisCopy = .false.
      call storage_new ('lsyssc_duplicateVector','vec-copy',rx%NEQ,&
                        ry%cdataType, ry%h_Ddata, ST_NEWBLOCK_NOINIT)
      
      ! Transfer the information of the new data source to the subvectors.
      ry%RvectorBlock(1:ry%nblocks)%h_Ddata = ry%h_Ddata
      ry%RvectorBlock(1:ry%nblocks)%cdataType = ry%cdataType
      ry%RvectorBlock(1:ry%nblocks)%bisCopy = .true.
            
      ! Create the starting positions of the subvectors.
      call updateIndex (ry,1)
    
    end subroutine

    ! ---------------------------------------------------------------

    subroutine updateIndex (ry,iidxFirstEntry)
    
    ! Updates the index position in ry based on an initial position iidxFirstEntry
    ! and the existing structure in ry.
    
    ! Vector whose subvectors should be corrected
    type(t_vectorBlock), intent(INOUT) :: ry
    
    ! New first position of the first subvector in the global data array.
    integer, intent(IN) :: iidxFirstEntry

      integer :: ioffset,i

      ! Create the starting positions of the subvectors.
      ry%iidxFirstEntry = iidxFirstEntry
      ry%RvectorBlock(1)%iidxFirstEntry = iidxFirstEntry
      do i=2,ry%nblocks
        ry%RvectorBlock(i)%iidxFirstEntry = ry%RvectorBlock(i-1)%iidxFirstEntry + &
          ry%RvectorBlock(i-1)%NEQ
      end do
    
    end subroutine
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_sortVectorInSitu (rvector,rtemp,bsort)
  
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
  logical, intent(IN) :: bsort
!</input>
  
!<inputoutput>
  ! The vector which is to be resorted
  type(t_vectorBlock), intent(INOUT) :: rvector

  ! A scalar temporary vector. Must be of the same data type as rvector.
  ! Must be at least as large as the longest subvector in rvector.
  type(t_vectorScalar), intent(INOUT) :: rtemp
!</inputoutput>

!</subroutine>

  integer :: iblock
  
  ! Loop over the blocks
  if (bsort) then
    do iblock = 1,rvector%nblocks
      call lsyssc_sortVectorInSitu (&
          rvector%RvectorBlock(iblock), rtemp,&
          abs(rvector%RvectorBlock(iblock)%isortStrategy))
    end do
  else
    do iblock = 1,rvector%nblocks
      call lsyssc_sortVectorInSitu (&
          rvector%RvectorBlock(iblock), rtemp,&
          -abs(rvector%RvectorBlock(iblock)%isortStrategy))
    end do
  end if

  end subroutine

  !****************************************************************************
  
!<subroutine>
  
  subroutine lsysbl_synchroniseSortVecVec (rvectorSrc,rvectorDst,rtemp)
  
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
  type(t_vectorBlock), intent(IN)               :: rvectorSrc
!</input>

!<inputoutput>
  ! Destination vector; is resorted according to the sort strategy in rvectorSrc
  ! or is unsorted, if rvectorSrc is unsorted.
  ! Must have the same size as rvectorSrc.
  type(t_vectorBlock), intent(INOUT)            :: rvectorDst

  ! A scalar temporary vector. Must be of the same data type as rvector.
  ! Must be at least as large as the longest subvector in rvector.
  type(t_vectorScalar), intent(INOUT) :: rtemp
!</inputoutput>

!</subroutine>

    integer :: iblock
    
    ! Loop over the blocks
    do iblock = 1,rvectorDst%nblocks
      ! Synchronise every subvector
      call lsyssc_synchroniseSortVecVec (rvectorSrc%RvectorBlock(iblock),&
          rvectorDst%RvectorBlock(iblock),rtemp)
    end do

  end subroutine

  !****************************************************************************
  
!<subroutine>
  
  subroutine lsysbl_synchroniseSortMatVec (rmatrixSrc,rvectorDst,rtemp)
  
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
  type(t_matrixBlock), intent(IN)               :: rmatrixSrc
!</input>

!<inputoutput>
  ! Destination vector; is resorted according to the sort strategy in rmatrixSrc
  ! or is unsorted, if rmatrixSrc is unsorted.
  ! Must have the same size (NEQ) as rmatrixSrc.
  type(t_vectorBlock), intent(INOUT)            :: rvectorDst

  ! A temporary vector. Must be of the same data type as rvector.
  ! Must be at least as large as rvectorDst.
  type(t_vectorScalar), intent(INOUT)            :: rtemp
!</inputoutput>

!</subroutine>

    integer :: i,j
    
    ! Loop over the columns of the block matrix
    do j = 1, rmatrixSrc%nblocksPerRow
      ! Loop through all rows in that column. Find the first matrix that can
      ! provide us with a sorting strategy we can use.
      ! We can assume that all submatrices in that matrix column have 
      ! the same sorting strategy, otherwise something like matrix-vector
      ! multiplication will quickly lead to a program failure...
      do i = 1,rmatrixSrc%nblocksPerCol
        if (rmatrixSrc%RmatrixBlock(i,j)%NEQ .ne. 0) then
          call lsyssc_synchroniseSortMatVec (rmatrixSrc%RmatrixBlock(i,j),&
              rvectorDst%RvectorBlock(i),rtemp)
          ! Next column / subvector
          exit
        end if
      end do
    end do

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_setSortStrategy (rvector,IsortStrategy,Hpermutations)
  
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
  type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>

!<input>
  ! An array of sorting strategy identifiers (SSTRAT_xxxx), each for one
  ! subvector of the global vector. 
  ! The negative value of this identifier is saved to the corresponding
  ! subvector.
  ! DIMENSION(rvector\%nblocks)
  integer, dimension(:), intent(IN) :: IsortStrategy
  
  ! An array of handles. Each handle corresponds to a subvector and
  ! defines a permutation how to resort the subvector.
  ! Each permutation associated to a handle must be of the form
  !    array [1..2*NEQ] of integer  (NEQ=NEQ(subvector))
  ! with entries (1..NEQ)       = permutation
  ! and  entries (NEQ+1..2*NEQ) = inverse permutation.
  ! DIMENSION(rvector\%nblocks)
  integer, dimension(:), intent(IN) :: Hpermutations

!</input>
  
!</subroutine>

  ! Install the sorting strategy in every block. 
  rvector%RvectorBlock(1:rvector%nblocks)%isortStrategy = &
    -abs(IsortStrategy(1:rvector%nblocks))
  rvector%RvectorBlock(1:rvector%nblocks)%h_IsortPermutation = &
    Hpermutations(1:rvector%nblocks)

  end subroutine

  !****************************************************************************

!<function>
  
  logical function lsysbl_isMatrixSorted (rmatrix)
  
!<description>
  ! Loops through all submatrices in a block matrix and checks if any of
  ! them is sorted. If yes, the block matrix is given the state 'sorted'.
  ! If not, the matrix is stated 'unsorted'
!</description>
  
!<input>
  ! Block matrix to check
  type(t_matrixBlock), intent(IN)                  :: rmatrix
!</input>

!<result>
  ! Whether the matrix is sorted or not.
!</result>

!</function>

  integer :: i,j,k
  
  ! Loop through all matrices and get the largest sort-identifier
  k = 0
  do j=1,rmatrix%nblocksPerRow
    do i=1,rmatrix%nblocksPerCol
      if (rmatrix%RmatrixBlock(i,j)%NEQ .ne. 0) then
        k = max(k,rmatrix%RmatrixBlock(i,j)%isortStrategy)
      end if
    end do
  end do
  
  ! If k is greater than 0, at least one submatrix is sorted
  lsysbl_isMatrixSorted = k .gt. 0

  end function

  !****************************************************************************

!<function>
  
  logical function lsysbl_isVectorSorted (rvector)
  
!<description>
  ! Loops through all subvectors in a block vector and checks if any of
  ! them is sorted. If yes, the block vector is given the state 'sorted'.
  ! If not, the vector is stated 'unsorted'
!</description>
  
!<input>
  ! Block vector to check
  type(t_vectorBlock), intent(IN)                  :: rvector
!</input>

!<result>
  ! Whether the vector is sorted or not.
!</result>

!</function>

  integer :: i,k
  
  ! Loop through all matrices and get the largest sort-identifier
  k = 0
  do i=1,rvector%nblocks
    if (rvector%RvectorBlock(i)%NEQ .ne. 0) then
      k = max(k,rvector%RvectorBlock(i)%isortStrategy)
    end if
  end do
  
  ! If k is greater than 0, at least one submatrix is sorted
  lsysbl_isVectorSorted = k .gt. 0

  end function

  !****************************************************************************

!<subroutine>

  subroutine lsysbl_swapVectors(rvector1,rvector2)

!<description>
    ! This subroutine swaps two vectors
!</description>

!<inputoutput>
    ! first block vector
    type(t_vectorBlock), intent(INOUT) :: rvector1

    ! second block vector
    type(t_vectorBlock), intent(INOUT) :: rvector2
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_vectorBlock) :: rvector
    integer :: iblock

    ! Vcetor1 -> Vector
    rvector = rvector1
    rvector%p_rblockDiscr => rvector1%p_rblockDiscr
    rvector%p_rdiscreteBC          => rvector1%p_rdiscreteBC
    rvector%p_rdiscreteBCfict      => rvector1%p_rdiscreteBCfict
    do iblock=1,rvector%nblocks
      rvector1%RvectorBlock(iblock)=rvector%RvectorBlock(iblock)
      rvector1%RvectorBlock(iblock)%p_rspatialDiscr =>&
          rvector%RvectorBlock(iblock)%p_rspatialDiscr
    end do

    ! Vector2 -> Vector1
    rvector1 = rvector2
    rvector1%p_rblockDiscr => rvector2%p_rblockDiscr
    rvector1%p_rdiscreteBC          => rvector2%p_rdiscreteBC
    rvector1%p_rdiscreteBCfict      => rvector2%p_rdiscreteBCfict
    do iblock=1,rvector1%nblocks
      rvector1%RvectorBlock(iblock)=rvector2%RvectorBlock(iblock)
      rvector1%RvectorBlock(iblock)%p_rspatialDiscr =>&
          rvector2%RvectorBlock(iblock)%p_rspatialDiscr
    end do

    ! Vector -> Vector2
    rvector2 = rvector
    rvector2%p_rblockDiscr => rvector%p_rblockDiscr
    rvector2%p_rdiscreteBC          => rvector%p_rdiscreteBC
    rvector2%p_rdiscreteBCfict      => rvector%p_rdiscreteBCfict
    do iblock=1,rvector2%nblocks
      rvector2%RvectorBlock(iblock)=rvector%RvectorBlock(iblock)
      rvector2%RvectorBlock(iblock)%p_rspatialDiscr =>&
          rvector%RvectorBlock(iblock)%p_rspatialDiscr
    end do
  end subroutine lsysbl_swapVectors

  !****************************************************************************

!<subroutine>

  subroutine lsysbl_deriveSubvector(rvectorSrc,rvectorDest, &
      ifirstSubvector,ilastSubvector,bshare)

!<description>
    ! This routine derives a block vector as a subset of another block vector.
    !
    ! ifirstSubvector is the index of a scalar subvector in rvectorSrc
    ! which should be used as the first vector in rvectorDest.
    ! ilastSubvector is the index of a scalar subvector in rvectorSrc
    ! which should be used as the last vector in rvectorDest.
    !
    ! rvectorDest will therefore contain the subvectors 
    ! ifirstSubvector..ilastSubvector of rvectorSrc.
    !
    ! If bshare=TRUE, the vector rvectorDest will be created by copying
    ! handles instead of memory from rvectorSrc. Every change in rvectorDest
    ! will therefore also affect rvectorSrc.
    !
    ! If bshare=FALSE (the standard setting), a new vector will be created
    ! and the content of the specified subvectors will be copied to that.
    !
    ! The newly created block vector will not have any block discretisation 
    ! structure attached!
    ! The caller may therefore want to use lsysbl_assignDiscrDirectVec
    ! to specify the correct discretisation structure after
    ! creating the vector.
!</description>

!<input>
  ! Source block vector
  type(t_vectorBlock), intent(IN) :: rvectorSrc

  ! OPTIONAL: Number of the subvector of rvectorSrc that should be used as first 
  ! subvector in rvectorDest. Default value is =1.
  integer, intent(IN), optional :: ifirstSubvector

  ! OPTIONAL: Number of the subvector of rvectorSrc that should be used as 
  ! last subvector in rvectorDest. Default value is the number of blocks
  ! in rvectorSrc.
  integer, intent(IN), optional :: ilastSubvector

  ! OPTIONAL: Whether to share the content between rvectorSrc and rvectorDest.
  ! = TRUE: Create a 'virtual' copy of rvectorSrc in rvectorDest that uses
  !         the same memory.
  ! =FALSE: Copy the sub-content of rvectorSrc to rvectorDest.
  logical, intent(IN), optional :: bshare
!</input>

!<inputoutput>
    ! Destination block vector. Any previous data is discarded, so the
    ! application must take care of that no allocated handles are inside here!
    type(t_vectorBlock), intent(OUT) :: rvectorDest
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ifirst, ilast, ncount, i, idupflag,n
    logical :: bshareContent
    integer, dimension(:), allocatable :: Isize
    
    ! Evaluate the optional parameters
    ifirst = 1
    ilast = rvectorSrc%nblocks
    bshareContent = .false.
    
    if (present(ifirstSubvector)) then
      ifirst = min(max(ifirst,ifirstSubvector),ilast)
    end if
    
    if (present(ilastSubvector)) then
      ilast = max(min(ilast,ilastSubvector),ifirst)
    end if
    
    if (present(bshare)) then
      bshareContent = bshare
    end if
    
    allocate(Isize(max(rvectorSrc%nblocks,1)))

    ! Let's start. At first, create a new vector based on the old, which contains
    ! only those subvectors specified in ifirst..ilast.
    ncount = ilast-ifirst+1
    Isize(1:ncount) = rvectorSrc%RvectorBlock(ifirst:ilast)%NEQ

    ! Data type
    rvectorDest%cdataType = rvectorDest%cdataType
    
    ! Initialise the sub-blocks. Save a pointer to the starting address of
    ! each sub-block.
    ! Denote in the subvector that the handle belongs to us - not to
    ! the subvector.
    allocate(rvectorDest%RvectorBlock(ncount))
    
    n=1
    do i = 1,ncount
      if (Isize(i) .gt. 0) then
        rvectorDest%RvectorBlock(i)%NEQ = Isize(i)
        rvectorDest%RvectorBlock(i)%iidxFirstEntry = n
        rvectorDest%RvectorBlock(i)%h_Ddata = rvectorSrc%h_Ddata
        rvectorDest%RvectorBlock(i)%cdataType = rvectorSrc%cdataType
        rvectorDest%RvectorBlock(i)%bisCopy = .true.
        n = n+Isize(i)
      else
        rvectorDest%RvectorBlock(i)%NEQ = 0
        rvectorDest%RvectorBlock(i)%iidxFirstEntry = 0
        rvectorDest%RvectorBlock(i)%h_Ddata = ST_NOHANDLE
      end if
    end do
    
    deallocate(Isize)
    
    rvectorDest%NEQ = n-1
    rvectorDest%nblocks = ncount
    
    ! Next question: should we share the vector content or not?
    if (.not. bshareContent) then
      ! Allocate a new large vector holding all data.
      call storage_new ('lsysbl_createVecBlockDirect', 'Vector', &
                          rvectorDest%NEQ, rvectorDest%cdataType, &
                          rvectorDest%h_Ddata, ST_NEWBLOCK_NOINIT)
      rvectorDest%RvectorBlock(1:ncount)%h_Ddata = rvectorDest%h_Ddata
    else
      ! The new vector should be a subvector of the old one. 
      ! That means, we have to set the index of the first entry 
      ! in rvectorDest to the first entry of the specified first 
      ! subvector in rvectorSrc.
      rvectorDest%iidxFirstEntry = rvectorSrc%RvectorBlock(ifirst)%iidxFirstEntry
      
      ! Copy the data handle
      rvectorDest%h_Ddata = rvectorSrc%h_Ddata
      
      ! Furthermore, we have to shift the starting addresses in all
      ! subvectors such that they point to the correct place in the
      ! global vector.
      rvectorDest%RvectorBlock(1:ncount)%iidxFirstEntry = &
          rvectorDest%RvectorBlock(1:ncount)%iidxFirstEntry + &
          rvectorDest%iidxFirstEntry - 1
    end if

    rvectorDest%bisCopy = bshareContent
    
    ! Basically copy the data, overwrite the allocated memory without
    ! reallocating.
    ! If bshareContent=true, just share the data.
    idupflag = LSYSSC_DUP_COPYOVERWRITE
    if (bshareContent) idupflag = LSYSSC_DUP_SHARE

    ! Everything is prepared, just copy -- either only the structure or
    ! the structure as well as the content.
    do i=1,ncount
      call lsyssc_duplicateVector (rvectorSrc%RvectorBlock(i+ifirst-1),&
          rvectorDest%RvectorBlock(i),&
          LSYSSC_DUP_COPY,idupflag)
    end do
    
  end subroutine 

  !****************************************************************************
  
!<subroutine>
  
  subroutine lsysbl_deriveSubmatrix (rsourceMatrix,rdestMatrix,&
                                     cdupStructure, cdupContent,&
                                     ifirstBlock,ilastBlock,&
                                     ifirstBlockCol,ilastBlockCol,&
                                     bignoreScaleFactors)
  
!<description>
  ! This routine derives a block matrix as a subset of another block matrix.
  !
  ! ifirstBlock is the number of the first diagonal block in rsourceMatrix
  ! which should be put to position (1,1) in rdestMatrix.
  ! ilastBlock is the index the last diagonal block in rdestMatrix
  ! which should be put to rdestMatrix.
  !
  ! The newly created block matrix will not have any block discretisation 
  ! structure attached!
  ! The caller may therefore want to use lsysbl_assignDiscrDirectMat
  ! to specify the correct discretisation structure after
  ! creating the matrix.
  !
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
  !  does not maintain it. Therefore, copying a permutation in one of
  !  the submatrices means copying the corresponding handle. 
  !  The application must keep track of the permutations.
  !
  ! Remark 2: When ifirstBlockCol,ilastBlockCol is not specified, the routine
  !  creates a matrix oriented at the diagonal:
  !     rdestmatrix = rsourceMatrix (ifirstBlock:ilastBlock, ifirstBlock:ilastBlock)
  !  If ifirstBlockCol,ilastBlockCol is specified, the routine creates
  !  a submatrix based on X/Y block coordinates. ifirstBlock,ilastBlock in this
  !  case specify the Y-coordinates in the block matrix while 
  !  ifirstBlockCol,ilastBlockCol specify the X-coordinates:
  !     rdestmatrix = rsourceMatrix (ifirstBlock:ilastBlock, ifirstBlockCol:ilastBlockCol)
!</description>
  
!<input>
  ! Source matrix.
  type(t_matrixBlock), intent(IN)               :: rsourceMatrix
  
  ! OPTIONAL: X-coordinate of the block in rsourceMatrix that should be put to 
  ! position (1,1) into rdestMatrix. Default value is =1.
  ! If ifirstBlockY is not specified, this also specifies the Y-coordinate,
  ! thus the diagonal block rsourceMatrix (ifirstBlock,ifirstBlockY) is put
  ! to (1,1) of rdestMatrix.
  integer, intent(IN), optional :: ifirstBlock

  ! OPTIONAL: X-coordinate of the last block in rsourceMatrix that should be put to 
  ! rdestMatrix. Default value is the number of blocks in rsourceMatrix.
  ! If ilastBlockY is not specified, this also specifies the Y-coordinate,
  ! thus the diagonal block rsourceMatrix (ilastBlock,ilastBlock) is put
  ! to (ndiagBlocks,ndiagblocks) of rdestMatrix.
  integer, intent(IN), optional :: ilastBlock

  ! OPTIONAL: Y-coordinate of the block in rsourceMatrix that should be put to 
  ! position (1,1) into rdestMatrix. Default value is ifirstBlock.
  integer, intent(IN), optional :: ifirstBlockCol

  ! OPTIONAL: Number of the last block in rsourceMatrix that should be put to 
  ! rdestMatrix. Default value is ilastBlock.
  integer, intent(IN), optional :: ilastBlockCol

  ! OPTIONAL: Ignore the scaling factors in the source matrix.
  ! TRUE: Submatrices of the source matrix will be copied regardless of
  !   whether dscaleFactor=0.0 or not. Standard.
  ! FALSE: Submatrices of the source matrix will only be copied if
  !   dscaleFacor <> 0.0.
  logical, intent(in), optional                  :: bignoreScaleFactors

  ! Duplication flag that decides on how to set up the structure
  ! of rdestMatrix. This duplication flag is applied to all submatrices
  ! of rsourceMatrix.
  !
  ! One of the LSYSSC_DUP_xxxx flags:
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
  integer, intent(IN)                            :: cdupStructure
  
  ! Duplication flag that decides on how to set up the content
  ! of rdestMatrix. This duplication flag is applied to all submatrices
  ! of rsourceMatrix.
  !
  ! One of the LSYSSC_DUP_xxxx flags:
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
  integer, intent(IN)                            :: cdupContent
!</input>

!<inputoutput>
  ! Destination matrix. 
  ! If the matrix exists and has exactly as many blocks as specified by
  ! ifirstBlock/ilastBlock, the destination matrix is kept and overwritten
  ! as specified by cdupStructure/cdubContent.
  ! If the matrix exists but the block count does not match, the matrix
  ! is released and recreated in the correct size.
  type(t_matrixBlock), intent(INOUT)            :: rdestMatrix
!</inputoutput>  

!</subroutine>

    ! local variables
    integer :: i,j
    integer :: ifirstX, ilastX, ifirstY, ilastY
    logical :: bignoreScale
    
    bignoreScale = .true.
    if (present(bignoreScaleFactors)) bignoreScale = bignoreScaleFactors
    
    ! Evaluate the optional parameters
    ifirstX = 1
    ilastX = rsourceMatrix%nblocksPerRow

    ifirstY = 1
    ilastY = rsourceMatrix%nblocksPerCol
    
    if (present(ifirstBlock)) then
      ifirstY = min(max(ifirstY,ifirstBlock),ilastX)
    end if
    
    if (present(ilastBlock)) then
      ilastY = max(min(ilastY,ilastBlock),ifirstX)
    end if
    
    if (present(ifirstBlockCol)) then
      ifirstX = min(max(ifirstX,ifirstBlockCol),ilastX)
    else
      ifirstX = ifirstY
    end if
    
    if (present(ilastBlockCol)) then
      ilastX = max(min(ilastX,ilastBlockCol),ifirstX)
    else
      ilastX = ilastY
    end if
    
    ! If the destination matrix has the correct size, leave it.
    ! if not, release it and create a new one.
    if ((rdestMatrix%NEQ .ne. 0) .and. &
        ((rdestMatrix%nblocksPerRow .ne. ilastX-ifirstX+1) .or. &
         (rdestMatrix%nblocksPerCol .ne. ilastY-ifirstY+1))) then
      call lsysbl_releaseMatrix (rdestMatrix)
    end if
    
    ! For every submatrix in the source matrix, call the 'scalar' variant
    ! of duplicateMatrix. 
    if (rdestMatrix%NEQ .eq. 0) &
      call lsysbl_createEmptyMatrix(rdestMatrix,(ilastY-ifirstY+1),(ilastX-ifirstX+1))
      
    do j=ifirstX,ilastX
      do i=ifirstY,ilastY
        if (lsysbl_isSubmatrixPresent(rsourceMatrix,i,j,bignoreScale)) then
        
          call lsyssc_duplicateMatrix ( &
              rsourceMatrix%RmatrixBlock(i,j), &
              rdestMatrix%RmatrixBlock(i-ifirstY+1,j-ifirstX+1),&
              cdupStructure, cdupContent)
               
        else if (lsysbl_isSubmatrixPresent(rdestMatrix,&
            i-ifirstY+1,j-ifirstX+1,.true.)) then
        
          ! Release the submatrix in the destination matrix if present.
          ! This is independent of the scale factor (we ignore it!) as
          ! we want to produce a destination matrix that looks like
          ! the source matrix without any skeletons in the cupboard
          ! that may be activated by switching the scaling factors...
          
          call lsyssc_releaseMatrix (&
              rdestMatrix%RmatrixBlock(i-ifirstY+1,j-ifirstX+1))
        
        end if
      end do
    end do
    
    ! Update the structural information of the block matrix for completeness.
    nullify(rdestMatrix%p_rblockDiscrTrial)
    nullify(rdestMatrix%p_rblockDiscrTest)
    call lsysbl_updateMatStrucInfo(rdestMatrix)
    
  end subroutine
    
!****************************************************************************
  
!<subroutine>
  
  subroutine lsysbl_extractSubmatrix (rsourceMatrix,rdestMatrix,&
                                      ifirstBlock,ilastBlock,&
                                      ifirstBlockCol,ilastBlockCol,&
                                      bignoreScaleFactors)
  
!<description>
  ! This routine extracts a submatrix from another source matrix.
  !
  ! ifirstBlock is the number of the first diagonal block in rsourceMatrix
  ! which should be put to position (1,1) in rdestMatrix.
  ! ilastBlock is the index the last diagonal block in rdestMatrix
  ! which should be put to rdestMatrix.
  !
  ! In contrast to lsysbl_deriveSubmatrix, the newly created submatrix receives
  ! all ownership information from the source matrix, leaving the source matrix
  ! as 'shared copy'. This allows to extract a submatrix, do some local modifications
  ! to it and move it back to the original matrix with lsysbl_moveToSubmatrix
  ! when finished:
  !
  ! a) Extract a submatrix:
  !      lsysbl_extractSubmatrix (rsource,rdest,a,b,a,b)
  ! b) Do some local modifications, e.g.
  !      rdest%RmatrixBlock(1,1)%dscaleFactor = 5
  !      ...
  ! c) Move the matrix back to its original position
  !      lsysbl_moveToSubmatrix (rdest,rsource,a,a)
  !
  ! During step b), the original matrix rsource is 'delayed', nothing should be done
  ! with it until reintegration with lsysbl_moveToSubmatrix.
  !
  ! Remark: There is never memory allocated on the heap for the sorting
  !  permutation. A matrix is never the 'owner' of a permutation, i.e.
  !  does not maintain it. Therefore, copying a permutation in one of
  !  the submatrices means copying the corresponding handle. 
  !  The application must keep track of the permutations.
  !
  ! Remark 2: When ifirstBlockCol,ilastBlockCol is not specified, the routine
  !  creates a matrix oriented at the diagonal:
  !     rdestmatrix = rsourceMatrix (ifirstBlock:ilastBlock, ifirstBlock:ilastBlock)
  !  If ifirstBlockCol,ilastBlockCol is specified, the routine creates
  !  a submatrix based on X/Y block coordinates. ifirstBlock,ilastBlock in this
  !  case specify the Y-coordinates in the block matrix while 
  !  ifirstBlockCol,ilastBlockCol specify the X-coordinates:
  !     rdestmatrix = rsourceMatrix (ifirstBlock:ilastBlock, ifirstBlockCol:ilastBlockCol)
  !
  ! Remark 3: If the destination matrix is completely empty, a block matrix
  !  with no discretisation structure is created. The 'local' discretisation
  !  structures in the scalar matrices will remain.
  !  If the destination matrix is based on a discretisation structure, the
  !  discretisation structures of the scalar matrices in the source matrix 
  !  is changed according to this 'guiding' block discretisation structure.
!</description>
  
!<input>
  ! OPTIONAL: X-coordinate of the block in rsourceMatrix that should be put to 
  ! position (1,1) into rdestMatrix. Default value is =1.
  ! If ifirstBlockY is not specified, this also specifies the Y-coordinate,
  ! thus the diagonal block rsourceMatrix (ifirstBlock,ifirstBlockY) is put
  ! to (1,1) of rdestMatrix.
  integer, intent(IN), optional :: ifirstBlock

  ! OPTIONAL: X-coordinate of the last block in rsourceMatrix that should be put to 
  ! rdestMatrix. Default value is the number of blocks in rsourceMatrix.
  ! If ilastBlockY is not specified, this also specifies the Y-coordinate,
  ! thus the diagonal block rsourceMatrix (ilastBlock,ilastBlock) is put
  ! to (ndiagBlocks,ndiagblocks) of rdestMatrix.
  integer, intent(IN), optional :: ilastBlock

  ! OPTIONAL: Y-coordinate of the block in rsourceMatrix that should be put to 
  ! position (1,1) into rdestMatrix. Default value is ifirstBlock.
  integer, intent(IN), optional :: ifirstBlockCol

  ! OPTIONAL: Number of the last block in rsourceMatrix that should be put to 
  ! rdestMatrix. Default value is ilastBlock.
  integer, intent(IN), optional :: ilastBlockCol

  ! OPTIONAL: Ignore the scaling factors in the source matrix.
  ! TRUE: Submatrices of the source matrix will be copied regardless of
  !   whether dscaleFactor=0.0 or not. Standard.
  ! FALSE: Submatrices of the source matrix will only be copied if
  !   dscaleFacor <> 0.0.
  logical, intent(in), optional                  :: bignoreScaleFactors

!</input>

!<inputoutput>
  ! Source matrix.
  type(t_matrixBlock), intent(inout)            :: rsourceMatrix

  ! Destination matrix. 
  ! If necessary, this matrix is regenerated based on the size of the matrix
  ! part to be extracted from rsourceMatrix
  type(t_matrixBlock), intent(inout)            :: rdestMatrix
!</inputoutput>  

!</subroutine>

    ! local variables
    integer :: i,j
    integer :: ifirstX, ilastX, ifirstY, ilastY
    logical :: bignoreScale
    
    bignoreScale = .true.
    if (present(bignoreScaleFactors)) bignoreScale = bignoreScaleFactors
    
    ! Evaluate the optional parameters
    ifirstX = 1
    ilastX = rsourceMatrix%nblocksPerRow

    ifirstY = 1
    ilastY = rsourceMatrix%nblocksPerCol
    
    if (present(ifirstBlock)) then
      ifirstY = min(max(ifirstY,ifirstBlock),ilastX)
    end if
    
    if (present(ilastBlock)) then
      ilastY = max(min(ilastY,ilastBlock),ifirstX)
    end if
    
    if (present(ifirstBlockCol)) then
      ifirstX = min(max(ifirstX,ifirstBlockCol),ilastX)
    else
      ifirstX = ifirstY
    end if
    
    if (present(ilastBlockCol)) then
      ilastX = max(min(ilastX,ilastBlockCol),ifirstX)
    else
      ilastX = ilastY
    end if
    
    ! If the destination matrix has the correct size, leave it.
    ! if not, release it and create a new one.
    if ((rdestMatrix%NEQ .ne. 0) .and. &
        ((rdestMatrix%nblocksPerRow .ne. ilastX-ifirstX+1) .or. &
         (rdestMatrix%nblocksPerCol .ne. ilastY-ifirstY+1))) then
      call lsysbl_releaseMatrix (rdestMatrix)
    end if
    
    ! For every submatrix in the source matrix, call the 'scalar' variant
    ! of duplicateMatrix. 
    if (rdestMatrix%NEQ .eq. 0) &
      call lsysbl_createEmptyMatrix(rdestMatrix,(ilastY-ifirstY+1),(ilastX-ifirstX+1))
    
    ! Move the source matrix to the destination matrix.
    !    
    ! loop over all columns and rows in the source matrix
    do j=ifirstX,ilastX
      do i=ifirstY,ilastY
        ! Copy the submatrix from the source matrix to the destination
        ! matrix.
        if (lsysbl_isSubmatrixPresent(rsourceMatrix,i,j,bignoreScale)) then
        
          call lsyssc_moveMatrix(rsourceMatrix%RmatrixBlock(i,j),&
              rdestMatrix%RmatrixBlock(i-ifirstY+1,j-ifirstX+1))
               
          ! Reassign the discretisation structures and change the spatial
          ! discretisation structures of the submatrices according to the
          ! block discretisation if there is one.
          if (associated(rdestMatrix%p_rblockDiscrTrial) .and. &
              associated(rdestMatrix%p_rblockDiscrTest)) then
            call lsyssc_assignDiscrDirectMat (&
                rdestMatrix%RmatrixBlock(i-ifirstY+1,j-ifirstX+1),&
                rdestMatrix%p_rblockDiscrTrial%RspatialDiscr(j-ifirstX+1),&
                rdestMatrix%p_rblockDiscrTest%RspatialDiscr(i-ifirstY+1))
                
            if (rdestMatrix%RmatrixBlock(i-ifirstY+1,j-ifirstX+1)%NCOLS .ne.&
                dof_igetNDofGlob(rdestMatrix%p_rblockDiscrTrial%&
                RspatialDiscr(j-ifirstX+1))) then
              call output_line('Matrix not compatible with discretisation (NCOLS)!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_extractSubmatrix')
              call sys_halt()
            end if

            if (rdestMatrix%RmatrixBlock(i-ifirstY+1,j-ifirstX+1)%NEQ .ne.&
                dof_igetNDofGlob(rdestMatrix%p_rblockDiscrTest%&
                RspatialDiscr(i-ifirstY+1))) then
              call output_line('Matrix not compatible with discretisation (NEQ)!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_extractSubmatrix')
              call sys_halt()
            end if
          end if

        else if (lsysbl_isSubmatrixPresent(rdestMatrix,&
            i-ifirstY+1,j-ifirstX+1,.true.)) then
        
          ! Release the submatrix in the destination matrix if present.
          ! This is independent of the scale factor (we ignore it!) as
          ! we want to produce a destination matrix that looks like
          ! the source matrix without any skeletons in the cupboard
          ! that may be activated by switching the scaling factors...
          
          call lsyssc_releaseMatrix (&
              rdestMatrix%RmatrixBlock(i-ifirstY+1,j-ifirstX+1))
        
        end if
            
      end do
    end do
    
    ! Update the structural information of the block matrix for completeness.
    nullify(rdestMatrix%p_rblockDiscrTrial)
    nullify(rdestMatrix%p_rblockDiscrTest)
    call lsysbl_updateMatStrucInfo(rdestMatrix)
    
  end subroutine

  !****************************************************************************
  
!<function>
  
  elemental logical function lsysbl_isSubmatrixPresent (rmatrix, &
      irow,icolumn,bignoreScaleFactor) result (bispresent)
  
!<description>
  ! This routine checks if the submatrix at position (irow,icolumn)
  ! is present in the matrix rmatrix or not.
!</description>
  
!<input>
  ! The block matrix.
  type(t_matrixBlock), intent(IN)               :: rmatrix
  
  ! Submatrix row to be checked
  integer, intent(IN) :: irow
  
  ! Submatrix column to be checked
  integer, intent(IN) :: icolumn
  
  ! OPTIONAL: Whether to check the scaling factor.
  ! FALSE: A scaling factor of 0.0 disables a submatrix. 
  !        This is the standard setting.
  ! TRUE: The scaling factor is ignored.
  logical, intent(IN), optional :: bignoreScaleFactor
!</input>

!<output>
  ! Whether the submatrix at position (irow,icolumn) exists or not.
!</output>  

!</function>
    logical :: bscale

    if (present(bignoreScaleFactor)) then
      bscale = bignoreScaleFactor
    else
      bscale = .false.
    end if

    bispresent = &
      (rmatrix%RmatrixBlock(irow,icolumn)%cmatrixFormat .ne. LSYSSC_MATRIXUNDEFINED) &
      .and. (bscale .or. &
             (rmatrix%RmatrixBlock(irow,icolumn)%dscaleFactor .ne. 0.0_DP))

  end function

  !****************************************************************************

!<subroutine>

  subroutine lsysbl_resizeVecBlockDirect (rx, Isize, bclear, bcopy, NEQMAX)

!<description>
  ! Resize the vector block structure rx. Isize is an array
  ! of integers containing the length of the individual blocks.
  ! 
  ! Remark: If the vector structure rx is not initialized, then
  ! it cannot be resized and the routine stops with an error.
  ! If it exists already, then the total memory on the heap is
  ! reallocated if "SUM(Isize) > SIZE(rx)", that is, the new
  ! dimension of rx exceeds that of the old vectors. If the old
  ! vector is larger than the new one, then the dimensions of the
  ! subvectors are only adjusted without reallocation.
  ! If the optional parameter NEQMAX is specified, then this value
  ! is taken as upper bound for the total memory.
  !
  ! If the parameter bclear=.TRUE., then the resized vector is cleared.
  ! Otherwise it is left 'as is'.
  !
  ! Remark: In contrast to scalar vectors where the content of the vector
  ! can be copied upon resize this is not possible for block vectors.
  ! Theoretically, this would be possible. However, some subvectors may
  ! have increased whereas other may have decreased, so that an additional
  ! vector is mandatory to reorganize the global array. This is not efficient.
!</description>

!<input>
    
    ! An array with desired length-tags for the different blocks
    integer, dimension(:), intent(IN) :: Isize

    ! Whether to fill the vector with zero initially
    logical, intent(IN)               :: bclear

    ! OPTIONAL: Whether to copy the content of the vector to the resized one
    logical, intent(IN), optional     :: bcopy

    ! OPTIONAL: Maximum length of the vector
    integer, intent(IN), optional     :: NEQMAX

!</input>

!<inputoutput>

    ! Block vector structure
    type(t_vectorBlock), intent(INOUT)             :: rx

!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorBlock)                 :: rxTmp
    real(DP), dimension(:), pointer     :: p_Ddata,p_DdataTmp
    real(SP), dimension(:), pointer     :: p_Fdata,p_FDataTmp
    integer, dimension(:), pointer      :: p_Idata,p_IdataTmp
    integer :: iNEQ,iisize,i,n
    logical              :: bdocopy

    ! Check, that the vector is not a copy of another (possibly larger) vector
    if (rx%bisCopy) then
      call output_line('A copied vector cannot be resized!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_resizeVecBlockDirect')
      call sys_halt()
    end if
    
    ! Check, if vector has been initialized before
    if (rx%NEQ == 0 .or. rx%h_Ddata == ST_NOHANDLE) then
      call output_line(' A vector can only be resized if it has been created correctly!',&
          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_resizeVecBlockDirect')
      call sys_halt()
    end if

    ! Set working dimensions
    iNEQ = max(0,sum(Isize))
    if (present(NEQMAX)) iNEQ = max(iNEQ,NEQMAX)

    ! Set copy/clear attributes
    bdocopy = (.not.bclear)
    if (present(bcopy)) bdocopy = (bdocopy .and. bcopy)

    ! If the vector should be cleared, then the sorting strategy (if any)
    ! can be ignored and reset. Otherwise, the vector needs to be unsorted
    ! prior to copying some part of it. Afterwards, no sorting strategy is
    ! available in any case.
    if (bdocopy) then
      do i = 1, rx%nblocks
        if (rx%RvectorBlock(i)%isortStrategy > 0) then
          call lsyssc_vectorActivateSorting(rx%RvectorBlock(i),.false.)
        end if
      end do
      ! Make a copy of the unsorted content
      call lsysbl_duplicateVector(rx, rxTmp,&
          LSYSSC_DUP_COPY, LSYSSC_DUP_COPY)
    end if

    ! Get current size of vector memory
    call storage_getsize(rx%h_Ddata, iisize)

    ! Update the global NEQ.
    rx%NEQ = max(0,sum(Isize))

    ! Do we really have to reallocate the vector physically?
    if (rx%NEQ > iisize) then

      ! Yes, so adopt the new size. Note that some extra memory is
      ! allocated if the optional argument NEQMAX is specified
      iisize = iNEQ
      
      ! Reallocate the memory for handle h_Ddata
      call storage_realloc('lsysbl_resizeVecBlockDirect', iisize, rx%h_Ddata, &
          ST_NEWBLOCK_NOINIT, .false.)
      
    elseif (present(NEQMAX)) then

      ! The available memory suffices for all componenets of the block vector.
      ! Let's check if the user supplied a new upper limit which makes it
      ! mandatory to "shrink" the allocated memory. Note that memory for at
      ! least NEQ=SUM(Isize) vector entries as allocated in any case.
      if (iisize > iNEQ) then
        
        ! Compute new size, i.e. MAX(0,NEQ,NEQMAX)
        iisize = iNEQ

        if (iisize == 0) then
          ! If nothing is left, then the vector can also be release
          call lsysbl_releaseVector(rx)
          call lsysbl_releaseVector(rxTmp)
          return
        else
          ! Reallocate the memory for handle h_Ddata
          call storage_realloc('lsysbl_resizeVecBlockDirect', iisize, rx%h_Ddata, &
              ST_NEWBLOCK_NOINIT, .false.)
        end if
      end if
    end if

    ! Should the vector be cleared?
    if (bclear) call storage_clear(rx%h_Ddata)

    ! Restore the structure of the scalar subvectors
    n=1
    do i=1,rx%nblocks

      ! Check that Isize(i) is a multiple of NVAR
      if (mod(Isize(i), rx%RvectorBlock(i)%NVAR) .ne. 0) then
        call output_line('Size of the scalar subvector is not a multiple of NVAR!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_resizeVecBlockDirect')
        call sys_halt()
      end if
      rx%RvectorBlock(i)%NEQ = Isize(i)/rx%RvectorBlock(i)%NVAR
      rx%RvectorBlock(i)%iidxFirstEntry = n
      n = n + rx%RvectorBlock(i)%NEQ

      ! Remove any sorting strategy 
      rx%RvectorBlock(i)%isortStrategy      = 0
      rx%RvectorBlock(i)%h_iSortPermutation = ST_NOHANDLE
    end do

    ! If the content should be copied use the temporal vector
    if (bdocopy) then
      select case(rx%cdataType)
      case (ST_DOUBLE)
        do i=1,rx%nblocks
          call lsyssc_getbase_double(rx%RvectorBlock(i), p_Ddata)
          call lsyssc_getbase_double(rxTmp%RvectorBlock(i), p_DdataTmp)
          n = min(size(p_Ddata), size(p_DdataTmp))
          call lalg_copyVectorDble(p_DdataTmp(1:n), p_Ddata(1:n))
        end do

      case (ST_SINGLE)
        do i=1,rx%nblocks
          call lsyssc_getbase_single(rx%RvectorBlock(i), p_Fdata)
          call lsyssc_getbase_single(rxTmp%RvectorBlock(i), p_FdataTmp)
          n = min(size(p_Fdata), size(p_FdataTmp))
          call lalg_copyVectorSngl(p_FdataTmp(1:n), p_Fdata(1:n))
        end do

      case (ST_INT)
        do i=1,rx%nblocks
          call lsyssc_getbase_int(rx%RvectorBlock(i), p_Idata)
          call lsyssc_getbase_int(rxTmp%RvectorBlock(i), p_IdataTmp)
          n = min(size(p_Idata), size(p_IdataTmp))
          call lalg_copyVectorInt(p_IdataTmp(1:n), p_Idata(1:n))
        end do

      case DEFAULT
        call output_line('Unsupported data format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_resizeVecBlockDirect')
        call sys_halt()
      end select

      ! Release temporal vector
      call lsysbl_releaseVector(rxTmp)
    end if

  end subroutine lsysbl_resizeVecBlockDirect

  !****************************************************************************

!<subroutine>

  subroutine lsysbl_resizeVecBlockDirectDims (rx, isize, bclear, bcopy, NEQMAX)

!<description>
  ! Resize the vector block structure rx. Isize is an integer value
  ! which denotes the new length of all blocks of the global block vector.
  !
  ! Remark: If the vector structure rx is not initialized, then
  ! it cannot be resized and the routine stops with an error.
  ! If it exists already, then the total memory on the heap is
  ! reallocated if "(rx%nblocks*isize) > SIZE(rx)", that is, 
  ! the new dimension of rx exceeds that of the old vectors. If the old
  ! vector is larger than the new one, then the dimensions of the
  ! subvectors are only adjusted without reallocation.
  ! If the optional parameter NEQMAX is specified, then this value
  ! is taken as upper bound for the total memory.
  !
  ! If the parameter bclear=.TRUE., then the resized vector is cleared.
  ! Otherwise it is left 'as is'.
  !
  ! Remark: In contrast to scalar vectors where the content of the vector
  ! can be copied upon resize this is not possible for block vectors.
  ! Theoretically, this would be possible. However, some subvectors may
  ! have increased whereas other may have decreased, so that an additional
  ! vector is mandatory to reorganize the global array. This is not efficient.
!</description>

!<input>
    
    ! Integer with desired length of all vector blocks
    integer, intent(IN)               :: isize

    ! Whether to fill the vector with zero initially
    logical, intent(IN)               :: bclear

    ! OPTIONAL: Whether to copy the content of the vector to the resized one
    logical, intent(IN), optional     :: bcopy

    ! OPTIONAL: Maximum length of the vector
    integer, intent(IN), optional     :: NEQMAX

!</input>

!<inputoutput>

    ! Block vector structure
    type(t_vectorBlock), intent(INOUT)             :: rx

!</inputoutput>

!</subroutine>
    
    ! local variabls
    integer, dimension(max(rx%nblocks,1)) :: Isubsize

    ! Fill auxiliary vector Iisize
    Isubsize(1:rx%nblocks) = isize

    ! Call the direct resize routine
    call lsysbl_resizeVecBlockDirect(rx, Isubsize(1:rx%nblocks), bclear, bcopy, NEQMAX)
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsysbl_resizeVecBlockIndirect (rx, rTemplate, bclear, bcopy)

!<description>
  ! Resize the vector block structure rx so that is resembles that of
  ! the template vector. If rx has not been initialized, than it is
  ! initialized adopting the same memory layout as the template vector.
  !
  ! Remark: If the vector exists already, then the total memory on 
  ! the heap is reallocated if "SIZE(rTemplate) > SIZE(rx)", that is, 
  ! the new dimension of rx exceeds that of the old vectors. If the old
  ! vector is larger than the new one, then the dimensions of the
  ! subvectors are only adjusted without reallocation.
  !
  ! If the parameter bclear=.TRUE., then the resized vector is cleared.
  ! Otherwise it is left 'as is'.
  !
  ! Remark: In contrast to scalar vectors where the content of the vector
  ! can be copied upon resize this is not possible for block vectors.
  ! Theoretically, this would be possible. However, some subvectors may
  ! have increased whereas other may have decreased, so that an additional
  ! vector is mandatory to reorganize the global array. This is not efficient.
!</description>

!<input>

    ! Template block vector structure
    type(t_vectorBlock), intent(IN)                :: rTemplate

    ! Whether to fill the vector with zero initially
    logical, intent(IN)                            :: bclear

    ! OPTIONAL: Whether to copy the content of the vector to the resized one
    logical, intent(IN), optional                  :: bcopy

!</input>

!<inputoutput>

    ! Block vector structure
    type(t_vectorBlock), intent(INOUT)             :: rx

!</inputoutput>

!</subroutine>

    ! local variabls
    integer, dimension(max(rx%nblocks,1)) :: Isize
    integer :: NEQMAX
    integer :: i

    ! Check if vector is initialized
    if (rx%NEQ == 0 .or. rx%h_Ddata == ST_NOHANDLE) then
      call lsysbl_createVectorBlock(rTemplate, rx, bclear)
      
    else
      
      ! Check if vectors are compatible
      if ((rx%cdataType .ne. rTemplate%cdataType) .or. &
          (rx%nblocks   .ne. rTemplate%nblocks  )) then
        call output_line('Vectors are incompatible!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_resizeVecBlockIndirect')
        call sys_halt()
      end if

      ! Fill auxiliary vector Iisize
      do i=1,rTemplate%nblocks

        if (rx%RvectorBlock(i)%NVAR .ne. rTemplate%RvectorBlock(i)%NVAR) then
          call output_line('Scalar subvectors are incompatible!',&
              OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_resizeVecBlockIndirect')
          call sys_halt()
        end if
        Isize(i) = rTemplate%RvectorBlock(i)%NEQ*rTemplate%RvectorBlock(i)%NVAR
      end do

      ! Get current size of global vector
      call storage_getsize(rTemplate%h_Ddata,NEQMAX)
      
      ! Resize vector directly
      call lsysbl_resizeVecBlockDirect(rx, Isize(1:rx%nblocks), bclear, bcopy, NEQMAX)
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsysbl_resizeVecBlockIndMat (rtemplateMat, rx, bclear, bcopy)

!<description>
    ! Resize the vector block structure rx so that it is compatible with
    ! the template matrix. If rx has not been initialized, than it is
    ! initialized adopting the same memory layout as the template matrix.
    !
    ! Remark: If the vector exists already, then the total memory on 
    ! the heap is reallocated if "rtemplateMat%NCOLS > SIZE(rx)", that is, 
    ! the new dimension of rx exceeds that of the old vector. If the old
    ! vector is larger than the new one, then the dimensions of the
    ! subvectors are only adjusted without reallocation.
    !
    ! If the parameter bclear=.TRUE., then the resized vector is cleared.
    ! Otherwise it is left 'as is'.
    !
    ! Remark: In contrast to scalar vectors where the content of the vector
    ! can be copied upon resize this is not possible for block vectors.
    ! Theoretically, this would be possible. However, some subvectors may
    ! have increased whereas other may have decreased, so that an additional
    ! vector is mandatory to reorganize the global array. This is not efficient.
!</description>

!<input>
    ! A template matrix structure
    type(t_matrixBlock), intent(IN)    :: rtemplateMat

    ! OPTIONAL: If set to TRUE, the vector will be filled with zero initially.
    ! Otherwise the content of rx is undefined.
    logical, intent(IN), optional      :: bclear

    ! OPTIONAL: Whether to copy the content of the vector to the resized one
    logical, intent(IN), optional      :: bcopy
!</input>

!<inputoutput>
    ! Block vector structure
    type(t_vectorBlock), intent(INOUT) :: rx
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorBlock)                 :: rxTmp
    real(DP), dimension(:), pointer     :: p_Ddata,p_DdataTmp
    real(SP), dimension(:), pointer     :: p_Fdata,p_FDataTmp
    integer, dimension(:), pointer :: p_Idata,p_IdataTmp
    integer :: isize,i,j,n
    logical              :: bdocopy

    ! Check, that the vector is not a copy of another (possibly larger) vector
    if (rx%bisCopy) then
      print *, "lsysbl_resizeVecBlockIndMat: A copied vector cannot be resized!"
      call sys_halt()
    end if

    ! Set copy/clear attributes
    bdocopy = (.not.bclear)
    if (present(bcopy)) bdocopy = (bdocopy .and. bcopy)

    ! Check if vector exists?
    if (rx%NEQ == 0 .or.&
        rx%h_DData == ST_NOHANDLE) then

      call lsysbl_createVecBlockIndMat(rtemplateMat,rx,bclear)

    else
      
      ! Check if vector/matrix are compatible
      if (rx%nblocks .ne. rtemplateMat%nblocksPerRow) then
        print *, "lsysbl_resizeVecBlockIndMat: Matrix/Vector incompatible!"
        call sys_halt()
      end if

      ! If the vector should be cleared, then the sorting strategy (if any)
      ! can be ignored and reset. Otherwise, the vector needs to be unsorted
      ! prior to copying some part of it. Afterwards, no sorting strategy is
      ! available in any case.
      if (bdocopy) then
        do i = 1, rx%nblocks
          if (rx%RvectorBlock(i)%isortStrategy > 0) then
            call lsyssc_vectorActivateSorting(rx%RvectorBlock(i),.false.)
          end if
        end do
        ! Make a copy of the unsorted content
        call lsysbl_duplicateVector(rx, rxTmp,&
            LSYSSC_DUP_COPY, LSYSSC_DUP_COPY)
      end if

      ! Get current size of vector memory
      call storage_getsize(rx%h_Ddata, isize)

      ! Update the global NEQ
      rx%NEQ = rtemplateMat%NCOLS
      
      ! Do we really have to reallocate the vector physically?
      if (rx%NEQ > isize) then

        ! Yes, reallocate memory for handle h_Ddata
        call storage_realloc('lsysbl_resizeVecBlockIndMat', rx%NEQ, rx%h_Ddata, &
            ST_NEWBLOCK_NOINIT, bdocopy)
      end if

      n=1
      do i = 1, rtemplateMat%nblocksPerRow
        ! Search for the first matrix in column i of the block matrix -
        ! this will give us information about the vector. Note that the
        ! diagonal blocks does not necessarily have to exist!
        do j = 1, rtemplateMat%nblocksPerCol
          
          ! Check if the matrix is not empty
          if (rtemplateMat%RmatrixBlock(j,i)%NCOLS .gt. 0) then
            
            ! Found a template matrix we can use :-)
            rx%RvectorBlock(i)%NEQ = rtemplateMat%RmatrixBlock(j,i)%NCOLS
            
            ! Take the handle of the complete-solution vector, but set the index of
            ! the first entry to a value >= 1 - so that it points to the first
            ! entry in the global solution vector!
            rx%RvectorBlock(i)%h_Ddata = rx%h_Ddata
            rx%RvectorBlock(i)%iidxFirstEntry = n
            
            ! Give the vector the discretisation of the matrix
            rx%RvectorBlock(i)%p_rspatialDiscr => &
                rtemplateMat%RmatrixBlock(j,i)%p_rspatialDiscrTrial
            
            ! Give the vector the same sorting strategy as the matrix, so that
            ! the matrix and vector get compatible. Otherwise, things
            ! like matrix vector multiplication won't work...
            rx%RvectorBlock(i)%isortStrategy = &
                rtemplateMat%RmatrixBlock(j,i)%isortStrategy
            rx%RvectorBlock(i)%h_IsortPermutation = &
                rtemplateMat%RmatrixBlock(j,i)%h_IsortPermutation
            
            ! Denote in the subvector that the handle belongs to us - not to
            ! the subvector.
            rx%RvectorBlock(i)%bisCopy = .true.
            
            n = n+rtemplateMat%RmatrixBlock(j,i)%NCOLS
            
            ! Finish this loop, continue with the next column
            exit
            
          end if
          
        end do

        if (j .gt. rtemplateMat%nblocksPerCol) then
          ! Let's hope this situation (an empty equation) never occurs - 
          ! might produce some errors elsewhere :)
          rx%RvectorBlock(i)%NEQ = 0
          rx%RvectorBlock(i)%iidxFirstEntry = 0
        end if

      end do
    end if

    ! Transfer the boundary conditions and block discretisation pointers
    ! from the matrix to the vector.
    rx%p_rblockDiscr => rtemplateMat%p_rblockDiscrTrial
    rx%p_rdiscreteBC     => rtemplateMat%p_rdiscreteBC
    rx%p_rdiscreteBCfict => rtemplateMat%p_rdiscreteBCfict
    
    ! Should the vector be cleared?
    if (bclear) then
      call lsysbl_clearVector (rx)
    end if

    ! If the content should be copied use the temporal vector
    if (bdocopy) then
      select case(rx%cdataType)
      case (ST_DOUBLE)
        do i=1,rx%nblocks
          call lsyssc_getbase_double(rx%RvectorBlock(i), p_Ddata)
          call lsyssc_getbase_double(rxTmp%RvectorBlock(i), p_DdataTmp)
          n = min(size(p_Ddata), size(p_DdataTmp))
          call lalg_copyVectorDble(p_DdataTmp(1:n), p_Ddata(1:n))
        end do

      case (ST_SINGLE)
        do i=1,rx%nblocks
          call lsyssc_getbase_single(rx%RvectorBlock(i), p_Fdata)
          call lsyssc_getbase_single(rxTmp%RvectorBlock(i), p_FdataTmp)
          n = min(size(p_Fdata), size(p_FdataTmp))
          call lalg_copyVectorSngl(p_FdataTmp(1:n), p_Fdata(1:n))
        end do

      case (ST_INT)
        do i=1,rx%nblocks
          call lsyssc_getbase_int(rx%RvectorBlock(i), p_Idata)
          call lsyssc_getbase_int(rxTmp%RvectorBlock(i), p_IdataTmp)
          n = min(size(p_Idata), size(p_IdataTmp))
          call lalg_copyVectorInt(p_IdataTmp(1:n), p_Idata(1:n))
        end do

      case DEFAULT
        call output_line('Unsupported data format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_resizeVecBlockIndMat')
        call sys_halt()
      end select

      ! Release temporal vector
      call lsysbl_releaseVector(rxTmp)
    end if
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsysbl_infoVector(rvector)

!<description>
    ! This subroutine outputs information about the block vector
!</description>

!<input>
    ! block vector
    type(t_vectorBlock), intent(IN) :: rvector
!</input>
!</subroutine>

    ! local variables
    integer :: iblock

    call output_lbrk()
    call output_line ('UUID: '//uuid_conv2String(rvector%ruuid))
    call output_line ('Vector is a ('&
        //trim(sys_siL(rvector%nblocks,3))//') vector.')

    do iblock = 1, rvector%nblocks
      call output_lbrk()
      call output_line ('Vector-block #'//trim(sys_siL(iblock,15)))
      call lsyssc_infoVector(rvector%RvectorBlock(iblock))
    end do

  end subroutine lsysbl_infoVector

  !****************************************************************************

!<subroutine>

  subroutine lsysbl_infoMatrix(rmatrix)

!<description>
    ! This subroutine outputs information about the block matrix
!</description>

!<input>
    ! block matrix
    type(t_matrixBlock) :: rmatrix
!</input>
!</subroutine>

    ! local variables
    integer :: iblock,jblock

    call output_lbrk()
    call output_line ('UUID: '//uuid_conv2String(rmatrix%ruuid))
    call output_line ('Matrix is a ('&
        //trim(sys_siL(rmatrix%nblocksPerCol,3))//','&
        //trim(sys_siL(rmatrix%nblocksPerRow,3))//') matrix.')

    do jblock=1,rmatrix%nblocksPerRow
      do iblock=1,rmatrix%nblocksPerCol
        if ((rmatrix%RmatrixBlock(iblock,jblock)%NEQ /= 0) .and.&
            (rmatrix%RmatrixBlock(iblock,jblock)%NCOLS /= 0)) then
          call output_lbrk()
          call output_line ('Matrix-block #('//trim(sys_siL(iblock,2))//','//&
              trim(sys_siL(jblock,2))//')')
          call output_lbrk()
          call lsyssc_infoMatrix(rmatrix%RmatrixBlock(iblock,jblock))
        end if
      end do
    end do
  end subroutine lsysbl_infoMatrix

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_convertVecFromScalar (rscalarVec,rvector,&
                                          rblockDiscretisation,nblocks)
  
!<description>
  ! This routine creates a 1-block vector rvector from a scalar vector
  ! rscalarVec. Both, rscalarVec and rvector will share the same handles,
  ! so changing the content of rvector will change rscalarVec, too.
  ! In contrast to lsysbl_createVecFromScalar, the 1-block vector rvector
  ! is marked as template vector, whereas the bisCopy flag of the original
  ! scalar vector rscalarVec is set to TRUE. This allows to release the
  ! scalar vector without releasing the 1-block vector.
!</description>
  
!<input>
  ! OPTIONAL: A block discretisation structure.
  ! A pointer to this will be saved to the vector.
  type(t_blockDiscretisation), intent(IN), optional, target :: rblockDiscretisation
  
  ! OPTIONAL: Number of blocks to reserve.
  ! Normally, the created scalar vector has only one block. If nblocks
  ! is specified, the resulting vector will have more blocks while only
  ! the first block vector is used.
  integer, intent(IN), optional :: nblocks
!</input>

!<inputoutput>
  ! The scalar vector which should provide the data
  type(t_vectorScalar), intent(INOUT) :: rscalarVec
!</inputoutput>

!<output>
  ! The block vector, created from rscalarVec.
  type(t_vectorBlock), intent(OUT) :: rvector
!</output>
  
!</subroutine>

    integer :: nactblocks
    
    nactblocks = 1
    if (present(nblocks)) nactblocks = max(nblocks,nactblocks)
    if (present(rblockDiscretisation)) &
      nactblocks = max(rblockDiscretisation%ncomponents,nactblocks)

    ! Fill the rvector structure with data.
    rvector%NEQ         = rscalarVec%NEQ * rscalarVec%NVAR
    rvector%cdataType   = rscalarVec%cdataType
    rvector%h_Ddata     = rscalarVec%h_Ddata
    rvector%nblocks     = nactblocks
    
    ! Copy the content of the scalar matrix structure into the
    ! first block of the block vector.
    allocate(rvector%RvectorBlock(nactblocks))
    rvector%RvectorBlock(1) = rscalarVec
    
    ! Copy the starting address of the scalar vector to our block vector.
    rvector%iidxFirstEntry  = rscalarVec%iidxFirstEntry

    ! Adopt the copyFlag from the template vector
    rvector%bisCopy = rscalarVec%bisCopy

    ! The handle of the scalar subvector belongs to the block vector.
    rvector%RvectorBlock(1)%bisCopy = .true.

    ! The handle of the template vector no longer belongs to it.
    rscalarVec%bisCopy = .true.
    
    ! Save a pointer to the discretisation structure if that one
    ! is specified.
    if (present(rblockDiscretisation)) then
      rvector%p_rblockDiscr => rblockDiscretisation
    else
      nullify(rvector%p_rblockDiscr)
    end if
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine lsysbl_convertMatFromScalar (rscalarMat,rmatrix,&
                                          rblockDiscrTest,rblockDiscrTrial)
  
!<description>
  ! This routine creates a 1x1 block matrix rmatrix from a scalar matrix
  ! rscalarMat. Both, rscalarMat and rmatrix will share the same handles,
  ! so changing the content of rmatrix will change rscalarMat, too.
  ! In contrast to lsysbl_createMatFromScalar, the 1x1 block matrix rmatrix
  ! is marked as template matrix, whereas the imatrixSpec flag of the original
  ! scalar matrix rscalarMat is set to LSYSSC_MSPEC_ISCOPY. This allows to 
  ! release the scalar matrix without releasing the 1x1 block matrix.
!</description>
  
!<input>
  ! OPTIONAL: A block discretisation structure specifying the test functions.
  ! A pointer to this will be saved to the matrix.
  type(t_blockDiscretisation), intent(IN), optional, target :: rblockDiscrTest

  ! OPTIONAL: A block discretisation structure specifying the trial functions.
  ! A pointer to this will be saved to the matrix.
  ! If not specified, the trial functions coincide with the test functions
  type(t_blockDiscretisation), intent(IN), optional, target :: rblockDiscrTrial
!</input>

!<inputoutput>
  ! The scalar matrix which should provide the data
  type(t_matrixScalar), intent(INOUT) :: rscalarMat
!</inputoutput>

!<output>
  ! The 1x1 block matrix, created from rscalarMat.
  type(t_matrixBlock), intent(OUT) :: rmatrix
!</output>
  
!</subroutine>

    ! Fill the rmatrix structure with data.
    rmatrix%NEQ           = rscalarMat%NEQ * rscalarMat%NVAR
    rmatrix%NCOLS         = rscalarMat%NCOLS * rscalarMat%NVAR
    rmatrix%nblocksPerCol = 1
    rmatrix%nblocksPerRow = 1
    rmatrix%imatrixSpec   = LSYSBS_MSPEC_SCALAR
    
    ! Copy the content of the scalar matrix structure into the
    ! first block of the block matrix
    allocate(rmatrix%RmatrixBlock(1,1))
    rmatrix%RmatrixBlock(1,1) = rscalarMat

    ! The handles of the template matrix no longer belong to it.
    rscalarMat%imatrixSpec = ior(rscalarMat%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
      
    ! Save a pointer to the discretisation structure if that one
    ! is specified.
    if (present(rblockDiscrTrial)) then
      rmatrix%p_rblockDiscrTrial => rblockDiscrTrial
      if (present(rblockDiscrTest)) then
        rmatrix%p_rblockDiscrTest => rblockDiscrTest
        rmatrix%bidenticalTrialAndTest = .false.
      else
        rmatrix%p_rblockDiscrTest => rblockDiscrTrial
        rmatrix%bidenticalTrialAndTest = .true.
      end if
    else
      nullify(rmatrix%p_rblockDiscrTest)
      nullify(rmatrix%p_rblockDiscrTrial)
      rmatrix%bidenticalTrialAndTest = .true.
    end if
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
 
  subroutine lsysbl_insertSubmatrix (rsourceMatrix,rdestMatrix,&
                                     cdupStructure, cdupContent,&
                                     iy,ix)
  
!<description>
  ! Inserts the block matrix as submatrix into another block matrix.
  !
  ! (iy,ix) is the upper left position in the destination matrix where
  ! rdestMatrix where rsourceMatrix should be inserted.
!</description>

  ! The source matrix to put into the destination matrix.
  type(t_matrixBlock), intent(IN) :: rsourceMatrix
  
  ! Duplication flag that decides on how to set up the structure
  ! of rdestMatrix. This duplication flag is applied to all submatrices
  ! of rsourceMatrix.
  !
  ! One of the LSYSSC_DUP_xxxx flags:
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
  integer, intent(IN)                            :: cdupStructure
  
  ! Duplication flag that decides on how to set up the content
  ! of rdestMatrix. This duplication flag is applied to all submatrices
  ! of rsourceMatrix.
  !
  ! One of the LSYSSC_DUP_xxxx flags:
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
  integer, intent(IN)                            :: cdupContent
  
  ! X- and Y-position in the destination matrix where rsourceMatrix
  ! should be put to.
  integer, intent(IN)                            :: iy,ix
!</input>

!<inputoutput>
  ! Destination matrix where the source matrix should be put into. 
  type(t_matrixBlock), intent(INOUT)             :: rdestMatrix
!</inputoutput>  
  
!</subroutine>

    ! local variables
    integer :: i,j
    
    ! loop over all columns and rows in the source matrix
    do j=1,rsourceMatrix%nblocksPerRow
      do i=1,rsourceMatrix%nblocksPerCol
        ! Copy the submatrix from the source matrix to the destination
        ! matrix.
        call lsyssc_duplicateMatrix (rsourceMatrix%RmatrixBlock(i,j),&
            rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1),&
            cdupStructure,cdupContent)
      end do
    end do

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine lsysbl_createFpdbObjectVec (rfpdbObjectItem, sname, rvector)

!<description>
    ! This subroutine creates an abstract ObjectItem from the block
    ! vector that can be stored in the persistence database
!</description>

!<input>
    ! The full qualified name of the object
    character(LEN=*), intent(IN) :: sname
!</input>

!<inputoutput>
    ! The ObjectItem  that is created
    type(t_fpdbObjectItem), intent(INOUT) :: rfpdbObjectItem

    ! The block vector
    type(t_vectorBlock), intent(INOUT), target :: rvector
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
    rfpdbObjectItem%stype = 't_vectorBlock'

    ! Allocate the array of data items: 10
    allocate(rfpdbObjectItem%p_RfpdbDataItem(10))
    
    ! Fill the array of data items: NEQ
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(1)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'NEQ'
    p_fpdbDataItem%iinteger = rvector%NEQ

    ! Fill the array of data items: h_Ddata
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(2)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'h_Ddata'
    p_fpdbDataItem%iinteger = rvector%h_Ddata

    ! Fill the array of data items: iidxFirstEntry
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(3)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'iidxFirstEntry'
    p_fpdbDataItem%iinteger = rvector%iidxFirstEntry

    ! Fill the array of data items: cdataType
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(4)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'cdataType'
    p_fpdbDataItem%iinteger = rvector%cdataType
    
    ! Fill the array of data items: bisCopy
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(5)
    p_fpdbDataItem%ctype    = FPDB_LOGICAL
    p_fpdbDataItem%sname    = 'bisCopy'
    p_fpdbDataItem%blogical = rvector%bisCopy

    ! Fill the array of data items: nblocks
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(6)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'nblocks'
    p_fpdbDataItem%iinteger = rvector%nblocks

    ! Fill the array of data items: p_rblockDiscr
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(7)
    p_fpdbDataItem%sname = 'p_rblockDiscr'
    if (associated(rvector%p_rblockDiscr)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if
    
    ! Fill the array of data items: p_rdiscreteBC
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(8)
    p_fpdbDataItem%sname = 'p_rdiscreteBC'
    if (associated(rvector%p_rdiscreteBC)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if  

    ! Fill the array of data items: p_rdiscreteBCfict
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(9)
    p_fpdbDataItem%sname = 'p_rdiscreteBCfict'
    if (associated(rvector%p_rdiscreteBCfict)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if

    ! Fill the array of data items: RvectorBlock
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(10)
    p_fpdbDataItem%sname = 'RvectorBlock'
    if (.not.associated(rvector%RvectorBlock)) then
      p_fpdbDataItem%ctype = FPDB_NULL
    else
      ! Creating the data item for the RvectorBlock is slightly more
      ! complicated. First, we need to create a new ObjectItem which
      ! is associated to the DataItem of the original block vector.
      p_fpdbDataItem%ctype = FPDB_OBJECT
      allocate(p_fpdbDataItem%p_fpdbObjectItem)
      call createFpdbObjectVectorBlock(p_fpdbDataItem%p_fpdbObjectItem)
    end if

  contains

    !**************************************************************
    ! Create a separate ObjectItem for the RvectorBlock

    subroutine createFpdbObjectVectorBlock (rfpdbObjectItem)

      ! The object item that represents the RvectorBlock
      type(t_fpdbObjectItem), intent(inout) :: rfpdbObjectItem

      ! local variables
      integer :: iblock

      ! Check if storage block has UUID; otherwise create one
      if (uuid_isNil(rfpdbObjectItem%ruuid)) then
        call uuid_createUUID(4, rfpdbObjectItem%ruuid)
      end if

      ! Set the name and type of the object
      rfpdbObjectItem%sname = 'RvectorBlock'
      rfpdbObjectItem%stype = 't_vectorScalar'
      
      ! Allocate the array of data items: size(rvector%RvectorBlock)
      allocate(rfpdbObjectItem%p_RfpdbDataItem(size(rvector%RvectorBlock)))

      ! For each scalar subvector a new ObjectItem is created and associated
      ! to the DataItem of the array of scalar vectors.
      do iblock = 1, rvector%nblocks
        rfpdbObjectItem%p_RfpdbDataItem(iblock)%ctype = FPDB_OBJECT
        rfpdbObjectItem%p_RfpdbDataItem(iblock)%sname = 't_vectorScalar'
        allocate(rfpdbObjectItem%p_RfpdbDataItem(iblock)%p_fpdbObjectItem)
        call lsyssc_createFpdbObjectVec(&
            rfpdbObjectItem%p_RfpdbDataItem(iblock)%p_fpdbObjectItem,&
            '', rvector%RvectorBlock(iblock))
      end do
      
    end subroutine createFpdbObjectVectorBlock

  end subroutine lsysbl_createFpdbObjectVec

  !************************************************************************

!<subroutine>

  subroutine lsysbl_createFpdbObjectMat (rfpdbObjectItem, sname, rmatrix)

!<description>
    ! This subroutine creates an abstract ObjectItem from the block
    ! matrix that can be stored in the persistence database
!</description>

!<input>
    ! The full qualified name of the object
    character(LEN=*), intent(IN) :: sname
!</input>

!<inputoutput>
    ! The ObjectItem  that is created
    type(t_fpdbObjectItem), intent(INOUT) :: rfpdbObjectItem

    ! The block matrix
    type(t_matrixBlock), intent(INOUT), target :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fpdbDataItem), pointer :: p_fpdbDataItem
    
    ! Check if vector has UUID; otherwise create one
    if (uuid_isNil(rmatrix%ruuid)) then
      call uuid_createUUID(4, rmatrix%ruuid)
    end if
    rfpdbObjectItem%ruuid = rmatrix%ruuid

    ! Set the name and type of the object
    rfpdbObjectItem%sname = sname
    rfpdbObjectItem%stype = 't_matrixBlock'

    ! Allocate the array of data items: 11
    allocate(rfpdbObjectItem%p_RfpdbDataItem(11))
    
    ! Fill the array of data items: NEQ
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(1)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'NEQ'
    p_fpdbDataItem%iinteger = rmatrix%NEQ

    ! Fill the array of data items: NCOLS
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(2)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'NCOLS'
    p_fpdbDataItem%iinteger = rmatrix%NCOLS

    ! Fill the array of data items: nblocksPerCol
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(3)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'nblocksPerCol'
    p_fpdbDataItem%iinteger = rmatrix%nblocksPerCol

    ! Fill the array of data items: nblocksPerRow
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(4)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'nblocksPerRow'
    p_fpdbDataItem%iinteger = rmatrix%nblocksPerRow

    ! Fill the array of data items: imatrixSpec
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(5)
    p_fpdbDataItem%ctype    = FPDB_INT
    p_fpdbDataItem%sname    = 'imatrixSpec'
    p_fpdbDataItem%iinteger = rmatrix%imatrixSpec

    ! Fill the array of data items: p_rblockDiscrTest
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(6)
    p_fpdbDataItem%sname = 'p_rblockDiscrTest'
    if (associated(rmatrix%p_rblockDiscrTest)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if

    ! Fill the array of data items: p_rblockDiscrTrial
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(7)
    p_fpdbDataItem%sname = 'p_rblockDiscrTrial'
    if (associated(rmatrix%p_rblockDiscrTrial)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if

    ! Fill the array of data items: bidenticalTrialAndTest
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(8)
    p_fpdbDataItem%ctype    = FPDB_LOGICAL
    p_fpdbDataItem%sname    = 'bidenticalTrialAndTest'
    p_fpdbDataItem%blogical = rmatrix%bidenticalTrialAndTest

    ! Fill the array of data items: p_rdiscreteBC
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(9)
    p_fpdbDataItem%sname = 'p_rdiscreteBC'
    if (associated(rmatrix%p_rdiscreteBC)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if
    
    ! Fill the array of data items: p_rdiscreteBCfict
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(10)
    p_fpdbDataItem%sname = 'p_rdiscreteBCfict'
    if (associated(rmatrix%p_rdiscreteBCfict)) then
      p_fpdbDataItem%ctype = FPDB_LINK
    else
      p_fpdbDataItem%ctype = FPDB_NULL
    end if

    ! Fill the array of data items: RmatrixBlock
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(11)
    p_fpdbDataItem%sname = 'RmatrixBlock'
    if (.not.associated(rmatrix%RmatrixBlock)) then
      p_fpdbDataItem%ctype = FPDB_NULL
    else
      ! Creating the data item for the RmatrixBlock is slightly more
      ! complicated. First, we need to create a new ObjectItem which
      ! is associated to the DataItem of the original block matrix.
      p_fpdbDataItem%ctype = FPDB_OBJECT
      allocate(p_fpdbDataItem%p_fpdbObjectItem)
      call createFpdbObjectMatrixBlock(p_fpdbDataItem%p_fpdbObjectItem)
    end if
    
  contains

    !**************************************************************
    ! Create a separate ObjectItem for the RmatrixBlock

    subroutine createFpdbObjectMatrixBlock (rfpdbObjectItem)

      ! The object item that represents the RmatrixBlock
      type(t_fpdbObjectItem), intent(inout) :: rfpdbObjectItem

      ! local variables
      integer :: i,j,idx

      ! Check if storage block has UUID; otherwise create one
      if (uuid_isNil(rfpdbObjectItem%ruuid)) then
        call uuid_createUUID(4, rfpdbObjectItem%ruuid)
      end if

      ! Set the name and type of the object
      rfpdbObjectItem%sname = 'RmatrixBlock'
      rfpdbObjectItem%stype = 't_matrixScalar'
    
      ! Allocate the array of data items
      allocate(rfpdbObjectItem%p_RfpdbDataItem(&
               size(rmatrix%RmatrixBlock,1)*&
               size(rmatrix%RmatrixBlock,2)))

      ! For each scalar submatrix a new ObjectItem is created and 
      ! associated to the DataItem of the array of scalar matrices.
      do j = 1, rmatrix%nblocksPerRow
        do i = 1, rmatrix%nblocksPerCol
          
          ! Compute index position
          idx = rmatrix%nblocksPerCol*(j-1)+i
          
          rfpdbObjectItem%p_RfpdbDataItem(idx)%ctype = FPDB_OBJECT
          rfpdbObjectItem%p_RfpdbDataItem(idx)%sname = 't_matrixScalar'
          allocate(rfpdbObjectItem%p_RfpdbDataItem(idx)%p_fpdbObjectItem)
          call lsyssc_createFpdbObjectMat(&
              rfpdbObjectItem%p_RfpdbDataItem(idx)%p_fpdbObjectItem,&
              '', rmatrix%RmatrixBlock(i,j))
        end do
      end do

    end subroutine createFpdbObjectMatrixBlock

  end subroutine lsysbl_createFpdbObjectMat

  !************************************************************************

!<subroutine>

  subroutine lsysbl_restoreFpdbObjectVec (rfpdbObjectItem, rvector)

!<description>
    ! This subroutine restores the block vector from the abstract ObjectItem 
!</description>

!<input>
    ! The object item that is created
    type(t_fpdbObjectItem), intent(IN) :: rfpdbObjectItem
!</input>

!<inputoutput>
    ! The block vector
    type(t_vectorBlock), intent(OUT), target :: rvector
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fpdbObjectItem), pointer :: p_fpdbObjectItem
    type(t_fpdbDataItem), pointer :: p_fpdbDataItem

    ! Check if ObjectItem has correct type
    if (trim(rfpdbObjectItem%stype) .ne. 't_vectorBlock') then
      call output_line ('Invalid object type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    end if

    ! Check if DataItems are associated
    if (.not.associated(rfpdbObjectItem%p_RfpdbDataItem)) then
      call output_line ('Missing data!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    end if

    ! Check if DataItems have correct size
    if (size(rfpdbObjectItem%p_RfpdbDataItem) .ne. 10) then
      call output_line ('Invalid data!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    end if

    ! Restore the UUID
    rvector%ruuid = rfpdbObjectItem%ruuid

    ! Restore the data from the DataItem: NEQ
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(1)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'NEQ')) then
      call output_line ('Invalid data: NEQ!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    else
      rvector%NEQ = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: h_Ddata
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(2)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'h_Ddata')) then
      call output_line ('Invalid data: h_Ddata!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    else
      rvector%h_Ddata = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: iidxFirstEntry
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(3)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'iidxFirstEntry')) then
      call output_line ('Invalid data: iidxFirstEntry!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    else
      rvector%iidxFirstEntry = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: cdataType
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(4)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'cdataType')) then
      call output_line ('Invalid data: cdataType!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    else
      rvector%cdataType = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: bisCopy
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(5)
    if ((p_fpdbDataItem%ctype .ne. FPDB_LOGICAL) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'bisCopy')) then
      call output_line ('Invalid data: bisCopy!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    else
      rvector%bisCopy = p_fpdbDataItem%blogical
    end if

    ! Restore the data from the DataItem: nblocks
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(6)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'nblocks')) then
      call output_line ('Invalid data: nblocks!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    else
      rvector%nblocks = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: p_rblockDiscr
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(7)
    if (trim(p_fpdbDataItem%sname) .ne. 'p_rblockDiscr') then
      call output_line ('Invalid data: p_rblockDiscr!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rvector%p_rblockDiscr)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ('Invalid data: p_rblockDiscr!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
        call sys_halt()
      end if
    end if

    ! Restore the data from the DataItem: p_rdiscreteBC
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(8)
    if (trim(p_fpdbDataItem%sname) .ne. 'p_rdiscreteBC') then
      call output_line ('Invalid data: p_rdiscreteBC!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rvector%p_rdiscreteBC)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ('Invalid data: p_rdiscreteBC!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
        call sys_halt()
      end if
    end if

    ! Restore the data from the DataItem: p_rdiscreteBCfict
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(9)
    if (trim(p_fpdbDataItem%sname) .ne. 'p_rdiscreteBCfict') then
      call output_line ('Invalid data: p_rdiscreteBCfict!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rvector%p_rdiscreteBCfict)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ('Invalid data: p_rdiscreteBCfict!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
        call sys_halt()
      end if
    end if

    ! Restore the data from the DataItem: RvectorBlock
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(10)
    if (trim(p_fpdbDataItem%sname) .ne. 'RvectorBlock') then
      call output_line ('Invalid data: RvectorBlock!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_OBJECT) then
        p_fpdbObjectItem => p_fpdbDataItem%p_fpdbObjectItem
        call restoreFpdbObjectVectorBlock(p_fpdbObjectItem)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rvector%Rvectorblock)
      else
        call output_line ('Invalid data: RvectorBlock!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectVec')
        call sys_halt()
      end if
    end if

  contains
    
    !**************************************************************
    ! Restore the RvectorBlock from the separate ObjectItem

    subroutine restoreFpdbObjectVectorBlock (rfpdbObjectItem)

      ! The object item that represents the RvectorBlock
      type(t_fpdbObjectItem), intent(inout) :: rfpdbObjectItem
      
      ! local variables
      integer :: iblock
      
      ! Check if ObjectItem has correct type
      if ((trim(rfpdbObjectItem%stype) .ne. 't_vectorScalar') .or.&
          (trim(rfpdbObjectItem%sname) .ne. 'RvectorBlock')) then
        call output_line ('Invalid object type!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectVectorBlock')
        call sys_halt()
      end if

      ! Check if DataItems are associated
      if (.not.associated(rfpdbObjectItem%p_RfpdbDataItem)) then
        call output_line ('Missing data!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectVectorBlock')
        call sys_halt()
      end if
      
      ! Check if DataItems have correct size
      if (size(rfpdbObjectItem%p_RfpdbDataItem) .ne. rvector%nblocks) then
        call output_line ('Invalid data!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectVectorBlock')
        call sys_halt()
      end if

      ! Allocate the array of scalar vectors
      allocate(rvector%RvectorBlock(rvector%nblocks))
      
      do iblock = 1, rvector%nblocks
        call lsyssc_restoreFpdbObjectVec(&
            rfpdbObjectItem%p_RfpdbDataItem(iblock)%p_fpdbObjectItem,&
            rvector%RvectorBlock(iblock))
      end do

    end subroutine restoreFpdbObjectVectorBlock

  end subroutine lsysbl_restoreFpdbObjectVec

  !************************************************************************

!<subroutine>

  subroutine lsysbl_restoreFpdbObjectMat (rfpdbObjectItem, rmatrix)

!<description>
    ! This subroutine restores the block matrix from the abstract ObjectItem 
!</description>

!<input>
    ! The object item that is created
    type(t_fpdbObjectItem), intent(IN) :: rfpdbObjectItem
!</input>

!<inputoutput>
    ! The block matrix
    type(t_matrixBlock), intent(OUT), target :: rmatrix
!</inputoutput>
!</subroutine>

    ! local variables
    type(t_fpdbObjectItem), pointer :: p_fpdbObjectItem
    type(t_fpdbDataItem), pointer :: p_fpdbDataItem

    ! Check if ObjectItem has correct type
    if (trim(rfpdbObjectItem%stype) .ne. 't_matrixBlock') then
      call output_line ('Invalid object type!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    end if

    ! Check if DataItems are associated
    if (.not.associated(rfpdbObjectItem%p_RfpdbDataItem)) then
      call output_line ('Missing data!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    end if

    ! Check if DataItems have correct size
    if (size(rfpdbObjectItem%p_RfpdbDataItem) .ne. 11) then
      call output_line ('Invalid data!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    end if

    ! Restore the UUID
    rmatrix%ruuid = rfpdbObjectItem%ruuid

    ! Restore the data from the DataItem: NEQ
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(1)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'NEQ')) then
      call output_line ('Invalid data: NEQ!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    else
      rmatrix%NEQ = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: NCOLS
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(2)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'NCOLS')) then
      call output_line ('Invalid data: NCOLS!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    else
      rmatrix%NCOLS = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: nblocksPerCol
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(3)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'nblocksPerCol')) then
      call output_line ('Invalid data: nblocksPerCol!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    else
      rmatrix%nblocksPerCol = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: nblocksPerRow
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(4)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'nblocksPerRow')) then
      call output_line ('Invalid data: nblocksPerRow!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    else
      rmatrix%nblocksPerRow = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: imatrixSpec
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(5)
    if ((p_fpdbDataItem%ctype .ne. FPDB_INT) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'imatrixSpec')) then
      call output_line ('Invalid data: imatrixSpec!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    else
      rmatrix%imatrixSpec = p_fpdbDataItem%iinteger
    end if

    ! Restore the data from the DataItem: p_rblockDiscrTest
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(6)
    if (trim(p_fpdbDataItem%sname) .ne. 'p_rblockDiscrTest') then
      call output_line ('Invalid data: p_rblockDiscrTest!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rmatrix%p_rblockDiscrTest)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ('Invalid data: p_rblockDiscrTest!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
        call sys_halt()
      end if
    end if

    ! Restore the data from the DataItem: p_rblockDiscrTrial
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(7)
    if (trim(p_fpdbDataItem%sname) .ne. 'p_rblockDiscrTrial') then
      call output_line ('Invalid data: p_rblockDiscrTrial!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rmatrix%p_rblockDiscrTrial)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ('Invalid data: p_rblockDiscrTrial!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
        call sys_halt()
      end if
    end if
    
    ! Restore the data from the DataItem: bidenticalTrialAndTest
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(8)
    if ((p_fpdbDataItem%ctype .ne. FPDB_LOGICAL) .or.&
        (trim(p_fpdbDataItem%sname) .ne. 'bidenticalTrialAndTest')) then
      call output_line ('Invalid data: bidenticalTrialAndTest!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    else
      rmatrix%bidenticalTrialAndTest = p_fpdbDataItem%blogical
    end if

    ! Restore the data from the DataItem: p_rdiscreteBC
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(9)
    if (trim(p_fpdbDataItem%sname) .ne. 'p_rdiscreteBC') then
      call output_line ('Invalid data: p_rdiscreteBC!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rmatrix%p_rdiscreteBC)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ('Invalid data: p_rdiscreteBC!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
        call sys_halt()
      end if
    end if
    
    ! Restore the data from the DataItem: p_rdiscreteBCfict
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(10)
    if (trim(p_fpdbDataItem%sname) .ne. 'p_rdiscreteBCfict') then
      call output_line ('Invalid data: p_rdiscreteBCfict!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rmatrix%p_rdiscreteBC)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_LINK) then
      else
        call output_line ('Invalid data: p_rdiscreteBCfict!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
        call sys_halt()
      end if
    end if
    
    ! Restore the data from the DataItem: RmatrixBlock
    p_fpdbDataItem => rfpdbObjectItem%p_RfpdbDataItem(11)
    if (trim(p_fpdbDataItem%sname) .ne. 'RmatrixBlock') then
      call output_line ('Invalid data: RmatrixBlock!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
      call sys_halt()
    else
      if (p_fpdbDataItem%ctype .eq. FPDB_OBJECT) then
        p_fpdbObjectItem => p_fpdbDataItem%p_fpdbObjectItem
        call restoreFpdbObjectMatrixBlock(p_fpdbObjectItem)
      elseif (p_fpdbDataItem%ctype .eq. FPDB_NULL) then
        nullify(rmatrix%RmatrixBlock)
      else
        call output_line ('Invalid data: RmatrixBlock!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_restoreFpdbObjectMat')
        call sys_halt()
      end if
    end if

    contains

      !**************************************************************
      ! Restore the RmatrixBlock from the separate ObjectItem
      
      subroutine restoreFpdbObjectMatrixBlock (rfpdbObjectItem)
        
        ! The object item that represents the RvectorBlock
        type(t_fpdbObjectItem), intent(inout) :: rfpdbObjectItem
        
        ! local variables
        integer :: i,j,idx

        ! Check if ObjectItem has correct type
        if ((trim(rfpdbObjectItem%stype) .ne. 't_matrixScalar') .or.&
            (trim(rfpdbObjectItem%sname) .ne. 'RmatrixBlock')) then
          call output_line ('Invalid object type!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectMatrixBlock')
          call sys_halt()
        end if
        
        ! Check if DataItems are associated
        if (.not.associated(rfpdbObjectItem%p_RfpdbDataItem)) then
          call output_line ('Missing data!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectMatrixBlock')
          call sys_halt()
        end if
      
        ! Check if DataItems have correct size
        if (size(rfpdbObjectItem%p_RfpdbDataItem) .ne.&
            rmatrix%nblocksPerCol*rmatrix%nblocksPerRow) then
          call output_line ('Invalid data!', &
                            OU_CLASS_ERROR,OU_MODE_STD,'restoreFpdbObjectMatrixBlock')
          call sys_halt()
        end if
        
        ! Allocate the array of scalar vectors
        allocate(rmatrix%RmatrixBlock(rmatrix%nblocksPerCol,&
                                      rmatrix%nblocksPerRow))

        do j = 1, rmatrix%nblocksPerRow
          do i = 1, rmatrix%nblocksPerCol
          
            ! Compute index position
            idx = rmatrix%nblocksPerCol*(j-1)+i
            call lsyssc_restoreFpdbObjectMat(&
                rfpdbObjectItem%p_RfpdbDataItem(idx)%p_fpdbObjectItem,&
                rmatrix%RmatrixBlock(i,j))
          end do
        end do

      end subroutine restoreFpdbObjectMatrixBlock
      
  end subroutine lsysbl_restoreFpdbObjectMat

  !****************************************************************************
  
!<subroutine>
  
  subroutine lsysbl_unshareMatrix (rmatrix,bstructure,bdata)
  
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
  type(t_matrixBlock), intent(inout)            :: rmatrix
!</inputoutput>  

!</subroutine>

    ! local variables
    integer :: i,j
    
    ! Loop over all blocks in rmatrix and 'unshare' the data.
    ! That's all.
    do j=1,rmatrix%nblocksPerRow
      do i=1,rmatrix%nblocksPerCol
        call lsyssc_unshareMatrix (rmatrix%RmatrixBlock(i,j),bstructure,bdata)
      end do
    end do    
    
  end subroutine
    
  !****************************************************************************
  
!<subroutine>
  
  subroutine lsysbl_unshareVector (rvector)
  
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
  type(t_vectorBlock), intent(inout)            :: rvector
!</output>  

!</subroutine>

    type(t_vectorBlock) :: rvector2

    ! If the vector is not the owner, make it the owner.
    ! Note that here, we cannot do a loop over the subblocks as a vector
    ! is (normally) realised as a large array with the subvectors
    ! being part of the large array. We have to treat the block vector
    ! in the same way as a scalar one...
    if (rvector%bisCopy) then
      ! Copy & unshare the data
      call lsysbl_duplicateVector(rvector,rvector2,LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      
      ! Release the original vector
      call lsysbl_releaseVector (rvector)
      
      ! Create a hard-copy of rvector.
      ! This is one of the very rare cases where we on purpose assign a matrix
      ! structure to another...
      rvector = rvector2
    end if
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>
 
  subroutine lsysbl_moveToSubmatrix (rsourceMatrix,rdestMatrix,itop,ileft,&
      bignoreScaleFactors)
  
!<description>
  ! Inserts the block matrix as submatrix into another block matrix.
  !
  ! (iy,ix) is the upper left position in the destination matrix where
  ! rdestMatrix where rsourceMatrix should be inserted.
  !
  ! The matrix information of the source matrix is moved to the destination
  ! matrix. The source matrix remains as 'shared copy' of the destination
  ! matrix and can (and should) be removed.
  !
  ! Remark: If no block discretisation is assigned to the destination matrix, 
  !  the 'local' discretisation structures in the source matrices will
  !  be transferred to the destination matrix.
  !  If a block discretisation is assigned to the destination matrix, the
  !  discretisation structures of the scalar matrices in the source matrix 
  !  is changed according to this 'guiding' block discretisation structure.
!</description>

!<input>
  ! X- and Y-position in the destination matrix where rsourceMatrix
  ! should be put to. If not specified, 1 is assumed.
  integer, intent(IN), optional :: itop,ileft

  ! OPTIONAL: Ignore the scaling factors in the source matrix.
  ! TRUE: Submatrices of the source matrix will be copied regardless of
  !   whether dscaleFactor=0.0 or not. Standard.
  ! FALSE: Submatrices of the source matrix will only be copied if
  !   dscaleFacor <> 0.0.
  logical, intent(in), optional                  :: bignoreScaleFactors
!</input>

!<inputoutput>
  ! The source matrix to put into the destination matrix.
  ! All ownership flags change to the destination matrix, the source
  ! matrix remains as 'shared copy' of the destination.
  type(t_matrixBlock), intent(INOUT) :: rsourceMatrix

  ! Destination matrix where the source matrix should be put into. 
  type(t_matrixBlock), intent(INOUT) :: rdestMatrix
!</inputoutput>  
  
!</subroutine>

    ! local variables
    integer :: i,j,ix,iy
    integer(I32) :: cflags1,cflags2
    logical :: bignoreScale
    
    bignoreScale = .true.
    if (present(bignoreScaleFactors)) bignoreScale = bignoreScaleFactors

    ix = 1
    iy = 1
    if (present(itop)) iy=itop
    if (present(ileft)) ix=ileft
    
    ! loop over all columns and rows in the source matrix
    do j=1,rsourceMatrix%nblocksPerRow
      do i=1,rsourceMatrix%nblocksPerCol
        ! Copy the submatrix from the source matrix to the destination
        ! matrix.
        if (lsysbl_isSubmatrixPresent(rsourceMatrix,i,j,bignoreScale)) then
          call lsyssc_duplicateMatrix (rsourceMatrix%RmatrixBlock(i,j),&
              rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
              
          ! Change the ownership from the source to the destination matrix
          cflags1 = iand (rsourceMatrix%RmatrixBlock(i,j)%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
          cflags2 = iand (rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1)%imatrixSpec,LSYSSC_MSPEC_ISCOPY)
          rsourceMatrix%RmatrixBlock(i,j)%imatrixSpec = &
            ior(iand(rsourceMatrix%RmatrixBlock(i,j)%imatrixSpec,&
                    not(LSYSSC_MSPEC_ISCOPY)),cflags2)
          rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1)%imatrixSpec = &
            ior(iand(rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1)%imatrixSpec,&
                    not(LSYSSC_MSPEC_ISCOPY)),cflags1)
          
          ! Reassign the discretisation structures and change the spatial
          ! discretisation structures of the submatrices according to the
          ! block discretisation if there is one.
          if (associated(rdestMatrix%p_rblockDiscrTrial) .and. &
              associated(rdestMatrix%p_rblockDiscrTest)) then
            call lsyssc_assignDiscrDirectMat (&
                rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1),&
                rdestMatrix%p_rblockDiscrTrial%RspatialDiscr(j+ix-1),&
                rdestMatrix%p_rblockDiscrTest%RspatialDiscr(i+iy-1))
                
            if (rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1)%NCOLS .ne.&
                dof_igetNDofGlob(rdestMatrix%p_rblockDiscrTrial%&
                RspatialDiscr(j+ix-1))) then
              call output_line('Matrix not compatible with discretisation (NCOLS)!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_moveToSubmatrix')
              call sys_halt()
            end if

            if (rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1)%NEQ .ne.&
                dof_igetNDofGlob(rdestMatrix%p_rblockDiscrTest%&
                RspatialDiscr(i+iy-1))) then
              call output_line('Matrix not compatible with discretisation (NEQ)!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_moveToSubmatrix')
              call sys_halt()
            end if
          end if

          ! Now, the destination matrix is the owner and the source matrix
          ! the copy...
        else
          ! Release the submatrix in the destination matrix if present.
          ! This is independent of the scale factor (we ignore it!) as
          ! we want to produce a destination matrix that looks like
          ! the source matrix without any skeletons in the cupboard
          ! that may be activated by switching the scaling factors...
          
          call lsyssc_releaseMatrix (&
              rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1))
        end if
      end do
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine lsysbl_allocEmptyMatrix (rmatrix,iclear,bignoreExisting,cdataType)
  
!<description>
  ! This routine allocates memory for the matrix entries of all existing
  ! submatrix without computing the entries. 
  ! This can be used to attach an 'empty' matrix to a matrix
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
  integer, intent(IN) :: iclear
  
  ! OPTIONAL: If set to TRUE, existing submatrices are ignored.
  ! Standard value is FALSE which stops the application with an error.
  logical, intent(in), optional :: bignoreExisting
  
  ! OPTIONAL: Data type of the matrix (ST_SINGLE, ST_DOUBLE)
  ! If not present, the standard data type cdataType-variable in the matrix
  ! is used (usually ST_DOUBLE).
  integer, intent(IN), optional :: cdataType
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixBlock), intent(INOUT) :: rmatrix
!</inputoutput>

!</subroutine>

    integer :: i,j

    ! loop over all columns and rows in the matrix and allocate
    do j=1,rmatrix%nblocksPerRow
      do i=1,rmatrix%nblocksPerCol
        if (lsysbl_isSubmatrixPresent(rmatrix,i,j)) then
          call lsyssc_allocEmptyMatrix (rmatrix%RmatrixBlock(i,j),iclear,&
              bignoreExisting,cdataType)
        end if
      end do
    end do

  end subroutine

  ! ***************************************************************************
  
!<subroutine>
 
  subroutine lsysbl_reintegrateSubmatrix (rsourceMatrix,rdestMatrix,iy,ix,&
      bignoreScaleFactors)
  
!<description>
  ! Reintegrates a submatrix to a larger block matrix.
  !
  ! When a matrix rsourceMatrix is derived from a larger matrix by using
  ! lsysbl_deriveSubmatrix, the matrix is a standalone matrix.
  ! lsysbl_reintegrateSubmatrix allows now to 'reintegrate' such a submatrix
  ! into the previous position in the larger matrix, taking care of
  ! all scaling factors, ownership stati and so on.
  !
  ! (iy,ix) is the upper left position in the destination matrix 
  ! rdestMatrix where rsourceMatrix should be reintegrated.
  !
  ! The matrix information of the source matrix is moved to the destination
  ! matrix. The source matrix remains as 'shared copy' of the destination
  ! matrix and can (and should) be removed.
!</description>

!<input>
  ! X- and Y-position in the destination matrix where rsourceMatrix
  ! should be put to.
  integer, intent(IN)                            :: iy,ix

  ! OPTIONAL: Ignore the scaling factors in the source matrix.
  ! TRUE: Submatrices of the source matrix will be copied regardless of
  !   whether dscaleFactor=0.0 or not. Standard.
  ! FALSE: Submatrices of the source matrix will only be copied if
  !   dscaleFacor <> 0.0.
  logical, intent(in), optional                  :: bignoreScaleFactors
!</input>

!<inputoutput>
  ! The source matrix to put into the destination matrix.
  ! All ownership flags change to the destination matrix, the source
  ! matrix remains as 'shared copy' of the destination.
  type(t_matrixBlock), intent(INOUT) :: rsourceMatrix

  ! Destination matrix where the source matrix should be put into. 
  type(t_matrixBlock), intent(INOUT) :: rdestMatrix
!</inputoutput>  
  
!</subroutine>

    ! local variables
    integer :: i,j
    integer :: idupflag1,idupflag2
    logical :: bignoreScale
    
    bignoreScale = .true.
    if (present(bignoreScaleFactors)) bignoreScale = bignoreScaleFactors
    
    ! loop over all columns and rows in the source matrix
    do j=1,rsourceMatrix%nblocksPerRow
      do i=1,rsourceMatrix%nblocksPerCol
        if (lsysbl_isSubmatrixPresent(rsourceMatrix,i,j,bignoreScale)) then
          ! By default, we don't have to do anything. Just if the data of the
          ! source and destination matrix do not coincide, we copy the data
          ! to the destination.
          idupflag1 = LSYSSC_DUP_IGNORE
          idupflag2 = LSYSSC_DUP_IGNORE
          
          if (.not. lsyssc_isMatrixStructureShared(&
              rsourceMatrix%RmatrixBlock(i,j),&
              rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1))) then
            idupflag1 = LSYSSC_DUP_COPYOVERWRITE
          end if

          if (.not. lsyssc_isMatrixContentShared(&
              rsourceMatrix%RmatrixBlock(i,j),&
              rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1))) then
            idupflag2 = LSYSSC_DUP_COPYOVERWRITE
          end if
          
          ! Copy -- or do nothing. 
          call lsyssc_duplicateMatrix (rsourceMatrix%RmatrixBlock(i,j),&
              rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1),idupflag1,idupflag2)
          
          ! Ensure that we copy the scale factor, that's important!
          rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1)%dscaleFactor = &
              rsourceMatrix%RmatrixBlock(i,j)%dscaleFactor
        else if (lsysbl_isSubmatrixPresent(rsourceMatrix,i,j,.true.) .and. &
                 .not. lsysbl_isSubmatrixPresent(rdestMatrix,i+iy-1,j+ix-1,.true.)) then
          ! The matrix in the source matrix is not present anymore. Mirror that
          ! by releasing the corresponding matrix in the destination.
          call lsyssc_releaseMatrix (rdestMatrix%RmatrixBlock(i+iy-1,j+ix-1))
        end if
            
      end do
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine lsysbl_getVectorMagnitude (rvectorBlock,rvectorVecMag,dumax)
  
!<description>
  ! Compute the vector magnitude. rvectorBlock identifies a vector field. The routine
  ! will calculate the length of each vector in the vector field and store it to
  ! rvectorVecMag (if present). The maximum vector length is returned in the parameter
  ! dunorm (if present). All subvectors must have the same length.
  !
  ! IMPORTANT NOTE: This routine works algebraically. It should not be used
  ! to calculate the maximum vector length of a vector field given by a
  ! finite element function (like a velocity field), as this would need the 
  ! evaluation of finite elements! This is NOT done here. Use routines from
  ! the feevaluation.f90 module instead!
!</description>
  
!<input>
  ! Block vector defining the vector field.
  type(t_vectorBlock), intent(in) :: rvectorBlock
!</input>
  
!<output>
  ! OPTIONAL: Scalar vector receiving the vector magnitude.
  type(t_vectorScalar), intent(inout), optional :: rvectorVecMag

  ! OPTIONAL: Returns the maximum vector magnitude.
  real(dp), intent(out), optional :: dumax
!</output>
  
!</subroutine>

    ! local variables
    real(dp), dimension(:), pointer :: p_Dx,p_Dy,p_Dx1,p_Dx2,p_Dx3
    type(t_vectorScalar) :: rvectorVecMagTemp
    integer :: icomp,ieq

    ! Probably nothing to do...
    if (.not. present(rvectorVecMag) .and. .not. present(dumax)) return
    
    do icomp = 2,rvectorBlock%nblocks
      if (rvectorBlock%RvectorBlock(icomp)%NEQ .ne. rvectorBlock%RvectorBlock(1)%NEQ) then
        call output_line ('Component '//trim(sys_siL(icomp,10))//' has invalid length!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'lsysbl_getVectorMagnitude')
        call sys_halt()
      end if
    end do
    
    ! Simple case: Compute the maximum vector magnitude.
    select case (rvectorBlock%nblocks)
    case (NDIM1D)
      if (present(dumax)) then
        dumax = lsyssc_vectorNorm(rvectorBlock%RvectorBlock(1),LINALG_NORMMAX)
      end if
      
      if (present(rvectorVecMag)) then
        call lsyssc_copyVector (rvectorBlock%RvectorBlock(1),rvectorVecMag)
      end if
      
      return
      
    case (NDIM2D)
      call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Dx1)
      call lsyssc_getbase_double (rvectorBlock%RvectorBlock(2),p_Dx2)
      if (present(rvectorVecMag)) then
        call lsyssc_getbase_double (rvectorVecMag,p_Dy)
        if (present(dumax)) then
          dumax = 0.0_DP
          do ieq = 1,size(p_Dy)
            p_Dy(ieq) = sqrt(p_Dx1(ieq)**2+p_Dx2(ieq)**2)
            dumax = max(dumax,p_Dy(ieq))
          end do
        else
          do ieq = 1,size(p_Dy)
            p_Dy(ieq) = sqrt(p_Dx1(ieq)**2+p_Dx2(ieq)**2)
          end do
        end if
      else
        ! dunorm must be present here because of the above IF clause...
        dumax = 0.0_DP
        do ieq = 1,size(p_Dx1)
          dumax = max(dumax,sqrt(p_Dx1(ieq)**2+p_Dx2(ieq)**2))
        end do
      end if
      
      return
      
    case (NDIM3D)
      call lsyssc_getbase_double (rvectorBlock%RvectorBlock(1),p_Dx1)
      call lsyssc_getbase_double (rvectorBlock%RvectorBlock(2),p_Dx2)
      call lsyssc_getbase_double (rvectorBlock%RvectorBlock(3),p_Dx3)
      if (present(rvectorVecMag)) then
        call lsyssc_getbase_double (rvectorVecMag,p_Dy)
        if (present(dumax)) then
          dumax = 0.0_DP
          do ieq = 1,size(p_Dy)
            p_Dy(ieq) = sqrt(p_Dx1(ieq)**2+p_Dx2(ieq)**2+p_Dx3(ieq)**2)
            dumax = max(dumax,p_Dy(ieq))
          end do
        else
          do ieq = 1,size(p_Dy)
            p_Dy(ieq) = sqrt(p_Dx1(ieq)**2+p_Dx2(ieq)**2+p_Dx3(ieq)**2)
          end do
        end if
      else
        ! dunorm must be present here because of the above IF clause...
        dumax = 0.0_DP
        do ieq = 1,size(p_Dx1)
          dumax = max(dumax,sqrt(p_Dx1(ieq)**2+p_Dx2(ieq)**2+p_Dx3(ieq)**2))
        end do
      end if
      
      return
      
    end select
    
    ! General case: Loop over the components and calculate...
    !
    ! We need a temporary vector for this operation
    if (.not. present(rvectorVecMag)) then
      call lsyssc_duplicateVector (rvectorBlock%RvectorBlock(1),rvectorVecMagTemp,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
    else
      call lsyssc_duplicateVector (rvectorVecMag,rvectorVecMagTemp,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end if
    
    ! Clear the output
    call lsyssc_clearVector (rvectorVecMagTemp)
    call lsyssc_getbase_double (rvectorVecMagTemp,p_Dy)
    
    ! Sum up the contributions
    do icomp = 1,rvectorBlock%nblocks
    
      call lsyssc_getbase_double (rvectorBlock%RvectorBlock(icomp),p_Dx)
      
      ! For each component compute sqrt(Dx1^2 + Dx2^2 + ...)
      do ieq = 1,size(p_Dy)
        p_Dy(ieq) = p_Dy(ieq) + p_Dx(ieq)**2
      end do

    end do

    dumax = 0.0_DP
    if (present(dumax)) then
      do ieq = 1,size(p_Dy)
        p_Dy(ieq) = sqrt(p_Dy(ieq))
        dumax = max(dumax,p_Dy(ieq))
      end do
    else
      do ieq = 1,size(p_Dy)
        p_Dy(ieq) = sqrt(p_Dy(ieq))
      end do
    end if
    
    ! Release the temp vector
    call lsyssc_releaseVector(rvectorVecMagTemp)

  end subroutine

end module

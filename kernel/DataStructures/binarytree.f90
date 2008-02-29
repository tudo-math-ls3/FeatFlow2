!##############################################################################
!# ****************************************************************************
!# <name> binarytree </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements a (linear) binary/AVL tree implemented as an array
!#
!# The following routines are available:
!#
!# 1.) btree_createTree
!#     -> Create an empty tree
!#
!# 2.) btree_releaseTree
!#     -> Release an existing tree
!#
!# 3.) btree_resizeTree
!#     -> Reallocate memory for an existing tree
!#
!# 4.) btree_copyToTree = btree_copyToTree_handle /
!#                        btree_copyToTree_arrayDble /
!#                        btree_copyToTree_arraySngl / 
!#                        btree_copyToTree_arrayInt
!#     -> Copy key and auxiliary data to tree
!#
!# 5.) btree_copyFromTreeKey = btree_copyFromTreeKey_handle /
!#                             btree_copyFromTreeKey_arrayDble /
!#                             btree_copyFromTreeKey_arraySngl /
!#                             btree_copyFromTreeKey_arrayInt
!#     -> Copy the key of the tree to handle/array
!#
!# 6.) btree_copyFromTreeDble = btree_copyFromTree_handle /
!#                              btree_copyFromTree_arrayDble
!#                              btree_copyFromTree_arrayDble2D
!#     -> Copy some auxiliary double data of the tree to handle/array
!#
!# 7.) btree_copyFromTreeSngl = btree_copyFromTree_handle /
!#                              btree_copyFromTree_arraySngl
!#                              btree_copyFromTree_arraySngl2D
!#     -> Copy some auxiliary single data of the tree to handle/array
!#
!# 8.) btree_copyFromTreeInt = btree_copyFromTree_handle /
!#                             btree_copyFromTree_arrayInt
!#                             btree_copyFromTree_arrayInt2D
!#     -> Copy some auxiliary integer data of the tree to handle/array
!#
!# 9.) btree_insertToTree = btree_insertToTreeDble /
!#                          btree_insertToTreeSngl /
!#                          btree_insertToTreeInt
!#     -> Insert key into tree
!#
!# 10.) btree_deleteFromTree = btree_deleteFromTreeDble /
!#                             btree_deleteFromTreeSngl /
!#                             btree_deleteFromTreeInt
!#      -> Delete key from tree
!#
!# 11.) btree_searchInTree = btree_searchInTreeDble /
!#                           btree_searchInTreeSngl /
!#                           btree_searchInTreeInt
!#      -> Search for key in tree
!#
!# 12.) btree_getItemInTree = btree_getItemInTreeDble /
!#                            btree_getItemInTreeSngl / 
!#                            btree_getItemInTreeInt
!#      -> Search for key in tree and return position of item directly
!#
!# 13.) btree_printTree
!#      -> Print out tree
!#
!# 14.) btree_getHeight
!#      -> Get height of the tree
!#
!# 15.) btree_infoTree
!#      -> Output statistical info about the tree
!#
!# 16.) btree_duplicateTree
!#      -> Create a duplicate / backup of a binary tree
!#
!# 17.) btree_restoreTree
!#      -> Restore a binary tree from a previous backup
!#
!# </purpose>
!##############################################################################
MODULE binarytree
  USE fsystem
  USE storage
  USE genoutput
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_btree
  PUBLIC :: btree_createTree
  PUBLIC :: btree_releaseTree
  PUBLIC :: btree_resizeTree
  PUBLIC :: btree_copyToTree
  PUBLIC :: btree_copyFromTreeKey
  PUBLIC :: btree_copyFromTreeDble
  PUBLIC :: btree_copyFromTreeSngl
  PUBLIC :: btree_copyFromTreeInt
  PUBLIC :: btree_insertIntoTree
  PUBLIC :: btree_deleteFromTree
  PUBLIC :: btree_searchInTree
  PUBLIC :: btree_getItemInTree
  PUBLIC :: btree_printTree
  PUBLIC :: btree_getHeight
  PUBLIC :: btree_infoTree
  PUBLIC :: btree_duplicateTree
  PUBLIC :: btree_restoreTree

!<constants>

!<constantblock description="KIND values for tree data">

  ! kind value for indices in tree
  INTEGER, PARAMETER, PUBLIC :: PREC_TREEIDX   = I32

!</constantblock>

!<constantblock description="Global flags for tree output">

  ! Tag for preorder traversal
  INTEGER, PARAMETER, PUBLIC :: BTREE_PREORDER  = 0

  ! Tag for inorder traversal
  INTEGER, PARAMETER, PUBLIC :: BTREE_INORDER   = 1

  ! Tag for postorder traversal
  INTEGER, PARAMETER, PUBLIC :: BTREE_POSTORDER = 2

!</constantblock>

!<constantblock description="Global flags for tree operations">

  ! Identifier for "not found in tree"
  INTEGER, PARAMETER, PUBLIC :: BTREE_NOT_FOUND = -1

  ! Identifier for "found in tree"
  INTEGER, PARAMETER, PUBLIC :: BTREE_FOUND     =  0
  
  ! Identifier for "copy to tree"
  INTEGER, PARAMETER, PUBLIC :: BTREE_COPY1     =  1

  ! Identifier for "copy from tree"
  INTEGER, PARAMETER, PUBLIC :: BTREE_COPY2     =  2

!</constantblock>

!<constantblock description="Internal tags for list status">

  ! Tag for empty tree
  INTEGER(PREC_TREEIDX), PARAMETER, PUBLIC :: TNULL = 0

  ! Tag for root of tree
  INTEGER(PREC_TREEIDX), PARAMETER, PUBLIC :: TROOT = 0
  
  ! Tag for next free item
  INTEGER(PREC_TREEIDX), PARAMETER :: TFREE = 0

  ! Tag for left child
  INTEGER, PARAMETER, PUBLIC :: TLEFT       = 0

  ! Tag for right child
  INTEGER, PARAMETER, PUBLIC :: TRIGHT      = 1

!</constantblock>
!</constants>

!<types>
!<typeblock>
  
  TYPE t_btree
    ! Format-Tag:
    INTEGER :: ctreeFormat         = ST_NOHANDLE
    
    ! Current depth of the tree
    INTEGER(PREC_TREEIDX) :: depth = 0

    ! Number of items that are currently stored in the tree
    INTEGER(PREC_TREEIDX) :: NA    = 0

    ! Total number of items that can be stored in the tree
    INTEGER(PREC_TREEIDX) :: NNA   = 0

    ! Total number of items that can initially be stored in the tree
    ! This information is needed to compute the growth of the tree
    ! after several resize operations
    INTEGER(PREC_TREEIDX) :: NNA0  = 0

    ! Total number of resize operations
    INTEGER :: NRESIZE             = 0
    
    ! Dimension of the auxiliary Integer values to be stored
    INTEGER :: isizeInt            = 0

    ! Dimension of the auxiliary Double values to be stored
    INTEGER :: isizeDble           = 0

    ! Dimension of the auxiliary Single values to be stored
    INTEGER :: isizeSngl           = 0

    ! Factor by which the list is enlarged if new storage is allocate
    REAL(DP) :: dfactor            = 1.5_DP

    ! Handle to the tree key data
    INTEGER :: h_Key               = ST_NOHANDLE

    ! Handle to the balance data
    INTEGER :: h_Kbal              = ST_NOHANDLE

    ! Handle to the path data
    INTEGER :: h_Kpath             = ST_NOHANDLE

    ! Handle to the children's data
    INTEGER :: h_Kchild            = ST_NOHANDLE

    ! Handle to the tree auxiliary Integer data
    INTEGER :: h_IData             = ST_NOHANDLE

    ! Handle to the tree auxiliary Double data
    INTEGER :: h_DData             = ST_NOHANDLE

    ! Handle to the tree auxiliary Single data
    INTEGER :: h_FData             = ST_NOHANDLE

    ! Tree child structure
    ! NOTE: This array is introduced to increase performance. It
    ! should not be touched by the user. However, if the handle would
    ! be dereferenced for each operation such as search, delete,
    ! performance would be very poor.
    INTEGER(PREC_TREEIDX), DIMENSION(:,:), POINTER :: p_Kchild => NULL()

    ! Tree path structure
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_TREEIDX), DIMENSION(:), POINTER ::   p_Kpath => NULL()

    ! Tree balance structure
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER, DIMENSION(:), POINTER ::                 p_Kbal => NULL()

    ! Tree key data (Double)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(DP), DIMENSION(:), POINTER ::                p_DKey => NULL()

    ! Tree key data (Single)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(SP), DIMENSION(:), POINTER ::                p_FKey => NULL()

    ! Tree key data (Integer)
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_TREEIDX), DIMENSION(:), POINTER ::   p_IKey => NULL()

    ! Tree data (Double)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(DP), DIMENSION(:,:), POINTER ::              p_DData => NULL()

    ! Tree data (Single)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(SP), DIMENSION(:,:), POINTER ::              p_FData => NULL()

    ! Tree data (Integer)
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_TREEIDX), DIMENSION(:,:), POINTER :: p_IData => NULL()
  END TYPE t_btree

!</typeblock>
!</types>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************
      
  INTERFACE btree_copyToTree
    MODULE PROCEDURE btree_copyToTree_handle
    MODULE PROCEDURE btree_copyToTree_arrayDble
    MODULE PROCEDURE btree_copyToTree_arraySngl
    MODULE PROCEDURE btree_copyToTree_arrayInt
  END INTERFACE

  INTERFACE btree_copyFromTreeKey
    MODULE PROCEDURE btree_copyFromTreeKey_handle
    MODULE PROCEDURE btree_copyFromTreeKey_arrayDble
    MODULE PROCEDURE btree_copyFromTreeKey_arraySngl
    MODULE PROCEDURE btree_copyFromTreeKey_arrayInt
  END INTERFACE

  INTERFACE btree_copyFromTreeDble
    MODULE PROCEDURE btree_copyFromTree_handle
    MODULE PROCEDURE btree_copyFromTree_arrayDble
    MODULE PROCEDURE btree_copyFromTree_arrayDble2D
  END INTERFACE

  INTERFACE btree_copyFromTreeSngl
    MODULE PROCEDURE btree_copyFromTree_handle
    MODULE PROCEDURE btree_copyFromTree_arraySngl
    MODULE PROCEDURE btree_copyFromTree_arraySngl2D
  END INTERFACE

  INTERFACE btree_copyFromTreeInt
    MODULE PROCEDURE btree_copyFromTree_handle
    MODULE PROCEDURE btree_copyFromTree_arrayInt
    MODULE PROCEDURE btree_copyFromTree_arrayInt2D
  END INTERFACE
  
  INTERFACE btree_insertIntoTree
    MODULE PROCEDURE btree_insertToTreeDble
    MODULE PROCEDURE btree_insertToTreeSngl
    MODULE PROCEDURE btree_insertToTreeInt
  END INTERFACE

  INTERFACE btree_deleteFromTree
    MODULE PROCEDURE btree_deleteFromTreeDble
    MODULE PROCEDURE btree_deleteFromTreeSngl
    MODULE PROCEDURE btree_deleteFromTreeInt
  END INTERFACE
  
  INTERFACE btree_searchInTree
    MODULE PROCEDURE btree_searchInTreeDble
    MODULE PROCEDURE btree_searchInTreeSngl
    MODULE PROCEDURE btree_searchInTreeInt
  END INTERFACE

  INTERFACE btree_getItemInTree
    MODULE PROCEDURE btree_getItemInTreeDble
    MODULE PROCEDURE btree_getItemInTreeSngl
    MODULE PROCEDURE btree_getItemInTreeInt
  END INTERFACE  

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_createTree(rtree,nna,ctreeFormat,isizeDble,isizeSngl,isizeInt,dfactor)

!<description>
    ! This subroutine creates a new tree
!</description>

!<input>
    ! Total number of items that can be stored in the tree
    INTEGER(PREC_TREEIDX), INTENT(IN) :: nna

    ! Format-tag. Type of tree format (Double,Single,Integer)
    INTEGER, INTENT(IN) :: ctreeFormat

    ! Dimension of the auxiliary Double values to be stored
    INTEGER, INTENT(IN) :: isizeDble

    ! Dimension of the auxiliary Single values to be stored
    INTEGER, INTENT(IN) :: isizeSngl

    ! Dimension of the auxiliary Integer values to be stored
    INTEGER, INTENT(IN) :: isizeInt

    ! OPTIONAL: Factor by which the list should be enlarged if memory
    ! needs to be reallocated
    REAL(DP), INTENT(IN), OPTIONAL :: dfactor
!</input>

!<output>
    ! tree
    TYPE(t_btree), INTENT(OUT) :: rtree
!</output>
!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(2) :: Isize1,Isize2
    
    ! Set factor
    IF (PRESENT(dfactor)) THEN
      IF (dfactor > 1_DP) rtree%dfactor=dfactor
    END IF

    ! Set tree format
    rtree%ctreeFormat=ctreeFormat

    ! Initialize tree
    rtree%isizeInt  = isizeInt
    rtree%isizeDble = isizeDble
    rtree%isizeSngl = isizeSngl
    rtree%depth     = 0
    rtree%NA        = 0
    rtree%NNA       = nna
    rtree%NNA0      = nna
    rtree%NRESIZE   = 0

    ! Allocate memory for tree structure
    Isize1 = (/TLEFT,TROOT/)
    Isize2 = (/TRIGHT,nna/)
    CALL storage_new('btree_createTree','Kchild',Isize1,Isize2,&
        ST_INT,rtree%h_Kchild,ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int2D(rtree%h_Kchild,rtree%p_Kchild)
    
    CALL storage_new('btree_createTree','Kpath',CEILING(1.441*LOG2(nna))+1,&
        ST_INT,rtree%h_Kpath,ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int(rtree%h_Kpath,rtree%p_Kpath)

    CALL storage_new('btree_createTree','Kbal',nna,&
        ST_INT,rtree%h_Kbal,ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int(rtree%h_Kbal,rtree%p_Kbal)

    ! Allocate memory for Key
    SELECT CASE(rtree%ctreeFormat)
    CASE (ST_DOUBLE)
      CALL storage_new('btree_createTree','Key',nna,&
          ST_DOUBLE,rtree%h_Key,ST_NEWBLOCK_ZERO)
      CALL storage_getbase_double(rtree%h_Key,rtree%p_DKey)
    CASE (ST_SINGLE)
      CALL storage_new('btree_createTree','Key',nna,&
          ST_SINGLE,rtree%h_Key,ST_NEWBLOCK_ZERO)
      CALL storage_getbase_single(rtree%h_Key,rtree%p_FKey)
    CASE (ST_INT)
      CALL storage_new('btree_createTree','Key',nna,&
          ST_INT,rtree%h_Key,ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(rtree%h_Key,rtree%p_IKey)
    CASE DEFAULT
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_createTree')
      CALL sys_halt()
    END SELECT

    ! Allocate memory fo auxiliary data
    IF (isizeDble > 0) THEN
      Isize1 = (/isizeDble,nna/)
      CALL storage_new('btree_createTree','DData',Isize1,&
          ST_DOUBLE,rtree%h_DData,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double2D(rtree%h_DData,rtree%p_DData)
    END IF

    IF (isizeSngl > 0) THEN
      Isize1 = (/isizeSngl,nna/)
      CALL storage_new('btree_createTree','FData',Isize1,&
          ST_SINGLE,rtree%h_FData,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_single2D(rtree%h_FData,rtree%p_FData)
    END IF

    IF (isizeInt > 0) THEN
      Isize1 = (/isizeInt,nna/)
      CALL storage_new('btree_createTree','IData',Isize1,&
          ST_INT,rtree%h_IData,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_int2D(rtree%h_IData,rtree%p_IData)
    END IF

    ! Initialize tree structure
    rtree%p_Kchild(TFREE,TROOT) = 1
  END SUBROUTINE btree_createTree

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE btree_releaseTree(rtree)

!<description>
    ! This subroutine releases an existing tree
!</description>

!<inputoutput>
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>
!</subroutine>

    ! Release memory
    IF (rtree%h_Key    .NE. ST_NOHANDLE) CALL storage_free(rtree%h_Key)
    IF (rtree%h_Kbal   .NE. ST_NOHANDLE) CALL storage_free(rtree%h_Kbal)
    IF (rtree%h_Kchild .NE. ST_NOHANDLE) CALL storage_free(rtree%h_Kchild)
    IF (rtree%h_Kpath  .NE. ST_NOHANDLE) CALL storage_free(rtree%h_Kpath)

    IF (rtree%h_DData  .NE. ST_NOHANDLE) CALL storage_free(rtree%h_DDATA)
    IF (rtree%h_FData  .NE. ST_NOHANDLE) CALL storage_free(rtree%h_FDATA)
    IF (rtree%h_IData  .NE. ST_NOHANDLE) CALL storage_free(rtree%h_IDATA)

    ! Reset tree
    rtree%ctreeFormat = ST_NOHANDLE
    rtree%isizeInt    = 0
    rtree%isizeDble   = 0
    rtree%isizeSngl   = 0
    rtree%depth       = 0
    rtree%NA          = 0
    rtree%NNA         = 0
    rtree%NNA         = 0
    rtree%NRESIZE     = 0

    NULLIFY(rtree%p_DKey,rtree%p_FKey,rtree%p_IKey)
    NULLIFY(rtree%p_Kbal,rtree%p_Kpath,rtree%p_Kchild)
    NULLIFY(rtree%p_DData,rtree%p_FData,rtree%p_IData)
  END SUBROUTINE btree_releaseTree

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE btree_resizeTree(rtree,nna)

!<description>
    ! This subroutine reallocates memory for an existing tree
!</description>

!<input>
    ! New number of total items that can be stored in the tree
    INTEGER(PREC_TREEIDX), INTENT(IN) :: nna
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>
!</subroutine>

    ! Set new size
    rtree%NNA=nna
    rtree%NRESIZE=rtree%NRESIZE+1

    CALL storage_realloc('btree_resizeTree',TROOT,nna,&
        rtree%h_Kchild,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_getbase_int2D(rtree%h_Kchild,rtree%p_Kchild)

    CALL storage_realloc('btree_resizeTree',nna,&
        rtree%h_Kbal,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_getbase_int(rtree%h_Kbal,rtree%p_Kbal)

    CALL storage_realloc('btree_resizeTree',&
        CEILING(1.441*LOG2(nna))+1,rtree%h_Kpath,&
        ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_getbase_int(rtree%h_Kpath,rtree%p_Kpath)

    CALL storage_realloc('btree_resizeTree',nna,&
        rtree%h_Key,ST_NEWBLOCK_ZERO,.TRUE.)
    SELECT CASE(rtree%ctreeFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(rtree%h_Key,rtree%p_DKey)
    CASE (ST_SINGLE)
      CALL storage_getbase_single(rtree%h_Key,rtree%p_FKey)
    CASE (ST_INT)
      CALL storage_getbase_int(rtree%h_Key,rtree%p_IKey)
    CASE DEFAULT
      CALL output_line('Unsupported key data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_resizeTree')
      CALL sys_halt()
    END SELECT

    ! Auxiliary data
    IF (rtree%isizeDble > 0) THEN
      CALL storage_realloc('btree_resizeTree',nna,&
          rtree%h_DData,ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_getbase_double2D(rtree%h_DData,rtree%p_DData)
    END IF

    IF (rtree%isizeSngl > 0) THEN
      CALL storage_realloc('btree_resizeTree',nna,&
          rtree%h_FData,ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_getbase_single2D(rtree%h_FData,rtree%p_FData)
    END IF

    IF (rtree%isizeInt > 0) THEN
      CALL storage_realloc('btree_resizeTree',nna,&
          rtree%h_IData,ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_getbase_int2D(rtree%h_IData,rtree%p_IData)
    END IF
  END SUBROUTINE btree_resizeTree

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyToTree_handle(h_Key,rtree,h_DData,h_FData,h_IData)

!<description>
    ! This subroutine copies the content of handles to the tree.
    ! The key handle is mandatory to build up the tree data
    ! structure. The optional handles h_XData can be used to attach
    ! auxiliary data to the tree which is "carried" through all operations.
!</description>

!<input>
    ! handle to the key
    INTEGER, INTENT(IN) :: h_Key

    ! OPTIONAL: handle to auxiliary Double data
    INTEGER, INTENT(IN), OPTIONAL :: h_DData

    ! OPTIONAL: handle to auxiliary Single data
    INTEGER, INTENT(IN), OPTIONAL :: h_FData

    ! OPTIONAL: handle to auxiliary Integer data
    INTEGER, INTENT(IN), OPTIONAL :: h_IData
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:), POINTER   :: p_DKey
    REAL(SP), DIMENSION(:), POINTER   :: p_FKey
    INTEGER,  DIMENSION(:), POINTER   :: p_IKey
    INTEGER(PREC_TREEIDX) :: j

    ! Transform the content of the key handle to the tree
    SELECT CASE (rtree%ctreeFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(h_Key,p_DKey)
      DO j=1,SIZE(p_DKey)
        CALL btree_insertIntoTree(rtree,p_DKey(j))
      END DO
      
    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_Key,p_FKey)
      DO j=1,SIZE(p_FKey)
        CALL btree_insertIntoTree(rtree,p_FKey(j))
      END DO
      
    CASE (ST_INT)
      CALL storage_getbase_int(h_Key,p_IKey)
      DO j=1,SIZE(p_IKey)
        CALL btree_insertIntoTree(rtree,p_IKey(j))
      END DO
      
    CASE DEFAULT
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_handle')
      CALL sys_halt()
    END SELECT
    
    ! Transform the content of the auxiliary handles to the tree.
    ! Note that the sorting key was inserted linearly. In principle,
    ! only the tree data structure needs to be generated but the
    ! actual key values are still stored in sequential order. Hence,
    ! it suffices to "copy" the handles of the auxiliary data to the
    ! handles of the tree structure
    
    IF (PRESENT(h_DData) .AND. rtree%isizeDble > 0) &
        CALL storage_copy(h_DData,rtree%h_DData)
    IF (PRESENT(h_FData) .AND. rtree%isizeSngl > 0) &
        CALL storage_copy(h_FData,rtree%h_FData)
    IF (PRESENT(h_IData) .AND. rtree%isizeInt  > 0) &
        CALL storage_copy(h_IData,rtree%h_IData)
  END SUBROUTINE btree_copyToTree_handle

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyToTree_arrayDble(p_DKey,rtree,p_DData,p_FData,p_IData)

!<description>
    ! This subroutine copies the content of arrays to the tree
    ! and makes use of a double-valued key. The optional arrays
    ! p_XData can be used to attach auxiliary data to the tree
    ! which is "carried" through all operations.
!</description>

!<input>
    ! double-valued key
    REAL(DP), DIMENSION(:), INTENT(IN) :: p_DKey

    ! OPTIONAL: auxiliary Double data
    REAL(DP), DIMENSION(:,:), OPTIONAL :: p_DData

    ! OPTIONAL: auxiliary Single data
    REAL(SP), DIMENSION(:,:), OPTIONAL :: p_FData

    ! OPTIONAL: auxiliary Integer data
    INTEGER,  DIMENSION(:,:), OPTIONAL :: p_IData
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_TREEIDX) :: j

    ! Check if tree has correct data format
    IF (rtree%ctreeFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Invalid data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayDble')
      CALL sys_halt()
    END IF

    ! Build-up treee structure
    DO j=1,SIZE(p_DKey)
      CALL btree_insertIntoTree(rtree,p_DKey(j))
    END DO

    ! Copy auxiliary arrays to the tree.
    IF (PRESENT(p_DData) .AND. rtree%isizeDble > 0) THEN
      CALL storage_getsize(rtree%h_DData,Isize)
      IF ((rtree%isizeDble .NE. SIZE(p_DData,1)) .OR.&
          (rtree%NA .NE. SIZE(p_DData,2))) THEN
        CALL output_line('Invalid auxiliary data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayDble')
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('btree_copyToTree_arrayDble',&
            rtree%NNA,rtree%h_DData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_double2D(rtree%h_DData,rtree%p_DData)
      END IF
      ! Copy data
      rtree%p_DData(:,1:rtree%NA) = p_DData
    END IF

    IF (PRESENT(p_FData) .AND. rtree%isizeSngl > 0) THEN
      CALL storage_getsize(rtree%h_FData,Isize)
      IF ((rtree%isizeSngl .NE. SIZE(p_FData,1)) .OR.&
          (rtree%NA .NE. SIZE(p_FData,2))) THEN
        CALL output_line('Invalid auxiliary data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayDble')
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('btree_copyToTree_arrayDble',&
            rtree%NNA,rtree%h_FData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_single2D(rtree%h_FData,rtree%p_FData)
      END IF
      ! Copy data
      rtree%p_FData(:,1:rtree%NA) = p_FData
    END IF

    IF (PRESENT(p_IData) .AND. rtree%isizeInt > 0) THEN
      CALL storage_getsize(rtree%h_IData,Isize)
      IF ((rtree%isizeInt .NE. SIZE(p_IData,1)) .OR.&
          (rtree%NA .NE. SIZE(p_IData,2))) THEN
        CALL output_line('Invalid auxiliary data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayDble')
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('btree_copyToTree_arrayDble',&
            rtree%NNA,rtree%h_IData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_int2D(rtree%h_IData,rtree%p_IData)
      END IF
      ! Copy data
      rtree%p_IData(:,1:rtree%NA) = p_IData
    END IF
    
  END SUBROUTINE btree_copyToTree_arrayDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyToTree_arraySngl(p_FKey,rtree,p_DData,p_FData,p_IData)

!<description>
    ! This subroutine copies the content of arrays to the tree
    ! and makes use of a single-valued key. The optional arrays
    ! p_XData can be used to attach auxiliary data to the tree
    ! which is "carried" through all operations.
!</description>

!<input>
    ! single-valued key
    REAL(SP), DIMENSION(:), INTENT(IN) :: p_FKey

    ! OPTIONAL: auxiliary Double data
    REAL(DP), DIMENSION(:,:), OPTIONAL :: p_DData

    ! OPTIONAL: auxiliary Single data
    REAL(SP), DIMENSION(:,:), OPTIONAL :: p_FData

    ! OPTIONAL: auxiliary Integer data
    INTEGER,  DIMENSION(:,:), OPTIONAL :: p_IData
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_TREEIDX) :: j

    ! Check if tree has correct data format
    IF (rtree%ctreeFormat .NE. ST_SINGLE) THEN
      CALL output_line('Invalid data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arraySngl')
      CALL sys_halt()
    END IF

    ! Build-up treee structure
    DO j=1,SIZE(p_FKey)
      CALL btree_insertIntoTree(rtree,p_FKey(j))
    END DO

    ! Copy auxiliary arrays to the tree.
    IF (PRESENT(p_DData) .AND. rtree%isizeDble > 0) THEN
      CALL storage_getsize(rtree%h_DData,Isize)
      IF ((rtree%isizeDble .NE. SIZE(p_DData,1)) .OR.&
          (rtree%NA .NE. SIZE(p_DData,2))) THEN
        CALL output_line('Invalid auxiliary data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arraySngl')
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('btree_copyToTree_arraySngl',&
            rtree%NNA,rtree%h_DData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_double2D(rtree%h_DData,rtree%p_DData)
      END IF
      ! Copy data
      rtree%p_DData(:,1:rtree%NA) = p_DData
    END IF
    
    IF (PRESENT(p_FData) .AND. rtree%isizeSngl > 0) THEN
      CALL storage_getsize(rtree%h_FData,Isize)
      IF ((rtree%isizeSngl .NE. SIZE(p_FData,1)) .OR.&
          (rtree%NA .NE. SIZE(p_FData,2))) THEN
        CALL output_line('Invalid auxiliary data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arraySngl')
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('btree_copyToTree_arraySngl',&
            rtree%NNA,rtree%h_FData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_single2D(rtree%h_FData,rtree%p_FData)
      END IF
      ! Copy data
      rtree%p_FData(:,1:rtree%NA) = p_FData
    END IF

    IF (PRESENT(p_IData) .AND. rtree%isizeInt > 0) THEN
      CALL storage_getsize(rtree%h_IData,Isize)
      IF ((rtree%isizeInt .NE. SIZE(p_IData,1)) .OR.&
          (rtree%NA .NE. SIZE(p_IData,2))) THEN
        CALL output_line('Invalid auxiliary data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arraySngl')
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('btree_copyToTree_arraySngl',&
            rtree%NNA,rtree%h_IData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_int2D(rtree%h_IData,rtree%p_IData)
      END IF
      ! Copy data
      rtree%p_IData(:,1:rtree%NA) = p_IData
    END IF

  END SUBROUTINE btree_copyToTree_arraySngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyToTree_arrayInt(p_IKey,rtree,p_DData,p_FData,p_IData)

!<description>
    ! This subroutine copies the content of arrays to the tree
    ! and makes use of an integer-valued key. The optional arrays
    ! p_XData can be used to attach auxiliary data to the tree
    ! which is "carried" through all operations.
!</description>

!<input>
    ! integer-valued key
    INTEGER, DIMENSION(:), INTENT(IN) :: p_IKey

    ! OPTIONAL: auxiliary Double data
    REAL(DP), DIMENSION(:,:), OPTIONAL :: p_DData

    ! OPTIONAL: auxiliary Single data
    REAL(SP), DIMENSION(:,:), OPTIONAL :: p_FData

    ! OPTIONAL: auxiliary Integer data
    INTEGER,  DIMENSION(:,:), OPTIONAL :: p_IData
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_TREEIDX) :: j

    ! Check if tree has correct data format
    IF (rtree%ctreeFormat .NE. ST_INT) THEN
      CALL output_line('Invalid data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayInt')
      CALL sys_halt()
    END IF

    ! Build-up treee structure
    DO j=1,SIZE(p_IKey)
      CALL btree_insertIntoTree(rtree,p_IKey(j))
    END DO

    ! Copy auxiliary arrays to the tree.
    IF (PRESENT(p_DData) .AND. rtree%isizeDble > 0) THEN
      CALL storage_getsize(rtree%h_DData,Isize)
      IF ((rtree%isizeDble .NE. SIZE(p_DData,1)) .OR.&
          (rtree%NA .NE. SIZE(p_DData,2))) THEN
        CALL output_line('Invalid auxiliary data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayInt')
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('btree_copyToTree_arrayInt',&
            rtree%NNA,rtree%h_DData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_double2D(rtree%h_DData,rtree%p_DData)
      END IF
      ! Copy data
      rtree%p_DData(:,1:rtree%NA) = p_DData
    END IF
    
    IF (PRESENT(p_FData) .AND. rtree%isizeSngl > 0) THEN
      CALL storage_getsize(rtree%h_FData,Isize)
      IF ((rtree%isizeSngl .NE. SIZE(p_FData,1)) .OR.&
          (rtree%NA .NE. SIZE(p_FData,2))) THEN
        CALL output_line('Invalid auxiliary data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayInt')
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('btree_copyToTree_arrayInt',&
            rtree%NNA,rtree%h_FData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_single2D(rtree%h_FData,rtree%p_FData)
      END IF
      ! Copy data
      rtree%p_FData(:,1:rtree%NA) = p_FData
    END IF

    IF (PRESENT(p_IData) .AND. rtree%isizeInt > 0) THEN
      CALL storage_getsize(rtree%h_IData,Isize)
      IF ((rtree%isizeInt .NE. SIZE(p_IData,1)) .OR.&
          (rtree%NA .NE. SIZE(p_IData,2))) THEN
        CALL output_line('Invalid auxiliary data!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyToTree_arrayInt')
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('btree_copyToTree_arrayInt',&
            rtree%NNA,rtree%h_IData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_int2D(rtree%h_IData,rtree%p_IData)
      END IF
      ! Copy data
      rtree%p_IData(:,1:rtree%NA) = p_IData
    END IF

  END SUBROUTINE btree_copyToTree_arrayInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyFromTreeKey_handle(rtree,h_Key)

!<description>
    ! This subroutine copies the key of the tree to the handle h_Key
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
!</input>

!<inputoutput>
    ! handle to the key
    INTEGER, INTENT(INOUT) :: h_Key
!</inputoutput>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:), POINTER   :: p_DKey
    REAL(SP), DIMENSION(:), POINTER   :: p_FKey
    INTEGER,  DIMENSION(:), POINTER   :: p_IKey
    INTEGER(I32) :: isize
    INTEGER(PREC_TREEIDX) :: j

    ! Check if handle needs to be (re-)allocated
    IF (h_Key .EQ. ST_NOHANDLE) THEN
      CALL storage_new('btree_copyFromTreeKey_handle','Key',rtree%NA,&
          rtree%ctreeFormat,h_Key,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_Key,isize)
      IF (isize < rtree%NA) THEN
        CALL storage_realloc('btree_copyFromTreeKey_handle',rtree%NA,&
            h_Key,ST_NEWBLOCK_ZERO,.FALSE.)
      END IF
    END IF

    ! What kind of key is used?
    SELECT CASE(rtree%ctreeFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(h_Key,p_DKey); j=0
      CALL inorderKeyDble(rtree%p_Kchild(TRIGHT,TROOT))
      
    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_Key,p_FKey); j=0
      CALL inorderKeySngl(rtree%p_Kchild(TRIGHT,TROOT))
      
    CASE (ST_INT)
      CALL storage_getbase_int(h_Key,p_IKey); j=0
      CALL inorderKeyInt(rtree%p_Kchild(TRIGHT,TROOT))
      
    CASE DEFAULT
      CALL output_line('Unsupported data format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
      CALL sys_halt()
    END SELECT
    
  CONTAINS

    ! Here, the real working routines follow.
    
    !**************************************************************
    ! Copy the content of the Double key to an array
    
    RECURSIVE SUBROUTINE inorderKeyDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderKeyDble(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_DKey(j)=rtree%p_DKey(i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderKeyDble(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderKeyDble
    
    !**************************************************************
    ! Copy the content of the Single key to an array

    RECURSIVE SUBROUTINE inorderKeySngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderKeySngl(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_FKey(j)=rtree%p_FKey(i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderKeySngl(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderKeySngl

    !**************************************************************
    ! Copy the content of the Integer key to an array

    RECURSIVE SUBROUTINE inorderKeyInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderKeyInt(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_IKey(j)=rtree%p_IKey(i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderKeyInt(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderKeyInt
  END SUBROUTINE btree_copyFromTreeKey_handle

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyFromTreeKey_arrayDble(rtree,p_DKey)

!<description>
    ! This subroutine copies the key of the tree to the double-valued array p_DKey
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
!</input>

!<inputoutput>
    ! double-valued array
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: p_DKey
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: j

    ! Check data format
    IF (rtree%ctreeFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Invalid data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTreeKey_arrayDble')
      CALL sys_halt()
    END IF

    ! Check size of array
    IF (SIZE(p_DKey) < rtree%NA) THEN
      CALL output_line('Array too small!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTreeKey_arrayDble')
      CALL sys_halt()
    END IF
    
    j=0
    CALL inorderKeyDble(rtree%p_Kchild(TRIGHT,TROOT))

  CONTAINS
    
    ! Here, the real working routine follows
    
    !**************************************************************
    ! Copy the content of the Double key to an array
    
    RECURSIVE SUBROUTINE inorderKeyDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderKeyDble(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_DKey(j)=rtree%p_DKey(i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderKeyDble(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderKeyDble
  END SUBROUTINE btree_copyFromTreeKey_arrayDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyFromTreeKey_arraySngl(rtree,p_FKey)

!<description>
    ! This subroutine copies the key of the tree to the single-valued array p_FKey
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
!</input>

!<inputoutput>
    ! single-valued array
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: p_FKey
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: j

    ! Check data format
    IF (rtree%ctreeFormat .NE. ST_SINGLE) THEN
      CALL output_line('Invalid data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTreeKey_arraySngl')
      CALL sys_halt()
    END IF

    ! Check size of array
    IF (SIZE(p_FKey) < rtree%NA) THEN
      CALL output_line('Array too small!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTreeKey_arraySngl')
      CALL sys_halt()
    END IF
    
    j=0
    CALL inorderKeySngl(rtree%p_Kchild(TRIGHT,TROOT))

  CONTAINS

    ! Here, the real working routine follows

    !**************************************************************
    ! Copy the content of the Single key to an array

    RECURSIVE SUBROUTINE inorderKeySngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderKeySngl(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_FKey(j)=rtree%p_FKey(i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderKeySngl(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderKeySngl
  END SUBROUTINE btree_copyFromTreeKey_arraySngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyFromTreeKey_arrayInt(rtree,p_IKey)

!<description>
    ! This subroutine copies the key of the tree to the integer-valued array p_IKey
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
!</input>

!<inputoutput>
    ! integer-valued array
    INTEGER, DIMENSION(:), INTENT(INOUT) :: p_IKey
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: j

    ! Check data format
    IF (rtree%ctreeFormat .NE. ST_INT) THEN
      CALL output_line('Invalid data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTreeKey_arrayInt')
      PRINT *, "btree_copyFromTreeKey_arrayInt: Invalid data format!"
      CALL sys_halt()
    END IF

    ! Check size of array
    IF (SIZE(p_IKey) < rtree%NA) THEN
      CALL output_line('Array too small!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTreeKey_arrayInt')
      CALL sys_halt()
    END IF
    
    j=0
    CALL inorderKeyInt(rtree%p_Kchild(TRIGHT,TROOT))

  CONTAINS
    
    ! Here, the real working routine follows

    !**************************************************************
    ! Copy the content of the Integer key to an array
    
    RECURSIVE SUBROUTINE inorderKeyInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderKeyInt(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_IKey(j)=rtree%p_IKey(i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderKeyInt(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderKeyInt
  END SUBROUTINE btree_copyFromTreeKey_arrayInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyFromTree_handle(rtree,ctype,h_Data,mask)

!<description>
    ! This subroutine copies some (part of) the auxiliary data
    ! attached to the tree to the given handle. The part of the
    ! data to be copied can be specified by the optional mask, e.g.
    ! mask=(/1,2,7/) only copies the first, second and seventh
    ! components of the auxiliary data.
    ! 
    ! The type of auxiliary data to be copied (ST_DOUBLE, ST_SINGLE,
    ! ST_INT) is specified by ctype. If the handle is already 
    ! associated, then it is reallocated if required but the type
    ! is not allowed to change.
    !
    ! Remark: There are quite many checks in this subroutine but it
    ! seems that there are some checks missing, e.g., if the mask
    ! matches the number of data present (or not) in the tree.
    ! Be informed that these checks are performed in the routines
    ! which are called for the arrays, and hece, can be omitted here.
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
    
    ! type of data
    INTEGER, INTENT(IN) :: ctype

    ! OPTIONAL: mask of components to be copied
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: mask
!</input>

!<inputoutput>
    ! handle to the data
    INTEGER, INTENT(INOUT) :: h_Data
!</inputoutput>
!</subroutine>

    ! local variables
    REAL(DP), DIMENSION(:), POINTER   :: p_DData
    REAL(SP), DIMENSION(:), POINTER   :: p_FData
    INTEGER,  DIMENSION(:), POINTER   :: p_IData
    REAL(DP), DIMENSION(:,:), POINTER :: p_DData2D
    REAL(SP), DIMENSION(:,:), POINTER :: p_FData2D
    INTEGER,  DIMENSION(:,:), POINTER :: p_IData2D
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(I32) :: iisize
    INTEGER :: idimension,idatatype
    
    ! Check if handle is already associated
    IF (h_Data .EQ. ST_NOHANDLE) THEN
      
      ! Get second working dimension
      Isize(2)=rtree%NA

      ! What kind of data should be copied?
      SELECT CASE(ctype)
      CASE (ST_DOUBLE)

        ! Get first working dimension
        Isize(1)=rtree%isizeDble
        IF (PRESENT(mask)) Isize(1) = SIZE(mask)
        
        IF (Isize(1) > 1) THEN
          CALL storage_new('btree_copyFromTree_handle','DData',Isize,&
              ST_DOUBLE,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_double2D(h_Data,p_DData2D)
          CALL btree_copyFromTree_arrayDble2D(rtree,p_DData2D,mask)
        ELSEIF (PRESENT(mask)) THEN
          CALL storage_new('btree_copyFromTree_handle','DData',Isize(2),&
              ST_DOUBLE,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_double(h_Data,p_DData)
          CALL btree_copyFromTree_arrayDble(rtree,p_DData,mask(1))
        ELSE
          CALL storage_new('btree_copyFromTree_handle','DData',Isize(2),&
              ST_DOUBLE,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_double(h_Data,p_DData)
          CALL btree_copyFromTree_arrayDble(rtree,p_DData,1)
        END IF

      CASE (ST_SINGLE)
        
        ! Get first working dimension
        Isize(1)=rtree%isizeSngl
        IF (PRESENT(mask)) Isize(1) = SIZE(mask)

        IF (Isize(1) > 1) THEN
          CALL storage_new('btree_copyFromTree_handle','FData',Isize,&
              ST_SINGLE,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_single2D(h_Data,p_FData2D)
          CALL btree_copyFromTree_arraySngl2D(rtree,p_FData2D,mask)
        ELSEIF (PRESENT(mask)) THEN
          CALL storage_new('btree_copyFromTree_handle','FData',Isize(2),&
              ST_SINGLE,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_single(h_Data,p_FData)
          CALL btree_copyFromTree_arraySngl(rtree,p_FData,mask(1))
        ELSE
          CALL storage_new('btree_copyFromTree_handle','FData',Isize(2),&
              ST_SINGLE,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_single(h_Data,p_FData)
          CALL btree_copyFromTree_arraySngl(rtree,p_FData,1)
        END IF

      CASE (ST_INT)
        
        ! Get first working dimension
        Isize(1)=rtree%isizeInt
        IF (PRESENT(mask)) Isize(1) = SIZE(mask)

        IF (Isize(1) > 1) THEN
          CALL storage_new('btree_copyFromTree_handle','IData',Isize,&
              ST_INT,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int2D(h_Data,p_IData2D)
          CALL btree_copyFromTree_arrayInt2D(rtree,p_IData2D,mask)
        ELSEIF (PRESENT(mask)) THEN
          CALL storage_new('btree_copyFromTree_handle','IData',Isize(2),&
              ST_INT,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int(h_Data,p_IData)
          CALL btree_copyFromTree_arrayInt(rtree,p_IData,mask(1))
        ELSE
          CALL storage_new('btree_copyFromTree_handle','IData',Isize(2),&
              ST_INT,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int(h_Data,p_IData)
          CALL btree_copyFromTree_arrayInt(rtree,p_IData,1)
        END IF
        
      CASE DEFAULT
        CALL output_line('Unsupported data format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
        CALL sys_halt()
      END SELECT

    ELSE   ! The handle is already associated

      ! Check if data type is valid
      CALL storage_getdatatype(h_Data,idatatype)
      IF (idatatype .NE. ctype) THEN
        CALL output_line('Data type mismatch!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
        CALL sys_halt()
      END IF

      ! Are we 1D- or 2D-array?
      CALL storage_getdimension(h_Data,idimension)
      SELECT CASE(idimension)

      CASE (1)   !!! 1D-ARRAY  !!!

        ! Do we have to reallocate the array?
        CALL storage_getsize(h_Data,iisize)
        IF (iisize < rtree%NA)&
            CALL storage_realloc('btree_copyFromTree_handle',rtree%NA,&
            h_Data,ST_NEWBLOCK_NOINIT,.FALSE.)

        ! What kind of data type are we?
        SELECT CASE(ctype)
        CASE (ST_DOUBLE)
          CALL storage_getbase_double(h_Data,p_DData)

          ! Which component should be copied
          IF (PRESENT(mask)) THEN
            IF (SIZE(mask) > 1) THEN
              CALL output_line('For 1D-array mask can only have one entry!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
              CALL sys_halt()
            END IF
            CALL btree_copyFromTree_arrayDble(rtree,p_DData,mask(1))
          ELSEIF (rtree%isizeDble .EQ. 1) THEN
            CALL btree_copyFromTree_arrayDble(rtree,p_DData,1)
          ELSE
            CALL output_line('A 1D-array was given but there are more than &
                &one components in the tree!',&
                OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
            CALL sys_halt()
          END IF
                    
        CASE (ST_SINGLE)
          CALL storage_getbase_single(h_Data,p_FData)

          ! Which component should be copied
          IF (PRESENT(mask)) THEN
            IF (SIZE(mask) > 1) THEN
              CALL output_line('For 1D-array mask can only have one entry!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
              CALL sys_halt()
            END IF
            CALL btree_copyFromTree_arraySngl(rtree,p_FData,mask(1))
          ELSEIF (rtree%isizeSngl .EQ. 1) THEN
            CALL btree_copyFromTree_arraySngl(rtree,p_FData,1)
          ELSE
            CALL output_line('A 1D-array was given but there are more than &
                &one components in the tree!',&
                OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
            CALL sys_halt()
          END IF
          
        CASE (ST_INT)
          CALL storage_getbase_int(h_Data,p_IData)

          ! Which component should be copied
          IF (PRESENT(mask)) THEN
            IF (SIZE(mask) > 1) THEN
              CALL output_line('For 1D-array mask can only have one entry!',&
                  OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
              CALL sys_halt()
            END IF
            CALL btree_copyFromTree_arrayInt(rtree,p_IData,mask(1))
          ELSEIF (rtree%isizeInt .EQ. 1) THEN
            CALL btree_copyFromTree_arrayInt(rtree,p_IData,1)
          ELSE
            CALL output_line('A 1D-array was given but there are more than &
                &one components in the tree!',&
                OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
            CALL sys_halt()
          END IF

        CASE DEFAULT
          CALL output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
          CALL sys_halt()
        END SELECT

      CASE (2)   !!! 2D-ARRAY !!!

        ! Do we have to reallocate the array?
        CALL storage_getsize(h_Data,Isize)
        IF (Isize(2) < rtree%NA)&
            CALL storage_realloc('btree_copyFromTree_handle',rtree%NA,&
            h_Data,ST_NEWBLOCK_NOINIT,.FALSE.)

        ! What kind of data type are we
        SELECT CASE(ctype)
        CASE (ST_DOUBLE)
          CALL storage_getbase_double2D(h_Data,p_DData2D)
          CALL btree_copyFromTree_arrayDble2D(rtree,p_DData2D,mask)

        CASE (ST_SINGLE)
          CALL storage_getbase_single2D(h_Data,p_FData2D)
          CALL btree_copyFromTree_arraySngl2D(rtree,p_FData2D,mask)

        CASE (ST_INT)
          CALL storage_getbase_int2D(h_Data,p_IData2D)
          CALL btree_copyFromTree_arrayInt2D(rtree,p_IData2D,mask)

        CASE DEFAULT
          CALL output_line('Unsupported data type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
          CALL sys_halt()
        END SELECT

      CASE DEFAULT
        CALL output_line('Unsupported data dimension!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_handle')
        CALL sys_halt()
      END SELECT
    END IF
  END SUBROUTINE btree_copyFromTree_handle

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyFromTree_arrayDble(rtree,p_Ddata,mask)

!<description>
    ! This subroutine copies one component of the auxiliary double data
    ! attached to the tree to the given array. The component must be specified
    ! by the mask, e.g. mask=4 copies the fourth component of the data.
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
    
    ! mask of component to be copied
    INTEGER, INTENT(IN) :: mask
!</input>

!<inputoutput>
    ! double data array
    REAL(DP), DIMENSION(:),INTENT(INOUT) :: p_DData
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: j

    ! Check if auxiliary double data is available
    IF (rtree%isizeDble .EQ. 0) THEN
      CALL output_line('No double data available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble')
      CALL sys_halt()
    END IF

    ! Check if given array is large enough in its second dimension
    IF (SIZE(p_DData) < rtree%NA) THEN
      CALL output_line('Array too small!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble')
      CALL sys_halt()
    END IF

    ! Check if mask is valid
    IF (mask < 1 .OR. mask > rtree%isizeDble) THEN
      CALL output_line('Invalid mask!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble')
      CALL sys_halt()
    END IF
    
    ! Do the work
    j=0; CALL inorderDataDble(rtree%p_Kchild(TRIGHT,TROOT))

  CONTAINS

    ! Here, the real working routine follows.
    
    !**************************************************************
    ! Copy the content of the Double data to an array
    
    RECURSIVE SUBROUTINE inorderDataDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderDataDble(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_DData(j)=rtree%p_DData(mask,i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderDataDble(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataDble

  END SUBROUTINE btree_copyFromTree_arrayDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyFromTree_arrayDble2D(rtree,p_Ddata,mask)

!<description>
    ! This subroutine copies some (part of) the auxiliary double data
    ! attached to the tree to the given array. The part of the data to
    ! be copied can be specified by the optional mask, e.g.
    ! mask=(/1,2,7/) only copies the first, second and seventh
    ! components of the auxiliary data.
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
    
    ! OPTIONAL: mask of components to be copied
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: mask
!</input>

!<inputoutput>
    ! double data array
    REAL(DP), DIMENSION(:,:),INTENT(INOUT) :: p_DData
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_TREEIDX) :: j

    ! Check if auxiliary double data is available
    IF (rtree%isizeDble .EQ. 0) THEN
      CALL output_line('No double data available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble2D')
      CALL sys_halt()
    END IF

    ! Check if given array is large enough in its second dimension
    Isize = SHAPE(p_Ddata)
    IF (Isize(2) < rtree%NA) THEN
      CALL output_line('Array too small!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble2D')
      CALL sys_halt()
    END IF

    ! Copy the content to the array
    IF (PRESENT(mask)) THEN
      ! Check if given array has correct size in its first dimension
      IF (Isize(1) .NE. SIZE(mask)) THEN
        CALL output_line('Array dimensions do not match mask!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble2D')
        CALL sys_halt()
      END IF

      ! Check if mask is valid
      IF (ANY(mask < 1) .OR. ANY(mask > rtree%isizeDble)) THEN
        CALL output_line('Invalid mask!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble2D')
        CALL sys_halt()
      END IF
      j=0; CALL inorderDataDbleMask(rtree%p_Kchild(TRIGHT,TROOT))
    ELSE
      ! Check if given array has correct size in its first dimension
      IF (Isize(1) .NE. rtree%isizeDble) THEN
        CALL output_line('Array too small!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayDble2D')
        CALL sys_halt()
      END IF
      j=0; CALL inorderDataDble(rtree%p_Kchild(TRIGHT,TROOT))
    END IF

  CONTAINS

    ! Here, the real working routines follow.

    !**************************************************************
    ! Copy the content of the Double data to an array

    RECURSIVE SUBROUTINE inorderDataDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderDataDble(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_DData(:,j)=rtree%p_DData(:,i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderDataDble(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataDble

    !**************************************************************
    ! Copy the content of the Double data to an array (masked)

    RECURSIVE SUBROUTINE inorderDataDbleMask(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderDataDbleMask(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_DData(:,j)=rtree%p_DData(mask,i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderDataDbleMask(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataDbleMask

  END SUBROUTINE btree_copyFromTree_arrayDble2D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyFromTree_arraySngl(rtree,p_FData,mask)

!<description>
    ! This subroutine copies one component of the auxiliary single data
    ! attached to the tree to the given array. The component must be specified
    ! by the mask, e.g. mask=4 copies the fourth component of the data.
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
    
    ! mask of component to be copied
    INTEGER, INTENT(IN) :: mask
!</input>

!<inputoutput>
    ! single data array
    REAL(SP), DIMENSION(:),INTENT(INOUT) :: p_FData
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: j

    ! Check if auxiliary single data is available
    IF (rtree%isizeSngl .EQ. 0) THEN
      CALL output_line('No single data available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl')
      CALL sys_halt()
    END IF

    ! Check if given array is large enough in its second dimension
    IF (SIZE(p_FData) < rtree%NA) THEN
      CALL output_line('Array too small!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl')
      CALL sys_halt()
    END IF

    ! Check if mask is valid
    IF (mask < 1 .OR. mask > rtree%isizeSngl) THEN
      CALL output_line('Invalid mask!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl')
      CALL sys_halt()
    END IF
    
    ! Do the work
    j=0; CALL inorderDataSngl(rtree%p_Kchild(TRIGHT,TROOT))

  CONTAINS

    ! Here, the real working routine follows.
    
    !**************************************************************
    ! Copy the content of the Single data to an array
    
    RECURSIVE SUBROUTINE inorderDataSngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderDataSngl(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_FData(j)=rtree%p_FData(mask,i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderDataSngl(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataSngl

  END SUBROUTINE btree_copyFromTree_arraySngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyFromTree_arraySngl2D(rtree,p_Fdata,mask)

!<description>
    ! This subroutine copies some (part of) the auxiliary single data
    ! attached to the tree to the given array. The part of the data to
    ! be copied can be specified by the optional mask, e.g.
    ! mask=(/1,2,7/) only copies the first, second and seventh
    ! components of the auxiliary data.
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
    
    ! OPTIONAL: mask of components to be copied
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: mask
!</input>

!<inputoutput>
    ! single data array
    REAL(SP), DIMENSION(:,:),INTENT(INOUT) :: p_FData
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_TREEIDX) :: j

    ! Check if auxiliary single data is available
    IF (rtree%isizeSngl .EQ. 0) THEN
      CALL output_line('No single data available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl2D')
      CALL sys_halt()
    END IF
    
    ! Check if given array is large enough in its second dimension
    Isize = SHAPE(p_Fdata)
    IF (Isize(2) < rtree%NA) THEN
      CALL output_line('Array too small!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl2D')
      CALL sys_halt()
    END IF

    ! Copy the content to the array
    IF (PRESENT(mask)) THEN
      ! Check if given array has correct size in its first dimension
      IF (Isize(1) .NE. SIZE(mask)) THEN
        CALL output_line('Array dimensions do not match mask!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl2D')
        CALL sys_halt()
      END IF

      ! Check if mask is valid
      IF (ANY(mask < 1) .OR. ANY(mask > rtree%isizeSngl)) THEN
        CALL output_line('Invalid mask!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl2D')
        CALL sys_halt()
      END IF
      j=0; CALL inorderDataSnglMask(rtree%p_Kchild(TRIGHT,TROOT))
    ELSE
      ! Check if given array has correct size in its first dimension
      IF (Isize(1) .NE. rtree%isizeSngl) THEN
        CALL output_line('Array too small!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arraySngl2D')
        CALL sys_halt()
      END IF
      j=0; CALL inorderDataSngl(rtree%p_Kchild(TRIGHT,TROOT))
    END IF

  CONTAINS

    ! Here, the real working routines follow.

    !**************************************************************
    ! Copy the content of the Single data to an array

    RECURSIVE SUBROUTINE inorderDataSngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderDataSngl(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_FData(:,j)=rtree%p_FData(:,i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderDataSngl(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataSngl

    !**************************************************************
    ! Copy the content of the Single data to an array (masked)

    RECURSIVE SUBROUTINE inorderDataSnglMask(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderDataSnglMask(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_FData(:,j)=rtree%p_FData(mask,i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderDataSnglMask(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataSnglMask

  END SUBROUTINE btree_copyFromTree_arraySngl2D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyFromTree_arrayInt(rtree,p_Idata,mask)

!<description>
    ! This subroutine copies one component of the auxiliary integer data
    ! attached to the tree to the given array. The component must be specified
    ! by the mask, e.g. mask=4 copies the fourth component of the data.
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
    
    ! mask of component to be copied
    INTEGER, INTENT(IN) :: mask
!</input>

!<inputoutput>
    ! integer data array
    INTEGER, DIMENSION(:),INTENT(INOUT) :: p_IData
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: j

    ! Check if auxiliary integer data is available
    IF (rtree%isizeInt .EQ. 0) THEN
      CALL output_line('No integer data available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt')
      CALL sys_halt()
    END IF

    ! Check if given array is large enough in its second dimension
    IF (SIZE(p_IData) < rtree%NA) THEN
      CALL output_line('Array too small!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt')
      CALL sys_halt()
    END IF

    ! Check if mask is valid
    IF (mask < 1 .OR. mask > rtree%isizeInt) THEN
      CALL output_line('Invalid mask!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt')
      CALL sys_halt()
    END IF
    
    ! Do the work
    j=0; CALL inorderDataInt(rtree%p_Kchild(TRIGHT,TROOT))

  CONTAINS

    ! Here, the real working routine follows.
    
    !**************************************************************
    ! Copy the content of the Integer data to an array
    
    RECURSIVE SUBROUTINE inorderDataInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderDataInt(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_IData(j)=rtree%p_IData(mask,i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderDataInt(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataInt

  END SUBROUTINE btree_copyFromTree_arrayInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_copyFromTree_arrayInt2D(rtree,p_Idata,mask)

!<description>
    ! This subroutine copies some (part of) the auxiliary integer data
    ! attached to the tree to the given array. The part of the data to
    ! be copied can be specified by the optional mask, e.g.
    ! mask=(/1,2,7/) only copies the first, second and seventh
    ! components of the auxiliary data.
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
    
    ! OPTIONAL: mask of components to be copied
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: mask
!</input>

!<inputoutput>
    ! integer data array
    INTEGER, DIMENSION(:,:),INTENT(INOUT) :: p_IData
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_TREEIDX) :: j

    ! Check if auxiliary integer data is available
    IF (rtree%isizeInt .EQ. 0) THEN
      CALL output_line('No integer data available!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt2D')
      CALL sys_halt()
    END IF

    ! Check if given array is large enough in its second dimension
    Isize = SHAPE(p_IData)
    IF (Isize(2) < rtree%NA) THEN
      CALL output_line('Invalid mask!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt2D')
      CALL sys_halt()
    END IF

    ! Copy the content to the array
    IF (PRESENT(mask)) THEN
      ! Check if given array has correct size in its first dimension
      IF (Isize(1) .NE. SIZE(mask)) THEN
        CALL output_line('Array dimensions do not match mask!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt2D')
        CALL sys_halt()
      END IF

      ! Check if mask is valid
      IF (ANY(mask < 1) .OR. ANY(mask > rtree%isizeInt)) THEN
        CALL output_line('Invalid mask!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt2D')
        CALL sys_halt()
      END IF
      j=0; CALL inorderDataIntMask(rtree%p_Kchild(TRIGHT,TROOT))
    ELSE
      ! Check if given array has correct size in its first dimension
      IF (Isize(1) .NE. rtree%isizeInt) THEN
        CALL output_line('Array too small!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_copyFromTree_arrayInt2D')
        CALL sys_halt()
      END IF
      j=0; CALL inorderDataInt(rtree%p_Kchild(TRIGHT,TROOT))
    END IF

  CONTAINS

    ! Here, the real working routines follow.

    !**************************************************************
    ! Copy the content of the Integer data to an array

    RECURSIVE SUBROUTINE inorderDataInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderDataInt(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_IData(:,j)=rtree%p_IData(:,i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderDataInt(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataInt

    !**************************************************************
    ! Copy the content of the Integer data to an array (masked)

    RECURSIVE SUBROUTINE inorderDataIntMask(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderDataIntMask(rtree%p_Kchild(TLEFT,i))
      j=j+1; p_IData(:,j)=rtree%p_IData(mask,i)
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderDataIntMask(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataIntMask

  END SUBROUTINE btree_copyFromTree_arrayInt2D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_insertToTreeDble(rtree,dkey,DData,FData,IData,iposOpt)

!<description>
    ! This subroutine inserts a new Double key into the tree and
    ! (possibly) some auxiliary data
!</description>

!<input>
    ! key
    REAL(DP), INTENT(IN) :: dkey
    
    ! OPTIONAL: Double data
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DData

    ! OPTIONAL: Single data
    REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: FData

    ! OPTIONAL: Integer data
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: IData
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>

!<output>
    ! OPTIONAL: Position of the new entry. If ipos < 0 then
    ! the entry already exists and ABS(ipos) is its position
    INTEGER(PREC_TREEIDX), INTENT(OUT), OPTIONAL :: iposOpt
!</output>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_TREEIDX) :: ipos,jpos
    
    ! Check if tree format is ok
    IF (rtree%ctreeFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_insertToTreeDble')
      CALL sys_halt()
    END IF

    ! Check if key is already stored in tree
    IF (btree_searchInTree(rtree,dkey,jpos) .EQ. BTREE_NOT_FOUND) THEN
      
      ! Adjust size and position
      rtree%na = rtree%na+1
      ipos     = rtree%p_Kchild(TFREE,TROOT)

      IF (PRESENT(iposOpt)) iposOpt=ABS(ipos)

      ! Check if memory needs to be expanded
      IF (ABS(ipos) > rtree%nna) &
          CALL btree_resizeTree(rtree,CEILING(rtree%dfactor*rtree%nna))
      
      ! Compute next free position in memory
      IF (ipos > 0) THEN
        rtree%p_Kchild(TFREE,TROOT) = ipos+1
      ELSE
        ipos                      = ABS(ipos)
        rtree%p_Kchild(TFREE,TROOT) = rtree%p_Kchild(TFREE,ipos)
      END IF
      
      ! Store given data in memory
      rtree%p_DKey(ipos)     = dkey
      rtree%p_Kbal(ipos)     = 0
      rtree%p_Kchild(:,ipos) = TNULL
      
      ! Store optional data in memory
      IF ((rtree%isizeInt > 0) .AND. &
          PRESENT(IData)) rtree%p_IData(:,ipos) = IData

      IF ((rtree%isizeDble > 0) .AND. &
          PRESENT(DData)) rtree%p_DData(:,ipos) = DData

      IF ((rtree%isizeSngl > 0) .AND. &
          PRESENT(FData)) rtree%p_FData(:,ipos) = FData
      
      ! Insert new node into the tree and into the search path
      rtree%p_Kchild(MERGE(TLEFT,TRIGHT,jpos < 0),ABS(jpos)) = ipos
      rtree%p_Kpath(rtree%depth+1) = MERGE(-ipos,ipos,jpos <= 0)
      
      ! Invoke rebalance procedure
      IF (rtree%depth > 0) CALL rebalanceAfterInsert(rtree,rtree%depth)

    ELSEIF(PRESENT(iposOpt)) THEN
      iposOpt=-rtree%p_Kchild(MERGE(TLEFT,TRIGHT,jpos < 0),ABS(jpos))
    END IF
  END SUBROUTINE btree_insertToTreeDble
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_insertToTreeSngl(rtree,skey,DData,FData,IData,iposOpt)

!<description>
    ! This subroutine inserts a new Single key into the tree and
    ! (possibly) some auxiliary data
!</description>

!<input>
    ! key
    REAL(SP), INTENT(IN) :: skey
    
    ! OPTIONAL: Double data
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DData

    ! OPTIONAL: Single data
    REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: FData

    ! OPTIONAL: Integer data
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: IData
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>

!<output>
    ! OPTIONAL: Position of the new entry. If ipos < 0 then
    ! the entry already exists and ABS(ipos) is its position
    INTEGER(PREC_TREEIDX), INTENT(OUT), OPTIONAL :: iposOpt
!</output>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_TREEIDX) :: ipos,jpos
    
    ! Check if tree format is ok
    IF (rtree%ctreeFormat .NE. ST_SINGLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_insertToTreeSngl')
      CALL sys_halt()
    END IF

    ! Check if key is already stored in tree
    IF (btree_searchInTree(rtree,skey,jpos) .EQ. BTREE_NOT_FOUND) THEN
      
      ! Adjust size and position
      rtree%na = rtree%na+1
      ipos     = rtree%p_Kchild(TFREE,TROOT)

      IF (PRESENT(iposOpt)) iposOpt=ABS(ipos)

      ! Check if memory needs to be expanded
      IF (ABS(ipos) > rtree%nna) &
          CALL btree_resizeTree(rtree,CEILING(rtree%dfactor*rtree%nna))
      
      ! Compute next free position in memory
      IF (ipos > 0) THEN
        rtree%p_Kchild(TFREE,TROOT) = ipos+1
      ELSE
        ipos                      = ABS(ipos)
        rtree%p_Kchild(TFREE,TROOT) = rtree%p_Kchild(TFREE,ipos)
      END IF
      
      ! Store given data in memory
      rtree%p_FKey(ipos)     = skey
      rtree%p_Kbal(ipos)     = 0
      rtree%p_Kchild(:,ipos) = TNULL
      
      ! Store optional data in memory
      IF ((rtree%isizeInt > 0) .AND. &
          PRESENT(IData)) rtree%p_IData(:,ipos) = IData

      IF ((rtree%isizeDble > 0) .AND. &
          PRESENT(DData)) rtree%p_DData(:,ipos) = DData

      IF ((rtree%isizeSngl > 0) .AND. &
          PRESENT(FData)) rtree%p_FData(:,ipos) = FData
      
      ! Insert new node into the tree and into the search path
      rtree%p_Kchild(MERGE(TLEFT,TRIGHT,jpos < 0),ABS(jpos)) = ipos
      rtree%p_Kpath(rtree%depth+1) = MERGE(-ipos,ipos,jpos <= 0)
      
      ! Invoke rebalance procedure
      IF (rtree%depth > 0) CALL rebalanceAfterInsert(rtree,rtree%depth)

    ELSEIF(PRESENT(iposOpt)) THEN
      iposOpt=-rtree%p_Kchild(MERGE(TLEFT,TRIGHT,jpos < 0),ABS(jpos))
    END IF
  END SUBROUTINE btree_insertToTreeSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_insertToTreeInt(rtree,ikey,DData,FData,IData,iposOpt)

!<description>
    ! This subroutine inserts a new Integer key into the tree and
    ! (possibly) some auxiliary data
!</description>

!<input>
    ! key
    INTEGER(PREC_TREEIDX), INTENT(IN) :: ikey
    
    ! OPTIONAL: Double data
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: DData

    ! OPTIONAL: Single data
    REAL(SP), DIMENSION(:), INTENT(IN), OPTIONAL :: FData

    ! OPTIONAL: Integer data
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: IData
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>

!<output>
    ! OPTIONAL: Position of the new entry. If ipos < 0 then
    ! the entry already exists and ABS(ipos) is its position
    INTEGER(PREC_TREEIDX), INTENT(OUT), OPTIONAL :: iposOpt
!</output>
!</subroutine>
    
    ! local variables
    INTEGER(PREC_TREEIDX) :: ipos,jpos
    
    ! Check if tree format is ok
    IF (rtree%ctreeFormat .NE. ST_INT) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_insertToTreeInt')
      CALL sys_halt()
    END IF

    ! Check if key is already stored in tree
    IF (btree_searchInTree(rtree,ikey,jpos) .EQ. BTREE_NOT_FOUND) THEN
      
      ! Adjust size and position
      rtree%na = rtree%na+1
      ipos     = rtree%p_Kchild(TFREE,TROOT)

      IF (PRESENT(iposOpt)) iposOpt=ABS(ipos)
      
      ! Check if memory needs to be expanded
      IF (ABS(ipos) > rtree%nna) &
          CALL btree_resizeTree(rtree,CEILING(rtree%dfactor*rtree%nna))
      
      ! Compute next free position in memory
      IF (ipos > 0) THEN
        rtree%p_Kchild(TFREE,TROOT) = ipos+1
      ELSE
        ipos                      = ABS(ipos)
        rtree%p_Kchild(TFREE,TROOT) = rtree%p_Kchild(TFREE,ipos)
      END IF
      
      ! Store given data in memory
      rtree%p_IKey(ipos)     = ikey
      rtree%p_Kbal(ipos)     = 0
      rtree%p_Kchild(:,ipos) = TNULL
      
      ! Store optional data in memory
      IF ((rtree%isizeInt > 0) .AND. &
          PRESENT(IData)) rtree%p_IData(:,ipos) = IData

      IF ((rtree%isizeDble > 0) .AND. &
          PRESENT(DData)) rtree%p_DData(:,ipos) = DData

      IF ((rtree%isizeSngl > 0) .AND. &
          PRESENT(FData)) rtree%p_FData(:,ipos) = FData
      
      ! Insert new node into the tree and into the search path
      rtree%p_Kchild(MERGE(TLEFT,TRIGHT,jpos < 0),ABS(jpos)) = ipos
      rtree%p_Kpath(rtree%depth+1) = MERGE(-ipos,ipos,jpos <= 0)
      
      ! Invoke rebalance procedure
      IF (rtree%depth > 0) CALL rebalanceAfterInsert(rtree,rtree%depth)
      
    ELSEIF(PRESENT(iposOpt)) THEN
      iposOpt=-rtree%p_Kchild(MERGE(TLEFT,TRIGHT,jpos < 0),ABS(jpos))
    END IF
  END SUBROUTINE btree_insertToTreeInt

  ! ***************************************************************************

!<subroutine>
  
  RECURSIVE SUBROUTINE rebalanceAfterInsert(rtree,i)

!<description>
    ! This subroutine rebalances the AVL tree after insertion
!</description>

!<input>
    ! starting node
    INTEGER(PREC_TREEIDX), INTENT(IN) :: i
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>
!</subroutine>    
    
    ! local variables
    INTEGER(PREC_TREEIDX) :: v,x,w
    INTEGER :: dir
    
    v   = rtree%p_Kpath(i)
    dir = MERGE(TLEFT,TRIGHT,v < 0)
    v   = ABS(v)
    
    ! Which rotation needs to be performed?
    SELECT CASE (dir)
      
    CASE (TLEFT)
      ! Node V has been reached from the left
      SELECT CASE(rtree%p_Kbal(v))
        
      CASE (1) ! bal(v)=1
        rtree%p_Kbal(v) = 0
        
      CASE (0) ! bal(v)=0
        rtree%p_Kbal(v) = -1
        IF (v .NE. rtree%p_Kchild(TRIGHT,TROOT)) &
            CALL rebalanceAfterInsert(rtree,i-1)
        
      CASE (-1) ! bal(v)=-1
        x = rtree%p_Kchild(TLEFT,v)
        w = rtree%p_Kchild(TRIGHT,x)
        
        SELECT CASE(rtree%p_Kbal(x))
          
        CASE (-1,0) ! bal(x)=-1 or bal(x)=0
          ! Single-Right-rotation
          rtree%p_Kchild(TLEFT,v)  = rtree%p_Kchild(TRIGHT,x)
          rtree%p_Kchild(TRIGHT,x) = v
          rtree%p_Kbal(x)          = rtree%p_Kbal(x)+1
          rtree%p_Kbal(v)          = -rtree%p_Kbal(x)
          
          IF (v .EQ. rtree%p_Kchild(TRIGHT,TROOT)) THEN
            rtree%p_Kchild(TRIGHT,TROOT) = x
          ELSE
            rtree%p_Kchild(MERGE(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                ABS(rtree%p_Kpath(i-1))) = x
          END IF
          
        CASE (1) ! bal(x)=1
          ! Double-Left-Right-rotation
          rtree%p_Kchild(TLEFT,v)  = rtree%p_Kchild(TRIGHT,w)
          rtree%p_Kchild(TRIGHT,x) = rtree%p_Kchild(TLEFT,w)
          rtree%p_Kchild(TLEFT,w)  = x
          rtree%p_Kchild(TRIGHT,w) = v
          rtree%p_Kbal(v)          = -MIN(0,rtree%p_Kbal(w))
          rtree%p_Kbal(x)          = -MAX(0,rtree%p_Kbal(w))
          rtree%p_Kbal(w)          = 0
          
          IF (v .EQ. rtree%p_Kchild(TRIGHT,TROOT)) THEN
            rtree%p_Kchild(TRIGHT,TROOT) = w
          ELSE
            rtree%p_Kchild(MERGE(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                ABS(rtree%p_Kpath(i-1))) = w
          END IF
          
        END SELECT
        
      END SELECT
      
    CASE (TRIGHT)
      ! Node V has been reached from the right
      SELECT CASE(rtree%p_Kbal(v))
        
      CASE (-1) ! bal(v)=-1
        rtree%p_Kbal(v) = 0
        
      CASE (0) ! bal(v)=0
        rtree%p_Kbal(v) = 1
        IF (v .NE. rtree%p_Kchild(TRIGHT,TROOT))&
            CALL rebalanceAfterInsert(rtree,i-1)
        
      CASE (1) ! bal(v)=1
        x = rtree%p_Kchild(TRIGHT,v)
        w = rtree%p_Kchild(TLEFT,x)
        
        SELECT CASE(rtree%p_Kbal(x))
          
        CASE (0,1) ! bal(x)=0 or bal(x)=1
          ! Single-Left-rotation
          rtree%p_Kchild(TRIGHT,v) = rtree%p_Kchild(TLEFT,x)
          rtree%p_Kchild(TLEFT,x)  = v
          rtree%p_Kbal(x)          = rtree%p_Kbal(x)-1
          rtree%p_Kbal(v)          = -rtree%p_Kbal(x)
          
          IF (v .EQ. rtree%p_Kchild(TRIGHT,TROOT)) THEN
            rtree%p_Kchild(TRIGHT,TROOT) = x
          ELSE
            rtree%p_Kchild(MERGE(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                ABS(rtree%p_Kpath(i-1))) = x
          END IF
          
        CASE (-1) ! bal(x)=-1
          ! Double-Right-Left-rotation
          rtree%p_Kchild(TRIGHT,v) = rtree%p_Kchild(TLEFT,w)
          rtree%p_Kchild(TLEFT,x)  = rtree%p_Kchild(TRIGHT,w)
          rtree%p_Kchild(TLEFT,w)  = v
          rtree%p_Kchild(TRIGHT,w) = x
          rtree%p_Kbal(v)          = -MAX(0,rtree%p_Kbal(w))
          rtree%p_Kbal(x)          = -MIN(0,rtree%p_Kbal(w))
          rtree%p_Kbal(w)          = 0
          
          IF (v .EQ. rtree%p_Kchild(TRIGHT,TROOT)) THEN
            rtree%p_Kchild(TRIGHT,TROOT) = w
          ELSE
            rtree%p_Kchild(MERGE(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                ABS(rtree%p_Kpath(i-1))) = w
          END IF
          
        END SELECT
        
      END SELECT
      
    END SELECT
  END SUBROUTINE rebalanceAfterInsert

  ! ***************************************************************************

!<function>

  FUNCTION btree_deleteFromTreeDble(rtree,dkey) RESULT(f)

!<description>
    ! This functions deletes a Double key from the tree
!</description>

!<input>
    ! Key
    REAL(DP), INTENT(IN) :: dkey
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>

!<result>
    ! Result of the deletion BTREE_NOT_FOUND / BTREE_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_TREEIDX) :: ipred,ipos,jpos

    ! Check if tree format is ok
    IF (rtree%ctreeFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_deleteFromTreeDble')
      CALL sys_halt()
    END IF

    ! Search for key
    f=btree_searchInTree(rtree,dkey,ipred)
    IF (f .EQ. BTREE_NOT_FOUND) RETURN

    ! Compute new dimensions
    rtree%na = rtree%na-1
    ipos     = rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))

    ! Check if node to be deleted has two children.
    IF ((rtree%p_Kchild(TLEFT,ipos) .NE. TNULL) .AND. &
        (rtree%p_Kchild(TRIGHT,ipos) .NE. TNULL)) THEN
      
      ! Descent in the left subtree of node IPOS always to the right
      ! until a node JPOS  without right child is found and
      ! interchange data.
      jpos  = rtree%p_Kchild(TLEFT,ipos)
      ipred = -ipos
      DO
        rtree%depth              = rtree%depth+1
        rtree%p_Kpath(rtree%depth) = ipred
        
        IF (rtree%p_Kchild(TRIGHT,jpos) .EQ. TNULL) THEN
          ! Change key
          rtree%p_DKey(ipos) = rtree%p_DKey(jpos)

          ! Change auxiliary data
          IF (rtree%isizeDble > 0) rtree%p_DData(:,ipos) = rtree%p_DData(:,jpos)
          IF (rtree%isizeSngl > 0) rtree%p_FData(:,ipos) = rtree%p_FData(:,jpos)
          IF (rtree%isizeInt  > 0) rtree%p_IData(:,ipos) = rtree%p_IData(:,jpos)
          EXIT
        END IF
        ipred = jpos
        jpos = rtree%p_Kchild(TRIGHT,jpos)
      END DO
    END IF

    ! Node to be deleted has less than two childen
    ipos = rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
    
    IF (rtree%p_Kchild(TLEFT,ipos) .EQ. TNULL) THEN
      IF (rtree%p_Kchild(TRIGHT,ipos) .EQ. TNULL) THEN
        ! Node to be deleted is leaf: Nullify pointer to node and
        ! mark position in array as deleted
        rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) = TNULL
        rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
        rtree%p_Kchild(TFREE,TROOT) = -ipos
      ELSE
        ! Node to be deleted has right child: Set pointer to right
        ! child and mark position in array as deleted
        rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) =&
            rtree%p_Kchild(TRIGHT,ipos)
        rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
        rtree%p_Kchild(TFREE,TROOT) = -ipos
      END IF
    ELSE
      ! Node to be deleted has left child: Set pointer to left child
      ! and mark position in array as deleted
      rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) =&
          rtree%p_Kchild(TLEFT,ipos)
      rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
      rtree%p_Kchild(TFREE,TROOT) = -ipos
    END IF
    
    ! Invoke rebalance procedure
    IF (rtree%depth > 0) CALL rebalanceAfterDeletion(rtree,rtree%depth)
  END FUNCTION btree_deleteFromTreeDble

  ! ***************************************************************************

!<function>

  FUNCTION btree_deleteFromTreeSngl(rtree,skey) RESULT(f)

!<description>
    ! This functions deletes a Single key from the tree
!</description>

!<input>
    ! Key
    REAL(SP), INTENT(IN) :: skey
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>

!<result>
    ! Result of the deletion BTREE_NOT_FOUND / BTREE_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_TREEIDX) :: ipred,ipos,jpos

    ! Check if tree format is ok
    IF (rtree%ctreeFormat .NE. ST_SINGLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_deleteFromTreeSngl')
      CALL sys_halt()
    END IF

    ! Search for key
    f=btree_searchInTree(rtree,skey,ipred)
    IF (f .EQ. BTREE_NOT_FOUND) RETURN

    ! Compute new dimensions
    rtree%na = rtree%na-1
    ipos     = rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))

    ! Check if node to be deleted has two children.
    IF ((rtree%p_Kchild(TLEFT,ipos) .NE. TNULL) .AND. &
        (rtree%p_Kchild(TRIGHT,ipos) .NE. TNULL)) THEN
      
      ! Descent in the left subtree of node IPOS always to the right
      ! until a node JPOS  without right child is found and
      ! interchange data.
      jpos  = rtree%p_Kchild(TLEFT,ipos)
      ipred = -ipos
      DO
        rtree%depth              = rtree%depth+1
        rtree%p_Kpath(rtree%depth) = ipred
        
        IF (rtree%p_Kchild(TRIGHT,jpos) .EQ. TNULL) THEN
          ! Change key
          rtree%p_FKey(ipos) = rtree%p_FKey(jpos)

          ! Change auxiliary data
          IF (rtree%isizeDble > 0) rtree%p_DData(:,ipos) = rtree%p_DData(:,jpos)
          IF (rtree%isizeSngl > 0) rtree%p_FData(:,ipos) = rtree%p_FData(:,jpos)
          IF (rtree%isizeInt  > 0) rtree%p_IData(:,ipos) = rtree%p_IData(:,jpos)
          EXIT
        END IF
        ipred = jpos
        jpos = rtree%p_Kchild(TRIGHT,jpos)
      END DO
    END IF

    ! Node to be deleted has less than two childen
    ipos = rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
    
    IF (rtree%p_Kchild(TLEFT,ipos) .EQ. TNULL) THEN
      IF (rtree%p_Kchild(TRIGHT,ipos) .EQ. TNULL) THEN
        ! Node to be deleted is leaf: Nullify pointer to node and
        ! mark position in array as deleted
        rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) = TNULL
        rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
        rtree%p_Kchild(TFREE,TROOT) = -ipos
      ELSE
        ! Node to be deleted has right child: Set pointer to right
        ! child and mark position in array as deleted
        rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) =&
            rtree%p_Kchild(TRIGHT,ipos)
        rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
        rtree%p_Kchild(TFREE,TROOT) = -ipos
      END IF
    ELSE
      ! Node to be deleted has left child: Set pointer to left child
      ! and mark position in array as deleted
      rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) =&
          rtree%p_Kchild(TLEFT,ipos)
      rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
      rtree%p_Kchild(TFREE,TROOT) = -ipos
    END IF
    
    ! Invoke rebalance procedure
    IF (rtree%depth > 0) CALL rebalanceAfterDeletion(rtree,rtree%depth)
  END FUNCTION btree_deleteFromTreeSngl

  ! ***************************************************************************

!<function>

  FUNCTION btree_deleteFromTreeInt(rtree,ikey) RESULT(f)

!<description>
    ! This functions deletes a Integer key from the tree
!</description>

!<input>
    ! Key
    INTEGER(PREC_TREEIDX), INTENT(IN) :: ikey
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>

!<result>
    ! Result of the deletion BTREE_NOT_FOUND / BTREE_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_TREEIDX) :: ipred,ipos,jpos

    ! Check if tree format is ok
    IF (rtree%ctreeFormat .NE. ST_INT) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_deleteFromTreeInt')
      CALL sys_halt()
    END IF

    ! Search for key
    f=btree_searchInTree(rtree,ikey,ipred)
    IF (f .EQ. BTREE_NOT_FOUND) RETURN

    ! Compute new dimensions
    rtree%na = rtree%na-1
    ipos     = rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))

    ! Check if node to be deleted has two children.
    IF ((rtree%p_Kchild(TLEFT,ipos) .NE. TNULL) .AND. &
        (rtree%p_Kchild(TRIGHT,ipos) .NE. TNULL)) THEN
      
      ! Descent in the left subtree of node IPOS always to the right
      ! until a node JPOS  without right child is found and
      ! interchange data.
      jpos  = rtree%p_Kchild(TLEFT,ipos)
      ipred = -ipos
      DO
        rtree%depth              = rtree%depth+1
        rtree%p_Kpath(rtree%depth) = ipred
        
        IF (rtree%p_Kchild(TRIGHT,jpos) .EQ. TNULL) THEN
          ! Change key
          rtree%p_IKey(ipos) = rtree%p_IKey(jpos)

          ! Change auxiliary data
          IF (rtree%isizeDble > 0) rtree%p_DData(:,ipos) = rtree%p_DData(:,jpos)
          IF (rtree%isizeSngl > 0) rtree%p_FData(:,ipos) = rtree%p_FData(:,jpos)
          IF (rtree%isizeInt  > 0) rtree%p_IData(:,ipos) = rtree%p_IData(:,jpos)
          EXIT
        END IF
        ipred = jpos
        jpos = rtree%p_Kchild(TRIGHT,jpos)
      END DO
    END IF

    ! Node to be deleted has less than two childen
    ipos = rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
    
    IF (rtree%p_Kchild(TLEFT,ipos) .EQ. TNULL) THEN
      IF (rtree%p_Kchild(TRIGHT,ipos) .EQ. TNULL) THEN
        ! Node to be deleted is leaf: Nullify pointer to node and
        ! mark position in array as deleted
        rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) = TNULL
        rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
        rtree%p_Kchild(TFREE,TROOT) = -ipos
      ELSE
        ! Node to be deleted has right child: Set pointer to right
        ! child and mark position in array as deleted
        rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) =&
            rtree%p_Kchild(TRIGHT,ipos)
        rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
        rtree%p_Kchild(TFREE,TROOT) = -ipos
      END IF
    ELSE
      ! Node to be deleted has left child: Set pointer to left child
      ! and mark position in array as deleted
      rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) =&
          rtree%p_Kchild(TLEFT,ipos)
      rtree%p_Kchild(TFREE,ipos)  = rtree%p_Kchild(TFREE,TROOT)
      rtree%p_Kchild(TFREE,TROOT) = -ipos
    END IF
    
    ! Invoke rebalance procedure
    IF (rtree%depth > 0) CALL rebalanceAfterDeletion(rtree,rtree%depth)
  END FUNCTION btree_deleteFromTreeInt

  ! ***************************************************************************

!<subroutine>  
  
  RECURSIVE SUBROUTINE rebalanceAfterDeletion(rtree,i)

!<description>
    ! This subroutine rebalances the AVL tree after deletion
!</description>

!<input>
    ! starting node
    INTEGER(PREC_TREEIDX), INTENT(IN) :: i
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: v,x,w
    INTEGER :: dir,xbal
    
    v   = rtree%p_Kpath(i)
    dir = MERGE(TLEFT,TRIGHT,v < 0)
    v   = ABS(v)
    
      ! Which rotation needs to be performed?
    SELECT CASE (dir)
      
    CASE (TLEFT)
      ! Node V has been reached from the left
      SELECT CASE (rtree%p_Kbal(v))
        
      CASE (0) ! bal(v)=0
        rtree%p_Kbal(v) = 1
        
      CASE (-1) ! bal(v)=-1
        rtree%p_Kbal(v) = 0
        IF (v .NE. rtree%p_Kchild(TRIGHT,TROOT))&
            CALL rebalanceAfterDeletion(rtree,i-1)
        
      CASE (1) ! bal(v)=1
        x = rtree%p_Kchild(TRIGHT,v)
        
        SELECT CASE(rtree%p_Kbal(x))
          
        CASE (-1) ! bal(x)=-1
          w = rtree%p_Kchild(TLEFT,x)
          ! Double-Right-Left-rotation
          rtree%p_Kchild(TRIGHT,v) = rtree%p_Kchild(TLEFT,w)
          rtree%p_Kchild(TLEFT,x)  = rtree%p_Kchild(TRIGHT,w)
          rtree%p_Kchild(TLEFT,w)  = v
          rtree%p_Kchild(TRIGHT,w) = x
          rtree%p_Kbal(v)          = -MAX(0,rtree%p_Kbal(w))
          rtree%p_Kbal(x)          = -MIN(0,rtree%p_Kbal(w))
          rtree%p_Kbal(w)          = 0
          
          IF (v .EQ. rtree%p_Kchild(TRIGHT,TROOT)) THEN
            rtree%p_Kchild(TRIGHT,TROOT) = w
          ELSE
            rtree%p_Kchild(MERGE(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                ABS(rtree%p_Kpath(i-1))) = w
            CALL rebalanceAfterDeletion(rtree,i-1)
          END IF
          
        CASE DEFAULT ! bal(x)=0 or bal(x)=1
          ! Single-Left-rotation
          rtree%p_Kchild(TRIGHT,v) = rtree%p_Kchild(TLEFT,x)
          rtree%p_Kchild(TLEFT,x)  = v
          xbal                   = rtree%p_Kbal(x)
          rtree%p_Kbal(x)          = rtree%p_Kbal(x)-1
          rtree%p_Kbal(v)          = -rtree%p_Kbal(x)
          
          IF (v .EQ. rtree%p_Kchild(TRIGHT,TROOT)) THEN
            rtree%p_Kchild(TRIGHT,TROOT) = x
          ELSE
            rtree%p_Kchild(MERGE(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                ABS(rtree%p_Kpath(i-1))) = x
            IF (xbal .EQ. 1) CALL rebalanceAfterDeletion(rtree,i-1)
          END IF
          
        END SELECT
        
      END SELECT
      
    CASE (TRIGHT)
      ! Node V has been reached from the right
      SELECT CASE (rtree%p_Kbal(v))
        
      CASE (0) ! bal(v)=0
        rtree%p_Kbal(v) = -1
        
      CASE (1) ! bal(v)=1
        rtree%p_Kbal(v) = 0
        IF (v .NE. rtree%p_Kchild(TRIGHT,TROOT))&
            CALL rebalanceAfterDeletion(rtree,i-1)
        
      CASE (-1) ! bal(v)=-1
        x=rtree%p_Kchild(TLEFT,v)
        
        SELECT CASE(rtree%p_Kbal(x))
          
        CASE (1) ! bal(x)=1
          w = rtree%p_Kchild(TRIGHT,x)
          ! Double-Left-Right-rotation
          rtree%p_Kchild(TLEFT,v ) = rtree%p_Kchild(TRIGHT,w)
          rtree%p_Kchild(TRIGHT,x) = rtree%p_Kchild(TLEFT,w)
          rtree%p_Kchild(TLEFT,w)  = x
          rtree%p_Kchild(TRIGHT,w) = v
          rtree%p_Kbal(v)          = -MIN(0,rtree%p_Kbal(w))
          rtree%p_Kbal(x)          = -MAX(0,rtree%p_Kbal(w))
          rtree%p_Kbal(w)          = 0
          
          IF (v .EQ. rtree%p_Kchild(TRIGHT,TROOT)) THEN
            rtree%p_Kchild(TRIGHT,TROOT) = w
          ELSE
            rtree%p_Kchild(MERGE(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                ABS(rtree%p_Kpath(i-1))) = w
            CALL rebalanceAfterDeletion(rtree,i-1)
          END IF
          
        CASE DEFAULT ! bal(x)=0 or bal(x)=-1
          ! Single-Right-rotation
          rtree%p_Kchild(TLEFT,v)  = rtree%p_Kchild(TRIGHT,x)
          rtree%p_Kchild(TRIGHT,x) = v
          xbal                   = rtree%p_Kbal(x)
          rtree%p_Kbal(x)          = rtree%p_Kbal(x)+1
          rtree%p_Kbal(v)          = -rtree%p_Kbal(x)
          
          IF (v .EQ. rtree%p_Kchild(TRIGHT,TROOT)) THEN
            rtree%p_Kchild(TRIGHT,TROOT) = x
          ELSE
            rtree%p_Kchild(MERGE(TLEFT,TRIGHT,rtree%p_Kpath(i-1) < 0),&
                ABS(rtree%p_Kpath(i-1))) = x
            IF (xbal .EQ. -1) CALL rebalanceAfterDeletion(rtree,i-1)
          END IF
          
        END SELECT
        
      END SELECT
      
    END SELECT
  END SUBROUTINE rebalanceAfterDeletion

  ! ***************************************************************************
  
!<function>

  FUNCTION btree_searchInTreeDble(rtree,dkey,ipos) RESULT(f)
!<description>
    ! This subroutine searches for a given Double key in the tree and
    ! returns the position of its predecessor.
!</description>

!<input>
    ! Key
    REAL(DP), INTENT(IN) :: dkey
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>

!<output>
    ! Position of the predecessor
    INTEGER(PREC_TREEIDX), INTENT(OUT) :: ipos
!</output>

!<result>
    ! Result of the search BTREE_NOT_FOUND / BTREE_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_TREEIDX) :: jpos
    INTEGER :: dir

    ! Check if list format is ok
    IF (rtree%ctreeFormat .NE. ST_DOUBLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_searchInTreeDble')
      CALL sys_halt()
    END IF
    
    f           = BTREE_NOT_FOUND
    ipos        = TROOT
    dir         = TRIGHT
    jpos        = rtree%p_Kchild(dir,ipos)
    rtree%depth = 0
    
    search: DO
      IF (jpos .EQ. TNULL) THEN
        ipos = MERGE(-ipos,ipos,dir .EQ. TLEFT)
        EXIT search
      END IF
      
      IF (rtree%p_DKey(jpos) .EQ. dkey) THEN
        f    = BTREE_FOUND
        ipos = MERGE(-ipos,ipos,dir .EQ. TLEFT)
        EXIT search
      END IF
      
      ipos                     = jpos
      dir                      = MERGE(TLEFT,TRIGHT,rtree%p_DKey(ipos) > dkey)
      jpos                     = rtree%p_Kchild(dir,ipos)
      rtree%depth              = rtree%depth+1
      rtree%p_Kpath(rtree%depth) = MERGE(-ipos,ipos,dir .EQ. TLEFT)
    END DO search
  END FUNCTION btree_searchInTreeDble

  ! ***************************************************************************
  
!<function>

  FUNCTION btree_searchInTreeSngl(rtree,skey,ipos) RESULT(f)
!<description>
    ! This subroutine searches for a given Single key in the tree and
    ! returns the position of its predecessor.
!</description>

!<input>
    ! Key
    REAL(SP), INTENT(IN) :: skey
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>

!<output>
    ! Position of the predecessor
    INTEGER(PREC_TREEIDX), INTENT(OUT) :: ipos
!</output>

!<result>
    ! Result of the search BTREE_NOT_FOUND / BTREE_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_TREEIDX) :: jpos
    INTEGER :: dir

    ! Check if list format is ok
    IF (rtree%ctreeFormat .NE. ST_SINGLE) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_searchInTreeSngl')
      CALL sys_halt()
    END IF
    
    f           = BTREE_NOT_FOUND
    ipos        = TROOT
    dir         = TRIGHT
    jpos        = rtree%p_Kchild(dir,ipos)
    rtree%depth = 0
    
    search: DO
      IF (jpos .EQ. TNULL) THEN
        ipos = MERGE(-ipos,ipos,dir .EQ. TLEFT)
        EXIT search
      END IF
      
      IF (rtree%p_FKey(jpos) .EQ. skey) THEN
        f    = BTREE_FOUND
        ipos = MERGE(-ipos,ipos,dir .EQ. TLEFT)
        EXIT search
      END IF
      
      ipos                     = jpos
      dir                      = MERGE(TLEFT,TRIGHT,rtree%p_FKey(ipos) > skey)
      jpos                     = rtree%p_Kchild(dir,ipos)
      rtree%depth              = rtree%depth+1
      rtree%p_Kpath(rtree%depth) = MERGE(-ipos,ipos,dir .EQ. TLEFT)
    END DO search
  END FUNCTION btree_searchInTreeSngl

  ! ***************************************************************************
  
!<function>

  FUNCTION btree_searchInTreeInt(rtree,ikey,ipos) RESULT(f)
!<description>
    ! This subroutine searches for a given Integer key in the tree and
    ! returns the position of its predecessor.
!</description>

!<input>
    ! Key
    INTEGER(PREC_TREEIDX), INTENT(IN) :: ikey
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>

!<output>
    ! Position of the predecessor
    INTEGER(PREC_TREEIDX), INTENT(OUT) :: ipos
!</output>

!<result>
    ! Result of the search BTREE_NOT_FOUND / BTREE_FOUND
    INTEGER :: f
!</result>
!</function>

    ! local variables
    INTEGER(PREC_TREEIDX) :: jpos
    INTEGER :: dir

    ! Check if list format is ok
    IF (rtree%ctreeFormat .NE. ST_INT) THEN
      CALL output_line('Unsupported data format!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_searchInTreeInt')
      CALL sys_halt()
    END IF
    
    f           = BTREE_NOT_FOUND
    ipos        = TROOT
    dir         = TRIGHT
    jpos        = rtree%p_Kchild(dir,ipos)
    rtree%depth = 0
    
    search: DO
      IF (jpos .EQ. TNULL) THEN
        ipos = MERGE(-ipos,ipos,dir .EQ. TLEFT)
        EXIT search
      END IF
      
      IF (rtree%p_IKey(jpos) .EQ. ikey) THEN
        f    = BTREE_FOUND
        ipos = MERGE(-ipos,ipos,dir .EQ. TLEFT)
        EXIT search
      END IF
      
      ipos                     = jpos
      dir                      = MERGE(TLEFT,TRIGHT,rtree%p_IKey(ipos) > ikey)
      jpos                     = rtree%p_Kchild(dir,ipos)
      rtree%depth              = rtree%depth+1
      rtree%p_Kpath(rtree%depth) = MERGE(-ipos,ipos,dir .EQ. TLEFT)
    END DO search
  END FUNCTION btree_searchInTreeInt

  ! ***************************************************************************

!<function>

  FUNCTION btree_getItemInTreeDble(rtree,dkey) RESULT(ipos)

!<description>
    ! This subroutine searches for a given Double key in the tree and
    ! returns its position. If the item cannot not be found than 
    ! program execution is terminated.
!</description>

!<input>
    ! Key
    REAL(DP), INTENT(IN) :: dkey
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>

!<result>
    ! Position of the item
    INTEGER(PREC_TREEIDX) :: ipos
!</result>
!</function>
  
    ! local variables
    INTEGER(PREC_TREEIDX) :: ipred
    
    IF (btree_searchInTree(rtree,dkey,ipred) .NE. BTREE_FOUND) THEN
      CALL output_line('Unable to find item in tree!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_getItemInTreeDble')
      CALL sys_halt()
    END IF
    ipos = rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
  END FUNCTION btree_getItemInTreeDble

  ! ***************************************************************************

!<function>

  FUNCTION btree_getItemInTreeSngl(rtree,skey) RESULT(ipos)

!<description>
    ! This subroutine searches for a given Single key in the tree and
    ! returns its position. If the item cannot not be found than 
    ! program execution is terminated.
!</description>

!<input>
    ! Key
    REAL(SP), INTENT(IN) :: skey
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>

!<result>
    ! Position of the item
    INTEGER(PREC_TREEIDX) :: ipos
!</result>
!</function>
  
    ! local variables
    INTEGER(PREC_TREEIDX) :: ipred
    
    IF (btree_searchInTree(rtree,skey,ipred) .NE. BTREE_FOUND) THEN
      CALL output_line('Unable to find item in tree!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_getItemInTreeSngl')
      CALL sys_halt()
    END IF
    ipos = rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
  END FUNCTION btree_getItemInTreeSngl
  
  ! ***************************************************************************

!<function>

  FUNCTION btree_getItemInTreeInt(rtree,ikey) RESULT(ipos)

!<description>
    ! This subroutine searches for a given Integer key in the tree and
    ! returns its position. If the item cannot not be found than 
    ! program execution is terminated.
!</description>

!<input>
    ! Key
    INTEGER(PREC_TREEIDX), INTENT(IN) :: ikey
!</input>

!<inputoutput>
    ! tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>

!<result>
    ! Position of the item
    INTEGER(PREC_TREEIDX) :: ipos
!</result>
!</function>
  
    ! local variables
    INTEGER(PREC_TREEIDX) :: ipred
    
    IF (btree_searchInTree(rtree,ikey,ipred) .NE. BTREE_FOUND) THEN
      CALL output_line('Unable to find item in tree!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_getItemInTreeInt')
      CALL sys_halt()
    END IF
    ipos = rtree%p_Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
  END FUNCTION btree_getItemInTreeInt

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE btree_printTree(rtree,op)

!<description>
    ! This subroutine prints the content of the tree
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree

    ! type of traversal
    INTEGER, INTENT(IN) :: op
!</input>
!</subroutine>

    ! Which kind of traversal should be applied
    SELECT CASE (op)
      
    CASE (BTREE_PREORDER)
      IF (rtree%p_Kchild(TRIGHT,TROOT) .NE. TNULL) THEN
        SELECT CASE(rtree%ctreeFormat)
        CASE (ST_DOUBLE)
          CALL preorderDble(rtree%p_Kchild(TRIGHT,TROOT))
        CASE (ST_SINGLE)
          CALL preorderSngl(rtree%p_Kchild(TRIGHT,TROOT))
        CASE (ST_INT)
          CALL preorderInt(rtree%p_Kchild(TRIGHT,TROOT))
        CASE DEFAULT
          CALL output_line('Unsupported data format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'btree_printTree')
          CALL sys_halt()
        END SELECT
      END IF
      
    CASE (BTREE_INORDER)
      IF (rtree%p_Kchild(TRIGHT,TROOT) .NE. TNULL) THEN
        SELECT CASE(rtree%ctreeFormat)
        CASE (ST_DOUBLE)
          CALL inorderDble(rtree%p_Kchild(TRIGHT,TROOT))
        CASE (ST_SINGLE)
          CALL inorderSngl(rtree%p_Kchild(TRIGHT,TROOT))
        CASE (ST_INT)
          CALL inorderInt(rtree%p_Kchild(TRIGHT,TROOT))
        CASE DEFAULT
          CALL output_line('Unsupported data format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'btree_printTree')
          CALL sys_halt()
        END SELECT
      END IF
      
    CASE (BTREE_POSTORDER)
      IF (rtree%p_Kchild(TRIGHT,TROOT) .NE. TNULL) THEN
        SELECT CASE(rtree%ctreeFormat)
        CASE (ST_DOUBLE)
          CALL postorderDble(rtree%p_Kchild(TRIGHT,TROOT))
        CASE (ST_SINGLE)
          CALL postorderSngl(rtree%p_Kchild(TRIGHT,TROOT))
        CASE (ST_INT)
          CALL postorderInt(rtree%p_Kchild(TRIGHT,TROOT))
        CASE DEFAULT
          CALL output_line('Unsupported data format!',&
              OU_CLASS_ERROR,OU_MODE_STD,'btree_printTree')
          CALL sys_halt()
        END SELECT
      END IF
    END SELECT
  CONTAINS

    ! Here, the real working routines follow.

    !**************************************************************
    ! Print the content of the tree in preorder for Double key
    
    RECURSIVE SUBROUTINE preorderDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      WRITE(*,FMT='(A)',ADVANCE='NO') rtree%p_DKey(i),','
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL preorderDble(rtree%p_Kchild(TLEFT,i))
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL preorderDble(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE preorderDble

    !**************************************************************
    ! Print the content of the tree in preorder for Single key
    
    RECURSIVE SUBROUTINE preorderSngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      WRITE(*,FMT='(A)',ADVANCE='NO') rtree%p_FKey(i),','
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL preorderSngl(rtree%p_Kchild(TLEFT,i))
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL preorderSngl(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE preorderSngl

    !**************************************************************
    ! Print the content of the tree in preorder for Integer key
    
    RECURSIVE SUBROUTINE preorderInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      WRITE(*,FMT='(A)',ADVANCE='NO') rtree%p_IKey(i),','
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL preorderInt(rtree%p_Kchild(TLEFT,i))
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL preorderInt(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE preorderInt

    !**************************************************************
    ! Print the content of the tree in postorder for Double key
    
    RECURSIVE SUBROUTINE postorderDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL postorderDble(rtree%p_Kchild(TLEFT,i))
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL postorderDble(rtree%p_Kchild(TRIGHT,i))
      WRITE(*,FMT='(A)',ADVANCE='NO') rtree%p_DKey(i),','
    END SUBROUTINE postorderDble

    !**************************************************************
    ! Print the content of the tree in postorder for Single key
    
    RECURSIVE SUBROUTINE postorderSngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL postorderSngl(rtree%p_Kchild(TLEFT,i))
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL postorderSngl(rtree%p_Kchild(TRIGHT,i))
      WRITE(*,FMT='(A)',ADVANCE='NO') rtree%p_FKey(i),','
    END SUBROUTINE postorderSngl

    !**************************************************************
    ! Print the content of the tree in postorder for Integer key
    
    RECURSIVE SUBROUTINE postorderInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL postorderInt(rtree%p_Kchild(TLEFT,i))
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL postorderInt(rtree%p_Kchild(TRIGHT,i))
      WRITE(*,FMT='(A)',ADVANCE='NO') rtree%p_IKey(i),','
    END SUBROUTINE postorderInt

    !**************************************************************
    ! Print the content of the tree in inorder for Double key
    
    RECURSIVE SUBROUTINE inorderDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderDble(rtree%p_Kchild(TLEFT,i))
      WRITE(*,FMT='(A)',ADVANCE='NO') rtree%p_DKey(i),','
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderDble(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderDble

    !**************************************************************
    ! Print the content of the tree in inorder for Single key
    
    RECURSIVE SUBROUTINE inorderSngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderSngl(rtree%p_Kchild(TLEFT,i))
      WRITE(*,FMT='(A)',ADVANCE='NO') rtree%p_FKey(i),','
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderSngl(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderSngl

    !**************************************************************
    ! Print the content of the tree in inorder for Integer key
    
    RECURSIVE SUBROUTINE inorderInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderInt(rtree%p_Kchild(TLEFT,i))
      WRITE(*,FMT='(A)',ADVANCE='NO') rtree%p_IKey(i),','
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderInt(rtree%p_Kchild(TRIGHT,i))
    END SUBROUTINE inorderInt
  END SUBROUTINE btree_printTree

  ! ***************************************************************************
  
!<function>
  
  PURE FUNCTION btree_getHeight(rtree) RESULT(h)

!<description>
    ! This function computes the height of a given tree
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
!</input>

!<result>
    ! height of the tree
    INTEGER(PREC_TREEIDX) :: h
!</result>
!</function>
    
    h = height(rtree%p_Kchild(TRIGHT,TROOT))
    
  CONTAINS
    
    PURE RECURSIVE FUNCTION height(i) RESULT(h)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      INTEGER(PREC_TREEIDX) :: h,hl,hr
      
      IF (rtree%p_Kchild(TLEFT,i) .NE. TNULL) THEN
        hl = height(rtree%p_Kchild(TLEFT,i))
      ELSE
        hl = 0
      END IF
      
      IF (rtree%p_Kchild(TRIGHT,i) .NE. TNULL) THEN
        hr = height(rtree%p_Kchild(TRIGHT,i))
      ELSE
        hr = 0
      END IF
      
      h = MAX(hl,hr)+1
    END FUNCTION height
  END FUNCTION btree_getHeight

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE btree_infoTree(rtree)

!<description>
    ! This subroutine outputs statistical info about the tree
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
!</input>
!</subroutine>

    CALL output_line('Tree statistics:')
    CALL output_line('----------------')
    CALL output_line('ctreeFormat: '//TRIM(sys_siL(rtree%ctreeFormat,5)))
    CALL output_line('NA:          '//TRIM(sys_siL(rtree%NA,15)))
    CALL output_line('NNA:         '//TRIM(sys_siL(rtree%NNA,15)))
    CALL output_line('NNA0:        '//TRIM(sys_siL(rtree%NNA0,15)))
    CALL output_line('NRESIZE:     '//TRIM(sys_siL(rtree%NRESIZE,15)))
    CALL output_line('h_Key:       '//TRIM(sys_siL(rtree%h_Key,15)))
    CALL output_line('h_Kbal:      '//TRIM(sys_siL(rtree%h_Kbal,15)))
    CALL output_line('h_Kpath:     '//TRIM(sys_siL(rtree%h_Kpath,15)))
    CALL output_line('h_Kchild:    '//TRIM(sys_siL(rtree%h_Kchild,15)))
    CALL output_line('h_DData:     '//TRIM(sys_siL(rtree%h_DData,15)))
    CALL output_line('h_FData:     '//TRIM(sys_siL(rtree%h_FData,15)))
    CALL output_line('h_IData:     '//TRIM(sys_siL(rtree%h_IData,15)))
    CALL output_lbrk()
    CALL output_line('Current data  memory usage: '//&
        TRIM(sys_sdL(100*rtree%NA/REAL(rtree%NNA,DP),2))//'%')

  END SUBROUTINE btree_infoTree

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_duplicateTree(rtree,rtreeBackup)

!<description>
    ! This subroutine makes a copy of a binary tree in memory.
    ! It does not make sense to share some information between binary trees,
    ! so each vectors is physically copied from the source tree to the destination tree.
!</description>

!<input>
    ! Source tree
    TYPE(t_btree), INTENT(IN) :: rtree
!</input>

!<inputoutput>
    ! Destination tree
    TYPE(t_btree), INTENT(INOUT) :: rtreeBackup
!</inputoutput>
!</subroutine>

    ! Release backup tree
    CALL btree_releaseTree(rtreeBackup)

    ! Copy all data
    rtreeBackup = rtree

    ! Reset Handles
    rtreeBackup%h_Key    = ST_NOHANDLE
    rtreeBackup%h_Kbal   = ST_NOHANDLE
    rtreeBackup%h_Kpath  = ST_NOHANDLE
    rtreeBackup%h_Kchild = ST_NOHANDLE
    rtreeBackup%h_DData  = ST_NOHANDLE
    rtreeBackup%h_FData  = ST_NOHANDLE
    rtreeBackup%h_IData  = ST_NOHANDLE

    ! Copy storage blocks
    IF (rtree%h_Key .NE. ST_NOHANDLE) THEN
      CALL storage_copy(rtree%h_Key,rtreeBackup%h_Key)
      SELECT CASE(rtreeBackup%ctreeFormat)
      CASE(ST_DOUBLE)
        CALL storage_getbase_double(rtreeBackup%h_Key,rtreeBackup%p_DKey)
      CASE(ST_SINGLE)
        CALL storage_getbase_single(rtreeBackup%h_Key,rtreeBackup%p_FKey)
      CASE(ST_INT)
        CALL storage_getbase_int(rtreeBackup%h_Key,rtreeBackup%p_IKey)
      CASE DEFAULT
        CALL output_line('Unsupported data format!',&
            OU_CLASS_ERROR,OU_MODE_STD,'btree_duplicateTree')
        CALL sys_halt()
      END SELECT
    END IF
          
    IF (rtree%h_Kbal .NE. ST_NOHANDLE) THEN
      CALL storage_copy(rtree%h_Kbal,rtreeBackup%h_Kbal)
      CALL storage_getbase_int(rtreeBackup%h_Kbal,rtreeBackup%p_Kbal)
    END IF

    IF (rtree%h_Kpath .NE. ST_NOHANDLE) THEN
      CALL storage_copy(rtree%h_Kpath,rtreeBackup%h_Kpath)
      CALL storage_getbase_int(rtreeBackup%h_Kpath,rtreeBackup%p_Kpath)
    END IF

    IF (rtree%h_Kchild .NE. ST_NOHANDLE) THEN
      CALL storage_copy(rtree%h_Kchild,rtreeBackup%h_Kchild)
      CALL storage_getbase_int2D(rtreeBackup%h_Kchild,rtreeBackup%p_Kchild)
    END IF

    IF (rtree%h_DData .NE. ST_NOHANDLE) THEN
      CALL storage_copy(rtree%h_DData,rtreeBackup%h_DData)
      CALL storage_getbase_double2D(rtreeBackup%h_DData,rtreeBackup%p_DData)
    END IF

    IF (rtree%h_FData .NE. ST_NOHANDLE) THEN
      CALL storage_copy(rtree%h_FData,rtreeBackup%h_Fdata)
      CALL storage_getbase_single2D(rtreeBackup%h_FData,rtreeBackup%p_FData)
    END IF

    IF (rtree%h_IData .NE. ST_NOHANDLE) THEN
      CALL storage_copy(rtree%h_IData,rtreeBackup%h_IData)
      CALL storage_getbase_int2D(rtreeBackup%h_IData,rtreeBackup%p_IData)
    END IF
  END SUBROUTINE btree_duplicateTree

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE btree_restoreTree(rtreeBackup,rtree)

!<description>
    ! This subroutine restores a binary tree from a previous backup.
    ! The format of both trees must be the same.
!</description>

!<input>
    ! Backup of a binary tree
    TYPE(t_btree), INTENT(IN) :: rtreeBackup
!</input>

!<inputoutput>
    ! Destination binary tree
    TYPE(t_btree), INTENT(INOUT) :: rtree
!</inputoutput>
!</subroutine>

    ! Check that both trees are compatible
    IF (rtree%ctreeFormat .NE. rtreeBackup%ctreeFormat .OR.&
        rtree%isizeDble   .NE. rtreeBackup%isizeDble   .OR.&
        rtree%isizeSngl   .NE. rtreeBackup%isizeSngl   .OR.&
        rtree%isizeInt    .NE. rtreeBackup%isizeInt) THEN
      CALL output_line('Incompatible binary trees!',&
          OU_CLASS_ERROR,OU_MODE_STD,'btree_restoreTree')
      CALL sys_halt()
    END IF

    ! Release binary tree
    CALL btree_releaseTree(rtree)

    ! Duplicate the backup
    CALL btree_duplicateTree(rtreeBackup,rtree)
  END SUBROUTINE btree_restoreTree

  ! ***************************************************************************
  
!<function>
  
  ELEMENTAL FUNCTION LOG2(i)

!<description>
    ! Compute log_2(i)
!</description>

!<input>
    ! integer value
    INTEGER, INTENT(IN) :: i
!</input>

!<result>
    ! log_2(i)
    REAL(DP) :: log2
!</result>
!</function>
        
    LOG2=LOG(REAL(i,DP))/LOG(2._DP)
  END FUNCTION LOG2
END MODULE binarytree

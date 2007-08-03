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
!# 4.) btree_copyToTree = t_tree_copyto_handle /
!#                       t_tree_copyto_arrayDble /
!#                       t_tree_copyto_arraySngl / 
!#                       t_tree_copyto_arrayInt
!#     -> Copy key and auxiliary data to tree
!#
!# 5.) btree_copyFromTreeKey = t_tree_copyfrom_key_handle /
!#                            t_tree_copyfrom_key_arrayDble /
!#                            t_tree_copyfrom_key_arraySngl /
!#                            t_tree_copyfrom_key_arrayInt
!#     -> Copy the key of the tree to handle/array
!#
!# 6.) btree_copyFromTreeDble = t_tree_copyfrom_handle /
!#                             t_tree_copyfrom_arrayDble
!#                             t_tree_copyfrom_arrayDble2D
!#     -> Copy some auxiliary double data of the tree to handle/array
!#
!# 7.) btree_copyFromTreeSngl = t_tree_copyfrom_handle /
!#                             t_tree_copyfrom_arraySngl
!#                             t_tree_copyfrom_arraySngl2D
!#     -> Copy some auxiliary single data of the tree to handle/array
!#
!# 8.) btree_copyFromTreeInt = t_tree_copyfrom_handle /
!#                            t_tree_copyfrom_arrayInt
!#                            t_tree_copyfrom_arrayInt2D
!#     -> Copy some auxiliary integer data of the tree to handle/array
!#
!# 9.) btree_insertToTree = t_tree_insertDble /
!#                         t_tree_insertSngl /
!#                         t_tree_insertInt
!#     -> Insert key into tree
!#
!# 10.) btree_deleteFromTree = t_tree_deleteDble /
!#                            t_tree_deleteSngl /
!#                            t_tree_deleteInt
!#      -> Delete key from tree
!#
!# 11.) btree_searchInTree = t_tree_searchDble /
!#                          t_tree_searchSngl /
!#                          t_tree_searchInt
!#      -> Search for key in tree
!#
!# 12.) btree_getItemInTree = t_tree_getitemDble /
!#                           t_tree_getitemSngl / 
!#                           t_tree_getitemInt
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
!# </purpose>
!##############################################################################
MODULE binarytree

  USE fsystem
  USE storage
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
    INTEGER(PREC_TREEIDX), DIMENSION(:,:), POINTER :: Kchild => NULL()

    ! Tree path structure
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_TREEIDX), DIMENSION(:), POINTER :: Kpath => NULL()

    ! Tree balance structure
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER, DIMENSION(:), POINTER :: Kbal => NULL()

    ! Tree key data (Double)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(DP), DIMENSION(:), POINTER :: DKey => NULL()

    ! Tree key data (Single)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(SP), DIMENSION(:), POINTER :: FKey => NULL()

    ! Tree key data (Integer)
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_TREEIDX), DIMENSION(:), POINTER :: IKey => NULL()

    ! Tree data (Double)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(DP), DIMENSION(:,:), POINTER :: DData => NULL()

    ! Tree data (Single)
    ! NOTE: This array is introduced to increase performance (see above).
    REAL(SP), DIMENSION(:,:), POINTER :: FData => NULL()

    ! Tree data (Integer)
    ! NOTE: This array is introduced to increase performance (see above).
    INTEGER(PREC_TREEIDX), DIMENSION(:,:), POINTER :: IData => NULL()
  END TYPE t_btree

!</typeblock>
!</types>

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************
  
  
  INTERFACE btree_createTree
    MODULE PROCEDURE t_tree_create
  END INTERFACE
  
  INTERFACE btree_releaseTree
    MODULE PROCEDURE t_tree_release
  END INTERFACE
  
  INTERFACE btree_resizeTree
    MODULE PROCEDURE t_tree_resize
  END INTERFACE
  
  INTERFACE resize   ! internal use
    MODULE PROCEDURE t_tree_resize
  END INTERFACE

  INTERFACE btree_copyToTree
    MODULE PROCEDURE t_tree_copyto_handle
    MODULE PROCEDURE t_tree_copyto_arrayDble
    MODULE PROCEDURE t_tree_copyto_arraySngl
    MODULE PROCEDURE t_tree_copyto_arrayInt
  END INTERFACE

  INTERFACE btree_copyFromTreeKey
    MODULE PROCEDURE t_tree_copyfrom_key_handle
    MODULE PROCEDURE t_tree_copyfrom_key_arrayDble
    MODULE PROCEDURE t_tree_copyfrom_key_arraySngl
    MODULE PROCEDURE t_tree_copyfrom_key_arrayInt
  END INTERFACE

  INTERFACE btree_copyFromTreeDble
    MODULE PROCEDURE t_tree_copyfrom_handle
    MODULE PROCEDURE t_tree_copyfrom_arrayDble
    MODULE PROCEDURE t_tree_copyfrom_arrayDble2D
  END INTERFACE

  INTERFACE btree_copyFromTreeSngl
    MODULE PROCEDURE t_tree_copyfrom_handle
    MODULE PROCEDURE t_tree_copyfrom_arraySngl
    MODULE PROCEDURE t_tree_copyfrom_arraySngl2D
  END INTERFACE

  INTERFACE btree_copyFromTreeInt
    MODULE PROCEDURE t_tree_copyfrom_handle
    MODULE PROCEDURE t_tree_copyfrom_arrayInt
    MODULE PROCEDURE t_tree_copyfrom_arrayInt2D
  END INTERFACE
  
  INTERFACE btree_insertIntoTree
    MODULE PROCEDURE t_tree_insertDble
    MODULE PROCEDURE t_tree_insertSngl
    MODULE PROCEDURE t_tree_insertInt
  END INTERFACE

  INTERFACE insert   ! internal use
    MODULE PROCEDURE t_tree_insertDble
    MODULE PROCEDURE t_tree_insertSngl
    MODULE PROCEDURE t_tree_insertInt
  END INTERFACE

  INTERFACE btree_deleteFromTree
    MODULE PROCEDURE t_tree_deleteDble
    MODULE PROCEDURE t_tree_deleteSngl
    MODULE PROCEDURE t_tree_deleteInt
  END INTERFACE

  INTERFACE delete   ! internal use
    MODULE PROCEDURE t_tree_deleteDble
    MODULE PROCEDURE t_tree_deleteSngl
    MODULE PROCEDURE t_tree_deleteInt
  END INTERFACE
  
  INTERFACE btree_searchInTree
    MODULE PROCEDURE t_tree_searchDble
    MODULE PROCEDURE t_tree_searchSngl
    MODULE PROCEDURE t_tree_searchInt
  END INTERFACE

  INTERFACE btree_getItemInTree
    MODULE PROCEDURE t_tree_getitemDble
    MODULE PROCEDURE t_tree_getitemSngl
    MODULE PROCEDURE t_tree_getitemInt
  END INTERFACE

  INTERFACE search   ! internal use
    MODULE PROCEDURE t_tree_searchDble
    MODULE PROCEDURE t_tree_searchSngl
    MODULE PROCEDURE t_tree_searchInt
  END INTERFACE
  
  INTERFACE btree_printTree
    MODULE PROCEDURE t_tree_print
  END INTERFACE
  
  INTERFACE btree_getHeight
    MODULE PROCEDURE t_tree_height
  END INTERFACE

  INTERFACE btree_infoTree
    MODULE PROCEDURE t_tree_info
  END INTERFACE
  

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_create(rtree,nna,ctreeFormat,isizeDble,isizeSngl,isizeInt,dfactor)

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
    CALL storage_new('t_tree_create','Kchild',Isize1,Isize2,&
        ST_INT,rtree%h_Kchild,ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int2D(rtree%h_Kchild,rtree%Kchild)
    
    CALL storage_new('t_tree_create','Kpath',CEILING(1.441*LOG2(nna))+1,&
        ST_INT,rtree%h_Kpath,ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int(rtree%h_Kpath,rtree%Kpath)

    CALL storage_new('t_tree_create','Kbal',nna,ST_INT,rtree%h_Kbal,ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int(rtree%h_Kbal,rtree%Kbal)

    ! Allocate memory for Key
    SELECT CASE(rtree%ctreeFormat)
    CASE (ST_DOUBLE)
      CALL storage_new('t_tree_create','Key',nna,ST_DOUBLE,rtree%h_Key,ST_NEWBLOCK_ZERO)
      CALL storage_getbase_double(rtree%h_Key,rtree%DKey)
    CASE (ST_SINGLE)
      CALL storage_new('t_tree_create','Key',nna,ST_SINGLE,rtree%h_Key,ST_NEWBLOCK_ZERO)
      CALL storage_getbase_single(rtree%h_Key,rtree%FKey)
    CASE (ST_INT)
      CALL storage_new('t_tree_create','Key',nna,ST_INT,rtree%h_Key,ST_NEWBLOCK_ZERO)
      CALL storage_getbase_int(rtree%h_Key,rtree%IKey)
    CASE DEFAULT
      PRINT *, 't_tree_create: Unsupported data format!'
      CALL sys_halt()
    END SELECT

    ! Allocate memory fo auxiliary data
    IF (isizeDble > 0) THEN
      Isize1 = (/isizeDble,nna/)
      CALL storage_new('t_tree_create','DData',Isize1,ST_DOUBLE,rtree%h_DData,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_double2D(rtree%h_DData,rtree%DData)
    END IF

    IF (isizeSngl > 0) THEN
      Isize1 = (/isizeSngl,nna/)
      CALL storage_new('t_tree_create','FData',Isize1,ST_SINGLE,rtree%h_FData,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_single2D(rtree%h_FData,rtree%FData)
    END IF

    IF (isizeInt > 0) THEN
      Isize1 = (/isizeInt,nna/)
      CALL storage_new('t_tree_create','IData',Isize1,ST_INT,rtree%h_IData,ST_NEWBLOCK_NOINIT)
      CALL storage_getbase_int2D(rtree%h_IData,rtree%IData)
    END IF

    ! Initialize tree structure
    rtree%Kchild(TFREE,TROOT) = 1
  END SUBROUTINE t_tree_create

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE t_tree_release(rtree)

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

    NULLIFY(rtree%Dkey,rtree%FKey,rtree%IKey)
    NULLIFY(rtree%Kbal,rtree%Kpath,rtree%Kchild)
    NULLIFY(rtree%DData,rtree%FData,rtree%IData)
  END SUBROUTINE t_tree_release

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE t_tree_resize(rtree,nna)

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

    CALL storage_realloc('t_tree_resize',TROOT,nna,rtree%h_Kchild,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_getbase_int2D(rtree%h_Kchild,rtree%Kchild)

    CALL storage_realloc('t_tree_resize',nna,rtree%h_Kbal,ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_getbase_int(rtree%h_Kbal,rtree%Kbal)

    CALL storage_realloc('t_tree_resize',CEILING(1.441*LOG2(nna))+1,rtree%h_Kpath,&
        ST_NEWBLOCK_ZERO,.TRUE.)
    CALL storage_getbase_int(rtree%h_Kpath,rtree%Kpath)

    CALL storage_realloc('t_tree_resize',nna,rtree%h_Key,ST_NEWBLOCK_ZERO,.TRUE.)
    SELECT CASE(rtree%ctreeFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(rtree%h_Key,rtree%DKey)
    CASE (ST_SINGLE)
      CALL storage_getbase_single(rtree%h_Key,rtree%FKey)
    CASE (ST_INT)
      CALL storage_getbase_int(rtree%h_Key,rtree%IKey)
    CASE DEFAULT
      PRINT *, 't_tree_resize: Unsupported key data format!'
      CALL sys_halt()
    END SELECT

    ! Auxiliary data
    IF (rtree%isizeDble > 0) THEN
      CALL storage_realloc('t_tree_resize',nna,rtree%h_DData,ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_getbase_double2D(rtree%h_DData,rtree%DData)
    END IF

    IF (rtree%isizeSngl > 0) THEN
      CALL storage_realloc('t_tree_resize',nna,rtree%h_FData,ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_getbase_single2D(rtree%h_FData,rtree%FData)
    END IF

    IF (rtree%isizeInt > 0) THEN
      CALL storage_realloc('t_tree_resize',nna,rtree%h_IData,ST_NEWBLOCK_NOINIT,.TRUE.)
      CALL storage_getbase_int2D(rtree%h_IData,rtree%IData)
    END IF
  END SUBROUTINE t_tree_resize

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyto_handle(h_Key,rtree,h_DData,h_FData,h_IData)

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
        CALL insert(rtree,p_DKey(j))
      END DO
      
    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_Key,p_FKey)
      DO j=1,SIZE(p_FKey)
        CALL insert(rtree,p_FKey(j))
      END DO
      
    CASE (ST_INT)
      CALL storage_getbase_int(h_Key,p_IKey)
      DO j=1,SIZE(p_IKey)
        CALL insert(rtree,p_IKey(j))
      END DO
      
    CASE DEFAULT
      PRINT *, 't_tree_copyto_handle: Unsupported data format!'
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
  END SUBROUTINE t_tree_copyto_handle

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyto_arrayDble(p_DKey,rtree,p_DData,p_FData,p_IData)

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
      PRINT *, "t_tree_copyto_arrayDble: Invalid data format!"
      CALL sys_halt()
    END IF

    ! Build-up treee structure
    DO j=1,SIZE(p_DKey)
      CALL insert(rtree,p_DKey(j))
    END DO

    ! Copy auxiliary arrays to the tree.
    IF (PRESENT(p_DData) .AND. rtree%isizeDble > 0) THEN
      CALL storage_getsize(rtree%h_DData,Isize)
      IF ((rtree%isizeDble .NE. SIZE(p_DData,1)) .OR. (rtree%NA .NE. SIZE(p_DData,2))) THEN
        PRINT *, "t_tree_copyto_arrayDble: Invalid auxiliary data!"
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('t_tree_copyto_arrayDble',rtree%NNA,rtree%h_DData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_double2D(rtree%h_DData,rtree%DData)
      END IF
      ! Copy data
      rtree%DData(:,1:rtree%NA) = p_DData
    END IF

    IF (PRESENT(p_FData) .AND. rtree%isizeSngl > 0) THEN
      CALL storage_getsize(rtree%h_FData,Isize)
      IF ((rtree%isizeSngl .NE. SIZE(p_FData,1)) .OR. (rtree%NA .NE. SIZE(p_FData,2))) THEN
        PRINT *, "t_tree_copyto_arrayDble: Invalid auxiliary data!"
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('t_tree_copyto_arrayDble',rtree%NNA,rtree%h_FData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_single2D(rtree%h_FData,rtree%FData)
      END IF
      ! Copy data
      rtree%FData(:,1:rtree%NA) = p_FData
    END IF

    IF (PRESENT(p_IData) .AND. rtree%isizeInt > 0) THEN
      CALL storage_getsize(rtree%h_IData,Isize)
      IF ((rtree%isizeInt .NE. SIZE(p_IData,1)) .OR. (rtree%NA .NE. SIZE(p_IData,2))) THEN
        PRINT *, "t_tree_copyto_arrayDble: Invalid auxiliary data!"
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('t_tree_copyto_arrayDble',rtree%NNA,rtree%h_IData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_int2D(rtree%h_IData,rtree%IData)
      END IF
      ! Copy data
      rtree%IData(:,1:rtree%NA) = p_IData
    END IF
    
  END SUBROUTINE t_tree_copyto_arrayDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyto_arraySngl(p_FKey,rtree,p_DData,p_FData,p_IData)

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
      PRINT *, "t_tree_copyto_arraySngl: Invalid data format!"
      CALL sys_halt()
    END IF

    ! Build-up treee structure
    DO j=1,SIZE(p_FKey)
      CALL insert(rtree,p_FKey(j))
    END DO

    ! Copy auxiliary arrays to the tree.
    IF (PRESENT(p_DData) .AND. rtree%isizeDble > 0) THEN
      CALL storage_getsize(rtree%h_DData,Isize)
      IF ((rtree%isizeDble .NE. SIZE(p_DData,1)) .OR. (rtree%NA .NE. SIZE(p_DData,2))) THEN
        PRINT *, "t_tree_copyto_arraySngl: Invalid auxiliary data!"
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('t_tree_copyto_arraySngl',rtree%NNA,rtree%h_DData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_double2D(rtree%h_DData,rtree%DData)
      END IF
      ! Copy data
      rtree%DData(:,1:rtree%NA) = p_DData
    END IF
    
    IF (PRESENT(p_FData) .AND. rtree%isizeSngl > 0) THEN
      CALL storage_getsize(rtree%h_FData,Isize)
      IF ((rtree%isizeSngl .NE. SIZE(p_FData,1)) .OR. (rtree%NA .NE. SIZE(p_FData,2))) THEN
        PRINT *, "t_tree_copyto_arraySngl: Invalid auxiliary data!"
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('t_tree_copyto_arraySngl',rtree%NNA,rtree%h_FData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_single2D(rtree%h_FData,rtree%FData)
      END IF
      ! Copy data
      rtree%FData(:,1:rtree%NA) = p_FData
    END IF

    IF (PRESENT(p_IData) .AND. rtree%isizeInt > 0) THEN
      CALL storage_getsize(rtree%h_IData,Isize)
      IF ((rtree%isizeInt .NE. SIZE(p_IData,1)) .OR. (rtree%NA .NE. SIZE(p_IData,2))) THEN
        PRINT *, "t_tree_copyto_arraySngl: Invalid auxiliary data!"
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('t_tree_copyto_arraySngl',rtree%NNA,rtree%h_IData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_int2D(rtree%h_IData,rtree%IData)
      END IF
      ! Copy data
      rtree%IData(:,1:rtree%NA) = p_IData
    END IF

  END SUBROUTINE t_tree_copyto_arraySngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyto_arrayInt(p_IKey,rtree,p_DData,p_FData,p_IData)

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
      PRINT *, "t_tree_copyto_arrayInt: Invalid data format!"
      CALL sys_halt()
    END IF

    ! Build-up treee structure
    DO j=1,SIZE(p_IKey)
      CALL insert(rtree,p_IKey(j))
    END DO

    ! Copy auxiliary arrays to the tree.
    IF (PRESENT(p_DData) .AND. rtree%isizeDble > 0) THEN
      CALL storage_getsize(rtree%h_DData,Isize)
      IF ((rtree%isizeDble .NE. SIZE(p_DData,1)) .OR. (rtree%NA .NE. SIZE(p_DData,2))) THEN
        PRINT *, "t_tree_copyto_arrayInt: Invalid auxiliary data!"
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('t_tree_copyto_arrayInt',rtree%NNA,rtree%h_DData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_double2D(rtree%h_DData,rtree%DData)
      END IF
      ! Copy data
      rtree%DData(:,1:rtree%NA) = p_DData
    END IF
    
    IF (PRESENT(p_FData) .AND. rtree%isizeSngl > 0) THEN
      CALL storage_getsize(rtree%h_FData,Isize)
      IF ((rtree%isizeSngl .NE. SIZE(p_FData,1)) .OR. (rtree%NA .NE. SIZE(p_FData,2))) THEN
        PRINT *, "t_tree_copyto_arrayInt: Invalid auxiliary data!"
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('t_tree_copyto_arrayInt',rtree%NNA,rtree%h_FData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_single2D(rtree%h_FData,rtree%FData)
      END IF
      ! Copy data
      rtree%FData(:,1:rtree%NA) = p_FData
    END IF

    IF (PRESENT(p_IData) .AND. rtree%isizeInt > 0) THEN
      CALL storage_getsize(rtree%h_IData,Isize)
      IF ((rtree%isizeInt .NE. SIZE(p_IData,1)) .OR. (rtree%NA .NE. SIZE(p_IData,2))) THEN
        PRINT *, "t_tree_copyto_arrayInt: Invalid auxiliary data!"
        CALL sys_halt()
      END IF
      ! Do we have to reallocate memory?
      IF (Isize(2) .NE. rtree%NNA) THEN
        CALL storage_realloc('t_tree_copyto_arrayInt',rtree%NNA,rtree%h_IData,ST_NEWBLOCK_NOINIT)
        CALL storage_getbase_int2D(rtree%h_IData,rtree%IData)
      END IF
      ! Copy data
      rtree%IData(:,1:rtree%NA) = p_IData
    END IF

  END SUBROUTINE t_tree_copyto_arrayInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyfrom_key_handle(rtree,h_Key)

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
!</inputoutpu>
!</subroutine>
    
    ! local variables
    REAL(DP), DIMENSION(:), POINTER   :: p_DKey
    REAL(SP), DIMENSION(:), POINTER   :: p_FKey
    INTEGER,  DIMENSION(:), POINTER   :: p_IKey
    INTEGER(I32) :: isize
    INTEGER(PREC_TREEIDX) :: j

    ! Check if handle needs to be (re-)allocated
    IF (h_Key .EQ. ST_NOHANDLE) THEN
      CALL storage_new('t_tree_copyfrom_key_handle','Key',rtree%NA,&
          rtree%ctreeFormat,h_Key,ST_NEWBLOCK_NOINIT)
    ELSE
      CALL storage_getsize(h_Key,isize)
      IF (isize < rtree%NA) THEN
        CALL storage_realloc('t_tree_copyfrom_key_handle',rtree%NA,&
            h_Key,ST_NEWBLOCK_ZERO,.FALSE.)
      END IF
    END IF

    ! What kind of key is used?
    SELECT CASE(rtree%ctreeFormat)
    CASE (ST_DOUBLE)
      CALL storage_getbase_double(h_Key,p_DKey); j=0
      CALL inorderKeyDble(rtree%Kchild(TRIGHT,TROOT))
      
    CASE (ST_SINGLE)
      CALL storage_getbase_single(h_Key,p_FKey); j=0
      CALL inorderKeySngl(rtree%Kchild(TRIGHT,TROOT))
      
    CASE (ST_INT)
      CALL storage_getbase_int(h_Key,p_IKey); j=0
      CALL inorderKeyInt(rtree%Kchild(TRIGHT,TROOT))
      
    CASE DEFAULT
      PRINT *, 't_tree_copyfrom_handle: Unsupported data format!'
      CALL sys_halt()
    END SELECT
    
  CONTAINS

    ! Here, the real working routines follow.
    
    !**************************************************************
    ! Copy the content of the Double key to an array
    
    RECURSIVE SUBROUTINE inorderKeyDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderKeyDble(rtree%Kchild(TLEFT,i))
      j=j+1; p_DKey(j)=rtree%DKey(i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderKeyDble(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderKeyDble
    
    !**************************************************************
    ! Copy the content of the Single key to an array

    RECURSIVE SUBROUTINE inorderKeySngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderKeySngl(rtree%Kchild(TLEFT,i))
      j=j+1; p_FKey(j)=rtree%FKey(i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderKeySngl(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderKeySngl

    !**************************************************************
    ! Copy the content of the Integer key to an array

    RECURSIVE SUBROUTINE inorderKeyInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderKeyInt(rtree%Kchild(TLEFT,i))
      j=j+1; p_IKey(j)=rtree%IKey(i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderKeyInt(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderKeyInt
  END SUBROUTINE t_tree_copyfrom_key_handle

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyfrom_key_arrayDble(rtree,p_DKey)

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
!</inputoutpu>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: j

    ! Check data format
    IF (rtree%ctreeFormat .NE. ST_DOUBLE) THEN
      PRINT *, "t_tree_copyfrom_key_arrayDble: Invalid data format!"
      CALL sys_halt()
    END IF

    ! Check size of array
    IF (SIZE(p_DKey) < rtree%NA) THEN
      PRINT *, "t_tree_copyfrom_key_arrayDble: Array too small!"
      CALL sys_halt()
    END IF
    
    j=0
    CALL inorderKeyDble(rtree%Kchild(TRIGHT,TROOT))

  CONTAINS
    
    ! Here, the real working routine follows
    
    !**************************************************************
    ! Copy the content of the Double key to an array
    
    RECURSIVE SUBROUTINE inorderKeyDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderKeyDble(rtree%Kchild(TLEFT,i))
      j=j+1; p_DKey(j)=rtree%DKey(i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderKeyDble(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderKeyDble
  END SUBROUTINE t_tree_copyfrom_key_arrayDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyfrom_key_arraySngl(rtree,p_FKey)

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
!</inputoutpu>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: j

    ! Check data format
    IF (rtree%ctreeFormat .NE. ST_SINGLE) THEN
      PRINT *, "t_tree_copyfrom_key_arraySngl: Invalid data format!"
      CALL sys_halt()
    END IF

    ! Check size of array
    IF (SIZE(p_FKey) < rtree%NA) THEN
      PRINT *, "t_tree_copyfrom_key_arraySngl: Array too small!"
      CALL sys_halt()
    END IF
    
    j=0
    CALL inorderKeySngl(rtree%Kchild(TRIGHT,TROOT))

  CONTAINS

    ! Here, the real working routine follows

    !**************************************************************
    ! Copy the content of the Single key to an array

    RECURSIVE SUBROUTINE inorderKeySngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderKeySngl(rtree%Kchild(TLEFT,i))
      j=j+1; p_FKey(j)=rtree%FKey(i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderKeySngl(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderKeySngl
  END SUBROUTINE t_tree_copyfrom_key_arraySngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyfrom_key_arrayInt(rtree,p_IKey)

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
!</inputoutpu>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: j

    ! Check data format
    IF (rtree%ctreeFormat .NE. ST_INT) THEN
      PRINT *, "t_tree_copyfrom_key_arrayInt: Invalid data format!"
      CALL sys_halt()
    END IF

    ! Check size of array
    IF (SIZE(p_IKey) < rtree%NA) THEN
      PRINT *, "t_tree_copyfrom_key_arrayInt: Array too small!"
      CALL sys_halt()
    END IF
    
    j=0
    CALL inorderKeyInt(rtree%Kchild(TRIGHT,TROOT))

  CONTAINS
    
    ! Here, the real working routine follows

    !**************************************************************
    ! Copy the content of the Integer key to an array
    
    RECURSIVE SUBROUTINE inorderKeyInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderKeyInt(rtree%Kchild(TLEFT,i))
      j=j+1; p_IKey(j)=rtree%IKey(i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderKeyInt(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderKeyInt
  END SUBROUTINE t_tree_copyfrom_key_arrayInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyfrom_handle(rtree,ctype,h_Data,mask)

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
!</inputoutpu>
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
          CALL storage_new('t_tree_copyfrom_handle','DData',Isize,&
              ST_DOUBLE,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_double2D(h_Data,p_DData2D)
          CALL t_tree_copyfrom_arrayDble2D(rtree,p_DData2D,mask)
        ELSEIF (PRESENT(mask)) THEN
          CALL storage_new('t_tree_copyfrom_handle','DData',Isize(2),&
              ST_DOUBLE,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_double(h_Data,p_DData)
          CALL t_tree_copyfrom_arrayDble(rtree,p_DData,mask(1))
        ELSE
          CALL storage_new('t_tree_copyfrom_handle','DData',Isize(2),&
              ST_DOUBLE,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_double(h_Data,p_DData)
          CALL t_tree_copyfrom_arrayDble(rtree,p_DData,1)
        END IF

      CASE (ST_SINGLE)
        
        ! Get first working dimension
        Isize(1)=rtree%isizeSngl
        IF (PRESENT(mask)) Isize(1) = SIZE(mask)

        IF (Isize(1) > 1) THEN
          CALL storage_new('t_tree_copyfrom_handle','FData',Isize,&
              ST_SINGLE,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_single2D(h_Data,p_FData2D)
          CALL t_tree_copyfrom_arraySngl2D(rtree,p_FData2D,mask)
        ELSEIF (PRESENT(mask)) THEN
          CALL storage_new('t_tree_copyfrom_handle','FData',Isize(2),&
              ST_SINGLE,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_single(h_Data,p_FData)
          CALL t_tree_copyfrom_arraySngl(rtree,p_FData,mask(1))
        ELSE
          CALL storage_new('t_tree_copyfrom_handle','FData',Isize(2),&
              ST_SINGLE,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_single(h_Data,p_FData)
          CALL t_tree_copyfrom_arraySngl(rtree,p_FData,1)
        END IF

      CASE (ST_INT)
        
        ! Get first working dimension
        Isize(1)=rtree%isizeInt
        IF (PRESENT(mask)) Isize(1) = SIZE(mask)

        IF (Isize(1) > 1) THEN
          CALL storage_new('t_tree_copyfrom_handle','IData',Isize,&
              ST_INT,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int2D(h_Data,p_IData2D)
          CALL t_tree_copyfrom_arrayInt2D(rtree,p_IData2D,mask)
        ELSEIF (PRESENT(mask)) THEN
          CALL storage_new('t_tree_copyfrom_handle','IData',Isize(2),&
              ST_INT,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int(h_Data,p_IData)
          CALL t_tree_copyfrom_arrayInt(rtree,p_IData,mask(1))
        ELSE
          CALL storage_new('t_tree_copyfrom_handle','IData',Isize(2),&
              ST_INT,h_Data,ST_NEWBLOCK_NOINIT)
          CALL storage_getbase_int(h_Data,p_IData)
          CALL t_tree_copyfrom_arrayInt(rtree,p_IData,1)
        END IF
        
      CASE DEFAULT
        PRINT *, "t_tree_copyfrom_handle: Unsupported data format!"
        CALL sys_halt()
      END SELECT

    ELSE   ! The handle is already associated

      ! Check if data type is valid
      CALL storage_getdatatype(h_Data,idatatype)
      IF (idatatype .NE. ctype) THEN
        PRINT *, "t_tree_copyfrom_handle: Data type mismatch!"
        CALL sys_halt()
      END IF

      ! Are we 1D- or 2D-array?
      CALL storage_getdimension(h_Data,idimension)
      SELECT CASE(idimension)

      CASE (1)   !!! 1D-ARRAY  !!!

        ! Do we have to reallocate the array?
        CALL storage_getsize(h_Data,iisize)
        IF (iisize < rtree%NA)&
            CALL storage_realloc('t_tree_copyfrom_handle',rtree%NA,&
            h_Data,ST_NEWBLOCK_NOINIT,.FALSE.)

        ! What kind of data type are we?
        SELECT CASE(ctype)
        CASE (ST_DOUBLE)
          CALL storage_getbase_double(h_Data,p_DData)

          ! Which component should be copied
          IF (PRESENT(mask)) THEN
            IF (SIZE(mask) > 1) THEN
              PRINT *, "t_tree_copyfrom_handle: For 1D-array mask can only&
                  & have one entry!"
              CALL sys_halt()
            END IF
            CALL t_tree_copyfrom_arrayDble(rtree,p_DData,mask(1))
          ELSEIF (rtree%isizeDble .EQ. 1) THEN
            CALL t_tree_copyfrom_arrayDble(rtree,p_DData,1)
          ELSE
            PRINT *, "t_tree_copyfrom_handle: A 1D-array was given but there&
                & are more than one components in the tree!"
            CALL sys_halt()
          END IF
                    
        CASE (ST_SINGLE)
          CALL storage_getbase_single(h_Data,p_FData)

          ! Which component should be copied
          IF (PRESENT(mask)) THEN
            IF (SIZE(mask) > 1) THEN
              PRINT *, "t_tree_copyfrom_handle: For 1D-array mask can only&
                  & have one entry!"
              CALL sys_halt()
            END IF
            CALL t_tree_copyfrom_arraySngl(rtree,p_FData,mask(1))
          ELSEIF (rtree%isizeSngl .EQ. 1) THEN
            CALL t_tree_copyfrom_arraySngl(rtree,p_FData,1)
          ELSE
            PRINT *, "t_tree_copyfrom_handle: A 1D-array was given but there&
                & are more than one components in the tree!"
            CALL sys_halt()
          END IF
          
        CASE (ST_INT)
          CALL storage_getbase_int(h_Data,p_IData)

          ! Which component should be copied
          IF (PRESENT(mask)) THEN
            IF (SIZE(mask) > 1) THEN
              PRINT *, "t_tree_copyfrom_handle: For 1D-array mask can only&
                  & have one entry!"
              CALL sys_halt()
            END IF
            CALL t_tree_copyfrom_arrayInt(rtree,p_IData,mask(1))
          ELSEIF (rtree%isizeInt .EQ. 1) THEN
            CALL t_tree_copyfrom_arrayInt(rtree,p_IData,1)
          ELSE
            PRINT *, "t_tree_copyfrom_handle: A 1D-array was given but there&
                & are more than one components in the tree!"
            CALL sys_halt()
          END IF

        CASE DEFAULT
          PRINT *, "t_tree_copyfrom_handle: Unsupported data type!"
          CALL sys_halt()
        END SELECT

      CASE (2)   !!! 2D-ARRAY !!!

        ! Do we have to reallocate the array?
        CALL storage_getsize(h_Data,Isize)
        IF (Isize(2) < rtree%NA)&
            CALL storage_realloc('t_tree_copyfrom_handle',rtree%NA,&
            h_Data,ST_NEWBLOCK_NOINIT,.FALSE.)

        ! What kind of data type are we
        SELECT CASE(ctype)
        CASE (ST_DOUBLE)
          CALL storage_getbase_double2D(h_Data,p_DData2D)
          CALL t_tree_copyfrom_arrayDble2D(rtree,p_DData2D,mask)

        CASE (ST_SINGLE)
          CALL storage_getbase_single2D(h_Data,p_FData2D)
          CALL t_tree_copyfrom_arraySngl2D(rtree,p_FData2D,mask)

        CASE (ST_INT)
          CALL storage_getbase_int2D(h_Data,p_IData2D)
          CALL t_tree_copyfrom_arrayInt2D(rtree,p_IData2D,mask)

        CASE DEFAULT
          PRINT *, "t_tree_copyfrom_handle: Unsupported data type!"
          CALL sys_halt()
        END SELECT

      CASE DEFAULT
        PRINT *, "t_tree_copyfrom_handle: Unsupported data dimension!"
        CALL sys_halt()
      END SELECT
    END IF

  END SUBROUTINE t_tree_copyfrom_handle

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyfrom_arrayDble(rtree,p_Ddata,mask)

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
!</inputoutpu>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: j

    ! Check if auxiliary double data is available
    IF (rtree%isizeDble .EQ. 0) THEN
      PRINT *, "t_tree_copyfrom_arrayDble: No double data available!"
      CALL sys_halt()
    END IF

    ! Check if given array is large enough in its second dimension
    IF (SIZE(p_DData) < rtree%NA) THEN
      PRINT *, "t_tree_copyfrom_arrayDble: Array too small!"
      CALL sys_halt()
    END IF

    ! Check if mask is valid
    IF (mask < 1 .OR. mask > rtree%isizeDble) THEN
      PRINT *, "t_tree_copyfrom_arrayDble: Invalid mask!"
      CALL sys_halt()
    END IF
    
    ! Do the work
    j=0; CALL inorderDataDble(rtree%Kchild(TRIGHT,TROOT))

  CONTAINS

    ! Here, the real working routine follows.
    
    !**************************************************************
    ! Copy the content of the Double data to an array
    
    RECURSIVE SUBROUTINE inorderDataDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderDataDble(rtree%Kchild(TLEFT,i))
      j=j+1; p_DData(j)=rtree%DData(mask,i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderDataDble(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataDble

  END SUBROUTINE t_tree_copyfrom_arrayDble

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyfrom_arrayDble2D(rtree,p_Ddata,mask)

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
!</inputoutpu>
!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_TREEIDX) :: j

    ! Check if auxiliary double data is available
    IF (rtree%isizeDble .EQ. 0) THEN
      PRINT *, "t_tree_copyfrom_arrayDble2D: No double data available!"
      CALL sys_halt()
    END IF

    ! Check if given array is large enough in its second dimension
    Isize = SHAPE(p_Ddata)
    IF (Isize(2) < rtree%NA) THEN
      PRINT *, "t_tree_copyfrom_arrayDble2D: Array too small!"
      CALL sys_halt()
    END IF

    ! Copy the content to the array
    IF (PRESENT(mask)) THEN
      ! Check if given array has correct size in its first dimension
      IF (Isize(1) .NE. SIZE(mask)) THEN
        PRINT *, "t_tree_copyfrom_arrayDble2D: Array dimensions do not match mask!"
        CALL sys_halt()
      END IF

      ! Check if mask is valid
      IF (ANY(mask < 1) .OR. ANY(mask > rtree%isizeDble)) THEN
        PRINT *, "t_tree_copyfrom_arrayDble2D: Invalid mask!"
        CALL sys_halt()
      END IF
      j=0; CALL inorderDataDbleMask(rtree%Kchild(TRIGHT,TROOT))
    ELSE
      ! Check if given array has correct size in its first dimension
      IF (Isize(1) .NE. rtree%isizeDble) THEN
        PRINT *, "t_tree_copyfrom_arrayDble2D: Array too smalle!"
        CALL sys_halt()
      END IF
      j=0; CALL inorderDataDble(rtree%Kchild(TRIGHT,TROOT))
    END IF

  CONTAINS

    ! Here, the real working routines follow.

    !**************************************************************
    ! Copy the content of the Double data to an array

    RECURSIVE SUBROUTINE inorderDataDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderDataDble(rtree%Kchild(TLEFT,i))
      j=j+1; p_DData(:,j)=rtree%DData(:,i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderDataDble(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataDble

    !**************************************************************
    ! Copy the content of the Double data to an array (masked)

    RECURSIVE SUBROUTINE inorderDataDbleMask(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderDataDbleMask(rtree%Kchild(TLEFT,i))
      j=j+1; p_DData(:,j)=rtree%DData(mask,i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderDataDbleMask(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataDbleMask

  END SUBROUTINE t_tree_copyfrom_arrayDble2D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyfrom_arraySngl(rtree,p_FData,mask)

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
!</inputoutpu>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: j

    ! Check if auxiliary single data is available
    IF (rtree%isizeSngl .EQ. 0) THEN
      PRINT *, "t_tree_copyfrom_arraySngl: No single data available!"
      CALL sys_halt()
    END IF

    ! Check if given array is large enough in its second dimension
    IF (SIZE(p_FData) < rtree%NA) THEN
      PRINT *, "t_tree_copyfrom_arraySngl: Array too small!"
      CALL sys_halt()
    END IF

    ! Check if mask is valid
    IF (mask < 1 .OR. mask > rtree%isizeSngl) THEN
      PRINT *, "t_tree_copyfrom_arraySngl: Invalid mask!"
      CALL sys_halt()
    END IF
    
    ! Do the work
    j=0; CALL inorderDataSngl(rtree%Kchild(TRIGHT,TROOT))

  CONTAINS

    ! Here, the real working routine follows.
    
    !**************************************************************
    ! Copy the content of the Single data to an array
    
    RECURSIVE SUBROUTINE inorderDataSngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderDataSngl(rtree%Kchild(TLEFT,i))
      j=j+1; p_FData(j)=rtree%FData(mask,i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderDataSngl(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataSngl

  END SUBROUTINE t_tree_copyfrom_arraySngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyfrom_arraySngl2D(rtree,p_Fdata,mask)

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
!</inputoutpu>
!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_TREEIDX) :: j

    ! Check if auxiliary single data is available
    IF (rtree%isizeSngl .EQ. 0) THEN
      PRINT *, "t_tree_copyfrom_arraySngl2D: No single data available!"
      CALL sys_halt()
    END IF

    ! Check if given array is large enough in its second dimension
    Isize = SHAPE(p_Fdata)
    IF (Isize(2) < rtree%NA) THEN
      PRINT *, "t_tree_copyfrom_arraySngl2D: Array too small!"
      CALL sys_halt()
    END IF

    ! Copy the content to the array
    IF (PRESENT(mask)) THEN
      ! Check if given array has correct size in its first dimension
      IF (Isize(1) .NE. SIZE(mask)) THEN
        PRINT *, "t_tree_copyfrom_arraySngl2D: Array dimensions do not match mask!"
        CALL sys_halt()
      END IF

      ! Check if mask is valid
      IF (ANY(mask < 1) .OR. ANY(mask > rtree%isizeSngl)) THEN
        PRINT *, "t_tree_copyfrom_arraySngl2D: Invalid mask!"
        CALL sys_halt()
      END IF
      j=0; CALL inorderDataSnglMask(rtree%Kchild(TRIGHT,TROOT))
    ELSE
      ! Check if given array has correct size in its first dimension
      IF (Isize(1) .NE. rtree%isizeSngl) THEN
        PRINT *, "t_tree_copyfrom_arraySngl2D: Array too smalle!"
        CALL sys_halt()
      END IF
      j=0; CALL inorderDataSngl(rtree%Kchild(TRIGHT,TROOT))
    END IF

  CONTAINS

    ! Here, the real working routines follow.

    !**************************************************************
    ! Copy the content of the Single data to an array

    RECURSIVE SUBROUTINE inorderDataSngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderDataSngl(rtree%Kchild(TLEFT,i))
      j=j+1; p_FData(:,j)=rtree%FData(:,i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderDataSngl(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataSngl

    !**************************************************************
    ! Copy the content of the Single data to an array (masked)

    RECURSIVE SUBROUTINE inorderDataSnglMask(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderDataSnglMask(rtree%Kchild(TLEFT,i))
      j=j+1; p_FData(:,j)=rtree%FData(mask,i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderDataSnglMask(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataSnglMask

  END SUBROUTINE t_tree_copyfrom_arraySngl2D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyfrom_arrayInt(rtree,p_Idata,mask)

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
!</inputoutpu>
!</subroutine>

    ! local variables
    INTEGER(PREC_TREEIDX) :: j

    ! Check if auxiliary integer data is available
    IF (rtree%isizeInt .EQ. 0) THEN
      PRINT *, "t_tree_copyfrom_arrayInt: No integer data available!"
      CALL sys_halt()
    END IF

    ! Check if given array is large enough in its second dimension
    IF (SIZE(p_IData) < rtree%NA) THEN
      PRINT *, "t_tree_copyfrom_arrayInt: Array too small!"
      CALL sys_halt()
    END IF

    ! Check if mask is valid
    IF (mask < 1 .OR. mask > rtree%isizeInt) THEN
      PRINT *, "t_tree_copyfrom_arrayInt: Invalid mask!"
      CALL sys_halt()
    END IF
    
    ! Do the work
    j=0; CALL inorderDataInt(rtree%Kchild(TRIGHT,TROOT))

  CONTAINS

    ! Here, the real working routine follows.
    
    !**************************************************************
    ! Copy the content of the Integer data to an array
    
    RECURSIVE SUBROUTINE inorderDataInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderDataInt(rtree%Kchild(TLEFT,i))
      j=j+1; p_IData(j)=rtree%IData(mask,i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderDataInt(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataInt

  END SUBROUTINE t_tree_copyfrom_arrayInt

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_copyfrom_arrayInt2D(rtree,p_Idata,mask)

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
!</inputoutpu>
!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(2) :: Isize
    INTEGER(PREC_TREEIDX) :: j

    ! Check if auxiliary integer data is available
    IF (rtree%isizeInt .EQ. 0) THEN
      PRINT *, "t_tree_copyfrom_arrayInt2D: No integer data available!"
      CALL sys_halt()
    END IF

    ! Check if given array is large enough in its second dimension
    Isize = SHAPE(p_IData)
    IF (Isize(2) < rtree%NA) THEN
      PRINT *, "t_tree_copyfrom_arrayInt2D: Array too small!"
      CALL sys_halt()
    END IF

    ! Copy the content to the array
    IF (PRESENT(mask)) THEN
      ! Check if given array has correct size in its first dimension
      IF (Isize(1) .NE. SIZE(mask)) THEN
        PRINT *, "t_tree_copyfrom_arrayInt2D: Array dimensions do not match mask!"
        CALL sys_halt()
      END IF

      ! Check if mask is valid
      IF (ANY(mask < 1) .OR. ANY(mask > rtree%isizeInt)) THEN
        PRINT *, "t_tree_copyfrom_arrayInt2D: Invalid mask!"
        CALL sys_halt()
      END IF
      j=0; CALL inorderDataIntMask(rtree%Kchild(TRIGHT,TROOT))
    ELSE
      ! Check if given array has correct size in its first dimension
      IF (Isize(1) .NE. rtree%isizeInt) THEN
        PRINT *, "t_tree_copyfrom_arrayInt2D: Array too smalle!"
        CALL sys_halt()
      END IF
      j=0; CALL inorderDataInt(rtree%Kchild(TRIGHT,TROOT))
    END IF

  CONTAINS

    ! Here, the real working routines follow.

    !**************************************************************
    ! Copy the content of the Integer data to an array

    RECURSIVE SUBROUTINE inorderDataInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderDataInt(rtree%Kchild(TLEFT,i))
      j=j+1; p_IData(:,j)=rtree%IData(:,i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderDataInt(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataInt

    !**************************************************************
    ! Copy the content of the Integer data to an array (masked)

    RECURSIVE SUBROUTINE inorderDataIntMask(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) CALL inorderDataIntMask(rtree%Kchild(TLEFT,i))
      j=j+1; p_IData(:,j)=rtree%IData(mask,i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) CALL inorderDataIntMask(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderDataIntMask

  END SUBROUTINE t_tree_copyfrom_arrayInt2D

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_insertDble(rtree,dkey,DData,FData,IData,iposOpt)

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
      PRINT *, 't_tree_insertDble: Unsupported data format!'
      CALL sys_halt()
    END IF

    ! Check if key is already stored in tree
    IF (search(rtree,dkey,jpos) .EQ. BTREE_NOT_FOUND) THEN
      
      ! Adjust size and position
      rtree%na = rtree%na+1
      ipos     = rtree%Kchild(TFREE,TROOT)

      IF (PRESENT(iposOpt)) iposOpt=ABS(ipos)

      ! Check if memory needs to be expanded
      IF (ABS(ipos) > rtree%nna) &
          CALL resize(rtree,CEILING(rtree%dfactor*rtree%nna))
      
      ! Compute next free position in memory
      IF (ipos > 0) THEN
        rtree%Kchild(TFREE,TROOT) = ipos+1
      ELSE
        ipos                      = ABS(ipos)
        rtree%Kchild(TFREE,TROOT) = rtree%Kchild(TFREE,ipos)
      END IF
      
      ! Store given data in memory
      rtree%DKey(ipos)     = dkey
      rtree%Kbal(ipos)     = 0
      rtree%Kchild(:,ipos) = TNULL
      
      ! Store optional data in memory
      IF ((rtree%isizeInt > 0) .AND. &
          PRESENT(IData)) rtree%IData(:,ipos) = IData

      IF ((rtree%isizeDble > 0) .AND. &
          PRESENT(DData)) rtree%DData(:,ipos) = DData

      IF ((rtree%isizeSngl > 0) .AND. &
          PRESENT(FData)) rtree%FData(:,ipos) = FData
      
      ! Insert new node into the tree and into the search path
      rtree%Kchild(MERGE(TLEFT,TRIGHT,jpos < 0),ABS(jpos)) = ipos
      rtree%Kpath(rtree%depth+1) = MERGE(-ipos,ipos,jpos <= 0)
      
      ! Invoke rebalance procedure
      IF (rtree%depth > 0) CALL rebalanceAfterInsert(rtree,rtree%depth)

    ELSEIF(PRESENT(iposOpt)) THEN
      iposOpt=-rtree%Kchild(MERGE(TLEFT,TRIGHT,jpos < 0),ABS(jpos))
    END IF
  END SUBROUTINE t_tree_insertDble
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_insertSngl(rtree,skey,DData,FData,IData,iposOpt)

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
      PRINT *, 't_tree_insertSngl: Unsupported data format!'
      CALL sys_halt()
    END IF

    ! Check if key is already stored in tree
    IF (search(rtree,skey,jpos) .EQ. BTREE_NOT_FOUND) THEN
      
      ! Adjust size and position
      rtree%na = rtree%na+1
      ipos     = rtree%Kchild(TFREE,TROOT)

      IF (PRESENT(iposOpt)) iposOpt=ABS(ipos)

      ! Check if memory needs to be expanded
      IF (ABS(ipos) > rtree%nna) &
          CALL resize(rtree,CEILING(rtree%dfactor*rtree%nna))
      
      ! Compute next free position in memory
      IF (ipos > 0) THEN
        rtree%Kchild(TFREE,TROOT) = ipos+1
      ELSE
        ipos                      = ABS(ipos)
        rtree%Kchild(TFREE,TROOT) = rtree%Kchild(TFREE,ipos)
      END IF
      
      ! Store given data in memory
      rtree%FKey(ipos)     = skey
      rtree%Kbal(ipos)     = 0
      rtree%Kchild(:,ipos) = TNULL
      
      ! Store optional data in memory
      IF ((rtree%isizeInt > 0) .AND. &
          PRESENT(IData)) rtree%IData(:,ipos) = IData

      IF ((rtree%isizeDble > 0) .AND. &
          PRESENT(DData)) rtree%DData(:,ipos) = DData

      IF ((rtree%isizeSngl > 0) .AND. &
          PRESENT(FData)) rtree%FData(:,ipos) = FData
      
      ! Insert new node into the tree and into the search path
      rtree%Kchild(MERGE(TLEFT,TRIGHT,jpos < 0),ABS(jpos)) = ipos
      rtree%Kpath(rtree%depth+1) = MERGE(-ipos,ipos,jpos <= 0)
      
      ! Invoke rebalance procedure
      IF (rtree%depth > 0) CALL rebalanceAfterInsert(rtree,rtree%depth)

    ELSEIF(PRESENT(iposOpt)) THEN
      iposOpt=-rtree%Kchild(MERGE(TLEFT,TRIGHT,jpos < 0),ABS(jpos))
    END IF
  END SUBROUTINE t_tree_insertSngl

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE t_tree_insertInt(rtree,ikey,DData,FData,IData,iposOpt)

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
      PRINT *, 't_tree_insertInt: Unsupported data format!'
      CALL sys_halt()
    END IF

    ! Check if key is already stored in tree
    IF (search(rtree,ikey,jpos) .EQ. BTREE_NOT_FOUND) THEN
      
      ! Adjust size and position
      rtree%na = rtree%na+1
      ipos     = rtree%Kchild(TFREE,TROOT)

      IF (PRESENT(iposOpt)) iposOpt=ABS(ipos)
      
      ! Check if memory needs to be expanded
      IF (ABS(ipos) > rtree%nna) &
          CALL resize(rtree,CEILING(rtree%dfactor*rtree%nna))
      
      ! Compute next free position in memory
      IF (ipos > 0) THEN
        rtree%Kchild(TFREE,TROOT) = ipos+1
      ELSE
        ipos                      = ABS(ipos)
        rtree%Kchild(TFREE,TROOT) = rtree%Kchild(TFREE,ipos)
      END IF
      
      ! Store given data in memory
      rtree%IKey(ipos)     = ikey
      rtree%Kbal(ipos)     = 0
      rtree%Kchild(:,ipos) = TNULL
      
      ! Store optional data in memory
      IF ((rtree%isizeInt > 0) .AND. &
          PRESENT(IData)) rtree%IData(:,ipos) = IData

      IF ((rtree%isizeDble > 0) .AND. &
          PRESENT(DData)) rtree%DData(:,ipos) = DData

      IF ((rtree%isizeSngl > 0) .AND. &
          PRESENT(FData)) rtree%FData(:,ipos) = FData
      
      ! Insert new node into the tree and into the search path
      rtree%Kchild(MERGE(TLEFT,TRIGHT,jpos < 0),ABS(jpos)) = ipos
      rtree%Kpath(rtree%depth+1) = MERGE(-ipos,ipos,jpos <= 0)
      
      ! Invoke rebalance procedure
      IF (rtree%depth > 0) CALL rebalanceAfterInsert(rtree,rtree%depth)
      
    ELSEIF(PRESENT(iposOpt)) THEN
      iposOpt=-rtree%Kchild(MERGE(TLEFT,TRIGHT,jpos < 0),ABS(jpos))
    END IF
  END SUBROUTINE t_tree_insertInt

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
    
    v   = rtree%Kpath(i)
    dir = MERGE(TLEFT,TRIGHT,v < 0)
    v   = ABS(v)
    
    ! Which rotation needs to be performed?
    SELECT CASE (dir)
      
    CASE (TLEFT)
      ! Node V has been reached from the left
      SELECT CASE(rtree%Kbal(v))
        
      CASE (1) ! bal(v)=1
        rtree%Kbal(v) = 0
        
      CASE (0) ! bal(v)=0
        rtree%Kbal(v) = -1
        IF (v .NE. rtree%Kchild(TRIGHT,TROOT)) &
            CALL rebalanceAfterInsert(rtree,i-1)
        
      CASE (-1) ! bal(v)=-1
        x = rtree%Kchild(TLEFT,v)
        w = rtree%Kchild(TRIGHT,x)
        
        SELECT CASE(rtree%Kbal(x))
          
        CASE (-1,0) ! bal(x)=-1 or bal(x)=0
          ! Single-Right-rotation
          rtree%Kchild(TLEFT,v)  = rtree%Kchild(TRIGHT,x)
          rtree%Kchild(TRIGHT,x) = v
          rtree%Kbal(x)          = rtree%Kbal(x)+1
          rtree%Kbal(v)          = -rtree%Kbal(x)
          
          IF (v .EQ. rtree%Kchild(TRIGHT,TROOT)) THEN
            rtree%Kchild(TRIGHT,TROOT) = x
          ELSE
            rtree%Kchild(MERGE(TLEFT,TRIGHT,rtree%Kpath(i-1) < 0),&
                ABS(rtree%Kpath(i-1))) = x
          END IF
          
        CASE (1) ! bal(x)=1
          ! Double-Left-Right-rotation
          rtree%Kchild(TLEFT,v)  = rtree%Kchild(TRIGHT,w)
          rtree%Kchild(TRIGHT,x) = rtree%Kchild(TLEFT,w)
          rtree%Kchild(TLEFT,w)  = x
          rtree%Kchild(TRIGHT,w) = v
          rtree%Kbal(v)          = -MIN(0,rtree%Kbal(w))
          rtree%Kbal(x)          = -MAX(0,rtree%Kbal(w))
          rtree%Kbal(w)          = 0
          
          IF (v .EQ. rtree%Kchild(TRIGHT,TROOT)) THEN
            rtree%Kchild(TRIGHT,TROOT) = w
          ELSE
            rtree%Kchild(MERGE(TLEFT,TRIGHT,rtree%Kpath(i-1) < 0),&
                ABS(rtree%Kpath(i-1))) = w
          END IF
          
        END SELECT
        
      END SELECT
      
    CASE (TRIGHT)
      ! Node V has been reached from the right
      SELECT CASE(rtree%Kbal(v))
        
      CASE (-1) ! bal(v)=-1
        rtree%Kbal(v) = 0
        
      CASE (0) ! bal(v)=0
        rtree%Kbal(v) = 1
        IF (v .NE. rtree%Kchild(TRIGHT,TROOT))&
            CALL rebalanceAfterInsert(rtree,i-1)
        
      CASE (1) ! bal(v)=1
        x = rtree%Kchild(TRIGHT,v)
        w = rtree%Kchild(TLEFT,x)
        
        SELECT CASE(rtree%Kbal(x))
          
        CASE (0,1) ! bal(x)=0 or bal(x)=1
          ! Single-Left-rotation
          rtree%Kchild(TRIGHT,v) = rtree%Kchild(TLEFT,x)
          rtree%Kchild(TLEFT,x)  = v
          rtree%Kbal(x)          = rtree%Kbal(x)-1
          rtree%Kbal(v)          = -rtree%Kbal(x)
          
          IF (v .EQ. rtree%Kchild(TRIGHT,TROOT)) THEN
            rtree%Kchild(TRIGHT,TROOT) = x
          ELSE
            rtree%Kchild(MERGE(TLEFT,TRIGHT,rtree%Kpath(i-1) < 0),&
                ABS(rtree%Kpath(i-1))) = x
          END IF
          
        CASE (-1) ! bal(x)=-1
          ! Double-Right-Left-rotation
          rtree%Kchild(TRIGHT,v) = rtree%Kchild(TLEFT,w)
          rtree%Kchild(TLEFT,x)  = rtree%Kchild(TRIGHT,w)
          rtree%Kchild(TLEFT,w)  = v
          rtree%Kchild(TRIGHT,w) = x
          rtree%Kbal(v)          = -MAX(0,rtree%Kbal(w))
          rtree%Kbal(x)          = -MIN(0,rtree%Kbal(w))
          rtree%Kbal(w)          = 0
          
          IF (v .EQ. rtree%Kchild(TRIGHT,TROOT)) THEN
            rtree%Kchild(TRIGHT,TROOT) = w
          ELSE
            rtree%Kchild(MERGE(TLEFT,TRIGHT,rtree%Kpath(i-1) < 0),&
                ABS(rtree%Kpath(i-1))) = w
          END IF
          
        END SELECT
        
      END SELECT
      
    END SELECT
  END SUBROUTINE rebalanceAfterInsert

  ! ***************************************************************************

!<function>

  FUNCTION t_tree_deleteDble(rtree,dkey) RESULT(f)

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
      PRINT *, 't_tree_deleteDble: Unsupported data format!'
      CALL sys_halt()
    END IF

    ! Search for key
    f=search(rtree,dkey,ipred)
    IF (f .EQ. BTREE_NOT_FOUND) RETURN

    ! Compute new dimensions
    rtree%na = rtree%na-1
    ipos     = rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))

    ! Check if node to be deleted has two children.
    IF ((rtree%Kchild(TLEFT,ipos) .NE. TNULL) .AND. &
        (rtree%Kchild(TRIGHT,ipos) .NE. TNULL)) THEN
      
      ! Descent in the left subtree of node IPOS always to the right
      ! until a node JPOS  without right child is found and
      ! interchange data.
      jpos  = rtree%Kchild(TLEFT,ipos)
      ipred = -ipos
      DO
        rtree%depth              = rtree%depth+1
        rtree%Kpath(rtree%depth) = ipred
        
        IF (rtree%Kchild(TRIGHT,jpos) .EQ. TNULL) THEN
          ! Change key
          rtree%DKey(ipos) = rtree%DKey(jpos)

          ! Change auxiliary data
          IF (rtree%isizeDble > 0) rtree%DData(:,ipos) = rtree%DData(:,jpos)
          IF (rtree%isizeSngl > 0) rtree%FData(:,ipos) = rtree%FData(:,jpos)
          IF (rtree%isizeInt  > 0) rtree%IData(:,ipos) = rtree%IData(:,jpos)
          EXIT
        END IF
        ipred = jpos
        jpos = rtree%Kchild(TRIGHT,jpos)
      END DO
    END IF

    ! Node to be deleted has less than two childen
    ipos = rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
    
    IF (rtree%Kchild(TLEFT,ipos) .EQ. TNULL) THEN
      IF (rtree%Kchild(TRIGHT,ipos) .EQ. TNULL) THEN
        ! Node to be deleted is leaf: Nullify pointer to node and
        ! mark position in array as deleted
        rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) = TNULL
        rtree%Kchild(TFREE,ipos)  = rtree%Kchild(TFREE,TROOT)
        rtree%Kchild(TFREE,TROOT) = -ipos
      ELSE
        ! Node to be deleted has right child: Set pointer to right
        ! child and mark position in array as deleted
        rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) =&
            rtree%Kchild(TRIGHT,ipos)
        rtree%Kchild(TFREE,ipos)  = rtree%Kchild(TFREE,TROOT)
        rtree%Kchild(TFREE,TROOT) = -ipos
      END IF
    ELSE
      ! Node to be deleted has left child: Set pointer to left child
      ! and mark position in array as deleted
      rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) =&
          rtree%Kchild(TLEFT,ipos)
      rtree%Kchild(TFREE,ipos)  = rtree%Kchild(TFREE,TROOT)
      rtree%Kchild(TFREE,TROOT) = -ipos
    END IF
    
    ! Invoke rebalance procedure
    IF (rtree%depth > 0) CALL rebalanceAfterDeletion(rtree,rtree%depth)
  END FUNCTION t_tree_deleteDble

  ! ***************************************************************************

!<function>

  FUNCTION t_tree_deleteSngl(rtree,skey) RESULT(f)

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
      PRINT *, 't_tree_deleteSngl: Unsupported data format!'
      CALL sys_halt()
    END IF

    ! Search for key
    f=search(rtree,skey,ipred)
    IF (f .EQ. BTREE_NOT_FOUND) RETURN

    ! Compute new dimensions
    rtree%na = rtree%na-1
    ipos     = rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))

    ! Check if node to be deleted has two children.
    IF ((rtree%Kchild(TLEFT,ipos) .NE. TNULL) .AND. &
        (rtree%Kchild(TRIGHT,ipos) .NE. TNULL)) THEN
      
      ! Descent in the left subtree of node IPOS always to the right
      ! until a node JPOS  without right child is found and
      ! interchange data.
      jpos  = rtree%Kchild(TLEFT,ipos)
      ipred = -ipos
      DO
        rtree%depth              = rtree%depth+1
        rtree%Kpath(rtree%depth) = ipred
        
        IF (rtree%Kchild(TRIGHT,jpos) .EQ. TNULL) THEN
          ! Change key
          rtree%FKey(ipos) = rtree%FKey(jpos)

          ! Change auxiliary data
          IF (rtree%isizeDble > 0) rtree%DData(:,ipos) = rtree%DData(:,jpos)
          IF (rtree%isizeSngl > 0) rtree%FData(:,ipos) = rtree%FData(:,jpos)
          IF (rtree%isizeInt  > 0) rtree%IData(:,ipos) = rtree%IData(:,jpos)
          EXIT
        END IF
        ipred = jpos
        jpos = rtree%Kchild(TRIGHT,jpos)
      END DO
    END IF

    ! Node to be deleted has less than two childen
    ipos = rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
    
    IF (rtree%Kchild(TLEFT,ipos) .EQ. TNULL) THEN
      IF (rtree%Kchild(TRIGHT,ipos) .EQ. TNULL) THEN
        ! Node to be deleted is leaf: Nullify pointer to node and
        ! mark position in array as deleted
        rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) = TNULL
        rtree%Kchild(TFREE,ipos)  = rtree%Kchild(TFREE,TROOT)
        rtree%Kchild(TFREE,TROOT) = -ipos
      ELSE
        ! Node to be deleted has right child: Set pointer to right
        ! child and mark position in array as deleted
        rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) =&
            rtree%Kchild(TRIGHT,ipos)
        rtree%Kchild(TFREE,ipos)  = rtree%Kchild(TFREE,TROOT)
        rtree%Kchild(TFREE,TROOT) = -ipos
      END IF
    ELSE
      ! Node to be deleted has left child: Set pointer to left child
      ! and mark position in array as deleted
      rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) =&
          rtree%Kchild(TLEFT,ipos)
      rtree%Kchild(TFREE,ipos)  = rtree%Kchild(TFREE,TROOT)
      rtree%Kchild(TFREE,TROOT) = -ipos
    END IF
    
    ! Invoke rebalance procedure
    IF (rtree%depth > 0) CALL rebalanceAfterDeletion(rtree,rtree%depth)
  END FUNCTION t_tree_deleteSngl

  ! ***************************************************************************

!<function>

  FUNCTION t_tree_deleteInt(rtree,ikey) RESULT(f)

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
      PRINT *, 't_tree_deleteInt: Unsupported data format!'
      CALL sys_halt()
    END IF

    ! Search for key
    f=search(rtree,ikey,ipred)
    IF (f .EQ. BTREE_NOT_FOUND) RETURN

    ! Compute new dimensions
    rtree%na = rtree%na-1
    ipos     = rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))

    ! Check if node to be deleted has two children.
    IF ((rtree%Kchild(TLEFT,ipos) .NE. TNULL) .AND. &
        (rtree%Kchild(TRIGHT,ipos) .NE. TNULL)) THEN
      
      ! Descent in the left subtree of node IPOS always to the right
      ! until a node JPOS  without right child is found and
      ! interchange data.
      jpos  = rtree%Kchild(TLEFT,ipos)
      ipred = -ipos
      DO
        rtree%depth              = rtree%depth+1
        rtree%Kpath(rtree%depth) = ipred
        
        IF (rtree%Kchild(TRIGHT,jpos) .EQ. TNULL) THEN
          ! Change key
          rtree%IKey(ipos) = rtree%IKey(jpos)

          ! Change auxiliary data
          IF (rtree%isizeDble > 0) rtree%DData(:,ipos) = rtree%DData(:,jpos)
          IF (rtree%isizeSngl > 0) rtree%FData(:,ipos) = rtree%FData(:,jpos)
          IF (rtree%isizeInt  > 0) rtree%IData(:,ipos) = rtree%IData(:,jpos)
          EXIT
        END IF
        ipred = jpos
        jpos = rtree%Kchild(TRIGHT,jpos)
      END DO
    END IF

    ! Node to be deleted has less than two childen
    ipos = rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
    
    IF (rtree%Kchild(TLEFT,ipos) .EQ. TNULL) THEN
      IF (rtree%Kchild(TRIGHT,ipos) .EQ. TNULL) THEN
        ! Node to be deleted is leaf: Nullify pointer to node and
        ! mark position in array as deleted
        rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) = TNULL
        rtree%Kchild(TFREE,ipos)  = rtree%Kchild(TFREE,TROOT)
        rtree%Kchild(TFREE,TROOT) = -ipos
      ELSE
        ! Node to be deleted has right child: Set pointer to right
        ! child and mark position in array as deleted
        rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) =&
            rtree%Kchild(TRIGHT,ipos)
        rtree%Kchild(TFREE,ipos)  = rtree%Kchild(TFREE,TROOT)
        rtree%Kchild(TFREE,TROOT) = -ipos
      END IF
    ELSE
      ! Node to be deleted has left child: Set pointer to left child
      ! and mark position in array as deleted
      rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred)) =&
          rtree%Kchild(TLEFT,ipos)
      rtree%Kchild(TFREE,ipos)  = rtree%Kchild(TFREE,TROOT)
      rtree%Kchild(TFREE,TROOT) = -ipos
    END IF
    
    ! Invoke rebalance procedure
    IF (rtree%depth > 0) CALL rebalanceAfterDeletion(rtree,rtree%depth)
  END FUNCTION t_tree_deleteInt

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
    
    v   = rtree%Kpath(i)
    dir = MERGE(TLEFT,TRIGHT,v < 0)
    v   = ABS(v)
    
      ! Which rotation needs to be performed?
    SELECT CASE (dir)
      
    CASE (TLEFT)
      ! Node V has been reached from the left
      SELECT CASE (rtree%Kbal(v))
        
      CASE (0) ! bal(v)=0
        rtree%Kbal(v) = 1
        
      CASE (-1) ! bal(v)=-1
        rtree%Kbal(v) = 0
        IF (v .NE. rtree%Kchild(TRIGHT,TROOT))&
            CALL rebalanceAfterDeletion(rtree,i-1)
        
      CASE (1) ! bal(v)=1
        x = rtree%Kchild(TRIGHT,v)
        
        SELECT CASE(rtree%Kbal(x))
          
        CASE (-1) ! bal(x)=-1
          w = rtree%Kchild(TLEFT,x)
          ! Double-Right-Left-rotation
          rtree%Kchild(TRIGHT,v) = rtree%Kchild(TLEFT,w)
          rtree%Kchild(TLEFT,x)  = rtree%Kchild(TRIGHT,w)
          rtree%Kchild(TLEFT,w)  = v
          rtree%Kchild(TRIGHT,w) = x
          rtree%Kbal(v)          = -MAX(0,rtree%Kbal(w))
          rtree%Kbal(x)          = -MIN(0,rtree%Kbal(w))
          rtree%Kbal(w)          = 0
          
          IF (v .EQ. rtree%Kchild(TRIGHT,TROOT)) THEN
            rtree%Kchild(TRIGHT,TROOT) = w
          ELSE
            rtree%Kchild(MERGE(TLEFT,TRIGHT,rtree%Kpath(i-1) < 0),&
                ABS(rtree%Kpath(i-1))) = w
            CALL rebalanceAfterDeletion(rtree,i-1)
          END IF
          
        CASE DEFAULT ! bal(x)=0 or bal(x)=1
          ! Single-Left-rotation
          rtree%Kchild(TRIGHT,v) = rtree%Kchild(TLEFT,x)
          rtree%Kchild(TLEFT,x)  = v
          xbal                   = rtree%Kbal(x)
          rtree%Kbal(x)          = rtree%Kbal(x)-1
          rtree%Kbal(v)          = -rtree%Kbal(x)
          
          IF (v .EQ. rtree%Kchild(TRIGHT,TROOT)) THEN
            rtree%Kchild(TRIGHT,TROOT) = x
          ELSE
            rtree%Kchild(MERGE(TLEFT,TRIGHT,rtree%Kpath(i-1) < 0),&
                ABS(rtree%Kpath(i-1))) = x
            IF (xbal .EQ. 1) CALL rebalanceAfterDeletion(rtree,i-1)
          END IF
          
        END SELECT
        
      END SELECT
      
    CASE (TRIGHT)
      ! Node V has been reached from the right
      SELECT CASE (rtree%Kbal(v))
        
      CASE (0) ! bal(v)=0
        rtree%Kbal(v) = -1
        
      CASE (1) ! bal(v)=1
        rtree%Kbal(v) = 0
        IF (v .NE. rtree%Kchild(TRIGHT,TROOT))&
            CALL rebalanceAfterDeletion(rtree,i-1)
        
      CASE (-1) ! bal(v)=-1
        x=rtree%Kchild(TLEFT,v)
        
        SELECT CASE(rtree%Kbal(x))
          
        CASE (1) ! bal(x)=1
          w = rtree%Kchild(TRIGHT,x)
          ! Double-Left-Right-rotation
          rtree%Kchild(TLEFT,v ) = rtree%Kchild(TRIGHT,w)
          rtree%Kchild(TRIGHT,x) = rtree%Kchild(TLEFT,w)
          rtree%Kchild(TLEFT,w)  = x
          rtree%Kchild(TRIGHT,w) = v
          rtree%Kbal(v)          = -MIN(0,rtree%Kbal(w))
          rtree%Kbal(x)          = -MAX(0,rtree%Kbal(w))
          rtree%Kbal(w)          = 0
          
          IF (v .EQ. rtree%Kchild(TRIGHT,TROOT)) THEN
            rtree%Kchild(TRIGHT,TROOT) = w
          ELSE
            rtree%Kchild(MERGE(TLEFT,TRIGHT,rtree%Kpath(i-1) < 0),&
                ABS(rtree%Kpath(i-1))) = w
            CALL rebalanceAfterDeletion(rtree,i-1)
          END IF
          
        CASE DEFAULT ! bal(x)=0 or bal(x)=-1
          ! Single-Right-rotation
          rtree%Kchild(TLEFT,v)  = rtree%Kchild(TRIGHT,x)
          rtree%Kchild(TRIGHT,x) = v
          xbal                   = rtree%Kbal(x)
          rtree%Kbal(x)          = rtree%Kbal(x)+1
          rtree%Kbal(v)          = -rtree%Kbal(x)
          
          IF (v .EQ. rtree%Kchild(TRIGHT,TROOT)) THEN
            rtree%Kchild(TRIGHT,TROOT) = x
          ELSE
            rtree%Kchild(MERGE(TLEFT,TRIGHT,rtree%Kpath(i-1) < 0),&
                ABS(rtree%Kpath(i-1))) = x
            IF (xbal .EQ. -1) CALL rebalanceAfterDeletion(rtree,i-1)
          END IF
          
        END SELECT
        
      END SELECT
      
    END SELECT
  END SUBROUTINE rebalanceAfterDeletion

  ! ***************************************************************************
  
!<function>

  FUNCTION t_tree_searchDble(rtree,dkey,ipos) RESULT(f)
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
      PRINT *, 't_tree_searchDble: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    f           = BTREE_NOT_FOUND
    ipos        = TROOT
    dir         = TRIGHT
    jpos        = rtree%Kchild(dir,ipos)
    rtree%depth = 0
    
    search: DO
      IF (jpos .EQ. TNULL) THEN
        ipos = MERGE(-ipos,ipos,dir .EQ. TLEFT)
        EXIT search
      END IF
      
      IF (rtree%DKey(jpos) .EQ. dkey) THEN
        f    = BTREE_FOUND
        ipos = MERGE(-ipos,ipos,dir .EQ. TLEFT)
        EXIT search
      END IF
      
      ipos                     = jpos
      dir                      = MERGE(TLEFT,TRIGHT,rtree%DKey(ipos) > dkey)
      jpos                     = rtree%Kchild(dir,ipos)
      rtree%depth              = rtree%depth+1
      rtree%Kpath(rtree%depth) = MERGE(-ipos,ipos,dir .EQ. TLEFT)
    END DO search
  END FUNCTION t_tree_searchDble

  ! ***************************************************************************
  
!<function>

  FUNCTION t_tree_searchSngl(rtree,skey,ipos) RESULT(f)
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
      PRINT *, 't_tree_searchSngl: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    f           = BTREE_NOT_FOUND
    ipos        = TROOT
    dir         = TRIGHT
    jpos        = rtree%Kchild(dir,ipos)
    rtree%depth = 0
    
    search: DO
      IF (jpos .EQ. TNULL) THEN
        ipos = MERGE(-ipos,ipos,dir .EQ. TLEFT)
        EXIT search
      END IF
      
      IF (rtree%FKey(jpos) .EQ. skey) THEN
        f    = BTREE_FOUND
        ipos = MERGE(-ipos,ipos,dir .EQ. TLEFT)
        EXIT search
      END IF
      
      ipos                     = jpos
      dir                      = MERGE(TLEFT,TRIGHT,rtree%FKey(ipos) > skey)
      jpos                     = rtree%Kchild(dir,ipos)
      rtree%depth              = rtree%depth+1
      rtree%Kpath(rtree%depth) = MERGE(-ipos,ipos,dir .EQ. TLEFT)
    END DO search
  END FUNCTION t_tree_searchSngl

  ! ***************************************************************************
  
!<function>

  FUNCTION t_tree_searchInt(rtree,ikey,ipos) RESULT(f)
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
      PRINT *, 't_tree_searchInt: Unsupported data format!'
      CALL sys_halt()
    END IF
    
    f           = BTREE_NOT_FOUND
    ipos        = TROOT
    dir         = TRIGHT
    jpos        = rtree%Kchild(dir,ipos)
    rtree%depth = 0
    
    search: DO
      IF (jpos .EQ. TNULL) THEN
        ipos = MERGE(-ipos,ipos,dir .EQ. TLEFT)
        EXIT search
      END IF
      
      IF (rtree%IKey(jpos) .EQ. ikey) THEN
        f    = BTREE_FOUND
        ipos = MERGE(-ipos,ipos,dir .EQ. TLEFT)
        EXIT search
      END IF
      
      ipos                     = jpos
      dir                      = MERGE(TLEFT,TRIGHT,rtree%IKey(ipos) > ikey)
      jpos                     = rtree%Kchild(dir,ipos)
      rtree%depth              = rtree%depth+1
      rtree%Kpath(rtree%depth) = MERGE(-ipos,ipos,dir .EQ. TLEFT)
    END DO search
  END FUNCTION t_tree_searchInt

  ! ***************************************************************************

!<function>

  FUNCTION t_tree_getitemDble(rtree,dkey) RESULT(ipos)

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
    
    IF (search(rtree,dkey,ipred) .NE. BTREE_FOUND) THEN
      PRINT *, "t_tree_getitemDble: Unable to find item in tree"
      CALL sys_halt()
    END IF
    ipos = rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
  END FUNCTION t_tree_getitemDble

  ! ***************************************************************************

!<function>

  FUNCTION t_tree_getitemSngl(rtree,skey) RESULT(ipos)

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
    
    IF (search(rtree,skey,ipred) .NE. BTREE_FOUND) THEN
      PRINT *, "t_tree_getitemSngl: Unable to find item in tree"
      CALL sys_halt()
    END IF
    ipos = rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))
  END FUNCTION t_tree_getitemSngl
  
  ! ***************************************************************************

!<function>

  FUNCTION t_tree_getitemInt(rtree,ikey) RESULT(ipos)

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
    
    IF (search(rtree,ikey,ipred) .NE. BTREE_FOUND) THEN
      PRINT *, "t_tree_getitemInt: Unable to find item in tree"
      CALL sys_halt()
    END IF
    ipos = rtree%Kchild(MERGE(TLEFT,TRIGHT,ipred < 0),ABS(ipred))

    IF (ipos .NE. ikey) THEN
      WRITE(*,'(A,I5,A,I5,A)') 'Position',ipos,'and key',ikey,'are different'
      CALL sys_halt()
    END IF
  END FUNCTION t_tree_getitemInt

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE t_tree_print(rtree,op)

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
      IF (rtree%Kchild(TRIGHT,TROOT) .NE. TNULL) THEN
        SELECT CASE(rtree%ctreeFormat)
        CASE (ST_DOUBLE)
          CALL preorderDble(rtree%Kchild(TRIGHT,TROOT))
        CASE (ST_SINGLE)
          CALL preorderSngl(rtree%Kchild(TRIGHT,TROOT))
        CASE (ST_INT)
          CALL preorderInt(rtree%Kchild(TRIGHT,TROOT))
        CASE DEFAULT
          PRINT *, 't_tree_print: Unsupported data format!'
          CALL sys_halt()
        END SELECT
      END IF
      
    CASE (BTREE_INORDER)
      IF (rtree%Kchild(TRIGHT,TROOT) .NE. TNULL) THEN
        SELECT CASE(rtree%ctreeFormat)
        CASE (ST_DOUBLE)
          CALL inorderDble(rtree%Kchild(TRIGHT,TROOT))
        CASE (ST_SINGLE)
          CALL inorderSngl(rtree%Kchild(TRIGHT,TROOT))
        CASE (ST_INT)
          CALL inorderInt(rtree%Kchild(TRIGHT,TROOT))
        CASE DEFAULT
          PRINT *, 't_tree_print: Unsupported data format!'
          CALL sys_halt()
        END SELECT
      END IF
      
    CASE (BTREE_POSTORDER)
      IF (rtree%Kchild(TRIGHT,TROOT) .NE. TNULL) THEN
        SELECT CASE(rtree%ctreeFormat)
        CASE (ST_DOUBLE)
          CALL postorderDble(rtree%Kchild(TRIGHT,TROOT))
        CASE (ST_SINGLE)
          CALL postorderSngl(rtree%Kchild(TRIGHT,TROOT))
        CASE (ST_INT)
          CALL postorderInt(rtree%Kchild(TRIGHT,TROOT))
        CASE DEFAULT
          PRINT *, 't_tree_print: Unsupported data format!'
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
      
      PRINT *, rtree%DKey(i)
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL)&
          CALL preorderDble(rtree%Kchild(TLEFT,i))
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL)&
          CALL preorderDble(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE preorderDble

    !**************************************************************
    ! Print the content of the tree in preorder for Single key
    
    RECURSIVE SUBROUTINE preorderSngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      PRINT *, rtree%FKey(i)
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL)&
          CALL preorderSngl(rtree%Kchild(TLEFT,i))
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL)&
          CALL preorderSngl(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE preorderSngl

    !**************************************************************
    ! Print the content of the tree in preorder for Integer key
    
    RECURSIVE SUBROUTINE preorderInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      PRINT *, rtree%IKey(i)
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL)&
          CALL preorderInt(rtree%Kchild(TLEFT,i))
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL)&
          CALL preorderInt(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE preorderInt

    !**************************************************************
    ! Print the content of the tree in postorder for Double key
    
    RECURSIVE SUBROUTINE postorderDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL)&
          CALL postorderDble(rtree%Kchild(TLEFT,i))
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL)&
          CALL postorderDble(rtree%Kchild(TRIGHT,i))
      PRINT *, rtree%DKey(i)
    END SUBROUTINE postorderDble

    !**************************************************************
    ! Print the content of the tree in postorder for Single key
    
    RECURSIVE SUBROUTINE postorderSngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL)&
          CALL postorderSngl(rtree%Kchild(TLEFT,i))
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL)&
          CALL postorderSngl(rtree%Kchild(TRIGHT,i))
      PRINT *, rtree%FKey(i)
    END SUBROUTINE postorderSngl

    !**************************************************************
    ! Print the content of the tree in postorder for Integer key
    
    RECURSIVE SUBROUTINE postorderInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL)&
          CALL postorderInt(rtree%Kchild(TLEFT,i))
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL)&
          CALL postorderInt(rtree%Kchild(TRIGHT,i))
      PRINT *, rtree%IKey(i)
    END SUBROUTINE postorderInt

    !**************************************************************
    ! Print the content of the tree in inorder for Double key
    
    RECURSIVE SUBROUTINE inorderDble(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderDble(rtree%Kchild(TLEFT,i))
      PRINT *, rtree%DKey(i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderDble(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderDble

    !**************************************************************
    ! Print the content of the tree in inorder for Single key
    
    RECURSIVE SUBROUTINE inorderSngl(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderSngl(rtree%Kchild(TLEFT,i))
      PRINT *, rtree%FKey(i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderSngl(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderSngl

    !**************************************************************
    ! Print the content of the tree in inorder for Integer key
    
    RECURSIVE SUBROUTINE inorderInt(i)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL)&
          CALL inorderInt(rtree%Kchild(TLEFT,i))
      PRINT *, rtree%IKey(i)
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL)&
          CALL inorderInt(rtree%Kchild(TRIGHT,i))
    END SUBROUTINE inorderInt
  END SUBROUTINE t_tree_print

  ! ***************************************************************************
  
!<function>
  
  PURE FUNCTION t_tree_height(rtree) RESULT(h)

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
    
    h = height(rtree%Kchild(TRIGHT,TROOT))
    
  CONTAINS
    
    PURE RECURSIVE FUNCTION height(i) RESULT(h)
      INTEGER(PREC_TREEIDX), INTENT(IN) :: i
      INTEGER(PREC_TREEIDX) :: h,hl,hr
      
      IF (rtree%Kchild(TLEFT,i) .NE. TNULL) THEN
        hl = height(rtree%Kchild(TLEFT,i))
      ELSE
        hl = 0
      END IF
      
      IF (rtree%Kchild(TRIGHT,i) .NE. TNULL) THEN
        hr = height(rtree%Kchild(TRIGHT,i))
      ELSE
        hr = 0
      END IF
      
      h = MAX(hl,hr)+1
    END FUNCTION height
  END FUNCTION t_tree_height

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE t_tree_info(rtree)

!<description>
    ! This subroutine outputs statistical info about the tree
!</description>

!<input>
    ! tree
    TYPE(t_btree), INTENT(IN) :: rtree
!</input>
!</subroutine>

    WRITE(*,FMT=*) ' Tree:'
    WRITE(*,FMT=*) ' ====='
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Key    =',rtree%h_Key
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Kbal   =',rtree%h_Kbal
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Kpath  =',rtree%h_Kpath
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_Kchild =',rtree%h_Kchild
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_DData  =',rtree%h_DData
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_FData  =',rtree%h_FData
    WRITE(*,FMT='(1X,A,1X,I5)') '  h_IData  =',rtree%h_IData
    WRITE(*,*)
    WRITE(*,FMT='(1X,A,1X,I8,3X,A,1X,I8,3X,A,F5.1,A)')&
        '  NA      =',rtree%NA,&
        'NNA      =',rtree%NNA,&
        'FILLING  =',100*rtree%NA/REAL(rtree%NNA,DP),'%'
    WRITE(*,FMT='(1X,A,1X,I8,A,2X,F7.1,A)') '  NRESIZE =',&
        rtree%NRESIZE, '   NNA´s    =',100*rtree%NNA/REAL(rtree%NNA0,DP),'%'
    WRITE(*,*)
  END SUBROUTINE t_tree_info

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
